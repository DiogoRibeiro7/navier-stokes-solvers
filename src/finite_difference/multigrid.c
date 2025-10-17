#include "multigrid.h"
#include "ns_fd_solver.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    int nx, ny;
    double dx, dy;
    double *phi;
    double *rhs;
    double *residual;
    double *scratch;
} NSMGLevel;

typedef struct {
    int num_levels;
    NSMGLevel *levels;
} NSMGHierarchy;

static int ns_mg_compute_levels(int nx, int ny) {
    int levels = 1;
    while (nx > 3 && ny > 3) {
        int interior_x = nx - 2;
        int interior_y = ny - 2;
        if (interior_x <= 1 || interior_y <= 1) {
            break;
        }
        interior_x = (interior_x + 1) / 2;
        interior_y = (interior_y + 1) / 2;
        nx = interior_x + 2;
        ny = interior_y + 2;
        levels++;
    }
    return levels;
}

static void ns_mg_free_level(NSMGLevel *level) {
    if (!level) {
        return;
    }
    free(level->phi);
    free(level->rhs);
    free(level->residual);
    free(level->scratch);
    level->phi = level->rhs = level->residual = level->scratch = NULL;
}

static void ns_mg_free_hierarchy(NSMGHierarchy *hier) {
    if (!hier || !hier->levels) {
        return;
    }
    for (int i = 0; i < hier->num_levels; ++i) {
        ns_mg_free_level(&hier->levels[i]);
    }
    free(hier->levels);
    hier->levels = NULL;
    hier->num_levels = 0;
}

static int ns_mg_allocate_level(NSMGLevel *level, int nx, int ny, double dx, double dy) {
    const size_t total = (size_t)nx * (size_t)ny;
    level->nx = nx;
    level->ny = ny;
    level->dx = dx;
    level->dy = dy;
    level->phi = (double *)calloc(total, sizeof(double));
    level->rhs = (double *)calloc(total, sizeof(double));
    level->residual = (double *)calloc(total, sizeof(double));
    level->scratch = (double *)calloc(total, sizeof(double));
    if (!level->phi || !level->rhs || !level->residual || !level->scratch) {
        return 0;
    }
    return 1;
}

static int ns_mg_build_hierarchy(NSMGHierarchy *hier,
                                 int finest_nx, int finest_ny,
                                 double finest_dx, double finest_dy) {
    const int levels = ns_mg_compute_levels(finest_nx, finest_ny);
    hier->levels = (NSMGLevel *)calloc((size_t)levels, sizeof(NSMGLevel));
    if (!hier->levels) {
        hier->num_levels = 0;
        return 0;
    }
    hier->num_levels = levels;

    int nx = finest_nx;
    int ny = finest_ny;
    double dx = finest_dx;
    double dy = finest_dy;

    for (int ell = 0; ell < levels; ++ell) {
        if (!ns_mg_allocate_level(&hier->levels[ell], nx, ny, dx, dy)) {
            ns_mg_free_hierarchy(hier);
            return 0;
        }

        if (ell < levels - 1) {
            const int interior_x = nx - 2;
            const int interior_y = ny - 2;
            const int coarse_interior_x = (interior_x + 1) / 2;
            const int coarse_interior_y = (interior_y + 1) / 2;
            nx = coarse_interior_x + 2;
            ny = coarse_interior_y + 2;
            dx *= 2.0;
            dy *= 2.0;
        }
    }
    return 1;
}

static inline int ns_mg_idx(int i, int j, int nx) {
    return i * nx + j;
}

static double ns_mg_eval_boundary_value(const NSBoundaryCondition *bc,
                                        double x, double y, double t,
                                        double fallback) {
    if (!bc) {
        return fallback;
    }

    if (bc->type == BC_CUSTOM && bc->func) {
        double u = 0.0;
        double v = 0.0;
        double p = fallback;
        bc->func(x, y, t, &u, &v, &p, bc->user_data);
        return p;
    }

    if (isfinite(bc->value_p)) {
        return bc->value_p;
    }

    if (bc->type == BC_DIRICHLET) {
        return fallback;
    }

    return fallback;
}

static void ns_fd_multigrid_apply_pressure_bc(const NSFiniteDiffData *data,
                                              NSMGLevel *level) {
    if (!data || !level) {
        return;
    }

    const int nx = level->nx;
    const int ny = level->ny;

    const NSBoundaryCondition *bottom = &data->boundaries[NS_BC_BOTTOM];
    const NSBoundaryCondition *top = &data->boundaries[NS_BC_TOP];
    const NSBoundaryCondition *left = &data->boundaries[NS_BC_LEFT];
    const NSBoundaryCondition *right = &data->boundaries[NS_BC_RIGHT];

    const double t = data->t;

    const int periodic_y = (bottom->type == BC_PERIODIC && top->type == BC_PERIODIC);
    const int periodic_x = (left->type == BC_PERIODIC && right->type == BC_PERIODIC);

    if (periodic_y) {
#pragma omp parallel for
        for (int j = 0; j < nx; ++j) {
            level->phi[ns_mg_idx(0, j, nx)] = level->phi[ns_mg_idx(ny - 2, j, nx)];
            level->phi[ns_mg_idx(ny - 1, j, nx)] = level->phi[ns_mg_idx(1, j, nx)];
        }
    } else {
#pragma omp parallel for
        for (int j = 0; j < nx; ++j) {
            const double x = j * level->dx;
            const int bottom_idx = ns_mg_idx(0, j, nx);
            const int top_idx = ns_mg_idx(ny - 1, j, nx);
            const double neighbor_bottom = level->phi[ns_mg_idx((ny > 1) ? 1 : 0, j, nx)];
            const double neighbor_top = level->phi[ns_mg_idx((ny > 1) ? ny - 2 : ny - 1, j, nx)];

            double bottom_val = neighbor_bottom;
            double top_val = neighbor_top;

            if (bottom->type == BC_CUSTOM || isfinite(bottom->value_p)) {
                bottom_val = ns_mg_eval_boundary_value(bottom, x, 0.0, t, neighbor_bottom);
            } else if (bottom->type == BC_DIRICHLET) {
                bottom_val = 0.0;
            }

            if (top->type == BC_CUSTOM || isfinite(top->value_p)) {
                const double y_top = (ny - 1) * level->dy;
                top_val = ns_mg_eval_boundary_value(top, x, y_top, t, neighbor_top);
            } else if (top->type == BC_DIRICHLET) {
                top_val = 0.0;
            }

            level->phi[bottom_idx] = bottom_val;
            level->phi[top_idx] = top_val;
        }
    }

    if (periodic_x) {
#pragma omp parallel for
        for (int i = 0; i < ny; ++i) {
            level->phi[ns_mg_idx(i, 0, nx)] = level->phi[ns_mg_idx(i, nx - 2, nx)];
            level->phi[ns_mg_idx(i, nx - 1, nx)] = level->phi[ns_mg_idx(i, 1, nx)];
        }
    } else {
#pragma omp parallel for
        for (int i = 0; i < ny; ++i) {
            const double y = i * level->dy;
            const int left_idx = ns_mg_idx(i, 0, nx);
            const int right_idx = ns_mg_idx(i, nx - 1, nx);
            const double neighbor_left = level->phi[ns_mg_idx(i, (nx > 1) ? 1 : 0, nx)];
            const double neighbor_right = level->phi[ns_mg_idx(i, (nx > 1) ? nx - 2 : nx - 1, nx)];

            double left_val = neighbor_left;
            double right_val = neighbor_right;

            if (left->type == BC_CUSTOM || isfinite(left->value_p)) {
                left_val = ns_mg_eval_boundary_value(left, 0.0, y, t, neighbor_left);
            } else if (left->type == BC_DIRICHLET) {
                left_val = 0.0;
            }

            if (right->type == BC_CUSTOM || isfinite(right->value_p)) {
                const double x_right = (nx - 1) * level->dx;
                right_val = ns_mg_eval_boundary_value(right, x_right, y, t, neighbor_right);
            } else if (right->type == BC_DIRICHLET) {
                right_val = 0.0;
            }

            level->phi[left_idx] = left_val;
            level->phi[right_idx] = right_val;
        }
    }

    level->phi[ns_mg_idx(0, 0, nx)] = 0.5 * (level->phi[ns_mg_idx(0, 1, nx)] + level->phi[ns_mg_idx(1, 0, nx)]);
    level->phi[ns_mg_idx(0, nx - 1, nx)] = 0.5 * (level->phi[ns_mg_idx(0, nx - 2, nx)] + level->phi[ns_mg_idx(1, nx - 1, nx)]);
    level->phi[ns_mg_idx(ny - 1, 0, nx)] = 0.5 * (level->phi[ns_mg_idx(ny - 1, 1, nx)] + level->phi[ns_mg_idx(ny - 2, 0, nx)]);
    level->phi[ns_mg_idx(ny - 1, nx - 1, nx)] = 0.5 * (level->phi[ns_mg_idx(ny - 1, nx - 2, nx)] + level->phi[ns_mg_idx(ny - 2, nx - 1, nx)]);
}

static void ns_mg_relax(NSMGLevel *level, int iterations) {
    if (!level || iterations <= 0) {
        return;
    }

    const int nx = level->nx;
    const int ny = level->ny;
    const double dx2 = level->dx * level->dx;
    const double dy2 = level->dy * level->dy;
    const double denom = 2.0 * (1.0 / dx2 + 1.0 / dy2);

    for (int iter = 0; iter < iterations; ++iter) {
        for (int color = 0; color < 2; ++color) {
#pragma omp parallel for
            for (int i = 1; i < ny - 1; ++i) {
                int j_start = 1 + ((i + color) & 1);
                for (int j = j_start; j < nx - 1; j += 2) {
                    const int idx = ns_mg_idx(i, j, nx);
                    const double rhs = level->rhs[idx];
                    const double phi_e = level->phi[ns_mg_idx(i, j + 1, nx)];
                    const double phi_w = level->phi[ns_mg_idx(i, j - 1, nx)];
                    const double phi_n = level->phi[ns_mg_idx(i + 1, j, nx)];
                    const double phi_s = level->phi[ns_mg_idx(i - 1, j, nx)];

                    const double sum = (phi_e + phi_w) / dx2 + (phi_n + phi_s) / dy2;
                    level->phi[idx] = (rhs + sum) / denom;
                }
            }
        }
    }
}

static void ns_mg_compute_residual(NSMGLevel *level) {
    if (!level) {
        return;
    }

    const int nx = level->nx;
    const int ny = level->ny;
    const double dx2 = level->dx * level->dx;
    const double dy2 = level->dy * level->dy;

#pragma omp parallel for
    for (int i = 1; i < ny - 1; ++i) {
        for (int j = 1; j < nx - 1; ++j) {
            const int idx = ns_mg_idx(i, j, nx);
            const double laplace =
                (level->phi[ns_mg_idx(i, j + 1, nx)] - 2.0 * level->phi[idx] + level->phi[ns_mg_idx(i, j - 1, nx)]) / dx2 +
                (level->phi[ns_mg_idx(i + 1, j, nx)] - 2.0 * level->phi[idx] + level->phi[ns_mg_idx(i - 1, j, nx)]) / dy2;
            level->residual[idx] = level->rhs[idx] - laplace;
        }
    }
}

static void ns_mg_restrict(const NSMGLevel *fine, NSMGLevel *coarse) {
    const int nx_f = fine->nx;
    const int nx_c = coarse->nx;
    const int ny_c = coarse->ny;

    memset(coarse->rhs, 0, (size_t)nx_c * (size_t)ny_c * sizeof(double));

    for (int I = 1; I < ny_c - 1; ++I) {
        for (int J = 1; J < nx_c - 1; ++J) {
            const int i = 2 * I;
            const int j = 2 * J;
            const int idx_f = ns_mg_idx(i, j, nx_f);
            const int idx_c = ns_mg_idx(I, J, nx_c);

            double sum = 0.25 * fine->residual[idx_f];
            sum += 0.125 * (fine->residual[ns_mg_idx(i, j - 1, nx_f)] +
                            fine->residual[ns_mg_idx(i, j + 1, nx_f)] +
                            fine->residual[ns_mg_idx(i - 1, j, nx_f)] +
                            fine->residual[ns_mg_idx(i + 1, j, nx_f)]);
            sum += 0.0625 * (fine->residual[ns_mg_idx(i - 1, j - 1, nx_f)] +
                             fine->residual[ns_mg_idx(i - 1, j + 1, nx_f)] +
                             fine->residual[ns_mg_idx(i + 1, j - 1, nx_f)] +
                             fine->residual[ns_mg_idx(i + 1, j + 1, nx_f)]);
            coarse->rhs[idx_c] = sum;
        }
    }
}

static void ns_mg_prolongate_add(const NSMGLevel *coarse, NSMGLevel *fine) {
    const int nx_f = fine->nx;
    const int ny_f = fine->ny;
    const int nx_c = coarse->nx;
    const int ny_c = coarse->ny;

    memset(fine->scratch, 0, (size_t)nx_f * (size_t)ny_f * sizeof(double));

    for (int i = 1; i < ny_f - 1; ++i) {
        const double fy = 0.5 * (double)(i % 2);
        const int I = i / 2;
        const int I1 = (I + 1 < ny_c) ? I + 1 : I;

        for (int j = 1; j < nx_f - 1; ++j) {
            const double fx = 0.5 * (double)(j % 2);
            const int J = j / 2;
            const int J1 = (J + 1 < nx_c) ? J + 1 : J;

            const double c00 = coarse->phi[ns_mg_idx(I, J, nx_c)];
            const double c10 = coarse->phi[ns_mg_idx(I, J1, nx_c)];
            const double c01 = coarse->phi[ns_mg_idx(I1, J, nx_c)];
            const double c11 = coarse->phi[ns_mg_idx(I1, J1, nx_c)];

            const double interp =
                (1.0 - fx) * (1.0 - fy) * c00 +
                fx * (1.0 - fy) * c10 +
                (1.0 - fx) * fy * c01 +
                fx * fy * c11;

            fine->scratch[ns_mg_idx(i, j, nx_f)] = interp;
        }
    }

    for (int i = 1; i < ny_f - 1; ++i) {
        for (int j = 1; j < nx_f - 1; ++j) {
            const int idx = ns_mg_idx(i, j, nx_f);
            fine->phi[idx] += fine->scratch[idx];
        }
    }
}

static double ns_mg_residual_norm(const NSMGLevel *level) {
    const int nx = level->nx;
    const int ny = level->ny;
    const int interior_x = nx - 2;
    const int interior_y = ny - 2;
    if (interior_x <= 0 || interior_y <= 0) {
        return 0.0;
    }

    double sum = 0.0;

#pragma omp parallel for reduction(+:sum)
    for (int i = 1; i < ny - 1; ++i) {
        for (int j = 1; j < nx - 1; ++j) {
            const double r = level->residual[ns_mg_idx(i, j, nx)];
            sum += r * r;
        }
    }

    const double norm = sqrt(sum / (double)(interior_x * interior_y));
    return norm;
}

static void ns_mg_zero_level(NSMGLevel *level) {
    if (!level) {
        return;
    }
    const size_t total = (size_t)level->nx * (size_t)level->ny;
    memset(level->phi, 0, total * sizeof(double));
}

static void ns_mg_cycle(NSMGHierarchy *hier,
                        const NSFiniteDiffData *data,
                        int level_idx,
                        const NSMultigridConfig *config) {
    NSMGLevel *level = &hier->levels[level_idx];
    ns_fd_multigrid_apply_pressure_bc(data, level);
    ns_mg_relax(level, config->pre_smoothing);
    ns_fd_multigrid_apply_pressure_bc(data, level);

    if (level_idx == hier->num_levels - 1) {
        ns_mg_relax(level, config->coarse_smoothing);
        ns_fd_multigrid_apply_pressure_bc(data, level);
        return;
    }

    ns_mg_compute_residual(level);

    NSMGLevel *coarse = &hier->levels[level_idx + 1];
    ns_mg_zero_level(coarse);
    memset(coarse->rhs, 0, (size_t)coarse->nx * (size_t)coarse->ny * sizeof(double));
    ns_mg_restrict(level, coarse);
    ns_fd_multigrid_apply_pressure_bc(data, coarse);

    const int repeats = (config->cycle_type == NS_MG_CYCLE_W) ? 2 : 1;
    for (int repeat = 0; repeat < repeats; ++repeat) {
        ns_mg_cycle(hier, data, level_idx + 1, config);
    }

    ns_mg_prolongate_add(coarse, level);
    ns_fd_multigrid_apply_pressure_bc(data, level);
    ns_mg_relax(level, config->post_smoothing);
    ns_fd_multigrid_apply_pressure_bc(data, level);
}

void ns_fd_multigrid_default_config(NSMultigridConfig *config) {
    if (!config) {
        return;
    }
    config->pre_smoothing = 2;
    config->post_smoothing = 2;
    config->coarse_smoothing = 20;
    config->max_cycles = 50;
    config->tolerance = 1e-6;
    config->cycle_type = NS_MG_CYCLE_V;
}

NSMultigridStats ns_fd_multigrid_solve_poisson(struct NSFiniteDiffData *data,
                                               const double *rhs,
                                               double *solution,
                                               const NSMultigridConfig *user_config,
                                               NSMultigridMonitor *monitor) {
    NSMultigridConfig local_config;
    if (user_config) {
        local_config = *user_config;
    } else {
        ns_fd_multigrid_default_config(&local_config);
    }

    NSMultigridStats stats = {0};
    if (!data || !solution || !rhs) {
        return stats;
    }

    if (monitor) {
        monitor->count = 0;
    }

    NSMGHierarchy hier = {0};
    if (!ns_mg_build_hierarchy(&hier, data->nx, data->ny, data->dx, data->dy)) {
        return stats;
    }

    NSMGLevel *fine = &hier.levels[0];
    const size_t total = (size_t)data->nx * (size_t)data->ny;
    memcpy(fine->phi, solution, total * sizeof(double));
    memcpy(fine->rhs, rhs, total * sizeof(double));

    ns_fd_multigrid_apply_pressure_bc(data, fine);
    ns_mg_compute_residual(fine);
    double residual = ns_mg_residual_norm(fine);
    stats.initial_residual = residual;
    stats.final_residual = residual;
    stats.cycles_performed = 0;
    stats.converged = (residual <= local_config.tolerance) ? 1 : 0;

    if (monitor && monitor->history && monitor->capacity > 0) {
        monitor->history[monitor->count++] = residual;
    }

    int cycle = 0;
    while (cycle < local_config.max_cycles && residual > local_config.tolerance) {
        ns_mg_cycle(&hier, data, 0, &local_config);
        ns_mg_compute_residual(fine);
        residual = ns_mg_residual_norm(fine);
        stats.final_residual = residual;
        ++cycle;
        stats.cycles_performed = cycle;

        if (monitor && monitor->history && monitor->count < monitor->capacity) {
            monitor->history[monitor->count++] = residual;
        }

        if (residual <= local_config.tolerance) {
            stats.converged = 1;
            break;
        }
    }

    memcpy(solution, fine->phi, total * sizeof(double));

    ns_mg_free_hierarchy(&hier);
    return stats;
}
