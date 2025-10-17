#include "ns_fd_solver.h"
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static void ns_fd_apply_periodic_y(NSFiniteDiffData *data) {
    const int nx = data->nx;
    const int ny = data->ny;
    if (ny < 2) {
        return;
    }

    const int lower_src = (ny >= 2) ? ny - 2 : 0;
    const int upper_src = (ny >= 2) ? 1 : 0;

#pragma omp parallel for
    for (int j = 0; j < nx; ++j) {
        const int bottom = IDX2D(0, j, nx);
        const int top = IDX2D(ny - 1, j, nx);
        const int bottom_src = IDX2D(lower_src, j, nx);
        const int top_src = IDX2D(upper_src, j, nx);

        data->u[bottom] = data->u[bottom_src];
        data->v[bottom] = data->v[bottom_src];
        data->p[bottom] = data->p[bottom_src];

        data->u[top] = data->u[top_src];
        data->v[top] = data->v[top_src];
        data->p[top] = data->p[top_src];
    }
}

static void ns_fd_apply_periodic_x(NSFiniteDiffData *data) {
    const int nx = data->nx;
    const int ny = data->ny;
    if (nx < 2) {
        return;
    }

    const int left_src = (nx >= 2) ? nx - 2 : 0;
    const int right_src = (nx >= 2) ? 1 : 0;

#pragma omp parallel for
    for (int i = 0; i < ny; ++i) {
        const int left = IDX2D(i, 0, nx);
        const int right = IDX2D(i, nx - 1, nx);
        const int left_neighbor = IDX2D(i, left_src, nx);
        const int right_neighbor = IDX2D(i, right_src, nx);

        data->u[left] = data->u[left_neighbor];
        data->v[left] = data->v[left_neighbor];
        data->p[left] = data->p[left_neighbor];

        data->u[right] = data->u[right_neighbor];
        data->v[right] = data->v[right_neighbor];
        data->p[right] = data->p[right_neighbor];
    }
}

static void ns_fd_apply_horizontal_boundary(NSFiniteDiffData *data, int side) {
    const NSBoundaryCondition *bc = &data->boundaries[side];
    const int nx = data->nx;
    const int ny = data->ny;

    if (nx == 0 || ny == 0) {
        return;
    }

    const int row = (side == NS_BC_BOTTOM) ? 0 : ny - 1;
    const int neighbor_row = (side == NS_BC_BOTTOM)
                                 ? ((ny >= 2) ? 1 : 0)
                                 : ((ny >= 2) ? ny - 2 : 0);

#pragma omp parallel for
    for (int j = 0; j < nx; ++j) {
        const double x = j * data->dx;
        const double y = row * data->dy;
        const int idx = IDX2D(row, j, nx);
        const int neighbor = IDX2D(neighbor_row, j, nx);
        int pressure_set = 0;

        switch (bc->type) {
            case BC_NO_SLIP:
            case BC_DIRICHLET:
                data->u[idx] = bc->value_u;
                data->v[idx] = bc->value_v;
                if (bc->type == BC_DIRICHLET || isfinite(bc->value_p)) {
                    data->p[idx] = bc->value_p;
                    pressure_set = 1;
                }
                break;
            case BC_FREE_SLIP:
                data->v[idx] = 0.0;
                data->u[idx] = data->u[neighbor];
                break;
            case BC_INFLOW:
                data->u[idx] = bc->value_u;
                data->v[idx] = bc->value_v;
                if (isfinite(bc->value_p)) {
                    data->p[idx] = bc->value_p;
                    pressure_set = 1;
                }
                break;
            case BC_OUTFLOW:
            case BC_NEUMANN:
            case BC_PERIODIC:
                data->u[idx] = data->u[neighbor];
                data->v[idx] = data->v[neighbor];
                if (isfinite(bc->value_p)) {
                    data->p[idx] = bc->value_p;
                    pressure_set = 1;
                }
                break;
            case BC_CUSTOM:
                if (bc->func) {
                    double u_val = data->u[idx];
                    double v_val = data->v[idx];
                    double p_val = data->p[idx];
                    bc->func(x, y, data->t, &u_val, &v_val, &p_val, bc->user_data);
                    data->u[idx] = u_val;
                    data->v[idx] = v_val;
                    data->p[idx] = p_val;
                    pressure_set = 1;
                } else {
                    data->u[idx] = bc->value_u;
                    data->v[idx] = bc->value_v;
                    if (isfinite(bc->value_p)) {
                        data->p[idx] = bc->value_p;
                        pressure_set = 1;
                    }
                }
                break;
        }

        if (!pressure_set && bc->type != BC_CUSTOM) {
            data->p[idx] = data->p[neighbor];
        }
    }
}

static void ns_fd_apply_vertical_boundary(NSFiniteDiffData *data, int side) {
    const NSBoundaryCondition *bc = &data->boundaries[side];
    const int nx = data->nx;
    const int ny = data->ny;

    if (nx == 0 || ny == 0) {
        return;
    }

    const int column = (side == NS_BC_LEFT) ? 0 : nx - 1;
    const int neighbor_col = (side == NS_BC_LEFT)
                                 ? ((nx >= 2) ? 1 : 0)
                                 : ((nx >= 2) ? nx - 2 : 0);

#pragma omp parallel for
    for (int i = 0; i < ny; ++i) {
        const double x = column * data->dx;
        const double y = i * data->dy;
        const int idx = IDX2D(i, column, nx);
        const int neighbor = IDX2D(i, neighbor_col, nx);
        int pressure_set = 0;

        switch (bc->type) {
            case BC_NO_SLIP:
            case BC_DIRICHLET:
                data->u[idx] = bc->value_u;
                data->v[idx] = bc->value_v;
                if (bc->type == BC_DIRICHLET || isfinite(bc->value_p)) {
                    data->p[idx] = bc->value_p;
                    pressure_set = 1;
                }
                break;
            case BC_FREE_SLIP:
                data->u[idx] = 0.0;
                data->v[idx] = data->v[neighbor];
                break;
            case BC_INFLOW:
                data->u[idx] = bc->value_u;
                data->v[idx] = bc->value_v;
                if (isfinite(bc->value_p)) {
                    data->p[idx] = bc->value_p;
                    pressure_set = 1;
                }
                break;
            case BC_OUTFLOW:
            case BC_NEUMANN:
            case BC_PERIODIC:
                data->u[idx] = data->u[neighbor];
                data->v[idx] = data->v[neighbor];
                if (isfinite(bc->value_p)) {
                    data->p[idx] = bc->value_p;
                    pressure_set = 1;
                }
                break;
            case BC_CUSTOM:
                if (bc->func) {
                    double u_val = data->u[idx];
                    double v_val = data->v[idx];
                    double p_val = data->p[idx];
                    bc->func(x, y, data->t, &u_val, &v_val, &p_val, bc->user_data);
                    data->u[idx] = u_val;
                    data->v[idx] = v_val;
                    data->p[idx] = p_val;
                    pressure_set = 1;
                } else {
                    data->u[idx] = bc->value_u;
                    data->v[idx] = bc->value_v;
                    if (isfinite(bc->value_p)) {
                        data->p[idx] = bc->value_p;
                        pressure_set = 1;
                    }
                }
                break;
        }

        if (!pressure_set && bc->type != BC_CUSTOM) {
            data->p[idx] = data->p[neighbor];
        }
    }
}

// Initialize lid-driven cavity flow
void ns_fd_initialize_lid_cavity(NSFiniteDiffData *data) {
    const int nx = data->nx;
    const int ny = data->ny;

    memset(data->u, 0, (size_t)nx * ny * sizeof(double));
    memset(data->v, 0, (size_t)nx * ny * sizeof(double));
    memset(data->p, 0, (size_t)nx * ny * sizeof(double));

    const double lid_u = data->boundaries[NS_BC_TOP].value_u;
    const double lid_v = data->boundaries[NS_BC_TOP].value_v;

#pragma omp parallel for
    for (int j = 0; j < nx; ++j) {
        const int idx = IDX2D(ny - 1, j, nx);
        data->u[idx] = lid_u;
        data->v[idx] = lid_v;
    }

    ns_fd_apply_boundary_conditions(data);
}

// Initialize Taylor-Green vortex
void ns_fd_initialize_taylor_green(NSFiniteDiffData *data) {
    const int nx = data->nx;
    const int ny = data->ny;
    const double L = data->L;
    const double H = data->H;

#pragma omp parallel for
    for (int i = 0; i < ny; ++i) {
#pragma omp simd
        for (int j = 0; j < nx; ++j) {
            const double x = j * data->dx;
            const double y = i * data->dy;
            const int idx = IDX2D(i, j, nx);

            data->u[idx] = sin(2.0 * M_PI * x / L) * cos(2.0 * M_PI * y / H);
            data->v[idx] = -cos(2.0 * M_PI * x / L) * sin(2.0 * M_PI * y / H);
            data->p[idx] = -0.25 * (cos(4.0 * M_PI * x / L) + cos(4.0 * M_PI * y / H));
        }
    }
}

// Apply boundary conditions for the finite-difference solver
void ns_fd_apply_boundary_conditions(NSFiniteDiffData *data) {
    if (!data) {
        return;
    }

    const int periodic_y = (data->boundaries[NS_BC_BOTTOM].type == BC_PERIODIC &&
                            data->boundaries[NS_BC_TOP].type == BC_PERIODIC);
    const int periodic_x = (data->boundaries[NS_BC_LEFT].type == BC_PERIODIC &&
                            data->boundaries[NS_BC_RIGHT].type == BC_PERIODIC);

    if (periodic_y) {
        ns_fd_apply_periodic_y(data);
    } else {
        ns_fd_apply_horizontal_boundary(data, NS_BC_BOTTOM);
        ns_fd_apply_horizontal_boundary(data, NS_BC_TOP);
    }

    if (periodic_x) {
        ns_fd_apply_periodic_x(data);
    } else {
        ns_fd_apply_vertical_boundary(data, NS_BC_LEFT);
        ns_fd_apply_vertical_boundary(data, NS_BC_RIGHT);
    }
}
