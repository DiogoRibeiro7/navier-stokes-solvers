// Lid-driven cavity flow at Re = 100
// Expected results: Steady single primary vortex with velocity profiles matching
//                   Ghia et al. (1982) benchmark data.
// References: Ghia, U., Ghia, K. N., & Shin, C. T. (1982). High-Re solutions for
//             incompressible flow using the Navier-Stokes equations and a
//             multigrid method. J. Comput. Phys., 48(3), 387-411.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <errno.h>

#ifdef _WIN32
#include <direct.h>
#define mkdir_p(path) _mkdir(path)
#else
#include <sys/stat.h>
#define mkdir_p(path) mkdir(path, 0777)
#endif

#include "../src/finite_difference/ns_fd_solver.h"

typedef struct {
    double position;
    double value;
} LineSample;

static const LineSample GHIA_RE100_VERTICAL[] = {
    {1.0000, 1.00000}, {0.9766, 0.84123}, {0.9688, 0.78871}, {0.9609, 0.73722},
    {0.9531, 0.68717}, {0.9453, 0.63866}, {0.9063, 0.47244}, {0.8594, 0.34682},
    {0.8047, 0.24299}, {0.5000, 0.00000}, {0.2344, -0.22445}, {0.2266, -0.23366},
    {0.1563, -0.28536}, {0.0938, -0.32066}, {0.0781, -0.32407}, {0.0703, -0.32726},
    {0.0625, -0.32995}, {0.0547, -0.33190}, {0.0469, -0.33257}, {0.0000, 0.00000}
};

static const LineSample GHIA_RE100_HORIZONTAL[] = {
    {1.0000, 0.00000}, {0.9688, -0.13641}, {0.9609, -0.14612}, {0.9531, -0.15660},
    {0.9453, -0.16730}, {0.9063, -0.20719}, {0.8594, -0.22847}, {0.8047, -0.23827},
    {0.5000, 0.00000}, {0.2344, 0.22445}, {0.2266, 0.23411}, {0.1563, 0.27669},
    {0.0938, 0.29093}, {0.0781, 0.29592}, {0.0703, 0.29916}, {0.0625, 0.30174},
    {0.0547, 0.30337}, {0.0469, 0.30413}, {0.0000, 0.00000}
};

static int ensure_directory(const char *path) {
    if (mkdir_p(path) == 0 || errno == EEXIST) {
        return 0;
    }
    fprintf(stderr, "Failed to create directory '%s': %s\n", path, strerror(errno));
    return -1;
}

static double interpolate_centerline_u(const NSFiniteDiffData *data, double y) {
    const int nx = data->nx;
    const int ny = data->ny;
    const int mid = nx / 2;
    const double dy = data->dy;
    double idxf = y / dy;
    if (idxf <= 0.0) {
        return data->u[IDX2D(0, mid, nx)];
    }
    if (idxf >= ny - 1) {
        return data->u[IDX2D(ny - 1, mid, nx)];
    }
    int i0 = (int)floor(idxf);
    double alpha = idxf - i0;
    double u0 = data->u[IDX2D(i0, mid, nx)];
    double u1 = data->u[IDX2D(i0 + 1, mid, nx)];
    return (1.0 - alpha) * u0 + alpha * u1;
}

static double interpolate_centerline_v(const NSFiniteDiffData *data, double x) {
    const int nx = data->nx;
    const int ny = data->ny;
    const int mid = ny / 2;
    const double dx = data->dx;
    double idxf = x / dx;
    if (idxf <= 0.0) {
        return data->v[IDX2D(mid, 0, nx)];
    }
    if (idxf >= nx - 1) {
        return data->v[IDX2D(mid, nx - 1, nx)];
    }
    int j0 = (int)floor(idxf);
    double alpha = idxf - j0;
    double v0 = data->v[IDX2D(mid, j0, nx)];
    double v1 = data->v[IDX2D(mid, j0 + 1, nx)];
    return (1.0 - alpha) * v0 + alpha * v1;
}

static void write_centerline_profiles(const NSFiniteDiffData *data, const char *filepath) {
    FILE *fp = fopen(filepath, "w");
    if (!fp) {
        fprintf(stderr, "Failed to write profiles '%s': %s\n", filepath, strerror(errno));
        return;
    }

    const int nx = data->nx;
    const int ny = data->ny;
    const int mid_x = nx / 2;
    const int mid_y = ny / 2;

    fprintf(fp, "# Vertical centreline (x = 0.5)\n");
    fprintf(fp, "# y u\n");
    for (int i = 0; i < ny; ++i) {
        double y = i * data->dy;
        double u_mid = data->u[IDX2D(i, mid_x, nx)];
        fprintf(fp, "%.6f %.8e\n", y, u_mid);
    }

    fprintf(fp, "\n# Horizontal centreline (y = 0.5)\n");
    fprintf(fp, "# x v\n");
    for (int j = 0; j < nx; ++j) {
        double x = j * data->dx;
        double v_mid = data->v[IDX2D(mid_y, j, nx)];
        fprintf(fp, "%.6f %.8e\n", x, v_mid);
    }

    fclose(fp);
}

static void evaluate_against_ghia(const NSFiniteDiffData *data) {
    const size_t vert_count = sizeof(GHIA_RE100_VERTICAL) / sizeof(GHIA_RE100_VERTICAL[0]);
    const size_t horiz_count = sizeof(GHIA_RE100_HORIZONTAL) / sizeof(GHIA_RE100_HORIZONTAL[0]);
    double err_vert_l2 = 0.0, err_vert_linf = 0.0;
    double err_horiz_l2 = 0.0, err_horiz_linf = 0.0;

    for (size_t i = 0; i < vert_count; ++i) {
        double sim = interpolate_centerline_u(data, GHIA_RE100_VERTICAL[i].position);
        double diff = sim - GHIA_RE100_VERTICAL[i].value;
        err_vert_l2 += diff * diff;
        err_vert_linf = fmax(err_vert_linf, fabs(diff));
    }
    err_vert_l2 = sqrt(err_vert_l2 / vert_count);

    for (size_t i = 0; i < horiz_count; ++i) {
        double sim = interpolate_centerline_v(data, GHIA_RE100_HORIZONTAL[i].position);
        double diff = sim - GHIA_RE100_HORIZONTAL[i].value;
        err_horiz_l2 += diff * diff;
        err_horiz_linf = fmax(err_horiz_linf, fabs(diff));
    }
    err_horiz_l2 = sqrt(err_horiz_l2 / horiz_count);

    printf("\nValidation against Ghia et al. (1982) Re=100 data:\n");
    printf("  Vertical centerline  L2 error = %.4e, Linf error = %.4e\n",
           err_vert_l2, err_vert_linf);
    printf("  Horizontal centerline L2 error = %.4e, Linf error = %.4e\n",
           err_horiz_l2, err_horiz_linf);
}

static double compute_velocity_change(const NSFiniteDiffData *data) {
    const int total = data->nx * data->ny;
    double max_diff = 0.0;
    for (int k = 0; k < total; ++k) {
        double du = fabs(data->u[k] - data->u_old[k]);
        double dv = fabs(data->v[k] - data->v_old[k]);
        max_diff = fmax(max_diff, du);
        max_diff = fmax(max_diff, dv);
    }
    return max_diff;
}

int main(void) {
    const int nx = 96;
    const int ny = 96;
    const double L = 1.0;
    const double H = 1.0;
    const double Re = 100.0;
    const double steady_residual_tol = 5e-7;
    const double steady_change_tol = 5e-7;
    const int max_steps = 12000;
    const int report_interval = 50;

    printf("============================================================\n");
    printf(" Lid-Driven Cavity Flow (Re = 100)\n");
    printf(" Finite Difference Solver Demo (ny=nx=%d)\n", nx);
    printf("============================================================\n\n");

    NSFiniteDiffData *solver = ns_fd_allocate(nx, ny, L, H, Re);
    if (!solver) {
        fprintf(stderr, "Failed to allocate solver memory.\n");
        return EXIT_FAILURE;
    }

    solver->bc_type = BC_NO_SLIP;
    solver->bc_values[0] = 0.0;
    solver->bc_values[1] = 1.0;  // moving lid
    solver->bc_values[2] = 0.0;
    solver->bc_values[3] = 0.0;

    ns_fd_initialize_lid_cavity(solver);

    const char *output_dir = "examples/results";
    if (ensure_directory("examples") || ensure_directory(output_dir)) {
        ns_fd_free(solver);
        return EXIT_FAILURE;
    }

    clock_t start_clock = clock();
    SolverStats stats = {0};
    int step;
    for (step = 1; step <= max_steps; ++step) {
        SolverStatus status = ns_fd_advance_timestep(solver);
        if (status != SOLVER_SUCCESS) {
            fprintf(stderr, "Solver step failed (code %d) at step %d.\n", status, step);
            break;
        }

        ns_fd_apply_boundary_conditions(solver);
        ns_fd_compute_stats(solver, &stats);

        double max_change = compute_velocity_change(solver);

        if (step % report_interval == 0) {
            printf("Step %5d | t = %.4f | dt = %.3e | res = %.3e | max|Δu| = %.3e\n",
                   step, solver->t, solver->dt, stats.residual, max_change);
        }

        if (stats.residual < steady_residual_tol && max_change < steady_change_tol) {
            printf("\nSteady state reached at step %d (residual %.3e, max|Δu| %.3e).\n",
                   step, stats.residual, max_change);
            break;
        }
    }

    clock_t end_clock = clock();
    double wall_time = (double)(end_clock - start_clock) / CLOCKS_PER_SEC;

    char solution_path[256];
    snprintf(solution_path, sizeof(solution_path),
             "%s/lid_driven_cavity_re100.dat", output_dir);
    ns_fd_output_solution(solver, solution_path);

    char profile_path[256];
    snprintf(profile_path, sizeof(profile_path),
             "%s/lid_driven_cavity_re100_profiles.dat", output_dir);
    write_centerline_profiles(solver, profile_path);

    evaluate_against_ghia(solver);

    printf("\nSummary:\n");
    printf("  Total steps: %d\n", step);
    printf("  Final time:  %.4f\n", solver->t);
    printf("  Final residual: %.3e\n", stats.residual);
    printf("  CPU time: %.2f s\n", wall_time);
    printf("  Results:\n");
    printf("    %s\n", solution_path);
    printf("    %s\n", profile_path);

    ns_fd_free(solver);
    return EXIT_SUCCESS;
}
