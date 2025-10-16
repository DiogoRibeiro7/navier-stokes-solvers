#include "ns_fd_solver.h"

#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

double ns_fd_compute_residual(NSFiniteDiffData *data) {
    const int nx = data->nx;
    const int ny = data->ny;
    const double dx = data->dx;
    const double dy = data->dy;
    const double dt = data->dt;
    const double inv_2dx = 0.5 / dx;
    const double inv_2dy = 0.5 / dy;
    const double inv_dx2 = 1.0 / (dx * dx);
    const double inv_dy2 = 1.0 / (dy * dy);
    const double inv_Re = 1.0 / data->Re;

    double residual = 0.0;

    double *NS_RESTRICT F = (double *NS_RESTRICT)NS_ASSUME_ALIGNED(data->F, 64);
    double *NS_RESTRICT G = (double *NS_RESTRICT)NS_ASSUME_ALIGNED(data->G, 64);
    const double *NS_RESTRICT u = (const double *NS_RESTRICT)NS_ASSUME_ALIGNED(data->u, 64);
    const double *NS_RESTRICT v = (const double *NS_RESTRICT)NS_ASSUME_ALIGNED(data->v, 64);
    const double *NS_RESTRICT p = (const double *NS_RESTRICT)NS_ASSUME_ALIGNED(data->p, 64);
    const double *NS_RESTRICT u_old = (const double *NS_RESTRICT)NS_ASSUME_ALIGNED(data->u_old, 64);
    const double *NS_RESTRICT v_old = (const double *NS_RESTRICT)NS_ASSUME_ALIGNED(data->v_old, 64);

#pragma omp parallel for collapse(2) reduction(+:residual)
    for (int i = 1; i < ny - 1; ++i) {
        for (int j = 1; j < nx - 1; ++j) {
            const int ij = IDX2D(i, j, nx);

            const double uij = u[ij];
            const double vij = v[ij];

            const double du_dx = (u[IDX2D(i, j + 1, nx)] - u[IDX2D(i, j - 1, nx)]) * inv_2dx;
            const double du_dy = (u[IDX2D(i + 1, j, nx)] - u[IDX2D(i - 1, j, nx)]) * inv_2dy;
            const double dv_dx = (v[IDX2D(i, j + 1, nx)] - v[IDX2D(i, j - 1, nx)]) * inv_2dx;
            const double dv_dy = (v[IDX2D(i + 1, j, nx)] - v[IDX2D(i - 1, j, nx)]) * inv_2dy;

            const double convection_u = uij * du_dx + vij * du_dy;
            const double convection_v = uij * dv_dx + vij * dv_dy;

            const double u_lap = (u[IDX2D(i, j + 1, nx)] - 2.0 * uij + u[IDX2D(i, j - 1, nx)]) * inv_dx2 +
                                 (u[IDX2D(i + 1, j, nx)] - 2.0 * uij + u[IDX2D(i - 1, j, nx)]) * inv_dy2;

            const double v_lap = (v[IDX2D(i, j + 1, nx)] - 2.0 * vij + v[IDX2D(i, j - 1, nx)]) * inv_dx2 +
                                 (v[IDX2D(i + 1, j, nx)] - 2.0 * vij + v[IDX2D(i - 1, j, nx)]) * inv_dy2;

            const double dp_dx = (p[IDX2D(i, j + 1, nx)] - p[IDX2D(i, j - 1, nx)]) * inv_2dx;
            const double dp_dy = (p[IDX2D(i + 1, j, nx)] - p[IDX2D(i - 1, j, nx)]) * inv_2dy;

            const double R_u = (uij - u_old[ij]) / dt + convection_u + dp_dx - inv_Re * u_lap;
            const double R_v = (vij - v_old[ij]) / dt + convection_v + dp_dy - inv_Re * v_lap;

            F[2 * ij] = R_u;
            F[2 * ij + 1] = R_v;

            const double divergence = du_dx + dv_dy;
            G[ij] = divergence;

            residual += R_u * R_u + R_v * R_v + divergence * divergence;
        }
    }

    return sqrt(residual);
}

void ns_fd_assemble_jacobian(NSFiniteDiffData *data) {
    const int nx = data->nx;
    const int ny = data->ny;
    const double dx = data->dx;
    const double dy = data->dy;
    const double dt = data->dt;
    const double Re = data->Re;
    const double diag = 1.0 / dt + 2.0 / (Re * dx * dx) + 2.0 / (Re * dy * dy);

#pragma omp parallel for collapse(2)
    for (int i = 1; i < ny - 1; ++i) {
        for (int j = 1; j < nx - 1; ++j) {
            const int ij = IDX2D(i, j, nx);
            data->J[ij] = diag;
        }
    }
}

void ns_fd_solve_linear_system(NSFiniteDiffData *data, double *delta) {
    const int nx = data->nx;
    const int ny = data->ny;
    const int total_points = nx * ny;
    const double inv_diag = 1.0 / (1.0 / data->dt + 1e-6);

#pragma omp parallel for simd
    for (int i = 0; i < 2 * total_points; ++i) {
        delta[i] = -data->F[i] * inv_diag;
    }

#pragma omp parallel for simd
    for (int i = 0; i < total_points; ++i) {
        delta[2 * total_points + i] = -data->G[i];
    }
}

SolverStatus ns_fd_newton_raphson_step(NSFiniteDiffData *data) {
    const int nx = data->nx;
    const int ny = data->ny;
    const int total_points = nx * ny;
    const int n = 3 * total_points;
    double *NS_RESTRICT delta = (double *NS_RESTRICT)NS_ASSUME_ALIGNED(data->delta, 64);

    if (!delta) {
        return SOLVER_MEM_ERROR;
    }

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        delta[i] = 0.0;
    }

    const double residual_norm = ns_fd_compute_residual(data);
    if (residual_norm < TOL) {
        return SOLVER_SUCCESS;
    }

    ns_fd_assemble_jacobian(data);
    ns_fd_solve_linear_system(data, delta);

#pragma omp parallel for collapse(2)
    for (int i = 1; i < ny - 1; ++i) {
        for (int j = 1; j < nx - 1; ++j) {
            const int ij = IDX2D(i, j, nx);
            data->u[ij] += delta[2 * ij];
            data->v[ij] += delta[2 * ij + 1];
            data->p[ij] += delta[2 * total_points + ij];
        }
    }

    ns_fd_apply_boundary_conditions(data);

    if (residual_norm > 1e6) {
        return SOLVER_DIVERGED;
    }
    return SOLVER_MAX_ITER;
}
