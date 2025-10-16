#include "ns_fd_solver.h"
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// Compute adaptive timestep
double ns_fd_adaptive_timestep(NSFiniteDiffData *data) {
    int nx = data->nx, ny = data->ny;
    double dx = data->dx, dy = data->dy;
    double max_u = 0.0, max_v = 0.0;
    
#pragma omp parallel for collapse(2) reduction(max:max_u, max_v)
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++) {
            int ij = IDX2D(i, j, nx);
            max_u = fmax(max_u, fabs(data->u[ij]));
            max_v = fmax(max_v, fabs(data->v[ij]));
        }
    }
    
    // CFL condition
    double dt_cfl = 0.5 / (max_u/dx + max_v/dy + 1e-10);
    
    // Viscous stability
    double dt_visc = 0.5 * data->Re * dx * dy / (dx + dy);
    
    double dt_new = fmin(dt_cfl, dt_visc);
    dt_new = fmax(dt_new, MIN_DT);
    dt_new = fmin(dt_new, MAX_DT);
    
    return dt_new;
}

// Advance one timestep with Newton-Raphson
SolverStatus ns_fd_advance_timestep(NSFiniteDiffData *data) {
    int nx = data->nx, ny = data->ny;
    
    // Store old values
    memcpy(data->u_old, data->u, nx * ny * sizeof(double));
    memcpy(data->v_old, data->v, nx * ny * sizeof(double));
    memcpy(data->p_old, data->p, nx * ny * sizeof(double));
    
    // Adaptive timestep
    data->dt = ns_fd_adaptive_timestep(data);
    
    // Newton-Raphson iteration
    SolverStatus status = SOLVER_MAX_ITER;
    for (int iter = 0; iter < MAX_ITER; iter++) {
        status = ns_fd_newton_raphson_step(data);
        
        if (status == SOLVER_SUCCESS) break;
        
        if (iter % 10 == 0) {
            double residual = ns_fd_compute_residual(data);
            printf("  NR iter %d: residual = %.2e\n", iter, residual);
        }
        
        if (status == SOLVER_DIVERGED) {
            printf("Newton-Raphson diverged!\n");
            return status;
        }
    }
    
    data->t += data->dt;
    return status;
}
