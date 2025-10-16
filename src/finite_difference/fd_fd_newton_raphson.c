#include "ns_fd_solver.h"

// Compute convective terms
static void compute_convection(NSFiniteDiffData *data, double *uu, double *uv) {
    int nx = data->nx, ny = data->ny;
    double dx = data->dx, dy = data->dy;
    
    for (int i = 1; i < ny - 1; i++) {
        for (int j = 1; j < nx - 1; j++) {
            int ij = IDX2D(i, j, nx);
            
            // u * du/dx + v * du/dy
            uu[ij] = data->u[ij] * (data->u[IDX2D(i, j+1, nx)] - data->u[IDX2D(i, j-1, nx)]) / (2*dx) +
                     data->v[ij] * (data->u[IDX2D(i+1, j, nx)] - data->u[IDX2D(i-1, j, nx)]) / (2*dy);
            
            // u * dv/dx + v * dv/dy
            uv[ij] = data->u[ij] * (data->v[IDX2D(i, j+1, nx)] - data->v[IDX2D(i, j-1, nx)]) / (2*dx) +
                     data->v[ij] * (data->v[IDX2D(i+1, j, nx)] - data->v[IDX2D(i-1, j, nx)]) / (2*dy);
        }
    }
}

// Compute residual for Newton-Raphson
double ns_fd_compute_residual(NSFiniteDiffData *data) {
    int nx = data->nx, ny = data->ny;
    double dx = data->dx, dy = data->dy;
    double dt = data->dt, Re = data->Re;
    double residual = 0.0;
    
    double *uu = calloc(nx * ny, sizeof(double));
    double *uv = calloc(nx * ny, sizeof(double));
    
    compute_convection(data, uu, uv);
    
    // Momentum equations residual
    for (int i = 1; i < ny - 1; i++) {
        for (int j = 1; j < nx - 1; j++) {
            int ij = IDX2D(i, j, nx);
            
            // u-momentum
            double u_xx = (data->u[IDX2D(i, j+1, nx)] - 2*data->u[ij] + data->u[IDX2D(i, j-1, nx)]) / (dx*dx);
            double u_yy = (data->u[IDX2D(i+1, j, nx)] - 2*data->u[ij] + data->u[IDX2D(i-1, j, nx)]) / (dy*dy);
            double dp_dx = (data->p[IDX2D(i, j+1, nx)] - data->p[IDX2D(i, j-1, nx)]) / (2*dx);
            
            double R_u = (data->u[ij] - data->u_old[ij])/dt + uu[ij] + dp_dx - (u_xx + u_yy)/Re;
            data->F[2*ij] = R_u;
            residual += R_u * R_u;
            
            // v-momentum
            double v_xx = (data->v[IDX2D(i, j+1, nx)] - 2*data->v[ij] + data->v[IDX2D(i, j-1, nx)]) / (dx*dx);
            double v_yy = (data->v[IDX2D(i+1, j, nx)] - 2*data->v[ij] + data->v[IDX2D(i-1, j, nx)]) / (dy*dy);
            double dp_dy = (data->p[IDX2D(i+1, j, nx)] - data->p[IDX2D(i-1, j, nx)]) / (2*dy);
            
            double R_v = (data->v[ij] - data->v_old[ij])/dt + uv[ij] + dp_dy - (v_xx + v_yy)/Re;
            data->F[2*ij + 1] = R_v;
            residual += R_v * R_v;
        }
    }
    
    // Continuity equation residual
    for (int i = 1; i < ny - 1; i++) {
        for (int j = 1; j < nx - 1; j++) {
            int ij = IDX2D(i, j, nx);
            double du_dx = (data->u[IDX2D(i, j+1, nx)] - data->u[IDX2D(i, j-1, nx)]) / (2*dx);
            double dv_dy = (data->v[IDX2D(i+1, j, nx)] - data->v[IDX2D(i-1, j, nx)]) / (2*dy);
            
            double R_p = du_dx + dv_dy;
            data->G[ij] = R_p;
            residual += R_p * R_p;
        }
    }
    
    free(uu); free(uv);
    return sqrt(residual);
}

// Simplified Jacobian assembly (diagonal approximation)
void ns_fd_assemble_jacobian(NSFiniteDiffData *data) {
    int nx = data->nx, ny = data->ny;
    double dx = data->dx, dy = data->dy;
    double dt = data->dt, Re = data->Re;
    
    // Diagonal dominant approximation for efficiency
    for (int i = 1; i < ny - 1; i++) {
        for (int j = 1; j < nx - 1; j++) {
            int ij = IDX2D(i, j, nx);
            data->J[ij] = 1.0/dt + 2.0/(Re*dx*dx) + 2.0/(Re*dy*dy);
        }
    }
}

// Solve linear system (simplified diagonal preconditioned)
void ns_fd_solve_linear_system(NSFiniteDiffData *data, double *delta) {
    int nx = data->nx, ny = data->ny;
    int n = 2 * nx * ny + nx * ny;
    
    for (int i = 0; i < n; i++) {
        if (i < 2 * nx * ny) {
            delta[i] = -data->F[i] / (1.0/data->dt + 1e-6);
        } else {
            int p_idx = i - 2 * nx * ny;
            delta[i] = -data->G[p_idx];
        }
        delta[i] *= 0.1; // Damping
    }
}

// Newton-Raphson iteration step
SolverStatus ns_fd_newton_raphson_step(NSFiniteDiffData *data) {
    int nx = data->nx, ny = data->ny;
    int n = 2 * nx * ny + nx * ny;
    double *delta = calloc(n, sizeof(double));
    
    if (!delta) return SOLVER_MEM_ERROR;
    
    double residual_norm = ns_fd_compute_residual(data);
    
    if (residual_norm < TOL) {
        free(delta);
        return SOLVER_SUCCESS;
    }
    
    ns_fd_assemble_jacobian(data);
    ns_fd_solve_linear_system(data, delta);
    
    // Update solution
    for (int i = 1; i < ny - 1; i++) {
        for (int j = 1; j < nx - 1; j++) {
            int ij = IDX2D(i, j, nx);
            data->u[ij] += delta[2*ij];
            data->v[ij] += delta[2*ij + 1];
            data->p[ij] += delta[2*nx*ny + ij];
        }
    }
    
    ns_fd_apply_boundary_conditions(data);
    
    free(delta);
    return (residual_norm > 1e6) ? SOLVER_DIVERGED : SOLVER_MAX_ITER;
}
