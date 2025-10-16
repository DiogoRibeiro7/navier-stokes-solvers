#include "ns_fd_solver.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// Initialize lid-driven cavity flow
void ns_fd_initialize_lid_cavity(NSFiniteDiffData *data) {
    int nx = data->nx, ny = data->ny;
    
    // Zero initial conditions
    memset(data->u, 0, nx * ny * sizeof(double));
    memset(data->v, 0, nx * ny * sizeof(double));
    memset(data->p, 0, nx * ny * sizeof(double));
    
    // Lid velocity at top boundary
#pragma omp parallel for
    for (int j = 0; j < nx; j++) {
        data->u[IDX2D(ny-1, j, nx)] = data->bc_values[1];
    }
    
    ns_fd_apply_boundary_conditions(data);
}

// Initialize Taylor-Green vortex
void ns_fd_initialize_taylor_green(NSFiniteDiffData *data) {
    int nx = data->nx, ny = data->ny;
    double L = data->L, H = data->H;
    
#pragma omp parallel for
    for (int i = 0; i < ny; i++) {
#pragma omp simd
        for (int j = 0; j < nx; j++) {
            const double x = j * data->dx;
            const double y = i * data->dy;
            const int idx = IDX2D(i, j, nx);

            data->u[idx] = sin(2 * M_PI * x / L) * cos(2 * M_PI * y / H);
            data->v[idx] = -cos(2 * M_PI * x / L) * sin(2 * M_PI * y / H);
            data->p[idx] = -0.25 * (cos(4 * M_PI * x / L) + cos(4 * M_PI * y / H));
        }
    }
}

// Apply boundary conditions (no-slip for lid-driven cavity)
void ns_fd_apply_boundary_conditions(NSFiniteDiffData *data) {
    int nx = data->nx, ny = data->ny;
    
    if (data->bc_type == BC_NO_SLIP) {
        // Bottom wall
#pragma omp parallel for
        for (int j = 0; j < nx; j++) {
            data->u[IDX2D(0, j, nx)] = data->bc_values[0];
            data->v[IDX2D(0, j, nx)] = 0.0;
        }
        
        // Top wall (lid)
#pragma omp parallel for
        for (int j = 0; j < nx; j++) {
            data->u[IDX2D(ny-1, j, nx)] = data->bc_values[1];
            data->v[IDX2D(ny-1, j, nx)] = 0.0;
        }
        
        // Left wall
#pragma omp parallel for
        for (int i = 0; i < ny; i++) {
            data->u[IDX2D(i, 0, nx)] = data->bc_values[2];
            data->v[IDX2D(i, 0, nx)] = 0.0;
        }
        
        // Right wall
#pragma omp parallel for
        for (int i = 0; i < ny; i++) {
            data->u[IDX2D(i, nx-1, nx)] = data->bc_values[3];
            data->v[IDX2D(i, nx-1, nx)] = 0.0;
        }
    }
    
    // Pressure: Neumann BC (dp/dn = 0)
#pragma omp parallel for
    for (int j = 1; j < nx-1; j++) {
        data->p[IDX2D(0, j, nx)] = data->p[IDX2D(1, j, nx)];
        data->p[IDX2D(ny-1, j, nx)] = data->p[IDX2D(ny-2, j, nx)];
    }
#pragma omp parallel for
    for (int i = 1; i < ny-1; i++) {
        data->p[IDX2D(i, 0, nx)] = data->p[IDX2D(i, 1, nx)];
        data->p[IDX2D(i, nx-1, nx)] = data->p[IDX2D(i, nx-2, nx)];
    }
}
