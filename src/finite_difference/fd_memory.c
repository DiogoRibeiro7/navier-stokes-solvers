#include "ns_fd_solver.h"

// Allocate memory for finite difference solver
NSFiniteDiffData* ns_fd_allocate(int nx, int ny, double L, double H, double Re) {
    NSFiniteDiffData *data = malloc(sizeof(NSFiniteDiffData));
    if (!data) return NULL;
    
    data->nx = nx;
    data->ny = ny;
    data->dx = L / (nx - 1);
    data->dy = H / (ny - 1);
    data->Re = Re;
    data->dt = 0.01;
    data->t = 0.0;
    data->L = L;
    data->H = H;
    
    int total_points = nx * ny;
    
    // Allocate field arrays
    data->u = calloc(total_points, sizeof(double));
    data->v = calloc(total_points, sizeof(double));
    data->p = calloc(total_points, sizeof(double));
    data->u_old = calloc(total_points, sizeof(double));
    data->v_old = calloc(total_points, sizeof(double));
    data->p_old = calloc(total_points, sizeof(double));
    data->F = calloc(2 * total_points, sizeof(double));
    data->G = calloc(total_points, sizeof(double));
    
    // Sparse Jacobian matrix (conservative estimate)
    int nnz = 15 * total_points;
    data->J = calloc(nnz, sizeof(double));
    data->ia = calloc(3 * total_points + 1, sizeof(int));
    data->ja = calloc(nnz, sizeof(int));
    
    // Check allocations
    if (!data->u || !data->v || !data->p || !data->u_old || !data->v_old || 
        !data->p_old || !data->F || !data->G || !data->J || !data->ia || !data->ja) {
        ns_fd_free(data);
        return NULL;
    }
    
    // Default boundary conditions
    data->bc_type = BC_NO_SLIP;
    data->bc_values[0] = 0.0; // bottom
    data->bc_values[1] = 1.0; // top (lid)
    data->bc_values[2] = 0.0; // left
    data->bc_values[3] = 0.0; // right
    
    return data;
}

// Free all allocated memory
void ns_fd_free(NSFiniteDiffData *data) {
    if (!data) return;
    
    free(data->u); free(data->v); free(data->p);
    free(data->u_old); free(data->v_old); free(data->p_old);
    free(data->F); free(data->G);
    free(data->J); free(data->ia); free(data->ja);
    free(data);
}
