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
    data->delta = calloc(3 * total_points, sizeof(double));
    
    // Sparse Jacobian matrix (conservative estimate)
    int nnz = 15 * total_points;
    data->J = calloc(nnz, sizeof(double));
    data->ia = calloc(3 * total_points + 1, sizeof(int));
    data->ja = calloc(nnz, sizeof(int));
    
    // Check allocations
    if (!data->u || !data->v || !data->p || !data->u_old || !data->v_old || 
        !data->p_old || !data->F || !data->G || !data->delta || !data->J || !data->ia || !data->ja) {
        ns_fd_free(data);
        return NULL;
    }

    // Default boundary conditions: no-slip walls with a moving lid at the top
    ns_fd_set_boundary_condition(data, NS_BC_BOTTOM, BC_NO_SLIP, 0.0, 0.0, NAN);
    ns_fd_set_boundary_condition(data, NS_BC_TOP, BC_NO_SLIP, 1.0, 0.0, NAN);
    ns_fd_set_boundary_condition(data, NS_BC_LEFT, BC_NO_SLIP, 0.0, 0.0, NAN);
    ns_fd_set_boundary_condition(data, NS_BC_RIGHT, BC_NO_SLIP, 0.0, 0.0, NAN);
    
    return data;
}

// Free all allocated memory
void ns_fd_free(NSFiniteDiffData *data) {
    if (!data) return;
    
    free(data->u); free(data->v); free(data->p);
    free(data->u_old); free(data->v_old); free(data->p_old);
    free(data->F); free(data->G); free(data->delta);
    free(data->J); free(data->ia); free(data->ja);
    free(data);
}

void ns_fd_set_boundary_condition(NSFiniteDiffData *data, int side,
                                  BoundaryCondition type,
                                  double value_u, double value_v, double value_p) {
    if (!data || side < 0 || side >= NS_BC_COUNT) {
        return;
    }

    NSBoundaryCondition *bc = &data->boundaries[side];
    bc->type = type;
    bc->value_u = value_u;
    bc->value_v = value_v;
    bc->value_p = value_p;

    if (type != BC_CUSTOM) {
        bc->func = NULL;
        bc->user_data = NULL;
    }
}

void ns_fd_set_boundary_function(NSFiniteDiffData *data, int side,
                                 NSBoundaryFunction func, void *user_data) {
    if (!data || side < 0 || side >= NS_BC_COUNT) {
        return;
    }

    data->boundaries[side].func = func;
    data->boundaries[side].user_data = user_data;
}
