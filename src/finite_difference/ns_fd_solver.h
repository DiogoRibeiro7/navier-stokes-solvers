#ifndef NS_FD_SOLVER_H
#define NS_FD_SOLVER_H

#include "../../include/common_types.h"

// Finite difference Navier-Stokes data structure
typedef struct {
    int nx, ny;                    // Grid dimensions
    double dx, dy;                 // Grid spacing
    double Re;                     // Reynolds number
    double dt;                     // Current time step
    double t;                      // Current time
    double L, H;                   // Domain dimensions
    
    // Flow fields
    double *NS_RESTRICT u;
    double *NS_RESTRICT v;
    double *NS_RESTRICT p;
    double *NS_RESTRICT u_old;
    double *NS_RESTRICT v_old;
    double *NS_RESTRICT p_old;
    double *NS_RESTRICT F;
    double *NS_RESTRICT G;
    double *NS_RESTRICT delta;
    
    // Sparse Jacobian (CSR format)
    double *J;                     // Matrix values
    int *ia, *ja;                  // Row/column indices
    
    // Boundary conditions
    BoundaryCondition bc_type;
    double bc_values[4];           // [bottom, top, left, right]
} NSFiniteDiffData;

// Memory management
NSFiniteDiffData* ns_fd_allocate(int nx, int ny, double L, double H, double Re);
void ns_fd_free(NSFiniteDiffData *data);

// Initialization
void ns_fd_initialize_lid_cavity(NSFiniteDiffData *data);
void ns_fd_initialize_taylor_green(NSFiniteDiffData *data);
void ns_fd_apply_boundary_conditions(NSFiniteDiffData *data);

// Newton-Raphson solver
SolverStatus ns_fd_newton_raphson_step(NSFiniteDiffData *data);
double ns_fd_compute_residual(NSFiniteDiffData *data);
void ns_fd_assemble_jacobian(NSFiniteDiffData *data);
void ns_fd_solve_linear_system(NSFiniteDiffData *data, double *delta);

// Time advancement
SolverStatus ns_fd_advance_timestep(NSFiniteDiffData *data);
double ns_fd_adaptive_timestep(NSFiniteDiffData *data);

// Analysis and output
void ns_fd_compute_stats(NSFiniteDiffData *data, SolverStats *stats);
void ns_fd_analyze_convergence(NSFiniteDiffData *data, int step);
void ns_fd_output_solution(NSFiniteDiffData *data, const char *filename);

#endif // NS_FD_SOLVER_H
