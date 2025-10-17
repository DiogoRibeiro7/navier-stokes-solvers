#ifndef NS_FD_SOLVER_H
#define NS_FD_SOLVER_H

#include "../../include/common_types.h"

typedef void (*NSBoundaryFunction)(double x, double y, double t,
                                   double *u, double *v, double *p,
                                   void *user_data);

typedef struct {
    BoundaryCondition type;
    double value_u;
    double value_v;
    double value_p;
    NSBoundaryFunction func;
    void *user_data;
} NSBoundaryCondition;

enum {
    NS_BC_BOTTOM = 0,
    NS_BC_TOP = 1,
    NS_BC_LEFT = 2,
    NS_BC_RIGHT = 3,
    NS_BC_COUNT = 4
};

// Finite difference Navier-Stokes data structure
typedef struct NSFiniteDiffData {

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
    NSBoundaryCondition boundaries[NS_BC_COUNT];
} NSFiniteDiffData;

// Memory management
NSFiniteDiffData* ns_fd_allocate(int nx, int ny, double L, double H, double Re);
void ns_fd_free(NSFiniteDiffData *data);
void ns_fd_set_boundary_condition(NSFiniteDiffData *data, int side,
                                  BoundaryCondition type,
                                  double value_u, double value_v, double value_p);
void ns_fd_set_boundary_function(NSFiniteDiffData *data, int side,
                                 NSBoundaryFunction func, void *user_data);

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
