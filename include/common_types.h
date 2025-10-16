#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Common constants
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_ITER 1000
#define TOL 1e-8
#define MIN_DT 1e-8
#define MAX_DT 0.1

// Grid indexing macro for 2D arrays
#define IDX2D(i, j, nx) ((i) * (nx) + (j))

// Boundary condition types
typedef enum {
    BC_DIRICHLET,
    BC_NEUMANN,
    BC_PERIODIC,
    BC_NO_SLIP,
    BC_FREE_SLIP
} BoundaryCondition;

// Solver status codes
typedef enum {
    SOLVER_SUCCESS = 0,
    SOLVER_MAX_ITER,
    SOLVER_DIVERGED,
    SOLVER_MEM_ERROR,
    SOLVER_INVALID_PARAM
} SolverStatus;

// Common solver statistics
typedef struct {
    double energy;
    double enstrophy;
    double max_velocity;
    double max_divergence;
    double residual;
    int iterations;
    double cpu_time;
} SolverStats;

// Output options
typedef struct {
    int output_interval;
    int write_vtk;
    int write_text;
    int write_binary;
    char output_dir[256];
} OutputOptions;

#endif // COMMON_TYPES_H
