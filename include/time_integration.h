#ifndef TIME_INTEGRATION_H
#define TIME_INTEGRATION_H

#include "common_types.h"

// Time integration schemes
typedef enum {
    TIME_EULER_EXPLICIT,
    TIME_EULER_IMPLICIT,
    TIME_CRANK_NICOLSON,
    TIME_RK4,
    TIME_ADAMS_BASHFORTH
} TimeScheme;

// Adaptive time stepping parameters
typedef struct {
    double cfl_factor;
    double dt_min;
    double dt_max;
    double safety_factor;
    int adapt_enabled;
} AdaptiveTimeParams;

// Function prototypes for time integration utilities
double compute_cfl_condition(double max_u, double max_v, double dx, double dy, double cfl_factor);
double compute_viscous_timestep(double dx, double dy, double Re);
double adaptive_timestep(double max_vel, double dx, double dy, double Re, AdaptiveTimeParams *params);

#endif // TIME_INTEGRATION_H
