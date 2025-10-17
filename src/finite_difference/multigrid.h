#ifndef NS_FD_MULTIGRID_H
#define NS_FD_MULTIGRID_H

#include "../../include/common_types.h"

struct NSFiniteDiffData;

typedef enum {
    NS_MG_CYCLE_V = 0,
    NS_MG_CYCLE_W = 1
} NSMultigridCycleType;

typedef struct {
    int pre_smoothing;
    int post_smoothing;
    int coarse_smoothing;
    int max_cycles;
    double tolerance;
    NSMultigridCycleType cycle_type;
} NSMultigridConfig;

typedef struct {
    double *history;
    int capacity;
    int count;
} NSMultigridMonitor;

typedef struct {
    double initial_residual;
    double final_residual;
    int cycles_performed;
    int converged;
} NSMultigridStats;

void ns_fd_multigrid_default_config(NSMultigridConfig *config);

NSMultigridStats ns_fd_multigrid_solve_poisson(struct NSFiniteDiffData *data,
                                               const double *rhs,
                                               double *solution,
                                               const NSMultigridConfig *config,
                                               NSMultigridMonitor *monitor);

#endif // NS_FD_MULTIGRID_H
