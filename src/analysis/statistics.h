#ifndef NS_ANALYSIS_STATISTICS_H
#define NS_ANALYSIS_STATISTICS_H

#include "../../include/common_types.h"

#ifdef __cplusplus
extern "C" {
#endif

// Central moments of a scalar field
typedef struct {
    double mean;
    double rms;
    double skewness;
    double kurtosis;
} NSFieldMoments;

// Collection of turbulence statistics for a 2D velocity field
typedef struct {
    NSFieldMoments u;
    NSFieldMoments v;
    double integral_length_scale_x;
    double integral_length_scale_y;
    double taylor_microscale;
    double kolmogorov_scale;
} NSTurbulenceStatistics;

/**
 * Compute turbulence statistics for a 2D velocity field.
 *
 * The optional correlation buffers (corr_x, corr_y) store the one-dimensional,
 * normalized two-point correlations of the streamwise velocity component along
 * the x and y directions respectively. The correlation arrays must have
 * max_sep_x and max_sep_y entries if provided.
 *
 * Returns 0 on success, non-zero on error (invalid input or allocation failure).
 */
int ns_stats_compute_velocity_statistics(const double *NS_RESTRICT u,
                                         const double *NS_RESTRICT v,
                                         int nx, int ny,
                                         double dx, double dy,
                                         double nu,
                                         int max_sep_x, int max_sep_y,
                                         NSTurbulenceStatistics *stats,
                                         double *corr_x, double *corr_y);

#ifdef __cplusplus
}
#endif

#endif // NS_ANALYSIS_STATISTICS_H
