#ifndef NS_COMMON_FORCING_H
#define NS_COMMON_FORCING_H

#include "../../include/common_types.h"
#include <complex.h>

typedef enum {
    NS_FORCING_WHITE_NOISE = 0,
    NS_FORCING_COLORED_NOISE = 1
} NSForcingType;

typedef struct {
    double kmin;
    double kmax;
    double energy_target;
    double correlation_time;
    unsigned long seed;
    NSForcingType type;
} NSForcingConfig;

typedef struct {
    NSForcingConfig config;
    int nkx;
    int nky;
    int spectral_size;
    int active_modes;
    double *weights;
    double complex *colored_state;
    unsigned int have_spare;
    double spare_gaussian;
    unsigned long long rng_state;
} NSForcingState;

void ns_forcing_config_default(NSForcingConfig *config);

int ns_forcing_state_init(NSForcingState *state,
                          const NSForcingConfig *config,
                          const double *k2,
                          int nkx, int nky);

void ns_forcing_state_free(NSForcingState *state);

void ns_forcing_generate(NSForcingState *state,
                         double dt,
                         const double *k2,
                         double complex *forcing_out);

#endif // NS_COMMON_FORCING_H
