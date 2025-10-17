#include "forcing.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define NS_FORCE_EPS 1e-16

static uint64_t ns_forcing_splitmix64(uint64_t *state) {
    uint64_t z = (*state += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
}

static double ns_forcing_uniform(NSForcingState *state) {
    uint64_t bits = ns_forcing_splitmix64(&state->rng_state);
    return (bits >> 11) * (1.0 / 9007199254740992.0); // 2^53
}

static double ns_forcing_gaussian(NSForcingState *state) {
    if (state->have_spare) {
        state->have_spare = 0;
        return state->spare_gaussian;
    }

    double u1 = ns_forcing_uniform(state);
    double u2 = ns_forcing_uniform(state);

    u1 = fmax(u1, DBL_MIN);
    const double radius = sqrt(-2.0 * log(u1));
    const double theta = 2.0 * M_PI * u2;

    state->spare_gaussian = radius * sin(theta);
    state->have_spare = 1;
    return radius * cos(theta);
}

void ns_forcing_config_default(NSForcingConfig *config) {
    if (!config) {
        return;
    }
    config->kmin = 0.0;
    config->kmax = 4.0;
    config->energy_target = 1.0;
    config->correlation_time = 0.5;
    config->seed = 0;
    config->type = NS_FORCING_WHITE_NOISE;
}

int ns_forcing_state_init(NSForcingState *state,
                          const NSForcingConfig *config,
                          const double *k2,
                          int nkx, int nky) {
    if (!state || !config || !k2 || nkx <= 0 || nky <= 0) {
        return -1;
    }

    memset(state, 0, sizeof(*state));
    state->config = *config;
    state->nkx = nkx;
    state->nky = nky;
    state->spectral_size = nkx * nky;
    state->weights = calloc((size_t)state->spectral_size, sizeof(double));
    if (!state->weights) {
        return -1;
    }

    if (config->type == NS_FORCING_COLORED_NOISE) {
        state->colored_state = calloc((size_t)state->spectral_size, sizeof(double complex));
        if (!state->colored_state) {
            ns_forcing_state_free(state);
            return -1;
        }
    }

    const double kmin = fmax(config->kmin, 0.0);
    const double kmax = fmax(config->kmax, kmin);

    int active = 0;
    for (int j = 0; j < nky; ++j) {
        for (int i = 0; i < nkx; ++i) {
            const int idx = j * nkx + i;
            const double kmag = sqrt(k2[idx]);
            if (kmag >= kmin && kmag <= kmax) {
                state->weights[idx] = 1.0;
                active++;
            } else {
                state->weights[idx] = 0.0;
            }
        }
    }

    state->active_modes = active;
    state->have_spare = 0;

    uint64_t seed = config->seed ? config->seed : (uint64_t)0x123456789abcdefULL;
    state->rng_state = seed;

    return (active == 0) ? -1 : 0;
}

void ns_forcing_state_free(NSForcingState *state) {
    if (!state) {
        return;
    }
    free(state->weights);
    state->weights = NULL;
    free(state->colored_state);
    state->colored_state = NULL;
    state->spectral_size = 0;
    state->active_modes = 0;
}

void ns_forcing_generate(NSForcingState *state,
                         double dt,
                         const double *k2,
                         double complex *forcing_out) {
    if (!state || !forcing_out || !k2) {
        return;
    }

    const int size = state->spectral_size;
    if (size == 0 || state->active_modes == 0) {
        memset(forcing_out, 0, (size_t)size * sizeof(double complex));
        return;
    }

    const NSForcingConfig *cfg = &state->config;
    if (cfg->energy_target <= 0.0) {
        memset(forcing_out, 0, (size_t)size * sizeof(double complex));
        if (state->colored_state) {
            memset(state->colored_state, 0, (size_t)size * sizeof(double complex));
        }
        return;
    }
    const double tau = fmax(cfg->correlation_time, 1e-6);
    const int colored = (cfg->type == NS_FORCING_COLORED_NOISE);

    double alpha = 0.0;
    double beta = 1.0;
    if (colored) {
        alpha = exp(-dt / tau);
        beta = sqrt(fmax(1.0 - alpha * alpha, 0.0));
    }

    double total_energy = 0.0;

    for (int idx = 0; idx < size; ++idx) {
        if (state->weights[idx] <= 0.0) {
            forcing_out[idx] = 0.0 + 0.0 * I;
            continue;
        }

        double g1 = ns_forcing_gaussian(state);
        double g2 = ns_forcing_gaussian(state);
        double complex noise = g1 + I * g2;

        double complex value;
        if (colored) {
            state->colored_state[idx] =
                alpha * state->colored_state[idx] + beta * noise;
            value = state->colored_state[idx];
        } else {
            value = noise;
        }

        forcing_out[idx] = state->weights[idx] * value;
        total_energy += pow(state->weights[idx], 2.0) *
                        (creal(value) * creal(value) + cimag(value) * cimag(value));
    }

    total_energy = fmax(total_energy, NS_FORCE_EPS);
    const double scale = sqrt(cfg->energy_target / total_energy);

    for (int idx = 0; idx < size; ++idx) {
        forcing_out[idx] *= scale;
        if (colored) {
            state->colored_state[idx] *= scale;
        }
    }
}
