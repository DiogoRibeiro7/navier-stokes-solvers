#include "ns_spectral_solver.h"
#include <math.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static inline double ns_hv_clamp(double value, double min_val, double max_val) {
    if (value < min_val) return min_val;
    if (value > max_val) return max_val;
    return value;
}

void ns_spectral_hyperviscosity_default(NSSpectralHyperviscosityConfig *config) {
    if (!config) {
        return;
    }

    config->enabled = 0;
    config->coefficient = 0.0;
    config->order = 8;
    config->cutoff_ratio = 0.65;
    config->sfd_strength = 0.0;
    config->sfd_timescale = 1.0;
    config->sfd_enabled = 0;
    config->kmax = 0.0;
    config->cutoff_wavenumber = 0.0;
    config->inv_bandwidth = 0.0;
}

void ns_spectral_hyperviscosity_set(NSSpectralData *data,
                                    const NSSpectralHyperviscosityConfig *config) {
    if (!data || !config) {
        return;
    }

    data->hyperviscosity = *config;

    if (!data->hyperviscosity.enabled) {
        memset(data->hyperviscosity_weights, 0,
               (size_t)data->nkx * (size_t)data->nky * sizeof(double));
    }

    if (!data->hyperviscosity.sfd_enabled || data->hyperviscosity.sfd_strength <= 0.0) {
        data->hyperviscosity.sfd_enabled = 0;
        data->hyperviscosity.sfd_strength = 0.0;
        memset(data->sfd_hat, 0,
               (size_t)data->nkx * (size_t)data->nky * sizeof(fftw_complex));
    }

    ns_spectral_hyperviscosity_update_weights(data);

    if (data->hyperviscosity.sfd_enabled && data->sfd_hat && data->omega_hat) {
        const size_t spectral_size = (size_t)data->nkx * (size_t)data->nky;
        memcpy(data->sfd_hat, data->omega_hat, spectral_size * sizeof(fftw_complex));
    }
}

void ns_spectral_hyperviscosity_update_weights(NSSpectralData *data) {
    if (!data || !data->hyperviscosity_weights || !data->k2) {
        return;
    }

    const int spectral_size = data->nkx * data->nky;
    NSSpectralHyperviscosityConfig *cfg = &data->hyperviscosity;

    double kmax = 0.0;
#pragma omp parallel for reduction(max:kmax)
    for (int idx = 0; idx < spectral_size; ++idx) {
        kmax = fmax(kmax, sqrt(data->k2[idx]));
    }
    cfg->kmax = kmax;

    const double cutoff_ratio = ns_hv_clamp(cfg->cutoff_ratio, 0.0, 1.0);
    cfg->cutoff_wavenumber = cutoff_ratio * kmax;
    const double bandwidth = fmax(kmax - cfg->cutoff_wavenumber, 1e-12);
    cfg->inv_bandwidth = (bandwidth > 0.0) ? 1.0 / bandwidth : 0.0;

    if (!cfg->enabled || cfg->coefficient <= 0.0 || kmax <= 0.0) {
#pragma omp parallel for
        for (int idx = 0; idx < spectral_size; ++idx) {
            data->hyperviscosity_weights[idx] = 0.0;
        }
        return;
    }

    const int order = (cfg->order > 0) ? cfg->order : 1;

#pragma omp parallel for
    for (int idx = 0; idx < spectral_size; ++idx) {
        const double k_mag = sqrt(data->k2[idx]);
        double weight = 0.0;

        if (k_mag > cfg->cutoff_wavenumber) {
            const double ratio = (k_mag - cfg->cutoff_wavenumber) * cfg->inv_bandwidth;
            weight = cfg->coefficient * pow(ratio, order);
        }

        data->hyperviscosity_weights[idx] = weight;
    }
}

void ns_spectral_hyperviscosity_apply(NSSpectralData *data,
                                      fftw_complex *rhs_hat,
                                      const fftw_complex *state_hat) {
    if (!data || !rhs_hat || !state_hat) {
        return;
    }

    const NSSpectralHyperviscosityConfig *cfg = &data->hyperviscosity;
    const double *weights = data->hyperviscosity_weights;
    const fftw_complex *sfd = data->sfd_hat;
    const int spectral_size = data->nkx * data->nky;

    const int apply_visc = (cfg->enabled && cfg->coefficient > 0.0);
    const int apply_sfd = (cfg->sfd_enabled && cfg->sfd_strength > 0.0);
    const int order = (cfg->order > 0) ? cfg->order : 1;

    if (!apply_visc && !apply_sfd) {
        return;
    }

#pragma omp parallel for simd
    for (int idx = 0; idx < spectral_size; ++idx) {
        fftw_complex damping = 0.0 + 0.0 * I;

        if (apply_visc) {
            damping += weights[idx] * state_hat[idx];
        }

        if (apply_sfd) {
            double ratio = 0.0;
            if (cfg->inv_bandwidth > 0.0) {
                const double k_mag = sqrt(data->k2[idx]);
                if (k_mag > cfg->cutoff_wavenumber) {
                    double scaled = (k_mag - cfg->cutoff_wavenumber) * cfg->inv_bandwidth;
                    scaled = ns_hv_clamp(scaled, 0.0, 1.0);
                    ratio = pow(scaled, order);
                }
            }
            const double sfd_weight = cfg->sfd_strength * ratio;
            if (sfd_weight > 0.0) {
                damping += sfd_weight * (state_hat[idx] - sfd[idx]);
            }
        }

        rhs_hat[idx] -= damping;
    }
}

void ns_spectral_hyperviscosity_post_step(NSSpectralData *data,
                                          const fftw_complex *state_hat,
                                          double dt) {
    if (!data || !state_hat || dt <= 0.0) {
        return;
    }

    NSSpectralHyperviscosityConfig *cfg = &data->hyperviscosity;
    if (!cfg->sfd_enabled || cfg->sfd_strength <= 0.0 || cfg->sfd_timescale <= 0.0) {
        return;
    }

    double alpha = dt / cfg->sfd_timescale;
    if (alpha > 1.0) {
        alpha = 1.0;
    } else if (alpha < 0.0) {
        alpha = 0.0;
    }

    const int spectral_size = data->nkx * data->nky;

#pragma omp parallel for simd
    for (int idx = 0; idx < spectral_size; ++idx) {
        const fftw_complex previous = data->sfd_hat[idx];
        data->sfd_hat[idx] = (1.0 - alpha) * previous + alpha * state_hat[idx];
    }
}
