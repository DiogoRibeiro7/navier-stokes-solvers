#include "ns_spectral_solver.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// Compute kinetic energy
double ns_spectral_compute_energy(NSSpectralData *data) {
    int nx = data->nx, ny = data->ny;
    double energy = 0.0;
    
    const int total = nx * ny;
#pragma omp parallel for reduction(+:energy)
    for (int k = 0; k < total; k++) {
        energy += 0.5 * (data->u[k] * data->u[k] + data->v[k] * data->v[k]);
    }
    
    return energy / (nx * ny);
}

// Compute enstrophy
double ns_spectral_compute_enstrophy(NSSpectralData *data) {
    int nx = data->nx, ny = data->ny;
    double enstrophy = 0.0;
    
    const int total = nx * ny;
#pragma omp parallel for reduction(+:enstrophy)
    for (int k = 0; k < total; k++) {
        enstrophy += 0.5 * data->omega[k] * data->omega[k];
    }
    
    return enstrophy / (nx * ny);
}

// Check spectral resolution
int ns_spectral_check_resolution(NSSpectralData *data) {
    int nkx = data->nkx, nky = data->nky;
    double high_k_energy = 0.0;
    double total_energy = 0.0;
    double threshold = 0.8;
    
    const int spectral_size = nkx * nky;
    const double max_k = sqrt(data->k2[spectral_size - 1]);
#pragma omp parallel for collapse(2) reduction(+:total_energy, high_k_energy)
    for (int j = 0; j < nky; j++) {
        for (int i = 0; i < nkx; i++) {
            const int idx = j * nkx + i;
            const double omega_mag2 = cabs(data->omega_hat[idx]) * cabs(data->omega_hat[idx]);
            total_energy += omega_mag2;

            const double k_normalized = (max_k > 0.0) ? sqrt(data->k2[idx]) / max_k : 0.0;
            if (k_normalized > threshold) {
                high_k_energy += omega_mag2;
            }
        }
    }
    
    double high_k_fraction = (total_energy > 1e-16) ? high_k_energy / total_energy : 0.0;
    
    if (high_k_fraction > 1e-6) {
        printf("Warning: %.2e%% of energy in high wavenumbers\n", high_k_fraction * 100);
        return 0;
    }
    
    return 1;
}

// Analyze convergence
void ns_spectral_analyze_convergence(NSSpectralData *data, int step) {
    if (step % 10 != 0) return;
    
    double energy = ns_spectral_compute_energy(data);
    double enstrophy = ns_spectral_compute_enstrophy(data);
    
    int spectral_size = data->nkx * data->nky;
    double max_omega_hat = 0.0;
    double energy_spectrum = 0.0;
    int high_modes = 0;
    
    for (int j = 0; j < data->nky; j++) {
#pragma omp parallel for reduction(max:max_omega_hat) reduction(+:energy_spectrum, high_modes)
        for (int i = 0; i < data->nkx; i++) {
            const int idx = j * data->nkx + i;
            const double omega_mag = cabs(data->omega_hat[idx]);
            if (omega_mag > max_omega_hat) {
                max_omega_hat = omega_mag;
            }

            if (data->k2[idx] > 0.25 * data->k2[spectral_size - 1]) {
                energy_spectrum += omega_mag * omega_mag;
                high_modes++;
            }
        }
    }
    
    double spectral_error = (high_modes > 0) ? sqrt(energy_spectrum / high_modes) : 0.0;
    
    printf("Step %d: t=%.4f, E=%.6e, Ω=%.6e, max|ω̂|=%.2e, err=%.2e, dt=%.2e\n",
           step, data->t, energy, enstrophy, max_omega_hat, spectral_error, data->dt);
}

// Analyze energy spectrum
void ns_spectral_analyze_spectrum(NSSpectralData *data) {
    int nkx = data->nkx, nky = data->nky;
    int max_k = fmin(nkx, nky/2);
    
    double *energy_k = calloc(max_k, sizeof(double));
    int *count_k = calloc(max_k, sizeof(int));
    
    for (int j = 0; j < nky; j++) {
        for (int i = 0; i < nkx; i++) {
            int idx = j * nkx + i;
            double k_mag = sqrt(data->k2[idx]);
            int k_bin = (int)(k_mag + 0.5);
            
            if (k_bin < max_k) {
                double u_mag2 = cabs(data->u_hat[idx]) * cabs(data->u_hat[idx]);
                double v_mag2 = cabs(data->v_hat[idx]) * cabs(data->v_hat[idx]);
                energy_k[k_bin] += 0.5 * (u_mag2 + v_mag2);
                count_k[k_bin]++;
            }
        }
    }
    
    printf("\nEnergy Spectrum (t=%.4f):\n", data->t);
    printf("k\tE(k)\t\tSlope\n");
    
    for (int k = 1; k < max_k && count_k[k] > 0; k++) {
        energy_k[k] /= count_k[k];
        
        double slope = 0.0;
        if (k > 1 && energy_k[k-1] > 1e-16 && energy_k[k] > 1e-16) {
            slope = log(energy_k[k] / energy_k[k-1]) / log((double)k / (k-1));
        }
        
        if (k <= 10 || k % 5 == 0) {
            printf("%d\t%.2e\t%.2f\n", k, energy_k[k], slope);
        }
    }
    
    free(energy_k); free(count_k);
}
