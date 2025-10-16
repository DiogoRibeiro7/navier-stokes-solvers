#include "ns_spectral_solver.h"
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// Physical to spectral transform
void ns_spectral_physical_to_spectral(NSSpectralData *data, double *field, fftw_complex *field_hat) {
    fftw_execute_dft_r2c(data->forward_plan, field, field_hat);
    
    // Normalize
    const double norm = 1.0 / (data->nx * data->ny);
    const int spectral_size = data->nkx * data->nky;

#pragma omp parallel for simd
    for (int k = 0; k < spectral_size; k++) {
        field_hat[k] *= norm;
    }
    
    ns_spectral_apply_dealiasing(data, field_hat);
}

// Spectral to physical transform
void ns_spectral_spectral_to_physical(NSSpectralData *data, fftw_complex *field_hat, double *field) {
    const int spectral_size = data->nkx * data->nky;
    fftw_complex *scratch = data->nonlinear_hat;
    memcpy(scratch, field_hat, spectral_size * sizeof(fftw_complex));

    ns_spectral_apply_dealiasing(data, scratch);
    fftw_execute_dft_c2r(data->backward_plan, scratch, field);
}

// Apply 2/3 dealiasing rule
void ns_spectral_apply_dealiasing(NSSpectralData *data, fftw_complex *field_hat) {
    const int spectral_size = data->nkx * data->nky;
#pragma omp parallel for simd
    for (int k = 0; k < spectral_size; k++) {
        if (!data->dealias_mask[k]) {
            field_hat[k] = 0.0 + 0.0 * I;
        }
    }
}
