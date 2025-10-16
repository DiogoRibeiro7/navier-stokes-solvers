#include "ns_spectral_solver.h"

// Physical to spectral transform
void ns_spectral_physical_to_spectral(NSSpectralData *data, double *field, fftw_complex *field_hat) {
    fftw_execute_dft_r2c(data->forward_plan, field, field_hat);
    
    // Normalize
    double norm = 1.0 / (data->nx * data->ny);
    int spectral_size = data->nkx * data->nky;
    for (int k = 0; k < spectral_size; k++) {
        field_hat[k] *= norm;
    }
    
    ns_spectral_apply_dealiasing(data, field_hat);
}

// Spectral to physical transform
void ns_spectral_spectral_to_physical(NSSpectralData *data, fftw_complex *field_hat, double *field) {
    int spectral_size = data->nkx * data->nky;
    fftw_complex *temp_hat = fftw_alloc_complex(spectral_size);
    memcpy(temp_hat, field_hat, spectral_size * sizeof(fftw_complex));
    
    ns_spectral_apply_dealiasing(data, temp_hat);
    fftw_execute_dft_c2r(data->backward_plan, temp_hat, field);
    
    fftw_free(temp_hat);
}

// Apply 2/3 dealiasing rule
void ns_spectral_apply_dealiasing(NSSpectralData *data, fftw_complex *field_hat) {
    int spectral_size = data->nkx * data->nky;
    for (int k = 0; k < spectral_size; k++) {
        if (!data->dealias_mask[k]) {
            field_hat[k] = 0.0 + 0.0*I;
        }
    }
}
