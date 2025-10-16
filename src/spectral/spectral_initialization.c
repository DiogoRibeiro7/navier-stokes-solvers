#include "ns_spectral_solver.h"

// Initialize grids and wavenumbers
void ns_spectral_initialize_grids(NSSpectralData *data) {
    int nx = data->nx, ny = data->ny;
    double Lx = data->Lx, Ly = data->Ly;
    
    // Physical coordinates (uniform for Fourier)
    for (int i = 0; i < nx; i++) {
        data->x[i] = i * Lx / nx;
    }
    for (int j = 0; j < ny; j++) {
        data->y[j] = j * Ly / ny;
    }
    
    // Wavenumbers for FFT
    for (int i = 0; i < data->nkx; i++) {
        data->kx[i] = 2.0 * M_PI * i / Lx;
    }
    
    for (int j = 0; j < data->nky; j++) {
        if (j <= ny/2) {
            data->ky[j] = 2.0 * M_PI * j / Ly;
        } else {
            data->ky[j] = 2.0 * M_PI * (j - ny) / Ly;
        }
    }
    
    // Compute k^2 and dealiasing mask
    for (int j = 0; j < data->nky; j++) {
        for (int i = 0; i < data->nkx; i++) {
            int idx = j * data->nkx + i;
            data->k2[idx] = data->kx[i] * data->kx[i] + data->ky[j] * data->ky[j];
            
            // 2/3 dealiasing rule
            int dealias_x = (i <= (2 * data->nkx) / 3);
            int dealias_y = (abs(j <= ny/2 ? j : j - ny) <= (2 * ny) / 6);
            data->dealias_mask[idx] = dealias_x && dealias_y;
        }
    }
}

// Initialize FFT plans
void ns_spectral_initialize_fft_plans(NSSpectralData *data) {
    data->forward_plan = fftw_plan_dft_r2c_2d(data->ny, data->nx,
                                              data->omega, data->omega_hat,
                                              FFTW_MEASURE);
    
    data->backward_plan = fftw_plan_dft_c2r_2d(data->ny, data->nx,
                                               data->omega_hat, data->omega,
                                               FFTW_MEASURE);
}

// Initialize Taylor-Green vortex
void ns_spectral_initialize_taylor_green(NSSpectralData *data) {
    int nx = data->nx, ny = data->ny;
    double Lx = data->Lx, Ly = data->Ly;
    
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = j * nx + i;
            double x = data->x[i];
            double y = data->y[j];
            
            data->u[idx] = sin(2*M_PI*x/Lx) * cos(2*M_PI*y/Ly);
            data->v[idx] = -cos(2*M_PI*x/Lx) * sin(2*M_PI*y/Ly);
            data->omega[idx] = 2.0 * (2*M_PI/Lx) * cos(2*M_PI*x/Lx) * cos(2*M_PI*y/Ly);
            data->p[idx] = 0.0;
        }
    }
    
    // Transform to spectral space
    ns_spectral_physical_to_spectral(data, data->omega, data->omega_hat);
    ns_spectral_physical_to_spectral(data, data->u, data->u_hat);
    ns_spectral_physical_to_spectral(data, data->v, data->v_hat);
}
