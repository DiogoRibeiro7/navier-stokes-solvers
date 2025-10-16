#include "ns_spectral_solver.h"

// Compute CFL timestep
double ns_spectral_compute_cfl_timestep(NSSpectralData *data) {
    int nx = data->nx, ny = data->ny;
    double max_vel = 0.0;
    
    for (int k = 0; k < nx * ny; k++) {
        double vel_mag = sqrt(data->u[k]*data->u[k] + data->v[k]*data->v[k]);
        max_vel = fmax(max_vel, vel_mag);
    }
    
    double dx_min = fmin(data->Lx/nx, data->Ly/ny);
    double dt_cfl = data->cfl_factor * dx_min / (max_vel + 1e-12);
    double dt_visc = 0.5 * data->Re * dx_min * dx_min;
    
    double dt_new = fmin(dt_cfl, dt_visc);
    dt_new = fmax(dt_new, MIN_DT);
    dt_new = fmin(dt_new, MAX_DT);
    
    return dt_new;
}

// 4th-order Runge-Kutta time advancement
void ns_spectral_advance_rk4(NSSpectralData *data) {
    int spectral_size = data->nkx * data->nky;
    double dt = data->dt;
    
    fftw_complex *k1 = fftw_alloc_complex(spectral_size);
    fftw_complex *k2 = fftw_alloc_complex(spectral_size);
    fftw_complex *k3 = fftw_alloc_complex(spectral_size);
    fftw_complex *k4 = fftw_alloc_complex(spectral_size);
    fftw_complex *omega_temp = fftw_alloc_complex(spectral_size);
    
    memcpy(omega_temp, data->omega_hat, spectral_size * sizeof(fftw_complex));
    
    // Stage 1
    ns_spectral_compute_nonlinear_terms(data);
    for (int k = 0; k < spectral_size; k++) {
        k1[k] = data->nonlinear_hat[k] - data->k2[k] * data->omega_hat[k] / data->Re;
    }
    
    // Stage 2
    for (int k = 0; k < spectral_size; k++) {
        data->omega_hat[k] = omega_temp[k] + 0.5 * dt * k1[k];
    }
    ns_spectral_compute_derivatives(data);
    ns_spectral_compute_nonlinear_terms(data);
    for (int k = 0; k < spectral_size; k++) {
        k2[k] = data->nonlinear_hat[k] - data->k2[k] * data->omega_hat[k] / data->Re;
    }
    
    // Stage 3
    for (int k = 0; k < spectral_size; k++) {
        data->omega_hat[k] = omega_temp[k] + 0.5 * dt * k2[k];
    }
    ns_spectral_compute_derivatives(data);
    ns_spectral_compute_nonlinear_terms(data);
    for (int k = 0; k < spectral_size; k++) {
        k3[k] = data->nonlinear_hat[k] - data->k2[k] * data->omega_hat[k] / data->Re;
    }
    
    // Stage 4
    for (int k = 0; k < spectral_size; k++) {
        data->omega_hat[k] = omega_temp[k] + dt * k3[k];
    }
    ns_spectral_compute_derivatives(data);
    ns_spectral_compute_nonlinear_terms(data);
    for (int k = 0; k < spectral_size; k++) {
        k4[k] = data->nonlinear_hat[k] - data->k2[k] * data->omega_hat[k] / data->Re;
    }
    
    // Final update
    for (int k = 0; k < spectral_size; k++) {
        data->omega_hat[k] = omega_temp[k] + dt/6.0 * (k1[k] + 2.0*k2[k] + 2.0*k3[k] + k4[k]);
    }
    
    data->t += dt;
    ns_spectral_compute_derivatives(data);
    
    fftw_free(k1); fftw_free(k2); fftw_free(k3); fftw_free(k4);
    fftw_free(omega_temp);
}
