#include "ns_spectral_solver.h"
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// Compute CFL timestep
double ns_spectral_compute_cfl_timestep(NSSpectralData *data) {
    const int nx = data->nx;
    const int ny = data->ny;
    double max_vel = 0.0;

#pragma omp parallel for reduction(max:max_vel)
    for (int k = 0; k < nx * ny; k++) {
        const double vel_mag = sqrt(data->u[k] * data->u[k] + data->v[k] * data->v[k]);
        if (vel_mag > max_vel) {
            max_vel = vel_mag;
        }
    }

    const double dx_min = fmin(data->Lx / nx, data->Ly / ny);
    const double dt_cfl = data->cfl_factor * dx_min / (max_vel + 1e-12);
    const double dt_visc = 0.5 * data->Re * dx_min * dx_min;

    double dt_new = fmin(dt_cfl, dt_visc);
    dt_new = fmax(dt_new, MIN_DT);
    dt_new = fmin(dt_new, MAX_DT);

    return dt_new;
}

// 4th-order Runge-Kutta time advancement
void ns_spectral_advance_rk4(NSSpectralData *data) {
    const int spectral_size = data->nkx * data->nky;
    const double dt = data->dt;

    fftw_complex *NS_RESTRICT k1 = data->rk1_hat;
    fftw_complex *NS_RESTRICT k2 = data->rk2_hat;
    fftw_complex *NS_RESTRICT k3 = data->rk3_hat;
    fftw_complex *NS_RESTRICT k4 = data->rk4_hat;
    fftw_complex *NS_RESTRICT omega_tmp = data->omega_tmp_hat;

    memcpy(omega_tmp, data->omega_hat, spectral_size * sizeof(fftw_complex));

    const double inv_Re = 1.0 / data->Re;

    // Stage 1
    ns_spectral_compute_nonlinear_terms(data);
#pragma omp parallel for simd
    for (int k = 0; k < spectral_size; k++) {
        k1[k] = data->nonlinear_hat[k] - data->k2[k] * data->omega_hat[k] * inv_Re;
    }
    ns_spectral_hyperviscosity_apply(data, k1, data->omega_hat);

    // Stage 2
#pragma omp parallel for simd
    for (int k = 0; k < spectral_size; k++) {
        data->omega_hat[k] = omega_tmp[k] + 0.5 * dt * k1[k];
    }
    ns_spectral_compute_derivatives(data);
    ns_spectral_compute_nonlinear_terms(data);
#pragma omp parallel for simd
    for (int k = 0; k < spectral_size; k++) {
        k2[k] = data->nonlinear_hat[k] - data->k2[k] * data->omega_hat[k] * inv_Re;
    }
    ns_spectral_hyperviscosity_apply(data, k2, data->omega_hat);

    // Stage 3
#pragma omp parallel for simd
    for (int k = 0; k < spectral_size; k++) {
        data->omega_hat[k] = omega_tmp[k] + 0.5 * dt * k2[k];
    }
    ns_spectral_compute_derivatives(data);
    ns_spectral_compute_nonlinear_terms(data);
#pragma omp parallel for simd
    for (int k = 0; k < spectral_size; k++) {
        k3[k] = data->nonlinear_hat[k] - data->k2[k] * data->omega_hat[k] * inv_Re;
    }
    ns_spectral_hyperviscosity_apply(data, k3, data->omega_hat);

    // Stage 4
#pragma omp parallel for simd
    for (int k = 0; k < spectral_size; k++) {
        data->omega_hat[k] = omega_tmp[k] + dt * k3[k];
    }
    ns_spectral_compute_derivatives(data);
    ns_spectral_compute_nonlinear_terms(data);
#pragma omp parallel for simd
    for (int k = 0; k < spectral_size; k++) {
        k4[k] = data->nonlinear_hat[k] - data->k2[k] * data->omega_hat[k] * inv_Re;
    }
    ns_spectral_hyperviscosity_apply(data, k4, data->omega_hat);

    // Final update
#pragma omp parallel for simd
    for (int k = 0; k < spectral_size; k++) {
        data->omega_hat[k] =
            omega_tmp[k] + dt / 6.0 * (k1[k] + 2.0 * k2[k] + 2.0 * k3[k] + k4[k]);
    }

    ns_spectral_hyperviscosity_post_step(data, data->omega_hat, dt);

    data->t += dt;
    ns_spectral_compute_derivatives(data);
}
