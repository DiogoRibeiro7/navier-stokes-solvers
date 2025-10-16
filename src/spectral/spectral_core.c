#include "ns_spectral_solver.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// Compute velocity from vorticity via stream function
void ns_spectral_compute_velocity_from_vorticity(NSSpectralData *data) {
    int spectral_size = data->nkx * data->nky;
    
    // Solve Poisson: ∇²ψ = -ω
#pragma omp parallel for collapse(2)
    for (int j = 0; j < data->nky; j++) {
        for (int i = 0; i < data->nkx; i++) {
            const int idx = j * data->nkx + i;

            if (data->k2[idx] > 1e-12) {
                data->psi_hat[idx] = -data->omega_hat[idx] / data->k2[idx];
            } else {
                data->psi_hat[idx] = 0.0 + 0.0 * I;
            }
        }
    }
    
    // Compute velocity: u = ∂ψ/∂y, v = -∂ψ/∂x
#pragma omp parallel for collapse(2)
    for (int j = 0; j < data->nky; j++) {
        for (int i = 0; i < data->nkx; i++) {
            const int idx = j * data->nkx + i;
            data->u_hat[idx] = I * data->ky[j] * data->psi_hat[idx];
            data->v_hat[idx] = -I * data->kx[i] * data->psi_hat[idx];
        }
    }
}

// Compute derivatives in spectral space
void ns_spectral_compute_derivatives(NSSpectralData *data) {
    ns_spectral_compute_velocity_from_vorticity(data);
    
    ns_spectral_spectral_to_physical(data, data->u_hat, data->u);
    ns_spectral_spectral_to_physical(data, data->v_hat, data->v);
    ns_spectral_spectral_to_physical(data, data->omega_hat, data->omega);
}

// Compute nonlinear terms: -J(ψ,ω)
void ns_spectral_compute_nonlinear_terms(NSSpectralData *data) {
    const int nx = data->nx;
    const int ny = data->ny;
    const int spectral_size = data->nkx * data->nky;

    ns_spectral_spectral_to_physical(data, data->psi_hat, data->psi);
    ns_spectral_spectral_to_physical(data, data->omega_hat, data->omega);
    
    double *psi_x = data->u;
    double *psi_y = data->v;
    double *omega_x = data->psi;
    double *omega_y = data->omega;
    double *jacobian = data->p;
    fftw_complex *scratch = data->nonlinear_hat;

    // ∂ψ/∂x
#pragma omp parallel for collapse(2)
    for (int j = 0; j < data->nky; j++) {
        for (int i = 0; i < data->nkx; i++) {
            const int idx = j * data->nkx + i;
            scratch[idx] = I * data->kx[i] * data->psi_hat[idx];
        }
    }
    ns_spectral_spectral_to_physical(data, scratch, psi_x);

    // ∂ψ/∂y
#pragma omp parallel for collapse(2)
    for (int j = 0; j < data->nky; j++) {
        for (int i = 0; i < data->nkx; i++) {
            const int idx = j * data->nkx + i;
            scratch[idx] = I * data->ky[j] * data->psi_hat[idx];
        }
    }
    ns_spectral_spectral_to_physical(data, scratch, psi_y);

    // ∂ω/∂x
#pragma omp parallel for collapse(2)
    for (int j = 0; j < data->nky; j++) {
        for (int i = 0; i < data->nkx; i++) {
            const int idx = j * data->nkx + i;
            scratch[idx] = I * data->kx[i] * data->omega_hat[idx];
        }
    }
    ns_spectral_spectral_to_physical(data, scratch, omega_x);

    // ∂ω/∂y
#pragma omp parallel for collapse(2)
    for (int j = 0; j < data->nky; j++) {
        for (int i = 0; i < data->nkx; i++) {
            const int idx = j * data->nkx + i;
            scratch[idx] = I * data->ky[j] * data->omega_hat[idx];
        }
    }
    ns_spectral_spectral_to_physical(data, scratch, omega_y);

#pragma omp parallel for simd
    for (int k = 0; k < nx * ny; k++) {
        jacobian[k] = psi_x[k] * omega_y[k] - psi_y[k] * omega_x[k];
    }

    ns_spectral_physical_to_spectral(data, jacobian, data->nonlinear_hat);

#pragma omp parallel for simd
    for (int k = 0; k < spectral_size; k++) {
        data->nonlinear_hat[k] = -data->nonlinear_hat[k];
    }
}
