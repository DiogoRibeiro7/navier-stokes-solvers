#include "ns_spectral_solver.h"

// Compute velocity from vorticity via stream function
void ns_spectral_compute_velocity_from_vorticity(NSSpectralData *data) {
    int spectral_size = data->nkx * data->nky;
    
    // Solve Poisson: ∇²ψ = -ω
    for (int j = 0; j < data->nky; j++) {
        for (int i = 0; i < data->nkx; i++) {
            int idx = j * data->nkx + i;
            
            if (data->k2[idx] > 1e-12) {
                data->psi_hat[idx] = -data->omega_hat[idx] / data->k2[idx];
            } else {
                data->psi_hat[idx] = 0.0 + 0.0*I;
            }
        }
    }
    
    // Compute velocity: u = ∂ψ/∂y, v = -∂ψ/∂x
    for (int j = 0; j < data->nky; j++) {
        for (int i = 0; i < data->nkx; i++) {
            int idx = j * data->nkx + i;
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
    int nx = data->nx, ny = data->ny;
    
    ns_spectral_spectral_to_physical(data, data->psi_hat, data->psi);
    ns_spectral_spectral_to_physical(data, data->omega_hat, data->omega);
    
    double *psi_x = fftw_alloc_real(nx * ny);
    double *psi_y = fftw_alloc_real(nx * ny);
    double *omega_x = fftw_alloc_real(nx * ny);
    double *omega_y = fftw_alloc_real(nx * ny);
    double *jacobian = fftw_alloc_real(nx * ny);
    
    fftw_complex *temp_hat = fftw_alloc_complex(data->nkx * data->nky);
    
    // ∂ψ/∂x
    for (int j = 0; j < data->nky; j++) {
        for (int i = 0; i < data->nkx; i++) {
            int idx = j * data->nkx + i;
            temp_hat[idx] = I * data->kx[i] * data->psi_hat[idx];
        }
    }
    ns_spectral_spectral_to_physical(data, temp_hat, psi_x);
    
    // ∂ψ/∂y
    for (int j = 0; j < data->nky; j++) {
        for (int i = 0; i < data->nkx; i++) {
            int idx = j * data->nkx + i;
            temp_hat[idx] = I * data->ky[j] * data->psi_hat[idx];
        }
    }
    ns_spectral_spectral_to_physical(data, temp_hat, psi_y);
    
    // ∂ω/∂x
    for (int j = 0; j < data->nky; j++) {
        for (int i = 0; i < data->nkx; i++) {
            int idx = j * data->nkx + i;
            temp_hat[idx] = I * data->kx[i] * data->omega_hat[idx];
        }
    }
    ns_spectral_spectral_to_physical(data, temp_hat, omega_x);
    
    // ∂ω/∂y
    for (int j = 0; j < data->nky; j++) {
        for (int i = 0; i < data->nkx; i++) {
            int idx = j * data->nkx + i;
            temp_hat[idx] = I * data->ky[j] * data->omega_hat[idx];
        }
    }
    ns_spectral_spectral_to_physical(data, temp_hat, omega_y);
    
    // J(ψ,ω) = ∂ψ/∂x * ∂ω/∂y - ∂ψ/∂y * ∂ω/∂x
    for (int k = 0; k < nx * ny; k++) {
        jacobian[k] = psi_x[k] * omega_y[k] - psi_y[k] * omega_x[k];
    }
    
    ns_spectral_physical_to_spectral(data, jacobian, data->nonlinear_hat);
    
    // Negate for RHS
    int spectral_size = data->nkx * data->nky;
    for (int k = 0; k < spectral_size; k++) {
        data->nonlinear_hat[k] = -data->nonlinear_hat[k];
    }
    
    fftw_free(psi_x); fftw_free(psi_y);
    fftw_free(omega_x); fftw_free(omega_y);
    fftw_free(jacobian); fftw_free(temp_hat);
}
