#include "ns_spectral_solver.h"
#include <time.h>

int main(int argc, char *argv[]) {
    printf("Spectral Navier-Stokes Solver\n");
    printf("Fourier Methods with RK4 Time Integration\n\n");
    
    // Problem parameters
    int nx = 128, ny = 128;
    double Lx = 2*M_PI, Ly = 2*M_PI;
    double Re = 1000.0;
    double t_final = 5.0;
    int output_interval = 50;
    
    printf("Resolution: %dx%d modes\n", nx, ny);
    printf("Domain: [0,%.2f] x [0,%.2f]\n", Lx, Ly);
    printf("Reynolds number: %.1f\n\n", Re);
    
    // Initialize solver
    NSSpectralData *data = ns_spectral_allocate(nx, ny, Lx, Ly, Re);
    if (!data) {
        printf("Memory allocation failed!\n");
        return 1;
    }
    
    ns_spectral_initialize_grids(data);
    ns_spectral_initialize_fft_plans(data);
    ns_spectral_initialize_taylor_green(data);
    
    printf("Initial conditions:\n");
    printf("Energy: %.6e\n", ns_spectral_compute_energy(data));
    printf("Enstrophy: %.6e\n\n", ns_spectral_compute_enstrophy(data));
    
    // Time integration
    clock_t start_time = clock();
    int step = 0;
    int resolution_adequate = 1;
    
    while (data->t < t_final && resolution_adequate) {
        data->dt = ns_spectral_compute_cfl_timestep(data);
        ns_spectral_advance_rk4(data);
        step++;
        
        ns_spectral_analyze_convergence(data, step);
        
        if (step % 100 == 0) {
            resolution_adequate = ns_spectral_check_resolution(data);
        }
        
        if (step % output_interval == 0) {
            char filename[100];
            sprintf(filename, "spectral_solution_step_%04d.dat", step);
            ns_spectral_output_solution(data, filename);
            
            if (step % (5 * output_interval) == 0) {
                ns_spectral_analyze_spectrum(data);
            }
        }
        
        if (data->dt < MIN_DT) {
            printf("Time step too small!\n");
            break;
        }
    }
    
    clock_t end_time = clock();
    double cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    
    printf("\n=== SIMULATION COMPLETED ===\n");
    printf("Total steps: %d\n", step);
    printf("Final time: %.6f\n", data->t);
    printf("Final energy: %.6e\n", ns_spectral_compute_energy(data));
    printf("Final enstrophy: %.6e\n", ns_spectral_compute_enstrophy(data));
    printf("CPU time: %.2f seconds\n", cpu_time);
    printf("Steps/second: %.1f\n\n", step / cpu_time);
    
    ns_spectral_output_solution(data, "spectral_final_solution.dat");
    ns_spectral_analyze_spectrum(data);
    
    printf("\nSpectral Accuracy:\n");
    if (ns_spectral_check_resolution(data)) {
        printf("✓ Resolution adequate\n");
        printf("✓ High-order accuracy maintained\n");
    } else {
        printf("⚠ Consider increasing resolution\n");
    }
    
    ns_spectral_free(data);
    
    return 0;
}
