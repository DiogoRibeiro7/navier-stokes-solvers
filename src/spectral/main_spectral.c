#include "ns_spectral_solver.h"
#include "../common/performance.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char *argv[]) {
    ns_perf_configure_threading("Spectral solver");

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
    ns_perf_counter wall_timer;
    ns_perf_counter_start(&wall_timer, 0, 0);
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
    
    wall_timer.bytes_moved = (size_t)step * (size_t)nx * (size_t)ny * sizeof(double) * 8;
    ns_perf_counter_stop(&wall_timer);
    double cpu_time = ns_perf_seconds(&wall_timer.start, &wall_timer.end);
    
    printf("\n=== SIMULATION COMPLETED ===\n");
    printf("Total steps: %d\n", step);
    printf("Final time: %.6f\n", data->t);
    printf("Final energy: %.6e\n", ns_spectral_compute_energy(data));
    printf("Final enstrophy: %.6e\n", ns_spectral_compute_enstrophy(data));
    printf("CPU time: %.2f seconds\n", cpu_time);
    printf("Steps/second: %.1f\n\n", step / cpu_time);
    printf("Resident memory: %.2f MiB (peak %.2f MiB)\n",
           ns_perf_resident_memory() / (1024.0 * 1024.0),
           ns_perf_peak_resident_memory() / (1024.0 * 1024.0));
    ns_perf_print_summary("Spectral solver", &wall_timer);
    
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
