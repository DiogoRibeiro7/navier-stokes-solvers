#include "ns_fd_solver.h"
#include "../common/performance.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char *argv[]) {
    ns_perf_configure_threading("Finite Difference solver");

    printf("2D Incompressible Navier-Stokes Solver\n");
    printf("Finite Difference Method with Newton-Raphson\n\n");
    
    // Problem parameters
    int nx = 32, ny = 32;
    double L = 1.0, H = 1.0;
    double Re = 100.0;
    double t_final = 1.0;
    
    printf("Grid: %dx%d, Re=%.1f, Domain=[%.1fx%.1f]\n", nx, ny, Re, L, H);
    
    // Allocate solver
    NSFiniteDiffData *data = ns_fd_allocate(nx, ny, L, H, Re);
    if (!data) {
        printf("Memory allocation failed!\n");
        return 1;
    }
    
    // Initialize flow
    ns_fd_initialize_lid_cavity(data);
    
    // Time integration
    ns_perf_counter wall_timer;
    ns_perf_counter_start(&wall_timer, 0, 0);
    int step = 0;
    
    while (data->t < t_final) {
        printf("\nTime step %d:\n", ++step);
        
        SolverStatus status = ns_fd_advance_timestep(data);
        if (status == SOLVER_DIVERGED) {
            printf("Solver diverged!\n");
            break;
        }
        
        ns_fd_analyze_convergence(data, step);
        
        if (step % 50 == 0) {
            char filename[100];
            sprintf(filename, "fd_solution_step_%04d.dat", step);
            ns_fd_output_solution(data, filename);
        }
    }
    
    wall_timer.bytes_moved = (size_t)step * (size_t)nx * (size_t)ny * sizeof(double) * 6;
    ns_perf_counter_stop(&wall_timer);
    double cpu_time = ns_perf_seconds(&wall_timer.start, &wall_timer.end);
    
    printf("\nSimulation completed!\n");
    printf("Total steps: %d\n", step);
    printf("Final time: %.4f\n", data->t);
    printf("CPU time: %.2f seconds\n", cpu_time);
    printf("Resident memory: %.2f MiB (peak %.2f MiB)\n",
           ns_perf_resident_memory() / (1024.0 * 1024.0),
           ns_perf_peak_resident_memory() / (1024.0 * 1024.0));
    ns_perf_print_summary("FD solver", &wall_timer);
    
    ns_fd_output_solution(data, "fd_final_solution.dat");
    ns_fd_free(data);
    
    return 0;
}
