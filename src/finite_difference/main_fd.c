#include "ns_fd_solver.h"
#include <time.h>

int main(int argc, char *argv[]) {
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
    clock_t start = clock();
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
    
    clock_t end = clock();
    double cpu_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    printf("\nSimulation completed!\n");
    printf("Total steps: %d\n", step);
    printf("Final time: %.4f\n", data->t);
    printf("CPU time: %.2f seconds\n", cpu_time);
    
    ns_fd_output_solution(data, "fd_final_solution.dat");
    ns_fd_free(data);
    
    return 0;
}
