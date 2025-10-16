# Navier-Stokes Solvers Makefile
# Author: Diogo Ribeiro (dfr@esmad.ipp.pt)

CC = gcc
CFLAGS = -O3 -Wall -Wextra -std=c11 -march=native -fopenmp
INCLUDES = -I./include
LDFLAGS = -lm
FFTW_FLAGS = -lfftw3 -lfftw3_omp

# Directories
SRC_FD = src/finite_difference
SRC_SPECTRAL = src/spectral
BIN_DIR = bin
OBJ_DIR = obj

# Create directories
$(shell mkdir -p $(BIN_DIR) $(OBJ_DIR)/fd $(OBJ_DIR)/spectral)

# Finite Difference solver objects
FD_OBJS = $(OBJ_DIR)/fd/fd_memory.o \
          $(OBJ_DIR)/fd/fd_initialization.o \
          $(OBJ_DIR)/fd/fd_newton_raphson.o \
          $(OBJ_DIR)/fd/fd_time_advance.o \
          $(OBJ_DIR)/fd/fd_analysis.o \
          $(OBJ_DIR)/fd/main_fd.o

# Spectral solver objects
SPECTRAL_OBJS = $(OBJ_DIR)/spectral/spectral_memory.o \
                $(OBJ_DIR)/spectral/spectral_initialization.o \
                $(OBJ_DIR)/spectral/spectral_transforms.o \
                $(OBJ_DIR)/spectral/spectral_core.o \
                $(OBJ_DIR)/spectral/spectral_time_integration.o \
                $(OBJ_DIR)/spectral/spectral_analysis.o \
                $(OBJ_DIR)/spectral/spectral_output.o \
                $(OBJ_DIR)/spectral/main_spectral.o

# Targets
.PHONY: all clean fd spectral help

all: fd spectral

fd: $(BIN_DIR)/ns_fd_solver

spectral: $(BIN_DIR)/ns_spectral_solver

# Finite Difference solver
$(BIN_DIR)/ns_fd_solver: $(FD_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Built Finite Difference solver: $@"

$(OBJ_DIR)/fd/%.o: $(SRC_FD)/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Spectral solver
$(BIN_DIR)/ns_spectral_solver: $(SPECTRAL_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(FFTW_FLAGS)
	@echo "Built Spectral solver: $@"

$(OBJ_DIR)/spectral/%.o: $(SRC_SPECTRAL)/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Clean
clean:
	rm -rf $(BIN_DIR) $(OBJ_DIR)
	rm -f *.dat
	@echo "Cleaned build artifacts"

# Help
help:
	@echo "Navier-Stokes Solvers Build System"
	@echo "=================================="
	@echo "Targets:"
	@echo "  all       - Build both solvers (default)"
	@echo "  fd        - Build finite difference solver"
	@echo "  spectral  - Build spectral solver"
	@echo "  clean     - Remove build artifacts"
	@echo "  help      - Show this message"
	@echo ""
	@echo "Requirements:"
	@echo "  - GCC compiler with C11 support"
	@echo "  - FFTW3 library (for spectral solver)"
	@echo "  - OpenMP support (optional, for parallelization)"
	@echo ""
	@echo "Usage:"
	@echo "  make all"
	@echo "  ./bin/ns_fd_solver"
	@echo "  ./bin/ns_spectral_solver"
