# Navier-Stokes Solvers Build System

CC              ?= gcc
CFLAGS_BASE      = -std=c11 -Wall -Wextra -fopenmp
OPT_FLAGS      ?= -O3 -march=native
PROFILE_FLAGS  ?=
CFLAGS          = $(CFLAGS_BASE) $(OPT_FLAGS) $(PROFILE_FLAGS)
INCLUDES        = -I./include -I./src/common
LDFLAGS_BASE    = -lm
FFTW_FLAGS      = -lfftw3 -lfftw3_omp -lfftw3_threads
LDFLAGS         = $(LDFLAGS_BASE) $(PROFILE_FLAGS)

BIN_DIR         = bin
OBJ_DIR         = obj
BENCH_DIR       = benchmarks

SRC_COMMON      = src/common
SRC_FD          = src/finite_difference
SRC_SPECTRAL    = src/spectral

COMMON_OBJS     = $(OBJ_DIR)/common/performance.o

FD_CORE_OBJS    = $(OBJ_DIR)/fd/fd_memory.o \
                  $(OBJ_DIR)/fd/fd_initialization.o \
                  $(OBJ_DIR)/fd/fd_fd_newton_raphson.o \
                  $(OBJ_DIR)/fd/fd_time_advance.o \
                  $(OBJ_DIR)/fd/fd_analysis.o
FD_OBJS         = $(FD_CORE_OBJS) $(OBJ_DIR)/fd/main_fd.o

SPECTRAL_CORE_OBJS = $(OBJ_DIR)/spectral/spectral_memory.o \
                     $(OBJ_DIR)/spectral/spectral_initialization.o \
                     $(OBJ_DIR)/spectral/spectral_transforms.o \
                     $(OBJ_DIR)/spectral/spectral_core.o \
                     $(OBJ_DIR)/spectral/spectral_time_integration.o \
                     $(OBJ_DIR)/spectral/spectral_analysis.o \
                     $(OBJ_DIR)/spectral/spectral_output.o
SPECTRAL_OBJS      = $(SPECTRAL_CORE_OBJS) $(OBJ_DIR)/spectral/main_spectral.o

BENCH_OBJS         = $(OBJ_DIR)/bench/performance_benchmark.o

.PHONY: all clean fd spectral performance bench pgo test help

all: fd spectral

fd: $(BIN_DIR)/ns_fd_solver

spectral: $(BIN_DIR)/ns_spectral_solver

performance:
	@echo "[build] Optimised performance build"
	$(MAKE) clean
	$(MAKE) OPT_FLAGS="-O3 -march=native -flto -funroll-loops -fomit-frame-pointer"

bench: $(BIN_DIR)/performance_benchmark
	@echo "[bench] Running performance benchmark..."
	$(BIN_DIR)/performance_benchmark > $(BENCH_DIR)/performance_results.csv
	@echo "[bench] Results stored in $(BENCH_DIR)/performance_results.csv"

pgo:
	@echo "[build] Profile-guided optimisation workflow"
	$(MAKE) clean
	$(MAKE) OPT_FLAGS="$(OPT_FLAGS)" PROFILE_FLAGS="-fprofile-generate"
	$(BIN_DIR)/performance_benchmark > /dev/null || true
	$(MAKE) clean
	$(MAKE) OPT_FLAGS="$(OPT_FLAGS)" PROFILE_FLAGS="-fprofile-use -fprofile-correction"

help:
	@echo "Navier-Stokes Solvers Build System"
	@echo "Targets:"
	@echo "  make            Build finite-difference and spectral solvers"
	@echo "  make fd         Build finite-difference solver only"
	@echo "  make spectral   Build spectral solver only"
	@echo "  make performance Rebuild with aggressive optimisation flags"
	@echo "  make bench      Build and run performance benchmarks"
	@echo "  make pgo        Execute profile-guided optimisation workflow"
	@echo "  make test       Build and execute regression tests"
	@echo "  make clean      Remove build artefacts"

test: all
	@$(MAKE) -C tests test

$(BIN_DIR)/ns_fd_solver: $(COMMON_OBJS) $(FD_OBJS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "[build] Built $@"

$(BIN_DIR)/ns_spectral_solver: $(COMMON_OBJS) $(SPECTRAL_OBJS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(FFTW_FLAGS)
	@echo "[build] Built $@"

$(BIN_DIR)/performance_benchmark: $(COMMON_OBJS) $(FD_CORE_OBJS) $(SPECTRAL_CORE_OBJS) $(BENCH_OBJS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(FFTW_FLAGS)
	@echo "[build] Built benchmark harness $@"

$(OBJ_DIR)/common/%.o: $(SRC_COMMON)/%.c
	@mkdir -p $(OBJ_DIR)/common
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/fd/%.o: $(SRC_FD)/%.c
	@mkdir -p $(OBJ_DIR)/fd
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/spectral/%.o: $(SRC_SPECTRAL)/%.c
	@mkdir -p $(OBJ_DIR)/spectral
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/bench/%.o: $(BENCH_DIR)/%.c
	@mkdir -p $(OBJ_DIR)/bench
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf $(BIN_DIR) $(OBJ_DIR) pgo-data *.dat benchmarks/performance_results.csv
	@echo "[clean] Removed build artefacts"
