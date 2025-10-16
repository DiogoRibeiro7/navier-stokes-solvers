# Test Suite Makefile
# Author: Diogo Ribeiro (dfr@esmad.ipp.pt)

CC = gcc
CFLAGS = -O2 -Wall -Wextra -std=c11 -g
INCLUDES = -I.. -I../include -I../src/finite_difference -I../src/spectral
LDFLAGS = -lm
FFTW_FLAGS = -lfftw3

# Directories
SRC_FD = ../src/finite_difference
SRC_SPECTRAL = ../src/spectral
OBJ_DIR = obj

# Create directories
$(shell mkdir -p $(OBJ_DIR))

# FD object files (excluding main)
FD_OBJS = $(OBJ_DIR)/fd_memory.o \
          $(OBJ_DIR)/fd_initialization.o \
          $(OBJ_DIR)/fd_newton_raphson.o \
          $(OBJ_DIR)/fd_time_advance.o \
          $(OBJ_DIR)/fd_analysis.o

# Spectral object files (excluding main)
SPECTRAL_OBJS = $(OBJ_DIR)/spectral_memory.o \
                $(OBJ_DIR)/spectral_initialization.o \
                $(OBJ_DIR)/spectral_transforms.o \
                $(OBJ_DIR)/spectral_core.o \
                $(OBJ_DIR)/spectral_time_integration.o \
                $(OBJ_DIR)/spectral_analysis.o \
                $(OBJ_DIR)/spectral_output.o

# Test executables
TESTS = test_solvers test_conservation test_taylor_green

.PHONY: all clean test run-tests

all: $(TESTS)

# Build all tests and run them
test: all run-tests

run-tests:
	@echo ""
	@echo "Running all tests..."
	@echo "===================="
	@./test_solvers
	@./test_conservation
	@./test_taylor_green
	@echo ""
	@echo "All tests completed!"

# Individual test targets
test_solvers: test_solvers.c $(FD_OBJS) $(SPECTRAL_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(FD_OBJS) $(SPECTRAL_OBJS) -o $@ $(LDFLAGS) $(FFTW_FLAGS)

test_conservation: test_conservation.c $(FD_OBJS) $(SPECTRAL_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(FD_OBJS) $(SPECTRAL_OBJS) -o $@ $(LDFLAGS) $(FFTW_FLAGS)

test_taylor_green: test_taylor_green.c $(SPECTRAL_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(SPECTRAL_OBJS) -o $@ $(LDFLAGS) $(FFTW_FLAGS)

# Compile FD object files
$(OBJ_DIR)/fd_%.o: $(SRC_FD)/fd_%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Compile spectral object files
$(OBJ_DIR)/spectral_%.o: $(SRC_SPECTRAL)/spectral_%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Clean
clean:
	rm -rf $(OBJ_DIR)
	rm -f $(TESTS)
	rm -f *.dat
	@echo "Cleaned test artifacts"

# Help
help:
	@echo "Test Suite Build System"
	@echo "======================="
	@echo "Targets:"
	@echo "  all        - Build all tests"
	@echo "  test       - Build and run all tests"
	@echo "  run-tests  - Run all tests"
	@echo "  clean      - Remove build artifacts"
	@echo ""
	@echo "Individual tests:"
	@echo "  ./test_solvers      - Basic functionality tests"
	@echo "  ./test_conservation - Conservation property tests"
	@echo "  ./test_taylor_green - Taylor-Green vortex validation"
