# Testing Documentation

## Quick Start

```bash
# Setup and run all tests
chmod +x setup_tests.sh
./setup_tests.sh
```

## Test Organization

### Test Suites

| Suite | File | Tests | Purpose | Runtime |
|-------|------|-------|---------|---------|
| **Functionality** | `test_solvers.c` | 12 | Memory, initialization, basic operations | ~2s |
| **Conservation** | `test_conservation.c` | 5 | Physical properties verification | ~10s |
| **Validation** | `test_taylor_green.c` | 4 | Analytical solution comparison | ~30s |

**Total: 21 comprehensive tests covering ~95% of code**

## Test Categories

### 1. Unit Tests (Functionality)

**Finite Difference Tests:**
- ✓ Memory allocation and deallocation
- ✓ Grid initialization (lid-driven cavity)
- ✓ Boundary condition enforcement
- ✓ Adaptive time step computation
- ✓ Residual computation

**Spectral Tests:**
- ✓ Memory allocation with FFTW
- ✓ Wavenumber grid initialization
- ✓ FFT forward/backward transforms
- ✓ 2/3 dealiasing rule
- ✓ Velocity from vorticity (Poisson solver)
- ✓ Energy and enstrophy computation

### 2. Conservation Tests

**Mass Conservation:**
- FD: max(∇·u) < 1e-6 ✓
- Spectral: max(∇·u) < 1e-12 (machine precision) ✓

**Energy Conservation:**
- Energy dissipates at correct rate
- Matches viscous theory prediction
- Monotonic decay verified

**Enstrophy Evolution:**
- Proper decay characteristics
- Physical bounds maintained

**Symmetry Preservation:**
- Taylor-Green symmetries maintained
- Numerical scheme preserves structure

### 3. Validation Tests

**Taylor-Green Vortex (Analytical Solution):**

| Test | Metric | Tolerance | Status |
|------|--------|-----------|--------|
| Short-time accuracy | L2 error < 1e-4 | t=0.1 | ✓ |
| Energy decay | < 1% error | t=1.0 | ✓ |
| Spectral convergence | Exponential | Grid refinement | ✓ |
| Mass conservation | Machine precision | All times | ✓ |

**Convergence Rates:**
```
Resolution    L2 Error      Convergence Rate
16×16         1.2e-3        -
32×32         8.5e-5        ~14.1× (exponential)
64×64         4.2e-6        ~20.2× (exponential)
```

## Running Tests

### All Tests
```bash
cd tests
make test
```

### Individual Test Suites
```bash
# Functionality tests
./test_solvers

# Conservation tests
./test_conservation

# Validation tests
./test_taylor_green
```

### With Detailed Output
```bash
cd tests
./run_all_tests.sh
```

### From Root Directory
```bash
make test
```

## Test Output Interpretation

### Passing Test
```
Running: test_spectral_taylor_green_short_time
  Testing spectral solver with Taylor-Green vortex (short time)...
    Time: 0.1000, Steps: 15, L2 error: 3.42e-05
✓ PASS: test_spectral_taylor_green_short_time
```

### Failing Test
```
Running: test_mass_conservation
  Testing mass conservation...
    Step 10: max_div=1.2e-05
✗ FAIL: Maximum divergence acceptable
  Expected: <1e-06, Got: 1.2e-05
  File: test_conservation.c, Line: 85
```

## Continuous Integration

Tests run automatically on:
- Push to `main` or `develop` branches
- Pull requests
- Weekly schedule (Sundays at midnight)

### CI Platforms
- ✓ Ubuntu 22.04 (GCC, Clang)
- ✓ macOS latest (Clang)
- ✓ Memory leak detection (Valgrind)

## Benchmarks

### Performance Targets

| Solver | Resolution | Time/Step | Target |
|--------|-----------|-----------|--------|
| FD | 32×32 | 0.02s | <0.05s |
| FD | 64×64 | 0.08s | <0.15s |
| Spectral | 64×64 | 0.06s | <0.10s |
| Spectral | 128×128 | 0.25s | <0.50s |

### Accuracy Targets

| Test | Metric | Target | Current |
|------|--------|--------|---------|
| Taylor-Green (64²) | L2 error | <1e-4 | 4.2e-6 ✓ |
| Mass conservation | max(∇·u) | <1e-10 | 1e-12 ✓ |
| Energy decay | Relative error | <1% | 0.07% ✓ |

## Debugging

### Enable Debug Mode
```bash
cd tests
make clean
make CFLAGS="-O0 -g -DDEBUG"
gdb ./test_solvers
```

### Check Memory Leaks
```bash
valgrind --leak-check=full --show-leak-kinds=all ./test_solvers
```

### Profile Performance
```bash
make CFLAGS="-O2 -pg"
./test_taylor_green
gprof test_taylor_green gmon.out > profile.txt
```

## Adding New Tests

### 1. Create Test File
```c
// tests/test_new_feature.c
#include "test_framework.h"
#include "../src/your_module.h"

int test_new_feature() {
    printf("  Testing new feature...\n");
    // Test implementation
    TEST_ASSERT(condition, "Description");
    return 1;
}

int main() {
    TEST_SUITE_START("New Feature Tests");
    RUN_TEST(test_new_feature);
    TEST_SUITE_END();
}
```

### 2. Update Makefile
```makefile
TESTS = test_solvers test_conservation test_taylor_green test_new_feature

test_new_feature: test_new_feature.c $(NEEDED_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(NEEDED_OBJS) -o $@ $(LDFLAGS)
```

### 3. Update Test Runner
```bash
# In tests/run_all_tests.sh
run_test "New Feature Tests" "test_new_feature"
```

## Test Framework API

### Available Assertions

```c
// Boolean assertions
TEST_ASSERT(condition, "message")

// Equality assertions
TEST_ASSERT_EQUAL(actual, expected, "message")
TEST_ASSERT_DOUBLE_EQUAL(actual, expected, tolerance, "message")

// Pointer assertions
TEST_ASSERT_NULL(ptr, "message")
TEST_ASSERT_NOT_NULL(ptr, "message")
```

### Test Macros

```c
RUN_TEST(test_function)         // Run single test
TEST_SUITE_START("Suite Name")  // Begin suite
TEST_SUITE_END()                // End with summary
```

## Coverage Report

```
File                               Lines    Executed    Coverage
================================================================
src/finite_difference/
  fd_memory.c                      45       45          100%
  fd_initialization.c              78       78          100%
  fd_newton_raphson.c              156      142         91%
  fd_time_advance.c                64       64          100%
  fd_analysis.c                    92       88          96%

src/spectral/
  spectral_memory.c                52       52          100%
  spectral_initialization.c        89       89          100%
  spectral_transforms.c            67       67          100%
  spectral_core.c                  134      128         95%
  spectral_time_integration.c      98       95          97%
  spectral_analysis.c              112      105         94%

Overall Coverage:                  987      953         96.6%
```

## Known Issues

None currently. All tests passing on supported platforms.

## FAQ

**Q: Tests fail with "FFTW3 not found"**  
A: Install FFTW3: `sudo apt-get install libfftw3-dev`

**Q: Tests pass but with warnings**  
A: Warnings are OK. Only failures matter for correctness.

**Q: How long should tests take?**  
A: Complete suite should finish in < 2 minutes on modern hardware.

**Q: Can I run tests in parallel?**  
A: Currently sequential only. Parallel support coming soon.

**Q: What if a test fails?**  
A: Check the error message, review recent changes, run with debug flags.

## References

- Ghia, U., et al. (1982). "High-Re solutions for incompressible flow"
- Taylor, G.I. & Green, A.E. (1937). "Mechanism of the production of small eddies"
- Canuto, C., et al. (2007). "Spectral Methods: Fundamentals in Single Domains"

---

**Maintained by:** Diogo Ribeiro (dfr@esmad.ipp.pt)  
**Last Updated:** 2025-01-16  
**Test Framework Version:** 1.0.0
