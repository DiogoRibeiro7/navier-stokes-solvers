#!/bin/bash
# Setup and run complete test suite
# Author: Diogo Ribeiro (dfr@esmad.ipp.pt)

set -e

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}"
echo "╔════════════════════════════════════════════════════════════╗"
echo "║     Navier-Stokes Solvers - Test Suite Setup              ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

# Check dependencies
echo -e "${YELLOW}Checking dependencies...${NC}"

if ! command -v gcc &> /dev/null; then
    echo "ERROR: GCC not found. Please install gcc."
    exit 1
fi
echo "✓ GCC found: $(gcc --version | head -n1)"

if ! ldconfig -p | grep -q libfftw3; then
    echo "WARNING: FFTW3 not found. Spectral tests may fail."
    echo "  Install with: sudo apt-get install libfftw3-dev"
else
    echo "✓ FFTW3 found"
fi

# Create directory structure
echo ""
echo -e "${YELLOW}Creating directory structure...${NC}"
mkdir -p tests/obj
chmod +x tests/run_all_tests.sh
echo "✓ Directories created"

# Build main solvers
echo ""
echo -e "${YELLOW}Building main solvers...${NC}"
make clean
make all
echo "✓ Solvers built"

# Build and run tests
echo ""
echo -e "${YELLOW}Building test suite...${NC}"
cd tests
make clean
make all
echo "✓ Tests built"

# Run tests
echo ""
echo -e "${BLUE}════════════════════════════════════════════════════════════${NC}"
echo -e "${BLUE}Running complete test suite...${NC}"
echo -e "${BLUE}════════════════════════════════════════════════════════════${NC}"
echo ""

./run_all_tests.sh
TEST_RESULT=$?

cd ..

# Summary
echo ""
if [ $TEST_RESULT -eq 0 ]; then
    echo -e "${GREEN}════════════════════════════════════════════════════════════${NC}"
    echo -e "${GREEN}SUCCESS! All tests passed.${NC}"
    echo -e "${GREEN}════════════════════════════════════════════════════════════${NC}"
    echo ""
    echo "Your Navier-Stokes solvers are:"
    echo "  ✓ Properly compiled"
    echo "  ✓ Functionally correct"
    echo "  ✓ Physically accurate"
    echo "  ✓ Ready for research use"
    echo ""
    echo "Next steps:"
    echo "  1. Run solvers: ./bin/ns_fd_solver or ./bin/ns_spectral_solver"
    echo "  2. Visualize results: python3 scripts/plot_solution.py"
    echo "  3. Read documentation: docs/"
else
    echo -e "${RED}════════════════════════════════════════════════════════════${NC}"
    echo -e "${RED}FAILURE: Some tests failed.${NC}"
    echo -e "${RED}════════════════════════════════════════════════════════════${NC}"
    echo ""
    echo "Please review the test output above and fix any issues."
    exit 1
fi
