#!/bin/bash
# Master test runner with detailed reporting
# Author: Diogo Ribeiro (dfr@esmad.ipp.pt)

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test statistics
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Header
echo ""
echo -e "${BLUE}╔════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║         Navier-Stokes Solvers Test Suite                  ║${NC}"
echo -e "${BLUE}║         Comprehensive Validation & Verification            ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Check if tests are built
if [ ! -f test_solvers ] || [ ! -f test_conservation ] || [ ! -f test_taylor_green ]; then
    echo -e "${YELLOW}Tests not built. Building now...${NC}"
    make clean
    make all
    echo ""
fi

# Function to run a test
run_test() {
    local test_name=$1
    local test_exec=$2
    
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${YELLOW}Running: ${test_name}${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    
    if ./$test_exec; then
        echo -e "${GREEN}✓ ${test_name} PASSED${NC}"
        ((PASSED_TESTS++))
    else
        echo -e "${RED}✗ ${test_name} FAILED${NC}"
        ((FAILED_TESTS++))
    fi
    ((TOTAL_TESTS++))
    echo ""
}

# Run all tests
run_test "Solver Functionality Tests" "test_solvers"
run_test "Conservation Property Tests" "test_conservation"
run_test "Taylor-Green Vortex Validation" "test_taylor_green"

# Summary
echo -e "${BLUE}╔════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║                      TEST SUMMARY                          ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Total Test Suites: $TOTAL_TESTS"
echo -e "${GREEN}Passed: $PASSED_TESTS${NC}"

if [ $FAILED_TESTS -gt 0 ]; then
    echo -e "${RED}Failed: $FAILED_TESTS${NC}"
    echo ""
    echo -e "${RED}Some tests failed! Please review the output above.${NC}"
    exit 1
else
    echo "Failed: 0"
    echo ""
    echo -e "${GREEN}╔════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║              ALL TESTS PASSED! ✓✓✓                        ║${NC}"
    echo -e "${GREEN}╚════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    echo -e "${GREEN}Your Navier-Stokes solvers are validated and ready!${NC}"
    exit 0
fi
