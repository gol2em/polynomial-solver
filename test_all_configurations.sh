#!/bin/bash
# Test all high-precision configurations

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test counter
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Function to print section header
print_header() {
    echo ""
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
    echo ""
}

# Function to run a test configuration
test_config() {
    local config_name="$1"
    shift
    local cmake_flags="$@"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    echo -e "${YELLOW}Testing: ${config_name}${NC}"
    echo "CMake flags: ${cmake_flags}"
    
    # Clean build directory
    cd build
    rm -rf *
    
    # Configure
    if ! cmake .. ${cmake_flags} > /dev/null 2>&1; then
        echo -e "${RED}✗ Configuration failed${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        cd ..
        return 1
    fi
    
    # Build test
    if ! make test_high_precision_types -j$(nproc) > /dev/null 2>&1; then
        echo -e "${RED}✗ Build failed${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        cd ..
        return 1
    fi
    
    # Run test
    if ./bin/test_high_precision_types > /tmp/test_output.txt 2>&1; then
        echo -e "${GREEN}✓ Test passed${NC}"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        cd ..
        return 0
    else
        echo -e "${RED}✗ Test failed${NC}"
        echo "Output:"
        cat /tmp/test_output.txt
        FAILED_TESTS=$((FAILED_TESTS + 1))
        cd ..
        return 1
    fi
}

# Main test suite
cd "$(dirname "$0")"

print_header "High-Precision Configuration Test Suite"

# Test 1: MPFR backend (default, all libraries)
test_config "MPFR Backend (Default)" \
    -DENABLE_HIGH_PRECISION=ON

# Test 2: cpp_dec_float backend (Boost only)
test_config "cpp_dec_float Backend (Boost only)" \
    -DENABLE_HIGH_PRECISION=ON \
    -DUSE_MPFR=OFF \
    -DUSE_GMP=OFF

# Test 3: Quadmath backend (no Boost)
test_config "Quadmath Backend (Standalone)" \
    -DENABLE_HIGH_PRECISION=ON \
    -DUSE_BOOST=OFF \
    -DUSE_MPFR=OFF \
    -DUSE_GMP=OFF

# Test 4: MPFR with templates disabled
test_config "MPFR Backend (No Templates)" \
    -DENABLE_HIGH_PRECISION=ON \
    -DDISABLE_TEMPLATES=ON

# Test 5: cpp_dec_float with templates disabled
test_config "cpp_dec_float Backend (No Templates)" \
    -DENABLE_HIGH_PRECISION=ON \
    -DUSE_MPFR=OFF \
    -DUSE_GMP=OFF \
    -DDISABLE_TEMPLATES=ON

# Test 6: Quadmath with templates disabled
test_config "Quadmath Backend (No Templates)" \
    -DENABLE_HIGH_PRECISION=ON \
    -DUSE_BOOST=OFF \
    -DUSE_MPFR=OFF \
    -DUSE_GMP=OFF \
    -DDISABLE_TEMPLATES=ON

# Test 7: MPFR without quadmath
test_config "MPFR Backend (No Quadmath)" \
    -DENABLE_HIGH_PRECISION=ON \
    -DUSE_QUADMATH=OFF

# Test 8: cpp_dec_float without quadmath
test_config "cpp_dec_float Backend (No Quadmath)" \
    -DENABLE_HIGH_PRECISION=ON \
    -DUSE_MPFR=OFF \
    -DUSE_GMP=OFF \
    -DUSE_QUADMATH=OFF

# Test 9: Quadmath only (no other libraries)
test_config "Quadmath Only (Minimal)" \
    -DENABLE_HIGH_PRECISION=ON \
    -DUSE_BOOST=OFF \
    -DUSE_MPFR=OFF \
    -DUSE_GMP=OFF \
    -DUSE_QUADMATH=ON

# Test 10: High-precision disabled (should skip test)
print_header "Testing with High-Precision Disabled"
TOTAL_TESTS=$((TOTAL_TESTS + 1))
cd build
rm -rf *
if cmake .. -DENABLE_HIGH_PRECISION=OFF > /dev/null 2>&1; then
    if make test_high_precision_types -j$(nproc) > /dev/null 2>&1; then
        if ./bin/test_high_precision_types > /tmp/test_output.txt 2>&1; then
            if grep -q "not enabled" /tmp/test_output.txt; then
                echo -e "${GREEN}✓ Correctly skipped when disabled${NC}"
                PASSED_TESTS=$((PASSED_TESTS + 1))
            else
                echo -e "${RED}✗ Should have skipped test${NC}"
                FAILED_TESTS=$((FAILED_TESTS + 1))
            fi
        else
            echo -e "${GREEN}✓ Test skipped (expected)${NC}"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        fi
    else
        echo -e "${YELLOW}⊗ Test not built (expected when HP disabled)${NC}"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    fi
else
    echo -e "${RED}✗ Configuration failed${NC}"
    FAILED_TESTS=$((FAILED_TESTS + 1))
fi
cd ..

# Print summary
print_header "Test Summary"
echo "Total tests:  ${TOTAL_TESTS}"
echo -e "Passed:       ${GREEN}${PASSED_TESTS}${NC}"
echo -e "Failed:       ${RED}${FAILED_TESTS}${NC}"
echo ""

if [ ${FAILED_TESTS} -eq 0 ]; then
    echo -e "${GREEN}All tests passed! ✓${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed! ✗${NC}"
    exit 1
fi

