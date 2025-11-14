#!/bin/bash

# Polynomial Solver - Build Script
# This script builds the project and optionally runs tests

set -e

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Default options
BUILD_TYPE="Release"
RUN_TESTS=false
CLEAN_BUILD=false
JOBS=4
VERBOSE=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        --test)
            RUN_TESTS=true
            shift
            ;;
        --clean)
            CLEAN_BUILD=true
            shift
            ;;
        --jobs|-j)
            JOBS="$2"
            shift 2
            ;;
        --verbose|-v)
            VERBOSE=true
            shift
            ;;
        --help|-h)
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --debug          Build in Debug mode (default: Release)"
            echo "  --test           Run tests after building"
            echo "  --clean          Clean build directory before building"
            echo "  --jobs, -j N     Use N parallel jobs (default: 4)"
            echo "  --verbose, -v    Verbose build output"
            echo "  --help, -h       Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0                    # Build in Release mode"
            echo "  $0 --debug --test     # Build in Debug mode and run tests"
            echo "  $0 --clean -j 8       # Clean build with 8 parallel jobs"
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Polynomial Solver - Build Script${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Check if we're in the right directory
if [ ! -f "CMakeLists.txt" ]; then
    echo -e "${RED}Error: CMakeLists.txt not found!${NC}"
    echo -e "${RED}Please run this script from the project root directory.${NC}"
    exit 1
fi

# Clean build if requested
if [ "$CLEAN_BUILD" = true ]; then
    echo -e "${YELLOW}Cleaning build directory...${NC}"
    rm -rf build
    echo -e "${GREEN}✓ Build directory cleaned${NC}"
    echo ""
fi

# Create build directory
echo -e "${BLUE}Creating build directory...${NC}"
mkdir -p build
cd build

# Configure with CMake
echo -e "${BLUE}Configuring with CMake...${NC}"
echo -e "  Build type: ${GREEN}$BUILD_TYPE${NC}"
echo -e "  Parallel jobs: ${GREEN}$JOBS${NC}"
echo ""

if [ "$VERBOSE" = true ]; then
    cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE ..
else
    cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE .. > /dev/null
fi

echo -e "${GREEN}✓ Configuration complete${NC}"
echo ""

# Build
echo -e "${BLUE}Building project...${NC}"
START_TIME=$(date +%s)

if [ "$VERBOSE" = true ]; then
    make -j$JOBS
else
    make -j$JOBS 2>&1 | grep -E "Built target|Linking|error:|warning:" || true
fi

END_TIME=$(date +%s)
BUILD_TIME=$((END_TIME - START_TIME))

echo -e "${GREEN}✓ Build complete in ${BUILD_TIME}s${NC}"
echo ""

# List built executables
echo -e "${BLUE}Built executables:${NC}"
if [ -d "bin" ]; then
    ls -lh bin/ | tail -n +2 | awk '{printf "  %-30s %10s\n", $9, $5}'
fi
echo ""

# Run tests if requested
if [ "$RUN_TESTS" = true ]; then
    echo -e "${BLUE}Running tests...${NC}"
    echo ""
    ctest --output-on-failure
    echo ""
    echo -e "${GREEN}✓ All tests completed${NC}"
    echo ""
fi

# Summary
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Build Summary${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "${GREEN}✓ Project built successfully!${NC}"
echo ""
echo -e "Build artifacts location: ${BLUE}$(pwd)${NC}"
echo ""
echo -e "Next steps:"
echo -e "  Run tests:        ${BLUE}cd build && ctest${NC}"
echo -e "  Run comparison:   ${BLUE}./build/bin/test_method_comparison${NC}"
echo -e "  Run application:  ${BLUE}./build/bin/polynomial_solver_app --help${NC}"
echo ""

