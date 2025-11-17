#!/bin/bash

# Polynomial Solver - Prerequisite Checker
# This script checks if all required tools and libraries are installed

set -e

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Counters
PASSED=0
FAILED=0
WARNINGS=0

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Polynomial Solver - Prerequisite Check${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Function to check if a command exists
check_command() {
    local cmd=$1
    local name=$2
    local required=$3
    local min_version=$4
    
    if command -v "$cmd" &> /dev/null; then
        local version=$($cmd --version 2>&1 | head -n 1)
        echo -e "${GREEN}✓${NC} $name found: $version"
        ((PASSED++))
        return 0
    else
        if [ "$required" = "yes" ]; then
            echo -e "${RED}✗${NC} $name NOT FOUND (required)"
            ((FAILED++))
        else
            echo -e "${YELLOW}⚠${NC} $name NOT FOUND (optional)"
            ((WARNINGS++))
        fi
        return 1
    fi
}

# Function to check C++ compiler
check_cpp_compiler() {
    echo -e "\n${BLUE}Checking C++ Compiler...${NC}"
    
    if check_command "g++" "GCC C++ Compiler" "yes"; then
        local version=$(g++ -dumpversion)
        local major=$(echo $version | cut -d. -f1)
        if [ "$major" -ge 4 ]; then
            echo -e "  ${GREEN}→${NC} GCC version $version supports C++11"
        else
            echo -e "  ${YELLOW}→${NC} GCC version $version may not fully support C++11"
            ((WARNINGS++))
        fi
    fi
}

# Function to check CMake
check_cmake() {
    echo -e "\n${BLUE}Checking Build System...${NC}"
    
    if check_command "cmake" "CMake" "yes"; then
        local version=$(cmake --version | head -n 1 | grep -oP '\d+\.\d+\.\d+')
        local major=$(echo $version | cut -d. -f1)
        local minor=$(echo $version | cut -d. -f2)
        
        if [ "$major" -gt 3 ] || ([ "$major" -eq 3 ] && [ "$minor" -ge 15 ]); then
            echo -e "  ${GREEN}→${NC} CMake version $version meets minimum requirement (3.15)"
        else
            echo -e "  ${RED}→${NC} CMake version $version is below minimum requirement (3.15)"
            ((FAILED++))
        fi
    fi
    
    check_command "make" "GNU Make" "yes"
}

# Function to check Python (optional, for visualization)
check_python() {
    echo -e "\n${BLUE}Checking Python (optional for visualization)...${NC}"
    
    if check_command "python3" "Python 3" "no"; then
        # Check for numpy
        if python3 -c "import numpy" 2>/dev/null; then
            echo -e "  ${GREEN}→${NC} NumPy is installed"
        else
            echo -e "  ${YELLOW}→${NC} NumPy not found (optional for visualization)"
        fi
        
        # Check for matplotlib
        if python3 -c "import matplotlib" 2>/dev/null; then
            echo -e "  ${GREEN}→${NC} Matplotlib is installed"
        else
            echo -e "  ${YELLOW}→${NC} Matplotlib not found (optional for visualization)"
        fi
    fi
    
    check_command "uv" "UV (Python package manager)" "no"
}

# Function to check Git
check_git() {
    echo -e "\n${BLUE}Checking Version Control...${NC}"
    check_command "git" "Git" "yes"
}

# Function to check system info
check_system() {
    echo -e "\n${BLUE}System Information...${NC}"
    echo -e "${GREEN}✓${NC} OS: $(uname -s)"
    echo -e "${GREEN}✓${NC} Architecture: $(uname -m)"
    echo -e "${GREEN}✓${NC} Kernel: $(uname -r)"
    
    if [ -f /etc/os-release ]; then
        . /etc/os-release
        echo -e "${GREEN}✓${NC} Distribution: $NAME $VERSION"
    fi
}

# Run all checks
check_system
check_cpp_compiler
check_cmake
check_python
check_git

# Summary
echo -e "\n${BLUE}========================================${NC}"
echo -e "${BLUE}Summary${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "${GREEN}Passed:${NC}   $PASSED"
echo -e "${YELLOW}Warnings:${NC} $WARNINGS"
echo -e "${RED}Failed:${NC}   $FAILED"
echo ""

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}✓ All required prerequisites are met!${NC}"
    echo -e "${GREEN}  You can proceed with building the project.${NC}"
    echo ""
    echo -e "Next steps:"
    echo -e "  1. Run: ${BLUE}./build.sh${NC}"
    echo -e "  2. Or manually: ${BLUE}mkdir -p build && cd build && cmake .. && make -j4${NC}"
    echo -e "  3. Run tests: ${BLUE}cd build && ctest${NC}"
    exit 0
else
    echo -e "${RED}✗ Some required prerequisites are missing!${NC}"
    echo -e "${RED}  Please install the missing components before building.${NC}"
    echo ""
    echo -e "Installation hints:"
    echo -e "  Ubuntu/Debian: ${BLUE}sudo apt-get install build-essential cmake git${NC}"
    echo -e "  Fedora/RHEL:   ${BLUE}sudo dnf install gcc-c++ cmake git${NC}"
    echo -e "  Arch Linux:    ${BLUE}sudo pacman -S base-devel cmake git${NC}"
    exit 1
fi

