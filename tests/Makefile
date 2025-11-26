# Tests Makefile for polynomial solver
#
# Usage:
#   make <test_name>      - Compile and run specific test
#   make all              - Build all tests
#   make run-all          - Build and run all tests
#   make clean            - Remove generated executables
#   make clean-all        - Clean executables AND library build
#   make rebuild          - Rebuild library and all tests
#   make lib              - Build the library only
#   make cmake            - Reconfigure CMake (use CMAKE_OPTS for options)
#
# Examples:
#   make test_polynomial_conversion
#   make test_result_refiner
#   make test_solver_subdivision
#
# Advanced:
#   make cmake CMAKE_OPTS="-DENABLE_GEOMETRY_DUMP=OFF"
#   make rebuild

# Project paths
PROJECT_ROOT := ..
BUILD_DIR := $(PROJECT_ROOT)/build
LIB_DIR := $(BUILD_DIR)/lib
INCLUDE_DIR := $(PROJECT_ROOT)/include

# Library
LIBRARY := $(LIB_DIR)/libpolynomial_solver.a

# Compiler settings
CXX := g++
CXXFLAGS := -std=c++11 -Wall -I$(INCLUDE_DIR)
LDFLAGS := -L$(LIB_DIR) -lpolynomial_solver -lm

# CMake options (can be overridden)
CMAKE_OPTS ?=

# Find all test sources
SOURCES := $(wildcard test_*.cpp)
TESTS := $(SOURCES:.cpp=_exe)

# Default target: show help
.PHONY: help
help:
	@echo "Tests Makefile - Build and run polynomial solver tests"
	@echo ""
	@echo "Usage:"
	@echo "  make <test>       Build and run specific test (without .cpp extension)"
	@echo "  make all          Build all tests"
	@echo "  make run-all      Build and run all tests"
	@echo "  make clean        Remove local executables"
	@echo "  make clean-all    Clean executables AND library build"
	@echo "  make rebuild      Rebuild library and all tests"
	@echo "  make lib          Build the library only"
	@echo "  make cmake        Reconfigure CMake (use CMAKE_OPTS for options)"
	@echo "  make list         List all available tests"
	@echo ""
	@echo "Examples:"
	@echo "  make test_polynomial_conversion"
	@echo "  make test_result_refiner"
	@echo "  make test_solver_subdivision"
	@echo ""
	@echo "Advanced:"
	@echo "  make cmake CMAKE_OPTS=\"-DENABLE_GEOMETRY_DUMP=OFF\""
	@echo "  make rebuild"
	@echo ""
	@echo "Note: This Makefile creates local executables for quick testing."
	@echo "      Use CMake build for production: cmake --build ../build"

# List all tests
.PHONY: list
list:
	@echo "Available tests:"
	@for src in $(SOURCES); do \
		name=$${src%.cpp}; \
		echo "  $$name"; \
	done

# Build all tests
.PHONY: all
all: $(TESTS)

# Run all tests
.PHONY: run-all
run-all: all
	@echo "========================================="
	@echo "Running all tests..."
	@echo "========================================="
	@for test in $(TESTS); do \
		echo ""; \
		echo "--- Running $$test ---"; \
		./$$test 2>&1 | head -50; \
		echo ""; \
	done

# Pattern rule: compile any test_*.cpp to test_*_exe
%_exe: %.cpp $(LIBRARY)
	@echo "Compiling $< ..."
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)
	@echo "✅ Created executable: $@"

# Pattern rule: run any test by name (without .cpp)
.PHONY: test_%
test_%: test_%_exe
	@echo "Running $< ..."
	@./$<

# Build the library
.PHONY: lib
lib:
	@echo "Building polynomial_solver library..."
	@cd $(PROJECT_ROOT) && cmake --build build --target polynomial_solver
	@echo "✅ Library built: $(LIBRARY)"

# Ensure library exists
$(LIBRARY):
	@echo "Library not found. Building..."
	@$(MAKE) lib

# Rebuild everything
.PHONY: rebuild
rebuild:
	@echo "Rebuilding library..."
	@cd $(PROJECT_ROOT) && cmake --build build --target polynomial_solver
	@echo "Recompiling all tests..."
	@$(MAKE) clean
	@$(MAKE) all

# Clean local executables
.PHONY: clean
clean:
	@echo "Cleaning local test executables..."
	@rm -f $(TESTS)
	@echo "✅ Clean complete"
	@echo "Note: CMake build directory not affected. Use 'make clean-all' for full clean."

# Clean everything including library build
.PHONY: clean-all
clean-all: clean
	@echo "Cleaning library build..."
	@cd $(PROJECT_ROOT) && rm -rf build
	@echo "✅ Full clean complete"

# Reconfigure CMake with options
.PHONY: cmake
cmake:
	@echo "Reconfiguring CMake..."
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake $(CMAKE_OPTS) ..
	@echo "✅ CMake configuration complete"

