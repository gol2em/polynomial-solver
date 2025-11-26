# Examples Makefile for polynomial solver
#
# Usage:
#   make <example_name>   - Compile and run specific example
#   make all              - Build all examples
#   make run-all          - Build and run all examples
#   make clean            - Remove generated executables
#   make clean-all        - Clean executables AND library build
#   make rebuild          - Rebuild library and all examples
#   make lib              - Build the library only
#   make cmake            - Reconfigure CMake (use CMAKE_OPTS for options)
#
# Examples:
#   make simple_cubic
#   make multiplicity_1d
#   make wilkinson_1d
#   make circle_ellipse
#
# Advanced:
#   make cmake CMAKE_OPTS="-DENABLE_GEOMETRY_DUMP=OFF"
#   make rebuild

# Project paths
PROJECT_ROOT := ..
BUILD_DIR := $(PROJECT_ROOT)/build
BIN_DIR := $(BUILD_DIR)/bin
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

# Example sources and targets
SOURCES := simple_cubic.cpp multiplicity_1d_roots.cpp wilkinson_1d_roots.cpp circle_ellipse_intersection.cpp
EXAMPLES := simple_cubic multiplicity_1d wilkinson_1d circle_ellipse

# Default target: show help
.PHONY: help
help:
	@echo "Examples Makefile - Build and run polynomial solver examples"
	@echo ""
	@echo "Usage:"
	@echo "  make <example>    Build and run specific example"
	@echo "  make all          Build all examples"
	@echo "  make run-all      Build and run all examples"
	@echo "  make clean        Remove local executables"
	@echo "  make clean-all    Clean executables AND library build"
	@echo "  make rebuild      Rebuild library and all examples"
	@echo "  make lib          Build the library only"
	@echo "  make cmake        Reconfigure CMake (use CMAKE_OPTS for options)"
	@echo ""
	@echo "Available examples:"
	@echo "  simple_cubic      - Simple cubic polynomial (x-0.2)(x-0.5)(x-0.8)"
	@echo "  multiplicity_1d   - High-multiplicity root (x-0.2)(x-0.6)^6"
	@echo "  wilkinson_1d      - Wilkinson polynomial (ill-conditioned)"
	@echo "  circle_ellipse    - 2D circle-ellipse intersection"
	@echo ""
	@echo "Advanced:"
	@echo "  make cmake CMAKE_OPTS=\"-DENABLE_GEOMETRY_DUMP=OFF\""
	@echo "  make rebuild"
	@echo ""
	@echo "Note: This Makefile creates local executables for quick testing."
	@echo "      Use CMake build for production: cmake --build ../build"

# Build all examples
.PHONY: all
all: $(EXAMPLES)

# Run all examples
.PHONY: run-all
run-all: all
	@echo "========================================="
	@echo "Running all examples..."
	@echo "========================================="
	@for ex in $(EXAMPLES); do \
		echo ""; \
		echo "--- Running $$ex ---"; \
		./$$ex | head -30; \
		echo ""; \
	done

# Individual example targets
.PHONY: simple_cubic
simple_cubic: simple_cubic_exe
	@echo "Running simple_cubic..."
	@./simple_cubic_exe

.PHONY: multiplicity_1d
multiplicity_1d: multiplicity_1d_exe
	@echo "Running multiplicity_1d..."
	@./multiplicity_1d_exe

.PHONY: wilkinson_1d
wilkinson_1d: wilkinson_1d_exe
	@echo "Running wilkinson_1d..."
	@./wilkinson_1d_exe

.PHONY: circle_ellipse
circle_ellipse: circle_ellipse_exe
	@echo "Running circle_ellipse..."
	@./circle_ellipse_exe

# Compile executables
simple_cubic_exe: simple_cubic.cpp $(LIBRARY)
	@echo "Compiling simple_cubic..."
	$(CXX) $(CXXFLAGS) simple_cubic.cpp -o simple_cubic_exe $(LDFLAGS)

multiplicity_1d_exe: multiplicity_1d_roots.cpp $(LIBRARY)
	@echo "Compiling multiplicity_1d..."
	$(CXX) $(CXXFLAGS) multiplicity_1d_roots.cpp -o multiplicity_1d_exe $(LDFLAGS)

wilkinson_1d_exe: wilkinson_1d_roots.cpp $(LIBRARY)
	@echo "Compiling wilkinson_1d..."
	$(CXX) $(CXXFLAGS) wilkinson_1d_roots.cpp -o wilkinson_1d_exe $(LDFLAGS)

circle_ellipse_exe: circle_ellipse_intersection.cpp $(LIBRARY)
	@echo "Compiling circle_ellipse..."
	$(CXX) $(CXXFLAGS) circle_ellipse_intersection.cpp -o circle_ellipse_exe $(LDFLAGS)

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
	@echo "Recompiling all examples..."
	@$(MAKE) clean
	@$(MAKE) all

# Clean local executables
.PHONY: clean
clean:
	@echo "Cleaning local example executables..."
	@rm -f simple_cubic_exe multiplicity_1d_exe wilkinson_1d_exe circle_ellipse_exe
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

