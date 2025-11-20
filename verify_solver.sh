#!/bin/bash

# Standard verification script for polynomial solver
# Tests circle-ellipse intersection with all strategies and generates visualizations

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=========================================="
echo "Polynomial Solver - Standard Verification"
echo "=========================================="
echo ""

# Check if build exists
if [ ! -d "build" ]; then
    echo "Error: build/ directory not found. Please run ./build.sh first."
    exit 1
fi

# Check if test_strategies exists
if [ ! -f "build/bin/test_strategies" ]; then
    echo "Error: test_strategies not found. Please run ./build.sh first."
    exit 1
fi

# Create necessary directories
mkdir -p dumps visualizations

# Clean old dumps and visualizations
echo "Cleaning old dumps and visualizations..."
rm -f dumps/strategy_*.txt
rm -rf visualizations/viz_*

# Run the test
echo ""
echo "Running circle-ellipse intersection test with all strategies..."
echo "--------------------------------------------------------------"
./build/bin/test_strategies

# Check if Python venv exists
if [ ! -d ".venv" ]; then
    echo ""
    echo "Warning: Python virtual environment not found."
    echo "Skipping visualization. To enable visualization:"
    echo "  1. Create venv: uv venv .venv"
    echo "  2. Install deps: uv pip install numpy matplotlib"
    exit 0
fi

# Check if visualization tool exists
if [ ! -f "tools/visualize_solver.py" ]; then
    echo ""
    echo "Warning: Visualization tool not found at tools/visualize_solver.py"
    echo "Skipping visualization."
    exit 0
fi

# Activate venv and run visualizations
echo ""
echo "Generating visualizations..."
echo "----------------------------"

source .venv/bin/activate

strategies=("ContractFirst" "SubdivideFirst" "Simultaneous")

for strategy in "${strategies[@]}"; do
    dump_file="dumps/strategy_${strategy}_geometry.txt"
    output_dir="visualizations/viz_${strategy}"
    
    if [ -f "$dump_file" ]; then
        echo "Visualizing $strategy strategy..."
        python tools/visualize_solver.py "$dump_file" --output-dir "$output_dir"
        
        # Count generated images
        num_images=$(find "$output_dir" -name "*.png" 2>/dev/null | wc -l)
        echo "  Generated $num_images visualization(s)"
    else
        echo "Warning: Dump file not found: $dump_file"
    fi
done

echo ""
echo "=========================================="
echo "Verification Complete!"
echo "=========================================="
echo ""
echo "Results:"
echo "  Dumps:          dumps/strategy_*.txt"
echo "  Visualizations: visualizations/viz_*/"
echo ""
echo "To view visualizations, open PNG files in visualizations/viz_*/"

