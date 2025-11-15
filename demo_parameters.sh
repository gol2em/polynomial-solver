#!/bin/bash

# Demonstration of command-line parameters for polynomial solver

echo "=========================================="
echo "Polynomial Solver Parameter Demo"
echo "=========================================="
echo ""

echo "1. Show help message:"
echo "   ./bin/polynomial_solver_app --help"
echo ""
./bin/polynomial_solver_app --help
echo ""

echo "=========================================="
echo ""

echo "2. Run with default parameters:"
echo "   ./bin/polynomial_solver_app"
echo ""
./bin/polynomial_solver_app
echo ""

echo "=========================================="
echo ""

echo "3. Run with custom parameters:"
echo "   ./bin/polynomial_solver_app --debug --crit 0.9 --tolerance 1e-4 --max-depth 25"
echo ""
./bin/polynomial_solver_app --debug --crit 0.9 --tolerance 1e-4 --max-depth 25
echo ""

echo "=========================================="
echo ""

echo "4. Run tests without debug output:"
echo "   ./bin/test_solver_subdivision"
echo ""
./bin/test_solver_subdivision
echo ""

echo "=========================================="
echo ""

echo "5. Run tests with debug output (showing first 50 lines):"
echo "   ./bin/test_solver_subdivision --debug | head -50"
echo ""
./bin/test_solver_subdivision --debug 2>&1 | head -50
echo ""

echo "=========================================="
echo "Demo complete!"
echo "=========================================="

