# Research Tests

This directory contains research and debug test files that were used during development to explore algorithm behavior, compare methods, and debug specific issues. These are **not** part of the main test suite (`ctest`) but are preserved for reference and future research.

## Categories

### Research/Comparison Tools (5 files)
Tools for comparing different algorithms and methods:

- **test_multiplicity_methods.cpp** - Compare Taylor, Ostrowski, and Sturm multiplicity detection methods
- **test_convergence_comparison.cpp** - Compare convergence rates using Taylor vs Ostrowski estimates
- **test_method_comparison.cpp** - Compare root bounding methods (GraphHull vs ProjectedPolyhedral)
- **test_strategies.cpp** - Compare subdivision strategies (ContractFirst vs SubdivideFirst)
- **test_hp_conversion_benefit.cpp** - Compare conversion paths (power→HP→Bernstein vs power→Bernstein)

### Debug Tools (6 files)
Tools for debugging specific issues:

- **test_multiplicity_debug.cpp** - Debug high multiplicity detection
- **test_ostrowski_debug.cpp** - Debug Ostrowski formula
- **test_sturm_debug.cpp** - Debug Sturm sequence computation
- **test_differentiation_debug.cpp** - Debug 2D polynomial differentiation
- **test_manual_mult6.cpp** - Manual test of multiplicity-6 root
- **test_halley_mult6.cpp** - Test Halley's method for multiplicity-6 root

### Precision Research (5 files)
Tools for studying precision requirements and limits:

- **test_detection_limits.cpp** - Find limits of multiplicity detection with 64-bit precision
- **test_double_precision_limits.cpp** - Test double precision limits for multiplicity detection
- **test_precision_requirements.cpp** - Systematically determine precision requirements
- **test_high_multiplicity.cpp** - Test limits of multiplicity detection methods
- **test_convergence_rate.cpp** - Measure convergence rate of modified Newton method

### Duplicate/Redundant (4 files)
Tests that duplicate functionality already covered elsewhere:

- **test_precision_escalation.cpp** - Duplicate of auto_precision_lift
- **test_hp_refiner_direct.cpp** - Covered by test_result_refiner_hp
- **test_auto_precision_lift.cpp** - Feature works, covered by workflow test
- **test_double_ostrowski.cpp** - Research on Ostrowski in double precision

## Building These Tests

These tests are **not** built by default. To build them, you would need to add them back to `tests/CMakeLists.txt` or create a separate CMakeLists.txt in this directory.

## Key Findings from Research

### Multiplicity Detection
- **Taylor series ratio test**: Works well for m ≤ 4 in double precision
- **Ostrowski method**: More robust, works for m ≤ 6 in double precision
- **Sturm sequence**: Most robust but computationally expensive

### Precision Requirements
- Double precision: Reliable for m ≤ 4
- 64-bit HP: Reliable for m ≤ 6
- Higher multiplicities require more precision

### Conversion Paths
- **Finding**: Power-basis differentiation is better than Bernstein conversion for HP
- **Reason**: Avoids Bernstein conversion which loses precision for high-degree polynomials

### Root Bounding Methods
- **GraphHull**: Exact for linear functions, good for low-degree polynomials
- **ProjectedPolyhedral**: More robust for general cases, handles degeneracies better

## Notes

These files are preserved for:
1. Future research on algorithm improvements
2. Reference for understanding design decisions
3. Benchmarking new methods against old approaches
4. Debugging similar issues in the future

They are **not maintained** and may not compile with the current codebase without modifications.

