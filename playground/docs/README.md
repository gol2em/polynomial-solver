# Research and Analysis Documentation

This directory contains research documents, analysis, and design explorations from the development of the polynomial solver.

## Purpose

These documents are **not user-facing documentation**. They are:
- Research notes and background information
- Analysis of different approaches and methods
- Design explorations and comparisons
- Development progress tracking
- Experimental results and findings

For user documentation, see the main `docs/` directory.

## Contents

### Background Research
- `polynomial_system_solvers_classification.md` - Classification of polynomial system solvers
- `univariate_polynomial_solvers.md` - Survey of 1D polynomial solvers
- `pp_vs_homotopy_comparison.md` - Comparison of PP method vs homotopy continuation
- `deflation_method.md` - Deflation methods for multiple roots
- `library_algorithm_summary.md` - Summary of algorithms in other libraries

### Multiplicity Detection Research
- `multiplicity_detection_integration.md` - Integration of Ostrowski method
- `multiple_roots_precision_handling.md` - Precision requirements for multiple roots
- `MULTIPLICITY_DETECTION_COMPARISON.md` - Comparison of detection methods
- `MULTIPLE_ROOT_HANDLING_SUMMARY.md` - Complete implementation summary

### Convergence Analysis
- `CONVERGENCE_ANALYSIS.md` - Convergence behavior analysis
- `CONVERGENCE_CRITERIA_ANALYSIS.md` - Analysis of convergence criteria
- `ERROR_BOUND_CONVERGENCE.md` - Error bound computation analysis
- `DETECTION_LIMITS.md` - Limits of detection methods

### Precision Requirements
- `PRECISION_REQUIREMENTS.md` - Systematic precision requirements analysis
- `CONDITIONING_AND_PRECISION.md` - Condition numbers and precision needs
- `1D_POLYNOMIAL_WORKFLOW.md` - Recommended workflow for 1D polynomials

### High-Precision Design (Tier 3)
- `HIGH_PRECISION_MINIMAL_DEPENDENCIES.md` - Minimal dependency design
- `HIGH_PRECISION_MULTIPLE_TYPES.md` - Multiple type support design
- `RUNTIME_PRECISION_DESIGN.md` - Runtime precision control design
- `TEMPLATE_CONTROL_LOGIC.md` - Template control logic explanation

### Configuration and Testing
- `CONFIGURATION_TEST_RESULTS.md` - Configuration testing results

## Organization

Documents are organized by topic rather than chronologically. Some documents may reference features or designs that are not yet implemented (e.g., Tier 3 templates).

## For Developers

If you're working on the polynomial solver:
1. Check these documents for background on design decisions
2. Review analysis documents before implementing new features
3. Add new research/analysis documents here (not in main `docs/`)
4. Update this README when adding new documents

## For Users

If you're using the polynomial solver library, you probably don't need these documents. See the main documentation in `docs/` instead:
- `docs/HIGH_PRECISION.md` - High-precision arithmetic guide
- `docs/PARAMETERS.md` - Configuration parameters
- `docs/WORKFLOW_DESIGN.md` - Recommended workflows

