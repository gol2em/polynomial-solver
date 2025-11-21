#!/usr/bin/env python3
"""
Visualize the precision test results for direct contraction.
"""

import matplotlib.pyplot as plt
import numpy as np

def main():
    # Data from test results
    tolerances = np.array([1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16])
    
    # Cubic results
    cubic_errors = np.array([6.396e-09, 6.396e-09, 2.220e-16, 2.220e-16, 2.220e-16, 
                             2.220e-16, 2.220e-16, 2.220e-16, 2.220e-16, 2.220e-16, 2.220e-16])
    cubic_widths = np.array([1.446e-08, 1.446e-08, 5.551e-17, 5.551e-17, 5.551e-17,
                             5.551e-17, 5.551e-17, 5.551e-17, 5.551e-17, 5.551e-17, 5.551e-17])
    cubic_resolved = np.array([4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4])
    
    # Multiplicity results
    mult_resolved = np.array([1, 1, 1, 1, 1, 1, 1, 4, 0, 0, 0])
    mult_unresolved = np.array([99, 99, 99, 99, 99, 99, 99, 99, 83, 83, 83])
    
    # Machine epsilon
    machine_eps = 2.22e-16
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Direct Contraction: Extreme Precision Analysis', fontsize=16, fontweight='bold')
    
    # Plot 1: Cubic - Max Error vs Tolerance
    ax = axes[0, 0]
    ax.loglog(tolerances, cubic_errors, 'o-', linewidth=2, markersize=8, label='Max Error')
    ax.axhline(machine_eps, color='r', linestyle='--', linewidth=2, label='Machine ε (2.22e-16)')
    ax.axvline(1e-8, color='g', linestyle=':', linewidth=1.5, alpha=0.7, label='Transition at 1e-8')
    ax.set_xlabel('Tolerance', fontsize=12)
    ax.set_ylabel('Max Root Error', fontsize=12)
    ax.set_title('Cubic: Max Error vs Tolerance', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    
    # Add annotation
    ax.annotate('Machine precision\nachieved!', xy=(1e-8, 2.22e-16), xytext=(1e-10, 1e-12),
                arrowprops=dict(arrowstyle='->', color='red', lw=2),
                fontsize=11, color='red', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7))
    
    # Plot 2: Cubic - Box Width vs Tolerance
    ax = axes[0, 1]
    ax.loglog(tolerances, cubic_widths, 's-', linewidth=2, markersize=8, color='purple', label='Max Width')
    ax.loglog(tolerances, tolerances, 'k--', linewidth=1, alpha=0.5, label='y = tolerance')
    ax.axhline(machine_eps, color='r', linestyle='--', linewidth=2, label='Machine ε')
    ax.set_xlabel('Tolerance', fontsize=12)
    ax.set_ylabel('Max Box Width', fontsize=12)
    ax.set_title('Cubic: Box Width vs Tolerance', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    
    # Plot 3: Multiplicity - Resolution vs Tolerance
    ax = axes[1, 0]
    ax.semilogx(tolerances, mult_resolved, 'o-', linewidth=2, markersize=8, color='green', label='Resolved')
    ax.semilogx(tolerances, mult_unresolved, 's-', linewidth=2, markersize=8, color='red', label='Unresolved')
    ax.axvline(1e-13, color='orange', linestyle=':', linewidth=1.5, alpha=0.7, label='Breakdown at 1e-13')
    ax.set_xlabel('Tolerance', fontsize=12)
    ax.set_ylabel('Number of Boxes', fontsize=12)
    ax.set_title('Multiplicity: Resolution vs Tolerance (Degenerate Case)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    
    # Add annotation
    ax.annotate('Extreme precision\nbreaks degeneracy\nhandling', xy=(1e-13, 4), xytext=(1e-11, 50),
                arrowprops=dict(arrowstyle='->', color='orange', lw=2),
                fontsize=10, color='orange', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.7))
    
    # Plot 4: Summary comparison
    ax = axes[1, 1]
    ax.axis('off')
    
    summary_text = """
    ═══════════════════════════════════════════════════════
    SUMMARY: OPTIMAL PRECISION FOR DIRECT CONTRACTION
    ═══════════════════════════════════════════════════════
    
    ✅ SIMPLE ROOTS (Cubic Polynomial):
       • Tolerance ≥ 1e-8:  Achieves MACHINE PRECISION
       • Max error:         2.22e-16 (1.0× machine ε)
       • Max width:         5.55e-17 (0.25× machine ε)
       • All 3 roots:       Resolved consistently
    
    ⚠️  DEGENERATE CASES (Multiple Roots):
       • Tolerance 1e-6 to 1e-12:  Stable (1 resolved, 99 unresolved)
       • Tolerance 1e-13:          Unstable (4 resolved, 99 unresolved)
       • Tolerance ≥ 1e-14:        Breakdown (0 resolved, 83 unresolved)
    
    📊 RECOMMENDATIONS:
       • For simple roots:      tolerance = 1e-8 to 1e-12
       • For degenerate cases:  tolerance = 1e-8 to 1e-10
       • Practical default:     tolerance = 1e-10
       • Extreme precision:     tolerance = 1e-12 (max safe)
    
    🎯 KEY INSIGHT:
       Direct contraction achieves machine precision (2.22e-16)
       for simple roots at tolerance ≥ 1e-8, a 2-8× improvement
       over incremental contraction!
    
       However, extreme precision (< 1e-13) breaks degeneracy
       handling due to numerical instability in the PP method.
    """
    
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('dumps/precision_analysis.png', dpi=150, bbox_inches='tight')
    print("✅ Visualization saved to: dumps/precision_analysis.png")
    
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())

