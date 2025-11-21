# Direct Contraction Implementation Plan

## Executive Summary

This document provides a detailed implementation plan for direct contraction from original [0,1] domain to eliminate error accumulation.

**Key Changes**:
1. Add `original_polys` and global interval tracking to `SubdivisionNode`
2. Modify contraction logic to update global interval and recompute from original
3. Modify subdivision logic to propagate original polynomials to children

**Estimated effort**: ~100 lines of code changes

---

## 1. Current Implementation Analysis

### 1.1 Data Structure: `SubdivisionNode`

**Location**: `src/solver.cpp:128-135`

```cpp
struct SubdivisionNode {
    std::vector<double> box_lower;  ///< Global box lower corner in [0,1]^n.
    std::vector<double> box_upper;  ///< Global box upper corner in [0,1]^n.
    unsigned int depth;             ///< Subdivision depth.

    // Polynomials restricted to this node's box and re-parameterised to [0,1]^n.
    std::vector<polynomial_solver::Polynomial> polys;
};
```

**Issue**: `box_lower` and `box_upper` are in global coordinates, but `polys` are always re-parameterized to local [0,1]^n. When we contract, we:
1. Get new bounds [a,b] in local [0,1] coordinates
2. Update global box: `new_low = old_low + a * old_width`
3. Restrict polynomials: `poly.restrictedToInterval(i, a, b)` (from current [0,1] to [a,b])

**Problem**: Each restriction accumulates errors.

### 1.2 Contraction Logic

**Location**: `src/solver.cpp:876-900`

```cpp
// Step 5: Contract the box (if strategy allows)
if (should_contract) {
    for (std::size_t i = 0; i < dim; ++i) {
        const double a = local_bound_lower[i];
        const double b = local_bound_upper[i];

        const double old_low = node.box_lower[i];
        const double old_high = node.box_upper[i];
        const double old_width = old_high - old_low;

        const double new_low = old_low + a * old_width;
        const double new_high = old_low + b * old_width;

        node.box_lower[i] = new_low;
        node.box_upper[i] = new_high;

        // Restrict polynomials to new interval
        for (std::size_t eq = 0; eq < node.polys.size(); ++eq) {
            node.polys[eq] = node.polys[eq].restrictedToInterval(i, a, b);
        }
    }
}
```

**Key observation**: 
- `local_bound_lower[i]` and `local_bound_upper[i]` are in local [0,1] coordinates
- We update global box correctly: `new_low = old_low + a * old_width`
- But we restrict from current [0,1] to [a,b], not from original [0,1] to global [new_low, new_high]

### 1.3 Subdivision Logic

**Location**: `src/solver.cpp:1059-1103`

```cpp
// Seed for children: start from the (possibly contracted) node.
std::vector<SubdivisionNode> children;
SubdivisionNode seed = node;
seed.depth = node.depth + 1u;
children.push_back(seed);

for (std::size_t axis = 0; axis < dim; ++axis) {
    if (!split_dim[axis]) {
        continue;
    }

    std::vector<SubdivisionNode> next_children;
    next_children.reserve(children.size() * 2u);

    for (const SubdivisionNode& child : children) {
        const double old_low = child.box_lower[axis];
        const double old_high = child.box_upper[axis];
        const double mid = 0.5 * (old_low + old_high);

        SubdivisionNode left_child = child;
        SubdivisionNode right_child = child;

        left_child.box_upper[axis] = mid;
        right_child.box_lower[axis] = mid;

        left_child.polys.clear();
        right_child.polys.clear();
        left_child.polys.reserve(child.polys.size());
        right_child.polys.reserve(child.polys.size());

        for (const Polynomial& poly : child.polys) {
            left_child.polys.push_back(poly.restrictedToInterval(axis, 0.0, 0.5));
            right_child.polys.push_back(poly.restrictedToInterval(axis, 0.5, 1.0));
        }

        next_children.push_back(std::move(left_child));
        next_children.push_back(std::move(right_child));
    }

    children.swap(next_children);
}
```

**Key observation**:
- Children inherit parent's (contracted) polynomials
- Subdivision restricts from current [0,1] to [0,0.5] and [0.5,1]
- Global box is updated correctly (midpoint in global coordinates)

### 1.4 Root Node Initialization

**Location**: `src/solver.cpp:609-614`

```cpp
// Root node: full box [0,1]^dim with the original equations.
SubdivisionNode root;
root.box_lower.assign(dim, 0.0);
root.box_upper.assign(dim, 1.0);
root.depth = 0u;
root.polys = system.equations();
```

**Key observation**: Root node has global box = [0,1]^n and polys = original equations.

---

## 2. Proposed Changes

### 2.1 Modify `SubdivisionNode` Structure

**Location**: `src/solver.cpp:128-135`

**Change**:
```cpp
struct SubdivisionNode {
    std::vector<double> box_lower;  ///< Global box lower corner in [0,1]^n.
    std::vector<double> box_upper;  ///< Global box upper corner in [0,1]^n.
    unsigned int depth;             ///< Subdivision depth.

    // Polynomials restricted to this node's box and re-parameterised to [0,1]^n.
    std::vector<polynomial_solver::Polynomial> polys;
    
    // NEW: Original polynomials (on [0,1]^n, never modified)
    std::vector<polynomial_solver::Polynomial> original_polys;
};
```

**Rationale**: Store original polynomials to enable direct contraction.

**Memory impact**: 2× polynomial storage per node. For degree 10, 2D: ~200 doubles per node. Typical solver has 100-1000 nodes, so ~20-200 KB total (negligible).

### 2.2 Modify Root Node Initialization

**Location**: `src/solver.cpp:609-614`

**Change**:
```cpp
// Root node: full box [0,1]^dim with the original equations.
SubdivisionNode root;
root.box_lower.assign(dim, 0.0);
root.box_upper.assign(dim, 1.0);
root.depth = 0u;
root.polys = system.equations();
root.original_polys = system.equations();  // NEW: Store original
```

**Rationale**: Initialize original_polys at root.

### 2.3 Modify Contraction Logic

**Location**: `src/solver.cpp:876-900`

**Change**:
```cpp
// Step 5: Contract the box (if strategy allows)
if (should_contract) {
    // First, update global box coordinates
    for (std::size_t i = 0; i < dim; ++i) {
        // For Simultaneous strategy, only contract directions that don't need subdivision
        if (config.strategy == SubdivisionStrategy::Simultaneous && needs_subdivision[i]) {
            continue;  // Skip contraction for this direction
        }

        const double a = local_bound_lower[i];
        const double b = local_bound_upper[i];

        const double old_low = node.box_lower[i];
        const double old_high = node.box_upper[i];
        const double old_width = old_high - old_low;

        const double new_low = old_low + a * old_width;
        const double new_high = old_low + b * old_width;

        node.box_lower[i] = new_low;
        node.box_upper[i] = new_high;
    }
    
    // NEW: Recompute polynomials from original using global box
    for (std::size_t eq = 0; eq < node.polys.size(); ++eq) {
        node.polys[eq] = node.original_polys[eq];
        for (std::size_t i = 0; i < dim; ++i) {
            // Skip if this dimension wasn't contracted (Simultaneous strategy)
            if (config.strategy == SubdivisionStrategy::Simultaneous && needs_subdivision[i]) {
                continue;
            }
            
            // Restrict from original [0,1] to global [box_lower[i], box_upper[i]]
            node.polys[eq] = node.polys[eq].restrictedToInterval(
                i, node.box_lower[i], node.box_upper[i]);
        }
    }
}
```

**Rationale**: 
- Update global box first (same as before)
- Then recompute polynomials from original using global box coordinates
- This eliminates error accumulation!

**Cost**: Same as before (2 subdivisions per dimension per contraction).

### 2.4 Modify Subdivision Logic

**Location**: `src/solver.cpp:1059-1103`

**Change**:
```cpp
// Seed for children: start from the (possibly contracted) node.
std::vector<SubdivisionNode> children;
SubdivisionNode seed = node;
seed.depth = node.depth + 1u;
children.push_back(seed);

for (std::size_t axis = 0; axis < dim; ++axis) {
    if (!split_dim[axis]) {
        continue;
    }

    std::vector<SubdivisionNode> next_children;
    next_children.reserve(children.size() * 2u);

    for (const SubdivisionNode& child : children) {
        const double old_low = child.box_lower[axis];
        const double old_high = child.box_upper[axis];
        const double mid = 0.5 * (old_low + old_high);

        SubdivisionNode left_child = child;
        SubdivisionNode right_child = child;

        left_child.box_upper[axis] = mid;
        right_child.box_lower[axis] = mid;

        // NEW: Recompute polynomials from original using global box
        left_child.polys.clear();
        right_child.polys.clear();
        left_child.polys.reserve(child.original_polys.size());
        right_child.polys.reserve(child.original_polys.size());

        for (const Polynomial& orig_poly : child.original_polys) {
            Polynomial left_poly = orig_poly;
            Polynomial right_poly = orig_poly;

            // Restrict from original [0,1] to global box for each child
            for (std::size_t i = 0; i < dim; ++i) {
                left_poly = left_poly.restrictedToInterval(
                    i, left_child.box_lower[i], left_child.box_upper[i]);
                right_poly = right_poly.restrictedToInterval(
                    i, right_child.box_lower[i], right_child.box_upper[i]);
            }

            left_child.polys.push_back(std::move(left_poly));
            right_child.polys.push_back(std::move(right_poly));
        }

        // NEW: Children inherit original polynomials
        left_child.original_polys = child.original_polys;
        right_child.original_polys = child.original_polys;

        next_children.push_back(std::move(left_child));
        next_children.push_back(std::move(right_child));
    }

    children.swap(next_children);
}
```

**Rationale**:
- Children inherit original polynomials from parent
- Recompute current polynomials from original using global box coordinates
- This ensures no error accumulation through subdivision either

**Cost**: More expensive than before! We now restrict from [0,1] to global box (2 subdivisions × dim dimensions), instead of just [0,1] to [0,0.5] or [0.5,1] (2 subdivisions × 1 dimension).

**Trade-off**: Subdivision becomes dim× more expensive, but contraction stays the same. Since we typically do many contractions per subdivision, this is acceptable.

---

## 3. Alternative Implementation (More Efficient)

The above implementation makes subdivision more expensive. Here's a more efficient alternative:

### 3.1 Keep Incremental Subdivision, Direct Contraction Only

**Idea**: Only use direct contraction for the contraction step, but keep incremental subdivision.

**Rationale**:
- Contraction happens many times per node (up to 100 iterations)
- Subdivision happens once per node
- Error accumulation is worse for contraction (many iterations)

**Modified Subdivision Logic**:
```cpp
// Keep the original incremental subdivision (more efficient)
for (const Polynomial& poly : child.polys) {
    left_child.polys.push_back(poly.restrictedToInterval(axis, 0.0, 0.5));
    right_child.polys.push_back(poly.restrictedToInterval(axis, 0.5, 1.0));
}

// But propagate original polynomials
left_child.original_polys = child.original_polys;
right_child.original_polys = child.original_polys;
```

**Trade-off**:
- Contraction: Direct from original (no error accumulation)
- Subdivision: Incremental (some error accumulation, but only once per node)
- Overall: Much better than current, but not perfect

**Recommendation**: Start with this approach (simpler, faster), then measure if full direct approach is needed.

---

## 4. Implementation Steps

### Step 1: Add `original_polys` to `SubdivisionNode`

**File**: `src/solver.cpp`
**Line**: 128-135

```cpp
struct SubdivisionNode {
    std::vector<double> box_lower;
    std::vector<double> box_upper;
    unsigned int depth;
    std::vector<polynomial_solver::Polynomial> polys;
    std::vector<polynomial_solver::Polynomial> original_polys;  // ADD THIS
};
```

### Step 2: Initialize `original_polys` at root

**File**: `src/solver.cpp`
**Line**: 609-614

```cpp
SubdivisionNode root;
root.box_lower.assign(dim, 0.0);
root.box_upper.assign(dim, 1.0);
root.depth = 0u;
root.polys = system.equations();
root.original_polys = system.equations();  // ADD THIS
```

### Step 3: Modify contraction logic (direct from original)

**File**: `src/solver.cpp`
**Line**: 876-900

Replace the contraction loop with:
```cpp
if (should_contract) {
    // Update global box
    for (std::size_t i = 0; i < dim; ++i) {
        if (config.strategy == SubdivisionStrategy::Simultaneous && needs_subdivision[i]) {
            continue;
        }
        const double a = local_bound_lower[i];
        const double b = local_bound_upper[i];
        const double old_low = node.box_lower[i];
        const double old_high = node.box_upper[i];
        const double old_width = old_high - old_low;
        node.box_lower[i] = old_low + a * old_width;
        node.box_upper[i] = old_low + b * old_width;
    }

    // Recompute from original
    for (std::size_t eq = 0; eq < node.polys.size(); ++eq) {
        node.polys[eq] = node.original_polys[eq];
        for (std::size_t i = 0; i < dim; ++i) {
            node.polys[eq] = node.polys[eq].restrictedToInterval(
                i, node.box_lower[i], node.box_upper[i]);
        }
    }
}
```

### Step 4: Propagate `original_polys` in subdivision

**File**: `src/solver.cpp`
**Line**: 1089-1092

Add after creating children:
```cpp
for (const Polynomial& poly : child.polys) {
    left_child.polys.push_back(poly.restrictedToInterval(axis, 0.0, 0.5));
    right_child.polys.push_back(poly.restrictedToInterval(axis, 0.5, 1.0));
}

// ADD THESE TWO LINES:
left_child.original_polys = child.original_polys;
right_child.original_polys = child.original_polys;
```

### Step 5: Test

Run all examples and verify:
1. Results are the same or better (more roots resolved)
2. Errors are smaller (2-8× improvement expected)
3. No performance regression

---

## 5. Testing Strategy

### 5.1 Correctness Tests

Run all existing examples:
```bash
./build/examples/cubic_1d_roots
./build/examples/multiplicity_1d_roots
./build/examples/wilkinson_1d_roots
./build/examples/circle_ellipse_intersection
```

**Expected**: Same or more roots resolved.

### 5.2 Accuracy Tests

Use `tools/test_extreme_precision_examples.py` to test at various tolerances:
```bash
python tools/test_extreme_precision_examples.py
```

**Expected**:
- At tolerance 10⁻⁸: Same results
- At tolerance 10⁻¹²: More roots resolved
- At tolerance 10⁻¹⁵: Fewer unresolved boxes (better degeneracy handling)

### 5.3 Performance Tests

Measure runtime for each example before and after.

**Expected**: Similar or slightly slower (due to recomputing from original), but acceptable.

---

## 6. Summary

### Changes Required

| File | Lines | Change |
|------|-------|--------|
| `src/solver.cpp` | 128-135 | Add `original_polys` field |
| `src/solver.cpp` | 614 | Initialize `original_polys` at root |
| `src/solver.cpp` | 876-900 | Modify contraction to use direct approach |
| `src/solver.cpp` | 1092 | Propagate `original_polys` in subdivision |

**Total**: ~30 lines of code changes

### Expected Benefits

1. ✅ **2-8× error reduction** in contraction step
2. ✅ **Same computational cost** for contraction
3. ✅ **Better degeneracy handling** at extreme precision
4. ✅ **More roots resolved** for ill-conditioned problems

### Recommendation

**Implement the efficient version** (Section 3.1):
- Direct contraction from original
- Incremental subdivision (keep current approach)
- Minimal code changes (~30 lines)
- Maximum benefit with minimal cost

This provides most of the benefit (contraction is where most iterations happen) with minimal implementation complexity.


