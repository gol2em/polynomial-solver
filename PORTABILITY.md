# Portability Guide

This document explains how to move the Polynomial Solver project to another environment.

## Quick Summary

The project is **fully portable** and **self-contained**. No external dependencies beyond standard C++11 libraries.

## What You Need

### Minimum Requirements
- C++11 compiler (GCC 4.8+, Clang 3.4+)
- CMake 3.15+
- Make
- Git

### Optional
- Python 3.11+ (for visualization only)
- NumPy, Matplotlib (for visualization only)

## Two Ways to Move

### Method 1: Git Clone (Recommended)

If the target environment has internet access:

```bash
git clone https://github.com/gol2em/polynomial-solver.git
cd polynomial-solver
./check_prerequisites.sh
./build.sh --test
```

### Method 2: Manual Transfer

If the target environment is offline:

**On source machine:**
```bash
cd polynomial-solver
tar -czf polynomial-solver.tar.gz .
# Creates: polynomial-solver.tar.gz
```

**Transfer the file** (USB, SCP, etc.)

**On target machine:**
```bash
mkdir polynomial-solver
tar -xzf polynomial-solver.tar.gz -C polynomial-solver
cd polynomial-solver
./check_prerequisites.sh
./build.sh --test
```

**Alternative: Git Bundle (With History)**

To preserve Git history without internet:

**On source machine:**
```bash
cd polynomial-solver
git bundle create polynomial-solver.bundle --all
```

**On target machine:**
```bash
git clone polynomial-solver.bundle polynomial-solver
cd polynomial-solver
./check_prerequisites.sh
./build.sh --test
```

## Verification Checklist

After moving to new environment:

- [ ] Run `./check_prerequisites.sh` - All checks pass
- [ ] Run `./build.sh --test` - Build succeeds
- [ ] Check test results - 100% tests passed (12/12)
- [ ] Run standard verification workflow (see below)

## Standard Verification Workflow

After setup, run the circle-ellipse intersection example:

```bash
# 1. Build and test
./build.sh --test

# 2. Run circle-ellipse intersection
./build/bin/test_strategies

# 3. Set up Python environment (one-time)
uv venv .venv
source .venv/bin/activate
uv pip install numpy matplotlib

# 4. Visualize results
python examples/visualize_circle_ellipse.py
```

**Expected Results:**
- All 12 tests pass
- 3 geometry dumps generated in `dumps/`
- Visualizations generated in `visualizations/`
- All strategies converge to root: (0.894427, 0.447214)

## Package Contents

### Essential Files
```
polynomial-solver/
├── include/              # Header files (4 files)
├── src/                  # Source files (5 files)
├── tests/                # Test suite (12 files)
├── examples/             # Example programs and scripts
├── tools/                # Visualization tools
├── CMakeLists.txt        # Build configuration
├── build.sh              # Build script
└── check_prerequisites.sh # Prerequisite checker
```

### Documentation
```
├── README.md             # Project overview
├── QUICKSTART.md         # Quick start guide
├── SETUP.md              # Complete setup guide
├── PORTABILITY.md        # This file
├── PROJECT_SUMMARY.md    # Project summary
├── examples/README.md    # Example documentation
├── tools/README.md       # Visualization API docs
└── docs/                 # Algorithm documentation
```

### Generated Files (Not Included)
```
├── build/                # Build artifacts (regenerated)
├── .venv/                # Python environment (regenerated)
├── dumps/                # Geometry dumps (regenerated)
└── visualizations/       # PNG visualizations (regenerated)
```

## Platform Support

### Tested Platforms
- ✅ Ubuntu 20.04+ (WSL2, Native)
- ✅ Debian 10+
- ✅ Fedora 30+
- ✅ Arch Linux

### Should Work On
- macOS 10.14+ (with Xcode Command Line Tools)
- CentOS 7+ / RHEL 7+
- Any Linux with GCC 4.8+ or Clang 3.4+

### Not Supported
- Windows native (use WSL2)
- Very old compilers without C++11 support

## Common Scenarios

### Scenario 1: Moving to HPC Cluster

```bash
# On local machine
tar -czf polynomial-solver.tar.gz polynomial-solver/

# Transfer to cluster
scp polynomial-solver.tar.gz user@cluster:~

# On cluster
ssh user@cluster
tar -xzf polynomial-solver.tar.gz
cd polynomial-solver
module load gcc/9.3.0 cmake/3.18.0  # Load required modules
./check_prerequisites.sh
./build.sh --test -j 16  # Use more cores
```

### Scenario 2: Moving to Air-Gapped System

```bash
# On internet-connected machine
tar -czf polynomial-solver.tar.gz polynomial-solver/

# Transfer via USB or approved method
# On air-gapped machine
tar -xzf polynomial-solver.tar.gz
cd polynomial-solver
./check_prerequisites.sh
# If prerequisites missing, install from offline repos
./build.sh --test
```

### Scenario 3: Moving to Docker Container

```bash
# Create Dockerfile
cat > Dockerfile << 'EOF'
FROM ubuntu:22.04
RUN apt-get update && apt-get install -y build-essential cmake git
COPY polynomial-solver /app/polynomial-solver
WORKDIR /app/polynomial-solver
RUN ./build.sh --test
CMD ["./build/bin/test_method_comparison"]
EOF

# Build and run
docker build -t polynomial-solver .
docker run polynomial-solver
```

## Troubleshooting

### Issue: Prerequisites missing on target

**Solution**: Use `check_prerequisites.sh` output to identify missing tools, then:
- Ubuntu/Debian: `sudo apt-get install build-essential cmake git`
- Fedora/RHEL: `sudo dnf install gcc-c++ cmake git`
- Arch: `sudo pacman -S base-devel cmake git`

### Issue: CMake too old

**Solution**: Download newer CMake from https://cmake.org/download/ or use pip:
```bash
pip3 install cmake --upgrade
```

### Issue: Compiler too old

**Solution**: Install newer GCC:
```bash
# Ubuntu/Debian
sudo apt-get install gcc-9 g++-9
export CC=gcc-9
export CXX=g++-9
```

## What Gets Regenerated

These directories/files are **not** in the package (will be regenerated):

- `build/` - Build artifacts (created by `./build.sh`)
- `.venv/` - Python virtual environment (created by `uv sync`)
- `*.png`, `*.csv` - Generated visualization files
- `*.o`, `*.a` - Compiled object files and libraries

## Summary

✅ **Fully portable** - No external dependencies
✅ **Self-contained** - All source included
✅ **Easy setup** - 2 commands to build
✅ **Well documented** - Multiple guides included
✅ **Tested** - Verified in isolated environment
✅ **Standard verification** - Circle-ellipse intersection workflow

## Next Steps

1. Choose your transfer method (Git or Manual)
2. Transfer to target environment
3. Run `./check_prerequisites.sh`
4. Run `./build.sh --test`
5. Run standard verification workflow
6. Start developing!

For more details, see:
- [QUICKSTART.md](QUICKSTART.md) - Quick start guide
- [SETUP.md](SETUP.md) - Complete setup guide
- [README.md](README.md) - Project overview

