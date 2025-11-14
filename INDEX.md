# Documentation Index

Complete guide to all documentation in the Polynomial Solver project.

## Quick Navigation

| Document | Purpose | When to Use |
|----------|---------|-------------|
| [README.md](README.md) | Project overview | First time visitors |
| [QUICKSTART.md](QUICKSTART.md) | 5-minute setup | Want to get started fast |
| [SETUP.md](SETUP.md) | Complete setup guide | Detailed installation |
| [PORTABILITY.md](PORTABILITY.md) | Moving to new environment | Transferring project |
| [INDEX.md](INDEX.md) | This file | Finding documentation |

## For New Users

### I want to get started quickly
→ Read [QUICKSTART.md](QUICKSTART.md)

### I want to understand the project
→ Read [README.md](README.md)

### I need detailed installation instructions
→ Read [SETUP.md](SETUP.md)

### I want to move to another environment
→ Read [PORTABILITY.md](PORTABILITY.md)

## Documentation Structure

### Main Documentation

#### [README.md](README.md)
**Purpose**: Project overview and introduction

**Contents**:
- Features and capabilities
- Project structure
- Quick start (brief)
- Algorithm overview
- Test suite summary
- Moving to another environment

**Read this if**: You're new to the project or want a high-level overview.

#### [QUICKSTART.md](QUICKSTART.md)
**Purpose**: Get up and running in 5 minutes

**Contents**:
- Prerequisites
- Installation steps (3 steps)
- Running examples
- Basic usage code
- Configuration options
- Troubleshooting

**Read this if**: You want to start using the project immediately.

#### [SETUP.md](SETUP.md)
**Purpose**: Comprehensive setup guide

**Contents**:
- System requirements (detailed table)
- Installation for all platforms
- Verification steps
- Troubleshooting (detailed)
- Development setup
- Moving to another environment

**Read this if**: You need detailed installation instructions or encounter issues.

#### [PORTABILITY.md](PORTABILITY.md)
**Purpose**: Moving project to another environment

**Contents**:
- Three transfer methods (Git, Package, Bundle)
- Package contents and sizes
- Platform support
- Common scenarios (HPC, air-gapped, Docker)
- Troubleshooting

**Read this if**: You need to move the project to a different machine or environment.

### Algorithm Documentation

Located in `docs/` directory:

#### [docs/GEOMETRY_ALGORITHMS.md](docs/GEOMETRY_ALGORITHMS.md)
**Purpose**: Geometric algorithm details

**Contents**:
- Convex hull algorithms (2D, 3D)
- Hyperplane intersection
- Algorithm complexity
- Implementation notes

**Read this if**: You want to understand the geometric algorithms used.

#### [docs/GEOMETRY_ROBUSTNESS.md](docs/GEOMETRY_ROBUSTNESS.md)
**Purpose**: Robustness improvements

**Contents**:
- Intrinsic vs ambient dimension
- Degenerate case handling
- Numerical stability
- Design decisions

**Read this if**: You want to understand how robustness is achieved.

#### [docs/DEGENERATE_BOXES.md](docs/DEGENERATE_BOXES.md)
**Purpose**: Degeneracy detection and handling

**Contents**:
- Degeneracy detection algorithm
- Threshold calculation
- Unresolved box handling
- Examples

**Read this if**: You want to understand degeneracy detection.

### Scripts and Tools

#### `check_prerequisites.sh`
**Purpose**: Verify all required tools are installed

**Usage**:
```bash
./check_prerequisites.sh
```

**Output**: Color-coded report of installed/missing tools with installation hints.

#### `build.sh`
**Purpose**: Automated build script

**Usage**:
```bash
./build.sh [options]
```

**Options**:
- `--debug` - Debug build
- `--test` - Run tests after build
- `--clean` - Clean build
- `-j N` - Use N parallel jobs
- `--verbose` - Verbose output
- `--help` - Show help

#### `package.sh`
**Purpose**: Create distribution package

**Usage**:
```bash
./package.sh [options]
```

**Options**:
- `--no-docs` - Exclude documentation
- `--no-python` - Exclude Python tools
- `--with-git` - Include Git history
- `-o NAME` - Output filename
- `--help` - Show help

## Common Tasks

### Task: First Time Setup

1. Read [QUICKSTART.md](QUICKSTART.md)
2. Run `./check_prerequisites.sh`
3. Run `./build.sh --test`
4. Run examples

### Task: Understanding Algorithms

1. Read [README.md](README.md) - Algorithm overview section
2. Read [docs/GEOMETRY_ALGORITHMS.md](docs/GEOMETRY_ALGORITHMS.md)
3. Read [docs/GEOMETRY_ROBUSTNESS.md](docs/GEOMETRY_ROBUSTNESS.md)
4. Look at test files in `tests/`

### Task: Moving to New Environment

1. Read [PORTABILITY.md](PORTABILITY.md)
2. Choose transfer method
3. Run `./package.sh` (if using package method)
4. Transfer and extract
5. Run `./check_prerequisites.sh`
6. Run `./build.sh --test`

### Task: Troubleshooting Build Issues

1. Check [QUICKSTART.md](QUICKSTART.md) - Troubleshooting section
2. Check [SETUP.md](SETUP.md) - Troubleshooting section
3. Run `./check_prerequisites.sh` to verify environment
4. Try clean build: `./build.sh --clean --verbose`

### Task: Development

1. Read [SETUP.md](SETUP.md) - Development setup section
2. Build in debug mode: `./build.sh --debug`
3. Run specific tests: `./build/bin/test_*`
4. Read algorithm docs in `docs/`

## File Organization

```
polynomial-solver/
├── README.md              # Project overview
├── QUICKSTART.md          # Quick start guide
├── SETUP.md               # Complete setup guide
├── PORTABILITY.md         # Portability guide
├── INDEX.md               # This file
├── docs/                  # Algorithm documentation
│   ├── GEOMETRY_ALGORITHMS.md
│   ├── GEOMETRY_ROBUSTNESS.md
│   └── DEGENERATE_BOXES.md
├── check_prerequisites.sh # Prerequisite checker
├── build.sh               # Build script
├── package.sh             # Packaging script
├── include/               # Header files
├── src/                   # Source files
├── tests/                 # Test suite
└── python/                # Visualization tools
```

## Support

- **Email**: 664862601@qq.com
- **GitHub**: https://github.com/gol2em/polynomial-solver
- **Issues**: https://github.com/gol2em/polynomial-solver/issues

