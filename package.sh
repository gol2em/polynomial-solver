#!/bin/bash

# Polynomial Solver - Packaging Script
# Creates a clean distribution package for moving to another environment

set -e

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Default options
INCLUDE_DOCS=true
INCLUDE_PYTHON=true
INCLUDE_GIT=false
OUTPUT_NAME="polynomial-solver-package"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --no-docs)
            INCLUDE_DOCS=false
            shift
            ;;
        --no-python)
            INCLUDE_PYTHON=false
            shift
            ;;
        --with-git)
            INCLUDE_GIT=true
            shift
            ;;
        --output|-o)
            OUTPUT_NAME="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --no-docs        Exclude documentation files"
            echo "  --no-python      Exclude Python visualization tools"
            echo "  --with-git       Include Git history"
            echo "  --output, -o     Output filename (without extension)"
            echo "  --help, -h       Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0                           # Full package"
            echo "  $0 --no-docs --no-python     # Minimal package"
            echo "  $0 --with-git -o my-package  # Include Git history"
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Polynomial Solver - Packaging Script${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Check if we're in the right directory
if [ ! -f "CMakeLists.txt" ]; then
    echo -e "${RED}Error: CMakeLists.txt not found!${NC}"
    echo -e "${RED}Please run this script from the project root directory.${NC}"
    exit 1
fi

# Create temporary directory
TEMP_DIR=$(mktemp -d)
PACKAGE_DIR="$TEMP_DIR/polynomial-solver"

echo -e "${BLUE}Creating package structure...${NC}"
mkdir -p "$PACKAGE_DIR"

# Copy essential files
echo -e "  Copying source files..."
cp -r include "$PACKAGE_DIR/"
cp -r src "$PACKAGE_DIR/"
cp -r tests "$PACKAGE_DIR/"
cp CMakeLists.txt "$PACKAGE_DIR/"

# Copy scripts
echo -e "  Copying scripts..."
cp build.sh "$PACKAGE_DIR/"
cp check_prerequisites.sh "$PACKAGE_DIR/"
chmod +x "$PACKAGE_DIR/build.sh"
chmod +x "$PACKAGE_DIR/check_prerequisites.sh"

# Copy documentation
if [ "$INCLUDE_DOCS" = true ]; then
    echo -e "  Copying documentation..."
    cp README.md "$PACKAGE_DIR/" 2>/dev/null || true
    cp QUICKSTART.md "$PACKAGE_DIR/" 2>/dev/null || true
    cp SETUP.md "$PACKAGE_DIR/" 2>/dev/null || true
    [ -d "docs" ] && cp -r docs "$PACKAGE_DIR/" || true
fi

# Copy Python tools
if [ "$INCLUDE_PYTHON" = true ]; then
    echo -e "  Copying Python tools..."
    [ -d "python" ] && cp -r python "$PACKAGE_DIR/" || true
fi

# Copy Git history
if [ "$INCLUDE_GIT" = true ]; then
    echo -e "  Copying Git history..."
    cp -r .git "$PACKAGE_DIR/" 2>/dev/null || echo -e "  ${YELLOW}Warning: No Git history found${NC}"
fi

# Create archive
echo -e "${BLUE}Creating archive...${NC}"
cd "$TEMP_DIR"
tar -czf "${OUTPUT_NAME}.tar.gz" polynomial-solver/
cd - > /dev/null

# Move archive to current directory
mv "$TEMP_DIR/${OUTPUT_NAME}.tar.gz" .

# Calculate size
SIZE=$(du -h "${OUTPUT_NAME}.tar.gz" | cut -f1)

# Cleanup
rm -rf "$TEMP_DIR"

# Summary
echo -e "${GREEN}✓ Package created successfully!${NC}"
echo ""
echo -e "${BLUE}Package Information:${NC}"
echo -e "  Filename: ${GREEN}${OUTPUT_NAME}.tar.gz${NC}"
echo -e "  Size: ${GREEN}${SIZE}${NC}"
echo -e "  Location: ${GREEN}$(pwd)/${OUTPUT_NAME}.tar.gz${NC}"
echo ""
echo -e "${BLUE}Package Contents:${NC}"
echo -e "  ✓ Source code (include/, src/)"
echo -e "  ✓ Test suite (tests/)"
echo -e "  ✓ Build scripts (build.sh, check_prerequisites.sh)"
[ "$INCLUDE_DOCS" = true ] && echo -e "  ✓ Documentation (README.md, QUICKSTART.md, SETUP.md, docs/)"
[ "$INCLUDE_PYTHON" = true ] && echo -e "  ✓ Python tools (python/)"
[ "$INCLUDE_GIT" = true ] && echo -e "  ✓ Git history (.git/)"
echo ""
echo -e "${BLUE}To extract and use:${NC}"
echo -e "  1. Copy ${OUTPUT_NAME}.tar.gz to target environment"
echo -e "  2. Extract: ${GREEN}tar -xzf ${OUTPUT_NAME}.tar.gz${NC}"
echo -e "  3. Navigate: ${GREEN}cd polynomial-solver${NC}"
echo -e "  4. Check prerequisites: ${GREEN}./check_prerequisites.sh${NC}"
echo -e "  5. Build: ${GREEN}./build.sh --test${NC}"
echo ""

