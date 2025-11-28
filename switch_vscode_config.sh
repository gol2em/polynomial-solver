#!/bin/bash
# Switch VS Code IntelliSense configuration for different high-precision backends

set -e

# Colors
BLUE='\033[0;34m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Function to show usage
show_usage() {
    echo -e "${BLUE}Usage: $0 <configuration>${NC}"
    echo ""
    echo "Available configurations:"
    echo "  1. mpfr          - MPFR backend (default, fastest)"
    echo "  2. cpp_dec_float - cpp_dec_float backend (Boost only)"
    echo "  3. quadmath      - Quadmath backend (standalone)"
    echo "  4. no_templates  - MPFR backend without templates"
    echo "  5. double        - Double precision only (no high-precision)"
    echo ""
    echo "Example:"
    echo "  $0 mpfr"
    echo "  $0 quadmath"
    echo ""
    echo "After switching, reload VS Code window:"
    echo "  Ctrl+Shift+P → 'Developer: Reload Window'"
}

# Check arguments
if [ $# -ne 1 ]; then
    show_usage
    exit 1
fi

CONFIG=$1
VSCODE_DIR=".vscode"
CONFIG_FILE="$VSCODE_DIR/c_cpp_properties.json"

# Function to write configuration
write_config() {
    local defines="$1"
    local include_paths="$2"
    
    cat > "$CONFIG_FILE" << EOF
{
    "configurations": [
        {
            "name": "Linux",
            "includePath": [
$include_paths
            ],
            "defines": [
$defines
            ],
            "compilerPath": "/usr/local/bin/c++",
            "cStandard": "c17",
            "cppStandard": "c++11",
            "intelliSenseMode": "linux-gcc-x64",
            "compileCommands": "\${workspaceFolder}/compile_commands.json",
            "configurationProvider": "ms-vscode.cmake-tools"
        }
    ],
    "version": 4
}
EOF
}

# Common include paths
COMMON_INCLUDES='                "${workspaceFolder}/include",
                "${workspaceFolder}/build/include"'

BOOST_INCLUDES='                "${workspaceFolder}/include",
                "${workspaceFolder}/build/include",
                "/usr/include",
                "/usr/include/x86_64-linux-gnu"'

# Switch configuration
case "$CONFIG" in
    mpfr)
        echo -e "${YELLOW}Switching to MPFR backend configuration...${NC}"
        DEFINES='                "ENABLE_HIGH_PRECISION",
                "USE_MPFR_BACKEND",
                "ENABLE_QUADMATH"'
        write_config "$DEFINES" "$BOOST_INCLUDES"
        echo -e "${GREEN}✓ Switched to MPFR backend${NC}"
        ;;
        
    cpp_dec_float)
        echo -e "${YELLOW}Switching to cpp_dec_float backend configuration...${NC}"
        DEFINES='                "ENABLE_HIGH_PRECISION",
                "USE_CPP_DEC_FLOAT_BACKEND",
                "ENABLE_QUADMATH"'
        write_config "$DEFINES" "$BOOST_INCLUDES"
        echo -e "${GREEN}✓ Switched to cpp_dec_float backend${NC}"
        ;;
        
    quadmath)
        echo -e "${YELLOW}Switching to Quadmath backend configuration...${NC}"
        DEFINES='                "ENABLE_HIGH_PRECISION",
                "USE_QUADMATH_BACKEND",
                "ENABLE_QUADMATH"'
        write_config "$DEFINES" "$COMMON_INCLUDES"
        echo -e "${GREEN}✓ Switched to Quadmath backend${NC}"
        ;;
        
    no_templates)
        echo -e "${YELLOW}Switching to No Templates configuration...${NC}"
        DEFINES='                "ENABLE_HIGH_PRECISION",
                "USE_MPFR_BACKEND",
                "DISABLE_TEMPLATES"'
        write_config "$DEFINES" "$BOOST_INCLUDES"
        echo -e "${GREEN}✓ Switched to No Templates (MPFR backend)${NC}"
        ;;
        
    double)
        echo -e "${YELLOW}Switching to Double Precision Only configuration...${NC}"
        DEFINES=''
        write_config "$DEFINES" "$COMMON_INCLUDES"
        echo -e "${GREEN}✓ Switched to Double Precision Only${NC}"
        ;;
        
    *)
        echo -e "${YELLOW}Unknown configuration: $CONFIG${NC}"
        show_usage
        exit 1
        ;;
esac

echo ""
echo -e "${BLUE}Next steps:${NC}"
echo "1. Reload VS Code window: Ctrl+Shift+P → 'Developer: Reload Window'"
echo "2. Wait for IntelliSense to reindex (watch bottom-right status bar)"
echo "3. Code blocks should now be highlighted correctly!"

