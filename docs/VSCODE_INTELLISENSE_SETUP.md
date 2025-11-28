# VS Code IntelliSense Setup for High-Precision Code

This document explains how to configure VS Code to correctly highlight code blocks based on preprocessor macros.

## Problem

When using conditional compilation with `#ifdef`, `#ifndef`, etc., VS Code needs to know which macros are defined to correctly highlight the active code blocks and gray out the inactive ones.

For example:
```cpp
#ifdef ENABLE_HIGH_PRECISION
    // This code should be highlighted when ENABLE_HIGH_PRECISION is defined
    using mpreal = boost::multiprecision::mpfr_float;
#else
    // This code should be grayed out when ENABLE_HIGH_PRECISION is defined
    using mpreal = double;
#endif
```

## Solution

There are **three methods** to configure IntelliSense. We've implemented all three for maximum flexibility.

---

## Method 1: Use `compile_commands.json` (Automatic) ✅

This is the **recommended** method because it automatically uses the actual build configuration.

### How it works:
1. CMake generates `compile_commands.json` with all compiler flags
2. VS Code reads this file to understand which macros are defined
3. IntelliSense automatically highlights the correct code blocks

### Setup:

**Step 1:** Generate `compile_commands.json` during CMake configuration:
```bash
cd build
cmake .. -DENABLE_HIGH_PRECISION=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
```

**Step 2:** Create a symlink in the project root (already done):
```bash
ln -sf build/compile_commands.json compile_commands.json
```

**Step 3:** VS Code configuration (already in `.vscode/c_cpp_properties.json`):
```json
{
    "compileCommands": "${workspaceFolder}/compile_commands.json"
}
```

**Step 4:** Reload VS Code:
- Press `Ctrl+Shift+P`
- Type "Developer: Reload Window"
- Wait for IntelliSense to reindex (watch bottom-right status bar)

### Advantages:
✅ Automatic - uses actual build configuration  
✅ Always in sync with CMake settings  
✅ No manual updates needed  

### Disadvantages:
❌ Requires rebuilding CMake when switching configurations  

---

## Method 2: Manual Defines in `c_cpp_properties.json` ✅

This method manually specifies which macros are defined.

### Setup:

Edit `.vscode/c_cpp_properties.json`:
```json
{
    "configurations": [{
        "name": "Linux",
        "defines": [
            "ENABLE_HIGH_PRECISION",
            "USE_MPFR_BACKEND",
            "ENABLE_QUADMATH"
        ],
        "includePath": [
            "${workspaceFolder}/include",
            "${workspaceFolder}/build/include",
            "/usr/include",
            "/usr/include/x86_64-linux-gnu"
        ]
    }]
}
```

### Switching Configurations:

We provide a helper script to easily switch between configurations:

```bash
# Switch to MPFR backend
./switch_vscode_config.sh mpfr

# Switch to cpp_dec_float backend
./switch_vscode_config.sh cpp_dec_float

# Switch to Quadmath backend
./switch_vscode_config.sh quadmath

# Switch to no templates
./switch_vscode_config.sh no_templates

# Switch to double precision only
./switch_vscode_config.sh double
```

After switching, reload VS Code window (Ctrl+Shift+P → "Developer: Reload Window").

### Advantages:
✅ Fast - no need to rebuild CMake  
✅ Easy to switch between configurations  
✅ Works even if build directory is empty  

### Disadvantages:
❌ Manual - must remember to switch when changing build config  
❌ Can get out of sync with actual build  

---

## Method 3: Configuration Selector in VS Code UI ✅

VS Code allows you to select different configurations from the UI.

### Setup:

Create multiple configurations in `.vscode/c_cpp_properties.json`:
```json
{
    "configurations": [
        {
            "name": "MPFR Backend",
            "defines": ["ENABLE_HIGH_PRECISION", "USE_MPFR_BACKEND"]
        },
        {
            "name": "Quadmath Backend",
            "defines": ["ENABLE_HIGH_PRECISION", "USE_QUADMATH_BACKEND"]
        }
    ]
}
```

### Usage:

1. Click on the configuration name in the bottom-right status bar
2. Select the desired configuration from the dropdown
3. IntelliSense will update automatically

### Advantages:
✅ Visual - easy to see which configuration is active  
✅ No command line needed  

### Disadvantages:
❌ More verbose configuration file  
❌ Still requires manual switching  

---

## Available Configurations

| Configuration | Defines | Use Case |
|---------------|---------|----------|
| **MPFR Backend** | `ENABLE_HIGH_PRECISION`, `USE_MPFR_BACKEND`, `ENABLE_QUADMATH` | Default, fastest, runtime precision |
| **cpp_dec_float** | `ENABLE_HIGH_PRECISION`, `USE_CPP_DEC_FLOAT_BACKEND`, `ENABLE_QUADMATH` | Boost only, fixed precision |
| **Quadmath** | `ENABLE_HIGH_PRECISION`, `USE_QUADMATH_BACKEND`, `ENABLE_QUADMATH` | Minimal deps, fixed precision |
| **No Templates** | `ENABLE_HIGH_PRECISION`, `USE_MPFR_BACKEND`, `DISABLE_TEMPLATES` | Code duplication instead of templates |
| **Double Only** | (none) | No high-precision support |

---

## Troubleshooting

### IntelliSense shows wrong code blocks

**Solution:** Reload VS Code window
- Press `Ctrl+Shift+P`
- Type "Developer: Reload Window"
- Wait for reindexing to complete

### Undefined references or red squiggles

**Solution 1:** Check include paths in `c_cpp_properties.json`
```json
"includePath": [
    "${workspaceFolder}/include",
    "${workspaceFolder}/build/include",
    "/usr/include",
    "/usr/include/x86_64-linux-gnu"
]
```

**Solution 2:** Regenerate `compile_commands.json`
```bash
cd build
rm -rf *
cmake .. -DENABLE_HIGH_PRECISION=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
```

### Configuration out of sync with build

**Solution:** Use Method 1 (`compile_commands.json`) for automatic sync

---

## Recommended Workflow

**For development:**
1. Use Method 1 (`compile_commands.json`) for automatic sync
2. Regenerate when switching build configurations:
   ```bash
   cd build
   cmake .. -DENABLE_HIGH_PRECISION=ON [flags] -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
   ```
3. Reload VS Code window

**For quick testing:**
1. Use Method 2 (manual defines) with the helper script:
   ```bash
   ./switch_vscode_config.sh mpfr
   ```
2. Reload VS Code window

---

## Files

- `.vscode/c_cpp_properties.json` - Main configuration file
- `.vscode/c_cpp_properties_templates.json` - Template configurations
- `switch_vscode_config.sh` - Helper script to switch configurations
- `compile_commands.json` - Symlink to `build/compile_commands.json`

