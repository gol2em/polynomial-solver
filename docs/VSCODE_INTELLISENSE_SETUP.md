# VS Code IntelliSense Setup for High-Precision Code

This document explains how VS Code IntelliSense automatically detects the correct preprocessor macros and highlights code blocks correctly.

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

## Solution: Use `compile_commands.json` ✅

**The project uses the standard CMake approach with `compile_commands.json`.**

This is IDE-agnostic and works for:
- ✅ VS Code
- ✅ CLion
- ✅ Qt Creator
- ✅ Any IDE that supports `compile_commands.json`

### How It Works

1. **CMake generates `compile_commands.json`** with all compiler flags and macros
2. **VS Code reads this file** to understand the build configuration
3. **IntelliSense automatically highlights** the correct code blocks

**No manual configuration needed!**

---

## Setup (One-Time)

### Step 1: Configure VS Code

The repository includes `.vscode/c_cpp_properties.json` that tells VS Code to use `compile_commands.json`:

```json
{
    "configurations": [{
        "name": "Linux",
        "compileCommands": "${workspaceFolder}/compile_commands.json"
    }],
    "version": 4
}
```

**This file is already committed to the repository, so you don't need to do anything!**

### Step 2: Run CMake

CMake automatically exports `compile_commands.json`:

```bash
cd build
cmake .. -DENABLE_HIGH_PRECISION=ON
```

CMake will:
- ✅ Export `build/compile_commands.json` with all compiler flags
- ✅ Create symlink `compile_commands.json` → `build/compile_commands.json`
- ✅ Include all preprocessor macros (ENABLE_HIGH_PRECISION, USE_MPFR_BACKEND, etc.)
- ✅ Include all include paths

### Step 3: Reload VS Code

After CMake finishes:
1. Press `Ctrl+Shift+P`
2. Type "Developer: Reload Window"
3. Wait for IntelliSense to reindex (watch bottom-right status bar)

**Done!** Code blocks are now correctly highlighted based on your build configuration.

---

## Switching Configurations

To switch between different backends or configurations, simply re-run CMake with different flags:

```bash
# Switch to MPFR backend
cd build && cmake .. -DENABLE_HIGH_PRECISION=ON

# Switch to cpp_dec_float backend
cd build && cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_MPFR=OFF -DUSE_GMP=OFF

# Switch to Quadmath backend
cd build && cmake .. -DENABLE_HIGH_PRECISION=ON -DUSE_BOOST=OFF -DUSE_MPFR=OFF -DUSE_GMP=OFF

# Disable templates
cd build && cmake .. -DENABLE_HIGH_PRECISION=ON -DDISABLE_TEMPLATES=ON

# Disable high-precision
cd build && cmake .. -DENABLE_HIGH_PRECISION=OFF
```

After each CMake run:
1. CMake automatically regenerates `compile_commands.json`
2. Reload VS Code: `Ctrl+Shift+P` → "Developer: Reload Window"
3. IntelliSense will automatically use the new configuration

---

## What's in `compile_commands.json`

The `compile_commands.json` file contains the exact compiler commands used to build each source file, including:

### Preprocessor Defines
All macros defined by CMake:
- `ENABLE_HIGH_PRECISION` (if enabled)
- `USE_MPFR_BACKEND` / `USE_CPP_DEC_FLOAT_BACKEND` / `USE_QUADMATH_BACKEND`
- `ENABLE_QUADMATH` (if available)
- `USE_TEMPLATES` (if templates enabled)
- `DISABLE_TEMPLATES` (if templates disabled)

### Include Paths
All include directories:
- `${workspaceFolder}/include`
- `${workspaceFolder}/build/include`
- Detected library paths (e.g., `/usr/include`, `/usr/include/x86_64-linux-gnu`)

### Compiler Settings
- Compiler path and version
- C++ standard (C++11)
- Optimization flags
- Warning flags

VS Code reads this file and automatically configures IntelliSense to match your build configuration.

---

## Example: What VS Code Sees

When you build with MPFR backend, `compile_commands.json` contains entries like:

```json
{
  "directory": "/home/user/polynomial-solver/build",
  "command": "/usr/bin/c++ -DENABLE_HIGH_PRECISION -DUSE_MPFR_BACKEND -DUSE_TEMPLATES -I/home/user/polynomial-solver/include -I/home/user/polynomial-solver/build/include -I/usr/include -std=c++11 -o polynomial.cpp.o -c /home/user/polynomial-solver/src/polynomial.cpp",
  "file": "/home/user/polynomial-solver/src/polynomial.cpp"
}
```

VS Code parses this and knows:
- ✅ `ENABLE_HIGH_PRECISION` is defined → highlight high-precision code
- ✅ `USE_MPFR_BACKEND` is defined → highlight MPFR-specific code
- ✅ `USE_TEMPLATES` is defined → highlight template code
- ✅ Include paths are `/home/user/polynomial-solver/include`, etc.

This is why code blocks are correctly highlighted!

---

## Troubleshooting

### IntelliSense shows wrong code blocks

**Solution:** Reload VS Code window
1. Press `Ctrl+Shift+P`
2. Type "Developer: Reload Window"
3. Wait for IntelliSense to reindex

### Configuration not updating after CMake

**Possible causes:**
1. CMake didn't finish successfully → Check CMake output for errors
2. VS Code hasn't reloaded → Reload window (`Ctrl+Shift+P` → "Developer: Reload Window")
3. IntelliSense is still indexing → Wait for indexing to complete (bottom-right status bar)

### Wrong backend detected

**Solution:** Clean build and reconfigure
```bash
cd build
rm -rf *
cmake .. -DENABLE_HIGH_PRECISION=ON [your flags]
```

Then reload VS Code.

### compile_commands.json not found

**Solution:** Symlink should be created automatically, but you can create it manually:
```bash
ln -sf build/compile_commands.json compile_commands.json
```

---

## Benefits of This Approach

✅ **IDE-agnostic** - Works with VS Code, CLion, Qt Creator, Vim, Emacs, etc.
✅ **Standard approach** - Uses CMake's built-in `CMAKE_EXPORT_COMPILE_COMMANDS`
✅ **Always in sync** - Configuration matches your build exactly
✅ **No manual editing** - CMake handles everything
✅ **Backend-aware** - Correct macros for MPFR/cpp_dec_float/quadmath
✅ **Template-aware** - Knows if templates are enabled/disabled
✅ **Include paths** - Automatically detects all library paths
✅ **Easy switching** - Just re-run CMake with different flags

---

## For Other IDEs

The same `compile_commands.json` file works with:

- **CLion** - Automatically detected
- **Qt Creator** - Import via "Open Project" → select `compile_commands.json`
- **Vim/Neovim** - Use with clangd LSP
- **Emacs** - Use with lsp-mode or eglot
- **Sublime Text** - Use with LSP plugin

All these IDEs will correctly highlight code blocks based on the same `compile_commands.json` file!
