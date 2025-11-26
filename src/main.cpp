#include <iostream>
#include <cstring>
#include <cstdlib>
#include "solver.h"
#include "geometry.h"

/**
 * @file main.cpp
 * @brief Main entry point for the polynomial solver application
 */

void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]" << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  --debug                     Enable debug output for geometric operations" << std::endl;
    std::cout << "  --dump [prefix]             Dump detailed geometric information (default prefix: 'dump')" << std::endl;
    std::cout << "  --tolerance <value>         Set box width tolerance (default: 1e-8)" << std::endl;
    std::cout << "  --max-depth <value>         Set maximum subdivision depth (default: 20)" << std::endl;
    std::cout << "  --degeneracy-mult <value>   Set degeneracy multiplier (default: 5.0)" << std::endl;
    std::cout << "  --help                      Show this help message" << std::endl;
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    std::cout << "Polynomial Solver v1.0.0" << std::endl;
    std::cout << "=========================" << std::endl;
    std::cout << std::endl;

    // Parse command-line arguments
    polynomial_solver::SubdivisionConfig config;
    polynomial_solver::GeometryConfig& geom_config = polynomial_solver::getGeometryConfig();

    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--debug") == 0) {
            geom_config.debug = true;
            std::cout << "Debug mode enabled" << std::endl;
        } else if (std::strcmp(argv[i], "--dump") == 0) {
#ifdef ENABLE_GEOMETRY_DUMP
            config.dump_geometry = true;
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                config.dump_prefix = argv[++i];
                std::cout << "Geometry dump enabled with prefix: " << config.dump_prefix << std::endl;
            } else {
                std::cout << "Geometry dump enabled with default prefix: " << config.dump_prefix << std::endl;
            }
#else
            std::cout << "Warning: Geometry dump is not available (compiled out in release mode)" << std::endl;
            // Skip the next argument if it looks like a prefix
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                ++i;
            }
#endif
        } else if (std::strcmp(argv[i], "--tolerance") == 0) {
            if (i + 1 < argc) {
                config.tolerance = std::atof(argv[++i]);
                std::cout << "Box width tolerance set to " << config.tolerance << std::endl;
            } else {
                std::cerr << "Error: --tolerance requires a value" << std::endl;
                return 1;
            }
        } else if (std::strcmp(argv[i], "--max-depth") == 0) {
            if (i + 1 < argc) {
                config.max_depth = static_cast<unsigned int>(std::atoi(argv[++i]));
                std::cout << "Maximum subdivision depth set to " << config.max_depth << std::endl;
            } else {
                std::cerr << "Error: --max-depth requires a value" << std::endl;
                return 1;
            }
        } else if (std::strcmp(argv[i], "--degeneracy-mult") == 0) {
            if (i + 1 < argc) {
                config.degeneracy_multiplier = std::atof(argv[++i]);
                std::cout << "Degeneracy multiplier set to " << config.degeneracy_multiplier << std::endl;
            } else {
                std::cerr << "Error: --degeneracy-mult requires a value" << std::endl;
                return 1;
            }
        } else if (std::strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Error: Unknown option '" << argv[i] << "'" << std::endl;
            print_usage(argv[0]);
            return 1;
        }
    }

    std::cout << std::endl;
    std::cout << "Configuration:" << std::endl;
    std::cout << "  Debug mode: " << (geom_config.debug ? "enabled" : "disabled") << std::endl;
#ifdef ENABLE_GEOMETRY_DUMP
    std::cout << "  Geometry dump: " << (config.dump_geometry ? "enabled" : "disabled") << std::endl;
    if (config.dump_geometry) {
        std::cout << "  Dump prefix: " << config.dump_prefix << std::endl;
    }
#else
    std::cout << "  Geometry dump: disabled (compiled out in release mode)" << std::endl;
#endif
    std::cout << "  Tolerance: " << config.tolerance << std::endl;
    std::cout << "  Max depth: " << config.max_depth << std::endl;
    std::cout << "  Degeneracy multiplier: " << config.degeneracy_multiplier << std::endl;
    std::cout << std::endl;

    // TODO: Implement command-line interface
    // - Read polynomial coefficients
    // - Call solver
    // - Display results

    std::cout << "Application initialized successfully." << std::endl;
    std::cout << "Implementation pending..." << std::endl;

    return 0;
}

