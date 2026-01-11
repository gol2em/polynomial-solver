#include "core/polynomial.h"
#include "solver/solver.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

using namespace polynomial_solver;

namespace {

bool approx_equal(double a, double b,
                  double eps = std::numeric_limits<double>::epsilon() * 100.0) {
    return std::fabs(a - b) <= eps;

}

void dump_graph_data_1d(const Polynomial& p,
                        const std::vector<double>& control_points,
                        unsigned int n)
{
    const std::size_t dim = p.dimension();
    const std::size_t stride = dim + 1u;

    // 1D control points: (x, y) = (k/n, Bernstein coefficient).
    std::ofstream ctrl("graph_1d_control_points.csv");
    ctrl << "x,y\n";
    for (unsigned int k = 0; k <= n; ++k) {
        const std::size_t base = static_cast<std::size_t>(k) * stride;
        const double x = control_points[base + 0u];
        const double y = control_points[base + 1u];
        ctrl << x << "," << y << "\n";
    }

    // Dense sampling of the curve using the Polynomial evaluation.
    std::ofstream curve("graph_1d_curve.csv");
    curve << "t,f\n";
    const int num_samples = 201;
    for (int j = 0; j < num_samples; ++j) {
        const double t = static_cast<double>(j) /
                         static_cast<double>(num_samples - 1);
        const double f = p.evaluate(t);
        curve << t << "," << f << "\n";
    }
}

void dump_graph_data_2d(const Polynomial& p,
                        const std::vector<double>& control_points,
                        unsigned int dx,
                        unsigned int dy)
{
    const std::size_t dim = p.dimension();
    const std::size_t stride = dim + 1u;

    auto index = [dy](unsigned int i, unsigned int j) {
        return static_cast<std::size_t>(i * (dy + 1u) + j);
    };

    // 2D control points: (x, y, z) where z is the Bernstein coefficient.
    std::ofstream ctrl("graph_2d_control_points.csv");
    ctrl << "x,y,z\n";
    for (unsigned int i = 0; i <= dx; ++i) {
        for (unsigned int j = 0; j <= dy; ++j) {
            const std::size_t idx = index(i, j);
            const std::size_t base = idx * stride;
            const double x = control_points[base + 0u];
            const double y = control_points[base + 1u];
            const double z = control_points[base + 2u];
            ctrl << x << "," << y << "," << z << "\n";
        }
    }

    // Dense sampling of the surface using the Polynomial evaluation.
    std::ofstream surf("graph_2d_surface.csv");
    surf << "x,y,z\n";
    const int grid_n = 41;
    for (int ix = 0; ix < grid_n; ++ix) {
        for (int iy = 0; iy < grid_n; ++iy) {
            const double tx = static_cast<double>(ix) /
                              static_cast<double>(grid_n - 1);
            const double ty = static_cast<double>(iy) /
                              static_cast<double>(grid_n - 1);
            std::vector<double> params{tx, ty};
            const double z = p.evaluate(params);
            surf << tx << "," << ty << "," << z << "\n";
        }
    }
}

int test_graph_control_points_1d() {
    const unsigned int n = 3u;
    std::vector<unsigned int> degrees{n};

    // p(t) = 1 + 2 t + 3 t^2 + 4 t^3
    std::vector<double> power_coeffs{1.0, 2.0, 3.0, 4.0};

    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

    std::vector<double> control_points;
    p.graphControlPoints(control_points);

    const std::vector<double>& b = p.bernsteinCoefficients();
    const std::size_t dim = p.dimension();
    const std::size_t stride = dim + 1u;

    if (dim != 1u) {
        std::cerr << "1D graph test: unexpected dimension " << dim << '\n';
        return 1;
    }

    if (control_points.size() != b.size() * stride) {
        std::cerr << "1D graph test: unexpected control_points size "
                  << control_points.size() << ", expected "
                  << b.size() * stride << '\n';
        return 1;
    }

    // Check that x-coordinates are k/n and y-coordinates equal Bernstein coeffs.
    for (unsigned int k = 0; k <= n; ++k) {
        const std::size_t base = static_cast<std::size_t>(k) * stride;
        const double x = control_points[base + 0u];
        const double y = control_points[base + 1u];

        const double expected_x = (n > 0u)
                                      ? static_cast<double>(k) /
                                            static_cast<double>(n)
                                      : 0.0;
        const double expected_y = b[static_cast<std::size_t>(k)];

        if (!approx_equal(x, expected_x)) {
            std::cerr << "1D graph test: x-coordinate mismatch at k=" << k
                      << ": expected " << expected_x << ", got " << x << '\n';
            return 1;
        }
        if (!approx_equal(y, expected_y)) {
            std::cerr << "1D graph test: y-coordinate mismatch at k=" << k
                      << ": expected " << expected_y << ", got " << y << '\n';
            return 1;
        }
    }

    // Build a polynomial from the x-coordinates and verify it is the identity map.
    std::vector<double> x_coeffs(b.size());
    for (std::size_t k = 0; k < b.size(); ++k) {
        x_coeffs[k] = control_points[k * stride + 0u];
    }

    Polynomial x_poly = Polynomial::fromBernstein(degrees, x_coeffs);

    const int num_points = 11;
    for (int j = 0; j < num_points; ++j) {
        const double t = static_cast<double>(j) /
                         static_cast<double>(num_points - 1);
        const double x_t = x_poly.evaluate(t);
        if (!approx_equal(x_t, t)) {
            std::cerr << "1D graph test: identity check failed at t=" << t
                      << ": expected " << t << ", got " << x_t << '\n';
            return 1;
        }
    }

    // Verify that the PolynomialSystem wrapper forwards correctly.
    PolynomialSystem system({p});
    const std::vector<GraphControlNet> nets = system.graphControlNets();
    if (nets.size() != 1u) {
        std::cerr << "1D graph test: system produced " << nets.size()
                  << " nets, expected 1\n";
        return 1;
    }
    if (nets[0].degrees != degrees) {
        std::cerr << "1D graph test: system degrees mismatch\n";
        return 1;
    }
    if (nets[0].control_points != control_points) {
        std::cerr << "1D graph test: system control_points mismatch\n";
        return 1;
    }

    dump_graph_data_1d(p, control_points, n);


    return 0;
}

int test_graph_control_points_2d() {
    const unsigned int dx = 2u;
    const unsigned int dy = 3u;
    std::vector<unsigned int> degrees{dx, dy};

    const std::size_t coeff_count =
        static_cast<std::size_t>((dx + 1u) * (dy + 1u));
    std::vector<double> power_coeffs(coeff_count, 0.0);

    auto index = [dy](unsigned int i, unsigned int j) {
        return static_cast<std::size_t>(i * (dy + 1u) + j);
    };

    // Simple pattern: a_{i,j} = i + 10*j.
    for (unsigned int i = 0; i <= dx; ++i) {
        for (unsigned int j = 0; j <= dy; ++j) {
            power_coeffs[index(i, j)] =
                static_cast<double>(i + 10u * j);
        }
    }

    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

    std::vector<double> control_points;
    p.graphControlPoints(control_points);

    const std::vector<double>& b = p.bernsteinCoefficients();
    const std::size_t dim = p.dimension();
    const std::size_t stride = dim + 1u;

    if (dim != 2u) {
        std::cerr << "2D graph test: unexpected dimension " << dim << '\n';
        return 1;
    }

    if (control_points.size() != b.size() * stride) {
        std::cerr << "2D graph test: unexpected control_points size "
                  << control_points.size() << ", expected "
                  << b.size() * stride << '\n';
        return 1;
    }

    // Check coordinates for each multi-index (i,j).
    for (unsigned int i = 0; i <= dx; ++i) {
        for (unsigned int j = 0; j <= dy; ++j) {
            const std::size_t idx = index(i, j);
            const std::size_t base = idx * stride;

            const double x = control_points[base + 0u];
            const double y = control_points[base + 1u];
            const double z = control_points[base + 2u];

            const double expected_x = (dx > 0u)
                                          ? static_cast<double>(i) /
                                                static_cast<double>(dx)
                                          : 0.0;
            const double expected_y = (dy > 0u)
                                          ? static_cast<double>(j) /
                                                static_cast<double>(dy)
                                          : 0.0;
            const double expected_z = b[idx];

            if (!approx_equal(x, expected_x)) {
                std::cerr << "2D graph test: x-coordinate mismatch at (" << i
                          << "," << j << ") expected " << expected_x
                          << ", got " << x << '\n';
                return 1;
            }
            if (!approx_equal(y, expected_y)) {
                std::cerr << "2D graph test: y-coordinate mismatch at (" << i
                          << "," << j << ") expected " << expected_y
                          << ", got " << y << '\n';
                return 1;
            }
            if (!approx_equal(z, expected_z)) {
                std::cerr << "2D graph test: z-coordinate mismatch at (" << i
                          << "," << j << ") expected " << expected_z
                          << ", got " << z << '\n';
                return 1;
            }
        }
    }

    // Build polynomials from x and y coordinates and from z, and verify
    // that they represent the identity maps and the original polynomial.
    std::vector<double> x_coeffs(b.size());
    std::vector<double> y_coeffs(b.size());
    std::vector<double> z_coeffs(b.size());

    for (unsigned int i = 0; i <= dx; ++i) {
        for (unsigned int j = 0; j <= dy; ++j) {
            const std::size_t idx = index(i, j);
            const std::size_t base = idx * stride;
            x_coeffs[idx] = control_points[base + 0u];
            y_coeffs[idx] = control_points[base + 1u];
            z_coeffs[idx] = control_points[base + 2u];
        }
    }

    Polynomial x_poly = Polynomial::fromBernstein(degrees, x_coeffs);
    Polynomial y_poly = Polynomial::fromBernstein(degrees, y_coeffs);
    Polynomial z_poly = Polynomial::fromBernstein(degrees, z_coeffs);

    const int grid_n = 5;
    for (int ix = 0; ix < grid_n; ++ix) {
        for (int iy = 0; iy < grid_n; ++iy) {
            const double tx = static_cast<double>(ix) /
                              static_cast<double>(grid_n - 1);
            const double ty = static_cast<double>(iy) /
                              static_cast<double>(grid_n - 1);
            std::vector<double> params{tx, ty};

            const double x_t = x_poly.evaluate(params);
            const double y_t = y_poly.evaluate(params);
            const double z_t = z_poly.evaluate(params);
            const double p_t = p.evaluate(params);

            if (!approx_equal(x_t, tx)) {
                std::cerr << "2D graph test: x identity failed at (" << tx
                          << "," << ty << ") expected " << tx
                          << ", got " << x_t << '\n';
                return 1;
            }
            if (!approx_equal(y_t, ty)) {
                std::cerr << "2D graph test: y identity failed at (" << tx
                          << "," << ty << ") expected " << ty
                          << ", got " << y_t << '\n';
                return 1;
            }
            if (!approx_equal(z_t, p_t)) {
                std::cerr << "2D graph test: value mismatch at (" << tx
                          << "," << ty << ") expected " << p_t
                          << ", got " << z_t << '\n';
                return 1;
            }
        }
    }
    dump_graph_data_2d(p, control_points, dx, dy);



    return 0;
}

} // namespace

int main() {
    int status = 0;
    status |= test_graph_control_points_1d();
    status |= test_graph_control_points_2d();
    return status;
}

