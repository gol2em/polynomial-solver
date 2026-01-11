#include "core/polynomial.h"
#include "solver/solver.h"
#include "core/de_casteljau.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using namespace polynomial_solver;

namespace {

bool approx_equal(double a, double b,
                  double eps = std::numeric_limits<double>::epsilon() * 200.0)
{
    return std::fabs(a - b) <= eps;
}


int test_de_casteljau_subdivide_1d()
{
    // p(t) = t on [0,1] with degree 3 has Bernstein coefficients b_k = k/3.
    std::vector<unsigned int> degrees{3u};
    std::vector<double> bernstein{0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0};

    Polynomial p = Polynomial::fromBernstein(degrees, bernstein);

    const double t_split = 0.3;
    std::vector<double> left_coeffs;
    std::vector<double> right_coeffs;

    DeCasteljau::subdivide1D(bernstein, t_split, left_coeffs, right_coeffs);

    if (left_coeffs.size() != bernstein.size() ||
        right_coeffs.size() != bernstein.size()) {
        std::cerr << "DeCasteljau 1D subdivision: size mismatch" << '\n';
        return 1;
    }

    Polynomial p_left = Polynomial::fromBernstein(degrees, left_coeffs);
    Polynomial p_right = Polynomial::fromBernstein(degrees, right_coeffs);

    const int num_points = 11;
    for (int i = 0; i < num_points; ++i) {
        const double u = static_cast<double>(i) /
                         static_cast<double>(num_points - 1);

        const double t_left = t_split * u;
        const double t_right = t_split + (1.0 - t_split) * u;

        const double expected_left = p.evaluate(t_left);
        const double expected_right = p.evaluate(t_right);

        const double value_left = p_left.evaluate(u);
        const double value_right = p_right.evaluate(u);

        if (!approx_equal(value_left, expected_left)) {
            std::cerr << "DeCasteljau 1D subdivision: left segment mismatch at u="
                      << u << '\n';
            return 1;
        }
        if (!approx_equal(value_right, expected_right)) {
            std::cerr << "DeCasteljau 1D subdivision: right segment mismatch at u="
                      << u << '\n';
            return 1;
        }
    }

    return 0;
}

int test_polynomial_interval_restriction()
{
    // 2D polynomial p(x,y) = 1 + x + 2 y on [0,1]^2.
    const unsigned int dx = 1u;
    const unsigned int dy = 1u;
    std::vector<unsigned int> degrees{dx, dy};

    auto idx = [dy](unsigned int ix, unsigned int iy) {
        return static_cast<std::size_t>(ix * (dy + 1u) + iy);
    };

    std::vector<double> power_coeffs((dx + 1u) * (dy + 1u), 0.0);
    power_coeffs[idx(0u, 0u)] = 1.0; // constant term
    power_coeffs[idx(1u, 0u)] = 1.0; // x term
    power_coeffs[idx(0u, 1u)] = 2.0; // y term

    Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

    // Ensure Bernstein representation is available before using restrictedToInterval
    p.ensureBernsteinPrimary();

    const double a = 0.2;
    const double b = 0.8;

    Polynomial p_x = p.restrictedToInterval(0u, a, b);
    Polynomial p_y = p.restrictedToInterval(1u, a, b);

    const int num_points = 5;

    // Check restriction along x (axis 0).
    for (int iu = 0; iu < num_points; ++iu) {
        const double u = static_cast<double>(iu) /
                         static_cast<double>(num_points - 1);
        for (int iy = 0; iy < num_points; ++iy) {
            const double y = static_cast<double>(iy) /
                             static_cast<double>(num_points - 1);

            const double x_full = a + (b - a) * u;
            std::vector<double> pt_full{x_full, y};
            std::vector<double> pt_restricted{u, y};

            const double expected = p.evaluate(pt_full);
            const double value = p_x.evaluate(pt_restricted);

            if (!approx_equal(value, expected)) {
                std::cerr << "Interval restriction (axis 0) mismatch at (u,y)=("
                          << u << "," << y << ")" << '\n';
                return 1;
            }
        }
    }

    // Check restriction along y (axis 1).
    for (int ix = 0; ix < num_points; ++ix) {
        const double x = static_cast<double>(ix) /
                         static_cast<double>(num_points - 1);
        for (int iu = 0; iu < num_points; ++iu) {
            const double u = static_cast<double>(iu) /
                             static_cast<double>(num_points - 1);

            const double y_full = a + (b - a) * u;
            std::vector<double> pt_full{x, y_full};
            std::vector<double> pt_restricted{x, u};

            const double expected = p.evaluate(pt_full);
            const double value = p_y.evaluate(pt_restricted);

            if (!approx_equal(value, expected)) {
                std::cerr << "Interval restriction (axis 1) mismatch at (x,u)=("
                          << x << "," << u << ")" << '\n';
                return 1;
            }
        }
    }

    return 0;
}

// Build the 1D polynomial corresponding (up to scaling and normalization)
// to (x - 1)(x - 2)...(x - 20) = 0 on x in [1, 20].
//
// We normalize x to t in [0,1] via x = 1 + 19 t, so the roots map to
// t_r = (r - 1) / 19 for r = 1,...,20. Up to a nonzero constant factor,
// this is represented by
//   p(t) = prod_{r=0}^{19} (t - r/19).
Polynomial make_1d_system_polynomial()
{
    const unsigned int degree = 20u;

    // Start from the constant polynomial 1.
    std::vector<double> coeffs(1u, 1.0);

    for (unsigned int r = 0; r < degree; ++r) {
        const double root = static_cast<double>(r) /
                            static_cast<double>(degree - 1u); // r / 19

        const std::size_t old_deg = coeffs.size() - 1u;
        std::vector<double> next(old_deg + 2u, 0.0);

        // Multiply by (t - root).
        for (std::size_t i = 0; i <= old_deg; ++i) {
            next[i + 1u] += coeffs[i];        // t * coeffs[i]
            next[i]      -= root * coeffs[i]; // -root * coeffs[i]
        }

        coeffs.swap(next);
    }

    std::vector<unsigned int> degrees{degree};
    return Polynomial::fromPower(degrees, coeffs);
}

int test_1d_system()
{
    Polynomial p = make_1d_system_polynomial();

    if (p.dimension() != 1u) {
        std::cerr << "1D system: unexpected dimension " << p.dimension() << '\n';
        return 1;
    }

    const std::vector<unsigned int>& degrees = p.degrees();
    if (degrees.size() != 1u || degrees[0] != 20u) {
        std::cerr << "1D system: unexpected degree "
                  << (degrees.empty() ? 0u : degrees[0]) << '\n';
        return 1;
    }

    // Check that the normalized roots t_r = r / 19 are zeros.
    const unsigned int degree_n = degrees[0];
    for (unsigned int r = 0; r < degree_n; ++r) {
        const double t = static_cast<double>(r) /
                         static_cast<double>(degree_n - 1u);
        const double value = p.evaluate(t);
        if (!approx_equal(value, 0.0)) {
            std::cerr << "1D system: root check failed at r=" << r
                      << " (t=" << t << ") value=" << value << '\n';
            return 1;
        }
    }

    // Wrap into a PolynomialSystem and sanity-check.
    PolynomialSystem system({p});
    if (system.dimension() != 1u || system.equationCount() != 1u) {
        std::cerr << "1D system: system shape mismatch" << '\n';
        return 1;
    }

    const std::vector<GraphControlNet> nets = system.graphControlNets();
    if (nets.size() != 1u) {
        std::cerr << "1D system: expected 1 graph control net, got "
                  << nets.size() << '\n';
        return 1;
    }

    if (nets[0].degrees != degrees) {
        std::cerr << "1D system: graph control net degrees mismatch" << '\n';
        return 1;
    }

    return 0;
}

int test_2d_system()
{
    // System in x,y in [0,1]^2:
    //   f1(x,y) = x^2 + y^2 - 1 = 0
    //   f2(x,y) = x^2/4 + 4 y^2 - 1 = 0
    // Their intersection in the first quadrant is at (sqrt(0.8), sqrt(0.2)).

    const unsigned int dx = 2u;
    const unsigned int dy = 2u;
    std::vector<unsigned int> degrees{dx, dy};

    const std::size_t coeff_count = static_cast<std::size_t>((dx + 1u) * (dy + 1u));

    auto idx = [dy](unsigned int ix, unsigned int iy) {
        return static_cast<std::size_t>(ix * (dy + 1u) + iy);
    };

    std::vector<double> power1(coeff_count, 0.0);
    std::vector<double> power2(coeff_count, 0.0);

    // f1(x,y) = x^2 + y^2 - 1
    power1[idx(0u, 0u)] = -1.0; // constant
    power1[idx(2u, 0u)] =  1.0; // x^2
    power1[idx(0u, 2u)] =  1.0; // y^2

    // f2(x,y) = x^2/4 + 4 y^2 - 1
    power2[idx(0u, 0u)] = -1.0;   // constant
    power2[idx(2u, 0u)] =  0.25;  // x^2 / 4
    power2[idx(0u, 2u)] =  4.0;   // 4 y^2

    Polynomial f1 = Polynomial::fromPower(degrees, power1);
    Polynomial f2 = Polynomial::fromPower(degrees, power2);

    PolynomialSystem system({f1, f2});

    if (system.dimension() != 2u || system.equationCount() != 2u) {
        std::cerr << "2D system: system shape mismatch" << '\n';
        return 1;
    }

    // Check the known intersection point in [0,1]^2.
    const double x = std::sqrt(0.8); // sqrt(4/5)
    const double y = std::sqrt(0.2); // sqrt(1/5)
    std::vector<double> pt{x, y};

    const double v1 = f1.evaluate(pt);
    const double v2 = f2.evaluate(pt);

    if (!approx_equal(v1, 0.0)) {
        std::cerr << "2D system: f1(x,y) != 0 at intersection, value=" << v1 << '\n';
        return 1;
    }
    if (!approx_equal(v2, 0.0)) {
        std::cerr << "2D system: f2(x,y) != 0 at intersection, value=" << v2 << '\n';
        return 1;
    }

    const std::vector<GraphControlNet> nets = system.graphControlNets();
    if (nets.size() != 2u) {
        std::cerr << "2D system: expected 2 graph control nets, got "
                  << nets.size() << '\n';
        return 1;
    }

    if (nets[0].degrees != degrees || nets[1].degrees != degrees) {
        std::cerr << "2D system: graph control net degrees mismatch" << '\n';
        return 1;
    }

    return 0;
}

} // namespace

int main()
{
    int status = 0;
    status |= test_de_casteljau_subdivide_1d();
    status |= test_polynomial_interval_restriction();
    status |= test_1d_system();
    status |= test_2d_system();

    if (status == 0) {
        std::cout << "Polynomial system example tests passed." << std::endl;
    }

    return status;
}

