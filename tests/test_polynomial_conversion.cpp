#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <random>


#include "polynomial.h"

using polynomial_solver::Polynomial;

namespace {

bool approx_equal(double a, double b,
                  double eps = std::numeric_limits<double>::epsilon() * 100.0)
{
    return std::fabs(a - b) <= eps;
}

} // anonymous namespace

int main()
{
    {
        // 1D: p(t) = 1 + 2 t + 3 t^2 on [0,1].
        std::vector<unsigned int> degrees{2u};
        std::vector<double> power_coeffs{1.0, 2.0, 3.0};

        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

        // Evaluate at a few points and compare with direct power evaluation.
        const double ts[] = {0.0, 0.25, 0.5, 0.75, 1.0};
        for (double t : ts) {
            double expected = power_coeffs[0]
                              + power_coeffs[1] * t
                              + power_coeffs[2] * t * t;
            double value = p.evaluate(t);
            if (!approx_equal(expected, value)) {
                std::cerr << "1D conversion test (values) failed at t=" << t
                          << ": expected " << expected
                          << ", got " << value << '\n';
                return 1;
            }
        }
    }

    {
        // 2D: p(x,y) = 1 + x + y on [0,1]^2 (value check).
        std::vector<unsigned int> degrees{1u, 1u};
        // Coefficients in tensor-product power basis with last dimension (y) fastest.
        std::vector<double> power_coeffs{1.0, 1.0, 1.0, 0.0};

        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

        const double xs[] = {0.0, 0.3, 0.7, 1.0};
        const double ys[] = {0.0, 0.2, 0.8, 1.0};

        for (double x : xs) {
            for (double y : ys) {
                double expected = 1.0 + x + y;
                std::vector<double> pt{x, y};
                double value = p.evaluate(pt);
                if (!approx_equal(expected, value)) {
                    std::cerr << "2D conversion test (values) failed at (" << x << "," << y
                              << "): expected " << expected
                              << ", got " << value << '\n';
                    return 1;
                }
            }
        }
    }

    {
        // 1D coefficient check: p(t) = 1 with degree 7.
        const unsigned int n = 7u;
        std::vector<unsigned int> degrees{n};
        std::vector<double> power_coeffs(n + 1u, 0.0);
        power_coeffs[0] = 1.0;

        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        p.ensureBernsteinPrimary();  // Explicitly convert to Bernstein
        const std::vector<double>& b = p.bernsteinCoefficients();

        if (b.size() != power_coeffs.size()) {
            std::cerr << "1D coefficient test (constant) size mismatch: expected "
                      << power_coeffs.size() << ", got " << b.size() << '\n';
            return 1;
        }

        for (std::size_t k = 0; k < b.size(); ++k) {
            if (!approx_equal(b[k], 1.0)) {
                std::cerr << "1D coefficient test (constant) failed at k=" << k
                          << ": expected 1.0, got " << b[k] << '\n';
                return 1;
            }
        }
    }

    {
        // 1D coefficient check: p(t) = t with degree 10.
        const unsigned int n = 10u;
        std::vector<unsigned int> degrees{n};
        std::vector<double> power_coeffs(n + 1u, 0.0);
        power_coeffs[1] = 1.0;

        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        p.ensureBernsteinPrimary();  // Explicitly convert to Bernstein
        const std::vector<double>& b = p.bernsteinCoefficients();

        if (b.size() != power_coeffs.size()) {
            std::cerr << "1D coefficient test (linear) size mismatch: expected "
                      << power_coeffs.size() << ", got " << b.size() << '\n';
            return 1;
        }

        for (unsigned int k = 0; k <= n; ++k) {
            double expected = static_cast<double>(k) / static_cast<double>(n);
            if (!approx_equal(b[k], expected)) {
                std::cerr << "1D coefficient test (linear) failed at k=" << k
                          << ": expected " << expected << ", got " << b[k] << '\n';
                return 1;
            }
        }
    }

    {
        // 1D coefficient check: p(t) = t^2 with degree 8.
        const unsigned int n = 8u;
        std::vector<unsigned int> degrees{n};
        std::vector<double> power_coeffs(n + 1u, 0.0);
        power_coeffs[2] = 1.0;

        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        p.ensureBernsteinPrimary();  // Explicitly convert to Bernstein
        const std::vector<double>& b = p.bernsteinCoefficients();

        if (b.size() != power_coeffs.size()) {
            std::cerr << "1D coefficient test (quadratic) size mismatch: expected "
                      << power_coeffs.size() << ", got " << b.size() << '\n';
            return 1;
        }

        for (unsigned int k = 0; k <= n; ++k) {
            double expected = 0.0;
            if (k >= 2u) {
                expected = (static_cast<double>(k) * static_cast<double>(k - 1)) /
                           (static_cast<double>(n) * static_cast<double>(n - 1));
            }
            if (!approx_equal(b[k], expected)) {
                std::cerr << "1D coefficient test (quadratic) failed at k=" << k
                          << ": expected " << expected << ", got " << b[k] << '\n';
                return 1;
            }
        }
    }

    {
        // 2D coefficient check: p(x,y) = 1 + x + y on [0,1]^2.
        const unsigned int dx = 3u;
        const unsigned int dy = 4u;
        std::vector<unsigned int> degrees{dx, dy};
        std::vector<double> power_coeffs((dx + 1u) * (dy + 1u), 0.0);

        auto idx = [dy](unsigned int ix, unsigned int iy) {
            return static_cast<std::size_t>(ix * (dy + 1u) + iy);
        };

        power_coeffs[idx(0u, 0u)] = 1.0; // constant term
        power_coeffs[idx(1u, 0u)] = 1.0; // x term
        power_coeffs[idx(0u, 1u)] = 1.0; // y term

        Polynomial p = Polynomial::fromPower(degrees, power_coeffs);
        p.ensureBernsteinPrimary();  // Explicitly convert to Bernstein
        const std::vector<double>& b = p.bernsteinCoefficients();

        if (b.size() != power_coeffs.size()) {
            std::cerr << "2D coefficient test size mismatch: expected "
                      << power_coeffs.size() << ", got " << b.size() << '\n';
            return 1;
        }

        for (unsigned int ix = 0; ix <= dx; ++ix) {
            for (unsigned int iy = 0; iy <= dy; ++iy) {
                const std::size_t linear = idx(ix, iy);
                double expected = 1.0
                                  + static_cast<double>(ix) / static_cast<double>(dx)
                                  + static_cast<double>(iy) / static_cast<double>(dy);
                if (!approx_equal(b[linear], expected)) {
                    std::cerr << "2D coefficient test failed at (" << ix << "," << iy
                              << "): expected " << expected << ", got " << b[linear]
                              << '\n';
                    return 1;
                }
            }
        }
    }

    {
        // 1D stress test: random polynomials up to degree 12.
        const unsigned int n = 12u;
        std::vector<unsigned int> degrees{n};

        std::mt19937 rng(123456u);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);

        const int num_polys = 10;
        const int num_points = 101;

        for (int p_idx = 0; p_idx < num_polys; ++p_idx) {
            std::vector<double> power_coeffs(n + 1u);
            for (unsigned int i = 0; i <= n; ++i) {
                power_coeffs[i] = dist(rng);
            }

            Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

            for (int j = 0; j < num_points; ++j) {
                double t = static_cast<double>(j) /
                           static_cast<double>(num_points - 1);

                // Horner evaluation in power basis.
                double expected = 0.0;
                for (int i = static_cast<int>(n); i >= 0; --i) {
                    expected = expected * t +
                               power_coeffs[static_cast<std::size_t>(i)];
                }

                double value = p.evaluate(t);
                if (!approx_equal(expected, value)) {
                    std::cerr << "1D random stress test failed for polynomial "
                              << p_idx << " at t=" << t
                              << ": expected " << expected
                              << ", got " << value << '\n';
                    return 1;
                }
            }
        }
    }

    {
        // 2D stress test: random polynomials with degrees (6,6).
        const unsigned int dx = 6u;
        const unsigned int dy = 6u;
        std::vector<unsigned int> degrees{dx, dy};

        std::mt19937 rng(7891011u);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);

        const int num_polys = 5;
        const int num_points = 21;

        auto idx2d = [dy](unsigned int ix, unsigned int iy) {
            return static_cast<std::size_t>(ix * (dy + 1u) + iy);
        };

        for (int p_idx = 0; p_idx < num_polys; ++p_idx) {
            std::vector<double> power_coeffs((dx + 1u) * (dy + 1u));
            for (unsigned int ix = 0; ix <= dx; ++ix) {
                for (unsigned int iy = 0; iy <= dy; ++iy) {
                    power_coeffs[idx2d(ix, iy)] = dist(rng);
                }
            }

            Polynomial p = Polynomial::fromPower(degrees, power_coeffs);

            for (int ix_sample = 0; ix_sample < num_points; ++ix_sample) {
                double x = static_cast<double>(ix_sample) /
                           static_cast<double>(num_points - 1);
                for (int iy_sample = 0; iy_sample < num_points; ++iy_sample) {
                    double y = static_cast<double>(iy_sample) /
                               static_cast<double>(num_points - 1);

                    // Horner evaluation in y, then x.
                    std::vector<double> row(dx + 1u);
                    for (unsigned int ix = 0; ix <= dx; ++ix) {
                        double v = 0.0;
                        for (int iy = static_cast<int>(dy); iy >= 0; --iy) {
                            v = v * y +
                                power_coeffs[idx2d(ix, static_cast<unsigned int>(iy))];
                        }
                        row[ix] = v;
                    }

                    double expected = 0.0;
                    for (int ix = static_cast<int>(dx); ix >= 0; --ix) {
                        expected = expected * x +
                                   row[static_cast<std::size_t>(ix)];
                    }

                    std::vector<double> pt{x, y};
                    double value = p.evaluate(pt);
                    if (!approx_equal(expected, value)) {
                        std::cerr << "2D random stress test failed for polynomial "
                                  << p_idx << " at (x,y) = (" << x << "," << y
                                  << "): expected " << expected
                                  << ", got " << value << '\n';
                        return 1;
                    }
                }
            }
        }
    }


    std::cout << "Polynomial power-to-Bernstein conversion tests passed." << std::endl;
    return 0;
}

