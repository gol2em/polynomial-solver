#ifdef ENABLE_HIGH_PRECISION

#include "polynomial_hp.h"
#include "precision_conversion.h"
#include <stdexcept>

namespace polynomial_solver {

// Helper: De Casteljau evaluation for 1D high-precision
static mpreal evaluate1D_HP(const std::vector<mpreal>& bernstein_coeffs, const mpreal& t) {
    const std::size_t n = bernstein_coeffs.size();
    if (n == 0) {
        return mpreal(0);
    }

    // Standard De Casteljau algorithm (univariate)
    std::vector<mpreal> tmp(bernstein_coeffs.begin(), bernstein_coeffs.end());

    mpreal one_minus_t = mpreal(1) - t;
    for (std::size_t r = 1; r < n; ++r) {
        const std::size_t upper = n - r;
        for (std::size_t i = 0; i < upper; ++i) {
            tmp[i] = one_minus_t * tmp[i] + t * tmp[i + 1];
        }
    }

    return tmp[0];
}

// Helper: De Casteljau evaluation for tensor-product high-precision
static mpreal evaluateTensorProduct_HP(
    const std::vector<unsigned int>& degrees,
    const std::vector<mpreal>& bernstein_coeffs,
    const std::vector<mpreal>& parameters)
{
    const std::size_t dim = degrees.size();

    // Univariate case
    if (dim == 1 && parameters.size() == 1) {
        const unsigned int deg = degrees[0];
        const std::size_t count = static_cast<std::size_t>(deg + 1);
        if (bernstein_coeffs.size() != count) {
            return mpreal(0);
        }
        return evaluate1D_HP(bernstein_coeffs, parameters[0]);
    }

    // Bivariate tensor-product case
    if (dim == 2 && parameters.size() == 2) {
        const unsigned int deg_x = degrees[0];
        const unsigned int deg_y = degrees[1];
        const std::size_t nx = static_cast<std::size_t>(deg_x + 1);
        const std::size_t ny = static_cast<std::size_t>(deg_y + 1);

        if (bernstein_coeffs.size() != nx * ny) {
            return mpreal(0);
        }

        const mpreal& t_x = parameters[0];
        const mpreal& t_y = parameters[1];

        // First, evaluate along y-direction for each x-slice
        std::vector<mpreal> y_vals(nx);
        for (std::size_t i = 0; i < nx; ++i) {
            std::vector<mpreal> y_coeffs(ny);
            for (std::size_t j = 0; j < ny; ++j) {
                y_coeffs[j] = bernstein_coeffs[i * ny + j];
            }
            y_vals[i] = evaluate1D_HP(y_coeffs, t_y);
        }

        // Then evaluate along x-direction
        return evaluate1D_HP(y_vals, t_x);
    }

    // Higher dimensions not implemented
    throw std::runtime_error("PolynomialHP: evaluateTensorProduct only supports 1D and 2D");
}

// PolynomialHP implementation

PolynomialHP::PolynomialHP()
    : degrees_(), bernstein_coeffs_()
{
}

PolynomialHP::PolynomialHP(const std::vector<unsigned int>& degrees,
                           const std::vector<mpreal>& bernstein_coeffs)
    : degrees_(degrees), bernstein_coeffs_(bernstein_coeffs)
{
}

PolynomialHP::PolynomialHP(const Polynomial& poly)
    : degrees_(poly.degrees())
{
    // Convert double coefficients to high-precision
    const std::vector<double>& double_coeffs = poly.bernsteinCoefficients();
    bernstein_coeffs_.reserve(double_coeffs.size());
    for (double coeff : double_coeffs) {
        bernstein_coeffs_.push_back(toHighPrecision(coeff));
    }
}

PolynomialHP::~PolynomialHP() {
}

std::size_t PolynomialHP::dimension() const {
    return degrees_.size();
}

const std::vector<unsigned int>& PolynomialHP::degrees() const {
    return degrees_;
}

std::size_t PolynomialHP::coefficientCount() const {
    return bernstein_coeffs_.size();
}

const std::vector<mpreal>& PolynomialHP::bernsteinCoefficients() const {
    return bernstein_coeffs_;
}

mpreal PolynomialHP::evaluate(const std::vector<mpreal>& parameters) const {
    return evaluateTensorProduct_HP(degrees_, bernstein_coeffs_, parameters);
}

mpreal PolynomialHP::evaluate(const mpreal& t) const {
    std::vector<mpreal> params(1, t);
    return evaluate(params);
}

bool PolynomialHP::empty() const {
    return bernstein_coeffs_.empty();
}

// Conversion utility

PolynomialHP convertToHighPrecision(const Polynomial& poly) {
    return PolynomialHP(poly);
}

} // namespace polynomial_solver

#endif // ENABLE_HIGH_PRECISION

