#include "paw/coulomb_correction.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace dft
{

namespace
{

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kSqrtFourPi = 3.5449077018110318;
constexpr double kOriginTolerance = 1e-14;

int dense_matrix_extent(std::size_t value)
{
    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("DenseMatrix extent exceeds int range");
    }

    return static_cast<int>(value);
}

DenseMatrix make_dense_matrix(std::size_t row_count, std::size_t column_count)
{
    return DenseMatrix(typename DenseMatrix::ShapeType{dense_matrix_extent(row_count), dense_matrix_extent(column_count)});
}

std::size_t dense_matrix_row_count(const DenseMatrix &matrix)
{
    return static_cast<std::size_t>(matrix.shape()[0]);
}

std::size_t dense_matrix_column_count(const DenseMatrix &matrix)
{
    return static_cast<std::size_t>(matrix.shape()[1]);
}

std::vector<double> trapezoid_weights(const std::vector<double> &radii)
{
    if (radii.size() < 2)
    {
        throw std::runtime_error("Need at least two radial grid points for quadrature");
    }

    std::vector<double> weights(radii.size(), 0.0);
    weights.front() = 0.5 * (radii[1] - radii[0]);

    for (std::size_t index = 1; index + 1 < radii.size(); ++index)
    {
        weights[index] = 0.5 * (radii[index + 1] - radii[index - 1]);
    }

    weights.back() = 0.5 * (radii.back() - radii[radii.size() - 2]);
    return weights;
}

double radial_moment_integral(const RadialFunction &radial_function)
{
    const std::vector<double> weights = trapezoid_weights(radial_function.radii);

    double integral = 0.0;
    for (std::size_t index = 0; index < radial_function.size(); ++index)
    {
        integral +=
            weights[index] * radial_function.radii[index] * radial_function.radii[index] * radial_function.values[index];
    }

    return integral;
}

RadialFunction make_normalized_shape_function_00(const XmlAttributeMap &shape_attributes,
                                                 const std::vector<double> &radial_grid)
{
    if (!shape_attributes.has("type") || !shape_attributes.has("rc"))
    {
        throw std::runtime_error("shape_function is missing required attributes");
    }

    const std::string type = shape_attributes.get_string("type");
    if (type != "sinc")
    {
        throw std::runtime_error("Only sinc shape_function is currently implemented for static Coulomb correction");
    }

    const double cutoff_radius = shape_attributes.get_double("rc");

    RadialFunction shape_function;
    shape_function.radii = radial_grid;
    shape_function.values.resize(radial_grid.size(), 0.0);

    for (std::size_t index = 0; index < radial_grid.size(); ++index)
    {
        const double radius = radial_grid[index];
        if (radius > cutoff_radius)
        {
            shape_function.values[index] = 0.0;
            continue;
        }

        const double x = kPi * radius / cutoff_radius;
        const double sinc_value = (std::abs(x) < kOriginTolerance) ? 1.0 : std::sin(x) / x;
        shape_function.values[index] = sinc_value;
    }

    const double moment = radial_moment_integral(shape_function);
    if (moment <= 0.0)
    {
        throw std::runtime_error("Failed to normalize shape_function");
    }

    const double normalization = 1.0 / (kSqrtFourPi * moment);
    for (double &value : shape_function.values)
    {
        value *= normalization;
    }

    return shape_function;
}

double coulomb_inner_product_spherical(const RadialFunction &left, const RadialFunction &right)
{
    if (left.size() != right.size())
    {
        throw std::runtime_error("Coulomb inner product requires matched radial grids");
    }

    const std::vector<double> left_weights = trapezoid_weights(left.radii);
    double integral = 0.0;

    for (std::size_t left_index = 0; left_index < left.size(); ++left_index)
    {
        for (std::size_t right_index = 0; right_index < right.size(); ++right_index)
        {
            const double max_radius = std::max(left.radii[left_index], right.radii[right_index]);
            if (max_radius == 0.0)
            {
                continue;
            }

            const double kernel = 1.0 / max_radius;
            integral += left_weights[left_index] * left_weights[right_index] * left.radii[left_index] *
                        left.radii[left_index] * right.radii[right_index] * right.radii[right_index] *
                        left.values[left_index] * right.values[right_index] * kernel;
        }
    }

    return integral;
}

RadialFunction multiply_radial_functions(const RadialFunction &left, const RadialFunction &right)
{
    if (left.size() != right.size())
    {
        throw std::runtime_error("multiply_radial_functions requires matched radial grids");
    }

    RadialFunction product;
    product.radii = left.radii;
    product.values.resize(left.size(), 0.0);

    for (std::size_t index = 0; index < left.size(); ++index)
    {
        product.values[index] = left.values[index] * right.values[index];
    }

    return product;
}

RadialFunction subtract_radial_functions(const RadialFunction &left, const RadialFunction &right)
{
    if (left.size() != right.size())
    {
        throw std::runtime_error("subtract_radial_functions requires matched radial grids");
    }

    RadialFunction difference;
    difference.radii = left.radii;
    difference.values.resize(left.size(), 0.0);

    for (std::size_t index = 0; index < left.size(); ++index)
    {
        difference.values[index] = left.values[index] - right.values[index];
    }

    return difference;
}

double integral_phi_phi_over_r(const RadialFunction &left, const RadialFunction &right)
{
    const std::vector<double> weights = trapezoid_weights(left.radii);

    double integral = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index)
    {
        integral += weights[index] * left.radii[index] * left.values[index] * right.values[index];
    }

    return integral;
}

} // namespace

DenseMatrix build_two_index_coulomb_correction(const PAWSetup &setup)
{
    const auto shape_function_iterator = setup.metadata_blocks().find("shape_function");
    if (shape_function_iterator == setup.metadata_blocks().end())
    {
        throw std::runtime_error("PAW XML is missing shape_function metadata");
    }

    const RadialFunction &ae_core_density = setup.named_radial_functions().at("ae_core_density");
    const RadialFunction &pseudo_core_density = setup.named_radial_functions().at("pseudo_core_density");
    const RadialFunction shape_function =
        make_normalized_shape_function_00(shape_function_iterator->second, ae_core_density.radii);
    const double delta_a =
        radial_moment_integral(subtract_radial_functions(ae_core_density, pseudo_core_density)) -
        setup.atomic_number() / kSqrtFourPi;
    const double nc_tilde_g00 = coulomb_inner_product_spherical(pseudo_core_density, shape_function);
    const double g00_self = coulomb_inner_product_spherical(shape_function, shape_function);

    DenseMatrix correction = make_dense_matrix(setup.states().size(), setup.states().size());
    for (std::size_t left_state_index = 0; left_state_index < setup.states().size(); ++left_state_index)
    {
        for (std::size_t right_state_index = 0; right_state_index < setup.states().size(); ++right_state_index)
        {
            if (setup.states()[left_state_index].l != setup.states()[right_state_index].l)
            {
                continue;
            }

            const RadialFunction ae_pair =
                multiply_radial_functions(setup.all_electron_partial_waves_by_state()[left_state_index],
                                          setup.all_electron_partial_waves_by_state()[right_state_index]);
            const RadialFunction pseudo_pair =
                multiply_radial_functions(setup.pseudo_partial_waves_by_state()[left_state_index],
                                          setup.pseudo_partial_waves_by_state()[right_state_index]);
            const RadialFunction delta_pair = subtract_radial_functions(ae_pair, pseudo_pair);

            const double delta_00_ij = radial_moment_integral(delta_pair) / kSqrtFourPi;
            const double ae_core_term = coulomb_inner_product_spherical(ae_pair, ae_core_density);
            const double pseudo_core_term = coulomb_inner_product_spherical(pseudo_pair, pseudo_core_density);
            const double nuclear_term =
                setup.atomic_number() *
                integral_phi_phi_over_r(setup.all_electron_partial_waves_by_state()[left_state_index],
                                        setup.all_electron_partial_waves_by_state()[right_state_index]);
            const double shape_term = delta_a * coulomb_inner_product_spherical(pseudo_pair, shape_function);
            const double delta_term = delta_00_ij * (delta_a * nc_tilde_g00 + g00_self);

            correction[dense_matrix_extent(left_state_index), dense_matrix_extent(right_state_index)] =
                ae_core_term - pseudo_core_term - nuclear_term - shape_term - delta_term;
        }
    }

    return correction;
}

} // namespace dft
