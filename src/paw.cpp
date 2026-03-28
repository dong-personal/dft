#include "paw.h"

#include "tinyxml2.h"

#include <cctype>
#include <cmath>
#include <sstream>
#include <unordered_map>

namespace dft
{

namespace
{

std::vector<double> ParseDoubles_(const std::string &text)
{
    std::vector<double> values;
    std::istringstream iss(text);
    double value = 0.0;
    while (iss >> value)
    {
        values.push_back(value);
    }
    return values;
}

const tinyxml2::XMLElement *RequireChild_(const tinyxml2::XMLElement *parent, const char *name)
{
    const tinyxml2::XMLElement *child = parent->FirstChildElement(name);
    if (child == nullptr)
    {
        throw std::runtime_error(std::string("Missing XML element: ") + name);
    }
    return child;
}

std::string ReadText_(const tinyxml2::XMLElement *element, const char *name)
{
    const char *text = element->GetText();
    if (text == nullptr)
    {
        throw std::runtime_error(std::string("Missing text for XML element: ") + name);
    }
    return std::string(text);
}

int ReadIntAttribute_(const tinyxml2::XMLElement *element, const char *attr_name)
{
    int value = 0;
    tinyxml2::XMLError err = element->QueryIntAttribute(attr_name, &value);
    if (err != tinyxml2::XML_SUCCESS)
    {
        throw std::runtime_error(std::string("Failed to read integer XML attribute: ") + attr_name);
    }
    return value;
}

double ReadDoubleAttribute_(const tinyxml2::XMLElement *element, const char *attr_name)
{
    double value = 0.0;
    tinyxml2::XMLError err = element->QueryDoubleAttribute(attr_name, &value);
    if (err != tinyxml2::XML_SUCCESS)
    {
        throw std::runtime_error(std::string("Failed to read double XML attribute: ") + attr_name);
    }
    return value;
}

std::string ReadStringAttribute_(const tinyxml2::XMLElement *element, const char *attr_name)
{
    const char *value = element->Attribute(attr_name);
    if (value == nullptr)
    {
        throw std::runtime_error(std::string("Missing XML attribute: ") + attr_name);
    }
    return std::string(value);
}

XMLAttributeMap ReadAttributes_(const tinyxml2::XMLElement *element)
{
    XMLAttributeMap attrs;
    for (const tinyxml2::XMLAttribute *attr = element->FirstAttribute(); attr != nullptr; attr = attr->Next())
    {
        attrs.values.emplace(attr->Name(), attr->Value());
    }
    return attrs;
}

RadialFunction LoadRadialFunction_(const tinyxml2::XMLElement *element, const std::vector<double> &grid,
                                   const char *name)
{
    RadialFunction rf;
    rf.r = grid;
    rf.values = ParseDoubles_(ReadText_(element, name));
    rf.Validate(name);
    if (rf.r.size() != rf.values.size())
    {
        throw std::runtime_error(std::string(name) + ": radial grid and sampled values size mismatch");
    }
    return rf;
}

const std::vector<double> &ResolveGrid_(const tinyxml2::XMLElement *element,
                                        const std::unordered_map<std::string, std::vector<double>> &grids)
{
    const std::string grid_id = ReadStringAttribute_(element, "grid");
    auto it = grids.find(grid_id);
    if (it == grids.end())
    {
        throw std::runtime_error("Unknown radial_grid id referenced by XML element: " + grid_id);
    }
    return it->second;
}

DenseMatrix ParseDenseSquareMatrix_(const tinyxml2::XMLElement *element, std::size_t size, const char *name)
{
    const std::vector<double> values = ParseDoubles_(ReadText_(element, name));
    if (values.size() != size * size)
    {
        throw std::runtime_error(std::string(name) + ": expected " + std::to_string(size * size) + " entries, got " +
                                 std::to_string(values.size()));
    }

    DenseMatrix matrix(size, std::vector<double>(size, 0.0));
    for (std::size_t i = 0; i < size; ++i)
    {
        for (std::size_t j = 0; j < size; ++j)
        {
            matrix[i][j] = values[i * size + j];
        }
    }
    return matrix;
}

double IntegrateRadialDifferenceWithPotential_(const RadialFunction &ae_i, const RadialFunction &ae_j,
                                               const RadialFunction &ps_i, const RadialFunction &ps_j,
                                               const RadialFunction &potential)
{
    const std::size_t n = ae_i.size();
    if (ae_j.size() != n || ps_i.size() != n || ps_j.size() != n || potential.size() != n)
    {
        throw std::runtime_error("Static Coulomb correction requires matched radial grids");
    }

    double integral = 0.0;
    for (std::size_t k = 0; k + 1 < n; ++k)
    {
        const double r0 = potential.r[k];
        const double r1 = potential.r[k + 1];
        const double dr = r1 - r0;

        const double f0 =
            r0 * r0 * potential.values[k] * (ae_i.values[k] * ae_j.values[k] - ps_i.values[k] * ps_j.values[k]);
        const double f1 = r1 * r1 * potential.values[k + 1] *
                          (ae_i.values[k + 1] * ae_j.values[k + 1] - ps_i.values[k + 1] * ps_j.values[k + 1]);

        integral += 0.5 * dr * (f0 + f1);
    }
    return integral;
}

std::vector<double> TrapezoidWeights_(const std::vector<double> &r)
{
    if (r.size() < 2)
    {
        throw std::runtime_error("Need at least two radial grid points for quadrature");
    }

    std::vector<double> weights(r.size(), 0.0);
    weights.front() = 0.5 * (r[1] - r[0]);
    for (std::size_t i = 1; i + 1 < r.size(); ++i)
    {
        weights[i] = 0.5 * (r[i + 1] - r[i - 1]);
    }
    weights.back() = 0.5 * (r.back() - r[r.size() - 2]);
    return weights;
}

double RadialMomentIntegral_(const RadialFunction &rf)
{
    const auto weights = TrapezoidWeights_(rf.r);
    double integral = 0.0;
    for (std::size_t i = 0; i < rf.size(); ++i)
    {
        integral += weights[i] * rf.r[i] * rf.r[i] * rf.values[i];
    }
    return integral;
}

// need to be fixed, normalized coef and specific from of the shape function
RadialFunction MakeNormalizedShapeFunction00_(const XMLAttributeMap &shape_attrs, const std::vector<double> &grid)
{
    if (!shape_attrs.Has("type") || !shape_attrs.Has("rc"))
    {
        throw std::runtime_error("shape_function is missing required attributes");
    }
    const std::string type = shape_attrs.GetString("type");
    if (type != "sinc")
    {
        throw std::runtime_error("Only sinc shape_function is currently implemented for static Coulomb correction");
    }

    const double rc = shape_attrs.GetDouble("rc");
    constexpr double sqrt4pi = 3.5449077018110318;

    RadialFunction g00;
    g00.r = grid;
    g00.values.resize(grid.size(), 0.0);
    for (std::size_t i = 0; i < grid.size(); ++i)
    {
        if (grid[i] > rc)
        {
            g00.values[i] = 0.0;
            continue;
        }
        const double x = M_PI * grid[i] / rc;
        const double sinc = (std::abs(x) < 1e-14) ? 1.0 : std::sin(x) / x;
        g00.values[i] = sinc;
    }

    const double moment = RadialMomentIntegral_(g00);
    if (moment <= 0.0)
    {
        throw std::runtime_error("Failed to normalize shape_function");
    }
    const double normalization = 1.0 / (sqrt4pi * moment);
    for (double &value : g00.values)
    {
        value *= normalization;
    }
    return g00;
}

double CoulombInnerProductSpherical_(const RadialFunction &lhs, const RadialFunction &rhs)
{
    if (lhs.size() != rhs.size())
    {
        throw std::runtime_error("Coulomb inner product requires matched radial grids");
    }

    const auto weights = TrapezoidWeights_(lhs.r);
    double integral = 0.0;
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        for (std::size_t j = 0; j < rhs.size(); ++j)
        {
            const double rmax = std::max(lhs.r[i], rhs.r[j]);
            if (rmax == 0.0)
            {
                continue;
            }
            const double kernel = 1.0 / rmax;
            integral += weights[i] * weights[j] * lhs.r[i] * lhs.r[i] * rhs.r[j] * rhs.r[j] * lhs.values[i] *
                        rhs.values[j] * kernel;
        }
    }
    return integral;
}

RadialFunction ProductRadialFunction_(const RadialFunction &lhs, const RadialFunction &rhs)
{
    if (lhs.size() != rhs.size())
    {
        throw std::runtime_error("ProductRadialFunction_ requires matched radial grids");
    }
    RadialFunction product;
    product.r = lhs.r;
    product.values.resize(lhs.size(), 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        product.values[i] = lhs.values[i] * rhs.values[i];
    }
    return product;
}

RadialFunction DifferenceRadialFunction_(const RadialFunction &lhs, const RadialFunction &rhs)
{
    if (lhs.size() != rhs.size())
    {
        throw std::runtime_error("DifferenceRadialFunction_ requires matched radial grids");
    }
    RadialFunction diff;
    diff.r = lhs.r;
    diff.values.resize(lhs.size(), 0.0);
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        diff.values[i] = lhs.values[i] - rhs.values[i];
    }
    return diff;
}

double IntegralPhiPhiOverR_(const RadialFunction &phi_i, const RadialFunction &phi_j)
{
    const auto weights = TrapezoidWeights_(phi_i.r);
    double integral = 0.0;
    for (std::size_t i = 0; i < phi_i.size(); ++i)
    {
        integral += weights[i] * phi_i.r[i] * phi_i.values[i] * phi_j.values[i];
    }
    return integral;
}

DenseMatrix AddDenseMatrices_(const DenseMatrix &lhs, const DenseMatrix &rhs, const std::string &name)
{
    if (lhs.size() != rhs.size())
    {
        throw std::runtime_error(name + ": matrix size mismatch");
    }

    DenseMatrix sum = lhs;
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        if (lhs[i].size() != rhs[i].size())
        {
            throw std::runtime_error(name + ": row size mismatch");
        }
        for (std::size_t j = 0; j < lhs[i].size(); ++j)
        {
            sum[i][j] += rhs[i][j];
        }
    }
    return sum;
}

DenseMatrix ExtractBlock_(const DenseMatrix &matrix, const std::vector<std::size_t> &indices)
{
    DenseMatrix block(indices.size(), std::vector<double>(indices.size(), 0.0));
    for (std::size_t i = 0; i < indices.size(); ++i)
    {
        for (std::size_t j = 0; j < indices.size(); ++j)
        {
            block[i][j] = matrix[indices[i]][indices[j]];
        }
    }
    return block;
}

} // namespace

bool XMLAttributeMap::Has(const std::string &key) const
{
    return values.find(key) != values.end();
}

const std::string &XMLAttributeMap::GetString(const std::string &key) const
{
    auto it = values.find(key);
    if (it == values.end())
    {
        throw std::runtime_error("Missing XML attribute key: " + key);
    }
    return it->second;
}

double XMLAttributeMap::GetDouble(const std::string &key) const
{
    return std::stod(GetString(key));
}

void RadialFunction::Validate(const std::string &name) const
{
    if (r.size() != values.size())
    {
        throw std::runtime_error(name + ": radial grid and values size mismatch");
    }
    if (r.empty())
    {
        throw std::runtime_error(name + ": radial function is empty");
    }
}

void PAWChannel::Validate(const std::string &setup_name) const
{
    if (projectors.empty())
    {
        throw std::runtime_error(setup_name + ": PAW channel has no projectors");
    }
    if (projectors.size() != all_electron_partial_waves.size() || projectors.size() != pseudo_partial_waves.size())
    {
        throw std::runtime_error(setup_name + ": inconsistent projector/partial-wave counts");
    }
    if (dij.size() != projectors.size())
    {
        throw std::runtime_error(setup_name + ": dij row count must match projector count");
    }
    if (kinetic_correction.size() != projectors.size() || static_coulomb_correction.size() != projectors.size() ||
        fixed_nonlocal_correction.size() != projectors.size())
    {
        throw std::runtime_error(setup_name + ": channel fixed nonlocal matrices must match projector count");
    }

    for (std::size_t i = 0; i < projectors.size(); ++i)
    {
        projectors[i].Validate(setup_name + ": projector[" + std::to_string(i) + "]");
        all_electron_partial_waves[i].Validate(setup_name + ": all-electron partial wave[" + std::to_string(i) + "]");
        pseudo_partial_waves[i].Validate(setup_name + ": pseudo partial wave[" + std::to_string(i) + "]");

        if (dij[i].size() != projectors.size())
        {
            throw std::runtime_error(setup_name + ": dij must be square");
        }
        if (kinetic_correction[i].size() != projectors.size() ||
            static_coulomb_correction[i].size() != projectors.size() ||
            fixed_nonlocal_correction[i].size() != projectors.size())
        {
            throw std::runtime_error(setup_name + ": channel fixed nonlocal matrices must be square");
        }
    }
}

std::size_t PAWSetup::NumProjectors() const
{
    std::size_t total = 0;
    for (const auto &channel : channels_)
    {
        total += channel.NumProjectors();
    }
    return total;
}

void PAWSetup::Validate() const
{
    if (symbol_.empty())
    {
        throw std::runtime_error("PAW setup symbol is empty");
    }
    if (valence_charge_ <= 0.0)
    {
        throw std::runtime_error(symbol_ + ": valence charge must be positive");
    }
    if (cutoff_radius_ <= 0.0)
    {
        throw std::runtime_error(symbol_ + ": cutoff radius must be positive");
    }

    local_potential_.Validate(symbol_ + ": local potential");
    core_density_.Validate(symbol_ + ": core density");

    if (channels_.empty())
    {
        throw std::runtime_error(symbol_ + ": no PAW channels defined");
    }
    if (states_.empty())
    {
        throw std::runtime_error(symbol_ + ": no PAW valence states defined");
    }
    if (all_electron_partial_waves_by_state_.size() != states_.size() ||
        pseudo_partial_waves_by_state_.size() != states_.size() || projectors_by_state_.size() != states_.size())
    {
        throw std::runtime_error(symbol_ + ": state-resolved PAW radial data is incomplete");
    }
    if (kinetic_energy_differences_.size() != states_.size() || static_coulomb_correction_.size() != states_.size() ||
        fixed_nonlocal_correction_.size() != states_.size())
    {
        throw std::runtime_error(symbol_ + ": fixed nonlocal matrices must match the number of states");
    }

    for (const auto &channel : channels_)
    {
        channel.Validate(symbol_);
    }
}

void PAWSetupRegistry::Add(PAWSetup setup)
{
    setup.Validate();
    setups_[setup.symbol()] = std::move(setup);
}

bool PAWSetupRegistry::Has(const std::string &symbol) const
{
    return setups_.find(symbol) != setups_.end();
}

const PAWSetup &PAWSetupRegistry::Get(const std::string &symbol) const
{
    auto it = setups_.find(symbol);
    if (it == setups_.end())
    {
        throw std::runtime_error("No PAW setup registered for symbol: " + symbol);
    }
    return it->second;
}

PAWSetup LoadPAWSetupXML(const std::string &filename)
{
    tinyxml2::XMLDocument doc;
    tinyxml2::XMLError load_err = doc.LoadFile(filename.c_str());
    if (load_err != tinyxml2::XML_SUCCESS)
    {
        throw std::runtime_error("Failed to load PAW XML file: " + filename);
    }

    const tinyxml2::XMLElement *root = RequireChild_(doc.RootElement(), "atom");
    const tinyxml2::XMLElement *dataset = doc.RootElement();
    if (dataset == nullptr || std::string(dataset->Name()) != "paw_dataset")
    {
        throw std::runtime_error("Unsupported PAW XML root element");
    }

    const tinyxml2::XMLElement *atom = RequireChild_(dataset, "atom");
    const tinyxml2::XMLElement *paw_radius = RequireChild_(dataset, "paw_radius");
    const tinyxml2::XMLElement *pseudo_core_density = RequireChild_(dataset, "pseudo_core_density");
    const tinyxml2::XMLElement *zero_potential = RequireChild_(dataset, "zero_potential");
    const tinyxml2::XMLElement *valence_states = RequireChild_(dataset, "valence_states");
    const tinyxml2::XMLElement *kinetic_energy_differences = RequireChild_(dataset, "kinetic_energy_differences");

    (void)root;

    PAWSetup setup(ReadStringAttribute_(atom, "symbol"), ReadDoubleAttribute_(atom, "valence"),
                   ReadDoubleAttribute_(paw_radius, "rc"));

    for (const tinyxml2::XMLElement *child = dataset->FirstChildElement(); child != nullptr;
         child = child->NextSiblingElement())
    {
        const std::string tag = child->Name();
        if (tag == "valence_states" || tag == "radial_grid" || tag == "ae_partial_wave" ||
            tag == "pseudo_partial_wave" || tag == "projector_function" || tag == "ae_core_density" ||
            tag == "pseudo_core_density" || tag == "ae_core_kinetic_energy_density" ||
            tag == "pseudo_core_kinetic_energy_density" || tag == "pseudo_valence_density" || tag == "zero_potential" ||
            tag == "blochl_local_ionic_potential")
        {
            continue;
        }
        setup.metadata_blocks()[tag] = ReadAttributes_(child);
    }

    std::unordered_map<std::string, std::vector<double>> grids;
    for (const tinyxml2::XMLElement *radial_grid = dataset->FirstChildElement("radial_grid"); radial_grid != nullptr;
         radial_grid = radial_grid->NextSiblingElement("radial_grid"))
    {
        const std::string grid_id = ReadStringAttribute_(radial_grid, "id");
        const tinyxml2::XMLElement *radial_values = RequireChild_(radial_grid, "values");
        std::vector<double> grid = ParseDoubles_(ReadText_(radial_values, "radial_grid/values"));
        if (grid.empty())
        {
            throw std::runtime_error("PAW XML radial grid is empty for id: " + grid_id);
        }
        grids.emplace(grid_id, std::move(grid));
    }

    if (grids.empty())
    {
        throw std::runtime_error("PAW XML does not contain any radial_grid definitions");
    }
    setup.radial_grids() = grids;

    setup.core_density() =
        LoadRadialFunction_(pseudo_core_density, ResolveGrid_(pseudo_core_density, grids), "pseudo_core_density");
    setup.local_potential() =
        LoadRadialFunction_(zero_potential, ResolveGrid_(zero_potential, grids), "zero_potential");
    setup.named_radial_functions()["pseudo_core_density"] = setup.core_density();
    setup.named_radial_functions()["zero_potential"] = setup.local_potential();

    for (const char *tag_name : {"ae_core_density", "pseudo_core_density", "ae_core_kinetic_energy_density",
                                 "pseudo_core_kinetic_energy_density", "pseudo_valence_density", "zero_potential",
                                 "blochl_local_ionic_potential"})
    {
        const tinyxml2::XMLElement *element = dataset->FirstChildElement(tag_name);
        if (element != nullptr)
        {
            setup.named_radial_functions()[tag_name] =
                LoadRadialFunction_(element, ResolveGrid_(element, grids), tag_name);
        }
    }

    std::unordered_map<std::string, std::size_t> state_to_channel;
    std::unordered_map<std::string, std::size_t> state_to_index;
    std::unordered_map<int, std::size_t> l_to_channel;

    for (const tinyxml2::XMLElement *state = valence_states->FirstChildElement("state"); state != nullptr;
         state = state->NextSiblingElement("state"))
    {
        const std::string state_id = ReadStringAttribute_(state, "id");
        const int l = ReadIntAttribute_(state, "l");

        PAWState paw_state;
        paw_state.id = state_id;
        paw_state.l = l;
        paw_state.attributes = ReadAttributes_(state);
        paw_state.has_n = paw_state.attributes.Has("n");
        if (paw_state.has_n)
        {
            paw_state.n = std::stoi(paw_state.attributes.GetString("n"));
        }
        if (paw_state.attributes.Has("f"))
        {
            paw_state.occupation = paw_state.attributes.GetDouble("f");
        }
        if (paw_state.attributes.Has("rc"))
        {
            paw_state.rc = paw_state.attributes.GetDouble("rc");
        }
        if (paw_state.attributes.Has("e"))
        {
            paw_state.energy = paw_state.attributes.GetDouble("e");
        }
        setup.states().push_back(std::move(paw_state));
        state_to_index.emplace(state_id, setup.states().size() - 1);

        std::size_t channel_index = 0;
        auto it = l_to_channel.find(l);
        if (it == l_to_channel.end())
        {
            channel_index = setup.channels().size();
            PAWChannel channel;
            channel.l = l;
            setup.channels().push_back(channel);
            l_to_channel.emplace(l, channel_index);
        }
        else
        {
            channel_index = it->second;
        }
        state_to_channel.emplace(state_id, channel_index);
        setup.channels()[channel_index].state_ids.push_back(state_id);
    }

    setup.all_electron_partial_waves_by_state().resize(setup.states().size());
    setup.pseudo_partial_waves_by_state().resize(setup.states().size());
    setup.projectors_by_state().resize(setup.states().size());

    for (const char *tag_name : {"ae_partial_wave", "pseudo_partial_wave", "projector_function"})
    {
        for (const tinyxml2::XMLElement *element = dataset->FirstChildElement(tag_name); element != nullptr;
             element = element->NextSiblingElement(tag_name))
        {
            const std::string state_id = ReadStringAttribute_(element, "state");
            auto channel_it = state_to_channel.find(state_id);
            if (channel_it == state_to_channel.end())
            {
                throw std::runtime_error("PAW XML references unknown state id: " + state_id);
            }
            auto state_it = state_to_index.find(state_id);
            if (state_it == state_to_index.end())
            {
                throw std::runtime_error("PAW XML references unknown state index: " + state_id);
            }

            PAWChannel &channel = setup.channels()[channel_it->second];
            RadialFunction rf = LoadRadialFunction_(element, ResolveGrid_(element, grids), tag_name);

            const std::string tag(tag_name);
            if (tag == "ae_partial_wave")
            {
                channel.all_electron_partial_waves.push_back(rf);
                setup.all_electron_partial_waves_by_state()[state_it->second] = std::move(rf);
            }
            else if (tag == "pseudo_partial_wave")
            {
                channel.pseudo_partial_waves.push_back(rf);
                setup.pseudo_partial_waves_by_state()[state_it->second] = std::move(rf);
            }
            else
            {
                channel.projectors.push_back(rf);
                setup.projectors_by_state()[state_it->second] = std::move(rf);
            }
        }
    }

    setup.kinetic_energy_differences() =
        ParseDenseSquareMatrix_(kinetic_energy_differences, setup.states().size(), "kinetic_energy_differences");

    const auto shape_it = setup.metadata_blocks().find("shape_function");
    if (shape_it == setup.metadata_blocks().end())
    {
        throw std::runtime_error("PAW XML is missing shape_function metadata");
    }
    const RadialFunction &ae_core_density = setup.named_radial_functions().at("ae_core_density");
    const RadialFunction &pseudo_core_density_rf = setup.named_radial_functions().at("pseudo_core_density");
    const RadialFunction g00 = MakeNormalizedShapeFunction00_(shape_it->second, ae_core_density.r);
    constexpr double sqrt4pi = 3.5449077018110318;
    const double delta_a = RadialMomentIntegral_(DifferenceRadialFunction_(ae_core_density, pseudo_core_density_rf)) -
                           ReadDoubleAttribute_(atom, "Z") / sqrt4pi;
    const double nc_tilde_g00 = CoulombInnerProductSpherical_(pseudo_core_density_rf, g00);
    const double g00_self = CoulombInnerProductSpherical_(g00, g00);

    setup.static_coulomb_correction().assign(setup.states().size(), std::vector<double>(setup.states().size(), 0.0));
    for (std::size_t i = 0; i < setup.states().size(); ++i)
    {
        for (std::size_t j = 0; j < setup.states().size(); ++j)
        {
            if (setup.states()[i].l != setup.states()[j].l)
            {
                continue;
            }
            const RadialFunction ae_pair = ProductRadialFunction_(setup.all_electron_partial_waves_by_state()[i],
                                                                  setup.all_electron_partial_waves_by_state()[j]);
            const RadialFunction ps_pair = ProductRadialFunction_(setup.pseudo_partial_waves_by_state()[i],
                                                                  setup.pseudo_partial_waves_by_state()[j]);
            const RadialFunction delta_pair = DifferenceRadialFunction_(ae_pair, ps_pair);

            const double delta00_ij = RadialMomentIntegral_(delta_pair) / sqrt4pi;
            const double term_ae_core = CoulombInnerProductSpherical_(ae_pair, ae_core_density);
            const double term_ps_core = CoulombInnerProductSpherical_(ps_pair, pseudo_core_density_rf);
            const double term_nuclear =
                ReadDoubleAttribute_(atom, "Z") * IntegralPhiPhiOverR_(setup.all_electron_partial_waves_by_state()[i],
                                                                       setup.all_electron_partial_waves_by_state()[j]);
            const double term_shape = delta_a * CoulombInnerProductSpherical_(ps_pair, g00);
            const double term_delta = delta00_ij * (delta_a * nc_tilde_g00 + g00_self);

            setup.static_coulomb_correction()[i][j] =
                term_ae_core - term_ps_core - term_nuclear - term_shape - term_delta;
        }
    }

    setup.fixed_nonlocal_correction() = AddDenseMatrices_(
        setup.kinetic_energy_differences(), setup.static_coulomb_correction(), "fixed_nonlocal_correction");

    for (auto &channel : setup.channels())
    {
        const std::size_t nproj = channel.projectors.size();
        channel.dij.assign(nproj, std::vector<double>(nproj, 0.0));
        std::vector<std::size_t> state_indices;
        state_indices.reserve(channel.state_ids.size());
        for (const std::string &state_id : channel.state_ids)
        {
            state_indices.push_back(state_to_index.at(state_id));
        }
        channel.kinetic_correction = ExtractBlock_(setup.kinetic_energy_differences(), state_indices);
        channel.static_coulomb_correction = ExtractBlock_(setup.static_coulomb_correction(), state_indices);
        channel.fixed_nonlocal_correction = ExtractBlock_(setup.fixed_nonlocal_correction(), state_indices);
    }

    setup.Validate();
    return setup;
}

} // namespace dft
