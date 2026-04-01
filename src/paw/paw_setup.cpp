#include "paw/paw_setup.hpp"

#include "tinyxml2.h"

#include <sstream>
#include <unordered_map>

namespace dft
{

namespace
{

std::vector<double> parse_doubles(const std::string &text)
{
    std::vector<double> values;
    std::istringstream input_stream(text);

    double value = 0.0;
    while (input_stream >> value)
    {
        values.push_back(value);
    }

    return values;
}

const tinyxml2::XMLElement *require_child(const tinyxml2::XMLElement *parent, const char *name)
{
    const tinyxml2::XMLElement *child = parent->FirstChildElement(name);
    if (child == nullptr)
    {
        throw std::runtime_error(std::string("Missing XML element: ") + name);
    }

    return child;
}

std::string read_text(const tinyxml2::XMLElement *element, const char *name)
{
    const char *text = element->GetText();
    if (text == nullptr)
    {
        throw std::runtime_error(std::string("Missing text for XML element: ") + name);
    }

    return std::string(text);
}

int read_int_attribute(const tinyxml2::XMLElement *element, const char *attribute_name)
{
    int value = 0;
    const tinyxml2::XMLError error = element->QueryIntAttribute(attribute_name, &value);
    if (error != tinyxml2::XML_SUCCESS)
    {
        throw std::runtime_error(std::string("Failed to read integer XML attribute: ") + attribute_name);
    }

    return value;
}

double read_double_attribute(const tinyxml2::XMLElement *element, const char *attribute_name)
{
    double value = 0.0;
    const tinyxml2::XMLError error = element->QueryDoubleAttribute(attribute_name, &value);
    if (error != tinyxml2::XML_SUCCESS)
    {
        throw std::runtime_error(std::string("Failed to read double XML attribute: ") + attribute_name);
    }

    return value;
}

std::string read_string_attribute(const tinyxml2::XMLElement *element, const char *attribute_name)
{
    const char *value = element->Attribute(attribute_name);
    if (value == nullptr)
    {
        throw std::runtime_error(std::string("Missing XML attribute: ") + attribute_name);
    }

    return std::string(value);
}

XmlAttributeMap read_attributes(const tinyxml2::XMLElement *element)
{
    XmlAttributeMap attributes;
    for (const tinyxml2::XMLAttribute *attribute = element->FirstAttribute(); attribute != nullptr;
         attribute = attribute->Next())
    {
        attributes.values.emplace(attribute->Name(), attribute->Value());
    }

    return attributes;
}

RadialFunction load_radial_function(const tinyxml2::XMLElement *element, const std::vector<double> &radial_grid,
                                    const char *name)
{
    RadialFunction radial_function;
    radial_function.radii = radial_grid;
    radial_function.values = parse_doubles(read_text(element, name));
    radial_function.validate(name);

    if (radial_function.radii.size() != radial_function.values.size())
    {
        throw std::runtime_error(std::string(name) + ": radial grid and sampled values size mismatch");
    }

    return radial_function;
}

const std::vector<double> &resolve_grid(const tinyxml2::XMLElement *element,
                                        const std::unordered_map<std::string, std::vector<double>> &radial_grids)
{
    const std::string grid_id = read_string_attribute(element, "grid");
    const auto grid_iterator = radial_grids.find(grid_id);
    if (grid_iterator == radial_grids.end())
    {
        throw std::runtime_error("Unknown radial_grid id referenced by XML element: " + grid_id);
    }

    return grid_iterator->second;
}

} // namespace

bool XmlAttributeMap::has(const std::string &key) const
{
    return values.find(key) != values.end();
}

const std::string &XmlAttributeMap::get_string(const std::string &key) const
{
    const auto value_iterator = values.find(key);
    if (value_iterator == values.end())
    {
        throw std::runtime_error("Missing XML attribute key: " + key);
    }

    return value_iterator->second;
}

double XmlAttributeMap::get_double(const std::string &key) const
{
    return std::stod(get_string(key));
}

void RadialFunction::validate(const std::string &name) const
{
    if (radii.size() != values.size())
    {
        throw std::runtime_error(name + ": radial grid and values size mismatch");
    }

    if (radii.empty())
    {
        throw std::runtime_error(name + ": radial function is empty");
    }
}

void PAWChannel::validate(const std::string &setup_name) const
{
    if (projectors.empty())
    {
        throw std::runtime_error(setup_name + ": PAW channel has no projectors");
    }

    if (projectors.size() != all_electron_partial_waves.size() || projectors.size() != pseudo_partial_waves.size())
    {
        throw std::runtime_error(setup_name + ": inconsistent projector/partial-wave counts");
    }

    for (std::size_t projector_index = 0; projector_index < projectors.size(); ++projector_index)
    {
        projectors[projector_index].validate(setup_name + ": projector[" + std::to_string(projector_index) + "]");
        all_electron_partial_waves[projector_index].validate(setup_name + ": all-electron partial wave[" +
                                                             std::to_string(projector_index) + "]");
        pseudo_partial_waves[projector_index].validate(setup_name + ": pseudo partial wave[" +
                                                       std::to_string(projector_index) + "]");
    }
}

std::size_t PAWSetup::num_projectors() const
{
    std::size_t total_projectors = 0;
    for (const PAWChannel &channel : m_channels)
    {
        total_projectors += channel.num_projectors();
    }

    return total_projectors;
}

void PAWSetup::validate() const
{
    if (m_symbol.empty())
    {
        throw std::runtime_error("PAW setup symbol is empty");
    }

    if (m_atomic_number <= 0.0)
    {
        throw std::runtime_error(m_symbol + ": atomic number must be positive");
    }

    if (m_valence_charge <= 0.0)
    {
        throw std::runtime_error(m_symbol + ": valence charge must be positive");
    }

    if (m_cutoff_radius <= 0.0)
    {
        throw std::runtime_error(m_symbol + ": cutoff radius must be positive");
    }

    m_local_potential.validate(m_symbol + ": local potential");
    m_core_density.validate(m_symbol + ": core density");

    if (m_channels.empty())
    {
        throw std::runtime_error(m_symbol + ": no PAW channels defined");
    }

    if (m_states.empty())
    {
        throw std::runtime_error(m_symbol + ": no PAW valence states defined");
    }

    if (m_all_electron_partial_waves_by_state.size() != m_states.size() ||
        m_pseudo_partial_waves_by_state.size() != m_states.size() || m_projectors_by_state.size() != m_states.size())
    {
        throw std::runtime_error(m_symbol + ": state-resolved PAW radial data is incomplete");
    }

    if (m_kinetic_difference_values.size() != m_states.size() * m_states.size())
    {
        throw std::runtime_error(m_symbol + ": kinetic-difference values must match the number of states");
    }

    for (const PAWChannel &channel : m_channels)
    {
        channel.validate(m_symbol);
    }
}

void PAWSetupRegistry::add(PAWSetup setup)
{
    setup.validate();
    m_setups[setup.symbol()] = std::move(setup);
}

bool PAWSetupRegistry::has(const std::string &symbol) const
{
    return m_setups.find(symbol) != m_setups.end();
}

const PAWSetup &PAWSetupRegistry::get(const std::string &symbol) const
{
    const auto setup_iterator = m_setups.find(symbol);
    if (setup_iterator == m_setups.end())
    {
        throw std::runtime_error("No PAW setup registered for symbol: " + symbol);
    }

    return setup_iterator->second;
}

PAWSetup load_paw_setup_xml(const std::string &filename)
{
    tinyxml2::XMLDocument document;
    const tinyxml2::XMLError load_error = document.LoadFile(filename.c_str());
    if (load_error != tinyxml2::XML_SUCCESS)
    {
        throw std::runtime_error("Failed to load PAW XML file: " + filename);
    }

    const tinyxml2::XMLElement *dataset = document.RootElement();
    if (dataset == nullptr || std::string(dataset->Name()) != "paw_dataset")
    {
        throw std::runtime_error("Unsupported PAW XML root element");
    }

    const tinyxml2::XMLElement *atom = require_child(dataset, "atom");
    const tinyxml2::XMLElement *paw_radius = require_child(dataset, "paw_radius");
    const tinyxml2::XMLElement *pseudo_core_density = require_child(dataset, "pseudo_core_density");
    const tinyxml2::XMLElement *zero_potential = require_child(dataset, "zero_potential");
    const tinyxml2::XMLElement *valence_states = require_child(dataset, "valence_states");
    const tinyxml2::XMLElement *kinetic_energy_differences = require_child(dataset, "kinetic_energy_differences");

    PAWSetup setup(read_string_attribute(atom, "symbol"), read_double_attribute(atom, "valence"),
                   read_double_attribute(paw_radius, "rc"));
    setup.set_atomic_number(read_double_attribute(atom, "Z"));

    for (const tinyxml2::XMLElement *child = dataset->FirstChildElement(); child != nullptr;
         child = child->NextSiblingElement())
    {
        const std::string tag_name = child->Name();
        if (tag_name == "valence_states" || tag_name == "radial_grid" || tag_name == "ae_partial_wave" ||
            tag_name == "pseudo_partial_wave" || tag_name == "projector_function" || tag_name == "ae_core_density" ||
            tag_name == "pseudo_core_density" || tag_name == "ae_core_kinetic_energy_density" ||
            tag_name == "pseudo_core_kinetic_energy_density" || tag_name == "pseudo_valence_density" ||
            tag_name == "zero_potential" || tag_name == "blochl_local_ionic_potential")
        {
            continue;
        }

        setup.metadata_blocks()[tag_name] = read_attributes(child);
    }

    std::unordered_map<std::string, std::vector<double>> radial_grids;
    for (const tinyxml2::XMLElement *radial_grid = dataset->FirstChildElement("radial_grid"); radial_grid != nullptr;
         radial_grid = radial_grid->NextSiblingElement("radial_grid"))
    {
        const std::string grid_id = read_string_attribute(radial_grid, "id");
        const tinyxml2::XMLElement *radial_values = require_child(radial_grid, "values");
        std::vector<double> grid_values = parse_doubles(read_text(radial_values, "radial_grid/values"));
        if (grid_values.empty())
        {
            throw std::runtime_error("PAW XML radial grid is empty for id: " + grid_id);
        }

        radial_grids.emplace(grid_id, std::move(grid_values));
    }

    if (radial_grids.empty())
    {
        throw std::runtime_error("PAW XML does not contain any radial_grid definitions");
    }

    setup.radial_grids() = radial_grids;

    setup.core_density() =
        load_radial_function(pseudo_core_density, resolve_grid(pseudo_core_density, radial_grids), "pseudo_core_density");
    setup.local_potential() =
        load_radial_function(zero_potential, resolve_grid(zero_potential, radial_grids), "zero_potential");
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
                load_radial_function(element, resolve_grid(element, radial_grids), tag_name);
        }
    }

    std::unordered_map<std::string, std::size_t> state_to_channel_index;
    std::unordered_map<std::string, std::size_t> state_to_state_index;
    std::unordered_map<int, std::size_t> angular_momentum_to_channel_index;

    for (const tinyxml2::XMLElement *state = valence_states->FirstChildElement("state"); state != nullptr;
         state = state->NextSiblingElement("state"))
    {
        const std::string state_id = read_string_attribute(state, "id");
        const int angular_momentum = read_int_attribute(state, "l");

        PAWState paw_state;
        paw_state.id = state_id;
        paw_state.l = angular_momentum;
        paw_state.attributes = read_attributes(state);
        paw_state.has_n = paw_state.attributes.has("n");
        if (paw_state.has_n)
        {
            paw_state.n = std::stoi(paw_state.attributes.get_string("n"));
        }
        if (paw_state.attributes.has("f"))
        {
            paw_state.occupation = paw_state.attributes.get_double("f");
        }
        if (paw_state.attributes.has("rc"))
        {
            paw_state.cutoff_radius = paw_state.attributes.get_double("rc");
        }
        if (paw_state.attributes.has("e"))
        {
            paw_state.energy = paw_state.attributes.get_double("e");
        }

        setup.states().push_back(std::move(paw_state));
        state_to_state_index.emplace(state_id, setup.states().size() - 1);

        std::size_t channel_index = 0;
        const auto channel_iterator = angular_momentum_to_channel_index.find(angular_momentum);
        if (channel_iterator == angular_momentum_to_channel_index.end())
        {
            channel_index = setup.channels().size();
            PAWChannel channel;
            channel.l = angular_momentum;
            setup.channels().push_back(channel);
            angular_momentum_to_channel_index.emplace(angular_momentum, channel_index);
        }
        else
        {
            channel_index = channel_iterator->second;
        }

        state_to_channel_index.emplace(state_id, channel_index);
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
            const std::string state_id = read_string_attribute(element, "state");
            const auto channel_iterator = state_to_channel_index.find(state_id);
            if (channel_iterator == state_to_channel_index.end())
            {
                throw std::runtime_error("PAW XML references unknown state id: " + state_id);
            }

            const auto state_iterator = state_to_state_index.find(state_id);
            if (state_iterator == state_to_state_index.end())
            {
                throw std::runtime_error("PAW XML references unknown state index: " + state_id);
            }

            PAWChannel &channel = setup.channels()[channel_iterator->second];
            RadialFunction radial_function = load_radial_function(element, resolve_grid(element, radial_grids), tag_name);
            const std::string tag_value(tag_name);

            if (tag_value == "ae_partial_wave")
            {
                channel.all_electron_partial_waves.push_back(radial_function);
                setup.all_electron_partial_waves_by_state()[state_iterator->second] = std::move(radial_function);
            }
            else if (tag_value == "pseudo_partial_wave")
            {
                channel.pseudo_partial_waves.push_back(radial_function);
                setup.pseudo_partial_waves_by_state()[state_iterator->second] = std::move(radial_function);
            }
            else
            {
                channel.projectors.push_back(radial_function);
                setup.projectors_by_state()[state_iterator->second] = std::move(radial_function);
            }
        }
    }

    setup.kinetic_difference_values() =
        parse_doubles(read_text(kinetic_energy_differences, "kinetic_energy_differences"));

    setup.validate();
    return setup;
}

} // namespace dft
