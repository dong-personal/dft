#include "PAW.h"

#include "tinyxml2.h"

#include <cctype>
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

RadialFunction LoadRadialFunction_(const tinyxml2::XMLElement *element, const std::vector<double> &grid, const char *name)
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
    if (projectors.size() != all_electron_partial_waves.size() ||
        projectors.size() != pseudo_partial_waves.size())
    {
        throw std::runtime_error(setup_name + ": inconsistent projector/partial-wave counts");
    }
    if (dij.size() != projectors.size())
    {
        throw std::runtime_error(setup_name + ": dij row count must match projector count");
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

    (void)root;

    PAWSetup setup(ReadStringAttribute_(atom, "symbol"),
                   ReadDoubleAttribute_(atom, "valence"),
                   ReadDoubleAttribute_(paw_radius, "rc"));

    for (const tinyxml2::XMLElement *child = dataset->FirstChildElement(); child != nullptr;
         child = child->NextSiblingElement())
    {
        const std::string tag = child->Name();
        if (tag == "valence_states" || tag == "radial_grid" || tag == "ae_partial_wave" || tag == "pseudo_partial_wave" ||
            tag == "projector_function" || tag == "ae_core_density" || tag == "pseudo_core_density" ||
            tag == "ae_core_kinetic_energy_density" || tag == "pseudo_core_kinetic_energy_density" ||
            tag == "pseudo_valence_density" || tag == "zero_potential" || tag == "blochl_local_ionic_potential")
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

    setup.core_density() = LoadRadialFunction_(pseudo_core_density, ResolveGrid_(pseudo_core_density, grids), "pseudo_core_density");
    setup.local_potential() = LoadRadialFunction_(zero_potential, ResolveGrid_(zero_potential, grids), "zero_potential");
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
    }

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

            PAWChannel &channel = setup.channels()[channel_it->second];
            RadialFunction rf = LoadRadialFunction_(element, ResolveGrid_(element, grids), tag_name);

            const std::string tag(tag_name);
            if (tag == "ae_partial_wave")
            {
                channel.all_electron_partial_waves.push_back(std::move(rf));
            }
            else if (tag == "pseudo_partial_wave")
            {
                channel.pseudo_partial_waves.push_back(std::move(rf));
            }
            else
            {
                channel.projectors.push_back(std::move(rf));
            }
        }
    }

    for (auto &channel : setup.channels())
    {
        const std::size_t nproj = channel.projectors.size();
        channel.dij.assign(nproj, std::vector<double>(nproj, 0.0));
    }

    setup.Validate();
    return setup;
}

} // namespace dft
