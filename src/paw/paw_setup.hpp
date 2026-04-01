#pragma once
#ifndef DFT_PAW_SETUP_HPP
#define DFT_PAW_SETUP_HPP
#include "ndarray/ndarray.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace dft
{

using DenseMatrix = NDArray<double, 2, int>;

struct XmlAttributeMap
{
    std::unordered_map<std::string, std::string> values;

    bool has(const std::string &key) const;
    const std::string &get_string(const std::string &key) const;
    double get_double(const std::string &key) const;
};

struct RadialFunction
{
    std::vector<double> radii;
    std::vector<double> values;

    bool empty() const { return radii.empty() || values.empty(); }
    std::size_t size() const { return values.size(); }

    void validate(const std::string &name) const;
};

struct PAWState
{
    std::string id;
    int n{-1};
    int l{0};
    bool has_n{false};
    double occupation{0.0};
    double cutoff_radius{0.0};
    double energy{0.0};

    XmlAttributeMap attributes;
};

struct PAWChannel
{
    int l{0};
    std::vector<std::string> state_ids;

    std::vector<RadialFunction> projectors;
    std::vector<RadialFunction> all_electron_partial_waves;
    std::vector<RadialFunction> pseudo_partial_waves;

    void validate(const std::string &setup_name) const;
    std::size_t num_projectors() const { return projectors.size(); }
};

class PAWSetup
{
  public:
    PAWSetup() = default;
    PAWSetup(std::string symbol, double valence_charge, double cutoff_radius)
        : m_symbol(std::move(symbol)), m_valence_charge(valence_charge), m_cutoff_radius(cutoff_radius)
    {
    }

    const std::string &symbol() const { return m_symbol; }
    double atomic_number() const { return m_atomic_number; }
    double valence_charge() const { return m_valence_charge; }
    double cutoff_radius() const { return m_cutoff_radius; }

    void set_symbol(std::string symbol) { m_symbol = std::move(symbol); }
    void set_atomic_number(double atomic_number) { m_atomic_number = atomic_number; }
    void set_valence_charge(double charge) { m_valence_charge = charge; }
    void set_cutoff_radius(double radius) { m_cutoff_radius = radius; }

    RadialFunction &local_potential() { return m_local_potential; }
    const RadialFunction &local_potential() const { return m_local_potential; }

    RadialFunction &core_density() { return m_core_density; }
    const RadialFunction &core_density() const { return m_core_density; }

    std::unordered_map<std::string, XmlAttributeMap> &metadata_blocks() { return m_metadata_blocks; }
    const std::unordered_map<std::string, XmlAttributeMap> &metadata_blocks() const { return m_metadata_blocks; }

    std::unordered_map<std::string, std::vector<double>> &radial_grids() { return m_radial_grids; }
    const std::unordered_map<std::string, std::vector<double>> &radial_grids() const { return m_radial_grids; }

    std::unordered_map<std::string, RadialFunction> &named_radial_functions() { return m_named_radial_functions; }
    const std::unordered_map<std::string, RadialFunction> &named_radial_functions() const
    {
        return m_named_radial_functions;
    }

    std::vector<PAWState> &states() { return m_states; }
    const std::vector<PAWState> &states() const { return m_states; }

    std::vector<RadialFunction> &all_electron_partial_waves_by_state() { return m_all_electron_partial_waves_by_state; }
    const std::vector<RadialFunction> &all_electron_partial_waves_by_state() const
    {
        return m_all_electron_partial_waves_by_state;
    }

    std::vector<RadialFunction> &pseudo_partial_waves_by_state() { return m_pseudo_partial_waves_by_state; }
    const std::vector<RadialFunction> &pseudo_partial_waves_by_state() const { return m_pseudo_partial_waves_by_state; }

    std::vector<RadialFunction> &projectors_by_state() { return m_projectors_by_state; }
    const std::vector<RadialFunction> &projectors_by_state() const { return m_projectors_by_state; }

    std::vector<double> &kinetic_difference_values() { return m_kinetic_difference_values; }
    const std::vector<double> &kinetic_difference_values() const { return m_kinetic_difference_values; }

    std::vector<PAWChannel> &channels() { return m_channels; }
    const std::vector<PAWChannel> &channels() const { return m_channels; }

    std::size_t num_channels() const { return m_channels.size(); }
    std::size_t num_projectors() const;

    void validate() const;

  private:
    std::string m_symbol;
    double m_atomic_number{0.0};
    double m_valence_charge{0.0};
    double m_cutoff_radius{0.0};

    RadialFunction m_local_potential;
    RadialFunction m_core_density;
    std::unordered_map<std::string, XmlAttributeMap> m_metadata_blocks;
    std::unordered_map<std::string, std::vector<double>> m_radial_grids;
    std::unordered_map<std::string, RadialFunction> m_named_radial_functions;
    std::vector<PAWState> m_states;
    std::vector<RadialFunction> m_all_electron_partial_waves_by_state;
    std::vector<RadialFunction> m_pseudo_partial_waves_by_state;
    std::vector<RadialFunction> m_projectors_by_state;
    std::vector<double> m_kinetic_difference_values;
    std::vector<PAWChannel> m_channels;
};

class PAWSetupRegistry
{
  public:
    void add(PAWSetup setup);
    bool has(const std::string &symbol) const;
    const PAWSetup &get(const std::string &symbol) const;

  private:
    std::unordered_map<std::string, PAWSetup> m_setups;
};

PAWSetup load_paw_setup_xml(const std::string &filename);

} // namespace dft
#endif
