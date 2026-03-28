#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace dft
{

using DenseMatrix = std::vector<std::vector<double>>;

struct XMLAttributeMap
{
    std::unordered_map<std::string, std::string> values;

    bool Has(const std::string &key) const;
    const std::string &GetString(const std::string &key) const;
    double GetDouble(const std::string &key) const;
};

struct RadialFunction
{
    std::vector<double> r;
    std::vector<double> values;

    bool empty() const { return r.empty() || values.empty(); }
    std::size_t size() const { return values.size(); }

    void Validate(const std::string &name) const;
};

struct PAWState
{
    std::string id;
    int n{-1};
    int l{0};
    bool has_n{false};
    double occupation{0.0};
    double rc{0.0};
    double energy{0.0};

    XMLAttributeMap attributes;
};

struct PAWChannel
{
    int l{0};
    std::vector<std::string> state_ids;

    std::vector<RadialFunction> projectors;
    std::vector<RadialFunction> all_electron_partial_waves;
    std::vector<RadialFunction> pseudo_partial_waves;

    DenseMatrix kinetic_correction;
    DenseMatrix static_coulomb_correction;
    DenseMatrix fixed_nonlocal_correction;
    std::vector<std::vector<double>> dij;

    void Validate(const std::string &setup_name) const;
    std::size_t NumProjectors() const { return projectors.size(); }
};

class PAWSetup
{
  public:
    PAWSetup() = default;
    PAWSetup(std::string symbol, double valence_charge, double cutoff_radius)
        : symbol_(std::move(symbol)), valence_charge_(valence_charge), cutoff_radius_(cutoff_radius)
    {
    }

    const std::string &symbol() const { return symbol_; }
    double valence_charge() const { return valence_charge_; }
    double cutoff_radius() const { return cutoff_radius_; }

    void set_symbol(std::string symbol) { symbol_ = std::move(symbol); }
    void set_valence_charge(double charge) { valence_charge_ = charge; }
    void set_cutoff_radius(double radius) { cutoff_radius_ = radius; }

    RadialFunction &local_potential() { return local_potential_; }
    const RadialFunction &local_potential() const { return local_potential_; }

    RadialFunction &core_density() { return core_density_; }
    const RadialFunction &core_density() const { return core_density_; }

    std::unordered_map<std::string, XMLAttributeMap> &metadata_blocks() { return metadata_blocks_; }
    const std::unordered_map<std::string, XMLAttributeMap> &metadata_blocks() const { return metadata_blocks_; }

    std::unordered_map<std::string, std::vector<double>> &radial_grids() { return radial_grids_; }
    const std::unordered_map<std::string, std::vector<double>> &radial_grids() const { return radial_grids_; }

    std::unordered_map<std::string, RadialFunction> &named_radial_functions() { return named_radial_functions_; }
    const std::unordered_map<std::string, RadialFunction> &named_radial_functions() const { return named_radial_functions_; }

    std::vector<PAWState> &states() { return states_; }
    const std::vector<PAWState> &states() const { return states_; }

    std::vector<RadialFunction> &all_electron_partial_waves_by_state() { return all_electron_partial_waves_by_state_; }
    const std::vector<RadialFunction> &all_electron_partial_waves_by_state() const { return all_electron_partial_waves_by_state_; }

    std::vector<RadialFunction> &pseudo_partial_waves_by_state() { return pseudo_partial_waves_by_state_; }
    const std::vector<RadialFunction> &pseudo_partial_waves_by_state() const { return pseudo_partial_waves_by_state_; }

    std::vector<RadialFunction> &projectors_by_state() { return projectors_by_state_; }
    const std::vector<RadialFunction> &projectors_by_state() const { return projectors_by_state_; }

    DenseMatrix &kinetic_energy_differences() { return kinetic_energy_differences_; }
    const DenseMatrix &kinetic_energy_differences() const { return kinetic_energy_differences_; }

    DenseMatrix &static_coulomb_correction() { return static_coulomb_correction_; }
    const DenseMatrix &static_coulomb_correction() const { return static_coulomb_correction_; }

    DenseMatrix &fixed_nonlocal_correction() { return fixed_nonlocal_correction_; }
    const DenseMatrix &fixed_nonlocal_correction() const { return fixed_nonlocal_correction_; }

    std::vector<PAWChannel> &channels() { return channels_; }
    const std::vector<PAWChannel> &channels() const { return channels_; }

    std::size_t NumChannels() const { return channels_.size(); }
    std::size_t NumProjectors() const;

    void Validate() const;

  private:
    std::string symbol_;
    double valence_charge_{0.0};
    double cutoff_radius_{0.0};

    RadialFunction local_potential_;
    RadialFunction core_density_;
    std::unordered_map<std::string, XMLAttributeMap> metadata_blocks_;
    std::unordered_map<std::string, std::vector<double>> radial_grids_;
    std::unordered_map<std::string, RadialFunction> named_radial_functions_;
    std::vector<PAWState> states_;
    std::vector<RadialFunction> all_electron_partial_waves_by_state_;
    std::vector<RadialFunction> pseudo_partial_waves_by_state_;
    std::vector<RadialFunction> projectors_by_state_;
    DenseMatrix kinetic_energy_differences_;
    DenseMatrix static_coulomb_correction_;
    DenseMatrix fixed_nonlocal_correction_;
    std::vector<PAWChannel> channels_;
};

class PAWSetupRegistry
{
  public:
    void Add(PAWSetup setup);
    bool Has(const std::string &symbol) const;
    const PAWSetup &Get(const std::string &symbol) const;

  private:
    std::unordered_map<std::string, PAWSetup> setups_;
};

PAWSetup LoadPAWSetupXML(const std::string &filename);

} // namespace dft
