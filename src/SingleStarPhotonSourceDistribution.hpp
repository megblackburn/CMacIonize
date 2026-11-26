/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CMacIonize is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file SingleStarPhotonSourceDistribution.hpp
 *
 * @brief PhotonSourceDistribution implementation containing a single stellar
 * source.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SINGLESTARPHOTONSOURCEDISTRIBUTION_HPP
#define SINGLESTARPHOTONSOURCEDISTRIBUTION_HPP

#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhotonSourceDistribution.hpp"
#include "RestartReader.hpp"
#include "RestartWriter.hpp"
#include "RecombinationRates.hpp"
#include "CrossSections.hpp"

/**
 * @brief General interface for photon source distribution functors.
 */
class SingleStarPhotonSourceDistribution : public PhotonSourceDistribution {
private:
  /*! @brief Position of the single stellar source (in m). */
  CoordinateVector<> _position;

  /*! @brief Luminosity of the single stellar source (in s^-1). */
  double _luminosity;
  double _initial_luminosity;

  // sources can move now
  double _velocity;

  bool _calculate_physical_luminosity_adjust;

  double _gamma;

  double _update_interval;

  double source_age = 0.0;

  double last_sf = 0.;

  double total_time = 0.;

  int_fast32_t iteration = 0;

  const double unit_Myr = 3.15576e13;


 // RecombinationRates *_recombination_rates;
 // CrossSections *_cross_sections;



  double get_mean_particle_mass(const double neutral_fraction) {
    return 0.5 *
          PhysicalConstants::get_physical_constant(
              PHYSICALCONSTANT_PROTON_MASS) *
          (1. + neutral_fraction);
  }


  double find_stromgren_radius(const RecombinationRates &recombination_rates, double cell_number_density, double source_ionizing_luminosity, double temperature, double source_age, double sound_speed){

    double alpha_B = recombination_rates.get_recombination_rate(ION_H_n, temperature);
    double recombination_time = 1.0/(alpha_B * cell_number_density);

    double R_ST = std::pow((3.0 * source_ionizing_luminosity)/(4.0 * M_PI * alpha_B * std::pow(cell_number_density, 2.0)), 1.0/3.0);


    // radiation phase of the stromgren volume
    double R_ST_rad_phase = R_ST * std::pow((1 - std::exp(-source_age/recombination_time)), 1.0/3.0);

    double transition_threshold = 3.5 * recombination_time; // time at which radius is 0.99 R_ST

    // D-type expansion phase of the stromgren volume
    if (source_age >= transition_threshold) {
      R_ST = R_ST * std::pow(1.0 + (7./4.) * ((sound_speed * (source_age - transition_threshold) /R_ST)), 4.0/7.0);
    } else {
      R_ST = R_ST_rad_phase;
    }

    return R_ST;

  }

  double get_luminosity_adjust(const RecombinationRates &recombination_rates, const CrossSections &cross_sections, const HydroVariables &state, double cell_number_density, double source_ionizing_luminosity, double temperature, double source_age, double cell_neutral_fraction, double cell_dx, double cell_volume) {

    double mean_molecular_mass = get_mean_particle_mass(cell_neutral_fraction);

    const double dvxdx = state.primitive_gradients(1).x();
    const double dvxdy = state.primitive_gradients(1).y();
    const double dvxdz = state.primitive_gradients(1).z();

    const double dvydx = state.primitive_gradients(2).x();
    const double dvydy = state.primitive_gradients(2).y();
    const double dvydz = state.primitive_gradients(2).z();

    const double dvzdx = state.primitive_gradients(3).x();
    const double dvzdy = state.primitive_gradients(3).y();
    const double dvzdz = state.primitive_gradients(3).z(); 

    double velocity_dispersion = cell_dx * std::sqrt(1.0/3.0 * 
                      (dvxdx*dvxdx + dvxdy*dvxdy + dvxdz*dvxdz + 
                        dvydx*dvydx + dvydy*dvydy + dvydz*dvydz +
                        dvzdx*dvzdx + dvzdy*dvzdy + dvzdz*dvzdz));  
                        

    double sound_speed = std::sqrt(PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN) * temperature * _gamma / mean_molecular_mass);

    double b = 0.35; // proportionality constant: Li et al. 2003, Federrath et al. 2008 - pure hydro sims 

    double Mach = velocity_dispersion/sound_speed;
    double turbulence_factor = std::sqrt(std::log(1.0 + b*b * Mach*Mach));

    double R_ST = find_stromgren_radius(recombination_rates,cell_number_density, source_ionizing_luminosity, temperature, source_age, sound_speed);

    double cross_section = cross_sections.get_cross_section(ION_H_n, 3.28847e15); // cross section for H at 13.6 eV

    double critical_density = 1 / (cross_section * R_ST);

    double density_threshold = std::log(critical_density/cell_number_density);

    double minimum_escape_frac;
    if (turbulence_factor != 0.0) {
      minimum_escape_frac = 0.5 * ( 1. + std::erf((density_threshold + (turbulence_factor*turbulence_factor/2.))/(std::sqrt(2.0) * turbulence_factor)));
    } else {
      minimum_escape_frac = 0.0;
    }
    double ionized_volume = 4./3. * M_PI * std::pow(R_ST, 3.0);
    double ionized_fraction = ionized_volume / cell_volume;       

    double escape_scale_factor;
    if (minimum_escape_frac <= 1.0) {
      escape_scale_factor = 1./(1. - minimum_escape_frac);
    } else {
       escape_scale_factor = 1.0;
    }
    
    double luminosity_adjust_factor = 1.0 - std::exp(-ionized_fraction * escape_scale_factor);

    return luminosity_adjust_factor;
  }


public:
  /**
   * @brief Constructor.
   *
   * @param position Position of the single stellar source (in m).
   * @param luminosity Luminosity of the single stellar source (in s^-1).
   * @param log Log to write logging information to.
   * @param calculate_physical_luminosity_adjust Whether to calculate the physical luminosity adjust.
   */
  SingleStarPhotonSourceDistribution(CoordinateVector<> position,
                                     double luminosity, double velocity, bool calculate_physical_luminosity_adjust, double gamma, double update_interval,
                                    Log *log = nullptr)
      : _position(position), _luminosity(luminosity), _velocity(velocity), _calculate_physical_luminosity_adjust(calculate_physical_luminosity_adjust), _gamma(gamma), _update_interval(update_interval) {

    if (log) {
      log->write_status(
          "Created SingleStarPhotonSourceDistribution at position [",
          _position.x(), " m, ", _position.y(), " m, ", _position.z(),
          " m], x velocity ", _velocity, " m/s with luminosity ", _luminosity, " s^-1.");
    }
    _initial_luminosity = _luminosity;
  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - position: Position of the stellar source (default: [0. pc, 0. pc, 0.
   *    pc])
   *  - luminosity: Ionizing luminosity of the source (default: 4.26e49 s^-1)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging information to.
   */
  SingleStarPhotonSourceDistribution(ParameterFile &params, Log *log = nullptr)
      : SingleStarPhotonSourceDistribution(
            params.get_physical_vector< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:position", "[0. pc, 0. pc, 0. pc]"),
            params.get_physical_value< QUANTITY_FREQUENCY >(
                "PhotonSourceDistribution:luminosity", "4.26e49 s^-1"),
            params.get_physical_value< QUANTITY_VELOCITY>(
                 "PhotonSourceDistribution:velocityx", "0.0 km s^-1"),
            params.get_value<bool>("PhotonSourceDistribution:calculate physical luminosity adjust", false),
            params.get_value<double>("Hydro:polytropic index", 5./3.),
            params.get_physical_value< QUANTITY_TIME >("PhotonSourceDistribution:update interval", "0.1 Myr"),
            log) {}

  /**
   * @brief Get the number of sources contained within this distribution.
   *
   * @return 1, as this distribution contains a single stellar source
   */
  virtual photonsourcenumber_t get_number_of_sources() const { return 1; }

  /**
   * @brief Get a valid position from the distribution.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return CoordinateVector of the single stellar source position (in m).
   */
  virtual CoordinateVector<> get_position(photonsourcenumber_t index) {
    return _position;
  };

  /**
   * @brief Get the weight of a photon source.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return Weight of the single photon source: 1.
   */
  virtual double get_weight(photonsourcenumber_t index) const { return 1.; }

  /**
   * @brief Get the luminosity of the single source.
   *
   * @return Luminosity (in s^-1).
   */
  virtual double get_total_luminosity() const { 
    return _luminosity;
  }


  virtual void float_sources(DensitySubGridCreator< HydroDensitySubGrid > *grid_creator, double timestep) {

     _position[0] += _velocity*timestep;


  }


  virtual bool update(DensitySubGridCreator< HydroDensitySubGrid > *grid_creator, const RecombinationRates *recombination_rates, const CrossSections *cross_sections, double actual_timestep) override {

    source_age += actual_timestep;
    total_time += actual_timestep;
    bool updated = false;
    
    std::cout<<"within update function" << std::endl;

    // update the luminosity adjust factor of all sources
    if (_calculate_physical_luminosity_adjust) {
      std::cout<< "calculating luminosity adjust factor at time " << total_time << std::endl;
      if (iteration == 0 || (total_time - last_sf >=  _update_interval || _update_interval < 0.)) {
        std::cout<< "Updating luminosity adjust factor for single star source at time " << total_time << std::endl;
        double luminosity_adjust_factor = 1.0;
        
            if (grid_creator->get_box().inside(_position)) {
            auto grid = grid_creator->get_subgrid(_position);
            auto cell = (*grid).get_hydro_cell(_position);
            const HydroVariables &state = cell.get_hydro_variables();

            double cell_num_dens = cell.get_ionization_variables().get_number_density();
            double cell_temp = cell.get_ionization_variables().get_temperature();
            double cell_neutral_frac = cell.get_ionization_variables().get_ionic_fraction(ION_H_n);

            double cell_volume = cell.get_volume();
            double cell_dx = std::pow(cell_volume, 1./3.);

            luminosity_adjust_factor = get_luminosity_adjust(*recombination_rates, *cross_sections, state, cell_num_dens, _initial_luminosity, cell_temp, source_age, cell_neutral_frac, cell_dx, cell_volume);
            _luminosity = _initial_luminosity * luminosity_adjust_factor;
          }
        
        iteration += 1;
        last_sf = total_time;
        updated = true;
        std::cout<< "Updated luminosity adjust factor for single star source to " << _luminosity << "with adjust factor " << luminosity_adjust_factor << std::endl;
      }
    }

    return updated;
  }

  /**
   * @brief Write the distribution to the given restart file.
   *
   * @param restart_writer RestartWriter to use.
   */
  virtual void write_restart_file(RestartWriter &restart_writer) const {

    _position.write_restart_file(restart_writer);
    restart_writer.write(_luminosity);
  }

  /**
   * @brief Restart constructor.
   *
   * @param restart_reader Restart file to read from.
   */
  inline SingleStarPhotonSourceDistribution(RestartReader &restart_reader)
      : _position(restart_reader),
        _luminosity(restart_reader.read< double >()) {}
};

#endif // PHOTONSOURCEDISTRIBUTION_HPP
