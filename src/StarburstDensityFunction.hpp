/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file StarburstDensityFunction.hpp
 *
 * @brief Starburst density function.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#ifndef STARBURSTDENSITYFUNCTION_HPP
#define STARBURSTDENSITYFUNCTION_HPP

#include "CoordinateVector.hpp"
#include "DensityFunction.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"
#include "RandomGenerator.hpp"

#include <cmath>
#include <string>
#include <stdexcept>

/**
 * @brief Starburst Galaxy Density Function.
 * matches CGOLS for full M82 like simulation or can simulate a patch of a starburst galaxy
 */
class StarburstDensityFunction : public DensityFunction {
public:
  enum ProfileMode {
      MODE_GLOBAL = 0, // Uses fully dynamic 3D coordinate wells from center (0,0,0)
      MODE_PATCH  = 1  //  stratified local tall-box patch radius
  };

private:
  ProfileMode _mode;
  double _patch_radius_kpc;

  const double _disc_z;
  const double _b_inv;
  const double _exponent;
  const double _density_norm;

  const double _temperature;
  const bool _trace_initial_neutral_flag;
  const double _temperature_to_trace;
  const double _neutral_fraction;
  const double _velocity_dispersion;

  mutable RandomGenerator _random_generator;

  static inline double get_mean_particle_mass(const double neutral_fraction) {
    return 0.5 *
           PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PROTON_MASS) *
           (1. + neutral_fraction);
  }

  static inline double get_turbulent_gas_disc_scale_height(const double surface_density,
                                                          const double velocity_dispersion) {
    return (velocity_dispersion * velocity_dispersion / (M_PI *
            PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_NEWTON_CONSTANT) *
            surface_density));
  }

  static inline double get_mass_fraction_factor(const double exponent) {
    const double x = std::log10(-0.5 * exponent);
    const double x2 = x * x;
    const double y = 0.01499337 * x2 * x - 0.08454788 * x2 + 0.63503798 * x - 0.01018254;
    return std::pow(10., y);
  }

public:
  /**
   * @brief construct
   */
  inline StarburstDensityFunction(const ProfileMode mode,
                                  const double patch_radius_kpc,
                                  const double disc_z,
                                  const double surface_density,
                                  const double scale_height,
                                  const double gas_fraction,
                                  const double temperature,
                                  const double velocity_dispersion,
                                  const bool trace_initial_neutral_flag,
                                  const double temperature_to_trace,
                                  const double neutral_fraction)
      : _mode(mode),
        _patch_radius_kpc(patch_radius_kpc),
        _disc_z(disc_z), _b_inv(1. / scale_height),
        _exponent(-2. * scale_height / get_turbulent_gas_disc_scale_height(surface_density, velocity_dispersion)),
        _density_norm(0.5 * gas_fraction * surface_density *
                      get_mass_fraction_factor(_exponent) * _b_inv /
                      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PROTON_MASS)),
        _temperature(temperature),
        _trace_initial_neutral_flag(trace_initial_neutral_flag),
        _temperature_to_trace(temperature_to_trace),
        _neutral_fraction(neutral_fraction),
        _velocity_dispersion(velocity_dispersion),
        _random_generator() {}

  /**
   * @brief Read parameters
   */
  inline StarburstDensityFunction(ParameterFile &params)
      : _disc_z(params.get_physical_value< QUANTITY_LENGTH >("DensityFunction:disc z", "0. m")),
        _b_inv(1. / params.get_physical_value< QUANTITY_LENGTH >("DensityFunction:scale height", "200. pc")),
        _exponent(-2. * params.get_physical_value< QUANTITY_LENGTH >("DensityFunction:scale height", "200. pc") /
                  get_turbulent_gas_disc_scale_height(
                      params.get_physical_value< QUANTITY_SURFACE_DENSITY >("DensityFunction:surface density", "30. Msol pc^-2"),
                      params.get_physical_value< QUANTITY_VELOCITY >("InitialConditions:velocity_dispersion", "10. km s^-1"))),
        _density_norm(0.5 * params.get_value< double >("DensityFunction:gas fraction", 0.1) *
                      params.get_physical_value< QUANTITY_SURFACE_DENSITY >("DensityFunction:surface density", "30. Msol pc^-2") *
                      get_mass_fraction_factor(_exponent) * _b_inv /
                      PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PROTON_MASS)),
        _temperature(params.get_physical_value< QUANTITY_TEMPERATURE >("DensityFunction:temperature", "1.e4 K")),
        _trace_initial_neutral_flag(params.get_value< bool >("DensityFunction:trace initial neutral flag", false)),
        _temperature_to_trace(params.get_physical_value< QUANTITY_TEMPERATURE >("DensityFunction:temperature to trace", "500. K")),
        _neutral_fraction(params.get_value< double >("DensityFunction:neutral fraction", 1.0)),
        _velocity_dispersion(params.get_physical_value< QUANTITY_VELOCITY >("InitialConditions:velocity_dispersion", "10. km s^-1")),
        _random_generator() 
  {
      std::string mode_str = params.get_value<std::string>("ExternalPotential:mode", "patch");
      if (mode_str == "global") {
          _mode = MODE_GLOBAL;
      } else if (mode_str == "patch") {
          _mode = MODE_PATCH;
      } else {
          throw std::runtime_error("StarburstDensityFunction: Invalid parameter mode! Use 'global' or 'patch'.");
      }

      double r_m = params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:patch radius", "1.5 kpc");
      _patch_radius_kpc = r_m / 3.086e16; 
  }

  virtual ~StarburstDensityFunction() {}

  double get_cgols_potential_value(const double x, const double y, const double z) const {
    const double G = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_NEWTON_CONSTANT);
    
    const double M_disk  = 6.0e9  * 1.989e30; 
    const double a_disk  = 1100.  * 3.086e16; 
    const double b_disk  = 150.   * 3.086e16; 
    const double M_bulge = 5.0e8  * 1.989e30; 
    const double c_bulge = 180.   * 3.086e16; 
    const double M_halo  = 5.0e10 * 1.989e30; 
    const double r_s     = 9000.  * 3.086e16;

    const double R2 = x*x + y*y;
    const double r  = std::sqrt(R2 + z*z + 1e-20);

    const double z_param = std::sqrt(z*z + b_disk*b_disk);
    const double disk_term = a_disk + z_param;
    const double phi_disk = -G * M_disk / std::sqrt(R2 + disk_term * disk_term);

    const double phi_bulge = -G * M_bulge / (r + c_bulge);

    const double x_halo = r / r_s;
    const double nfw_norm = std::log(1.0 + (90000.0/9000.0)) - (10.0 / 11.0); 
    const double phi_halo = (-G * M_halo / (r * nfw_norm)) * std::log(1.0 + x_halo);

    return phi_disk + phi_bulge + phi_halo;
  }
  /**
   * @brief Set initial fields
   */
  virtual DensityValues operator()(const Cell &cell) override {
    double x = cell.get_cell_midpoint().x();
    double y = cell.get_cell_midpoint().y();
    const double z = cell.get_cell_midpoint().z() - _disc_z;

    if (_mode == MODE_PATCH) {
        const double R_fixed = _patch_radius_kpc * 3.086e16;
        x = R_fixed;
        y = 0.0;
    }

    // hydrostatic equilibrium
    const double phi_midplane = get_cgols_potential_value(x, y, 0.0);
    const double phi_local    = get_cgols_potential_value(x, y, z);
    const double delta_phi    = phi_local - phi_midplane;

    const double thermal_cs2 = (PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN) * _temperature) /
                               get_mean_particle_mass(_neutral_fraction);
    
    const double total_support_sigma2 = thermal_cs2 + (_velocity_dispersion * _velocity_dispersion);
    const double equilibrium_exponent = -delta_phi / total_support_sigma2;
    
    double number_density;
    if (equilibrium_exponent < -50.0) {
        number_density = 1.0e-10 * 1.e6; 
    } else {
        number_density = _density_norm * std::exp(equilibrium_exponent);
    }

    // Gaussian turbulent velocity vectors (Box-Muller)
    double u1 = _random_generator.get_uniform_random_double();
    double u2 = _random_generator.get_uniform_random_double();
    double u3 = _random_generator.get_uniform_random_double();
    double u4 = _random_generator.get_uniform_random_double();

    double g1 = std::sqrt(-2.0 * std::log(u1 + 1e-20)) * std::cos(2.0 * M_PI * u2);
    double g2 = std::sqrt(-2.0 * std::log(u1 + 1e-20)) * std::sin(2.0 * M_PI * u2);
    double g3 = std::sqrt(-2.0 * std::log(u3 + 1e-20)) * std::cos(2.0 * M_PI * u4);

    double vx = g1 * _velocity_dispersion;
    double vy = g2 * _velocity_dispersion;
    double vz = g3 * _velocity_dispersion;

    DensityValues values;
    values.set_number_density(number_density);
    values.set_temperature(_temperature);
    values.set_ionic_fraction(ION_H_n, _neutral_fraction);
    values.set_velocity(CoordinateVector<>(vx, vy, vz));

    if (_trace_initial_neutral_flag == true) {
        if (values.get_temperature() <= _temperature_to_trace) {
            values.set_initial_neutral_scalar_field(1.0);
            values.set_remaining_initial_neutral_scalar_field(1.0);
        } else {
            values.set_initial_neutral_scalar_field(0.0);
            values.set_remaining_initial_neutral_scalar_field(0.0);
        }
        values.set_cooled_neutral_scalar_field(0.0);
        values.set_remaining_cooled_neutral_scalar_field(0.0);
    } else {
        values.set_initial_neutral_scalar_field(0.0);
        values.set_cooled_neutral_scalar_field(0.0);
        values.set_remaining_initial_neutral_scalar_field(0.0);
        values.set_remaining_cooled_neutral_scalar_field(0.0);
    }

    return values;
  }
};

#endif // STARBURSTDENSITYFUNCTION_HPP
