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

class StarburstDensityFunction : public DensityFunction {
public:
  enum ProfileMode {
      MODE_GLOBAL = 0,
      MODE_PATCH  = 1
  };

private:
  ProfileMode _mode;
  double _patch_radius_kpc;

  const double _disc_z;
  const double _density_norm;
  const double _temperature;
  const bool _trace_initial_neutral_flag;
  const double _temperature_to_trace;
  const double _neutral_fraction;
  const double _velocity_dispersion;

  // --- Dynamic Settable Parameters Parsed From ParameterFile ---
  double _M_disk_stars;
  double _a_disk_stars;
  double _b_disk_stars;
  double _M_halo;
  double _r_s_halo;
  double _c_nfw;
  double _R_disk_gas;
  double _n_peak_disk;
  double _n_halo_floor;
  double _T_halo_floor;

  mutable RandomGenerator _random_generator;

  static inline double get_mean_particle_mass(const double neutral_fraction) {
    return 0.5 * PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_PROTON_MASS) * (1. + neutral_fraction);
  }

public:
  /**
   * @brief ParameterFile configuration runtime parsing constructor.
   */
  inline StarburstDensityFunction(ParameterFile &params)
      : _disc_z(params.get_physical_value< QUANTITY_LENGTH >("DensityFunction:disc z", "0. m")),
        _density_norm(1.0), // Redundant placeholder for base structure mapping
        _temperature(params.get_physical_value< QUANTITY_TEMPERATURE >("DensityFunction:temperature", "1.e4 K")),
        _trace_initial_neutral_flag(params.get_value< bool >("DensityFunction:trace initial neutral flag", false)),
        _temperature_to_trace(params.get_physical_value< QUANTITY_TEMPERATURE >("DensityFunction:temperature to trace", "500. K")),
        _neutral_fraction(params.get_value< double >("DensityFunction:neutral fraction", 1.0)),
        _velocity_dispersion(params.get_physical_value< QUANTITY_VELOCITY >("InitialConditions:velocity_dispersion", "0. km s^-1")),
        _random_generator() 
  {
      // 1. Parse Simulation Runtime Operational Mode
      std::string mode_str = params.get_value<std::string>("ExternalPotential:mode", "global");
      if (mode_str == "global") _mode = MODE_GLOBAL;
      else if (mode_str == "patch") _mode = MODE_PATCH;
      else throw std::runtime_error("StarburstDensityFunction: Invalid mode! Use 'global' or 'patch'.");

      double r_m = params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:patch radius", "1.5 kpc");
      _patch_radius_kpc = r_m / 3.08568e19; 

      // 2. Fetch ALL Gravitational Potential Parameters directly from the Param File (No Hardcoding)
      _M_disk_stars = params.get_physical_value< QUANTITY_MASS >("ExternalPotential:disk mass", "1.0e10 Msol");
      _a_disk_stars = params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:disk scale length", "800. pc");
      _b_disk_stars = params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:disk scale height", "150. pc");
      
      _M_halo       = params.get_physical_value< QUANTITY_MASS >("ExternalPotential:halo mass", "5.0e10 Msol");
      _r_s_halo     = params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:halo scale radius", "5300. pc");
      _c_nfw        = params.get_value< double >("ExternalPotential:concentration", 10.0);

      // 3. Fetch Gas Disk & Hot Halo initial state parameters
      _R_disk_gas   = params.get_physical_value< QUANTITY_LENGTH >("InitialConditions:gas_scale_radius", "1.6 kpc");
      _n_peak_disk  = params.get_value< double >("InitialConditions:peak_midplane_density", 200.0) * 1.0e6; // cm^-3 to m^-3
      _n_halo_floor = params.get_value< double >("InitialConditions:halo_density_floor", 1.0e-3) * 1.0e6;   // cm^-3 to m^-3
      _T_halo_floor = params.get_physical_value< QUANTITY_TEMPERATURE >("InitialConditions:halo_temperature_floor", "2.0e6 K");
  }

  /**
   * @brief Computes CGOLS V potential dynamically using parameter class fields.
   */
  double get_cgols_potential_value(const double x, const double y, const double z) const {
    const double G = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_NEWTON_CONSTANT);
    
    const double R2 = x*x + y*y;
    const double r  = std::sqrt(R2 + z*z + 1e-20);

    // Miyamoto-Nagai Stellar Disk Profile
    const double z_param = std::sqrt(z*z + _b_disk_stars*_b_disk_stars);
    const double disk_term = _a_disk_stars + z_param;
    const double phi_disk = -G * _M_disk_stars / std::sqrt(R2 + disk_term * disk_term);

    // NFW Dark Matter Halo Profile
    const double x_halo = r / _r_s_halo;
    const double nfw_norm = std::log(1.0 + _c_nfw) - (_c_nfw / (1.0 + _c_nfw)); 
    const double phi_halo = (-G * _M_halo / (r * nfw_norm)) * std::log(1.0 + x_halo);

    return phi_disk + phi_halo;
  }

  /**
   * @brief Calculates rotation curve velocities derived from class fields.
   */
  double get_circular_velocity(const double R) const {
    if (R < 1.0) return 0.0;
    const double G = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_NEWTON_CONSTANT);
    
    const double disk_term = _a_disk_stars + _b_disk_stars;
    const double disk_accel = (G * _M_disk_stars * R) / std::pow(R*R + disk_term*disk_term, 1.5);

    const double nfw_norm = std::log(1.0 + _c_nfw) - (_c_nfw / (1.0 + _c_nfw));
    const double x_halo = R / _r_s_halo;
    const double halo_accel = (G * _M_halo / (R*R * nfw_norm)) * (std::log(1.0 + x_halo) - (x_halo / (1.0 + x_halo)));

    return std::sqrt(R * (disk_accel + halo_accel));
  }

  /**
   * @brief Set initial fields
   */
    /**
   * @brief Initialises cellular conditions.
   */
  virtual DensityValues operator()(const Cell &cell) override {
    double x = cell.get_cell_midpoint().x();
    double y = cell.get_cell_midpoint().y();
    const double z = cell.get_cell_midpoint().z() - _disc_z;

    double R = std::sqrt(x*x + y*y);

    if (_mode == MODE_PATCH) {
        const double R_fixed = _patch_radius_kpc * 3.08568e19; // Safe 1 kpc conversion scale
        x = R_fixed;
        y = 0.0;
        R = R_fixed;
    }

    const double k_B = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);

    // A. Evaluate Exponential Gas Disk Hydrostatic Equilibrium Component
    const double phi_midplane = get_cgols_potential_value(x, y, 0.0);
    const double phi_local    = get_cgols_potential_value(x, y, z);
    const double delta_phi    = phi_local - phi_midplane;

    const double thermal_cs2 = (k_B * _temperature) / get_mean_particle_mass(_neutral_fraction);
    const double scale_exponent = -delta_phi / thermal_cs2;

    double n_disk = 0.0;
    if (scale_exponent >= -50.0) {
        // Core tracking: Midplane normalization scales exponentially with radius R
        double n_midplane_R = _n_peak_disk * std::exp(-R / _R_disk_gas);
        n_disk = n_midplane_R * std::exp(scale_exponent);
    }

    double final_density = 0.0;
    double final_temperature = _temperature;

    // Superimpose gas disk inside the background hot adiabatic halo envelope
    if (n_disk >= _n_halo_floor) {
        final_density = n_disk;
        final_temperature = _temperature; // Cool gas disk layer (1.0e4 K)
    } else {
        final_density = _n_halo_floor;
        final_temperature = _T_halo_floor; // Hot Halo envelope (2.0e6 K)
    }

    // C. Isotropic Micro-Turbulent Velocity Dispersion Component (Box-Muller)
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

    // D. Compute Large-Scale Galactic Rotation Component in GLOBAL mode
    if (_mode == MODE_GLOBAL && R > 1.0 && final_temperature == _temperature) {
        double v_phi = get_circular_velocity(R); // Dynamically matches CGOLS centrifugal curve (~130 km/s)
        double pure_x = (std::abs(x) < 1.0e-5) ? 0.0 : x;
        double pure_y = (std::abs(y) < 1.0e-5) ? 0.0 : y;

        // Correct Cartesian projection for counter-clockwise rotation:
        // v_x = -v_phi * sin(theta) = -v_phi * (y / R)
        // v_y =  v_phi * cos(theta) =  v_phi * (x / R)
        vx += -v_phi * (pure_y / R);
        vy +=  v_phi * (pure_x / R);
    }

    DensityValues values;
    values.set_number_density(final_density);
    values.set_temperature(final_temperature); 
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
