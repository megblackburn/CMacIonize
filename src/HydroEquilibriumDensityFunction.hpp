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
 * @file HYDROEQUILIBRIUMDensityFunction.hpp
 *
 * @brief Disc patch density function.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#ifndef HYDROEQUILIBRIUMDENSITYFUNCTION_HPP
#define HYDROEQUILIBRIUMDENSITYFUNCTION_HPP

#include "CoordinateVector.hpp"
#include "DensityFunction.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"

#include <cmath>

/**
 * @brief Disc patch density function.
 * model for hydro equilibrium with the SILCC potential - used in Li, Bryan & Ostriker 2017, ApJ, 841, 101
 * not true equilibrium due to the DM potential and density floor
 */
class HydroEquilibriumDensityFunction : public DensityFunction {
private:
    const double _stellar_surface_density;
    const double _gas_surface_density;
    const double _midplane_number_density;
    const double _scale_height;
    const double _neutral_fraction;
    const double _density_floor;
    const double _temperature;
    const double _gamma;
    const bool _trace_initial_neutral_flag;
    const double _temperature_to_trace;

  static inline double get_mean_particle_mass(const double neutral_fraction) {
    return 0.5 *
           PhysicalConstants::get_physical_constant(
               PHYSICALCONSTANT_PROTON_MASS) *
           (1. + neutral_fraction);
  }

  static inline double sech(const double x) {
    return 1.0 / std::cosh(x);
  }

public:
  /**
   * @brief Constructor.
   *
   * @param disc_z Vertical position of the disc (in m).
   * @param surface_density Surface density of the disc, @f$\Sigma{}_M@f$
   * (in kg m^-2).
   * @param scale_height Scale height of the disc, @f$b_M@f$ (in m).
   * @param gas_fraction Fraction of the total mass content of the disc that is
   * in gas, @f$f_g@f$.
   * @param temperature Constant initial temperature, @f$T@f$ (in K).
   * @param neutral_fraction Constant initial neutral fraction for hydrogen,
   * @f$x_{\rm{}H}@f$.
   * @param observational_disc Use the observational vertical density profile?
   */
  inline HydroEquilibriumDensityFunction(const double stellar_surface_density,
                                  const double gas_surface_density,
                                  const double midplane_number_density,
                                  const double scale_height,
                                  const double neutral_fraction,
                                  const double density_floor,
                                  const double temperature,
                                  const double gamma,
                                  const bool trace_initial_neutral_flag,
                                  const double temperature_to_trace)
      : _stellar_surface_density(stellar_surface_density), _gas_surface_density(gas_surface_density), _midplane_number_density(midplane_number_density), _scale_height(scale_height), _neutral_fraction(neutral_fraction), _density_floor(density_floor), _temperature(temperature), _gamma(gamma), _trace_initial_neutral_flag(trace_initial_neutral_flag), _temperature_to_trace(temperature_to_trace) {}

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   */
  inline HydroEquilibriumDensityFunction(ParameterFile &params)
      : HydroEquilibriumDensityFunction(
            params.get_physical_value< QUANTITY_SURFACE_DENSITY >(
                "DensityFunction:stellar surface density", "10 Msol pc^-2"),
                params.get_physical_value< QUANTITY_SURFACE_DENSITY >("DensityFunction:gas surface density", "10. Msol pc^-2"),
                params.get_physical_value< QUANTITY_NUMBER_DENSITY >("DensityFunction:midplane number density", "0.822 cm^-3"),
                params.get_physical_value< QUANTITY_LENGTH >("DensityFunction:scale height", "300. pc"),
                params.get_value< double >("DensityFunction:neutral fraction", 0.99999),
                params.get_physical_value< QUANTITY_DENSITY >("DensityFunction:density floor", "3e-28 g cm^-3"),
                params.get_physical_value< QUANTITY_TEMPERATURE >("DensityFunction:temperature", "1.e4 K"),
                params.get_value< double >("Hydro:polytropic index", 5./3.),
                params.get_value< bool >("DensityFunction:trace initial neutral flag", false),
                params.get_physical_value< QUANTITY_TEMPERATURE >(
                    "DensityFunction:temperature to trace", "500. K")
            ) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~HydroEquilibriumDensityFunction() {} 
  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) {

    const double z = cell.get_cell_midpoint()[2];

    const double G = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_NEWTON_CONSTANT);
    const double velocity_dispersion = std::sqrt(M_PI * _scale_height * G * _stellar_surface_density);

    const double f_star = _stellar_surface_density / (_gas_surface_density + _stellar_surface_density);

    const double kB = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);
    const double sound_speed = std::sqrt(kB * _temperature / get_mean_particle_mass(_neutral_fraction));
    const double alpha = _gamma * velocity_dispersion * velocity_dispersion / ( f_star * sound_speed * sound_speed);


    double number_density = _midplane_number_density * std::pow((sech(z/_scale_height)), 2 * alpha);
    double number_density_floor = _density_floor / get_mean_particle_mass(_neutral_fraction);

    number_density = std::max(number_density, number_density_floor);
    

    DensityValues values;
    values.set_number_density(number_density);
    values.set_temperature(_temperature);
    values.set_ionic_fraction(ION_H_n, _neutral_fraction);

    if (_trace_initial_neutral_flag == true){
        if (values.get_temperature() <= _temperature_to_trace){
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

#endif // HYDROEQUILIBRIUMDENSITYFUNCTION_HPP
