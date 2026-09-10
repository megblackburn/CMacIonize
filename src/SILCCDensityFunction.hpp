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
 * @file SILCCDensityFunction.hpp
 *
 * @brief Disc patch density function.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SILCCDENSITYFUNCTION_HPP
#define SILCCDENSITYFUNCTION_HPP

#include "CoordinateVector.hpp"
#include "DensityFunction.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"

#include <cmath>

/**
 * @brief Disc patch density function.
 *
 * Represents a gas density profile that is initially in hydrostatic equilibrium
 * with a SILCCExternalPotential corresponding to a matter density profile
 * of the form
 * \f[
 *   \rho{}_M = \frac{\Sigma{}_M}{2b_M} \left(\cosh\left(\frac{z}{b_M}\right)
 *   \right)^{-2},
 * \f]
 * with @f$z@f$ the third component of the position, @f$\Sigma{}_M@f$ the
 * surface density of matter in the plane @f$z=0@f$ and @f$b_M@f$ a vertical
 * scale height for the density profile.
 *
 * The gas density profile itself has the general form
 * \f[
 *   \rho{}_g = \frac{\Sigma{}_g}{2b_M} \left(\cosh\left(\frac{z}{b_M}\right)
 *   \right)^{-\frac{2b_M}{b_g}},
 * \f]
 * where @f$\Sigma{}_g@f$ and @f$b_g@f$ are the surface density and scale height
 * for the gas density profile. The latter is given by
 * \f[
 *   b_g = \frac{k_B T}{\mu{} m_p \pi{} G \Sigma{}_M},
 * \f]
 * with @f$T@f$ the hydrostatic equilibrium temperature and @f$\mu{}@f$ the
 * mean molecular weight,
 * \f[
 *   \mu{} = \frac{1}{2} (1 + x_{\rm{}H}),
 * \f]
 * with @f$x_{\rm{}H}@f$ the hydrogen neutral fraction (we assume a hydrogen
 * only gas). @f$k_B@f$, @f$m_p@f$ and @f$G@f$ are respectively Bolzmann's
 * constant, the proton mass and Newton's constant.
 *
 * The gas surface density @f$\Sigma{}_g@f$ can be related to the total surface
 * density @f$\Sigma{}_M@f$ by imposing a fixed mass ratio
 * @f$f_g = \frac{M_g}{M_M}@f$, with
 * \f[
 *   M_X = \int_{-\infty{}}^{+\infty{}} \rho{}_X (z)~{\rm{}d}z, X = [M, g].
 * \f]
 * The corresponding expression is
 * \f[
 *   \Sigma{}_g = \frac{2}{I\left(-\frac{2b_M}{b_g}\right)} f_g \Sigma{}_M,
 * \f]
 * with
 * \f[
 *   I(d) = \int_{-\infty{}}^{+\infty{}} (\cosh(x))^d~{\rm{}d}x.
 * \f]
 * This integral has to be evaluated numerically. We use a third order
 * polynomial fit in log-log space to approximate it.
 */
class SILCCDensityFunction : public DensityFunction {
private:
    const double _surface_density;
    const double _scale_height;
    const double _neutral_fraction;
    const double _density_floor;
    const double _temperature;
    const bool _trace_initial_neutral_flag;
    const double _temperature_to_trace;

  static inline double get_mean_particle_mass(const double neutral_fraction) {
    return 0.5 *
           PhysicalConstants::get_physical_constant(
               PHYSICALCONSTANT_PROTON_MASS) *
           (1. + neutral_fraction);
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
  inline SILCCDensityFunction(const double surface_density,
                                  const double scale_height,
                                  const double neutral_fraction,
                                  const double density_floor,
                                  const double temperature,
                                  const bool trace_initial_neutral_flag,
                                  const double temperature_to_trace)
      : _surface_density(surface_density), _scale_height(scale_height), _neutral_fraction(neutral_fraction), _density_floor(density_floor), _temperature(temperature), _trace_initial_neutral_flag(trace_initial_neutral_flag), _temperature_to_trace(temperature_to_trace) {}

  /**
   * @brief ParameterFile constructor.
   *
   * We accept the following parameters:
   *  - disc z: Vertical position of the disc (default: 0. pc)
   *  - surface density: Surface density of the disc (default: 30. Msol pc^-2)
   *  - scale height: Scale height of the disc (default: 200. pc)
   *  - gas fraction: Fraction of the total mass content of the disc that is in
   *    gas (default: 0.1)
   *  - temperature: Constant initial temperature (default: 1.e4 K)
   *  - neutral fraction: Constant initial neutral fraction for hydrogen
   *    (default: 1.e-6)
   *  - observationaldisc: Use the observational vertical density profile. If
   *    false, use the hydrostatic-equilibrium profile (default: true)
   *
   * @param params ParameterFile to read from.
   */
  inline SILCCDensityFunction(ParameterFile &params)
      : SILCCDensityFunction(
            params.get_physical_value< QUANTITY_SURFACE_DENSITY >(
                "DensityFunction:surface density", "10 Msol pc^-2"),
                params.get_physical_value< QUANTITY_LENGTH >("DensityFunction:scale height", "30. pc"),
                params.get_value< double >("DensityFunction:neutral fraction", 0.9999),
                params.get_physical_value< QUANTITY_DENSITY >("DensityFunction:density floor", "1e-27 g cm^-3"),
                params.get_physical_value< QUANTITY_TEMPERATURE >("DensityFunction:temperature", "5000. K"),
                params.get_value< bool >("DensityFunction:trace initial neutral flag", false),
                params.get_physical_value< QUANTITY_TEMPERATURE >(
                    "DensityFunction:temperature to trace", "500. K")
            ) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~SILCCDensityFunction() {} 
  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) {

    const double z = cell.get_cell_midpoint()[2];

    const double midplane_density = _surface_density/ (std::sqrt(2. * M_PI) * _scale_height);
    double density = midplane_density * std::exp(-0.5 * std::pow(z / _scale_height, 2));

    density = std::max(density, _density_floor);
    
    const double number_density = density / get_mean_particle_mass(_neutral_fraction);

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

#endif // SILCCDENSITYFUNCTION_HPP
