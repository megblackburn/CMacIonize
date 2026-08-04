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
 * @file OrszagTanghDensityFunction.hpp
 *
 * @brief Disc patch density function.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef ORSZAGTANGDENSITYFUNCTIONMHD_HPP
#define ORSZAGTANGDENSITYFUNCTIONMHD_HPP

#include "CoordinateVector.hpp"
#include "DensityFunctionMHD.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"

#include <cmath>

/**
 * @brief Disc patch density function.
 *
 * for field loop mhd test case
 */
class OrszagTangDensityFunctionMHD : public DensityFunctionMHD {
private:
  /*! @brief Vertical position of the disc (in m). */
  const double _density;

  const double _delta_T;

  const double _temperature;

  const double _B_0;

  const double _v_0;

  const double _init_B_scalar;

  const double _L;

  const double _neutral_fraction;

  

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
   */
  inline OrszagTangDensityFunctionMHD(const double density,
                                    const double delta_T,
                                    const double temperature,
                                    const double B_0,
                                    const double v_0,
                                    const double init_B_scalar,
                                    const double L,
                                  const double neutral_fraction)
      : _density(density), _delta_T(delta_T), _temperature(temperature), _B_0(B_0), _v_0(v_0), _init_B_scalar(init_B_scalar), _L(L), _neutral_fraction(neutral_fraction) {}

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
   *
   * @param params ParameterFile to read from.
   */
  inline OrszagTangDensityFunctionMHD(ParameterFile &params)
      : OrszagTangDensityFunctionMHD(
            params.get_physical_value< QUANTITY_NUMBER_DENSITY >(
                "DensityFunction:density", "1.e19 m^-3"),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:delta T", "600. K"),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:temperature", "10000. K"),
            params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                "DensityFunction:B_0", "4.3e-5 T"),
            params.get_physical_value< QUANTITY_VELOCITY >(
                "DensityFunction:v_0", "15000 m s^-1"),
            params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                "DensityFunction:init B_scalar", "0. T"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "DensityFunction:L", "1.0 m"),
            params.get_value< double >("DensityFunction:neutral fraction",
                                       1e-6)) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~OrszagTangDensityFunctionMHD() {} /// Lewis's edited density function: mgb note 30.10.2025

  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) {

    double cell_x = cell.get_cell_midpoint().x();
    double cell_y = cell.get_cell_midpoint().y();

    DensityValues values;

    
    double density = _density;

    double Bx = _B_0 * std::sin(2. * M_PI * cell_y/_L);
    double By = _B_0 * std::sin(4. * M_PI * cell_x/_L);
    double Bz = 0.0;

    double vx = - _v_0 * std::sin(2. * M_PI * cell_y/_L);
    double vy = _v_0 * std::sin(4. * M_PI * cell_x/_L);
    double vz = 0.;

    double temperature = _temperature + _delta_T * (std::cos(4 * M_PI * cell_x/_L)-1.0);


    values.set_number_density(density);
    values.set_temperature(temperature);
    values.set_ionic_fraction(ION_H_n, _neutral_fraction);
    values.set_magnetic_field(CoordinateVector<>(Bx,By,Bz));
    values.set_B_scalar(_init_B_scalar);
    values.set_velocity(CoordinateVector<>(vx,vy,vz));

return values;
  }

};

#endif // ORSZAGTANGDENSITYFUNCTION_HPP
