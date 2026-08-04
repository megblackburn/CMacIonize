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
 * @file SodShockhDensityFunction.hpp
 *
 * @brief Disc patch density function.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef SODSHOCKDENSITYFUNCTIONMHD_HPP
#define SODSHOCKDENSITYFUNCTIONMHD_HPP

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
class SodShockDensityFunctionMHD : public DensityFunctionMHD {
private:
  /*! @brief Vertical position of the disc (in m). */
  const double _left_density;

  /*! @brief Inverse vertical scale height of the disc, @f$\frac{1}{b_M}@f$
   *  (in m^-1). */
  const double _right_density;

  /*! @brief Exponent of the density profile, @f$-\frac{2b_M}{b_g}@f$. */
  const double _left_temperature;

  /*! @brief Norm of the density profile, @f$\frac{\Sigma{}_g}{2b_Mm_p}@f$
   *  (in m^-3). */
  const double _right_temperature;

  /*! @brief Constant initial temperature, @f$T@f$ (in K). */
  const double _left_By;

  /*! @brief Constant initial temperature, @f$T@f$ (in K). */
  const double _right_By;

  /*! @brief Constant initial B_x, @f$T@f$ (in K). */
  const double _init_Bx;

  /*! @brief Constant initial B_y, @f$T@f$ (in K). */
  const double _init_By;

  /*! @brief Constant initial B_z, @f$T@f$ (in K). */
  const double _init_Bz;

  /*! @brief Constant initial B_scalar, @f$T@f$ (in K). */
  const double _init_B_scalar;

  /*! @brief Constant initial B_scalar, @f$T@f$ (in K). */
  const double _discontinuity;

  /*! @brief Constant initial neutral fraction for hydrogen,
   *  @f$x_{\rm{}H}@f$. */
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
  inline SodShockDensityFunctionMHD(const double left_density,
                                  const double right_density,
                                  const double left_temperature,
                                  const double right_temperature,
                                  const double left_By,
                                  const double right_By,
                                  const double init_Bx,
                                  const double init_By,
                                  const double init_Bz,
                                  const double init_B_scalar,
                                  const double discontinuity,
                                  const double neutral_fraction)
      : _left_density(left_density), _right_density(right_density), _left_temperature(left_temperature), _right_temperature(right_temperature),
        _left_By(left_By), _right_By(right_By), _init_Bx(init_Bx), _init_By(init_By), _init_Bz(init_Bz), _init_B_scalar(init_B_scalar), _discontinuity(discontinuity), _neutral_fraction(neutral_fraction) {}

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
  inline SodShockDensityFunctionMHD(ParameterFile &params)
      : SodShockDensityFunctionMHD(
            params.get_physical_value< QUANTITY_NUMBER_DENSITY >(
                "DensityFunction:left density", "1. cm^-3"),
            params.get_physical_value< QUANTITY_NUMBER_DENSITY >(
                "DensityFunction:right density", "1. cm^-3"),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:left temperature", "0. K"),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:right temperature", "0. K"),
            params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                "DensityFunction:left By", "0. T"),
            params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                "DensityFunction:right By", "0. T"),
            params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                "DensityFunction:init Bx", "0. T"),
            params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                "DensityFunction:init By", "0. T"),
            params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                "DensityFunction:init Bz", "0. T"),
            params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                "DensityFunction:init B_scalar", "0. T"),
            params.get_value< double >("DensityFunction:discontinuity", 0.5),
            params.get_value< double >("DensityFunction:neutral fraction",
                                       1e-6)) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~SodShockDensityFunctionMHD() {} /// Lewis's edited density function: mgb note 30.10.2025

  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) {

    double cell_x = cell.get_cell_midpoint().x();

    DensityValues values;

    double density = 0.;
    // double pressure = 0.;
    
    double temperature = 0.;
    double Bx = _init_Bx;
    double By = _init_By;
    double Bz = _init_Bz;

    if (cell_x < 0.5){
        density = _left_density;
        temperature = _left_temperature;
        By = _left_By;
    } else {
        density = _right_density;
        temperature = _right_temperature;
        By = _left_By;
    }
    
    values.set_number_density(density);
    values.set_temperature(temperature);
    values.set_ionic_fraction(ION_H_n, _neutral_fraction);
    values.set_magnetic_field(CoordinateVector<>(Bx,By,Bz));
    values.set_B_scalar(_init_B_scalar);

return values;
  }

};

#endif // SODSHOCKDENSITYFUNCTION_HPP
