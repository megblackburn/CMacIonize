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
 * @file KelvinHelmholtzDensityFunction.hpp
 *
 * @brief KelvinHelmholtz DensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef KELVINHELMHOLTZDENSITYFUNCTION_HPP
#define KELVINHELMHOLTZDENSITYFUNCTION_HPP

#include "DensityFunction.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include <cmath>

class KelvinHelmholtzDensityFunction : public DensityFunction {
private:
  const double _density_high;      // Dense fluid at the bottom half
  const double _density_low;       // Less dense fluid at the top half
  const double _temperature_cold;  // Temperature of the dense bottom fluid
  const double _neutral_fraction_H;
  const double _dust_to_gas;
  const double _fraction_silicates;
  const double _box_size_x;
  const double _box_size_y;
  const double _interface_width;

  // Integrated velocity properties
  const double _v_shear;
  const double _v_perturbation;

public:
  KelvinHelmholtzDensityFunction(double density_high = 2., double density_low = 1.,
                                 double temp_cold = 5000., double neutral_fraction_H = 1.e-6,
                                 double dust_to_gas = 0.0, double fraction_silicates = 0.6,
                                 double box_size_x = 10.0, double box_size_y = 10.0,
                                 double interface_width = 0.1, double v_shear = 1.0e4,
                                 double v_perturbation = 1.0e2, Log *log = nullptr)
      : _density_high(density_high), _density_low(density_low),
        _temperature_cold(temp_cold), _neutral_fraction_H(neutral_fraction_H),
        _dust_to_gas(dust_to_gas), _fraction_silicates(fraction_silicates),
        _box_size_x(box_size_x), _box_size_y(box_size_y),
        _interface_width(interface_width), _v_shear(v_shear),
        _v_perturbation(v_perturbation) {

    if (log) {
      log->write_status("Created KelvinHelmholtzDensityFunction with single stratified shear layer.");
    }
  }

  KelvinHelmholtzDensityFunction(ParameterFile &params, Log *log = nullptr)
      : KelvinHelmholtzDensityFunction(
            params.get_physical_value< QUANTITY_NUMBER_DENSITY >("DensityFunction:density high", "2 cm^-3"),
            params.get_physical_value< QUANTITY_NUMBER_DENSITY >("DensityFunction:density low", "1 cm^-3"),
            params.get_physical_value< QUANTITY_TEMPERATURE >("DensityFunction:temperature cold", "5000. K"),
            params.get_value< double >("DensityFunction:neutral fraction H", 1.e-06),
            params.get_value< double >("DensityFunction:dust to gas", 0.0),
            params.get_value< double >("DensityFunction:fraction silicates", 0.6),
            params.get_physical_vector< QUANTITY_LENGTH >("SimulationBox:sides").x(),
            params.get_physical_vector< QUANTITY_LENGTH >("SimulationBox:sides").y(),
            params.get_physical_value< QUANTITY_LENGTH >("DensityFunction:interface width", "0.1 pc"),
            params.get_physical_value< QUANTITY_VELOCITY >("DensityFunction:shear velocity", "10.0 km/s"),
            params.get_physical_value< QUANTITY_VELOCITY >("DensityFunction:perturbation amplitude", "0.1 km/s"),
            log) {}

  virtual DensityValues operator()(const Cell &cell) {
    DensityValues values;

    double x = cell.get_cell_midpoint().x();
    double y = cell.get_cell_midpoint().y();

    double centered_x = x + (0.5 * _box_size_x);
    double centered_y = y + (0.5 * _box_size_y);

    double profile = std::tanh((centered_y - 0.5 * _box_size_y) / _interface_width);

    double norm_fraction = 0.5 * (profile + 1.0);
    double final_density = _density_high + (_density_low - _density_high) * norm_fraction;

    double balanced_temperature = (_density_high * _temperature_cold) / final_density;

    double vx = _v_shear * profile;

    double k_mode = 2.0 * M_PI / _box_size_x;
    double vy = _v_perturbation * std::sin(k_mode * centered_x) * 
                std::exp(-std::pow((centered_y - 0.5 * _box_size_y) / _interface_width, 2));
    double vz = 0.0;

    // 7. Assign primitive values and the velocity vector block
    values.set_number_density(final_density);
    values.set_temperature(balanced_temperature);
    values.set_velocity(CoordinateVector<>(vx, vy, vz));
    
    values.set_ionic_fraction(ION_H_n, _neutral_fraction_H);
    values.set_dust_gas_ratio(_dust_to_gas);
    values.set_fraction_silicates(_fraction_silicates);

#ifdef HAS_HELIUM
    values.set_ionic_fraction(ION_He_n, 0.999);
    values.set_ionic_fraction(ION_He_p1, 0.001);
#endif

    return values;
  }
};

#endif // KELVINHELMHOLTZDENSITYFUNCTION_HPP