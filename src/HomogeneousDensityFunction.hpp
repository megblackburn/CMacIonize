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
 * @file HomogeneousDensityFunction.hpp
 *
 * @brief Homogeneous DensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HOMOGENEOUSDENSITYFUNCTION_HPP
#define HOMOGENEOUSDENSITYFUNCTION_HPP

#include "DensityFunction.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

/**
 * @brief DensityFunction that returns a constant value for all coordinates,
 * corresponding to a homogeneous density field.
 */
class HomogeneousDensityFunction : public DensityFunction {
private:
  double _density;
  double _temperature;
  double _neutral_fraction_H;
  double _dust_to_gas;
  double _fraction_silicates;
  bool _seed_turbulence;
  double _target_mach;
  double _gamma;
  double _box_L;

  struct WaveMode { 
      double kx, ky, kz; 
      double phase_x, phase_y, phase_z; 
      double amp; 
  };

  std::vector<WaveMode> _modes;
  
public:
  /**
   * @brief Constructor.
   */
   HomogeneousDensityFunction(double density = 1., double temperature = 8000.,
                             double neutral_fraction_H = 1.e-6,
                             double dust_to_gas = 0.0, double fraction_silicates = 0.6,
                             Log *log = nullptr)
      : _density(density), _temperature(temperature),
        _neutral_fraction_H(neutral_fraction_H),
        _dust_to_gas(dust_to_gas), _fraction_silicates(fraction_silicates),
        _seed_turbulence(false), _target_mach(0.0), _gamma(1.6666666666666667), _box_L(0.0) {

    if (log) {
      log->write_status(
          "Created HomogeneousDensityFunction with constant density ", _density,
          " m^-3 and constant temperature ", _temperature, " K. Dust gas ratio of ",
           _dust_to_gas, " and fraction silicates ", _fraction_silicates);
    }
  }

  /**
   * @brief ParameterFile constructor.
   */
  HomogeneousDensityFunction(ParameterFile &params, Log *log = nullptr)
      : HomogeneousDensityFunction(
            params.get_physical_value< QUANTITY_NUMBER_DENSITY >(
                "DensityFunction:density", "100. cm^-3"),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:temperature", "8000. K"),
            params.get_value< double >("DensityFunction:neutral fraction H",
                                       1.e-6),
            params.get_value<double>("DensityFunction:dust to gas",0.0),
            params.get_value<double>("DensityFunction:fraction silicates",0.6),
            log) {
    
    // Read optional configuration variables from the ParameterFile safely inside the constructor block
    _seed_turbulence = params.get_value<bool>("DensityFunction:seed_turbulence", false);
    
    if (_seed_turbulence) {
      _target_mach = params.get_value<double>("DensityFunction:target_mach", 5.0);
      _gamma = params.get_value<double>("DensityFunction:gamma", 5.0 / 3.0);
      
      double cell_dx = params.get_physical_value< QUANTITY_LENGTH >("DensityFunction:cell_dx", "1. pc");
      _box_L = params.get_value<double>("DensityFunction:box_cells", 32.0) * cell_dx;
      
      std::mt19937 generator(42); 
      std::uniform_real_distribution<double> dist(0.0, 2.0 * M_PI);
      std::uniform_real_distribution<double> amp_dist(0.5, 1.5);
      
      double running_rms_sq = 0.0;
      size_t sampled_modes = 0;
      
      for (int kx = -2; kx <= 2; ++kx) {
        for (int ky = -2; ky <= 2; ++ky) {
            for (int kz = -2; kz <= 2; ++kz) {
                double k_sq = kx*kx + ky*ky + kz*kz;
                if (k_sq >= 1.0 && k_sq <= 9.0) {
                    WaveMode mode;
                    mode.kx = (2.0 * M_PI * kx) / _box_L;
                    mode.ky = (2.0 * M_PI * ky) / _box_L;
                    mode.kz = (2.0 * M_PI * kz) / _box_L;
                    mode.phase_x = dist(generator);
                    mode.phase_y = dist(generator);
                    mode.phase_z = dist(generator);
                    mode.amp = amp_dist(generator) / k_sq; 
                    _modes.push_back(mode);
                    
                    running_rms_sq += (mode.amp * mode.amp * 1.5); 
                    sampled_modes++;
                }
            }
        }
      }
      
      double theoretical_sigma = std::sqrt(running_rms_sq / 3.0);
      double m_p = 1.672622e-27;
      double k_B = 1.380649e-23;
      double m_mean = 0.5 * m_p * (1.0 + _neutral_fraction_H);
      double c_s = std::sqrt((_gamma * k_B * _temperature) / m_mean);
      double target_sigma = _target_mach * c_s;
      double norm_factor = (theoretical_sigma > 0.0) ? (target_sigma / theoretical_sigma) : 1.0;
      
      for (auto &mode : _modes) {
          mode.amp *= norm_factor;
      }
    }
  }

  /**
   * @brief Function that gives the density for a given cell.
   */
  virtual DensityValues operator()(const Cell &cell) override {
    DensityValues values;
    values.set_number_density(_density);
    values.set_temperature(_temperature);
    values.set_ionic_fraction(ION_H_n, _neutral_fraction_H);
    values.set_dust_gas_ratio(_dust_to_gas);
    values.set_fraction_silicates(_fraction_silicates);
    
    if (_seed_turbulence) {
        CoordinateVector<> pos = cell.get_cell_midpoint();
        double vx = 0.0, vy = 0.0, vz = 0.0;
        
        for (const auto &mode : _modes) {
            vx += mode.amp * std::sin(mode.kx * pos.x() + mode.phase_x);
            vy += mode.amp * std::sin(mode.ky * pos.y() + mode.phase_y);
            vz += mode.amp * std::sin(mode.kz * pos.z() + mode.phase_z);
        }
        
        CoordinateVector<> seeded_vel(vx, vy, vz);
        values.set_velocity(seeded_vel); 
    }

#ifdef HAS_HELIUM
    values.set_ionic_fraction(ION_He_n, 0.999);
    values.set_ionic_fraction(ION_He_p1, 0.001);
#endif
#ifdef HAS_CARBON
    values.set_ionic_fraction(ION_C_p1, 0.999);
    values.set_ionic_fraction(ION_C_p2, 0.001);
#endif
#ifdef HAS_NITROGEN
    values.set_ionic_fraction(ION_N_n, 0.999);
    values.set_ionic_fraction(ION_N_p1, 0.001);
    values.set_ionic_fraction(ION_N_p2, 0.000);
#endif
#ifdef HAS_OXYGEN
    values.set_ionic_fraction(ION_O_n, 0.999);
    values.set_ionic_fraction(ION_O_p1, 0.001);
    values.set_ionic_fraction(ION_O_p2, 0.000);
    values.set_ionic_fraction(ION_O_p3, 0.000);
#endif
#ifdef HAS_NEON
    values.set_ionic_fraction(ION_Ne_n, 0.999);
    values.set_ionic_fraction(ION_Ne_p1, 0.001);
    values.set_ionic_fraction(ION_Ne_p2, 0.000);
    values.set_ionic_fraction(ION_Ne_p3, 0.000);
#endif
#ifdef HAS_SULPHUR
    values.set_ionic_fraction(ION_S_p1, 0.999);
    values.set_ionic_fraction(ION_S_p2, 0.001);
    values.set_ionic_fraction(ION_S_p3, 0.000);
#endif
    return values;
  }
};


#endif