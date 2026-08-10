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
 * @file InhomogeneousDensityFunction.hpp
 *
 * @brief Inhomogeneous DensityFunction implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef INHOMOGENEOUSDENSITYFUNCTION_HPP
#define INHOMOGENEOUSDENSITYFUNCTION_HPP

#include "DensityFunction.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include <cmath>
#include <algorithm>


class InhomogeneousDensityFunction : public DensityFunction {
private:
  const double _base_density;
  const double _temperature;
  const double _neutral_fraction_H;
  const double _dust_to_gas;
  const double _fraction_silicates;

  // Custom parameters for clump layout
  const int _seed;
  const double _frequency; // Restored: acts as the direct clump count across the box width
  const double _amplitude;

  // Cached box dimensions to prevent integer overflow and scope issues
  const double _anchor_x;
  const double _anchor_y;
  const double _anchor_z;
  const double _side_x;
  const double _side_y;
  const double _side_z;

  /**
   * @brief A completely deterministic 3D hash function.
   * Given integer coordinates and a seed, it returns a value between 0.0 and 1.0.
   */
  inline double hash3D(int x, int y, int z) const {
    long long n = x * 1619 + y * 31337 + z * 6971 + _seed * 1013;
    n = (n >> 13) ^ n;
    long long nn = (n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff;
    return static_cast<double>(nn) / 2147483647.0;
  }

  /**
   * @brief Smooth fade function (S-curve) to blend boundaries smoothly.
   */
  inline double fade(double t) const {
    return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
  }

  /**
   * @brief Manual 3D Value Noise generator.
   */
  double get_noise(double x, double y, double z) const {
    int ix = static_cast<int>(std::floor(x));
    int iy = static_cast<int>(std::floor(y));
    int iz = static_cast<int>(std::floor(z));

    double fx = x - std::floor(x);
    double fy = y - std::floor(y);
    double fz = z - std::floor(z);

    double u = fade(fx);
    double v = fade(fy);
    double w = fade(fz);

    double c000 = hash3D(ix,     iy,     iz);
    double c100 = hash3D(ix + 1, iy,     iz);
    double c010 = hash3D(ix,     iy + 1, iz);
    double c110 = hash3D(ix + 1, iy + 1, iz);
    double c001 = hash3D(ix,     iy,     iz + 1);
    double c101 = hash3D(ix + 1, iy,     iz + 1);
    double c011 = hash3D(ix,     iy + 1, iz + 1);
    double c111 = hash3D(ix + 1, iy + 1, iz + 1);

    double x00 = c000 + u * (c100 - c000);
    double x10 = c010 + u * (c110 - c010);
    double x01 = c001 + u * (c101 - c001);
    double x11 = c011 + u * (c111 - c011);

    double y0 = x00 + v * (x10 - x00);
    double y1 = x01 + v * (x11 - x01);

    return y0 + w * (y1 - y0);
  }

public:
  // Primary Constructor
  InhomogeneousDensityFunction(double density, double temperature, double neutral_fraction_H,
                               double dust_to_gas, double fraction_silicates, int seed,
                               double frequency, double amplitude,
                               double ax, double ay, double az, double sx, double sy, double sz,
                               Log *log = nullptr)
      : _base_density(density), _temperature(temperature),
        _neutral_fraction_H(neutral_fraction_H), _dust_to_gas(dust_to_gas),
        _fraction_silicates(fraction_silicates), _seed(seed),
        _frequency(frequency), _amplitude(amplitude),
        _anchor_x(ax), _anchor_y(ay), _anchor_z(az),
        _side_x(sx), _side_y(sy), _side_z(sz) {

    if (log) {
      log->write_status("Created InhomogeneousDensityFunction with adjustable frequency configuration.");
    }
  }

  // ParameterFile Constructor (Fixes array subscripts by mapping direct scalar vector extraction)
  InhomogeneousDensityFunction(ParameterFile &params, Log *log = nullptr)
      : InhomogeneousDensityFunction(
            params.get_physical_value< QUANTITY_NUMBER_DENSITY >("DensityFunction:density", "100. cm^-3"),
            params.get_physical_value< QUANTITY_TEMPERATURE >("DensityFunction:temperature", "8000. K"),
            params.get_value< double >("DensityFunction:neutral fraction H", 1.e-6),
            params.get_value< double >("DensityFunction:dust to gas", 0.0),
            params.get_value< double >("DensityFunction:fraction silicates", 0.6),
            params.get_value< int >("DensityFunction:seed", 42),
            params.get_value< double >("DensityFunction:frequency", 3.0), // Maps direct default of 3 clumps
            params.get_value< double >("DensityFunction:amplitude", 3.0),
            // Resolve bracket-subscript mismatch by passing native component parameters directly
            params.get_physical_vector< QUANTITY_LENGTH >("SimulationBox:anchor").x(),
            params.get_physical_vector< QUANTITY_LENGTH >("SimulationBox:anchor").y(),
            params.get_physical_vector< QUANTITY_LENGTH >("SimulationBox:anchor").z(),
            params.get_physical_vector< QUANTITY_LENGTH >("SimulationBox:sides").x(),
            params.get_physical_vector< QUANTITY_LENGTH >("SimulationBox:sides").y(),
            params.get_physical_vector< QUANTITY_LENGTH >("SimulationBox:sides").z(),
            log) {}

  virtual DensityValues operator()(const Cell &cell) {
    DensityValues values;

    // 1. Fetch cell coordinates from the physical midpoint vector
    double raw_x = cell.get_cell_midpoint().x();
    double raw_y = cell.get_cell_midpoint().y();
    double raw_z = cell.get_cell_midpoint().z();

    // 2. Normalize coordinates strictly to a safe [0.0, 1.0] range using stored members
    double x_norm = (raw_x - _anchor_x) / _side_x;
    double y_norm = (raw_y - _anchor_y) / _side_y;
    double z_norm = (raw_z - _anchor_z) / _side_z;

    // 3. Scale the 0-1 normalized space by frequency (acting as target clump count)
    double noise = get_noise(x_norm * _frequency, y_norm * _frequency, z_norm * _frequency); 
    double centered_noise = 2.0 * noise - 1.0; 

    // 4. Apply an exponential scaling factor for high-amplitude clumps and voids
    double extreme_multiplier = std::exp(_amplitude * centered_noise);
    double final_density = _base_density * extreme_multiplier;
    final_density = std::max(1e-10, final_density); // Safety floor

    // 5. Assign fields
    values.set_number_density(final_density);
    values.set_temperature(_temperature);
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


#endif // INHOMOGENEOUSDENSITYFUNCTION_HPP
