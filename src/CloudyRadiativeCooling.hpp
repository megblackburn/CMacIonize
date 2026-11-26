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
 * @file CloudyRadiativeCooling.hpp
 *
 * @brief Radiative cooling rates from De Rijcke et al. (2013).
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#ifndef CLOUDYRADIATIVECOOLING_HPP
#define CLOUDYRADIATIVECOOLING_HPP

#include "Error.hpp"
#include "RadiativeCoolingHandler.hpp"

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <iostream>

/*! @brief Number of temperature bins. */
#define CLOUDYRADIATIVECOOLING_NUMTEMPERATURE 351

/*! @brief Number of frequency bins. */
#define CLOUDYRADIATIVECOOLING_NUMFREQ 250

class CloudyRadiativeCooling : public RadiativeCoolingHandler {
private:
  double _temperatures[CLOUDYRADIATIVECOOLING_NUMTEMPERATURE];
  double _cooling_rates[CLOUDYRADIATIVECOOLING_NUMTEMPERATURE];
  double _ionising_emissivity[CLOUDYRADIATIVECOOLING_NUMTEMPERATURE];

  double _cdf_frequencies[CLOUDYRADIATIVECOOLING_NUMFREQ];
  double _cdf_values[CLOUDYRADIATIVECOOLING_NUMTEMPERATURE][CLOUDYRADIATIVECOOLING_NUMFREQ];

  double _min_logT;
  double _inverse_avg_dlogT;

  inline void find_temperature_interval(const double temperature,
                                        uint_fast32_t &ilow,
                                        uint_fast32_t &ihigh,
                                        double &fac) const {
    if (temperature < _temperatures[0]) {
      ilow = 0;
      ihigh = 1;
    } else if (temperature >=
               _temperatures[CLOUDYRADIATIVECOOLING_NUMTEMPERATURE - 1]) {
      ilow = CLOUDYRADIATIVECOOLING_NUMTEMPERATURE - 2;
      ihigh = CLOUDYRADIATIVECOOLING_NUMTEMPERATURE - 1;
    } else {
      ilow = static_cast< uint_fast32_t >((std::log(temperature) - _min_logT) *
                                          _inverse_avg_dlogT);
      if (ilow >= CLOUDYRADIATIVECOOLING_NUMTEMPERATURE - 1) {
        ilow = CLOUDYRADIATIVECOOLING_NUMTEMPERATURE - 2;
      }
      if (temperature < _temperatures[ilow]) {
        ihigh = ilow;
        ilow = 0;
      } else {
        ihigh = CLOUDYRADIATIVECOOLING_NUMTEMPERATURE - 1;
      }
      while ((ihigh - ilow) != 1) {
        uint_fast32_t imid = (ilow + ihigh) >> 1;
        if (temperature >= _temperatures[imid]) {
          ilow = imid;
        } else {
          ihigh = imid;
        }
      }
    }
    fac = (temperature - _temperatures[ilow]) /
          (_temperatures[ihigh] - _temperatures[ilow]);
  }

public:
  CloudyRadiativeCooling();

  inline double get_cooling_rate(const double temperature) const override {
    uint_fast32_t ilow, ihigh;
    double fac;
    find_temperature_interval(temperature, ilow, ihigh, fac);

    const double cooling_rate =
        (1. - fac) * _cooling_rates[ilow] + fac * _cooling_rates[ihigh];
    cmac_assert(cooling_rate == cooling_rate);
    return std::max(cooling_rate, 0.);
  }

  /**
   * @brief Get the ionising photon emissivity (13.6 - 54.4 eV) for the given temperature.
   * @return Emissivity in m^3 s^-1 (multiply by n_H * n_e * V to get photons/s).
   */
  inline double get_ionizing_emissivity(const double temperature) const {
    if (temperature < _temperatures[0]) return 0.0;
    uint_fast32_t ilow, ihigh;
    double fac;
    find_temperature_interval(temperature, ilow, ihigh, fac);

    const double emis =
        (1. - fac) * _ionising_emissivity[ilow] + fac * _ionising_emissivity[ihigh];
    return std::max(emis, 0.0);
  }

  /**
   * @brief Sample a random ionising photon frequency (in Hz) from the 2D CDF(T, nu).
   */
  template < typename RNGType >
  inline double get_random_frequency(const double temperature,
                                     RNGType &local_rng) const {
    uint_fast32_t ilow, ihigh;
    double fac;
    find_temperature_interval(temperature, ilow, ihigh, fac);

    // Clamp temperature interpolation weight to [0, 1] for T > T_max
    const double w = std::min(1.0, std::max(0.0, fac));
    const double u = local_rng.get_uniform_random_double();

    // Binary search on the blended CDF row at temperature T
    uint_fast32_t flow = 0;
    uint_fast32_t fhigh = CLOUDYRADIATIVECOOLING_NUMFREQ - 1;
    while ((fhigh - flow) > 1) {
      uint_fast32_t fmid = (flow + fhigh) >> 1;
      const double cdf_mid =
          (1.0 - w) * _cdf_values[ilow][fmid] + w * _cdf_values[ihigh][fmid];
      if (u >= cdf_mid) {
        flow = fmid;
      } else {
        fhigh = fmid;
      }
    }

    const double cdf_low =
        (1.0 - w) * _cdf_values[ilow][flow] + w * _cdf_values[ihigh][flow];
    const double cdf_high =
        (1.0 - w) * _cdf_values[ilow][fhigh] + w * _cdf_values[ihigh][fhigh];
    const double denom = cdf_high - cdf_low;
    const double f_frac = (denom > 0.0) ? (u - cdf_low) / denom : 0.0;

    return (1.0 - f_frac) * _cdf_frequencies[flow] +
           f_frac * _cdf_frequencies[fhigh];
  }

  inline bool is_metal_only() const override { return true; }
  double get_minimum_temperature() const { return _temperatures[0]; }
  double get_maximum_temperature() const {
    return _temperatures[CLOUDYRADIATIVECOOLING_NUMTEMPERATURE - 1];
  }
};

#endif // CLOUDYRADIATIVECOOLING_HPP