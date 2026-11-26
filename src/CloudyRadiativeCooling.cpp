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
 * @file CloudyRadiativeCooling.cpp
 *
 * @brief CloudyRadiativeCooling implementation.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#include "CloudyRadiativeCooling.hpp"
#include "CloudyDataLocation.hpp"
#include "UnitConverter.hpp"
#include <fstream>
#include <sstream>
#include <string>

/**
 * @brief Constructor.
 */
CloudyRadiativeCooling::CloudyRadiativeCooling() {
  const double nH =
      UnitConverter::to_SI< QUANTITY_NUMBER_DENSITY >(1., "cm^-3");
  const double nH2inv = 1. / (nH * nH);

#ifdef HAS_HELIUM
  const std::string cool_name = "RadLoss_MetalsOnly_Cloudy.rates";
  const std::string emis_name = "Metals_Only_Ionising_Emissivity.rates";
  const std::string cdf_name  = "Metals_Only_Ionising_Spectrum_CDF.dat";
#else
  const std::string cool_name = "RadLoss_Metals_and_He_Only_Cloudy.rates";
  const std::string emis_name = "Metals_and_He_Ionising_Emissivity.rates";
  const std::string cdf_name  = "Metals_and_He_Ionising_Spectrum_CDF.dat";
#endif

  // 1. Load Cooling Rates
  {
    std::ifstream drfile(std::string(CLOUDYDATALOCATION) + cool_name);
    cmac_assert_message(drfile.is_open(), "Could not open %s", cool_name.c_str());
    std::string line;
    std::getline(drfile, line);
    std::getline(drfile, line);
    for (uint_fast32_t i = 0; i < CLOUDYRADIATIVECOOLING_NUMTEMPERATURE; ++i) {
      std::getline(drfile, line);
      std::istringstream linestream(line);
      linestream >> _temperatures[i] >> _cooling_rates[i];
      _cooling_rates[i] = UnitConverter::to_SI< QUANTITY_ENERGY_CHANGE_RATE >(
                              _cooling_rates[i], "erg cm^-3 s^-1") *
                          nH2inv;
    }
  }

  // 2. Load Ionising Photon Emissivity (cm^3 s^-1 -> m^3 s^-1)
  {
    std::ifstream emfile(std::string(CLOUDYDATALOCATION) + emis_name);
    cmac_assert_message(emfile.is_open(), "Could not open %s", emis_name.c_str());
    std::string line;
    std::getline(emfile, line);
    std::getline(emfile, line);
    for (uint_fast32_t i = 0; i < CLOUDYRADIATIVECOOLING_NUMTEMPERATURE; ++i) {
      double t_dummy = 0.0;
      std::getline(emfile, line);
      std::istringstream linestream(line);
      linestream >> t_dummy >> _ionising_emissivity[i];
      _ionising_emissivity[i] *= 1.e-6; // cm^3 s^-1 -> m^3 s^-1
    }
  }

  // 3. Load 2D Spectral CDF(T, nu)
  {
    std::ifstream cdffile(std::string(CLOUDYDATALOCATION) + cdf_name);
    cmac_assert_message(cdffile.is_open(), "Could not open %s", cdf_name.c_str());
    std::string line;
    // Skip 3 comment header lines
    std::getline(cdffile, line);
    std::getline(cdffile, line);
    std::getline(cdffile, line);

    // Row 0: dummy 0.0 followed by 250 frequency values in Hz
    {
      std::getline(cdffile, line);
      std::istringstream linestream(line);
      double dummy_zero = 0.0;
      linestream >> dummy_zero;
      for (uint_fast32_t j = 0; j < CLOUDYRADIATIVECOOLING_NUMFREQ; ++j) {
        linestream >> _cdf_frequencies[j];
      }
    }

    // Subsequent 351 rows: T followed by 250 CDF values in [0, 1]
    for (uint_fast32_t i = 0; i < CLOUDYRADIATIVECOOLING_NUMTEMPERATURE; ++i) {
      std::getline(cdffile, line);
      std::istringstream linestream(line);
      double t_dummy = 0.0;
      linestream >> t_dummy;
      for (uint_fast32_t j = 0; j < CLOUDYRADIATIVECOOLING_NUMFREQ; ++j) {
        linestream >> _cdf_values[i][j];
      }
    }
  }

  _min_logT = std::log(_temperatures[0]);
  _inverse_avg_dlogT =
      static_cast< double >(CLOUDYRADIATIVECOOLING_NUMTEMPERATURE - 1) /
      (std::log(_temperatures[CLOUDYRADIATIVECOOLING_NUMTEMPERATURE - 1]) -
       _min_logT);
}
