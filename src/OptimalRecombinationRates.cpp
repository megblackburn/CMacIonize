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
 * @file OptimalRecombinationRates.cpp
 *
 * @brief RecombinationRates implementation with Verner's recombination rates:
 * implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "OptimalRecombinationRates.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#include "ChiantiRecombinationRates.hpp"
#include "Error.hpp"
#include "ParameterFile.hpp"
#include "VernerRecombinationRates.hpp"

/**
 * @brief Constructor.
 */
OptimalRecombinationRates::OptimalRecombinationRates(ParameterFile &params) {

  _apply_diffuse_field =
      params.get_value<bool>("DiffuseReemissionHandler:diffuse field", true) ||
      params.get_value<bool>("DiffuseReemissionHandler:time dependent diffuse field", true);

  std::string metal_type =
      params.get_value<std::string>("RecombinationRates:metal rate type", "Verner");

  std::cout << "Setting up OptimalRecombinationRates (H/HeII: Hui & Gnedin 1997, "
            << "HeI: Verner + Benjamin/Wood, Metals: " << metal_type << ")"
            << std::endl;

  std::cout << "Diffuse field reemission is "
            << (_apply_diffuse_field ? "enabled" : "disabled")
            << ", resorting to " << (_apply_diffuse_field ? "case A" : "case B")
            << std::endl;
  std::cout << "HAS " << NUMBER_OF_IONNAMES << " different ions." << std::endl;

  if (metal_type == "Verner") {
    _metal_recombination_rates =
        std::make_unique<VernerRecombinationRates>(params);
  } else if (metal_type == "Chianti") {
    _metal_recombination_rates =
        std::make_unique<ChiantiRecombinationRates>(params);
  } else {
    cmac_error("Unknown RecombinationRates:metal rate type: \"%s\"! Use \"Verner\" or \"Chianti\".",
               metal_type.c_str());
  }
}

double OptimalRecombinationRates::get_recombination_rate(
    const int_fast32_t ion, const double temperature) const {

  double rate = 0.;

  switch (ion) {

  case ION_H_n: {
    // Hui & Gnedin (1997, MNRAS, 292, 27, Appendix A), valid for 1 K < T < 10^9 K
    const double T_TR = 157807.0; // ionization threshold for H I in K
    const double lambda = 2.0 * T_TR / temperature;
    if (_apply_diffuse_field) {
      // Case A
      rate = 1.269e-13 * std::pow(lambda, 1.503) *
             std::pow(1.0 + std::pow(lambda / 0.522, 0.470), -1.923);
    } else {
      // Case B
      rate = 2.753e-14 * std::pow(lambda, 1.500) *
             std::pow(1.0 + std::pow(lambda / 2.740, 0.407), -2.242);
    }
    return 1.e-6 * rate; // convert cm^3 s^-1 to m^3 s^-1
  }

#ifdef HAS_HELIUM
  case ION_He_n: {
    // Verner & Ferland (1996) formula (4) with values from Table 1 (HeIa) — Case A
    // Valid in the range [3 K, 10^10 K]
    const double T1 = temperature / 4.266e-2;
    const double T2 = temperature / 4.677e6;
    rate = 9.356e-10 / (std::sqrt(T1) * std::pow(1. + std::sqrt(T1), 0.2108) *
                        std::pow(1. + std::sqrt(T2), 1.7892));
    if (!_apply_diffuse_field) {
      // Revert to Case B -> alpha_A - alpha_1_He using Wood, Mathis & Ercolano (2004)
      // / Benjamin et al. (1999) so it remains consistent with He I diffuse emission
      const double alpha_1_He = get_helium_ground_state_recombination_rate(temperature);
      rate = std::max(0.0, rate - alpha_1_He);
    }
    return 1.e-6 * rate; // convert cm^3 s^-1 to m^3 s^-1
  }

  case ION_He_p1: {
    // He II diffuse emission is NOT transported, so always use Case B
    // Hui & Gnedin (1997, MNRAS, 292, 27, Appendix A) Z = 2 scaling, valid [1 K, 10^9 K]
    const double T_TR = 631515.0; // ionization threshold for He II in K
    const double lambda = 2.0 * T_TR / temperature;
    rate = 2.0 * 2.753e-14 * std::pow(lambda, 1.500) *
           std::pow(1.0 + std::pow(lambda / 2.740, 0.407), -2.242);
    return 1.e-6 * rate; // convert cm^3 s^-1 to m^3 s^-1
  }
#endif

  default:
    // Delegate all heavy metals to Verner or CHIANTI 
    return _metal_recombination_rates->get_recombination_rate(ion, temperature);
  }
}
