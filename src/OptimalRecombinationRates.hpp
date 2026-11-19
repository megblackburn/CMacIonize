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
#ifndef OPTIMALRECOMBINATIONRATES_HPP
#define OPTIMALRECOMBINATIONRATES_HPP

#include <memory>
#include <string>

#include "ParameterFile.hpp"
#include "RecombinationRates.hpp"

class OptimalRecombinationRates : public RecombinationRates {
private:
  bool _apply_diffuse_field;
  std::unique_ptr<RecombinationRates> _metal_recombination_rates;

public:
  OptimalRecombinationRates(ParameterFile &params);

  double get_recombination_rate(const int_fast32_t ion,
                                const double temperature) const override;

  double get_hydrogen_ground_state_recombination_rate(
      const double temperature) const;

    // Wood, Mathis & Ercolano (2004), sections 3.3 and 7
  // equation (25)
  double get_helium_ground_state_recombination_rate(const double temperature) const {
    const double T4 = temperature * 1.e-4;
    const double alpha_1_He = 1.54e-13 * std::pow(T4, -0.486);
    return alpha_1_He;
  } 
};



#endif // OPTIMALRECOMBINATIONRATES_HPP
