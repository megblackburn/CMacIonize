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
 * @file ChiantiRecombinationRates.hpp
 *
 * @brief RecombinationRates implementation with Chianti's recombination rates:
 * header.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef CHIANTIRECOMBINATIONRATES_HPP
#define CHIANTIRECOMBINATIONRATES_HPP

#include "RecombinationRates.hpp"
#include "ParameterFile.hpp"
#include <vector>
#include <cmath>
/**
 * @brief RecombinationRates implementation with Chianti's recombination rates.
 *

 */
class ChiantiRecombinationRates : public RecombinationRates {
private:

  /*! @brief flag to indicate the inclusion or exclusion of the diffuse field reemission*/
  const bool _apply_diffuse_field;
  /*! @brief Temperature values (in K). */
  std::vector<double> _temperatures;

  /*! @brief Collisional rate values (in m^-3 s^-1). */
  //double _collisional_rates[NUMBER_OF_IONNAMES][250];


  std::vector<std::vector<double>> _recomb_rates;
  /*! @brief Logarithm of the minimum tabulated `temperature
   *  (in log(T / K)). */
  double _min_logT;

  /*! @brief Inverse of the average logarithmic distance between tabulated
   *  temperature values (in 1 / log(T / K)). */
  double _inverse_avg_dlogT;


public:
  ChiantiRecombinationRates(const bool apply_diffuse_field);

  double get_recombination_rate_chianti(const int ion, const double temperature) const;

  virtual double get_recombination_rate(const int_fast32_t ion,
                                        const double temperature) const;

   // Wood, Mathis & Ercolano (2004), sections 3.3 and 7
  // equation (24)
  virtual double get_hydrogen_ground_state_recombination_rate(const double temperature) const {
    const double T4 = temperature * 1.e-4;
    const double alpha_1_H = 1.58e-13 * std::pow(T4, -0.53);
    return alpha_1_H;
  }

  // Wood, Mathis & Ercolano (2004), sections 3.3 and 7
  // equation (25)
  virtual double get_helium_ground_state_recombination_rate(const double temperature) const {
    const double T4 = temperature * 1.e-4;
    const double alpha_1_He = 1.54e-13 * std::pow(T4, -0.486);
    return alpha_1_He;
  } 
  
   /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   */
  ChiantiRecombinationRates(ParameterFile &params)
      : ChiantiRecombinationRates(
        params.get_value< bool >(
          "TaskBasedRadiationHydrodynamicsSimulation:diffuse field", true)
      ) {}
};

#endif // CHIANTIRECOMBINATIONRATES_HPP
