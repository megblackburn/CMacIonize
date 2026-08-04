/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file ReadFUVLumToMass.hpp
 *
 * @brief PhotonSourceSpectrum implementation for the Pauldrach, Hoffmann &
 * Lennon (2001) stellar model spectra.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#ifndef READFUVLUMTOMASS_HPP
#define READFUVLUMTOMASS_HPP

#include "PhotonSourceSpectrum.hpp"

#include <string>
#include <vector>

class Log;
class ParameterFile;
class RandomGenerator;

/**
 * @brief Number of frequency bins used in the internal table.
 */

/**
 * @brief Gets the star formation rate at a specific time from the Swiggum et al profile
 */
class ReadFUVLumToMass : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency bins. */
  std::vector< double > _file_age;

  /*! @brief Cumulative distribution of the spectrum. */
  std::vector< double > _file_l_fuv_per_mass;


  double interpolate_value(double target_time, const std::vector< double > &time_axis, const std::vector< double > &data_axis) const;

public:
  ReadFUVLumToMass(std::string input_filename, Log *log = nullptr);

  ReadFUVLumToMass(std::string role, ParameterFile &params,
                              Log *log = nullptr);

  static std::string get_filename(std::string input_filename);
  double get_l_fuv_per_mass_at_time(double target_time) const;

  /**
   * @brief Virtual destructor.
   */
  virtual ~ReadFUVLumToMass() {}

  double get_total_flux() const;
  double get_random_frequency(RandomGenerator &random_generator, double temperature) const;


};

#endif // READFUVLUMTOMASS_HPP
