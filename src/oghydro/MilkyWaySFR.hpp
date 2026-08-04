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
 * @file MilkyWaySFR.hpp
 *
 * @brief PhotonSourceSpectrum implementation for the Pauldrach, Hoffmann &
 * Lennon (2001) stellar model spectra.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MILKYWAYSFR_HPP
#define MILKYWAYSFR_HPP

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
class MilkyWaySFR : public PhotonSourceSpectrum {
private:
  /*! @brief Frequency bins. */
  std::vector< double > _file_time;

  /*! @brief Cumulative distribution of the spectrum. */
  std::vector< double > _file_SFR;

  /*! @brief Cumulative distribution of the spectrum. */
  std::vector< double > _file_SFR_error;

  double interpolate_value(double target_time, const std::vector< double > &time_axis, const std::vector< double > &data_axis) const;

public:
  MilkyWaySFR(std::string input_filename, Log *log = nullptr);

  MilkyWaySFR(std::string role, ParameterFile &params,
                              Log *log = nullptr);

  static std::string get_filename(std::string input_filename);
  double get_sfr_at_time(double target_time) const;

  /**
   * @brief Virtual destructor.
   */
  virtual ~MilkyWaySFR() {}

  virtual double get_random_frequency(RandomGenerator &random_generator, double temperature = 0.) const override;
  virtual double get_total_flux() const override;


};

#endif // MILKYWAYSFR_HPP
