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
 * @file ReadFUVLumToMass.cpp
 *
 * @brief ReadFUVLumToMass implementation.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */

#include "ReadFUVLumToMass.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"
#include "RandomGenerator.hpp"
#include "UnitConverter.hpp"
#include "Utilities.hpp"
#include "FUVLumToMassDataLocation.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

/**
 * @brief Constructor.
 *
 * Reads in the correct model spectrum and resamples it on the 1000 bin internal
 * frequency grid.
 *
 * @param temperature Effective temperature of the star (in K).
 * @param surface_gravity Surface gravity of the star (in m s^-2).
 * @param log Log to write logging info to.
 */
ReadFUVLumToMass::ReadFUVLumToMass(std::string filename, Log *log) {
  std::string path_plus_filename = get_filename(filename);

  std::ifstream file(path_plus_filename);
  if (!file) {
    cmac_error("Error while opening file '%s'. FUV luminosity to mass conversion data does not exist!",
               path_plus_filename.c_str());
  }

  if (log) {
    log->write_status("Reading FUV luminosity to mass conversion from ", filename, "...");
  }
  std::string line;
  // skip the first three lines from the file, they contain column headers
  std::getline(file, line);


  while (std::getline(file, line)) {

    if (line.empty()) continue;


    std::stringstream linestream(line);
    std::string age, l_fuv_per_mass;

    if (std::getline(linestream, age, ',') && std::getline(linestream, l_fuv_per_mass, ',')) {
      double t =  std::stod(age);
      double l_fuv_per_mass_val = std::stod(l_fuv_per_mass);

      _file_age.push_back(t);
      _file_l_fuv_per_mass.push_back(l_fuv_per_mass_val);

    }
  }
   uint_fast32_t num_frequency = _file_age.size();
  if (num_frequency < 2) {
    cmac_error("Error: Not enough data points found in CSV file '%s' to interpolate!", filename.c_str());
  }


}

double ReadFUVLumToMass::get_l_fuv_per_mass_at_time(double target_time) const {
  return interpolate_value(target_time, _file_age, _file_l_fuv_per_mass);
}

/**
 * @brief ParameterFile constructor.
 *
 * Parameters are:
 *  - temperature: Temperature of the star (default: 4.e4 K)
 *  - surface gravity: Surface gravity of the star (default: 25. m s^-2)
 *
 * @param role Role the spectrum will fulfil in the simulation. Parameters are
 * read from the corresponding block in the parameter file.
 * @param params ParameterFile to read from.
 * @param log Log to write logging info to.
 */
ReadFUVLumToMass::ReadFUVLumToMass(std::string role, ParameterFile &params, Log *log)
    : ReadFUVLumToMass(
          params.get_value< std::string >(
              role + ":filename", "Kroupa_IMF_L_FUV_per_mass.csv"),
          log) {}


/**
 * @brief Get the name of the table containing the spectrum for the given
 * effective temperature and surface gravity.
 *
 * @param filename
 * @return Name of the file containing the requested SFH.
 */
std::string ReadFUVLumToMass::get_filename(std::string input_filename) {
  std::stringstream ss;
  ss << FUVLUMTOMASSDATALOCATION << input_filename;
  return ss.str();
}



double ReadFUVLumToMass::interpolate_value(double target_time, const std::vector< double > &time_axis, const std::vector< double > &data_axis) const {
  
  if (target_time <= time_axis.front()) {
    return data_axis.front();
  }
  if (target_time >= time_axis.back()){
    return data_axis.back();
  }

  auto it = std::upper_bound(time_axis.begin(), time_axis.end(), target_time);

  size_t idx2 = std::distance(time_axis.begin(), it);
  size_t idx1 = idx2 - 1;

  double t1 = time_axis[idx1];
  double t2 = time_axis[idx2];
  double f = (target_time-t1)/(t2 - t1);

  return data_axis[idx1] + f * (data_axis[idx2] - data_axis[idx1]);
}


double ReadFUVLumToMass::get_random_frequency(RandomGenerator &random_generator, double temperature) const {
    return 0.0; 
}

double ReadFUVLumToMass::get_total_flux() const {
    return 0.0; 
}