/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file RadiativeCoolingHandlerFactory.hpp
 *
 * @brief Factory for RadiativeCoolingHandler instances.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef RADIATIVECOOLINGHANDLERFACTORY_HPP
#define RADIATIVECOOLINGHANDLERFACTORY_HPP

#include "RadiativeCoolingHandler.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"

// implementations
#include "DeRijckeRadiativeCooling.hpp"
#include "CloudyRadiativeCooling.hpp"

class RadiativeCoolingHandlerFactory {
public:
  /**
   * @brief Generate a RadiativeCoolingHandler instance of the type found in
   * the given ParameterFile.
   *
   * Supported types are (default: None):
   *  - FixedValue: Fixed reemission probability and frequency handler.
   *  - Physical: Physical diffuse reemission handler.
   *
   * @param cross_sections Photoionization cross sections to use.
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   * @return Pointer to a newly created RadiativeCoolingHandler instance (or a
   * null pointer if no RadiativeCoolingHandler was selected in the
   * ParameterFile). Memory management for this pointer should be done by the
   * calling routine.
   */
  inline static RadiativeCoolingHandler *
  generate(ParameterFile &params,
           Log *log = nullptr) {

    // 1. Fetch the type from the parameters
    const std::string type = params.get_value< std::string >(
        "RadiativeCoolingHandler:type", "Cloudy");

    if (log) {
      log->write_info("Requested RadiativeCoolingHandler type: ", type);
    }
    
    std::cout << "Radiative Cooling type: " << type << std::endl;


    // 3. Factory returns - allocating the actual handler classes based on type
    if (type == "DeRijcke") {
      return new DeRijckeRadiativeCooling();
    } else if (type == "Cloudy") {
      return new CloudyRadiativeCooling();
    } else {
      cmac_error("Unknown RadiativeCoolingHandler type: \"%s\"!",
                 type.c_str());
      return nullptr;
    }
  }
};

#endif // RADIATIVECOOLINGHANDLERFACTORY_HPP
