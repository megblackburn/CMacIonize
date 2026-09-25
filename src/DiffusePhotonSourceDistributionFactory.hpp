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
 * @file DiffusePhotonSourceDistributionFactory.hpp
 *
 * @brief Factory class for DiffusePhotonSourceDistribution instances.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#ifndef DIFFUSEPHOTONSOURCEDISTRIBUTIONFACTORY_HPP
#define DIFFUSEPHOTONSOURCEDISTRIBUTIONFACTORY_HPP

#include "Configuration.hpp"
#include "Error.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "RecombinationRates.hpp"
#include "CrossSections.hpp"

// non library dependent implementations
#include "TimeDependentDiffusePhotonSourceDistribution.hpp"

#include <string>
#include <typeinfo>

/**
 * @brief Factory class for DiffusePhotonSourceDistribution instances.
 */
class DiffusePhotonSourceDistributionFactory {
private:
  

public:


  static DiffusePhotonSourceDistribution *generate(ParameterFile &params,
                                            RecombinationRates &recombination_rates,
                                            CrossSections &cross_sections,
                                            Log *log = nullptr) {

    const std::string type = params.get_value< std::string >(
        "DiffusePhotonSourceDistribution:type", "TimeDependent");
    if (log) {
      log->write_info("Requested DiffusePhotonSourceDistribution type: ", type, ".");
    }

#ifndef HAVE_HDF5
    check_hdf5(type, log);
#endif
    if (type == "None") {
      return nullptr;
    } else if (type == "TimeDependent") {
      return new TimeDependentDiffusePhotonSourceDistribution(params, recombination_rates, cross_sections, log);
    } else {
      cmac_error("Unknown DiffusePhotonSourceDistribution type: \"%s\".",
                 type.c_str());
      return nullptr;
    }
    std::cout << "Diffuse Photon Source type: " << type << std::endl;
  }

  


  /**
   * @brief Restart the distribution from the given restart file.
   *
   * @param restart_reader Restart file to read from.
   * @param log Log to write logging info to.
   * @param use_tigress_like_injection True for TIGRESS-like supernova
   * injection, false for the legacy SILCC-like prescription.
   * @return Pointer to a newly created DiffusePhotonSourceDistribution implementation.
   * Memory management for the pointer needs to be done by the calling routine.
   */
  inline static DiffusePhotonSourceDistribution *restart(
      RestartReader &restart_reader, RecombinationRates &recombination_rates,
      CrossSections &cross_sections, Log *log = nullptr,
      const bool use_tigress_like_injection = true) {

    const std::string tag = restart_reader.read< std::string >();
    if (tag == typeid(TimeDependentDiffusePhotonSourceDistribution).name()) {
      return new TimeDependentDiffusePhotonSourceDistribution(restart_reader, recombination_rates, cross_sections, log);
    } else {
      cmac_error("Restarting is not supported for distribution type: \"%s\".",
                 tag.c_str());
      return nullptr;
    }
  }
};

#endif // DIFFUSEPHOTONSOURCEDISTRIBUTIONFACTORY_HPP
