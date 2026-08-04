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
 * @file MagnetoHydroBoundaryManager.hpp
 *
 * @brief MagnetoHydro boundary management.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#ifndef MAGNETOHYDROBOUNDARYMANAGER_HPP
#define MAGNETOHYDROBOUNDARYMANAGER_HPP

#include "BondiMagnetoHydroBoundary.hpp"
#include "MagnetoHydroBoundary.hpp"
#include "ParameterFile.hpp"
#include "TravelDirections.hpp"

/**
 * @brief MagnetoHydro boundary management.
 */
class MagnetoHydroBoundaryManager {
private:
  /*! @brief MagnetoHydroBoundary for each boundary direction. */
  MagnetoHydroBoundary *_boundaries[6];

  /**
   * @brief Get a MagnetoHydroBoundary with the given type name.
   *
   * @param type Type of boundary condition.
   * @param params ParameterFile to read from.
   * @return Pointer to a new MagnetoHydroBoundary object.
   */
  inline static MagnetoHydroBoundary *get_boundary(const std::string type,
                                                   ParameterFile &params) {
    if (type == "bondi") {
      return new BondiMagnetoHydroBoundary(params);
    } else if (type == "inflow") {
      return new InflowMagnetoHydroBoundary();
    } else if (type == "outflow") {
      return new OutflowMagnetoHydroBoundary();
    } else if (type == "outflowBfield") {
      return new OutflowBfieldMagnetoHydroBoundary();
    } else if (type == "periodic") {
      return nullptr;
    } else if (type == "reflective") {
      return new ReflectiveMagnetoHydroBoundary();
    } else {
      cmac_error("Unknown magneto-hydro boundary type: \"%s\"!", type.c_str());
      return nullptr;
    }
  }

public:
  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   */
  inline MagnetoHydroBoundaryManager(ParameterFile &params)
      : _boundaries{
            get_boundary(params.get_value< std::string >(
                             "MagnetoHydroBoundaryManager:boundary x high", "inflow"),
                         params),
            get_boundary(params.get_value< std::string >(
                             "MagnetoHydroBoundaryManager:boundary x low", "inflow"),
                         params),
            get_boundary(params.get_value< std::string >(
                             "MagnetoHydroBoundaryManager:boundary y high", "inflow"),
                         params),
            get_boundary(params.get_value< std::string >(
                             "MagnetoHydroBoundaryManager:boundary y low", "inflow"),
                         params),
            get_boundary(params.get_value< std::string >(
                             "MagnetoHydroBoundaryManager:boundary z high", "inflow"),
                         params),
            get_boundary(params.get_value< std::string >(
                             "MagnetoHydroBoundaryManager:boundary z low", "inflow"),
                         params)} {}

  /**
   * @brief Destructor.
   */
  inline ~MagnetoHydroBoundaryManager() {
    for (uint_fast8_t i = 0; i < 6; ++i) {
      delete _boundaries[i];
    }
  }

  /**
   * @brief Get the MagnetoHydroBoundary for the given boundary direction.
   *
   * @param direction Boundary direction.
   * @return Corresponding MagnetoHydroBoundary.
   */
  inline MagnetoHydroBoundary &get_boundary_condition(int_fast8_t direction) const {
    cmac_assert_message(direction >= TRAVELDIRECTION_FACE_X_P &&
                            direction <= TRAVELDIRECTION_FACE_Z_N,
                        "Invalid boundary direction: %" PRIiFAST8, direction);
    MagnetoHydroBoundary *boundary = _boundaries[direction - TRAVELDIRECTION_FACE_X_P];
    if (boundary == nullptr) {
      cmac_error("Periodic boundaries are not properly linked!");
    }
    return *boundary;
  }
};

#endif // MAGNETOHYDROBOUNDARYMANAGER_HPP
