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
 * @file MagnetoHydroBoundaryConditions.hpp
 *
 * @brief Types of boundary conditions for the magneto-hydro.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#ifndef MAGNETOHYDROBOUNDARYCONDITIONS_HPP
#define MAGNETOHYDROBOUNDARYCONDITIONS_HPP

/**
 * @brief Types of boundary conditions implemented for the boundaries of the
 * box.
 */
enum MagnetoHydroBoundaryConditionType {
  /*! @brief A periodic boundary (only works if the grid is also periodic). */
  MAGNETOHYDRO_BOUNDARY_PERIODIC = 0,
  /*! @brief Reflective boundaries (elastic collisions are assumed at the
   *  boundaries). */
  MAGNETOHYDRO_BOUNDARY_REFLECTIVE,
  /*! @brief Inflow boundaries (material is assumed to flow in or out of the box
   *  at the same rate it flows near the boundary). */
  MAGNETOHYDRO_BOUNDARY_INFLOW,
  /*! @brief Outflow boundaries (material is allowed to leave the box, but
   *  cannot enter it). */
  MAGNETOHYDRO_BOUNDARY_OUTFLOW,
  /*! @brief Bondi inflow boundary conditions. */
  MAGNETOHYDRO_BOUNDARY_BONDI,
  /*! @brief Invalid boundaries selected. */
  MAGNETOHYDRO_BOUNDARY_INVALID
};

#endif // MAGNETOHYDROBOUNDARYCONDITIONS_HPP
