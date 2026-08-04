/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file HydroBoundary.hpp
 *
 * @brief Hydro boundary conditions related functionality.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#ifndef MAGNETOHYDROBOUNDARY_HPP
#define MAGNETOHYDROBOUNDARY_HPP

#include "MagnetoHydroVariables.hpp"

/**
 * @brief General interface for hydrodynamical boundary conditions.
 */
class MagnetoHydroBoundary {
public:
  /**
   * @brief Virtual destructor.
   */
  inline virtual ~MagnetoHydroBoundary() {}

  /**
   * @brief Get the right state primitive variables corresponding to the given
   * left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables).
   */
  virtual MagnetoHydroVariables get_right_state_gradient_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR,
      const MagnetoHydroVariables &left_state) const = 0;

  /**
   * @brief Get the right state primitive variables and gradients corresponding
   * to the given left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables and
   * gradients).
   */
  virtual MagnetoHydroVariables
  get_right_state_flux_variables(const int_fast8_t i,
                                 const int_fast8_t orientation,
                                 const CoordinateVector<> posR,
                                 const MagnetoHydroVariables &left_state) const = 0;
};

/**
 * @brief Inflow hydro boundary.
 */
class InflowMagnetoHydroBoundary : public MagnetoHydroBoundary {
public:
  /**
   * @brief Get the right state primitive variables corresponding to the given
   * left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables).
   */
  virtual MagnetoHydroVariables get_right_state_gradient_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR, const MagnetoHydroVariables &left_state) const {

    MagnetoHydroVariables right_state;
    for (uint_fast8_t j = 0; j < 10; ++j) {
      right_state.primitives(j) = left_state.primitives(j);
    }

    return right_state;
  }

  /**
   * @brief Get the right state primitive variables and gradients corresponding
   * to the given left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables and
   * gradients).
   */
  virtual MagnetoHydroVariables get_right_state_flux_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR, const MagnetoHydroVariables &left_state) const {

    MagnetoHydroVariables right_state;

    for (uint_fast8_t j = 0; j < 10; ++j) {
      right_state.primitives(j) = left_state.primitives(j);
      right_state.primitive_gradients(j) = left_state.primitive_gradients(j);
    }

    return right_state;
  }
};

/**
 * @brief Outflow magneto-hydro boundary.
 */
class OutflowMagnetoHydroBoundary : public MagnetoHydroBoundary {
public:
  /**
   * @brief Get the right state primitive variables corresponding to the given
   * left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state magneto-hydro variables.
   * @return Corresponding right state (only containing primitive variables).
   */
  virtual MagnetoHydroVariables get_right_state_gradient_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR, const MagnetoHydroVariables &left_state) const {

    MagnetoHydroVariables right_state;

    for (uint_fast8_t j = 0; j < 10; ++j) {
      right_state.primitives(j) = left_state.primitives(j);
    }
    // we need to reverse the velocity component aligned with the surface
    // normal, but only if it is entering the box
    if (orientation * right_state.primitives(1 + i) < 0.) {
      right_state.primitives(1 + i) = -right_state.primitives(1 + i);
    }

    return right_state;
  }

  /**
   * @brief Get the right state primitive variables and gradients corresponding
   * to the given left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables and
   * gradients).
   */
  virtual MagnetoHydroVariables get_right_state_flux_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR, const MagnetoHydroVariables &left_state) const {

    MagnetoHydroVariables right_state;

    for (uint_fast8_t j = 0; j < 10; ++j) {
      right_state.primitives(j) = left_state.primitives(j);
      right_state.primitive_gradients(j) = left_state.primitive_gradients(j);
    }
    // we need to reverse the velocity component aligned with the surface
    // normal, but only if it is entering the box
    if (orientation * right_state.primitives(1 + i) < 0.) {
      right_state.primitives(1 + i) = -right_state.primitives(1 + i);
      right_state.primitive_gradients(1 + i) = CoordinateVector<>(0.);
    }

    return right_state;
  }
};

/**
 * @brief TIGRESS-Equivalent Stratified Vertical Boundary Condition.
 * Combines zero-gradient continuous transport for fluid fields and transverse B,
 * with strict solenoidal zero-derivative tracking for the normal field.
 */
class OutflowBfieldMagnetoHydroBoundary : public MagnetoHydroBoundary {
public:
  virtual MagnetoHydroVariables get_right_state_gradient_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR, const MagnetoHydroVariables &left_state) const {

    MagnetoHydroVariables right_state;
    
    for (uint_fast8_t j = 0; j < 5; ++j) {
      right_state.primitives(j) = left_state.primitives(j);
    }


    right_state.primitives(5) = left_state.primitives(5); // Bx
    right_state.primitives(6) = left_state.primitives(6); // By
    right_state.primitives(7) = left_state.primitives(7); // Bz


    right_state.primitives(8) = left_state.primitives(8); 
    right_state.primitives(9) = left_state.primitives(9); // Eint

    return right_state;
  }

  virtual MagnetoHydroVariables get_right_state_flux_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR, const MagnetoHydroVariables &left_state) const {

    MagnetoHydroVariables right_state;

    // Symmetrically copy all primitive tracks over to the face flux calculator
    for (uint_fast8_t j = 0; j < 10; ++j) {
      right_state.primitives(j) = left_state.primitives(j);
      right_state.primitive_gradients(j) = left_state.primitive_gradients(j);
    }

    // Prevent back-reflection of gradients at the outer edge boundary plane
    right_state.primitive_gradients(1 + i) = CoordinateVector<>(0.);

    return right_state;
  }
};

/**
 * @brief Reflective magneto-hydro boundary.
 */
class ReflectiveMagnetoHydroBoundary : public MagnetoHydroBoundary {
public:
  /**
   * @brief Get the right state primitive variables corresponding to the given
   * left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @return Corresponding right state (only containing primitive variables).
   */
  virtual MagnetoHydroVariables get_right_state_gradient_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR, const MagnetoHydroVariables &left_state) const {

    MagnetoHydroVariables right_state;
    for (uint_fast8_t j = 0; j < 10; ++j) {
      right_state.primitives(j) = left_state.primitives(j);
    }
    // we need to reverse the velocity component aligned with the surface normal
    right_state.primitives(1 + i) = -right_state.primitives(1 + i);

    // mgb edit 21.04.2026: we also need to reverse the magnetic field component aligned with the surface normal
    right_state.primitives(5 + i) = -right_state.primitives(5 + i);

    // mgb edit 21.04.2026: we also need to reverse the magnetic field scalar
    right_state.primitives(8) = -right_state.primitives(8);

    return right_state;
  }

  /**
   * @brief Get the right state primitive variables and gradients corresponding
   * to the given left state boundary ghost.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param orientation Interface orientation: negative (-1) or positive (1).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state magneto-hydro variables.
   * @return Corresponding right state (only containing primitive variables and
   * gradients).
   */
  virtual MagnetoHydroVariables get_right_state_flux_variables(
      const int_fast8_t i, const int_fast8_t orientation,
      const CoordinateVector<> posR, const MagnetoHydroVariables &left_state) const {

    MagnetoHydroVariables right_state;
    
    for (uint_fast8_t j = 0; j < 10; ++j) {
      right_state.primitives(j) = left_state.primitives(j);
      right_state.primitive_gradients(j) = left_state.primitive_gradients(j);
    }
    // we need to reverse the velocity component aligned with the surface normal
    right_state.primitives(1 + i) = -right_state.primitives(1 + i);
    // mgb edit 21.04.2026: we also need to reverse the magnetic field component aligned with the surface normal
    right_state.primitives(5 + i) = -right_state.primitives(5 + i);
    // mgb edit 21.04.2026: we also need to reverse the magnetic field scalar
    right_state.primitives(8) = -right_state.primitives(8);

    // we need to invert all gradients that are aligned with the interface
    // normal: idx
    // however, we do not invert the gradient of the velocity aligned with the
    // interface normal, as the velocity itself also changes sign
    // mgb edit 21.04.2026: we want to flip all gradients but then unflip the velocity, B field and B scalar -following original logic
    for (uint_fast8_t j = 0; j < 10; ++j) {
      right_state.primitive_gradients(j)[i] =
          -right_state.primitive_gradients(j)[i];
    }
    
    right_state.primitive_gradients(1 + i)[i] =
        -right_state.primitive_gradients(1 + i)[i];
    right_state.primitive_gradients(5 + i)[i] =
        -right_state.primitive_gradients(5 + i)[i];
    right_state.primitive_gradients(8)[i] =
        -right_state.primitive_gradients(8)[i];

    return right_state;
  }
};

#endif // MAGNETOHYDROBOUNDARYWIP_HPP
