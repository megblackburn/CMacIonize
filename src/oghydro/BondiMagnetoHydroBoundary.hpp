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
 * @file BondiMagnetoHydroBoundary.hpp
 *
 * @brief Bondi magneto-hydro boundaries.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#ifndef BONDIMAGNETOHYDROBOUNDARY_HPP
#define BONDIMAGNETOHYDROBOUNDARY_HPP

#include "BondiProfileMHD.hpp"
#include "MagnetoHydroBoundary.hpp"
/**
 * @brief Bondi magneto-hydro boundaries.
 */
class BondiMagnetoHydroBoundary : public MagnetoHydroBoundary {
  /*! @brief BondiProfile to use. */
  const BondiProfileMHD _bondi_profile;

public:
  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   */
  inline BondiMagnetoHydroBoundary(ParameterFile &params) : _bondi_profile(params) {}

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

    double rhoR, PR, nfrac;
    CoordinateVector<> uR;
    _bondi_profile.get_magnetohydrodynamic_variables(posR, rhoR, uR, PR, nfrac);
    MagnetoHydroVariables right_state;
    right_state.set_primitives_density(rhoR);
    right_state.set_primitives_velocity(uR);
    right_state.set_primitives_pressure(PR);

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

    double rhoR, PR, nfrac;
    CoordinateVector<> uR;
    _bondi_profile.get_magnetohydrodynamic_variables(posR, rhoR, uR, PR, nfrac);
    MagnetoHydroVariables right_state;
    right_state.set_primitives_density(rhoR);
    right_state.set_primitives_velocity(uR);
    right_state.set_primitives_pressure(PR);

    return right_state;
  }
};

#endif // BONDIMAGNETOHYDROBOUNDARY_HPP
