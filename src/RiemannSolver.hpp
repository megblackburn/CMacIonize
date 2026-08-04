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
 * @file RiemannSolver.hpp
 *
 * @brief General interface for Riemann solvers.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef RIEMANNSOLVER_HPP
#define RIEMANNSOLVER_HPP

#include "CoordinateVector.hpp"

#include <cinttypes>

/**
 * @brief General interface for Riemann solvers.
 */
class RiemannSolver {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~RiemannSolver() {}

  /**
   * @brief Solve the Riemann problem with the given left and right state and
   * get the resulting flux accross an interface.
   *
   * @param rhoL Left state density.
   * @param uL Left state velocity.
   * @param PL Left state pressure.
   * @param rhoR Right state density.
   * @param uR Right state velocity.
   * @param PR Right state pressure.
   * @param mflux Mass flux solution.
   * @param pflux Momentum flux solution.
   * @param Eflux Energy flux solution.
   * @param normal Surface normal of the interface.
   * @param vface Velocity of the interface, used to boost the fluxes.
   * @param BL Left state magnetic field.
   * @param BR Right state magnetic field.
   * @param BL_scalar Left state scalar magnetic field.
   * @param BR_scalar Right state scalar magnetic field.
   * @param Bflux Magnetic field flux solution.
   * @param B_scalar_flux Scalar magnetic field flux solution.

   */
  virtual void solve_for_flux(const double rhoL, const CoordinateVector<> uL,
                              const double PL, const double rhoR,
                              const CoordinateVector<> uR, const double PR,
                              double &mflux, CoordinateVector<> &pflux,
                              double &Eflux, double &EintFlux, const CoordinateVector<> normal,
                              const CoordinateVector<> vface = 0.) const = 0;

 
 virtual void solve_for_flux_MHD(const double rhoL, 
                              const CoordinateVector<> uL, const double PL, 
                              const CoordinateVector<> BL, const double BL_scalar,
                              const double rhoR,
                              const CoordinateVector<> uR, const double PR,
                              const CoordinateVector<> BR, const double BR_scalar, 
                              double &mflux, CoordinateVector<> &pflux,
                              double &Eflux, double &EintFlux, CoordinateVector<> &Bflux,
                              double &B_scalar_flux,
                              const CoordinateVector<> normal,
                              double current_etot,
                              double _mach_limit,
                              double current_eint,
                              const CoordinateVector<> vface = 0.) const = 0;                             
};

#endif // RIEMANNSOLVER_HPP
