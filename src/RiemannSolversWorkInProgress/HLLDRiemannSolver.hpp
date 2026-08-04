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
 * @file HLLDRiemannSolver.hpp
 *
 * @brief HLLD Riemann solver.
 *
 * Based on the HLLD Riemann solver in Shadowfax (Vandenbroucke & De Rijcke,
 * 2016).
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#ifndef HLLDRIEMANNSOLVER_HPP
#define HLLDRIEMANNSOLVER_HPP

#include "Error.hpp"
#include "RiemannSolver.hpp"
#include "Utilities.hpp"

#include <algorithm>
#include <cfloat>
#include <cinttypes>
#include <cmath>

/**
 * @brief HLLD Riemann solver.
 */
class HLLDRiemannSolver : public RiemannSolver {
private:
  /*! @brief Polytropic index of the gas, @f$\gamma{}@f$. */
  const double _gamma;


  /*! @brief @f${\gamma{}-1}@f$. */
  const double _gm1;

  /*! @brief @f$\frac{\gamma{}+1}{2\gamma{}}@f$. */
  const double _gp1d2g;

  /*! @brief @f$\frac{1}{\gamma{}-1}@f$. */
  const double _odgm1;

  /*! @brief @f$\frac{\gamma{}-1}{2}@f$. */
  const double _gm1d2;

  /*! @brief @f$\frac{\gamma{}-1}{\gamma{}+1}@f$. */
  const double _gm1dgp1;

  /*! @brief @f$\frac{2\gamma{}}{\gamma{}-1}@f$. */
  const double _tgdgm1;

  /*! @brief @f$\frac{2}{\gamma{}+1}@f$. */
  const double _tdgp1;

  

  

  inline double get_fast_magnetosonic_speed(const double rho, const double p, const CoordinateVector<> &B) const { // mgb edit 21.04.2026 - from Minoshima et al.
    if (rho <= 0.) {
      return 0.;
    }
    const double cs_squared = _gamma * p / rho;
    const double B_squared = CoordinateVector<>::dot_product(B, B);
    const double va_squared = B_squared * _inv_mu0 / rho;
    const double vax_squared = (B.x() * B.x()) * _inv_mu0 / rho;

    const double cf_squared = 0.5 * ((cs_squared + va_squared) + std::sqrt(std::pow(cs_squared + va_squared, 2.) - 4. * cs_squared * vax_squared));

    return std::sqrt(cf_squared);
  }


  /**
  

public:
  /**
   * @brief Constructor.
   *
   * @param gamma Polytropic index of the gas.
   */
  inline HLLDRiemannSolver(const double gamma, const double mu0, const double cleaning_speed) // mgb edit 21.04.2026 - add mu0 and cleaning speed parameters for MHD
      : _gamma(std::max(gamma, 1.00000001)), _gm1(_gamma - 1.), // mgb edit 21.04.2026 from 2/(gamma-1) to (gamma-1) for HLLD
        _odgm1(1. / (_gamma - 1.)),
        _mu0(mu0), _inv_mu0(1/mu0), _inv_sqrt_mu0(1./std::sqrt(mu0)), _cleaning_speed(cleaning_speed), _cleaning_speed_squared(cleaning_speed * cleaning_speed) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~HLLDRiemannSolver() {}

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
   */
  virtual void solve_for_flux_MHD(const double rhoL, const CoordinateVector<> uL,
                              const double PL, const CoordinateVector<>BL, const double BL_scalar, 
                              const double rhoR, const CoordinateVector<> uR, 
                              const double PR, const CoordinateVector<> BR, const double BR_scalar,
                              double &mflux, CoordinateVector<> &pflux,
                              double &Eflux, CoordinateVector<> &Bflux, double &B_scalarflux, const CoordinateVector<> normal, // mgbb edit 21.04.2026 - add Bflux and B_scalar_flux parameters for MHD
                              const CoordinateVector<> vface = 0.) const {

    
    // check input values
    cmac_assert_message(
        rhoL == rhoL && rhoL >= 0. && uL.x() == uL.x() && uL.y() == uL.y() &&
            uL.z() == uL.z() && PL == PL && PL >= 0. && rhoR == rhoR &&
            rhoR >= 0. && uR.x() == uR.x() && uR.y() == uR.y() &&
            uR.z() == uR.z() && PR == PR && PR >= 0.,
        "Invalid Riemann solver input: rhoL: %g, uL: %g %g %g, PL: %g, rhoR: "
        "%g, uR: %g %g %g, PR: %g",
        rhoL, uL.x(), uL.y(), uL.z(), PL, rhoR, uR.x(), uR.y(), uR.z(), PR);


    

    const double B2L = CoordinateVector<>::dot_product(BL, BL);
    const double B2R = CoordinateVector<>::dot_product(BR, BR);
    const double PtotL = PL + 0.5 * B2L * _inv_mu0;
    const double PtotR = PR + 0.5 * B2R * _inv_mu0;

    const double Bnorm = 0.5 * (CoordinateVector<>::dot_product(BL, normal) + CoordinateVector<>::dot_product(BR, normal));

    const double cfL = get_fast_magnetosonic_speed(rhoL, PL, BL);
    const double cfR = get_fast_magnetosonic_speed(rhoR, PR, BR);
    const cf_max = std::max(cfL, cfR);

    const double max_u = std::max(std::abs(uL.x()), std::abs(uR.x()));
    const double min_u = std::min(std::abs(uL.x()), std::abs(uR.x()));
    const double S_R = std::max(0., max_u + cf_max); // right fast magnetosonic
    const double S_L = std::min(0., min_u - cf_max); // left fast magnetosonic 


    const double S_M = ((S_R-uR.x())*rhoR*uR.x()-(S_L-uL.x())*rhoL*uL.x()-(PtotR-PtotL))/((S_R - uR.x())*rhoR-(S_L - uL.x())*rhoL); // middle speed (discontinuity)

    const double total_P = ((S_R-uR.x())*rhoR*PtotL - (S_L-uL.x())*rhoL*PtotR + (rhoL*rhoR*(S_R-uR.x())*(S_L-uL.x())*(uR.x()-uL.x())))/((S_R - uR.x())*rhoR-(S_L - uL.x())*rhoL);

    const double rhoL_star = rhoL * (S_L - uL.x()) / (S_L - S_M);
    const double rhoR_star = rhoR * (S_R - uR.x()) / (S_R - S_M);

    const double S_L_star = S_M - std::abs(Bnorm) * _inv_sqrt_mu0 / std::sqrt(rhoL_star); // left alfven speed
    const double S_R_star = S_M + std::abs(Bnorm) * _inv_sqrt_mu0 / std::sqrt(rhoR_star); // right alfven speed
      // UstarL = [rhoL_star, rhoL_star*S_M, rhoL_star*uL.y(), rhoL_star*uL.z(), BL.x(), BL.y(), BL.z(), PL/(gamma-1.) + 0.5*rhoL_star*S_M*S_M + 0.5*B2L*_inv_mu0] check this
    
    if (S_L > 0.) { // right-going waves only
      mflux = rhoL * uL.x();
      pflux = rhoL * uL.x() * uL + PtotL * normal - BL * (CoordinateVector<>::dot_product(BL, normal) * _inv_mu0);
      Eflux = (PL / _gm1 + 0.5 * rhoL * CoordinateVector<>::dot_product(uL, uL) + 0.5 * B2L * _inv_mu0)
      Bflux.x() = 0.0;
      Bflux.y() = uL.x() * BL.y() - uL.y() * BL.x();
      Bflux.z() = uL.x() * BL.z() - uL.z() * BL.x();
      //B_scalarflux = 0.;
    } else if (S_R > 0.) { // left-going waves only
      double factor = S_R * ()
      mflux = rhoR * uR.x();
      pflux = rhoR * uR.x() * uR + PtotR * normal - BR * (CoordinateVector<>::dot_product(BR, normal) * _inv_mu0);
      Eflux = (PR / _gm1 + 0.5 * rhoR * CoordinateVector<>::dot_product(uR, uR) + 0.5 * B2R * _inv_mu0);
      Bflux.x() = 0.0;
      Bflux.y() = uR.x() * BR.y() - uR.y() * BR.x();
      Bflux.z() = uR.x() * BR.z() - uR.z() * BR.x();
      //B_scalarflux = 0.;
    } else if (S_L_star > 0.) { // left-going waves and middle discontinuity
      
    } else if (S_R_star > 0.) { // right-going waves and middle discontinuity
      
    } else if (S_M > 0.) { // middle discontinuity
     
    } else { // all waves

    }
    
    else { // all waves

  }
};

#endif // HLLDRIEMANNSOLVER_HPP
