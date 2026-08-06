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
#include "ParameterFile.hpp"
#include "ExactRiemannSolver.hpp"

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

  /*! @brief Magnetic permeability of free space. */
  const double _mu0;

  /*! @brief Maximum velocity for the solver. */
  const double _max_velocity;

  /*! @brief Factor for the cleaning speed. */
  const double _cleaning_speed_factor;

  /*! @brief @f${\gamma{}-1}@f$. */
  const double _gm1;


  /*! @brief @f$\frac{1}{\gamma{}-1}@f$. */
  const double _odgm1;

  /*! @brief Inverse of the magnetic permeability of free space. */
  const double _inv_mu0;

  /*! @brief Inverse of the square root of the magnetic permeability of free space. */
  const double _inv_sqrt_mu0;

  /*! @brief the cleaning speed. */
  const double _cleaning_speed;

  /*! @brief the cleaning speed squared. */
  const double _cleaning_speed_squared;

  /*! @brief Riemann Solver to revert to if no B fields*/
  const ExactRiemannSolver _riemann_solver_hd;

  
  inline double get_limited_alfven_speed(double B_squared, double rho, double v_limit) const {
   // v_limit *= 10. // take v_limit to be 10 * max velocity
    double alfspeed_squared = B_squared/(_mu0*rho);
    double alfspeed = std::sqrt(alfspeed_squared);
    double alfspeed_limited = alfspeed/std::sqrt(1.0+std::pow(alfspeed/v_limit,2.));

    return alfspeed_limited * alfspeed_limited; 
  }
  
  /**
   * @brief limit the fast magnetosonic speed 
   *
   * @param B_mag Magnitude of the magnetic field.
   * @param rho Density.
   * @param v_limit Maximum velocity.
   * @return Limited Alfven speed.
   */

  inline double get_fast_magnetosonic_speed(const double rho, const double p, const CoordinateVector<> &B, const double max_velocity) const { // mgb edit 21.04.2026 - from Minoshima et al.
    if (rho <= 0.) {
      return 0.;
    }
    const double cs_squared = _gamma * p / rho;
    const double B_squared = CoordinateVector<>::dot_product(B, B);
    const double Bx_squared = B.x() * B.x();
    double va_squared = B_squared /(_mu0 * rho);
    double vax_squared = (Bx_squared) /(_mu0 * rho);

    if (va_squared > _max_velocity*_max_velocity) {
        va_squared = get_limited_alfven_speed(B_squared, rho, max_velocity);
    }
    if (vax_squared > _max_velocity*_max_velocity) {
        vax_squared = get_limited_alfven_speed(B_squared, rho, max_velocity);
    }


    const double cf_squared = 0.5 * ((cs_squared + va_squared) + std::sqrt(std::pow(cs_squared + va_squared, 2.) - 4. * cs_squared * vax_squared));

  //  //std::cout << "get_fast_magnetosonic_speed: cs_squared = " << cs_squared << ", va_squared = " << va_squared << ", vax_squared = " << vax_squared << ", cf_squared = " << cf_squared << std::endl;

    return std::sqrt(cf_squared);
  }



  

public:
  /**
   * @brief Constructor.
   *
   * @param gamma Polytropic index of the gas.
   */
  inline HLLDRiemannSolver(const double gamma, const double max_velocity, const double cleaning_speed_factor) // mgb edit 21.04.2026 - add mu0 and cleaning speed parameters for MHD
      :  _gamma(std::max(gamma, 1.00000001)), _mu0(PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_MAGNETC_PERMEABILITY)), _max_velocity(max_velocity), _cleaning_speed_factor(cleaning_speed_factor), _gm1(_gamma - 1.), // mgb edit 21.04.2026 from 2/(gamma-1) to (gamma-1) for HLLD
        _odgm1(1. / (_gamma - 1.)),
         _inv_mu0(1/_mu0), _inv_sqrt_mu0(1./std::sqrt(_mu0)), _cleaning_speed(max_velocity * cleaning_speed_factor), _cleaning_speed_squared(_cleaning_speed * _cleaning_speed), _riemann_solver_hd(gamma) {}

  /**
  * @brief ParameterFile constructor.
   *
   * The following parameters are read:
   *  - polytropic index: Polytropic index of the gas (default: 5. / 3.)
   *  - neutral temperature: Assumed neutral temperature for the gas (default:
   *    100. K)
   *  - ionised temperature: Assumed ionised temperature for the gas (default:
   *    1.e4 K)
   *  - maximum velocity: Maximum allowed velocity for the gas. The gas velocity
   *    is capped at this value (default: 1.e99 m s^-1)
   *  - do explicit heating: Enable explicit radiation heating? (default: false)
   *
   * @param params ParameterFile to read from.
   */
  inline HLLDRiemannSolver(ParameterFile &params)
      : HLLDRiemannSolver(params.get_value< double >("Hydro:polytropic index", 5. / 3.),
              params.get_physical_value< QUANTITY_VELOCITY >(
                  "Hydro:maximum velocity", "1.e99 m s^-1"),
              
              params.get_value<double>("Hydro:cleaning speed factor",1.0) // mgb edit 21.04.2026 - factor to multiply the cleaning speed by (default: 1.0)
              ) {}
  /**
   * @brief Virtual destructor.
   */
  virtual ~HLLDRiemannSolver() {}

  virtual void solve_for_flux(const double rhoL, const CoordinateVector<> uL,
                              const double PL, const double rhoR,
                              const CoordinateVector<> uR, const double PR,
                              double &mflux, CoordinateVector<> &pflux,
                              double &Eflux, double &EintFlux, const CoordinateVector<> normal,
                              const CoordinateVector<> vface) const override {
    // This solver is for MHD. If the code calls this hydro version, 
    // it's likely a logic error elsewhere.
    cmac_error("HLLDRiemannSolver::solve_for_flux (Hydro) called! Use solve_for_flux_MHD.");
  }

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
   * @param B_scalarL Left state scalar magnetic field.
   * @param B_scalarR Right state scalar magnetic field.
   * @param Bflux Magnetic field flux solution.
   * @param B_scalar_flux Scalar magnetic field flux solution.
   */
  virtual void solve_for_flux_MHD(const double rhoL, const CoordinateVector<> uL,
                              const double PL, const CoordinateVector<>BL, const double BL_scalar, 
                              const double rhoR, const CoordinateVector<> uR, 
                              const double PR, const CoordinateVector<> BR, const double BR_scalar,
                              double &mflux, CoordinateVector<> &pflux,
                              double &Eflux, double &EintFlux, CoordinateVector<> &Bflux, double &B_scalarflux, const CoordinateVector<> normal, double current_etot, double _mach_limit, double current_eint,// mgbb edit 21.04.2026 - add Bflux and B_scalar_flux parameters for MHD
                              const CoordinateVector<> vface=0.) const { // mgb note: vface is currently ignored but if nonzero then need to do uL - vface and uR - vface for all velocity components in the following code

    // check input values
    cmac_assert_message(
        rhoL == rhoL && rhoL >= 0. && uL.x() == uL.x() && uL.y() == uL.y() &&
            uL.z() == uL.z() && PL == PL && PL >= 0. && rhoR == rhoR &&
            rhoR >= 0. && uR.x() == uR.x() && uR.y() == uR.y() &&
            uR.z() == uR.z() && PR == PR && PR >= 0.,
        "Invalid Riemann solver input: rhoL: %g, uL: %g %g %g, PL: %g, rhoR: "
        "%g, uR: %g %g %g, PR: %g",
        rhoL, uL.x(), uL.y(), uL.z(), PL, rhoR, uR.x(), uR.y(), uR.z(), PR);


    const bool debug_flag = false;

    double left_conserved[10];
    double right_conserved[10];
    double left_flux[10];
    double right_flux[10];
    double star_conserved[10];
    double starstar_conserved[10];

    CoordinateVector<> t1,t2;
    CoordinateVector<>temp;

    if (std::abs(normal.x()) < 0.7) {
      temp = CoordinateVector<>(1.,0.,0.);
    } else {
      temp = CoordinateVector<>(0.,1.,0.);
    }

    t1 = CoordinateVector<>::cross_product(temp, normal);
    t1/=t1.norm();
    t2 = CoordinateVector<>::cross_product(normal, t1);

    const double BFLOAT_EPSILON = 1e-5; // small number to prevent division by zero in case of very low magnetic fields, should be small enough to not affect results when magnetic field is significant

    const double B2L = CoordinateVector<>::dot_product(BL, BL);
    const double B2R = CoordinateVector<>::dot_product(BR, BR);
    const double Bu_L = CoordinateVector<>::dot_product(BL,uL);
    const double Bu_R = CoordinateVector<>::dot_product(BR,uR);
    
    const double max_B2 = std::max(B2L, B2R);
    const double P_gas_max = std::max(PL,PR);
    const double magnetic_vacuum_threshold = 1e-20 * P_gas_max * _mu0;

    const double v2L = CoordinateVector<>::dot_product(uL, uL);
    const double v2R = CoordinateVector<>::dot_product(uR, uR);

    const double uLn = CoordinateVector<>::dot_product(uL, normal);
    const double uRn = CoordinateVector<>::dot_product(uR, normal);
    const double uLt1 = CoordinateVector<>::dot_product(uL, t1);
    const double uRt1 = CoordinateVector<>::dot_product(uR, t1);
    const double uLt2 = CoordinateVector<>::dot_product(uL, t2);
    const double uRt2 = CoordinateVector<>::dot_product(uR, t2);

    const double BLn = CoordinateVector<>::dot_product(BL, normal);
    const double BRn = CoordinateVector<>::dot_product(BR, normal);
    const double BLt1 = CoordinateVector<>::dot_product(BL, t1);
    const double BRt1 = CoordinateVector<>::dot_product(BR, t1);
    const double BLt2 = CoordinateVector<>::dot_product(BL, t2);
    const double BRt2 = CoordinateVector<>::dot_product(BR, t2);

    const double EtotL =  (PL / _gm1) + (0.5 * rhoL * v2L) + 0.5 * B2L * _inv_mu0;
    const double EtotR = (PR / _gm1) + (0.5 * rhoR * v2R) + 0.5 * B2R * _inv_mu0; // double check this uses Einternal

    const double EintL = PL / _gm1;
    const double EintR = PR / _gm1;

    const double PtotL = PL + 0.5 * B2L * _inv_mu0;
    const double PtotR = PR + 0.5 * B2R * _inv_mu0;

    const double Bnorm = 0.5 * (CoordinateVector<>::dot_product(BL, normal) + CoordinateVector<>::dot_product(BR, normal));

    const double cfL = get_fast_magnetosonic_speed(rhoL, PL, BL, _max_velocity);
    const double cfR = get_fast_magnetosonic_speed(rhoR, PR, BR, _max_velocity);
    const double cf_max = std::max(cfL, cfR);
    const double v_max = std::max(std::abs(uLn), std::abs(uRn));
    const double _cleaning_speed_max = std::min(v_max + cf_max, _max_velocity * _cleaning_speed_factor);

    left_conserved[0] = rhoL;
    left_conserved[1] = rhoL * uLn;
    left_conserved[2] = rhoL * uLt1;
    left_conserved[3] = rhoL * uLt2;
    left_conserved[4] = EtotL;
    left_conserved[5] = BLn;
    left_conserved[6] = BLt1;
    left_conserved[7] = BLt2;
    left_conserved[8] = BL_scalar;
    left_conserved[9] = EintL;

    left_flux[0] = rhoL * uLn;
    left_flux[1] = left_conserved[1] * uLn + PtotL - BLn * BLn * _inv_mu0;
    left_flux[2] = left_conserved[2] * uLn - BLn * BLt1 * _inv_mu0;
    left_flux[3] = left_conserved[3] * uLn - BLn * BLt2 * _inv_mu0;
    left_flux[4] = (EtotL + PtotL) * uLn - BLn * Bu_L * _inv_mu0;
    left_flux[5] = 0.0;
    left_flux[6] = uLn * BLt1 - uLt1*BLn;
    left_flux[7] = -uLt2 * BLn + uLn * BLt2;
    left_flux[9] = (EintL) * uLn; 

    right_conserved[0] = rhoR;
    right_conserved[1] = rhoR * uRn;
    right_conserved[2] = rhoR * uRt1;
    right_conserved[3] = rhoR * uRt2;
    right_conserved[4] = EtotR;
    right_conserved[5] = BRn;
    right_conserved[6] = BRt1;
    right_conserved[7] = BRt2;
    right_conserved[8] = BR_scalar;
    right_conserved[9] = EintR;

    right_flux[0] = rhoR * uRn;
    right_flux[1] = right_conserved[1] * uRn + PtotR - BRn * BRn * _inv_mu0;
    right_flux[2] = right_conserved[2] * uRn - BRn * BRt1 * _inv_mu0;
    right_flux[3] = right_conserved[3] * uRn - BRn * BRt2 * _inv_mu0;
    right_flux[4] = (EtotR + PtotR) * uRn - BRn * Bu_R * _inv_mu0;
    right_flux[5] = 0.0;
    right_flux[6] = uRn * BRt1 - uRt1 * BRn;
    right_flux[7] = -uRt2 * BRn + uRn * BRt2;
    right_flux[9] = (EintR) * uRn;

    const double Bx = 0.5 * (BLn + BRn); // normal magnetic field component, shouLnd be the same on both sides but take average to reduce numerical errors

    const double max_u = std::max(uLn, uRn); // (std::abs(uLn), std::abs(uRn));
    const double min_u = std::min(uLn, uRn); // (std::abs(uLn), std::abs(uRn));
    const double S_R =  max_u + cf_max; // right fast magnetosonic
    const double S_L =  min_u - cf_max; // left fast magnetosonic 

  //  double Bx_flux = 0.5 * (BL_scalar + BR_scalar) - 0.5 * _cleaning_speed * (BRn - BLn);
  //  double Bscalar_flux = 0.5 * (_cleaning_speed * _cleaning_speed) * (BLn + BRn) - 0.5 * _cleaning_speed * (BR_scalar - BL_scalar);
//    double Bx_flux = 0.5 * (BL_scalar + BR_scalar) - 0.5 * (BRn - BLn);
  //  double Bscalar_flux = 0.5 *  (BLn + BRn) - 0.5 * (BR_scalar - BL_scalar);
    //double Bx_flux = 0.5 * (left_conserved[8] + right_conserved[8]) 
                   //- 0.5 * _cleaning_speed * (right_conserved[5] - left_conserved[5]);

    //double Bscalar_flux = 0.5 * (_cleaning_speed *_cleaning_speed) * (left_conserved[5] + right_conserved[5]) 
                    //    - 0.5/_cleaning_speed  * (right_conserved[8] - left_conserved[8]);

   //std::cout<< "B flux before HLLD: B_flux = " << Bflux.norm() << std::endl;
   //std::cout<< "Cleaning speed = " << _cleaning_speed << ", cleaning speed factor = " << _cleaning_speed_factor << ", Max velocity = " << _max_velocity << std::endl;

    if (std::isnan(Bflux.x()) || std::isnan(Bflux.y()) || std::isnan(Bflux.z()) || std::isnan(B_scalarflux)) {
        std::cout << "NaN detected in B flux calculation before HLLD solver!" << std::endl;
        std::cout << "B flux before HLLD: B_flux = " << Bflux.norm() << ", Bx_flux " << 0.5 * (BL_scalar + BR_scalar) - 0.5 * _cleaning_speed * (BRn - BLn) << ", By_flux = " << 0.5 * (uLn * BRt1 - uRt1 * BRn) << ", Bz_flux = " << 0.5 * (-uLt2 * BLn + uLn * BLt2) << ", Bscalar_flux = " << 0.5 * (_cleaning_speed_max * _cleaning_speed_max) * (BLn + BRn) - 0.5 * _cleaning_speed_max * (BR_scalar - BL_scalar) << std::endl;
        std::cout << "P flux before HLLD: p_flux = " << pflux.norm() << std::endl;
    };

    if (max_B2 < magnetic_vacuum_threshold){
    //  std::cout<< "Warning: within HD solution as B is negligible!" << std::endl;
      _riemann_solver_hd.solve_for_flux(rhoL, uL, PL, rhoR, uR, PR, mflux, pflux,
                                   Eflux, EintFlux, normal);
      
      /*if (std::abs(Eflux) > current_etot) {
        std::cout << "Warning: within HD solution with Neglibible B field." << std::endl;
        std::cout << "HLLD: Mflux: " << mflux << ", EintFlux: " << EintFlux << ", Eflux: " << Eflux << ", Current Eint: " << current_eint << std::endl;
        std::cout << "HLLD: MASS left_flux: " << left_flux[0] << ", right_flux: " << right_flux[0] << ", right_conserved: " << right_conserved[0] << ", left_conserved: " << left_conserved[0] << std::endl;
        std::cout << "HLLD: ENERGY left_flux: " << left_flux[4] << ", right_flux: " << right_flux[4] << ", right_conserved: " << right_conserved[4] << ", left_conserved: " << left_conserved[4] << std::endl;
        std::cout << "HLLD: S_R: " << S_R << ", S_L: " << S_L << std::endl;
      }*/

      
      // 4. Clear out magnetic field fluxes for the hydro fallback region
      Bflux = CoordinateVector<double>(0.0, 0.0, 0.0);
      B_scalarflux = 0.0;
      if (S_L>0.){
        EintFlux = left_flux[9];
      }
      if (S_R<0.){
        EintFlux = right_flux[9];
      }

     // if (EintFlux < 0. && mflux > 0) {
      //  std::cout << "Warning: EintFlux < 0: " << EintFlux << ", mflux: " << mflux << ", mass right_conserved: " << right_conserved[0] << ", Eint right_conserved:" << right_conserved[9] << std::endl;
     // }

      /*if (std::abs(EintFlux) > current_eint) {
        std::cout << "Warning: within HD solution with Neglibible B field." << std::endl;
        std::cout << "HLLD: Mflux: " << mflux << ", EintFlux: " << EintFlux << ", Eflux: " << Eflux << ", Current Eint: " << current_eint << std::endl;
        std::cout << "HLLD: MASS left_flux: " << left_flux[0] << ", right_flux: " << right_flux[0] << ", right_conserved: " << right_conserved[0] << ", left_conserved: " << left_conserved[0] << std::endl;
        std::cout << "HLLD: ENERGY left_flux: " << left_flux[4] << ", right_flux: " << right_flux[4] << ", right_conserved: " << right_conserved[4] << ", left_conserved: " << left_conserved[4] << std::endl;
        std::cout << "HLLD: INTERNAL_ENERGY left_flux: " << left_flux[9] << ", right_flux: " << right_flux[9] << std::endl;
        std::cout << "HLLD: S_R: " << S_R << ", S_L: " << S_L << std::endl;
      }*/
      
      return;
    };

    if (S_L > 0.0) {
     //std::cout << "HLLD: using left fluxes, S_L = " << S_L << std::endl;

      const CoordinateVector<> Bfluxog = Bflux;
      const CoordinateVector<> pfluxog = pflux;

      mflux = left_flux[0];
      pflux = (left_flux[1] * normal) + (left_flux[2] * t1) + (left_flux[3] * t2); //CoordinateVector<>(left_flux[1], left_flux[2], left_flux[3]);
      Eflux = left_flux[4];
      EintFlux = left_flux[9];
      double Bx_flux = left_conserved[8] + 0.5 * (right_conserved[8] - left_conserved[8]) - 0.5 * _cleaning_speed_max * (right_conserved[5] - left_conserved[5]) ; // HLLD Bx flux
      double By_flux = left_flux[6];
      double Bz_flux = left_flux[7];
      double Bscalar_flux = (_cleaning_speed_max) * (left_conserved[5] + 0.5*(right_conserved[5]-left_conserved[5]) - 0.5 * (right_conserved[8] - left_conserved[8])); // HLLD scalar B flux
      Bflux = (Bx_flux * normal) + (By_flux * t1) + (Bz_flux * t2); //CoordinateVector<>(Bx_flux, left_flux[6], left_flux[7]);
      B_scalarflux = Bscalar_flux;
      //if (right_conserved[0]>1.e-25){
      //  EintFlux = mflux * (right_conserved[9]/(right_conserved[0]*1.e10));
     //} else{
     //   EintFlux = 0.0;
     // }
     // if (EintFlux < 0.&& mflux > 0) {
      //  std::cout << "Warning: EintFlux < 0: " << EintFlux << ", mflux: " << mflux << ", mass right_conserved: " << right_conserved[0] << ", Eint right_conserved:" << right_conserved[9] << std::endl;
     // }
     if (std::isnan(Bflux.x()) || std::isnan(Bflux.y()) || std::isnan(Bflux.z()) || std::isnan(B_scalarflux)) {
        std::cout << "NaN detected in B flux calculation in HLLD S_L>0 solver!" << std::endl;
        std::cout << "B flux before HLLD: B_flux = " << Bfluxog.norm() << "Bx_flux " << Bx_flux << ", By_flux = " << By_flux << ", Bz_flux = " << Bz_flux << ", Bscalar_flux = " << Bscalar_flux << std::endl;
        std::cout << "P flux before HLLD: p_flux = " << pfluxog.norm() << std::endl;
        std::cout<< "B flux after HLLD: B_flux = " << Bflux.norm() << ", Bx_flux " << Bx_flux << ", By_flux = " << By_flux << ", Bz_flux = " << Bz_flux << ", Bscalar_flux = " << Bscalar_flux << std::endl;
        std::cout<< "P flux after HLLD: p_flux = " << pflux.norm() << std::endl;

      
    };
    if (std::abs(Eflux) > current_etot && debug_flag == true) {
        std::cout << "Warning: S_L > 0." << std::endl;
        std::cout << "HLLD: Mflux: " << mflux << ", EintFlux: " << EintFlux << ", Eflux: " << Eflux << ", Current Eint: " << current_eint << std::endl;
        std::cout << "HLLD: MASS left_flux: " << left_flux[0] << ", right_flux: " << right_flux[0] << ", right_conserved: " << right_conserved[0] << ", left_conserved: " << left_conserved[0] << std::endl;
        std::cout << "HLLD: ENERGY left_flux: " << left_flux[4] << ", right_flux: " << right_flux[4] << ", right_conserved: " << right_conserved[4] << ", left_conserved: " << left_conserved[4] << std::endl;
        std::cout << "HLLD: S_R: " << S_R << ", S_L: " << S_L << std::endl;
      };
    if (std::abs(EintFlux) > current_eint && debug_flag == true) {
        std::cout << "Warning: S_L > 0." << std::endl;
        std::cout << "HLLD: Mflux: " << mflux << ", EintFlux: " << EintFlux << ", Eflux: " << Eflux << ", Current Eint: " << current_eint << std::endl;
        std::cout << "HLLD: MASS left_flux: " << left_flux[0] << ", right_flux: " << right_flux[0] << ", right_conserved: " << right_conserved[0] << ", left_conserved: " << left_conserved[0] << std::endl;
        std::cout << "HLLD: ENERGY left_flux: " << left_flux[4] << ", right_flux: " << right_flux[4] << ", right_conserved: " << right_conserved[4] << ", left_conserved: " << left_conserved[4] << std::endl;
        std::cout << "HLLD: INTERNAL_ENERGY left_flux: " << left_flux[9] << ", right_flux: " << right_flux[9] << std::endl;
        std::cout << "HLLD: S_R: " << S_R << ", S_L: " << S_L << std::endl;
      }

      return;
    } 
    
    if (S_R < 0.0) {
   //  std::cout << "HLLD: using right fluxes, S_R = " << S_R << std::endl;

      const  CoordinateVector<> Bfluxog = Bflux;
      const  CoordinateVector<> pfluxog = pflux;

      mflux = right_flux[0];
      pflux = (right_flux[1] * normal) + (right_flux[2] * t1) + (right_flux[3] * t2); //CoordinateVector<>(right_flux[1], right_flux[2], right_flux[3]);
      Eflux = right_flux[4];
      EintFlux = right_flux[9];
      double Bx_flux = left_conserved[8] + 0.5 * (right_conserved[8] - left_conserved[8]) - 0.5 * _cleaning_speed_max * (right_conserved[5] - left_conserved[5]) ; // HLLD Bx flux
      double By_flux = right_flux[6];
      double Bz_flux = right_flux[7];
      double Bscalar_flux = (_cleaning_speed_max) * (left_conserved[5] + 0.5*(right_conserved[5]-left_conserved[5]) - 0.5 * (right_conserved[8] - left_conserved[8])); //(_cleaning_speed_max * _cleaning_speed_max) * (left_conserved[5] + 0.5*(right_conserved[5]-left_conserved[5]) - 0.5/_cleaning_speed_max * (right_conserved[8] - left_conserved[8])); // HLLD scalar B flux
      Bflux = (Bx_flux * normal) + (By_flux * t1) + (Bz_flux * t2); //CoordinateVector<>(Bx_flux, right_flux[6], right_flux[7]);
      B_scalarflux = Bscalar_flux;

      

     if (std::isnan(Bflux.x()) || std::isnan(Bflux.y()) || std::isnan(Bflux.z()) || std::isnan(B_scalarflux) || debug_flag == true) {
        std::cout << "NaN detected in B flux calculation in HLLD S_R<0 solver!" << std::endl;
        std::cout << "B flux before HLLD: B_flux = " << Bfluxog.norm() << "Bx_flux " << Bx_flux << ", By_flux = " << By_flux << ", Bz_flux = " << Bz_flux << ", Bscalar_flux = " << Bscalar_flux << std::endl;
        std::cout << "P flux before HLLD: p_flux = " << pfluxog.norm() << std::endl;
        std::cout<< "B flux after HLLD: B_flux = " << Bflux.norm() << ", Bx_flux " << Bx_flux << ", By_flux = " << By_flux << ", Bz_flux = " << Bz_flux << ", Bscalar_flux = " << Bscalar_flux << std::endl;
        std::cout<< "P flux after HLLD: p_flux = " << pflux.norm() << std::endl;

      
    };
    if (std::abs(Eflux) > current_etot && debug_flag == true) {
        std::cout << "Warning: S_R < 0." << std::endl;
        std::cout << "HLLD: Mflux: " << mflux << ", EintFlux: " << EintFlux << ", Eflux: " << Eflux << ", Current Eint: " << current_eint << std::endl;
        std::cout << "HLLD: MASS left_flux: " << left_flux[0] << ", right_flux: " << right_flux[0] << ", right_conserved: " << right_conserved[0] << ", left_conserved: " << left_conserved[0] << std::endl;
        std::cout << "HLLD: ENERGY left_flux: " << left_flux[4] << ", right_flux: " << right_flux[4] << ", right_conserved: " << right_conserved[4] << ", left_conserved: " << left_conserved[4] << std::endl;
        std::cout << "HLLD: S_R: " << S_R << ", S_L: " << S_L << std::endl;
      };

    if (std::abs(EintFlux) > current_eint && debug_flag == true) {
        std::cout << "Warning: S_R < 0." << std::endl;
        std::cout << "HLLD: Mflux: " << mflux << ", EintFlux: " << EintFlux << ", Eflux: " << Eflux << ", Current Eint: " << current_eint << std::endl;
        std::cout << "HLLD: MASS left_flux: " << left_flux[0] << ", right_flux: " << right_flux[0] << ", right_conserved: " << right_conserved[0] << ", left_conserved: " << left_conserved[0] << std::endl;
        std::cout << "HLLD: ENERGY left_flux: " << left_flux[4] << ", right_flux: " << right_flux[4] << ", right_conserved: " << right_conserved[4] << ", left_conserved: " << left_conserved[4] << std::endl;
        std::cout << "HLLD: INTERNAL_ENERGY left_flux: " << left_flux[9] << ", right_flux: " << right_flux[9] << std::endl;
        std::cout << "HLLD: S_R: " << S_R << ", S_L: " << S_L << std::endl;
      }
      return;
    }
    

    double S_M = ((S_R-uRn)*rhoR*uRn-(S_L-uLn)*rhoL*uLn-(PtotR-PtotL))/((S_R - uRn)*rhoR-(S_L - uLn)*rhoL); // middle speed (discontinuity)
    if (!std::isfinite(S_M)) {
      std::cout << "Warning: S_M is not finite, which may cause numerical issues. S_L: " << S_L << ", S_R: " << S_R << ", uLn: " << uLn << ", uRn: " << uRn << ", rhoL: " << rhoL << ", rhoR: " << rhoR << ", PtotL: " << PtotL << ", PtotR: " << PtotR << std::endl;
      S_M = 0.5 * (uLn + uRn);
    }
    double denL = S_L - S_M;
    if (std::abs(denL) < 1e-12)
    {
     //std::cout << "Warning: S_L - S_M is very small, which may cause numerical issues. S_L: " << S_L << ", S_M: " << S_M << std::endl;
      double denL = std::copysign(1e-12, denL);
    }
    double denR = S_R - S_M;
    if (std::abs(denR) < 1e-12)    {
     //std::cout << "Warning: S_R - S_M is very small, which may cause numerical issues. S_R: " << S_R << ", S_M: " << S_M << std::endl;
      double denR = std::copysign(1e-12, denR);
    }

    double rhoL_star = rhoL * (S_L - uLn) / denL;
    double rhoR_star = rhoR * (S_R - uRn) / denR;
    rhoL_star = std::max(rhoL_star, 1e-25); 
    rhoR_star = std::max(rhoR_star, 1e-25); // 

    const double S_L_star = S_M - std::abs(Bnorm) * _inv_sqrt_mu0 / std::sqrt(rhoL_star); // left alfven speed
    const double S_R_star = S_M + std::abs(Bnorm) * _inv_sqrt_mu0 / std::sqrt(rhoR_star); // right alfven speed
    double p_den = (S_R - uRn) * rhoR - (S_L - uLn) * rhoL;
    if (std::abs(p_den) < 1e-15) {
      p_den = std::copysign(1e-15, p_den);
    }

    const double Ptot_star = ((S_R-uRn)*rhoR*PtotL - (S_L-uLn)*rhoL*PtotR + (rhoL*rhoR*(S_R-uRn)*(S_L-uLn)*(uRn-uLn)))/p_den;

    const double samL = uLn - cfL;
    const double sapL = uLn + cfL; 
    const double samR = uRn - cfR; 
    const double sapR = uRn + cfR;

    double uyLs, uzLs, ByLs, BzLs;
    double uyRs, uzRs, ByRs, BzRs;
    double BBLs, BBRs, vvLs, vvRs;

    if ((std::abs(S_M-uLn) < BFLOAT_EPSILON) and (std::abs(BLt1) < BFLOAT_EPSILON) and (std::abs(BLt2) < BFLOAT_EPSILON) and (Bx*Bx*_inv_mu0 >= _gamma * PL) and ((std::abs(S_L - samL) <= BFLOAT_EPSILON) or (std::abs(S_L-sapL) <= BFLOAT_EPSILON))) { // check for left star region being very close to left state and magnetic field being very weak, in which case we can just use the left fluxes to avoid numerical issues with calcuLnating the star region
       uyLs = uLt1;
       uzLs = uLt2;
       ByLs = BLt1;
       BzLs = BLt2;
    
    } else {
       double den = rhoL * (S_L-uLn) * (S_L-S_M) - Bx*Bx*_inv_mu0;
      // if (std::abs(den) < 1e-15) {
       //  den = std::copysign(1e-15, den);
      // }
       vvLs = (S_M - uLn) / den;
       BBLs = (rhoL * (S_L-uLn)*(S_L-uLn)-Bx*Bx*_inv_mu0) / den;
       uyLs = uLt1 - Bx * BLt1 * vvLs;
       uzLs = uLt2 - Bx * BLt2 * vvLs;
       ByLs = BLt1 * BBLs;
       BzLs = BLt2 * BBLs;
    }

    if (std::abs(S_M-uRn) < BFLOAT_EPSILON and std::abs(BRt1) < BFLOAT_EPSILON and std::abs(BRt2) < BFLOAT_EPSILON and Bx*Bx*_inv_mu0 >= _gamma * PR and (std::abs(S_R - samR) <= BFLOAT_EPSILON or std::abs(S_R-sapR) <= BFLOAT_EPSILON)) { // check for right star region being very close to right state and magnetic field being very weak, in which case we can just use the right fluxes to avoid numerical issues with calcuLnating the star region
       uyRs = uRt1;
       uzRs = uRt2;
       ByRs = BRt1;
       BzRs = BRt2;
    
    } else {
       double den = rhoR * (S_R-uRn) * (S_R-S_M) - Bx*Bx*_inv_mu0;
      // if (std::abs(den) < 1e-15) {
      //   den = std::copysign(1e-15, den);
     //  }
       vvRs = (S_M - uRn) / den;
       BBRs = (rhoR * (S_R-uRn)*(S_R-uRn)-Bx*Bx*_inv_mu0) / den;
       uyRs = uRt1 - Bx * BRt1 * vvRs;
       uzRs = uRt2 - Bx * BRt2 * vvRs;
       ByRs = BRt1 * BBRs;
       BzRs = BRt2 * BBRs;
    }

    const double Bu_Ls = S_M * Bx + uyLs * ByLs + uzLs * BzLs;
    const double Bu_Rs = S_M * Bx + uyRs * ByRs + uzRs * BzRs;

    const double EtotLS = ((S_L - uLn) * EtotL - PtotL * uLn + Ptot_star * S_M + Bx*(Bu_L - Bu_Ls)*_inv_mu0) / (S_L - S_M);
    const double EtotRS = ((S_R - uRn) * EtotR - PtotR * uRn + Ptot_star * S_M + Bx*(Bu_R - Bu_Rs)*_inv_mu0) / (S_R - S_M);

//    const double c_sound = std::sqrt(_gamma * PL/rhoL);
 //   double mach = std::abs(uLn)/c_sound;


    if (S_L <= 0 && S_L_star >=0){
    // std::cout << "HLLD: using left star fluxes, S_L = " << S_L << ", S_L_star = " << S_L_star << std::endl;
      star_conserved[0] = rhoL_star;
      star_conserved[1] = rhoL_star * S_M;
      star_conserved[2] = rhoL_star * uyLs;
      star_conserved[3] = rhoL_star * uzLs;
      star_conserved[4] = EtotLS;
      star_conserved[5] = Bx;
      star_conserved[6] = ByLs;
      star_conserved[7] = BzLs;
      star_conserved[8] = BL_scalar;

      double proposed_Eint_star = EintL * (rhoL_star / rhoL);

      double e_kin_star = 0.5 * rhoL_star * (S_M*S_M + uyLs*uyLs + uzLs*uzLs);
      double e_mag_star = (Bx*Bx + ByLs*ByLs + BzLs*BzLs) * 0.5 * _inv_mu0;

      double max_allowable_Eint_star = EtotLS - e_kin_star - e_mag_star;

      if (proposed_Eint_star >= max_allowable_Eint_star || proposed_Eint_star <= 0.0) {
          proposed_Eint_star = max_allowable_Eint_star;
      }
      star_conserved[9] = proposed_Eint_star;

      const CoordinateVector<> Bfluxog = Bflux;
      const CoordinateVector<> pfluxog = pflux;
    
      mflux = left_flux[0] + S_L * (star_conserved[0] - left_conserved[0]);
      double pflux_x = left_flux[1] + S_L * (star_conserved[1] - left_conserved[1]);
      double pflux_y = left_flux[2] + S_L * (star_conserved[2] - left_conserved[2]);
      double pflux_z = left_flux[3] + S_L * (star_conserved[3] - left_conserved[3]);
      pflux = (pflux_x * normal) + (pflux_y * t1) + (pflux_z * t2); //CoordinateVector<>(pflux_x, pflux_y, pflux_z);
      Eflux = left_flux[4] + S_L * (star_conserved[4] - left_conserved[4]);
      EintFlux = left_flux[9] + S_L * (star_conserved[9] - left_conserved[9]);
      double Bx_flux = left_conserved[8] + 0.5 * (right_conserved[8] - left_conserved[8]) - 0.5 * _cleaning_speed_max * (right_conserved[5] - left_conserved[5]) ; // HLLD Bx flux
      double Bscalar_flux = (_cleaning_speed_max ) * (left_conserved[5] + 0.5*(right_conserved[5]-left_conserved[5]) - 0.5 * (right_conserved[8] - left_conserved[8])); // HLLD scalar B flux
      double By_flux = left_flux[6] + S_L * (star_conserved[6] - left_conserved[6]);
      double Bz_flux = left_flux[7] + S_L * (star_conserved[7] - left_conserved[7]);
      Bflux = (Bx_flux * normal) + (By_flux * t1) + (Bz_flux * t2); // HLLD B flux
      B_scalarflux = Bscalar_flux;
      
      
    
     if (std::isnan(Bflux.x()) || std::isnan(Bflux.y()) || std::isnan(Bflux.z()) || std::isnan(B_scalarflux) ) {
        std::cout << "NaN detected in B flux calculation in HLLD S_L<=0, S_L_star>=0 solver!" << std::endl;
        std::cout << "B flux before HLLD: B_flux = " << Bfluxog.norm() << ", Bx_flux " << Bx_flux << ", By_flux = " << By_flux << ", Bz_flux = " << Bz_flux << ", Bscalar_flux = " << Bscalar_flux << std::endl;
        std::cout << "P flux before HLLD: p_flux = " << pfluxog.norm() << ", px_flux " << pfluxog.x() << ", py_flux " << pfluxog.y() << ", pz_flux " << pfluxog.z() << std::endl;
        std::cout<< "B flux after HLLD: B_flux = " << Bflux.norm() << ", Bx_flux " << Bx_flux << ", By_flux = " << By_flux << ", Bz_flux = " << Bz_flux << ", Bscalar_flux = " << Bscalar_flux << std::endl;
        std::cout<< "P flux after HLLD: p_flux = " << pflux.norm() << ", px_flux " << pflux_x << ", py_flux " << pflux_y << ", pz_flux " << pflux_z << std::endl;
        std::cout<< "variables used Bflux: star_conserved[6] = " << star_conserved[6] << ", star_conserved[7] = " << star_conserved[7] << ", left_conserved[6] = " << left_conserved[6] << ", left_conserved[7] = " << left_conserved[7] << ", left_flux[6] = " << left_flux[6] << ", left_flux[7] = " << left_flux[7] << std::endl;
        std::cout<< "variables used Pflux: star_conserved[1] = " <<star_conserved[1] << ", star_conserved[2] = " << star_conserved[2] << ", star_conserved[3] = " << star_conserved[3] << ", left_conserved[1] = " << left_conserved[1] << ", left_conserved[2] = " << left_conserved[2] << ", left_conserved[3] = " << left_conserved[3] << ", left_flux[1] = " << left_flux[1] << ", left_flux[2] = " << left_flux[2] << ", left_flux[3] = " << left_flux[3] << std::endl;
        std::cout<< "primary vars: rhoL_star = " << rhoL_star << ", S_M = " << S_M << ", uyLs = " << uyLs << ", uzLs = " << uzLs << ", ByLs = " << ByLs << ", BzLs = " << BzLs << std::endl;
        
    };
    if (std::abs(Eflux) > current_etot && debug_flag == true) {
        std::cout << "Warning: S_L <= 0 && S_L_star >=0." << std::endl;
        std::cout << "HLLD: Mflux: " << mflux << ", EintFlux: " << EintFlux << ", Eflux: " << Eflux << ", Current Eint: " << current_eint << std::endl;
        std::cout << "HLLD: MASS left_flux: " << left_flux[0] << ", right_flux: " << right_flux[0] << ", right_conserved: " << right_conserved[0] << ", left_conserved: " << left_conserved[0] << std::endl;
        std::cout << "HLLD: ENERGY left_flux: " << left_flux[4] << ", right_flux: " << right_flux[4] << ", right_conserved: " << right_conserved[4] << ", left_conserved: " << left_conserved[4] << std::endl;
        std::cout << "HLLD: S_R: " << S_R << ", S_L: " << S_L << std::endl;
      };
    if (std::abs(EintFlux) > current_eint && debug_flag == true) {
        std::cout << "Warning: S_L <= 0 && S_L_star >=0." << std::endl;
        std::cout << "HLLD: Mflux: " << mflux << ", EintFlux: " << EintFlux << ", Eflux: " << Eflux << ", Current Eint: " << current_eint << std::endl;
        std::cout << "HLLD: MASS left_flux: " << left_flux[0] << ", right_flux: " << right_flux[0] << ", right_conserved: " << right_conserved[0] << ", left_conserved: " << left_conserved[0] << std::endl;
        std::cout << "HLLD: ENERGY left_flux: " << left_flux[4] << ", right_flux: " << right_flux[4] << ", right_conserved: " << right_conserved[4] << ", left_conserved: " << left_conserved[4] << std::endl;
        std::cout << "HLLD: INTERNAL_ENERGY left_flux: " << left_flux[9] << ", right_flux: " << right_flux[9] << std::endl;
        std::cout << "HLLD: S_R: " << S_R << ", S_L: " << S_L << std::endl;
      }
      return;

    }
    if (S_R_star <= 0 && S_R >= 0) {
    // std::cout << "HLLD: using right star fluxes, S_R = " << S_R << ", S_R_star = " << S_R_star << std::endl;
      star_conserved[0] = rhoR_star;
      star_conserved[1] = rhoR_star * S_M;
      star_conserved[2] = rhoR_star * uyRs;
      star_conserved[3] = rhoR_star * uzRs;
      star_conserved[4] = EtotRS;
      star_conserved[5] = Bx;
      star_conserved[6] = ByRs;
      star_conserved[7] = BzRs;
      star_conserved[8] = BR_scalar;
      double proposed_Eint_star = EintR * (rhoR_star / rhoR);

      double e_kin_star = 0.5 * rhoR_star * (S_M*S_M + uyRs*uyRs + uzRs*uzRs);
      double e_mag_star = (Bx*Bx + ByRs*ByRs + BzRs*BzRs) * 0.5 * _inv_mu0;

      double max_allowable_Eint_star = EtotRS - e_kin_star - e_mag_star;

      if (proposed_Eint_star >= max_allowable_Eint_star || proposed_Eint_star <= 0.0) {
          proposed_Eint_star = max_allowable_Eint_star;
      }
      star_conserved[9] = proposed_Eint_star;

      const CoordinateVector<> Bfluxog = Bflux;
      const CoordinateVector<> pfluxog = pflux;

      mflux = right_flux[0] + S_R * (star_conserved[0] - right_conserved[0]);
      double pflux_x = right_flux[1] + S_R * (star_conserved[1] - right_conserved[1]);
      double pflux_y = right_flux[2] + S_R * (star_conserved[2] - right_conserved[2]);
      double pflux_z = right_flux[3] + S_R * (star_conserved[3] - right_conserved[3]);
      pflux = (pflux_x * normal) + (pflux_y * t1) + (pflux_z * t2); //CoordinateVector<>(pflux_x, pflux_y, pflux_z);
      Eflux = right_flux[4] + S_R * (star_conserved[4] - right_conserved[4]);
      EintFlux = right_flux[9]+ S_R * (star_conserved[9] - right_conserved[9]);
      double Bx_flux = left_conserved[8] + 0.5 * (right_conserved[8] - left_conserved[8]) - 0.5 * _cleaning_speed_max * (right_conserved[5] - left_conserved[5]) ; // HLLD Bx flux
      double By_flux = right_flux[6] + S_R * (star_conserved[6] - right_conserved[6]);
      double Bz_flux = right_flux[7] + S_R * (star_conserved[7] - right_conserved[7]);
      double Bscalar_flux = (_cleaning_speed_max) * (left_conserved[5] + 0.5*(right_conserved[5]-left_conserved[5]) - 0.5 * ( right_conserved[8] - left_conserved[8])); // HLLD scalar B flux
      Bflux = (Bx_flux * normal) + (By_flux * t1) + (Bz_flux * t2); // HLLD B flux
      B_scalarflux = Bscalar_flux;

      
      
      
     if (std::isnan(Bflux.x()) || std::isnan(Bflux.y()) || std::isnan(Bflux.z()) || std::isnan(B_scalarflux)) {
        std::cout << "NaN detected in B flux calculation in HLLD S_R_star<=0, S_R>=0 solver!" << std::endl;
        std::cout << "B flux before HLLD: B_flux = " << Bfluxog.norm() << "Bx_flux " << Bx_flux << ", By_flux = " << By_flux << ", Bz_flux = " << Bz_flux << ", Bscalar_flux = " << Bscalar_flux << std::endl;
        std::cout << "P flux before HLLD: p_flux = " << pfluxog.norm() << std::endl;
        std::cout<< "B flux after HLLD: B_flux = " << Bflux.norm() << ", Bx_flux " << Bx_flux << ", By_flux = " << By_flux << ", Bz_flux = " << Bz_flux << ", Bscalar_flux = " << Bscalar_flux << std::endl;
        std::cout<< "P flux after HLLD: p_flux = " << pflux.norm() << std::endl;

      
    };
    if (std::abs(Eflux) > current_etot && debug_flag == true) {
        std::cout << "Warning: S_R_star <= 0 && S_R >= 0." << std::endl;
        std::cout << "HLLD: Mflux: " << mflux << ", EintFlux: " << EintFlux << ", Eflux: " << Eflux << ", Current Eint: " << current_eint << std::endl;
        std::cout << "HLLD: MASS left_flux: " << left_flux[0] << ", right_flux: " << right_flux[0] << ", right_conserved: " << right_conserved[0] << ", left_conserved: " << left_conserved[0] << std::endl;
        std::cout << "HLLD: ENERGY left_flux: " << left_flux[4] << ", right_flux: " << right_flux[4] << ", right_conserved: " << right_conserved[4] << ", left_conserved: " << left_conserved[4] << std::endl;
        std::cout << "HLLD: S_R: " << S_R << ", S_L: " << S_L << std::endl;
      };
    if (std::abs(EintFlux) > current_eint && debug_flag == true) {
        std::cout << "Warning: S_R_star <= 0 && S_R >= 0." << std::endl;
        std::cout << "HLLD: Mflux: " << mflux << ", EintFlux: " << EintFlux << ", Eflux: " << Eflux << ", Current Eint: " << current_eint << std::endl;
        std::cout << "HLLD: MASS left_flux: " << left_flux[0] << ", right_flux: " << right_flux[0] << ", right_conserved: " << right_conserved[0] << ", left_conserved: " << left_conserved[0] << std::endl;
        std::cout << "HLLD: ENERGY left_flux: " << left_flux[4] << ", right_flux: " << right_flux[4] << ", right_conserved: " << right_conserved[4] << ", left_conserved: " << left_conserved[4] << std::endl;
        std::cout << "HLLD: INTERNAL_ENERGY left_flux: " << left_flux[9] << ", right_flux: " << right_flux[9] << std::endl;
        std::cout << "HLLD: S_R: " << S_R << ", S_L: " << S_L << std::endl;
      }

      return;
    } 
    

    const double rho_star_average = std::sqrt(rhoL_star) + std::sqrt(rhoR_star);
    const double uy_starstar = (std::sqrt(rhoL_star) * uyLs + std::sqrt(rhoR_star) * uyRs + (ByRs - ByLs) * std::copysign(1.0, Bx) * _inv_sqrt_mu0) / rho_star_average;
    const double uz_starstar = (std::sqrt(rhoL_star) * uzLs + std::sqrt(rhoR_star) * uzRs + (BzRs - BzLs) * std::copysign(1.0, Bx) * _inv_sqrt_mu0) / rho_star_average;
    const double By_starstar = (std::sqrt(rhoL_star) * ByRs + std::sqrt(rhoR_star) * ByLs + std::sqrt(rhoL_star * rhoR_star) * (uyRs - uyLs) * std::copysign(1.0, Bx) / _inv_sqrt_mu0) / rho_star_average;
    const double Bz_starstar = (std::sqrt(rhoL_star) * BzRs + std::sqrt(rhoR_star) * BzLs + std::sqrt(rhoL_star * rhoR_star) * (uzRs - uzLs) * std::copysign(1.0, Bx) / _inv_sqrt_mu0) / rho_star_average;
    const double Bu_starstar = S_M * Bx + uy_starstar * By_starstar + uz_starstar * Bz_starstar;
    const double EtotL_starstar = EtotLS + std::sqrt(rhoL_star) * (Bu_Ls - Bu_starstar) * std::copysign(1.0, Bx)*_inv_sqrt_mu0;
    const double EtotR_starstar = EtotRS + std::sqrt(rhoR_star) * (Bu_Rs - Bu_starstar) * std::copysign(1.0, Bx)*_inv_sqrt_mu0;

   //std::cout<< " used params above HLLD: rhoL_star = " << rhoL_star << ", uyLs = " << uyLs << ", uzLs = " << uzLs << ", ByLs = " << ByLs << ", BzLs = " << BzLs << ", rho_star_average = " << rho_star_average << std::endl;
    //std::cout<< " used params above HLLD: rhoR_star = " << rhoR_star << ", uyRs = " << uyRs << ", uzRs = " << uzRs << ", ByRs = " << ByRs << ", BzRs = " << BzRs << ", rho_star_average = " << rho_star_average << std::endl;

    if (S_L_star <= 0. && S_M >= 0.0) {
     // std::cout << "HLLD: using left star-star fluxes, S_L_star = " << S_L_star << ", S_M = " << S_M << std::endl;
      star_conserved[0] = rhoL_star; 
      star_conserved[1] = rhoL_star * S_M;
      star_conserved[2] = rhoL_star * uyLs;
      star_conserved[3] = rhoL_star * uzLs;
      star_conserved[4] = EtotLS;
      star_conserved[5] = Bx;
      star_conserved[6] = ByLs;
      star_conserved[7] = BzLs;
      star_conserved[8] = BL_scalar;
      double proposed_Eint_star = EintL * (rhoL_star / rhoL);

      double e_kin_star = 0.5 * rhoL_star * (S_M*S_M + uyLs*uyLs + uzLs*uzLs);
      double e_mag_star = (Bx*Bx + ByLs*ByLs + BzLs*BzLs) * 0.5 * _inv_mu0;

      double max_allowable_Eint_star = EtotLS - e_kin_star - e_mag_star;

      if (proposed_Eint_star >= max_allowable_Eint_star || proposed_Eint_star <= 0.0) {
          proposed_Eint_star = max_allowable_Eint_star;
      }
      star_conserved[9] = proposed_Eint_star;
      
      starstar_conserved[0] = rhoL_star;
      starstar_conserved[1] = rhoL_star * S_M;
      starstar_conserved[2] = rhoL_star * uy_starstar;
      starstar_conserved[3] = rhoL_star * uz_starstar;
      starstar_conserved[4] = EtotL_starstar;
      starstar_conserved[5] = Bx;
      starstar_conserved[6] = By_starstar;
      starstar_conserved[7] = Bz_starstar;
      starstar_conserved[8] = BL_scalar;
      double proposed_Eint_starstar = EintL * (rhoL_star / rhoL);

      double e_kin_starstar = 0.5 * rhoL_star * (S_M*S_M + uy_starstar*uy_starstar + uz_starstar*uz_starstar);
      double e_mag_starstar = (Bx*Bx + By_starstar*By_starstar + Bz_starstar*Bz_starstar) * 0.5 * _inv_mu0;

      double max_allowable_Eint_starstar = EtotL_starstar - e_kin_starstar - e_mag_starstar;

      if (proposed_Eint_starstar >= max_allowable_Eint_starstar || proposed_Eint_starstar <= 0.0) {
          proposed_Eint_starstar = max_allowable_Eint_starstar;
      }
      starstar_conserved[9] = proposed_Eint_starstar;
      //std::cout << "star stars used: 5 = " << starstar_conserved[5] << ", 6 = " << starstar_conserved[6] << ", 7 = " << starstar_conserved[7] << std::endl; 
      //std::cout << "star used: 5 = " << star_conserved[5] << ", 6 = " << star_conserved[6] << ", 7 = " << star_conserved[7] << std::endl;
      //std::cout << "left used: 5 = " << left_conserved[5] << ", 6 = " << left_conserved[6] << ", 7 = " << left_conserved[7] << std::endl;
      //std::cout << "right used: 5 = " << right_conserved[5] << ", 6 = " << right_conserved[6] << ", 7 = " << right_conserved[7] << std::endl;
      //std::cout << "inv mu0 = " << _inv_sqrt_mu0 << std::endl;
      const CoordinateVector<> Bfluxog = Bflux;
      const CoordinateVector<> pfluxog = pflux;
      mflux = left_flux[0] + S_L_star * starstar_conserved[0] - (S_L_star - S_L) * star_conserved[0] - S_L * left_conserved[0];
      double pflux_x = left_flux[1] + S_L_star * starstar_conserved[1] - (S_L_star - S_L) * star_conserved[1] - S_L * left_conserved[1];
      double pflux_y = left_flux[2] + S_L_star * starstar_conserved[2] - (S_L_star - S_L) * star_conserved[2] - S_L * left_conserved[2];
      double pflux_z = left_flux[3] + S_L_star * starstar_conserved[3] - (S_L_star - S_L) * star_conserved[3] - S_L * left_conserved[3];
      pflux = (pflux_x * normal) + (pflux_y * t1) + (pflux_z * t2); //CoordinateVector<>(pflux_x, pflux_y, pflux_z);
      Eflux = left_flux[4] + S_L_star * starstar_conserved[4] - (S_L_star - S_L) * star_conserved[4] - S_L * left_conserved[4];
      EintFlux = left_flux[9] + S_L_star * starstar_conserved[9] - (S_L_star - S_L) * star_conserved[9] - S_L * left_conserved[9];
      double Bx_flux = left_conserved[8] + 0.5 * (right_conserved[8] - left_conserved[8]) - 0.5 * _cleaning_speed_max * (right_conserved[5] - left_conserved[5]) ; // HLLD Bx flux
      double Bscalar_flux = (_cleaning_speed_max) * (left_conserved[5] + 0.5*(right_conserved[5]-left_conserved[5]) - 0.5 * (right_conserved[8] - left_conserved[8])); // HLLD scalar B flux
      double By_flux = left_flux[6] + S_L_star * starstar_conserved[6] - (S_L_star - S_L) * star_conserved[6] - S_L * left_conserved[6];
      double Bz_flux = left_flux[7] + S_L_star * starstar_conserved[7] - (S_L_star - S_L) * star_conserved[7] - S_L * left_conserved[7];
      Bflux = (Bx_flux * normal) + (By_flux * t1) + (Bz_flux * t2); //CoordinateVector<>(Bx_flux, By_flux, Bz_flux); // HLLD B flux
      B_scalarflux = Bscalar_flux;
      
    
      
    
      if (std::isnan(Bflux.x()) || std::isnan(Bflux.y()) || std::isnan(Bflux.z()) || std::isnan(B_scalarflux)) {
        std::cout << "NaN detected in B flux calculation in HLLD S_L_star<=0, S_M>0 solver!" << std::endl;
        std::cout << "B flux before HLLD: B_flux = " << Bfluxog.norm() << "Bx_flux " << Bx_flux << ", By_flux = " << By_flux << ", Bz_flux = " << Bz_flux << ", Bscalar_flux = " << Bscalar_flux << std::endl;
        std::cout << "P flux before HLLD: p_flux = " << pfluxog.norm() << std::endl;
        std::cout<< "B flux after HLLD: B_flux = " << Bflux.norm() << ", Bx_flux " << Bx_flux << ", By_flux = " << By_flux << ", Bz_flux = " << Bz_flux << ", Bscalar_flux = " << Bscalar_flux << std::endl;
        std::cout<< "P flux after HLLD: p_flux = " << pflux.norm() << std::endl;

      
      };
      if (std::abs(Eflux) > current_etot && debug_flag == true) {
        std::cout << "Warning: S_L_star <= 0. && S_M >= 0.0." << std::endl;
        std::cout << "HLLD: Mflux: " << mflux << ", EintFlux: " << EintFlux << ", Eflux: " << Eflux << ", Current Eint: " << current_eint << std::endl;
        std::cout << "HLLD: MASS left_flux: " << left_flux[0] << ", right_flux: " << right_flux[0] << ", right_conserved: " << right_conserved[0] << ", left_conserved: " << left_conserved[0] << std::endl;
        std::cout << "HLLD: ENERGY left_flux: " << left_flux[4] << ", right_flux: " << right_flux[4] << ", right_conserved: " << right_conserved[4] << ", left_conserved: " << left_conserved[4] << std::endl;
        std::cout << "HLLD: S_R: " << S_R << ", S_L: " << S_L << std::endl;
      };
    if (std::abs(EintFlux) > current_eint && debug_flag == true) {
        std::cout << "Warning: S_L_star <= 0. && S_M >= 0.0." << std::endl;
        std::cout << "HLLD: Mflux: " << mflux << ", EintFlux: " << EintFlux << ", Eflux: " << Eflux << ", Current Eint: " << current_eint << std::endl;
        std::cout << "HLLD: MASS left_flux: " << left_flux[0] << ", right_flux: " << right_flux[0] << ", right_conserved: " << right_conserved[0] << ", left_conserved: " << left_conserved[0] << std::endl;
        std::cout << "HLLD: ENERGY left_flux: " << left_flux[4] << ", right_flux: " << right_flux[4] << ", right_conserved: " << right_conserved[4] << ", left_conserved: " << left_conserved[4] << std::endl;
        std::cout << "HLLD: INTERNAL_ENERGY left_flux: " << left_flux[9] << ", right_flux: " << right_flux[9] << std::endl;
        std::cout << "HLLD: S_R: " << S_R << ", S_L: " << S_L << std::endl;
      }
      return;
    }
    if (S_M <= 0.0 && S_R_star >= 0.0) {
     // std::cout << "HLLD: using right star-star fluxes, S_M = " << S_M << ", S_R_star = " << S_R_star << std::endl;
      star_conserved[0] = rhoR_star; 
      star_conserved[1] = rhoR_star * S_M;
      star_conserved[2] = rhoR_star * uyRs;
      star_conserved[3] = rhoR_star * uzRs;
      star_conserved[4] = EtotRS;
      star_conserved[5] = Bx;
      star_conserved[6] = ByRs;
      star_conserved[7] = BzRs;
      star_conserved[8] = BR_scalar;
      double proposed_Eint_star = EintR * (rhoR_star / rhoR);

      double e_kin_star = 0.5 * rhoR_star * (S_M*S_M + uyRs*uyRs + uzRs*uzRs);
      double e_mag_star = (Bx*Bx + ByRs*ByRs + BzRs*BzRs) * 0.5 * _inv_mu0;

      double max_allowable_Eint_star = EtotRS - e_kin_star - e_mag_star;

      if (proposed_Eint_star >= max_allowable_Eint_star || proposed_Eint_star <= 0.0) {
          proposed_Eint_star = max_allowable_Eint_star;
      }
      star_conserved[9] = proposed_Eint_star;

      starstar_conserved[0] = rhoR_star;
      starstar_conserved[1] = rhoR_star * S_M;
      starstar_conserved[2] = rhoR_star * uy_starstar;
      starstar_conserved[3] = rhoR_star * uz_starstar;
      starstar_conserved[4] = EtotR_starstar;
      starstar_conserved[5] = Bx;
      starstar_conserved[6] = By_starstar;
      starstar_conserved[7] = Bz_starstar;
      starstar_conserved[8] = BR_scalar;
      double proposed_Eint_starstar = EintR * (rhoR_star / rhoR);

      double e_kin_starstar = 0.5 * rhoR_star * (S_M*S_M + uy_starstar*uy_starstar + uz_starstar*uz_starstar);
      double e_mag_starstar = (Bx*Bx + By_starstar*By_starstar + Bz_starstar*Bz_starstar) * 0.5 * _inv_mu0;

      double max_allowable_Eint_starstar = EtotR_starstar - e_kin_starstar - e_mag_starstar;

      if (proposed_Eint_starstar >= max_allowable_Eint_starstar || proposed_Eint_starstar <= 0.0) {
          proposed_Eint_starstar = max_allowable_Eint_starstar;
      }
      starstar_conserved[9] = proposed_Eint_starstar;

      const CoordinateVector<> Bfluxog = Bflux;
      const CoordinateVector<> pfluxog = pflux;

      mflux = right_flux[0] + S_R_star * starstar_conserved[0] - (S_R_star - S_R) * star_conserved[0] - S_R * right_conserved[0];
      double pflux_x = (right_flux[1] + S_R_star * starstar_conserved[1] - (S_R_star - S_R) * star_conserved[1] - S_R * right_conserved[1]);
      double pflux_y = (right_flux[2] + S_R_star * starstar_conserved[2] - (S_R_star - S_R) * star_conserved[2] - S_R * right_conserved[2]);
      double pflux_z = (right_flux[3] + S_R_star * starstar_conserved[3] - (S_R_star - S_R) * star_conserved[3] - S_R * right_conserved[3]);
      pflux = (pflux_x * normal) + (pflux_y * t1) + (pflux_z * t2); //CoordinateVector<>(pflux_x, pflux_y, pflux_z);
      Eflux = (right_flux[4] + S_R_star * starstar_conserved[4] - (S_R_star - S_R) * star_conserved[4] - S_R * right_conserved[4]);
      EintFlux = (right_flux[9] + S_R_star * starstar_conserved[9] - (S_R_star - S_R) * star_conserved[9] - S_R * right_conserved[9]);
      double Bx_flux = left_conserved[8] + 0.5 * (right_conserved[8] - left_conserved[8]) - 0.5 * _cleaning_speed_max * (right_conserved[5] - left_conserved[5]) ; // HLLD Bx flux
      double Bscalar_flux = (_cleaning_speed_max ) * (left_conserved[5] + 0.5*(right_conserved[5]-left_conserved[5]) - 0.5 * (right_conserved[8] - left_conserved[8])); // HLLD scalar B flux
      double By_flux = right_flux[6] + S_R_star * starstar_conserved[6] - (S_R_star - S_R) * star_conserved[6] - S_R * right_conserved[6];
      double Bz_flux = right_flux[7] + S_R_star * starstar_conserved[7] - (S_R_star - S_R) * star_conserved[7] - S_R * right_conserved[7];
      Bflux = (Bx_flux * normal) + (By_flux * t1) + (Bz_flux * t2); //CoordinateVector<>(Bx_flux, By_flux, Bz_flux); // HLLD B flux
      B_scalarflux = Bscalar_flux;

      /*if (mach < _mach_limit){
        double interface_u = 0.5 * (uLn + uRn);
        double interface_rho = 0.5 * (rhoL + rhoR);
        double ekin_flux = 0.5 * interface_rho * interface_u * interface_u * interface_u;
        double B_squared = Bflux.x()*Bflux.x() + Bflux.y()*Bflux.y() + Bflux.z()*Bflux.z();
        double B_dot_v = Bflux.x()*interface_u; 
        double emag_flux = (B_squared / (2.0 * _mu0)) * interface_u - (Bflux.x() * B_dot_v) / _mu0;
        
        EintFlux = Eflux - ekin_flux - emag_flux;
      }*/

      
    
     // if (std::abs(Eflux) > 1e20) Eflux = std::copysign(1e20, Eflux); // debug mgb 03.05.2026
     if (std::isnan(Bflux.x()) || std::isnan(Bflux.y()) || std::isnan(Bflux.z()) || std::isnan(B_scalarflux)) {
        std::cout << "NaN detected in B flux calculation in HLLD S_M<=0, S_R_star>=0 solver!" << std::endl;
        std::cout << "B flux before HLLD: B_flux = " << Bfluxog.norm() << "Bx_flux " << Bx_flux << ", By_flux = " << By_flux << ", Bz_flux = " << Bz_flux << ", Bscalar_flux = " << Bscalar_flux << std::endl;
        std::cout << "P flux before HLLD: p_flux = " << pfluxog.norm() << std::endl;
        std::cout<< "B flux after HLLD: B_flux = " << Bflux.norm() << ", Bx_flux " << Bx_flux << ", By_flux = " << By_flux << ", Bz_flux = " << Bz_flux << ", Bscalar_flux = " << Bscalar_flux << std::endl;
        std::cout<< "P flux after HLLD: p_flux = " << pflux.norm() << std::endl;

      
    };
    if (std::abs(Eflux) > current_etot && debug_flag == true) {
        std::cout << "Warning: S_M <= 0.0 && S_R_star >= 0.0." << std::endl;
        std::cout << "HLLD: Mflux: " << mflux << ", EintFlux: " << EintFlux << ", Eflux: " << Eflux << ", Current Eint: " << current_eint << std::endl;
        std::cout << "HLLD: MASS left_flux: " << left_flux[0] << ", right_flux: " << right_flux[0] << ", right_conserved: " << right_conserved[0] << ", left_conserved: " << left_conserved[0] << std::endl;
        std::cout << "HLLD: ENERGY left_flux: " << left_flux[4] << ", right_flux: " << right_flux[4] << ", right_conserved: " << right_conserved[4] << ", left_conserved: " << left_conserved[4] << std::endl;
        std::cout << "HLLD: S_R: " << S_R << ", S_L: " << S_L << std::endl;
      };
    if (std::abs(EintFlux) > current_eint && debug_flag == true) {
        std::cout << "Warning: S_M <= 0.0 && S_R_star >= 0.0." << std::endl;
        std::cout << "HLLD: Mflux: " << mflux << ", EintFlux: " << EintFlux << ", Eflux: " << Eflux << ", Current Eint: " << current_eint << std::endl;
        std::cout << "HLLD: MASS left_flux: " << left_flux[0] << ", right_flux: " << right_flux[0] << ", right_conserved: " << right_conserved[0] << ", left_conserved: " << left_conserved[0] << std::endl;
        std::cout << "HLLD: ENERGY left_flux: " << left_flux[4] << ", right_flux: " << right_flux[4] << ", right_conserved: " << right_conserved[4] << ", left_conserved: " << left_conserved[4] << std::endl;
        std::cout << "HLLD: INTERNAL_ENERGY left_flux: " << left_flux[9] << ", right_flux: " << right_flux[9] << std::endl;
        std::cout << "HLLD: S_R: " << S_R << ", S_L: " << S_L << std::endl;
      }
    return;
    
    
  };
};
};


#endif // HLLDRIEMANNSOLVER_HPP
