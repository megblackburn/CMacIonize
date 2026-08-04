/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018, 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file Hydro.hpp
 *
 * @brief Hydro related functionality.
 *
 * @author Meg Blackburn mgb27@st-andrews.ac.uk
 */
#ifndef MAGNETOHYDRO_HPP
#define MAGNETOHYDRO_HPP

#include "HLLDRiemannSolver.hpp"
#include "MagnetoHydroBoundary.hpp" // mgb edit 21.04.2026 from HydroBoundary.hpp to MagnetoHydroBoundary.hpp
#include "MagnetoHydroVariables.hpp" // mgb edit 21.04.2026 from HydroVariables.hpp to MagnetoHydroVariables.hpp
#include "IonizationVariables.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"
#include "Abundances.hpp"

#include <cfloat>
#include <iomanip>

/*! @brief Uncomment this to enable hard resets for unphysical hydro
 *  variables. */
#define SAFE_HYDRO_VARIABLES

/*! @brief Uncomment this to activate the flux limiter. */
#define FLUX_LIMITER 2.

/**
 * @brief MagnetoHydro related functionality.
 */
class MagnetoHydro {
private:
  /*! @brief Polytropic index @f$\gamma{}@f$ of the gas. */
  const double _gamma;

  /*! @brief Assumed temperature for neutral gas (in K). */
  const double _neutral_temperature;

  /*! @brief Assumed temperature for ionised gas (in K). */
  const double _ionised_temperature;

  /*! @brief Velocity limit. Gas velocities higher than this value are capped
   *  (in m s^-1). */
  const double _max_velocity;

  /*! @brief Enable explicit radiation heating? */
  const bool _do_explicit_heating;

  /*! @brief Initial magnetic field strength (in T) (mgb edit 21.04.2026) */  
  const double _init_Bx; // mgb edit 21.04.2026
  const double _init_By; // mgb edit 21.04.2026
  const double _init_Bz; // mgb edit 21.04.2026
  const double _init_B_scalar; // mgb edit 21.04.2026 - initial scalar magnetic field strength (in T)

  /*! @brief Cleaning speed factor for divergence cleaning (default: 1.0) (mgb edit 21.04.2026) */
  const double _cleaning_speed_factor; // mgb edit 21.04.2026

  /*! @brief Mach limit (default: 1.0) (mgb edit 18.05.2026) */
  const double _mach_limit; // mgb edit 21.04.2026

  const Abundances &_abundances;

  /*! @brief @f$\gamma{}-1@f$. */
  const double _gamma_minus_one;

  /*! @brief @f$\frac{1}{\gamma{}-1}@f$. */
  const double _one_over_gamma_minus_one;

  /*! @brief mu0 */
  const double _mu0; // mgb edit 21.04.2026 - magnetic permeability of free space (in H m^-1)

  /*! @brief Conversion factor from number density to density (in kg). */
  const double _density_conversion_factor;

  /*! @brief Conversion factor from temperature to pressure,
   *  @f$P_{fac} = \frac{k}{m_{\rm{}H}}@f$ (in m^2 K^-1 s^-2). */
  const double _pressure_conversion_factor;

  /*! @brief Conversion factor from density to number density,
   *  @f$n_{fac} = \frac{1}{m_{\rm{}H}}@f$ (in kg^-1). */
  const double _n_conversion_factor;

  /*! @brief Conversion factor from pressure to temperature,
   *  @f$T_{fac} = \frac{m_{\rm{}H}}{k}@f$ (in K s^2 m^-2). */
  const double _T_conversion_factor;

  /*! @brief Conversion factor from temperature to internal energy,
   *  @f$u_{fac} = \frac{2k}{(\gamma{}-1)m_{\rm{}H}}@f$ (in m^2 K^-1 s^-2). */
  const double _u_conversion_factor;

  /*! @brief Conversion factor from magnetic field strength to energy density,
   *  @f$B_{fac} = \frac{1}{2\mu_0}@f$ (in kg m^-1 s^-2 T^-2). */
  const double _B_conversion_factor; // mgb edit 21.04.2026

  /*! @brief Riemann solver used to solve the Riemann problem. */
  const HLLDRiemannSolver _riemann_solver; // mgb edit 21.04.2026 - HLLD Riemann solver for MHD (changed from HLLCRiemannSolver)

  /**
   * @brief Per face slope limiter for a single quantity.
   *
   * Based on the slope limiter described in one of the appendices of Hopkins'
   * GIZMO paper.
   *
   * @param phimid0 Reconstructed value of the quantity at the interface.
   * @param phiL Value at the left of the interface.
   * @param phiR Value at the right of the interface.
   * @param dnrm_over_r Ratio of the distance between the left cell midpoint
   * to the face midpoint and the distances between left and right cell
   * midpoint.
   * @return Limited value of the quantity at the interface.
   */
  inline static double limit(const double phimid0, const double phiL,
                             const double phiR, const double dnrm_over_r) {

    const static double psi1 = 0.5;
    const static double psi2 = 0.25;

    const double delta1 = psi1 * std::abs(phiL - phiR);
    const double delta2 = psi2 * std::abs(phiL - phiR);

    const double phimin = std::min(phiL, phiR);
    const double phimax = std::max(phiL, phiR);

    const double phibar = phiL + dnrm_over_r * (phiR - phiL);

    // if sign(phimax+delta1) == sign(phimax)
    double phiplus;
    if ((phimax + delta1) * phimax > 0.) {
      phiplus = phimax + delta1;
    } else {
      const double absphimax = std::abs(phimax);
      phiplus = phimax * absphimax / (absphimax + delta1 + DBL_MIN);
    }

    // if sign(phimin-delta1) == sign(phimin)
    double phiminus;
    if ((phimin - delta1) * phimin > 0.) {
      phiminus = phimin - delta1;
    } else {
      const double absphimin = std::abs(phimin);
      phiminus = phimin * absphimin / (absphimin + delta1 + DBL_MIN);
    }

    double phimid;
    if (phiL == phiR) {
      phimid = phiL;
    } else {
      if (phiL < phiR) {
        phimid = std::max(phiminus, std::min(phibar + delta2, phimid0));
      } else {
        phimid = std::min(phiplus, std::max(phibar - delta2, phimid0));
      }
    }
    return phimid;
  }

public:
  /**
   * @brief Constructor.
   *
   * @param gamma Polytropic index @f$\gamma{}@f$ of the gas.
   * @param neutral_temperature Assumed neutral temperature for the gas (in K).
   * @param ionised_temperature Assumed ionised temperature for the gas (in K).
   * @param max_velocity Maximum allowed velocity for the gas (in m s^-1).
   * @param do_explicit_heating Enable explicit radiation heating?
   * @param init_Bx Initial x-component of the magnetic field (in T) (mgb edit 21.04.2026)
   * @param init_By Initial y-component of the magnetic field (in T) (mgb edit 21.04.2026)
   * @param init_Bz Initial z-component of the magnetic field (in T) (mgb edit 21.04.2026)
   * @param init_B_scalar Initial scalar magnetic field strength (in T) (mgb edit 21.04.2026)
   * @param cleaning_speed_factor Factor to multiply the cleaning speed by (default: 1.0) (mgb edit 21.04.2026)
   */
  inline MagnetoHydro(const double gamma, const double neutral_temperature,
               const double ionised_temperature, const double max_velocity,
               const bool do_explicit_heating, const double init_Bx, const double init_By, const double init_Bz, const double init_B_scalar, const double cleaning_speed_factor, const double mach_limit, const Abundances &abundances)
      : _gamma(gamma), _neutral_temperature(neutral_temperature),
        _ionised_temperature(ionised_temperature), _max_velocity(max_velocity),
        _do_explicit_heating(do_explicit_heating), _init_Bx(init_Bx), _init_By(init_By), _init_Bz(init_Bz), _init_B_scalar(init_B_scalar), _cleaning_speed_factor(cleaning_speed_factor), _mach_limit(mach_limit), _abundances(abundances),
        _gamma_minus_one(_gamma - 1.),
        _one_over_gamma_minus_one(1. / _gamma_minus_one),
        _mu0(PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_MAGNETC_PERMEABILITY)), // mgb edit 21.04.2026
        _density_conversion_factor(PhysicalConstants::get_physical_constant(
            PHYSICALCONSTANT_PROTON_MASS)),
        _pressure_conversion_factor(PhysicalConstants::get_physical_constant(
                                        PHYSICALCONSTANT_BOLTZMANN) /
                                    PhysicalConstants::get_physical_constant(
                                        PHYSICALCONSTANT_PROTON_MASS)),
        _n_conversion_factor(1. / PhysicalConstants::get_physical_constant(
                                      PHYSICALCONSTANT_PROTON_MASS)),
        _T_conversion_factor(PhysicalConstants::get_physical_constant(
                                 PHYSICALCONSTANT_PROTON_MASS) /
                             PhysicalConstants::get_physical_constant(
                                 PHYSICALCONSTANT_BOLTZMANN)),
        _u_conversion_factor(2. *
                             PhysicalConstants::get_physical_constant(
                                 PHYSICALCONSTANT_BOLTZMANN) *
                             _one_over_gamma_minus_one /
                             PhysicalConstants::get_physical_constant(
                                 PHYSICALCONSTANT_PROTON_MASS)),
        _B_conversion_factor(1./(2. * _mu0)), // mgb edit 21.04.2026
        _riemann_solver(gamma, _max_velocity, _cleaning_speed_factor) {}

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
  inline MagnetoHydro(const Abundances &abundances, ParameterFile &params)
      : MagnetoHydro(params.get_value< double >("MagnetoHydro:polytropic index", 5. / 3.),
              params.get_physical_value< QUANTITY_TEMPERATURE >(
                  "MagnetoHydro:neutral temperature", "100. K"),
              params.get_physical_value< QUANTITY_TEMPERATURE >(
                  "MagnetoHydro:ionised temperature", "1.e4 K"),
              params.get_physical_value< QUANTITY_VELOCITY >(
                  "MagnetoHydro:maximum velocity", "1.e99 m s^-1"),
              params.get_value< bool >("MagnetoHydro:do explicit heating", false),
              params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                  "MagnetoHydro:init Bx", "0.0 T"), // mgb edit 21.04.2026 - initial x-component of the magnetic field (in T)
              params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                  "MagnetoHydro:init By", "0.0 T"), // mgb edit 21.04.2026 - initial y-component of the magnetic field (in T)
              params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                  "MagnetoHydro:init Bz", "0.0 T"), // mgb edit 21.04.2026 - initial z-component of the magnetic field (in T)
              params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                  "MagnetoHydro:init B_scalar", "0.0 T"), // mgb edit 21.04.2026 - initial scalar magnetic field strength (in T)
              params.get_value<double>("MagnetoHydro:cleaning speed factor",1.0), // mgb edit 21.04.2026 - factor to multiply the cleaning speed by (default: 1.0)
              params.get_value<double>("MagnetoHydro: mach limit", 3.0),
              abundances) {}

  /**
   * @brief Get the soundspeed for the given hydrodynamic variables.
   *
   * @param magnetohydro_variables Hydro variables.
   * @param ionization_variables IonizationVariables for the same cell.
   * @return Soundspeed (in m s^-1).
   */
  inline double
  get_soundspeed(const MagnetoHydroVariables &magnetohydro_variables,
                 const IonizationVariables &ionization_variables) const {

    if (_gamma > 1.) {
      const double rho = magnetohydro_variables.get_primitives_density();
      const double P = magnetohydro_variables.get_primitives_pressure();
      if (rho > 0. && P > 0.) {
        const double rho_inv = 1. / rho;
        if (!std::isinf(rho_inv)) {
          const double cs = std::sqrt(_gamma * P * rho_inv);
          cmac_assert(cs == cs);
          cmac_assert_message(cs > 0., "gamma: %g, rho: %g, rho_inv: %g, P: %g",
                              _gamma, rho, rho_inv, P);
          return cs;
        } else {
          return DBL_MIN;
        }
      } else {
        return DBL_MIN;
      }
    } else {
      const double mean_molecular_mass =
          0.5 * (1. + ionization_variables.get_ionic_fraction(ION_H_n));
      const double temperature = ionization_variables.get_temperature();
      const double cs = std::sqrt(_pressure_conversion_factor * temperature /
                                  mean_molecular_mass);
      cmac_assert(cs == cs);
      cmac_assert_message(cs > 0., "xH: %g, T: %g",
                          ionization_variables.get_ionic_fraction(ION_H_n),
                          temperature);
      return cs;
    }
  }

  /**
   * @brief Get the maximum MHD wave speed for the given hydrodynamic variables.
   *
   * @param magnetohydro_variables Hydro variables.
   * @param ionization_variables IonizationVariables for the same cell.
   * @return Soundspeed (in m s^-1).
   */
  inline double
  get_vMHD(const MagnetoHydroVariables &magnetohydro_variables,
                 const IonizationVariables &ionization_variables) const {

      const double cs = get_soundspeed(magnetohydro_variables, ionization_variables);
      double rho = magnetohydro_variables.get_primitives_density();
      if (rho <= 0.) { 
        return cs;
      }
      double Bx = magnetohydro_variables.get_primitives_magnetic_field().x();
      double By = magnetohydro_variables.get_primitives_magnetic_field().y();
      double Bz = magnetohydro_variables.get_primitives_magnetic_field().z();
      double B_squared = Bx*Bx + By*By + Bz*Bz;

      double vA_squared = B_squared/(_mu0 * rho); // mgb edit 21.04.2026 - Alfvén speed squared + sound speed squared
      double vMHD = std::sqrt(cs * cs + vA_squared);
      cmac_assert(vMHD == vMHD);
      return vMHD;
  }
  

  /**
   * @brief Set the primitive variables for the given state and inverse volume.
   *
   * @param magnetohydro_state Hydrodynamical state variables.
   * @param ionization_state Ionization variables.
   * @param inverse_volume Inverse of the volume (in m^-3).
   */
  inline void set_primitive_variables(MagnetoHydroVariables &magnetohydro_state,
                                      IonizationVariables &ionization_state,
                                      const double inverse_volume) const { 

    double density = 0.;
    CoordinateVector<> velocity;
    double pressure = 0.;
    CoordinateVector<>B_field = magnetohydro_state.get_conserved_magnetic_field() * inverse_volume;
    double Bx = B_field.x(); // mgb edit 21.04.2026 - initial x-component of the magnetic field (in T)
    double By = B_field.y(); // mgb edit 21.04.2026 - initial y-component of the magnetic field (in T)
    double Bz = B_field.z(); // mgb edit 21.04.2026 - initial z-component of the magnetic field (in T)
    double B_scalar = magnetohydro_state.get_conserved_B_scalar() * inverse_volume; // mgb edit 21.04.2026 - initial scalar magnetic field strength (in T)
    
    double initial_P = magnetohydro_state.get_primitives_pressure();
    double initial_rho = magnetohydro_state.get_primitives_density();

    if (std::isnan(Bx) || std::isnan(By) || std::isnan(Bz)) {
        std::cout << "Warning: Magnetic field is NaN before primitive variable conversion for pressure:" << initial_P << ", density:" << initial_rho << std::endl;
    }
    if (std::isnan(initial_P) || std::isnan(initial_rho)) {
        std::cout << "Warning: Pressure or density is NaN before primitive variable conversion for pressure:" << initial_P << ", density:" << initial_rho << std::endl;
    }
    double internal_energy;

    if (magnetohydro_state.get_conserved_mass() > 0.) {
      const double inverse_mass = 1. / magnetohydro_state.get_conserved_mass();
      if (!std::isinf(inverse_mass)) {
        density = magnetohydro_state.get_conserved_mass() * inverse_volume; 
        velocity = inverse_mass * magnetohydro_state.get_conserved_momentum();

        double e_int = magnetohydro_state.get_conserved_internal_energy(); // mgb edit 06.05.2026 - get internal energy density
        double test_pressure = _gamma_minus_one * e_int * inverse_volume; // mgb edit 06.05.2026 - calculate pressure from internal energy density
        double cs_squared = _gamma * test_pressure / density; // mgb edit 06.05.2026 - calculate sound speed squared from test pressure and density
        double va_squared = (Bx*Bx + By*By + Bz*Bz) / (_mu0 * density); // mgb edit 06.05.2026 - calculate Alfvén speed squared from magnetic field strength and density
        double v_fast = std::sqrt(cs_squared + va_squared);
        double v_mag = velocity.norm();

        double mach_fast = v_mag/(v_fast + 1.e-15); // mgb edit 06.05.2026 - calculate fast magnetosonic Mach number (add small number to avoid division by zero

        double ekin = 0.5 * CoordinateVector<>::dot_product(
                         velocity, magnetohydro_state.get_conserved_momentum());
        double emag = (1./inverse_volume) *0.5 * (Bx*Bx + By*By + Bz*Bz)/_mu0; // mgb edit 21.04.2026 - magnetic energy density
        double etot = emag + ekin + e_int;
        if (_gamma > 1.) {
          

          if (std::isnan(ekin) || std::isnan(emag)) {
            std::cout << "Warning: Kinetic or magnetic energy is NaN before primitive variable conversion, emag:" << emag << ", ekin:" << ekin << std::endl;
          }

          if (mach_fast > _mach_limit || (emag + ekin > 0.95 * magnetohydro_state.get_conserved_total_energy()) || (test_pressure-magnetohydro_state.get_primitives_pressure())/(std::min(test_pressure, magnetohydro_state.get_primitives_pressure())+1e-10) > 0.25) { // mgb edit 06.05.2026 - if the fast magnetosonic Mach number is very low or magnetic energy contribution is high, ignore magnetic energy contribution to internal energy to avoid negative pressures due to numerical errors in total energy when kinetic and magnetic energy contributions are close to total energy
         //   std::cout << "Warning: dual energy triggered, emag:" << emag << ", ekin:" << ekin << ", total energy:" << magnetohydro_state.get_conserved_total_energy() << std::endl;
            internal_energy = magnetohydro_state.get_conserved_internal_energy(); // mgb edit 21.04.2026 - subtract kinetic and magnetic energy contributions from total energy to get internal energy
            magnetohydro_state.set_primitives_internal_energy(internal_energy * inverse_volume); // mgb edit 21.04.2026 - update conserved internal energy with new internal energy value (after subtracting kinetic and magnetic energy contributions)  
            if (internal_energy > etot) {
              std::cout << "Warning: internal energy is GREATER than total energy in set_primitive_variables dual energy! Eint: " << internal_energy << ", Etot: " << magnetohydro_state.get_conserved_total_energy() << ", Emag: " << emag << ", Ekin: " << ekin << std::endl;
            }
            if (internal_energy > etot) {
              std::cout << "Warning: tracked internal energy is GREATER than internal energy in set_primitive_variables dual energy! Eint: " << internal_energy << ", Etot: " << magnetohydro_state.get_conserved_total_energy() << ", Emag: " << emag << ", Ekin: " << ekin << std::endl;
            }
            if (internal_energy < 0) {
              std::cout << "Warning: set_primitive_variables dual energy, Eint < 0: " << internal_energy << std::endl;
            }

            //  std::cout << "Warning: dual energy, Bx = " << magnetohydro_state.get_primitives_magnetic_field().x() << ", By = " << magnetohydro_state.get_primitives_magnetic_field().y() << ", Bz = " << magnetohydro_state.get_primitives_magnetic_field().z() << ", 1/inverse volume = " << 1./inverse_volume  << std::endl;
          //  std::cout << "Warning: dual energy triggered, emag:" << emag << ", ekin:" << ekin << ", eint:" << internal_energy << ", total energy:" << magnetohydro_state.get_conserved_total_energy() << std::endl;
          } else {
            internal_energy = magnetohydro_state.get_conserved_total_energy() - ekin - emag; // mgb edit 21.04.2026 - subtract kinetic energy contribution from total energy to get internal energy (ignore magnetic energy contribution for very low internal energies to avoid negative pressures)
            magnetohydro_state.set_conserved_internal_energy(internal_energy); // mgb edit 21.04.2026 - update conserved internal energy with new internal energy value (after subtracting kinetic and magnetic energy contributions)
          //  magnetohydro_state.set_primitives_internal_energy(internal_energy*inverse_volume); // mgb edit 21.04.2026 - update primitive internal energy with new internal energy value (after subtracting kinetic and magnetic energy contributions)
          //  std::cout << "Warning: dual energy not triggered, Bx = " << magnetohydro_state.get_primitives_magnetic_field().x() << ", By = " << magnetohydro_state.get_primitives_magnetic_field().y() << ", Bz = " << magnetohydro_state.get_primitives_magnetic_field().z() << ", 1/inverse volume = " << 1./inverse_volume  << std::endl;
          //  std::cout << "Warning: dual energy not triggered, emag:" << emag << ", ekin:" << ekin << ", eint:" << internal_energy << ", total energy:" << magnetohydro_state.get_conserved_total_energy() << std::endl;
            if (magnetohydro_state.get_conserved_internal_energy() > magnetohydro_state.get_conserved_total_energy()) {
              std::cout << "Warning: internal energy is GREATER than total energy in set_primitive_variables! Eint: " << internal_energy << ", Etot: " << magnetohydro_state.get_conserved_total_energy() << ", Emag: " << emag << ", Ekin: " << ekin << std::endl;
            }
            if (internal_energy < 0) {
              std::cout << "Warning: set_primitive_variables, Eint < 0: " << internal_energy << std::endl;
            }
          }

          pressure =
              _gamma_minus_one * inverse_volume *
              internal_energy; // mgb edit 21.04.2026 - subtract magnetic energy density from total energy to get thermal energy density
        } else {
          std::cout<< "Warning: Temperature is being used in set_primitive_variables - this temp will now be out of sync with the internal energy!" << std::endl;
          const double mean_molecular_mass =
              0.5 * (1. + ionization_state.get_ionic_fraction(ION_H_n));
          const double temperature = ionization_state.get_temperature();
          pressure = _pressure_conversion_factor * density * temperature /
                     mean_molecular_mass;
        }

        // apply velocity limiter
        const double vnrm = velocity.norm();
        if (vnrm > _max_velocity) {
          velocity *= (_max_velocity / vnrm);
        }
        if (density > 0.) {
          const double inverse_density = 1. / density;
          if (!std::isinf(inverse_density)) {
            const double cs = std::sqrt(_gamma * pressure * inverse_density);
            const double B_squared = Bx*Bx + By*By + Bz*Bz;
            const double vA_squared = B_squared/(_mu0 * density); // mgb edit 21.04.2026 - Alfvén speed squared
            const double vMHD = std::sqrt(cs * cs + vA_squared); // mgb edit 21.04.2026 - MHD wave speed
            if (vMHD > _max_velocity) {
              const double factor = _max_velocity / vMHD;
              pressure *= factor * factor;

              if (std::isnan(pressure)) {
                std::cout << "Warning: Pressure is NaN, P: " << pressure << std::endl;
              }
            }
          }
        }

        cmac_assert(density == density);
        cmac_assert(velocity.x() == velocity.x());
        cmac_assert(velocity.y() == velocity.y());
        cmac_assert(velocity.z() == velocity.z());
        cmac_assert(_gamma == 1. || pressure == pressure);


#ifdef SAFE_HYDRO_VARIABLES
        density = std::max(density, 0.);
        pressure = std::max(pressure, 0.);
        internal_energy = std::max(internal_energy, 0.);
#else
        cmac_assert(density >= 0.);
        cmac_assert(_gamma == 1. || pressure >= 0.);
#endif
      }
    }
    magnetohydro_state.set_primitives_density(density);
    magnetohydro_state.set_primitives_velocity(velocity);
    cmac_assert_message(pressure > 0.0,"ABOUT TO SET P=0, energy=%g,kin=%g,mag=%g",magnetohydro_state.get_conserved_total_energy(),0.5 * CoordinateVector<>::dot_product(
                         velocity, magnetohydro_state.get_conserved_momentum()), emag);
    magnetohydro_state.set_primitives_pressure(pressure);
    magnetohydro_state.set_primitives_magnetic_field(CoordinateVector<>(Bx, By, Bz)); // mgb edit 21.04.2026 - set magnetic field components
    magnetohydro_state.set_primitives_B_scalar(B_scalar); // mgb edit 21.04.2026 - set B scalar field
    magnetohydro_state.set_primitives_internal_energy(magnetohydro_state.get_conserved_internal_energy()* inverse_volume); // mgb edit 21.04.2026 - set primitive internal energy from conserved internal energy (after subtracting kinetic and magnetic energy contributions if necessary)
  
    if (std::isnan(magnetohydro_state.get_conserved_total_energy())) {
      std::cerr << "Total energy is NaN before conserved update in set_primitive_variables for density " << density << ", pressure " << pressure << ", B field " << Bx << ", " << By << ", " << Bz << std::endl;
    }
    if (magnetohydro_state.get_conserved_total_energy() < magnetohydro_state.get_conserved_internal_energy()){
        std::cout<<"Warning! Etot < Eint before update in set primitive variables! Etot: " << magnetohydro_state.get_conserved_total_energy() << ", Eint: " << magnetohydro_state.get_conserved_internal_energy() << std::endl;
      }

    set_conserved_variables(magnetohydro_state, 1./inverse_volume);

    if (magnetohydro_state.get_conserved_total_energy() < magnetohydro_state.get_conserved_internal_energy()){
        std::cout<<"Warning! Etot < Eint after update in set primitive variables! Etot: " << magnetohydro_state.get_conserved_total_energy() << ", Eint: " << magnetohydro_state.get_conserved_internal_energy() << std::endl;
      }

    
    if (std::isnan(magnetohydro_state.get_conserved_total_energy())) {
      std::cerr << "Total energy is NaN after conserved update in set_primitive_variables for density " << density << ", pressure " << pressure << ", B field " << Bx << ", " << By << ", " << Bz << std::endl;
    }

    if (magnetohydro_state.get_conserved_internal_energy() > magnetohydro_state.get_conserved_total_energy()) {
      std::cout << "Warning: internal energy is GREATER than total energy after set_primitive_variables! Eint: " << magnetohydro_state.get_conserved_internal_energy() << ", Etot: " << magnetohydro_state.get_conserved_total_energy() << std::endl;
    }

  //  std::cout << "Primitive and conserved B field set Line 407, B field = " << Bx << ", " << By << ", " << Bz << std::endl;
  }

  /**
   * @brief Set the conserved variables for the given state and volume.
   *
   * @param state State variables.
   * @param volume Volume (in m^-3).
   */
  inline void set_conserved_variables(MagnetoHydroVariables &state,
                                      const double volume) const {

    double mass = state.get_primitives_density() * volume;
    const CoordinateVector<> momentum = mass * state.get_primitives_velocity();
    double total_energy = 0.;
    double internal_energy = 0.;
    const CoordinateVector<> magnetic_field = state.get_primitives_magnetic_field(); // mgb edit 21.04.2026 - get magnetic field strength from primitive magnetic field
    if (_gamma > 1.) {
      internal_energy = _one_over_gamma_minus_one * state.get_primitives_pressure() * volume; // mgb edit 06.05.2026 - calculate internal energy from primitive pressure
      double emag = 0.5 * CoordinateVector<>::dot_product(state.get_primitives_magnetic_field(), state.get_primitives_magnetic_field())/ _mu0;
      double ekin = 0.5 * CoordinateVector<>::dot_product(
                    momentum, state.get_primitives_velocity());
      total_energy =
          internal_energy +
          ekin + emag*volume; // mgb edit 21.04.2026 - add magnetic energy to total energy
     // if (internal_energy > state.get_conserved_total_energy()) {
       // std::cout << "Warning: internal energy is GREATER than total energy before set_conserved_variables! Eint: " << internal_energy << ", Etot: " << state.get_conserved_total_energy() << std::endl;
     // }
    }
    

    cmac_assert(mass == mass);
    cmac_assert(momentum.x() == momentum.x());
    cmac_assert(momentum.y() == momentum.y());
    cmac_assert(momentum.z() == momentum.z());
    cmac_assert(_gamma == 1. || total_energy == total_energy);

    if (internal_energy < 0.0) {
      std::cout << "Warning: set_conserved_variables, Eint < 0: " << internal_energy << std::endl;
    }

#ifdef SAFE_HYDRO_VARIABLES
    mass = std::max(mass, 0.);
    total_energy = std::max(total_energy, 0.);
    internal_energy = std::max(internal_energy, 0.);
#else
    cmac_assert(mass >= 0.);
    cmac_assert(_gamma == 1. || total_energy >= 0.);
#endif

    state.set_conserved_mass(mass);
    state.set_conserved_momentum(momentum);
    state.set_conserved_total_energy(total_energy);
    state.set_conserved_magnetic_field(magnetic_field*volume); // mgb edit 21.04.2026 - set conserved magnetic field from primitive magnetic field
    state.set_conserved_B_scalar(state.get_primitives_B_scalar() * volume); // mgb edit 21.04.2026 - set conserved B scalar field from primitive B scalar field
    state.set_conserved_internal_energy(internal_energy)  ; // mgb edit 21.04.2026 - set conserved internal energy from calculated internal energy

    if (state.get_conserved_internal_energy() > state.get_conserved_total_energy()) {
      std::cout << "Warning: internal energy is GREATER than total energy after set_conserved_variables! Eint: " << internal_energy << ", Etot: " << state.get_conserved_total_energy()  << std::endl;
    }
  }

  /**
   * @brief Do the flux calculation for the given interface.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param left_state Left state hydro variables.
   * @param right_state Right state hydro variables.
   * @param dx Distance between left and right state midpoint (in m).
   * @param A Surface area of the interface (in m^2).
   * @param dt Current system time step, used for flux limiter (in s).
   */
  inline void do_flux_calculation(const uint_fast8_t i,
                                  MagnetoHydroVariables &left_state,
                                  MagnetoHydroVariables &right_state, const double dx,
                                  const double A, const double dt) const {

    const double halfdx = 0.5 * dx;
    double rhoL = left_state.get_primitives_density() +
                  halfdx * left_state.primitive_gradients(0)[i];
    CoordinateVector<> vL(left_state.primitives(1) +
                              halfdx * left_state.primitive_gradients(1)[i],
                          left_state.primitives(2) +
                              halfdx * left_state.primitive_gradients(2)[i],
                          left_state.primitives(3) +
                              halfdx * left_state.primitive_gradients(3)[i]);
    double PL = left_state.get_primitives_pressure() + // mgb edit 21.04.2026 start
                halfdx * left_state.primitive_gradients(4)[i];
    CoordinateVector<> BL(left_state.primitives(5) +
                              halfdx * left_state.primitive_gradients(5)[i],
                          left_state.primitives(6) +
                              halfdx * left_state.primitive_gradients(6)[i],
                          left_state.primitives(7) +
                              halfdx * left_state.primitive_gradients(7)[i]); // mgb edit 21.04.2026 - reconstructed magnetic field components
    double B_scalarL = left_state.primitives(8) + // mgb edit 21.04.2026 - reconstructed B scalar field
                       halfdx * left_state.primitive_gradients(8)[i];
    double rhoR = right_state.get_primitives_density() -
                  halfdx * right_state.primitive_gradients(0)[i];
    CoordinateVector<> vR(right_state.primitives(1) -
                              halfdx * right_state.primitive_gradients(1)[i],
                          right_state.primitives(2) -
                              halfdx * right_state.primitive_gradients(2)[i],
                          right_state.primitives(3) -
                              halfdx * right_state.primitive_gradients(3)[i]);
    double PR = right_state.get_primitives_pressure() -
                halfdx * right_state.primitive_gradients(4)[i];
    CoordinateVector<> BR(right_state.primitives(5) - // mgb edit 21.04.2026 start
                              halfdx * right_state.primitive_gradients(5)[i],
                          right_state.primitives(6) -
                              halfdx * right_state.primitive_gradients(6)[i],
                          right_state.primitives(7) -
                              halfdx * right_state.primitive_gradients(7)[i]); // mgb edit 21.04.2026 - reconstructed magnetic field components
    double B_scalarR = right_state.primitives(8) - // mgb edit 21.04.2026 - reconstructed B scalar field
                       halfdx * right_state.primitive_gradients(8)[i];
  //  std::cout << "initial B scalar set: " << left_state.primitives(8) << ", " << right_state.primitives(8) << std::endl;
  //  std::cout << "initial B scalar values: B_scalarL: " << B_scalarL << ", B_scalarR: " << B_scalarR << std::endl;
    rhoL = limit(rhoL, left_state.get_primitives_density(),
                 right_state.get_primitives_density(), 0.5);
    vL[0] = limit(vL.x(), left_state.get_primitives_velocity().x(),
                  right_state.get_primitives_velocity().x(), 0.5);
    vL[1] = limit(vL.y(), left_state.get_primitives_velocity().y(),
                  right_state.get_primitives_velocity().y(), 0.5);
    vL[2] = limit(vL.z(), left_state.get_primitives_velocity().z(),
                  right_state.get_primitives_velocity().z(), 0.5);
    PL = limit(PL, left_state.get_primitives_pressure(),
               right_state.get_primitives_pressure(), 0.5);
    BL[0] = limit(BL.x(), left_state.get_primitives_magnetic_field().x(),
                  right_state.get_primitives_magnetic_field().x(), 0.5); // mgb edit 21.04.2026 - limit reconstructed magnetic field components
    BL[1] = limit(BL.y(), left_state.get_primitives_magnetic_field().y(),
                  right_state.get_primitives_magnetic_field().y(), 0.5);
    BL[2] = limit(BL.z(), left_state.get_primitives_magnetic_field().z(),
                  right_state.get_primitives_magnetic_field().z(), 0.5);
    B_scalarL = limit(B_scalarL, left_state.get_primitives_B_scalar(),
                      right_state.get_primitives_B_scalar(), 0.5); // mgb edit 21

    rhoR = limit(rhoR, right_state.get_primitives_density(),
                 left_state.get_primitives_density(), 0.5);
    vR[0] = limit(vR.x(), right_state.get_primitives_velocity().x(),
                  left_state.get_primitives_velocity().x(), 0.5);
    vR[1] = limit(vR.y(), right_state.get_primitives_velocity().y(),
                  left_state.get_primitives_velocity().y(), 0.5);
    vR[2] = limit(vR.z(), right_state.get_primitives_velocity().z(),
                  left_state.get_primitives_velocity().z(), 0.5);
    PR = limit(PR, right_state.get_primitives_pressure(),
               left_state.get_primitives_pressure(), 0.5);
    BR[0] = limit(BR.x(), right_state.get_primitives_magnetic_field().x(),
                  left_state.get_primitives_magnetic_field().x(), 0.5); // mgb edit 21.04.2026 - limit reconstructed magnetic field components
    BR[1] = limit(BR.y(), right_state.get_primitives_magnetic_field().y(),
                  left_state.get_primitives_magnetic_field().y(), 0.5);
    BR[2] = limit(BR.z(), right_state.get_primitives_magnetic_field().z(),
                  left_state.get_primitives_magnetic_field().z(), 0.5);
    B_scalarR = limit(B_scalarR, right_state.get_primitives_B_scalar(),
                      left_state.get_primitives_B_scalar(), 0.5); // mgb edit 21.04.2026 - limit reconstructed B scalar field

  //  std::cout << "limited B scalar values: B_scalarL: " << B_scalarL << ", B_scalarR: " << B_scalarR << std::endl;
    cmac_assert(rhoL == rhoL);
    cmac_assert(vL.x() == vL.x());
    cmac_assert(vL.y() == vL.y());
    cmac_assert(vL.z() == vL.z());
    cmac_assert(PL == PL);
    cmac_assert(BL.x() == BL.x());
    cmac_assert(BL.y() == BL.y());
    cmac_assert(BL.z() == BL.z());
    cmac_assert(B_scalarL == B_scalarL);

    cmac_assert(rhoR == rhoR);
    cmac_assert(vR.x() == vR.x());
    cmac_assert(vR.y() == vR.y());
    cmac_assert(vR.z() == vR.z());
    cmac_assert(PR == PR);
    cmac_assert(BR.x() == BR.x());
    cmac_assert(BR.y() == BR.y());
    cmac_assert(BR.z() == BR.z());
    cmac_assert(B_scalarR == B_scalarR);

    // make sure all densities and pressures are physical
#ifdef SAFE_HYDRO_VARIABLES
    rhoL = std::max(rhoL, 0.);
    PL = std::max(PL, 0.);
    rhoR = std::max(rhoR, 0.);
    PR = std::max(PR, 0.);
#else
    cmac_assert(rhoL >= 0.);
    cmac_assert(PL >= 0.);
    cmac_assert(rhoR >= 0.);
    cmac_assert(PR >= 0.);
#endif

    double mflux = 0.;
    CoordinateVector<> pflux;
    double Eflux = 0.;
    double EintFlux = 0.; // mgb edit 21.04.2026 - internal energy flux
    CoordinateVector<> Bflux; // mgb edit 21.04.2026 - magnetic field flux
    double B_scalarflux = 0.; // mgb edit 21.04.2026
    CoordinateVector<> normal;
    normal[i] = 1.;
    double current_etot = left_state.get_conserved_total_energy();
    double current_eint = left_state.get_conserved_internal_energy();
   // std::cout << "Pass to Riemann - " << "rhoL: " << rhoL << ", vL: " << vL.x() << ", " << vL.y() << ", " << vL.z() << ", PL: " << PL << ", BL: " << BL.x() << ", " << BL.y() << ", " << BL.z() << ", B_scalarL: " << B_scalarL
           //   << "; rhoR: " << rhoR << ", vR: " << vR.x() << ", " << vR.y() << ", " << vR.z() << ", PR: " << PR << ", BR: " << BR.x() << ", " << BR.y() << ", " << BR.z() << ", B_scalarR: " << B_scalarR
            //  << std::endl;
    if (current_eint == 0 || right_state.get_conserved_internal_energy()==0){
      std::cout << "Warning: Internal Energy is zero before entering riemann flux solver!" << std::endl;
    }
    _riemann_solver.solve_for_flux_MHD(rhoL, vL, PL, BL, B_scalarL, rhoR, vR, PR, BR, B_scalarR, mflux, pflux, Eflux, EintFlux, Bflux, B_scalarflux,
                                    normal, current_etot, _mach_limit, current_eint);
  //  std::cout << "Fluxes calculated Line 583 MagnetoHydro Bflux = " << Bflux.x() << ", " << Bflux.y() << ", " << Bflux.z() << std::endl;

    cmac_assert(mflux == mflux);
    cmac_assert(pflux.x() == pflux.x());
    cmac_assert(pflux.y() == pflux.y());
    cmac_assert(pflux.z() == pflux.z());
    cmac_assert(Eflux == Eflux);

    mflux *= A;
    pflux[0] *= A;
    pflux[1] *= A;
    pflux[2] *= A;
    Eflux *= A;
    EintFlux *= A; // mgb edit 21.04.2026 - scale internal energy flux by area
    Bflux[0] *= A; // mgb edit 21.04.2026 - scale magnetic field flux by area
    Bflux[1] *= A;
    Bflux[2] *= A;
    B_scalarflux *= A; // mgb edit 21.04.2026 - scale B scalar flux by area

#ifdef FLUX_LIMITER
    // limit the flux
    double fluxfac = 1.;
    const double absmflux = mflux * dt;
    if (absmflux > FLUX_LIMITER * left_state.get_conserved_mass()) {
      double safe_mass = std::max(left_state.get_conserved_mass(), 1e-30);
      fluxfac = FLUX_LIMITER * safe_mass / absmflux;
    }
    if (-absmflux > FLUX_LIMITER * right_state.get_conserved_mass()) {
      double safe_massr = std::max(right_state.get_conserved_mass(), 1e-30);
      fluxfac = std::min(
          fluxfac, -FLUX_LIMITER * safe_massr / absmflux);
    }
    if (_gamma > 1.) {
      const double absEflux = Eflux * dt;
      if (absEflux > FLUX_LIMITER * left_state.get_conserved_total_energy()) {
        fluxfac = std::min(
            fluxfac,
            FLUX_LIMITER * left_state.get_conserved_total_energy() / absEflux);
            cmac_assert_message(fluxfac >= 0. && fluxfac <= 1., "fluxfac1a: %g , absEflux = %g, eflux = %g, e_r = %g",
                                        fluxfac, absEflux, Eflux, left_state.get_conserved_total_energy());

      }
      if (-absEflux > FLUX_LIMITER * right_state.get_conserved_total_energy()) {
        fluxfac = std::min(
            fluxfac, -FLUX_LIMITER * right_state.get_conserved_total_energy() /
                         absEflux);
        cmac_assert_message(fluxfac >= 0. && fluxfac <= 1., "fluxfac1b: %g , absEflux = %g, eflux = %g, e_r = %g",
                                    fluxfac, absEflux, Eflux, right_state.get_conserved_total_energy());
      }
    }
    if (_gamma > 1.) {
      const double absEintFlux = EintFlux * dt;
      if (absEintFlux > FLUX_LIMITER * left_state.get_conserved_internal_energy()) {
        fluxfac = std::min(
            fluxfac,
            FLUX_LIMITER * left_state.get_conserved_internal_energy() / absEintFlux);
            cmac_assert_message(fluxfac >= 0. && fluxfac <= 1., "fluxfac1a: %g , absEflux = %g, eflux = %g, e_r = %g",
                                        fluxfac, absEintFlux, EintFlux, left_state.get_conserved_internal_energy());
      }
      
      // Check if we are "over-stealing" from the right cell (negative flux)
      if (-absEintFlux > FLUX_LIMITER * right_state.get_conserved_internal_energy()) {
        fluxfac = std::min(
            fluxfac, 
            -FLUX_LIMITER * right_state.get_conserved_internal_energy() / absEintFlux);
        cmac_assert_message(fluxfac >= 0. && fluxfac <= 1., "fluxfac1b: %g , absEflux = %g, eflux = %g, e_r = %g",
                                    fluxfac, absEintFlux, EintFlux, right_state.get_conserved_internal_energy());
      }
    }

    

    cmac_assert_message(fluxfac >= 0. && fluxfac <= 1., "fluxfac2: %g", fluxfac);
    // momentum flux limiter
    // note that we only apply this for cells that have high momentum, i.e.
    // whose momentum is higher than the thermal momentum of the cell
    // without this condition, cells with zero momentum would never be able
    // to gain momentum...
    const double p2 = left_state.get_conserved_momentum().norm2();
    const double m2 =
        left_state.get_conserved_mass() * left_state.get_conserved_mass();
    if (p2 * left_state.get_primitives_density() >
        _gamma * m2 * left_state.get_primitives_pressure()) {
      const double pflux2 = pflux.norm2() * dt * dt;
      if (pflux2 > (FLUX_LIMITER * FLUX_LIMITER) * p2) {
        fluxfac = std::min(
            fluxfac, std::sqrt((FLUX_LIMITER * FLUX_LIMITER) * p2 / pflux2));
      }
    }
    cmac_assert_message(fluxfac >= 0. && fluxfac <= 1., "fluxfac3: %g", fluxfac);
    {
      const double pn2 = right_state.get_conserved_momentum().norm2();
      const double mn2 =
          right_state.get_conserved_mass() * right_state.get_conserved_mass();
      if (p2 * right_state.get_primitives_density() >
          _gamma * mn2 * right_state.get_primitives_pressure()) {
        const double pflux2 = pflux.norm2() * dt * dt;
        if (pflux2 > (FLUX_LIMITER * FLUX_LIMITER) * pn2) {
          fluxfac = std::min(
              fluxfac, std::sqrt((FLUX_LIMITER * FLUX_LIMITER) * pn2 / pflux2));
        }
      }
    }
    cmac_assert_message(fluxfac >= 0. && fluxfac <= 1., "fluxfac4: %g", fluxfac);
    if (std::isnan(pflux.x()) || std::isnan(pflux.y()) || std::isnan(pflux.z())) {
      std::cout << "Warning: Momentum flux is NaN before flux limiter, pflux: " << pflux.x() << ", " << pflux.y() << ", " << pflux.z() << std::endl;
    }
    if (std::isnan(Bflux.x()) || std::isnan(Bflux.y()) || std::isnan(Bflux.z())) {
      std::cout << "Warning: Magnetic field flux is NaN before flux limiter, Bflux: " << Bflux.x() << ", " << Bflux.y() << ", " << Bflux.z() << std::endl;
    }

    
    mflux *= fluxfac;
    pflux *= fluxfac;
    Eflux *= fluxfac;
    EintFlux *= fluxfac; // mgb edit 21.04.2026 - apply flux limiter to internal energy flux
    Bflux[0] *= fluxfac; // mgb edit 21.04.2026 - apply flux limiter to magnetic field flux
    Bflux[1] *= fluxfac;
    Bflux[2] *= fluxfac;
    B_scalarflux *= fluxfac; // mgb edit 21.04.2026 - apply flux limiter to B scalar flux
#endif

    if (std::isnan(pflux.x()) || std::isnan(pflux.y()) || std::isnan(pflux.z())) {
      std::cout << "Warning: Momentum flux is NaN after flux limiter, pflux: " << pflux.x() << ", " << pflux.y() << ", " << pflux.z() << std::endl;
    }
    if (std::isnan(Bflux.x()) || std::isnan(Bflux.y()) || std::isnan(Bflux.z())) {
      std::cout << "Warning: Magnetic field flux is NaN after flux limiter, Bflux: " << Bflux.x() << ", " << Bflux.y() << ", " << Bflux.z() << std::endl;
    }

    const double safe_fraction = 0.8;
    double axis_scale_factor = 1.0;

    if (Eflux > 0.0) {
        double net_E_lossl = std::abs(left_state.delta_conserved(4)- Eflux);
        double max_allowed_E_lossl = safe_fraction * left_state.get_conserved_total_energy();
        if (net_E_lossl > max_allowed_E_lossl) {
            axis_scale_factor = std::min(axis_scale_factor, max_allowed_E_lossl / net_E_lossl);
        }
    }
    if (mflux > 0.0) {
        double existing_loss = (left_state.delta_conserved(0) < 0.0) ? -left_state.delta_conserved(0) : 0.0;
        double net_M_loss_l = existing_loss + (mflux * dt);
        
        double max_allowed_M_lossl = safe_fraction * left_state.get_conserved_mass();
        if (net_M_loss_l > max_allowed_M_lossl) {
            axis_scale_factor = std::min(axis_scale_factor, max_allowed_M_lossl / (net_M_loss_l + 1e-20));
        }
    }
    if (EintFlux > 0.0) {
        double net_Eint_lossl = std::abs(left_state.delta_conserved(9) - EintFlux);
        double max_allowed_Eint_lossl = safe_fraction * left_state.get_conserved_internal_energy();
        if (net_Eint_lossl > max_allowed_Eint_lossl) {
            axis_scale_factor = std::min(axis_scale_factor, max_allowed_Eint_lossl / net_Eint_lossl);
        }
    }

    // --- RIGHT CELL INTEGRATED CHECKS ---
    // If Eflux is negative, energy is moving left (OUT of the right cell).
    if (Eflux < 0.0) {
        double net_E_lossr = std::abs(right_state.delta_conserved(4) + Eflux);
        double abs_Eflux_r = -net_E_lossr; // Convert to absolute loss
        double max_allowed_E_lossr = safe_fraction * right_state.get_conserved_total_energy();
        if (abs_Eflux_r > max_allowed_E_lossr) {
            axis_scale_factor = std::min(axis_scale_factor, max_allowed_E_lossr / abs_Eflux_r);
        }
    }
    if (mflux < 0.0) {
        double existing_loss = (right_state.delta_conserved(0) < 0.0) ? -right_state.delta_conserved(0) : 0.0;
        double net_M_loss_r = existing_loss + (-mflux * dt);
        
        double max_allowed_M_lossr = safe_fraction * right_state.get_conserved_mass();
        if (net_M_loss_r > max_allowed_M_lossr) {
            axis_scale_factor = std::min(axis_scale_factor, max_allowed_M_lossr / (net_M_loss_r + 1e-20));
        }
    }
    if (EintFlux < 0.0) {
        double net_Eint_lossr = std::abs(right_state.delta_conserved(9) + EintFlux);
        double abs_eintflux_r = -net_Eint_lossr;
        double max_allowed_Eint_lossr = safe_fraction * right_state.get_conserved_internal_energy();
        if (abs_eintflux_r > max_allowed_Eint_lossr) {
            axis_scale_factor = std::min(axis_scale_factor, max_allowed_Eint_lossr / abs_eintflux_r);
        }
    }

    mflux *= axis_scale_factor;
    pflux *= axis_scale_factor;
    Eflux *= axis_scale_factor;
    EintFlux *= axis_scale_factor;
    Bflux *= axis_scale_factor;
    B_scalarflux *= axis_scale_factor;

    left_state.delta_conserved(0) -= mflux;
    left_state.delta_conserved(1) -= pflux.x();
    left_state.delta_conserved(2) -= pflux.y();
    left_state.delta_conserved(3) -= pflux.z();
    left_state.delta_conserved(4) -= Eflux;
    left_state.delta_conserved(5) -= Bflux.x(); // mgb edit 21.04.2026 - update magnetic field flux in left state
    left_state.delta_conserved(6) -= Bflux.y(); 
    left_state.delta_conserved(7) -= Bflux.z();
    left_state.delta_conserved(8) -= B_scalarflux; // mgb edit 21.04.2026 - update B scalar flux in left state
    left_state.delta_conserved(9) -= EintFlux; // mgb edit 21.04.2026 - update internal energy flux in left state

    right_state.delta_conserved(0) += mflux;
    right_state.delta_conserved(1) += pflux.x();
    right_state.delta_conserved(2) += pflux.y();
    right_state.delta_conserved(3) += pflux.z();
    right_state.delta_conserved(4) += Eflux;
    right_state.delta_conserved(5) += Bflux.x(); // mgb edit 21.04.2026 - update magnetic field flux in right state
    right_state.delta_conserved(6) += Bflux.y();
    right_state.delta_conserved(7) += Bflux.z();
    right_state.delta_conserved(8) += B_scalarflux; // mgb edit 21.04.2026 - update B scalar flux in right state
    right_state.delta_conserved(9) += EintFlux; // mgb edit 21.04.2026 - update internal energy flux in right state

  }

  /**
   * @brief Do the flux calculation across a box boundary.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state hydro variables.
   * @param boundary HydroBoundary that sets the right state variables.
   * @param dx Distance between left and right state midpoint (in m).
   * @param A Surface area of the interface (in m^2).
   * @param dt Current system time step, used for flux limiter (in s).
   */
  inline void do_ghost_flux_calculation(const uint_fast8_t i,
                                        const CoordinateVector<> posR,
                                        MagnetoHydroVariables &left_state,
                                        const MagnetoHydroBoundary &boundary,
                                        const double dx, const double A,
                                        const double dt) const {

    // the sign bit is set (1) for negative values
    int_fast8_t orientation = 1 - 2 * std::signbit(dx);
    MagnetoHydroVariables right_state = boundary.get_right_state_flux_variables(
        i, orientation, posR, left_state);

    const double halfdx = 0.5 * dx;
    double rhoL = left_state.get_primitives_density() +
                  halfdx * left_state.primitive_gradients(0)[i];
    CoordinateVector<> vL(left_state.primitives(1) +
                              halfdx * left_state.primitive_gradients(1)[i],
                          left_state.primitives(2) +
                              halfdx * left_state.primitive_gradients(2)[i],
                          left_state.primitives(3) +
                              halfdx * left_state.primitive_gradients(3)[i]);
    double PL = left_state.get_primitives_pressure() +
                halfdx * left_state.primitive_gradients(4)[i];
    CoordinateVector<> BL(left_state.primitives(5) +
                              halfdx * left_state.primitive_gradients(5)[i],
                          left_state.primitives(6) +
                              halfdx * left_state.primitive_gradients(6)[i],
                          left_state.primitives(7) +
                              halfdx * left_state.primitive_gradients(7)[i]); // mgb edit 21.04.2026 - reconstructed magnetic field components
    double B_scalarL = left_state.primitives(8) + // mgb edit 21.04.2026 - reconstructed B scalar field
                       halfdx * left_state.primitive_gradients(8)[i];
    double rhoR = right_state.get_primitives_density() -
                  halfdx * right_state.primitive_gradients(0)[i];
    CoordinateVector<> vR(right_state.primitives(1) -
                              halfdx * right_state.primitive_gradients(1)[i],
                          right_state.primitives(2) -
                              halfdx * right_state.primitive_gradients(2)[i],
                          right_state.primitives(3) -
                              halfdx * right_state.primitive_gradients(3)[i]);
    double PR = right_state.get_primitives_pressure() -
                halfdx * right_state.primitive_gradients(4)[i];
    CoordinateVector<> BR(right_state.primitives(5) - // mgb edit 21.04.2026 start
                              halfdx * right_state.primitive_gradients(5)[i],
                          right_state.primitives(6) -
                              halfdx * right_state.primitive_gradients(6)[i],
                          right_state.primitives(7) -
                              halfdx * right_state.primitive_gradients(7)[i]); // mgb edit 21.04.2026 - reconstructed magnetic field components
    double B_scalarR = right_state.primitives(8) - // mgb edit 21.04.2026 - reconstructed B scalar field
                       halfdx * right_state.primitive_gradients(8)[i];

    rhoL = limit(rhoL, left_state.get_primitives_density(),
                 right_state.get_primitives_density(), 0.5);
    vL[0] = limit(vL.x(), left_state.get_primitives_velocity().x(),
                  right_state.get_primitives_velocity().x(), 0.5);
    vL[1] = limit(vL.y(), left_state.get_primitives_velocity().y(),
                  right_state.get_primitives_velocity().y(), 0.5);
    vL[2] = limit(vL.z(), left_state.get_primitives_velocity().z(),
                  right_state.get_primitives_velocity().z(), 0.5);
    PL = limit(PL, left_state.get_primitives_pressure(),
               right_state.get_primitives_pressure(), 0.5);
    BL[0] = limit(BL.x(), left_state.get_primitives_magnetic_field().x(),
                  right_state.get_primitives_magnetic_field().x(), 0.5); // mgb edit 21.04.2026 - limit reconstructed magnetic field components
    BL[1] = limit(BL.y(), left_state.get_primitives_magnetic_field().y(),
                  right_state.get_primitives_magnetic_field().y(), 0.5);
    BL[2] = limit(BL.z(), left_state.get_primitives_magnetic_field().z(),
                  right_state.get_primitives_magnetic_field().z(), 0.5);
    B_scalarL = limit(B_scalarL, left_state.get_primitives_B_scalar(),
                      right_state.get_primitives_B_scalar(), 0.5); // mgb edit 21.04.2026 - limit reconstructed B scalar field

    rhoR = limit(rhoR, right_state.get_primitives_density(),
                 left_state.get_primitives_density(), 0.5);
    vR[0] = limit(vR.x(), right_state.get_primitives_velocity().x(),
                  left_state.get_primitives_velocity().x(), 0.5);
    vR[1] = limit(vR.y(), right_state.get_primitives_velocity().y(),
                  left_state.get_primitives_velocity().y(), 0.5);
    vR[2] = limit(vR.z(), right_state.get_primitives_velocity().z(),
                  left_state.get_primitives_velocity().z(), 0.5);
    PR = limit(PR, right_state.get_primitives_pressure(),
               left_state.get_primitives_pressure(), 0.5);
    BR[0] = limit(BR.x(), right_state.get_primitives_magnetic_field().x(),
                  left_state.get_primitives_magnetic_field().x(), 0.5); // mgb edit 21.04.2026 - limit reconstructed magnetic field components
    BR[1] = limit(BR.y(), right_state.get_primitives_magnetic_field().y(),
                  left_state.get_primitives_magnetic_field().y(), 0.5);
    BR[2] = limit(BR.z(), right_state.get_primitives_magnetic_field().z(),
                  left_state.get_primitives_magnetic_field().z(), 0.5);
    B_scalarR = limit(B_scalarR, right_state.get_primitives_B_scalar(),
                      left_state.get_primitives_B_scalar(), 0.5); // mgb edit 21.04.2026 - limit reconstructed B scalar field

    cmac_assert(rhoL == rhoL);
    cmac_assert(vL.x() == vL.x());
    cmac_assert(vL.y() == vL.y());
    cmac_assert(vL.z() == vL.z());
    cmac_assert(PL == PL);
    cmac_assert(BL.x() == BL.x());
    cmac_assert(BL.y() == BL.y());
    cmac_assert(BL.z() == BL.z());
    cmac_assert(B_scalarL == B_scalarL);

    cmac_assert(rhoR == rhoR);
    cmac_assert(vR.x() == vR.x());
    cmac_assert(vR.y() == vR.y());
    cmac_assert(vR.z() == vR.z());
    cmac_assert(PR == PR);
    cmac_assert(BR.x() == BR.x());
    cmac_assert(BR.y() == BR.y());
    cmac_assert(BR.z() == BR.z());
    cmac_assert(B_scalarR == B_scalarR);

    // make sure all densities and pressures are physical
#ifdef SAFE_HYDRO_VARIABLES
    rhoL = std::max(rhoL, 0.);
    PL = std::max(PL, 0.);
    rhoR = std::max(rhoR, 0.);
    PR = std::max(PR, 0.);
#else
    cmac_assert(rhoL >= 0.);
    cmac_assert(PL >= 0.);
    cmac_assert(rhoR >= 0.);
    cmac_assert(PR >= 0.);
#endif

    double mflux = 0.;
    CoordinateVector<> pflux;
    double Eflux = 0.;
    double EintFlux = 0.; // mgb edit 21.04.2026 - internal energy flux
    CoordinateVector<> Bflux; // mgb edit 21.04.2026 - magnetic field flux
    double B_scalarflux = 0.; // mgb edit 21.04.2026 - B scalar flux
    CoordinateVector<> normal;
    double current_etot = left_state.get_conserved_total_energy();
    double current_eint = left_state.get_conserved_internal_energy();
    normal[i] = orientation;
    _riemann_solver.solve_for_flux_MHD(rhoL, vL, PL, BL, B_scalarL, rhoR, vR, PR, BR, B_scalarR, mflux, pflux,
                                   Eflux, EintFlux, Bflux, B_scalarflux, normal, current_etot, _mach_limit, current_eint);
  //  std::cout << "Fluxes calculated Line 828 MagnetoHydro ghost cells Bflux = " << Bflux.x() << ", " << Bflux.y() << ", " << Bflux.z() << std::endl; // mgb 04.05.2026

    cmac_assert(mflux == mflux);
    cmac_assert(pflux.x() == pflux.x());
    cmac_assert(pflux.y() == pflux.y());
    cmac_assert(pflux.z() == pflux.z());
    cmac_assert(Eflux == Eflux);
    cmac_assert(Bflux.x() == Bflux.x());
    cmac_assert(Bflux.y() == Bflux.y());
    cmac_assert(Bflux.z() == Bflux.z());
    cmac_assert(B_scalarflux == B_scalarflux);

    mflux *= A;
    pflux[0] *= A;
    pflux[1] *= A;
    pflux[2] *= A;
    Eflux *= A;
    EintFlux *= A; // mgb edit 21.04.2026 - scale internal energy flux by area
    Bflux[0] *= A; // mgb edit 21.04.2026 - scale magnetic field flux by area
    Bflux[1] *= A;
    Bflux[2] *= A;
    B_scalarflux *= A; // mgb edit 21.04.2026 - scale B scalar flux by area

#ifdef FLUX_LIMITER
    // limit the flux
    double fluxfac = 1.;
    const double absmflux = mflux * dt;
    if (absmflux > FLUX_LIMITER * left_state.get_conserved_mass()) {
      fluxfac = FLUX_LIMITER * left_state.get_conserved_mass() / absmflux;
    }
    if (_gamma > 1.) {
      const double absEflux = Eflux * dt;
      if (absEflux > FLUX_LIMITER * left_state.get_conserved_total_energy()) {
        fluxfac = std::min(
            fluxfac,
            FLUX_LIMITER * left_state.get_conserved_total_energy() / absEflux);
      }
    }
    // momentum flux limiter
    // note that we only apply this for cells that have high momentum, i.e.
    // whose momentum is higher than the thermal momentum of the cell
    // without this condition, cells with zero momentum would never be able
    // to gain momentum...
    const double p2 = left_state.get_conserved_momentum().norm2();
    const double m2 =
        left_state.get_conserved_mass() * left_state.get_conserved_mass();
    if (p2 * left_state.get_primitives_density() >
        _gamma * m2 * left_state.get_primitives_pressure()) {
      const double pflux2 = pflux.norm2() * dt;
      if (pflux2 > (FLUX_LIMITER * FLUX_LIMITER) * p2) {
        fluxfac = std::min(
            fluxfac, std::sqrt((FLUX_LIMITER * FLUX_LIMITER) * p2 / pflux2));
      }
    }
    if (std::isnan(pflux.x()) || std::isnan(pflux.y()) || std::isnan(pflux.z())) {
      std::cout << "Warning: Momentum flux is NaN before ghost cell flux limiter, pflux: " << pflux.x() << ", " << pflux.y() << ", " << pflux.z() << std::endl;
    }
    if (std::isnan(Bflux.x()) || std::isnan(Bflux.y()) || std::isnan(Bflux.z())) {
      std::cout << "Warning: Magnetic field flux is NaN before ghost cell flux limiter, Bflux: " << Bflux.x() << ", " << Bflux.y() << ", " << Bflux.z() << std::endl;
    }
    cmac_assert_message(fluxfac >= 0. && fluxfac <= 1., "fluxfac: %g", fluxfac);
    mflux *= fluxfac;
    pflux *= fluxfac;
    Eflux *= fluxfac;
    EintFlux *= fluxfac; // mgb edit 21.04.2026 - apply flux limiter to internal energy flux
    Bflux[0] *= fluxfac; // mgb edit 21.04.2026 - apply flux limiter to magnetic field flux
    Bflux[1] *= fluxfac;
    Bflux[2] *= fluxfac;
    B_scalarflux *= fluxfac; // mgb edit 21.04.2026 - apply flux limiter to B scalar flux
#endif

    if (std::isnan(pflux.x()) || std::isnan(pflux.y()) || std::isnan(pflux.z())) {
      std::cout << "Warning: Momentum flux is NaN after ghost cell flux limiter, pflux: " << pflux.x() << ", " << pflux.y() << ", " << pflux.z() << std::endl;
    }
    if (std::isnan(Bflux.x()) || std::isnan(Bflux.y()) || std::isnan(Bflux.z())) {
      std::cout << "Warning: Magnetic field flux is NaN after ghost cell flux limiter, Bflux: " << Bflux.x() << ", " << Bflux.y() << ", " << Bflux.z() << std::endl;
    }
    const double safe_fraction = 0.8;
    double axis_scale_factor = 1.0;


    if (Eflux > 0.0) {
        double max_allowed_E_lossl = safe_fraction * left_state.get_conserved_total_energy();
        if (Eflux > max_allowed_E_lossl) {
            axis_scale_factor = std::min(axis_scale_factor, max_allowed_E_lossl / Eflux);
        }
    }
    if (mflux > 0.0) {
        double max_allowed_M_lossl = safe_fraction * left_state.get_conserved_mass();
        if (mflux > max_allowed_M_lossl) {
            axis_scale_factor = std::min(axis_scale_factor, max_allowed_M_lossl / mflux);
        }
    }
    if (EintFlux > 0.0) {
        double max_allowed_Eint_lossl = safe_fraction * left_state.get_conserved_internal_energy();
        if (EintFlux > max_allowed_Eint_lossl) {
            axis_scale_factor = std::min(axis_scale_factor, max_allowed_Eint_lossl / EintFlux);
        }
    }

    // --- RIGHT CELL INTEGRATED CHECKS ---
    // If Eflux is negative, energy is moving left (OUT of the right cell).
    if (Eflux < 0.0) {
        double abs_Eflux_r = -Eflux; // Convert to absolute loss
        double max_allowed_E_lossr = safe_fraction * right_state.get_conserved_total_energy();
        if (abs_Eflux_r > max_allowed_E_lossr) {
            axis_scale_factor = std::min(axis_scale_factor, max_allowed_E_lossr / abs_Eflux_r);
        }
    }
    if (mflux < 0.0) {
        double abs_mflux_r = -mflux;
        double max_allowed_M_lossr = safe_fraction * right_state.get_conserved_mass();
        if (abs_mflux_r > max_allowed_M_lossr) {
            axis_scale_factor = std::min(axis_scale_factor, max_allowed_M_lossr / abs_mflux_r);
        }
    }
    if (EintFlux < 0.0) {
        double abs_eintflux_r = -EintFlux;
        double max_allowed_Eint_lossr = safe_fraction * right_state.get_conserved_internal_energy();
        if (abs_eintflux_r > max_allowed_Eint_lossr) {
            axis_scale_factor = std::min(axis_scale_factor, max_allowed_Eint_lossr / abs_eintflux_r);
        }
    }

    mflux *= axis_scale_factor;
    pflux *= axis_scale_factor;
    Eflux *= axis_scale_factor;
    EintFlux *= axis_scale_factor;
    Bflux *= axis_scale_factor;
    B_scalarflux *= axis_scale_factor;

    left_state.delta_conserved(0) -= mflux;
    left_state.delta_conserved(1) -= pflux.x();
    left_state.delta_conserved(2) -= pflux.y();
    left_state.delta_conserved(3) -= pflux.z();
    left_state.delta_conserved(4) -= Eflux;
    left_state.delta_conserved(5) -= Bflux.x(); // mgb edit 21.04.2026 - update magnetic field flux in left state
    left_state.delta_conserved(6) -= Bflux.y(); 
    left_state.delta_conserved(7) -= Bflux.z();
    left_state.delta_conserved(8) -= B_scalarflux; // mgb edit 21.04.2026 - update B scalar flux in left state
    left_state.delta_conserved(9) -= EintFlux; // mgb edit 21.04.2026 - update internal energy flux in left state

  
  }

  /**
   * @brief Do the gradient calculation for the given interface.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param left_state Left state variables.
   * @param right_state Right state variables.
   * @param dxinv Inverse distance between left and right state midpoint (in m).
   * @param WLlim Left state primitive variable limiters (updated; density -
   * kg m^-3, velocity - m s^-1, pressure - kg m^-1 s^-2).
   * @param WRlim Right state primitive variable limiters (updated; density -
   * kg m^-3, velocity - m s^-1, pressure - kg m^-1 s^-2).
   */
  inline void do_gradient_calculation(const int i, MagnetoHydroVariables &left_state,
                                      MagnetoHydroVariables &right_state,
                                      const double dxinv, double WLlim[20], // mgb edit 21.04.2026 - 10 -> 18
                                      double WRlim[20]) const { // mgb edit 21.04.2026 - 10 -> 18

    for (int_fast32_t j = 0; j < 10; ++j) { // mgb edit 21.04.2026 - 5 -> 10
      cmac_assert_message(left_state.primitives(j) == left_state.primitives(j),
                          "j: %" PRIiFAST32, j);
      cmac_assert_message(right_state.primitives(j) ==
                              right_state.primitives(j),
                          "j: %" PRIiFAST32, j);

      const double dwdx =
          0.5 * (right_state.primitives(j) - left_state.primitives(j)) * dxinv; // mgb edit 15.05.2026 (was +) 

      
      cmac_assert_message(
          dwdx == dwdx, "j: %" PRIiFAST32 ", left: %g, right: %g, dxinv: %g", j,
          left_state.primitives(j), right_state.primitives(j), dxinv);

      CoordinateVector<> gradL = left_state.primitive_gradients(j);
      CoordinateVector<> gradR = right_state.primitive_gradients(j);

      gradL[i] += dwdx;
      gradR[i] += dwdx; // mgb edit 04.08.2026 changed from - to + as gradient is true dw/dx now

      left_state.primitive_gradients(j)  = gradL[i];
      right_state.primitive_gradients(j) = gradR[i];


     // if (dwdx != 0){
     /*
      if (std::abs(left_state.primitive_gradients(j)[0]) >0 || std::abs(left_state.primitive_gradients(j)[1]) >0 ||std::abs(left_state.primitive_gradients(j)[2]) >0) {
        std::cout << "Warning: non zero gradient has been added! state gradients x[0]: "  << left_state.primitive_gradients(5)[0] << ", x[1]: " << left_state.primitive_gradients(5)[1] << ", x[2]: " << left_state.primitive_gradients(5)[1]<< std::endl;
        std::cout << "dwdx = " << dwdx << std::endl;
        std::cout<< "j: " << j << ", i: " << i << std::endl;
      }
        if (left_state.primitive_gradients(5)[0] ==0 || left_state.primitive_gradients(5)[1] ==0 || left_state.primitive_gradients(5)[2] ==0) {
          std::cout << "Warning: Primitive Gradients cross terms are not being set. x[0]: " << left_state.primitive_gradients(5)[0] << ", x[1]: " << left_state.primitive_gradients(5)[1] << ", x[2]: " << left_state.primitive_gradients(5)[1] << std::endl; 
          std::cout << "dwdx = " << dwdx << std::endl;
          std::cout<< "j: " << j << ", i: " << i << std::endl;
          
        }*/
      //}

      WLlim[2 * j] = std::min(WLlim[2 * j], right_state.primitives(j));
      WLlim[2 * j + 1] = std::max(WLlim[2 * j + 1], right_state.primitives(j));
      WRlim[2 * j] = std::min(WRlim[2 * j], left_state.primitives(j));
      WRlim[2 * j + 1] = std::max(WRlim[2 * j + 1], left_state.primitives(j));
    }
  }

  /**
   * @brief Do the gradient calculation across a box boundary.
   *
   * @param i Interface direction: x (0), y (1) or z (2).
   * @param posR Midpoint position of the ghost cell (in m).
   * @param left_state Left state variables.
   * @param boundary HydroBoundary that sets the right state variables.
   * @param dxinv Inverse distance between left and right state midpoint (in m).
   * @param WLlim Left state primitive variable limiters (updated; density -
   * kg m^-3, velocity - m s^-1, pressure - kg m^-1 s^-2).
   */
  inline void do_ghost_gradient_calculation(const int_fast32_t i,
                                            const CoordinateVector<> posR,
                                            MagnetoHydroVariables &left_state,
                                            const MagnetoHydroBoundary &boundary,
                                            const double dxinv,
                                            double WLlim[20]) const { // mgb edit 21.04.2026 - 10 -> 18

    // the sign bit is set (1) for negative values
    int_fast8_t orientation = 1 - 2 * std::signbit(dxinv);
    MagnetoHydroVariables right_state = boundary.get_right_state_gradient_variables(
        i, orientation, posR, left_state);
    for (int_fast32_t j = 0; j < 10; ++j) {  // mgb edit 21.04.2026 - 5 -> 10
      cmac_assert_message(left_state.primitives(j) == left_state.primitives(j),
                          "j: %" PRIiFAST32, j);
      cmac_assert_message(right_state.primitives(j) ==
                              right_state.primitives(j),
                          "j: %" PRIiFAST32, j);

      const double dwdx =
          0.5 * (left_state.primitives(j) + right_state.primitives(j)) * dxinv;

      cmac_assert_message(
          dwdx == dwdx, "j: %" PRIiFAST32 ", left: %g, right: %g, dxinv: %g", j,
          left_state.primitives(j), right_state.primitives(j), dxinv);

      left_state.primitive_gradients(j)[i] += dwdx;
      WLlim[2 * j] = std::min(WLlim[2 * j], right_state.primitives(j));
      WLlim[2 * j + 1] = std::max(WLlim[2 * j + 1], right_state.primitives(j));
    }
  }

  /**
   * @brief Apply the slope limiter for the given variables.
   *
   * @param state Hydro state.
   * @param Wlim Primitive variable limiters (density - kg m^-3, velocity -
   * m s^-1, pressure - kg m^-1 s^-2).
   * @param dx Distance between the cell and the neighbouring cells in all
   * directions (in m).
   */
  inline void apply_slope_limiter(MagnetoHydroVariables &state, const double Wlim[20], // mgb edit 21.04.2026 - 10 -> 18
                                  const CoordinateVector<> dx) const {

    for (int_fast8_t i = 0; i < 10; ++i) { // mgb edit 21.04.2026 - 5 -> 10
      cmac_assert_message(
          state.primitive_gradients(i)[0] == state.primitive_gradients(i)[0],
          "%" PRIiFAST8 ", (%g %g %g)", i, state.primitive_gradients(i)[0],
          state.primitive_gradients(i)[1], state.primitive_gradients(i)[2]);
      cmac_assert_message(
          state.primitive_gradients(i)[1] == state.primitive_gradients(i)[1],
          "%" PRIiFAST8 ", (%g %g %g)", i, state.primitive_gradients(i)[0],
          state.primitive_gradients(i)[1], state.primitive_gradients(i)[2]);
      cmac_assert_message(
          state.primitive_gradients(i)[2] == state.primitive_gradients(i)[2],
          "%" PRIiFAST8 ", (%g %g %g)", i, state.primitive_gradients(i)[0],
          state.primitive_gradients(i)[1], state.primitive_gradients(i)[2]);

      const double dwext[3] = {state.primitive_gradients(i)[0] * 0.5 * dx[0],
                               state.primitive_gradients(i)[1] * 0.5 * dx[1],
                               state.primitive_gradients(i)[2] * 0.5 * dx[2]};
      double dwmax = std::max(state.primitives(i) + dwext[0],
                              state.primitives(i) - dwext[0]);
      double dwmin = std::min(state.primitives(i) + dwext[0],
                              state.primitives(i) - dwext[0]);
      for (int_fast8_t j = 1; j < 3; ++j) {
        dwmax = std::max(dwmax, state.primitives(i) + dwext[j]);
        dwmin = std::min(dwmin, state.primitives(i) + dwext[j]);
        dwmax = std::max(dwmax, state.primitives(i) - dwext[j]);
        dwmin = std::min(dwmin, state.primitives(i) - dwext[j]);
      }
      dwmax -= state.primitives(i);
      dwmin -= state.primitives(i);
      double maxfac = DBL_MAX;
      if (dwmax != 0.) {
        const double dwngbmax = Wlim[2 * i + 1] - state.primitives(i);
        maxfac = dwngbmax / dwmax;
      }
      double minfac = DBL_MAX;
      if (dwmin != 0.) {
        const double dwngbmin = Wlim[2 * i] - state.primitives(i);
        minfac = dwngbmin / dwmin;
      }
      const double alpha = std::min(1., 0.5 * std::min(maxfac, minfac));
      state.primitive_gradients(i) *= alpha;

      cmac_assert_message(
          state.primitive_gradients(i)[0] == state.primitive_gradients(i)[0],
          "%" PRIiFAST8 ", (%g %g %g), alpha: %g", i,
          state.primitive_gradients(i)[0], state.primitive_gradients(i)[1],
          state.primitive_gradients(i)[2], alpha);
      cmac_assert_message(
          state.primitive_gradients(i)[1] == state.primitive_gradients(i)[1],
          "%" PRIiFAST8 ", (%g %g %g), alpha: %g", i,
          state.primitive_gradients(i)[0], state.primitive_gradients(i)[1],
          state.primitive_gradients(i)[2], alpha);
      cmac_assert_message(
          state.primitive_gradients(i)[2] == state.primitive_gradients(i)[2],
          "%" PRIiFAST8 ", (%g %g %g), alpha: %g", i,
          state.primitive_gradients(i)[0], state.primitive_gradients(i)[1],
          state.primitive_gradients(i)[2], alpha);
    }
  }

  /**
   * @brief Predict the primitive variables forward in time with the given time
   * step.
   *
   * @param state Hydro variables of the cell.
   * @param dt Time step (in s).
   */
  inline void predict_primitive_variables(MagnetoHydroVariables &state,
                                          const double dt) const {

    const double rho = state.get_primitives_density();

    if (rho == 0.) {
      return;
    }
    const double rhoinv = 1. / rho;
    if (std::isinf(rhoinv)) {
      return;
    }

    const CoordinateVector<> B = state.get_primitives_magnetic_field();

    const double vx = state.get_primitives_velocity().x();
    const double vy = state.get_primitives_velocity().y();
    const double vz = state.get_primitives_velocity().z();
    const double P = state.get_primitives_pressure();
    const double Bx = state.get_primitives_magnetic_field().x(); // mgb edit 21.04.2026 - magnetic field components
    const double By = state.get_primitives_magnetic_field().y();
    const double Bz = state.get_primitives_magnetic_field().z();
    const double B_scalar = state.get_primitives_B_scalar(); // mgb edit 21.04.2026 - B scalar field
    const double ax = state.get_gravitational_acceleration().x();
    const double ay = state.get_gravitational_acceleration().y();
    const double az = state.get_gravitational_acceleration().z();

    const double drhodx = state.primitive_gradients(0).x();
    const double drhody = state.primitive_gradients(0).y();
    const double drhodz = state.primitive_gradients(0).z();

    const double dvxdx = state.primitive_gradients(1).x();
    const double dvxdy = state.primitive_gradients(1).y();
    const double dvxdz = state.primitive_gradients(1).z();
    const double dvydx = state.primitive_gradients(2).x();
    const double dvydy = state.primitive_gradients(2).y();
    const double dvydz = state.primitive_gradients(2).z();
    const double dvzdx = state.primitive_gradients(3).x();
    const double dvzdy = state.primitive_gradients(3).y();
    const double dvzdz = state.primitive_gradients(3).z();

    const double dPdx = state.primitive_gradients(4).x();
    const double dPdy = state.primitive_gradients(4).y();
    const double dPdz = state.primitive_gradients(4).z();

    const double dBxdx = state.primitive_gradients(5).x(); // mgb edit 21.04.2026 - magnetic field gradients
    const double dBydy = state.primitive_gradients(6).y();
    const double dBzdz = state.primitive_gradients(7).z();
    const double dBxdy = state.primitive_gradients(5).y();
    const double dBxdz = state.primitive_gradients(5).z();
    const double dBydx = state.primitive_gradients(6).x();
    const double dBydz = state.primitive_gradients(6).z();
    const double dBzdx = state.primitive_gradients(7).x();
    const double dBzdy = state.primitive_gradients(7).y();

    const double dB_scalardx = state.primitive_gradients(8).x(); // mgb edit 21.04.2026 - B scalar field gradient
    const double dB_scalardy = state.primitive_gradients(8).y();
    const double dB_scalardz = state.primitive_gradients(8).z();

    const double dBsdx = dB_scalardx * Bx/_mu0;
    const double dBsdy = dB_scalardy * By/_mu0;
    const double dBsdz = dB_scalardz * Bz/_mu0;

    const double Tx = (Bx * dBxdx + By * dBxdy + Bz * dBxdz) / _mu0; // mgb edit 21.04.2026 - magnetic tension terms = B . grad B / mu0
    const double Ty = (Bx * dBydx + By * dBydy + Bz * dBydz) / _mu0; //. mgb comment 27.05.2026: to be consistent with vx * divv etc - don't include cross terms?
    const double Tz = (Bx * dBzdx + By * dBzdy + Bz * dBzdz) / _mu0;

    const double dPmagdx = (Bx * dBxdx + By * dBydx + Bz * dBzdx)/ _mu0; // mgb edit 21.04.2026 - magnetic pressure gradient terms = grad (B^2 / 2 mu0) (= (B . grad B) / mu0 for ideal MHD where div B = 0)
    const double dPmagdy = (Bx * dBxdy + By * dBydy + Bz * dBzdy) / _mu0;
    const double dPmagdz = (Bx * dBxdz + By * dBydz + Bz * dBzdz) / _mu0;

    const double divv = dvxdx + dvydy + dvzdz;
    const double divB = dBxdx + dBydy + dBzdz;


   // double rho_new =
     //   rho - dt * (rho * divv + vx * drhodx + vy * drhody + vz * drhodz);
    
    double lnRho = std::log(rho);

    double dlnRhodx = drhodx / rho;
    double dlnRhody = drhody / rho;
    double dlnRhodz = drhodz / rho;

    double lnRho_new = lnRho - dt * (divv + vx * dlnRhodx + vy * dlnRhody + vz * dlnRhodz);

    double rho_new = std::exp(lnRho_new);

    const double vxdiv = vx * dvxdx + vy * dvxdy + vz * dvxdz; //vx * dvxdx + vy * dvxdy + vz * dvxdz;
    const double vydiv = vx * dvydx + vy * dvydy + vz * dvydz;  //vx * dvydx + vy * dvydy + vz * dvydz;
    const double vzdiv = vx * dvzdx + vy * dvzdy + vz * dvzdz; //vx * dvzdx + vy * dvzdy + vz * dvzdz; 

    const CoordinateVector<> v_new(vx - dt * (vxdiv + rhoinv * (dPdx + dPmagdx - Tx + dBsdx) - ax), // mgb edit 21.04.2026 - include magnetic pressure gradient and tension terms in velocity update
                                   vy - dt * (vydiv + rhoinv * (dPdy + dPmagdy - Ty + dBsdy) - ay),
                                   vz - dt * (vzdiv + rhoinv * (dPdz + dPmagdz - Tz + dBsdz) - az));
    double P_new =
        P - dt * (_gamma * P * divv + vx * dPdx + vy * dPdy + vz * dPdz);

   /* if (dPmagdx - Tx != 0){
      std::cout<< "Warning: magnetic tension and pressure are equal and gas dynamics are not affected by B field! Tx: " << Tx << ", dPmagdx: " << dPmagdx << std::endl;
    }*/
        
    // mgb comment 27.05.2026: might need
    //double lnP = std::log(P);
    //double lnP_new = lnP - dt * (_gamma * divv + vx * dlnPdx + vy * dlnPdy + vz * dlnPdz);

    //double P_new = std::exp(lnP_new);

    const double cf = get_fast_magnetosonic_speed(rho, P, B, _max_velocity);

    CoordinateVector<> normal;
    normal[0] = 1.;

    const double v_new_norm = std::sqrt(vx*vx + vy*vy + vz*vz);//CoordinateVector<>::dot_product(v_new, normal); // mgb note 21.05.2026: picks out the x-velocity
    
    const double cleaning_speed_max = v_new_norm + cf;

    const double cleaning_speed = std::min(cleaning_speed_max, _max_velocity * _cleaning_speed_factor); // mgb edit 21.04.2026 - cleaning speed for divergence cleaning, set to maximum velocity in the simulation
    double B_scalar_new = B_scalar - dt * (cleaning_speed * divB); // mgb edit 21.04.2026 - B scalar field update with divergence cleaning term

     const CoordinateVector<> B_new(Bx - dt * (vy * dBxdy + vz * dBxdz - Bx * (dvydy + dvzdz) + By * dvxdy + Bz * dvxdz + cleaning_speed * dB_scalardx),
        By - dt * (vx * dBydx + vz * dBydz - By * (dvxdx + dvzdz) + Bx * dvydx + Bz * dvydz + cleaning_speed * dB_scalardy),
        Bz - dt * (vx * dBzdx + vy * dBzdy - Bz * (dvxdx + dvydy) + Bx * dvzdx + By * dvzdy + cleaning_speed * dB_scalardz));


    /*
    // Internal energy update mgb edit 15.05.2026
    const double Eint = state.get_primitives_internal_energy();
    //const double dEintdx = state.primitive_gradients(9).x();
    //const double dEintdy = state.primitive_gradients(9).y();
    //const double dEintdz = state.primitive_gradients(9).z();

    // mgb comment 21.05.2026: may want this?
   // double Eint_new = 0.;
    //if (Eint == 0) {
     double Eint_new = P_new/(_gamma_minus_one * rho_new);
    //} else{
      // Eint_new = Eint - dt * (_gamma_minus_one * Eint * divv +  vx * dEintdx + vy * dEintdy + vz * dEintdz);
   // }
   //double Eint_new = P_new/(_gamma_minus_one * rho_new);
   if (Eint_new<=0){
      std::cout <<  "Warning: predict_primitive_variables, Eint_new: " << Eint_new <<", Eint: " << Eint << ", P: " << P_new << ", rho: " << rho_new  << ", Eint from P: " << P_new/(_gamma_minus_one * rho_new) << std::endl;
    }*/

    double lnP = std::log(P);
    
    double dlnPdx = dPdx / P;
    double dlnPdy = dPdy / P;
    double dlnPdz = dPdz / P;

    double lnP_new = lnP - dt * (_gamma * divv + vx * dlnPdx + vy * dlnPdy + vz * dlnPdz);

     P_new = std::exp(lnP_new);

    double Eint_new = P_new / (_gamma_minus_one * rho_new);

    cmac_assert_message(rho_new == rho_new,
                        "rho: %g, divv: %g, v: %g %g %g, drho: %g %g %g", rho,
                        divv, vx, vy, vz, drhodx, drhody, drhodz);
    cmac_assert_message(
        v_new.x() == v_new.x(),
        "v: %g %g %g, divv: %g, rhoinv: %g, dP: %g %g %g, a: %g %g %g", vx, vy,
        vz, divv, rhoinv, dPdx, dPdy, dPdz, ax, ay, az);
    cmac_assert_message(
        v_new.y() == v_new.y(),
        "v: %g %g %g, divv: %g, rhoinv: %g, dP: %g %g %g, a: %g %g %g", vx, vy,
        vz, divv, rhoinv, dPdx, dPdy, dPdz, ax, ay, az);
    cmac_assert_message(
        v_new.z() == v_new.z(),
        "v: %g %g %g, divv: %g, rhoinv: %g, dP: %g %g %g, a: %g %g %g", vx, vy,
        vz, divv, rhoinv, dPdx, dPdy, dPdz, ax, ay, az);
    cmac_assert_message(P_new == P_new,
                        "P: %g, divv: %g, v: %g %g %g, dP: %g %g %g", P, divv,
                        vx, vy, vz, dPdx, dPdy, dPdz);

#ifdef SAFE_HYDRO_VARIABLES
    rho_new = std::max(rho_new, 0.);
    P_new = std::max(P_new, 0.);
    Eint_new = std::max(Eint_new, 0.);
#else
    cmac_assert(rho_new >= 0.);
    cmac_assert(P_new >= 0.);
    cmac_assert(Eint_new >= 0.);
#endif

    

    state.set_primitives_density(rho_new);
    state.set_primitives_velocity(v_new);
    state.set_primitives_pressure(P_new);
    state.set_primitives_magnetic_field(B_new); // mgb edit 21.04.2026 - update magnetic field primitives
    state.set_primitives_B_scalar(B_scalar_new); // mgb edit 21.04.2026 - update B scalar primitive
    state.set_primitives_internal_energy(Eint_new);
  //  std::cout << "Primitive B field set Line 1168 " << std::endl;
//    if (std::abs(B_new.x()) > 100 * std::abs(Bx)) {
     // std::cout << "Warning: Large increase in Bx detected: " << std::setprecision(15)  << Bx << " -> " << B_new.x() << std::endl;
   // }
   //  if (std::abs(B_new.y()) > 100 * std::abs(By)) {
   //   std::cout << "Warning: Large increase in By detected: " << std::setprecision(15) << By << " -> " << B_new.y() << std::endl;
   // }
  //   if (std::abs(B_new.z()) > 100 * std::abs(Bz)) {
   //   std::cout << "Warning: Large increase in Bz detected: " << std::setprecision(15)  << Bz << " -> " << B_new.z() << std::endl;
   // }

   // can uncomment probably mgb note 02.06.2026
   // if (std::abs(B_scalar_new) > 1e-5) {
    //    std::cout << "Warning: Large B_scalar detected: " << std::setprecision(15) << B_scalar << " -> " << B_scalar_new << std::endl;
    //  }
    if (state.get_conserved_total_energy() < state.get_conserved_internal_energy()) {
      std::cout << "Warning: Etot < Eint within predict_primitive_variables!"<<std::endl;
    }

  }

  /**
   * @brief Set the hydrodynamic variables based on the given ionization
   * variables.
   *
   * @param ionization_variables IonizationVariables.
   * @param magnetohydro_variables HydroVariables.
   */
  inline void
  ionization_to_hydro(const IonizationVariables &ionization_variables,
                      MagnetoHydroVariables &magnetohydro_variables) const { // mgb comment 26.05.2026: this is only called within initialize hydro vars

    double density =
        _density_conversion_factor * ionization_variables.get_number_density();
  //  const double mean_molecular_mass =
  //      0.5 * (1. + ionization_variables.get_ionic_fraction(ION_H_n));
    const double h0 = ionization_variables.get_ionic_fraction(ION_H_n);
#ifdef HAS_HELIUM
    const double hep = ionization_variables.get_ionic_fraction(ION_He_p1);
    const double hepp = 1 - hep - ionization_variables.get_ionic_fraction(ION_He_n);
    const double AHe = _abundances.get_abundance(ELEMENT_He);
    const double mean_molecular_mass = (1.0 + 4.0*AHe)/(2 + AHe - h0 + AHe*hep + 2.0*AHe*hepp);
#else
    const double mean_molecular_mass = 1.0/(2.0-h0);
#endif
    double pressure = _pressure_conversion_factor * density *
                      ionization_variables.get_temperature() /
                      mean_molecular_mass;

    // apply velocity limiter
    if (density > 0.) {
      const double inverse_density = 1. / density;
      if (!std::isinf(inverse_density)) {
        const double cs = std::sqrt(_gamma * pressure * inverse_density);
        const double B2 = magnetohydro_variables.get_primitives_magnetic_field().norm2();
        const double va_squared = B2/(_mu0 * density); // mgb edit 21.04.2026 - Alfvén speed squared = B^2 / (mu0 * rho)
        const double vMHD = std::sqrt(cs*cs + va_squared); // mgb edit 21.04.2026 - magnetohydrodynamic wave speed = sqrt(sound speed squared + Alfvén speed squared)
        if (vMHD > _max_velocity) {
          const double factor = _max_velocity / vMHD;
          pressure *= factor * factor;
        }
      }
    }


    cmac_assert_message(pressure > 0.0,"p=%g,t=%g,mmm=%g",pressure, ionization_variables.get_temperature(),mean_molecular_mass);

    cmac_assert(density == density);
    cmac_assert(pressure == pressure);

#ifdef SAFE_HYDRO_VARIABLES
    density = std::max(density, 0.);
    pressure = std::max(pressure, 0.);
#else
    cmac_assert(density >= 0.);
    cmac_assert(pressure >= 0.);
#endif

    // the velocity is directly set from the initial condition
    magnetohydro_variables.set_primitives_density(density);
    magnetohydro_variables.set_primitives_pressure(pressure);
    //magnetohydro_variables.set_primitives_magnetic_field(CoordinateVector<> (_init_Bx, _init_By, _init_Bz)); // mgb edit 21.04.2026 - set initial magnetic field from problem definition unsure if this is needed!
   // magnetohydro_variables.set_primitives_B_scalar(_init_B_scalar); // mgb edit 21.04.2026 - set initial B scalar field from problem definition unsure if this is needed!
    if (std::isnan(magnetohydro_variables.get_conserved_total_energy())) {
        std::cout << "Warning: Total energy is NaN before ionization_to_hydro conversion for pressure:" << pressure << ", density:" << density << std::endl;
    }
    if (magnetohydro_variables.get_conserved_total_energy() < magnetohydro_variables.get_conserved_internal_energy()){
        std::cout<<"Warning! Etot < Eint before update in ionization to hydro! Etot: " << magnetohydro_variables.get_conserved_total_energy() << ", Eint: " << magnetohydro_variables.get_conserved_internal_energy() << std::endl;
      }
    set_conserved_variables(magnetohydro_variables, magnetohydro_variables.get_conserved_mass()/density);

    if (magnetohydro_variables.get_conserved_total_energy() < magnetohydro_variables.get_conserved_internal_energy()){
        std::cout<<"Warning! Etot < Eint after update in ionization to hydro! Etot: " << magnetohydro_variables.get_conserved_total_energy() << ", Eint: " << magnetohydro_variables.get_conserved_internal_energy() << std::endl;
      } 
    if (std::isnan(magnetohydro_variables.get_conserved_total_energy())) {
        std::cout << "Warning: Total energy is NaN after ionization_to_hydro conversion for pressure:" << pressure << ", density:" << density << std::endl;
    }
  }

  /**
   * @brief Set all temperature related variables to the given temperature.
   *
   * @param ionization_variables IonizationVariables of the cell.
   * @param magnetohydro_variables MagnetoHydroVariables of the cell.
   * @param volume Volume of the cell (in m^3).
   * @param temperature Temperature value to enforce (in K).
   */
  inline void set_temperature(IonizationVariables &ionization_variables,
                              MagnetoHydroVariables &magnetohydro_variables,
                              const double volume,
                              const double temperature) const {

    const double density = magnetohydro_variables.get_primitives_density();
    if (density > 0.) {
      // const double xH = ionization_variables.get_ionic_fraction(ION_H_n);
      // const double thermal_energy =
      //     _u_conversion_factor * temperature / (1. + xH);
    const double h0 = ionization_variables.get_ionic_fraction(ION_H_n);
#ifdef HAS_HELIUM
    const double hep = ionization_variables.get_ionic_fraction(ION_He_p1);
    const double hepp = 1 - hep - ionization_variables.get_ionic_fraction(ION_He_n);
    const double AHe = _abundances.get_abundance(ELEMENT_He);
    const double mean_molecular_mass = (1.0 + 4.0*AHe)/(2 + AHe - h0 + AHe*hep + 2.0*AHe*hepp);
#else
    const double mean_molecular_mass = 1.0/(2.0-h0);
#endif
      const double thermal_energy = _u_conversion_factor*temperature/mean_molecular_mass/2.0;
      const double internal_energy = volume * density * thermal_energy;

      const double kinetic_energy =
          0.5 * CoordinateVector<>::dot_product(
                    magnetohydro_variables.get_primitives_velocity(),
                    magnetohydro_variables.get_conserved_momentum());
      const double magnetic_energy_density = (CoordinateVector<>::dot_product(magnetohydro_variables.get_primitives_magnetic_field(), magnetohydro_variables.get_primitives_magnetic_field()))/(2.0*_mu0); // mgb edit 21.04.2026 - magnetic energy density = B^2 / (2 mu0)
      
      const double total_energy =
          internal_energy + kinetic_energy + volume * magnetic_energy_density;
      const double pressure = _gamma_minus_one * density * thermal_energy;

      if ((magnetohydro_variables.get_conserved_total_energy() < magnetohydro_variables.get_conserved_internal_energy()) || (total_energy < internal_energy)){
        std::cout << "Warning! Etot < Eint before update in set_temperature! Etot: " << magnetohydro_variables.get_conserved_total_energy() << ", Eint: " << magnetohydro_variables.get_conserved_internal_energy() << ", total_energy: " << total_energy << ", internal_energy: " << internal_energy << std::endl;
      }
      if (internal_energy < 0){
        std::cout << "Warning: set_temperature Eint < 0: " << internal_energy << std::endl;
      }
      ionization_variables.set_temperature(temperature);
      magnetohydro_variables.set_primitives_pressure(pressure);
      magnetohydro_variables.set_conserved_total_energy(total_energy);
      magnetohydro_variables.set_conserved_internal_energy(internal_energy);
      
      if ((magnetohydro_variables.get_conserved_total_energy() < magnetohydro_variables.get_conserved_internal_energy()) || (total_energy < internal_energy)){
        std::cout << "Warning! Etot < Eint after update in set_temperature! Etot: " << magnetohydro_variables.get_conserved_total_energy() << ", Eint: " << magnetohydro_variables.get_conserved_internal_energy() << ", total_energy: " << total_energy << ", internal_energy: " << internal_energy << std::endl;
      }

    }
  }

  /**
   * @brief Get the temperature difference corresponding to the given change in
   * total energy.
   *
   * We assume that the change in total energy is entirely due to a change in
   * thermal energy.
   *
   * @param ionization_variables IonizationVariables of the cell.
   * @param magnetohydro_variables MagnetoHydroVariables of the cell.
   * @param inverse_volume Inverse volume of the cell (in m^-3).
   * @param delta_energy Change in energy (in J).
   * @return Corresponding temperature difference (in K).
   */
  inline double
  get_temperature_difference(const IonizationVariables &ionization_variables,
                             const MagnetoHydroVariables &magnetohydro_variables,
                             const double inverse_volume,
                             const double delta_energy) const {

    const double density = magnetohydro_variables.get_primitives_density();
    if (density == 0.) {
      return 0.;
    }
    const double inverse_density = 1. / density;
    if (std::isinf(inverse_density)) {
      return 0.;
    }

  //  const double mean_molecular_weight =
//        0.5 * (1. + ionization_variables.get_ionic_fraction(ION_H_n));
    const double h0 = ionization_variables.get_ionic_fraction(ION_H_n);
#ifdef HAS_HELIUM
    const double hep = ionization_variables.get_ionic_fraction(ION_He_p1);
    const double hepp = 1 - hep - ionization_variables.get_ionic_fraction(ION_He_n);
    const double AHe = _abundances.get_abundance(ELEMENT_He);
    const double mean_molecular_weight = (1.0 + 4.0*AHe)/(2 + AHe - h0 + AHe*hep + 2.0*AHe*hepp);
#else
    const double mean_molecular_weight = 1.0/(2.0-h0);
#endif
    return _T_conversion_factor * mean_molecular_weight * _gamma_minus_one *
           inverse_density * inverse_volume * delta_energy;
  }

  /**
   * @brief Update all variables that depend on the energy after a change in
   * total energy of the given cell by the given amount.
   *
   * @param ionization_variables IonizationVariables of the cell.
   * @param magnetohydro_variables MagnetoHydroVariables of the cell.
   * @param inverse_volume Inverse volume of the cell (in m^-3).
   * @param delta_energy Change in energy (in J).
   */
  inline void update_energy_variables(IonizationVariables &ionization_variables,
                                      MagnetoHydroVariables &magnetohydro_variables,
                                      const double inverse_volume,
                                      const double delta_energy) const {

    const double density = magnetohydro_variables.get_primitives_density();
    if (density == 0.) {
      return;
    }
    const double inverse_density = 1. / density;
    if (std::isinf(inverse_density)) {
      return;
    }

    const double old_energy = magnetohydro_variables.get_conserved_total_energy();
    const double old_eint = magnetohydro_variables.get_conserved_internal_energy();
    const double kinetic_energy =
        0.5 * CoordinateVector<>::dot_product(
                  magnetohydro_variables.get_primitives_velocity(),
                  magnetohydro_variables.get_conserved_momentum());
    const double magnetic_energy_density = (CoordinateVector<>::dot_product(magnetohydro_variables.get_primitives_magnetic_field(), magnetohydro_variables.get_primitives_magnetic_field()))/(2.0*_mu0); // mgb edit 21.04.2026 - magnetic energy density = B^2 / (2 mu0)
    //comment out berts version, trying mine with helium for now.... see what happens
    // const double mean_molecular_mass =
    //     0.5 * (1. + ionization_variables.get_ionic_fraction(ION_H_n));


    const double h0 = ionization_variables.get_ionic_fraction(ION_H_n);
#ifdef HAS_HELIUM
    const double hep = ionization_variables.get_ionic_fraction(ION_He_p1);
    const double hepp = 1 - hep - ionization_variables.get_ionic_fraction(ION_He_n);
    const double AHe = _abundances.get_abundance(ELEMENT_He);
    const double mean_molecular_mass = (1.0 + 4.0*AHe)/(2 + AHe - h0 + AHe*hep + 2.0*AHe*hepp);
#else
    const double mean_molecular_mass = 1.0/(2.0-h0);
#endif


    

    double new_energy = old_energy + delta_energy;
    double new_eint = old_eint + delta_energy;
    double internal_energy;

    if (new_eint < 0.0){
      std::cout << "Warning: update_energy_variables, Eint < 0: " << new_eint << ", dE: " << delta_energy << ", old Eint:" << old_eint << std::endl;
      std::cout << "Warning: update_energy_variables, Etot: " << new_energy << ", dE: " << delta_energy << ", old Etot:" << old_energy << std::endl;
      std::cout << "Warning: update_energy_variables, Ekin: " << kinetic_energy << ", Emag: " << magnetic_energy_density / inverse_volume << ", Eint2:" << old_energy-kinetic_energy-(magnetic_energy_density / inverse_volume) << std::endl;
    }

    const CoordinateVector<> velocity = magnetohydro_variables.get_primitives_velocity();
    double test_pressure = magnetohydro_variables.get_primitives_pressure(); // mgb edit 06.05.2026 - calculate pressure from internal energy density
    double cs_squared = _gamma * test_pressure / density; // mgb edit 06.05.2026 - calculate sound speed squared from test pressure and density
    double va_squared = (CoordinateVector<>::dot_product(magnetohydro_variables.get_primitives_magnetic_field(), magnetohydro_variables.get_primitives_magnetic_field())) / (_mu0 * density); // mgb edit 06.05.2026 - calculate Alfvén speed squared from magnetic field strength and density
    double v_fast = std::sqrt(cs_squared + va_squared);
    double v_mag = velocity.norm();

    double mach_fast = v_mag/(v_fast + 1.e-15); // mgb edit 06.05.2026 - calculate fast magnetosonic Mach number (add small number to avoid division by zero

    if (new_eint > 0.){
      if (mach_fast > _mach_limit || (kinetic_energy + magnetic_energy_density * 1/inverse_volume > 0.95 * new_energy) || (test_pressure-magnetohydro_variables.get_primitives_pressure())/(std::min(test_pressure, magnetohydro_variables.get_primitives_pressure())+1.e-10) > 0.25) {
      // std::cout << "Warning: dual energy triggered, emag:" << magnetic_energy_density * 1/inverse_volume << ", ekin:" << kinetic_energy << ", total energy:" << magnetohydro_variables.get_conserved_total_energy() << std::endl;
        internal_energy = new_eint;
      // std::cout << "Warning: dual energy triggered, Bx:" << magnetohydro_variables.get_primitives_magnetic_field().x() << ", By:" << magnetohydro_variables.get_primitives_magnetic_field().y() << ", Bz:" << magnetohydro_variables.get_primitives_magnetic_field().z() << ", 1/inverse volume:" << 1./inverse_volume << std::endl;
      // std::cout << "Warning: dual energy triggered, emag:" << magnetic_energy_density * 1/inverse_volume << ", ekin:" << kinetic_energy << ", eint:" << internal_energy << ", total energy:" << magnetohydro_variables.get_conserved_total_energy() << std::endl;
        if (std::isnan(internal_energy)) {
          std::cout << "NaN detected in internal energy tracker! old_energy: " << old_eint << " delta_energy: " << delta_energy << " new_energy: " << new_eint << std::endl;
          //exit(1);
        }

        new_energy = kinetic_energy + magnetic_energy_density/inverse_volume + internal_energy; // mgb edit 010.05.2026: update energy with new internal energy.
        if (new_energy < 0.){
          std::cout<<"Warning: New Etot < 0! Etot: " << new_energy << ", Emag: " << magnetic_energy_density/inverse_volume << ", Ekin: " << kinetic_energy << ", Eint: " << internal_energy << std::endl;
          new_energy = std::max(new_energy, 0.);
        }
        
        magnetohydro_variables.set_conserved_total_energy(new_energy);
        magnetohydro_variables.set_conserved_internal_energy(internal_energy);
        magnetohydro_variables.set_primitives_internal_energy(internal_energy * inverse_volume);

        if (internal_energy > magnetohydro_variables.get_conserved_total_energy()) {
          std::cout << "Warning: internal energy is GREATER than total energy in update_energy_variables dual energy! Eint: " << internal_energy << ", Etot: " << new_energy << ", Emag: " << magnetic_energy_density/inverse_volume << ", Ekin: " << kinetic_energy << std::endl;
        }
        if (internal_energy < 0) {
          std::cout << "Warning: update_energy_variables dual energy, Eint < 0: " << internal_energy << std::endl;
        }
        
        
      } else {
        internal_energy = (new_energy - kinetic_energy - magnetic_energy_density/inverse_volume);
        if (std::isnan(internal_energy)) {
          std::cout << "NaN detected in internal energy calculation! new_energy: " << new_energy  << ", old_energy: " << old_energy << " delta_energy: " << delta_energy << " kinetic_energy: " << kinetic_energy << " magnetic_energy_density: " << magnetic_energy_density/inverse_volume << std::endl;
          //exit(1);
        }
        if (internal_energy < 0) {
          std::cout << "Warning: update_energy_variables no dual energy, Eint < 0: " << internal_energy << std::endl;
        }
        new_energy = std::max(new_energy, 0.);
        internal_energy = std::max(internal_energy, 0.);
        magnetohydro_variables.set_conserved_total_energy(new_energy);
        magnetohydro_variables.set_conserved_internal_energy(internal_energy);
        magnetohydro_variables.set_primitives_internal_energy(internal_energy * inverse_volume);
    //   std::cout << "Warning: dual energy not triggered, Bx:" << magnetohydro_variables.get_primitives_magnetic_field().x() << ", By:" << magnetohydro_variables.get_primitives_magnetic_field().y() << ", Bz:" << magnetohydro_variables.get_primitives_magnetic_field().z() << ", 1/inverse volume:" << 1./inverse_volume << std::endl;
      //  std::cout << "Warning: dual energy not triggered, emag:" << magnetic_energy_density * 1/inverse_volume << ", ekin:" << kinetic_energy << ", eint:" << internal_energy << ", total energy:" << magnetohydro_variables.get_conserved_total_energy() << std::endl;
        if (internal_energy > magnetohydro_variables.get_conserved_total_energy()) {
          std::cout << "Warning: internal energy is GREATER than total energy in update_energy_variables! Eint: " << internal_energy << ", Etot: " << new_energy << ", Emag: " << magnetic_energy_density/inverse_volume << ", Ekin: " << kinetic_energy << std::endl;
        }
    }
  } else {
    internal_energy = (new_energy - kinetic_energy - magnetic_energy_density/inverse_volume);
    if (std::isnan(internal_energy)) {
      std::cout << "NaN detected in internal energy calculation! new_energy: " << new_energy  << ", old_energy: " << old_energy << " delta_energy: " << delta_energy << " kinetic_energy: " << kinetic_energy << " magnetic_energy_density: " << magnetic_energy_density/inverse_volume << std::endl;
      //exit(1);
    }
    if (internal_energy < 0) {
      std::cout << "Warning: update_energy_variables, no dual as internal tracker <0, Eint < 0: " << internal_energy << std::endl;
    }
    new_energy = std::max(new_energy, 0.);
    internal_energy = std::max(internal_energy, 0.);
    magnetohydro_variables.set_conserved_total_energy(new_energy);
    magnetohydro_variables.set_conserved_internal_energy(internal_energy);
    magnetohydro_variables.set_primitives_internal_energy(internal_energy * inverse_volume);
//   std::cout << "Warning: dual energy not triggered, Bx:" << magnetohydro_variables.get_primitives_magnetic_field().x() << ", By:" << magnetohydro_variables.get_primitives_magnetic_field().y() << ", Bz:" << magnetohydro_variables.get_primitives_magnetic_field().z() << ", 1/inverse volume:" << 1./inverse_volume << std::endl;
  //  std::cout << "Warning: dual energy not triggered, emag:" << magnetic_energy_density * 1/inverse_volume << ", ekin:" << kinetic_energy << ", eint:" << internal_energy << ", total energy:" << magnetohydro_variables.get_conserved_total_energy() << std::endl;
    if (internal_energy > magnetohydro_variables.get_conserved_total_energy()) {
      std::cout << "Warning: internal energy is GREATER than total energy in update_energy_variables! Eint: " << internal_energy << ", Etot: " << new_energy << ", Emag: " << magnetic_energy_density/inverse_volume << ", Ekin: " << kinetic_energy << std::endl;
    }
  }
    double new_pressure =
        _gamma_minus_one * inverse_volume * internal_energy;
    double new_temperature = _T_conversion_factor * mean_molecular_mass *
                             new_pressure * inverse_density;

    if (std::isnan(new_temperature)) {
        std::cout << "NaN detected! rho: " << density << " P: " << new_pressure << " Etot: " << new_energy << std::endl;
       // exit(1); 
    }

    cmac_assert_message(new_pressure == new_pressure,
                        "old_energy: %g, delta_energy: %g, new_energy: %g",
                        old_energy, delta_energy, new_energy);
    cmac_assert_message(new_energy == new_energy,
                        "old_energy: %g, delta_energy: %g, new_energy: %g",
                        old_energy, delta_energy, new_energy);
    cmac_assert_message(new_temperature == new_temperature,
                        "old_energy: %g, delta_energy: %g, new_energy: %g",
                        old_energy, delta_energy, new_energy);

    cmac_assert_message(new_temperature > 0.0,
                        "old_temp: %g, new_temp: %g, de=%g",
                        ionization_variables.get_temperature(), new_temperature,delta_energy);

#ifdef SAFE_HYDRO_VARIABLES
    new_pressure = std::max(new_pressure, 0.);
    new_energy = std::max(new_energy, 0.);
    new_temperature = std::max(new_temperature, 0.);
#else
    cmac_assert(new_pressure > 0.);
    cmac_assert(new_energy > 0.);
    cmac_assert(new_temperature > 0.);
#endif

 //   magnetohydro_variables.set_conserved_total_energy(new_energy);
    magnetohydro_variables.set_primitives_pressure(new_pressure);
    ionization_variables.set_temperature(new_temperature);

    if (magnetohydro_variables.get_conserved_internal_energy() == 0) {
      std::cout<< "Warning: internal energy is 0 within update_energy_variables!" << std::endl;
    }
  }

  /**
   * @brief Add the energy due to ionization to the given hydrodynamic
   * variables.
   *
   * @param ionization_variables IonizationVariables.
   * @param magnetohydro_variables MagnetoHydroVariables.
   * @param inverse_volume Inverse volume of the cell (in m^-3).
   * @param timestep Integration timestep (in s).
   */
  inline void add_ionization_energy(IonizationVariables &ionization_variables,
                                    MagnetoHydroVariables &magnetohydro_variables,
                                    const double inverse_volume,
                                    const double timestep) const {

    std::cout<<"Warning: add_ionization_energy has been called!" << std::endl;

    if (_do_explicit_heating) {
      double dE = ionization_variables.get_heating(HEATINGTERM_H) *
                        timestep * ionization_variables.get_number_density() /
                        inverse_volume *
                        ionization_variables.get_ionic_fraction(ION_H_n);
#ifdef HAS_HELIUM
   //calculate heating due to helium photoionization
    const double AHe = _abundances.get_abundance(ELEMENT_He);
    const double n = ionization_variables.get_number_density();
    const double h0 = ionization_variables.get_ionic_fraction(ION_H_n);
    const double he0 = ionization_variables.get_ionic_fraction(ION_He_n);
    const double hep = ionization_variables.get_ionic_fraction(ION_He_p1);
     dE += AHe*ionization_variables.get_heating(HEATINGTERM_He) *
                        timestep * ionization_variables.get_number_density() /
                        inverse_volume *
                        ionization_variables.get_ionic_fraction(ION_He_n);
// calculate heating due to hydrogen OTS absorption of HeLyAlpha photon
    const double T4 = 1.e-4*ionization_variables.get_temperature();
    const double sqrtT = std::sqrt(ionization_variables.get_temperature());
    const double alpha_e_2sP = 4.17e-20 * std::pow(T4, -0.861);
    const double pHots = 1. / (1. + 77. * he0 / (sqrtT * h0));
    const double ne = n * (1. - h0 + AHe * hep + 2*AHe*(1. - hep - he0));
    const double nenhep = ne * hep * n * AHe;
    dE += pHots * 1.21765423e-18 * alpha_e_2sP * nenhep/inverse_volume*timestep;
#endif
      if (std::abs(dE) > magnetohydro_variables.get_conserved_internal_energy()){
        std::cout << "Warning: add_ionization_energy, Eint < dE: " << magnetohydro_variables.get_conserved_internal_energy() << ", dE: " << dE << std::endl;
      }
      update_energy_variables(ionization_variables, magnetohydro_variables,
                              inverse_volume, dE);
      return;
    }

    if (_gamma == 1.) {
      const double xH = ionization_variables.get_ionic_fraction(ION_H_n);
      const double mean_molecular_mass = 0.5 * (1. + xH);
      double temperature = 
          _ionised_temperature * (1. - xH) + _neutral_temperature * xH;
      double pressure = _pressure_conversion_factor *
                        magnetohydro_variables.get_primitives_density() * temperature /
                        mean_molecular_mass;

      cmac_assert(temperature == temperature);
      cmac_assert(pressure == pressure);

#ifdef SAFE_HYDRO_VARIABLES
      temperature = std::max(temperature, 0.);
      pressure = std::max(pressure, 0.);
#else
      cmac_assert(temperature >= 0.);
      cmac_assert(pressure >= 0.);
#endif

      ionization_variables.set_temperature(temperature);
      magnetohydro_variables.set_primitives_pressure(pressure);
    } else {
      const double rho = magnetohydro_variables.get_primitives_density();
      if (rho <= 0.) {
        return;
      }
      const double rho_inv = 1. / rho;
      if (std::isinf(rho_inv)) {
        return;
      }
      const double P = magnetohydro_variables.get_primitives_pressure();
      const double xH = ionization_variables.get_ionic_fraction(ION_H_n);
      const double m = magnetohydro_variables.get_conserved_mass();
      const double Tgas_new =
          _ionised_temperature * (1. - xH) + _neutral_temperature * xH;
      const double ufac = _u_conversion_factor / (1. + xH);
      const double ugas_new = ufac * Tgas_new;
      const double ugas_old = _one_over_gamma_minus_one * P * rho_inv;
      const double du = ugas_new - ugas_old;
      const double dE = m * du;

      const CoordinateVector<> velocity = magnetohydro_variables.get_primitives_velocity();
      double cs_squared = _gamma * P / rho; // mgb edit 06.05.2026 - calculate sound speed squared from test pressure and density
      double va_squared = (CoordinateVector<>::dot_product(magnetohydro_variables.get_primitives_magnetic_field(), magnetohydro_variables.get_primitives_magnetic_field()))/ (_mu0 * rho); // mgb edit 06.05.2026 - calculate Alfvén speed squared from magnetic field strength and density
      double v_fast = std::sqrt(cs_squared + va_squared);
      double v_mag = velocity.norm();
      double test_pressure = _gamma_minus_one * rho * ugas_new; // mgb edit 06.05.2026 - calculate new pressure from updated internal energy
     
      double mach_fast = v_mag/(v_fast + 1.e-15); // mgb edit 06.05.2026 - calculate fast magnetosonic Mach number (add small number to avoid division by zero

      if (dE > 0.) {
        magnetohydro_variables.conserved(4) += dE;
        magnetohydro_variables.conserved(9) += dE;

        double new_energy = magnetohydro_variables.get_conserved_total_energy();
        double kinetic_energy =
            0.5 * CoordinateVector<>::dot_product(
                      magnetohydro_variables.get_primitives_velocity(),
                      magnetohydro_variables.get_conserved_momentum());
        double magnetic_energy_density = (CoordinateVector<>::dot_product(magnetohydro_variables.get_primitives_magnetic_field(), magnetohydro_variables.get_primitives_magnetic_field()))/(2.0*_mu0); // mgb edit 21.04.2026 - magnetic energy density =
        double new_eint = magnetohydro_variables.get_conserved_internal_energy();

        double internal_energy;

        if (mach_fast > _mach_limit || (kinetic_energy + magnetic_energy_density * 1/inverse_volume >  0.95 * new_energy) || (test_pressure-magnetohydro_variables.get_primitives_pressure())/(std::min(test_pressure, magnetohydro_variables.get_primitives_pressure())+1.e-10) > 0.25 ) {
        //  std::cout << "Warning: dual energy triggered, emag:" << magnetic_energy_density * 1/inverse_volume << ", ekin:" << kinetic_energy << ", total energy:" << magnetohydro_variables.get_conserved_total_energy() << std::endl;
          internal_energy = new_eint;
          magnetohydro_variables.set_conserved_internal_energy(internal_energy);
          if (internal_energy > new_energy) {
            std::cout << "Warning: internal energy is GREATER than total energy in add_ionization_energy dual energy! Eint: " << internal_energy << ", Etot: " << new_energy << ", Emag: " << magnetic_energy_density/inverse_volume << ", Ekin: " << kinetic_energy << std::endl;
          }
          if (internal_energy < 0) {
            std::cout << "Warning: add_ionization_energy dual energy, Eint < 0: " << internal_energy << std::endl;
          }
          
         // std::cout << "Warning: dual energy, Bx = " << magnetohydro_variables.get_primitives_magnetic_field().x() << ", By = " << magnetohydro_variables.get_primitives_magnetic_field().y() << ", Bz = " << magnetohydro_variables.get_primitives_magnetic_field().z() << std::endl;
        //  std::cout << "Warning: dual energy triggered, emag:" << magnetic_energy_density * 1/inverse_volume << ", ekin:" << kinetic_energy << ", eint:" << internal_energy << ", total energy:" << magnetohydro_variables.get_conserved_total_energy() << std::endl;
        } else {
          internal_energy = new_energy - kinetic_energy - magnetic_energy_density/inverse_volume;
          magnetohydro_variables.set_conserved_internal_energy(internal_energy);
          if (internal_energy > new_energy) {
            std::cout << "Warning: internal energy is GREATER than total energy in add_ionization_energy! Eint: " << internal_energy << ", Etot: " << new_energy << ", Emag: " << magnetic_energy_density/inverse_volume << ", Ekin: " << kinetic_energy << std::endl;
          }
          if (internal_energy < 0) {
            std::cout << "Warning: add_ionization_energy, Eint < 0: " << internal_energy << std::endl;
          }
        //  std::cout << "Warning: dual energy not triggered, Bx = " << magnetohydro_variables.get_primitives_magnetic_field().x() << ", By = " << magnetohydro_variables.get_primitives_magnetic_field().y() << ", Bz = " << magnetohydro_variables.get_primitives_magnetic_field().z() << std::endl;
         // std::cout << "Warning: dual energy not triggered, emag:" << magnetic_energy_density * 1/inverse_volume << ", ekin:" << kinetic_energy << ", eint:" << internal_energy << ", total energy:" << magnetohydro_variables.get_conserved_total_energy() << std::endl;
        
        }

        const double h0 = ionization_variables.get_ionic_fraction(ION_H_n);
#ifdef HAS_HELIUM
        const double hep = ionization_variables.get_ionic_fraction(ION_He_p1);
        const double hepp = 1 - hep - ionization_variables.get_ionic_fraction(ION_He_n);
        const double AHe = _abundances.get_abundance(ELEMENT_He);
        const double mean_molecular_mass = (1.0 + 4.0*AHe)/(2 + AHe - h0 + AHe*hep + 2.0*AHe*hepp);
#else
        const double mean_molecular_mass = 1.0/(2.0-h0);
#endif
        const double pressure =
            _gamma_minus_one * inverse_volume * internal_energy; // mgb edit 21.04.2026 - calculate pressure from total energy by subtracting kinetic and magnetic energy contributions
        double new_temperature = _T_conversion_factor * mean_molecular_mass *
                             pressure / rho;
            magnetohydro_variables.set_primitives_pressure(pressure);
        ionization_variables.set_temperature(new_temperature);
      }
    }

    // apply velocity limiter
    if (magnetohydro_variables.get_primitives_density() > 0.) {
      const double rho = magnetohydro_variables.get_primitives_density();
      const double cs = get_soundspeed(magnetohydro_variables, ionization_variables);
      double Bx = magnetohydro_variables.get_primitives_magnetic_field().x();
      double By = magnetohydro_variables.get_primitives_magnetic_field().y();
      double Bz = magnetohydro_variables.get_primitives_magnetic_field().z();
      double B_squared = Bx*Bx + By*By + Bz*Bz;

      double vA_squared = B_squared/(_mu0 * rho); // mgb edit 21.04.2026 - Alfvén speed squared + sound speed squared
      double vMHD = std::sqrt(cs * cs + vA_squared);

      if (vMHD > _max_velocity) {
        const double factor = _max_velocity / vMHD;
        magnetohydro_variables.primitives(4) *= factor * factor;
      }
    }
    if (std::isnan(magnetohydro_variables.get_conserved_total_energy())) {
        std::cout << "Warning: Total energy is NaN after ionization energy addition for pressure:" << magnetohydro_variables.get_primitives_pressure() << ", density:" << magnetohydro_variables.get_primitives_density() << std::endl;
    }
    if (magnetohydro_variables.get_conserved_total_energy() < magnetohydro_variables.get_conserved_internal_energy()){
        std::cout<<"Warning! Etot < Eint before update in add ionization energy! Etot: " << magnetohydro_variables.get_conserved_total_energy() << ", Eint: " << magnetohydro_variables.get_conserved_internal_energy() << std::endl;
      }
    set_conserved_variables(magnetohydro_variables, 1./inverse_volume);
    if (magnetohydro_variables.get_conserved_total_energy() < magnetohydro_variables.get_conserved_internal_energy()){
        std::cout<<"Warning! Etot < Eint before update in add ionization energy! Etot: " << magnetohydro_variables.get_conserved_total_energy() << ", Eint: " << magnetohydro_variables.get_conserved_internal_energy() << std::endl;
      }

    if (std::isnan(magnetohydro_variables.get_conserved_total_energy())) {
        std::cout << "Warning: Total energy is NaN after ionization energy addition and conserved variable update for pressure:" << magnetohydro_variables.get_primitives_pressure() << ", density:" << magnetohydro_variables.get_primitives_density() << std::endl;
    }
  }

  /**
   * @brief Set the ionization variables based on the given hydrodynamic
   * variables.
   *
   * @param magnetohydro_variables MagnetoHydroVariables.
   * @param ionization_variables IonizationVariables.
   */
  inline void
  hydro_to_ionization(const MagnetoHydroVariables &magnetohydro_variables,
                      IonizationVariables &ionization_variables) const {

    double number_density =
        _n_conversion_factor * magnetohydro_variables.get_primitives_density();

    cmac_assert(number_density == number_density);

#ifdef SAFE_HYDRO_VARIABLES
    number_density = std::max(number_density,0.);
#else
    cmac_assert(number_density >= 0.);
#endif

    ionization_variables.set_number_density(number_density);
  }


  inline void align_temp_to_p(const MagnetoHydroVariables &magnetohydro_variables, // mgb comment 19.05.2026: not used
                                IonizationVariables &ionization_variables) const {


          double inverse_density = 1./magnetohydro_variables.get_primitives_density();

          double pressure = magnetohydro_variables.get_primitives_pressure();
          //double mean_molecular_mass = 0.5*(1 +  ionization_variables.get_ionic_fraction(ION_H_n));
          const double h0 = ionization_variables.get_ionic_fraction(ION_H_n);
#ifdef HAS_HELIUM
          const double hep = ionization_variables.get_ionic_fraction(ION_He_p1);
          const double hepp = 1 - hep - ionization_variables.get_ionic_fraction(ION_He_n);
          const double AHe = _abundances.get_abundance(ELEMENT_He);
          const double mean_molecular_mass = (1.0 + 4.0*AHe)/(2 + AHe - h0 + AHe*hep + 2.0*AHe*hepp);
#else
          const double mean_molecular_mass = 1.0/(2.0-h0);
#endif
          double temperature = _T_conversion_factor*mean_molecular_mass*pressure*inverse_density;
          cmac_assert_message(temperature > 0.0,"t=%g,p=%g,mmm=%g",temperature, pressure,mean_molecular_mass);
          ionization_variables.set_temperature(temperature);


                                }

  inline void align_p_to_temp(MagnetoHydroVariables &magnetohydro_variables, // mgb edit 26.05.2026: potentially helpful function to realign the temperature with pressure and internal energy
                                IonizationVariables &ionization_variables, const double volume) const {


          const double inverse_density = 1./magnetohydro_variables.get_primitives_density();

          const double temperature = ionization_variables.get_temperature();
          //double mean_molecular_mass = 0.5*(1 +  ionization_variables.get_ionic_fraction(ION_H_n));
          const double h0 = ionization_variables.get_ionic_fraction(ION_H_n);
#ifdef HAS_HELIUM
          const double hep = ionization_variables.get_ionic_fraction(ION_He_p1);
          const double hepp = 1 - hep - ionization_variables.get_ionic_fraction(ION_He_n);
          const double AHe = _abundances.get_abundance(ELEMENT_He);
          const double mean_molecular_mass = (1.0 + 4.0*AHe)/(2 + AHe - h0 + AHe*hep + 2.0*AHe*hepp);
#else
          const double mean_molecular_mass = 1.0/(2.0-h0);
#endif
       //   double temperature = _T_conversion_factor*mean_molecular_mass*pressure*inverse_density;
          const double pressure = temperature/(_T_conversion_factor * mean_molecular_mass * inverse_density);
          cmac_assert_message(temperature > 0.0,"t=%g,p=%g,mmm=%g",temperature, pressure,mean_molecular_mass);
          

          // update internal energy too!
          const double internal_energy = pressure*volume/(_gamma_minus_one);
          const double kinetic_energy =
          0.5 * CoordinateVector<>::dot_product(
                    magnetohydro_variables.get_primitives_velocity(),
                    magnetohydro_variables.get_conserved_momentum());
          const double magnetic_energy_density = (CoordinateVector<>::dot_product(magnetohydro_variables.get_primitives_magnetic_field(), magnetohydro_variables.get_primitives_magnetic_field()))/(2.0*_mu0); // mgb edit 21.04.2026 - magnetic energy density = B^2 / (2 mu0)
          const double total_energy = kinetic_energy + magnetic_energy_density * volume + internal_energy;
          
          magnetohydro_variables.set_conserved_internal_energy(internal_energy);
          magnetohydro_variables.set_conserved_total_energy(total_energy);
          magnetohydro_variables.set_primitives_internal_energy(internal_energy/volume);
          magnetohydro_variables.set_primitives_pressure(pressure);


                                }

  /**
   * @brief limit the alfven speed if needed
   *
   * @param B_mag Magnitude of the magnetic field.
   * @param rho Density.
   * @param v_limit Maximum velocity.
   * @return Limited Alfven speed.
   */

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
  
  /**
   * @brief Get the hydrodynamical timestep for the given cell.
   *
   * @param magnetohydro_variables Hydro variables.
   * @param ionization_variables IonizationVariables for the cell.
   * @param volume Volume of the cell (in m^3).
   * @return Corresponding timestep.
   */
  inline double get_timestep(const MagnetoHydroVariables &magnetohydro_variables, // mgb note - may need to edit this to use a different v
                             const IonizationVariables &ionization_variables,
                             const double volume) const {

    double cs = get_soundspeed(magnetohydro_variables, ionization_variables);

    double rho = magnetohydro_variables.get_primitives_density();

    
    double cf = 0.0;
    if (rho > 0.0) {
      const double p = magnetohydro_variables.get_primitives_pressure();
      const CoordinateVector<> B = magnetohydro_variables.get_primitives_magnetic_field();

      cf = get_fast_magnetosonic_speed(rho, p, B, _max_velocity);
    } else {
      cf = cs; // Fallback for vacuum cells
    } // mgb edit 21.04.2026 - magnetohydrodynamic wave speed = sqrt(sound speed squared + Alfvén speed squared)
  
    const double v = magnetohydro_variables.get_primitives_velocity().norm();
    const double R = std::cbrt(0.75 * volume * M_1_PI);

    const double dt = R / (cf + v);
    if (cf > _max_velocity*_max_velocity ){ // || (dt < 1.e10
      std::cout << "Timestep calculation: R=" << R << ", cs=" << cs << ", cf=" << cf << ", v=" << v << ", dt=" << dt << std::endl;  
    }
    cmac_assert(dt > 0.);
    return dt;
  }
};

#endif // MAGNETOHYDRO_HPP
