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
 * @file SILCCExternalPotential.hpp
 *
 * @brief Disc patch external potential.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#ifndef SILCCEXTERNALPOTENTIAL_HPP
#define SILCCEXTERNALPOTENTIAL_HPP

#include "CoordinateVector.hpp"
#include "ExternalPotential.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"

#include <cmath>

/**
 * @brief Disc patch external potential.
 *
 * The potential is based on Rathjen et al. 2021 and Li et al. 2017 - consists of a stellar component and NFW dark matter component.
 */
class SILCCExternalPotential : public ExternalPotential {
private:
  /*! @brief Params for stellar component */
  const double _disc_z;
  const double _stellar_surface_density;
  const double _gas_surface_density;
  const double _stellar_scale_height;

  /*! @brief Params for NFW dark matter component */
  const double _R_vir;
  const double _rho_dm;
  const double _concentration;
  const double _R_sol;

  
  CoordinateVector<double> get_stellar_acceleration(const CoordinateVector<> position) const {

    const double G = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_NEWTON_CONSTANT);

    const double dz = position.z() - _disc_z;
    const double f_star = _stellar_surface_density / (_stellar_surface_density + _gas_surface_density);

    const double az = - 1/f_star * 2 * M_PI * G * _stellar_surface_density * std::tanh(dz / _stellar_scale_height);
    return CoordinateVector<>(0., 0., az);
  }

  CoordinateVector<double> get_nfw_acceleration(const CoordinateVector<> position) const {

    const double G = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_NEWTON_CONSTANT);

    const double z = position.z();

    const double r  = std::sqrt(_R_sol*_R_sol + z*z);

    const double R_S = _R_vir / _concentration;
    const double rRs = r / R_S;

    const double nfw_numerator = G * 4 * M_PI * _rho_dm * R_S*R_S * (std::log(1.0 + rRs) - (r/(r+R_S)));
    const double nfw_denominator = std::pow(r, 3.0);

    const double nfw_acceleration_factor = nfw_numerator / nfw_denominator;

    // we do not require x and y components in disc patch
   // ax = x * nfw_acceleration_factor;
   // ay = y * nfw_acceleration_factor;
    const double az = - z * nfw_acceleration_factor;

    return CoordinateVector<double>(0., 0., az);
  }



public:
  /**
   * @brief Constructor.
   *
   * @param disc_z Vertical position of the disc (in m).
   * @param surface_density Surface density of the disc (in kg m^-2).
   * @param scale_height Scale height of the disc (in m).
   */
  inline SILCCExternalPotential(const double disc_z,
                                    const double stellar_surface_density,
                                    const double gas_surface_density,
                                    const double stellar_scale_height,
                                    const double R_vir,
                                    const double rho_dm,
                                    const double concentration,
                                    const double R_sol)
      : _disc_z(disc_z), _stellar_surface_density(stellar_surface_density), _gas_surface_density(gas_surface_density), _stellar_scale_height(stellar_scale_height), _R_vir(R_vir), _rho_dm(rho_dm), _concentration(concentration), _R_sol(R_sol) {}

  /**
   * @brief ParameterFile constructor.
   *
   * We accept the following parameters (defaults based on Creasey, Theuns &
   * Bower, 2013):
   *  - disc z: Vertical position of the disc (default: 0. pc)
   *  - surface density: Surface density of the disc (default: 12. Msol pc^-2)
   *  - scale height: Scale height of the disc (default: 200. pc)
   *
   * @param params ParameterFile to read from.
   */
  inline SILCCExternalPotential(ParameterFile &params)
      : SILCCExternalPotential(
            params.get_physical_value< QUANTITY_LENGTH >(
                "ExternalPotential:disc z", "0. m"),
            params.get_physical_value< QUANTITY_SURFACE_DENSITY >(
                "ExternalPotential:stellar surface density", "30.0 Msol pc^-2"),
            params.get_physical_value< QUANTITY_SURFACE_DENSITY >(
                "ExternalPotential:gas surface density", "10.0 Msol pc^-2"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "ExternalPotential:stellar scale height", "300. pc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "ExternalPotential:virial radius", "200. kpc"),
            params.get_physical_value< QUANTITY_DENSITY >(
                "ExternalPotential:dark matter density", "0.013 Msol pc^-3"),
            params.get_value< double >(
                "ExternalPotential:concentration", 12.0), 
            params.get_physical_value< QUANTITY_LENGTH >(
                "ExternalPotential:solar radius", "8.0 kpc")) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~SILCCExternalPotential() {}

  /**
   * @brief Get the acceleration caused by the external disc on a mass at the
   * given position.
   *
   * @param position Position (in m).
   * @return Acceleration (in m s^-2).
   */
  virtual CoordinateVector<>
  get_acceleration(const CoordinateVector<> position) const {
    const CoordinateVector<> stellar_acceleration =
        get_stellar_acceleration(position);
    const CoordinateVector<> nfw_acceleration = get_nfw_acceleration(position);
    const CoordinateVector<> total_acceleration =
        stellar_acceleration + nfw_acceleration;
    return total_acceleration;
  }
};

#endif // SILCCEXTERNALPOTENTIAL_HPP
