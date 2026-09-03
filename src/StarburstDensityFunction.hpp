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
 * @file StarburstDensityFunction.hpp
 *
 * @brief Starburst CGOLS model density function.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef STARBURSTDENSITYFUNCTION_HPP
#define STARBURSTDENSITYFUNCTION_HPP

#include "CoordinateVector.hpp"
#include "DensityFunction.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"

#include <cmath>

/**
 * @brief Starburst CGOLS model density function.
 *
 *
 */
class StarburstDensityFunction : public DensityFunction {
private:
  
    const double _neutral_fraction;
    const double _temperature;
    const bool _trace_initial_neutral_flag;
    const double _temperature_to_trace;
    const double _R_gas;
    const double _initial_surface_density;
    const double _z_min;
    const double _z_max;
    const double _nz;
    const double _gamma;
    const double _M_stars;
    const double _R_stars;
    const double _z_stars;
    const double _M_halo;
    const double _R_halo;
    const double _concentration;

  /**
   * @brief Get the mean particle mass @f$\mu{} m_p@f$ corresponding to the
   * given neutral fraction.
   *
   * @param neutral_fraction Neutral fraction of hydrogen, @f$x_{\rm{}H}@f$.
   * @return Mean particle mass, @f$\mu{}m_p@f$ (in kg).
   */
  static inline double get_mean_particle_mass(const double neutral_fraction) {
    return 0.5 *
           PhysicalConstants::get_physical_constant(
               PHYSICALCONSTANT_PROTON_MASS) *
           (1. + neutral_fraction);
  }

/**
   * @brief Computes CGOLS V potential
   */
  double get_starburst_potential_value(const double x, const double y, const double z) const {
    const double G = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_NEWTON_CONSTANT);
    
    const double r_cyl_2 = x*x + y*y;
    const double r_sphere  = std::sqrt(r_cyl_2 + z*z);

    // Miyamoto-Nagai Stellar Disk Profile
    const double z_par = std::sqrt(z * z + _z_stars * _z_stars);
    const double stellar_disk_denominator = std::sqrt(r_cyl_2 + std::pow(_R_stars + z_par, 2.));
    const double stellar_disk_numerator =  G * _M_stars;
    const double stellar_disk_potential = -stellar_disk_numerator / stellar_disk_denominator;

    // 2. NFW Dark Matter Halo Potential
    const double rRhalo = r_sphere / _R_halo;
    const double nfw_concentration = std::log(1.0 + _concentration) - (_concentration / (1.0 + _concentration));
    const double nfw_r_factor = std::log(1.0 + rRhalo);

    const double nfw_numerator = G * _M_halo * nfw_r_factor;
    const double nfw_denominator = r_sphere * nfw_concentration;
    const double nfw_potential = -nfw_numerator / nfw_denominator;

    return stellar_disk_potential + nfw_potential;
  }


  double get_normalisation(const double x, const double y, const double z, const double c_s2, const double surface_density, const double zmin, const double zmax, const double nz) const {

    double integral_sum = 0.0;
    const double phi_0_midplane = get_starburst_potential_value(x,y,0.0);

    const double dz = (zmax - zmin) / static_cast<double>(nz);


    for (int k = 0; k < nz; ++k){
        double z_0 = zmin + k * dz;
        double z_1 = zmin + (k+1) * dz;

        double phi_0 = get_starburst_potential_value(x,y,z_0);
        double phi_1 = get_starburst_potential_value(x,y,z_1);

        double f_0 = std::exp((phi_0 - phi_0_midplane) / (c_s2));
        double f_1 = std::exp((phi_1 - phi_0_midplane) / (c_s2));

        integral_sum += 0.5 * (f_0 + f_1) * dz;
    }

    const double rho_0 = surface_density / integral_sum;
    return rho_0;
  }


  CoordinateVector<double> get_starburst_acceleration(const double x, const double y, const double z) const {

    const double G = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_NEWTON_CONSTANT);

    double r_cyl_2 = 0.0;
    double r_sphere = 0.0;

    r_cyl_2 = x * x + y * y; // cylindrical radius
    r_sphere = std::sqrt(x * x + y * y + z * z); // spherical radius


    double ax = 0.0; double ay = 0.0; double az = 0.0;

    const double z_par = std::sqrt(z * z + _z_stars * _z_stars);
    const double stellar_disk_denominator = std::sqrt(r_cyl_2 + std::pow(_R_stars + z_par, 2.));
    const double stellar_disk_numerator =  G * _M_stars;

    ax -= stellar_disk_numerator * x / std::pow(stellar_disk_denominator, 3.);
    ay -= stellar_disk_numerator * y / std::pow(stellar_disk_denominator, 3.);

    az -=  stellar_disk_numerator * z * (_R_stars + z_par) / (z_par * std::pow(stellar_disk_denominator, 3.));


    const double rRhalo = r_sphere / _R_halo;
    const double nfw_concentration = std::log(1.0 + _concentration) - (_concentration / (1.0 + _concentration));
    const double nfw_r_factor = std::log(1.0 + rRhalo);

    const double nfw_acceleration_num = G * _M_halo * (nfw_r_factor * _R_halo * (1.0 + rRhalo) - r_sphere);
    const double nfw_acceleration_den = std::pow(r_sphere, 3.) * _R_halo * nfw_concentration * (1.0 + rRhalo);

    const double nfw_acceleration_factor = nfw_acceleration_num / nfw_acceleration_den;

    ax -= x * nfw_acceleration_factor;
    ay -= y * nfw_acceleration_factor;
    az -= z * nfw_acceleration_factor;

    return CoordinateVector<double>(ax, ay, az);
  }

public:
  /**
   * @brief Constructor.
   *
   * @param disc_z Vertical position of the disc (in m).
   * @param surface_density Surface density of the disc, @f$\Sigma{}_M@f$
   * (in kg m^-2).
   * @param scale_height Scale height of the disc, @f$b_M@f$ (in m).
   * @param gas_fraction Fraction of the total mass content of the disc that is
   * in gas, @f$f_g@f$.
   * @param temperature Constant initial temperature, @f$T@f$ (in K).
   * @param neutral_fraction Constant initial neutral fraction for hydrogen,
   * @f$x_{\rm{}H}@f$.
   * @param observational_disc Use the observational vertical density profile?
   */
  inline StarburstDensityFunction(const double neutral_fraction,
                                    const double temperature,
                                    const bool trace_initial_neutral_flag,
                                    const double temperature_to_trace,
                                    const double R_gas,
                                    const double initial_surface_density,
                                    const double z_min,
                                    const double z_max,
                                    const int nz,
                                    const double gamma, 
                                    const double M_stars, 
                                    const double R_stars, 
                                    const double z_stars,
                                    const double M_halo, 
                                    const double R_halo, 
                                    const double concentration
                                  )
      : _neutral_fraction(neutral_fraction),
        _temperature(temperature),
        _trace_initial_neutral_flag(trace_initial_neutral_flag),
        _temperature_to_trace(temperature_to_trace),
        _R_gas(R_gas),
        _initial_surface_density(initial_surface_density),
        _z_min(z_min),
        _z_max(z_max),
        _nz(nz),
        _gamma(gamma),
        _M_stars(M_stars),
        _R_stars(R_stars),
        _z_stars(z_stars),
        _M_halo(M_halo),
        _R_halo(R_halo),
        _concentration(concentration)
      {}

  /**
   * @brief ParameterFile constructor.
   *
   * We accept the following parameters:
   *  - disc z: Vertical position of the disc (default: 0. pc)
   *  - surface density: Surface density of the disc (default: 30. Msol pc^-2)
   *  - scale height: Scale height of the disc (default: 200. pc)
   *  - gas fraction: Fraction of the total mass content of the disc that is in
   *    gas (default: 0.1)
   *  - temperature: Constant initial temperature (default: 1.e4 K)
   *  - neutral fraction: Constant initial neutral fraction for hydrogen
   *    (default: 1.e-6)
   *  - observationaldisc: Use the observational vertical density profile. If
   *    false, use the hydrostatic-equilibrium profile (default: true)
   *
   * @param params ParameterFile to read from.
   */
  inline StarburstDensityFunction(ParameterFile &params)
      : StarburstDensityFunction(
            params.get_value< double >(
                "DensityFunction:neutral fraction", 1.0),
            params.get_physical_value< QUANTITY_TEMPERATURE >("DensityFunction:temperature", "1.e4 K"),
             params.get_value< bool >("DensityFunction:trace initial neutral flag", false),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:temperature to trace", "500. K"),
            params.get_physical_value< QUANTITY_LENGTH >("DensityFunction:gas radius", "1600 pc"),
            params.get_physical_value< QUANTITY_SURFACE_DENSITY >("DensityFunction:initial surface density", "30. Msol pc^2"),
            params.get_physical_value< QUANTITY_LENGTH >("DensityFunction:z min", "-5. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >("DensityFunction:z max", "5. kpc"),
            params.get_value< int >("DensityFunction:nz", 256),
            params.get_value< double >("Hydro:polytropic index", 1.0),
            params.get_physical_value< QUANTITY_MASS >("ExternalPotential:stellar mass", "1.e10 Msol"),
            params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:stellar scale radius", "800. pc"),
            params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:stellar scale height", "150. pc"),
            params.get_physical_value< QUANTITY_MASS >("ExternalPotential:halo mass", "5.e10 Msol"),
            params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:halo scale radius", "5300. pc"),
            params.get_value< double >("ExternalPotential:halo concentration", 10.0)

            ) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~StarburstDensityFunction() {} /// Lewis's edited density function: mgb note 30.10.2025

  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) {

    const double z = cell.get_cell_midpoint()[2];
    const double y = cell.get_cell_midpoint()[1];
    const double x = cell.get_cell_midpoint()[0];

    const double r = std::sqrt(x*x + y*y);

    const double kB = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_BOLTZMANN);
    const double mu = get_mean_particle_mass(_neutral_fraction);

    const double surface_density = _initial_surface_density * std::exp(-r / _R_gas);
    const double sound_speed = std::sqrt(kB * _temperature / mu); 
    const double c_s2 = sound_speed * sound_speed;


    const double rho_normalisation = get_normalisation(x, y, z, c_s2, surface_density, _z_min, _z_max, _nz);

    const double phi_z = get_starburst_potential_value(x, y, z);
    const double phi_0 = get_starburst_potential_value(x, y, 0.0);

    const double mass_dens = rho_normalisation * std::exp((phi_z - phi_0) / c_s2);
    const double number_density = mass_dens / mu;



    DensityValues values;
    values.set_number_density(number_density);
    values.set_temperature(_temperature);
    values.set_ionic_fraction(ION_H_n, _neutral_fraction);

    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;

    if (r > 0.0) {
        const double dr = 1.0e-3 * _R_gas;
        const double r_shift = r + dr;

        const double x_shift = x * (r_shift/r);
        const double y_shift = y * (r_shift/r);

        const double surface_dens_shift = _initial_surface_density * std::exp(-r_shift / _R_gas);
        const double rho_normalisation_shift = get_normalisation(x_shift, y_shift, z, c_s2, surface_dens_shift, _z_min, _z_max, _nz);
        const double phi_z_shift = get_starburst_potential_value(x_shift, y_shift, z);
        const double phi_0_shift = get_starburst_potential_value(x_shift, y_shift, 0.0);

        const double mass_dens_shift = rho_normalisation_shift * std::exp((phi_z_shift - phi_0_shift) / c_s2);
        const double dln_rho_dr = (std::log(mass_dens_shift) - std::log(mass_dens)) / dr;

        CoordinateVector< double > acceleration = get_starburst_acceleration(x,y,z);
        const double a_r = (acceleration.x() * x + acceleration.y() * y) / r;
        const double v_phi_square = r * (-a_r + (c_s2 / _gamma) * dln_rho_dr);

        double vphi = 0.0;
        if (v_phi_square > 0.0) {
            vphi = std::sqrt(v_phi_square);
        }
        vx = -vphi * (y/r);
        vy = vphi * (x/r);
    }

    CoordinateVector< double > velocity(vx,vy,vz);
    values.set_velocity(velocity);

    if (_trace_initial_neutral_flag == true){
        if (values.get_temperature() <= _temperature_to_trace){
            values.set_initial_neutral_scalar_field(1.0);
            values.set_remaining_initial_neutral_scalar_field(1.0);
        } else {
            values.set_initial_neutral_scalar_field(0.0);
            values.set_remaining_initial_neutral_scalar_field(0.0);
        }
        values.set_cooled_neutral_scalar_field(0.0);
        values.set_remaining_cooled_neutral_scalar_field(0.0);
    } else {
        values.set_initial_neutral_scalar_field(0.0);
        values.set_cooled_neutral_scalar_field(0.0);
        values.set_remaining_initial_neutral_scalar_field(0.0);
        values.set_remaining_cooled_neutral_scalar_field(0.0);
    }

    return values;
  }
};

#endif // STARBURSTDENSITYFUNCTION_HPP
