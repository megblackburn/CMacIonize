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
 * @file StarburstExternalPotential.hpp
 *
 * @brief Starburst external potential.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */

#ifndef STARBURSTEXTERNALPOTENTIAL_HPP
#define STARBURSTEXTERNALPOTENTIAL_HPP

#include <cmath>
#include <string>
#include <stdexcept>

/**
 * @brief Starburst External Potential matched perfectly to CGOLS V profiles.
 */
class StarburstExternalPotential : public ExternalPotential {
public:
    enum PotentialMode {
        POTENTIAL_GLOBAL = 0, // (0,0,0) is galactic center; complete 3D vectors
        POTENTIAL_PATCH  = 1  // (0,0) is box center; vertical gravity (0, 0, az) only
    };

private:
    double _G;
    std::string _mode_str;
    PotentialMode _mode;  
    const double _patch_radius;      
    
    // --- Miyamoto-Nagai Disk Parameters ---
    const double _M_stars;   
    const double _R_stars;   
    const double _z_stars;   

    // --- NFW Dark Matter Halo Parameters ---
    const double _M_halo;  
    const double _R_halo;
    const double _concentration;  


public:
    /**
     * @brief Explicit constructor matched exactly to your CGOLS V / M82 parameters.
     */
    inline StarburstExternalPotential(
        const std::string mode,
        const double patch_radius,
        const double M_stars, 
        const double R_stars, 
        const double z_stars,
        const double M_halo, 
        const double R_halo, 
        const double concentration
    ) : _mode_str(mode), _patch_radius(patch_radius), _M_stars(M_stars), _R_stars(R_stars), _z_stars(z_stars), _M_halo(M_halo), _R_halo(R_halo), _concentration(concentration)
    {
        _G = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_NEWTON_CONSTANT);

        if (_mode_str == "global") {
            _mode = POTENTIAL_GLOBAL;
        } else if (_mode_str == "patch") {
            _mode = POTENTIAL_PATCH;
        } else {
            throw std::runtime_error("StarburstExternalPotential: Invalid mode string! Use 'global' or 'patch'.");
        }   
           
    }

    /**
     * @brief ParameterFile constructor reading directly from runtime script.
     */
    StarburstExternalPotential(ParameterFile &params) : StarburstExternalPotential(
       params.get_value<std::string>("ExternalPotential:mode", "global"),
        params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:patch radius", "1.5 kpc"),
        params.get_physical_value< QUANTITY_MASS >("ExternalPotential:stellar mass", "1.e10 Msol"),
        params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:stellar scale radius", "800. pc"),
        params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:stellar scale height", "150. pc"),
        params.get_physical_value< QUANTITY_MASS >("ExternalPotential:halo mass", "5.e10 Msol"),
        params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:halo scale radius", "5300. pc"),
        params.get_value< double >("ExternalPotential:halo concentration", 10.0)
    ) {  }

    virtual ~StarburstExternalPotential() {}
    /**
     * @brief Computes combined 3D fluid acceleration vectors.
     */
    virtual CoordinateVector<> get_acceleration(const CoordinateVector<> position) const override {
        const double x = position.x();
        const double y = position.y();
        const double z = position.z();

        double r_cyl_2 = 0.0;
        double r_sphere = 0.0;

        if (_mode == POTENTIAL_GLOBAL) {
            r_cyl_2 = x * x + y * y; // cylindrical radius
            r_sphere = std::sqrt(x * x + y * y + z * z); // spherical radius
        } else { // POTENTIAL_PATCH
            r_cyl_2 = _patch_radius * _patch_radius;
            r_sphere = std::sqrt(_patch_radius * _patch_radius + z * z);
        }

        double ax = 0.0; double ay = 0.0; double az = 0.0;

        // ====================================================================
        // A. MIYAMOTO-NAGAI DISK ACCELERATION COMPONENT
        // ====================================================================
        const double z_par = std::sqrt(z * z + _z_stars * _z_stars);
        const double stellar_disk_denominator = std::sqrt(r_cyl_2 + std::pow(_R_stars + z_par, 2.));
        const double stellar_disk_numerator =  _G * _M_stars;

        if (_mode == POTENTIAL_GLOBAL) {
            ax -= stellar_disk_numerator * x / std::pow(stellar_disk_denominator, 3.);
            ay -= stellar_disk_numerator * y / std::pow(stellar_disk_denominator, 3.);
        }
        az -=  stellar_disk_numerator * z * (_R_stars + z_par) / (z_par * std::pow(stellar_disk_denominator, 3.));


        // ====================================================================
        // C. NAVARRO-FRENK-WHITE (NFW) DARK MATTER HALO COMPONENT
        // ====================================================================

        const double rRhalo = r_sphere / _R_halo;
        const double nfw_concentration = std::log(1.0 + _concentration) - (_concentration / (1.0 + _concentration));
        const double nfw_r_factor = std::log(1.0 + rRhalo);

        const double nfw_acceleration_num = _G * _M_halo * (nfw_r_factor * _R_halo * (1.0 + rRhalo) - r_sphere);
        const double nfw_acceleration_den = std::pow(r_sphere, 3.) * _R_halo * nfw_concentration * (1.0 + rRhalo);

        const double nfw_acceleration_factor = nfw_acceleration_num / nfw_acceleration_den;

        if (_mode == POTENTIAL_GLOBAL) {
            ax -= x * nfw_acceleration_factor;
            ay -= y * nfw_acceleration_factor;
        }

        az -= z * nfw_acceleration_factor;

        if (_mode == POTENTIAL_PATCH) {
            return CoordinateVector<>(0.0, 0.0, az);
        }
        
        return CoordinateVector<>(ax, ay, az);
    }
};

#endif // STARBURSTEXTERNALPOTENTIAL_HPP
