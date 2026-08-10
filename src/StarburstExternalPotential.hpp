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
 * @brief Starburst External Potential based on CGOLS (Even E. Schnieder and S. Alwin Mao 2024).
 * 
 * Supports both full 3D global galaxy tracking and local stratified tall-box patch configurations
 * via a runtime execution flag parameter.
 */
class StarburstExternalPotential : public ExternalPotential {
public:
    /**
     * @brief Runtime operational flags determining the geometric coordinate mapping system.
     */
    enum PotentialMode {
        POTENTIAL_GLOBAL = 0, // (0,0,0) is galactic center; calculates complete 3D vectors (ax, ay, az).
        POTENTIAL_PATCH  = 1  // (0,0) is box center; locks radial radius to R_patch, returns vertical (0, 0, az).
    };

private:
    double _G;
    PotentialMode _mode;  // Runtime structural toggle switch
    double _R_patch;      // Localised orbital distance from galaxy center (m) [Used only in PATCH mode]
    
    // --- Miyamoto-Nagai Disk Parameters ---
    double _M_disk;   // Total Mass of the Stellar Disk (kg)
    double _a_disk;   // Radial Scale Length of Disk (m)
    double _b_disk;   // Vertical Scale Thickness of Disk (m)

    // --- Hernquist Bulge Parameters ---
    double _M_bulge;  // Total Mass of Central Bulge Core (kg)
    double _c_bulge;  // Scale Radius of Bulge Core (m)

    // --- NFW Dark Matter Halo Parameters ---
    double _M_halo;   // Scale Characteristic Mass of Dark Matter Halo (kg)
    double _r_s;      // Core Scale Radius of Halo Profile (m)

public:
    /**
     * @brief Explicit C++ parameter constructor.
     */
    inline StarburstExternalPotential(
        const PotentialMode mode = POTENTIAL_PATCH,
        const double patch_radius_kpc = 1.5,
        const double M_disk  = 6.0e9  * 1.989e30, 
        const double a_disk  = 1100.  * 3.086e16, 
        const double b_disk  = 150.   * 3.086e16, 
        const double M_bulge = 5.0e8  * 1.989e30, 
        const double c_bulge = 180.   * 3.086e16, 
        const double M_halo  = 5.0e10 * 1.989e30, 
        const double r_s     = 9000.  * 3.086e16  
    ) : _mode(mode), _M_disk(M_disk), _a_disk(a_disk), _b_disk(b_disk),
        _M_bulge(M_bulge), _c_bulge(c_bulge), _M_halo(M_halo), _r_s(r_s) 
    {
        _G = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_NEWTON_CONSTANT);
        _R_patch = patch_radius_kpc * 3.086e16; 
    }

    /**
     * @brief ParameterFile configuration runtime parsing constructor.
     */
    inline StarburstExternalPotential(ParameterFile &params) {
        _G = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_NEWTON_CONSTANT);
        
        std::string mode_str = params.get_value<std::string>("ExternalPotential:mode", "patch");
        if (mode_str == "global") {
            _mode = POTENTIAL_GLOBAL;
        } else if (mode_str == "patch") {
            _mode = POTENTIAL_PATCH;
        } else {
            throw std::runtime_error("StarburstExternalPotential: Invalid mode string! Use 'global' or 'patch'.");
        }

        _R_patch = params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:patch radius", "1.5 kpc");
        _M_disk  = params.get_physical_value< QUANTITY_MASS >("ExternalPotential:disk mass", "6.0e9 Msol");
        _a_disk  = params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:disk scale length", "1100. pc");
        _b_disk  = params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:disk scale height", "150. pc");
        _M_bulge = params.get_physical_value< QUANTITY_MASS >("ExternalPotential:bulge mass", "5.0e8 Msol");
        _c_bulge = params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:bulge scale radius", "180. pc");
        _M_halo  = params.get_physical_value< QUANTITY_MASS >("ExternalPotential:halo mass", "5.0e10 Msol");
        _r_s     = params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:halo scale radius", "9000. pc");
    }

    virtual ~StarburstExternalPotential() {}

    /**
     * @brief Evaluates spatial parameters and calculates combined multi-component gravity.
     */
    virtual CoordinateVector<> get_acceleration(const CoordinateVector<> position) const override {
        const double x = position.x();
        const double y = position.y();
        const double z = position.z();

        double R2 = 0.0;
        double r  = 0.0;

        if (_mode == POTENTIAL_GLOBAL) {
            R2 = x * x + y * y;
            r  = std::sqrt(R2 + z * z + 1e-20);
        } else { // POTENTIAL_PATCH
            R2 = _R_patch * _R_patch;
            r  = std::sqrt(R2 + z * z + 1e-20);
        }

        double ax = 0.0; double ay = 0.0; double az = 0.0;

        // ====================================================================
        // A. MIYAMOTO-NAGAI DISK ACCELERATION COMPONENT
        // ====================================================================
        const double z_param = std::sqrt(z*z + _b_disk*_b_disk);
        const double disk_denom_inner = _a_disk + z_param;
        const double disk_denom = std::sqrt(R2 + disk_denom_inner*disk_denom_inner);
        const double disk_factor = -_G * _M_disk / (disk_denom * disk_denom * disk_denom);

        if (_mode == POTENTIAL_GLOBAL) {
            ax += disk_factor * x;
            ay += disk_factor * y;
        }
        az += disk_factor * z * (disk_denom_inner / z_param);

        // ====================================================================
        // B. HERNQUIST CENTRAL BULGE ACCELERATION COMPONENT
        // ====================================================================
        const double bulge_factor = -_G * _M_bulge / (r * (r + _c_bulge) * (r + _c_bulge));
        
        if (_mode == POTENTIAL_GLOBAL) {
            ax += bulge_factor * x;
            ay += bulge_factor * y;
        }
        az += bulge_factor * z;

        // ====================================================================
        // C. NFW DARK MATTER HALO ACCELERATION COMPONENT
        // ====================================================================
        const double x_halo = r / _r_s;
        const double g_x = std::log(1.0 + x_halo) - (x_halo / (1.0 + x_halo));
        const double halo_factor = -_G * _M_halo * g_x / (r * r * r);

        if (_mode == POTENTIAL_GLOBAL) {
            ax += halo_factor * x;
            ay += halo_factor * y;
        }
        az += halo_factor * z;

        if (_mode == POTENTIAL_PATCH) {
            return CoordinateVector<>(0.0, 0.0, az);
        }
        
        return CoordinateVector<>(ax, ay, az);
    }
};

#endif // STARBURSTEXTERNALPOTENTIAL_HPP
