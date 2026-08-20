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
    PotentialMode _mode;  
    double _R_patch;      
    
    // --- Miyamoto-Nagai Disk Parameters ---
    double _M_disk;   
    double _a_disk;   
    double _b_disk;   

    // --- Hernquist Bulge Parameters (Set to 0 for CGOLS V) ---
    double _M_bulge;  
    double _c_bulge;  

    // --- NFW Dark Matter Halo Parameters ---
    double _M_halo;   
    double _r_s;      
    double _c_nfw;

public:
    /**
     * @brief Explicit constructor matched exactly to your CGOLS V / M82 parameters.
     */
    inline StarburstExternalPotential(
        const PotentialMode mode = POTENTIAL_GLOBAL,
        const double patch_radius_kpc = 1.5,
        const double M_disk  = 1.0e10 * 1.98847e30, // 10^10 Msun
        const double a_disk  = 800.0  * 3.08568e16, // 0.8 kpc
        const double b_disk  = 150.0  * 3.08568e16, // 0.15 kpc
        const double M_bulge = 0.0    * 1.98847e30, // Removed for CGOLS V
        const double c_bulge = 0.0    * 3.08568e16,
        const double M_halo  = 5.0e10 * 1.98847e30, // 5x10^10 Msun
        const double r_s     = 5300.0 * 3.08568e16, // 5.3 kpc
        const double c_nfw   = 10.0
    ) : _mode(mode), _M_disk(M_disk), _a_disk(a_disk), _b_disk(b_disk),
        _M_bulge(M_bulge), _c_bulge(c_bulge), _M_halo(M_halo), _r_s(r_s), _c_nfw(c_nfw) 
    {
        _G = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_NEWTON_CONSTANT);
        _R_patch = patch_radius_kpc * 3.08568e19; // Corrected to 1e19 for kpc conversion scale
    }

    /**
     * @brief ParameterFile constructor reading directly from runtime script.
     */
    inline StarburstExternalPotential(ParameterFile &params) {
        _G = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_NEWTON_CONSTANT);
        
        std::string mode_str = params.get_value<std::string>("ExternalPotential:mode", "global");
        if (mode_str == "global") {
            _mode = POTENTIAL_GLOBAL;
        } else if (mode_str == "patch") {
            _mode = POTENTIAL_PATCH;
        } else {
            throw std::runtime_error("StarburstExternalPotential: Invalid mode string! Use 'global' or 'patch'.");
        }

        // Dynamically parsed from parameter script - completely eliminating hardcoding
        _R_patch = params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:patch radius", "1.5 kpc");
        _M_disk  = params.get_physical_value< QUANTITY_MASS >("ExternalPotential:disk mass", "1.0e10 Msol");
        _a_disk  = params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:disk scale length", "800. pc");
        _b_disk  = params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:disk scale height", "150. pc");
        _M_bulge = params.get_physical_value< QUANTITY_MASS >("ExternalPotential:bulge mass", "0.0 Msol");
        _c_bulge = params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:bulge scale radius", "0.0 pc");
        _M_halo  = params.get_physical_value< QUANTITY_MASS >("ExternalPotential:halo mass", "5.0e10 Msol");
        _r_s     = params.get_physical_value< QUANTITY_LENGTH >("ExternalPotential:halo scale radius", "5300. pc");
        _c_nfw   = params.get_value< double >("ExternalPotential:concentration", 10.0);
    }

    virtual ~StarburstExternalPotential() {}
    /**
     * @brief Computes combined 3D fluid acceleration vectors.
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
        // B. HERNQUIST CENTRAL BULGE COMPONENT (Automatically skips if mass is 0)
        // ====================================================================
        if (_M_bulge > 0.0) {
            const double bulge_factor = -_G * _M_bulge / (r * (r + _c_bulge) * (r + _c_bulge));
            if (_mode == POTENTIAL_GLOBAL) {
                ax += bulge_factor * x;
                ay += bulge_factor * y;
            }
            az += bulge_factor * z;
        }

        // ====================================================================
        // C. NAVARRO-FRENK-WHITE (NFW) DARK MATTER HALO COMPONENT
        // ====================================================================
        const double x_halo = r / _r_s;
        const double nfw_norm = std::log(1.0 + _c_nfw) - (_c_nfw / (1.0 + _c_nfw));
        const double g_x = std::log(1.0 + x_halo) - (x_halo / (1.0 + x_halo));
        const double halo_factor = -_G * _M_halo * g_x / (r * r * r * nfw_norm);

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
