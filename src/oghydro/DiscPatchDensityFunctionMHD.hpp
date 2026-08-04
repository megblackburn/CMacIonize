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
 * @file DiscPatchDensityFunction.hpp
 *
 * @brief Disc patch density function.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DISCPATCHDENSITYFUNCTIONMHD_HPP
#define DISCPATCHDENSITYFUNCTIONMHD_HPP

#include "CoordinateVector.hpp"
#include "DensityFunctionMHD.hpp"
#include "ParameterFile.hpp"
#include "PhysicalConstants.hpp"

#include <cmath>

/**
 * @brief Disc patch density function.
 *
 * Represents a gas density profile that is initially in hydrostatic equilibrium
 * with a DiscPatchExternalPotential corresponding to a matter density profile
 * of the form
 * \f[
 *   \rho{}_M = \frac{\Sigma{}_M}{2b_M} \left(\cosh\left(\frac{z}{b_M}\right)
 *   \right)^{-2},
 * \f]
 * with @f$z@f$ the third component of the position, @f$\Sigma{}_M@f$ the
 * surface density of matter in the plane @f$z=0@f$ and @f$b_M@f$ a vertical
 * scale height for the density profile.
 *
 * The gas density profile itself has the general form
 * \f[
 *   \rho{}_g = \frac{\Sigma{}_g}{2b_M} \left(\cosh\left(\frac{z}{b_M}\right)
 *   \right)^{-\frac{2b_M}{b_g}},
 * \f]
 * where @f$\Sigma{}_g@f$ and @f$b_g@f$ are the surface density and scale height
 * for the gas density profile. The latter is given by
 * \f[
 *   b_g = \frac{k_B T}{\mu{} m_p \pi{} G \Sigma{}_M},
 * \f]
 * with @f$T@f$ the hydrostatic equilibrium temperature and @f$\mu{}@f$ the
 * mean molecular weight,
 * \f[
 *   \mu{} = \frac{1}{2} (1 + x_{\rm{}H}),
 * \f]
 * with @f$x_{\rm{}H}@f$ the hydrogen neutral fraction (we assume a hydrogen
 * only gas). @f$k_B@f$, @f$m_p@f$ and @f$G@f$ are respectively Bolzmann's
 * constant, the proton mass and Newton's constant.
 *
 * The gas surface density @f$\Sigma{}_g@f$ can be related to the total surface
 * density @f$\Sigma{}_M@f$ by imposing a fixed mass ratio
 * @f$f_g = \frac{M_g}{M_M}@f$, with
 * \f[
 *   M_X = \int_{-\infty{}}^{+\infty{}} \rho{}_X (z)~{\rm{}d}z, X = [M, g].
 * \f]
 * The corresponding expression is
 * \f[
 *   \Sigma{}_g = \frac{2}{I\left(-\frac{2b_M}{b_g}\right)} f_g \Sigma{}_M,
 * \f]
 * with
 * \f[
 *   I(d) = \int_{-\infty{}}^{+\infty{}} (\cosh(x))^d~{\rm{}d}x.
 * \f]
 * This integral has to be evaluated numerically. We use a third order
 * polynomial fit in log-log space to approximate it.
 */
class DiscPatchDensityFunctionMHD : public DensityFunctionMHD {
private:
  /*! @brief Vertical position of the disc (in m). */
  const double _disc_z;

  /*! @brief Inverse vertical scale height of the disc, @f$\frac{1}{b_M}@f$
   *  (in m^-1). */
  const double _b_inv;

  /*! @brief Exponent of the density profile, @f$-\frac{2b_M}{b_g}@f$. */
  const double _exponent;

  /*! @brief Norm of the density profile, @f$\frac{\Sigma{}_g}{2b_Mm_p}@f$
   *  (in m^-3). */
  const double _density_norm;

  /*! @brief Constant initial temperature, @f$T@f$ (in K). */
  const double _temperature;

  /*! @brief Constant initial B_x, @f$T@f$ (in K). */
  const double _init_Bx;

  /*! @brief Constant initial B_y, @f$T@f$ (in K). */
  const double _init_By;

  /*! @brief Constant initial B_z, @f$T@f$ (in K). */
  const double _init_Bz;

  /*! @brief Constant initial B_scalar, @f$T@f$ (in K). */
  const double _init_B_scalar;

  /*! @brief Constant initial plasma_beta, @f$T@f$ (in K). */
  const double _plasma_beta;


  /*! @brief Constant initial neutral fraction for hydrogen,
   *  @f$x_{\rm{}H}@f$. */
  const double _neutral_fraction;

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
   * @brief Get the scale height of a gas disc in hydrostatic equilibrium with
   * a potential with the given surface density.
   *
   * @param surface_density Surface density of the disc, @f$\Sigma{}_M@f$
   * (in kg m^-2).
   * @param temperature Hydrostatic equilibrium temperature of the gas, @f$T@f$
   * (in K).
   * @param neutral_fraction Neutral fraction for hydrogen, @f$x_{\rm{}H}@f$.
   * @return Scale height for the corresponding equilibrium gas disc, @f$b_g@f$.
   */
  static inline double
  get_gas_disc_scale_height(const double surface_density,
                            const double temperature,
                            const double neutral_fraction) {
    return (PhysicalConstants::get_physical_constant(
                PHYSICALCONSTANT_BOLTZMANN) *
            temperature) /
           (get_mean_particle_mass(neutral_fraction) * M_PI *
            PhysicalConstants::get_physical_constant(
                PHYSICALCONSTANT_NEWTON_CONSTANT) *
            surface_density);
  }

  /**
   * @brief Get the mass fraction factor for the given density profile exponent
   * @f$d@f$.
   *
   * This is the numerical integration of
   * \f[
   *   \int_{-\infty{}}^{+\infty{}} \left(\cosh(x)\right)^d~{\rm{}d}x.
   * \f]
   *
   * We made a fit to this expression for values of the scale height ratio
   * @f$r = \frac{b_M}{b_g} = -\frac{1}{2} d@f$ in log-log space.
   *
   * @param exponent Exponent @f$d@f$.
   * @return Mass fraction factor: extra factor needed to convert from total
   * surface density @f$\Sigma{}_M@f$ into gas surface density @f$\Sigma{}_g@f$
   * due to the different slope of the gas density profile.
   */
  static inline double get_mass_fraction_factor(const double exponent) {
    const double x = std::log10(-0.5 * exponent);
    const double x2 = x * x;
    // polynomial fit to actual integral
    const double y =
        0.01499337 * x2 * x - 0.08454788 * x2 + 0.63503798 * x - 0.01018254;
    return std::pow(10., y);
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
   */
  inline DiscPatchDensityFunctionMHD(const double disc_z,
                                  const double surface_density,
                                  const double scale_height,
                                  const double gas_fraction,
                                  const double temperature,
                                  const double init_Bx,
                                  const double init_By,
                                  const double init_Bz,
                                  const double init_B_scalar,
                                  const double plasma_beta,
                                  const double neutral_fraction)
      : _disc_z(disc_z), _b_inv(1. / scale_height),
        _exponent(-2. * scale_height /
                  get_gas_disc_scale_height(surface_density, temperature,
                                            neutral_fraction)),
        _density_norm(0.5 * gas_fraction * surface_density *
                      get_mass_fraction_factor(_exponent) * _b_inv /
                      PhysicalConstants::get_physical_constant(
                          PHYSICALCONSTANT_PROTON_MASS)),
        _temperature(temperature), _init_Bx(init_Bx), _init_By(init_By), _init_Bz(init_Bz), _init_B_scalar(init_B_scalar), _plasma_beta(plasma_beta), _neutral_fraction(neutral_fraction) {}

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
   *
   * @param params ParameterFile to read from.
   */
  inline DiscPatchDensityFunctionMHD(ParameterFile &params)
      : DiscPatchDensityFunctionMHD(
            params.get_physical_value< QUANTITY_LENGTH >(
                "DensityFunction:disc z", "0. m"),
            params.get_physical_value< QUANTITY_SURFACE_DENSITY >(
                "DensityFunction:surface density", "30. Msol pc^-2"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "DensityFunction:scale height", "200. pc"),
            params.get_value< double >("DensityFunction:gas fraction", 0.1),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:temperature", "1.e4 K"),
            params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                "DensityFunction:init Bx", "0. T"),
            params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                "DensityFunction:init By", "0. T"),
            params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                "DensityFunction:init Bz", "0. T"),
            params.get_physical_value< QUANTITY_MAGNETIC_FIELD >(
                "DensityFunction:init B_scalar", "0. T"),
            params.get_value< double > ("DensityFunction:plasma beta", 1),
            params.get_value< double >("DensityFunction:neutral fraction",
                                       1e-6)) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~DiscPatchDensityFunctionMHD() {} /// Lewis's edited density function: mgb note 30.10.2025

  /**
   * @brief Function that gives the density for a given cell.
   *
   * @param cell Geometrical information about the cell.
   * @return Initial physical field values for that cell.
   */
  virtual DensityValues operator()(const Cell &cell) {

   // HEY FUTURE ME! If you're gonna swap this back remember this equation is in cm^-3 not m^-3.
    const double abs_z = std::abs(cell.get_cell_midpoint().z() - _disc_z)/3.086e+19;
    const double nH = 0.47*std::exp(-0.5*std::pow(abs_z/0.09,2)) + 0.13*std::exp(-0.5*std::pow(abs_z/0.225,2))
         + 0.077*std::exp(-1*(abs_z/0.403)) +
           0.025*std::exp(-abs_z);

  //  const double dz = cell.get_cell_midpoint().z() - _disc_z;
  //  const double cosh = std::cosh(dz * _b_inv);
  //  const double nH = _density_norm * std::pow(cosh, _exponent);
    const double dens_profile = 0.47*std::exp(-0.5*std::pow(abs_z/0.09,2)) + 0.13*std::exp(-0.5*std::pow(abs_z/0.225,2))
         + 0.077*std::exp(-1*(abs_z/0.403)) +
           0.025*std::exp(-abs_z);

    const double temperature = _temperature/dens_profile;

    DensityValues values;
    values.set_number_density(nH*1e6);
    values.set_temperature(temperature);
    values.set_ionic_fraction(ION_H_n, _neutral_fraction);
    /*
    const double dens_profile = 0.47*std::exp(-0.5*std::pow(abs_z/0.09,2)) + 0.13*std::exp(-0.5*std::pow(abs_z/0.225,2))
         + 0.077*std::exp(-1*(abs_z/0.403)) +
           0.025*std::exp(-abs_z);*/

    const double _mu0 = PhysicalConstants::get_physical_constant(PHYSICALCONSTANT_MAGNETC_PERMEABILITY);
    
    const double k_B = 1.380649e-23; 
    const double local_number_density_m3 = values.get_number_density(); 
    const double local_gas_pressure = local_number_density_m3 * k_B * _temperature;

    const double target_beta = _plasma_beta;
    
    double beta_scaling = std::sqrt((2.0 * _mu0 * local_gas_pressure) / target_beta);


    double Bx = _init_Bx * beta_scaling; //* dens_profile;
    double By = _init_By * beta_scaling; //* dens_profile;
    double Bz = _init_Bz * beta_scaling;// * dens_profile;

    values.set_magnetic_field(CoordinateVector<>(Bx,By,Bz));
    values.set_B_scalar(_init_B_scalar);

    return values;
  }
};

#endif // DISCPATCHDENSITYFUNCTION_HPP
