/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file GradientCalculator.hpp
 *
 * @brief Visitor that calculates the primitive variable gradients for each cell
 * in the grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef GRADIENTCALCULATOR_HPP
#define GRADIENTCALCULATOR_HPP

#include "DensityGridMHD.hpp"
#include "MagnetoHydroBoundaryConditions.hpp"

#include <cfloat>

/**
 * @brief Visitor that calculates the primitive variable gradients for each cell
 * in the grid.
 */
class GradientCalculatorMHD {
public:
  /**
   * @brief Compute the gradient for a single cell of the grid.
   *
   * @param cell DensityGrid::iterator pointing to a single cell of the grid.
   * @param grid_end DensityGrid::iterator pointing to the end of the grid.
   * @param boundaries Boundary condition flags.
   * @param inverse_length_unit_in_SI Conversion factor from SI lengths to
   * internal length unit (in m^-1).
   * @param inverse_surface_area_unit_in_SI Conversion factor from SI surface
   * area to internal surface area unit (in m^-2).
   * @param volume_unit_in_SI Conversion factor from internal volume units to
   * SI units (in m^3).
   */
  static inline void
  compute_gradient(DensityGridMHD::iterator &cell,
                   const DensityGridMHD::iterator &grid_end,
                   const MagnetoHydroBoundaryConditionType *boundaries,
                   const double inverse_length_unit_in_SI = 1.,
                   const double inverse_surface_area_unit_in_SI = 1.,
                   const double volume_unit_in_SI = 1.) {

    // get the cell variables
    const double WL[10] = {cell.get_hydro_variables().primitives(0),
                          cell.get_hydro_variables().primitives(1),
                          cell.get_hydro_variables().primitives(2),
                          cell.get_hydro_variables().primitives(3),
                          cell.get_hydro_variables().primitives(4),
                          cell.get_hydro_variables().primitives(5),
                          cell.get_hydro_variables().primitives(6),
                          cell.get_hydro_variables().primitives(7),
                          cell.get_hydro_variables().primitives(8),
                          cell.get_hydro_variables().primitives(9)}; // mgb edit 21.04.2026 - internal energy scalar field
    const CoordinateVector<> position_L = cell.get_cell_midpoint();

    // first loop over the neighbours: compute gradients
    auto ngbs = cell.get_neighbours();
    double phi_ngb_max[10] = {-DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX};
    double phi_ngb_min[10] = {DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX};
    for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
      // get the neighbour variables
      DensityGridMHD::iterator ngb = std::get< 0 >(*ngbit);
      const CoordinateVector<> midpoint = std::get< 1 >(*ngbit);
      const CoordinateVector<> normal = std::get< 2 >(*ngbit);
      const double surface_area =
          std::get< 3 >(*ngbit) * inverse_surface_area_unit_in_SI;
      const CoordinateVector<> position_R = position_L + std::get< 4 >(*ngbit);

      double WR[10];
      if (ngb != grid_end) {
        WR[0] = ngb.get_hydro_variables().primitives(0);
        WR[1] = ngb.get_hydro_variables().primitives(1);
        WR[2] = ngb.get_hydro_variables().primitives(2);
        WR[3] = ngb.get_hydro_variables().primitives(3);
        WR[4] = ngb.get_hydro_variables().primitives(4);
        WR[5] = ngb.get_hydro_variables().primitives(5);
        WR[6] = ngb.get_hydro_variables().primitives(6);
        WR[7] = ngb.get_hydro_variables().primitives(7);
        WR[8] = ngb.get_hydro_variables().primitives(8);
        WR[9] = ngb.get_hydro_variables().primitives(9); // mgb edit 21.04.2026 - internal energy scalar field
      } else {
        // apply boundary conditions
        WR[0] = WL[0];
        WR[1] = WL[1];
        if (normal[0] < 0. && boundaries[0] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[1] = -WR[1];
        }
        if (normal[0] > 0. && boundaries[1] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[1] = -WR[1];
        }
        WR[2] = WL[2];
        if (normal[1] < 0. && boundaries[2] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[2] = -WR[2];
        }
        if (normal[1] > 0. && boundaries[3] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[2] = -WR[2];
        }
        WR[3] = WL[3];
        if (normal[2] < 0. && boundaries[4] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[3] = -WR[3];
        }
        if (normal[2] > 0. && boundaries[5] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[3] = -WR[3];
        }
        WR[4] = WL[4];
        WR[5] = WL[5];
        if (normal[0] < 0. && boundaries[0] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[5] = -WR[5];
        }
        if (normal[0] > 0. && boundaries[1] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[5] = -WR[5];
        }
        WR[6] = WL[6];
        if (normal[1] < 0. && boundaries[2] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[6] = -WR[6];
        }
        if (normal[1] > 0. && boundaries[3] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[6] = -WR[6];
        }
        WR[7] = WL[7];
        if (normal[2] < 0. && boundaries[4] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[7] = -WR[7];
        }
        if (normal[2] > 0. && boundaries[5] == HYDRO_BOUNDARY_REFLECTIVE) {
          WR[7] = -WR[7];
        }
        WR[8] = WL[8];
        if ((boundaries[0] == HYDRO_BOUNDARY_REFLECTIVE || boundaries[1] == HYDRO_BOUNDARY_REFLECTIVE ||
             boundaries[2] == HYDRO_BOUNDARY_REFLECTIVE || boundaries[3] == HYDRO_BOUNDARY_REFLECTIVE ||
             boundaries[4] == HYDRO_BOUNDARY_REFLECTIVE || boundaries[5] == HYDRO_BOUNDARY_REFLECTIVE))  {
          WR[8] = -WR[8];
        }
        WR[9] = WL[9]; // mgb edit 21.04.2026 - internal energy scalar field
      }

      const CoordinateVector<> halfpoint = 0.5 * (position_L + position_R);
      const CoordinateVector<> rLR =
          (position_L - position_R) * inverse_length_unit_in_SI;
      const CoordinateVector<> cLR =
          (midpoint - halfpoint) * inverse_length_unit_in_SI;
      const double rLR_inv = 1. / rLR.norm();
      const double fac = surface_area * rLR_inv;

      for (uint_fast8_t i = 0; i < 10; ++i) {
        phi_ngb_max[i] = std::max(phi_ngb_max[i], WR[i]);
        phi_ngb_min[i] = std::min(phi_ngb_min[i], WR[i]);
        for (uint_fast8_t j = 0; j < 3; ++j) {
          cell.get_hydro_variables().primitive_gradients(i)[j] +=
              fac * ((WR[i] - WL[i]) * cLR[j] - 0.5 * (WR[i] + WL[i]) * rLR[j]);
        }
      }
    }

    // normalize the gradients
    const double Vinv = volume_unit_in_SI / cell.get_volume();
    for (uint_fast8_t i = 0; i < 10; ++i) {
      for (uint_fast8_t j = 0; j < 3; ++j) {
        cell.get_hydro_variables().primitive_gradients(i)[j] *= Vinv;
      }
    }

    double phi_ext_max[10] = {-DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX};
    double phi_ext_min[10] = {DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX};
    // second loop over the neighbours: slope limit
    for (auto ngbit = ngbs.begin(); ngbit != ngbs.end(); ++ngbit) {
      // get the neighbour variables
      const CoordinateVector<> midpoint = std::get< 1 >(*ngbit);

      const CoordinateVector<> deltaLR =
          (midpoint - position_L) * inverse_length_unit_in_SI;
      for (uint_fast8_t i = 0; i < 10; ++i) {
        const double phi_ext = CoordinateVector<>::dot_product(
            cell.get_hydro_variables().primitive_gradients(i), deltaLR);
        phi_ext_max[i] = std::max(phi_ext_max[i], phi_ext);
        phi_ext_min[i] = std::min(phi_ext_min[i], phi_ext);
      }
    }

    // slope limiting
    for (uint_fast8_t i = 0; i < 10; ++i) {
      double alpha = 0.;
      if (phi_ext_max[i] != 0. && phi_ext_min[i] != 0.) {
        alpha = std::min(
            1., 0.5 * std::min((phi_ngb_max[i] - WL[i]) / phi_ext_max[i],
                               (phi_ngb_min[i] - WL[i]) / phi_ext_min[i]));
      }
      cell.get_hydro_variables().primitive_gradients(i) *= alpha;
    }
  }

  /**
   * @brief Functor that does the gradient computation for a single cell.
   */
  class GradientComputation {
  private:
    /*! @brief Boundary condition flags. */
    const MagnetoHydroBoundaryConditionType *_boundaries;

    /*! @brief Iterator to the end of the DensityGrid. */
    const DensityGridMHD::iterator &_grid_end;

    /*! @brief Conversion factor from SI lengths to internal length unit
     *  (in m^-1).*/
    const double _inverse_length_unit_in_SI;

    /*! @brief Conversion factor from SI surface area to internal surface area
     *  unit (in m^-2). */
    const double _inverse_surface_area_unit_in_SI;

    /*! @brief Conversion factor from internal volume units to SI units
     *  (in m^3).*/
    const double _volume_unit_in_SI;

  public:
    /**
     * @brief Constructor.
     *
     * @param boundaries Boundary condition flags.
     * @param grid_end Iterator to the end of the DensityGrid.
     * @param inverse_length_unit_in_SI Conversion factor from SI lengths to
     * internal length unit (in m^-1).
     * @param inverse_surface_area_unit_in_SI Conversion factor from SI surface
     * area to internal surface area unit (in m^-2).
     * @param volume_unit_in_SI Conversion factor from internal volume units to
     * SI units (in m^3).
     */
    inline GradientComputation(
        const MagnetoHydroBoundaryConditionType *boundaries,
        const DensityGridMHD::iterator &grid_end,
        const double inverse_length_unit_in_SI = 1.,
        const double inverse_surface_area_unit_in_SI = 1.,
        const double volume_unit_in_SI = 1.)
        : _boundaries(boundaries), _grid_end(grid_end),
          _inverse_length_unit_in_SI(inverse_length_unit_in_SI),
          _inverse_surface_area_unit_in_SI(inverse_surface_area_unit_in_SI),
          _volume_unit_in_SI(volume_unit_in_SI) {}

    /**
     * @brief Perform the gradient computation for a single cell of the grid.
     *
     * @param cell DensityGridMHD::iterator pointing to a grid cell.
     */
    inline void operator()(DensityGridMHD::iterator &cell) {
      compute_gradient(cell, _grid_end, _boundaries, _inverse_length_unit_in_SI,
                       _inverse_surface_area_unit_in_SI, _volume_unit_in_SI);
    }
  };
};

#endif // GRADIENTCALCULATORMHD_HPP

