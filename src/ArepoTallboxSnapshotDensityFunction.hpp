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
#ifndef AREPOTALLBOXSNAPSHOTDENSITYFUNCTION_HPP
#define AREPOTALLBOXSNAPSHOTDENSITYFUNCTION_HPP

#include "CoordinateVector.hpp"
#include "DensityFunction.hpp"
#include "PhysicalConstants.hpp"
#include "HDF5Tools.hpp"
#include "Box.hpp"
#include <cinttypes>
#include <random>
#include <vector>
#include <algorithm>

#include <cmath>

class Log;
class ParameterFile;

/**
 * @brief Disc patch density function.
 *
 * Represents a gas density profile that is initialised from an arepo snapshot of tallbox region from a full galaxy simulation.
 */
class ArepoTallboxSnapshotDensityFunction : public DensityFunction {
private:

  /*! @brief filename */
  std::string _filename;

  const bool _trace_initial_neutral_flag;

  const double _temperature_to_trace;

  const double _dust_gas_ratio;

  const double _fraction_silicates;

  const bool _zero_initial_xy_velocity;

  /*! @brief Box containing the grid. */
  const Box<> _box;

   /*! @brief Number of subgrids in each coordinate direction. */
  const CoordinateVector< uint_fast32_t > _number_of_subgrids;

   /*! @brief Number of cells. */
  const CoordinateVector< uint_fast32_t > _ncell;

  /*! @brief Log to write logging info to. */
  Log *_log;

  
  
  /*! @brief Cartesian density grid (if applicable). */
  DensityValues ***_cartesian_grid;

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

public:
  /**
   * @brief Constructor.
   *
   * 
   */
  ArepoTallboxSnapshotDensityFunction(std::string filename, 
                                            const bool trace_initial_neutral_flag,
                                            const double temperature_to_trace,
                                            const double dust_gas_ratio,
                                            const double fraction_silicates,
                                            const bool zero_initial_xy_velocity,
                                            const Box<> box, 
                                            const CoordinateVector< uint_fast32_t > number_of_subgrids,
                                            const CoordinateVector< uint_fast32_t > number_of_cells,
                                            Log *log
                                  );
    
  ArepoTallboxSnapshotDensityFunction(ParameterFile &params, Log *log = nullptr);

  virtual ~ArepoTallboxSnapshotDensityFunction();

  virtual void initialize();

  virtual void free();

  virtual DensityValues operator()(const Cell &cell);

};

#endif // AREPOTALLBOXSNAPSHOTDENSITYFUNCTION_HPP