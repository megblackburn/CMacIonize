/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file DensityGridWriter.hpp
 *
 * @brief Snapshot file writer for the DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYGRIDWRITERMHD_HPP
#define DENSITYGRIDWRITERMHD_HPP

#include "DensityGridWriterFieldsMHD.hpp"
#include "Log.hpp"

#include <cstdlib>
#include <string>

class DensityGridMHD;
class DensitySubGridMHD;
template < class _subgrid_type_ > class DensitySubGridCreatorMHD;
class MagnetoHydroDensitySubGrid;
class ParameterFile;
class InternalMagnetoHydroUnits;

/**
 * @brief Snapshot file writer for the DensityGrid.
 */
class DensityGridWriterMHD {
protected:
  /*! @brief Name of the folder where output files should be placed. */
  std::string _output_folder;

  /*! @brief Flag specifying whether or not hydro is active. */
  const bool _hydro;

  /*! @brief DensityGridWriterFields containing information about which
   *  output fields are active. */
  const DensityGridWriterFieldsMHD _fields;

  /*! @brief Log to write logging information to. */
  Log *_log;

public:
  /**
   * @brief Constructor.
   *
   * @param output_folder Name of the folder where output files should be
   * placed.
   * @param hydro Flag specifying whether or not hydro is active.
   * @param fields DensityGridWriterFields containing information about which
   * output fields are active.
   * @param log Log to write logging information to.
   */
  DensityGridWriterMHD(std::string output_folder, const bool hydro,
                    const DensityGridWriterFieldsMHD fields, Log *log)
      : _output_folder(output_folder), _hydro(hydro), _fields(fields),
        _log(log) {

    if (_log) {
      _log->write_status("Output will be written to ", _output_folder, "/");
    }
  }

  /**
   * @brief Virtual destructor.
   */
  virtual ~DensityGridWriterMHD() {}

  /**
   * @brief Write a snapshot.
   *
   * @param grid DensityGrid to write.
   * @param iteration Iteration number to use in the snapshot file name(s).
   * @param params ParameterFile containing the run parameters that should be
   * written to the file.
   * @param time Simulation time (in s).
   * @param hydro_units Internal unit system for the hydro.
   */
  virtual void write(DensityGridMHD &grid, uint_fast32_t iteration,
                     ParameterFile &params, double time = 0.,
                     const InternalMagnetoHydroUnits *hydro_units = nullptr) = 0;

  /**
   * @brief Write a snapshot for a split grid.
   *
   * @param grid_creator Grid.
   * @param counter Counter value to add to the snapshot file name.
   * @param params ParameterFile containing the run parameters that should be
   * written to the file.
   * @param time Simulation time (in s).
   */
  virtual void write(DensitySubGridCreatorMHD< DensitySubGridMHD > &grid_creator,
                     const uint_fast32_t counter, ParameterFile &params,
                     double time = 0.) = 0;

  /**
   * @brief Write a snapshot for a split grid with hydro.
   *
   * @param grid_creator Grid.
   * @param counter Counter value to add to the snapshot file name.
   * @param params ParameterFile containing the run parameters that should be
   * written to the file.
   * @param time Simulation time (in s).
   */
  virtual void write(DensitySubGridCreatorMHD< MagnetoHydroDensitySubGrid > &grid_creator,
                     const uint_fast32_t counter, ParameterFile &params,
                     double time = 0.) = 0;
};

#endif // DENSITYGRIDWRITER_HPP
