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
 * @file GadgetDensityGridWriter.hpp
 *
 * @brief HDF5-file writer for the DensityGrid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef GADGETDENSITYGRIDWRITERMHD_HPP
#define GADGETDENSITYGRIDWRITERMHD_HPP

#include "DensityGridWriterMHD.hpp"


#include <string>

class ParameterFile;

/**
 * @brief HDF5-file writer for the DensityGrid.
 */
class GadgetDensityGridWriterMHD : public DensityGridWriterMHD {
private:
  /*! @brief Prefix of the name for the file to write. */
  const std::string _prefix;

  /*! @brief Number of digits used for the counter in the filenames. */
  const uint_fast8_t _padding;

  /*! @brief Compress the HDF5 output? */
  const bool _compression;

public:
  GadgetDensityGridWriterMHD(
      std::string prefix, std::string output_folder = std::string("."),
      const bool hydro = false,
      const DensityGridWriterFieldsMHD fields = DensityGridWriterFieldsMHD(false),
      Log *log = nullptr, uint_fast8_t padding = 3,
      const bool compression = false);
  GadgetDensityGridWriterMHD(std::string output_folder, ParameterFile &params,
                          const bool hydro, Log *log = nullptr);

  virtual void write(DensityGridMHD &grid, uint_fast32_t iteration,
                     ParameterFile &params, double time = 0.,
                     const InternalMagnetoHydroUnits *hydro_units = nullptr);

  virtual void write(DensitySubGridCreatorMHD< DensitySubGridMHD > &grid_creator,
                     const uint_fast32_t counter, ParameterFile &params,
                     double time = 0.);

  virtual void write(DensitySubGridCreatorMHD< MagnetoHydroDensitySubGrid > &grid_creator,
                     const uint_fast32_t counter, ParameterFile &params,
                     double time = 0.);
};

#endif // GADGETDENSITYGRIDWRITERMHD_HPP
