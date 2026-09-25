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
 * @file DiffusePhotonSourceDistribution.hpp
 *
 * @brief Distribution functor for photon sources.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DiffusePhotonSourceDistribution_HPP
#define DiffusePhotonSourceDistribution_HPP

#include "CoordinateVector.hpp"
#include "DensityGrid.hpp"
#include "HydroDensitySubGrid.hpp"
#include "DensitySubGridCreator.hpp"

/*! @brief Size of a variable that stores the number of photon sources. */
typedef uint_fast32_t photonsourcenumber_t;
class RandomGenerator;
class HydrogenLymanContinuumSpectrum;
class RecombinationRates;
class CrossSections;

/**
 * @brief General interface for photon source distribution functors.
 */
class DiffusePhotonSourceDistribution {
public:
  /**
   * @brief Virtual destructor.
   */
  virtual ~DiffusePhotonSourceDistribution() {}


  /**
   * @brief Get the total luminosity of all sources together.
   *
   * @return Total luminosity (in s^-1).
   */
  virtual double get_total_diffuse_luminosity() const = 0;

  virtual double get_total_diffuse_photons() const = 0;

  virtual photonsourcenumber_t get_number_of_diffuse_sources() const = 0;

  virtual double get_diffuse_photon_frequency(photonsourcenumber_t index) {
    return 0.0;
  };

  virtual double get_diffuse_weight(photonsourcenumber_t index) const = 0;

  virtual size_t get_diffuse_subgrid_index(photonsourcenumber_t index) const = 0;

  virtual bool has_active_diffuse_sources() const {
    return false;
  }

  

  /**
   * @brief Get the position of a diffuse photon source.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return CoordinateVector of a valid and photon source position (in m).
   */
  virtual CoordinateVector<> get_diffuse_position(photonsourcenumber_t index) = 0;

  /**
   * @brief Update the distribution after the system moved to the given time.
   *
   * @param simulation_time Current simulation time (in s).
   * @return True if the distribution changed, false otherwise.
   */

  virtual bool update(DensitySubGridCreator< HydroDensitySubGrid > *grid_creator, double actual_timestep) {return false;}


  /**
   * @brief Append distribution-specific metadata to a completed snapshot.
   *
   * Most distributions have no extra snapshot data. Implementations that do
   * write metadata should open and close the file entirely within this call.
   */
  virtual void write_snapshot_diffuse_metadata(const std::string &filename,
                                       const double simulation_time) {}

  /**
   * @brief Write the distribution to the given restart file.
   *
   * @param restart_writer RestartWriter to use.
   */
  virtual void write_restart_file(RestartWriter &restart_writer) const {
    cmac_error(
        "Restarting is not supported for this DiffusePhotonSourceDistribution!");
  }

};

#endif // DiffusePhotonSourceDistribution_HPP
