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
 * @file TimeDependentDiffusePhotonSourceDistribution.hpp
 *
 * @brief TimeDependentDiffuse PhotonSourceDistribution.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#ifndef TIMEDEPENDENTDIFFUSEPHOTONSOURCEDISTRIBUTION_HPP
#define TIMEDEPENDENTDIFFUSEPHOTONSOURCEDISTRIBUTION_HPP

#include "PhysicalConstants.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "DiffusePhotonSourceDistribution.hpp"
#include "RandomGenerator.hpp"
#include "HydrogenLymanContinuumSpectrum.hpp"
#include "RecombinationRates.hpp"
#include "CrossSections.hpp"

#include <algorithm>
#include <cinttypes>
#include <fstream>
#include <unistd.h>
#include <vector>
#include <sys/stat.h>


/**
 * @brief Disc patch PhotonSourceDistribution.
 */
class TimeDependentDiffusePhotonSourceDistribution : public DiffusePhotonSourceDistribution {
private:

  std::vector< double> _snap_diffuse_source_luminosities;
  std::vector< CoordinateVector<> > _snap_diffuse_source_positions;
  std::vector< double> _snap_diffuse_photons;
  std::vector< double> _snap_diffuse_source_frequencies;
  
  std::vector< double> _diffuse_source_luminosities;
  std::vector< CoordinateVector<> > _diffuse_source_positions;
  std::vector< size_t > _diffuse_source_subgrid_indices;
  std::vector< double> _diffuse_photons;
  std::vector< double> _diffuse_source_frequencies;

  double _total_time = 0.;

  double _last_sf = 0.;

  double _total_diffuse_luminosity = 0.;
  double _total_diffuse_photons = 0.;

  /*! @brief RecombinationRates used to calculate ionic fractions. */
  const RecombinationRates &_recombination_rates;

  const CrossSections &_cross_sections;

  /*! @brief Hydrogen Lyman continuum spectrum. */
  const HydrogenLymanContinuumSpectrum _HLyc_spectrum;

  const bool _restart_flag;
  const double _restart_time;
  const double _update_interval;


  RandomGenerator _random_generator;







  Log *_log;
 

public:
  /**
   * @brief Constructor.
   */
  inline TimeDependentDiffusePhotonSourceDistribution(
      const RecombinationRates &recombination_rates,
      const CrossSections &cross_sections,
      const bool restart_flag, const double restart_time,
      const double update_interval,
      Log *log = nullptr)
      : _recombination_rates(recombination_rates), _cross_sections(cross_sections), _HLyc_spectrum(cross_sections), _restart_flag(restart_flag), _restart_time(restart_time), _update_interval(update_interval), _log(log){

    
    if (_restart_flag == true) {
      _total_time += _restart_time;
      _last_sf = _total_time;
      std::cout<< "total time: " << _total_time << " | restart time: " << _restart_time 
       << std::endl;
    }
  }



  // -----------------------------------------------

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  TimeDependentDiffusePhotonSourceDistribution(ParameterFile &params, RecombinationRates &recombination_rates, CrossSections &cross_sections, Log *log = nullptr)
      : TimeDependentDiffusePhotonSourceDistribution(recombination_rates, cross_sections, 
        params.get_value< bool >("TaskBasedRadiationHydrodynamicsSimulation:restart flag", false),
        params.get_physical_value< QUANTITY_TIME >("PhotonSourceDistribution:restart time", "0. Myr"),
        params.get_physical_value< QUANTITY_TIME >("PhotonSourceDistribution:update interval", "0.1 Myr"),
        log) {}


  /**
   * @brief Virtual destructor.
   */
  virtual ~TimeDependentDiffusePhotonSourceDistribution() {}

  


  /**
   * @brief Append live stars and supernovae since the previous snapshot.
   */
  virtual void write_snapshot_diffuse_metadata(const std::string &filename,
                                      double simulation_time) override {
#ifdef HAVE_HDF5
    if (filename.empty()) {
      return;
    }
    cmac_assert(_snap_diffuse_source_positions.size() == _snap_diffuse_source_luminosities.size());
    HDF5Tools::HDF5File file =
        HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_APPEND);

    HDF5Tools::HDF5Group diffuse_sources =
        HDF5Tools::create_group(file, "DiffusePhotonSources");
    uint32_t number_of_sources = _snap_diffuse_source_positions.size();
    std::string coordinate_units = "m";
    std::string luminosity_units = "s^-1";
    std::string photons_units = "";
    std::string frequency_units = "Hz";

    double total_diffuse_luminosity = get_total_diffuse_luminosity();

    HDF5Tools::write_attribute< uint32_t >(
      diffuse_sources, "NumberOfDiffuseSources", number_of_sources);
    HDF5Tools::write_attribute< double >(
      diffuse_sources, "SimulationTime", simulation_time);
    HDF5Tools::write_attribute< double >(
      diffuse_sources, "TotalDiffuseLuminosity", total_diffuse_luminosity);
    HDF5Tools::write_attribute< std::string >(
        diffuse_sources, "CoordinateUnits", coordinate_units);
    HDF5Tools::write_attribute< std::string >(
        diffuse_sources, "IonizingLuminosityUnits", luminosity_units);
    HDF5Tools::write_attribute< std::string >(
        diffuse_sources, "NumberOfPhotonsUnits", photons_units);
    HDF5Tools::write_attribute< std::string >(
        diffuse_sources, "FrequencyUnits", frequency_units);
    if (!_snap_diffuse_source_positions.empty()) {
      HDF5Tools::write_dataset< CoordinateVector<> >(
          diffuse_sources, "DiffuseCoordinates", _snap_diffuse_source_positions);
      HDF5Tools::write_dataset< double >(
          diffuse_sources, "DiffuseIonizingLuminosity", _snap_diffuse_source_luminosities);
      HDF5Tools::write_dataset< double >(
          diffuse_sources, "DiffuseNumberOfPhotons", _snap_diffuse_photons);
      HDF5Tools::write_dataset< double >(
          diffuse_sources, "DiffuseFrequency", _snap_diffuse_source_frequencies);
    }
    HDF5Tools::close_group(diffuse_sources);
    HDF5Tools::close_file(file);

    // Only clear after the HDF5 file was closed successfully.
    _snap_diffuse_source_positions.clear();
    _snap_diffuse_source_luminosities.clear();
    _snap_diffuse_photons.clear();
    _snap_diffuse_source_frequencies.clear();
#else
    (void)filename;
    (void)simulation_time;
#endif
  }


  /**
   * @brief Get the total luminosity of all sources together.
   *
   * @return Total luminosity (in s^-1).
   */


  virtual double get_total_diffuse_luminosity() const {
    return _total_diffuse_luminosity;
  }

    /**
   * @brief Get the total number of diffuse photons.
   *
   * @return Total number of diffuse photons.
   */


  virtual double get_total_diffuse_photons() const {
    return _total_diffuse_photons;
  }

    /**
   * @brief Get a valid position from the distribution.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return CoordinateVector of a valid and photon source position (in m).
   */
  virtual CoordinateVector<> get_diffuse_position(photonsourcenumber_t index) {
    return _diffuse_source_positions[index];
  }

    /**
   * @brief Get the number of sources contained within this distribution.
   *
   * The PhotonSourceDistribution will return exactly this number of valid
   * and unique positions by successive application of operator().
   *
   * @return Number of sources.
   */
  virtual photonsourcenumber_t get_number_of_diffuse_sources() const {
    return _diffuse_source_positions.size();
  }


  double get_diffuse_photon_frequency(photonsourcenumber_t index) {
    return _diffuse_source_frequencies[index];
  }

  /**
   * @brief Get the weight of a photon source.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return Weight of the photon source, used to determine how many photons are
   * emitted from this particular source.
   */
  virtual double get_diffuse_weight(photonsourcenumber_t index) const {
    return _total_diffuse_luminosity > 0.0 ?
        _diffuse_source_luminosities[index] / _total_diffuse_luminosity : 0.0;
  }

  virtual size_t get_diffuse_subgrid_index(photonsourcenumber_t index) const override {
    return _diffuse_source_subgrid_indices[index];
  }

  /**
   * @brief Update the distribution after the system moved to the given time.
   *
   * @param simulation_time Current simulation time (in s).
   * @return True if the distribution changed, false otherwise.
   */
       virtual bool update(DensitySubGridCreator< HydroDensitySubGrid > *grid_creator, 
                      double actual_timestep) override {

    _total_time += actual_timestep;
    bool updated = false;

    if (_total_time - _last_sf >= _update_interval) {
      const double delta_time = _total_time - _last_sf;

      _total_diffuse_luminosity = 0.0;
      _total_diffuse_photons = 0.0;

      _diffuse_source_luminosities.clear();
      _diffuse_source_positions.clear();
      _diffuse_source_subgrid_indices.clear();
      _diffuse_photons.clear();
      _diffuse_source_frequencies.clear();  

      _snap_diffuse_source_luminosities.clear();
      _snap_diffuse_source_positions.clear();
      _snap_diffuse_photons.clear();
      _snap_diffuse_source_frequencies.clear();

      const size_t num_subgrids = grid_creator->number_of_original_subgrids();
      const double EPSILON_PHOTONS = 1e-5;

      #pragma omp parallel
      {
        std::vector<double> local_lums;
        std::vector<CoordinateVector<>> local_positions;
        std::vector<size_t> local_subgrid_indices;
        std::vector<double> local_photons;
        std::vector<double> local_frequencies;
        double local_total_luminosity = 0.0;
        double local_total_photons = 0.0;

        unsigned int thread_seed = 42 + omp_get_thread_num();
        RandomGenerator local_rng;
        local_rng.set_seed(thread_seed);

        #pragma omp for schedule(dynamic, 1)
        for (size_t this_igrid = 0; this_igrid < num_subgrids; ++this_igrid) {
          HydroDensitySubGrid &subgrid = *grid_creator->get_subgrid(this_igrid);
          
          for (auto it = subgrid.hydro_begin(); it != subgrid.hydro_end(); ++it) {
            const double neutral_fraction =
                it.get_ionization_variables().get_ionic_fraction(ION_H_n);
            
            if (neutral_fraction >= 1.0) continue; 

            const double ionic_fraction = 1.0 - neutral_fraction;
            const double temperature = it.get_ionization_variables().get_temperature();
            const double number_density = it.get_ionization_variables().get_number_density();

            const double n_e = number_density * ionic_fraction;
            const double n_Hp = number_density * ionic_fraction; 

            const double alphaH = _recombination_rates.get_recombination_rate(ION_H_n, temperature);
            const double cell_volume = it.get_volume();

            const double diffuse_photon_rate = n_e * n_Hp * alphaH * cell_volume; 
            const double number_of_diffuse_photons = diffuse_photon_rate * delta_time;

            if (number_of_diffuse_photons < EPSILON_PHOTONS) continue;

            const double photon_frequency =
                _HLyc_spectrum.get_random_frequency(local_rng, temperature);
            const CoordinateVector<> cell_position = it.get_cell_midpoint();

            local_lums.push_back(diffuse_photon_rate);
            local_positions.push_back(cell_position);
            local_subgrid_indices.push_back(this_igrid);
            local_photons.push_back(number_of_diffuse_photons);
            local_frequencies.push_back(photon_frequency);
            local_total_luminosity += diffuse_photon_rate;
            local_total_photons += number_of_diffuse_photons;
          }
        }

        #pragma omp critical
        {
          _diffuse_source_luminosities.insert(_diffuse_source_luminosities.end(), local_lums.begin(), local_lums.end());
          _diffuse_source_positions.insert(_diffuse_source_positions.end(), local_positions.begin(), local_positions.end());
          _diffuse_source_subgrid_indices.insert(_diffuse_source_subgrid_indices.end(), local_subgrid_indices.begin(), local_subgrid_indices.end());
          _diffuse_photons.insert(_diffuse_photons.end(), local_photons.begin(), local_photons.end());
          _diffuse_source_frequencies.insert(_diffuse_source_frequencies.end(), local_frequencies.begin(), local_frequencies.end());
          _total_diffuse_luminosity += local_total_luminosity;
          _total_diffuse_photons += local_total_photons;
        }
      }
      std::cout<< "Time Dependent Diffuse photon source updated with " << _diffuse_source_positions.size() << " sources, total luminosity: " << _total_diffuse_luminosity << " s^-1, total photons: " << _total_diffuse_photons << std::endl;

      _snap_diffuse_source_luminosities = _diffuse_source_luminosities;
      _snap_diffuse_source_positions = _diffuse_source_positions;
      _snap_diffuse_photons = _diffuse_photons;
      _snap_diffuse_source_frequencies = _diffuse_source_frequencies;

      _last_sf = _total_time;
      updated = true;
    }

    updated = updated && get_total_diffuse_luminosity() > 0.0;
    return updated;
  }



  





// --------------------------------------

   /**
   * @brief Write the distribution variables to the restart file binary stream.
   */
  virtual void write_restart_file(RestartWriter &restart_writer) const override {
    restart_writer.write< double >(_total_time);
    restart_writer.write< double >(_last_sf);
    restart_writer.write< bool >(_restart_flag);
    restart_writer.write< double >(_restart_time);
    restart_writer.write< double >(_update_interval);
  }

  /**
   * @brief Restart constructor.
   *        
   */
  inline TimeDependentDiffusePhotonSourceDistribution(
      RestartReader &restart_reader,
      const RecombinationRates &recombination_rates,
      const CrossSections &cross_sections,
      Log *log = nullptr)
      : _recombination_rates(recombination_rates),
        _cross_sections(cross_sections),
        _HLyc_spectrum(cross_sections),
        _restart_flag(restart_reader.read< bool >()),      // FIXED: const initialized here!
        _restart_time(restart_reader.read< double >()),    // FIXED: const initialized here!
        _update_interval(restart_reader.read< double >()),  // FIXED: const initialized here!
        _log(log) {
    
    _total_time = restart_reader.read< double >();
    _last_sf    = restart_reader.read< double >();
  }


};
#endif // TIMEDEPENDENTDIFFUSEPHOTONSOURCEDISTRIBUTION_HPP

