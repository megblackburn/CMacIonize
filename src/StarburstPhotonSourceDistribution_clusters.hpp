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
 * @file StarburstPhotonSourceDistribution.hpp
 *
 * @brief Disc patch PhotonSourceDistribution.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#ifndef STARBURSTPHOTONSOURCEDISTRIBUTION_HPP
#define STARBURSTPHOTONSOURCEDISTRIBUTION_HPP

#include "Log.hpp"
#include "ParameterFile.hpp"
#include "PhotonSourceDistribution.hpp"
#include "RandomGenerator.hpp"
#include "DensitySubGridCreator.hpp"
#include "SupernovaHandler.hpp"
#include "WMBasicPhotonSourceSpectrum.hpp"
#include "PhotonSourceSpectrum.hpp"


#include <algorithm>
#include <cinttypes>
#include <fstream>
#include <unistd.h>
#include <vector>

/**
 * @brief Disc patch PhotonSourceDistribution.
 */
class StarburstPhotonSourceDistribution : public PhotonSourceDistribution {
private:
  /*! @brief Lifetime of a source (in s). */
  const double _source_lifetime;

  /*! @brief Ionising luminosity of a single source (in s^-1). */
  const double _source_luminosity;

  /*! @brief Probability of a new source being created (in s^-1). */
  const double _source_probability;

  /*! @brief Average number of sources at any given time. */
  const uint_fast32_t _average_number_of_sources;

  /*! @brief x component of the anchor of the rectangular disk (in m). */
  const double _anchor_x;

  /*! @brief y component of the anchor of the rectangular disk (in m). */
  const double _anchor_y;

  /*! @brief x side length of the rectangular disk (in m). */
  const double _sides_x;

  /*! @brief y side length of the rectangular disk (in m). */
  const double _sides_y;

  /*! @brief Origin of the Gaussian disk height distribution (in m). */
  const double _origin_z;

  /*! @brief Scale height of the Gaussian disk height distribution (in m). */
  const double _scaleheight_z;

  /*! @brief Update time interval (in s). */
  const double _update_interval;


  const double _cluster_luminosity;

  const double _cluster_mass;

  const double _length_of_burst;

  const double _time_of_burst_peak;

  const double _restart_time_myr;

  const bool _restart_flag;

  const double _star_formation_rate;

  const double _burst_amplitude;

  const bool _output_sources;

  std::string _source_filename;
  std::string _total_luminosity_filename;

  /*! @brief Pseudo-random number generator. */
  RandomGenerator _random_generator;

  /*! @brief Positions of the sources (in m). */
  std::vector< CoordinateVector<> > _source_positions;

  /*! @brief Remaining lifetime of the sources (in s). */
  std::vector< double > _source_lifetimes;

  /*! @brief Source Luminosity. */
  std::vector< double > _source_luminosities;

  /*! @brief Spectrum Index. */
  std::vector< double > _spectrum_index;


  /*! @brief Output file for the sources (if applicable). */
  std::ofstream *_output_file;

  /*! @brief Output file for the sources (if applicable). */
  std::ofstream *_output_file2;

  /*! @brief Number of updates since the start of the simulation. */
  uint_fast32_t _number_of_updates;

  /*! @brief Indices of the sources (if output is enabled). */
  std::vector< uint_fast32_t > _source_indices;

  /*! @brief spectrum index to use */
  size_t _starburst_spec_index;

  std::vector<PhotonSourceSpectrum*> _all_spectra;


  /*! @brief Index of the next source to add (if output is enabled). */
  uint_fast32_t _next_index;

  std::vector< CoordinateVector<> > _to_do_feedback;

  std::vector< double > _r_inj;
  std::vector< double > _r_st;
  std::vector< double > _num_cells_injected;
  std::vector< double > _nbar;


  const double _sne_energy = 1.e44;

  // defining constants: mgb edit 12.11.2025 
  const double unit_Myr = 3.15576e13;
  const double unit_yr = 3.154e+7;
  const double unit_kpc = 3.086e+19;
  const double unit_Msol = 1.988e30;

  SupernovaHandler *novahandler;

  /**
   * @brief Generate a new source position.
   *
   * @return New source position (in m).
   */
  inline CoordinateVector<> generate_source_position() {
    const double x =
        _anchor_x + _random_generator.get_uniform_random_double() * _sides_x;
    const double y =
        _anchor_y + _random_generator.get_uniform_random_double() * _sides_y;
    // we use the Box-Muller method to sample the Gaussian
    const double z =
        _scaleheight_z *
            std::sqrt(-2. *
                      std::log(_random_generator.get_uniform_random_double())) *
            std::cos(2. * M_PI *
                     _random_generator.get_uniform_random_double()) +
        _origin_z;

    return CoordinateVector<>(x, y, z);
  }




public:
  /**
   * @brief Constructor.
   *
   * @param source_lifetime Lifetime of a source (in s).
   * @param source_luminosity Ionising luminosity of a single source (in s^-1).
   * @param average_number Average number of sources at any given time.
   * @param anchor_x x component of the anchor of the rectangular disk (in m).
   * @param sides_x x side length of the rectangular disk (in m).
   * @param anchor_y  y component of the anchor of the rectangular disk (in m).
   * @param sides_y y side length of the rectangular disk (in m).
   * @param origin_z Origin of the Gaussian disk height distribution (in m).
   * @param scaleheight_z Scale height of the Gaussian disk height distribution
   * (in m).
   * @param seed Seed for the pseudo-random number generator.
   * @param update_interval Time interval in between successive source
   * distribution updates (in s).
   * @param starting_time Start time of the simulation. The distribution is
   * evolved forward in time to this point before it is used (in s).
   * @param output_sources Should the source positions be written to a file?
   */
  inline StarburstPhotonSourceDistribution(
      const double source_lifetime, const double source_luminosity,
      const uint_fast32_t average_number, const double anchor_x,
      const double sides_x, const double anchor_y, const double sides_y,
      const double origin_z, const double scaleheight_z,
      const int_fast32_t seed, const double update_interval,
      const double cluster_luminosity, const double cluster_mass,
      const double length_of_burst, const double time_of_burst_peak, const double restart_time_myr, const bool restart_flag,
      const double star_formation_rate, const double burst_amplitude,
      bool output_sources,
      const std::string source_filename="Starburst_source_positions.txt",
      const std::string total_luminosity_filename="TotalLuminosity.txt",
      const double sne_energy = 1.e44
      )
      : _source_lifetime(source_lifetime),
        _source_luminosity(source_luminosity),
        _source_probability(update_interval / source_lifetime),
        _average_number_of_sources(average_number), _anchor_x(anchor_x),
        _anchor_y(anchor_y), _sides_x(sides_x), _sides_y(sides_y),
        _origin_z(origin_z), _scaleheight_z(scaleheight_z),
        _update_interval(update_interval), _cluster_luminosity(cluster_luminosity), _cluster_mass(cluster_mass), 
        _length_of_burst(length_of_burst), _time_of_burst_peak(time_of_burst_peak), _restart_time_myr(restart_time_myr), _restart_flag(restart_flag),
        _star_formation_rate(star_formation_rate), _burst_amplitude(burst_amplitude), _output_sources(output_sources),
        _source_filename(source_filename), _total_luminosity_filename(total_luminosity_filename),
        _random_generator(seed),
        _output_file(nullptr), _output_file2(nullptr), _number_of_updates(1), _next_index(0),
        _sne_energy(sne_energy) {


    novahandler = new SupernovaHandler(_sne_energy);

    _all_spectra.push_back(new WMBasicPhotonSourceSpectrum(39000,25));
    _starburst_spec_index = _all_spectra.size()-1;
  

    if (_output_sources) {
          _output_file = new std::ofstream(_source_filename, std::ios_base::out | std::ios_base::app);
        if (_output_file->tellp() == 0) {
          *_output_file << "#time (s)\tx (m)\ty (m)\tz (m)\tevent\tindex\tluminosity\tMass\ttype\n";
        }
        _output_file2 = new std::ofstream(_total_luminosity_filename, std::ios_base::out | std::ios_base::app);
        if (_output_file2->tellp() == 0) {
          *_output_file2 << "simulation time (s)\tlum (s^-1)\tnumsne\tSFR (Msol Myr-1)\tSFR (Msol Myr^-1 kpc^-2)\n";
        }
    }


  }

  /**
   * @brief ParameterFile constructor.
   *
   * Parameters are:
   *  - source lifetime: Lifetime of a source (default: 20. Myr)
   *  - source luminosity: Ionising luminosity of a single source
   *    (default: 1.e48 s^-1)
   *  - average number of sources: Average number of sources (default: 24)
   *  - anchor x: X position of the anchor of the 2D disc (default: -1. kpc)
   *  - sides x: X side length of the 2D disc (default: 2. kpc)
   *  - anchor y: Y position of the anchor of the 2D disc (default: -1. kpc)
   *  - sides y: Y side length of the 2D disc (default: 2. kpc)
   *  - origin z: Origin of the exponential disc profile in the z direction
   *    (default: 0. pc)
   *  - scaleheight z: Vertical scale height of the exponential disc profile
   *    (default: 63. pc)
   *  - random seed: Random seed used to initialize the random generator that
   *    is used to sample the individual positions (default: 42)
   *  - update interval: Time interval in between successive distribution
   *    updates (default: 0.1 Myr)
   *  - starting time: Starting time of the simulation. The distribution is
   *    evolved forward in time to this point before it is used
   *    (default: 0. Myr)
   *  - output sources: Whether or not to write the source positions to a file
   *    (default: false)
   *
   * @param params ParameterFile to read from.
   * @param log Log to write logging info to.
   */
  StarburstPhotonSourceDistribution(ParameterFile &params, Log *log = nullptr)
      : StarburstPhotonSourceDistribution(
            params.get_physical_value< QUANTITY_TIME >(
                "PhotonSourceDistribution:source lifetime", "20. Myr"),
            params.get_physical_value< QUANTITY_FREQUENCY >(
                "PhotonSourceDistribution:source luminosity", "3.125e49 s^-1"),
            params.get_value< photonsourcenumber_t >(
                "PhotonSourceDistribution:average number of sources", 24),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:anchor x", "-1. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:sides x", "2. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:anchor y", "-1. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:sides y", "2. kpc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:origin z", "0. pc"),
            params.get_physical_value< QUANTITY_LENGTH >(
                "PhotonSourceDistribution:scaleheight z", "63. pc"),
            params.get_value< int_fast32_t >(
                "PhotonSourceDistribution:random seed", 42),
            params.get_physical_value< QUANTITY_TIME >(
                "PhotonSourceDistribution:update interval", "0.1 Myr"),
            params.get_physical_value< QUANTITY_FREQUENCY > (
                "PhotonSourceDistribution:cluster luminosity", "1.5e50 s^-1"
            ),
            params.get_physical_value< QUANTITY_MASS> (
                "PhotonSourceDistribution:cluster mass", "10000.0 Msol"
            ),
            params.get_physical_value< QUANTITY_TIME> (
                "PhotonSourceDistribution:length of burst", "0.0 Myr"
            ),
            params.get_physical_value< QUANTITY_TIME> (
                "PhotonSourceDistribution:time of burst peak", "0.0 Myr"
            ),
            params.get_physical_value< QUANTITY_TIME> (
                "PhotonSourceDistribution:restart time myr", "0.0 Myr"
            ),
            params.get_value< bool > (
                "PhotonSourceDistribution:restart flag", false
            ),
            params.get_physical_value< QUANTITY_MASS_RATE> (
                "PhotonSourceDistribution:star formation rate", "0.1 Msol yr^-1" // technically kpc^-2
            ),
            params.get_value< double> (
                "PhotonSourceDistribution:burst amplitude", 0.0
            ),
            params.get_value< bool >("PhotonSourceDistribution:output sources",
                                     false),
            params.get_value<std::string>("PhotonSourceDistribution:filename","SourceFile.txt"),
            params.get_value<std::string>("PhotonSourceDistribution:source filename","Bursty_source_positions.txt"),
            params.get_physical_value< QUANTITY_ENERGY > (
                "PhotonSourceDistribution:supernova energy", "1.e51 erg")) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~StarburstPhotonSourceDistribution() {}

  /**
   * @brief Get the number of sources contained within this distribution.
   *
   * The PhotonSourceDistribution will return exactly this number of valid
   * and unique positions by successive application of operator().
   *
   * @return Number of sources.
   */
  virtual photonsourcenumber_t get_number_of_sources() const {
    return _source_positions.size();
  }


  /**
   * @brief Will the distribution do stellar feedback at the given time?
   *
   * @param current_time Current simulation time (in s).
   * @return True if the star has not exploded yet and its lifetime has been
   * exceeded.
   */
  virtual bool do_stellar_feedback(const double current_time) const {
    return (_to_do_feedback.size() > 0);
  }



  virtual void get_sne_radii(DensitySubGridCreator< HydroDensitySubGrid > &grid_creator) {

       for (uint_fast32_t i = 0; i < _to_do_feedback.size(); ++i) {

        double r_inj,r_st,nbar,num_inj;


        std::tie(r_inj,r_st,nbar,num_inj) = novahandler->get_r_inj(&grid_creator,_to_do_feedback[i]);

         _r_inj.push_back(r_inj);
         _r_st.push_back(r_st);
         _nbar.push_back(nbar);
         _num_cells_injected.push_back(num_inj);

       }
  }



  virtual void add_stellar_feedback(HydroDensitySubGrid &subgrid, Hydro &hydro) {



    for (uint_fast32_t i = 0; i < _to_do_feedback.size(); ++i) {

      novahandler->inject_sne(subgrid, hydro, _to_do_feedback[i], _r_inj[i],_r_st[i],_nbar[i],_num_cells_injected[i]);

    }
  }

  virtual void done_stellar_feedback() {

    for (uint_fast32_t i=0; i<_to_do_feedback.size();i++) {

    std::cout << "\n SNe INJECTION HERE: R_inj = " << _r_inj[i] << " R_st = " <<  _r_st[i]
       << " num_cells = " <<  _num_cells_injected[i] << " nbar = "  << _nbar[i] << "\n";
    }



    _to_do_feedback.clear();
    _r_inj.clear();
    _r_st.clear();
    _num_cells_injected.clear();
    _nbar.clear();

  }



  /**
   * @brief Get a valid position from the distribution.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return CoordinateVector of a valid and photon source position (in m).
   */
  virtual CoordinateVector<> get_position(photonsourcenumber_t index) {
    return _source_positions[index];
  }

  
  /**
   * @brief Get the weight of a photon source.
   *
   * @param index Index of the photon source, must be in between 0 and
   * get_number_of_sources().
   * @return Weight of the photon source, used to determine how many photons are
   * emitted from this particular source.
   */
  virtual double get_weight(photonsourcenumber_t index) const {
    return 1. / get_number_of_sources();
  }

  /**
   * @brief Get the total luminosity of all sources together.
   *
   * @return Total luminosity (in s^-1).
   */
  virtual double get_total_luminosity() const {
    return _source_luminosity * get_number_of_sources();
  }

  double get_photon_frequency(RandomGenerator &random_generator,
    photonsourcenumber_t index) {

      return _all_spectra[_spectrum_index[index]]->get_random_frequency(random_generator,0.0);

}


  /**
   * @brief Update the distribution after the system moved to the given time.
   *
   * @param simulation_time Current simulation time (in s).
   * @return True if the distribution changed, false otherwise.
   */
  virtual bool update(DensitySubGridCreator< HydroDensitySubGrid > *grid_creator, const double simulation_time) {

    bool changed = false;
    while (_number_of_updates * _update_interval <= simulation_time) {

      const double total_time = _number_of_updates * _update_interval;
      // first clear out sources that do no longer exist
      size_t i = 0;
      while (i < _source_lifetimes.size()) {
        _source_lifetimes[i] -= _update_interval;
        if (_source_lifetimes[i] <= 0.) {
          // remove the element
          if (_output_file != nullptr) {
            *_output_file << total_time << "\t0.\t0.\t0.\t2\t"
                          << _source_indices[i] << "\n";
            _source_indices.erase(_source_indices.begin() + i);
          }

          _to_do_feedback.push_back(_source_positions[i]);
          _source_positions.erase(_source_positions.begin() + i);
          _source_lifetimes.erase(_source_lifetimes.begin() + i);
          _source_luminosities.erase(_source_luminosities.begin() + i);


          changed = true;
        } else {
          // check the next element
          ++i;
        }
      }

      double _bursty_star_formation_rate;


      double _sides_x = grid_creator->get_box().get_sides()[0];
      double _sides_y = grid_creator->get_box().get_sides()[1];


      double area_kpc = _sides_x*_sides_y/(unit_kpc*unit_kpc);
    
      double _total_time_for_sfr_burst;
      double _star_formation_rate_area = _star_formation_rate * area_kpc; // need to convert from Myr^-1 kpc^-2 to Myr^-1 mgb edit 11.11.2025 
      if (_restart_flag == true) {
          _total_time_for_sfr_burst = total_time + _restart_time_myr;
      }
      else {
          _total_time_for_sfr_burst = total_time;
      }

      _bursty_star_formation_rate = _star_formation_rate_area + (_burst_amplitude*_star_formation_rate_area - _star_formation_rate_area) * exp(-pow(((_total_time_for_sfr_burst - _time_of_burst_peak)/_length_of_burst), 2)); // edit mgb 19.09.2025 should be in units kg s^-1 kpc^-2
      std::cout << "SFR Amplitude = " << _bursty_star_formation_rate << _bursty_star_formation_rate*unit_Myr/unit_Msol << std::endl;
      std::cout << "Burst amplitude = " << _burst_amplitude << std::endl;
      std::cout << "Star formation rate = " << _star_formation_rate << std::endl;
      
        

        // 0.073 (now 0.207) factor is to take into account we only form stars over 8Msol
        // mass_to_generate in units of Msol to match IMF
        double mass_to_generate = 0.;
        
        mass_to_generate = _update_interval*_bursty_star_formation_rate/unit_Msol*0.207; // mgb edit 25.05.2026

        if (_output_file2 != nullptr) {
          double totallum = get_total_luminosity();
          *_output_file2 << total_time << "\t" << _total_time_for_sfr_burst << "\t" << totallum  << "\t" << ((_bursty_star_formation_rate*unit_Myr)/(unit_Msol)) << "\t" << ((_bursty_star_formation_rate*unit_Myr)/(area_kpc*unit_Msol)) << "\n"; // output the SFR in Msol Myr^-1 kpc^-2 and in kg s^-1 kpc^-2
          _output_file2->flush();

        }

        double clusters_required = mass_to_generate / _cluster_mass;

        while (clusters_required > 0.0){

          double sample = _random_generator.get_uniform_random_double();
          double cluster_creation_threshold = std::min(1.0, clusters_required);

          if (sample <= cluster_creation_threshold){

            double a0z = 9.955209529401348;
            double a1z = -3.3370109454102326;
            double a2z = 0.8116654874025604;

            const double lifetime_offset = _random_generator.get_uniform_random_double() * _update_interval;
            double mcut = 8.0;
            double lifetime = a0z + a1z*std::log10(mcut) + a2z*(std::log10(mcut)*std::log10(mcut));
            lifetime = std::pow(10.0,lifetime);
            lifetime = lifetime*3.154e+7;

            _source_lifetimes.push_back(lifetime - lifetime_offset);
            _source_positions.push_back(generate_source_position());
            _source_luminosities.push_back(_cluster_luminosity);
            _spectrum_index.push_back(_starburst_spec_index);

            if (_output_file != nullptr) {
              const CoordinateVector<> &pos = _source_positions.back();
              *_output_file << total_time << "\t" << pos.x() << "\t" << pos.y()
                        << "\t" << pos.z() << "\t1\t"
                        << _source_indices.back() << "\t"
                        << _source_luminosities.back() << "\t"
                        << _cluster_mass << "\t"
                        << "cluster\n";
        
                  ++_next_index;
              }
              changed = true;

          }
          clusters_required -= 1.0;

        }
        
      if (_output_file != nullptr) {
        _output_file->flush();
      }

      ++_number_of_updates;
    }

    return changed;
  }

  /**
   * @brief Write the distribution to the given restart file.
   *
   * @param restart_writer RestartWriter to use.
   */
  virtual void write_restart_file(RestartWriter &restart_writer) const {

    restart_writer.write(_source_lifetime);
    restart_writer.write(_source_luminosity);
    restart_writer.write(_source_probability);
    restart_writer.write(_average_number_of_sources);
    restart_writer.write(_anchor_x);
    restart_writer.write(_anchor_y);
    restart_writer.write(_sides_x);
    restart_writer.write(_sides_y);
    restart_writer.write(_origin_z);
    restart_writer.write(_scaleheight_z);
    restart_writer.write(_update_interval);
    restart_writer.write(_cluster_luminosity);
    restart_writer.write(_cluster_mass);
    restart_writer.write(_length_of_burst);
    restart_writer.write(_time_of_burst_peak);
    restart_writer.write(_restart_time_myr);
    restart_writer.write(_restart_flag);
    restart_writer.write(_output_sources);

    _random_generator.write_restart_file(restart_writer);
    {
      const auto size = _source_positions.size();
      restart_writer.write(size);
      for (std::vector< CoordinateVector<> >::size_type i = 0; i < size; ++i) {
        _source_positions[i].write_restart_file(restart_writer);
      }
    }
    {
      const auto size = _source_lifetimes.size();
      restart_writer.write(size);
      for (std::vector< double >::size_type i = 0; i < size; ++i) {
        restart_writer.write(_source_lifetimes[i]);
      }
    }
    restart_writer.write(_number_of_updates);
    const bool has_output = (_output_file != nullptr);
    restart_writer.write(has_output);
    if (has_output) {
      // store current position in the std::ofstream
      // we want to be able to continue writing from that point
      const auto filepos = _output_file->tellp();
      restart_writer.write(filepos);
      {
        const auto size = _source_indices.size();
        restart_writer.write(size);
        for (std::vector< uint_fast32_t >::size_type i = 0; i < size; ++i) {
          restart_writer.write(_source_indices[i]);
        }
      }
      restart_writer.write(_next_index);
    }
  }

  /**
   * @brief Restart constructor.
   *
   * @param restart_reader Restart file to read from.
   */
  inline StarburstPhotonSourceDistribution(RestartReader &restart_reader)
      : _source_lifetime(restart_reader.read< double >()),
        _source_luminosity(restart_reader.read< double >()),
        _source_probability(restart_reader.read< double >()),
        _average_number_of_sources(restart_reader.read< uint_fast32_t >()),
        _anchor_x(restart_reader.read< double >()),
        _anchor_y(restart_reader.read< double >()),
        _sides_x(restart_reader.read< double >()),
        _sides_y(restart_reader.read< double >()),
        _origin_z(restart_reader.read< double >()),
        _scaleheight_z(restart_reader.read< double >()),
        _update_interval(restart_reader.read< double >()),
        _cluster_luminosity(restart_reader.read< double >()),
        _cluster_mass(restart_reader.read< double >()),
        _length_of_burst(restart_reader.read< double >()),
        _time_of_burst_peak(restart_reader.read< double >()),
        _restart_time_myr(restart_reader.read< double >()),
        _restart_flag(restart_reader.read< double >()),
        _star_formation_rate(restart_reader.read< double>()),
        _burst_amplitude(restart_reader.read< double >()),
        _output_sources(restart_reader.read< double >()),
        _random_generator(restart_reader) {

    {
      const std::vector< CoordinateVector<> >::size_type size =
          restart_reader.read< std::vector< CoordinateVector<> >::size_type >();
      _source_positions.resize(size);
      for (std::vector< CoordinateVector<> >::size_type i = 0; i < size; ++i) {
        _source_positions[i] = CoordinateVector<>(restart_reader);
      }
    }
    {
      const std::vector< double >::size_type size =
          restart_reader.read< std::vector< double >::size_type >();
      _source_lifetimes.resize(size);
      for (std::vector< double >::size_type i = 0; i < size; ++i) {
        _source_lifetimes[i] = restart_reader.read< double >();
      }
    }
    _number_of_updates = restart_reader.read< uint_fast32_t >();
    const bool has_output = restart_reader.read< bool >();
    if (has_output) {
      const std::streampos filepos = restart_reader.read< std::streampos >();
      // truncate the original file to the size we were at
      if (truncate("Starburst_source_positions.txt", filepos) != 0) {
        cmac_error("Error while truncating output file!");
      }
      // now open the file in append mode
      _output_file = new std::ofstream("Starburst_source_positions.txt",
                                       std::ios_base::app);
      {
        const std::vector< uint_fast32_t >::size_type size =
            restart_reader.read< std::vector< uint_fast32_t >::size_type >();
        _source_indices.resize(size);
        for (std::vector< uint_fast32_t >::size_type i = 0; i < size; ++i) {
          _source_indices[i] = restart_reader.read< uint_fast32_t >();
        }
      }
      _next_index = restart_reader.read< uint_fast32_t >();
    }
    novahandler = new SupernovaHandler(_sne_energy);
  }
};

#endif // StarburstPHOTONSOURCEDISTRIBUTION_HPP
