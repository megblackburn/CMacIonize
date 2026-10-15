/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file PhotonPacketStatistics.hpp
 *
 * @brief Statistical information about the photon packets and re-emission
 * events.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef PHOTONPACKETSTATISTICS_HPP
#define PHOTONPACKETSTATISTICS_HPP

#include "AtomicValue.hpp"
#include "ParameterFile.hpp"
#include "PhotonPacket.hpp"
#include "HDF5Tools.hpp"
#include <fstream>
#include <vector>

/**
 * @brief Statistical information about the photon packets and re-emission
 * events.
 */
class PhotonPacketStatistics {
private:
  /**
   * @brief variable to store the number of photons on bins according to number
   * of scatterings.
   */
  std::vector< AtomicValue< uint_fast64_t > > _scatter_histogram;


  AtomicValue<uint_fast32_t> _num_abs;
  AtomicValue<uint_fast32_t> _num_escape;

  // mgb 22.09.2026 - add counters for escaping ionizing and non-ionizing photons
  AtomicValue<uint_fast32_t> _num_escape_ionizing;
  AtomicValue<uint_fast32_t> _num_escape_nonionizing;

  AtomicValue<uint_fast32_t> _num_abs_source;
  AtomicValue<uint_fast32_t> _num_escape_source;
  AtomicValue<uint_fast32_t> _num_abs_reemitted_H;
  AtomicValue<uint_fast32_t> _num_abs_reemitted_He;

  AtomicValue<uint_fast32_t> _num_escape_ionizing_source;
  AtomicValue<uint_fast32_t> _num_escape_nonionizing_source;
  AtomicValue<uint_fast32_t> _num_escape_reemitted_ionizing_H;
  AtomicValue<uint_fast32_t> _num_escape_reemitted_nonionizing_H;
  AtomicValue<uint_fast32_t> _num_escape_reemitted_ionizing_He;
  AtomicValue<uint_fast32_t> _num_escape_reemitted_nonionizing_He;

  AtomicValue<uint_fast32_t> _num_abs_dust;

  AtomicValue<uint_fast32_t> _num_abs_dens;
  AtomicValue<uint_fast32_t> _num_abs_dif;

  AtomicValue<uint_fast32_t> _num_reemitted_H;
  AtomicValue<uint_fast32_t> _num_reemitted_He;

  const uint_fast32_t _track_source_num;

  const uint_fast32_t numbins = 1000;

  std::vector<double> _frequencies;
  std::vector<AtomicValue<uint_fast32_t>> _ingoing_spectrum;
  std::vector<AtomicValue<uint_fast32_t>> _outgoing_spectrum;
  const double min_frequency = 3.289e15;
  const double max_frequency = 4. * min_frequency;
  const double _bin_width = (max_frequency - min_frequency) / numbins;


  // mgb 15.10.2026 - record each radiation step

  std::vector< uint_fast32_t > _history_absorbed;
  std::vector< uint_fast32_t > _history_escaped;
  std::vector< uint_fast32_t > _history_escaped_source;
  std::vector< uint_fast32_t > _history_escaped_ionizing;
  std::vector< uint_fast32_t > _history_escaped_nonionizing;
  std::vector< uint_fast32_t > _history_escaped_ionizing_source;
  std::vector< uint_fast32_t > _history_escaped_nonionizing_source;
  std::vector< uint_fast32_t > _history_reemitted_H;
  std::vector< uint_fast32_t > _history_reemitted_He;
  std::vector< uint_fast32_t > _history_reemitted_ionizing_H;
  std::vector< uint_fast32_t > _history_reemitted_nonionizing_H;
  std::vector< uint_fast32_t > _history_reemitted_ionizing_He;
  std::vector< uint_fast32_t > _history_reemitted_nonionizing_He;
  std::vector< uint_fast32_t > _history_absorbed_source;
  std::vector< uint_fast32_t > _history_absorbed_reemitted_H;
  std::vector< uint_fast32_t > _history_absorbed_reemitted_He;
  std::vector< uint_fast32_t > _history_absorbed_dust;
  std::vector< uint_fast32_t > _history_absorbed_dens;
  std::vector< uint_fast32_t > _history_absorbed_dif;

  std::vector<uint_fast32_t> _history_ingoing_spectra;
  std::vector<uint_fast32_t> _history_outgoing_spectra;
  
  uint_fast32_t _recorded_steps_count = 0;

public:


  inline void record_and_reset_radiation_step(){

    _history_absorbed.push_back(_num_abs.value());
    _history_escaped.push_back(_num_escape.value());
    _history_escaped_ionizing.push_back(_num_escape_ionizing.value());
    _history_escaped_nonionizing.push_back(_num_escape_nonionizing.value());
    _history_absorbed_source.push_back(_num_abs_source.value());
    _history_absorbed_reemitted_H.push_back(_num_abs_reemitted_H.value());
    _history_absorbed_reemitted_He.push_back(_num_abs_reemitted_He.value());
    _history_escaped_source.push_back(_num_escape_source.value());
    _history_escaped_ionizing_source.push_back(_num_escape_ionizing_source.value());
    _history_escaped_nonionizing_source.push_back(_num_escape_nonionizing_source.value());
    _history_reemitted_H.push_back(_num_reemitted_H.value());
    _history_reemitted_ionizing_H.push_back(_num_escape_reemitted_ionizing_H.value());
    _history_reemitted_nonionizing_H.push_back(_num_escape_reemitted_nonionizing_H.value());
    _history_reemitted_ionizing_He.push_back(_num_escape_reemitted_ionizing_He.value());
    _history_reemitted_nonionizing_He.push_back(_num_escape_reemitted_nonionizing_He.value());
    _history_reemitted_He.push_back(_num_reemitted_He.value());
    _history_absorbed_dust.push_back(_num_abs_dust.value());
    _history_absorbed_dens.push_back(_num_abs_dens.value());
    _history_absorbed_dif.push_back(_num_abs_dif.value());


     for (uint_fast32_t i = 0; i < numbins; ++i) {
      _history_ingoing_spectra.push_back(_ingoing_spectrum[i].value());
      _history_outgoing_spectra.push_back(_outgoing_spectrum[i].value());
    }
    
    _recorded_steps_count++;
    
    reset_counters();
  }

  /**
   * @brief constructor
   *
   * @param max_scatter maximum number of scatters recorded by the statistics
   */
  inline PhotonPacketStatistics(uint_fast32_t max_scatter, uint_fast32_t track_source_num)
      : _scatter_histogram(max_scatter + 2), _track_source_num(track_source_num),
       _frequencies(numbins),_ingoing_spectrum(numbins),
      _outgoing_spectrum(numbins){
        _num_abs.set(0);
        _num_escape.set(0);
        _num_escape_ionizing.set(0);
        _num_escape_nonionizing.set(0);
        _num_abs_dens.set(0);
        _num_abs_dif.set(0);
        _num_reemitted_H.set(0);
        _num_reemitted_He.set(0);

        for (uint_fast32_t i = 0; i < numbins; ++i) {
          _frequencies[i] =
            min_frequency + i * (max_frequency - min_frequency) /
                            (numbins - 1.);
          _ingoing_spectrum[i].set(0.0);
          _outgoing_spectrum[i].set(0.0);
        }
      }
  /**
   * @brief parameter file constructor
   *
   * These are the parameters that are used by this function:
   *   - maximum number of scatters; (default: 5)
   * @param params reference to the parameter file
   */
  inline PhotonPacketStatistics(ParameterFile &params)
      : PhotonPacketStatistics(params.get_value< uint_fast32_t >(
            "PhotonPacketStatistics:maximum number of scatters", 5),
            params.get_value< uint_fast32_t >(
            "PhotonPacketStatistics:track source index", 0)) {}
  /**
   * @brief function that implements photon absorption termination
   *
   * @param packet photon packet retrieved to be terminated
   */
  inline void absorb_photon(bool dense) {
    _num_abs.pre_increment();
    if (dense){
      _num_abs_dens.pre_increment();
    } else {
      _num_abs_dif.pre_increment();
    }
  }

  inline void absorb_photon_dust() {
      _num_abs_dust.pre_increment();
  }

  inline void absorb_photon(const PhotonPacket &packet) {
    _num_abs.pre_increment();
    size_t scatter_counter = packet.get_scatter_counter();
    _scatter_histogram[std::min(scatter_counter, _scatter_histogram.size() - 1)]
        .pre_increment();
  }

  inline void absorb_source_photon(const PhotonPacket &packet) {
    _num_abs.pre_increment();
    _num_abs_source.pre_increment();
    size_t scatter_counter = packet.get_scatter_counter();
    _scatter_histogram[std::min(scatter_counter, _scatter_histogram.size() - 1)]
        .pre_increment();
  }
  
  inline void absorb_reemitted_H_photon(const PhotonPacket &packet) {
    _num_abs.pre_increment();
    _num_abs_reemitted_H.pre_increment();
    size_t scatter_counter = packet.get_scatter_counter();
    _scatter_histogram[std::min(scatter_counter, _scatter_histogram.size() - 1)]
        .pre_increment();
  }

  inline void absorb_reemitted_He_photon(const PhotonPacket &packet) {
    _num_abs.pre_increment();
    _num_abs_reemitted_He.pre_increment();
    size_t scatter_counter = packet.get_scatter_counter();
    _scatter_histogram[std::min(scatter_counter, _scatter_histogram.size() - 1)]
        .pre_increment();
  }

  inline void injected_photon(const PhotonPacket &packet) {
    double photon_energy = packet.get_energy();
    if (photon_energy >= min_frequency && photon_energy < max_frequency) {
      // Determine the bin index for this energy
      uint_fast32_t bin_index = static_cast<uint_fast32_t>((photon_energy - min_frequency) / _bin_width);
      // Increment the count in the corresponding bin
      _ingoing_spectrum[bin_index].pre_increment();
    }
  }

  inline uint_fast32_t get_tracked_source() {return _track_source_num;}


  /**
   * @brief Function that implements photon escape termination
   *
   * @param packet Photon packet.
   */
  inline void escape_photon(const PhotonPacket &packet) {
    _num_escape.pre_increment();
    size_t scatter_counter = packet.get_scatter_counter();
    _scatter_histogram[std::min(scatter_counter, _scatter_histogram.size() - 1)]
        .pre_increment();
    double photon_energy = packet.get_energy();
    if (photon_energy >= min_frequency) {

      _num_escape_ionizing.pre_increment();
    
      if (photon_energy < max_frequency) {
        // Determine the bin index for this energy
        uint_fast32_t bin_index = static_cast<uint_fast32_t>((photon_energy - min_frequency) / _bin_width);
        // Increment the count in the corresponding bin
        _outgoing_spectrum[bin_index].pre_increment();
      }
    } else {
      _num_escape_nonionizing.pre_increment();
    }
  }

  /**
   * @brief Function that implements photon escape termination
   *
   * @param packet Photon packet.
   */
  inline void escape_reemitted_H_photon(const PhotonPacket &packet) {
    _num_escape.pre_increment();
    size_t scatter_counter = packet.get_scatter_counter();
    _scatter_histogram[std::min(scatter_counter, _scatter_histogram.size() - 1)]
        .pre_increment();
    double photon_energy = packet.get_energy();
    if (photon_energy >= min_frequency) {

      _num_escape_reemitted_ionizing_H.pre_increment();
    
      if (photon_energy < max_frequency) {
        // Determine the bin index for this energy
        uint_fast32_t bin_index = static_cast<uint_fast32_t>((photon_energy - min_frequency) / _bin_width);
        // Increment the count in the corresponding bin
        _outgoing_spectrum[bin_index].pre_increment();
      }
    } else {
      _num_escape_reemitted_nonionizing_H.pre_increment();
    }
  }

    /**
   * @brief Function that implements photon escape termination
   *
   * @param packet Photon packet.
   */
  inline void escape_reemitted_He_photon(const PhotonPacket &packet) {
    _num_escape.pre_increment();
    size_t scatter_counter = packet.get_scatter_counter();
    _scatter_histogram[std::min(scatter_counter, _scatter_histogram.size() - 1)]
        .pre_increment();
    double photon_energy = packet.get_energy();
    if (photon_energy >= min_frequency) {

      _num_escape_reemitted_ionizing_He.pre_increment();
    
      if (photon_energy < max_frequency) {
        // Determine the bin index for this energy
        uint_fast32_t bin_index = static_cast<uint_fast32_t>((photon_energy - min_frequency) / _bin_width);
        // Increment the count in the corresponding bin
        _outgoing_spectrum[bin_index].pre_increment();
      }
    } else {
      _num_escape_reemitted_nonionizing_He.pre_increment();
    }
  }

  /**
   * @brief Function that implements photon escape termination
   *
   * @param packet Photon packet.
   */
  inline void escape_source_photon(const PhotonPacket &packet) {
    _num_escape.pre_increment();
    _num_escape_source.pre_increment();
    size_t scatter_counter = packet.get_scatter_counter();
    _scatter_histogram[std::min(scatter_counter, _scatter_histogram.size() - 1)]
        .pre_increment();
    double photon_energy = packet.get_energy();
    if (photon_energy >= min_frequency) {

      _num_escape_ionizing_source.pre_increment();
    
      if (photon_energy < max_frequency) {
        // Determine the bin index for this energy
        uint_fast32_t bin_index = static_cast<uint_fast32_t>((photon_energy - min_frequency) / _bin_width);
        // Increment the count in the corresponding bin
        _outgoing_spectrum[bin_index].pre_increment();
      }
    } else {
      _num_escape_nonionizing_source.pre_increment();
    }
  }

  inline void reemit_H_photon() {
    _num_reemitted_H.pre_increment();
  }

  inline void reemit_He_photon() {
    _num_reemitted_He.pre_increment();
  }



  inline uint_fast32_t get_num_escaped() {
    return _num_escape.value();
  }

  inline uint_fast32_t get_num_escaped_source() {
    return _num_escape_source.value();
  }

  inline uint_fast32_t get_num_escaped_ionizing() {
    return _num_escape_ionizing.value();
  }

  inline uint_fast32_t get_num_escaped_nonionizing() {
    return _num_escape_nonionizing.value();
  }

  inline uint_fast32_t get_num_escaped_ionizing_source() {
    return _num_escape_ionizing_source.value();
  }

  inline uint_fast32_t get_num_escaped_nonionizing_source() {
    return _num_escape_nonionizing_source.value();
  }

  inline uint_fast32_t get_num_escaped_ionizing_H() {
    return _num_escape_reemitted_ionizing_H.value();
  }

  inline uint_fast32_t get_num_escaped_nonionizing_H() {
    return _num_escape_reemitted_nonionizing_H.value();
  }

  inline uint_fast32_t get_num_escaped_ionizing_He() {
    return _num_escape_reemitted_ionizing_He.value();
  }

  inline uint_fast32_t get_num_escaped_nonionizing_He() {
    return _num_escape_reemitted_nonionizing_He.value();
  }

  inline uint_fast32_t get_num_absorbed() {
    return _num_abs.value();
  }

  inline uint_fast32_t get_num_abs_source() {
    return _num_abs_source.value();
  }

  inline uint_fast32_t get_num_abs_reemitted_H() {
    return _num_abs_reemitted_H.value();
  }

  inline uint_fast32_t get_num_abs_reemitted_He() {
    return _num_abs_reemitted_He.value();
  }

  inline uint_fast32_t get_num_abs_dens() {
    return _num_abs_dens.value();
  }

  inline uint_fast32_t get_num_abs_dif() {
    return _num_abs_dif.value();
  }

  inline uint_fast32_t get_num_abs_dust() {
    return _num_abs_dust.value();
  }

  inline uint_fast32_t get_num_reemitted_H() {
    return _num_reemitted_H.value();
  }

  inline uint_fast32_t get_num_reemitted_He() {
    return _num_reemitted_He.value();
  }

  inline void reset_counters() {
    _num_abs.set(0);
    _num_escape.set(0);
    _num_escape_ionizing.set(0);
    _num_escape_nonionizing.set(0);
    _num_abs_dust.set(0);
    _num_abs_dens.set(0);
    _num_abs_dif.set(0);
    _num_reemitted_H.set(0);
    _num_reemitted_He.set(0);
    _num_abs_source.set(0);
    _num_abs_reemitted_H.set(0);
    _num_abs_reemitted_He.set(0);
    _num_escape_ionizing_source.set(0);
    _num_escape_nonionizing_source.set(0);
    _num_escape_reemitted_ionizing_H.set(0);
    _num_escape_reemitted_nonionizing_H.set(0);
    _num_escape_reemitted_ionizing_He.set(0);
    _num_escape_reemitted_nonionizing_He.set(0);
    for (uint_fast32_t i = 0; i < numbins; ++i) {
      _ingoing_spectrum[i].set(0.0);
      _outgoing_spectrum[i].set(0.0);
    }
  }
  /**
   * @brief function that outputs re-emission statistics of photons
   */
  inline void print_stats() {
// histogram thing
  {
    std::ofstream output_stats("photon_statistics.txt");
    output_stats << "# Scattering statistics for photons\n";
    output_stats << "# Nscatter\t BinCount  \n";
    for (uint_fast32_t i = 0; i < _scatter_histogram.size(); i++) {
      output_stats << i << "\t" << _scatter_histogram[i].value() << "\n";
    }
  }
// outgoing spectrum
  {
    std::ofstream output_stats("leaking_spectrum.txt");
    output_stats << "# Frequency\t Spectrum \n";
    for (uint_fast32_t i = 0; i < _frequencies.size(); i++) {
      output_stats << _frequencies[i] << "\t" << _outgoing_spectrum[i].value() << "\n";
    }
  }  
//ingoing spectrum
  {
    std::ofstream output_stats("input_spectrum.txt");
    output_stats << "# Frequency\t Spectrum \n";
    for (uint_fast32_t i = 0; i < _frequencies.size(); i++) {
      output_stats << _frequencies[i] << "\t" << _ingoing_spectrum[i].value() << "\n";
    }
  }

}

  
inline void write_snapshot_photon_statistics(const std::string &filename) {
#ifdef HAVE_HDF5
    if (filename.empty()) {
      return;
    }

    HDF5Tools::HDF5File file = HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_APPEND);

    HDF5Tools::HDF5Group photon_stats = HDF5Tools::create_group(file, "PhotonStatistics");

    // Write metadata attributes describing the dataset properties
    double min_freq_val = min_frequency;
    double max_freq_val = max_frequency;
    uint_fast32_t num_bins_val = numbins;
    uint_fast32_t steps_in_snapshot = _recorded_steps_count;

    HDF5Tools::write_attribute< double >(photon_stats, "MinFrequency", min_freq_val);
    HDF5Tools::write_attribute< double >(photon_stats, "MaxFrequency", max_freq_val);
    HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "NumBins", num_bins_val);
    HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "StepsInSnapshot", steps_in_snapshot);

    std::string freq_units = "Hz";
    HDF5Tools::write_attribute< std::string >(photon_stats, "FrequencyUnits", freq_units);

    // Save frequency bin axis definition
    HDF5Tools::write_dataset< double >(photon_stats, "Frequencies", _frequencies);

    if (steps_in_snapshot == 0) {
      _history_absorbed.push_back(0);
      _history_escaped.push_back(0);
      _history_escaped_ionizing.push_back(0);
      _history_escaped_nonionizing.push_back(0);
      _history_absorbed_source.push_back(0);
      _history_absorbed_reemitted_H.push_back(0);
      _history_absorbed_reemitted_He.push_back(0);
      _history_escaped_source.push_back(0);
      _history_escaped_ionizing_source.push_back(0);
      _history_escaped_nonionizing_source.push_back(0);
      _history_reemitted_H.push_back(0);
      _history_reemitted_ionizing_H.push_back(0);
      _history_reemitted_nonionizing_H.push_back(0);
      _history_reemitted_ionizing_He.push_back(0);
      _history_reemitted_nonionizing_He.push_back(0);
      _history_reemitted_He.push_back(0);
      _history_absorbed_dust.push_back(0);
      _history_absorbed_dens.push_back(0);
      _history_absorbed_dif.push_back(0);

      // Pad out spectra history vectors with zeros for 1 full dummy step block
      for (uint_fast32_t i = 0; i < numbins; ++i) {
        _history_ingoing_spectra.push_back(0);
        _history_outgoing_spectra.push_back(0);
      }
    }

    // Write out the step-resolved scalar histories securely
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalAbsorbed", _history_absorbed);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalEscaped", _history_escaped);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalEscapedIonizing", _history_escaped_ionizing);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalEscapedNonIonizing", _history_escaped_nonionizing);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalAbsorbedSource", _history_absorbed_source);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalAbsorbedReemittedHydrogen", _history_absorbed_reemitted_H);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalAbsorbedReemittedHelium", _history_absorbed_reemitted_He);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalEscapedSource", _history_escaped_source);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalEscapedIonizingSource", _history_escaped_ionizing_source);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalEscapedNonIonizingSource", _history_escaped_nonionizing_source);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalReemittedHydrogen", _history_reemitted_H);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalEscapedReemittedIonizingHydrogen", _history_reemitted_ionizing_H);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalEscapedReemittedNonIonizingHydrogen", _history_reemitted_nonionizing_H);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalEscapedReemittedIonizingHelium", _history_reemitted_ionizing_He);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalEscapedReemittedNonIonizingHelium", _history_reemitted_nonionizing_He);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalReemittedHelium", _history_reemitted_He);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalAbsorbedDust", _history_absorbed_dust);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalAbsorbedDense", _history_absorbed_dens);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryTotalAbsorbedDiffuse", _history_absorbed_dif);

    // Write out the flat step-resolved spectra histories
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryIngoingSpectrum", _history_ingoing_spectra);
    HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "HistoryOutgoingSpectrum", _history_outgoing_spectra);

    HDF5Tools::close_group(photon_stats);
    HDF5Tools::close_file(file);

    // Clear history cache containers so the next block doesn't accumulate leftovers
    _history_absorbed.clear();
    _history_escaped.clear();
    _history_escaped_ionizing.clear();
    _history_escaped_nonionizing.clear();
    _history_absorbed_source.clear();
    _history_absorbed_reemitted_H.clear();
    _history_absorbed_reemitted_He.clear();
    _history_escaped_source.clear();
    _history_escaped_ionizing_source.clear();
    _history_escaped_nonionizing_source.clear();
    _history_reemitted_H.clear();
    _history_reemitted_ionizing_H.clear();
    _history_reemitted_nonionizing_H.clear();
    _history_reemitted_ionizing_He.clear();
    _history_reemitted_nonionizing_He.clear();
    _history_reemitted_He.clear();
    _history_absorbed_dust.clear();
    _history_absorbed_dens.clear();
    _history_absorbed_dif.clear();
    _history_ingoing_spectra.clear();
    _history_outgoing_spectra.clear();
    _recorded_steps_count = 0;
#else
    (void)filename;
#endif
  }




};






#endif // PHOTONPACKETSTATISTICS_HPP

