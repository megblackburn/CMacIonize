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


public:
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

  uint_fast32_t total_abs = _num_abs.value();
  uint_fast32_t total_escape = _num_escape.value();
  uint_fast32_t escape_ion = _num_escape_ionizing.value();
  uint_fast32_t escape_nonion = _num_escape_nonionizing.value();
  uint_fast32_t abs_dens = _num_abs_dens.value();
  uint_fast32_t abs_dif = _num_abs_dif.value();
  uint_fast32_t abs_dust = _num_abs_dust.value();
  uint_fast32_t reemitted_H = _num_reemitted_H.value();
  uint_fast32_t reemitted_He = _num_reemitted_He.value();
  uint_fast32_t abs_source = _num_abs_source.value();
  uint_fast32_t abs_reemitted_H = _num_abs_reemitted_H.value();
  uint_fast32_t abs_reemitted_He = _num_abs_reemitted_He.value();
  uint_fast32_t escape_ionizing_source = _num_escape_ionizing_source.value();
  uint_fast32_t escape_nonionizing_source = _num_escape_nonionizing_source.value();
  uint_fast32_t escape_reemitted_ionizing_H = _num_escape_reemitted_ionizing_H.value();
  uint_fast32_t escape_reemitted_nonionizing_H = _num_escape_reemitted_nonionizing_H.value();
  uint_fast32_t escape_reemitted_ionizing_He = _num_escape_reemitted_ionizing_He.value();
  uint_fast32_t escape_reemitted_nonionizing_He = _num_escape_reemitted_nonionizing_He.value();

  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalAbsorbed", total_abs);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalEscaped", total_escape);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalEscapedIonizing", escape_ion);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalEscapedNonIonizing", escape_nonion);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalAbsorbedDense", abs_dens);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalAbsorbedDiffuse", abs_dif);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalAbsorbedDust", abs_dust);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalReemittedHydrogen", reemitted_H);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalReemittedHelium", reemitted_He);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalAbsorbedSource", abs_source);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalAbsorbedReemittedHydrogen", abs_reemitted_H);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalAbsorbedReemittedHelium", abs_reemitted_He);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalEscapedIonizingSource", escape_ionizing_source);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalEscapedNonIonizingSource", escape_nonionizing_source);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalEscapedReemittedIonizingHydrogen", escape_reemitted_ionizing_H);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalEscapedReemittedNonIonizingHydrogen", escape_reemitted_nonionizing_H);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalEscapedReemittedIonizingHelium", escape_reemitted_ionizing_He);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "TotalEscapedReemittedNonIonizingHelium", escape_reemitted_nonionizing_He);

  double min_freq_val = min_frequency;
  double max_freq_val = max_frequency;
  uint_fast32_t num_bins_val = numbins;

  HDF5Tools::write_attribute< double >(photon_stats, "MinFrequency", min_freq_val);
  HDF5Tools::write_attribute< double >(photon_stats, "MaxFrequency", max_freq_val);
  HDF5Tools::write_attribute< uint_fast32_t >(photon_stats, "NumBins", num_bins_val);

  std::string freq_units = "Hz";
  HDF5Tools::write_attribute< std::string >(photon_stats, "FrequencyUnits", freq_units);

  std::vector< uint_fast32_t > ingoing_write_spectrum(numbins);
  std::vector< uint_fast32_t > outgoing_write_spectrum(numbins);
  for (uint_fast32_t i = 0; i < numbins; ++i) {
    ingoing_write_spectrum[i] = _ingoing_spectrum[i].value();
    outgoing_write_spectrum[i] = _outgoing_spectrum[i].value();
  }

  HDF5Tools::write_dataset< double >(photon_stats, "Frequencies", _frequencies);
  HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "IngoingSpectrum", ingoing_write_spectrum);
  HDF5Tools::write_dataset< uint_fast32_t >(photon_stats, "OutgoingSpectrum", outgoing_write_spectrum);

  HDF5Tools::close_group(photon_stats);
  HDF5Tools::close_file(file);
#else
    (void)filename;
#endif
}

};






#endif // PHOTONPACKETSTATISTICS_HPP

