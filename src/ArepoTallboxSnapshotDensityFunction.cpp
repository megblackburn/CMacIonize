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
 * @file ArepoTallboxSnapshotDensityFunction.cpp
 *
 * @brief Arepo Snapshot patch density function.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */

#include "ArepoTallboxSnapshotDensityFunction.hpp"
#include "CoordinateVector.hpp"
#include "DensityFunction.hpp"
#include "PhysicalConstants.hpp"
#include "HDF5Tools.hpp"
#include "Log.hpp"
#include "ParameterFile.hpp"
#include "Box.hpp"
#include <cinttypes>
#include <random>
#include <vector>
#include <algorithm>
#include <iostream>

#include <cmath>

/**
 * @brief Arepo Snapshot patch density function.
 *
 * Represents a gas density profile that is initialised from an arepo snapshot of tallbox region from a full galaxy simulation.
 */


  ArepoTallboxSnapshotDensityFunction::ArepoTallboxSnapshotDensityFunction(std::string filename, 
                                            const bool trace_initial_neutral_flag,
                                            const double temperature_to_trace,
                                            const double dust_gas_ratio,
                                            const double fraction_silicates,
                                            const bool zero_initial_xy_velocity,
                                            const Box<> box, 
                                            const CoordinateVector< uint_fast32_t > number_of_subgrids,
                                            const CoordinateVector< uint_fast32_t > number_of_cells,
                                            Log *log
                                  )
      : _filename(filename), _trace_initial_neutral_flag(trace_initial_neutral_flag), _temperature_to_trace(temperature_to_trace), _dust_gas_ratio(dust_gas_ratio), _fraction_silicates(fraction_silicates), _zero_initial_xy_velocity(zero_initial_xy_velocity), _box(box), _number_of_subgrids(number_of_subgrids), _ncell(number_of_cells), _log(log), _cartesian_grid(nullptr) {
        
        std::ifstream file(filename);
        
        if (!file.is_open()) {
            throw std::runtime_error("Could not open custom snapshot file: " + filename);
        }

      }

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
    ArepoTallboxSnapshotDensityFunction::ArepoTallboxSnapshotDensityFunction(ParameterFile &params, Log *log)
      : ArepoTallboxSnapshotDensityFunction(
            params.get_filename("DensityFunction:filename"),
            params.get_value< bool >("DensityFunction:trace initial neutral flag", false),
            params.get_physical_value< QUANTITY_TEMPERATURE >(
                "DensityFunction:temperature to trace", "500. K"),
            params.get_value< double >("DensityFunction:dust to gas",
                                      0.0),
            params.get_value< double >("DensityFunction:fraction silicates",
                                                    1.e-6),
            params.get_value< bool >("DensityFunction:zero intial xy velocity", false),
            Box<>(params.get_physical_vector< QUANTITY_LENGTH >(
                      "SimulationBox:anchor", "[-5. pc, -5. pc, -5. pc]"),
                  params.get_physical_vector< QUANTITY_LENGTH >(
                      "SimulationBox:sides", "[10. pc, 10. pc, 10. pc]")),
            params.get_value< CoordinateVector< uint_fast32_t > >(
                "DensitySubGridCreator:number of subgrids",
                CoordinateVector< uint_fast32_t >(8)),
            params.get_value< CoordinateVector< uint_fast32_t > >(
          "DensityGrid:number of cells", CoordinateVector< uint_fast32_t >(-1)),
            log) {} 

  /**
   * @brief Virtual destructor.
   */
   ArepoTallboxSnapshotDensityFunction::~ArepoTallboxSnapshotDensityFunction() {
        if (_cartesian_grid) {
          for (uint_fast32_t ix = 0; ix < _ncell.x(); ++ix) {
            for (uint_fast32_t iy = 0; iy < _ncell.y(); ++iy) {
              delete[] _cartesian_grid[ix][iy];
        }
        delete[] _cartesian_grid[ix];
        }
        delete[] _cartesian_grid;
    }
  } 


   void ArepoTallboxSnapshotDensityFunction::initialize() {

    if (_log) {
        _log->write_info("Opening the interpolated arepo tallbox snapshot: ", _filename);
    }

    std::cout << "Start of the Arepo Tallbox Snapshot Density Function" << std::endl;
    
    uint_fast32_t _nx = _ncell.x();
    uint_fast32_t _ny = _ncell.y();
    uint_fast32_t _nz = _ncell.z();
    

    HDF5Tools::HDF5File file = HDF5Tools::open_file(_filename, HDF5Tools::HDF5FILEMODE_READ);

   

    if (_log) {
            _log->write_info("Allocating task-based Cartesian grid structure dimensions...");
        }

   
    std::cout<<"number density set"<<std::endl;

    if (_log) {
            _log->write_info("Loading interpolated particle data layers from /PartType0...");
        }

    HDF5Tools::HDF5Group gas_group = HDF5Tools::open_group(file, "/PartType0");
    std::cout<<"parttype0 group opened"<<std::endl;

   // size_t total_elements = static_cast<size_t>(_nx) * _ny * _nz;
    //std::cout<<"total elements"<<std::endl;
      std::cout << "Allocating 3D nested data vectors..." << std::endl;

    // 1. Properly dimensioned 3D vectors to handle the [ix][iy][iz] syntax
    std::vector<std::vector<std::vector<double>>> cell_densities(_nx, std::vector<std::vector<double>>(_ny, std::vector<double>(_nz)));
    std::vector<std::vector<std::vector<double>>> cell_temperatures(_nx, std::vector<std::vector<double>>(_ny, std::vector<double>(_nz)));
    std::vector<std::vector<std::vector<double>>> cell_neutral_fractions(_nx, std::vector<std::vector<double>>(_ny, std::vector<double>(_nz)));
    std::vector<std::vector<std::vector<double>>> cell_velocity_x(_nx, std::vector<std::vector<double>>(_ny, std::vector<double>(_nz)));
    std::vector<std::vector<std::vector<double>>> cell_velocity_y(_nx, std::vector<std::vector<double>>(_ny, std::vector<double>(_nz)));
    std::vector<std::vector<std::vector<double>>> cell_velocity_z(_nx, std::vector<std::vector<double>>(_ny, std::vector<double>(_nz)));

    // 2. Temporarily read the flat stream from the file safely
    size_t total_elements = static_cast<size_t>(_nx) * _ny * _nz;
    std::vector<float> flat_dens(total_elements);
    std::vector<float> flat_temp(total_elements);
    std::vector<float> flat_nf(total_elements);
    std::vector<float> flat_vx(total_elements);
    std::vector<float> flat_vy(total_elements);
    std::vector<float> flat_vz(total_elements);

    hid_t raw_group_id = gas_group; 
    hid_t dset_dens = H5Dopen2(raw_group_id, "Density", H5P_DEFAULT);
    hid_t dset_temp = H5Dopen2(raw_group_id, "Temperature", H5P_DEFAULT);
    hid_t dset_nf   = H5Dopen2(raw_group_id, "NeutralFraction", H5P_DEFAULT);
    hid_t dset_vx   = H5Dopen2(raw_group_id, "Velocity_X", H5P_DEFAULT);
    hid_t dset_vy   = H5Dopen2(raw_group_id, "Velocity_Y", H5P_DEFAULT);
    hid_t dset_vz   = H5Dopen2(raw_group_id, "Velocity_Z", H5P_DEFAULT);

    // 2. USE H5T_NATIVE_FLOAT: This ensures a perfect byte-stride alignment
    H5Dread(dset_dens, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat_dens.data());
    H5Dclose(dset_dens);

    H5Dread(dset_temp, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat_temp.data());
    H5Dclose(dset_temp);

    H5Dread(dset_nf, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat_nf.data());
    H5Dclose(dset_nf);

    H5Dread(dset_vx, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat_vx.data());
    H5Dclose(dset_vx);

    H5Dread(dset_vy, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat_vy.data());
    H5Dclose(dset_vy);

    H5Dread(dset_vz, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat_vz.data());
    H5Dclose(dset_vz);
    HDF5Tools::close_group(gas_group);
    HDF5Tools::close_file(file);

    auto min_dens_it = std::min_element(flat_dens.begin(), flat_dens.end());
    auto max_dens_it = std::max_element(flat_dens.begin(), flat_dens.end());
    
    auto min_temp_it = std::min_element(flat_temp.begin(), flat_temp.end());
    auto max_temp_it = std::max_element(flat_temp.begin(), flat_temp.end());

    std::cout << "\n================ RAW FLAT HDF5 STREAM PROFILE ================" << std::endl;
    if (min_dens_it != flat_dens.end() && max_dens_it != flat_dens.end()) {
        std::cout << "Raw HDF5 Density  | Min: " << *min_dens_it << " | Max: " << *max_dens_it << " kg/m^3" << std::endl;
    }
    if (min_temp_it != flat_temp.end() && max_temp_it != flat_temp.end()) {
        std::cout << "Raw HDF5 Temp (K) | Min: " << *min_temp_it << " | Max: " << *max_temp_it << " K" << std::endl;
    }

    std::cout << "Populating 3D matrix variables..." << std::endl;

     for (uint_fast32_t ix = 0; ix < _nx; ++ix) {
        for (uint_fast32_t iy = 0; iy < _ny; ++iy) {
            for (uint_fast32_t iz = 0; iz < _nz; ++iz) {
                
                // Map the flat row-major array index sequentially: ix * (ny * nz) + iy * nz + iz
                const uint_fast32_t file_index = ix * (_ny * _nz) + iy * _nz + iz;

                // Sort the flat stream into your un-jumbled 3D arrays
                cell_densities[ix][iy][iz]          = flat_dens[file_index];
                cell_temperatures[ix][iy][iz]       = flat_temp[file_index];
                cell_neutral_fractions[ix][iy][iz]  = flat_nf[file_index];
                cell_velocity_x[ix][iy][iz]         = flat_vx[file_index];
                cell_velocity_y[ix][iy][iz]         = flat_vy[file_index];
                cell_velocity_z[ix][iy][iz]         = flat_vz[file_index];
            }
        }
    }

    const uint_fast32_t numblockx = _ncell.x() / _number_of_subgrids.x();
    const uint_fast32_t numblocky = _ncell.y() / _number_of_subgrids.y();
    const uint_fast32_t numblockz = _ncell.z() / _number_of_subgrids.z();
   // const uint_fast32_t numblocktot = numblockx * numblocky * numblockz;

    CoordinateVector< uint_fast32_t > numsubgrid = _number_of_subgrids;

    std::cout << "Data matrices ready for grid assignment loop." << std::endl;

    
  
   // const uint_fast32_t numblockx = _ncell.x() / numsubgrid.x(); // 128 / 16 = 8
    //const uint_fast32_t numblocky = _ncell.y() / numsubgrid.y(); // 128 / 16 = 8
    //const uint_fast32_t numblockz = _ncell.z() / numsubgrid.z(); // 768 / 16 = 48

    std::cout << "Assigning data using true 6-nested task-based geometry loops..." << std::endl;
    
    _cartesian_grid = new DensityValues **[_ncell.x()];
    for (uint_fast32_t ix = 0; ix < _ncell.x(); ++ix) {
      _cartesian_grid[ix] = new DensityValues *[_ncell.y()];
      for (uint_fast32_t iy = 0; iy < _ncell.y(); ++iy) {
        _cartesian_grid[ix][iy] = new DensityValues[_ncell.z()];
        for (uint_fast32_t iz = 0; iz < _ncell.z(); ++iz) {
          _cartesian_grid[ix][iy][iz].set_number_density(-1.);
        }
      }
    }

   

    for (uint_fast32_t six = 0; six < numsubgrid.x(); ++six) {
      for (uint_fast32_t siy = 0; siy < numsubgrid.y(); ++siy) {
        for (uint_fast32_t siz = 0; siz < numsubgrid.z(); ++siz) {
          
          for (uint_fast32_t cix = 0; cix < numblockx; ++cix) {
            for (uint_fast32_t ciy = 0; ciy < numblocky; ++ciy) {
              for (uint_fast32_t ciz = 0; ciz < numblockz; ++ciz) {
                
                const uint_fast32_t ix = six * numblockx + cix;
                const uint_fast32_t iy = siy * numblocky + ciy;
                const uint_fast32_t iz = siz * numblockz + ciz;

                // CHANGE THIS: Pull from your un-jumbled 3D arrays using spatial coordinates
                double mean_particle_mass = get_mean_particle_mass(cell_neutral_fractions[ix][iy][iz]);
                double number_dens = cell_densities[ix][iy][iz] / mean_particle_mass;

                _cartesian_grid[ix][iy][iz].set_number_density(number_dens);
                _cartesian_grid[ix][iy][iz].set_temperature(cell_temperatures[ix][iy][iz]);
                
                // mgb edit 12.06.2026: add neutral gas scalar field counter
                if (_trace_initial_neutral_flag == true){
                    if (_cartesian_grid[ix][iy][iz].get_temperature() <= _temperature_to_trace){
                        _cartesian_grid[ix][iy][iz].set_initial_neutral_scalar_field(1.0);
                        _cartesian_grid[ix][iy][iz].set_remaining_initial_neutral_scalar_field(1.0);
                    } else {
                       _cartesian_grid[ix][iy][iz].set_initial_neutral_scalar_field(0.0);
                       _cartesian_grid[ix][iy][iz].set_remaining_initial_neutral_scalar_field(0.0);
                    }
                  _cartesian_grid[ix][iy][iz].set_cooled_neutral_scalar_field(0.0);
                  _cartesian_grid[ix][iy][iz].set_remaining_cooled_neutral_scalar_field(0.0);
                  
                } else {
                    _cartesian_grid[ix][iy][iz].set_initial_neutral_scalar_field(0.);
                    _cartesian_grid[ix][iy][iz].set_remaining_initial_neutral_scalar_field(0.);
                    _cartesian_grid[ix][iy][iz].set_cooled_neutral_scalar_field(0.);
                    _cartesian_grid[ix][iy][iz].set_remaining_cooled_neutral_scalar_field(0.);

                } // end of mgb edit 12.06.2026
                _cartesian_grid[ix][iy][iz].set_fraction_silicates(_fraction_silicates);
                _cartesian_grid[ix][iy][iz].set_dust_gas_ratio(_dust_gas_ratio);
                
                _cartesian_grid[ix][iy][iz].set_ionic_fraction(ION_H_n, cell_neutral_fractions[ix][iy][iz]);

                if (flat_vx.size() > 0) {
                  if (_zero_initial_xy_velocity == true){
                    CoordinateVector<double> cell_velocity(
                        0.0, 
                        0.0,
                        cell_velocity_z[ix][iy][iz]
                    );
                    _cartesian_grid[ix][iy][iz].set_velocity(cell_velocity);
                  } else {
                    // Pull velocities from your un-jumbled 3D arrays
                    CoordinateVector<double> cell_velocity(
                        cell_velocity_x[ix][iy][iz],
                        cell_velocity_y[ix][iy][iz],
                        cell_velocity_z[ix][iy][iz]
                    );
                    _cartesian_grid[ix][iy][iz].set_velocity(cell_velocity);
                  }
                }
              }
            }
          }
        }
      }
    }

    for (uint_fast32_t ix = 0; ix < _ncell.x(); ++ix) {
      for (uint_fast32_t iy = 0; iy < _ncell.y(); ++iy) {
        for (uint_fast32_t iz = 0; iz < _ncell.z(); ++iz) {
          if (_cartesian_grid[ix][iy][iz].get_number_density() < 0.) {
            cmac_error("No values found for cell (%" PRIuFAST32 ", %" PRIuFAST32
                       ", %" PRIuFAST32 ")!",
                       ix, iy, iz);
          }
        }
      }
    }

    double min_final_temp = _cartesian_grid[0][0][0].get_temperature();
    double max_final_temp = min_final_temp;
    for (uint_fast32_t ix = 0; ix < _nx; ++ix) {
        for (uint_fast32_t iy = 0; iy < _ny; ++iy) {
            for (uint_fast32_t iz = 0; iz < _nz; ++iz) {
                double t = _cartesian_grid[ix][iy][iz].get_temperature();
                if (t < min_final_temp) min_final_temp = t;
                if (t > max_final_temp) max_final_temp = t;
            }
        }
    }
    std::cout << "[TRACE STAGE 3] FINAL ASSIGNED CARTESIAN GRID" << std::endl;
    std::cout << " -> _cartesian_grid   | Min Temp: " << min_final_temp << " K | Max Temp: " << max_final_temp << " K\n" << std::endl;
        
    size_t valid_cells = 0;
    size_t unassigned_cells = 0;
    double sample_dens = -1.0;
    double sample_temp = -1.0;

    for (uint_fast32_t ix = 0; ix < _nx; ++ix) {
        for (uint_fast32_t iy = 0; iy < _ny; ++iy) {
            for (uint_fast32_t iz = 0; iz < _nz; ++iz) {
                
                // 1. Check if the physical cell entry pointer structure is broken
                if (&_cartesian_grid[ix][iy][iz] == nullptr) {
                    std::cerr << "CRITICAL ERROR: Null memory block pointer discovered at (" 
                              << ix << "," << iy << "," << iz << ")!" << std::endl;
                    continue;
                }

                double d = _cartesian_grid[ix][iy][iz].get_number_density();
                double t = _cartesian_grid[ix][iy][iz].get_temperature();

                // 2. Check if a cell was skipped (retains the initial fallback -1.0 value)
                if (d == -1.0) {
                    unassigned_cells++;
                } else {
                    valid_cells++;
                    // Grab one real point as a sanity indicator
                    if (sample_dens == -1.0 && d > 0.0) {
                        sample_dens = d;
                        sample_temp = t;
                    }
                }
            }
        }
    }

    std::cout << "Grid Verification -> Total Valid Cells: " << valid_cells 
              << " | Mapped Skipped/Unassigned: " << unassigned_cells << std::endl;
    if (valid_cells > 0) {
        std::cout << "Sample Cell Density: " << sample_dens 
                  << " | Temp: " << sample_temp << " K" << std::endl;
    }



  }

/**
 * @brief Free up the memory used by the density function. After this,
 * operator() will no longer work.
 */
 void ArepoTallboxSnapshotDensityFunction::free() {
  if (_cartesian_grid) {
    for (uint_fast32_t ix = 0; ix < _ncell.x(); ++ix) {
      for (uint_fast32_t iy = 0; iy < _ncell.y(); ++iy) {
        delete[] _cartesian_grid[ix][iy];
      }
      delete[] _cartesian_grid[ix];
    }
    delete[] _cartesian_grid;
    _cartesian_grid = nullptr;
  }
}

/**
 * @brief Function that gives the density for a given cell.
 *
 * @param cell Geometrical information about the cell.
 * @return Initial physical field values for that cell.
 */
 DensityValues ArepoTallboxSnapshotDensityFunction::operator()(const Cell &cell) {

  const CoordinateVector<> position = cell.get_cell_midpoint();

  if (_cartesian_grid) {
    // get the indices of the cell containing the position
     int_fast32_t ix = _ncell.x() *
                             (position.x() - _box.get_anchor().x()) /
                             _box.get_sides().x();
     int_fast32_t iy = _ncell.y() *
                             (position.y() - _box.get_anchor().y()) /
                             _box.get_sides().y();
     int_fast32_t iz = _ncell.z() *
                             (position.z() - _box.get_anchor().z()) /
                             _box.get_sides().z();

    if (ix < 0) ix = 0;
    if (ix >= static_cast<int_fast32_t>(_ncell.x())) ix = _ncell.x() - 1;

    if (iy < 0) iy = 0;
    if (iy >= static_cast<int_fast32_t>(_ncell.y())) iy = _ncell.y() - 1;

    if (iz < 0) iz = 0;
    if (iz >= static_cast<int_fast32_t>(_ncell.z())) iz = _ncell.z() - 1;

    cmac_assert_message(ix < _ncell.x(), "%" PRIuFAST32, ix);
    cmac_assert_message(iy < _ncell.y(), "%" PRIuFAST32, iy);
    cmac_assert_message(iz < _ncell.z(), "%" PRIuFAST32, iz);

 //   std::cout<<"Assigning Cartesian Grid!"<<std::endl;
    return _cartesian_grid[ix][iy][iz];
  } else {
    cmac_error("This grid type is not supported (and you should never see this "
               "error)!");
    return DensityValues();
  }

};

