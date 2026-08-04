/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file MagnetoHydroDensitySubGrid.hpp
 *
 * @brief Extension of DensitySubGrid that adds magneto-hydro variables.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#ifndef MAGNETOHYDRODENSITYSUBGRID_HPP
#define MAGNETOHYDRODENSITYSUBGRID_HPP

#include "DensitySubGridMHD.hpp"
#include "DensityValues.hpp"
#include "MagnetoHydro.hpp"
#include "MagnetoHydroVariables.hpp"
#include "TaskQueue.hpp"

/**
 * @brief Extension of DensitySubGrid that adds magneto-hydro variables.
 */
class alignas(16) MagnetoHydroDensitySubGrid : public DensitySubGridMHD {
private:
  /*! @brief Volume of a single cell (in m^3). */
  double _cell_volume;

  /*! @brief Inverse volume of a single cell (in m^-3). */
  double _inverse_cell_volume;

  /*! @brief Surface areas of a single cell (in m^2). */
  double _cell_areas[3];

  /*! @brief Hydrodynamical variables. */
 // MagnetoHydroVariables *_hydro_variables;

  /*! @brief Gradient limiters for the primitive hydrodynamical variables. */
  double *_primitive_variable_limiters;

  /*! @brief Indices of the hydro tasks associated with this subgrid. */
  size_t _hydro_tasks[18];

public:
  /**
   * @brief Constructor.
   *
   * @param box Dimensions of the box that contains the grid (in m; first 3
   * elements are the anchor of the box, 3 last elements are the side lengths
   * of the box).
   * @param ncell Number of cells in each dimension.
   */
  inline MagnetoHydroDensitySubGrid(const double *box,
                             const CoordinateVector< int_fast32_t > ncell)
      : DensitySubGridMHD(box, ncell),
        _cell_volume(_cell_size[0] * _cell_size[1] * _cell_size[2]),
        _inverse_cell_volume(1. / _cell_volume),
        _cell_areas{_cell_size[1] * _cell_size[2],
                    _cell_size[0] * _cell_size[2],
                    _cell_size[0] * _cell_size[1]} {

    // allocate memory for data arrays
    const int_fast32_t local_tot_ncell = _number_of_cells[3] * ncell[0];
   // const int_fast32_t local_tot_ncell = (size_t)ncell[0] * ncell[1] * ncell[2];
    _hydro_variables = new MagnetoHydroVariables[local_tot_ncell];
    _primitive_variable_limiters = new double[local_tot_ncell * 18]; // mgb edit 21.04.2026 from 10 to 18
        
    for (int_fast32_t i = 0; i < 9 * local_tot_ncell; ++i) { // mgb edit 21.04.2026 from 10 to 18
      _primitive_variable_limiters[2 * i] = DBL_MAX;
      _primitive_variable_limiters[2 * i + 1] = -DBL_MAX;
    }
    for (int i = 0; i < 18; ++i) {
        _hydro_tasks[i] = NO_TASK; 
    }
  }

  /**
   * @brief Copy constructor.
   *
   * @param original DensitySubGrid to copy.
   */
  inline MagnetoHydroDensitySubGrid(const MagnetoHydroDensitySubGrid &original)
      : DensitySubGridMHD(original), _cell_volume(original._cell_volume),
        _inverse_cell_volume(original._inverse_cell_volume),
        _cell_areas{original._cell_areas[0], original._cell_areas[1],
                    original._cell_areas[2]} {

    const int_fast32_t local_tot_ncell = _number_of_cells[3] * _number_of_cells[0];
   // const int_fast32_t local_tot_ncell = (size_t)_number_of_cells[0] * _number_of_cells[1] * _number_of_cells[2];
    _hydro_variables = new MagnetoHydroVariables[local_tot_ncell];
    _primitive_variable_limiters = new double[local_tot_ncell * 18];

    // copy data arrays
    for (int_fast32_t i = 0; i < local_tot_ncell; ++i) {
      _hydro_variables[i].copy_all(original._hydro_variables[i]);
    }

    for (int_fast32_t i = 0; i < 18 * local_tot_ncell; ++i) {
      _primitive_variable_limiters[i] =
          original._primitive_variable_limiters[i];
    }
    for (int i = 0; i < 18; ++i) {
        _hydro_tasks[i] = original._hydro_tasks[i];
    }
  }

  /**
   * @brief Destructor.
   */
  virtual ~MagnetoHydroDensitySubGrid() {
    // deallocate data arrays
    if (_hydro_variables != nullptr) {
      delete[] _hydro_variables;
      _hydro_variables = nullptr;
    }
    if (_primitive_variable_limiters != nullptr) {
      delete[] _primitive_variable_limiters;
      _primitive_variable_limiters = nullptr;
    }
  }

  /**
   * @brief Initialize the conserved variables for the grid.
   *
   * @param magnetohydro Hydro instance to use.
   * @param do_primitives Initialize the primitive variables based on the
   * ionization variables?
   * @return Minimum initial time step for the cells in the grid.
   */
  inline double initialize_hydrodynamic_variables(const MagnetoHydro &magnetohydro,
                                                  const bool do_primitives) {

    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    double timestep = DBL_MAX;
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      if (do_primitives) {
        magnetohydro.ionization_to_hydro(_ionization_variables[i],
                                  _hydro_variables[i]);
      }
      magnetohydro.set_conserved_variables(_hydro_variables[i], _cell_volume);
      timestep = std::min(timestep, magnetohydro.get_timestep(_hydro_variables[i],
                                                              _ionization_variables[i],
                                                              _cell_volume));
    }
    return timestep;
  }

  /**
   * @brief Update the conserved variables for all cells in the grid.
   *
   * @param timestep Integration time step size (in s).
   */
  inline void update_conserved_variables(const double timestep) {

    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      const CoordinateVector<> a =
          _hydro_variables[i].get_gravitational_acceleration();
      const CoordinateVector<> p = _hydro_variables[i].get_conserved_momentum();
      const double mdt = _hydro_variables[i].get_conserved_mass() * timestep;
      _hydro_variables[i].conserved(1) += mdt * a.x(); // mgb note - don't need edited as only gravitational effects not present on magnetic fields
      _hydro_variables[i].conserved(2) += mdt * a.y();
      _hydro_variables[i].conserved(3) += mdt * a.z();
      _hydro_variables[i].conserved(4) +=
          timestep * CoordinateVector<>::dot_product(p, a);
      _hydro_variables[i].conserved(4) += _hydro_variables[i].get_energy_term();
      _hydro_variables[i].set_energy_term(0.);
      for (int_fast8_t j = 0; j < 9; ++j) {
        _hydro_variables[i].conserved(j) +=
            _hydro_variables[i].delta_conserved(j) * timestep;

        // reset magnetohydro variables
        _hydro_variables[i].delta_conserved(j) = 0;
        _hydro_variables[i].primitive_gradients(j) = CoordinateVector<>(0.);
        _primitive_variable_limiters[18 * i + 2 * j] = DBL_MAX;
        _primitive_variable_limiters[18 * i + 2 * j + 1] = -DBL_MAX;
      }

      // mgb edit 21.04.2026 - check total energy + magnetic 
     // double pi = M_PI;
    //  const double mu0 = 4. * pi * 1e-7;
//      double Bx = _hydro_variables[i].get_primitives_magnetic_field().x();
   //   double By = _hydro_variables[i].get_primitives_magnetic_field().y();
   //   double Bz = _hydro_variables[i].get_primitives_magnetic_field().z();
    //  double B_squared = Bx*Bx + By*By + Bz*Bz;
    //  double magnetic_energy_density = B_squared/(2.*mu0);

//      double p2 = CoordinateVector<>::dot_product(_hydro_variables[i].get_conserved_momentum(), _hydro_variables[i].get_conserved_momentum());
      //double kinetic_energy_density = 0.5*p2/_hydro_variables[i].get_conserved_mass();
      // end of mgb edit 21.04.2026

      cmac_assert(_hydro_variables[i].get_conserved_mass() ==
                  _hydro_variables[i].get_conserved_mass());
      cmac_assert(_hydro_variables[i].get_conserved_momentum().x() ==
                  _hydro_variables[i].get_conserved_momentum().x());
      cmac_assert(_hydro_variables[i].get_conserved_momentum().y() ==
                  _hydro_variables[i].get_conserved_momentum().y());
      cmac_assert(_hydro_variables[i].get_conserved_momentum().z() ==
                  _hydro_variables[i].get_conserved_momentum().z());
      cmac_assert(_hydro_variables[i].get_conserved_total_energy() ==
                  _hydro_variables[i].get_conserved_total_energy());
      cmac_assert(_hydro_variables[i].get_conserved_total_energy() > 0.5*(CoordinateVector<>::dot_product(_hydro_variables[i].get_conserved_momentum(), _hydro_variables[i].get_conserved_momentum()))/_hydro_variables[i].get_conserved_mass() + (_hydro_variables[i].get_primitives_magnetic_field().x()*_hydro_variables[i].get_primitives_magnetic_field().x()+_hydro_variables[i].get_primitives_magnetic_field().y()*_hydro_variables[i].get_primitives_magnetic_field().y()+_hydro_variables[i].get_primitives_magnetic_field().z()*_hydro_variables[i].get_primitives_magnetic_field().z())/(2.*4. * M_PI * 1e-7)); // mgb edit 21.04.2026 - check that total energy is greater than kinetic + magnetic energy

#ifdef SAFE_HYDRO_VARIABLES
      _hydro_variables[i].conserved(0) =
          std::max(_hydro_variables[i].get_conserved_mass(), 0.);
      _hydro_variables[i].conserved(4) =
          std::max(_hydro_variables[i].get_conserved_total_energy(), 0.);
#else
      cmac_assert(_hydro_variables[i].get_conserved_mass() >= 0.);
      cmac_assert(_hydro_variables[i].get_conserved_total_energy() >= 0.);
#endif
    cmac_assert_message(_hydro_variables[i].get_conserved_mass() > 0.0, "about to set mass = 0");
    cmac_assert_message(_hydro_variables[i].get_conserved_total_energy() > 0.5*CoordinateVector<>::dot_product(_hydro_variables[i].get_conserved_momentum(), _hydro_variables[i].get_conserved_momentum())/_hydro_variables[i].get_conserved_mass(),
    "about to set kinetic greater than total energy...."); // mgb note - this may not be needed anymore as checking if > ekin + emag
    }
  }

  /**
   * @brief Update the primitive variables for the grid.
   *
   * @param magnetohydro Magnetohydro instance to use.
   */
  inline void update_primitive_variables(const MagnetoHydro &magnetohydro) {

    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      magnetohydro.set_primitive_variables(
          _hydro_variables[i], _ionization_variables[i], _inverse_cell_volume);
    }
  }

  /**
   * @brief Update the ionization variables for all cells in the subgrid using
   * their hydrodynamical variables.
   *
   * @param magnetohydro Magnetohydro instance to use.
   * @param maximum_neutral_fraction Maximum neutral fraction for hydrogen.
   */
  inline void
  update_ionization_variables(const MagnetoHydro &magnetohydro,
                              const double maximum_neutral_fraction) {

    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      magnetohydro.hydro_to_ionization(_hydro_variables[i], _ionization_variables[i]);
      if (maximum_neutral_fraction > 0. &&
          _ionization_variables[i].get_ionic_fraction(ION_H_n) >
              maximum_neutral_fraction) {
        _ionization_variables[i].set_ionic_fraction(ION_H_n,
                                                    maximum_neutral_fraction);
      }
    }
  }

  /**
   * @brief Add the energy because of photoionization to the hydro variables.
   *
   * @param magnetohydro Magnetohydro instance to use.
   * @param timestep Integration time step (in s).
   */
  inline void add_ionization_energy(const MagnetoHydro &magnetohydro, const double timestep) {

    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      magnetohydro.add_ionization_energy(_ionization_variables[i], _hydro_variables[i],
                                  _inverse_cell_volume, timestep);
    }
  }

  /**
   * @brief Half time step prediction for the primitive variables.
   *
   * @param magnetohydro Magnetohydro instance to use.
   * @param timestep Half system time step (in s).
   */
  inline void predict_primitive_variables(const MagnetoHydro &magnetohydro,
                                          const double timestep) {
    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      magnetohydro.predict_primitive_variables(_hydro_variables[i], timestep);
    }
  }

  /**
   * @brief Apply the slope limiter to all primitive variable gradients.
   *
   * @param magnetohydro Magnetohydro instance to use.
   */
  inline void apply_slope_limiter(const MagnetoHydro &magnetohydro) {
    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      magnetohydro.apply_slope_limiter(_hydro_variables[i],
                                        &_primitive_variable_limiters[18 * i], // mgb edit 21.04.2026 from 10 to 18
                                _cell_size);
    }
  }

  /**
   * @brief Compute the hydrodynamical fluxes for all interfaces inside the
   * subgrid.
   *
   * @param magnetohydro Magnetohydro instance to use.
   * @param dt Current system time step (in s).
   */
  inline void inner_flux_sweep(const MagnetoHydro &magnetohydro, const double dt) {

    // we do three separate sweeps: one for every coordinate direction
    for (int_fast32_t ix = 0; ix < _number_of_cells[0] - 1; ++ix) {
      for (int_fast32_t iy = 0; iy < _number_of_cells[1]; ++iy) {
        for (int_fast32_t iz = 0; iz < _number_of_cells[2]; ++iz) {
          const int_fast32_t index000 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          const int_fast32_t index100 =
              (ix + 1) * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          // x direction
          magnetohydro.do_flux_calculation(0, _hydro_variables[index000],
                                    _hydro_variables[index100], _cell_size[0],
                                    _cell_areas[0], dt);
        }
      }
    }
    for (int_fast32_t ix = 0; ix < _number_of_cells[0]; ++ix) {
      for (int_fast32_t iy = 0; iy < _number_of_cells[1] - 1; ++iy) {
        for (int_fast32_t iz = 0; iz < _number_of_cells[2]; ++iz) {
          const int_fast32_t index000 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          const int_fast32_t index010 =
              ix * _number_of_cells[3] + (iy + 1) * _number_of_cells[2] + iz;
          // y direction
          magnetohydro.do_flux_calculation(1, _hydro_variables[index000],
                                    _hydro_variables[index010], _cell_size[1],
                                    _cell_areas[1], dt);
        }
      }
    }
    for (int_fast32_t ix = 0; ix < _number_of_cells[0]; ++ix) {
      for (int_fast32_t iy = 0; iy < _number_of_cells[1]; ++iy) {
        for (int_fast32_t iz = 0; iz < _number_of_cells[2] - 1; ++iz) {
          const int_fast32_t index000 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          const int_fast32_t index001 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz + 1;
          // z direction
          magnetohydro.do_flux_calculation(2, _hydro_variables[index000],
                                    _hydro_variables[index001], _cell_size[2],
                                    _cell_areas[2], dt);
        }
      }
    }
  }

  /**
   * @brief Compute the hydrodynamical fluxes for all interfaces at the boundary
   * between this subgrid and the given neighbouring subgrid.
   *
   * @param direction TravelDirection of the neighbour.
   * @param magnetohydro Magnetohydro instance to use.
   * @param neighbour Neighbouring DensitySubGrid.
   * @param dt Current system time step (in s).
   */
  inline void outer_flux_sweep(const int_fast32_t direction, const MagnetoHydro &magnetohydro,
                               MagnetoHydroDensitySubGrid &neighbour,
                               const double dt) {

    int_fast32_t i, start_index_left, start_index_right, row_increment,
        row_length, column_increment, column_length;
    double dx, A;
    MagnetoHydroDensitySubGrid *left_grid, *right_grid;
    switch (direction) {
    case TRAVELDIRECTION_FACE_X_P:
      i = 0;
      left_grid = this;
      right_grid = &neighbour;
      start_index_left = (_number_of_cells[0] - 1) * _number_of_cells[3];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dx = _cell_size[0];
      A = _cell_areas[0];
      break;
    case TRAVELDIRECTION_FACE_X_N:
      i = 0;
      left_grid = &neighbour;
      right_grid = this;
      start_index_left = (_number_of_cells[0] - 1) * _number_of_cells[3];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dx = _cell_size[0];
      A = _cell_areas[0];
      break;
    case TRAVELDIRECTION_FACE_Y_P:
      i = 1;
      left_grid = this;
      right_grid = &neighbour;
      start_index_left = (_number_of_cells[1] - 1) * _number_of_cells[2];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = _cell_size[1];
      A = _cell_areas[1];
      break;
    case TRAVELDIRECTION_FACE_Y_N:
      i = 1;
      left_grid = &neighbour;
      right_grid = this;
      start_index_left = (_number_of_cells[1] - 1) * _number_of_cells[2];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = _cell_size[1];
      A = _cell_areas[1];
      break;
    case TRAVELDIRECTION_FACE_Z_P:
      i = 2;
      left_grid = this;
      right_grid = &neighbour;
      start_index_left = _number_of_cells[2] - 1;
      start_index_right = 0;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = _cell_size[2];
      A = _cell_areas[2];
      break;
    case TRAVELDIRECTION_FACE_Z_N:
      i = 2;
      left_grid = &neighbour;
      right_grid = this;
      start_index_left = _number_of_cells[2] - 1;
      start_index_right = 0;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = _cell_size[2];
      A = _cell_areas[2];
      break;
    default:
      cmac_error("Unknown magneto-hydro neighbour: %" PRIiFAST32, direction);
      break;
    }

    // using the index computation below is (much) faster than setting the
    // increment correctly and summing the indices manually
    for (int_fast32_t ic = 0; ic < column_length; ++ic) {
      for (int_fast32_t ir = 0; ir < row_length; ++ir) {
        const int_fast32_t index_left =
            start_index_left + ic * column_increment + ir * row_increment;
        const int_fast32_t index_right =
            start_index_right + ic * column_increment + ir * row_increment;
        magnetohydro.do_flux_calculation(i, left_grid->_hydro_variables[index_left],
                                  right_grid->_hydro_variables[index_right], dx,
                                  A, dt);
      }
    }
  }

  /**
   * @brief Compute the hydrodynamical fluxes for all interfaces at the boundary
   * between this subgrid and the given box boundary.
   *
   * @param direction TravelDirection of the neighbour.
   * @param magnetohydro Magnetohydro instance to use.
   * @param boundary MagnetohydroBoundary that sets the right state primitive
   * variables.
   * @param dt Current system time step (in s).
   */
  inline void outer_ghost_flux_sweep(const int_fast32_t direction,
                                     const MagnetoHydro &magnetohydro,
                                     const MagnetoHydroBoundary &boundary,
                                     const double dt) {

    int_fast32_t i, start_index_left, row_increment, row_length,
        column_increment, column_length;
    double dx, A;
    CoordinateVector<> offset;
    switch (direction) {
    case TRAVELDIRECTION_FACE_X_P:
      i = 0;
      start_index_left = (_number_of_cells[0] - 1) * _number_of_cells[3];
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dx = _cell_size[0];
      A = _cell_areas[0];
      offset = CoordinateVector<>(_cell_size[0], 0., 0.);
      break;
    case TRAVELDIRECTION_FACE_X_N:
      i = 0;
      start_index_left = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dx = -_cell_size[0];
      A = _cell_areas[0];
      offset = CoordinateVector<>(-_cell_size[0], 0., 0.);
      break;
    case TRAVELDIRECTION_FACE_Y_P:
      i = 1;
      start_index_left = (_number_of_cells[1] - 1) * _number_of_cells[2];
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = _cell_size[1];
      A = _cell_areas[1];
      offset = CoordinateVector<>(0., _cell_size[1], 0.);
      break;
    case TRAVELDIRECTION_FACE_Y_N:
      i = 1;
      start_index_left = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = -_cell_size[1];
      A = _cell_areas[1];
      offset = CoordinateVector<>(0., -_cell_size[1], 0.);
      break;
    case TRAVELDIRECTION_FACE_Z_P:
      i = 2;
      start_index_left = _number_of_cells[2] - 1;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = _cell_size[2];
      A = _cell_areas[2];
      offset = CoordinateVector<>(0., 0., _cell_size[2]);
      break;
    case TRAVELDIRECTION_FACE_Z_N:
      i = 2;
      start_index_left = 0;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dx = -_cell_size[2];
      A = _cell_areas[2];
      offset = CoordinateVector<>(0., 0., -_cell_size[2]);
      break;
    default:
      cmac_error("Unknown magneto-hydro neighbour: %" PRIiFAST32, direction);
      break;
    }

    // using the index computation below is (much) faster than setting the
    // increment correctly and summing the indices manually
    for (int_fast32_t ic = 0; ic < column_length; ++ic) {
      for (int_fast32_t ir = 0; ir < row_length; ++ir) {
        const int_fast32_t index_left =
            start_index_left + ic * column_increment + ir * row_increment;
        magnetohydro.do_ghost_flux_calculation(
            i, get_cell_midpoint(index_left) + offset,
            _hydro_variables[index_left], boundary, dx, A, dt);
      }
    }
  }

  /**
   * @brief Compute the hydrodynamical gradients for all interfaces inside the
   * subgrid.
   *
   * @param magnetohydro Magnetohydro instance to use.
   */
  inline void inner_gradient_sweep(const MagnetoHydro &magnetohydro) {

    // we do three separate sweeps: one for every coordinate direction
    for (int_fast32_t ix = 0; ix < _number_of_cells[0] - 1; ++ix) {
      for (int_fast32_t iy = 0; iy < _number_of_cells[1]; ++iy) {
        for (int_fast32_t iz = 0; iz < _number_of_cells[2]; ++iz) {
          const int_fast32_t index000 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          const int_fast32_t index100 =
              (ix + 1) * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          // x direction
          magnetohydro.do_gradient_calculation(
              0, _hydro_variables[index000], _hydro_variables[index100],
              _inv_cell_size[0], &_primitive_variable_limiters[18 * index000], // mgb edit 21.04.2026 from 10 to 18
              &_primitive_variable_limiters[18 * index100]); // mgb edit 21.04.2026 from 10 to 18
        }
      }
    }
    for (int_fast32_t ix = 0; ix < _number_of_cells[0]; ++ix) {
      for (int_fast32_t iy = 0; iy < _number_of_cells[1] - 1; ++iy) {
        for (int_fast32_t iz = 0; iz < _number_of_cells[2]; ++iz) {
          const int_fast32_t index000 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          const int_fast32_t index010 =
              ix * _number_of_cells[3] + (iy + 1) * _number_of_cells[2] + iz;
          // y direction
          magnetohydro.do_gradient_calculation(
              1, _hydro_variables[index000], _hydro_variables[index010],
              _inv_cell_size[1], &_primitive_variable_limiters[18 * index000], // mgb edit 21.04.2026 from 10 to 18
              &_primitive_variable_limiters[18 * index010]); // mgb edit 21.04.2026 from 10 to 18
        }
      }
    }
    for (int_fast32_t ix = 0; ix < _number_of_cells[0]; ++ix) {
      for (int_fast32_t iy = 0; iy < _number_of_cells[1]; ++iy) {
        for (int_fast32_t iz = 0; iz < _number_of_cells[2] - 1; ++iz) {
          const int_fast32_t index000 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          const int_fast32_t index001 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz + 1;
          // z direction
          magnetohydro.do_gradient_calculation(
              2, _hydro_variables[index000], _hydro_variables[index001],
              _inv_cell_size[2], &_primitive_variable_limiters[18 * index000], // mgb edit 21.04.2026 from 10 to 18
              &_primitive_variable_limiters[18 * index001]); // mgb edit 21.04.2026 from 10 to 18
        }
      }
    }
  }

  /**
   * @brief Compute the hydrodynamical gradients for all interfaces at the
   * boundary between this subgrid and the given neighbouring subgrid.
   *
   * @param direction TravelDirection of the neighbour.
   * @param magnetohydro Magnetohydro instance to use.
   * @param neighbour Neighbouring DensitySubGrid.
   */
  inline void outer_gradient_sweep(const int_fast32_t direction,
                                   const MagnetoHydro &magnetohydro,
                                   MagnetoHydroDensitySubGrid &neighbour) {

    int_fast32_t i, start_index_left, start_index_right, row_increment,
        row_length, column_increment, column_length;
    double dxinv;
    MagnetoHydroDensitySubGrid *left_grid, *right_grid;
    switch (direction) {
    case TRAVELDIRECTION_FACE_X_P:
      i = 0;
      left_grid = this;
      right_grid = &neighbour;
      start_index_left = (_number_of_cells[0] - 1) * _number_of_cells[3];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dxinv = _inv_cell_size[0];
      break;
    case TRAVELDIRECTION_FACE_X_N:
      i = 0;
      left_grid = &neighbour;
      right_grid = this;
      start_index_left = (_number_of_cells[0] - 1) * _number_of_cells[3];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dxinv = _inv_cell_size[0];
      break;
    case TRAVELDIRECTION_FACE_Y_P:
      i = 1;
      left_grid = this;
      right_grid = &neighbour;
      start_index_left = (_number_of_cells[1] - 1) * _number_of_cells[2];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = _inv_cell_size[1];
      break;
    case TRAVELDIRECTION_FACE_Y_N:
      i = 1;
      left_grid = &neighbour;
      right_grid = this;
      start_index_left = (_number_of_cells[1] - 1) * _number_of_cells[2];
      start_index_right = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = _inv_cell_size[1];
      break;
    case TRAVELDIRECTION_FACE_Z_P:
      i = 2;
      left_grid = this;
      right_grid = &neighbour;
      start_index_left = _number_of_cells[2] - 1;
      start_index_right = 0;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = _inv_cell_size[2];
      break;
    case TRAVELDIRECTION_FACE_Z_N:
      i = 2;
      left_grid = &neighbour;
      right_grid = this;
      start_index_left = _number_of_cells[2] - 1;
      start_index_right = 0;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = _inv_cell_size[2];
      break;
    default:
      cmac_error("Unknown magneto-hydro neighbour: %" PRIiFAST32, direction);
      break;
    }

    // using the index computation below is (much) faster than setting the
    // increment correctly and summing the indices manually
    for (int_fast32_t ic = 0; ic < column_length; ++ic) {
      for (int_fast32_t ir = 0; ir < row_length; ++ir) {
        const int_fast32_t index_left =
            start_index_left + ic * column_increment + ir * row_increment;
        const int_fast32_t index_right =
            start_index_right + ic * column_increment + ir * row_increment;
        magnetohydro.do_gradient_calculation(
            i, left_grid->_hydro_variables[index_left],
            right_grid->_hydro_variables[index_right], dxinv,
            &left_grid->_primitive_variable_limiters[18 * index_left], // mgb edit 21.04.2026 from 10 to 18
            &right_grid->_primitive_variable_limiters[18 * index_right]); // mgb edit 21.04.2026 from 10 to 18
      }
    }
  }

  /**
   * @brief Compute the hydrodynamical gradients for all interfaces at the
   * boundary between this subgrid and the given box boundary with boundary
   * condition.
   *
   * @param direction TravelDirection of the neighbour.
   * @param magnetohydro Magnetohydro instance to use.
   * @param boundary MagnetohydroBoundary that sets the right state primitive
   * variables.
   */
  inline void outer_ghost_gradient_sweep(const int_fast32_t direction,
                                         const MagnetoHydro &magnetohydro,
                                         const MagnetoHydroBoundary &boundary) {

    int_fast32_t i, start_index_left, row_increment, row_length,
        column_increment, column_length;
    double dxinv;
    MagnetoHydroDensitySubGrid *left_grid;
    CoordinateVector<> offset;
    switch (direction) {
    case TRAVELDIRECTION_FACE_X_P:
      i = 0;
      left_grid = this;
      start_index_left = (_number_of_cells[0] - 1) * _number_of_cells[3];
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dxinv = _inv_cell_size[0];
      offset = CoordinateVector<>(_cell_size[0], 0., 0.);
      break;
    case TRAVELDIRECTION_FACE_X_N:
      i = 0;
      left_grid = this;
      start_index_left = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[2];
      column_length = _number_of_cells[1];
      dxinv = -_inv_cell_size[0];
      offset = CoordinateVector<>(-_cell_size[0], 0., 0.);
      break;
    case TRAVELDIRECTION_FACE_Y_P:
      i = 1;
      left_grid = this;
      start_index_left = (_number_of_cells[1] - 1) * _number_of_cells[2];
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = _inv_cell_size[1];
      offset = CoordinateVector<>(0., _cell_size[1], 0.);
      break;
    case TRAVELDIRECTION_FACE_Y_N:
      i = 1;
      left_grid = this;
      start_index_left = 0;
      row_increment = 1;
      row_length = _number_of_cells[2];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = -_inv_cell_size[1];
      offset = CoordinateVector<>(0., -_cell_size[1], 0.);
      break;
    case TRAVELDIRECTION_FACE_Z_P:
      i = 2;
      left_grid = this;
      start_index_left = _number_of_cells[2] - 1;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = _inv_cell_size[2];
      offset = CoordinateVector<>(0., 0., _cell_size[2]);
      break;
    case TRAVELDIRECTION_FACE_Z_N:
      i = 2;
      left_grid = this;
      start_index_left = 0;
      row_increment = _number_of_cells[2];
      row_length = _number_of_cells[1];
      column_increment = _number_of_cells[3];
      column_length = _number_of_cells[0];
      dxinv = -_inv_cell_size[2];
      offset = CoordinateVector<>(0., 0., -_cell_size[2]);
      break;
    default:
      cmac_error("Unknown magneto-hydro neighbour: %" PRIiFAST32, direction);
      break;
    }

    // using the index computation below is (much) faster than setting the
    // increment correctly and summing the indices manually
    for (int_fast32_t ic = 0; ic < column_length; ++ic) {
      for (int_fast32_t ir = 0; ir < row_length; ++ir) {
        const int_fast32_t index_left =
            start_index_left + ic * column_increment + ir * row_increment;
        magnetohydro.do_ghost_gradient_calculation(
            i, get_cell_midpoint(index_left) + offset,
            left_grid->_hydro_variables[index_left], boundary, dxinv,
            &left_grid->_primitive_variable_limiters[18 * index_left]); // mgb edit 21.04.2026 from 10 to 18
      }
    }
  }

  /**
   * @brief Set the hydro task with the given index.
   *
   * @param i Index.
   * @param task Task.
   */
  inline void set_hydro_task(const int_fast32_t i, const size_t task) {
    _hydro_tasks[i] = task;
  }

  /**
   * @brief Get the magnetohydro task with the given index.
   *
   * @param i Index.
   * @return Task.
   */
  inline size_t get_hydro_task(const int_fast32_t i) const {
    return _hydro_tasks[i];
  }

  /**
   * @brief Initialize the hydrodynamic variables for a cell in this subgrid.
   *
   * @param index Index of a cell in the subgrid.
   * @param values DensityValues to use.
   */
  virtual void initialize_hydro(const uint_fast32_t index,
                                const DensityValues &values) {
     if (index >= DensitySubGridMHD::get_number_of_cells()) { 
        std::cerr << "OUT OF BOUNDS!" << std::endl;
        return;
    }

    _hydro_variables[index].set_primitives_velocity(values.get_velocity());
    _hydro_variables[index].set_primitives_magnetic_field(values.get_magnetic_field());
  }

  /**
   * @brief Iterator to loop over the cells in the subgrid.
   */
  class hydroiterator : public Cell {
  private:
    /*! @brief Index of the cell the iterator is currently pointing to. */
    uint_fast32_t _index;

    /*! @brief Pointer to the underlying subgrid (we cannot use a reference,
     *  since then things like it = it would not work). */
    MagnetoHydroDensitySubGrid *_subgrid;

  public:
    /**
     * @brief Constructor.
     *
     * @param index Index of the cell the iterator is currently pointing to.
     * @param subgrid MagnetoHydroDensitySubGrid over which we iterate.
     */
    inline hydroiterator(const uint_fast32_t index,
                         MagnetoHydroDensitySubGrid &subgrid)
        : _index(index), _subgrid(&subgrid) {}

    // Cell interface

    /**
     * @brief Get the midpoint of the cell the iterator is pointing to.
     *
     * @return Cell midpoint (in m).
     */
    virtual CoordinateVector<> get_cell_midpoint() const {
      return _subgrid->get_cell_midpoint(_index);
    }

    

    /**
     * @brief Get the volume of the cell the iterator is pointing to.
     *
     * @return Cell volume (in m^3).
     */
    virtual double get_volume() const { return _subgrid->_cell_volume; }

    /**
     * @brief Get the faces of the cell.
     *
     * @return Faces of the cell.
     */
    virtual std::vector< Face > get_faces() const {
      return std::vector< Face >();
    }

    // HydroDensitySubGrid access functionality

    /**
     * @brief Get read only access to the magnetohydro variables stored in this
     * cell.
     *
     * @return Read only access to the magnetohydro variables.
     */
    inline const MagnetoHydroVariables &get_hydro_variables() const {
      return _subgrid->_hydro_variables[_index];
    }

    /**
     * @brief Get read/write access to the magnetohydro variables stored in this
     * cell.
     *
     * @return Read/write access to the magnetohydro variables.
     */
    inline MagnetoHydroVariables &get_hydro_variables() {
      return _subgrid->_hydro_variables[_index];
    }

    /**
     * @brief Get read only access to the ionization variables stored in this
     * cell.
     *
     * @return Read only access to the ionization variables.
     */
    inline const IonizationVariables &get_ionization_variables() const {
      return _subgrid->_ionization_variables[_index];
    }

    /**
     * @brief Get read/write access to the ionization variables stored in this
     * cell.
     *
     * @return Read/write access to the ionization variables.
     */
    inline IonizationVariables &get_ionization_variables() {
      return _subgrid->_ionization_variables[_index];
    }

    // Iterator functionality

    /**
     * @brief Increment operator.
     *
     * We only implemented the pre-increment version, since the post-increment
     * version creates a new object and is computationally more expensive.
     *
     * @return Reference to the incremented iterator.
     */
    inline hydroiterator &operator++() {
      ++_index;
      return *this;
    }

    /**
     * @brief Increment operator.
     *
     * @param increment Increment to add.
     * @return Reference to the incremented iterator.
     */
    inline hydroiterator &operator+=(const uint_fast32_t increment) {
      _index += increment;
      return *this;
    }

    /**
     * @brief Free addition operator.
     *
     * @param increment Increment to add to the iterator.
     * @return Incremented iterator.
     */
    inline hydroiterator operator+(const uint_fast32_t increment) const {
      hydroiterator it(*this);
      it += increment;
      return it;
    }

    /**
     * @brief Get the index of the cell the iterator is currently pointing to.
     *
     * @return Index of the current cell.
     */
    inline uint_fast32_t get_index() const { return _index; }

    /**
     * @brief Compare iterators.
     *
     * @param it Iterator to compare with.
     * @return True if the iterators point to the same cell of the same grid.
     */
    inline bool operator==(hydroiterator it) const {
      return (_subgrid == it._subgrid && _index == it._index);
    }

    /**
     * @brief Compare iterators.
     *
     * @param it Iterator to compare with.
     * @return True if the iterators do not point to the same cell of the same
     * grid.
     */
    inline bool operator!=(hydroiterator it) const { return !(*this == it); }
  };

  /**
   * @brief Get an iterator to the first cell in the subgrid.
   *
   * @return Iterator to the first cell in the subgrid.
   */
  inline hydroiterator hydro_begin() { return hydroiterator(0, *this); }

  /**
   * @brief Get an iterator to the beyond last cell in the subgrid.
   *
   * @return Iterator to the beyond last cell in the subgrid.
   */
  inline hydroiterator hydro_end() {
    return hydroiterator(
        _number_of_cells[0] * _number_of_cells[1] * _number_of_cells[2], *this);
  }

  /**
   * @brief Get an iterator to the cell that contains the given position.
   *
   * @param position Position (in m).
   * @return Iterator to the corresponding cell.
   */
  inline hydroiterator get_hydro_cell(const CoordinateVector<> position) {
    CoordinateVector< int_fast32_t > three_index;
    return hydroiterator(get_start_index(position - _anchor,
                                                TRAVELDIRECTION_INSIDE, three_index),
                                *this);
  }

  /**
   * @brief Dump the subgrid to the given restart file.
   *
   * @param restart_writer RestartWriter to write to.
   */
  virtual void write_restart_file(RestartWriter &restart_writer) const {

    DensitySubGridMHD::write_restart_file(restart_writer);

    restart_writer.write(_cell_volume);
    restart_writer.write(_cell_areas[0]);
    restart_writer.write(_cell_areas[1]);
    restart_writer.write(_cell_areas[2]);
    const int_fast32_t number_of_cells =
        _number_of_cells[0] * _number_of_cells[1] * _number_of_cells[2];
    for (int_fast32_t i = 0; i < number_of_cells; ++i) {
      _hydro_variables[i].write_restart_file(restart_writer);
      _ionization_variables[i].write_restart_file(restart_writer);
    }
  }

  /**
   * @brief Restart constructor.
   *
   * @param restart_reader Restart file to read from.
   */
  inline MagnetoHydroDensitySubGrid(RestartReader &restart_reader)
      : DensitySubGridMHD(restart_reader) {

    _cell_volume = restart_reader.read< double >();
    _inverse_cell_volume = 1. / _cell_volume;
    _cell_areas[0] = restart_reader.read< double >();
    _cell_areas[1] = restart_reader.read< double >();
    _cell_areas[2] = restart_reader.read< double >();
    const int_fast32_t number_of_cells =
        _number_of_cells[0] * _number_of_cells[1] * _number_of_cells[2];
    _hydro_variables = new MagnetoHydroVariables[number_of_cells];
    for (int_fast32_t i = 0; i < number_of_cells; ++i) {
      _hydro_variables[i] = MagnetoHydroVariables(restart_reader);
      _ionization_variables[i] = IonizationVariables(restart_reader);
    }
    _primitive_variable_limiters = new double[18 * number_of_cells]; // mgb edit 21.04.2026 from 10 to 18
    for (int_fast32_t i = 0; i < 9 * number_of_cells; ++i) { // mgb edit 21.04.2026 from 5 to 9
      _primitive_variable_limiters[2 * i] = DBL_MAX;
      _primitive_variable_limiters[2 * i + 1] = -DBL_MAX;
    }
  }
};

#endif // MAGNETOHYDRODENSITYSUBGRID_HPP
