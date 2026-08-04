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
 * @file HydroDensitySubGrid.hpp
 *
 * @brief Extension of DensitySubGrid that adds hydro variables.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef MAGNETOHYDRODENSITYSUBGRID_HPP
#define MAGNETOHYDRODENSITYSUBGRID_HPP

#include "DensitySubGridMHD.hpp"
#include "DensityValues.hpp"
#include "MagnetoHydro.hpp"
#include "MagnetoHydroVariables.hpp"

/**
 * @brief Extension of DensitySubGrid that adds hydro variables.
 */
class MagnetoHydroDensitySubGrid : public DensitySubGridMHD {
private:
  /*! @brief Volume of a single cell (in m^3). */
  double _cell_volume;

  /*! @brief Inverse volume of a single cell (in m^-3). */
  double _inverse_cell_volume;

  /*! @brief Surface areas of a single cell (in m^2). */
  double _cell_areas[3];

  /*! @brief Hydrodynamical variables. */
  MagnetoHydroVariables *_hydro_variables;

  /*! @brief Gradient limiters for the primitive hydrodynamical variables. */
  double *_primitive_variable_limiters;

  /*! @brief Indices of the hydro tasks associated with this subgrid. */
  size_t _hydro_tasks[18]; //18

 // size_t _padding[100];

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
    const int_fast32_t tot_ncell = (size_t)ncell[0] * ncell[1] * ncell[2]; //_number_of_cells[3] * ncell[0];
    _hydro_variables = new MagnetoHydroVariables[tot_ncell];
    _primitive_variable_limiters = new double[tot_ncell * 20];

    for (int_fast32_t i = 0; i < 10 * tot_ncell; ++i) {
      _primitive_variable_limiters[2 * i] = DBL_MAX;
      _primitive_variable_limiters[2 * i + 1] = -DBL_MAX;
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

    const int_fast32_t tot_ncell = (size_t)original._number_of_cells[0] * original._number_of_cells[1] * original._number_of_cells[2];
    _hydro_variables = new MagnetoHydroVariables[tot_ncell];
    _primitive_variable_limiters = new double[tot_ncell * 20];

    // copy data arrays
    for (int_fast32_t i = 0; i < tot_ncell; ++i) {
      _hydro_variables[i].copy_all(original._hydro_variables[i]);
    }

    for (int_fast32_t i = 0; i < 20 * tot_ncell; ++i) {
      _primitive_variable_limiters[i] =
          original._primitive_variable_limiters[i];
    }
  
    
  }

  /**
   * @brief Destructor.
   */
  virtual ~MagnetoHydroDensitySubGrid() {
    // deallocate data arrays
    delete[] _hydro_variables;
    delete[] _primitive_variable_limiters;
  }

  /**
   * @brief Initialize the conserved variables for the grid.
   *
   * @param hydro Hydro instance to use.
   * @param do_primitives Initialize the primitive variables based on the
   * ionization variables?
   * @return Minimum initial time step for the cells in the grid.
   */
  inline double initialize_hydrodynamic_variables(const MagnetoHydro &hydro,
                                                  const bool do_primitives) {

    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    double timestep = DBL_MAX;
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      if (do_primitives) {
        hydro.ionization_to_hydro(_ionization_variables[i],
                                  _hydro_variables[i]);
      }
    //  if (hydro.get_hydro_variables().get_conserved_total_energy() < hydro.get_hydro_variables().get_conserved_internal_energy()){
      //  std::cout<<"Within Initialize hydrodynamic variables" << std::endl;
    //  }
      hydro.set_conserved_variables(_hydro_variables[i], _cell_volume);
   //   if (hydro.get_hydro_variables().get_conserved_total_energy() < hydro.get_hydro_variables().get_conserved_internal_energy()){
   //     std::cout<<"Warning! Etot < Eint after update in initialize hydro variables! Etot: " << hydro.get_hydro_variables().get_conserved_total_energy() << ", Eint: ", hydro.get_hydro_variables().get_conserved_internal_energy() << std::endl;
   //   }
      
      timestep = std::min(timestep, hydro.get_timestep(_hydro_variables[i],
                                                       _ionization_variables[i],
                                                       _cell_volume));
    //  std::cout << "timestep calculated in MagnetoHydroDensitySubGrid Line 151, dt = " << timestep << std::endl;
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
      if ((std::isnan(_hydro_variables[i].get_conserved_total_energy())) || (_hydro_variables[i].get_conserved_total_energy() == 0)) {
        std::cerr << "Total energy is NaN or 0. before conserved variable update for cell in MagnetoHydroDensitySubGrid" << i << ", mass: " << _hydro_variables[i].get_conserved_mass() << ", momentum: (" << _hydro_variables[i].get_conserved_momentum().x() << ", " << _hydro_variables[i].get_conserved_momentum().y() << ", " << _hydro_variables[i].get_conserved_momentum().z() << ")" << ", Etot: " << _hydro_variables[i].get_conserved_total_energy() << std::endl;
      }
      const CoordinateVector<> a =
          _hydro_variables[i].get_gravitational_acceleration();
      const CoordinateVector<> p = _hydro_variables[i].get_conserved_momentum();
      const CoordinateVector<> v = _hydro_variables[i].get_primitives_velocity();
      if (std::isnan(p.x()) || std::isnan(p.y()) || std::isnan(p.z())) {
        std::cerr << "Warning: Momentum is NaN before conserved variable update for cell in MagnetoHydroDensitySubGrid" << i << ", mass: " << _hydro_variables[i].get_conserved_mass() << ", velocity: (" << v.x() << ", " << v.y() << ", " << v.z() << ")" << std::endl;
      }
      const double mdt = _hydro_variables[i].get_conserved_mass() * timestep;
      if (std::isnan(mdt)) {
        std::cerr << "Warning: Mass is NaN before conserved variable update for cell in MagnetoHydroDensitySubGrid" << i << ", mass: " << _hydro_variables[i].get_conserved_mass() << std::endl;
      }
      
      double current_Etot = _hydro_variables[i].get_conserved_total_energy();
      double grav_work    = timestep * CoordinateVector<>::dot_product(p, a);
      double energy_term  = _hydro_variables[i].get_energy_term();
      double flux_delta   = _hydro_variables[i].delta_conserved(4) * timestep;
      double current_Eint = _hydro_variables[i].get_conserved_internal_energy();

      double final_Etot = current_Etot + grav_work + energy_term + flux_delta;
      if (final_Etot <= 0.0 || std::isnan(final_Etot)) {
          std::cout << "[FATAL ENERGY LOSS] Cell Index: " << i << "\n"
                    << "  Initial Etot:  " << current_Etot << "\n"
                    << " Initial Eint: " << current_Eint << "\n"
                    << "  Grav Work:     " << grav_work << "\n"
                    << "  Energy Term:   " << energy_term << "\n"
                    << "  Flux Delta:    " << flux_delta << "\n"
                    << "  Mass:          " << _hydro_variables[i].get_conserved_mass() << "\n"
                    << "  Result:        " << final_Etot << std::endl;
      }

      _hydro_variables[i].conserved(1) += mdt * a.x();
      _hydro_variables[i].conserved(2) += mdt * a.y();
      _hydro_variables[i].conserved(3) += mdt * a.z();
      _hydro_variables[i].conserved(4) +=
          timestep * CoordinateVector<>::dot_product(p, a);
      //_hydro_variables[i].conserved(9) +=
          //timestep * CoordinateVector<>::dot_product(p, a); // mgb edit 21.04.2026 - add gravitational energy change to conserved energy variable for better tracking of energy conservation in the presence of gravity
      _hydro_variables[i].conserved(4) += _hydro_variables[i].get_energy_term();
      _hydro_variables[i].conserved(9) += _hydro_variables[i].get_energy_term();
      _hydro_variables[i].set_energy_term(0.);
      for (int_fast8_t j = 0; j < 10; ++j) {
        _hydro_variables[i].conserved(j) +=
            _hydro_variables[i].delta_conserved(j) * timestep;
            
            

        // reset hydro variables
        _hydro_variables[i].delta_conserved(j) = 0;
        _hydro_variables[i].primitive_gradients(j) = CoordinateVector<>(0.);
        _primitive_variable_limiters[20 * i + 2 * j] = DBL_MAX;
        _primitive_variable_limiters[20 * i + 2 * j + 1] = -DBL_MAX;
      }

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
    "about to set kinetic greater than total energy....");

    if ((std::isnan(_hydro_variables[i].get_conserved_total_energy())) || (_hydro_variables[i].get_conserved_total_energy() == 0)) {
        std::cerr << "Total energy is NaN or 0. after conserved variable update for cell in MagnetoHydroDensitySubGrid" << i << ", mass: " << _hydro_variables[i].get_conserved_mass() << ", momentum: (" << _hydro_variables[i].get_conserved_momentum().x() << ", " << _hydro_variables[i].get_conserved_momentum().y() << ", " << _hydro_variables[i].get_conserved_momentum().z() << ")" << ", Etot: " << _hydro_variables[i].get_conserved_total_energy() << std::endl;
        
      }
    }
  }

  /**
   * @brief Update the primitive variables for the grid.
   *
   * @param hydro Hydro instance to use.
   */
  inline void update_primitive_variables(const MagnetoHydro &hydro) {

    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      if (std::isnan(_hydro_variables[i].get_conserved_total_energy()) || _hydro_variables[i].get_conserved_total_energy() == 0.) {
        std::cerr << "Total energy is NaN or zero before primitive variable update for cell in MagnetoHydroDensitySubGrid" << i << ", mass: " << _hydro_variables[i].get_conserved_mass() << ", momentum: (" << _hydro_variables[i].get_conserved_momentum().x() << ", " << _hydro_variables[i].get_conserved_momentum().y() << ", " << _hydro_variables[i].get_conserved_momentum().z() << ")" << std::endl;
      }
      hydro.set_primitive_variables(
          _hydro_variables[i], _ionization_variables[i], _inverse_cell_volume);
      
    }
  }

  /**
   * @brief Update the ionization variables for all cells in the subgrid using
   * their hydrodynamical variables.
   *
   * @param hydro Hydro instance to use.
   * @param maximum_neutral_fraction Maximum neutral fraction for hydrogen.
   */
  inline void
  update_ionization_variables(const MagnetoHydro &hydro,
                              const double maximum_neutral_fraction) {

    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      hydro.hydro_to_ionization(_hydro_variables[i], _ionization_variables[i]);
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
   * @param hydro Hydro instance to use.
   * @param timestep Integration time step (in s).
   */
  inline void add_ionization_energy(const MagnetoHydro &hydro, const double timestep) {

    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      hydro.add_ionization_energy(_ionization_variables[i], _hydro_variables[i],
                                  _inverse_cell_volume, timestep);
    }
  }

  /**
   * @brief Half time step prediction for the primitive variables.
   *
   * @param hydro Hydro instance to use.
   * @param timestep Half system time step (in s).
   */
  inline void predict_primitive_variables(const MagnetoHydro &hydro,
                                          const double timestep) {
    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      hydro.predict_primitive_variables(_hydro_variables[i], timestep);
    }
  }

  /**
   * @brief Apply the slope limiter to all primitive variable gradients.
   *
   * @param hydro Hydro instance to use.
   */
  inline void apply_slope_limiter(const MagnetoHydro &hydro) {
    const int_fast32_t tot_num_cells =
        _number_of_cells[0] * _number_of_cells[3];
    for (int_fast32_t i = 0; i < tot_num_cells; ++i) {
      hydro.apply_slope_limiter(_hydro_variables[i],
                                &_primitive_variable_limiters[20 * i],
                                _cell_size);
    }
  }

  /**
   * @brief Compute the hydrodynamical fluxes for all interfaces inside the
   * subgrid.
   *
   * @param hydro Hydro instance to use.
   * @param dt Current system time step (in s).
   */
  inline void inner_flux_sweep(const MagnetoHydro &hydro, const double dt) {

    // we do three separate sweeps: one for every coordinate direction
    for (int_fast32_t ix = 0; ix < _number_of_cells[0] - 1; ++ix) {
      for (int_fast32_t iy = 0; iy < _number_of_cells[1]; ++iy) {
        for (int_fast32_t iz = 0; iz < _number_of_cells[2]; ++iz) {
          const int_fast32_t index000 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          const int_fast32_t index100 =
              (ix + 1) * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          // x direction
          hydro.do_flux_calculation(0, _hydro_variables[index000],
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
          hydro.do_flux_calculation(1, _hydro_variables[index000],
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
          hydro.do_flux_calculation(2, _hydro_variables[index000],
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
   * @param hydro Hydro instance to use.
   * @param neighbour Neighbouring DensitySubGrid.
   * @param dt Current system time step (in s).
   */
  inline void outer_flux_sweep(const int_fast32_t direction, const MagnetoHydro &hydro,
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
      cmac_error("Unknown hydro neighbour: %" PRIiFAST32, direction);
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
        hydro.do_flux_calculation(i, left_grid->_hydro_variables[index_left],
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
   * @param hydro Hydro instance to use.
   * @param boundary HydroBoundary that sets the right state primitive
   * variables.
   * @param dt Current system time step (in s).
   */
  inline void outer_ghost_flux_sweep(const int_fast32_t direction,
                                     const MagnetoHydro &hydro,
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
      cmac_error("Unknown hydro neighbour: %" PRIiFAST32, direction);
      break;
    }
    if (sizeof(MagnetoHydro) < 10 * sizeof(double)) {
        cmac_error("CRITICAL: MagnetoHydro size mismatch! Expected at least 100 bytes, got %lu", sizeof(MagnetoHydro));
    }

    // using the index computation below is (much) faster than setting the
    // increment correctly and summing the indices manually
    for (int_fast32_t ic = 0; ic < column_length; ++ic) {
      for (int_fast32_t ir = 0; ir < row_length; ++ir) {
        const int_fast32_t index_left =
            start_index_left + ic * column_increment + ir * row_increment;
        hydro.do_ghost_flux_calculation(
            i, get_cell_midpoint(index_left) + offset,
            _hydro_variables[index_left], boundary, dx, A, dt);
      }
    }
  }

  /**
   * @brief Compute the hydrodynamical gradients for all interfaces inside the
   * subgrid.
   *
   * @param hydro Hydro instance to use.
   */
  inline void inner_gradient_sweep(const MagnetoHydro &hydro) {

    // we do three separate sweeps: one for every coordinate direction
    for (int_fast32_t ix = 0; ix < _number_of_cells[0] - 1; ++ix) {
      for (int_fast32_t iy = 0; iy < _number_of_cells[1]; ++iy) {
        for (int_fast32_t iz = 0; iz < _number_of_cells[2]; ++iz) {
          const int_fast32_t index000 =
              ix * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          const int_fast32_t index100 =
              (ix + 1) * _number_of_cells[3] + iy * _number_of_cells[2] + iz;
          // x direction
          hydro.do_gradient_calculation(
              0, _hydro_variables[index000], _hydro_variables[index100],
              _inv_cell_size[0], &_primitive_variable_limiters[20 * index000],
              &_primitive_variable_limiters[20 * index100]);
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
          hydro.do_gradient_calculation(
              1, _hydro_variables[index000], _hydro_variables[index010],
              _inv_cell_size[1], &_primitive_variable_limiters[20 * index000],
              &_primitive_variable_limiters[20 * index010]);
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
          hydro.do_gradient_calculation(
              2, _hydro_variables[index000], _hydro_variables[index001],
              _inv_cell_size[2], &_primitive_variable_limiters[20 * index000],
              &_primitive_variable_limiters[20 * index001]);
        }
      }
    }
  }

  /**
   * @brief Compute the hydrodynamical gradients for all interfaces at the
   * boundary between this subgrid and the given neighbouring subgrid.
   *
   * @param direction TravelDirection of the neighbour.
   * @param hydro Hydro instance to use.
   * @param neighbour Neighbouring DensitySubGrid.
   */
  inline void outer_gradient_sweep(const int_fast32_t direction,
                                   const MagnetoHydro &hydro,
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
      cmac_error("Unknown hydro neighbour: %" PRIiFAST32, direction);
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
        hydro.do_gradient_calculation(
            i, left_grid->_hydro_variables[index_left],
            right_grid->_hydro_variables[index_right], dxinv,
            &left_grid->_primitive_variable_limiters[20 * index_left],
            &right_grid->_primitive_variable_limiters[20 * index_right]);
      }
    }
  }

  /**
   * @brief Compute the hydrodynamical gradients for all interfaces at the
   * boundary between this subgrid and the given box boundary with boundary
   * condition.
   *
   * @param direction TravelDirection of the neighbour.
   * @param hydro Hydro instance to use.
   * @param boundary HydroBoundary that sets the right state primitive
   * variables.
   */
  inline void outer_ghost_gradient_sweep(const int_fast32_t direction,
                                         const MagnetoHydro &hydro,
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
      cmac_error("Unknown hydro neighbour: %" PRIiFAST32, direction);
      break;
    }
    if (sizeof(MagnetoHydro) < 10 * sizeof(double)) {
        cmac_error("CRITICAL: MagnetoHydro size mismatch! Expected at least 100 bytes, got %lu", sizeof(MagnetoHydro));
    }

    // using the index computation below is (much) faster than setting the
    // increment correctly and summing the indices manually
    for (int_fast32_t ic = 0; ic < column_length; ++ic) {
      for (int_fast32_t ir = 0; ir < row_length; ++ir) {
        const int_fast32_t index_left =
            start_index_left + ic * column_increment + ir * row_increment;
        hydro.do_ghost_gradient_calculation(
            i, get_cell_midpoint(index_left) + offset,
            left_grid->_hydro_variables[index_left], boundary, dxinv,
            &left_grid->_primitive_variable_limiters[20 * index_left]); 
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
   * @brief Get the hydro task with the given index.
   *
   * @param i Index.
   * @return Task.
   */
  inline size_t get_hydro_task(const int_fast32_t i) const {
   // std::cout << "Inside get hydro task, itask = " << _hydro_tasks[i] << std::endl; // mgb comment 28.04.2026
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
    _hydro_variables[index].set_primitives_velocity(values.get_velocity());
    _hydro_variables[index].set_primitives_magnetic_field(values.get_magnetic_field());
    _hydro_variables[index].set_primitives_B_scalar(values.get_B_scalar());
    // std::cout << "Initializing hydro variables for cell " << std::endl;
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
     * @brief Get read only access to the hydro variables stored in this
     * cell.
     *
     * @return Read only access to the hydro variables.
     */
    inline const MagnetoHydroVariables &get_hydro_variables() const {
      return _subgrid->_hydro_variables[_index];
    }

    /**
     * @brief Get read/write access to the hydro variables stored in this
     * cell.
     *
     * @return Read/write access to the hydro variables.
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
    _primitive_variable_limiters = new double[20 * number_of_cells];
    for (int_fast32_t i = 0; i < 10 * number_of_cells; ++i) {
      _primitive_variable_limiters[2 * i] = DBL_MAX;
      _primitive_variable_limiters[2 * i + 1] = -DBL_MAX;
    }
  }
};

#endif // HYDRODENSITYSUBGRID_HPP
