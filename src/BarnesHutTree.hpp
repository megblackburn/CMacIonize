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
 * @file BarnesHutTree.hpp
 *
 * @brief Barns Hut Tree
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk) ??
 */

#ifndef BARNESHUTTREE_HPP
#define BARNESHUTTREE_HPP

#include "PhysicalConstants.hpp"

#include <cmath>
#include <cstdint>
#include <memory>

class BarnesHutTree {
public:
  BarnesHutTree(const CoordinateVector<double> &grid_dimensions,
                const CoordinateVector<double> box_anchor, const double theta)
      : _grid_dimensions(grid_dimensions), _box_anchor(box_anchor),
        _theta(theta),
        _newton_G(PhysicalConstants::get_physical_constant(
            PHYSICALCONSTANT_NEWTON_CONSTANT)) {

    // One cubic root covers the whole domain, also for non-cubic boxes.
    const double root_dimension = _grid_dimensions.max();
    const CoordinateVector<double> root_midpoint =
        _box_anchor + 0.5 * _grid_dimensions;
    _root = std::unique_ptr< TreeNode >(
        new TreeNode(root_dimension, root_midpoint));
  }

  void insert(const HydroDensitySubGrid::hydroiterator &cell) {
    const double cell_mass = cell.get_hydro_variables().get_conserved_mass();
    if (!(cell_mass > 0.) || !std::isfinite(cell_mass)) {
      return;
    }
    if (_root) {
      _root->insert(cell);
    }
  }

  CoordinateVector<double>
  compute_gravity(HydroDensitySubGrid::hydroiterator &cell) const {
    CoordinateVector<double> acceleration(0.0, 0.0, 0.0);
    if (_root && _root->mass > 0.) {
      acceleration = _root->compute_gravity_recursive(cell, _theta, _newton_G);
    }
    return acceleration;
  }

  ~BarnesHutTree() = default;

  class TreeNode {
  public:
    double mass;
    CoordinateVector<double> center_of_mass;
    double dimensions;
    CoordinateVector<double> node_midpoint;
    std::unique_ptr<TreeNode> children[8];

    TreeNode(const double node_dimensions,
             const CoordinateVector<double> node_midpoint)
        : mass(0.), center_of_mass(CoordinateVector<double>()),
          dimensions(node_dimensions), node_midpoint(node_midpoint) {
      for (int i = 0; i < 8; ++i) {
        children[i] = nullptr;
      }
    }

    void insert(const HydroDensitySubGrid::hydroiterator &cell) {
      const double cell_mass = cell.get_hydro_variables().get_conserved_mass();
      if (!(cell_mass > 0.) || !std::isfinite(cell_mass)) {
        return;
      }

      const CoordinateVector<double> cell_position = cell.get_cell_midpoint();

      // Stable COM update that also handles previously empty nodes.
      if (mass > 0.) {
        const double total_mass = mass + cell_mass;
        center_of_mass =
            (center_of_mass * mass + cell_position * cell_mass) / total_mass;
        mass = total_mass;
      } else {
        center_of_mass = cell_position;
        mass = cell_mass;
      }

      const double resolution = std::cbrt(cell.get_volume());
      if (!(resolution > 0.) || dimensions <= 1.01 * resolution) {
        // Leaf-scale node: no further subdivision needed.
        return;
      }

      if (children[0] == nullptr) {
        subdivide();
      }

      children[child_index(cell_position)]->insert(cell);
    }

    int child_index(const CoordinateVector<double> &position) const {
      int index = 0;
      if (position.x() >= node_midpoint.x()) index |= 1;
      if (position.y() >= node_midpoint.y()) index |= 2;
      if (position.z() >= node_midpoint.z()) index |= 4;
      return index;
    }

    void subdivide() {
      for (int i = 0; i < 8; ++i) {
        const double child_dimensions = dimensions / 2.;
        CoordinateVector<double> child_midpoint =
            node_midpoint +
            CoordinateVector<double>((i & 1 ? 1. : -1.) * child_dimensions / 2.,
                                     (i & 2 ? 1. : -1.) * child_dimensions / 2.,
                                     (i & 4 ? 1. : -1.) * child_dimensions / 2.);
        children[i].reset(new TreeNode(child_dimensions, child_midpoint));
      }
    }

    CoordinateVector<double>
    compute_gravity_recursive(HydroDensitySubGrid::hydroiterator &cell,
                              const double theta,
                              const double newton_G) const {
      if (!(mass > 0.)) {
        return CoordinateVector<double>(0.0, 0.0, 0.0);
      }

      const CoordinateVector<double> position = cell.get_cell_midpoint();
      const CoordinateVector<double> r = center_of_mass - position;
      const double r2 = r.norm2();
      const bool is_leaf = (children[0] == nullptr);

      // Avoid singular self-force; tie the cut-off to local cell size.
      const double cell_size = std::cbrt(cell.get_volume());
      const double softening = 0.25 * cell_size;
      const double softening2 = softening * softening;

      if (is_leaf || (r2 > 0. && dimensions * dimensions <= theta * theta * r2)) {
        if (r2 <= softening2) {
          return CoordinateVector<double>(0.0, 0.0, 0.0);
        }

        const double inv_r = 1. / std::sqrt(r2);
        const double inv_r3 = inv_r * inv_r * inv_r;
        const double cell_mass = cell.get_hydro_variables().get_conserved_mass();

        // Store pair-wise potential energy with 1/2 factor to avoid double counting.
        const double potential = -0.5 * newton_G * mass * cell_mass * inv_r;
        cell.get_hydro_variables().set_gravitational_potential(
            cell.get_hydro_variables().get_gravitational_potential() +
            potential);

        return newton_G * mass * r * inv_r3;
      } else {
        CoordinateVector<double> acceleration(0., 0., 0.);
        for (const auto &child : children) {
          if (child) {
            acceleration +=
                child->compute_gravity_recursive(cell, theta, newton_G);
          }
        }
        return acceleration;
      }
    }

  };

private:
  std::unique_ptr< TreeNode > _root;
  CoordinateVector<double> _grid_dimensions;
  CoordinateVector<double> _box_anchor;
  double _theta;
  double _newton_G;
};

#endif
