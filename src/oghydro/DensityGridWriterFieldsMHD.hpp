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
 * @file DensityGridWriterFields.hpp
 *
 * @brief Variable field selection for DensityGrid output.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef DENSITYGRIDWRITERFIELDSMHD_HPP
#define DENSITYGRIDWRITERFIELDSMHD_HPP

#include "DensityGridMHD.hpp"
#include "Error.hpp"
#include "InternalMagnetoHydroUnits.hpp"
#include "ParameterFile.hpp"

#include <cinttypes>
#include <string>
#include <iostream>

/**
 * @brief Convenient indices for all supported output fields.
 */
enum DensityGridFieldMHD {
  DENSITYGRIDFIELDMHD_COORDINATES = 0,
  DENSITYGRIDFIELDMHD_NUMBER_DENSITY,
  DENSITYGRIDFIELDMHD_TEMPERATURE,
  DENSITYGRIDFIELDMHD_NEUTRAL_FRACTION,
#ifdef DO_OUTPUT_COOLING
  DENSITYGRIDFIELDMHD_COOLING,
#endif
#ifdef DO_OUTPUT_PHOTOIONIZATION_RATES
  DENSITYGRIDFIELDMHD_PHOTOIONIZATION_RATE,
#endif
#ifdef DO_OUTPUT_HEATING
  DENSITYGRIDFIELDMHD_HEATING_RATE,
#endif
  DENSITYGRIDFIELDMHD_COSMIC_RAY_FACTOR,
  DENSITYGRIDFIELDMHD_DENSITY,
  DENSITYGRIDFIELDMHD_VELOCITIES,
  DENSITYGRIDFIELDMHD_PRESSURE,
  DENSITYGRIDFIELDMHD_MASS,
  DENSITYGRIDFIELDMHD_MOMENTUM,
  DENSITYGRIDFIELDMHD_TOTAL_ENERGY,
  DENSITYGRIDFIELDMHD_ACCELERATION,
  DENSITYGRIDFIELDMHD_PHOTONCOUNTER,
  DENSITYGRIDFIELDMHD_MAGNETIC_FIELD,
  DENSITYGRIDFIELDMHD_B_SCALAR,
  DENSITYGRIDFIELDMHD_NUMBER
};

/**
 * @brief Convenient indices for all supported output field types.
 */
enum DensityGridFieldTypeMHD {
  DENSITYGRIDFIELDTYPEMHD_SCALAR_DOUBLE = 0,
  DENSITYGRIDFIELDTYPEMHD_VECTOR_DOUBLE,
  DENSITYGRIDFIELDTYPEMHD_NUMBER
};

/**
 * @brief Functionality to handle DensityGridFields and DensityGridFieldTypes.
 */
class DensityGridWriterFieldsMHD {
  /// STATIC FUNCTIONS
public:
  /**
   * @brief Get the DensityGridFieldType corresponding to the given
   * DensityGridField.
   *
   * @param field_name DensityGridField.
   * @return Corresponding DensityGridFieldType.
   */
  inline static int_fast32_t get_type(const int_fast32_t field_name) {
    switch (field_name) {
    case DENSITYGRIDFIELDMHD_COORDINATES:
      return DENSITYGRIDFIELDTYPEMHD_VECTOR_DOUBLE;
    case DENSITYGRIDFIELDMHD_NUMBER_DENSITY:
      return DENSITYGRIDFIELDTYPEMHD_SCALAR_DOUBLE;
    case DENSITYGRIDFIELDMHD_TEMPERATURE:
      return DENSITYGRIDFIELDTYPEMHD_SCALAR_DOUBLE;
    case DENSITYGRIDFIELDMHD_NEUTRAL_FRACTION:
      return DENSITYGRIDFIELDTYPEMHD_SCALAR_DOUBLE;
#ifdef DO_OUTPUT_COOLING
    case DENSITYGRIDFIELDMHD_COOLING:
      return DENSITYGRIDFIELDTYPEMHD_SCALAR_DOUBLE;
#endif
#ifdef DO_OUTPUT_PHOTOIONIZATION_RATES
    case DENSITYGRIDFIELDMHD_PHOTOIONIZATION_RATE:
      return DENSITYGRIDFIELDTYPEMHD_SCALAR_DOUBLE;
#endif
#ifdef DO_OUTPUT_HEATING
    case DENSITYGRIDFIELDMHD_HEATING_RATE:
      return DENSITYGRIDFIELDTYPEMHD_SCALAR_DOUBLE;
#endif
    case DENSITYGRIDFIELDMHD_COSMIC_RAY_FACTOR:
      return DENSITYGRIDFIELDTYPEMHD_SCALAR_DOUBLE;
    case DENSITYGRIDFIELDMHD_DENSITY:
      return DENSITYGRIDFIELDTYPEMHD_SCALAR_DOUBLE;
    case DENSITYGRIDFIELDMHD_VELOCITIES:
      return DENSITYGRIDFIELDTYPEMHD_VECTOR_DOUBLE;
    case DENSITYGRIDFIELDMHD_PRESSURE:
      return DENSITYGRIDFIELDTYPEMHD_SCALAR_DOUBLE;
    case DENSITYGRIDFIELDMHD_MASS:
      return DENSITYGRIDFIELDTYPEMHD_SCALAR_DOUBLE;
    case DENSITYGRIDFIELDMHD_MOMENTUM:
      return DENSITYGRIDFIELDTYPEMHD_VECTOR_DOUBLE;
    case DENSITYGRIDFIELDMHD_TOTAL_ENERGY:
      return DENSITYGRIDFIELDTYPEMHD_SCALAR_DOUBLE;
    case DENSITYGRIDFIELDMHD_ACCELERATION:
      return DENSITYGRIDFIELDTYPEMHD_VECTOR_DOUBLE;
    case DENSITYGRIDFIELDMHD_PHOTONCOUNTER:
      return DENSITYGRIDFIELDTYPEMHD_SCALAR_DOUBLE;
    case DENSITYGRIDFIELDMHD_MAGNETIC_FIELD:
      return DENSITYGRIDFIELDTYPEMHD_VECTOR_DOUBLE;
    case DENSITYGRIDFIELDMHD_B_SCALAR:
      return DENSITYGRIDFIELDTYPEMHD_SCALAR_DOUBLE;
    default:
      cmac_error("Unknown DensityGridField: %" PRIiFAST32, field_name);
      return DENSITYGRIDFIELDTYPEMHD_NUMBER;
    }
  }

  /**
   * @brief Get the name corresponding to the given DensityGridField.
   *
   * @param field_name DensityGridField.
   * @return Corresponding human readable name.
   */
  inline static std::string get_name(const int_fast32_t field_name) {
    switch (field_name) {
    case DENSITYGRIDFIELDMHD_COORDINATES:
      return "Coordinates";
    case DENSITYGRIDFIELDMHD_NUMBER_DENSITY:
      return "NumberDensity";
    case DENSITYGRIDFIELDMHD_TEMPERATURE:
      return "Temperature";
    case DENSITYGRIDFIELDMHD_NEUTRAL_FRACTION:
      return "NeutralFraction";
#ifdef DO_OUTPUT_COOLING
    case DENSITYGRIDFIELDMHD_COOLING:
      return "Cooling";
#endif
#ifdef DO_OUTPUT_PHOTOIONIZATION_RATES
    case DENSITYGRIDFIELDMHD_PHOTOIONIZATION_RATE:
      return "PhotoionizationRate";
#endif
#ifdef DO_OUTPUT_HEATING
    case DENSITYGRIDFIELDMHD_HEATING_RATE:
      return "HeatingRate";
#endif
    case DENSITYGRIDFIELDMHD_COSMIC_RAY_FACTOR:
      return "CosmicRayFactor";
    case DENSITYGRIDFIELDMHD_DENSITY:
      return "Density";
    case DENSITYGRIDFIELDMHD_VELOCITIES:
      return "Velocities";
    case DENSITYGRIDFIELDMHD_PRESSURE:
      return "Pressure";
    case DENSITYGRIDFIELDMHD_MASS:
      return "Mass";
    case DENSITYGRIDFIELDMHD_MOMENTUM:
      return "Momentum";
    case DENSITYGRIDFIELDMHD_TOTAL_ENERGY:
      return "TotalEnergy";
    case DENSITYGRIDFIELDMHD_ACCELERATION:
      return "Acceleration";
    case DENSITYGRIDFIELDMHD_PHOTONCOUNTER:
      return "PhotonCounter";
    case DENSITYGRIDFIELDMHD_MAGNETIC_FIELD:
      return "MagneticField";
    case DENSITYGRIDFIELDMHD_B_SCALAR:
      return "B_Scalar";
    default:
      cmac_error("Unknown DensityGridField: %" PRIiFAST32, field_name);
      return "";
    }
  }

  /**
   * @brief Default output flag for the given field.
   *
   * @param field_name DensityGridField.
   * @param hydro Flag specifying if hydro is active or not.
   * @return Whether or not the field is present by default.
   */
  inline static uint_fast32_t default_flag(const int_fast32_t field_name,
                                           const bool hydro) {
    switch (field_name) {
    case DENSITYGRIDFIELDMHD_COORDINATES:
      return true;
    case DENSITYGRIDFIELDMHD_NUMBER_DENSITY:
      // only output if hydro is not active
      return !hydro;
    case DENSITYGRIDFIELDMHD_TEMPERATURE:
      return hydro;
    case DENSITYGRIDFIELDMHD_NEUTRAL_FRACTION:
      // only output hydrogen neutral fraction
      return true;
#ifdef DO_OUTPUT_COOLING
    case DENSITYGRIDFIELDMHD_COOLING:
      return false;
#endif
#ifdef DO_OUTPUT_PHOTOIONIZATION_RATES
    case DENSITYGRIDFIELDMHD_PHOTOIONIZATION_RATE:
      return false;
#endif
#ifdef DO_OUTPUT_HEATING
    case DENSITYGRIDFIELDMHD_HEATING_RATE:
      return false;
#endif
    case DENSITYGRIDFIELDMHD_COSMIC_RAY_FACTOR:
      return false;
    case DENSITYGRIDFIELDMHD_DENSITY:
      return hydro;
    case DENSITYGRIDFIELDMHD_VELOCITIES:
      return hydro;
    case DENSITYGRIDFIELDMHD_PRESSURE:
      return hydro;
    case DENSITYGRIDFIELDMHD_MASS:
      return false;
    case DENSITYGRIDFIELDMHD_MOMENTUM:
      return false;
    case DENSITYGRIDFIELDMHD_TOTAL_ENERGY:
      return false;
    case DENSITYGRIDFIELDMHD_ACCELERATION:
      return false;
    case DENSITYGRIDFIELDMHD_PHOTONCOUNTER:
      return false;
    case DENSITYGRIDFIELDMHD_MAGNETIC_FIELD:
      return false;
    case DENSITYGRIDFIELDMHD_B_SCALAR:
      return false;
    default:
      cmac_error("Unknown DensityGridField: %" PRIiFAST32, field_name);
      return false;
    }
  }

  /**
   * @brief Is the given DensityGridField an ion property?
   *
   * @param field_name DensityGridField.
   * @return True if the given field has a value for each ion.
   */
  inline static bool is_ion_property(const int_fast32_t field_name) {
    switch (field_name) {
    case DENSITYGRIDFIELDMHD_COORDINATES:
      return false;
    case DENSITYGRIDFIELDMHD_NUMBER_DENSITY:
      return false;
    case DENSITYGRIDFIELDMHD_TEMPERATURE:
      return false;
    case DENSITYGRIDFIELDMHD_NEUTRAL_FRACTION:
      return true;
#ifdef DO_OUTPUT_COOLING
    case DENSITYGRIDFIELDMHD_COOLING:
      return true;
#endif
#ifdef DO_OUTPUT_PHOTOIONIZATION_RATES
    case DENSITYGRIDFIELDMHD_PHOTOIONIZATION_RATE:
      return true;
#endif
#ifdef DO_OUTPUT_HEATING
    case DENSITYGRIDFIELDMHD_HEATING_RATE:
      return false;
#endif
    case DENSITYGRIDFIELDMHD_COSMIC_RAY_FACTOR:
      return false;
    case DENSITYGRIDFIELDMHD_DENSITY:
      return false;
    case DENSITYGRIDFIELDMHD_VELOCITIES:
      return false;
    case DENSITYGRIDFIELDMHD_PRESSURE:
      return false;
    case DENSITYGRIDFIELDMHD_MASS:
      return false;
    case DENSITYGRIDFIELDMHD_MOMENTUM:
      return false;
    case DENSITYGRIDFIELDMHD_TOTAL_ENERGY:
      return false;
    case DENSITYGRIDFIELDMHD_ACCELERATION:
      return false;
    case DENSITYGRIDFIELDMHD_PHOTONCOUNTER:
      return false;
    case DENSITYGRIDFIELDMHD_MAGNETIC_FIELD:
      return false;
    case DENSITYGRIDFIELDMHD_B_SCALAR:
      return false;
    default:
      cmac_error("Unknown DensityGridField: %" PRIiFAST32, field_name);
      return false;
    }
  }

  /**
   * @brief Is the given DensityGridField a heating term property?
   *
   * @param field_name DensityGridField.
   * @return True if the given field has a value for each heating term.
   */
  inline static bool is_heating_property(const int_fast32_t field_name) {
    switch (field_name) {
    case DENSITYGRIDFIELDMHD_COORDINATES:
      return false;
    case DENSITYGRIDFIELDMHD_NUMBER_DENSITY:
      return false;
    case DENSITYGRIDFIELDMHD_TEMPERATURE:
      return false;
    case DENSITYGRIDFIELDMHD_NEUTRAL_FRACTION:
      return false;
#ifdef DO_OUTPUT_COOLING
    case DENSITYGRIDFIELDMHD_COOLING:
      return false;
#endif
#ifdef DO_OUTPUT_PHOTOIONIZATION_RATES
    case DENSITYGRIDFIELDMHD_PHOTOIONIZATION_RATE:
      return false;
#endif
#ifdef DO_OUTPUT_HEATING
    case DENSITYGRIDFIELDMHD_HEATING_RATE:
      return true;
#endif
    case DENSITYGRIDFIELDMHD_COSMIC_RAY_FACTOR:
      return false;
    case DENSITYGRIDFIELDMHD_DENSITY:
      return false;
    case DENSITYGRIDFIELDMHD_VELOCITIES:
      return false;
    case DENSITYGRIDFIELDMHD_PRESSURE:
      return false;
    case DENSITYGRIDFIELDMHD_MASS:
      return false;
    case DENSITYGRIDFIELDMHD_MOMENTUM:
      return false;
    case DENSITYGRIDFIELDMHD_TOTAL_ENERGY:
      return false;
    case DENSITYGRIDFIELDMHD_ACCELERATION:
      return false;
    case DENSITYGRIDFIELDMHD_PHOTONCOUNTER:
      return false;
    case DENSITYGRIDFIELDMHD_MAGNETIC_FIELD:
      return false;
    case DENSITYGRIDFIELDMHD_B_SCALAR:
      return false;
    default:
      cmac_error("Unknown DensityGridField: %" PRIiFAST32, field_name);
      return false;
    }
  }

  /**
   * @brief Is the given DensityGridField a hydro property?
   *
   * @param field_name DensityGridField.
   * @return True if the given field is a hydro only property.
   */
  inline static bool is_hydro_property(const int_fast32_t field_name) {
    switch (field_name) {
    case DENSITYGRIDFIELDMHD_COORDINATES:
      return false;
    case DENSITYGRIDFIELDMHD_NUMBER_DENSITY:
      return false;
    case DENSITYGRIDFIELDMHD_TEMPERATURE:
      return false;
    case DENSITYGRIDFIELDMHD_NEUTRAL_FRACTION:
      return false;
#ifdef DO_OUTPUT_COOLING
    case DENSITYGRIDFIELDMHD_COOLING:
      return false;
#endif
#ifdef DO_OUTPUT_PHOTOIONIZATION_RATES
    case DENSITYGRIDFIELDMHD_PHOTOIONIZATION_RATE:
      return false;
#endif
#ifdef DO_OUTPUT_HEATING
    case DENSITYGRIDFIELDMHD_HEATING_RATE:
      return false;
#endif
    case DENSITYGRIDFIELDMHD_COSMIC_RAY_FACTOR:
      return false;
    case DENSITYGRIDFIELDMHD_DENSITY:
      return true;
    case DENSITYGRIDFIELDMHD_VELOCITIES:
      return true;
    case DENSITYGRIDFIELDMHD_PRESSURE:
      return true;
    case DENSITYGRIDFIELDMHD_MASS:
      return true;
    case DENSITYGRIDFIELDMHD_MOMENTUM:
      return true;
    case DENSITYGRIDFIELDMHD_TOTAL_ENERGY:
      return true;
    case DENSITYGRIDFIELDMHD_ACCELERATION:
      return true;
    case DENSITYGRIDFIELDMHD_PHOTONCOUNTER:
      return false;
    case DENSITYGRIDFIELDMHD_MAGNETIC_FIELD:
      return false;
    case DENSITYGRIDFIELDMHD_B_SCALAR:
      return false;
    default:
      cmac_error("Unknown DensityGridField: %" PRIiFAST32, field_name);
      return false;
    }
  }

  /**
   * @brief Get the value of the double scalar field with the given name.
   *
   * @param field_name DensityGridField.
   * @param it DensityGrid::iterator to a cell.
   * @param hydro_units Internal unit system used for hydro.
   * @return Value of the double scalar field.
   */
  inline static double
  get_scalar_double_value(const int_fast32_t field_name,
                          const DensityGridMHD::iterator &it,
                          const InternalMagnetoHydroUnits *hydro_units = nullptr) {
    switch (field_name) {
    case DENSITYGRIDFIELDMHD_NUMBER_DENSITY:
      return it.get_ionization_variables().get_number_density();
    case DENSITYGRIDFIELDMHD_TEMPERATURE:
      return it.get_ionization_variables().get_temperature();
    case DENSITYGRIDFIELDMHD_COSMIC_RAY_FACTOR:
      return it.get_ionization_variables().get_cosmic_ray_factor();
    case DENSITYGRIDFIELDMHD_DENSITY: {
      if (hydro_units != nullptr) {
        return hydro_units->convert_to_SI_units< QUANTITY_DENSITY >(
            it.get_hydro_variables().get_primitives_density());
      } else {
        return it.get_hydro_variables().get_primitives_density();
      }
    }
    case DENSITYGRIDFIELDMHD_PRESSURE: {
      if (hydro_units != nullptr) {
        return hydro_units->convert_to_SI_units< QUANTITY_PRESSURE >(
            it.get_hydro_variables().get_primitives_pressure());
      } else {
        return it.get_hydro_variables().get_primitives_pressure();
      }
    }
    case DENSITYGRIDFIELDMHD_MASS: {
      if (hydro_units != nullptr) {
        return hydro_units->convert_to_SI_units< QUANTITY_MASS >(
            it.get_hydro_variables().get_conserved_mass());
      } else {
        return it.get_hydro_variables().get_conserved_mass();
      }
    }
    case DENSITYGRIDFIELDMHD_TOTAL_ENERGY: {
      if (hydro_units != nullptr) {
        return hydro_units->convert_to_SI_units< QUANTITY_ENERGY >(
            it.get_hydro_variables().get_conserved_total_energy());
      } else {
        return it.get_hydro_variables().get_conserved_total_energy();
      }
    }
    case DENSITYGRIDFIELDMHD_PHOTONCOUNTER: {
      return it.get_ionization_variables().get_counter(false);
    }
    case DENSITYGRIDFIELDMHD_B_SCALAR: {
      if (hydro_units != nullptr) {
        return hydro_units->convert_to_SI_units<QUANTITY_MAGNETIC_FIELD>(
          it.get_hydro_variables().get_primitives_B_scalar());
      } else {
        return it.get_hydro_variables().get_primitives_B_scalar();
      }
    }
    default:
      cmac_error("Not a scalar DensityGridField: %" PRIiFAST32, field_name);
      return 0.;
    }
  }

  /**
   * @brief Get the value of the double scalar field with the given name.
   *
   * @param field_name DensityGridField.
   * @param it DensitySubGrid::iterator to a cell.
   * @param hydro_units Internal unit system used for hydro.
   * @return Value of the double scalar field.
   */
  template < typename _cell_iterator_ >
  inline static double
  get_scalar_double_value(const int_fast32_t field_name,
                          const _cell_iterator_ &it,
                          const InternalMagnetoHydroUnits *hydro_units = nullptr) {
    switch (field_name) {
    case DENSITYGRIDFIELDMHD_NUMBER_DENSITY:
      return it.get_ionization_variables().get_number_density();
    case DENSITYGRIDFIELDMHD_TEMPERATURE:
      return it.get_ionization_variables().get_temperature();
    case DENSITYGRIDFIELDMHD_COSMIC_RAY_FACTOR:
      return it.get_ionization_variables().get_cosmic_ray_factor();
    case DENSITYGRIDFIELDMHD_DENSITY: {
      if (hydro_units != nullptr) {
        cmac_error("Not implemented yet!");
        return 0.;
      } else {
        return it.get_hydro_variables().get_primitives_density();
      }
    }
    case DENSITYGRIDFIELDMHD_PRESSURE: {
      if (hydro_units != nullptr) {
        cmac_error("Not implemented yet!");
        return 0.;
      } else {
        return it.get_hydro_variables().get_primitives_pressure();
      }
    }
    case DENSITYGRIDFIELDMHD_MASS: {
      if (hydro_units != nullptr) {
        cmac_error("Not implemented yet!");
        return 0.;
      } else {
        return it.get_hydro_variables().get_conserved_mass();
      }
    }
    case DENSITYGRIDFIELDMHD_TOTAL_ENERGY: {
      if (hydro_units != nullptr) {
        cmac_error("Not implemented yet!");
        return 0.;
      } else {
        return it.get_hydro_variables().get_conserved_total_energy();
      }
    }
    case DENSITYGRIDFIELDMHD_PHOTONCOUNTER: {
      return it.get_ionization_variables().get_counter(false);
    }
    case DENSITYGRIDFIELDMHD_B_SCALAR: {
      if (hydro_units != nullptr) {
        cmac_error("Not implemented yet!");
        return 0.;
      } else {
        return it.get_hydro_variables().get_primitives_B_scalar();
      }
    }
    default:
      cmac_error("Not a scalar DensityGridField: %" PRIiFAST32, field_name);
      return 0.;
    }
  }

  /**
   * @brief Get the value of the double vector field with the given name.
   *
   * @param field_name DensityGridField.
   * @param it DensityGrid::iterator to a cell.
   * @param box_anchor Anchor of the simulation box, used for vector corrections
   * (in m).
   * @param hydro_units Internal unit system used for hydro.
   * @return Value of the double vector field.
   */
  inline static CoordinateVector<>
  get_vector_double_value(const int_fast32_t field_name,
                          const DensityGridMHD::iterator &it,
                          const CoordinateVector<> box_anchor,
                          const InternalMagnetoHydroUnits *hydro_units = nullptr) {
    switch (field_name) {
    case DENSITYGRIDFIELDMHD_COORDINATES:
      return it.get_cell_midpoint() - box_anchor;
    case DENSITYGRIDFIELDMHD_VELOCITIES: {
      if (hydro_units != nullptr) {
        return hydro_units->convert_to_SI_units< QUANTITY_VELOCITY >(
            it.get_hydro_variables().get_primitives_velocity());
      } else {
        return it.get_hydro_variables().get_primitives_velocity();
      }
    }
    case DENSITYGRIDFIELDMHD_MOMENTUM: {
      if (hydro_units != nullptr) {
        return hydro_units->convert_to_SI_units< QUANTITY_MOMENTUM >(
            it.get_hydro_variables().get_conserved_momentum());
      } else {
        return it.get_hydro_variables().get_conserved_momentum();
      }
    }
    case DENSITYGRIDFIELDMHD_ACCELERATION: {
      if (hydro_units != nullptr) {
        return hydro_units->convert_to_SI_units< QUANTITY_ACCELERATION >(
            it.get_hydro_variables().get_gravitational_acceleration());
      } else {
        return it.get_hydro_variables().get_gravitational_acceleration();
      }
    }
    case DENSITYGRIDFIELDMHD_MAGNETIC_FIELD: {
      if (hydro_units != nullptr) {
        return hydro_units->convert_to_SI_units<QUANTITY_MAGNETIC_FIELD>(
          it.get_hydro_variables().get_conserved_magnetic_field());
      } else {
        return it.get_hydro_variables().get_conserved_magnetic_field();
      }
    }
    default:
      cmac_error("Not a vector DensityGridField: %" PRIiFAST32, field_name);
      return CoordinateVector<>();
    }
  }

  /**
   * @brief Get the value of the double vector field with the given name.
   *
   * @param field_name DensityGridField.
   * @param it DensitySubGrid::iterator to a cell.
   * @param box_anchor Anchor of the simulation box, used for vector corrections
   * (in m).
   * @param hydro_units Internal unit system used for hydro.
   * @return Value of the double vector field.
   */
  template < typename _cell_iterator_ >
  inline static CoordinateVector<>
  get_vector_double_value(const int_fast32_t field_name,
                          const _cell_iterator_ &it,
                          const CoordinateVector<> box_anchor,
                          const InternalMagnetoHydroUnits *hydro_units = nullptr) {
    switch (field_name) {
    case DENSITYGRIDFIELDMHD_COORDINATES:
      return it.get_cell_midpoint() - box_anchor;
    case DENSITYGRIDFIELDMHD_VELOCITIES: {
      if (hydro_units != nullptr) {
        cmac_error("Not implemented yet!");
        return 0.;
      } else {
        return it.get_hydro_variables().get_primitives_velocity();
      }
    }
    case DENSITYGRIDFIELDMHD_MOMENTUM: {
      if (hydro_units != nullptr) {
        cmac_error("Not implemented yet!");
        return 0.;
      } else {
        return it.get_hydro_variables().get_conserved_momentum();
      }
    }
    case DENSITYGRIDFIELDMHD_ACCELERATION: {
      if (hydro_units != nullptr) {
        cmac_error("Not implemented yet!");
        return 0.;
      } else {
        return it.get_hydro_variables().get_gravitational_acceleration();
      }
    }
    case DENSITYGRIDFIELDMHD_MAGNETIC_FIELD: {
      if (hydro_units != nullptr) {
        cmac_error("Not implemented yet!");
        return 0.;
      } else {
        return it.get_hydro_variables().get_primitives_magnetic_field();
      }
    }

    default:
      cmac_error("Not a vector DensityGridField: %" PRIiFAST32, field_name);
      return CoordinateVector<>();
    }
  }

  /**
   * @brief Get the value of the double scalar ion field with the given name.
   *
   * @param field_name DensityGridField.
   * @param ion_name IonName.
   * @param it DensityGrid::iterator to a cell.
   * @return Value of the double scalar ion field.
   */
  inline static double
  get_scalar_double_ion_value(const int_fast32_t field_name,
                              const int_fast32_t ion_name,
                              const DensityGridMHD::iterator &it) {
    switch (field_name) {
    case DENSITYGRIDFIELDMHD_NEUTRAL_FRACTION:
      return it.get_ionization_variables().get_ionic_fraction(ion_name);
#ifdef DO_OUTPUT_COOLING
    case DENSITYGRIDFIELDMHD_COOLING:
      return it.get_ionization_variables().get_cooling(ion_name);
#endif
#ifdef DO_OUTPUT_PHOTOIONIZATION_RATES
    case DENSITYGRIDFIELDMHD_PHOTOIONIZATION_RATE:
      return it.get_ionization_variables().get_mean_intensity(ion_name);
#endif
    default:
      cmac_error("Not a scalar ion DensityGridField: %" PRIiFAST32, field_name);
      return 0.;
    }
  }

  /**
   * @brief Get the value of the double scalar ion field with the given name.
   *
   * @param field_name DensityGridField.
   * @param ion_name IonName.
   * @param it DensitySubGrid::iterator to a cell.
   * @return Value of the double scalar ion field.
   */
  template < typename _cell_iterator_ >
  inline static double
  get_scalar_double_ion_value(const int_fast32_t field_name,
                              const int_fast32_t ion_name,
                              const _cell_iterator_ &it) {

    switch (field_name) {
    case DENSITYGRIDFIELDMHD_NEUTRAL_FRACTION:
      return it.get_ionization_variables().get_ionic_fraction(ion_name);
#ifdef DO_OUTPUT_COOLING
    case DENSITYGRIDFIELDMHD_COOLING:
      return it.get_ionization_variables().get_cooling(ion_name);
#endif
#ifdef DO_OUTPUT_PHOTOIONIZATION_RATES
    case DENSITYGRIDFIELDMHD_PHOTOIONIZATION_RATE:
      return it.get_ionization_variables().get_mean_intensity(ion_name);
#endif
    default:
      cmac_error("Not a scalar ion DensityGridField: %" PRIiFAST32, field_name);
      return 0.;
    }
  }

  /**
   * @brief Get the value of the double scalar heating property field with the
   * given name.
   *
   * @param field_name DensityGridField.
   * @param heating_property_name HeatingTermName.
   * @param it DensityGrid::iterator to a cell.
   * @return Value fo the double scalar heating property field.
   */
  inline static double
  get_scalar_double_heating_value(const int_fast32_t field_name,
                                  const int_fast32_t heating_property_name,
                                  const DensityGridMHD::iterator &it) {
    switch (field_name) {
#ifdef DO_OUTPUT_HEATING
    case DENSITYGRIDFIELDMHD_HEATING_RATE:
      return it.get_ionization_variables().get_heating(heating_property_name);
#endif
    default:
      cmac_error("Not a scalar heating property DensityGridField: %" PRIiFAST32,
                 field_name);
      return false;
    }
  }

  /**
   * @brief Get the value of the double scalar heating property field with the
   * given name.
   *
   * @param field_name DensityGridField.
   * @param heating_property_name HeatingTermName.
   * @param it DensitySubGrid::iterator to a cell.
   * @return Value fo the double scalar heating property field.
   */
  template < typename _cell_iterator_ >
  inline static double
  get_scalar_double_heating_value(const int_fast32_t field_name,
                                  const int_fast32_t heating_property_name,
                                  const _cell_iterator_ &it) {
    switch (field_name) {
#ifdef DO_OUTPUT_HEATING
    case DENSITYGRIDFIELDMHD_HEATING_RATE:
      return it.get_ionization_variables().get_heating(heating_property_name);
#endif
    default:
      cmac_error("Not a scalar heating property DensityGridField: %" PRIiFAST32,
                 field_name);
      return false;
    }
  }

private:
  /*! @brief Field flag for each DensityGridField. For fields with multiple
   *  variables this flag encodes the variables that should be written to the
   *  output. */
  uint_least64_t _field_flag[DENSITYGRIDFIELDMHD_NUMBER];

  /*! @brief Number of active fields of each type. */
  uint_least8_t _field_count[DENSITYGRIDFIELDTYPEMHD_NUMBER];

  /**
   * @brief Return the number of 1 bits in the given sequence.
   *
   * Note that some cpus have built in instructions for this. We do not
   * particularly care about speed, so we just use a naive loop.
   *
   * @param sequence 32-bit sequence.
   * @return Number of 1 bits.
   */
  inline static uint_fast8_t bit_count(const uint_fast32_t sequence) {
    uint_fast8_t count = 0;
    for (uint_fast8_t i = 0; i < 32; ++i) {
      count += (sequence >> i) & 1;
    }
    return count;
  }

public:
  /**
   * @brief Default constructor.
   *
   * @param hydro Output hydro variables?
   */
  inline DensityGridWriterFieldsMHD(const bool hydro) {

    for (int_fast32_t type = 0; type < DENSITYGRIDFIELDTYPEMHD_NUMBER; ++type) {
      _field_count[type] = 0;
    }

    for (int_fast32_t property = 0; property < DENSITYGRIDFIELDMHD_NUMBER;
         ++property) {
      _field_flag[property] = default_flag(property, hydro);
      _field_count[get_type(property)] += bit_count(_field_flag[property]);
    }
  }

  /**
   * @brief Test constructor.
   *
   * Set up a DensityGridWriterFields instance with hand selected output fields.
   *
   * @param flags Field flags.
   */
  inline DensityGridWriterFieldsMHD(
      const uint_fast32_t flags[DENSITYGRIDFIELDMHD_NUMBER]) {

    for (int_fast32_t type = 0; type < DENSITYGRIDFIELDTYPEMHD_NUMBER; ++type) {
      _field_count[type] = 0;
    }

    for (int_fast32_t property = 0; property < DENSITYGRIDFIELDMHD_NUMBER;
         ++property) {
      _field_flag[property] = flags[property];
      _field_count[get_type(property)] += bit_count(_field_flag[property]);
    }
  }

  /**
   * @brief Copy constructor.
   *
   * @param copy DensityGridWriterFields instance to copy from.
   */
  inline DensityGridWriterFieldsMHD(const DensityGridWriterFieldsMHD &copy) {

    for (int_fast32_t property = 0; property < DENSITYGRIDFIELDMHD_NUMBER;
         ++property) {
      _field_flag[property] = copy._field_flag[property];
    }

    for (int_fast32_t type = 0; type < DENSITYGRIDFIELDTYPEMHD_NUMBER; ++type) {
      _field_count[type] = copy._field_count[type];
    }
  }

  /**
   * @brief ParameterFile constructor.
   *
   * @param params ParameterFile to read from.
   * @param hydro Flag specifying whether or not hydro is active.
   */
  inline DensityGridWriterFieldsMHD(ParameterFile &params, const bool hydro) {

    for (int_fast32_t type = 0; type < DENSITYGRIDFIELDTYPEMHD_NUMBER; ++type) {
      _field_count[type] = 0;
    }

    for (int_fast32_t property = 0; property < DENSITYGRIDFIELDMHD_NUMBER;
         ++property) {
      if (hydro || !is_hydro_property(property)) {
        const std::string prop_name = get_name(property);
        if (is_ion_property(property)) {
          _field_flag[property] = 0;
          for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
            const std::string name = prop_name + get_ion_name(ion);
            _field_flag[property] +=
                params.get_value< uint_fast32_t >(
                    "DensityGridWriterFields:" + name,
                    default_flag(property, hydro) & (1u << ion))
                << ion;
          }
        } else if (is_heating_property(property)) {
          _field_flag[property] = 0;
          for (int_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
               ++heating) {
            const std::string name = prop_name + get_ion_name(heating);
            _field_flag[property] +=
                params.get_value< uint_fast32_t >(
                    "DensityGridWriterFields:" + name,
                    default_flag(property, hydro) & (1u << heating))
                << heating;
          }
        } else {
          // normal single variable property
          _field_flag[property] = params.get_value< uint_fast32_t >(
              "DensityGridWriterFields:" + prop_name,
              default_flag(property, hydro));
        }
        _field_count[get_type(property)] += bit_count(_field_flag[property]);
      } else {
        _field_flag[property] = false;
      }
    }
  }

  /**
   * @brief Check if the given field should be output.
   *
   * @param field_name DensityGridField.
   * @return True if the field should be output.
   */
  inline bool field_present(const int_fast32_t field_name) const {
    return _field_flag[field_name] > 0;
  }

  /**
   * @brief Check if the given ion for the given field should be output.
   *
   * @param field_name DensityGridField.
   * @param ion IonName.
   * @return True if the given ion for the given field should be output.
   */
  inline bool ion_present(const int_fast32_t field_name,
                          const int_fast32_t ion) const {
    return ((_field_flag[field_name] >> ion) & 1) != 0;
  }

  /**
   * @brief Check if the given heating term for the given field should be
   * output.
   *
   * @param field_name DensityGridField.
   * @param heating HeatingTermName.
   * @return True if the given ion for the given field should be output.
   */
  inline bool heatingterm_present(const int_fast32_t field_name,
                                  const int_fast32_t heating) const {
    return ((_field_flag[field_name] >> heating) & 1) != 0;
  }

  /**
   * @brief Get the number of fields of the given type that should be output.
   *
   * @param type DensityGridFieldType.
   * @return Number of fields of the given type that should be output.
   */
  inline uint_fast32_t get_field_count(const int_fast32_t type) const {
    return _field_count[type];
  }
};

#endif // DENSITYGRIDWRITERFIELDS_HPP
