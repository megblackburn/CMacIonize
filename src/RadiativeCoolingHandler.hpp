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
 * @file RadiativeCooling.hpp
 *
 * @brief Radiative cooling rates parent class.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */


 #ifndef RADIATIVECOOLINGHANDLER_HPP
 #define RADIATIVECOOLINGHANDLER_HPP


 class RadiativeCoolingHandler {
 public:

    virtual ~RadiativeCoolingHandler() {};

    virtual double get_cooling_rate(double temperature) const = 0;

    virtual bool is_metal_only() const = 0;

    virtual double get_minimum_temperature() const = 0;

    virtual double get_maximum_temperature() const = 0;


 };

#endif // RADIATIVECOOLINGHANDLER_HPP