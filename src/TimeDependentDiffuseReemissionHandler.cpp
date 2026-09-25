/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file TimeDependentDiffuseReemissionHandler.cpp
 *
 * @brief TimeDependentDiffuseReemissionHandler implementation.
 *
 * @author Meg Blackburn (mgb27@st-andrews.ac.uk)
 */
#include "TimeDependentDiffuseReemissionHandler.hpp"
#include "PhotonPacket.hpp"

/**
 * @brief Constructor.
 *
 * @param cross_sections Cross sections for photoionization.
 */
TimeDependentDiffuseReemissionHandler::TimeDependentDiffuseReemissionHandler(
    const CrossSections &cross_sections)
    : _HLyc_spectrum(cross_sections), _HeLyc_spectrum(cross_sections) {
    }


/**
 * @brief Reemit the given Photon.
 *
 * This routine randomly chooses if the photon is absorbed by hydrogen or
 * helium, and then determines a new frequency for the photon.
 *
 * @param photon Photon to reemit.
 * @param helium_abundance Abundance of helium.
 * @param ionization_variables IonizationVariables of the cell that contains the
 * current location of the Photon.
 * @param random_generator RandomGenerator to use.
 * @param type New type of the reemitted photon.
 * @return New frequency for the photon, or zero if the photon is absorbed.
 */
double TimeDependentDiffuseReemissionHandler::reemit(
    const PhotonPacket &photon, const double helium_abundance,
    IonizationVariables &ionization_variables,
    RandomGenerator &random_generator, PhotonType &type,
    PhotonPacketStatistics *statistics) const {

  double new_frequency = 0.;

  // Wood, Mathis & Ercolano (2004), section 3.3

  // determine whether the photon is absorbed by hydrogen or by helium
  const double nH0anuH0 = ionization_variables.get_ionic_fraction(ION_H_n) *
                          photon.get_photoionization_cross_section(ION_H_n)*
                          ionization_variables.get_number_density();
#ifdef HAS_HELIUM
  const double nHe0anuHe0 = ionization_variables.get_ionic_fraction(ION_He_n) *
                            helium_abundance *
                            photon.get_photoionization_cross_section(ION_He_n)*
                            ionization_variables.get_number_density();
#else
  const double nHe0anuHe0 = 0.;
#endif

    double dust_opacity = photon.get_dust_opacity();


    double ndustdust = dust_opacity * ionization_variables.get_dust_density();

  

    const double pDustabs = ndustdust / (nH0anuH0 + nHe0anuHe0 + ndustdust);

   


    double x = random_generator.get_uniform_random_double();

    if (x <= pDustabs) {
      // interaction wuth dust, lets do stuff

      x = random_generator.get_uniform_random_double();
      if (x <= ionization_variables.get_reemission_probability(
                   REEMISSIONPROBABILITY_DUST_ALBEDO)) {

        // keep same frequency as incoming photon
        new_frequency = photon.get_energy();
        type = PHOTONTYPE_SCATTERED;


      } else {

        // photon absorbed
        type = PHOTONTYPE_ABSORBED;
        //num_abs_dust.pre_increment();
        if (statistics != nullptr){
          if (photon.get_source_index() == statistics->get_tracked_source()){
            ionization_variables.increment_counter(true);
          }
          statistics->absorb_photon_dust();
        }

      }



    } else {


    type = PHOTONTYPE_ABSORBED;

    if (statistics != nullptr) {
      if (photon.get_source_index() == statistics->get_tracked_source()){
        ionization_variables.increment_counter(false);
      } 
    if (ionization_variables.get_number_density() > dens_thresh){
        statistics->absorb_photon(true);
      } else {
        statistics->absorb_photon(false);
      }
    }
  
}

  return new_frequency;

}

