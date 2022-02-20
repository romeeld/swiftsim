/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Tsang Keung Chan (chantsangkeung@gmail.com)
 *               2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_RT_PROPERTIES_SPHM1RT_H
#define SWIFT_RT_PROPERTIES_SPHM1RT_H

#include "rt_parameters.h"

/**
 * @file src/rt/SPHM1RT/rt_properties.h
 * @brief Main header file for the 'SPHM1RT' radiative transfer scheme
 * properties. SPHM1RT method described in Chan+21: 2102.08404
 */

/**
 * @brief Properties of the 'SPHM1RT' radiative transfer model
 */
struct rt_props {

  /* CFL condition */
  float CFL_condition;

  /* reduced speed of light in code unit */
  float cred;

  /*! initial opacity */
  float initialchi[RT_NGROUPS];

  /* Frequency bin edges for photon groups
   * Includes 0 as leftmost edge, doesn't include infinity as
   * rightmost bin edge*/
  float photon_groups[RT_NGROUPS];

  /* Are we running with hydro or star controlled injection?
   * This is added to avoid #ifdef macros as far as possible */
  int hydro_controlled_injection;

  /* Do we need to run a conversion after the zeroth
   * step, but before the first step? */
  int convert_stars_after_zeroth_step;
  int convert_parts_after_zeroth_step;

  /* Are we using constant stellar emission rates? */
  int use_const_emission_rates;

  /* Global constant stellar emission rates */
  double stellar_const_emission_rates[RT_NGROUPS];

  /* Optionally restrict maximal timestep for stars */
  float stars_max_timestep;

  /* Which stellar spectrum type to use? */
  int stellar_spectrum_type;
  /* If constant: get max frequency */
  double const_stellar_spectrum_max_frequency;
  /* If blackbody: get temperature */
  double stellar_spectrum_blackbody_T;

  /*! Fraction of the particle mass in given elements at the start of the run */
  float initial_metal_mass_fraction[rt_chemistry_element_count];

  /*! Fraction of the particle mass in *all* metals at the start of the run */
  float initial_metal_mass_fraction_total;

  /* switch for the on the spot approximation? */
  int onthespot;

  /*! The energy of an ionizing photon in cgs units */
  double ionizing_photon_energy_cgs[3];  

  /* switch for cooling (and photoheating) */
  /* however, there is still photo-ionization even if the switch is off */
  int coolingon;

  /* switch for fixing the photo-density */
  int fixphotondensity;

  /* the photo-density values if the photo density is fixed  */
  float Fgamma_fixed_cgs[RT_NGROUPS];

  /*! switch to use thermo-chemistry parameters from the parameter file */
  float useparams;

  /* the following only use when useparam = 1 */
  /*! The case A recombination coefficient for hydrogen (cgs) */
  double alphaA_cgs_H;

  /*! The case B recombination coefficient for hydrogen (cgs) */
  double alphaB_cgs_H;

  /*! The collisional ionization coefficient for hydrogen (cgs) */
  double beta_cgs_H;

  /*! The cross section of ionizing photons for hydrogen (cgs) */
  double sigma_cross_cgs_H[3];   

  /*** end of useparams ***/

  /*! tolerance of relative change to shift from explicit solver to CVODE */
  double explicitRelTolerance;

  /*! CVODE absolute tolerance */
  double absoluteTolerance;

  /*! CVODE relative tolerance */
  double relativeTolerance;



};


/**
 * @brief Return a string containing the name of a given #rt_cooling_species.
 */
__attribute__((always_inline)) INLINE static const char*
rt_cooling_get_species_name(enum rt_cooling_species spec) {

  static const char* rt_cooling_species_names[rt_species_count] = {
      "e","HI","HII","HeI","HeII","HeIII"};

  return rt_cooling_species_names[spec];
}


/**
 * @brief Return a string containing the name of a given #rt_chemistry_element.
 */
__attribute__((always_inline)) INLINE static const char*
rt_chemistry_get_element_name(enum rt_chemistry_element elem) {

  static const char* rt_chemistry_element_names[rt_chemistry_element_count] = {
      "Hydrogen", "Helium"};

  return rt_chemistry_element_names[elem];
}

/**
 * @brief Print the RT model.
 *
 * @param rtp The #rt_props
 */
__attribute__((always_inline)) INLINE static void rt_props_print(
    const struct rt_props* rtp) {

  /* Only the master print */
  if (engine_rank != 0) return;

  message("Radiative transfer scheme: '%s'", "SPH M1closure");

  char messagestring[200] = "Using photon frequency bins: [0.";
  char freqstring[20];
  for (int g = 1; g < RT_NGROUPS; g++) {
    sprintf(freqstring, ", %.3g", rtp->photon_groups[g]);
    strcat(messagestring, freqstring);
  }
  strcat(messagestring, "]");
  message("%s", messagestring);

  if (rtp->use_const_emission_rates) {
    strcpy(messagestring, "Using constant stellar emission rates: [ ");
    for (int g = 0; g < RT_NGROUPS; g++) {
      sprintf(freqstring, "%.3g ", rtp->stellar_const_emission_rates[g]);
      strcat(messagestring, freqstring);
    }
    strcat(messagestring, "]");
    message("%s", messagestring);
  }
}

/**
 * @brief Initialize the global properties of the RT scheme.
 *
 * @param rtp The #rt_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void rt_props_init(
    struct rt_props* rtp, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params,
    struct cosmology* cosmo) {

#ifdef RT_HYDRO_CONTROLLED_INJECTION
  rtp->hydro_controlled_injection = 1;
#else
  rtp->hydro_controlled_injection = 0;
#endif

  rtp->convert_parts_after_zeroth_step = 0;
  rtp->convert_stars_after_zeroth_step = 0;

  if (RT_NGROUPS <= 0) {
    error(
        "You need to run SPHM1RT with at least 1 photon group, "
        "you have %d",
        RT_NGROUPS);
  } else if (RT_NGROUPS == 1) {
    rtp->photon_groups[0] = 0.f;
  } else {
    /* Read in parameters */
    float frequencies[RT_NGROUPS - 1];
    parser_get_param_float_array(params, "SPHM1RT:photon_groups_Hz",
                                 RT_NGROUPS - 1, frequencies);
    float Hz_internal = units_cgs_conversion_factor(us, UNIT_CONV_INV_TIME);
    float Hz_internal_inv = 1.f / Hz_internal;
    for (int g = 0; g < RT_NGROUPS - 1; g++) {
      rtp->photon_groups[g + 1] = frequencies[g] * Hz_internal_inv;
    }
    rtp->photon_groups[0] = 0.f;
  }

  /* Are we using constant emission rates? */
  /* ------------------------------------- */
  rtp->use_const_emission_rates = parser_get_opt_param_int(
      params, "SPHM1RT:use_const_emission_rates", /* default = */ 0);

  if (rtp->use_const_emission_rates) {
    double emission_rates[RT_NGROUPS];
    parser_get_param_double_array(params, "SPHM1RT:star_emission_rates_LSol",
                                  RT_NGROUPS, emission_rates);
    const double unit_power = units_cgs_conversion_factor(us, UNIT_CONV_POWER);
    const double unit_power_inv = 1. / unit_power;
    for (int g = 0; g < RT_NGROUPS; g++) {
      rtp->stellar_const_emission_rates[g] = emission_rates[g] * unit_power_inv;
    }
  } else {
    /* kill the run for now */
    error("SPHM1RT can't run without constant stellar emission rates for now.");
  }

  /* get reduced speed of light in code unit */
  const float cred = parser_get_param_float(params, "SPHM1RT:cred");
  rtp->cred = cred;

  /* get initial opacity in code unit */
  float initialchi[RT_NGROUPS];
  int errorint = parser_get_opt_param_float_array(params, "SPHM1RT:chi",
                                                  RT_NGROUPS, initialchi);

  if (errorint == 0) {
    message("SPHM1RT:chi not found in params");
    for (int g = 0; g < RT_NGROUPS; g++) {
      rtp->initialchi[g] = 0.0f;
    }
  }

  /* get CFL condition */
  const float CFL = parser_get_param_float(params, "SPHM1RT:CFL_condition");
  rtp->CFL_condition = CFL;

  /* Read the total metallicity */
  rtp->initial_metal_mass_fraction_total = parser_get_opt_param_float(
      params, "SPHM1RT:init_abundance_metal", -1.f);

  if (rtp->initial_metal_mass_fraction_total != -1.f) {
    /* Read the individual mass fractions */
    for (int elem = 0; elem < rt_chemistry_element_count; ++elem) {
      char buffer[50];
      sprintf(buffer, "SPHM1RT:init_abundance_%s",
              rt_chemistry_get_element_name((enum rt_chemistry_element)elem));
      rtp->initial_metal_mass_fraction[elem] =
          parser_get_param_float(params, buffer);
    }
  }



  /* Initialize conditional parameters to bogus values */
  rtp->const_stellar_spectrum_max_frequency = -1.;
  rtp->stellar_spectrum_blackbody_T = -1.;

  rtp->stellar_spectrum_type =
      parser_get_param_int(params, "SPHM1RT:stellar_spectrum_type");
  if (rtp->stellar_spectrum_type == 0) {
    /* Constant spectrum: Read additional parameter */
    /* TODO: also translate back to internal units at later. For now, keep it in
     * Hz */
    rtp->const_stellar_spectrum_max_frequency = parser_get_param_float(
        params, "SPHM1RT:stellar_spectrum_const_max_frequency_Hz");
  } else if (rtp->stellar_spectrum_type == 1) {
    /* Blackbody spectrum: Read additional parameter */
    rtp->stellar_spectrum_blackbody_T = parser_get_param_float(
        params, "SPHM1RT:stellar_spectrum_blackbody_temperature_K");
    rtp->stellar_spectrum_blackbody_T /=
        units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  } else {
    error("Selected unknown stellar spectrum type %d",
          rtp->stellar_spectrum_type);
  }

  /* Maximal Star Timestep? */
  /* ---------------------- */
  rtp->stars_max_timestep = parser_get_opt_param_float(
      params, "SPHM1RT:stars_max_timestep", /*default=*/FLT_MAX);
  /* Turn off if negative value given */
  if (rtp->stars_max_timestep < 0.f) rtp->stars_max_timestep = FLT_MAX;
  /* Better safe than sorry */
  if (rtp->stars_max_timestep == 0.f)
    error("You are restricting star time step to 0. That's a no-no.");

  /* After initialisation, print params to screen */
  rt_props_print(rtp);

  /* Print a final message. */
  if (engine_rank == 0) {
    message("Radiative transfer initialized");
  }
}

/**
 * @brief Write an RT properties struct to the given FILE as a
 * stream of bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
__attribute__((always_inline)) INLINE static void rt_struct_dump(
    const struct rt_props* props, FILE* stream) {

  restart_write_blocks((void*)props, sizeof(struct rt_props), 1, stream,
                       "RT props", "RT properties struct");
}

/**
 * @brief Restore an RT properties struct from the given FILE as
 * a stream of bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
__attribute__((always_inline)) INLINE static void rt_struct_restore(
    struct rt_props* props, FILE* stream) {

  restart_read_blocks((void*)props, sizeof(struct rt_props), 1, stream, NULL,
                      "RT properties struct");
}

#endif /* SWIFT_RT_PROPERTIES_SPHM1RT_H */
