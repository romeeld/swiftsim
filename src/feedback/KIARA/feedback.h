/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2022 Doug Rennehan (douglas.rennehan@gmail.com)
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
#ifndef SWIFT_FEEDBACK_KIARA_H
#define SWIFT_FEEDBACK_KIARA_H

#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"
#include "timers.h"
#include "timestep_sync.h"
#include "timestep_sync_part.h"
#include "star_formation.h"

#include <strings.h>

double feedback_get_lum_from_star_particle(const struct spart *sp, 
				           double age,
                                           const struct feedback_props* fb_props);
void feedback_get_ejecta_from_star_particle(const struct spart* sp,
                                            double age,
                                            const struct feedback_props* fb_props,
                                            double dt,
                                            float *N_SNe,
                                            float *ejecta_energy,
                                            float *ejecta_mass,
                                            float *ejecta_unprocessed,
                                            float ejecta_metal_mass[chem5_element_count]);
void feedback_dust_production_condensation(struct spart* sp,
                                           double star_age,
                                           const struct feedback_props* fb_props,
                                           float delta_metal_mass[chemistry_element_count]);
float feedback_life_time(const struct feedback_props* fb_props,
                         const float m, 
                         const float z);
float feedback_imf(const struct feedback_props* fb_props, 
                   const float m);
void feedback_set_turnover_mass(const struct feedback_props* fb_props, 
                                const float z, double* LFLT2);
float feedback_get_turnover_mass(const struct feedback_props* fb_props, 
                                 const float t, const float z);
void feedback_prepare_interpolation_tables(const struct feedback_props* fb_props);

/**
 * @brief Determine the probability of a gas particle being kicked
 *        due to stellar feedback in star forming gas.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param fb_props The feedback properties.
 * @param ti_current The current timestep.
 * @param dt_part The time step of the particle.
 * @param rand_for_sf_wind The random number for the wind generation.
 * @param wind_mass The amount of mass in the wind (code units).
 */
__attribute__((always_inline)) INLINE static double feedback_wind_probability(
    struct part* p, struct xpart* xp, const struct engine* e, 
    const struct cosmology* cosmo,
    const struct feedback_props* fb_props, 
    const integertime_t ti_current, 
    const double dt_part,
    double *rand_for_sf_wind,
    double *wind_mass) {

  return 0.f;
}


/**
 * @brief Kick a gas particle selected for stellar feedback.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param fb_props The feedback properties.
 * @param ti_current The current timestep.
 * @param with_cosmology Is cosmological integration on?
 * @param dt_part The time step of the particle.
 * @param wind_mass The amount of mass in the wind (code units).
 */
__attribute__((always_inline)) INLINE static void feedback_kick_and_decouple_part(
    struct part* p, struct xpart* xp, 
    const struct engine* e, 
    const struct cosmology* cosmo,
    const struct feedback_props* fb_props, 
    const integertime_t ti_current,
    const int with_cosmology,
    const double dt_part,
    const double wind_mass) {};

/**
 * @brief Recouple wind particles.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param with_cosmology Is this a cosmological simulation?
 */
__attribute__((always_inline)) INLINE static void feedback_recouple_part(
    struct part* p, struct xpart* xp, const struct engine* e,
    const int with_cosmology, 
    const struct cosmology* cosmo, const struct unit_system* us,
    const struct feedback_props* fb_props) {

  /* No reason to do this if the decoupling time is zero */
  if (p->feedback_data.decoupling_delay_time > 0.f) {
    const integertime_t ti_step = get_integer_timestep(p->time_bin);
    const integertime_t ti_begin =
        get_integer_time_begin(e->ti_current - 1, p->time_bin);

    /* Get particle time-step */
    double dt_part;
    if (with_cosmology) {
      dt_part =
          cosmology_get_delta_time(e->cosmology, ti_begin, ti_begin + ti_step);
    } else {
      dt_part = get_timestep(p->time_bin, e->time_base);
    }

    p->feedback_data.decoupling_delay_time -= dt_part;

    /* Recouple if below recoupling density */
    const double rho_cgs = hydro_get_physical_density(p, cosmo) *
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
    if (rho_cgs < fb_props->recouple_density_factor * 
      fb_props->recouple_ism_density_cgs) {
      p->feedback_data.decoupling_delay_time = 0.f;
    }

    /* Recouple if it rejoined ISM after some time away */
    const float current_delay_time = fb_props->wind_decouple_time_factor *
      cosmology_get_time_since_big_bang(cosmo, cosmo->a);
    if (p->feedback_data.decoupling_delay_time > 0.2 * current_delay_time && 
        rho_cgs > fb_props->recouple_ism_density_cgs) {
      p->feedback_data.decoupling_delay_time = 0.f;
    }

    /* Firehose wind model: This variable being negative signifies particle should recouple, if it's been a wind long enough */
    if (p->chemistry_data.radius_stream < 0.f) {
      p->feedback_data.decoupling_delay_time = 0.f;
      p->chemistry_data.radius_stream = 0.f;
    }

    /* Don't recouple jets in the ISM (causes issues with SPH) */
    //if (p->feedback_data.decoupling_delay_time <= 0.f && p->feedback_data.number_of_times_decoupled >= 100000 && rho_cgs > fb_props->recouple_ism_density_cgs) {
      //p->feedback_data.decoupling_delay_time += dt_part;
    //}

    /* Here we recouple if needed, and if necessary we also allow cooling again */ 
    if (p->feedback_data.decoupling_delay_time <= 0.f) {
      p->feedback_data.decoupling_delay_time = 0.f;
      p->feedback_data.cooling_shutoff_delay_time = 0.f;

      /* Reset subgrid properties */
      p->cooling_data.subgrid_temp = 0.f;
      p->cooling_data.subgrid_dens = hydro_get_physical_density(p, cosmo);

      /* Make sure to sync the newly coupled part on the timeline */
      timestep_sync_part(p);
    }
    else {
      /* Reset subgrid properties */
      p->cooling_data.subgrid_temp = 0.f;
      p->cooling_data.subgrid_dens = hydro_get_physical_density(p, cosmo);
    }
  } else {
    /* Because we are using floats, always make sure to set exactly zero */
    p->feedback_data.decoupling_delay_time = 0.f;
  }
}

/**
 * @brief Determine if particles that ignore cooling should start cooling again.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param with_cosmology Is this a cosmological simulation?
 */
__attribute__((always_inline)) INLINE static void feedback_ready_to_cool(
    struct part* p, struct xpart* xp, const struct engine* e,
    const struct cosmology* restrict cosmo, const int with_cosmology) {

  /* No reason to do this if the decoupling time is zero */
  if (p->feedback_data.cooling_shutoff_delay_time > 0.f) {
    /* Reset subgrid properties */
    p->cooling_data.subgrid_temp = 0.f;
    //p->cooling_data.subgrid_dens = hydro_get_physical_density(p, cosmo);
    
    const integertime_t ti_step = get_integer_timestep(p->time_bin);
    const integertime_t ti_begin =
        get_integer_time_begin(e->ti_current - 1, p->time_bin);

    /* Get particle time-step */
    double dt_part;
    if (with_cosmology) {
      dt_part =
          cosmology_get_delta_time(e->cosmology, ti_begin, ti_begin + ti_step);
    } else {
      dt_part = get_timestep(p->time_bin, e->time_base);
    }

    p->feedback_data.cooling_shutoff_delay_time -= dt_part;
    if (p->feedback_data.cooling_shutoff_delay_time < 0.f) {
      p->feedback_data.cooling_shutoff_delay_time = 0.f;

      /* Make sure to sync the newly coupled part on the timeline */
      timestep_sync_part(p);
    }
  } else {
    /* Because we are using floats, always make sure to set exactly zero */
    p->feedback_data.cooling_shutoff_delay_time = 0.f;
  }
}

/**
 * @brief Update the properties of a particle due to feedback effects after
 * the cooling was applied.
 *
 * Nothing to do here in the KIARA model.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param with_cosmology Is this a cosmological simulation?
 */
__attribute__((always_inline)) INLINE static void feedback_update_part(
    struct part* p, struct xpart* xp, const struct engine* e,
    const int with_cosmology) {}

/**
 * @brief Reset the gas particle-carried fields related to feedback at the
 * start of a step.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
__attribute__((always_inline)) INLINE static void feedback_reset_part(
    struct part* p, struct xpart* xp) {}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * @param sp The #spart.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static int feedback_is_active(
    const struct spart* sp, const struct engine* e) {

  return 1; /* always */
}

/**
 * @brief Should this particle be doing any DM looping?
 *
 * @param sp The #spart.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static int stars_dm_loop_is_active(
    const struct spart* sp, const struct engine* e) {
  /* Active stars always do the DM loop for the KIARA model */
  return 0;
}

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void feedback_init_spart(
    struct spart* sp) {

  sp->feedback_data.enrichment_weight_inv = 0.f;
  sp->feedback_data.ngb_N = 0;
  sp->feedback_data.ngb_mass = 0.f;
  sp->feedback_data.ngb_rho = 0.f;
  sp->feedback_data.ngb_Z = 0.f;
#ifdef SWIFT_STARS_DENSITY_CHECKS
  sp->has_done_feedback = 0;
#endif
}

/**
 * @brief Returns the length of time since the particle last did
 * enrichment/feedback.
 *
 * @param sp The #spart.
 * @param with_cosmology Are we running with cosmological time integration on?
 * @param cosmo The cosmological model.
 * @param time The current time (since the Big Bang / start of the run) in
 * internal units.
 * @param dt_star the length of this particle's time-step in internal units.
 * @return The length of the enrichment step in internal units.
 */
INLINE static double feedback_get_enrichment_timestep(
    const struct spart* sp, const int with_cosmology,
    const struct cosmology* cosmo, const double time, const double dt_star) {
  /* D. Rennehan: We always use the timestep of the star particle */
  return dt_star;
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 */
__attribute__((always_inline)) INLINE static void feedback_reset_feedback(
    struct spart* sp, const struct feedback_props* feedback_props) {

  /* Zero the distribution weights */
  sp->feedback_data.enrichment_weight = 0.f;

  /* Zero the amount of mass that is distributed */
  sp->feedback_data.mass = 0.f;

  /* Zero the metal enrichment quantities */
  for (int i = 0; i < chemistry_element_count; i++) {
    sp->feedback_data.metal_mass[i] = 0.f;
#if COOLING_GRACKLE_MODE >= 2
    sp->feedback_data.delta_dust_mass[i] = 0.f;
#endif
  }
  sp->feedback_data.total_metal_mass = 0.f;

  /* Zero the energy to inject */
  sp->feedback_data.energy = 0.f;

}

/**
 * @brief Initialises the s-particles feedback props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param feedback_props The properties of the feedback model.
 */
__attribute__((always_inline)) INLINE static void feedback_first_init_spart(
    struct spart* sp, const struct feedback_props* feedback_props) {

  feedback_init_spart(sp);
  sp->feedback_data.SNe_ThisTimeStep = 0.f;
  sp->feedback_data.firehose_radius_stream = 0.f;

}

/**
 * @brief Initialises the s-particles feedback props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param feedback_props The properties of the feedback model.
 */
__attribute__((always_inline)) INLINE static void feedback_prepare_spart(
    struct spart* sp, const struct feedback_props* feedback_props) {}

/**
 * @brief Compute kick velocity for particle sp based on host galaxy properties, in code units
 *
 * @param sp The #spart to consider
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time used value for seeding random number
 */
__attribute__((always_inline)) INLINE static double
feedback_compute_kick_velocity(struct spart* sp, const struct cosmology* cosmo,
    const struct feedback_props *fb_props, const integertime_t ti_current) {

  /* Compute galaxy masses. This is done in the RUNNER files.
   * Therefore, we have access to the gpart */
  double galaxy_stellar_mass_Msun =
      sp->gpart->fof_data.group_stellar_mass;
  if (galaxy_stellar_mass_Msun < fb_props->minimum_galaxy_stellar_mass) 
    galaxy_stellar_mass_Msun = fb_props->minimum_galaxy_stellar_mass;
  galaxy_stellar_mass_Msun *= fb_props->mass_to_solar_mass;

  /* Physical circular velocity km/s from z=0-2 DEEP2 measurements by Dutton+11 */
  /* Dutton+11 eq 6: log (M* / 1e10) = -0.61 + 4.51 log (vdisk / 100) */
  double v_circ_km_s =
      100. * pow(4.0738 * galaxy_stellar_mass_Msun * 1.e-10, 0.221729) *
      pow(cosmo->H / cosmo->H0, 1. / 3.);

  /* Compute using vcirc=sqrt(GM/R), with R(M*) from Ward+24 using CEERS data 
  double galaxy_mass_cgs = sp->gpart->fof_data.group_mass * fb_props->mass_to_solar_mass;
  galaxy_mass_cgs = 1.99e33 * fmax(galaxy_mass_cgs, galaxy_stellar_mass_Msun);
  double reff = 7.1 * pow(1.f + cosmo->z, -0.63) * pow(galaxy_stellar_mass_Msun / 5.e10, 0.19) * 3.086e21;
  double v_circ_km_s = sqrtf(6.67e-8 * galaxy_mass_cgs / reff) * 1.e-5;  // to physical km/s */

  //message("VCIRC: z=%g mg=%g reff=%g ms=%g vc1=%g vc2=%g\n",cosmo->z, galaxy_mass_cgs / 1.99e33, reff/3.086e21, galaxy_stellar_mass_Msun, v_circ_km_s, v_circ_2);

  /* Physical circular velocity km/s from z=0 baryonic tully-fisher relation 
  double galaxy_gas_stellar_mass_Msun =
      sp->gpart->fof_data.group_mass;
  if (galaxy_gas_stellar_mass_Msun <= fb_props->minimum_galaxy_stellar_mass)
      galaxy_gas_stellar_mass_Msun = fb_props->minimum_galaxy_stellar_mass;
  galaxy_gas_stellar_mass_Msun *= fb_props->mass_to_solar_mass;

  const double v_circ_km_s =
      pow(galaxy_gas_stellar_mass_Msun / 102.329, 0.26178) *
      pow(cosmo->H / cosmo->H0, 1. / 3.);*/

  const double rand_for_scatter = random_unit_interval(sp->id, ti_current,
                                      random_number_stellar_feedback_2);

  /* The wind velocity in internal units from FIRE scalings */
  double wind_velocity =
      fb_props->FIRE_velocity_normalization *
      pow(v_circ_km_s / 200., fb_props->FIRE_velocity_slope) *
      ( 1. - fb_props->kick_velocity_scatter +
        2. * fb_props->kick_velocity_scatter * rand_for_scatter) *
      v_circ_km_s *
      fb_props->kms_to_internal *
      /* Note that pj->v_full = a^2 * dx/dt, with x the comoving coordinate. 
       * Thus a physical kick, dv, gets translated into a code velocity kick, a * dv */
      cosmo->a;
 
  const float a_suppress_inv = (1.f + fabs(fb_props->early_wind_suppression_redshift)); 
  if (fb_props->early_wind_suppression_redshift > 0 && cosmo->z > fb_props->early_wind_suppression_redshift) {
    wind_velocity *= cosmo->a * cosmo->a * a_suppress_inv * a_suppress_inv;
  }
  if (fb_props->early_wind_suppression_redshift < 0) {
    wind_velocity *= exp(-pow(cosmo->a * a_suppress_inv, -3.f));
  }

  /* Boost wind speed based on metallicity which governs photon energy output */
  if (fb_props->metal_dependent_vwind == 1 || fb_props->metal_dependent_vwind == 3) {
    float Z_fac = 2.61634;
    const float Z_met = sp->chemistry_data.metal_mass_fraction_total;
    if (Z_met > 1.e-9) Z_fac = pow(10., -0.0029 * pow(log10(Z_met)+9, 2.5) + 0.417694);
    if (Z_fac < 1.f) Z_fac = 1.f;
    wind_velocity *= sqrtf(Z_fac);
  }

  return wind_velocity;
}

/**
 * @brief Prepare a #spart for the feedback task.
 *
 * In KIARA, this function does the stellar evolution for a #spart.
 *
 * @param sp The particle to act upon
 * @param feedback_props The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param phys_const The physical constants in internal units.
 * @param star_age_beg_step The age of the star at the star of the time-step in
 * internal units.
 * @param dt The time-step size of this star in internal units.
 * @param time The physical time in internal units.
 * @param ti_begin The integer time at the beginning of the step.
 * @param with_cosmology Are we running with cosmology on?
 */
__attribute__((always_inline)) INLINE static void feedback_prepare_feedback(
    struct spart* restrict sp, const struct feedback_props* feedback_props,
    const struct cosmology* cosmo, const struct unit_system* us,
    const struct phys_const* phys_const, const double star_age_beg_step,
    const double dt, const double time, const integertime_t ti_begin,
    const int with_cosmology) {

  if (sp->feedback_data.ngb_rho <= 0.) {
    warning("Star %lld has zero neighbor gas density.", sp->id);
    return;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->birth_time == -1.) error("Evolving a star particle that should not!");
#endif

  /* Start by finishing the loops over neighbours */
  const float h = sp->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  sp->feedback_data.ngb_rho *= h_inv_dim;

  const float rho_inv = 1.f / sp->feedback_data.ngb_rho;
  sp->feedback_data.ngb_Z *= h_inv_dim * rho_inv;

  TIMER_TIC;

#ifdef SIMBA_DEBUG_CHECKS
  if (sp->feedback_data.ngb_rho <= 0) {
    warning("Star %lld with mass %g has no neighbors!",
            sp->id, sp->mass);
    return;
  }
#endif

#ifdef SWIFT_DEBUG_CHECKS
  if (star_age_beg_step < 0.f) message("Negative age for a star %g.",star_age_beg_step);
  if (sp->feedback_data.ngb_rho <= 0)
    error("Star %lld with mass %g has no neighbors!",
            sp->id, sp->mass);
  if (sp->count_since_last_enrichment != 0 && engine_current_step > 0)
    error("Computing feedback on a star that should not");
#endif

  /* Properties collected in the stellar density loop. */
  const float ngb_gas_mass = sp->feedback_data.ngb_mass;

  /* Check if there are neighbours, otherwise exit */
  if (ngb_gas_mass == 0.f || sp->density.wcount * pow_dimension(sp->h) < 1e-4) {
    feedback_reset_feedback(sp, feedback_props);
    return;
  }

  /* Update the enrichment weights */
  const float enrichment_weight_inv =
      sp->feedback_data.enrichment_weight_inv;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.enrichment_weight_inv < 0.)
    error("Negative inverse weight!");
#endif

  /* Update the weights used for distribution */
  const float enrichment_weight =
      (enrichment_weight_inv != 0.f) ? 1.f / enrichment_weight_inv : 0.f;
  sp->feedback_data.enrichment_weight = enrichment_weight;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.enrichment_weight < 0.)
    error("Negative weight!");
#endif

  /* Check if this is a newly formed star; if so, set mass to be ejected in wind */
  double eta = feedback_mass_loading_factor(sp->gpart->fof_data.group_stellar_mass, feedback_props->minimum_galaxy_stellar_mass, feedback_props->FIRE_eta_normalization, feedback_props->FIRE_eta_break, feedback_props->FIRE_eta_lower_slope, feedback_props->FIRE_eta_upper_slope);
  /* Boost wind speed based on metallicity which governs photon energy output */
  if (feedback_props->metal_dependent_vwind == 2 || feedback_props->metal_dependent_vwind == 3) {
    float Z_fac = 2.61634;
    float Z_met = sp->chemistry_data.metal_mass_fraction_total;
    if (Z_met > 1.e-9) Z_fac = pow(10., -0.0029 * pow(log10(Z_met)+9, 2.5) + 0.417694);
    if (Z_fac < 1.f) Z_fac = 1.f;
    eta *= Z_fac;
  }

  /* Set total mass to launch and kick velocity for this star */
  if( star_age_beg_step <= dt ) {
     sp->feedback_data.feedback_mass_to_launch = eta * sp->mass;
     sp->feedback_data.feedback_wind_velocity = feedback_compute_kick_velocity(sp, cosmo, feedback_props, ti_begin);
  }
  for (int i=0; i<FEEDBACK_N_KICK_MAX; i++) {
     sp->feedback_data.id_gas_to_be_kicked[i] = -1;
     sp->feedback_data.r2_gas_to_be_kicked[i] = FLT_MAX;
  }

#if COOLING_GRACKLE_MODE >= 2
  /* Compute Habing luminosity of star for use in ISRF (only with Grackle subgrid ISM model) */
  //sp->feedback_data.lum_habing = feedback_get_lum_from_star_particle(sp, star_age_beg_step, feedback_props);
  //message("G0: age %g  Lhabing %g\n",star_age_beg_step*feedback_props->time_to_Myr, sp->feedback_data.lum_habing);
#endif

  /* Zero out mass return for this step */
  sp->feedback_data.mass = 0.f;

  /* Compute the time since the last chemical evolution step was done for this star */
  assert(sp->last_enrichment_time <= cosmo->a);
  double t_since_last = cosmology_get_delta_time_from_scale_factors(
        cosmo, (double)sp->last_enrichment_time, cosmo->a);

  /* If not enough time has passed since last enrichment, then skip it and return */
  if (t_since_last < feedback_props->stellar_enrichment_frequency * star_age_beg_step) return;

  /* Do chem5 chemical evolution model */
  int elem;
  float N_SNe = 0.f;
  float ejecta_energy = 0.f;
  float ejecta_mass = 0.f;
  float ejecta_unprocessed = 0.f;
  float ejecta_metal_mass[chem5_element_count];
  for (elem = 0; elem < chem5_element_count; elem++) ejecta_metal_mass[elem] = 0.f;

  feedback_get_ejecta_from_star_particle(sp, star_age_beg_step, feedback_props, t_since_last, 
		  			 &N_SNe,
                                         &ejecta_energy, 
                                         &ejecta_mass, 
                                         &ejecta_unprocessed, 
                                         ejecta_metal_mass);

  if (isnan(ejecta_mass)) {
    for (elem = 0; elem < chem5_element_count; elem++) {
      message("ejecta_metal_mass[%d]=%g", elem, ejecta_metal_mass[elem]);
    }

    message("[Fe/H] = %g", sp->chemistry_data.metal_mass_fraction[chemistry_element_Fe] / sp->chemistry_data.metal_mass_fraction[chemistry_element_H]);
    message("Z = %g", sp->chemistry_data.metal_mass_fraction_total);

    error("Star particle %lld with mass %g (init_mass %g) is trying to give away NaN mass (Mejecta=%g, Energy=%g, Unprocessed=%g)!",
          sp->id, sp->mass, sp->mass_init, ejecta_mass, ejecta_energy, ejecta_unprocessed);
  }

  if (ejecta_energy < 0.f) {
    warning("Star particle %lld with mass %g (init_mass %g) is trying to give away negative energy (Mejecta=%g, Energy=%g, Unprocessed=%g)!",
          sp->id, sp->mass, sp->mass_init, ejecta_mass, ejecta_energy, ejecta_unprocessed);
    feedback_reset_feedback(sp, feedback_props);
    return;
  }

  if (sp->mass-ejecta_mass < 0.2 * sp->mass_init) {
    warning("Star particle %lld with mass %g is trying to lower its mass past 0.2 of initial (Mejecta=%g)!",
          sp->id, sp->mass, ejecta_mass);
    feedback_reset_feedback(sp, feedback_props);
    return;
  }

  /* D. Rennehan: Do some magic that I still don't understand 
   * https://www.youtube.com/watch?v=cY2xBNWrBZ4
   */
  float dum = 0.f;
  int flag_negative = 0;
  /* Here we can loop over Swift metals because metal_mass_fraction 
   * would be zero for the unique Chem5 metals anyway, and would
   * not activate the condition.
   */
  for (elem = 0; elem < chemistry_element_count; elem++) {
    dum = ejecta_unprocessed * sp->chemistry_data.metal_mass_fraction[elem];
    ejecta_metal_mass[feedback_props->element_index_conversions[elem]] += dum;
    if (ejecta_metal_mass[feedback_props->element_index_conversions[elem]] < 0.f) {
      ejecta_metal_mass[feedback_props->element_index_conversions[elem]] = 0.f;
      flag_negative = 1;
      /* Do not break here, we need the zeroed elements where negative */
    }
  }
  
  /* Check for any remaining that have negative mass after adding unprocessed */
  for (elem = 0; elem < chem5_element_count; elem++) {
    if (ejecta_metal_mass[elem] < 0.f) {
      ejecta_metal_mass[elem] = 0.f;
      flag_negative = 1;
    }
  }

  /* If ANY element ended up negative we recompute everything */
  if (flag_negative) {
    ejecta_mass = 0.f;
    ejecta_metal_mass[chem5_element_Z] = 0.f;
    for (elem = chem5_element_H; elem < chem5_element_Zn; elem++) ejecta_mass += ejecta_metal_mass[elem];
    for (elem = chem5_element_C; elem < chem5_element_Zn; elem++) ejecta_metal_mass[chem5_element_Z] += ejecta_metal_mass[elem];
  }

  /* Now we loop over the Swift metals and set the proper values using the conversion map */
  sp->feedback_data.total_metal_mass = 0.f;
  for (elem = 0; elem < chemistry_element_count; elem++) {
    sp->feedback_data.metal_mass[elem] = ejecta_metal_mass[feedback_props->element_index_conversions[elem]];

    /* Only count real metals in total_metal_mass */
    if (elem != chemistry_element_H && elem != chemistry_element_He) {
      sp->feedback_data.total_metal_mass += ejecta_metal_mass[feedback_props->element_index_conversions[elem]];
    }
  }

#ifdef SIMBA_DEBUG_CHECKS
    if (sp->mass/sp->mass_init<0.2) message("Star particle %lld with mass %g (init %g) is giving away %g Msun and %g erg (%g Msun metals).",
          sp->id, 
          sp->mass, 
          sp->mass_init, 
          ejecta_mass * feedback_props->mass_to_solar_mass, 
          ejecta_energy * feedback_props->energy_to_cgs,
          sp->feedback_data.total_metal_mass * feedback_props->mass_to_solar_mass);
#endif

#if COOLING_GRACKLE_MODE >= 2
  /* Put some of the ejecta metals into dust.  Must be done after chem5->chemistry conversion map is applied */
  if (sp->feedback_data.total_metal_mass > 0.f) {
    feedback_dust_production_condensation(sp, star_age_beg_step, feedback_props, 
                                         sp->feedback_data.metal_mass);
  }
#endif

  /* Compute the total mass to distribute */
  sp->feedback_data.mass = ejecta_mass;

  /* Energy from SNII feedback goes into launching winds; convert units to physical from comoving */
  const float SN_energy_multiplier = 1.f;  // factor above SNII energy to get total energy output by stellar pop 
  //sp->feedback_data.feedback_energy_reservoir += SN_energy_multiplier * ejecta_energy;
  const double E_SN = SN_energy_multiplier * N_SNe * 1.e51 / feedback_props->energy_to_cgs;
  sp->feedback_data.feedback_energy_reservoir += E_SN;

  /* Add early stellar feedback for recently born stellar pops */
  if (feedback_props->early_stellar_feedback_alpha > 0.f && star_age_beg_step > 0.f && star_age_beg_step < 1. / feedback_props->early_stellar_feedback_tfb_inv) { 
      const float alpha_factor = fmax(4. * feedback_props->early_stellar_feedback_alpha - 1.f, 0.f);
      const float dt_factor = pow((star_age_beg_step + dt) * feedback_props->early_stellar_feedback_tfb_inv, alpha_factor) - pow(star_age_beg_step * feedback_props->early_stellar_feedback_tfb_inv, alpha_factor);
      /* Momentum input from ESF from eq. 10 in Keller+22, multiplied by 0.5*v_kick to convert to energy */
      const double E_esf = feedback_props->early_stellar_feedback_alpha * feedback_props->early_stellar_feedback_p0 * sp->mass * dt_factor * 0.5 * sp->feedback_data.feedback_wind_velocity;
      sp->feedback_data.feedback_energy_reservoir += E_esf;
      //message("ESF: z=%g Esn=%g Esf=%g age=%g dtfac=%g vw=%g", cosmo->z, E_SN, E_esf, star_age_beg_step * feedback_props->time_to_Myr, dt_factor, sp->feedback_data.feedback_wind_velocity);
  }

  //if (sp->feedback_data.feedback_mass_to_launch > 0. && sp->feedback_data.feedback_energy_reservoir > 0.) printf("WIND_STAR z=%.5f %lld age=%g mw=%g vw=%g Ew=%g E*=%g ESN=%g Eres=%g\n", cosmo->z, sp->id, star_age_beg_step * feedback_props->time_to_yr * 1.e-6, sp->feedback_data.feedback_mass_to_launch * feedback_props->mass_to_solar_mass, sp->feedback_data.feedback_wind_velocity, 0.5 * sp->feedback_data.feedback_mass_to_launch * sp->feedback_data.feedback_wind_velocity * sp->feedback_data.feedback_wind_velocity, SN_energy_multiplier * ejecta_energy * feedback_props->energy_to_cgs, N_SNe * 1.e51,  sp->feedback_data.feedback_energy_reservoir * feedback_props->energy_to_cgs);
  /* Decrease star mass by amount of mass distributed to gas neighbours */
  sp->mass -= ejecta_mass;

  /* Update last enrichment time */
  sp->last_enrichment_time = cosmo->a;

  /* Set stream radius for firehose particles kicked by this star */
  const float stream_init_density = 0.1;  // in n_cgs units
  const float rho_volumefilling = stream_init_density / feedback_props->rho_to_n_cgs;
  double galaxy_stellar_mass_Msun =
      sp->gpart->fof_data.group_stellar_mass;
  if (galaxy_stellar_mass_Msun < feedback_props->minimum_galaxy_stellar_mass) 
      galaxy_stellar_mass_Msun = feedback_props->minimum_galaxy_stellar_mass;
  galaxy_stellar_mass_Msun *= feedback_props->mass_to_solar_mass;
  const float redge_obs = pow(10.f, 0.34 * log10(galaxy_stellar_mass_Msun) - 2.26) / ((1.f + cosmo->z) * feedback_props->length_to_kpc);  // observed size out to edge of disk galaxies, Buitrago+Trujillo 2024
  sp->feedback_data.firehose_radius_stream = redge_obs;  
  if (sp->group_data.stellar_mass > 0.f && sp->group_data.ssfr > 0.f && eta > 0.f) {
    sp->feedback_data.firehose_radius_stream = min(sqrtf(sp->group_data.ssfr * sp->group_data.stellar_mass * eta / (M_PI * rho_volumefilling * fabs(sp->feedback_data.feedback_wind_velocity))), redge_obs);
  }
  if (sp->feedback_data.firehose_radius_stream <= 1.e-20) sp->feedback_data.firehose_radius_stream = redge_obs;
  if (sp->feedback_data.firehose_radius_stream <= 0.f) error("FIREHOSE stream radius=0! %lld m*=%g ssfr=%g eta=%g robs=%g r=%g\n",sp->id, galaxy_stellar_mass_Msun, sp->group_data.ssfr, eta, redge_obs, sp->feedback_data.firehose_radius_stream);

#if COOLING_GRACKLE_MODE >= 2
  /* Update the number of SNe that have gone off, used in Grackle dust model */
  sp->feedback_data.SNe_ThisTimeStep = N_SNe / t_since_last;  // actually stores SNe rate
#endif

#ifdef SWIFT_STARS_DENSITY_CHECKS
  sp->has_done_feedback = 1;
#endif

  TIMER_TOC(timer_do_star_evol);

}

/**
 * @brief Will this star particle want to do feedback during the next time-step?
 *
 * This is called in the time step task and increases counters of time-steps
 * that have been performed.
 *
 * @param sp The particle to act upon
 * @param feedback_props The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param phys_const The #phys_const.
 * @param time The physical time in internal units.
 * @param with_cosmology Are we running with cosmology on?
 * @param ti_current The current time (in integer)
 * @param time_base The time base.
 */
__attribute__((always_inline)) INLINE static void feedback_will_do_feedback(
    struct spart* sp, const struct feedback_props* feedback_props,
    const int with_cosmology, const struct cosmology* cosmo, const double time,
    const struct unit_system* us, const struct phys_const* phys_const,
    const integertime_t ti_current, const double time_base) { }

void feedback_clean(struct feedback_props* fp);

void feedback_struct_dump(const struct feedback_props* feedback, FILE* stream);

void feedback_struct_restore(struct feedback_props* feedback, FILE* stream);

#ifdef HAVE_HDF5
/**
 * @brief Writes the current model of feedback to the file
 *
 * @param feedback The properties of the feedback scheme.
 * @param h_grp The HDF5 group in which to write.
 */
INLINE static void feedback_write_flavour(struct feedback_props* feedback,
                                          hid_t h_grp) {

  io_write_attribute_s(h_grp, "Feedback Model", "KIARA (decoupled kinetic + chem5 enrichment)");
}
#endif  // HAVE_HDF5

#endif /* SWIFT_FEEDBACK_KIARA_H */
