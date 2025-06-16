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
                                double *N_SNe,
                                double *ejecta_energy,
                                double *ejecta_mass,
                                double *ejecta_unprocessed,
                                double ejecta_metal_mass[chem5_element_count]);
void feedback_dust_production_condensation(struct spart* sp,
                              double star_age,
                              const struct feedback_props* fb_props,
                              double delta_metal_mass[chemistry_element_count]);
double feedback_life_time(const struct feedback_props* fb_props,
                          const double m, 
                          const double z);
double feedback_imf(const struct feedback_props* fb_props, 
                    const double m);
void feedback_set_turnover_mass(const struct feedback_props* fb_props, 
                                const double z, double* LFLT2);
double feedback_get_turnover_mass(const struct feedback_props* fb_props, 
                                  const double t, const double z);
void feedback_prepare_interpolation_tables(
                                  const struct feedback_props* fb_props);

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
__attribute__((always_inline)) INLINE static 
void feedback_kick_and_decouple_part(
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
__attribute__((always_inline)) INLINE static 
void feedback_recouple_set_flags(struct part* p, 
                                 const struct cosmology* cosmo) {

  p->feedback_data.decoupling_delay_time = 0.f;
  p->chemistry_data.radius_stream = 0.f;

  /* Flag for recoupling */
  p->to_be_recoupled = 1;

  /* Reset subgrid properties */
  p->cooling_data.subgrid_temp = 0.f;
  p->cooling_data.subgrid_dens = hydro_get_physical_density(p, cosmo);

  /* Make sure to sync the newly coupled part on the timeline */
  timestep_sync_part(p);
}

/**
 * @brief Recouple wind particles.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param with_cosmology Is this a cosmological simulation?
 * @param cosmo The cosmology of the simulation.
 * @param fb_props The #feedback_props feedback parameters.
 */
__attribute__((always_inline)) INLINE static void feedback_recouple_part(
    struct part* p, struct xpart* xp, const struct engine* e,
    const int with_cosmology, 
    const struct cosmology* cosmo,
    const struct feedback_props* fb_props) {

  if (p->decoupled) {
    const integertime_t ti_step = get_integer_timestep(p->time_bin);
    const integertime_t ti_begin =
        get_integer_time_begin(e->ti_current - 1, p->time_bin);

    /* Get particle time-step */
    double dt_part;
    if (with_cosmology) {
      dt_part =
          cosmology_get_delta_time(e->cosmology, ti_begin, ti_begin + ti_step);
    } 
    else {
      dt_part = get_timestep(p->time_bin, e->time_base);
    }

    /* Decrement the counter */
    p->feedback_data.decoupling_delay_time -= dt_part;

    /**
     * Recouple under 3 conditions:
     * (1) Below the density threshold.
     * (2) If the stream radius is negative.
     * (3) If the timer has run out.
     */
    const double rho_nH_cgs = 
        hydro_get_physical_density(p, cosmo) * fb_props->rho_to_n_cgs;
    const double rho_recouple_cgs = 
        fb_props->recouple_density_factor * 
            fb_props->recouple_ism_density_nH_cgs;

    const int recouple = (p->feedback_data.decoupling_delay_time <= 0.f ||
                          p->chemistry_data.radius_stream < 0.f ||
                          rho_nH_cgs < rho_recouple_cgs);

    if (recouple) {
      feedback_recouple_set_flags(p, cosmo);
    }
    else {
      /* Reset subgrid properties if decoupled for safety */
      p->cooling_data.subgrid_temp = 0.f;
      p->cooling_data.subgrid_dens = hydro_get_physical_density(p, cosmo);
    }
  }
}

/**
 * @brief Sets the wind direction vector for feedback kicks
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param with_cosmology Is this a cosmological simulation?
 * @param cosmo The cosmology of the simulation.
 * @param fb_props The #feedback_props feedback parameters.
 */
__attribute__((always_inline)) INLINE static void feedback_set_wind_direction(
    struct part* p, struct xpart* xp, const struct engine* e,
    const int with_cosmology, 
    const struct cosmology* cosmo,
    const struct feedback_props* fb_props) {

  p->feedback_data.wind_direction[0] = 
      p->gpart->a_grav[1] * p->gpart->v_full[2] -
      p->gpart->a_grav[2] * p->gpart->v_full[1];
  p->feedback_data.wind_direction[1] = 
      p->gpart->a_grav[2] * p->gpart->v_full[0] -
      p->gpart->a_grav[0] * p->gpart->v_full[2];
  p->feedback_data.wind_direction[2] = 
      p->gpart->a_grav[0] * p->gpart->v_full[1] -
      p->gpart->a_grav[1] * p->gpart->v_full[0];
}

/**
 * @brief Determine if particles that ignore cooling should start cooling again.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param cosmo The cosmological information of the simulation.
 * @param with_cosmology Is this a cosmological simulation?
 */
__attribute__((always_inline)) INLINE static void feedback_ready_to_cool(
    struct part* p, struct xpart* xp, const struct engine* e,
    const struct cosmology* restrict cosmo, const int with_cosmology) {

  /* No reason to do this if the decoupling time is zero */
  if (p->feedback_data.cooling_shutoff_delay_time > 0.f) {

    /* Reset subgrid properties */
    p->cooling_data.subgrid_temp = 0.f;
    p->cooling_data.subgrid_dens = hydro_get_physical_density(p, cosmo);
    
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
  } 
  else {
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

  return e->step <= 0 ||
         ((sp->birth_time != -1.) && (sp->count_since_last_enrichment == 0));
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

  /* Default to not suppression the mass loading in the winds */
  sp->feedback_data.eta_suppression_factor = 1.f;
  sp->feedback_data.kernel_wt_sum = 0.f;
  sp->feedback_data.wind_wt_sum = 0.f;
  sp->feedback_data.ngb_mass = 0.f;
  sp->feedback_data.wind_ngb_mass = 0.f;

  /* Check reservoirs each time-step for out-of-bounds values */
  if (sp->feedback_data.mass_to_launch < 0.f) {
    sp->feedback_data.mass_to_launch = 0.f;
  }

  if (sp->feedback_data.physical_energy_reservoir < 0.) {
    sp->feedback_data.physical_energy_reservoir = 0.;
  }

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
  
  if (with_cosmology) {
    return cosmology_get_delta_time_from_scale_factors(
        cosmo, (double)sp->last_enrichment_time, cosmo->a);
  } 
  else {
    return time - (double)sp->last_enrichment_time;
  }
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 */
__attribute__((always_inline)) INLINE static void feedback_reset_feedback(
    struct spart* sp, const struct feedback_props* feedback_props) {

  /* Zero the amount of mass that is distributed */
  sp->feedback_data.mass = 0.;

  /* Zero the metal enrichment quantities */
  for (int i = 0; i < chemistry_element_count; i++) {
    sp->feedback_data.metal_mass[i] = 0.;
#if COOLING_GRACKLE_MODE >= 2
    sp->feedback_data.delta_dust_mass[i] = 0.;
#endif
  }
  sp->feedback_data.total_metal_mass = 0.;

  /* Zero the energy to inject */
  sp->feedback_data.energy = 0.;

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
  sp->feedback_data.SNe_ThisTimeStep = 0.;
  sp->feedback_data.SNe_Total = 0.;
  sp->feedback_data.firehose_radius_stream = 0.f;
  sp->feedback_data.mass_to_launch = 0.f;
  sp->feedback_data.total_mass_kicked = 0.f;
  sp->feedback_data.wind_velocity = 0.f;
  sp->feedback_data.physical_energy_reservoir = 0.;
  sp->feedback_data.N_launched = 0;
}

/**
 * @brief Initialises the particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions or assignments between the particle
 * and extended particle fields.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void feedback_first_init_part(
    struct part *restrict p, struct xpart *restrict xp) {

  p->feedback_data.decoupling_delay_time = 0.f;
  p->feedback_data.number_of_times_decoupled = 0;
  p->feedback_data.cooling_shutoff_delay_time = 0.f;
  p->feedback_data.kick_id = -1;
  for (int i = 0; i < 3; i++) p->feedback_data.wind_direction[i] = 0.f;
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
 * @brief Compute kick velocity for particle sp based on host galaxy properties, 
 *        in code units
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
  double galaxy_stellar_mass_Msun = sp->galaxy_data.stellar_mass;
  if (galaxy_stellar_mass_Msun < fb_props->minimum_galaxy_stellar_mass) {
    galaxy_stellar_mass_Msun = fb_props->minimum_galaxy_stellar_mass;
  }

  galaxy_stellar_mass_Msun *= fb_props->mass_to_solar_mass;

  /* Physical circular velocity km/s from z=0-2 DEEP2 
     measurements by Dutton+11 */

  /* Dutton+11 eq 6: log (M* / 1e10) = -0.61 + 4.51 log (vdisk / 100) */
  const float v_circ_km_s =
      100.f * powf(4.0738f * galaxy_stellar_mass_Msun * 1.e-10f, 0.221729f) *
          pow(cosmo->H / cosmo->H0, 1.f / 3.f);

  const float rand_for_scatter = 
      random_unit_interval(sp->id, ti_current, 
                           random_number_stellar_feedback_2);

  /* The wind velocity in internal units and COMOVING from FIRE scalings */
  float wind_velocity =
      fb_props->FIRE_velocity_normalization *
      powf(v_circ_km_s / 200.f, fb_props->FIRE_velocity_slope) *
      ( 1.f - fb_props->kick_velocity_scatter +
        2.f * fb_props->kick_velocity_scatter * rand_for_scatter) *
      v_circ_km_s *
      fb_props->kms_to_internal *
      /* Note that pj->v_full = a^2 * dx/dt, with x the comoving coordinate. 
       * Thus a physical kick, dv, gets translated into a code velocity kick, 
       * a * dv */
      cosmo->a;
 
  const float a_suppress_inv = 
      (1.f + fabs(fb_props->wind_velocity_suppression_redshift)); 
  if (fb_props->wind_velocity_suppression_redshift > 0 && 
          cosmo->z > fb_props->wind_velocity_suppression_redshift) {
    wind_velocity *= cosmo->a * cosmo->a * a_suppress_inv * a_suppress_inv;
  }
  else if (fb_props->wind_velocity_suppression_redshift < 0) {
    wind_velocity *= expf(-powf(cosmo->a * a_suppress_inv, -3.f));
  }

  /* internal, COMOVING units */
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

  if (sp->feedback_data.ngb_mass <= 0.f) {
    warning("Star %lld has zero neighbor gas density.", sp->id);
    return;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->birth_time == -1.) error("Evolving a star particle that should not!");
#endif

  TIMER_TIC;

#if COOLING_GRACKLE_MODE >= 2
  /* Compute Habing luminosity of star for use in ISRF 
     (only with Grackle subgrid ISM model) */
  /*sp->feedback_data.lum_habing = 
      feedback_get_lum_from_star_particle(sp, star_age_beg_step, feedback_props);
    message("G0: age %g  Lhabing %g\n",
            star_age_beg_step * feedback_props->time_to_Myr, 
            sp->feedback_data.lum_habing);
    */
#endif

  /* Zero out mass and energy return for this step */
  sp->feedback_data.mass = 0.;
  sp->feedback_data.energy = 0.;

  /* Do chem5 chemical evolution model */
  int elem;
  double N_SNe = 0.;
  double ejecta_energy = 0.;
  double ejecta_mass = 0.;
  double ejecta_unprocessed = 0.;
  double ejecta_metal_mass[chem5_element_count];
  for (elem = 0; elem < chem5_element_count; elem++) {
    ejecta_metal_mass[elem] = 0.;
  }

  feedback_get_ejecta_from_star_particle(sp, 
                                         star_age_beg_step, 
                                         feedback_props, 
                                         dt, 
                                         &N_SNe,
                                         &ejecta_energy, 
                                         &ejecta_mass, 
                                         &ejecta_unprocessed, 
                                         ejecta_metal_mass);

  if (isnan(ejecta_mass)) {
    for (elem = 0; elem < chem5_element_count; elem++) {
      message("ejecta_metal_mass[%d]=%g", elem, ejecta_metal_mass[elem]);
    }

    message("[Fe/H] = %g", 
            sp->chemistry_data.metal_mass_fraction[chemistry_element_Fe] / 
                sp->chemistry_data.metal_mass_fraction[chemistry_element_H]);
    message("Z = %g", sp->chemistry_data.metal_mass_fraction_total);

    error("Star particle %lld with mass %g (init_mass %g) is trying to give "
          "away NaN mass (Mejecta=%g, Energy=%g, Unprocessed=%g)!",
          sp->id, sp->mass, sp->mass_init, ejecta_mass, ejecta_energy, 
          ejecta_unprocessed);
  }

  if (ejecta_energy < 0.f) {
    warning("Star particle %lld with mass %g (init_mass %g) is trying to give "
            "away negative energy (Mejecta=%g, Energy=%g, Unprocessed=%g)!",
          sp->id, sp->mass, sp->mass_init, ejecta_mass, ejecta_energy, 
          ejecta_unprocessed);
    feedback_reset_feedback(sp, feedback_props);
    return;
  }

  if (sp->mass-ejecta_mass < 0.2 * sp->mass_init) {
    warning("Star particle %lld with mass %g is trying to lower its mass "
            "past 0.2 of initial (Mejecta=%g)!",
            sp->id, sp->mass, ejecta_mass);
    feedback_reset_feedback(sp, feedback_props);
    return;
  }

  /* Now set up mass to launch in wind */
  /* Check if this is a newly formed star; if so, set mass to be ejected in 
   * wind */
  float eta = 
      feedback_mass_loading_factor(sp->galaxy_data.stellar_mass, 
                                   feedback_props->minimum_galaxy_stellar_mass, 
                                   feedback_props->FIRE_eta_normalization, 
                                   feedback_props->FIRE_eta_break, 
                                   feedback_props->FIRE_eta_lower_slope, 
                                   feedback_props->FIRE_eta_upper_slope);

  /* Suppress eta based on the redshift, if desired */
  const float z_suppress = feedback_props->wind_eta_suppression_redshift;
  if (z_suppress > 0.f && cosmo->z > z_suppress) {
    const float a_suppress = 1.f / (1.f + z_suppress);
    sp->feedback_data.eta_suppression_factor = 
        (cosmo->a * cosmo->a) / (a_suppress * a_suppress);
  }
  else {
    /* no supression */
    sp->feedback_data.eta_suppression_factor = 1.f;
  }


  /**
   * Newly born star particles only ever kick out gas once in their lifetime
   * since the mass loading is set based on the stellar mass formed. Therefore,
   * check if the launched flag has not been set. It has to do at least one
   * wind kick event!
   * 
   * If true, then set up a wind launch event with a fixed mass loading, and
   * wind mass based on the initial stellar mass. 
   * 
   * That wind may launch over several time-steps if the amount of mass to kick
   * if larger than half of the gas reservoir within the kernel. If the amount
   * to kick is above 50% of the gas mass in the kernel, the code limits the 
   * amount kicked to 50% of the gas mass. The next time-step, the star should
   * then kick the remaining amount of gas, if possible. 
   */
  //if (!sp->feedback_data.launched && eta > 0.f) {
  const float wind_mass = eta * sp->mass_init;
  const float total_mass_kicked = sp->feedback_data.total_mass_kicked;
  if (N_SNe > 0.f && total_mass_kicked < wind_mass && eta > 0.f) {
    /* velocity in internal units which is a^2*comoving, or a*physical */
    float v_internal = 
        feedback_compute_kick_velocity(sp, cosmo, feedback_props, ti_begin);

    /* Boost wind speed based on metallicity which governs 
      * photon energy output */
    float Z_fac = 1.f;
    const int vwind_boost_flag = feedback_props->metal_dependent_vwind;
    if (vwind_boost_flag != kiara_metal_boosting_off) {
      Z_fac = 2.61634f;
      const float Z_met = sp->chemistry_data.metal_mass_fraction_total;
      if (Z_met > 1.e-9f) {
        Z_fac = powf(10.f, 
                    -0.0029f * powf(log10f(Z_met) + 9.f, 2.5f) + 0.417694f);
      }
  
      Z_fac = sqrtf(max(Z_fac, 1.f));
    }

    switch (vwind_boost_flag) {
      case kiara_metal_boosting_vwind:
        v_internal *= Z_fac;
        break;
      case kiara_metal_boosting_eta:
        eta *= Z_fac *  Z_fac;
        break;
      case kiara_metal_boosting_both:
        v_internal *= Z_fac;
        eta *= Z_fac * Z_fac;
        break;
    }

    /* Add SNII energy to the energy reservoir available for launching winds */
    const float scaling = 
        sqrtf(feedback_props->SNII_energy_multiplier) * Z_fac;
    const double E_SNe = 1.e51 * N_SNe / feedback_props->energy_to_cgs;

    /* Comoving energy reservoir */
    sp->feedback_data.physical_energy_reservoir += E_SNe * scaling;

    /* Set kick for this star */
    sp->feedback_data.wind_velocity = v_internal;
  
    /* Compute mass to launch in this timestep, based on 
     * currently available SNII energy */
    const double wind_energy_phys = 
        0.5 * sp->mass_init * v_internal * v_internal * cosmo->a2_inv;
    const double eta_max_this_timestep = 
        sp->feedback_data.physical_energy_reservoir / wind_energy_phys;

    /* How much can the star still kick? */
    const double wind_mass_max = eta_max_this_timestep * sp->mass_init;
    const double wind_mass_left = max(wind_mass - total_mass_kicked, 0.);

    const double mass_to_launch = min(wind_mass_max, wind_mass_left);

    sp->feedback_data.mass_to_launch = mass_to_launch;

#ifdef KIARA_DEBUG_CHECKS
    message("ETA: z=%g id=%lld age=%g Eres=%g dE=%g NSNe=%g NSNtot=%g eta=%g "
            "max=%g tot=%g mlaunch=%g Ntot=%d",
            cosmo->z, 
            sp->id, 
            star_age_beg_step * feedback_props->time_to_Myr, 
            sp->feedback_data.physical_energy_reservoir *
                feedback_props->energy_to_cgs, 
            1.e51 * N_SNe * scaling, 
            N_SNe, 
	          sp->mass_init * feedback_props->mass_to_solar_mass / 80.f,  
            /* 1 SNII for ~80 Mo for Kroupa/Chabrier IMF */
            mass_to_launch / sp->mass_init, 
            eta_max_this_timestep, 
            eta, 
            sp->feedback_data.mass_to_launch,
	          sp->feedback_data.N_launched);
#endif

    /* Set stream radius for firehose particles kicked by this star */
    const double stream_init_density = 0.1; /* n_H units CGS */
    const double rho_volumefilling = 
        stream_init_density / feedback_props->rho_to_n_cgs;
    float galaxy_stellar_mass_Msun = sp->galaxy_data.stellar_mass;
    const float min_gal_mass = feedback_props->minimum_galaxy_stellar_mass;
    if (galaxy_stellar_mass_Msun < min_gal_mass) {
      galaxy_stellar_mass_Msun = min_gal_mass;
    }
    galaxy_stellar_mass_Msun *= feedback_props->mass_to_solar_mass;

    /* stream size = 2 * comoving effective size of disk galaxies 
     * (Ward+2024 CEERS) */
    const float redge_obs = 2.f * 
        7.1f * pow(cosmo->a, 0.63f) * 
          pow(galaxy_stellar_mass_Msun / 5.e10, 0.16f) /
            feedback_props->length_to_kpc * cosmo->a_inv; 
    sp->feedback_data.firehose_radius_stream = redge_obs;

    if (sp->galaxy_data.stellar_mass > 0.f && sp->galaxy_data.ssfr > 0.f && 
          eta > 0.f) {
      sp->feedback_data.firehose_radius_stream = 
          min(
            sqrtf(sp->galaxy_data.ssfr * sp->galaxy_data.stellar_mass * eta / 
                  (M_PI * rho_volumefilling * 
                    fabs(sp->feedback_data.wind_velocity))), redge_obs);
    }

    /* Stream cannot be smaller than the smoothing length */
    sp->feedback_data.firehose_radius_stream = 
        fmax(sp->feedback_data.firehose_radius_stream, sp->h);
  }

  /* D. Rennehan: Do some magic that I still don't understand 
   */
  double dum = 0.;
  int flag_negative = 0;
  /* Here we can loop over Swift metals because metal_mass_fraction 
   * would be zero for the unique Chem5 metals anyway, and would
   * not activate the condition.
   */
  for (elem = 0; elem < chemistry_element_count; elem++) {
    dum = ejecta_unprocessed * sp->chemistry_data.metal_mass_fraction[elem];
    const int elem_conv = feedback_props->element_index_conversions[elem];
    ejecta_metal_mass[elem_conv] += dum;
    if (ejecta_metal_mass[elem_conv] < 0.) {
      ejecta_metal_mass[elem_conv] = 0.;
      flag_negative = 1;
      /* Do not break here, we need the zeroed elements where negative */
    }
  }
  
  /* Check for any remaining that have negative mass after adding unprocessed */
  for (elem = 0; elem < chem5_element_count; elem++) {
    if (ejecta_metal_mass[elem] < 0.) {
      ejecta_metal_mass[elem] = 0.;
      flag_negative = 1;
    }
  }

  /* If ANY element ended up negative we recompute everything */
  if (flag_negative) {
    ejecta_mass = 0.;
    ejecta_metal_mass[chem5_element_Z] = 0.;
    for (elem = chem5_element_H; elem < chem5_element_Zn; elem++) {
      ejecta_mass += ejecta_metal_mass[elem];
    }

    for (elem = chem5_element_C; elem < chem5_element_Zn; elem++) {
      ejecta_metal_mass[chem5_element_Z] += ejecta_metal_mass[elem];
    }
  }

  /* Now we loop over the Swift metals and set the proper values using the 
     conversion map */
  sp->feedback_data.total_metal_mass = ejecta_metal_mass[chem5_element_Z];
  for (elem = 0; elem < chemistry_element_count; elem++) {
    sp->feedback_data.metal_mass[elem] = 
        ejecta_metal_mass[feedback_props->element_index_conversions[elem]];
  }

#if COOLING_GRACKLE_MODE >= 2
  /* Put some of the ejecta metals into dust.  Must be done after 
     chem5->chemistry conversion map is applied */
  if (sp->feedback_data.total_metal_mass > 0.) {
    feedback_dust_production_condensation(sp, star_age_beg_step, feedback_props, 
                                          sp->feedback_data.metal_mass);
  }
#endif

  /* Compute the total mass to distribute */
  sp->feedback_data.mass = ejecta_mass;
  sp->feedback_data.energy = ejecta_energy;

  /* Decrease star mass by amount of mass distributed to gas neighbours */
  sp->mass -= ejecta_mass;

  /* Mark this is the last time we did enrichment */
  sp->last_enrichment_time = (with_cosmology) ? cosmo->a : time;

#if COOLING_GRACKLE_MODE >= 2
  /* Update the number of SNe that have gone off, used in Grackle dust model. 
     Actually stores SNe rate */
  sp->feedback_data.SNe_ThisTimeStep = N_SNe / dt; 
  sp->feedback_data.SNe_Total += N_SNe;
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
    const integertime_t ti_current, const double time_base) { 

  /* Special case for new-born stars */
  if (with_cosmology) {
    if (sp->birth_scale_factor == (float)cosmo->a) {

      /* Set the counter to "let's do enrichment" */
      sp->count_since_last_enrichment = 0;

      /* Ok, we are done. */
      return;
    }
  } else {
    if (sp->birth_time == (float)time) {

      /* Set the counter to "let's do enrichment" */
      sp->count_since_last_enrichment = 0;

      /* Ok, we are done. */
      return;
    }
  }

  /* Calculate age of the star at current time */
  double age_of_star;
  if (with_cosmology) {
    age_of_star = cosmology_get_delta_time_from_scale_factors(
        cosmo, (double)sp->birth_scale_factor, cosmo->a);
  } else {
    age_of_star = time - (double)sp->birth_time;
  }

  /* Is the star still young? */
  if (age_of_star < feedback_props->stellar_evolution_age_cut) {

    /* Set the counter to "let's do enrichment" */
    sp->count_since_last_enrichment = 0;

  } else {

    /* Increment counter */
    sp->count_since_last_enrichment++;

    if ((sp->count_since_last_enrichment %
         feedback_props->stellar_evolution_sampling_rate) == 0) {

      /* Reset counter */
      sp->count_since_last_enrichment = 0;
    }
  }
}

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

  io_write_attribute_s(h_grp, "Feedback Model", "KIARA "
                       "(decoupled kinetic + chem5 enrichment)");
}
#endif  // HAVE_HDF5

#endif /* SWIFT_FEEDBACK_KIARA_H */
