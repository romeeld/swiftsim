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
#ifndef SWIFT_FEEDBACK_SIMBA_H
#define SWIFT_FEEDBACK_SIMBA_H

#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"
#include "timestep_sync.h"
#include "timestep_sync_part.h"

#include <strings.h>


void compute_stellar_evolution(const struct feedback_props* feedback_props,
                               const struct phys_const* phys_const,
                               const struct cosmology* cosmo, struct spart* sp,
                               const struct unit_system* us, const double age,
                               const double dt, const integertime_t ti_begin);


/**
 * @brief Determine the probability of a gas particle being kicked
 *        due to stellar feedback in star forming gas.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param cosmo The cosmological model.
 * @param fb_props The feedback properties.
 * @param ti_current The current timestep.
 * @param dt_part The time step of the particle.
 * @param rand_for_sf_wind The random number for the wind generation.
 * @param wind_mass The amount of mass in the wind (code units).
 */
__attribute__((always_inline)) INLINE static double feedback_wind_probability(
    struct part* p, struct xpart* xp, 
    const struct engine* e, 
    const struct cosmology* cosmo,
    const struct feedback_props* fb_props, 
    const integertime_t ti_current, 
    const double dt_part,
    double *rand_for_sf_wind,
    double *wind_mass) {

  /* First thing we will do is generate a random number */
  *rand_for_sf_wind = random_unit_interval(p->id, ti_current,
                                           random_number_stellar_feedback_1);

  /* This is done in the RUNNER files. Therefore, we have access
   * to the gpart. */
  double galaxy_stellar_mass = p->galaxy_data.stellar_mass;
  if (galaxy_stellar_mass <= 0.) return 0.;

  const double stellar_mass_this_step = p->sf_data.SFR * dt_part;
  if (stellar_mass_this_step <= 0.) return 0.;

  /* If M* is non-zero, make sure it is at least resolved in the
   * following calculations.
   */
  if (galaxy_stellar_mass < fb_props->minimum_galaxy_stellar_mass) {
    galaxy_stellar_mass = fb_props->minimum_galaxy_stellar_mass;
  }

  /* When early wind suppression is enabled, we alter the minimum
   * stellar mass to be safe.
   */
  if (fb_props->early_wind_suppression_enabled) {
    const double early_minimum_stellar_mass =
        fb_props->early_stellar_mass_norm *
        exp(
          -1. *
          (
            (cosmo->a / fb_props->early_wind_suppression_scale_factor) *
            (cosmo->a / fb_props->early_wind_suppression_scale_factor)
          )
        );
    if (cosmo->a < fb_props->early_wind_suppression_scale_factor) {
      galaxy_stellar_mass = early_minimum_stellar_mass;
    }
  }

  *wind_mass = 
      fb_props->FIRE_eta_normalization * stellar_mass_this_step;
  if (galaxy_stellar_mass < fb_props->FIRE_eta_break) {
    (*wind_mass) *= pow(
      galaxy_stellar_mass / fb_props->FIRE_eta_break, 
      fb_props->FIRE_eta_lower_slope /*-0.317*/
    );
  } else {
    (*wind_mass) *= pow(
      galaxy_stellar_mass / fb_props->FIRE_eta_break, 
      fb_props->FIRE_eta_upper_slope /*-0.761*/
    );
  }

  /* Suppress stellar feedback in the early universe when galaxies are
   * too small. Star formation can destroy unresolved galaxies, so
   * we must suppress the stellar feedback.
   */
  if (fb_props->early_wind_suppression_enabled) {
    if (cosmo->a < fb_props->early_wind_suppression_scale_factor) {
      (*wind_mass) *= pow(cosmo->a / fb_props->early_wind_suppression_scale_factor, 
                       fb_props->early_wind_suppression_slope);
    }
  }

  return 1. - exp(-(*wind_mass) / hydro_get_mass(p));
}

/**
 * @brief Kick a gas particle selected for stellar feedback.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param cosmo The cosmological model.
 * @param fb_props The feedback properties.
 * @param ti_current The current timestep.
 * @param with_cosmology Is cosmological integration on?
 * @param dt_part The time step of the particle.
 * @param wind_mass The amount of mass in the wind (code units).
 */
__attribute__((always_inline)) INLINE static void 
feedback_kick_and_decouple_part(
   struct part* p, struct xpart* xp, 
   const struct engine* e, 
   const struct cosmology* cosmo,
   const struct feedback_props* fb_props, 
   const integertime_t ti_current,
   const int with_cosmology,
   const double dt_part,
   const double wind_mass) {

  const double galaxy_stellar_mass = 
      p->galaxy_data.stellar_mass;
  const double galaxy_stellar_mass_Msun =
      galaxy_stellar_mass * fb_props->mass_to_solar_mass;
  /* This is done in the RUNNER files. Therefore, we have
   * access to the gpart */
  const double galaxy_gas_stellar_mass_Msun = 
      (p->galaxy_data.gas_mass + p->galaxy_data.stellar_mass) * 
          fb_props->mass_to_solar_mass;
  if (galaxy_gas_stellar_mass_Msun <= 0. || galaxy_stellar_mass <= 0.) return;

  /* Physical circular velocity km/s */
  const double v_circ_km_s = 
      pow(galaxy_gas_stellar_mass_Msun / 102.329, 0.26178) *
      pow(cosmo->H / cosmo->H0, 1. / 3.);
  const double rand_for_scatter = random_unit_interval(p->id, ti_current,
                                      random_number_stellar_feedback_2);

  /* The wind velocity in internal units */
  double wind_velocity =
      fb_props->FIRE_velocity_normalization *
      pow(v_circ_km_s / 200., fb_props->FIRE_velocity_slope) *
      (
        1. - fb_props->kick_velocity_scatter + 
        2. * fb_props->kick_velocity_scatter * rand_for_scatter
      ) *
      v_circ_km_s *
      fb_props->kms_to_internal;

  /* Now we have wind_velocity in internal units, determine how much 
   * should go to heating */
  const double u_wind = 0.5 * wind_velocity * wind_velocity;
  
  /* Metal mass fraction (Z) of the gas particle */
  const double Z = p->chemistry_data.metal_mass_fraction_total;

  /* Supernova energy in internal units */
  double u_SN = ((1.e51 * (0.0102778 / fb_props->solar_mass_in_g) * 
                    (p->sf_data.SFR * dt_part / wind_mass)) /
                    (fb_props->kms_to_cms * fb_props->kms_to_cms)) *
		                (fb_props->kms_to_internal * fb_props->kms_to_internal);
  if (Z > 1.e-9) {
    u_SN *= pow(10., -0.0029 * pow(log10(Z) + 9., 2.5) + 0.417694);
  } else {
    u_SN *= 2.61634;
  }

  /* Limit the kinetic energy in the winds to the available SN energy */
  if (u_wind > u_SN) wind_velocity *= sqrt(u_SN / u_wind);

  /* 0.2511886 = pow(10., -0.6) */
  float pandya_slope = 0.f;
  if (galaxy_stellar_mass_Msun > 3.16e10) {
    pandya_slope = -2.1f;
  } else {
    pandya_slope = -0.1f;
  }

  const double f_warm = 
      0.2511886 * pow(galaxy_stellar_mass_Msun / 3.16e10, pandya_slope);
  /* additional 10% removed for cold phase */
  const double hot_wind_fraction = max(0., 0.9 - f_warm);
  const double rand_for_hot = 
      random_unit_interval(p->id, ti_current, random_number_stellar_feedback_3);
  const double rand_for_spread = 
      random_unit_interval(p->id, ti_current, random_number_stellar_feedback);

  /* We want these for logging purposes */
  const double u_init = hydro_get_physical_internal_energy(p, xp, cosmo);
  double u_new = 0.;
  if (u_SN > u_wind && rand_for_hot < hot_wind_fraction) {
    u_new = u_init + (u_SN - u_wind) * (0.5 + rand_for_spread);
  } else {
    u_new = fb_props->cold_wind_internal_energy;
  }

  if (u_new / u_init > 1000) {
    warning("Wind heating too large! T0=%g Tnew=%g fw=%g hwf=%g "
            " TSN=%g Tw=%g vw=%g ms=%g mwind=%g", 
            u_init / fb_props->temp_to_u_factor, 
            u_new / fb_props->temp_to_u_factor, 
            f_warm, 
            hot_wind_fraction, 
            u_SN / fb_props->temp_to_u_factor, 
            u_wind / fb_props->temp_to_u_factor, 
            wind_velocity, 
            p->sf_data.SFR * dt_part, 
            wind_mass);

    u_new = u_init * 1000;
  }

  hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
  hydro_set_drifted_physical_internal_energy(p, cosmo, NULL, u_new);

  const double dir[3] = {
    p->gpart->a_grav[1] * p->gpart->v_full[2] - 
        p->gpart->a_grav[2] * p->gpart->v_full[1],
    p->gpart->a_grav[2] * p->gpart->v_full[0] - 
        p->gpart->a_grav[0] * p->gpart->v_full[2],
    p->gpart->a_grav[0] * p->gpart->v_full[1] - 
        p->gpart->a_grav[1] * p->gpart->v_full[0]
  };
  const double norm = sqrt(
    dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]
  );
  
  /* No norm, no wind */
  if (norm <= 0.) return;
  const double prefactor = cosmo->a * wind_velocity / norm;

  xp->v_full[0] += prefactor * dir[0];
  xp->v_full[1] += prefactor * dir[1];
  xp->v_full[2] += prefactor * dir[2];

  /* Update the signal velocity of the particle based on the velocity kick. */
  hydro_set_v_sig_based_on_velocity_kick(p, cosmo, wind_velocity);

  /* Impose maximal viscosity */
  hydro_diffusive_feedback_reset(p);

  /* Synchronize the particle on the timeline */
  timestep_sync_part(p);

  /* Need time-step for decoupling */
  double dt;
  if (with_cosmology) { 
    const integertime_t ti_step = get_integer_timestep(p->time_bin);
    const integertime_t ti_begin =
      get_integer_time_begin(ti_current - 1, p->time_bin);

    dt = cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
  } 
  else {
    dt = get_timestep(p->time_bin, e->time_base);
  }

  /* Mark to be decoupled */
  p->to_be_decoupled = 1;
  p->to_be_recoupled = 0;
  
  /* Decouple the particles from the hydrodynamics */
  p->feedback_data.decoupling_delay_time = 
      dt + fb_props->wind_decouple_time_factor * 
           cosmology_get_time_since_big_bang(cosmo, cosmo->a);

#ifdef WITH_FOF_GALAXIES
  /* Wind particles are never grouppable. This is done in the
   * RUNNER files. Therefore, we have access to the gpart. */
  p->gpart->fof_data.is_grouppable = 0;
#endif

  /* Wind cannot be star forming */
  if (p->sf_data.SFR > 0.f) {

    /* Record the current time as an indicator of when this particle was last
       star-forming. */
    if (with_cosmology) {
      p->sf_data.SFR = -e->cosmology->a;
    } else {
      p->sf_data.SFR = -e->time;
    }

  }

  /**
   * z pid dt M* Mb vkick vkx vky vkz h x y z vx vy vz T rho v_sig decoupletime 
   * Ndecouple
   */
  const float length_convert = cosmo->a * fb_props->length_to_kpc;
  const float velocity_convert = cosmo->a_inv / fb_props->kms_to_internal;
  const float rho_convert = cosmo->a3_inv * fb_props->rho_to_n_cgs;
  const float u_convert = 
      cosmo->a_factor_internal_energy / fb_props->temp_to_u_factor;

  printf("WIND_LOG %.3f %lld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %g\n",
          cosmo->z,
          p->id, 
          dt_part * fb_props->time_to_Myr,
          galaxy_stellar_mass * fb_props->mass_to_solar_mass,
          galaxy_gas_stellar_mass_Msun,
          wind_velocity / fb_props->kms_to_internal,
          prefactor * dir[0] * velocity_convert,
          prefactor * dir[1] * velocity_convert,
          prefactor * dir[2] * velocity_convert,
          p->h * cosmo->a * fb_props->length_to_kpc,
          p->x[0] * length_convert, 
          p->x[1] * length_convert, 
          p->x[2] * length_convert,
          xp->v_full[0] * velocity_convert, 
          xp->v_full[1] * velocity_convert, 
          xp->v_full[2] * velocity_convert,
          p->u * u_convert, 
          p->rho * rho_convert, 
          p->viscosity.v_sig * velocity_convert,
          p->feedback_data.decoupling_delay_time * fb_props->time_to_Myr, 
          p->feedback_data.number_of_times_decoupled,
          u_new / u_init);
}

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
    const struct cosmology* cosmo,
    const struct feedback_props* fb_props) {

  /* No reason to do this is the decoupling time is zero */
  if (p->decoupled) {
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
    if (p->feedback_data.decoupling_delay_time < 0.f) {
      p->feedback_data.decoupling_delay_time = 0.f;

      /* Make sure to sync the newly coupled part on the timeline */
      timestep_sync_part(p);
    }
  } else {
    /* Because we are using floats, always make sure to set exactly zero */
    p->feedback_data.decoupling_delay_time = 0.f;
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
 * @param with_cosmology Is this a cosmological simulation?
 */
__attribute__((always_inline)) INLINE static void feedback_ready_to_cool(
    struct part* p, struct xpart* xp, const struct engine* e,
    const struct cosmology* restrict cosmo, const int with_cosmology) {

  /* No reason to do this is the decoupling time is zero */
  if (p->feedback_data.cooling_shutoff_delay_time > 0.f) {
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
 * Nothing to do here in the SIMBA model.
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
 * Nothing to do here in the SIMBA model.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
__attribute__((always_inline)) INLINE static void feedback_reset_part(
    struct part* p, struct xpart* xp) {

}

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
  /* Active stars always do the DM loop for the SIMBA model */
  return 0;
}

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void feedback_init_spart(
    struct spart* sp) {

  sp->feedback_data.to_collect.enrichment_weight_inv = 0.f;
  sp->feedback_data.to_collect.ngb_N = 0;
  sp->feedback_data.to_collect.ngb_mass = 0.f;
  sp->feedback_data.to_collect.ngb_rho = 0.f;
  sp->feedback_data.to_collect.ngb_Z = 0.f;
#ifdef SWIFT_STARS_DENSITY_CHECKS
  sp->has_done_feedback = 0;
#endif
}

/**
 * @brief Returns the length of time since the particle last did
 * enrichment/feedback.
 *
 * @param sp The #spart.
 * @param dm_ngb_N the integer number of neighbours from the previous loop
 * @param dm_mean_velocity the mass-weighted (unnormalized) three components of
 * velocity
 */
INLINE static void feedback_intermediate_density_normalize(
    struct spart* sp, const int dm_ngb_N, float dm_mean_velocity[3]) {

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
  } else {
    return time - (double)sp->last_enrichment_time;
  }
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 */
__attribute__((always_inline)) INLINE static void feedback_reset_feedback(
    struct spart* sp, const struct feedback_props* feedback_props) {

  /* Zero the distribution weights */
  sp->feedback_data.to_distribute.enrichment_weight = 0.f;

  /* Zero the amount of mass that is distributed */
  sp->feedback_data.to_distribute.mass = 0.f;

  /* Zero the metal enrichment quantities */
  for (int i = 0; i < chemistry_element_count; i++) {
    sp->feedback_data.to_distribute.metal_mass[i] = 0.f;
  }
  sp->feedback_data.to_distribute.total_metal_mass = 0.f;
  sp->feedback_data.to_distribute.mass_from_AGB = 0.f;
  sp->feedback_data.to_distribute.metal_mass_from_AGB = 0.f;
  sp->feedback_data.to_distribute.mass_from_SNII = 0.f;
  sp->feedback_data.to_distribute.metal_mass_from_SNII = 0.f;
  sp->feedback_data.to_distribute.mass_from_SNIa = 0.f;
  sp->feedback_data.to_distribute.metal_mass_from_SNIa = 0.f;
  sp->feedback_data.to_distribute.Fe_mass_from_SNIa = 0.f;

  /* Zero the energy to inject */
  sp->feedback_data.to_distribute.energy = 0.f;

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
 * @brief Prepare a #spart for the feedback task.
 *
 * In SIMBA, this function evolves the stellar properties of a #spart.
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

  if (sp->feedback_data.to_collect.ngb_rho <= 0.) {
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

  sp->feedback_data.to_collect.ngb_rho *= h_inv_dim;

  const float rho_inv = 1.f / sp->feedback_data.to_collect.ngb_rho;
  sp->feedback_data.to_collect.ngb_Z *= h_inv_dim * rho_inv;

  /* Compute amount of enrichment and feedback that needs to be done in this
   * step */
  compute_stellar_evolution(feedback_props, phys_const, cosmo, sp, us,
                            star_age_beg_step, dt, ti_begin);

  /* Decrease star mass by amount of mass distributed to gas neighbours */
  sp->mass -= sp->feedback_data.to_distribute.mass;

  /* Mark this is the last time we did enrichment */
  if (with_cosmology)
    sp->last_enrichment_time = cosmo->a;
  else
    sp->last_enrichment_time = time;

#ifdef SWIFT_STARS_DENSITY_CHECKS
  sp->has_done_feedback = 1;
#endif
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

  io_write_attribute_s(h_grp, "Feedback Model", "SIMBA (decoupled kinetic)");
}
#endif  // HAVE_HDF5

#endif /* SWIFT_FEEDBACK_SIMBA_H */
