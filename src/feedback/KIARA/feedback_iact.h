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
#ifndef SWIFT_KIARA_FEEDBACK_IACT_H
#define SWIFT_KIARA_FEEDBACK_IACT_H

/* Local includes */
#include "random.h"
#include "timestep_sync_part.h"
#include "tools.h"
#include "tracers.h"
#include <assert.h>

#define KICK_RADIUS_OVER_H 0.5f
#define MAX_FRAC_OF_KERNEL_TO_LAUNCH 0.5f

/**
 * @brief Compute the mean DM velocity around a star. (non-symmetric).
 *
 * @param si First sparticle.
 * @param gj Second particle (not updated).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_dm_vel_sum(struct spart *si, const struct gpart *gj,
                                       int *dm_ngb_N,
                                       float dm_mean_velocity[3]) {}

/**
 * @brief Compute the DM velocity dispersion around a star. (non-symmetric).
 *
 * @param si First sparticle.
 * @param gj Second particle.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_dm_vel_disp(struct spart *si,
                                        const struct gpart *gj,
                                        const float dm_mean_velocity[3]) {}

/**
 * @brief Compute customized kernel weight for feedback
 *
 * @param pj gas particle.
 * @param wi SPH kernel weight at location of pj
 */
__attribute__((always_inline)) INLINE static float
feedback_kernel_weight(const struct part *pj, const float wi, const float ui)
{
  /* If it's beyond the kick radius, then the weighting is zero */
  if (ui >= KICK_RADIUS_OVER_H) return 0.f;

  /* Weight towards higher SFR particles. As SFR->0, SFR_wi->wi and
  * then radial weighting returns to normal. */
  float weight = (pj->sf_data.SFR > 0.f) ? wi + pj->sf_data.SFR : wi;
  weight *= hydro_get_mass(pj);
  return weight;
}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param xpj Extra particle data (not updated).
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time value
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_density(const float r2, const float dx[3],
                                    const float hi, const float hj,
                                    struct spart *si, const struct part *pj,
                                    const struct xpart *xpj,
                                    const struct cosmology *cosmo,
                                    const struct feedback_props *fb_props,
                                    const integertime_t ti_current) {

  /* Do not count winds in the density */
  if (pj->decoupled) return;

  const float rho = hydro_get_comoving_density(pj);
  if (rho <= 0.f) return;

  /* Get the gas mass. */
  const float mj = hydro_get_mass(pj);

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;

  float wi;
  kernel_eval(ui, &wi);

  /* Add mass of pj to neighbour mass of si  */
  si->feedback_data.ngb_mass += mj;

  /* sum(mj * wj) */
  si->feedback_data.kernel_wt_sum += mj * wi;

  /* If pj is being kicked in this step, don't kick again */
#ifndef BLACK_HOLES_NONE
  if (pj->black_holes_data.swallow_id > -1) return;
#endif
  if (pj->feedback_data.kick_id > -1) return;

  /* Sum up the weights for normalizing the kernel later */
  const float wt = feedback_kernel_weight(pj, wi, ui);
  si->feedback_data.wind_wt_sum += wt;
  if (wt > 0.f) si->feedback_data.wind_ngb_mass += mj;

}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep1(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  const struct spart *si, struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const integertime_t ti_current) {

  /* No need to even check anything else if there is no mass to launch */
  if (si->feedback_data.mass_to_launch <= 0.f) return;

  /* No mass surrounding the star, no kick */
  if (si->feedback_data.wind_wt_sum <= 0.f) return;

#ifndef BLACK_HOLES_NONE
  /* If pj is being swallowed by a black hole, don't kick again */
  if (pj->black_holes_data.swallow_id > -1) return;
#endif

  /* If pj is being kicked by a star particle, don't kick again */
  if (pj->feedback_data.kick_id > -1) return;

  /* If pj is already a wind particle, don't kick again */
  if (pj->decoupled) return; 

  /* Get the gas mass. */
  const float mj = hydro_get_mass(pj);

  /* Get r. */
  const float r = sqrtf(r2);

  /* No kicks far away from the star */
  if (r >= KICK_RADIUS_OVER_H * hi) return;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;

  float wi;
  kernel_eval(ui, &wi);

  /* Bias towards the center of the kernel and to high SFR */
  const float wt = feedback_kernel_weight(pj, wi, ui);

  /* No kick if weight is zero */
  if (wt <= 0.f) return;

  /* Estimated number of particles to kick out of the kernel */
  float N_to_launch = si->feedback_data.mass_to_launch / mj;

  /* Make sure that stars do not kick too much mass out of the kernel */
  /* The rest of the mass will be kicked out later */
  if (N_to_launch > MAX_FRAC_OF_KERNEL_TO_LAUNCH * si->feedback_data.wind_ngb_mass / mj) {
    N_to_launch =  MAX_FRAC_OF_KERNEL_TO_LAUNCH * si->feedback_data.wind_ngb_mass / mj;
  }

  /* Apply redshift correction */
  N_to_launch *= si->feedback_data.eta_suppression_factor;

  /* Probability to swallow this particle */
  const float prob = N_to_launch * wt / si->feedback_data.wind_wt_sum;

#ifdef KIARA_DEBUG_CHECKS
  message("STAR_PROB: sid=%lld, gid=%lld, prob=%g, eta=%g, mlaunch=%g, "
          "m*=%g, N_to_launch=%g, mgas=%g, wt_sum=%g",
          si->id,
          pj->id,
          prob,
          si->feedback_data.mass_to_launch / si->mass_init,
          si->feedback_data.mass_to_launch,
          si->mass_init,
          N_to_launch,
          mj,
          si->feedback_data.wind_wt_sum);
#endif

  /* Draw a random number (Note mixing both IDs) */
  const float rand = random_unit_interval(si->id + pj->id, ti_current,
                                          random_number_stellar_feedback_1);

  /* We kick! */
  if (rand < prob) {
    pj->feedback_data.kick_id = si->id;
  }
}

/**
 * @brief Compile gas particles to be kicked by stellar feedback in this step,
 * These are the (up to) FEEDBACK_N_KICK_MAX nearest gas particles to the star.
 *
 * @param r2 Distance squared from star i to gas j.
 * @param dx[3] x,y,z distance from star i to gas j.
 * @param hi Smoothing length of star.
 * @param hj Smoothing length of gas.
 * @param si First (star) particle.
 * @param pj Second (gas) particle (not updated).
 * @param xpj Extra particle data.
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep2(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  struct spart *si, const struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const integertime_t ti_current) {

  /* Remove mass from the mass_to_launch reservoir */
  if (pj->feedback_data.kick_id == si->id) {
    si->feedback_data.mass_to_launch -= hydro_get_mass(pj);
    si->feedback_data.total_mass_kicked += hydro_get_mass(pj);

    /* Reservoir is comoving */
    const float v2 = 
        si->feedback_data.wind_velocity * si->feedback_data.wind_velocity;
    const double energy_phys = 0.5 * hydro_get_mass(pj) * v2 * cosmo->a2_inv;

    /* Remove energy used to kick particle from the SNII energy reservoir */
    si->feedback_data.physical_energy_reservoir -= energy_phys;

    /* Keep track of how many particles launched */
    si->feedback_data.N_launched += 1;

  }

}

/**
 * @brief Kick and sometimes heat gas particle near a star, 
 * if star has enough mass and energy for an ejection event.
 *
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 */
__attribute__((always_inline)) INLINE static void
feedback_kick_gas_around_star(
    const struct spart *si, struct part *pj, struct xpart *xpj,
    const struct cosmology *cosmo, 
    const struct feedback_props *fb_props, 
    const integertime_t ti_current) {

  if (pj->feedback_data.kick_id == si->id) {
    const float rand_for_kick_dir = 
        random_unit_interval(pj->id, ti_current, 
                             random_number_stellar_feedback);
    float kick_dir = 1.f;
    if (rand_for_kick_dir < 0.5f) {
      kick_dir = -1.f;
    }

    /* Compute velocity and KE of wind event.
    * Note that pj->v_full = a^2 * dx/dt, with x the comoving
    * coordinate. Therefore, a physical kick, dv, gets translated into a
    * code velocity kick, a * dv.
    */
    const float wind_velocity = 
        kick_dir * si->feedback_data.wind_velocity;
    const float wind_velocity_phys = 
        fabs(wind_velocity * cosmo->a_inv);

    /* Set kick direction as v x a */
    const float dir[3] = {pj->feedback_data.wind_direction[0],
                          pj->feedback_data.wind_direction[1],
                          pj->feedback_data.wind_direction[2]};
    const float norm = 
        sqrtf(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);

    /* No normalization, no wind (should basically never happen) */
    if (norm <= 0.f) {
      error("Normalization of wind direction is zero!\n(x, y, z) "
            "= (%g, %g, %g); vw=%g",
            dir[0], dir[1], dir[2], fabs(wind_velocity * cosmo->a_inv));
    }

    const float prefactor = wind_velocity / norm;

    /* Do the kicks by updating the particle velocity. */
    pj->v_full[0] += dir[0] * prefactor;
    pj->v_full[1] += dir[1] * prefactor;
    pj->v_full[2] += dir[2] * prefactor;

    /* DO WIND HEATING */
    float u_new = fb_props->cold_wind_internal_energy;
    if (fb_props->cold_wind_internal_energy < 
            fb_props->hot_wind_internal_energy) {
      float galaxy_stellar_mass =
            pj->galaxy_data.stellar_mass;
      if (galaxy_stellar_mass < fb_props->minimum_galaxy_stellar_mass) {
        galaxy_stellar_mass = fb_props->minimum_galaxy_stellar_mass;
      }
      const float galaxy_stellar_mass_Msun = galaxy_stellar_mass * 
      fb_props->mass_to_solar_mass;

      /* Based on Pandya et al 2022 FIRE results */
      float pandya_slope = 0.f;
      if (galaxy_stellar_mass_Msun > 3.16e10) {
        pandya_slope = -2.1f;
      } else {
        pandya_slope = -0.1f;
      }

      /* 0.2511886 = pow(10., -0.6) */
      const float f_warm = 
          0.2511886f * pow(galaxy_stellar_mass_Msun / 3.16e10f, pandya_slope);
      /* additional 10% removed for cold phase */
      const float hot_wind_fraction = max(0.f, 0.9f - f_warm);
      const float rand_for_hot = 
          random_unit_interval(pj->id, ti_current, 
              random_number_stellar_feedback_3);
      const float rand_for_spread = 
          random_unit_interval(pj->id, ti_current, 
              random_number_stellar_feedback);

      /* If selected, heat the particle */
      const float u_wind = 0.5 * wind_velocity_phys * wind_velocity_phys;
      if (rand_for_hot < hot_wind_fraction && 
              fb_props->hot_wind_internal_energy > u_wind) {
          u_new = (fb_props->hot_wind_internal_energy - u_wind) * 
                      (0.5 + rand_for_spread);
          u_new += hydro_get_drifted_physical_internal_energy(pj, cosmo);
      }
    }

    /* Set the wind particle internal energy */
    hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new);
    hydro_set_drifted_physical_internal_energy(pj, cosmo, NULL, u_new);

#ifdef FIREHOSE_DEBUG_CHECKS
    /* For firehose model, set initial radius of stream */
    if (si->feedback_data.firehose_radius_stream <= 0.f) {
      error("Firehose error: firehose_radius_stream <= 0. sid=%lld "
            "pid=%lld Rstream=%g",
            si->id,
            pj->id,
            si->feedback_data.firehose_radius_stream);
    }
#endif

    pj->chemistry_data.radius_stream = si->feedback_data.firehose_radius_stream;
    pj->chemistry_data.exchanged_mass = 0.f;

    /* FINISH UP FEEDBACK */
    /* Turn off any star formation in wind particle.
    * Record exp factor of when this particle was last ejected as -SFR. */
    pj->sf_data.SFR = -cosmo->a;

    /* Update the signal velocity of the particle based on the velocity kick,
      wind_velocity must be PHYSICAL passed into this function */
    hydro_set_v_sig_based_on_velocity_kick(pj, cosmo, wind_velocity_phys);

    /* Impose maximal viscosity */
    hydro_diffusive_feedback_reset(pj);

    /* Synchronize the particle on the timeline */
    timestep_sync_part(pj);

    /* Need time-step for decoupling */
    const integertime_t ti_step = get_integer_timestep(pj->time_bin);
    const integertime_t ti_begin =
      get_integer_time_begin(ti_current - 1, pj->time_bin);

    /* TODO: Requires always having with_cosmology! */
    const double dt = 
        cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);

    /* Mark to be decoupled */
    pj->to_be_decoupled = 1;
    pj->to_be_recoupled = 0;

    /* Decouple the particles from the hydrodynamics */
    pj->feedback_data.decoupling_delay_time =
        dt + fb_props->wind_decouple_time_factor *
             cosmology_get_time_since_big_bang(cosmo, cosmo->a);

    /* TODO: Move to chemistry module */
    pj->chemistry_data.diffusion_coefficient = 0.f;

    /* Take particle out of subgrid ISM mode */
    pj->cooling_data.subgrid_temp = 0.f;
    pj->cooling_data.subgrid_dens = hydro_get_physical_density(pj, cosmo);
    pj->cooling_data.subgrid_fcold = 0.f;

    pj->feedback_data.number_of_times_decoupled += 1;

    /* Kicked and handled */
    pj->feedback_data.kick_id = -1;

    /** Log the wind event.
     * z starid gasid dt M* vkick vkx vky vkz h x y z vx vy vz T rho v_sig tdec 
     * Ndec Z
     */
    const float length_convert = cosmo->a * fb_props->length_to_kpc;
    const float velocity_convert = cosmo->a_inv / fb_props->kms_to_internal;
    const float rho_convert = cosmo->a3_inv * fb_props->rho_to_n_cgs;
    const float u_convert =
        cosmo->a_factor_internal_energy / fb_props->temp_to_u_factor;

    printf("WIND_LOG %.5f %lld %g %g %g %g %g %lld %g %g %g %g %g %g "
           "%g %g %g %g %g "
           "%g %g %g %g %d %g\n",
            cosmo->z,
            si->id,
            si->feedback_data.mass_to_launch * 
                fb_props->mass_to_solar_mass,
            si->feedback_data.total_mass_kicked *
                fb_props->mass_to_solar_mass,
            si->feedback_data.total_mass_kicked / 
                si->mass,
            1.f/si->birth_scale_factor - 1.f,
            pj->galaxy_data.stellar_mass * 
                fb_props->mass_to_solar_mass,
            pj->id,
            fabs(wind_velocity) * velocity_convert,
            prefactor * dir[0] * velocity_convert,
            prefactor * dir[1] * velocity_convert,
            prefactor * dir[2] * velocity_convert,
            pj->h * length_convert, 
            pj->x[0] * length_convert,
            pj->x[1] * length_convert,
            pj->x[2] * length_convert,
            pj->v_full[0] * velocity_convert,
            pj->v_full[1] * velocity_convert,
            pj->v_full[2] * velocity_convert,
            pj->u * u_convert,
            pj->rho * rho_convert,
            pj->viscosity.v_sig * velocity_convert,
            pj->feedback_data.decoupling_delay_time * fb_props->time_to_Myr,
            pj->feedback_data.number_of_times_decoupled,
            pj->chemistry_data.metal_mass_fraction_total);

  }
}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 * Used for updating properties of gas particles neighbouring a star particle
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 */
__attribute__((always_inline)) INLINE static void
feedback_do_chemical_enrichment_of_gas_around_star(
    const float r2, const float dx[3], const float hi, const float hj,
    const struct spart *si, struct part *pj, struct xpart *xpj,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props, 
    const integertime_t ti_current) {

  /* Nothing to distribute */
  if (si->feedback_data.mass <= 0.f || 
      si->feedback_data.kernel_wt_sum <= 0.f) return;

  /* Gas particle density */
  const float rho_j = hydro_get_comoving_density(pj);
  if (rho_j <= 0.f) return;

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  const double current_mass = hydro_get_mass(pj);
  /* Compute weighting for distributing feedback quantities.
   * f = (mi * wi) / sum(mj * wj) */
  float Omega_frac = current_mass * wi / si->feedback_data.kernel_wt_sum;

  /* Never apply feedback if Omega_frac is bigger than or equal to unity */
  if (Omega_frac < 0.f || (Omega_frac > 1.f && ui < 1.f)) {
    warning(
        "Invalid fraction of material to distribute for star ID=%lld "
        "Omega_frac=%e count since last enrich=%d kernel_wt_sum=%g "
        "wi=%g rho_j=%g",
        si->id, Omega_frac, si->count_since_last_enrichment,
	      si->feedback_data.kernel_wt_sum, wi , rho_j);
    if (Omega_frac < 0.f || (Omega_frac > 1.01f && ui < 1.f)) {
      error("Omega_frac negative or too large! aborting");
    }
    
    Omega_frac = fmin(Omega_frac, 1.f);
  }

  /* Update particle mass */
  const double delta_mass = si->feedback_data.mass * Omega_frac;
  const double new_mass = current_mass + delta_mass;

  hydro_set_mass(pj, new_mass);

  /* Inverse of the new mass */
  const double new_mass_inv = 1. / new_mass;

  /* Update particle energy */
  const double injected_energy = si->feedback_data.energy * Omega_frac;

  /* Compute the current kinetic energy */
  const double current_v2 = pj->v_full[0] * pj->v_full[0] +
                            pj->v_full[1] * pj->v_full[1] +
                            pj->v_full[2] * pj->v_full[2];
  const double current_kinetic_energy_gas =
      0.5 * cosmo->a2_inv * current_mass * current_v2;

  /* Compute the current thermal energy */
  const double current_thermal_energy =
      current_mass * hydro_get_physical_internal_energy(pj, xpj, cosmo);

  /* Update velocity following change in gas mass */
  pj->v_full[0] *= current_mass * new_mass_inv;
  pj->v_full[1] *= current_mass * new_mass_inv;
  pj->v_full[2] *= current_mass * new_mass_inv;

  /* Update velocity following addition of mass with different momentum */
  pj->v_full[0] += delta_mass * new_mass_inv * si->v[0];
  pj->v_full[1] += delta_mass * new_mass_inv * si->v[1];
  pj->v_full[2] += delta_mass * new_mass_inv * si->v[2];

  /* Compute the new kinetic energy */
  const double new_v2 = pj->v_full[0] * pj->v_full[0] +
                        pj->v_full[1] * pj->v_full[1] +
                        pj->v_full[2] * pj->v_full[2];
  const double new_kinetic_energy_gas = 0.5 * cosmo->a2_inv * new_mass * new_v2;

  const double delta_KE = new_kinetic_energy_gas - current_kinetic_energy_gas;
  double new_thermal_energy = 
      current_thermal_energy + injected_energy - delta_KE;

  /* Following SPHENIX, don't decrease energy by more than 2x */
  new_thermal_energy = max(0.5f * current_thermal_energy, new_thermal_energy);

  /* Never go below the absolute minimum */
  const double min_u = hydro_props->minimal_internal_energy * new_mass;
  new_thermal_energy = max(new_thermal_energy, min_u);

  const double new_u = new_thermal_energy * new_mass_inv;

  hydro_set_physical_internal_energy(pj, xpj, cosmo, new_u);
  hydro_set_drifted_physical_internal_energy(pj, cosmo, /*pfloor=*/NULL,
                                             new_u);

  /* Recompute Z since we do not track all of the metals from Chem5 */
  pj->chemistry_data.metal_mass_fraction_total = 0.f;

  /* Update mass fraction of each tracked element  */
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    const double current_metal_mass =
        pj->chemistry_data.metal_mass_fraction[elem] * current_mass;
    const double delta_metal_mass =
        si->feedback_data.metal_mass[elem] * Omega_frac;
    const double new_metal_mass = current_metal_mass + delta_metal_mass;

    pj->chemistry_data.metal_mass_fraction[elem] =
        new_metal_mass * new_mass_inv;

    if (elem != chemistry_element_H && elem != chemistry_element_He) {
      pj->chemistry_data.metal_mass_fraction_total +=
          pj->chemistry_data.metal_mass_fraction[elem];
    }
  }

  /* Make sure that X + Y + Z = 1 */
  const float Y_He = 
      pj->chemistry_data.metal_mass_fraction[chemistry_element_He];
  const float Z = pj->chemistry_data.metal_mass_fraction_total;
  const float X_H = 1.f - Y_He - Z;

  if (X_H < 0.f || X_H > 1.f) {
    for (int elem = 0; elem < chemistry_element_count; elem++) {
      warning("\telem[%d] is %g",
              elem, pj->chemistry_data.metal_mass_fraction[elem]);
    }

    error("Hydrogen fraction exeeds unity or is negative for"
          " particle id=%lld due to stellar feedback.", pj->id);
  }

  pj->chemistry_data.metal_mass_fraction[chemistry_element_H] = X_H;

  /* Compute kernel-smoothed contribution to number of SNe going off 
   * this timestep */
  pj->feedback_data.SNe_ThisTimeStep += 
      si->feedback_data.SNe_ThisTimeStep * Omega_frac;
  pj->feedback_data.SNe_ThisTimeStep = 
      fmax(pj->feedback_data.SNe_ThisTimeStep, 0.);

  /* Spread dust ejecta to gas */
  for (int elem = chemistry_element_He; 
          elem < chemistry_element_count; elem++) {
    const double current_dust_mass =
        pj->cooling_data.dust_mass_fraction[elem] * pj->cooling_data.dust_mass;
    const double delta_dust_mass =
        si->feedback_data.delta_dust_mass[elem] * Omega_frac;

    /* at the moment this stores the mass (not mass frac) in each elem */
    pj->cooling_data.dust_mass_fraction[elem] =
        (current_dust_mass + delta_dust_mass);
  }

  /* Sum up each element to get total dust mass */
  pj->cooling_data.dust_mass = 0.;
  for (int elem = chemistry_element_He; 
          elem < chemistry_element_count; elem++) {
    pj->cooling_data.dust_mass += pj->cooling_data.dust_mass_fraction[elem];
  }

  if (pj->cooling_data.dust_mass > 0.) {
    const double dust_mass_inv = 1. / pj->cooling_data.dust_mass;

    /* Divide by new dust mass to get the fractions */
    for (int elem = chemistry_element_He; 
            elem < chemistry_element_count; elem++) {
      pj->cooling_data.dust_mass_fraction[elem] *= dust_mass_inv;
    }

    /* Check for inconsistency */
    if (pj->cooling_data.dust_mass > pj->mass) {
      for (int elem = chemistry_element_He; 
              elem < chemistry_element_count; elem++) {
        message("DUST EXCEEDS MASS elem=%d md=%g delta=%g \n",
                elem, 
                pj->cooling_data.dust_mass_fraction[elem] * 
                    pj->cooling_data.dust_mass,
                si->feedback_data.delta_dust_mass[elem] * Omega_frac);
      }

      error("DUST EXCEEDS MASS mgas=%g  mdust=%g\n",
            pj->mass, 
            pj->cooling_data.dust_mass);
    }
  }
  else {
    /* For some reason the dust mass is zero or negative */
    for (int elem = 0; elem < chemistry_element_count; elem++) {
      pj->cooling_data.dust_mass_fraction[elem] = 0.f;
    }
  }

}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 * Used for updating properties of gas particles neighbouring a star particle
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(
    const float r2, const float dx[3], const float hi, const float hj,
    const struct spart *si, struct part *pj, struct xpart *xpj,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props, 
    const integertime_t ti_current) {

  /* Ignore decoupled particles */
  if (pj->decoupled) return;

  /* Do chemical enrichment of gas, metals and dust from star */
  feedback_do_chemical_enrichment_of_gas_around_star(
    r2, dx, hi, hj, si, pj, xpj, cosmo, hydro_props,
    fb_props, ti_current);

  /* Do kinetic wind feedback */
  feedback_kick_gas_around_star(si, pj, xpj, cosmo, fb_props, ti_current);

#if COOLING_GRACKLE_MODE >= 2
  /* NOT USED: Compute G0 contribution from star to the gas particle in Habing units of 
   * 1.6e-3 erg/s/cm^2. Note that this value includes the 4*pi geometric factor 
  if (0) {
    const float length_to_physical_cm = cosmo->a * fb_props->length_to_kpc * 3.08567758e21f;
    // Compute a softened distance from star to gas particle 
    const double r2_in_cm = (r2 + 0.01*hi*hi) * length_to_physical_cm * length_to_physical_cm;
    const double r_in_cm = sqrt(r2_in_cm);
  
    // Compute self-shielding from H2, from Schauer et al. 2015 eq 8,9
    // H attenuation factor 
    const double NH_cgs = hydro_get_physical_density(pj, cosmo) * fb_props->rho_to_n_cgs * r_in_cm;
    const double xH = NH_cgs / 2.85e23;
    const double fH_shield = pow(1.f+xH,-1.62) * exp(-0.149*xH);
    // H2 attenuation factor
    const double NH2_cgs = pj->sf_data.H2_fraction * NH_cgs;
    const double DH2_cgs = 1.e-5 * sqrt(2.*1.38e-16*cooling_get_subgrid_temperature(pj, xpj) / 3.346e-24);
    const double xH2 = NH2_cgs / 8.465e13;
    const double fH2_shield = 0.9379/pow(1.f+xH2/DH2_cgs,1.879) + 0.03465/pow(1.f+xH2,0.473) * exp(-2.293e-4*sqrt(1+xH2));
    //message("G0 shield: r=%g xH2=%g xH=%g fH2=%g fH=%g\n",r_in_cm/3.086e21,xH2,xH,fH2_shield,fH_shield);
  
    if (si->feedback_data.lum_habing > -10.) {  
      pj->chemistry_data.G0 += fH2_shield * fH_shield * pow(10.,si->feedback_data.lum_habing) / (1.6e-3 * r2_in_cm);
    }
  }*/
#endif

}

#endif /* SWIFT_KIARA_FEEDBACK_IACT_H */
