/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_SIMBA_BLACK_HOLES_H
#define SWIFT_SIMBA_BLACK_HOLES_H

/* Local includes */
#include "black_holes_properties.h"
#include "black_holes_struct.h"
#include "star_formation.h"
#include "cooling.h"
#include "cosmology.h"
#include "dimension.h"
#include "exp10.h"
#include "gravity.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "physical_constants.h"
#include "random.h"

/* Standard includes */
#include <float.h>
#include <math.h>

/**
 * @brief Determine the mass loading factor for
 *        a given star particle based on its host galaxy.
 *
 * @param group_stellar_mass The stellar mass of the host galaxy
 * @param minimum_galaxy_stellar_mass Floor for stellar mass in eta calculation
 * @param FIRE_eta_normalization Normalization of eta at FIRE_eta_break
 * @param FIRE_eta_break M* at which eta(M*) slope changes
 * @param FIRE_eta_lower_slope Slope below break
 * @param FIRE_eta_upper_slope Slope above break
 */
__attribute__((always_inline)) INLINE static float feedback_mass_loading_factor(
    const double group_stellar_mass,
    const float minimum_galaxy_stellar_mass,
    const float FIRE_eta_normalization,
    const float FIRE_eta_break,
    const float FIRE_eta_lower_slope,
    const float FIRE_eta_upper_slope) {

  const float m_star =
      max(group_stellar_mass, minimum_galaxy_stellar_mass);

  float slope = FIRE_eta_lower_slope;
  if (m_star > FIRE_eta_break) {
    slope = FIRE_eta_upper_slope;
  }

  float eta =
      FIRE_eta_normalization * powf(m_star / FIRE_eta_break, slope);

  return eta;
}

/**
 * @brief Computes the time-step of a given black hole particle.
 *
 * @param bp Pointer to the s-particle data.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 */
__attribute__((always_inline)) INLINE static float black_holes_compute_timestep(
    const struct bpart* const bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo) {

  /* Allow for finer timestepping if necessary! */
  float dt_accr = FLT_MAX;
  /* Only need to the limit in the resolved feedback regime */
  if (bp->total_accretion_rate > 0.f && bp->subgrid_mass > bp->mass) {
    dt_accr = props->dt_accretion_factor * bp->mass / bp->total_accretion_rate;
  }

  if (dt_accr < props->time_step_min) {
    message(
        "Warning! BH_TIMESTEP_LOW: id=%lld (%g Myr) is below time_step_min (%g "
        "Myr).",
        bp->id, dt_accr * props->time_to_Myr,
        props->time_step_min * props->time_to_Myr);
  }

  return max(dt_accr, props->time_step_min);
}

/**
 * @brief Initialises the b-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param bp The particle to act upon
 * @param props The properties of the black holes model.
 */
__attribute__((always_inline)) INLINE static void black_holes_first_init_bpart(
    struct bpart* bp, const struct black_holes_props* props) {

  bp->time_bin = 0;
  if (props->use_subgrid_mass_from_ics == 0) {
    bp->subgrid_mass = bp->mass;
  } else if (props->with_subgrid_mass_check && bp->subgrid_mass <= 0) {
    error(
        "Black hole %lld has a subgrid mass of %f (internal units).\n"
        "If this is because the ICs do not contain a 'SubgridMass' data "
        "set, you should set the parameter "
        "'SIMBAAGN:use_subgrid_mass_from_ics' to 0 to initialize the "
        "black hole subgrid masses to the corresponding dynamical masses.\n"
        "If the subgrid mass is intentionally set to this value, you can "
        "disable this error by setting 'SIMBAAGN:with_subgrid_mass_check' "
        "to 0.",
        bp->id, bp->subgrid_mass);
  }
  bp->total_accreted_mass = 0.f;
  bp->accretion_disk_mass = 0.f; 
  bp->accretion_rate = 0.f;
  bp->total_accretion_rate = 0.f;
  bp->formation_time = -1.f;
  bp->cumulative_number_seeds = 1;
  bp->number_of_mergers = 0;
  bp->number_of_gas_swallows = 0;
  bp->number_of_direct_gas_swallows = 0;
  bp->number_of_repositions = 0;
  bp->number_of_reposition_attempts = 0;
  bp->number_of_time_steps = 0;
  bp->last_minor_merger_time = -1.;
  bp->last_major_merger_time = -1.;
  bp->swallowed_angular_momentum[0] = 0.f;
  bp->swallowed_angular_momentum[1] = 0.f;
  bp->swallowed_angular_momentum[2] = 0.f;
  bp->last_repos_vel = 0.f;
  bp->dt_heat = FLT_MAX;
  bp->dt_accr = FLT_MAX;
  bp->radiative_luminosity = 0.f;


}

/**
 * @brief Prepares a b-particle for its interactions
 *
 * @param bp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void black_holes_init_bpart(
    struct bpart* bp) {

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    bp->ids_ngbs_density[i] = -1;
  bp->num_ngb_density = 0;
#endif

  bp->density.wcount = 0.f;
  bp->density.wcount_dh = 0.f;
  bp->rho_gas = 0.f;
  bp->sound_speed_gas = 0.f;
  bp->internal_energy_gas = 0.f;
  bp->hot_gas_mass = 0.f;
  bp->cold_gas_mass = 0.f;
  bp->gas_SFR = 0.f;
  bp->hot_gas_internal_energy = 0.f;
  bp->rho_subgrid_gas = -1.f;
  bp->sound_speed_subgrid_gas = -1.f;
  bp->circular_velocity_gas[0] = 0.f;
  bp->circular_velocity_gas[1] = 0.f;
  bp->circular_velocity_gas[2] = 0.f;
  bp->angular_momentum_gas[0] = 0.f;
  bp->angular_momentum_gas[1] = 0.f;
  bp->angular_momentum_gas[2] = 0.f;
  bp->specific_angular_momentum_stars[0] = 0.f;
  bp->specific_angular_momentum_stars[1] = 0.f;
  bp->specific_angular_momentum_stars[2] = 0.f;
  bp->stellar_mass = 0.f;
  bp->radiative_luminosity = 0.f;
  bp->ngb_mass = 0.f;
  bp->num_ngbs = 0;
  bp->reposition.delta_x[0] = -FLT_MAX;
  bp->reposition.delta_x[1] = -FLT_MAX;
  bp->reposition.delta_x[2] = -FLT_MAX;
  bp->reposition.min_potential = FLT_MAX;
  bp->reposition.potential = FLT_MAX;
  bp->accretion_rate = 0.f; /* Optionally accumulated ngb-by-ngb */
  bp->total_accretion_rate = 0.f;
  bp->accretion_boost_factor = 0.f;
  bp->mass_at_start_of_step = bp->mass; /* bp->mass may grow in nibbling mode */
  bp->cold_disk_mass = 0.f;
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * The fields do not get predicted but we move the BH to its new position
 * if a new one was calculated in the repositioning loop.
 *
 * @param bp The particle
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void black_holes_predict_extra(
    struct bpart* restrict bp, float dt_drift) {

  /* Are we doing some repositioning? */
  if (bp->reposition.min_potential != FLT_MAX) {

#ifdef SWIFT_DEBUG_CHECKS
    if (bp->reposition.delta_x[0] == -FLT_MAX ||
        bp->reposition.delta_x[1] == -FLT_MAX ||
        bp->reposition.delta_x[2] == -FLT_MAX) {
      error("Something went wrong with the new repositioning position");
    }

    const double dx = bp->reposition.delta_x[0];
    const double dy = bp->reposition.delta_x[1];
    const double dz = bp->reposition.delta_x[2];
    const double d = sqrt(dx * dx + dy * dy + dz * dz);
    if (d > 1.01 * kernel_gamma * bp->h)
      error("Repositioning BH beyond the kernel support!");
#endif

    /* Move the black hole */
    bp->x[0] += bp->reposition.delta_x[0];
    bp->x[1] += bp->reposition.delta_x[1];
    bp->x[2] += bp->reposition.delta_x[2];

    /* Move its gravity properties as well */
    bp->gpart->x[0] += bp->reposition.delta_x[0];
    bp->gpart->x[1] += bp->reposition.delta_x[1];
    bp->gpart->x[2] += bp->reposition.delta_x[2];

    /* Store the delta position */
    bp->x_diff[0] -= bp->reposition.delta_x[0];
    bp->x_diff[1] -= bp->reposition.delta_x[1];
    bp->x_diff[2] -= bp->reposition.delta_x[2];

    /* Reset the reposition variables */
    bp->reposition.delta_x[0] = -FLT_MAX;
    bp->reposition.delta_x[1] = -FLT_MAX;
    bp->reposition.delta_x[2] = -FLT_MAX;
    bp->reposition.min_potential = FLT_MAX;

    /* Count the jump */
    bp->number_of_repositions++;
  }
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param bp The particle.
 */
__attribute__((always_inline)) INLINE static void
black_holes_reset_predicted_values(struct bpart* bp) {}

/**
 * @brief Kick the additional variables
 *
 * @param bp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void black_holes_kick_extra(
    struct bpart* bp, float dt) {}

/**
 * @brief Finishes the calculation of density on black holes
 *
 * @param bp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void black_holes_end_density(
    struct bpart* bp, const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = bp->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* --- Finish the calculation by inserting the missing h factors --- */
  bp->density.wcount *= h_inv_dim;
  bp->density.wcount_dh *= h_inv_dim_plus_one;
  bp->rho_gas *= h_inv_dim;
  float rho_inv = 1.f;
  if (bp->rho_gas > 0.f) rho_inv = 1.f / bp->rho_gas;
  /* All mass-weighted quantities are for the hot & cold gas */
  float m_hot_inv = 1.f;
  if (bp->hot_gas_mass > 0.f) m_hot_inv /= bp->hot_gas_mass;

  /* For the following, we also have to undo the mass smoothing */
  bp->sound_speed_gas *= h_inv_dim * rho_inv;
  bp->internal_energy_gas *= h_inv_dim * rho_inv;
  bp->hot_gas_internal_energy *= m_hot_inv;
  bp->circular_velocity_gas[0] *=
      h_inv_dim * rho_inv; /* h_inv_dim * rho_inv; */
  bp->circular_velocity_gas[1] *=
      h_inv_dim * rho_inv; /* h_inv_dim * rho_inv; */
  bp->circular_velocity_gas[2] *=
      h_inv_dim * rho_inv; /* h_inv_dim * rho_inv; */

  /* Calculate circular velocity at the smoothing radius from specific
   * angular momentum (extra h_inv) */
  bp->circular_velocity_gas[0] *= h_inv;
  bp->circular_velocity_gas[1] *= h_inv;
  bp->circular_velocity_gas[2] *= h_inv;
}

/**
 * @brief Sets all particle fields to sensible values when the #spart has 0
 * ngbs.
 *
 * @param bp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
black_holes_bpart_has_no_neighbours(struct bpart* bp,
                                    const struct cosmology* cosmo) {

#ifdef SWIFT_DEBUG_CHECKS
  warning(
      "BH particle with ID %lld treated as having no neighbours (h: %g, "
      "wcount: %g).",
      bp->id, bp->h, bp->density.wcount);
#endif

  /* Some smoothing length multiples. */
  const float h = bp->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  bp->density.wcount = kernel_root * h_inv_dim;
  bp->density.wcount_dh = 0.f;

  bp->internal_energy_gas = 0.f;
  bp->hot_gas_internal_energy = 0.f;
}

/**
 * @brief Return the current instantaneous accretion rate of the BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_accretion_rate(const struct bpart* bp) {
  return bp->accretion_rate;
}

/**
 * @brief Return the total accreted gas mass of this BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_accreted_mass(const struct bpart* bp) {
  return bp->total_accreted_mass;
}

/**
 * @brief Return the subgrid mass of this BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_subgrid_mass(const struct bpart* bp) {
  return bp->subgrid_mass;
}

/**
 *  * @brief Return the current bolometric luminosity of the BH.
 *   *
 *    * @param bp the #bpart.
 *     */
__attribute__((always_inline)) INLINE static double
black_holes_get_bolometric_luminosity(const struct bpart* bp,
                                      const struct phys_const* constants) {
  const double c = constants->const_speed_light_c;
  return bp->accretion_rate * 0.1 * c * c;
}

/**
 *  * @brief Return the current kinetic jet power of the BH.
 *   *
 *    * @param bp the #bpart.
 *     */
__attribute__((always_inline)) INLINE static double black_holes_get_jet_power(
    const struct bpart* bp, const struct phys_const* constants,
    const struct black_holes_props* props) {
  /* Now that we have v_kick we can determine the accretion fraction f_acc */
  const double psi = 
      props->wind_momentum_flux * props->epsilon_r * 
        constants->const_speed_light_c / bp->v_kick;
  return 0.5 * psi * bp->accretion_rate * bp->v_kick * bp->v_kick;
}

/**
 * @brief Update the properties of a black hole particles by swallowing
 * a gas particle.
 *
 * @param bp The #bpart to update.
 * @param p The #part that is swallowed.
 * @param xp The #xpart that is swallowed.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void black_holes_swallow_part(
    struct bpart* bp, const struct part* p, const struct xpart* xp,
    const struct cosmology* cosmo) {

  /* Get the current dynamical masses */
  const float gas_mass = hydro_get_mass(p);
  const float BH_mass = bp->mass;

  /* Increase the dynamical mass of the BH. */
  bp->mass += gas_mass;
  bp->gpart->mass += gas_mass;

  /* Physical velocity difference between the particles */
  const float dv[3] = {(bp->v[0] - p->v[0]) * cosmo->a_inv,
                       (bp->v[1] - p->v[1]) * cosmo->a_inv,
                       (bp->v[2] - p->v[2]) * cosmo->a_inv};

  /* Physical distance between the particles */
  const float dx[3] = {(bp->x[0] - p->x[0]) * cosmo->a,
                       (bp->x[1] - p->x[1]) * cosmo->a,
                       (bp->x[2] - p->x[2]) * cosmo->a};

  /* Collect the swallowed angular momentum */
  bp->swallowed_angular_momentum[0] +=
      gas_mass * (dx[1] * dv[2] - dx[2] * dv[1]);
  bp->swallowed_angular_momentum[1] +=
      gas_mass * (dx[2] * dv[0] - dx[0] * dv[2]);
  bp->swallowed_angular_momentum[2] +=
      gas_mass * (dx[0] * dv[1] - dx[1] * dv[0]);

  /* Update the BH momentum */
  const float BH_mom[3] = {BH_mass * bp->v[0] + gas_mass * p->v[0],
                           BH_mass * bp->v[1] + gas_mass * p->v[1],
                           BH_mass * bp->v[2] + gas_mass * p->v[2]};

  bp->v[0] = BH_mom[0] / bp->mass;
  bp->v[1] = BH_mom[1] / bp->mass;
  bp->v[2] = BH_mom[2] / bp->mass;
  bp->gpart->v_full[0] = bp->v[0];
  bp->gpart->v_full[1] = bp->v[1];
  bp->gpart->v_full[2] = bp->v[2];

  const float dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  message(
      "BH %lld swallowing gas particle %lld "
      "(Delta_v = [%f, %f, %f] U_V, "
      "Delta_x = [%f, %f, %f] U_L, "
      "Delta_v_rad = %f)",
      bp->id, p->id, -dv[0], -dv[1], -dv[2], -dx[0], -dx[1], -dx[2],
      (dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2]) / dr);

  /* Update the BH metal masses */
  struct chemistry_bpart_data* bp_chem = &bp->chemistry_data;
  const struct chemistry_part_data* p_chem = &p->chemistry_data;
  chemistry_add_part_to_bpart(bp_chem, p_chem, gas_mass);

  /* This BH swallowed a gas particle */
  bp->number_of_gas_swallows++;
  bp->number_of_direct_gas_swallows++;

  /* This BH lost a neighbour */
  bp->num_ngbs--;
  bp->ngb_mass -= gas_mass;
}

/**
 * @brief Update the properties of a black hole particles by swallowing
 * a BH particle.
 *
 * @param bpi The #bpart to update.
 * @param bpj The #bpart that is swallowed.
 * @param cosmo The current cosmological model.
 * @param time Time since the start of the simulation (non-cosmo mode).
 * @param with_cosmology Are we running with cosmology?
 * @param props The properties of the black hole scheme.
 */
__attribute__((always_inline)) INLINE static void black_holes_swallow_bpart(
    struct bpart* bpi, const struct bpart* bpj, const struct cosmology* cosmo,
    const double time, const int with_cosmology,
    const struct black_holes_props* props, const struct phys_const* constants) {

  /* Get the current dynamical masses */
  const float bpi_dyn_mass = bpi->mass;
  const float bpj_dyn_mass = bpj->mass;

  /* Is this merger ratio above the threshold for recording? */
  const double merger_ratio = bpj->subgrid_mass / bpi->subgrid_mass;
  if (merger_ratio > props->major_merger_threshold) {
    if (with_cosmology) {
      bpi->last_major_merger_scale_factor = cosmo->a;
    } else {
      bpi->last_major_merger_time = time;
    }
  } else if (merger_ratio > props->minor_merger_threshold) {
    if (with_cosmology) {
      bpi->last_minor_merger_scale_factor = cosmo->a;
    } else {
      bpi->last_minor_merger_time = time;
    }
  }

  /* Increase the masses of the BH. */
  bpi->mass += bpj->mass;
  bpi->gpart->mass += bpj->mass;
  bpi->subgrid_mass += bpj->subgrid_mass;

  /* Collect the swallowed angular momentum */
  bpi->swallowed_angular_momentum[0] += bpj->swallowed_angular_momentum[0];
  bpi->swallowed_angular_momentum[1] += bpj->swallowed_angular_momentum[1];
  bpi->swallowed_angular_momentum[2] += bpj->swallowed_angular_momentum[2];

  /* Update the BH momentum */
  const float BH_mom[3] = {bpi_dyn_mass * bpi->v[0] + bpj_dyn_mass * bpj->v[0],
                           bpi_dyn_mass * bpi->v[1] + bpj_dyn_mass * bpj->v[1],
                           bpi_dyn_mass * bpi->v[2] + bpj_dyn_mass * bpj->v[2]};

  bpi->v[0] = BH_mom[0] / bpi->mass;
  bpi->v[1] = BH_mom[1] / bpi->mass;
  bpi->v[2] = BH_mom[2] / bpi->mass;
  bpi->gpart->v_full[0] = bpi->v[0];
  bpi->gpart->v_full[1] = bpi->v[1];
  bpi->gpart->v_full[2] = bpi->v[2];

  /* Update the BH metal masses */
  struct chemistry_bpart_data* bpi_chem = &bpi->chemistry_data;
  const struct chemistry_bpart_data* bpj_chem = &bpj->chemistry_data;
  chemistry_add_bpart_to_bpart(bpi_chem, bpj_chem);

  /* Add up all the BH seeds */
  bpi->cumulative_number_seeds += bpj->cumulative_number_seeds;

  /* Add up all the gas particles we swallowed */
  bpi->number_of_gas_swallows += bpj->number_of_gas_swallows;

  /* We had another merger */
  bpi->number_of_mergers++;
}

/**
 * @brief Compute the accretion rate of the black hole and all the quantities
 * required for the feedback loop.
 *
 * @param bp The black hole particle.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 * @param cosmo The cosmological model.
 * @param cooling Properties of the cooling model.
 * @param floor_props Properties of the entropy fllor.
 * @param time Time since the start of the simulation (non-cosmo mode).
 * @param with_cosmology Are we running with cosmology?
 * @param dt The time-step size (in physical internal units).
 * @param ti_begin The time at which the step begun (ti_current).
 */
__attribute__((always_inline)) INLINE static void black_holes_prepare_feedback(
    struct bpart* restrict bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo,
    const struct cooling_function_data* cooling,
    const struct entropy_floor_properties* floor_props, const double time,
    const int with_cosmology, const double dt, const integertime_t ti_begin) {

  /* Record that the black hole has another active time step */
  bp->number_of_time_steps++;

  if (dt == 0. || bp->rho_gas == 0.) return;

  /* A black hole should never accrete/feedback if it is not in a galaxy */
  if (bp->galaxy_data.stellar_mass <= 0.f) return;

  /* Gather some physical constants (all in internal units) */
  const double G = constants->const_newton_G;
  const double c = constants->const_speed_light_c;
  const double proton_mass = constants->const_proton_mass;
  const double sigma_Thomson = constants->const_thomson_cross_section;

  /* Gather the parameters of the model */
  const double f_Edd_limit = props->f_Edd;
  const double epsilon_r = props->epsilon_r;

  /* (Subgrid) mass of the BH (internal units) */
  const double BH_mass = bp->subgrid_mass;

  /* Convert the quantities we gathered to physical frame (internal units). */

  const float h_inv = 1.f / (kernel_gamma * bp->h);
  const double volume_inv = (3. / (4. * M_PI)) * pow_dimension(h_inv);
  const double gas_rho = bp->hot_gas_mass * volume_inv;
  const double gas_rho_phys = gas_rho * cosmo->a3_inv;

  double gas_c_phys = 0.;
  if (bp->hot_gas_internal_energy > 0.) {
    gas_c_phys = gas_soundspeed_from_internal_energy(
                     gas_rho, bp->hot_gas_internal_energy) *
                 cosmo->a_factor_sound_speed;
  }
  double gas_c_phys2 = gas_c_phys * gas_c_phys;

  /* We can now compute the Bondi accretion rate (internal units)
   * D. Rennehan: In Simba, we only consider the hot gas within
   * the kernel for the Bondi rate, and the cold gas using the
   * torque accretion estimator.
   */
  double Bondi_rate = 0.;

  const double denominator2 = gas_c_phys2;
#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure that the denominator is strictly positive */
  if (denominator2 <= 0)
    warning(
        "Invalid denominator for black hole particle %lld in Bondi rate "
        "calculation.",
        bp->id);
#endif

  if (denominator2 > 0.) {
    const double denominator_inv = 1. / sqrt(denominator2);
    Bondi_rate = 4. * M_PI * G * G * BH_mass * BH_mass * gas_rho_phys *
                 denominator_inv * denominator_inv * denominator_inv;
  }

  /* Hot gas can only be boosted by an alpha factor, no density dependence */
  if (props->with_boost_factor) {
    Bondi_rate *= props->boost_alpha;
    bp->accretion_boost_factor = props->boost_alpha;
  } else {
    bp->accretion_boost_factor = 1.f;
  }

  /* Compute the Eddington rate (internal units) */
  const double Eddington_rate =
      4. * M_PI * G * BH_mass * proton_mass / (epsilon_r * c * sigma_Thomson);

  /* The accretion rate estimators give Mdot,inflow  (Mdot,BH = f_acc *
   * Mdot,inflow) */
  double accr_rate = props->f_accretion * Bondi_rate;

  /* Limit Bondi to Eddington_rate, as well as the Eddington rate for a specified,
   * custom Msun BH */
  const double Eddington_rate_custom_mass =
      Eddington_rate * (props->bondi_rate_limiting_bh_mass / BH_mass);

  accr_rate = min3(accr_rate, Eddington_rate, Eddington_rate_custom_mass);

#ifdef SIMBA_DEBUG_CHECKS
  double bondi_accr_rate = accr_rate;
#endif

  /* Now we compute the torque-limited accretion model */
  float torque_accr_rate = 0.f;
  const float disk_gas_mass = bp->cold_gas_mass;
  const float corot_gas_mass = 
      fmax(bp->cold_gas_mass - 2. * (bp->cold_gas_mass - bp->cold_disk_mass), 
           0.f); 
  const float galaxy_gas_mass = bp->galaxy_data.gas_mass;

  /* corrects from gas density to total density; set this to max value allowed */
  float f_corr_stellar = 10.f;
  if (galaxy_gas_mass > 0) {
    f_corr_stellar = 
        min(1. + bp->galaxy_data.stellar_mass / galaxy_gas_mass, 
            f_corr_stellar);
  }

  const float rho_bh = bp->subgrid_mass * volume_inv;
  /* Correct gas density for subgrid */
  float rho_subgrid_factor = 1.f;
  if (bp->rho_subgrid_gas > 0.f && bp->rho_gas > 0.f) {
    rho_subgrid_factor = bp->rho_subgrid_gas / bp->rho_gas;
  }

  /* Compute dynamical time */
  const float rho_kernel = 
      f_corr_stellar * bp->rho_gas * rho_subgrid_factor + rho_bh;
  const float tdyn_inv = 
      sqrtf(G * rho_kernel * cosmo->a3_inv);

  if (props->torque_accretion_method == 2) {
    /* Let's compute the torque-limited accretion rate from the cold gas as in Simba  */
    if (corot_gas_mass > 0.f &&  bp->cold_gas_mass > 0.f) {
      const float m_disk = bp->cold_gas_mass * f_corr_stellar;
      const float f_disk = corot_gas_mass / bp->cold_gas_mass; 
  
      const float r0 = bp->h * cosmo->a * (props->length_to_parsec / 100.f);
  
      const float alpha = 5.f;
      const float mass_to_1e9solar = props->mass_to_solar_mass / 1.0e9f;
      const float mass_to_1e8solar = props->mass_to_solar_mass / 1.0e8f;
  
      const float f0 =
          0.31f * f_disk * f_disk * 
          pow(m_disk * mass_to_1e9solar, -1.f / 3.f);
      const float f_gas = corot_gas_mass / m_disk;
  
      torque_accr_rate = props->torque_accretion_norm * alpha *
                          corot_gas_mass * mass_to_1e9solar *
                          powf(f_disk, 5.f / 2.f) *
                          powf(bp->subgrid_mass * mass_to_1e8solar, 1.f / 6.f) *
                          powf(r0, -3.f / 2.f) / (1 + f0 / f_gas);

      /* correct for subgrid density tdyn */
      torque_accr_rate *= sqrtf(rho_subgrid_factor);
      torque_accr_rate *=
          props->f_accretion * (props->time_to_yr / props->mass_to_solar_mass);
  
    } 
  }
  else if (props->torque_accretion_method == 1 || 
              props->torque_accretion_method == 0) {
    /* Here the accretion rate is only based on Mgas / tdyn.
     * We do not use the DM mass to compute tdyn since it probably
     * doesn't contribute much near the core of the system. We also
     * assume that the gas fraction is constant in the galaxy in
     * order to compute Mstar within the kernel of the black hole.
     * Therefore, Mdot = Mgas / tdyn = Mgas / sqrt(3pi/(32 G rho))
     * and rho = (Mgas + Mstar + Mdm) / (4pi h^3 / 3) where
     * Mstar = Mgas / fgas, Mdm = 0. Therefore,
     * rho = 3 * ((1 + fgas) / fgas) * Mgas / (4 * pi * h^3)
     * and
     * Mdot = Mgas * sqrt(32 * G * 3 * ((1 + fgas) / fgas) * Mgas)) /
     *    sqrt(3 * pi * 4 * pi * h^3)
     *      = sqrt(96 * G * ((1 + fgas) / fgas) * Mgas^3) /
     *    sqrt(12 * pi^2 * h^3)
     *      = (1 / pi) * sqrt(8 * G * ((1 + fgas) / fgas) * (Mgas / h)^3))
     */
    bp->stellar_mass = (bp->cold_gas_mass + bp->hot_gas_mass) * f_corr_stellar;
  
    /* Disk mass is total mass minus twice the counter-rotating part */
    if (props->torque_accretion_method == 1) {
      torque_accr_rate = props->torque_accretion_norm * corot_gas_mass * tdyn_inv;
    }
    /* Torque rate is based on disk mass falling in on a dynamical time */
    if (props->torque_accretion_method == 0) {
      torque_accr_rate = props->torque_accretion_norm * disk_gas_mass * tdyn_inv;
    }

    torque_accr_rate *= props->f_accretion;
  }

  /* Do suppression of BH growth */
  /* Mode 1 is based on the Hopkins+22 model for how SF suppresses BH accretion */
  double f_suppress = 1.f;
  if (props->suppress_growth == BH_suppress_Hopkins22) {
    const float r0 = bp->h * cosmo->a * props->length_to_parsec;
    const float sigma_eff = f_corr_stellar * bp->ngb_mass *
                            props->mass_to_solar_mass /
                            (M_PI * r0 * r0);

    f_suppress = sigma_eff / (sigma_eff + props->sigma_crit_resolution_factor *
                                     props->sigma_crit_Msun_pc2);
  }

  /* Mode 2 or 6 is exponential suppression like Simba-C, 
   * 2 uses all cold gas and 6 accretes only SF gas */
  if (props->suppress_growth == BH_suppress_ExponentialOnTorque || 
          props->suppress_growth == BH_suppress_ExponentialOnSFGas) {
    float m_suppress = fabs(props->characteristic_suppression_mass);
    if (props->characteristic_suppression_mass < 0) {
      m_suppress *= cosmo->a;
    }

    f_suppress = 
        1.f - exp(-bp->subgrid_mass * props->mass_to_solar_mass / m_suppress);
  }

  /* Mode 4 or 5 is based on amount of gas blown out in FIRE SF-driven winds 
   * during infall; 5 accretes only SF gas */
  if (props->suppress_growth == BH_suppress_OutflowsOnAllGas || 
          props->suppress_growth == BH_suppress_OutflowsOnSFGas || 
          props->suppress_growth == BH_suppress_OutflowsOrExponential) {
    /* compute mass loading factor from SF feedback, should be same as used 
     * in feedback_mass_loading_factor() */
    const float eta = 
        feedback_mass_loading_factor(bp->galaxy_data.stellar_mass, 
                                     props->minimum_galaxy_stellar_mass, 
                                     props->FIRE_eta_normalization, 
                                     props->FIRE_eta_break, 
                                     props->FIRE_eta_lower_slope, 
                                     props->FIRE_eta_upper_slope);

    /* Get star formation efficiency */
    const float sf_eff = props->star_formation_efficiency;

    f_suppress = expf(-eta * sf_eff);

    /* Mode 7 suppresses accretion by factor accounting for mass 
     * lost in outflow over dynamical time */
    if (props->suppress_growth == BH_suppress_OutflowsOrExponential) {
      float m_suppress = fabs(props->characteristic_suppression_mass);
      if (props->characteristic_suppression_mass < 0) {
        m_suppress *= cosmo->a;
      }

      f_suppress = 
        fmin(f_suppress, 
             1. - exp(-bp->subgrid_mass * props->mass_to_solar_mass / m_suppress));
    }
  }

  /* For these modes, only torque accretion is suppressed */
  torque_accr_rate *= f_suppress;

  if (isnan(torque_accr_rate) || torque_accr_rate < 0.) {
    error("torque_accr_rate is incorrect!\n"
          "f_corr_stellar = %g\n"
          "cold_disk_mass = %g Msun\n"
          "rho_gas = %g cm^-3\n"
          "torque_accretion_norm = %g\n"
          "f_accretion = %g\n"
          "torque_accretion_rate = %g\n",
          f_corr_stellar, 
          bp->cold_disk_mass * props->mass_to_solar_mass,
          gas_rho_phys * props->rho_to_n_cgs,
          props->torque_accretion_norm,
          props->f_accretion, torque_accr_rate);
  }

  /* Total accretion rate is Bondi + torque */
  accr_rate += torque_accr_rate;

  /* In Mode 3, suppression applies to all accretion not just torque */
  if (props->suppress_growth == BH_suppress_ExponentialOnTotal || 
        props->suppress_growth == BH_suppress_ExpOutflowsOnTotal) {
    float m_suppress = fabs(props->characteristic_suppression_mass);
    if (props->characteristic_suppression_mass < 0.f) {
      m_suppress *= cosmo->a;
    }

    f_suppress = 
        1. - exp(-bp->subgrid_mass * props->mass_to_solar_mass / m_suppress);

    /* Mode 8 additionally limits by exponential suppression applied 
     * to total accretion */
    if (props->suppress_growth == BH_suppress_ExpOutflowsOnTotal) {
      /* compute mass loading factor from SF feedback, should be same as used in 
       * feedback_mass_loading_factor() */
      const float eta = 
          feedback_mass_loading_factor(bp->galaxy_data.stellar_mass, 
                                       props->minimum_galaxy_stellar_mass, 
                                       props->FIRE_eta_normalization, 
                                       props->FIRE_eta_break, 
                                       props->FIRE_eta_lower_slope, 
                                       props->FIRE_eta_upper_slope);

      /* Get star formation efficiency */
      const float sf_eff = props->star_formation_efficiency;

      /* suppress accretion by factor accounting for mass lost in outflow 
       * over dynamical time */
      f_suppress = fmin(f_suppress, exp(-eta * sf_eff));
    }

    accr_rate *= f_suppress;
  }

  /* Limit overall accretion rate to some factor time the Eddington rate */
  accr_rate = min(accr_rate, f_Edd_limit * Eddington_rate);
  bp->accretion_rate = accr_rate;
  bp->total_accretion_rate = accr_rate;

  /* No accretion, so no feedback or anything else */
  if (accr_rate <= 0.) return;

  /* Factor in the radiative efficiency */
  const double mass_rate = (1. - epsilon_r) * accr_rate;
  const double luminosity = epsilon_r * accr_rate * c * c;

  /* This is used for X-ray feedback later */
  bp->radiative_luminosity = luminosity;

  /* Integrate mass forward in time */
  double delta_mass = mass_rate * dt;

  /* Compute infall times to BH at this redshift */
  const double E_inv = cosmo->H0 / cosmo->H;
  double t_infall = props->accretion_disk_timescale * E_inv;
  if (props->accretion_disk_timescale <= 10.) {
    /* assume a factor times the dynamical time */
    t_infall = min(props->accretion_disk_timescale / tdyn_inv, 100.f * E_inv);
  }

  /* If desired we put mass into accretion disk which feeds 
   * BH on some frac of tdyn 
   */
  if (t_infall > 0.) {
    /* Add accreted mass into a reservoir representing BH accretion disk */
    bp->accretion_disk_mass += delta_mass;

    /* Compute mass that will actually go into BH */
    delta_mass = bp->accretion_disk_mass * (1.f - exp(-dt / t_infall));

    /* This mass gets removed from the accretion disk */
    bp->accretion_disk_mass -= delta_mass;

    /* Recompute accretion rate based on the reservoir change */
    bp->accretion_rate = delta_mass / ((1. - epsilon_r) * dt);
  }

  bp->subgrid_mass += delta_mass;
  bp->total_accreted_mass += delta_mass;

  /* Calculate Eddington ratio, using true accretion rate */
  bp->eddington_fraction = bp->accretion_rate / Eddington_rate;

  if (props->use_nibbling && bp->subgrid_mass < bp->mass) {
    /* In this case, the BH is still accreting from its (assumed) subgrid gas
     * mass reservoir left over when it was formed. There is some loss in this
     * due to radiative losses, so we must decrease the particle mass
     * in proportion to its current accretion rate. We do not account for this
     * in the swallowing approach, however. */
    bp->mass -= epsilon_r * accr_rate * dt;
    if (bp->mass < 0)
      error("Black hole %lld reached negative mass (%g). Trouble ahead...",
            bp->id, bp->mass);
  }

  /* Compute the properties that don't change over the timestep,
   * i.e. feedback v_kick, f_acc, etc.
   */
  const float subgrid_mass_Msun = bp->subgrid_mass * props->mass_to_solar_mass;
  const float bondi_fraction = 1.f - (torque_accr_rate / accr_rate);

  /* Compute non-jet kick velocity */
  float v_kick = 0.f;
  if (bp->subgrid_mass > props->subgrid_seed_mass) {
    v_kick = 500.f + (500.f / 3.f) * (log10f(subgrid_mass_Msun) - 6.f);
    v_kick *= props->kms_to_internal;

    if (v_kick < 0.f) v_kick = 0.f;
  }

  /* Decide whether it is in jet mode */
  const int max_jet_mode = 
      (bondi_fraction > props->bondi_fraction_thresh_always_jet
	      || subgrid_mass_Msun > props->mass_thresh_always_jet
        || bp->accretion_rate > props->accr_rate_thresh_always_jet);

  const int jet_mode = 
      (bp->eddington_fraction < props->eddington_fraction_lower_boundary 
		    || max_jet_mode);

  /* Here we compute the jet mode feedback */
  if (v_kick > 0.f) {
    if (bp->eddington_fraction < 1.e-10f) {
      bp->eddington_fraction = 1.e-10f;
    }

    if (jet_mode) {
      /* Generate a spread in jet mass threshold among BH particles, based on id */
      const float mass_min = props->jet_mass_min_Msun;
      const float mass_max = 
          (props->jet_mass_min_Msun + props->jet_mass_spread_Msun);
      const float jet_mass_thresh_Msun = 
          mass_min + 0.01f * (float)(bp->id % 100) * (mass_max - mass_min);

      /* If the BH is massive enough, we compute the jet mode boost to v_kick */
      if (subgrid_mass_Msun > jet_mass_thresh_Msun) {
        /* jet_velocity is already in internal velocity units here
         * Determine max velocity of jet for this BH; max at 
         * max_multiplier*vjet @ MBH=1e8 */
        float jet_vmax = 
            props->jet_velocity * 
                powf(subgrid_mass_Msun * 1.e-8, 
                     props->jet_velocity_scaling_with_mass);
	      jet_vmax = max(jet_vmax, props->jet_velocity);
	      jet_vmax = 
            min(jet_vmax, 
                props->jet_velocity_max_multiplier * props->jet_velocity);

        /* Add some spread around the maximum velocity */
        const double random_number =
            random_unit_interval(bp->id, ti_begin, random_number_BH_feedback);
        jet_vmax *= props->jet_velocity_spread_alpha 
                    + props->jet_velocity_spread_beta * random_number;

        if (max_jet_mode) {
          v_kick += jet_vmax;
        }
      	else {
          const float scaled_jet_vel = 
                props->jet_velocity 
                * log10f(props->eddington_fraction_lower_boundary / 
                          bp->eddington_fraction);
          v_kick += min(scaled_jet_vel, jet_vmax);
        }
      }
    }

    /* Set v_kick with the jet velocity */
    bp->v_kick = v_kick;

    /* Now that we have v_kick we can determine the accretion fraction f_acc */
    const float momentum_scaling = epsilon_r * c / bp->v_kick;
    bp->f_accretion =
        1.f / (1.f + props->wind_momentum_flux * momentum_scaling);
  } 
  else {
    bp->f_accretion = props->f_accretion;
  }

#ifdef SIMBA_DEBUG_CHECKS
  double f_gas = bp->galaxy_data.gas_mass / 
      (bp->galaxy_data.stellar_mass + bp->galaxy_data.gas_mass);

  if (bp->subgrid_mass * props->mass_to_solar_mass > 1.e6) {
    message("BH_ACC: z=%g bid=%lld ms=%g mbh=%g ssfr=%g sfr=%g state=%d "
            "torque=%g bondi=%g fEdd=%g facc=%g fsupp=%g mcold=%g mhot=%g "
            "mdisk=%g tin=%g vkick=%g dmass=%g radeff=%g mres=%g fsub=%g fgas=%g",
            cosmo->z, 
            bp->id,
            bp->galaxy_data.stellar_mass * props->mass_to_solar_mass,
            bp->subgrid_mass * props->mass_to_solar_mass,
            bp->galaxy_data.ssfr / props->time_to_yr,
            bp->galaxy_data.ssfr * props->mass_to_solar_mass / props->time_to_yr,
            (bp->eddington_fraction > props->eddington_fraction_lower_boundary) 
                ? 1 : 0,
            torque_accr_rate * props->mass_to_solar_mass / props->time_to_yr,
            bondi_accr_rate * props->mass_to_solar_mass / props->time_to_yr,
            bp->eddington_fraction,
            bp->f_accretion,
            f_suppress,
            bp->cold_gas_mass * props->mass_to_solar_mass,
            bp->hot_gas_mass * props->mass_to_solar_mass,
            corot_gas_mass * props->mass_to_solar_mass,
            t_infall * props->time_to_Myr,
            bp->v_kick / props->kms_to_internal,
            delta_mass, 
            props->epsilon_r, 
            bp->accretion_disk_mass, 
            rho_subgrid_factor, 
            f_gas);
  }
#endif

  /**
   * 0 = scale factor
   * 1 = id
   * 2 = dynamical mass
   * 3 = subgrid (actual) BH mass
   * 4 = accretion rate
   * 5 = Bondi rate
   * 6 = torque rate
   * 7 = time step
   * 8 = accretion region mean density
   * 9 = hot gas energy
   * 10 = accretion region SFR
   * 11 = accretion region gas mass
   * 12 = kernel hot gas mass
   * 13 = kernel stellar mass
   * 14 = kernel kinematic disk mass
   * 15 = cumulative accretion disk mass
   * 16 = dynamical time
   * 17 = host galaxy stellar mass
   * 18 = host galaxy sSFR
   * 19 = size of the accretion region
   * 20 = x
   * 21 = y
   * 22 = z
   * 23 = vx
   * 24 = vy
   * 25 = vz
   * 26 = Lx kernel gas
   * 27 = Ly kernel gas
   * 28 = Lz kernel gas
   * 29 = Lx kernel stars
   * 30 = Ly kernel stars
   * 31 = Lz kernel stars
   * 32 = f_Edd
   * 33 = Mass accretion this timestep
   * 34 = f_Bondi
   * 35 = v_kick
   */
  printf("BH_DETAILS "
          "%2.12f %lld "
          " %g %g %g %g %g %g %g "
          " %g %g %g %g "
          " %g %g %g %g %g %g %g "
          " %2.10f %2.10f %2.10f "
          " %2.7f %2.7f %2.7f "
          " %g %g %g  %g %g %g"
          " %g %g %g %f\n",
          cosmo->a, 
          bp->id, 
          bp->mass * props->mass_to_solar_mass,
          bp->subgrid_mass * props->mass_to_solar_mass,
          bp->accretion_rate * props->mass_to_solar_mass / props->time_to_yr,
          Bondi_rate * props->mass_to_solar_mass / props->time_to_yr,
          torque_accr_rate * props->mass_to_solar_mass / props->time_to_yr,
          dt * props->time_to_Myr,
          (bp->rho_gas * cosmo->a3_inv) * props->rho_to_n_cgs,
          bp->hot_gas_internal_energy * cosmo->a_factor_internal_energy * 
              props->conv_factor_specific_energy_to_cgs,
          /* index 10 below */
          bp->gas_SFR, 
          (bp->hot_gas_mass + bp->cold_gas_mass) * props->mass_to_solar_mass,
          bp->hot_gas_mass * props->mass_to_solar_mass,
          bp->stellar_mass * props->mass_to_solar_mass, 
          disk_gas_mass * props->mass_to_solar_mass,
          bp->accretion_disk_mass * props->mass_to_solar_mass,
          1. / tdyn_inv * props->time_to_Myr,
          bp->galaxy_data.stellar_mass * props->mass_to_solar_mass, 
          bp->galaxy_data.ssfr / props->time_to_yr,
          bp->h * props->length_to_parsec / 1.0e3f, 
          /* index 20 below */
          bp->x[0] * cosmo->a * props->length_to_parsec / 1.0e3f, 
          bp->x[1] * cosmo->a * props->length_to_parsec / 1.0e3f,
          bp->x[2] * cosmo->a * props->length_to_parsec / 1.0e3f,
          bp->v[0] * cosmo->a_inv / props->kms_to_internal,
          bp->v[1] * cosmo->a_inv / props->kms_to_internal,
          bp->v[2] * cosmo->a_inv / props->kms_to_internal,
          bp->angular_momentum_gas[0], 
          bp->angular_momentum_gas[1],
          bp->angular_momentum_gas[2], 
          bp->specific_angular_momentum_stars[0], 
          /* index 30 below */
          bp->specific_angular_momentum_stars[1],
          bp->specific_angular_momentum_stars[2], 
          bp->eddington_fraction, 
          delta_mass * props->mass_to_solar_mass, 
          bondi_fraction, 
          v_kick / props->kms_to_internal);
}

/**
 * @brief Computes the (maximal) repositioning speed for a black hole.
 *f_Edd
 * Calculated as upsilon * (m_BH / m_ref) ^ beta_m * (n_H_BH / n_ref) ^ beta_n
 * where m_BH = BH subgrid mass, n_H_BH = physical gas density around BH
 * and upsilon, m_ref, beta_m, n_ref, and beta_n are parameters.
 *
 * @param bp The #bpart.
 * @param props The properties of the black hole model.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_repositioning_speed(const struct bpart* restrict bp,
                                    const struct black_holes_props* props,
                                    const struct cosmology* cosmo) {

  const double n_gas_phys = bp->rho_gas * cosmo->a3_inv * props->rho_to_n_cgs;
  const double v_repos =
      props->reposition_coefficient_upsilon *
      pow(bp->subgrid_mass / props->reposition_reference_mass,
          props->reposition_exponent_mass) *
      pow(n_gas_phys / props->reposition_reference_n_H,
          props->reposition_exponent_n_H);

  /* Make sure the repositioning is not back-firing... */
  if (v_repos < 0)
    error(
        "BH %lld wants to reposition at negative speed (%g U_V). Do you "
        "think you are being funny? No-one is laughing.",
        bp->id, v_repos);

  return v_repos;
}

/**
 * @brief Finish the calculation of the new BH position.
 *
 * Here, we check that the BH should indeed be moved in the next drift.
 *
 * @param bp The black hole particle.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 * @param cosmo The cosmological model.
 * @param dt The black hole particle's time step.
 * @param ti_begin The time at the start of the temp
 */
__attribute__((always_inline)) INLINE static void black_holes_end_reposition(
    struct bpart* restrict bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo,
    const double dt, const integertime_t ti_begin) {

  /* First check: did we find any eligible neighbour particle to jump to? */
  if (bp->reposition.min_potential != FLT_MAX) {

    /* Record that we have a (possible) repositioning situation */
    bp->number_of_reposition_attempts++;

    /* Is the potential lower (i.e. the BH is at the bottom already)
     * OR is the BH massive enough that we don't reposition? */
    const float potential = gravity_get_comoving_potential(bp->gpart);
    if (potential < bp->reposition.min_potential ||
        bp->subgrid_mass > props->max_reposition_mass) {

      /* No need to reposition */
      bp->reposition.min_potential = FLT_MAX;
      bp->reposition.delta_x[0] = -FLT_MAX;
      bp->reposition.delta_x[1] = -FLT_MAX;
      bp->reposition.delta_x[2] = -FLT_MAX;

    } else if (props->set_reposition_speed) {

      /* If we are re-positioning, move the BH a fraction of delta_x, so
       * that we have a well-defined re-positioning velocity (repos_vel
       * cannot be negative). */
      double repos_vel = black_holes_get_repositioning_speed(bp, props, cosmo);

      /* Convert target reposition velocity to a fractional reposition
       * along reposition.delta_x */
      const double dx = bp->reposition.delta_x[0];
      const double dy = bp->reposition.delta_x[1];
      const double dz = bp->reposition.delta_x[2];
      const double d = sqrt(dx * dx + dy * dy + dz * dz);

      /* Exclude the pathological case of repositioning by zero distance */
      if (d > 0) {
        double repos_frac = repos_vel * dt / d;

        /* We should never get negative repositioning fractions... */
        if (repos_frac < 0)
          error("Wanting to reposition by negative fraction (%g)?", repos_frac);

        /* ... but fractions > 1 can occur if the target velocity is high.
         * We do not want this, because it could lead to overshooting the
         * actual potential minimum. */
        if (repos_frac > 1) {
          repos_frac = 1.;
          repos_vel = repos_frac * d / dt;
        }

        bp->last_repos_vel = (float)repos_vel;
        bp->reposition.delta_x[0] *= repos_frac;
        bp->reposition.delta_x[1] *= repos_frac;
        bp->reposition.delta_x[2] *= repos_frac;
      }

      /* ends section for fractional repositioning */
    } else {

      /* We _should_ reposition, but not fractionally. Here, we will
       * reposition exactly on top of another gas particle - which
       * could cause issues, so we add on a small fractional offset
       * of magnitude 0.001 h in the reposition delta. */

      /* Generate three random numbers in the interval [-0.5, 0.5[; id,
       * id**2, and id**3 are required to give unique random numbers (as
       * random_unit_interval is completely reproducible). */
      const float offset_dx =
          random_unit_interval(bp->id, ti_begin, random_number_BH_reposition) -
          0.5f;
      const float offset_dy =
          random_unit_interval(bp->id * bp->id, ti_begin,
                               random_number_BH_reposition) -
          0.5f;
      const float offset_dz =
          random_unit_interval(bp->id * bp->id * bp->id, ti_begin,
                               random_number_BH_reposition) -
          0.5f;

      const float length_inv =
          1.0f / sqrtf(offset_dx * offset_dx + offset_dy * offset_dy +
                       offset_dz * offset_dz);

      const float norm = 0.001f * bp->h * length_inv;

      bp->reposition.delta_x[0] += offset_dx * norm;
      bp->reposition.delta_x[1] += offset_dy * norm;
      bp->reposition.delta_x[2] += offset_dz * norm;
    }
  } /* ends section if we found eligible repositioning target(s) */
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * This is the equivalent of hydro_reset_acceleration.
 * We do not compute the acceleration on black hole, therefore no need to use
 * it.
 *
 * @param bp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void black_holes_reset_feedback(
    struct bpart* restrict bp) {

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    bp->ids_ngbs_force[i] = -1;
  bp->num_ngb_force = 0;
#endif
}

/**
 * @brief Store the gravitational potential of a black hole by copying it from
 * its #gpart friend.
 *
 * @param bp The black hole particle.
 * @param gp The black hole's #gpart.
 */
__attribute__((always_inline)) INLINE static void
black_holes_store_potential_in_bpart(struct bpart* bp, const struct gpart* gp) {

#ifdef SWIFT_DEBUG_CHECKS
  if (bp->gpart != gp) error("Copying potential to the wrong black hole!");
#endif

  bp->reposition.potential = gp->potential;
}

/**
 * @brief Store the gravitational potential of a particle by copying it from
 * its #gpart friend.
 *
 * @param p_data The black hole data of a gas particle.
 * @param gp The black hole's #gpart.
 */
__attribute__((always_inline)) INLINE static void
black_holes_store_potential_in_part(struct black_holes_part_data* p_data,
                                    const struct gpart* gp) {
  p_data->potential = gp->potential;
}

/**
 * @brief Initialise a BH particle that has just been seeded.
 *
 * @param bp The #bpart to initialise.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param p The #part that became a black hole.
 * @param xp The #xpart that became a black hole.
 */
INLINE static void black_holes_create_from_gas(
    struct bpart* bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo,
    const struct part* p, const struct xpart* xp,
    const integertime_t ti_current) {

  /* All the non-basic properties of the black hole have been zeroed
   * in the FOF code. We update them here.
   * (i.e. position, velocity, mass, time-step have been set) */

  /* Birth time and density */
  bp->formation_scale_factor = cosmo->a;
  //bp->formation_gas_density = hydro_get_physical_density(p, cosmo);
  bp->formation_gas_density = cooling_get_subgrid_density(p, xp);

  /* Initial seed mass */
  bp->subgrid_mass = props->subgrid_seed_mass;

  /* We haven't accreted anything yet */
  bp->total_accreted_mass = 0.f;
  bp->cumulative_number_seeds = 1;
  bp->number_of_mergers = 0;
  bp->number_of_gas_swallows = 0;
  bp->number_of_direct_gas_swallows = 0;
  bp->number_of_time_steps = 0;

  /* We haven't repositioned yet, nor attempted it */
  bp->number_of_repositions = 0;
  bp->number_of_reposition_attempts = 0;
  bp->last_repos_vel = 0.f;

  /* Copy over the splitting struct */
  bp->split_data = xp->split_data;

  /* Initial metal masses */
  const float gas_mass = hydro_get_mass(p);
  struct chemistry_bpart_data* bp_chem = &bp->chemistry_data;
  const struct chemistry_part_data* p_chem = &p->chemistry_data;
  chemistry_bpart_from_part(bp_chem, p_chem, gas_mass);

  /* No swallowed angular momentum */
  bp->swallowed_angular_momentum[0] = 0.f;
  bp->swallowed_angular_momentum[1] = 0.f;
  bp->swallowed_angular_momentum[2] = 0.f;

  /* Last time of mergers */
  bp->last_minor_merger_time = -1.;
  bp->last_major_merger_time = -1.;

  /* First initialisation */
  black_holes_init_bpart(bp);

  black_holes_mark_bpart_as_not_swallowed(&bp->merger_data);
}

#endif /* SWIFT_SIMBA_BLACK_HOLES_H */
