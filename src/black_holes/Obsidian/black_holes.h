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
#ifndef SWIFT_OBSIDIAN_BLACK_HOLES_H
#define SWIFT_OBSIDIAN_BLACK_HOLES_H

/* Local includes */
#include "black_holes_properties.h"
#include "black_holes_struct.h"
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
#include <gsl/gsl_poly.h>


/**
 * @brief How much of the feedback actually couples to the medium?
 *
 * @param props The properties of the black hole scheme.
 * @param BH_state The current state of the black hole.
 */
__attribute__((always_inline)) INLINE static double get_black_hole_coupling(
    const struct black_holes_props* props, const struct cosmology* cosmo, 
    const int BH_state) {
  switch (BH_state) {
    case BH_states_adaf:
      const double scaling =
          min(pow(1. + cosmo->z, props->adaf_z_scaling), 1.);
      return props->adaf_coupling * scaling;
      break;
    case BH_states_quasar:
      return props->quasar_coupling;
      break;
    case BH_states_slim_disk:
      return props->slim_disk_coupling;
      break;
    default:
      error("Invalid black hole state.");
      return 0.;
      break;
  }
}

/**
 * @brief Computes the radiative efficiency in the slim disk mode.
 *
 * @param props The properties of the black hole scheme.
 * @param f_Edd M_dot,BH / M_dot,Edd
 */
__attribute__((always_inline)) INLINE static 
double get_black_hole_slim_disk_efficiency(
    const struct black_holes_props* props, 
    const double f_Edd) {
  if (f_Edd <= 0.) return 0.;
  const double R = 1. / f_Edd;
  /* Efficiency from Lupi et al. (2014), 
   * super eddington accretion and feedback */
  return (R / 16.) * props->A_lupi * 
         (0.985 / (R + (5. / 8.) * props->B_lupi) + 0.015 / 
            (R + (5. / 8.) * props->C_lupi));
}

/**
 * @brief Computes the radiative efficiency in the ADAF mode.
 *
 * @param props The properties of the black hole scheme.
 * @param f_Edd M_dot,BH / M_dot,Edd
 */
__attribute__((always_inline)) INLINE static 
double get_black_hole_adaf_efficiency(
    const struct black_holes_props* props, const double f_Edd) {
  return props->epsilon_r * f_Edd;  /* scales with M_dot,BH */
}

/**
 * @brief Chooses and calls the proper radiative efficiency function for the state.
 *
 * @param props The properties of the black hole scheme.
 * @param f_Edd The accretion rate over the Eddington rate.
 * @param BH_state The current state of the BH.
 */
__attribute__((always_inline)) INLINE static 
double get_black_hole_radiative_efficiency(
    const struct black_holes_props* props, 
    const double f_Edd, const int BH_state) {
  switch(BH_state) {
    case BH_states_adaf:
      return get_black_hole_adaf_efficiency(props, f_Edd);
    case BH_states_quasar:
      return props->epsilon_r;
    case BH_states_slim_disk:
      return get_black_hole_slim_disk_efficiency(props, f_Edd);
    default:
      error("Invalid black hole state.");
      break;
  }

  return 0.;
}

/**
 * @brief Computes the fraction of M_dot,inflow that should go into the BH.
 *
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 * @param m_dot_inflow_m_dot_edd M_dot,inflow scaled to M_dot,Edd for the BH.
 */
__attribute__((always_inline)) INLINE static double 
get_black_hole_upper_mdot_medd(
    const struct black_holes_props* props, 
    const struct phys_const* constants,
    const double m_dot_inflow_m_dot_edd) {

  if (m_dot_inflow_m_dot_edd <= 0.) return 0.;

  double x1, x2, x3;
  double a3, a2, a1, a0;
  const double phi = props->slim_disk_phi;

  int num_roots;

  a3 = ((5. * 5.) / (8. * 8.)) * props->B_lupi * props->C_lupi;
  a2 = (5. / 8.) * ((props->B_lupi + props->C_lupi) + (phi / 16.) * 
       props->A_lupi * (0.015 * props->B_lupi + 0.985 * props->C_lupi) - 
       (5. / 8.) * props->B_lupi * props->C_lupi * m_dot_inflow_m_dot_edd);
  a1 = 1. + (phi / 16.) * props->A_lupi - (5. / 8.) * 
       (props->B_lupi + props->C_lupi) * m_dot_inflow_m_dot_edd;
  a0 = -m_dot_inflow_m_dot_edd;

  a2 /= a3;
  a1 /= a3;
  a0 /= a3;

  num_roots = gsl_poly_solve_cubic(a2, a1, a0, &x1, &x2, &x3);
  if (num_roots == 1) {
    if (x1 >= 0.) {
      return x1;
    } 
    else {
      warning("num_roots=1 m_dot_inflow_m_dot_edd=%g phi=%g a3=%g a2=%g "
              "a1=%g a0=%g",
              m_dot_inflow_m_dot_edd, phi, a3, a2, a1, a0);
      return 0.;
    }
  }
  if (x3 >= 0.) {
    return x3;
  } 
  else {
    warning("num_roots=0 m_dot_inflow_m_dot_edd=%g phi=%g a3=%g a2=%g a1=%g a0=%g",
            m_dot_inflow_m_dot_edd, phi, a3, a2, a1, a0);
    return 0.;
  }

  return 0.;
}

/**
 * @brief Computes the fraction of M_dot,inflow that should go into the BH.
 *
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 * @param m_dot_inflow M_dot,inflow in internal units.
 * @param BH_mass The subgrid mass of the BH in internal units.
 * @param BH_state The current state of the BH.
 * @param Eddington_rate M_dot,Edd in internal units.
 */
__attribute__((always_inline)) INLINE static 
double get_black_hole_accretion_factor(
    const struct black_holes_props* props, 
    const struct phys_const* constants,
    const double m_dot_inflow, const double BH_mass, 
    const int BH_state, 
    const double Eddington_rate) {
  
  if (m_dot_inflow <= 0. || BH_mass <= 0.) return 0.;

  
  switch (BH_state) {
    case BH_states_adaf:
      return props->adaf_f_accretion;
      break;
    case BH_states_quasar:
      return props->quasar_f_accretion;
      break;
    case BH_states_slim_disk:
      /* This is the FRACTION of the total so divide by M_dot,inflow */
      const double f_edd = m_dot_inflow / Eddington_rate;
      double mdot_medd = 
          get_black_hole_upper_mdot_medd(props, constants, f_edd);
      return mdot_medd * Eddington_rate / m_dot_inflow;
      break;
    default:
      error("Invalid black hole state.");
      return 0.;
      break;
  }
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
  float dt_overall = FLT_MAX;
  float dt_kick = FLT_MAX;
  const float min_subgrid_mass = props->minimum_black_hole_mass_unresolved;

  /* Only limit when in the resolved feedback regime */
  if (bp->accretion_rate > 0.f && bp->subgrid_mass > min_subgrid_mass) {
    dt_accr = props->dt_accretion_factor * bp->mass / bp->accretion_rate;

    if (bp->state == BH_states_adaf) {
      dt_kick = bp->ngb_mass / (props->jet_mass_loading * bp->accretion_rate);
    }
    else {
      if (bp->f_accretion > 0.f) {
        /* Make sure that the wind mass does not exceed the kernel gas mass */
        const float psi = (1.f - bp->f_accretion) / bp->f_accretion;
        dt_kick = bp->ngb_mass / (psi * bp->accretion_rate);
      }
    }

    dt_overall = min(dt_kick, dt_accr);
  }

  if (dt_overall < props->time_step_min) {
    message(
        "Warning! BH_TIMESTEP_LOW: id=%lld (%g Myr) is below time_step_min (%g "
        "Myr).",
        bp->id, dt_overall * props->time_to_Myr,
        props->time_step_min * props->time_to_Myr);
  }

  return max(dt_overall, props->time_step_min);
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
        "'ObsidianAGN:use_subgrid_mass_from_ics' to 0 to initialize the "
        "black hole subgrid masses to the corresponding dynamical masses.\n"
        "If the subgrid mass is intentionally set to this value, you can "
        "disable this error by setting 'ObsidianAGN:with_subgrid_mass_check' "
        "to 0.",
        bp->id, bp->subgrid_mass);
  }
  bp->total_accreted_mass = 0.f;
  bp->accretion_disk_mass = 0.f;
  bp->gas_SFR = 0.f;
  bp->accretion_rate = 0.f;
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
  bp->accreted_angular_momentum[0] = 0.f;
  bp->accreted_angular_momentum[1] = 0.f;
  bp->accreted_angular_momentum[2] = 0.f;
  bp->last_repos_vel = 0.f;
  bp->radiative_luminosity = 0.f;
  bp->delta_energy_this_timestep = 0.f;
  bp->state = BH_states_slim_disk;
  bp->radiative_efficiency = 0.f;
  bp->f_accretion = 0.f;
  bp->m_dot_inflow = 0.f;
  bp->cold_disk_mass = 0.f;
  bp->jet_mass_reservoir = 0.f;
  bp->jet_mass_kicked_this_step = 0.f;
  bp->adaf_energy_to_dump = 0.f;

#ifdef WITH_FOF_GALAXIES
  bp->group_data.mass = 0.f;
  bp->group_data.stellar_mass = 0.f;
  bp->group_data.ssfr = 0.f;
#endif
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
  bp->hot_gas_internal_energy = 0.f;
  bp->sound_speed_subgrid_gas = -1.f;
  bp->velocity_gas[0] = 0.f;
  bp->velocity_gas[1] = 0.f;
  bp->velocity_gas[2] = 0.f;
  bp->circular_velocity_gas[0] = 0.f;
  bp->circular_velocity_gas[1] = 0.f;
  bp->circular_velocity_gas[2] = 0.f;
  bp->angular_momentum_gas[0] = 0.f;
  bp->angular_momentum_gas[1] = 0.f;
  bp->angular_momentum_gas[2] = 0.f;
  bp->stellar_mass = 0.f;
  bp->stellar_bulge_mass = 0.f;
  bp->radiative_luminosity = 0.f;
  bp->ngb_mass = 0.f;
  bp->gravitational_ngb_mass = 0.f;
  bp->num_ngbs = 0;
  bp->num_gravitational_ngbs = 0;
  bp->reposition.delta_x[0] = -FLT_MAX;
  bp->reposition.delta_x[1] = -FLT_MAX;
  bp->reposition.delta_x[2] = -FLT_MAX;
  bp->reposition.min_potential = FLT_MAX;
  bp->reposition.potential = FLT_MAX;
  bp->accretion_rate = 0.f; /* Optionally accumulated ngb-by-ngb */
  bp->cold_disk_mass = 0.f;
  bp->mass_at_start_of_step = bp->mass; /* bp->mass may grow in nibbling mode */
  bp->m_dot_inflow = 0.f; /* reset accretion rate */
  bp->adaf_energy_to_dump = 0.f;
  bp->adaf_wt_sum = 0.f;
  bp->kernel_wt_sum = 0.f;
  /* update the reservoir */
  bp->jet_mass_reservoir -= bp->jet_mass_kicked_this_step;
  bp->jet_mass_kicked_this_step = 0.f;
  if (bp->jet_mass_reservoir < 0.f) {
    bp->jet_mass_reservoir = 0.f; /* reset reservoir if used up */
  }
  /* update the unresolved reservoir */
  bp->unresolved_mass_reservoir -= bp->unresolved_mass_kicked_this_step;
  bp->unresolved_mass_kicked_this_step = 0.f;
  if (bp->unresolved_mass_reservoir < 0.f) {
    bp->unresolved_mass_reservoir = 0.f; 
  }
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

  /* For the following, we also have to undo the mass smoothing
   * (N.B.: bp->velocity_gas is in BH frame, in internal units). */
  bp->sound_speed_gas *= h_inv_dim * rho_inv;
  bp->internal_energy_gas *= h_inv_dim * rho_inv;

  /* Non-weighted (no decoupled winds) properties below.
   * All mass-weighted quantities are for the hot & cold gas */
  float m_hot_inv = 1.f;
  if (bp->hot_gas_mass > 0.f) m_hot_inv /= bp->hot_gas_mass;
  /* Or the total mass */
  float m_tot_inv = 1.f;
  if (bp->ngb_mass > 0.f) m_tot_inv /= bp->ngb_mass;

  bp->hot_gas_internal_energy *= m_hot_inv;
  bp->velocity_gas[0] *= m_tot_inv;
  bp->velocity_gas[1] *= m_tot_inv;
  bp->velocity_gas[2] *= m_tot_inv;
  bp->circular_velocity_gas[0] *= m_tot_inv;
  bp->circular_velocity_gas[1] *= m_tot_inv;
  bp->circular_velocity_gas[2] *= m_tot_inv;

  /* Calculate circular velocity at the smoothing radius from specific
   * angular momentum (extra h_inv). It is now a VELOCITY.
   */
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
__attribute__((always_inline)) INLINE static void black_holes_bpart_has_no_neighbours(struct bpart* bp,
                                    const struct cosmology* cosmo) {

  //warning(
  //    "BH particle with ID %lld treated as having no neighbours (h: %g, "
  //    "wcount: %g).",
  //    bp->id, bp->h, bp->density.wcount);

  /* Some smoothing length multiples. */
  const float h = bp->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  bp->density.wcount = kernel_root * h_inv_dim;
  bp->density.wcount_dh = 0.f;

  bp->velocity_gas[0] = FLT_MAX;
  bp->velocity_gas[1] = FLT_MAX;
  bp->velocity_gas[2] = FLT_MAX;

  bp->internal_energy_gas = -FLT_MAX;
  bp->hot_gas_internal_energy = -FLT_MAX;
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
 * @brief Return the current bolometric luminosity of the BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_bolometric_luminosity(const struct bpart* bp,
                                      const struct phys_const* constants) {
  const double c = constants->const_speed_light_c;
  return bp->accretion_rate * bp->radiative_efficiency * c * c;
}

/**
 * @brief Return the current kinetic jet power of the BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double black_holes_get_jet_power(
    const struct bpart* bp, const struct phys_const* constants,
    const struct black_holes_props* props) {
  const double c = constants->const_speed_light_c;
  double eta_jet = props->jet_efficiency;
  if (bp->state != BH_states_adaf) {
    eta_jet = FLT_MIN;
  }
  /* accretion_rate is M_dot,acc from the paper */
  return eta_jet * bp->accretion_rate * c * c;
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
  bp->num_gravitational_ngbs--;
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

  /* Update the energy reservoir */
  bpi->jet_mass_reservoir += bpj->jet_mass_reservoir;

  /* Add up all the BH seeds */
  bpi->cumulative_number_seeds += bpj->cumulative_number_seeds;

  /* Add up all the gas particles we swallowed */
  bpi->number_of_gas_swallows += bpj->number_of_gas_swallows;

  /* Add the subgrid angular momentum that we swallowed */
  bpi->accreted_angular_momentum[0] += bpj->accreted_angular_momentum[0];
  bpi->accreted_angular_momentum[1] += bpj->accreted_angular_momentum[1];
  bpi->accreted_angular_momentum[2] += bpj->accreted_angular_momentum[2];

  /* We had another merger */
  bpi->number_of_mergers++;
}

/**
 * @brief Compute the wind launch speed for this feedback step.
 *
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 * @param bp The black hole particle.
 * @param Eddington_rate M_dot,Edd for the black hole.
 */
__attribute__((always_inline)) INLINE static double get_black_hole_wind_speed(
    const struct black_holes_props* props,
    const struct phys_const* constants,
    const struct bpart *bp,
    const double Eddington_rate) {

  if (bp->accretion_rate < 0.f || bp->m_dot_inflow < 0.f) return 0.f;

  switch (bp->state) {   
    case BH_states_adaf:
      return props->adaf_wind_speed;
      break;
    case BH_states_quasar:
      return props->quasar_wind_speed;
      break;
    case BH_states_slim_disk:
      return props->slim_disk_wind_speed;
      break;
    default:
      error("Invalid black hole state.");
      return 0.f;
      break;
  }
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

  if (dt == 0. || bp->rho_gas == 0. || bp->h == 0.) return;

  /* A black hole should never accrete/feedback if it is not in a galaxy */
  if (bp->group_data.mass <= 0.f) return;

  /* Gather some physical constants (all in internal units) */
  const double G = constants->const_newton_G;
  const double c = constants->const_speed_light_c;
  const double proton_mass = constants->const_proton_mass;
  const double sigma_Thomson = constants->const_thomson_cross_section;

  /* Gather the parameters of the model */
  const double f_Edd_maximum = props->f_Edd_maximum;

#ifdef OBSIDIAN_DEBUG_CHECKS
  if (isnan(bp->subgrid_mass)) error("subgrid_mass nan");
#endif
  /* (Subgrid) mass of the BH (internal units) */
  const double BH_mass = bp->subgrid_mass;

  /* Compute the Eddington rate (internal units).
   * IMPORTANT: epsilon_r = 0.1 is the SET value for the Eddington rate.
   * It is assumed in all of the derivations for the state model.
   */
  const double Eddington_rate =
      4. * M_PI * G * BH_mass * proton_mass / (0.1 * c * sigma_Thomson);

  const double bh_h = kernel_gamma * bp->h;
  const double bh_h_inv = 1. / bh_h;
  const double volume_bh_inv = 
      (3. / (4. * M_PI)) * bh_h_inv * bh_h_inv * bh_h_inv;
  double gas_rho = 0.;
  if (props->bondi_use_all_gas) {
    gas_rho = bp->rho_gas;
  }
  else {
    gas_rho = bp->hot_gas_mass * volume_bh_inv;
  }

  const double gas_rho_phys = gas_rho * cosmo->a3_inv;

  /* We can now compute the Bondi accretion rate (internal units) 
   * D. Rennehan: In Simba, we only consider the hot gas within
   * the kernel for the Bondi rate, and the cold gas using the
   * torque accretion estimator.
   */
  double Bondi_rate = 0.;

  /* Use all gas within the kernel or only the hot gas*/
  double gas_internal_energy = 0.;
  if (props->bondi_use_all_gas) {
    if (bp->internal_energy_gas > 0.) {
      gas_internal_energy = bp->internal_energy_gas;
    }
  }
  else {
    if (bp->hot_gas_internal_energy > 0.) {
      gas_internal_energy = bp->hot_gas_internal_energy;
    }
  }

  /* Check if there is hot/any gas in the kernel */
  if (gas_internal_energy > 0.) {
    double gas_c = 
        gas_soundspeed_from_internal_energy(gas_rho, gas_internal_energy);

    if (gas_c > 0.) {
      const double gas_c_phys_inv = 
          1. / (cosmo->a_factor_sound_speed * gas_c);

      Bondi_rate = 4. * M_PI * G * G * BH_mass * BH_mass * gas_rho_phys *
                    gas_c_phys_inv * gas_c_phys_inv * gas_c_phys_inv;

      /* In the case of standard Bondi, we limit it to the Eddington rate */
      Bondi_rate = fmin(Bondi_rate, props->f_Edd_Bondi_maximum * Eddington_rate);
    }
  }

  /* The accretion rate estimators give Mdot,inflow  
    * (Mdot,BH = f_acc * Mdot,inflow) */
  const double bondi_accr_rate = props->bondi_alpha * Bondi_rate;

  /* Compute the torque-limited accretion rate */
  double torque_accr_rate = 0.;

  double f_corr_stellar = 10.;
  const double galaxy_gas_mass = 
      bp->group_data.mass - bp->group_data.stellar_mass;
  if (galaxy_gas_mass > 0.) {
    f_corr_stellar = 
        min(bp->group_data.stellar_mass / galaxy_gas_mass, f_corr_stellar);
  }

  /* Torque accretion rate based on some fraction of gas near BH 
   * falling in on dynamical time. (This is the default.)
   * Here the accretion rate is only based on Mgas / tdyn.
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
  double tdyn_inv = FLT_MAX;
  const float potential = fabs(gravity_get_comoving_potential(bp->gpart));
  /* Includes dynamical mass of the BH */
  double total_mass = gravity_get_total_mass(bp->gpart);
  switch (props->dynamical_time_calculation_method) {
    /* Assume gas fraction is the same in the kernel and outside */
    case 0:
      /* Compute correction to total dynamical mass around 
       * BH contributed by stars */
      const double m_star_gal = bp->group_data.stellar_mass;
      const double m_gas_cold_gal = 
          bp->group_data.mass - bp->group_data.stellar_mass;
      const double m_gas_bh = bp->gravitational_ngb_mass;
      const double m_bh = bp->mass;

      /* Compute stellar mass assuming a constant cold gas fraction in the
       * entire galaxy. If m_gas_cold_bh is zero it doesn't matter since
       * the BH won't acrrete in the torque mode anyway. */
      const double m_gas_cold_bh = bp->cold_gas_mass;
      double m_star_bh = 0.;
      if (m_gas_cold_gal > 0.) {
        m_star_bh = (m_star_gal / m_gas_cold_gal) * m_gas_cold_bh;
      }

      /* Have to assume baryon dominance within the kernel */
      const double rho_est = (m_star_bh + m_gas_bh + m_bh) * volume_bh_inv;

      /* Inverse physical dynamical time */
      tdyn_inv = sqrt(32. * G  * rho_est * cosmo->a3_inv / (3. * M_PI));
      break;

    /* Assume BH potential */
    case 1:
      if (potential >= 0.f) {
        tdyn_inv = (sqrt(potential) / bh_h) * cosmo->a2_inv;
      }
      break;

    /* Assume dynamical time from the kernel mass */
    case 2:
      /* do not have gravity_props here */
      const double eps = gravity_get_softening(bp->gpart, NULL);
      const double volume = (4.f * M_PI / 3.f) * eps * eps * eps;
      const double rho = total_mass / volume;
      tdyn_inv = sqrt(32. * G * rho * cosmo->a3_inv / (3. * M_PI));
      break;

    default:
      error("Unknown dynamical time calculation method %d", 
            props->dynamical_time_calculation_method);
      break;
  }

  /* Compute infall times to BH at this redshift */
  double t_infall = props->bh_accr_dyn_time_fac / tdyn_inv;

  /* If the input value is above 10, assume it is constant in Myr, 
   * scaled with 1/H */
  if (props->bh_accr_dyn_time_fac > 10.) {
    t_infall = 
        props->bh_accr_dyn_time_fac * cosmo->H0 / 
          (props->time_to_Myr * cosmo->H);
  }

  const double corot_gas_mass = 
        bp->cold_gas_mass - 2. * (bp->cold_gas_mass - bp->cold_disk_mass);
  if (props->torque_accretion_norm > 0.f) {
    switch (props->torque_accretion_method) {
      case 0:
        if (galaxy_gas_mass > 0.) {
          torque_accr_rate = 
              props->torque_accretion_norm * bp->cold_disk_mass * tdyn_inv;
        }
        break;

      case 1:
        if (corot_gas_mass > 0. && bp->cold_gas_mass > 0.) {
          torque_accr_rate = 
              props->torque_accretion_norm * corot_gas_mass * tdyn_inv;
        }
        break;

      case 2:
        if (corot_gas_mass > 0. && bp->cold_gas_mass > 0.) {
          const double m_disk = bp->cold_gas_mass * f_corr_stellar;
          const double f_disk = corot_gas_mass / bp->cold_gas_mass;

          const double r0 = bh_h * cosmo->a * (props->length_to_parsec * 0.01);

          const double alpha = 5.;
          const double mass_to_1e9solar = props->mass_to_solar_mass * 1.0e-9;
          const double mass_to_1e8solar = props->mass_to_solar_mass * 1.0e-8;

          const double f0 =
              0.31 * f_disk * f_disk * 
                pow(m_disk * mass_to_1e9solar, -1. / 3.);
          const double f_gas = corot_gas_mass / m_disk;
          const double mass_in_1e8solar = bp->subgrid_mass * mass_to_1e8solar;

          torque_accr_rate = props->torque_accretion_norm * alpha *
                             corot_gas_mass * mass_to_1e9solar *
                             powf(f_disk, 5. / 2.) *
                             powf(mass_in_1e8solar, 1. / 6.) *
                             powf(r0, -3. / 2.) / (1. + f0 / f_gas);
          torque_accr_rate *= 
              (props->time_to_yr / props->mass_to_solar_mass);
        }
        break;

        case 3:
          if (galaxy_gas_mass > 0.) {
            torque_accr_rate = 
              props->torque_accretion_norm * bp->cold_gas_mass * tdyn_inv;
          }
          break;

        default:
          error("Unknown torque_accretion_method=%d", 
                props->torque_accretion_method);
          break;
    }

    /* Do suppression of BH growth */
    switch (props->suppress_growth) {
      case 1:
        const double r0 = bh_h * cosmo->a * props->length_to_parsec;
        const double sigma_eff = f_corr_stellar * bp->ngb_mass *
                                 props->mass_to_solar_mass /
                                 (M_PI * r0 * r0);

        torque_accr_rate *= sigma_eff / (sigma_eff + 3000.);
        break;

      case 2:
      case 6:
      case 7:
        double m_suppress = fabs(props->bh_characteristic_suppression_mass);
        if (props->bh_characteristic_suppression_mass < 0) {
          m_suppress *= cosmo->a;
        }

        torque_accr_rate *= 
            1. - exp(-bp->subgrid_mass * props->mass_to_solar_mass / m_suppress);
        break;

      case 4:
      case 5:
        /* compute mass loading factor from SF feedback, 
         * should be same as used in feedback_mass_loading_factor() 
         */
        const double galaxy_stellar_mass = 
            max(bp->group_data.stellar_mass, 5.8e8 / props->mass_to_solar_mass);
        double slope = -0.317;

        const double min_mass = 5.2e9 / props->mass_to_solar_mass;
        if (galaxy_stellar_mass > min_mass) slope = -0.716;

        const double eta = 
            12. * pow(galaxy_stellar_mass / min_mass, slope);

        /* compute fraction of mass within kernel going into outflows 
         * over accretion time 
         */
        double sfr = 0.;
        if (bp->cold_gas_mass > 0.f) {
          sfr = bp->group_data.ssfr * bp->group_data.stellar_mass;
          sfr *= props->torque_accretion_norm;
        }

        /* suppress accretion by factor accounting for mass lost
         * in SF-driven outflow 
         */
        if (sfr > 0.) {
          torque_accr_rate *= torque_accr_rate / (torque_accr_rate + eta * sfr);
        }
        break;

      default:
        break;
    }
  }

#ifdef OBSIDIAN_DEBUG_CHECKS
  if (isnan(bondi_accr_rate)) error("bondi_accr_rate nan");
  if (isnan(torque_accr_rate)) error("torque_accr_rate nan");
#endif

  /* Right now this is M_dot,inflow. We will multiply by 
   * f_accretion later to make it M_dot,acc */
  bp->accretion_rate = bondi_accr_rate + torque_accr_rate;

  /* We will use eddington_fraction_lower_boundary and 
   * eddington_fraction_upper_boundary to divide up the accretion rate 
   * in three regimes.
   * 
   * In order to switch out of a regime (i.e. a state), it is necessary 
   * for the true accretion rate (compared to Eddington rate) to switch 
   * over the boundary. Therefore, before we switch a state we must calculate 
   * what the previous state predicts the true accretion rate onto the SMBH is, 
   * and then update the state if it crosses a boundary.
   */
  double predicted_mdot_medd = 0.;
  switch (bp->state) {
    case BH_states_adaf:
      predicted_mdot_medd 
            = bp->accretion_rate * props->adaf_f_accretion / Eddington_rate;

      if (predicted_mdot_medd > props->eddington_fraction_upper_boundary) {
        bp->state = BH_states_slim_disk;
        break;
      }
      if (predicted_mdot_medd > props->eddington_fraction_lower_boundary) {
        bp->state = BH_states_quasar;
      }

      break; /* end case ADAF */
    case BH_states_quasar:
      predicted_mdot_medd 
            = bp->accretion_rate * props->quasar_f_accretion / Eddington_rate;

      if (BH_mass > props->adaf_mass_limit &&
          predicted_mdot_medd < props->eddington_fraction_lower_boundary) {
        bp->state = BH_states_adaf;
        break;
      }

      if (predicted_mdot_medd > props->eddington_fraction_upper_boundary) {
        bp->state = BH_states_slim_disk;
      }
  
      break; /* end case quasar */
    case BH_states_slim_disk:
      predicted_mdot_medd = 
          get_black_hole_upper_mdot_medd(props, constants, 
                                         bp->accretion_rate / Eddington_rate);

      if (BH_mass > props->adaf_mass_limit &&
          predicted_mdot_medd < props->eddington_fraction_lower_boundary) {
        bp->state = BH_states_adaf;
        break;
      }

      if (predicted_mdot_medd < props->eddington_fraction_upper_boundary) {
        bp->state = BH_states_quasar;
      }

      break; /* end case slim disk */
    default:
      error("Invalid black hole state.");
      break;
  }

  /* We need to store the full M_dot,inflow rate to calculate the 
   * fraction at high accretion rate */
  bp->m_dot_inflow = bp->accretion_rate;

  /* This depends on the new state */
  bp->f_accretion = 
      get_black_hole_accretion_factor(props, constants, bp->m_dot_inflow,
                                      BH_mass, bp->state, Eddington_rate);
#ifdef OBSIDIAN_DEBUG_CHECKS
  if (isnan(bp->f_accretion)) error("f_accretion nan");
#endif
  if (bp->f_accretion <= 1.e-10f) bp->f_accretion = 0.f;

  /* bp->accretion_rate is M_dot,acc in Rennehan+24 */
  bp->accretion_rate *= bp->f_accretion;

  if (!props->bondi_use_all_gas) {
    /* Now we can Eddington limit. */
    bp->accretion_rate = 
        min(bp->accretion_rate, f_Edd_maximum * Eddington_rate);
  }

  /* All accretion is done, now we can set the eddington fraction */
  bp->eddington_fraction = bp->accretion_rate / Eddington_rate;

  /* Get the new radiative efficiency based on the new state */
  bp->radiative_efficiency = 
      get_black_hole_radiative_efficiency(props, bp->eddington_fraction, 
                                          bp->state);
#ifdef OBSIDIAN_DEBUG_CHECKS
  if (isnan(bp->radiative_efficiency)) error("radiative_efficiency nan");
#endif
  if (bp->radiative_efficiency < 1.e-10f) bp->radiative_efficiency = 0.f;

  double mass_rate = 0.;
  const double luminosity = 
      bp->radiative_efficiency * bp->accretion_rate * c * c;

  /* Factor in the radiative efficiency, don't subtract 
   * jet BZ efficiency (spin is fixed) */
  mass_rate = (1. - bp->radiative_efficiency) * bp->accretion_rate;

  /* This is used for X-ray feedback later */
  bp->radiative_luminosity = luminosity;

#ifdef OBSIDIAN_DEBUG_CHECKS
  if (isnan(dt)) error("dt nan");
#endif
  /* Integrate forward in time */
  double delta_mass = mass_rate * dt;

  /* If desired we put mass into accretion disk which feeds BH on some 
   * frac of tdyn 
   */
  if (t_infall > 0.f) {
    /* Add accreted mass into a reservoir representing BH accretion disk */
    bp->accretion_disk_mass += delta_mass;

    /* Compute mass that will actually go into BH */
    delta_mass = bp->accretion_disk_mass * (1. - exp(-dt / t_infall));

    /* This mass gets removed from the accretion disk */
    if (bp->accretion_disk_mass > delta_mass) {
      bp->accretion_disk_mass -= delta_mass;
    }
    else {
      delta_mass = bp->accretion_disk_mass;
      bp->accretion_disk_mass = 0.;
    }

    /* Recompute accretion rate based on the reservoir change */
    bp->accretion_rate = delta_mass / (dt * (1. - bp->radiative_efficiency));
  }

  bp->subgrid_mass += delta_mass;
  bp->total_accreted_mass += delta_mass;

  if (bp->state == BH_states_adaf) {
    /* ergs to dump in a kernel-weighted fashion */
    if (props->adaf_wind_mass_loading == 0.f) {
      bp->adaf_energy_to_dump = 
          get_black_hole_coupling(props, cosmo, bp->state) *
            props->adaf_disk_efficiency * bp->accretion_rate * c * c;
    }
    else {
      const float adaf_v2 = props->adaf_wind_speed * props->adaf_wind_speed;
      const float mass_this_step = 
          props->adaf_wind_mass_loading * bp->accretion_rate * dt;
      bp->adaf_energy_to_dump = 0.5f * mass_this_step * adaf_v2;
    }
  }

  if (bp->state == BH_states_adaf || 
        (props->slim_disk_jet_active && bp->state == BH_states_slim_disk)) { 
    /* Psi_jet*M_dot,acc*dt is the total mass expected in the jet this step */
    bp->jet_mass_reservoir += props->jet_mass_loading * bp->accretion_rate * dt;
  }

  if (bp->subgrid_mass < bp->mass) {
    /* In this case, the BH is still accreting from its (assumed) subgrid gas
     * mass reservoir left over when it was formed. There is some loss in this
     * due to radiative losses, so we must decrease the particle mass
     * in proportion to its current accretion rate. We do not account for this
     * in the swallowing approach, however. */
    bp->mass -= bp->radiative_efficiency * bp->accretion_rate * dt;

    if (bp->mass < 0) {
      error("Black hole %lld reached negative mass (%g). Trouble ahead...",
            bp->id, bp->mass);
    }

    /* Make sure not to destroy low mass galaxies */
    if (bp->subgrid_mass > props->minimum_black_hole_mass_unresolved &&
        bp->state != BH_states_adaf) {
      /* Make sure if many mergers have driven up the dynamical mass at low
      * subgrid mass, that we still kick out particles! */
      const float psi = (1.f - bp->f_accretion) / bp->f_accretion;
      bp->unresolved_mass_reservoir += psi * bp->accretion_rate * dt;
    }
  }

  /* Increase the subgrid angular momentum according to what we accreted
   * Note that this is already in physical units, a factors from velocity and
   * radius cancel each other. Also, the circular velocity contains an extra
   * smoothing length factor that we undo here. */
  const double m_times_r = (mass_rate * dt) * bp->h;
  /* L = m * r * v */
  bp->accreted_angular_momentum[0] += m_times_r * bp->circular_velocity_gas[0];
  bp->accreted_angular_momentum[1] += m_times_r * bp->circular_velocity_gas[1];
  bp->accreted_angular_momentum[2] += m_times_r * bp->circular_velocity_gas[2];

  /* Keep v_kick physical, there are a lot of comparisons */
  bp->v_kick = 
      get_black_hole_wind_speed(props, constants, bp, Eddington_rate);

  /* This is always true in the ADAF mode; only heating happens */
  if (bp->state == BH_states_adaf) {
    bp->v_kick = 0.f;
  }

#ifdef OBSIDIAN_DEBUG_CHECKS
  tdyn_inv = (tdyn_inv > 0.f) ? tdyn_inv : FLT_MIN;
  message("BH_ACC: z=%g bid=%lld ms=%g dms=%g sfr=%g mbh=%g dmbh=%g state=%d "
          "torque=%g bondi=%g fEdd=%g facc=%g fsupp=%g mcold=%g mhot=%g mdisk=%g"
          " tin=%g vkick=%g dmass=%g radeff=%g mres=%g tdyn=%g", 
          cosmo->z, 
          bp->id, 
          bp->group_data.stellar_mass * props->mass_to_solar_mass,
          bp->group_data.ssfr * bp->group_data.stellar_mass * dt * 
              props->mass_to_solar_mass,
          bp->group_data.ssfr * bp->group_data.stellar_mass * 
              props->mass_to_solar_mass / props->time_to_yr, 
          bp->subgrid_mass * props->mass_to_solar_mass,
          delta_mass * props->mass_to_solar_mass,
          bp->state,
          torque_accr_rate * props->mass_to_solar_mass / props->time_to_yr,
          bondi_accr_rate * props->mass_to_solar_mass / props->time_to_yr,
          bp->eddington_fraction,
          bp->f_accretion,
          1. - exp(-bp->subgrid_mass * props->mass_to_solar_mass / 
                   fabs(props->bh_characteristic_suppression_mass) * cosmo->a),
          bp->cold_gas_mass * props->mass_to_solar_mass,
          bp->hot_gas_mass * props->mass_to_solar_mass,
          corot_gas_mass * props->mass_to_solar_mass,
          t_infall * props->time_to_Myr,
          bp->v_kick / props->kms_to_internal,
          delta_mass, 
          bp->radiative_efficiency, 
          bp->accretion_disk_mass, 
          (1.f / tdyn_inv) * props->time_to_Myr);

  message("BH_STATES: id=%lld, new_state=%d, predicted_mdot_medd=%g, "
          "eps_r=%g, f_Edd=%g, f_acc=%g, "
          "luminosity=%g, accr_rate=%g Msun/yr, coupling=%g, v_kick=%g km/s, "
          "jet_mass_reservoir=%g Msun unresolved_reservoir=%g Msun",
          bp->id,
          bp->state, 
          predicted_mdot_medd, 
          bp->radiative_efficiency,
          bp->eddington_fraction,
          bp->f_accretion, 
          bp->radiative_luminosity * props->conv_factor_energy_rate_to_cgs, 
          bp->accretion_rate * props->mass_to_solar_mass / props->time_to_yr,  
          get_black_hole_coupling(props, cosmo, bp->state), 
          bp->v_kick / props->kms_to_internal,
          bp->jet_mass_reservoir * props->mass_to_solar_mass,
          bp->unresolved_mass_reservoir * props->mass_to_solar_mass);
#endif

  printf("BH_DETAILS "
         "%2.12f %lld "
         " %g %g %g %g %g %g %g "
         " %g %g %g %g " 
         " %g %g %g %g %g "
         " %2.10f %2.10f %2.10f "
         " %2.7f %2.7f %2.7f "
         " %g %g %g  %g %g %g"
         " %g %d %g %g"
         " %g %g\n",
         cosmo->a,
         bp->id,
         bp->mass * props->mass_to_solar_mass, 
         bp->subgrid_mass * props->mass_to_solar_mass, 
         total_mass * props->mass_to_solar_mass, 
         bp->accretion_rate * props->mass_to_solar_mass / props->time_to_yr, 
         Bondi_rate * props->mass_to_solar_mass / props->time_to_yr, 
         torque_accr_rate * props->mass_to_solar_mass / props->time_to_yr, 
         dt * props->time_to_Myr,
         (bp->rho_gas * cosmo->a3_inv) * props->rho_to_n_cgs, 
         bp->hot_gas_internal_energy * cosmo->a_factor_internal_energy * 
             props->conv_factor_specific_energy_to_cgs, 
         bp->gas_SFR * props->mass_to_solar_mass / props->time_to_yr, 
         bp->ngb_mass * props->mass_to_solar_mass, 
         bp->hot_gas_mass * props->mass_to_solar_mass, 
         bp->stellar_mass * props->mass_to_solar_mass, 
         0.f /* Mgas,bulge */, 
         bp->stellar_bulge_mass * props->mass_to_solar_mass, 
         0.f,
         bp->x[0] * cosmo->a * props->length_to_parsec / 1.0e3f, 
         bp->x[1] * cosmo->a * props->length_to_parsec / 1.0e3f, 
         bp->x[2] * cosmo->a * props->length_to_parsec / 1.0e3f, 
         bp->v[0] * cosmo->a_inv / props->kms_to_internal, 
         bp->v[1] * cosmo->a_inv / props->kms_to_internal, 
         bp->v[2] * cosmo->a_inv / props->kms_to_internal,
         bp->angular_momentum_gas[0], 
         bp->angular_momentum_gas[1], 
         bp->angular_momentum_gas[2],
         0.f,  /* specific angular momentum of the stars */
         0.f,  /* specific angular momentum of the stars */
         0.f,  /* specific angular momentum of the stars */
         bp->radiative_luminosity * props->conv_factor_energy_rate_to_cgs,
         bp->state,
         bp->f_accretion,
         bp->radiative_efficiency,
         bp->eddington_fraction,
         bp->gravitational_ngb_mass * props->mass_to_solar_mass);
}

/**
 * @brief Computes the (maximal) repositioning speed for a black hole.
 *
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
  bp->formation_gas_density = hydro_get_physical_density(p, cosmo);

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

  bp->state = BH_states_slim_disk;
  
  black_holes_mark_bpart_as_not_swallowed(&bp->merger_data);
}

/**
 * @brief Should this bh particle be doing any stars looping?
 *
 * @param bp The #bpart.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static int bh_stars_loop_is_active(
    const struct bpart* bp, const struct engine* e) {
  /* Active bhs never do the stars loop for the Obsidian model */
  return 0;
}

#endif /* SWIFT_OBSIDIAN_BLACK_HOLES_H */
