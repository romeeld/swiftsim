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

  /* Ignore wind in density computation */
  /*if (pj->feedback_data.decoupling_delay_time > 0.f) return;*/

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

  /* We found a neighbour! */
  si->feedback_data.ngb_N++;

  /* Add mass of pj to neighbour mass of si  */
  si->feedback_data.ngb_mass += mj;

  /* Update counter of total (integer) neighbours */
  si->feedback_data.num_ngbs++;

  /* Contribution to the star's surrounding gas density */
  si->feedback_data.ngb_rho += mj * wi;

  const float Zj = chemistry_get_total_metal_mass_fraction_for_feedback(pj);

  /* Contribution to the star's surrounding metallicity (metal mass fraction */
  si->feedback_data.ngb_Z += mj * Zj * wi;

  /* Add contribution of pj to normalisation of density weighted fraction
   * which determines how much mass to distribute to neighbouring
   * gas particles */
  si->feedback_data.enrichment_weight_inv += wi / rho;

}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep1(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  const struct spart *si, struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const integertime_t ti_current) {}

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

  /* If pj is already a wind particle, don't kick again */
  if (pj->feedback_data.decoupling_delay_time > 0.f) return;

  /* If more nearby particles already fill the kick list, then pj is not needed */
  if (r2 >= si->feedback_data.r2_gas_to_be_kicked[FEEDBACK_N_KICK_MAX-1]) return;

  /* Check if there is enough mass to launch */
  if (si->feedback_data.feedback_mass_to_launch < pj->mass) return;

  /* Check if there is enough energy to launch */
  const float energy = 0.5 * pj->mass * si->feedback_data.feedback_wind_velocity 
	  * si->feedback_data.feedback_wind_velocity;
  if (si->feedback_data.feedback_energy_reservoir < energy) return;

  /* Find location within kick list to insert this gas particle */
  int i;
  for (i = FEEDBACK_N_KICK_MAX; i > 0; i--) {
    if (r2 > si->feedback_data.r2_gas_to_be_kicked[i-1]) break;
    if (i < FEEDBACK_N_KICK_MAX) {
      si->feedback_data.id_gas_to_be_kicked[i] = 
          si->feedback_data.id_gas_to_be_kicked[i-1];
      si->feedback_data.r2_gas_to_be_kicked[i] = 
          si->feedback_data.r2_gas_to_be_kicked[i-1];
    }
  }

  /* Insert pj into kick list */
  si->feedback_data.id_gas_to_be_kicked[i] = pj->id;
  si->feedback_data.r2_gas_to_be_kicked[i] = r2;

  /* Remove mass and energy needed to launch this particle */
  si->feedback_data.feedback_mass_to_launch -= pj->mass;
  si->feedback_data.feedback_energy_reservoir -= energy;
  //message("KICKLIST: %lld %g %g %g %g", si->id, si->feedback_data.r2_gas_to_be_kicked[0], si->feedback_data.r2_gas_to_be_kicked[1], si->feedback_data.r2_gas_to_be_kicked[2], si->feedback_data.r2_gas_to_be_kicked[3]);
}

/**
 * @brief Kick and sometimes heat gas particle near a star, if star has enough mass
 * and energy for an ejection event.
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

  /* DO KINETIC FEEDBACK */
  /* Check if particle pj has been selected for a kick, and set alternating kick direction */
  int i;
  float kick_dir = 1.f;
  if (si->feedback_data.id_gas_to_be_kicked[0] % 2 == 1) kick_dir = -1.f;
  for (i = 0; i < FEEDBACK_N_KICK_MAX; i++) {
    if (pj->id == si->feedback_data.id_gas_to_be_kicked[i]) break;
    kick_dir *= -1.f;
  }
  /* This particle was not in the kick list, so nothing to do */
  if (i == FEEDBACK_N_KICK_MAX) return;

  /* Compute velocity and KE of wind event.
  * Note that pj->v_full = a^2 * dx/dt, with x the comoving
  * coordinate. Therefore, a physical kick, dv, gets translated into a
  * code velocity kick, a * dv */

  const double wind_velocity = 
      kick_dir * si->feedback_data.feedback_wind_velocity;

  //message("FEEDBACK %lld %lld E_sn=%g Ew=%g %g   M_ej=%g Mp=%g %g",si->id, pj->id, si->feedback_data.feedback_energy_reservoir, wind_energy, si->feedback_data.feedback_energy_reservoir/wind_energy, si->feedback_data.feedback_mass_to_launch, pj->mass, si->feedback_data.feedback_mass_to_launch/pj->mass);

  /* Set kick direction as v x a */
  const double dir[3] = {
    pj->gpart->a_grav[1] * pj->gpart->v_full[2] -
        pj->gpart->a_grav[2] * pj->gpart->v_full[1],
    pj->gpart->a_grav[2] * pj->gpart->v_full[0] -
        pj->gpart->a_grav[0] * pj->gpart->v_full[2],
    pj->gpart->a_grav[0] * pj->gpart->v_full[1] -
        pj->gpart->a_grav[1] * pj->gpart->v_full[0]
  };
  const double norm = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);

  /* No normalization, no wind (should basically never happen) */
  if (norm <= 0.) {
    warning("Normalization of wind direction is zero!\n(x, y, z) = (%g, %g, %g)",
            dir[0], dir[1], dir[2]);
    return;
  }

  const double prefactor = wind_velocity / norm;

  /* Do the kicks by updating the particle velocity. */
  pj->v_full[0] += dir[0] * prefactor;
  pj->v_full[1] += dir[1] * prefactor;
  pj->v_full[2] += dir[2] * prefactor;

  /* DO WIND HEATING */
  double u_new = fb_props->cold_wind_internal_energy;
  if (fb_props->cold_wind_internal_energy < fb_props->hot_wind_internal_energy) {
    double galaxy_stellar_mass =
          pj->gpart->fof_data.group_stellar_mass;
    if (galaxy_stellar_mass < fb_props->minimum_galaxy_stellar_mass) {
      galaxy_stellar_mass = fb_props->minimum_galaxy_stellar_mass;
    }
    const double galaxy_stellar_mass_Msun = galaxy_stellar_mass * 
	  fb_props->mass_to_solar_mass;

    /* Based on Pandya et al 2022 FIRE results */
    float pandya_slope = 0.f;
    if (galaxy_stellar_mass_Msun > 3.16e10) {
      pandya_slope = -2.1f;
    } else {
      pandya_slope = -0.1f;
    }

    /* 0.2511886 = pow(10., -0.6) */
    const double f_warm = 
        0.2511886 * pow(galaxy_stellar_mass_Msun / 3.16e10, pandya_slope);
    /* additional 10% removed for cold phase */
    const double hot_wind_fraction = max(0., 0.9 - f_warm);
    const double rand_for_hot = 
        random_unit_interval(pj->id, ti_current, random_number_stellar_feedback_3);
    const double rand_for_spread = 
        random_unit_interval(pj->id, ti_current, random_number_stellar_feedback);

    /* If selected, heat the particle */
    const double u_wind = 0.5 * wind_velocity * wind_velocity;
    if (rand_for_hot < hot_wind_fraction && 
            fb_props->hot_wind_internal_energy > u_wind) {
        u_new = (fb_props->hot_wind_internal_energy - u_wind) * 
                    (0.5 + rand_for_spread);
    }
  }

  /* Set the wind particle temperature */
  hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new);
  hydro_set_drifted_physical_internal_energy(pj, cosmo, NULL, u_new);

  /* For firehose model, set initial radius of stream */
  assert(si->feedback_data.firehose_radius_stream > 0.f);
  pj->chemistry_data.radius_stream = si->feedback_data.firehose_radius_stream;
  pj->chemistry_data.exchanged_mass = 0.f;

  /* FINISH UP FEEDBACK */
  /* Turn off any star formation in wind particle.
   * Record exp factor of when this particle was last ejected as -SFR. */
  pj->sf_data.SFR = -cosmo->a;

  /* Update the signal velocity of the particle based on the velocity kick */
  hydro_set_v_sig_based_on_velocity_kick(pj, cosmo, wind_velocity);

  /* Impose maximal viscosity */
  hydro_diffusive_feedback_reset(pj);

  /* Synchronize the particle on the timeline */
  timestep_sync_part(pj);

  /* Decouple the particles from the hydrodynamics */
  pj->feedback_data.decoupling_delay_time =
      fb_props->wind_decouple_time_factor *
      cosmology_get_time_since_big_bang(cosmo, cosmo->a);
  pj->chemistry_data.diffusion_coefficient = 0.f;

  /** Log the wind event.
   * z starid gasid dt M* vkick vkx vky vkz h x y z vx vy vz T rho v_sig tdec Ndec Z
   */
  const float length_convert = cosmo->a * fb_props->length_to_kpc;
  const float velocity_convert = cosmo->a_inv / fb_props->kms_to_internal;
  const float rho_convert = cosmo->a3_inv * fb_props->rho_to_n_cgs;
  const float u_convert =
      cosmo->a_factor_internal_energy / fb_props->temp_to_u_factor;
  printf("WIND_LOG %.5f %lld %g %g %g %g %lld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %g\n",
          cosmo->z,
          si->id,
	  si->feedback_data.feedback_mass_to_launch * fb_props->mass_to_solar_mass,
	  si->feedback_data.feedback_energy_reservoir * fb_props->energy_to_cgs / 1.e51,
	  1./si->birth_scale_factor - 1.,
          pj->gpart->fof_data.group_stellar_mass * fb_props->mass_to_solar_mass,
          pj->id,
          wind_velocity,
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

  /* If no mass to distribute, nothing to do */
  if (si->feedback_data.mass <= 0.f) return;

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

  /* Compute weighting for distributing feedback quantities */
  float Omega_frac = si->feedback_data.enrichment_weight * wi / rho_j;

  /* Never apply feedback if Omega_frac is bigger than or equal to unity */
  if (Omega_frac < 0. || Omega_frac > 1.01) {
    warning(
        "Invalid fraction of material to distribute for star ID=%lld "
        "Omega_frac=%e count since last enrich=%d enrichment_weight=%g wi=%g rho_j=%g",
        si->id, Omega_frac, si->count_since_last_enrichment,
	  si->feedback_data.enrichment_weight , wi , rho_j);
    if (Omega_frac > 1.0) Omega_frac = 1.0;
  }

  /* Update particle mass */
  const double current_mass = hydro_get_mass(pj);
  const double delta_mass = si->feedback_data.mass * Omega_frac;
  const double new_mass = current_mass + delta_mass;

  hydro_set_mass(pj, new_mass);

  /* Inverse of the new mass */
  const double new_mass_inv = 1. / new_mass;

  /* Update total metallicity */
  const double current_metal_mass_total =
      pj->chemistry_data.metal_mass_fraction_total * current_mass;
  const double delta_metal_mass_total =
      si->feedback_data.total_metal_mass * Omega_frac;
  const double new_metal_mass_total =
      current_metal_mass_total + delta_metal_mass_total;

  pj->chemistry_data.metal_mass_fraction_total =
      new_metal_mass_total * new_mass_inv;

  if (pj->chemistry_data.metal_mass_fraction_total < 0.f) {
    error("Stellar feedback led to negative metallicity!"
          "\tpid=%lld\n\tMnew=%g\n\tZ=%g\n\tOmega=%g\n"
          "\tMZtot=%g\n\tMcurr=%g",
          pj->id,
          new_mass,
          pj->chemistry_data.metal_mass_fraction_total,
          Omega_frac,
          si->feedback_data.total_metal_mass,
          current_mass);
    return;
  }
  /* Update mass fraction of each tracked element  */
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    const double current_metal_mass =
        pj->chemistry_data.metal_mass_fraction[elem] * current_mass;
    const double delta_metal_mass =
        si->feedback_data.metal_mass[elem] * Omega_frac;
    const double new_metal_mass = current_metal_mass + delta_metal_mass;

    pj->chemistry_data.metal_mass_fraction[elem] =
        new_metal_mass * new_mass_inv;
  }

  /* Compute kernel-smoothed contribution to number of SNe going off this timestep */
  pj->feedback_data.SNe_ThisTimeStep += 
      si->feedback_data.SNe_ThisTimeStep * Omega_frac;
  pj->feedback_data.SNe_ThisTimeStep = 
      fmax(pj->feedback_data.SNe_ThisTimeStep, 0.);

  /* Spread dust ejecta to gas */
  for (int elem = chemistry_element_He; elem < chemistry_element_count; elem++) {
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
  for (int elem = chemistry_element_He; elem < chemistry_element_count; elem++) {
    pj->cooling_data.dust_mass += pj->cooling_data.dust_mass_fraction[elem];
  }

  if (pj->cooling_data.dust_mass > 0.) {
    const double dust_mass_inv = 1. / pj->cooling_data.dust_mass;

    /* Divide by new dust mass to get the fractions */
    for (int elem = chemistry_element_He; elem < chemistry_element_count; elem++) {
      pj->cooling_data.dust_mass_fraction[elem] *= dust_mass_inv;
    }

    /* Check for inconsistency */
    if (pj->cooling_data.dust_mass > pj->mass) {
      for (int elem = chemistry_element_He; elem < chemistry_element_count; elem++) {
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
  if (pj->feedback_data.decoupling_delay_time > 0.f) return;

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
