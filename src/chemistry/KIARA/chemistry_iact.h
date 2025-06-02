/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_KIARA_CHEMISTRY_IACT_H
#define SWIFT_KIARA_CHEMISTRY_IACT_H

/**
 * @file KIARA/chemistry_iact.h
 * @brief Smooth metal interaction functions following the KIARA model.
 */

#include "timestep_sync_part.h"
#include <assert.h>

/**
 * @brief Sums ambient quantities for the firehose wind model
 *
 * This is called from runner_iact_chemistry, which is called during the density loop
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 */
__attribute__((always_inline)) INLINE static void firehose_compute_ambient_sym(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj) {

  /* Only decoupled particles need the ambient properties since they
   * are in the stream. */
  const int decoupled_i = pi->decoupled;
  const int decoupled_j = pj->decoupled;
  if (decoupled_i && decoupled_j) return;

  const float r = sqrtf(r2);

  /* Accumulate ambient neighbour quantities with an SPH gather operation */
  if (decoupled_i && !decoupled_j) {
    struct chemistry_part_data* chi = &pi->chemistry_data;
    const float hi_inv = 1. / hi;
    const float ui = r * hi_inv;
    const float mj = hydro_get_mass(pj);
    float wi;
    kernel_eval(ui, &wi);

#ifdef FIREHOSE_DEBUG_CHECKS
    if (!isfinite(mj * pj->u * wi)) {
      message("FIREHOSE_BAD pi=%lld ui=%g neighbour pj=%lld uj=%g  mj=%g  hi=%g"
              "wi=%g\n",
              pi->id, 
              pi->u,
              pj->id, 
              pj->u, 
              mj, 
              hi, 
              wi);
    }
#endif

    chi->u_ambient += mj * pj->u * wi;
    chi->rho_ambient += mj * wi;
    chi->w_ambient += wi;
    chi->v_sig_ambient += mj * pj->viscosity.v_sig * pj->viscosity.v_sig * pj->viscosity.v_sig * wi;
  }

  if (!decoupled_i && decoupled_j) {
    struct chemistry_part_data* chj = &pj->chemistry_data;
    const float hj_inv = 1. / hj;
    const float uj = r * hj_inv;
    const float mi = hydro_get_mass(pi);
    float wj;
    kernel_eval(uj, &wj);

#ifdef FIREHOSE_DEBUG_CHECKS
    if (!isfinite(mi * pi->u * wj)) {
      message("FIREHOSE_BAD pj=%lld uj=%g neighbour pi=%lld ui=%g  mi=%g  hj=%g"
              "wj=%g\n",
              pj->id,
              pj->u,
              pi->id, 
              pi->u,
              mi, 
              hj, 
              wj);
    }
#endif

    chj->u_ambient += mi * pi->u * wj;
    chj->rho_ambient += mi * wj;
    chj->w_ambient += wj;
    chj->v_sig_ambient += mi * pi->viscosity.v_sig * pi->viscosity.v_sig * pi->viscosity.v_sig * wj;
  }
}

/**
 * @brief Sums ambient quantities for the firehose wind model
 *
 * This is called from runner_iact_chemistry, which is called during the density loop
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 */
__attribute__((always_inline)) INLINE static void 
firehose_compute_ambient_nonsym(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj) {

  /* Only decoupled winds need the ambient quantities computed because
   * they are in the stream. */
  const int decoupled_i = pi->decoupled;
  const int decoupled_j = pj->decoupled;

  /* A wind particle cannot be in the ambient medium */
  if (!decoupled_i || decoupled_j) return;

  struct chemistry_part_data* chi = &pi->chemistry_data;

  /* Do accumulation of ambient quantities */

  /* Compute the kernel function for pi */
  const float r = sqrtf(r2);
  const float h_inv = 1. / hi;
  const float ui = r * h_inv;
  const float mj = hydro_get_mass(pj);
  float wi;
  kernel_eval(ui, &wi);

#ifdef FIREHOSE_DEBUG_CHECKS
  if (!isfinite(mj * pj->u * wi)) {
    message("FIREHOSE_BAD pi=%lld ui=%g neighbour pj=%lld uj=%g  mj=%g  hi=%g"
            "wi=%g\n",
            pi->id, 
            pi->u,
            pj->id, 
            pj->u, 
            mj, 
            hi, 
            wi);
  }
#endif

  /* Accumulate ambient neighbour quantities with an SPH gather operation */
  chi->u_ambient += mj * pj->u * wi;
  chi->rho_ambient += mj * wi;
  chi->w_ambient += wi;
  chi->v_sig_ambient += mj * pj->viscosity.v_sig * pj->viscosity.v_sig * pj->viscosity.v_sig * wi;
}

/**
 * @brief do chemistry computation after the runner_iact_density (symmetric
 * version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_chemistry(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  /* If in wind mode, compute ambient quantities for firehose wind diffusion */
  firehose_compute_ambient_sym(r2, dx, hi, hj, pi, pj);

  /* Do not need diffusion properties for wind particles */
  if (pi->decoupled || pj->decoupled) return;
          
  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  float wi, wi_dx;
  float wj, wj_dx;

  /* Get the masses. */
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  /* Get r */
  const float r = sqrtf(r2);
  const float r_inv = 1.f / r;

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute the kernel function for pj */
  const float uj = r / hj;
  kernel_deval(uj, &wj, &wj_dx);

  const float wi_dr = wi_dx * r_inv;
  const float mj_wi_dr = mj * wi_dr;

  const float wj_dr = wj_dx * r_inv;
  const float mi_wj_dr = mi * wj_dr;

  /* dx[i] is from i -> j, so should dv_ij be from i -> j */
  const float dv_ij[3] = {
    pi->v[0] - pj->v[0],
    pi->v[1] - pj->v[1],
    pi->v[2] - pj->v[2]
  };

  /* Compute the shear tensor, i = spatial direction */
  for (int i = 0; i < 3; i++) {
    const float dxi_mj_wi_dr = dx[i] * mj_wi_dr;
    const float dxi_mi_wj_dr = dx[i] * mi_wj_dr;

    chi->shear_tensor[0][i] += dv_ij[0] * dxi_mj_wi_dr;
    chi->shear_tensor[1][i] += dv_ij[1] * dxi_mj_wi_dr;
    chi->shear_tensor[2][i] += dv_ij[2] * dxi_mj_wi_dr;

    /* Sign must be positive since dx is always i -> j */
    chj->shear_tensor[0][i] += dv_ij[0] * dxi_mi_wj_dr;
    chj->shear_tensor[1][i] += dv_ij[1] * dxi_mi_wj_dr;
    chj->shear_tensor[2][i] += dv_ij[2] * dxi_mi_wj_dr;
  }
}

/**
 * @brief do chemistry computation after the runner_iact_density (non symmetric
 * version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_chemistry(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H) {

  /* If in wind mode, compute ambient quantities for firehose wind diffusion */
  firehose_compute_ambient_nonsym(r2, dx, hi, hj, pi, pj);

  /* Do not need diffusion properties for wind particles */
  if (pi->decoupled || pj->decoupled) return;

  struct chemistry_part_data *chi = &pi->chemistry_data;

  float wi, wi_dx;

  /* Get the masses. */
  const float mj = hydro_get_mass(pj);

  /* Get r */
  const float r = sqrtf(r2);
  const float r_inv = 1.f / r;

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_deval(ui, &wi, &wi_dx);

  const float wi_dr = wi_dx * r_inv;
  const float mj_wi_dr = mj * wi_dr;

  const float dv_ij[3] = {
    pi->v[0] - pj->v[0],
    pi->v[1] - pj->v[1],
    pi->v[2] - pj->v[2]
  };

  /* Compute the shear tensor */
  for (int i = 0; i < 3; i++) {
    const float dxi_mj_wi_dr = dx[i] * mj_wi_dr;
    chi->shear_tensor[0][i] += dv_ij[0] * dxi_mj_wi_dr;
    chi->shear_tensor[1][i] += dv_ij[1] * dxi_mj_wi_dr;
    chi->shear_tensor[2][i] += dv_ij[2] * dxi_mj_wi_dr;
  }
}

/**
 * @brief Computes the mass exchanged between the firehose stream 
 * and the ambient medium. Note that either i or j could be the stream 
 * particle, with j or i being ambient.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Wind particle.
 * @param pj Ambient particle.
 * @param time_base The time base used to convert integer to float time.
 * @param ti_current Current integer time for seeding random number generator   
 * @param phys_const Physical constants
 * @param cd #chemistry_global_data containing chemistry information.
 * @param v2 velocity difference squared between i and j.
 *
 */
__attribute__((always_inline)) INLINE static float 
firehose_compute_mass_exchange(
    const float r2, const float dx[3], const float hi, const float hj,
    const struct part *pi, const struct part *pj,
    const float time_base, const integertime_t ti_current,
    const struct phys_const* phys_const, const struct chemistry_global_data* cd, 
    float *v2,
    const struct cosmology *cosmo) {

  const int i_stream = pi->decoupled;
  const int j_stream = pj->decoupled;

  /* Both particles cannot be in the stream. The one with >0 delay time is 
   * the stream particle. At least one particle must be in the stream. */
  if (i_stream && j_stream) return 0.f;
  if (!i_stream && !j_stream) return 0.f;

  const double dt_i = get_timestep(pi->time_bin, time_base);
  const double dt_j = get_timestep(pj->time_bin, time_base);

  /* Must be a timestep for both */
  if (dt_i <= 0. || dt_j <= 0.) return 0.f;

  /* For stream particle, make sure the stream radius > 0 */
  if (i_stream && pi->chemistry_data.radius_stream <= 0.f) return 0.f;
  if (j_stream && pj->chemistry_data.radius_stream <= 0.f) return 0.f;
    
  /* Stream particle cannot be an AGN jet particle */
  if (i_stream && 
          pi->feedback_data.number_of_times_decoupled >= 10000) return 0.f;
  if (j_stream && 
          pj->feedback_data.number_of_times_decoupled >= 10000) return 0.f;

  /* Compute the velocity of the stream relative to ambient gas.
   * Order does not matter here because it is symmetric. */
  *v2 = 0.f;
  for (int i = 0; i < 3; i++) {
    *v2 += (pj->v_full[i] - pi->v_full[i]) * (pj->v_full[i] - pi->v_full[i]);
  }

  /* Don't apply above some velocity to avoid jets */
  const float v2_phys = *v2 * cosmo->a_inv * cosmo->a_inv;
  const float max_v2_phys = 
      cd->firehose_max_velocity * cd->firehose_max_velocity;
  if (v2_phys > max_v2_phys) return 0.f;

  /* Get mixing layer cooling time, which is negative if cooling */
  const float mixing_layer_time_i = pi->cooling_data.mixing_layer_cool_time;
  const float mixing_layer_time_j = pj->cooling_data.mixing_layer_cool_time;

  const float r = sqrtf(r2);
  const float gamma_gamma_minus_1 = hydro_gamma * hydro_gamma_minus_one;

  /* Set these based on which is in the stream */
  float mi;
  float wi;
  float sum_wi;
  float chi;
  float c_stream;
  float c_amb;
  float radius_stream;
  float mixing_layer_time;
  double dt = 0.;

  if (i_stream) {
    dt = dt_i;
    mi = hydro_get_mass(pi);
    const float h_inv = 1. / hi;
    const float ui = r * h_inv;
    kernel_eval(ui, &wi);

    /* Normalization of the ambient medium */
    sum_wi = pi->chemistry_data.w_ambient;

    /* Compute thermal energy ratio for stream and ambient */
    chi = pi->chemistry_data.u_ambient / pi->u;
    c_stream = sqrtf(pi->u * gamma_gamma_minus_1);
    c_amb = sqrtf(pi->chemistry_data.u_ambient * gamma_gamma_minus_1);
    radius_stream = pi->chemistry_data.radius_stream;
    mixing_layer_time = mixing_layer_time_i;
  }
  
  if (j_stream) { /* j must be the stream here */
    dt = dt_j;
    mi = hydro_get_mass(pj);
    const float h_inv = 1. / hj;
    const float ui = r * h_inv;
    kernel_eval(ui, &wi);

    /* Normalization of the ambient medium */
    sum_wi = pj->chemistry_data.w_ambient;

    /* Compute thermal energy ratio for stream and ambient */
    chi = pj->chemistry_data.u_ambient / pj->u;
    c_stream = sqrtf(pj->u * gamma_gamma_minus_1);
    c_amb = sqrtf(pj->chemistry_data.u_ambient * gamma_gamma_minus_1);
    radius_stream = pj->chemistry_data.radius_stream;
    mixing_layer_time = mixing_layer_time_j;
  }

  /* This should never happen. */
  if (dt == 0.) return 0.f;

  double dm = 0.;
  double t_cool_mix = 1.e10 * dt;

  /* Mass change is growth due to cooling minus loss due to shearing, 
     kernel-weighted. */
  const float v_stream = sqrtf(*v2);
  const float Mach = v_stream / (c_stream + c_amb);
  const float alpha = 0.21f * (0.8f * exp(-3.f * Mach * Mach) + 0.2f);

  t_cool_mix = 
      (mixing_layer_time < 0.f) ? fabs(mixing_layer_time) : t_cool_mix;

  const double t_shear = radius_stream / (alpha * v_stream);
  const double t_sound = 2.f * radius_stream / c_stream;

  const double delta_growth = 
      (4. / (chi * t_sound)) * pow(t_cool_mix / t_sound, -0.25) * dt;

  double delta_shear = 0.;
  if (t_shear < t_cool_mix) delta_shear = (1. - exp(-dt / t_shear));

  dm = mi * (delta_growth - delta_shear) * (wi / sum_wi);

  /* If stream is growing, don't mix */
  if (dm > 0.f) dm = 0.f;

#ifdef FIREHOSE_DEBUG_CHECKS
  if (dm < 0.f && i_stream && pj->cooling_data.subgrid_temp > 0.f) {
    message("FIREHOSE: z=%g %lld %lld m=%g nHamb=%g rhoamb/rhoi=%g rhoamb/rhoj=%g"
            " Tamb=%g Tj/Tamb=%g cstr/camb=%g M=%g r=%g grow=%g shear=%g tshear=%g"
            " tcool=%g fexch=%g", 
	    cosmo->z,
            pi->id, 
            pj->id, 
            pi->mass, 
            pi->chemistry_data.rho_ambient * cosmo->a3_inv * cd->rho_to_n_cgs,
            pi->chemistry_data.rho_ambient / pi->rho, 
            pi->chemistry_data.rho_ambient / pj->rho, 
            pi->chemistry_data.u_ambient * cosmo->a_factor_internal_energy / cd->temp_to_u_factor, 
            pj->u/pi->chemistry_data.u_ambient, 
	    c_stream / c_amb,
	    Mach,
	    radius_stream * cd->length_to_kpc * cosmo->a,
            delta_growth, 
            delta_shear, 
            t_shear, 
            pi->cooling_data.mixing_layer_cool_time, 
            dm/pi->mass);
  }
#endif

  return dm;
}


/**
 * @brief Computes the particle interaction via the firehose stream model
 *
 * This is called from runner_iact_diffusion, which is called during the force 
 * loop
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Wind particle.
 * @param pj Ambient particle.
 * @param time_base The time base used to convert integer to float time.
 * @param ti_current Current integer time for seeding random number generator.
 * @param phys_const Physical constants
 * @param cd #chemistry_global_data containing chemistry information.
 *
 */
__attribute__((always_inline)) INLINE static void firehose_evolve_particle_sym(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *pi, struct part *pj,
    struct xpart *xpi, struct xpart *xpj,
    const float time_base, const integertime_t ti_current,
    const struct phys_const* phys_const, const struct chemistry_global_data* cd,
    const struct cosmology *cosmo) {

  /* Both particles must be within each others smoothing lengths */
  const float Hi = kernel_gamma * hi;
  const float Hj = kernel_gamma * hj;
  const int r_in_Hi = (r2 < Hi * Hi) ? 1 : 0;
  const int r_in_Hj = (r2 < Hj * Hj) ? 1 : 0;
  if (!r_in_Hi || !r_in_Hj) return;

  const int i_stream = pi->decoupled;
  const int j_stream = pj->decoupled;

  /* Only one of the particles must be in the stream. */
  if (i_stream && j_stream) return;
  if (!i_stream && !j_stream) return;

  if (r2 <= 0.f) return;

  struct chemistry_part_data* chi = &pi->chemistry_data;
  struct chemistry_part_data* chj = &pj->chemistry_data;

  /* Compute the amount of mass mixed between stream particle and ambient gas */
  float v2 = 0.f;
  const float dm = firehose_compute_mass_exchange(r2, dx, hi, hj, pi, pj, 
                                                  time_base, ti_current, 
                                                  phys_const, cd, &v2,
                                                  cosmo);
  float delta_m = fabs(dm);
  if (delta_m <= 0.f) return; 

  /* Limit mass exchange to some fraction of particles' mass */
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  const float max_fmix_this_step = cd->firehose_max_fmix_per_step;
  if (delta_m > max_fmix_this_step * mi) delta_m = max_fmix_this_step * mi;
  if (delta_m > max_fmix_this_step * mj) delta_m = max_fmix_this_step * mj;

  /* Track amount of gas mixed in stream particle */
  if (i_stream) chi->exchanged_mass += delta_m;
  if (j_stream) chj->exchanged_mass += delta_m;

  chi->dm += delta_m;
  chj->dm += delta_m;

  /* set weights for averaging i and j */
  const float pii_weight = (mi - delta_m) / mi;
  const float pij_weight = delta_m / mi;
  const float pji_weight = delta_m / mj;
  const float pjj_weight = (mj - delta_m) / mj;

  /* Mixing is negligibly small, avoid underflows */
  if (pij_weight < 1.e-10f || pji_weight < 1.e-10f) return;

  const float wt_ii = pii_weight * mi;
  const float wt_ij = pij_weight * mi;
  const float wt_jj = pjj_weight * mj;
  const float wt_ji = pji_weight * mj;

  /* Spread dust masses between particles */
  const float pi_dust_mass = pi->cooling_data.dust_mass;
  const float pj_dust_mass = pj->cooling_data.dust_mass;

  const float dust_wt_ii = pii_weight * pi_dust_mass;
  const float dust_wt_ij = pij_weight * pi_dust_mass;
  const float dust_wt_ji = pji_weight * pj_dust_mass;
  const float dust_wt_jj = pjj_weight * pj_dust_mass;

  const float new_pi_dust_mass = dust_wt_ii + dust_wt_ij;
  const float new_pj_dust_mass = dust_wt_ji + dust_wt_jj;
  chi->dm_dust += new_pi_dust_mass - pi_dust_mass;
  chj->dm_dust += new_pj_dust_mass - pj_dust_mass;

  /* Spread individual dust elements */
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    /* Exchange metals */
    const float term_ii = wt_ii * chi->metal_mass_fraction[elem];
    const float term_ij = wt_ij * chj->metal_mass_fraction[elem];
    const float term_jj = wt_jj * chj->metal_mass_fraction[elem];
    const float term_ji = wt_ji * chi->metal_mass_fraction[elem];

    const float old_pi_Z_mass = mi * chi->metal_mass_fraction[elem];
    const float new_pi_Z_mass = term_ii + term_ij;
    const float old_pj_Z_mass = mj * chj->metal_mass_fraction[elem];
    const float new_pj_Z_mass = term_jj + term_ji;

    chi->dm_Z[elem] += new_pi_Z_mass - old_pi_Z_mass;
    chj->dm_Z[elem] += new_pj_Z_mass - old_pj_Z_mass;

    /* Exchange dust metals */
    const float pi_dust_frac = pi->cooling_data.dust_mass_fraction[elem];
    const float pj_dust_frac = pj->cooling_data.dust_mass_fraction[elem];

    /* Particle i */
    const float old_pi_dust_mass_Z = pi_dust_frac * pi_dust_mass;

    const float dust_term_ii = dust_wt_ii * pi_dust_frac;
    const float dust_term_ij = dust_wt_ij * pj_dust_frac;
    const float new_pi_dust_mass_Z = dust_term_ii + dust_term_ij;

    chi->dm_dust_Z[elem] += new_pi_dust_mass_Z - old_pi_dust_mass_Z;

    /* Particle j */
    const float old_pj_dust_mass_Z = pj_dust_frac * pj_dust_mass;

    const float dust_term_jj = dust_wt_jj * pj_dust_frac;
    const float dust_term_ji = dust_wt_ji * pi_dust_frac;
    const float new_pj_dust_mass_Z = dust_term_jj + dust_term_ji;

    chj->dm_dust_Z[elem] += new_pj_dust_mass_Z - old_pj_dust_mass_Z;
  }

  /* Update particles' internal energy per unit mass */
  const float old_pi_u = pi->u;
  const float old_pj_u = pj->u;

  float new_pi_u = (wt_ii * old_pi_u + wt_ij * old_pj_u) / mi;
  float new_pj_u = (wt_ji * old_pi_u + wt_jj * old_pj_u) / mj;

  /* Accumulate the velocity exchange due to the mass, conserving the 
   * momentum */
  float new_v2 = 0.f;
  for (int i = 0; i < 3; i++) {
    const float new_pi_v_full_i = 
        (wt_ii * pi->v_full[i] + wt_ij * pj->v_full[i]) / mi;
    /* Keep track of the final new velocity */
    chi->dv[i] += new_pi_v_full_i - pi->v_full[i];
    
    const float new_pj_v_full_i = 
        (wt_ji * pi->v_full[i] + wt_jj * pj->v_full[i]) / mj;
    chj->dv[i] += new_pj_v_full_i - pj->v_full[i];

    const float dv_i = new_pi_v_full_i - new_pj_v_full_i;
    new_v2 += dv_i * dv_i;
  }

  /* 4) Split excess energy between stream and ambient particle */
  float delta_KE = 0.5f * delta_m * (v2 - new_v2);
  const float min_KE = min(mi * new_pi_u, mj * new_pj_u);

  /* Limit to minimum of the two particles */
  if (delta_KE > min_KE) delta_KE = min_KE;

  const float dKE_split = 0.5f * delta_KE;

  /* Add in the delta KE difference */
  new_pi_u += dKE_split / mi;
  new_pj_u += dKE_split / mj;

  /* Accumulate the changes in energy in all interactions */
  chi->du += new_pi_u - old_pi_u;
  chj->du += new_pj_u - old_pj_u;

#ifdef FIREHOSE_DEBUG_CHECKS
  message("FIREHOSE_EXCHANGE: pi=%lld pj=%lld si=%d sj=%d npi_u=%g npj_u=%g "
          "opi_u=%g opj_u=%g dKE=%g nv2=%g ov2=%g",
          pi->id,
          pj->id,
          i_stream,
          j_stream,
          new_pi_u,
          new_pj_u,
          old_pi_u,
          old_pj_u,
          delta_KE,
          new_v2,
          v2);
#endif

}

/**
 * @brief do metal diffusion computation in the <FORCE LOOP>
 * (symmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param time_base The time base used to convert integer to float time.
 * @param ti_current The current time (in integer)
 * @param cosmo The cosmology information.
 * @param with_cosmology Are we running with cosmology?
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_diffusion(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, 
    struct xpart *restrict xpi, struct xpart *restrict xpj,
    const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology, 
    const struct phys_const* phys_const, const struct chemistry_global_data *cd) {

  if (pi->decoupled || pj->decoupled) {
    if (cd->use_firehose_wind_model) {
      /* If in wind mode, do firehose wind diffusion */
      firehose_evolve_particle_sym(r2, dx, hi, hj, pi, pj, xpi, xpj, time_base, 
                                    t_current, phys_const, cd, cosmo);
    }

    return;
  }

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  /* No need to diffuse if both particles are not diffusing. */
  if (chj->diffusion_coefficient > 0.f && chi->diffusion_coefficient > 0.f) {

    /* Get mass */
    const float mj = hydro_get_mass(pj);
    const float mi = hydro_get_mass(pi);
    const float rhoj = hydro_get_physical_density(pj, cosmo);
    const float rhoi = hydro_get_physical_density(pi, cosmo);

    float wi, wj, dwi_dx, dwj_dx;

    /* Get r */
    const float r = sqrtf(r2);

    /* part j */
    /* Get the kernel for hj */
    const float hj_inv = 1.0f / hj;

    /* Compute the kernel function for pj */
    const float xj = r * hj_inv;
    kernel_deval(xj, &wj, &dwj_dx);

    /* part i */
    /* Get the kernel for hi */
    const float hi_inv = 1.0f / hi;

    /* Compute the kernel function for pi */
    const float xi = r * hi_inv;
    kernel_deval(xi, &wi, &dwi_dx);

    /* Get 1/r */
    const float r_inv = r > 0.f ? 1.f / r : 0.f;

    const float wi_dr = dwi_dx * r_inv;
    const float wj_dr = dwj_dx * r_inv;

    const float mj_dw_r = mj * wi_dr;
    const float mi_dw_r = mi * wj_dr;

    const float rhoij_inv = 1.f / (rhoi * rhoj);

    /**
     * Compute the diffusion following Eq. 2.14
     * from Monaghan, Huppert, & Worster (2006).
     */
    float coef = 4.f * chi->diffusion_coefficient * chj->diffusion_coefficient;
    coef /= chi->diffusion_coefficient + chj->diffusion_coefficient;

    const float coef_i = coef * mj_dw_r * rhoij_inv;
    const float coef_j = coef * mi_dw_r * rhoij_inv;

    /* Compute the time derivative of metals due to diffusion */
    const float dZ_ij_tot = 
        chi->metal_mass_fraction_total - chj->metal_mass_fraction_total;
    chi->dZ_dt_total += coef_i * dZ_ij_tot;
    chj->dZ_dt_total -= coef_j * dZ_ij_tot;

    for (int elem = 0; elem < chemistry_element_count; elem++) {
      const float dZ_ij = 
          chi->metal_mass_fraction[elem] - chj->metal_mass_fraction[elem];
      chi->dZ_dt[elem] += coef_i * dZ_ij;
      chj->dZ_dt[elem] -= coef_j * dZ_ij;
    }
  }
}

/**
 * @brief do metal diffusion computation in the <FORCE LOOP>
 * (nonsymmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle. (not updated)
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param time_base The time base used to convert integer to float time.
 * @param ti_current The current time (in integer)
 * @param cosmo The #cosmology.
 * @param with_cosmology Are we running with cosmology?
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_diffusion(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct phys_const* phys_const, 
    const struct chemistry_global_data *cd) {

  /* In nonsym case, two cases: depending on whether i is stream or ambient */
  if (pi->decoupled || pj->decoupled) return;

  struct chemistry_part_data *chi = &pi->chemistry_data;
  const struct chemistry_part_data *chj = &pj->chemistry_data;

  if (chj->diffusion_coefficient > 0 && chi->diffusion_coefficient > 0) {

    /* Get mass */
    const float mj = hydro_get_mass(pj);
    const float rhoj = hydro_get_physical_density(pj, cosmo);
    const float rhoi = hydro_get_physical_density(pi, cosmo);

    float wi, dwi_dx;

    /* Get r */
    const float r = sqrtf(r2);

    /* part i */
    /* Get the kernel for hi */
    const float hi_inv = 1.0f / hi;

    /* Compute the kernel function for pi */
    const float xi = r * hi_inv;
    kernel_deval(xi, &wi, &dwi_dx);

    /* Get 1/r */
    const float r_inv = 1.f / sqrtf(r2);
    const float wi_dr = dwi_dx * r_inv;

    const float mj_dw_r = mj * wi_dr;

    const float rhoij_inv = 1.f / (rhoi * rhoj);

    /**
     * Compute the diffusion following Eq. 2.14
     * from Monaghan, Huppert, & Worster (2006).
     */
    float coef = 4.f * chi->diffusion_coefficient * chj->diffusion_coefficient;
    coef /= chi->diffusion_coefficient + chj->diffusion_coefficient;

    const float coef_i = coef * mj_dw_r * rhoij_inv;

    /* Compute the time derivative */
    const float dZ_ij_tot = 
        chi->metal_mass_fraction_total - chj->metal_mass_fraction_total;
    chi->dZ_dt_total += coef_i * dZ_ij_tot;

    for (int elem = 0; elem < chemistry_element_count; elem++) {
      const float dZ_ij = 
          chi->metal_mass_fraction[elem] - chj->metal_mass_fraction[elem];
      chi->dZ_dt[elem] += coef_i * dZ_ij;
    }
  }
}


#endif /* SWIFT_KIARA_CHEMISTRY_IACT_H */
