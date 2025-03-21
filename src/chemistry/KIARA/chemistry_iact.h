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

  if (pi->feedback_data.decoupling_delay_time <= 0.f 
        && pj->feedback_data.decoupling_delay_time <= 0.f) return;

  struct chemistry_part_data* chi = &pi->chemistry_data;
  struct chemistry_part_data* chj = &pj->chemistry_data;

  /* Do accumulation of ambient quantities */
  if (r2 > 0.f) {
    /* Compute the kernel function for pi */
    const float r = sqrtf(r2);
    const float hi_inv = 1. / hi;
    const float hj_inv = 1. / hj;
    const float ui = r * hi_inv;
    const float uj = r * hj_inv;
    float wi, wj;
    kernel_eval(ui, &wi);
    kernel_eval(uj, &wj);
    const float mj_wi = pj->mass * wi;
    const float mi_wj = pi->mass * wj;

    /* Accumulate ambient neighbour quantities with an SPH gather operation */
    if (wi > 0.f) {
      chi->u_ambient += pj->u * mj_wi;
      chi->rho_ambient += mj_wi * pow_dimension(hi_inv);
      chi->w_ambient += mj_wi;

      chj->u_ambient += pi->u * mi_wj;
      chj->rho_ambient += mi_wj * pow_dimension(hj_inv);
      chj->w_ambient += mi_wj;
    }
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
__attribute__((always_inline)) INLINE static void firehose_compute_ambient_nonsym(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj) {

  if (pi->feedback_data.decoupling_delay_time <= 0.f) return;

  struct chemistry_part_data* chi = &pi->chemistry_data;

  /* Do accumulation of ambient quantities */
  if (r2 > 0.f) {
    /* Compute the kernel function for pi */
    const float r = sqrtf(r2);
    const float h_inv = 1. / hi;
    const float ui = r * h_inv;
    float wi;
    kernel_eval(ui, &wi);
    const float mj_wi = pj->mass * wi;

    /* Accumulate ambient neighbour quantities with an SPH gather operation */
    if (wi > 0.f) {
      chi->u_ambient += pj->u * mj_wi;
      chi->rho_ambient += mj_wi * pow_dimension(h_inv);
      chi->w_ambient += mj_wi;
    }
  }
}

/**
 * @file KIARA/chemistry_iact.h
 * @brief Smooth metal interaction functions following the KIARA model.
 */

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
 * @brief Computes the mass exchanged between the firehose stream and the ambient medium.
 * Note that either i or j could be the stream particle, with j or i being ambient.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Wind particle.
 * @param pj Ambient particle.
 * @param time_base The time base used to convert integer to float time.
 * @param ti_current Current integer time used value for seeding random number generator   
 * @param phys_const Physical constants
 * @param cd #chemistry_global_data containing chemistry information.
 * @param v2 velocity difference squared between i and j.
 *
 */
__attribute__((always_inline)) INLINE static float firehose_compute_mass_exchange(
    const float r2, const float dx[3], const float hi, const float hj,
    const struct part *pi, const struct part *pj,
    const float time_base, const integertime_t ti_current,
    const struct phys_const* phys_const, const struct chemistry_global_data* cd, 
    float *v2,
    const struct cosmology *cosmo) {

  /* never do diffusion between a cooling shutoff particle */
  if (pi->feedback_data.cooling_shutoff_delay_time > 0.f ||
        pj->feedback_data.cooling_shutoff_delay_time > 0.f) return 0.f;

  const float decouple_time_i = pi->feedback_data.decoupling_delay_time;
  const float decouple_time_j = pj->feedback_data.decoupling_delay_time;

  /* Both particles cannot be in the stream.  
     The one with >0 delay time is the stream particle */
  if (decouple_time_i * decouple_time_j > 0.f) return 0.;

  /* For stream particle, make sure the stream radius > 0 */
  if (decouple_time_i > 0.f && 
        pi->chemistry_data.radius_stream <= 0.f) return 0.;
  if (decouple_time_j > 0.f && 
        pj->chemistry_data.radius_stream <= 0.f) return 0.;
    
  /* Compute the kernel function */
  const float r = sqrtf(r2);
  const float h_inv = 1. / hi;
  const float ui = r * h_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* too far out to interact */
  if (wi <= 0.f) return 0.f;

  /* Get the timestep */ 
  const float dt = get_timestep(pi->time_bin, time_base);
  
  /* zero time step so no mixing */
  if (dt <= 0.f) return 0.f;

  /* Compute the velocity of the stream relative to ambient gas */
  *v2 = 0.f;
  for (int i = 0; i < 3; i++) {
    *v2 += (pi->v_full[i] - pj->v_full[i]) * (pi->v_full[i] - pj->v_full[i]);
  }

  /* Don't apply above some velocity to avoid jets */
  const float v2_phys = *v2 * cosmo->a_inv * cosmo->a_inv;
  const float max_v2_phys = cd->firehose_max_velocity * cd->firehose_max_velocity;
  if (v2_phys > max_v2_phys) return 0.f;

  /* Compute thermal energy ratio for stream and ambient */
  const float chi = pi->chemistry_data.u_ambient / pi->u;
  const float v_stream = sqrtf(*v2);
  const float c_stream = sqrtf(pi->u * hydro_gamma * hydro_gamma_minus_one);
  const float c_amb = 
      sqrtf(pi->chemistry_data.u_ambient * hydro_gamma * hydro_gamma_minus_one);
  const float Mach = v_stream / (c_stream + c_amb);
  const float alpha = 0.21f * (0.8f * exp(-3.f * Mach * Mach) + 0.2f);

  /* Define timescales for streams */

  /* Get mixing layer cooling time, which is negative if cooling */
  const float mixing_layer_time_i = pi->cooling_data.mixing_layer_cool_time;
  const float mixing_layer_time_j = pj->cooling_data.mixing_layer_cool_time;

  /* some large number that will result in no mixing */
  float t_cool_mix = 1.e10f * dt;
  if (decouple_time_i > 0.f && mixing_layer_time_i < 0.f) {
    t_cool_mix = fabs(mixing_layer_time_i);
  }

  if (decouple_time_j > 0.f && mixing_layer_time_j < 0.f) {
    t_cool_mix = fabs(mixing_layer_time_j);
  }

  double t_shear = pi->chemistry_data.radius_stream / (alpha * v_stream);
  float t_sound = 2.f * pi->chemistry_data.radius_stream / c_stream;
  if (decouple_time_j > 0.f) {
    t_shear = pj->chemistry_data.radius_stream / (alpha * v_stream);
    t_sound = 2.f * pj->chemistry_data.radius_stream / c_stream;
  }

  /* Mass change is growth due to cooling minus loss due to shearing, 
     kernel-weighted */
  float dm = 0.f;
  float delta_shear = 0.f;
  float delta_growth = 0.f;
  if (t_shear < t_cool_mix) delta_shear = (1.f - exp(-dt / t_shear));
  if (decouple_time_i > 0.f && mixing_layer_time_i < 0.f) {
    delta_growth = 
        (4.f / (chi * t_sound)) * pow(t_cool_mix / t_sound, -0.25f) * dt;
  }

  if (decouple_time_j > 0.f && mixing_layer_time_j < 0.f) {
    delta_growth = 
        (4.f / (chi * t_sound)) * pow(t_cool_mix / t_sound, -0.25f) * dt;
  }

  dm = wi * pi->mass * (delta_growth - delta_shear);

  /* Limit amount of mixing per neighbor,
     with ~50 neighbours, this limits total loss/gain to a particle's mass 
     in a single step. */
  const float fmix_max = 0.02f;
  if (dm > fmix_max * pi->mass) dm = fmix_max * pi->mass;
  if (dm < -fmix_max * pi->mass) dm = -fmix_max * pi->mass;

  /* If stream is growing, don't mix */
  if (dm > 0.f) dm = 0.f;

#ifdef FIREHOSE_DEBUG_CHECKS
  if (dm < 0.f) {
    message("FIREHOSE: %lld %lld m=%g rhoi=%g rhoamb=%g rhoj=%g"
            " ui=%g uamb=%g uj=%g Ti/Tamb=%g grow=%g shear=%g tshear=%g tcmix=%g"
            " tcool=%g fexch=%g", 
            pi->id, 
            pj->id, 
            pi->mass, 
            pi->rho, 
            pi->chemistry_data.rho_ambient, 
            pj->rho, 
            pi->u, 
            pi->chemistry_data.u_ambient, 
            pj->u, 
            pi->u/pi->chemistry_data.u_ambient, 
            delta_growth, 
            delta_shear, 
            tshear, 
            t_cool_mix/tshear, 
            pi->cooling_data.mixing_layer_cool_time, 
            dm/pi->mass);
  }
#endif

  return dm;
}


/**
 * @brief Check recoupling criterion for firehose stream particle .
 * Returns negative value if it should recouple.
 * Actual recoupling is done in feedback.h.
 *
 * @param pi Wind particle (not updated).
 * @param Mach Stream Mach number vs ambient
 * @param r_stream Current radius of stream 
 * @param cd #chemistry_global_data containing chemistry information.
 *
 */
__attribute__((always_inline)) INLINE static float firehose_recoupling_criterion(
	  struct part *pi, const float Mach, const float r_stream, 
    const struct chemistry_global_data* cd) {

  const float u_max = max(pi->u, pi->chemistry_data.u_ambient);
  const float u_diff = fabs(pi->u - pi->chemistry_data.u_ambient) / u_max;
  if (Mach < cd->firehose_recoupling_mach && 
        u_diff < cd->firehose_recoupling_u_factor) return -1.f;

  const float exchanged_mass_frac = pi->chemistry_data.exchanged_mass / pi->mass;
  if (exchanged_mass_frac > cd->firehose_recoupling_fmix) return -1.f;  
  if (r_stream == 0.f) return -1.f;

  return pi->chemistry_data.radius_stream;
}


/**
 * @brief Computes the particle interaction via the firehose stream model
 *
 * This is called from runner_iact_diffusion, which is called during the force loop
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Wind particle.
 * @param pj Ambient particle.
 * @param time_base The time base used to convert integer to float time.
 * @param ti_current Current integer time used value for seeding random number generator.
 * @param phys_const Physical constants
 * @param cd #chemistry_global_data containing chemistry information.
 *
 */
__attribute__((always_inline)) INLINE static void firehose_evolve_particle_sym(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *pi, struct part *pj,
    const float time_base, const integertime_t ti_current,
    const struct phys_const* phys_const, const struct chemistry_global_data* cd,
    const struct cosmology *cosmo) {

  /* never do diffusion between a cooling shutoff particle */
  if (pi->feedback_data.cooling_shutoff_delay_time > 0.f ||
        pj->feedback_data.cooling_shutoff_delay_time > 0.f) return;

  if (r2 <= 0.f) return;

  /* Compute the amount of mass mixed between stream particle and ambient gas */
  float v2 = 0.f;

  const float dm = firehose_compute_mass_exchange(r2, dx, hi, hj, pi, pj, 
                                                  time_base, ti_current, 
                                                  phys_const, cd, &v2,
                                                  cosmo);
  const float delta_m = fabs(dm);
  if (delta_m <= 0.f) return; 

  struct chemistry_part_data* chi = &pi->chemistry_data;
  struct chemistry_part_data* chj = &pj->chemistry_data;

  /* Track amount of gas mixed in stream particle */
  if (pi->feedback_data.decoupling_delay_time > 0.f) {
    chi->exchanged_mass += delta_m;
  }

  if (pj->feedback_data.decoupling_delay_time > 0.f) {
    chj->exchanged_mass += delta_m;
  }

  /* Mach number */ 
  const float Mach = 
      sqrtf(v2 / (chi->u_ambient * hydro_gamma * hydro_gamma_minus_one));  

  /* set weights for averaging i and j */
  const float pii_weight = (pi->mass - delta_m) / pi->mass;
  const float pij_weight = delta_m / pi->mass;
  const float pji_weight = delta_m / pj->mass;
  const float pjj_weight = (pj->mass - delta_m) / pj->mass;

  /* Mixing is negligibly small, avoid underflows */
  if (pij_weight < 1.e-10f || pji_weight < 1.e-10f) return;

  /* Mixing is erroneously large */
  if (pij_weight > 0.25f || pji_weight > 0.25f) return;

  /* 1) Update chemistry */
  chi->metal_mass_fraction_total = 0.f;
  chj->metal_mass_fraction_total = 0.f;

  const float wt_ii = pii_weight * pi->mass;
  const float wt_ij = pij_weight * pj->mass;
  const float wt_jj = pjj_weight * pj->mass;
  const float wt_ji = pji_weight * pi->mass;

  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    const float term_ii = wt_ii * chi->metal_mass_fraction[elem];
    const float term_ij = wt_ij * chj->metal_mass_fraction[elem];
    const float term_jj = wt_jj * chj->metal_mass_fraction[elem];
    const float term_ji = wt_ji * chi->metal_mass_fraction[elem];

    chi->metal_mass_fraction[elem] = (term_ii + term_ij) / pi->mass;
    chj->metal_mass_fraction[elem] = (term_jj + term_ji) / pj->mass;

    /* Keep track of the total "metallicity" Z */
    if (elem > chemistry_element_He) {
      chi->metal_mass_fraction_total += chi->metal_mass_fraction[elem];
      chj->metal_mass_fraction_total += chj->metal_mass_fraction[elem];
    }
  }

  if (chi->metal_mass_fraction_total < 0.f) {
    warning("FIREHOSE led to negative metallicity in particle i!\n"
            "\tpid=%lld\n\tZ_tot=%g\n\twt_ii=%g\n\twt_ij=%g\n"
            "\twt_jj=%g\n\twt_ji=%g\n\tdelta_m=%g\n",
            pi->id,
            chi->metal_mass_fraction_total,
            wt_ii, wt_ij, wt_jj, wt_ji,
            delta_m);
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      warning("\telem[%d]=%g\n", elem, chi->metal_mass_fraction[elem]);
    }
    error("Did not check particle j.");
    return;
  }

  if (chj->metal_mass_fraction_total < 0.f) {
    warning("FIREHOSE led to negative metallicity in particle j!\n"
            "\tpid=%lld\n\tZ_tot=%g\n\twt_ii=%g\n\twt_ij=%g\n"
            "\twt_jj=%g\n\twt_ji=%g\n\tdelta_m=%g\n",
            pj->id,
            chj->metal_mass_fraction_total,
            wt_ii, wt_ij, wt_jj, wt_ji,
            delta_m);
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      warning("\telem[%d]=%g\n", elem, chj->metal_mass_fraction[elem]);
    }
    error("Already checked particle i.");
    return;
  }

  const float total_dust_mass_ij =
          pi->cooling_data.dust_mass + pj->cooling_data.dust_mass;
  if (total_dust_mass_ij > 0.f) {
    /* Spread dust mass between particles */
    const float pi_dust_mass = pi->cooling_data.dust_mass;
    const float pj_dust_mass = pj->cooling_data.dust_mass;

    const float dust_wt_ii = pii_weight * pi_dust_mass;
    const float dust_wt_ij = pij_weight * pj_dust_mass;
    const float dust_wt_ji = pji_weight * pi_dust_mass;
    const float dust_wt_jj = pjj_weight * pj_dust_mass;

    pi->cooling_data.dust_mass = dust_wt_ii + dust_wt_ij;
    pj->cooling_data.dust_mass = dust_wt_ji + dust_wt_jj;

    /* Spread individual dust elements */
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      const float pi_dust_frac = pi->cooling_data.dust_mass_fraction[elem];
      const float pj_dust_frac = pj->cooling_data.dust_mass_fraction[elem];

      /* These will be reset if there is dust mass */
      pi->cooling_data.dust_mass_fraction[elem] = 0.f;
      pj->cooling_data.dust_mass_fraction[elem] = 0.f;

      if (pi->cooling_data.dust_mass > 0.f) {
        const float dust_term_ii = dust_wt_ii * pi_dust_frac;
        const float dust_term_ij = dust_wt_ij * pj_dust_frac;

        pi->cooling_data.dust_mass_fraction[elem] = dust_term_ii + dust_term_ij;
        pi->cooling_data.dust_mass_fraction[elem] /= pi->cooling_data.dust_mass;
      }

      if (pj->cooling_data.dust_mass > 0.f) {
        const float dust_term_jj = dust_wt_jj * pj_dust_frac;
        const float dust_term_ji = dust_wt_ji * pi_dust_frac;

        pj->cooling_data.dust_mass_fraction[elem] = dust_term_jj + dust_term_ji;
        pj->cooling_data.dust_mass_fraction[elem] /= pj->cooling_data.dust_mass;
      }
    }
  }

  /* 2) Update particles' internal energy per unit mass */
  const float pi_u = pi->u;
  const float pj_u = pj->u;

  pi->u = (wt_ii * pi_u + wt_ij * pj_u) / pi->mass;
  pj->u = (wt_ji * pi_u + wt_jj * pj_u) / pj->mass;

  /* 3) Update particles' velocities, conserving momentum */
  float pi_vfull;
  float new_v2 = 0.f;
  for (int i = 0; i < 3; i++) {
    pi_vfull = pi->v_full[i];
    pi->v_full[i] = (wt_ii * pi_vfull + wt_ij * pj->v_full[i]) / pi->mass;
    pj->v_full[i] = (wt_ji * pi_vfull + wt_jj * pj->v_full[i]) / pj->mass;
    new_v2 += (pi->v_full[i] - pj->v_full[i]) * (pi->v_full[i] - pj->v_full[i]);
  }

   /* 4) Split excess energy between stream and ambient particle */
  float delta_KE = 0.5f * delta_m * (v2 - new_v2);
  if (delta_KE > min(pi->mass * pi->u, pj->mass * pj->u)) {
    delta_KE = min(pi->mass * pi->u, pj->mass * pj->u);
  }

  /* Check extreme energies for pi->u and update */
  float new_pi_u = pi->u + 0.5f * delta_KE / pi->mass;
  float energy_fraction = pi->u > 0.f ? new_pi_u / pi->u : 0.f;

  if (energy_fraction > FIREHOSE_HEATLIM) {
    new_pi_u = FIREHOSE_HEATLIM * pi->u;
  }

  if (energy_fraction < FIREHOSE_COOLLIM || new_pi_u < 0.f) {
    new_pi_u = FIREHOSE_COOLLIM * pi->u;
  }

  pi->u = new_pi_u;

  /* Check extreme energies for pj->u and update */
  float new_pj_u = pj->u + 0.5f * delta_KE / pj->mass;
  energy_fraction = pj->u > 0.f ? new_pj_u / pj->u : 0.f;

  if (energy_fraction > FIREHOSE_HEATLIM) {
    new_pj_u = FIREHOSE_HEATLIM * pj->u;
  }

  if (energy_fraction < FIREHOSE_COOLLIM || new_pj_u < 0.f) {
    new_pj_u = FIREHOSE_COOLLIM * pj->u;
  }

  pj->u = new_pj_u;

  /* Update stream radius */
  float stream_growth_factor;
  if (pi->feedback_data.decoupling_delay_time > 0.f) {
    stream_growth_factor = (pi->mass + dm) / pi->mass;

    if (stream_growth_factor > 0.f) {
      chi->radius_stream *= sqrtf(stream_growth_factor);
    }
    else {
      chi->radius_stream = 0.f;
    }
  }
  else if (pj->feedback_data.decoupling_delay_time > 0.f) {
    stream_growth_factor = (pj->mass + dm) / pj->mass;

    if (stream_growth_factor > 0.f) {
      chj->radius_stream *= sqrtf(stream_growth_factor);
    }
    else {
      chj->radius_stream = 0.f;
    }
  }
  else {
    error("In firehose model, both i and j have negative delay times %g %g",
          pi->feedback_data.decoupling_delay_time, 
          pj->feedback_data.decoupling_delay_time);
  }

#ifdef FIREHOSE_DEBUG_CHECKS
  if (pi->feedback_data.decoupling_delay_time > 0.f) {
    message("FIREHOSE: %lld %lld dv=%g tdi=%g ui=%g uj=%g ua=%g dm=%g dE=%g "
            "Ri=%g", 
            pi->id, 
            pj->id, 
            sqrtf(v2), 
            pi->feedback_data.decoupling_delay_time, 
            pi->u, 
            pj->u, 
            pi->chemistry_data.u_ambient, 
            delta_m/pi->mass, 
            delta_KE/pi->u, 
            pi->chemistry_data.radius_stream);
  }
#endif

  /* Check if particle should recouple. 
     Negative value signifies particle should be recoupled 
     (which happens in feedback.h) */
  if (pi->feedback_data.decoupling_delay_time > 0.f) {
    chi->radius_stream = 
        firehose_recoupling_criterion(pi, Mach, 
                                      chi->radius_stream, cd);
  }  
  else if (pj->feedback_data.decoupling_delay_time > 0.f) {
    chj->radius_stream = 
        firehose_recoupling_criterion(pj, Mach, 
                                      chj->radius_stream, cd); 
  } 

  return;
}

/**
 * @brief Computes non-symmetric (single) particle interaction via the firehose stream model
 * Here particle i is the stream and j the ambient.
 *
 * This is called from runner_iact_nonsym_diffusion, which is called during the <FORCE> loop
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Stream particle.
 * @param pj Gas particle (not updated)..
 * @param time_base The time base used to convert integer to float time.
 * @param ti_current Current integer time used value for seeding random number generator   
 * @param phys_const Physical constants
 * @param cd #chemistry_global_data containing chemistry information.
 *
 */
__attribute__((always_inline)) INLINE static void firehose_evolve_stream_particle_nonsym(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *pi, const struct part *pj,
    const float time_base, const integertime_t ti_current,
    const struct phys_const* phys_const, const struct chemistry_global_data* cd,
    const struct cosmology *cosmo) {

  /* never do diffusion between a cooling shutoff particle */
  if (pi->feedback_data.cooling_shutoff_delay_time > 0.f ||
        pj->feedback_data.cooling_shutoff_delay_time > 0.f) return;

  if (r2 <= 0.f) return;

  /* Compute the interaction terms between stream particle and ambient gas */
  float v2 = 0.f;
  float dm = firehose_compute_mass_exchange(r2, dx, hi, hj, pi, pj, time_base, 
                                            ti_current, phys_const, cd, &v2,
                                            cosmo);
  float delta_m = fabs(dm);
  if (delta_m <= 0.) return;

  struct chemistry_part_data* chi = &pi->chemistry_data;
  const struct chemistry_part_data* chj = &pj->chemistry_data;

  /* Track the amount of gas mixed */
  chi->exchanged_mass += delta_m;

  /* Constants */ 
  const float Mach = 
      sqrtf(v2 / (chi->u_ambient * hydro_gamma * hydro_gamma_minus_one));  

  /* set weights for averaging i and j */
  float pii_weight = (pi->mass - delta_m) / pi->mass;
  float pij_weight = delta_m / pi->mass;
  float pji_weight = delta_m / pj->mass;
  float pjj_weight = (pj->mass - delta_m) / pj->mass;

  /* Mixing is negligibly small */
  if (pij_weight < 1.e-10f || pji_weight < 1.e-10f) return;

  /* Mixing is erroneously large */
  if (pij_weight > 0.25f || pji_weight > 0.25f) return;

  const float wt_ii = pii_weight * pi->mass;
  const float wt_ij = pij_weight * pj->mass;
  const float wt_jj = pjj_weight * pj->mass;
  const float wt_ji = pji_weight * pi->mass;

  /* 1) Update chemistry */
  chi->metal_mass_fraction_total = 0.f;
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    const float term_ii = wt_ii * chi->metal_mass_fraction[elem];
    const float term_ij = wt_ij * chj->metal_mass_fraction[elem];

    chi->metal_mass_fraction[elem] = (term_ii + term_ij) / pi->mass;

    /* Track the new "metallicity" Z */
    if (elem > chemistry_element_He) {
      chi->metal_mass_fraction_total += chi->metal_mass_fraction[elem];
    }
  }

  if (chi->metal_mass_fraction_total < 0.f) {
    warning("FIREHOSE led to negative metallicity in particle i!\n"
          "\tpid=%lld\n\tZ_tot=%g\n\twt_ii=%g\n\twt_ij=%g\n"
          "\twt_jj=%g\n\twt_ji=%g\n\tdelta_m=%g\n",
          pi->id,
          chi->metal_mass_fraction_total,
          wt_ii, wt_ij, wt_jj, wt_ji,
          delta_m);
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      warning("\telem[%d]=%g\n", elem, chi->metal_mass_fraction[elem]);
    }
    error("This is the stream nonsym case.");
    return;
  }

  /* Spread dust for testing */
  const float total_dust_mass_ij = 
      pi->cooling_data.dust_mass + pj->cooling_data.dust_mass;
  if (total_dust_mass_ij > 0.f) {
    /* Spread dust mass between particles */
    const float pi_dust_mass = pi->cooling_data.dust_mass;
    const float pj_dust_mass = pj->cooling_data.dust_mass;

    const float dust_wt_ii = pii_weight * pi_dust_mass;
    const float dust_wt_ij = pij_weight * pj_dust_mass;

    pi->cooling_data.dust_mass = dust_wt_ii + dust_wt_ij;

    if (pi->cooling_data.dust_mass > 0.f) {
      /* Spread individual dust elements */
      for (int elem = 0; elem < chemistry_element_count; ++elem) {
        const float dust_term_ii = 
          dust_wt_ii * pi->cooling_data.dust_mass_fraction[elem];
        const float dust_term_ij = 
          dust_wt_ij * pj->cooling_data.dust_mass_fraction[elem];

        pi->cooling_data.dust_mass_fraction[elem] = dust_term_ii + dust_term_ij;
        pi->cooling_data.dust_mass_fraction[elem] /= pi->cooling_data.dust_mass;
      }
    }
  }

  /* 2) Update particles' internal energy per unit mass */
  const float pi_u = pi->u;
  const float pj_u = pj->u;

  pi->u = (wt_ii * pi_u + wt_ij * pj_u) / pi->mass;

  /* 3) Update particles' velocities, conserving momentum */
  float new_v2 = 0.f;
  for (int i = 0; i < 3; i++) {
    const float pi_vfull = pi->v_full[i];
    pi->v_full[i] = (wt_ii * pi_vfull + wt_ij * pj->v_full[i]) / pi->mass;
    const float pj_vfull = (wt_ji * pi_vfull + wt_jj * pj->v_full[i]) / pj->mass;
    new_v2 += (pi->v_full[i] - pj_vfull) * (pi->v_full[i] - pj_vfull);
  }

   /* 4) Deposit excess energy into stream */
  float delta_KE = 0.5f * delta_m * (v2 - new_v2);
  if (delta_KE > min(pi->mass * pi->u, pj->mass * pj->u)) {
    delta_KE = min(pi->mass * pi->u, pj->mass * pj->u);
  }

  float new_pi_u = pi->u + 0.5f * delta_KE / pi->mass;
  const float energy_fraction = pi->u > 0.f ? new_pi_u / pi->u : 0.f;

  if (energy_fraction > FIREHOSE_HEATLIM) {
    new_pi_u = FIREHOSE_HEATLIM * pi->u;
  }

  if (energy_fraction < FIREHOSE_COOLLIM || new_pi_u < 0.f) {
    new_pi_u = FIREHOSE_COOLLIM * pi->u;
  }

  pi->u = new_pi_u;

  /* Update stream radius */
  float stream_growth_factor;
  if (pi->feedback_data.decoupling_delay_time > 0.f) {
    stream_growth_factor = (pi->mass + dm) / pi->mass;
    if (stream_growth_factor > 0.f) {
      chi->radius_stream *= sqrtf(stream_growth_factor);
    }
    else {
      chi->radius_stream = 0.f;
    }
  }

  /* Check if particle should recouple. Negative value signifies particle 
     should be recoupled (which happens in feedback.h)*/
  chi->radius_stream = 
      firehose_recoupling_criterion(pi, Mach, 
                                    chi->radius_stream, cd);

  return;
}

/**
 * @brief Computes non-symmetric (single) particle interaction via the 
 * firehose stream model.
 * 
 * Here particle j is the ambient and i the stream.
 *
 * This is called from runner_iact_nonsym_diffusion, 
 * which is called during the <FORCE> loop
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Stream particle (not updated).
 * @param pj Ambient particle.
 * @param time_base The time base used to convert integer to float time.
 * @param ti_current Current integer time used value for seeding random 
 *                   number generator   
 * @param phys_const Physical constants
 * @param us Unit system
 *
 */
__attribute__((always_inline)) INLINE static 
void firehose_evolve_ambient_particle_nonsym(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *pj, const struct part *pi,
    const float time_base, const integertime_t ti_current,
    const struct phys_const* phys_const, 
    const struct chemistry_global_data* cd,
    const struct cosmology *cosmo) {

  /* never do diffusion between a cooling shutoff particle */
  if (pi->feedback_data.cooling_shutoff_delay_time > 0.f ||
        pj->feedback_data.cooling_shutoff_delay_time > 0.f) return;

  if (r2 <= 0.f) return;

  /* Compute the interaction terms between stream particle and ambient gas */
  float v2 = 0.f;
  float dm = firehose_compute_mass_exchange(r2, dx, hi, hj, pi, pj, 
                                            time_base, ti_current, 
                                            phys_const, cd, &v2,
                                            cosmo);
  float delta_m = fabs(dm);
  if (delta_m <= 0.) return;

  const struct chemistry_part_data* chi = &pi->chemistry_data;
  struct chemistry_part_data* chj = &pj->chemistry_data;

  /* set weights for averaging i and j */
  float pii_weight = (pi->mass - delta_m) / pi->mass;
  float pij_weight = delta_m / pi->mass;
  float pji_weight = delta_m / pj->mass;
  float pjj_weight = (pj->mass - delta_m) / pj->mass;

  /* Mixing is negligbly small */
  if (pij_weight < 1.e-10f || pji_weight < 1.e-10f) return;

  /* Mixing is erroneously large */
  if (pij_weight > 0.25f || pji_weight > 0.25f) return;

  const float wt_ii = pii_weight * pi->mass;
  const float wt_ij = pij_weight * pj->mass;
  const float wt_jj = pjj_weight * pj->mass;
  const float wt_ji = pji_weight * pi->mass;

  /* 1) Update chemistry */
  chj->metal_mass_fraction_total = 0.f;
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    const float term_ji = wt_ji * chi->metal_mass_fraction[elem];
    const float term_jj = wt_jj * chj->metal_mass_fraction[elem];

    chj->metal_mass_fraction[elem] = (term_ji + term_jj) / pj->mass;

    /* Check the total "metallicity" Z */
    if (elem > chemistry_element_He) {
      chj->metal_mass_fraction_total += chj->metal_mass_fraction[elem];
    }
  }

  if (chj->metal_mass_fraction_total < 0.f) {
    warning("FIREHOSE led to negative metallicity in particle j!\n"
            "\tpid=%lld\n\tZ_tot=%g\n\twt_ii=%g\n\twt_ij=%g\n"
            "\twt_jj=%g\n\twt_ji=%g\n\tdelta_m=%g\n",
            pj->id,
            chj->metal_mass_fraction_total,
            wt_ii, wt_ij, wt_jj, wt_ji,
            delta_m);
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      warning("\telem[%d]=%g\n", elem, chj->metal_mass_fraction[elem]);
    }
    error("This is the ambient nonsym case.");
    return;
  }
  
  /* Spread dust for testing */
  const float total_dust_mass_ij = 
      pi->cooling_data.dust_mass + pj->cooling_data.dust_mass;
  if (total_dust_mass_ij > 0.f) {
    /* Spread dust mass between particles */
    const float pi_dust_mass = pi->cooling_data.dust_mass;
    const float pj_dust_mass = pj->cooling_data.dust_mass;

    const float dust_wt_ji = pji_weight * pi_dust_mass;
    const float dust_wt_jj = pjj_weight * pj_dust_mass;

    pj->cooling_data.dust_mass = dust_wt_ji + dust_wt_jj;

    if (pj->cooling_data.dust_mass > 0.f) {
      /* Spread individual dust elements */
      for (int elem = 0; elem < chemistry_element_count; ++elem) {
        const float dust_term_ji = 
            dust_wt_ji * pi->cooling_data.dust_mass_fraction[elem];
        const float dust_term_jj =
            dust_wt_jj * pj->cooling_data.dust_mass_fraction[elem];

        pj->cooling_data.dust_mass_fraction[elem] = dust_term_ji + dust_term_jj;
        pj->cooling_data.dust_mass_fraction[elem] /= pj->cooling_data.dust_mass;
      }
    }
  }

  const float pi_u = pi->u;
  const float pj_u = pj->u;

  /* 2) Update particles' internal energy per unit mass */
  pj->u = (wt_ii * pi_u + wt_ij * pj_u) / pi->mass;

  /* 3) Update particles' velocities, conserving momentum */
  float new_v2 = 0.f;
  for (int i = 0; i < 3; i++){
    pj->v_full[i] = (wt_ji * pi->v_full[i] + wt_jj * pj->v_full[i]) / pj->mass;
    new_v2 += (pi->v_full[i] - pj->v_full[i]) * (pi->v_full[i] - pj->v_full[i]);
  }

   /* 4) Deposit excess energy into ambient */
  float delta_KE = 0.5f * delta_m * (v2 - new_v2);
  if (delta_KE > min(pi->mass * pi->u, pj->mass * pj->u)) {
    delta_KE = min(pi->mass * pi->u, pj->mass * pj->u);
  }

  float new_pj_u = pj->u + 0.5f * delta_KE / pj->mass;
  const float energy_fraction = pj->u > 0.f ? new_pj_u / pj->u : 0.f;

  if (energy_fraction > FIREHOSE_HEATLIM) {
    new_pj_u = FIREHOSE_HEATLIM * pj->u;
  }

  if (energy_fraction < FIREHOSE_COOLLIM || new_pj_u < 0.f) {
    new_pj_u = FIREHOSE_COOLLIM * pj->u;
  }

  pj->u = new_pj_u;

  return;
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
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology, 
    const struct phys_const* phys_const, const struct chemistry_global_data *cd) {

  /* never do diffusion between a cooling shutoff particle */
  if (pi->feedback_data.cooling_shutoff_delay_time > 0.f ||
        pj->feedback_data.cooling_shutoff_delay_time > 0.f) return;

  if (pi->feedback_data.decoupling_delay_time > 0.f || 
        pj->feedback_data.decoupling_delay_time > 0.f) {
    /* If in wind mode, do firehose wind diffusion */
    firehose_evolve_particle_sym(r2, dx, hi, hj, pi, pj, time_base, 
                                 t_current, phys_const, cd, cosmo);
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
    const struct phys_const* phys_const, const struct chemistry_global_data *cd) {

  /* never do diffusion between a cooling shutoff particle */
  if (pi->feedback_data.cooling_shutoff_delay_time > 0.f ||
        pj->feedback_data.cooling_shutoff_delay_time > 0.f) return;

  /* In nonsym case, two cases: depending on whether i is stream or ambient */
  if (pi->feedback_data.decoupling_delay_time > 0.f && 
        pj->feedback_data.decoupling_delay_time <= 0.f) {
    /* Here, i is the stream particle, j is the ambient */
    firehose_evolve_stream_particle_nonsym(r2, dx, hi, hj, pi, pj, 
                                           time_base, t_current, phys_const, cd,
                                           cosmo);
    return;
  }

  if (pi->feedback_data.decoupling_delay_time <= 0.f && 
        pj->feedback_data.decoupling_delay_time > 0.f) {
    /* Here, i is ambient. note that inside the called routine, 
       i and j are switched, so j is ambient there */
    firehose_evolve_ambient_particle_nonsym(r2, dx, hi, hj, pi, pj, 
                                            time_base, t_current, phys_const, cd,
                                            cosmo);
    return;
  }

  /* If we are here and one of them is still positive then skip diffusion */
  if (pi->feedback_data.decoupling_delay_time > 0.f ||
        pj->feedback_data.decoupling_delay_time > 0.f) return;

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
