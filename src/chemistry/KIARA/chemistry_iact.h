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
 * @param p The particle to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void firehose_compute_ambient_quantities(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj) {

  struct chemistry_part_data* cpd = &pi->chemistry_data;

  if (cpd->weight_ambient < 0.f) return;  // this signified that firehose model is not on

  /* Get r */
  if (r2 > 0.f) {
    //const float X_H = chemistry_get_metal_mass_fraction_for_cooling(p)[chemistry_element_H];
    //const float yhelium = (1. - X_H) / (4. * X_H);
    //const float mu = (1. + yhelium) / (1. + ne + 4. * yhelium);
    //const float Tj = pj->u / (mu * cooling->temp_to_u_factor);

    /* Compute the kernel function for pi */
    const float r = sqrtf(r2);
    const float ui = r / hi;
    float wi, wi_dx;
    kernel_deval(ui, &wi, &wi_dx);

    /* Sum ambient neighbour quantities */
    cpd->u_ambient += pj->u * wi;
    cpd->rho_ambient += pj->rho * wi;
    cpd->weight_ambient += wi;
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

  if (pi->feedback_data.decoupling_delay_time > 0.f) {
    /* If in wind mode, compute ambient quantities for firehose wind diffusion */
    firehose_compute_ambient_quantities(r2, dx, hi, hj, pi, pj);
  }

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

    /* Compute the shear tensor */
  for (int i = 0; i < 3; i++) {
    const float dxi_mj_wi_dr = dx[i] * mj_wi_dr;
    const float dxi_mi_wj_dr = dx[i] * mi_wj_dr;

    chi->shear_tensor[i][0] += (pj->v[0] - pi->v[0]) * dxi_mj_wi_dr;
    chi->shear_tensor[i][1] += (pj->v[1] - pi->v[1]) * dxi_mj_wi_dr;
    chi->shear_tensor[i][2] += (pj->v[2] - pi->v[2]) * dxi_mj_wi_dr;

    chj->shear_tensor[i][0] -= (pj->v[0] - pi->v[0]) * dxi_mi_wj_dr;
    chj->shear_tensor[i][1] -= (pj->v[1] - pi->v[1]) * dxi_mi_wj_dr;
    chj->shear_tensor[i][2] -= (pj->v[2] - pi->v[2]) * dxi_mi_wj_dr;
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

  /* Compute the shear tensor */
  for (int i = 0; i < 3; i++) {
    const float dxi_mj_wi_dr = dx[i] * mj_wi_dr;
    chi->shear_tensor[i][0] += (pj->v[0] - pi->v[0]) * dxi_mj_wi_dr;
    chi->shear_tensor[i][1] += (pj->v[1] - pi->v[1]) * dxi_mj_wi_dr;
    chi->shear_tensor[i][2] += (pj->v[2] - pi->v[2]) * dxi_mj_wi_dr;
  }
}


/**
 * @brief Computes the mass exchange between the firehose stream and the ambient medium
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Wind particle (not updated).
 * @param pj Gas particle.
 * @param xpi Extra particle data (wind)
 * @param xpj Extra particle data (gas)
 * @param ti_current Current integer time used value for seeding random number generator   
 * 
 * @param phys_const Physical constants
 * @param us Unit system
 *
 */
__attribute__((always_inline)) INLINE static float firehose_compute_mass_exchange(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *pi, struct part *pj,
    const float time_base, const integertime_t ti_current,
    const struct phys_const* phys_const, float *wi, float *chi, float *v2) {

  /* Are we using firehose model? If not, return */
  if (pi->chemistry_data.weight_ambient <= 0.f) return 0.;
  
  /* Ignore COUPLED particles */
  if (pi->feedback_data.decoupling_delay_time <= 0.f) return 0.;
  
  /* No wind-wind interaction */
  if (pj->feedback_data.decoupling_delay_time > 0.f) return 0.;

  //message("FIREHOSE: %g %g %g\n", pi->chemistry_data.weight_ambient, pi->feedback_data.decoupling_delay_time, pj->feedback_data.decoupling_delay_time);

  /* Gas particle density */
  float rho_j = hydro_get_comoving_density(pj);

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  kernel_eval(r/hi, wi);
  if (*wi <= 0.) return 0.;  // Too far out to interact
  /* This is the fraction of each quantity which particle j will exchange, from among all ambient gas */
  *wi = *wi / pi->chemistry_data.weight_ambient;

  /* Get the timestep */ 
  const float dt = get_timestep(pi->time_bin, time_base);

  /* Compute temperatures for stream (i) and ambient (j), with a floor at 10^4 K */
  const float T_floor = 1.e4;
  const float temp_to_u_factor = phys_const->const_boltzmann_k / (hydro_gamma_minus_one * phys_const->const_proton_mass); //assumes units of T are Kelvin
  const float X_Hi = pi->chemistry_data.metal_mass_fraction[chemistry_element_H];
  const float X_Hj = pj->chemistry_data.metal_mass_fraction[chemistry_element_H];
  const float X_Hei = pi->chemistry_data.metal_mass_fraction[chemistry_element_He];
  const float X_Hej = pj->chemistry_data.metal_mass_fraction[chemistry_element_He];
  const float yhelium_i = (1. - X_Hi) / (4. * X_Hi);
  const float yhelium_j = (1. - X_Hj) / (4. * X_Hj);
  const float e_fraci = X_Hi + 0.5 * X_Hei; // assume fully ionised
  const float e_fracj = X_Hj + 0.5 * X_Hej; // assume fully ionised
  const float mu_i = (1. + yhelium_i) / (1. + e_fraci + 4. * yhelium_i);
  const float mu_j = (1. + yhelium_j) / (1. + e_fracj + 4. * yhelium_j);

  const float T_amb =  max(pi->chemistry_data.u_ambient / (mu_i * temp_to_u_factor), T_floor);
  const float Tstream =  max(pi->u / (mu_i * temp_to_u_factor), T_floor);
  *chi = T_amb * mu_i / (Tstream * mu_j); //assuming collisional equilibrium

  /* Compute the velocity of the stream relative to ambient gas */
  *v2 = 0.f;
  for (int i=0; i<3; i++){
    *v2 += (pi->v_full[i] - pj->v_full[i]);
  }
  //const float Mach = sqrtf(*v2/(pi->chemistry_data.u_ambient * hydro_gamma * hydro_gamma_minus_one));  

  /* Define cooling and destruction timescales for streams*/
  const float Lambda_mix = pi->cooling_data.mixing_layer_cool_rate;
  float n_mix = sqrtf(pi->chemistry_data.rho_ambient * pi->rho) / phys_const->const_proton_mass;

  float tcoolmix = phys_const->const_boltzmann_k * sqrt(Tstream*T_amb) / (hydro_gamma_minus_one * n_mix * Lambda_mix);
  float tshear = pi->chemistry_data.radius_stream / sqrt(*v2);

  //float virtual_mass = *chi * rho_j * pow(pi->chemistry_data.radius_stream,2) * M_PI * sqrtf(*v2) * dt;
  float mdot = 0.f;

  /* If the stream is losing mass, updated destruction time and compute mdot */
  if (tcoolmix / tshear >= 1.f){

    /* If cloud just began to be destroyed, update initial mass */
    //if (pi->chemistry_data.destruction_time == 0.f || Mach < 1){
    //  pi->chemistry_data.initial_mass = virtual_mass;
    //}

    pi->chemistry_data.destruction_time += dt;
    mdot = -1. / tshear * pi->mass * exp(-pi->chemistry_data.destruction_time / tshear);
    message("FIREHOSE shear: %g %g %g %g\n", pi->chemistry_data.initial_mass, pi->chemistry_data.radius_stream, pi->chemistry_data.destruction_time, tshear);
  }
  else {
    /* If stream is growing, cancel destruction time and update mdot*/

    pi->chemistry_data.destruction_time = 0.f;
    /* Compute sound crossing time */
    float tsc = 2 * pi->chemistry_data.radius_stream / sqrt(pi->u / *chi / rho_j * hydro_gamma / hydro_gamma_minus_one);
    mdot = 4. / *chi * pi->mass / tsc * pow(tcoolmix / tsc, -0.25);
  
    message("FIREHOSE cool: %lld %g %g %g %g\n", pi->id, pi->chemistry_data.initial_mass, pi->chemistry_data.radius_stream, tsc, tcoolmix);
  } 

  pi->chemistry_data.initial_mass += fabs(mdot);  // track amount of gas mixed

  return mdot;
}


/**
 * @brief Check recoupling criterion for firehose stream particle 
 *
 * @param pi Wind particle (not updated).
 * @param Mach Stream Mach number vs ambient
 * @param r_stream Current radius of stream 
 * 
 * @param phys_const Physical constants
 * @param us Unit system
 *
 */
__attribute__((always_inline)) INLINE static float firehose_recoupling_criterion(
	struct part *pi, const float Mach, const float r_stream) {

  // negative value indicates it should recouple
  if (Mach < 1.f) return -1.;  
  if (pi->chemistry_data.initial_mass / pi->mass > 0.9f) return -1.;  
  if (r_stream == 0.f) return -1.;
  return pi->chemistry_data.weight_ambient;
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
 * @param pi Wind particle (not updated).
 * @param pj Gas particle.
 * @param xpi Extra particle data (wind)
 * @param xpj Extra particle data (gas)
 * @param ti_current Current integer time used value for seeding random number generator   
 * 
 * @param phys_const Physical constants
 * @param us Unit system
 *
 */
__attribute__((always_inline)) INLINE static void firehose_evolve_particle_sym(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *pi, struct part *pj,
    const float time_base, const integertime_t ti_current,
    const struct phys_const* phys_const) {

  /* Compute the interaction terms between stream particle and ambient gas */
  float wi, chi, v2;
  float mdot = firehose_compute_mass_exchange(r2, dx, hi, hj, pi, pj, time_base, ti_current, phys_const, &wi, &chi, &v2);
  float delta_m = fabs(mdot);
  if (delta_m <= 0.) return; // no interaction, nothing to do

  /* Constants */ 
  const float dt = get_timestep(pi->time_bin, time_base);
  const float Mach = sqrtf(v2 / (pi->chemistry_data.u_ambient * hydro_gamma * hydro_gamma_minus_one));  

  /* set weights for averaging i and j */
  if (delta_m > min(0.5 * pi->mass, 0.5 * pj->mass)) delta_m = min(0.5 * pi->mass, 0.5 * pj->mass);  // limit single-step mixing
  float pii_weight = wi * (pi->mass - delta_m) / pi->mass;
  float pij_weight = wi * delta_m / pi->mass;
  float pji_weight = wi * delta_m / pj->mass;
  float pjj_weight = wi * (pj->mass - delta_m) / pj->mass;

  /* 1) Update chemistry */
  pi->chemistry_data.metal_mass_fraction_total = 0.f;
  pj->chemistry_data.metal_mass_fraction_total = 0.f;
  float pi_metal_frac;
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    pi_metal_frac = pi->chemistry_data.metal_mass_fraction[elem];
    pi->chemistry_data.metal_mass_fraction[elem] = (pii_weight * pi->mass * pi->chemistry_data.metal_mass_fraction[elem] + pij_weight * pj->mass * pj->chemistry_data.metal_mass_fraction[elem]) / pi->mass;
    pj->chemistry_data.metal_mass_fraction[elem] = (pjj_weight * pj->mass * pj->chemistry_data.metal_mass_fraction[elem] + pji_weight * pi->mass * pi_metal_frac) / pj->mass;
    if (elem > chemistry_element_He) {
      pi->chemistry_data.metal_mass_fraction_total += pi->chemistry_data.metal_mass_fraction[elem];
      pj->chemistry_data.metal_mass_fraction_total += pj->chemistry_data.metal_mass_fraction[elem];
    }
  }

  int spread_dust=1; // flag (for testing) whether to spread dust in addition to metals
  message("FIREHOSE before: %g %g %g %g %g %g %g %g\n", pi->cooling_data.dust_mass, pj->cooling_data.dust_mass, pi->cooling_data.dust_mass_fraction[2], pj->cooling_data.dust_mass_fraction[2], delta_m, mdot, chi, v2);
  if (spread_dust && (pi->cooling_data.dust_mass+pj->cooling_data.dust_mass) > 0.f) {
    /* Spread dust mass between particles */
    float pi_dust_mass = pi->cooling_data.dust_mass;
    float pj_dust_mass = pj->cooling_data.dust_mass;
    pi->cooling_data.dust_mass = pii_weight * pi_dust_mass + pij_weight * pj_dust_mass;
    pj->cooling_data.dust_mass = pji_weight * pi_dust_mass + pjj_weight * pj_dust_mass;
    assert(pi->cooling_data.dust_mass*pj->cooling_data.dust_mass > 0.);
    /* Spread individual dust elements */
    float pi_dust_frac;
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      pi_dust_frac = pi->cooling_data.dust_mass_fraction[elem];
      pi->cooling_data.dust_mass_fraction[elem] = (pii_weight * pi_dust_mass * pi->cooling_data.dust_mass_fraction[elem] + pij_weight * pj_dust_mass * pj->cooling_data.dust_mass_fraction[elem]) / pi->cooling_data.dust_mass;
      pj->cooling_data.dust_mass_fraction[elem] = (pjj_weight * pj_dust_mass * pj->cooling_data.dust_mass_fraction[elem] + pji_weight * pi_dust_mass * pi_dust_frac) / pj->cooling_data.dust_mass;
    }
  }

  /* 2) Update particles' internal energy per unit mass */
  float pi_u = pi->u;
  pi->u = (pii_weight * pi->mass * pi_u + pij_weight * pj->mass * pj->u) / pi->mass;
  pj->u = (pji_weight * pi->mass * pi_u + pjj_weight * pj->mass * pj->u) / pj->mass;


    /* 3) Update particles' velocities, conserving momentum */
  float pi_vfull, new_v2 = 0.f;
  for (int i=0; i<3; i++){
    pi_vfull = pi->v_full[i];
    pi->v_full[i] = (pii_weight * pi->mass * pi_vfull + pij_weight * pj->mass * pj->v_full[i]) / pi->mass;
    pj->v_full[i] = (pji_weight * pi->mass * pi_vfull + pjj_weight * pj->mass * pj->v_full[i]) / pj->mass;
    new_v2 += (pi->v_full[i]-pj->v_full[i]) * (pi->v_full[i]-pj->v_full[i]);
  }

   /* 4) Deposit excess energy onto stream */
  float delE = 0.5f * delta_m * (v2 - new_v2);
  pi->u += wi * delE / pi->mass;

  /* Update stream radius */
  float stream_growth_factor = (pi->mass + wi * mdot * dt) / pi->mass;
  if (stream_growth_factor > 0.f) pi->chemistry_data.radius_stream *= sqrtf(stream_growth_factor);
  else pi->chemistry_data.radius_stream = 0.;

  /* Check if particle should recouple */
  pi->chemistry_data.weight_ambient = firehose_recoupling_criterion(pi, Mach, pi->chemistry_data.radius_stream);  // Negative value signifies particle should be recoupled (which happens in feedback.h)

  return;
}

/**
 * @brief Computes non-symmetric (single) particle interaction via the firehose stream model
 *
 * This is called from runner_iact_nonsym_diffusion, which is called during the force loop
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Wind particle (not updated).
 * @param pj Gas particle.
 * @param xpi Extra particle data (wind)
 * @param xpj Extra particle data (gas)
 * @param ti_current Current integer time used value for seeding random number generator   
 * 
 * @param phys_const Physical constants
 * @param us Unit system
 *
 */
__attribute__((always_inline)) INLINE static void firehose_evolve_particle_nonsym(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *pi, struct part *pj,
    const float time_base, const integertime_t ti_current,
    const struct phys_const* phys_const) {

  /* Compute the interaction terms between stream particle and ambient gas */
  float wi, chi, v2;
  float mdot = firehose_compute_mass_exchange(r2, dx, hi, hj, pi, pj, time_base, ti_current, phys_const, &wi, &chi, &v2);
  float delta_m = fabs(mdot);
  if (delta_m <= 0.) return; // no interaction, nothing to do

  /* Constants */ 
  const float dt = get_timestep(pi->time_bin, time_base);
  const float Mach = sqrtf(v2/(pi->chemistry_data.u_ambient * hydro_gamma * hydro_gamma_minus_one));  

  /* set weights for averaging i and j */
  float pii_weight = wi * (pi->mass - delta_m) / pi->mass;
  float pij_weight = wi * delta_m / pi->mass;
  float pji_weight = wi * delta_m / pj->mass;
  float pjj_weight = wi * (pj->mass - delta_m) / pj->mass;

  /* 1) Update chemistry */
  pi->chemistry_data.metal_mass_fraction_total = 0.f;
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    pi->chemistry_data.metal_mass_fraction[elem] = (pii_weight * pi->mass * pi->chemistry_data.metal_mass_fraction[elem] + pij_weight * pj->mass * pj->chemistry_data.metal_mass_fraction[elem]) / pi->mass;
    if (elem > chemistry_element_He) {
      pi->chemistry_data.metal_mass_fraction_total += pi->chemistry_data.metal_mass_fraction[elem];
    }
  }

  int spread_dust=1;
  if (spread_dust && (pi->cooling_data.dust_mass+pj->cooling_data.dust_mass) > 0.f) {
    /* Spread dust mass between particles */
    float pi_dust_mass = pi->cooling_data.dust_mass;
    float pj_dust_mass = pj->cooling_data.dust_mass;
    pi->cooling_data.dust_mass = pii_weight * pi_dust_mass + pij_weight * pj_dust_mass;
    /* Spread individual dust elements */
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      pi->cooling_data.dust_mass_fraction[elem] = (pii_weight * pi_dust_mass * pi->cooling_data.dust_mass_fraction[elem] + pij_weight * pj_dust_mass * pj->cooling_data.dust_mass_fraction[elem]) / pi->cooling_data.dust_mass;
    }
  }
  //message("FIREHOSE after: %g %g %g %g\n", pi->cooling_data.dust_mass, pj->cooling_data.dust_mass, pi->cooling_data.dust_mass_fraction[2], pj->cooling_data.dust_mass_fraction[2]);

  /* 2) Update particles' internal energy per unit mass */
  float pi_u = pi->u;
  pi->u = (pii_weight * pi->mass * pi_u + pij_weight * pj->mass * pj->u) / pi->mass;


    /* 3) Update particles' velocities, conserving momentum */
  float pi_vfull, pj_vfull, new_v2 = 0.f;
  for (int i=0; i<3; i++){
    pi_vfull = pi->v_full[i];
    pj_vfull = pj->v_full[i];
    pi->v_full[i] = (pii_weight * pi->mass * pi_vfull + pij_weight * pj->mass * pj->v_full[i]) / pi->mass;
    pj_vfull = (pji_weight * pi->mass * pi_vfull + pjj_weight * pj->mass * pj->v_full[i]) / pj->mass;
    new_v2 += (pi->v_full[i]-pj_vfull) * (pi->v_full[i]-pj_vfull);
  }

   /* 4) Deposit excess energy onto stream */
  float delE = 0.5f * delta_m * (v2 - new_v2);
  pi->u += wi * delE / pi->mass;

  /* Update stream radius */
  float stream_growth_factor = (pi->mass + wi * mdot * dt) / pi->mass;
  if (stream_growth_factor > 0.f) pi->chemistry_data.radius_stream *= sqrtf(stream_growth_factor);
  else pi->chemistry_data.radius_stream = 0.;

  /* Check if particle should recouple */
  pi->chemistry_data.weight_ambient = firehose_recoupling_criterion(pi, Mach, pi->chemistry_data.radius_stream);  // Negative value signifies particle should be recoupled (which happens in feedback.h)

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
 * @param time_base The time base used in order to convert integer to float
 * time.
 * @param ti_current The current time (in integer)
 * @param cosmo The #cosmology.
 * @param with_cosmology Are we running with cosmology?
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_diffusion(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology, 
    const struct phys_const* phys_const) {

  if (pi->feedback_data.decoupling_delay_time > 0.f) {
    /* If in wind mode, do firehose wind diffusion */
    firehose_evolve_particle_sym(r2, dx, hi, hj, pi, pj, time_base, t_current, phys_const);
    return;
  }

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  /* No need to diffuse if both particles are not diffusing. */
  if (chj->diffusion_coefficient > 0 && chi->diffusion_coefficient > 0) {

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
    const float r_inv = 1.f / sqrtf(r2);

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
    const float dZ_ij_tot = chi->metal_mass_fraction_total - chj->metal_mass_fraction_total;
    chi->dZ_dt_total += coef_i * dZ_ij_tot;
    chj->dZ_dt_total -= coef_j * dZ_ij_tot;

    for (int elem = 0; elem < chemistry_element_count; elem++) {
      const float dZ_ij = chi->metal_mass_fraction[elem] - chj->metal_mass_fraction[elem];
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
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param time_base The time base used in order to convert integer to float
 * time.
 * @param ti_current The current time (in integer)
 * @param cosmo The #cosmology.
 * @param with_cosmology Are we running with cosmology?
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_diffusion(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct phys_const* phys_const) {

  if (pi->feedback_data.decoupling_delay_time > 0.f) {
    /* If in wind mode, do firehose wind diffusion */
    firehose_evolve_particle_nonsym(r2, dx, hi, hj, pi, pj, time_base, t_current, phys_const);
    return;
  }

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

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
    const float dZ_ij_tot = chi->metal_mass_fraction_total - chj->metal_mass_fraction_total;
    chi->dZ_dt_total += coef_i * dZ_ij_tot;

    for (int elem = 0; elem < chemistry_element_count; elem++) {
      const float dZ_ij = chi->metal_mass_fraction[elem] - chj->metal_mass_fraction[elem];
      chi->dZ_dt[elem] += coef_i * dZ_ij;
    }
  }
}


#endif /* SWIFT_KIARA_CHEMISTRY_IACT_H */
