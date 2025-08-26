/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk) &
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_GASOLINE_HYDRO_IACT_H
#define SWIFT_GASOLINE_HYDRO_IACT_H

/**
 * @file Gasoline/hydro_iact.h
 * @brief Density-Energy conservative implementation of SPH,
 *        with added Gasoline physics (Wadsley+ 2017) (interaction routines)
 */

#include "adaptive_softening_iact.h"
#include "adiabatic_index.h"
#include "hydro_parameters.h"
#include "minmax.h"
#include "signal_velocity.h"

/**
 * @brief Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {
  float wi, wj, wi_dx, wj_dx;

  const int decoupled_i = pi->decoupled;
  const int decoupled_j = pj->decoupled;

  if (decoupled_i && decoupled_j) return;

  const float r = sqrtf(r2);

  /* Compute density of pi. */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;
  
  if (!decoupled_j) {
    kernel_deval(ui, &wi, &wi_dx);
  }
  else {
    wi = 0.f;
    wi_dx = 0.f;
  }

  /* Compute density of pj. */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;

  if (!decoupled_i) {
    kernel_deval(uj, &wj, &wj_dx);
  }
  else {
    wj = 0.f;
    wj_dx = 0.f;
  }

  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  adaptive_softening_add_correction_term(pi, ui, hi_inv, mj);

  pj->rho += mi * wj;
  pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);

  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);

  adaptive_softening_add_correction_term(pj, uj, hj_inv, mi);

  /* Now we need to compute the div terms */
  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  /* Smooth pressure gradient */
  pi->smooth_pressure_gradient[0] += faci * pj->u * dx[0];
  pi->smooth_pressure_gradient[1] += faci * pj->u * dx[1];
  pi->smooth_pressure_gradient[2] += faci * pj->u * dx[2];

  pj->smooth_pressure_gradient[0] -= facj * pi->u * dx[0];
  pj->smooth_pressure_gradient[1] -= facj * pi->u * dx[1];
  pj->smooth_pressure_gradient[2] -= facj * pi->u * dx[2];

  /* Finally, the big boy; the velocity gradient tensor. Note that the
   * loops here are over the coordinates, i=0 -> x, and so on. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      const float dv_ij = pi->v[i] - pj->v[i];
      const float dx_ij = pi->x[j] - pj->x[j];

      pi->viscosity.velocity_gradient[i][j] +=
          mj * dv_ij * dx_ij * wi_dx * r_inv;
      pj->viscosity.velocity_gradient[i][j] +=
          mi * dv_ij * dx_ij * wj_dx * r_inv;
    }
  }

  /* Correction factors for kernel gradients, and norm for the velocity
   * gradient. */

  pi->weighted_wcount += mj * r2 * wi_dx * r_inv;
  pj->weighted_wcount += mi * r2 * wj_dx * r_inv;
  if (r < hi && r < hj) {
    pi->weighted_self_wcount += mj * r2 * wi_dx * r_inv;
    pj->weighted_self_wcount += mi * r2 * wj_dx * r_inv;
  }

}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, const struct part* restrict pj, const float a,
    const float H) {
  float wi, wi_dx;

  /* In the non-sym case only the neighbor matters */
  const int decoupled_j = pj->decoupled;
  if (decoupled_j) return;

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);

  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  adaptive_softening_add_correction_term(pi, ui, h_inv, mj);

  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = mj * wi_dx * r_inv;

  /* Compute pressure gradient */
  pi->smooth_pressure_gradient[0] += faci * pj->u * dx[0];
  pi->smooth_pressure_gradient[1] += faci * pj->u * dx[1];
  pi->smooth_pressure_gradient[2] += faci * pj->u * dx[2];

  /* Finally, the big boy; the velocity gradient tensor. Note that the
   * loops here are over the coordinates, i=0 -> x, and so on. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      const float dv_ij = pi->v[i] - pj->v[i];
      const float dx_ij = pi->x[j] - pj->x[j];

      pi->viscosity.velocity_gradient[i][j] +=
          mj * dv_ij * dx_ij * wi_dx * r_inv;
    }
  }

  /* Correction factors for kernel gradients, and norm for the velocity
   * gradient. */
  pi->weighted_wcount += mj * r2 * wi_dx * r_inv;
  if (r < hi && r < hj) {
    pi->weighted_self_wcount += mj * r2 * wi_dx * r_inv;
  }

}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * This method wraps around hydro_gradients_collect, which can be an empty
 * method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {

  const int decoupled_i = pi->decoupled;
  const int decoupled_j = pj->decoupled;

  if (decoupled_i && decoupled_j) return;

  /* We need to construct the maximal signal velocity between our particle
   * and all of it's neighbours */

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Cosmology terms for the signal velocity */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */

  const float dvdr_Hubble = dvdr + a2_Hubble * r2;
  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float new_v_sig =
      signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

  /* Update if we need to */
  if (!decoupled_j) {
    pi->viscosity.v_sig = max(pi->viscosity.v_sig, new_v_sig);
  }
  if (!decoupled_i) {
    pj->viscosity.v_sig = max(pj->viscosity.v_sig, new_v_sig);
  }
  if (pj->viscosity.v_sig > 3.e5) message("V_SIG: sym z=%g idj=%lld idi=%lld vsigj=%g vsigi=%g dec=%d tdec=%g dvdr=%g omij=%g muij=%g r=%g", 1./a - 1., pj->id, pi->id, pj->viscosity.v_sig, pi->viscosity.v_sig, pj->decoupled, pj->feedback_data.decoupling_delay_time, dvdr, omega_ij, mu_ij, 1./r_inv);

  /* Calculate Del^2 u for the thermal diffusion coefficient. */
  /* Need to get some kernel values F_ij = wi_dx */
  float wi, wi_dx, wj, wj_dx;

  if (!decoupled_j) {
    const float ui = r / hi;
    kernel_deval(ui, &wi, &wi_dx);
  }
  else {
    wi = 0.f;
    wi_dx = 0.f;
  }
  if (!decoupled_i) {
    const float uj = r / hj;
    kernel_deval(uj, &wj, &wj_dx);
  }
  else {
    wj = 0.f;
    wj_dx = 0.f;
  }

  /* Calculate the shock limiter component */
  const float shock_ratio_i =
      pj->viscosity.tensor_norm > 0.f
          ? pj->viscosity.shock_indicator / pj->viscosity.tensor_norm
          : 0.f;

  const float shock_ratio_j =
      pi->viscosity.tensor_norm > 0.f
          ? pi->viscosity.shock_indicator / pi->viscosity.tensor_norm
          : 0.f;

  /* Rennehan: Add in complicated shock weighting to reduce noise
   * Eq 29 Wadsley+'17 */
  const float distance_norm_i = 2.f * kernel_gamma * hi;
  const float distance_norm_j = 2.f * kernel_gamma * hj;
  float w_R_i = 0.f;
  float w_R_j = 0.f;
  if (r < distance_norm_i) {
    const float w_R_i_core = 1.f - (r / distance_norm_i);
    const float w_R_i_core2 = w_R_i_core * w_R_i_core;
    w_R_i = w_R_i_core2 * w_R_i_core2;
  }
  if (r < distance_norm_j) {
    const float w_R_j_core = 1.f - (r / distance_norm_j);
    const float w_R_j_core2 = w_R_j_core * w_R_j_core;
    w_R_j = w_R_j_core2 * w_R_j_core2;
  }

  pi->viscosity.shock_limiter += pj->mass * shock_ratio_i * w_R_i;
  pj->viscosity.shock_limiter += pi->mass * shock_ratio_j * w_R_j;

  /* Rennehan: Collect the norm of the shock limiter */
  pi->viscosity.shock_limiter_norm += pj->mass * w_R_i;
  pj->viscosity.shock_limiter_norm += pi->mass * w_R_j;

  /* Correction factors for kernel gradients */
  const float rho_inv_i = 1.f / pi->rho;
  const float rho_inv_j = 1.f / pj->rho;

  if (r < hi && r < hj) {
    pi->weighted_neighbour_wcount += pj->mass * r2 * wi_dx * rho_inv_j * r_inv;
    pj->weighted_neighbour_wcount += pi->mass * r2 * wj_dx * rho_inv_i * r_inv;
  }

  /* Gradient of the density field */
  for (int j = 0; j < 3; j++) {
    const float drho_ij = pi->rho - pj->rho;
    const float dx_ij = pi->x[j] - pj->x[j];

    pi->rho_gradient[j] +=
        pj->mass * drho_ij * dx_ij * wi_dx * r_inv;
    pj->rho_gradient[j] +=
        pi->mass * drho_ij * dx_ij * wj_dx * r_inv;
  }
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * This method wraps around hydro_gradients_nonsym_collect, which can be an
 * empty method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {
  
  /* In the non-sym case only the neighbor matters */
  const int decoupled_j = pj->decoupled;
  if (decoupled_j) return;

  /* We need to construct the maximal signal velocity between our particle
   * and all of it's neighbours */

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Cosmology terms for the signal velocity */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */

  const float dvdr_Hubble = dvdr + a2_Hubble * r2;
  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float new_v_sig =
      signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

  /* Update if we need to */
  pi->viscosity.v_sig = max(pi->viscosity.v_sig, new_v_sig);

  if (pi->viscosity.v_sig > 3.e5) message("V_SIG: nonsym z=%g idi=%lld idj=%lld vsigi=%g vsigj=%g dec=%d tdec=%g dvdr=%g omij=%g muij=%g r=%g", 1./a - 1., pi->id, pj->id, pi->viscosity.v_sig, pj->viscosity.v_sig, pi->decoupled, pi->feedback_data.decoupling_delay_time, dvdr, omega_ij, mu_ij, 1./r_inv);

  /* Need to get some kernel values F_ij = wi_dx */
  float wi, wi_dx;

  const float ui = r / hi;

  kernel_deval(ui, &wi, &wi_dx);

  /* Calculate the shock limiter component */
  const float shock_ratio_i =
      pj->viscosity.tensor_norm > 0.f
          ? pj->viscosity.shock_indicator / pj->viscosity.tensor_norm
          : 0.f;

  /* Rennehan: Add in complicated shock weighting to reduce noise
   * Eq 29 Wadsley+'17 */
  const float distance_norm_j = 2.f * kernel_gamma * hj;
  float w_R_i = 0.f;
  if (r < distance_norm_j) {
    const float w_R_i_core = 1.f - (r / distance_norm_j);
    const float w_R_i_core2 = w_R_i_core * w_R_i_core;
    w_R_i = w_R_i_core2 * w_R_i_core2;
  }

  pi->viscosity.shock_limiter += pj->mass * shock_ratio_i * w_R_i;
  /* Rennehan: Normalize the complicated weights */
  pi->viscosity.shock_limiter_norm += pj->mass * w_R_i;

  /* Correction factors for kernel gradients */

  const float rho_inv_j = 1.f / pj->rho;

  if (r < hi && r < hj) {
    pi->weighted_neighbour_wcount += pj->mass * r2 * wi_dx * rho_inv_j * r_inv;
  }

  /* Gradient of the density field */
  for (int j = 0; j < 3; j++) {
    const float drho_ij = pi->rho - pj->rho;
    const float dx_ij = pi->x[j] - pj->x[j];

    pi->rho_gradient[j] +=
        pj->mass * drho_ij * dx_ij * wi_dx * r_inv;
  }
}

/**
 * @brief Force interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {

  const int decoupled_i = pi->decoupled;
  const int decoupled_j = pj->decoupled;

  if (decoupled_i && decoupled_j) return;

  /* Cosmological factors entering the EoMs */
  const float a2_Hubble = a * a * H;

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mj = pj->mass;
  const float mi = pi->mass;

  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  if (!decoupled_j) {
    kernel_deval(xi, &wi, &wi_dx);
  }
  else {
    wi = 0.f;
    wi_dx = 0.f;
  }
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  if (!decoupled_i) {
    kernel_deval(xj, &wj, &wj_dx);
  }
  else {
    wj = 0.f;
    wj_dx = 0.f;
  }
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;
  pj->force.h_dt -= mi * dvdr * r_inv / rhoi * wj_dr;

  /* Decoupled winds do not need any further calculations. */
  if (decoupled_i || decoupled_j) return;

  /* Includes the hubble flow term; not used for du/dt */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;
  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);

#ifndef hydro_props_default_mu_ij_softening
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */
#else
  /* Equation 4.2 Monaghan 1992 */
  const float h_ij = 0.5f * (hi + hj);
  const float eta_ij = hydro_props_default_mu_ij_softening * h_ij;
  const float mu_ij = h_ij * omega_ij / (r2 + eta_ij * eta_ij);
#endif

  /* Variable smoothing length term */
  const float kernel_gradient =
      0.5f * (wi_dr * pi->force.f + wj_dr * pj->force.f);

  /* Construct the full viscosity term */
  const float rho_ij = rhoi + rhoj;
  const float alpha = pi->viscosity.alpha + pj->viscosity.alpha;
  const float visc =
      omega_ij < 0.f
          ? (-0.25f * alpha * (pi->force.soundspeed + pj->force.soundspeed) *
                 mu_ij +
             const_viscosity_beta_mu * mu_ij * mu_ij) /
                (0.5f * rho_ij)
          : 0.f;

  /* Convolve with the kernel */
  const float visc_acc_term = visc * kernel_gradient * r_inv;

  /* SPH acceleration term */
  const float sph_acc_term =
      (pressurei + pressurej) * r_inv * kernel_gradient / (pi->rho * pj->rho);

  /* Adaptive softening acceleration term */
  const float adapt_soft_acc_term = adaptive_softening_get_acc_term(
      pi, pj, wi_dr, wj_dr, pi->force.f, pj->force.f, r_inv);

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term + adapt_soft_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  pj->a_hydro[0] += mi * acc * dx[0];
  pj->a_hydro[1] += mi * acc * dx[1];
  pj->a_hydro[2] += mi * acc * dx[2];

  /* Get the time derivative for u. */
  const float sph_du_term_i =
      pressurei * dvdr * r_inv * kernel_gradient / (pi->rho * pj->rho);
  const float sph_du_term_j =
      pressurej * dvdr * r_inv * kernel_gradient / (pi->rho * pj->rho);

  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;

  /* Diffusion term */
  /*const float diff_du_term = 2.f * (pi->diffusion.rate + pj->diffusion.rate) *
                             (pi->u - pj->u) * kernel_gradient / rho_ij;*/
  /* Rennehan: replacing with Monaghan, Huppert, & Worster (2006) Eq 2.14 */
  float diff_du_term = 0.f;
  if (pi->diffusion.rate > 0.f && pj->diffusion.rate > 0.f && 
      pi->rho > 0.f && pj->rho > 0.f) {
    const float D_i_weighted = pi->rho * pi->diffusion.rate;
    const float D_j_weighted = pj->rho * pj->diffusion.rate;
    const float rho_ij_multiplied = pi->rho * pj->rho;
    const float diff_rate_ij = 
        2.f * (D_i_weighted * D_j_weighted) / (D_i_weighted + D_j_weighted);
    const float diff_weight = 2.f * diff_rate_ij / rho_ij_multiplied;
    diff_du_term = diff_weight * (pi->u - pj->u) * kernel_gradient;
  }

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term + diff_du_term;
  const float du_dt_j = sph_du_term_j + visc_du_term - diff_du_term;

  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;

  assert(pi->u_dt == pi->u_dt);
  assert(pj->u_dt == pj->u_dt);
  assert(pi->a_hydro[0] == pi->a_hydro[0]);
  assert(pj->a_hydro[0] == pj->a_hydro[0]);
}

/**
 * @brief Force interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, const struct part* restrict pj, const float a,
    const float H) {

  const int decoupled_i = pi->decoupled;
  const int decoupled_j = pj->decoupled;

  if (decoupled_i && decoupled_j) return;

  /* Cosmological factors entering the EoMs */
  const float a2_Hubble = a * a * H;

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mj = pj->mass;

  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  if (!decoupled_j) {
    kernel_deval(xi, &wi, &wi_dx);
  }
  else {
    wi = 0.f;
    wi_dx = 0.f;
  }
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  if (!decoupled_i) {
    kernel_deval(xj, &wj, &wj_dx);
  }
  else {
    wj = 0.f;
    wj_dx = 0.f;
  }
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;
  
  /* Decoupled winds do not need any further calculations. */
  if (decoupled_i || decoupled_j) return;

  /* Includes the hubble flow term; not used for du/dt */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);

#ifndef hydro_props_default_mu_ij_softening
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */
#else
  /* Equation 4.2 Monaghan 1992 */
  const float h_ij = 0.5f * (hi + hj);
  const float eta_ij = hydro_props_default_mu_ij_softening * h_ij;
  const float mu_ij = h_ij * omega_ij / (r2 + eta_ij * eta_ij);
#endif

  /* Variable smoothing length term */
  const float kernel_gradient =
      0.5f * (wi_dr * pi->force.f + wj_dr * pj->force.f);

  /* Construct the full viscosity term */
  const float rho_ij = rhoi + rhoj;
  const float alpha = pi->viscosity.alpha + pj->viscosity.alpha;
  const float visc =
      omega_ij < 0.f
          ? (-0.25f * alpha * (pi->force.soundspeed + pj->force.soundspeed) *
                 mu_ij +
             const_viscosity_beta_mu * mu_ij * mu_ij) /
                (0.5f * rho_ij)
          : 0.f;

  /* Convolve with the kernel */
  const float visc_acc_term = visc * kernel_gradient * r_inv;

  /* SPH acceleration term */
  const float sph_acc_term =
      (pressurei + pressurej) * r_inv * kernel_gradient / (pi->rho * pj->rho);

  /* Adaptive softening acceleration term */
  const float adapt_soft_acc_term = adaptive_softening_get_acc_term(
      pi, pj, wi_dr, wj_dr, pi->force.f, pj->force.f, r_inv);

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term + adapt_soft_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  /* Get the time derivative for u. */
  const float sph_du_term_i =
      pressurei * dvdr * r_inv * kernel_gradient / (pi->rho * pj->rho);

  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;

  /* Diffusion term */
  /*const float diff_du_term = 2.f * (pi->diffusion.rate + pj->diffusion.rate) *
                             (pi->u - pj->u) * kernel_gradient / rho_ij;*/
  /* Rennehan: replacing with Monaghan, Huppert, & Worster (2006) Eq 2.14 */
  float diff_du_term = 0.f;
  if (pi->diffusion.rate > 0.f && pj->diffusion.rate > 0.f && 
      pi->rho > 0.f && pj->rho > 0.f) {
    const float D_i_weighted = pi->rho * pi->diffusion.rate;
    const float D_j_weighted = pj->rho * pj->diffusion.rate;
    const float rho_ij_multiplied = pi->rho * pj->rho;
    const float diff_rate_ij = 
        2.f * (D_i_weighted * D_j_weighted) / (D_i_weighted + D_j_weighted);
    const float diff_weight = 2.f * diff_rate_ij / rho_ij_multiplied;
    diff_du_term = diff_weight * (pi->u - pj->u) * kernel_gradient;
  }

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term + diff_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;

  assert(pi->u_dt == pi->u_dt);
  assert(pj->u_dt == pj->u_dt);
  assert(pi->a_hydro[0] == pi->a_hydro[0]);
  assert(pj->a_hydro[0] == pj->a_hydro[0]);
  assert(pi->weighted_wcount!=0.f);
  assert(pj->weighted_wcount!=0.f);
}

#endif /* SWIFT_GASOLINE_HYDRO_IACT_H */
