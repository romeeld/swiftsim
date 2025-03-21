/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_CHEMISTRY_KIARA_H
#define SWIFT_CHEMISTRY_KIARA_H

/**
 * @file src/chemistry/KIARA/chemistry.h
 * @brief Empty infrastructure for the cases without chemistry function
 */

/* Some standard headers. */
#include <float.h>
#include <math.h>
#include <signal.h>

/* Local includes. */
#include "chemistry_struct.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/**
 * @brief Initializes summed particle quantities for the firehose wind model
 *
 * This is called from chemistry_init_part.
 *
 * @param p The particle to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void firehose_init_ambient_quantities(
    struct part* restrict p, const struct chemistry_global_data* cd) {

  struct chemistry_part_data* cpd = &p->chemistry_data;

  cpd->rho_ambient = 0.f;
  cpd->u_ambient = 0.f;
  cpd->w_ambient = 0.f;
}

/**
 * @brief Finishes up ambient quantity calculation for the firehose wind model
 *
 * This is called from chemistry_end_density
 *
 * @param p The particle to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void firehose_end_ambient_quantities(
    struct part* restrict p, const struct chemistry_global_data* cd,
    const struct cosmology* cosmo) {

  /* Limit ambient properties */
  p->chemistry_data.rho_ambient = min(
    p->chemistry_data.rho_ambient, 
    cd->firehose_ambient_rho_max * cosmo->a * cosmo->a * cosmo->a
  );

  if (p->chemistry_data.w_ambient >= 0.f) {
    p->chemistry_data.u_ambient = max(
      p->chemistry_data.u_ambient / p->chemistry_data.w_ambient, 
      cd->firehose_u_floor / cosmo->a_factor_internal_energy
    );
  }
  else {
    p->chemistry_data.u_ambient = 
        cd->firehose_u_floor / cosmo->a_factor_internal_energy;
  }
  //message("FIREHOSE_lim: %lld %g %g %g %g\n",p->id, p->chemistry_data.rho_ambient, rho_amb_limit, p->chemistry_data.u_ambient, T_floor * cd->temp_to_u_factor);
}


__attribute__((always_inline)) INLINE static void
logger_windprops_printprops(
    struct part *pi,
    const struct cosmology *cosmo, const struct chemistry_global_data* cd,
    FILE *fp) {

#ifdef FIREHOSE_DEBUG_CHECKS
  // Ignore COUPLED particles 
  if (pi->feedback_data.decoupling_delay_time <= 0.f) return;
  
  if (pi->id % 100 != 0) return;
  
  // Print wind properties
  const float length_convert = cosmo->a * cd->length_to_kpc;
  const float velocity_convert = cosmo->a_inv / cd->kms_to_internal;
  const float rho_convert = cosmo->a3_inv * cd->rho_to_n_cgs;
  const float u_convert =
      cosmo->a_factor_internal_energy / cd->temp_to_u_factor;
  /*fprintf(fp, "%.3f %lld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d\n",*/
  message("FIREHOSE: %.3f %lld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d\n",
        cosmo->z,
        pi->id,
        pi->gpart->fof_data.group_mass * cd->mass_to_solar_mass, 
        pi->h * cosmo->a * cd->length_to_kpc,
        pi->x[0] * length_convert,
        pi->x[1] * length_convert,
        pi->x[2] * length_convert,
        pi->v_full[0] * velocity_convert,
        pi->v_full[1] * velocity_convert,
        pi->v_full[2] * velocity_convert,
        pi->u * u_convert,
        pi->rho * rho_convert,
        pi->chemistry_data.radius_stream * length_convert,
        pi->chemistry_data.metal_mass_fraction_total,
        pi->viscosity.v_sig * velocity_convert,
        pi->feedback_data.decoupling_delay_time * cd->time_to_Myr,
        pi->feedback_data.number_of_times_decoupled);
#endif

  return;
}


/**
 * @brief Return a string containing the name of a given #chemistry_element.
 */
__attribute__((always_inline)) INLINE static const char*
chemistry_get_element_name(enum chemistry_element elem) {

  static const char* chemistry_element_names[chemistry_element_count] = {
      "Hydrogen", "Helium",    "Carbon",  "Nitrogen", "Oxygen",
      "Neon",     "Magnesium", "Silicon", "Sulfur", "Calcium", "Iron"};

  return chemistry_element_names[elem];
}

/**
 * @brief Prepares a particle for the smooth metal calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various smooth metallicity tasks
 *
 * @param p The particle to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void chemistry_init_part(
    struct part* restrict p, const struct chemistry_global_data* cd) {

  struct chemistry_part_data* cpd = &p->chemistry_data;

  /* Reset the shear tensor */
  for (int i = 0; i < 3; i++) {
    cpd->shear_tensor[i][0] = 0.f;
    cpd->shear_tensor[i][1] = 0.f;
    cpd->shear_tensor[i][2] = 0.f;
  }

  /* Reset the diffusion. */
  cpd->diffusion_coefficient = 0.f;

  /* Reset the change in metallicity */
  cpd->dZ_dt_total = 0.f;
  for (int elem = 0; elem < chemistry_element_count; ++elem) cpd->dZ_dt[elem] = 0.f;

#if COOLING_GRACKLE_MODE >= 2
  cpd->local_sfr_density = 0.f;
#endif
  firehose_init_ambient_quantities(p, cd);
}

/**
 * @brief Finishes the smooth metal calculation.
 *
 * Multiplies the metallicity and number of neighbours by the
 * appropiate constants and add the self-contribution term.
 *
 * This function requires the #hydro_end_density to have been called.
 *
 * @param p The particle to act upon.
 * @param cd #chemistry_global_data containing chemistry informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void chemistry_end_density(
    struct part* restrict p, const struct chemistry_global_data* cd,
    const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h; /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */

  struct chemistry_part_data* cpd = &p->chemistry_data;

  /* If diffusion is on, finish up shear tensor & particle's diffusion coeff */
  if (cd->diffusion_flag == 1 && cd->C_Smagorinsky > 0.f) {
    const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */
    const float rho = hydro_get_comoving_density(p);

    /* convert the shear factor into physical */
    const float factor_shear = h_inv_dim_plus_one * cosmo->a2_inv / rho;
    for (int k = 0; k < 3; k++) {
      cpd->shear_tensor[k][0] *= factor_shear;
      cpd->shear_tensor[k][1] *= factor_shear;
      cpd->shear_tensor[k][2] *= factor_shear;
    }

    /* Compute the trace over 3 and add the hubble flow. */
    float trace_3 = 0.f;
    for (int i = 0; i < 3; i++) {
      cpd->shear_tensor[i][i] += cosmo->H;
      trace_3 += cpd->shear_tensor[i][i];
    }
    trace_3 /= 3.f;

    float shear_tensor[3][3] = {{0.f}};
    for (int i = 0; i < 3; i++) {
      /* Make the tensor symmetric. */
      float avg = 0.5f * (cpd->shear_tensor[i][0] + cpd->shear_tensor[0][i]);
      shear_tensor[i][0] = avg;
      shear_tensor[0][i] = avg;

      avg = 0.5f * (cpd->shear_tensor[i][1] + cpd->shear_tensor[1][i]);
      shear_tensor[i][1] = avg;
      shear_tensor[1][i] = avg;

      avg = 0.5f * (cpd->shear_tensor[i][2] + cpd->shear_tensor[2][i]);
      shear_tensor[i][2] = avg;
      shear_tensor[2][i] = avg;

      /* Remove the trace. */
      shear_tensor[i][i] -= trace_3;
    }

    /* Compute the norm. */
    float velocity_gradient_norm = 0.f;
    for (int i = 0; i < 3; i++) {
      velocity_gradient_norm += shear_tensor[i][0] * shear_tensor[i][0];
      velocity_gradient_norm += shear_tensor[i][1] * shear_tensor[i][1];
      velocity_gradient_norm += shear_tensor[i][2] * shear_tensor[i][2];

      /* Copy the final values into the particle quantity */
      cpd->shear_tensor[i][0] = shear_tensor[i][0];
      cpd->shear_tensor[i][1] = shear_tensor[i][1];
      cpd->shear_tensor[i][2] = shear_tensor[i][2];
    }

    velocity_gradient_norm = sqrtf(2.f * velocity_gradient_norm);

    /* Never set D for cooling shutoff, wind, or ISM particles */
    if (!(p->feedback_data.cooling_shutoff_delay_time > 0.f) &&
        !(p->feedback_data.decoupling_delay_time > 0.f) &&
        !(p->cooling_data.subgrid_temp > 0.f)) {

      /* Compute the diffusion coefficient in physical coordinates.
      * The norm is already in physical coordinates.
      * kernel_gamma is necessary (see Rennehan 2021)
      */
      const float rho_phys = hydro_get_physical_density(p, cosmo);
      const float h_phys = cosmo->a * p->h * kernel_gamma;
      const float smag_length_scale = cd->C_Smagorinsky * h_phys;
      float D_phys = rho_phys * smag_length_scale * smag_length_scale * 
                         velocity_gradient_norm;

      /* Sometimes, the diffusion coefficient can be quite large
       * and actually be the limiting time step for the particle. That is 
       * because the diffusion coefficient can be something like
       * D ~ 1e8 Msun / (Myr * kpc) which is basically like transporting more
       * than a particle's mass in a single standard timestep.
       * For that reason, we impose an upper limit on the diffusion coefficient.
       * 
       * The mass transfer estimate is D * h * dt, so we can simplify
       * since dt = Beta * rho * h^2 / D and D = rho * (C_s * h)^2 * |S|
       * to:
       *    dM = D * h * dt
       *    dM = D * h * (Beta * rho * h^2 / D)
       *    dM = Beta * rho * h^3
       * 
       * If dM > dM_max, then we need to reset D such that the mass transfer
       * rate is satisfied. Compute dM / M with a maximum value of of 1. Then,
       * multiply by Beta to reduce the diffusion rate.
       * 
       */
      const float dm = 
          cd->diffusion_beta * rho_phys * h_phys * h_phys * h_phys;
      const float m = hydro_get_mass(p);
      const float dm_m = m > 0.f ? fmin(dm / m, 1.f) : 0.f;

      if (dm_m > cd->diffusion_beta) {
        D_phys *= cd->diffusion_beta;
      }

      cpd->diffusion_coefficient = D_phys;
    }
  }

#if COOLING_GRACKLE_MODE >= 2
  /* Add self contribution to SFR */
  cpd->local_sfr_density += max(0.f, p->sf_data.SFR);
  const float vol_factor = h_inv_dim / (4.f * M_PI / 3.f);
  /* Convert to physical density from comoving */
  cpd->local_sfr_density *= vol_factor * cosmo->a3_inv;
#endif
  if (cd->use_firehose_wind_model) {
    firehose_end_ambient_quantities(p, cd, cosmo);
  }
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_has_no_neighbours(struct part* restrict p,
                                 struct xpart* restrict xp,
                                 const struct chemistry_global_data* cd,
                                 const struct cosmology* cosmo) {

  /* Just make all the smoothed fields default to the un-smoothed values */
  struct chemistry_part_data* cpd = &p->chemistry_data;
  /* Reset the shear tensor */
  for (int i = 0; i < 3; i++) {
    cpd->shear_tensor[i][0] = 0.f;
    cpd->shear_tensor[i][1] = 0.f;
    cpd->shear_tensor[i][2] = 0.f;
  }

  /* Reset the diffusion. */
  cpd->diffusion_coefficient = 0.f;

  /* Reset the change in metallicity */
  cpd->dZ_dt_total = 0.f;
  for (int elem = 0; elem < chemistry_element_count; ++elem) cpd->dZ_dt[elem] = 0.f;

#if COOLING_GRACKLE_MODE >= 2
  /* If there is no nearby SF, set to zero */
  cpd->local_sfr_density = 0.f;
#endif
}

/**
 * @brief Sets the chemistry properties of the (x-)particles to a valid start
 * state.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param data The global chemistry information.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_part(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct chemistry_global_data* data, struct part* restrict p,
    struct xpart* restrict xp) {

  /* Initialize mass fractions for total metals and each metal individually */
  if (data->initial_metal_mass_fraction_total != -1) {
    p->chemistry_data.metal_mass_fraction_total =
        data->initial_metal_mass_fraction_total;

    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      p->chemistry_data.metal_mass_fraction[elem] =
          data->initial_metal_mass_fraction[elem];
    }
  }
  chemistry_init_part(p, data);
  firehose_init_ambient_quantities(p, data);
}

/**
 * @brief Sets the chemistry properties of the sparticles to a valid start
 * state.
 *
 * @param data The global chemistry information.
 * @param sp Pointer to the sparticle data.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_spart(
    const struct chemistry_global_data* data, struct spart* restrict sp) {

  /* Initialize mass fractions for total metals and each metal individually */
  if (data->initial_metal_mass_fraction_total != -1) {
    sp->chemistry_data.metal_mass_fraction_total =
        data->initial_metal_mass_fraction_total;

    for (int elem = 0; elem < chemistry_element_count; ++elem)
      sp->chemistry_data.metal_mass_fraction[elem] =
          data->initial_metal_mass_fraction[elem];
  }
}

/**
 * @brief Sets the chemistry properties of the sink particles to a valid start
 * state.
 *
 * @param data The global chemistry information.
 * @param sink Pointer to the sink particle data.
 * Required by space_first_init.c
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_sink(
    const struct chemistry_global_data* data, struct sink* restrict sink) {}


/**
 * @brief Initialises the chemistry properties.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param data The properties to initialise.
 */
static INLINE void chemistry_init_backend(struct swift_params* parameter_file,
                                          const struct unit_system* us,
                                          const struct phys_const* phys_const,
                                          struct chemistry_global_data* data) {

  /* Set some useful unit conversions */
  const double Msun_cgs = phys_const->const_solar_mass *
                          units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  const double unit_mass_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  data->mass_to_solar_mass = unit_mass_cgs / Msun_cgs;
  data->temp_to_u_factor = 
      phys_const->const_boltzmann_k / 
        (hydro_gamma_minus_one * phys_const->const_proton_mass *
          units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE));
  const double X_H = 0.75;
  data->rho_to_n_cgs =
      (X_H / phys_const->const_proton_mass) * 
        units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);
  data->kms_to_internal = 1.0e5f / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);
  data->time_to_Myr = units_cgs_conversion_factor(us, UNIT_CONV_TIME) /
      (1.e6f * 365.25f * 24.f * 60.f * 60.f);
  data->length_to_kpc =
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH) / 3.08567758e21f;

  /* Is metal diffusion turned on? */
  data->diffusion_flag = parser_get_param_int(parameter_file,
                                              "KIARAChemistry:diffusion_on");

  /* Read the diffusion coefficient */
  data->C_Smagorinsky = 
      parser_get_opt_param_float(parameter_file,
                                 "KIARAChemistry:diffusion_coefficient",
                                 0.23f);

  /* Time-step restriction to <= 0.15*rho*h^2 / D from 
     Parshikov & Medin 2002 equation 41 */
  data->diffusion_beta = 
      parser_get_opt_param_float(parameter_file,
                                 "KIARAChemistry:diffusion_beta",
                                 0.1f);
  if (data->diffusion_beta < 0.f || data->diffusion_beta > 0.1f) {
    error("diffusion_beta must be >= 0 and <= 0.1");
  }

  data->time_step_min =
      parser_get_opt_param_float(parameter_file,
                                 "KIARAChemistry:minimum_timestep_Myr",
                                 0.1f);
  data->time_step_min /= data->time_to_Myr;
  if (data->time_step_min < 0.f) {
    error("time_step_min must be > 0");
  }

  data->max_fractional_Z_transfer = 
      parser_get_opt_param_float(parameter_file,
                                 "KIARAChemistry:max_fractional_Z_transfer",
                                 0.25f);
  if (data->max_fractional_Z_transfer < 0.f || 
          data->max_fractional_Z_transfer > 1.f) {
    error("diffusion_beta must be >= 0 and <= 1");
  }

  /* Are we using the firehose wind model? */
  data->use_firehose_wind_model = 
      parser_get_opt_param_int(parameter_file,
                               "KIARAChemistry:use_firehose_wind_model", 
                               0);

  /* Firehose model parameters */
  data->firehose_ambient_rho_max = 
      parser_get_opt_param_float(parameter_file,
                                 "KIARAChemistry:firehose_nh_ambient_max_cgs", 
                                 0.1f);
  data->firehose_ambient_rho_max /= data->rho_to_n_cgs;

  data->firehose_u_floor = 
      parser_get_opt_param_float(parameter_file, 
                                 "KIARAChemistry:firehose_temp_floor", 
                                 1.e4f);
  data->firehose_u_floor *= data->temp_to_u_factor;

  /* Firehose recoupling parameters */
  data->firehose_recoupling_mach = 
      parser_get_opt_param_float(parameter_file,
                                 "KIARAChemistry:firehose_recoupling_mach", 
                                 0.5f);

  data->firehose_recoupling_u_factor = 
      parser_get_opt_param_float(parameter_file,
                                 "KIARAChemistry:firehose_recoupling_u_factor", 
                                 0.5f);

  data->firehose_recoupling_fmix = 
      parser_get_opt_param_float(parameter_file,
                                 "KIARAChemistry:firehose_recoupling_fmix", 
                                 0.9f);

  data->firehose_max_velocity = 
      parser_get_opt_param_float(parameter_file,
                                 "KIARAChemistry:firehose_max_velocity", 
                                 2000.f);
  data->firehose_max_velocity *= data->kms_to_internal;

  /* Read the total metallicity */
  data->initial_metal_mass_fraction_total = parser_get_opt_param_float(
      parameter_file, "KIARAChemistry:init_abundance_metal", -1);

  if (data->initial_metal_mass_fraction_total != -1) {
    /* Read the individual mass fractions */
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      char buffer[50];
      sprintf(buffer, "KIARAChemistry:init_abundance_%s",
              chemistry_get_element_name((enum chemistry_element)elem));

      data->initial_metal_mass_fraction[elem] =
          parser_get_param_float(parameter_file, buffer);
    }

    /* Let's check that things make sense (broadly) */

    /* H + He + Z should be ~1 */
    float total_frac = data->initial_metal_mass_fraction[chemistry_element_H] +
                       data->initial_metal_mass_fraction[chemistry_element_He] +
                       data->initial_metal_mass_fraction_total;

    if (total_frac < 0.98 || total_frac > 1.02)
      error("The abundances provided seem odd! H + He + Z = %f =/= 1.",
            total_frac);

    /* Sum of metal elements should be <= Z */
    total_frac = 0.f;
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      if (elem != chemistry_element_H && elem != chemistry_element_He) {
        total_frac += data->initial_metal_mass_fraction[elem];
      }
    }

    if (!data->diffusion_flag) {
      if (total_frac > 1.02 * data->initial_metal_mass_fraction_total) {
        error(
            "The abundances provided seem odd! \\sum metal elements (%f) > Z "
            "(%f)",
            total_frac, data->initial_metal_mass_fraction_total);
      }
    }
    else {
      /* If diffusion is on, need a metallicity floor so reset Z total */
      if (total_frac > 1.02 * data->initial_metal_mass_fraction_total) {
        warning("Resetting total Z to the sum of all available metals.");
        data->initial_metal_mass_fraction_total = total_frac;

        /* H + He + Z should be ~1 */
        float total_frac_check = 
            data->initial_metal_mass_fraction[chemistry_element_H] +
            data->initial_metal_mass_fraction[chemistry_element_He] +
            data->initial_metal_mass_fraction_total;

        if (total_frac_check < 0.98 || total_frac_check > 1.02) {
          error("After resetting, the abundances provided seem odd! "
                "H + He + Z = %f =/= 1.",
                total_frac_check);
        }
      }
    }

    /* Sum of all elements should be <= 1 */
    total_frac = 0.f;
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      total_frac += data->initial_metal_mass_fraction[elem];
    }

    if (total_frac > 1.02) {
      error("The abundances provided seem odd! \\sum elements (%f) > 1",
            total_frac);
    }
  }

}

/**
 * @brief Prints the properties of the chemistry model to stdout.
 *
 * @brief The #chemistry_global_data containing information about the current
 * model.
 */
static INLINE void chemistry_print_backend(
    const struct chemistry_global_data* data) {

  if (data->use_firehose_wind_model) {
    if (data->diffusion_flag) {
      message("Chemistry model is 'KIARA' tracking %d elements with the "
              "Firehose wind model and metal diffusion.",
              chemistry_element_count);
    }
    else {
      message("Chemistry model is 'KIARA' tracking %d elements with "
              "the Firehose wind model.",
              chemistry_element_count);
    }
  }
  else {
    if (data->diffusion_flag) {
      message("Chemistry model is 'KIARA' tracking %d elements with "
              " metal diffusion on.",
              chemistry_element_count);
    }
    else {
      message("Chemistry model is 'KIARA' tracking %d elements.",
            chemistry_element_count);
    }
  }
}

/**
 * @brief Updates to the chemistry data after the hydro force loop.
 *
 * Finish off the diffusion by actually exchanging the metals
 *
 * @param p The particle to act upon.
 * @param cosmo The current cosmological model.
 * @param with_cosmology Are we running with the cosmology?
 * @param time Current time of the simulation.
 * @param dt Time step (in physical units).
 */
__attribute__((always_inline)) INLINE static void chemistry_end_force(
    struct part* restrict p, const struct cosmology* cosmo,
    const int with_cosmology, const double time, 
    const struct chemistry_global_data* cd, const double dt) {

  /* never do diffusion for a cooling shutoff or wind particle */
  if (p->feedback_data.cooling_shutoff_delay_time > 0.f ||
        p->feedback_data.decoupling_delay_time > 0.f) return;

  if (dt == 0.) return;

  /* Check if we are hypersonic*/
  /* Reset dZ_dt and return? */
  bool reset_time_derivatives = false;

  struct chemistry_part_data* ch = &p->chemistry_data;

  const float h_inv = 1.f / p->h;
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */
  /* Missing factors in iact. */
  const float factor = h_inv_dim * h_inv;

  /* Add diffused metals to particle */
  const float dZ_tot = ch->dZ_dt_total * dt * factor;
  const float new_metal_mass_fraction_total 
                  = ch->metal_mass_fraction_total + dZ_tot;
  if (ch->metal_mass_fraction_total > 0.f) {
    const float abs_fractional_change = 
        fabs(dZ_tot) / ch->metal_mass_fraction_total;
    /* Check if dZ is bigger than 1/4 of the Z */
    if (abs_fractional_change > cd->max_fractional_Z_transfer) {
      reset_time_derivatives = true;
    }
  }

  /* Handle edge case where diffusion leads to negative metallicity */
  if (new_metal_mass_fraction_total < 0.f) {
    warning("Metal diffusion led to negative metallicity!\n"
            "\tpid=%lld\n\tdt=%g\n\tZ=%g\n\tdZ_dt=%g\n"
            "\tdZtot=%g\n\tZnewtot=%g\n\tfactor=%g",
            p->id,
            dt,
            ch->metal_mass_fraction_total,
            ch->dZ_dt_total,
            dZ_tot,
            new_metal_mass_fraction_total,
            factor);
    reset_time_derivatives = true;
  }

  /* Handle edge case where diffusion leads to super-unity metallicity */
  if (new_metal_mass_fraction_total > 1.f) {
    warning("Metal diffusion led to metal fractions above unity!\n"
            "pid=%lld\n\tdt=%g\n\tZ=%g\n\tdZ_dt=%g\n"
            "\tdZtot=%g\n\tZnewtot=%g\n\tfactor=%g",
            p->id,
            dt,
            ch->metal_mass_fraction_total,
            ch->dZ_dt_total,
            dZ_tot,
            new_metal_mass_fraction_total,
            factor);
    reset_time_derivatives = true;
  }


  /* Add individual element contributions from diffusion */
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    const float dZ = ch->dZ_dt[elem] * dt * factor;
    const float new_metal_fraction_elem
                    = ch->metal_mass_fraction[elem] + dZ;

    if (ch->metal_mass_fraction[elem] > 0.f) {
      const float abs_fractional_change = 
          fabs(dZ) / ch->metal_mass_fraction[elem];
      if (abs_fractional_change > cd->max_fractional_Z_transfer) {
        reset_time_derivatives = true;
      }
    }

    /* Make sure that the metallicity is 0 <= x <= 1 */
    if (new_metal_fraction_elem < 0.f) {
      warning("Z[elem] < 0! pid=%lld, dt=%g, elem=%d, Z=%g, dZ_dt=%g, dZ=%g, "
              "dZtot=%g Ztot=%g Zdust=%g.", 
              p->id, 
              dt, 
              elem, 
              ch->metal_mass_fraction[elem], 
              ch->dZ_dt[elem], 
              dZ,
              dZ_tot, 
              ch->metal_mass_fraction_total, 
              p->cooling_data.dust_mass_fraction[elem]);
      reset_time_derivatives = true;
    }

    if (new_metal_fraction_elem > 1.f) {
      warning("Z[elem] > 1! pid=%lld, dt=%g, elem=%d, Z=%g, dZ_dt=%g, "
              "dZ=%g, dZtot=%g Ztot=%g.", 
              p->id, 
              dt, 
              elem, 
              ch->metal_mass_fraction[elem], 
              ch->dZ_dt[elem], 
              dZ,
              dZ_tot, 
              ch->metal_mass_fraction_total);
      reset_time_derivatives = true;
    }
  }

  /* Found weird dZ_dt values so we should reset everything and exit */
  if (reset_time_derivatives) {
    ch->dZ_dt_total = 0.f;
    for (int elem = 0; elem < chemistry_element_count; elem++) {
      ch->dZ_dt[elem] = 0.f;
    }
    return;
  }
  else {
#if COOLING_GRACKLE_MODE >= 2
    if (ch->metal_mass_fraction_total > 0.f) {
      /* Add diffused dust to particle, in proportion to added metals */
      p->cooling_data.dust_mass *= 1.f + dZ_tot / ch->metal_mass_fraction_total;
    }
#endif

    /* Reset the total metallicity Z */
    ch->metal_mass_fraction_total = new_metal_mass_fraction_total;

    /* Add individual element contributions from diffusion */
    for (int elem = 0; elem < chemistry_element_count; elem++) {
      const float dZ = ch->dZ_dt[elem] * dt * factor;
      const float new_metal_fraction_elem = 
          ch->metal_mass_fraction[elem] + dZ;
                      
  #if COOLING_GRACKLE_MODE >= 2
    /* Add diffused dust to particle, in proportion to added metals */
      if (ch->metal_mass_fraction[elem] > 0.f) {
        p->cooling_data.dust_mass_fraction[elem] *= 
            1.f + dZ / ch->metal_mass_fraction[elem];
      }
  #endif

      /* Treating Z like a passive scalar */
      ch->metal_mass_fraction[elem] = new_metal_fraction_elem;
    }
  }
}

/**
 * @brief Computes the chemistry-related time-step constraint.
 *
 * Only constraint in KIARA is the diffusion time-step.
 *
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param us The internal system of units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cd The global properties of the chemistry scheme.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float chemistry_timestep(
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_props,
    const struct chemistry_global_data* cd, const struct part* restrict p) {

  float dt_chem = FLT_MAX;
  if (cd->diffusion_flag) {
    if (p->chemistry_data.diffusion_coefficient > 0.f) {
      const struct chemistry_part_data *ch = &p->chemistry_data;

      /* Parshikov & Medin 2002 equation 41 */
      const float h_phys = p->h * cosmo->a * kernel_gamma;
      const float D_phys = ch->diffusion_coefficient;
      const float rho_phys = hydro_get_physical_density(p, cosmo);
      dt_chem = cd->diffusion_beta * rho_phys * h_phys * h_phys / D_phys;
      if (dt_chem < cd->time_step_min) {
        message(
          "Warning! dZ_dt timestep low: id=%lld (%g Myr) is below "
          "time_step_min (%g Myr).",
          p->id, dt_chem * cd->time_to_Myr,
          cd->time_step_min * cd->time_to_Myr);
      }
    }
  }

  return max(cd->time_step_min, dt_chem);
}

/**
 * @brief Initialise the chemistry properties of a black hole with
 * the chemistry properties of the gas it is born from.
 *
 * Black holes don't store fractions so we need to use element masses.
 *
 * @param bp_data The black hole data to initialise.
 * @param p_data The gas data to use.
 * @param gas_mass The mass of the gas particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_bpart_from_part(
    struct chemistry_bpart_data* bp_data,
    const struct chemistry_part_data* p_data, const double gas_mass) {

  bp_data->metal_mass_total = p_data->metal_mass_fraction_total * gas_mass;
  for (int i = 0; i < chemistry_element_count; ++i) {
    bp_data->metal_mass[i] = p_data->metal_mass_fraction[i] * gas_mass;
  }

  bp_data->formation_metallicity = p_data->metal_mass_fraction_total;
}

/**
 * @brief Add the chemistry data of a gas particle to a black hole.
 *
 * Black holes don't store fractions so we need to add element masses.
 *
 * @param bp_data The black hole data to add to.
 * @param p_data The gas data to use.
 * @param gas_mass The mass of the gas particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_add_part_to_bpart(
    struct chemistry_bpart_data* bp_data,
    const struct chemistry_part_data* p_data, const double gas_mass) {

  bp_data->metal_mass_total += p_data->metal_mass_fraction_total * gas_mass;
  for (int i = 0; i < chemistry_element_count; ++i) {
    bp_data->metal_mass[i] += p_data->metal_mass_fraction[i] * gas_mass;
  }
}

/**
 * @brief Transfer chemistry data of a gas particle to a black hole.
 *
 * In contrast to `chemistry_add_part_to_bpart`, only a fraction of the
 * masses stored in the gas particle are transferred here. Absolute masses
 * of the gas particle are adjusted as well.
 * Black holes don't store fractions so we need to add element masses.
 *
 * We expect the nibble_mass to be the gas particle mass multiplied by the
 * nibble_fraction.
 *
 * @param bp_data The black hole data to add to.
 * @param p_data The gas data to use.
 * @param nibble_mass The mass to be removed from the gas particle.
 * @param nibble_fraction The fraction of the (original) mass of the gas
 *        particle that is removed.
 */
__attribute__((always_inline)) INLINE static void
chemistry_transfer_part_to_bpart(struct chemistry_bpart_data* bp_data,
                                 struct chemistry_part_data* p_data,
                                 const double nibble_mass,
                                 const double nibble_fraction) {

  bp_data->metal_mass_total += p_data->metal_mass_fraction_total * nibble_mass;
  for (int i = 0; i < chemistry_element_count; ++i)
    bp_data->metal_mass[i] += p_data->metal_mass_fraction[i] * nibble_mass;

}

/**
 * @brief Add the chemistry data of a black hole to another one.
 *
 * @param bp_data The black hole data to add to.
 * @param swallowed_data The black hole data to use.
 */
__attribute__((always_inline)) INLINE static void chemistry_add_bpart_to_bpart(
    struct chemistry_bpart_data* bp_data,
    const struct chemistry_bpart_data* swallowed_data) {

  bp_data->metal_mass_total += swallowed_data->metal_mass_total;
  for (int i = 0; i < chemistry_element_count; ++i) {
    bp_data->metal_mass[i] += swallowed_data->metal_mass[i];
  }
}

/**
 * @brief Split the metal content of a particle into n pieces
 *
 * We only need to split the fields that are not fractions.
 *
 * @param p The #part.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void chemistry_split_part(
    struct part* p, const double n) { }

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in feedback/enrichment related routines.
 *
 * We return the un-smoothed quantity here as the star will smooth
 * over its gas neighbours.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_fraction_for_feedback(
    const struct part* restrict p) {

  return p->chemistry_data.metal_mass_fraction_total;
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in feedback/enrichment related routines.
 *
 * We return the un-smoothed quantity here as the star will smooth
 * over its gas neighbours.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float const*
chemistry_get_metal_mass_fraction_for_feedback(const struct part* restrict p) {

  return p->chemistry_data.metal_mass_fraction;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * star particle to be used in feedback/enrichment related routines.
 *
 * KIARA uses smooth abundances for everything.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_star_total_metal_mass_fraction_for_feedback(
    const struct spart* restrict sp) {

  return sp->chemistry_data.metal_mass_fraction_total;
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * star particle to be used in feedback/enrichment related routines.
 *
 * KIARA uses smooth abundances for everything.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float const*
chemistry_get_star_metal_mass_fraction_for_feedback(
    const struct spart* restrict sp) {

  return sp->chemistry_data.metal_mass_fraction;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in cooling related routines.
 *
 * KIARA uses smooth abundances for everything.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_fraction_for_cooling(
    const struct part* restrict p) {

  return p->chemistry_data.metal_mass_fraction_total;
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in cooling related routines.
 *
 * KIARA uses smooth abundances for everything.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float const*
chemistry_get_metal_mass_fraction_for_cooling(const struct part* restrict p) {

  return p->chemistry_data.metal_mass_fraction;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in star formation related routines.
 *
 * KIARA uses smooth abundances for everything.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_fraction_for_star_formation(
    const struct part* restrict p) {

  return p->chemistry_data.metal_mass_fraction_total;
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in star formation related routines.
 *
 * KIARA uses smooth abundances for everything.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float const*
chemistry_get_metal_mass_fraction_for_star_formation(
    const struct part* restrict p) {

  return p->chemistry_data.metal_mass_fraction;
}

/**
 * @brief Returns the total metal mass of the
 * gas particle to be used in the stats related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_for_stats(const struct part* restrict p) {

  return p->chemistry_data.metal_mass_fraction_total * hydro_get_mass(p);
}

/**
 * @brief Returns the total metal mass of the
 * star particle to be used in the stats related routines.
 *
 * @param sp Pointer to the star particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_star_total_metal_mass_for_stats(const struct spart* restrict sp) {

  return sp->chemistry_data.metal_mass_fraction_total * sp->mass;
}

/**
 * @brief Returns the total metal mass of the
 * black hole particle to be used in the stats related routines.
 *
 * @param bp Pointer to the BH particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_bh_total_metal_mass_for_stats(const struct bpart* restrict bp) {

  return bp->chemistry_data.metal_mass_total;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * star particle to be used in the luminosity calculations.
 *
 * @param sp Pointer to the star particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_star_total_metal_mass_fraction_for_luminosity(
    const struct spart* restrict sp) {

  return sp->chemistry_data.metal_mass_fraction_total;
}


#endif /* SWIFT_CHEMISTRY_KIARA_H */
