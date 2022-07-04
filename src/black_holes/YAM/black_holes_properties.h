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
#ifndef SWIFT_YAM_BLACK_HOLES_PROPERTIES_H
#define SWIFT_YAM_BLACK_HOLES_PROPERTIES_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "chemistry.h"
#include "exp10.h"
#include "hydro_properties.h"

/* Includes. */
#include <string.h>

/**
 * @brief Modes of energy injection for AGN feedback
 */
enum AGN_feedback_models {
  AGN_random_ngb_model,       /*< Random neighbour model for AGN feedback */
  AGN_isotropic_model,        /*< Isotropic model of AGN feedback */
  AGN_minimum_distance_model, /*< Minimum-distance model of AGN feedback */
  AGN_minimum_density_model   /*< Minimum-density model of AGN feedback */
};

enum BH_merger_thresholds {
  BH_mergers_circular_velocity,        /*< v_circ at separation, as in EAGLE */
  BH_mergers_escape_velocity,          /*< v_esc at separation */
  BH_mergers_dynamical_escape_velocity /*< combined v_esc_dyn at separation */
};

/**
 * @brief Properties of black holes and AGN feedback in the EAGEL model.
 */
struct black_holes_props {

  /* ----- Basic neighbour search properties ------ */

  /*! Resolution parameter */
  float eta_neighbours;

  /*! Target weightd number of neighbours (for info only)*/
  float target_neighbours;

  /*! Smoothing length tolerance */
  float h_tolerance;

  /*! Tolerance on neighbour number  (for info only)*/
  float delta_neighbours;

  /*! Maximal number of iterations to converge h */
  int max_smoothing_iterations;

  /*! Maximal change of h over one time-step */
  float log_max_h_change;

  /* ----- Initialisation properties  ------ */

  /*! Are we seeding from star formation? */
  int seed_during_star_formation;

  /*! If above this threshold, seed a SMBH from SF gas */
  float seed_n_H_threshold_cgs;

  /*! If below this threshold, seed a SMBH from SF gas */
  float seed_temperature_threshold_cgs;

  /*! If below this threshold, seed a SMBH from SF gas */
  float seed_metallicity_threshold;

  /*! Mass of a BH seed at creation time */
  float subgrid_seed_mass;

  /*! Should we use the subgrid mass specified in ICs? */
  int use_subgrid_mass_from_ics;

  /*! Should we enforce positive subgrid masses initially? */
  int with_subgrid_mass_check;
  
  /* ----- Properties of the accretion model ------ */

  /*! Radiative efficiency of the black holes. */
  float epsilon_r;

  /*! Maximal fraction of the Eddington rate allowed. */
  float f_Edd_maximum;

  /*! Lowest value of the boost of the Booth & Schaye 2009 model */
  float boost_alpha;

  /*! Minimum gas particle mass in nibbling mode */
  float min_gas_mass_for_nibbling;

  /*! Switch to calculate the sound speed with a fixed T near the EoS */
  int with_fixed_T_near_EoS;

  /*! Factor above EoS below which fixed T applies for sound speed */
  float fixed_T_above_EoS_factor;

  /*! Fixed T (expressed as internal energy) for sound speed near EoS */
  float fixed_u_for_soundspeed;

  /*! Where do we distinguish between hot gas for Bondi? */
  float environment_temperature_cut;

  /*! Normalization of the torque accretion rate */
  float torque_accretion_norm;

  /*! Factor in front of M/(dM/dt) for timestepping */
  float dt_accretion_factor;

  /*! A from Lupi+17 */
  float A_lupi;

  /*! B from Lupi+17 */
  float B_lupi;

  /*! C from Lupi+17 */
  float C_lupi;

  /*! The spin of EVERY black hole */
  float fixed_spin;

  /* ---- Properties of the feedback model ------- */

  /*! AGN feedback model: random, isotropic or minimum distance */
  enum AGN_feedback_models feedback_model;

  /*! Is the AGN feedback model deterministic or stochastic? */
  int AGN_deterministic;

  /*! Feedback coupling efficiency of the black holes. */
  float epsilon_f;

  /*! When does the jet start heating? (km/s) */
  float jet_heating_velocity_threshold;

  /*! At what v_kick does the X-ray heating start? */
  float xray_heating_velocity_threshold;

  /*! How many times the particle's u value can we heat? */
  float xray_maximum_heating_factor;

  /*! How much of the X-ray energy should go into velocity for dense gas? */
  float xray_kinetic_fraction;

  /*! Above this density we should split X-ray energy into kinetic */ 
  float xray_heating_n_H_threshold_cgs;

  /*! Below this temperature we should split X-ray energy into kinetic */
  float xray_heating_T_threshold_cgs;

  /*! What is the physical max. velocity of the jet? (km/s) */
  float jet_velocity;

  /*! The temperature of the jet. Set < 0.f for halo virial temperature */
  float jet_temperature;

  /*! What lower Mdot,BH/Mdot,Edd boundary does the jet activate? */
  float eddington_fraction_lower_boundary;

  /*! What upper Mdot,BH/Mdot,Edd boundary does the slim disk mode activate? */
  float eddington_fraction_upper_boundary;

  /*! Minimum mass for starting the jet (Msun) */
  float jet_mass_min_Msun;

  /*! Maximum mass for starting the jet (Msun) */
  float jet_mass_max_Msun;

  /*! Constrains momentum of outflowing wind to p = F * L / c */
  float quasar_wind_momentum_flux;

  /*! The mass loading of the quasar outflow */
  float quasar_wind_mass_loading;

  /*! The wind speed of the quasar outflow */
  float quasar_wind_speed;

  /*! f_acc for the quasar mode */
  float quasar_f_accretion;

  /*! eps_f for the quasar mode */
  float quasar_coupling;

  /*! eps_f for the ADAF mode */
  float adaf_coupling;

  /*! f_acc for the ADAF mode */
  float adaf_f_accretion;

  /*! Wind momnetum flux for the ADAF mode */
  float adaf_wind_momentum_flux;

  /*! eps_f for the slim disk mode */
  float slim_disk_coupling;

  /*! Momentum flux for the slim disk mode */
  float slim_disk_wind_momentum_flux;

  /*! wind speed in the slim disk mode */
  float slim_disk_wind_speed;

  /*! Factor in front of E/(dE/dt) for timestepping. */
  float dt_feedback_factor;

  /*! The efficiency of the jet */
  float jet_efficiency;

  /*! The mass loading in the jet */
  float jet_loading;

  /*! The quadratic term (see paper) for the jet */
  float jet_quadratic_term;

  /* ---- Properties of the repositioning model --- */

  /*! Maximal mass of BH to reposition */
  float max_reposition_mass;

  /*! Maximal distance to reposition, in units of softening length */
  float max_reposition_distance_ratio;

  /*! Switch to enable a relative velocity limit for particles to which the
   * black holes can reposition */
  int with_reposition_velocity_threshold;

  /*! Maximal velocity offset of particles to which the black hole can
   * reposition, in units of the ambient sound speed of the black hole */
  float max_reposition_velocity_ratio;

  /*! Minimum value of the velocity repositioning threshold */
  float min_reposition_velocity_threshold;

  /*! Switch to enable repositioning at fixed (maximum) speed */
  int set_reposition_speed;

  /*! Normalisation factor for repositioning velocity */
  float reposition_coefficient_upsilon;

  /*! Reference black hole mass for repositioning scaling */
  float reposition_reference_mass;

  /*! Repositioning velocity scaling with black hole mass */
  float reposition_exponent_mass;

  /*! Reference gas density for repositioning scaling */
  float reposition_reference_n_H;

  /*! Repositioning velocity scaling with gas density */
  float reposition_exponent_n_H;

  /*! Correct potential of BH? */
  int correct_bh_potential_for_repositioning;

  /* ---- Properties of the merger model ---------- */

  /*! Mass ratio above which a merger is considered 'minor' */
  float minor_merger_threshold;

  /*! Mass ratio above which a merger is considered 'major' */
  float major_merger_threshold;

  /*! Type of merger threshold */
  enum BH_merger_thresholds merger_threshold_type;

  /*! Maximal distance over which BHs merge, in units of softening length */
  float max_merging_distance_ratio;

  /* ---- Black hole time-step properties ---------- */

  /*! Minimum allowed time-step of BH (internal units) */
  float time_step_min;

  /* ---- Common conversion factors --------------- */

  /*! Conversion factor from temperature to internal energy */
  float temp_to_u_factor;

  /*! Conversion factor from physical density to n_H [cgs] */
  float rho_to_n_cgs;

  /*! Conversion factor from internal mass to solar masses */
  float mass_to_solar_mass;

  /*! Conversion factor from km/s to internal velocity units (without a-factor) */
  float kms_to_internal;

  /*! Conversion factor from internal length to parsec */
  float length_to_parsec;

  /*! Conversion factor from internal time to yr */
  float time_to_yr;

  /*! Conversion factor from internal time to Myr */
  float time_to_Myr;

  /*! Conversion factor from density to cgs */
  double conv_factor_density_to_cgs;

  /*! Conversion factor from luminosity to cgs */
  double conv_factor_energy_rate_to_cgs;

  /*! Conversion factor from length to cgs */
  double conv_factor_length_to_cgs;

  /*! Conversion factor from mass to cgs */
  double conv_factor_mass_to_cgs;

  /*! Conversion factor from time to cgs */
  double conv_factor_time_to_cgs;

  /*! Conversion factor from specific energy to cgs */
  double conv_factor_specific_energy_to_cgs;

  /*! Proton mass */
  double proton_mass_cgs_inv;

  /*! Convert Kelvin to internal temperature */
  double T_K_to_int;
};

/**
 * @brief Initialise the black hole properties from the parameter file.
 *
 * For the basic black holes neighbour finding properties we use the
 * defaults from the hydro scheme if the users did not provide specific
 * values.
 *
 * @param bp The #black_holes_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
INLINE static void black_holes_props_init(struct black_holes_props *bp,
                                          const struct phys_const *phys_const,
                                          const struct unit_system *us,
                                          struct swift_params *params,
                                          const struct hydro_props *hydro_props,
                                          const struct cosmology *cosmo) {

  /* Calculate temperature to internal energy conversion factor (all internal
   * units) */
  const double k_B = phys_const->const_boltzmann_k;
  const double m_p = phys_const->const_proton_mass;
  const double mu = hydro_props->mu_ionised;
  bp->temp_to_u_factor = k_B / (mu * hydro_gamma_minus_one * m_p);

  /* Read in the basic neighbour search properties or default to the hydro
     ones if the user did not provide any different values */

  /* Kernel properties */
  bp->eta_neighbours = parser_get_opt_param_float(
      params, "BlackHoles:resolution_eta", hydro_props->eta_neighbours);

  /* Tolerance for the smoothing length Newton-Raphson scheme */
  bp->h_tolerance = parser_get_opt_param_float(params, "BlackHoles:h_tolerance",
                                               hydro_props->h_tolerance);

  /* Get derived properties */
  bp->target_neighbours = pow_dimension(bp->eta_neighbours) * kernel_norm;
  const float delta_eta = bp->eta_neighbours * (1.f + bp->h_tolerance);
  bp->delta_neighbours =
      (pow_dimension(delta_eta) - pow_dimension(bp->eta_neighbours)) *
      kernel_norm;

  /* Number of iterations to converge h */
  bp->max_smoothing_iterations =
      parser_get_opt_param_int(params, "BlackHoles:max_ghost_iterations",
                               hydro_props->max_smoothing_iterations);

  /* Time integration properties */
  const float max_volume_change =
      parser_get_opt_param_float(params, "BlackHoles:max_volume_change", -1);
  if (max_volume_change == -1)
    bp->log_max_h_change = hydro_props->log_max_h_change;
  else
    bp->log_max_h_change = logf(powf(max_volume_change, hydro_dimension_inv));

  /* Initialisation properties  ---------------------------- */

  bp->subgrid_seed_mass =
      parser_get_param_float(params, "YAMAGN:subgrid_seed_mass_Msun");

  bp->seed_n_H_threshold_cgs =
      parser_get_param_float(params, "YAMAGN:seed_n_H_threshold_cgs");

  bp->seed_temperature_threshold_cgs =
      parser_get_param_float(params, "YAMAGN:seed_temperature_threshold_cgs");

  bp->seed_metallicity_threshold =
      parser_get_param_float(params, "YAMAGN:seed_metallicity_threshold");

  bp->seed_during_star_formation =
      parser_get_param_int(params, "YAMAGN:seed_during_star_formation");

  /* Convert to internal units */
  bp->subgrid_seed_mass *= phys_const->const_solar_mass;

  bp->use_subgrid_mass_from_ics =
      parser_get_opt_param_int(params, "YAMAGN:use_subgrid_mass_from_ics", 1);
  if (bp->use_subgrid_mass_from_ics)
    bp->with_subgrid_mass_check =
        parser_get_opt_param_int(params, "YAMAGN:with_subgrid_mass_check", 1);

  /* Accretion parameters ---------------------------------- */

  bp->environment_temperature_cut =
      parser_get_opt_param_float(params, "YAMAGN:environment_temperature_cut", 1.0e5f);

  bp->torque_accretion_norm =
      parser_get_param_float(params, "YAMAGN:torque_accretion_norm");

  bp->dt_accretion_factor =
      parser_get_opt_param_float(params, "YAMAGN:dt_accretion_factor", 1.f);
  if (bp->dt_accretion_factor > 1.f || bp->dt_accretion_factor < 0.f) {
    error("YAMAGN:dt_accretion_factor must be between 0 and 1");
  }

  bp->f_Edd_maximum = parser_get_param_float(params, "YAMAGN:max_eddington_fraction");

  bp->boost_alpha = parser_get_param_float(params, "YAMAGN:boost_alpha");

  bp->with_fixed_T_near_EoS =
      parser_get_param_int(params, "YAMAGN:with_fixed_T_near_EoS");
  if (bp->with_fixed_T_near_EoS) {
    bp->fixed_T_above_EoS_factor =
        exp10(parser_get_param_float(params, "YAMAGN:fixed_T_above_EoS_dex"));
    bp->fixed_u_for_soundspeed =
        parser_get_param_float(params, "YAMAGN:fixed_T_near_EoS_K") /
        units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
    bp->fixed_u_for_soundspeed *= bp->temp_to_u_factor;
  }

  /* Feedback parameters ---------------------------------- */

  char temp[40];
  parser_get_param_string(params, "YAMAGN:AGN_feedback_model", temp);
  if (strcmp(temp, "Random") == 0)
    bp->feedback_model = AGN_random_ngb_model;
  else if (strcmp(temp, "Isotropic") == 0)
    bp->feedback_model = AGN_isotropic_model;
  else if (strcmp(temp, "MinimumDistance") == 0)
    bp->feedback_model = AGN_minimum_distance_model;
  else if (strcmp(temp, "MinimumDensity") == 0)
    bp->feedback_model = AGN_minimum_density_model;
  else
    error(
        "The AGN feedback model must be either 'Random', 'MinimumDistance', "
        "'MinimumDensity' or 'Isotropic', not %s",
        temp);

  bp->AGN_deterministic =
      parser_get_param_int(params, "YAMAGN:AGN_use_deterministic_feedback");

  bp->epsilon_f =
      parser_get_param_float(params, "YAMAGN:coupling_efficiency");

  const double T_K_to_int =
      1. / units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  bp->jet_heating_velocity_threshold = 
      parser_get_param_float(params, "YAMAGN:jet_heating_velocity_threshold");

  /* Convert to internal units */
  const float jet_heating_velocity_threshold = 
      bp->jet_heating_velocity_threshold * 1.0e5f;
  bp->jet_heating_velocity_threshold =
      jet_heating_velocity_threshold / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);

  bp->xray_heating_velocity_threshold =
      parser_get_param_float(params, "YAMAGN:xray_heating_velocity_threshold");

  /* Convert to internal units */
  const float xray_heating_velocity_threshold = 
      bp->xray_heating_velocity_threshold * 1.0e5f;
  bp->xray_heating_velocity_threshold =
      xray_heating_velocity_threshold / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);

  bp->xray_maximum_heating_factor = 
      parser_get_opt_param_float(params, "YAMAGN:xray_maximum_heating_factor",
            1000.0f);

  bp->xray_kinetic_fraction =
      parser_get_opt_param_float(params, "YAMAGN:xray_kinetic_fraction",
            0.5f);

  bp->xray_heating_n_H_threshold_cgs = 
      parser_get_opt_param_float(params, "YAMAGN:xray_heating_n_H_threshold_cgs",
        0.13f);

  bp->xray_heating_T_threshold_cgs =
      parser_get_opt_param_float(params, "YAMAGN:xray_heating_T_threshold_cgs",
        5.0e5f);

  bp->jet_velocity = 
      parser_get_param_float(params, "YAMAGN:jet_velocity");

  /* Convert to internal units */
  const float jet_velocity = bp->jet_velocity * 1.0e5f;
  bp->jet_velocity =
      jet_velocity / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);

  bp->jet_temperature = 
      parser_get_param_float(params, "YAMAGN:jet_temperature");

  bp->eddington_fraction_lower_boundary = 
      parser_get_param_float(params, "YAMAGN:eddington_fraction_lower_boundary");

  bp->jet_mass_min_Msun = 
      parser_get_opt_param_float(params, "YAMAGN:jet_mass_min_Msun", 4.5e7f);

  bp->jet_mass_max_Msun = 
      parser_get_opt_param_float(params, "YAMAGN:jet_mass_max_Msun", 5.0e7f); 

  bp->quasar_wind_momentum_flux =
      parser_get_param_float(params, "YAMAGN:quasar_wind_momentum_flux");

  bp->dt_feedback_factor =
      parser_get_opt_param_float(params, "YAMAGN:dt_feedback_factor", 1.f);
  if (bp->dt_feedback_factor > 1.f || bp->dt_feedback_factor < 0.f) {
    error("YAMAGN:dt_feedback_factor must be between 0 and 1");
  }

  bp->fixed_spin = 
        parser_get_param_float(params, "YAMAGN:fixed_spin");
  if (bp->fixed_spin >= 1.f || bp->fixed_spin <= 0.f) {
    error("Black hole must have spin > 0.0 and < 1.0");
  }

  bp->A_lupi = powf(0.9663f - 0.9292f * bp->fixed_spin, -0.5639f);
  bp->B_lupi = powf(4.627f - 4.445f * bp->fixed_spin, -0.5524f);
  bp->C_lupi = powf(827.3f - 718.1f * bp->fixed_spin, -0.7060f);

  /* Fit to Benson & Babul (2009) disk+wind efficiency */
  const float first = (0.83f / (1.085f - bp->fixed_spin)) - 3.455f;
  const float second = (0.14f / (1.085f - bp->fixed_spin)) - 1.9118f + 
                       (0.003f / (1.085f - powf(bp->fixed_spin, 2.1f))) + 
                       (0.0002f / (1.085f - powf(bp->fixed_spin, 3.f)));
  bp->jet_efficiency = powf(1.f / (powf(1.f / powf(10.f, first), 4.f) + 
                            powf(1.f / powf(10.f, second), 4.f)), 0.25f);

  bp->jet_loading = 2.f * bp->jet_efficiency * 
                    phys_const->const_speed_light_c * phys_const->const_speed_light_c /
                    (bp->jet_velocity * bp->jet_velocity + 
                        3.f * (phys_const->const_boltzmann_k * bp->jet_temperature) / 
                        (0.59 * phys_const->const_proton_mass));

  const float R = 1.f / bp->eddington_fraction_upper_boundary; 
  const float eta_at_slim_disk_boundary = 
        (R / 16.f) * bp->A_lupi * ((0.985f / (R + (5.f / 8.f) * bp->B_lupi)) + 
                                (0.015f / (R + (5.f / 8.f) * bp->C_lupi)));
  /* Divide C_factor by M_dot,Edd later */
  const float C_factor = eta_at_slim_disk_boundary / bp->eddington_fraction_lower_boundary;

  bp->jet_quadratic_term = (1.f + bp->jet_efficiency + bp->jet_loading) / C_factor;
  /* Overwrite the value, we need to keep it continuous over all M_dot,BH/M_dot,Edd */
  bp->epsilon_r = eta_at_slim_disk_boundary;
  if (bp->epsilon_r > 1.f) error("Somehow epsilon_r is greater than 1.0.");

  /* These are for momentum constrained winds */
  bp->quasar_wind_mass_loading = bp->quasar_wind_momentum_flux * bp->epsilon_f * bp->epsilon_r *
          (phys_const->const_speed_light_c / bp->quasar_wind_speed);
  bp->quasar_f_accretion = 1.f / (1.f + bp->quasar_wind_mass_loading);

  bp->adaf_coupling = parser_get_param_float(params, "YAMAGN:adaf_coupling");
  bp->slim_disk_coupling = parser_get_param_float(params, "YAMAGN:slim_disk_coupling");
  bp->quasar_coupling = parser_get_param_float(params, "YAMAGN:quasar_coupling");

  bp->slim_disk_wind_momentum_flux = 
        parser_get_param_float(params, "YAMAGN:slim_disk_wind_momentum_flux");
  bp->adaf_wind_momentum_flux =
        parser_get_param_float(params, "YAMAGN:adaf_wind_momentum_flux");

  bp->slim_disk_wind_speed =
        parser_get_param_float(params, "YAMAGN:slim_disk_wind_speed");
  bp->adaf_f_accretion = 
        parser_get_param_float(params, "YAMAGN:adaf_f_accretion");

  /* Reposition parameters --------------------------------- */

  bp->max_reposition_mass =
      parser_get_param_float(params, "YAMAGN:max_reposition_mass") *
      phys_const->const_solar_mass;
  bp->max_reposition_distance_ratio =
      parser_get_param_float(params, "YAMAGN:max_reposition_distance_ratio");

  bp->with_reposition_velocity_threshold = parser_get_param_int(
      params, "YAMAGN:with_reposition_velocity_threshold");

  if (bp->with_reposition_velocity_threshold) {
    bp->max_reposition_velocity_ratio = parser_get_param_float(
        params, "YAMAGN:max_reposition_velocity_ratio");

    /* Prevent nonsensical input */
    if (bp->max_reposition_velocity_ratio <= 0)
      error("max_reposition_velocity_ratio must be positive, not %f.",
            bp->max_reposition_velocity_ratio);

    bp->min_reposition_velocity_threshold = parser_get_param_float(
        params, "YAMAGN:min_reposition_velocity_threshold");
    /* Convert from km/s to internal units */
    bp->min_reposition_velocity_threshold *=
        (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));
  }

  bp->set_reposition_speed =
      parser_get_param_int(params, "YAMAGN:set_reposition_speed");

  if (bp->set_reposition_speed) {
    bp->reposition_coefficient_upsilon = parser_get_param_float(
        params, "YAMAGN:reposition_coefficient_upsilon");

    /* Prevent the user from making silly wishes */
    if (bp->reposition_coefficient_upsilon <= 0)
      error(
          "reposition_coefficient_upsilon must be positive, not %f "
          "km/s/M_sun.",
          bp->reposition_coefficient_upsilon);

    /* Convert from km/s to internal units */
    bp->reposition_coefficient_upsilon *=
        (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));

    /* Scaling parameters with BH mass and gas density */
    bp->reposition_reference_mass =
        parser_get_param_float(params, "YAMAGN:reposition_reference_mass") *
        phys_const->const_solar_mass;
    bp->reposition_exponent_mass = parser_get_opt_param_float(
        params, "YAMAGN:reposition_exponent_mass", 2.0);
    bp->reposition_reference_n_H =
        parser_get_param_float(params, "YAMAGN:reposition_reference_n_H");
    bp->reposition_exponent_n_H = parser_get_opt_param_float(
        params, "YAMAGN:reposition_exponent_n_H", 1.0);
  }

  bp->correct_bh_potential_for_repositioning =
      parser_get_param_int(params, "YAMAGN:with_potential_correction");

  /* Merger parameters ------------------------------------- */

  bp->minor_merger_threshold =
      parser_get_param_float(params, "YAMAGN:threshold_minor_merger");

  bp->major_merger_threshold =
      parser_get_param_float(params, "YAMAGN:threshold_major_merger");

  char temp2[40];
  parser_get_param_string(params, "YAMAGN:merger_threshold_type", temp2);
  if (strcmp(temp2, "CircularVelocity") == 0)
    bp->merger_threshold_type = BH_mergers_circular_velocity;
  else if (strcmp(temp2, "EscapeVelocity") == 0)
    bp->merger_threshold_type = BH_mergers_escape_velocity;
  else if (strcmp(temp2, "DynamicalEscapeVelocity") == 0)
    bp->merger_threshold_type = BH_mergers_dynamical_escape_velocity;
  else
    error(
        "The BH merger model must be either CircularVelocity, EscapeVelocity, "
        "or DynamicalEscapeVelocity, not %s",
        temp2);

  bp->max_merging_distance_ratio =
      parser_get_param_float(params, "YAMAGN:merger_max_distance_ratio");

  /* ---- Black hole time-step properties ------------------ */

  const double Myr_in_cgs = 1e6 * 365.25 * 24. * 60. * 60.;

  const double time_step_min_Myr = parser_get_opt_param_float(
      params, "YAMAGN:minimum_timestep_Myr", FLT_MAX);

  bp->time_step_min = time_step_min_Myr * Myr_in_cgs /
                      units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Common conversion factors ----------------------------- */

  /* Calculate conversion factor from rho to n_H.
   * Note this assumes primoridal abundance */
  const double X_H = hydro_props->hydrogen_mass_fraction;
  bp->rho_to_n_cgs =
      (X_H / m_p) * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  /* Conversion factor for internal mass to M_solar */
  bp->mass_to_solar_mass = 1.f / phys_const->const_solar_mass;

  bp->kms_to_internal = 1.0e5f / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);

  bp->length_to_parsec = 1.f / phys_const->const_parsec;

  bp->time_to_yr = 1.f / phys_const->const_year;

  bp->time_to_Myr = units_cgs_conversion_factor(us, UNIT_CONV_TIME) /
      (1.e6f * 365.25f * 24.f * 60.f * 60.f);

  /* Some useful conversion values */
  bp->conv_factor_density_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  bp->conv_factor_energy_rate_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY) /
      units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  bp->conv_factor_length_to_cgs = 
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
  bp->conv_factor_mass_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  bp->conv_factor_time_to_cgs = 
      units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  bp->conv_factor_specific_energy_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);

  /* Useful constants */
  bp->proton_mass_cgs_inv =
      1. / (phys_const->const_proton_mass *
            units_cgs_conversion_factor(us, UNIT_CONV_MASS));

  bp->T_K_to_int = T_K_to_int;
}

/**
 * @brief Write a black_holes_props struct to the given FILE as a stream of
 * bytes.
 *
 * @param props the black hole properties struct
 * @param stream the file stream
 */
INLINE static void black_holes_struct_dump(
    const struct black_holes_props *props, FILE *stream) {
  restart_write_blocks((void *)props, sizeof(struct black_holes_props), 1,
                       stream, "black_holes props", "black holes props");
}

/**
 * @brief Restore a black_holes_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param props the black hole properties struct
 * @param stream the file stream
 */
INLINE static void black_holes_struct_restore(
    const struct black_holes_props *props, FILE *stream) {
  restart_read_blocks((void *)props, sizeof(struct black_holes_props), 1,
                      stream, NULL, "black holes props");
}

#endif /* SWIFT_YAM_BLACK_HOLES_PROPERTIES_H */
