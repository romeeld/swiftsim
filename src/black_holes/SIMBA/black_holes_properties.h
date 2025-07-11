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
#ifndef SWIFT_SIMBA_BLACK_HOLES_PROPERTIES_H
#define SWIFT_SIMBA_BLACK_HOLES_PROPERTIES_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "chemistry.h"
#include "exp10.h"
#include "hydro_properties.h"

/* Includes. */
#include <string.h>

enum BH_merger_thresholds {
  BH_mergers_circular_velocity,        /*< v_circ at separation, as in EAGLE */
  BH_mergers_escape_velocity,          /*< v_esc at separation */
  BH_mergers_dynamical_escape_velocity /*< combined v_esc_dyn at separation */
};

enum BH_suppress_models {
  BH_suppress_None,
  BH_suppress_Hopkins22,
  BH_suppress_ExponentialOnTorque,
  BH_suppress_ExponentialOnTotal,
  BH_suppress_OutflowsOnAllGas,
  BH_suppress_OutflowsOnSFGas,
  BH_suppress_ExponentialOnSFGas,
  BH_suppress_OutflowsOrExponential,
  BH_suppress_ExpOutflowsOnTotal
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

  /*! Maximum smoothing length */
  float h_max;

  /*! Minimum smoothing length */
  float h_min;

  /*! Maximal number of iterations to converge h */
  int max_smoothing_iterations;

  /*! Maximal change of h over one time-step */
  float log_max_h_change;

  /* ----- Initialisation properties  ------ */

  /*! Mass of a BH seed at creation time */
  float subgrid_seed_mass;

  /*! Should we use the subgrid mass specified in ICs? */
  int use_subgrid_mass_from_ics;

  /*! Should we enforce positive subgrid masses initially? */
  int with_subgrid_mass_check;

  /* ----- Properties of the accretion model ------ */

  /*! Calculate Bondi accretion rate for individual neighbours? */
  int use_multi_phase_bondi;

  /*! Are we using the subgrid gas properties in the Bondi model? */
  int use_subgrid_bondi;

  /*! Normalisation of the viscuous angular momentum accretion reduction */
  float alpha_visc;

  /*! Radiative efficiency of the black holes. */
  float epsilon_r;

  /*! Maximal fraction of the Eddington rate allowed. */
  float f_Edd;

  /*! Eddington fraction threshold for recording */
  float f_Edd_recording;

  /*! Switch for the Booth & Schaye 2009 model */
  int with_boost_factor;

  /*! value of the boost of the Booth & Schaye 2009 model */
  float boost_alpha;

  /*! Switch for nibbling mode */
  int use_nibbling;

  /*! Minimum gas particle mass in nibbling mode */
  float min_gas_mass_for_nibbling;

  /*! Factor above EoS below which fixed T applies for sound speed */
  float fixed_T_above_EoS_factor;

  /*! Where do we distinguish between hot gas for Bondi? */
  float environment_temperature_cut;

  /*! How much of Mdot,inflow should we accrete? The rest is an outflowing wind
   */
  float f_accretion;

  /*! Number of dynamical times over which gas is accreted from accretion disk */
  float accretion_disk_timescale;

  /*! Method to compute torque accretion rate: 0=Mgas/tdyn on cold gas; 
   * 1=Mgas/tdyn on corot gas; 2=Simba-style(HQ11) */
  int torque_accretion_method;

  /*! Normalization of the torque accretion rate */
  float torque_accretion_norm;

  /*! Factor in front of M/(dM/dt) for timestepping */
  float dt_accretion_factor;

  /*! Suppress growth via Hopkins+'21 */
  int suppress_growth;

  /*! sigma_crit from Hopkins+'21 */
  float sigma_crit_Msun_pc2;

  /*! A factor to account for not being able to resolve sigma_crit */
  float sigma_crit_resolution_factor;

  /*! The factor for exponentially limiting black hole growth in early stages.*/
  float characteristic_suppression_mass;

  /*! Bondi rate cannot exceed the rate for BH of this mass. */
  float bondi_rate_limiting_bh_mass;

  /* ---- Properties of the feedback model ------- */

  /*! Maximum heating temperature */
  float jet_maximum_temperature;

  /*! When does the jet start heating? (internal velocity) */
  float jet_heating_velocity_threshold;

  /*! At what v_kick does the X-ray heating start? (internal velocity) */
  float xray_heating_velocity_threshold;

  /*! How many times the particle's u value can we heat? */
  float xray_maximum_heating_factor;

  /*! How much of the X-ray energy should go into velocity for dense gas? */
  float xray_kinetic_fraction;

  /*! Above this density we should split X-ray energy into kinetic */
  float xray_heating_n_H_threshold_cgs;

  /*! Below this temperature we should split X-ray energy into kinetic */
  float xray_heating_T_threshold_cgs;

  /*! Amount of radiation lost to X-rays */
  float xray_radiation_loss;

  /*! The limit for X-ray feedback, below this use X-ray feedback */
  float xray_f_gas_limit;

  /*! Should cooling be shut off for X-ray impacted particles? */
  int xray_shutoff_cooling;

  /*! Physical max. velocity of the jet (internal velocity), for MBH=1e8 */
  float jet_velocity;

  /*! factor setting maximum jet velcoity as multiplier of jet_velocity */
  float jet_velocity_max_multiplier;

  /*! Power law scaling of v_jet with MBH; set =0 for constant value at jet_velocity */
  float jet_velocity_scaling_with_mass;

  /*! The temperature of the jet. Set < 0.f for halo virial temperature */
  float jet_temperature;

  /*! Should we scale the jet temperature with the BH mass? */
  int scale_jet_temperature;

  /*! What lower Mdot,BH/Mdot,Edd boundary does the jet activate? */
  float eddington_fraction_lower_boundary;

  /*! Full jets are always on if Bondi accretion fraction above this */
  float bondi_fraction_thresh_always_jet;

  /*! Always use maximum jet speed above this black hole mass */
  float mass_thresh_always_jet;

  /*! Always use maximum jet speed above this black hole accretion rate */
  float accr_rate_thresh_always_jet;

  /*! Add spread around the jet velocity alpha + beta * random */
  float jet_velocity_spread_alpha;

  /*! Add spread around the jet velocity alpha + beta * random */
  float jet_velocity_spread_beta;

  /*! Minimum mass for starting the jet (Msun) */
  float jet_mass_min_Msun;

  /*! The maximum mass for the jet, to be added to jet_mass_min_Msun (Msun) */
  float jet_mass_spread_Msun;

  /*! Constrains momentum of outflowing wind to p = F * L / c */
  float wind_momentum_flux;

  /*! This times 1/H(z) for decoupling time */
  float wind_decouple_time_factor;

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
  double temp_to_u_factor;

  /*! Conversion factor from physical density to n_H [cgs] */
  double rho_to_n_cgs;

  /*! Conversion factor from internal mass to solar masses */
  double mass_to_solar_mass;

  /*! Conversion factor from km/s to internal velocity units (without a-factor)
   */
  double kms_to_internal;

  /*! Conversion factor from internal length to parsec */
  double length_to_parsec;

  /*! Conversion factor from internal time to yr */
  double time_to_yr;

  /*! Conversion factor from internal time to Myr */
  double time_to_Myr;

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

  /* ------------ Stellar feedback properties for eta computation --------------- */

  /*! Normalization for the mass loading curve */
  float FIRE_eta_normalization;

  /*! The location (in internal mass units) where the break in the
   * mass loading curve occurs */
  float FIRE_eta_break;

  /*! The power-law slope of eta below FIRE_eta_break */
  float FIRE_eta_lower_slope;

  /*! The power-law slope of eta above FIRE_eta_break */
  float FIRE_eta_upper_slope;

  /*! The minimum galaxy stellar mass in internal units */
  float minimum_galaxy_stellar_mass;

  /*! Star formation efficiency */
  float star_formation_efficiency;
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

  /* Common conversion factors ----------------------------- */

  /* Calculate temperature to internal energy conversion factor (all internal
   * units) */
  const double k_B = phys_const->const_boltzmann_k;
  const double m_p = phys_const->const_proton_mass;
  const double mu = hydro_props->mu_ionised;
  bp->temp_to_u_factor = k_B / (mu * hydro_gamma_minus_one * m_p);

  /* Calculate conversion factor from rho to n_H.
   * Note this assumes primoridal abundance */
  const double X_H = hydro_props->hydrogen_mass_fraction;
  bp->rho_to_n_cgs =
      (X_H / m_p) * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  bp->kms_to_internal =
      1.0e5 / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);

  bp->length_to_parsec = 1. / phys_const->const_parsec;

  bp->time_to_yr = 1. / phys_const->const_year;

  bp->time_to_Myr = bp->time_to_yr * 1.e-6;

  if (engine_rank == 0) {
    message("TIME_TO_YR CONVERSION: %g %g %g", bp->time_to_yr, bp->time_to_Myr, 
	          units_cgs_conversion_factor(us, UNIT_CONV_TIME) / 
                (1.e6 * 365.25 * 24. * 60. * 60.) );
  }

  /* Some useful conversion values */
  bp->conv_factor_density_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  bp->conv_factor_energy_rate_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY) /
      units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  bp->conv_factor_length_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
  bp->conv_factor_mass_to_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  bp->conv_factor_time_to_cgs = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  bp->conv_factor_specific_energy_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  /* Conversion factor for internal mass to M_solar */
  bp->mass_to_solar_mass = 1. / phys_const->const_solar_mass;

  /* Useful constants */
  bp->proton_mass_cgs_inv =
      1. / (phys_const->const_proton_mass *
            units_cgs_conversion_factor(us, UNIT_CONV_MASS));

  const double T_K_to_int =
      1. / units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  bp->T_K_to_int = T_K_to_int;

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

  /* Maximal smoothing length */
  bp->h_max = parser_get_opt_param_float(params, "BlackHoles:h_max",
                                         hydro_props->h_max);

  /* Maximal smoothing length */
  bp->h_min = parser_get_opt_param_float(params, "BlackHoles:h_min",
                                         hydro_props->h_min);

  /* Time integration properties */
  const float max_volume_change =
      parser_get_opt_param_float(params, "BlackHoles:max_volume_change", -1);
  if (max_volume_change == -1)
    bp->log_max_h_change = hydro_props->log_max_h_change;
  else
    bp->log_max_h_change = logf(powf(max_volume_change, hydro_dimension_inv));

  /* Initialisation properties  ---------------------------- */

  bp->subgrid_seed_mass =
      parser_get_param_float(params, "SIMBAAGN:subgrid_seed_mass_Msun");

  /* Convert to internal units */
  bp->subgrid_seed_mass *= phys_const->const_solar_mass;

  bp->use_subgrid_mass_from_ics =
      parser_get_opt_param_int(params, "SIMBAAGN:use_subgrid_mass_from_ics", 1);
  if (bp->use_subgrid_mass_from_ics)
    bp->with_subgrid_mass_check =
        parser_get_opt_param_int(params, "SIMBAAGN:with_subgrid_mass_check", 1);

  /* Accretion parameters ---------------------------------- */

  bp->environment_temperature_cut = parser_get_opt_param_float(
      params, "SIMBAAGN:environment_temperature_cut", 1.0e5f);

  bp->f_accretion = parser_get_param_float(params, "SIMBAAGN:f_accretion");

  bp->accretion_disk_timescale = parser_get_opt_param_float(
      params, "SIMBAAGN:accretion_disk_timescale", 0.f);
  bp->accretion_disk_timescale /= bp->time_to_Myr; /* to internal units */

  bp->torque_accretion_method =
      parser_get_opt_param_int(params, "SIMBAAGN:torque_accretion_method", 0);

  bp->torque_accretion_norm =
      parser_get_param_float(params, "SIMBAAGN:torque_accretion_norm");

  bp->dt_accretion_factor =
      parser_get_opt_param_float(params, "SIMBAAGN:dt_accretion_factor", 1.f);
  if (bp->dt_accretion_factor > 1.f || bp->dt_accretion_factor < 0.f) {
    error("SIMBAAGN:dt_accretion_factor must be between 0 and 1");
  }

  bp->suppress_growth =
      parser_get_param_int(params, "SIMBAAGN:suppress_growth");

  if (bp->suppress_growth == BH_suppress_OutflowsOnAllGas || 
          bp->suppress_growth == BH_suppress_OutflowsOnSFGas || 
          BH_suppress_ExpOutflowsOnTotal > 0) {
    bp->FIRE_eta_normalization =
        parser_get_param_float(params, "KIARAFeedback:FIRE_eta_normalization");
    bp->FIRE_eta_break =
        parser_get_param_float(params, "KIARAFeedback:FIRE_eta_break_Msun");
    bp->FIRE_eta_break /= bp->mass_to_solar_mass;
    bp->FIRE_eta_lower_slope =
        parser_get_param_float(params, "KIARAFeedback:FIRE_eta_lower_slope");
    bp->FIRE_eta_upper_slope =
        parser_get_param_float(params, "KIARAFeedback:FIRE_eta_upper_slope");
    bp->minimum_galaxy_stellar_mass =
      parser_get_param_float(params, 
          "KIARAFeedback:minimum_galaxy_stellar_mass_Msun");
    bp->minimum_galaxy_stellar_mass /= bp->mass_to_solar_mass;
     /* Get the star formation efficiency */
    bp->star_formation_efficiency = parser_get_param_float(
        params, "KIARAStarFormation:star_formation_efficiency");
  }

  bp->sigma_crit_Msun_pc2 =
      parser_get_param_float(params, "SIMBAAGN:sigma_crit_Msun_pc2");

  bp->sigma_crit_resolution_factor =
      parser_get_param_float(params, "SIMBAAGN:sigma_crit_resolution_factor");

  bp->characteristic_suppression_mass = parser_get_param_float(
      params, "SIMBAAGN:bh_characteristic_suppression_mass");

  bp->bondi_rate_limiting_bh_mass =
      parser_get_param_float(params, "SIMBAAGN:bondi_rate_limiting_bh_mass");
  bp->bondi_rate_limiting_bh_mass /=
      bp->mass_to_solar_mass; /* Convert to internal units */

  bp->use_multi_phase_bondi =
      parser_get_param_int(params, "SIMBAAGN:use_multi_phase_bondi");

  bp->use_subgrid_bondi =
      parser_get_param_int(params, "SIMBAAGN:use_subgrid_bondi");

  if (bp->use_multi_phase_bondi && bp->use_subgrid_bondi)
    error(
        "Cannot run with both the multi-phase Bondi and subgrid Bondi models "
        "at the same time!");

  /*  Booth & Schaye (2009) Parameters */
  bp->with_boost_factor =
      parser_get_param_int(params, "SIMBAAGN:with_boost_factor");

  if (bp->with_boost_factor) {
    bp->boost_alpha = parser_get_param_float(params, "SIMBAAGN:boost_alpha");
  }

  bp->fixed_T_above_EoS_factor =
      exp10(parser_get_param_float(params, "SIMBAAGN:fixed_T_above_EoS_dex"));

  bp->epsilon_r =
      parser_get_param_float(params, "SIMBAAGN:radiative_efficiency");
  if (bp->epsilon_r > 1.f)
    error("SIMBAAGN:radiative_efficiency must be <= 1, not %f.", bp->epsilon_r);

  bp->use_nibbling =
      parser_get_opt_param_int(params, "SIMBAAGN:use_nibbling", 1);
  if (bp->use_nibbling) {
    bp->min_gas_mass_for_nibbling = parser_get_param_float(
        params, "SIMBAAGN:min_gas_mass_for_nibbling_Msun");
    bp->min_gas_mass_for_nibbling *= phys_const->const_solar_mass;

    if ((bp->min_gas_mass_for_nibbling < 1e-5 * bp->subgrid_seed_mass) ||
        (bp->min_gas_mass_for_nibbling > 1e5 * bp->subgrid_seed_mass)) {
      error(
          "The BH seeding mass and minimal gas mass for nibbling differ by "
          "more "
          "than 10^5. That is probably indicating a typo in the parameter "
          "file.");
    }
  } else {
    error("It is impossible to use SIMBA without nibbling.");
  }

  bp->f_Edd = parser_get_param_float(params, "SIMBAAGN:max_eddington_fraction");

  bp->bondi_fraction_thresh_always_jet = 
      parser_get_opt_param_float(params, 
          "SIMBAAGN:bondi_fraction_thresh_always_jet", FLT_MAX);

  bp->mass_thresh_always_jet = 
      parser_get_opt_param_float(params, 
          "SIMBAAGN:mass_thresh_always_jet", FLT_MAX);

  bp->accr_rate_thresh_always_jet = 
      parser_get_opt_param_float(params, 
          "SIMBAAGN:accr_rate_thresh_always_jet", FLT_MAX);
  /* Convert from Mo/yr to internal units */
  bp->accr_rate_thresh_always_jet *= 
      bp->time_to_yr / bp->mass_to_solar_mass;

  bp->jet_velocity_spread_alpha =
      parser_get_opt_param_float(params, 
          "SIMBAAGN:jet_velocity_spread_alpha", 0.8f);

  bp->jet_velocity_spread_beta =
      parser_get_opt_param_float(params, 
          "SIMBAAGN:jet_velocity_spread_beta", 0.4f);

  /* Feedback parameters ---------------------------------- */

  bp->jet_maximum_temperature =
      parser_get_opt_param_float(params, "SIMBAAGN:jet_maximum_temperature", 
                                 1.e8f);

  bp->jet_heating_velocity_threshold =
      parser_get_param_float(params, "SIMBAAGN:jet_heating_velocity_threshold");
  bp->jet_heating_velocity_threshold *= bp->kms_to_internal;

  bp->xray_heating_velocity_threshold = parser_get_param_float(
      params, "SIMBAAGN:xray_heating_velocity_threshold");
  bp->xray_heating_velocity_threshold *= bp->kms_to_internal;

  bp->xray_maximum_heating_factor = parser_get_opt_param_float(
      params, "SIMBAAGN:xray_maximum_heating_factor", 1000.0f);

  bp->xray_kinetic_fraction = parser_get_opt_param_float(
      params, "SIMBAAGN:xray_kinetic_fraction", 0.5f);
  if (bp->xray_kinetic_fraction < 0.f || bp->xray_kinetic_fraction > 1.f) {
    error("xray_kinetic_fraction must be between 0 <= f_xray <= 1.");
  }
  
  bp->xray_heating_n_H_threshold_cgs = parser_get_opt_param_float(
      params, "SIMBAAGN:xray_heating_n_H_threshold_cgs", 0.13f);

  bp->xray_heating_T_threshold_cgs = parser_get_opt_param_float(
      params, "SIMBAAGN:xray_heating_T_threshold_cgs", 5.0e5f);

  bp->xray_radiation_loss =
      parser_get_param_float(params, "SIMBAAGN:xray_radiation_loss");

  bp->xray_f_gas_limit =
      parser_get_param_float(params, "SIMBAAGN:xray_f_gas_limit");

  bp->xray_shutoff_cooling =
      parser_get_param_int(params, "SIMBAAGN:xray_shutoff_cooling");

  bp->jet_velocity = parser_get_param_float(params, "SIMBAAGN:jet_velocity");
  bp->jet_velocity *= bp->kms_to_internal;

  bp->jet_velocity_scaling_with_mass = 
      parser_get_param_float(params, "SIMBAAGN:jet_velocity_scaling_with_mass");

  bp->jet_velocity_max_multiplier = 
      parser_get_param_float(params, "SIMBAAGN:jet_velocity_max_multiplier");

  bp->scale_jet_temperature =
      parser_get_param_int(params, "SIMBAAGN:scale_jet_temperature_with_mass");

  bp->jet_temperature =
      parser_get_param_float(params, "SIMBAAGN:jet_temperature");
  if ((bp->scale_jet_temperature <= 1 && bp->jet_temperature < 100.) || 
          (bp->scale_jet_temperature == 2 && bp->jet_temperature > 100.)) {
    error("Parameter SIMBAAGN:jet_temperature=%g doesn't look right!  "
          " for scale_jet_temperature==1 it should be a temperature like 1e7, "
          "whereas if scale_jet_temperature==2 it should be a multiplier for "
          "Tvir like ~1-10 or so", 
          bp->jet_temperature);
  }
      
  bp->eddington_fraction_lower_boundary = parser_get_param_float(
      params, "SIMBAAGN:eddington_fraction_lower_boundary");

  bp->jet_mass_min_Msun =
      parser_get_opt_param_float(params, "SIMBAAGN:jet_mass_min_Msun", 4.5e7f);

  bp->jet_mass_spread_Msun =
      parser_get_opt_param_float(params, "SIMBAAGN:jet_mass_spread_Msun", 5.0e6f);

  bp->wind_momentum_flux =
      parser_get_param_float(params, "SIMBAAGN:wind_momentum_flux");

  bp->wind_decouple_time_factor =
      parser_get_param_float(params, "SIMBAAGN:wind_decouple_time_factor");

  /* Reposition parameters --------------------------------- */

  bp->max_reposition_mass =
      parser_get_param_float(params, "SIMBAAGN:max_reposition_mass") *
      phys_const->const_solar_mass;
  bp->max_reposition_distance_ratio =
      parser_get_param_float(params, "SIMBAAGN:max_reposition_distance_ratio");

  bp->with_reposition_velocity_threshold = parser_get_param_int(
      params, "SIMBAAGN:with_reposition_velocity_threshold");

  if (bp->with_reposition_velocity_threshold) {
    bp->max_reposition_velocity_ratio = parser_get_param_float(
        params, "SIMBAAGN:max_reposition_velocity_ratio");

    /* Prevent nonsensical input */
    if (bp->max_reposition_velocity_ratio <= 0)
      error("max_reposition_velocity_ratio must be positive, not %f.",
            bp->max_reposition_velocity_ratio);

    bp->min_reposition_velocity_threshold = parser_get_param_float(
        params, "SIMBAAGN:min_reposition_velocity_threshold");
    /* Convert from km/s to internal units */
    bp->min_reposition_velocity_threshold *=
        (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));
  }

  bp->set_reposition_speed =
      parser_get_param_int(params, "SIMBAAGN:set_reposition_speed");

  if (bp->set_reposition_speed) {
    bp->reposition_coefficient_upsilon = parser_get_param_float(
        params, "SIMBAAGN:reposition_coefficient_upsilon");

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
        parser_get_param_float(params, "SIMBAAGN:reposition_reference_mass") *
        phys_const->const_solar_mass;
    bp->reposition_exponent_mass = parser_get_opt_param_float(
        params, "SIMBAAGN:reposition_exponent_mass", 2.0);
    bp->reposition_reference_n_H =
        parser_get_param_float(params, "SIMBAAGN:reposition_reference_n_H");
    bp->reposition_exponent_n_H = parser_get_opt_param_float(
        params, "SIMBAAGN:reposition_exponent_n_H", 1.0);
  }

  bp->correct_bh_potential_for_repositioning =
      parser_get_param_int(params, "SIMBAAGN:with_potential_correction");

  /* Merger parameters ------------------------------------- */

  bp->minor_merger_threshold =
      parser_get_param_float(params, "SIMBAAGN:threshold_minor_merger");

  bp->major_merger_threshold =
      parser_get_param_float(params, "SIMBAAGN:threshold_major_merger");

  char temp2[40];
  parser_get_param_string(params, "SIMBAAGN:merger_threshold_type", temp2);
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
      parser_get_param_float(params, "SIMBAAGN:merger_max_distance_ratio");

  /* ---- Black hole time-step properties ------------------ */

  const double Myr_in_cgs = 1e6 * 365.25 * 24. * 60. * 60.;

  const double time_step_min_Myr = parser_get_opt_param_float(
      params, "SIMBAAGN:minimum_timestep_Myr", FLT_MAX);

  bp->time_step_min = time_step_min_Myr * Myr_in_cgs /
                      units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  if (engine_rank == 0) {
    message("Black holes kernel: %s with eta=%f (%.2f neighbours).", 
            kernel_name, bp->eta_neighbours, bp->target_neighbours);

    message("Black holes relative tolerance in h: %.5f (+/- %.4f neighbours).",
            bp->h_tolerance, bp->delta_neighbours);

    message("Black hole model is SIMBA");
    message("Black hole h_max is %g", bp->h_max);
    message("Black hole h_min is %g", bp->h_min);
  }
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

#endif /* SWIFT_SIMBA_BLACK_HOLES_PROPERTIES_H */
