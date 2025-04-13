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
#ifndef SWIFT_OBSIDIAN_BLACK_HOLES_PROPERTIES_H
#define SWIFT_OBSIDIAN_BLACK_HOLES_PROPERTIES_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "chemistry.h"
#include "exp10.h"
#include "hydro_properties.h"

/* Includes. */
#include <string.h>


enum BH_states {
  BH_states_adaf = 0, /* < 0.03 Mdot,BH / Mdot,Edd */
  BH_states_quasar,   /* 0.03 < Mdot,BH / Mdot,Edd < 0.3 */
  BH_states_slim_disk /* Mdot,BH/Mdot,Edd > 0.3 */
};

enum BH_merger_thresholds {
  BH_mergers_circular_velocity,        /*< v_circ at separation, as in EAGLE */
  BH_mergers_escape_velocity,          /*< v_esc at separation */
  BH_mergers_dynamical_escape_velocity /*< combined v_esc_dyn at separation */
};

enum BH_loading_types {
  BH_jet_momentum_loaded,  /* Momentum loaded jet, with subgrid energy loading */
  BH_jet_energy_loaded     /* Energy loaded jet, no subgrid */
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

  /*! Maximum smoothing length */
  float h_max;

  /*! Minimum smoothing length */
  float h_min;

  /*! Tolerance on neighbour number  (for info only)*/
  float delta_neighbours;

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

  /*! Radiative efficiency of the black holes. */
  float epsilon_r;

  /*! Maximal fraction of the Eddington rate allowed for total accretion. */
  float f_Edd_maximum;

  /*! Maximal fraction of the Eddington rate allowed for Bondi accretion alone. */
  float f_Edd_Bondi_maximum;

  /*! Minimum gas particle mass in nibbling mode */
  float min_gas_mass_for_nibbling;

  /*! Switch to calculate the sound speed with a fixed T near the EoS */
  int with_fixed_T_near_EoS;

  /*! Factor above EoS below which fixed T applies for sound speed */
  float fixed_T_above_EoS_factor;

  /*! Where do we distinguish between hot gas for Bondi? */
  float environment_temperature_cut;

  /*! Number of dynamical times over which gas is accreted from accretion disk */
  float bh_accr_dyn_time_fac;

  /*! Method to compute torque accretion rate: 0=Mgas/tdyn on all gas; 1=Mgas/tdyn on disk gas; 2=Simba-style(HQ11) */
  int torque_accretion_method;

  /*! Normalization of the torque accretion rate */
  float torque_accretion_norm;

  /*! Factor in front of M/(dM/dt) for timestepping */
  float dt_accretion_factor;

  /*! Factor for exponentially limiting black hole growth in early stages. */
  float bh_characteristic_suppression_mass;

  /*! Method to suppress early growth of BH */
  int suppress_growth;

  /*! A from Lupi+17 */
  float A_lupi;

  /*! B from Lupi+17 */
  float B_lupi;

  /*! C from Lupi+17 */
  float C_lupi;

  /*! The spin of EVERY black hole */
  float fixed_spin;

  /*! Method to compute the dynamical time within the kernel */
  int dynamical_time_calculation_method;

  /* ---- Properties of the feedback model ------- */

  /*! The loading for the jet: momentum or energy */
  enum BH_loading_types jet_loading_type;

  /*! What is the physical max. velocity of the jet? (km/s) */
  float jet_velocity;

  /*! The temperature of the jet. Set < 0.f for halo virial temperature */
  float jet_temperature;

  /*! What lower Mdot,BH/Mdot,Edd boundary does the jet activate? */
  float eddington_fraction_lower_boundary;

  /*! What upper Mdot,BH/Mdot,Edd boundary does the slim disk mode activate? */
  float eddington_fraction_upper_boundary;

  /*! How long to decouple black hole winds? */
  float wind_decouple_time_factor;

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

  /*! The disk wind efficiency from Benson & Babul 2009 */
  float adaf_disk_efficiency;

  /*! The wind speed of the ADAF outflow */
  float adaf_wind_speed;

  /*! eps_f for the ADAF mode */
  float adaf_coupling;

  /*! Power-law scaling of adaf_coupling with (1+z) */
  float adaf_z_scaling;

  /*! f_acc for the ADAF mode */
  float adaf_f_accretion;

  /*! The maximum temperature to heat a gas particle */
  float adaf_maximum_temperature;

  /*! Above this density we should shut off cooling for heated particles */
  double adaf_heating_n_H_threshold_cgs;

  /*! Below this temperature we should shut off cooling for particles */
  double adaf_heating_T_threshold_cgs;

  /*! Lower mass limit (internal units) for BH to enter ADAF mode */
  float adaf_mass_limit;

  /*! A multiplicative factor for delaying cooling on a particle */
  float adaf_cooling_shutoff_factor;

  /*! A multiplicative factor 0. < f < 1. to multiply E_inject in the ADAF mode */
  float adaf_kick_factor;

  /*! Should we use nibbling */
  int use_nibbling;

  /*! Use all of the gas in the kernel to compute Bondi */
  int bondi_use_all_gas;

  /*! Multiplicative factor in front of Bondi rate */
  float bondi_alpha;

  /*! Minimum BH mass for unresolved feedback (internal units) */
  float minimum_black_hole_mass_unresolved;

  /*! The phi term for the slim disk mode (Eq. 9 from Rennehan+24) */
  float slim_disk_phi;

  /*! eps_f for the slim disk mode */
  float slim_disk_coupling;

  /*! wind speed in the slim disk mode */
  float slim_disk_wind_speed;

  /*! Is the slim disk jet model active? */
  int slim_disk_jet_active;

  /*! The efficiency of the jet */
  float jet_efficiency;

  /*! The mass loading in the jet */
  float jet_mass_loading;
  
  /*! The subgrid jet speed to set the accretion fraction */
  float jet_subgrid_velocity;

  /*! Is the jet isotropic (=1) or along angular momentum (=0) */
  int jet_is_isotropic;

  /*! The minimum mass required before the jet will launch */
  float jet_minimum_reservoir_mass;

  /*! Direction flag to kick the winds by default */
  int default_dir_flag;

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
  double rho_to_n_cgs;

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

  /* Maximum smoothing length */
  bp->h_max = parser_get_opt_param_float(params, "BlackHoles:h_max",
                                         hydro_props->h_max);

  /* Minimum smoothing length */
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
      parser_get_param_float(params, "ObsidianAGN:subgrid_seed_mass_Msun");
      
  /* Convert to internal units */
  bp->subgrid_seed_mass *= phys_const->const_solar_mass;

  bp->use_subgrid_mass_from_ics = parser_get_opt_param_int(params, 
                                  "ObsidianAGN:use_subgrid_mass_from_ics", 1);
  if (bp->use_subgrid_mass_from_ics) {
    bp->with_subgrid_mass_check = parser_get_opt_param_int(params, 
                                  "ObsidianAGN:with_subgrid_mass_check", 1);
  }

  /* Accretion parameters ---------------------------------- */

 /* Conversion factor for internal mass to M_solar */
  bp->mass_to_solar_mass = 1.f / phys_const->const_solar_mass;

  bp->min_gas_mass_for_nibbling =
      parser_get_param_float(params, 
                             "ObsidianAGN:min_gas_mass_for_nibbling_Msun");
  bp->min_gas_mass_for_nibbling /= bp->mass_to_solar_mass;

  const double T_K_to_int =
      1. / units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  bp->environment_temperature_cut =
      parser_get_opt_param_float(params, 
                                 "ObsidianAGN:environment_temperature_cut_K", 
                                 1.0e5f);
  bp->environment_temperature_cut *= T_K_to_int;

  bp->bh_accr_dyn_time_fac = parser_get_opt_param_float(
      params, "ObsidianAGN:bh_accr_dyn_time_fac", 0.f);
  
  bp->torque_accretion_method =
      parser_get_opt_param_int(params, "ObsidianAGN:torque_accretion_method", 0);

  bp->torque_accretion_norm =
      parser_get_param_float(params, "ObsidianAGN:torque_accretion_norm");

  bp->suppress_growth =
      parser_get_opt_param_int(params, "ObsidianAGN:suppress_growth", 0);


  if (bp->suppress_growth == 5 && bp->torque_accretion_method == 2) {
    error("SF-based suppression of BH mass will not work correctly with "
          "Simba-style torque-limited BH growth -- use tdyn method");
  }

  bp->dt_accretion_factor =
      parser_get_opt_param_float(params, "ObsidianAGN:dt_accretion_factor", 1.f);
  if (bp->dt_accretion_factor > 1.f || bp->dt_accretion_factor < 0.f) {
    error("ObsidianAGN:dt_accretion_factor must be between 0 and 1");
  }

  bp->bh_characteristic_suppression_mass = parser_get_opt_param_float(
      params, "ObsidianAGN:bh_characteristic_suppression_mass", 0.f);

  bp->f_Edd_maximum = 
        parser_get_param_float(params, "ObsidianAGN:max_eddington_fraction");

  bp->f_Edd_Bondi_maximum = parser_get_opt_param_float(params, 
                            "ObsidianAGN:max_bondi_eddington_fraction", 1.f);

  bp->fixed_T_above_EoS_factor =
      exp10(parser_get_param_float(params, "ObsidianAGN:fixed_T_above_EoS_dex"));

  bp->dynamical_time_calculation_method = parser_get_opt_param_int(params, 
      "ObsidianAGN:dynamical_time_calculation_method", 1);

  /* Feedback parameters ---------------------------------- */

  bp->kms_to_internal = 1.0e5f / 
        units_cgs_conversion_factor(us, UNIT_CONV_SPEED);

  char temp3[40];
  parser_get_param_string(params, "ObsidianAGN:jet_loading_type", temp3);
  if (strcmp(temp3, "EnergyLoaded") == 0) {
    bp->jet_loading_type = BH_jet_energy_loaded;
  }
  else if (strcmp(temp3, "MomentumLoaded") == 0) {
    bp->jet_loading_type = BH_jet_momentum_loaded;
  }
  else {
    error(
        "The BH jet loading must be either EnergyLoaded or "
        "MomentumLoaded, not %s",
        temp3);
  }

  bp->jet_velocity = 
      parser_get_param_float(params, "ObsidianAGN:jet_velocity_km_s");
  bp->jet_velocity *= bp->kms_to_internal;

  bp->jet_temperature = 
      parser_get_param_float(params, "ObsidianAGN:jet_temperature_K");
  bp->jet_temperature *= T_K_to_int;

  bp->eddington_fraction_lower_boundary = 
      parser_get_param_float(params, 
            "ObsidianAGN:eddington_fraction_lower_boundary");

  bp->eddington_fraction_upper_boundary =
      parser_get_param_float(params, 
            "ObsidianAGN:eddington_fraction_upper_boundary");

  bp->wind_decouple_time_factor =
      parser_get_param_float(params, 
            "ObsidianAGN:wind_decouple_time_factor");

  bp->fixed_spin = 
        parser_get_param_float(params, "ObsidianAGN:fixed_spin");
  if (bp->fixed_spin >= 1.f || bp->fixed_spin <= 0.f) {
    error("Black hole must have spin > 0.0 and < 1.0");
  }

  bp->A_lupi = powf(0.9663f - 0.9292f * bp->fixed_spin, -0.5639f);
  bp->B_lupi = powf(4.627f - 4.445f * bp->fixed_spin, -0.5524f);
  bp->C_lupi = powf(827.3f - 718.1f * bp->fixed_spin, -0.7060f);

  const float phi_bh = -20.2f * powf(bp->fixed_spin, 3.f)
                       -14.9f * powf(bp->fixed_spin, 2.f)
                       + 34.f * bp->fixed_spin
                       + 52.6f;
  const float big_J = bp->fixed_spin / 
        (2.f * (1.f + sqrtf(1.f - pow(bp->fixed_spin, 2.f))));
  const float f_j = powf(big_J, 2.f) + 
        1.38f * powf(big_J, 4.f) - 9.2f * powf(big_J, 6.f);

  bp->jet_efficiency = (1.f / (24.f * M_PI * M_PI)) * powf(phi_bh, 2.f) * f_j;

  if (bp->jet_loading_type == BH_jet_momentum_loaded) {
    bp->jet_mass_loading = 
        bp->jet_efficiency * (phys_const->const_speed_light_c / bp->jet_velocity);
  } else {
    bp->jet_mass_loading =
        2.f * bp->jet_efficiency * powf(
          phys_const->const_speed_light_c / bp->jet_velocity,
          2.f
        );
  }

  bp->jet_subgrid_velocity = 
      parser_get_param_float(params, "ObsidianAGN:jet_subgrid_velocity_km_s");
  bp->jet_subgrid_velocity *= bp->kms_to_internal;

  const float R = 1.f / bp->eddington_fraction_upper_boundary; 
  const float eta_at_slim_disk_boundary = 
        (R / 16.f) * bp->A_lupi * ((0.985f / (R + (5.f / 8.f) * bp->B_lupi)) + 
                                (0.015f / (R + (5.f / 8.f) * bp->C_lupi)));

  bp->jet_is_isotropic = 
    parser_get_param_int(params, "ObsidianAGN:jet_is_isotropic");

  bp->jet_minimum_reservoir_mass
      = parser_get_param_float(params, 
                               "ObsidianAGN:jet_minimum_reservoir_mass_Msun");
  bp->jet_minimum_reservoir_mass /= bp->mass_to_solar_mass;

  /* 2 is along accreted angular momentum direction */
  bp->default_dir_flag =
      parser_get_opt_param_int(params, "ObsidianAGN:default_dir_flag", 2);

  /* We need to keep epsilon_r continuous over all M_dot,BH/M_dot,Edd */
  bp->epsilon_r = eta_at_slim_disk_boundary;
  if (bp->epsilon_r > 1.f) error("Somehow epsilon_r is greater than 1.0.");

  bp->adaf_coupling = 
      parser_get_param_float(params, "ObsidianAGN:adaf_coupling");
  bp->adaf_z_scaling = 
      parser_get_opt_param_float(params, "ObsidianAGN:adaf_z_scaling", 0.f);
  bp->quasar_coupling = 
      parser_get_param_float(params, "ObsidianAGN:quasar_coupling");
  bp->slim_disk_coupling = parser_get_opt_param_float(params, 
      "ObsidianAGN:slim_disk_coupling", bp->quasar_coupling);

  /* These are for momentum constrained winds */
  bp->quasar_wind_momentum_flux = parser_get_opt_param_float(params, 
      "ObsidianAGN:quasar_wind_momentum_flux", 20.f);
  bp->quasar_wind_speed = 
      parser_get_param_float(params, "ObsidianAGN:quasar_wind_speed_km_s");
  bp->quasar_wind_speed *= bp->kms_to_internal;

  bp->quasar_wind_mass_loading = bp->quasar_wind_momentum_flux * 
        bp->quasar_coupling * bp->epsilon_r *
        (phys_const->const_speed_light_c / bp->quasar_wind_speed);
  bp->quasar_f_accretion = 1.f / (1.f + bp->quasar_wind_mass_loading);

  const double slim_disk_wind_momentum_flux = parser_get_opt_param_float(params, 
      "ObsidianAGN:slim_disk_wind_momentum_flux", bp->quasar_wind_momentum_flux);

  bp->slim_disk_wind_speed = parser_get_opt_param_float(params, 
      "ObsidianAGN:slim_disk_wind_speed_km_s", 
      bp->quasar_wind_speed / bp->kms_to_internal);
  bp->slim_disk_wind_speed *= bp->kms_to_internal;

  /* Set the slim disk mass loading to be continuous at the
   * eta upper boundary. Compute the phi term to solve for the 
   * accretion fraction */
  bp->slim_disk_phi = 
      slim_disk_wind_momentum_flux * bp->slim_disk_coupling * 
          (phys_const->const_speed_light_c / bp->slim_disk_wind_speed);
  const double slim_disk_wind_mass_loading = bp->slim_disk_phi * bp->epsilon_r;

  bp->slim_disk_jet_active =
        parser_get_param_int(params, "ObsidianAGN:slim_disk_jet_active");

  bp->adaf_disk_efficiency =
        parser_get_param_float(params, "ObsidianAGN:adaf_disk_efficiency");

  bp->adaf_kick_factor =
      parser_get_opt_param_float(params,
                                 "ObsidianAGN:adaf_kick_factor",
                                 0.5f);
  if (bp->adaf_kick_factor < 0.f || bp->adaf_kick_factor > 1.f) {
    error("adaf_kick_factor must be >= 0 and <= 1.");
  }

  bp->adaf_wind_speed = parser_get_opt_param_float(params, 
        "ObsidianAGN:adaf_wind_speed_km_s", -1.f);
  bp->adaf_wind_speed *= bp->kms_to_internal;

  const float jet_subgrid_mass_loading
      = 2.f * bp->jet_efficiency *
          (phys_const->const_speed_light_c / bp->jet_subgrid_velocity) *
          (phys_const->const_speed_light_c / bp->jet_subgrid_velocity);
  float adaf_wind_mass_loading = 0.f;

  if (bp->adaf_wind_speed < 0.f) {
    bp->adaf_f_accretion = bp->quasar_f_accretion;

    /* Ignore subgrid mass loading in this case */
    double adaf_wind_speed_inv = (1. / bp->adaf_f_accretion) - 1.;
    if (adaf_wind_speed_inv < 0.) {
      error("adaf_wind_speed_inv is negative! adaf_f_accretion=%g",
            bp->adaf_f_accretion);
    }
    adaf_wind_speed_inv /= 2. * bp->adaf_coupling * bp->adaf_disk_efficiency;
    adaf_wind_speed_inv = sqrt(adaf_wind_speed_inv);
    adaf_wind_speed_inv /= phys_const->const_speed_light_c;
    bp->adaf_wind_speed = 1. / adaf_wind_speed_inv;
  }
  else {
    adaf_wind_mass_loading = 2.f * bp->adaf_coupling * bp->adaf_disk_efficiency;
    adaf_wind_mass_loading *= pow(
      phys_const->const_speed_light_c / bp->adaf_wind_speed,
      2.f
    );
    bp->adaf_f_accretion = 1.f /
          (1.f + jet_subgrid_mass_loading + adaf_wind_mass_loading);
  }

  bp->adaf_maximum_temperature =
        parser_get_opt_param_float(params, 
                                   "ObsidianAGN:adaf_maximum_temperature_K", 
                                   5.e7f);
  bp->adaf_maximum_temperature *= T_K_to_int;

  bp->adaf_heating_n_H_threshold_cgs = parser_get_opt_param_float(
      params, "ObsidianAGN:adaf_heating_n_H_threshold_cgs", 0.13f);

  bp->adaf_heating_T_threshold_cgs = parser_get_opt_param_float(
      params, "ObsidianAGN:adaf_heating_T_threshold_cgs", 5.0e5f);

  bp->adaf_mass_limit = 
      parser_get_param_float(params, "ObsidianAGN:adaf_mass_limit_Msun");
  bp->adaf_mass_limit /= bp->mass_to_solar_mass;

  bp->adaf_cooling_shutoff_factor =
      parser_get_opt_param_float(params, 
                                 "ObsidianAGN:adaf_cooling_shutoff_factor",
                                 -1.f);

  /* Always use nibbling in Obsidian */
  bp->use_nibbling = 1;

  bp->bondi_use_all_gas =
        parser_get_opt_param_int(params,
                                 "ObsidianAGN:bondi_use_all_gas",
                                 0);

  bp->bondi_alpha =
        parser_get_opt_param_float(params,
                                   "ObsidianAGN:bondi_alpha",
                                   1.f);

  bp->minimum_black_hole_mass_unresolved = 
      parser_get_param_float(params, 
          "ObsidianAGN:minimum_black_hole_mass_unresolved_Msun");
  bp->minimum_black_hole_mass_unresolved /= bp->mass_to_solar_mass;

  /* Reposition parameters --------------------------------- */

  bp->max_reposition_mass =
      parser_get_param_float(params, "ObsidianAGN:max_reposition_mass") *
      phys_const->const_solar_mass;
  bp->max_reposition_distance_ratio =
      parser_get_param_float(params, "ObsidianAGN:max_reposition_distance_ratio");

  bp->with_reposition_velocity_threshold = parser_get_param_int(
      params, "ObsidianAGN:with_reposition_velocity_threshold");

  if (bp->with_reposition_velocity_threshold) {
    bp->max_reposition_velocity_ratio = parser_get_param_float(
        params, "ObsidianAGN:max_reposition_velocity_ratio");

    /* Prevent nonsensical input */
    if (bp->max_reposition_velocity_ratio <= 0)
      error("max_reposition_velocity_ratio must be positive, not %f.",
            bp->max_reposition_velocity_ratio);

    bp->min_reposition_velocity_threshold = parser_get_param_float(
        params, "ObsidianAGN:min_reposition_velocity_threshold");
    /* Convert from km/s to internal units */
    bp->min_reposition_velocity_threshold *=
        (1e5 / (us->UnitLength_in_cgs / us->UnitTime_in_cgs));
  }

  bp->set_reposition_speed =
      parser_get_param_int(params, "ObsidianAGN:set_reposition_speed");

  if (bp->set_reposition_speed) {
    bp->reposition_coefficient_upsilon = parser_get_param_float(
        params, "ObsidianAGN:reposition_coefficient_upsilon");

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
        parser_get_param_float(params, "ObsidianAGN:reposition_reference_mass") *
        phys_const->const_solar_mass;
    bp->reposition_exponent_mass = parser_get_opt_param_float(
        params, "ObsidianAGN:reposition_exponent_mass", 2.0);
    bp->reposition_reference_n_H =
        parser_get_param_float(params, "ObsidianAGN:reposition_reference_n_H");
    bp->reposition_exponent_n_H = parser_get_opt_param_float(
        params, "ObsidianAGN:reposition_exponent_n_H", 1.0);
  }

  bp->correct_bh_potential_for_repositioning =
      parser_get_param_int(params, "ObsidianAGN:with_potential_correction");

  /* Merger parameters ------------------------------------- */

  bp->minor_merger_threshold =
      parser_get_param_float(params, "ObsidianAGN:threshold_minor_merger");

  bp->major_merger_threshold =
      parser_get_param_float(params, "ObsidianAGN:threshold_major_merger");

  char temp2[40];
  parser_get_param_string(params, "ObsidianAGN:merger_threshold_type", temp2);
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
      parser_get_param_float(params, "ObsidianAGN:merger_max_distance_ratio");

  /* ---- Black hole time-step properties ------------------ */

  const double Myr_in_cgs = 1e6 * 365.25 * 24. * 60. * 60.;

  const double time_step_min_Myr = parser_get_opt_param_float(
      params, "ObsidianAGN:minimum_timestep_Myr", FLT_MAX);

  bp->time_step_min = time_step_min_Myr * Myr_in_cgs /
                      units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Common conversion factors ----------------------------- */

  /* Calculate conversion factor from rho to n_H.
   * Note this assumes primoridal abundance */
  const double X_H = hydro_props->hydrogen_mass_fraction;
  bp->rho_to_n_cgs =
      (X_H / m_p) * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

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

  if (engine_rank == 0) {
    message("Black holes kernel: %s with eta=%f (%.2f neighbours).", 
            kernel_name, bp->eta_neighbours, bp->target_neighbours);

    message("Black holes relative tolerance in h: %.5f (+/- %.4f neighbours).",
            bp->h_tolerance, bp->delta_neighbours);

    message("Black hole model is Rennehan+24");
    message("Black hole jet velocity is %g km/s",
            bp->jet_velocity / bp->kms_to_internal);
    if (bp->jet_loading_type == BH_jet_momentum_loaded) {
      message("Black hole jet loading (momentum) is %g", 
              bp->jet_mass_loading);
    } else {
      message("Black hole jet loading (energy) is %g",
              bp->jet_mass_loading);
    }
    message("Black hole subgrid jet velocity is %g km/s",
              bp->jet_subgrid_velocity / bp->kms_to_internal);
    message("Black hole subgrid jet loading (energy) is %g", 
              jet_subgrid_mass_loading);
    message("Black hole jet efficiency is %g", 
            bp->jet_efficiency);
    message("Black hole quasar radiative efficiency is %g",
            bp->epsilon_r);
    message("Black hole quasar wind speed is %g km/s",
            bp->quasar_wind_speed / bp->kms_to_internal);
    message("Black hole quasar mass loading (momentum) is %g", 
            bp->quasar_wind_mass_loading);
    message("Black hole quasar f_accretion is %g", 
            bp->quasar_f_accretion);
    message("Black hole slim disk wind speed is %g km/s",
            bp->slim_disk_wind_speed / bp->kms_to_internal);
    message("Black hole slim disk mass loading (momentum) is %g "
            "(at the eta=%g boundary)", 
            slim_disk_wind_mass_loading, bp->epsilon_r);
    message("Black hole ADAF mass loading (energy) is %g",
            adaf_wind_mass_loading);
    message("Black hole ADAF f_accretion is %g",
            bp->adaf_f_accretion);
    message("Black hole ADAF v_wind is %g km/s (thermally dumped)",
            bp->adaf_wind_speed / bp->kms_to_internal);
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

#endif /* SWIFT_OBSIDIAN_BLACK_HOLES_PROPERTIES_H */
