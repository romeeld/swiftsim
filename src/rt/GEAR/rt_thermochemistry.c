/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#include "rt_thermochemistry.h"
#include "rt_grackle_utils.h"
#include "rt_interaction_cross_sections.h"
#include "rt_interaction_rates.h"
#include "rt_ionization_equilibrium.h"
#include "rt_unphysical.h"
#include "rt_getters.h"

/**
 * @file src/rt/GEAR/rt_thermochemistry.h
 * @brief Main header file for the GEAR M1 closure radiative transfer scheme
 * thermochemistry related functions.
 */

/**
 * @brief initialize particle quantities relevant for the thermochemistry.
 *
 * @param p part to work with
 * @param rt_props rt_properties struct
 * @param hydro_props hydro properties struct
 * @param phys_const physical constants struct
 * @param us unit system struct
 * @param cosmo cosmology struct
 */
void rt_tchem_first_init_part(
    struct part* restrict p, struct xpart* restrict xp, const struct rt_props* rt_props,
    const struct hydro_props* hydro_props,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cooling_function_data* cooling,
    const struct cosmology* restrict cosmo) {

  if (rt_props->set_equilibrium_initial_ionization_mass_fractions) {
    float XHI, XHII, XHeI, XHeII, XHeIII;
    rt_ion_equil_get_mass_fractions(&XHI, &XHII, &XHeI, &XHeII, &XHeIII, p,
                                    rt_props, hydro_props, phys_const, us,
                                    cosmo);
    p->rt_data.tchem.mass_fraction_HI = XHI;
    p->rt_data.tchem.mass_fraction_HII = XHII;
    p->rt_data.tchem.mass_fraction_HeI = XHeI;
    p->rt_data.tchem.mass_fraction_HeII = XHeII;
    p->rt_data.tchem.mass_fraction_HeIII = XHeIII;
  } else if (rt_props->set_initial_ionization_mass_fractions) {
    p->rt_data.tchem.mass_fraction_HI = rt_props->mass_fraction_HI_init;
    p->rt_data.tchem.mass_fraction_HII = rt_props->mass_fraction_HII_init;
    p->rt_data.tchem.mass_fraction_HeI = rt_props->mass_fraction_HeI_init;
    p->rt_data.tchem.mass_fraction_HeII = rt_props->mass_fraction_HeII_init;
    p->rt_data.tchem.mass_fraction_HeIII = rt_props->mass_fraction_HeIII_init;
  }

    /* Initialize the cooling initial particle quantities. */
    cooling_first_init_part(phys_const, us, hydro_props, cosmo, cooling, p, xp);

    /* Here for GEARRT test. We need to reset the value for cooling data.
     * TODO :And we need to restructure it when we do the cosmological run. */    
    xp->cooling_data.HI_frac = p->rt_data.tchem.mass_fraction_HI;
    xp->cooling_data.HII_frac = p->rt_data.tchem.mass_fraction_HII;
    xp->cooling_data.HeI_frac = p->rt_data.tchem.mass_fraction_HeI;
    xp->cooling_data.HeII_frac = p->rt_data.tchem.mass_fraction_HeII;
    xp->cooling_data.HeIII_frac = p->rt_data.tchem.mass_fraction_HeIII;
    xp->cooling_data.e_frac = xp->cooling_data.HII_frac +
                            0.25 * xp->cooling_data.HeII_frac +
                            0.5 * xp->cooling_data.HeIII_frac;
}

/**
 * @brief Main function for the thermochemistry step.
 *
 * @param p Particle to work on.
 * @param xp Pointer to the particle' extended data.
 * @param rt_props RT properties struct
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param dt The time-step of this particle.
 * @param depth recursion depth
 */
INLINE void rt_do_thermochemistry(
    struct part* restrict p, struct xpart* restrict xp,
    struct rt_props* rt_props, const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct phys_const* restrict phys_const,
    const struct cooling_function_data* restrict cooling,
    const struct unit_system* restrict us, const double dt, 
    const double dt_therm, int depth) {
  /* Note: Can't pass rt_props as const struct because of grackle
   * accessinging its properties there */

  /* Nothing to do here? */
  if (rt_props->skip_thermochemistry) return;
  if (dt == 0.) return;

  /* This is where the fun begins */
  /* ---------------------------- */

  /* initialize data so it'll be in scope */
  //grackle_field_data particle_grackle_data;

  gr_float density = hydro_get_physical_density(p, cosmo);

  /* In rare cases, unphysical solutions can arise with negative densities
   * which won't be fixed in the hydro part until further down the dependency
   * graph. Also, we can have vacuum, in which case we have nothing to do here.
   * So exit early if that is the case. */
  if (density <= 0.) return;

  const float u_minimal = hydro_props->minimal_internal_energy;
#ifdef GIZMO_MFV_SPH
  /* Physical internal energy */
  gr_float internal_energy_phys =
      hydro_get_physical_internal_energy(p, xp, cosmo);
  gr_float internal_energy = max(internal_energy_phys, u_minimal);

  const float u_old = internal_energy;
#else
  /* Get physical internal energy */
  const float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);

  gr_float internal_energy_phys = hydro_get_physical_internal_energy(p, xp, cosmo);

  gr_float internal_energy = max(internal_energy_phys, u_minimal);

  const float u_old = internal_energy;
#endif
  
  /* initialize data to send to grackle */
  gr_float *species_densities;
  species_densities = (gr_float *)calloc(N_SPECIES, sizeof(gr_float));
  grackle_field_data data;

  /* load particle information from particle to grackle data */
  cooling_copy_to_grackle(&data, us, cosmo, cooling, p, xp, dt, 0., species_densities);

  float radiation_energy_density[RT_NGROUPS];
  rt_part_get_physical_radiation_energy_density(p, radiation_energy_density, cosmo);

  /* TODO: put the iact_rates to the cooling grackle data. */
  gr_float iact_rates[5];
  rt_get_interaction_rates_for_grackle(
      iact_rates, radiation_energy_density, species_densities,
      rt_props->average_photon_energy, rt_props->energy_weighted_cross_sections,
      rt_props->number_weighted_cross_sections, phys_const, us);

  /* TODO: currently manually add iact_rates in grackle data field here. */
  data.RT_heating_rate = &iact_rates[0];
  data.RT_HI_ionization_rate = &iact_rates[1];
  data.RT_HeI_ionization_rate = &iact_rates[2];
  data.RT_HeII_ionization_rate = &iact_rates[3];
  data.RT_H2_dissociation_rate = &iact_rates[4];

  /* solve chemistry */
  /* Note: `grackle_rates` is a global variable defined by grackle itself.
   * Using a manually allocd and initialized variable here fails with MPI
   * for some reason. */
  if (local_solve_chemistry(
          &rt_props->grackle_chemistry_data, &rt_props->grackle_chemistry_rates,
          &rt_props->grackle_units, &data, dt) == 0)
    error("Error in solve_chemistry.");
  
  /* copy from grackle data to particle */
  cooling_copy_from_grackle(&data, p, xp, cooling, species_densities[12]);

  /* copy updated grackle data to particle */
  /* update particle internal energy. Grackle had access by reference
   * to internal_energy */
  internal_energy_phys = data.internal_energy[0];
  const float u_new = max(internal_energy_phys, u_minimal);

  /* Re-do thermochemistry? */
  if ((rt_props->max_tchem_recursion > depth) &&
      (fabsf(u_old - u_new) > 0.1 * u_old)) {
    /* Note that grackle already has internal "10% rules". But sometimes, they
     * may not suffice. */
    //rt_clean_grackle_fields(&particle_grackle_data);
    cooling_grackle_free_data(&data);
    free(species_densities);
    rt_do_thermochemistry(p, xp, rt_props, cosmo, hydro_props, phys_const, cooling, us,
                          0.5 * dt, 0.5 * dt_therm, depth + 1);
    rt_do_thermochemistry(p, xp, rt_props, cosmo, hydro_props, phys_const, cooling, us,
                          0.5 * dt, 0.5 * dt_therm, depth + 1);
    return;
  }

  /* If we're good, update the particle data from grackle results */
#ifdef GIZMO_MFV_SPH
  hydro_set_physical_internal_energy(p, u_new);
#else
  /* compute the heating/cooling due to the thermochemistry */
  float cool_du_dt = (u_new - u_old) / dt_therm;
  
  /* check whether the the thermochemistry heating/cooling is larger
   * than du/dt of the particle. If it is, directly set the new internal energy 
   * of the particle, and set du/dt = 0.*/
  if (fabsf(cool_du_dt) > fabsf(hydro_du_dt)){
    hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
  
    hydro_set_physical_internal_energy_dt(p, cosmo, 0.);
  } else {
  /* If it isn't, ignore the radiative cooling and apply only hydro du/dt. */
    hydro_set_physical_internal_energy_dt(p, cosmo, hydro_du_dt);
   }
#endif

  /* Update mass fractions */
  const gr_float one_over_rho = 1. / density;
  p->rt_data.tchem.mass_fraction_HI =
      data.HI_density[0] * one_over_rho;
  p->rt_data.tchem.mass_fraction_HII =
      data.HII_density[0] * one_over_rho;
  p->rt_data.tchem.mass_fraction_HeI =
      data.HeI_density[0] * one_over_rho;
  p->rt_data.tchem.mass_fraction_HeII =
      data.HeII_density[0] * one_over_rho;
  p->rt_data.tchem.mass_fraction_HeIII =
      data.HeIII_density[0] * one_over_rho;

  /* Update radiation fields */
  /* First get absorption rates at the start and the end of the step */
  double absorption_rates[RT_NGROUPS];
  rt_get_absorption_rates(
      absorption_rates, species_densities, rt_props->average_photon_energy,
      rt_props->number_weighted_cross_sections, phys_const, us);

  gr_float species_densities_new[6];
  species_densities_new[0] = data.HI_density[0];
  species_densities_new[1] = data.HII_density[0];
  species_densities_new[2] = data.HeI_density[0];
  species_densities_new[3] = data.HeII_density[0];
  species_densities_new[4] = data.HeIII_density[0];
  species_densities_new[5] = data.e_density[0];
  double absorption_rates_new[RT_NGROUPS];
  rt_get_absorption_rates(absorption_rates_new, species_densities_new,
                          rt_props->average_photon_energy,
                          rt_props->number_weighted_cross_sections, phys_const,
                          us);

  /* Now remove absorbed radiation */
  for (int g = 0; g < RT_NGROUPS; g++) {
    const float E_old = p->rt_data.radiation[g].energy_density;
    double f = dt * 0.5 * (absorption_rates[g] + absorption_rates_new[g]);
    f = min(1., f);
    f = max(0., f);
    p->rt_data.radiation[g].energy_density *= (1. - f);
    for (int i = 0; i < 3; i++) {
      p->rt_data.radiation[g].flux[i] *= (1. - f);
    }

    rt_check_unphysical_state(&p->rt_data.radiation[g].energy_density,
                              p->rt_data.radiation[g].flux, E_old,
                              /*callloc=*/2);
  }

  /* Clean up after yourself. */
  cooling_grackle_free_data(&data);
  free(species_densities);
}

/**
 * @brief Computes an upper boundary for the thermochemistry/cooling
 * time.
 *
 * @param p Particle to work on.
 * @param xp Pointer to the particle' extended data.
 * @param rt_props RT properties struct
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 */
float rt_tchem_get_tchem_time(
    const struct part* restrict p, const struct xpart* restrict xp,
    struct rt_props* rt_props, const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct phys_const* restrict phys_const,
    const struct cooling_function_data* restrict cooling,
    const struct unit_system* restrict us) {
  /* Note: Can't pass rt_props as const struct because of grackle
   * accessinging its properties there */


  /* initialize data to send to grackle */
  gr_float *species_densities;
  species_densities = (gr_float *)calloc(N_SPECIES, sizeof(gr_float));
  grackle_field_data data;

  float radiation_energy_density[RT_NGROUPS];
  rt_part_get_physical_radiation_energy_density(p, radiation_energy_density, cosmo);

  gr_float iact_rates[5];
  rt_get_interaction_rates_for_grackle(
      iact_rates, radiation_energy_density, species_densities,
      rt_props->average_photon_energy, rt_props->energy_weighted_cross_sections,
      rt_props->number_weighted_cross_sections, phys_const, us);

  /* load particle information from particle to grackle data */
  cooling_copy_to_grackle(&data, us, cosmo, cooling, p, xp, 0., 0., species_densities);

  /* TODO: currently manually add iact_rates in grackle data field here. */
  data.RT_heating_rate = &iact_rates[0];
  data.RT_HI_ionization_rate = &iact_rates[1];
  data.RT_HeI_ionization_rate = &iact_rates[2];
  data.RT_HeII_ionization_rate = &iact_rates[3];
  data.RT_H2_dissociation_rate = &iact_rates[4];

  /* Compute 'cooling' time */
  /* Note: grackle_rates is a global variable defined by grackle itself.
   * Using a manually allocd and initialized variable here fails with MPI
   * for some reason. */
  gr_float tchem_time;
  if (local_calculate_cooling_time(
          &rt_props->grackle_chemistry_data, &rt_props->grackle_chemistry_rates,
          &rt_props->grackle_units, &data, &tchem_time) == 0)
    error("Error in calculate_cooling_time.");

  /* Clean up after yourself. */
  cooling_grackle_free_data(&data);
  free(species_densities);

  return (float)tchem_time;
}

/**
 * @brief Main function for the thermochemistry step when coupling with 
 * subgrid physics.
 *
 * @param p Particle to work on.
 * @param xp Pointer to the particle' extended data.
 * @param rt_props RT properties struct
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param dt The time-step of this particle.
 * @param depth recursion depth
 */
INLINE void rt_do_thermochemistry_with_subgrid(
    struct part* restrict p, struct xpart* restrict xp,
    struct rt_props* rt_props, const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct entropy_floor_properties* floor_props,
    const struct phys_const* restrict phys_const,
    const struct cooling_function_data* restrict cooling,
    const struct unit_system* restrict us, const double dt,
    const double dt_therm, int depth) {
  /* Note: Can't pass rt_props as const struct because of grackle
   * accessinging its properties there */

  /* Nothing to do here? */
  if (rt_props->skip_thermochemistry) return;
  if (dt == 0.f || dt_therm == 0.f) return;

  /* Never cool if there is a cooling shut off, let the hydro do its magic */
  if (p->feedback_data.cooling_shutoff_delay_time > 0.f) {
    /* These need to be set for the shut off particles */
    p->cooling_data.subgrid_dens = hydro_get_physical_density(p, cosmo);
    p->cooling_data.subgrid_temp = 0.;

    return;
  }

  /* No cooling if particle is decoupled, but do the firehose model */
  if (p->decoupled) {
    /* Make sure these are always set for the wind particles */
    p->cooling_data.subgrid_dens = hydro_get_physical_density(p, cosmo);
    p->cooling_data.subgrid_temp = 0.;

    firehose_cooling_and_dust(phys_const, us, cosmo, hydro_props, cooling, p, xp, dt);
    return;
  }
  
  /* Update the subgrid properties */
  cooling_set_particle_subgrid_properties( phys_const, us,
	  cosmo, hydro_props, floor_props, cooling, p, xp);

  /* Compute the entropy floor */
  const double T_floor = entropy_floor_temperature(p, cosmo, floor_props);
  const double u_floor = cooling_convert_temp_to_u(T_floor, xp->cooling_data.e_frac, cooling, p);

  /* If it's eligible for SF and metal-free, crudely self-enrich to very small level; needed to kick-start Grackle dust */
  double init_dust_to_gas=0.2;  // arbitrary choice to kick-start dust model
  if (p->cooling_data.subgrid_temp > 0.f && chemistry_get_total_metal_mass_fraction_for_cooling(p) < (1.-init_dust_to_gas) * cooling->self_enrichment_metallicity) {
    /* Fraction of mass in stars going SNe */
    //const float yield = 0.02;
    /* Compute increase in metal mass fraction due to self-enrichment from own SFR */
    //const double delta_metal_frac = p->sf_data.SFR * dt * yield / p->mass;
    /* Set metal fraction to floor value when star-forming */
    p->chemistry_data.metal_mass_fraction_total += (1.-init_dust_to_gas) * cooling->self_enrichment_metallicity;
    p->cooling_data.dust_mass += init_dust_to_gas * cooling->self_enrichment_metallicity * p->mass;
    /* SolarAbundances has He as element 0, while chemistry_element struct has H as element 0, hence an offset of 1 */
    double solar_met_total=0.f;
    for (int i = 1; i < chemistry_element_count; i++)
       if (i > chemistry_element_He) solar_met_total += cooling->chemistry.SolarAbundances[i-1];
    /* Distribute the self-enrichment metallicity among elements assuming solar abundance ratios*/
    for (int i = 1; i < chemistry_element_count; i++) {
      if (i > chemistry_element_He) {
  	 p->chemistry_data.metal_mass_fraction[i] +=
	    (1.-init_dust_to_gas) * cooling->chemistry.SolarAbundances[i-1] * cooling->self_enrichment_metallicity / solar_met_total;  // fraction of gas mass in each element
  	 p->cooling_data.dust_mass_fraction[i] +=
	    init_dust_to_gas * cooling->chemistry.SolarAbundances[i-1] / solar_met_total;  // fraction of dust mass in each element
      }
    }
  }

  /* This is where the fun begins */
  /* ---------------------------- */

  /* initialize data so it'll be in scope */
  //grackle_field_data particle_grackle_data;

  gr_float density = hydro_get_physical_density(p, cosmo);

  /* In rare cases, unphysical solutions can arise with negative densities
   * which won't be fixed in the hydro part until further down the dependency
   * graph. Also, we can have vacuum, in which case we have nothing to do here.
   * So exit early if that is the case. */
  if (density <= 0.) return;

  const float u_minimal = hydro_props->minimal_internal_energy;
#ifdef GIZMO_MFV_SPH
  /* Physical internal energy */
  gr_float internal_energy_phys =
      hydro_get_physical_internal_energy(p, xp, cosmo);
  gr_float internal_energy = max(internal_energy_phys, u_minimal);

  const float u_old = internal_energy;
#else
  /* Get physical internal energy */
  const float hydro_du_dt = hydro_get_physical_internal_energy_dt(p, cosmo);

  gr_float internal_energy_phys = hydro_get_physical_internal_energy(p, xp, cosmo);

  gr_float internal_energy = max(internal_energy_phys, u_minimal);

  const float u_old = internal_energy;
#endif

  /* initialize data to send to grackle */
  gr_float *species_densities;
  species_densities = (gr_float *)calloc(N_SPECIES, sizeof(gr_float));
  grackle_field_data data;

  /* load particle information from particle to grackle data */
  cooling_copy_to_grackle(&data, us, cosmo, cooling, p, xp, dt, 0., species_densities);

  float radiation_energy_density[RT_NGROUPS];
  rt_part_get_physical_radiation_energy_density(p, radiation_energy_density, cosmo);

  /* TODO: put the iact_rates to the cooling grackle data. */
  gr_float iact_rates[5];
  rt_get_interaction_rates_for_grackle(
      iact_rates, radiation_energy_density, species_densities,
      rt_props->average_photon_energy, rt_props->energy_weighted_cross_sections,
      rt_props->number_weighted_cross_sections, phys_const, us);

  /* TODO: currently manually add iact_rates in grackle data field here. */
  data.RT_heating_rate = &iact_rates[0];
  data.RT_HI_ionization_rate = &iact_rates[1];
  data.RT_HeI_ionization_rate = &iact_rates[2];
  data.RT_HeII_ionization_rate = &iact_rates[3];
  data.RT_H2_dissociation_rate = &iact_rates[4];

  /* solve chemistry */
  /* Note: `grackle_rates` is a global variable defined by grackle itself.
   * Using a manually allocd and initialized variable here fails with MPI
   * for some reason. */
  if (local_solve_chemistry(
          &rt_props->grackle_chemistry_data, &rt_props->grackle_chemistry_rates,
          &rt_props->grackle_units, &data, dt) == 0)
    error("Error in solve_chemistry.");

  /* copy from grackle data to particle */
  cooling_copy_from_grackle(&data, p, xp, cooling, species_densities[12]);

  /* copy updated grackle data to particle */
  /* update particle internal energy. Grackle had access by reference
   * to internal_energy */
  internal_energy_phys = data.internal_energy[0];
  float u_new = max(internal_energy_phys, u_minimal);

  /* Re-do thermochemistry? */
  if ((rt_props->max_tchem_recursion > depth) &&
      (fabsf(u_old - u_new) > 0.1 * u_old)) {
    /* Note that grackle already has internal "10% rules". But sometimes, they
     * may not suffice. */
    //rt_clean_grackle_fields(&particle_grackle_data);
    cooling_grackle_free_data(&data);
    free(species_densities);
    rt_do_thermochemistry_with_subgrid(p, xp, rt_props, cosmo, hydro_props, floor_props, phys_const, cooling, us,
                          0.5 * dt, 0.5 * dt_therm, depth + 1);
    rt_do_thermochemistry_with_subgrid(p, xp, rt_props, cosmo, hydro_props, floor_props, phys_const, cooling, us,
                          0.5 * dt, 0.5 * dt_therm, depth + 1);
    return;
  }

  /* Assign new thermal energy to particle */
  float cool_du_dt = 0.;

  if (p->cooling_data.subgrid_temp == 0.) {
    /* Normal cooling; check that we are not going to go below any of the limits */
    if (u_new > GRACKLE_HEATLIM * u_old) u_new = GRACKLE_HEATLIM * u_old;
    if (u_new < GRACKLE_COOLLIM * u_old) u_new = GRACKLE_COOLLIM * u_old;
    u_new = max(u_new, u_floor);

    /* If we're good, update the particle data from grackle results */
#ifdef GIZMO_MFV_SPH
    hydro_set_physical_internal_energy(p, u_new);
#else
    /* compute the heating/cooling due to the thermochemistry */
    cool_du_dt = (u_new - u_old) / dt_therm;

    /* check whether the the thermochemistry heating/cooling is larger
     * than du/dt of the particle. If it is, directly set the new internal energy
     * of the particle, and set du/dt = 0.*/
    if (fabsf(cool_du_dt) > fabsf(hydro_du_dt)){
      hydro_set_physical_internal_energy(p, xp, cosmo, u_new);

      hydro_set_physical_internal_energy_dt(p, cosmo, 0.);
    } else {
    /* If it isn't, ignore the radiative cooling and apply only hydro du/dt. */
      hydro_set_physical_internal_energy_dt(p, cosmo, hydro_du_dt);
     }
#endif
  }
  else {
    /* Particle is in subgrid mode; result is stored in subgrid_temp */
    p->cooling_data.subgrid_temp = cooling_convert_u_to_temp(u_new, xp->cooling_data.e_frac, cooling, p);

    /* Set internal energy time derivative to 0 for overall particle */
    hydro_set_physical_internal_energy_dt(p, cosmo, 0.f);

    /* Force the overall particle to lie on the equation of state */
    hydro_set_physical_internal_energy(p, xp, cosmo, u_floor);
  }

  /* Store the radiated energy */
  xp->cooling_data.radiated_energy -= hydro_get_mass(p) * cool_du_dt * dt_therm;

  /* set subgrid properties for use in SF routine */
  cooling_set_particle_subgrid_properties(
      phys_const, us, cosmo, hydro_props, floor_props, cooling, p, xp);

  /* Update mass fractions */
  const gr_float one_over_rho = 1. / density;
  p->rt_data.tchem.mass_fraction_HI =
      data.HI_density[0] * one_over_rho;
  p->rt_data.tchem.mass_fraction_HII =
      data.HII_density[0] * one_over_rho;
  p->rt_data.tchem.mass_fraction_HeI =
      data.HeI_density[0] * one_over_rho;
  p->rt_data.tchem.mass_fraction_HeII =
      data.HeII_density[0] * one_over_rho;
  p->rt_data.tchem.mass_fraction_HeIII =
      data.HeIII_density[0] * one_over_rho;

  /* Update radiation fields */
  /* First get absorption rates at the start and the end of the step */
  double absorption_rates[RT_NGROUPS];
  rt_get_absorption_rates(
      absorption_rates, species_densities, rt_props->average_photon_energy,
      rt_props->number_weighted_cross_sections, phys_const, us);

  gr_float species_densities_new[6];
  species_densities_new[0] = data.HI_density[0];
  species_densities_new[1] = data.HII_density[0];
  species_densities_new[2] = data.HeI_density[0];
  species_densities_new[3] = data.HeII_density[0];
  species_densities_new[4] = data.HeIII_density[0];
  species_densities_new[5] = data.e_density[0];
  double absorption_rates_new[RT_NGROUPS];
  rt_get_absorption_rates(absorption_rates_new, species_densities_new,
                          rt_props->average_photon_energy,
                          rt_props->number_weighted_cross_sections, phys_const,
                          us);

  /* Now remove absorbed radiation */
  for (int g = 0; g < RT_NGROUPS; g++) {
    const float E_old = p->rt_data.radiation[g].energy_density;
    double f = dt * 0.5 * (absorption_rates[g] + absorption_rates_new[g]);
    f = min(1., f);
    f = max(0., f);
    p->rt_data.radiation[g].energy_density *= (1. - f);
    for (int i = 0; i < 3; i++) {
      p->rt_data.radiation[g].flux[i] *= (1. - f);
    }

    rt_check_unphysical_state(&p->rt_data.radiation[g].energy_density,
                              p->rt_data.radiation[g].flux, E_old,
                              /*callloc=*/2);
  }

  /* Clean up after yourself. */
  cooling_grackle_free_data(&data);
  free(species_densities);
}

