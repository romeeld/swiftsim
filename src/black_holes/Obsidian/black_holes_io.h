/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_OBSIDIAN_BLACK_HOLES_IO_H
#define SWIFT_OBSIDIAN_BLACK_HOLES_IO_H

#include "adiabatic_index.h"
#include "black_holes_part.h"
#include "black_holes_properties.h"
#include "io_properties.h"
#include "kick.h"

/**
 * @brief Specifies which b-particle fields to read from a dataset
 *
 * @param bparts The b-particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void black_holes_read_particles(struct bpart* bparts,
                                              struct io_props* list,
                                              int* num_fields) {

  int num = 0;

  /* List what we want to read */
  list[num] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, bparts, x);
  num++;

  list[num] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, bparts, v);
  num++;

  list[num] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                bparts, mass);
  num++;

  list[num] = io_make_input_field("ParticleIDs", LONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, bparts, id);
  num++;

  list[num] = io_make_input_field("SmoothingLengths", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_LENGTH, bparts, h);
  num++;

  list[num] = io_make_input_field("JetMassReservoirs", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_MASS, bparts, jet_mass_reservoir);
  num++;

  list[num] = io_make_input_field("SubgridMasses", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_MASS, bparts, subgrid_mass);
  num++;

  *num_fields = num;

}

INLINE static void convert_bpart_pos(const struct engine* e,
                                     const struct bpart* bp, double* ret) {

  const struct space* s = e->s;
  if (s->periodic) {
    ret[0] = box_wrap(bp->x[0], 0.0, s->dim[0]);
    ret[1] = box_wrap(bp->x[1], 0.0, s->dim[1]);
    ret[2] = box_wrap(bp->x[2], 0.0, s->dim[2]);
  } else {
    ret[0] = bp->x[0];
    ret[1] = bp->x[1];
    ret[2] = bp->x[2];
  }
  if (e->snapshot_use_delta_from_edge) {
    ret[0] = min(ret[0], s->dim[0] - e->snapshot_delta_from_edge);
    ret[1] = min(ret[1], s->dim[1] - e->snapshot_delta_from_edge);
    ret[2] = min(ret[2], s->dim[2] - e->snapshot_delta_from_edge);
  }
}

INLINE static void convert_bpart_vel(const struct engine* e,
                                     const struct bpart* bp, float* ret) {

  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const struct cosmology* cosmo = e->cosmology;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;
  const float dt_kick_grav_mesh = e->dt_kick_grav_mesh_for_io;

  const integertime_t ti_beg = get_integer_time_begin(ti_current, bp->time_bin);
  const integertime_t ti_end = get_integer_time_end(ti_current, bp->time_bin);

  /* Get time-step since the last kick */
  const float dt_kick_grav =
      kick_get_grav_kick_dt(ti_beg, ti_current, time_base, with_cosmology,
                            cosmo) -
      kick_get_grav_kick_dt(ti_beg, (ti_beg + ti_end) / 2, time_base,
                            with_cosmology, cosmo);

  /* Extrapolate the velocites to the current time */
  const struct gpart* gp = bp->gpart;
  ret[0] = gp->v_full[0] + gp->a_grav[0] * dt_kick_grav;
  ret[1] = gp->v_full[1] + gp->a_grav[1] * dt_kick_grav;
  ret[2] = gp->v_full[2] + gp->a_grav[2] * dt_kick_grav;

  /* Extrapolate the velocites to the current time (mesh forces) */
  ret[0] += gp->a_grav_mesh[0] * dt_kick_grav_mesh;
  ret[1] += gp->a_grav_mesh[1] * dt_kick_grav_mesh;
  ret[2] += gp->a_grav_mesh[2] * dt_kick_grav_mesh;

  /* Conversion from internal to physical units */
  ret[0] *= cosmo->a_inv;
  ret[1] *= cosmo->a_inv;
  ret[2] *= cosmo->a_inv;
}

INLINE static void convert_bpart_potential(const struct engine* e,
                                           const struct bpart* bp, float* ret) {

  if (bp->gpart != NULL)
    ret[0] = gravity_get_comoving_potential(bp->gpart);
  else
    ret[0] = 0.f;
}

INLINE static void convert_bpart_gas_vel(const struct engine* e,
                                         const struct bpart* bp, float* ret) {

  const struct cosmology* cosmo = e->cosmology;

  /* Convert relative velocities to physical units */
  ret[0] = bp->velocity_gas[0] * cosmo->a_inv;
  ret[1] = bp->velocity_gas[1] * cosmo->a_inv;
  ret[2] = bp->velocity_gas[2] * cosmo->a_inv;
}

INLINE static void convert_bpart_gas_circular_vel(const struct engine* e,
                                                  const struct bpart* bp,
                                                  float* ret) {

  const struct cosmology* cosmo = e->cosmology;

  /* Conversion from internal to physical units */
  ret[0] = bp->circular_velocity_gas[0] * cosmo->a_inv;
  ret[1] = bp->circular_velocity_gas[1] * cosmo->a_inv;
  ret[2] = bp->circular_velocity_gas[2] * cosmo->a_inv;
}

INLINE static void convert_bpart_gas_temperatures(const struct engine* e,
                                                  const struct bpart* bp,
                                                  float* ret) {

  const struct black_holes_props* props = e->black_holes_properties;
  const struct cosmology* cosmo = e->cosmology;

  /* Conversion from specific internal energy to temperature */
  ret[0] = bp->internal_energy_gas * cosmo->a_factor_internal_energy /
           props->temp_to_u_factor;
}

/**
 * @brief Specifies which b-particle fields to write to a dataset
 *
 * @param bparts The b-particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 * @param with_cosmology Are we running a cosmological simulation?
 */
INLINE static void black_holes_write_particles(const struct bpart* bparts,
                                               struct io_props* list,
                                               int* num_fields,
                                               int with_cosmology) {

  int num = 0;

  /* List what we want to write */
  list[num] = io_make_output_field_convert_bpart(
      "Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH, 1.f, bparts,
      convert_bpart_pos, "Co-moving position of the particles");
  num++;

  list[num] = io_make_output_field_convert_bpart(
      "Velocities", FLOAT, 3, UNIT_CONV_SPEED, 0.f, bparts, convert_bpart_vel,
      "Peculiar velocities of the particles. This is a * dx/dt where x is the "
      "co-moving position of the particles.");
  num++;

  list[num] =
      io_make_output_field("DynamicalMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                           bparts, mass, "Dynamical masses of the particles");
  num++;

  list[num] =
      io_make_output_field("ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
                           bparts, id, "Unique ID of the particles");
  num++;

  list[num] = io_make_output_field(
      "SmoothingLengths", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, bparts, h,
      "Co-moving smoothing lengths (FWHM of the kernel) of the particles");
  num++;

  list[num] = io_make_output_field("SubgridMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                                 bparts, subgrid_mass,
                                 "Subgrid masses of the particles");
  num++;

  if (with_cosmology) {
    list[num] = io_make_output_field(
        "FormationScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
        formation_scale_factor, "Scale-factors at which the BHs were formed");
  } else {
    list[num] = io_make_output_field("FormationTimes", FLOAT, 1, UNIT_CONV_TIME,
                                   0.f, bparts, formation_time,
                                   "Times at which the BHs were formed");
  }
  num++;

  list[num] = io_make_output_field(
      "GasDensities", FLOAT, 1, UNIT_CONV_DENSITY, -3.f, bparts, rho_gas,
      "Co-moving densities of the gas around the particles");
  num++;

  list[num] = io_make_output_field(
      "GasSoundSpeeds", FLOAT, 1, UNIT_CONV_SPEED,
      -1.5f * hydro_gamma_minus_one, bparts, sound_speed_gas,
      "Co-moving sound-speeds of the gas around the particles");
  num++;

  list[num] = io_make_output_field(
      "JetMassReservoirs", FLOAT, 1, UNIT_CONV_MASS, 0.f, bparts,
      jet_mass_reservoir,
      "Physcial mass contained in the jet reservoir of the particles");
  num++;

  list[num] = io_make_output_field(
      "UnresolvedMassReservoirs", FLOAT, 1, UNIT_CONV_MASS, 0.f, bparts,
      unresolved_mass_reservoir,
      "Physcial mass contained in the unresolved reservoir of the particles");
  num++;

  list[num] = io_make_output_field(
      "AccretionRates", FLOAT, 1, UNIT_CONV_MASS_PER_UNIT_TIME, 0.f, bparts,
      accretion_rate,
      "Physical instantaneous accretion rates of the particles");
  num++;

  list[num] = io_make_output_field(
      "TotalAccretedMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, bparts,
      total_accreted_mass,
      "Total mass accreted onto the particles since its birth");
  num++;

  list[num] = io_make_output_field(
      "CumulativeNumberOfSeeds", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      cumulative_number_seeds,
      "Total number of BH seeds that have merged into this black hole");
  num++;

  list[num] =
      io_make_output_field("NumberOfMergers", INT, 1, UNIT_CONV_NO_UNITS, 0.f,
                           bparts, number_of_mergers,
                           "Number of mergers the black holes went through. "
                           "This does not include the number of mergers "
                           "accumulated by any merged black hole.");
  num++;

  if (with_cosmology) {
    list[num] = io_make_output_field(
        "LastMinorMergerScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
        bparts, last_minor_merger_scale_factor,
        "Scale-factors at which the black holes last had a minor merger.");
  } else {
    list[num] = io_make_output_field(
        "LastMinorMergerTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, bparts,
        last_minor_merger_time,
        "Times at which the black holes last had a minor merger.");
  }
  num++;

  if (with_cosmology) {
    list[num] = io_make_output_field(
        "LastMajorMergerScaleFactors", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
        bparts, last_major_merger_scale_factor,
        "Scale-factors at which the black holes last had a major merger.");
  } else {
    list[num] = io_make_output_field(
        "LastMajorMergerTimes", FLOAT, 1, UNIT_CONV_TIME, 0.f, bparts,
        last_major_merger_time,
        "Times at which the black holes last had a major merger.");
  }
  num++;

  list[num] = io_make_output_field(
      "SwallowedAngularMomenta", FLOAT, 3, UNIT_CONV_ANGULAR_MOMENTUM, 0.f,
      bparts, swallowed_angular_momentum,
      "Physical angular momenta that the black holes have accumulated by "
      "swallowing gas particles.");
  num++;

  list[num] = io_make_output_field_convert_bpart(
      "GasRelativeVelocities", FLOAT, 3, UNIT_CONV_SPEED, 0.f, bparts,
      convert_bpart_gas_vel,
      "Peculiar relative velocities of the gas particles around the black "
      "holes. This is a * dx/dt where x is the co-moving position of the "
      "particles.");
  num++;

  list[num] = io_make_output_field_convert_bpart(
      "GasCircularVelocities", FLOAT, 3, UNIT_CONV_SPEED, 0.f, bparts,
      convert_bpart_gas_circular_vel,
      "Circular velocities of the gas around the black hole at the "
      "smoothing radius. This is j / h_BH, where j is the smoothed, peculiar "
      "specific angular momentum of gas around the black holes, and h_BH is "
      "the smoothing length of each black hole.");
  num++;

  list[num] =
      io_make_output_field("TimeBins", CHAR, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
                           time_bin, "Time-bins of the particles");
  num++;

  list[num] = io_make_output_field(
      "NumberOfSwallows", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      number_of_gas_swallows,
      "Number of gas particles the black holes have swallowed. "
      "This includes the particles swallowed by any of the black holes that "
      "merged into this one.");
  num++;

  list[num] = io_make_output_field(
      "NumberOfDirectSwallows", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      number_of_direct_gas_swallows,
      "Number of gas particles the black holes have swallowed. "
      "This does not include any particles swallowed by any of the black holes "
      "that merged into this one.");
  num++;

  list[num] = io_make_output_field(
      "NumberOfRepositions", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      number_of_repositions,
      "Number of repositioning events the black holes went through. This does "
      "not include the number of reposition events accumulated by any merged "
      "black holes.");
  num++;

  list[num] = io_make_output_field(
      "NumberOfRepositionAttempts", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      number_of_reposition_attempts,
      "Number of time steps in which the black holes had an eligible particle "
      "to reposition to. They may or may not have ended up moving there, "
      "depending on their subgrid mass and on whether these particles were at "
      "a lower or higher potential than the black holes themselves. It does "
      "not include attempted repositioning events accumulated by any merged "
      "black holes.");
  num++;

  list[num] = io_make_output_field(
      "NumberOfTimeSteps", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      number_of_time_steps,
      "Total number of time steps at which the black holes were active.");
  num++;

  list[num] = io_make_output_field(
      "SubgridSoundSpeeds", FLOAT, 1, UNIT_CONV_SPEED, 0.f, bparts,
      sound_speed_subgrid_gas,
      "Physical subgrid sound-speeds used in the subgrid-Bondi model.");
  num++;

  list[num] = io_make_output_field(
      "BirthGasDensities", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, bparts,
      formation_gas_density,
      "Physical densities of the converted part at the time of birth. "
      "We store the physical density at the birth redshift, no conversion is "
      "needed.");
  num++;

  list[num] = io_make_output_field(
      "AccretedAngularMomenta", FLOAT, 3, UNIT_CONV_ANGULAR_MOMENTUM, 0.f,
      bparts, accreted_angular_momentum,
      "Physical angular momenta that the black holes have accumulated through "
      "subgrid accretion.");
  num++;

  list[num] = io_make_output_field(
      "NumberOfGasNeighbours", INT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      num_ngbs,
      "Integer number of gas neighbour particles within the black hole "
      "kernels.");
  num++;

  list[num] = io_make_output_field(
      "LastRepositionVelocities", FLOAT, 1, UNIT_CONV_SPEED, 0.f, bparts,
      last_repos_vel,
      "Physical speeds at which the black holes repositioned most recently. "
      "This is 0 for black holes that have never repositioned, or if the "
      "simulation has been run without prescribed repositioning speed.");
  num++;

  list[num] = io_make_output_field_convert_bpart(
      "GasTemperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, bparts,
      convert_bpart_gas_temperatures,
      "Temperature of the gas surrounding the black holes.");
  num++;

  list[num] = io_make_output_field(
      "EddingtonFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, bparts,
      eddington_fraction,
      "Accretion rates of black holes in units of their Eddington rates. "
      "This is based on the unlimited accretion rates, so these fractions "
      "can be above the limiting fEdd.");
  num++;

  list[num] = io_make_output_field_convert_bpart(
      "Potentials", FLOAT, 1, UNIT_CONV_POTENTIAL, -1.f, bparts,
      convert_bpart_potential, "Gravitational potentials of the particles");
  num++;

  *num_fields = num;

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES

  list += *num_fields;
  *num_fields += 4;

  list[0] = io_make_output_field("Num_ngb_density", INT, 1, UNIT_CONV_NO_UNITS,
                                 bparts, num_ngb_density);
  list[1] = io_make_output_field("Num_ngb_force", INT, 1, UNIT_CONV_NO_UNITS,
                                 bparts, num_ngb_force);
  list[2] = io_make_output_field("Ids_ngb_density", LONGLONG,
                                 MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES,
                                 UNIT_CONV_NO_UNITS, bparts, ids_ngbs_density);
  list[3] = io_make_output_field("Ids_ngb_force", LONGLONG,
                                 MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES,
                                 UNIT_CONV_NO_UNITS, bparts, ids_ngbs_force);
#endif
}

#endif /* SWIFT_OBSIDIAN_BLACK_HOLES_IO_H */
