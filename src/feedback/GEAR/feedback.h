/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_FEEDBACK_GEAR_H
#define SWIFT_FEEDBACK_GEAR_H

#include "cosmology.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "stars.h"
#include "stellar_evolution.h"
#include "units.h"

#include <strings.h>

void feedback_update_part(struct part* p, struct xpart* xp,
                          const struct engine* e, const int with_cosmology);

/**
 * @brief Determine the probability of a gas particle being kicked
 *        due to stellar feedback in star forming gas.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param fb_props The feedback properties.
 * @param ti_current The current timestep.
 * @param dt_part The time step of the particle.
 * @param rand_for_sf_wind The random number for the wind generation.
 * @param wind_mass The amount of mass in the wind (code units).
 */
__attribute__((always_inline)) INLINE static double feedback_wind_probability(
    struct part* p, struct xpart* xp, const struct engine* e,
    const struct cosmology* cosmo,
    const struct feedback_props* fb_props,
    const integertime_t ti_current,
    const double dt_part,
    double *rand_for_sf_wind,
    double *wind_mass) {

  return 0.f;
}


/**
 * @brief Kick a gas particle selected for stellar feedback.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param fb_props The feedback properties.
 * @param ti_current The current timestep.
 * @param with_cosmology Is cosmological integration on?
 * @param dt_part The time step of the particle.
 * @param wind_mass The amount of mass in the wind (code units).
 */
__attribute__((always_inline)) INLINE static void feedback_kick_and_decouple_part(
    struct part* p, struct xpart* xp,
    const struct engine* e,
    const struct cosmology* cosmo,
    const struct feedback_props* fb_props,
    const integertime_t ti_current,
    const int with_cosmology,
    const double dt_part,
    const double wind_mass) {};


/**
 * @brief Recouple wind particles.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param with_cosmology Is this a cosmological simulation?
 */
__attribute__((always_inline)) INLINE static void feedback_recouple_part(
    struct part* p, struct xpart* xp, const struct engine* e,
    const int with_cosmology, 
    const struct cosmology* cosmo,
    const struct feedback_props* fb_props) {}

/**
 * @brief Sets the wind direction vector for feedback kicks
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param with_cosmology Is this a cosmological simulation?
 * @param cosmo The cosmology of the simulation.
 * @param fb_props The #feedback_props feedback parameters.
 */
__attribute__((always_inline)) INLINE static void feedback_set_wind_direction(
    struct part* p, struct xpart* xp, const struct engine* e,
    const int with_cosmology, 
    const struct cosmology* cosmo,
    const struct feedback_props* fb_props) {}

void feedback_reset_part(struct part* p, struct xpart* xp);

void feedback_will_do_feedback(
    struct spart* sp, const struct feedback_props* feedback_props,
    const int with_cosmology, const struct cosmology* cosmo, const double time,
    const struct unit_system* us, const struct phys_const* phys_const,
    const integertime_t ti_current, const double time_base);

int feedback_is_active(const struct spart* sp, const struct engine* e);

/**
 * @brief Should this particle be doing any DM looping?
 *
 * @param sp The #spart.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static int stars_dm_loop_is_active(
    const struct spart* sp, const struct engine* e) {
  /* No */
  return 0;
}

double feedback_get_enrichment_timestep(const struct spart* sp,
                                        const int with_cosmology,
                                        const struct cosmology* cosmo,
                                        const double time,
                                        const double dt_star);
void feedback_init_spart(struct spart* sp);

void feedback_init_after_star_formation(
    struct spart* sp, const struct feedback_props* feedback_props,
    enum stellar_type star_type);
void feedback_reset_feedback(struct spart* sp,
                             const struct feedback_props* feedback_props);
void feedback_first_init_spart(struct spart* sp,
                               const struct feedback_props* feedback_props);
void feedback_first_init_part(struct part *restrict p,
                              struct xpart *restrict xp);
void feedback_prepare_spart(struct spart* sp,
                            const struct feedback_props* feedback_props);
void feedback_prepare_feedback(struct spart* restrict sp,
                               const struct feedback_props* feedback_props,
                               const struct cosmology* cosmo,
                               const struct unit_system* us,
                               const struct phys_const* phys_const,
                               const double star_age_beg_step, const double dt,
                               const double time, const integertime_t ti_begin,
                               const int with_cosmology);
void feedback_struct_dump(const struct feedback_props* feedback, FILE* stream);
void feedback_struct_restore(struct feedback_props* feedback, FILE* stream);
void feedback_clean(struct feedback_props* feedback);

/**
 * @brief Writes the current model of feedback to the file
 *
 * @param feedback The #feedback_props.
 * @param h_grp The HDF5 group in which to write
 */
INLINE static void feedback_write_flavour(struct feedback_props* feedback,
                                          hid_t h_grp) {

  io_write_attribute_s(h_grp, "Feedback Model", "GEAR");
};

#endif /* SWIFT_FEEDBACK_GEAR_H */
