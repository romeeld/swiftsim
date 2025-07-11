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
#ifndef SWIFT_NONE_FEEDBACK_IACT_H
#define SWIFT_NONE_FEEDBACK_IACT_H

/* Local includes */
#include "feedback.h"

/**
 * @brief Compute the mean DM velocity around a star. (non-symmetric).
 *
 * @param si First sparticle.
 * @param gj Second particle (not updated).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_dm_vel_sum(struct spart *si, const struct gpart *gj,
                                       int *dm_ngb_N,
                                       float dm_mean_velocity[3]) {}

/**
 * @brief Compute the DM velocity dispersion around a star. (non-symmetric).
 *
 * @param si First sparticle.
 * @param gj Second particle.
 * @param fb_props Properties of the feedback scheme.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_dm_vel_disp(struct spart *si,
                                        const struct gpart *gj,
                                        const float dm_mean_velocity[3]) {}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * Nothing to do here for the no-feedback model.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param xp Extra particle data (not updated).
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time value
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_density(const float r2, const float dx[3],
                                    const float hi, const float hj,
                                    struct spart *si, const struct part *pj,
                                    const struct xpart *xp,
                                    const struct cosmology *cosmo,
                                    const struct feedback_props *fb_props,
                                    const integertime_t ti_current) {}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 * Used for updating properties of gas particles neighbouring a star particle
 *
 * Nothing to do here for the no-feedback model.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xp Extra particle data
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(
    const float r2, const float dx[3], const float hi, const float hj,
    const struct spart *si, struct part *pj, struct xpart *xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props, const integertime_t ti_current) {}

#endif /* SWIFT_NONE_FEEDBACK_IACT_H */
