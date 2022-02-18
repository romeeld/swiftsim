/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_ACTIVE_H
#define SWIFT_RT_ACTIVE_H

/* Local includes. */
#include "active.h"

/**
 * @file rt_active.h
 * @brief This file contains the checks for whether radiative transfer
 * (injection) tasks should be activated for given cells, parts, or sparts.
 * The functions differ from the checks in `src/active.h` in that not only
 * time-step data is being checked, but more properties as well. So for
 * consistency, they get their own file. Finally, the functions are gathered
 * in one RT related file to concentrate all #ifdef macro shenanigans in a
 * single place as far as possible.
 */

/**
 * @brief Does a cell contain particles that should do RT this step?
 * This function is for a self-type interaction, where we need a cell
 * to have active hydro particles and star particles in any state.
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int rt_should_iact_cell(
    const struct cell *c, const struct engine *e) {

  return ((cell_is_active_stars(c, e) && (c->stars.count > 0)) &&
          (c->hydro.count > 0));
}

/**
 * @brief Does a cell contain particles that should do RT this step?
 * This function is for a pair-type interaction, where we take stars from
 * cell ci and hydro particles from cell cj.
 *
 * @param ci First #cell.
 * @param cj Second #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int rt_should_iact_cell_pair(
    const struct cell *ci, const struct cell *cj, const struct engine *e) {

  if (cj == NULL) return 0;

  return (cell_is_active_stars(ci, e) && (ci->stars.count > 0) &&
          (cj->hydro.count > 0));
}

/**
 * @brief Do we need to unskip this cell's RT (injection) related tasks?
 * For unskipping (recursively), don't check about the star's count here:
 * Pair-type interactions don't require stars in every cell. This way, we
 * can include the check for star count > 0 there on a pair by pair basis.
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell needs to activate tasks
 */
__attribute__((always_inline)) INLINE static int rt_should_do_unskip_cell(
    const struct cell *c, const struct engine *e) {
  /* whether it's hydro controlled or not, we need to check for hydro
   * activity at the top level. We also need to check for star activity
   * so we can activate the rt_in implicit tasks to catch dependencies
   * before the injection and not be overwritten by work in star density
   * ghosts. */
  return ((cell_is_active_hydro(c, e) && (c->hydro.count > 0)) ||
          cell_is_active_stars(c, e));
}

#endif /* SWIFT_RT_ACTIVE_H */
