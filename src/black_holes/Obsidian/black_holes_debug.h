/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_BLACK_HOLES_OBSIDIAN_DEBUG_H
#define SWIFT_BLACK_HOLES_OBSIDIAN_DEBUG_H

__attribute__((always_inline)) INLINE static void black_holes_debug_particle(
    const struct part* p, const struct xpart* xp) {

  warning("[PID%lld] black_holes_part_data:", p->id);
  warning("[PID%lld] swallow_id = %lld, potential = %.3e", p->id,
          p->black_holes_data.swallow_id, p->black_holes_data.potential);
}

#endif /* SWIFT_BLACK_HOLES_OBSIDIAN_DEBUG_H */
