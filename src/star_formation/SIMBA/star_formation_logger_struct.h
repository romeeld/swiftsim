/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_SIMBA_STAR_FORMATION_LOGGER_STRUCT_H
#define SWIFT_SIMBA_STAR_FORMATION_LOGGER_STRUCT_H

/* Starformation history struct */
struct star_formation_history {
  /*! Total new stellar mass */
  float new_stellar_mass;

  /*! SFR of all particles */
  float SFR_inactive;

  /*! SFR of active particles */
  float SFR_active;

  /*! SFR*dt of active particles */
  float SFRdt_active;
};

/* Starformation history struct for the engine.
 Allows to integrate in time some values.
 Nothing to do in SIMBA => copy of star_formation_history */
struct star_formation_history_accumulator {
  /*! Total new stellar mass */
  float new_stellar_mass;

  /*! SFR of all particles */
  float SFR_inactive;

  /*! SFR of active particles */
  float SFR_active;

  /*! SFR*dt of active particles */
  float SFRdt_active;
};

#endif /* SWIFT_SIMBA_STAR_FORMATION_LOGGER_STRUCT_H */
