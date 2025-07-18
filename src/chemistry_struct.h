/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_CHEMISTRY_STRUCT_H
#define SWIFT_CHEMISTRY_STRUCT_H

/**
 * @file src/chemistry_struct.h
 * @brief Branches between the different chemistry functions.
 */

/* Config parameters. */
#include <config.h>

/* Import the right chemistry definition */
#if defined(CHEMISTRY_NONE)
#include "./chemistry/none/chemistry_struct.h"
#elif defined(CHEMISTRY_GEAR)
#include "./chemistry/GEAR/chemistry_struct.h"
#elif defined(CHEMISTRY_GEAR_DIFFUSION)
#include "./chemistry/GEAR_DIFFUSION/chemistry_struct.h"
#elif defined(CHEMISTRY_AGORA)
#include "./chemistry/AGORA/chemistry_struct.h"
#elif defined(CHEMISTRY_QLA)
#include "./chemistry/QLA/chemistry_struct.h"
#elif defined(CHEMISTRY_EAGLE)
#include "./chemistry/EAGLE/chemistry_struct.h"
#elif defined(CHEMISTRY_SIMBA)
#include "./chemistry/SIMBA/chemistry_struct.h"
#elif defined(CHEMISTRY_KIARA)
#include "./chemistry/KIARA/chemistry_struct.h"
#else
#error "Invalid choice of chemistry function."
#endif

#endif /* SWIFT_CHEMISTRY_STRUCT_H */
