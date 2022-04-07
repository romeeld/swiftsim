/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c)    2022 Tsang Keung Chan (chantsangkeung@gmail.com)
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
#ifndef SWIFT_RT_SPHM1RT_PROTO_H
#define SWIFT_RT_SPHM1RT_PROTO_H
/**
 * @file src/rt/SPHM1RT/rt_proto.h
 * @brief Function proto-type for the SPHM1RT radiative transfer scheme
 * thermochemistry related functions.
 */


#include "rt_properties.h"
#include "rt_struct.h"


/* Local includes. */
#include <cvode/cvode.h>
#include <math.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <sys/types.h>
#include <time.h>

static int frateeq(realtype t, N_Vector y, N_Vector ydot, void *user_data);


#endif /* SWIFT_RT_SPHM1RT_PROTO_H */