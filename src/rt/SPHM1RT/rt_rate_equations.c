/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Tsang Keung Chan (chantsangkeung@gmail.com)
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


/**
 * @brief Defines the right-hand side function.
 *
 * Defines the system of differential equations that make
 * up the right-hand side function, which will be integrated
 * by CVode.
 *
 * @param t Current time.
 * @param y Vector containing the variables to be integrated.
 * @param ydot Vector containing the time derivatives of the variables.
 * @param user_data The #UserData struct containing the input data.
 */
int f(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
  struct UserData *data;

  data = (struct UserData *)user_data;

  /* First, loop through the enum types of all
   * non-eq species. If they are included in
   * the network then their abundance is in
   * the vector y. */
  int icount = 0; /* We use this to keep track of where we are in the vector y */
  int aindex[3];
  for (int i = 0; i < 3; i++) {
    aindex[i] = data->aindex[i];
  }
  for (int i = 0; i < 3; i++) {  
    data->abundances[aindex[i]] = (double)NV_Ith_S(y, icount);
    icount += 1;
  }

  /* Update the species not in the network */
  double finish_abundances[rt_species_count];
  enforce_constraint_equations(data->abundances, data->metal_mass_fraction, finish_abundances);
  for (int j = 0; j < rt_species_count; j++) {
    data->abundances[j] = finish_abundances[j];
  }

  /* If Thermal Evolution is switched on, the element in the
   * vector y is the internal energy (per unit volume). Use this
   * to update the temperature, and also the rates that depend on T */
  double u_cgs; 
  if (data->coolingon==1) {
    u_cgs = (double)NV_Ith_S(y,icount);
    icount += 1;
  } else {
    u_cgs = data->u_cgs;
  }

  /* the final element in the
   * vector y is the photon density.
   * */
  double ngamma_cgs[3]; 
  if (data->fixphotondensity==0) {
    for (int i = 0; i < 3; i++) {  
      ngamma_cgs[i] = (double)NV_Ith_S(y, icount);
      icount += 1;
    }
  } else {
    for (int i = 0; i < 3; i++) {  
      ngamma_cgs[i] = data->ngamma_cgs[i];
    }
  }



  double log_u_cgs = log10(u_cgs);
  double log_T_cgs = convert_u_to_temp(data->k_B_cgs, data->m_H_cgs, data->metal_mass_fraction[rt_chemistry_element_H], log_u_cgs, data->abundances);

  // Update rates
  double alphalist[rt_species_count], betalist[rt_species_count], Gammalist[rt_species_count], sigmalist[3][3], epsilonlist[3][3];

  compute_rate_coefficients(log_T_cgs, data->onthespot, alphalist, betalist, Gammalist, sigmalist, epsilonlist, aindex);

  if (data->useparams == 1) {
    betalist[rt_sp_elec] = 0.0;
    betalist[rt_sp_HI] = data->beta_cgs_H;
    betalist[rt_sp_HII] = 0.0;
    betalist[rt_sp_HeI] = 0.0;
    betalist[rt_sp_HeII] = 0.0;
    betalist[rt_sp_HeIII] = 0.0;
    alphalist[rt_sp_elec] = 0.0;
    alphalist[rt_sp_HI] = 0.0;
    alphalist[rt_sp_HeI] = 0.0; 
    if (data->onthespot==1) {
      alphalist[rt_sp_HII] = data->alphaB_cgs_H;
      alphalist[rt_sp_HeII] = 0.0;
      alphalist[rt_sp_HeIII] = 0.0;
    } else {
      alphalist[rt_sp_HII] = data->alphaA_cgs_H;
      alphalist[rt_sp_HeII] = 0.0;
      alphalist[rt_sp_HeIII] = 0.0;
    }

    sigmalist[0][0] = data->sigma_cross_cgs_H[0];
    sigmalist[1][0] = data->sigma_cross_cgs_H[1];
    sigmalist[2][0] = data->sigma_cross_cgs_H[2];
    sigmalist[0][1] = 0.0;
    sigmalist[1][1] = 0.0;
    sigmalist[2][1] = 0.0;
    sigmalist[0][2] = 0.0;
    sigmalist[1][2] = 0.0;
    sigmalist[2][2] = 0.0;    
  }

  // Compute creation and destruction rates
  double absorption_rate[3], chemistry_rates[rt_species_count];

  compute_radiation_rate(data->n_H_cgs, data->cred_cgs, data->abundances, ngamma_cgs, sigmalist, aindex, absorption_rate);

  compute_chemistry_rate(data->n_H_cgs, data->cred_cgs, data->abundances, ngamma_cgs, alphalist, betalist, sigmalist, aindex, chemistry_rates);

  double Lambda_net_cgs;
  Lambda_net_cgs = compute_cooling_rate(data->n_H_cgs, data->cred_cgs, data->abundances, ngamma_cgs, Gammalist, sigmalist, epsilonlist, aindex);

  int jcount = 0;
  /* Now set the output ydot vector for the chemical abundances */
  for (int i = 0; i < 3; i++) {  
    NV_Ith_S(ydot, jcount) = (realtype)(chemistry_rates[aindex[i]] / data->n_H_cgs);
    jcount += 1;
  }
  /* Now set the output ydot vector for the internal energy */
  if (data->coolingon == 1) {
    NV_Ith_S(ydot, jcount) = (realtype)(Lambda_net_cgs / data->rho_cgs);
    jcount += 1;
  }

  /* Now set the output ydot vector for the radiation density */
  if (data->fixphotondensity==0) {
    for (int i = 0; i < 3; i++) { 
      NV_Ith_S(ydot, jcount) = (realtype) (-absorption_rate[i]);
      jcount += 1;
    }
  }
  return(0);
}