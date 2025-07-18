/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2014 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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
#ifndef SWIFT_GIZMO_MFV_HYDRO_PART_H
#define SWIFT_GIZMO_MFV_HYDRO_PART_H

/* Data of a single particle. */
struct part {

  /*! Particle ID. */
  long long id;

  /*! Associated gravitas. */
  struct gpart *gpart;

  /*! Particle position. */
  double x[3];

  /*! Particle predicted velocity. */
  float v[3];

  /*! Particle velocity for drift */
  float v_full[3];

  /*! Particle acceleration. */
  float a_hydro[3];

  /*! Particle smoothing length. */
  float h;

  /*! Density. */
  float rho;

  /*! Fluid velocity. */
  float fluid_v[3];

  /*! Pressure. */
  float P;

  /*! Gradients of the primitive variables. */
  struct {

    /* Density gradients. */
    float rho[3];

    /* Fluid velocity gradients. */
    float v[3][3];

    /* Pressure gradients. */
    float P[3];

  } gradients;

  /*! Quantities needed by the slope limiter. */
  struct {

    /* Extreme values of the density among the neighbours. */
    float rho[2];

    /* Extreme values of the fluid velocity among the neighbours. */
    float v[3][2];

    /* Extreme values of the pressure among the neighbours. */
    float P[2];

    /* Maximal distance to all neighbouring faces. */
    float maxr;

  } limiter;

  /*! The conserved hydrodynamical variables. */
  struct {

    /* Fluid mass */
    float mass;

    /* Fluid momentum. */
    float momentum[3];

    /* Fluid thermal energy (not per unit mass!). */
    float energy;

  } conserved;

  /*! Fluxes. */
  struct {

    /* Mass flux. */
    float mass;

    /* Momentum flux. */
    float momentum[3];

    /* Energy flux. */
    float energy;

    /* Particle time step. Used to compute time-integrated fluxes. */
    float dt;

  } flux;

  /*! Geometrical quantities used for hydro. */
  struct fvpm_geometry_struct geometry;

  /*! Variables used for timestep calculation. */
  struct {

    /* Maximum signal velocity among all the neighbours of the particle. The
     * signal velocity encodes information about the relative fluid velocities
     * AND particle velocities of the neighbour and this particle, as well as
     * the sound speed of both particles. */
    float vmax;

  } timestepvars;

  /*! Quantities used during the volume (=density) loop. */
  struct {

    /* Derivative of particle number density. */
    float wcount_dh;

    /* Particle number density. */
    float wcount;

  } density;

  /*! Quantities used during the force loop. */
  struct {

    /* Needed to drift the primitive variables. */
    float h_dt;

  } force;

  /*! Specific stuff for the gravity-hydro coupling. */
  struct {

    /* Current value of the mass flux vector. */
    float mflux[3];

  } gravity;

  /*! Flag for decoupling from the hydrodynamics/feedback routines */
  unsigned char decoupled;

  /*! Flag to indicate that the decoupling task will run */
  unsigned char to_be_decoupled;
  
  /*! Flag to indicate that the recoupling task will run */
  unsigned char to_be_recoupled;
  
  /*! Chemistry information */
  struct chemistry_part_data chemistry_data;

  /*! Cooling information */
  struct cooling_part_data cooling_data;

  /*! Additional data used by the feedback */
  struct feedback_part_data feedback_data;

#ifdef WITH_FOF_GALAXIES
  /*! Additional data used by the FoF */
  struct galaxy_data galaxy_data;
#endif

  /*! Black holes information (e.g. swallowing ID) */
  struct black_holes_part_data black_holes_data;

  /*! Sink information (e.g. swallowing ID) */
  struct sink_part_data sink_data;

  /*! Additional Radiative Transfer Data */
  struct rt_part_data rt_data;

  /*! RT sub-cycling time stepping data */
  struct rt_timestepping_data rt_time_data;

  /*! Time-step length */
  timebin_t time_bin;

  /*! Tree-depth at which size / 2 <= h * gamma < size */
  char depth_h;

  /*! Time-step limiter information */
  struct timestep_limiter_data limiter_data;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_GIZMO_MFV_HYDRO_PART_H */
