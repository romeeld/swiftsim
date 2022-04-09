/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Tsang Keung Chan (chantsangkeung@gmail.com)
 *               2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_STRUCT_SPHM1RT_H
#define SWIFT_RT_STRUCT_SPHM1RT_H

/**
 * @file src/rt/SPHM1RT/rt_struct.h
 * @brief Main header file for no radiative transfer struct.
 * SPHM1RT method described in Chan+21: 2102.08404
 */

/**
 * @brief The individual elements traced in the SPHM1RT model.
 */
enum rt_chemistry_element {
  rt_chemistry_element_H = 0,
  rt_chemistry_element_He,
  rt_chemistry_element_count
};

/**
 * @brief The individual species traced.
 */
enum rt_cooling_species {
  rt_sp_elec = 0, /* 0 */
  rt_sp_HI,       /* 1 */
  rt_sp_HII,      /* 2 */
  rt_sp_HeI,      /* 3 */
  rt_sp_HeII,     /* 4 */
  rt_sp_HeIII,    /* 5 */
  rt_species_count
};

/* Additional RT data in hydro particle struct */
struct rt_part_data {

  /*! time step of the gas particle */
  float dt;

  /* conserved state vector in comoving units */
  /* but comoving and physical urad and frad are the same in our convention here
   */
  /* urad: radiation energy per mass */
  /* frad: radiation flux per gas density */
  /* (they are conserved in the sense of energy/mass; assuming mass is equal) */
  struct {
    float urad;
    float frad[3];
  } conserved[RT_NGROUPS];

  /* rate of change of the conserved state vector */
  struct {
    float urad;
    float frad[3];
  } dconserved_dt[RT_NGROUPS];

  /* Store viscosity information in a separate struct. */
  struct {

    /*! Particle radiation flux divergence */
    float divf;

    /*! Particle radiation flux divergence from previous step */
    float divf_previous_step;

    /* parameter to control dissipation */
    float alpha;

  } viscosity[RT_NGROUPS];

  /* Store artificial diffusion information in a separate struct. */
  struct {

    /*! gradient of radiation energy density per gas density */
    float graduradc[3];

    /* parameter to control dissipation */
    float alpha;

  } diffusion[RT_NGROUPS];

  /* Store radiation parameter in a separate struct. */
  struct {

    /*! initial mean opacity */
    float chi[RT_NGROUPS];

    /*! reduced speed of light */
    float cred;

    /*! switch to reinject radiation from injection radius */
    int reinject;

  } params;

  /* Store hydro information in a separate struct. */
  struct {

    /*! "Grad h" term */
    float f;

  } force;

  struct {

    /*! Fraction of the particle mass in a given element */
    float metal_mass_fraction[rt_chemistry_element_count];

    /*! Fraction of the particle mass in *all* metals */
    float metal_mass_fraction_total;

    /*! abundances of species i, i.e. n_i/nH */
    /* note that we use hydrogen density in the denominator */
    float abundances[rt_species_count];

  } tchem;
};

/* Additional RT data in star particle struct */
struct rt_spart_data {

  /* Stellar energy emission that will be injected in to gas.
   * Total energy, not density, not rate! */
  float emission_this_step[RT_NGROUPS];

  /*! normalisation factor used for the radiation injection */
  float injection_weight;

  /*! radiation energy within injection radius */
  float emission_reinject[RT_NGROUPS];
};

struct UserData {

  // void *cvode_mem;           /*!< Pointer to the CVODE memory. */

  /* switch for on the spot approximation */
  int onthespot;

  /* switch for gas cooling */
  int coolingon;

  /* switch for not changing photon density */
  int fixphotondensity;

  /* 1: to use the input parameters; 0: calculate with temperature. */
  /* (H coefficient only; no heating or cooling) */
  int useparams;

  /*! Fraction of the particle mass in a given element */
  double metal_mass_fraction[rt_chemistry_element_count];

  double m_H_cgs;

  double k_B_cgs;

  double cred_cgs;

  double rho_cgs;

  double n_H_cgs;

  double ngamma_cgs[3];

  double u_cgs;

  /*! abundances of species i, i.e. n_i/nH */
  /* note that we use hydrogen density in the denominator */
  double abundances[rt_species_count];

  double u_min_cgs;

  int aindex[3];

  /* only use when useparam = 1 */
  /*! The case A recombination coefficient for hydrogen (cgs) */
  double alphaA_cgs_H;

  /*! The case B recombination coefficient for hydrogen (cgs) */
  double alphaB_cgs_H;

  /*! The collisional ionization coefficient for hydrogen (cgs) */
  double beta_cgs_H;

  /*! The cross section of ionizing photons for hydrogen (cgs) */
  double sigma_cross_cgs_H[3];
};

#endif /* SWIFT_RT_STRUCT_SPHM1RT_H */
