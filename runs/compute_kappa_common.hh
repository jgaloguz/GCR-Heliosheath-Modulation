#ifndef COMPUTE_KAPPA_COMMON_HH
#define COMPUTE_KAPPA_COMMON_HH

#include "common/physics.hh"
#include <iostream>
#include <iomanip>
#include <fstream>

//! Number of (k,PSD) points in the file to be read
extern int Nk;

//! Number of (k,PSD) points to use for fits of ends of spectrum (for QLT)
const int Nkfit = fmax(Nk / 10, 2);

//! Number of Bmag points in the file to be read
extern int NB;

//! Number of pitch-angle points to use for plotting and integration
const int Nmu = 100;

//! Minimum pitch-angle to use for plotting and integration
const double mu0 = 0.001;

//! Pitch-angle bin width in linear space
const double dmu = (1.0 - mu0) / (Nmu-1);

//! Number of kinetic energy values to use in computations
const int NE = 50;

//! Minimum kinetic energy
const double E1 = 100.0 * SPC_CONST_CGSM_KILO_ELECTRON_VOLT / unit_energy_particle;

//! Maximum kinetic energy
const double E2 = 100.0 * SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle;

//! Kinetic energy bin width in logarithmic space 
const double dlnE = (log(E2) - log(E1)) / (NE-1);

//! Number of points per period to resolve in the integral of resonance function
const int Ntpc = 25;

//! Threshold value of resonance envelope function to limit improper integral (upper bound is infinity)
const double exp_lim = -log(0.0001);

//! Arrays with k and PSD values
extern double *k_vals, *PSD_vals;

//! Array with magnitude of radius and B along V2 trajectory
extern double *R_V2, *Bmag_V2;

//! Array with parallel diffusion coefficients
extern double *kappa_para_rig, *kappa_para_V2;

//! Function to read turbulence spectrum from file
void ReadPSD(std::string PSD_fp, int Nk_in, double norm_const);

//! Function to read V2 magnetic field measurements from file
void ReadBmagV2(std::string V2_fp, int NB_in);

//! Function to de-allocate memory
void FreeMemory();

#endif