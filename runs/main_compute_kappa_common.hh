#ifndef MAIN_COMPUTE_KAPPA_COMMON_HH
#define MAIN_COMPUTE_KAPPA_COMMON_HH

#include "common/physics.hh"
#include <iostream>
#include <iomanip>
#include <fstream>

//! Number of (k,PSD) points in the file to be read
const int Nk = 40;

//! Number of (k,PSD) points to use for fits of ends of spectrum (for QLT)
const int Nkfit = Nk / 10;

//! Number of Bmag points in the file to be read
const int NB = 365;

//! Lower bound to plot spectrum
const double k1 = 1.0e-12 * unit_length_fluid;

//! Upper bound to plot spectrum
const double k2 = 1.0e-5 * unit_length_fluid;

//! Wavenumber bin width in logarithmic space
const double dlnk = (log(k2) - log(k1)) / (Nk-1);

//! Number of pitch-angle points to use for plotting and integration
const int Nmu = 100;

//! Minimum pitch-angle to use for plotting and integration
const double mu0 = 0.001;

//! Pitch-angle bin width in linear space
const double dmu = (1.0 - mu0) / (Nmu-1);

//! Number of kinetic energy values to use in computations
const int NE = 10;

//! Minimum kinetic energy
const double E1 = 100.0 * SPC_CONST_CGSM_KILO_ELECTRON_VOLT / unit_energy_particle;

//! Maximum kinetic energy
const double E2 = 10.0 * SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle;

//! Kinetic energy bin width in logarithmic space 
const double dlnE = (log(E2) - log(E1)) / (NE-1);

//! Number of points per period to resolve in the integral of resonance function
const int Ntpc = 25;

//! Threshold value of resonance envelope function to limit improper integral (upper bound is infinity)
const double exp_lim = -log(0.0001);

//! Arrays with k and PSD values
extern double k_vals[Nk], PSD_vals[Nk];

//! Array with magnitude of radius and B along V2 trajectory
extern double R_V2[NB], Bmag_V2[NB];

//! Function to read turbulence spectrum from file
void ReadPSD(std::string PSD_filename);

//! Function to read V2 magnetic field measurements from file
void ReadBmagV2(std::string V2_filename);

//! Function to compute pitch-angle scattering coefficient
double Dmumu(double v, double mu, double B0, unsigned int isp = 0);

//! Function to compute parallel diffusion coefficient
double KappaPara(double v, double B0, unsigned int isp = 0);

//! Function to "plot" parallel diffusion coefficient vs rigidity
void KappaParaVsRigidity(double B0, unsigned int isp = 0);

//! Function to "plot" parallel diffusion coefficient along V2 trajectory
void KappaParaVsRadius(double v, unsigned int isp = 0);

#endif