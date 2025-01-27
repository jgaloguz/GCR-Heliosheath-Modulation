#ifndef COMPUTE_KAPPA_PARALLEL_HH
#define COMPUTE_KAPPA_PARALLEL_HH

#include "compute_kappa_common.hh"
#include <iostream>
#include <iomanip>
#include <fstream>

//! Function to compute pitch-angle scattering coefficient
double Dmumu(double v, double mu, double B0, unsigned int isp = 0);

//! Function to compute the parallel diffusion coefficient
double KappaPara(double v, double B0, unsigned int isp = 0);

//! Function to "plot" parallel diffusion coefficient vs rigidity
void PlotKappaParaVsRigidity(double B0, unsigned int isp = 0);

//! Function to "plot" parallel diffusion coefficient along V2 trajectory
void PlotKappaParaVsRadius(double v, unsigned int isp = 0);

#endif