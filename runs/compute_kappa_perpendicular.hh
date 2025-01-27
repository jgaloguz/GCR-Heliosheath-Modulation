#ifndef COMPUTE_KAPPA_PERPENDICULAR_HH
#define COMPUTE_KAPPA_PERPENDICULAR_HH

#include "compute_kappa_common.hh"
#include <iostream>
#include <iomanip>
#include <fstream>

//! Function to compute the perpendicular diffusion coefficient
double KappaPerp(double v, double B0, double Kpara, double Kperp);

//! Function to "plot" perpendicular diffusion coefficient vs rigidity
void PlotKappaPerpVsRigidity(double B0, double Kpara, double Kperp, unsigned int isp = 0);

//! Function to "plot" perpendicular diffusion coefficient along V2 trajectory
void PlotKappaPerpVsRadius(double v, double Kperp);

#endif