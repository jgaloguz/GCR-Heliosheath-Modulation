#ifndef COMPUTE_KAPPA_PERPENDICULAR_UNLT_HH
#define COMPUTE_KAPPA_PERPENDICULAR_UNLT_HH

#include "compute_kappa_perpendicular.hh"

//! Function to compute the next fixed-point iteration of the perpendicular diffusion coefficient
double FixedPointUNLT(double v, double B0, double Kpara, double Kperp);

#endif