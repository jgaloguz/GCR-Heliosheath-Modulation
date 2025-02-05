#ifndef COMPUTE_KAPPA_PARALLEL_QLT_HH
#define COMPUTE_KAPPA_PARALLEL_QLT_HH

#include "compute_kappa_parallel.hh"

//! How to handle wavenumbers outside of data range: 0 = use edge values, 1 = extrapolate with linear fit
#define K_OUTSIDE_RANGE 1

//! Function to find linear fits to continue high and low wave-number ends of finite spectrum from data
void LinearFit(double* x, double* y, double& a, double& b);

//! Function to find linear fits to continue high and low wave-number ends of finite spectrum from data
void FitEndsPSD(void);

//! Function to interpolate discrete spectrum from data
double g_slab(double k);

#endif