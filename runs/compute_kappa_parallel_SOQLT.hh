#ifndef COMPUTE_KAPPA_PARALLEL_SOQLT_HH
#define COMPUTE_KAPPA_PARALLEL_SOQLT_HH

#include "compute_kappa_parallel.hh"

//! Type of Gaussian width: 0 = 90-deg, large-time approx, 1 = 90-deg approx, 2 = no approx
#define GAUSSIAN_WIDTH_TYPE 2

//! Function to compute dB^2 / B0^2 based on PSD read from file.
void ComputeA2(double B0);

//! Function to compute the positive or negative resonance function
double Resonance(double v, double mu, double Omega, double k, double t, bool sign);

//! Function to integrate the positive + negative resonance functions over wavenumber using the trapezoid rule
double IntegralPSD_R(double v, double mu, double Omega, double t);

//! Function to compute the Gaussian width (squared)
double GaussianWidth2(double v, double mu, double B0, double Omega, double t);

//! Function to compute the derivative of Gaussian width (squared)
double dGaussianWidth2(double v, double mu, double Omega, double t);

//! Function to compute the characteristic function
double Characteristic(double v, double mu, double B0, double Omega, double k, double t);

//! Function to ompute the maximum time to integrate Gaussian resonance function
double ExpLimTime(double v, double mu, double B0, double Omega, double k);

//! Function to compute the analytic, improper (0 to infinity) time integral of the characteristic function for the 90-deg, large-time approximation
double IntegralCharacteristic90DegLargeTime(double v, double mu, double Omega, double k);

//! Function to integrate the characteristic function over time using the trapezoid rule
double IntegralCharacteristic(double v, double mu, double B0, double Omega, double k);

//! Function to integrate the PSD multiplied by the integral of the characteristic function over wavenumber using the trapezoid rule
double IntegralPSD_C(double v, double mu, double B0, double Omega);

#endif