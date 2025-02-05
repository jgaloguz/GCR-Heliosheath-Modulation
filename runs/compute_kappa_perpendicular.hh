#ifndef COMPUTE_KAPPA_PERPENDICULAR_HH
#define COMPUTE_KAPPA_PERPENDICULAR_HH

#include "compute_kappa_common.hh"

//! Function to read the computed values of parallel diffusion coefficient vs rigidity from file
void ReadKappaParaVsRigidity(std::string k_para_fp);

//! Function to read the computed values of parallel diffusion coefficient vs radius from file
void ReadKappaParaVsRadius(std::string k_para_fp);

//! Function to compute the perpendicular diffusion coefficient
double KappaPerp(double v, double B0, double Kpara, double Kperp);

//! Function to "plot" perpendicular diffusion coefficient vs rigidity
void PlotKappaPerpVsRigidity(std::string out_fp, double B0, double Kperp, unsigned int isp = 0);

//! Function to "plot" perpendicular diffusion coefficient along V2 trajectory
void PlotKappaPerpVsRadius(std::string out_fp, double v, double Kperp, unsigned int isp = 0);

#endif