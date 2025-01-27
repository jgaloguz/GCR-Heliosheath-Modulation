#include "compute_kappa_perpendicular.hh"

using namespace Spectrum;

//! Maximum number of fixed point iterations
const int FP_maxiters = 100;

//! Convergence threshold for fixed-point iterations
const double FP_threshold = 1.0e-8;

//! Phenomenological constant
const double a_UNLT = 0.4;

/*!
\brief Compute the next fixed-point iteration of the perpendicular diffusion coefficient
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v       particle speed
\param[in] B0      mean magnetic field
\param[in] Kpara   parallel diffusion coefficient
\param[in] Kperp   perpendicular diffusion coefficient from previous fixed-point iteration
\return next fixed-point iteration of the perpendicular diffusion coefficient according to UNLT
*/
double FixedPointUNLT(double v, double B0, double Kpara, double Kperp)
{
   int i;
   double S = 0.5 * (k_vals[1] - k_vals[0]) * ( PSD_vals[1] / (4.0 * Kperp * Sqr(k_vals[1]) + Sqr(v) / Kpara)
                                              + PSD_vals[0] / (4.0 * Kperp * Sqr(k_vals[0]) + Sqr(v) / Kpara) );
// Integrate PSD using trapezoid rule
   for (i = 1; i < Nk-1; i++) S += 0.5 * ( PSD_vals[i+1] / (4.0 * Kperp * Sqr(k_vals[i+1]) + Sqr(v) / Kpara)
                                         + PSD_vals[i] / (4.0 * Kperp * Sqr(k_vals[i]) + Sqr(v) / Kpara) );

   return M_PI * Sqr(a_UNLT * v / B0) * S;
};

/*!
\brief Compute the perpendicular diffusion coefficient
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v       particle speed
\param[in] B0      mean magnetic field
\param[in] Kpara   parallel diffusion coefficient
\param[in] Kperp   perpendicular diffusion coefficient (to start iterations if necessary)
\return perpendicular diffusion coefficient according to UNLT
*/
double KappaPerp(double v, double B0, double Kpara, double Kperp)
{
   int i;
   double Kperp_old;

// Perform fixed-point iterations until convergence
   for (i = 0; i < FP_maxiters; i++) {
      Kperp_old = Kperp;
      Kperp = FixedPointUNLT(v, B0, Kpara, Kperp_old);
      if (fabs(Kperp - Kperp_old) / (fabs(Kperp) + fabs(Kperp_old)) < FP_threshold) break;
   };
   if (i == FP_maxiters) std::cerr << "WARNING: Maximum number of FP iterations reached." << std::endl;

   return Kperp;
};

int main(int argc, char** argv)
{
// Specie
   int specie = Specie::alpha_particle;
   double vel = Vel(Mom(1.0 * SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle, specie), specie);

// Read PSD from file
   ReadPSD("data/k_spectra_perpendicular_SHS.dat", M_2PI);
   ReadBmagV2("data/V2_Bmag_2013_303_2014_365.dat");
   ReadKappaParaVsRadius("../results/kappa_SOQLT_He.dat");

// Compute kappa_perp vs R_V2
// Initializing iterations with Kperp = 0 will yield the FLRW limit after the first iteration
   PlotKappaPerpVsRadius(vel, 0.0);

   return 0;
};