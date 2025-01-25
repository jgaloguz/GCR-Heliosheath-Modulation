#include "main_compute_kappa_common.hh"

using namespace Spectrum;

//! How to handle wavenumbers outside of data range: 0 = use edge values, 1 = extrapolate with linear fit
#define K_OUTSIDE_RANGE 1

//! Parameters for linear fits to the left and right of spectrum
double a_l, b_l, a_r, b_r;

/*!
\brief Find linear fits to continue high and low wave-number ends of finite spectrum from data
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] x array of independent values (abscissa)
\param[in] y array of dependent values (ordinate)
\param[in] a constant coefficient
\param[in] b linear coefficient
*/
void LinearFit(double* x, double* y, double& a, double& b)
{
   int i;
   double xm = 0.0, ym = 0.0, denom = 0.0;

// Find averages
   for(i = 0; i < Nkfit; i++) {
      xm += x[i];
      ym += y[i];
   };
   xm /= Nkfit;
   ym /= Nkfit;

// Find optimal parameters
   for(i = 0; i < Nkfit; i++) denom += Sqr(x[i] - xm);
   for(i = 0; i < Nkfit; i++) b += (x[i] - xm) * (y[i] - ym);
   b /= denom;
   a = ym - b * xm;
};

/*!
\brief Find linear fits to continue high and low wave-number ends of finite spectrum from data
\author Juan G Alonso Guzman
\date 01/24/2025
*/
void FitEndsPSD(void)
{
   int i;
   double lnk_vals[Nkfit], lnPSD_vals[Nkfit];

// Fit left end of PSD with constant power slope
   for(i = 0; i < Nkfit; i++) {
      lnk_vals[i] = log(k_vals[i]);
      lnPSD_vals[i] = log(PSD_vals[i]);
   };
   LinearFit(lnk_vals, lnPSD_vals, a_l, b_l);

// Fit right end of PSD with constant power slope
   for(i = 0; i < Nkfit; i++) {
      lnk_vals[i] = log(k_vals[Nk-Nkfit+i]);
      lnPSD_vals[i] = log(PSD_vals[Nk-Nkfit+i]);
   };
   LinearFit(lnk_vals, lnPSD_vals, a_r, b_r);
};

/*!
\brief Interpolate discrete spectrum from data
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] k wave-number
\return magnitude of turbulence spectrum at given wave-number
*/
double g_slab(double k)
{
   int k_idx;
// Locate interval in PSD array containing value "k"
   k_idx = LocateInArray(0, Nk-1, k_vals, k, true);
// Interpolate PSD
#if K_OUTSIDE_RANGE == 0
   if(k_idx == -1) return PSD_vals[0];
   else if(k_idx == Nk-1) return PSD_vals[Nk-1];
#elif K_OUTSIDE_RANGE == 1
   if(k_idx == -1) return exp(a_l + b_l*log(k));
   else if(k_idx == Nk-1) return exp(a_r + b_r*log(k));
#endif
   else return PSD_vals[k_idx] + (k - k_vals[k_idx]) * (PSD_vals[k_idx+1] - PSD_vals[k_idx]) / (k_vals[k_idx+1] - k_vals[k_idx]);
};

/*!
\brief Compute pitch-angle scattering coefficient according to QLT
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v   particle speed
\param[in] mu  partickle pitch-angle
\param[in] B0  mean magnetic field
\param[in] isp specie index
\return pitch-angle scattering coefficient
*/
double Dmumu(double v, double mu, double B0, unsigned int isp)
{
   double Omega = CyclotronFrequency(v, B0, isp);
   double k = Omega / (v * mu);
   return 2.0 * Sqr(M_PI) * Omega * (1.0 - Sqr(mu)) * k * g_slab(k) / Sqr(B0);
};

int main(int argc, char** argv)
{
// Specie
   int specie = Specie::alpha_particle;
   double vel = Vel(Mom(1.0 * SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle, specie), specie);

// Read PSD and data from files
   ReadPSD("data/k_spectra_parallel_SHS.dat");
   FitEndsPSD();
   ReadBmagV2("data/V2_Bmag_2013_303_2014_365.dat");

// Compute kappa_parallel vs R_V2
   KappaParaVsRadius(vel);

   return 0;
};