#include "main_compute_kappa_common.hh"

using namespace Spectrum;

//! Type of Gaussian width: 0 = 90-deg, large-time approx, 1 = 90-deg approx, 2 = no approx
#define GAUSSIAN_WIDTH_TYPE 2

//! Ratios dB / B0 and dB^2 / B0^2
double A, A2;

/*!
\brief Compute dB^2 / B0^2 based on PSD read from file.
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] B0 mean magnetic field
*/
void ComputeA2(double B0)
{
   int i;
// Integrate PSD using the trapezoidal rule to find delta B^2
   A2 = 0.5 * (k_vals[1] - k_vals[0]) * (PSD_vals[1] + PSD_vals[0]);
   for(i = 1; i < Nk-1; i++) A2 += 0.5 * (k_vals[i+1] - k_vals[i]) * (PSD_vals[i+1] + PSD_vals[i]);
   A2 *= M_8PI;
   A2 /= Sqr(B0);
   A = sqrt(A2);
};

/*!
\brief Compute the positive or negative resonance function
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v     particle speed
\param[in] mu    partickle pitch-angle
\param[in] Omega cyclotron frequency
\param[in] k     wave-number
\param[in] t     time
\param[in] sign  whether to compute positive (True) or negative (False) resonance
\return positive (sign == True) or negative (sign == False) resonance function
*/
double Resonance(double v, double mu, double Omega, double k, double t, bool sign)
{
   double beta = k * v * mu + (sign ? 1.0 : -1.0) * Omega;
   double betat = beta * t;
// Theoretical limit of resonance function as beta -> 0 is t^4 / 8
   if(betat <= 1.0e-3) return 0.125 * Sqr(Sqr(t));
// Empirically, the resonance formula is reliable for beta * t > 10^-3. Below this threshold, round-off errors give garbage because of finite precision.
   else return (1.0 - cos(betat) - betat * sin(betat) + 0.5 * Sqr(betat)) / Sqr(Sqr(beta));
};

/*!
\brief Integrate the positive + negative resonance functions over wavenumber using the trapezoid rule
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v     particle speed
\param[in] mu    partickle pitch-angle
\param[in] Omega cyclotron frequency
\param[in] t     time
\return integral of the positive + negative resonance functions over wavenumber
*/
double IntegralPSD_R(double v, double mu, double Omega, double t)
{
   int i;
   double Res[Nk];

// Pre-compute integral of resonance function
   for(i = 0; i < Nk; i++) Res[i] = Resonance(v, mu, Omega, k_vals[i], t, true) + Resonance(v, mu, Omega, k_vals[i], t, false);

   double S = 0.5 * (k_vals[1] - k_vals[0]) * (PSD_vals[1] * Res[1] + PSD_vals[0] * Res[0]);
// Integrate PSD using trapezoid rule
   for(i = 1; i < Nk-1; i++) S += 0.5 * (k_vals[i+1] - k_vals[i]) * (PSD_vals[i+1] * Res[i+1] + PSD_vals[i] * Res[i]);
   return S;
};

/*!
\brief Compute the Gaussian width (squared)
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v     particle speed
\param[in] mu    partickle pitch-angle
\param[in] B0    mean magnetic field
\param[in] Omega cyclotron frequency
\param[in] t     time
\return Gaussian width (squared)
*/
double GaussianWidth2(double v, double mu, double B0, double Omega, double t)
{
#if GAUSSIAN_WIDTH_TYPE == 0
   return 0.5 * A2 * Sqr(v * t);
#elif GAUSSIAN_WIDTH_TYPE == 1
   double Omegat = Omega * t;
   return A2 * Sqr(v / Omega) * (1.0 - cos(Omegat) - Omegat * sin(Omegat) + 0.5 * Sqr(Omegat));
#elif GAUSSIAN_WIDTH_TYPE == 2
   return M_4PI * Sqr(Omega * v / B0) * (1.0 - Sqr(mu)) * IntegralPSD_R(v, mu, Omega, t);
#endif
};

/*!
\brief Compute the derivative of Gaussian width (squared)
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v     particle speed
\param[in] mu    partickle pitch-angle
\param[in] Omega cyclotron frequency
\param[in] t     time
\return derivative of Gaussian width (squared)
*/
double dGaussianWidth2(double v, double mu, double Omega, double t)
{
#if GAUSSIAN_WIDTH_TYPE == 0
   return A2 * Sqr(v) * t;
#elif GAUSSIAN_WIDTH_TYPE == 1
   double Omegat = Omega * t;
   return A2 * Sqr(v) / Omega * (sin(Omegat) - sin(Omegat) - Omegat * cos(Omegat) + Omega * t);
#elif GAUSSIAN_WIDTH_TYPE == 2
   return 0.0; //TODO
#endif
};

/*!
\brief Compute the characteristic function
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v     particle speed
\param[in] mu    partickle pitch-angle
\param[in] B0    mean magnetic field
\param[in] Omega cyclotron frequency
\param[in] k     wave-number
\param[in] t     time
\return characteristic function
*/
double Characteristic(double v, double mu, double B0, double Omega, double k, double t)
{
   return (cos((k * v * mu + Omega) * t) + cos((k * v * mu - Omega) * t)) * exp(-0.5 * Sqr(k) * GaussianWidth2(v, mu, B0, Omega, t));
};

/*!
\brief Compute the maximum time to integrate Gaussian resonance function
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v     particle speed
\param[in] mu    partickle pitch-angle
\param[in] B0    mean magnetic field
\param[in] Omega cyclotron frequency
\param[in] k     wave-number
\return maximum time to integrate Gaussian resonance function
*/
double ExpLimTime(double v, double mu, double B0, double Omega, double k)
{
   int i = 0, imax = 100;
   double t0, t1 = 2.0 * sqrt(exp_lim) / (v * k * A);
   double C = -2.0 * exp_lim / Sqr(k);
   double eps = 1.0, tol = 1.0e-8;

// Use Newton's method to find root. We use the 90-deg, large-time approximation as an initial guess
   // while(i < imax && eps > tol) {
   //    t0 = t1;
   //    t1 = t0 - (GaussianWidth2(v, mu, B0, Omega, t0) + C) / dGaussianWidth2(v, mu, Omega, t0);
   //    eps = fabs(t1 - t0) / (fabs(t1) + fabs(t0));
   //    i++;
   // };

// Output warning if max number of interations was reached
   if(i == imax) std::cerr << "Max number of iterations reached with eps = " << eps << std::endl;

   return t1;
};

/*!
\brief Compute the analytic, improper (0 to infinity) time integral of the characteristic function for the 90-deg, large-time approximation
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v     particle speed
\param[in] mu    partickle pitch-angle
\param[in] Omega cyclotron frequency
\param[in] k     wave-number
\return analytic, improper (0 to infinity) time integral of the characteristic function for the 90-deg, large-time approximation
*/
double IntegralCharacteristic90DegLargeTime(double v, double mu, double Omega, double k)
{
   double tau_A = k * v * A;
   double tau_mu = k * v * mu;
   return 0.5 * sqrt(M_PI) / tau_A * (exp(-Sqr((tau_mu + Omega)/tau_A)) + exp(-Sqr((tau_mu - Omega)/tau_A)));
};

/*!
\brief Integrate the characteristic function over time using the trapezoid rule
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v     particle speed
\param[in] mu    partickle pitch-angle
\param[in] B0    mean magnetic field
\param[in] Omega cyclotron frequency
\param[in] k     wave-number
\return integral of the characteristic function over time
*/
double IntegralCharacteristic(double v, double mu, double B0, double Omega, double k)
{
   int i;
   double freq = fmax(fabs(k * v * mu + Omega), fabs(k * v * mu - Omega)); // Largest frequency to integrate
   double Tf = ExpLimTime(v, mu, B0, Omega, k); // Final time of integration
   double dt = M_2PI / (Ntpc * freq); // Increment resolution for integral

// Integrate resonance function using trapezoid rule
   double t = 0.0;
   double S = 0.5 * Characteristic(v, mu, B0, Omega, k, t);
   while(t < Tf) {
      t += dt;
      S += Characteristic(v, mu, B0, Omega, k, t);
   };
   t += dt;
   S += 0.5 * Characteristic(v, mu, B0, Omega, k, t);

// The 0.5 factor is because the characteristic function should actually have a 0.5 factor but it's faster to multiply it after the integral
   return 0.5 * S * dt;
};

/*!
\brief Integrate the PSD multiplied by the integral of the characteristic function over wavenumber using the trapezoid rule
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v     particle speed
\param[in] mu    partickle pitch-angle
\param[in] B0    mean magnetic field
\param[in] Omega cyclotron frequency
\return integral of the characteristic function over time
*/
double IntegralPSD_C(double v, double mu, double B0, double Omega)
{
   int i;
   double IntRes[Nk];

// Pre-compute integral of resonance function
   for(i = 0; i < Nk; i++) IntRes[i] = IntegralCharacteristic(v, mu, B0, Omega, k_vals[i]);
   // for(i = 0; i < Nk; i++) IntRes[i] = IntegralCharacteristic90DegLargeTime(v, mu, Omega, k_vals[i]);

   double S = 0.5 * (k_vals[1] - k_vals[0]) * (PSD_vals[1] * IntRes[1] + PSD_vals[0] * IntRes[0]);
// Integrate PSD using trapezoid rule
   for(i = 1; i < Nk-1; i++) S += 0.5 * (k_vals[i+1] - k_vals[i]) * (PSD_vals[i+1] * IntRes[i+1] + PSD_vals[i] * IntRes[i]);
   return S;
};

/*!
\brief Compute pitch-angle scattering coefficient according to SOQLT
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v   particle speed
\param[in] mu  partickle pitch-angle
\param[in] B0  mean magnetic field
\param[in] isp specie index
\return SOQLT pitch-angle scattering coefficient
*/
double Dmumu(double v, double mu, double B0, unsigned int isp)
{
   double Omega = CyclotronFrequency(v, B0, isp);
   ComputeA2(B0);
   return M_4PI * Sqr(Omega) * (1.0 - Sqr(mu)) * IntegralPSD_C(v, mu, B0, Omega) / Sqr(B0);
};

int main(int argc, char** argv)
{
// Specie
   int specie = Specie::alpha_particle;
   double vel = Vel(Mom(1.0 * SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle, specie), specie);

// Read PSD from file
   ReadPSD("data/k_spectra_parallel_SHS.dat");
   ReadBmagV2("data/V2_Bmag_2013_303_2014_365.dat");

// Compute kappa_parallel vs R_V2
   KappaParaVsRadius(vel);

   return 0;
};