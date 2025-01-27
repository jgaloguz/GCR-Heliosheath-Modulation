#include "compute_kappa_parallel.hh"

using namespace Spectrum;

double KappaPara(double v, double B0, unsigned int isp)
{
   int i;
   double mu = mu0;
   double S = 0.5 * Sqr(1.0 - Sqr(mu)) / Dmumu(v, mu, B0, isp);

// Integrate using trapezoid rule
   for (i = 1; i < Nmu-1; i++) {
      mu += dmu;
      S += Sqr(1.0 - Sqr(mu)) / Dmumu(v, mu, B0, isp);
   };
   // S += 0.0;
   return 0.25 * Sqr(v) * dmu * S;
};

/*!
\brief "Plot" the parallel diffusion coefficient vs rigidity
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] B0  mean magnetic field
\param[in] isp specie index
*/
void PlotKappaParaVsRigidity(double B0, unsigned int isp)
{
   int i;
   double E;

// Iterate over E values
   for (i = 0; i < NE; i++) {
      E = exp(log(E1) + i * dlnE);
      std::cout << std::setw(18) << Rigidity(Mom(E, isp), isp) / 1.0e9
                << std::setw(18) << KappaPara(Vel(Mom(E, isp), isp), B0, isp) * unit_diffusion_fluid
                << std::endl;
   };
};

/*!
\brief "Plot" the parallel diffusion coefficient along V2 trajectory
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v   particle speed
\param[in] isp specie index
*/
void PlotKappaParaVsRadius(double v, unsigned int isp)
{
   int i;

// Iterate over E values
   for (i = 0; i < NB; i++) {
      std::cout << std::setw(18) << R_V2[i]
                << std::setw(18) << KappaPara(v, Bmag_V2[i], isp) * unit_diffusion_fluid
                << std::endl;
   };
};