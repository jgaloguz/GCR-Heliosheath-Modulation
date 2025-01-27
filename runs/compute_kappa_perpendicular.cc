#include "compute_kappa_perpendicular.hh"

using namespace Spectrum;

/*!
\brief "Plot" the perpendicular diffusion coefficient vs rigidity
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] B0    mean magnetic field
\param[in] Kpara parallel diffusion coefficient
\param[in] Kperp perpendicular diffusion coefficient (to start iterations if necessary)
\param[in] isp   specie index
*/
void PlotKappaPerpVsRigidity(double B0, double Kpara, double Kperp, unsigned int isp)
{
   int i;
   double E;

// Iterate over E values
   for (i = 0; i < NE; i++) {
      E = exp(log(E1) + i * dlnE);
      Kperp = KappaPerp(Vel(Mom(E, isp), isp), B0, Kpara, Kperp);
      std::cout << std::setw(18) << Rigidity(Mom(E, isp), isp) / 1.0e9
                << std::setw(18) << Kperp * unit_diffusion_fluid
                << std::endl;
   };
};

/*!
\brief "Plot" the perpendicular diffusion coefficient along V2 trajectory
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v     particle speed
\param[in] Kperp perpendicular diffusion coefficient (to start iterations if necessary)
*/
void PlotKappaPerpVsRadius(double v, double Kperp)
{
   int i;

// Iterate over E values
   for (i = 0; i < NB; i++) {
      Kperp = KappaPerp(v, Bmag_V2[i], kappa_para_V2[i], Kperp);
      std::cout << std::setw(18) << R_V2[i]
                << std::setw(18) << Kperp * unit_diffusion_fluid
                << std::endl;
   };
};