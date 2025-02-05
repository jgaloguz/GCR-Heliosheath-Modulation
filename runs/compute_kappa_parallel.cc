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
\param[in] out_fp path to output file
\param[in] B0     mean magnetic field
\param[in] isp    specie index
*/
void PlotKappaParaVsRigidity(std::string out_fp, double B0, unsigned int isp)
{
   int i;
   double vel, mom, E;
   double Kpara;
   std::ofstream file;

// Status message
   std::cout << "Writing parallel diffusion coefficient vs rigidity..." << std::endl;

// Iterate over E values
   file.open(out_fp); 
   for (i = 0; i < NE; i++) {
      E = exp(log(E1) + i * dlnE);
      mom = Mom(E, isp);
      vel = Vel(mom, isp);
      Kpara = KappaPara(vel, B0, isp);
      std::cout << "E = " << E * unit_energy_particle / SPC_CONST_CGSM_MEGA_ELECTRON_VOLT << " MeV" << std::endl;
      file << std::setw(18) << Rigidity(mom, isp) * unit_rigidity_particle * 300.0 / 1.0e9 // GV
           << std::setw(18) << Kpara * unit_diffusion_fluid    // diffusion coefficient (cm^2 s^-1)
           << std::setw(18) << Kpara * 3.0 / vel               // mean free path (code units = au)
           << std::endl;
   };
   file.close();
};

/*!
\brief "Plot" the parallel diffusion coefficient along V2 trajectory
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] out_fp path to output file
\param[in] v      particle speed
\param[in] isp    specie index
*/
void PlotKappaParaVsRadius(std::string out_fp, double v, unsigned int isp)
{
   int i;
   std::ofstream file;

// Status message
   std::cout << "Writing parallel diffusion coefficient along V2 trajectory..." << std::endl;

// Iterate over E values
   file.open(out_fp);
   for (i = 0; i < NB; i++) {
      std::cout << "R = " << R_V2[i] << " au" << std::endl;
      file << std::setw(18) << R_V2[i]
           << std::setw(18) << KappaPara(v, Bmag_V2[i], isp) * unit_diffusion_fluid
           << std::endl;
   };
   file.close();
};