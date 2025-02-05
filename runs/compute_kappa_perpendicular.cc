#include "compute_kappa_perpendicular.hh"

using namespace Spectrum;

/*!
\brief Read the computed values of parallel diffusion coefficient vs rigidity from file
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] k_para_fp path to file containing coefficients
*/
void ReadKappaParaVsRigidity(std::string k_para_fp)
{
   int i;
   double rig, mfp;
   std::ifstream file;

// Status message
   std::cout << "Reading parallel diffusion coefficient vs rigidity..." << std::endl;

// Open file with PSD
   file.open(k_para_fp);
// Iterate over k values
   for (i = 0; i < NE; i++) {
      file >> rig;
      file >> kappa_para_rig[i];
      file >> mfp;
      kappa_para_rig[i] /= unit_diffusion_fluid;
   };
// Close file
   file.close();
};

/*!
\brief Read the computed values of parallel diffusion coefficient vs radius from file
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] k_para_fp path to file containing coefficients
*/
void ReadKappaParaVsRadius(std::string k_para_fp)
{
   int i;
   std::ifstream file;

// Status message
   std::cout << "Reading parallel diffusion coefficient along V2 trajectory..." << std::endl;

// Open file with PSD
   file.open(k_para_fp);
// Iterate over k values
   for (i = 0; i < NB; i++) {
      file >> R_V2[i];
      file >> kappa_para_V2[i];
      kappa_para_V2[i] /= unit_diffusion_fluid;
   };
// Close file
   file.close();
};

/*!
\brief "Plot" the perpendicular diffusion coefficient vs rigidity
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] out_fp path to output file
\param[in] B0     mean magnetic field
\param[in] Kperp  perpendicular diffusion coefficient (to start iterations if necessary)
\param[in] isp    specie index
*/
void PlotKappaPerpVsRigidity(std::string out_fp, double B0, double Kperp, unsigned int isp)
{
   int i;
   double vel, mom, E;
   std::ofstream file;

// Status message
   std::cout << "Writing perpendicular diffusion coefficient vs rigidity..." << std::endl;

// Iterate over E values
   file.open(out_fp);
   for (i = 0; i < NE; i++) {
      E = exp(log(E1) + i * dlnE);
      mom = Mom(E, isp);
      vel = Vel(mom, isp);
      Kperp = KappaPerp(vel, B0, kappa_para_rig[i], Kperp);
      std::cout << "E = " << E * unit_energy_particle / SPC_CONST_CGSM_MEGA_ELECTRON_VOLT << " MeV" << std::endl;
      file << std::setw(18) << Rigidity(mom, isp) * unit_rigidity_particle * 300.0 / 1.0e9 // GV
           << std::setw(18) << Kperp * unit_diffusion_fluid    // diffusion coefficient (cm^2 s^-1)
           << std::setw(18) << Kperp * 3.0 / vel               // mean free path (code units = au)
           << std::endl;
   };
   file.close();
};

/*!
\brief "Plot" the perpendicular diffusion coefficient along V2 trajectory
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] out_fp path to output file
\param[in] v      particle speed
\param[in] Kperp  perpendicular diffusion coefficient (to start iterations if necessary)
\param[in] isp    specie index
*/
void PlotKappaPerpVsRadius(std::string out_fp, double v, double Kperp, unsigned int isp)
{
   int i;
   std::ofstream file;

// Status message
   std::cout << "Writing perpendicular diffusion coefficient along V2 trajectory..." << std::endl;

// Iterate over E values
   file.open(out_fp);
   for (i = 0; i < NB; i++) {
      Kperp = KappaPerp(v, Bmag_V2[i], kappa_para_V2[i], Kperp);
      std::cout << "R = " << R_V2[i] << " au" << std::endl;
      file << std::setw(18) << R_V2[i]
           << std::setw(18) << Kperp * unit_diffusion_fluid
           << std::endl;
   };
   file.close();
};