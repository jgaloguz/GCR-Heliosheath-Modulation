#include "main_compute_kappa_common.hh"

using namespace Spectrum;

//! Arrays with k and PSD values
double k_vals[Nk], PSD_vals[Nk];

//! Array with magnitude of radius and B along V2 trajectory
double R_V2[NB], Bmag_V2[NB];

/*!
\brief Read the turbulence spectrum from file
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] PSD_filename file containing spectrum
\note the input k spectrum is assumed to be in units of nT^2 * km, with k given in units of km^-1
*/
void ReadPSD(std::string PSD_filename)
{
   int i;
   std::ifstream file;

// Open file with PSD
   file.open(PSD_filename);
// Iterate over k values
   for(i = 0; i < Nk; i++) {
      file >> k_vals[i];
      file >> PSD_vals[i];
      k_vals[i] *= (unit_length_fluid * 1.0e-5); // km^-1 to cm^-1
      PSD_vals[i] /= (unit_length_fluid * 1.0e-5) * Sqr(unit_magnetic_fluid * 1.0e5); // nT^2 km to G^2 cm
      PSD_vals[i] /= M_8PI;
   };
// Close file
   file.close();
};

/*!
\brief Read the V2 magnetic field measurements from file
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] V2_filename file containing data
\note the magnetic field measurements are in nT, and the radius in au
*/
void ReadBmagV2(std::string V2_filename)
{
   int i;
   std::ifstream file;

// Open file with PSD
   file.open(V2_filename);
// Iterate over k values
   for(i = 0; i < NB; i++) {
      file >> R_V2[i];
      file >> Bmag_V2[i];
      Bmag_V2[i] /= unit_magnetic_fluid * 1.0e5; // nT to G
   };
// Close file
   file.close();
};

/*!
\brief Compute the parallel diffusion coefficient
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v   particle speed
\param[in] B0  mean magnetic field
\param[in] isp specie index
\return parallel diffusion coefficient
*/
double KappaPara(double v, double B0, unsigned int isp)
{
   int i;
   double mu = mu0;
   double S = 0.5 * Sqr(1.0 - Sqr(mu)) / Dmumu(v, mu, B0, isp);

// Integrate using trapezoid rule
   for(i = 1; i < Nmu-1; i++) {
      mu += dmu;
      S += Sqr(1.0 - Sqr(mu)) / Dmumu(v, mu, B0, isp);
   };
   // S += 0.0;
   return 0.25 * Sqr(v) * dmu * S;
};

/*!
\brief "Plot" the QLT parallel diffusion coefficient vs rigidity
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] B0  mean magnetic field
\param[in] isp specie index
*/
void KappaParaVsRigidity(double B0, unsigned int isp)
{
   int i;
   double E;

// Iterate over E values
   for(i = 0; i < NE; i++) {
      E = exp(log(E1) + i * dlnE);
      std::cout << std::setw(18) << Rigidity(Mom(E, isp), isp) / 1.0e9
                << std::setw(18) << KappaPara(Vel(Mom(E, isp), isp), B0, isp) * unit_diffusion_fluid
                << std::endl;
   };
};

/*!
\brief "Plot" the QLT parallel diffusion coefficient along V2 trajectory
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] v   particle speed
\param[in] isp specie index
*/
void KappaParaVsRadius(double v, unsigned int isp)
{
   int i;

// Iterate over E values
   for(i = 0; i < NB; i++) {
      std::cout << std::setw(18) << R_V2[i]
                << std::setw(18) << KappaPara(v, Bmag_V2[i], isp) * unit_diffusion_fluid
                << std::endl;
   };
};