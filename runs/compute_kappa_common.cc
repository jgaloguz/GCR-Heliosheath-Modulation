#include "compute_kappa_common.hh"

using namespace Spectrum;

//! Arrays with k and PSD values
double k_vals[Nk], PSD_vals[Nk];

//! Array with magnitude of radius and B along V2 trajectory
double R_V2[NB], Bmag_V2[NB];

//! Array with parallel diffusion coefficient along V2 trajectory
double kappa_para_rig[NE], kappa_para_V2[NB];

/*!
\brief Read the turbulence spectrum from file
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] PSD_fp     path to file containing spectrum
\param[in] norm_const normalization constant relative to dB^2
\note the input k spectrum is assumed to be in units of nT^2 * km, with k given in units of km^-1
*/
void ReadPSD(std::string PSD_fp, double norm_const)
{
   int i;
   std::ifstream file;

// Status message
   std::cout << "Reading power spectrum..." << std::endl;

// Open file with PSD
   file.open(PSD_fp);
// Iterate over k values
   for (i = 0; i < Nk; i++) {
      file >> k_vals[i];
      file >> PSD_vals[i];
      k_vals[i] *= (unit_length_fluid * 1.0e-5); // km^-1 to cm^-1
      PSD_vals[i] /= (unit_length_fluid * 1.0e-5) * Sqr(unit_magnetic_fluid * 1.0e5); // nT^2 km to G^2 cm
      PSD_vals[i] /= norm_const;
   };
// Close file
   file.close();
};

/*!
\brief Read the V2 magnetic field measurements from file
\author Juan G Alonso Guzman
\date 01/24/2025
\param[in] V2_fp path to file containing data
\note the magnetic field measurements are in nT, and the radius in au
*/
void ReadBmagV2(std::string V2_fp)
{
   int i;
   std::ifstream file;

// Status message
   std::cout << "Reading V2 data..." << std::endl;

// Open file with PSD
   file.open(V2_fp);
// Iterate over k values
   for (i = 0; i < NB; i++) {
      file >> R_V2[i];
      file >> Bmag_V2[i];
      Bmag_V2[i] /= unit_magnetic_fluid * 1.0e5; // nT to G
   };
// Close file
   file.close();
};
