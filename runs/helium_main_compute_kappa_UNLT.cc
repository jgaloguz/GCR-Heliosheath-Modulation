#include "compute_kappa_perpendicular_UNLT.hh"

using namespace Spectrum;

int main(int argc, char** argv)
{
// Specie
   int specie = Specie::alpha_particle;
   double vel = Vel(Mom(1.0 * SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle, specie), specie);

// Read PSD, V2 data, and parallel coefficients from files
   ReadPSD("data/k_spectra_perpendicular_SHS.dat", M_2PI);
   ReadBmagV2("data/V2_Bmag_2013_303_2014_365.dat");
   ReadKappaParaVsRigidity("../results/kappa_SOQLT_rig_He.dat");
   ReadKappaParaVsRadius("../results/kappa_SOQLT_V2_He.dat");

// Iterate over voyager trajectory and average over values
   double B0 = Average(NB, Bmag_V2, true);
   std::cout << "B0 = " << B0 * unit_magnetic_fluid * 1.0e5 << " nT" << std::endl;

// Compute kappa_perp vs rigidity
// Initializing iterations with Kperp = 0 will yield the FLRW limit after the first iteration
   PlotKappaPerpVsRigidity("../results/kappa_UNLT_rig_He.dat", B0, 0.0, specie);

// Compute kappa_perp vs R_V2
// Initializing iterations with Kperp = 0 will yield the FLRW limit after the first iteration
   PlotKappaPerpVsRadius("../results/kappa_UNLT_V2_He.dat", vel, 0.0, specie);

   return 0;
};