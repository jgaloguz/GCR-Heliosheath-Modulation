#include "compute_kappa_perpendicular_UNLT.hh"

using namespace Spectrum;

int main(int argc, char** argv)
{
// Specie
   int specie = Specie::proton;
   double vel = Vel(Mom(100.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie), specie);

// Find region
   std::string region = "SHS";
   std::string Bmag_file = "data/V2_Bmag_2013_303_2014_365.dat";
   int Nk_reg = 40;
   int NB_reg = 428;
   if (argc > 1) {
      if ((std::string)argv[1] == "UHS") {
         region = "UHS";
         Bmag_file = "data/V2_Bmag_2010_152_2011_210.dat";
         Nk_reg = 16;
         NB_reg = 424;
      };
   };

// Read PSD, V2 data, and parallel coefficients from files
   ReadPSD("data/k_spectra_perpendicular_" + region + ".dat", Nk_reg, M_2PI);
   ReadBmagV2(Bmag_file, NB_reg);
   ReadKappaParaVsRigidity("../results/kappa_SOQLT_rig_H_" + region + ".dat");
   ReadKappaParaVsRadius("../results/kappa_SOQLT_V2_H_" + region + ".dat");

// Iterate over voyager trajectory and average over values
   double B0 = Average(NB, Bmag_V2, true);
   std::cout << "B0 = " << B0 * unit_magnetic_fluid * 1.0e5 << " nT" << std::endl;

// Compute kappa_perp vs rigidity
// Initializing iterations with Kperp = 0 will yield the FLRW limit after the first iteration
   PlotKappaPerpVsRigidity("../results/kappa_UNLT_rig_H_" + region + ".dat", B0, 0.0, specie);

// Compute kappa_perp vs R_V2
// Initializing iterations with Kperp = 0 will yield the FLRW limit after the first iteration
   PlotKappaPerpVsRadius("../results/kappa_UNLT_V2_H_" + region + ".dat", vel, 0.0, specie);

// De-allocate memory
   FreeMemory();

   return 0;
};