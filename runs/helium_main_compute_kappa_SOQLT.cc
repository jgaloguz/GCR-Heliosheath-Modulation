#include "compute_kappa_parallel_SOQLT.hh"

using namespace Spectrum;

int main(int argc, char** argv)
{
// Specie
   int specie = Specie::alpha_particle;
   double vel = Vel(Mom(1.0 * SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle, specie), specie);

// Read PSD and V2 data from files
   ReadPSD("data/k_spectra_parallel_SHS.dat", M_8PI);
   ReadBmagV2("data/V2_Bmag_2013_303_2014_365.dat");

// Iterate over voyager trajectory and average over values
   double B0 = Average(NB, Bmag_V2, true);
   std::cout << "B0 = " << B0 * unit_magnetic_fluid * 1.0e5 << " nT" << std::endl;

// Compute kappa_parallel vs rigidity
   PlotKappaParaVsRigidity("../results/kappa_SOQLT_rig_He.dat", B0, specie);

// Compute kappa_parallel vs R_V2
   PlotKappaParaVsRadius("../results/kappa_SOQLT_V2_He.dat", vel, specie);

   return 0;
};