#include "src/distribution_other.hh"
#include "common/physics.hh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Spectrum;

const double J0 = 100.0;
const double mu0 = 0.335;
const double T0 = SPC_CONST_CGSM_GIGA_ELECTRON_VOLT;
const double Tb = 0.438 * T0;
const double mu1 = -2.52;
const double ds = 3.76;
const int specie = Specie::proton;

// Unmodulated spectrum
inline double unmod_spectrum(double T) {return J0 * pow(T / T0, mu0) / pow(1.0 + pow(T / Tb, (mu0-mu1)/ds), ds);};

int main(int argc, char** argv)
{
   int percent = 0;
   if(argc > 1) percent = atoi(argv[1]);
   std::string distroname1 = "../results/HS_mod_spec_H/HS_mod_parker_" + std::to_string(percent) + "_pct_0.out";
   std::string infilename1 = "../results/HS_mod_spec_H/HS_mod_parker_" + std::to_string(percent) + "_pct_spec.dat";
   std::string outfilename = "../results/HS_mod_spec_H/HS_mod_parker_" + std::to_string(percent) + "_pct_spec_comp.dat";
   std::string line;
   int i, N = 100;
   int sum_c1[N];
   double energy1[N], distro1[N], sum_w1[N];
   double p2, energy0 = SPC_CONST_CGSM_GIGA_ELECTRON_VOLT;

// Restore distribution
   DistributionSpectrumKineticEnergyBentPowerLaw distro_obj;
   distro_obj.Restore(distroname1);
   distro_obj.Print1D(0, infilename1, true);

// Open input analytic distro file
   std::ifstream input_spectrum_file1(infilename1);

// Read first two lines of distro file
   std::getline(input_spectrum_file1, line);
   std::getline(input_spectrum_file1, line);

// Read data
   for(i = 0; i < N; i++) {
      input_spectrum_file1 >> energy1[i];
      input_spectrum_file1 >> distro1[i];
      input_spectrum_file1 >> sum_w1[i];
      input_spectrum_file1 >> sum_c1[i];
   };

// Close input cartesian distro file
   input_spectrum_file1.close();

// Open output distro file
   std::ofstream output_spectrum_file(outfilename);

// Output data
   output_spectrum_file << std::setprecision(8);
   for(i = 0; i < N; i++) {
      p2 = Sqr(Mom(energy1[i] / unit_energy_particle, specie));
      output_spectrum_file << std::setw(20) << energy1[i] / energy0;
      output_spectrum_file << std::setw(20) << unmod_spectrum(energy1[i]);
      output_spectrum_file << std::setw(20) << p2 * distro1[i];
      output_spectrum_file << std::endl;
   };

// Close output distro file
   output_spectrum_file.close();

   return 0;
};
