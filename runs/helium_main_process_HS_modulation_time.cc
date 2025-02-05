#include "src/distribution_other.hh"
#include "common/physics.hh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Spectrum;

int main(int argc, char** argv)
{
   int percent = 0;
   if(argc > 1) percent = atoi(argv[1]);
   std::string distroname1 = "../results/HS_mod_spec_He/HS_mod_parker_" + std::to_string(percent) + "_pct_1.out";
   std::string infilename1 = "../results/HS_mod_spec_He/HS_mod_parker_" + std::to_string(percent) + "_pct_time.dat";
   std::string outfilename = "../results/HS_mod_spec_He/HS_mod_parker_" + std::to_string(percent) + "_pct_time_pp.dat";
   std::string line;
   int i, N = 100;
   int sum_c1[N], total_c = 0;
   double time1[N], distro1[N], sum_w1[N];

// Restore distribution
   DistributionTimeUniform distro_obj;
   distro_obj.Restore(distroname1);
   distro_obj.Print1D(0, infilename1, true);

// Open input analytic distro file
   std::ifstream input_spectrum_file1(infilename1);

// Read first two lines of distro file
   std::getline(input_spectrum_file1, line);
   std::getline(input_spectrum_file1, line);

// Read data
   for(i = 0; i < N; i++) {
      input_spectrum_file1 >> time1[i];
      input_spectrum_file1 >> distro1[i];
      input_spectrum_file1 >> sum_w1[i];
      input_spectrum_file1 >> sum_c1[i];
      total_c += sum_c1[i];
   };

// Close input cartesian distro file
   input_spectrum_file1.close();

// Open output distro file
   std::ofstream output_spectrum_file(outfilename);

// Output data
   output_spectrum_file << std::setprecision(8);
   for(i = 0; i < N; i++) {
      output_spectrum_file << std::setw(20) << YEARS(time1[i]);
      output_spectrum_file << std::setw(20) << (double)sum_c1[i] / (double)total_c;
      output_spectrum_file << std::endl;
   };

// Close output distro file
   output_spectrum_file.close();

   return 0;
};
