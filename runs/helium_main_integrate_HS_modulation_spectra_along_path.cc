#include "common/physics.hh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Spectrum;

const int nuc = 4;
const double J0 = 4.1 / nuc;
const double mu0 = 0.28;
const double T0 = nuc * 1000.0;
const double Tb = 0.428 * T0;
const double mu1 = -2.54;
const double ds = 2.57;

// Unmodulated spectrum
inline double unmod_spectrum(double T) {return J0 * pow(T / T0, mu0) / pow(1.0 + pow(T / Tb, (mu0-mu1)/ds), ds);};

// Integrate differential intensity (trapezoid rule)
double integrate_diff_int(double *J, double *T, double *R, int N)
{
   double S = 0.0;
   for(int i = 0; i < N - 1; i++) S += 0.5 * (J[i]*R[i] + J[i+1]*R[i+1]) * (T[i+1] - T[i]);
   return S;
};

// Interpolate response function
double interp_response(double energy, double* energy_full, double* response_full, int N)
{
   int idx = LocateInArray(0, N-1, energy_full, energy, true);
   if(idx == -1) return response_full[0];
   else if(idx == N-1) return response_full[N-1];
   else return response_full[idx] + (energy - energy_full[idx]) * (response_full[idx+1] - response_full[idx]) / (energy_full[idx+1] - energy_full[idx]);
};

int main(int argc, char** argv)
{
   std::ifstream input_spectrum_file;
   std::ofstream output_spectrum_file;
   std::ifstream response_file;
   std::string infilename;
   std::string outfilename = "../results/HS_mod_spec_He/HS_mod_parker_integ_spec.dat";
   std::string resfilename = "data/TET_D123_response_m2sr.dat";
   int eng, n_eng = 100;
   int seg, n_seg = 1;
   if(argc > 1) n_seg = atoi(argv[1]);
   int percent, d_percent = 100 / n_seg;
   double r_Tshock = 81.6 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   double r_HP = 118.9 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   double r_init, d_r = r_HP - r_Tshock;
   double energy[n_eng], diff_int[n_eng], sum_w1[n_eng];
   int n_res = 1361;
   double energy_full[n_res], response_full[n_res], response[n_eng];
   double skip;

// Get response function
   response_file.open(resfilename);
   for(eng = 0; eng <= n_res; eng++) {
      response_file >> energy_full[eng];
      response_file >> response_full[eng];
   };
   response_file.close();

// Open output distro file
   output_spectrum_file.open(outfilename);

// Iterate over segments endpoints
   for(seg = 0; seg <= n_seg; seg++) {

// Status message
      std::cerr << "\tSegment endpoint " << seg << std::endl;

      percent = seg * d_percent;
      r_init = r_Tshock + ((double)percent / 100.0) * d_r;

// Open input differential intensity file
      infilename = "../results/HS_mod_spec/HS_mod_parker_" + std::to_string(percent) + "_pct_spec_comp.dat";
      input_spectrum_file.open(infilename);

// Read differential intensity
      for(eng = 0; eng < n_eng; eng++) {
         input_spectrum_file >> energy[eng];
         energy[eng] *= 1000;
         response[eng] = interp_response(energy[eng], energy_full, response_full, n_res);
         input_spectrum_file >> skip;
         input_spectrum_file >> diff_int[eng];
      };

// Output integrated differential intensity
      output_spectrum_file << std::setw(20) << r_init
                           << std::setw(20) << integrate_diff_int(diff_int, energy, response, n_eng)
                           << std::endl;

// Close input differential intensity file
      input_spectrum_file.close();
   };

// Close output distro file
      output_spectrum_file.close();

// Output LISM integrated differential intensity for reference
   for(eng = 0; eng < n_eng; eng++) diff_int[eng] = unmod_spectrum(energy[eng]);
   std::cout << std::setw(20) << integrate_diff_int(diff_int, energy, response, n_eng) << std::endl;

   return 0;
};
