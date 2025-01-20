#include "src/simulation.hh"
#include "src/distribution_other.hh"
#include "src/background_solarwind_termshock.hh"
#include "src/diffusion_other.hh"
#include "src/boundary_time.hh"
#include "src/boundary_space.hh"
#include "src/initial_time.hh"
#include "src/initial_space.hh"
#include "src/initial_momentum.hh"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <filesystem>

using namespace Spectrum;

int main(int argc, char** argv)
{

   DataContainer container;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Create a simulation object
//----------------------------------------------------------------------------------------------------------------------------------------------------

   std::unique_ptr<SimulationWorker> simulation;
   simulation = CreateSimulation(argc, argv);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Particle type
//----------------------------------------------------------------------------------------------------------------------------------------------------

   int specie = Specie::alpha_particle;
   int nuclei = 4;
   simulation->SetSpecie(specie);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Background
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Initial time
   double t0 = 60.0 * 60.0 * 24.0 * 365.0 * 2001.0 / unit_time_fluid;
   container.Insert(t0);

// Origin
   container.Insert(gv_zeros);

// Velocity
   double umag = 4.0e7 / unit_velocity_fluid;
   GeoVector u0(umag, 0.0, 0.0);
   container.Insert(u0);

// Magnetic field
   double RS = 6.957e10 / unit_length_fluid;
   double r_ref = 3.0 * RS;
   double BmagE = 6.0e-5 / unit_magnetic_fluid;
   double dBmag_E = 1.5e-5 / unit_magnetic_fluid;
   double Bmag_ref = 0.71 * BmagE * Sqr((GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid) / r_ref);
   double dBmag_ref = Bmag_ref * (dBmag_E / BmagE);
   GeoVector B0(-Bmag_ref, -dBmag_ref, 0.0);
   container.Insert(B0);

// Effective "mesh" resolution
   double dmax = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(dmax);

// solar rotation vector
   double w0 = M_2PI / (27.0 * 24.0 * 3600.0) / unit_frequency_fluid;
   GeoVector Omega(0.0, 0.0, w0);
   container.Insert(Omega);

// Reference equatorial distance
   container.Insert(r_ref);

// dmax fraction for distances closer to the Sun
   double dmax_fraction = 0.1;
   container.Insert(dmax_fraction);

// Termination shock radius
   double r_TS = 83.1 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(r_TS);

// Termination shock width
   double w_TS = 1.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(w_TS);

// Termination shock strength
   double s_TS = 2.5;
   container.Insert(s_TS);
   simulation->AddBackground(BackgroundSolarWindTermShock(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

   int perc_seg = 0;
   if(argc > 3) perc_seg = atoi(argv[3]);
   if(perc_seg < 0) perc_seg = 0;
   else if(perc_seg > 100) perc_seg = 100;

   double t_min = 60.0 * 60.0 * 24.0 * 365.0 * 2007.00 / unit_time_fluid;
   double t_max = 60.0 * 60.0 * 24.0 * 365.0 * 2018.85 / unit_time_fluid;
   double t_init = t_min + ((double)perc_seg / 100.0) * (t_max - t_min);
   container.Insert(t_init);

   simulation->AddInitial(InitialTimeFixed(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Spatial initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

   GeoVector r_min(-58.67, -42.78, -37.16);
   GeoVector r_max(-79.02, -62.40, -63.41);
   GeoVector init_pos = r_min + ((double)perc_seg / 100.0) * (r_max - r_min);
   container.Insert(init_pos);

   simulation->AddInitial(InitialSpaceFixed(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Momentum initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Lower bound for momentum
   double momentum1 = Mom(nuclei * 130.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
   container.Insert(momentum1);

// Upper bound for momentum
   double momentum2 = Mom(nuclei * 460.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
   container.Insert(momentum2);

// Log bias
   bool log_bias = true;
   container.Insert(log_bias);

   simulation->AddInitial(InitialMomentumThickShell(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Inner boundary
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings_Sun = -1;
   container.Insert(max_crossings_Sun);

// Action
   std::vector<int> actions_Sun;
   actions_Sun.push_back(-1);
   actions_Sun.push_back(-1);
   container.Insert(actions_Sun);

// Origin
   container.Insert(gv_zeros);

// Radius
   double inner_boundary = 1.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(inner_boundary);

   simulation->AddBoundary(BoundarySphereReflect(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Outer boundary
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Max crossings
   int max_crossings_outer = 1;
   container.Insert(max_crossings_outer);

// Action
   std::vector<int> actions_outer;
   actions_outer.push_back(0);
   actions_outer.push_back(0);
   container.Insert(actions_outer);

// Origin
   container.Insert(gv_zeros);

// Radius
   double outer_boundary = 119.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(outer_boundary);

   simulation->AddBoundary(BoundarySphereAbsorb(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time limit
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Not needed because this class sets the value to -1
   int max_crossings_time = 1;
   container.Insert(max_crossings_time);

// Action
   std::vector<int> actions_time;
   actions_time.push_back(-1);
   actions_time.push_back(-1);
   container.Insert(actions_time);
   
// Max duration of the trajectory
   double maxtime = 60.0 * 60.0 * 24.0 * 365.0 * 10.0 / unit_time_fluid;
   container.Insert(maxtime);

   simulation->AddBoundary(BoundaryTimeExpire(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Diffusion model
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Import diffusion parameters
   int n_diff_params = 3;
   double diff_params[n_diff_params];
   std::ifstream diff_params_file("params.dat");
   for(int i = 0; i < n_diff_params; i++) diff_params_file >> diff_params[i];
   diff_params_file.close();

   container.Clear();

// LISM indicator variable index
   int LISM_idx = 0;
   container.Insert(LISM_idx);

// Parallel inner mean free path
   double lam_in = diff_params[0] * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(lam_in);

// Parallel outer mean free path
   double lam_inc_fac = 1.0e6;
   double lam_out = lam_in * lam_inc_fac;
   container.Insert(lam_out);

// Rigidity normalization factor
   double R0 = 1.0e9 / unit_rigidity_particle;
   container.Insert(R0);

// Magnetic field normalization factor
   container.Insert(BmagE);

// Ratio of perp to para diffusion inner
   double kap_rat_in = diff_params[1];
   container.Insert(kap_rat_in);

// Ratio of perp to para diffusion outer
   double kap_rat_out = kap_rat_in / lam_inc_fac / 1.0e4;
   container.Insert(kap_rat_out);

// Bmix indicator variable index
   int Bmix_idx = 1;
   container.Insert(Bmix_idx);

// Ratio reduction factor in unipolar regions
   double kap_red_fac = diff_params[2];
   container.Insert(kap_red_fac);

// Pass ownership of "diffusion" to simulation
   simulation->AddDiffusion(DiffusionStraussEtAl2013(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 1 (spectrum)
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   MultiIndex n_bins1(100, 0, 0);
   container.Insert(n_bins1);
   
// Smallest value
   GeoVector minval1(EnrKin(momentum1, specie), 0.0, 0.0);
   container.Insert(minval1);

// Largest value
   GeoVector maxval1(EnrKin(momentum2, specie), 0.0, 0.0);
   container.Insert(maxval1);

// Linear or logarithmic bins
   MultiIndex log_bins1(1, 0, 0);
   container.Insert(log_bins1);

// Add outlying events to the end bins
   MultiIndex bin_outside1(0, 0, 0);
   container.Insert(bin_outside1);

// Physical units of the distro variable
   double unit_distro1 = 1.0;
   container.Insert(unit_distro1);

// Physical units of the bin variable
   GeoVector unit_val1 = {unit_energy_particle, 1.0, 1.0};
   container.Insert(unit_val1);

// Don't keep records
   bool keep_records1 = false;
   container.Insert(keep_records1);

//! Normalization for the "hot" boundary
   double J0 = 4.1 / nuclei;
   container.Insert(J0);

//! Characteristic energy
   double T0 = nuclei * SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle;
   container.Insert(T0);

//! Spectral power law
   double pow_law_T = 0.28;
   container.Insert(pow_law_T);

//! Constant value for the "cold" condition
   double val_cold1 = 0.0;
   container.Insert(val_cold1);

//! Bendover energy
   double Tb = 0.428 * T0;
   container.Insert(Tb);

//! Spectral power law after bend
   double pow_law_Tb = -2.54;
   container.Insert(pow_law_Tb);

//! Smoothness of bend
   double bend_smooth = 2.57;
   container.Insert(bend_smooth);

   simulation->AddDistribution(DistributionSpectrumKineticEnergyBentPowerLaw(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution 2 (exit time)
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   MultiIndex n_bins2(100, 0, 0);
   container.Insert(n_bins2);
   
// Smallest value
   GeoVector minval2(0.0, 0.0, 0.0);
   container.Insert(minval2);

// Largest value
   GeoVector maxval2(maxtime, 0.0, 0.0);
   container.Insert(maxval2);

// Linear or logarithmic bins
   MultiIndex log_bins2(0, 0, 0);
   container.Insert(log_bins2);

// Add outlying events to the end bins
   MultiIndex bin_outside2(1, 0, 0);
   container.Insert(bin_outside2);

// Physical units of the distro variable
   double unit_distro2 = 1.0;
   container.Insert(unit_distro2);

// Physical units of the bin variable
   GeoVector unit_val2 = {unit_time_fluid, 1.0, 1.0};
   container.Insert(unit_val2);

// Keep records
   bool keep_records2 = false;
   container.Insert(keep_records2);

// Value for the "hot" condition
   double val_hot2 = 1.0;
   container.Insert(val_hot2);

// Value for the "cold" condition
   double val_cold2 = 0.0;
   container.Insert(val_cold2);

// Which time to take (initial or final)
   double val_time2 = 1;
   container.Insert(val_time2);
   
   simulation->AddDistribution(DistributionTimeUniform(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Run the simulation
//----------------------------------------------------------------------------------------------------------------------------------------------------

   int n_traj;
   int batch_size;

   batch_size = n_traj = 1;
   if(argc > 1) n_traj = atoi(argv[1]);
   if(argc > 2) batch_size = atoi(argv[2]);

   std::filesystem::create_directory("../results/HS_mod_spec_He");
   std::string simulation_files_prefix = "../results/HS_mod_spec_He/HS_mod_parker_" + std::to_string(perc_seg) + "_pct_";

   simulation->DistroFileName(simulation_files_prefix);
   simulation->SetTasks(n_traj, batch_size);
   simulation->MainLoop();

   return 0;
};
