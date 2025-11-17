#include "src/server_config.hh"
#include "src/background_solarwind_termshock.hh"
#include "src/diffusion_other.hh"
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace Spectrum;

int main(int argc, char** argv)
{
   BackgroundSolarWindTermShock background;
   DiffusionEmpiricalSOQLTandUNLT diffusion;
   std::ofstream diffusion_file;

   SpatialData spdata;
   int specie = Specie::alpha_particle;
   double t;
   int i,j,k;
   GeoVector pos, vel = gv_zeros, mom = gv_zeros;

   spdata._mask = BACKGROUND_U | BACKGROUND_B;

   DataContainer container;

// Import simulation parameters
   int n_sim_params = 17;
   double sim_params[n_sim_params];
   std::ifstream diff_sim_s_file("params_He.txt");
   for(i = 0; i < n_sim_params; i++) diff_sim_s_file >> sim_params[i];
   diff_sim_s_file.close();

   container.Clear();

// Initial time
   double one_year = 60.0 * 60.0 * 24.0 * 365.0 / unit_time_fluid;
   double t0 = 2000.5 * one_year;
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
   double dBmag_E = sim_params[0] / unit_magnetic_fluid;
   double one_au = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   double Bmag_ref = 0.71 * BmagE * Sqr(one_au / r_ref);
   double dBmag_ref = Bmag_ref * (dBmag_E / BmagE);
   GeoVector B0(-Bmag_ref, -dBmag_ref, 0.0);
   container.Insert(B0);

// Effective "mesh" resolution
   double dmax = one_au;
   container.Insert(dmax);

// solar rotation vector
   double w0 = M_2PI / (30.0 * 24.0 * 3600.0) / unit_frequency_fluid;
   GeoVector Omega(0.0, 0.0, w0);
   container.Insert(Omega);

// Reference equatorial distance
   container.Insert(r_ref);

// dmax fraction for distances closer to the Sun
   double dmax_fraction = 0.1;
   container.Insert(dmax_fraction);

// WSO datafile
#if SOLARWIND_CURRENT_SHEET == 4
   std::string WSO_datafile = "data/WSO_tilt_angle_slice_LRs.dat";
   container.Insert(WSO_datafile);
#endif

// Turbulence structures
#if N_ADV_TUR_STR > 0
   double t0_trb[N_ADV_TUR_STR];
   double sig_trb[N_ADV_TUR_STR];
   double amp_trb[N_ADV_TUR_STR];
   for (i = 0; i < N_ADV_TUR_STR; i++) {
      t0_trb[i] = sim_params[5 + 3 * i] * one_year;
      sig_trb[i] = sim_params[6 + 3 * i] * one_au;
      amp_trb[i] = sim_params[7 + 3 * i] / sqrt(M_2PI) / sig_trb[i];
   };
   for (i = 0; i < N_ADV_TUR_STR; i++) {
      container.Insert(t0_trb[i]);
      container.Insert(sig_trb[i]);
      container.Insert(amp_trb[i]);
   };
#endif

// Termination shock radius
   double r_TS = 83.1 * one_au;
   container.Insert(r_TS);

// Termination shock width
   double w_TS = 1.0 * one_au;
   container.Insert(w_TS);

// Termination shock strength
   double s_TS = 2.5;
   container.Insert(s_TS);

   background.SetupObject(container);

//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Parallel mean free path
   double lam_para = sim_params[1] * one_au;
   container.Insert(lam_para);

// Perpendicular mean free path
   double lam_perp = sim_params[2] * one_au;
   container.Insert(lam_perp);

// Rigidity normalization factor
   double R0 = 3.33e7 / unit_rigidity_particle;
   container.Insert(R0);

// Magnetic field normalization factor
   container.Insert(BmagE);

// Bmix indicator variable index
   int Bmix_idx = 1;
   container.Insert(Bmix_idx);

// Ratio reduction factor in unipolar regions
   double kap_red_fac = sim_params[3];
   container.Insert(kap_red_fac);

// Solar cycle indicator variable index
   int solar_cycle_idx = 2;
   container.Insert(solar_cycle_idx);

// Magnitude of solar cycle effect
   double solar_cycle_effect = sim_params[4];
   container.Insert(solar_cycle_effect);

// Pass ownership of "diffusion" to simulation
   diffusion.SetupObject(container);
   diffusion.SetSpecie(specie);

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Number of frames
   int NE = 50;
   int Nt = 1000;
   if(argc > 1) NE = atoi(argv[1]);
   if(argc > 2) Nt = atoi(argv[2]);

// Time bounds for simulation
   double t_min = 2007.00 * one_year;
   double t_max = 2018.85 * one_year;
   double dt = (t_max - t_min) / (Nt-1);
   t = t_min;

// Radial bounds for V2 trajectory
   GeoVector r_min(-58.67, -42.78, -37.16);
   GeoVector r_max(-79.02, -62.40, -63.41);
   double ang = -atan2(r_max[1],r_max[0]);
   r_min.Rotate(gv_nz, ang);
   r_max.Rotate(gv_nz, ang);
   GeoVector dr = (r_max - r_min) / (Nt-1);

// Energy bounds for rigidity plot
   double Emin = 100.0 * SPC_CONST_CGSM_KILO_ELECTRON_VOLT / unit_energy_particle;
   double Emax = 100.0 * SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle;
   double dlnE = (log(Emax) - log(Emin)) / (NE-1);

// Output
   double E, Kpara, Kperp;
   double perc = 0.63;                                      // fraction of V2 path
   t = t_min + perc * (Nt-1) * dt;
   pos = r_min + perc * (Nt-1) * dr;
   background.GetFields(t, pos, mom, spdata);               // get region information from model
   diffusion_file.open("../results/kappa_SIM_rig_He.dat");
   std::cout << "r = " << pos.Norm() << " au" << std::endl;
   for(i = 0; i < NE; i++) {
      E = exp(log(Emin) + i * dlnE);
      mom[0] = Mom(E, specie);
      Kpara = diffusion.GetComponent(1, t, pos, mom, spdata);
      Kperp = diffusion.GetComponent(0, t, pos, mom, spdata);
      diffusion_file << std::setw(18) << Rigidity(mom[0], specie) * unit_rigidity_particle * 300.0 / 1.0e9
                     << std::setw(18) << Kpara * unit_diffusion_fluid
                     << std::setw(18) << Kpara * 3.0 / Vel(mom[0], specie)
                     << std::setw(18) << Kperp * unit_diffusion_fluid
                     << std::setw(18) << Kperp * 3.0 / Vel(mom[0], specie)
                     << std::endl;
   };
   diffusion_file.close();

// Output diffusion along V2 trajectory
   t = t_min;
   pos = r_min;
   diffusion_file.open("../results/kappa_SIM_V2_He.dat");
   mom[0] = Mom(1.0 * SPC_CONST_CGSM_GIGA_ELECTRON_VOLT / unit_energy_particle, specie);
   std::cout << "R = " << Rigidity(mom[0], specie) * unit_rigidity_particle * 300.0 / 1.0e9 << " GV" << std::endl;
   for(i = 0; i < Nt; i++) {
      background.GetFields(t, pos, mom, spdata);
      diffusion_file << std::setw(18) << pos.Norm()
                     << std::setw(18) << diffusion.GetComponent(1, t, pos, mom, spdata) * unit_diffusion_fluid
                     << std::setw(18) << diffusion.GetComponent(0, t, pos, mom, spdata) * unit_diffusion_fluid
                     << std::endl;
      t += dt;
      pos += dr;
   };
   diffusion_file.close();

// Print diffusion coefficients at 1 au for a sanity check
   pos = gv_nx;
   t = t0;
   background.GetFields(t, pos, mom, spdata);
   std::cout << "Solar Max @ 1 au:" << std::endl
             << "\tk_para = " << diffusion.GetComponent(1, t, pos, mom, spdata) * unit_diffusion_fluid
             << " cm^2 s^-1" << std::endl
             << "\tk_perp = " << diffusion.GetComponent(0, t, pos, mom, spdata) * unit_diffusion_fluid
             << " cm^2 s^-1" << std::endl;
   t = t0 + 5.5 * one_year;
   background.GetFields(t, pos, mom, spdata);
   std::cout << "Solar Min @ 1 au:" << std::endl
             << "\tk_para = " << diffusion.GetComponent(1, t, pos, mom, spdata) * unit_diffusion_fluid
             << " cm^2 s^-1" << std::endl
             << "\tk_perp = " << diffusion.GetComponent(0, t, pos, mom, spdata) * unit_diffusion_fluid
             << " cm^2 s^-1"<< std::endl;

// Print normalized DSA parameters
   t = t_min;
   pos = r_min;
   background.GetFields(t, pos, mom, spdata);
   Kperp = diffusion.GetComponent(0, t, pos, mom, spdata);
   std::cout << "normalized shock width = " << w_TS / (Kperp * umag) << std::endl;
   std::cout << "shock planarity parameter = " << umag * r_TS / Kperp << std::endl;

   return 0;
};
