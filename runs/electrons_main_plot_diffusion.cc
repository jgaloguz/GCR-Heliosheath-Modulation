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
   DiffusionPotgieterEtAl2015 diffusion;
   std::ofstream diffusion_file;

   SpatialData spdata;
   int specie = Specie::electron;
   double t;
   int i,j,k;
   GeoVector pos, vel = gv_zeros, mom = gv_zeros;

   spdata._mask = BACKGROUND_U | BACKGROUND_B;

   DataContainer container;
   container.Clear();

// Initial time
   double t0 = 60.0 * 60.0 * 24.0 * 365.0 * 2002.0 / unit_time_fluid;
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
   double w0 = M_2PI / (30.0 * 24.0 * 3600.0) / unit_frequency_fluid;
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

   background.SetupObject(container);

//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// LISM indicator variable index
   int LISM_idx = 0;
   container.Insert(LISM_idx);

// Parallel inner mean free path
   double lam_in = 0.23 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(lam_in);

// Parallel outer mean free path
   double lam_inc_fac = 1.0e6;
   double lam_out = lam_in * lam_inc_fac;
   container.Insert(lam_out);

// Rigidity normalization factor
   double R0 = 3.33e7 / unit_rigidity_particle;    // 10 GeV
   container.Insert(R0);

// Magnetic field normalization factor
   container.Insert(BmagE);

// Ratio of perp to para diffusion inner
   double kap_rat_in = 0.03;
   container.Insert(kap_rat_in);

// Ratio of perp to para diffusion outer
   double kap_rat_out = kap_rat_in / lam_inc_fac / 1.0e4;
   container.Insert(kap_rat_out);

// Bmix indicator variable index
   int Bmix_idx = 1;
   container.Insert(Bmix_idx);

// Ratio reduction factor in unipolar regions
   double kap_red_fac = 0.05;
   container.Insert(kap_red_fac);

// Pass ownership of "diffusion" to simulation
   diffusion.SetupObject(container);
   diffusion.SetSpecie(specie);

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Number of frames
   int NE = 50;
   int Nt = 100;
   if(argc > 1) NE = atoi(argv[1]);
   if(argc > 2) Nt = atoi(argv[2]);

// Time bounds for simulation
   double t_min = 60.0 * 60.0 * 24.0 * 365.0 * 2007.00 / unit_time_fluid;
   double t_max = 60.0 * 60.0 * 24.0 * 365.0 * 2018.85 / unit_time_fluid;
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
   diffusion_file.open("../results/kappa_SIM_rig_e.dat");
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
   diffusion_file.open("../results/kappa_SIM_V2_e.dat");
   mom[0] = Mom(10.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
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
   t = 60.0 * 60.0 * 24.0 * 365.0 * 2002.00 / unit_time_fluid;
   background.GetFields(t, pos, mom, spdata);
   std::cout << "Solar Max @ 1 au:" << std::endl
             << "\tk_para = " << diffusion.GetComponent(1, t, pos, mom, spdata) * unit_diffusion_fluid
             << " cm^2 s^-1" << std::endl
             << "\tk_perp = " << diffusion.GetComponent(0, t, pos, mom, spdata) * unit_diffusion_fluid
             << " cm^2 s^-1" << std::endl;
   t = 60.0 * 60.0 * 24.0 * 365.0 * 2007.00 / unit_time_fluid;
   background.GetFields(t, pos, mom, spdata);
   std::cout << "Solar Min @ 1 au:" << std::endl
             << "\tk_para = " << diffusion.GetComponent(1, t, pos, mom, spdata) * unit_diffusion_fluid
             << " cm^2 s^-1" << std::endl
             << "\tk_perp = " << diffusion.GetComponent(0, t, pos, mom, spdata) * unit_diffusion_fluid
             << " cm^2 s^-1" << std::endl;

   return 0;
};


