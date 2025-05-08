#include "src/server_config.hh"
#include "src/background_solarwind_termshock.hh"
#include "src/diffusion_other.hh"
#include <iostream>
#include <iomanip>
#include <fstream>

#define COMPUTE_DIVK 1

using namespace Spectrum;

// Drift functions
inline double drift_theory(double x)
{
   return 0.457 - 0.412 * fabs(x) + 0.0915 * Sqr(x);
};

inline GeoVector drift_numer(double r_L, double vel, SpatialData spdata)
{
   return (r_L * vel / 3.0) * (spdata.curlB() - 2.0 * (spdata.gradBmag ^ spdata.bhat)) / spdata.Bmag; 
};

int main(int argc, char** argv)
{
   BackgroundSolarWindTermShock background;
   DiffusionEmpiricalSOQLTandUNLT diffusion;

   SpatialData spdata, spdata_forw, spdata_back;
   double t;
   int i,j,k,l;
   GeoVector pos, vel = gv_zeros, mom = gv_zeros;

   spdata._mask = BACKGROUND_U | BACKGROUND_B | BACKGROUND_gradB;
   spdata_forw._mask = spdata._mask;
   spdata_back._mask = spdata._mask;

// Import paramrs
   int n_sim_params = 13;
   double sim_params[n_sim_params];
   std::ifstream sim_params_file("params_drifts.txt");
   for (int i = 0; i < n_sim_params; i++) sim_params_file >> sim_params[i];
   sim_params_file.close();

   int specie;
   if (sim_params[0] == 0) specie = Specie::proton;
   else if (sim_params[0] == 1) specie = Specie::alpha_particle;
   else if (sim_params[0] == 2) specie = Specie::electron;
   else std::cout << "Specie index not recognized." << std::endl;

   DataContainer container;
   container.Clear();

// Initial time
   double t0 = 60.0 * 60.0 * 24.0 * 365.0 * 2000.5 / unit_time_fluid;
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
   double dBmag_E = sim_params[1] / unit_magnetic_fluid;
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

// WSO datafile
#if SOLARWIND_CURRENT_SHEET == 4
   std::string WSO_datafile = "data/WSO_tilt_angle_slice_LRs.dat";
   container.Insert(WSO_datafile);
#endif

// Termination shock radius
   double r_TS = 83.5 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(r_TS);

// Termination shock width
   double w_TS = 0.1 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(w_TS);

// Termination shock strength
   double s_TS = 2.5;
   container.Insert(s_TS);

   background.SetupObject(container);
   background.SetSpecie(specie);

// Clear container
   container.Clear();

// Parallel mean free path
   double lam_para = sim_params[2] * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(lam_para);

// Perpendicular mean free path
   double lam_perp = sim_params[3] * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
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
   double kap_red_fac = sim_params[4];
   container.Insert(kap_red_fac);

// Solar cycle indicator variable index
   int solar_cycle_idx = 2;
   container.Insert(solar_cycle_idx);

// Magnitude of solar cycle effect
   double solar_cycle_effect = sim_params[5];
   container.Insert(solar_cycle_effect);

   diffusion.SetupObject(container);
   diffusion.SetSpecie(specie);

//----------------------------------------------------------------------------------------------------------------------------------------------------

// 2D colormap
   int Nx = sim_params[6], Nz = sim_params[7];
   double AU = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   GeoVector sim_vel, divK, pos_tmp;
   double r_L, delta;
   double Kperp_forw, Kperp_back, Kpara_forw, Kpara_back, Kappa_forw, Kappa_back;
   double Kperp, Kpara;
   GeoVector gradKpara, gradKperp;
   GeoMatrix bhatbhat;
   double x0  = sim_params[8] * AU;
   double z0  = sim_params[9] * AU;
   double drx = sim_params[10] * AU / (double)Nx;
   double drz = sim_params[11] * AU / (double)Nz;
   mom[0] = Mom(sim_params[12] * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
   vel[0] = Vel(mom[0], specie);

   double t_min = 60.0 * 60.0 * 24.0 * 365.0 * 2007.00 / unit_time_fluid;
   double t_max = 60.0 * 60.0 * 24.0 * 365.0 * 2018.85 / unit_time_fluid;
   t = t_min + 0.9 * (t_max - t_min);

   std::ofstream drifts_file("../results/gcr_drifts.dat");
   drifts_file << COMPUTE_DIVK << std::endl;
   drifts_file << std::scientific;

   pos[0] = x0 + 0.5 * drx;
   pos[1] = 0.0;
   for (k = 0; k < Nx; k++) {
      pos[2] = z0 + 0.5 * drz;
      for (l = 0; l < Nz; l++) {
         background.GetFields(t, pos, mom, spdata);
         r_L = LarmorRadius(mom[0], spdata.Bmag, specie);

#if COMPUTE_DIVK == 0
// Drift velocity
         sim_vel = drift_numer(r_L, vel[0], spdata);
// Correct magnitude if necessary
         if (sim_vel.Norm() > 0.5 * vel[0]) {
            sim_vel.Normalize();
            sim_vel *= 0.5 * vel[0];
         };

         drifts_file << std::setw(16) << sim_vel.x / c_code
                     << std::setw(16) << sim_vel.y / c_code
                     << std::setw(16) << sim_vel.z / c_code;

#elif COMPUTE_DIVK == 1
// Divergence of diffusion using finite difference
         delta = fmin(r_L, spdata.dmax);
         divK = gv_zeros;
         for (j = 0; j < 3; j++) {
            pos_tmp = pos + delta * cart_unit_vec[j];
            background.GetFields(t, pos_tmp, mom, spdata_forw);
            Kperp_forw = diffusion.GetComponent(0, t, pos_tmp, mom, spdata_forw);
            Kpara_forw = diffusion.GetComponent(1, t, pos_tmp, mom, spdata_forw);
            pos_tmp[j] -= 2.0 * delta;
            background.GetFields(t, pos_tmp, mom, spdata_back);
            Kperp_back = diffusion.GetComponent(0, t, pos_tmp, mom, spdata_back);
            Kpara_back = diffusion.GetComponent(1, t, pos_tmp, mom, spdata_back);
            for (i = 0; i < 3; i++) {
               Kappa_forw = Kperp_forw * (i == j ? 1.0 : 0.0) + (Kpara_forw - Kperp_forw) * spdata_forw.bhat[j] * spdata_forw.bhat[i];
               Kappa_back = Kperp_back * (i == j ? 1.0 : 0.0) + (Kpara_back - Kperp_back) * spdata_back.bhat[j] * spdata_back.bhat[i];
               divK[i] += 0.5 * (Kappa_forw - Kappa_back) / delta;
            };
         };

#elif COMPUTE_DIVK == 2
// Divergence of diffusion using chain rule
         Kperp = diffusion.GetComponent(0, t, pos, mom, spdata);
         gradKperp[0] = diffusion.GetDirectionalDerivative(0);
         gradKperp[1] = diffusion.GetDirectionalDerivative(1);
         gradKperp[2] = diffusion.GetDirectionalDerivative(2);
         Kpara = diffusion.GetComponent(1, t, pos, mom, spdata);
         gradKpara[0] = diffusion.GetDirectionalDerivative(0);
         gradKpara[1] = diffusion.GetDirectionalDerivative(1);
         gradKpara[2] = diffusion.GetDirectionalDerivative(2);
         bhatbhat.Dyadic(spdata.bhat);
         divK = gradKperp + bhatbhat * (gradKpara - gradKperp)
              + (Kpara - Kperp) * (spdata.divbhat() * spdata.bhat + spdata.bhat * spdata.gradbhat());

#endif

#if COMPUTE_DIVK > 0
// Correct magnitude if necessary
         if (divK.Norm() > 0.1 * vel[0]) {
            divK.Normalize();
            divK *= 0.1 * vel[0];
         };

         drifts_file << std::setw(16) << divK.x / c_code
                     << std::setw(16) << divK.y / c_code
                     << std::setw(16) << divK.z / c_code;
#endif

         drifts_file << std::setw(16) << r_L << std::endl;
         pos[2] += drz;
      };
      pos[0] += drx;
   };

   drifts_file.close();

   return 0;
};
