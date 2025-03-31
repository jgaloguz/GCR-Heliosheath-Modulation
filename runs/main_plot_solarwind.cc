#include "src/server_config.hh"
#include "src/background_solarwind_termshock.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>

using namespace Spectrum;

int main(int argc, char** argv)
{
   BackgroundSolarWindTermShock background;
   std::ofstream trajectory_file;

   SpatialData spdata;
   double t;
   int i,j,k;
   GeoVector pos, vel = gv_zeros, mom = gv_zeros;

   spdata._mask = BACKGROUND_U | BACKGROUND_B;

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

// WSO datafile
#if SOLARWIND_CURRENT_SHEET == 4
   std::string WSO_datafile = "data/WSO_tilt_angle_slice_LRs.dat";
   container.Insert(WSO_datafile);
#endif

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

   container.Clear();

//----------------------------------------------------------------------------------------------------------------------------------------------------

// Number of frames
   int Nt = 25;
   if(argc > 1) Nt = atoi(argv[1]);

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
   pos = r_min;

// 2D box plot with SILO functions
   GeoVector voyager[Nt];
   GeoVector xyz_min(   0.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid,
                     -120.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid,
                        0.0);
   GeoVector xyz_max( 120.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid,
                        0.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid,
                        0.0);
   MultiIndex dims_z(1000, 1000, 1);
   GeoVector normal(0.0, 1.0, 0.0);
   GeoVector right(1.0, 0.0, 0.0);
   background.SetBox(xyz_min, xyz_max, dims_z, normal, right);

// Check orientation (Bz along z-axis < 0 @ tmin)
   background.GetFields(t_min, gv_nz, mom, spdata);
   std::cerr << "Bz = " << spdata.Bvec[2] << std::endl;
   std::cerr << "Bz should be negative." << std::endl;

// Output files for visualization
   std::string frame;
   int it, it2;
   // std::filesystem::create_directory("../results/solarwind");
   // for(it = 0; it < Nt; it++) {
   //    frame = std::to_string(it);
   //    frame.insert(0, 5 - frame.size(), '0');
   //    std::cerr << "frame " << frame << std::endl;
   //    background.BoxPlot2DMesh("../results/solarwind/sw_" + frame + ".silo", false);
   //    background.BoxPlot2DScalar("By", false, t);
   //    background.BoxPlot2DScalar("Region2", false, t);
   //    background.BoxPlot2DScalar("Umag", false, t);
   //    background.BoxPlotFinalize();

   //    trajectory_file.open("../results/solarwind/traj_" + frame + ".lines");
   //    for(it2 = 0; it2 < it; it2++) {
   //       trajectory_file << std::setw(18) << voyager[it2][0] << ","
   //                       << std::setw(18) << voyager[it2][2] << ","
   //                       << std::setw(18) << 0.0
   //                       << std::endl;
   //    };
   //    voyager[it] = r_min + it * dr;
   //    trajectory_file << std::setw(18) << voyager[it][0] << ","
   //                    << std::setw(18) << voyager[it][2] << ","
   //                    << std::setw(18) << 0.0
   //                    << std::endl;
   //    trajectory_file.close();

   //    t += dt;
   // };

// Get tilt angle at Voyager trajectory
   double V2_lat, CS_lat;
   GeoVector pos_voy;
   t = t_min;
   pos = r_min;
   trajectory_file.open("../results/V2_vs_HCS_lat.dat");
   for (it = 0; it < Nt; it++) {
      background.GetFields(t, pos, mom, spdata);
      pos_voy = pos;
      pos_voy.XYZ_RTP();
      V2_lat = RadToDeg(pos_voy[1] - M_PI_2);
      CS_lat = RadToDeg(spdata.region[0]);
      trajectory_file << std::setw(18) << pos_voy[0]
                      << std::setw(18) << t * unit_time_fluid / (60.0 * 60.0 * 24.0 * 365.0)
                      << std::setw(18) << V2_lat
                      << std::setw(18) << CS_lat
                      << std::setw(18) << spdata.region[1]
                      << std::setw(18) << cos(spdata.region[2])
                      << std::endl;
      pos += dr;
      t += dt;
   };
   trajectory_file.close();


   return 0;
};


