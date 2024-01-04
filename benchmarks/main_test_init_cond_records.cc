#include "src/simulation.hh"
#include "src/distribution_other.hh"
#include "src/background_uniform.hh"
#include "src/boundary_time.hh"
#include "src/initial_space.hh"
#include "src/initial_momentum.hh"
#include <iostream>
#include <iomanip>

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

   int specie = Specie::proton;
   simulation->SetSpecie(specie);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Background
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Origin
   container.Insert(gv_zeros);

// Velocity
   container.Insert(gv_zeros);

// Magnetic field
   double Bmag = 1.0 / unit_magnetic_fluid;
   GeoVector B0(0.0, 0.0, Bmag);
   container.Insert(B0);

// Effective "mesh" resolution
   double dmax = 0.1;
   container.Insert(dmax);

   simulation->AddBackground(BackgroundUniform(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Spatial initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

   container.Insert(gv_zeros);

   double radius = GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
   container.Insert(radius);

   simulation->AddInitial(InitialSpaceSphere(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Momentum initial condition
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Initial momentum
   double momentum = Mom(100.0 * SPC_CONST_CGSM_MEGA_ELECTRON_VOLT / unit_energy_particle, specie);
   container.Insert(momentum);

   double theta = DegToRad(90.0);
   container.Insert(theta);

   simulation->AddInitial(InitialMomentumRing(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Time boundary condition (terminal)
//----------------------------------------------------------------------------------------------------------------------------------------------------

   container.Clear();

// Not needed because this class sets the value to -1
   int max_crossings_time = 1;
   container.Insert(max_crossings_time);

// Action
   std::vector<int> actions_time;
   actions_time.push_back(0);
   container.Insert(actions_time);
   
// Duration of the trajectory
   double maxtime = 0.1;
   container.Insert(maxtime);

   simulation->AddBoundary(BoundaryTimeExpire(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Distribution
//----------------------------------------------------------------------------------------------------------------------------------------------------

// Parameters for distribution
   container.Clear();

// Number of bins
   MultiIndex n_bins(0, 0, 0);
   container.Insert(n_bins);
   
// Smallest value
   GeoVector minval(0.0, 0.0, 0.0);
   container.Insert(minval);

// Largest value
   GeoVector maxval(0.0, 0.0, 0.0);
   container.Insert(maxval);

// Linear or logarithmic bins
   MultiIndex log_bins(0, 0, 0);
   container.Insert(log_bins);

// Add outlying events to the end bins
   MultiIndex bin_outside(0, 0, 0);
   container.Insert(bin_outside);

// Physical units of the distro variable
   double unit_distro = 1.0;
   container.Insert(unit_distro);

// Physical units of the bin variable
   GeoVector unit_val = {1.0, 1.0, 1.0};
   container.Insert(unit_val);

// Don't keep records
   bool keep_records = true;
   container.Insert(keep_records);

// Value for the "hot" condition
   double val_hot = 1.0;
   container.Insert(val_hot);

// Value for the "cold" condition
   double val_cold = 0.0;
   container.Insert(val_cold);

// Coordinates to use (initial or final)
   int val_time = 0;
   container.Insert(val_time);

// Coordinate representation to use
   int val_coord = 0;
   container.Insert(val_coord);

   simulation->AddDistribution(DistributionPositionUniform(), container);

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Run the simulation
//----------------------------------------------------------------------------------------------------------------------------------------------------

   int n_traj;
   int batch_size;

   batch_size = n_traj = 1;
   if(argc > 1) n_traj = atoi(argv[1]);
   if(argc > 2) batch_size = atoi(argv[2]);

   std::string simulation_files_prefix = "output_data/main_test_init_cond_records_" + simulation->GetTrajectoryName() + "_";
   simulation->DistroFileName(simulation_files_prefix);
   simulation->SetTasks(n_traj, batch_size);
   simulation->MainLoop();
   simulation->PrintRecords(0, simulation_files_prefix + "records.dat", true);

   if(simulation->IsMaster()) {
      std::cout << std::endl;
      std::cout << "INITIAL CONDITION RECORDS" << std::endl;
      std::cout << "=========================================================" << std::endl;
      std::cout << "Trajectory type: " << simulation->GetTrajectoryName() << std::endl;
      std::cout << "=========================================================" << std::endl;
      std::cout << "Distribution files outputed to " << simulation_files_prefix << std::endl;
      std::cout << std::endl;
   };
   
   return 0;
};
