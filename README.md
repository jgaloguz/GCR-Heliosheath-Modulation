# GCR Heliosheath Modulation

This is a specialization of the SPECTRUM software applied to modeling galactic cosmic rays (GCR) in the heliosphere. Specifically, the codes in this repository simulate GCR modulation in the heliosheath, i.e. the region of space between the supersonic solar wind and the local interstellar medium, analyze the results, and generate useful figures for a future scientific publication.

For this application, the model of the heliosphere is analytic with the following main features:
  - the termination shock (TS) and heliopause (HP), which bound the heliosheath, are spherical,
  - the HP is the outer boundary of the simulation domain, i.e. there is no explicit model for the interstellar medium,
  - the solar wind is radially uniform up to the TS, then decreases squared inverse of radial distance,
  - the heliospheric current sheet is warped and advected with the plasma flow with a solar cycle dependent tilt angle,
  - the magnetic field has a Parker spiral configuration which adapts to the radially changing solar wind speed and has a polar correction to the azimuthal component,
  - the temporal behavior of the solar wind is data-informed, i.e. calibrated to be consistent with remote and in-situ (1 au) measurements.

The boundary conditions for GCR quantities are obtained from measurements taken by the Voyager spacecrafts, and the diffusion coefficients dictating the particle transport are empirical.
The simulations are run along a simulated trajectory of Voyager 2 (V2) within the simplified heliospheric model.

## Using the code

Before first use, the code must be configured with autotools. After cloning this repository, execute the configure script in the working directory
```
git clone https://github.com/jgaloguz/GCR-Heliosheath-Modulation
cd GCR-Heliosheath-Modulation
./configure.sh <mpi-option>
```
where `<mpi-option>` is either `openmpi` or `mpich`, whichever is installed in your system. You may have to change the permissions of `configure.sh` before you can execute it. You will know the configuration stage ran successfully if a `config.h` file was generated in the working directory.

If the configuration stage runs successfully, you can compile and run codes in the `runs` folder. Simply,
```
cd runs
make <name-of-code>
./<name-of-code> <inputs>
```
or
```
mpirun -np <N> <name-of-code> <inputs>
```
where `<name-of-code>` is the name of the C++ file containing the program you wish to compile and execute (without the `.cc` extension), `<N>` is the number of processors to use when running the code in parallel (with MPI), and `<inputs>` are any arguments fed to the programs from the terminal (separated by a space).

For convenience, the script `compile_code.sh` compiles the relevant code and generates/cleans the folders in the `results` directory that will house the simulation output data. It can be used by
```
./compile_code.sh <specie>
```
where `<specie>` is either `helium`, `hydrogen`, or `electrons`, depending on which specie you want to simulate.

**Typical simulation pipeline**
This section assumes you have compiled all the relevant codes according to the procedure described previously.
To run a modulation simulation, run the code named `<specie>_main_sim_HS_modulation_parker` where `<specie>` is either `helium`, `hydrogen`, or `electrons`.
More precisely,
```
mpirun -np <N> <specie>_main_sim_HS_modulation_parker <number-of-trajectories> <batch-size> <position-percent>
```
where `<number-of-trajectories>` is the total number of trajectories to average over during this run, `<batch-size>` is the number of trajectories assigned per batch as each processor becomes available to perform work, and `<position-percent>` is an integer between 0 and 100 representing the V2 position at which to compute the modulated spectrum of GCRs as a percentage of distance from a location just before the TS to the HP.

To postprocess the results of this simulation, runs the following codes
```
./<specie>_main_process_HS_modulation_spectrum <position-percent>
./<specie>_main_process_HS_modulation_time <position-percent>
```
This generates the files `results/HS_mod_spec_<sp>/HS_mod_parker_<position-percent>_*.dat`, where `<sp>` is either `He`, `H`, or `e`.
The final postprocessing step is to integrate the spectra at each position, folding them through the detector's response function, to obtain a simulated counts per second at each position
Assuming that `M` modulation simulations have been performed and postprocessed at equally spaced values within [0,100] for the parameter `<position-percent>`, this is done with the command
```
<specie>_main_integrate_HS_modulation_spectra_along_path M
```
The output of this step is in `results/HS_mod_spec_<sp>/HS_mod_parker_integ_spec.dat`.
After running these codes, the results can be visualized using
```
python plot_rate.py <specie>
```
In addition,
```
python plot_modulated_flux.py <specie>
```
can plot the modulated flux at a particular location, specified by the `perc` variable in the code (line 81).

**Parameters File**
The parameters of the modulation simulations are in the corresponding `params_<sp>.txt` files (units in parentheses).

 - *Parameter 1*: Amplitude of the solar cycle variation of the magnetic field at 1 au (G).
 - *Parameter 2*: parallel mean free path in sectored region for a 10 GV particle near 1 au (au).
 - *Parameter 3*: perpendicular mean free path in sectored region for a 10 GV particle near 1 au (au).
 - *Parameter 4*: reduction factor of perpendicular mean free path in unipolar region (unitless).
 - *Parameter 5*: factor controlling amplitude of solar cycle dependence of diffusion coefficients (unitless). This coefficient should be greater than 0. A value of 1 means no solar cycle dependance, whereas a value greater/smaller than 1 gives a proportional/inverse dependence with solar cycle (higher/lower diffusion during solar maximum). See `DiffusionEmpiricalSOQLTandUNLT` source code for details.

**Auxiliary Plotting Routines**
The code `main_plot_solarwind.cc` will plot relevant solar wind quantities along V2's trajectory during its passage through the heliosheath, e.g. magnetic field magnitude, plasma flow speed, heliospheric current sheet latitude, etc.
To run it, execute
```
./main_plot_solarwind <Nt>
```
where `<Nt>` is the number of frames desired (25 by default).
It can also output 2D slices of solar wind quantities on a plane containing V2's trajectory.
To output these, make sure that the macro `PLOT_2D` is defined near the top of the file (uncomment line 10) before compilation, and also configure the code with silo, i.e. add `--enable-silo` to the configuration command in `configure.sh`.
The scripts `plot_V2_and_HCS_lat.py`, `plot_V2_flow_and_field.py`, and `plot_WSO_tilt_angle.py` plot the outputs of this code.

The code `<specie>_main_compute_kappa_<diffusion>.cc` estimates the diffusion coefficient along V2's trajectory for particles of type `<specie>` according to `<diffusion>` theory, which can be `QLT`, `SOQLT`, and `UNLT`.
This uses in-situ observations already downloaded and placed in the `runs/data` folder.
To plot these estimated coefficients, along with the simulated coefficients, use code `<specie>_main_plot_diffusion.cc`.
After running these codes, the results can be visualized using
```
python plot_diffusion.py <specie>
```

The code `main_plot_gcr_drifts.cc` calculates the magnitude of the GCR drifts within a specified window.
The parameters within `params_drifts.txt` specify the particle (specie, energy, etc), background, plotting window, and type of drift to compute.
After running this code, the results can be visualized as 2D color maps using
```
python plot_drifts.py <specie>
```

The remaining Python scripts plot different quantities of interest.
`plot_fit_spectrum_LISM.py` plots the interstellar spectrum used in the simulations, which is fitted to in-situ observations.
`plot_TET_response.py` plots the TET instrument response function.

## Important note

**This is NOT the official SPECTRUM repository.** For information about SPECTRUM, go to https://github.com/vflorins/SPECTRUM.
