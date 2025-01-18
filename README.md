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
./<name-of-code> <input>
```
or
```
mpirun -np <N> <name-of-code> <input>
```
where `<name-of-code>` is the name of the C++ file containing the program you wish to compile and execute (without the `.cc` extension), `<N>` is the number of processors to use when running the code in parallel (with MPI), and `<input>` are any arguments fed to the programs from the terminal (separated by a space).

## Important note

**This is NOT the official SPECTRUM repository.** For information about SPECTRUM, go to https://github.com/vflorins/SPECTRUM.
