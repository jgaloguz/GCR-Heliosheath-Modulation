#!/bin/bash

if [ "$1" == "helium" ]; then
   mkdir -p ../results/HS_mod_spec_He
   rm ../results/HS_mod_spec_He/*
elif [ "$1" == "electrons" ]; then
   mkdir -p ../results/HS_mod_spec_e
   rm ../results/HS_mod_spec_e/*
else
   echo "Unrecognized argument for specie."
   echo "Use 'helium' or 'electrons'."
   exit 1
fi

make "main_sim_HS_modulation_parker_$1"
make "main_process_HS_modulation_spectrum_$1"
make "main_process_HS_modulation_time_$1"
make "main_integrate_HS_modulation_spectra_along_path_$1"