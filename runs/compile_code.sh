#!/bin/bash

# Terminal arguments:
# $1: specie to simulate (helium, hydrogen, or electrons)

# Check input and create directories
if [ "${1}" == "helium" ]; then
   mkdir -p ../results/HS_mod_spec_He
   rm ../results/HS_mod_spec_He/*
elif [ "${1}" == "hydrogen" ]; then
   mkdir -p ../results/HS_mod_spec_H
   rm ../results/HS_mod_spec_H/*
elif [ "${1}" == "electrons" ]; then
   mkdir -p ../results/HS_mod_spec_e
   rm ../results/HS_mod_spec_e/*
else
   echo "Unrecognized argument for specie."
   echo "Use 'helium', 'hydrogen', or 'electrons'."
   exit 1
fi
mkdir -p ../results/solarwind
rm ../results/solarwind/*

# Compile code
make "${1}_main_sim_HS_modulation_parker"
make "${1}_main_process_HS_modulation_spectrum"
make "${1}_main_process_HS_modulation_time"
make "${1}_main_integrate_HS_modulation_spectra_along_path"
make "${1}_main_compute_kappa_QLT"
make "${1}_main_compute_kappa_SOQLT"
make "${1}_main_compute_kappa_UNLT"
make "${1}_main_plot_diffusion"
