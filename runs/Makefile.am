## Process this file with automake to produce Makefile.in

AM_CXXFLAGS = $(MPI_CFLAGS)
AM_LDFLAGS = $(MPI_LIBS)

bin_PROGRAMS = helium_main_sim_HS_modulation_parker \
               helium_main_process_HS_modulation_spectrum \
               helium_main_process_HS_modulation_time \
               helium_main_integrate_HS_modulation_spectra_along_path \
               helium_main_compute_kappa_QLT \
               helium_main_compute_kappa_SOQLT \
               helium_main_compute_kappa_UNLT \
               helium_main_plot_diffusion \
               hydrogen_main_sim_HS_modulation_parker \
               hydrogen_main_process_HS_modulation_spectrum \
               hydrogen_main_process_HS_modulation_time \
               hydrogen_main_integrate_HS_modulation_spectra_along_path \
               hydrogen_main_compute_kappa_QLT \
               hydrogen_main_compute_kappa_SOQLT \
               hydrogen_main_compute_kappa_UNLT \
               hydrogen_main_plot_diffusion \
               electrons_main_plot_diffusion \
               electrons_main_compute_kappa_UNLT \
               electrons_main_compute_kappa_SOQLT \
               electrons_main_compute_kappa_QLT \
               electrons_main_process_HS_modulation_time \
               electrons_main_process_HS_modulation_spectrum \
               electrons_main_sim_HS_modulation_parker \
               electrons_main_integrate_HS_modulation_spectra_along_path \
               main_plot_solarwind main_plot_gcr_drifts

SPBL_COMMON_DIR = ../common
SPBL_SOURCE_DIR = ../src

helium_main_sim_HS_modulation_parker_SOURCES = helium_main_sim_HS_modulation_parker.cc \
   $(SPBL_SOURCE_DIR)/simulation.cc \
   $(SPBL_SOURCE_DIR)/simulation.hh \
   $(SPBL_SOURCE_DIR)/trajectory_parker.cc \
   $(SPBL_SOURCE_DIR)/trajectory_parker.hh \
   $(SPBL_SOURCE_DIR)/trajectory_base.cc \
   $(SPBL_SOURCE_DIR)/trajectory_base.hh \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.hh \
   $(SPBL_SOURCE_DIR)/background_solarwind.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/reader_cartesian.cc \
   $(SPBL_SOURCE_DIR)/reader_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/boundary_momentum.cc \
   $(SPBL_SOURCE_DIR)/boundary_momentum.hh \
   $(SPBL_SOURCE_DIR)/boundary_space.cc \
   $(SPBL_SOURCE_DIR)/boundary_space.hh \
   $(SPBL_SOURCE_DIR)/boundary_time.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.hh \
   $(SPBL_SOURCE_DIR)/boundary_base.cc \
   $(SPBL_SOURCE_DIR)/boundary_base.hh \
   $(SPBL_SOURCE_DIR)/initial_time.cc \
   $(SPBL_SOURCE_DIR)/initial_time.hh \
   $(SPBL_SOURCE_DIR)/initial_space.cc \
   $(SPBL_SOURCE_DIR)/initial_space.hh \
   $(SPBL_SOURCE_DIR)/initial_momentum.cc \
   $(SPBL_SOURCE_DIR)/initial_momentum.hh \
   $(SPBL_SOURCE_DIR)/initial_base.cc \
   $(SPBL_SOURCE_DIR)/initial_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_other.cc \
   $(SPBL_SOURCE_DIR)/diffusion_other.hh \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_SOURCE_DIR)/distribution_other.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.hh \
   $(SPBL_SOURCE_DIR)/distribution_templated.cc \
   $(SPBL_SOURCE_DIR)/distribution_templated.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_SOURCE_DIR)/traj_config.hh \
   $(SPBL_SOURCE_DIR)/server_config.hh \
   $(SPBL_COMMON_DIR)/workload_manager.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

helium_main_sim_HS_modulation_parker_LDADD = $(MPI_LIBS) $(GSL_LIBS)

hydrogen_main_sim_HS_modulation_parker_SOURCES = hydrogen_main_sim_HS_modulation_parker.cc \
   $(SPBL_SOURCE_DIR)/simulation.cc \
   $(SPBL_SOURCE_DIR)/simulation.hh \
   $(SPBL_SOURCE_DIR)/trajectory_parker.cc \
   $(SPBL_SOURCE_DIR)/trajectory_parker.hh \
   $(SPBL_SOURCE_DIR)/trajectory_base.cc \
   $(SPBL_SOURCE_DIR)/trajectory_base.hh \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.hh \
   $(SPBL_SOURCE_DIR)/background_solarwind.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/reader_cartesian.cc \
   $(SPBL_SOURCE_DIR)/reader_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/boundary_momentum.cc \
   $(SPBL_SOURCE_DIR)/boundary_momentum.hh \
   $(SPBL_SOURCE_DIR)/boundary_space.cc \
   $(SPBL_SOURCE_DIR)/boundary_space.hh \
   $(SPBL_SOURCE_DIR)/boundary_time.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.hh \
   $(SPBL_SOURCE_DIR)/boundary_base.cc \
   $(SPBL_SOURCE_DIR)/boundary_base.hh \
   $(SPBL_SOURCE_DIR)/initial_time.cc \
   $(SPBL_SOURCE_DIR)/initial_time.hh \
   $(SPBL_SOURCE_DIR)/initial_space.cc \
   $(SPBL_SOURCE_DIR)/initial_space.hh \
   $(SPBL_SOURCE_DIR)/initial_momentum.cc \
   $(SPBL_SOURCE_DIR)/initial_momentum.hh \
   $(SPBL_SOURCE_DIR)/initial_base.cc \
   $(SPBL_SOURCE_DIR)/initial_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_other.cc \
   $(SPBL_SOURCE_DIR)/diffusion_other.hh \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_SOURCE_DIR)/distribution_other.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.hh \
   $(SPBL_SOURCE_DIR)/distribution_templated.cc \
   $(SPBL_SOURCE_DIR)/distribution_templated.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_SOURCE_DIR)/traj_config.hh \
   $(SPBL_SOURCE_DIR)/server_config.hh \
   $(SPBL_COMMON_DIR)/workload_manager.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

hydrogen_main_sim_HS_modulation_parker_LDADD = $(MPI_LIBS) $(GSL_LIBS)

electrons_main_sim_HS_modulation_parker_SOURCES = electrons_main_sim_HS_modulation_parker.cc \
   $(SPBL_SOURCE_DIR)/simulation.cc \
   $(SPBL_SOURCE_DIR)/simulation.hh \
   $(SPBL_SOURCE_DIR)/trajectory_parker.cc \
   $(SPBL_SOURCE_DIR)/trajectory_parker.hh \
   $(SPBL_SOURCE_DIR)/trajectory_base.cc \
   $(SPBL_SOURCE_DIR)/trajectory_base.hh \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.hh \
   $(SPBL_SOURCE_DIR)/background_solarwind.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/reader_cartesian.cc \
   $(SPBL_SOURCE_DIR)/reader_cartesian.hh \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/boundary_momentum.cc \
   $(SPBL_SOURCE_DIR)/boundary_momentum.hh \
   $(SPBL_SOURCE_DIR)/boundary_space.cc \
   $(SPBL_SOURCE_DIR)/boundary_space.hh \
   $(SPBL_SOURCE_DIR)/boundary_time.cc \
   $(SPBL_SOURCE_DIR)/boundary_time.hh \
   $(SPBL_SOURCE_DIR)/boundary_base.cc \
   $(SPBL_SOURCE_DIR)/boundary_base.hh \
   $(SPBL_SOURCE_DIR)/initial_time.cc \
   $(SPBL_SOURCE_DIR)/initial_time.hh \
   $(SPBL_SOURCE_DIR)/initial_space.cc \
   $(SPBL_SOURCE_DIR)/initial_space.hh \
   $(SPBL_SOURCE_DIR)/initial_momentum.cc \
   $(SPBL_SOURCE_DIR)/initial_momentum.hh \
   $(SPBL_SOURCE_DIR)/initial_base.cc \
   $(SPBL_SOURCE_DIR)/initial_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_other.cc \
   $(SPBL_SOURCE_DIR)/diffusion_other.hh \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_SOURCE_DIR)/distribution_other.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.hh \
   $(SPBL_SOURCE_DIR)/distribution_templated.cc \
   $(SPBL_SOURCE_DIR)/distribution_templated.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_SOURCE_DIR)/traj_config.hh \
   $(SPBL_SOURCE_DIR)/server_config.hh \
   $(SPBL_COMMON_DIR)/workload_manager.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

electrons_main_sim_HS_modulation_parker_LDADD = $(MPI_LIBS) $(GSL_LIBS)

helium_main_process_HS_modulation_spectrum_SOURCES = helium_main_process_HS_modulation_spectrum.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.hh \
   $(SPBL_SOURCE_DIR)/distribution_templated.cc \
   $(SPBL_SOURCE_DIR)/distribution_templated.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

helium_main_process_HS_modulation_spectrum_LDADD = $(MPI_LIBS) $(GSL_LIBS)

hydrogen_main_process_HS_modulation_spectrum_SOURCES = hydrogen_main_process_HS_modulation_spectrum.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.hh \
   $(SPBL_SOURCE_DIR)/distribution_templated.cc \
   $(SPBL_SOURCE_DIR)/distribution_templated.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

hydrogen_main_process_HS_modulation_spectrum_LDADD = $(MPI_LIBS) $(GSL_LIBS)

electrons_main_process_HS_modulation_spectrum_SOURCES = electrons_main_process_HS_modulation_spectrum.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.hh \
   $(SPBL_SOURCE_DIR)/distribution_templated.cc \
   $(SPBL_SOURCE_DIR)/distribution_templated.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

electrons_main_process_HS_modulation_spectrum_LDADD = $(MPI_LIBS) $(GSL_LIBS)

helium_main_process_HS_modulation_time_SOURCES = helium_main_process_HS_modulation_time.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.hh \
   $(SPBL_SOURCE_DIR)/distribution_templated.cc \
   $(SPBL_SOURCE_DIR)/distribution_templated.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

helium_main_process_HS_modulation_time_LDADD = $(MPI_LIBS) $(GSL_LIBS)

hydrogen_main_process_HS_modulation_time_SOURCES = hydrogen_main_process_HS_modulation_time.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.hh \
   $(SPBL_SOURCE_DIR)/distribution_templated.cc \
   $(SPBL_SOURCE_DIR)/distribution_templated.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

hydrogen_main_process_HS_modulation_time_LDADD = $(MPI_LIBS) $(GSL_LIBS)

electrons_main_process_HS_modulation_time_SOURCES = electrons_main_process_HS_modulation_time.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.cc \
   $(SPBL_SOURCE_DIR)/distribution_other.hh \
   $(SPBL_SOURCE_DIR)/distribution_templated.cc \
   $(SPBL_SOURCE_DIR)/distribution_templated.hh \
   $(SPBL_SOURCE_DIR)/distribution_base.cc \
   $(SPBL_SOURCE_DIR)/distribution_base.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

electrons_main_process_HS_modulation_time_LDADD = $(MPI_LIBS) $(GSL_LIBS)

helium_main_integrate_HS_modulation_spectra_along_path_SOURCES = helium_main_integrate_HS_modulation_spectra_along_path.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

helium_main_integrate_HS_modulation_spectra_along_path_LDADD = $(MPI_LIBS) $(GSL_LIBS)

hydrogen_main_integrate_HS_modulation_spectra_along_path_SOURCES = hydrogen_main_integrate_HS_modulation_spectra_along_path.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

hydrogen_main_integrate_HS_modulation_spectra_along_path_LDADD = $(MPI_LIBS) $(GSL_LIBS)

electrons_main_integrate_HS_modulation_spectra_along_path_SOURCES = electrons_main_integrate_HS_modulation_spectra_along_path.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

electrons_main_integrate_HS_modulation_spectra_along_path_LDADD = $(MPI_LIBS) $(GSL_LIBS)

helium_main_compute_kappa_QLT_SOURCES = helium_main_compute_kappa_QLT.cc \
   compute_kappa_common.cc \
   compute_kappa_common.hh \
   compute_kappa_parallel.cc \
   compute_kappa_parallel.hh \
   compute_kappa_parallel_QLT.cc \
   compute_kappa_parallel_QLT.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

helium_main_compute_kappa_QLT_LDADD = $(MPI_LIBS) $(GSL_LIBS)

hydrogen_main_compute_kappa_QLT_SOURCES = hydrogen_main_compute_kappa_QLT.cc \
   compute_kappa_common.cc \
   compute_kappa_common.hh \
   compute_kappa_parallel.cc \
   compute_kappa_parallel.hh \
   compute_kappa_parallel_QLT.cc \
   compute_kappa_parallel_QLT.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

hydrogen_main_compute_kappa_QLT_LDADD = $(MPI_LIBS) $(GSL_LIBS)

electrons_main_compute_kappa_QLT_SOURCES = electrons_main_compute_kappa_QLT.cc \
   compute_kappa_common.cc \
   compute_kappa_common.hh \
   compute_kappa_parallel.cc \
   compute_kappa_parallel.hh \
   compute_kappa_parallel_QLT.cc \
   compute_kappa_parallel_QLT.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

electrons_main_compute_kappa_QLT_LDADD = $(MPI_LIBS) $(GSL_LIBS)

helium_main_compute_kappa_SOQLT_SOURCES = helium_main_compute_kappa_SOQLT.cc \
   compute_kappa_common.cc \
   compute_kappa_common.hh \
   compute_kappa_parallel.cc \
   compute_kappa_parallel.hh \
   compute_kappa_parallel_SOQLT.cc \
   compute_kappa_parallel_SOQLT.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

helium_main_compute_kappa_SOQLT_LDADD = $(MPI_LIBS) $(GSL_LIBS)

hydrogen_main_compute_kappa_SOQLT_SOURCES = hydrogen_main_compute_kappa_SOQLT.cc \
   compute_kappa_common.cc \
   compute_kappa_common.hh \
   compute_kappa_parallel.cc \
   compute_kappa_parallel.hh \
   compute_kappa_parallel_SOQLT.cc \
   compute_kappa_parallel_SOQLT.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

hydrogen_main_compute_kappa_SOQLT_LDADD = $(MPI_LIBS) $(GSL_LIBS)

electrons_main_compute_kappa_SOQLT_SOURCES = electrons_main_compute_kappa_SOQLT.cc \
   compute_kappa_common.cc \
   compute_kappa_common.hh \
   compute_kappa_parallel.cc \
   compute_kappa_parallel.hh \
   compute_kappa_parallel_SOQLT.cc \
   compute_kappa_parallel_SOQLT.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

electrons_main_compute_kappa_SOQLT_LDADD = $(MPI_LIBS) $(GSL_LIBS)

helium_main_compute_kappa_UNLT_SOURCES = helium_main_compute_kappa_UNLT.cc \
   compute_kappa_common.cc \
   compute_kappa_common.hh \
   compute_kappa_perpendicular.cc \
   compute_kappa_perpendicular.hh \
   compute_kappa_perpendicular_UNLT.cc \
   compute_kappa_perpendicular_UNLT.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

helium_main_compute_kappa_UNLT_LDADD = $(MPI_LIBS) $(GSL_LIBS)

hydrogen_main_compute_kappa_UNLT_SOURCES = hydrogen_main_compute_kappa_UNLT.cc \
   compute_kappa_common.cc \
   compute_kappa_common.hh \
   compute_kappa_perpendicular.cc \
   compute_kappa_perpendicular.hh \
   compute_kappa_perpendicular_UNLT.cc \
   compute_kappa_perpendicular_UNLT.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

hydrogen_main_compute_kappa_UNLT_LDADD = $(MPI_LIBS) $(GSL_LIBS)

electrons_main_compute_kappa_UNLT_SOURCES = electrons_main_compute_kappa_UNLT.cc \
   compute_kappa_common.cc \
   compute_kappa_common.hh \
   compute_kappa_perpendicular.cc \
   compute_kappa_perpendicular.hh \
   compute_kappa_perpendicular_UNLT.cc \
   compute_kappa_perpendicular_UNLT.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

electrons_main_compute_kappa_UNLT_LDADD = $(MPI_LIBS) $(GSL_LIBS)

helium_main_plot_diffusion_SOURCES = helium_main_plot_diffusion.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.hh \
   $(SPBL_SOURCE_DIR)/background_solarwind.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_other.cc \
   $(SPBL_SOURCE_DIR)/diffusion_other.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

helium_main_plot_diffusion_LDADD = $(MPI_LIBS) $(GSL_LIBS)

hydrogen_main_plot_diffusion_SOURCES = hydrogen_main_plot_diffusion.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.hh \
   $(SPBL_SOURCE_DIR)/background_solarwind.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_other.cc \
   $(SPBL_SOURCE_DIR)/diffusion_other.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

hydrogen_main_plot_diffusion_LDADD = $(MPI_LIBS) $(GSL_LIBS)

electrons_main_plot_diffusion_SOURCES = electrons_main_plot_diffusion.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.hh \
   $(SPBL_SOURCE_DIR)/background_solarwind.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_other.cc \
   $(SPBL_SOURCE_DIR)/diffusion_other.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

electrons_main_plot_diffusion_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_plot_solarwind_SOURCES = main_plot_solarwind.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.hh \
   $(SPBL_SOURCE_DIR)/background_solarwind.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

main_plot_solarwind_LDADD = $(MPI_LIBS) $(GSL_LIBS)

main_plot_gcr_drifts_SOURCES = main_plot_gcr_drifts.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind_termshock.hh \
   $(SPBL_SOURCE_DIR)/background_solarwind.cc \
   $(SPBL_SOURCE_DIR)/background_solarwind.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base_visual.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_other.cc \
   $(SPBL_SOURCE_DIR)/diffusion_other.hh \
   $(SPBL_SOURCE_DIR)/server_base.cc \
   $(SPBL_SOURCE_DIR)/server_base.hh \
   $(SPBL_SOURCE_DIR)/block_base.cc \
   $(SPBL_SOURCE_DIR)/block_base.hh \
   $(SPBL_SOURCE_DIR)/cache_lru.cc \
   $(SPBL_SOURCE_DIR)/cache_lru.hh \
   $(SPBL_COMMON_DIR)/mpi_config.cc \
   $(SPBL_COMMON_DIR)/mpi_config.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/rk_config.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

main_plot_gcr_drifts_LDADD = $(MPI_LIBS) $(GSL_LIBS)