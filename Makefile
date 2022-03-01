FORTRAN=gfortran
FORTRAN_FLAGS= -O3 -fopenmp
SUB_FLAGS=-c

UTILS= mod_utils.f90
PROCESSES= mod_IC.f90 mod_synchrotron.f90 mod_particle_distribution.f90 mod_pair_production.f90

MODULES=$(subst f90,o,$(UTILS) $(PROCESSES)) 

#optical data from RS Oph 
OPTICAL_DATA= observations/RSOph_Bmag.dat observations/RSOph_Vmag.dat observations/RSOph_Rmag.dat observations/RSOph_Imag.dat
HESS_DATA=  observations/sed_stereo_090821_bias2p5.dat observations/sed_stereo_100821_bias2p5.dat observations/sed_stereo_110821_bias2p5.dat observations/sed_stereo_120821_bias2p5.dat observations/sed_stereo_130821_bias2p5.dat
HESSII_DATA= observations/sed_mono_090821_low.dat observations/sed_mono_100821_low.dat observations/sed_mono_110821_low.dat observations/sed_mono_120821_low.dat observations/sed_mono_130821_low_noshow.dat 
HE_DATA= observations/sed_LAT_090821_HC_REBIN.dat observations/sed_LAT_100821_HC.dat observations/sed_LAT_110821_HC.dat observations/sed_LAT_120821_HC.dat observations/sed_LAT_130821_HC_REBIN.dat


figures: dat/table_input.tex pp_emission.jpg ic_emission.jpg

%.jpg: %.eps
	convert -density 300x $^ $@

%.o: %.f90
	$(FORTRAN) -c $^

pp_emission.eps ic_emission.eps: emission.gnuplot multi_pannel.cfg add_data.cfg dat/pp_examples.dat dat/ic_examples.dat
	gnuplot emission.gnuplot

dat/ic_examples.dat  dat/pp_examples.dat dat/gammagamma_examples.dat dat/gammagamma_pp_examples.dat dat/gammagamma_ic_examples.dat: dat/physical_conditions_unform.dat dat/protons.dat dat/electrons.dat module_pp.f90 nonthermal_emission.f90
	$(FORTRAN) -O2 $(UTILS) $(PROCESSES) myparameters.f90 setup.f90 module_pp.f90 nonthermal_emission.f90 -o nonthermal_emission.out
	./nonthermal_emission.out

dat/table_input.tex dat/physical_conditions.dat dat/physical_conditions_unform.dat dat/protons.dat dat/electrons.dat dat/electrons_examples.dat dat/protons_examples.dat: $(MODULES) myparameters.f90 setup.f90 nonthermal_particles.f90 $(OPTICAL_DATA)
	mkdir -p dat/
	$(FORTRAN) -O2 $(UTILS) $(PROCESSES) myparameters.f90 setup.f90 nonthermal_particles.f90 -o nonthermal_particles.out
	./nonthermal_particles.out

