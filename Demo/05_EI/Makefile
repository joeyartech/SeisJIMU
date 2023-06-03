


fwd : ../../exe/fwd_BornEnv_reconstruction_FD_O2  setup_default
	makevel nz=50 nx=101 v000=1500 > c1
	makevel nz=1  nx=101 v000=2000 > c2
	cat c1 c2 c1 > tmp && transp < tmp n1=101 > model
	rm c1 c2 tmp
	#
	cp setup_default setup
	echo 'DIR_OUT   results_fwd' >> setup
	$< setup

pfei_wave : ../../exe/pfei_BornEnv_reconstruction_FD_O2_Wolfe_SD  setup_default
	suphase < results_fwd/wavelet.su a=-90 > ricker_shift.su
	#
	#makevel nz=50 nx=101 v000=1600 > c1
	#makevel nz=1  nx=101 v000=1900 > c2
	#cat c1 c2 c1 > tmp && transp < tmp n1=101 > model
	#rm c1 c2 tmp
	#
	makevel nz=101 nx=101 v000=1600 > model
	#
	cp setup_default setup
	echo 'DIR_OUT   results_wave' >> setup
	echo 'MODE   wave ' >> setup
	echo 'FILE_WAVELET_ENV   ricker_shift.su' >> setup
	$< setup > out 2>&1

plot_wave :
	suximage < results_wave/dfault_Shot0001.su perc=99 legend=1 &
	suximage < results_wave/RdE_Shot0001.su    perc=99 legend=1 title=RdE &
	#
	ximage < results_wave/snap_surD_rev n1=101  perc=99 legend=1 title=surD_rev &
	ximage < results_wave/snap_fld_dE%p_fwdscat n1=149  perc=99 legend=1 title=dE_fwdscat &
	ximage < results_wave/snap_fld_F0%p_adj n1=149  perc=99 legend=1 title=F0_adj &
	ximage < results_wave/snap_fld_dF%p_adj n1=149  clip=3e10 legend=1 title=dF_adj &
	#
	ximage < results_wave/F0_star_E0 n1=101 perc=99 legend=1 title=F0_star_E0 &
	ximage < results_wave/F0_star_dE n1=101 perc=99 legend=1 title=F0_star_dE &
	ximage < results_wave/dF_star_E0 n1=101 perc=99 legend=1 title=dF_star_E0 &

pfei_env : ../../exe/pfei_BornEnv_reconstruction_FD_O2_Wolfe_SD  setup_default results_fwd/dsyn_Shot0001.su
	suenv < results_fwd/dsyn_Shot0001.su | sushw key=trid a=11 > denv_Shot0001.su
	suenv < results_fwd/wavelet.su | sushw key=trid a=11 |  suphase a=-45 | sumute key=trid mode=1  xmute=1 tmute=1 ntaper=50 > ricker_env_shift.su
	#
	#makevel nz=50 nx=101 v000=1600 > c1
	#makevel nz=1  nx=101 v000=1900 > c2
	#cat c1 c2 c1 > tmp && transp < tmp n1=101 > model
	#rm c1 c2 tmp
	#
	makevel nz=101 nx=101 v000=1600 > model
	#
	cp setup_default setup
	echo 'DIR_OUT   results_env' >> setup
	echo 'MODE   env ' >> setup
	echo 'FILE_DATA_PREFIX   denv_Shot' >> setup
	echo 'FILE_WAVELET_ENV   ricker_env_shift.su' >> setup
	$< setup > out 2>&1
	#
	suximage < results_env/dobs-RE0_Shot0001.su perc=99 legend=1 &
	suximage < results_env/RdE_fwd_Shot0001.su perc=99 legend=1 & #note the small artefacts near source


cleanall :
	-rm setup
	-rm -r results*
