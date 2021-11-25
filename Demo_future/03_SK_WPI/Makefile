n=301
d=10
freq=6

vtrue=1000
#vinit=1000
vinit=1100

model : 
	unisam xin=0,2450,2460,2500 yin=1000,1000,1000,2000 nout=$(n) dxout=$(d) | suaddhead ns=$(n) | sustack key=tracl repeat=1 nrepeat=$(n) | sustrip > rho
	makevel nx=$(n) nz=$(n) v000=$(vtrue) > vp
	cat vp rho > true
	#
	makevel nx=$(n) nz=$(n) v000=$(vinit) > vp
	makevel nx=$(n) nz=$(n) v000=1000     > rho
	cat vp rho > init
	#
	rm vp rho

obs :
	cp setup.in_bkp setup.in
	../../exe/FWD setup.in
	suwind < results/dsyn_0001.su key=gx min=2000 > faroff_0001.su
	suwind < results/dsyn_0002.su key=gx max=1500 > nearoff_0001.su

image :
	echo "0    3000    #offsets of vertices" > fwei
	#echo "0    0       #time instants of vertices" >> fwei
	#echo "0    0       #weight=0 below this polygon" >> fwei
	#echo "0    3000    #offsets of vertices" >> fwei
	echo "1    1       #time instants of vertices" >> fwei
	echo "1    0       #weight=1 below this polygon" >> fwei
	echo "END" >> fwei
	#
	cp setup.in_bkp setup.in
	echo 'FILE_MODEL	init'  >> setup.in
	echo 'FILE_DATA_PREFIX	nearoff_' >> setup.in
	echo 'PARAMETER		vp:1000:1200' >> setup.in
	echo 'UPDATE_WAVELET	no' >> setup.in
	echo 'RDT	0.' >> setup.in
	echo 'WEI_IMAGING	polygon:fwei' >> setup.in
	echo 'SMOOTHING	none' >> setup.in
	echo 'JOB	imaging' >> setup.in
	echo 'DIR_OUT   "./results_imaging/"' >> setup.in
	../../exe/WPI setup.in
	cp results_imaging/I .
	#suaddhead < results_imaging/I ns=301 | sushw key=dt a=1e6 | suwind key=tracl min=150 max=150 | sumute key=tracl xmute=1 tmute=3.5 | sustack key=tracl repeat=0 nrepeat=301 | sustrip > I

plot_image :
	cat nearoff_0001.su  results_imaging/Ru_imag_0001.su | suxwigb title=obs_nearoff+Ru &
	ximage < results_imaging/I n1=301 perc=99 &

image_multioff :
	#fake multioffset I
	suaddhead < results_imaging/I ns=$(n) | sustack key=cdp repeat=1 nrepeat=$(n) | sustrip > multioff_I
	
gradient :
	cp setup.in_bkp setup.in
	echo 'FILE_MODEL	init'  >> setup.in
	echo 'FILE_DATA_PREFIX	faroff_' >> setup.in
	echo 'PARAMETER		vp:1000:1200' >> setup.in
	echo 'UPDATE_WAVELET	no' >> setup.in
	echo 'RDT	0.' >> setup.in
	#echo "IMAGE 'I'" >> setup.in
	echo "IMAGE 'multioff_I'" >> setup.in
	echo 'SMOOTHING	none' >> setup.in
	echo 'JOB	gradient' >> setup.in
	../../exe/WPI setup.in

plot_gradient :
	cat faroff_0001.su  results/Ru_0001.su  results/Rdu_0001.su | suxwigb title=obs_faroff+Ru+Rdu &
	ximage < results/u_star_da n1=301 perc=99 title=u_star_da &
	ximage < results/u_star_a n1=301 perc=99 title=u_star_a &

movie_gradient :
	xmovie < results/snap_fld_du%p n1=349 n2=349 clip=0.006 loop=2 title=%g
