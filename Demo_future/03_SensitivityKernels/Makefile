n=301
d=10
freq=6

vtrue=1000
#vinit=1000
vinit=1100

model : 
	unisam xin=0,2490,2500,2510 yin=1000,1000,2000,1000 nout=$(n) dxout=$(d) | suaddhead ns=$(n) | sustack key=tracl repeat=1 nrepeat=$(n) | sustrip > rho
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
	mv results/dsyn_0001.su obs_0001.su

gradient :
	cp setup.in_bkp setup.in
	echo 'FILE_MODEL	init'  >> setup.in
	echo 'FILE_DATA_PREFIX	obs_' >> setup.in
	echo 'PARAMETER		vp:1000:1200' >> setup.in
	echo 'UPDATE_WAVELET	no' >> setup.in
	echo 'RDT	0.' >> setup.in
	echo 'SMOOTHING	none' >> setup.in
	echo 'JOB	gradient' >> setup.in
	../../exe/WPI setup.in
	
