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
	cp setup.in_bkp setup.in
	echo 'FILE_MODEL	init'  >> setup.in
	echo 'FILE_DATA_PREFIX	nearoff_' >> setup.in
	echo 'PARAMETER		vp:1000:1200' >> setup.in
	echo 'UPDATE_WAVELET	no' >> setup.in
	echo 'RDT	0.' >> setup.in
	echo 'WEI	polygon:fwei' >> setup.in
	echo 'SMOOTHING	none' >> setup.in
	echo 'JOB	imaging' >> setup.in
	echo 'DIR_OUT   results_imaging' >> setup.in
	../../exe/WPI setup.in
	cp results_imaging/I .
	#suaddhead < results_imaging/I ns=301 | sushw key=dt a=1e6 | suwind key=tracl min=150 max=150 | sumute key=tracl xmute=1 tmute=3.5 | sustack key=tracl repeat=0 nrepeat=301 | sustrip > I
	
gradient :
	cp setup.in_bkp setup.in
	echo 'FILE_MODEL	init'  >> setup.in
	echo 'FILE_DATA_PREFIX	faroff_' >> setup.in
	echo 'PARAMETER		vp:1000:1200' >> setup.in
	echo 'UPDATE_WAVELET	no' >> setup.in
	echo 'RDT	0.' >> setup.in
	echo "I 'I'" >> setup.in
	#echo 'WEI	polygon:fwei' >> setup.in
	echo 'SMOOTHING	none' >> setup.in
	echo 'JOB	gradient' >> setup.in
	../../exe/WPI setup.in
	
