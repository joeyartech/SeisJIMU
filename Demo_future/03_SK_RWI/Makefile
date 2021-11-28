n=301
d=10
freq=6

vtrue=1000
#vinit=1000
vinit=1100

model : 
	unisam xin=0,2450,2460,2500 yin=1,1,1,2 nout=$(n) dxout=$(d) | suaddhead ns=$(n) | sustack key=tracl repeat=1 nrepeat=$(n) | sustrip > true_rho
	makevel nx=$(n) nz=$(n) v000=$(vtrue) > true_vp
	
	makevel nx=$(n) nz=$(n) v000=1        > init_rho
	makevel nx=$(n) nz=$(n) v000=1        > init_rhobckg
	makevel nx=$(n) nz=$(n) v000=$(vinit) > init_vp

obs :
	cp setup.in_bkp setup.in
	../../../CPP/SeisJIMU/exe/FWD setup.in
	#
	sustrip < ../03_SK_WPI/results/dsyn_0001.su > temp1 head=head1
	sustrip < ../03_SK_WPI/results/dsyn_0002.su > temp2 head=head2
	supaste < synth_data_0001 ns=1000 head=head1 > dsyn_0001.su
	supaste < synth_data_0002 ns=1000 head=head2 > dsyn_0002.su
	#
	suwind < dsyn_0001.su key=gx min=2000 > faroff_0001.su
	suwind < dsyn_0002.su key=gx max=1500 > nearoff_0001.su
	#
	rm synth_data_000? dsyn_000?.su head? temp?


IpWI :
	echo "0    3000    #offsets of vertices" 	> fwei_ipwi
	echo "0    0       #time instants of vertices" 	>> fwei_ipwi
	echo "0    0       #weight=0 below this polygon" >> fwei_ipwi
	echo "0    3000    #offsets of vertices" >> fwei_ipwi
	echo "1.5    1.5   #time instants of vertices" >> fwei_ipwi
	echo "1    0       #weight=1 below this polygon" >> fwei_ipwi
	echo "END" >> fwei_ipwi
	#
	cp setup.in_bkp setup.in
	echo 'FILE_MODELS	"../init"'  >> setup.in
	echo 'FILE_DATA		"../nearoff_"' >> setup.in
	echo 'ACTIVE_PARAMETER	ip1000:1200' >> setup.in
	echo 'UPDATE_WAVELET	no' >> setup.in
	echo 'FILE_WEIGHT_POLYGON	"../fwei_ipwi"' >> setup.in
	echo 'ZPOWER	0' >> setup.in
	echo 'SMOOTHING	none' >> setup.in
	#echo 'JOB	gradient ' >> setup.in
	#
	-mkdir results_IpWI
	(cd results_IpWI; ../../../../CPP/SeisJIMU/exe/FWI ../setup.in > out)

plot_IpWI :
	sustrip < nearoff_0001.su > temp
	cat temp results_IpWI/synth_data_0001 results_IpWI/residual_0001 | xgraph n1=1000 d1=1 nplot=3 title=obs_nearoff+Ru+residu &
	ximage < results_IpWI/model_update n1=301 perc=99 &
	rm temp

#k=1
k=5
#k=10
JFWI :
	echo "0    1000    #offsets of vertices" > fsepa
	echo "0    1.75    #time instants of vertices" >> fsepa
	echo "END" >> fsepa
	#
	echo "0    1000    #offsets of vertices" 	> fwei_rwi
	echo "0    0       #time instants of vertices" 	>> fwei_rwi
	echo "0    0       #weight=0 below this polygon" >> fwei_rwi
	echo "0    1000    #offsets of vertices" >> fwei_rwi
	echo "1.75 1.75   #time instants of vertices" >> fwei_rwi
	echo "1    0       #weight=1 below this polygon" >> fwei_rwi
	echo "END" >> fwei_rwi
	#
	gfortran fake_init_rho.f90
	./a.out 'results_IpWI/model_update' $(k)
	#
	cp setup.in_bkp setup.in
	echo 'FILE_MODELS	"../init"'    >> setup.in
	echo 'FILE_DATA		"../faroff_"' >> setup.in
	echo 'FILE_SEPARATER    "../fsepa"'  >> setup.in
	echo 'FILE_WEIGHT_POLYGON	"../fwei_rwi"' >> setup.in
	echo 'ACTIVE_PARAMETER	vp1000:1200' >> setup.in
	echo 'UPDATE_WAVELET	no' >> setup.in
	echo 'ZPOWER	0' >> setup.in
	echo 'SMOOTHING	none' >> setup.in
	echo 'JOB	gradient ' >> setup.in
	-mkdir results_JFWI
	(cd results_JFWI; ../../../../CPP/SeisJIMU/exe/RWI ../setup.in > out)
	(cd results_JFWI; cp ../init_rho init_rho-$(k); cp cb%gradient cb%gradient_init_rho-$(k))
	

plot_JFWI :
	sustrip < faroff_0001.su >temp
	cat temp results_JFWI/synth_data_0001 results_JFWI/synth_data_bckg_0001 | xgraph title=obs_faroff+Ru+Ru0-$(k) n1=1000 d1=0.004 nplot=3 &
	cat results_JFWI/residual_0001 results_JFWI/residual2_0001 | xgraph title=Residu1+2-$(k) n1=1000 d1=0.004 nplot=2 &
	ximage < results_JFWI/cb%gradient_init_rho-$(k) n1=301 perc=95 title=cb%gradient_init_rho-$(k) &

movie_JFWI :
	xmovie < results/snap_fld_du%p n1=349 n2=349 clip=0.006 loop=2 title=%g
