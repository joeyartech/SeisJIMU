#!/bin/bash

n=101

# # FWD ##
# makevel nz=51 nx=$n v000=1800 > vp1
# makevel nz=50 nx=$n v000=1500 > vp2
# cat vp1  vp2  > tmp1 && transp < tmp1 n1=$n > model
# rm vp1 vp2 rho1 rho2 tmp*
#
# ../../exe/FWD setup.in > out_fwd
#
# rm -r results_fwd; mv results  results_fwd
# (cd results_fwd; suenv < dsyn_Shot0001.su > denv_Shot0001.su)
#
#
# # FWD no Refl ##
# makevel nz=$n nx=$n v000=1800 > model
#
# ../../exe/FWD setup.in > out_fwd_norefl
#
# rm -r results_fwd_norefl; mv results  results_fwd_norefl
# (cd results_fwd_norefl; suenv < dsyn_Shot0001.su > denv_Shot0001.su)


# ## PFEI estimate wavelet for envelope modeling ##
# cp setup.in setup_wl.in
# echo "WAVELET_TYPE   'ricker envelope'" >> setup_wl.in
#echo "FILE_DATA_PREFIX    'results_fwd_norefl/denv_Shot'" >> setup_grad.in
# echo "JOB    'estimate wavelet'" >> setup_wl.in
#
# ../../exe/PFEI setup_wl.in > out_wl
# rm -r results_wl; mv results results_wl
# (cd results_wl; sumute < updated_wavelet.su key=tracl xmute=1   ntaper=50 mode=1  tmute=1.8 > muted_updated_wavelet.su)



# ## PFEI gradient ##
# makevel nz=51 nx=$n v000=1800 > vp1
# #makevel nz=50 nx=$n v000=1500 > vp2
# makevel nz=50 nx=$n v000=1800 > vp2
# cat vp1  vp2 > tmp1 && transp < tmp1 n1=$n > model
#
# # makevel nz=$n nx=$n v000=1750 > model
# #
# # makevel nz=49 nx=$n v000=1750 > vp1
# # makevel nz=52 nx=$n v000=1500 > vp2
# # cat vp1  vp2 > tmp1 && transp < tmp1 n1=$n > model
#
# makevel nz=1 nx=$n v000=100 > topo
# #makevel nz=1 nx=$n v000=400 > topo
#
# rm vp1 vp2 rho1 rho2 tmp*
#
# cp setup.in setup_grad.in
# echo "FILE_WAVELET   './results_wl/muted_updated_wavelet.su'" >> setup_grad.in
#
# ../../exe/PFEI setup_grad.in  > out_grad
# rm -r results_grad; mv results results_grad
#
# suximage < results_grad/RE0_Shot0001.su clip=1e-5 &
# suximage < results_grad/RdE_Shot0001.su clip=1e-5 &
# susum  results_grad/RE0_Shot0001.su  results_grad/RdE_Shot0001.su | suximage clip=1e-5 &
#
# ximage < results_grad/gvp2_F1_star_E0 title=gvp2_F1_star_E0 n1=101 clip=1e-19 &
# ximage < results_grad/gvp2_F2_star_dE title=gvp2_F2_star_dE n1=101 clip=1e-19 &
# ximage < results_grad/gvp2_F2_star_E0 title=gvp2_F2_star_E0 n1=101 clip=1e-19 &



## GradientTest ##
cp setup.in setup_test.in
echo "FILE_WAVELET   './results_wl/muted_updated_wavelet.su'" >> setup_test.in
echo "JOB   'build tilD'" >> setup_test.in

rm -r results
../../exe/GradientTest  setup_test.in > out

echo '            alpha    pert%f    curr%f    (pert%f-curr%f)/alpha    curr%g_dot_d    if_1st_cond'
grep '1st cond' out


