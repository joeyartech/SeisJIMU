#!/bin/bash -x

#compile
#(cd ../FWD; make)


run() {
        echo "MODEL_SIZE              '$1    $1   1' " >  setup.in
        echo "MODEL_SPACING           '20    20   1' " >> setup.in
        echo "FILE_MODEL              model " >> setup.in
        echo "MODEL_ATTRIBUTES        vp " >> setup.in
        echo " " >> setup.in
        echo "IS_FREESURFACE          F " >> setup.in
        echo "IF_BLOOM                F " >> setup.in
        echo "IF_HICKS                  F " >> setup.in
        echo " " >> setup.in
        echo "NUMBER_SOURCE       $2 " >> setup.in
        echo " " >> setup.in
        echo "NT           $4 " >> setup.in
        echo "DT           0.004 " >> setup.in
        echo "FPEAK        7 " >> setup.in
        echo " " >> setup.in
        echo "IF_USE_CHECKPOINT   F " >> setup.in

        echo "Writing setup.in OK"

        mpirun -np $3 ../../exe/FWD setup.in
}


#single OMP thread
export OMP_NUM_THREADS=1

#repeating times
REPEAT=5

#number of time steps
NT=5000

gfortran gen_vel.f90

for N in 10 20 40 80 160 320 640 1280 2560; do
#N=1280
for NPROC in 1 2 4 8; do
#NPROC=8
  echo "##################################################"
  ./a.out $N
  run $N $((NPROC*REPEAT)) $NPROC $NT
  echo "##################################################"
done
done

