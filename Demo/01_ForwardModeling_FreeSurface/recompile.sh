#!/bin/bash

cd ../../

sed -i 's/WaveEq=DAS/WaveEq=PSV/' modules.inc
#cat modules.inc
(cd Modules/Modeling; make clean)
(cd FWD; rm *.o)
make fwd

sed -i 's/WaveEq=PSV/WaveEq=DAS/' modules.inc
#cat modules.inc
(cd Modules/Modeling; make clean)
(cd FWD; rm *.o)
make fwd

cd -
