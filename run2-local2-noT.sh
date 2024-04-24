#!/bin/bash
#set -x

./clean.sh
rm  template01-test
rm -rf dump-i*
rm log.dat
rm -rf log-*
rm -rf dumpfaterdiff3-i*
rm time-performance.dat
rm volume.dat

rm template01-test.c
tpvalue=(15.0) #contact angle
#maxlevelvalue=(256 512 1024 2048) #MESHN
maxlevelvalue=(11) #(7,8,9) #MESHN:maxl-level of grid
Rbubvalue=(30) #(30) #Ja number 
for tp in "${tpvalue[@]}"; 
   do
   for rbubv in "${Rbubvalue[@]}";
    do   
      for maxlval in "${maxlevelvalue[@]}";
        do
        StandardIN=template01-template-2-2-restart-2-noT.c
     sed s/CONTACTANGLE/$tp/g $StandardIN | sed s/MESHN/$maxlval/g | sed s/JACOBNUMBER/$rbubv/ > template01-test.c
    done
done

done


#qcc -O2 -Wall case1-3-direct.c -o case1-3-direct-test1 -L$BASILISK/gl/ -lOSMesa -lglutils -lfb_osmesa -lGLU -lm
#CC99='mpicc -std=c99' qcc -O2 -events -Wall -D_MPI=1 axi_rising_bubble.c -o axi_rising_bubble01 -L$BASILISK/gl/ -lOSMesa -lglutils -lfb_osmesa -lGLU -lm #-lmpi
# CC99='mpicc -std=c99' qcc -O2 -Wall -D_MPI=1 -events template01-test.c  -o template01-test -L$BASILISK/gl/ -lOSMesa -lglutils -lfb_osmesa -lGLU -lm #-lmpi
CC99='mpicc -std=c99' qcc -O2 -Wall -D_MPI=1 template01-test.c  -o template01-test -L$BASILISK/gl/ -lOSMesa -lglutils -lfb_osmesa -lGLU -lm #-lmpi

#qcc -O2 -Wall  heat-test.c -o heat-test1 -L$BASILISK/gl/ -lOSMesa -lglutils -lfb_osmesa -lGLU -lm #-lmpi
#CC99='mpicc -std=c99' qcc -O2 -Wall heat-test.c -o heat-test1 -L$BASILISK/gl/ -lOSMesa -lglutils -lfb_osmesa -lGLU -lm #-lmpi

DIR=./outfacets
if [ -d "$DIR" ];
then
    echo "$DIR directory exists. then delete and re-create\n"
    rm -rf outfacets
    mkdir outfacets
else
	echo "$DIR directory does not exist. then create"
	mkdir outfacets
fi

#./case1-3-direct-test1 > log.dat
#./heat-test1 > log.dat
# mpirun --prefix /usr/bin/mpicc -np 8 ./template01-template-2-2-restart-noT > log.dat
mpirun --host localhost:24 -np 24 ./template01-test > log.dat
#mpirun -np 8 ./heat-test1 > log.dat
