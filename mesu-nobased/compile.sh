#!/bin/bash
#set -x

tpvalue=(15)  # Update with multiple values if needed
maxlevelvalue=(11)  # Update with multiple values if needed
Rbubvalue=(30)  # Update with multiple values if needed

cd ../
for tp in "${tpvalue[@]}"; do
   for rbubv in "${Rbubvalue[@]}"; do
      for maxlval in "${maxlevelvalue[@]}"; do
          rm -f DDestIN DestIN.c DestIN
          StandardIN="template01-template-2-2-restart-2-noT.c"
          DestIN='compile'
          sed s/CONTACTANGLE/$tp/g $StandardIN | sed s/MESHN/$maxlval/g | sed s/JACOBNUMBER/$rbubv/ > $DestIN.c
        #   mpicc -std=c99 -O2 -Wall -D_MPI=1 -events $DestIN.c -o $DestIN -L$BASILISK/gl/ -lOSMesa -lglutils -lfb_osmesa -lGLU -lm
          CC99='mpicc -std=c99' qcc -O2 -Wall -D_MPI=1 -events $DestIN.c  -o $DestIN -L$BASILISK/gl/ -lOSMesa -lglutils -lfb_osmesa -lGLU -lm 
          
          if [ $? -eq 0 ]; then
              echo -e "\033[0;32mCompile success\033[0m"
          else
              echo -e "\033[0;31mCompile failed\033[0m"
          fi
      done
   done
done
