#!/bin/bash
#set -x
module purge
module load mpt
# module load gcc/11.2
# module load gcc/8.2
#module load openMPI/1.10.6
# module load mpt/2.18-ga
# module load mpt/2.18
# module load openMPI/4.1.2-gcc82
# module load mpi
#module load openMPI/3.1.5
#module load openMPI/1.10.6
echo PBS: working directory is $PWD
Exefilename='MESUPROJECTNAME'
Cfile='_template01.c'
[ -f ./$Cfile ] && echo "Found" || echo "Not found"
#PATH=$PATH:/home/basilisk/bin
#export PATH="$HOME:$PATH"
which mpicc
# mpicc -Wall -O2 -std=c99 $Cfile -o $Exefilename -I$HOME -L$HOME/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lm

# mpicc -O2 -std=c99 -Wl,--verbose $Cfile -o $Exefilename -L$HOME/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lm

# mpicc -Wall -O2 -std=c99 -D_XOPEN_SOURCE=700 $Cfile  -o $Exefilename \
#         -L$HOME/gl -lglutils -lfb_tiny -lm

# mpicc -Wall -O2 -std=c99 -D_XOPEN_SOURCE=700 $Cfile  -o $Exefilename \
#         -L$HOME/gl -lglutils -lm

mpicc -Wall -O2 -std=c99 -D_XOPEN_SOURCE=700 $Cfile  -o $Exefilename \
        -I$HOME -L$HOME/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lm
# ldd $Exefilename
qsub run2-sub.sh
#cd $PBS_O_WORKDIR
#mpiexec_mpt -n 64 ./$Exefilename -m WALLTIME 2>> log >> out2
