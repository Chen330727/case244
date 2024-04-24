#!/bin/bash
#set -x
#PBS -l select=1:ncpus=24:mpiprocs=24:mem=64gb
#PBS -l walltime=24:00:00
#PBS -N MESUPROJECTNAME
#PBS -M xiangbin.chen@sorbonne-universite.fr
#PBS -j oe  
# load modules 
#module load openMPI/4.1.2
module load mpt
# module load mpt/2.18
#module load gcc/11.2
#module load openMPI/3.1.5
#module load openMPI/1.10.6
Exefilename='MESUPROJECTNAME'
cd $PBS_O_WORKDIR
echo PBS: working directory is $PBS_O_WORKDIR

DIR=./outfacets
if [ -d "$DIR" ];
then
    echo "$DIR directory exists. then delete files inside\n"
    rm -rf outfacets/*
    #mkdir outfacets
else
        echo "$DIR directory does not exist. then create"
        mkdir outfacets
fi

DIR2=./fieldout
if [ -d "$DIR2" ];
then
    echo "$DIR2 directory exists. then delete files inside\n"
    rm -rf fieldout/*
    #mkdir outfacets
else
        echo "$DIR2 directory does not exist. then create"
        mkdir fieldout
fi

mpiexec_mpt -n 24 ./MESUPROJECTNAME -m WALLTIME 2>> log.dat #PBS -l select=8:ncpus=8:mpiprocs=8

#mpiexec_mpt -n 64 ./MESUPROJECTNAME -m WALLTIME 2>> log.dat #PBS -l select=8:ncpus=8:mpiprocs=8
#mpirun -n 16 ./MESUPROJECTNAME # 12 2>> log >> out # -m WALLTIME 2>> log >> out2 #PBS -l pmem=128gb
