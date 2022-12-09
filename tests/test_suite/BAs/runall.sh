#!/bin/sh
#QSUB2 queue qM
#QSUB2 core 192 
#QSUB2 mpi 32
#QSUB2 smp 6
#QSUB2 wtime 48:00:00
#PBS -N whzjob
cd $PBS_O_WORKDIR
source /etc/profile.d/modules.sh

module load comp1
module load intel-mkl/20.0.0

export PATH=$PATH:/home/whzhang/work/packages/qe-6.4.1/bin
export PATH=$PATH:/home/whzhang/work/packages/my_bin
export LD_LIBRARY_PATH=/home/whzhang/work/packages/alamode/spglib/lib:$LD_LIBRARY_PATH

anphon=/home/whzhang/work/calculations/4ph/codes/cccad8/anphon

jbs=(pure_rta pure_3ph pure_4ph natu_rta natu_3ph natu_4ph)

for folder in "${jbs[@]}" 
do
  cd $folder
  ln -s ../BAs-paw-40-480-N60-rcut12-8.xml BAs.xml
  ln -s ../BAs.born .
  mpijob $anphon anphon.in > anphon.out
  rm BAs.born
  rm BAs.xml
  cd - 
done
