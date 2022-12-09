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

jbs=(natu_rta natu_3ph natu_4ph)

for folder in "${jbs[@]}" 
do
  cd $folder
  ln -s ../AlSb.100DFSET.xml AlSb.xml
  ln -s ../AlSb.born .
  mpijob $anphon anphon.in > anphon.out
  rm AlSb.born
  rm AlSb.xml
  cd - 
done
