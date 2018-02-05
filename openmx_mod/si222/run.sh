#! /usr/bin/bash

alm=si
openmx=si222
mag=0.02
temperature=300

# generate pattern file
../alm/alm ${alm}_alm.in > ${alm_alm}.log1
#../alm/alm ${alm}_alm2.in > ${alm}_alm2.log1

# generate openmx input file for all displacement patterns
python ../tools/displace.py --OpenMX=${openmx}.dat --mag=${mag} ${openmx}.pattern_HARMONIC
#python ../tools/displace.py --OpenMX=${openmx}.dat --mag=${mag} ${openmx}_cubic.pattern_ANHARM3
mv disp*.dat disp_cubic


python ../tools/extract.py --OpenMX=${openmx}.dat --get=disp md/*.md > disp.dat
python ../tools/extract.py --OpenMX=${openmx}.dat --get=force md/*.md > force.dat

#python ../tools/extract.py --OpenMX=${openmx}.dat --get=disp md_cubic/*.md > disp3.dat
#python ../tools/extract.py --OpenMX=${openmx}.dat --get=force md_cubic/*.md > force3.dat


# generate fcs file and result file
../alm/alm ${alm}_alm.in > ${alm}_alm.log2
#../alm/alm ${alm}_alm2.in > ${alm}_alm2.log2

# calculate band and DOS
#../anphon/anphon ${alm}_phband.in > ${alm}_phband.log
#../anphon/anphon ${alm}_phdos.in > ${alm}_phdos.log

# calcualte thermal conductivity
#../anphon/anphon ${alm}_RTA.in > ${alm}_RTA.log

# calcualte life time (tau)
#python ../tools/analyze_phonons.y --calc tau --temp ${temperature} ${openmx}.result > tau${temperature}K.dat

