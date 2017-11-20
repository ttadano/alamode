#!/bin/bash

#ALAMODE_root = $HOME/src/alamode
ALAMODE_root=$HOME/Work/alamode
echo ${ALAMODE_root}

for ((i = 0; i <= 1000; i = i + 100))
do
${ALAMODE_root}/tools/dfc2 <<TextForInput
  PbTe_harm444.xml
  PbTe_scph_${i}K.xml
  PbTe_scph4-4.scph_fc2_correction
  ${i}
TextForInput
done
