set terminal pdf
set output "ZnO_thermal_strain.pdf"

set grid
set xlabel "Temperature [K]"
set ylabel "thermal strain"

plot "ZnO_qha.umn_tensor" using 1:2 with lines title "u_{xx} = u_{yy}",\
"ZnO_qha.umn_tensor" using 1:10 with lines title "u_{zz}",\
# "ZnO_qha_zsisa.umn_tensor" using 1:2 with lines title "u_{xx} = u_{yy}, ZSISA",\
"ZnO_qha_zsisa.umn_tensor" using 1:10 with lines title "u_{zz}, ZSISA",\
"ZnO_qha_vzsisa.umn_tensor" using 1:2 with lines title "u_{xx} = u_{yy}, v-ZSISA",\
"ZnO_qha_vzsisa.umn_tensor" using 1:10 with lines title "u_{zz}, v-ZSISA",\
