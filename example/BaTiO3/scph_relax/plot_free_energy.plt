set terminal pdf size 4,2.8
set output "BaTiO3_scph_relax_free_energy.pdf"

set grid
set xlabel "Temperature [K]"
set ylabel "Free energy [Ry]"

plot "cBTO222_scph.scph_thermo" using 1:6 title "total" with lines,\
"cBTO222_scph.scph_thermo" using 1:5 title "U_0" with lines,\
"cBTO222_scph.scph_thermo" using 1:($3+$4) title "F_{vib}" with lines,\
