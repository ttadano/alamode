set terminal pdf size 4,2.8
set output "BaTiO3_scph_relax.pdf"

set grid
set xlabel "Temperature [K]"
set ylabel "atomic displacements [Bohr]"

plot "cBTO222_scph.atom_disp" using 1:4 title "Ba(z)" with lines,\
"cBTO222_scph.atom_disp" using 1:7 title "Ti(z)" with lines,\
"cBTO222_scph.atom_disp" using 1:10 title "O(1,2,z)" with lines,\
"cBTO222_scph.atom_disp" using 1:16 title "O(3,z)" with lines,\