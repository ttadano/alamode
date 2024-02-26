set logscale x
set format x "%.1e"
set grid

set terminal png
set output "BTO_IFC_cv.png"

set xlabel "CV alpha"
set ylabel "error"

plot "cBTO222.cvscore" using 1:2 title "Fitting error (mean)" with linespoints,\
"cBTO222.cvscore" using 1:4 title "Validation error (mean)" with linespoints


