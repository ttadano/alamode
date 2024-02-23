set logscale x
set grid

set terminal pdf
set output "cv.pdf"

set xlabel "CV alpha"
set ylabel "error"

plot "cBTO222.cvscore" using 1:2 title "Fitting error (mean)" with linespoints,\
"cBTO222.cvscore" using 1:4 title "Validation error (mean)" with linespoints


