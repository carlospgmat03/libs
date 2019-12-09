unset label
unset xrange
unset yrange

set xlabel "Purity"
set ylabel "Concurrence"

set xrange [0.25:1.0]
set yrange [0.0:1.0]

set key left top

set xlabel font ", 12"
set ylabel font ", 12"

plot '../../tests/q13_t4500/test_bpe1.4_bpa1.4_1.dat' using 2:3 with l title 'bx=bz=1.4' lc 'green' lw 2, \
     '../../tests/q13_t4500/test_bpe1.55_bpa0.0_1.dat' using 2:3 with l title 'bx=1.55,bz=0' lc 'red', \
     '../../tests/q13_t4500/test_bpe1.88_bpa0.58_1.dat' using 2:3 with l title 'bx=1.88,bz=0.58' lc 'blue'
