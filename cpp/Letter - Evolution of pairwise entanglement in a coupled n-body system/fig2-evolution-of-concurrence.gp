unset label
unset xrange
unset yrange

set xlabel "Time"
set ylabel "Concurrence"

set xrange [0.0:20]
set yrange [0.69:1.0]

set xlabel font ", 12"
set ylabel font ", 12"

plot for [i=1:60] '../../tests/q13_t20_b1.4/test_'.i.'.dat' using 1:3 title '' with l lc 'forest-green' lw 1, \
     for [j=1:60] '../../tests/q13_t20_b1.55/test_'.j.'.dat' using 1:3 title '' with l lc 'red' lw 1
