set termopt enhanced

unset label
unset xrange
unset yrange
unset title

set size 500,300

set multiplot layout 1,2 title "b_{parallel}=1.4, J_{e}=1, J_{ce}=0.03, J_{cs}=2, alpha=0.5\n" font ", 18"
#set multiplot layout 1,2 title "B=0.5, W=0, J_{cs}=1, J_{z}=1, alpha=0.5\n" font ", 18"

set xlabel "Site"
set ylabel "Time"

set autoscale xfix
set autoscale yfix
set autoscale cbfix

set xlabel font ", 16"
set ylabel font ", 16"

set xtics font ", 12"
set ytics font ", 12"


plot '../../tests/epr1/test_epr1_fig1_bpa1.40_Je1.00_Jce0.03_Jcs2.00_a0.50_q9_t50__26353241.dat' every :::0::10 matrix with image notitle


set xlabel "Site"
set ylabel ""

plot '../../tests/epr1/test_epr1_fig1_bpa1.40_Je1.00_Jce0.03_Jcs2.00_a0.50_q9_t50__26353241.dat' matrix with image notitle




unset multiplot
