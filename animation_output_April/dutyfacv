# Axes label
unset key
set xlabel 'Velocity'
set ylabel 'Duty Factor(%)'
set key left top
# Axes ranges
set xrange [0:140]
set yrange [0:100]
set xtics 0,20,140
set ytics 0,10,100
set grid
# Functions
dfv(x)=(x>=20&&x<=110) ? 79.51-3.489*sqrt(x) : NaN
# Plot
plot dfv(x) title 'Experimental Data' with lines linestyle 1, "strideinfo" u 2:($5!=-1 && $3==1? 100*($5 * $6):NaN) title 'Model'