# Axes label
set xlabel 'Velocity'
set ylabel 'Stride Length'
set key left top
# Axes ranges
set xrange [0:140]
set yrange [0:20]
set grid
# Functions
slv(x)=(x>=20&&x<=120) ? -3.64+3.048*log(x) : NaN
# Plot
plot slv(x) title 'Experimental Data' with lines linestyle 1, "strideinfo" u 2:($4!=-1 && $3==1? $2/$6 : NaN) title 'Model'
