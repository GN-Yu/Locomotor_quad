# Axes label
set xlabel 'Velocity'
set ylabel 'Stride Frequency'
set key left top
# Axes ranges
set xrange [0:140]
set yrange [0:16]
set xtics 0,10,140
set ytics 0,2,16
set grid
# Functions
sfv(x)=(x>=20&&x<=120) ? 1.004+0.3978*x/log(x) : NaN
# Plot
plot sfv(x) title 'Experimental Data' with lines linestyle 1, "strideinfo" u 2:($4!=-1 && $3==1? $6:NaN) title 'Model'