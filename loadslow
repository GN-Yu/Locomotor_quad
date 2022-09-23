# Axes label
set xlabel 'Time'
set ylabel 'Total Load'
set key right top
# Axes ranges
set xrange [1:1.5]
set yrange [0:1]
set xtics 1,0.1,1.5
set ytics 0,0.2,1
set grid
# Plot
plot "dat" u ($1>=1 && $1<=1.5 ? $1:NaN):6 w l title '20cm/s'
