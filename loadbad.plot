# Axes label
set xlabel 'Time'
set ylabel 'Total Load'
set key right top
# Axes ranges
set xrange [19.5:20]
set yrange [0:1]
set xtics 19.5,0.1,20
set ytics 0,0.2,1
set grid
# Plot
plot "dat" u ($1>=19.5 && $1<=20 ? $1:NaN):6 w l title '125cm/s'