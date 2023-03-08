# Axes label
set xlabel 'Time'
set ylabel 'Total Load'
set key right top
# Axes ranges
set xrange [17:17.5]
set yrange [0:1]
set xtics 17,0.1,17.5
set ytics 0,0.2,1
set grid
# Plot
plot "dat" u ($1>=17 && $1<=17.5 ? $1:NaN):6 w l title '110cm/s'