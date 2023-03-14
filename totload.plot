# Axes label
set size ratio 0
set xlabel 'time'
set ylabel 'total load'
#set key noautotitle

# Axes ranges
set xrange [0:10]
set yrange [-1:1]
set xtics 0,0.1,20
set ytics -1,0.1,1
set grid
# Plot
plot "dat" u 1:6 w lp