# Axes label
set xlabel 'time'
set ylabel 'total load'
#set key noautotitle

# Axes ranges
set xrange [0:15]
set yrange [-1:1]
set xtics 0,0.5,20
set ytics -1,0.1,1
set grid
# Plot
plot "dat" u 1:6 w lp