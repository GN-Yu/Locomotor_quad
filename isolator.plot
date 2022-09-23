# Axes label
set xlabel 'Time'
set ylabel 'phi'
set key left top
# Axes ranges
set xrange [0:10]
set yrange [-3.2:3.2]
set xtics 0,1,10
set ytics -3.2,0.5,3.2
set grid
set size ratio 0
# Functions
stdOmega(x)=62.8
# Plot
plot "dat" u 1:15 title 'FL' with lines, "dat" u 1:16 title 'FR' with lines, "dat" u 1:17 title 'HL' with lines, "dat" u 1:18 title 'HR' with lines