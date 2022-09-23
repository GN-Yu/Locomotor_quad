# Axes label
set xlabel 'Time'
set ylabel 'Omega'
set key left top
# Axes ranges
set xrange [0:10]
set yrange [0:200]
set xtics 0,1,10
set ytics 0,30,300
set grid
set size ratio 0
# Functions
stdOmega(x)=62.8
# Plot
plot stdOmega(x) title 'stdOmega' with lines linestyle 1 linecolor 1, "dat" u 1:19 title 'FL' with lines, "dat" u 1:20 title 'FR' with lines, "dat" u 1:21 title 'HL' with lines, "dat" u 1:22 title 'HR' with lines