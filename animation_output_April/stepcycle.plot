# Axes label
set xlabel 'time'
set ylabel 'limb number'
set key noautotitle
# Axes ranges
set xrange [0:20]
set yrange [0:5]
set xtics 0,0.1,20
set ytics ("FL" 1, "FR" 2, "HL" 3, "HR" 4)
set grid
# Plot
plot "dat" u 1:($11==0? 1:NaN), "dat" u 1:($12==0? 2:NaN), "dat" u 1:($13==0? 3:NaN), "dat" u 1:($14==0? 4:NaN)