# Axes label
set terminal qt size 1000,200
set xlabel 'time'
set ylabel 'limb'
set key noautotitle
set size ratio 0
# Axes ranges
set size ratio 0
set xrange [0:20]
set yrange [.9:4.1]
set xtics 0,0.1,20
set ytics ("FL" 1, "HR" 2, "FR" 3, "HL" 4)
set grid
# Plot
plot "dat" u 1:($11==0? 1:NaN), "dat" u 1:($14==0? 2:NaN), "dat" u 1:($12==0? 3:NaN), "dat" u 1:($13==0? 4:NaN)