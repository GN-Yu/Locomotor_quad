## Limb (k) is starts swing at time (t) with parameters (F0), (kv), (Gu) and (Tsw);
## it has velocity (vv), stance time (ttst[k]), previous swing time (ttsw[k]) and duty factor (dutyfactor).

## simulation data saved as: 1 (k),  2 (t),  3 (F0),  4 (kv),  5 (Gu),  6 (Tsw), 7(vv), 8 (ttst[k]), 9 (ttsw[k]), 10 (dutyfactor)

# Axes label
set xlabel 'time'
set ylabel 'limb number'
set key noautotitle
set size ratio 0
# Axes ranges
set size ratio 0
set xrange [0:20]
set yrange [0:5]
set xtics 0,0.1,20
set ytics ("FL" 1, "FR" 2, "HL" 3, "HR" 4)
set grid
# Plot
plot "dat" u 1:($11==0? 1:NaN), "dat" u 1:($12==0? 2:NaN), "dat" u 1:($13==0? 3:NaN), "dat" u 1:($14==0? 4:NaN)