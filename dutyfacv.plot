## Limb (k) is starts swing at time (t) with parameters (F0), (kv), (Gu) and (Tsw);
## it has velocity (vv), stance time (ttst[k]), previous swing time (ttsw[k]) and duty factor (dutyfactor).

## simulation data saved as: 1 (k),  2 (t),  3 (F0),  4 (kv),  5 (Gu),  6 (Tsw), 7(vv), 8 (ttst[k]), 9 (ttsw[k]), 10 (dutyfactor)


# Axes label
unset key
set xlabel 'Velocity'
set ylabel 'Duty Factor(%)'
set key left top
set size ratio 0
# Axes ranges
set xrange [0:140]
set yrange [0:100]
set xtics 0,20,140
set ytics 0,10,100
set grid
# Functions
dfv(x)=(x>=20&&x<=110) ? 79.51-3.489*sqrt(x) : NaN
# Plot
plot dfv(x) title 'Experimental Data' with lines linestyle 1, "dur" u 7:10 title 'Model'