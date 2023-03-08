## Limb (k) is starts swing at time (t) with parameters (F0), (kv), (Gu) and (Tsw);
## it has velocity (vv), stance time (ttst[k]), previous swing time (ttsw[k]) and duty factor (dutyfactor).

## simulation data saved as: 1 (k),  2 (t),  3 (F0),  4 (kv),  5 (Gu),  6 (Tsw), 7(vv), 8 (ttst[k]), 9 (ttsw[k]), 10 (dutyfactor)


# Axes label
set xlabel 'Velocity'
set ylabel 'Stride Length'
set key left top
set size ratio 0
# Axes ranges
set xrange [0:140]
set yrange [0:20]
set grid
# Functions
slv(x)=(x>=20&&x<=120) ? -3.64+3.048*log(x) : NaN
# Plot
plot slv(x) title 'Experimental Data' with lines linestyle 1, "dur" u 2:($4!=-1 && $3==1? $2/$6 : NaN) title 'Model'
