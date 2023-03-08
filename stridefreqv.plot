## Limb (k) is starts swing at time (t) with parameters (F0), (kv), (Gu) and (Tsw);
## it has velocity (vv), stance time (ttst[k]), previous swing time (ttsw[k]) and duty factor (dutyfactor).

## simulation data saved as: 1 (k),  2 (t),  3 (F0),  4 (kv),  5 (Gu),  6 (Tsw), 7(vv), 8 (ttst[k]), 9 (ttsw[k]), 10 (dutyfactor)


# Axes label
set xlabel 'Velocity'
set ylabel 'Stride Frequency'
set key left top
set size ratio 0
# Axes ranges
set xrange [0:140]
set yrange [0:16]
set xtics 0,10,140
set ytics 0,2,16
set grid
# Functions
sfv(x)=(x>=20&&x<=120) ? 1.004+0.3978*x/log(x) : NaN
# Plot
plot sfv(x) title 'Experimental Data' with lines linestyle 1, "dur" u 7:(1/($8+$9)) title 'Model'