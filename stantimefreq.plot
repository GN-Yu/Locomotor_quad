## Limb (k) is starts swing at time (t) with parameters (F0), (kv), (Gu) and (Tsw);
## it has velocity (vv), stance time (ttst[k]), previous swing time (ttsw[k]) and duty factor (dutyfactor).

## simulation data saved as: 1 (k),  2 (t),  3 (F0),  4 (kv),  5 (Gu),  6 (Tsw), 7(vv), 8 (ttst[k]), 9 (ttsw[k]), 10 (dutyfactor)

# Axes label
set xlabel 'stride frequency'
set ylabel 'stance time'
set key left top
set size ratio 0
# Axes ranges
set xrange [0:14]
set yrange [0:0.4]
# Functions
stsf(x)=(x>=2 && x<=11) ? -0.03+0.72/x : NaN
# Plot
plot stsf(x) title 'Experimental Data' with lines linestyle 1, "dur" u (1/($8+$9)):8 title 'Model'