## Limb (k) is starts swing at time (t) with parameters (F0), (kv), (Gu) and (Tsw);
## it has velocity (vv), stance time (ttst[k]), previous swing time (ttsw[k]) and duty factor (dutyfactor).

## simulation data saved as: 1 (k),  2 (t),  3 (F0),  4 (kv),  5 (Gu),  6 (Tsw), 7(vv), 8 (ttst[k]), 9 (ttsw[k]), 10 (dutyfactor)

# Axes label
set xlabel 'Stride Frequency'
set ylabel 'Swing/Stance Time'
set key left top
set size ratio 0
# Axes ranges
set xrange [0:12]
set yrange [0:0.4]
set xtics 0,1,12
set ytics 0,0.05,0.4
set grid
# Functions
swsf(x)=(x>=3 && x<=10) ? 0.121-0.006*x : NaN
stsf(x)=(x>=3 && x<=10) ? -0.03+0.72/x : NaN
# Plot
plot swsf(x) title 'Experimental Swing' with lines linestyle 1 linecolor 1, stsf(x) title 'Experimental Stance' with lines linestyle 2 linecolor 2, "dur" u (1/($8+$9)):($1==1? $9:NaN) title 'Model FL Swing' linecolor 1, "dur" u (1/($8+$9)):($1==1? $8:NaN) title 'Model FL Stance' linecolor 2,  "dur" u (1/($8+$9)):($1==3? $9:NaN) title 'Model HL Swing' linecolor 1, "dur" u (1/($8+$9)):($1==3? $8:NaN) title 'Model HL Stance' linecolor 2