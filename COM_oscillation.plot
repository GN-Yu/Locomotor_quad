# ocilliation of the COM in a stride

# Limb (k) is starts swing at time (t) with parameters (F0), (kv), (Gu) and (Tsw);
# it has velocity (vv), stance time (ttst[k]), previous swing time (ttsw[k]) and duty factor (dutyfactor).

# simulation data saved as: 1 (k),  2 (t),  3 (F0),  4 (kv),  5 (Gu),  6 (Tsw), 7(vv), 8 (ttst[k]), 9 (ttsw[k]), 10 (dutyfactor) 11 (COM x coordinate) 12 (COM y coordinate) 13 (body angle) 14 (number of limbs in swing)

# Axes label
set xlabel 'Velocity'
set ylabel 'Stride Length'
set key left top
set size ratio 0
# Axes ranges
set xrange [0:10]
set yrange [0:5]
set grid


# Plot
plot "dur" u 2:($4!=-1 && $3==1? $2/$6 : NaN) title 'Model' with lines linestyle 1