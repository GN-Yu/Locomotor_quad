## output settings
set term gif size 1134,204

# Axes label
set xlabel ''
set ylabel ''
set title ''
set key noautotitle

# Axes ranges
set size ratio 0
set xrange [0:20]
set yrange [.9:4.1]
set xtics 0,0.1,20
set ytics ("FL" 1, "HR" 2, "FR" 3, "HL" 4)
set grid


# Plot Model 1 point A in heatmap
set out "single_sims_fig/StepCyc_model1_A.gif"
set xrange [18:18.5]
set xtics ("0" 18, "0.1" 18.1, "0.2" 18.2, "0.3" 18.3, "0.4" 18.4, "0.5" 18.5)
!g++ -O2 quad_loco.cc -o quad_loco.o && ./quad_loco.o -T 20 -fps 800 -v0 11.4 -v1 11.4 -Gu0 -10 -Gu1 -10 -kv0 1e10 -kv1 1e10 -tsw0 .02 -tsw1 .02 -kr 0.5 -presetTROT >dat 2>ppp
plot "dat" u 1:($11==0? 1:NaN) lw 3 pt 5, "dat" u 1:($14==0? 2:NaN) lw 3 pt 5, "dat" u 1:($12==0? 3:NaN) lw 3 pt 5, "dat" u 1:($13==0? 4:NaN) lw 3 pt 5

# Plot Model 1 point B in heatmap
set out "single_sims_fig/StepCyc_model1_B.gif"
set xrange [18:18.5]
set xtics ("0" 18, "0.1" 18.1, "0.2" 18.2, "0.3" 18.3, "0.4" 18.4, "0.5" 18.5)
!g++ -O2 quad_loco.cc -o quad_loco.o && ./quad_loco.o -T 20 -fps 800 -v0 13 -v1 13 -Gu0 -10 -Gu1 -10 -kv0 1e10 -kv1 1e10 -tsw0 .01 -tsw1 .05 -kr 0.5 -presetTROT >dat 2>ppp
plot "dat" u 1:($11==0? 1:NaN) lw 3 pt 5, "dat" u 1:($14==0? 2:NaN) lw 3 pt 5, "dat" u 1:($12==0? 3:NaN) lw 3 pt 5, "dat" u 1:($13==0? 4:NaN) lw 3 pt 5


# Plot Model 2 point A in heatmap
set out "single_sims_fig/StepCyc_model2_A.gif"
set xrange [2:2.5]
set xtics ("0" 2, "0.1" 2.1, "0.2" 2.2, "0.3" 2.3, "0.4" 2.4, "0.5" 2.5)
!g++ -O2 quad_loco.cc -o quad_loco.o && ./quad_loco.o -T 20 -fps 800 -v0 14 -v1 14 -Gu0 -10 -Gu1 -10 -kv0 0 -kv1 5 -tsw0 1000 -tsw1 1000 -kr 0.2 -presetWALK >dat 2>ppp
plot "dat" u 1:($11==0? 1:NaN) lw 3 pt 5, "dat" u 1:($14==0? 2:NaN) lw 3 pt 5, "dat" u 1:($12==0? 3:NaN) lw 3 pt 5, "dat" u 1:($13==0? 4:NaN) lw 3 pt 5

# Plot Model 2 point B in heatmap
set out "single_sims_fig/StepCyc_model2_B.gif"
set xrange [8.5:9]
set xtics ("0" 8.5, "0.1" 8.6, "0.2" 8.7, "0.3" 8.8, "0.4" 8.9, "0.5" 9)
!g++ -O2 quad_loco.cc -o quad_loco.o && ./quad_loco.o -T 20 -fps 800 -v0 18 -v1 18 -Gu0 -10 -Gu1 -10 -kv0 0 -kv1 5 -tsw0 1000 -tsw1 1000 -kr 0.2 -presetWALK >dat 2>ppp
plot "dat" u 1:($11==0? 1:NaN) lw 3 pt 5, "dat" u 1:($14==0? 2:NaN) lw 3 pt 5, "dat" u 1:($12==0? 3:NaN) lw 3 pt 5, "dat" u 1:($13==0? 4:NaN) lw 3 pt 5


# Plot Model 3 point A in heatmap
set out "single_sims_fig/StepCyc_model3_A.gif"
set xrange [1.5:2]
set xtics ("0" 1.5, "0.1" 1.6, "0.2" 1.7, "0.3" 1.8, "0.4" 1.9, "0.5" 2)
!g++ -O2 quad_loco.cc -o quad_loco.o && ./quad_loco.o -T 20 -fps 800 -v0 14 -v1 14 -Gu0 -10 -Gu1 -10 -kv0 0 -kv1 10 -tsw0 1000 -tsw1 1000 -kr 0.3 -inhib -presetWALK >dat 2>ppp
plot "dat" u 1:($11==0? 1:NaN) lw 3 pt 5, "dat" u 1:($14==0? 2:NaN) lw 3 pt 5, "dat" u 1:($12==0? 3:NaN) lw 3 pt 5, "dat" u 1:($13==0? 4:NaN) lw 3 pt 5

# Plot Model 3 point B in heatmap
set out "single_sims_fig/StepCyc_model3_B.gif"
set xrange [4.7:5.2]
set xtics ("0" 4.7, "0.1" 4.8, "0.2" 4.9, "0.3" 5, "0.4" 5.1, "0.5" 5.2)
!g++ -O2 quad_loco.cc -o quad_loco.o && ./quad_loco.o -T 20 -fps 800 -v0 16 -v1 16 -Gu0 -10 -Gu1 -10 -kv0 0 -kv1 10 -tsw0 1000 -tsw1 1000 -kr 0.3 -inhib -presetWALK >dat 2>ppp
plot "dat" u 1:($11==0? 1:NaN) lw 3 pt 5, "dat" u 1:($14==0? 2:NaN) lw 3 pt 5, "dat" u 1:($12==0? 3:NaN) lw 3 pt 5, "dat" u 1:($13==0? 4:NaN) lw 3 pt 5