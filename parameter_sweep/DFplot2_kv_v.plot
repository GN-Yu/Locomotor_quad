set term gif
set out "DFmodel2_kv_v.gif"
set key 
set cbrange [50:100]
set xlabel "Balance Parameter"
set ylabel "COM Velocity"
set xrange [0:5]
set yrange [0:80]
set xtics 0,.5,5
set ytics 0,20,80
set title "Duty Factor (for FL)"

#data saved as: 1 time Tswc,  2 force F0,  3 actual velocity vv,  4 leg indicator k,  6 stance duration ttst[k],  7 duty factor (100*ttst[k]/stridetime)
plot "<cat data/*" u 1:($4==1? $3:NaN):7 w p pt 5 lc palette z notitle