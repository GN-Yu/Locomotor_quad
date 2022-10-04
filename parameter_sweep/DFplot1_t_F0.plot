set term gif
set out "DFmodel1_t_F0.gif"
set key 
set cbrange [50:100]
set xlabel "Swing Duration"
set ylabel "Nomalized Force"
set xrange [0:0.1]
set yrange [0:50]
set xtics 0,0.02,0.1
set ytics 0,5,50 
set title "Duty Factor (for FL)"

#data saved as: 1 time Tswc,  2 force F0,  3 actual velocity vv,  4 leg indicator k,  6 stance duration ttst[k],  7 duty factor (100*ttst[k]/stridetime)
plot "<cat data/*" u 1:($4==1? $2:NaN):7 w p pt 5 lc palette z notitle