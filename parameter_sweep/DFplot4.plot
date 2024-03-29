## Limb (k) is starts swing at time (t) with parameters (F0), (kv), (Gu) and (Tsw);
## it has velocity (vv), stance time (ttst[k]), previous swing time (ttsw[k]) and duty factor (dutyfactor).

## simulation data saved as: 1 (k),  2 (t),  3 (F0),  4 (kv),  5 (Gu),  6 (Tsw), 7(vv), 8 (ttst[k]), 9 (ttsw[k]), 10 (dutyfactor)

## save actual experimental data in table
set table $realdata
set param
set trange [20:80]
plot .121-.006*(1.004+.3978*t/log(t)),t
unset table
unset param


## output settings
set term gif
set key
set cbrange [40:100]


## model 4 plot 1
set out "DFmodel4_kv_F0.gif"
set xlabel "Balance Parameter"
set ylabel "Limb Horizontal Force"
set xrange [0:10]
set yrange [0:40]
set xtics 0,1,10
set ytics 0,5,40 
set title "Duty Factor Heat Map"

plot "<cat data/*" u 4:3:(($10<100 && $10>45)? $10:NaN) w p pt 5 lc palette z notitle

#plot "<cat data/*" u 4:($1==1? $3:NaN):10 w p pt 5 lc palette z notitle


## model 4 plot 2
set out "DFmodel4_kv_v.gif"
set xlabel "Balance Parameter"
set ylabel "COM Velocity"
set xrange [0:10]
set yrange [0:80]
set xtics 0,1,10
set ytics 0,10,80
set title "Duty Factor Heat Map"

plot "<cat data/*" u 4:7:(($10<100 && $10>45)? $10:NaN) w p pt 5 lc palette z notitle

#plot "<cat data/*" u 4:($1==1? $7:NaN):10 w p pt 5 lc palette z notitle


## model 4 plot 3
set out "DFmodel4_tsw_F0.gif"
set xlabel "Swing Duration"
set ylabel "Limb Horizontal Force"
set xrange [0:0.2]
set yrange [0:40]
set xtics 0,0.02,0.2
set ytics 0,5,40 
set title "Duty Factor Heat Map"

plot "<cat data/*" u 9:3:(($10<100 && $10>45)? $10:NaN) w p pt 5 lc palette z notitle

#plot "<cat data/*" u 9:($1==1? $3:NaN):10 w p pt 5 lc palette z notitle


## model 4 plot 4
set out "DFmodel4_tsw_v.gif"
set xlabel "Swing Duration"
set ylabel "COM Velocity"
set xrange [0:0.2]
set yrange [0:80]
set xtics 0,0.02,0.2
set ytics 0,10,80
set title "Duty Factor Heat Map"

plot "<cat data/*" u 9:7:(($10<100 && $10>45)? $10:NaN) w p pt 5 lc palette z notitle, $realdata u 1:2:(79.51-3.489*sqrt($2)) w l lw 5 lc -1 title "Experimental data"

#plot "<cat data/*" u 9:($1==1? $7:NaN):10 w p pt 5 lc palette z notitle, $realdata u 1:2:(79.51-3.489*sqrt($2)) w l lw 5 lc palette z notitle