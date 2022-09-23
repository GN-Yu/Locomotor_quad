set term gif
set out "Duty_factor_temp.gif"
set key 
set cbrange [50:100]
set xlabel "Swing Duration"
set ylabel "COM Velocity"
set xrange [0:0.12]
set yrange [0:80]
set xtics 0,0.02,0.12
set ytics 0,20,80
set title "Duty Factor (for FL)"

velocity(x)=(x>0.05 && x<=0.1) ? ((0.121-x)/0.006 - 1.66338)/0.08251 : NaN
dfv(x)=79.51-3.489*sqrt(velocity(x))

#set parametric
#Tsw(t)=(t>0 && t<=100) ? 0.121-0.006*(1.004+0.3978*t/log(t)) : NaN
#velocity(t)=(t>0.05 && t<=0.1) ? ((0.121-t)/0.006 - 1.66338)/0.08251 : NaN
#dfv(t)=79.51-3.489*sqrt(velocity(t))

#0.114976-0.0023868*t/log(t) : NaN

#data saved as: 1 time Tswc,  2 force F0,  3 actual velocity vv,  4 leg indicator k,  6 stance duration ttst[k],  7 duty factor (100*ttst[k]/stridetime)
plot velocity(x) title 'Experimental Data' w l lw 2 lc palette, "<cat data/*" u 1:($4==1? $3:NaN):7 w p pt 5 lc palette z notitle
#plot x:velocity(x):dfv(x) title 'Experimental Data' w l lw 7 lc palette z