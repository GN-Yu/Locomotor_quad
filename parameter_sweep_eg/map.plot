set term gif
set out "map.gif"
set cbrange [0:.5]
set xlabel "time"
set ylabel "F"
set title "stance duration"
plot "<cat data/*" u 1:($4==1?$2:NaN):5 w p pt 5 lc palette z