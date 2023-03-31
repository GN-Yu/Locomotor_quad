reset session
set title ""

$DATA << EOD
0  1
90 1
320 1
230 1
EOD

set polar
set theta clockwise top
set angles degrees
set grid polar 45. lt -1 dt 0 lw 1
unset xlabel
unset ylabel
unset border
set ttics 0,90,360
unset xtics
unset ytics
unset raxis
set rtics 0,1,1
set format r ''
set rrange [0:1.04]
set size ratio -1

plot 1.03 lt -1 title "", $DATA lw 1.5 title ""