# Axes label
set xlabel 'Stride Frequency'
set ylabel 'Swing/Stance Time'
set key left top
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
plot swsf(x) title 'Experimental Swing' with lines linestyle 1 linecolor 1, stsf(x) title 'Experimental Stance' with lines linestyle 2 linecolor 2, "strideinfo" u ($4!=-1 && $3==1? $6:NaN):4 title 'Model FL Swing' linecolor 1, "strideinfo" u ($5!=-1 && $3==1? $6:NaN):5 title 'Model FL Stance' linecolor 2,  "strideinfo" u ($4!=-1 && $3==3? $6:NaN):4 title 'Model HL Swing' linecolor 1, "strideinfo" u ($5!=-1 && $3==3? $6:NaN):5 title 'Model HL Stance' linecolor 2