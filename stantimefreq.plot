# Axes label
set xlabel 'stride frequency'
set ylabel 'stance time'
set key left top
# Axes ranges
set xrange [0:14]
set yrange [0:0.4]
# Functions
stsf(x)=(x>=2 && x<=11) ? -0.03+0.72/x : NaN
# Plot
plot stsf(x) title 'Experimental Data' with lines linestyle 1, "strideinfo" u ($5!=-1 && $3==1? $6:NaN):5 title 'Model'