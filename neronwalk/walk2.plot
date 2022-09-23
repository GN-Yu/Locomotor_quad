!cp ini20 ini
!g++ -O2 walk2.cc -o walk2 && ./walk2 -v 20 20 -T 1 -dt .001 -w .99 -el -56 -i >dat 2>ppp
set term x11; set size ratio -1; set xrange [*:*]; set yrange [*:*]; load 'ppp'
pause -1
