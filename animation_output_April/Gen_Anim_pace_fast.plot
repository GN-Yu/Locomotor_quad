!g++ -O2 Mouse_trot.cc -o Mouse_trot && ./Mouse_trot -T 2 -dt .001 -v0 30 -v1 30 -Gu0 .25 -Gu1 .25 -kv0 12 -kv1 12 -fasti >dat 2>ppp
set term png size 600,600
set size ratio -1
unset key; unset xtics; unset ytics; unset label; unset arrow; unset object;
load "ppp"
!ffmpeg -r 50 -i tmp%d.png -y -hide_banner -loglevel error Trot.mp4
!rm tmp*