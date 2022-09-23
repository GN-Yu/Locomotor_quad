!g++ -O2 Mouse_pace.cc -o Mouse_pace && ./Mouse_pace -T 15 -dt .001 -v0 5 -v1 30 -Gu0 .01 -Gu1 .25 -kv0 .5 -kv1 12 -slowi >dat 2>ppp
set term png size 600,600
set size ratio -1
unset key
set grid
unset label; unset arrow; unset object;
load "ppp"
!ffmpeg -r 500 -i tmp%d.png -y -hide_banner -loglevel error pace.mp4
!rm tmp*