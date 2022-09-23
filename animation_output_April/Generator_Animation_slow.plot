!g++ -O2 Mouse_trot.cc -o Mouse_trot && ./Mouse_trot -T 3 -dt .001 -v0 8 -v1 8 -Gu0 .14 -Gu1 .14 -kv0 1.3 -kv1 1.3 -slowi >dat 2>ppp
set term png size 600,600
set size ratio -1
unset key
set grid
unset label; unset arrow; unset object;
load "ppp"
!ffmpeg -r 100 -i tmp%d.png -y -hide_banner -loglevel error Walk15.mp4
!rm tmp*