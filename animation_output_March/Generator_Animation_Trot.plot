!g++ -O2 walk2.cc -o walk2 && ./walk2 -T 20 -dt .001 -v0 20 -v1 50 -Gu0 .1 -Gu1 .2 -Gv0 .98 -Gv1 .95 -newi >dat 2>ppp
!g++ -O2 walk2.cc -o walk2 && ./walk2 -T 0.6 -dt .001 -v0 50 -v1 50 -Gu0 .2 -Gu1 .2 -Gv0 .95 -Gv1 .95 -i >dat 2>ppp
set term png size 600,600
set size ratio -1
unset key
set grid
unset label; unset arrow; unset object;
load "ppp"
# !ffmpeg -r 100 -i tmp%d.png -y -hide_banner -loglevel error Trot110.mp4
# !rm tmp*