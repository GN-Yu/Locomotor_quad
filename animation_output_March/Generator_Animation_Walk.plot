!g++ -O2 walk2.cc -o walk2 && ./walk2 -T 2 -dt .001 -v0 20 -v1 20 -Gu0 .1 -Gu1 .1 -Gv0 .98 -Gv1 .98 -newi >dat 2>ppp
!g++ -O2 walk2.cc -o walk2 && ./walk2 -T 15 -dt .001 -v0 20 -v1 5 -Gu0 .1 -Gu1 .01 -Gv0 .98 -Gv1 .99999 -i -newi >dat 2>ppp
!g++ -O2 walk2.cc -o walk2 && ./walk2 -T 2 -dt .001 -v0 5 -v1 5 -Gu0 .01 -Gu1 .01 -Gv0 .99999 -Gv1 .99999 -i >dat 2>ppp
set term png size 600,600
set size ratio -1
unset key
set grid
unset label; unset arrow; unset object;
load "ppp"
# !ffmpeg -r 500 -i tmp%d.png -y -hide_banner -loglevel error Walk15.mp4
# !rm tmp*