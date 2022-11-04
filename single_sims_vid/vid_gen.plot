#!g++ -O2 quad_loco.cc -o quad_loco.o && ./quad_loco.o -T 20 -fps 500 -v0 14 -v1 14 -Gu0 -10 -Gu1 -10 -kv0 0 -kv1 5 -tsw0 1000 -tsw1 1000 -kr 0.2 -presetWALK >dat 2>ppp

load "ppp"

!ffmpeg -r 100 -i pics/tmp%d.png -y -hide_banner -loglevel error model1_A.mp4
# !rm pics/tmp*





# !g++ -O2 quad_loco.cc -o quad_loco.o && ./quad_loco.o -T 20 -fps 500 -v0 13 -v1 13 -Gu0 -10 -Gu1 -10 -kv0 1e10 -kv1 1e10 -tsw0 .01 -tsw1 .05 -kr 0.5 -presetTROT >dat 2>ppp



# !g++ -O2 quad_loco.cc -o quad_loco.o && ./quad_loco.o -T 20 -fps 500 -v0 14 -v1 14 -Gu0 -10 -Gu1 -10 -kv0 0 -kv1 5 -tsw0 1000 -tsw1 1000 -kr 0.2 -presetWALK >dat 2>ppp


# !g++ -O2 quad_loco.cc -o quad_loco.o && ./quad_loco.o -T 20 -fps 500 -v0 18 -v1 18 -Gu0 -10 -Gu1 -10 -kv0 0 -kv1 5 -tsw0 1000 -tsw1 1000 -kr 0.2 -presetWALK >dat 2>ppp



# !g++ -O2 quad_loco.cc -o quad_loco.o && ./quad_loco.o -T 20 -fps 500 -v0 14 -v1 14 -Gu0 -10 -Gu1 -10 -kv0 0 -kv1 10 -tsw0 1000 -tsw1 1000 -kr 0.3 -inhib -presetWALK >dat 2>ppp


# !g++ -O2 quad_loco.cc -o quad_loco.o && ./quad_loco.o -T 20 -fps 500 -v0 16 -v1 16 -Gu0 -10 -Gu1 -10 -kv0 0 -kv1 10 -tsw0 1000 -tsw1 1000 -kr 0.3 -inhib -presetWALK >dat 2>ppp
