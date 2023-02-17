set encoding iso_8859_1
set terminal png size 1280, 720 font "Helvetica,25"
set output 'cp.png'
#set style fill transparent solid 0.25 noborder
set title 'Charge patch, K_0 = 8 keV, {/Symbol q} = 8º'
set xlabel 'z [a_0]'
set ylabel 'x [a_0]'
set size ratio -1
#set size square
set grid front
#set colorbox# horiz user origin .1,.02 size .8,.04

lz = 1.77*75*1.88973#In x-axis
lx = 1.77*30*1.88973#In y-axis
lzr = lz*1.1
#lxr = lx*1.5
#set style fill transparent solid 0.25 noborder
set style function filledcurves above y1=-lx
set style function filledcurves below y2=lx
set xrange [-lzr:lzr]#z range in data
set yrange [-lx*1.3:lx*1.5]#x range in data
#set xtics axis 1e5
#set ytics axis
#set format x '%.0tx10^{%T}'
#set format y '%.1tx10^{%T}'
set clip two
#set key left top
set view map
set palette defined (-1 "red", 0 "blue")
#set palette rgbformulae 3,2,2#3,11,6
#set palette rgb 

f = 'emb.dat'

f1(x) = lx
f2(x) = -lx
f3(x) = (x > -lz && x < lz) ? lx : 1/0
f4(x) = -f3(x)

plot f3(x) lc rgb "gray" title 'Material region',\
	f4(x) lc rgb "gray" notitle, \
	f using 3:1:2 palette ps 2 pt 7 with points title 'Embedded electrons'

#	f2(x) lc rgb "gray" notitle, \

#pause mouse button2
exit