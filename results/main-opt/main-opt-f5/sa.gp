set encoding iso_8859_1
set terminal png size 1280, 720 font "Helvetica,25"
set output 'sa.png'
set title 'Scattering Angles, K_0 = 9 keV, {/Symbol q} = 8º'
set xlabel '{/Symbol b} [º]'
set ylabel '{/Symbol a} [º]'
#set size square
set grid
set xrange [-25:25]
set yrange [0:30]
#set xtics axis 1e5
#set ytics axis
#set format x '%.0tx10^{%T}'
#set format y '%.1tx10^{%T}'
set clip two
#set key left top
set palette rgbformulae 3,11,6#Hot palette
f = 'sea.dat'
plot f using 2:1:3 notitle palette ps 2 pt 7 with points
#pause mouse button2
exit