set encoding iso_8859_1
set terminal png size 1280, 720 font "Helvetica,25"
set output 'td.png'
set title "Electrons' final state, K_0 = 10 keV, {/Symbol q} = 6º"
set xlabel 'Number of simulated trajectories'
set ylabel 'Number of electrons'
#set size square
set grid
#set xtics axis 1e5
#set ytics axis
#set format x '%.0tx10^{%T}'
#set format y '%.1tx10^{%T}'
set clip two
set key left top
f = 'td.dat'
t1 = 'Embedded'
t2 = 'Scattered'
plot f using 3:1 lc rgb 'green' lt 1 lw 6 pt 7 ps 0.5 title t1 with lines, \
f using 3:2 lc rgb 'blue' lt 1 lw 6 pt 7 ps 0.5 title t2 with lines

#pause mouse button2
exit