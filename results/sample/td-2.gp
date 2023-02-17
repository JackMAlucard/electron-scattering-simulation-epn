set encoding iso_8859_1
#set terminal pngcairo size 14 in, 14 in
#set output '1.d+1_1.d+5_1.d-3_1.d-0_trajectories_fit.png'
set title 'Time distribution, K_0 = 10 keV, {/Symbol q} = 6º'
set xlabel 'Trajectories simulated (time)'
set ylabel 'Fraction of electrons'
#set size square
set grid
#set xtics axis 1e5
#set ytics axis
#set format x '%.0tx10^{%T}'
#set format y '%.1tx10^{%T}'
set clip two
set key left top
f = 'td.dat'
t1 = 'Fraction of electrons embedded'
t2 = 'Fraction of electrons scattered'
c0 = '#e63431'
c1 = '#2c3b55'
c2 = '#445168'
c3 = '#5b667b'
c4 = '#8a92a1'
plot f using 3:4 lc rgb 'green' lt 1 lw 4 pt 7 ps 0.5 title t1 with lines, \
f using 3:5 lc rgb 'blue' lt 1 lw 4 pt 7 ps 0.5 title t2 with lines

pause mouse button2
exit