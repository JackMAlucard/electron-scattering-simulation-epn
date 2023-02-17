set encoding iso_8859_1
#set terminal pngcairo size 14 in, 14 in
#set output 'cb.png'
set title 'Simple Silica (SiO_2) Model'
set xlabel 'z'
set ylabel 'x'
set zlabel 'y'
set view equal xyz
#set size square

#set xrange [-0.5:0.5]
#set yrange [-0.5:0.5]

f = 'ssm_Si.dat'
g = 'ssm_O.dat'
plot g using 3:1 notitle with points pointsize 2 pointtype 7 lc rgb 'blue', \
	f using 3:1 notitle with points pointsize 2 pointtype 7 lc rgb 'gray'

pause mouse button2

exit


#set key off
#splot 'haz.dat' using 3:2:1 with points palette pointsize 3 pointtype 7

#pause mouse button2

#exit
