set terminal wxt size 900,600 enhanced font 'Verdana,20' persist
#set output 'Wendland_kernel.png'

set xlabel "step"
set ylabel "water level"
set title "Dynamic buffer - water level"
set key top left

file = 'waterLevel.dat'

p file u 1:2  lc "black" ps 2 title "max", file u 1:3 lc "green" title "-3h", file u 1:4 lc "red" title "-2h", file u 1:5 lc "blue" title "-h"
