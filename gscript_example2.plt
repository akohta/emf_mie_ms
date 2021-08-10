# Gnuplot script file
set terminal postscript eps color enhanced "Arial" 16 size 7in,7in
set output "I_example2.eps"

set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')

set pm3d map

set size square
set multiplot layout 2,2

set xrange [-4e-6 : 4e-6]
set yrange [-4e-6 : 4e-6]
set xtics -3e-6, 2.0e-6, 3e-6
set ytics -3e-6, 2.0e-6, 3e-6

set title "y=0 plane"
set xlabel "{/Arial-Italic x} (m)"
set ylabel "{/Arial-Italic z} (m)"
splot "Ie_xz.txt"

set title "x=0 plane"
set xlabel "{/Arial-Italic y} (m)"
splot "Ie_yz.txt"

set title "z=0 plane"
set xlabel "{/Arial-Italic x} (m)"
set ylabel "{/Arial-Italic y} (m)"
splot "Ie_xy.txt"

unset multiplot
set terminal x11
