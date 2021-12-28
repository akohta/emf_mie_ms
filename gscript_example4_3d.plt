# Gnuplot script file
set terminal postscript eps color enhanced "Arial" 24 size 7in,10in
set output "I_example4_3d.eps"

file="fsIe_3d.txt"
angle1=45
angle2=45

set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')

set multiplot layout 2,1

set pm3d depthorder
set angles degree
set view equal xyz
set ticslevel 0

set xrange [-1.2:1.2]
set xlabel "{/Arial-Italic x} [a.u.]"
set yrange [-1.2:1.2]
set ylabel "{/Arial-Italic y} [a.u.]"
set zrange [-1.2:1.2] 
set zlabel "{/Arial-Italic z} [a.u.]"

set title "scattered field intensity distribution in far-field"
set view angle1,angle2
set logscale cb
splot file using (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):3 with pm3d title ""

set title ""
set view angle1+180,angle2+90
set logscale cb
splot file using (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):3 with pm3d title ""

