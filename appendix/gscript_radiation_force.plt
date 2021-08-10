#Gnuplot script

set terminal postscript size 7in,7in eps color enhanced "Arial" 20 
set output "radiation_force.eps"


set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')

set size square
set multiplot layout 2,2

set xrange[-1.0e-6 : 1.0e-6]
set yrange[-1.0e-6 : 1.0e-6]

scale=0.8e-7;

set title "y-z plane (x=0) radiation force"
set xlabel "{/Arial-Italic y} (m)"
set ylabel "{/Arial-Italic z} (m)"
plot "force_yz.txt" u 1:2:(scale*$4/sqrt($4*$4+$5*$5)):(scale*$5/sqrt($4*$4+$5*$5)):(sqrt($4*$4+$5*$5)) w vector lc palette ti ""

set title "x-z plane (y=0) radiation force"
set xlabel "{/Arial-Italic x} (m)"
set ylabel "{/Arial-Italic z} (m)"
plot "force_xz.txt" u 1:2:(scale*$3/sqrt($3*$3+$5*$5)):(scale*$5/sqrt($3*$3+$5*$5)):(sqrt($3*$3+$5*$5)) w vector lc palette ti ""

set title "x-y plane (z=0) radiation force"
set xlabel "{/Arial-Italic x} (m)"
set ylabel "{/Arial-Italic y} (m)"
plot "force_xy.txt" u 1:2:(scale*$3/sqrt($3*$3+$4*$4)):(scale*$4/sqrt($3*$3+$4*$4)):(sqrt($3*$3+$4*$4)) w vector lc palette ti ""

