#Gnuplot script

set terminal postscript size 7in,7in eps color enhanced "Arial" 20 
set output "radiation_force.eps"


set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')

set size square
set multiplot layout 2,2

set xrange[-1.0 : 1.0]
set yrange[-1.0 : 1.0]

scale=0.8e-1;

set title "radiation force on x=0 plane"
set xlabel "{/Arial-Italic y}"
set ylabel "{/Arial-Italic z}"
plot "force_yz.txt" u 1:2:(scale*$4/sqrt($4*$4+$5*$5)):(scale*$5/sqrt($4*$4+$5*$5)):(sqrt($4*$4+$5*$5)) w vector lc palette ti ""

set title "radiation force on y=0 plane"
set xlabel "{/Arial-Italic x}"
set ylabel "{/Arial-Italic z}"
plot "force_xz.txt" u 1:2:(scale*$3/sqrt($3*$3+$5*$5)):(scale*$5/sqrt($3*$3+$5*$5)):(sqrt($3*$3+$5*$5)) w vector lc palette ti ""

set title "radiation force on z=0 plane"
set xlabel "{/Arial-Italic x}"
set ylabel "{/Arial-Italic y}"
plot "force_xy.txt" u 1:2:(scale*$3/sqrt($3*$3+$4*$4)):(scale*$4/sqrt($3*$3+$4*$4)):(sqrt($3*$3+$4*$4)) w vector lc palette ti ""

