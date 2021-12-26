# Gnuplot script file
# Usage : gnuplot -e 'file="particles_filename"' gscript_particles.plt

if ( !exists("file") ) file="ex.particles"

set terminal postscript eps color enhanced "Arial" 20 size 7in,7in
set output "particles.eps"

#set palette model RGB rgbformulae 35,13,10
set palette model RGB rgbformulae 22,13,-31
set view equal xyz
set ticslevel 0
set cbtics 1
set cblabel "Sphere ID"
set grid x y z vertical

set xlabel "{/Arial-Italic x}"
set ylabel "{/Arial-Italic y}" 
set zlabel "{/Arial-Italic z}"

splot file using 1:2:3:4 w d palette title "nodes for surface integral"
