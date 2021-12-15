# Gnuplot script file for visualization of "radiatoin_force_a.out" output 
# usage on command line : gnuplot -e 'file="output_filename"' gscript_force_a.plt

set terminal postscript eps color enhanced "Arial" 22 size 7in,4in
set output "force_a_out.eps"
set encoding iso

if ( !exists("file") ) file="force_a.txt"

set title "radiation force and torque acting on the sphere"

set xlabel "sphere radius [\265m]"
set ylabel  "z component of radiation force [N]"
set y2label "radiation torque around the y-axis [Nm]"
set ytics nomirror
set y2tics nomirror

plot file using 1:4 axis x1y1 with line title "force", file using 1:6 axis x1y2 with line title "torque"
