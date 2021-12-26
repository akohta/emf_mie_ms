# emf_mie_ms
This is the electromagnetic field analysis program for spherical particles. This is based on Mie scattering theory. 
This program can analyze multiple scattering between spheres by using iterative solution. 
The electromagnetic field analysis program "multi_fbeam" is used to analyze incident field. 
GNU Scientific Library and libpng are required.


## Usage of example code

1. type 'make' command to compile.  
   The executable mie_ms_solver, example1.out, example2.out, example3.out are created. 
   The mie_ms_solver is the solver of Mie coefficients. 
   The example1.out is the executable of source code example1.c, it shows a simplest example using "emf_mie_ms". 
   The example2.out is the executable of source code example2.c, it shows a example of electromagnetic field intensity analysis. 
   The example3.out is the executable of source code example3.c, 
   it shows a example of outputting the instantaneous value of the electromagnetic field as an image.
   The example4.out is the executable of source code example4.c, it shows a example of far-field intensity analysis.  
   
2. type './mie_ms_solver' with an argument of output datafile name.  
   For example, './mie_ms_solver ex.dat'. 
   The beam datafile "ipw.txt" (plane wave is defined) and the sphere datafile "msphr.txt" (radius 2.5 water droplet is defined) are used. 
   This executable calcluates Mie coefficients, outputs them to binary file with the specified file name.
   This program searches for a sphere datafile in current directory using the default datafile name "msphr.txt". 
   As a simple representation of the analysis model, the nodes used for the surface integral are output as point cloud data. 
   In this case, the file "ex.particles" is output, and the visualization result is "ex_particle.png".
   The image was created using gnuplot script "gscript_particles.plt" and converted eps to png by using ImageMagick.  
   
3. type './example1.out' with an argument of datafile name.   
   For example, './example1.out ex.dat'. 
   This executable calculates electromagnetic field, radiation force and torque. This is the simplest example using this code.   
   
4. type './example2.out' with an argument of datafile name.  
   For example, './example2.out ex.dat'. 
   This executable calculates electromagnetic field intensity distributions, outputs them to text files.
   The I_example2.png is the visualization result of electromagnetic field intensity distributions, created by gnuplot script gscript_example2.plt
   (converted eps to png by using ImageMagick).  
   
5. type './example3.out' with an argument of datafile name.  
   For example, './example3.out ex.dat'.
   This executable calculates instantaneous value of the electromagnetic fields, outputs them to png image files.
   The image files are output to the folder which has a name adding "images" to the datafile name specified in the argument (file-extension is excluded). 
   Each image file has a name that indicates the cross section, field component, and number of time steps (ex. xz_Ex_014.png). 
   The color bar is output as color_bar.png in the same folder.
   The range of color bar in each cross section is output to the info.txt file (ex. xy_info.txt for z=0 plane).
   The xz_Ex.gif, yz_Ex.gif and xy_Ex.gif are animated gifs that concatenate the png files created by using the shell script gif_animation.sh.  
   
6. type './example4.out' with an argument of datafile name.  
   For example, './example4.out ex.dat'. 
   This executable calculates far-field intensity distributions and outputs them to text files. 
   The I_example4.png is the visualization result of electric field intensity distributions, created by gnuplot script gscript_example4.plt.  
   
Please see exmie_src/emf_mie_ms.h for detail of functions, mfb_src/multi_fbeam.h for detail of incident fields. 
The mie_ms_solver, example2.out and example3.out are parallelized using OpenMP. 
The number of threads is controlled by the environment variable OMP_NUM_THREADS.
The file named make_icx is the makefile for Intel compiler. 
For this code, Intel compiler is about 1.5 times faster than gcc (using gcc version 9.3.0, icx version 2021.4.0). 
The additional analysis example of single sphere is in the folder analysis_sample1.  

![point cloud data](ex_particles.png "nodes for surface integral (ex_particles.png)")  
![intensity distributions](I_example2.png "intensity distributions (I_example2.png)")  
![xz_Ex.gif](xz_Ex.gif "instantaneous value of the E_x on y=0 plane (xz_Ex.gif)")![yz_Ex.gif](yz_Ex.gif "instantaneous value of the E_x on x=0 plane (yz_Ex.gif)")  
![xy_Ex.gif](xy_Ex.gif "instantaneous value of the E_x on z=0 plane (xy_Ex.gif)")  
![far-field intensity](I_example4.png "far-field intensity distributions (I_example4.png)")  


## Analysis sample of multi-spheres (in the folder analysis_sample2)

The analysis sample of multi-spheres is in the folder analysis_sample2. 
The multiple scattering between spheres can be analyzed by the executable mie_ms_solver.
The usage is the same as example code. 
For exampe, '../mie_ms_solver ex2.dat', '../example2.out ex2.dat' (run in analysis_sample2 folder). 
The I_example2.png in this folder is the visualization result of electromagnetic field intensity distributions.  

![point cloud data 2](analysis_sample2/ex2_particles.png "nodes for surface integral (analysis_sample2/ex_particles.png)")  
![intensity distributions 2](analysis_sample2/I_example2.png "intensity distributions (analysis_sample2/I_example2.png)")  
![xz_Ex.gif 2](analysis_sample2/xz_Ex.gif "instantaneous value of the E_x on y=0 plane (analysis_sample2/xz_Ex.gif)")![yz_Ex.gif 2](analysis_sample2/yz_Ex.gif "instantaneous value of the E_x on x=0 plane (analysis_sample2/yz_Ex.gif)")  
![xy_Ex.gif 2](analysis_sample2/xy_Ex.gif "instantaneous value of the E_x on z=0 plane (analysis_sample2/xy_Ex.gif)")  
![far-field intensity 2](analysis_sample2/I_example4.png "far-field intensity distributions (analysis_sample2/I_example4.png)")  

## Analysis sample of radiation force (in the folder analysis_sample3)  

The code in the folder analysis_sample3 is the radiation force analysis program for optical trapping. 
The usage is 'make' and './radiation_force.out'.
The radiation_force.png is the visualization result of radiation forces, created by Gnuplot script gscript_radiation_force.plt.  

![radiation force analysis](analysis_sample3/radiation_force.png "vector plot of radiation force (analysis_sample3/radiation_force.png)")


## System of units

This program use the own defined system of units (OSU), optimized for optics. 
The system of units is defined as <img src="https://latex.codecogs.com/gif.latex?c_0=1"> ( speed of light in vacuum ), 
<img src="https://latex.codecogs.com/gif.latex?\mu_0=1"> ( permeability of vacuum ). 
For the conversion from OSU to MKSA system of units, the unit of length in OSU is defined as 
<img src="https://latex.codecogs.com/gif.latex?1\times10^{-6}"> [m] in MKSA, the unit of power in OSU is defined as
<img src="https://latex.codecogs.com/gif.latex?1\times10^{-3}"> [W] in MKSA. The conversions of base unit are follows.  
<img src="https://latex.codecogs.com/gif.latex?a=1\times10^{-6}">,  
<img src="https://latex.codecogs.com/gif.latex?b=1\times10^{-3}">,  
<img src="https://latex.codecogs.com/gif.latex?a\,\mathrm{[m]}=1\,\mathrm{[L]}">,  
<img src="https://latex.codecogs.com/gif.latex?\frac{ab}{c_0^3}\,\mathrm{[kg]}=1\,\mathrm{[M]}">,  
<img src="https://latex.codecogs.com/gif.latex?\frac{a}{c_0}\,\mathrm{[s]}=1\,\mathrm{[T]}">,  
<img src="https://latex.codecogs.com/gif.latex?\sqrt{\frac{b}{c_0\mu_0}}\,\mathrm{[A]}=1\,\mathrm{[I]}">.  
Please see com_src/osu_mksa.h and com_src/osu_mksa.c for detail of conversions.


## References
1. Barton, J. P., D. R. Alexander, and S. A. Schaub. "Internal and near‐surface electromagnetic fields for a spherical particle irradiated by a focused laser beam." Journal of Applied Physics 64.4 (1988): 1632-1639.  
2. Barton, J. P., D. R. Alexander, and S. A. Schaub. "Theoretical determination of net radiation force and torque for a spherical particle illuminated by a focused laser beam." Journal of Applied Physics 66.10 (1989): 4594-4602.  
3. Barton, John P., et al. "Electromagnetic field for a beam incident on two adjacent spherical particles." Applied optics 30.33 (1991): 4706-4715.  
4. Abramowitz, Milton, and Irene A. Stegun, eds. Handbook of mathematical functions with formulas, graphs, and mathematical tables. Vol. 55. US Government printing office, 1948.  
5. GNU Scientific Library [GSL](https://www.gnu.org/software/gsl/)
6. The command-line driven graphing utility [gnuplot](http://www.gnuplot.info/)  
7. The utilities for manipulating images [ImageMagick](https://imagemagick.org/)  
8. The official PNG reference library [libpng](http://www.libpng.org/pub/png/libpng.html)  
9. The electromagnetic field analysis program [multi_fbeam](https://github.com/akohta/multi_fbeam/)  


The formula (12) in the Reference 2 ( z-component of radiation torque ) is misprinted. The following formula is correct.  
<img src="https://latex.codecogs.com/gif.latex?\frac{\left<N_z\right>}{a^3E_0^2}=-\frac{a}{8\pi}\sum_{l=1}^{\infty}\sum_{m=-l}^{l}l(l+1)m\left[\epsilon_{\mathrm{ext}}|a_{lm}|^2+|b_{lm}|^2+\Re(\epsilon_{\mathrm{ext}}a_{lm}A_{lm}^*+b_{lm}B_{lm}^*)\right]">.  
The first letter <img src="https://latex.codecogs.com/gif.latex?l"> in the term of sum is missed.
