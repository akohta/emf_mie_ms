# emf_mie_ms
This is the electromagnetic field analysis program for spherical particles. This is based on Mie scattering theory. 
This program can analyze multiple scattering between spheres by using iterative solution. 
The electromagnetic field analysis program "multi_fbeam" is used to analyze incident field. 
GNU Scientific Library is required.

## Usage of example code

1. type 'make' command to compile.  
   The executable mie_ms_solver, example1.out, example2.out are created. 
   The mie_ms_solver is the solver of Mie coefficients. 
   The example1.out is the executable of source code example1.c, it shows a simplest example using "emf_mie_ms". 
   The example2.out is the executable of source code example2.c, it shows a example of electromagnetic field intensity analysis.  
   
2. type './mie_ms_solver' with an argument of output datafile name.  
   For example, './mie_ms_solver ex.dat'. 
   The beam datafile "ipw.txt" (plane wave is defined) and the sphere datafile "msphr.txt" (radius 2.5 water droplet is defined) are used. 
   This executable calcluates Mie coefficients, outputs them to binary file with the specified file name.
   This program searches for a sphere datafile in current directory using the default datafile name "msphr.txt". 
   As a simple representation of the analysis model, the nodes used for the surface integral are output as point cloud data. 
   In this case, the file "ex.particles" is output, and the visualization result is "ex_particle.png" ( using ParaView ). 
   
3. type './example1.out' with an argument of datafile name.   
   For example, './example1.out ex.dat'. 
   This executable calculates electromagnetic field, radiation force and torque. This is the simplest example using this code.   
   
4. type './example2.out' with an argument of datafile name.  
   For example, './example2.out ex.dat'. 
   This executable calculates electromagnetic field intensity distributions, outputs them to text files.
   The I_example2.png is the visualization result of electromagnetic field intensity distributions, created by Gnuplot script gscript_example2.plt
   (converted eps to png by using ImageMagick).  

![point cloud data](ex_particles.png "nodes for surface integral (ex_particles.png)") 
![intensity distributions](I_example2.png "intensity distributions (I_example2.png)")

Please see exmie_src/emf_mie_ms.h for detail of functions, mfb_src/multi_fbeam.h for detail of incident fields. 
The mie_ms_solver is parallelized using OpenMP. The number of threads is controlled by the environment variable OMP_NUM_THREADS.
The additional analysis example of single sphere is in the folder analysis_sample1.


## Analysis sample of multi-spheres 

The analysis sample of multi-spheres is in the folder analysis_sample2. 
The multiple scattering between spheres can be analyzed by the executable mie_ms_solver.
The usage is the same as example code. 
For exampe, '../mie_ms_solver ex2.dat', '../example2.out ex2.dat' (run in analysis_sample2 folder). 
The I_example2.png in this folder is the visualization result of electromagnetic field intensity distributions.  

![point cloud data 2](analysis_sample2/ex2_particles.png "nodes for surface integral (analysis_sample2/ex_particles.png)") 
![intensity distributions 2](analysis_sample2/I_example2.png "intensity distributions (analysis_sample2/I_example2.png)")


## Analysis sample of radiation force  

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
6. Command-line driven graphing utility [gnuplot](http://www.gnuplot.info/)  
7. The electromagnetic field analysis program [multi_fbeam](https://github.com/akohta/multi_fbeam/)  
8. The data analysis and visualization application [ParaView](https://www.paraview.org/)  

The formula (12) in the Reference 2 ( z-component of radiation torque ) is misprinted. The following formula is correct.  
<img src="https://latex.codecogs.com/gif.latex?\frac{\left<N_z\right>}{a^3E_0^2}=-\frac{a}{8\pi}\sum_{l=1}^{\infty}\sum_{m=-l}^{l}l(l+1)m\left[\epsilon_{\mathrm{ext}}|a_{lm}|^2+|b_{lm}|^2+\Re(\epsilon_{\mathrm{ext}}a_{lm}A_{lm}^*+b_{lm}B_{lm}^*)\right]">.  
The first letter <img src="https://latex.codecogs.com/gif.latex?l"> in the term of sum is missed.
