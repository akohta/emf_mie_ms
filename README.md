# emf_mie_ms
This is the electromagnetic field analysis program for spherical particles. This is based on Mie scattering theory. 
This program can analyze multiple scattering between spheres by using iterative solutions. 
The electromagnetic field analysis program "multi_fbeam" is used to analyze incident field. 

## Usage of example code

1. type 'make' command to compile.  
   The executable mie_ms_solver, example1.out, example2.out are created.  
   
2. type './mie_ms_solver' with a argument of output datafile name. For example, './mie_ms_solver ex.dat'.  
   This executable calcluates coefficients, outputs them to binary file with specified name.  
   The sphere datafile 'msphr.txt' ( defined by default ) is searched in current directory and loaded automatically. 
   
3. type './example1.out' with a argument of datafile name. For example, './example1.out ex.dat'.  
   This executable ( source code is example1.c ) calculates electromagnetic field and radiation force.  
   This is the simplest example using this code. 
   
4. type './example2.out' with a argument of datafile name. For example, './example2.out ex.dat'.  
   This executable ( source code is example2.c ) calculates electric field intensity distributions, outputs them to text files.  
  The I_example2.pdf is the visualization result of intensity distributions, created by the Gnuplot script 'gscript_example2.plt'.
   
Please see 'exmie_src/emf_mie_ms.h' for detail of functions, 'mfb_src/multi_fbeam.h' for detail of incident fields. 


### Example of multi-spheres
The example of multi-spheres is in the folder 'multi_sphere_sample'. 
The multiple scattering between spheres can be analyzed by the executable 'mie_ms_solver'.
The usage is same as example code. 
For exampe, '../mie_ms_solver exms.dat', '../example2.out exms.dat' ( run in 'multi_sphere_sample' ). 
The I_example2.pdf in this folder is the visualization result of intensity distributions. 


### Appendix for radiation analysis  
The code in the folder 'appendix' is the radiation force analysis program for optical trapping. 
Usage is 'make' and './radiation_force.out'.
The radiation_force.pdf is the visualization result of radiation forces, created by the Gnuplot script 'gscript_radiation_force.plt'.


## References
1. Barton, J. P., D. R. Alexander, and S. A. Schaub. "Internal and near‐surface electromagnetic fields for a spherical particle irradiated by a focused laser beam." Journal of Applied Physics 64.4 (1988): 1632-1639.  
2. Barton, J. P., D. R. Alexander, and S. A. Schaub. "Theoretical determination of net radiation force and torque for a spherical particle illuminated by a focused laser beam." Journal of Applied Physics 66.10 (1989): 4594-4602.  
3. Barton, John P., et al. "Electromagnetic field for a beam incident on two adjacent spherical particles." Applied optics 30.33 (1991): 4706-4715.  
4. Abramowitz, Milton, and Irene A. Stegun, eds. Handbook of mathematical functions with formulas, graphs, and mathematical tables. Vol. 55. US Government printing office, 1948.  

The formula (12) in the Reference 2 ( z-component of radiation torque ) is misprinted. The following formula is correct.  
<img src="https://latex.codecogs.com/gif.latex?\frac{\left<N_z\right>}{a^3E_0^2}=-\frac{a}{8\pi}\sum_{l=1}^{\infty}\sum_{m=-l}^{l}l(l+1)m\left[\epsilon_{\mathrm{ext}}|a_{lm}|^2+|b_{lm}|^2+\Re(\epsilon_{\mathrm{ext}}a_{lm}A_{lm}^*+b_{lm}B_{lm}^*)\right]">.  
The first letter <img src="https://latex.codecogs.com/gif.latex?l"> in the term of sum is missed.



