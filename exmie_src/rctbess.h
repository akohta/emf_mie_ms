#if !defined RTCBESS_H
#define RTCBESS_H
// Calculate Riccati-Bessel functions and their derivatives
void rctjd(int n,double x,int *mn,double *rj,double *dj); 
/*Calculate Riccati-Bessel functions of first kind and their derivatives
        Input:   x     --- Argument of Riccati-Bessel function (Real number)
                 n     --- Order of Riccati-Bessel function
        Output:  rj[n] --- j_n(x)
                 dj[n] --- j_n(x)'  (derivative)
                 mn    --- Highest order calculated */
void rctyd(int n,double x,int *mn,double *ry,double *dy); 
/*Calculate Riccati-Bessel functions of second kind and their derivatives
       Input:   x     --- Argument of Riccati-Bessel function (Real number)
                n     --- Order of Riccati-Bessel function 
       Output:  ry[n] --- y_n(x)
                dy[n] --- y_n(x)' (derivative)
                mn    --- Highest order calculated */
void rctjc(int n,double _Complex z,int *mn,double _Complex *cj,double _Complex *dcj);
/*Calculate Riccati-Bessel functions of first kind and their derivatives 
  Input : z   --- Complex number */
void rctyc(int n,double _Complex z,int *mn,double _Complex *cy,double _Complex *dcy);
/*Calculate Riccati-Bessel functions of second kind and their derivatives
  Input : z   --- Complex number */
void rcth1d(int n,double x,int *mn,double _Complex *ch,double _Complex *dch);
/*Calculate Riccati-Bessel hankel type functions of first kind and their derivatives 
  xi_l^(1)(x)=psi_l(x)-I*chi_l(x) */
void rcth2d(int n,double x,int *mn,double _Complex *ch,double _Complex *dch);
/*Calculate Riccati-Bessel hankel type functions of second kind and their derivatives
  xi_l^(1)(x)=psi_l(x)+I*chi_l(x)  */
void rcth1c(int n,double _Complex z,int *mn,double _Complex *ch,double _Complex *dch);
/*Calculate Riccati-Bessel hankel type functions of first kind and their derivatives
  Input : z  --- Complex argument */
void rcth2c(int n,double _Complex z,int *mn,double _Complex *ch,double _Complex *dch);
/*Calculate Riccati-Bessel hankel type functions of second kind and their derivatives
  Input : z  --- Complex argument */

int msta1(double x,int mp);
int msta2(double x,int n,int mp);

#endif
