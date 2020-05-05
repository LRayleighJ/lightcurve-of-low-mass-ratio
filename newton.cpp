/*
Use the Newton-Raphson method to find images close to the source and lenses

the initial guess in newton's method for this problem was discussed at
https://ui.adsabs.harvard.edu/abs/1993ApJ...403..530W/abstract

-------------------------------------
variables introduction and example:
zlens[0]=complex(0.7, -0.05);
zlens[1]=complex(-0.4, 0.1);
zlens[2]=complex( 0.1, 0.5);
rs = 1e-6;
mlens[0]= 1e-6;
mlens[1]= 1e-6;
mlens[2]= 1-mlens[0]-mlens[1];

xs, ys: the position of the source

NLENS: number of lenses
-------------------------------------

add following lines to the VBBinaryLensingLibrary.h :

void findCloseImages(double mlens[], complex zlens[], double xs, double ys, int NLENS, complex close_zr[], int iter[]);
void newtonSolver(double mlens[], complex zlens[], double xs, double ys, int NLENS, complex *z,int *iter);

*/


#define ITER_MAX 50
#define rootEPS 1e-15
void newtonSolver(double mlens[], complex zlens[], double xs, double ys, int NLENS, complex *z, int *iter)
{
/* 
*z is the initial guess of the solution, after this function is done, *z is the corresponding solution
*/

  int j, k;
  int success = 0;
  complex delta;
  complex fz, dfdz, dfdzc, dz;
  complex zsc, zs;
  zs = complex(xs, ys);
  zsc = conj(zs);
  /* now refine it */
  for (k=0; k<ITER_MAX; k++) {
    fz = *z;
    dfdz = 1.0;
    dfdzc = 0.0;
    
    for (j=0; j<NLENS; j++) {
      fz = fz - mlens[j]/conj(*z-zlens[j]);
      // dfdzc += mlens[j]/cpow(conj(*z-zlens[j]), 2);
      dfdzc = dfdzc + mlens[j]/(conj(*z-zlens[j]) * conj(*z-zlens[j]));
    }
    
    delta = zs - fz;
    
    dz = (delta - dfdzc*conj(delta))/(1.0-pow(abs(dfdzc), 2));
    *z = *z + dz;
    
    if (abs(dz) < rootEPS) {
      success = 1;
      // fprintf(stderr, "<><><><> find solution, iteration = %d, x1, x2 = %.5f, %.5f\n", *k+1, (*z).re, (*z).im);
      break;
    }
  }
  *iter = k;

}


/* Use the Newton-Raphson method to find images close to the source and lenses */
void findCloseImages(double mlens[], complex zlens[], double xs, double ys, int NLENS, complex close_zr[], int iter[])
{
/*

close_zr[] is an  NLENS+1 dimensions arrry for saving the solved close solutions
close_zr[0] save the primary image
close_zr[1] ~ close_zr[NLENS] save the solved image corresponding to the NLENS masses

iter[] is an NLENS+1 dimensions array for saving how many iterations used to solve out the corresponding solution.

*/
  int i, k, it;
  complex sum;
  complex zs, zsc;
  complex z;

  zs = complex(xs, ys);
  zsc = conj(zs);

  /* primary image, we assume this has always been found */
  z = zs; 
  newtonSolver(mlens, zlens, xs, ys, NLENS, &z, &it);

  close_zr[0] = z;
  iter[0] = it;

  /* find close images to point lenses */
  for (i=0; i<NLENS; i++) {
    /* initial guess */
    sum = 0.0; z = zlens[i];
    for (k=0; k<NLENS; k++) {
      if (k == i) continue;
      sum = sum + mlens[k]/(zlens[i]-zlens[k]);
    }
    z = z + mlens[i]/(zlens[i]-zsc-sum);
    // z = z + mlens[i]/(conj(zlens[i])-zsc-sum);
    newtonSolver(mlens, zlens, xs, ys, NLENS, &z, &it);
    close_zr[i+1] = z;
    iter[i+1] = it;
  }
}
