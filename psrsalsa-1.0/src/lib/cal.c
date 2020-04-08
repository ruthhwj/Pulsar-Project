/*
If gsl_linalg_LU_invert crashes, this might be because inversion of
Mueller matrix failed. This happens for instance if the ellipticities
are opposite and or0-or1=90 deg.
 */

#include <math.h>
#include <string.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_randist.h"
#include <gsl/gsl_linalg.h>
#include "psrsalsa.h"

int output_FittedProfile_Of_Itterations;  // Variable to enable debug information at specific points in code

double refmodel_phase(double freq, long refmodel, verbose_definition verbose);
double refmodel_diffgain(double freq, long refmodel, verbose_definition verbose);
long phasemodels_names(int index);
long diffgainmodels_names(int index);

/* Calculates the Mueller matrix corresponding to the leakage terms
   (ellipticities). feedtype indicates if the receiver is linear or
   circular. This Mueller matrix is to fix the Stokes parameters if
   inverse is set. A distorted (described by the receiver parameters)
   observation will be restored to what it should look like. If
   inverse set to zero this Mueller matrix is to distort the Stokes
   parameters. A intrinsic signal will be distorted to what it will
   look like after propagation through (described by the receiver
   parameters). If ignore_dodgy_leakage is set, it makes the
   Mueller matrix zero (i.e. it will zap data) if the normalisation
   factor in the leakage matrix is abnormally high. This will only be
   in effect together with the reconstruct flag.

   Return:
   0 = dodgy solution, M=0
   1 = solution accepted, M is set
*/
int leakage_Mueller_matrix(float M[4][4], double el0, double el1, double or0, double or1, int feedtype, int inverse, int ignore_dodgy_leakage)
{
  double norm;
  double a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p;
  int i1, i2;

  a = cos(2.0*el0);
  b = sin(2.0*or0);
  c = cos(2.0*el1);
  d = sin(2.0*or1);
  e = sin(2.0*el1);
  f = sin(2.0*el0);
  g = cos(2.0*or0);
  h = cos(2.0*or1);

  i = sin(or0-or1);
  j = cos(el0+el1);
  k = cos(or0+or1);
  l = cos(el0-el1);
  m = sin(el1-el0);
  n = sin(or0+or1);
  o = cos(or0-or1);
  p = sin(el0+el1);

  if(inverse) {
    norm = 4.0/(2.0 + cos(2*(el0 - el1)) - cos(2*(el0 + el1)) +  2*cos(2*el0)*cos(2*el1)*cos(2*(or0 - or1)));
    //    printf("XXXXX %e (2 + %e - %e + %e) (%e %e %e %e)\n", norm, cos(2*(el0 - el1)), cos(2*(el0 + el1)),  2*cos(2*el0)*cos(2*el1)*cos(2*(or0 - or1)), el0*180.0/M_PI, el1*180.0/M_PI, or0*180.0/M_PI, or1*180.0/M_PI);
    // Check the normalisation constant if we want to ignore dodgy channels.
    if(ignore_dodgy_leakage && norm > 3.0) {
      for(i1 = 0; i1 < 4; i1++) {
	for(i2 = 0; i2 < 4; i2++) {
	  M[i1][i2] = 0.0;
	}
      }
      //      printwarning(0, "WARNING leakage_Mueller_matrix: Ignoring dodgy solution with norm=%f", norm);
      return 0;
    }
    //    printf("XXXXX norm = %lf\n", norm);
  }

  if(feedtype == FEEDTYPE_CIRCULAR) {  // Circular basis
    if(inverse == 0) {
      /*
1                            0.5*(cos(2*el0)*sin(2*or0) - cos(2*el1)*sin(2*or1))      0.5*(sin(2*el1)-sin(2*el0))    0.5*(cos(2*el0)*cos(2*or0)-cos(2*el1)*cos(2*or1))
sin(or0-or1)*cos(el0+el1)    cos(or0+or1)*cos(el0-el1)                                sin(or0-or1)*sin(el1-el0)      -sin(or0+or1)*cos(el1-el0)
cos(or0-or1)*sin(el1-el0)    sin(or0+or1)*sin(el0+el1)                                cos(or0-or1)*cos(el0+el1)      cos(or0+or1)*sin(el0+el1)
0                            0.5*(sin(2*or0)*cos(2*el0)+sin(2*or1)*cos(2*el1))        -0.5*(sin(2*el0)+sin(2*el1))   0.5*(cos(2*el0)*cos(2*or0)+cos(2*el1)*cos(2*or1))
     */

      M[0][0] = 1;
      M[0][1] = 0.5*(a*b - c*d);
      M[0][2] = 0.5*(e-f);
      M[0][3] = 0.5*(a*g-c*h);
      
      M[1][0] = i*j;
      M[1][1] = k*l;
      M[1][2] = i*m;
      M[1][3] = -n*l;

      M[2][0] = o*m;
      M[2][1] = n*p;
      M[2][2] = o*j;
      M[2][3] = k*p;
      
      M[3][0] = 0;
      M[3][1] = 0.5*(b*a+d*c);
      M[3][2] = -0.5*(f+e);
      M[3][3] = 0.5*(a*g+c*h);
    }else {
      M[0][0] = norm;
      M[0][1] = -norm*j*i;
      M[0][2] = -norm*m*o;
      M[0][3] = 0;
      
      M[1][0] = norm*0.5*(c*d - a*b);
      M[1][1] = norm*l*k;
      M[1][2] = norm*p*n;
      M[1][3] = norm*0.5*(a*b + c*d);
      
      M[2][0] = norm*0.5*(f - e);
      M[2][1] = norm*m*i;
      M[2][2] = norm*j*o;
      M[2][3] = -norm*0.5*(f + e);
      
      M[3][0] = norm*0.25*(2*c*h - cos(2.0*(el0-or0)) - cos(2.0*(el0+or0)));
      M[3][1] = -norm*l*n;
      M[3][2] = norm*p*k;
      M[3][3] = norm*0.5*(a*g + c*h);
    }
  }else if(feedtype == FEEDTYPE_LINEAR) {  // Linear basis
    if(inverse == 0) {
    /*
1                                0.5*(cos(2*el0)*cos(2*or0)-cos(2*el1)*cos(2*or1))  0.5*(sin(2*or0)*cos(2*el0)-sin(2*or1)*cos(2*el1))    0.5*(sin(2*el1)-sin(2*el0))
0                                0.5*(cos(2*el0)*cos(2*or0)+cos(2*el1)*cos(2*or1))  0.5*(sin(2*or0)*cos(2*el0)+sin(2*or1)*cos(2*el1))  -0.5*(sin(2*el1)+sin(2*el0))
sin(or0-or1)*cos(el0+el1)        -sin(or0+or1)*cos(el0-el1)                         cos(el0-el1)*cos(or0+or1)                          sin(or0-or1)*sin(el1-el0)
sin(el1-el0)*cos(or0-or1)        sin(el1+el0)*cos(or0+or1)                          sin(or0+or1)*sin(el0+el1)                          cos(el1+el0)*cos(or0-or1)
     */

      M[0][0] = 1;
      M[0][1] = 0.5*(a*g-c*h);
      M[0][2] = 0.5*(b*a-d*c);
      M[0][3] = 0.5*(e-f);

      M[1][0] = 0;
      M[1][1] = 0.5*(a*g+c*h);
      M[1][2] = 0.5*(b*a+d*c);
      M[1][3] = -0.5*(e+f);
      
      M[2][0] = i*j;
      M[2][1] = -n*l;
      M[2][2] = l*k;
      M[2][3] = i*m;
      
      M[3][0] = m*o;
      M[3][1] = p*k;
      M[3][2] = n*p;
      M[3][3] = j*o;
    }else {
      M[0][0] = norm;
      M[0][1] = 0;
      M[0][2] = -norm*j*i;
      M[0][3] = -norm*m*o;
      
      M[1][0] = norm*0.25*(2.0*c*h - cos(2.0*(el0-or0))- cos(2.0*(el0+or0)));
      M[1][1] = norm*0.5*(a*g + c*h);
      M[1][2] = -norm*l*n;
      M[1][3] = norm*p*k;
      
      M[2][0] = norm*0.5*(c*d - a*b);
      M[2][1] = norm*0.5*(a*b + c*d);
      M[2][2] = norm*l*k;
      M[2][3] = norm*p*n;
      
      M[3][0] = norm*0.5*(f - e);
      M[3][1] = -norm*0.5*(f + e);
      M[3][2] = norm*m*i;
      M[3][3] = norm*j*o;
    }
  }
  /*
  for(i1 = 0; i1 < 4; i1++) {
    printf("XXXXX ");
    for(i2 = 0; i2 < 4; i2++) {
      printf("%e ", M[i1][i2]);
    }
    printf("\n");
  }
  printf("\n");
  */
  return 1;
}

/* Calculates the Mueller matrix corresponding to the leakage terms
   (ellipticities). feedtype indicates if the receiver is linear or
   circular. This Mueller matrix is to fix the Stokes parameters if
   inverse is set. A distorted (described by the receiver parameters)
   observation will be restored to what it should look like. If
   inverse set to zero this Mueller matrix is to distort the Stokes
   parameters. A intrinsic signal will be distorted to what it will
   look like after propagation through (described by the receiver
   parameters). 
*/
void leakage_Mueller_matrix_old(float M[4][4], float el0, float el1, float or0, float or1, int feedtype, int inverse)
{
  float norm;
  //  int j, k;

  if(inverse) {
    norm = 4.0/(2.0 + cos(2*(el0 - el1)) - cos(2*(el0 + el1)) +  2*cos(2*el0)*cos(2*el1)*cos(2*(or0 - or1)));
  }

  if(feedtype == FEEDTYPE_CIRCULAR) {  // Circular basis
    if(inverse == 0) {
      /*
1                            0.5*(cos(2*el0)*sin(2*or0) - cos(2*el1)*sin(2*or1))      0.5*(sin(2*el1)-sin(2*el0))    0.5*(cos(2*el0)*cos(2*or0)-cos(2*el1)*cos(2*or1))
sin(or0-or1)*cos(el0+el1)    cos(or0+or1)*cos(el0-el1)                                sin(or0-or1)*sin(el1-el0)      -sin(or0+or1)*cos(el1-el0)
cos(or0-or1)*sin(el1-el0)    sin(or0+or1)*sin(el0+el1)                                cos(or0-or1)*cos(el0+el1)      cos(or0+or1)*sin(el0+el1)
0                            0.5*(sin(2*or0)*cos(2*el0)+sin(2*or1)*cos(2*el1))        -0.5*(sin(2*el0)+sin(2*el1))   0.5*(cos(2*el0)*cos(2*or0)+cos(2*el1)*cos(2*or1))
     */

      M[0][0] = 1;
      M[0][1] = 0.5*(cos(2.0*el0)*sin(2.0*or0) - cos(2.0*el1)*sin(2.0*or1));
      M[0][2] = 0.5*(sin(2.0*el1)-sin(2.0*el0));
      M[0][3] = 0.5*(cos(2.0*el0)*cos(2.0*or0)-cos(2.0*el1)*cos(2.0*or1));
      
      M[1][0] = sin(or0-or1)*cos(el0+el1);
      M[1][1] = cos(or0+or1)*cos(el0-el1);
      M[1][2] = sin(or0-or1)*sin(el1-el0);
      M[1][3] = -sin(or0+or1)*cos(el1-el0);
      
      M[2][0] = cos(or0-or1)*sin(el1-el0);
      M[2][1] = sin(or0+or1)*sin(el0+el1);
      M[2][2] = cos(or0-or1)*cos(el0+el1);
      M[2][3] = cos(or0+or1)*sin(el0+el1);
      
      M[3][0] = 0;
      M[3][1] = 0.5*(sin(2.0*or0)*cos(2.0*el0)+sin(2.0*or1)*cos(2.0*el1));
      M[3][2] = -0.5*(sin(2.0*el0)+sin(2.0*el1));
      M[3][3] = 0.5*(cos(2.0*el0)*cos(2.0*or0)+cos(2.0*el1)*cos(2.0*or1));
    }else {
      M[0][0] = norm;
      M[0][1] = norm*cos(el0+el1)*sin(or1-or0);
      M[0][2] = norm*sin(el0-el1)*cos(or0-or1);
      M[0][3] = 0;
      
      M[1][0] = norm*0.5*(cos(2.0*el1)*sin(2.0*or1) - cos(2.0*el0)*sin(2.0*or0));
      M[1][1] = norm*cos(el0-el1)*cos(or0+or1);
      M[1][2] = norm*sin(el0+el1)*sin(or0+or1);
      M[1][3] = norm*0.5*(cos(2.0*el0)*sin(2.0*or0) + cos(2.0*el1)*sin(2.0*or1));
      
      M[2][0] = norm*0.5*(sin(2.0*el0) - sin(2.0*el1));
      M[2][1] = norm*sin(el1-el0)*sin(or0-or1);
      M[2][2] = norm*cos(el0+el1)*cos(or0-or1);
      M[2][3] = -norm*0.5*(sin(2.0*el0) + sin(2.0*el1));
      
      M[3][0] = norm*0.25*(2*cos(2.0*el1)*cos(2.0*or1) - cos(2.0*(el0-or0))- cos(2.0*(el0+or0)));
      M[3][1] = -norm*cos(el0-el1)*sin(or0+or1);
      M[3][2] = norm*sin(el0+el1)*cos(or0+or1);
      M[3][3] = norm*0.5*(cos(2.0*el0)*cos(2.0*or0) + cos(2.0*el1)*cos(2.0*or1));
    }
  }else if(feedtype == FEEDTYPE_LINEAR) {  // Linear basis
    if(inverse == 0) {
    /*
1                                0.5*(cos(2*el0)*cos(2*or0)-cos(2*el1)*cos(2*or1))  0.5*(sin(2*or0)*cos(2*el0)-sin(2*or1)*cos(2*el1))    0.5*(sin(2*el1)-sin(2*el0))
0                                0.5*(cos(2*el0)*cos(2*or0)+cos(2*el1)*cos(2*or1))  0.5*(sin(2*or0)*cos(2*el0)+sin(2*or1)*cos(2*el1))  -0.5*(sin(2*el1)+sin(2*el0))
sin(or0-or1)*cos(el0+el1)        -sin(or0+or1)*cos(el0-el1)                         cos(el0-el1)*cos(or0+or1)                          sin(or0-or1)*sin(el1-el0)
sin(el1-el0)*cos(or0-or1)        sin(el1+el0)*cos(or0+or1)                          sin(or0+or1)*sin(el0+el1)                          cos(el1+el0)*cos(or0-or1)
     */

      M[0][0] = 1;
      M[0][1] = 0.5*(cos(2.0*el0)*cos(2.0*or0)-cos(2.0*el1)*cos(2.0*or1));
      M[0][2] = 0.5*(sin(2.0*or0)*cos(2.0*el0)-sin(2.0*or1)*cos(2.0*el1));
      M[0][3] = 0.5*(sin(2.0*el1)-sin(2.0*el0));

      M[1][0] = 0;
      M[1][1] = 0.5*(cos(2.0*el0)*cos(2.0*or0)+cos(2.0*el1)*cos(2.0*or1));
      M[1][2] = 0.5*(sin(2.0*or0)*cos(2.0*el0)+sin(2.0*or1)*cos(2.0*el1));
      M[1][3] = -0.5*(sin(2.0*el1)+sin(2.0*el0));
      
      M[2][0] = sin(or0-or1)*cos(el0+el1);
      M[2][1] = -sin(or0+or1)*cos(el0-el1);
      M[2][2] = cos(el0-el1)*cos(or0+or1);
      M[2][3] = sin(or0-or1)*sin(el1-el0);
      
      M[3][0] = sin(el1-el0)*cos(or0-or1);
      M[3][1] = sin(el1+el0)*cos(or0+or1);
      M[3][2] = sin(or0+or1)*sin(el0+el1);
      M[3][3] = cos(el1+el0)*cos(or0-or1);
    }else {
      M[0][0] = norm;
      M[0][1] = 0;
      M[0][2] = norm*cos(el0+el1)*sin(or1-or0);
      M[0][3] = norm*sin(el0-el1)*cos(or0-or1);
      
      M[1][0] = norm*0.25*(2.0*cos(2.0*el1)*cos(2.0*or1) - cos(2.0*(el0-or0))- cos(2.0*(el0+or0)));
      M[1][1] = norm*0.5*(cos(2.0*el0)*cos(2.0*or0) + cos(2.0*el1)*cos(2.0*or1));
      M[1][2] = -norm*cos(el0-el1)*sin(or0+or1);
      M[1][3] = norm*sin(el0+el1)*cos(or0+or1);
      
      M[2][0] = norm*0.5*(cos(2.0*el1)*sin(2.0*or1) - cos(2.0*el0)*sin(2.0*or0));
      M[2][1] = norm*0.5*(cos(2.0*el0)*sin(2.0*or0) + cos(2.0*el1)*sin(2.0*or1));
      M[2][2] = norm*cos(el0-el1)*cos(or0+or1);
      M[2][3] = norm*sin(el0+el1)*sin(or0+or1);
      
      M[3][0] = norm*0.5*(sin(2.0*el0) - sin(2.0*el1));
      M[3][1] = -norm*0.5*(sin(2.0*el0) + sin(2.0*el1));
      M[3][2] = norm*sin(el1-el0)*sin(or0-or1);
      M[3][3] = norm*cos(el0+el1)*cos(or0-or1);
    }
  }
  /*
  if(inverse) {
    norm = 4.0/(2.0 + cos(2*(el0 - el1)) - cos(2*(el0 + el1)) +  2*cos(2*el0)*cos(2*el1)*cos(2*(or0 - or1)));
    for(j = 0; j < 4; j++) {
      for(k = 0; k < 4; k++) {
	M[j][k] *= norm;
      }
    }
  }
  */
}

/* Calculates the Mueller matrix corresponding to a gain
   imbalance. feedtype indicates if the receiver is linear or
   circular. This Mueller matrix is to fix the Stokes parameters if
   inverse is set. A distorted (described by the receiver parameters)
   observation will be restored to what it should look like. If
   inverse set to zero this Mueller matrix is to distort the Stokes
   parameters. A intrinsic signal will be distorted to what it will
   look like after propagation through (described by the receiver
   parameters). If the gains of the individual polarization channels
   are g0 and g1 (in the Jones matrix), then G=sqrt(g0*g1) and
   gamma=g0/g1-1.  */
void gainImbalance_Mueller_matrix(float M[4][4], double G, double gamma, double phase, int feedtype, int inverse)
{
  double a, b, c, d;  
  double g0, g1;
  /*
    Calculate the gains of the two channels
    
    G     = sqrt(g0*g1)
    gamma = g0/g1-1
    
    g0 = g1*(gamma+1)
    G     = g1*sqrt(gamma+1) -> g0 = g1*(gamma+1.0)
  */
  g1 = G/sqrt(gamma+1.0);
  g0 = g1*(gamma+1.0);

  /* To get the inverse Mueller matrix, take the inverse of the
     gains. We also need the opposite sign of the differential phase,
     and the phase always appears with a factor two in the
     equations. */
  if(inverse) {
    g0 = 1.0/g0;
    g1 = 1.0/g1;
    phase *= 2.0;
  }else {
    phase *= -2.0;
  }

  /*
    printf("XXXXX G=%f    gamma=%f  phase=%f\n", G, gamma, phase);
    printf("XXXXX g0=%f   g1=%f\n", g0, g1);
  */

  /* Calculate the values of the non-zero elements of the matrix first */
  c = g0*g0;
  d = g1*g1;
  a = 0.5*(c+d);
  b = 0.5*(c-d);
  c = g0*g1;
  d = c*cos(phase);
  c *= sin(phase);

  if(feedtype == FEEDTYPE_CIRCULAR) {  // Circular feed
    M[0][0] = a;
    M[0][1] = 0;
    M[0][2] = 0;
    M[0][3] = b;

    M[1][0] = 0;
    M[1][1] = d;
    M[1][2] = -c;
    M[1][3] = 0;

    M[2][0] = 0;
    M[2][1] = c;
    M[2][2] = d;
    M[2][3] = 0;

    M[3][0] = b;
    M[3][1] = 0;
    M[3][2] = 0;
    M[3][3] = a;
  }else if(feedtype == FEEDTYPE_LINEAR) { // Linear feed
    M[0][0] = a;
    M[0][1] = b;
    M[0][2] = 0;
    M[0][3] = 0;

    M[1][0] = b;
    M[1][1] = a;
    M[1][2] = 0;
    M[1][3] = 0;

    M[2][0] = 0;
    M[2][1] = 0;
    M[2][2] = d;
    M[2][3] = -c;

    M[3][0] = 0;
    M[3][1] = 0;
    M[3][2] = c;
    M[3][3] = d;
  }
  /*
  int i1, i2;
  for(i1 = 0; i1 < 4; i1++) {
    printf("XXXXX: ");
    for(i2 = 0; i2 < 4; i2++) {
      printf("%e ", M[i1][i2]);
    }
    printf("\n");
  }
  printf("\n");
  */
}

/* Construct the Mueller matrix for the given receiver parameters. If
   reconstruct is set, the Mueller matrix will recover the intrinsic
   Stokes parameters from the recorded Stokes parameters. If
   reconstruct is not set it will distord the Stokes parameters in the
   same way the receiver would do. If ignore_leakage is set, the
   leakage (ellipticities/orientations) are ignored, which makes this
   function faster. If ignore_dodgy_leakage is set, it makes the
   Mueller matrix zero (i.e. it will zap data) if the normalisation
   factor in the leakage matrix is unrealisticly high. This will only
   be in effect together with the reconstruct flag. If normalise is
   set, the matrix is divided by the value of the first element,
   i.e. after applying the Mueller matrix to the data the rms of
   Stokes I in the noise is unchanged, assuming that the rms of the
   four Stokes parameters are equal to start off with. Normalise is
   only identical to setting the gain to 1 for an ideal receiver.

   Return:
   0 = dodgy solution, M=0
   1 = solution accepted, M is set   
*/
int constructMueller(float M[4][4], double G, double gamma, double phase, double el0, double el1, double or0, double or1, int feedtype, int reconstruct, int ignore_leakage, int ignore_dodgy_leakage, int normalise)
{
  int ret = 1;
  if(ignore_leakage) {
    gainImbalance_Mueller_matrix(M, G, gamma, phase, feedtype, reconstruct);
  }else {
    /* First calculate the Mueller matrix that distords the signal */
    float M1[4][4], M2[4][4];
    int j, k, l;
    gainImbalance_Mueller_matrix(M1, G, gamma, phase, feedtype, reconstruct);
    if(leakage_Mueller_matrix(M2, el0, el1, or0, or1, feedtype, reconstruct, ignore_dodgy_leakage) == 0) {
      //      printwarning(0, "WARNING constructMueller: solution is ignored");
      ret = 0;
    }
    // Multiply matrices. Order depends on if we calculating the inverse.
    if(reconstruct == 0) {
      for(j = 0; j < 4; j++) {
	for(k = 0; k < 4; k++) {
	  M[j][k] = 0;
	  for(l = 0; l < 4; l++) {
	    M[j][k] += M1[j][l]*M2[l][k];
	  }
	}
      }
    }else {
      for(j = 0; j < 4; j++) {
	for(k = 0; k < 4; k++) {
	  M[j][k] = 0;
	  for(l = 0; l < 4; l++) {
	    M[j][k] += M2[j][l]*M1[l][k];
	  }
	}
      }
    }
    /*
    int i1, i2;
    for(i1 = 0; i1 < 4; i1++) {
      printf("XXXXX: ");
      for(i2 = 0; i2 < 4; i2++) {
	printf("%e ", M[i1][i2]);
      }
      printf("     ");
      for(i2 = 0; i2 < 4; i2++) {
	printf("%e ", M1[i1][i2]);
      }
      printf("     ");
      for(i2 = 0; i2 < 4; i2++) {
	printf("%e ", M2[i1][i2]);
      }
      printf("\n");
    }
    printf("\n");
    */
  }
  if(normalise) {
    double varI, fac;
    int i, j;
    //Error propagation: varI is the expected variance of Stokes I after applying the Mueller matrix on data with unit variance for all four Stokes parameters.
    varI = M[0][0]*M[0][0] + M[0][1]*M[0][1] + M[0][2]*M[0][2] + M[0][3]*M[0][3];
    // Want this variance to remain unity, which means that the initial variance in Stokes I is unchanged. This can be achieved by multiplying the Mueller matrix with constant fac:
    //    1 = fac**2*(M[0][0]*M[0][0] + M[0][1]*M[0][1] + M[0][2]*M[0][2] + M[0][3]*M[0][3]);
    //    1 = fac**2*varI;
    fac = sqrt(1.0/varI);
    for(i = 0; i < 4; i++) {
      for(j = 0; j < 4; j++) {
	M[i][j] *= fac;
      }
    }
  }
  return ret;
}

/* S must be 4 Stokes parameters. The vector S is multiplied with the
   Mueller matrix and S is replaced with the result. */
void applyMueller(float M[4][4], float *S)
{
  float Inew, Qnew, Unew, Vnew;

  Inew = M[0][0]*S[0] + M[0][1]*S[1] + M[0][2]*S[2] + M[0][3]*S[3];
  Qnew = M[1][0]*S[0] + M[1][1]*S[1] + M[1][2]*S[2] + M[1][3]*S[3];
  Unew = M[2][0]*S[0] + M[2][1]*S[1] + M[2][2]*S[2] + M[2][3]*S[3];
  Vnew = M[3][0]*S[0] + M[3][1]*S[1] + M[3][2]*S[2] + M[3][3]*S[3];

  /*
  printf("applyMueller:\n");
  printf("%f * %f + %f * %f + %f * %f + %f * %f = %f\n", M[0][0],S[0] , M[0][1],S[1] , M[0][2],S[2] , M[0][3],S[3], Inew);
  printf("%f * %f + %f * %f + %f * %f + %f * %f = %f\n", M[1][0],S[0] , M[1][1],S[1] , M[1][2],S[2] , M[1][3],S[3], Qnew);
  printf("%f * %f + %f * %f + %f * %f + %f * %f = %f\n", M[2][0],S[0] , M[2][1],S[1] , M[2][2],S[2] , M[2][3],S[3], Unew);
  printf("%f * %f + %f * %f + %f * %f + %f * %f = %f\n", M[3][0],S[0] , M[3][1],S[1] , M[3][2],S[2] , M[3][3],S[3], Vnew);
  */

  S[0] = Inew;
  S[1] = Qnew;
  S[2] = Unew;
  S[3] = Vnew;
}

/*
   The Stokes parameters are passed through an imperfect quarter wave
   plate. The system is characterized with qw_phi (rotation of
   basis) and the phase delay qw_delay, with -45,90 representing an
   ideal quarter wave plate. This is assuming a circular basis.
 */
void applyQW(float *S, float qw_phi, float qw_delay)
{
  float M[4][4];
  int i, j;

  qw_phi *= M_PI/180.0;
  qw_delay *= M_PI/180.0;
  for(i = 0; i < 4; i++)
    for(j = 0; j < 4; j++)
      M[i][j] = 0;
  M[0][0] = 1;
  M[1][1] = -sin(2*qw_phi);
  M[1][2] = cos(qw_delay)*cos(2*qw_phi);
  M[1][3] = -sin(qw_delay)*cos(2*qw_phi);
  M[2][2] = sin(qw_delay);
  M[2][3] = cos(qw_delay);
  M[3][1] = cos(2*qw_phi);
  M[3][2] = cos(qw_delay)*sin(2*qw_phi);
  M[3][3] = -sin(qw_delay)*sin(2*qw_phi);

  /*
  printf("applyMueller:\n");
  printf("%f * %f + %f * %f + %f * %f + %f * %f\n", M[0][0],S[0] , M[0][1],S[1] , M[0][2],S[2] , M[0][3],S[3]);
  printf("%f * %f + %f * %f + %f * %f + %f * %f\n", M[1][0],S[0] , M[1][1],S[1] , M[1][2],S[2] , M[1][3],S[3]);
  printf("%f * %f + %f * %f + %f * %f + %f * %f\n", M[2][0],S[0] , M[2][1],S[1] , M[2][2],S[2] , M[2][3],S[3]);
  printf("%f * %f + %f * %f + %f * %f + %f * %f\n", M[3][0],S[0] , M[3][1],S[1] , M[3][2],S[2] , M[3][3],S[3]);
  */

  applyMueller(M, S);
}



/* S must be 4 Stokes parameters. The applied rotation (to Q and U)
   are to correct for rotations befor the receiver: for instance
   parallactic angle. */
void calibrate_applyQURotation(float *S, float phase)
{
  float Qnew, Unew;
  Qnew = S[1]*cos(2.0*phase) - S[2]*sin(2.0*phase);
  Unew = S[1]*sin(2.0*phase) + S[2]*cos(2.0*phase);
  S[1] = Qnew;
  S[2] = Unew;
}

/*
   0 - Gain
   1 - Differential gain [%%]
   2 - Differential phase [deg]
   3 - Elipticity 0 [deg]
   4 - Orientation 0 [deg]
   5 - Elipticity 1 [deg]
   6 - Orientation 1 [deg]
   7 - Chi-square
   8 - Nfree

   Returns the polarization nr in which the parameter is
   stored (value >= 0). On error the return value is negative:
      -1 if parameter is not part of the solution (consistent with gentype).
      -2 According to gentype the parameter exist, but it doens't match the nr of stored parameters
      -3 The input parameter is not recognized (beyond range?)
*/
int get_polnr_in_receiversolution(datafile_definition solution, int parameter)
{
  int polnr;
  if(parameter == 0) {           // Gain
    polnr = 0;
  }else if(parameter == 1) {     // Differential gain [%%]
    polnr = 1;
  }else if(parameter == 2) {     // Differential phase [deg]
    polnr = 2;
  }else if(parameter == 3) {     // Elipticity 0 [deg]
    polnr = 3;
  }else if(parameter == 4) {     // Orientation 0 [deg]
    polnr = 4;
  }else if(parameter == 5) {     // Elipticity 1 [deg]
    polnr = 5;
  }else if(parameter == 6) {     // Orientation 1 [deg]
    polnr = 6;
  }else if(parameter == 7) {     // Chi-square
    polnr = solution.NrPols-2;
    if(solution.gentype != GENTYPE_RECEIVERMODEL2) {
      //      fflush(stdout);
      //      printwarning(application.debug, "WARNING get_polnr_in_receiversolution: Requested parameter is not part of this type of receiver solution");
      return -1;
    }
  }else if(parameter == 8) {     // Nfree
    polnr = solution.NrPols-1;
    if(solution.gentype != GENTYPE_RECEIVERMODEL2) {
      //      fflush(stdout);
      //      printwarning(application.debug, "WARNING get_polnr_in_receiversolution: Requested parameter is not part of this type of receiver solution");
      return -1;
    }
  }else {
    fflush(stdout);
    printwarning(0, "WARNING get_polnr_in_receiversolution: Requested parameter not recognized");
    return -3;
  }
  if(polnr >= solution.NrPols) {
    //    fflush(stdout);
    //    printerror(debug, "ERROR get_polnr_in_receiversolution: Requested parameter does not appear to be part of the receiver solution");
    return -2;
  }
  return polnr;
}

/* Extract a parameter from the receiver model. Provide an
   analyticsolution != NULL if you want to subtract this from the
   data. The data is returned in xdata and ydata and yerrorbar, for
   which an array of nrpoints is allocated (so you should free
   it). xdata will contain the frequency, ydata the requested
   parameter. Only frequency channels with a valid solution are
   returned. The parameters are defined to be:

   0 - Gain
   1 - Differential gain [%%]
   2 - Differential phase [deg]
   3 - Elipticity 0 [deg]
   4 - Orientation 0 [deg]
   5 - Elipticity 1 [deg]
   6 - Orientation 1 [deg]
   7 - Chi-square
   8 - Nfree

   verbose determines the number of spaces before output.
   Return 1 = success, 0 on error 
*/
int extractReceiverParameter(int parameter, datafile_definition solution, long subintnr, analyticReceiverSolution_def *analyticsolution, float **xdata, float **ydata, float **yerrorbar, long *nrpoints, verbose_definition verbose)
{
  long i, f2, f;
  int polnr;
  double ysolution;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Extract parameter %d from receiver model\n", parameter);
  }
  if(solution.gentype != GENTYPE_RECEIVERMODEL2 && solution.gentype != GENTYPE_RECEIVERMODEL) {
    fflush(stdout);
    printerror(0, "ERROR extractReceiverParameter: File does not appear to be a receiver solution.");
    return 0;
  }
  if(solution.NrPols < 3) {
    fflush(stdout);
    printerror(0, "ERROR extractReceiverParameter: Expecteded at least 3 receiver parameters, only %ld present.", solution.NrPols);
    return 0;
  }
  if(solution.format != MEMORY_format) {
    fflush(stdout);
    printerror(0, "ERROR extractReceiverParameter: receiver model should be read into memory first");
    return 0;
  }

  *xdata = (float *)malloc(solution.NrFreqChan*sizeof(float));
  *ydata = (float *)malloc(solution.NrFreqChan*sizeof(float));
  *yerrorbar = (float *)malloc(solution.NrFreqChan*sizeof(float));
  if(*xdata == NULL || *ydata == NULL || *yerrorbar == NULL) {
    fflush(stdout);
    printerror(0, "ERROR extractReceiverParameter: Memory allocation error.");
    return 0;
  }

  polnr = get_polnr_in_receiversolution(solution, parameter);
  //  printf("XXXXXX paramater %d = %d\n", parameter, polnr);
  if(polnr == -1 || polnr == -2) {
    fflush(stdout);
    printerror(0, "ERROR extractReceiverParameter: Requested parameter is not part of this type of receiver solution");
    return 0;
  }
  if(polnr < 0) {
    fflush(stdout);
    printerror(0, "ERROR extractReceiverParameter: Requested parameter not recognized");
    return 0;
  }

  f2 = 0;
  for(f = 0; f < solution.NrFreqChan; f++) {
    /* If gain is zero, skip data point. */
    if(solution.data[solution.NrBins*(0+solution.NrPols*(f+subintnr*solution.NrFreqChan))+0] > 0) {
      (*ydata)[f2] = solution.data[solution.NrBins*(polnr+solution.NrPols*(f+subintnr*solution.NrFreqChan))+0];
      if(isnan((*ydata)[f2]) == 0 && isinf((*ydata)[f2]) == 0) {
	//	(*xdata)[f2] = solution.freq_cent-0.5*(solution.NrFreqChan-1.0)*get_channelbw(solution) + f*get_channelbw(solution);
	(*xdata)[f2] = get_weighted_channel_freq(solution, subintnr, f, verbose);
	if(solution.NrBins == 2)
	  (*yerrorbar)[f2] = solution.data[solution.NrBins*(polnr+solution.NrPols*(f+subintnr*solution.NrFreqChan))+1];
	else
	  (*yerrorbar)[f2] = 0;
	f2++;
      }
    }
  }
  *nrpoints = f2;

  for(f = 0; f < *nrpoints; f++) {
    if(parameter == 1) {
      (*ydata)[f] = 100.0*(exp(2.0*(*ydata)[f])-1.0);
      (*yerrorbar)[f] = 0;                                   // Should do some error propagation, but ignored for now
    }else if(parameter >= 2 && parameter <= 6) {
      (*ydata)[f] *= 180.0/M_PI;
      (*yerrorbar)[f] *= 180.0/M_PI;
    }else if(parameter == 7 || parameter == 8) {
      (*yerrorbar)[f] = 0;
    }
  }

  if(analyticsolution != NULL) {
    for(f = 0; f < *nrpoints; f++) {
      if(parameter == 0) {
	if(analyticsolution->gainmodel_defined) {
	  if(analyticReceiverSolution_gain(*analyticsolution, (*xdata)[f], &ysolution, verbose) == 0) {
	    fflush(stdout);
	    printerror(0, "ERROR extractReceiverParameter: Model calculation failed.");
	    return 0;
	  }
	  (*ydata)[f] -= ysolution;
	}
      }else if(parameter == 1) {
	if(analyticsolution->diffgainmodel_defined) {
	  if(analyticReceiverSolution_diffgain(*analyticsolution, (*xdata)[f], &ysolution, verbose) == 0) {
	    fflush(stdout);
	    printerror(0, "ERROR extractReceiverParameter: Model calculation failed.");
	    return 0;
	  }
	  (*ydata)[f] -= ysolution;
	}
      }else if(parameter == 2) {
	if(analyticsolution->phasemodel_defined) {
	  if(analyticReceiverSolution_diffphase(*analyticsolution, (*xdata)[f], &ysolution, verbose) == 0) {
	    fflush(stdout);
	    printerror(0, "ERROR extractReceiverParameter: Model calculation failed.");
	    return 0;
	  }
	  (*ydata)[f] -= ysolution;
	}
      }else if(parameter == 3 || parameter == 5) {
	polnr = 1;
	if(parameter == 5)
	  polnr = 2;
	if((polnr == 1 && analyticsolution->ell1model_defined) || (polnr == 2 && analyticsolution->ell2model_defined)) {
	  if(analyticReceiverSolution_ellipticity(*analyticsolution, polnr, (*xdata)[f], &ysolution) == 0) {
	    fflush(stdout);
	    printerror(0, "ERROR extractReceiverParameter: Model calculation failed.");
	    return 0;
	  }
	  (*ydata)[f] -= ysolution;
	}
      }else if(parameter == 4 || parameter == 6) {
	polnr = 1;
	if(parameter == 6)
	  polnr = 2;
	if((polnr == 1 && analyticsolution->or1model_defined) || (polnr == 2 && analyticsolution->or2model_defined)) {
	  if(analyticReceiverSolution_orientation(*analyticsolution, polnr, (*xdata)[f], &ysolution) == 0) {
	    fflush(stdout);
	    printerror(0, "ERROR extractReceiverParameter: Model calculation failed.");
	    return 0;
	  }
	  (*ydata)[f] -= ysolution;
	}
      }
    }
  }

  return 1;
}

/* Shows the receiver model on a pgplot device. The solution should be
   read into memory. Neat sets the style of the output:

   0 = Each parameter is shown seperately
   1 = More like a psrchive style output
   2 = No chi^2 and nfree, and all parameters on a single page

   If fixscale is set, the differential phase, orientation and
   ellipticity range is set to -90..90 and the differential gain range
   to -100..100. Set analyticsolution = NULL if you don't want an
   analytic solution to be superimposed. Set device to a pgplot
   device, or set to NULL if you want the user to specify the
   device. Set first when you want to open the device and last to
   close it. Multiple files can be plotted in this way. When plotting
   one solution, set both first and last. Multiple solutions can be
   show. If overlay is non-zero, the receiver solutions are
   overplotted in the same graph. If the title is set to NULL, the
   filename is used instead. The title supports header keywords set by
   str_replace_header_params(). If showmodel is set, the type of model
   used to fit the parameters is identified in the plots via
   text. Return 1 = success */
int showReceiverModel(datafile_definition *solution, long nrsolutions, int neat, int fixscale, analyticReceiverSolution_def *analyticsolution, char *device, int first, int last, int overlay, char *title, int showmodel, verbose_definition verbose)
{
  int firsttime, color, ok, parameter;
  long n, f, nrpoints, nrpoints_total, p, i, solnr;
  float *ydata, *xdata, *yerrorbar, *ydata_all, *xdata_all, *yerrorbar_all, ymin, ymax;
  char label[1000];
  pgplot_options_definition *pgplot_options;

  firsttime = 1;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Show receiver model\n");
  }
  if(nrsolutions > 1) {
    neat = 2;
    first = 1;
    last = 1;
    overlay = 1;
  }
  for(solnr = 0; solnr < nrsolutions; solnr++) {
    if(solution[solnr].gentype != GENTYPE_RECEIVERMODEL2 && solution[solnr].gentype != GENTYPE_RECEIVERMODEL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR showReceiverModel: File does not appear to be a receiver solution.");
      return 0;
    }
    if(solution[solnr].NrPols < 3) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR showReceiverModel: Expecteded at least 3 receiver parameters, only %ld present.", solution[solnr].NrPols);
      return 0;
    }
    if(solution[solnr].format != MEMORY_format) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR showReceiverModel: receiver model should be read into memory first");
      return 0;
    }
    if(solution[solnr].NrPols != solution[0].NrPols) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR showReceiverModel: Nr of parameters different in different input files.");
      return 0;
    }
    if(solution[solnr].NrSubints != solution[0].NrSubints) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR showReceiverModel: Nr of parameters different in different input files.");
      return 0;
    }
    if(solution[solnr].NrFreqChan != solution[0].NrFreqChan) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR showReceiverModel: Nr of frequency channels is different in different input files.");
      return 0;
    }
    if(solution[solnr].gentype != solution[0].gentype) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR showReceiverModel: Input files appear to be of different types.");
      return 0;
    }
  }
    
  pgplot_options = (pgplot_options_definition *)malloc(sizeof(pgplot_options_definition));
  ydata_all = malloc(nrsolutions*solution[0].NrFreqChan*sizeof(float));
  xdata_all = malloc(nrsolutions*solution[0].NrFreqChan*sizeof(float));
  yerrorbar_all = malloc(nrsolutions*solution[0].NrFreqChan*sizeof(float));
  if(pgplot_options == NULL || xdata_all == NULL || ydata_all == NULL || yerrorbar_all == NULL) {
    printerror(verbose.debug, "ERROR showReceiverModel: Memory allocation error");
    return 0;
  }
  pgplot_clear_options(pgplot_options);
    

  /* Loop over pulses and frequency channels etc. */
  for(n = 0; n < solution[0].NrSubints; n++) {
    for(p = 0; p < solution[0].NrPols; p++) {
      nrpoints_total = 0;
      if(p == 0) {
	parameter = 0;
      }else if(solution[0].gentype == GENTYPE_RECEIVERMODEL2 && p == solution[0].NrPols-2) {
	parameter = 7;
      }else if(solution[0].gentype == GENTYPE_RECEIVERMODEL2 && p == solution[0].NrPols-1) {
	parameter = 8;
      }else if(p == 1) {
	parameter = 1;
      }else {
	parameter = p;
      }
      if(overlay) {
	for(solnr = 0; solnr < nrsolutions; solnr++) {
	  if(extractReceiverParameter(parameter, solution[solnr], n, NULL, &xdata, &ydata, &yerrorbar, &nrpoints, verbose) == 0) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR showReceiverModel: extracting parameter failed");
	    return 0;
	  }
	  memcpy(&(xdata_all[nrpoints_total]), xdata, nrpoints*sizeof(float));
	  memcpy(&(ydata_all[nrpoints_total]), ydata, nrpoints*sizeof(float));	
	  memcpy(&(yerrorbar_all[nrpoints_total]), yerrorbar, nrpoints*sizeof(float));
	  nrpoints_total += nrpoints;
	  if(solnr != nrsolutions -1) {  // Don't free last solution, as it will be done again below
	    free(xdata);
	    free(ydata);
	    free(yerrorbar);
	  }
	}
      }else {
	if(extractReceiverParameter(parameter, solution[0], n, NULL, &xdata, &ydata, &yerrorbar, &nrpoints, verbose) == 0) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR showReceiverModel: extracting parameter failed");
	  return 0;
	}
      }

      if(parameter == 0) {
	sprintf(label, "Gain");
      }else if(parameter == 7) {
	sprintf(label, "Chi-square");
      }else if(parameter == 8) {
	sprintf(label, "Nfree");
      }else if(parameter == 1) {
	if(neat == 1)
	  sprintf(label, "Differential gain [%%]");
	else
	  sprintf(label, "Differential gain (gamma) [%%]");
      }else if(parameter == 2) {
	sprintf(label, "Differential phase [deg]");
      }else if(parameter == 3) {
	if(neat) sprintf(label, "Elipticity [deg]"); else sprintf(label, "Elipticity 0 [deg]");
      }else if(parameter == 4) {
	if(neat) sprintf(label, "Orientation [deg]"); else sprintf(label, "Orientation 0 [deg]");
      }else if(parameter == 5) {
	if(neat) label[0] = 0; else sprintf(label, "Elipticity 1 [deg]");
      }else if(parameter == 6) {
	if(neat) label[0] = 0; else sprintf(label, "Orientation 1 [deg]"); 
      }else {
	fflush(stdout); 
	printerror(verbose.debug, "ERROR showReceiverModel: Something wrong with parameter number.");
	return 0;
      }

      pgplot_clear_options(pgplot_options);
      strcpy(pgplot_options->box.xlabel, "Frequency [MHz]");
      strcpy(pgplot_options->box.ylabel, label);
      if(title == NULL) {
	strcpy(pgplot_options->box.title, solution[0].filename);
      }else {
	char *newtext;
	newtext = str_replace_header_params(solution[0], title, verbose);
	if(newtext != NULL) {
	  strcpy(pgplot_options->box.title, newtext);
	  free(newtext);
	}else {
	  pgplot_options->box.title[0] = 0;
	  fflush(stdout);
	  printwarning(0, "WARNING showReceiverModel: Error formatting title");
	}
      }
      if((parameter == 5 || parameter == 6) && neat)  // Prevent title to be printed twice in different colours
	pgplot_options->box.title[0] = 0;
      color = 1;
      if(device != NULL)
	sprintf(pgplot_options->viewport.plotDevice, device);
      pgplot_options->viewport.dontclose = 1;
      if(firsttime == 0 || first == 0)
	pgplot_options->viewport.dontopen = 1;
      ymin = ymax = 0;
      if(neat == 2) {
	pgplot_options->viewport.xsize = 0.53;
	if(p > 2)
	  pgplot_options->viewport.dxplot = 0.44;
	else
	  pgplot_options->viewport.dxplot = -0.05;
	if(p != 0)
	  pgplot_options->viewport.noclear = 1;
	pgplot_options->box.label_ch = 0.75;
	pgplot_options->box.box_labelsize = 0.75;
      }
      if(neat) {
	if(p < 3)
	  pgplot_options->viewport.ysize = 0.35;
	else
	  pgplot_options->viewport.ysize = 0.525;
	if(parameter == 0 || parameter == 3 || parameter == 5 || parameter == 7) {
	  pgplot_options->viewport.dyplot = 0;
	  if(parameter == 5)
	    pgplot_options->viewport.noclear = 1;
	  else if(parameter == 0 || neat == 1)
	    pgplot_options->viewport.noclear = 0;
	  pgplot_options->box.title[0] = 0;
	  pgplot_options->box.drawtitle = 0;
	  sprintf(pgplot_options->box.box_xopt, "bcnsti");
	  if(parameter == 3) {
	    sprintf(pgplot_options->box.box_yopt, "bnsti");
	    if(fixscale) {
	      ymin = -90;
	      ymax = 90;
	    }
	  }else if(parameter == 5) {
	    color = 2;
	    //	    sprintf(pgplotbox.box_yopt, "cmsti");
	    sprintf(pgplot_options->box.box_yopt, "csti");
	    if(fixscale) {
	      ymin = -90;
	      ymax = 90;
	    }
	  }
	}else if(parameter == 1 || parameter == 4 || parameter == 6 || parameter == 8) {
	  pgplot_options->viewport.dyplot = 0.8*pgplot_options->viewport.ysize;
	  pgplot_options->viewport.noclear = 1;
	  if(parameter == 1)
	    pgplot_options->box.drawtitle = 0;
	  pgplot_options->box.xlabel[0] = 0;
	  if(parameter == 1) {
	    if(fixscale) {
	      ymin = -100;
	      ymax = 100;
	    }
	  }else if(parameter == 4) {
	    sprintf(pgplot_options->box.box_yopt, "bnsti");
	    if(fixscale) {
	      ymin = -90;
	      ymax = 90;
	    }
	  }else if(parameter == 6) {
	    color = 2;
	    //	    sprintf(pgplotbox.box_yopt, "cmsti");
	    sprintf(pgplot_options->box.box_yopt, "csti");
	    if(fixscale) {
	      ymin = -90;
	      ymax = 90;
	    }
	  }
	  sprintf(pgplot_options->box.box_xopt, "bcsti");
	}else if(parameter == 2) {
	  pgplot_options->viewport.dyplot = 0.8*pgplot_options->viewport.ysize*2;
	  pgplot_options->viewport.noclear = 1;
	  pgplot_options->box.xlabel[0] = 0;
	  if(fixscale) {
	    ymin = -90;
	    ymax = 90;
	  }
	  sprintf(pgplot_options->box.box_xopt, "bcsti");
	}
      }

      // If overlaying data points: never do a clearpage except for the first receiver solution
      if(overlay && first == 0) {
	pgplot_options->viewport.noclear = 1;
      }
      if(solution[0].NrBins == 2) {
	if(overlay && nrsolutions > 1) {
	  //	  pgplotGraph1(pgplot_options, ydata_all, xdata_all, yerrorbar_all, nrpoints_total, xdata[0], xdata[nrpoints-1], 0, xdata[0], xdata[nrpoints-1], ymin, ymax, 0, 0, 1, -5, color, 1, NULL, -1, verbose);
	  pgplotGraph1(pgplot_options, ydata_all, xdata_all, yerrorbar_all, nrpoints_total, 0, 0, 0, 0, 0, ymin, ymax, 0, 0, 1, 1, -5, color, 1, NULL, -1, verbose);
	}else {
	  //	  pgplotGraph1(pgplot_options, ydata, xdata, yerrorbar, nrpoints, xdata[0], xdata[nrpoints-1], 0, xdata[0], xdata[nrpoints-1], ymin, ymax, 0, 0, 1, -5, color, 1, NULL, -1, verbose);
	  pgplotGraph1(pgplot_options, ydata, xdata, yerrorbar, nrpoints, 0, 0, 0, 0, 0, ymin, ymax, 0, 0, 1, 1, -5, color, 1, NULL, -1, verbose);
	}
      }else {
	if(overlay && nrsolutions > 1) {
	  //	  pgplotGraph1(pgplot_options, ydata_all, xdata_all, NULL, nrpoints_total, xdata[0], xdata[nrpoints-1], 0, xdata[0], xdata[nrpoints-1], ymin, ymax, 0, 0, 1, -5, color, 1, NULL, -1, verbose);
	  pgplotGraph1(pgplot_options, ydata_all, xdata_all, NULL, nrpoints_total, 0, 0, 0, 0, 0, ymin, ymax, 0, 0, 1, 1, -5, color, 1, NULL, -1, verbose);
	}else {
	  //	  pgplotGraph1(pgplot_options, ydata, xdata, NULL, nrpoints, xdata[0], xdata[nrpoints-1], 0, xdata[0], xdata[nrpoints-1], ymin, ymax, 0, 0, 1, -5, color, 1, NULL, -1, verbose);
	  pgplotGraph1(pgplot_options, ydata, xdata, NULL, nrpoints, 0, 0, 0, 0, 0, ymin, ymax, 0, 0, 1, 1, -5, color, 1, NULL, -1, verbose);
	}
      }

      if(analyticsolution != NULL) {
	double value;
	int havegraph, polnr;
	char text[1000];
	havegraph = 0;
	for(f = 0; f < nrpoints; f++) {
	  if(parameter == 0) {
	    if(analyticsolution->gainmodel_defined) {
	      if(analyticReceiverSolution_gain(*analyticsolution, xdata[f], &value, verbose) == 0) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR showReceiverModel: Model calculation failed.");
		return 0;
	      }
	      ydata[f] = value;
	      havegraph = 1;
	    }
	  }else if(parameter == 1) {
	    if(analyticsolution->diffgainmodel_defined) {
	      if(analyticReceiverSolution_diffgain(*analyticsolution, xdata[f], &value, verbose) == 0) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR showReceiverModel: Model calculation failed.");
		return 0;
	      }
	      ydata[f] = value;
	      havegraph = 1;
	      if(showmodel && f == 0) {
		sprintf(text, "ref=%ld", analyticsolution->diffgainRefmodel);
		ppgmtxt("r", 1.5, 0.05, 0, text);
	      }
	    }
	  }else if(parameter == 2) {
	    if(analyticsolution->phasemodel_defined) {
	      if(analyticReceiverSolution_diffphase(*analyticsolution, xdata[f], &value, verbose) == 0) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR showReceiverModel: Model calculation failed.");
		return 0;
	      }
	      ydata[f] = value;
	      havegraph = 1;
	      if(showmodel && f == 0) {
		sprintf(text, "ref=%ld", analyticsolution->phaseRefmodel);
		//		ppgmtxt("t", -2, 1.02, 0, text);
		ppgmtxt("r", 1.5, 0.05, 0, text);
	      }
	    }
	  }else if(parameter == 3 || parameter == 5) {
	    polnr = 1;
	    if(parameter == 5)
	      polnr = 2;
	    if((polnr == 1 && analyticsolution->ell1model_defined) || (polnr == 2 && analyticsolution->ell2model_defined)) {
	      if(analyticReceiverSolution_ellipticity(*analyticsolution, polnr, xdata[f], &value) == 0) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR showReceiverModel: Model calculation failed.");
		return 0;
	      }
	      ydata[f] = value;
	      havegraph = 1;
	    }
	  }else if(parameter == 4 || parameter == 6) {
	    polnr = 1;
	    if(parameter == 6)
	      polnr = 2;
	    if((polnr == 1 && analyticsolution->or1model_defined) || (polnr == 2 && analyticsolution->or2model_defined)) {
	      if(analyticReceiverSolution_orientation(*analyticsolution, polnr, xdata[f], &value) == 0) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR showReceiverModel: Model calculation failed.");
		return 0;
	      }
	      ydata[f] = value;
	      havegraph = 1;
	    }
	  }
	}
	if(havegraph) {
	  pgplot_options->viewport.dontopen = 1; // If plotting gain, this option is still set to 0
	  pgplot_options->viewport.noclear = 1;
	  pgplotGraph1(pgplot_options, ydata, xdata, NULL, nrpoints, xdata[0], xdata[nrpoints-1], 1, xdata[0], xdata[nrpoints-1], ymin, ymax, 0, 0, 0, 1, 0, color+2, 1, NULL, -1, verbose);
	}
      }
    
      if(parameter == 5 || parameter == 6) {
	    //	    sprintf(pgplotbox.box_yopt, "cmsti");
	sprintf(pgplot_options->box.box_yopt, "mti");
	pgplot_options->box.box_xopt[0] = 0;
	pgplot_options->box.xlabel[0] = 0;
	ppgsci(color);
	ppgslw(pgplot_options->box.box_lw);
	ppgsch(pgplot_options->box.box_labelsize);
	pgplot_drawbox(&(pgplot_options->box));
	ppgslw(1);
	ppgsch(1);
	ppgsci(1);
      }

      firsttime = 0;

      /* Find out if we need to wait for a key press */
      ok = 0;
      /* If we show all parameters one by one, wait for key press */
      if(neat == 0) {
	ok = 1;
      }else if(neat) {
	/* After gain, diff gain and diff phase, wait for key. */
	if(parameter == 2 && neat == 1) {
	  ok = 1;
	}
	/* And after the orientation/ellipticity plots */
	if((solution[0].gentype == GENTYPE_RECEIVERMODEL2 && parameter == solution[0].NrPols-3) || (solution[0].gentype == GENTYPE_RECEIVERMODEL && parameter == 6)) {
	  ok = 1;
	}
      }
      
      free(xdata);
      free(ydata);
      free(yerrorbar);

      /* Always wait after last parameter of the last file */
      if((p == solution[0].NrPols - 1 && (neat == 1 || neat == 2)) || (p == 6 && neat == 2)) {
	/* Unless for the last parameter of the last file */
	if(last)
	  ok = 0;
	else
	  ok = 1;
      }

      /* If device is selected, assume no key presses are desired */
      //      printf("ok=%d  device=%s\n", ok, device);
      if(ok && device == NULL) {
	printf("Press a key to continue\n");
	pgetch();
      }
      // Skip nfree and chi^2 plot if neat == 2
      if(neat == 2 && p == 6)
	break;

      }
  }

  free(xdata_all);
  free(ydata_all);
  free(yerrorbar_all);

  free(pgplot_options);
  if(last)
    ppgclos();
  return 1;
}



void clear_analyticReceiverSolution(analyticReceiverSolution_def *solution)
{
  solution->phasemodel = 0;
  solution->phaseRefmodel = 0;
  solution->diffgainRefmodel = 0;
  solution->phasemodel_defined = 0;
  solution->phasemodel_redchi2 = 0;
  solution->diffgainmodel = 0;
  solution->diffgainmodel_defined = 0;
  solution->diffgainmodel_redchi2 = 0;
  solution->gainmodel = 0;
  solution->gainmodel_defined = 0;
  solution->gainmodel_redchi2 = 0;
  solution->ell1model = 0;
  solution->ell1model_defined = 0;
  solution->ell1model_redchi2 = 0;
  solution->ell2model = 0;
  solution->ell2model_defined = 0;
  solution->ell2model_redchi2 = 0;
  solution->or1model = 0;
  solution->or1model_defined = 0;
  solution->or1model_redchi2 = 0;
  solution->or2model = 0;
  solution->or2model_defined = 0;
  solution->or2model_redchi2 = 0;
  solution->mjd = 0;
}

void printAnalyticReceiverSolution(analyticReceiverSolution_def solution, int indent)
{
  int i;

  for(i = 0; i < indent; i++)
    printf(" ");
  printf("Solutions in expressed as function of x=(freq[GHz]-%e)\n", 0.001*solution.reffreq);
  if(solution.gainmodel_defined) {
    if(solution.gainmodel == 1) {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Gain = %e  (red. chi2=%e)\n", solution.gainmodel_values[0], solution.gainmodel_redchi2);
    }else if(solution.gainmodel == 2) {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Gain = %e + %e * x  (red. chi2=%e)\n", solution.gainmodel_values[0], solution.gainmodel_values[1], solution.gainmodel_redchi2);
    }else {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Gain = %e + %e * x", solution.gainmodel_values[0], solution.gainmodel_values[1]);
      for(i = 0; i < solution.gainmodel - 2; i++) {
	printf(" + %e * x^%d", solution.gainmodel_values[2+i], i+2);
      }
      printf("  (red. chi2=%e)\n", solution.gainmodel_redchi2);
    }
  }

  if(solution.diffgainmodel_defined) {
    if(solution.diffgainmodel == 1) {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Differential gain [%%] = %e  (red. chi2=%e)\n", solution.diffgainmodel_values[0], solution.diffgainmodel_redchi2);
    }else if(solution.diffgainmodel == 2) {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Differential gain [%%] = %e + %e * x  (red. chi2=%e)\n", solution.diffgainmodel_values[0], solution.diffgainmodel_values[1], solution.diffgainmodel_redchi2);
    }else {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Differential gain [%%] = %e + %e * x", solution.diffgainmodel_values[0], solution.diffgainmodel_values[1]);
      for(i = 0; i < solution.diffgainmodel - 2; i++) {
	printf(" + %e * x^%d", solution.diffgainmodel_values[2+i], i+2);
      }
      printf("  (red. chi2=%e)\n", solution.diffgainmodel_redchi2);
    }
    for(i = 0; i < indent; i++)
      printf(" ");
    printf("Differential gain reference model = %ld\n", solution.diffgainRefmodel);
  }

  if(solution.phasemodel_defined) {
    if(solution.phasemodel == 1) {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Differential phase [deg] = %e (red. chi2=%e)\n", solution.phasemodel_values[0], solution.phasemodel_redchi2);
    }else if(solution.phasemodel == 2) {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Differential phase [deg] = %e + %e * x  (red. chi2=%e)\n", solution.phasemodel_values[0], solution.phasemodel_values[1], solution.phasemodel_redchi2);
    }else {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Differential phase [deg] = %e + %e * x ", solution.phasemodel_values[0], solution.phasemodel_values[1]);
      for(i = 0; i < solution.phasemodel - 2; i++) {
	printf(" + %e * x^%d", solution.phasemodel_values[2+i], i+2);
      }
      printf("  (red. chi2=%e)\n", solution.phasemodel_redchi2);
    }
    for(i = 0; i < indent; i++)
      printf(" ");
    printf("Differential phase reference model = %ld\n", solution.phaseRefmodel);
  }

  if(solution.ell1model_defined) {
    if(solution.ell1model == 1) {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Ellipticity 1 [deg] = %e  (red. chi2=%e)\n", solution.ell1model_values[0], solution.ell1model_redchi2);
    }else if(solution.ell1model == 2) {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Ellipticity 1 [deg] = %e + %e * x  (red. chi2=%e)\n", solution.ell1model_values[0], solution.ell1model_values[1], solution.ell1model_redchi2);
    }else {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Ellipticity 1 [deg] = %e + %e * x", solution.ell1model_values[0], solution.ell1model_values[1]);
      for(i = 0; i < solution.ell1model - 2; i++) {
	printf(" + %e * x^%d", solution.ell1model_values[2+i], i+2);
      }
      printf("  (red. chi2=%e)\n", solution.ell1model_redchi2);
    }
  }

  if(solution.ell2model_defined) {
    if(solution.ell2model == 1) {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Ellipticity 2 [deg] = %e  (red. chi2=%e)\n", solution.ell2model_values[0], solution.ell2model_redchi2);
    }else if(solution.ell2model == 2) {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Ellipticity 2 [deg] = %e + %e * x  (red. chi2=%e)\n", solution.ell2model_values[0], solution.ell2model_values[1], solution.ell2model_redchi2);
    }else {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Ellipticity 2 [deg] = %e + %e * x", solution.ell2model_values[0], solution.ell2model_values[1]);
      for(i = 0; i < solution.ell2model - 2; i++) {
	printf(" + %e * x^%d", solution.ell2model_values[2+i], i+2);
      }
      printf("  (red. chi2=%e)\n", solution.ell2model_redchi2);
    }
  }

  if(solution.or1model_defined) {
    if(solution.or1model == 1) {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Orientation 1 [deg] = %e  (red. chi2=%e)\n", solution.or1model_values[0], solution.or1model_redchi2);
    }else if(solution.or1model == 2) {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Orientation 1 [deg] = %e + %e * x  (red. chi2=%e)\n", solution.or1model_values[0], solution.or1model_values[1], solution.or1model_redchi2);
    }else {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Orientation 1 [deg] = %e + %e * x", solution.or1model_values[0], solution.or1model_values[1]);
      for(i = 0; i < solution.or1model - 2; i++) {
	printf(" + %e * x^%d", solution.or1model_values[2+i], i+2);
      }
      printf("  (red. chi2=%e)\n", solution.or1model_redchi2);
    }
  }

  if(solution.or2model_defined) {
    if(solution.or2model == 1) {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Orientation 2 [deg] = %e  (red. chi2=%e)\n", solution.or2model_values[0], solution.or2model_redchi2);
    }else if(solution.or2model == 2) {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Orientation 2 [deg] = %e + %e * x  (red. chi2=%e)\n", solution.or2model_values[0], solution.or2model_values[1], solution.or2model_redchi2);
    }else {
      for(i = 0; i < indent; i++)
	printf(" ");
      printf("Orientation 2 [deg] = %e + %e * x", solution.or2model_values[0], solution.or2model_values[1]);
      for(i = 0; i < solution.or2model - 2; i++) {
	printf(" + %e * x^%d", solution.or2model_values[2+i], i+2);
      }
      printf("  (red. chi2=%e)\n", solution.or2model_redchi2);
    }
  }

}


/*
  Writes analytic template to file

  Returns 1 on success, 0 on error
 */
int writeAnalyticReceiverSolution(char *filename, analyticReceiverSolution_def solution)
{
  FILE *fout;
  int version, i;

  version = 1;

  fout = fopen(filename, "wb");
  if(fout == NULL) {
    fflush(stdout);
    printerror(0, "ERROR writeAnalyticReceiverSolution: Cannot open %s", filename);
    return 0;
  }
  fprintf(fout, "AnalyticReceiverSolution\n");
  fprintf(fout, "version %d\n", version);
  fprintf(fout, "mjd %lf\n", solution.mjd);

  if(solution.gainmodel_defined) {
    fprintf(fout, "total_gain_model %d =", solution.gainmodel);
    for(i = 0; i < solution.gainmodel; i++) {
      fprintf(fout, " %e", solution.gainmodel_values[i]);
    }
    fprintf(fout, " %e ", solution.gainmodel_redchi2);
    if(solution.gainmodel == 1) {
      fprintf(fout, "(value red_chi2)\n");
    }else if(solution.gainmodel == 2) {
      fprintf(fout, "(offset slope_per_GHz red_chi2)\n");
    }else {
      fprintf(fout, "(offset slope_per_GHz");
      for(i = 0; i < solution.gainmodel - 2; i++) {
	fprintf(fout, " coeff[GHz^-%d]", i+2);
      }
      fprintf(fout, " red_chi2)\n");
    }
  }

  if(solution.diffgainmodel_defined) {
    fprintf(fout, "differential_gain_refmodel %ld\n", solution.diffgainRefmodel);
    fprintf(fout, "differential_gain_model %d =", solution.diffgainmodel);
    for(i = 0; i < solution.diffgainmodel; i++) {
      fprintf(fout, " %e", solution.diffgainmodel_values[i]);
    }
    fprintf(fout, " %e ", solution.diffgainmodel_redchi2);
    if(solution.diffgainmodel == 1) {
      fprintf(fout, "(value[%%] red_chi2)\n");
    }else if(solution.diffgainmodel == 2) {
      fprintf(fout, "(offset[%%] slope[%%/GHz] red_chi2)\n");
    }else {
      fprintf(fout, "(offset[%%] slope[%%/GHz]");
      for(i = 0; i < solution.diffgainmodel - 2; i++) {
	fprintf(fout, " coeff[%%/GHz^%d]", i+2);
      }
      fprintf(fout, " red_chi2)\n");
    }
  }

  if(solution.phasemodel_defined) {
    fprintf(fout, "phaserefmodel %ld\n", solution.phaseRefmodel);
    fprintf(fout, "phasemodel %d =", solution.phasemodel);
    for(i = 0; i < solution.phasemodel; i++) {
      fprintf(fout, " %e", solution.phasemodel_values[i]);
    }
    fprintf(fout, " %e ", solution.phasemodel_redchi2);
    if(solution.phasemodel == 1) {
      fprintf(fout, "(value[deg] red_chi2)\n");
    }else if(solution.phasemodel == 2) {
      fprintf(fout, "(offset[deg] slope[deg/GHz] red_chi2)\n");
    }else {
      fprintf(fout, "(offset[deg] slope[deg/GHz]");
      for(i = 0; i < solution.diffgainmodel - 2; i++) {
	fprintf(fout, " coeff[deg/GHz^%d]", i+2);
      }
      fprintf(fout, " red_chi2)\n");
    }
  }

  if(solution.ell1model_defined) {
    fprintf(fout, "ellipticity1_model %d =", solution.ell1model);
    for(i = 0; i < solution.ell1model; i++) {
      fprintf(fout, " %e", solution.ell1model_values[i]);
    }
    fprintf(fout, " %e ", solution.ell1model_redchi2);
    if(solution.ell1model == 1) {
      fprintf(fout, "(value[deg] red_chi2)\n");
    }else if(solution.ell1model == 2) {
      fprintf(fout, "(offset[deg] slope[deg/GHz] red_chi2)\n");
    }else {
      fprintf(fout, "(offset[deg] slope[deg/GHz]");
      for(i = 0; i < solution.ell1model - 2; i++) {
	fprintf(fout, " coeff[deg/GHz^%d]", i+2);
      }
      fprintf(fout, " red_chi2)\n");
    }
  }

  if(solution.ell2model_defined) {
    fprintf(fout, "ellipticity2_model %d =", solution.ell2model);
    for(i = 0; i < solution.ell2model; i++) {
      fprintf(fout, " %e", solution.ell2model_values[i]);
    }
    fprintf(fout, " %e ", solution.ell2model_redchi2);
    if(solution.ell2model == 1) {
      fprintf(fout, "(value[deg] red_chi2)\n");
    }else if(solution.ell2model == 2) {
      fprintf(fout, "(offset[deg] slope[deg/GHz] red_chi2)\n");
    }else {
      fprintf(fout, "(offset[deg] slope[deg/GHz]");
      for(i = 0; i < solution.ell2model - 2; i++) {
	fprintf(fout, " coeff[deg/GHz^%d]", i+2);
      }
      fprintf(fout, " red_chi2)\n");
    }
  }

  if(solution.or1model_defined) {
    fprintf(fout, "orientation1_model %d =", solution.or1model);
    for(i = 0; i < solution.or1model; i++) {
      fprintf(fout, " %e", solution.or1model_values[i]);
    }
    fprintf(fout, " %e ", solution.or1model_redchi2);
    if(solution.or1model == 1) {
      fprintf(fout, "(value[deg] red_chi2)\n");
    }else if(solution.or1model == 2) {
      fprintf(fout, "(offset[deg] slope[deg/GHz] red_chi2)\n");
    }else {
      fprintf(fout, "(offset[deg] slope[deg/GHz]");
      for(i = 0; i < solution.or1model - 2; i++) {
	fprintf(fout, " coeff[deg/GHz^%d]", i+2);
      }
      fprintf(fout, " red_chi2)\n");
    }
  }

  if(solution.or2model_defined) {
    fprintf(fout, "orientation2_model %d =", solution.or2model);
    for(i = 0; i < solution.or2model; i++) {
      fprintf(fout, " %e", solution.or2model_values[i]);
    }
    fprintf(fout, " %e ", solution.or2model_redchi2);
    if(solution.or2model == 1) {
      fprintf(fout, "(value[deg] red_chi2)\n");
    }else if(solution.or2model == 2) {
      fprintf(fout, "(offset[deg] slope[deg/GHz] red_chi2)\n");
    }else {
      fprintf(fout, "(offset[deg] slope[deg/GHz]");
      for(i = 0; i < solution.or2model - 2; i++) {
	fprintf(fout, " coeff[deg/GHz^%d]", i+2);
      }
      fprintf(fout, " red_chi2)\n");
    }
  }

  fclose(fout);
  return 1;
}

/*
  Calculates the differential phase (degrees) at frequency freq (MHz)
  according to the analytic receiver solution.

  Return 0 on error, 1 if success
 */
int analyticReceiverSolution_diffphase(analyticReceiverSolution_def solution, double freq, double *phase, verbose_definition verbose)
{
  freq = 0.001*(freq-solution.reffreq);
  if(solution.phasemodel_defined == 0) {
    fflush(stdout);
    printerror(0, "ERROR analyticReceiverSolution_diffphase: No differential phase model appears to be defined");
    return 0;
  }
  if(solution.phasemodel == 1) { // constant: phase[deg] = [0]
    *phase = solution.phasemodel_values[0];
    *phase = derotate_180_small_double(*phase);
  }else if(solution.phasemodel == 2) { // linear dep. with freq: phase[deg] = [0] + [1]*freq
    *phase = solution.phasemodel_values[0] + freq*solution.phasemodel_values[1];
    *phase = derotate_180_small_double(*phase);
  }else if(solution.phasemodel == 3) {
    double freq2;
    freq2 = freq*freq;
    *phase = solution.phasemodel_values[0] + freq*solution.phasemodel_values[1] + freq2*solution.phasemodel_values[2];
    *phase = derotate_180_small_double(*phase);
  }else if(solution.phasemodel == 4) {
    double freq2, freq3;
    freq2 = freq*freq;
    freq3 = freq2*freq;
    *phase = solution.phasemodel_values[0] + freq*solution.phasemodel_values[1] + freq2*solution.phasemodel_values[2] + freq3*solution.phasemodel_values[3];
    *phase = derotate_180_small_double(*phase);
  }else if(solution.phasemodel == 5) {
    double freq2, freq3, freq4;
    freq2 = freq*freq;
    freq3 = freq2*freq;
    freq4 = freq2*freq2;
    *phase = solution.phasemodel_values[0] + freq*solution.phasemodel_values[1] + freq2*solution.phasemodel_values[2] + freq3*solution.phasemodel_values[3] + freq4*solution.phasemodel_values[4];
    *phase = derotate_180_small_double(*phase);
  }else if(solution.phasemodel == 6) {
    double freq2, freq3, freq4, freq5;
    freq2 = freq*freq;
    freq3 = freq2*freq;
    freq4 = freq2*freq2;
    freq5 = freq3*freq2;
    *phase = solution.phasemodel_values[0] + freq*solution.phasemodel_values[1] + freq2*solution.phasemodel_values[2] + freq3*solution.phasemodel_values[3] + freq4*solution.phasemodel_values[4] + freq5*solution.phasemodel_values[5];
    *phase = derotate_180_small_double(*phase);
  }else if(solution.phasemodel == 7) {
    double freq2, freq3, freq4, freq5, freq6;
    freq2 = freq*freq;
    freq3 = freq2*freq;
    freq4 = freq2*freq2;
    freq5 = freq3*freq2;
    freq6 = freq3*freq3;
    *phase = solution.phasemodel_values[0] + freq*solution.phasemodel_values[1] + freq2*solution.phasemodel_values[2] + freq3*solution.phasemodel_values[3] + freq4*solution.phasemodel_values[4] + freq5*solution.phasemodel_values[5] + freq6*solution.phasemodel_values[6];
    *phase = derotate_180_small_double(*phase);
  }else if(solution.phasemodel == 8) {
    double freq2, freq3, freq4, freq5, freq6;
    freq2 = freq*freq;
    freq3 = freq2*freq;
    freq4 = freq2*freq2;
    freq5 = freq3*freq2;
    freq6 = freq3*freq3;
    *phase = solution.phasemodel_values[0] + freq*solution.phasemodel_values[1] + freq2*solution.phasemodel_values[2] + freq3*solution.phasemodel_values[3] + freq4*solution.phasemodel_values[4] + freq5*solution.phasemodel_values[5] + freq6*solution.phasemodel_values[6] + freq6*freq*solution.phasemodel_values[7];
    *phase = derotate_180_small_double(*phase);
  }else if(solution.phasemodel == 9) {
    double freq2, freq3, freq4, freq5, freq6;
    freq2 = freq*freq;
    freq3 = freq2*freq;
    freq4 = freq2*freq2;
    freq5 = freq3*freq2;
    freq6 = freq3*freq3;
    *phase = solution.phasemodel_values[0] + freq*solution.phasemodel_values[1] + freq2*solution.phasemodel_values[2] + freq3*solution.phasemodel_values[3] + freq4*solution.phasemodel_values[4] + freq5*solution.phasemodel_values[5] + freq6*solution.phasemodel_values[6] + freq6*freq*solution.phasemodel_values[7] + freq6*freq2*solution.phasemodel_values[8];
    *phase = derotate_180_small_double(*phase);
  }else if(solution.phasemodel == 10) {
    double freq2, freq3, freq4, freq5, freq6;
    freq2 = freq*freq;
    freq3 = freq2*freq;
    freq4 = freq2*freq2;
    freq5 = freq3*freq2;
    freq6 = freq3*freq3;
    *phase = solution.phasemodel_values[0] + freq*solution.phasemodel_values[1] + freq2*solution.phasemodel_values[2] + freq3*solution.phasemodel_values[3] + freq4*solution.phasemodel_values[4] + freq5*solution.phasemodel_values[5] + freq6*solution.phasemodel_values[6] + freq6*freq*solution.phasemodel_values[7] + freq6*freq2*solution.phasemodel_values[8] + freq6*freq3*solution.phasemodel_values[9];
    *phase = derotate_180_small_double(*phase);
  }else if(solution.phasemodel == 11) {
    double freq2, freq3, freq4, freq5, freq6;
    freq2 = freq*freq;
    freq3 = freq2*freq;
    freq4 = freq2*freq2;
    freq5 = freq3*freq2;
    freq6 = freq3*freq3;
    *phase = solution.phasemodel_values[0] + freq*solution.phasemodel_values[1] + freq2*solution.phasemodel_values[2] + freq3*solution.phasemodel_values[3] + freq4*solution.phasemodel_values[4] + freq5*solution.phasemodel_values[5] + freq6*solution.phasemodel_values[6] + freq6*freq*solution.phasemodel_values[7] + freq6*freq2*solution.phasemodel_values[8] + freq6*freq3*solution.phasemodel_values[9] + freq6*freq4*solution.phasemodel_values[10];
    *phase = derotate_180_small_double(*phase);
  }else {
    fflush(stdout);
    printerror(0, "ERROR analyticReceiverSolution_diffphase: Phase model %d is not implemented", solution.phasemodel);
    return 0;
  }

  if(solution.phaseRefmodel != 0) {
    // Check if phase model is defined.
    int i, ok;
    long model;
    ok = 0;
    for(i = 1; i < 1000; i++) {
      model = phasemodels_names(i);
      if(model == -1)                   // Reached end of list
	break;
      if(model == solution.phaseRefmodel)
	ok = 1;
    }
    if(ok) {
      *phase += refmodel_phase(freq, solution.phaseRefmodel, verbose);
      *phase = derotate_180_small_double(*phase);
    }else {
      fflush(stdout);
      printerror(0, "ERROR analyticReceiverSolution_diffphase: Phase reference model %ld is not implemented", solution.phaseRefmodel);
      return 0;
    }
  }
  return 1;
}


/*
  Calculates the differential gain (percentage) at frequency freq (MHz)
  according to the analytic receiver solution.

  Return 0 on error, 1 if success
 */
int analyticReceiverSolution_diffgain(analyticReceiverSolution_def solution, double freq, double *diffgain, verbose_definition verbose)
{
  freq = 0.001*(freq-solution.reffreq);
  if(solution.diffgainmodel_defined == 0) {
    fflush(stdout);
    printerror(0, "ERROR analyticReceiverSolution_diffgain: No differential gain model appears to be defined");
    return 0;
  }
  if(solution.diffgainmodel == 1) { // Constant: diffgain = [0]
    *diffgain = solution.diffgainmodel_values[0];
  }else if(solution.diffgainmodel == 2) { // linear dep. with freq: diffgain = [0] + [1]*freq
    *diffgain = solution.diffgainmodel_values[0] + freq*solution.diffgainmodel_values[1];
  }else if(solution.diffgainmodel == 3) { // quadratic dep. with freq: diffgain = [0] + [1]*freq + [2]*freq^2
    *diffgain = solution.diffgainmodel_values[0] + freq*solution.diffgainmodel_values[1] + freq*freq*solution.diffgainmodel_values[2];
  }else if(solution.diffgainmodel == 4) {
    double freq2;
    freq2 = freq*freq;
    *diffgain = solution.diffgainmodel_values[0] + freq*solution.diffgainmodel_values[1] + freq2*solution.diffgainmodel_values[2] + freq*freq2*solution.diffgainmodel_values[3];
  }else if(solution.diffgainmodel == 5) {
    double freq2;
    freq2 = freq*freq;
    *diffgain = solution.diffgainmodel_values[0] + freq*solution.diffgainmodel_values[1] + freq2*solution.diffgainmodel_values[2] + freq*freq2*solution.diffgainmodel_values[3] + freq2*freq2*solution.diffgainmodel_values[4];
  }else if(solution.diffgainmodel == 6) {
    double freq2, freq3;
    freq2 = freq*freq;
    freq3 = freq2*freq;
    *diffgain = solution.diffgainmodel_values[0] + freq*solution.diffgainmodel_values[1] + freq2*solution.diffgainmodel_values[2] + freq3*solution.diffgainmodel_values[3] + freq2*freq2*solution.diffgainmodel_values[4] + freq2*freq3*solution.diffgainmodel_values[5];
  }else if(solution.diffgainmodel == 7) {
    double freq2, freq3, freq4;
    freq2 = freq*freq;
    freq3 = freq2*freq;
    freq4 = freq2*freq2;
    *diffgain = solution.diffgainmodel_values[0] + freq*solution.diffgainmodel_values[1] + freq2*solution.diffgainmodel_values[2] + freq3*solution.diffgainmodel_values[3] + freq4*solution.diffgainmodel_values[4] + freq2*freq3*solution.diffgainmodel_values[5] + freq3*freq3*solution.diffgainmodel_values[6];
  }else if(solution.diffgainmodel == 8) {
    double freq2, freq3, freq4;
    freq2 = freq*freq;
    freq3 = freq2*freq;
    freq4 = freq2*freq2;
    *diffgain = solution.diffgainmodel_values[0] + freq*solution.diffgainmodel_values[1] + freq2*solution.diffgainmodel_values[2] + freq3*solution.diffgainmodel_values[3] + freq4*solution.diffgainmodel_values[4] + freq2*freq3*solution.diffgainmodel_values[5] + freq3*freq3*solution.diffgainmodel_values[6] + freq3*freq4*solution.diffgainmodel_values[7];
  }else if(solution.diffgainmodel == 9) {
    double freq2, freq3, freq4, freq5;
    freq2 = freq*freq;
    freq3 = freq2*freq;
    freq4 = freq2*freq2;
    freq5 = freq3*freq2;
    *diffgain = solution.diffgainmodel_values[0] + freq*solution.diffgainmodel_values[1] + freq2*solution.diffgainmodel_values[2] + freq3*solution.diffgainmodel_values[3] + freq4*solution.diffgainmodel_values[4] + freq5*solution.diffgainmodel_values[5] + freq3*freq3*solution.diffgainmodel_values[6] + freq3*freq4*solution.diffgainmodel_values[7] + freq4*freq4*solution.diffgainmodel_values[8];
  }else if(solution.diffgainmodel == 10) {
    double freq2, freq3, freq4, freq5;
    freq2 = freq*freq;
    freq3 = freq2*freq;
    freq4 = freq2*freq2;
    freq5 = freq3*freq2;
    *diffgain = solution.diffgainmodel_values[0] + freq*solution.diffgainmodel_values[1] + freq2*solution.diffgainmodel_values[2] + freq3*solution.diffgainmodel_values[3] + freq4*solution.diffgainmodel_values[4] + freq5*solution.diffgainmodel_values[5] + freq3*freq3*solution.diffgainmodel_values[6] + freq3*freq4*solution.diffgainmodel_values[7] + freq4*freq4*solution.diffgainmodel_values[8] + freq4*freq5*solution.diffgainmodel_values[9];
  }else if(solution.diffgainmodel == 11) {
    double freq2, freq3, freq4, freq5, freq6;
    freq2 = freq*freq;
    freq3 = freq2*freq;
    freq4 = freq2*freq2;
    freq5 = freq3*freq2;
    freq6 = freq3*freq3;
    *diffgain = solution.diffgainmodel_values[0] + freq*solution.diffgainmodel_values[1] + freq2*solution.diffgainmodel_values[2] + freq3*solution.diffgainmodel_values[3] + freq4*solution.diffgainmodel_values[4] + freq5*solution.diffgainmodel_values[5] + freq6*solution.diffgainmodel_values[6] + freq3*freq4*solution.diffgainmodel_values[7] + freq4*freq4*solution.diffgainmodel_values[8] + freq4*freq5*solution.diffgainmodel_values[9] + freq5*freq5*solution.diffgainmodel_values[10];
  }else {
    fflush(stdout);
    printerror(0, "ERROR analyticReceiverSolution_diffgain: Differential gain model %d is not implemented", solution.diffgainmodel);
    return 0;
  }
  if(solution.diffgainRefmodel != 0) {
    // Check if phase model is defined.
    int i, ok;
    long model;
    ok = 0;
    for(i = 1; i < 1000; i++) {
      model = diffgainmodels_names(i);
      if(model == -1)                   // Reached end of list
	break;
      if(model == solution.diffgainRefmodel)
	ok = 1;
    }
    if(ok) {
      *diffgain += refmodel_diffgain(freq, solution.diffgainRefmodel, verbose);
    }else {
      fflush(stdout);
      printerror(0, "ERROR analyticReceiverSolution_diffgain: Differential gain reference model %ld is not implemented", solution.diffgainRefmodel);
      return 0;
    }
  }
  return 1;
}

/*
  Calculates the gain (factor) at frequency freq (MHz)
  according to the analytic receiver solution.

  Return 0 on error, 1 if success
 */
int analyticReceiverSolution_gain(analyticReceiverSolution_def solution, double freq, double *gain, verbose_definition verbose)
{
  freq = 0.001*(freq-solution.reffreq);
  if(solution.gainmodel_defined == 0) {
    fflush(stdout);
    printerror(0, "ERROR analyticReceiverSolution_gain: No gain model appears to be defined");
    return 0;
  }
  if(solution.gainmodel == 1) { // Constant: gain = [0]
    *gain = solution.gainmodel_values[0];
    return 1;
  }else if(solution.gainmodel == 2) { // linear dep. with freq: gain = [0] + [1]*freq
    *gain = solution.gainmodel_values[0] + freq*solution.gainmodel_values[1];
    return 1;
  }else if(solution.gainmodel == 3) { 
    double freq2 = freq*freq;
    *gain = solution.gainmodel_values[0] + freq*solution.gainmodel_values[1] + freq2*solution.gainmodel_values[2];
    return 1;
  }else if(solution.gainmodel == 4) { 
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    *gain = solution.gainmodel_values[0] + freq*solution.gainmodel_values[1] + freq2*solution.gainmodel_values[2] + freq3*solution.gainmodel_values[3];
    return 1;
  }else if(solution.gainmodel == 5) { 
    double freq2 = freq*freq;
    double freq3 = freq2*freq; 
    double freq4 = freq2*freq2;
    *gain = solution.gainmodel_values[0] + freq*solution.gainmodel_values[1] + freq2*solution.gainmodel_values[2] + freq3*solution.gainmodel_values[3] + freq4*solution.gainmodel_values[4];
    return 1;
  }else if(solution.gainmodel == 6) { 
    double freq2 = freq*freq;
    double freq3 = freq2*freq; 
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    *gain = solution.gainmodel_values[0] + freq*solution.gainmodel_values[1] + freq2*solution.gainmodel_values[2] + freq3*solution.gainmodel_values[3] + freq4*solution.gainmodel_values[4] + freq5*solution.gainmodel_values[5];
    return 1;
  }else if(solution.gainmodel == 7) { 
    double freq2 = freq*freq;
    double freq3 = freq2*freq; 
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    double freq6 = freq3*freq3;
    *gain = solution.gainmodel_values[0] + freq*solution.gainmodel_values[1] + freq2*solution.gainmodel_values[2] + freq3*solution.gainmodel_values[3] + freq4*solution.gainmodel_values[4] + freq5*solution.gainmodel_values[5] + freq6*solution.gainmodel_values[6];
    return 1;
  }else if(solution.gainmodel == 8) { 
    double freq2 = freq*freq;
    double freq3 = freq2*freq; 
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    double freq6 = freq3*freq3;
    double freq7 = freq4*freq3;
    *gain = solution.gainmodel_values[0] + freq*solution.gainmodel_values[1] + freq2*solution.gainmodel_values[2] + freq3*solution.gainmodel_values[3] + freq4*solution.gainmodel_values[4] + freq5*solution.gainmodel_values[5] + freq6*solution.gainmodel_values[6] + freq7*solution.gainmodel_values[7];
    return 1;
  }else if(solution.gainmodel == 9) { 
    double freq2 = freq*freq;
    double freq3 = freq2*freq; 
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    double freq6 = freq3*freq3;
    double freq7 = freq4*freq3;
    double freq8 = freq4*freq4;
    *gain = solution.gainmodel_values[0] + freq*solution.gainmodel_values[1] + freq2*solution.gainmodel_values[2] + freq3*solution.gainmodel_values[3] + freq4*solution.gainmodel_values[4] + freq5*solution.gainmodel_values[5] + freq6*solution.gainmodel_values[6] + freq7*solution.gainmodel_values[7] + freq8*solution.gainmodel_values[8];
    return 1;
  }else if(solution.gainmodel == 10) { 
    double freq2 = freq*freq;
    double freq3 = freq2*freq; 
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    double freq6 = freq3*freq3;
    double freq7 = freq4*freq3;
    double freq8 = freq4*freq4;
    double freq9 = freq5*freq4;
    *gain = solution.gainmodel_values[0] + freq*solution.gainmodel_values[1] + freq2*solution.gainmodel_values[2] + freq3*solution.gainmodel_values[3] + freq4*solution.gainmodel_values[4] + freq5*solution.gainmodel_values[5] + freq6*solution.gainmodel_values[6] + freq7*solution.gainmodel_values[7] + freq8*solution.gainmodel_values[8] + freq9*solution.gainmodel_values[9];
    return 1;
  }else if(solution.gainmodel == 11) { 
    double freq2 = freq*freq;
    double freq3 = freq2*freq; 
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    double freq6 = freq3*freq3;
    double freq7 = freq4*freq3;
    double freq8 = freq4*freq4;
    double freq9 = freq5*freq4;
    double freq10 = freq5*freq5;
    *gain = solution.gainmodel_values[0] + freq*solution.gainmodel_values[1] + freq2*solution.gainmodel_values[2] + freq3*solution.gainmodel_values[3] + freq4*solution.gainmodel_values[4] + freq5*solution.gainmodel_values[5] + freq6*solution.gainmodel_values[6] + freq7*solution.gainmodel_values[7] + freq8*solution.gainmodel_values[8] + freq9*solution.gainmodel_values[9] + freq10*solution.gainmodel_values[10];
    return 1;
  }else {
    fflush(stdout);
    printerror(0, "ERROR analyticReceiverSolution_gain: Gain model %d is not implemented", solution.gainmodel);
    return 0;
  }
  return 1;
}

/*
  Calculates the ellipticity (degrees) at frequency freq (MHz)
  according to the analytic receiver solution. polnr sets which of the
  two polarizations (1 or 2) should be evaluated.

  Return 0 on error, 1 if success
 */
int analyticReceiverSolution_ellipticity(analyticReceiverSolution_def solution, int polnr, double freq, double *ell)
{
  freq = 0.001*(freq-solution.reffreq);
  if((polnr == 1 && solution.ell1model_defined == 0) || (polnr == 2 && solution.ell2model_defined == 0)) {
    fflush(stdout);
    printerror(0, "ERROR analyticReceiverSolution_ellipticity: No ellipticity model appears to be defined for channel %d", polnr);
    return 0;
  }
  if((polnr == 1 && solution.ell1model == 1) || (polnr == 2 && solution.ell2model == 1)) { // Constant = [0]
    if(polnr == 1) {
      *ell = solution.ell1model_values[0];
      *ell = derotate_180_small_double(*ell);
    }else {
      *ell = solution.ell2model_values[0];
      *ell = derotate_180_small_double(*ell);
    }
    return 1;
  }else if((polnr == 1 && solution.ell1model == 2) || (polnr == 2 && solution.ell2model == 2)) { // linear dep. with freq: ell[deg] = [0] + [1]*freq
    if(polnr == 1) {
      *ell = solution.ell1model_values[0] + freq*solution.ell1model_values[1];
      *ell = derotate_180_small_double(*ell);
    }else {
      *ell = solution.ell2model_values[0] + freq*solution.ell2model_values[1];
      *ell = derotate_180_small_double(*ell);
    }
    return 1;
  }else if((polnr == 1 && solution.ell1model == 3) || (polnr == 2 && solution.ell2model == 3)) {
    double freq2 = freq*freq;
    if(polnr == 1) {
      *ell = solution.ell1model_values[0] + freq*solution.ell1model_values[1] + freq2*solution.ell1model_values[2];
      *ell = derotate_180_small_double(*ell);
    }else {
      *ell = solution.ell2model_values[0] + freq*solution.ell2model_values[1] + freq2*solution.ell2model_values[2];
      *ell = derotate_180_small_double(*ell);
    }
    return 1;
  }else if((polnr == 1 && solution.ell1model == 4) || (polnr == 2 && solution.ell2model == 4)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    if(polnr == 1) {
      *ell = solution.ell1model_values[0] + freq*solution.ell1model_values[1] + freq2*solution.ell1model_values[2] + freq3*solution.ell1model_values[3];
      *ell = derotate_180_small_double(*ell);
    }else {
      *ell = solution.ell2model_values[0] + freq*solution.ell2model_values[1] + freq2*solution.ell2model_values[2] + freq3*solution.ell2model_values[3];
      *ell = derotate_180_small_double(*ell);
    }
    return 1;
  }else if((polnr == 1 && solution.ell1model == 5) || (polnr == 2 && solution.ell2model == 5)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    double freq4 = freq2*freq2;
    if(polnr == 1) {
      *ell = solution.ell1model_values[0] + freq*solution.ell1model_values[1] + freq2*solution.ell1model_values[2] + freq3*solution.ell1model_values[3] + freq4*solution.ell1model_values[4];
      *ell = derotate_180_small_double(*ell);
    }else {
      *ell = solution.ell2model_values[0] + freq*solution.ell2model_values[1] + freq2*solution.ell2model_values[2] + freq3*solution.ell2model_values[3] + freq4*solution.ell2model_values[4];
      *ell = derotate_180_small_double(*ell);
    }
    return 1;
  }else if((polnr == 1 && solution.ell1model == 6) || (polnr == 2 && solution.ell2model == 6)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    if(polnr == 1) {
      *ell = solution.ell1model_values[0] + freq*solution.ell1model_values[1] + freq2*solution.ell1model_values[2] + freq3*solution.ell1model_values[3] + freq4*solution.ell1model_values[4] + freq5*solution.ell1model_values[5];
      *ell = derotate_180_small_double(*ell);
    }else {
      *ell = solution.ell2model_values[0] + freq*solution.ell2model_values[1] + freq2*solution.ell2model_values[2] + freq3*solution.ell2model_values[3] + freq4*solution.ell2model_values[4] + freq5*solution.ell2model_values[5];
      *ell = derotate_180_small_double(*ell);
    }
    return 1;
  }else if((polnr == 1 && solution.ell1model == 7) || (polnr == 2 && solution.ell2model == 7)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    double freq6 = freq3*freq3;
    if(polnr == 1) {
      *ell = solution.ell1model_values[0] + freq*solution.ell1model_values[1] + freq2*solution.ell1model_values[2] + freq3*solution.ell1model_values[3] + freq4*solution.ell1model_values[4] + freq5*solution.ell1model_values[5] + freq6*solution.ell1model_values[6];
      *ell = derotate_180_small_double(*ell);
    }else {
      *ell = solution.ell2model_values[0] + freq*solution.ell2model_values[1] + freq2*solution.ell2model_values[2] + freq3*solution.ell2model_values[3] + freq4*solution.ell2model_values[4] + freq5*solution.ell2model_values[5] + freq6*solution.ell2model_values[6];
      *ell = derotate_180_small_double(*ell);
    }
    return 1;
  }else if((polnr == 1 && solution.ell1model == 8) || (polnr == 2 && solution.ell2model == 8)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    double freq6 = freq3*freq3;
    double freq7 = freq4*freq3;
    if(polnr == 1) {
      *ell = solution.ell1model_values[0] + freq*solution.ell1model_values[1] + freq2*solution.ell1model_values[2] + freq3*solution.ell1model_values[3] + freq4*solution.ell1model_values[4] + freq5*solution.ell1model_values[5] + freq6*solution.ell1model_values[6] + freq7*solution.ell1model_values[7];
      *ell = derotate_180_small_double(*ell);
    }else {
      *ell = solution.ell2model_values[0] + freq*solution.ell2model_values[1] + freq2*solution.ell2model_values[2] + freq3*solution.ell2model_values[3] + freq4*solution.ell2model_values[4] + freq5*solution.ell2model_values[5] + freq6*solution.ell2model_values[6] + freq7*solution.ell2model_values[7];
      *ell = derotate_180_small_double(*ell);
    }
    return 1;
  }else if((polnr == 1 && solution.ell1model == 9) || (polnr == 2 && solution.ell2model == 9)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    double freq6 = freq3*freq3;
    double freq7 = freq4*freq3;
    double freq8 = freq4*freq4;
    if(polnr == 1) {
      *ell = solution.ell1model_values[0] + freq*solution.ell1model_values[1] + freq2*solution.ell1model_values[2] + freq3*solution.ell1model_values[3] + freq4*solution.ell1model_values[4] + freq5*solution.ell1model_values[5] + freq6*solution.ell1model_values[6] + freq7*solution.ell1model_values[7] + freq8*solution.ell1model_values[8];
      *ell = derotate_180_small_double(*ell);
    }else {
      *ell = solution.ell2model_values[0] + freq*solution.ell2model_values[1] + freq2*solution.ell2model_values[2] + freq3*solution.ell2model_values[3] + freq4*solution.ell2model_values[4] + freq5*solution.ell2model_values[5] + freq6*solution.ell2model_values[6] + freq7*solution.ell2model_values[7] + freq8*solution.ell2model_values[8];
      *ell = derotate_180_small_double(*ell);
    }
    return 1;
  }else if((polnr == 1 && solution.ell1model == 10) || (polnr == 2 && solution.ell2model == 10)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    double freq6 = freq3*freq3;
    double freq7 = freq4*freq3;
    double freq8 = freq4*freq4;
    double freq9 = freq5*freq4;
    if(polnr == 1) {
      *ell = solution.ell1model_values[0] + freq*solution.ell1model_values[1] + freq2*solution.ell1model_values[2] + freq3*solution.ell1model_values[3] + freq4*solution.ell1model_values[4] + freq5*solution.ell1model_values[5] + freq6*solution.ell1model_values[6] + freq7*solution.ell1model_values[7] + freq8*solution.ell1model_values[8] + freq9*solution.ell1model_values[9];
      *ell = derotate_180_small_double(*ell);
    }else {
      *ell = solution.ell2model_values[0] + freq*solution.ell2model_values[1] + freq2*solution.ell2model_values[2] + freq3*solution.ell2model_values[3] + freq4*solution.ell2model_values[4] + freq5*solution.ell2model_values[5] + freq6*solution.ell2model_values[6] + freq7*solution.ell2model_values[7] + freq8*solution.ell2model_values[8] + freq9*solution.ell2model_values[9];
      *ell = derotate_180_small_double(*ell);
    }
    return 1;
  }else if((polnr == 1 && solution.ell1model == 11) || (polnr == 2 && solution.ell2model == 11)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    double freq6 = freq3*freq3;
    double freq7 = freq4*freq3;
    double freq8 = freq4*freq4;
    double freq9 = freq5*freq4;
    double freq10 = freq5*freq5;
    if(polnr == 1) {
      *ell = solution.ell1model_values[0] + freq*solution.ell1model_values[1] + freq2*solution.ell1model_values[2] + freq3*solution.ell1model_values[3] + freq4*solution.ell1model_values[4] + freq5*solution.ell1model_values[5] + freq6*solution.ell1model_values[6] + freq7*solution.ell1model_values[7] + freq8*solution.ell1model_values[8] + freq9*solution.ell1model_values[9] + freq10*solution.ell1model_values[10];
      *ell = derotate_180_small_double(*ell);
    }else {
      *ell = solution.ell2model_values[0] + freq*solution.ell2model_values[1] + freq2*solution.ell2model_values[2] + freq3*solution.ell2model_values[3] + freq4*solution.ell2model_values[4] + freq5*solution.ell2model_values[5] + freq6*solution.ell2model_values[6] + freq7*solution.ell2model_values[7] + freq8*solution.ell2model_values[8] + freq9*solution.ell2model_values[9] + freq10*solution.ell2model_values[10];
      *ell = derotate_180_small_double(*ell);
    }
    return 1;
  }else {
    fflush(stdout);
    if(polnr == 1) {
      printerror(0, "ERROR analyticReceiverSolution_ellipticity: Ellipticity model %d is not implemented", solution.ell1model);
    }else {
      printerror(0, "ERROR analyticReceiverSolution_ellipticity: Ellipticity model %d is not implemented", solution.ell2model);
    }
    return 0;
  }
  return 1;
}

/*
  Calculates the orientation (degrees) at frequency freq (MHz)
  according to the analytic receiver solution. polnr sets which of the
  two polarizations (1 or 2) should be evaluated.

  Return 0 on error, 1 if success
 */
int analyticReceiverSolution_orientation(analyticReceiverSolution_def solution, int polnr, double freq, double *or)
{
  freq = 0.001*(freq-solution.reffreq);
  if((polnr == 1 && solution.or1model_defined == 0) || (polnr == 2 && solution.or2model_defined == 0)) {
    fflush(stdout);
    printerror(0, "ERROR analyticReceiverSolution_orientation: No orientation model appears to be defined for channel %d", polnr);
    return 0;
  }
  if((polnr == 1 && solution.or1model == 1) || (polnr == 2 && solution.or2model == 1)) { // Constant = [0]
    if(polnr == 1) {
      *or = solution.or1model_values[0];
      *or = derotate_180_small_double(*or);
    }else {
      *or = solution.or2model_values[0];
      *or = derotate_180_small_double(*or);
    }
    return 1;
  }else if((polnr == 1 && solution.or1model == 2) || (polnr == 2 && solution.or2model == 2)) { // linear dep. with freq: or[deg] = [0] + [1]*freq
    if(polnr == 1) {
      *or = solution.or1model_values[0] + freq*solution.or1model_values[1];
      *or = derotate_180_small_double(*or);
    }else {
      *or = solution.or2model_values[0] + freq*solution.or2model_values[1];
      *or = derotate_180_small_double(*or);
    }
    return 1;
  }else if((polnr == 1 && solution.or1model == 3) || (polnr == 2 && solution.or2model == 3)) {
    double freq2 = freq*freq;
    if(polnr == 1) {
      *or = solution.or1model_values[0] + freq*solution.or1model_values[1] + freq2*solution.or1model_values[2];
      *or = derotate_180_small_double(*or);
    }else {
      *or = solution.or2model_values[0] + freq*solution.or2model_values[1] + freq2*solution.or2model_values[2];
      *or = derotate_180_small_double(*or);
    }
    return 1;
  }else if((polnr == 1 && solution.or1model == 4) || (polnr == 2 && solution.or2model == 4)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    if(polnr == 1) {
      *or = solution.or1model_values[0] + freq*solution.or1model_values[1] + freq2*solution.or1model_values[2] + freq3*solution.or1model_values[3];
      *or = derotate_180_small_double(*or);
    }else {
      *or = solution.or2model_values[0] + freq*solution.or2model_values[1] + freq2*solution.or2model_values[2] + freq3*solution.or2model_values[3];
      *or = derotate_180_small_double(*or);
    }
    return 1;
  }else if((polnr == 1 && solution.or1model == 5) || (polnr == 2 && solution.or2model == 5)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    double freq4 = freq2*freq2;
    if(polnr == 1) {
      *or = solution.or1model_values[0] + freq*solution.or1model_values[1] + freq2*solution.or1model_values[2] + freq3*solution.or1model_values[3] + freq4*solution.or1model_values[4];
      *or = derotate_180_small_double(*or);
    }else {
      *or = solution.or2model_values[0] + freq*solution.or2model_values[1] + freq2*solution.or2model_values[2] + freq3*solution.or2model_values[3] + freq4*solution.or2model_values[4];
      *or = derotate_180_small_double(*or);
    }
    return 1;
  }else if((polnr == 1 && solution.or1model == 6) || (polnr == 2 && solution.or2model == 6)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    if(polnr == 1) {
      *or = solution.or1model_values[0] + freq*solution.or1model_values[1] + freq2*solution.or1model_values[2] + freq3*solution.or1model_values[3] + freq4*solution.or1model_values[4] + freq5*solution.or1model_values[5];
      *or = derotate_180_small_double(*or);
    }else {
      *or = solution.or2model_values[0] + freq*solution.or2model_values[1] + freq2*solution.or2model_values[2] + freq3*solution.or2model_values[3] + freq4*solution.or2model_values[4] + freq5*solution.or2model_values[5];
      *or = derotate_180_small_double(*or);
    }
    return 1;
  }else if((polnr == 1 && solution.or1model == 7) || (polnr == 2 && solution.or2model == 7)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    double freq6 = freq3*freq3;
    if(polnr == 1) {
      *or = solution.or1model_values[0] + freq*solution.or1model_values[1] + freq2*solution.or1model_values[2] + freq3*solution.or1model_values[3] + freq4*solution.or1model_values[4] + freq5*solution.or1model_values[5] + freq6*solution.or1model_values[6];
      *or = derotate_180_small_double(*or);
    }else {
      *or = solution.or2model_values[0] + freq*solution.or2model_values[1] + freq2*solution.or2model_values[2] + freq3*solution.or2model_values[3] + freq4*solution.or2model_values[4] + freq5*solution.or2model_values[5] + freq6*solution.or2model_values[6];
      *or = derotate_180_small_double(*or);
    }
    return 1;
  }else if((polnr == 1 && solution.or1model == 8) || (polnr == 2 && solution.or2model == 8)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    double freq6 = freq3*freq3;
    double freq7 = freq4*freq3;
    if(polnr == 1) {
      *or = solution.or1model_values[0] + freq*solution.or1model_values[1] + freq2*solution.or1model_values[2] + freq3*solution.or1model_values[3] + freq4*solution.or1model_values[4] + freq5*solution.or1model_values[5] + freq6*solution.or1model_values[6] + freq7*solution.or1model_values[7];
      *or = derotate_180_small_double(*or);
    }else {
      *or = solution.or2model_values[0] + freq*solution.or2model_values[1] + freq2*solution.or2model_values[2] + freq3*solution.or2model_values[3] + freq4*solution.or2model_values[4] + freq5*solution.or2model_values[5] + freq6*solution.or2model_values[6] + freq7*solution.or2model_values[7];
      *or = derotate_180_small_double(*or);
    }
    return 1;
  }else if((polnr == 1 && solution.or1model == 9) || (polnr == 2 && solution.or2model == 9)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    double freq6 = freq3*freq3;
    double freq7 = freq4*freq3;
    double freq8 = freq4*freq4;
    if(polnr == 1) {
      *or = solution.or1model_values[0] + freq*solution.or1model_values[1] + freq2*solution.or1model_values[2] + freq3*solution.or1model_values[3] + freq4*solution.or1model_values[4] + freq5*solution.or1model_values[5] + freq6*solution.or1model_values[6] + freq7*solution.or1model_values[7] + freq8*solution.or1model_values[8];
      *or = derotate_180_small_double(*or);
    }else {
      *or = solution.or2model_values[0] + freq*solution.or2model_values[1] + freq2*solution.or2model_values[2] + freq3*solution.or2model_values[3] + freq4*solution.or2model_values[4] + freq5*solution.or2model_values[5] + freq6*solution.or2model_values[6] + freq7*solution.or2model_values[7] + freq8*solution.or2model_values[8];
      *or = derotate_180_small_double(*or);
    }
    return 1;
  }else if((polnr == 1 && solution.or1model == 10) || (polnr == 2 && solution.or2model == 10)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    double freq6 = freq3*freq3;
    double freq7 = freq4*freq3;
    double freq8 = freq4*freq4;
    double freq9 = freq5*freq4;
    if(polnr == 1) {
      *or = solution.or1model_values[0] + freq*solution.or1model_values[1] + freq2*solution.or1model_values[2] + freq3*solution.or1model_values[3] + freq4*solution.or1model_values[4] + freq5*solution.or1model_values[5] + freq6*solution.or1model_values[6] + freq7*solution.or1model_values[7] + freq8*solution.or1model_values[8] + freq9*solution.or1model_values[9];
      *or = derotate_180_small_double(*or);
    }else {
      *or = solution.or2model_values[0] + freq*solution.or2model_values[1] + freq2*solution.or2model_values[2] + freq3*solution.or2model_values[3] + freq4*solution.or2model_values[4] + freq5*solution.or2model_values[5] + freq6*solution.or2model_values[6] + freq7*solution.or2model_values[7] + freq8*solution.or2model_values[8] + freq9*solution.or2model_values[9];
      *or = derotate_180_small_double(*or);
    }
    return 1;
  }else if((polnr == 1 && solution.or1model == 11) || (polnr == 2 && solution.or2model == 11)) {
    double freq2 = freq*freq;
    double freq3 = freq2*freq;
    double freq4 = freq2*freq2;
    double freq5 = freq3*freq2;
    double freq6 = freq3*freq3;
    double freq7 = freq4*freq3;
    double freq8 = freq4*freq4;
    double freq9 = freq5*freq4;
    double freq10 = freq5*freq5;
    if(polnr == 1) {
      *or = solution.or1model_values[0] + freq*solution.or1model_values[1] + freq2*solution.or1model_values[2] + freq3*solution.or1model_values[3] + freq4*solution.or1model_values[4] + freq5*solution.or1model_values[5] + freq6*solution.or1model_values[6] + freq7*solution.or1model_values[7] + freq8*solution.or1model_values[8] + freq9*solution.or1model_values[9] + freq10*solution.or1model_values[10];
      *or = derotate_180_small_double(*or);
    }else {
      *or = solution.or2model_values[0] + freq*solution.or2model_values[1] + freq2*solution.or2model_values[2] + freq3*solution.or2model_values[3] + freq4*solution.or2model_values[4] + freq5*solution.or2model_values[5] + freq6*solution.or2model_values[6] + freq7*solution.or2model_values[7] + freq8*solution.or2model_values[8] + freq9*solution.or2model_values[9] + freq10*solution.or2model_values[10];
      *or = derotate_180_small_double(*or);
    }
    return 1;
  }else {
    fflush(stdout);
    if(polnr == 1) {
      printerror(0, "ERROR analyticReceiverSolution_orientation: Orientation model %d is not implemented", solution.or1model);
    }else {
      printerror(0, "ERROR analyticReceiverSolution_orientation: Orientation model %d is not implemented", solution.or2model);
    }
    return 0;
  }
  return 1;
}

struct {
  int param;                                  // Parameter we want to fit for
  analyticReceiverSolution_def solution;
  int nrPoints;
  double *xdata;
  double *ydata;
  double *yerror;
  int yerror_defined;
}fitReceiverModel_funk_info;

double fitReceiverModel_funk(double x[])
{
  double y, dy, chi2, xmhz;
  int i, j, polnr;
  verbose_definition verbose;
  cleanVerboseState(&verbose);

  chi2 = 0;
  for(i = 0; i < fitReceiverModel_funk_info.nrPoints; i++) {
    dy = 0;
    xmhz = fitReceiverModel_funk_info.xdata[i]*1000.0+fitReceiverModel_funk_info.solution.reffreq;
    if(fitReceiverModel_funk_info.param == 0) {   // gain
      //      fitReceiverModel_funk_info.solution.gainmodel_values[0] = x[0];
      //      fitReceiverModel_funk_info.solution.gainmodel_values[1] = x[1];
      int jmax;
      jmax = fitReceiverModel_funk_info.solution.gainmodel;
      if(jmax == 1)
	jmax = 2;
      for(j = 0; j < jmax; j++) {
	fitReceiverModel_funk_info.solution.gainmodel_values[j] = x[j];
      }
      if(analyticReceiverSolution_gain(fitReceiverModel_funk_info.solution, xmhz, &y, verbose) == 0) {
	fflush(stdout);
	printerror(0, "ERROR fitReceiverModel_funk: Cannot obtain value from model");
	exit(0);
      }
      dy = fitReceiverModel_funk_info.ydata[i] - y;
      if(fitReceiverModel_funk_info.yerror_defined) {
	chi2 += dy*dy/(fitReceiverModel_funk_info.yerror[i]*fitReceiverModel_funk_info.yerror[i]);
      }else {
	chi2 += dy*dy;
      }
    }else if(fitReceiverModel_funk_info.param == 1) {   // Differential gain
      //      fitReceiverModel_funk_info.solution.diffgainmodel_values[0] = x[0];
      //      fitReceiverModel_funk_info.solution.diffgainmodel_values[1] = x[1];
      int jmax;
      jmax = fitReceiverModel_funk_info.solution.diffgainmodel;
      if(jmax == 1)
	jmax = 2;
      for(j = 0; j < jmax; j++) {
	fitReceiverModel_funk_info.solution.diffgainmodel_values[j] = x[j];
      }
      if(analyticReceiverSolution_diffgain(fitReceiverModel_funk_info.solution, xmhz, &y, verbose) == 0) {
	fflush(stdout);
	printerror(0, "ERROR fitReceiverModel_funk: Cannot obtain value from model");
	exit(0);
      }
      dy = fitReceiverModel_funk_info.ydata[i] - y;
      //      printerror(0, "XXXX entered funk: dy=%f errordefine=%d error=%f", dy, fitReceiverModel_funk_info.yerror_defined, fitReceiverModel_funk_info.yerror[i]);
      if(fitReceiverModel_funk_info.yerror_defined) {
	chi2 += dy*dy/(fitReceiverModel_funk_info.yerror[i]*fitReceiverModel_funk_info.yerror[i]);
      }else {
	chi2 += dy*dy;
      }
      //      printerror(0, "XXXX entered funk: chi2=%f", chi2);
    }else if(fitReceiverModel_funk_info.param == 2) {  // Differential phase
      //      fitReceiverModel_funk_info.solution.phasemodel_values[0] = x[0];
      //      fitReceiverModel_funk_info.solution.phasemodel_values[1] = x[1];
      int jmax;
      jmax = fitReceiverModel_funk_info.solution.phasemodel;
      if(jmax == 1)
	jmax = 2;
      for(j = 0; j < jmax; j++) {
	fitReceiverModel_funk_info.solution.phasemodel_values[j] = x[j];
      }
      if(analyticReceiverSolution_diffphase(fitReceiverModel_funk_info.solution, xmhz, &y, verbose) == 0) {
	fflush(stdout);
	printerror(0, "ERROR fitReceiverModel_funk: Cannot obtain value from model");
	exit(0);
      }
      dy = fitReceiverModel_funk_info.ydata[i] - y;
      dy = derotate_180_small_double(dy);   // Resolve wraps in difference between -90 and +90 deg
      //      printerror(0, "XXXX entered funk: dy=%f errordefine=%d error=%f", dy, fitReceiverModel_funk_info.yerror_defined, fitReceiverModel_funk_info.yerror[i]);
      if(fitReceiverModel_funk_info.yerror_defined) {
	chi2 += dy*dy/(fitReceiverModel_funk_info.yerror[i]*fitReceiverModel_funk_info.yerror[i]);
      }else {
	chi2 += dy*dy;
      }
    }else if(fitReceiverModel_funk_info.param == 3 || fitReceiverModel_funk_info.param == 5) {   // Ellipticity
      if(fitReceiverModel_funk_info.param == 3) {
	//	fitReceiverModel_funk_info.solution.ell1model_values[0] = x[0];
	//	fitReceiverModel_funk_info.solution.ell1model_values[1] = x[1];
	int jmax;
	jmax = fitReceiverModel_funk_info.solution.ell1model;
	if(jmax == 1)
	  jmax = 2;
	for(j = 0; j < jmax; j++) {
	  fitReceiverModel_funk_info.solution.ell1model_values[j] = x[j];
	}
	polnr = 1;
      }else {
	//	fitReceiverModel_funk_info.solution.ell2model_values[0] = x[0];
	//	fitReceiverModel_funk_info.solution.ell2model_values[1] = x[1];
	int jmax;
	jmax = fitReceiverModel_funk_info.solution.ell2model;
	if(jmax == 1)
	  jmax = 2;
	for(j = 0; j < jmax; j++) {
	  fitReceiverModel_funk_info.solution.ell2model_values[j] = x[j];
	}
	polnr = 2;
      }
      if(analyticReceiverSolution_ellipticity(fitReceiverModel_funk_info.solution, polnr, xmhz, &y) == 0) {
	fflush(stdout);
	printerror(0, "ERROR fitReceiverModel_funk: Cannot obtain value from model");
	exit(0);
      }
      dy = fitReceiverModel_funk_info.ydata[i] - y;
      dy = derotate_180_small_double(dy);   // Resolve wraps in difference between -90 and +90 deg
      if(fitReceiverModel_funk_info.yerror_defined) {
	chi2 += dy*dy/(fitReceiverModel_funk_info.yerror[i]*fitReceiverModel_funk_info.yerror[i]);
      }else {
	chi2 += dy*dy;
      }
    }else if(fitReceiverModel_funk_info.param == 4 || fitReceiverModel_funk_info.param == 6) {   // Orientation
      if(fitReceiverModel_funk_info.param == 4) {
	//	if(i == 0) {
	//	  fprintf(stderr, "Orientation 1: x[0]=%e  x[1]=%e\n", x[0], x[1]);
	//	}
	//	fitReceiverModel_funk_info.solution.or1model_values[0] = x[0];
	//	fitReceiverModel_funk_info.solution.or1model_values[1] = x[1];
	int jmax;
	jmax = fitReceiverModel_funk_info.solution.or1model;
	if(jmax == 1)
	  jmax = 2;
	for(j = 0; j < jmax; j++) {
	  fitReceiverModel_funk_info.solution.or1model_values[j] = x[j];
	}
	polnr = 1;
      }else {
	//	fitReceiverModel_funk_info.solution.or2model_values[0] = x[0];
	//	fitReceiverModel_funk_info.solution.or2model_values[1] = x[1];
	int jmax;
	jmax = fitReceiverModel_funk_info.solution.or2model;
	if(jmax == 1)
	  jmax = 2;
	for(j = 0; j < jmax; j++) {
	  fitReceiverModel_funk_info.solution.or2model_values[j] = x[j];
	}
	polnr = 2;
      }
      if(analyticReceiverSolution_orientation(fitReceiverModel_funk_info.solution, polnr, xmhz, &y) == 0) {
	fflush(stdout);
	printerror(0, "ERROR fitReceiverModel_funk: Cannot obtain value from model");
	exit(0);
      }
      dy = fitReceiverModel_funk_info.ydata[i] - y;
      dy = derotate_180_small_double(dy);   // Resolve wraps in difference between -90 and +90 deg
      if(fitReceiverModel_funk_info.yerror_defined) {
	chi2 += dy*dy/(fitReceiverModel_funk_info.yerror[i]*fitReceiverModel_funk_info.yerror[i]);
      }else {
	chi2 += dy*dy;
      }
    }
  }
  //  printerror(0, "      chi2 = %f", chi2);
  return chi2+1e-5;  // Make the chi2 a bit larger. If the datapoints are zero, and the chi2=0, the algorithm can get stuck, possibly only if the initial guess is perfect (not tested).
}

/* Fits a analytic solution to a receiver model. phasemodel defines
   the type of model to fit for the differential phase solution,
   diffgainmodel is for the differential gain, gainmodel, ell1model
   and ell2model is for the ellipticity and or1model.  0 = not fitted,
   otherwise it is a n-1 order polynomial, so 1 = constant, 2 = linear
   function etc.

  In addition to the fit models, a reference models can be specified,
  i.e. a specific possibly complicated frequency dependent shape. See
  also the function refmodellist()

  If rmsOutliers is positive, all points rmsOutliers times the rms
  will be thrown away while deriving the solution. 

  Referencemodel models can be specified, which will be subtracted
  from the solution before fitting is done. If a reference model is
  set to -1, all reference models are explored.

  Verbose determines the number of spaces in front of the output.
  Return 1 = success, 0 on error */
int fitReceiverModel(datafile_definition solution, int gainmodel, int diffgainmodel, int phasemodel, long phasereferencemodel, long diffgainreferencemodel, int ell1model, int ell2model, int or1model, int or2model, double rmsOutliers, analyticReceiverSolution_def *analyticsolution, verbose_definition verbose)
{
  double *xdata, *ydata, *yerrorbar, *ypred, xstart[MaxNrfitReceiverModelFitParameters], dx[MaxNrfitReceiverModelFitParameters], xfit[MaxNrfitReceiverModelFitParameters], chi2, chi2best, ftol, ydata_min, ydata_max, rms_solution, rms_point, currefmodeltrial_chi2best;
  int i, subint, pol, bin, freq, freq2, dofit, ret, fixed[MaxNrfitReceiverModelFitParameters], nrparams, nfunk, trial, trial2, throwOutOutliersItt;
  int currefmodeltrial;
  long referencemodel_best;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Fit receiver model\n");
  }

  currefmodeltrial_chi2best = 0;  // Not necessary, but to get rid of valgrind errors

  clear_analyticReceiverSolution(&(fitReceiverModel_funk_info.solution));
  clear_analyticReceiverSolution(analyticsolution);
  fitReceiverModel_funk_info.solution.reffreq = get_centre_frequency(solution, verbose);
  analyticsolution->reffreq = get_centre_frequency(solution, verbose);
  for(i = 0; i < MaxNrfitReceiverModelFitParameters; i++) {
    fixed[i] = 0;
  }
  ftol = 1e-4;
  analyticsolution->mjd = solution.mjd_start;

  if(solution.gentype != GENTYPE_RECEIVERMODEL2 && solution.gentype != GENTYPE_RECEIVERMODEL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fitReceiverModel: File does not appear to be a receiver solution.");
    return 0;
  }
  if(solution.NrPols < 3) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fitReceiverModel: Expecteded at least 3 receiver parameters, only %ld present.", solution.NrPols);
    return 0;
  }
  if(solution.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fitReceiverModel: receiver model should be read into memory first");
    return 0;
  }

  xdata = (double *)malloc(solution.NrFreqChan*sizeof(double));
  ydata = (double *)malloc(solution.NrFreqChan*sizeof(double));
  ypred = (double *)malloc(solution.NrFreqChan*sizeof(double));
  yerrorbar = (double *)malloc(solution.NrFreqChan*sizeof(double));
  if(xdata == NULL || ydata == NULL || yerrorbar == NULL || ypred == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fitReceiverModel: Memory allocation error.");
    return 0;
  }
  fitReceiverModel_funk_info.xdata = xdata;
  fitReceiverModel_funk_info.ydata = ydata;
  fitReceiverModel_funk_info.yerror = yerrorbar;

  /* Loop over pulses and frequency channels etc. */
  for(subint = 0; subint < solution.NrSubints; subint++) {
    for(pol = 0; pol < solution.NrPols; pol++) {
      referencemodel_best = 0;
      // Loop over all possible reference models. Only actually used when no refmodel is specified. 0 = no model subtracted.
      for(currefmodeltrial = 0; currefmodeltrial <= 1000; currefmodeltrial++) {
	if(pol != 1 && pol != 2) // No subtraction of reference model implemented for pols other than diff gain and phase
	  currefmodeltrial = 1000;
	for(throwOutOutliersItt = 0; throwOutOutliersItt < 4; throwOutOutliersItt++) {
	  if(rmsOutliers <= 0 && throwOutOutliersItt == 1)   // If not interested in throwing away outliers, quit loop.
	    break;
	  if(throwOutOutliersItt >= 1) {  // Determine the rms of the solution
	    for(freq = 0; freq < fitReceiverModel_funk_info.nrPoints; freq++) {
	      if(pol == 0) {
		if(analyticReceiverSolution_gain(*analyticsolution, xdata[freq], &(ypred[freq]), verbose) == 0) {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR fitReceiverModel_funk: Cannot obtain value from model");
		  return 0;
		}
	      }else if(pol == 1) {
		if(analyticReceiverSolution_diffgain(*analyticsolution, xdata[freq], &(ypred[freq]), verbose) == 0) {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR fitReceiverModel_funk: Cannot obtain value from model");
		  return 0;
		}
	      }else if(pol == 2) {
		if(analyticReceiverSolution_diffphase(*analyticsolution, xdata[freq], &(ypred[freq]), verbose) == 0) {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR fitReceiverModel_funk: Cannot obtain value from model");
		  return 0;
		}
	      }else if(pol == 3) {
		if(analyticReceiverSolution_ellipticity(*analyticsolution, 1, xdata[freq], &(ypred[freq])) == 0) {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR fitReceiverModel_funk: Cannot obtain value from model");
		  return 0;
		}
	      }else if(pol == 5) {
		if(analyticReceiverSolution_ellipticity(*analyticsolution, 2, xdata[freq], &(ypred[freq])) == 0) {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR fitReceiverModel_funk: Cannot obtain value from model");
		  return 0;
		}
	      }else if(pol == 4) {
		if(analyticReceiverSolution_orientation(*analyticsolution, 1, xdata[freq], &(ypred[freq])) == 0) {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR fitReceiverModel_funk: Cannot obtain value from model");
		  return 0;
		}
	      }else if(pol == 6) {
		if(analyticReceiverSolution_orientation(*analyticsolution, 2, xdata[freq], &(ypred[freq])) == 0) {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR fitReceiverModel_funk: Cannot obtain value from model");
		  return 0;
		}
	      }
	      if(pol == 0 || pol == 1) {  // Not an angle
		rms_solution += (ydata[freq]-ypred[freq])*(ydata[freq]-ypred[freq]);
	      }else {
		rms_solution += derotate_180_small_double(ydata[freq]-ypred[freq])*derotate_180_small_double(ydata[freq]-ypred[freq]);
	      }
	    }
	    rms_solution /= (double)fitReceiverModel_funk_info.nrPoints;
	    rms_solution = sqrt(rms_solution);
	  } // End of throwOutOutliersItt if statement
	  if(solution.NrBins == 2) {           // By default, errors are used if they are available
	    fitReceiverModel_funk_info.yerror_defined = 1;
	  }else {
	    fitReceiverModel_funk_info.yerror_defined = 0;
	  }
	  fitReceiverModel_funk_info.param = pol;  // The current parameter which might be fitted for below

	  // Fill arrays with receiver parameters
	  for(bin = 0; bin < solution.NrBins; bin++) {
	    if(bin == 0) {
	      fitReceiverModel_funk_info.nrPoints = 0;
	      for(freq = 0; freq < solution.NrFreqChan; freq++) {
		/* If gain is zero, skip data point. */
		if(solution.data[solution.NrBins*(0+solution.NrPols*(freq+subint*solution.NrFreqChan))+0] > 0) {
		  ydata[fitReceiverModel_funk_info.nrPoints] = solution.data[solution.NrBins*(pol+solution.NrPols*(freq+subint*solution.NrFreqChan))+0];
		  //		  xdata[fitReceiverModel_funk_info.nrPoints] = solution.freq_cent-0.5*(solution.NrFreqChan-1.0)*get_channelbw(solution) + freq*get_channelbw(solution);
		  //		  xdata[fitReceiverModel_funk_info.nrPoints] = get_nonweighted_channel_freq(solution, freq, verbose);
		  xdata[fitReceiverModel_funk_info.nrPoints] = 0.001*(get_weighted_channel_freq(solution, subint, freq, verbose) - analyticsolution->reffreq);
		  fitReceiverModel_funk_info.nrPoints++;
		}
	      }
	    }else if(bin == 1) {
	      fitReceiverModel_funk_info.nrPoints = 0;
	      for(freq = 0; freq < solution.NrFreqChan; freq++) {
		/* If gain is zero, skip data point. */
		if(solution.data[solution.NrBins*(0+solution.NrPols*(freq+subint*solution.NrFreqChan))+0] > 0) {
		  yerrorbar[fitReceiverModel_funk_info.nrPoints] = solution.data[solution.NrBins*(pol+solution.NrPols*(freq+subint*solution.NrFreqChan))+1];
		  fitReceiverModel_funk_info.nrPoints++;
		}
	      }
	    }
	  }
	  freq2 = 0;
	  for(freq = 0; freq < fitReceiverModel_funk_info.nrPoints; freq++) {
	    if(isnan(ydata[freq]) || isinf(ydata[freq])) {
	      ydata[freq] = 0;
	      yerrorbar[freq] = 0;
	    }else {
	      if(solution.gentype == GENTYPE_RECEIVERMODEL2 && pol == solution.NrPols-2) {  //chi^2
		yerrorbar[freq] = 0;
	      }else if(solution.gentype == GENTYPE_RECEIVERMODEL2 && pol == solution.NrPols-1) {  // Nfree
		yerrorbar[freq] = 0;
	      }else if(pol == 1) {  // Differential gain
		ydata[freq] = 100.0*(exp(2.0*ydata[freq])-1.0);
		//	  yerrorbar[freq] *= 200.0*exp(2.0*ydata[freq]); // Assume error-bar is very small 
		// didn't work properly, as lower points get far too much weight. Simply set to uniform errorbars
		yerrorbar[freq] = 1;
	      }else if(pol >= 2 && pol <= 6) {   // Angles
		ydata[freq] *= 180.0/M_PI;
		yerrorbar[freq] *= 180.0/M_PI;
	      }
	      if(pol == 0 || pol == 1) {  // Not an angle
		rms_point = fabs(ydata[freq] - ypred[freq]);
	      }else {
		//	      printerror(debug, "freq %d: %f %f", freq, ydata[freq], ypred[freq]);
		if(throwOutOutliersItt == 0)
		  rms_point = 0;                   // Shouldn't be used in the first itteration
		else
		  rms_point = fabs(derotate_180_small_double(ydata[freq] - ypred[freq]));
	      }
	      if(throwOutOutliersItt == 0 || rms_point < rms_solution *rmsOutliers) {
		if(freq == 0 || ydata[freq] < ydata_min) {
		  ydata_min = ydata[freq];
		}
		if(freq == 0 || ydata[freq] > ydata_max) {
		  ydata_max = ydata[freq];
		}
		xdata[freq2] = xdata[freq];
		ydata[freq2] = ydata[freq];
		yerrorbar[freq2] = yerrorbar[freq];
		freq2++;
	      }else {
		if(verbose.verbose) {
		  for(i = 0; i < verbose.indent; i++)
		    printf(" ");
		  printf("Itteration %d: Throw away freq %f MHz\n", throwOutOutliersItt, 1000*xdata[freq]+analyticsolution->reffreq);
		}
	      }
	    }
	  }
	  fitReceiverModel_funk_info.nrPoints = freq2;
	
	  chi2best = 0;
	  // Each trial is a new start point of the fit (if defined)
	  // 90 = 9*10 trials if 2D search
	  for(trial = 0; trial < 90; trial++) {   
	    dofit = 0;
	    if(pol == 0) { //Gain
	      if(gainmodel > 0) {   // If fitting of gain is requested
		fitReceiverModel_funk_info.solution.gainmodel = gainmodel;
		analyticsolution->gainmodel                   = gainmodel;
		fitReceiverModel_funk_info.solution.gainmodel_defined = 1;
		analyticsolution->gainmodel_defined                   = 1;
		if(gainmodel >= 1 || gainmodel <= 11) {
		  if(gainmodel == 1)
		    nrparams = 2;  // Should be 1, but then amouba doesn't work. Setting it to two still finds solution
		  else
		    nrparams = gainmodel; 
		  dx[0] = 0.1*(ydata_max-ydata_min);     // 10% variation
		  // y = a*x^n
		  // dy = a*(xmax^n-xcent^n)   but xcent = 0 since freq = 0.001*(freq-solution.reffreq);
		  // dy = dx*xmax^n    
		  double dy = 0.01*(ydata_max-ydata_min);    // dy = max extension param, take to be 1%.
		  double xmax = 0.5*get_bandwidth(solution, verbose)*0.001;
		  double xmax2 = xmax*xmax;
		  double xmax3 = xmax2*xmax;
		  double xmax4 = xmax2*xmax2;
		  double xmax5 = xmax3*xmax2;
		  dx[1] = 10*dy/xmax;    // 10% per 100 MHz
		  dx[2] = dy/(xmax2);   // 1% per 100 MHz
		  dx[3] = dy/(xmax3);   // 1% per 100 MHz
		  dx[4] = dy/(xmax4);   // 1% per 100 MHz
		  dx[5] = dy/(xmax5);   // 1% per 100 MHz
		  dx[6] = dy/(xmax5*xmax);   // 1% per 100 MHz
		  dx[7] = dy/(xmax5*xmax2);   // 1% per 100 MHz
		  dx[8] = dy/(xmax5*xmax3);   // 1% per 100 MHz
		  dx[9] = dy/(xmax5*xmax4);   // 1% per 100 MHz
		  dx[10] = dy/(xmax5*xmax5);   // 1% per 100 MHz
		  //		  dx[11] = dy/(xmax5*xmax5*xmax);   // 1% per 100 MHz
		  //		  dx[1] = 0.1*(ydata_max-ydata_min)/100.0;   // 10% per 100 MHz
		
		  xstart[0] = (-0.5*(10-1) + trial % 10) * 0.15*(ydata_max-ydata_min);   // Offset - 10 trials
		  trial2 = trial/10;
		  xstart[1] = dx[1]*(trial2-0.5*(10-2));   // Gradient - 9 trials
		  xstart[2] = 0;
		  xstart[3] = 0;
		  xstart[4] = 0;
		  xstart[5] = 0;
		  xstart[6] = 0;
		  xstart[7] = 0;
		  xstart[8] = 0;
		  xstart[9] = 0;
		  xstart[10] = 0;
		  //		  xstart[11] = 0;
		  
		  //		printf("xstart = %f %f\n", xstart[0], xstart[1]);
		  
		  /*		  if(gainmodel == 2) {
		    //		    xstart[0] -= solution.freq_cent*xstart[1]; // Make sure that the defined offset is that at the centre frequency rather then zero frequency
		    xstart[0] -= get_centre_frequency(solution, verbose)*xstart[1]; // Make sure that the defined offset is that at the centre frequency rather then zero frequency
		    }*/
		  dofit = 1;
		}else {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR fitReceiverModel: gain model %d not implemented", gainmodel);
		  return 0;
		}
	      }
	    }else if(solution.gentype == GENTYPE_RECEIVERMODEL2 && pol == solution.NrPols-2) {  //chi^2
	    }else if(solution.gentype == GENTYPE_RECEIVERMODEL2 && pol == solution.NrPols-1) {  // Nfree
	    }else if(pol == 1) {  // Differential gain
	      if(diffgainmodel <= 0) {   // If no fitting of differential gain is requested
		currefmodeltrial = 20000;
	      }else {
		if(diffgainreferencemodel != -1) { // Use specified model
		  fitReceiverModel_funk_info.solution.diffgainRefmodel = diffgainreferencemodel;
		  analyticsolution->diffgainRefmodel                   = diffgainreferencemodel;
		    currefmodeltrial = 20000;  // We're done after re-calculating best reference model. Make sure this is last time we go through loop.
		}else {                 // Trying out all available models
		  if(currefmodeltrial == 0) {  // First trial: no subtraction
		    fitReceiverModel_funk_info.solution.diffgainRefmodel = 0;
		    analyticsolution->diffgainRefmodel                   = 0;
		  }else {
		    if(diffgainmodels_names(currefmodeltrial) != -1) {    // A unexplored model exists
		      fitReceiverModel_funk_info.solution.diffgainRefmodel = diffgainmodels_names(currefmodeltrial);
		      analyticsolution->diffgainRefmodel                   = diffgainmodels_names(currefmodeltrial);
		    }else {
		      if(referencemodel_best == 0) {
			fitReceiverModel_funk_info.solution.diffgainRefmodel = 0;
			analyticsolution->diffgainRefmodel                   = 0;
		      }else {
			fitReceiverModel_funk_info.solution.diffgainRefmodel = diffgainmodels_names(referencemodel_best);
			analyticsolution->diffgainRefmodel                   = diffgainmodels_names(referencemodel_best);
		      }
		      currefmodeltrial = 20000;  // We're done after re-calculating best reference model. Make sure this is last time we go through loop.
		    }
		  }
		}
		fitReceiverModel_funk_info.solution.diffgainmodel = diffgainmodel;
		analyticsolution->diffgainmodel                   = diffgainmodel;
		fitReceiverModel_funk_info.solution.diffgainmodel_defined = 1;
		analyticsolution->diffgainmodel_defined                   = 1;
		if(diffgainmodel >= 1 && diffgainmodel <= 11) {
		  if(diffgainmodel == 1) {
		    nrparams = 2;  // Should be 1, but then amouba doesn't work. Setting it to two still finds solution
		  }else {
		    nrparams = diffgainmodel;
		  }
		  dx[0] = 10;        // 10% stepsize in diff gain as a constant term
		  // y = a*x^n
		  // dy = a*(xmax^n-xcent^n)   but xcent = 0 since freq = 0.001*(freq-solution.reffreq);
		  // dy = dx*xmax^n    
		  double dy = 1.0;    // dy = max extension param, take to be 1%.
		  double xmax = 0.5*get_bandwidth(solution, verbose)*0.001;
		  double xmax2 = xmax*xmax;
		  double xmax3 = xmax2*xmax;
		  double xmax4 = xmax2*xmax2;
		  double xmax5 = xmax3*xmax2;
		  dx[1] = 10/xmax;    // 10% per 100 MHz
		  dx[2] = dy/(xmax2);   // 1% per 100 MHz
		  dx[3] = dy/(xmax3);   // 1% per 100 MHz
		  dx[4] = dy/(xmax4);   // 1% per 100 MHz
		  dx[5] = dy/(xmax5);   // 1% per 100 MHz
		  dx[6] = dy/(xmax5*xmax);   // 1% per 100 MHz
		  dx[7] = dy/(xmax5*xmax2);   // 1% per 100 MHz
		  dx[8] = dy/(xmax5*xmax3);   // 1% per 100 MHz
		  dx[9] = dy/(xmax5*xmax4);   // 1% per 100 MHz
		  dx[10] = dy/(xmax5*xmax5);   // 1% per 100 MHz
		  //		  dx[11] = dy/(xmax5*xmax5*xmax);   // 1% per 100 MHz
		  
		  xstart[0] = (-0.5*(10-1) + trial % 10) * 100/(float)(10);   // Offset - 10 trials
		  trial2 = trial/10;
		  xstart[1] = dx[1]*(trial2-0.5*(10-2))*7.0/(float)10;   // Gradient - 9 trials
		  xstart[2] = 0;
		  xstart[3] = 0;
		  xstart[4] = 0;
		  xstart[5] = 0;
		  xstart[6] = 0;
		  xstart[7] = 0;
		  xstart[8] = 0;
		  xstart[9] = 0;
		  xstart[10] = 0;
		  //		  xstart[11] = 0;
		  //		  xstart[1] = 1.03;
		  //		  xstart[2] = -3.5e-4;

		  
		  /*		  if(diffgainmodel >= 2) {
		    //		    xstart[0] -= solution.freq_cent*xstart[1]; // Make sure that the defined offset is that at the centre frequency rather then zero frequency
		    xstart[0] -= get_centre_frequency(solution, verbose)*xstart[1]; // Make sure that the defined offset is that at the centre frequency rather then zero frequency
		    }*/
		  //		  printf("xstart = %e + %e * 0.001*(x-%e)\n", xstart[0], xstart[1], analyticsolution->reffreq);
		  //		  printf("dx = %e,  %e\n", dx[0], dx[1]);
		  dofit = 1;
		}else {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR fitReceiverModel: differential gain model %d not implemented", diffgainmodel);
		  return 0;
		}
	      }
	    }else if(pol == 2) {  // Differential phase
	      if(phasemodel <= 0) {   // If no fitting of phase is requested
		currefmodeltrial = 20000;
	      }else {
		if(phasereferencemodel != -1) {
		  fitReceiverModel_funk_info.solution.phaseRefmodel = phasereferencemodel;
		  analyticsolution->phaseRefmodel                   = phasereferencemodel;
		  currefmodeltrial = 20000;  // We're done after re-calculating best reference model. Make sure this is last time we go through loop.
		}else {            // Trying out all available models
		  if(currefmodeltrial == 0) {  // First trial: no subtraction
		    fitReceiverModel_funk_info.solution.phaseRefmodel = 0;
		    analyticsolution->phaseRefmodel                   = 0;
		  }else {
		    if(phasemodels_names(currefmodeltrial) != -1) {
		      fitReceiverModel_funk_info.solution.phaseRefmodel = phasemodels_names(currefmodeltrial);
		      analyticsolution->phaseRefmodel                   = phasemodels_names(currefmodeltrial);
		    }else {
		      if(referencemodel_best == 0) {
			fitReceiverModel_funk_info.solution.phaseRefmodel = 0;
			analyticsolution->phaseRefmodel                   = 0;
		      }else {
			fitReceiverModel_funk_info.solution.phaseRefmodel = phasemodels_names(referencemodel_best);
			analyticsolution->phaseRefmodel                   = phasemodels_names(referencemodel_best);
		      }
		      currefmodeltrial = 20000;  // We're done after re-calculating best reference model. Make sure this is last time we go through loop.
		    }
		  }
		}
		fitReceiverModel_funk_info.solution.phasemodel = phasemodel;
		analyticsolution->phasemodel                   = phasemodel;
		fitReceiverModel_funk_info.solution.phasemodel_defined = 1;
		analyticsolution->phasemodel_defined                   = 1;
		if(phasemodel >= 1 && phasemodel <= 11) {
		  if(phasemodel == 1) {
		    nrparams = 2; // Should be 1, but then amouba doesn't work. Setting it to two still finds solution
		  }else {
		    nrparams = phasemodel;
		  }
		  dx[0] = 45;           // 45 deg in phase as a constant term
		  // y = a*x^n
		  // dy = a*(xmax^n-xcent^n)   but xcent = 0 since freq = 0.001*(freq-solution.reffreq);
		  // dy = dx*xmax^n    
		  double dy = 10.0;    // dy = max extension param, take to be 10 deg
		  double xmax = 0.5*get_bandwidth(solution, verbose)*0.001;
		  double xmax2 = xmax*xmax;
		  double xmax3 = xmax2*xmax;
		  double xmax4 = xmax2*xmax2;
		  double xmax5 = xmax3*xmax2;
		  dx[1] = 90.0/xmax;   // 90 deg over BW
		  dx[2] = dy/(xmax2);   // 10 deg over BW
		  dx[3] = dy/(xmax3);   // 10 deg over BW
		  dx[4] = dy/(xmax4);   // 10 deg over BW
		  dx[5] = dy/(xmax5);   // 10 deg over BW
		  dx[6] = dy/(xmax5*xmax);   // 10 deg over BW
		  dx[7] = dy/(xmax5*xmax2);   // 10 deg over BW
		  dx[8] = dy/(xmax5*xmax3);   // 10 deg over BW
		  dx[9] = dy/(xmax5*xmax4);   // 10 deg over BW
		  dx[10] = dy/(xmax5*xmax5);   // 10 deg over BW
		  //		  dx[11] = dy/(xmax5*xmax5*xmax);   // 10 deg over BW
		  
		  xstart[0] = (trial % 10) * 180.0/(float)(10) - 90;   // Offset - 10 trials
		  trial2 = trial/10;
		  xstart[1] = dx[1]*(trial2-4);   // Gradient - 9 trials
		  xstart[2] = 0;
		  xstart[3] = 0;
		  xstart[4] = 0;
		  xstart[5] = 0;
		  xstart[6] = 0;
		  xstart[7] = 0;
		  xstart[8] = 0;
		  xstart[9] = 0;
		  xstart[10] = 0;
		  //		  xstart[11] = 0;
		  
		  //		printf("xstart = %f %f\n", xstart[0], xstart[1]);
		  
		  dofit = 1;
		}else {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR fitReceiverModel: phase model %d not implemented", phasemodel);
		  return 0;
		}
	      }
	    }else if(pol == 3 || pol == 5) {  // el0 or el1
	      if((pol == 3 && ell1model > 0) || (pol == 5 && ell2model > 0)) {   // If fitting of ellipticity is requested
		if(pol == 3) {
		  fitReceiverModel_funk_info.solution.ell1model = ell1model;
		  analyticsolution->ell1model                   = ell1model;
		  fitReceiverModel_funk_info.solution.ell1model_defined = 1;
		  analyticsolution->ell1model_defined                   = 1;
		}else {
		  fitReceiverModel_funk_info.solution.ell2model = ell2model;
		  analyticsolution->ell2model                   = ell2model;
		  fitReceiverModel_funk_info.solution.ell2model_defined = 1;
		  analyticsolution->ell2model_defined                   = 1;
		}
		if((pol == 3 && (ell1model >= 1 && ell1model <= 11)) || (pol == 5 && (ell2model >= 1 && ell2model <= 11))) {
		  if((pol == 3 && ell1model == 1) || (pol == 5 && ell2model == 1)) {
		    nrparams = 2; // Should be 1, but then amouba doesn't work. Setting it to two still finds solution
		  }else {
		    if(pol == 3)
		      nrparams = ell1model;
		    else
		      nrparams = ell2model;
		  }
		  dx[0] = 45;
		  // y = a*x^n
		  // dy = a*(xmax^n-xcent^n)   but xcent = 0 since freq = 0.001*(freq-solution.reffreq);
		  // dy = dx*xmax^n    
		  double dy = 1.0;    // dy = max extension param, take to be 1 deg
		  double xmax = 0.5*get_bandwidth(solution, verbose)*0.001;
		  double xmax2 = xmax*xmax;
		  double xmax3 = xmax2*xmax;
		  double xmax4 = xmax2*xmax2;
		  double xmax5 = xmax3*xmax2;
		  dx[1] = 90.0/xmax;   // 90 deg over BW
		  dx[2] = dy/(xmax2);   // 1 deg over BW
		  dx[3] = dy/(xmax3);   // 1 deg over BW
		  dx[4] = dy/(xmax4);   // 1 deg over BW
		  dx[5] = dy/(xmax5);   // 1 deg over BW
		  dx[6] = dy/(xmax5*xmax);   // 1 deg over BW
		  dx[7] = dy/(xmax5*xmax2);   // 1 deg over BW
		  dx[8] = dy/(xmax5*xmax3);   // 1 deg over BW
		  dx[9] = dy/(xmax5*xmax4);   // 1 deg over BW
		  dx[10] = dy/(xmax5*xmax5);   // 1 deg over BW
		  
		  xstart[0] = (trial % 10) * 180.0/(float)(10) - 90;   // Offset - 10 trials
		  
		  trial2 = trial/10;
		  xstart[1] = dx[1]*(trial2-0.5*(10-2))*0.5/(float)10;   // Gradient - 9 trials
		  xstart[2] = 0;
		  xstart[3] = 0;
		  xstart[4] = 0;
		  xstart[5] = 0;
		  xstart[6] = 0;
		  xstart[7] = 0;
		  xstart[8] = 0;
		  xstart[9] = 0;
		  xstart[10] = 0;
		  
		  //		printf("xstart = %f %f\n", xstart[0], xstart[1]);
		  
		  dofit = 1;
		}else {
		  fflush(stdout);
		  if(pol == 3) {
		    printerror(verbose.debug, "ERROR fitReceiverModel: ellipticity model %d not implemented", ell1model);
		  }else {
		    printerror(verbose.debug, "ERROR fitReceiverModel: ellipticity model %d not implemented", ell2model);
		  }
		  return 0;
		}
	      }
	    }else if(pol == 4 || pol == 6) {  // or0 or or1
	      if((pol == 4 && or1model > 0) || (pol == 6 && or2model > 0)) {   // If fitting of orientation is requested
		if(pol == 4) {
		  fitReceiverModel_funk_info.solution.or1model = or1model;
		  analyticsolution->or1model                   = or1model;
		  fitReceiverModel_funk_info.solution.or1model_defined = 1;
		  analyticsolution->or1model_defined                   = 1;
		}else {
		  fitReceiverModel_funk_info.solution.or2model = or2model;
		  analyticsolution->or2model                   = or2model;
		  fitReceiverModel_funk_info.solution.or2model_defined = 1;
		  analyticsolution->or2model_defined                   = 1;
		}
		if((pol == 4 && (or1model >= 1 && or1model <= 11)) || (pol == 6 && (or2model >= 1 || or2model <= 11))) {
		  if((pol == 4 && or1model == 1) || (pol == 6 && or2model == 1)) {
		    nrparams = 2; // Should be 1, but then amouba doesn't work. Setting it to two still finds solution
		  }else {
		    if(pol == 4)
		      nrparams = or1model;
		    else
		      nrparams = or2model;
		  }
		  dx[0] = 45;
		  // y = a*x^n
		  // dy = a*(xmax^n-xcent^n)   but xcent = 0 since freq = 0.001*(freq-solution.reffreq);
		  // dy = dx*xmax^n    
		  double dy = 1.0;    // dy = max extension param, take to be 1 deg
		  double xmax = 0.5*get_bandwidth(solution, verbose)*0.001;
		  double xmax2 = xmax*xmax;
		  double xmax3 = xmax2*xmax;
		  double xmax4 = xmax2*xmax2;
		  double xmax5 = xmax3*xmax2;
		  dx[1] = 90.0/xmax;   // 90 deg over BW
		  dx[2] = dy/(xmax2);   // 1 deg over BW
		  dx[3] = dy/(xmax3);   // 1 deg over BW
		  dx[4] = dy/(xmax4);   // 1 deg over BW
		  dx[5] = dy/(xmax5);   // 1 deg over BW
		  dx[6] = dy/(xmax5*xmax);   // 1 deg over BW
		  dx[7] = dy/(xmax5*xmax2);   // 1 deg over BW
		  dx[8] = dy/(xmax5*xmax3);   // 1 deg over BW
		  dx[9] = dy/(xmax5*xmax4);   // 1 deg over BW
		  dx[10] = dy/(xmax5*xmax5);   // 1 deg over BW
		  
		  xstart[0] = (trial % 10) * 180.0/(float)(10) - 90;   // Offset - 10 trials
		  trial2 = trial/10;
		  xstart[1] = dx[1]*(trial2-0.5*(10-2))*0.5/(float)10;   // Gradient - 9 trials
		  xstart[2] = 0;
		  xstart[3] = 0;
		  xstart[4] = 0;
		  xstart[5] = 0;
		  xstart[6] = 0;
		  xstart[7] = 0;
		  xstart[8] = 0;
		  xstart[9] = 0;
		  xstart[10] = 0;

		  //		printf("xstart = %f %f\n", xstart[0], xstart[1]);
		  dofit = 1;
		}else {
		  fflush(stdout);
		  if(pol == 4) {
		    printerror(verbose.debug, "ERROR fitReceiverModel: orientation model %d not implemented", or1model);
		  }else {
		    printerror(verbose.debug, "ERROR fitReceiverModel: orientation model %d not implemented", or2model);
		  }
		  return 0;
		}
	      }
	    }else {
	      fflush(stdout); 
	      printerror(verbose.debug, "ERROR fitReceiverModel: Something wrong.");
	    }
	    
	    if(dofit) {
	      ret = doAmoeba_d(0, xstart, dx, fixed, xfit, &chi2, nrparams, fitReceiverModel_funk, ftol, &nfunk, 1, 0, 3, NULL, NULL);
	      if(ret == 0) {
		if(pol == 0) { //Gain
		  if(trial == 0 || chi2 < chi2best) {
		    chi2best = chi2;	
		    int j, jmax;
		    jmax = fitReceiverModel_funk_info.solution.gainmodel;
		    if(jmax == 1)
		      jmax = 2;
		    for(j = 0; j < jmax; j++) {
		      analyticsolution->gainmodel_values[j] = xfit[j];
		    }
		    //		    analyticsolution->gainmodel_values[0] = xfit[0];
		    //		    analyticsolution->gainmodel_values[1] = xfit[1];
		    analyticsolution->gainmodel_redchi2 = chi2/(double)fitReceiverModel_funk_info.nrPoints;
		    if(currefmodeltrial == 0 || chi2best < currefmodeltrial_chi2best) {
		      currefmodeltrial_chi2best = chi2best;
		      referencemodel_best = currefmodeltrial;
		    }
		  }
		}else if(solution.gentype == GENTYPE_RECEIVERMODEL2 && pol == solution.NrPols-2) {  //chi^2
		}else if(solution.gentype == GENTYPE_RECEIVERMODEL2 && pol == solution.NrPols-1) {  // Nfree
		}else if(pol == 1) {  // Differential gain
		  if(trial == 0 || chi2 < chi2best) {
		    chi2best = chi2;
		    int j, jmax;
		    jmax = fitReceiverModel_funk_info.solution.diffgainmodel;
		    if(jmax == 1)
		      jmax = 2;
		    for(j = 0; j < jmax; j++) {
		      analyticsolution->diffgainmodel_values[j] = xfit[j];
		    }
		    //		    analyticsolution->diffgainmodel_values[0] = xfit[0];
		    //		    analyticsolution->diffgainmodel_values[1] = xfit[1];

		    analyticsolution->diffgainmodel_redchi2 = chi2/(double)fitReceiverModel_funk_info.nrPoints;
		    if(currefmodeltrial == 0 || chi2best < currefmodeltrial_chi2best) {
		      currefmodeltrial_chi2best = chi2best;
		      referencemodel_best = currefmodeltrial;
		    }
		  }
		}else if(pol == 2) {  // Differential phase
		  if(trial == 0 || chi2 < chi2best) {
		    chi2best = chi2;
		    analyticsolution->phasemodel_values[0] = derotate_180_small_double(xfit[0]);
		    if(fitReceiverModel_funk_info.solution.phasemodel > 1) {
		      int j, jmax;
		      jmax = fitReceiverModel_funk_info.solution.phasemodel;
		      for(j = 0; j < jmax; j++) {
			analyticsolution->phasemodel_values[j] = xfit[j];
		      }
		    }
		    //		    analyticsolution->phasemodel_values[1] = xfit[1];
		    analyticsolution->phasemodel_redchi2 = chi2/(double)fitReceiverModel_funk_info.nrPoints;
		    //		    printf("XXXXXX %d %lf %lf", currefmodeltrial, chi2best, currefmodeltrial_chi2best);
		    if(currefmodeltrial == 0 || chi2best < currefmodeltrial_chi2best) {
		      currefmodeltrial_chi2best = chi2best;
		      referencemodel_best = currefmodeltrial;
		    }
		  }
		}else if(pol == 3) {  // el0
		  if(trial == 0 || chi2 < chi2best) {
		    chi2best = chi2;
		    int j, jmax;
		    jmax = fitReceiverModel_funk_info.solution.ell1model;
		    if(jmax == 1)
		      jmax = 2;
		    for(j = 0; j < jmax; j++) {
		      analyticsolution->ell1model_values[j] = xfit[j];
		    }
		    //		    analyticsolution->ell1model_values[0] = xfit[0];
		    //		    analyticsolution->ell1model_values[1] = xfit[1];
		    analyticsolution->ell1model_redchi2 = chi2/(double)fitReceiverModel_funk_info.nrPoints;
		    if(currefmodeltrial == 0 || chi2best < currefmodeltrial_chi2best) {
		      currefmodeltrial_chi2best = chi2best;
		      referencemodel_best = currefmodeltrial;
		    }
		  }
		}else if(pol == 4) {  // or0
		  if(trial == 0 || chi2 < chi2best) {
		    chi2best = chi2;
		    int j, jmax;
		    jmax = fitReceiverModel_funk_info.solution.or1model;
		    if(jmax == 1)
		      jmax = 2;
		    for(j = 0; j < jmax; j++) {
		      analyticsolution->or1model_values[j] = xfit[j];
		    }
		    //		    analyticsolution->or1model_values[0] = xfit[0];
		    //		    analyticsolution->or1model_values[1] = xfit[1];
		    analyticsolution->or1model_redchi2 = chi2/(double)fitReceiverModel_funk_info.nrPoints;
		    if(currefmodeltrial == 0 || chi2best < currefmodeltrial_chi2best) {
		      currefmodeltrial_chi2best = chi2best;
		      referencemodel_best = currefmodeltrial;
		    }
		  }
		}else if(pol == 5) {  // el1
		  if(trial == 0 || chi2 < chi2best) {
		    chi2best = chi2;
		    int j, jmax;
		    jmax = fitReceiverModel_funk_info.solution.ell2model;
		    if(jmax == 1)
		      jmax = 2;
		    for(j = 0; j < jmax; j++) {
		      analyticsolution->ell2model_values[j] = xfit[j];
		    }
		    //		    analyticsolution->ell2model_values[0] = xfit[0];
		    //		    analyticsolution->ell2model_values[1] = xfit[1];
		    analyticsolution->ell2model_redchi2 = chi2/(double)fitReceiverModel_funk_info.nrPoints;
		    if(currefmodeltrial == 0 || chi2best < currefmodeltrial_chi2best) {
		      currefmodeltrial_chi2best = chi2best;
		      referencemodel_best = currefmodeltrial;
		    }
		  }
		}else if(pol == 6) {  // or1
		  if(trial == 0 || chi2 < chi2best) {
		    chi2best = chi2;
		    int j, jmax;
		    jmax = fitReceiverModel_funk_info.solution.or2model;
		    if(jmax == 1)
		      jmax = 2;
		    for(j = 0; j < jmax; j++) {
		      analyticsolution->or2model_values[j] = xfit[j];
		    }
		    //		    analyticsolution->or2model_values[0] = xfit[0];
		    //		    analyticsolution->or2model_values[1] = xfit[1];
		    analyticsolution->or2model_redchi2 = chi2/(double)fitReceiverModel_funk_info.nrPoints;
		    if(currefmodeltrial == 0 || chi2best < currefmodeltrial_chi2best) {
		      currefmodeltrial_chi2best = chi2best;
		      referencemodel_best = currefmodeltrial;
		    }
		  }
		}else {
		  fflush(stdout); 
		  printerror(verbose.debug, "ERROR fitReceiverModel: Something wrong.");
		}
	      }else {
		if(ret == 1) {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR : Maximum number of itterations exceeded (maybe try smaller ftol)");
		}else if(ret == 2) {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR : Memory allocation error");
		}else if(ret == 3) {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR : Nr of parameters to fit for should be at least 1");
		}else if(ret == 4) {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR: Specified fit algorithm is not available");
		}
	      }
	    }
	  }
	  if(dofit == 0)  // If not fitting, we can quit the throwOutOutliersItt loop
	    break;
	} // End of throwOutOutliersItt loop
	if(verbose.verbose) {
	  for(i = 0; i < verbose.indent; i++)      
	    printf(" ");
	  printf("  For pol %d, reference model trial %d, chi2 = %e\n", pol, currefmodeltrial, chi2best);
	}
      } // End of currefmodeltrial loop
    }  // End of pol loop
  } // End of subint loop
  free(xdata);
  free(ydata);
  free(ypred);
  return 1;
}

/* Replaces the data with a fake cal signal, so memory should already
   be allocated in original. The data are created as coherency
   parameters rather than Stokes parameters. A list of gains of the
   two polarization channels and phases (in radians) should be
   provided (one value per frequency channel). In addition the
   orientations and ellipticity of leakage can be given (if set to
   NULL these parameters will be ignored). If slope is set, a phase
   gradient (PA-slope) is introduced in the cal signal as function of
   pulse longitude. The default will produce a cal signal with a fixed
   PA as function of pulse longitude. rms sets the rms level of the
   noise, with the signal strength being of the order of the gains
   squared. Two levels are provided which defines a gradient over the
   frequency channels. Set rms_low=rms_high if no gradient is desired.
   If fixseed is set, the noise is reproducable. Set type=1 for the
   standard square-wave signal, type=2 for a gaussian peak. Set
   circular to the fraction of the signal that will be Stokes V. Set
   unpolarised to the fraction of the signal that will be
   unpolarised. By default a pure Q signal is generated (refpa=0). If
   refpa=pi/4 (45 deg), a pure U signal is generated instead. If qw is
   set, the signal is first passed through an imperfect quarter wave
   plate before the receiver parameters are applied. The system is
   characterized with qw_phi (rotation of basis) and the phase delay
   qw_delay, with -45,90 representing an ideal quarter wave plate.
   Verbose level detirmines nr of spaces before output.  Return 1 =
   successx */
int generateCal(datafile_definition *original, float *gain1, float *gain2, float *phases, int slope, float *el0, float *el1, float *or0, float *or1, float rms_low, float rms_high, int fixseed, int type, float circular, float unpolarised, float refpa, int qw, float qw_phi, float qw_delay, verbose_definition verbose)
{
  long f, n, b, i, idnum;
  float S[4], M[4][4], G, gamma, steeperslope, rms;
  int basis;
  gsl_rng *rand_num_gen;  /* global generator */
  const gsl_rng_type *rand_num_gen_type;

  gsl_rng_env_setup();    /* Set the default generators, can be influenced by environment variables */
  rand_num_gen_type = gsl_rng_default;
  rand_num_gen = gsl_rng_alloc (rand_num_gen_type);
  if(fixseed)
    idnum = 2;
  else
    randomize_idnum(&idnum);
  //  idnum = 0;
  //  printf("XXXXXX %ld %d\n", idnum, fixseed);
  gsl_rng_set(rand_num_gen, idnum);


  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("generate cal signal\n");
  }
  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR generateCal: no memory allocated?");
    return 0;
  }
  if(original->NrPols != 4) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR generateCal: Expected 4 polarization channels, only %ld present.", original->NrPols);
    return 0;
  }
  if(original->poltype != POLTYPE_COHERENCY) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR generateCal: Expected poltype=%d (coherency parameters), got %d.", POLTYPE_COHERENCY, original->poltype);
    return 0;
  }
  if(original->feedtype == FEEDTYPE_CIRCULAR || original->feedtype == FEEDTYPE_INV_CIRCULAR) {
    basis = 2;
  }else if(original->feedtype == FEEDTYPE_LINEAR || original->feedtype == FEEDTYPE_INV_LINEAR) {
    basis = 1;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR generateCal: Got feedtype=%d, which is not implemented.", original->feedtype);
    return 0;
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(basis == 1)
      printf("  using a linear basis\n");
    else if(basis == 2)
      printf("  using a circular basis\n");
  }
  if(qw && original->feedtype != FEEDTYPE_CIRCULAR) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR generateCal: Simulation of a quarter wave plate not implemented for a linear receiver.");
    return 0;
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA || original->poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR generateCal: Cannot handle PA data.");
    return 0;
  }

  if(unpolarised < 0 || unpolarised > 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR generateCal: Requested amount of unpolarised emission should be a fraction, %f is not a fraction.", unpolarised);
    return 0;
  }
  if(fabs(circular/(1.0-unpolarised)) > 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR generateCal: Requested amount fraction of Stokes V (%f) should be smaller than the total amount of polarised emission %f.", circular, 1.0-unpolarised);
    return 0;
  }

  circular = asin(circular/(1.0-unpolarised));
  unpolarised = acos(1.0-unpolarised);


  /* If generating gaussian peak, make slope pa-swing steeper */
  if(type == 2)
    steeperslope = 5;
  else
    steeperslope = 1;
  if(slope == 0)      /* If no slope, use rotation multiplication factor 0. */
    steeperslope = 0;

  /* Loop over pulses and frequency channels etc. */
  for(f = 0; f < original->NrFreqChan; f++) {
    if(original->NrFreqChan > 1) 
      rms = rms_low + (rms_high-rms_low)*f/(float)(original->NrFreqChan-1);
    else
      rms = rms_low;
    gamma = gain1[f]/gain2[f]-1.0;
    G = sqrt(gain1[f]*gain2[f]);
    if(el0 != NULL && el1 != NULL && or0 != NULL && or1 != NULL)
      constructMueller(M, G, gamma, phases[f], el0[f], el1[f], or0[f], or1[f], original->feedtype, 0, 0, 0, 0);
    else
      constructMueller(M, G, gamma, phases[f], 0, 0, 0, 0, original->feedtype, 0, 1, 0, 0);
    for(n = 0; n < original->NrSubints; n++) {
      for(b = 0; b < original->NrBins; b++) {
	/* Start with a 100% linearly polarized source. These are Stokes parameters at this point. */
	S[0] = 1;
	S[1] = cos(unpolarised)*cos(circular)*cos(2.0*refpa+2.0*M_PI*steeperslope*b/(float)original->NrBins);
	S[2] = cos(unpolarised)*cos(circular)*sin(2.0*refpa+2.0*M_PI*steeperslope*b/(float)original->NrBins);
	S[3] = cos(unpolarised)*sin(circular);

	/* 	if(b == 0) printf("XXXX chan=%ld: %f %f %f %f   ->  ", f, S[0], S[1], S[2], S[3]); */
	if(qw)
	  applyQW(S, qw_phi, qw_delay);
	applyMueller(M, S);
	/* 	if(b == 0) printf("%f %f %f %f\n", S[0], S[1], S[2], S[3]); */

	/* If type is 2, multiply coherency parameters with a Gaussian shaped curve. */
	if(type == 2) {
	  G = exp(-(b-0.5*original->NrBins)*(b-0.5*original->NrBins)/(0.5*0.01*original->NrBins*original->NrBins));
	  S[0] *= G;
	  S[1] *= G;
	  S[2] *= G;
	  S[3] *= G;
	}else {  /* Otherwise, by default make a square wave by setting half of the profile to zero. */
	  if(b < original->NrBins/2) {
	    S[0] = 0;
	    S[1] = 0;
	    S[2] = 0;
	    S[3] = 0;
	  }
	}
	/* Add random noise in Stokes space */
	S[0] += gsl_ran_gaussian(rand_num_gen, rms);
	S[1] += gsl_ran_gaussian(rand_num_gen, rms);
	S[2] += gsl_ran_gaussian(rand_num_gen, rms);
	S[3] += gsl_ran_gaussian(rand_num_gen, rms);

	/* Convert the Stokes parameters into coherency parameters. */
	convertStokesToCoherency(S, 2, basis);

	/* Add the white noise and store the result in the output array. */
	/* NR code:
	original->data[original->NrBins*(0+original->NrPols*(f+n*original->NrFreqChan))+b] = c1+rms*gasdev(&idnum);
	original->data[original->NrBins*(1+original->NrPols*(f+n*original->NrFreqChan))+b] = c2+rms*gasdev(&idnum);
	original->data[original->NrBins*(2+original->NrPols*(f+n*original->NrFreqChan))+b] = c3+rms*gasdev(&idnum);
	original->data[original->NrBins*(3+original->NrPols*(f+n*original->NrFreqChan))+b] = c4+rms*gasdev(&idnum);
	*/

	original->data[original->NrBins*(0+original->NrPols*(f+n*original->NrFreqChan))+b] = S[0];
	original->data[original->NrBins*(1+original->NrPols*(f+n*original->NrFreqChan))+b] = S[1];
	original->data[original->NrBins*(2+original->NrPols*(f+n*original->NrFreqChan))+b] = S[2];
	original->data[original->NrBins*(3+original->NrPols*(f+n*original->NrFreqChan))+b] = S[3];
      }
    }
  }
  gsl_rng_free (rand_num_gen);
  /*  if(verbose) printf(" done                            \n");*/
  return 1;
}


/* The conversion will only be done if poltype == POLTYPE_COHERENCY (i.e. coherency
   parameters). feedtype indicates if the receiver is linear or
   circular. */
void convertCoherencyToStokes(float *C, int poltype, int fd_type)
{
  float I, Q, U, V;
  /* If coherency parameters, they need to be converted */
  if(poltype == POLTYPE_COHERENCY) {
    if(fd_type == FEEDTYPE_CIRCULAR) {
      I = C[0]+C[1];
      Q = 2.0*C[2];
      U = 2.0*C[3];
      V = C[0]-C[1];
      C[0] = I;
      C[1] = Q;
      C[2] = U;
      C[3] = V;
    }else if(fd_type == FEEDTYPE_LINEAR) {
      I = C[0]+C[1];
      Q = C[0]-C[1];
      U = 2.0*C[2];
      V = 2.0*C[3];
      C[0] = I;
      C[1] = Q;
      C[2] = U;
      C[3] = V;
    }
  }
}

/* The conversion will only be done if poltype == POLTYPE_COHERENCY (i.e. set to
   coherency parameters). fd_type indicates if the receiver is linear
   or circular. */
void convertStokesToCoherency(float *C, int poltype, int fd_type)
{
  float c1, c2, c3, c4;
  /* If coherency parameters, they need to be converted */
  if(poltype == POLTYPE_COHERENCY) {
    if(fd_type == FEEDTYPE_CIRCULAR) {
      c1 = 0.5*(C[0]+C[3]);
      c2 = 0.5*(C[0]-C[3]);
      c3 = 0.5*C[1];
      c4 = 0.5*C[2];

      C[0] = c1;
      C[1] = c2;
      C[2] = c3;
      C[3] = c4;
    }else if(fd_type == FEEDTYPE_LINEAR) {
      c1 = 0.5*(C[0]+C[1]);
      c2 = 0.5*(C[0]-C[1]);
      c3 = 0.5*C[2];
      c4 = 0.5*C[3];
      C[0] = c1;
      C[1] = c2;
      C[2] = c3;
      C[3] = c4;
    }
  }
}


/* Applies the receiver model to the datafile. The datafile and the
   solution should be read into memory. If reconstruct is set, the
   intrinsic Stokes parameters are reconstructed using the receiver
   parameters. If not, the data is distorted. If ignoreGain is set,
   the gain parameters are ignored, which might be useful to prevent
   the data to be dominated by a few channels. If ignore_dodgy_leakage
   is set, it makes the Mueller matrix zero (i.e. it will zap data) if
   the normalisation factor in the leakage matrix is abnormally
   high. This will only be in effect together with the reconstruct
   flag. If normalise is set, the matrix is divided by the value of
   the first element, i.e. after applying the Mueller matrix to the
   data the rms of Stokes I in the noise is unchanged, assuming that
   the rms of the four Stokes parameters are equal to start off
   with. Normalise is only identical to setting the gain to 1 for an
   ideal receiver. Return 1 = success */
int applyReceiverModel(datafile_definition *datafile, datafile_definition solution, int reconstruct, int ignoreGain, int ignore_dodgy_leakage, int normalise, verbose_definition verbose)
{
  long f, n, b, i, index, index2;
  int leakage;
  float S[4], M[4][4], phase, gamma, G, el0, el1, or0, or1;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Applying receiver model to data\n");
  }
  if(datafile->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR applyReceiverModel: datafile should be read into memory first");
    return 0;
  }
  if(datafile->NrPols != 4) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR applyReceiverModel: Expected 4 polarization channels, only %ld present.", datafile->NrPols);
    return 0;
  }
  if(datafile->poltype != POLTYPE_COHERENCY && datafile->poltype != POLTYPE_STOKES) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR applyReceiverModel: Expected poltype=%d or %d (Stokes or coherency parameters), got %d.", POLTYPE_STOKES, POLTYPE_COHERENCY, datafile->poltype);
    return 0;
  }
  if(verbose.verbose) {
    if(datafile->poltype == POLTYPE_STOKES) {
      printf("  Got Stokes parameters\n");
    }else if(datafile->poltype == POLTYPE_COHERENCY) {
      if(verbose.verbose) printf("  Got Coherency parameters\n");
    }
  }
  if (datafile->feedtype != FEEDTYPE_LINEAR && datafile->feedtype != FEEDTYPE_CIRCULAR) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR applyReceiverModel: Don't know how to process feedtype=%d.", datafile->feedtype);
    return 0;
  }
  if(datafile->NrFreqChan != solution.NrFreqChan) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR applyReceiverModel: The number of frequency channels are different in data file and receiver solution.");    
    return 0;
  }
  if(solution.NrPols < 3) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR applyReceiverModel: Expecteded at least 3 receiver parameters, only %ld present.", solution.NrPols);
    return 0;
  }
  if(solution.NrPols > 3 && ((solution.NrPols != 7 && solution.gentype == GENTYPE_RECEIVERMODEL) || (solution.NrPols != 9 && solution.gentype == GENTYPE_RECEIVERMODEL2))) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING applyReceiverModel: Only gain, differential gain and phase difference will be applied.");
  }


  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(datafile->feedtype == FEEDTYPE_LINEAR || datafile->feedtype == FEEDTYPE_INV_LINEAR)
      printf("  using a linear basis\n");
    else if(datafile->feedtype == FEEDTYPE_CIRCULAR || datafile->feedtype == FEEDTYPE_INV_CIRCULAR)
      printf("  using a circular basis\n");
  }
  if(datafile->poltype == POLTYPE_ILVPAdPA || datafile->poltype == POLTYPE_PAdPA || datafile->poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR applyReceiverModel: Cannot handle PA data.");
    return 0;
  }

  if((solution.NrPols >= 7 && solution.gentype == GENTYPE_RECEIVERMODEL) || (solution.NrPols >= 9 && solution.gentype == GENTYPE_RECEIVERMODEL2)) {
    leakage = 1;
  }else {
    leakage = 0;
  }     

  /* Loop over pulses and frequency channels etc. */
  for(f = 0; f < datafile->NrFreqChan; f++) {
    /* Get PSRCHIVE receiver solution parameters and convert gamma = exp(2*beta)-1 */
    index = solution.NrBins*(0+solution.NrPols*(f+0*solution.NrFreqChan))+0;
    if(ignoreGain)
      G = 1;
    else
      G     =       solution.data[index];
    index += solution.NrBins;
    gamma = exp(2*solution.data[index])-1;
    index += solution.NrBins;
    phase =       solution.data[index];
    if(leakage) {
      index += solution.NrBins;
      el0 =     solution.data[index];
      index += solution.NrBins;
      or0 =     solution.data[index];
      index += solution.NrBins;
      el1 =     solution.data[index];
      index += solution.NrBins;
      or1 =     solution.data[index];
      if(constructMueller(M, G, gamma, phase, el0, el1, or0, or1, datafile->feedtype, reconstruct, 0, ignore_dodgy_leakage, normalise) == 0) {
	printwarning(verbose.debug, "WARNING applyReceiverModel: Dodgy solution results in frequency channel %d is set to zero", f);
      }
    }else {
      constructMueller(M, G, gamma, phase, el0, el1, or0, or1, datafile->feedtype, reconstruct, 1, ignore_dodgy_leakage, normalise);
    }
    for(n = 0; n < datafile->NrSubints; n++) { 
      index = datafile->NrBins*(0+datafile->NrPols*(f+n*datafile->NrFreqChan));
      for(b = 0; b < datafile->NrBins; b++) {
	index2 = index;
	S[0] = datafile->data[index2];
	index2 += datafile->NrBins;
	S[1] = datafile->data[index2];
	index2 += datafile->NrBins;
	S[2] = datafile->data[index2];
	index2 += datafile->NrBins;
	S[3] = datafile->data[index2];
	//	if(b==0) printf("XXXXXXXX %f\n", S[0]);

	/* Convert from coherency parameters into Stokes parameters if necessary */
	convertCoherencyToStokes(S, datafile->poltype, datafile->feedtype);

      
	/* Apply the corrections */
	applyMueller(M, S);
	/*	if(b==0) printf("XXXXXXXX %f\n", S[0]); */

	/* Not sure why, but for a linear feed we have to add an extra rotation in Q and U. Needs some more testing, so took it out again. */
	/*	  if(datafile->feedtype == FEEDTYPE_LINEAR) 
		  calibrate_applyQURotation(S, -45*M_PI/180.0);
	*/

	/* Convert back into coherency parameters (if necessary)*/
	convertStokesToCoherency(S, datafile->poltype, datafile->feedtype);

	/* If result is nan or infinity (most likely because there is no valid receiver model for that frequency channel), set intensity Stokes parameters to 0. */
	if(isnan(S[0]) || isnan(S[1]) || isnan(S[2]) || isnan(S[3]) || isinf(S[0]) || isinf(S[1]) || isinf(S[2]) || isinf(S[3])) {
	  S[0] = S[1] = S[2] = S[3] = 0;
	}
	index2 = index;
	datafile->data[index2] = S[0];
	index2 += datafile->NrBins;
	datafile->data[index2] = S[1];
	index2 += datafile->NrBins;
	datafile->data[index2] = S[2];
	index2 += datafile->NrBins;
	datafile->data[index2] = S[3];

	index++;

	/*	if(b==0) printf("XXXXXXXX %f\n", S[0]);   */
      }
    }
  }
  /*  if(verbose) printf(" done                            \n");*/
  return 1;
}

/* Prints the Mueller matrix constructed from the receiver model.  Set
   fd_type=1 for a linear basis, 2 for a circular basis. If singlerow
   is set, the matrix elements appear on a single row. The solution
   should be read into memory. If reconstruct is set to 1, the
   intrinsic Stokes parameters are reconstructed using the receiver
   parameters. If set to zero, the data is distorted. If reconstruct
   is set to 2, the Mueller matrix is multiplied with its inverse, to
   test if the unity matrix is obtained. If ignoreGain is set, the
   gain parameters are ignored. If normalise is set, the matrix is
   divided by the value of the first element, i.e. after applying the
   Mueller matrix to the data the rms of Stokes I in the noise is
   unchanged, assuming that the rms of the four Stokes parameters are
   equal to start off with. Normalise is only identical to setting the
   gain to 1 for an ideal receiver. If stat is set, some statistics
   are returned as well. Return 1 = success, 0 = error. */
int printMuellerFromReceiverModel(datafile_definition solution, int reconstruct, int ignoreGain, int normalise, int fd_type, int singlerow, int stat, verbose_definition verbose)
{
  long f, i, index;
  int leakage, ignore_dodgy_leakage, j, k;
  float M[4][4], M2[4][4], phase, gamma, G, el0, el1, or0, or1;
  double overall_channelmax_leakage, overall_channelmax_crosspol;

  // Show all solutions, even if looking dodgy
  ignore_dodgy_leakage = 0;

  if(reconstruct < 0 || reconstruct > 2) {
    printerror(verbose.debug, "ERROR printMuellerFromReceiverModel: The parameter reconstruct is an invalid number: %d.", reconstruct);
    return 0;
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Output construncted Mueller matrix from receiver model. ");
    if(reconstruct == 1)
      printf("This is the matrix to correct the data\n");
    else if(reconstruct == 0)
      printf("This is the matrix which has distorted the data\n");
    else
      printf("This should be the unity matrix\n");
  }
  if(solution.NrPols < 3) {
    printerror(verbose.debug, "ERROR printMuellerFromReceiverModel: Expecteded at least 3 receiver parameters, only %ld present.", solution.NrPols);
    return 0;
  }
  if(solution.NrPols > 3 && ((solution.NrPols != 7 && solution.gentype == GENTYPE_RECEIVERMODEL) || (solution.NrPols != 9 && solution.gentype == GENTYPE_RECEIVERMODEL2))) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING printMuellerFromReceiverModel: Only gain, differential gain and phase difference will be used.");
  }

  if((solution.NrPols >= 7 && solution.gentype == GENTYPE_RECEIVERMODEL) || (solution.NrPols >= 9 && solution.gentype == GENTYPE_RECEIVERMODEL2)) {
    leakage = 1;
  }else {
    leakage = 0;
  }     

  overall_channelmax_leakage = 0;
  overall_channelmax_crosspol = 0;
  /* Loop over pulses and frequency channels etc. */
  for(f = 0; f < solution.NrFreqChan; f++) {
    /* Get PSRCHIVE receiver solution parameters and convert gamma = exp(2*beta)-1 */
    index = solution.NrBins*(0+solution.NrPols*(f+0*solution.NrFreqChan))+0;
    if(ignoreGain)
      G = 1;
    else
      G     =       solution.data[index];
    index += solution.NrBins;
    gamma = exp(2*solution.data[index])-1;
    index += solution.NrBins;
    phase =       solution.data[index];
    if(leakage) {
      index += solution.NrBins;
      el0 =     solution.data[index];
      index += solution.NrBins;
      or0 =     solution.data[index];
      index += solution.NrBins;
      el1 =     solution.data[index];
      index += solution.NrBins;
      or1 =     solution.data[index];
      if(constructMueller(M, G, gamma, phase, el0, el1, or0, or1, fd_type, reconstruct, 0, ignore_dodgy_leakage, normalise) == 0) {
	printwarning(verbose.debug, "WARNING applyReceiverModel: Dodgy solution results in frequency channel %d is set to zero", f);
      }
      if(reconstruct == 2) { // Also need the inverse
	//	printf("XXXXX %d -> %d\n", reconstruct, !reconstruct);
	if(constructMueller(M2, G, gamma, phase, el0, el1, or0, or1, fd_type, 0, 0, ignore_dodgy_leakage, normalise) == 0) {
	  printwarning(verbose.debug, "WARNING applyReceiverModel: Dodgy solution results in frequency channel %d is set to zero", f);
	}
      }
    }else {
      constructMueller(M, G, gamma, phase, el0, el1, or0, or1, fd_type, reconstruct, 1, ignore_dodgy_leakage, normalise);
      if(reconstruct == 2) { // Also need the inverse
	constructMueller(M2, G, gamma, phase, el0, el1, or0, or1, fd_type, 0, 1, ignore_dodgy_leakage, normalise);
      }
    }
    if(reconstruct == 2) {
      int l;
      float M3[4][4];
      for(j = 0; j < 4; j++) {
	for(k = 0; k < 4; k++) {
	  M3[j][k] = 0;
	  for(l = 0; l < 4; l++) {
	    M3[j][k] += M[j][l]*M2[l][k];
	  }
	}
      }
      for(j = 0; j < 4; j++) {
	for(k = 0; k < 4; k++) {
	  M[j][k] = M3[j][k];
	}
      }
    }
    printf("Mueller matrix for channel %ld", f);
    if(singlerow)
      printf(": ");
    else
      printf("\n");
    for(j = 0; j < 4; j++) {
      for(k = 0; k < 4; k++) {
	printf("%e   ", M[j][k]);
      }
      if(singlerow == 0)
	printf("\n");
    }
    if(singlerow)
      printf("\n");
    if(stat) {
      double channelmax_leakage, channelmax_crosspol;
      channelmax_leakage = 0;
      channelmax_crosspol = 0;
      for(j = 1; j < 4; j++) {
	if(fabs(M[0][j]) > channelmax_leakage)
	  channelmax_leakage = fabs(M[0][j]);
	if(fabs(M[j][0]) > channelmax_leakage)
	  channelmax_leakage = fabs(M[j][0]);
      }
      for(j = 1; j < 4; j++) {
	for(k = 1; k < 4; k++) {
	  if(j != k) {
	    if(fabs(M[j][k]) > channelmax_crosspol) {
	      channelmax_crosspol = fabs(M[j][k]);
	    }
	  }
	}
      }
      if(channelmax_leakage > overall_channelmax_leakage)
	overall_channelmax_leakage = channelmax_leakage;
      if(channelmax_crosspol > overall_channelmax_crosspol)
	overall_channelmax_crosspol = channelmax_crosspol;
      printf("Maximum conversion factor to/from Stokes I (leakage): %lf\n", channelmax_leakage);
      printf("Maximum conversion factor to/from Stokes Q/U/V: %lf\n", channelmax_crosspol);
    }
  }
  if(stat) {
    printf("\n");
    printf("Summary of all frequency channels:\n");
    printf("Maximum conversion factor to/from Stokes I (leakage) when considering the first row/column of the matrix: %lf\n", overall_channelmax_leakage);
    printf("Maximum conversion factor to/from Stokes Q/U/V when considering the other elements of the matrix: %lf\n", overall_channelmax_crosspol);
  }

  return 1;
}

void cal_internal_drawTwoProfiles(float *profile1, float *profile2, long NrBins, char *title)
{
  long bin;
  float miny, maxy;
#ifdef HAVEPGPLOT2 
  ppgopen("/pw");
#else
  ppgopen("/xs");
#endif
  ppgask(0);
  ppgslw(1);
  ppgpage();
  ppgsvp(0.1, 0.9, 0.1, 0.9);
  miny = maxy = profile1[0];
  for(bin=0; bin < NrBins; bin++) {
    if(profile1[bin] > maxy)
      maxy = profile1[bin];
    if(profile1[bin] < miny)
      miny = profile1[bin];
    if(profile2[bin] > maxy)
      maxy = profile2[bin];
    if(profile2[bin] < miny)
      miny = profile2[bin];
  }
  miny -= 0.05*(maxy-miny);
  maxy += 0.05*(maxy-miny);
  ppgswin(0, NrBins, miny, maxy);
  ppgsci(1);
  ppgbox("bcnsti",0.0,0,"bcntsi",0.0,0);
  ppglab("Bins", "Stokes I", title);
  ppgmove(0, profile1[0]);
  for(bin=0; bin < NrBins; bin++)
    ppgdraw(bin, profile1[bin]);
  ppgsci(2);
  ppgmove(0, profile2[0]);
  for(bin=0; bin < NrBins; bin++)
    ppgdraw(bin, profile2[bin]);
  ppgsci(1);
  ppgend();
}


// Set this flag if you want to find the solution by correcting the observation rather than 
// distorting the template. The latter is the best if the S/N of the observation < template.
int internal_mtm_inverse;
// Solution containing the solution of the current frequency channel being processed
datafile_definition internal_solution_cur_channel;
// Template data, onpulse region only
datafile_definition internal_template_channel_onpulse;
// Uncalibrated data, onpulse region only 
datafile_definition internal_data_channel_onpulse;
// Calibrated data, onpulse region only 
datafile_definition internal_data_channel_onpulse_calib;
// The distorted template.
datafile_definition internal_template_channel_onpulse_calib;
// Set if the leakage is fitted for
int internal_leakage_fit_flag;
// Set if the gain will be analytically determined rather than treat it as a fit parameter
int internal_analytic_gain_calc_flag;

/* Position is number between -1 and 3. If device != NULL, it will be the pgplot device to use. If onlyone is set, a wider box is generated. */
void cal_internal_drawBeforeAfter(datafile_definition profile, int position, int onlyone, char *title, char *device, int profileresid_device)
{
  long bin, pol;
  static float miny, maxy;

  if(position == -1 || position == 0 || position == 2) {
    if(position == 0) {
      if(device == NULL) {
	if(profileresid_device == 0) {
#ifdef HAVEPGPLOT2 
	  ppgopen("/pw");
#else
	  ppgopen("/xs");
#endif
	}else {
	  ppgslct(profileresid_device);
	}
      }else {
	ppgopen(device);
      }
      ppgask(0);
      ppgpage();
    }
    if(position != -1) {
      ppgslw(1);
      if(position == 0) {
	if(onlyone) {
	  ppgsvp(0.1, 0.9, 0.1, 0.9);
	}else {
	  ppgsvp(0.1, 0.45, 0.1, 0.9);
	}
      }else {
	ppgsvp(0.55, 0.9, 0.1, 0.9);      
      }
    }
    if(position == -1 || position == 2) {
      miny = maxy = profile.data[0];
    }
    for(pol = 0; pol < 4; pol++) {
      for(bin=0; bin < profile.NrBins; bin++) {
	if(profile.data[bin+profile.NrBins*pol] > maxy)
	  maxy = profile.data[bin+profile.NrBins*pol];
	if(profile.data[bin+profile.NrBins*pol] < miny)
	  miny = profile.data[bin+profile.NrBins*pol];
      }
    }
    if(position == -1)
      return;
    miny -= 0.1*(maxy-miny);
    maxy += 0.3*(maxy-miny);
    ppgswin(0, profile.NrBins-1, miny, maxy);
    ppgsci(1);
    ppgbox("bcnsti",0.0,0,"bcntsi",0.0,0);
    ppglab("Bins", "Stokes I/Q/U/V", title);
  }

  for(pol = 0; pol < 4; pol++) {
    switch(position) {
    case 0: ppgsci(1+pol); ppgslw(5); break;
    case 1: ppgsci(1+pol); ppgslw(15); break;
    case 2: ppgsci(1+pol); ppgslw(5); break;
    case 3: ppgsci(1+pol); ppgslw(15); break;
    }
    for(bin=0; bin < profile.NrBins; bin++) {
      if(position == 0 || position == 2) {
	if(bin == 0)
	  ppgmove(bin, profile.data[bin+profile.NrBins*pol]);
	else
	  ppgdraw(bin, profile.data[bin+profile.NrBins*pol]);
      }else {
	ppgpt1(bin, profile.data[bin+profile.NrBins*pol], -1);
      }
    }
  }

  if(position == 3) {
    if(profileresid_device == 0) {
      ppgend();
    }
  }
}


/* This function takes the calibrated solution and rescales the
   intensity scale such that the integral over the onpulse region of
   Stokes I is the same for the template and data set. */
void normalise_solution_with_template(void)
{
  long b;
  double intensity_template, intensity_data;
  float norm;

  intensity_template = 0;
  intensity_data = 0;
  for(b = 0; b < internal_data_channel_onpulse_calib.NrBins; b++) {
    intensity_template += internal_template_channel_onpulse.data[b];
    intensity_data += internal_data_channel_onpulse_calib.data[b];
  }

  norm = intensity_template/intensity_data;
  for(b = 0; b < 4*internal_data_channel_onpulse_calib.NrBins; b++) {
    internal_data_channel_onpulse_calib.data[b] *= norm;
  }
}

/* This function returns a figure of merrit which will be minimised by
   changing the fitparameters (receiver solution). */
double funk_optimizeMTM(double *fitparams)
{
  long b, p, i;
  double chi2, delta; 
  verbose_definition noverbose;
  cleanVerboseState(&noverbose);
  noverbose.nocounters = 1;

  // Prepare a single channel solution. If the flag is set, we're not
  // fitting for the gain, but it should be determined analytically.
  if(internal_analytic_gain_calc_flag)
    internal_solution_cur_channel.data[0] = 1.0;
  else
    internal_solution_cur_channel.data[0] = fitparams[0];
  internal_solution_cur_channel.data[1] = fitparams[1];
  internal_solution_cur_channel.data[2] = fitparams[2];
  if(internal_leakage_fit_flag) {
    internal_solution_cur_channel.data[3] = fitparams[3];
    internal_solution_cur_channel.data[4] = fitparams[4];
    internal_solution_cur_channel.data[5] = fitparams[5];
    internal_solution_cur_channel.data[6] = fitparams[6];
  }
 
  // Make a copy of the uncalibrated data. The copy will be replaced with the calibrated result.
  if(internal_mtm_inverse) {
    memcpy(internal_data_channel_onpulse_calib.data, internal_data_channel_onpulse.data, (internal_data_channel_onpulse.NrBins)*(internal_data_channel_onpulse.NrPols)*(internal_data_channel_onpulse.NrFreqChan)*(internal_data_channel_onpulse.NrSubints)*sizeof(float));
    if(output_FittedProfile_Of_Itterations) {
      i = 0;
      for(p = 0; p < 4; p++) { 
	for(b = 0; b < internal_data_channel_onpulse_calib.NrBins; b++) {
	  printf("FITITTERATION_DATA_BEFORE: %f\n", internal_data_channel_onpulse_calib.data[i]);
	  i++;
	}
      }
    }
    if(applyReceiverModel(&internal_data_channel_onpulse_calib, internal_solution_cur_channel, 1, 0, 0, 0, noverbose) == 0) {
    fflush(stdout);
      printerror(0, "ERROR cal_mtm: Cannot apply solution in present itteration.");
      exit(0);
    }
    if(internal_analytic_gain_calc_flag)
      normalise_solution_with_template();
  }else {
    memcpy(internal_template_channel_onpulse_calib.data, internal_template_channel_onpulse.data, (internal_data_channel_onpulse.NrBins)*(internal_data_channel_onpulse.NrPols)*(internal_data_channel_onpulse.NrFreqChan)*(internal_data_channel_onpulse.NrSubints)*sizeof(float));
    if(output_FittedProfile_Of_Itterations) {
      i = 0;
      for(p = 0; p < 4; p++) { 
	for(b = 0; b < internal_data_channel_onpulse_calib.NrBins; b++) {
	  printf("FITITTERATION_TEMPLATE_BEFORE: %f\n", internal_template_channel_onpulse_calib.data[i]);
	  i++;
	}
      }
    }
    if(applyReceiverModel(&internal_template_channel_onpulse_calib, internal_solution_cur_channel, 0, 0, 0, 0, noverbose) == 0) {
      fflush(stdout);
      printerror(0, "ERROR cal_mtm: Cannot apply solution in present itteration.");
      exit(0);
    }
  }


  /* Calculate a figure of merit, based on a classic chi^2
calculation. However, this FOM does not includes a normalisation by
the variance. Including this normalisation will bias the gain such
that the calibrated intensities are smaller. This is because this
makes the variance smaller, hence the chi^2 better. So instead use the
unnormalised chi^2 as a figure of merrit. This should be fine, as long
as the variance on the Stokes parameters is the same. Later in the
code a proper chi^2 is determined for the best solution. Note that the
absence of this normalisation also makes the code faster, especially
because only the onpulse region needs to be processed in each
itteration.
   */
  chi2 = 0;
  i = 0;
  for(p = 0; p < 4; p++) { 
    for(b = 0; b < internal_data_channel_onpulse_calib.NrBins; b++) {
      if(internal_mtm_inverse) {
	delta = (double)internal_data_channel_onpulse_calib.data[i]-(double)internal_template_channel_onpulse.data[i];
      }else {
	delta = (double)internal_data_channel_onpulse.data[i]-(double)internal_template_channel_onpulse_calib.data[i];
      }
            chi2 += delta*delta;
      //      chi2 += fabs(delta*delta*delta);
      //      chi2 += fabs(delta);
      i++;
    }
  }

  /* Avoid that the differential phase crosses -90 or +90 boundary, as it wraps every 180 degrees. */
  if(internal_solution_cur_channel.data[2] > 0.5*M_PI || internal_solution_cur_channel.data[2] < -0.5*M_PI)
    chi2 *= 1e10;
  /* Ellipticities/orientations do not wrap after 180 deg, but there
     is an equivalent solution with a different differential phase
     when 180 deg are added. In addition if we take or -> or + 90 and
     at same time el -> 90-el and dphase -> dphase - 45 an covarient
     solution is found with the differential phase. Hence, we can take
     ellipticities in limit +-45 deg. */
  if(internal_leakage_fit_flag) {
    if(internal_solution_cur_channel.data[3] > 0.25*M_PI || internal_solution_cur_channel.data[3] < -0.25*M_PI)
      chi2 *= 1e10;
    if(internal_solution_cur_channel.data[4] > 0.5*M_PI || internal_solution_cur_channel.data[4] < -0.5*M_PI)
      chi2 *= 1e10;
    if(internal_solution_cur_channel.data[5] > 0.25*M_PI || internal_solution_cur_channel.data[5] < -0.25*M_PI)
      chi2 *= 1e10;
    if(internal_solution_cur_channel.data[6] > 0.5*M_PI || internal_solution_cur_channel.data[6] < -0.5*M_PI)
      chi2 *= 1e10;
  }

  if(output_FittedProfile_Of_Itterations) {
    printf("FITITTERATION: FOM of itteration = %lf\n", chi2);
    i = 0;
    for(p = 0; p < 4; p++) { 
      for(b = 0; b < internal_data_channel_onpulse_calib.NrBins; b++) {
	if(internal_mtm_inverse) {
	  printf("FITITTERATION_DATA_TEMPLATE: %f %f\n", internal_data_channel_onpulse_calib.data[i], internal_template_channel_onpulse.data[i]);
	}else {
	  printf("FITITTERATION_DATA_TEMPLATE: %f %f\n", internal_data_channel_onpulse.data[i], internal_template_channel_onpulse_calib.data[i]);
	}
	i++;
      }
    }
    for(p = 0; p <= 6; p++) {
      printf("FITITTERATION_PARAMS: %f\n", fitparams[p]);
    }
    printf("FITITTERATION_MUELLER_DISTORT: ");
    printMuellerFromReceiverModel(internal_solution_cur_channel, 0, 0, 0, internal_template_channel_onpulse_calib.feedtype, 1, 0, noverbose);
    printf("FITITTERATION_MUELLER_CORRECT: ");
    printMuellerFromReceiverModel(internal_solution_cur_channel, 1, 0, 0, internal_template_channel_onpulse_calib.feedtype, 1, 0, noverbose);
  }

  return chi2;
}

int getinitialguess_internal(datafile_definition *initialguess, long pol, long freq, double *value, verbose_definition verbose)
{
  int nrpols;
  float valuef;
  if(initialguess != NULL) {
    nrpols = initialguess->NrPols;
    if(initialguess->gentype == GENTYPE_RECEIVERMODEL2)
      nrpols -= 2;
    if(pol < nrpols) {
      if(readPulsePSRData(initialguess, 0, pol, freq, 0, 1, &valuef, verbose) != 1) {
	fflush(stdout);
	printerror(0, "ERROR cal_mtm: Cannot read from initial guess receiver model.");
	return 0;
      }
      *value = valuef;
      return 1;
    }
  }
  *value = 0;
  return 1;
}

/* This function fits the data with the template. Note that both
datasets could be destroyed, as data might need to be rebinned and
shifted and dedispersed etc. The onpulse region can be set with
onpulseregion for noise calculation. If set to NULL, the user will be
asked to select it. This will also happen when the number of regions
is 0, but in that case the region will be returned (useful when
processing multiple files). onpulseregion_fit is similar, but then for
the region that is used to do the actual fit. If aligned is set, the
template and observation are not aligned using a cross-correlation.
If initialguess != NULL, this receiver model is used to set the
initial guess. The solution is written as a fits file to
outputfilename. argc and argv should be provided to allow writing of a
history line to the output file. If debug is set, a lot of extra
information is produced, including plots (the before/after plots are
sent to beforeafter_device). Instead, if profileresid_device != 0, the
plots are sent to this pgplot id instead, which allows plots from
subsequent calls to cal_mtm() to be appendend. If leakage is set, the
leakage terms will be fitted for as well, otherwise just the gain
imbalance. The non-zero elements in the array called fixed (which has
to contain 7 elements) indicates which parameters have to be fixed to
the initial guess value. If mtm_inverse is set, the observation is
corrected during the fitting process to match the template. This is in
general a bad idea, because if the S/N of the data is low, the
amplitude of the Stokes parameters is underestimated. By default the
template is distorted to match the observation. If analytic_gain is
set, the gain is not treated as a free parameter, but it is fixed to
whatever is needed to make the average onpulse Stokes I match the
template. This resolves the bias which makes the amplitude of the
observation lower than the template (in low S/N situations), but there
is still a bias in the determined receiver solution. If bruteforce is
zero, only 1 initial value for each receiver parameter is
tried. Otherwise the number of trials are (1+2*bruteforce)**(nr fit
parameters). Hence this variable takes bruteforce extra trials at each
side of the central value, which is the ideal receiver value. Note
that the number of trials are increasing extremely rapidly, so be
careful! A somewhat quicker way is to set extremes, in which case for
all angles 0 and the extreme values are tried, without any
intermediate values. If sequencenr is set to a negative number, a
dumpfile will be generated, containing all the required information to
run this process in multiple instances, each instance processing part
of the band. If totsequencenr is set to 0, cal_mtm processes the band
in one go without generating a dump file. Set sequencenr to
1...totsequencenr to process part of the band. When part of the band
is processed, the datafile, template and initialguess should not be
allocated yet, as data will be read from dumpfile and memory will be
allocated. If nocounters is set, no counters are shown. If
ignoreParallacticAngle is set, all parallactic angle corrections are
ignored. This function returns 0 when there is an error. */
int cal_mtm(datafile_definition *datafile, datafile_definition *template, pulselongitude_regions_definition *onpulseregion, pulselongitude_regions_definition *onpulseregion_fit, int aligned, datafile_definition *initialguess, int leakage, int *fixed, int analytic_gain, int mtm_inverse, int bruteforce, int extremes, char *outputfilename, int argc, char **argv, char *dumpfilename, int sequencenr, int totsequencenr, int ignoreParallacticAngle, verbose_definition verbose, char *beforeafter_device, int profileresid_device)
{
  int iszapped, ok, failed, ret, dummyi, internal_nfree;
  long i, n, f, p, b, ftemplate, np1, np2, np3, np4, np5, np6, np7, np1tot, np2tot, np3tot, np4tot, np5tot, np6tot, np7tot;
  long nrbinsNew_template, nrbinsNew_datafile, nractivebins_offpulse;
  double delay;
  float *template_channel_full;
  datafile_definition template_fscr, clone, datafile_fscr, data_channel_full;
  datafile_definition solutionall, fscrunched_calib_obs, fscrunched_template;
  datafile_definition data_channel_calib_best;
  pulselongitude_regions_definition onpulseregion_tmp, onpulseregion_tmp2;
  FILE *dumpfile;
  int dumpfile_version;
  char txt[1000];
  float fittedparamsbest[7], ftol, chisqbest_norm;
  double startx[7], dx[7], fittedparams[7], chisq, chisqbest, data_fom_best, newFOM, chisq_norm, dplus[7], dmin[7];
  int nfunk, finderrors;
  pgplot_options_definition *pgplot_options;
  verbose_definition verbose2, verbose3, noverbose, noverbose_nocounters;

  /* To make messages from other functions to appear with more spaces in front. */
  cleanVerboseState(&noverbose);
  noverbose.debug = verbose.debug;
  cleanVerboseState(&noverbose_nocounters);
  noverbose_nocounters.nocounters = 1;
  noverbose_nocounters.debug = verbose.debug;
  copyVerboseState(verbose, &verbose2);
  copyVerboseState(verbose, &verbose3);
  verbose2.indent = verbose.indent + 4;
  verbose3.indent = verbose.indent + 6;
  output_FittedProfile_Of_Itterations = 0;  // Default: don't output fitted profiles during amoeba search, can be triggered for specific conditions in code below
  if(verbose.verbose) {
    printf("\n");
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("MTM: Fit template to observation\n");
  }

  if(sequencenr == 0 && totsequencenr > 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm: instance index starts counting at 1, not 0.");
    return 0;
  }

  if(initPulselongitudeRegion(&onpulseregion_tmp, verbose) == 0) {
    printerror(verbose.debug, "ERROR cal_mtm: Initialising onpulse region failed.");
    return 0;
  }
  if(initPulselongitudeRegion(&onpulseregion_tmp2, verbose) == 0) {
    printerror(verbose.debug, "ERROR cal_mtm: Initialising onpulse region failed.");
    return 0;
  }

  internal_analytic_gain_calc_flag = analytic_gain;
  internal_mtm_inverse = mtm_inverse;

  /* Start writing dumpfile to prepare running as multiple processes */
  if(sequencenr < 0 && totsequencenr > 0) {
    dumpfile = fopen(dumpfilename, "wb");
    if(dumpfile == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Cannot open %s", dumpfilename);
      return 0;
    }
    memset(txt, 0, 9);
    sprintf(txt, "pcaldump");
    ret = fwrite(txt, 1, 8, dumpfile);
    if(ret != 8) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
      return 0;
    }
    dumpfile_version = 1;
    ret = fwrite(&dumpfile_version, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
      return 0;
    }
    if(initialguess == NULL)
      dummyi = 0;
    else
      dummyi = 1;
    ret = fwrite(&dummyi, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
      return 0;
    }
    if(dummyi) {
      initialguess->fptr_hdr = dumpfile;
      if(writePSRSALSAHeader(initialguess, verbose) == 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
	return 0;
      }
      // SHOULD DATA NOT BE WRITTEN OUT?????
    }
    ret = fwrite(&leakage, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
      return 0;
    }
    ret = fwrite(fixed, sizeof(int), 7, dumpfile);
    if(ret != 7) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
      return 0;
    }
    ret = fwrite(&bruteforce, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
      return 0;
    }
    ret = fwrite(&extremes, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
      return 0;
    }
    ret = fwrite(&analytic_gain, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
      return 0;
    }
    ret = fwrite(&internal_mtm_inverse, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
      return 0;
    }
    ret = fwrite(&totsequencenr, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
      return 0;
    }
  }

  /* Start reading dumpfile */
  if(sequencenr > 0 && totsequencenr > 0) {
    dumpfile = fopen(dumpfilename, "rb");
    if(dumpfile == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Cannot open %s", dumpfilename);
      return 0;
    }
    memset(txt, 0, 9);
    ret = fread(txt, 1, 8, dumpfile);
    if(ret != 8) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s", dumpfilename);
      return 0;
    }
    if(strcmp(txt, "pcaldump") != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s - not a dump file?", dumpfilename);
      return 0;
    }
    ret = fread(&dumpfile_version, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s", dumpfilename);
      return 0;
    }
    if(dumpfile_version != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s - not a valid dump file version?", dumpfilename);
      return 0;
    }
    ret = fread(&dummyi, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s", dumpfilename);
      return 0;
    }
    if(dummyi) {
      cleanPSRData(initialguess, verbose);
      initialguess->fptr_hdr = dumpfile;
      if(readPSRSALSAHeader(initialguess, 1, verbose) == 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: Read error from %s", dumpfilename);
	return 0;
      }
      initialguess->format = MEMORY_format;
      // SHOULD DATA NOT BE READ IN?????
    }else {
      initialguess = NULL;
    }
    ret = fread(&leakage, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s", dumpfilename);
      return 0;
    }
    ret = fread(fixed, sizeof(int), 7, dumpfile);
    if(ret != 7) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s", dumpfilename);
      return 0;
    }
    ret = fread(&bruteforce, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s", dumpfilename);
      return 0;
    }
    ret = fread(&extremes, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s", dumpfilename);
      return 0;
    }
    ret = fread(&analytic_gain, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s", dumpfilename);
      return 0;
    }
    internal_analytic_gain_calc_flag = analytic_gain;
    ret = fread(&internal_mtm_inverse, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s", dumpfilename);
      return 0;
    }
    internal_mtm_inverse = mtm_inverse;
    ret = fread(&totsequencenr, sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s", dumpfilename);
      return 0;
    }
    if(sequencenr > totsequencenr && totsequencenr > 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: instance index cannot exceed total number of instances.");
      return 0;
    }
  }

  if(extremes && bruteforce) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm: Cannot use brute force and search extremes approach at the same time.");
    return 0;
  }

  /* The gain cannot be fixed to allow for, for instance, scintillation. */
  fixed[0] = 0;
  // Unless we force the gain to be whatever fits Stokes I, in which case the gain no longer is a free parameter
  if(internal_analytic_gain_calc_flag)
    fixed[0] = 1;
  if(verbose.verbose) {
    if(leakage == 0) {
      for(p = 0; p < 3; p++) {
	if(fixed[p]) {
	  for(i = 0; i < verbose.indent; i++)      
	    printf(" ");
	  printf("  Parameter %ld is fixed\n", p);
	}
      }
    }else {
      for(p = 0; p < 7; p++) {
	if(fixed[p]) {
	  for(i = 0; i < verbose.indent; i++)      
	    printf(" ");
	  printf("  Parameter %ld is fixed\n", p);
	}
      }
    }
  }

  pgplot_options = (pgplot_options_definition *)malloc(sizeof(pgplot_options_definition));
  if(pgplot_options == NULL) {
    printerror(verbose.debug, "ERROR cal_mtm: Memory allocation error");
    return 0;
  }

  /* Do some checking, rebinning etc, unless reading in dump file */
  if(sequencenr <= 0) {
    if(datafile->format != MEMORY_format) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: datafile should be read into memory first");
      return 0;
    }
    if(template->format != MEMORY_format) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: datafile should be read into memory first");
      return 0;
    }
    if(datafile->NrPols != 4) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Expected 4 polarization channels, only %ld present in data.", datafile->NrPols);
      return 0;
    }
    if(template->NrPols != 4) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Expected 4 polarization channels, only %ld present in template", template->NrPols);
      return 0;
    }
    if(datafile->poltype != POLTYPE_STOKES) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Expected poltype=%d (Stokes parameters), got %d in datafile", POLTYPE_STOKES, datafile->poltype);
      return 0;
    }
    if(template->poltype != POLTYPE_STOKES) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Expected poltype=%d (Stokes parameters), got %d in template", POLTYPE_STOKES, template->poltype);
      return 0;
    }
    if (datafile->feedtype != FEEDTYPE_LINEAR && datafile->feedtype != FEEDTYPE_CIRCULAR) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Don't know how to process datafile with feedtype=%d.", datafile->feedtype);
      return 0;
    }
    if(datafile->freqMode != FREQMODE_UNIFORM) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Expected frequency channels to be equally separated.");
      return 0;
    }
    if(template->freqMode != FREQMODE_UNIFORM) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Expected frequency channels to be equally separated.");
      return 0;
    }
    /* The template is already assumed to be Stokes, so it doesn't matter if the feed type is set. BUT HOW RECEIVER PARAMETERS ARE APPLIED MAKES A DIFFERENCE!!!
       if (template->feedtype != FEEDTYPE_LINEAR && template->feedtype != FEEDTYPE_CIRCULAR) {
      fflush(stdout);
       printerror(verbose.debug, "ERROR cal_mtm: Don't know how to process template with fd_type=%d.", template->fd_type);
       return 0;
       }*/
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      if(datafile->feedtype == FEEDTYPE_LINEAR || datafile->feedtype == FEEDTYPE_INV_LINEAR)
	printf("  Using a linear basis for datafile\n");
      else if(datafile->feedtype == FEEDTYPE_CIRCULAR || datafile->feedtype == FEEDTYPE_INV_CIRCULAR)
	printf("  Using a circular basis for datafile\n");
    }
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      if(template->feedtype == FEEDTYPE_LINEAR || template->feedtype == FEEDTYPE_INV_LINEAR)
	printf("  Using a linear basis for template\n");
      else if(template->feedtype == FEEDTYPE_CIRCULAR || template->feedtype == FEEDTYPE_INV_LINEAR)
	printf("  Using a circular basis for template\n");
    }
    if(datafile->NrFreqChan != template->NrFreqChan && template->NrFreqChan != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: The number of frequency channels are different in data file and template (%ld != %ld). Note that the nr of frequency channels in principle is allowed to be 1.", datafile->NrFreqChan, template->NrFreqChan);
      return 0;
    }
    if(datafile->NrFreqChan != template->NrFreqChan && template->NrFreqChan == 1) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING cal_mtm: The template does not have frequency resolution. Although this is allowed, it is probably not smart.");    
    }
    if(datafile->poltype == POLTYPE_ILVPAdPA || datafile->poltype == POLTYPE_PAdPA || datafile->poltype == POLTYPE_ILVPAdPATEldEl) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Cannot handle PA data for datafile.");
      return 0;
    }
    if(template->poltype == POLTYPE_ILVPAdPA || template->poltype == POLTYPE_PAdPA || datafile->poltype == POLTYPE_ILVPAdPATEldEl) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Cannot handle PA data for template.");
      return 0;
    }
    if(initialguess != NULL) {
      if(initialguess->format != MEMORY_format) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: initial guess for receiver model should be read into memory first");
	return 0;
      }
      if(initialguess->gentype != GENTYPE_RECEIVERMODEL && initialguess->gentype != GENTYPE_RECEIVERMODEL2) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: initial guess for receiver model doesn't appear to be a receiver model (gentype=%s)", returnGenType_str(initialguess->gentype));
	return 0;
      }
      if(datafile->NrFreqChan != initialguess->NrFreqChan) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: The number of frequency channels are different in the data file and the initial guess receiver model.");    
	return 0;
      }
    }

    // Set the reference frequency of both the data and template to infinite frequency.
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Set reference frequency of both data and template to infinite frequency\n");
      if(preprocess_changeRefFreq(datafile, -1.0, verbose2) == 0) {
	printerror(verbose.debug, "ERROR cal_mtm: Changing reference frequency of observation failed");    
	return 0;
      }
      if(preprocess_changeRefFreq(template, -1.0, verbose2) == 0) {
	printerror(verbose.debug, "ERROR cal_mtm: Changing reference frequency of template failed");    
	return 0;
      }
    }


    /* Dedisperse both the datafile and template. If the template and
       data contain the same number of frequency channels etc this
       should strictly speaking no necessary, but it doesn't hurt. We
       now ensured that all frequency channels are aligned. */
    if(fabs(datafile->dm - template->dm) > 1e-6) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING cal_mtm: The template and data file appear to have different DM's (%f != %f). In principle this is allowed.", datafile->dm, template->dm);    
    }
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Dedispersing observation\n");
    }
    if(preprocess_dedisperse(datafile, 0, 0, 0, verbose2) == 0)
      return 0;
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Dedispersing frequency channels of the template\n");
    }
    if(preprocess_dedisperse(template, 0, 0, 0, verbose2) == 0)
      return 0;

    /* To correct for parallactic angle effect, the data shouldn't be
       corrected here, because this correction should be applied after
       correcting for receiver effects. Instead, apply opposite effect
       to the template. */

    if(datafile->isDePar == 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: The data file appears to be parallactic angle corrected. This is not allowed.");    
      return 0;
    }else if(datafile->isDePar == -1) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING cal_mtm: Parallactic angle correction state of the data is unknown. It is assumed NOT to be corrected for.");    
    }

    if(ignoreParallacticAngle == 0) {
      //First correct template
      if(verbose.verbose) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("  Remove parallactic angle effect from template\n");
      }
      if(preprocess_corrParAng(template, NULL, 0, verbose2) == 0)
	return 0;
      //Then undo effect using source and location parameters in data file
      //      template->telescope_long = datafile->telescope_long;
      //      template->telescope_lat = datafile->telescope_lat;
      template->telescope_X = datafile->telescope_X;
      template->telescope_Y = datafile->telescope_Y;
      template->telescope_Z = datafile->telescope_Z;
      template->ra = datafile->ra;
      template->dec = datafile->dec;
      template->mjd_start = datafile->mjd_start;
      template->tsubMode = datafile->tsubMode;
      //      template->tsub_list = datafile->tsub_list;
      if(template->tsub_list != NULL)
	free(template->tsub_list);
      template->tsub_list = malloc(datafile->NrSubints*sizeof(double));
      if(template->tsub_list == NULL) {
	printerror(verbose.debug, "ERROR cal_mtm: Memory allocation error.");
	return 0;
      }
      memcpy(template->tsub_list, datafile->tsub_list, datafile->NrSubints*sizeof(double));
      if(template->NrSubints != 1 || datafile->NrSubints != 1) {
	printerror(verbose.debug, "ERROR cal_mtm: At the moment multiple subints in either the template or the data file is not handled.");
	if(verbose.debug) {
	  printf("WANT TO CHANGE BOTH THE START MJD AND THE SUBINT DURATION TO SOMETHING SENSIBLE FOR EACH DATAFILE SUBINT TO BE PROCESSED\n");
	}
	exit(0);
      }
      // If correcting of template failed, at least try to apply the
      // opposite rotation, which is useful if we use the same template on
      // multiple observations.
      template->isDePar = 1;
      if(verbose.verbose) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("  Make parallactic angle of template identical to data\n");
      }
      if(preprocess_corrParAng(template, NULL, 1, verbose2) == 0)
	return 0;
    }

    /* Take out Faraday rotation in the template. Later in the code
       the template will be Faraday rotated back w.r.t. the data
       frequency and RM. */
    if(fabs(datafile->rm - template->rm) > 1e-6) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING cal_mtm: The template and data file appear to have different RM's (%f != %f). In principle this is allowed.", datafile->rm, template->rm);    
    }
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  De-Faraday rotating template\n");
    }

    if(template->isDeFarad == -1) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING cal_mtm: The de-Faraday rotation state of the template is unknown. It is assumed to be NOT de-Faraday rotated.");
      template->isDeFarad = 0;
    }else {
      if(fabs(template->rm) < 1e-6) {
	fflush(stdout);
	printwarning(verbose.debug, "WARNING cal_mtm: The template does not appear to have a defined RM. This sounds like a problem that you would like to fix before using the MTM method.");    
      }
      if(!preprocess_deFaraday(template, 0, 0, 0, NULL, verbose2))
	return 0;
    }
    // Always allow the Faraday rotation back w.r.t. the data
    // frequency and RM.
    template->isDeFarad = 1;
    template->rm = datafile->rm;


    // The datafile should still have Faraday rotation in there. Faraday rotation will be put back into template to make it match the data.
    if(datafile->isDeFarad == -1) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING cal_mtm: The de-Faraday rotation state of the data is unknown. It is assumed to be NOT de-Faraday rotated.");
      datafile->isDeFarad = 0;
    }
    if(fabs(datafile->rm) < 1e-6) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING cal_mtm: The datafile does not appear to have a defined RM. This sounds like a problem that you would like to fix before using the MTM method.");
    }
    if(!preprocess_deFaraday(datafile, 1, 0, 0, NULL, verbose2))
      return 0;


    nrbinsNew_template = template->NrBins;
    nrbinsNew_datafile = datafile->NrBins;

    if(template->NrBins != datafile->NrBins) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING cal_mtm: template and observation have different amount of bins.");
      if(template->NrBins > datafile->NrBins) {
	nrbinsNew_template = datafile->NrBins;
      }else {
	nrbinsNew_datafile = template->NrBins;
      }
    }
    i = log10(nrbinsNew_template)/log10(2);
    n = pow(2.0,(i));
    ok = 0;
    if(pow(2.0, i+1) == nrbinsNew_template)
      ok = 1;
    if(n != nrbinsNew_template && ok == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING cal_mtm: Rebinning to a power of 2 bins (for cross-correlation).");
      nrbinsNew_template = n;
      nrbinsNew_datafile = n;
    }
    if(template->NrBins != nrbinsNew_template) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING cal_mtm: Rebinning template to %ld bins.", nrbinsNew_template);
      if(!preprocess_rebin(*template, &clone, nrbinsNew_template, verbose2))
	return 0;
      swap_orig_clone(template, &clone, verbose); 
    }
    if(datafile->NrBins != nrbinsNew_datafile) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING cal_mtm: Rebinning observation to %ld bins.", nrbinsNew_datafile);
      if(!preprocess_rebin(*datafile, &clone, template->NrBins, verbose2))
	return 0;
      swap_orig_clone(datafile, &clone, verbose); 
    }
 
    /* Make a frequency scrunched clone of the template, unless we read in dump file */
    if(!(sequencenr > 0 && totsequencenr > 0)) {
      if(verbose.verbose) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("  Making frequency scrunched template and data profile for aligning the data and template\n");
      }
      cleanPSRData(&template_fscr, verbose);
      if(copy_params_PSRData(*template, &template_fscr, verbose) == 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: Error making copy of data parameters.");
	return 0;
      }
      template_fscr.data = (float *)malloc((template_fscr.NrBins)*(template_fscr.NrPols)*(template_fscr.NrFreqChan)*(template_fscr.NrSubints)*sizeof(float));
      if(template_fscr.data == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: Memory allocation error.");
	return 0;
      }
      memcpy(template_fscr.data, template->data, (template_fscr.NrBins)*(template_fscr.NrPols)*(template_fscr.NrFreqChan)*(template_fscr.NrSubints)*sizeof(float));
      if(verbose.verbose) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("    Making template profile (for aligning data and template)\n");
      }
      if(!preprocess_addsuccessiveFreqChans(template_fscr, &clone, template_fscr.NrFreqChan, NULL, verbose3))
	return 0;
      swap_orig_clone(&template_fscr, &clone, verbose); 
      
      if(verbose.verbose) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("    Making observation profile (for aligning data and template)\n");
      }
      if(!preprocess_addsuccessiveFreqChans(*datafile, &datafile_fscr, datafile->NrFreqChan, NULL, verbose3))
	return 0;
      if(verbose.verbose) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("    done\n");
      }

      /* Show template allowing user to define what the noise region is. */
      clearPulselongitudeRegion(&onpulseregion_tmp);
      pgplot_clear_options(pgplot_options);
      if(onpulseregion != NULL) {
	copyPulselongitudeRegion(*onpulseregion, &onpulseregion_tmp);
	//	if(onpulseregion->nrRegions > 0) {
	//	  memcpy(&onpulseregion_tmp, onpulseregion, sizeof(regions_definition));
	//	}
      }
      /* If onpulseregion is defined, show it, otherwise ask user for input. */
      if(onpulseregion_tmp.nrRegions > 0) {
	/* Convert onpulse region as a fraction to bins (after rebinning is done) is fractions are defined. */
	region_frac_to_int(&(onpulseregion_tmp), datafile->NrBins, 0);
	sprintf(pgplot_options->box.xlabel, "Bin");
	sprintf(pgplot_options->box.ylabel, "Intensity");
	sprintf(pgplot_options->box.title, "Excluded region(s) in noise calculation");
#ifdef HAVEPGPLOT2 
	sprintf(pgplot_options->viewport.plotDevice, "/pw");
#else
	sprintf(pgplot_options->viewport.plotDevice, "/xs");
#endif
	if(pgplotGraph1(pgplot_options, template_fscr.data, NULL, NULL, template_fscr.NrBins, 0, template_fscr.NrBins-1, 0, 0, template_fscr.NrBins-1, 0, 0, 0, 1, 0, 1, 0, 1, 1, &onpulseregion_tmp, -1, verbose) == 0) {
	  fflush(stdout);
	  printwarning(verbose.debug, "cal_mtm: WARNING plotting template failed, but continuing processing.");
	  //	  return 0;
	}
      }else {
	sprintf(pgplot_options->viewport.plotDevice, "?");
	strcpy(pgplot_options->box.xlabel, "Bin");
	strcpy(pgplot_options->box.ylabel, "Intensity");
	strcpy(pgplot_options->box.title, "Choose onpulse region(s) in the template to exclude from noise calculation");
	while(selectRegions(template_fscr.data, template_fscr.NrBins, pgplot_options, 0, 0, 0, &onpulseregion_tmp, verbose) == 0) {
	}
	if(onpulseregion != NULL) {
	  copyPulselongitudeRegion(onpulseregion_tmp, onpulseregion);
	  //	  memcpy(onpulseregion, &onpulseregion_tmp, sizeof(regions_definition));
	}
      }      
      region_int_to_frac(&(onpulseregion_tmp), 1.0/(float)datafile->NrBins, 0);
      regionShowNextTimeUse(onpulseregion_tmp, "-onpulse", "-onpulsef", stdout);
      
      /* Similar to above, but now for the high S/N part of profile. */
      clearPulselongitudeRegion(&onpulseregion_tmp2);
      if(onpulseregion_fit != NULL) {
	if(onpulseregion_fit->nrRegions > 0) {
	  copyPulselongitudeRegion(*onpulseregion_fit, &onpulseregion_tmp2);
	  //	  memcpy(&onpulseregion_tmp2, onpulseregion_fit, sizeof(regions_definition));
	}
      }
      if(onpulseregion_tmp2.nrRegions > 0) {
	/* Convert onpulse region as a fraction to bins (after rebinning is done) is fractions are defined. */
	region_frac_to_int(&(onpulseregion_tmp2), datafile->NrBins, 0);
	sprintf(pgplot_options->box.xlabel, "Bin");
	sprintf(pgplot_options->box.ylabel, "Intensity");
	sprintf(pgplot_options->box.title, "High S/N bins used in fit");
#ifdef HAVEPGPLOT2 
	sprintf(pgplot_options->viewport.plotDevice, "/pw");
#else
	sprintf(pgplot_options->viewport.plotDevice, "/xs");
#endif
	if(pgplotGraph1(pgplot_options, template_fscr.data, NULL, NULL, template_fscr.NrBins, 0, template_fscr.NrBins-1, 0, 0, template_fscr.NrBins-1, 0, 0, 0, 1, 0, 1, 0, 1, 1, &onpulseregion_tmp2, -1, verbose) == 0) {
	  fflush(stdout);
	  printerror(verbose.debug, "cal_mtm: ERROR plotting template");
	  //	  return 0;
	}
      }else {
	sprintf(pgplot_options->viewport.plotDevice, "?");
	strcpy(pgplot_options->box.xlabel, "Bin");
	strcpy(pgplot_options->box.ylabel, "Intensity");
	strcpy(pgplot_options->box.title, "Choose high S/N bins in the template to use in fit");
	while(selectRegions(template_fscr.data, template_fscr.NrBins, pgplot_options, 0, 0, 0, &onpulseregion_tmp2, verbose) == 0) {
	}
	if(onpulseregion_fit != NULL) {
	  copyPulselongitudeRegion(onpulseregion_tmp2, onpulseregion_fit);
	  //	  memcpy(onpulseregion_fit, &onpulseregion_tmp2, sizeof(regions_definition));
	}
      }
      region_int_to_frac(&(onpulseregion_tmp2), 1.0/(float)datafile->NrBins, 0);
      regionShowNextTimeUse(onpulseregion_tmp2, "-onpulse2", "-onpulsef2", stdout);
    
      /* Do cross-correlation to determine offset in pulse longitude
	 between template and data, but first debase and normalise
	 profiles to reduce numerical noise. */
      if(preprocess_debase(&datafile_fscr, &onpulseregion_tmp, NULL, 0, noverbose) == 0)
	return 0;
      if(preprocess_debase(&template_fscr, &onpulseregion_tmp, NULL, 0, noverbose) == 0)
	return 0;
      if(preprocess_norm(datafile_fscr, 1, NULL, 0, noverbose) == 0)
	return 0;
      if(preprocess_norm(template_fscr, 1, NULL, 0, noverbose) == 0)
	return 0;
      if(verbose.verbose) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("  Determining phase delay between data and template\n");
      }
      if(delay_Fourier_phase_gradient(datafile_fscr.data, template_fscr.data, datafile_fscr.NrBins, &delay, verbose2) != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: Finding delay failed");
	return 0;
      }
      //      printf("  Found a delay of %f in phase between template and observation\n", delay);

      /*
      if(verbose.debug) {
	cal_internal_drawTwoProfiles(cc_function, cc_function, cclength, "Cross-correlation function. Data (white) and template (red)");
	printf("Press a key in terminal to continue\n");
	pgetch();
      }
      */
      if(verbose.verbose) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("  Rotate data to match template\n");
      }
      if(aligned == 0) {
	/* Rotate data to match template */
	if(preprocess_fftshift(*datafile, delay, 0, 0, verbose2) == 0)
	  return 0;
      }else {
	printf("  Ignoring delay, assuming profiles are already aligned.\n");
      }
      if(verbose.debug) {
	printf("  Showing aligned profile\n");
	/* Now the profiles are aligned, debase and normalise again as
	   onpulse region should now match the aligned data. */
	if(preprocess_fftshift(datafile_fscr, delay, 0, 0, verbose2) == 0)
	  return 0;
	if(preprocess_debase(&datafile_fscr, &onpulseregion_tmp, NULL, 0, noverbose) == 0)
	  return 0;
	if(preprocess_norm(datafile_fscr, 1, NULL, 0, noverbose_nocounters) == 0)
	  return 0;
	cal_internal_drawTwoProfiles(datafile_fscr.data, template_fscr.data, datafile_fscr.NrBins, "Aligned profiles. Data (white) and template (red)");
	printf("Press a key in terminal to continue\n");
	pgetch();
      }
      closePSRData(&datafile_fscr, 0, verbose);
      closePSRData(&template_fscr, 0, verbose);
    }
  }

  /* Do reading/writing of dumpfile */
  /* Start writing dumpfile */
  if(sequencenr < 0 && totsequencenr > 0) {
    datafile->fptr_hdr = dumpfile;
    if(writePSRSALSAHeader(datafile, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
      return 0;
    }
    ret = fwrite(datafile->data, 1, (datafile->NrBins)*(datafile->NrPols)*(datafile->NrFreqChan)*(datafile->NrSubints)*sizeof(float), dumpfile);
    if(ret != (datafile->NrBins)*(datafile->NrPols)*(datafile->NrFreqChan)*(datafile->NrSubints)*sizeof(float)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
      return 0;
    }
  
    template->fptr_hdr = dumpfile;
    if(writePSRSALSAHeader(template, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
      return 0;
    }
    ret = fwrite(template->data, 1, (template->NrBins)*(template->NrPols)*(template->NrFreqChan)*(template->NrSubints)*sizeof(float), dumpfile);
    if(ret != (template->NrBins)*(template->NrPols)*(template->NrFreqChan)*(template->NrSubints)*sizeof(float)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
      return 0;
    }
    ret = fwrite(&(onpulseregion_tmp.nrRegions), sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
      return 0;
    }
    for(i = 0; i < onpulseregion_tmp.nrRegions; i++) {
      ret = fwrite(&(onpulseregion_tmp.left_bin[i]), sizeof(int), 1, dumpfile);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
	return 0;
      }
      ret = fwrite(&(onpulseregion_tmp.right_bin[i]), sizeof(int), 1, dumpfile);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
	return 0;
      }
    }
    ret = fwrite(&(onpulseregion_tmp2.nrRegions), sizeof(int), 1, dumpfile);
    for(i = 0; i < onpulseregion_tmp2.nrRegions; i++) {
      ret = fwrite(&(onpulseregion_tmp2.left_bin[i]), sizeof(int), 1, dumpfile);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
	return 0;
      }
      ret = fwrite(&(onpulseregion_tmp2.right_bin[i]), sizeof(int), 1, dumpfile);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: Write error to %s", dumpfilename);
	return 0;
      }
    }
    //      destination->bins_defined[i] = source.bins_defined[i];

    fclose(dumpfile);

    freePulselongitudeRegion(&onpulseregion_tmp);
    freePulselongitudeRegion(&onpulseregion_tmp2);
    free(pgplot_options);

    return 1;
  }

  /* Start reading dumpfile */
  if(sequencenr > 0 && totsequencenr > 0) {
    cleanPSRData(datafile, verbose);
    datafile->fptr_hdr = dumpfile;
    if(readPSRSALSAHeader(datafile, 1, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s", dumpfilename);
      return 0;
    }
    datafile->format = MEMORY_format;
    datafile->data = (float *)malloc((datafile->NrBins)*(datafile->NrPols)*(datafile->NrFreqChan)*(datafile->NrSubints)*sizeof(float));
    if(datafile->data == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Memory allocation error");
      return 0;
    }
    ret = fread(datafile->data, 1, (datafile->NrBins)*(datafile->NrPols)*(datafile->NrFreqChan)*(datafile->NrSubints)*sizeof(float), dumpfile);
    if(ret != (datafile->NrBins)*(datafile->NrPols)*(datafile->NrFreqChan)*(datafile->NrSubints)*sizeof(float)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s", dumpfilename);
      return 0;
    }

    cleanPSRData(template, verbose);
    template->fptr_hdr = dumpfile;
    if(readPSRSALSAHeader(template, 1, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s", dumpfilename);
      return 0;
    }
    template->format = MEMORY_format;
    template->data = (float *)malloc((template->NrBins)*(template->NrPols)*(template->NrFreqChan)*(template->NrSubints)*sizeof(float));
    if(template->data == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Memory allocation error");
      return 0;
    }
    ret = fread(template->data, 1, (template->NrBins)*(template->NrPols)*(template->NrFreqChan)*(template->NrSubints)*sizeof(float), dumpfile);
    if(ret != (template->NrBins)*(template->NrPols)*(template->NrFreqChan)*(template->NrSubints)*sizeof(float)) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error from %s", dumpfilename);
      return 0;
    }
    ret = fread(&(onpulseregion_tmp.nrRegions), sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error to %s", dumpfilename);
      return 0;
    }
    for(i = 0; i < onpulseregion_tmp.nrRegions; i++) {
      ret = fread(&(onpulseregion_tmp.left_bin[i]), sizeof(int), 1, dumpfile);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: Read error to %s", dumpfilename);
	return 0;
      }
      ret = fread(&(onpulseregion_tmp.right_bin[i]), sizeof(int), 1, dumpfile);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: Read error to %s", dumpfilename);
	return 0;
      }
      onpulseregion_tmp.bins_defined[i] = 1;
      onpulseregion_tmp.frac_defined[i] = 0;
    }
    ret = fread(&(onpulseregion_tmp2.nrRegions), sizeof(int), 1, dumpfile);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Read error to %s", dumpfilename);
      return 0;
    }
    for(i = 0; i < onpulseregion_tmp2.nrRegions; i++) {
      ret = fread(&(onpulseregion_tmp2.left_bin[i]), sizeof(int), 1, dumpfile);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: Read error to %s", dumpfilename);
	return 0;
      }
      ret = fread(&(onpulseregion_tmp2.right_bin[i]), sizeof(int), 1, dumpfile);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR cal_mtm: Read error to %s", dumpfilename);
	return 0;
      }
      onpulseregion_tmp2.bins_defined[i] = 1;
      onpulseregion_tmp2.frac_defined[i] = 0;
    }
  }

  /* Remove baseline from data and template */
  if(preprocess_debase(datafile, &onpulseregion_tmp, NULL, 0, noverbose) == 0)
    return 0;
  if(preprocess_debase(template, &onpulseregion_tmp, NULL, 0, noverbose) == 0)
    return 0;


  /* Get the number of bins in the high-S/N selected range and store
     the value in internal_nfree for now. Will be converted into the
     number of degrees of freedom in a moment.  */
  internal_nfree = 0;
  for(b = 0; b < datafile->NrBins; b++) {
    if(checkRegions(b, &onpulseregion_tmp2, 0, verbose) != 0)
      internal_nfree++;
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  Data will be fit over %d on-pulse bins\n", internal_nfree);
  }

  nractivebins_offpulse = 0;
  for(b = 0; b < datafile->NrBins; b++) {
    if(checkRegions(b, &onpulseregion_tmp, 0, verbose) == 0)
      nractivebins_offpulse++;
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  There are %ld off-pulse bins for noise calculation\n", nractivebins_offpulse);
  }

  /* Allocate data to hold one frequency channel and one data channel
     and a calibrated channel. Only the high S/N part is stored as
     that is the only relevant bit for the fitting process. */
  cleanPSRData(&internal_data_channel_onpulse, verbose);
  if(copy_params_PSRData(*datafile, &internal_data_channel_onpulse, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm:  Making copy of data parameters failed.");
    return 0;
  }
  internal_data_channel_onpulse.NrBins = internal_nfree;
  internal_data_channel_onpulse.NrFreqChan = 1;
  internal_data_channel_onpulse.NrSubints = 1;
  internal_data_channel_onpulse.data = (float *)malloc((internal_data_channel_onpulse.NrBins)*(internal_data_channel_onpulse.NrPols)*(internal_data_channel_onpulse.NrFreqChan)*(internal_data_channel_onpulse.NrSubints)*sizeof(float));
  /* printf("XXXX allocating %ld bytes (channel) \n", (internal_data_channel_onpulse.NrBins)*(internal_data_channel_onpulse.NrPols)*(internal_data_channel_onpulse.NrFreqChan)*(internal_data_channel_onpulse.NrSubints)*sizeof(float)); */

  cleanPSRData(&internal_data_channel_onpulse_calib, verbose);
  if(copy_params_PSRData(internal_data_channel_onpulse, &internal_data_channel_onpulse_calib, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm:  Making copy of data parameters failed.");
    return 0;
  }
  internal_data_channel_onpulse_calib.data = (float *)malloc((internal_data_channel_onpulse_calib.NrBins)*(internal_data_channel_onpulse_calib.NrPols)*(internal_data_channel_onpulse_calib.NrFreqChan)*(internal_data_channel_onpulse_calib.NrSubints)*sizeof(float));
  /*  printf("XXXX allocating %ld bytes (channel) \n", (internal_data_channel_onpulse_calib.NrBins)*(internal_data_channel_onpulse_calib.NrPols)*(internal_data_channel_onpulse_calib.NrFreqChan)*(internal_data_channel_onpulse_calib.NrSubints)*sizeof(float)); */

  cleanPSRData(&data_channel_calib_best, verbose);
  if(copy_params_PSRData(internal_data_channel_onpulse, &data_channel_calib_best, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm:  Making copy of data parameters failed.");
    return 0;
  }
  data_channel_calib_best.data = (float *)malloc((data_channel_calib_best.NrBins)*(data_channel_calib_best.NrPols)*(data_channel_calib_best.NrFreqChan)*(data_channel_calib_best.NrSubints)*sizeof(float));
  /*  printf("XXXX allocating %ld bytes (channel) \n", (data_channel_calib_best.NrBins)*(data_channel_calib_best.NrPols)*(data_channel_calib_best.NrFreqChan)*(data_channel_calib_best.NrSubints)*sizeof(float)); */


  cleanPSRData(&internal_template_channel_onpulse, verbose);
  if(copy_params_PSRData(*template, &internal_template_channel_onpulse, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm:  Making copy of data parameters failed.");
    return 0;
  }
  internal_template_channel_onpulse.NrBins = internal_nfree;
  internal_template_channel_onpulse.NrFreqChan = 1;
  internal_template_channel_onpulse.NrSubints = 1;
  internal_template_channel_onpulse.data = (float *)malloc((internal_template_channel_onpulse.NrBins)*(internal_template_channel_onpulse.NrPols)*(internal_template_channel_onpulse.NrFreqChan)*(internal_template_channel_onpulse.NrSubints)*sizeof(float));
  /*  printf("XXXX allocating %ld bytes (channel) \n", (internal_template_channel_onpulse.NrBins)*(internal_template_channel_onpulse.NrPols)*(internal_template_channel_onpulse.NrFreqChan)*(internal_template_channel_onpulse.NrSubints)*sizeof(float)); */

  cleanPSRData(&internal_template_channel_onpulse_calib, verbose);
  if(copy_params_PSRData(internal_template_channel_onpulse, &internal_template_channel_onpulse_calib, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm:  Making copy of data parameters failed.");
    return 0;
  }
  // Set the polarization type (lin/cir) to be equal to that of the data. Because we normally fit the template to the data we should find the Mueller matrix in the basis of the data.
  internal_template_channel_onpulse_calib.feedtype = datafile->feedtype;
  internal_template_channel_onpulse_calib.data = (float *)malloc((internal_template_channel_onpulse_calib.NrBins)*(internal_template_channel_onpulse_calib.NrPols)*(internal_template_channel_onpulse_calib.NrFreqChan)*(internal_template_channel_onpulse_calib.NrSubints)*sizeof(float));
  /*  printf("XXXX allocating %ld bytes (channel) \n", (internal_template_channel_onpulse_calib.NrBins)*(internal_template_channel_onpulse_calib.NrPols)*(internal_template_channel_onpulse_calib.NrFreqChan)*(internal_template_channel_onpulse_calib.NrSubints)*sizeof(float)); */

  cleanPSRData(&data_channel_full, verbose);
  if(copy_params_PSRData(*datafile, &data_channel_full, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm:  Making copy of data parameters failed.");
    return 0;
  }

  data_channel_full.NrFreqChan = 1;
  data_channel_full.NrSubints = 1;
  data_channel_full.data = (float *)malloc((data_channel_full.NrBins)*(data_channel_full.NrPols)*(data_channel_full.NrFreqChan)*(data_channel_full.NrSubints)*sizeof(float));

  /* For verbose purposes: Allocate data to hold frequency scrunched
     result. Only the high S/N part is stored as that is the only
     relevant bit for the fitting process. */
  cleanPSRData(&fscrunched_calib_obs, verbose);
  if(copy_params_PSRData(*datafile, &fscrunched_calib_obs, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm:  Making copy of data parameters failed.");
    return 0;
  }
  fscrunched_calib_obs.NrBins = internal_nfree;
  fscrunched_calib_obs.NrFreqChan = 1;
  fscrunched_calib_obs.NrSubints = 1;
  fscrunched_calib_obs.data = (float *)malloc((fscrunched_calib_obs.NrBins)*(fscrunched_calib_obs.NrPols)*(fscrunched_calib_obs.NrFreqChan)*(fscrunched_calib_obs.NrSubints)*sizeof(float));

  /* For verbose purposes: Allocate data to hold frequency scrunched
     result. Only the high S/N part is stored as that is the only
     relevant bit for the fitting process. */
  cleanPSRData(&fscrunched_template, verbose);
  if(copy_params_PSRData(*datafile, &fscrunched_template, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm:  Making copy of data parameters failed.");
    return 0;
  }
  fscrunched_template.NrBins = internal_nfree;
  fscrunched_template.NrFreqChan = 1;
  fscrunched_template.NrSubints = 1;
  fscrunched_template.data = (float *)malloc((fscrunched_template.NrBins)*(fscrunched_template.NrPols)*(fscrunched_template.NrFreqChan)*(fscrunched_template.NrSubints)*sizeof(float));

  // Clear memory for frequency scrunched end result
  if(verbose.verbose) {
    long i1, i2, i3;
    for(i1 = 0; i1 < 4; i1++) {
      for(i2 = 0; i2 < fscrunched_calib_obs.NrBins; i2++) {
	i3 = i1*fscrunched_calib_obs.NrBins+i2;
	fscrunched_calib_obs.data[i3] = 0;
	fscrunched_template.data[i3] = 0;
      }
    }
  }

  /* Also allocate memory to hold a single data and template channel,
     which will be used to calculate rms */
  template_channel_full = (float *)malloc((datafile->NrBins)*(datafile->NrPols)*sizeof(float));

  if(internal_data_channel_onpulse.data == NULL || internal_data_channel_onpulse_calib.data == NULL || internal_template_channel_onpulse_calib.data == NULL || data_channel_calib_best.data == NULL || internal_template_channel_onpulse.data == NULL || data_channel_full.data == NULL || template_channel_full == NULL || fscrunched_calib_obs.data == NULL || fscrunched_template.data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm: Memory allocation error.");
    return 0;
  }

  /* Calculate the degrees of freedom in fitting process */
  internal_leakage_fit_flag = leakage;
  internal_nfree = 4*internal_nfree;
  if(leakage)
    internal_nfree -= 7;
  else
    internal_nfree -= 3;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  There are %d degrees of freedom\n", internal_nfree);
  }

  /* Prepare the receiver solution file */
  cleanPSRData(&solutionall, verbose);
  if(copy_params_PSRData(*datafile, &solutionall, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm:  Making copy of data parameters failed.");
    return 0;
  }
  /* Data contains npol parameters. If nbin=2 then there is also an errorbar */
  if(leakage == 0) {
    solutionall.NrPols = 3;
    solutionall.gentype = GENTYPE_RECEIVERMODEL;  /* set to GENTYPE_RECEIVERMODEL2 if storing nfree and chi^2 */
  }else {
    solutionall.NrPols = 9;
    solutionall.gentype = GENTYPE_RECEIVERMODEL2;
  }
  solutionall.NrBins = 1;
  solutionall.poltype = POLTYPE_UNKNOWN;
  /*  printf("XXXXX allocating %ld bytes\n", (solutionall.NrBins)*(solutionall.NrPols)*(solutionall.NrFreqChan)*(solutionall.NrSubints)*sizeof(float)); */
  solutionall.data = (float *)malloc((solutionall.NrBins)*(solutionall.NrPols)*(solutionall.NrFreqChan)*(solutionall.NrSubints)*sizeof(float));
  /*  printf("XXXX allocating %ld bytes (total solution) \n", (solutionall.NrBins)*(solutionall.NrPols)*(solutionall.NrFreqChan)*(solutionall.NrSubints)*sizeof(float)); */

  cleanPSRData(&internal_solution_cur_channel, verbose);
  if(copy_params_PSRData(solutionall, &internal_solution_cur_channel, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm:  Making copy of data parameters failed.");
    return 0;
  }
  internal_solution_cur_channel.NrFreqChan = 1;
  internal_solution_cur_channel.NrSubints = 1;
  internal_solution_cur_channel.data = (float *)malloc((internal_solution_cur_channel.NrBins)*(internal_solution_cur_channel.NrPols)*(internal_solution_cur_channel.NrFreqChan)*(internal_solution_cur_channel.NrSubints)*sizeof(float));
  if(internal_solution_cur_channel.data == NULL || solutionall.data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm: Memory allocation error.");
    return 0;
  }
  /* printf("XXXX allocating %ld bytes (single solution) \n", (internal_solution_cur_channel.NrBins)*(internal_solution_cur_channel.NrPols)*(internal_solution_cur_channel.NrFreqChan)*(internal_solution_cur_channel.NrSubints)*sizeof(float)); */

  /* Allocate memory to store rms values */
  internal_data_channel_onpulse_calib.offpulse_rms = (float *)malloc(internal_data_channel_onpulse.NrPols*sizeof(float));
  //  internal_data_channel_onpulse_calib.rmspoints_defined = 1;
  internal_template_channel_onpulse.offpulse_rms = (float *)malloc(internal_template_channel_onpulse.NrPols*sizeof(float));
  //  internal_template_channel_onpulse.rmspoints_defined = 1;
  if(internal_data_channel_onpulse_calib.offpulse_rms == NULL || internal_template_channel_onpulse.offpulse_rms == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm: Memory allocation error.");
    return 0;
  }
  /*  printerror(verbose.debug, "XXXX %ld %ld", internal_data_channel_onpulse.NrPols*sizeof(float), internal_template_channel_onpulse.NrPols*sizeof(float)); */

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  Start fitting process\n");
  }

  /* Loop over pulses and frequency channels etc and start the hard work. */
  for(n = 0; n < datafile->NrSubints; n++) {
    for(f = 0; f < datafile->NrFreqChan; f++) {
      /*
      if(f == 160)
	debug_flag = 1;
      else
	debug_flag = 0;
      */
      int oktoprocesschannel;
      /* Determine if this particular channel needs to be processed */
      oktoprocesschannel = 1;
      if(sequencenr > 0 && totsequencenr > 0) {
	oktoprocesschannel = 0;
	if(f%totsequencenr == sequencenr-1)
	  oktoprocesschannel = 1;
      }
      if(oktoprocesschannel) {
	/* Compare same frequency channel in data and template, unless the template only has one frequency channel. */
	ftemplate = f;
	if(template->NrFreqChan == 1)
	  ftemplate = 0;
	/* Read a channel from template and data stream */
	int curbin;
	double avrg1, avrg2, rms_temp;
	for(p = 0; p < internal_data_channel_onpulse.NrPols; p++) {
	  if(readPulsePSRData(datafile, n, p, f, 0, datafile->NrBins, &data_channel_full.data[p*data_channel_full.NrBins], verbose) == 0)
	    return 0;
	  if(readPulsePSRData(template, n, p, ftemplate, 0, datafile->NrBins, &(template_channel_full[p*template->NrBins]), verbose) == 0)
	    return 0;
	  /* Calculate the rms in the template */
	  rms_temp = 0;
	  for(b = 0; b < datafile->NrBins; b++) {
	    if(checkRegions(b, &onpulseregion_tmp, 0, verbose) == 0) {
	      /* Note: baseline is already subtracted */
	      rms_temp += (double)(template_channel_full[b+p*template->NrBins])*(double)(template_channel_full[b+p*template->NrBins]);
	    }
	  }
	  rms_temp /= (double)nractivebins_offpulse;
	  internal_template_channel_onpulse.offpulse_rms[p] = sqrt(rms_temp);
	  internal_data_channel_onpulse_calib.offpulse_rms[p] = 0;  /* Not sure if used, but old code had it set to zero */
	  if(verbose.debug) printf("  RMS polarization %ld: %f (template)\n", p+1, internal_template_channel_onpulse.offpulse_rms[p]);

	  /* Copy high-S/N part of profile to single-pulse single-freq
	     data-sets to be passed on to fitter. Also get scaling of
	     Stokes I to get initial guess for required gain*/
	  curbin = 0;
	  if(p == 0)
	    avrg1 = avrg2 = 0;
	  for(b = 0; b < datafile->NrBins; b++) {
	    if(checkRegions(b, &onpulseregion_tmp2, 0, verbose) != 0) {
	      internal_data_channel_onpulse.data[p*internal_data_channel_onpulse.NrBins+curbin] = data_channel_full.data[b+p*data_channel_full.NrBins];
	      internal_template_channel_onpulse.data[p*internal_template_channel_onpulse.NrBins+curbin] = template_channel_full[b+p*template->NrBins];
	      if(p == 0) {
		avrg1 += data_channel_full.data[b+p*data_channel_full.NrBins];
		avrg2 += template_channel_full[b+p*template->NrBins];
	      }
	      curbin++;
	    }
	  }
	}

	//Faraday rotate template back to that expected for the current data channel (relative to inf freq)
	double dphi, L, phi;
	dphi = -calcRMAngle(get_nonweighted_channel_freq(*datafile, f, verbose), 0.0, 1, datafile->rm);
	for(b = 0; b < internal_template_channel_onpulse.NrBins; b++) {
	  L = sqrt(internal_template_channel_onpulse.data[internal_template_channel_onpulse.NrBins+b]*internal_template_channel_onpulse.data[internal_template_channel_onpulse.NrBins+b]+internal_template_channel_onpulse.data[2*internal_template_channel_onpulse.NrBins+b]*internal_template_channel_onpulse.data[2*internal_template_channel_onpulse.NrBins+b]);
	  phi = atan2(internal_template_channel_onpulse.data[2*internal_template_channel_onpulse.NrBins+b], internal_template_channel_onpulse.data[internal_template_channel_onpulse.NrBins+b]);
	  phi -= 2.0*dphi;
	  internal_template_channel_onpulse.data[internal_template_channel_onpulse.NrBins+b] = L*cos(phi);
	  internal_template_channel_onpulse.data[2*internal_template_channel_onpulse.NrBins+b] = L*sin(phi);
	}

  /*  Test to see if initial estimate for gain worked
      internal_solution_cur_channel.data[0] = sqrt(fabs(avrg1/avrg2));
      internal_solution_cur_channel.data[1] = 0;
      internal_solution_cur_channel.data[2] = 0;
      if(debug) printf("  Initial guess for gain: %f\n", internal_solution_cur_channel.data[0]);
      
      if(applyReceiverModel(&internal_data_channel_onpulse, internal_solution_cur_channel, 1, 0, verbose) == 0) {
      fflush(stdout);
      printerror(debug, "ERROR cal_mtm: Cannot apply initional solution.");
      return 0;
      }*/
	
  /* Calculate the number of trials requested for each parameter */
	np1tot = (1+2*bruteforce);
	np2tot = (1+2*bruteforce);
	np3tot = (1+2*bruteforce);
	np4tot = (1+2*bruteforce);
	np5tot = (1+2*bruteforce);
	np6tot = (1+2*bruteforce);
	np7tot = (1+2*bruteforce);
	if(extremes) {
	  np3tot = 3;
	  np4tot = 3;
	  np5tot = 3;
	  np6tot = 3;
	  np7tot = 3;
	}
	if(fixed[0])
	  np1tot = 1;
	if(fixed[1])
	  np2tot = 1;
	if(fixed[2])
	  np3tot = 1;
	if(fixed[3])
	  np4tot = 1;
	if(fixed[4])
	  np5tot = 1;
	if(fixed[5])
	  np6tot = 1;
	if(fixed[6])
	  np7tot = 1;
	if(leakage == 0) {
	  np4tot = 1;
	  np5tot = 1;
	  np6tot = 1;
	  np7tot = 1;
	}
	
	double spacing;
	chisqbest = -1;
	data_fom_best = -1;
	for(np1 = 0; np1 < np1tot; np1++) {
	  if(getinitialguess_internal(initialguess, 0, f, &startx[0], verbose) == 0)
	    return 0;
	  /* Ignore gain from receiver model */
	  startx[0] = sqrt(fabs(avrg1/avrg2));
	  if(np1tot > 1) {
	    /* Let the gain vary with a factor between 0.5 to 1.5 */
	    startx[0] *= 0.5+1.0*np1/(double)(np1tot-1);
	    /* Take initial step size 0.3 times separation between grid points */
	    //	    dx[0] = 0.3*1.0/(double)(np1tot-1);
	    dx[0] = 0.1*startx[0];
	  }else {
	    dx[0] = 0.1*startx[0];
	  }
	  for(np2 = 0; np2 < np2tot; np2++) {
	    if(getinitialguess_internal(initialguess, 1, f, &startx[1], verbose) == 0)
	      return 0;
	    /* Ignore parameter from receiver model if applying brute force method */
	    if(np2tot > 1) {
	      /* Let the differential gain vary between -1 and 1 */
	      spacing = 2.0/(double)(2.0*(1.0+bruteforce));
	      startx[1] = 0 + bruteforce*spacing*(2.0*np2/(double)(np2tot-1) -1);
	      //	      printf("XXXXXX: gamma=%lf ", startx[1]);
	      /* Convert from gamma to beta */
	      startx[1] = 0.5*log(startx[1]+1.0);
	      /* Take initial step size 0.3 times separation between grid points */
	      dx[1] = 0.5*log(startx[1]+0.3*spacing+1.0)-startx[1];
	      if(isnan(dx[1])) {
		dx[1] = 0.1;
	      }
	      //	      printf("beta=%lf dx=%lf spacing=lf\n", startx[1], dx[1], spacing);
	    }else {
	      dx[1] = 0.1;
	    }
	    for(np3 = 0; np3 < np3tot; np3++) {
	      if(getinitialguess_internal(initialguess, 2, f, &startx[2], verbose) == 0)
		return 0;
	      /* Ignore parameter from receiver model if applying brute force method */
	      if(np3tot > 1) {
		/* The differential phase can vary between -pi/2 and pi/2, after which it repeats, so place points at -pi/4, 0 and pi/4. */
		if(extremes) {
		  if(np3 == 0)
		    startx[2] = 0;
		  else if(np3 == 1)
		    startx[2] = M_PI/4.0;
		  else
		    startx[2] = -M_PI/4.0;
		  dx[2] = 0.1;
		}else {
		  spacing = M_PI/(double)(2*(1+bruteforce));
		  startx[2] = 0 + bruteforce*spacing*(2.0*np3/(double)(np3tot-1) -1);
		  /* Take initial step size 0.3 times separation between grid points */
		  dx[2] = 0.3*spacing;
		}
	      }else {
		dx[2] = 0.1;
	      }
	      for(np4 = 0; np4 < np4tot; np4++) {
		if(getinitialguess_internal(initialguess, 3, f, &startx[3], verbose) == 0)
		  return 0;
		/* Ignore parameter from receiver model if applying brute force method */
		if(np4tot > 1) {
		  /* Let the ellipticities vary between -pi/4 and pi/4, so place points at -pi/8, 0 and pi/8. The 0.95 is to avoid "nice" angles that cause the Mueller matrix to contain rows of zero's. This happens for instance if the ellipticities are opposite and or0-or1=90 deg. */
		  if(extremes) {
		    if(np4 == 0)
		      startx[3] = 0;
		    else if(np4 == 1)
		      startx[3] = 0.95*M_PI/8.0;		
		    else
		      startx[3] = -0.95*M_PI/8.0;		
		    dx[3] = 0.1;
		  }else {
		    spacing = 0.95*M_PI/(double)(4.0*(1+bruteforce));
		    startx[3] = 0 + bruteforce*spacing*(2.0*np4/(double)(np4tot-1) -1);
		    /* Take initial step size 0.3 times separation between grid points */
		    dx[3] = 0.3*spacing;
		  }
		}else {
		  dx[3] = 5*M_PI/180.0;
		}
		for(np5 = 0; np5 < np5tot; np5++) {
		  if(getinitialguess_internal(initialguess, 4, f, &startx[4], verbose) == 0)
		    return 0;
		  /* Ignore parameter from receiver model if applying brute force method */
		  if(np5tot > 1) {
		    /* Let the orientations vary between -pi/2 and pi/2, so place points at -pi/4, 0 and pi/4. The 0.95 is to avoid "nice" angles that cause the Mueller matrix to contain rows of zero's. This happens for instance if the ellipticities are opposite and or0-or1=90 deg.*/
		    if(extremes) {
		      if(np5 == 0)
			startx[4] = 0;
		      else if(np5 == 1)
			startx[4] = 0.95*M_PI/4.0;		
		      else
			startx[4] = -0.95*M_PI/4.0;		
		      dx[4] = 0.1;
		    }else {
		      spacing = 0.95*M_PI/(double)(2*(1+bruteforce));
		      startx[4] = 0 + bruteforce*spacing*(2.0*np5/(double)(np5tot-1) -1);
		      /* Take initial step size 0.3 times separation between grid points */
		      dx[4] = 0.3*spacing;
		    }
		  }else {
		    dx[4] = 5*M_PI/180.0;
		  }
		  for(np6 = 0; np6 < np6tot; np6++) {
		    if(getinitialguess_internal(initialguess, 5, f, &startx[5], verbose) == 0)
		      return 0;
		    /* Ignore parameter from receiver model if applying brute force method */
		    if(np6tot > 1) {
		      /* Let the ellipticities vary between -pi/4 and pi/4, so place points at -pi/8, 0 and pi/8. The 0.95 is to avoid "nice" angles that cause the Mueller matrix to contain rows of zero's. This happens for instance if the ellipticities are opposite and or0-or1=90 deg. */
		      if(extremes) {
			if(np6 == 0)
			  startx[5] = 0;
			else if(np6 == 1)
			  startx[5] = 0.95*M_PI/8.0;		
			else
			  startx[5] = -0.95*M_PI/8.0;		
			dx[5] = 0.1;
		      }else {
			spacing = 0.95*M_PI/(double)(4.0*(1+bruteforce));
			startx[5] = 0 + bruteforce*spacing*(2.0*np6/(double)(np6tot-1) -1);
			/* Take initial step size 0.3 times separation between grid points */
			dx[5] = 0.3*spacing;
		      }
		    }else {
		      dx[5] = 5*M_PI/180.0;
		    }
		    for(np7 = 0; np7 < np7tot; np7++) {
		      if(getinitialguess_internal(initialguess, 6, f, &startx[6], verbose) == 0)
			return 0;
		      /* Ignore parameter from receiver model if applying brute force method */
		      if(np7tot > 1) {
			/* Let the orientations vary between -pi/2 and pi/2, so place points at -pi/4, 0 and pi/4. The 0.95 is to avoid "nice" angles that cause the Mueller matrix to contain rows of zero's. This happens for instance if the ellipticities are opposite and or0-or1=90 deg. */
			if(extremes) {
			  if(np7 == 0)
			    startx[6] = 0;
			  else if(np7 == 1)
			    startx[6] = 0.95*M_PI/4.0;		
			  else
			    startx[6] = -0.95*M_PI/4.0;		
			  dx[6] = 0.1;
			}else {
			  spacing = 0.95*M_PI/(double)(2*(1+bruteforce));
			  startx[6] = 0 + bruteforce*spacing*(2.0*np7/(double)(np7tot-1) -1);
			  /* Take initial step size 0.3 times separation between grid points */
			  dx[6] = 0.3*spacing;
			}
		      }else {
			dx[6] = 5*M_PI/180.0;
		      }
		      /* Now we now initial parameters, do the actual fitting. */
		      /*
		      if((verbose && f == 0) || verbose.debug) {
			if(leakage == 0) {
			  printf("Initial guess for channel n=%04ld freq=%04ld: G=%f  diff=%f   phase=%f\n", n, f, startx[0], startx[1], startx[2]);
			  printf("                                  stepsize: G=%f  diff=%f   phase=%f\n", dx[0], dx[1], dx[2]);
			}else {
			  printf("Initial guess for channel n=%04ld freq=%04ld: G=%f  diff=%f   phase=%f   %f   %f   %f   %f\n", n, f, startx[0], startx[1], startx[2], startx[3], startx[4], startx[5], startx[6]);
			  printf("                                  stepsize: G=%f  diff=%f   phase=%f   %f   %f   %f   %f\n", dx[0], dx[1], dx[2], dx[3], dx[4], dx[5], dx[6]);
			}
		      }
		      */

/* Far too big!!!! Is in radians!!!! */
/*
      dx[3] = 5;   
      dx[4] = 5;
      dx[5] = 5;
      dx[6] = 5;
*/
/*
      if(initialguess != NULL) {
	int nrpols;
	nrpols = initialguess->NrPols;
	if(initialguess->gentype == GENTYPE_RECEIVERMODEL2)
	  nrpols -= 2;
	if(nrpols > 1) {
	  for(p = 1; p < nrpols; p++) {
	    if(readPulsePSRData(initialguess, 0, p, f, 0, 1, &startx[p], verbose.debug) != 1) {
      fflush(stdout);
	      printerror(verbose.debug, "ERROR cal_mtm: Cannot read from initial guess receiver model.");
	      return 0;
	    }
	  }
	}
	}*/
      /*
      fixed[0] = 0;
      fixed[1] = 0;
      fixed[2] = 0;
      fixed[3] = 0;
      fixed[4] = 0;
      fixed[5] = 0;
      fixed[6] = 0;
      */
		      ftol = 1e-4;
		      finderrors = 0;
		      //		      if(verbose.debug) printf("Initial guess for gain: %lf\n", startx[0]);
		      iszapped = 0;
		      if(isnan(avrg1) || isnan(avrg2) || isinf(avrg1) || isinf(avrg2) || avrg1 == 0 || avrg2 == 0) {
			if(verbose.debug) printf("  Channel %ld appears to be zapped\n", f);
			fittedparams[0] = -1e-6;
			fittedparams[1] = 0;
			fittedparams[2] = 0;
			fittedparams[3] = 0;
			fittedparams[4] = 0;
			fittedparams[5] = 0;
			fittedparams[6] = 0;
			chisq = 1;
			chisq_norm = 1;
			iszapped = 1;
			/*	printerror(verbose.debug, "XXXXX startx[0] = %f", startx[0]); */
		      }
		      failed = 0;
		      
		      if(iszapped == 0) {
			//			gsl_error_handler_t *oldhandler;
			//			oldhandler = gsl_set_error_handler_off();
			if(leakage == 0) {
			  if(doAmoeba_d(0, startx, dx, fixed, fittedparams, &chisq, 3, funk_optimizeMTM, ftol, &nfunk,  verbose.debug, finderrors, 1, dplus, dmin) != 0) {
			    fflush(stdout);
			    printwarning(verbose.debug, "WARNING cal_mtm: Fit algorithm didn't converge (f=%ld sub=%ld), trying lowering tollerance.", f, n);
			    if(doAmoeba_d(0, startx, dx, fixed, fittedparams, &chisq, 3, funk_optimizeMTM, 100*ftol, &nfunk,  verbose.debug, finderrors, 1, dplus, dmin) != 0) {

			      if(bruteforce) {
				fflush(stdout);
				printwarning(verbose.debug, "WARNING cal_mtm: Fit algorithm didn't converge (f=%ld sub=%ld).", f, n);
				failed = 1;
			      }else {
				fflush(stdout);
				printerror(verbose.debug, "ERROR cal_mtm: Fit algorithm didn't converge (f=%ld sub=%ld).", f, n);
				return 0;
			      }
			    }else {
			      fflush(stdout);
			      printwarning(verbose.debug, "WARNING cal_mtm: Success!");
			    }
			  }
			}else {
			  if(doAmoeba_d(0, startx, dx, fixed, fittedparams, &chisq, 7, funk_optimizeMTM, ftol, &nfunk,  verbose.debug, finderrors, 1, dplus, dmin) != 0) {
			    fflush(stdout);
			    printwarning(verbose.debug, "WARNING cal_mtm: Fit algorithm didn't converge (f=%ld sub=%ld), trying lowering tollerance.", f, n);
			    if(doAmoeba_d(0, startx, dx, fixed, fittedparams, &chisq, 7, funk_optimizeMTM, 100*ftol, &nfunk,  verbose.debug, finderrors, 1, dplus, dmin) != 0) {
			      if(bruteforce) {
				fflush(stdout);
				printwarning(verbose.debug, "WARNING cal_mtm: Fit algorithm didn't converge (f=%ld sub=%ld).", f, n);
				failed = 1;
			      }else {
				fflush(stdout);
				printerror(verbose.debug, "ERROR cal_mtm: Fit algorithm didn't converge (f=%ld sub=%ld).", f, n);
				return 0;
			      }
			    }else {
			      fflush(stdout);
			      printwarning(verbose.debug, "WARNING cal_mtm: Success!");
			    }
			  }
			}
			//			gsl_set_error_handler(oldhandler); 

			if(failed == 0) {
			  // Now calculate the chisq of the best
			  // solution. I will assume the rms of I,Q,U,V
			  // are the same (they should), so take the
			  // average of the rms'es. For reasons
			  // explained in the funk_optimizeMTM, no
			  // normalisation of the chi^2 was done. We
			  // will now calculate a proper chi^2, but
			  // first we have to know the rms'es of the
			  // calibrated data. In addition, by default
			  // the template was distorted to find the
			  // solution, in which case we want
			  // to do the opposite to the data here to
			  // see the result.
			  if(internal_analytic_gain_calc_flag)
			    fittedparams[0] = 1.0;
			  internal_solution_cur_channel.data[0] = fittedparams[0];
			  internal_solution_cur_channel.data[1] = fittedparams[1];
			  internal_solution_cur_channel.data[2] = fittedparams[2];
			  if(internal_leakage_fit_flag) {
			    internal_solution_cur_channel.data[3] = fittedparams[3];
			    internal_solution_cur_channel.data[4] = fittedparams[4];
			    internal_solution_cur_channel.data[5] = fittedparams[5];
			    internal_solution_cur_channel.data[6] = fittedparams[6];
			  }
			  // Re-read data_channel_full, as it will be
			  // overwritten with calibration applied in next step, so it will need to
			  // be restored first before applying the calibration of the next trial.
			  for(p = 0; p < internal_data_channel_onpulse.NrPols; p++) {
			    if(readPulsePSRData(datafile, n, p, f, 0, datafile->NrBins, &data_channel_full.data[p*data_channel_full.NrBins], verbose) == 0)
			      return 0;
			  }
			  if(applyReceiverModel(&data_channel_full, internal_solution_cur_channel, 1, 0, 0, 0, noverbose) == 0) {
			    fflush(stdout);
			    printerror(verbose.debug, "ERROR cal_mtm: Cannot apply initional solution.");
			    exit(0);
			  }

			  if(internal_analytic_gain_calc_flag) {
			    double intensity_template, intensity_data;
			    float norm;

			    intensity_template = 0;
			    for(b = 0; b < internal_template_channel_onpulse.NrBins; b++) {
			      intensity_template += internal_template_channel_onpulse.data[b];
			    }

			    intensity_data = 0;
			    for(b = 0; b < data_channel_full.NrBins; b++) {
			      if(checkRegions(b, &onpulseregion_tmp2, 0, verbose) != 0) {
				intensity_data += data_channel_full.data[b];
			      }
			    }

			    norm = intensity_template/intensity_data;
			    for(b = 0; b < 4*data_channel_full.NrBins; b++) {
			      data_channel_full.data[b] *= norm;
			    }
			    fittedparams[0] = 1.0/sqrt(norm);
			  }

			  chisq_norm = 0;
			  // A new figure of merrit needs to be calculated, as we now correct observation rather than distort the template, i.e. if we made template smaller during fitting, we now make data bigger.
			  double value;
			  newFOM = 0;
			  //			  int junk_counter = 0;
			  for(p = 0; p < 4; p++) {
			    rms_temp = 0;
			    for(b = 0; b < data_channel_full.NrBins; b++) {
			      //			      if(p == 0 && f == 161) {
			      //				printf("  calib profile %f\n", data_channel_full.data[b+p*data_channel_full.NrBins]);
			      //			      }
			      if(checkRegions(b, &onpulseregion_tmp, 0, verbose) == 0) {
				// Note: baseline is already subtracted 
				rms_temp += (double)(data_channel_full.data[b+p*data_channel_full.NrBins])*(double)(data_channel_full.data[b+p*data_channel_full.NrBins]);
				//				printf("  Added sample %f\n", data_channel_full.data[b+p*data_channel_full.NrBins]);
			      }
			      // If in fit region, add to FOM
			      if(checkRegions(b, &onpulseregion_tmp2, 0, verbose) != 0) {
				value = template_channel_full[b+p*template->NrBins] - data_channel_full.data[b+p*data_channel_full.NrBins];
				newFOM += value*value;
				//				junk_counter++;
			      }
			    }
			    //			    printf("XXXXXX rms_temp=%f, nrbins=%ld\n", rms_temp, nractivebins_offpulse);
			    rms_temp /= (double)nractivebins_offpulse;
			    internal_data_channel_onpulse_calib.offpulse_rms[p] = sqrt(rms_temp);
			    //			    printf("XXXXXX RMS pol %ld: template %f, obs %f\n", p, internal_template_channel_onpulse.offpulse_rms[p], internal_data_channel_onpulse_calib.offpulse_rms[p]);
			    chisq_norm += internal_template_channel_onpulse.offpulse_rms[p]*internal_template_channel_onpulse.offpulse_rms[p]+internal_data_channel_onpulse_calib.offpulse_rms[p]*internal_data_channel_onpulse_calib.offpulse_rms[p];
			  }
			  //			  printerror(verbose.debug, "XXXXX %d", junk_counter);
			  // Normalise by average of rms I,Q,U,V. 
			  //			  printf("XXXXXX G=%f\n", fittedparams[0]);
			  //			  printf("XXXXXX chisq = %f/%f ", chisq, 0.25*chisq_norm);
			  //			  printf("XXXXXX: chisq=%f, chisq_norm=%f\n", chisq, chisq_norm);
			  chisq_norm = newFOM/(0.25*chisq_norm);
			  //			printf("= %f\n", chisq_norm);
			  // Not sure you always want to do this. When you combine receiver solutions you don't 
			  // want the result to be dominated by crappy data with low chisq. Also, the RMS'es 
			  // are influenced by solution, so now you pick
			  // the "best" solution, which might be a crappy one with large rms'es
			  //			chisq = chisq_norm;
			}
		      } // end of if(iszapped == 0)
		      if(failed == 0) {
			if(chisqbest < 0 || chisq < chisqbest) {
			  /*
			  if(f == 5 && chisq < 7402) {
			    output_FittedProfile_Of_Itterations = 1;
			    funk_optimizeMTM(fittedparams);
			    output_FittedProfile_Of_Itterations = 0;
			    }*/
			  //			  printf("XXXXXX found new best FOM: %lf\n", chisq);
			  fittedparamsbest[0] = fittedparams[0];
			  fittedparamsbest[1] = fittedparams[1];
			  fittedparamsbest[2] = fittedparams[2];
			  fittedparamsbest[3] = fittedparams[3];
			  fittedparamsbest[4] = fittedparams[4];
			  fittedparamsbest[5] = fittedparams[5];
			  fittedparamsbest[6] = fittedparams[6];
			  chisqbest = chisq;
			  chisqbest_norm = chisq_norm;
			  data_fom_best = newFOM;
			  // Store result of best solution
			  for(p = 0; p < internal_data_channel_onpulse.NrPols; p++) {
			    curbin = 0;
			    for(b = 0; b < datafile->NrBins; b++) {
			      if(checkRegions(b, &onpulseregion_tmp2, 0, verbose) != 0) {
				data_channel_calib_best.data[p*internal_data_channel_onpulse.NrBins+curbin] = data_channel_full.data[b+p*data_channel_full.NrBins];
				curbin++;
			      }
			    }
			  }
			  //			  if(f == 161)
			  //			    printf("freq %ld: Accepted solution with chisq=%f\n", f, chisqbest);
			}
			//			else {
			  //			  if(f == 161)
			  //			    printf("freq %ld: Rejected solution with chisq=%f\n", f, chisqbest);
			//			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }else { /* End oktoprocesschannel if statement */ // This is never reached? Similar if statement exist above
	fittedparamsbest[0] = 1;
	fittedparamsbest[1] = 0;
	fittedparamsbest[2] = 0;
	fittedparamsbest[3] = 0;
	fittedparamsbest[4] = 0;
	fittedparamsbest[5] = 0;
	fittedparamsbest[6] = 0;
	chisqbest = -1;
	chisqbest_norm = -1;
	iszapped = 1;
      }

      if(iszapped) {
	chisqbest = -1;
	chisqbest_norm = -1;
	data_fom_best = -1;
      }

      if(verbose.verbose) {
	if(f == 0) {
	  for(i = 0; i < verbose.indent; i++)      
	    printf(" ");
	  printf("G=gain, differential gain gamma=exp(2dG)-1, dph=differential phase [rad], el1/2=ellipticity [rad], or1/2=orientation [rad] FOMt=Fig. of merrit fitting template to data (used to optimise solution), FOMd=FOM when applying solution to data, chisq=reduced chi2 corresponding to FOMd, nfree=Degr. of freedom in fit, iszapped=1 if channel is zapped.\n");
	}
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	if(leakage == 0)
	  printf("  Solution nsub=%04ld freq=%04ld: G=%6.3f  dG=%6.3f   dph=%6.3f   (FOMt=%10.3e  FOMd=%10.3e  chisq=%10.3e  nfree=%d  iszapped=%d)\n", n, f, fittedparamsbest[0], fittedparamsbest[1], fittedparamsbest[2], chisqbest, data_fom_best, chisqbest_norm/(float)internal_nfree, internal_nfree, iszapped);
	else
	  printf("  Solution nsub=%04ld freq=%04ld: G=%6.3f  dG=%6.3f   dph=%6.3f   el1=%6.3f   or1=%6.3f   el2=%6.3f   or2=%6.3f   (FOMt=%10.3e  FOMd=%10.3e  chisq=%10.3e  nfree=%d  iszapped=%d)\n", n, f, fittedparamsbest[0], fittedparamsbest[1], fittedparamsbest[2], fittedparamsbest[3], fittedparamsbest[4], fittedparamsbest[5], fittedparamsbest[6], chisqbest, data_fom_best, chisqbest_norm/(float)internal_nfree, internal_nfree, iszapped);
      }

      if(writePulsePSRData(&solutionall, n, 0, f, 0, 1, &fittedparamsbest[0], verbose) != 1) {
	fflush(stdout);
 	printerror(verbose.debug, "ERROR cal_mtm: Error writing solution.");
	return 0;
      }
      if(writePulsePSRData(&solutionall, n, 1, f, 0, 1, &fittedparamsbest[1], verbose) != 1) {
	fflush(stdout);
 	printerror(verbose.debug, "ERROR cal_mtm: Error writing solution.");
	return 0;
      }
      if(writePulsePSRData(&solutionall, n, 2, f, 0, 1, &fittedparamsbest[2], verbose) != 1) {
	fflush(stdout);
 	printerror(verbose.debug, "ERROR cal_mtm: Error writing solution.");
	return 0;
      }
      if(leakage) {
	if(writePulsePSRData(&solutionall, n, 3, f, 0, 1, &fittedparamsbest[3], verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR cal_mtm: Error writing solution.");
	  return 0;
	}
	if(writePulsePSRData(&solutionall, n, 4, f, 0, 1, &fittedparamsbest[4], verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR cal_mtm: Error writing solution.");
	  return 0;
	}
	if(writePulsePSRData(&solutionall, n, 5, f, 0, 1, &fittedparamsbest[5], verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR cal_mtm: Error writing solution.");
	  return 0;
	}
	if(writePulsePSRData(&solutionall, n, 6, f, 0, 1, &fittedparamsbest[6], verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR cal_mtm: Error writing solution.");
	  return 0;
	}
	if(writePulsePSRData(&solutionall, n, 7, f, 0, 1, &chisqbest_norm, verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR cal_mtm: Error writing solution.");
	  return 0;
	}
	float valuef;
	valuef = internal_nfree;
	if(writePulsePSRData(&solutionall, n, 8, f, 0, 1, &valuef, verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR cal_mtm: Error writing solution.");
	  return 0;
	}
      }

      if((verbose.debug && iszapped == 0)) {
	printf("  After %d steps I found a solution with a chi^2 of %lf\n", nfunk, chisq);
	//	cal_internal_drawTwoProfiles(internal_data_channel_onpulse_calib.data, internal_template_channel_onpulse.data, internal_data_channel_onpulse.NrBins, "Data (white) and template (red) (prefit; initial guess gain)");

	cal_internal_drawBeforeAfter(internal_data_channel_onpulse, -1, 0, "", beforeafter_device, 0);
	cal_internal_drawBeforeAfter(internal_template_channel_onpulse, 0, 0, "IQUV=WRGB before: (dots are data)", beforeafter_device, 0);
	cal_internal_drawBeforeAfter(internal_data_channel_onpulse, 1, 0, "", beforeafter_device, 0);
	cal_internal_drawBeforeAfter(internal_template_channel_onpulse, 2, 0, "IQUV=WRGB after: (dots are data)", beforeafter_device, 0);
	cal_internal_drawBeforeAfter(data_channel_calib_best, 3, 0, "", beforeafter_device, 0);
	//	for(b = 0; b < data_channel_calib_best.NrBins*data_channel_calib_best.NrPols; b++) {
	  //	  printf("XXXXX %ld %f %f\n", b, data_channel_calib_best.data[b], internal_template_channel_onpulse.data[b]);
	//	}	
	printf("Press a key in terminal to continue\n");
	pgetch();
      }

      // Add solution to the overall frequency-scrunched addition
      if(verbose.verbose && iszapped == 0) {
	long i1, i2, i3;

	//De-Faraday rotate data to get a nice plot without loosing the linear polarization. Undo the applied Faraday rotation to template and apply the same de-rotation to the data (relative to infinite freq).
	double dphi, L, phi;
	dphi = calcRMAngle(get_nonweighted_channel_freq(*datafile, f, verbose), 0.0, 1, datafile->rm);
	for(b = 0; b < internal_template_channel_onpulse.NrBins; b++) {
	  L = sqrt(internal_template_channel_onpulse.data[internal_template_channel_onpulse.NrBins+b]*internal_template_channel_onpulse.data[internal_template_channel_onpulse.NrBins+b]+internal_template_channel_onpulse.data[2*internal_template_channel_onpulse.NrBins+b]*internal_template_channel_onpulse.data[2*internal_template_channel_onpulse.NrBins+b]);
	  phi = atan2(internal_template_channel_onpulse.data[2*internal_template_channel_onpulse.NrBins+b], internal_template_channel_onpulse.data[internal_template_channel_onpulse.NrBins+b]);
	  phi -= 2.0*dphi;
	  internal_template_channel_onpulse.data[internal_template_channel_onpulse.NrBins+b] = L*cos(phi);
	  internal_template_channel_onpulse.data[2*internal_template_channel_onpulse.NrBins+b] = L*sin(phi);

	  L = sqrt(data_channel_calib_best.data[data_channel_calib_best.NrBins+b]*data_channel_calib_best.data[data_channel_calib_best.NrBins+b]+data_channel_calib_best.data[2*data_channel_calib_best.NrBins+b]*data_channel_calib_best.data[2*data_channel_calib_best.NrBins+b]);
	  phi = atan2(data_channel_calib_best.data[2*data_channel_calib_best.NrBins+b], data_channel_calib_best.data[data_channel_calib_best.NrBins+b]);
	  phi -= 2.0*dphi;
	  data_channel_calib_best.data[data_channel_calib_best.NrBins+b] = L*cos(phi);
	  data_channel_calib_best.data[2*data_channel_calib_best.NrBins+b] = L*sin(phi);
	}

	for(i1 = 0; i1 < 4; i1++) {
	  for(i2 = 0; i2 < fscrunched_calib_obs.NrBins; i2++) {
	    //	    if(f == 150) {
	      i3 = i1*fscrunched_calib_obs.NrBins+i2;
	      fscrunched_calib_obs.data[i3] += data_channel_calib_best.data[i3];
	      fscrunched_template.data[i3] += internal_template_channel_onpulse.data[i3];
	      //	    }
	  }
	}
      }

    } // End loop over frequency channels
  } // End loop over subints
  if(sequencenr > 0 && totsequencenr > 0) {
    fclose(dumpfile);
  }

  /* Open output file and write data */
  char realoutputfilename[1000];
  if(sequencenr > 0 && totsequencenr > 0) {
    sprintf(realoutputfilename, "%s.%d", outputfilename, sequencenr);
  }else {
    sprintf(realoutputfilename, "%s", outputfilename);
  }
  // Set the reference frequency to be the centre frequency of the observation.
  // Required for the output fits file to be handled by psrchive, as otherwise it is set to infinite frequency (a very large number).
  solutionall.freq_ref = get_centre_frequency(*datafile, verbose);
  if(!openPSRData(&solutionall, realoutputfilename, FITS_format, 1, 0, 0, verbose))
    return 0;
  int cmdOnly = 0;
  if(!writeHeaderPSRData(&solutionall, argc, argv, cmdOnly, verbose))
    return 0;
  //  writeHistory(&solutionall, argc, argv, verbose);
  if(writePSRData(&solutionall, solutionall.data, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm: Cannot write data");
    return 0;
  }
  closePSRData(&solutionall, 0, verbose);

  if(verbose.verbose || beforeafter_device != NULL || profileresid_device) {
    /*
    if(preprocess_norm(fscrunched_calib_obs, 1.0, NULL, 0, verbose) != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cal_mtm: Normalisation failed");
      return 0;
    }
    if(preprocess_norm(fscrunched_template, 1.0, NULL, 0, verbose) != 1) {
	  fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Normalisation failed");
      return 0;
    }
    */

    char *title;
    if(datafile->filename != NULL) {
      title = malloc(100+strlen(datafile->filename));
    }else {
      title = malloc(100);
    }
    if(title == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cal_mtm: Memory allocation error.");
      return 0;
    }
    if(datafile->filename != NULL) {
      sprintf(title, "%s IQUV=WRGB freq-scrunched result: (dots are data)", datafile->filename);
    }else {
      sprintf(title, "IQUV=WRGB freq-scrunched result: (dots are data)");
    }

    cal_internal_drawBeforeAfter(fscrunched_calib_obs, -1, 1, "", beforeafter_device, profileresid_device);
    cal_internal_drawBeforeAfter(fscrunched_template, 0, 1, title, beforeafter_device, profileresid_device);
    cal_internal_drawBeforeAfter(fscrunched_calib_obs, 1, 1, "", beforeafter_device, profileresid_device);
    if(profileresid_device == 0) {
      ppgend();
    }
    free(title);

    double FOM = 0;
    double delta;
    i = 0;
    for(p = 0; p < 4; p++) { 
      for(b = 0; b < fscrunched_calib_obs.NrBins; b++) {
	delta = fscrunched_calib_obs.data[i]-fscrunched_template.data[i];
	i++;
	FOM += delta*delta;
      }
    }

    //    printf("Overal FOM: approx. %e (peak of template/calibrated data were normalised first)\n", FOM);
    printf("  Overall FOM: %e (based on frequency summed profile, lower is better)\n", FOM);
  }

  //  ppgend();   // Took this out to have multiple profile plots work. I don't think it is necessary, but maybe something broke.

  freePulselongitudeRegion(&onpulseregion_tmp);
  freePulselongitudeRegion(&onpulseregion_tmp2);
  free(pgplot_options);
  free(template_channel_full);
  closePSRData(&data_channel_full, 0, verbose);
  closePSRData(&fscrunched_template, 0, verbose);
  closePSRData(&fscrunched_calib_obs, 0, verbose);
  closePSRData(&data_channel_calib_best, 0, verbose);
  closePSRData(&internal_template_channel_onpulse_calib, 0, verbose);
  closePSRData(&internal_template_channel_onpulse, 0, verbose);
  closePSRData(&internal_data_channel_onpulse_calib, 0, verbose);
  closePSRData(&internal_data_channel_onpulse, 0, verbose);
  closePSRData(&internal_solution_cur_channel, 0, verbose);
  return 1;
}  // End of cal_mtm


double refmodel_phase(double freq, long refmodel, verbose_definition verbose)
{
  splineFit_def splineFit;
  double ans;
  int ret;
  if(refmodel == 20120126) {
    splineFit.order = 4;
    splineFit.nbreak = 23;
    splineFit.x_start = 1.323410e+03;
    splineFit.x_end = 1.726890e+03;
    splineFit.coefficients = (double *)malloc(25*sizeof(double));
    splineFit.breakpoints = (double *)malloc(23*sizeof(double));
    if(splineFit.coefficients == NULL || splineFit.breakpoints == NULL) {
      fflush(stdout);
      printerror(0, "ERROR refmodel_phase: Cannot allocate memory for model %ld", refmodel);
      exit(0);
    }
    splineFit.coefficients[0] = 4.859596e+01; splineFit.coefficients[1] = 1.288068e+00; splineFit.coefficients[2] = 1.129416e+01; splineFit.coefficients[3] = -3.082969e+01; splineFit.coefficients[4] = -1.412796e+02; splineFit.coefficients[5] = 5.107659e+02; splineFit.coefficients[6] = 4.123721e+07; splineFit.coefficients[7] = -2.006381e+05; splineFit.coefficients[8] = 5.388855e+03; splineFit.coefficients[9] = -9.544992e+01; splineFit.coefficients[10] = 1.290026e+02; splineFit.coefficients[11] = -1.560441e+01; splineFit.coefficients[12] = 5.128615e+00; splineFit.coefficients[13] = -2.273511e+02; splineFit.coefficients[14] = 7.740572e+00; splineFit.coefficients[15] = -7.423763e+01; splineFit.coefficients[16] = -4.577323e+01; splineFit.coefficients[17] = -1.141388e+02; splineFit.coefficients[18] = -1.426467e+02; splineFit.coefficients[19] = -1.816209e+02; splineFit.coefficients[20] = -1.777312e+02; splineFit.coefficients[21] = -1.879414e+02; splineFit.coefficients[22] = -1.843166e+02; splineFit.coefficients[23] = -1.877673e+02; splineFit.coefficients[24] = -1.542505e+02;
    splineFit.breakpoints[0] = 1.323410e+03; splineFit.breakpoints[1] = 1.425183e+03; splineFit.breakpoints[2] = 1.425705e+03; splineFit.breakpoints[3] = 1.443955e+03; splineFit.breakpoints[4] = 1.453464e+03; splineFit.breakpoints[5] = 1.482323e+03; splineFit.breakpoints[6] = 1.482828e+03; splineFit.breakpoints[7] = 1.526768e+03; splineFit.breakpoints[8] = 1.528129e+03; splineFit.breakpoints[9] = 1.528387e+03; splineFit.breakpoints[10] = 1.529755e+03; splineFit.breakpoints[11] = 1.530184e+03; splineFit.breakpoints[12] = 1.530962e+03; splineFit.breakpoints[13] = 1.531683e+03; splineFit.breakpoints[14] = 1.531871e+03; splineFit.breakpoints[15] = 1.535590e+03; splineFit.breakpoints[16] = 1.539895e+03; splineFit.breakpoints[17] = 1.545857e+03; splineFit.breakpoints[18] = 1.560766e+03; splineFit.breakpoints[19] = 1.621085e+03; splineFit.breakpoints[20] = 1.639320e+03; splineFit.breakpoints[21] = 1.655775e+03; splineFit.breakpoints[22] = 1.726890e+03;
  }else if(refmodel == 20101009) {
    splineFit.order = 4;
    splineFit.nbreak = 23;
    splineFit.x_start = 1.319250e+03;
    splineFit.x_end = 1.726250e+03;
    splineFit.coefficients = malloc(25*sizeof(double));
    splineFit.breakpoints = malloc(23*sizeof(double));
    if(splineFit.coefficients == NULL || splineFit.breakpoints == NULL) {
      fflush(stdout);
      printerror(0, "ERROR refmodel_phase: Cannot allocate memory for model %ld", refmodel);
      exit(0);
    }
    splineFit.coefficients[0] = 4.015991e+01; splineFit.coefficients[1] = 3.227456e+01; splineFit.coefficients[2] = 3.498294e+01; splineFit.coefficients[3] = 4.921152e+00; splineFit.coefficients[4] = 1.066882e+01; splineFit.coefficients[5] = 9.509486e+00; splineFit.coefficients[6] = -1.506720e+00; splineFit.coefficients[7] = -7.624645e+00; splineFit.coefficients[8] = -1.003923e+01; splineFit.coefficients[9] = -3.265026e+01; splineFit.coefficients[10] = -2.665463e+01; splineFit.coefficients[11] = -3.234767e+03; splineFit.coefficients[12] = 3.528427e+03; splineFit.coefficients[13] = -3.374283e+03; splineFit.coefficients[14] = 3.614096e+02; splineFit.coefficients[15] = 5.883754e+01; splineFit.coefficients[16] = 7.833478e+00; splineFit.coefficients[17] = 9.517329e+00; splineFit.coefficients[18] = -1.127469e+01; splineFit.coefficients[19] = 3.241603e+00; splineFit.coefficients[20] = -5.780397e+00; splineFit.coefficients[21] = 9.394108e-01; splineFit.coefficients[22] = -4.590822e+00; splineFit.coefficients[23] = 1.394599e+01; splineFit.coefficients[24] = 1.992430e+01;
    splineFit.breakpoints[0] = 1.319250e+03; splineFit.breakpoints[1] = 1.336746e+03; splineFit.breakpoints[2] = 1.369457e+03; splineFit.breakpoints[3] = 1.371033e+03; splineFit.breakpoints[4] = 1.374033e+03; splineFit.breakpoints[5] = 1.393250e+03; splineFit.breakpoints[6] = 1.411750e+03; splineFit.breakpoints[7] = 1.421752e+03; splineFit.breakpoints[8] = 1.422900e+03; splineFit.breakpoints[9] = 1.463257e+03; splineFit.breakpoints[10] = 1.475648e+03; splineFit.breakpoints[11] = 1.496193e+03; splineFit.breakpoints[12] = 1.504515e+03; splineFit.breakpoints[13] = 1.522750e+03; splineFit.breakpoints[14] = 1.541250e+03; splineFit.breakpoints[15] = 1.547634e+03; splineFit.breakpoints[16] = 1.580088e+03; splineFit.breakpoints[17] = 1.604216e+03; splineFit.breakpoints[18] = 1.631902e+03; splineFit.breakpoints[19] = 1.634590e+03; splineFit.breakpoints[20] = 1.652250e+03; splineFit.breakpoints[21] = 1.692053e+03; splineFit.breakpoints[22] = 1.726250e+03;
  }else if(refmodel == 20090816) {
    splineFit.order = 4;
    splineFit.nbreak = 23;
    splineFit.x_start = 1.300000e+03;
    splineFit.x_end = 1.750000e+03;
    splineFit.coefficients = malloc(25*sizeof(double));
    splineFit.breakpoints = malloc(23*sizeof(double));
    if(splineFit.coefficients == NULL || splineFit.breakpoints == NULL) {
      fflush(stdout);
      printerror(0, "ERROR refmodel_phase: Cannot allocate memory for model %ld", refmodel);
      exit(0);
    }
    splineFit.coefficients[0] = 7.910186e+01; splineFit.coefficients[1] = 8.895594e+01; splineFit.coefficients[2] = 3.103228e+01; splineFit.coefficients[3] = 2.537357e+01; splineFit.coefficients[4] = 1.233020e-01; splineFit.coefficients[5] = -5.061013e+01; splineFit.coefficients[6] = -1.199611e+02; splineFit.coefficients[7] = -1.287806e+03; splineFit.coefficients[8] = 6.346463e+03; splineFit.coefficients[9] = -4.321480e+02; splineFit.coefficients[10] = 6.470116e+00; splineFit.coefficients[11] = 1.573901e+01; splineFit.coefficients[12] = 2.716809e+01; splineFit.coefficients[13] = 1.532589e+01; splineFit.coefficients[14] = 1.725437e+01; splineFit.coefficients[15] = 9.871985e+00; splineFit.coefficients[16] = 1.015997e+01; splineFit.coefficients[17] = 1.654141e+01; splineFit.coefficients[18] = -1.301896e+01; splineFit.coefficients[19] = 5.259659e+00; splineFit.coefficients[20] = -5.764387e+00; splineFit.coefficients[21] = -9.766593e+00; splineFit.coefficients[22] = -1.480330e+00; splineFit.coefficients[23] = -3.910798e+01; splineFit.coefficients[24] = -1.754534e+01;
    splineFit.breakpoints[0] = 1.300000e+03; splineFit.breakpoints[1] = 1.357772e+03; splineFit.breakpoints[2] = 1.361364e+03; splineFit.breakpoints[3] = 1.404021e+03; splineFit.breakpoints[4] = 1.418892e+03; splineFit.breakpoints[5] = 1.461673e+03; splineFit.breakpoints[6] = 1.484091e+03; splineFit.breakpoints[7] = 1.486745e+03; splineFit.breakpoints[8] = 1.525000e+03; splineFit.breakpoints[9] = 1.541677e+03; splineFit.breakpoints[10] = 1.544423e+03; splineFit.breakpoints[11] = 1.565909e+03; splineFit.breakpoints[12] = 1.586364e+03; splineFit.breakpoints[13] = 1.606818e+03; splineFit.breakpoints[14] = 1.611747e+03; splineFit.breakpoints[15] = 1.617246e+03; splineFit.breakpoints[16] = 1.622603e+03; splineFit.breakpoints[17] = 1.627273e+03; splineFit.breakpoints[18] = 1.662363e+03; splineFit.breakpoints[19] = 1.668182e+03; splineFit.breakpoints[20] = 1.692180e+03; splineFit.breakpoints[21] = 1.692915e+03; splineFit.breakpoints[22] = 1.750000e+03;
  }else if(refmodel == 20120819) {
    splineFit.order = 4;
    splineFit.nbreak = 28;
    splineFit.x_start = 1.358050e+03;
    splineFit.x_end = 1.725450e+03;
    splineFit.coefficients = malloc(30*sizeof(double));
    splineFit.breakpoints = malloc(28*sizeof(double));
    if(splineFit.coefficients == NULL || splineFit.breakpoints == NULL) {
      fflush(stdout);
      printerror(0, "ERROR refmodel_phase: Cannot allocate memory for model %ld", refmodel);
      exit(0);
    }
    splineFit.coefficients[0] = 8.119092e+00; splineFit.coefficients[1] = 6.688840e+04; splineFit.coefficients[2] = -2.924567e+02; splineFit.coefficients[3] = 6.610589e+01; splineFit.coefficients[4] = -8.476670e+00; splineFit.coefficients[5] = 1.852434e+01; splineFit.coefficients[6] = -1.521222e+00; splineFit.coefficients[7] = 9.046425e+01; splineFit.coefficients[8] = -1.786008e+01; splineFit.coefficients[9] = 8.597220e+00; splineFit.coefficients[10] = -7.962144e+00; splineFit.coefficients[11] = 2.807497e+00; splineFit.coefficients[12] = -1.256141e+01; splineFit.coefficients[13] = 3.126480e+01; splineFit.coefficients[14] = -4.585746e+02; splineFit.coefficients[15] = 3.013664e+05; splineFit.coefficients[16] = -5.750043e+04; splineFit.coefficients[17] = 3.530249e+02; splineFit.coefficients[18] = -1.243538e+01; splineFit.coefficients[19] = -7.520005e+00; splineFit.coefficients[20] = -5.340787e+00; splineFit.coefficients[21] = -1.296171e+00; splineFit.coefficients[22] = -4.866165e+00; splineFit.coefficients[23] = -3.015304e-01; splineFit.coefficients[24] = 1.847827e+00; splineFit.coefficients[25] = 2.649273e+00; splineFit.coefficients[26] = -1.173209e+00; splineFit.coefficients[27] = 8.511918e+00; splineFit.coefficients[28] = -5.142319e+00; splineFit.coefficients[29] = 1.743925e+00;
    splineFit.breakpoints[0] = 1.358050e+03; splineFit.breakpoints[1] = 1.371393e+03; splineFit.breakpoints[2] = 1.375649e+03; splineFit.breakpoints[3] = 1.382521e+03; splineFit.breakpoints[4] = 1.384551e+03; splineFit.breakpoints[5] = 1.386398e+03; splineFit.breakpoints[6] = 1.387311e+03; splineFit.breakpoints[7] = 1.406328e+03; splineFit.breakpoints[8] = 1.406460e+03; splineFit.breakpoints[9] = 1.408924e+03; splineFit.breakpoints[10] = 1.411154e+03; splineFit.breakpoints[11] = 1.413752e+03; splineFit.breakpoints[12] = 1.433360e+03; splineFit.breakpoints[13] = 1.469725e+03; splineFit.breakpoints[14] = 1.477386e+03; splineFit.breakpoints[15] = 1.481242e+03; splineFit.breakpoints[16] = 1.517608e+03; splineFit.breakpoints[17] = 1.522620e+03; splineFit.breakpoints[18] = 1.538557e+03; splineFit.breakpoints[19] = 1.552142e+03; splineFit.breakpoints[20] = 1.556472e+03; splineFit.breakpoints[21] = 1.562095e+03; splineFit.breakpoints[22] = 1.575084e+03; splineFit.breakpoints[23] = 1.639390e+03; splineFit.breakpoints[24] = 1.648146e+03; splineFit.breakpoints[25] = 1.650904e+03; splineFit.breakpoints[26] = 1.666191e+03; splineFit.breakpoints[27] = 1.725450e+03;
  }else {
    printerror(0, "ERROR refmodel_phase: Unknown reference model: %ld", refmodel);
    exit(0);
  }
  ret = cubicBspline_eval_double(splineFit, freq, &ans, verbose);
  cubicBspline_free(&splineFit);
  if(ret != 1) {
    fflush(stdout);
    printerror(0, "ERROR refmodel_phase: interpolation failed");
    exit(0);
  }
  return ans;
}

double refmodel_diffgain(double freq, long refmodel, verbose_definition verbose)
{
  splineFit_def splineFit;
  double ans;
  int ret;
  if(refmodel == 20120126) {
    splineFit.order = 4;
    splineFit.nbreak = 18;
    splineFit.x_start = 1.300000e+03;
    splineFit.x_end = 1.750000e+03;
    splineFit.coefficients = malloc(20*sizeof(double));
    splineFit.breakpoints = malloc(18*sizeof(double));
    if(splineFit.coefficients == NULL || splineFit.breakpoints == NULL) {
      fflush(stdout);
      printerror(0, "ERROR refmodel_diffgain: Cannot allocate memory for model %ld", refmodel);
      exit(0);
    }
    splineFit.coefficients[0] = 4.493388e+00; splineFit.coefficients[1] = -2.017354e+01; splineFit.coefficients[2] = 2.019120e+01; splineFit.coefficients[3] = 1.018483e+01; splineFit.coefficients[4] = 1.471308e+01; splineFit.coefficients[5] = -1.702146e+00; splineFit.coefficients[6] = 1.257317e+01; splineFit.coefficients[7] = -3.868719e+02; splineFit.coefficients[8] = 1.084216e+03; splineFit.coefficients[9] = -9.450899e+01; splineFit.coefficients[10] = -8.091449e+01; splineFit.coefficients[11] = -8.319187e+01; splineFit.coefficients[12] = -1.148812e+01; splineFit.coefficients[13] = 4.095772e+00; splineFit.coefficients[14] = -2.020078e+01; splineFit.coefficients[15] = 1.416660e+00; splineFit.coefficients[16] = -1.755330e+00; splineFit.coefficients[17] = 3.258445e+01; splineFit.coefficients[18] = 4.218565e+00; splineFit.coefficients[19] = 2.570999e+01;
    splineFit.breakpoints[0] = 1.300000e+03; splineFit.breakpoints[1] = 1.341811e+03; splineFit.breakpoints[2] = 1.414866e+03; splineFit.breakpoints[3] = 1.416833e+03; splineFit.breakpoints[4] = 1.426632e+03; splineFit.breakpoints[5] = 1.431300e+03; splineFit.breakpoints[6] = 1.432353e+03; splineFit.breakpoints[7] = 1.506894e+03; splineFit.breakpoints[8] = 1.523077e+03; splineFit.breakpoints[9] = 1.532344e+03; splineFit.breakpoints[10] = 1.535770e+03; splineFit.breakpoints[11] = 1.538235e+03; splineFit.breakpoints[12] = 1.542393e+03; splineFit.breakpoints[13] = 1.607104e+03; splineFit.breakpoints[14] = 1.640970e+03; splineFit.breakpoints[15] = 1.644267e+03; splineFit.breakpoints[16] = 1.682920e+03; splineFit.breakpoints[17] = 1.750000e+03;
  }else if(refmodel == 20111003) {
    splineFit.order = 4;
    splineFit.nbreak = 23;
    splineFit.x_start = 1.300000e+03;
    splineFit.x_end = 1.800000e+03;
    splineFit.coefficients = malloc(25*sizeof(double));
    splineFit.breakpoints = malloc(23*sizeof(double));
    if(splineFit.coefficients == NULL || splineFit.breakpoints == NULL) {
      fflush(stdout);
      printerror(0, "ERROR refmodel_diffgain: Cannot allocate memory for model %ld", refmodel);
      exit(0);
    }
    splineFit.coefficients[0] = -3.988807e+01; splineFit.coefficients[1] = -4.848314e+01; splineFit.coefficients[2] = -7.029112e+00; splineFit.coefficients[3] = 1.272247e+01; splineFit.coefficients[4] = 6.473987e+01; splineFit.coefficients[5] = -3.224320e+01; splineFit.coefficients[6] = 1.908303e+02; splineFit.coefficients[7] = -2.488187e+03; splineFit.coefficients[8] = 8.644564e+03; splineFit.coefficients[9] = -3.848955e+02; splineFit.coefficients[10] = -1.452885e+01; splineFit.coefficients[11] = -3.019950e+00; splineFit.coefficients[12] = -1.350683e+01; splineFit.coefficients[13] = 7.443214e+00; splineFit.coefficients[14] = -1.300427e+01; splineFit.coefficients[15] = 5.856877e+00; splineFit.coefficients[16] = -1.540511e+01; splineFit.coefficients[17] = 1.058880e+01; splineFit.coefficients[18] = -1.133736e+01; splineFit.coefficients[19] = -1.170213e+01; splineFit.coefficients[20] = -2.515949e+01; splineFit.coefficients[21] = 7.711589e-01; splineFit.coefficients[22] = -3.619288e+01; splineFit.coefficients[23] = 1.894317e+01; splineFit.coefficients[24] = -2.035879e+01;
    splineFit.breakpoints[0] = 1.300000e+03; splineFit.breakpoints[1] = 1.359215e+03; splineFit.breakpoints[2] = 1.370561e+03; splineFit.breakpoints[3] = 1.399268e+03; splineFit.breakpoints[4] = 1.413636e+03; splineFit.breakpoints[5] = 1.451299e+03; splineFit.breakpoints[6] = 1.460187e+03; splineFit.breakpoints[7] = 1.498853e+03; splineFit.breakpoints[8] = 1.526238e+03; splineFit.breakpoints[9] = 1.537634e+03; splineFit.breakpoints[10] = 1.550000e+03; splineFit.breakpoints[11] = 1.586010e+03; splineFit.breakpoints[12] = 1.586588e+03; splineFit.breakpoints[13] = 1.590407e+03; splineFit.breakpoints[14] = 1.609120e+03; splineFit.breakpoints[15] = 1.626079e+03; splineFit.breakpoints[16] = 1.638907e+03; splineFit.breakpoints[17] = 1.650336e+03; splineFit.breakpoints[18] = 1.656506e+03; splineFit.breakpoints[19] = 1.664893e+03; splineFit.breakpoints[20] = 1.666493e+03; splineFit.breakpoints[21] = 1.680302e+03; splineFit.breakpoints[22] = 1.800000e+03;
  }else if(refmodel == 20111029) {
    splineFit.order = 4;
    splineFit.nbreak = 28;
    splineFit.x_start = 1.300000e+03;
    splineFit.x_end = 1.800000e+03;
    splineFit.coefficients = malloc(30*sizeof(double));
    splineFit.breakpoints = malloc(28*sizeof(double));
    if(splineFit.coefficients == NULL || splineFit.breakpoints == NULL) {
      fflush(stdout);
      printerror(0, "ERROR refmodel_diffgain: Cannot allocate memory for model %ld", refmodel);
      exit(0);
    }
    splineFit.coefficients[0] = -2.940648e+01; splineFit.coefficients[1] = -4.046898e-01; splineFit.coefficients[2] = -1.454796e+01; splineFit.coefficients[3] = 4.231329e+01; splineFit.coefficients[4] = 2.130416e+01; splineFit.coefficients[5] = 4.894985e+01; splineFit.coefficients[6] = 2.073708e+00; splineFit.coefficients[7] = 1.161830e+01; splineFit.coefficients[8] = -1.239231e+00; splineFit.coefficients[9] = 4.440942e+01; splineFit.coefficients[10] = -2.349959e+03; splineFit.coefficients[11] = 1.247811e+05; splineFit.coefficients[12] = -3.338316e+04; splineFit.coefficients[13] = 1.205427e+03; splineFit.coefficients[14] = -2.340035e+02; splineFit.coefficients[15] = -4.343715e+01; splineFit.coefficients[16] = -2.328265e+01; splineFit.coefficients[17] = -1.311590e+01; splineFit.coefficients[18] = -5.779602e+00; splineFit.coefficients[19] = -1.058578e+01; splineFit.coefficients[20] = 1.083812e+01; splineFit.coefficients[21] = -1.610325e+01; splineFit.coefficients[22] = 4.499790e+00; splineFit.coefficients[23] = -6.638076e+00; splineFit.coefficients[24] = 1.092505e+01; splineFit.coefficients[25] = -1.458500e+01; splineFit.coefficients[26] = 2.640306e+01; splineFit.coefficients[27] = -9.255648e+01; splineFit.coefficients[28] = 2.508970e+02; splineFit.coefficients[29] = -5.000000e+01;
    splineFit.breakpoints[0] = 1.300000e+03; splineFit.breakpoints[1] = 1.380308e+03; splineFit.breakpoints[2] = 1.386918e+03; splineFit.breakpoints[3] = 1.393067e+03; splineFit.breakpoints[4] = 1.402015e+03; splineFit.breakpoints[5] = 1.411111e+03; splineFit.breakpoints[6] = 1.423321e+03; splineFit.breakpoints[7] = 1.426402e+03; splineFit.breakpoints[8] = 1.432012e+03; splineFit.breakpoints[9] = 1.462299e+03; splineFit.breakpoints[10] = 1.483774e+03; splineFit.breakpoints[11] = 1.504005e+03; splineFit.breakpoints[12] = 1.511048e+03; splineFit.breakpoints[13] = 1.535393e+03; splineFit.breakpoints[14] = 1.537071e+03; splineFit.breakpoints[15] = 1.542278e+03; splineFit.breakpoints[16] = 1.555212e+03; splineFit.breakpoints[17] = 1.607899e+03; splineFit.breakpoints[18] = 1.619902e+03; splineFit.breakpoints[19] = 1.625640e+03; splineFit.breakpoints[20] = 1.628050e+03; splineFit.breakpoints[21] = 1.633333e+03; splineFit.breakpoints[22] = 1.645512e+03; splineFit.breakpoints[23] = 1.659378e+03; splineFit.breakpoints[24] = 1.666805e+03; splineFit.breakpoints[25] = 1.690266e+03; splineFit.breakpoints[26] = 1.741602e+03; splineFit.breakpoints[27] = 1.800000e+03;
  }else if(refmodel == 20130529) {
    splineFit.order = 4;
    splineFit.nbreak = 18;
    splineFit.x_start = 1.300000e+03;
    splineFit.x_end = 1.800000e+03;
    splineFit.coefficients = malloc(20*sizeof(double));
    splineFit.breakpoints = malloc(18*sizeof(double));
    if(splineFit.coefficients == NULL || splineFit.breakpoints == NULL) {
      fflush(stdout);
      printerror(0, "ERROR refmodel_diffgain: Cannot allocate memory for model %ld", refmodel);
      exit(0);
    }
    splineFit.coefficients[0] = 4.995084e+01; splineFit.coefficients[1] = 1.798997e+02; splineFit.coefficients[2] = -2.438269e+01; splineFit.coefficients[3] = 4.463499e+00; splineFit.coefficients[4] = -9.567145e+00; splineFit.coefficients[5] = 1.065858e+01; splineFit.coefficients[6] = -5.230518e+00; splineFit.coefficients[7] = -1.714824e+01; splineFit.coefficients[8] = 9.393130e+00; splineFit.coefficients[9] = -5.308989e+00; splineFit.coefficients[10] = 6.701936e+00; splineFit.coefficients[11] = -5.225284e+00; splineFit.coefficients[12] = 3.662463e+00; splineFit.coefficients[13] = -1.385745e+01; splineFit.coefficients[14] = 2.766813e+00; splineFit.coefficients[15] = -4.653257e+00; splineFit.coefficients[16] = 4.852618e+00; splineFit.coefficients[17] = -4.289351e+01; splineFit.coefficients[18] = 2.217414e+02; splineFit.coefficients[19] = -6.671080e-05;
    splineFit.breakpoints[0] = 1.300000e+03; splineFit.breakpoints[1] = 1.446600e+03; splineFit.breakpoints[2] = 1.450597e+03; splineFit.breakpoints[3] = 1.470759e+03; splineFit.breakpoints[4] = 1.493895e+03; splineFit.breakpoints[5] = 1.506086e+03; splineFit.breakpoints[6] = 1.545649e+03; splineFit.breakpoints[7] = 1.591238e+03; splineFit.breakpoints[8] = 1.594118e+03; splineFit.breakpoints[9] = 1.605051e+03; splineFit.breakpoints[10] = 1.646284e+03; splineFit.breakpoints[11] = 1.654109e+03; splineFit.breakpoints[12] = 1.662318e+03; splineFit.breakpoints[13] = 1.676718e+03; splineFit.breakpoints[14] = 1.685770e+03; splineFit.breakpoints[15] = 1.689433e+03; splineFit.breakpoints[16] = 1.699456e+03; splineFit.breakpoints[17] = 1.800000e+03;
  }else if(refmodel == 20130307) {
    splineFit.order = 4;
    splineFit.nbreak = 13;
    splineFit.x_start = 1.300000e+03;
    splineFit.x_end = 1.800000e+03;
    splineFit.coefficients = malloc(15*sizeof(double));
    splineFit.breakpoints = malloc(13*sizeof(double));
    if(splineFit.coefficients == NULL || splineFit.breakpoints == NULL) {
      fflush(stdout);
      printerror(0, "ERROR refmodel_diffgain: Cannot allocate memory for model %ld", refmodel);
      exit(0);
    }
    splineFit.coefficients[0] = -2.748254e-03; splineFit.coefficients[1] = 4.402792e+01; splineFit.coefficients[2] = -4.916294e+01; splineFit.coefficients[3] = -6.625215e+00; splineFit.coefficients[4] = -1.354236e+01; splineFit.coefficients[5] = 8.410877e+00; splineFit.coefficients[6] = 2.723172e+00; splineFit.coefficients[7] = -1.474274e+01; splineFit.coefficients[8] = 1.326340e+01; splineFit.coefficients[9] = -9.101807e+00; splineFit.coefficients[10] = 6.083985e+00; splineFit.coefficients[11] = -1.334314e+00; splineFit.coefficients[12] = 3.843177e+00; splineFit.coefficients[13] = -1.441683e+01; splineFit.coefficients[14] = -3.634655e-02;
    splineFit.breakpoints[0] = 1.300000e+03; splineFit.breakpoints[1] = 1.402068e+03; splineFit.breakpoints[2] = 1.411184e+03; splineFit.breakpoints[3] = 1.413229e+03; splineFit.breakpoints[4] = 1.451578e+03; splineFit.breakpoints[5] = 1.456427e+03; splineFit.breakpoints[6] = 1.527625e+03; splineFit.breakpoints[7] = 1.533484e+03; splineFit.breakpoints[8] = 1.536695e+03; splineFit.breakpoints[9] = 1.612799e+03; splineFit.breakpoints[10] = 1.628987e+03; splineFit.breakpoints[11] = 1.629588e+03; splineFit.breakpoints[12] = 1.800000e+03;
  }else if(refmodel == 20120823) {
    splineFit.order = 4;
    splineFit.nbreak = 13;
    splineFit.x_start = 1.300000e+03;
    splineFit.x_end = 1.800000e+03;
    splineFit.coefficients = malloc(15*sizeof(double));
    splineFit.breakpoints = malloc(13*sizeof(double));
    if(splineFit.coefficients == NULL || splineFit.breakpoints == NULL) {
      fflush(stdout);
      printerror(0, "ERROR refmodel_diffgain: Cannot allocate memory for model %ld", refmodel);
      exit(0);
    }
    splineFit.coefficients[0] = 1.090005e+00; splineFit.coefficients[1] = -1.666307e+01; splineFit.coefficients[2] = 2.400624e+01; splineFit.coefficients[3] = -3.882833e+00; splineFit.coefficients[4] = 4.800518e+01; splineFit.coefficients[5] = -4.628771e+01; splineFit.coefficients[6] = 2.453389e+01; splineFit.coefficients[7] = -1.284603e+01; splineFit.coefficients[8] = 5.628643e+00; splineFit.coefficients[9] = -1.199362e+00; splineFit.coefficients[10] = 4.819366e+00; splineFit.coefficients[11] = -1.151882e+01; splineFit.coefficients[12] = 2.561169e+01; splineFit.coefficients[13] = -2.184949e+01; splineFit.coefficients[14] = 2.741658e+00;
    splineFit.breakpoints[0] = 1.300000e+03; splineFit.breakpoints[1] = 1.382188e+03; splineFit.breakpoints[2] = 1.382826e+03; splineFit.breakpoints[3] = 1.384668e+03; splineFit.breakpoints[4] = 1.395310e+03; splineFit.breakpoints[5] = 1.458712e+03; splineFit.breakpoints[6] = 1.458745e+03; splineFit.breakpoints[7] = 1.505330e+03; splineFit.breakpoints[8] = 1.642995e+03; splineFit.breakpoints[9] = 1.651817e+03; splineFit.breakpoints[10] = 1.661119e+03; splineFit.breakpoints[11] = 1.668643e+03; splineFit.breakpoints[12] = 1.800000e+03;
  }else {
    printerror(0, "ERROR refmodel_diffgain: Unknown reference model: %ld", refmodel);
    exit(0);
  }
  ret = cubicBspline_eval_double(splineFit, freq, &ans, verbose);
  cubicBspline_free(&splineFit);
  if(ret != 1) {
    fflush(stdout);
    printerror(0, "ERROR refmodel_diffgain: interpolation failed");
    exit(0);
  }
  return ans;
}

// Shows the available system respondses which can be subtracted from polarization calibration receiver responds
void refmodellist()
{
  int i;
  long id;
  printf("Defined differential phase models:");
  for(i = 1; i < 1000; i++) {
    id = phasemodels_names(i);
    if(i == 1)
      printf(" %ld", id);
    else
      printf(", %ld", id);
  }
  printf("\n");
  printf("Defined differential gain models:");
  for(i = 1; i < 1000; i++) {
    id = diffgainmodels_names(i);
    if(i == 1)
      printf(" %ld", id);
    else
      printf(", %ld", id);
  }
  printf("\n");
}

/* DONE: Next steps: calculate initial guess for gain from averages. */
/* Do amoeba thing to get parameters. I have to make a fake 1 channel/subint wide solution file so I can apply solution. */
/* I guess also long one to store whole solution, so I can write out the pacv style solution when implemented in fits writer */
/* Didn't do rms yet */
/* For highest S/N points, calculate residual: Sum_bins (correcteddata-template)^2/(rmsdata^2+rmstemplate^2). I guess I have to add those chi^2 together? */
/* Maybe it works better when I divide by I in each bin first, so you're less suseptable for variations?*/


/*
  Starting from index 1, return the ID's of the defined phase
  models. Returns -1 if reached the end of the list.
 */
long phasemodels_names(int index)
{
  if(index == 1)
    return 20120126;
  if(index == 2)
    return 20101009;
  if(index == 3)
    return 20090816;
  if(index == 4)
    return 20120819;
  //  if(index == 5)
  //    return ;
  return -1;
}

/*
  Starting from index 1, return the ID's of the defined differential
  gain models. Returns -1 if reached the end of the list.
 */
long diffgainmodels_names(int index)
{
  if(index == 1)
    return 20120126;
  if(index == 2)
    return 20111003;
  if(index == 3)
    return 20111029;
  if(index == 4)
    return 20130529;
  if(index == 5)
    return 20130307;
  if(index == 6)
    return 20120823;
  //  if(index == 7)
  //    return ;
  return -1;
}


/* Given nrsolutions input solutions, construct a single overall
   average receiver solution (combined, not yet initialised with
   cleanPSRData).

   type:
   COMBINERECEIVERMODELS_CHI2: For each channel, store that solution with the smallest chi2
   COMBINERECEIVERMODELS_APPEND: All receiver solutions are combined in a single file with probably multiply defined frequency channels

   Return 1 = success */
int combineReceiverModels(datafile_definition *solution, long nrsolutions, int type, datafile_definition *combined, verbose_definition verbose)
{
  long solnr;
  cleanPSRData(combined, verbose);
  if(nrsolutions < 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR combineReceiverModels: Need at least one input solution.");
    return 0;	
  }
  if(type != COMBINERECEIVERMODELS_CHI2 && type != COMBINERECEIVERMODELS_APPEND) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR combineReceiverModels: Unknown combination algorithm specified.");
    return 0;	
  }
  // Do some sanity checks.
  for(solnr = 0; solnr < nrsolutions; solnr++) {
    if(solution[0].NrSubints != solution[solnr].NrSubints) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR combineReceiverModels: Number of subints mismatch between receiver files that are combined.");
      return 0;
    }
    if(solution[0].NrBins != solution[solnr].NrBins) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR combineReceiverModels: Number of bins mismatch between receiver files that are combined.");
      return 0;
    }
    if(solution[solnr].NrBins < 0 || solution[solnr].NrBins > 2) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR combineReceiverModels: Number of bins was expected to be either 1 or 2.");
      return 0;
    }
    if(solution[0].NrPols != solution[solnr].NrPols) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR combineReceiverModels: Number of receiver parameter (=\"polarizations\") mismatch between receiver files that are combined.");
      return 0;
    }
    
    if(type == COMBINERECEIVERMODELS_CHI2 || type == COMBINERECEIVERMODELS_APPEND) {
      if(solution[0].NrFreqChan != solution[solnr].NrFreqChan) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR combineReceiverModels: Number of subints mismatch between receiver files that are combined.");
	return 0;
      }
    }

    if(type == COMBINERECEIVERMODELS_CHI2) {
      if(solution[solnr].gentype != GENTYPE_RECEIVERMODEL2) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR combineReceiverModels: Cannot find chi^2 values, and therefore the file cannot be combined.");
	return 0;	
      }
      if(solution[solnr].NrPols != 9) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR combineReceiverModels: Expected 9 parameter (\"polarization\") channels. Cannot combine data as chi^2 values cannot be identified.");
	return 0;	
      }
    }
  }

  copy_params_PSRData(solution[0], combined, verbose);
  if(type == COMBINERECEIVERMODELS_CHI2) {
    combined->data = malloc(combined->NrSubints*combined->NrBins*combined->NrPols*combined->NrFreqChan*sizeof(float));
    if(combined->data == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR combineReceiverModels: Memory allocation error.");
      return 0;
    }
  }else if(type == COMBINERECEIVERMODELS_APPEND) {
    combined->NrFreqChan = nrsolutions*solution[0].NrFreqChan;
    //    printf("XXXXX there are %ld freq channels\n", combined->NrFreqChan);
    combined->data = malloc(combined->NrSubints*combined->NrBins*combined->NrPols*combined->NrFreqChan*sizeof(float));
    if(combined->data == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR combineReceiverModels: Memory allocation error.");
      return 0;
    }
    // Allocate memory to store frequency for each solution
    if(combined->freqlabel_list != NULL) {
      free(combined->freqlabel_list);
    }
    combined->freqlabel_list = malloc(combined->NrSubints*combined->NrFreqChan*sizeof(double));
    combined->freqMode = FREQMODE_FREQTABLE;
  }

  if(type == COMBINERECEIVERMODELS_CHI2) {
    for(solnr = 0; solnr < nrsolutions; solnr++) {
      /* Check if this is the first receiver solution to be read in */
      if(solnr == 0) {
	// Start of taking the first solution as the optimum
	memcpy(combined->data, solution[solnr].data, combined->NrSubints*combined->NrBins*combined->NrPols*combined->NrFreqChan*sizeof(float));
      }else {
	long n, b, pol, f;
	float value, value2;
	for(n = 0; n < combined->NrSubints; n++) {
	  for(f = 0; f < combined->NrFreqChan; f++) {
	    /* Get chi^2 */
	    if(readPulsePSRData(&(solution[solnr]), n, 7, f, 0, 1, &value, verbose) != 1) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR combineReceiverModels: Error reading solution.");
	      return 0;
	    }
	    if(readPulsePSRData(combined, n, 7, f, 0, 1, &value2, verbose) != 1) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR combineReceiverModels: Error reading solution.");
	      return 0;
	    }
	    /*	    fprintf(stderr, "XXXXX %f %f\n", value, value2); */
	    // Keep largest chi2 solution
	    if(value > 0 && (value < value2 || value2 < 0)) {
	      for(pol = 0; pol < combined->NrPols; pol++) {
		for(b = 0; b < combined->NrBins; b++) {
		  if(readPulsePSRData(&(solution[solnr]), n, pol, f, b, 1, &value, verbose) != 1) {
		    fflush(stdout);
		    printerror(verbose.debug, "ERROR combineReceiverModels: Error reading solution.");
		    return 0;
		  }
		  if(writePulsePSRData(combined, n, pol, f, b, 1, &value, verbose) != 1) {
		    fflush(stdout);
		    printerror(verbose.debug, "ERROR combineReceiverModels: Error writing solution.");
		    return 0;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }else if(type == COMBINERECEIVERMODELS_APPEND) {
    long freq_in_output = 0;
    for(solnr = 0; solnr < nrsolutions; solnr++) {
      long n, b, pol, f;
      float value;
      for(f = 0; f < solution[solnr].NrFreqChan; f++) {
	for(n = 0; n < combined->NrSubints; n++) {
	  double freqvalue;
	  freqvalue = get_weighted_channel_freq(solution[solnr], n, f, verbose);
	  if(set_weighted_channel_freq(combined, n, freq_in_output, freqvalue, verbose) == 0) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR combineReceiverModels: Error setting frequency.");
	    return 0;
	  }
	  for(pol = 0; pol < combined->NrPols; pol++) {
	    for(b = 0; b < combined->NrBins; b++) {
	      if(readPulsePSRData(&(solution[solnr]), n, pol, f, b, 1, &value, verbose) != 1) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR combineReceiverModels: Error reading solution.");
		return 0;
	      }
	      if(writePulsePSRData(combined, n, pol, freq_in_output, b, 1, &value, verbose) != 1) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR combineReceiverModels: Error writing solution.");
		return 0;
	      }
	    }
	  }
	}
	freq_in_output++;
      }
    }
  }
  return 1;
}



/* Take the input receiver model and replace the parameters with the
   analytic solution.

   Return 0 = Error, 1 = success 
*/
int makeReceiverModelAnalytic(datafile_definition *solution, analyticReceiverSolution_def *analyticsolution, verbose_definition verbose)
{
  long i, subintnr, freqnr, parameter;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Make receiver model analytic\n");
  }
  if(solution->gentype != GENTYPE_RECEIVERMODEL2 && solution->gentype != GENTYPE_RECEIVERMODEL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR makeReceiverModelAnalytic: File does not appear to be a receiver solution->");
    return 0;
  }
  if(solution->NrPols < 3) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR makeReceiverModelAnalytic: Expecteded at least 3 receiver parameters, only %ld present.", solution->NrPols);
    return 0;
  }
  if(solution->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR makeReceiverModelAnalytic: receiver model should be read into memory first");
    return 0;
  }


  // Check the frequency range in use: i.e. ignore band edges.
  // Should test if all parameters are defined (finite, and != nan)
  long freqmin = -1;
  for(freqnr = 0; freqnr < solution->NrFreqChan; freqnr++) {
    for(subintnr = 0; subintnr < solution->NrSubints; subintnr++) {
      parameter = 0;
      float gain = solution->data[solution->NrBins*(parameter+solution->NrPols*(freqnr+subintnr*solution->NrFreqChan))+0];
      if(gain > 0) {
	freqmin = freqnr;
	break;
      }
    }
    if(freqmin != -1)
      break;
  }
  long freqmax = -1;
  for(freqnr = solution->NrFreqChan-1; freqnr >= 0; freqnr--) {
    for(subintnr = 0; subintnr < solution->NrSubints; subintnr++) {
      parameter = 0;
      float gain = solution->data[solution->NrBins*(parameter+solution->NrPols*(freqnr+subintnr*solution->NrFreqChan))+0];
      if(gain > 0) {
	freqmax = freqnr;
	break;
      }
    }
    if(freqmax != -1)
      break;
  }

  for(parameter = 0; parameter < solution->NrPols; parameter++) {
    if((parameter == 0 && analyticsolution->gainmodel_defined)
       || (parameter == 1 && analyticsolution->diffgainmodel_defined)
       || (parameter == 2 && analyticsolution->phasemodel_defined)
       || (parameter == 3 && analyticsolution->ell1model_defined)
       || (parameter == 5 && analyticsolution->ell2model_defined)
       || (parameter == 4 && analyticsolution->or1model_defined)
       || (parameter == 6 && analyticsolution->or2model_defined)) {

      /* Loop over pulses and frequency channels etc. */
      for(subintnr = 0; subintnr < solution->NrSubints; subintnr++) {
	for(freqnr = 0; freqnr < solution->NrFreqChan; freqnr++) {
	  double freq, value;
	  freq = get_weighted_channel_freq(*solution, subintnr, freqnr, verbose);
	  value = 0;

	  if(parameter == 0) {
	    if(analyticReceiverSolution_gain(*analyticsolution, freq, &value, verbose) == 0) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR makeReceiverModelAnalytic: Model calculation failed.");
	      return 0;
	    }
	  }else if(parameter == 1) {
	    if(analyticReceiverSolution_diffgain(*analyticsolution, freq, &value, verbose) == 0) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR makeReceiverModelAnalytic: Model calculation failed.");
	      return 0;
	    }
	  }else if(parameter == 2) {
	    if(analyticReceiverSolution_diffphase(*analyticsolution, freq, &value, verbose) == 0) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR makeReceiverModelAnalytic: Model calculation failed.");
	      return 0;
	    }
	  }else if(parameter == 3 || parameter == 5) {
	    int polnr = 1;
	    if(parameter == 5)
	      polnr = 2;
	    if(analyticReceiverSolution_ellipticity(*analyticsolution, polnr, freq, &value) == 0) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR makeReceiverModelAnalytic: Model calculation failed.");
	      return 0;
	    }
	  }else if(parameter == 4 || parameter == 6) {
	    int polnr = 1;
	    if(parameter == 6)
	      polnr = 2;
	    if(analyticReceiverSolution_orientation(*analyticsolution, polnr, freq, &value) == 0) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR makeReceiverModelAnalytic: Model calculation failed.");
	      return 0;
	    }
	  }

	  float yerrorbar = 1e-5;
	  if(parameter == 1) {
	    //	    (*ydata)[f] = 100.0*(exp(2.0*(*ydata)[f])-1.0);
	    //	    actual = 100.0*(exp(2.0*stored)-1.0);
	    //	    0.01*actual + 1.0 = exp(2.0*stored)
	    //	    0.5*ln(0.01*actual + 1.0) = stored
	    value = 0.5*log(0.01*value + 1.0);
	    yerrorbar = 0;
	  }else if(parameter >= 2 && parameter <= 6) {
	    value *= M_PI/180.0;
	    yerrorbar = -1.0*M_PI/180.0;
	  }else if(parameter == 7 || parameter == 8) {
	    yerrorbar = 0;
	  }

	  if(freqnr < freqmin || freqnr > freqmax) {
	    yerrorbar = -1;
	    if(parameter == 1) {
	      value = -1;
	    }else {
	      value = sqrt(-1);
	    }
	  }

	  //	  printf("Replacing %e with %e\n", solution->data[solution->NrBins*(parameter+solution->NrPols*(freqnr+subintnr*solution->NrFreqChan))+0], value);
	  solution->data[solution->NrBins*(parameter+solution->NrPols*(freqnr+subintnr*solution->NrFreqChan))+0] = value;
	  if(solution->NrBins == 2) {
	    solution->data[solution->NrBins*(parameter+solution->NrPols*(freqnr+subintnr*solution->NrFreqChan))+1] = yerrorbar;
	  }
	  
	}
      }
    }
  }
  return 1;
}

