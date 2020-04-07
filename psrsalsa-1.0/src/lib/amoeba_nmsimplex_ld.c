/* This is an exact copy of amoeba_nmsimplex.c (float version) to make
   the long double version. This file should be regeneratable by
   applying the following substitutions:

   float     ->  long double 
   _nmsimplex -> _nmsimplex_ld
   %f        -> %Lf
   %e        -> %Le
   fabs      -> fabsl
   sqrt      -> sqrtl
*/

/* 
 * Program: nmsimplex.c
 * Author : Michael F. Hutt
 * http://www.mikehutt.com
 * 11/3/97
 *
 * An implementation of the Nelder-Mead simplex method.
 *
 * Copyright (c) 1997-2011 <Michael F. Hutt>
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *
 * Jan. 6, 1999 
 * Modified to conform to the algorithm presented
 * in Margaret H. Wright's paper on Direct Search Methods.
 *
 * Jul. 23, 2007
 * Fixed memory leak.
 *
 * Mar. 1, 2011
 * Added constraints.
 *
 * Patrick Weltevrede added the following functionality:
 * A step size should be provided, rather than a scale
 * A verbose option is added and output format is different
 * Removed the constrain option
 * Added the reachedEpsilon parameter
 * Changed convergence condition
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define MAX_IT      5000      /* maximum number of iterations */
#define ALPHA       1.0       /* reflection coefficient */
#define BETA        0.5       /* contraction coefficient */
#define GAMMA       2.0       /* expansion coefficient */

long double amoeba_nmsimplex_ld(long double (*objfunc)(long double[]), long double start[], long double dx[], int n, long double EPSILON, int *nfunc, long double *reachedEpsilon, int verbose)
{
  //  void (*constrain)(long double[],int n);

  int vs;         /* vertex with smallest value */
  int vh;         /* vertex with next smallest value */
  int vg;         /* vertex with largest value */
  
  int i,j,m,row;
  int k;   	      /* track the number of function evaluations */
  
  long double **v;     /* holds vertices of simplex */
  //	long double pn,qn;   /* values used to create initial simplex */
  long double *f;      /* value of function at each vertex */
  long double fr;      /* value of function at reflection point */
  long double fe;      /* value of function at expansion point */
  long double fc;      /* value of function at contraction point */
  long double *vr;     /* reflection - coordinates */
  long double *ve;     /* expansion - coordinates */
  long double *vc;     /* contraction - coordinates */
  long double *vm;     /* centroid - coordinates */
  long double min;
	
  long double s,cent;
  //  long double fsum,favg;

  /* dynamically allocate arrays */
  
  /* allocate the rows of the arrays */
  v =  (long double **) malloc ((n+1) * sizeof(long double *));
  f =  (long double *) malloc ((n+1) * sizeof(long double));
  vr = (long double *) malloc (n * sizeof(long double));
  ve = (long double *) malloc (n * sizeof(long double));  
  vc = (long double *) malloc (n * sizeof(long double));  
  vm = (long double *) malloc (n * sizeof(long double));  
	
  /* allocate the columns of the arrays */
  for (i=0;i<=n;i++) {
    v[i] = (long double *) malloc (n * sizeof(long double));
  }
  
  /* create the initial simplex */
  /* assume one of the vertices is 0,0 */
  /*	
	pn = scale*(sqrtl(n+1)-1+n)/(n*sqrtl(2));
	qn = scale*(sqrtl(n+1)-1)/(n*sqrtl(2));
  */
  
  for (i=0;i<n;i++) {
    v[0][i] = start[i];
  }
  
  for (i=1;i<=n;i++) {
    for (j=0;j<n;j++) {
      if (i-1 == j) {
	//				v[i][j] = pn + start[j];
	v[i][j] = start[j]+dx[j];
      }
      else {
	//				v[i][j] = qn + start[j];
	v[i][j] = start[j];
      }
    }
  }
  
  //	if (constrain != NULL) {
  //    constrain(v[j],n);
  //  } 
  
  /* find the initial function values */
  for (j=0;j<=n;j++) {
    f[j] = objfunc(v[j]);
  }
  
  k = n+1;
  
  /* print out the initial values */
  if(verbose) {
    printf("Initial Values\n");
    for (j=0;j<=n;j++) {
      printf("Vertex %d:", j+1);
      for (i=0;i<n;i++) {
	printf(" %Lf",v[j][i]);
      }
      printf(": %Lf\n",f[j]);
    }
  }
	
	
  /* begin the main loop of the minimization */
  int nritterations;
  for (nritterations=1;nritterations<=MAX_IT;(nritterations)++) {     
    /* find the index of the largest value */
    vg=0;
    for (j=0;j<=n;j++) {
      if (f[j] > f[vg]) {
	vg = j;
      }
    }
    
    /* find the index of the smallest value */
    vs=0;
    for (j=0;j<=n;j++) {
      if (f[j] < f[vs]) {
	vs = j;
      }
    }
    
    /* find the index of the second largest value */
    vh=vs;
    for (j=0;j<=n;j++) {
      if (f[j] > f[vh] && f[j] < f[vg]) {
	vh = j;
      }
    }
    
    /* calculate the centroid */
    for (j=0;j<=n-1;j++) {
      cent=0.0;
      for (m=0;m<=n;m++) {
	if (m!=vg) {
	  cent += v[m][j];
	}
      }
      vm[j] = cent/n;
    }
    
    /* reflect vg to new vertex vr */
    for (j=0;j<=n-1;j++) {
      /*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
      vr[j] = vm[j]+ALPHA*(vm[j]-v[vg][j]);
    }
    //		if (constrain != NULL) {
    //      constrain(vr,n);
    //    }
    fr = objfunc(vr);
    k++;
    
    if (fr < f[vh] && fr >= f[vs]) {
      for (j=0;j<=n-1;j++) {
	v[vg][j] = vr[j];
      }
      f[vg] = fr;
    }
    
    /* investigate a step further in this direction */
    if ( fr <  f[vs]) {
      for (j=0;j<=n-1;j++) {
	/*ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];*/
	ve[j] = vm[j]+GAMMA*(vr[j]-vm[j]);
      }
      //			if (constrain != NULL) {
      //        constrain(ve,n);
      //      }
      fe = objfunc(ve);
      k++;
      
      /* by making fe < fr as opposed to fe < f[vs], 			   
	 Rosenbrocks function takes 63 iterations as opposed 
	 to 64 when using long double variables. */
      
      if (fe < fr) {
	for (j=0;j<=n-1;j++) {
	  v[vg][j] = ve[j];
	}
	f[vg] = fe;
      }
      else {
	for (j=0;j<=n-1;j++) {
	  v[vg][j] = vr[j];
	}
	f[vg] = fr;
      }
    }
    
    /* check to see if a contraction is necessary */
    if (fr >= f[vh]) {
      if (fr < f[vg] && fr >= f[vh]) {
	/* perform outside contraction */
	for (j=0;j<=n-1;j++) {
	  /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
	  vc[j] = vm[j]+BETA*(vr[j]-vm[j]);
	}
	//				if (constrain != NULL) {
	//          constrain(vc,n);
	//        }
	fc = objfunc(vc);
	k++;
      }
      else {
	/* perform inside contraction */
	for (j=0;j<=n-1;j++) {
	  /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
	  vc[j] = vm[j]-BETA*(vm[j]-v[vg][j]);
	}
	//				if (constrain != NULL) {
	//          constrain(vc,n);
	//        }
	fc = objfunc(vc);
	k++;
      }
      
      
      if (fc < f[vg]) {
	for (j=0;j<=n-1;j++) {
	  v[vg][j] = vc[j];
	}
	f[vg] = fc;
      }
      /* at this point the contraction is not successful,
	 we must halve the distance from vs to all the 
	 vertices of the simplex and then continue.
	 10/31/97 - modified to account for ALL vertices. 
      */
      else {
	for (row=0;row<=n;row++) {
	  if (row != vs) {
	    for (j=0;j<=n-1;j++) {
	      v[row][j] = v[vs][j]+(v[row][j]-v[vs][j])/2.0;
	    }
	  }
	}
	//				if (constrain != NULL) {
	//          constrain(v[vg],n);
	//        }
	f[vg] = objfunc(v[vg]);
	k++;
	//				if (constrain != NULL) {
	//          constrain(v[vh],n);
	//        }
	f[vh] = objfunc(v[vh]);
	k++;
	
	
      }
    }
    
    /* print out the value at each iteration */
    if(verbose) {
      printf("Iteration %d\n",nritterations);
      for (j=0;j<=n;j++) {
	printf("Vertex %d:", j+1);
	for (i=0;i<n;i++) {
	  printf(" %Lf",v[j][i]);
	}
	printf(": %Lf\n",f[j]);
      }
    }
    
    /* test for convergence */
    /*
    fsum = 0.0;
    for (j=0;j<=n;j++) {
      fsum += f[j];
    }
    favg = fsum/(n+1);
    s = 0.0;
    for (j=0;j<=n;j++) {
      s += pow((f[j]-favg),2.0)/(n);
    }
    s = sqrtl(s);*/
    long double yhi, ylo;

    yhi = ylo = f[0];
    for (j=0;j<=n;j++) {
      if(f[j] > yhi)
	yhi = f[j];
      if(f[j] < ylo)
	ylo = f[j];
    }
    s = 2*fabsl(yhi-ylo)/(fabsl(yhi)+fabsl(ylo));

    if (s < EPSILON) break;
  }
  /* end main loop of the minimization */
  
  *reachedEpsilon = s;

  /* find the index of the smallest value */
  vs=0;
  for (j=0;j<=n;j++) {
    if (f[j] < f[vs]) {
      vs = j;
    }
  }
  
  if(verbose) {
    printf("The minimum was found at epsilon=%Le:\n  ", *reachedEpsilon); 
    for (j=0;j<n;j++) {
      printf(" %Le",v[vs][j]);
    }
    printf("\n");
  }
  for (j=0;j<n;j++) {
    start[j] = v[vs][j];
  }
  min=objfunc(v[vs]);
  k++;
  if(verbose) {
    printf("%d Function Evaluations\n",k);
    printf("%d Iterations through program\n",nritterations);
  }
  *nfunc = k;

  free(f);
  free(vr);
  free(ve);
  free(vc);
  free(vm);
  for (i=0;i<=n;i++) {
    free (v[i]);
  }
  free(v);
  return min;
}
