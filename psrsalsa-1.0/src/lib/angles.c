//START REGION RELEASE
#include <math.h>
//START REGION DEVELOP
#include <stdio.h>
#include <string.h>
#include "psrsalsa.h"


//START REGION RELEASE

/* Derotates an angle between 0 and 360 degrees */
float derotate_deg(float a)
{
  int i;
  i = fabs(a)/360.0;
  if(a > 0)
    a -= 360*i;
  else
    a += 360*i;
  if(a < 0)
    a += 360;
  return a;
}

//START REGION DEVELOP

/* Derotates an angle between 0 and 360 degrees */
double derotate_deg_double(double a)
{
  int i;
  i = fabs(a)/360.0;
  if(a > 0)
    a -= 360*i;
  else
    a += 360*i;
  if(a < 0)
    a += 360;
  return a;
}
/* Derotates an angle between 0 and 360 degrees */
long double derotate_deg_longdouble(long double a)
{
  int i;
  i = fabsl(a)/360.0;
  if(a > 0)
    a -= 360*i;
  else
    a += 360*i;
  if(a < 0)
    a += 360;
  return a;
}

/* Derotates an angle between -180 and 180 degrees */
double derotate_deg_small_double(double a)
{
  double y;
  y = derotate_deg_double(a);
  if(y >= 180.0)
    y = y-360.0;
  return y;
}

/* Derotates an angle between 0 and 2 pi radians */
long double derotate_rad_longdouble(long double a)
{
  int i;
  i = fabsl(a)/(2.0*M_PI);
  if(a > 0)
    a -= 2.0*M_PI*i;
  else
    a += 2.0*M_PI*i;
  if(a < 0)
    a += 2.0*M_PI;
  return a;
}

/* Derotates an angle between 0 and 2 pi radians */
double derotate_rad_double(double a)
{
  int i;
  i = fabs(a)/(2.0*M_PI);
  if(a > 0)
    a -= 2.0*M_PI*i;
  else
    a += 2.0*M_PI*i;
  if(a < 0)
    a += 2.0*M_PI;
  return a;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Derotates an angle between 0 and 180 degrees */
float derotate_180(float a)
{
  int i;
  i = fabs(a)/180.0;
  if(a > 0)
    a -= 180*i;
  else
    a += 180*i;
  if(a < 0)
    a += 180;
  return a;
}


/* Derotates an angle between 0 and 180 degrees */
double derotate_180_double(double a)
{
  int i;
  i = fabs(a)/180.0;
  if(a > 0)
    a -= 180*i;
  else
    a += 180*i;
  if(a < 0)
    a += 180;
  return a;
}


/* Derotates an angle between 0 and pi radians */
double derotate_180_rad_double(double a)
{
  int i;
  i = fabs(a)/M_PI;
  if(a > 0)
    a -= M_PI*i;
  else
    a += M_PI*i;
  if(a < 0)
    a += M_PI;
  return a;
}

/* Derotates an angle between 0 and 90 degrees */
float derotate_90(float a)
{
  int i;
  i = fabs(a)/90.0;
  if(a > 0)
    a -= 90*i;
  else
    a += 90*i;
  if(a < 0)
    a += 90;
  return a;
}


/* Derotates an angle between 0 and 90 degrees */
double derotate_90_double(double a)
{
  int i;
  i = fabs(a)/90.0;
  if(a > 0)
    a -= 90*i;
  else
    a += 90*i;
  if(a < 0)
    a += 90;
  return a;
}


//START REGION DEVELOP
//START REGION RELEASE

/* Derotates an angle between -90 and 90 degrees */
double derotate_180_small_double(double a)
{
  double y;
  y = derotate_180_double(a);
  if(y >= 90.0)
    y = y-180.0;
  return y;
}

// Derotates an angle between -45 and 45 degrees
double derotate_90_small_double(double a)
{
  double y;
  y = derotate_90_double(a);
  if(y >= 45.0)
    y = y-90.0;
  return y;
}


//START REGION DEVELOP
//START REGION RELEASE

/* Calculates the angle of the vector (x,y) in radians. (1,0)
   corresponds with 0 rad and counter clockwise corresponds to an
   increasing angle. */
float polar_angle_rad(float x, float y)
{
  float alpha;
  if(x == 0) {
    if(y <= 0)
      return 1.5*M_PI;
    else
      return 0.5*M_PI;
  }else {
    alpha = atan(y/x);
    if(x > 0 && y >= 0)
      return alpha;
    if(x > 0)
      return alpha + 2*M_PI;
    return alpha + M_PI;
  }
}

//START REGION DEVELOP

double polar_angle_rad_d(double x, double y)
{
  double alpha;
  if(x == 0) {
    if(y <= 0)
      return 1.5*M_PI;
    else
      return 0.5*M_PI;
  }else {
    alpha = atan(y/x);
    if(x > 0 && y >= 0)
      return alpha;
    if(x > 0)
      return alpha + 2*M_PI;
    return alpha + M_PI;
  }
}

/*
  Calculates the azimuth and zenith angle for a source at a given lst
  and latitude. All arguments are angles in radians.
 */
void calculate_az_zen(double raj, double decj, double lst, double latitude, double *az, double *zen)
{
  double x, y;
  *zen = acos(sin(latitude)*sin(decj)+cos(latitude)*cos(decj)*cos(lst-raj));
  x = -sin(lst-raj);
  y = cos(latitude)*tan(decj)-sin(latitude)*cos(lst-raj);
  *az    = polar_angle_rad(y,x);
}


/* For a given right ascention and declination calculate the rise and
   set times (LST) for a given elevation limit. All parameters are
   angles in radians).

   Return values:
     0 - Never visible
     1 - Does rise and set
     2 - Always visible
 */
int riseset(double alpha, double delta, double elevation_limit, double latitude, double *rise, double *set, int verbose)
{
  double u1, u2, u, t1, t2, s1, s2;
  int ret;
  u2 = sin(elevation_limit)-sin(delta)*sin(latitude);
  u2 /= cos(delta)*cos(latitude);

  t1 = alpha*12/M_PI - 6;
  if(t1 >= 24)
    t1 -= 24;
  if(t1 < 0)
    t1 += 24;
  t2 = alpha*12/M_PI + 6;
  if(t2 > 24)
    t2 -= 24;
  if(t2 < 0)
    t2 += 24;

  if(t2 < t1)
    t2 += 24;

  if(u2 > 1 || u2 < -1) {
    if(verbose) printf("\nAbove elevation limit: ");
    if(sin(delta)*sin(latitude)+cos(delta)*cos(latitude) > sin(elevation_limit)) {
      if(verbose) printf("Source does not set.\n");
      ret = 2;
      s1 = -1;
      s2 = 25;
    }else {
      if(verbose) printf("Source does not rise.\n");
      s1 = 25;
      s2 = -1;
      ret = 0;
    }

    /*
    printf("in HA range                       : ");
    printhour(t1*M_PI/12.0);
    printf(" - ");
    printhour(t2*M_PI/12.0);
    printf(" LST\n");
    s1 = t1*M_PI/12.0;
    s2 = t2*M_PI/12.0;

    printf("Visible with WSRT                 : ");
    printhour(s1);
    printf(" - ");
    printhour(s2);
    printf(" LST\n");
    */
  }else {
    ret = 1;
    u1 = acos(u2);
    u2 = 2*M_PI - u1;

    if(sin(delta)*sin(latitude)+cos(delta)*cos(latitude)*cos(u1+0.01) < sin(elevation_limit)) {
      u = u2;
      u2 = u1;
      u1 = u;
    }
    
    s1 = u1 + alpha;
    s2 = u2 + alpha;

    if(s1 > 2*M_PI)
      s1 -= 2*M_PI;
    if(s1 < 0)
      s1 += 2*M_PI;
    if(s2 > 2*M_PI)
      s2 -= 2*M_PI;
    if(s2 < 0)
      s2 += 2*M_PI;
    
    /*
    if(s2 < s1)
      s2 += 2*M_PI;
    */
    if(verbose) {
      printf("s1 = %lf, s2 = %lf\n", 12*s1/M_PI, 12*s2/M_PI);
      printf("\nAbove elevation limit: ");
      printhour(s1);
      printf(" - ");
      printhour(s2);
      printf(" LST\n");

      /*
      printf("in HA range                       : ");
      printhour(t1*M_PI/12.0);
      printf(" - ");
      printhour(t2*M_PI/12.0);
      printf(" LST\n");
      if(s1 < t1*M_PI/12.0)
	s1 = t1*M_PI/12.0;
      if(s2 > t2*M_PI/12.0)
	s2 = t2*M_PI/12.0;
      printf("Visible with WSRT                 : ");
      printhour(s1);
      printf(" - ");
      printhour(s2);
      printf(" LST\n");
      */
    }
  }
  *rise = s1;
  *set = s2;
  return ret;
}

/* Get the position of the vector from the pulsar to Earth with
   respect to the magnetic axis (in polar coordinates, radians) at
   time t for a pulse period period (should have the same
   units). alpha is the angle between the magnetic axis and the
   rototation axis and zeta between the rotation axis and the line of
   sight. All angles are in radians. */

void pulsar_polar_coordinates(float t, float period, float alpha, float zeta, float *phi, float *theta)
{
  float y = sin(2*M_PI*t/period);
  float x = cos(2*M_PI*t/period)*sin(0.5*M_PI - alpha) - 
            tan(0.5*M_PI - zeta)*cos(0.5*M_PI - alpha);
  *phi = polar_angle_rad(x, y);
  *theta = 0.5*M_PI - (asin(sin(0.5*M_PI - alpha)*sin(0.5*M_PI - zeta) +
                                 cos(0.5*M_PI - alpha)*cos(0.5*M_PI -
				 zeta)*cos(2*M_PI*t/period)));
}

void pulsar_polar_coordinates_d(double t, double period, double alpha, double zeta, double *phi, double *theta)
{
  double y = sin(2*M_PI*t/period);
  double x = cos(2*M_PI*t/period)*sin(0.5*M_PI - alpha) - 
            tan(0.5*M_PI - zeta)*cos(0.5*M_PI - alpha);
  *phi = polar_angle_rad_d(x, y);
  *theta = 0.5*M_PI - (asin(sin(0.5*M_PI - alpha)*sin(0.5*M_PI - zeta) +
                                 cos(0.5*M_PI - alpha)*cos(0.5*M_PI -
				 zeta)*cos(2*M_PI*t/period)));
}

/* If rmin or rmax are negative it will calculate the whole field line
   instead of section. theta0 (in deg) is the footprint polar angle of
   the field line from the magnetic axis. alpha (in deg) is the angle
   between the rotation axis (positive z-axis) and the magnetic
   axis. Set s to +1 or -1 to indicate which side of the polar cap
   you're considering. */
void fieldLine(float theta0, float alpha, int s, float *x, float *z, long nrpoints, float Rns, float rmin, float rmax)
{
  /* For a dipole fieldline sin(theta)^2/r is constant */

  long n;
  float theta, r, c, theta1, theta2;
  theta0 *= M_PI/180.0;
  c = sin(theta0)*sin(theta0)/Rns;
  if(rmin < 0 || rmax < 0) {
    theta1 = theta0;
    theta2 = M_PI-theta0;
  }else {
    theta1 = asin(sqrt(c*rmin));
    theta2 = asin(sqrt(c*rmax));
  }
  for(n = 0; n < nrpoints; n++) {
    if(nrpoints == 1)
      theta = theta1;
    else
      theta = theta1+(theta2-theta1)*n/(float)(nrpoints-1);
    r = sin(theta)*sin(theta)/c;
    x[n] = r*sin(s*theta+M_PI*alpha/180.0);
    z[n] = r*cos(s*theta+M_PI*alpha/180.0);
  }
}

void fieldLine_d(double theta0, double alpha, int s, double *x, double *z, long nrpoints, double Rns, double rmin, double rmax)
{
  long n;
  double theta, r, c, theta1, theta2;
  theta0 *= M_PI/180.0;
  c = sin(theta0)*sin(theta0)/Rns;
  if(rmin < 0 || rmax < 0) {
    theta1 = theta0;
    theta2 = M_PI-theta0;
  }else {
    theta1 = asin(sqrt(c*rmin));
    theta2 = asin(sqrt(c*rmax));
  }
  for(n = 0; n < nrpoints; n++) {
    if(nrpoints == 1)
      theta = theta1;
    else
      theta = theta1+(theta2-theta1)*n/(double)(nrpoints-1);
    r = sin(theta)*sin(theta)/c;
    x[n] = r*sin(s*theta+M_PI*alpha/180.0);
    z[n] = r*cos(s*theta+M_PI*alpha/180.0);
  }
}

/* This is equation 39 of Michel and Li 1999, Physics Reports 318, 227
This calculates the field strength at a position indicated by the
three spherical coordinates (r, theta, phi), normalized by the field
strength at the equator of the star. The neutron star radius must be
profided in meters, as well as the rotation period P. All angles are
in radians and all distances in meters. theta is the polar angle
measured from the rotation axis and alpha is the angle between the
rotation axis and the magnetic axis.
*/
void fieldStrength_Deutsch_d(double r, double theta, double phi, double alpha, double Rns, double P, double *Br, double *Btheta, double *Bphi)
{
  double psi, zeta, zeta2, rz, rho2, zeta4, rho, d1, d2, d3, d4, q1, q2, constant;
  double cost, sint, cosa, sina, cospsi, sinpsi;

  constant = 2.0*M_PI/(299792458.0*P);
  zeta = Rns*constant;
  rho  = r*constant;
  psi = phi+rho-zeta;


  zeta2 = zeta*zeta;
  zeta4 = zeta2*zeta2;
  rho2 = rho*rho;
  rz = rho*zeta;

  constant = 1.0/(zeta2+1.0);
  d1 = (rz+1.0)*constant;
  d2 = (rho - zeta)*constant;
  d3 = (1.0+rz-rho2)*constant;
  d4 = ((rho2 - 1.0)*zeta+rho)*constant;


  constant = zeta4*(zeta2-3.0) + 36.0;
  q1  = 3.0*rz*(6.0*zeta2-zeta4)+3.0*(3.0-rho2)*zeta2*(2.0-zeta2);
  q1 /= constant;

  q2  = (3.0-rho2)*zeta*(zeta4-6.0*zeta2)+9.0*rho*zeta2*(2.0-zeta2);
  q2 /= constant;


  cosa = cos(alpha);
  sina = sin(alpha);
  cost = cos(theta);
  sint = sin(theta);
  cospsi = cos(psi);
  sinpsi = sin(psi);
  constant = Rns*Rns*Rns/(r*r*r);
  *Br = 2.0*constant*(cosa*cost+sina*sint*(d1*cospsi+d2*sinpsi));
  *Btheta = constant*(cosa*sint-sina*cost*((q1+d3)*cospsi+(q2+d4)*sinpsi));
  *Bphi = constant*sina*(-1.0*(q2*cos(2.0*theta)+d4)*cospsi + (q1*cos(2.0*theta)+d3)*sinpsi);

  /* Between brackets is the predicted value in the near field */
  /*  fprintf(stderr, "d1=%e (1) d2=%e (0) d3=%e (1) d4=%e (0)\n", d1, d2, d3, d4);
      fprintf(stderr, "q1=%e (%e) q2=%e (0)\n", q1, zeta2/2.0, q2);*/
}


/* Converts the (vector dr, dtheta, dphi) at position (theta, phi)
   into the Cartesian vector (dx, dy, dz) */
void dSpherical2dCartesian(double theta, double phi, double dr, double dtheta, double dphi, double *dx, double *dy, double *dz)
{
  *dx = (dr*sin(theta)*cos(phi) + dtheta*cos(theta)*cos(phi) - dphi*sin(phi));
  *dy = (dr*sin(theta)*sin(phi) + dtheta*cos(theta)*sin(phi) + dphi*cos(phi));
  *dz = (dr*cos(theta) - dtheta*sin(theta));
}

/* theta0 and phi0 (in deg) are the polar (from the magnetic axis) and
   azimuthal angle of the field line at rmin. alpha (in deg) is the
   angle between the rotation axis (positive z-axis) and the magnetic
   axis. All distances should be in meters and the period P in
   seconds. If rmin is negative, it is considered to be the neutron
   star radius. If rmax is negative, it will be set to 3 times the
   light cylinder radius. */
void fieldLine_Deutsch_d_old(double theta0, double phi0, double alpha, double *x, double *y, double *z, long nrpoints, double Rns, double P, double rmin, double rmax) 
{
  /* This is the number of points along the field line that is stored,
     which will resampled to nrpoints. If this number is too small,
     the field line has to be calculated twice. Doesn't really appear
     to speed up things. */
  //  #define fieldLine_Deutsch_d_maxBufferSize  100000
  //  double xbuffer[fieldLine_Deutsch_d_maxBufferSize], ybuffer[fieldLine_Deutsch_d_maxBufferSize], zbuffer[fieldLine_Deutsch_d_maxBufferSize];
  //  long nrbufferpoints, nrbufferpointsskip;

  double r, theta, phi, Br, Btheta, Bphi, Bx, By, Bz, length, factor, sign, xpos, ypos, zpos, totlength, totlengthsave;
  double Bxold, Byold, Bzold, improd;
  long i;
  int firsttime, step;

  /* By default, start at the NS surface */
  if(rmin < 0)  
    rmin = Rns;
  if(rmax < 0)
    rmax = 3.0*299792458*P/(2.0*M_PI);
  
  //  nrbufferpoints = 0;
  //  nrbufferpointsskip = 0;

  /* First step: calculate total length of field line piece.
     Second step: save the requested number of points along field line */
  fprintf(stderr, "Input: r=%e theta=%.2lf phi=%.2lf nrpoints=%ld\n", rmin, theta0, phi0, nrpoints); 
  for(step = 1; step <= 2; step++) { 
    r = rmin;
    theta = M_PI*theta0/180.0;
    phi = M_PI*phi0/180.0;

    /* Do a coordinate transformation from magnetic coordinates to rotational coordinates. */
    double p, t;
    pulsar_polar_coordinates_d(phi, 2.0*M_PI, alpha*M_PI/180.0, theta, &p, &t);
    phi = p+M_PI;
    theta = t;

    /* So z-axis is at theta=0, y-axis is at phi=90 and x-axis at phi=0 */
    xpos = r*sin(theta)*cos(phi);
    ypos = r*sin(theta)*sin(phi);
    zpos = r*cos(theta);
    //    fprintf(stderr, "Starting with:  r=%e theta=%e phi=%e\n", r, theta*180.0/M_PI, phi*180.0/M_PI); 
    /*    printf("%e %e %e\n", xpos, ypos, zpos); */
    //    exit(0);

    if(step == 2) {
      x[0] = xpos;
      y[0] = ypos;
      z[0] = zpos;
    }
    //else if(step == 1) {
      //      xbuffer[0] = xpos;
    //      ybuffer[0] = ypos;
    //      zbuffer[0] = zpos;
    //      nrbufferpoints++;
    //      nrbufferpointsskip++;
    //    }
    firsttime = 1;
    if(step == 2) {
      totlengthsave = totlength;
      i = 1;
    }
    totlength = 0;
    //    fprintf(stderr, "Check: r=%e rmin=%e rmax=%e\n", r, rmin, rmax);
    while(r >= rmin && (r <= rmax || rmax < 0)) {
      //      if(phi0 == 0 && theta0 == 270)
	//	fprintf(stderr, "calculate r=%e/%.1e theta=%.2lf phi=%.2lf alpha=%.2f\n", r, rmax, theta*180/M_PI, phi*180/M_PI, alpha);
      fieldStrength_Deutsch_d(r, theta, phi, alpha*M_PI/180.0, Rns, P, &Br, &Btheta, &Bphi);
      if(firsttime == 1) {
	if(Br >= 0) {
	  sign = 1;
	  //	  fprintf(stderr, "Br=%e, so choose sign %f\n", Br, sign);
	}else {
	  sign = -1;
	  //	  fprintf(stderr, "Br=%e, so choose sign %f\n", Br, sign);
	}
      }
      dSpherical2dCartesian(theta, phi, Br, Btheta, Bphi, &Bx, &By, &Bz);
      Bx *= sign;
      By *= sign;
      Bz *= sign;

      /*
      Bx = sign*(Br*sin(theta)*cos(phi) + Btheta*cos(theta)*cos(phi) - Bphi*sin(phi));
      By = sign*(Br*sin(theta)*sin(phi) + Btheta*cos(theta)*sin(phi) + Bphi*cos(phi));
      Bz = sign*(Br*cos(theta) - Btheta*sin(theta));
      */
      /* Caclulate the length of the vector */
      length = sqrt(Bx*Bx+By*By+Bz*Bz);
      /* Normalize vector to unit length */
      factor = 1.0/length;
      Bx *= factor;
      By *= factor;
      Bz *= factor;
      Br *= factor;
      Btheta *= factor;
      Bphi *= factor;
      if(firsttime == 1) {
	firsttime = 0;
	improd = 0;
      }else {
	improd = Bx*Bxold + By*Byold + Bz*Bzold;
      }
      improd = 1.0-improd;
      if(improd != 0)
	improd = 1e-10/improd;
      if(improd > 10000)
	improd = 10000;
      if(improd < 10)
	improd = 10;
      Bxold = Bx;
      Byold = By;
      Bzold = Bz;

      //      improd = 100;
      /* Calculate the factor to normalize the vector */
      factor = (rmax-rmin)/((double)nrpoints-1.0);
      /* Do calculation at most at a 100m scale */
      if(factor > improd || rmin < 0 || rmax < 0)
	factor = improd; 
      totlength += factor;
      Bx *= factor;
      By *= factor;
      Bz *= factor;
      Br *= factor;
      Btheta *= factor;
      Bphi *= factor;
      
      xpos += Bx;
      ypos += By;
      zpos += Bz;
      
      /*
      if(step == 1) {
	if(nrbufferpoints < fieldLine_Deutsch_d_maxBufferSize && nrbufferpointsskip%100 == 0) {
	  xbuffer[nrbufferpoints] = xpos;
	  ybuffer[nrbufferpoints] = ypos;
	  zbuffer[nrbufferpoints] = zpos;
	  nrbufferpoints++;
	  nrbufferpointsskip = 0;
	}
	nrbufferpointsskip++;
      }
      */

      r = sqrt(xpos*xpos+ypos*ypos+zpos*zpos);
      theta = polar_angle_rad_d(zpos, sqrt(xpos*xpos+ypos*ypos));
      phi = polar_angle_rad_d(xpos, ypos);

      //      fprintf(stderr, "Br=%e Btheta=%e Bphi=%e\n", sign*Br, sign*Btheta, sign*Bphi);

      if(step == 2) {
	if(totlength > totlengthsave/(double)(nrpoints-1)) {
	  /*	  printf("%e %e %e\n", xpos, ypos, zpos); */
	  totlength -= totlengthsave/(double)(nrpoints);
	  if(i >= nrpoints) {
	    fprintf(stderr, "fieldLine_Deutsch_d: Array overflow\n");
	  }else {
	    x[i] = xpos;
	    y[i] = ypos;
	    z[i] = zpos;
	    /*
	    fprintf(stderr, "written point %ld (%e/%e r=%e)\n", i, totlength, totlengthsave, r); 
	    fprintf(stderr, "  x=%e y=%e z=%e\n", xpos, ypos, zpos);
	    fprintf(stderr, "  r=%e theta=%e phi=%e\n", r, theta*180.0/M_PI, phi*180.0/M_PI);
	    fprintf(stderr, "  Br=%e Btheta=%e Bphi=%e\n", sign*Br, sign*Btheta, sign*Bphi);
	    */
	    i++;
	    //	    exit(0);
	  }
	}
      }
      //      fprintf(stderr, "Check: r=%e rmin=%e rmax=%e\n", r, rmin, rmax);
    }
    if(i != nrpoints && step == 2) {
      fprintf(stderr, "fieldLine_Deutsch_d: Array underflow (%ld/%ld)\n", i, nrpoints);
    }

    /* All the points could be stored, so use the buffer rather than recalculate things */
    /*
    if(step == 1 && nrbufferpoints < fieldLine_Deutsch_d_maxBufferSize && nrbufferpoints > nrpoints) {
      step = 9999;
      totlengthsave = totlength;
      for(n = 0; n < nrbufferpoints; n++) {
	if(n > 0) {
	  totlength += sqrt((xbuffer[n] - xbuffer[n-1])*(xbuffer[n] - xbuffer[n-1])+(ybuffer[n] - ybuffer[n-1])*(ybuffer[n] - ybuffer[n-1])+(zbuffer[n] - zbuffer[n-1])*(zbuffer[n] - zbuffer[n-1]));
	}else {
	  totlength = 0;
	  i = 0;
	}
	//	fprintf(stderr, "%ld/%ld: %e/%e\n", n, nrbufferpoints, totlength, totlengthsave);     
	
	if(totlength > totlengthsave/(double)(nrpoints-1) || i == 0) {
	  if(i != 0)
	    totlength -= totlengthsave/(double)(nrpoints);
	  if(i >= nrpoints) {
	    fprintf(stderr, "fieldLine_Deutsch_d: Array overflow\n");
	  }else {
	    x[i] = xbuffer[n];
	    y[i] = ybuffer[n];
	    z[i] = zbuffer[n];
	    i++;
	  }
	}

      }
      if(i != nrpoints) {
	fprintf(stderr, "fieldLine_Deutsch_d: Array underflow (%ld/%ld)\n", i, nrpoints);
      }
      }*/
  }
  //  exit(0);
  fprintf(stderr, "done\n"); 
}


/*
  y[1] = r/rLC
  y[2] = theta (rad)
  y[3] = phi   (rad)
 */

static double internal_psrsalsa_Deutsch_derivs_alpha;
static double internal_psrsalsa_Deutsch_derivs_P;
static double internal_psrsalsa_Deutsch_derivs_Rns;

void internal_psrsalsa_Deutsch_derivs(double x, double *y, double *dydx)
{
  double constant, Br, Btheta, Bphi;
  constant = 299792458.0*internal_psrsalsa_Deutsch_derivs_P / (2.0*M_PI);

  fieldStrength_Deutsch_d(constant*y[1], y[2], y[3], internal_psrsalsa_Deutsch_derivs_alpha, internal_psrsalsa_Deutsch_derivs_Rns, internal_psrsalsa_Deutsch_derivs_P, &Br, &Btheta, &Bphi);
  constant = 1.0/sqrt(Br*Br+Btheta*Btheta+Bphi*Bphi);
  dydx[1] = Br*constant;
  dydx[2] = Btheta*constant;
  dydx[3] = Bphi*constant;

  dydx[2] /= y[1];
  dydx[3] /= y[1]*sin(y[2]);
  //Should make it x,y,z coordinates
  //      dSpherical2dCartesian(yrk[1], yrk[2], Br*rlc, Btheta, Bphi, &Bx, &By, &Bz);
}

/* theta0 and phi0 (in deg) are the polar (from the magnetic axis) and
   azimuthal angle of the field line at rmin. alpha (in deg) is the
   angle between the rotation axis (positive z-axis) and the magnetic
   axis. All distances should be in meters and the period P in
   seconds. If rmin is negative, it is considered to be the neutron
   star radius. If rmax is negative, it will be set to 3 times the
   light cylinder radius. eps is the precision required. If too small,
   it takes forever and you might get too small stepsize errors. If
   too large not enough points will be calculated to fill the array of
   points along the field line. openfieldline will be set to 1 if the
   fieldline penetrates the light cylinder. The number of points that
   are succesfully calculated is returned. If nounderflow_error there
   are no underflow errors generated. */
long fieldLine_Deutsch_d(double theta0, double phi0, double alpha, double *x, double *y, double *z, int *openfieldline, long nrpoints, double Rns, double P, double rmin, double rmax, double eps, int nounderflow_error) 
{
#ifdef NRAVAIL
  /* This is the number of points along the field line that is stored,
     which will resampled to nrpoints. If this number is too small,
     the field line has to be calculated twice. Doesn't really appear
     to speed up things. */
  //  #define fieldLine_Deutsch_d_maxBufferSize  100000
  //  double xbuffer[fieldLine_Deutsch_d_maxBufferSize], ybuffer[fieldLine_Deutsch_d_maxBufferSize], zbuffer[fieldLine_Deutsch_d_maxBufferSize];
  //  long nrbufferpoints, nrbufferpointsskip;

  double yrk[3], yrkold[3], dydxrk[3], yscal[3], hdid, htry, hnext, totlength;
  int step, firsttime;
  double xpos, ypos, zpos, totlengthsave, rlc, rlc_inverse, Br, Btheta, Bphi, sign, Bx, By, Bz;
  long i, totalNrPointsCalculated;

  rlc = 299792458.0*P/(2.0*M_PI);
  rlc_inverse = 1.0/rlc;
  *openfieldline = 0;

  /* By default, start at the NS surface */
  if(rmin < 0)  
    rmin = Rns;
  rmin *= rlc_inverse;
  if(rmax < 0)
    rmax = 3.0;
  else
    rmax *= rlc_inverse;
  
  internal_psrsalsa_Deutsch_derivs_alpha = alpha*M_PI/180.0;
  internal_psrsalsa_Deutsch_derivs_P = P;
  internal_psrsalsa_Deutsch_derivs_Rns = Rns;

  //  nrbufferpoints = 0;
  //  nrbufferpointsskip = 0;

  /* First step: calculate total length of field line piece.
     Second step: save the requested number of points along field line */
  //  fprintf(stderr, "Input: r=%e theta=%.2lf phi=%.2lf nrpoints=%ld\n", rmin*rlc, theta0, phi0, nrpoints); 
  for(step = 1; step <= 2; step++) { 
    yrk[0] = rmin;
    yrk[1] = M_PI*theta0/180.0;
    yrk[2] = M_PI*phi0/180.0;

    /* Do a coordinate transformation from magnetic coordinates to rotational coordinates. */
    double p, t;
    pulsar_polar_coordinates_d(yrk[2], 2.0*M_PI, alpha*M_PI/180.0, yrk[1], &p, &t);
    yrk[2] = p+M_PI;
    yrk[1] = t;

    /* So z-axis is at theta=0, y-axis is at phi=90 and x-axis at phi=0 */
    xpos = rlc*yrk[0]*sin(yrk[1])*cos(yrk[2]);
    ypos = rlc*yrk[0]*sin(yrk[1])*sin(yrk[2]);
    zpos = rlc*yrk[0]*cos(yrk[1]);
    //    fprintf(stderr, "Starting with:  r=%e theta=%e phi=%e\n", yrk[0]*rlc, yrk[1]*180.0/M_PI, yrk[2]*180.0/M_PI); 
    /*    printf("%e %e %e\n", xpos, ypos, zpos); */
    //    exit(0);

    if(step == 2) {
      x[0] = xpos;
      y[0] = ypos;
      z[0] = zpos;
    }
    //else if(step == 1) {
      //      xbuffer[0] = xpos;
    //      ybuffer[0] = ypos;
    //      zbuffer[0] = zpos;
    //      nrbufferpoints++;
    //      nrbufferpointsskip++;
    //    }
    firsttime = 1;
    if(step == 2) {
      totlengthsave = totlength;
      i = 1;
    }
    totlength = 0;
    htry = 10000*rlc_inverse;
    yscal[0] = 1;
    yscal[1] = 1;
    yscal[2] = 1;
    totalNrPointsCalculated = 0;
    //    fprintf(stderr, "Check: r=%e rmin=%e rmax=%e\n", yrk[0], rmin, rmax);
    while(yrk[0] >= 0.999*rmin && (yrk[0] <= rmax || rmax < 0)) {
      //      if(phi0 == 0 && theta0 == 270)
      //      fprintf(stderr, "calculate r=%e/%.1e theta=%.2lf phi=%.2lf alpha=%.2f\n", yrk[0], rmax, yrk[1]*180/M_PI, yrk[2]*180/M_PI, alpha);
      totalNrPointsCalculated++;
      if(totalNrPointsCalculated == 1e8) {
	fprintf(stderr, "fieldLine_Deutsch_d: Interupted calculation after %ld points\n", totalNrPointsCalculated);
	return 0;
      }
	

      //      yrk[0] = r*rlc_inverse;
      //      yrk[1] = theta;
      //      yrk[2] = phi;
      internal_psrsalsa_Deutsch_derivs(totlength, yrk-1, dydxrk-1);
      yrkold[0] = yrk[0];
      yrkold[1] = yrk[1];
      yrkold[2] = yrk[2];
      if(rkqs_d(yrk-1, dydxrk-1, 3, &totlength, htry, eps, yscal-1, &hdid, &hnext, internal_psrsalsa_Deutsch_derivs) == 0) {
	fprintf(stderr, "fieldLine_Deutsch_d: trying less precision.\n");
	if(rkqs_d(yrk-1, dydxrk-1, 3, &totlength, htry, eps*10, yscal-1, &hdid, &hnext, internal_psrsalsa_Deutsch_derivs) == 0) {
	  fprintf(stderr, "fieldLine_Deutsch_d: trying even less precision.\n");
	  if(rkqs_d(yrk-1, dydxrk-1, 3, &totlength, htry, eps*100, yscal-1, &hdid, &hnext, internal_psrsalsa_Deutsch_derivs) == 0) {
	    fprintf(stderr, "fieldLine_Deutsch_d: less precision did not work, I'm giving up.\n");
	    //	    exit(0);
	  }
	}
      }
      //      fprintf(stderr, "Stepsize %e\n", hdid);
      htry = hnext;

      Br     = yrk[0] - yrkold[0];
      Btheta = yrk[1] - yrkold[1];
      Bphi   = yrk[2] - yrkold[2];
      if(firsttime == 1) {
	if(Br >= 0) {
	  sign = 1;
	  //	  fprintf(stderr, "Br=%e, so choose sign %f\n", Br, sign);
	}else {
	  sign = -1;
	  //	  fprintf(stderr, "Br=%e, so choose sign %f\n", Br, sign);
	}
	firsttime = 0;
      }
      if(sign  < 0) {
	Br     *= sign;
	Btheta *= sign;
	Bphi   *= sign;
	yrk[0] = yrkold[0] + Br;
	yrk[1] = yrkold[1] + Btheta;
	yrk[2] = yrkold[2] + Bphi;
      }
      dSpherical2dCartesian(yrk[1], yrk[2], Br*rlc, Btheta*rlc*yrk[0], Bphi*rlc*yrk[0]*sin(yrk[1]), &Bx, &By, &Bz);
      //      Bx *= sign;
      //      By *= sign;
      //      Bz *= sign;

      /*
      Bx = sign*(Br*sin(theta)*cos(phi) + Btheta*cos(theta)*cos(phi) - Bphi*sin(phi));
      By = sign*(Br*sin(theta)*sin(phi) + Btheta*cos(theta)*sin(phi) + Bphi*cos(phi));
      Bz = sign*(Br*cos(theta) - Btheta*sin(theta));
      */
      /* Caclulate the length of the vector */
      //      length = sqrt(Bx*Bx+By*By+Bz*Bz);
      /* Normalize vector to unit length */
      /*      factor = 1.0/length;
      Bx *= factor;
      By *= factor;
      Bz *= factor;
      Br *= factor;
      Btheta *= factor;
      Bphi *= factor;
      if(firsttime == 1) {
	firsttime = 0;
	improd = 0;
      }else {
	improd = Bx*Bxold + By*Byold + Bz*Bzold;
      }
      improd = 1.0-improd;
      if(improd != 0)
	improd = 1e-10/improd;
      if(improd > 10000)
	improd = 10000;
      if(improd < 10)
	improd = 10;
      Bxold = Bx;
      Byold = By;
      Bzold = Bz;
      */

      //      improd = 100;
      /* Calculate the factor to normalize the vector */
      //      factor = (rmax-rmin)/((double)nrpoints-1.0);
      /* Do calculation at most at a 100m scale */
      /*
      if(factor > improd || rmin < 0 || rmax < 0)
	factor = improd; 
      totlength += factor;
      Bx *= factor;
      By *= factor;
      Bz *= factor;
      Br *= factor;
      Btheta *= factor;
      Bphi *= factor; */
      
      xpos += Bx;
      ypos += By;
      zpos += Bz;
      if(*openfieldline == 0) {
	if(sqrt(xpos*xpos+ypos*ypos) >= rlc)
	  *openfieldline = 1;
      }
      /*
      if(step == 1) {
	if(nrbufferpoints < fieldLine_Deutsch_d_maxBufferSize && nrbufferpointsskip%100 == 0) {
	  xbuffer[nrbufferpoints] = xpos;
	  ybuffer[nrbufferpoints] = ypos;
	  zbuffer[nrbufferpoints] = zpos;
	  nrbufferpoints++;
	  nrbufferpointsskip = 0;
	}
	nrbufferpointsskip++;
      }
      */

      //      r = sqrt(xpos*xpos+ypos*ypos+zpos*zpos);
      // theta = polar_angle_rad_d(zpos, sqrt(xpos*xpos+ypos*ypos));
      //phi = polar_angle_rad_d(xpos, ypos);

      //      fprintf(stderr, "Br=%e Btheta=%e Bphi=%e\n", sign*Br, sign*Btheta, sign*Bphi);

      if(step == 2) {
	if(totlength > totlengthsave/(double)(nrpoints-1)) {
	  /*	  printf("%e %e %e\n", xpos, ypos, zpos); */
	  totlength -= totlengthsave/(double)(nrpoints);
	  if(i >= nrpoints) {
	    fprintf(stderr, "fieldLine_Deutsch_d: Array overflow\n");
	  }else {
	    x[i] = xpos;
	    y[i] = ypos;
	    z[i] = zpos;
	    /*
	    fprintf(stderr, "written point %ld (%e/%e r=%e)\n", i, totlength, totlengthsave, r); 
	    fprintf(stderr, "  x=%e y=%e z=%e\n", xpos, ypos, zpos);
	    fprintf(stderr, "  r=%e theta=%e phi=%e\n", r, theta*180.0/M_PI, phi*180.0/M_PI);
	    fprintf(stderr, "  Br=%e Btheta=%e Bphi=%e\n", sign*Br, sign*Btheta, sign*Bphi);
	    */
	    i++;
	    //	    exit(0);
	  }
	}
      }
      //      fprintf(stderr, "Check: r=%e rmin=%e rmax=%e\n", r, rmin, rmax);
    }
    if(i != nrpoints && step == 2 && nounderflow_error == 0) {
      fprintf(stderr, "fieldLine_Deutsch_d: Array underflow (%ld/%ld)\n", i, nrpoints);
    }

    /* All the points could be stored, so use the buffer rather than recalculate things */
    /*
    if(step == 1 && nrbufferpoints < fieldLine_Deutsch_d_maxBufferSize && nrbufferpoints > nrpoints) {
      step = 9999;
      totlengthsave = totlength;
      for(n = 0; n < nrbufferpoints; n++) {
	if(n > 0) {
	  totlength += sqrt((xbuffer[n] - xbuffer[n-1])*(xbuffer[n] - xbuffer[n-1])+(ybuffer[n] - ybuffer[n-1])*(ybuffer[n] - ybuffer[n-1])+(zbuffer[n] - zbuffer[n-1])*(zbuffer[n] - zbuffer[n-1]));
	}else {
	  totlength = 0;
	  i = 0;
	}
	//	fprintf(stderr, "%ld/%ld: %e/%e\n", n, nrbufferpoints, totlength, totlengthsave);     
	
	if(totlength > totlengthsave/(double)(nrpoints-1) || i == 0) {
	  if(i != 0)
	    totlength -= totlengthsave/(double)(nrpoints);
	  if(i >= nrpoints) {
	    fprintf(stderr, "fieldLine_Deutsch_d: Array overflow\n");
	  }else {
	    x[i] = xbuffer[n];
	    y[i] = ybuffer[n];
	    z[i] = zbuffer[n];
	    i++;
	  }
	}

      }
      if(i != nrpoints) {
	fprintf(stderr, "fieldLine_Deutsch_d: Array underflow (%ld/%ld)\n", i, nrpoints);
      }
      }*/
  }
  //  exit(0);
  //  fprintf(stderr, "done\n"); 
  return i;
#else
  fprintf(stderr, "CODE NEEDS TO BE RECOMPILED WITH NR SUPPORT BEFORE FUNCTION fieldLine_Deutsch_d WILL WORK.\n");
  exit(0);
#endif
}


/* Z is up, X is right, Y is inside. Rotates X towards Y. Angle is in degrees.  */
void rotateZ(long nrpoints, float angle, float *x, float *y, float *z)
{
  long n;
  float x2, y2;
  angle *= M_PI/180.0;
  for(n = 0; n < nrpoints; n++) {
    x2 = x[n]*cos(angle);
    y2 = x[n]*sin(angle);
    x2 -= y[n]*sin(angle);
    y2 += y[n]*cos(angle);
    x[n] = x2;
    y[n] = y2;
  }
}

void rotateZ_d(long nrpoints, double angle, double *x, double *y, double *z)
{
  long n;
  double x2, y2;
  angle *= M_PI/180.0;
  for(n = 0; n < nrpoints; n++) {
    x2 = x[n]*cos(angle);
    y2 = x[n]*sin(angle);
    x2 -= y[n]*sin(angle);
    y2 += y[n]*cos(angle);
    x[n] = x2;
    y[n] = y2;
  }
}

/* Z is up, X is right, Y is inside. Rotates Z towards X. Angle is in degrees.  */
void rotateY(long nrpoints, float angle, float *x, float *y, float *z)
{
  long n;
  float x2, z2;
  angle *= M_PI/180.0;
  for(n = 0; n < nrpoints; n++) {
    x2 = x[n]*cos(angle);
    z2 = -x[n]*sin(angle);
    x2 += z[n]*sin(angle);
    z2 += z[n]*cos(angle);
    x[n] = x2;
    z[n] = z2;
  }
}

void rotateY_d(long nrpoints, double angle, double *x, double *y, double *z)
{
  long n;
  double x2, z2;
  angle *= M_PI/180.0;
  for(n = 0; n < nrpoints; n++) {
    x2 = x[n]*cos(angle);
    z2 = -x[n]*sin(angle);
    x2 += z[n]*sin(angle);
    z2 += z[n]*cos(angle);
    x[n] = x2;
    z[n] = z2;
  }
}

/* Z is up, X is right, Y is inside. Rotates Z towards Y. Angle is in degrees.  */
void rotateX(long nrpoints, float angle, float *x, float *y, float *z)
{
  long n;
  float y2, z2;
  angle *= M_PI/180.0;
  for(n = 0; n < nrpoints; n++) {
    y2 = y[n]*cos(angle);
    z2 = -y[n]*sin(angle);
    y2 += z[n]*sin(angle);
    z2 += z[n]*cos(angle);
    y[n] = y2;
    z[n] = z2;
  }
}

void rotateX_d(long nrpoints, double angle, double *x, double *y, double *z)
{
  long n;
  double y2, z2;
  angle *= M_PI/180.0;
  for(n = 0; n < nrpoints; n++) {
    y2 = y[n]*cos(angle);
    z2 = -y[n]*sin(angle);
    y2 += z[n]*sin(angle);
    z2 += z[n]*cos(angle);
    y[n] = y2;
    z[n] = z2;
  }
}

/* Lorentz transform (normalized) propagation direction (dx,dy,dz) on
   position (x,y,z) in the corotating frame to the observer frame. rlc
   is the light cylinder distance. All distances should be in the same
   units. */
void aberration(float x, float y, float z, float dx, float dy, float dz, float *dx2, float *dy2, float *dz2, float rlc)
{
  float phi, gamma_lorentz, c1, c2;
  float beta_lorentz, betax, betay, betaz, improd;
  /* Angle in x-y plane (z is rotation axis) */
  phi = atan2(y,x);
  /* The corotation velocity (km/s) */
  beta_lorentz = sqrt(x*x+y*y)/rlc;
  /* Beam rotates from x to y */
  betax = -beta_lorentz*sin(phi);
  betay = beta_lorentz*cos(phi);
  betaz = 0;
  gamma_lorentz = 1.0/sqrt(1.0-beta_lorentz*beta_lorentz);

  improd = dx*betax + dy*betay + dz*betaz;
  c1 = gamma_lorentz+(gamma_lorentz-1.0)*improd/(beta_lorentz*beta_lorentz);
  c2 = 1.0/(gamma_lorentz*(1.0+improd));
  *dx2 = c2*(dx + betax*c1);
  *dy2 = c2*(dy + betay*c1);
  *dz2 = c2*(dz + betaz*c1);
}

void aberration_d(double x, double y, double z, double dx, double dy, double dz, double *dx2, double *dy2, double *dz2, double rlc)
{
  double phi, gamma_lorentz, c1, c2;
  double beta_lorentz, betax, betay, betaz, improd;
  /* Angle in x-y plane (z is rotation axis) */
  phi = atan2(y,x);
  /* The corotation velocity (km/s) */
  beta_lorentz = sqrt(x*x+y*y)/rlc;
  /* Beam rotates from x to y */
  betax = -beta_lorentz*sin(phi);
  betay = beta_lorentz*cos(phi);
  betaz = 0;
  gamma_lorentz = 1.0/sqrt(1.0-beta_lorentz*beta_lorentz);

  improd = dx*betax + dy*betay + dz*betaz;
  c1 = gamma_lorentz+(gamma_lorentz-1.0)*improd/(beta_lorentz*beta_lorentz);
  c2 = 1.0/(gamma_lorentz*(1.0+improd));
  *dx2 = c2*(dx + betax*c1);
  *dy2 = c2*(dy + betay*c1);
  *dz2 = c2*(dz + betaz*c1);
}

/* Formula's from: Dyks & Rudak 2003, ApJ 598 1201: Two pole caustic model 

Takes as input: position (x,y,z) and normalized direction (dx, dy, dz)
and the light cylinder radius rlc. Distances should be in the same
units. Returns pulse longitude and angle between rotation axis and los
(in degrees). dt_ret is the time travel delay in degrees.
*/
void timedelays(float x, float y, float z, float dx, float dy, float dz, float rlc, float *phi, float *zeta, float *dt_ret, int noaberration, int noretardation)
{
  float dx2, dy2, dz2, dt;

  /* Get aberrated propagation direction */
  if(noaberration) {
    dx2 = dx;
    dy2 = dy;
    dz2 = dz;
  }else {
    aberration(x, y, z, dx, dy, dz, &dx2, &dy2, &dz2, rlc); 
  }

  /* This is the time delay part (in radians) */
  if(noretardation == 0) {
    *dt_ret = -x*dx2;
    *dt_ret -= y*dy2;
    *dt_ret -= z*dz2;
    *dt_ret /= rlc;
  }else {
    *dt_ret = 0;
  }

  /* This is the retardated/aberrated propagation direction angle (in radians) */
  dt = *dt_ret - atan2(dy2,dx2);
  *phi = 180.0*dt/M_PI;
  *dt_ret *= 180.0/M_PI;

  /* This is the angle between the propagation direction and the rotation axis */
  *zeta = sqrt(dx2*dx2 + dy2*dy2);
  *zeta = 180*atan2(*zeta, dz2)/M_PI;
}

void timedelays_d(double x, double y, double z, double dx, double dy, double dz, double rlc, double *phi, double *zeta, double *dt_ret, int noaberration, int noretardation)
{
  double dx2, dy2, dz2, dt;
  if(noaberration) {
    dx2 = dx;
    dy2 = dy;
    dz2 = dz;
  }else {
    aberration_d(x, y, z, dx, dy, dz, &dx2, &dy2, &dz2, rlc); 
  }
  if(noretardation == 0) {
    *dt_ret = -x*dx2;
    *dt_ret -= y*dy2;
    *dt_ret -= z*dz2;
    *dt_ret /= rlc;
  }else {
    *dt_ret = 0;
  }
  dt = *dt_ret - atan2(dy2,dx2);
  *phi = 180.0*dt/M_PI;
  *dt_ret *= 180.0/M_PI;
  *zeta = sqrt(dx2*dx2 + dy2*dy2);
  *zeta = 180*atan2(*zeta, dz2)/M_PI;
}

void crossprod(float a1, float a2, float a3, float b1, float b2, float b3, float *c1, float *c2, float *c3)
{
  *c1 = a2*b3-a3*b2;
  *c2 = a3*b1-a1*b3;
  *c3 = a1*b2-a2*b1;
}

// answer = axb
void crossprod_array(float *a, float *b, float *answer)
{
  answer[0] = a[1]*b[2]-a[2]*b[1];
  answer[1] = a[2]*b[0]-a[0]*b[2];
  answer[2] = a[0]*b[1]-a[1]*b[0];
}

double dotproduct_array(float *a, float *b)
{
  double answer;
  answer = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  return answer;
}

double dotproduct_array_d(double *a, double *b)
{
  double answer;
  answer = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  return answer;
}

double length_vector_array(float *a)
{
  double length;
  length = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  return length;
}

double length_vector_array_d(double *a)
{
  double length;
  length = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  return length;
}

void crossprod_array_d(double *a, double *b, double *answer)
{
  answer[0] = a[1]*b[2]-a[2]*b[1];
  answer[1] = a[2]*b[0]-a[0]*b[2];
  answer[2] = a[0]*b[1]-a[1]*b[0];
}

void crossprod_d(double a1, double a2, double a3, double b1, double b2, double b3, double *c1, double *c2, double *c3)
{
  *c1 = a2*b3-a3*b2;
  *c2 = a3*b1-a1*b3;
  *c3 = a1*b2-a2*b1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* All angles are in degrees. PA is an angle between 0 and 180. Note
   there is also a double version below. */
float paswing(float alpha, float beta, float l, float pa0, float l0, int nrJumps, float *jump_longitude, float *jump_offset, float add_height_longitude, float add_height_shift)
{
  float x1, y1, sa, dl, pa, dbcw, dha;
  int i;
  alpha *= M_PI/180.0;
  beta *= M_PI/180.0;
  if(l >= add_height_longitude) {
    dbcw = 2.0*add_height_shift*180.0/M_PI;
    // see Hibschman and  Arons, 2001, ApJ, 546, 382
    dha = 10.0*add_height_shift*cos(alpha)*180.0/(3.0*M_PI); 
    /*    printf("%f %f %f\n", add_height_longitude, add_height_shift, dbcw); */
  }else {
    dbcw = 0;
    dha = 0;
  }
  sa = sin(alpha);
  dl = (l-(l0+dbcw))*M_PI/180.0;
  y1 = sa*sin(dl);
  x1 = sin(alpha+beta)*cos(alpha);
  x1 -= cos(alpha+beta)*sa*cos(dl);
  pa = pa0+dha+atan2(y1,x1)*180.0/M_PI;
  if(nrJumps > 0) {
    for(i = 0; i < nrJumps; i++)
      if(l > jump_longitude[i])
	pa += jump_offset[i];
  }
  pa = derotate_180(pa);
  return pa;
}

/* All angles are in degrees. PA is an angle between 0 and 180. Note there is also a float version above */
double paswing_double(double alpha, double beta, double l, double pa0, double l0, int nrJumps, double *jump_longitude, double *jump_offset, double add_height_longitude, double add_height_shift)
{
  double x1, y1, sa, dl, pa, dbcw, dha;
  int i;
  alpha *= M_PI/180.0;
  beta *= M_PI/180.0;
  if(l >= add_height_longitude) {
    dbcw = 2.0*add_height_shift*180.0/M_PI; 
    dha = 10.0*add_height_shift*cos(alpha)*180.0/(3.0*M_PI); 
    /*    printf("%f %f %f\n", add_height_longitude, add_height_shift, dbcw); */
  }else {
    dbcw = 0;
    dha = 0;
  }
  sa = sin(alpha);
  dl = (l-(l0+dbcw))*M_PI/180.0;
  y1 = sa*sin(dl);
  x1 = sin(alpha+beta)*cos(alpha);
  x1 -= cos(alpha+beta)*sa*cos(dl);
  pa = pa0+dha+atan2(y1,x1)*180.0/M_PI;
  if(nrJumps > 0) {
    for(i = 0; i < nrJumps; i++)
      if(l > jump_longitude[i])
	pa += jump_offset[i];
  }
  pa = derotate_180_double(pa);
  return pa;
}

//START REGION DEVELOP

/* Prints an angle in radians as a time */
void printhour(double t)
{
  int th, tm, ts;
  t *= 12/M_PI;
  if(t < 0) t += 24;
  if(t > 24) t -= 24;
  th = t;
  tm = 60*t - 60*th;
  ts = 60*60*t - 60*60*th - 60*tm;
  printf("%02d:%02d:%02d", th, tm, ts);
}

/* Converts the (vector r, theta, phi) into the Cartesian vector (x,
   y, z). The z-axis is the polar axis (theta=0), the x-axis to phi=0
   and the y-axis to phi=pi/2. All angles are in radians. */
void spherical2Cartesian(double r, double theta, double phi, double *x, double *y, double *z)
{
  *x = r*sin(theta)*cos(phi);
  *y = r*sin(theta)*sin(phi);
  *z = r*cos(theta);
}

/* Reverse of spherical2Cartesian. */
void cartesian2spherical(double *r, double *theta, double *phi, double x, double y, double z)
{
  float r2;
  *r = sqrt(x*x+y*y+z*z);
  r2 = sqrt(x*x+y*y);
  *theta = atan2(r2, z);
  *phi = atan2(y, x);
}

/*
  Given the longitude and latitude (in radians), calculate the
  Hammer-Aitoff projection Cartesian coordinates x and y. The range of
  x is between -2 and 2, and that of y between -1 and 1. The angles
  dlongitude and dlatitude correspond to a rotation of the sphere such
  that positive angles correspond with the equator moving towards the
  Z-pole and a rotation in the positive longitude direction. 
*/
void projectionHammerAitoff_xy(float longitude, float latitude, float dlongitude, float dlatitude, float *x, float *y)
{
  float f;
  double x2, y2, z2, theta, phi, r;

  dlongitude *= 180.0/M_PI;
  dlatitude  *= 180.0/M_PI;

  spherical2Cartesian(1, latitude+0.5*M_PI, longitude, &x2, &y2, &z2);
  //  if(ok) printf("XXXX %f %f -> %lf %lf %lf\n", longitude*180.0/M_PI, latitude*180.0/M_PI, x2, y2, z2);
  rotateZ_d(1, dlongitude, &x2, &y2, &z2);
  //  if(ok) printf("XXXX       -> %lf %lf %lf   (%f)\n", x2, y2, z2, dlongitude*180.0/M_PI);
  rotateY_d(1, -dlatitude, &x2, &y2, &z2);
  //  if(ok) printf("XXXX       -> %lf %lf %lf   (%f)\n", x2, y2, z2, dlatitude*180.0/M_PI);
  cartesian2spherical(&r, &theta, &phi, x2, y2, z2);
  latitude = theta - 0.5*M_PI;
  longitude = phi;
  //  if(ok) printf("XXXX       -> %lf %lf -> %f %f\n", theta*180.0/M_PI, phi*180.0/M_PI, longitude*180.0/M_PI, latitude*180.0/M_PI);

  f = 1.0/sqrt(1.0+cos(latitude)*cos(0.5*longitude));
  *x = f*2.0*sqrt(2.0)*cos(latitude)*sin(0.5*longitude);
  *y = sqrt(2.0)*sin(latitude)*f;
  /* rescale x to be between -2 and 2 and the vertical coordinate
     between -1 and 1. */
  *x *= 0.5*sqrt(2);
  *y *= 0.5*sqrt(2);
}

/* Given the Cartesian coordinates x and y of the Hammer-Aitoff
   projection, calculate the corresponding longitude and latitude (in
   radians). The range of x is between -2 and 2, and that of y between
   -1 and 1. The angles dlongitude and dlatitude correspond to a
   rotation of the sphere such that positive angles correspond with
   the equator moving towards the Z-pole and a rotation in the
   positive longitude direction.

   Return value:
     1 = ok
     0 = (x,y) outside projection
*/
int projectionHammerAitoff_longlat(float x, float y, float dlongitude, float dlatitude, float *longitude, float *latitude)
{
  double x2, y2, z2, theta, phi, r;
  float z;
  x *= 1.0*sqrt(2);
  y *= 1.0*sqrt(2);
  z = 1.0-(0.25*x)*(0.25*x)-(0.5*y)*(0.5*y);
  if(z < 0.5)
    return 0;
  z = sqrt(z);
  *longitude = 2.0*atan(0.5*z*x/(2.0*z*z-1.0));
  *latitude = asin(z*y);

  dlongitude *= 180.0/M_PI;
  dlatitude  *= 180.0/M_PI;

  spherical2Cartesian(1, *latitude+0.5*M_PI, *longitude, &x2, &y2, &z2);
  //  if(ok) printf("XXXX %f %f -> %lf %lf %lf\n", longitude*180.0/M_PI, latitude*180.0/M_PI, x2, y2, z2);
  //  if(ok) printf("XXXX       -> %lf %lf %lf   (%f)\n", x2, y2, z2, dlongitude*180.0/M_PI);
  rotateY_d(1, dlatitude, &x2, &y2, &z2);
  rotateZ_d(1, -dlongitude, &x2, &y2, &z2);
  //  if(ok) printf("XXXX       -> %lf %lf %lf   (%f)\n", x2, y2, z2, dlatitude*180.0/M_PI);
  cartesian2spherical(&r, &theta, &phi, x2, y2, z2);
  *latitude = theta - 0.5*M_PI;
  *longitude = phi;

  return 1;
}


/*
  Given the longitude and latitude (in radians), calculate the
  projection Cartesian coordinates x and y, which looks like a
  sphere. The sphere has a radius of 1. If the point lies on the
  front-side, the sphere is centred at x=-1.125, y=0. If the point
  lies on the far side the point ends on a sphere with an x-offset in
  the opposite direction. The angles dlongitude and dlatitude
  correspond to a rotation of the sphere such that positive angles
  correspond with the equator moving towards the Z-pole and a rotation
  in the positive longitude direction. The weight is set to the cosine
  of the angle between the normal and the line of sight.

  Return value:
     0 = point in left-hand sphere (front of sphere)
     1 = point in right-hand sphere (back of sphere)
*/
int projection_sphere_xy(float longitude, float latitude, float dlongitude, float dlatitude, float *x, float *y, float *weight)
{
  float xE, yE, zE;
  int side;

  dlongitude *= 180.0/M_PI;
  dlatitude  *= 180.0/M_PI;

  /* Calculate point on sphere with radius 1 centred at origin. Take
     x=right, y=up, z=away from observer */
  yE = sin(latitude);
  xE = cos(latitude)*sin(longitude);
  zE = -cos(latitude)*cos(longitude);

  rotateY(1, -dlongitude, &xE, &yE, &zE);
  rotateX(1, dlatitude, &xE, &yE, &zE);

  /* Dot product with l.o.s. gives cos(angle), where angle is the
     angle with the line of sight. */
  *weight = fabs(zE);

  /* If positive z -> point is obscured by other points in first half sphere. */
  if(zE > 0)
    side = 1;
  else
    side = 0;

  /* A line connecting the observer (in origin) to the point on the
  sphere located at distance z=distE as function of z is:

  x = z*xE/(distE+zE);
  y = z*yE/(distE+zE);

  The projection is where those lines intersect a screen at distance z=distS

  x = distS*xE/(distE+zE);
  y = distS*yE/(distE+zE);

  Take distE to be effectively infinity (>> distS)

  */

  *x = xE;
  *y = yE;

  if(side == 0)
    *x -= 1.125;
  else
    *x = (-*x)+1.125;
  
  return side;
}

/* Given the Cartesian coordinates x and y of a spherical projection,
   calculate the corresponding longitude and latitude (in
   radians). The sphere has a radius of 1. If the point lies on the
   front-side, the sphere is centred at x=-1.125, y=0. If the point
   lies on the far side the point ends on a sphere with an x-offset in
   the opposite direction. The angles dlongitude and dlatitude
   correspond to a rotation of the sphere such that positive angles
   correspond with the equator moving towards the Z-pole and a
   rotation in the positive longitude direction.

   Return value:
     1 = ok
     0 = (x,y) outside projection
*/
int projection_sphere_longlat(float x, float y, float dlongitude, float dlatitude, float *longitude, float *latitude)
{
  float z, r2;
  int side;

  if(x < 0) {
    x += 1.125;
    side = 0;
  }else {
    x = -(x - 1.125);
    side = 1;
  }

  if(x*x+y*y > 1)
    return 0;

  if(side == 0)
    z = -sqrt(1-x*x-y*y);
  else
    z = sqrt(1-x*x-y*y);

  rotateX(1, -dlatitude*180.0/M_PI, &x, &y, &z);
  rotateY(1, dlongitude*180.0/M_PI, &x, &y, &z);

  r2 = sqrt(z*z+x*x);

  *latitude = atan2(r2, y);
  *latitude = 0.5*M_PI-(*latitude);

  *longitude = atan2(x, -z);

  return 1;
}

/* Given the Cartesian coordinates x and y of a simple linear
   longitude-latitude map projection, calculate the corresponding
   longitude and latitude (in radians). The x and y values are
   expected to be in units of 80 degrees degrees (to make map size
   compatible with make_projection_map()). The angles dlongitude and
   dlatitude correspond to a rotation of the sphere such that positive
   angles correspond with the equator moving towards the Z-pole and a
   rotation in the positive longitude direction.

   Return value:
     1 = ok
     0 = (x,y) outside projection
*/
int projection_longlat_longlat(float x, float y, float dlongitude, float dlatitude, float *longitude, float *latitude)
{
  float xE, yE, zE, r2;

  x *= 80.0;
  y *= 80.0;
  
  if(x < -180 || x > 180 || y < -90 || y > 90)
    return 0;

  *latitude = y*M_PI/180.0;
  *longitude = x*M_PI/180.0;
  
  /* Calculate point on sphere with radius 1 centred at origin. Take
     x=right, y=up, z=away from observer */
  yE = sin(*latitude);
  xE = cos(*latitude)*sin(*longitude);
  zE = -cos(*latitude)*cos(*longitude);

  rotateY(1, dlongitude*180.0/M_PI, &xE, &yE, &zE);
  rotateX(1, -dlatitude*180.0/M_PI, &xE, &yE, &zE);

  r2 = sqrt(zE*zE+xE*xE);

  *latitude = atan2(r2, yE);
  *latitude = 0.5*M_PI-(*latitude);
  *longitude = atan2(xE, -zE);

  return 1;
}

/* Given the longitude and latitude (in radians), calculate the
   corresponding Cartesian coordinates x and y of a simple linear
   longitude-latitude map projection (these will be in degrees to be
   compatable with drawSphericalGrid()). This projection is very
   simple. The angles dlongitude and dlatitude correspond to a
   rotation of the sphere such that positive angles correspond with
   the equator moving towards the Z-pole and a rotation in the
   positive longitude direction.

*/
void projection_longlat_xy(float longitude, float latitude, float dlongitude, float dlatitude, float *x, float *y)
{
  float xE, yE, zE, r2;

  dlongitude *= 180.0/M_PI;
  dlatitude  *= 180.0/M_PI;

  /* Calculate point on sphere with radius 1 centred at origin. Take
     x=right, y=up, z=away from observer */
  yE = sin(latitude);
  xE = cos(latitude)*sin(longitude);
  zE = -cos(latitude)*cos(longitude);

  rotateX(1, dlatitude, &xE, &yE, &zE);
  rotateY(1, -dlongitude, &xE, &yE, &zE);

  r2 = sqrt(zE*zE+xE*xE);

  *y = atan2(r2, yE);
  *y = 0.5*M_PI-(*y);
  *y *= 180.0/M_PI;
  *x = atan2(xE, -zE)*180.0/M_PI;
}

/* 
   This function fills the map (of nrx by nry pixels, idealy with a
   ration 2:1) using the a spherical projection. The area outside the
   spherical surface will be filled with value background. For
   projection 1 and 2 the horizontal axis runs from -2.25 to 2.25 and
   the vertical axis from -1.125 to 1.125. For projection 1 the poles
   are at (0, +1) and (0, -1) and the equator extends from (-2, 0) to
   (+2, 0). For projection 2 the sphere has a radius of 1. If a point
   lies on the front-side, the sphere is centred at x=-1.125, y=0. If
   the point lies on the far side the point ends on a sphere with an
   x-offset in the opposite direction. For projection 3 the horizontal
   axis runs from -180 to 180 degrees and the vertical axis from -90
   to 90 degrees. The angles dlongitude and dlatitude correspond to a
   rotation of the sphere such that positive angles correspond with
   the equator moving towards the Z-pole and a rotation in the
   positive longitude direction.

   projection = 1 Hammer-Aitoff projection. 
   projection = 2 3D sphere
   projection = 3 Simple linear longitude-latitude map

   funk is a function that accepts two floats (longitude and latitude
   in radians) and returns the value of the map.
*/
void make_projection_map(float *map, int nrx, int nry, float background, float dlongitude, float dlatitude, float (*funk)(float, float), int projection)
{
  int i, j;
  float x, y, l, b;
  for(i = 0; i < nrx; i++) {
    for(j = 0; j < nry; j++) {
      x = 4.5*(i-0.5*nrx)/(float)nrx;
      y = 2.25*(j-0.5*nry)/(float)nry;
      if(projection == 1) {
	if(projectionHammerAitoff_longlat(x, y, dlongitude, dlatitude,  &l, &b))
	  map[i+nrx*j] = funk(l, b);
	else
	  map[i+nrx*j] = background;
      }else if(projection == 2) {
	if(projection_sphere_longlat(x, y, dlongitude, dlatitude,  &l, &b))
	  map[i+nrx*j] = funk(l, b);
	else
	  map[i+nrx*j] = background;
      }else if(projection == 3) {
	if(projection_longlat_longlat(x, y, dlongitude, dlatitude,  &l, &b))
	  map[i+nrx*j] = funk(l, b);
	else
	  map[i+nrx*j] = background;
      }
    }
  }
}
