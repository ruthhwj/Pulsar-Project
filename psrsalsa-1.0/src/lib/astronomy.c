//START REGION RELEASE
#include <stdio.h>
#include <math.h>
#include "psrsalsa.h"
#include "sla_wrap.h"

//START REGION DEVELOP
#include <string.h>

//START REGION RELEASE

/* Returns the DM delay in seconds at freq freq relative to freq_ref
   (so returns a positive value is freq < freq_ref). If inffrq is
   nonzero, the reference frequency is set to infinity. The
   frequencies are in MHz and the DM in pc/cm^3. */
long double calcDMDelay(long double freq, long double freq_ref, int inffrq, long double dm)
{
  long double dt;
  if(freq <= 0.01)
    return 0;

  if(inffrq)
    dt = DM_CONST*dm/(freq*freq);
  else
    dt = DM_CONST*dm*(1.0/(freq*freq)-1.0/(freq_ref*freq_ref));
  return dt;
}

/* Returns the amount of Faraday rotation (of the position angle in
   radians) according to the RM at frequency freq relative to
   freq_ref. If inffrq is nonzero, the reference frequency is set to
   infinity. The frequency is in MHz and the RM in radians/meter^2. */
float calcRMAngle(float freq, float freq_ref, int inffrq, float rm)
{
  float lambda, lambda2, phi;
  lambda = 299792458.0*(1.0e-6)/freq;
  if(freq <= 0.01)
    return 0;
  if(inffrq) {
    phi = rm*lambda*lambda;
  }else {
    lambda2 = 299792458.0*(1.0e-6)/freq_ref;
    phi = rm*(lambda*lambda-lambda2*lambda2);
  }
  return phi;
}

//START REGION DEVELOP


/***********************
TOA = 35 (name=t090410_201838.SFTC flags=)
FREQ = 1362.176000
MJD = 54931.851562
wiki:  UT  = 2009 04 09 20:24:25
Meeus: UT  = 2009 04 10 08:24:25

Note: SAT = 54931.8471065206587 = 6.42 min earlier than BAT (i.e calculated UT exactly 0.5 day less)

 ***********************/

/*
Astronomical Formulae for Calculators, Jean Meeus: H3 Julian Day and Calendar Date, p26
*/

void mjd2date_old(long double mjd, int *year, int *month, int *day, int *hour, int *minute, float *seconds)
{
  long A, B, C, D, E, Z, a1, a2;
  double F;
  /*
    1957 4.81 oct = 2436116.31
    333 27.5 jan = 1842713.00
  */
  double JD = mjd + 2400000.5;
  /* Actually, this doesn't make sense. MJD XXXX.0 should be midnight,
     XXXX.5 afternoon, but maybe only in astronomy? Different
     definitions of MJD around? */
  JD += 0.5;
  int sign = 1;
  if(JD < 0)
    sign = -1;
  Z = sign*floor(sign*JD);
  F = JD - Z;
  A = Z;
  if(Z >= 2299161)
    {
      a1 = (Z - 1867216.25)/36524.25;
      a2 = a1/4.0;
      A += 1 + a1 - a2;
    }
  B = A + 1524;
  C = (B - 122.13)/365.25;
  D = 365.25*C;
  E = (B-D)/30.6001;
  a1 = 30.6001*E;
  *day = B - D - a1;
  *hour = F*24;
  *minute = (F-(*hour)/24.0)*24*60.0;
  *seconds = (F-(*hour+(*minute)/60.0)/24.0)*24*60.0*60.0;
  if(E < 14)
    *month=E - 1;
  else
    *month=E - 13;
  if(*month > 2)
    *year = C - 4716;
  else
    *year = C - 4715;
  if(*year < 1)
    *year -= 1;
}


//START REGION DEVELOP
//START REGION RELEASE

void internal_mjd_A_to_year_month_day(long Z, int *year, int *month, int *day)
{
  long A, B, C, D, E, a1, a2;
  A = Z;
  if(Z >= 2299161)
    {
      a1 = (Z - 1867216.25)/36524.25;
      a2 = a1/4.0;
      A += 1 + a1 - a2;
    }
  B = A + 1524;
  C = (B - 122.13)/365.25;
  D = 365.25*C;
  E = (B-D)/30.6001;
  a1 = 30.6001*E;
  *day = B - D - a1;
  if(E < 14)
    *month=E - 1;
  else
    *month=E - 13;
  if(*month > 2)
    *year = C - 4716;
  else
    *year = C - 4715;
  if(*year < 1)
    *year -= 1;
}

//START REGION DEVELOP

void mjd2date(long double mjd, int *year, int *month, int *day, int *hour, int *minute, float *seconds)
{
  long Z;
  long double F;
  /*
    1957 4.81 oct = 2436116.31
    333 27.5 jan = 1842713.00
  */
  long double JD = mjd + 2400000.5;
  /* Actually, this doesn't make sense. MJD XXXX.0 should be midnight,
     XXXX.5 afternoon, but maybe only in astronomy? Different
     definitions of MJD around? */
  JD += 0.5;
  int sign = 1;
  if(JD < 0)
    sign = -1;
  Z = sign*floor(sign*JD);
  F = JD - Z;

  internal_mjd_A_to_year_month_day(Z, year, month, day);

  *hour = F*24;
  *minute = (F-(*hour)/24.0)*24*60.0;
  *seconds = (F-(*hour+(*minute)/60.0)/24.0)*24*60.0*60.0;
}

//START REGION DEVELOP
//START REGION RELEASE

/* 
  Type 1: YYYY-MM-DDXHH:MM:SS

  X = separator
  precision gives the number of decimal places in the seconds.
*/
void mjd2dateString(long double mjd, char *string, int precision, int type, char *separator)
{
  long Z;
  long double F;
  char timestr[100];
  int year, month, day;
  /*
    1957 4.81 oct = 2436116.31
    333 27.5 jan = 1842713.00
  */
  long double JD = mjd + 2400000.5;
  /* Actually, this doesn't make sense. MJD XXXX.0 should be midnight,
     XXXX.5 afternoon, but maybe only in astronomy? Different
     definitions of MJD around? */
  JD += 0.5;
  int sign = 1;
  if(JD < 0)
    sign = -1;
  Z = sign*floor(sign*JD);
  F = JD - Z;

  converthms_string(timestr, F*24.0, precision, 1);
  if(timestr[0] == '2' && timestr[1] == '4') {
    /* Time after rounding is 24:00:00, set to 0 hours and add a day. */
    timestr[0] = '0';
    timestr[1] = '0';
    internal_mjd_A_to_year_month_day(Z+1, &year, &month, &day);
  } else {
    internal_mjd_A_to_year_month_day(Z, &year, &month, &day);
  }


  if(type == 1) {
    sprintf(string, "%04d-%02d-%02d%s%s", year, month, day, separator, timestr);
  }else {
    fflush(stdout);
    fprintf(stderr, "ERROR: Unknown type selected in mjd2dateString!");
    return;
  }
}

//START REGION DEVELOP
//START REGION RELEASE

/* Converts floating point number into a string with hours (or
   degrees), minutes and seconds. The number should be in hours (for
   instance 12.0 -> 12:00:00) or in degrees (-30.5 ->
   -30:30:00). precision states the number of decimals in the
   seconds. If precision is set to a negative number, the output is in
   minute precision only. Type controls the formatting:
   type 1 = 12:34:56.000
   type 2 = 12h34m56.000s
   type 3 = 12d34'56.000"
   type 4 = hhmmss
*/
void converthms_string(char *hms, long double number, int precision, int type)
{
  char dummy_str[100], dummy_str2[100];
  int hour, minute;
  if(number < 0) {
    sprintf(hms, "-");
    number *= -1;
  }else {
    hms[0] = 0;
  }
  hour = number;
  number -= hour;
  number *= 60.0;
  minute = number;
  number -= minute;
  number *= 60.0;

  if(precision < 0) {
    if(number >= 30) {
      minute += 1;
      if(minute == 60) {
	minute = 0;
	hour += 1;
      }
    }
    if(type == 2)
      sprintf(hms, "%02dh%02dm", hour, minute);
    else if(type == 3)
      sprintf(hms, "%02dd%02d'", hour, minute);
    else if(type == 4)
      sprintf(hms, "%02d%02d'", hour, minute);
    else
      sprintf(hms, "%02d:%02d", hour, minute);
  }else {
    if(precision == 0) {
      sprintf(dummy_str2, "%%0%d.%dLf", 2, precision);
    }else {
      sprintf(dummy_str2, "%%0%d.%dLf", precision+3, precision);
    }
    sprintf(dummy_str, dummy_str2, number);
    if(dummy_str[0] == '6' && dummy_str[1] == '0') {
      dummy_str[0] = '0';
      dummy_str[1] = '0';
      minute += 1;
      if(minute == 60) {
	minute = 0;
	hour += 1;
      }
    }
    // account for negative dec
    if (hms[0] == '-')  {
      if(type == 2)
	sprintf(hms, "-%02dh%02dm%ss", hour, minute, dummy_str);
      else if(type == 3)
	sprintf(hms, "-%02dd%02d'%s\"", hour, minute, dummy_str);
      else if(type == 4)
	sprintf(hms, "-%02d%02d%s", hour, minute, dummy_str);
      else
	sprintf(hms, "-%02d:%02d:%s", hour, minute, dummy_str);
    }else {
      if(type == 2)
	sprintf(hms, "%02dh%02dm%ss", hour, minute, dummy_str);    
      else if(type == 3)
	sprintf(hms, "%02dd%02d'%s\"", hour, minute, dummy_str);    
      else if(type == 4)
	sprintf(hms, "%02d%02d%s", hour, minute, dummy_str);
      else
	sprintf(hms, "%02d:%02d:%s", hour, minute, dummy_str);    
    }
  }
}

//START REGION DEVELOP

/* Converts a time in hours into a string with hours, minutes and
   seconds. If precision is set to zero only minutes are produced. */
/*
void gethourstring(double t, char *txt, int precision)
{
  int th, tm, ts;
  if(t < 0) t += 24;
  if(t > 24) t -= 24;
  th = t;
  tm = 60*t - 60*th;
  ts = 60*60*t - 60*60*th - 60*tm;
  if(precision == 0) {
    if(ts >= 30)
      tm++;
    if(tm == 60) {
      tm = 0;
      th++;
    }
    if(th == 24)
      th = 0;
    sprintf(txt, "%02d:%02d", th, tm);
  }else {
    sprintf(txt, "%02d:%02d:%02d", th, tm, ts);
  }
}

*/

//START REGION DEVELOP
//START REGION RELEASE

/* Converts a string with hours, minutes and seconds into a floating
   point number (in hours).*/
void converthms(char *hms, double *h)
{
  int ret, i1, i2;
  double f1;
  ret = sscanf(hms, "%d:%d:%lf", &i1, &i2, &f1);

  *h = 0;
  if(ret >= 1) {
    *h += abs(i1);
    if(ret >= 2)
      *h += +abs(i2)/60.0;
    if(ret >= 3)
      *h += fabs(f1)/3600.0;
    if(hms[0] == '-')
      *h *= -1.0;
  }
}

//START REGION DEVELOP

/* Converts a string with hours, minutes and seconds into a floating
   point number (in hours).*/
void converthms_ld(char *hms, long double *h)
{
  int ret, i1, i2;
  long double f1;
  ret = sscanf(hms, "%d:%d:%Lf", &i1, &i2, &f1);

  *h = 0;
  if(ret >= 1) {
    *h += abs(i1);
    if(ret >= 2)
      *h += +abs(i2)/60.0;
    if(ret >= 3)
      *h += fabsl(f1)/3600.0;
    if(hms[0] == '-')
      *h *= -1.0;
  }
}

/* Converts a string with hours, minutes and seconds into a floating
   point number (in hours).*/
void converthms_f(char *hms, float *h)
{
  int ret, i1, i2;
  double f1;
  ret = sscanf(hms, "%d:%d:%lf", &i1, &i2, &f1);
  *h = 0;
  if(ret >= 1) {
    *h += abs(i1);
    if(ret >= 2)
      *h += +abs(i2)/60.0;
    if(ret >= 3)
      *h += fabs(f1)/3600.0;
    if(hms[0] == '-')
      *h *= -1.0;
  }
}


//START REGION DEVELOP
//START REGION RELEASE

/* Given the geodetic longitude and latitude of the observatory
   (longitude and latitude) and the position of the source (right
   ascension and declination) and the mjd, return the parallactic
   angle. All angles are in radians. If precess is set, the ra and dec
   will be corrected for precession, nutation and aberration (using
   J2000). */
double calc_parang(double longitude, double latitude, double ra, double dec, double mjd, int precess)
{
  double curlst, ha, parang;

  /*  csla_preces("FK5", 2000.0, csla_epj(mjd), &ra, &dec); */
  if(precess)
    calc_precess_nut_ab('J', mjd, &ra, &dec, 1, 1, 0);

  curlst = csla_gmst(mjd);
  ha = (curlst+longitude-ra);
  parang = csla_pa(ha, dec, latitude);
  return parang;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Set system to B or J for B1950 or J2000 coordinates. The
   coordinates are defined at 1950.0 or 2000.0 depending on the
   specified system and will be transformed to epoch mjd1. The input
   ra and dec (radians) will be transformed. If nutation is set, the
   additional effect of nutation is corrected for. Likewise aberration
   sets the correction for Annual aberration. Returns 0 if system was
   not recognized. */
int calc_precess_nut_ab(char system, double mjd, double *ra, double *dec, int nutation, int aberration, int verbose)
{
  double pm[3][3], v1[3], v2[3], ep0, ep1;
  double ra_new, dec_new, ra_new2, dec_new2, ra_new3, dec_new3;
  int useB;
  char txt[100];

  /* Validate sys */
  if(system == 'J' || system == 'j') {
    useB = 0;
  }else if(system == 'B' || system == 'b') {
    useB = 1;
  }else {
    return 0;
  }

  if(useB) {
    ep0 = 1950.0;
    ep1 = csla_epb(mjd);
  }else {
    ep0 = 2000.0;
    ep1 = csla_epj(mjd);
  }

  if(verbose) {
    printf("calc_precess_nut: Correct from epoch %lf to %lf\n", ep0, ep1);
  }

  /* Generate appropriate precession matrix for precession only*/
  if(useB) {
    csla_prebn(ep0, ep1, pm);
  }else {
    /* Use precl instead???? */
    csla_prec(ep0, ep1, pm);
  }

  /* Convert RA,Dec to x,y,z */
  csla_dcs2c(*ra, *dec, v1);

  /* Precess */
  csla_dmxv(pm, v1, v2);

  if(verbose) {
    /* Back to RA, Dec */
    csla_dcc2s(v2, &ra_new, &dec_new);
    ra_new = csla_dranrm(ra_new);
    converthms_string(txt, (ra_new-*ra)*12.0/M_PI, 2, 2);
    printf("  Delta RA  precession:        %9.6lf rad = %9.6lf deg = %s\n", ra_new-*ra, (ra_new-*ra)*180.0/M_PI, txt);
    converthms_string(txt, (dec_new-*dec)*180.0/M_PI, 2, 3);
    printf("  Delta DEC precession:        %9.6lf rad = %9.6lf deg = %s\n", dec_new-*dec, (dec_new-*dec)*180.0/M_PI, txt);
  }

  /* In J coordinates, precession and nutation can be done together. */
  if(useB == 0 && nutation) {
    csla_prenut(ep0, mjd, pm);
  }
  /* If B coordinates, we have to do nutation correction separately in sla. */
  if(useB && nutation) {
    v1[0] = v2[0];
    v1[1] = v2[1];
    v1[2] = v2[2];
    csla_nut(mjd, pm);
  }
  csla_dmxv(pm, v1, v2);

  /* Back to RA, Dec */
  csla_dcc2s(v2, &ra_new2, &dec_new2);
  ra_new2 = csla_dranrm(ra_new2);

  if(verbose && nutation) {
    converthms_string(txt, (ra_new2-ra_new)*12.0/M_PI, 2, 2);
    printf("  Delta RA  nutation:          %9.6lf rad = %9.6lf deg = %s\n", ra_new2-ra_new, (ra_new2-ra_new)*180.0/M_PI, txt);
    converthms_string(txt, (dec_new2-dec_new)*180.0/M_PI, 2, 3);
    printf("  Delta DEC nutation:          %9.6lf rad = %9.6lf deg = %s\n", dec_new2-dec_new, (dec_new2-dec_new)*180.0/M_PI, txt);
  }

  if(aberration) {
    double ebd[3], dpb[3], edh[3], eh[3], abv[3], p1dv, p1dvp1, w, ab1, vn[3], vm;
    v1[0] = v2[0];
    v1[1] = v2[1];
    v1[2] = v2[2];
    csla_evp(mjd, 2000.0, ebd, dpb, edh, eh);
    abv[0] = ebd[0]*499.004782;
    abv[1] = ebd[1]*499.004782;
    abv[2] = ebd[2]*499.004782;
    p1dv = csla_dvdv(v1, abv);
    p1dvp1 = p1dv+1.0;
    csla_dvn(abv, vn, &vm);
    ab1 = sqrt(1.0-vm*vm);
    w = 1.0+p1dv/(ab1+1.0);
    v2[0] = (ab1*v1[0]+w*abv[0])/p1dvp1;
    v2[1] = (ab1*v1[1]+w*abv[1])/p1dvp1;
    v2[2] = (ab1*v1[2]+w*abv[2])/p1dvp1;

    csla_dcc2s(v2, &ra_new3, &dec_new3);
    ra_new3 = csla_dranrm(ra_new3);

    if(verbose) {
      converthms_string(txt, (ra_new3-ra_new2)*12.0/M_PI, 2, 2);
      printf("  Delta RA  Annual aberration: %9.6lf rad = %9.6lf deg = %s\n", ra_new3-ra_new2, (ra_new3-ra_new2)*180.0/M_PI, txt);
      converthms_string(txt, (dec_new3-dec_new2)*180.0/M_PI, 2, 3);
      printf("  Delta DEC Annual aberration: %9.6lf rad = %9.6lf deg = %s\n", dec_new3-dec_new2, (dec_new3-dec_new2)*180.0/M_PI, txt);
    }
  }else {
    ra_new3 = ra_new2;
    dec_new3 = dec_new2;
  }

  *ra = ra_new3;
  *dec = dec_new3;
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE


// The following code came from TEMPO2:
// redwards functions for handling geodetic coordinates, which are
// needed for atomspheric corrections, and also may be inadvertantly
// provided on input
#define GRS80_A 6378137.0           /* semi-major axis (m) */
#define GRS80_F 1.0/298.257222101   /* flattening */
// Geocentric to geodetic.
// Uses Vermeille (2004)'s method:
//http://www.springerlink.com/app/home/contribution.asp?wasp=08ea5d2c4c62464789a7961196d84ab5&referrer=parent&backto=issue,11,18;journal,9,85;linkingpublicationresults,1:100435,1
//
// Longitude/latitude are GRS80 geodetic coordinates in radians
// Height and ITRF X,Y,Z coordinates are in meters
void tempo2_ITRF_to_GRS80(double obs_X, double obs_Y, double obs_Z, double *longitude, double *latitude, double *height)
{
  double p = (obs_X*obs_X + obs_Y*obs_Y)/ (GRS80_A*GRS80_A);
  double esq = GRS80_F*(2.0-GRS80_F);
  double q = (1.0-esq)/(GRS80_A*GRS80_A)*obs_Z*obs_Z;
  double r = (p+q-esq*esq)/6.0;
  double s = esq*esq*p*q/(4*r*r*r);
  double t = pow(1.0+s+sqrt(s*(2.0+s)), 1.0/3.0);
  double u = r*(1.0+t+1.0/t);
  double v = sqrt(u*u+esq*esq*q);
  double w = esq*(u+v-q)/(2.0*v);
  double k = sqrt(u+v+w*w)-w;
  double D = k*sqrt(obs_X*obs_X+obs_Y*obs_Y)/(k+esq);
  
  *height = (k+esq-1.0)/k * sqrt(D*D+obs_Z*obs_Z);
  *latitude = 2.0*atan2(obs_Z, D+sqrt(D*D+obs_Z*obs_Z));
  if (obs_Y >= 0.0)
    *longitude =
      0.5*M_PI - 2.0*atan2(obs_X, sqrt(obs_X*obs_X+obs_Y*obs_Y)+obs_Y);
  else
    *longitude = 
      -0.5*M_PI + 2.0*atan2(obs_X, sqrt(obs_X*obs_X+obs_Y*obs_Y)-obs_Y);
}

// Geodetic to geocentric: standard formula
//
// Longitude/latitude are GRS80 geodetic coordinates in radians
// Height and ITRF X,Y,Z coordinates are in meters
void
tempo2_GRS80_to_ITRF(double longitude, double latitude, double height, double *obs_X, double *obs_Y, double *obs_Z)
{
  double esq = GRS80_F * (2.0 - GRS80_F);
  double N = GRS80_A / sqrt(1.0-esq*sin(latitude)*sin(latitude));
  *obs_X = (N+height)*cos(latitude)*cos(longitude);
  *obs_Y = (N+height)*cos(latitude)*sin(longitude);
  *obs_Z = (N*(1.0-esq)+height)*sin(latitude);
}
 
// returns the geocentric (ITRF based) longitude of the observatory in radians
double observatory_long_geocentric(datafile_definition datafile)
{
  double lat;
  lat = atan2(datafile.telescope_Y, datafile.telescope_X);
  return lat;
}

// returns the geocentric (ITRF based) latitude of the observatory in radians
double observatory_lat_geocentric(datafile_definition datafile)
{
  double r, longitude;
  r = sqrt(datafile.telescope_X*datafile.telescope_X+datafile.telescope_Y*datafile.telescope_Y+datafile.telescope_Z*datafile.telescope_Z);
  longitude = asin(datafile.telescope_Z/r);
  return longitude;
}

// returns the geodetic (GRS80 based) longitude of the observatory in radians
double observatory_long_geodetic(datafile_definition datafile)
{
  double longitude, latitude, height;
  tempo2_ITRF_to_GRS80(datafile.telescope_X, datafile.telescope_Y, datafile.telescope_Z, &longitude, &latitude, &height);
  return longitude;
}

// returns the geodetic (GRS80 based) latitude of the observatory in radians
double observatory_lat_geodetic(datafile_definition datafile)
{
  double longitude, latitude, height;
  tempo2_ITRF_to_GRS80(datafile.telescope_X, datafile.telescope_Y, datafile.telescope_Z, &longitude, &latitude, &height);
  return latitude;
}

// returns the geodetic (GRS80 based) height of the observatory in meters
double observatory_height_geodetic(datafile_definition datafile)
{
  double longitude, latitude, height;
  tempo2_ITRF_to_GRS80(datafile.telescope_X, datafile.telescope_Y, datafile.telescope_Z, &longitude, &latitude, &height);
  return height;
}


//START REGION DEVELOP
