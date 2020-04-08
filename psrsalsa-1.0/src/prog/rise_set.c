#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "sla_wrap.h"
#include "psrsalsa.h"

int verbose;

void extractRADEC(char *psrname, double *alpha, double *delta, int switch_galactic)
{
  double u1;
  char psr[200];
  if(strcasecmp(psrname, "crab") == 0) {
    *alpha = 5 + (34.0+31.973/60.0)/60.0;
    *delta = 22+(0.0+52.06/60.0)/60.0;
    return;
  }else if(strcasecmp(psrname, "vela") == 0) {
    *alpha = 8 + (35.0+20.61149/60.0)/60.0;
    *delta = -45-(10.0+34.8751/60.0)/60.0;
    return;
  }else if(strcasecmp(psrname, "M31") == 0) {
    *alpha = 0 + (42+44/60.0)/60.0;
    *delta = 41+(16+9/60.0)/60.0;
    return;
  }else if(strcasecmp(psrname, "M33") == 0) {
    *alpha = 1 + (33+50/60.0)/60.0;
    *delta = 30+(39+37/60.0)/60.0;
    return;
  }else if(strcasecmp(psrname, "CASA") == 0) {
    *alpha = 23 + (23+26.4/60.0)/60.0;
    *delta = 58+(49+34.0/60.0)/60.0;
    return;
  }else if(strcasecmp(psrname, "CYGA") == 0) {
    *alpha = 19 + (59+29.2/60.0)/60.0;
    *delta = 40+(44+16.0/60.0)/60.0;
    return;
  }else if(strcasecmp(psrname, "TAUA") == 0) {
    *alpha = 05 + (34+31.0/60.0)/60.0;
    *delta = 22+(00+52.1/60.0)/60.0;
    return;
  }else if(strcasecmp(psrname, "VIRA") == 0) {
    *alpha = 12 + (30+49.6/60.0)/60.0;
    *delta = 12+(23+21.0/60.0)/60.0;
    return;
  }else if(strcasecmp(psrname, "HERA") == 0) {
    *alpha = 16 + (51+08.3/60.0)/60.0;
    *delta = 04+(59+26.0/60.0)/60.0;
    return;
  }else if(strcasecmp(psrname, "3C353") == 0) {
    *alpha = 17 + (20+29.5/60.0)/60.0;
    *delta = 00-(58+52.0/60.0)/60.0;
    return;
  }
  strcpy(psr, psrname);
  if(psr[0] == 'B' || psr[0] == 'J')
    strcpy(psr, &psr[1]);
  if(verbose) printf("Pulsar: %s\n", psr);
  if(strlen(psr) == 9) {
    sscanf(&psr[7], "%lf", delta);
    (*delta) /= 60.0;
  }else
    (*delta) = 0;
  psr[7] = 0;
  sscanf(&psr[5], "%lf", &u1);
  (*delta) += u1;
  if(psr[4] == '-')
    (*delta) *= -1;
  psr[4] = 0;
  
  sscanf(&psr[2], "%lf", alpha);
  if(switch_galactic == 0)
    (*alpha) /= 60;
  psr[2] = 0;
  sscanf(psr, "%lf", &u1);
  if(switch_galactic == 0)
    *alpha += u1;
  else
    *alpha += 100*u1;
  if(switch_galactic)
    *alpha *= 12.0/180.0;
}

int catalogRADEC(char *psrname, char *RAJ, char *DECJ)
{
  char command[200], junk[100];
  FILE *fin;
  int i;

  sprintf(command, "psrcat -nonumber -all -c \"RAJ DECJ P0\" %s", psrname);
  strcat(command, " > junk");
  if(verbose)
    printf("%s\n", command);
  system("rm -f junk");
  system(command);
  /*    if(verbose)  system("more junk"); */
  if((fin = fopen("junk", "r")) == NULL) {
    printf("Error: no file 'junk'\n");
    return 0;
  }
  for(i = 0; i < 8; i++) {
    fscanf(fin, "%s", junk);
  }
  fscanf(fin, "%s", RAJ);
  fscanf(fin, "%s", junk);
  fscanf(fin, "%s", junk);
  fscanf(fin, "%s", DECJ);
  fclose(fin);
  
  if(verbose) {
    printf("RAJ:                  %s\n", RAJ);
    printf("DECJ:                 %s\n", DECJ);
  }
  system("rm -f junk");
  return 1;
}

int main(int argc, char* argv[])
{
  /* Default location is WSRT, take East to be positive */
  double latitude = M_PI*(52.731622)/180.0;
  double elevation_limit = M_PI*(5)/180.0;
  double longitude = M_PI*(6.604457)/180.0;
  double alpha, delta, rise, set, curmjd, curlst, gall, galb, ha, az, el, parallacticAngle; 
  double alpha_prec, delta_prec, dtime;
  int i, ret, UseCatalog, switch_galactic, pos_set, elevation_limit_manual;
  char RAJ[100], DECJ[100], telescopename[100], txt[100];
  time_t rawtime;
  int elcurve;
  struct tm * ptm;

  UseCatalog = 0;
  verbose = 0;
  curmjd = sqrt(-1);
  switch_galactic = 0;
  pos_set = 0;
  elevation_limit_manual = 0;
  elcurve = 0;
  strcpy(telescopename, "WSRT");

  /* Take current time/date as the UT time */

  time ( &rawtime );
  ptm = gmtime ( &rawtime );
  /*  printf("XXX %d %d %d\n", ptm->tm_year, ptm->tm_mon, ptm->tm_mday); */
  csla_cldj(1900+ptm->tm_year, ptm->tm_mon+1, ptm->tm_mday, &curmjd, &ret);

  //  printf("XXXXXX %lf\n", curmjd);
  if(ret != 0) {
    fprintf(stderr, "ERROR rise_set: cannot calculate mjd.\n");
    return 0;
  }
  curmjd += (ptm->tm_hour + ptm->tm_min/60.0 + ptm->tm_sec/3600.0)/24.0;


  if(argc >= 1) {
    for(i = 0; i < argc; i++) {
      if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-H") == 0 || argc == 1) {
	fprintf(stdout, "Usage: rise_set psrname\n");
	fprintf(stdout, "\n-p        Use sky location of prscat rather than pulsar name to determine RA and DEC.\n");
	fprintf(stdout, "-v        Verbose.\n");
	fprintf(stdout, "-jb       Assume Jodrell Bank location instead of default location (WSRT).\n");
	fprintf(stdout, "-pks      Assume Parkes location instead of default location (WSRT).\n");
	fprintf(stdout, "-fast     Assume FAST location instead of default location (WSRT).\n");
	fprintf(stdout, "-ao       Assume Arecibo location instead of default location (WSRT).\n");
	fprintf(stdout, "-chil     Assume Chilbolton LOFAR location instead of default location (WSRT).\n");
	fprintf(stdout, "-nenu     Assume NenuFAR location at Nancay observatory instead of default location (WSRT).\n");
	fprintf(stdout, "-loc      Give location \"longitude latitude\" in degrees (geodetic).\n");
	fprintf(stdout, "-mjd      Set the current mjd.\n");
	fprintf(stdout, "-cal      Set the current time \"year month day hour:minute:second\" (UT).\n");
	fprintf(stdout, "-uttime   Set the current time \"hour:minute:second\" (UT), but take the current date.\n");
	fprintf(stdout, "-pos      Set the source position in degrees (source name will be ignored).\n");
	fprintf(stdout, "-sun      Set the source position to that of the Sun.\n");
	fprintf(stdout, "-el       Set the elevation limit of the telescope.\n");
	fprintf(stdout, "-gal      The pulsar name or -pos coordinates are in galactic coordinates rather than RA and DEC.\n");
	fprintf(stdout, "-elcurve  Show the elevation as function of UT.\n");
	fprintf(stdout, "-h        This help.\n");
	fprintf(stdout, "\n\nCalculates when source is above elevation limit\nand within slew range.\n");
	return 0;
      }else if(strcmp(argv[i], "-v") == 0) {
	verbose = 1;   
      }else if(strcmp(argv[i], "-jb") == 0) {
	latitude = M_PI*(53.233901)/180.0;
	longitude = -M_PI*(2+18/60.0 + 25.74/3600.0)/180.0;
	strcpy(telescopename, "JB");
      }else if(strcmp(argv[i], "-chil") == 0) {
	latitude = M_PI*(51.143699)/180.0;    // From Google maps
	longitude = -M_PI*(1.433990)/180.0;
	strcpy(telescopename, "Chilbolton");
      }else if(strcmp(argv[i], "-nenu") == 0) {
	latitude = M_PI*(47.376676)/180.0;    // From Google maps
	longitude = M_PI*(2.192193)/180.0;
	strcpy(telescopename, "NenuFAR");
      }else if(strcmp(argv[i], "-pks") == 0) {
	latitude = M_PI*(-32.99994444444444444444)/180.0;
	longitude = M_PI*(148+15/60.0+48.636/3600.0)/180.0;
	if(elevation_limit_manual == 0)
	  elevation_limit = M_PI*(30.5)/180.0;
	strcpy(telescopename, "Parkes");
      }else if(strcmp(argv[i], "-fast") == 0) {
	latitude = M_PI*(25+39/60.0+10.6/3600.0)/180.0;
	longitude = M_PI*(106+51/60.0+24/3600.0)/180.0;
	if(elevation_limit_manual == 0)
	  elevation_limit = M_PI*(30.5)/180.0;
	strcpy(telescopename, "FAST");
      }else if(strcmp(argv[i], "-ao") == 0) {
	latitude = M_PI*(18+20/60.0+36.6/3600.0)/180.0;
	longitude = -M_PI*(66+45/60.0+11.1/3600.0)/180.0;
	if(elevation_limit_manual == 0)
	  elevation_limit = M_PI*(90.0-19.7)/180.0;
	strcpy(telescopename, "Arecibo");
      }else if(strcmp(argv[i], "-loc") == 0) {
	ret = sscanf(argv[i+1], "%lf %lf", &longitude, &latitude);
	if(ret != 2) {
	  fprintf(stderr, "ERROR rise_set: cannot parse %s option. Expected 2 floating point values.\n", argv[i]);
	  return 0;
	}
	i++;
	latitude *= M_PI/180.0;
	longitude *= M_PI/180.0;
	strcpy(telescopename, "User defined location");
      }else if(strcmp(argv[i], "-mjd") == 0) {
	ret = sscanf(argv[i+1], "%lf", &curmjd);
	if(ret != 1) {
	  fprintf(stderr, "ERROR rise_set: cannot parse %s option. Expected 1 floating point value.\n", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-pos") == 0) {
	ret = sscanf(argv[i+1], "%lf %lf", &alpha, &delta);
	if(ret != 2) {
	  fprintf(stderr, "ERROR rise_set: cannot parse %s option. Expected 2 floating point value.\n", argv[i]);
	  return 0;
	}
	pos_set = 1;
	i++;
      }else if(strcmp(argv[i], "-sun") == 0) {
	pos_set = 2;
      }else if(strcmp(argv[i], "-cal") == 0) {
	int year, month, day;
	double hour, minute, second;
	ret = sscanf(argv[i+1], "%d %d %d %lf:%lf:%lf", &year, &month, &day, &hour, &minute, &second);
	if(ret != 6) {
	  fprintf(stderr, "ERROR rise_set: cannot parse %s option. Expected 6 floating point values.\n", argv[i]);
	  return 0;
	}
	/*	fprintf(stderr, "%d %d %d %lf:%lf:%lf\n", year, month, day, hour, minute, second); */
	csla_cldj(year, month, day, &curmjd, &ret);
	if(ret != 0) {
	  fprintf(stderr, "ERROR rise_set: cannot calculate mjd.\n");
	  return 0;
	}
	curmjd += (hour + minute/60.0 + second/3600.0)/24.0;
	i++;
      }else if(strcmp(argv[i], "-uttime") == 0 || strcmp(argv[i], "-time") == 0) {
	double hour, minute, second;
	ret = sscanf(argv[i+1], "%lf:%lf:%lf", &hour, &minute, &second);
	if(ret != 3) {
	  fprintf(stderr, "ERROR rise_set: cannot parse %s option. Expected 3 floating point values.\n", argv[i]);
	  return 0;
	}
	/*	fprintf(stderr, "%lf:%lf:%lf\n", hour, minute, second); */
	curmjd = floor(curmjd);
	curmjd += (hour + minute/60.0 + second/3600.0)/24.0;
	i++;
      }else if(strcmp(argv[i], "-el") == 0) {
	ret = sscanf(argv[i+1], "%lf", &elevation_limit);
	if(ret != 1) {
	  fprintf(stderr, "ERROR rise_set: cannot parse %s option. Expected 1 floating point value.\n", argv[i]);
	  return 0;
	}
	elevation_limit *= M_PI/180.0;
	elevation_limit_manual = 1;
	i++;
      }else if(strcmp(argv[i], "-p") == 0) {
	UseCatalog = 1;   
      }else if(strcmp(argv[i], "-gal") == 0) {
	switch_galactic = 1;   
      }else if(strcmp(argv[i], "-elcurve") == 0) {
	elcurve = 1;
      }else {
	if(i != argc - 1 && i != 0) {
	  fprintf(stderr, "ERROR rise_set: Unrecognized option: %s\n", argv[i]);
	  return 0;
	}
      }
    }
  }

  if(verbose) {
    printf("Geodetic location of %s: %lf deg ", telescopename, fabs(longitude*180.0/M_PI));
    if(longitude >= 0)
      printf("East");
    else
      printf("West");
    printf(" and %lf deg ", fabs(latitude*180.0/M_PI));
    if(latitude >= 0)
      printf("North\n");
    else
      printf("South\n");
  }
  if(pos_set == 0) {
    extractRADEC(argv[argc-1], &alpha, &delta, switch_galactic);
    if(UseCatalog) {
      if(!catalogRADEC(argv[argc-1], RAJ, DECJ)) {
	return 0;
      }
      converthms(RAJ, &alpha);
      converthms(DECJ, &delta);
    }
    alpha *= M_PI/12.0;
    delta *= M_PI/180.0;
  }else if(pos_set == 2) {
    if(isnan(curmjd)) {
      fprintf(stderr, "Error rise_set: Date and time are not set\n");
      return 0;
    }
    double dvb[3], dpb[3], dvh[3], dph[3], ecl_1, ecl_2;
    csla_evp(curmjd, csla_epj(curmjd), dvb, dpb, dvh, dph);
    printf("Barycentric  X Y Z: %lf %lf %lf AU (epoch %lf)\n", dpb[0], dpb[1], dpb[2], csla_epj(curmjd));
    printf("Heliocentric X Y Z: %lf %lf %lf AU (epoch %lf)\n", dph[0], dph[1], dph[2], csla_epj(curmjd));
    csla_evp(curmjd, 2000.0, dvb, dpb, dvh, dph);
    printf("Barycentric  X Y Z: %lf %lf %lf AU (epoch 2000)\n", dpb[0], dpb[1], dpb[2]);
    printf("Heliocentric X Y Z: %lf %lf %lf AU (epoch 2000)\n", dph[0], dph[1], dph[2]);
    csla_dcc2s(dph, &ecl_1, &ecl_2);
    /*    printf("XXXXXX %lf %lf\n", ecl_1*180.0/M_PI, ecl_2*180.0/M_PI); */
    alpha = csla_dranrm(ecl_1+M_PI);
    delta = -ecl_2;
  }else {
    alpha *= M_PI/180.0;
    delta *= M_PI/180.0;
  }

  if(switch_galactic) {
    /*    if(verbose) printf("GALL  = %lf h (%lf deg)\n", alpha*12.0/M_PI, alpha*180.0/M_PI);
	  if(verbose) printf("GALB  = %lf deg\n", delta*180/M_PI);*/
    csla_galeq(alpha, delta, &alpha, &delta);
  }

  converthms_string(txt, alpha*12.0/M_PI, 2, 2);
  printf("RAJ  = %s", txt);
  if(verbose) 
    printf(" (%lf h or %lf deg or %lf rad)\n", alpha*12.0/M_PI, alpha*180.0/M_PI, alpha);
  else
    printf("     ");
  converthms_string(txt, delta*180.0/M_PI, 2, 3);
  printf("DECJ = %s", txt);
  if(verbose)
    printf(" (%lf deg or %lf rad)\n", delta*180/M_PI, delta);
  else
    printf("\n");

  if(verbose) {
    csla_eqgal(alpha, delta, &gall, &galb);
    converthms_string(txt, gall*180.0/M_PI, 2, 3);
    printf("Galactic l = %s (%lf deg or %lf h or %lf rad)\n", txt, gall*180.0/M_PI, gall*12.0/M_PI, gall);
    converthms_string(txt, galb*180.0/M_PI, 2, 3);
    printf("Galactic b = %s (%lf deg or %lf rad)\n", txt, galb*180/M_PI, galb);
  }

  if(elcurve) {
    long curmjd_int;
    long nr_points_curve = 1000;
    curmjd_int = curmjd;
    float *mjd_curve, *el_curve;
    mjd_curve = malloc(nr_points_curve*sizeof(float));
    el_curve = malloc(nr_points_curve*sizeof(float));
    if(mjd_curve == NULL || el_curve == NULL) {
      fprintf(stderr, "Memory allocation error\n");
      return 0;
    }
    double alpha_orig, delta_orig;
    alpha_orig = alpha;
    delta_orig = delta;
    for(i = 0; i < nr_points_curve; i++) {
      curmjd = curmjd_int + i*2.0/(double)(nr_points_curve-1); 
      alpha = alpha_orig;
      delta = delta_orig;
      alpha_prec = alpha_orig;
      delta_prec = delta_orig;
      if(calc_precess_nut_ab('j', curmjd, &alpha_prec, &delta_prec, 1, 1, 0) == 0) {
	fprintf(stderr, "Precession & nutation correction failed\n");
      }else {
	alpha = alpha_prec;
	delta = delta_prec;
      }

      curlst = csla_gmst(curmjd);
      ha = (curlst+longitude-alpha);
      csla_de2h(ha, delta, latitude, &az, &el);

      mjd_curve[i] = curmjd-curmjd_int;
      el_curve[i] = el*180.0/M_PI;
      //      printf("%lf %lf\n", curmjd, el*180.0/M_PI);
    }
    pgplot_options_definition pgplot;
    pgplot_clear_options(&pgplot);
    pgplot.viewport.dontclose = 1;
    sprintf(pgplot.box.xlabel, "MJD-%ld", curmjd_int);
    strcpy(pgplot.box.ylabel, "Elevation");
    sprintf(pgplot.box.title, "Elevation for %s", telescopename);
    verbose_definition verbose;
    cleanVerboseState(&verbose);
    pgplotGraph1(&pgplot, el_curve, mjd_curve, NULL, nr_points_curve, -1, -1, 0, -1, -1, -1, -1, 0, 0, 0, 1, 0, 1, 1, NULL, -1, verbose);
    ppgmove(0, elevation_limit*180.0/M_PI);
    ppgsci(2);
    ppgsls(2);
    ppgdraw(2, elevation_limit*180.0/M_PI);
    ppgend();
    free(mjd_curve);
    free(el_curve);
  }else if(!isnan(curmjd)) {
    mjd2dateString(curmjd, txt, 2, 1, " ");
    printf("MJD:  %lf (JD %lf)\n", curmjd, curmjd+2400000.5);
    printf("UT:   %s (epoch %lf)\n", txt, csla_epj(curmjd));
    alpha_prec = alpha;
    delta_prec = delta;
    //    csla_preces("FK5", 2000.0, csla_epj(curmjd), &alpha_prec, &delta_prec);
    if(calc_precess_nut_ab('j', curmjd, &alpha_prec, &delta_prec, 1, 1, verbose) == 0) {
      fprintf(stderr, "Precession & nutation correction failed\n");
    }else {
      if(verbose) {
	printf("After precession, nutation and abberation: \n");
	converthms_string(txt, (alpha_prec-alpha)*12.0/M_PI, 2, 2);
	printf("  Delta RAJ  = %s (%lf deg)\n", txt, (alpha_prec-alpha)*180.0/M_PI);
	converthms_string(txt, (delta_prec-delta)*180.0/M_PI, 2, 3);
	printf("  Delta DECJ = %s (%lf deg)\n", txt, (delta_prec-delta)*180.0/M_PI);
	converthms_string(txt, alpha_prec*12.0/M_PI, 2, 2);
	printf("  RAJ  = %s", txt);
	printf(" (%lf h or %lf deg or %lf rad)\n", alpha_prec*12.0/M_PI, alpha_prec*180.0/M_PI, alpha_prec);
	converthms_string(txt, delta_prec*180.0/M_PI, 2, 3);
	printf("  DECJ = %s", txt);
	printf(" (%lf deg or %lf rad)\n", delta_prec*180/M_PI, delta_prec);
      }

      /*
      csla_map(alpha, delta, 0, 0, 0, 0, 2000, curmjd, &alpha_prec, &delta_prec);
      converthms_string(txt, alpha_prec*12.0/M_PI, 2, 2);
      printf("  RAJ  = %s", txt);
      printf(" (%lf h or %lf deg or %lf rad)\n", alpha_prec*12.0/M_PI, alpha_prec*180.0/M_PI, alpha_prec);
      converthms_string(txt, delta_prec*180.0/M_PI, 2, 3);
      printf("  DECJ = %s", txt);
      printf(" (%lf deg or %lf rad)\n", delta_prec*180/M_PI, delta_prec);
      */

      alpha = alpha_prec;
      delta = delta_prec;
    }

    /* The MJD should be UT1, not UTC */
    curlst = csla_gmst(curmjd);
    ha = (curlst+longitude-alpha);
    curlst *= 12.0/M_PI;
    converthms_string(txt, curlst, 2, 1);
    printf("GMST: %s\n", txt);

    /* Current hour angle */

    curlst += longitude*12.0/M_PI;
    converthms_string(txt, curlst, 2, 1);
    dtime = curlst - (curmjd - floor(curmjd))*24.0;
    printf("LST:  %s (%s) = UT ", txt, telescopename);
    if(dtime > 12)
      dtime -= 24.0;
    if(dtime < 0) {
      printf("-");
    }else {
      printf("+");
    }
    converthms_string(txt, fabs(dtime), 2, 1);
    printf(" %s hours\n", txt);

    csla_de2h(ha, delta, latitude, &az, &el);
    parallacticAngle = calc_parang(longitude, latitude, alpha, delta, curmjd, 0);

    printf("Current Hour angle:        %lf deg\n", ha*180/M_PI);
    printf("Current parallactic angle: %lf deg\n", parallacticAngle*180/M_PI);
    printf("Current Azimuth:           %lf deg\n", az*180/M_PI);
    printf("Current Elevation:         %lf deg\n", el*180/M_PI);
  }
  if(verbose) printf("Elevation limit: %.1lf deg\n", elevation_limit*180.0/M_PI);
  i = riseset(alpha, delta, elevation_limit, latitude, &rise, &set, verbose);
  if(i == 1) {
    int time_type;
    for(time_type = 0; time_type < 2; time_type++) {
      if(time_type == 1) {
	rise -= dtime*M_PI/12.0;
	set -= dtime*M_PI/12.0;
      }
      printf("Visible from %s between: ", telescopename);
      printhour(rise);
      printf(" and ");
      printhour(set);
      switch(time_type) {
      case 0: printf(" LST\n"); break;
      case 1: printf(" UT\n"); break;
      }
    }
  }else if(i == 0) {
    printf("Source is not visible from %s\n", telescopename);
  }else {
    printf("Source is always visible from %s\n", telescopename);
  }

  return 0;
}
