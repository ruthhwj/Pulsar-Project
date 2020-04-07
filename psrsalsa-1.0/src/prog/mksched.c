#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "psrsalsa.h"

#define MaxNrSources 1000
#define MaxNrSchedules 5
#define MaxNrLookupEntries 3000

double latitude = M_PI*(-32.99994444444444444444)/180.0;
double elevation_limit = M_PI*(30.5)/180.0;
double slewtime =     120.0;

typedef struct {
  char name[25], commend[500], risesetstr[30], lst0str[10], linenrstr[10];
  double raj, decj, rise, set, az, zen, lst0;
  long duration, waitaftersource;
  int docal, dolevcal, fbsamp, calaftersource;
  float edot, dist, rms, mjd;
}source_type;

typedef struct {
  char name[25];
  double raj, decj;
  float edot, dist, rms, mjd;
  int fbsamp;
}lookuptable_entry;


void extractRADEC(char *psr, double *alpha, double *delta);
int catalogRADEC(char *psrname, char *RAJ, char *DECJ);
/* void converthms(char *hms, double *h);
   int riseset(double alpha, double delta, double *rise, double *set); */
/*
void calculateAZ(double raj, double decj, double h, double *az, double *zen);
void printhour(double t);
void gethourstring(double t, char *txt, int precision);
*/
int readSourceList(source_type *sources, int *nrsources, char *filename);
int catalogEdotDist(char *psrname, float *edot, float *dist);
int readlookuptable(char *filename, lookuptable_entry *entries, int *nrlookup_entries, float *latestmjd);
void lookup(source_type *source, lookuptable_entry *entries, int nrlookup_entries);
int readrmstable(char *filename, lookuptable_entry *entries, int nrlookup_entries);
int readfbtable(char *filename, lookuptable_entry *entries, int nrlookup_entries);
int FindLastObservation(char *psrname, float *mjd);



void help()
{
  fprintf(stderr, "Usage: mksched sourcelist\n\n");
  fprintf(stderr, "-p     Use sky location of prscat.\n");
  fprintf(stderr, "-d     Dump lookup table.\n");
  fprintf(stderr, "-l     Use this lookup table instead of using psrcat.\n");
  fprintf(stderr, "-rms   Use this rms lookup table.\n");
  fprintf(stderr, "-fbdt  Use sampling time lookup table.\n");
  fprintf(stderr, "-slew  Calculate slew times. Seems to work best with -t 100 option.\n");
  fprintf(stderr, "-t     Specify extra time per source (def=%.2f sec).\n", slewtime);
  fprintf(stderr, "-v     Verbose.\n");
  fprintf(stderr, "-h     Help.\n");
  fprintf(stderr, "-psr   Calculate position of this pulsar at specified lst\n");
  fprintf(stderr, "-short Print out the short version.\n");
  fprintf(stderr, "-skip  Skip this number of entries in sourcelist.\n");
  fprintf(stderr, "-old   Run program in old tcs mode instead of the next generation version.\n");
  fprintf(stderr, "-gnuplot  Generate a file called sky.gnu.\n");
}

int verbose;

int main(int argc, char* argv[])
{
  int i, j, UseCatalog, psrmode, nrsources, linenr, dumpmode, calc_slew, fb_dt_to_use, gnuplot;
  long cal_separation;
  char RAJ[100], DECJ[100], txt[100], txtcal[100], lookuptable[1000], rmstable[1000], fbsamptable[1000];
  double lst, lst_lastcal, timemultiplier, lst_lastlevcal, t1, t2;
  source_type sources[MaxNrSources];
  int nrschedules, levcal_wait, schedule_levcal[MaxNrSchedules], dfbnr[MaxNrSchedules], srchmode[MaxNrSchedules];
  char recv[MaxNrSchedules][100], config_file[MaxNrSchedules][100];
  float freq[MaxNrSchedules], timeconst[MaxNrSchedules];
  int wbtsub, wbtsub_skip, wbtsub_cal, tobs_cal, tobs_srch_levcal, feed_rot, usefilterbank, fbsamp, usefbsamptable;
  /* int nrfbsamptable*/
  char fbconfig[100];
  FILE *fout[MaxNrSchedules], *foutf, *gnuplotf;
  int nrlookup_entries, extendedprintout, nrentriesskip, nextgen;
  float dAz, dZen, tAz, tZen, dt, latestmjd;
  lookuptable_entry lookuptable_entrys[MaxNrLookupEntries];
  char gnuplot_buffer[10000];

  nextgen = 1;
  UseCatalog = 0;
  verbose = 0;
  psrmode = 0;
  dumpmode = 0;
  lookuptable[0] = 0;
  rmstable[0] = 0;
  calc_slew = 0;
  usefbsamptable = 0;
  extendedprintout = 1;
  nrentriesskip = 0;
  gnuplot = 0;
  gnuplot_buffer[0] = 0;
  if(argc < 2) {
    help();
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-H") == 0) {
	help();
	return 0;
      }else if(strcmp(argv[i], "-v") == 0) {
	verbose = 1;   
      }else if(strcmp(argv[i], "-old") == 0) {
	nextgen = 0;
      }else if(strcmp(argv[i], "-d") == 0) {
	dumpmode = 1;   
      }else if(strcmp(argv[i], "-p") == 0) {
	UseCatalog = 1;   
      }else if(strcmp(argv[i], "-slew") == 0) {
	calc_slew = 1;   
      }else if(strcmp(argv[i], "-short") == 0) {
	extendedprintout = 0;
      }else if(strcmp(argv[i], "-psr") == 0) {
	psrmode = 1;
	j = sscanf(argv[++i], "%s", txt);
	if(j != 1) {
	  printf("Cannot parse %s option.\n", argv[i-1]);
	  return 0;
	}
	converthms(txt, &lst);
      }else if(strcmp(argv[i], "-l") == 0) {
	j = sscanf(argv[++i], "%s", lookuptable);
	if(j != 1) {
	  printf("Cannot parse %s option.\n", argv[i-1]);
	  return 0;
	}
      }else if(strcmp(argv[i], "-rms") == 0) {
	j = sscanf(argv[++i], "%s", rmstable);
	if(j != 1) {
	  printf("Cannot parse %s option.\n", argv[i-1]);
	  return 0;
	}
      }else if(strcmp(argv[i], "-fbdt") == 0) {
	j = sscanf(argv[++i], "%s", fbsamptable);
	if(j != 1) {
	  printf("Cannot parse %s option.\n", argv[i-1]);
	  return 0;
	}
	usefbsamptable = 1;
      }else if(strcmp(argv[i], "-t") == 0) {
	j = sscanf(argv[++i], "%lf", &slewtime);
	if(j != 1) {
	  printf("Cannot parse %s option.\n", argv[i-1]);
	  return 0;
	}
      }else if(strcmp(argv[i], "-skip") == 0) {
	j = sscanf(argv[++i], "%d", &nrentriesskip);
	if(j != 1) {
	  printf("Cannot parse %s option.\n", argv[i-1]);
	  return 0;
	}
      }else if(strcmp(argv[i], "-gnuplot") == 0) {
	gnuplot = 1;
      }else {
	if(i == argc - 1)
	  break;
	fprintf(stderr, "Unrecognized option: %s\n", argv[i]);
	return 0;
      }
    }
  }

  if(psrmode) {
    /* Get RA and DEC */
    extractRADEC(argv[argc-1], &sources[0].raj, &sources[0].decj);
    if(UseCatalog) {
      if(!catalogRADEC(argv[argc-1], RAJ, DECJ)) {
	return 0;
      }
      converthms(RAJ, &sources[0].raj);
      converthms(DECJ, &sources[0].decj);
    }
    if(verbose) printf("alpha = %lf h\n", sources[0].raj);
    if(verbose) printf("delta = %lf deg\n", sources[0].decj);
    sources[0].raj *= M_PI/12;
    sources[0].decj *= M_PI/180.0;

    riseset(sources[0].raj, sources[0].decj, elevation_limit, latitude, &sources[0].rise, &sources[0].set, verbose);
    printf("Visible for Parkes between: ");
    printhour(sources[0].rise);
    printf(" and ");
    printhour(sources[0].set);
    printf("\n");

    lst *= M_PI/12.0;
    calculate_az_zen(sources[0].raj, sources[0].decj, lst, latitude, &sources[0].az, &sources[0].zen);
    printf("A = %lf Z=%lf\n", sources[0].az*180.0/M_PI, sources[0].zen*180.0/M_PI);
    return 0;
  }

  if(dumpmode) {
    if(!readSourceList(sources, &nrsources, argv[argc-1])) {
      return 0;
    }else {
      for(i = 0; i < nrsources; i++) {
	extractRADEC(sources[i].name, &sources[i].raj, &sources[i].decj);
	if(!catalogRADEC(sources[i].name, RAJ, DECJ)) {
	  return 0;
	}
	converthms(RAJ, &sources[i].raj);
	converthms(DECJ, &sources[i].decj);
	if(!catalogEdotDist(sources[i].name, &sources[i].edot, &sources[i].dist)) {
	  return 0;
	}
	if(!FindLastObservation(sources[i].name, &sources[i].mjd)) {
	  return 0;
	}
	printf("%s %e %e %e %e %f\n", sources[i].name, sources[i].raj, sources[i].decj, sources[i].edot, sources[i].dist, sources[i].mjd);
      }
    }
    return 0;
  }

  latestmjd = 0;
  if(lookuptable[0] != 0) {
    if(!readlookuptable(lookuptable, lookuptable_entrys, &nrlookup_entries, &latestmjd))
      return 0;
    fprintf(stderr, "Found %d entries in lookup table\n", nrlookup_entries);
  }
  if(rmstable[0] != 0) {
    if(!readrmstable(rmstable, lookuptable_entrys, nrlookup_entries))
      return 0;
    fprintf(stderr, "Found %d entries in rms lookup table\n", nrlookup_entries);
  }


  fprintf(stderr, "Start LST: ");
  scanf("%s", txt);
  converthms(txt, &lst);
  fprintf(stderr, "lst = %lf\n", lst);
  fprintf(stderr, "Multiply durations with this factor: ");
  scanf("%lf", &timemultiplier);
  fprintf(stderr, "Time between cals (0 = always): ");
  scanf("%ld", &cal_separation);
  fprintf(stderr, "Specify time between levcal's (0 = always, -1 = never): ");
  scanf("%d", &levcal_wait);
  fprintf(stderr, "How many schedules files: ");
  scanf("%d", &nrschedules);
  if(nrschedules > 0) {
    for(i = 0; i < nrschedules; i++) {
      fprintf(stderr, "For schedule %d: ", i+1);
      fprintf(stderr, "Specify reciever (e.g. MULTI): ");
      scanf("%s", recv[i]);
      fprintf(stderr, "Specify config file (e.g. wbb1024_256_512_3p_c): ");
      scanf("%s", config_file[i]);
      if(config_file[i][0] == 'd') {
	dfbnr[i] = 1;
      }else {
	if(config_file[i][0] == 'p') {
	  if(config_file[i][4] == '2') {
	    dfbnr[i] = 2;
	  }else if(config_file[i][4] == '3') {
	    dfbnr[i] = 3;
	  }else if(config_file[i][4] == '4') {
	    dfbnr[i] = 4;
	  }else {
	    dfbnr[i] = -1;
	  }
	}else {
	  dfbnr[i] = -1;
	}
      }
      if(dfbnr[i] < 1 || dfbnr[i] > 4) {
	fprintf(stdout, "Cannot guess dfb nr from specified configuration, give dfb number: ");
	fflush(stdout);
	scanf("%d", &dfbnr[i]);
      }
      srchmode[i] = 0;
      if(config_file[i][0] == 's' && config_file[i][1] == 'r' && config_file[i][2] == 'c' && config_file[i][3] == 'h') {
	printf("The config file suggests this is a seach mode observation\n");
	srchmode[i] = 1; 
	fprintf(stderr, "Baseline time constant (e.g. 1.0): ");
	scanf("%f", &timeconst[i]);
	fprintf(stderr, "tobs_levcal (e.g. 10 sec): ");
	scanf("%d", &tobs_srch_levcal);
      }
      fprintf(stderr, "Specify frequecy (e.g. 1369.0): ");
      scanf("%f", &freq[i]);
      fprintf(stderr, "Schedule levcal's (0 = no 1 = yes): ");
      scanf("%d", &schedule_levcal[i]);
    }
  }
  fprintf(stderr, "Specify wbtsub (e.g. 30): ");
  scanf("%d", &wbtsub);
  fprintf(stderr, "Specify wbtsub's to skip (e.g. 0): ");
  scanf("%d", &wbtsub_skip);
  fprintf(stderr, "Specify wbtsub_cal (e.g. 20): ");
  scanf("%d", &wbtsub_cal);
  fprintf(stderr, "tobs_cal (e.g. 60 sec): ");
  scanf("%d", &tobs_cal);
  fprintf(stderr, "Do you want to use feed rotation? (0 = no, 1 = yes): ");
  scanf("%d", &feed_rot);
  fprintf(stderr, "Do you want to use the filterbank? (0 = no, 1 = yes): ");
  scanf("%d", &usefilterbank);
  if(usefilterbank) {
    fprintf(stderr, "Filterbank config file (e.g. P501_MB.cfg): ");
    scanf("%s", fbconfig);
    fprintf(stderr, "Filterbank sampling time (e.g. 250 microseconds): ");
    scanf("%d", &fbsamp);
    if(lookuptable[0] != 0) {
      for(i = 0; i < nrlookup_entries; i++)
	lookuptable_entrys[i].fbsamp = fbsamp;
    }
  }


  printf("==================================================================================\nSpecified options:\n");
  /*  gethourstring(lst, txt, 1); */
  converthms_string(txt, lst, 0, 1);
  printf("Start LST: %s\n", txt);
  printf("Multiply durations with factor: %.1lf\n", timemultiplier);
  printf("Time between cals : ");
  if(cal_separation == 0)
    printf("always\n");
  else
    printf("%ld sec\n", cal_separation);
  printf("Time between levcal's: ");
  if(levcal_wait == 0)
    printf("always\n");
  else if(levcal_wait == -1)
    printf("never\n");
  else
    printf("%d sec\n", levcal_wait);
  printf("Number schedules files: %d\n", nrschedules);
  if(nrschedules > 0) {
    for(i = 0; i < nrschedules; i++) {
      printf("For schedule %d:\n", i+1);
      printf("  receiver:       %s\n", recv[i]);
      printf("  config file:    %s\n", config_file[i]);
      if(srchmode[i]) {
	printf("  search mode:    yes\n");
	printf("  time constant:  %.1f s\n", timeconst[i]);
	printf("  tobs_srch_levcal: %d s\n", tobs_srch_levcal);
      }else {
	printf("  search mode:    no\n");
      }
      printf("  guessed dfb nr: %d\n", dfbnr[i]);
      if(dfbnr[i] < 1 || dfbnr[i] > 4) {
	fprintf(stderr, "The dfb nr doesn't appear to make sense\n");
	return 0;
      }
      printf("  frequency:      %.1f MHz\n", freq[i]);
      printf("  levcal's:       ");
      if(schedule_levcal[i])
	printf("yes\n");
      else
	printf("no\n");
    }
  }
  printf("wbtsub:           %d sec\n", wbtsub);
  printf("wbtsub's to skip: %d\n", wbtsub_skip);
  printf("wbtsub_cal:       %d sec\n", wbtsub_cal);
  printf("tobs_cal:         %d sec\n", tobs_cal);
  printf("Feed rotation: ");
  if(feed_rot)
    printf("yes\n");
  else
    printf("no\n");
  printf("Use filterbank: ");
  if(usefilterbank)
    printf("yes\n");
  else
    printf("no\n");
  if(usefilterbank) {
    printf("Filterbank config file  : %s\n", fbconfig);
    printf("Filterbank sampling time: %d microseconds\n", fbsamp);
  }
  printf("==================================================================================\n");
  if(lookuptable[0] != 0 && usefbsamptable) {
    if(!readfbtable(fbsamptable, lookuptable_entrys, nrlookup_entries))
      return 0;
    fprintf(stderr, "Found %d entries in fb lookup table\n", nrlookup_entries);
  }

  if(nrschedules > 0) {
    for(i = 0; i < nrschedules; i++) {
      sprintf(txt, "sched%d", i+1);
      if(verbose) printf("Opening %s\n", txt);
      fout[i] = fopen(txt, "w");
      if(fout[i] == NULL) {
	fprintf(stderr, "Cannot open %s\n", txt);
	return 0;
      }
    }
  }
  if(usefilterbank) {
    if(verbose) printf("Opening sched_f\n");
    foutf = fopen("sched_f", "w");
    if(foutf == NULL) {
      fprintf(stderr, "Cannot open sched_f\n");
      return 0;
    }
  }

  if(gnuplot) {
    if(verbose) printf("Opening sky.gnu\n");
    gnuplotf = fopen("sky.gnu", "w");
    if(gnuplotf == NULL) {
      fprintf(stderr, "Cannot open sky.gnu\n");
      return 0;
    }
    fprintf(gnuplotf, "set term postscript color\n");
    fprintf(gnuplotf, "set output 'sky.ps'\n");
    fprintf(gnuplotf, "set size square\n");
    fprintf(gnuplotf, "set parametric\n");
    fprintf(gnuplotf, "unset key\n");
    fprintf(gnuplotf, "set xrange[-1.2:1.2]\n");
    fprintf(gnuplotf, "set yrange[-1.2:1.2]\n");

    //    fprintf(gnuplotf, "plot cos(t)*sin(90.0*3.14159265359/180.0),sin(t)*sin(90.0*3.14159265359/180.0) lt 1");
    //    fprintf(gnuplotf, ", cos(t)*sin(60.0*3.14159265359/180.0),sin(t)*sin(60.0*3.14159265359/180.0) lt 1");
    //    fprintf(gnuplotf, ", cos(t)*sin(30.0*3.14159265359/180.0),sin(t)*sin(30.0*3.14159265359/180.0) lt 1");


    fprintf(gnuplotf, "set label \"N\" at 0, 1.05 center\n");
    fprintf(gnuplotf, "set label \"S\" at 0, -1.05 center\n");
    fprintf(gnuplotf, "set label \"W\" at -1.05,0 center\n");
    fprintf(gnuplotf, "set label \"E\" at 1.05,0 center\n");

    fprintf(gnuplotf, "set arrow from sin(%e),cos(%e) to sin(%e)*(2./3.),cos(%e)*(2./3.) nohead\n", 205.0*M_PI/180.0, 205.0*M_PI/180.0, 205.0*M_PI/180.0, 205.0*M_PI/180.0);
    fprintf(gnuplotf, "set arrow from sin(%e),cos(%e) to sin(%e)*(2./3.),cos(%e)*(2./3.) nohead\n", 295.0*M_PI/180.0, 295.0*M_PI/180.0, 295.0*M_PI/180.0, 295.0*M_PI/180.0);
    fprintf(gnuplotf, "\n");
    fprintf(gnuplotf, "\n");
    fprintf(gnuplotf, "\n");

    fprintf(gnuplotf, "plot cos(t)*1.0,sin(t)*1.0 lt 1");
    fprintf(gnuplotf, ", cos(t)*(2.0/3.0),sin(t)*(2.0/3.0) lt 1");
    fprintf(gnuplotf, ", cos(t)*(1.0/3.0),sin(t)*(1.0/3.0) lt 1");
  }

  if(!readSourceList(sources, &nrsources, argv[argc-1])) {
    return 0;
  }else {
    linenr = 1;
    for(i = nrentriesskip; i < nrsources; i++) {
      extractRADEC(sources[i].name, &sources[i].raj, &sources[i].decj);
      //      printf("XXXX read declination as %lf deg\n", sources[i].decj);
      if(lookuptable[0] != 0) {
	lookup(&sources[i], lookuptable_entrys, nrlookup_entries);
      }else {
	sources[i].rms = sqrt(-1);
	sources[i].fbsamp = sqrt(-1);
	sources[i].mjd = sqrt(-1);
      }

      sources[i].edot = sqrt(-1);
      sources[i].dist = sqrt(-1);

      if(UseCatalog) {
	if(!catalogRADEC(sources[i].name, RAJ, DECJ)) {
	  return 0;
	}
	converthms(RAJ, &sources[i].raj);
	converthms(DECJ, &sources[i].decj);
	if(!catalogEdotDist(sources[i].name, &sources[i].edot, &sources[i].dist)) {
	  return 0;
	}
      }
      sources[i].raj *= M_PI/12.0;
      sources[i].decj *= M_PI/180.0;
      riseset(sources[i].raj, sources[i].decj, elevation_limit, latitude, &sources[i].rise, &sources[i].set, verbose);
      /*      gethourstring(12.0*sources[i].rise/M_PI, txt, 0); */
      converthms_string(txt, 12.0*sources[i].rise/M_PI, -1, 1);
      sprintf(sources[i].risesetstr, "rise %s - set ", txt);
      /*      gethourstring(12.0*sources[i].set/M_PI, txt, 0); */
      converthms_string(txt, 12.0*sources[i].set/M_PI, -1, 1);
      strcat(sources[i].risesetstr, txt);
      sources[i].lst0 = lst;
      /*      gethourstring(sources[i].lst0, sources[i].lst0str, 0); */
      if(sources[i].lst0 > 24)
	converthms_string(sources[i].lst0str, sources[i].lst0-24, -1, 1);
      else
	converthms_string(sources[i].lst0str, sources[i].lst0, -1, 1);
      calculate_az_zen(sources[i].raj, sources[i].decj, sources[i].lst0*M_PI/12.0, latitude, &sources[i].az, &sources[i].zen);

      /* Calculate drive time */
      if(calc_slew && i != nrentriesskip) {
	double diff;
	diff = fabs(sources[i].az - sources[i-1].az);
	if(fabs(sources[i].az - sources[i-1].az+2.0*M_PI) < diff)
	  diff = fabs(sources[i].az - sources[i-1].az+2.0*M_PI);
	if(fabs(sources[i].az - sources[i-1].az-2.0*M_PI) < diff)
	  diff = fabs(sources[i].az - sources[i-1].az-2.0*M_PI);
	dAz = diff*180.0/M_PI;
	dZen = (sources[i].zen - sources[i-1].zen)*180.0/M_PI;
	tAz = 60.0*dAz/24.0;
	if(dZen > 0)
	  tZen = 60.0*fabs(dZen)/12.0;
	else
	  tZen = 60.0*fabs(dZen)/10.0;
	if(tZen > tAz)
	  dt = tZen;
	else
	  dt = tAz;
	/*	printf("dAz = %f, dZen = %f, tAz = %f, tZen = %f, dt = %f\n", dAz, dZen, tAz, tZen, dt); */
	lst += dt/3600.0;
      }
      lst += slewtime/3600.0;

      /* Update time and position */
      sources[i].lst0 = lst;
      /*      gethourstring(sources[i].lst0, sources[i].lst0str, 0); */
      if(sources[i].lst0 > 24)
	converthms_string(sources[i].lst0str, sources[i].lst0-24, -1, 1);
      else
	converthms_string(sources[i].lst0str, sources[i].lst0, -1, 1);
      calculate_az_zen(sources[i].raj, sources[i].decj, sources[i].lst0*M_PI/12.0, latitude, &sources[i].az, &sources[i].zen);

      if(gnuplot) {
	int itt;
	double curlst, az, zen;
	char txt[100];
	for(itt = 0; itt < 21; itt++) {
	  curlst = sources[i].lst0 + ((double)itt/20.0)*sources[i].duration*timemultiplier/3600.0;
	  calculate_az_zen(sources[i].raj, sources[i].decj, curlst*M_PI/12.0, latitude, &az, &zen);
	  //	  printf("XXXX %f %f\n", zen, zen*180.0/M_PI);
	  if(zen < M_PI/2.0)
	    fprintf(gnuplotf, ",sin(%e)*(%e),cos(%e)*(%e) w p pt 7 lc 3", az, zen/(0.5*M_PI), az, zen/(0.5*M_PI));
	  sprintf(txt, "set label \"%s\" at sin(%e)*(%e),cos(%e)*(%e) left font \"Helvetica,8\"\n", sources[i].name, az, zen/(0.5*M_PI), az, zen/(0.5*M_PI));
	  if(itt == 0)
	    strcat(gnuplot_buffer, txt);
	}
      }


      sources[i].docal = 0;
      sources[i].dolevcal = 0;
      if(i == nrentriesskip) {
	sources[i].docal = 1;
	lst_lastcal = lst;
        if(levcal_wait >= 0) {
	  sources[i].dolevcal = 1;
	  lst_lastlevcal = lst;
        }
      }else if((lst - lst_lastcal)*3600.0 > cal_separation) {
	sources[i].docal = 1;
	lst_lastcal = lst;
        if(levcal_wait >= 0) {
  	  if((lst - lst_lastlevcal)*3600.0 > levcal_wait) {
	    sources[i].dolevcal = 1;
	    lst_lastlevcal = lst;
          }
	}
      }else if(i > nrentriesskip) {
	if(sources[i-1].calaftersource != 0) {
	  lst_lastcal = lst;
	  sources[i].docal = 1;
	  if(sources[i-1].calaftersource == 2) {
	    sources[i].dolevcal = 1;
	    lst_lastlevcal = lst;
	  }
	}
      }
      if(sources[i].zen*180.0/M_PI < 4.5) {
	printf("WARNING: Zenith angle is too small!!\n");
      }
      t1 = M_PI*lst/12.0;
      t2 = lst+sources[i].duration*timemultiplier/3600.0;
      if(sources[i].docal)
	t2 += tobs_cal/3600;
      t2 *= M_PI/12.0;
      if(t1 > 2.0*M_PI)
	t1 -= 2.0*M_PI;
      if(t2 > 2.0*M_PI)
	t2 -= 2.0*M_PI;
      //            printf("XXXXX rise=%lf, set=%lf, start=%lf, end=%lf\n", sources[i].rise, sources[i].set, t1, t2);
      if(sources[i].rise < sources[i].set) {
	if(t1 < t2) {
	  if(t1 < sources[i].rise) {
	    printf("WARNING: Source is below horizon at start of observation!!!\n");
	  }
	  if(t2 > sources[i].set) {
	    printf("WARNING: Source is below horizon at end of observation!!!\n");
	  }
	}else {
	  printf("WARNING: Source is below horizon during observation!!!\n");
	}
      }else {
	if(t1 < t2) {
	  if(!((t1 <= sources[i].set && t2 <= sources[i].set) || (t1 >= sources[i].rise && t2 >= sources[i].rise))) {
	    printf("WARNING: Source is below horizon during observation!!!\n");
	  }
	}else {
	  if(!(t1 >= sources[i].rise && t2 <= sources[i].set)) {
	    printf("WARNING: Source is below horizon during observation!!!\n");
	  }
	}
      }

      /*      gethourstring(12.0*sources[i].raj/M_PI, RAJ, 0); */
      converthms_string(RAJ, 12.0*sources[i].raj/M_PI, -1, 1);

      strcat(RAJ, ":00");
      if(sources[i].decj < 0)
	sprintf(DECJ, "-");
      else
	DECJ[0] = 0;
      sprintf(txt, "%02d:00:00", (int)fabs(180.0*sources[i].decj/M_PI));
      //      printf("XXXXX dec was set to %lf\n", sources[i].decj);
      strcat(DECJ, txt); 
      /* printf("%s - %s: ", RAJ, DECJ); */
      /*if(strcmp(sources[i].name, "J0543+2329") == 0) fprintf(stderr, "\nBlaf: %s %f %f --- '%s' '%s'\n", sources[i].name, sources[i].raj, sources[i].decj, RAJ, DECJ);*/

      if(sources[i].docal) {
	sprintf(sources[i].linenrstr, "%03d-%03d", linenr, linenr+1);
	linenr += 2;
      }else {
	sprintf(sources[i].linenrstr, "%03d-%03d", linenr, linenr);
	linenr += 1;
      }

      if(sources[i].docal == 0) {
	txt[0] = 'X';
      }else {
	txt[0] = 'C';
	if(sources[i].dolevcal)
	  txt[0] = 'L';
      }
      /*      printf("%3d %s %c %-10s%s dur=%4ld (%02ldm) raj=%lf decj = %lf %s A = %lf Z=%lf %s\n", i+1, sources[i].linenrstr, txt[0], sources[i].name, sources[i].lst0str, sources[i].duration, sources[i].duration/60, sources[i].raj, sources[i].decj, sources[i].risesetstr, sources[i].az*180.0/M_PI, sources[i].zen*180.0/M_PI, sources[i].commend);*/
      if(i == nrentriesskip) {
	printf(" NR  PSR NAME   LST    ");
	if(extendedprintout)
	  printf("        ");
	printf("DUR | LINENR CAL (RISE LST   - SET LST  ) | ");
	if(extendedprintout)
	  printf("AZMUT ZENIT  ");
	if(extendedprintout)
	  printf("EDOT         GAMMA-FLUX  ");
	printf("RMS_RESIDS DAYS_SINCE_LAST_OBS\n");
      }
      printf("%3d %-10s %s ", i+1, sources[i].name, sources[i].lst0str);
      if(extendedprintout)
	printf("%5lds = ", sources[i].duration);
      printf("%3ldm | %s %c  (%s) | ", sources[i].duration/60, sources[i].linenrstr, txt[0], sources[i].risesetstr);
      if(extendedprintout)
	printf("A=%3.0lf Z=%4.1lf ", sources[i].az*180.0/M_PI, sources[i].zen*180.0/M_PI);
      if(extendedprintout)
	printf("edot=%.1e flux=%6.1f ", sources[i].edot, 1e-15*sqrt(sources[i].edot)/(sources[i].dist*sources[i].dist));
      printf("rms=%5.1f days=%.0f %-10s ", sources[i].rms, latestmjd-sources[i].mjd, sources[i].name);
      if(strlen(sources[i].commend) > 0)
	printf(" ### %s\n", sources[i].commend);
      else
	printf("\n");


      if(nrschedules > 0) {
	for(j = 0; j < nrschedules; j++) {
	  fprintf(fout[j], "#\n");
	  fprintf(fout[j], "# PSR %s Tobs = %ld LST start %s:00 (%s)\n", sources[i].name, sources[i].duration, sources[i].lst0str, sources[i].risesetstr);
	  fprintf(fout[j], "#                Azimuth = %.2f degrees and Zenith angle = %.2f\n", sources[i].az*180.0/M_PI, sources[i].zen*180.0/M_PI);

	  /*
J0437-4715_R; P 04:37:00 -47:35:00; d3_cfg pdfb3_2048_1024_1024; d3_frq 3100.0; d3_tcyc 20; d3_tsub 60; CPFRQ1 3256.0; CPFRQ2 685.0; rcvr 1050CM; tobs 180; mode levcal; cal freq; fdmode FA; fdang 45.0
J0437-4715; J J0437-4715; mode psr;  tobs 3840

	  */
	  //	  printf("XxXXXx %d %d\n", sources[i].docal, srchmode[i]);
	  if(sources[i].docal) {
	    if(srchmode[j] == 0) {
	      if(sources[i].dolevcal && schedule_levcal[j])
		sprintf(txtcal, "levcal");
	      else
		sprintf(txtcal, "cal");
	      if(nextgen)
		fprintf(fout[j], "%s_R; P %s %s; rcvr %s; d%d_cfg %s; d%d_frq %.1f; d%d_tcyc 10; mode %s; cal freq; d%d_tsub %d; tobs %d\n", sources[i].name, RAJ, DECJ, recv[j], dfbnr[j], config_file[j], dfbnr[j], freq[j], dfbnr[j], txtcal, dfbnr[j], wbtsub_cal, tobs_cal);
	      else
		fprintf(fout[j], "%s_R; P %s %s; rcvr %s; wbcfg %s; wbfrq %.1f; mode %s; cal sync; wbtsub %d; tobs %d\n", sources[i].name, RAJ, DECJ, recv[j], config_file[j], freq[j], txtcal, wbtsub_cal, tobs_cal);
	    }else {
	      if(sources[i].dolevcal && schedule_levcal[j]) {
		fprintf(fout[j], "%s_R; P %s %s; rcvr %s; d%d_cfg %s; d%d_frq %.1f; d%d_tcyc 10; mode srchset; sr%d_nsblk 2048; sr%d_ftmx 3600; sr%d_btc %.1f; sr%d_nbit 8; sr%d_nprod 4; sr%d_tsmp 256; tobs %d\n", sources[i].name, RAJ, DECJ, recv[j], dfbnr[j], config_file[j], dfbnr[j], freq[j], dfbnr[j], dfbnr[j], dfbnr[j], dfbnr[j], timeconst[j], dfbnr[j], dfbnr[j], dfbnr[j], tobs_srch_levcal);
	      }
	      fprintf(fout[j], "%s_R; P %s %s; mode srch; tobs %d\n", sources[i].name, RAJ, DECJ, tobs_cal);
	    }
	  }
	  if(!feed_rot) {
	    if(srchmode[j] == 0) {
	      if(nextgen)
		fprintf(fout[j], "%s; J %s; mode psr; d%d_tsub %d; tobs %ld; fdmode FA; fdang 0.0\n", sources[i].name, sources[i].name, dfbnr[j], wbtsub, sources[i].duration);
	      else
		fprintf(fout[j], "%s; J %s; mode psr; cal off; wbtsub %d; tobs %ld; fd_mode FA; fd_ang 0.0\n", sources[i].name, sources[i].name, wbtsub, sources[i].duration);
	    }else {
	      fprintf(fout[j], "%s; J %s; mode srch; tobs %ld; fdmode FA; fdang 0.0\n", sources[i].name, sources[i].name, sources[i].duration);
	    }
	  }else {
	    if(srchmode[j]) {
	      fprintf(stderr, "Feed rotation in search mode is not implemented\n");
	      return 0;
	    }
	    if(nextgen) {
	      fprintf(fout[j], "%s; J %s; mode psr; d%d_tsub %d; tobs %ld; fdmode FA; fdang -45.0\n", sources[i].name, sources[i].name, dfbnr[j], wbtsub, sources[i].duration/2);
	      fprintf(fout[j], "%s; J %s; mode psr; d%d_tsub %d; tobs %ld; fdmode FA; fdang 45.0\n", sources[i].name, sources[i].name, dfbnr[j], wbtsub, sources[i].duration/2);
	    }else {
	      fprintf(fout[j], "%s; J %s; mode psr; cal off; wbtsub %d; tobs %ld; fd_mode FA; fd_ang -45.0\n", sources[i].name, sources[i].name, wbtsub, sources[i].duration/2);
	      fprintf(fout[j], "%s; J %s; mode psr; cal off; wbtsub %d; tobs %ld; fd_mode FA; fd_ang 45.0\n", sources[i].name, sources[i].name, wbtsub, sources[i].duration/2);
	    }
	  }
	}
	if(usefilterbank) {
	  if(sources[i].fbsamp)
	    fb_dt_to_use = sources[i].fbsamp;
	  else
	    fb_dt_to_use = fbsamp;
	  fprintf(foutf, "%s %s t%ld s%d\n", sources[i].name, fbconfig, sources[i].duration -5, fb_dt_to_use);
	}
      }
      lst += sources[i].duration*timemultiplier/3600.0;
      if(sources[i].docal)
	lst += tobs_cal/3600.0;

      lst += sources[i].waitaftersource/3600.0;
      if(sources[i].waitaftersource) {
	printf("############# Waiting %ld seconds.\n", sources[i].waitaftersource);
      }
    }
  }
  lst -= slewtime/3600.0;
  printf("############# Schedule ends at: ");
  printhour(M_PI*lst/12.0);
  printf("\n");
  
  if(nrschedules > 0) {
    for(i = 0; i < nrschedules; i++) {
      fclose(fout[i]);
    }
  }
  if(usefilterbank) {
    fclose(foutf);
  }
  if(gnuplot) {
    fprintf(gnuplotf, "\n"); 
    fprintf(gnuplotf, "%s", gnuplot_buffer);
   fprintf(gnuplotf, "replot\n");
   fprintf(gnuplotf, "set output\n");
    //    fprintf(gnuplotf, "pause -1\n");
    fclose(gnuplotf);
  }

  return 0;
}



void extractRADEC(char *psrname, double *alpha, double *delta)
{
  double u1;
  char psr[200];
  //  printf("XXXXX psr=%s\n", psrname);
  strcpy(psr, psrname);
  if(psr[0] == 'B' || psr[0] == 'J')
    strcpy(psr, &psrname[1]);
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
  (*alpha) /= 60;
  psr[2] = 0;
  sscanf(psr, "%lf", &u1);
  *alpha += u1;
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

int catalogEdotDist(char *psrname, float *edot, float *dist)
{
  char command[200], junk[100];
  FILE *fin;
  int i;

  sprintf(command, "psrcat -nonumber -all -c \"EDOT DIST\" %s", psrname);
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
  for(i = 0; i < 6; i++) {
    fscanf(fin, "%s", junk);
  }
  fscanf(fin, "%e", edot);
  fscanf(fin, "%f", dist);
  fclose(fin);
  
  if(verbose) {
    printf("EDOT:                  %e\n", *edot);
    printf("DIST:                 %f\n", *dist);
  }
  system("rm -f junk");
  return 1;
}




int readSourceList(source_type *sources, int *nrsources, char *filename)
{
  FILE *fin;
  int i, ret, commendmode;
  char txt[200];

  if(verbose) printf("Opening %s\n", filename);
  fin = fopen(filename, "r");
  if(fin == NULL) {
    fprintf(stderr, "Cannot open %s\n", filename);
    return 0;
  }
  *nrsources = 0;
  do {
    sources[*nrsources].waitaftersource = 0;
    ret = fscanf(fin, "%s", txt);
    if(ret != 1) 
      break;
    if(strcmp(txt, "end") == 0 || strcmp(txt, "END") == 0) {
      ret = EOF;
      break;
    }
    if(txt[0] == '#') {
      if(verbose) printf("Skipping line starting with %s\n", txt);
      fseek(fin, -1, SEEK_CUR);
      do {
	ret = fgetc(fin);
	if(ret == EOF)
	  break;
      }while(ret != '\n' && ret != 0 && ret != '\r');
    }else if(strcmp(txt, "wait") == 0 || strcmp(txt, "WAIT") == 0) {
      ret = fscanf(fin, "%ld", &sources[*nrsources - 1].waitaftersource);
      if(ret != 1) 
	break;
    }else if(strcmp(txt, "cal") == 0 || strcmp(txt, "CAL") == 0) {
      /* blaat */
      sources[*nrsources - 1].calaftersource = 1;
    }else if(strcmp(txt, "levcal") == 0 || strcmp(txt, "LEVCAL") == 0) {
      sources[*nrsources - 1].calaftersource = 2;
    }else {
      strcpy(sources[*nrsources].name, txt);
      ret = fscanf(fin, "%ld", &sources[*nrsources].duration);
      if(ret != 1) {
	printf("Error reading duration of source %s\n", sources[*nrsources].name);
	break;
      }
      sources[*nrsources].commend[0] = 0;
      fseek(fin, -1, SEEK_CUR);
      commendmode = 0;
      i = 0;
      do {
	ret = fgetc(fin);
	if(ret == EOF)
	  break;
	else if(ret == '#') {
	  if(commendmode == 0)
	    commendmode = 1;
	  else
	    commendmode = 0;
	}else if(commendmode && ret != '\n' && ret != 0 && ret != '\r') {
          sources[*nrsources].commend[i++] = ret;
          sources[*nrsources].commend[i] = 0;
	}
      }while(ret != '\n' && ret != 0 && ret != '\r');
      *nrsources += 1;
    }
  }while(*nrsources < MaxNrSources);
  fclose(fin);
  if(*nrsources == MaxNrSources) {
    printf("Too many sources\n");
    return 0;
  }
  if(ret != EOF) {
    printf("Error reading file\n");
    return 0;
  }
  return 1;
}

int readlookuptable(char *filename, lookuptable_entry *entries, int *nrlookup_entries, float *latestmjd)
{
  int i, j;
  FILE *fin;

  *latestmjd = -1;
  if((fin = fopen(filename, "r")) == NULL) {
    printf("Error: cannot open %s\n", filename);
    return 0;
  }
  *nrlookup_entries = 0;
  for(i = 0; i < MaxNrLookupEntries; i++) {
    j = fscanf(fin, "%s %lf %lf %f %f %f", entries[i].name, &(entries[i].raj), &(entries[i].decj), &(entries[i].edot), &(entries[i].dist), &(entries[i].mjd));
    entries[i].rms = 0;
    entries[i].fbsamp = 0;
    if(j == 6) {
      (*nrlookup_entries)++;
      if(entries[i].mjd > *latestmjd)
	*latestmjd = entries[i].mjd;
      /*      printf("Added %s\n", entries[i].name);  */
    }else
      break;
  }
  fclose(fin);
  return 1;
}

void lookup(source_type *source, lookuptable_entry *entries, int nrlookup_entries)
{
  int i, found;
  found = 0;
  for(i = 0; i < nrlookup_entries; i++) {
    if(strcmp(source->name, entries[i].name) == 0) {
      source->raj = entries[i].raj;
      source->decj = entries[i].decj;
      source->edot = entries[i].edot;
      source->dist = entries[i].dist;
      source->rms = entries[i].rms;
      source->fbsamp = entries[i].fbsamp;
      source->mjd = entries[i].mjd;
      found = 1;
      break;
    }
  }
  if(found == 0)
    printf("Cannot find %s in lookuptable.\n", source->name);
}

int readrmstable(char *filename, lookuptable_entry *entries, int nrlookup_entries)
{
  int i, j, k, found;
  FILE *fin;
  char psrname[100];
  float rms;

  if((fin = fopen(filename, "r")) == NULL) {
    printf("Error: cannot open %s\n", filename);
    return 0;
  }
  for(i = 0; i < MaxNrLookupEntries; i++) {
    j = fscanf(fin, "%s %f", psrname, &rms);
    if(j == 2) {
      found = 0;
      for(k = 0; k < nrlookup_entries; k++) {
	if(strcmp(psrname, entries[k].name) == 0) {
	  entries[k].rms = rms;
	  found = 1;
	  break;
	}
      }
      if(found == 0)
	fprintf(stderr, "PSR %s of the rms list is not in the lookuptable, so this line is ignored.\n", psrname);
    }else
      break;
  }
  fclose(fin);
  return 1;
}

int readfbtable(char *filename, lookuptable_entry *entries, int nrlookup_entries)
{
  int i, j, k, found, fbsamp;
  FILE *fin;
  char psrname[100];

  if((fin = fopen(filename, "r")) == NULL) {
    printf("Error: cannot open %s\n", filename);
    return 0;
  }
  for(i = 0; i < MaxNrLookupEntries; i++) {
    j = fscanf(fin, "%s s%d", psrname, &fbsamp);
    if(j == 2) {
      found = 0;
      for(k = 0; k < nrlookup_entries; k++) {
	if(strcmp(psrname, entries[k].name) == 0) {
	  entries[k].fbsamp = fbsamp;
	  found = 1;
	  break;
	}
      }
      if(found == 0)
	fprintf(stderr, "PSR %s of the filterbank lookuptable is not in the lookuptable, so this line is ignored.\n", psrname);
    }else
      break;
  }
  fclose(fin);
  return 1;
}



int FindLastObservation(char *psrname, float *mjd)
{
  char command[300], junk[100];
  FILE *fin;

  sprintf(command, "vap -c 'mjd' /pulsar/archive17/wel161/glast/pulsars/%s/*.SFTC | sort -n -k 2 | tail -1 > junk", psrname);
  if(verbose)
    printf("%s\n", command);
  system("rm -f junk");
  system(command);
  /*    if(verbose)  system("more junk"); */
  if((fin = fopen("junk", "r")) == NULL) {
    printf("Error: no file 'junk'\n");
    return 0;
  }
  fscanf(fin, "%s", junk);
  /*  printf("strlen=%d %s %c\n", strlen(junk), junk, junk[8]); */
  if(strlen(junk) == 19 && junk[7]=='_') {
    fscanf(fin, "%f", mjd);
  }else {
    *mjd = -1;
  }
  fclose(fin);
  
  if(verbose) {
    printf("MJD:                  %f\n", *mjd);
  }
  return 1;
}


/*
float polar_alfa(float x, float y)
{
  float alfa;
  if(x == 0) {
    if(y == 0)
      return 1.5*M_PI;
    else
      return 0.5*M_PI;
  }else {
    alfa = atan(y/x);
    if(x > 0 && y >= 0)
      return alfa;
    if(x > 0)
      return alfa + 2*M_PI;
    return alfa + M_PI;
  }
}

void calculateAZ(double raj, double decj, double lst, double *az, double *zen)
{
  double x, y;
  *zen = acos(sin(latitude)*sin(decj)+cos(latitude)*cos(decj)*cos(lst-raj));
  x = -sin(lst-raj);
  y = cos(latitude)*tan(decj)-sin(latitude)*cos(lst-raj);
  *az    = polar_angle_rad(y,x);
}
*/



