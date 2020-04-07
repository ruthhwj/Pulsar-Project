#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "psrsalsa.h"

void SHOWREVISIONINFO_prog() {
#include "pmod.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);

int main(int argc, char **argv)
{
  int skiplines, nrColumns, colnumX, colnumY, smooth, ret, firsttime, output_c, index;
  int notascending, notdescending;
  long i, j, nrdatapoints, order, smoothRND;
  float *dataX, *dataY, x, xstart, xend, dx, y, smoothRPLC, mindata, maxdata;
  FILE *fout;
  splineFit_def splineFit;
  psrsalsaApplication application;

  initApplication(&application, "interpol", "[options] file XVALUE (or XSTART,XEND,DX)");
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_output = 1;
  sprintf(application.outputname, "stdout");

  skiplines = 0;
  nrColumns = 2;
  colnumX = 1;
  colnumY = 2;
  order = 2;
  smooth = 0;
  output_c = 0;
  smoothRND = 0;
  smoothRPLC = 1;  // Default: replace all grid points in each itteration
  fout = stdout; 

  if(argc < 3) {
    fprintf(stderr, "This program is designed to read a 2 column ascii file (x f(x)), and outputs\nthe interpolated value f(XVALUE). The x coordinate is assumed to be monotonic.\n\n");
    printApplicationHelp(&application);
    printf("Other options:\n\n");
    fprintf(stderr, "-format     Specify file format: \"NrLinesToSkip NrColumns ColumnNrX ColumnNrY\".\n");
    fprintf(stderr, "            default is \"0 2 1 2\".\n");
    fprintf(stderr, "-n          Use a n-point interpolation function, default is %ld.\n", order);
    fprintf(stderr, "-smooth     Instead if interpolating, use a fit of cubic B-spline functions\n");
    fprintf(stderr, "            with the specified nr of coefficients.\n");
    fprintf(stderr, "-smoothRND  Instead of only using a uniform gridding with the -smooth option,\n");
    fprintf(stderr, "            assign the node position at random. Repeat this process the\n");
    fprintf(stderr, "            given nr of times to determine the optimum.\n");
    fprintf(stderr, "-smoothRPLC Only replace this fraction of grid-points in each itteration of the\n");
    fprintf(stderr, "            -smoothRND option. Default is 1.\n");
    fprintf(stderr, "-smooth_c   Output b-spline as C code to the specified filename (can be set to\n");
    fprintf(stderr, "            stdout). Use together with the -smooth options.\n");
    return 0;
  }else if(argc > 3) {
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcmp(argv[i], "-smooth_c") == 0) {
	output_c = i+1;
	i++;
      }else if(strcmp(argv[i], "-n") == 0) {
	j = sscanf(argv[i+1], "%ld", &order);
	if(j != 1) {
	  printerror(application.verbose_state.debug, "Cannot parse -n option, need 1 value.");
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-smooth") == 0) {
	j = sscanf(argv[i+1], "%ld", &order);
	if(j != 1) {
	  printerror(application.verbose_state.debug, "Cannot parse option %s, need 1 value.", argv[i+1]);
	  return 0;
	}
	smooth = 1;
	i++;
      }else if(strcmp(argv[i], "-smoothRND") == 0) {
	j = sscanf(argv[i+1], "%ld", &smoothRND);
	if(j != 1) {
	  printerror(application.verbose_state.debug, "Cannot parse option %s, need 1 value.", argv[i+1]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-smoothRPLC") == 0) {
	j = sscanf(argv[i+1], "%f", &smoothRPLC);
	if(j != 1) {
	  printerror(application.verbose_state.debug, "Cannot parse option %s, need 1 value.", argv[i+1]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-format") == 0) {
	j = sscanf(argv[i+1], "%d %d %d %d", &skiplines, &nrColumns, &colnumX, &colnumY);
	if(j != 4) {
	  printerror(application.verbose_state.debug, "Cannot parse -format option, need 4 value.");
	  return 0;
	}
	i++;
      }else {
	/* If the option is not recognized, assume it is a filename */
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "ERROR interpol: Unknown option: %s", argv[i]);
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0)
	    return 0;
	}
      }
    }
  }

  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) < 2) {
    printerror(application.verbose_state.debug, "ERROR interpol: A file and a value (or range) need to be specified on the command line.");
    return 0;
  }

  // last command line option should be the value/range command
  ret = sscanf(argv[argc-1], "%f,%f,%f", &xstart, &xend, &dx);
  if(ret == 1) {
    xend = xstart;
    dx = 1;
  }else if(ret != 3) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR interpol: Cannot parse XVALUE (or XSTART,XEND,DX) parameter on command line. This needs to be specified on command line. See usage help.");
    return 0;
  }

  //  if(read_ascii_file_float_column(argv[argc-2], skiplines, nrColumns, colnumX, &dataX, &nrdatapoints, NULL, NULL, NULL, 0, verbose, 1) == 0)
  //    return 0;
  //  if(read_ascii_file_float_column(argv[argc-2], skiplines, nrColumns, colnumY, &dataY, &nrdatapoints, NULL, NULL, NULL, 0, verbose, 1) == 0)
  //    return 0;

  if(read_ascii_column(argv[argc-2], skiplines, '#', nrColumns, 0, &nrdatapoints, colnumX, 1.0, 0, &dataX, &mindata, &maxdata, NULL, application.verbose_state, 1) == 0)
    return 0;
  if(read_ascii_column(argv[argc-2], skiplines, '#', nrColumns, 0, &nrdatapoints, colnumY, 1.0, 0, &dataY, NULL, NULL, NULL, application.verbose_state, 1) == 0)
    return 0;

  if(xstart < mindata || xend > maxdata) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR interpol: Specified interpolation range (from %f to %f) exceeds the range in the data (from %f to %f).", xstart, xend, mindata, maxdata);
    return 0;
  }
  if(nrdatapoints < 2) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR interpol: Need at least two data points.");
    return 0;
  }
  notascending = 0;
  notdescending = 0;
  for(i = 1; i < nrdatapoints; i++) {
    if(dataX[i] < dataX[i-1]) {
      notascending = 1;
    }
    if(dataX[i] > dataX[i-1]) {
      notdescending = 1;
    }
  }
  if(notascending && notdescending) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR interpol: Input data does not appear to be sorted. Sort input first.");
    // Could use gsl_sort_index
    return 0;
  }

  if(strcmp(application.outputname, "stdout") != 0) {
    fout = fopen(application.outputname, "w");
    if(fout == NULL) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "ERROR interpol: Cannot open output file %s.", application.outputname);
      return 0;
    }
  }

  firsttime = 1;
  for(x = xstart; x <= xend; x += dx) {
    if(x <= xend) {
      if(smooth == 0) {
	if(interpolate(dataX, dataY, nrdatapoints, x, &y, order, application.verbose_state) == 0)
	  return 0;
      }else {
	if(firsttime) {
	  if(cubicBspline_fit(dataX, dataY, NULL, nrdatapoints, order, x, &y, NULL, &splineFit, smoothRND, smoothRPLC, application.verbose_state) == 0)
	    return 0;
	  if(output_c) {
	    FILE *fout_code;
	    fprintf(stderr, "interpol: Output C code to file %s\n", argv[output_c]);
	    if(strcasecmp(argv[output_c], "stdout") == 0) {
	      fout_code = stdout;
	    }else {
	      fout_code = fopen(argv[output_c], "w");
	      if(fout_code == NULL) {
		printerror(application.verbose_state.debug, "ERROR interpol: Cannot open file %s", argv[output_c]);
		return 0;
	      }
	    }
	    // Insert C-code
	    fprintf(fout_code, "// Change function name AND error message\n");
	    fprintf(fout_code, "// Returns 1 on success, 0 on error\n");
	    fprintf(fout_code, "int define_spline(splineFit_def *splineFit)\n{\n");
	    fprintf(fout_code, "  splineFit->order = %ld;\n", splineFit.order);
	    fprintf(fout_code, "  splineFit->nbreak = %ld;\n", splineFit.nbreak);
	    fprintf(fout_code, "  splineFit->x_start = %e;\n", splineFit.x_start);
	    fprintf(fout_code, "  splineFit->x_end = %e;\n", splineFit.x_end);
	    fprintf(fout_code, "  splineFit->coefficients = malloc(%ld*sizeof(double));\n", splineFit.nbreak + splineFit.order - 2);
	    fprintf(fout_code, "  splineFit->breakpoints = malloc(%ld*sizeof(double));\n", splineFit.nbreak);
	    fprintf(fout_code, "  if(splineFit->coefficients == NULL || splineFit->breakpoints == NULL) {\n");
	    fprintf(fout_code, "    fprintf(stderr, \"ERROR define_spline: Cannot allocate memory\\n\");\n");
	    fprintf(fout_code, "    return 0;\n");
	    fprintf(fout_code, "  }\n");
	    fprintf(fout_code, " ");
	    for(i = 0; i < splineFit.nbreak + splineFit.order - 2; i++) {
	      fprintf(fout_code, " splineFit->coefficients[%ld] = %e;", i, splineFit.coefficients[i]);
	    }
	    fprintf(fout_code, "\n");
	    for(i = 0; i < splineFit.nbreak; i++) {
	      fprintf(fout_code, " splineFit->breakpoints[%ld] = %e;", i, splineFit.breakpoints[i]);
	    }
	    fprintf(fout_code, "\n");
	    fprintf(fout_code, "  return 1;\n");
	    fprintf(fout_code, "}\n");
	    if(strcasecmp(argv[output_c], "stdout"))
	      fclose(fout_code);
	  }
	  firsttime = 0;
	}else {
	  //	  fprintf(stderr, "Unsing old solution!\n");
	  if(cubicBspline_eval(splineFit, x, &y, application.verbose_state) == 0)
	    return 0;
	}
      }
      
      fprintf(fout, "x = %e = %f  ,  y = %e = %f\n", x, x, y, y);
    }
  }

  if(strcmp(application.outputname, "stdout") == 0) {
    fclose(fout);
  }

  free(dataX);
  free(dataY);

  terminateApplication(&application);
  return 0;
}


 
