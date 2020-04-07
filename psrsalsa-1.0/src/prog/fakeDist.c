//START REGION RELEASE
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "psrsalsa.h"

//START REGION DEVELOP
void SHOWREVISIONINFO_prog() {
#include "fakeDist.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);
//START REGION RELEASE

int main(int argc, char **argv)
{
  psrsalsaApplication application;
  int nrfunctions, nrfunctions2, outputfile, AddNulls, cmd_line_end_first_distr;
  int noisefile_id, randomize_seed, quiet;
  long NumberPoints, NumberPoints2, i, j, loopnr, nrloops, pointnr, n_noisedata, idnum;
  double noisesigma, average_value, *data_noise;
  long double total;
  /*
  int bin, functionnr, penergyflag;
  float Powerlaw_a,Powerlaw_b,Powerlaw_c, sin_a, totalE, average_value, di;
  float emin, emax, *distribution, *cdistribution, *chance, *chance2, E, x;
  long l, i_start, i_end, i_start2, i_end2, idnum, idnum2, nrebins, nrNulls;
  FILE *fout;
  char tmpstr[10000];
  */

//START REGION RELEASE
  initApplication(&application, "fakeDist", "[options]");
//START REGION DEVELOP
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
//START REGION RELEASE

  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_fixseed = 1;

  NumberPoints = 10000;
  NumberPoints2 = 0;   // Points in second distribution
  nrfunctions = 0;   
  nrfunctions2 = 0;    // Nr of functions after -N2
  outputfile = 0;
  noisesigma = 0;
  AddNulls = 0;
  average_value = 1;
  noisefile_id = 0;
  nrloops = 1;
  idnum = 1;
  randomize_seed = 1;
  cmd_line_end_first_distr = 0;
  quiet = 0;
  /*
  emin = -10;
  emax = 10;
  nrebins = 1000;
  idnum = 0;
  idnum2 = 0;
  functionnr = 0;
  penergyflag = 0;
  verbose = 0;
  */

  if(argc < 2) {
    printf("Program to generate a random list of values drawn from various (combinations of)\n");
    printf("distribution functions. When multiple distributions are specified, they are\n");
    printf("convolved. Functions specied after the -N2 (convoved with each other) are\n");
    printf("combined with the previous distribution. Examples:\n");
    printf("\n");
    printf("To generate a lognormal distribution convolved with a normal distribution:\n");
    printf("  fakeDist -N 100 -lognorm \"1.0 0.2\" -norm \"0 0.5\"\n");
    printf("To generate the same 100 values, followed by 1000 values drawn from a powerlaw\n");
    printf("distribution convolved with a normal distribution:\n");
    printf("  fakeDist -N 100 -lognorm \"1.0 0.2\" -norm \"0 0.5\" -N2 1000 -pwrlaw \"-2 0.5\" -norm \"0 0.5\"\n");
    printf("\n");
    printApplicationHelp(&application);
    printf("Distribution functions:\n\n");
    printf("-flat      \"min max\", specify a flat/uniform distribution between min and max.\n");
    printf("-gamma     \"k theta\", specify gamma function\n");
    printf("-null av   add zero's to distribution to make average equal to av. This implies\n");
    printf("           the nr of generated values is larger than what is specied with -N.\n");
    printf("-norm      \"mu sigma\", specify normal (Gaussian) function  with mean mu and\n");
    printf("           standard devitation sigma. dmu and dsigma\n");
    printf("-lognorm   \"mu sigma\", specify lognormal function\n");
    printf("-pwrlaw    \"idx min\", specify function f(x)=x**idx for x>=min, where idx should\n");
    printf("           be a negative number.\n");
    printf("-Rayleigh  \"sigma\", specify Rayleigh function.\n");
//START REGION DEVELOP
    printf("-rndwalk   \"stepsize n\", each generated value is obtained after doing n steps\n");
    printf("           with the given stepsize.\n");
//START REGION RELEASE
    printf("-sin       \"a b min max\", specify |sin(ax+b)| function, where x is in degrees\n");
    printf("           and between min and max.\n");
    printf("\nOther options:\n\n");
    printf("-N nr      Generate nr points, default is %ld.\n", NumberPoints);
    printf("-N2 nr     Generate nr points in second distribution (example given at top).\n");
    printf("-sigma s   Convolve the generated distribution with a Gaussian with sigma s.\n");
    printf("           This is identical to using -norm \"0 s\".\n");
    printf("-noisefile Use specified file (list of values) instead of Gaussian noise.\n");
    printf("           The values are read in from the first column.\n");
    printf("-output    Output to this file instead to stdout.\n");
    printf("-loop N    Run program N times (separate output files are generated).\n");
    printf("-seed S    Initialise the random number seed with this value, rather than\n");
    printf("           assigning it 'randomly' based on the clock time. -fixseed is\n");
    printf("           identical to -seed 1.\n");
    printf("\nPlease use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about fitting distributions (in the context of pulse energies) can be found in:\n");
    printf(" - Weltevrede et al. 2006, A&A, 458, 269\n\n");
    printCitationInfo();
  /*
    printf("Usage: fake_distribution\n\nOptional options:\n");
    printf("-E    Output in penergy style\n");
    printf("-r    \"min max points\", specify energy binning (def. is \"%.0f %.0f %ld\")\n", emin, emax, nrebins);
    printf("-1    add nulls to get average energy 1.\n");
    printf("-v    Turn on verbose.\n");
  */
    terminateApplication(&application);
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      int index;
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcasecmp(argv[i], "-Rayleigh") == 0) {
	float dummy_float;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &dummy_float, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakeDist: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	if(NumberPoints2 == 0)
	  nrfunctions++;
	else
	  nrfunctions2++;
        i++;
      }else if(strcmp(argv[i], "-gamma") == 0 || strcmp(argv[i], "-flat") == 0 || strcmp(argv[i], "-norm") == 0 || strcmp(argv[i], "-lognorm") == 0 || strcmp(argv[i], "-pwrlaw") == 0
//START REGION DEVELOP
	       || strcmp(argv[i], "-rndwalk") == 0
//START REGION RELEASE
	       ) {
	float dummy_float;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &dummy_float, &dummy_float, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakeDist: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	if(NumberPoints2 == 0)
	  nrfunctions++;
	else
	  nrfunctions2++;
        i++;
      }else if(strcmp(argv[i], "-sin") == 0) {
	float dummy_float;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f %f %f", &dummy_float, &dummy_float, &dummy_float, &dummy_float, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakeDist: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	if(NumberPoints2 == 0)
	  nrfunctions++;
	else
	  nrfunctions2++;
        i++;
      }else if(strcmp(argv[i], "-sigma") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &noisesigma, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakeDist: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-N") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &NumberPoints, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakeDist: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-N2") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &NumberPoints2, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakeDist: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	cmd_line_end_first_distr = i;
        i++;
      }else if(strcmp(argv[i], "-output") == 0) {
	outputfile = i+1;
	i++;
      }else if(strcmp(argv[i], "-null") == 0) {
	AddNulls = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &average_value, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakeDist: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-noisefile") == 0) {
	noisefile_id = i+1;
        i++;
      }else if(strcmp(argv[i], "-loop") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &nrloops, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakeDist: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-seed") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &idnum, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakeDist: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	randomize_seed = 0;
        i++;
      }else if(strcmp(argv[i], "-quiet") == 0) {  // Option to reduce the output, which is used by pdistFit
	quiet = 1;
      }else {
	printerror(application.verbose_state.debug, "fakeDist: Unknown option: %s\n\nRun fakeDist without command line arguments to show help", argv[i]);
	terminateApplication(&application);
	return 0;
      }
    }
    /*
      }else if(strcmp(argv[i], "-1") == 0) {
	AddNulls = 1;
	average_value = 1;
      }else if(strcmp(argv[i], "-A") == 0) {
	AddNulls = 1;
	average_value = atof(argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-E") == 0) {
	penergyflag = 1;
    */
  }

  if(nrfunctions == 0) {
    printerror(application.verbose_state.debug, "fakeDist: Specify at least one distribution function on the command line.\n");
    return 0;
  }
  if(NumberPoints2 > 0 && nrfunctions2 == 0) {
    printerror(application.verbose_state.debug, "fakeDist: Specify at least one distribution function after the -N2 command line option.\n");
    return 0;
  }

  if(noisefile_id != 0) {
    if(application.verbose_state.verbose)
      fprintf(stdout, "Loading noise values from ascii file\n");
    if(read_ascii_column_double(argv[noisefile_id], 0, '#', -1, 1, &n_noisedata, 1, 1.0, 0, &data_noise, NULL, NULL, NULL, application.verbose_state, 1) == 0) {
      printerror(application.verbose_state.debug, "ERROR fakeDist: Cannot load file with noise values.\n");
      return 0;
    }
    //    Readfile(argv[noisefile_id], &n_noisedata, &data_noise);
  }

  gsl_rng *rand_num_gen;
  const gsl_rng_type *rand_num_gen_type;
  gsl_rng_env_setup();    /* Set the default generators, can be influenced by environment variables */
  rand_num_gen_type = gsl_rng_default;
  rand_num_gen = gsl_rng_alloc(rand_num_gen_type);
  if(application.fixseed)
    idnum = 1;
  else if(randomize_seed)
    randomize_idnum(&idnum);
  gsl_rng_set(rand_num_gen, idnum);

  int distr_number;

  for(loopnr = 0; loopnr < nrloops; loopnr++) {
    FILE *fout;
    if(!outputfile) {
      fout = stdout;
    }else {
      char *tmpstr;
      tmpstr = malloc(strlen(argv[outputfile])+9+2);
      if(tmpstr == NULL) {
	printerror(application.verbose_state.debug, "ERROR fakeDist: fakeDist: Cannot allocate memory.\n");
	return 0;
      }
      if(nrloops == 1) {
	strcpy(tmpstr, argv[outputfile]);
      }else {
	sprintf(tmpstr, "%09ld%s", loopnr, argv[outputfile]);
      }
      fout = fopen(tmpstr, "w");
      if(fout == NULL) {
	printerror(application.verbose_state.debug, "ERROR fakeDist: Cannot open '%s'\n", tmpstr);
	return 0;
      }
      if(application.verbose_state.verbose) {
	fprintf(stderr, "file %s opened for output.\n", tmpstr);
      }
      free(tmpstr);
    }

    total = 0;
    // Loop over the two (possible) input distributions.
    for(distr_number = 0; distr_number < 2; distr_number++) {
      double sample;
      long max_pointnr;
      // Loop over the number of points that need to be generated.
      if(distr_number == 0)
	max_pointnr = NumberPoints;
      else
	max_pointnr = NumberPoints2;
      for(pointnr = 0; pointnr < max_pointnr; pointnr++) {
	sample = 0; // Add the different distributions as specified on command line
	int cmd_line_start;
	cmd_line_start = 1;
	if(distr_number == 1) {  // Skip the first distribution functions
	  cmd_line_start = cmd_line_end_first_distr;
	}
	for(i = cmd_line_start; i < argc; i++) { // Note: sscanf's are fine here as reading in already tested above
	  if(strcmp(argv[i], "-norm") == 0) {
	    double mu, sigma;
	    j = sscanf(argv[i+1], "%lf %lf", &mu, &sigma);
	    i++;
	    if(pointnr == 0 && application.verbose_state.verbose) {
	      fprintf(stderr, "Using distribution: exp(-(x-%f)^2/(2*%f^2))/(sqrt(2*pi)*%f)\n", mu, sigma, sigma);
	    }
	    sample += mu + gsl_ran_gaussian(rand_num_gen, sigma); 
	  }else if(strcmp(argv[i], "-lognorm") == 0) {
	    double mu, sigma;
	    j = sscanf(argv[i+1], "%lf %lf", &mu, &sigma);
	    i++;
	    if(pointnr == 0 && application.verbose_state.verbose) {
	      fprintf(stderr, "Using distribution: exp(-(log(x)-%f)^2/(2*%f^2))/(%f*x*sqrt(2*pi))\n", mu, sigma, sigma);
	    }
	    sample += gsl_ran_lognormal(rand_num_gen, mu, sigma); 
	    //	sample += exp(mu+gsl_ran_gaussian(rand_num_gen, sigma));
	  }else if(strcmp(argv[i], "-pwrlaw") == 0) {
	    // The Pareto Distribution is like powerlaw (257)
	    // (A/B)/(x/B)^(A+1), for x >= B
	    // (A/B)*(x/B)^(-A-1), for x >= B
	    // A*B^(-1)*x^(-A-1)*B^(A+1), for x >= B
	    // A*B^A*x^(-A-1), for x >= B
	    //
	    // Specified powerlaw is of form: x^a for x>=b
	    //    => B = b
	    //    => -A-1 = a, so 
	    //       A = -a -1
	    // Should be straight line distribution if plotting N vs log x
	    double a, b, pareto_A, pareto_B;
	    j = sscanf(argv[i+1], "%lf %lf", &a, &b);
	    i++;
	    pareto_A = -a -1.0;
	    pareto_B = b;
	    if(pointnr == 0 && application.verbose_state.verbose) {
	      fprintf(stderr, "Using distribution: %f*x^%f for x >= %f\n", pareto_A*pow(pareto_B, pareto_A), a, b);
	    }
	    sample += gsl_ran_pareto(rand_num_gen, pareto_A, pareto_B); 
	  }else if(strcmp(argv[i], "-flat") == 0) {
	    double min, max;
	    j = sscanf(argv[i+1], "%lf %lf", &min, &max);
	    i++;
	    if(pointnr == 0 && application.verbose_state.verbose) {
	      fprintf(stderr, "Using distribution: 1/(%f-%f)\n", max, min);
	    }
	    sample += gsl_ran_flat(rand_num_gen, min, max); 
	  }else if(strcmp(argv[i], "-gamma") == 0) {
	    double k, theta;
	    j = sscanf(argv[i+1], "%lf %lf", &k, &theta);
	    i++;
	    if(pointnr == 0 && application.verbose_state.verbose) {
	      fprintf(stderr, "Using distribution: %f*x^%f*exp(-x/%f)\n", 1.0/(tgamma(k)*pow(theta, k)), k-1.0, theta);
	    }
	    sample += gsl_ran_gamma(rand_num_gen, k, theta); 
//START REGION DEVELOP
	  }else if(strcmp(argv[i], "-rndwalk") == 0) {
	    double stepsize;
	    long int nrsteps, sum, itt;
	    j = sscanf(argv[i+1], "%lf %ld", &stepsize, &nrsteps);
	    i++;
	    if(pointnr == 0 && application.verbose_state.verbose) {
	      fprintf(stderr, "Using random walk with stepsize %lf and %ld steps\n", stepsize, nrsteps);
	    }
	    sum = 0;
	    for(itt = 0; itt < nrsteps; itt++) {
	      if(gsl_rng_uniform_int(rand_num_gen, 2) == 0)
		sum += 1;
	      else
		sum -= 1;
	      if(sum > 2147483640 || sum < -2147483640) {
		printerror(application.verbose_state.debug, "Too many steps in random walk: overflow happened.\n");
		return 0;
	      }
	    }
	    sample += sum*stepsize;
//START REGION RELEASE
	  }else if(strcmp(argv[i], "-sin") == 0) {
	    double a, b, min, max, angle;
	    j = sscanf(argv[i+1], "%lf %lf %lf %lf", &a, &b, &min, &max);
	    i++;
	    if(pointnr == 0 && application.verbose_state.verbose) {
	      fprintf(stderr, "Using distribution: |sin(%f*x+%f)| with %f <= x <= %f\n", a, b, min, max);
	    }
	    do {
	      // Sinusoidal distribution between 0 and 180.0
	      angle = acos((2.0*gsl_rng_uniform(rand_num_gen)-1.0))*180.0/M_PI;
	      // Sinusoidal distribution between -b and 180.0-b
	      angle -= b;
	      // Sinusoidal distribution between -b/a and (180.0-b)/a
	      angle /= a;
	      // Could add any multiple of 180/a
	      // min = -b/a+nmin*180/a
	      // (min + b/a)*a/180.0 = nmin
	      double nmin, nmax, n;
	      nmin = (min + b/a)*a/180.0;
	      nmax = (max + b/a)*a/180.0;
	      nmin -= 1;
	      nmax += 1;
	      nmin = floor(nmin);
	      nmax = floor(nmax);
	      n = gsl_rng_uniform_int(rand_num_gen, nmax+1-nmin)+nmin;
	      angle += n*180.0/a;
	    }while(angle <= min || angle >= max);
	    sample += angle;
	  }else if(strcasecmp(argv[i], "-Rayleigh") == 0) {
	    double sigma;
	    j = sscanf(argv[i+1], "%lf", &sigma);
	    i++;
	    if(pointnr == 0 && application.verbose_state.verbose) {
	      fprintf(stderr, "Using distribution: x*exp(-x^2/(2*%f^2))/(%f^2)\n", sigma, sigma);
	    }
	    sample += gsl_ran_rayleigh(rand_num_gen, sigma); 
	  }
	  // P = a*b**(-x*c)
	  // P = a*exp(ln(b**(-x*c)))
	  // P = a*exp(-x*c*ln(b))
	  // 1/mu * exp(-x/mu)
	  // Start at pg 230 of gsl ref.
	  if(distr_number == 0) {  // Skip the last distribution functions
	    if(NumberPoints2 > 0) {
	      if(i == cmd_line_end_first_distr)
		break;
	    }
	  }

	}

	// Convolve generated distribution with this additional Gaussian
	if(noisesigma > 0) {
	  sample += gsl_ran_gaussian(rand_num_gen, noisesigma); 
	}
	if(noisefile_id != 0) {
	  long n;
	  n = gsl_rng_uniform_int(rand_num_gen, n_noisedata);
	  sample += data_noise[n];
	}
      
	total += sample;
	fprintf(fout, "%e\n", sample);
      }  // End of point number loop
    }  // End of loop over first and second distribution
    total /= (long double)(NumberPoints+NumberPoints2);
    if(application.verbose_state.verbose) {
      fprintf(stderr, "Average = %Le\n", total);
    }
    

    if(AddNulls) {
      long nrNulls;
      if(total > average_value) {
	//      Etot = NumberPoints*total = average_value*(NumberPoints+nrNulls);
	//      NumberPoints*total/average_value - NumberPoints = nrNulls;
	nrNulls = (NumberPoints+NumberPoints2)*(total/average_value - 1.0);
      }else {
	nrNulls = 0;
      }
      if(nrNulls > 0) {
	double sample;
	if(quiet == 0)
	  fprintf(stderr, "Adding %ld nulls to distribution\n", nrNulls);
	total *= NumberPoints+NumberPoints2;
	for(j = 0; j < nrNulls; j++) {
	  sample = 0;
	  if(noisesigma > 0) {
	    sample += gsl_ran_gaussian(rand_num_gen, noisesigma); 
	  }
	  if(noisefile_id != 0) {
	    long n;
	    n = gsl_rng_uniform_int(rand_num_gen, n_noisedata);
	    sample += data_noise[n];
	  }
	  total += sample;
	  fprintf(fout, "%e\n", sample);
	}
	total /= (long double)(NumberPoints+NumberPoints2+nrNulls);
	if(application.verbose_state.verbose) {
	  fprintf(stderr, "New average = %Le\n", total);
	}
      }else {
	if(quiet == 0) {
	  printwarning(application.verbose_state.debug, "Average energy already less than %e, so cannot make average match the value specified with the -null option", average_value);
	}
      }
    }


    if(outputfile)
      fclose(fout);
  } // End of loopnr for loop
  gsl_rng_free(rand_num_gen); 
  if(noisefile_id != 0) {
    free(data_noise);
  }
  terminateApplication(&application);
  return 0;
}

