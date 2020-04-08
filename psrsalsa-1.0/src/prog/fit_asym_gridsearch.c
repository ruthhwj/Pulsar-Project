#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "psrsalsa.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

int dumpdebuginfo = 0;

int NrSubints = -1;      // Effectively the number of P3 bins
int NrBins;              // These are the number of longitude bins to process
double true_P3 = -1;     // Apparent P3 as observed, in units of pulse periods
double P3;               // This is the number of P3 bins in the P3 fold, value will be the same as NrSubints
double P4;               // Circulation time, in pulse periods, if the observed P3 would be the number of P3 bins in the P3 fold
double alpha;            // Magnetic inclination, in radians
double sin_alpha, cos_alpha;
double beta;             // Impact parameter of LOS, in radians
int fiducial_bin;        // Bin number of the fiducial plane, such that if this is 0, the fiducial plane happens directly after the end of the first bin.
int N_sparks;            // Number of sparks
int alias_order;         // That is the n as the Gupta et al. stuff, i.e. increases by 1 everytime the apparent drift direction reverses.
double sigma1, sigma2;   // Off pulse variance (so really sigma squared) of the two polarizations
int long_bin_start = -1; // An onpulse region which lives within the NrBins longitude bins (start bin)
int long_bin_end = -1;   // An onpulse region which lives within the NrBins longitude bins (end bin)


//double **data;           // The data to be fitted. There are NrBins arrays with length 4*NrSubints, which is the first polarization (O1 leading) at a given pulse longitude, followed by those corresponding to O1 trailing, O2 leading, O2 trailing.
float *data_O1l, *data_O2l, *data_O1t, *data_O2t;  // Data to be fitted. Their length is NrSubints.
long long_bin_nr;        // Current longitude bin which is being fitted for, as required by chi2_funk()
int pol_eq_nr;           // Current polarization equation number (0, 1) which is being fitted for, as required by chi2_funk()

int optimise_without_denominator; // If set nonzero, the denominator of the chi2 is not included in the chi2 optimisation.
int forcephysical;                // If set nonzero, a penalty is applied to try to ensure found solutions are physical


int is_physical(double a, double b, double c, double d);

double chi2_funk(double *x)
{
  double chi2;
  double denomenator, denomenator2;
  chi2 = 0;

  if(forcephysical == 0) {
    if(optimise_without_denominator) {
      denomenator = 1;
    }else {
      if(pol_eq_nr == 0) {
	denomenator = sigma1*(1+x[0]*x[0]) + x[1]*x[1]*sigma2;
      }else {
	denomenator = x[0]*x[0]*sigma1 + sigma2*(1+x[1]*x[1]);
      }
    }
  }else {
    if(optimise_without_denominator) {
      denomenator = 1;
      denomenator2 = 1;
    }else {
      denomenator  = sqrt(sigma1*(1+x[0]*x[0]) + x[1]*x[1]*sigma2);
      denomenator2 = sqrt(x[2]*x[2]*sigma1 + sigma2*(1+x[3]*x[3]));
    }
  }

  int i;
  for(i = 0; i < NrSubints; i++) {
    //    double term;
    //    float Ol1, Ol2, Ot1, Ot2;
    /*
    Ol1 = (data[long_bin_nr])[i+0*NrSubints];
    Ol2 = (data[long_bin_nr])[i+2*NrSubints];
    Ot1 = (data[long_bin_nr])[i+1*NrSubints];
    Ot2 = (data[long_bin_nr])[i+3*NrSubints];
    */
    /*
    Ol1 = data_O1l[i];
    Ol2 = data_O2l[i];
    Ot1 = data_O1t[i];
    Ot2 = data_O2t[i];
    if(pol_eq_nr == 0) {
      term = Ol1-x[0]*Ot1-x[1]*Ot2;
    }else {
      term = Ol2-x[0]*Ot1-x[1]*Ot2;
    }
    */
    if(forcephysical == 0) {
      double term;
      if(pol_eq_nr == 0) {
	term = data_O1l[i]-x[0]*data_O1t[i]-x[1]*data_O2t[i];
      }else {
	term = data_O2l[i]-x[0]*data_O1t[i]-x[1]*data_O2t[i];
      }
      chi2 += term*term;
    }else {
      double term1, term2;
      term1 = (data_O1l[i]-x[0]*data_O1t[i]-x[1]*data_O2t[i])/denomenator;
      term2 = (data_O2l[i]-x[2]*data_O1t[i]-x[3]*data_O2t[i])/denomenator2;
      chi2 += term1*term1 + term2*term2;
      //      printf("XXXXXX %e %e\n", term1, term2);
    }
  }

  //  if(dumpdebuginfo) {
    //    printf("c=%e d=%e chi2=%e", x[0], x[1], chi2/denomenator);
  //  printf("a=%e b=%e c=%e d=%e chi2=%e\n", x[0], x[1], x[2], x[3], chi2);
    /*
    for(i = 0; i < NrSubints; i++) {
      printf(" %e", data_O2l[i]);
    }
    for(i = 0; i < NrSubints; i++) {
      printf(" %e", x[0]*data_O1t[i]+x[1]*data_O2t[i]);
      }*/
    //    printf("\n");
  //  }
  if(forcephysical == 0) {
    return chi2/denomenator;
  }else {
    int physical = is_physical(x[0], x[1], x[2], x[3]);
    if(physical == 0) {
      return chi2*1e10+1e10;
    }
    return chi2;   // Normalisation already applied
  }
}


double calculate_theta_trans(double phi)
{
  // phi is defined with respect to the fiducial plane,
  // i.e. -pi < phi <= pi, where the phi_fid = 0
  //  R = calculate_R(phi);
  double R, sttrans, cttrans;
  double sin_ab = sin(alpha + beta);
  double cos_ab = cos(alpha + beta);
  R = acos( (cos_alpha*cos_ab) + (sin_alpha*sin_ab*cos(phi)));
  double sin_R = sin(R);
  sttrans = (sin_ab*sin(phi)) / sin_R;
  cttrans = ((cos_alpha*cos(R))-cos_ab)/(sin_alpha*sin_R);
  return atan2(sttrans, cttrans);
}


double floating_mod(double a, double b)
{
  // returns the remainder of a - n*b, where n is an integer.
  // Maps a to the range 0 - b.
  if(b < 0.0) {
    fprintf(stderr, "ERROR (floating_mod): the second value must be greater than 0.");
    exit(0);
  }
  /* 
     If a is positive something like this should work
     long n = a/b;
     a -= n*b;
     return a;

     Need to possibly do a +1 or -1 if a is negative

     Following is something I use to get a value between 0 and 360 deg
  int i;
  i = fabs(a)/360.0;
  if(a > 0)
    a -= 360*i;
  else
    a += 360*i;
  if(a < 0)
    a += 360;
  return a;


   */
  while(a >= b) {
    a -= b;
  }
  while(a < 0.0) {
    a += b;
  }
  return a;
}

void find_opposite_subint(int sub_nr, int bin_nr, double sign_modifier, int *sub1, double *weight1, int *sub2, double *weight2, double *true_opposite)
{
  // k_2 = \bigg( k_1 + \frac{P_3}{2\pi}\Delta\Phi\bigg)\mod P_3
  // \Delta\Phi = -N\bigg( \frac{\Delta\phi}{P_4}\mp\Delta\theta_\mathrm{trans}\bigg)
 
  // Calculate the pulse longitude wrt the fiducial plane
  // (\phi - \phi_\mathrm{fid})
  double phi;
  phi = (((double)bin_nr) + 0.5)*2.0*M_PI/((double)NrBins)-M_PI;

  //Remember that ttrans is symmetric, (ttrans2 = -ttrans1)
  // so Delta ttrans = -ttrans1 - ttrans 1 = -2ttrans
  double ttrans;
  ttrans = calculate_theta_trans(phi);

  double Dphi, Dttrans, DPhase;
  Dphi = -2.0*phi;
  Dttrans = -2.0*ttrans;
  double D;
  if(alias_order % 2 == 0) {
    D = 1.0;
  }else {
    D = -1.0;
  }
  // Changed the +sign_modifier to -sign_modifier, and changed it back again
  DPhase = -((double)N_sparks)*((Dphi/P4) + sign_modifier*Dttrans)*D;

  // find the opposite number.
  // Remember, this is the true opposite number + 0.5 because of the way bins work.
  double opposite_sub_nr;
  opposite_sub_nr = floating_mod(((double)sub_nr + DPhase*P3/(2.0*M_PI)), P3);
  *true_opposite = opposite_sub_nr - 0.5;
  // Return 2 subints with their associated weights, and the true opposite sub for a sanity check.

  double below;
  below = floor(opposite_sub_nr);

  *weight1 = 1.0 - (opposite_sub_nr - below);
  *weight2 = (opposite_sub_nr - below);
    // But actually the centres of the bins are non-integer...
  *sub1 = (int)(floating_mod(floor(*true_opposite),NrSubints));
  *sub2 = (int)(floating_mod(ceil(*true_opposite), NrSubints));
}

void find_initial_guess(double *guess_a, double *guess_b, double *guess_c, double *guess_d)
{
  int n;
  double Sxx, Sxy1, Sxy2, Szz, Szy1, Szy2, Sxz;
  Sxx = Sxy1 = Sxy2 = Szz = Szy1 = Szy2 = Sxz = 0;
  for(n = 0; n < NrSubints; n++) {
    double O1t = data_O1t[n];
    double O1l = data_O1l[n];
    double O2t = data_O2t[n];
    double O2l = data_O2l[n];
    //    printf("find_initial_guess: %e\t%e\t%e\t%e\n", O1l, O1t, O2l, O2t);
    Sxx  += O1t*O1t;
    Sxy1 += O1t*O1l;
    Sxy2 += O1t*O2l;
    Szz  += O2t*O2t;
    Szy1 += O2t*O1l; 
    Szy2 += O2t*O2l;
    Sxz  += O1t*O2t;
  }
  double det;
  det = Sxx*Szz-Sxz*Sxz;
  *guess_a = (Sxy1*Szz-Sxz*Szy1)/det;
  *guess_b = (Sxx*Szy1-Sxz*Sxy1)/det;
  *guess_c = (Sxy2*Szz-Sxz*Szy2)/det;
  *guess_d = (Sxx*Szy2-Sxz*Sxy2)/det;
}

int is_physical(double a, double b, double c, double d)
{
  int PHYSICAL = 1;
  if(!(((a < 0.0) && (b < 0.0)) || ((c < 0.0) && (d < 0.0)))) {
    //    if((a*fabs(c)) == (c*fabs(a))) {
    if(a*c >= 0) {
      PHYSICAL = 1.0;
    }else if((fabs(a)*d)+(b*fabs(c)) >= 0) {
      PHYSICAL = 1.0;
    }else {
      PHYSICAL = 0.0;
    }
  }else {
    PHYSICAL = 0.0;
  }
  return PHYSICAL;
}

int main(int argc, char **argv)
{
  long i;
  int profiledump, solutiondump;
  psrsalsaApplication application;
  double posmatrix_level;

  //  sigma1 = sigma2 = -1;

  initApplication(&application, "fit_asym", "[options] inputfile1 inputfile2");
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_history_cmd_only = 1;
  profiledump = 0;
  solutiondump = 0;
  optimise_without_denominator = 0;
  forcephysical = 0;
  posmatrix_level = -1;

  if(argc < 2) {
    printf("Program to fit the Asymmetry matrix, to do Crispin's 0031 stuff. The first input file is the data (P3-fold with three polarization channels, which are Stokes I and the two mode separated polarizations). The second input file defines the gridsearch parameters.\n");
    printApplicationHelp(&application);
    printf("Input options:\n\n");
    //    printf("-sigma   \"sigma1 sigma2\"    Define the two RMSes of the two polarizations.\n");
    //    printf("-np3bins n                    Set the expected the number of P3 bins to be n.\n");
    printf("-P3 value                 The apparent P3 as observed, in units of pulse periods.\n");
    printf("-onpulse \"bin1 bin2\"      The range of bins considered as onpulse (used to calculate the off-pulse rms).\n");
    printf("-profiledump              Produce dump files corresponding to the pulse profile for each itteration (it is actually the standard deviation of each column in the P3-fold).\n");
    printf("-solutiondump             Produce dump files corresponding to the solution of the fitting process in terms of the data. Two polarizations are produced, corresponding to the two observed mode-splitted polarizations. The leading half of the profile corresponds to the actually observed data, the trailing half to the by the asymetry modified trailing half data, which therefore should match what is observed in the leading half.\n");
    printf("-simplechi2               If specified, the denominator of the chi2 function is ignored during optimisation, so it effectively becomes an unweighted fit.\n");
    printf("-forcephysical            Ensure that only physical solutions are found\n");
    printf("-posmatrix level          Rather than doing a gridsearch, read the solutions from the matrix.DUMP file, and make a posmatrix.DUMP, which should only have positive matrix elements. The input files are expected to be the same as for the gridsearch which made the matrix.DUMP file. The variable level (between 0 and 1) indicates what solution is picked out of the range of possibilities. Each column in the mixing matrix is treated independently, and the sum of the squares of the elements in the trailing half mixing matrix is defined to be 1.\n");
    printf("\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else {

    for(i = 1; i < argc; i++) {
      int index;
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	      i = index;
      }else if(strcmp(argv[i], "-P3") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &true_P3, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fit_asym: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-onpulse") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &long_bin_start, &long_bin_end, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fit_asym: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
	/*
      }else if(strcmp(argv[i], "-sigma") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &sigma1, &sigma2, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fit_asym: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
	}else if(strcmp(argv[i], "-np3bins") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &NrSubints, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fit_asym: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
	*/
      }else if(strcmp(argv[i], "-profiledump") == 0) {
	profiledump = 1;
      }else if(strcmp(argv[i], "-solutiondump") == 0) {
	solutiondump = 1;
      }else if(strcmp(argv[i], "-simplechi2") == 0) {
	optimise_without_denominator = 1;
      }else if(strcmp(argv[i], "-forcephysical") == 0) {
	forcephysical = 1;
      }else if(strcmp(argv[i], "-posmatrix") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &posmatrix_level, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fit_asym: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	if(posmatrix_level < 0.0 || posmatrix_level > 1.0) {
	  printerror(application.verbose_state.debug, "ERROR fit_asym: Please use specify a number between 0 and 1 for option %s.", argv[i]);
	  return 0;
	}
	i++;
      }else {
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "ERROR fit_asym: Unknown option: %s\n\nRun fit_asym without command line arguments to show help", argv[i]);
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0)
	    return 0;
	}
      }
    }
  }

  if(true_P3 <= 0 && posmatrix_level < 0.0) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Please use the -P3 option to provide a positive number.");
    return 0;
  }
  if((long_bin_start < 0 || long_bin_start < 0) && posmatrix_level < 0.0) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Please use the -onpulse option to provide a positive number.");
    return 0;
  }


  /*
  if(sigma1 <= 0 || sigma2 <= 0) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Please use the -sigma option to provide two positive numbers.");
    return 0;
  }
  if(NrSubints <= 0) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Please use the -np3bins option to provide a positive number.");
    return 0;
  }
  */
  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Two input files were expected at the end of the command line.");
    return 0;
  }


  char *filename_ptr;
  filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state);
  if(application.verbose_state.verbose) {
    fprintf(stdout, "Loading pulsar data from file %s\n", filename_ptr);
  }
  datafile_definition datain;
  // Open a file in a non-specified format (will be hopefully figured out automatically), writing is not enabled, and the file is read into memory.
  if(openPSRData(&datain, filename_ptr, 0, 0, 1, 0, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR: Error opening data");
    return 0;
  }
  NrSubints = datain.NrSubints;      // Effectively the number of P3 bins
  NrBins = datain.NrBins;            // These are the number of longitude bins to process
  if(datain.NrPols != 3) {
    printerror(application.verbose_state.debug, "ERROR: Expected that the data file contains 3 polarizations.");
    return 0;
  }
  if(datain.NrBins % 2 == 1) {
    printerror(application.verbose_state.debug, "ERROR: Expected that the data file contains an even number of pulse longitude bins.");
    return 0;
  }
  
  filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state);
  if(application.verbose_state.verbose) {
    fprintf(stdout, "Loading grid search information from %s\n", filename_ptr);
  }
  FILE *fin;
  fin = fopen(filename_ptr, "r");
  if(fin == NULL) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: cannot open file %s.\n", filename_ptr);
    return 0;
  }
  char txt[1000];
  int ret, linenr;
  double alpha_start, alpha_end;    // Start/end alpha in degrees
  int alpha_steps;                 // Nr of steps in alpha gridsearch
  double drift_grad_start, drift_grad_end;   // Gradient of the drift bands in P3 bins per degree????
  int drift_grad_steps;
  int nrsparks_start, nrsparks_end, nrsparks_steps;
  int fid_bin_start, fid_bin_end, fid_bin_steps;
  int aliasorder_start, aliasorder_end, aliasorder_steps;


  linenr = 1;
  ret = fscanf(fin, "%lf %s", &alpha_start, txt);
  if(ret != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Reading line %d failed in %s.\n", linenr, filename_ptr);
    return 0;    
  }
  linenr++;

  ret = fscanf(fin, "%lf %s", &alpha_end, txt);
  if(ret != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Reading line %d failed in %s.\n", linenr, filename_ptr);
    return 0;    
  }
  linenr++;

  ret = fscanf(fin, "%d %s", &alpha_steps, txt);
  if(ret != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Reading line %d failed in %s.\n", linenr, filename_ptr);
    return 0;    
  }
  linenr++;

  ret = fscanf(fin, "%lf %s", &drift_grad_start, txt);
  if(ret != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Reading line %d failed in %s.\n", linenr, filename_ptr);
    return 0;    
  }
  linenr++;

  ret = fscanf(fin, "%lf %s", &drift_grad_end, txt);
  if(ret != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Reading line %d failed in %s.\n", linenr, filename_ptr);
    return 0;    
  }
  linenr++;

  ret = fscanf(fin, "%d %s", &drift_grad_steps, txt);
  if(ret != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Reading line %d failed in %s.\n", linenr, filename_ptr);
    return 0;    
  }
  linenr++;

  ret = fscanf(fin, "%d %s", &nrsparks_start, txt);
  if(ret != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Reading line %d failed in %s.\n", linenr, filename_ptr);
    return 0;    
  }
  linenr++;

  ret = fscanf(fin, "%d %s", &nrsparks_end, txt);
  if(ret != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Reading line %d failed in %s.\n", linenr, filename_ptr);
    return 0;    
  }
  linenr++;

  ret = fscanf(fin, "%d %s", &nrsparks_steps, txt);
  if(ret != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Reading line %d failed in %s.\n", linenr, filename_ptr);
    return 0;    
  }
  linenr++;

  if((nrsparks_end-nrsparks_start+1)%nrsparks_steps != 0) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: The grid search over the number of sparks (start=%d, end=%d and steps=%d) are not compatible since the number of sparks should be an integer.\n", nrsparks_start, nrsparks_end, nrsparks_steps);
    return 0;        
  }


  ret = fscanf(fin, "%d %s", &fid_bin_start, txt);
  if(ret != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Reading line %d failed in %s.\n", linenr, filename_ptr);
    return 0;    
  }
  linenr++;

  ret = fscanf(fin, "%d %s", &fid_bin_end, txt);
  if(ret != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Reading line %d failed in %s.\n", linenr, filename_ptr);
    return 0;    
  }
  linenr++;

  ret = fscanf(fin, "%d %s", &fid_bin_steps, txt);
  if(ret != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Reading line %d failed in %s.\n", linenr, filename_ptr);
    return 0;    
  }
  linenr++;

  ret = fscanf(fin, "%d %s", &aliasorder_start, txt);
  if(ret != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Reading line %d failed in %s.\n", linenr, filename_ptr);
    return 0;    
  }
  linenr++;

  ret = fscanf(fin, "%d %s", &aliasorder_end, txt);
  if(ret != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Reading line %d failed in %s.\n", linenr, filename_ptr);
    return 0;    
  }
  linenr++;

  ret = fscanf(fin, "%d %s", &aliasorder_steps, txt);
  if(ret != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Reading line %d failed in %s.\n", linenr, filename_ptr);
    return 0;    
  }
  linenr++;

  if((aliasorder_end-aliasorder_start+1)%aliasorder_steps != 0) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: The grid search over the alias order (start=%d, end=%d and steps=%d) are not compatible since the alias order should be an integer.\n", aliasorder_start, aliasorder_end, aliasorder_steps);
    return 0;        
  }

  fclose(fin);
  if(application.verbose_state.verbose) {
    printf("alpha (deg): %lf %lf in %d steps\n", alpha_start, alpha_end, alpha_steps);
    printf("drift gradient (?): %lf %lf in %d steps\n", drift_grad_start, drift_grad_end, drift_grad_steps);
    printf("Nr sparks: %d %d in %d steps\n", nrsparks_start, nrsparks_end, nrsparks_steps);
    printf("Fiducial bin: %d %d in %d steps\n", fid_bin_start, fid_bin_end, fid_bin_steps);
    printf("Alias order:  %d %d in %d steps\n", aliasorder_start, aliasorder_end, aliasorder_steps);
  }
  P3 = (double)NrSubints;

  long p3bin, polnr, longbin;
  for(polnr = 1; polnr <= 2; polnr++) {
    double variance;
    long nrpoints_used;
    variance = 0;
    nrpoints_used = 0;
    for(p3bin = 0; p3bin < NrSubints; p3bin++) {
      for(longbin = 0; longbin < NrBins; longbin++) {
	if(longbin < long_bin_start || longbin > long_bin_end) {
	  nrpoints_used++;
	  float sample;
	  //	  fprintf(stderr, "%ld %ld %ld\n", p3bin, longbin, polnr);
	  if(readPulsePSRData(&datain, p3bin, polnr, 0, longbin, 1, &sample, application.verbose_state) != 1) {
	    printerror(application.verbose_state.debug, "ERROR: read error");
	    return 0;
	  }
	  variance += sample*sample;
	}
      }
    }
    if(polnr == 1) {
      sigma1 = variance/(double)nrpoints_used;
    }else {
      sigma2 = variance/(double)nrpoints_used;
    }
  }
  if(application.verbose_state.verbose) {
    printf("Found offpulse RMS for the two polariations: sqrt(sigma1)=%e sqrt(sigma2)=%e\n", sqrt(sigma1), sqrt(sigma2));
  }  

  alpha_start *= M_PI/180.0;
  alpha_end *= M_PI/180.0;

  double alpha_stepsize, drift_grad_stepsize, drift_grad;
  int nrsparks_stepsize, fid_bin_stepsize, aliasorder_stepsize;

  if(alpha_steps >= 2) {
    alpha_stepsize = (alpha_end-alpha_start)/(double)(alpha_steps-1);
  }else {
    alpha_stepsize = 1;
    alpha_end = alpha_start;
  }

  if(drift_grad_steps >= 2) {
    drift_grad_stepsize = (drift_grad_end-drift_grad_start)/(double)(drift_grad_steps-1);
  }else {
    drift_grad_stepsize = 1;
    drift_grad_end = drift_grad_start;
  }

  if(nrsparks_steps >= 2) {
    nrsparks_stepsize = (nrsparks_end-nrsparks_start)/(nrsparks_steps-1);
  }else {
    nrsparks_stepsize = 1;
    nrsparks_end = nrsparks_start;
  }

  if(fid_bin_steps >= 2) {
    fid_bin_stepsize = (fid_bin_end-fid_bin_start)/(fid_bin_steps-1);
  }else {
    fid_bin_stepsize = 1;
    fid_bin_end = fid_bin_start;
  }

  if(aliasorder_steps >= 2) {
    aliasorder_stepsize = (aliasorder_end-aliasorder_start)/(aliasorder_steps-1);
  }else {
    aliasorder_stepsize = 1;
    aliasorder_end = aliasorder_start;
  }

  datafile_definition *rotated_data;
  rotated_data = malloc(fid_bin_steps*sizeof(datafile_definition));
  if(rotated_data == NULL) {
    printerror(application.verbose_state.debug, "ERROR fit_asym_gridsearch: Memory allocation error.\n");
    return 0;
  }
  float **rotated_profiles;
  if(profiledump) {
    rotated_profiles = malloc(fid_bin_steps*sizeof(float *));
    if(rotated_profiles == NULL) {
      printerror(application.verbose_state.debug, "ERROR fit_asym_gridsearch: Memory allocation error.\n");
      return 0;
    }
    for(i = 0; i < fid_bin_steps; i++) {
      rotated_profiles[i] = malloc(NrBins*sizeof(float));
      if(rotated_profiles[i] == NULL) {
	printerror(application.verbose_state.debug, "ERROR fit_asym_gridsearch: Memory allocation error.\n");
	return 0;
      }
    }
  }


  int fiducial_bin_index = 0;
  if(posmatrix_level < 0.0) {
    printf("Setting up some stuff before fitting starts\n");

    // Make different copies of the data, which has the fiducial plane in the centre
    for(fiducial_bin = fid_bin_start; fiducial_bin <= fid_bin_end; fiducial_bin += fid_bin_stepsize) {
      
      // Clone the original data
      if(make_clone(datain, &(rotated_data[fiducial_bin_index]), application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR fit_asym_gridsearch: Cloning data failed.\n");
	return 0;
      }
      
      float shiftPhase = (NrBins/2.0 - 1.0 - fiducial_bin)/(float)NrBins;
      if(preprocess_fftshift(rotated_data[fiducial_bin_index], shiftPhase, 0, 0.0, application.verbose_state) != 1) {
	printerror(application.verbose_state.debug, "ERROR fit_asym_gridsearch: Rotating data failed.\n");
	return 0;
      }
      //    printf("index = %d - Apply shift = %e bins\n", fiducial_bin_index, shiftPhase);
      if(profiledump) {
	int bin_nr;
	for(bin_nr = 0; bin_nr < NrBins; bin_nr++) {
	  (rotated_profiles[fiducial_bin_index])[bin_nr] = 0;
	  for(i = 0; i < NrSubints; i++) {
	    float sample;
	    if(readPulsePSRData(&(rotated_data[fiducial_bin_index]), i, 0, 0, bin_nr, 1, &sample, application.verbose_state) != 1) {
	      printerror(application.verbose_state.debug, "ERROR: read error");
	      return 0;
	    }
	    (rotated_profiles[fiducial_bin_index])[bin_nr] += sample*sample;
	  }
	  (rotated_profiles[fiducial_bin_index])[bin_nr] = sqrt((rotated_profiles[fiducial_bin_index])[bin_nr])/(float)NrSubints;
	}
      }
      fiducial_bin_index++;
    }
  }

  datafile_definition output_data_file_solution;
  if(solutiondump) {
    // When writing out data, start with initialisng some variables in the output data file struct
    cleanPSRData(&output_data_file_solution, application.verbose_state);
    copy_params_PSRData(datain, &output_data_file_solution, application.verbose_state);
    // Since in this example we only write out a single subintegration, make sure that information is correct in the header
    output_data_file_solution.NrFreqChan = 1;
    output_data_file_solution.NrPols = 2;
    output_data_file_solution.data = malloc(datain.NrBins*output_data_file_solution.NrPols*output_data_file_solution.NrSubints*sizeof(float));
    if(output_data_file_solution.data == NULL) {
      printerror(application.verbose_state.debug, "ERROR: Cannot allocate memory.");
      return 0;
    }
  }


  // Close the input file
  closePSRData(&datain, 0, application.verbose_state);

  FILE *file_grid, *file_matrix, *file_physical, *file_chisq, *file_profile;
  if(posmatrix_level < 0) {
    file_grid = fopen("gridsearch.DUMP", "w");
    file_matrix = fopen("matrix.DUMP", "w");
    file_physical = fopen("physical.DUMP", "w");
    file_chisq = fopen("chisq.DUMP", "w");
    if(file_grid == NULL || file_matrix == NULL || file_physical == NULL || file_chisq == NULL) {
      printerror(application.verbose_state.debug, "ERROR fit_asym_gridsearch: Cannot open output ascii file.\n");
      perror("");
      return 0;
    }
    if(profiledump) {
      file_profile = fopen("profile.DUMP", "w");
      if(file_profile == NULL) {
	printerror(application.verbose_state.debug, "ERROR fit_asym_gridsearch: Cannot open output ascii file.\n");
	perror("");
	return 0;
      }
    }
    fprintf(file_grid, "%f %f %d %f %f %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", alpha_start*180.0/M_PI, alpha_end*180.0/M_PI, alpha_steps, drift_grad_start, drift_grad_end, drift_grad_steps, nrsparks_start, nrsparks_end, nrsparks_steps, fid_bin_start, fid_bin_end, fid_bin_steps, aliasorder_start, aliasorder_end, aliasorder_steps, long_bin_start, long_bin_end, NrSubints, NrBins);
    fprintf(file_matrix, "%f %f %d %f %f %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", alpha_start*180.0/M_PI, alpha_end*180.0/M_PI, alpha_steps, drift_grad_start, drift_grad_end, drift_grad_steps, nrsparks_start, nrsparks_end, nrsparks_steps, fid_bin_start, fid_bin_end, fid_bin_steps, aliasorder_start, aliasorder_end, aliasorder_steps, long_bin_start, long_bin_end, NrSubints, NrBins);
    fprintf(file_physical, "%f %f %d %f %f %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", alpha_start*180.0/M_PI, alpha_end*180.0/M_PI, alpha_steps, drift_grad_start, drift_grad_end, drift_grad_steps, nrsparks_start, nrsparks_end, nrsparks_steps, fid_bin_start, fid_bin_end, fid_bin_steps, aliasorder_start, aliasorder_end, aliasorder_steps, long_bin_start, long_bin_end, NrSubints, NrBins);
    fprintf(file_chisq, "%f %f %d %f %f %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", alpha_start*180.0/M_PI, alpha_end*180.0/M_PI, alpha_steps, drift_grad_start, drift_grad_end, drift_grad_steps, nrsparks_start, nrsparks_end, nrsparks_steps, fid_bin_start, fid_bin_end, fid_bin_steps, aliasorder_start, aliasorder_end, aliasorder_steps, long_bin_start, long_bin_end, NrSubints, NrBins);
    if(profiledump) {
      fprintf(file_profile, "%f %f %d %f %f %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", alpha_start*180.0/M_PI, alpha_end*180.0/M_PI, alpha_steps, drift_grad_start, drift_grad_end, drift_grad_steps, nrsparks_start, nrsparks_end, nrsparks_steps, fid_bin_start, fid_bin_end, fid_bin_steps, aliasorder_start, aliasorder_end, aliasorder_steps, long_bin_start, long_bin_end, NrSubints, NrBins);
    }
  }

  float *data_O1t_notrotated, *data_O2t_notrotated;
  double *solution_a, *solution_b, *solution_c, *solution_d, *chi2_1, *chi2_2;
  int *physical, *fit_status_1, *fit_status_2;
  data_O1l = malloc(NrSubints*sizeof(float));
  data_O2l = malloc(NrSubints*sizeof(float));
  data_O1t = malloc(NrSubints*sizeof(float));
  data_O2t = malloc(NrSubints*sizeof(float));
  solution_a = malloc(NrBins*sizeof(double));
  solution_b = malloc(NrBins*sizeof(double));
  solution_c = malloc(NrBins*sizeof(double));
  solution_d = malloc(NrBins*sizeof(double));
  chi2_1 = malloc(NrBins*sizeof(double));
  chi2_2 = malloc(NrBins*sizeof(double));
  physical = malloc(NrBins*sizeof(int));
  fit_status_1 = malloc(NrBins*sizeof(int));
  fit_status_2 = malloc(NrBins*sizeof(int));
  data_O1t_notrotated = malloc(NrSubints*sizeof(float));
  data_O2t_notrotated = malloc(NrSubints*sizeof(float));
  if(data_O1l == NULL || data_O2l == NULL  || data_O1t == NULL || data_O2t == NULL || data_O1t_notrotated == NULL || data_O2t_notrotated == NULL || solution_a == NULL || solution_b == NULL || solution_c == NULL || solution_c == NULL || chi2_1 == NULL || chi2_1 == NULL || physical == NULL || fit_status_1 == NULL || fit_status_2 == NULL) {
    printerror(application.verbose_state.debug, "ERROR fit_asym_gridsearch: Memory allocation error.\n");
    return 0;
  }
  for(i = NrBins/2; i < NrBins; i++) {
    solution_a[i] = 1;
    solution_b[i] = 0;
    solution_c[i] = 0;
    solution_d[i] = 1;
    chi2_1[i] = 0;
    chi2_2[i] = 0;
    physical[i] = 1;
    fit_status_1[i] = 0;
    fit_status_2[i] = 0;
  }




  if(posmatrix_level >= 0.0) {
    FILE *file_posmatrix;
    file_matrix = fopen("matrix.DUMP", "r");
    file_posmatrix = fopen("posmatrix.DUMP", "w");
    if(file_matrix == NULL || file_posmatrix == NULL) {
      printerror(application.verbose_state.debug, "ERROR fit_asym_gridsearch: Cannot open output/input ascii file.\n");
      perror("");
      return 0;
    }
    char line[10001];
    fgets(line, 10000, file_matrix);
    fputs(line, file_posmatrix);

    printf("Start finding positive matrices process\n");
    fflush(stdout);



    gsl_rng *rand_num_gen;
    const gsl_rng_type *rand_num_gen_type;
    int repeatable = 1;    // If set, lines that contain random numbers etc are not produced
    gsl_rng_env_setup();    /* Set the default generators, can be influenced by environment variables */
    rand_num_gen_type = gsl_rng_default;
    rand_num_gen = gsl_rng_alloc(rand_num_gen_type);
    long idnum;
    if(repeatable == 0)
      randomize_idnum(&idnum);
    else
      idnum = 123;
    gsl_rng_set(rand_num_gen, idnum);

    double *best_a, *best_b, *best_c, *best_d;
    best_a = malloc(NrBins*sizeof(double));
    best_b = malloc(NrBins*sizeof(double));
    best_c = malloc(NrBins*sizeof(double));
    best_d = malloc(NrBins*sizeof(double));
    if(best_a == NULL || best_b == NULL || best_c == NULL || best_c == NULL) {
      printerror(application.verbose_state.debug, "ERROR fit_asym_gridsearch: Memory allocation error.\n");
      return 0;
    }


    //    double alpha_iterator = 0.0;
    int alpha_index;
    int drift_grad_index;

    //    for(alpha_iterator = alpha_start; alpha_iterator <= alpha_end; alpha_iterator += alpha_stepsize) {
    for(alpha_index = 0; alpha_index < alpha_steps; alpha_index++) {
      //      for(drift_grad = drift_grad_start; drift_grad <= drift_grad_end; drift_grad += drift_grad_stepsize) {
      for(drift_grad_index = 0; drift_grad_index < drift_grad_steps; drift_grad_index++) {
	for(N_sparks = nrsparks_start; N_sparks <= nrsparks_end; N_sparks += nrsparks_stepsize) {
	  for(fiducial_bin = fid_bin_start; fiducial_bin <= fid_bin_end; fiducial_bin += fid_bin_stepsize) {
	    for(alias_order = aliasorder_start; alias_order <= aliasorder_end; alias_order += aliasorder_stepsize) {
	      int bin_nr, ret;
	      for(bin_nr = 0; bin_nr < NrBins; bin_nr++) {
		ret = fscanf(file_matrix, "%lf", &(solution_a[bin_nr]));
		if(ret != 1) {
		  printerror(application.verbose_state.debug, "ERROR fit_asym_gridsearch: Cannot read input ascii file.\n");
		  return 0;
		}
	      }
	      for(bin_nr = 0; bin_nr < NrBins; bin_nr++) {
		ret = fscanf(file_matrix, "%lf", &(solution_b[bin_nr]));
		if(ret != 1) {
		  printerror(application.verbose_state.debug, "ERROR fit_asym_gridsearch: Cannot read input ascii file.\n");
		  return 0;
		}
	      }
	      for(bin_nr = 0; bin_nr < NrBins; bin_nr++) {
		ret = fscanf(file_matrix, "%lf", &(solution_c[bin_nr]));
		if(ret != 1) {
		  printerror(application.verbose_state.debug, "ERROR fit_asym_gridsearch: Cannot read input ascii file.\n");
		  return 0;
		}
	      }
	      for(bin_nr = 0; bin_nr < NrBins; bin_nr++) {
		ret = fscanf(file_matrix, "%lf", &(solution_d[bin_nr]));
		if(ret != 1) {
		  printerror(application.verbose_state.debug, "ERROR fit_asym_gridsearch: Cannot read input ascii file.\n");
		  return 0;
		}
	      }

	      // The default solution is what it was originally
	      for(bin_nr = 0; bin_nr < NrBins; bin_nr++) {
		best_a[bin_nr] = solution_a[bin_nr];
		best_b[bin_nr] = solution_b[bin_nr];
		best_c[bin_nr] = solution_c[bin_nr];
		best_d[bin_nr] = solution_d[bin_nr];
	      }

	      for(bin_nr = 0; bin_nr < NrBins/2; bin_nr++) {
		double la, lb, lc, ld;
		double ta, tb, tc, td;

		int opposite_bin_nr;
		opposite_bin_nr = NrBins - 1 - bin_nr;  // opposite bin is trivial
		
		int found_sol = 0;
		double xmin, xmin1, xmin2;
		double xmax, xmax1, xmax2;
		if(solution_a[bin_nr] >= 0.0 && solution_c[bin_nr] >= 0.0) {
		  xmin1 = -solution_b[bin_nr]/solution_a[bin_nr];
		  xmin2 = -solution_d[bin_nr]/solution_c[bin_nr];
		  if(xmin1 < 0) {
		    xmin1 = 0;
		  }
		  if(xmin2 < 0) {
		    xmin2 = 0;
		  }
		  xmax1 = -1;
		  xmax2 = -1;
		}else if(solution_a[bin_nr] >= 0.0 && solution_c[bin_nr] < 0.0) {
		  xmin1 = xmin2 = -solution_b[bin_nr]/solution_a[bin_nr];
		  xmax1 = xmax2 = -solution_d[bin_nr]/solution_c[bin_nr];
		  if(xmin1 < 0) {
		    xmin1 = 0;
		    xmin2 = 0;
		  }
		}else if(solution_a[bin_nr] < 0.0 && solution_c[bin_nr] >= 0.0) {
		  xmax1 = xmax2 = -solution_b[bin_nr]/solution_a[bin_nr];
		  xmin1 = xmin2 = -solution_d[bin_nr]/solution_c[bin_nr];
		  if(xmin1 < 0) {
		    xmin1 = 0;
		    xmin2 = 0;
		  }
		}else {
		  xmax1 = -solution_b[bin_nr]/solution_a[bin_nr];
		  xmax2 = -solution_d[bin_nr]/solution_c[bin_nr];
		  xmin1 = xmin2 = -1;
		  //		  found_sol = 1;
		}
		if(xmin2 > xmin1) {
		  xmin = xmin2;
		}else {
		  xmin = xmin1;
		}

		if(xmax2 > xmax1) {
		  xmax = xmax1;
		}else {
		  xmax = xmax2;
		}

		double xmin_norm, xmax_norm;
		xmin_norm = xmin/sqrt(xmin*xmin+1.0);
		xmax_norm = xmax/sqrt(xmax*xmax+1.0);

		if(xmin_norm < 0.0) {
		  xmin_norm = 0.0;
		}
		if(xmax_norm < 0.0) {
		  xmax_norm = 1.0;
		}

		ta = xmin_norm + posmatrix_level*(xmax_norm-xmin_norm);
		tb = xmin_norm + posmatrix_level*(xmax_norm-xmin_norm);
		tc = sqrt(1.0-ta*ta);
		td = sqrt(1.0-tb*tb);

		la = solution_a[bin_nr]*ta + solution_b[bin_nr]*tc;
		lb = solution_a[bin_nr]*tb + solution_b[bin_nr]*td;
		lc = solution_c[bin_nr]*ta + solution_d[bin_nr]*tc;
		ld = solution_c[bin_nr]*tb + solution_d[bin_nr]*td;

		best_a[bin_nr] = la;
		best_b[bin_nr] = lb;
		best_c[bin_nr] = lc;
		best_d[bin_nr] = ld;
		best_a[opposite_bin_nr] = ta;
		best_b[opposite_bin_nr] = tb;
		best_c[opposite_bin_nr] = tc;
		best_d[opposite_bin_nr] = td;
		if(la >= 0 && lb >= 0 && lc >= 0 && ld >= 0) {
		  found_sol = 1;
		}
		if(found_sol == 0) {
		  printf("Finding positive matrix failed for bin %d: %e %e %e %e physical=%d\n", bin_nr, solution_a[bin_nr], solution_b[bin_nr], solution_c[bin_nr], solution_d[bin_nr], is_physical(solution_a[bin_nr], solution_b[bin_nr], solution_c[bin_nr], solution_d[bin_nr]));
		}

		/*
		long itt;
		double ampl = 1e-0;
		int found_sol = 0;
		for(itt = 0; itt < 100000; itt++) {
		  ta = ampl*gsl_rng_uniform(rand_num_gen);
		  tb = ampl*gsl_rng_uniform(rand_num_gen);
		  tc = ampl*gsl_rng_uniform(rand_num_gen);
		  td = ampl*gsl_rng_uniform(rand_num_gen);
		  la = solution_a[bin_nr]*ta + solution_b[bin_nr]*tc;
		  lb = solution_a[bin_nr]*tb + solution_b[bin_nr]*td;
		  lc = solution_c[bin_nr]*ta + solution_d[bin_nr]*tc;
		  ld = solution_c[bin_nr]*tb + solution_d[bin_nr]*td;
		  if(la >= 0 && lb >= 0 && lc >= 0 && ld >= 0) {
		    best_a[bin_nr] = la;
		    best_b[bin_nr] = lb;
		    best_c[bin_nr] = lc;
		    best_d[bin_nr] = ld;
		    best_a[opposite_bin_nr] = ta;
		    best_b[opposite_bin_nr] = tb;
		    best_c[opposite_bin_nr] = tc;
		    best_d[opposite_bin_nr] = td;
		    found_sol = 1;
		  }
		}
		if(found_sol == 0) {
		  printf("Finding positive matrix failed for bin %d: %e %e %e %e physical=%d\n", bin_nr, solution_a[bin_nr], solution_b[bin_nr], solution_c[bin_nr], solution_d[bin_nr], is_physical(solution_a[bin_nr], solution_b[bin_nr], solution_c[bin_nr], solution_d[bin_nr]));
		}
		*/

	      }   // End of bin_nr loop

	      for(bin_nr = 0; bin_nr < NrBins; bin_nr++) {
		fprintf(file_posmatrix, "%e ", best_a[bin_nr]);
	      }
	      for(bin_nr = 0; bin_nr < NrBins; bin_nr++) {
		fprintf(file_posmatrix, "%e ", best_b[bin_nr]);
	      }
	      for(bin_nr = 0; bin_nr < NrBins; bin_nr++) {
		fprintf(file_posmatrix, "%e ", best_c[bin_nr]);
	      }
	      for(bin_nr = 0; bin_nr < NrBins; bin_nr++) {
		if(bin_nr == NrBins - 1) {
		  fprintf(file_posmatrix, "%e\n", best_d[bin_nr]);
		}else {
		  fprintf(file_posmatrix, "%e ", best_d[bin_nr]);
		}
	      }
	    }
	  }
	}
      }
    }

    fclose(file_posmatrix);
    fclose(file_matrix);
    return 0;
  }    // End of posmatrix stuff



  printf("Start fitting process\n");
  fflush(stdout);

  long current_itteration = 0;
  //  double alpha_iterator = 0.0;

  int alpha_index;
  //  for(alpha_iterator = alpha_start; alpha_iterator <= alpha_end; alpha_iterator += alpha_stepsize) {
  for(alpha_index = 0; alpha_index < alpha_steps; alpha_index++) {
    alpha = alpha_start + alpha_index*alpha_stepsize;
    // sanity check: alpha can never be zero (for the maths to work), so if it is
    // found to be 0, set it to something small but non-zero.
    //    alpha = alpha_iterator;
    if (alpha <= 0.0){
      alpha = 1.0E-5;
    }

    double sign_modifier;
    //    sign_modifier = -abs(np.cos(alpha))/np.cos(alpha)
    if(alpha < 0.5*M_PI) {
      sign_modifier = -1.0;
    }else {
      sign_modifier = 1.0;
    }

    sin_alpha = sin(alpha);
    cos_alpha = cos(alpha);

    int drift_grad_index;
    //    for(drift_grad = drift_grad_start; drift_grad <= drift_grad_end; drift_grad += drift_grad_stepsize) {
    for(drift_grad_index = 0; drift_grad_index < drift_grad_steps; drift_grad_index++) {
      drift_grad = drift_grad_start + drift_grad_stepsize*drift_grad_index;
      for(N_sparks = nrsparks_start; N_sparks <= nrsparks_end; N_sparks += nrsparks_stepsize) {
        fiducial_bin_index = 0;
        for(fiducial_bin = fid_bin_start; fiducial_bin <= fid_bin_end; fiducial_bin += fid_bin_stepsize) {
          for(alias_order = aliasorder_start; alias_order <= aliasorder_end; alias_order += aliasorder_stepsize) {

            // Update P4
            double k, D;
            k = floor((alias_order + 1.0)/2.0);
            if(alias_order % 2 == 0) {
              D = 1.0;
            }else {
              D = -1.0;
            }
            P4 = (N_sparks*P3) / (k*P3 + D);

            double true_P4;
            // Update true P4
            true_P4 = (N_sparks*true_P3) / (k*true_P3 + D);

            // Update beta
            double Q;
	    // Changed - sign
            Q = -sign_modifier*((1.0/true_P4)+(2.0*M_PI*drift_grad*D)/(N_sparks*true_P3));
            //	    printf("Q = %e\nsign=%e\ntrue_P3=%e\n", Q, sign_modifier, true_P4);
            beta = atan(sin_alpha/(Q-cos_alpha));

            int shift_value = NrBins/2.0 - 1.0 - fiducial_bin;
            int onpulse_bin_start_shifted = long_bin_start + shift_value;
            int onpulse_bin_end_shifted = long_bin_end + shift_value;

            int bin_nr;
            for(bin_nr = 0; bin_nr < NrBins/2; bin_nr++) {
              int opposite_bin_nr;
              opposite_bin_nr = NrBins - 1 - bin_nr;  // opposite bin is trivial
              int sub1, sub2;
              double weight1, weight2, true_opposite;
              // Weights and stuff should not be required in the code below, so code could be simplified
              find_opposite_subint(0, bin_nr, sign_modifier, &sub1, &weight1, &sub2, &weight2, &true_opposite);

              //	      printf("find_opposite_subint: %d %e %d %e %d %e %e\n", bin_nr, sign_modifier, sub1, weight1, sub2, weight2, true_opposite);

              long row_nr;
              // Extract the relevant data to use for the pulse longitude bin considered.
              for(row_nr = 0; row_nr < NrSubints; row_nr++) {

                if(readPulsePSRData(&(rotated_data[fiducial_bin_index]), row_nr, 1, 0, bin_nr, 1, &(data_O1l[row_nr]), application.verbose_state) != 1) {
                  printerror(application.verbose_state.debug, "ERROR: read error");
                  return 0;
                }

                if(readPulsePSRData(&(rotated_data[fiducial_bin_index]), row_nr, 2, 0, bin_nr, 1, &(data_O2l[row_nr]), application.verbose_state) != 1) {
                  printerror(application.verbose_state.debug, "ERROR: read error");
                  return 0;
                }
                if(readPulsePSRData(&(rotated_data[fiducial_bin_index]), row_nr, 1, 0, opposite_bin_nr, 1, &(data_O1t_notrotated[row_nr]), application.verbose_state) != 1) {
                  printerror(application.verbose_state.debug, "ERROR: read error");
                  return 0;
                }
                if(readPulsePSRData(&(rotated_data[fiducial_bin_index]), row_nr, 2, 0, opposite_bin_nr, 1, &(data_O2t_notrotated[row_nr]), application.verbose_state) != 1) {
                  printerror(application.verbose_state.debug, "ERROR: read error");
                  return 0;
                }
              //		printf("data_O1l=%e data_O2l=%e data_O1l=%e data_O2l=%e for bin_nr=%d\n", data_O1l[row_nr], data_O2l[row_nr], data_O1t_notrotated[row_nr], data_O2t_notrotated[row_nr], bin_nr);
              }

              // Certainly can speed this up. In any case linear
              // interpolation might be quicker???? Also if it isn't a
              // power of 2 it might fail??? Otherwise it might be
              // slower as well.
              //	      if(rotateSinglepulse(float *data, int npts, float epsilon, verbose_definition verbose) == 0) {
              //		printerror(application.verbose_state.debug, "ERROR: read error");
              //		return 0;
              //	      }

              // Populate the trailing side array with a rotated version using linear interpolation
              for(row_nr = 0; row_nr < NrSubints; row_nr++) {
                int row_nr_new1 = sub1 + row_nr;
                int row_nr_new2 = sub2 + row_nr;
                if(row_nr_new1 >= NrSubints) {
                  row_nr_new1 -= NrSubints;
                }
                if(row_nr_new2 >= NrSubints) {
                  row_nr_new2 -= NrSubints;
                }
                //		int offending_line;
                data_O1t[row_nr] = (data_O1t_notrotated[row_nr_new1]*weight1 + data_O1t_notrotated[row_nr_new2]*weight2); ///2.0;
                data_O2t[row_nr] = (data_O2t_notrotated[row_nr_new1]*weight1 + data_O2t_notrotated[row_nr_new2]*weight2); ///2.0;
              }

              // Determine guess parameters using the "wrong" linear interpolation
              double guess_a, guess_b, guess_c, guess_d;
              find_initial_guess(&guess_a, &guess_b, &guess_c, &guess_d);

              // Do the fitting, if forcephysical == 0, the two polarizations can be fitted independently
              for(pol_eq_nr = 0; pol_eq_nr < 2; pol_eq_nr++) {

                double xstart[4];
                double dx[4];
                int fixed[4];
                double solution[4];
                int nritt;

		int initial_physical;
		if(forcephysical == 0) {
		  if(pol_eq_nr == 0) {
		    xstart[0] = guess_a;
		    xstart[1] = guess_b;
		  }else {
		    xstart[0] = guess_c;
		    xstart[1] = guess_d;
		  }
		  dx[0] = 0.1;
		  dx[1] = 0.1;
		  fixed[0] = 0;
		  fixed[1] = 0;
		}else {
		  initial_physical = is_physical(guess_a, guess_b, guess_c, guess_d);
		  if(initial_physical) {
		    xstart[0] = guess_a;
		    xstart[1] = guess_b;
		    xstart[2] = guess_c;
		    xstart[3] = guess_d;
		  }else {
		    xstart[0] = 1;
		    xstart[1] = 0;
		    xstart[2] = 0;
		    xstart[3] = 1;
		  }
		  dx[0] = 0.1;
		  dx[1] = 0.1;
		  dx[2] = 0.1;
		  dx[3] = 0.1;
		  fixed[0] = 0;
		  fixed[1] = 0;
		  fixed[2] = 0;
		  fixed[3] = 0;
		}

                double chi2;
		int fit_status;
		//		dumpdebuginfo = 0;
		//		if(opposite_bin_nr == 550 && pol_eq_nr == 1) {
		//		  dumpdebuginfo = 1;
		//		}

		if(forcephysical == 0) {
		  if(optimise_without_denominator == 0) {
		    //		printf("xstart = %e %e\n", xstart[0], xstart[1]);
		    fit_status = doAmoeba_d(0, xstart, dx, fixed, solution, &chi2, 2, &chi2_funk, 1e-6, &nritt, application.verbose_state.verbose, 0, 0.0, NULL, NULL);
		    if(fit_status == 1) {
		      printerror(application.verbose_state.debug, "Error fit_asym: Downhill-Simplex method did not converge. You can try lowering the tolerance with -ftol.");
		      //		  return 0;
		    }else if(fit_status != 0) {
		      printerror(application.verbose_state.debug, "Error fit_asym: Downhill-Simplex method failed.");
		      //		  return 0;
		    }
		  }else {
		    fit_status = 0;
		    solution[0] = xstart[0];
		    solution[1] = xstart[1];
		  }
		}else {   // Fit all 4 parameters at the same time if forcephysical is set
		  if(optimise_without_denominator == 0 || initial_physical == 0) {
		    fit_status = doAmoeba_d(0, xstart, dx, fixed, solution, &chi2, 4, &chi2_funk, 1e-6, &nritt, application.verbose_state.verbose, 0, 0.0, NULL, NULL);
		    if(fit_status == 1) {
		      printerror(application.verbose_state.debug, "Error fit_asym: Downhill-Simplex method did not converge. You can try lowering the tolerance with -ftol.");
		      //		  return 0;
		    }else if(fit_status != 0) {
		      printerror(application.verbose_state.debug, "Error fit_asym: Downhill-Simplex method failed.");
		      //		  return 0;
		    }
		  }else {
		    fit_status = 0;
		    solution[0] = xstart[0];
		    solution[1] = xstart[1];
		    solution[2] = xstart[2];
		    solution[3] = xstart[3];
		  }
		}

		if(optimise_without_denominator) {    // Re-calculate chi2 with dominator included
		  optimise_without_denominator = 0;
		  chi2 = chi2_funk(solution);
		  optimise_without_denominator = 1;
		}

                //		printf("solution = %e %e\n", solution[0], solution[1]);
                //		printf("chi2 = %e\n", chi2);
                //		return 0;

		if(forcephysical == 0) {
		  if(pol_eq_nr == 0) {
		    solution_a[bin_nr] = solution[0];
		    solution_b[bin_nr] = solution[1];
		    chi2_1[bin_nr] = chi2;
		    fit_status_1[bin_nr] = fit_status;
		  }else {
		    solution_c[bin_nr] = solution[0];
		    solution_d[bin_nr] = solution[1];
		    chi2_2[bin_nr] = chi2;
		    fit_status_2[bin_nr] = fit_status;
		  }
		}else {
		  solution_a[bin_nr] = solution[0];
		  solution_b[bin_nr] = solution[1];
		  solution_c[bin_nr] = solution[2];
		  solution_d[bin_nr] = solution[3];
		  chi2_1[bin_nr] = chi2;
		  chi2_2[bin_nr] = chi2;
		  fit_status_1[bin_nr] = fit_status;
		  fit_status_2[bin_nr] = fit_status;
		}


		if(solutiondump) {
		  //int writePulsePSRData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose);
		  int subint;
		  for(subint = 0; subint < NrSubints; subint++) {
		    // Leading half is just the observed polarization
		    if(forcephysical == 0) {
		      float sample;
		      if(pol_eq_nr == 0) {
			sample = data_O1l[subint];
		      }else {
			sample = data_O2l[subint];
		      }
		      if(writePulsePSRData(&output_data_file_solution, subint, pol_eq_nr, 0, bin_nr, 1, &sample, application.verbose_state) != 1) {
			fflush(stdout);
			printerror(application.verbose_state.debug, "ERROR preprocess_removenan: Cannot write data.");
			return 0;
		      }
		    }else {
		      float sample;
		      sample = data_O1l[subint];
		      if(writePulsePSRData(&output_data_file_solution, subint, 0, 0, bin_nr, 1, &sample, application.verbose_state) != 1) {
			fflush(stdout);
			printerror(application.verbose_state.debug, "ERROR preprocess_removenan: Cannot write data.");
			return 0;
		      }
		      sample = data_O2l[subint];
		      if(writePulsePSRData(&output_data_file_solution, subint, 1, 0, bin_nr, 1, &sample, application.verbose_state) != 1) {
			fflush(stdout);
			printerror(application.verbose_state.debug, "ERROR preprocess_removenan: Cannot write data.");
			return 0;
		      }
		    }
		    // Trailing half should be the observed polarization made to match the first half
		    if(forcephysical == 0) {
		      float sample;
		      sample = solution[0]*data_O1t_notrotated[subint]+solution[1]*data_O2t_notrotated[subint];
		      if(writePulsePSRData(&output_data_file_solution, subint, pol_eq_nr, 0, NrBins-bin_nr-1, 1, &sample, application.verbose_state) != 1) {
			fflush(stdout);
			printerror(application.verbose_state.debug, "ERROR preprocess_removenan: Cannot write data.");
			return 0;
		      }
		    }else {
		      float sample;
		      sample = solution[0]*data_O1t_notrotated[subint]+solution[1]*data_O2t_notrotated[subint];
		      if(writePulsePSRData(&output_data_file_solution, subint, 0, 0, NrBins-bin_nr-1, 1, &sample, application.verbose_state) != 1) {
			fflush(stdout);
			printerror(application.verbose_state.debug, "ERROR preprocess_removenan: Cannot write data.");
			return 0;
		      }

		      sample = solution[2]*data_O1t_notrotated[subint]+solution[3]*data_O2t_notrotated[subint];
		      if(writePulsePSRData(&output_data_file_solution, subint, 1, 0, NrBins-bin_nr-1, 1, &sample, application.verbose_state) != 1) {
			fflush(stdout);
			printerror(application.verbose_state.debug, "ERROR preprocess_removenan: Cannot write data.");
			return 0;
		      }
		    }
		  }
		}


                // if(pol_eq_nr == 0) {
                //   solution_a[bin_nr] = guess_a;
                //   solution_b[bin_nr] = guess_b;
                //   chi2_1[bin_nr] = chi2;
                // }else {
                //   solution_c[bin_nr] = guess_c;
                //   solution_d[bin_nr] = guess_d;
                //   chi2_2[bin_nr] = chi2;
                // }
		if(forcephysical) {
		  break;
		}
              } // End of pol_eq_nr loop	   
              physical[bin_nr] = is_physical(solution_a[bin_nr], solution_b[bin_nr], solution_c[bin_nr], solution_d[bin_nr]);
            }  // End of bin loop - Fitting of all longitude bins/polarizations done

            //	    fprintf(stderr, "XXXX %e %e %d %d %d\n", alpha, drift_grad, N_sparks, fiducial_bin, alias_order);
            for(i = 0; i < NrBins; i++) {
	      fprintf(file_matrix, "%e ", solution_a[i]);
            }
            for(i = 0; i < NrBins; i++) {
              fprintf(file_matrix, "%e ", solution_b[i]);
            }
            for(i = 0; i < NrBins; i++) {
              fprintf(file_matrix, "%e ", solution_c[i]);
            }
            for(i = 0; i < NrBins; i++) {
              if(i == NrBins - 1) {
                fprintf(file_matrix, "%e\n", solution_d[i]);
              }else {
                fprintf(file_matrix, "%e ", solution_d[i]);
              }
            }

            fprintf(file_grid, "%f %f %d %d %d %d %d\n", alpha*180.0/M_PI, beta*180.0/M_PI, N_sparks, fiducial_bin, alias_order, onpulse_bin_start_shifted, onpulse_bin_end_shifted);
          
            // ======================= Write out physical file ===============================
            for(i = 0; i < NrBins; i++) {
              if(i == NrBins - 1) {
                fprintf(file_physical, "%d\n", physical[i]);
              }else {
                fprintf(file_physical, "%d ", physical[i]);
              }
            }

            // ======================= Write out chisq file ===============================
            for(i = 0; i < NrBins; i++) {
              fprintf(file_chisq, "%e ", chi2_1[i]);
            }
            for(i = 0; i < NrBins; i++) {
              fprintf(file_chisq, "%e ", chi2_2[i]);
            }
            for(i = 0; i < NrBins; i++) {
              fprintf(file_chisq, "%d ", fit_status_1[i]);
            }
            for(i = 0; i < NrBins; i++) {
              if(i == NrBins - 1) {
                fprintf(file_chisq, "%d\n", fit_status_2[i]);
              }else {
                fprintf(file_chisq, "%d ", fit_status_2[i]);
              }
            }

	    if(profiledump) {
	      for(i = 0; i < NrBins; i++) {
		if(i == NrBins - 1) {
		  fprintf(file_profile, "%e\n", (rotated_profiles[fiducial_bin_index])[i]);
		}else {
		  fprintf(file_profile, "%e ", (rotated_profiles[fiducial_bin_index])[i]);
		}
	      }
	    }

	    if(solutiondump) {
	      // Change the input filename and replaces the extension with something new
	      char outputname[1000];
	      char extension[1000];
	      sprintf(extension, "solutiondump_%05ld", current_itteration);
	      if(change_filename_extension(filename_ptr, outputname, extension, 1000, application.verbose_state) == 0) {
		return 0;
	      }

	      //	      fprintf(stderr, "Going to open %s\n", outputname);
	      // Open an actual file, enabling writing, on the harddisk with the given name and the given format
	      if(!openPSRData(&output_data_file_solution, outputname, FITS_format, 1, 0, 0, application.verbose_state))
		return 0;
	      
	      // Write out the header information
	      if(writeHeaderPSRData(&output_data_file_solution, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
		printerror(application.verbose_state.debug, "ERROR: Unable to write header.\n");
		return 0;
	      }

	      // Writes out all data which was previously already was put in memory to disk
	      if(writePSRData(&output_data_file_solution, output_data_file_solution.data, application.verbose_state) == 0) {
		printerror(application.verbose_state.debug, "ERROR : Unable to write data.\n");
		return 0;
	      }

	      // Close the output file
	      //	      fprintf(stderr, "Going to close %s\n", outputname);
	      closePSRData(&output_data_file_solution, 2, application.verbose_state);
	      // Necessary, so the next writepulse will write to the memory, rather than a closed fits file
	      output_data_file_solution.format = MEMORY_format;
	      //	      fprintf(stderr, "Closing done\n");
	    }


            current_itteration++;
            printf("Processing: %.2lf%%    %ld/%d     %.2f %.2f %d %d  %d  \r", 100.0*(double)current_itteration/(double)(alpha_steps*drift_grad_steps*nrsparks_steps*fid_bin_steps*aliasorder_steps), current_itteration, alpha_steps*drift_grad_steps*nrsparks_steps*fid_bin_steps*aliasorder_steps, alpha*180.0/M_PI, drift_grad, N_sparks, fiducial_bin, alias_order);
            fflush(stdout);
          } // End of alias_order loop (level 5)
          fiducial_bin_index++;
        } // End of fiducial_bin loop (level 4)
      } // End of N_sparks loop (level 3)
    } // End of drift_grad loop (level 2)
  } // End of alpha loop (level 1)
  printf("\n");

  //  fprintf(stderr, "See code at line with offending_line, which shouldn't divide by 2.\n");


  fclose(file_grid);
  fclose(file_matrix);
  fclose(file_physical);
  fclose(file_chisq);
  if(profiledump) {
    fclose(file_profile);
  }

  return 0;
}
