#define TEMPO3_LIB 1   // Needs to be set to something to make tempo3.h to work
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sort.h>
#include "psrsalsa.h"
#include "psrsalsa_tempo3.h"

int GTKinitialised = 0;

/* Some constants */
int TEMPO3_MaxNrWaves =                 25;          /* Maximum number of WAVE parameters supported. Can be changed on command line */
int TEMPO3_MaxNrGlitches =              20;          /* Maximum number of glitches which are supported. Can be changed on command line*/
int TEMPO3_MaxNrJumps =                 10;          /* Maximum number of jumps which are supported. Can be changed on command line */
int TEMPO3_MaxNrCompanions =             2;          // Max nr of companions that can be defined

// The fitting process requires all parameters to be available in a
// vector. The following list of indices identifies where each
// parameter will end up. These indices will automaticly be assingned
// by the tempo3_parfile_descr_add_parameter() function.
// If it is an EPOCH -> it should be EXCLUDED from perturb_parfile
// THESE NEED TO BE IN TEMPO3_LIB.c AND tempo3.h
int TEMPO3_PARAM_PSRJNAME    = -1;  // The JNAME, not a fit parameter, so need to be dealt with differently in places
int TEMPO3_PARAM_POSEPOCH    = -1;  // Epoch (MJD) at which the RAJ and DECJ are defined
int TEMPO3_PARAM_RAJ         = -1;  // The RAJ in radians
int TEMPO3_PARAM_DECJ        = -1;  // The DECJ in radians
int TEMPO3_PARAM_PMRA        = -1;  // The proper motion in RA in milli arcseconds per year (difference with original value)
int TEMPO3_PARAM_PMDEC       = -1;  // The proper motion in DEC in milli arcseconds per year (difference with original value)
int TEMPO3_PARAM_PHASE0      = -1;  // A phase offset (in units phase) applied to all residuals
int TEMPO3_PARAM_DMEPOCH     = -1;  // Epoch (MJD) at which DM and time derivatives are defined
int TEMPO3_PARAM_DM          = -1;  // The actual fit parameter is the DM difference from the actual DM, so there are some non-standard exceptions to be handled. There are TEMPO3_MaxNrDMderivatives time derivatives which follow.
int TEMPO3_PARAM_PEPOCH      = -1;  // Epoch (MJD) at which F0 and higher derivatives are defined
int TEMPO3_PARAM_F0          = -1;  // Points to F0 (rotational frequency in Hz), but it is followed by TEMPO3_MaxNrFderivatives derivatives
int TEMPO3_PARAM_NBRAKE      = -1;  // Braking index fit parameter
int TEMPO3_PARAM_WAVEEPOCH   = -1;
int TEMPO3_PARAM_WAVEOM      = -1;  // Defines the period of the fitwaves, not a fit parameter
int TEMPO3_PARAM_WAVESIN     = -1;  // Points to the first wavesin term, followed by the TEMPO3_MaxNrWaves-1 others
int TEMPO3_PARAM_WAVECOS     = -1;  // Points to the first wavecos term, followed by the TEMPO3_MaxNrWaves-1 others
int TEMPO3_PARAM_GLEP        = -1;  // Glitch epochs for TEMPO3_MaxNrGlitches glitches
int TEMPO3_PARAM_GLPH        = -1;  // Phase step for TEMPO3_MaxNrGlitches glitches
int TEMPO3_PARAM_GLF0        = -1;  // F0 for TEMPO3_MaxNrGlitches glitches, followed by TEMPO3_MaxNrGlitchFderivatives times (TEMPO3_MaxNrGlitches-1) other parameters
int TEMPO3_PARAM_GLTD        = -1;
int TEMPO3_PARAM_GLF0D       = -1;
//int TEMPO3_PARAM_GLF1D       = -1;  // This is not an tempo2 parameter, and it is degenerate with glf0d as is implemented?
int TEMPO3_PARAM_GLNBRAKE    = -1;  // Not yet implemented. Change of braking index after glitch
int TEMPO3_PARAM_GLEPMIN     = -1;  // Min bound (mjd) on glitch epochs for TEMPO3_MaxNrGlitches glitches
int TEMPO3_PARAM_GLEPMAX     = -1;  // Max bound (mjd) on glitch epochs for TEMPO3_MaxNrGlitches glitches
int TEMPO3_PARAM_JUMPS       = -1;  // The jump for different flagged toa's
int TEMPO3_PARAM_VALUE       = -1;  // The VALUE parameters set
int TEMPO3_PARAM_BINARY      = -1;  // The type of binary model to use, not a fit parameter and stored in a different place
int TEMPO3_PARAM_T0          = -1;  // TEMPO3_MaxNrCompanions T0 parameters
int TEMPO3_PARAM_TASC        = -1;  // TEMPO3_MaxNrCompanions parameters
int TEMPO3_PARAM_PB          = -1;  // TEMPO3_MaxNrCompanions parameters
int TEMPO3_PARAM_A1          = -1;  // TEMPO3_MaxNrCompanions parameters
int TEMPO3_PARAM_OM          = -1;  // TEMPO3_MaxNrCompanions parameters
int TEMPO3_PARAM_ECC         = -1;  // TEMPO3_MaxNrCompanions parameters
int TEMPO3_PARAM_EPS1        = -1;  // TEMPO3_MaxNrCompanions parameters
int TEMPO3_PARAM_EPS2        = -1;  // TEMPO3_MaxNrCompanions parameters
int TEMPO3_PARAM_OMDOT       = -1;  // TEMPO3_MaxNrCompanions parameters
int TEMPO3_PARAM_PBDOT       = -1;  // TEMPO3_MaxNrCompanions parameters
int TEMPO3_PARAM_A1DOT       = -1;  // TEMPO3_MaxNrCompanions parameters
int TEMPO3_PARAM_GAMMA       = -1;  // TEMPO3_MaxNrCompanions parameters
int TEMPO3_PARAM_TZRMJD      = -1;
int TEMPO3_PARAM_TZRFRQ      = -1;
int TEMPO3_PARAM_TZRSITE     = -1;  // not a fit parameter, so need to be dealt with differently in places
int TEMPO3_PARAM_MODE        = -1;  // Determines if a weighted fit is used, not a fit parameter and stored in a different place
int TEMPO3_PARAM_TRACK       = -1;  // Determines if pulse numbering is used, not a fit parameter and stored in a different place
int TEMPO3_PARAM_TRES        = -1;  // The RMS of the residual in micro-seconds
int TEMPO3_PARAM_TOTAL       = -1;  // Total nr of parameters to be determined automatically

parfile_description_def tempo3_parfile_descr;   // This holds information about the meaning of all possible parameters definable in a ephemeris
int tempo3_parfile_descr_initialised = 0;

// Defined below, but needs to be defined earlier.
long double evaluate_ephemeris(long double mjd, long double freq, long double ssb1, long double ssb2, long double ssb3, int showWaves, int noNu, char *site, char *flags, tempo3_parameters_def *parfile, long double dphase, int *paramset, int mode, int ignore_binary);


// Should really only be used internally in library, but is declared in tempo3 as well.
//
// If allowerrors is nonzero, space is allocated to store errors on parameters.
//
// Return 0: Error
// Return 1: Success
int allocate_ephemeris_par_only(tempo3_parameters_def *parfile, int allowerrors, verbose_definition verbose)
{
  int i;
  //If changing the allocations done, update the statement about memory per parfile printed out with the -listparameters option
  parfile->parameter = malloc(sizeof(long double)*TEMPO3_PARAM_TOTAL);
  parfile->companion_defined = malloc(sizeof(int)*TEMPO3_MaxNrCompanions);
  parfile->jumpflag = malloc(sizeof(char *)*TEMPO3_MaxNrJumps);
  parfile->jumpmeaning = malloc(sizeof(int)*TEMPO3_MaxNrJumps);
  if(parfile->parameter == NULL || parfile->companion_defined == NULL || parfile->jumpflag == NULL || parfile->jumpmeaning == NULL) {
    printerror(verbose.debug, "ERROR allocate_ephemeris_par_only: Cannot allocate memory for %ld parameters, %d companions and %d jumps\n", TEMPO3_PARAM_TOTAL, TEMPO3_MaxNrCompanions, TEMPO3_MaxNrJumps);
    return 0;
  }
  if(allowerrors) {
    parfile->dplus = malloc(TEMPO3_PARAM_TOTAL*sizeof(long double));
    parfile->dmin = malloc(TEMPO3_PARAM_TOTAL*sizeof(long double));
    if(parfile->dplus == NULL || parfile->dmin == NULL) {
      printerror(verbose.debug, "ERROR allocate_ephemeris_par_only: Cannot allocate memory for errors for %ld parameters, %d companions and %d jumps\n", TEMPO3_PARAM_TOTAL, TEMPO3_MaxNrCompanions, TEMPO3_MaxNrJumps);
      return 0;
    }
  }else {
    parfile->dplus = NULL;
    parfile->dmin = NULL;
  }
  //If changing the allocations done, update the statement about memory per parfile printed out with the -listparameters option
  for(i = 0; i < TEMPO3_MaxNrJumps; i++) {
    parfile->jumpflag[i] = malloc(TEMPO3_MaxFlagsLength);
    if(parfile->jumpflag[i] == NULL) {
      printerror(verbose.debug, "ERROR allocate_ephemeris_par_only: Cannot allocate memory for %d jumps\n", TEMPO3_MaxNrJumps);
      return 0;
    }
  }
  return 1;
}

// Should really only be used internally in library, but is declared in tempo3 as well.
void free_ephemeris_par_only(tempo3_parameters_def *parfile)
{
  int i;
  free(parfile->companion_defined);
  free(parfile->parameter);
  for(i = 0; i < TEMPO3_MaxNrJumps; i++) {
    free(parfile->jumpflag[i]);
  }
  free(parfile->jumpflag);
  free(parfile->jumpmeaning);
  if(parfile->dplus != NULL) {
    free(parfile->dplus);
    parfile->dplus = NULL;
  }
  if(parfile->dmin != NULL) {
    free(parfile->dmin);
    parfile->dmin = NULL;
  }
}

// Should really only be used internally in library, but is declared in tempo3 as well.
void initialise_ephemeris_par_only(tempo3_parameters_def *parfile)
{
  long i;
  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
    parfile->parameter[i] = 0;
  }
  parfile->psrjname[0] = 0;
  for(i = 0; i < TEMPO3_MaxNrDMderivatives + 1; i++) {
    parfile->dm[i] = 0;
  }
  parfile->raj = 0;
  parfile->decj = 0;
  parfile->pmra = 0;
  parfile->pmdec = 0;

  for(i = 0; i < TEMPO3_MaxNrGlitches; i++) {
    parfile->parameter[TEMPO3_PARAM_GLEPMIN+i] = -1;
    parfile->parameter[TEMPO3_PARAM_GLEPMAX+i] = -1;
  }

  parfile->nrglitches = 0;
  parfile->nrwaves = 0;
  for(i = 0; i < TEMPO3_MaxNrJumps; i++) {
    parfile->jumpflag[i][0] = 0;
    parfile->jumpmeaning[i] = 0;
  }
  parfile->nrjumps = 0;

  parfile->binarymodel = 0;
  for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
    parfile->companion_defined[i] = 0;
  }
  parfile->mode = 0;
  parfile->track = 0;
  parfile->tzrsite[0] = 0;
  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
    if(parfile->dplus != NULL) {
      parfile->dplus[i] = 0;
    }
    if(parfile->dmin != NULL) {
      parfile->dmin[i] = 0;
    }
  }
}



// Return 0: error
// Return 1: Success
// If countonly is set, only work out the number of parameters, otherwise assign pointers etc
// print_format pointers are assumed to remain available, but par_ident not.
int tempo3_parfile_descr_add_parameter(int *TEMPO3_PARAM, char *par_ident, char *print_format, long print_format_as_this_param, char *par_unit, long par_unit_as_this_param, char *par_descr, long par_descr_as_this_param, int countonly, char ****parfile_identifiers, char ***parfile_print_format, char ***parfile_units, char ***parfile_description, verbose_definition verbose)
{
  static int cur_paramter = 0;  // Keeps track of how many parameters were added
  int i, j;
  if(countonly) {
    *TEMPO3_PARAM = cur_paramter;         // Assign a fit parameter number
    TEMPO3_PARAM_TOTAL = cur_paramter+1;  // Keep track of how many were assigned already
    cur_paramter++;
  }else {
    // If we're dealing with the first parameter, do some memory allocations
    if(*TEMPO3_PARAM == 0) {
      // Pointer to TEMPO3_PARAM_TOTAL pointers to 2 pointers, giving the parfile names of the parameters (2nd NULL if not used)
      //      printf("XXXX allocating mem for %d pointers\n", TEMPO3_PARAM_TOTAL);
      *parfile_identifiers  = (char ***)malloc(TEMPO3_PARAM_TOTAL*sizeof(char ***));
      *parfile_print_format = (char **)malloc(TEMPO3_PARAM_TOTAL*sizeof(char **));
      *parfile_units        = (char **)malloc(TEMPO3_PARAM_TOTAL*sizeof(char **));
      *parfile_description  = (char **)malloc(TEMPO3_PARAM_TOTAL*sizeof(char **));
      if(*parfile_identifiers == NULL || *parfile_print_format == NULL || *parfile_units == NULL || *parfile_description == NULL) {
	printerror(verbose.debug, "ERROR tempo3_parfile_descr_add_parameter(): Cannot allocate memory\n");
	return 0;
      }
      for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
	(*parfile_identifiers)[i]  = (char **)malloc(2*sizeof(char **));
	//	parfile_print_format[i] = (char *)malloc(sizeof(char *));
	//	parfile_units[i]        = (char *)malloc(sizeof(char *));
	//	parfile_description[i]  = (char *)malloc(sizeof(char *));
	// || parfile_description[i] == NULL  || parfile_units[i] == NULL  || parfile_print_format[i] == NULL
	if((*parfile_identifiers)[i] == NULL) {
	  printerror(verbose.debug, "ERROR tempo3_parfile_descr_add_parameter(): Cannot allocate memory\n");
	  return 0;
	}
	for(j = 0; j < 2; j++) {
	  (*parfile_identifiers)[i][j] = NULL;
	}
	(*parfile_print_format)[i] = NULL;
	(*parfile_units)[i] = NULL;
	(*parfile_description)[i] = NULL;
      }
    }

    (*parfile_identifiers)[*TEMPO3_PARAM][0] = malloc(strlen(par_ident)+1);
    if((*parfile_identifiers)[*TEMPO3_PARAM][0] == NULL) {
      printerror(verbose.debug, "ERROR tempo3_parfile_descr_add_parameter(): Cannot allocate memory\n");
      return 0;
    }
    strcpy((*parfile_identifiers)[*TEMPO3_PARAM][0], par_ident);
    (*parfile_print_format)[*TEMPO3_PARAM] = print_format;

    if(par_unit_as_this_param >= 0) {
      (*parfile_units)[*TEMPO3_PARAM] = (*parfile_units)[par_unit_as_this_param];
    }else {
      (*parfile_units)[*TEMPO3_PARAM] = malloc(strlen(par_unit)+1);
      if((*parfile_units)[*TEMPO3_PARAM] == NULL) {
	printerror(verbose.debug, "ERROR tempo3_parfile_descr_add_parameter(): Cannot allocate memory\n");
	return 0;
      }
      strcpy((*parfile_units)[*TEMPO3_PARAM], par_unit);
    }

    if(par_descr_as_this_param >= 0) {
      (*parfile_description)[*TEMPO3_PARAM] = (*parfile_description)[par_descr_as_this_param];
    }else {
      (*parfile_description)[*TEMPO3_PARAM] = malloc(strlen(par_descr)+1);
      if((*parfile_description)[*TEMPO3_PARAM] == NULL) {
	printerror(verbose.debug, "ERROR tempo3 - tempo3_parfile_descr_add_parameter(): Cannot allocate memory\n");
	return 0;
      }
      strcpy((*parfile_description)[*TEMPO3_PARAM], par_descr);
    }
  }
  return 1;
}



// Should be called as soon as TEMPO3_MaxNrWaves,
// TEMPO3_MaxNrGlitches, TEMPO3_MaxNrJumps, TEMPO3_MaxNrCompanions are
// defined (or immediately if not changing these values). This
// initialises some global variables used in the library, such as
// descriptions of the different ephemeris parameters and their units
// which are used when parameters are reported. The associated memory
// will be released after calling cleanup_tempo3_lib(), which should
// be done when the program is finished using the tempo3 library.
// 
// Return 0: Error
// Return 1: Success
int initialise_tempo3_lib(verbose_definition verbose)
{
  int max_description_string_length = 1000;
  char *txt, *txt2, *txt3;
  int i, j, k, countonly;

  if(tempo3_parfile_descr_initialised) {
    printerror(verbose.verbose, "ERROR initialise_tempo3_lib: The library already appears to be initialised.");
    return 0;
  }

  txt = malloc(max_description_string_length+1);
  txt2 = malloc(max_description_string_length+1);
  txt3 = malloc(max_description_string_length+1);
  if(txt == NULL || txt2 == NULL ||txt3 == NULL) {
    printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Memory allocation error.");
    return 0;
  }

  // Note: if order changes, the order of the steps in the fitting process might be different, leading to (slightly) different results.
  for(countonly = 1; countonly >= 0; countonly--) {
    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_PSRJNAME, "PSRJ",   "%s", -1,        "", -1, "J2000 name of the pulsar", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }

    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_POSEPOCH, "POSEPOCH", "%-30.23Lf", -1, "MJD", -1, "Epoch at which RAJ/DECJ are defined", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_RAJ,      "RAJ",      "%-34s", -1,        "h:m:s", -1, "Right ascension", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_DECJ,     "DECJ",     "%-34s", -1,        "d:m:s", -1, "Declination", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_PMRA,     "PMRA",     "%-30.23Lf", -1,    "mas/yr", -1, "Proper motion in RA", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_PMDEC,    "PMDEC",    "%-30.23Lf", -1,    "mas/yr", -1, "Proper motion in DEC", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_DMEPOCH,       "DMEPOCH",     "%-30.23Lf", -1, "MJD", -1, "Epoch at which DM1 etc are defined", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_DM,       "DM",       "%-30.23Lf", -1, "cm^-3 pc", -1, "Dispersion measure", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    // Time derivatives of DM
    for(i = 0; i < TEMPO3_MaxNrDMderivatives; i++) {
      sprintf(txt, "DM%d", i+1);
      j = TEMPO3_PARAM_DM + i+1;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      sprintf(txt2, "cm^-3 pc yr^%d", -(i+1));
      if(i == 0) {
	sprintf(txt3, "The %dst DM derivative", i+1);
      }else if(i == 1) {
	sprintf(txt3, "The %dnd DM derivative", i+1);
      }else {
	sprintf(txt3, "The %dth DM derivative", i+1);
      }
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Le", -1, txt2, -1, txt3, -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
	printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
	return 0;
      }
      TEMPO3_PARAM_DM = j - i -1;      // Only store postition of DM itself, others follow this parameter
    }

    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_PHASE0,   "PHASE0",   "%-30.23Le", -1, "", -1, "Phase offset applied to all residuals to minimise goodness of fit", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }

    for(i = 0; i < TEMPO3_MaxNrValues; i++) {
      sprintf(txt, "VALUE%d", i+1);
      j = TEMPO3_PARAM_VALUE + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Le", -1, "UNDEFINED", -1, "User supplied value that can be assigned to other parameters (i.e. \"GLF0_1 VALUE1\" and \"GLF0_1 -VALUE1\")", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_VALUE = j - i;   // Only store postition of VALUE1, others follow this parameter
    }

    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_PEPOCH,       "PEPOCH",     "%-30.23Lf", -1, "MJD", -1, "Epoch at which F0 etc are defined", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    for(i = 0; i <= TEMPO3_MaxNrFderivatives; i++) {
      sprintf(txt, "F%d", i);
      j = TEMPO3_PARAM_F0 + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      sprintf(txt2, "s^%d", -1-i);
      if(i == 0) {
	sprintf(txt3, "Spin frequency");
      }else if(i == 1) {
	sprintf(txt3, "The %dst spin frequency derivative", i);
      }else {
	sprintf(txt3, "The %dth spin frequency derivative", i);
      }
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Le", -1, txt2, -1, txt3, -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_F0 = j - i;      // Only store postition of F0, others follow this parameter
    }
    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_NBRAKE,   "NBRAKE", "%-30.23Lf", -1, "", -1, "Braking index", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }

    long par_descr_as_this_param = -1;
    long par_unit_as_this_param = -1;
    for(i = 0; i < TEMPO3_MaxNrGlitches; i++) {
      sprintf(txt, "GLEP_%d", i+1);
      j = TEMPO3_PARAM_GLEP + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      if(i == 0 || i == TEMPO3_MaxNrGlitches-1) {
	sprintf(txt3, "Epoch at which GLF0_%d etc are defined", i+1);
	par_descr_as_this_param = -1;
      }else if(i == 1) {
	sprintf(txt3, "Similar to GLF0_1");
	par_unit_as_this_param = TEMPO3_PARAM_GLEP;
      }else {
	par_descr_as_this_param = TEMPO3_PARAM_GLEP + 1;
      }
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Lf", par_unit_as_this_param, "MJD", par_unit_as_this_param, txt3, par_descr_as_this_param, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_GLEP = j - i;   // Only store postition of GLEP_1, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrGlitches; i++) {
      sprintf(txt, "GLEPMIN_%d", i+1);
      j = TEMPO3_PARAM_GLEPMIN + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      if(i == 0) {
	sprintf(txt3, "Lower epoch bound allowed for GLEP_%d fitting", i+1);
	par_descr_as_this_param = -1;
	par_unit_as_this_param = -1;
      }else if(i == 1) {
	sprintf(txt3, "Similar to GLEPMIN_1");
	par_unit_as_this_param = TEMPO3_PARAM_GLEPMIN;
      }else {
	par_descr_as_this_param = TEMPO3_PARAM_GLEPMIN + 1;
      }
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Lf", par_unit_as_this_param, "MJD", par_unit_as_this_param, txt3, par_descr_as_this_param, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_GLEPMIN = j - i;   // Only store postition of GLEPMIN_1, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrGlitches; i++) {
      sprintf(txt, "GLEPMAX_%d", i+1);
      j = TEMPO3_PARAM_GLEPMAX + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      if(i == 0) {
	sprintf(txt3, "Upped epoch bound allowed for GLEP_%d fitting", i+1);
	par_descr_as_this_param = -1;
	par_unit_as_this_param = -1;
      }else if(i == 1) {
	sprintf(txt3, "Similar to GLEPMAX_1");
	par_unit_as_this_param = TEMPO3_PARAM_GLEPMAX;
      }else {
	par_descr_as_this_param = TEMPO3_PARAM_GLEPMAX + 1;
      }
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Lf", par_unit_as_this_param, "MJD", par_unit_as_this_param, txt3, par_descr_as_this_param, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_GLEPMAX = j - i;   // Only store postition of GLEPMAX_1, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrGlitches; i++) {
      sprintf(txt, "GLPH_%d", i+1);
      j = TEMPO3_PARAM_GLPH + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      if(i == 0) {
	sprintf(txt3, "Jump in phase at glitch %d", i+1);
	par_descr_as_this_param = -1;
	par_unit_as_this_param = -1;
      }else if(i == 1) {
	sprintf(txt3, "Similar to GLPH_1");
	par_unit_as_this_param = TEMPO3_PARAM_GLPH;
      }else {
	par_descr_as_this_param = TEMPO3_PARAM_GLPH + 1;
      }
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Le", par_unit_as_this_param, "", par_unit_as_this_param, txt3, par_descr_as_this_param, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_GLPH = j - i;   // Only store postition of GLPH_1, others follow this parameter
    }
    for(j = 0; j <= TEMPO3_MaxNrGlitchFderivatives; j++) {
      for(i = 0; i < TEMPO3_MaxNrGlitches; i++) {
	sprintf(txt, "GLF%d_%d", j, i+1);
	k = TEMPO3_PARAM_GLF0 + j*TEMPO3_MaxNrGlitches + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
	sprintf(txt2, "s^%d", -1-j);
	if(i == 0 || (i == TEMPO3_MaxNrGlitches-1 && j == TEMPO3_MaxNrGlitchFderivatives)) {
	  sprintf(txt3, "Jump in F%d at glitch %d", j, i+1);
	  par_descr_as_this_param = -1;
	  par_unit_as_this_param = -1;
	}else if(i == 1) {
	  sprintf(txt3, "Similar to GLF%d_1", j);
	  par_unit_as_this_param = TEMPO3_PARAM_GLF0 + j*TEMPO3_MaxNrGlitches;
	}else {
	  par_descr_as_this_param = TEMPO3_PARAM_GLF0 + j*TEMPO3_MaxNrGlitches + 1;
	}
	if(tempo3_parfile_descr_add_parameter(&k,       txt,   "%-30.23Le", par_unit_as_this_param, txt2, par_unit_as_this_param, txt3, par_descr_as_this_param, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
	TEMPO3_PARAM_GLF0 = k - (j*TEMPO3_MaxNrGlitches + i);   // Only store postition of GLF0_1, other derivatives/glitches follow this parameter
      }
    }
    for(i = 0; i < TEMPO3_MaxNrGlitches; i++) {
      sprintf(txt, "GLNBRAKE_%d", i+1);
      j = TEMPO3_PARAM_GLNBRAKE + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      if(i == 0) {
	sprintf(txt3, "Jump in the braking index at glitch %d", i+1);
	par_descr_as_this_param = -1;
	par_unit_as_this_param = -1;
      }else if(i == 1) {
	sprintf(txt3, "Similar to GLNBRAKE_1");
	par_unit_as_this_param = TEMPO3_PARAM_GLNBRAKE;
      }else {
	par_descr_as_this_param = TEMPO3_PARAM_GLNBRAKE + 1;
      }
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Lf", par_unit_as_this_param, "", par_unit_as_this_param, txt3, par_descr_as_this_param, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_GLNBRAKE = j - i;   // Only store postition of GLNBRAKE_1, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrGlitches; i++) {
      sprintf(txt, "GLTD_%d", i+1);
      j = TEMPO3_PARAM_GLTD + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      if(i == 0) {
	sprintf(txt3, "Timescale of exponential glitch recovery of glitch %d", i+1);
	par_descr_as_this_param = -1;
	par_unit_as_this_param = -1;
      }else if(i == 1) {
	sprintf(txt3, "Similar to GLTD_1");
	par_unit_as_this_param = TEMPO3_PARAM_GLTD;
      }else {
	par_descr_as_this_param = TEMPO3_PARAM_GLTD + 1;
      }
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Lf", par_unit_as_this_param, "days", par_unit_as_this_param, txt3, par_descr_as_this_param, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_GLTD = j - i;   // Only store postition of GLTD_1, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrGlitches; i++) {
      sprintf(txt, "GLF0D_%d", i+1);
      j = TEMPO3_PARAM_GLF0D + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      if(i == 0) {
	sprintf(txt3, "Jump in F0 which is exponentially recovering for glitch %d", i+1);
	par_descr_as_this_param = -1;
	par_unit_as_this_param = -1;
      }else if(i == 1) {
	sprintf(txt3, "Similar to GLF0D_1");
	par_unit_as_this_param = TEMPO3_PARAM_GLF0D;
      }else {
	par_descr_as_this_param = TEMPO3_PARAM_GLF0D + 1;
      }
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Le", par_unit_as_this_param, "Hz", par_unit_as_this_param, txt3, par_descr_as_this_param, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_GLF0D = j - i;   // Only store postition of GLF0D_1, others follow this parameter
    }
    /*
    for(i = 0; i < TEMPO3_MaxNrGlitches; i++) {
      sprintf(txt, "GLF1D_%d", i+1);
      j = TEMPO3_PARAM_GLF1D + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Le", "", "", countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_GLF1D = j - i;   // Only store postition of GLF1D_1, others follow this parameter
    }
    */

    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_BINARY,       "BINARY",     "%NOTUSED", -1,  "", -1, "Name of binary model to use", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
      if(i == 0)
	sprintf(txt, "T0");
      else
	sprintf(txt, "T0_%d", i+1);
      j = TEMPO3_PARAM_T0 + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      sprintf(txt3, "Epoch of periastron passage of companion %d (or ascending node passage for circular orbits with OM paramter set to zero)", i+1);
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Lf", -1, "MJD", -1, txt3, -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_T0 = j - i;   // Only store postition of T0, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
      if(i == 0)
	sprintf(txt, "TASC");
      else
	sprintf(txt, "TASC_%d", i+1);
      j = TEMPO3_PARAM_TASC + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      sprintf(txt3, "Epoch of ascending node passage of companion %d (can be used in T2 model)", i+1);
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Lf", -1, "MJD", -1, txt3, -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_TASC = j - i;   // Only store postition of TASC, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
      if(i == 0)
	sprintf(txt, "PB");
      else
	sprintf(txt, "PB_%d", i+1);
      j = TEMPO3_PARAM_PB + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      sprintf(txt3, "Orbital period of companion %d", i+1);
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Lf", -1, "days", -1, txt3, -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_PB = j - i;   // Only store postition of PB, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
      if(i == 0)
	sprintf(txt, "PBDOT");
      else
	sprintf(txt, "PBDOT_%d", i+1);
      j = TEMPO3_PARAM_PBDOT + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      sprintf(txt3, "Orbital period decay of companion %d", i+1);
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Le", -1, "s/s", -1, txt3, -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_PBDOT = j - i;   // Only store postition of PBDOT, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
      if(i == 0)
	sprintf(txt, "A1");
      else
	sprintf(txt, "A1_%d", i+1);
      j = TEMPO3_PARAM_A1 + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      sprintf(txt3, "Projected semimajor axis of orbit of companion %d", i+1);
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Lf", -1, "lt-s", -1, txt3, -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_A1 = j - i;   // Only store postition of A1, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
      if(i == 0)
	sprintf(txt, "A1DOT");
      else
	sprintf(txt, "A1DOT_%d", i+1);
      j = TEMPO3_PARAM_A1DOT + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      sprintf(txt3, "Rate of change of semimajor axisof companion %d", i+1);
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Le", -1, "lt-s/s", -1, txt3, -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_A1DOT = j - i;   // Only store postition of A1DOT, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
      if(i == 0)
	sprintf(txt, "OM");
      else
	sprintf(txt, "OM_%d", i+1);
      j = TEMPO3_PARAM_OM + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      sprintf(txt3, "Longitude of periastron of companion %d", i+1);
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Lf", -1, "deg", -1, txt3, -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_OM = j - i;   // Only store postition of OM, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
      if(i == 0)
	sprintf(txt, "OMDOT");
      else
	sprintf(txt, "OMDOT_%d", i+1);
      j = TEMPO3_PARAM_OMDOT + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      sprintf(txt3, "Periastron advance of companion %d", i+1);
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Le", -1, "deg/yr", -1, txt3, -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_OMDOT = j - i;   // Only store postition of OMDOT, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
      if(i == 0)
	sprintf(txt, "ECC");
      else
	sprintf(txt, "ECC_%d", i+1);
      j = TEMPO3_PARAM_ECC + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      sprintf(txt3, "Orbital eccentricity of companion %d", i+1);
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Lf", -1, "", -1, txt3, -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_ECC = j - i;   // Only store postition of ECC, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
      if(i == 0)
	sprintf(txt, "EPS1");
      else
	sprintf(txt, "EPS1_%d", i+1);
      j = TEMPO3_PARAM_EPS1 + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      sprintf(txt3, "ECC x sin(OM) for ELL1/T2 model of companion %d", i+1);
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Lf", -1, "", -1, txt3, -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_EPS1 = j - i;   // Only store postition of EPS1, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
      if(i == 0)
	sprintf(txt, "EPS2");
      else
	sprintf(txt, "EPS2_%d", i+1);
      j = TEMPO3_PARAM_EPS2 + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      sprintf(txt3, "ECC x cos(OM) for ELL1/T2 model of companion %d", i+1);
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Lf", -1, "", -1, txt3, -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_EPS2 = j - i;   // Only store postition of EPS2, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
      if(i == 0)
	sprintf(txt, "GAMMA");
      else
	sprintf(txt, "GAMMA_%d", i+1);
      j = TEMPO3_PARAM_GAMMA + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      sprintf(txt3, "Post-Keplerian \"gamma\" term of companion %d", i+1);
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Le", -1, "s", -1, txt3, -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_GAMMA = j - i;   // Only store postition of GAMMA, others follow this parameter
    }


    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_WAVEEPOCH,    "WAVEEPOCH", "%-30.23Lf", -1, "MJD", -1, "Epoch at which WAVE parameters are defined", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_WAVEOM,   "WAVE_OM", "%-30.23Lf", -1, "rad/day", -1, "Angular frequency of first WAVE parameter", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    for(i = 0; i < TEMPO3_MaxNrWaves; i++) {
      sprintf(txt, "WAVE%d", i+1);
      j = TEMPO3_PARAM_WAVESIN + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      sprintf(txt3, "Sine and cosine amplitudes of harmonic %d", i+1);
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Le", -1, "s", -1, txt3, -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_WAVESIN = j - i;   // Only store postition of WAVE1, others follow this parameter
    }
    for(i = 0; i < TEMPO3_MaxNrWaves; i++) {
      sprintf(txt, "WAVE%d", i+1);
      j = TEMPO3_PARAM_WAVECOS + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Le", -1, "", -1, "", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_WAVECOS = j - i;   // Only store postition of WAVE1, others follow this parameter
    }


    for(i = 0; i < TEMPO3_MaxNrJumps; i++) {
      sprintf(txt, "JUMP");
      j = TEMPO3_PARAM_JUMPS + i;   // tempo3_parfile_descr_add_parameter needs to know what parameter you talking about, if countonly != 0
      if(tempo3_parfile_descr_add_parameter(&j,       txt,   "%-30.23Lf", -1, "s", -1, "Example: JUMP -t JB XXXX adds XXXX seconds to all tao's with flag -t JB", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
      TEMPO3_PARAM_JUMPS = j - i;   // Only store postition of the first JUMP, others follow this parameter
    }

    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_TZRMJD,     "TZRMJD",   "%-30.23Lf", -1, "MJD", -1, "TOA corresponding to phase 0 definition.", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_TZRFRQ,     "TZRFRQ",   "%-30.23Lf", -1, "MHz", -1, "Freq of TOA defining phase 0", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_TZRSITE,    "TZRSITE",  "%s", -1,        "", -1, "Site code of TOA defining phase 0", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_MODE,       "MODE",     "%d", -1,        "", -1, "Weighted fit (=1) or non-weighted fit (=0)", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_TRACK,      "TRACK",    "%d", -1,        "", -1, "Enable pulse numbering (=-2) or not (=0)", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }
    if(tempo3_parfile_descr_add_parameter(&TEMPO3_PARAM_TRES,       "TRES",     "%-30.23Lf", -1, "micro-s", -1, "RMS of residuals", -1, countonly, &(tempo3_parfile_descr.identifiers), &(tempo3_parfile_descr.print_format), &(tempo3_parfile_descr.units), &(tempo3_parfile_descr.description), verbose) == 0) {
      printerror(verbose.verbose, "ERROR initialise_tempo3_lib: Adding parameter description failed.");
      return 0;
    }

    //    if(HARDCODED_NR_PARAMETERS != TEMPO3_PARAM_TOTAL) {
    //      printerror(0, "ERROR tempo3: HARDCODED_NR_PARAMETERS (%d) set incorrectly. The number of supported parameters is currently %d. I know it is very ugly coding, but please change this number and re-compile.", HARDCODED_NR_PARAMETERS, TEMPO3_PARAM_TOTAL);
    //      return 0;
    //    }
  }
  tempo3_parfile_descr_initialised = 1;
  free(txt);
  free(txt2);
  free(txt3);
  return 1;
}


// The memory used after calling initialise_tempo3_lib() will be
// released, which should be done when the program is finished using
// the tempo3 library.
//
// Return 0: Error
// Return 1: Ok
int cleanup_tempo3_lib(verbose_definition verbose)
{
  long i, j;
  void *ptr;

  if(tempo3_parfile_descr_initialised == 0) {
    printerror(verbose.verbose, "ERROR cleanup_tempo3_lib: The library does not appear to be initialised.");
    return 0;
  }

  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
    if(tempo3_parfile_descr.identifiers[i][0] != NULL)
      free(tempo3_parfile_descr.identifiers[i][0]);
    if(tempo3_parfile_descr.identifiers[i][1] != NULL)
      free(tempo3_parfile_descr.identifiers[i][1]);
    free(tempo3_parfile_descr.identifiers[i]);
    //    if(tempo3_parfile_descr.print_format[i] != NULL)
    //      free(tempo3_parfile_descr.print_format[i]);
    if(tempo3_parfile_descr.description[i] != NULL ) {
      ptr = tempo3_parfile_descr.description[i];
      free(tempo3_parfile_descr.description[i]);
      for(j = 0; j < TEMPO3_PARAM_TOTAL; j++) {  // Same pointer can be used multiple times to reduce memory use
	if(tempo3_parfile_descr.description[j] == ptr) {
	  tempo3_parfile_descr.description[j] = NULL;
	}
      }
    }
    if(tempo3_parfile_descr.units[i] != NULL ) {
      ptr = tempo3_parfile_descr.units[i];
      free(tempo3_parfile_descr.units[i]);
      for(j = 0; j < TEMPO3_PARAM_TOTAL; j++) {  // Same pointer can be used multiple times to reduce memory use
	if(tempo3_parfile_descr.units[j] == ptr) {
	  tempo3_parfile_descr.units[j] = NULL;
	}
      }
    }
  }
  free(tempo3_parfile_descr.identifiers);
  free(tempo3_parfile_descr.print_format);
  free(tempo3_parfile_descr.units);
  free(tempo3_parfile_descr.description);

  tempo3_parfile_descr_initialised = 0;
  return 1;
}

// If VALUE parameters are used, insert the value in the parameter for evaluation of the ephemeris.
// Should really only be used internally in library, but is declared in tempo3 as well.
void ephemeris_setLinkedValues_par_version(tempo3_parameters_def *parfile, int *paramlinked)
{
  int i, param, value_id, nrentries;
  long double value;

  if(paramlinked == NULL) {
    return;
  }

  nrentries = *paramlinked;               // First entry is the number of entries
  for(i = 0; i < nrentries; i++) {
    paramlinked++;
    param = *paramlinked;                     // The parameter number that needs to be replaced
    paramlinked++;
    value_id = *paramlinked;                  // The value ID to be used (counting from 1)
    if(value_id > 0) {
      value = parfile->parameter[TEMPO3_PARAM_VALUE + value_id - 1];
    }else {
      value = -parfile->parameter[TEMPO3_PARAM_VALUE - value_id - 1];
    }
    parfile->parameter[param] = value;
  }
}

// The ephemeris might contain VALUE parameters. The call to this
// function copies their values to the appropriate parameters such
// that the ephemeris can be evaluated.
void ephemeris_setLinkedValues(ephemeris_def *eph)
{
  ephemeris_setLinkedValues_par_version(&(eph->par), eph->paramlinked);
}

// Should really only be used internally in library, but is declared in tempo3 as well.
//
// If VALUE parameters are used, find out if, and which VALUE parameter is associated with a given ephemeris variable.
//
// Return:  0 = no
// Return: >0 = id (counting from 1, or neg if neg value is used).
int ephemeris_check_if_parameter_is_linked(int param_id, int *paramlinked)
{
  int i, param, nrentries, value_id;

  nrentries = *paramlinked;               // First entry is the number of entries
  for(i = 0; i < nrentries; i++) {
    paramlinked++;
    param = *paramlinked;                     // The parameter number that needs to be replaced
    paramlinked++;
    value_id = *paramlinked;                  // The value ID to be used (counting from 1)
    if(param == param_id)
      return value_id;
  }
  return 0;
}


/* Output the ephemeris eph to stream. A call to
  ephemeris_setLinkedValues() will be made.  If eph2 != NULL, then the
  difference (eph-eph2) is shown as well. No checks are done if the
  same parameters are set in eph and eph2. If toas != NULL, then TOA
  related information is given, such as the dataspan used and the
  number of TOAs. It also allows a warning to be generated if
  position-related (or DM) fitting was based on already barycentred
  toa's. If ignorewraps is nonzero, the wrap parameters are not
  shown. If ignoreunrecognizedlines is nonzero, the lines from the
  ephemeris which were read in, but not recognized/used, will not be
  shown. If ignorefitted is nonzero, the information about if the
  parameter is fitted for is omitted from the output. If ignoreerrors
  is nonzero, the information about if the fit errors is omitted from
  the output. autoquiet is used in tempo3 (with the -auto modes), but
  set it to zero otherwise. If hex is nonzero, parameters are written
  as hex encoded bytes values, which avoids small changes in
  parameters before/after writing out a parameter file. The indent
  parameter of verbose can be used to set the number of spaces before
  the output. The use of this function requires the following function
  to be defined by the user:

  // Not used function, but needs to be defined to use library
  void print_ephemeris_gtk()
  {
  }

  Return 0: Error
  Return 1: Success
*/
// Make gotBarycentreTOAs non-zero if the fitting was based on already barycentred toa's (obsolete, so delete this line). 
int print_ephemeris(FILE *stream, ephemeris_def *eph, ephemeris_def *eph2, toa_data_def *toas, int ignorewraps, int ignoreunrecognizedlines, int ignorefitted, int ignoreerrors, int autoquiet, int hex, verbose_definition verbose)
{
  char txt[100], textBuffer[1000];
  int i, curparam; //, linenr
  tempo3_parameters_def *parfile = &(eph->par);
  tempo3_parameters_def *parfile2 = &(eph2->par);

  // Insert VALUE values in right parameter
  ephemeris_setLinkedValues(eph);

  //  linenr = 0;
  txt[0] = 0;
  for(i = 0; i < verbose.indent; i++)
    strcat(txt, " ");

  if(autoquiet) {   // This bit is a special output mode? Normal output is generated below.
    if(ignoreerrors == 0) {  // Only output stuff when errors are determined, otherwise it is an intermediate result
      if(eph->paramset[TEMPO3_PARAM_PEPOCH]) {
	fprintf(stream, "%Lf ", parfile->parameter[TEMPO3_PARAM_PEPOCH]);	
      }
      for(curparam = 0; curparam < TEMPO3_PARAM_TOTAL; curparam++) {
	if(eph->fixed[curparam] == 0 && curparam != TEMPO3_PARAM_PHASE0) { // Show values + errors of all fitted parameters
	  fprintf(stream, "%.10Le ", parfile->parameter[curparam]);
	  fprintf(stream, "%.4Le ", parfile->dplus[curparam]);
	}
      }
      if(eph->paramset[TEMPO3_PARAM_TRES]) {
	fprintf(stream, "%.4Le ", parfile->parameter[TEMPO3_PARAM_TRES]);	
      }
      fprintf(stream, "\n");
    }
    return 1;
  }

  if(toas != NULL && stream != stdout) {
    if(toas->gotBarycentreTOAs) {  // Par file is being written out and the TOA's in TIM file are defined at the SSB
      // Check if position/DM was fitted for
      if(
	 (eph->paramset[TEMPO3_PARAM_DM] && parfile->parameter[TEMPO3_PARAM_DM] != 0.0)
	 || (eph->paramset[TEMPO3_PARAM_RAJ] && parfile->parameter[TEMPO3_PARAM_RAJ] != 0.0)
	 || (eph->paramset[TEMPO3_PARAM_DECJ] && parfile->parameter[TEMPO3_PARAM_DECJ] != 0.0)
	 || (eph->paramset[TEMPO3_PARAM_PMRA] && parfile->parameter[TEMPO3_PARAM_PMRA] != 0.0)
	 || (eph->paramset[TEMPO3_PARAM_PMDEC] && parfile->parameter[TEMPO3_PARAM_PMDEC] != 0.0)
	 ) {
	printwarning(verbose.debug, "WARNING: You have fitted for DM and/or position related parameters. Since the input toa's appear to be already barycentered, the ephemeris that is going to be written out will be incompatible with the toa's that were used.");
      }
    }
  }


  for(curparam = 0; curparam < TEMPO3_PARAM_TOTAL; curparam++) {
    if(eph->paramset[curparam] &&
       (curparam < TEMPO3_PARAM_WAVECOS || curparam >= TEMPO3_PARAM_WAVECOS + TEMPO3_MaxNrWaves)  // Ignore cos terms, as they were already dealt with when outputting the sin terms
       ) {
      int value_id = ephemeris_check_if_parameter_is_linked(curparam, eph->paramlinked);
      if(curparam >= TEMPO3_PARAM_JUMPS && curparam < TEMPO3_PARAM_JUMPS + TEMPO3_MaxNrJumps) {
	int maxlengthflag = 0;
	for(i = TEMPO3_PARAM_JUMPS; i < TEMPO3_PARAM_JUMPS + TEMPO3_MaxNrJumps; i++) {
	  if(eph->paramset[i]) {
	    if(strlen(parfile->jumpflag[i-TEMPO3_PARAM_JUMPS]) > maxlengthflag) {
	      maxlengthflag = strlen(parfile->jumpflag[i-TEMPO3_PARAM_JUMPS]);
	      if(parfile->jumpmeaning[i-TEMPO3_PARAM_JUMPS] == 1) { // TEL XXXXX option
		maxlengthflag += 4;
	      }
	    }
	  }
	}
	if(parfile->jumpmeaning[curparam-TEMPO3_PARAM_JUMPS] == 0) {
	  fprintf(stream, "%s%-9s %s", txt, tempo3_parfile_descr.identifiers[curparam][0], parfile->jumpflag[curparam-TEMPO3_PARAM_JUMPS]);
	}else if(parfile->jumpmeaning[curparam-TEMPO3_PARAM_JUMPS] == 1) {
	  fprintf(stream, "%s%-9s TEL %s", txt, tempo3_parfile_descr.identifiers[curparam][0], parfile->jumpflag[curparam-TEMPO3_PARAM_JUMPS]);
	}else {
	  printerror(verbose.debug, "ERROR print_ephemeris: jumpmeaning=%d is not defined. This is a bug.", parfile->jumpmeaning[i-TEMPO3_PARAM_JUMPS]);
	  return 0;
	}
	for(i = 0; i < maxlengthflag - strlen(parfile->jumpflag[curparam-TEMPO3_PARAM_JUMPS]) + 2; i++)
	  fprintf(stream, " ");
      }else {
	fprintf(stream, "%s%-9s  ", txt, tempo3_parfile_descr.identifiers[curparam][0]);
      }
      if(curparam == TEMPO3_PARAM_PSRJNAME) {
	fprintf(stream, " ");
	fprintf(stream, tempo3_parfile_descr.print_format[curparam], parfile->psrjname);
      }else if(curparam == TEMPO3_PARAM_RAJ) {
	fprintf(stream, " ");
	converthms_string(textBuffer, (parfile->raj+parfile->parameter[curparam])*12.0/M_PI, 20, 1);
	fprintf(stream, tempo3_parfile_descr.print_format[curparam], textBuffer);
      }else if(curparam == TEMPO3_PARAM_DECJ) {
	fprintf(stream, " ");
	converthms_string(textBuffer, (parfile->decj+parfile->parameter[curparam])*180.0/M_PI, 20, 1);
	fprintf(stream, tempo3_parfile_descr.print_format[curparam], textBuffer);
      }else if(curparam == TEMPO3_PARAM_TZRSITE) {
	fprintf(stream, " "); // Print a space instead of a minus sign for better looking alignment
	fprintf(stream, tempo3_parfile_descr.print_format[curparam], parfile->tzrsite);
      }else if(curparam == TEMPO3_PARAM_BINARY) {
	fprintf(stream, " "); // Print a space instead of a minus sign for better looking alignment
	if(parfile->binarymodel == 1) {
	  fprintf(stream, "BT");
	}else if(parfile->binarymodel == 2) {
	  fprintf(stream, "T2");
	}else if(parfile->binarymodel == 3) {
	  fprintf(stream, "HYPER");
	}else {
	  fprintf(stream, "NOT SUPPORTED");
	}
      }else if(curparam == TEMPO3_PARAM_MODE) {
	fprintf(stream, " "); // Print a space instead of a minus sign for better looking alignment
	fprintf(stream, tempo3_parfile_descr.print_format[curparam], parfile->mode);
      }else if(curparam == TEMPO3_PARAM_TRACK) {
	fprintf(stream, " "); // Print a space instead of a minus sign for better looking alignment
	fprintf(stream, tempo3_parfile_descr.print_format[curparam], parfile->track);
      }else if(curparam >= TEMPO3_PARAM_WAVESIN && curparam < TEMPO3_PARAM_WAVESIN + TEMPO3_MaxNrWaves) {
	fprintf(stream, " "); // Print a space instead of a minus sign for better looking alignment
	fprintf(stream, tempo3_parfile_descr.print_format[curparam], parfile->parameter[curparam]);
	fprintf(stream, " ");
	fprintf(stream, tempo3_parfile_descr.print_format[curparam], parfile->parameter[curparam-TEMPO3_PARAM_WAVESIN+TEMPO3_PARAM_WAVECOS]);
      }else {
	int counter, nrspacesadded;
	long double value;
	if(value_id == 0) {
	  value = parfile->parameter[curparam];
	  if(curparam >= TEMPO3_PARAM_DM && curparam <= TEMPO3_PARAM_DM+TEMPO3_MaxNrDMderivatives) {
	    value += parfile->dm[curparam - TEMPO3_PARAM_DM];
	  }else if(curparam == TEMPO3_PARAM_PMRA) {
	    value += parfile->pmra;
	  }else if(curparam == TEMPO3_PARAM_PMDEC) {
	    value += parfile->pmdec;
	  }
	  nrspacesadded = 0;
	  if(tempo3_parfile_descr.print_format[curparam][strlen(tempo3_parfile_descr.print_format[curparam])-1] == 'f' && hex == 0) {
	    //	    printf("PARAM=%d value=%Lf\n", curparam, value);
	    fflush(stdout);
	    if(fabsl(value) < 10000) { // To align the dots
	      fprintf(stream, " ");
	      nrspacesadded++;
	    }
	    if(fabsl(value) < 1000) {  // To align the dots
	      fprintf(stream, " ");
	      nrspacesadded++;
	    }
	    if(fabsl(value) < 100) {  // To align the dots
	      fprintf(stream, " ");
	      nrspacesadded++;
	    }
	    if(fabsl(value) < 10) { // To align the dots
	      fprintf(stream, " ");
	      nrspacesadded++;
	    }
	  }else {
	    fprintf(stream, "    ");
	  }
	  if(value >= 0)
	    fprintf(stream, " "); // Print a space instead of a minus sign for better looking alignment
	  if(hex == 0) {
	    fprintf(stream, tempo3_parfile_descr.print_format[curparam], value);
	  }else {
	    print_longdouble_as_hexstring(stream, value);
	  }
	  if(value < 0)
	    fprintf(stream, " "); // Print a space instead of a minus sign for better looking alignment
	  if(tempo3_parfile_descr.print_format[curparam][strlen(tempo3_parfile_descr.print_format[curparam])-1] == 'f') {
	    for(counter = 0; counter < 4-nrspacesadded; counter++)
	      fprintf(stream, " ");
	  }
	}else {
	  if(value_id > 0) {
	    fprintf(stream, " VALUE%d", value_id);
	  }else {
	    fprintf(stream, "-VALUE%d", value_id);
	  }
	}
      }
      if(value_id == 0 && curparam != TEMPO3_PARAM_PSRJNAME && curparam != TEMPO3_PARAM_WAVEOM && curparam != TEMPO3_PARAM_PEPOCH && curparam != TEMPO3_PARAM_WAVEEPOCH && curparam != TEMPO3_PARAM_TRES && curparam != TEMPO3_PARAM_BINARY && (curparam < TEMPO3_PARAM_GLEPMIN || curparam >= TEMPO3_PARAM_GLEPMIN+TEMPO3_MaxNrGlitches) && (curparam < TEMPO3_PARAM_GLEPMAX || curparam >= TEMPO3_PARAM_GLEPMAX+TEMPO3_MaxNrGlitches) && curparam != TEMPO3_PARAM_MODE && curparam != TEMPO3_PARAM_TRACK && curparam != TEMPO3_PARAM_TZRMJD && curparam != TEMPO3_PARAM_TZRFRQ && curparam != TEMPO3_PARAM_TZRSITE && curparam != TEMPO3_PARAM_POSEPOCH && curparam != TEMPO3_PARAM_DMEPOCH) {
	if(parfile2 != NULL) {
	  fprintf(stream, "%-10s ", tempo3_parfile_descr.units[curparam]);	  
	}
	if(ignorefitted == 0) {
	  if(eph->fixed[curparam])
	    fprintf(stream, "fixed  ");
	  else
	    fprintf(stream, "FIT    ");
	}
	if(parfile2 != NULL) {
	  if(curparam >= TEMPO3_PARAM_DM && curparam <= TEMPO3_PARAM_DM + TEMPO3_MaxNrDMderivatives) {
	    fprintf(stream, "%11.1Le", parfile->dm[curparam - TEMPO3_PARAM_DM] + parfile->parameter[curparam] - parfile2->dm[curparam - TEMPO3_PARAM_DM] - parfile2->parameter[curparam]);
	  }else if(curparam == TEMPO3_PARAM_PMRA) {
	    fprintf(stream, "%11.1Le", parfile->pmra + parfile->parameter[curparam] - parfile2->pmra - parfile2->parameter[curparam]);
	  }else if(curparam == TEMPO3_PARAM_PMDEC) {
	    fprintf(stream, "%11.1Le", parfile->pmdec + parfile->parameter[curparam] - parfile2->pmdec - parfile2->parameter[curparam]);
	  }else if(curparam >= TEMPO3_PARAM_WAVESIN && curparam < TEMPO3_PARAM_WAVESIN + TEMPO3_MaxNrWaves) {
	    fprintf(stream, "%11.1Le %11.1Le", parfile->parameter[curparam]-parfile2->parameter[curparam], parfile->parameter[curparam-TEMPO3_PARAM_WAVESIN+TEMPO3_PARAM_WAVECOS]-parfile2->parameter[curparam-TEMPO3_PARAM_WAVESIN+TEMPO3_PARAM_WAVECOS]);   
	  }else {
	    fprintf(stream, "%11.1Le", parfile->parameter[curparam]-parfile2->parameter[curparam]);
	  }
	}
	if(ignoreerrors == 0) {
	  if(curparam >= TEMPO3_PARAM_DM && curparam <= TEMPO3_PARAM_DM + TEMPO3_MaxNrDMderivatives) {
	    fprintf(stream, "%s%12.2Le (%.2Le)", txt, parfile->dplus[curparam], fabsl((parfile->dm[curparam-TEMPO3_PARAM_DM] + parfile->parameter[curparam])/parfile->dplus[curparam]));
	  }else if(curparam == TEMPO3_PARAM_PMRA) {
	    fprintf(stream, "%s%12.2Le (%.2Le)", txt, parfile->dplus[curparam], fabsl((parfile->pmra + parfile->parameter[curparam])/parfile->dplus[curparam]));
	  }else if(curparam == TEMPO3_PARAM_PMDEC) {
	    fprintf(stream, "%s%12.2Le (%.2Le)", txt, parfile->dplus[curparam], fabsl((parfile->pmdec + parfile->parameter[curparam])/parfile->dplus[curparam]));
	  }else if(curparam == TEMPO3_PARAM_RAJ) {
	    fprintf(stream, "%s%12.2Le (sec)", txt, parfile->dplus[curparam]*12.0*60.0*60.0/M_PI);
	  }else if(curparam == TEMPO3_PARAM_DECJ) {
	    fprintf(stream, "%s%12.2Le (arcsec)", txt, parfile->dplus[curparam]*180.0*60.0*60.0/M_PI);
	  }else if(curparam >= TEMPO3_PARAM_WAVESIN && curparam < TEMPO3_PARAM_WAVESIN + TEMPO3_MaxNrWaves) {
	    fprintf(stream, "%s%12.2Le  %12.2Le", txt, parfile->dplus[curparam], parfile->dplus[curparam-TEMPO3_PARAM_WAVESIN+TEMPO3_PARAM_WAVECOS]);
	  }else {
	    fprintf(stream, "%s%12.2Le (%.2Le)", txt, parfile->dplus[curparam], fabsl(parfile->parameter[curparam]/parfile->dplus[curparam]));
	  }
	}
      }
      fprintf(stream, "\n");
    }
  }

  if(toas != NULL) {
    fprintf(stream, "%sSTART       %Lf\n", txt, toas->data_mjdmin);
    fprintf(stream, "%sFINISH      %Lf\n", txt, toas->data_mjdmax);
    /* HaveALookHereIfFixingChi2BugZoomed */
    fprintf(stream, "%sNTOA        %d\n", txt, toas->nrtoas);
  }

  if(ignorewraps == 0) {
    for(i = 0; i < eph->wraps->nrwraps; i++) {
      fprintf(stream, "%sPHASE      %.0Lf %Lf\n", txt, eph->wraps->phase[i], eph->wraps->mjd[i]);
    }
  }

  if(ignoreunrecognizedlines == 0) {
    for(i = 0; i < eph->nrunused_parfile_lines; i++) {
      fprintf(stream, "%s", eph->unused_parfile_lines[i]);
    }
  }

#ifdef TEMPO3_EnableGTK
  if(GTKinitialised) {
    void print_ephemeris_gtk(FILE *stream, tempo3_parameters_def parfile, tempo3_parameters_def *parfile2, int *fixed, int ident, int otherparams, toa_data_def *toas, tempo3_wraps_def *wraps, char **unused_parfile_lines, int nrunused_parfile_lines, int *paramset, int finderrors, long double *dplus, long double *dmin, int *paramlinked);
    int otherparams, finderrors;
    if(ignoreunrecognizedlines) {
      otherparams = 0;
    }else {
      otherparams = 1;
    }
    if(ignoreerrors) {
      finderrors = 0;
    }else {
      finderrors = 1;
    }
    int *gtkfixed;
    if(ignorefitted) {
      gtkfixed = NULL;
    }else {
      gtkfixed = eph->fixed;
    }
    print_ephemeris_gtk(stream, *parfile, parfile2, gtkfixed, verbose.indent, otherparams, toas, eph->wraps, eph->unused_parfile_lines, eph->nrunused_parfile_lines, eph->paramset,  finderrors, parfile->dplus, parfile->dmin, eph->paramlinked);
  }
#endif
  return 1;
}

// Allocates memory in the ephemeris to store information, and
// initialise their values to be unset. If allowerrors is nonzero,
// space is allocated to store errors on parameters.
//
// Return 0: error
// Return 1: ok
int initialise_ephemeris(ephemeris_def *eph, int allowerrors, verbose_definition verbose)
{
  long i;
  //    char unused_parfile_lines[TEMPO3_MaxNrParLines][TEMPO3_MaxNrParLineLength];
  //    int nrunused_parfile_lines;
  //    int *fixed, *paramset, *paramlinked;
  eph->fixed = malloc(TEMPO3_PARAM_TOTAL*sizeof(int));
  eph->paramset = malloc(TEMPO3_PARAM_TOTAL*sizeof(int));
  eph->paramlinked = malloc((1+2*TEMPO3_PARAM_TOTAL)*sizeof(int));
  eph->wraps = malloc(sizeof(tempo3_wraps_def));
  eph->unused_parfile_lines = malloc(TEMPO3_MaxNrParLines*sizeof(char *));
  if(eph->paramset == NULL || eph->paramlinked == NULL || eph->fixed == NULL || eph->wraps == NULL || eph->unused_parfile_lines == NULL) {
    printerror(verbose.debug, "ERROR initialise_ephemeris: Cannot allocate memory\n");
    return 0;
  }
  /*
  for(i = 0; i < TEMPO3_MaxNrParLines; i++) {
    eph->unused_parfile_lines[i] = malloc(TEMPO3_MaxNrParLineLength+1);
    if(eph->unused_parfile_lines[i] == NULL) {
      printerror(verbose.debug, "ERROR initialise_ephemeris: Cannot allocate memory\n");
      return 0;
    }
  }
  */
  eph->nrunused_parfile_lines = 0;

  eph->wraps->nrwraps = 0;
  eph->paramlinked[0] = 0;  // Means: no parfile parameters are linked to VALUE parameters
  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
    eph->fixed[i] = 1;
    eph->paramset[i] = 0;
  }
  eph->fixed[TEMPO3_PARAM_PHASE0] = 0;
  eph->paramset[TEMPO3_PARAM_PHASE0] = 1;         // Could leave this out, and it will not show up in the output
  if(allocate_ephemeris_par_only(&(eph->par), allowerrors, verbose) == 0) {
    printerror(verbose.debug, "ERROR initialise_ephemeris: Cannot allocate memory\n");
    return 0;
  }
  initialise_ephemeris_par_only(&(eph->par));
  return 1;
}

// Releases the memory allocated by initialise_ephemeris().
//
// Return 0: error
// Return 1: ok
int free_ephemeris(ephemeris_def *eph, verbose_definition verbose)
{
  long i;
  free(eph->fixed);
  eph->fixed = NULL;
  free(eph->paramset);
  eph->paramset = NULL;
  free(eph->paramlinked);
  eph->paramlinked = NULL;
  free(eph->wraps);
  eph->wraps = NULL;

  for(i = 0; i < eph->nrunused_parfile_lines; i++) {
    free(eph->unused_parfile_lines[i]);
  }
  free(eph->unused_parfile_lines);
  eph->nrunused_parfile_lines = 0;

  free_ephemeris_par_only(&(eph->par));
  return 1;
}

// Reads in an ephemeris from a file. The ephemeris needs to be
// initialised first with the function initialise_ephemeris(). If hex
// is nonzero the values in the parameter file are assumed to be a hex
// representation of the bytes of the floating point value, rather
// than decimal numbers.
//
// Return 0 on error, 1 on success
int read_ephemeris(char *filename, ephemeris_def *eph, int hex, int nowarning, verbose_definition verbose)
{
  int i, ret, found, curparam, altid, jumpmeaning;
  char txt[1000], param[1000], value[1000], value2[1000], value3[1000], value4[1000], value5[1000], jumpstr[1000];
  char valueD[1000];
  long double jumpvalue;
  FILE *fin;
  char *ignorewarning;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Reading %s\n", filename);
  }

  if(eph->nrunused_parfile_lines != 0) {
    printwarning(verbose.debug, "WARNING read_ephemeris: memory slot where information will be stored is either uninitialised, or already in use.");
  }
  eph->nrunused_parfile_lines = 0;

  ignorewarning = malloc(10000);
  if(ignorewarning == NULL) {
    printerror(verbose.debug, "ERROR read_ephemeris: Memory allocation error");
    return 0;
  }
  ignorewarning[0] = 0;

  fin = fopen(filename, "r");
  if(fin == NULL) {
    printerror(verbose.debug, "ERROR read_ephemeris: Cannot open %s\n", filename);
    return 0;
  }
  eph->par.nrglitches = 0;
  while(fgets(txt, 1000, fin) != NULL) {
    /* First break up string in PARAM STRINGS VALUE for jump parameter */
    jumpstr[0] = 0;
    param[0] = 0;
    ret = sscanf(txt, "%s %s %s %s %s %s", param, value, value2, value3, value4, value5);
    if(ret >= 1) {   // Skip empty lines
    /*
    if(ret >= 3) {
      sscanf(value2, "%Lf", &jumpvalue);
      strcat(jumpstr, value);
    }
    if(ret >= 4) {
      sscanf(value3, "%Lf", &jumpvalue);
      strcat(jumpstr, " ");
      strcat(jumpstr, value2);
    }
    if(ret >= 5) {
      sscanf(value4, "%Lf", &jumpvalue);
      strcat(jumpstr, " ");
      strcat(jumpstr, value3);
    }
    if(ret >= 6) {
      sscanf(value5, "%Lf", &jumpvalue);
      strcat(jumpstr, " ");
      strcat(jumpstr, value4);
    }
    */
    // JUMPS are assumed to be in the format JUMP ID1 ID2 VALUE OPTIONAL_ERROR
    // NOT ALWAYS TRUE, AS ACCORDING TO TEMPO2 DOCUMENTATION:
    // JUMP FLAGID FLAGSTRING   // default (i.e JUMP -t lovell)
    // JUMP NAME str   // str part of filename
    // JUMP TEL ID
    // JUMP MJD v1 v2   // In MJD
    // JUMP FREQ v1 v2   // In MHz
    if(strcasecmp(value, "TEL") == 0) {
      strcpy(jumpstr, value2);
      if(hex == 0) {
	sscanf(value3, "%Lf", &jumpvalue);
      }else {
	if(convert_hexstring_to_longdouble(value3, &jumpvalue, verbose) == 0) {
	  printerror(verbose.debug, "ERROR read_ephemeris: Hex conversion failed.", value);
	  return 0;
	}
      }
      jumpmeaning = 1;
    }else if(strcasecmp(value, "NAME") == 0 || strcasecmp(value, "MJD") == 0 || strcasecmp(value, "FREQ") == 0) {
      printerror(verbose.debug, "ERROR read_ephemeris: The JUMP %s parameters are currently not supported in tempo3.", value);
      return 0;
    }else {
      if(ret >= 4) {
	strcpy(jumpstr, value);
	strcat(jumpstr, " ");
	strcat(jumpstr, value2);
	if(hex == 0) {
	  sscanf(value3, "%Lf", &jumpvalue);
	}else {
	  if(convert_hexstring_to_longdouble(value3, &jumpvalue, verbose) == 0) {
	    printerror(verbose.debug, "ERROR read_ephemeris: Hex conversion failed.", value);
	    return 0;
	  }
	}
	jumpmeaning = 0;
      }
    }

    /* Now make the normal PARAM VALUE parts. Value can actually be two values (for fitwaves)*/
    ret = sscanf(txt, "%s %s %s", param, value, value2);
    if(ret == 0 || ret == EOF) {
      param[0] = 0;
    }
    if(ret == 3) {
      strcat(value, " ");
      strcat(value, value2);
    }
    strcpy(valueD, value); 
    /* Tempo1 uses D-10 instead of e-10 */
    if(hex == 0) {
      int i;
      for(i = 0; i < strlen(valueD); i++) {
	if(valueD[i] == 'D') {
	  valueD[i] = 'e';
	}
      }
    }
    
    found = 0;

    if(strcmp(param, "JUMP") == 0) {   // Dealt with separately, because all JUMP parameters have the same identifier
      eph->par.parameter[TEMPO3_PARAM_JUMPS+eph->par.nrjumps] = jumpvalue;
      strcpy(eph->par.jumpflag[eph->par.nrjumps], jumpstr);
      eph->paramset[TEMPO3_PARAM_JUMPS + eph->par.nrjumps] = 1;
      eph->par.jumpmeaning[eph->par.nrjumps] = jumpmeaning;
      eph->par.nrjumps += 1;
      found = 1;
    }else if(strcmp(param, "BINARY") == 0) { // Dealt with separately, since code needs to be converted in an integer id, and it is stored in a different place
      if(strcmp(value, "BT") == 0) {
	eph->par.binarymodel = 1;
	eph->paramset[TEMPO3_PARAM_BINARY] = 1;
      }else if(strcmp(value, "T2") == 0) {
	eph->par.binarymodel = 2;
	eph->paramset[TEMPO3_PARAM_BINARY] = 1;
      }else if(strcmp(value, "HYPER") == 0) {
	eph->par.binarymodel = 3;
	eph->paramset[TEMPO3_PARAM_BINARY] = 1;
      }else {
	if(nowarning == 0) {
	  printwarning(verbose.debug, "WARNING read_ephemeris: BINARY model %s is not supported and the binary parameters are ignored.", value);
	}
	eph->par.binarymodel = 0;
      }
      found = 1;
    }else {
      for(curparam = 0; curparam < TEMPO3_PARAM_TOTAL; curparam++) {
	for(altid = 0; altid < 2; altid++) {
	  if(tempo3_parfile_descr.identifiers[curparam][altid] != NULL) {
	    if(strcasecmp(param, tempo3_parfile_descr.identifiers[curparam][altid]) == 0) {
	      eph->paramset[curparam] = 1;
	      
	      // Check first if value of parameter is set to the variable VALUE
	      int id;
	      if(strncmp(value, "VALUE", 5) == 0) {
		sscanf(value, "VALUE%d", &id);
		eph->paramlinked[2*eph->paramlinked[0] + 1] = curparam;
		eph->paramlinked[2*eph->paramlinked[0] + 2] = id;
		eph->paramlinked[0]++;
	      }else if(strncmp(value, "-VALUE", 6) == 0) {
		sscanf(value, "-VALUE%d", &id);
		eph->paramlinked[2*eph->paramlinked[0] + 1] = curparam;
		eph->paramlinked[2*eph->paramlinked[0] + 2] = -id;
		eph->paramlinked[0]++;
	      }else {                               // Read in value itself rather than it being VALUE
		if(curparam == TEMPO3_PARAM_PSRJNAME) {
		  sscanf(value, "%s", eph->par.psrjname);
		}else if(curparam == TEMPO3_PARAM_RAJ) {
		  converthms_ld(value, &(eph->par.raj));
		  eph->par.raj *= M_PI/12.0;          // RAJ stored here, the difference is fitted for
		}else if(curparam == TEMPO3_PARAM_DECJ) {
		  converthms_ld(value, &(eph->par.decj));
		  eph->par.decj *= M_PI/180.0;        // DECJ stored here, the difference is fitted for
		}else if(curparam == TEMPO3_PARAM_PMRA) {
		  if(hex == 0) {
		    sscanf(valueD, "%Lf", &(eph->par.pmra));                 // PMRA stored in .pmra, but difference is fitted
		  }else {
		    if(convert_hexstring_to_longdouble(valueD, &(eph->par.pmra), verbose) == 0) {
		      printerror(verbose.debug, "ERROR read_ephemeris: Hex conversion failed.", value);
		      return 0;
		    }
		  }
		}else if(curparam == TEMPO3_PARAM_PMDEC) {
		  if(hex == 0) {
		    sscanf(valueD, "%Lf", &(eph->par.pmdec));                 // PMRA stored in .pmra, but difference is fitted
		  }else {
		    if(convert_hexstring_to_longdouble(valueD, &(eph->par.pmdec), verbose) == 0) {
		      printerror(verbose.debug, "ERROR read_ephemeris: Hex conversion failed.", value);
		      return 0;
		    }
		  }
		}else if(curparam == TEMPO3_PARAM_TZRSITE) {
		  sscanf(value, "%s", eph->par.tzrsite);
		}else if(curparam >= TEMPO3_PARAM_DM && curparam <= TEMPO3_PARAM_DM+TEMPO3_MaxNrDMderivatives) {
		  if(hex == 0) {
		    sscanf(valueD, "%Lf", &(eph->par.dm[curparam-TEMPO3_PARAM_DM]));                 // DM stored in .dm, but .dmdiff is fitted
		  }else {
		    if(convert_hexstring_to_longdouble(valueD, &(eph->par.dm[curparam-TEMPO3_PARAM_DM]), verbose) == 0) {
		      printerror(verbose.debug, "ERROR read_ephemeris: Hex conversion failed.", value);
		      return 0;
		    }
		  }
		}else if(curparam == TEMPO3_PARAM_MODE) {
		  sscanf(value, "%d", &(eph->par.mode));                 // The value of mode is stored in a non-standard location
		}else if(curparam == TEMPO3_PARAM_TRACK) {
		  sscanf(value, "%d", &(eph->par.track));                // The value of track is stored in a non-standard location
		  if(eph->par.track != 0 && eph->par.track != -2) {
		    if(nowarning == 0) {
		      printwarning(verbose.debug, "WARNING read_ephemeris: TRACK %d is not recognized as a valid parameter and is set to zero.", eph->par.track);
		    }
		    eph->par.track = 0;
		  }
		}else {
		  if(hex == 0) {
		    sscanf(valueD, "%Le", &(eph->par.parameter[curparam]));
		  }else {
		    if(convert_hexstring_to_longdouble(valueD, &(eph->par.parameter[curparam]), verbose) == 0) {
		      printerror(verbose.debug, "ERROR read_ephemeris: Hex conversion failed.", value);
		      return 0;
		    }
		  }
		}
	      }
	      if(curparam >= TEMPO3_PARAM_GLEP && curparam < TEMPO3_PARAM_GLEP + TEMPO3_MaxNrGlitches) {
		if(curparam - TEMPO3_PARAM_GLEP + 1 > eph->par.nrglitches)
		  eph->par.nrglitches = curparam - TEMPO3_PARAM_GLEP + 1;
	      }else if(curparam >= TEMPO3_PARAM_T0 && curparam < TEMPO3_PARAM_T0 + TEMPO3_MaxNrCompanions) {
		eph->par.companion_defined[curparam-TEMPO3_PARAM_T0] = 1;
	      }else if(curparam >= TEMPO3_PARAM_TASC && curparam < TEMPO3_PARAM_TASC + TEMPO3_MaxNrCompanions) {
		eph->par.companion_defined[curparam-TEMPO3_PARAM_TASC] = 1;
	      }
	      //	      printf("XXXX read in parameter %s\n", parfile_identifiers[curparam][altid]);
	      found = 1;
	      break;
	    }
	  }else {
	    //	    printf("XXXX didn't match parameter %s\n", parfile_identifiers[curparam][altid]);
	  }
	}
	if(found)
	  break;
      }
    }

    // TRES parameter is always calculated, so it should always be active
    eph->paramset[TEMPO3_PARAM_TRES] = 1;
    // MODE parameter is always set
    eph->paramset[TEMPO3_PARAM_MODE] = 1;
    // TRACK parameter is always set
    eph->paramset[TEMPO3_PARAM_TRACK] = 1;

    //    if(found == 0) {
    //      printf("XXXX didn't recognized parameter '%s'\n", param);
    //    }
    if(strcmp(param, "PHASE") == 0) {
      sscanf(valueD, "%Lf %Lf", &(eph->wraps->phase[eph->wraps->nrwraps]), &(eph->wraps->mjd[eph->wraps->nrwraps]));
      eph->wraps->nrwraps = eph->wraps->nrwraps + 1;
      /* Just ignore start and finish entries and NTOA, regenerage them anyway */
    }else if(strcmp(param, "START") == 0) {
    }else if(strcmp(param, "FINISH") == 0) {
    }else if(strcmp(param, "NTOA") == 0) {
    }else {
      int i;
      for(i = 0; i < TEMPO3_MaxNrWaves; i++) {
	sprintf(value2, "WAVE%d", i+1);
	if(strcmp(param, value2) == 0) {
	  eph->paramset[TEMPO3_PARAM_WAVESIN+i] = 1;
	  eph->paramset[TEMPO3_PARAM_WAVECOS+i] = 1;
	  sscanf(valueD, "%Lf %Lf", &(eph->par.parameter[TEMPO3_PARAM_WAVESIN+i]), &(eph->par.parameter[TEMPO3_PARAM_WAVECOS+i]));
	  if(eph->par.nrwaves < i+1)
	    eph->par.nrwaves = i+1;
	  found = 1;
	}	
      }
      
      if(found == 0) {
	if(ignorewarning[0] == 0) {
	  strcpy(ignorewarning, "WARNING: Following parameters are ignored from the input par file: ");
	}
	strcat(ignorewarning, param);
	strcat(ignorewarning, " ");
	if(eph->nrunused_parfile_lines < TEMPO3_MaxNrParLines) {
	  eph->unused_parfile_lines[eph->nrunused_parfile_lines] = malloc(strlen(txt)+1);
	  if(eph->unused_parfile_lines[eph->nrunused_parfile_lines] == NULL) {
	    printerror(verbose.debug, "ERROR read_ephemeris: Cannot allocate memory\n");
	    return 0;
	  }
	  strcpy(eph->unused_parfile_lines[eph->nrunused_parfile_lines], txt);
	  (eph->nrunused_parfile_lines) += 1;
	}else {
	  if(nowarning == 0) {
	    printwarning(verbose.debug, "WARNING read_ephemeris: Too many lines in par file which could not be interpreted");
	  }
	}
      }
    }
    }
  }   /* End of while loop reading in file */

  for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
    if((eph->paramset[TEMPO3_PARAM_EPS1+i] || eph->paramset[TEMPO3_PARAM_EPS2+i]) && (eph->paramset[TEMPO3_PARAM_ECC+i] || eph->paramset[TEMPO3_PARAM_OM+i])) {
      printerror(verbose.debug, "ERROR read_ephemeris: Cannot define EPS1/EPS2 and ECC/OM parameters simultaneously for companion %d.", i+1);
      return 0;
    }
    if(eph->par.parameter[TEMPO3_PARAM_T0] != 0 && eph->par.parameter[TEMPO3_PARAM_TASC] != 0) {
      printerror(verbose.debug, "ERROR read_ephemeris: Cannot define T0 and TASC parameters simultaneously for companion %d.", i+1);
      return 0;
    }
    if(eph->paramset[TEMPO3_PARAM_PBDOT+i]) {
      if(fabsl(eph->par.parameter[TEMPO3_PARAM_PBDOT+i]) > 1e-7) {
	if(nowarning == 0) {
	  printwarning(verbose.debug, "WARNING read_ephemeris: PBDOT_%d appears to be abnormally large. Assume the \"tempo\" factor 1e-12 is missing.", i+1);
	}
	eph->par.parameter[TEMPO3_PARAM_PBDOT+i] *= 1e-12;
      }
    }
  }


  if(ignorewarning[0] != 0) {
    if(nowarning == 0) {
      printwarning(verbose.debug, ignorewarning);
    }
  }

  /* Assume wave epoch is period epoch if not set */
  if(eph->paramset[TEMPO3_PARAM_WAVEEPOCH] == 0) {
    eph->par.parameter[TEMPO3_PARAM_WAVEEPOCH] = eph->par.parameter[TEMPO3_PARAM_PEPOCH]; 
  }
  // WAVEEPOCH should always be enabled if waves are defined
  if(eph->paramset[TEMPO3_PARAM_WAVESIN]) {
    eph->paramset[TEMPO3_PARAM_WAVEEPOCH] = 1;
  }
  
  // Assume pos epoch is period epoch if not set if position defined
  if(eph->paramset[TEMPO3_PARAM_POSEPOCH] == 0) {
    eph->par.parameter[TEMPO3_PARAM_POSEPOCH] = eph->par.parameter[TEMPO3_PARAM_PEPOCH]; 
    //    printf("XXXX setting POSEPOCH = %Lf\n", eph->par.parameter[TEMPO3_PARAM_POSEPOCH]);
  }
  //  printf("XXXX PEPOCH = %Lf (param %d)\n", eph->par.parameter[TEMPO3_PARAM_PEPOCH], TEMPO3_PARAM_PEPOCH);
  //  printf("XXXX POSEPOCH = %Lf (param %d)\n", eph->par.parameter[TEMPO3_PARAM_POSEPOCH], TEMPO3_PARAM_POSEPOCH);
  if(eph->paramset[TEMPO3_PARAM_RAJ] != 0 || eph->paramset[TEMPO3_PARAM_DECJ] != 0) {
    eph->paramset[TEMPO3_PARAM_POSEPOCH] = 1;
  }

  // Assume DM epoch is period epoch if not set if DM1 etc defined
  if(eph->paramset[TEMPO3_PARAM_DMEPOCH] == 0) {
    eph->par.parameter[TEMPO3_PARAM_DMEPOCH] = eph->par.parameter[TEMPO3_PARAM_PEPOCH]; 
    //    printf("XXXX setting DMEPOCH = %Lf\n", eph->par.parameter[TEMPO3_PARAM_DMEPOCH]);
  }
  //  printf("XXXX DMEPOCH = %Lf (param %d)\n", eph->par.parameter[TEMPO3_PARAM_DMEPOCH], TEMPO3_PARAM_DMEPOCH);
  for(i = 0; i < TEMPO3_MaxNrDMderivatives; i++) {
    if(eph->paramset[TEMPO3_PARAM_DM+i+1] != 0) {
      eph->paramset[TEMPO3_PARAM_DMEPOCH] = 1;
    }
  }

  fclose(fin);

  /* Check if glitch epochs are not outside allowed range */

  for(i = 0; i < TEMPO3_MaxNrGlitches; i++) {
    if(eph->par.parameter[TEMPO3_PARAM_GLEPMIN+i] > 0 && eph->par.parameter[TEMPO3_PARAM_GLEP+i] > 0 && (eph->par.parameter[TEMPO3_PARAM_GLEP+i] < eph->par.parameter[TEMPO3_PARAM_GLEPMIN+i])) {
      printerror(verbose.debug, "ERROR read_ephemeris: Glitch epoch out of range: %Lf < %Lf\n", eph->par.parameter[TEMPO3_PARAM_GLEP+i], eph->par.parameter[TEMPO3_PARAM_GLEPMIN+i]);
      return 0;
    }
    if(eph->par.parameter[TEMPO3_PARAM_GLEPMAX+i] > 0 && eph->par.parameter[TEMPO3_PARAM_GLEP+i] > 0 && (eph->par.parameter[TEMPO3_PARAM_GLEP+i] > eph->par.parameter[TEMPO3_PARAM_GLEPMAX+i])) {
      printerror(verbose.debug, "ERROR read_ephemeris: Glitch epoch out of range: %Lf > %Lf\n", eph->par.parameter[TEMPO3_PARAM_GLEP+i], eph->par.parameter[TEMPO3_PARAM_GLEPMAX+i]);
      return 0;
    }
  }

  free(ignorewarning);

  return 1;
}

// Should really be abandoned. It is also declared in tempo3 as well.
//
// It calculates the predicted phase of the pulsar at mjd and 0.5/F0
// appart to derive the actual perceived spin frequency.
long double ephemeris_estimateNu_numerically_avoid_using(long double mjd, long double freq, long double ssb1, long double ssb2, long double ssb3, char *site, char *flags, tempo3_parameters_def *parfile, int *paramset)
{
  long double phase, dt, phase2;
  phase = evaluate_ephemeris(mjd, freq, ssb1, ssb2, ssb3, 0, 0, site, flags, parfile, 0, paramset, 0, 0);
  dt = 0.5/parfile->parameter[TEMPO3_PARAM_F0];
  /* Calculate phases ~0.5 periods appart */
  phase2 = evaluate_ephemeris(mjd+dt/TEMPO3_SecPerDay, freq, ssb1, ssb2, ssb3, 0, 0, site, flags, parfile, 0, paramset, 0, 0);
  /* Calculate actual period */
  phase2 -= phase;
  phase2 = phase2/dt;
  return phase2;
}


// Should really be abandoned. It is also declared in tempo3 as well.
//
// This function could be made more precise by calculating the
// difference in frequency, rather than phase, which is possible with
// evaluate_ephemeris(). In any case, it needs to be avoided.
long double ephemeris_estimateNudot_numerically_avoid_using(long double mjd, long double freq, long double ssb1, long double ssb2, long double ssb3, char *site, char *flags, tempo3_parameters_def parfile, int *paramset)
{
  long double daysUsedForNdotEstimate;
  long double phase, dt, phase2, nudot;
  /* Used for 1119 paper */
  /*  daysUsedForNdotEstimate = 19.0*1; */  
  /* After subtracting nu part of phase, much shorter time intervals is required to get significant effects. By including the nu part, the phases are much larger number, and therefore you run into numerical problems when deriving such small numbers. */
  daysUsedForNdotEstimate = 0.01;
  //  fprintf(stderr, "XXXX pepoch in parfile=%Lf phase0=%Lf\n", parfile.parameter[TEMPO3_PARAM_PEPOCH], parfile.parameter[TEMPO3_PARAM_PHASE0]);
  phase = evaluate_ephemeris(mjd+daysUsedForNdotEstimate*0.5, freq, ssb1, ssb2, ssb3, 0, 1, site, flags, &parfile, 0, paramset, 0, 0);
  //  fprintf(stderr, "XXXX phase=%Lf\n", phase);
  dt = 0.5/parfile.parameter[TEMPO3_PARAM_F0];
  /* Calculate phases ~0.5 periods appart */
  phase2 = evaluate_ephemeris(mjd+daysUsedForNdotEstimate*0.5+dt/TEMPO3_SecPerDay, freq, ssb1, ssb2, ssb3, 0, 1, site, flags, &parfile, 0, paramset, 0, 0);
  /* Calculate actual period */
  phase2 -= phase;
  phase2 = phase2/dt;
  nudot = phase2;
  //  fprintf(stderr, "XXXX nu(t1)=%Lf\n", nudot);
      
  /* Now do the same one day later to determine F1 */
  phase = evaluate_ephemeris(mjd-daysUsedForNdotEstimate*0.5, freq, ssb1, ssb2, ssb3, 0, 1, site, flags, &parfile, 0, paramset, 0, 0);
  dt = 0.5/parfile.parameter[TEMPO3_PARAM_F0];
  /* Calculate phases ~0.5 periods appart */
  phase2 = evaluate_ephemeris(mjd-daysUsedForNdotEstimate*0.5+dt/TEMPO3_SecPerDay, freq, ssb1, ssb2, ssb3, 0, 1, site, flags, &parfile, 0, paramset, 0, 0);
  /* Calculate actual period */
  phase2 -= phase;
  phase2 = phase2/dt;
  //  fprintf(stderr, "XXXX nu(t2)=%Lf\n", phase2);
  /* This is delta F per day */
  nudot -= phase2;
  nudot /= daysUsedForNdotEstimate*TEMPO3_SecPerDay;
  //  fprintf(stderr, "XXXX nudot=%Le\n", nudot);
  return nudot;
}

// Should really only be used internally in library, but is declared in tempo3 as well.
// Could easily make a ephemeris version with as a wrapper function.
/* Returns fitwave function in seconds. 

If mode    = 0, it calculates the phase (in seconds)
           = 1, it calculates the frequency  (in seconds per mjd)
           = 2, it calculates the frequency-devivative (in seconds per mjd**2)
           = 3, it calculates the second frequency-derivative (in seconds per mjd**3)

If positive, showWaves is the first fitwave to be excluded from the sum
If zero, all terms are included
If negative, -showWaves it is the first fitwave to be included
*/
long double fitwave_function(tempo3_parameters_def *parfile, long double mjd, int mode, int showWaves)
{
  int j, ok;
  long double y;
  y = 0;
  for(j = 0; j < parfile->nrwaves; j++) {
    ok = 0;
    if(showWaves == 0) {
      ok = 1;
    }else if(showWaves > 0) {
      if((j+1) < showWaves)
	ok = 1;
    }else if(showWaves < 0) {
      if((j+1) >= -showWaves)
	ok = 1;
    }
    if(ok) {
      if(mode == 0) {
	y += parfile->parameter[TEMPO3_PARAM_WAVESIN+j]*sinl(parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1)*(mjd-parfile->parameter[TEMPO3_PARAM_WAVEEPOCH]));
	y += parfile->parameter[TEMPO3_PARAM_WAVECOS+j]*cosl(parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1)*(mjd-parfile->parameter[TEMPO3_PARAM_WAVEEPOCH]));
      }else if(mode == 1) {
	y += parfile->parameter[TEMPO3_PARAM_WAVESIN+j]*cosl(parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1)*(mjd-parfile->parameter[TEMPO3_PARAM_WAVEEPOCH]))*parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1);
	y += -parfile->parameter[TEMPO3_PARAM_WAVECOS+j]*sinl(parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1)*(mjd-parfile->parameter[TEMPO3_PARAM_WAVEEPOCH]))*parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1);
      }else if(mode == 2) {
	y += -parfile->parameter[TEMPO3_PARAM_WAVESIN+j]*sinl(parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1)*(mjd-parfile->parameter[TEMPO3_PARAM_WAVEEPOCH]))*parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1)*parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1);
	y += -parfile->parameter[TEMPO3_PARAM_WAVECOS+j]*cosl(parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1)*(mjd-parfile->parameter[TEMPO3_PARAM_WAVEEPOCH]))*parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1)*parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1);
      }else if(mode == 3) {
	y += -parfile->parameter[TEMPO3_PARAM_WAVESIN+j]*cosl(parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1)*(mjd-parfile->parameter[TEMPO3_PARAM_WAVEEPOCH]))*parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1)*parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1)*parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1);
	y += parfile->parameter[TEMPO3_PARAM_WAVECOS+j]*sinl(parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1)*(mjd-parfile->parameter[TEMPO3_PARAM_WAVEEPOCH]))*parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1)*parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1)*parfile->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(j+1);
      }
    /*    printf("WAVE%d: sin=%Lf cos=%Lf (ep=%Lf om=%Lf)\n",j+1,  parfile->parameter[TEMPO3_PARAM_WAVESIN+j], parfile->parameter[TEMPO3_PARAM_WAVECOS+j], parfile->parameter[TEMPO3_PARAM_WAVEEPOCH], parfile->parameter[TEMPO3_PARAM_WAVEOM]);  */
    }
  }

  /*  printf("XXXXXXXXXXXX %d\n", parfile->nrwaves); */
  /*  y *= -1; */
  return y;
}

// MJD in seconds
// whichcompanion = -1 = all effects combined, otherwise only that of one companion (counting from 0)
// calc_deriv = 0: return time correction
// calc_deriv = 1: return derivative of time correction to a
// calc_deriv = 2: return derivative of time correction to Pb
// calc_deriv = 3: return derivative of time correction to T0
// calc_deriv = 4: return derivative of time correction to omega
long double ephemeris_get_torb(long double mjd, tempo3_parameters_def *parfile, int *paramset, int calc_deriv, int whichcompanion)
{
  int i, idx;
  long double torb;
  binary_def companion[TEMPO3_MaxNrCompanions];
  int nrCompanions;

  nrCompanions = 0;
  if(parfile->binarymodel != 0) {
    if(parfile->companion_defined[0] == 0) {
      printwarning(0, "WARNING ephemeris_get_torb: First companion is undefined, which is required for the any binary model. Ignoring any further binary calculations");
      parfile->binarymodel = 0;
    }

    for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
      if(whichcompanion < 0) {
	idx = i;
      }else {
	idx = whichcompanion;
      }
      companion[i].defined       = parfile->companion_defined[idx];
      companion[i].ecc_defined   = paramset[TEMPO3_PARAM_ECC+idx];
      companion[i].eps1_defined  = paramset[TEMPO3_PARAM_EPS1+idx];
      companion[i].t0            = parfile->parameter[TEMPO3_PARAM_T0+idx];
      companion[i].tasc          = parfile->parameter[TEMPO3_PARAM_TASC+idx];
      companion[i].pb            = parfile->parameter[TEMPO3_PARAM_PB+idx];
      companion[i].a1            = parfile->parameter[TEMPO3_PARAM_A1+idx];
      companion[i].om            = parfile->parameter[TEMPO3_PARAM_OM+idx];
      companion[i].ecc           = parfile->parameter[TEMPO3_PARAM_ECC+idx];
      companion[i].eps1          = parfile->parameter[TEMPO3_PARAM_EPS1+idx];
      companion[i].eps2          = parfile->parameter[TEMPO3_PARAM_EPS2+idx];
      companion[i].omdot         = parfile->parameter[TEMPO3_PARAM_OMDOT+idx];
      companion[i].pbdot         = parfile->parameter[TEMPO3_PARAM_PBDOT+idx];
      companion[i].a1dot         = parfile->parameter[TEMPO3_PARAM_A1DOT+idx];
      companion[i].gamma         = parfile->parameter[TEMPO3_PARAM_GAMMA+idx];
    }
    if(whichcompanion >= 0) {
      for(i = 1; i < TEMPO3_MaxNrCompanions; i++) {
	companion[i].defined = 0;
      }
    }
    for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
      if(companion[i].defined != 0) {
	nrCompanions ++;
      }
    }    
  }


  torb = 0;
  if(parfile->binarymodel == 1) {  // BT model
    // bbat = mjd as a bat (in days)
    // t0 = Epoch of periastron (MJD)
    // pb = Orbital period (days)
    // ecc = Orbital eccentricity
    // pbdot = The first time derivative of binary period, same unit as tempo2 (s/s?)
    // a1 = Projected semimajor axis of orbit (lt-s)
    // om = Longitude of periastron (deg)
    // omdot = Periastron advance (deg/yr) 
    // a1dot = Rate of change of semimajor axis (lt-s / s) 
    // gamma = Post-Keplerian "gamma" term (s) 
    for(i = 1; i < TEMPO3_MaxNrCompanions; i++) {
      if(companion[i].defined != 0) {
	printwarning(0, "WARNING ephemeris_get_torb: Only one companion is allowed in the BT model, ignoring this companion from now on.");
	companion[i].defined = 0;
      }
    }
      
    if(calc_deriv > 0) {
      printerror(0, "ERROR ephemeris_get_torb: This binary model does currently not support the calculation of derivatives.");
      exit(0);
    }


    torb = BTmodel(mjd/TEMPO3_SecPerDay, companion[0]);
  }else if(parfile->binarymodel == 2) {
    if(calc_deriv > 0) {
      printerror(0, "ERROR ephemeris_get_torb: This binary model does currently not support the calculation of derivatives.");
      exit(0);
    }
    //      printf("XXXXX %Lf %Lf\n", companion[0].t0, companion[1].t0);
    //      fprintf(stderr, "XXXXXX entering T2model()\n");
    torb = T2model(mjd/TEMPO3_SecPerDay, nrCompanions, companion);
    //      fprintf(stderr, "XXXXXX exiting T2model()\n");
  }else if(parfile->binarymodel == 3) {  // HYPER model
    // bbat = mjd as a bat (in days)
    // t0 = Epoch of periastron (MJD)
    // pb = Orbital period (days)
    // ecc = Orbital eccentricity
    // pbdot = The first time derivative of binary period, same unit as tempo2 (s/s?)
    // a1 = Projected semimajor axis of orbit (lt-s)
    // om = Longitude of periastron (deg)
    // omdot = Periastron advance (deg/yr) 
    // a1dot = Rate of change of semimajor axis (lt-s / s) 
    // gamma = Post-Keplerian "gamma" term (s) 
    
    torb = T2model_hyper(mjd/TEMPO3_SecPerDay, nrCompanions, calc_deriv, -1, companion);
  }else if(parfile->binarymodel != 0) {
    printwarning(0, "WARNING ephemeris_get_torb(): Binary model not recognized.");
    torb = 0;
  }
  return torb;
}

// Should really only be used internally in library, but is declared in tempo3 as well.
// Note that the function definition is also appearing at the top of the library.
/* If mode = 0, it calculates the phase in turns
           = 1, it calculates the frequency (time derivative of phase)
           = 2, it calculates the frequency-devivative
           = 3, it calculates the second frequency-derivative

   For a given mjd, at a given freq in MHz (dm delay) evaluate the
   ephemeris provided with parfile (and paramset to identify which
   parameters are defined). ssb1, ssb2, ssb3 give the position of the
   barycentre at mjd. ephemeris_setLinkedValues() is assumed to be
   already applied.

   showWaves = 0 - all terms are included
   showWaves > 0 - showWaves is the first fitwave to be excluded from the sum
   showWaves < 0 - -showWaves it is the first fitwave to be included

   If noNu is nonzero, the phase increase because F0 is excluded from
   the calculation. I assume this is implemented for numerical
   derivatives or something? site and flags (for jumps) refer to other
   properties of the toa which can affect the ephemeris
   evaluation. dphase allows a phase offset to be applied (in mode ==
   0). If ignore_binary is non-zero, the binary effects are ignored in
   the calculation.

 */
long double evaluate_ephemeris(long double mjd, long double freq, long double ssb1, long double ssb2, long double ssb3, int showWaves, int noNu, char *site, char *flags, tempo3_parameters_def *parfile, long double dphase, int *paramset, int mode, int ignore_binary)
{
  long double phase, spinfreq, spinfreqderiv, spinfreqderiv2;
  long double dt, dtfitwaves, pepoch;
  int i;

  spinfreq = 0;
  spinfreqderiv = 0;
  spinfreqderiv2 = 0;


  /*
  printf("XXXXX F0 = %Le\n", parfile->parameter[TEMPO3_PARAM_F0+0]);
  printf("XXXXX F1 = %Le\n", parfile->parameter[TEMPO3_PARAM_F0+1]);
  printf("XXXXX F2 = %Le\n", parfile->parameter[TEMPO3_PARAM_F0+2]);
  printf("XXXXX F3 = %Le\n", parfile->parameter[TEMPO3_PARAM_F0+3]);
  printf("XXXXX F4 = %Le\n", parfile->parameter[TEMPO3_PARAM_F0+4]);
  printf("XXXXX F5 = %Le\n", parfile->parameter[TEMPO3_PARAM_F0+5]);
  printf("XXXXX F6 = %Le\n", parfile->parameter[TEMPO3_PARAM_F0+6]);
  printf("XXXXX F7 = %Le\n", parfile->parameter[TEMPO3_PARAM_F0+7]);
  printf("XXXXX F8 = %Le\n", parfile->parameter[TEMPO3_PARAM_F0+8]);
  exit(0);
  */

  /*  fprintf(stderr, "ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ 1 %p\n", paramlinked); */
  // NO LONGER DO IT HERE, AS IT IS A SIGNIFICANT OVERHEAD
  //  ephemeris_setLinkedValues_par_version(&parfile, paramlinked);
  /*  fprintf(stderr, "ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ 2\n"); */

  if(mode) {
    long double estimatedperiod = 1.0/ephemeris_estimateNu_numerically_avoid_using(mjd, freq, ssb1, ssb2, ssb3, site, flags, parfile, paramset);
    if(mode == 1)
      spinfreq += fitwave_function(parfile, mjd, mode, showWaves)/(TEMPO3_SecPerDay*estimatedperiod);
    else if(mode == 2)
      spinfreqderiv += fitwave_function(parfile, mjd, mode, showWaves)/(TEMPO3_SecPerDay*TEMPO3_SecPerDay*estimatedperiod);
    else if(mode == 3)
      spinfreqderiv2 += fitwave_function(parfile, mjd, mode, showWaves)/(TEMPO3_SecPerDay*TEMPO3_SecPerDay*TEMPO3_SecPerDay*estimatedperiod);
  }

  /* Always take time a bit later, because fitwaves function is in seconds rather than phase */
  /*  fprintf(stderr, "XXXXX mjd =%Lf\n", mjd); */
  if(parfile->nrwaves > 0) {
    dtfitwaves = fitwave_function(parfile, mjd, 0, showWaves);
    mjd += dtfitwaves/TEMPO3_SecPerDay; /* Divide by TEMPO3_SecPerDay as it is multiplied with same constant at the next line */
  }
  /*  fprintf(stderr, "XXXXX mjd2=%Lf\n", mjd); */
  /* Convert in seconds */
  mjd            *= TEMPO3_SecPerDay; 
  pepoch  = parfile->parameter[TEMPO3_PARAM_PEPOCH]  * TEMPO3_SecPerDay;
  long double *glep, *gltd;
  if(parfile->nrglitches > 0) {
    glep = malloc(parfile->nrglitches*sizeof(long double));
    gltd = malloc(parfile->nrglitches*sizeof(long double));
    if(glep == NULL || gltd == NULL) {
      printerror(0, "ERROR evaluate_ephemeris: Cannot allocate memory\n");
      exit(0);
    }
    /*  fprintf(stderr, "XXXXX mjd3=%Lf\n", mjd); */
    for(i = 0; i < parfile->nrglitches; i++) {
      glep[i] = parfile->parameter[TEMPO3_PARAM_GLEP+i]  * TEMPO3_SecPerDay;
      gltd[i] = parfile->parameter[TEMPO3_PARAM_GLTD+i]  * TEMPO3_SecPerDay; 
    }
    
  }


  /* DM correction not necessary, because barycentred toa's are
     already corrected for this. For DM fitting, the fitted dm is the
     difference with the original DM from the par file.*/
  if(parfile->parameter[TEMPO3_PARAM_DM] != 0.0) {
    // The following should be equivalent of calling calcDMDelay(), but faster? Not really it seems, so left as it was.
    // dt = DM_CONST*parfile->parameter[TEMPO3_PARAM_DM]/(freq*freq);
    long double t_minus_tepoch, t_minus_tepoch2, t_minus_tepoch3, t_minus_tepoch4, t_minus_tepoch5;
    t_minus_tepoch = (mjd/TEMPO3_SecPerDay-parfile->parameter[TEMPO3_PARAM_DMEPOCH])/365.25;
    t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
    t_minus_tepoch3 = t_minus_tepoch2*t_minus_tepoch;
    t_minus_tepoch4 = t_minus_tepoch2*t_minus_tepoch2;
    t_minus_tepoch5 = t_minus_tepoch3*t_minus_tepoch2;
    long double dm;
    dm = parfile->parameter[TEMPO3_PARAM_DM];
    if(TEMPO3_MaxNrDMderivatives >= 1) {
      if(parfile->parameter[TEMPO3_PARAM_DM+1] != 0) {
	dm += parfile->parameter[TEMPO3_PARAM_DM+1]*t_minus_tepoch;
      }
    }
    if(TEMPO3_MaxNrDMderivatives >= 2) {
      if(parfile->parameter[TEMPO3_PARAM_DM+2] != 0) {
	dm += parfile->parameter[TEMPO3_PARAM_DM+2]*t_minus_tepoch2/2.0;   // 2!
      }
    }
    if(TEMPO3_MaxNrDMderivatives >= 3) {
      if(parfile->parameter[TEMPO3_PARAM_DM+3] != 0) {
	dm += parfile->parameter[TEMPO3_PARAM_DM+3]*t_minus_tepoch3/6.0;   // 3!
      }
    }
    if(TEMPO3_MaxNrDMderivatives >= 4) {
      if(parfile->parameter[TEMPO3_PARAM_DM+4] != 0) {
	dm += parfile->parameter[TEMPO3_PARAM_DM+4]*t_minus_tepoch4/24.0;   // 4!
      }
    }
    if(TEMPO3_MaxNrDMderivatives >= 5) {
      if(parfile->parameter[TEMPO3_PARAM_DM+5] != 0) {
	dm += parfile->parameter[TEMPO3_PARAM_DM+5]*t_minus_tepoch5/120.0;  // 5!
      }
    }
    if(TEMPO3_MaxNrDMderivatives >= 6) {
      if(parfile->parameter[TEMPO3_PARAM_DM+6] != 0) {
	dm += parfile->parameter[TEMPO3_PARAM_DM+6]*t_minus_tepoch5*t_minus_tepoch/720.0;  // 6!
      }
    }
    if(TEMPO3_MaxNrDMderivatives >= 7) {
      if(parfile->parameter[TEMPO3_PARAM_DM+7] != 0) {
	dm += parfile->parameter[TEMPO3_PARAM_DM+7]*t_minus_tepoch5*t_minus_tepoch2/5040.0;  // 7!
      }
    }
    if(TEMPO3_MaxNrDMderivatives >= 8) {
      if(parfile->parameter[TEMPO3_PARAM_DM+8] != 0) {
	dm += parfile->parameter[TEMPO3_PARAM_DM+8]*t_minus_tepoch5*t_minus_tepoch3/40320.0;  // 8!
      }
    }
    if(TEMPO3_MaxNrDMderivatives >= 9) {
      if(parfile->parameter[TEMPO3_PARAM_DM+9] != 0) {
	dm += parfile->parameter[TEMPO3_PARAM_DM+9]*t_minus_tepoch5*t_minus_tepoch4/362880.0;  // 9!
      }
    }
    if(TEMPO3_MaxNrDMderivatives >= 10) {
      printerror(0, "ERROR tempo3: there are more DM derivatives defined than implemented\n");
      exit(0);
    }
    dt = calcDMDelay(freq, -1, 1, dm);
    mjd -= dt;
  }

  /* Correct for jump if flag is set. */
  //  printf("TOA %Lf APPLIED JUMPS:", mjd/TEMPO3_SecPerDay);
  for(i = 0; i < parfile->nrjumps; i++) {
    if(parfile->jumpmeaning[i] == 0) {  // Search for flag
      if(strstr(flags, parfile->jumpflag[i]) != NULL) {
	mjd += parfile->parameter[TEMPO3_PARAM_JUMPS+i];
	//	printf(" %Lf (flag=%s)", parfile->parameter[TEMPO3_PARAM_JUMPS+i], parfile->jumpflag[i]);
      }
    }else if(parfile->jumpmeaning[i] == 1) {  // Search for site
      if(strcmp(site, parfile->jumpflag[i]) == 0) {
	mjd -= parfile->parameter[TEMPO3_PARAM_JUMPS+i];  // CONFUSINGLY, THE TELESCOPE JUMP IS SUBTRACTED, BUT THE FLAG JUMP IS ADDED?????
	//	printf(" %Lf (site=%s)", parfile->parameter[TEMPO3_PARAM_JUMPS+i], parfile->jumpflag[i]);
	//      }else {
	//	printf("site=%s (%s) rejected", site, parfile->jumpflag[i]);	
      }
    }
  }
  //  printf("\n");

  if(parfile->parameter[TEMPO3_PARAM_RAJ] != 0.0 || parfile->parameter[TEMPO3_PARAM_DECJ] != 0.0 || parfile->parameter[TEMPO3_PARAM_PMRA] != 0.0 || parfile->parameter[TEMPO3_PARAM_PMDEC] != 0.0) {
    if(fabsl(ssb1) > 1 || fabsl(ssb2) > 1 || fabsl(ssb3) > 1) {
      long double x[3], ra, dec, constant;

      // This bit only needs to be once: undo the barycentre/proper motion bit. It is time dependent, so it needs to be done for each toa. In principle the mjd is changed by jumps etc. I assume those effects can be neglected for any reasonable situation.
      ra = parfile->raj;
      dec = parfile->decj;
      if(parfile->parameter[TEMPO3_PARAM_PMRA] != 0.0 || parfile->parameter[TEMPO3_PARAM_PMDEC] != 0.0 || parfile->pmra != 0.0 || parfile->pmdec != 0) {
	// pmra [mas/yr] * dt [yr] * (1e-3 / 3600.0) [deg/mas]
      	constant = ((mjd/TEMPO3_SecPerDay - parfile->parameter[TEMPO3_PARAM_POSEPOCH])/365.25 )*(1e-3/3600.0)*(M_PI/180.0);
      	ra  += parfile->pmra  * constant;
      	dec += parfile->pmdec * constant;
      }

      x[0] = cosl(dec)*cosl(ra);
      x[1] = cosl(dec)*sinl(ra);
      x[2] = sinl(dec);

      dt = ssb1*x[0] + ssb2*x[1] + ssb3*x[2];

      ra  += parfile->parameter[TEMPO3_PARAM_RAJ];
      dec += parfile->parameter[TEMPO3_PARAM_DECJ];
      if(parfile->parameter[TEMPO3_PARAM_PMRA] != 0.0 || parfile->parameter[TEMPO3_PARAM_PMDEC] != 0.0 || parfile->pmra != 0.0 || parfile->pmdec != 0) {
      	ra  += parfile->parameter[TEMPO3_PARAM_PMRA]  * constant;
      	dec += parfile->parameter[TEMPO3_PARAM_PMDEC] * constant;
      }

      x[0] = cosl(dec)*cosl(ra);
      x[1] = cosl(dec)*sinl(ra);
      x[2] = sinl(dec);

      dt -= ssb1*x[0] + ssb2*x[1] + ssb3*x[2];

      mjd -= dt;
    }else {
      printwarning(0, "WARNING: Position Earth w.r.t. barycenter not known, cannot fit for position/proper motion.");
    }
  }

  if(ignore_binary == 0) {
    mjd += ephemeris_get_torb(mjd, parfile, paramset, 0, -1);
  }

  /* Phase is the predicted pulse phase for the TOA */
  phase = parfile->parameter[TEMPO3_PARAM_PHASE0];

  //  fprintf(stderr, "YYYY1 phase=%Le\n", phase);

  phase += dphase;
  //  fprintf(stderr, "YYYY2 phase=%Le\n", phase);

  long double t_minus_tepoch, t_minus_tepoch2, t_minus_tepoch3, t_minus_tepoch4, t_minus_tepoch5, t_minus_tepoch10;
  t_minus_tepoch = mjd-pepoch;
  t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
  t_minus_tepoch3 = t_minus_tepoch2*t_minus_tepoch;
  t_minus_tepoch4 = t_minus_tepoch2*t_minus_tepoch2;
  t_minus_tepoch5 = t_minus_tepoch3*t_minus_tepoch2;
  t_minus_tepoch10 = t_minus_tepoch5*t_minus_tepoch5;
  
  // Not sure something like the following is effective to avoid unecessary calculations above.
  //  int t_minus_tepoch5_defined;
  //  for(i = TEMPO3_MaxNrFderivatives; i > 2; i--) {
  //    if(parfile->parameter[TEMPO3_PARAM_F0+i] != 0) {
  //      break;
  //    }
  //  }
  //  if(i >= 5) {
  //    t_minus_tepoch5_defined = 1;
  //  }


  // If no braking index, take into account F0 and F1 contribution. Otherwise another equation is used to calculate these constributions
  if(parfile->parameter[TEMPO3_PARAM_NBRAKE] == 0) {
    /* Spin frequency contribution */
    if(noNu == 0) {
      if(mode == 0) {
	//	fprintf(stderr, "YYYY3 phase=%Le  (noNu=%d F0=%Lf t_minus_tepoch=%Lf)\n", phase, noNu, parfile->parameter[TEMPO3_PARAM_F0], t_minus_tepoch);
	phase += parfile->parameter[TEMPO3_PARAM_F0]*t_minus_tepoch;
      }else if(mode == 1)
	spinfreq += parfile->parameter[TEMPO3_PARAM_F0];
      else if(mode == 2)
	spinfreqderiv += 0;
      else if(mode == 3)
	spinfreqderiv2 += 0;
    }else {
      if(parfile->nrwaves > 0) {
	phase += parfile->parameter[TEMPO3_PARAM_F0]*dtfitwaves;   // If you don't include the spin-frequency contribution itself, make sure that at least the fitwave contribution is included.
      }
    }
    //    fprintf(stderr, "YYYY3 phase=%Le\n", phase);
    /* Spin frequency derivative contribution */
    if(mode == 0)
      phase += 0.5*parfile->parameter[TEMPO3_PARAM_F0+1]*t_minus_tepoch2;
    else if(mode == 1)
      spinfreq += parfile->parameter[TEMPO3_PARAM_F0+1]*t_minus_tepoch;
    else if(mode == 2)
      spinfreqderiv += parfile->parameter[TEMPO3_PARAM_F0+1];
    else if(mode == 3)
      spinfreqderiv2 += 0;
  }else { // Braking index is specified
    if(parfile->parameter[TEMPO3_PARAM_F0+2] != 0) {
      printwarning(0, "WARNING tempo3: Setting breaking index and F2 at the same time is not likely to be sensible");
    }
    // When using braking index parameter, first undo the effect of F0 and F1, as they are both included in the equation below. NO LONGER NECESSARY.
    //    printf("XXXXX       = %Le\n", phase);
    /*    if(mode == 0) {
      phase -= parfile->parameter[TEMPO3_PARAM_F0]*(mjd-pepoch);
      phase -= parfile->parameter[TEMPO3_PARAM_F0+1]*(mjd-pepoch)*(mjd-pepoch)/2.0;
    }else if(mode == 1) {
      spinfreq -= parfile->parameter[TEMPO3_PARAM_F0];
      spinfreq -= parfile->parameter[TEMPO3_PARAM_F0+1]*(mjd-pepoch);
    }else if(mode == 2) {
      spinfreqderiv -= parfile->parameter[TEMPO3_PARAM_F0+1];
    }
    */
    //    printf("XXXXX       = %Le\n", phase);
    long double expr1, expr2;
    //        printf("XXXXX       nbreak =%Le\n", parfile->parameter[TEMPO3_PARAM_NBRAKE]);
    //    printf("XXXXX       mjd    =%Le\n", mjd);
    //    printf("XXXXX       F0     =%Le\n", parfile->parameter[TEMPO3_PARAM_F0]);
    //    printf("XXXXX       F1     =%Le\n", parfile->parameter[TEMPO3_PARAM_F0+1]);
    //    printf("XXXXX       pepoch =%Le\n", pepoch);
    long double f0_pwr_n = powl(parfile->parameter[TEMPO3_PARAM_F0], parfile->parameter[TEMPO3_PARAM_NBRAKE]);
    //    expr2 = -parfile->parameter[TEMPO3_PARAM_F0+1]/f0_pwr_n;
    //    expr2 *= (parfile->parameter[TEMPO3_PARAM_NBRAKE]-1.0)*t_minus_tepoch;
    //    //    expr2 += powl(parfile->parameter[TEMPO3_PARAM_F0], 1.0-parfile->parameter[TEMPO3_PARAM_NBRAKE]);
    //    expr2 += parfile->parameter[TEMPO3_PARAM_F0]/f0_pwr_n;
    expr2 = parfile->parameter[TEMPO3_PARAM_F0]-parfile->parameter[TEMPO3_PARAM_F0+1]*(parfile->parameter[TEMPO3_PARAM_NBRAKE]-1.0)*t_minus_tepoch;
    expr2 /= f0_pwr_n;

    if(mode == 0) {
      //      expr1 = powl(parfile->parameter[TEMPO3_PARAM_F0], parfile->parameter[TEMPO3_PARAM_NBRAKE])/((2.0-parfile->parameter[TEMPO3_PARAM_NBRAKE])*parfile->parameter[TEMPO3_PARAM_F0+1]);
      expr1 = f0_pwr_n/((2.0-parfile->parameter[TEMPO3_PARAM_NBRAKE])*parfile->parameter[TEMPO3_PARAM_F0+1]);
      phase += expr1*powl(expr2,(2.0-parfile->parameter[TEMPO3_PARAM_NBRAKE])/(1.0-parfile->parameter[TEMPO3_PARAM_NBRAKE]));
      phase -= parfile->parameter[TEMPO3_PARAM_F0]*parfile->parameter[TEMPO3_PARAM_F0]/(parfile->parameter[TEMPO3_PARAM_F0+1]*(2.0-parfile->parameter[TEMPO3_PARAM_NBRAKE]));
    }else if(mode == 1) {
      expr1 = 1.0;
      spinfreq += expr1*powl(expr2,1.0/(1.0-parfile->parameter[TEMPO3_PARAM_NBRAKE]));
    }else if(mode == 2) {
      expr1 = parfile->parameter[TEMPO3_PARAM_F0+1]/powl(parfile->parameter[TEMPO3_PARAM_F0], parfile->parameter[TEMPO3_PARAM_NBRAKE]);
      spinfreqderiv += expr1*powl(expr2,parfile->parameter[TEMPO3_PARAM_NBRAKE]/(1.0-parfile->parameter[TEMPO3_PARAM_NBRAKE]));
    }else if(mode == 3) {
      expr1 = parfile->parameter[TEMPO3_PARAM_NBRAKE]*parfile->parameter[TEMPO3_PARAM_F0+1]*parfile->parameter[TEMPO3_PARAM_F0+1]/powl(parfile->parameter[TEMPO3_PARAM_F0], 2.0*parfile->parameter[TEMPO3_PARAM_NBRAKE]);
      spinfreqderiv2 += expr1*powl(expr2,(2.0*parfile->parameter[TEMPO3_PARAM_NBRAKE]-1.0)/(1.0-parfile->parameter[TEMPO3_PARAM_NBRAKE]));
    }
    //    printf("XXXXX       = %Le\n", phase);
  }  // End of braking index or F0/F1 contribution


  /* F2 constribution */
  if(TEMPO3_MaxNrFderivatives >= 2) {
    if(parfile->parameter[TEMPO3_PARAM_F0+2] != 0) {
      if(mode == 0)
	phase += parfile->parameter[TEMPO3_PARAM_F0+2]*t_minus_tepoch3/6.0;
      else if(mode == 1)
	spinfreq += parfile->parameter[TEMPO3_PARAM_F0+2]*t_minus_tepoch2/2.0;
      else if(mode == 2)
	spinfreqderiv += parfile->parameter[TEMPO3_PARAM_F0+2]*t_minus_tepoch;
      else if(mode == 3)
	spinfreqderiv2 += parfile->parameter[TEMPO3_PARAM_F0+2];
    }
  }
  /* F3 constribution */
  if(TEMPO3_MaxNrFderivatives >= 3) {
    if(parfile->parameter[TEMPO3_PARAM_F0+3] != 0) {
      if(mode == 0)
	phase += parfile->parameter[TEMPO3_PARAM_F0+3]*t_minus_tepoch4/24.0;
      else if(mode == 1)
	spinfreq += parfile->parameter[TEMPO3_PARAM_F0+3]*t_minus_tepoch3/6.0;
      else if(mode == 2)
	spinfreqderiv += parfile->parameter[TEMPO3_PARAM_F0+3]*t_minus_tepoch2/2.0;
      else if(mode == 3)
	spinfreqderiv2 += parfile->parameter[TEMPO3_PARAM_F0+3]*t_minus_tepoch;
    }
  }
  /* F4 constribution */
  if(TEMPO3_MaxNrFderivatives >= 4) {
    if(parfile->parameter[TEMPO3_PARAM_F0+4] != 0) {
      if(mode == 0)
	phase += parfile->parameter[TEMPO3_PARAM_F0+4]*t_minus_tepoch5/120.0;
      else if(mode == 1)
	spinfreq += parfile->parameter[TEMPO3_PARAM_F0+4]*t_minus_tepoch4/24.0;
      else if(mode == 2)
	spinfreqderiv += parfile->parameter[TEMPO3_PARAM_F0+4]*t_minus_tepoch3/6.0;
      else if(mode == 3)
	spinfreqderiv2 += parfile->parameter[TEMPO3_PARAM_F0+4]*t_minus_tepoch2/2.0;
    }
  }
  /* F5 constribution */
  if(TEMPO3_MaxNrFderivatives >= 5) {
    if(parfile->parameter[TEMPO3_PARAM_F0+5] != 0) {
      if(mode == 0)
	phase += parfile->parameter[TEMPO3_PARAM_F0+5]*t_minus_tepoch5*t_minus_tepoch/720.0;
      else if(mode == 1)
	spinfreq += parfile->parameter[TEMPO3_PARAM_F0+5]*t_minus_tepoch5/120.0;
      else if(mode == 2)
	spinfreqderiv += parfile->parameter[TEMPO3_PARAM_F0+5]*t_minus_tepoch4/24.0;
      else if(mode == 3)
	spinfreqderiv2 += parfile->parameter[TEMPO3_PARAM_F0+5]*t_minus_tepoch3/6.0;
    }
  }
  /* F6 constribution */
  if(TEMPO3_MaxNrFderivatives >= 6) {
    if(parfile->parameter[TEMPO3_PARAM_F0+6] != 0) {
      if(mode == 0)
	phase += parfile->parameter[TEMPO3_PARAM_F0+6]*t_minus_tepoch5*t_minus_tepoch2/5040.0;
      else if(mode == 1)
	spinfreq += parfile->parameter[TEMPO3_PARAM_F0+6]*t_minus_tepoch5*t_minus_tepoch/720.0;
      else if(mode == 2)
	spinfreqderiv += parfile->parameter[TEMPO3_PARAM_F0+6]*t_minus_tepoch5/120.0;
      else if(mode == 3)
	spinfreqderiv2 += parfile->parameter[TEMPO3_PARAM_F0+6]*t_minus_tepoch4/24.0;
    }
  }
  /* F7 constribution */
  if(TEMPO3_MaxNrFderivatives >= 7) {
    if(parfile->parameter[TEMPO3_PARAM_F0+7] != 0) {
      if(mode == 0)
	phase += parfile->parameter[TEMPO3_PARAM_F0+7]*t_minus_tepoch5*t_minus_tepoch3/40320.0;
      else if(mode == 1)
	spinfreq += parfile->parameter[TEMPO3_PARAM_F0+7]*t_minus_tepoch5*t_minus_tepoch2/5040.0;
      else if(mode == 2)
	spinfreqderiv += parfile->parameter[TEMPO3_PARAM_F0+7]*t_minus_tepoch5*t_minus_tepoch/720.0;
      else if(mode == 3)
	spinfreqderiv2 += parfile->parameter[TEMPO3_PARAM_F0+7]*t_minus_tepoch5/120.0;
    }
  }
  /* F8 constribution */
  if(TEMPO3_MaxNrFderivatives >= 8) {
    if(parfile->parameter[TEMPO3_PARAM_F0+8] != 0) {
      if(mode == 0)
	phase += parfile->parameter[TEMPO3_PARAM_F0+8]*t_minus_tepoch5*t_minus_tepoch4/362880.0;
      else if(mode == 1)
	spinfreq += parfile->parameter[TEMPO3_PARAM_F0+8]*t_minus_tepoch5*t_minus_tepoch3/40320.0;
      else if(mode == 2)
	spinfreqderiv += parfile->parameter[TEMPO3_PARAM_F0+8]*t_minus_tepoch5*t_minus_tepoch2/5040.0;
      else if(mode == 3)
	spinfreqderiv2 += parfile->parameter[TEMPO3_PARAM_F0+8]*t_minus_tepoch5*t_minus_tepoch/720.0;
    }
  }
  /* F9 constribution */
  if(TEMPO3_MaxNrFderivatives >= 9) {
    if(parfile->parameter[TEMPO3_PARAM_F0+9] != 0) {
      if(mode == 0)
	phase += parfile->parameter[TEMPO3_PARAM_F0+9]*t_minus_tepoch5*t_minus_tepoch5/3628800.0; //10!
      else if(mode == 1)
	spinfreq += parfile->parameter[TEMPO3_PARAM_F0+9]*t_minus_tepoch5*t_minus_tepoch4/362880.0;  //9!
      else if(mode == 2)
	spinfreqderiv += parfile->parameter[TEMPO3_PARAM_F0+9]*t_minus_tepoch5*t_minus_tepoch3/40320.0; //8!
      else if(mode == 3)
	spinfreqderiv2 += parfile->parameter[TEMPO3_PARAM_F0+9]*t_minus_tepoch5*t_minus_tepoch2/5040.0; //7!
    }
  }
  /* F10 constribution */
  if(TEMPO3_MaxNrFderivatives >= 10) {
    if(parfile->parameter[TEMPO3_PARAM_F0+10] != 0) {
      if(mode == 0)
	phase += parfile->parameter[TEMPO3_PARAM_F0+10]*t_minus_tepoch10*t_minus_tepoch/39916800.0; //11!
      else if(mode == 1)
	spinfreq += parfile->parameter[TEMPO3_PARAM_F0+10]*t_minus_tepoch10/3628800.0;  //10!
      else if(mode == 2)
	spinfreqderiv += parfile->parameter[TEMPO3_PARAM_F0+10]*t_minus_tepoch5*t_minus_tepoch4/362880.0; //9!
      else if(mode == 3)
	spinfreqderiv2 += parfile->parameter[TEMPO3_PARAM_F0+10]*t_minus_tepoch5*t_minus_tepoch3/40320.0; //8!
    }
  }
  // F11 constribution
  if(TEMPO3_MaxNrFderivatives >= 11) {
    if(parfile->parameter[TEMPO3_PARAM_F0+11] != 0) {
      if(mode == 0)
	phase += parfile->parameter[TEMPO3_PARAM_F0+11]*t_minus_tepoch10*t_minus_tepoch2/479001600.0; //12!
      else if(mode == 1)
	spinfreq += parfile->parameter[TEMPO3_PARAM_F0+11]*t_minus_tepoch10*t_minus_tepoch/39916800.0;  //11!
      else if(mode == 2)
	spinfreqderiv += parfile->parameter[TEMPO3_PARAM_F0+11]*t_minus_tepoch10/3628800.0; //10!
      else if(mode == 3)
	spinfreqderiv2 += parfile->parameter[TEMPO3_PARAM_F0+11]*t_minus_tepoch5*t_minus_tepoch4/362880.0; //9!
    }
  }
  // F12 constribution
  if(TEMPO3_MaxNrFderivatives >= 12) {
    if(parfile->parameter[TEMPO3_PARAM_F0+12] != 0) {
      if(mode == 0)
	phase += parfile->parameter[TEMPO3_PARAM_F0+12]*t_minus_tepoch10*t_minus_tepoch3/6227020800.0; //13!
      else if(mode == 1)
	spinfreq += parfile->parameter[TEMPO3_PARAM_F0+12]*t_minus_tepoch10*t_minus_tepoch2/479001600.0;  //12!
      else if(mode == 2)
	spinfreqderiv += parfile->parameter[TEMPO3_PARAM_F0+12]*t_minus_tepoch10*t_minus_tepoch/39916800.0; //11!
      else if(mode == 3)
	spinfreqderiv2 += parfile->parameter[TEMPO3_PARAM_F0+12]*t_minus_tepoch10/3628800.0; //10!
    }
    //    printf("XXXXX phase = %Le\n", parfile->parameter[TEMPO3_PARAM_F0+12]);
  }
  if(TEMPO3_MaxNrFderivatives >= 13) {
    printerror(0, "ERROR tempo3: there are more spin-frequency derivatives defined than implemented\n");
    exit(0);    
  }
  //  printf("XXXXX phase = %Le\n", phase);


  if(parfile->nrglitches != 0) {
    for(i = 0; i < parfile->nrglitches; i++) {
      if(mjd > glep[i]) {
	phase += parfile->parameter[TEMPO3_PARAM_GLPH+i];

	long double t_minus_tglep, t_minus_tglep2, t_minus_tglep3;
	t_minus_tglep = mjd-glep[i];
	t_minus_tglep2 = t_minus_tglep*t_minus_tglep;
	t_minus_tglep3 = t_minus_tglep2*t_minus_tglep;

	/* GLF0 contribution */
	if(noNu == 0) {
	  if(mode == 0)
	    phase += parfile->parameter[TEMPO3_PARAM_GLF0+i]*t_minus_tglep;
	  else if(mode == 1)
	    spinfreq += parfile->parameter[TEMPO3_PARAM_GLF0+i];
	  else if(mode == 2)
	    spinfreqderiv += 0;
	  else if(mode == 3)
	    spinfreqderiv2 += 0;
	}
	/* GLF1 contrubution */
	if(mode == 0)
	  phase += parfile->parameter[TEMPO3_PARAM_GLF0+1*TEMPO3_MaxNrGlitches+i]*t_minus_tglep2/2.0;
	else if(mode == 1)
	  spinfreq += parfile->parameter[TEMPO3_PARAM_GLF0+1*TEMPO3_MaxNrGlitches+i]*t_minus_tglep;
	else if(mode == 2)
	  spinfreqderiv += parfile->parameter[TEMPO3_PARAM_GLF0+1*TEMPO3_MaxNrGlitches+i];
	else if(mode == 3)
	  spinfreqderiv2 += 0;

	/* GLF2 contrubution */
	if(TEMPO3_MaxNrGlitchFderivatives >= 2) {
	  if(parfile->parameter[TEMPO3_PARAM_GLF0+2*TEMPO3_MaxNrGlitches+i] != 0) {
	    if(mode == 0) 
	      phase += parfile->parameter[TEMPO3_PARAM_GLF0+2*TEMPO3_MaxNrGlitches+i]*t_minus_tglep3/6.0;
	    else if(mode == 1) 
	      spinfreq += parfile->parameter[TEMPO3_PARAM_GLF0+2*TEMPO3_MaxNrGlitches+i]*t_minus_tglep2/2.0;
	    else if(mode == 2) 
	      spinfreqderiv += parfile->parameter[TEMPO3_PARAM_GLF0+2*TEMPO3_MaxNrGlitches+i]*t_minus_tglep;
	    else if(mode == 3) 
	      spinfreqderiv2 += parfile->parameter[TEMPO3_PARAM_GLF0+2*TEMPO3_MaxNrGlitches+i];
	  }
	}
	/* GLF3 contrubution */
	if(TEMPO3_MaxNrGlitchFderivatives >= 3) {
	  if(parfile->parameter[TEMPO3_PARAM_GLF0+3*TEMPO3_MaxNrGlitches+i] != 0) {
	    if(mode == 0)
	      phase += parfile->parameter[TEMPO3_PARAM_GLF0+3*TEMPO3_MaxNrGlitches+i]*t_minus_tglep3*t_minus_tglep/24.0;
	    else if(mode == 1)
	      spinfreq += parfile->parameter[TEMPO3_PARAM_GLF0+3*TEMPO3_MaxNrGlitches+i]*t_minus_tglep3/6.0;
	    else if(mode == 2)
	      spinfreqderiv += parfile->parameter[TEMPO3_PARAM_GLF0+3*TEMPO3_MaxNrGlitches+i]*t_minus_tglep2/2.0;
	    else if(mode == 3)
	      spinfreqderiv2 += parfile->parameter[TEMPO3_PARAM_GLF0+3*TEMPO3_MaxNrGlitches+i]*t_minus_tglep;
	  }
	}
	/* GLF4 contrubution */
	if(TEMPO3_MaxNrGlitchFderivatives >= 4) {
	  if(parfile->parameter[TEMPO3_PARAM_GLF0+4*TEMPO3_MaxNrGlitches+i] != 0) {
	    if(mode == 0)
	      phase += parfile->parameter[TEMPO3_PARAM_GLF0+4*TEMPO3_MaxNrGlitches+i]*t_minus_tglep3*t_minus_tglep2/120.0;
	    else if(mode == 1)
	      spinfreq  += parfile->parameter[TEMPO3_PARAM_GLF0+4*TEMPO3_MaxNrGlitches+i]*t_minus_tglep3*t_minus_tglep/24.0;
	    else if(mode == 2)
	      spinfreqderiv  += parfile->parameter[TEMPO3_PARAM_GLF0+4*TEMPO3_MaxNrGlitches+i]*t_minus_tglep3/6.0;
	    else if(mode == 3)
	      spinfreqderiv2  += parfile->parameter[TEMPO3_PARAM_GLF0+4*TEMPO3_MaxNrGlitches+i]*t_minus_tglep2/2.0;
	  }
	}
	/* GLF5 contrubution */
	if(TEMPO3_MaxNrGlitchFderivatives >= 5) {
	  if(parfile->parameter[TEMPO3_PARAM_GLF0+5*TEMPO3_MaxNrGlitches+i] != 0) {
	    if(mode == 0)
	      phase += parfile->parameter[TEMPO3_PARAM_GLF0+5*TEMPO3_MaxNrGlitches+i]*t_minus_tglep3*t_minus_tglep3/720.0;
	    else if(mode == 1)
	      spinfreq  += parfile->parameter[TEMPO3_PARAM_GLF0+5*TEMPO3_MaxNrGlitches+i]*t_minus_tglep3*t_minus_tglep2/120.0;
	    else if(mode == 2)
	      spinfreqderiv  += parfile->parameter[TEMPO3_PARAM_GLF0+5*TEMPO3_MaxNrGlitches+i]*t_minus_tglep3*t_minus_tglep/24.0;
	    else if(mode == 3)
	      spinfreqderiv2  += parfile->parameter[TEMPO3_PARAM_GLF0+5*TEMPO3_MaxNrGlitches+i]*t_minus_tglep3/6.0;
	  }
	}
      
	if(glep[i] > 0) {
	  /* GLF0D contribution */
	  /* I wouldn't do 1 - exp, but just -exp, but that is how it works in tempo2. */
	  if(gltd[i] > 0) {
	    if(mode == 0)
	      phase += parfile->parameter[TEMPO3_PARAM_GLF0D+i]*(gltd[i])*(1.0-expl(-t_minus_tglep/gltd[i]));
	    else if(mode == 1)
	      spinfreq += parfile->parameter[TEMPO3_PARAM_GLF0D+i]*expl(-t_minus_tglep/gltd[i]);
	    else if(mode == 2) {
	      /*	    fprintf(stderr, "XXXXXX mjd=%Lf td=%Lf glep=%Lf %Lf\n", mjd, gltd[i], glep[i], spinfreqderiv); */
	      spinfreqderiv += -parfile->parameter[TEMPO3_PARAM_GLF0D+i]*expl(-t_minus_tglep/gltd[i])/gltd[i];
	      /*	    fprintf(stderr, "XXXXXX mjd=%Lf td=%Lf glep=%Lf %Lf\n", mjd, gltd[i], glep[i], spinfreqderiv); */
	    }else if(mode == 3) {
	      spinfreqderiv2 += parfile->parameter[TEMPO3_PARAM_GLF0D+i]*expl(-t_minus_tglep/gltd[i])/(gltd[i]*gltd[i]);
	    }
	  }

	  /* GLF1D contribution */
	  /*
	  if(gltd[i] > 0) {
	    if(mode == 0)
	      phase += gltd[i]*gltd[i]*parfile->parameter[TEMPO3_PARAM_GLF1D+i]*expl(-(mjd-glep[i])/gltd[i]);
	    else if(mode == 1)
	      spinfreq += -gltd[i]*parfile->parameter[TEMPO3_PARAM_GLF1D+i]*expl(-(mjd-glep[i])/gltd[i]);
	    else if(mode == 2)
	      spinfreqderiv += parfile->parameter[TEMPO3_PARAM_GLF1D+i]*expl(-(mjd-glep[i])/gltd[i]);
	    else if(mode == 3)
	      spinfreqderiv2 += -parfile->parameter[TEMPO3_PARAM_GLF1D+i]*expl(-(mjd-glep[i])/gltd[i])/gltd[i];
	      }*/
	}

	if(parfile->parameter[TEMPO3_PARAM_GLNBRAKE+i] != 0) {
	  printerror(0, "ERROR tempo3: the GLNBRAKE parameter is not (yet?) implemented\n");
	  exit(0);
	}	  
	
      }
    }
  }
  
  if(parfile->nrglitches > 0) {
    free(glep);
    free(gltd);
  }


  if(mode == 1)
    return spinfreq;
  else if(mode == 2)
    return spinfreqderiv;
  else if(mode == 3)
    return spinfreqderiv2;
  return phase;
}


// Allocates memory to hold nrtoas active toas plus nrtoasdeleted
// deleted toas. Which ones are deleted, or how many, can be changed
// later but the total is the maximum capacity. The TOA values are not
// set or initialised, except pointers to strings which are set to
// NULL. If ssb is nonzero, space is allocated to store errors solar
// system barycentre coordinates. If floatversion is nonzero, float
// versions of the mjd and frequency is allocated, which can be useful
// for plotting purposes for example. If freqevol is nonzero, there is
// enough space allocated to hold information about the frequency and
// frequency derivative as function of time. A call to free_toas()
// will release the memory.
//
// Return 0: error
// Return 1: ok
int initialise_toas(toa_data_def *toas, long nrtoas, long nrtoasdeleted, int ssb, int floatversion, int freqevol, verbose_definition verbose)
{
  long i, total;
  toas->nrtoas = nrtoas;
  toas->nrtoasdeleted = nrtoasdeleted;
  total = nrtoas + nrtoasdeleted;

  toas->mjd = malloc(total*sizeof(long double));
  toas->err = malloc(total*sizeof(float));
  toas->freqSite = malloc(total*sizeof(long double));
  toas->freqSSB = malloc(total*sizeof(long double));
  toas->filenames = malloc(total*sizeof(char *));
  toas->flags = malloc(total*sizeof(char *));
  toas->site = malloc(total*sizeof(char *));
  toas->deleted = malloc(total*sizeof(int));
  if(toas->mjd == NULL || toas->err == NULL || toas->freqSite == NULL || toas->freqSSB == NULL || toas->filenames == NULL || toas->flags == NULL  || toas->site == NULL || toas->deleted == NULL) {
    printerror(verbose.debug, "ERROR initialise_toas: Cannot allocate memory\n");
    return 0;
  }
  if(ssb) {
    toas->ssb1 = malloc(total*sizeof(long double));
    toas->ssb2 = malloc(total*sizeof(long double));
    toas->ssb3 = malloc(total*sizeof(long double));
    if(toas->ssb1 == NULL || toas->ssb2 == NULL || toas->ssb3 == NULL) {
      printerror(verbose.debug, "ERROR initialise_toas: Cannot allocate memory\n");
      return 0;
    }
  }else {
    toas->ssb1 = NULL;
    toas->ssb2 = NULL;
    toas->ssb3 = NULL;
  }
  if(floatversion) {
    toas->mjdx = malloc(total*sizeof(float));
    toas->freqx = malloc(total*sizeof(float));
    //    printf("XXXX freqx allocated %p\n", toas->freqx);
    if(toas->mjdx == NULL || toas->freqx == NULL) {
      printerror(verbose.debug, "ERROR initialise_toas: Cannot allocate memory\n");
      return 0;
    }
  }else {
    toas->mjdx  = NULL;
    toas->freqx = NULL;
  }
  if(freqevol) {
    toas->nu = malloc(total*sizeof(long double));
    toas->nudot = malloc(total*sizeof(long double));
    toas->nu_mjd = malloc(total*sizeof(long double));
    if(toas->nu == NULL || toas->nudot == NULL || toas->nu_mjd == NULL) {
      printerror(verbose.debug, "ERROR initialise_toas: Cannot allocate memory\n");
      return 0;
    }
  }else {
    toas->nu     = NULL;
    toas->nudot  = NULL;
    toas->nu_mjd = NULL;
  }
  for(i = 0; i < total; i++) {
    toas->filenames[i] = NULL;
    toas->flags[i] = NULL;
    toas->site[i] = NULL;
  }

  return 1;
}

// Release the memory previously allocated with initialise_toas()
void free_toas(toa_data_def *toas)
{
  long i;
  for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
    if(toas->filenames[i] != NULL) {
      free(toas->filenames[i]);
      toas->filenames[i] = NULL;
    }
    if(toas->flags[i] != NULL) {
      free(toas->flags[i]);
      toas->flags[i] = NULL;
    }
    if(toas->site[i] != NULL) {
      free(toas->site[i]);
      toas->site[i] = NULL;
    }
  }
  free(toas->mjd);
  free(toas->err);
  free(toas->freqSite);
  free(toas->freqSSB);
  free(toas->filenames);
  free(toas->flags);
  free(toas->site);
  free(toas->deleted);
  if(toas->ssb1 != NULL && toas->ssb2 != NULL && toas->ssb3 != NULL) {
    free(toas->ssb1);
    free(toas->ssb2);
    free(toas->ssb3);
  }
  if(toas->mjdx != NULL) {
    free(toas->mjdx);
    toas->mjdx = NULL;
  }
  if(toas->freqx != NULL) {
    //    printf("XXXX freqx freed %p\n", toas->freqx);
    free(toas->freqx);
    toas->freqx = NULL;
  }
  if(toas->nu != NULL) {
    free(toas->nu);
    toas->nu = NULL;
  }
  if(toas->nudot != NULL) {
    free(toas->nudot);
    toas->nudot = NULL;
  }
  if(toas->nu_mjd != NULL) {
    free(toas->nu_mjd);
    toas->nu_mjd = NULL;
  }
  toas->nrtoas = 0;
  toas->nrtoasdeleted = 0;
}

// Add a TOA to toas at position toa_index. Its properties are given
// by toa_mjd, toa_err (in microseconds), toa_freq , if deleted != 0
// then the TOA is flagged to be deleted, toa_filename, toa_flags,
// toa_site, toa_ssb (X,Y,Z = position earth w.r.t. solar system
// barycentre in lt-s). nrtoas and nrtoasdeleted are not updated.
//
// Return 0: error
// Return 1: ok
int add_toa(toa_data_def *toas, long toa_index, long double toa_mjd, float toa_err, long double toa_freq, long double toa_freqSSB, int deleted, char *toa_filename, char *toa_flags, char *toa_site, long double *toa_ssb, verbose_definition verbose)
{
  toas->mjd[toa_index] = toa_mjd;
  toas->err[toa_index] = toa_err;
  toas->freqSite[toa_index] = toa_freq;
  toas->freqSSB[toa_index] = toa_freqSSB;

  if(toas->mjdx != NULL) {
    toas->mjdx[toa_index] = toa_mjd;
  }
  if(toas->freqx != NULL) {
    toas->freqx[toa_index] = toa_freq;
  }

  if(toas->filenames[toa_index] != NULL) {
    free(toas->filenames[toa_index]);
  }
  if(toa_filename != NULL) {
    toas->filenames[toa_index] = malloc(strlen(toa_filename)+1);
    if(toas->filenames[toa_index] == NULL) {
      printerror(verbose.debug, "ERROR add_toa: Memory allocation error\n");
      return 0;
    }
    strcpy(toas->filenames[toa_index], toa_filename);
  }

  if(toas->flags[toa_index] != NULL) {
    free(toas->flags[toa_index]);
  }
  if(toa_flags != NULL) {
    toas->flags[toa_index] = malloc(strlen(toa_flags)+1);
    if(toas->flags[toa_index] == NULL) {
      printerror(verbose.debug, "ERROR add_toa: Cannot allocate memory\n");
      return 0;
    }
    strcpy(toas->flags[toa_index], toa_flags);
  }

  if(toas->site[toa_index] != NULL) {
    free(toas->site[toa_index]);
  }
  if(toa_site != NULL) {
    toas->site[toa_index] = malloc(strlen(toa_site)+1);
    if(toas->site[toa_index] == NULL) {
      printerror(verbose.debug, "ERROR add_toa: Cannot allocate memory\n");
      return 0;
    }
    strcpy(toas->site[toa_index], toa_site);
  }

  if(toa_ssb != NULL) {
    toas->ssb1[toa_index] = toa_ssb[0];
    toas->ssb2[toa_index] = toa_ssb[1];
    toas->ssb3[toa_index] = toa_ssb[2];
  }
  toas->deleted[toa_index] = deleted;
  return 1;
}


// Given flag, a newflag is (assumed to be able to hold at least
// TEMPO3_MaxFlagsLength + 1 bytes) get filled in with flag without
// the pulse number flag.
//
// Return 0: error
// Return 1: ok
int stripPNflag(char *flag, char *newflag, verbose_definition verbose)
{
  int nrwords, i, index;
  char *txt;
  char *tmptxt;
  //  verbose_definition noverbose;
  //  cleanVerboseState(&noverbose);
  if(verbose.debug == 0) {
    verbose.verbose = 0;
  }

  tmptxt = malloc(TEMPO3_MaxFlagsLength + 1);
  if(tmptxt == NULL) {
    printerror(verbose.debug, "ERROR stripPNflag: Memory allocation error");
    return 0;
  }

  newflag[0] = 0;
  txt = pickWordFromString(flag, 1, &nrwords, 1, ' ', verbose);
  //  printf("XXXX flag=%s (nrwords=%d)\n", flag, nrwords);
  index = 0;
  for(i = 1; i <= nrwords; i++) {
    txt = pickWordFromString(flag, i, &nrwords, 1, ' ', verbose);
    //    printf("XXXX pick=%s (i=%d/%d)\n", txt, i, nrwords);
    sscanf(txt, "%s", tmptxt);
    //    printf("XXXX word=%s\n", tmptxt);
    if(strcmp(tmptxt, "-pn") == 0) {
      i++;  // Also skip value
    }else {
      if(index)
	strcat(newflag, " ");
      strcat(newflag, tmptxt);
      index++;
      //      printf("XXXX newflag=%s\n", newflag);
    }
  }
  free(tmptxt);
  return 1;
}

// Write out the information from toas to a file with name
// filename. If nturn != NULL, the pulse numbering information is
// included, corrected for wraps if != NULL. If isbarycentred != 0,
// then all the TOA' are assumed to be defined at the barycentre and
// the site codes are written out accordingly independent of what the
// original site was set to.
//
// Return 0: error
// Return 1: ok
int writetimfile(char *filename, toa_data_def toas, long long *nturn, tempo3_wraps_def *wraps, int isbarycentred, verbose_definition verbose)
{
  long toaindex, wrapindex;
  char *newflag;
  FILE *fin;
  long long toa_nturn;

  fin = fopen(filename, "w");
  if(fin == NULL) {
    printerror(verbose.debug, "ERROR writetimfile: Cannot open %s", filename);
    return 0;
  }
  newflag = malloc(TEMPO3_MaxFlagsLength + 1);
  if(newflag == NULL) {
    printerror(verbose.debug, "ERROR writetimfile: Memory allocation error");
    return 0;
  }

  toa_nturn = 0;
  fprintf(fin, "FORMAT 1\n");
  for(toaindex = 0; toaindex < toas.nrtoas + toas.nrtoasdeleted; toaindex++) {
    if(nturn != NULL) {
      toa_nturn = nturn[toaindex];
    }
    if(wraps != NULL) {
      for(wrapindex = 0; wrapindex < wraps->nrwraps; wrapindex++) {
	if(toas.mjdx[toaindex] >= wraps->mjd[wrapindex]) {
	  toa_nturn -= wraps->phase[wrapindex];
	}
      }
    }
    if(toas.deleted[toaindex] == 0) {
      long double freq;
      if(toas.freqSSB[toaindex] > 0 && isbarycentred) {
	freq = toas.freqSSB[toaindex];
      }else {
	freq = toas.freqSite[toaindex];
      }
      //      printf("XXXXXX FREQ: %Lf %Lf %Lf\n", freq, toas.freqSSB[toaindex], toas.freqSite[toaindex]);
      fprintf(fin, " %s %15.8Lf  %.13Lf   %.2lf", toas.filenames[toaindex], freq, toas.mjd[toaindex], toas.err[toaindex]);
      if(isbarycentred) {
	fprintf(fin, "  %s", "@");
      }else {
	if(toas.site == NULL) {
	  printerror(verbose.debug, "ERROR writetimfile: Site information is not defined for TOA %ld", toaindex);
	  return 0;
	}
	if(toas.site[toaindex] == NULL) {
	  printerror(verbose.debug, "ERROR writetimfile: Site information is not defined for TOA %ld", toaindex);
	  return 0;
	}
	fprintf(fin, "  %s", toas.site[toaindex]);
      }
      if(nturn != NULL) {
	fprintf(fin, " -pn %lld", toa_nturn);   // Used to be %Ld, which is incorrect.
      }
      if(toas.ssb1 != NULL && toas.ssb2 != NULL && toas.ssb3 != NULL) {
	if(fabsl(toas.ssb1[toaindex]) > 1 || fabsl(toas.ssb2[toaindex]) > 1 || fabsl(toas.ssb3[toaindex]) > 1) {
	  fprintf(fin, " -ssb1 %.20Lf -ssb2 %.20Lf -ssb3 %.20Lf -freqSSB %.20Lf", toas.ssb1[toaindex], toas.ssb2[toaindex], toas.ssb3[toaindex], toas.freqSSB[toaindex]);	    
	}
      }
      if(toas.flags[toaindex] != NULL) {
	stripPNflag(toas.flags[toaindex], newflag, verbose);
	fprintf(fin, " %s", newflag);
      }
      fprintf(fin, "\n");
    }
  }
  fclose(fin);
  free(newflag);
  if(verbose.verbose) {
    printf("Writing done to %s\n", filename);
  }
  return 1;
}



// Should really only be used internally in library, but is declared in tempo3 as well.
// Note that the function definition is also appearing at the top of the library.

/* If mode = 0, it calculates the phase in turns
           = 1, it calculates the frequency (time derivative of phase)
           = 2, it calculates the frequency-devivative
           = 3, it calculates the second frequency-derivative

   For a given mjd, at a given freq in MHz (dm delay) evaluate the
   ephemeris provided with parfile (and paramset to identify which
   parameters are defined). ssb1, ssb2, ssb3 give the position of the
   barycentre at mjd. ephemeris_setLinkedValues() is assumed to be
   already applied.

   showWaves = 0 - all terms are included
   showWaves > 0 - showWaves is the first fitwave to be excluded from the sum
   showWaves < 0 - -showWaves it is the first fitwave to be included

   If noNu is nonzero, the phase increase because F0 is excluded from
   the calculation. I assume this is implemented for numerical
   derivatives or something? site and flags (for jumps) refer to other
   properties of the toa which can affect the ephemeris
   evaluation. dphase allows a phase offset to be applied (in mode ==
   0). If ignore_binary is non-zero, the binary effects are ignored in
   the calculation.

 */
// long double evaluate_ephemeris(long double mjd, long double freq, long double ssb1, long double ssb2, long double ssb3, int showWaves, int noNu, char *site, char *flags, tempo3_parameters_def *parfile, long double dphase, int *paramset, int mode, int ignore_binary)

// Evaluate the spin-frequency at the given barycentric mjd.
//
// Return 0: error
// Return 1: ok
int evaluate_ephemeris_spinfrequency(ephemeris_def eph, long double bary_mjd, long double *spinfreq, verbose_definition verbose)
{
  // Observing frequency of 1400 MHz is only used for DM delay, not important for spin-frequency
  // SSB position also not relevant for spin-frequency (assumed to be a barycentric mjd anyway)
  // Site/flags not important as well (jumps only)
  *spinfreq = evaluate_ephemeris(bary_mjd, 1400.0, 0.0, 0.0, 0.0, 0, 0, "notimportant", "notimportant", &(eph.par), 0.0, eph.paramset, 1, 0);
  return 1;
}


// Filter some non-important parameters out (such as hyper binary
// model parameters), which are not important when determining the SSB
// positions.
//
// Return 0: error
// Return 1: ok
int filter_parfile_for_tempo2(char *inputfile, char *outputfile, tempo3_parameters_def *parfile, verbose_definition verbose)
{
  char *line, *word;
  FILE *fin, *fout;
  fin = fopen(inputfile, "r");
  if(fin == NULL) {
    printerror(verbose.debug, "ERROR filter_parfile_for_tempo2: Cannot open %s for reading", inputfile);	
    return 0;
  }
  fout = fopen(outputfile, "w");
  if(fout == NULL) {
    printerror(verbose.debug, "ERROR filter_parfile_for_tempo2: Cannot open %s for writing", outputfile);	
    return 0;
  }
  if(verbose.verbose) {
    printf("Generating %s which is identical to %s\n", outputfile, inputfile);
  }
  line = malloc(TEMPO3_MaxNrParLineLength+1);
  word = malloc(TEMPO3_MaxNrParLineLength+1);
  if(line == NULL || word == NULL) {
    printerror(verbose.debug, "ERROR filter_parfile_for_tempo2: Memory allocation error");
    return 0;
  }
  while(fgets(line, TEMPO3_MaxNrParLineLength, fin) != NULL) {
    if(strlen(line) > 0) {
      int ret;
      ret = sscanf(line, "%s", word);
      if(ret == 1) {
	if(strlen(word) > 0) {
	  if(parfile->binarymodel == 3) {  // Need to filter out all binary parameters if the HYPER model is used
	    if(strncmp(word, "BINARY", 6) != 0 && strncmp(word, "PB", 2) != 0) {
	      fputs(line, fout);
	    }else if(verbose.verbose) {
	      printf("  except line: %s", line);
	    }
	  }else {
	    fputs(line, fout);
	  }
	}
      }
    }
  }
  fclose(fin);
  fclose(fout);
  free(line);
  free(word);
  return 1;
}

/*
  After running tempo2 with the general2 plugin to extract the relevant information from tempo2 about each toa, filter the output (only consider lines with the list of toa's, ignore all other statements made by tempo2. indx will be initialised with the a list. The first element gives the toa number (in input file) corresponding to the smallest mjd etc. This is used when linking the site code as read in separately to the correct toa.
 */
int filter_tempo2_general2_output(char *filename_in, char *filename_out, int nosort, unsigned long **indx, int repeatable, verbose_definition verbose)
{
  long nrlines, firsttime, i, curoutputline;
  char *line, **output;
  double *mjd;
  FILE *fin, *fout;
  fin = fopen(filename_in, "r");
  if(fin == NULL) {
    printerror(verbose.debug, "ERROR filter_tempo2_general2_output(): Cannot open %s\n", filename_in);
    return 0;
  }
  fout = fopen(filename_out, "w");
  if(fout == NULL) {
    printerror(verbose.debug, "ERROR filter_tempo2_general2_output(): Cannot open %s\n", filename_out);
    return 0;
  }
  line = malloc(TEMPO3_MaxNrParLineLength+1);
  if(line == NULL) {
    printerror(verbose.debug, "ERROR filter_tempo2_general2_output(): Memory allocation error\n");
    return 0;
  }


  // Determine the number of lines
  nrlines = 0;
  while(fgets(line, TEMPO3_MaxNrParLineLength, fin) != NULL) {
    //    if(nrlines == 0) {
    //    fprintf(stderr, "First line of tim file: %s", line);
    //    }
    if(strncmp(line, "grepthislineout", 15) == 0) {
      nrlines++;
    }else if(verbose.verbose) {
      //repeatable && 
      if(strncmp(line, "Finishing off: time taken", 25) != 0) {
	if(strncmp(line, "*** Adding phase", 16) != 0 && strncmp(line, "Have 5", 6) != 0 && strncmp(line, "Have 4", 6) != 0 && strncmp(line, "WARNING [MISC1]: Unknown parameter in par file:  PHASE0", 55) != 0 && strncmp(line, "WARNING [MISC1]: Unknown parameter in par file:  NBRAKE", 55) != 0) {
	  // repeatable == 0 || 
	  if(verbose.verbose || (strncmp(line, "WARNING [PAR2]: Have not set a DM epoch", 39) != 0)) {
	    printf("%s", line);
	  }
	}
      }
    }
  }

  if(nrlines == 0) {
    printerror(verbose.debug, "ERROR filter_tempo2_general2_output(): No TOA's were found in %s. Tempo2 reports:\n", filename_in);
    rewind(fin);
    line[0] = ' ';
    line[1] = ' ';
    line[2] = 0;
    while(fgets(line+2, TEMPO3_MaxNrParLineLength, fin) != NULL) {
      line[strlen(line)-1] = 0;
      printerror(verbose.debug, line);
    }					     
    return 0;
  }

  // repeatable == 0
  if(verbose.verbose)
    printf("  Filtering %ld toa's from %s\n", nrlines, filename_in);
  output = malloc(nrlines*sizeof(char *)); // List of nrlines pointers that will point to the strings to be sent to output file
  *indx = (unsigned long *)malloc(nrlines*sizeof(unsigned long));
  mjd = (double *)malloc(nrlines*sizeof(double));
  if(output == NULL || *indx == NULL || mjd == NULL) {
    printerror(verbose.debug, "ERROR filter_tempo2_general2_output(): Memory allocation error\n");
    return 0;
  }

  rewind(fin);
  firsttime = 1;
  curoutputline = 0;
  while(fgets(line, 1000, fin) != NULL) {
    if(strncmp(line, "grepthislineout", 15) == 0) {  // Is a TOA to be stored
      output[curoutputline] = malloc(strlen(line)-14); // Don't store the "grepthislineout" part of line, but do store the line terminator
      if(output[curoutputline] == NULL) {
	printerror(verbose.debug, "ERROR filter_tempo2_general2_output(): Memory allocation error\n");
	return 0;
      }
      strcpy(output[curoutputline], line+15);
      curoutputline++;
    }else if(verbose.verbose) {  //  && repeatable == 0
      if(firsttime) {
	printf("    tempo2 reported:\n");
	firsttime = 0;
      }
      printf("      %s", line);
    }
  }

  if(nosort == 0) {
    printf("  Sorting %ld toa's\n", nrlines);
    for(i = 0; i < nrlines; i++) {
      if(sscanf(output[i], "%lf", &mjd[i]) != 1) {
	printerror(verbose.debug, "ERROR filter_tempo2_general2_output(): Cannot interpret first column as mjd's\n");
	return 0;	
      }
    }
    gsl_sort_index(*indx, mjd, 1, nrlines);
  }else {
    for(i = 0; i < nrlines; i++) {
      (*indx)[i] = i;
    }
  }


  for(i = 0; i < nrlines; i++) {
    if(fputs(output[(*indx)[i]], fout) == EOF) { 
      printerror(verbose.debug, "ERROR filter_tempo2_general2_output(): write error to %s\n", filename_out);
      return 0;
    }
    //    printf("TOA %ld: ", (*indx)[i]);
    //    puts(output[(*indx)[i]]);
  }


  fclose(fin);
  fclose(fout);
  free(line);
  for(i = 0; i < nrlines; i++)
    free(output[i]);
  free(output);
  free(mjd);
  return nrlines;
}

/*
  I THINK this is just used to read in the TOAs after being written
  out in the tempo3 style TIM file format as produced by the general2 plugin.

  If readturns is set, this will be the only things that is being read
  in, so do this after calling readtimfile_custom without this flag enabled
  (allowing you to allocate memory first for the residuals structure.

  If shortformat is set, the solar system barycentre coordinates are not expected to be set.
 */
int readtimfile_custom(char *filename, toa_data_def *toas, int shortformat, verbose_definition verbose)
{
  int nrlines, i, j;
  char *line, *toa_filename, *toa_flags;
  long double toa_mjd, toa_freq, toa_ssb[3], toa_freqSSB;
  float toa_err;

  if(verbose.verbose) {
    printf("Reading in %s\n", filename);
  }

  FILE *fin;
  fin = fopen(filename, "r");
  if(fin == NULL) {
    printerror(verbose.debug, "ERROR readtimfile_custom: Cannot open %s\n", filename);
    return 0;
  }
  line = malloc(10000);
  toa_filename = malloc(MaxFilenameLength+1);
  toa_flags = malloc(TEMPO3_MaxFlagsLength+1);
  if(line == NULL || toa_filename == NULL || toa_flags == NULL) {
    printerror(verbose.debug, "ERROR readtimfile_custom: Memory allocation error\n");
    return 0;
  }

  nrlines = 0;
  while(fgets(line, 1000, fin) != NULL) {
    nrlines++;
  }
  //  if(repeatable == 0 || verbose.debug) {
  if(verbose.verbose) {
    printf("  Reading in %d toa's from %s\n", nrlines, filename);
  }
  if(initialise_toas(toas, nrlines, 0, 1, 1, 0, verbose) == 0) {
    printerror(verbose.debug, "ERROR readtimfile_custom: Cannot initialise TOAs\n");
    free(line);
    free(toa_filename);
    free(toa_flags);
    return 0;
  }
  for(i = 0; i < nrlines; i++) {
    toas->deleted[i] = 0;
  }
  toas->nrtoasdeleted = 0;
  rewind(fin);

  for(i = 0; i < nrlines; i++) {
    if(fgets(line, 1000, fin) == NULL) {
      printerror(verbose.debug, "ERROR readtimfile_custom: Error reading line %d\n", i+1);
      free(line);
      free(toa_filename);
      free(toa_flags);
      return 0;
    }
    //    if(i == 0) {
    //      fprintf(stderr, "First line of tim file: %s", line);
    //      fprintf(stderr, "shortformat = %d\n", shortformat);
    //    }
    if(shortformat == 0) {
      if((j = sscanf(line, "%Lf %f %Lf %Lf %Lf %Lf %Lf %s", &toa_mjd, &toa_err, &(toa_ssb[0]), &(toa_ssb[1]), &(toa_ssb[2]), &toa_freq, &toa_freqSSB, toa_filename)) != 8) {
	printerror(verbose.debug, "ERROR readtimfile_custom: Error parsing line %d (only %d elements read in from '%s')\n", i+1, j, line);
	free(line);
	free(toa_filename);
	free(toa_flags);
	return 0;
      }
      //      printf("XXXXXX READ FREQ SSB: %Lf\n", toa_freqSSB);
      if(toa_freqSSB/toa_freq > 1000.0) {   // It turns out the general2 plugin writes out the SSB freq in Hz rather than MHz, UNLESS the TOA's are already barycentred, in which case general2 plugin writes the Site frequency in MHz.
	toa_freqSSB /= 1.0e6;
      }
    }else {
      if((j = sscanf(line, "%Lf %f %Lf %s", &toa_mjd, &toa_err, &toa_freq, toa_filename)) != 4) {
	printerror(verbose.debug, "ERROR readtimfile_custom: Error parsing line %d (only %d elements read in from '%s')\n", i+1, j, line);
	free(line);
	free(toa_filename);
	free(toa_flags);
	return 0;
      }
      toa_ssb[0] = 0;
      toa_ssb[1] = 0;
      toa_ssb[2] = 0;
      toa_freqSSB = 0;
    }

    toa_flags[0] = 0;
    for(j = 0; j < strlen(line); j++) {
      if(line[j] == '|') {
	char *str_ptr, word[1000];
	int nrwords, counter;
	verbose_definition verbose;
	cleanVerboseState(&verbose);
	for(counter = 1; counter < 1000; counter++) {
	  /* Parse the string for the nth word (counting from 1). The return
	     value is the ptr to the start of the nth word (or NULL if not
	     found). The words are separated by spaces. The total number of
	     words is returned as nrwords. The individual words cannot be larger
	     than 1000 bytes, or else there will be a memory overflow. If
	     replacetabs is set, all tabs in the input string are replaced by a
	     space. Any trailing spaces are ignored. The input string is not
	     altered. The return pointer is just a pointer to the start of the
	     requested word. The string is not null terminated after the
	     word. */
	  str_ptr = pickWordFromString(line+j+1, counter, &nrwords, 1, ' ', verbose);
	  if(str_ptr != NULL) {
	    sscanf(str_ptr, "%s", word);
	    if(strcmp(word, "-ssb1") == 0) {
	      sscanf(str_ptr, "%s %Lf", word, &(toa_ssb[0]));
	      counter++;
	    }else if(strcmp(word, "-ssb2") == 0) {
	      sscanf(str_ptr, "%s %Lf", word, &(toa_ssb[1]));
	      counter++;
	    }else if(strcmp(word, "-ssb3") == 0) {
	      sscanf(str_ptr, "%s %Lf", word, &(toa_ssb[2]));
	      counter++;
	    }else if(strcmp(word, "-freqSSB") == 0) {
	      sscanf(str_ptr, "%s %Lf", word, &(toa_freqSSB));
	      counter++;
	    }else {
	      if(strlen(toa_flags) + strlen(word) + 1 >= TEMPO3_MaxFlagsLength) {
		printerror(verbose.debug, "ERROR readtimfile_custom: The string length of the flags exceeds the maximum (%d)\nProcessing line '%s'", TEMPO3_MaxFlagsLength, line);
		free(line);
		free(toa_filename);
		free(toa_flags);
		return 0;
	      }
	      if(counter != 1)
		strcat(toa_flags, " ");
	      strcat(toa_flags, word);
	    }
	  }
	}
      }
    }
    for(j = 0; j < strlen(toa_flags); j++) {
      if(toa_flags[j] == '\n' || toa_flags[j] == '\r') {
	toa_flags[j] = 0;
      }
    }
    if(strlen(toa_flags) > 0) {
      if(toa_flags[strlen(toa_flags)-1] == ' ')
	toa_flags[strlen(toa_flags)-1] = 0;
    }
    /*    if(*(toas->mjd+i) > 54268.0 && *(toas->mjd+i) < 54268.2) {
	  printf("XXXXXXXXX FIDLING WITH TOA!!!!!!\n");
	  *(toas->mjd+i) += (25.0/360.0)*0.41/TEMPO3_SecPerDay;
	  }*/
    if(i == 0) {
      toas->data_mjdmin = toa_mjd;
      toas->data_mjdmax = toa_mjd;
    }else {
      if(toa_mjd > toas->data_mjdmax)
	toas->data_mjdmax = toa_mjd;
      if(toa_mjd < toas->data_mjdmin)
	toas->data_mjdmin = toa_mjd;     
    }

    //    printf("Adding TOA %Lf at %Lf MHz\n", toa_mjd, toa_freq);

    //int add_toa(toa_data_def *toas, long toa_index, long double toa_mjd, float toa_err, long double toa_freq, int deleted, char *toa_filename, char *toa_flags, char *toa_site, long double *toa_ssb, verbose_definition verbose)


    if(add_toa(toas, i, toa_mjd, toa_err, toa_freq, toa_freqSSB, 0, toa_filename, toa_flags, NULL, toa_ssb, verbose) == 0) {
      printerror(verbose.debug, "ERROR readtimfile_custom: Adding TOA failed\n");
      free(line);
      free(toa_filename);
      free(toa_flags);
      return 0;
    }
  }
  fclose(fin);
  if(verbose.verbose)
    printf("  MJD range: %Lf - %Lf\n", toas->data_mjdmin, toas->data_mjdmax);
  toas->data_mjdmin -= 1;
  toas->data_mjdmax += 1;
  //  residuals->mjdmin = toas->data_mjdmin;
  //  residuals->mjdmax = toas->data_mjdmax;
  free(line);
  free(toa_filename);
  free(toa_flags);
  if(verbose.verbose) {
    printf("  Reading in TOAs done\n");
  }
  return nrlines;
}


/* The site information is actually not required, or used, unless you
want to write out the timfile and/or if jumps are specified for
specific telescopes.  If force_colsite > 0, that column number is used
rather than guessing the column number.

Return 0: error
Return > 0: column number where site was found (starting at 1)
*/
int tempo3_loadSites(char *timfilename, toa_data_def *toas, unsigned long *indx_input_toas, int force_colsite, verbose_definition verbose)
{
  FILE *fin;
  char *txt, *txt2, *ret, *txtptr, *toa_site;
  int i, j, colsite;

  toas->gotBarycentreTOAs = 0; // Is potentially set to 1 when reading site codes.

  fin = fopen(timfilename, "r");
  if(fin == NULL) {
    printerror(verbose.debug, "ERROR tempo3_loadSites: Cannot open '%s'", timfilename);
    return 0;
  }
  txt = malloc(TEMPO3_MaxNrParLineLength+1);
  txt2 = malloc(TEMPO3_MaxNrParLineLength+1);
  toa_site = malloc(TEMPO3_MaxSiteStringLength+1);
  if(txt == NULL || txt2 == NULL || toa_site == NULL) {
    printerror(verbose.debug, "ERROR tempo3_loadSites: Memory allocation error.");
    fclose(fin);
    return 0;
  }


  /* Now try to guess in which column the site code is, as the general2 plugin does not allow you to output this information */
  ret = fgets(txt, TEMPO3_MaxNrParLineLength, fin);
  sscanf(txt, "%s", txt2);
  if(strcmp(txt2, "FORMAT") == 0 || strcmp(txt2, "MODE") == 0) {
    ret = fgets(txt, 1000, fin);
  }

  if(force_colsite > 0) {
    colsite = force_colsite - 1;
  }else {
    colsite = -1;
    if(ret != NULL) {
      for(i = 0; i < 100; i++) {
	if(i == 0) {
	  txtptr = strtok(txt, " ");
	}else {
	  txtptr = strtok(NULL, " ");
	}
	if(txtptr != NULL && i != 0) {  // Do not allow first column to be site, as sometimes the "filename" is the TOA number, which then gets interpreted as the site code.
	  sscanf(txtptr, "%s", txt2);
	  if(strcmp(txt2, "1") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "2") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "3") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "5") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "6") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "7") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "8") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "c") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "f") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "g") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "i") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "a") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "4") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "q") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "gbt") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "atca") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "ao") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "nanshan") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "tid43") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "pks") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "jb") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "jbdfb") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "vla") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "ncy") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "eff") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "jbm4") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "gb300") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "gb140") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "gb853") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "lap") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "hob") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "hart") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "wsrt") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "j") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "k") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "l") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "m") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "n") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "coe") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "jbmk2") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "jb42") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "p") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "@") == 0) {
	    colsite = i;
	    break;
	  }else if(strcmp(txt2, "BAT") == 0) {
	    colsite = i;
	    break;
	  }
	}
      }
    }else {
      fclose(fin);
      free(toa_site);
      free(txt);
      free(txt2);
      return 0;
    }
  }

  if(colsite != -1) {
    printf("  I guess column %d contains the site\n", colsite+1);

    rewind(fin);
    long linenr = 0;
    for(j = 0; j < toas->nrtoas; j++) {
      do {
	ret = fgets(txt, TEMPO3_MaxNrParLineLength, fin);
	linenr++;
	sscanf(txt, "%s", txt2);
      }while(strcmp(txt2, "C") == 0 || (j == 0 && strcmp(txt2, "FORMAT") == 0) || (j == 0 && strcmp(txt2, "MODE") == 0));  // Skip deleted TOA's as they were not loaded since general2 plugin didn't output them. Also skip FORMAT/MODE on first lines.
      for(i = 0; i <= colsite; i++) {
	if(i == 0) {
	  txtptr = strtok(txt, " ");
	}else {
	  txtptr = strtok(NULL, " ");
	}
      }
      if(txtptr == NULL) {
	printerror(verbose.debug, "ERROR tempo3_loadSites: Cannot read site code from column %d while reading in line %d.", colsite+1, linenr);
	fclose(fin);
	free(toa_site);
	free(txt);
	free(txt2);
	return 0;
      }
      //      fprintf(stderr, "XXXXX reading word from '%s' from line '%s'\n", txtptr, txt);
      long actual_indx = -1;  // Initialise to something to avoid compiler error. Note that this should give a segfault if the following loop doesn't exit via the break, but I think this should always happen unless something is very wrong.
      for(i = 0; i < toas->nrtoas; i++) {
	if(indx_input_toas[i] == j) {  // indx_input_toas[i]=Line number in TOA that ended up stored as index i. j=current line number that is read in.
	  actual_indx = i;  // So this is the location in the array of the current line number
	  break;
	}
      }
      if(strcmp(txt2, toas->filenames[actual_indx]) != 0) { 
	printwarning(verbose.debug, "WARNING tempo3_loadSites: When reading in the site codes, the TOA identifiers didn't match what was previously read in. The site codes might be incorrectly matched to tht TOAs.");
	printwarning(verbose.debug, "WARNING tempo3_loadSites: Expected ID '%s', but got '%s'", toas->filenames[actual_indx], txt2);
      }
      //      strncpy(toas->site[indx_input_toas[j]], txtptr, TEMPO3_MaxSiteStringLength-1);
      strncpy(toa_site, txtptr, TEMPO3_MaxSiteStringLength-1);
      for(i = 0; i < strlen(toa_site); i++) {
	if(toa_site[i] == '\n' || toa_site[i] == '\r')
	  toa_site[i] = 0;
      }
      if(toas->site[actual_indx] != NULL) {
	free(toas->site[actual_indx]);
      }
      toas->site[actual_indx] = malloc(strlen(toa_site)+1);
      if(toas->site[actual_indx] == NULL) {
	printerror(verbose.debug, "ERROR tempo3_loadSites: Memory allocation error\n");
	fclose(fin);
	free(toa_site);
	free(txt);
	free(txt2);
	return 0;
      }
      strcpy(toas->site[actual_indx], toa_site);
      //      printf("TOA %Lf (%s): index=%ld site=%s\n", toas->mjd[actual_indx], toas->filenames[actual_indx], indx_input_toas[j], toas->site[actual_indx]);
      //      if(strcmp(toas->site[indx_input_toas[j]], "@") == 0 || strcmp(toas->site[indx_input_toas[j]], "BAT") == 0) {
      if(strcmp(toas->site[actual_indx], "@") == 0 || strcmp(toas->site[actual_indx], "BAT") == 0) {
	toas->gotBarycentreTOAs = 1;
      }
    }
  }else {
    printwarning(verbose.debug, "WARNING tempo3_loadSites: determining in which column the site name is stored failed. Your site code might not be hardcoded in the code. The site information will not be available after barycentring. You can use the -sitecolumn option as well to avoid this problem.");
  }



  fclose(fin);
  free(toa_site);
  free(txt);
  free(txt2);
  return colsite+1;
}

// Load timfile with name timfilename. If loadoriginalsites != 0, the
// site codes are set to whatever the site codes were in timfile,
// despite the TOAs are barycentred. This information is required, for
// example, when dealing with jumps. When loading in this site
// information, force_colsite indicates which column the site code is
// in. If zero, an attempt is made to determine this
// automatically. The TOA's are described by a parfile with name
// parfilename (used to do the barycentring). The ephemeris should
// also be in memory (eph). A system call to tempo2 will be
// made. After running tempo2 with the general2 plugin the relevant
// information from tempo2 about each toa is extracted. indx will be
// initialised with the a list. The first element gives the toa number
// (in input file) corresponding to the smallest mjd etc. This is used
// when linking the site code as read in separately to the correct
// toa. If nosort != 0, no sorting is done and indx is simply an list
// such that indx[i] = i. If baryssbdumpfilename != NULL, this is the
// filename used to dump the barycentre information obtained from
// tempo2 to a file which is kept for the user.
//
// You probably should have a line like:   
// residuals->mjdmin = toas->data_mjdmin;
// residuals->mjdmax = toas->data_mjdmax;
//
// Return 0: error
// Return 1: ok
int loadtimfile(char *timfilename, int loadoriginalsites, int force_colsite, char *parfilename, ephemeris_def eph, toa_data_def *toas, int nosort, unsigned long **indx, int repeatable, char *baryssbdumpfilename, verbose_definition verbose)
{
  char *username, *cmd, *tmpfilename_tempo2_output, *tmpfilename_filtered_par, *tmpfilename_filtered_tim;
  if(getUsername(&username, verbose) == 0) {
    username = malloc(8);
    if(username == NULL) {
      printerror(verbose.debug, "ERROR loadtimfile: Memory allocation error");
      return 0;
    }
    sprintf(username, "Unknown");
  }

  tmpfilename_tempo2_output = malloc(MaxFilenameLength+1);
  tmpfilename_filtered_par = malloc(MaxFilenameLength+1);
  tmpfilename_filtered_tim = malloc(MaxFilenameLength+1);
  cmd = malloc(MaxFilenameLength + 1);
  if(tmpfilename_tempo2_output == NULL || tmpfilename_filtered_par == NULL || cmd == NULL || tmpfilename_filtered_tim == NULL) {
    printerror(verbose.debug, "ERROR loadtimfile: Memory allocation error");
    free(username);
    return 0;
  }
  sprintf(tmpfilename_filtered_par, "%sjunk_%s_%ld.par", TEMPO3_TMP_DIR, username, randomUnsignedInt());
  sprintf(tmpfilename_tempo2_output, "%sjunk_%s_%ld.tempo2", TEMPO3_TMP_DIR, username, randomUnsignedInt());
  if(baryssbdumpfilename != NULL) {
    strcpy(tmpfilename_filtered_tim, baryssbdumpfilename);
  }else {
    sprintf(tmpfilename_filtered_tim,  "%sjunk_%s_%ld.toa", TEMPO3_TMP_DIR, username, randomUnsignedInt());
  }
  free(username);

  // Make a new par file with some potential fatal parameters filtered out, which are not necessary to generate the BATs
  if(verbose.debug == 0)
    verbose.verbose = 0;
  if(filter_parfile_for_tempo2(parfilename, tmpfilename_filtered_par, &(eph.par), verbose) == 0) {
    printerror(verbose.debug, "ERROR loadtimfile: Making a tempo2 suitable filtered parameter file failed");
    free(cmd);
    free(tmpfilename_filtered_par);
    free(tmpfilename_tempo2_output);
    free(tmpfilename_filtered_tim);
    return 0;
  }

  sprintf(cmd, "tempo2 -nofit -output general2 -s \"grepthislineout\\t{bat}\\t{err}\\t{earth_ssb1}\\t{earth_ssb2}\\t{earth_ssb3}\\t{freq}\\t{freqSSB}\\t{FILE}\\t|{FLAGS}\\n\" -f %s %s >& %s\n", tmpfilename_filtered_par, timfilename, tmpfilename_tempo2_output);
  if(verbose.debug) {
    fprintf(stderr, "  tempo2 is called as: '%s'\n", cmd); 
  }
  system(cmd);
  //  int ret = system(cmd);
  //  printf("XXXXX: %d\n", ret);

  //  sprintf(cmd, "ls -la %s\n", tmpfilename_tempo2_output);
  //  ret = system(cmd);
  //  printf("XXXXX: %d\n", ret);

  //  sprintf(cmd, "cat %s\n", tmpfilename_tempo2_output);
  //  ret = system(cmd);
  //  printf("XXXXX: %d\n", ret);

  free(cmd);
  remove(tmpfilename_filtered_par);
  free(tmpfilename_filtered_par);
  //    sprintf(cmd, "grep 'grepthislineout' %s | sed '{s/grepthislineout/ /g}' > %s\n", tmpfilename_tempo2_output, tmpfilename_filtered_tim);
  //    system(cmd);
  if(verbose.debug) {
    fprintf(stderr, "  Calling filter_tempo2_general2_output\n"); 
  }
  if(filter_tempo2_general2_output(tmpfilename_tempo2_output, tmpfilename_filtered_tim, nosort, indx, repeatable, verbose) == 0) {
    printerror(verbose.debug, "ERROR loadtimfile: Filtering of tempo2 TOA output failed\n");
    //      sprintf(cmd, "rm -f %s %s %s %s\n", tmpfilename, tmpfilename_tempo2_output, tmpfilename_filtered_tim, tmpfilename_filtered_par);
    //      system(cmd);
    remove(tmpfilename_tempo2_output);
    free(tmpfilename_tempo2_output);
    return 0;
  }
  remove(tmpfilename_tempo2_output);
  free(tmpfilename_tempo2_output);
  //    if(verbose.debug) {
  //      fprintf(stderr, "  Calling readtimfile_custom\n"); 
  //    }
  // 0 was timfile_shortformat, but after system call the ssb information should always be there
  toas->nrtoas = readtimfile_custom(tmpfilename_filtered_tim, toas, 0, verbose);
    //    if(verbose_state.debug) {
    //      fprintf(stderr, "  Calling readtimfile_custom done\n"); 
    //    }
    //    sprintf(cmd, "rm -f %s %s %s %s\n", tmpfilename, tmpfilename_tempo2_output, tmpfilename_filtered_tim, tmpfilename_filtered_par);
    //    system(cmd);

  
  if(baryssbdumpfilename == NULL) {   // Keep this information if requested by the user, otherwise remove the temporary file
    remove(tmpfilename_filtered_tim);
  }
  free(tmpfilename_filtered_tim);

  if(loadoriginalsites) {
    if(tempo3_loadSites(timfilename, toas, *indx, force_colsite, verbose) == 0) {
      printerror(verbose.debug, "ERROR loadtimfile: Loading site information failed");
      return 0;
    }
  }
  return 1;
}
