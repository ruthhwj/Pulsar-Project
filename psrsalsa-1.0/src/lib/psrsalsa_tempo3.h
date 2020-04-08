#include "psrsalsa.h"

#define TEMPO3_TMP_DIR   "/tmp/"      // Don't forget the trailing /

/* Comment following line out to disable GTK stuff (setting to 0 doesn't help). This disables gtk during compilation.  */
//#define TEMPO3_EnableGTK    1

#define TEMPO3_MaxNrFderivatives        12           /* Maximum number of frequency derivatives that are supported beyond F0. Setting this lower than 2 is probably a bad idea and beyond 9 is not supported at the moment (calculations and toggling). */
#define TEMPO3_MaxNrDMderivatives        9           // I.e. 1 means DM and DM1 are defined, but DM2 isn't
#define TEMPO3_MaxNrGlitchFderivatives   5           /* Maximum number of frequency derivatives that are supported beyond GLF0 */
#define TEMPO3_MaxNrValues              10           /* Maximum number of VALUE%d parameters that can be specified in par file */
#define TEMPO3_MaxNrParLines           100
#define TEMPO3_MaxNrParLineLength     1000           // Actually also for tim file
#define TEMPO3_MaxFlagsLength         1000           // Total maximum string length of the flags
#define TEMPO3_MaxSiteStringLength     100           // Maximum length of the site code



#define TEMPO3_SecPerDay                 86400.0         /* 24.0*3600.0 */
#define TEMPO3_SPEEDOFLIGHT              299792458.0    // m/s
#define TEMPO3_AU                        149597870700.0 // m

/* Some constants defined in the library */
extern int TEMPO3_MaxNrWaves;          /* Maximum number of WAVE parameters supported. Can be changed on command line, but shouldn't be changed after calling initialise_tempo3_lib() */
extern int TEMPO3_MaxNrGlitches;       /* Maximum number of glitches which are supported. Can be changed on command line, but shouldn't be changed after calling initialise_tempo3_lib() */
extern int TEMPO3_MaxNrJumps;          /* Maximum number of jumps which are supported. Can be changed on command line, but shouldn't be changed after calling initialise_tempo3_lib()  */
extern int TEMPO3_MaxNrCompanions;     // Max nr of companions that can be defined, but shouldn't be changed after calling initialise_tempo3_lib()

extern int TEMPO3_PARAM_PSRJNAME;  // The JNAME, not a fit parameter, so need to be dealt with differently in places
extern int TEMPO3_PARAM_POSEPOCH;  // Epoch (MJD) at which the RAJ and DECJ are defined
extern int TEMPO3_PARAM_RAJ     ;  // The RAJ in radians
extern int TEMPO3_PARAM_DECJ    ;  // The DECJ in radians
extern int TEMPO3_PARAM_PMRA    ;  // The proper motion in RA in milli arcseconds per year (difference with original value)
extern int TEMPO3_PARAM_PMDEC   ;  // The proper motion in DEC in milli arcseconds per year (difference with original value)
extern int TEMPO3_PARAM_PHASE0  ;  // A phase offset (in units phase) applied to all residuals
extern int TEMPO3_PARAM_DMEPOCH ;  // Epoch (MJD) at which DM and time derivatives are defined
extern int TEMPO3_PARAM_DM      ;  // The actual fit parameter is the DM difference from the actual DM, so there are some non-standard exceptions to be handled. The DM is followed by TEMPO3_MaxNrDMderivatives time derivatives
extern int TEMPO3_PARAM_PEPOCH  ;
extern int TEMPO3_PARAM_F0      ;  // Points to F0 (rotational frequency in Hz), but it is followed by TEMPO3_MaxNrFderivatives derivatives
extern int TEMPO3_PARAM_NBRAKE  ;  // Braking index fit parameter
extern int TEMPO3_PARAM_WAVEEPOCH;
extern int TEMPO3_PARAM_WAVEOM  ;  // Defines the period of the fitwaves, not a fit parameter
extern int TEMPO3_PARAM_WAVESIN ;  // Points to the first wavesin term, followed by the TEMPO3_MaxNrWaves-1 others
extern int TEMPO3_PARAM_WAVECOS ;  // Points to the first wavecos term, followed by the TEMPO3_MaxNrWaves-1 others
extern int TEMPO3_PARAM_GLEP    ;  // Glitch epochs for TEMPO3_MaxNrGlitches glitches
extern int TEMPO3_PARAM_GLPH    ;  // Phase step for TEMPO3_MaxNrGlitches glitches
extern int TEMPO3_PARAM_GLF0    ;  // F0 for TEMPO3_MaxNrGlitches glitches, followed by TEMPO3_MaxNrGlitchFderivatives times (TEMPO3_MaxNrGlitches-1) other parameters
extern int TEMPO3_PARAM_GLTD    ;
extern int TEMPO3_PARAM_GLF0D   ;
//extern int TEMPO3_PARAM_GLF1D   ;  // This is not an tempo2 parameter, and it is degenerate with glf0d as is implemented?
extern int TEMPO3_PARAM_GLNBRAKE;  // Not yet implemented. Change of braking index after glitch
extern int TEMPO3_PARAM_GLEPMIN ;  // Min bound (mjd) on glitch epochs for TEMPO3_MaxNrGlitches glitches
extern int TEMPO3_PARAM_GLEPMAX ;  // Max bound (mjd) on glitch epochs for TEMPO3_MaxNrGlitches glitches
extern int TEMPO3_PARAM_JUMPS   ;  // The jump for different flagged toa's
extern int TEMPO3_PARAM_VALUE   ;  // The VALUE parameters set
extern int TEMPO3_PARAM_BINARY  ;  // The type of binary model to use, not a fit parameter and stored in a different place
extern int TEMPO3_PARAM_T0      ;  // TEMPO3_MaxNrCompanions T0 parameters
extern int TEMPO3_PARAM_TASC    ;  // TEMPO3_MaxNrCompanions parameters
extern int TEMPO3_PARAM_PB      ;  // TEMPO3_MaxNrCompanions parameters
extern int TEMPO3_PARAM_A1      ;  // TEMPO3_MaxNrCompanions parameters
extern int TEMPO3_PARAM_OM      ;  // TEMPO3_MaxNrCompanions parameters
extern int TEMPO3_PARAM_ECC     ;  // TEMPO3_MaxNrCompanions parameters
extern int TEMPO3_PARAM_EPS1    ;  // TEMPO3_MaxNrCompanions parameters
extern int TEMPO3_PARAM_EPS2    ;  // TEMPO3_MaxNrCompanions parameters
extern int TEMPO3_PARAM_OMDOT   ;  // TEMPO3_MaxNrCompanions parameters
extern int TEMPO3_PARAM_PBDOT   ;  // TEMPO3_MaxNrCompanions parameters
extern int TEMPO3_PARAM_A1DOT   ;  // TEMPO3_MaxNrCompanions parameters
extern int TEMPO3_PARAM_GAMMA   ;  // TEMPO3_MaxNrCompanions parameters
extern int TEMPO3_PARAM_TZRMJD  ;
extern int TEMPO3_PARAM_TZRFRQ  ;
extern int TEMPO3_PARAM_TZRSITE ;  // not a fit parameter, so need to be dealt with differently in places
extern int TEMPO3_PARAM_MODE    ;  // Determines if a weighted fit is used, not a fit parameter and stored in a different place
extern int TEMPO3_PARAM_TRACK   ;  // Determines if pulse numbering is used, not a fit parameter and stored in a different place
extern int TEMPO3_PARAM_TRES    ;  // The RMS of the residual in micro-seconds
extern int TEMPO3_PARAM_TOTAL   ;  // Total nr of parameters to be determined automatically

// This defines a structure which describes all identifiers used in a
// .par file, its meaning etc. This information is provided by the
// function initialise_tempo3_lib()

typedef struct {
  char ***identifiers;    // Pointer to TEMPO3_PARAM_TOTAL pointers to 2 pointers, giving the parfile names of the parameters (2nd NULL if not used)
  char **print_format;    // Pointer to TEMPO3_PARAM_TOTAL pointers, giving the the print format, i.e. "%30.20Le"
  char **units;           // Pointer to TEMPO3_PARAM_TOTAL pointers, giving the units of the parameter
  char **description;     // Pointer to TEMPO3_PARAM_TOTAL pointers, giving a description of the parameter
}parfile_description_def;


/* 
These are parameters which are part of an ephemeris, but there are
others as well. See ephemeris_def. When add something, search for
gltd.
*/
typedef struct {
  long double *parameter;      // The floating point parameters, identified via TEMPO3_PARAM_XXXXX
  long double *dplus, *dmin;   // The associated errors, or NULL if not stored.
  char psrjname[100];                                // Jname of the pulsar
  long double dm[1+TEMPO3_MaxNrDMderivatives];       // The DM and derivatives as specified in the par file. Their differences are fit parameters.
  long double raj, decj;                             // The RAJ and DECJ  as specified in the par file. The difference is one of the fit parameters.
  long double pmra, pmdec;                           // The proper motion as specified in the par file. The difference is one of the fit parameters.
  char **jumpflag; //  char jumpflag[TEMPO3_MaxNrJumps][100];           /* The flag corresonding to the jump */
  int *jumpmeaning; // 0 = jumpflag is the "FLAG ID" combination, 1 = jumpflag is the site code
  int binarymodel;                          // 1 = BT model, 2 = T2, 3 = HYPER
  int *companion_defined;                   // Will be allocated as TEMPO3_MaxNrCompanions ints
  int nrglitches, nrwaves, nrjumps;
  int mode;                                 /* 0 = no errorbars, 1 = weighted fit */
  int track;                                /* 0 = off (default), -2 is pulse numbering*/
  char tzrsite[100];
}tempo3_parameters_def;
// Don't forget to change:
//  clean_parfile(tempo3_parameters_def *parfile)
//  allocate_parfile(tempo3_parameters_def *parfile)
//  free_parfile(tempo3_parameters_def *parfile)
//  copy_parfile(tempo3_parameters_def *copy, tempo3_parameters_def *orig)

/* 
These are parameters which are part of an ephemeris, but there are
others as well. See ephemeris_def.
*/
typedef struct {
  int nrwraps;
  long double mjd[1000];
  long double phase[1000];
}tempo3_wraps_def;

typedef struct {
  tempo3_parameters_def par;
  tempo3_wraps_def *wraps;     // Set to NULL if not defined
  char **unused_parfile_lines;   // Has dimensions [TEMPO3_MaxNrParLines][TEMPO3_MaxNrParLineLength] after allocation;
  int nrunused_parfile_lines;
  int *fixed, *paramset, *paramlinked;
}ephemeris_def;

typedef struct {
  int nrtoas;                               /* The number of TOAs that are active */
  int nrtoasdeleted;                        /* The number of deleted TOAs. Total=nrtoas+nrtoasdeleted */
  long double *mjd;                         /* The TOAs as read in from TIM file */
  long double *freqSite, *freqSSB;          // The observing frequency at the observatory and at the SSB in MHz
  long double *ssb1, *ssb2, *ssb3;          /* Position earth w.r.t. solar system barycentre (SSB) in lt-s (if != NULL) */
  float *err;                               /* The TOAs errorbar in microsec */
  float *mjdx, *freqx;                      /* The TOAs and frequencies converted in floats */
  char **filenames;                         /* The filenames corresponding to the TOAs */
  char **flags;                             /* The flags corresponding to the TOAs */
  char **site;                              /* The sites corresponding to the TOAs */
  int gotBarycentreTOAs;                    /* Set to 1 if one or more toa's have their site set to the SSB, otherwise 0 */ 
  int *deleted;                             /* Set to 1 if deleted */
  long double data_mjdmin, data_mjdmax;     /* The time-span covered */
  long double *nu, *nudot, *nu_mjd;         /* Measured values for nu and nudot, if calculated */
}toa_data_def;


// FUCTION DESCRIPTIONS

// Should be called as soon as TEMPO3_MaxNrWaves,
// TEMPO3_MaxNrGlitches, TEMPO3_MaxNrJumps, TEMPO3_MaxNrCompanions are
// defined (or immediately if not changing these values). This
// initialises some global variables used in the library, such as
// descriptions of the different ephemeris parameters and their units
// which are used when parameters are reported. The associated memory
// will be released after calling cleanup_tempo3_lib(), which sassociatedhould
// be done when the program is finished using the tempo3 library.
// 
// Return 0: Error
// Return 1: Success
int initialise_tempo3_lib(verbose_definition verbose);

// The memory used after calling initialise_tempo3_lib() will be
// released, which should be done when the program is finished using
// the tempo3 library.
//
// Return 0: Error
// Return 1: Ok
int cleanup_tempo3_lib(verbose_definition verbose);

// Allocates memory in the ephemeris to store information, and
// initialise their values to be unset. If allowerrors is nonzero,
// space is allocated to store errors on parameters.
//
// Return 0: error
// Return 1: ok
int initialise_ephemeris(ephemeris_def *eph, int allowerrors, verbose_definition verbose);

// Releases the memory allocated by initialise_ephemeris().
//
// Return 0: error
// Return 1: ok
int free_ephemeris(ephemeris_def *eph, verbose_definition verbose);

// Reads in an ephemeris from a file. The ephemeris needs to be
// initialised first with the function initialise_ephemeris(). If hex
// is nonzero the values in the parameter file are assumed to be a hex
// representation of the bytes of the floating point value, rather
// than decimal numbers.
//
// Return 0 on error, 1 on success
int read_ephemeris(char *filename, ephemeris_def *eph, int hex, int nowarning, verbose_definition verbose);

// Evaluate the spin-frequency at the given barycentric mjd.
//
// Return 0: error
// Return 1: ok
int evaluate_ephemeris_spinfrequency(ephemeris_def eph, long double bary_mjd, long double *spinfreq, verbose_definition verbose);

// The ephemeris might contain VALUE parameters. The call to this
// function copies their values to the appropriate parameters such
// that the ephemeris can be evaluated.
void ephemeris_setLinkedValues(ephemeris_def *eph);

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
int print_ephemeris(FILE *stream, ephemeris_def *eph, ephemeris_def *eph2, toa_data_def *toas, int ignorewraps, int ignoreunrecognizedlines, int ignorefitted, int ignoreerrors, int autoquiet, int hex, verbose_definition verbose);

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
int initialise_toas(toa_data_def *toas, long nrtoas, long nrtoasdeleted, int ssb, int floatversion, int freqevol, verbose_definition verbose);

// Release the memory previously allocated with initialise_toas()
void free_toas(toa_data_def *tim);

// Add a TOA to toas at position toa_index. Its properties are given
// by toa_mjd, toa_err (in microseconds), toa_freq , if deleted != 0
// then the TOA is flagged to be deleted, toa_filename, toa_flags,
// toa_site, toa_ssb (X,Y,Z = position earth w.r.t. solar system
// barycentre in lt-s). nrtoas and nrtoasdeleted are not updated.
//
// Return 0: error
// Return 1: ok
int add_toa(toa_data_def *toas, long toa_index, long double toa_mjd, float toa_err, long double toa_freq, long double toa_freqSSB, int deleted, char *toa_filename, char *toa_flags, char *toa_site, long double *toa_ssb, verbose_definition verbose);

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
int loadtimfile(char *timfilename, int loadoriginalsites, int force_colsite, char *parfilename, ephemeris_def eph, toa_data_def *toas, int nosort, unsigned long **indx, int repeatable, char *baryssbdumpfilename, verbose_definition verbose);

// Write out the information from toas to a file with name
// filename. If nturn != NULL, the pulse numbering information is
// included, corrected for wraps if != NULL. If isbarycentred != 0,
// then all the TOA' are assumed to be defined at the barycentre and
// the site codes are written out accordingly independent of what the
// original site was set to.
//
// Return 0: error
// Return 1: ok
int writetimfile(char *filename, toa_data_def toas, long long *nturn, tempo3_wraps_def *wraps, int isbarycentred, verbose_definition verbose);
