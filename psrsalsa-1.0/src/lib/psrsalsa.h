//START REGION DEVELOP
//#include <math.h>
//#include "pumadata.h"
//START REGION RELEASE
#include <fitsio.h>
#include <psrsalsa_defines.h>
#include <psrsalsa_typedefs.h>
//START REGION DEVELOP
#ifdef HAVEPGPLOT2 
  #include "pgplot2.h"
#endif
//START REGION RELEASE

/* ***********************************
   Functions in psrio.c 
   *********************************** */

/*
  Return 1 if data format is recognized, or zero otherwise.
*/
int isValidPSRDATA_format(int format);

/* Print the available data formats to device (could be for instance
   stdio). nrspaces defines the number of spaces before each line.*/
void printPSRDataFormats(FILE *printdevice, int nrspaces);

/* Parse a string to find the data format for the -oformat/-iformat
   command line option. Returns 0 on error, otherwise the data
   format. */
int parsePSRDataFormats(char *cmd);

/* Clean the struct which is probably filled with random junk at the
   time of declaration. */
void cleanPSRData(datafile_definition *datafile, verbose_definition verbose);

/* Copy the struct to another, except the file pointers. No memory is
   allocated to hold data and the data pointer is not set. The
   destination should have been initialised with cleanPSRData() at
   some point before calling this function.

If return value = 0, then there is an error. Returns 1 on success.
*/
int copy_params_PSRData(datafile_definition datafile_source, datafile_definition *datafile_dest, verbose_definition verbose);

/* Copy the string filename in the data structure. If no memory is
allocated to hold the filename it is allocated. If already momory is
allocated, it will be freed first. If filename is the NULL pointer,
the pointer in the data struct will point to an empty string.

If return value = 0, then there is a memory allocation
error. Returns 1 on success.  */
int set_filename_PSRData(datafile_definition *datafile_dest, char *filename, verbose_definition verbose);
int set_psrname_PSRData(datafile_definition *datafile_dest, char *psrname, verbose_definition verbose);
int set_observatory_PSRData(datafile_definition *datafile_dest, char *observatory, verbose_definition verbose);
int set_institute_PSRData(datafile_definition *datafile_dest, char *institute, verbose_definition verbose);
int set_instrument_PSRData(datafile_definition *datafile_dest, char *instrument, verbose_definition verbose);
int set_scanID_PSRData(datafile_definition *datafile_dest, char *scanID, verbose_definition verbose);
int set_observer_PSRData(datafile_definition *datafile_dest, char *observer, verbose_definition verbose);
int set_projectID_PSRData(datafile_definition *datafile_dest, char *projectID, verbose_definition verbose);

//START REGION DEVELOP
/* 
Takes the command line and stores it in the scanid field of the datafile struct.
This function is now obsolete.
*/
//void putCmdlineInScanid(datafile_definition *datafile, int argc, char **argv, verbose_definition verbose);
//START REGION RELEASE

/* Tries to guess what format the file is. Verbose level determines 
   nr of spaces before output. If debug is set, more information is
   printed. If noerror is set, no error is generated if data format
   cannot be determined. Returns format:
     -3 is an other error
     -2 is a known type which cannot be read
     -1 is file opening/read error
      0 is unknown type

//START REGION DEVELOP
   Formats that should be recognized are: 
         PUMA_format, 
         PSRCHIVE_ASCII_format, 
         EPN_format,
         FITS_format, 
	 PSRSALSA_BINARY_format
         AO_ASCII_1_format
         AO_ASCII_2_format
	 PPOL_format
	 PPOL_SHORT_format
	 SIGPROC_format
//START REGION RELEASE
*/
int guessPSRData_format(char *filename, int noerror, verbose_definition verbose);

/* Open filename which is in the specified format. Set the
   enable_write flag to open the file with write permissions. If
   enable_write is not set, cleanPSRData() is called, meaning that if
   this is done in the main progamme, there will be a memory leak
   since the struct will be re-initialised.  read_in_memory is set the
   file is read into memory and the file is closed. You can acces the
   data as if the file is still open. If the format is set to zero it
   will try to guess the format. If enable_write is not set, the
   datafile structure is always cleared before any reading/writing is
   done. If nowarnings == 0, warnings are shown. If 1, no warnings are
   generated for incomplete headers (the checks done AFTER the
   file-format specific reading has been done), which can be useful
   for applications which do not rely on the header information. If 2,
   all warnings are ignored while reading in the header. If nocounters
   is set it will not write out any counters (useful for log
   files). Verbose level determines nr of spaces before output. If
   debug is set, more information is printed. Returns 0 on error. */
int openPSRData(datafile_definition *datafile, char *filename, int format, int enable_write, int read_in_memory, int nowarnings, verbose_definition verbose);

/* Close the files and free the memory if required. Data is written to
   file if writing is buffered (useful for ascii formats). To be
   specific:

   - If requested by the dumpOnClose flag, the data is written to the file.
   - If the file is opened it is closed.
   - If the format is in MEMORY_format the allocated memory (data) is released
   - If perserve_header_info is 1, memory associated with things like the pulsar name will not be freed

Returns 0 on success. 
*/
int closePSRData(datafile_definition *datafile, int perserve_header_info, verbose_definition verbose);

/* Automatically done in readHeaderPSRData if verbose is set. If
   update is set, the first line indicates these are updated header
   commands (used after using -header command line). Indent indicates
   the nr of spaces before output. If debug is set, the output is
   more precise, but less concise. */
void printHeaderPSRData(datafile_definition datafile, int update, verbose_definition verbose);

/* Read the header. If readnoscales is set, the scales/offsets/weights
  of the subints are not read in (if supported/used by file
  format). Useful when only interested in header information.  Verbose
  level determines nr of spaces before output. If debug is set, more
  information is printed. If nowarnings == 0, warnings are shown. If
  1, no warnings are generated for incomplete headers (the checks done
  AFTER the file-format specific reading has been done), which can be
  useful for applications which do not rely on the header
  information. If 2, all warnings are ignored while reading in the
  header. Returns 0 on error, 1 on success. */
int readHeaderPSRData(datafile_definition *datafile, int readnoscales, int nowarnings, verbose_definition verbose);

/* Write the header, including the present command line (if argc >
   0). If cmdOnly is set, only the command is writen out. This makes
   the history table without timestamps and therefore results are
   exactly reproducable.  Verbose level determines nr of spaces before
   output. 
   Returns 1 if successful or 0 on error. */
int writeHeaderPSRData(datafile_definition *datafile, int argc, char **argv, int cmdOnly, verbose_definition verbose);

/* Obtains a pointer pulse_ptr which points to a single pulse (or
   subint) starting at binnr etc. The behaviour is unpredictable if
   you read in beyond a subint/polarization/frequency channel. An
   error will be generated if the data was not yet read into memory.

  Returns 1 if successful, 0 on error. */
int get_pointer_PulsePSRData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, float **pulse_ptr, verbose_definition verbose);

/* Read a single pulse (or subint) starting at binnr and with length
   nrSamples. The behaviour of this function is unpredictable if you
   read in beyond a subint/polarization/frequency channel. The
   reading should not nessesarily be from start to end. If the data is
   psrfits, and psrfits_set_use_weighted_freq function is used to set
   the weighted frequency flag, the centre frequency is set to the
   wieghted frequency of the first pulse read in. Returns 1 if
   successful. */
int readPulsePSRData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose);

/* Write a single pulse/subint starting at binnr and with length
   nrSamples. The behaviour of this function is unpredictable if you
   read in beyond a subint/polarization/frequency channel. The
   writing should not nessesarily be from start to end. Return 1 if
   successful. */
int writePulsePSRData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Read the whole observation at once in memory. The data loops over
   pulse number, frequency channel, polarization channel and bin nr.
   If the data is psrfits, and psrfits_set_use_weighted_freq function
   is used to set the weighted frequency flag, the centre frequency is
   set to the wieghted frequency of the first pulse read in. If
   nocounters is set, no counters are shown. Return 0 if there is an
   error. */
int readPSRData(datafile_definition *datafile, float *data, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Write the whole dataset at once from memory. The data loops over
   pulse number, frequency channel, polarization channel and bin
   nr. Set pa data pointers to something random if not used. Returns 1
   if successful. If nocounters is set, no counters are shown. */
int writePSRData(datafile_definition *datafile, float *data, verbose_definition verbose);

/* Read file and produce a pulse profile. Set zapMask to NULL is you
   don't want to zap pulses. Set polchan to select the polarization
   channel. Returns 1 if successful. */
int read_profilePSRData(datafile_definition datafile, float *profileI, int *zapMask, int polchan, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Identical to read_profilePSRData, but now you can set nread and nskip. */
int read_partprofilePSRData(datafile_definition datafile, float *profileI, int *zapMask, int polchan, long nskip, long nread, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Read data file and produce a list of rms'es and averages. Set
   zapMask to NULL is you don't want to zap pulses. The defined
   regions are NOT included, unless invert is set. If onlyI is set
   only the total intensity is generated. Set polchan to select the
   polarization channel. Set freqchan to -1 if you want to get average
   of all channels instead of a specific channel. Set rms to NULL if
   you're not interested in the rms'es, likewise for avrg. Return 0 if
   there is an error. This function is very similar to rmsPSRData(). */
int read_rmsPSRData(datafile_definition datafile, float *rms, float *avrg, int *zapMask, pulselongitude_regions_definition *regions, int invert, int polchan, int freqchan, verbose_definition verbose);

//START REGION DEVELOP

/* Read data file and obtain an rms and average for a given subint,
   polchan and freqchan. The defined regions are NOT included, unless
   invert is set. Set rms to NULL if you're not interested in the rms,
   likewise for avrg. Return 0 if there is an error. This function is
   very similar to read_rmsPSRData().*/
int rmsPSRData(datafile_definition datafile, float *rms, float *avrg, pulselongitude_regions_definition *regions, int invert, long subint, long polchan, long freqchan, verbose_definition verbose);

//START REGION DEVELOP

/* Read data file and produce a list of correlations and averages. The
   file pointer should be set to the start of the data. Set zapMask to
   NULL is you don't want to zap pulses. The defined regions are NOT
   included, unless invert is set. Set polchan to select the
   polarization channel.  */
int read_correlPSRData(datafile_definition datafile, float *correl, float *avrg, int *zapMask, pulselongitude_regions_definition *regions, int invert, int polchan, verbose_definition verbose);

/* Read data file and produce a list of peak rms'es and averages
   (after adding addbins together). Set zapMask to NULL is you don't
   want to zap pulses. The defined regions are NOT included, unless
   invert is set. Set polchan to select the polarization channel. The
   defined regions are NOT included, unless invert is set.

   If return value = 0, then there is a memory allocation error.
*/
int read_peakrmsPSRData(datafile_definition datafile, float *rms, float *avrg, int *zapMask, pulselongitude_regions_definition *regions, int invert, int addbins, int polchan, verbose_definition verbose);

//START REGION RELEASE

/* If the data is in FREQMODE_FREQTABLE mode, check if the data really
   makes use of this feature. If all the channels are equally spaced
   in line with the given bandwidth and centre-frequency, then change
   the data to FREQMODE_UNIFORM, which is compatible with more
   functionality. 

   Return -
   0: Error
   1: Not converted, because input data not in FREQMODE_FREQTABLE mode
   2: Not converted, because frequency channels are not uniformly separated
   3: Converted
*/
int convert_if_uniform_frequency_spacing(datafile_definition *datafile, int nowarnings, verbose_definition verbose);

/*
  If the data is in FREQMODE_FREQTABLE mode, force all channels to
   become equally spaced with a single bandwidth and
   centre-frequency. Then change the data to FREQMODE_UNIFORM, which
   is compatible with more functionality. NO CHECK IS DONE IF THIS IS
   REASONALBLE!

   Return -
   1: Not converted, because input data not in FREQMODE_FREQTABLE mode
   3: Converted
 */
int force_uniform_frequency_spacing(datafile_definition *datafile, verbose_definition verbose);

// Set the struct to the default state (no verbose)
void cleanVerboseState(verbose_definition *verbose_state);

// Copy verbose state
void copyVerboseState(verbose_definition verbose_state_src, verbose_definition *verbose_state_dst);

/* Parse the command line for -header or -headerUFL options and applies the changes
   to psrdata. Returns 0 if parse error */
int PSRDataHeader_parse_commandline(datafile_definition *psrdata, int argc, char **argv, verbose_definition verbose);

/* Show command line options for the -header option */
void printHeaderCommandlineOptions(FILE *printdevice);

// Shows the options of the -header gentype option
void printHeaderGentypeOptions(FILE *printdevice);

//START REGION DEVELOP

/* Allocates memory for paswing, pulse longitudes and bins and RMS
   values. Returns 0 if not successful */
//int alloc_paswing_PSRdata(datafile_definition *psrdata, int paswing, int longitudes, int bins, int rmsValues, verbose_definition verbose);

//START REGION RELEASE

/* This function makes a clone of the data and puts it in a (as yet
   non-existent) clone. So the clone will contain have all data
   available in memory. Returns 1 on success, 0 on error. The data
   must have been read into memory first. */
int make_clone(datafile_definition original, datafile_definition *clone, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Moves the clone to original and free up the original. */
void swap_orig_clone(datafile_definition *original, datafile_definition *clone, verbose_definition verbose);


/* Print the description of the gentype to the destination (e.g. stdout) */
void printGenType(int gentype, FILE *destination);

/* return a string with the description of the gentype */
char *returnGenType_str(int gentype);

/* return a string with the description of the file format */
char *returnFileFormat_str(int format);

//START REGION DEVELOP
//START REGION RELEASE

/* Show the content of the history that is stored in the file.
   Returns 1 on success, 0 on error, for instance no history supported. */
int showHistory(datafile_definition datafile, verbose_definition verbose);


//START REGION DEVELOP
//START REGION RELEASE

/* Rewinds the file and skil lines starting with a #. The start
   datastart field in the datafile struct is adjusted. */
int skipallhashedlines(datafile_definition *datafile);


//START REGION DEVELOP
//START REGION RELEASE

// Set the ITRF telescope location XYZ by name of the observatory.
// If verbose is set, the location is printed, with verbose-1 spaces in front of output.
// Return 1 if observatory was recognized, or 0 if not.
int setITRFlocation_by_name(datafile_definition *datafile, char *observatory, verbose_definition verbose);

//START REGION DEVELOP


/* No seeking is done (important to allow for instance multiple files
   are dumped in a single dump file), header is written out
   immediately to the header file pointer defined in the datafile
   struct, so file should already be opened. The datastart variable in
   datafile is set to the position where the data can be found.

   Return 1 on success, 0 on error.
 */
int writePSRSALSAHeader(datafile_definition *datafile, verbose_definition verbose);

//START REGION DEVELOP

/* No seeking is done (important to allow for instance multiple files
   are dumped in a single dump file), header is read immediately from
   the header file pointer defined in the datafile struct, so file
   should already be opened. If nohistory_expected is nonzero, no
   warning is generated if there is no history (but a warning is
   generated if there is).

   Return 1 on success, 0 on error.
 */
int readPSRSALSAHeader(datafile_definition *datafile, int nohistory_expected, verbose_definition verbose);

//START REGION DEVELOP
/* If the data is constant and zero for all polarizations in a given
   single pulse (or subint) in a given frequency channel, then
   iszapped is set to 1, otherwise it is set to 0.

   Returns 1 if successful, or zero if an error occured.
 */
int pulse_iszapped(datafile_definition datafile, long subintnr, int freq, int *iszapped, verbose_definition verbose);

/* Using pulse_iszapped(), find the number of zapped frequency
   channels for a given subint.

   Returns 1 if successful, or zero if an error occured.
 */
int nrchannels_iszapped(datafile_definition datafile, long subintnr, long *nrzapped, verbose_definition verbose);


//START REGION RELEASE

/* ***********************************
   Functions in psrio_opperations.c 
   *********************************** */

// Returns the fold period for subint (counting from zero) in datafile.
// At the moment, only a fixed period folding model is supported, 
// so at the moment, subint is not used.
//
// Returns: 0 = ok
//          1 = not folded, which also results in a warning. (period = 0)
//          2 = error  (period = 0)
int get_period(datafile_definition datafile, long subint, double *period, verbose_definition verbose);

// Returns the fold period for subint (counting from zero) in datafile.
// At the moment, only a fixed period folding model is supported, 
// so at the moment, subint is not used.
double get_tsamp(datafile_definition datafile, long subint, verbose_definition verbose);

// Attempt to convert for instance a pulse longitude list in a fixed sampling time
// Return: 1 = successful, 0 = No changes could be made
int convert_to_fixed_tsamp(datafile_definition *datafile, verbose_definition verbose);

// Get the pulse longitude (in degrees)
double get_pulse_longitude(datafile_definition datafile, long subint, long binnr, verbose_definition verbose);

// Returns the length (in seconds) of the subint (counting from zero)
// in datafile. At the moment, only a fixed subint durations are
// supported, subint is not used.
double get_tsub(datafile_definition datafile, long subint, verbose_definition verbose);

// Returns the length of the observation in seconds.
double get_tobs(datafile_definition datafile, verbose_definition verbose);

// Returns the mid-mjd of the subint (counting from zero)
// in datafile.
long double get_mjd_subint(datafile_definition datafile, long subint, verbose_definition verbose);

// Returns the channel bandwidth (can be negative) for datafile
//
// Return 1 = ok, 0 = failed
int get_channelbandwidth(datafile_definition datafile, double *channelbw, verbose_definition verbose);

// Sets the channel bandwidth (can be negative) for datafile
void set_channelbandwidth(datafile_definition *datafile, double channelbw, verbose_definition verbose);

// Returns the overall bandwidth (can be negative) for datafile, which is the difference between the first and last channel, including their widths. The channels are therefore assumed to be in order.
double get_bandwidth(datafile_definition datafile, verbose_definition verbose);

// Sets the channel bandwidth (can be negative) for datafile
//
// Return 1 = ok, 0 = failed
int set_bandwidth(datafile_definition *datafile, double bw, verbose_definition verbose);

// Returns the centre frequency, i.e. (nonweighted) middle of the overall band.
double get_centre_frequency(datafile_definition datafile, verbose_definition verbose);

// Set the centre frequency, i.e. (nonweighted) middle of the overall band.
void set_centre_frequency(datafile_definition *datafile, double freq, verbose_definition verbose);

/* Returns the midpoint of the frequency  channel number in MHz. */
double get_nonweighted_channel_freq(datafile_definition psrdata, long channel, verbose_definition verbose);

/* Returns the frequency label (if FREQMODE_FREQTABLE) subint/channel number in MHz. Otherwise the same as get_nonweighted_channel_freq(). */
double get_weighted_channel_freq(datafile_definition psrdata, long subint, long channel, verbose_definition verbose);

/* Sets the frequency label (if FREQMODE_FREQTABLE) subint/channel number in MHz. Return 0 on error. */
int set_weighted_channel_freq(datafile_definition *psrdata, long subint, long channel, double freq, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Input Ipulse (NrBins long), output Ipulse2 (NrBins2 long). 
   If noDependencyWarning is set, there is no warning generated 
   when the rebinning is not done with an integer amount (this
   is used when warning is already generated in wrapper function,
   and no warning should be generated for each subint for instance).
   The resulting intensities are weighted such that the mean remains the same.
   Returns 1 if succesfull. */
int rebinPulse(float *Ipulse, long NrBins, float *Ipulse2, long NrBins2, int noDependencyWarning, verbose_definition verbose);

/* Takes data from fin and write it to fout (which will be opened and
   the header should be written out, as well as a history if argc !=
   0). The data is shifted with shift bins and the shift is carried
   over to the next pulse. If circularShift is set, no pulses are lost
   because the last pulse is used to fill in the first pulse. oformat
   can be set to MEMORY_format in which case the data will be written
   to memory. fout should only be declared, but it should be closed
   afterwards.  If verbose2 is set it lets you know what parts of the
   pulse are written out. Return 1 = success */
int continuous_shift(datafile_definition fin, datafile_definition *fout, int shift, int circularShift, char *output_name, int oformat, int argc, char **argv, verbose_definition verbose, int verbose2);

//START REGION DEVELOP
//START REGION RELEASE

/* Get the parallactic angle in radians for the midpoint of the
   specified subintnr (counting from zero and assumed to be of equal
   length). Set subintnr to -1 to obtain the parallactic angle at the
   midpoint of the observation. The ra and dec are assumed to be not
   yet precessed and to be in J2000. If successful (returns 1,
   otherwise 0).*/
int data_parang(datafile_definition data, long subintnr, double *parang, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* By looking at Stokes I to be either just positive or negative, it
   is determined if the baseline is subtracted or not. Returns 0 if
   the baseline appears not to be subtracted, and 1 if the baseline is
   (potentially) subtracted. */
int check_baseline_subtracted(datafile_definition data, verbose_definition verbose);

/* The string text will be parsed and certain keywords will replaced
   with header parameters. The return value is a pointer of the new
   string (memory is allocated, so it should be freed after use), or
   NULL if an error occured. Use the function str_list_replace_keys()
   to list all possible keywords. */
char *str_replace_header_params(datafile_definition data, char *text, verbose_definition verbose);

/* Lists the available keywords for the str_replace_header_params()
   function. nrspaces spaces are printed in front of the output. */
void str_list_replace_keys(int nrspaces);


//START REGION DEVELOP



/* ***********************************
   Functions in rmsynth.c 
   *********************************** */

//START REGION RELEASE

/* Apply the RM synthesis technique to the data for a range of RM's
   (from rm_low to rm_high) divided in nrrmsteps steps. The data is
   written out to rmsynth_array (memory will be allocated when set to
   the NULL pointer). The data is written out as two floats per bin
   per rm trial. The first being the amplitude and the second being
   the phase (in radians). If onpulse != NULL, a check will be done
   such that only selected pulse longitude bins are analysed. Return 0
   = error, 1 success. Verbose level determines nr of spaces before
   output. */
int rmSynthesis(datafile_definition data, float rm_low, float rm_high, float **rmsynth_array, int nrrmsteps, pulselongitude_regions_definition *onpulse, verbose_definition verbose);

/*
  This function takes the output of rmSynthesis(), i.e. rmsynth_array,
  and integrates the results for the pulse longitude bins selected in
  onpulse. The variables nrrmsteps and nrBins should be the same as
  passed to the rmSynthesis command. The output is stored in
  singlespectrum, which should be able to store nrrmsteps floats.
 */
void collapseRMSynthesisArray(float *rmsynth_array, int nrrmsteps, int nrBins, pulselongitude_regions_definition onpulse, float *singlespectrum, verbose_definition verbose);



/* Calculates the expected instrumental responds of the RM synthesis
   method (so intrinsic RM spectrum is convolved with this
   function). rmshift shifts this responds in RM (i.e. a delta
   function intrinsic Faraday depth at this RM will give this observed
   Faraday depth). Determine this function for a range of RM's (from
   rm_low to rm_high) divided in nrrmsteps steps. The data is written
   out to rmsynth_responds (memory will be allocated when set to the
   NULL pointer). If usedouble is set, the data is written in
   rmsynth_responds_double instead. The data is written out as one
   float per rm trial (corresponding to the amplitude). The data is
   assumed to be recorded in nrFreqChan channels each with a bandwidth
   chanbw MHz. The centre frequency of channel 0 is cfreq0.  Return 0
   = error, 1 success. If callmulti=1, previously determined
   intermediate results are stored for successive calls in where only
   the rmshift is different (not checked for). If zero, nothing is
   stored. If set to two, the values are reset and not stored. Verbose
   level determines nr of spaces before output. */
int rmSynthesis_instrument_responds(int nrFreqChan, double chanbw, double cfreq0, double rm_low, double rm_high, float **rmsynth_responds, double **rmsynth_responds_double, int usedouble, int nrrmsteps, double rmshift, int callmulti, verbose_definition verbose);


/*
  For a spectrum returned by collapseRMSynthesisArray (with parameters
  rmmin, rmmax, nrrmsteps), do a fit of the expected instrumental
  responds function rmSynthesis_instrument_responds, which is peaking
  at RM value rm and is offset (vertically) by offset and scaled by
  scale. The data is assumed to be recorded in nrFreqChan channels
  each with a bandwidth chanbw MHz. The centre frequency of channel 0
  is cfreq0.  ftol is the tolerance of the doAmoeba_d() function.
  Returns 0 on success.  */
int rmSynthesis_fitInstrumentalResponds(float *singlespectrum, float rmmin, float rmmax, int nrrmsteps, float *rm, float *offset, float *scale, int nrFreqChan, float chanbw, float cfreq0, float ftol, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* ***********************************
   Functions in psrio_fits.c 
   *********************************** */

void print_fitsio_version_used(FILE *stream);

/* Set a flag that will ignore the weights when reading in PSRFITS files. */
void psrfits_set_noweights(int val);

/* Set a flag that will take the absolute value of the weights when reading in PSRFITS files. */
void psrfits_set_absweights(int val);

/* Set a flag that will ignore read in the weighted frequency in
   PSRFITS files. Note that the centre frequency is changed each time
   you read a pulse from the fits file. */
void psrfits_set_use_weighted_freq(int val);

//START REGION DEVELOP

/* ***********************************
   Functions in psrio_sigproc.c 
   *********************************** */

/* Add the Ingrid flavour to the sigproc format. This adds the pulse
   phase of the start mjd to the header. Set parfile to the location
   of the par file. The name of the parfile should be stay in the same
   memory location, even after calling this funtion. */
void add_IngridFlavour_SigprocAsciiFormat(int val, char *parfile);

/* ***********************************
   Functions in psrio_paswing.c 
   *********************************** */

//START REGION DEVELOP
//START REGION RELEASE
/* Filters out all non-significant PA-points, i.e. those with an PA
   error which is negative. Returns the total number of significant
   points. */
int filterPApoints(datafile_definition *datafile, verbose_definition verbose);


//START REGION DEVELOP
//START REGION RELEASE

/*
It will ignore all points after 360 degrees longitude. This is because
it is important for the pa-fitting that there is only one profile
loaded, and not two rotation periods. If onlysignificantPA is set, all points
without a significant PA are ignored.
*/
int readPPOLfile(datafile_definition *datafile, float *data, int extended, float add_longitude_shift, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/*
Set fptr in datafile to to stdout is possible. If onlysignificantPA is
set only the significant PA's are written out. Can decide to write out
two profiles and or only the significant PA values. If extended == 0,
only the PA and its errorbar are written out. Return value of 1 is OK.
*/
int writePPOLfile(datafile_definition datafile, float *data, int extended, int onlysignificantPA, int twoprofiles, float PAoffset, verbose_definition verbose);


//START REGION DEVELOP
//START REGION RELEASE

/* Provide Q and U data and the onpulse region. If extended is
   non-zero, the total polarization sqrt(Q^2+U^2+V^2) and ellipticity
   + error will be computed as well. Each subint and frequency channel
   is calculated independent of each other. Set normalize if you want
   to normalize intensities to the peak value. Set correctQV to -1 if
   you want to swap sign of Q and V (or to +1 if you don't want to
   make a change). Stokes V is divided by correctV. The pulse
   longitudes are calculated (i.e. data will be in
   TSAMPMODE_LONGITUDELIST, unless nolongitudes in nonzero). loffset
   can be used to shift the data (has no effect in nolongitudes is
   used). Whatever extended is set to, it will always calculate total
   degree of linear polarization and calculates pulse longitudes. By
   default (correctLbias = 0) the median value of the offpulse L is
   subtracted from the data to do the "L bias" correction. If
   correctLbias is set to 1 the much better Wardle & Kronberg debias
   method is used. If set to -1 no bias is subtracted at all. paoffset
   (in degrees) is an offset that will be applied to all PA values. If
   rms_file is not the NULL pointer, the off-pulse rmses are
   determined from this dataset rather than from datafile. In both
   cases onpulse is used to identify the offpulse region in the
   relevant dataset. When applying the rmses on the actual dataset to
   be processed, a scaling can be applied via the variable
   rebin_factor, which indicates how many bins should be summed to get
   a sampling time of datafile. Returns 0 if
   not successful. */
int make_paswing_fromIQUV(datafile_definition *datafile, int extended, pulselongitude_regions_definition onpulse, int normalize, int correctLbias, float correctQV, float correctV, int nolongitudes, float loffset, float paoffset, datafile_definition *rms_file, float rebin_factor, verbose_definition verbose);

//START REGION DEVELOP

int decompose_polarization_modes(datafile_definition datafile, datafile_definition *mode1, datafile_definition *mode2, int method, verbose_definition verbose);

/* 
   Makes a projection of the the Poincare sphere on a
   map of nrx by nry points (idealy with a ration 2:1). The area
   outside the spherical surface will be filled with value
   background. 

   projection = 1 Hammer-Aitoff projection. 
   projection = 2 3D sphere
   projection = 3 Simple linear longitude-latitude map

   For projection 1 and 2 the horizontal axis runs from -2.25 to 2.25
   and the vertical axis from -1.125 to 1.125. For projection 1 the
   poles are at (0, +1) and (0, -1) and the equator extends from (-2,
   0) to (+2, 0). For projection 2 the sphere has a radius of 1. If a
   point lies on the front-side, the sphere is centred at x=-1.125,
   y=0. If the point lies on the far side the point ends on a sphere
   with an x-offset in the opposite direction. For projection 3 the
   horizontal axis runs from -180 to 180 degrees and the vertical axis
   from -90 to 90 degrees.


   The baseline of the data is expected to be subtracted. If weighting
   = 0, the map is basically a count map. If set to 1, weighting
   w.r.t. the polarized power is applied.

   If binnr is non-negative, only the selected pulse longitude bin is
   considered. If binnr == -1, the onpulse region is
   considered. Otherwise all bins are considered.

   Normally longitude 0 and latitude 0 defines the centre of the
   projection, but this grid can be rotated by using rot_long and
   rot_lat (both in degrees), which are added to the longitude and
   latitude lines which are drawn.

   If subtract_pa_data != NULL, the PA-swing from this datafile from
   the data is subtracted. So this should be 1 subit, 1 freq channel,
   converted in PA etc, and the number of bins should match.

   Return 0 on error, 1 on success

*/
int make_projection_map_fromIQUV(datafile_definition datafile, float *map, int nrx, int nry, float background, int binnr, pulselongitude_regions_definition onpulse, int weighting, int projection, float rot_long, float rot_lat, datafile_definition *subtract_pa_data, verbose_definition verbose);

/*
  Makes a PA (or ellipticity) distribution with nrbins pa bins. If
  ellipticity != 0, the ellipticities are used rather than the pa. The
  input data should already be in the form of I, L, V, PA and its
  error bar. The memory for the output data will be allocated and
  dataout can be completely uninitialized. If normalise is set, the
  distribution will be normalised such that the pulse longitude bin
  column with most counts will have a sum of 1.

  In addition, a "PA-mask" can be supplied. This is an array of data,
  with the dimensions as the pa-distribution to be produced, i.e. the
  number of "subint" should be equal to the nrbins variable and the
  number of pulse longitude bins should match those of datain. If this
  pointer to pamask is set to NULL, it will be ignored. Otherwise, the
  input data will be masked in the following way:

  For each PA bin (i.e. subint)/pulse-longitude bin combination,
  pamask will have a recorded value in the first polarization
  channel. If the PA value of a given sample in the input data
  (datain, the pulse stack) matches the value given by pamask_value
  (+- 0.01), the input data will not be altered. Otherwise, that
  particular sample is set to zero. If pamask_value is NaN, the output
  is either the corresponding value as specified in the mask file, or
  zero when the PA is not significant.

  Returns 0 on error
 */
int make_pa_distribution(datafile_definition datain, datafile_definition *dataout, int nrbins, int normalise, datafile_definition *pamask, float pamask_value, int ellipticity, verbose_definition verbose);

/*
  Fit a PA distribution with one or more von Mises functions. The
  memory will be allocated and fit (containing the fit) and resid
  (containing the residual) can be completely uninitialized.

  maxNrComp = The maximum allowed number of functions to fit.
  sigmalimit = the limit up to where components are fitted (5 seems to work reasonably good). 
  precision = The larger precision is set to the finer trials are tried. Something like 2 seems reasonable.

  Returns 0 on error
 */
int fit_pa_distribution(datafile_definition datain, datafile_definition *fit, datafile_definition *resid, int maxNrComp, float sigmalimit, int precision, verbose_definition verbose);


/* P0 = L/sigma (the std. deviation of any of the Stokes parameters,
   they should be equal). L is sqrt(Q^2+U^2) for the pa point, after
   subtracting the bias in L. The 1 sigma errorbar of the pa-point is
   returned in radians. Formula is from Naghizadeh-Khouei & Clarke
   1993.*/
double realisticPAerrorbar(double P0, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* ***********************************
   Functions in psrio_preprocess.c 
   *********************************** */

//START REGION RELEASE

/* Using the original data, construct a profile. Dedispersion,
   de-Faraday rotation etc is taken care off. If stokesI is set,
   Stokes parameters are formed and only StokesI is stored in the
   profile.  The profile should not exist yet. Return 1 = success, 0
   on error. Verbose level determines nr of spaces before output. If
   nocounters is not set, the progresss is shown. */
int preprocess_make_profile(datafile_definition original, datafile_definition *profile, int stokesI, verbose_definition verbose);

//START REGION DEVELOP

/* Estimate the onpulse region automatically of the DEDISPERSED
   profile. It actually finds the offpulse region, so the onpulse
   region is probably a bit too large, so it is best used to use the
   rest as baseline. Set smooth_width to the width (in rotational
   phase) of the smoothing function to determine the minimum in the
   profile, or to a negative value to take the default value (might
   change in future). sigma is the nr of times the profile intensity
   has to be deviating from the mean with the rms of the offpulse
   region to count as a detection (set to negative to take the
   default). minwidthbins removes not selected regions of a maximum of
   this width from the not-selected regions (set to negative to take
   the default). sigma2 is similar to sigma, but after applying the
   minwidthbins option, the not-selected regions are expanded using
   the sigma2 limit (which should be lower than sigma) to ensure there
   is no low-level emission included in not-selected regions. The same
   process is repeated after smoothing the profile with
   smooth_width2. Set nritts to the number of itterations to reach
   final result (set to negative to take the default). Return 1 =
   success, 0 on error. Verbose level determines nr of spaces before
   output. If nocounters is not set, the progresss is shown. */
int preprocess_autoOnpulse(datafile_definition original, float smooth_width, float smooth_width2, float sigma, float sigma2, int minwidthbins, int nritts, pulselongitude_regions_definition *onpulse, verbose_definition verbose);

/* Smooth the data by convolving each subint/freq. channel in Fourier
   space with a tophat function with a full width width (in
   bins). Return 0 = error, 1 = success. Verbose level determines nr
   of spaces before output. */
int preprocess_smooth(datafile_definition original, float width, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Makes a clone with Gaussian noise with the given standard deviation
   added to the data. The clone should not exist yet. Verbose level
   determines nr of spaces before output. If nocounters is not set,
   the progresss is shown. Return 1 = success. */
int preprocess_addNoise(datafile_definition original, datafile_definition *clone, float rms, verbose_definition verbose);

/* Makes a clone with the subints placed in a random order. The clone
   should not exist yet. If fixseed is set, the random number
   generator seed is initialised with a fixed seed to make results
   reproducable. Verbose level determines nr of spaces before
   output. If nocounters is not set, the progresss is shown.

   Return 1 = success, 0 = error. */
int preprocess_shuffle(datafile_definition original, datafile_definition *clone, int fixseed, verbose_definition verbose);

/* Convert coherency parameters into stokes parameters. Verbose level
   determines nr of spaces before output.  Return 1 = success, 0 =
   error. */
int preprocess_stokes(datafile_definition *original, verbose_definition verbose);

/* Convert Stokes parameters into coherency parameters. Verbose level
   determines nr of spaces before output.  Return 1 = success */
int preprocess_coherency(datafile_definition *original, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Rotate stokes parameter stokes1 into stokes parameter stokes2. So
   if stokes1 = 1 and stokes2 = 2 Stokes Q is rotated into U by a
   certain angle in degrees. If angle_array is set to NULL, the value
   angle is used for all pulse longitude bins. angle should be in
   degrees, angle_array in radians. Otherwise, angle_array should
   point to an array specifying the angle for each pulse longitude
   bin. The clone should not exist yet. If inplace is set, the
   variable clone is ignored, and original is modified instead. Set
   subint to the subint you want to modify (counting from zero), or to
   -1 if you want to correct all subints with the specified
   angle. Verbose level determines nr of spaces before output. If
   nocounters is not set, the progresss is shown.

   Return 1 = success, 0 = error. */
int preprocess_rotateStokes(datafile_definition *original, datafile_definition *clone, int inplace, int subint, float angle, float *angle_array, int stokes1, int stokes2, verbose_definition verbose);

//START REGION DEVELOP

/* Subtract a given RVM model (defined by alpha, beta, pa0 and l0 in
   degrees) from the data by rotating the Q & U parameters. The clone
   should not exist yet. If inplace is set, the variable clone is
   ignored, and original is modified instead. Verbose level determines
   nr of spaces before output. If nocounters is not set, the progresss
   is shown.

   Return 1 = success, 0 = error. */
int preprocess_subtractRVM(datafile_definition *original, datafile_definition *clone, int inplace, float alpha, float beta, float pa0, float l0, verbose_definition verbose);

/* Subtract the average (over pulse number) PA from each sample of
   each bin from the data by rotating the Q & U parameters. The clone
   should not exist yet. If inplace is set, the variable clone is
   ignored, and original is modified instead. Verbose level determines
   nr of spaces before output. If nocounters is not set, the progresss
   is shown.

   Return 1 = success, 0 = error. */
int preprocess_subtractAveragePA(datafile_definition *original, datafile_definition *clone, int inplace, verbose_definition verbose);

//START REGION RELEASE

/* Create a clone with nrpulses succesive pulses added together
   (i.e. if set to the total nr of pulses you get a pulse profile). If
   nrpulses is negative each subint is written out -nrpulses times. If
   complete is set, each subint are thrown away to make each subint
   the sum of the same number of input subints. The clone should not
   exist yet. Return 1 = success. Verbose level determines nr of
   spaces before output. If nocounters is set the progress is not
   shown. */
int preprocess_addsuccessivepulses(datafile_definition original, datafile_definition *clone, long nrpulses, int complete, verbose_definition verbose);

/* Rotates each pulse independent of each other in Fourier space based
   on the expected dispersive delay w.r.t. the reference frequency
   specified in datafile_definition. If undo is set, the dispersion
   delay is re-introduced into the data if the data was de-dispersed
   previously. Otherwise, if update is zero, the nothing is done if
   the data is already dedispersed. If update is nonzero, dedispersion
   will be applied again to change the reference frequency from
   freq_ref (-1 or 1e10 = inf freq) to the reference frequency
   specified in datafile_definition.

   Return 0 = error, 1 = command was ignored because of warnings (such
   as data already being dedispersed), 2 = dispersive delay was taken
   out. Verbose level determines nr of spaces before output. */
int preprocess_dedisperse(datafile_definition *original, int undo, int update, double freq_ref, verbose_definition verbose);

/* Compensate for the effect of Faraday rotation with the reference
   frequency determined in datafile_definition. If undo is set, a
   Faraday rotation is applied rather than corrected for. If update is
   set, the reference frequency is changed from freq_ref (-1 = inf
   freq) to what is set in datafile_definition. Cannot do an undo and
   update simultaneously. Verbose determines nr of spaces before
   output.

   If rm_table is set to NULL, the RM from the header is
   used. Otherwise rm_table is interpretted as a list of RM's for each
   pulse longitude bin.

   Return:
     0 = error, 
     1 = ignored because of warning (i.e. already de-Faraday rotated),
     2 = applied de-Faraday rotation. 
*/
int preprocess_deFaraday(datafile_definition *original, int undo, int update, double freq_ref, double *rm_table, verbose_definition verbose);

/*
  Replace the reference frequency with the new frequency new_freq (set
  to -1 or 1e10 if infinite frequency). If the data is not
  dedispersed/de-Faraday rotated, nothing is done except changing the
  header parameter. If the data is dedispersed/de-Faraday rotated, it
  will be adjusted to make it consistent with the new reference
  frequency.

  Return
    0 = error.
    1 = success.
*/
int preprocess_changeRefFreq(datafile_definition *original, double freq_ref_new, verbose_definition verbose);

/* Create a clone with nrfreq succesive frequency channels added
   together (i.e. if set to the total nr of frequency channels you get
   a dedispersed time series). The data should already be dedispersed
   and de-Faraday rotated. If nrfreq is negative each frequency
   channel is written out -nrfreq times. The clone should not exist
   yet. Return 1 = success, 0 on error. Verbose level determines nr of
   spaces before output. This function is almost identical to
   preprocess_addsuccessivepulses. If fzapMask is != NULL, the
   frequency channels that are != 0 in the array are ignored. If
   nocounters is set the progress is not shown. */
int preprocess_addsuccessiveFreqChans(datafile_definition original, datafile_definition *clone, long nrfreq, int *fzapMask, verbose_definition verbose);

/* Create a rebinned clone. The clone should not exist yet. Return 1 =
   success. Verbose level determines nr of spaces before output. If
   nocounters is not set, the progresss is shown. */
int preprocess_rebin(datafile_definition original, datafile_definition *clone, long NrBins, verbose_definition verbose);

/* Subtracts baseline from data. If onpulse != NULL, the onpulse
   region will be excluded from the calculation of the averages. If
   baseline != NULL, the subtracted baseline values are stored to be
   used in preprocess_restore_debase(). If remove_slope = 0, a simple
   offset is subtracted. If 1, a linear gradient is removed from the
   offpulse region. Memory will be allocated to store this
   information.

   Returns 0 on error. */
int preprocess_debase(datafile_definition *original, pulselongitude_regions_definition *onpulse, float **baseline, int remove_shape, verbose_definition verbose);


//START REGION DEVELOP

/* Restores the baseline as was previously subtracted with preprocess_debase(). 
   Returns 0 on error. */
int preprocess_restore_debase(datafile_definition *original, float *baseline, verbose_definition verbose);

//START REGION RELEASE

/* Create a clone with only the selected frequency channel. The clone
   should not exist yet. Return 1 = success */
int preprocess_channelselect(datafile_definition original, datafile_definition *clone, long chanelnr, verbose_definition verbose);

/* Create a clone with only the selected pulses. The clone should not
   exist yet. Verbose level determines nr of spaces before
   output. Return 1 = success */
int preprocess_pulsesselect(datafile_definition original, datafile_definition *clone, long nskip, long nread, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Create a clone with only a multiple of blocksize
   subintegrations. The clone should not exist yet. Verbose level
   determines nr of spaces before output. Return 1 = success */
int preprocess_blocksize(datafile_definition original, datafile_definition *clone, int blocksize, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Rotates each pulse independent of each other using an fft
   algorithm. shiftPhase sets a constant phase offset to be applied to
   each subint. If addslope is set, an additional offset is applied
   which is zero for the first subint, but for following subints it is
   increasing by slope.

   Return 1 = success. Verbose level determines nr of spaces before
   output. If nocounters is set, no progress is shown. */
int preprocess_fftshift(datafile_definition original, float shiftPhase, int addslope, float slope, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Create a clone with only the selected polarization channel. The
   clone should not exist yet. Return 1 = success */
int preprocess_polselect(datafile_definition original, datafile_definition *clone, long polnr, verbose_definition verbose);

//START REGION DEVELOP

/* Create a clone with an inverted bin-axis and subint-axis. The clone
   should not exist yet. Verbose level determines nr of spaces before
   output. Return 1 = success */
int preprocess_invertXY(datafile_definition original, datafile_definition *clone, verbose_definition verbose);

/* Create a clone with an inverted frequency-axis (the order of the
   channels is inverted). Verbose level determines nr of spaces before
   output. The clone should not exist yet. Return 1 = success */
int preprocess_invertF(datafile_definition original, datafile_definition *clone, verbose_definition verbose);

/* Create a clone with an inverted frequency-axis and subint-axis. In
   fact this is just a header opperation. Verbose level determines nr
   of spaces before output. Return 1 = success */
int preprocess_invertFP(datafile_definition *original, verbose_definition verbose);

/* Create a clone with an inverted bin-axis and frequency channel
   axis. The clone should not exist yet. Verbose level determines nr
   of spaces before output. Return 1 = success */
int preprocess_invertFX(datafile_definition original, datafile_definition *clone, verbose_definition verbose);

//START REGION RELEASE
/* Create a clone with a different order of the axis that makes it
   possible to directly plot raw filterbank data. It will have subints
   on the horizontal axis (each nrbins long) and the frequency
   channels on the vertical axis. This block of data for the first
   polarization is followed by the next polarization.

   Return 1 = success. Verbose
   level of verbose determines nr of spaces before output. */
int preprocess_transposeRawFBdata(datafile_definition original, datafile_definition *clone, verbose_definition verbose);
//START REGION DEVELOP

/* Create a clone in which the data outside range defined by onpulse
is thrown away.  The clone should not exist yet. Verbose level
determines nr of spaces before output. If nocounters is not set, the
progresss is shown. Return 1 = success */
int preprocess_gate(datafile_definition original, datafile_definition *clone, pulselongitude_regions_definition onpulse, verbose_definition verbose);

/* Smooth each pulse independent of each other using the fft
   algorithm. If highpass is set, the low frequencies are zapped
   rathern than the low frequencies. Verbose level determines nr of
   spaces before output. If nocounters is not set, the progresss is
   shown. Return 1 = success, 0 on error. */
int preprocess_fftSmooth(datafile_definition original, long nrZapFreqs, int highpass, verbose_definition verbose);

/* Zap each pulse independent of each other using fft
   algorithm. Verbose level determines nr of spaces before output. If
   nocounters is not set, the progresss is shown. Return 1 =
   success */
int preprocess_fftZap(datafile_definition original, long freq1, long freq2, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* If global=0: Normalizes each subint/channel independent of each
   other (based on first polarization channel).

   If global is set: Normalizes each subint/channel with the same
   scale factor based on the peak value found in the first
   polarization channel.

   If onpulse == NULL the whole pulse longitude range is searched for the
   peak value, otherwise only the onpulse region. Return 1 =
   success */
int preprocess_norm(datafile_definition original, float normvalue, pulselongitude_regions_definition *onpulse, int global, verbose_definition verbose);

/* 
   All values > clipvalue are set to the clipvalue. 
   All values < -clipvalue are set to -clipvalue. 
   Return 1 = success */
int preprocess_clip(datafile_definition original, float clipvalue, verbose_definition verbose);

//START REGION DEVELOP

/* Normalizes each subint/channel independent of each other (based on
   first polarization channel) by making the off-pulse rms =
   normvalue. If onpulse == NULL, the whole pulse longitude range is
   used to determine the rms, otherwise only the not selected region
   is used. Return 1 = success */
int preprocess_normRMS(datafile_definition original, float normvalue, pulselongitude_regions_definition *onpulse, verbose_definition verbose);


//START REGION DEVELOP
//START REGION RELEASE

/* Multiplies data values with a certain value after applying an
   offset. Verbose level determines nr of spaces before output. If
   nocounters is not set, the progresss is shown.  Return 1 =
   success */
int preprocess_scale(datafile_definition original, float factor, float offset, verbose_definition verbose);

//START REGION DEVELOP

/* Replaces data with fft. Verbose level determines nr of spaces
   before output. If nocounters is not set, the progresss is
   shown. Return 1 = success */
int preprocess_fft(datafile_definition *original, verbose_definition verbose);

/* Replaces signal with a test signal. Return 1 = success. Verbose
   level determines nr of spaces before output. If nocounters is not
   set, the progresss is shown. */
int preprocess_testsignal(datafile_definition original, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Returns 1 if there are NaN's in the data, or zero otherwise. If
   generate_warning is set, a warning is generated. Verbose level
   determines nr of spaces before output. */
int preprocess_checknan(datafile_definition original, int generate_warning, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Function sets all NaN's to zero. Returns 1 if there are NaN's in
   the data, or zero otherwise. Verbose level determines nr of spaces
   before output. */
int preprocess_removenan(datafile_definition original, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Returns 1 if there are INF's in the data, or zero otherwise. If
   generate_warning is set, a warning is generated. Verbose level
   determines nr of spaces before output. */
int preprocess_checkinf(datafile_definition original, int generate_warning, verbose_definition verbose);

//START REGION DEVELOP

/* Function sets all INF's (or -INF's) to zero. Returns 1 if there are INF's in
   the data, or zero otherwise. Verbose level determines nr of spaces
   before output. */
int preprocess_removeinf(datafile_definition original, verbose_definition verbose);


/* Change the coherency parameters in such a way that a cable-swap is
   undone. Return 1 = success */
int preprocess_swapcables(datafile_definition *original, verbose_definition verbose);

//START REGION RELEASE

/* Compensate for the effect of parallactic angle. If clone is set to
   NULL, the original is modified (if read into memory). Otherwise a
   clone is generated. Verbose level determines nr of spaces before
   output.  If nocounters is not set, the progresss is shown. If undo
   is set, the data is distorted rather than corrected.

   Return 0 = error, 
   1 = command is ignored because already parallactic angle corrected, 
   2 = ignored because of warnings (such as unknown location of telescope), 
   3 = parallactic angle got changed. 

   Note that a clone is generated or data is modified only when the
   return code is 3.
  */
int preprocess_corrParAng(datafile_definition *original, datafile_definition *clone, int undo, verbose_definition verbose);

//START REGION DEVELOP

/* This function takes a pulse of data, npts bins long and convolves
   it with an exponential tail (fast-rise-slow-decay) function. The
   timeconst (e-fold time) is in seconds, as defined at frequency
   reffreq (MHz). The timescale varies with frequency with a
   powerlaw-index pwrlawindex (typically -4 for scattering in
   ISM). Each pulse is treated seperately, so the end of scatter tails
   end up at the start of the same pulse/subint. Return 1 =
   success. Verbose level determines nr of spaces before output. If
   nocounters is set, no progress is shown. */
int convolveScatterTail(datafile_definition original, float timescale, float reffreq, float pwrlawindex, verbose_definition verbose);

/* This function converts the timeconst (e-fold time) in seconds, as
   defined at frequency reffreq (MHz), varying with frequency with a
   powerlaw-index pwrlawindex (typically -4 for scattering in ISM)
   into a scatter time in bins for frequency channel freqchan. */
double scatter_efoldtime_bins(datafile_definition original, long subint, long freqchan, double timescale, double reffreq, double pwrlawindex, verbose_definition verbose);

/* The 0-DM profile is subtracted from the individual non-zapped
   frequency channels. Any dedispersion is removed from the data
   first. This default behaviour is obtained if both
   posdips_stddev_thresh and negdips_stddev_thresh are set to a
   negative value.

   If float negdips_stddev_thresh is set to a positive number, then
   only pulse phases where there is a dip in the intensity exceeding
   negdips_stddev_thresh times the standard deviation will be
   affected.

   posdips_stddev_thresh does something similar to positive spikes.

   if either or both posdips_stddev_thresh and negdips_stddev_thresh
   are positive (or zero), then mode refines the behaviour:

   mode=1: Subtract freq-avrg from each channel
   mode=2: Set each channel to its avrg (zero per def)

   if either or both posdips_stddev_thresh and negdips_stddev_thresh
   are positive (or zero), then maxw refines the behaviour:

   When searching for peaks, if maxw is a positive integer, it will do
   a boxcar search for the most significant peak of various widths in
   the 0-DM profile. Only if the most significant peak is more narrow
   than maxw bins, it will be removed. At either side of the found
   signal extrabins extra bins are removed. In addition, if
   checkdedispersed is non-zero, the width of a signal will be
   compared with the width after dedispersion. Actual pulsar signals
   should become narrower, undispersed signals wider. So only if the
   width after dedispersion becomes wider by minwidening bins (allowed
   to be negative), the signal is regarded to be RFI signal.

   Return 1 = success, 0 on error. 
*/
int preprocess_zero_dming(datafile_definition *original, float posdips_stddev_thresh, float negdips_stddev_thresh, int mode, int maxw, int extrabins, int checkdedispersed, int minwidening, verbose_definition verbose);


/* ***********************************
   Functions in cal.c 
   *********************************** */



void clear_analyticReceiverSolution(analyticReceiverSolution_def *solution);
void printAnalyticReceiverSolution(analyticReceiverSolution_def solution, int indent);

/*
  Writes analytic template to file

  Returns 1 on success, 0 on error
 */
int writeAnalyticReceiverSolution(char *filename, analyticReceiverSolution_def solution);

/*
  Calculates the differential phase (degrees) at frequency freq (MHz)
  according to the analytic receiver solution.

  Return 0 on error, 1 if success
 */
int analyticReceiverSolution_diffphase(analyticReceiverSolution_def solution, double freq, double *phase, verbose_definition verbose);

/*
  Calculates the gain (factor) at frequency freq (MHz)
  according to the analytic receiver solution.

  Return 0 on error, 1 if success
 */
int analyticReceiverSolution_gain(analyticReceiverSolution_def solution, double freq, double *gain, verbose_definition verbose);

/*
  Calculates the differential gain (percentage) at frequency freq (MHz)
  according to the analytic receiver solution.

  Return 0 on error, 1 if success
 */
int analyticReceiverSolution_diffgain(analyticReceiverSolution_def solution, double freq, double *diffgain, verbose_definition verbose);

/*
  Calculates the ellipticity (degrees) at frequency freq (MHz)
  according to the analytic receiver solution. polnr sets which of the
  two polarizations (1 or 2) should be evaluated.

  Return 0 on error, 1 if success
 */
int analyticReceiverSolution_ellipticity(analyticReceiverSolution_def solution, int polnr, double freq, double *ell);

/*
  Calculates the orientation (degrees) at frequency freq (MHz)
  according to the analytic receiver solution. polnr sets which of the
  two polarizations (1 or 2) should be evaluated.

  Return 0 on error, 1 if success
 */
int analyticReceiverSolution_orientation(analyticReceiverSolution_def solution, int polnr, double freq, double *or);

/* Fits a analytic solution to a receiver model. phasemodel defines
   the type of model to fit for the differential phase solution,
   diffgainmodel is for the differential gain, gainmodel, ell1model
   and ell2model is for the ellipticity and or1model.  0 = not fitted,
   otherwise it is a n-1 order polynomial, so 1 = constant, 2 = linear
   function etc.

  In addition to the fit models, a reference models can be specified,
  i.e. a specific possibly complicated frequency dependent shape. See
  also the function refmodellist()

  If rmsOutliers is positive, all points rmsOutliers times the rms
  will be thrown away while deriving the solution. 

  Referencemodel models can be specified, which will be subtracted
  from the solution before fitting is done. If a reference model is
  set to -1, all reference models are explored.

  Verbose determines the number of spaces in front of the output.
  Return 1 = success, 0 on error */
int fitReceiverModel(datafile_definition solution, int gainmodel, int diffgainmodel, int phasemodel, long phasereferencemodel, long diffgainreferencemodel, int ell1model, int ell2model, int or1model, int or2model, double rmsOutliers, analyticReceiverSolution_def *analyticsolution, verbose_definition verbose);


/* Replaces the data with a fake cal signal, so memory should already
   be allocated in original. The data are created as coherency
   parameters rather than Stokes parameters. A list of gains of the
   two polarization channels and phases (in radians) should be
   provided (one value per frequency channel). In addition the
   orientations and ellipticity of leakage can be given (if set to
   NULL these parameters will be ignored). If slope is set, a phase
   gradient (PA-slope) is introduced in the cal signal as function of
   pulse longitude. The default will produce a cal signal with a fixed
   PA as function of pulse longitude. rms sets the rms level of the
   noise, with the signal strength being of the order of the gains
   squared. Two levels are provided which defines a gradient over the
   frequency channels. Set rms_low=rms_high if no gradient is desired.
   If fixseed is set, the noise is reproducable. Set type=1 for the
   standard square-wave signal, type=2 for a gaussian peak. Set
   circular to the fraction of the signal that will be Stokes V. Set
   unpolarised to the fraction of the signal that will be
   unpolarised. By default a pure Q signal is generated (refpa=0). If
   refpa=pi/4 (45 deg), a pure U signal is generated instead. If qw is
   set, the signal is first passed through an imperfect quarter wave
   plate before the receiver parameters are applied. The system is
   characterized with qw_phi (rotation of basis) and the phase delay
   qw_delay, with -45,90 representing an ideal quarter wave plate.
   Verbose level determines nr of spaces before output.  Return 1 =
   successx */
int generateCal(datafile_definition *original, float *gain1, float *gain2, float *phases, int slope, float *el0, float *el1, float *or0, float *or1, float rms_low, float rms_high, int fixseed, int type, float circular, float unpolarised, float refpa, int qw, float qw_phi, float qw_delay, verbose_definition verbose);

/* Applies the receiver model to the datafile. The datafile and the
   solution should be read into memory. If reconstruct is set, the
   intrinsic Stokes parameters are reconstructed using the receiver
   parameters. If not, the data is distorted. If ignoreGain is set,
   the gain parameters are ignored, which might be useful to prevent
   the data to be dominated by a few channels. If ignore_dodgy_leakage
   is set, it makes the Mueller matrix zero (i.e. it will zap data) if
   the normalisation factor in the leakage matrix is abnormally
   high. This will only be in effect together with the reconstruct
   flag. If normalise is set, the matrix is divided by the value of
   the first element, i.e. after applying the Mueller matrix to the
   data the rms of Stokes I in the noise is unchanged, assuming that
   the rms of the four Stokes parameters are equal to start off
   with. Normalise is only identical to setting the gain to 1 for an
   ideal receiver. Return 1 = success */
int applyReceiverModel(datafile_definition *datafile, datafile_definition solution, int reconstruct, int ignoreGain, int ignore_dodgy_leakage, int normalise, verbose_definition verbose);

/* Prints the Mueller matrix constructed from the receiver model.  Set
   fd_type=1 for a linear basis, 2 for a circular basis. If singlerow
   is set, the matrix elements appear on a single row. The solution
   should be read into memory. If reconstruct is set to 1, the
   intrinsic Stokes parameters are reconstructed using the receiver
   parameters. If set to zero, the data is distorted. If reconstruct
   is set to 2, the Mueller matrix is multiplied with its inverse, to
   test if the unity matrix is obtained. If ignoreGain is set, the
   gain parameters are ignored. If normalise is set, the matrix is
   divided by the value of the first element, i.e. after applying the
   Mueller matrix to the data the rms of Stokes I in the noise is
   unchanged, assuming that the rms of the four Stokes parameters are
   equal to start off with. Normalise is only identical to setting the
   gain to 1 for an ideal receiver. If stat is set, some statistics
   are returned as well. Return 1 = success, 0 = error. */
int printMuellerFromReceiverModel(datafile_definition solution, int reconstruct, int ignoreGain, int normalise, int fd_type, int singlerow, int stat, verbose_definition verbose);

/*
   0 - Gain
   1 - Differential gain [%%]
   2 - Differential phase [deg]
   3 - Elipticity 0 [deg]
   4 - Orientation 0 [deg]
   5 - Elipticity 1 [deg]
   6 - Orientation 1 [deg]
   7 - Chi-square
   8 - Nfree

   Returns the polarization nr in which the parameter is
   stored (value >= 0). On error the return value is negative:
      -1 if parameter is not part of the solution (consistent with gentype).
      -2 According to gentype the parameter exist, but it doens't match the nr of stored parameters
      -3 The input parameter is not recognized (beyond range?)
*/
int get_polnr_in_receiversolution(datafile_definition solution, int parameter);

/* Extract a parameter from the receiver model. Provide an
   analyticsolution != NULL if you want to subtract this from the
   data. The data is returned in xdata and ydata and yerrorbar, for
   which an array of nrpoints is allocated (so you should free
   it). xdata will contain the frequency, ydata the requested
   parameter. Only frequency channels with a valid solution are
   returned. The parameters are defined to be:

   0 - Gain
   1 - Differential gain [%%]
   2 - Differential phase [deg]
   3 - Elipticity 0 [deg]
   4 - Orientation 0 [deg]
   5 - Elipticity 1 [deg]
   6 - Orientation 1 [deg]
   7 - Chi-square
   8 - Nfree

   verbose determines the number of spaces before output.
   Return 1 = success, 0 on error 
*/
int extractReceiverParameter(int parameter, datafile_definition solution, long subintnr, analyticReceiverSolution_def *analyticsolution, float **xdata, float **ydata, float **yerrorbar, long *nrpoints, verbose_definition verbose);

// Shows the available system respondses which can be subtracted from polarization calibration receiver responds
void refmodellist();



/* Shows the receiver model on a pgplot device. The solution should be
   read into memory. Neat sets the style of the output:

   0 = Each parameter is shown seperately
   1 = More like a psrchive style output
   2 = No chi^2 and nfree, and all parameters on a single page

   If fixscale is set, the differential phase, orientation and
   ellipticity range is set to -90..90 and the differential gain range
   to -100..100. Set analyticsolution = NULL if you don't want an
   analytic solution to be superimposed. Set device to a pgplot
   device, or set to NULL if you want the user to specify the
   device. Set first when you want to open the device and last to
   close it. Multiple files can be plotted in this way. When plotting
   one solution, set both first and last. Multiple solutions can be
   show. If overlay is non-zero, the receiver solutions are
   overplotted in the same graph. If the title is set to NULL, the
   filename is used instead. The title supports header keywords set by
   str_replace_header_params(). If showmodel is set, the type of model
   used to fit the parameters is identified in the plots via
   text. Return 1 = success */
int showReceiverModel(datafile_definition *solution, long nrsolutions, int neat, int fixscale, analyticReceiverSolution_def *analyticsolution, char *device, int first, int last, int overlay, char *title, int showmodel, verbose_definition verbose);

/* Given nrsolutions input solutions, construct a single overall
   average receiver solution (combined, not yet initialised with
   cleanPSRData).

   type:
   1: For each channel, store that solution with the smallest chi2

   Return 1 = success */
int combineReceiverModels(datafile_definition *solution, long nrsolutions, int type, datafile_definition *combined, verbose_definition verbose);

/* Take the input receiver model and replace the parameters with the
   analytic solution.

   Return 0 = Error, 1 = success 
*/
int makeReceiverModelAnalytic(datafile_definition *solution, analyticReceiverSolution_def *analyticsolution, verbose_definition verbose);

/* The conversion will only be done if poltype == 2 (i.e. coherency
   parameters). fd_type indicates if the receiver is linear or
   circular. */
void convertCoherencyToStokes(float *C, int poltype, int fd_type);

/* The conversion will only be done if poltype == 2 (i.e. set to
   coherency parameters). fd_type indicates if the receiver is linear
   or circular. */
void convertStokesToCoherency(float *C, int poltype, int fd_type);

/* This function fits the data with the template. Note that both
datasets could be destroyed, as data might need to be rebinned and
shifted and dedispersed etc. The onpulse region can be set with
onpulseregion for noise calculation. If set to NULL, the user will be
asked to select it. This will also happen when the number of regions
is 0, but in that case the region will be returned (useful when
processing multiple files). onpulseregion_fit is similar, but then for
the region that is used to do the actual fit. If aligned is set, the
template and observation are not aligned using a cross-correlation.
If initialguess != NULL, this receiver model is used to set the
initial guess. The solution is written as a fits file to
outputfilename. argc and argv should be provided to allow writing of a
history line to the output file. If debug is set, a lot of extra
information is produced, including plots (the before/after plots are
sent to beforeafter_device). Instead, if profileresid_device != 0, the
plots are sent to this pgplot id instead, which allows plots from
subsequent calls to cal_mtm() to be appendend. If leakage is set, the
leakage terms will be fitted for as well, otherwise just the gain
imbalance. The non-zero elements in the array called fixed (which has
to contain 7 elements) indicates which parameters have to be fixed to
the initial guess value. If mtm_inverse is set, the observation is
corrected during the fitting process to match the template. This is in
general a bad idea, because if the S/N of the data is low, the
amplitude of the Stokes parameters is underestimated. By default the
template is distorted to match the observation. If analytic_gain is
set, the gain is not treated as a free parameter, but it is fixed to
whatever is needed to make the average onpulse Stokes I match the
template. This resolves the bias which makes the amplitude of the
observation lower than the template (in low S/N situations), but there
is still a bias in the determined receiver solution. If bruteforce is
zero, only 1 initial value for each receiver parameter is
tried. Otherwise the number of trials are (1+2*bruteforce)**(nr fit
parameters). Hence this variable takes bruteforce extra trials at each
side of the central value, which is the ideal receiver value. Note
that the number of trials are increasing extremely rapidly, so be
careful! A somewhat quicker way is to set extremes, in which case for
all angles 0 and the extreme values are tried, without any
intermediate values. If sequencenr is set to a negative number, a
dumpfile will be generated, containing all the required information to
run this process in multiple instances, each instance processing part
of the band. If totsequencenr is set to 0, cal_mtm processes the band
in one go without generating a dump file. Set sequencenr to
1...totsequencenr to process part of the band. When part of the band
is processed, the datafile, template and initialguess should not be
allocated yet, as data will be read from dumpfile and memory will be
allocated. If nocounters is set, no counters are shown. If
ignoreParallacticAngle is set, all parallactic angle corrections are
ignored. This function returns 0 when there is an error. */
int cal_mtm(datafile_definition *datafile, datafile_definition *template, pulselongitude_regions_definition *onpulseregion, pulselongitude_regions_definition *onpulseregion_fit, int aligned, datafile_definition *initialguess, int leakage, int *fixed, int analytic_gain, int mtm_inverse, int bruteforce, int extremes, char *outputfilename, int argc, char **argv, char *dumpfilename, int sequencenr, int totsequencenr, int ignoreParallacticAngle, verbose_definition verbose, char *beforeafter_device, int profileresid_device);

/* Construct the Mueller matrix for the given receiver parameters. If
   reconstruct is set, the Mueller matrix will recover the intrinsic
   Stokes parameters from the recorded Stokes parameters. If
   reconstruct is not set it will distord the Stokes parameters in the
   same way the receiver would do. If ignore_leakage is set, the
   leakage (ellipticities/orientations) are ignored, which makes this
   function faster. If ignore_dodgy_leakage is set, it makes the
   Mueller matrix zero (i.e. it will zap data) if the normalisation
   factor in the leakage matrix is unrealisticly high. This will only
   be in effect together with the reconstruct flag. If normalise is
   set, the matrix is divided by the value of the first element,
   i.e. after applying the Mueller matrix to the data the rms of
   Stokes I in the noise is unchanged, assuming that the rms of the
   four Stokes parameters are equal to start off with. Normalise is
   only identical to setting the gain to 1 for an ideal receiver.

   Return:
   0 = dodgy solution, M=0
   1 = solution accepted, M is set   
*/
int constructMueller(float M[4][4], double G, double gamma, double phase, double el0, double el1, double or0, double or1, int fd_type, int reconstruct, int ignore_leakage, int ignore_dodgy_leakage, int normalise);

/* S must be 4 Stokes parameters. The vector S is multiplied with the
   Mueller matrix and S is replaced with the result. */
void applyMueller(float M[4][4], float *S);

//START REGION DEVELOP
//START REGION RELEASE

/* ***********************************
   Functions in fft.c 
   *********************************** */

void print_fftw_version_used(FILE *stream);

//START REGION DEVELOP
/* Based on Russells arbitrary_bin_shift function
   (bin_shift.cc). Returns 0 if failed. */
//START REGION RELEASE
int rotateSinglepulse(float *data, int npts, float epsilon, verbose_definition verbose);

//START REGION DEVELOP

/* Smooth each pulse independent of each other using fft algorithm by
   zapping nrSmoothFreq frequency bins. If highpass is set, the low
   frequencies are zapped rathern than the low frequencies.  Returns 0
   if failed, otherwise 1. Verbose determines the number of spaces
   before output. */
int fftSmooth(float *data, long npts, int highpass, long nrSmoothFreq, verbose_definition verbose);

/* Zaps given range of frequency bins.  Returns 0 if failed, 
otherwise 1.
*/
int fftZap(float *data, long npts, long freq1, long freq2, verbose_definition verbose);

/* Replace data by its (real) fft with length npts/2. Last sample is
   ignored if the npts is odd. Returns 0 if failed, otherwise 1. */
int replacefft(float *data, long npts, verbose_definition verbose);

//START REGION RELEASE
/* Calculates the cross-correlation using fft's. The two data-sets
   should have the same length (ndata), not necessarily a power of
   two. You could use the function crosscorrelation_fft_padding
   instead, which allows you to zero-pad the data first. The
   cross-correlation function is returned in cc. Return value is 0 if
   there is an error. */
int crosscorrelation_fft(float *data1, float *data2, int ndata, float *cc, verbose_definition verbose);

/* Calculates the length of the output of the function
crosscorrelation_fft_padding (i.e. what cclength will be set to). This
includes the zero-padding in the calculation.  */
int crosscorrelation_fft_padding_cclength(int ndata, int extrazeropad);

//START REGION DEVELOP



//START REGION RELEASE
/* Calculates the cross-correlation using fft's. The two data-sets
   should have the same length (ndata). If extrazeropad is set, at
   least this amount of zero's are pasted after the data set. If the
   length of the data (plus any extra zeropadding) is not a power of
   two, it will be zero-padded further to a length that is a power of
   two. The cross-correlation function is returned in cc (memory will
   be allocated), which is cclength points long. Return value is 0 if
   there is an error. Verbose determines the number of spaces before
   output. */
int crosscorrelation_fft_padding(float *data1, float *data2, int ndata, int extrazeropad, float **cc, int *cclength, verbose_definition verbose);
//START REGION DEVELOP

/*
  Determine delay between data1 and data2, which should each be nrbins
  long arrays. The delay is found by determining the location of the
  maximum in the cross-correlation function, which will have a
  precision of 1 bin. The delay is returned in phase. Verbose
  determines the number of spaces before output.

  Return 0 on error, 1 if successful.
 */
int crosscorrelation_find_delay(float *data1, float *data2, long nrbins, double *delay, verbose_definition verbose);

/* Calculates the delay using the Fourier phase gradient method
   descibed in Taylor 1992. The two data-sets should have the same
   length (ndata), not necessarily a power of two. The delay is
   returned as a delay in phase. Verbose determines the nr of spaces
   before output.

   Return value:
      0 - Critical error occured
      1 - Converged
      2 - Not converged, not no critical error occured
*/
int delay_Fourier_phase_gradient(float *data1, float *data2, int ndata, double *delay, verbose_definition verbose);

/* This function takes a pulse of data, npts bins long and convolves
   it with an exponential tail (fast-rise-slow-decay) function. The
   timeconst is in bins. Note that the end of the scatter tail ends up
   at the start of the pulse.

   Returns 0 if failed, otherwise 1. */
int convolveScatterTail_singlepulse(float *data, int npts, float timeconst, verbose_definition verbose);

/* This function takes a pulse of data, npts bins long and convolves
   it with tophat function of full width timeconst (in bins). Verbose
   level determines nr of spaces before output.

   Returns 0 if failed, otherwise 1. */
int tophat_smooth(float *data, int npts, float timeconst, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* ***********************************
   Functions in spectra.c 
   *********************************** */


/* 
   Note there is a wrapper function calcLRFS_padding() which allows
   you to do zero-padding. Calculate the LRFS of data, which is an
   array of nrx by nry points. The powers of the fft's of length
   fft_size are added together in lrfs, which should be an array of
   nrx by fft_size/2+1 points (which corresponds to the frequencies
   0,1/fft_size,2/fft_size, ..., 0.5). If subtractDC is specified, the
   DC component is subtracted from the LRFS. If
   avrg_offpulse_lrfs_power != NULL, it will be set to the average
   offpulse power in the obtained lrfs (so the average power per
   spectral/longitude bin in *lrfs. phase_track should be nrx points
   long (if calcPhaseTrack is specified), which contains the
   vertically collapsed subpulse phase for each bin in degrees. The
   range in fluctuation frequencies specified by freq_min and freq_max
   is used for the track calculation. If phase_track_phases is not
   NULL, the phase offsets between the phase_tracks in subsequent
   fft_blocks are stored (so it should be an array of nry/fft_size
   floats). If track_only_first_region is zero, it means that the full
   onpulse region is used to align the phases obtained from the
   individual blocks. If set to 1, only the first on-pulse region is
   used to do the correlations, while the rest of the regions are
   still used for the noise subtraction of the LRFS. The noise
   subtraction is based on all not-selected pulse-longitude bins. If
   calcsubpulseAmplitude is set (subpulseAmplitude should have a
   length nrx otherwise), the power of the drifting subpulses is
   stored for each longitude bin (as the stddev in selected frequency
   range divided by the maximum stddev in the profile). These are the
   averages of the different spectral channels for a given block, not
   weighted by anything. The final subpulse phase track is a weighted
   average however. Set the mask_freqs flag to mask out everything
   outside the frequency range. var_rms is the rms (after subtracting
   the mean, using the separate blocks of data rather than the final
   lrfs for which the different blocks are summed) of the samples in
   the lrfs in the off-pulse region.
//START REGION DEVELOP
   If inverseFFT is set it calculates the inverse transform. In
   that case the data is overwritten. I'm not sure if the
   normalization is done correctly. You might be better of to select a
   single bin.
//START REGION RELEASE 
   If regions is defined (on-pulse), than the off-pulse noise
   contribution is subtracted and the rms of the offpulse rms is
   calculated (used by calcModindex). If verbose is set it shows the
   number of blocks that is calculated. If argc > 0 and argv != NULL,
   it will search the command line for the -p3zap option and zap these
   regions from the lrfs (can be specified in cpp or in bins). The 
   debug option gives more output.

   0 = nothing done 
   1 = ok
*/
int calcLRFS(float *data, long nry, long nrx, unsigned long fft_size, float *lrfs, int subtractDC, float *avrg_offpulse_lrfs_power, float *phase_track, float *phase_track_phases, int calcPhaseTrack, float freq_min, float freq_max, int track_only_first_region, float *subpulseAmplitude, int calcsubpulseAmplitude, int mask_freqs, int inverseFFT, pulselongitude_regions_definition *regions, float *var_rms, int argc, char **argv, verbose_definition verbose);

//START REGION DEVELOP

// This function is the same as calcLRFS, but it takes as an extra argument the number of subintegrations to add to the end of the dataset filled with zeros. The only extra parameter is nr_zero_pad_subints_to_add. So this is a wrapper function that does zero-padding.
int calcLRFS_padding(float *data, long nry, long nrx, unsigned long fft_size, float *lrfs, int subtractDC, float *avrg_offpulse_lrfs_power, float *phase_track, float *phase_track_phases, int calcPhaseTrack, float freq_min, float freq_max, int track_only_first_region, float *subpulseAmplitude, int calcsubpulseAmplitude, int mask_freqs, int inverseFFT, pulselongitude_regions_definition *regions, float *var_rms, int argc, char **argv, verbose_definition verbose, int nr_zero_pad_subints_to_add);

//START REGION DEVELOP
//START REGION RELEASE

/* Calculate std deviation (sigma) and the modulation index. The
   profile will be normalized in the process. If regions is set it
   will calculate error-bars as well (an array for rms_sigma and
   rms_modind). Calculate the profile for a integer number of
   fft-blocks and set nrpulses to a multiple of fft_size. You can pass
   on var_rms from calcLRFS to obtain a higher precision (and smaller
   numbers) in the error-bar calculation if there are baseline
   variations. Not sure if it is better, but it is more consistent
   with Russell's software. Not sure if it is really working well. In
   fact, I don't believe the errorbars are correct in Russells
   software either (I simulated a 2D sinusoid convolved with a top-hat
   function as a test). You should do bootstrapping to obtain reliable
   errorbars. If *avrg_offpulse_lrfs_power is != NULL, then this
   value as provided by calcLRFS is used to report the offpulse
   normalised variance as spectrally determined.
*/
void calcModindex(float *lrfs, float *profile, long nrx, unsigned long fft_size, unsigned long nrpulses, float *sigma, float *rms_sigma, float *modind, float *rms_modind, pulselongitude_regions_definition *regions, float var_rms, float *avrg_offpulse_lrfs_power, verbose_definition verbose);

//START REGION DEVELOP

/* Calculate std deviation (sigma) and the modulation index directly
   from the pulse stack. No errors are calculated and no baseline is
   subtracted. If regions is defined the variance of the off-pulse
   region is subtracted before calculating the modulation index. */
void calcModindexSimple(float *pulsestack, long nrx, unsigned long nrpulses, float *sigma, float *modind, pulselongitude_regions_definition *regions, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* 
   Calculate the 2DFS of data, which is an array of nrx by nry
   points. The on-pulse region should be defined and should have a
   length which is a power of two. The region number to calculate the
   2DFS for can be specified. Only multiples of fft_size are used in
   the y-direction). The 2DFS's of length fft_size are added together
   in 2dfs, which should be an array of nrx by 1+fft_size/2
   points. All the defined (on-pulse) regions are ignored for the
   off-pulse noise subtraction. The off-pulse subtraction will be done
   using an off-pulse region with a with equal to the selected
   region. This off-pulse region start at bin 0, if it fits. Otherwise
   the start bin of the off-pulse region is increased until a suitable
   off-pulse region is identified. If there is no suitable off-pulse
   region, no off-pulse noise subtraction will be applied. If verbose
   is set it shows the number of blocks that is calculated (unless
   dontshowprogress is set).

   0 = nothing done 
   1 = ok
*/
int calc2DFS(float *data, long nry, long nrx, unsigned long fft_size, float *twodfs, pulselongitude_regions_definition *onpulse, int region, verbose_definition verbose);

//START REGION DEVELOP
/* 
   Calculate the HRFS of data, which is an array of nrx by nry
   points. The data is analysed per block of pulses_per_block pulses and the
   resulting power spectra are added. The resulting spectrum is
   hrfs_length points in length (memory will be allocated). If regions
   is defined (on-pulse), than the off-pulse bins are set to zero
   before doing calculating the hrfs.

   0 = nothing done 
   1 = ok
*/
int calcHRFS(float *data, long nry, long nrx, long pulses_per_block, float **hrfs, long *hrfs_length, pulselongitude_regions_definition *regions, verbose_definition verbose);
//START REGION DEVELOP

/* 
   Calculate the LRAC (longitude resolved auto-correlation) of
   data. If lrac_keepbothhalves is set, both positive and negative
   lags are stored, which is an array of nrx by cclength points (as
   returned by crosscorrelation_fft_padding_cclength() when adding at
   least nry extra zero-padding zero's). If lrac_keepbothhalves is
   zero, only zero lag + positive lags are stored (the return value of
   crosscorrelation_fft_padding_cclength()/2+1 points). This length is
   returned in cclength. The data will be zero-padded in the nry
   direction to make it a power of 2. If remove_wf is set, the
   triangular shaped effect because of the window function is removed
   from the AC. If remove_zerolag the zero lag coefficients (often
   respondsible for a large spike) is set to zero. The debug option
   gives more output.

   0 = error
   1 = ok
*/
int calcLRAC(float *data, long nry, long nrx, float *lrac, int remove_wf, int remove_spike, int lrac_keepbothhalves, long *cclength, verbose_definition verbose);

//START REGION DEVELOP

/*
  If argc > 0 and argv != NULL,
   it will search the command line for the -p3zap option and zap these
   regions from the lrfs (can be specified in cpp or in bins). The 
   debug option gives more output.

 Return 0 on error
*/
int p3classify(float *data, long nry, long nrx, long padding_up_to, long nblock, long min_length, float cpp1, float cpp2, float *output, pulselongitude_regions_definition *regions, int argc, char **argv, verbose_definition verbose);




//START REGION DEVELOP
/* 
   Calculate the longitude resolved cross-correlation map of the data
   for the specified lag. The map is written out to lrcc, which should
   be nrx*nrx floats in size. If noSubtractMean is specified, the mean
   profiles is not subtracted from the LRCC. 

   lrcc_ij = sqrt(sum_i sum_j sum_t stack_i(t)*stack_j(t+lag))

   0 = nothing done 
   1 = ok
*/
int calcLRCC(float *data, long nry, long nrx, int lag, float *lrcc, int noSubtractMean, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE


//START REGION RELEASE

/* ***********************************
   Functions in fold.c 
   *********************************** */

/* Folds the data (an array of nrx phase bins by nry pulses) on a
   given p3 value (expressed in pulse periods). The resulting map
   (which should already be allocated) has is nrx by nr_p3_bins points
   in size. The resulting map is normalised by the number of times a
   specific bin in the P3 cycle is filled in. If refine is set the
   routine will align the blocks first by doing a cross-correlation in
   order to try to correct for P3 changes etc. If refine is larger
   than 1, it will first produce a template, which then will be used
   to calculate more accurate offsets etc. The larger the value of
   refine, the better the result should be. Set cyclesperblock to set
   the number of P3 cycles that are folded first before attempting to
   align the blocks (number of pulses folded first is
   foldp3*cyclesperblock rounded down). If onpulse is set, it will
   only use the onpulse region to align the blocks. By default each
   pulse is split over two p3 phase bins, depending how accurate it
   folds on the integer number of p3 phase bins. However, if
   smoothWidth is positive, then the for a given pulse the weight for
   each p3 phase is set by a gaussian with this width (in p3 phase
   bins). As another alternative, if noSmooth is set, all power ends
   up in a single bin. If slope is nonzero (in degrees per bin) a
   subpulse phase slope is subtracted. Offset (in degrees) is
   added. If offsetFileName is set, the found offsets are written out
   to a file (readwriteOffsets = 1) or read in and used
   (readwriteOffsets = 2)

   0 = nothing done 
   1 = ok
*/
int foldP3(float *data, long nry, long nrx, float *map, int nr_p3_bins, float foldp3, int refine, int cyclesperblock, int noSmooth, float smoothWidth, float slope, float subpulse_offset, pulselongitude_regions_definition *onpulse
//START REGION DEVELOP
	   , char *offsetFileName, int readwriteOffsets, char *pulse_list_filename, int writeP3foldPulseList
//START REGION RELEASE
	   , verbose_definition verbose);


//START REGION DEVELOP
float plotCarousel(carousel_parameter_list plist, char *filename_ptr, float *output_image, datafile_definition datain);
//START REGION RELEASE


/* ***********************************
   Functions in astronomy.c 
   *********************************** */

/* Returns the DM delay in seconds at freq freq relative to freq_ref
   (so returns a positive value is freq < freq_ref). If inffrq is
   nonzero, the reference frequency is set to infinity. The
   frequencies are in MHz and the DM in pc/cm^3. */
long double calcDMDelay(long double freq, long double freq_ref, int inffrq, long double dm);

/* Returns the amount of Faraday rotation (of the position angle in
   radians) according to the RM at frequency freq relative to
   freq_ref. If inffrq is nonzero, the reference frequency is set to
   infinity. The frequency is in MHz and the RM in radians/meter^2. */
float calcRMAngle(float freq, float freq_ref, int inffrq, float rm);

// Longitude/latitude are GRS80 geodetic coordinates in radians
// Height and ITRF X,Y,Z coordinates are in meters
void tempo2_ITRF_to_GRS80(double obs_X, double obs_Y, double obs_Z, double *longitude, double *latitude, double *height);

// Longitude/latitude are GRS80 geodetic coordinates in radians
// Height and ITRF X,Y,Z coordinates are in meters
void
tempo2_GRS80_to_ITRF(double longitude, double latitude, double height, double *obs_X, double *obs_Y, double *obs_Z);

// returns the geocentric (ITRF based) longitude of the observatory in radians
double observatory_long_geocentric(datafile_definition datafile);

// returns the geocentric (ITRF based) latitude of the observatory in radians
double observatory_lat_geocentric(datafile_definition datafile);

// returns the geodetic (GRS80 based) longitude of the observatory in radians
double observatory_long_geodetic(datafile_definition datafile);

// returns the geodetic (GRS80 based) latitude of the observatory in radians
double observatory_lat_geodetic(datafile_definition datafile);

// returns the geodetic (GRS80 based) height of the observatory in meters
double observatory_height_geodetic(datafile_definition datafile);

//START REGION DEVELOP

void mjd2date(long double mjd, int *year, int *month, int *day, int *hour, int *minute, float *seconds);

//START REGION DEVELOP
//START REGION RELEASE

/* 
  Type 1: YYYY-MM-DDXHH:MM:SS

  X = separator
  precision gives the number of decimal places in the seconds.
*/
void mjd2dateString(long double mjd, char *string, int precision, int type, char *separator);

//START REGION DEVELOP
//START REGION RELEASE

/* Converts a string with hours, minutes and seconds into a floating
   point number (in hours).*/
void converthms(char *hms, double *h);
//START REGION DEVELOP
void converthms_ld(char *hms, long double *h);
void converthms_f(char *hms, float *h);

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
void converthms_string(char *hms, long double number, int precision, int type);
//START REGION DEVELOP
//START REGION RELEASE

/* Given the geodetic longitude and latitude of the observatory
   (longitude and latitude) and the position of the source (right
   ascension and declination) and the mjd, return the parallactic
   angle. All angles are in radians. If precess is set, the ra and dec
   will be corrected for precession, nutation and aberration (using
   J2000). */
double calc_parang(double longitude, double latitude, double ra, double dec, double mjd, int precess);

//START REGION DEVELOP
//START REGION RELEASE

/* Set system to B or J for B1950 or J2000 coordinates. The
   coordinates are defined at 1950.0 or 2000.0 depending on the
   specified system and will be transformed to epoch mjd1. The input
   ra and dec (radians) will be transformed. If nutation is set, the
   additional effect of nutation is corrected for. Likewise aberration
   sets the correction for Annual aberration. Returns 0 if system was
   not recognized. */
int calc_precess_nut_ab(char system, double mjd, double *ra, double *dec, int nutation, int aberration, int verbose);

//START REGION DEVELOP

/* Converts a time in hours into a string with hours, minutes and
   seconds. If precision is set to zero only minutes are produced. */
/*void gethourstring(double t, char *txt, int precision);*/

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
void make_projection_map(float *map, int nrx, int nry, float background, float dlongitude, float dlatitude, float (*funk)(float, float), int projection);

/*
  Given the longitude and latitude (in radians), calculate the
  Hammer-Aitoff projection Cartesian coordinates x and y. The range of
  x is between -2 and 2, and that of y between -1 and 1. The angles
  dlongitude and dlatitude correspond to a rotation of the sphere such
  that positive angles correspond with the equator moving towards the
  Z-pole and a rotation in the positive longitude direction. 
*/
void projectionHammerAitoff_xy(float longitude, float latitude, float dlongitude, float dlatitude, float *x, float *y);

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
int projectionHammerAitoff_longlat(float x, float y, float dlongitude, float dlatitude, float *longitude, float *latitude);

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
int projection_sphere_xy(float longitude, float latitude, float dlongitude, float dlatitude, float *x, float *y, float *weight);

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
int projection_sphere_longlat(float x, float y, float dlongitude, float dlatitude, float *longitude, float *latitude);

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
int projection_longlat_longlat(float x, float y, float dlongitude, float dlatitude, float *longitude, float *latitude);

/* Given the longitude and latitude (in radians), calculate the
   corresponding Cartesian coordinates x and y of a simple linear
   longitude-latitude map projection (these will be in degrees to be
   compatable with drawSphericalGrid()). This projection is very
   simple. The angles dlongitude and dlatitude correspond to a
   rotation of the sphere such that positive angles correspond with
   the equator moving towards the Z-pole and a rotation in the
   positive longitude direction.

*/
void projection_longlat_xy(float longitude, float latitude, float dlongitude, float dlatitude, float *x, float *y);

//START REGION DEVELOP
//START REGION RELEASE

/* ***********************************
   Functions in angles.c 
   *********************************** */

//START REGION DEVELOP

/* Calculates the cross-product */
void crossprod(float a1, float a2, float a3, float b1, float b2, float b3, float *c1, float *c2, float *c3);
void crossprod_d(double a1, double a2, double a3, double b1, double b2, double b3, double *c1, double *c2, double *c3);

// answer = axb
void crossprod_array(float *a, float *b, float *answer);
void crossprod_array_d(double *a, double *b, double *answer);

double length_vector_array(float *a);
double length_vector_array_d(double *a);

double dotproduct_array(float *a, float *b);
double dotproduct_array_d(double *a, double *b);

//START REGION RELEASE
/* Derotates an angle between 0 and 360 degrees */
float derotate_deg(float a);
//START REGION DEVELOP
double derotate_deg_double(double a);
long double derotate_deg_longdouble(long double a);

/* Derotates an angle between 0 and 2 pi radians */
double derotate_rad_double(double a);
long double derotate_rad_longdouble(long double a);

//START REGION RELEASE
/* Derotates an angle between 0 and 180 degrees */
float derotate_180(float a);
double derotate_180_double(double a);
/* Derotates an angle between 0 and pi radians */
double derotate_180_rad_double(double a);
/* Derotates an angle between 0 and 90 degrees */
float derotate_90(float a);
double derotate_90_double(double a);

//START REGION DEVELOP

/* Derotates an angle between -180 and 180 degrees */
double derotate_deg_small_double(double a);

//START REGION RELEASE
/* Derotates an angle between -90 and 90 degrees */
double derotate_180_small_double(double a);


// Derotates an angle between -45 and 45 degrees
double derotate_90_small_double(double a);

//START REGION DEVELOP

//START REGION RELEASE

/* Calculates the angle of the vector (x,y) in radians. (1,0)
   corresponds with 0 rad and counter clockwise corresponds to an
   increasing angle. */
float polar_angle_rad(float x, float y);
//START REGION DEVELOP
double polar_angle_rad_d(double x, double y);

/*
  Calculates the azimuth and zenith angle for a source at a given lst
  and latitude. All arguments are angles in radians.
 */
void calculate_az_zen(double raj, double decj, double lst, double latitude, double *az, double *zen);

/* For a given right ascention and declination calculate the rise and
   set times (LST) for a given elevation limit. All parameters are
   angles in radians).

   Return values:
     0 - Never visible
     1 - Does rise and set
     2 - Always visible
 */
int riseset(double alpha, double delta, double elevation_limit, double latitude, double *rise, double *set, int verbose);


/* Prints a angle in radians as a time */
void printhour(double t);

/* Get the position of the vector from the pulsar to Earth with
   respect to the magnetic axis (in polar coordinates, radians) at
   time t for a pulse period period (should have the same
   units). alpha is the angle between the magnetic axis and the
   rototation axis and zeta between the rotation axis and the line of
   sight. All angles are in radians. */

void pulsar_polar_coordinates(float t, float period, float alpha, float zeta, float *phi, float *theta);
void pulsar_polar_coordinates_d(double t, double period, double alpha, double zeta, double *phi, double *theta);

/* Converts the (vector dr, dtheta, dphi) at position (theta, phi)
   into the Cartesian vector (dx, dy, dz) */
void dSpherical2dCartesian(double theta, double phi, double dr, double dtheta, double dphi, double *dx, double *dy, double *dz);

/* Converts the (vector r, theta, phi) into the Cartesian vector (x,
   y, z). The z-axis is the polar axis (theta=0), the x-axis to phi=0
   and the y-axis to phi=pi/2. All angles are in radians. */
void spherical2Cartesian(double r, double theta, double phi, double *x, double *y, double *z);
/* Reverse of spherical2Cartesian. */
void cartesian2spherical(double *r, double *theta, double *phi, double x, double y, double z);

/* If rmin or rmax are negative it will calculate the whole field line
   instead of section. theta0 (in deg) is the footprint polar angle of
   the field line from the magnetic axis. alpha (in deg) is the angle
   between the rotation axis (positive z-axis) and the magnetic
   axis. Set s to +1 or -1 to indicate which side of the polar cap
   you're considering. */
void fieldLine(float theta0, float alpha, int s, float *x, float *z, long nrpoints, float Rns, float rmin, float rmax);
void fieldLine_d(double theta0, double alpha, int s, double *x, double *z, long nrpoints, double Rns, double rmin, double rmax);

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
long fieldLine_Deutsch_d(double theta0, double phi0, double alpha, double *x, double *y, double *z, int *openfieldline, long nrpoints, double Rns, double P, double rmin, double rmax, double eps, int nounderflow_error);

/* Same function, but not using NR. theta0 and phi0 (in deg) are the
   polar (from the magnetic axis) and azimuthal angle of the field
   line at rmin. alpha (in deg) is the angle between the rotation axis
   (positive z-axis) and the magnetic axis. All distances should be in
   meters and the period P in seconds. */
void fieldLine_Deutsch_d_old(double theta0, double phi0, double alpha, double *x, double *y, double *z, long nrpoints, double Rns, double P, double rmin, double rmax);


/* Formula's from: Dyks & Rudak 2003, ApJ 598 1201: Two pole caustic model 

Takes as input: position (x,y,z) and normalized direction (dx, dy, dz)
and the light cylinder radius rlc. Distances should be in the same
units. Returns pulse longitude and angle between rotation axis and los
(in degrees). dt_ret is the time travel delay in degrees.
*/
void timedelays(float x, float y, float z, float dx, float dy, float dz, float rlc, float *phi, float *zeta, float *dt_ret, int noaberration, int noretardation);
void timedelays_d(double x, double y, double z, double dx, double dy, double dz, double rlc, double *phi, double *zeta, double *dt_ret, int noaberration, int noretardation);

/* Lorentz transform (normalized) propagation direction (dx,dy,dz) on
   position (x,y,z) in the corotating frame to the observer frame. rlc
   is the light cylinder distance. All distances should be in the same
   units. */
void aberration(float x, float y, float z, float dx, float dy, float dz, float *dx2, float *dy2, float *dz2, float rlc);
void aberration_d(double x, double y, double z, double dx, double dy, double dz, double *dx2, double *dy2, double *dz2, double rlc);

/* Z is up, X is right, Y is inside. Rotates Z towards Y. Angle is in degrees.  */
void rotateX(long nrpoints, float angle, float *x, float *y, float *z);
/* Z is up, X is right, Y is inside. Rotates Z towards X. Angle is in degrees.  */
void rotateY(long nrpoints, float angle, float *x, float *y, float *z);
/* Z is up, X is right, Y is inside. Rotates X towards Y. Angle is in degrees.  */
void rotateZ(long nrpoints, float angle, float *x, float *y, float *z);
void rotateX_d(long nrpoints, double angle, double *x, double *y, double *z);
void rotateY_d(long nrpoints, double angle, double *x, double *y, double *z);
void rotateZ_d(long nrpoints, double angle, double *x, double *y, double *z);

//START REGION RELEASE

/* All angles are in degrees. */
float paswing(float alpha, float beta, float l, float pa0, float l0, int nrJumps, float *jump_longitude, float *jump_offset, float add_height_longitude, float add_height_shift);
double paswing_double(double alpha, double beta, double l, double pa0, double l0, int nrJumps, double *jump_longitude, double *jump_offset, double add_height_longitude, double add_height_shift);

//START REGION DEVELOP
//START REGION RELEASE

/* ***********************************
   Functions in pgplot.c 
   *********************************** */

void print_pgplot_version_used(FILE *stream);

/* Sets viewport parameters to default values. */
void pgplot_clear_viewport_def(pgplot_viewport_definition *viewport);

/* Sets box options of pgplot plot to default values. */
void clear_pgplot_box(pgplot_box_definition *box);

void pgplot_clear_options(pgplot_options_definition *pgplot);

// Return 1 = success, 0 = fail
int initPulselongitudeRegion(pulselongitude_regions_definition *region, verbose_definition verbose);
// This function is automatically called from within initPulselongitudeRegion
void clearPulselongitudeRegion(pulselongitude_regions_definition *region);
// Release the memory associated with an onpulse region
void freePulselongitudeRegion(pulselongitude_regions_definition *region);
// Copy a region definition to an already initialised other onpulse region
void copyPulselongitudeRegion(pulselongitude_regions_definition source, pulselongitude_regions_definition *destination);

/* windowwidth & windowheight = (if specified) sets the resolution of the device. 
otherwise aspectratio (if set) sets the aspect ratio of the device, without specifying the resolution. */
void pgplot_setWindowsize(int windowwidth, int windowheight, float aspectratio);

/* Clear frame parameters */
void clear_pgplot_frame(pgplot_frame_def_internal *frame);

//START REGION DEVELOP
//START REGION RELEASE

/* Draw a box and labels.  */
void pgplot_drawbox(pgplot_box_definition *box);

/* If you use regions, make sure to clear it first */
void clearRegion(pulselongitude_regions_definition *region);

//START REGION DEVELOP

/* Given the array selection of nrbins integers, construct regions
   based on those bins where selection != 0. Returns 1 on success, 0
   on error. */
int selection_to_region(int *selection, int nrbins, pulselongitude_regions_definition *onpulse, verbose_definition verbose);

/* Prints the regions defined in regions. Se the nr of spaces before
   the output by indent. */
void printRegions(pulselongitude_regions_definition *regions, verbose_definition indent);

//START REGION RELEASE

/* Convert all regions defined as fracions to bin numbers via: 
   int = frac*scale + offset
*/
void region_frac_to_int(pulselongitude_regions_definition *region, float scale, float offset);

/* Convert all regions defined as bin numbers to pulse phase via: 
   frac = bin*scale + offset
*/
void region_int_to_frac(pulselongitude_regions_definition *region, float scale, float offset);

//START REGION DEVELOP
//START REGION RELEASE

/* Show the selected regions as command line options (option). Where
   could be set to stdout for instance. If optionFrac != NULL and
   region contains a definition in terms of pulse phase, that will be
   outputted as well. */
void regionShowNextTimeUse(pulselongitude_regions_definition region, char *option, char *optionFrac, FILE *where);

//START REGION DEVELOP
//START REGION RELEASE

/* Opens a pgplot device (deviceID is set) and sets the resolution if
   provided. Returns 0 if not successful. */
int pgplot_opendevice(pgplot_viewport_definition *viewport, int *deviceID, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/*
  Plot 'data' with nrx points. Properties of the device is set with
  viewport and properties of the box and labels with pgplotbox. If
  datax is set, those x values will be used. If datax is set to NULL,
  the data should be uniformely sampled between xmin and xmax, so data
  is an one dimensional array. If sigma != NULL, it will be used to
  plot errorbars. If ymin_show == ymax_show the minumum and maximum
  are obtained from the data (if forceMinZero is set than the yrange
  is set from 0 to the maximum in the data). If ymin_show !=
  ymax_show, these will used for the vertical range. Set xmin_show and
  xmax_show to zoom (if they are equal the full range is used
  instead). If dontsetranges is set, the ranges are expected to be
  already by previous pgplot calls. If hist is set, a histogram is
  produced rather than connected points (if forceMinZero is set, it
  means the outer edges of the histogram will go to zero.). If noline
  is set there is no line drawn in between the points. The width of
  the curve can be set with linewidth. If pointtype is set the points
  are drawn with this pgplot symbol and color sets the color. boxcolor
  sets the color of the axis (default = 1). You can specify regions by
  using the regions struct, or you can set it to NULL if not
  used. This will increase the color index from color (offpulse) by 1
  for each subsequent selected onpulse region. Alternatively, if
  onpulsecolor is non-negative, the color is onpulsecolor inside the
  selected region(s) and color outside the regions.

  Return 1: Succesful
         0: Error
*/

int pgplotGraph1(pgplot_options_definition *pgplot, float *data, float *datax, float *sigma, int nrx, float xmin, float xmax, int dontsetranges, float xmin_show, float xmax_show, float ymin_show, float ymax_show, int forceMinZero, int hist, int noline, int linewidth, int pointtype, int color, int boxcolor, pulselongitude_regions_definition *regions, int onpulsecolor, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/*
  Plot 2D data (cmap) with nrx by nry points. viewport sets the size
  of the device and the placement of the graph within the device. The
  data should be uniformely sampled between xmin and xmax and ymin and
  ymax. You can zoom in by choosing xmingshow, xmaxshow
  differently. You can plot contours as well. In that case you can
  either provide nrcontours levels via contours, or set contours=NULL,
  in which case the contours are equally separated across the range of
  the data. When plotting contours, you might want to use the nogray
  option to only show the contours. The image transfer function can be
  set with itf (0=linear, 1=logarithmic, 2=square-root). You can use
  levelInversion to invert the colors (usefull in combination with
  itf). If forceMinZero is set than the black corresponds to 0, or
  else black is set to the minimum in the data. White is set to the
  maximum of the data if saturize is 1. Larger values of saturize mean
  that the peaks become clipped. Instead, if levelset is set, you can
  set the levels by levelmin and levelmax. The look of the labels etc
  is set with pgplotbox. If onlyData is set no labels/boxes are drawn
  (if onlyData is set to 2, the side panel boxes are drawn and the
  title if set). If sideright is set a collapsed panel is shown at the
  right. If set to 2, the numbers are not plotted. If set to 3 or 4,
  it is the same as 1 and 2, except that no tickmarkes are placed on
  axis along the map. If forceMinZeroRight is set the minimum of the
  scale in the side panel is set to zero. sidelw sets the linewidth
  used to draw the curve in the side panels. Showwedge indicates if
  there should be a color-range bar in the plot. If plotSubset is set,
  only the relevant bit of the map will be sent to pgplot, thereby
  making ps files potentially a lot smaller. If debug is set info
  about the color scale is printed to the stdout. The xmin,xmax and
  ymin,ymax corresponds to the centres of the relevant bins. If
  showTwice is set the map is plotted twice above each other.

  The color maps are:
    PPGPLOT_GRAYSCALE
    PPGPLOT_INVERTED_GRAYSCALE
    etc.

  Return 1: Succesful
         0: Plot device couldn't be opened or no memory could be allocated
*/
int pgplotMap(pgplot_options_definition *pgplot, float *cmap, int nrx, int nry, float xmin, float xmax, float xminshow, float xmaxshow, float ymin, float ymax, float yminshow, float ymaxshow, int maptype, int itf, int nogray, int nrcontours, float *contours, int contourlw, int forceMinZero, float saturize, int levelset, float levelmin, float levelmax, int levelInversion, int onlyData, int sideright, int forceMinZeroRight, int sidetop, int forceMinZeroTop, int sidelw, int showwedge, int plotSubset, int showTwice, verbose_definition verbose);

/* Function sets nx and ny to the array ellements in the data that has
   been plotted with pgplotMap at coordinates (x, y). If a point is
   selected outside the array the point is set at the border. In that
   case the return value is 1 instead of 0. */
int pgplotMapCoordinate(float x, float y, int *nx, int *ny);
int pgplotMapCoordinate_dbl(double x, double y, int *nx, int *ny);

//START REGION RELEASE

/* Inverse of pgplotMapCoordinate. Returns coordinates of the centre of bin. */
void pgplotMapCoordinateInverse(float *x, float *y, int nx, int ny);
void pgplotMapCoordinateInverse_dbl(double *x, double *y, int nx, int ny);

//START REGION DEVELOP
//START REGION RELEASE

/* Returns the size of one bin in the map plotted with pgplotMap. */
void pgplotMapCoordinateBinSize(float *dx, float *dy);

//START REGION DEVELOP
//START REGION RELEASE

/* Ask user to define pulse longitude ranges in a graph with nrBins
   bins. If onlyOne is specified the user can keep on selecting one
   component until he/she is satisfied. If powerTwo is selected the
   region is forced to be a power of two. If evenNumber is set, the
   region is forced to be an even number. In viewport, if dontclose is
   set the plotdevice is not closed, so you can do additional
   drawing. If dontopen is set it assumes a plotting device is opened
   and selected. It returns the number of selected regions, which can
   be zero if nothing is selected. -1 is returned if device doesn't
   support a cursor. */
int selectRegions(float *profileI, int nrBins, pgplot_options_definition *pgplot, int onlyOne, int powerTwo, int evenNumber, pulselongitude_regions_definition *regions, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Check if bin is in the region specified by whichregion (counting
   from 1). Zero is returned when bin is outside the specified
   whichregion. If whichregion is 0, then it will return the first
   region (counting from 1) in which bin falls (or 0 if it is in none
   of the regions). If regions set to NULL it will return 0.*/
int checkRegions(int bin, pulselongitude_regions_definition *regions, int whichregion, verbose_definition verbose);

//START REGION DEVELOP

/* Counts the uniqe bins specified by regions. It does not double
   count bins which appear in multiple selections. */
int count_nrbins_in_Regions(pulselongitude_regions_definition *regions, int debug);

//START REGION RELEASE
/*
  Strips the device type from given device name, i.e. it returns:

  0  = undefined device type
  1  = /xwindow
  1  = /xw
  2  = /xserve
  2  = /xs
  3  = /ps
  4  = /vps
  5  = /cps
  6  = /vcps
  7  = /latex
  8  = /null
  9  = /png
  10 = /tpng
  100 = /pw

  If devicename = NULL, the name of the current opened device is used
  instead.

 */
int pgplot_device_type(char *devicename, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE


/*
Make a PA plot
data = data file
showtotpol = nonzero: plot the total degree of polarization
nopaswing = nonzero: do not produce PA-swing panel
showEll = nonzero: plot the ellipticity angle
viewport = The name of the device etc
xlabel, ylabel, ylabel_pa = labels (set title in pgplotbox (keywords are substituted with header parameters with the str_replace_header_params() function), xlabel and ylabel are ignored in pgplotbox)
longitude_left & longitude_right = pulse longitude range in degrees
xunit_type: 0=degrees, 1=phase
Imin, Imax = yrange of intensities. Set both to zero for autoscale
pa_bottom & pa_top = PA range
PAoffset = extra vertical offset PA's in plot
sigma_limit = required S/N to plot PA's
datalinewidth = line width of data
ysize2 = relative ysize of PA plot [1=default]
dashed = dashed lines instead (usefull for B/W plots)
noynumbers = Don't plot numbers on vertical axis
textoption = the option in command line for the texts, i.e. -text
herrorbaroption = the option in command line for the horizontal errorbars, i.e. -herrorbar
herrorbaroption2 = the option in command line for the horizontal errorbars in the pa plot, i.e. -herrorbar2
verrorbaroption = the option in command line for the vertical errorbars, i.e. -verrorbar
verrorbaroption2 = the option in command line for the vertical errorbars in the pa plot, i.e. -verrorbar2
argc, argv = number of command line options and the options themselves. 
             It parses for command lines like: 
             -text "x y ch lw font color" "text" 
outline_txt = The text appears as an outline
outline_lw, outline_color = line width and color of the outline
overlayPA = Overplot an RVM model
overlayalpha, overlaybeta, overlaypa0, overlayl0 = parameters of RVM curve
overlayPAfine: If set, a much finer resolution than the resolution of the data is used to draw the RVM curve
nrJumps: number of OPM jumps to be included in the RVM curve
jump_longitudes & jump_offsets: the list of longitudes and offsets where the OPMs occur
dontclose = plotdevice is not closed
dontopen = it assumes a plotting device is opened and selected
padist = the pa distribution to be attached underneath the plot, if not NULL pointer and if the data is already read in successfully
padist_pamin, padist_pamax = range in PA covered by PA distribution in degrees
padist_saturize: If > 1, the count rate will get saturated
elldist = same for ellipticity angle distribution
elldist_saturize: If > 1, the count rate will get saturated
nowedge: if non-zero, no wedge is shown next to PA/Ellipticity distributions
  Return 1: Succesful
         0: Plot device couldn't be opened or data couldn't be plotted
*/

int pgplotPAplot(datafile_definition data, int showtotpol, int nopaswing, int showEll, pgplot_options_definition *pgplot, char *xlabel, char *ylabel, char *ylabel_pa, char *ylabel_ell, float longitude_left, float longitude_right, int xunit_type, float Imin, float Imax, float pa_bottom, float pa_top, float PAoffset, float sigma_limit, float datalinewidth, float ysize2, int dashed, int noynumbers, char *textoption, char *herrorbaroption, char *herrorbaroption2, char *verrorbaroption, char *verrorbaroption2, int argc, char **argv, int outline_txt, int outline_lw, int outline_color, int overlayPA, float overlayalpha, float overlaybeta, float overlaypa0, float overlayl0, int overlayPAfine, int nrJumps, float *jump_longitudes, float *jump_offsets, datafile_definition *padist, float padist_pamin, float padist_pamax, float padist_saturize, datafile_definition *elldist, float elldist_saturize, int nowedge, verbose_definition verbose);

//START REGION DEVELOP

/* 
   The angles dlat and dlong set the step size in latitude and
   longitude in degrees. lw sets the line width. Normally longitude 0
   and latitude 0 defines the centre of the projection, but this grid
   can be rotated by using rot_long and rot_lat, which are added to
   the longitude and latitude lines which are drawn.

   projection = 1 Draws a grid of a Hammer-Aitoff projection. 
   projection = 2 Draws a grid of a 3D sphere. 
   projection = 3 Draws a grid of a simple linear longitude-latitude map

   For projection 1 the poles are assumed to be at (0, +1) and (0, -1)
   and the equator extends from (-2, 0) to (+2, 0). For projection 2
   the sphere has a radius of 1. If a point lies on the front-side,
   the sphere is centred at x=-1.125, y=0. If the point lies on the
   far side the point ends on a sphere with an x-offset in the
   opposite direction. For projection 3 the horizontal axis runs from
   -180 to 180 degrees and the vertical axis from -90 to 90 degrees.
   */
void drawSphericalGrid(float dlat, float dlong, float rot_long, float rot_lat, int lw, int projection);

/* Print all supported pgplot commands */
void printPgplotInterpratorHelp();
void printSubjectHelp(char *substr);

/* Executes pgplotinterprator command.
Return values:
 0   = error
 1   = processed ok
 100 = Reached quit command
*/
int doPgplotInterprator(char *txt, int debug);

/* Parses commandline for -pgplot options. Returns 1 if successful. */
int commandline_PgplotInterprator(int argc, char **argv, int debug);

//START REGION RELEASE

/* Show the -cmap options */
void printCMAPCommandlineOptions(FILE *printdevice);

/* Returns 0 if parse error, or else it returns the map type.  */
int cmap_parse_commandline(int argc, char **argv, int debug);

//START REGION DEVELOP

/* Returns 1 if succesfull. */
int pgopenoutputfile(char *filename);
void pgcloseoutputfile();
//START REGION DEVELOP
//START REGION RELEASE
int ppgopen(const char *device);
void ppgclos(void);
void ppgend(void);
void ppgsvp(float xleft, float xright, float ybot, float ytop);
void ppgswin(float x1, float x2, float y1, float y2);
void ppgsci(int ci);
void ppgbbuf(void);
void ppgebuf(void);
void ppgmove(float x, float y);
void ppgdraw(float x, float y);
void ppgpt1(float xpt, float ypt, int symbol);
void ppgerr1(int dir, float x, float y, float e, float t);
void ppgask(int flag);
void ppgpage(void);
void ppgslw(int lw);
void ppgsls(int ls);
void ppgsfs(int fs);
void ppgsch(float size);
void ppgqcs(int units, float *xch, float *ych);
void ppglab(const char *xlbl, const char *ylbl, const char *toplbl);
void ppgbox(const char *xopt, float xtick, int nxsub, const char *yopt, float ytick, int nysub);
void ppgaxis(const char *opt, float x1, float y1, float x2, float y2, float v1, float v2, float step, int nsub, float dmajl, float dmajr, float fmin, float disp, float orient);
void ppgmtxt(const char *side, float disp, float coord, float fjust, const char *text);
void ppgtick(float x1, float y1, float x2, float y2, float v, float tikl, float tikr, float disp, float orient, const char *str);
void ppgscr(int ci, float cr, float cg, float cb);
void ppgsitf(int itf);
void ppgarro(float x1, float y1, float x2, float y2);
void ppgcirc(float xcent, float ycent, float radius);
void ppgtext(float x, float y, const char *text);
void ppgptxt(float x, float y, float angle, float fjust, const char *text);
void ppgscf(int font);
void ppgshs(float angle, float sepn, float phase);
void ppgrect(float x1, float x2, float y1, float y2);
void ppgerry(int n, const float *x, const float *y1, const float *y2, float t);
void ppggray(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float fg, float bg, const float *tr);
void ppgimag(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float a1, float a2, const float *tr);
void ppgcont(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, const float *c, int nc, const float *tr);
void ppgpoly(int n, const float *xpts, const float *ypts);
void ppgconl(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, float c, const float *tr, const char *label, int intval, int minint);
int ppgband(int mode, int posn, float xref, float yref, float *x, float *y, char *ch_scalar);
int ppgcurs(float *x, float *y, char *ch);
void ppgscir(int icilo, int icihi);
void ppgctab(float *l, float *r, float *g, float *b, int nc, float contra, float bright);
void ppgwedg(const char *side, float disp, float width, float fg, float bg, const char *label);
int ppgqid(int *id);
int ppgslct(int id);
void ppgpap(float width, float aspect);
/* Returns 1 if successful */
int  ppgqvp(int units, float *xleft, float *xright, float *ybot, float *ytop);
/* Returns 1 if successful */
int  ppgqwin(float *xleft, float *xright, float *ybot, float *ytop);
// Note: value_length must be initialised with length of string value (ie. the available space)
void ppgqinf(const char *item, char *value, int *value_length);
//START REGION DEVELOP


/* ***********************************
   Functions in psrio_puma.c
   *********************************** */

/* Return the number of polarization channels in header */
//int puma_nrpol(Header_type hdr);

/* Read gated puma file and produce a pulse stack. The file pointer
   should be set to the start of the data and the profile is written
   to profileI. Set zapMask to NULL is you don't want to zap pulses. */
//void puma_read_profile(FILE *fin, long nrPulses, int nrBins, int nrFreqChan, float *profileI, int *zapMask, int verbose);

/* Read gated puma file and produce a list of rms'es and averages. The
   file pointer should be set to the start of the data and the profile
   is written to profileI. Set zapMask to NULL is you don't want to
   zap pulses. The defined regions are NOT included, unless invert is
   set. */
//void puma_read_rms(FILE *fin, long nrPulses, int nrBins, int nrFreqChan, int nrPol, float *rms, float *avrg, int *zapMask, pulselongitude_regions_definition *regions, int invert, int verbose);

/* Read gated puma file and produce a list of correlations and
   averages. The file pointer should be set to the start of the data
   and the profile is written to profileI. Set zapMask to NULL is you
   don't want to zap pulses. The defined regions are NOT included,
   unless invert is set. */
//void puma_read_correl(FILE *fin, long nrPulses, int nrBins, int nrFreqChan, int nrPol, float *correl, float *avrg, int *zapMask, pulselongitude_regions_definition *regions, int invert, int verbose);


/* Read gated puma file and produce a list of peak rms'es and averages
   (after adding addbins together). The file pointer should be set to
   the start of the data and the profile is written to profileI. Set
   zapMask to NULL is you don't want to zap pulses. The defined
   regions are NOT included, unless invert is set.

   If return value = 0, then there is a memory allocation error.
*/
//int puma_read_peakrms(FILE *fin, long nrPulses, int nrBins, int nrFreqChan, int nrPol, float *rms, float *avrg, int *zapMask, pulselongitude_regions_definition *regions, int invert, int addbins, int verbose);

/* ***********************************
   Functions in vonMises.c 
   *********************************** */

//START REGION RELEASE
/* Read a text file with the vonMises functions. 

  Return 1: OK
         0: Cannot open file
*/
int readVonMisesModel(char *filename, vonMises_collection_definition *components, verbose_definition verbose);
//START REGION DEVELOP

/* Write a text file with the vonMises functions.

   If filename = NULL, the stdout is used as input instead.

  Return 1: OK
         0: Cannot open/write file
*/
int writeVonMisesModel(char *filename, vonMises_collection_definition *components, verbose_definition verbose);

//START REGION RELEASE
/* Calculate the profile out of the model by at a certain phase by
   applying a phase shift. */
double calcVonMisesFunction(vonMises_collection_definition *components, double phase, double shift);

/* Calculate a vonMises function at a certain phase by applying a
   phase shift. */
double calcVonMisesFunction2(double centre, double concentration, double height, double phase, double shift);

//START REGION DEVELOP
/* Integrate a vonMises function over its entire domain. 
Here the 2pi in the formula corresponds to a full rotation in radians, so the units are intensity*radians. 
Replace the 2pi with 360 or the number of bins across the profile to get it in different units.
*/
double integrateVonMisesFunction2(double concentration, double height);

/* Returns the width of a single vonMises function given its
   concentration in radians. The width is measured at a level ampfrac,
   which is a number between 0 and 1. A fwhm corresponds to ampfrac =
   0.5 and a W10 to an ampfrac = 0.1. For example the fwhm is defined
   at half the peak amplitude. So the actual minimum of the function
   (which is only asymptotically approaching zero for narrow
   distributions) is not considered. Returns sqrt(-1) if the fwhm doesn't
   exists. This happens if the concentration is very low, which means
   that the function never goes to zero. */
double widthVonMisesFunction2(double concentration, double ampfrac);

/* Calculate the integral of the model over its entire domain.  Here
the 2pi in the formula corresponds to a full rotation in radians, so
the units are intensity*radians.  Replace the 2pi with 360 or the
number of bins across the profile to get it in different units.*/
double integrateVonMisesFunction(vonMises_collection_definition *components);

//START REGION RELEASE

/* Calculate the profile out of the model by applying a phase
   shift. If the normalize flag is set the profile will be
   normalized. */
void calcVonMisesProfile(vonMises_collection_definition *components, int nrbins, float *profile, double shift, int normalize);

//START REGION DEVELOP

/* For a given profile with nrbins bins, calculate the rms of the
   profile - model after applying a phase shift to the model. If the
   normalize flag is set the profile will be normalized. */
void calcVonMisesProfile_resid_rms(vonMises_collection_definition *components, int nrbins, float *profile, double shift, double *rms);

/* Calculate a shape parameter based on the profiled model after
   applying a phase shift. phaseprecision, a number << 1, indicates
   how precise the shape parameter needs to be. More precision = more
   computations. The possibilities for shapepar are defined in
   psrsalsa_defines.h, and start with SHAPEPAR_. The shape parameter
   measure might depend on one or mor additional input values (a range
   for example), which are provided as shapepar_aux. The measurement
   is returned as measurement.
*/
void calcVonMisesProfile_shape_parameter(vonMises_collection_definition *components, double shift, double phaseprecision, int shapepar, double *shapepar_aux, double *measurement, verbose_definition verbose);

/* Print the shape parameter, plus a description if showdescr != 0. The error is quoted when >= 0. */
void print_shape_par(FILE *fout, int showdescr, int shapepar, double measurement, double error);

//START REGION RELEASE
/* Find the best shift in phase to match the profile (the return
   value). verbose-1 is number of spaces before output. */
float correlateVonMisesFunction(vonMises_collection_definition *components, int nrbins, float *profile, verbose_definition verbose);
//START REGION DEVELOP


void find_boundaries(float *profile, int nrbins, float y, pulselongitude_regions_definition *regions);

/* Fit a vonMises function through the pulse profile in
   psrdata. MaxNrComp defines the maximum nr of components which are
   allowed to be fitted and sigmalimit the limit up to where
   components are fitted (5 seems to work reasonably good). If
   educatedguess it will use the boxcarFindpeak algorithm to find
   initial guesses (this option seems to work best, so you should use
   it by default). If avoid_neg_components is non-zero, then
   components with negative amplitudes are avoided. If fitbaseline is
   set it is not assumed that the baseline is removed and its level is
   stored in baseline. The larger precision is set to the finer trials
   are tried. Something like 2 seems reasonable. If onpulse != NULL (and the number of regions is non-zero),
   the components must be centered within the specified range.
   verbose-1 is number of spaces before output.

   superverbose = 1: extra print statements
                  2: plot best intermediate steps
                  3: plot all intermediate steps

Return 1 is success 
*/
int fitvonmises(datafile_definition psrdata, vonMises_collection_definition *components, int MaxNrComp, float sigmalimit, int educatedguess, int fitbaseline, int precision, int avoid_neg_components, float *baseline, pulselongitude_regions_definition *onpulse, verbose_definition verbose, int superverbose, char *plotdevice);

/* Fit a given vonMises function through the pulse profile in psrdata
   and refine its parameters. If fitbaseline is set it is not assumed
   that the baseline is removed and its level is stored in
   baseline. If avoid_neg_components is non-zero, then components with
   negative amplitudes are avoided.  If fixamp, fixwidth and fixphase
   are set, the amplitudes, widths and/or phases are held fixed. If
   fixrelphase is set, the components cannot shift relative to each
   other. With fixrelamp, the components cannot change height relative
   to each other. If there is only one fit parameter (for instance
   when fitting an overall scaling, without a baseline variation), a
   dummy fit parameter is introduced to avoid the downhill-simplex
   method to fail. Not elegant, but it works.

   return 0=error,    1=ok */
int fitvonmises_refine_model(datafile_definition psrdata, vonMises_collection_definition *components, int fitbaseline, int avoid_neg_components, float *baseline, int fixamp, int fixwidth, int fixphase, int fixrelamp, int fixrelphase, verbose_definition verbose);

/* ***********************************
   Functions in statistics.c 
   *********************************** */

//START REGION RELEASE
/* Add zeropad zero's before and after data and then add extra zero's
   to make the data length a power of two. If circularpad is set the
   padding is not done with zero's but with the edges of the data
   itself. When duplicate is set the profile is first duplicated,
   which can be usefull when correlating profiles which occur at the
   edge. Then it correlates the data and returns the lag number where
   the maximum was found. correl_max indices how many times the
   maximum of the correlation function was higher than the
   minimum. Returns 1 if succesful. */
int find_peak_correlation(float *data1, float *data2, int ndata, int zeropad, int circularpad, int duplicate, int *lag, float *correl_max, verbose_definition verbose);
//START REGION DEVELOP

//START REGION RELEASE
/* Generates a random idnum to be used in Numerical Recipies random
   number generator. There are 60,000,000 unique values that can be
   genetated. */
void randomize_idnum(long *idnum);
//START REGION DEVELOP

/* Generates a random unsigned int. There are 60,000,000 unique values
   that can be generated. Can be used for instance to generate unique
   temporary file names */
long randomUnsignedInt();

/*
  Converts a number from a uniform distribution (between 0 and 1) to a
  number from a sinusoidal distribution (between 0 and 180, peaking at
  90).
 */
double uniform_to_sinusoidal_180(double uniform);

/*
  Converts a number from a uniform distribution (between 0 and 1) to a
  number from a sinusoidal distribution (between 0 and 90, peaking at
  90).
 */
double uniform_to_sinusoidal_90(double uniform);

/*
  Interpolate the functions to get dataY(x). The value y is found via
  linear interpolation. Order indicates what type of interpolation is
  done. 2 indicates a 2-point interpolation (linear
  interpolation). The data points are assumed to be stored in
  increasing order.

  returns 1 on success, 0 if out of bounds.

  interpolate() will call interpolate_double() 
 */
int interpolate_double(double *dataX, double *dataY, long nrdatapoints, double x, double *y, long int order, verbose_definition verbose);
int interpolate(float *dataX, float *dataY, long nrdatapoints, float x, float *y, long int order, verbose_definition verbose);


/*
  Fits cubic B-spline functions to get dataY(x), and dataYerr defines
  the errorbar on the data points (set to NULL if not used). The value
  y is found fitting cubic B-spline basis functions with uniform
  breakpoints. The number of fit coefficients is set by
  nrcoefficients. Set y=NULL and/or yerr=NULL if you're not interested
  in the error on the y estimate. The fit is returned in splineFit
  (unless it is set to NULL, memory will be allocated) for future use,
  such as with cubicBspline_eval(). Use cubicBspline_free() to free up
  this memory. dataX is not required to be sorted. If nrtrials is set,
  instead of assuming uniform gridding of the nodes (which is first
  trial), random distributions of node positions are explored
  (nrtrials nr of times) to find a better fit. The variable
  frac_replace defines the fraction of grid-points which are replaced
  by random values in each itteration. verbose-1 is number of spaces
  before output. If debug is set, more information is printed.

  returns 1 on success, 0 on error.

  cubicBspline_fit() will call cubicBspline_fit_double() 
*/
int cubicBspline_fit_double(double *dataX, double *dataY, double *dataYerr, long nrdatapoints, int nrcoefficients, double x, double *y, double *yerr, splineFit_def *splineFit, long nrtrials, double frac_replace, verbose_definition verbose);
int cubicBspline_fit(float *dataX, float *dataY, float *dataYerr, long nrdatapoints, int nrcoefficients, float x, float *y, float *yerr, splineFit_def *splineFit, long nrtrials, float frac_replace, verbose_definition verbose);

// Free up allocated space after spline was allocated with cubicBspline_fit_double
void cubicBspline_free(splineFit_def *splineFit);

/*
  Using a fit made by cubicBspline_fit(), calculate the value of the
  cubic B-spline functions at x.

  returns 1 on success, 0 if out of bounds.

  cubicBspline_eval() will call cubicBspline_eval_double() 
*/
int cubicBspline_eval_double(splineFit_def splineFit, double x, double *y, verbose_definition verbose);
int cubicBspline_eval(splineFit_def splineFit, double x, float *y, verbose_definition verbose);

/* This wrapper function will call smoothFitFunction_double. */
int smoothFitFunction(float *dataX, float *dataY, float *dataYerr, long nrdatapoints, int nrcoefficients, float x, float *y, float *yerr, int verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Calculates in which bin number the value x belongs, given the bin
   width dx and the minimum value of the distribution min_x. Set
   centered_at_zero to one to make the centre of one bin equal to
   zero. exta_phase will shift the bin_centre by a certain fraction of
   a bin. */
long calculate_bin_number(double x, double dx, double min_x, int centered_at_zero, double extra_phase);

// Like calculate_bin_number(), but returns to central value of bin binnr
double calculate_bin_location(long binnr, double dx, double min_x, int centered_at_zero, double extra_phase);

// Like calculate_bin_number(), but returns to bin width dx required to make x fall in the centre of bin binnr. Can be useful to determine required 
double calculate_required_bin_width(double x, long binnr, double min_x, int centered_at_zero, double extra_phase, verbose_definition verbose);

/*
  Given data in the range [min_x_data, max_x_data], determine the
  required binning. If rangex_set is non-zero, the range is adjusted
  to [rangex_min, rangex_max]. nrbins is the requested number of bins,
  which is used when nrbins_specified is set. Otherwise *dx is
  assumed. If centered_at_zero is set, the mid-point of one bin
  (possibly outside range) will coincide with zero. Otherwise an edge
  will coincide with zero. extra_phase changes where zero falls with
  respect to the bin-centre (cyclic behaviour with period of 1).

  The return parameters are the actual range that will be covered
  [min_x, max_x], which might be different from what is specified if
  randex_set is set (to avoid rounding errors). Otherwise, it will be
  set to [min_x_data, max_x_data]. The bin-size is returned in dx.
  
  Return value:
  0: ok
  1: Determination of bin width failed.
  2: Fist bin not in correct location
 */
int set_binning_histogram(double min_x_data, double max_x_data, int rangex_set, double rangex_min, double rangex_max, int nrbins_specified, long nrbins, int centered_at_zero, double extra_phase, double *min_x, double *max_x, double *dx, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/*
  If data2 != NULL, n2 > 0: Given two data-sets, calculate the maximum
  difference between their cummulative distribution, and the
  corresponding probability that the distributions are not sampled
  from the same underlying distribution. A small probability means the
  distributions are statistically shown to be different. The input
  arrays data1 and data2 will be sorted after the call to kstest.

  If data2 == NULL, n2 > 0:
    - If cdf_type == 0: a pointer to the function cdf should be
      provided, which is the cummulative distribution of the distribution
      which is assumed to be the parent distribution from which the points
      data1 are drawn.
    - If cdf_type == 1: Instead use a flat distribution starting/ending 
      at maximum value of data1.
    - If cdf_type == 2: Instead use a sin distribution from 0 to 90 deg.
    - If cdf_type == 3: Instead use a flat distribution starting/ending 
      at values input_value1 and input_value2.

  This implementation should be similar to that implemented in
  Numerical Recipes in C, second edition.

  */
void kstest(double *data1, long n1, double *data2, long n2, int cdf_type, double input_value1, double input_value2, double (*cdf)(double), double *max_diff, double *prob, verbose_definition verbose);

//START REGION DEVELOP

/* ***********************************
   Functions in fitting.c 
   *********************************** */

//START REGION RELEASE
void print_gsl_version_used(FILE *stream);
//START REGION DEVELOP

void show_fitfunctions_commandline_options(FILE *stream, char *cmd_flag, int nrspaces);
int parse_commandline_fitfunctions(int argc, char **argv, char *cmd_flag, fitfunc_collection_type *function, verbose_definition verbose);

//START REGION RELEASE
/* 
This function will show the specified function. If novalue is set,
the functional form is shown, but not the actual values of the
parameters. If showerror is set the error is shown as well as the
value.  A negative index means show all parameters, otherwise the
function index of the function collection is shown.  
*/
int print_fitfunctions(fitfunc_collection_type *function, int novalue, int showerror, int index, verbose_definition verbose);
//START REGION DEVELOP

//START REGION RELEASE
/*
Use the Levenberg-Marquardt algorithm to fit a function specified by
function to ndata datapoints (set data_sigma to NULL to assume equal
errors, which implies force_chi2_1). If oneatatime is speciefied one
parameter is fit at a time in a loop. This makes the fitting slower,
but potentially more robust. If force_chi2_1 is set, the errorbars of
the datapoints are rescaled such that the reduced chi2 is equal to
unity and the errorbars on the fit parameters are probably more
sensible. You can set either epsabs or epsrel (set one to zero),
otherwise the stringents condition of the two will be used. If the
step size |dx_i| is smaller than epsabs or epsrel|x|, then the
function is taken to be converged. maxiter sets the maximum number of
iterations. If showresults is set, or if verbose is set, the fitted
function + chi2 is reported. If showcovariance is set, the covariance
matrix is reported as well. The return value is the same value as
returned via status.

Return values:
  0 = Converged
  1 = Memory allocation error or problem with the initial parameters
  2 = Max number of itterations reached
  3 = Machine precision limit reached
  4 = Cannot determine suitable trial step.
 10 = Shouldn't happen, not documented error code of gsl

algorithm = 1 = GSL
algorithm = 2 = PSRSALSA version
algorithm = 3 = NR double-version hack

 */
int fit_levmar(int algorithm, fitfunc_collection_type *function, double *data_x, double *data_y, double *data_sigma, long ndata, int oneatatime, int force_chi2_1, double epsabs, double epsrel, int maxiter, int *status, int showresults, int showcovariance, verbose_definition verbose);
//START REGION DEVELOP

/* This function is similar to fit_levmar2_ld(), so make changes in both
   if required. 

   This function is designed to look like doAmoeba(), allowing the
   user to switch easily between the two optimisation strategies. Only
   the psrsalsa fitter is being supported at the moment (not the gsl
   or NR equivalent). Switching between algorithms would be
   complicated because the user specified function needs to be quite
   different to be used by the different algorithms.

   This function fits a user-defined function to a data-set. This
   function has nrparams parameters and the initial values are defined
   by the array xstart[0]..xstart[nrparams-1]. The nonzero values of
   fixed indicate if the parameter is fixed and not fitted for. xfit
   is filled with the found solution at which the chi2 has a minimum
   with value called "yfit" to use the same language as the doAmoeba()
   function. ftol (a small number) indicates the precision of the fit
   (in this function it corresponds to the maximum required fractional
   change in the fit parameters before calling the solution converged)
   and nfunc is set to the number of iteration steps it took to
   converge. If finderrors is set it will use the determined
   covariance matrix to find the sigma-level errors (as set by sigma)
   for all nonfixed parameters.  This is stored in dplus and/or dmin
   if != NULL. The values in dplus and dmin will be identical per
   definition, but both can be set to make the function more
   comparable with doAmoeba().

   A function func needs to be provided which is of the form:

   void func(func_params, alpha, beta, &chisq, info);

   So for a set of func_params (nrparams long), is should calculate
   alpha (nxn matrix), beta (n vector) and the chisq. Here n is the
   number of parameters which are fitted for. info is a (void *),
   which can point to anything that might be relevant for that
   function. It is directly passed on from the argument list of
   fit_levmar2_ld().



   Returns 0 on success.
   Returns 1 on maximum number of itterations exceeded error (try smaller ftol)
   Returns 2 on memory error
   Returns 3 if less than 1 fit parameter are not fixed
   Returns 4 if algorithm is not available  NOT RELEVANT
   Returns 5 Other fit fail error (such as singular matrix)
*/
int fit_levmar2_ld(long double *xstart, int *fixed, long double *xfit, long double *yfit, int nrparams, 
		   void (*func)(long double *, long double *, long double *, long double *, void *), void *info, long double ftol, int *nfunk, int finderrors, long double sigma, long double *dplus, long double *dmin, verbose_definition verbose);

//START REGION RELEASE
/* Calculates collection of function at value x. */
double evaluate_fitfunc_collection(fitfunc_collection_type *function, double x, verbose_definition verbose);
//START REGION DEVELOP

// See, for instance, the general least squares discussion in NR
// Input:
// - The solution = the optimum values for the parameters for which the chi2 is minimised (must be done elsewhere). These are nrparams parameters passed to the function in the array optimum_params.
// - dydparam() is a function that takes for argument 1 the parameter number and for argument 2 the optimum_params array. The return value is the derivative of the fit function w.r.t. the parameter number (argument 1) at x=argument number 3
// - The chi2 of the solution, so NOT the reduced chi2
// - If rescale_errors is non-zero, the error-bars will be rescaled such that the REDUCED chi2=1, hence to incorporate a estimation for systematic errors. The reduced chi2 = chi2 / (nrdatapoints - nrfitparameters).
// - Information about the points that are fitted:
//     - Nr of data points = nrdatapoints
//     - The x values
//     - The y values are NOT required
//     - The y-errorbars = sigma
//     - If sigma = NULL, the errorbars are assumed to be equal and will be choosen such that the reduced chi2=1. The chi2 should be calculated by effectively taking all sigma values to be 1. This implies that rescale_errors is effectively always enabled in this case.
// - The covariance matrix is returned as the array covariance_matrix (nrparams*nrparams values). Memory should already be allocated
//
// NOTE: the error's on the parameters (to be procise the variance, so take the the sqrt if you want sigma, the std. dev.) can be found on the diagonal of the covariance matrix.
//
// Return values
//    0 = ok
//    1 = Error (such as memory allocation error)
//    2 = Not enough data points (nrdatapoints <= nrparams)
//    3 = Covariance matrix cannot be computed (matrix cannot be inverted, probably because there are fully covariant parameters?)
int compute_covariance_matrix(long nrparams, double *optimum_params, double (*dydparam)(long, double [], double), double chi2, int rescale_errors, long nrdatapoints, double *x, double *sigma, double *covariance_matrix, verbose_definition verbose_state);


/* ***********************************
   Functions in minimize.c 
   *********************************** */

//START REGION DEVELOP
//START REGION RELEASE

/*
  Given the function: double funk(double x, void *params)

  Minimize the function as function of x (if findroot == 0), otherwise
  find the root of the function. The variable params can be used to pass
  on other (fixed) parameters. The minimization is done using gsl with
  the Brent method.
  
  As input, the search range x_lower and x_upper must be known. The
  result is returned as x_minumum.  epsabs sets the minimum required
  uncertainty in x_minimum, while epsrel set the minimum required
  uncertainty as a fraction of x_minimum. It is allowed to set either
  epsabs or epsrel to zero. max_iter is the maximum nr of iterrations
  (maybe 100).
  
  If gridsearch is > 1, the specified x range in searched over in
  gridsearch points to find the rough location of the minimum first,
  before bracketing the minimum further. This initial gridsearch is
  repeated nested + 1 times. If not set correctly, the minimum can be
  missed. If investigateLocalMinima is set, 3 local minima are
  explored in a higher resolution in order to potentially discriminate
  between them and zoom in on what is hopefully the global
  minimum. This is probably not a good idea if the fit parameter is
  circular (like an angle), as the best solution can be at the lower
  and upper border of the parameter range simultaneously, which
  results in the routine not being able to distinguish between the two
  minima, resulting in an error.
  
  If set, verbose determines the nr of spaces before the output. If
  debug_verbose is set, function values are outputted when doing the
  nested gridsearch.
  
  Return values:
  0 = success, converged
  1 = Maximum nr of itterations reached, did not converge fully
  2 = Did not find root in specified range
  3 = Lower and upper limit do not bracket a root
  4 = memory allocation error
  5 = Other unspecified error
*/
int minimize_1D_double(int findroot, double (*funk)(double, void *), void *params, double x_lower, double x_upper, int gridsearch, int investigateLocalMinima, int nested, double *x_minimum, int max_iter, double epsabs, double epsrel, int verbose, int debug_verbose);

//START REGION DEVELOP
//START REGION RELEASE
/*
  You provide the function: double funk(double *x, void *params),
  which returns the chi2 as function of an array of parameters
  x[nrparameters] (set by xminimum) and potenital additional (fixed)
  parameter params. Of these parameters x[], parameter number paramnr
  (counting from zero) is varied to find an errorbar on that
  parameter, which is defined where the chi2 is (1+sigma) times higher
  than the minimum value chi2min. The chi2 at parameters xminimum
  should be within this limit. Set sigma to a positive value to find
  the errorbar in the positive direction, or to a negative value to
  find the errorbar in the negative direction. The size of the
  errorbar is returned as errorbar.

  The minimization is done using gsl with the Brent method.
  
  As input, xminimum is provided, which is the x value for which funk
  has a minimum (potentially found with function
  minimize_1D_double). Also dx must be provided, which is the stepsize
  used to find the specified sigma point. If separation exceeds dxmax
  (set to negative value to disable this feature), the errorbar is set
  to this value and no further searching is done.

  epsabs sets the minimum required uncertainty in the determined errorbar, while
  epsrel set the minimum required uncertainty as a fraction of
  the value. It is allowed to set either epsabs or epsrel to
  zero. max_iter is the maximum nr of iterrations (maybe 100).
  
  If set, verbose determines the nr of spaces before the output.
  
  Return values:
  0 = success, converged
  1 = Maximum nr of itterations reached, did not converge fully
  2 = Did not find root in specified range
  3 = Lower and upper limit do not bracket a root
  4 = Memory allocation error
  5 = Other unspecified error
*/
int find_1D_error(double (*funk)(double *, void *), double *xminimum, int paramnr, int nrparameters, double dx, double dxmax, void *params, double sigma, double chi2min, int max_itr, double epsabs, double epsrel, double *errorbar, int verbose);
//START REGION DEVELOP

/* ***********************************
   Functions in amoeba.c 
   *********************************** */

/* Does an amoeba search to minimize function funk, which has nrparams
   and is called with the array x[0]..x[nrparams-1]. The start point
   is given by xstart and the inital search step by dx. The nonzero
   values of fixed indicate if the parameter is fixed and not fitted
   for. xfit is filled with the position at which funk has a minimum
   with value yfit. All arrays start at entry 0, not 1. ftol (a small
   number) indicates the precision of the fit and nfunc is set to the
   number of iteration steps it took to converge. If finderrors is set
   it will search the parameterspace for all nonfixed parameters to
   find the sigma errors (assuming funk is a chi^2 surface, dplus &
   dmin are the two errorbars). GSL is not implemented. I cannot
   remember exactly why, but probably because it only supports double
   format.

   algorithm speciefies which code to use
     0: Modified code originally written by Michael F. Hutt
     1: Modified code originally from Numerical Recipes
     2: GSL (unsupported)

   Returns 0 on success.
   Returns 1 on maximum number of itterations exceeded error (try smaller ftol)
   Returns 2 on memory error
   Returns 3 if less than 2 fit parameter are not fixed
   Returns 4 if algorithm is not available
*/
int doAmoeba(int algorithm, float *xstart, float *dx, int *fixed, float *xfit, float *yfit, int nrparams, float (*funk)(float []), float ftol, int *nfunk, int verbose, int finderrors, float sigma, float *dplus, float *dmin);
//START REGION RELEASE
int doAmoeba_d(int algorithm, double *xstart, double *dx, int *fixed, double *xfit, double *yfit, int nrparams, double (*funk)(double []), double ftol, int *nfunk, int verbose, int finderrors, double sigma, double *dplus, double *dmin);
//START REGION DEVELOP
int doAmoeba_ld(int algorithm, long double *xstart, long double *dx, int *fixed, long double *xfit, long double *yfit, int nrparams, long double (*funk)(long double []), long double ftol, int *nfunk, int verbose, int finderrors, long double sigma, long double *dplus, long double *dmin);

/* Copy parameters of doAmoeba, but also give paramnr for which to
   find the errorbar at sigma times the minimum. 
   Returns 0 on success.
   Returns 1 on NMAX exceed error (smaller ftol doesn't work)
   Returns 2 on memory error
   Returns 3 if less than 2 fit parameter are not fixed
*/
int find_errors_amoeba(int algorithm, float *dx, int *fixed, float *xfit, float yfit, int nrparams, float (*funk)(float []), float ftol, int paramnr, float *dplus, float *dmin, float sigma);

//START REGION RELEASE
int find_errors_amoeba_d(int algorithm, double *dx, int *fixed, double *xfit, double yfit, int nrparams, double (*funk)(double []), double ftol, int paramnr, double *dplus, double *dmin, double sigma);
//START REGION DEVELOP

int find_errors_amoeba_ld(int algorithm, long double *dx, int *fixed, long double *xfit, long double yfit, int nrparams, long double (*funk)(long double []), long double ftol, int paramnr, long double *dplus, long double *dmin, long double sigma);


/* ***********************************
   Functions in calc.c 
   *********************************** */

/* Print the mathematics functions to device (could be for instance
   stdio). nrspaces defines the number of spaces before each line.*/
void printCalcFunctions(FILE *printdevice, int nrspaces);

/* Calculates equation and puts it in answer. Set verbose to 1 to see
   the answer and set it to 2 to show all intermediate steps. All
   calculations are done in double format, but the answer can both be
   a float or a double. Returns 0 on error. */
int calc_expression(char *equation, int nrvariables, char variables[][100], double *varables_values, double *answer, int verbose);
int calc_expressionf(char *equation, int nrvariables, char variables[][100], float *varables_values, float *answer, int verbose);

/* ***********************************
   Functions in pulseenergy.c 
   *********************************** */

/*
This function adds the bins between (and including) bin1 and bin2 by
taking into acount the specified baseline. If squared is set the
signal is squared first.
 */
float integratePulseEnergy(float *pulse, int bin1, int bin2, float baseline, int squared);

/*
This function determines the baseline, rms of the offpulse region
OUTSIDE the region defined by onpulse (Set to NULL if not used). If
nodebase is set, the baseline will be assumed to be already
subtracted, which can be useful when a running baseline has been used.
 */
void offpulseStats(float *pulse, int nrBins, float *baseline, float *rms, pulselongitude_regions_definition *onpulse, int nodebase, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Does a boxcar search for the highest s/n peak in pulse. onpulse
   (set to NULL if not used) is used to determine off-pulse statistics
   (baseline and rms). Returns the first bin of the box and its
   pulsewidth, the s/n and the integral over the box (if E_best !=
   NULL). If squared is set the signal is squared. posOrNeg defines if
   the algorithm searches for both positive or negative peaks.

   posOrNeg: 0 = positive peaks only
   posOrNeg: 1 = positive & negative peaks
   posOrNeg: 2 = negative peaks only

   If allwidths is set all bin widths are tried instead of a fixed nr
   of widths (extremely slow). If maxwidth >= 0, it indicates that
   only widths smaller than that number should be tried. If refine is
   set, try all bin-widths around the optimum of the finite nr of
   boxcars search if s2n > 2. The only_onpulse parameter is used to
   limit the search range. If it is set to zero the full pulse phase
   range is searched. If it is set to 1 only the onpulse region is
   searched. If it is set to 2 only the first onpulse region is
   considered, while the other onpulse regions are define the region
   excluded from the off-pulse region. If nodebase is set, the
   baseline will be assumed to be already subtracted, which can be
   useful when a running baseline has been used. Returns 1 if
   successful.  */
int boxcarFindpeak(float *pulse, int nrBins, pulselongitude_regions_definition *onpulse, int *bin, int *pulsewidth, float *snrbest, float *E_best, int squared, int posOrNeg, int allwidths, int refine, int maxwidth, int only_onpulse, int nodebase, verbose_definition verbose);

//START REGION DEVELOP

/*
  Determines the single pulse rms when excluded the onpulse region (if
  not set to NULL). You can select the polarization channel as
  well. If nodebase is set, the baseline will be assumed to be already
  subtracted, which can be useful when a running baseline has been
  used. This is the average determined from analysing each
  channel/subint individually.
 */
float determine_singlepulse_rms(datafile_definition psrdata, int polchan, pulselongitude_regions_definition *onpulse, int nodebase, verbose_definition verbose);


//START REGION RELEASE

/* ***********************************
   Functions in application.c 
   *********************************** */

void initApplication(psrsalsaApplication *application, char *name, char *genusage);
// Releases some memory
void terminateApplication(psrsalsaApplication *application);
void printApplicationHelp(psrsalsaApplication *application);
int processCommandLine(psrsalsaApplication *application, int argc, char **argv, int *index);
void printCitationInfo();

// argv[argv_index] = string to be parsed, so argv[argv_index-1] is the name of the option
// If check_only != 0, no errors messages are generated, so you can "try out" if a given format would work or not
// When parsing a string, the max length (plus termination character) should be provides: e.g. %100s
// minrequestedparameters = minimum number of parameters that should be successfully parsed. This is useful if the user can either provide for example 1 or 2 parameters. If <= 0, all parameters should be read in.
// Return value: nr of parameters successfully parsed if no errors are generated
//               0 is returned when an error occured
int parse_command_string(verbose_definition verbose, int argc, char **argv, int argv_index, int check_only, int minrequestedparameters, char *format, ...);

/*
Applies any preprocess options specified on the command line. Returns 0 on error, 1 on success.
 */
int preprocessApplication(psrsalsaApplication *application, datafile_definition *psrdata);

/* Call when you want to add string argi from argv to the list of
filenames.  returns 1 on success, 0 on error */
int applicationAddFilename(int argi, verbose_definition verbose);

/* returns 1 if the filenames appear consecutive, or else returns
   0. If the list is empty returns 1. If argv != NULL, an error
   message is produced stating which options failed. */
int applicationFilenameList_checkConsecutive(char **argv, verbose_definition verbose);

/* Returns the number of filenames in the commandline + those
   specified in the file specified with -filelist. */
int numberInApplicationFilenameList(psrsalsaApplication *application, char **argv, verbose_definition verbose);

/* Returns NULL if the last file has been processed. */
char *getNextFilenameFromList(psrsalsaApplication *application, char **argv, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/*
  Resets current position in the input filename list (either on the
  command-line and/or in the filelist specified with the -filelist
  option).
 */
void rewindFilenameList(psrsalsaApplication *application);
//START REGION DEVELOP

//START REGION RELEASE
/*
It is assume outputname has a length MaxFilenameLength.

The output name is either from the outputname set with -output or the
extension set with -ext

Returns 1 on success, 0 on error.
 */
int getOutputName(psrsalsaApplication *application, char *filename, char *outputname, verbose_definition verbose);
//START REGION DEVELOP

//START REGION DEVELOP
//START REGION RELEASE

// Show the version information about libraries that are used
void showlibraryversioninformation(FILE *stream);


/* ***********************************
   Functions in myio.c 
   *********************************** */

/* Define own getch function, which doesn't exist in all versions of
   C. */
int pgetch(void);

/*
 Prints out txt to destination with the specified colour (if
 destination points to terminal rather than a file). The rest of the
 function works like printf.

Colours are:
1 = normal
2 = red
3 = green
4 = yellow
5 = blue
6 = magenta
7 = cyan
8 = white
*/
void fprintf_color(FILE *destination, int color, const char *format, ...);

/* Like function pgetch(), except that it uses the on the command line
   defined macro file if available instead of the keyboard. A ^ is
   interpreted as a ctrl-key. Line feeds/returns are ignored in
   macro's. If end of macro is reached, it will be closed and input
   will happen from keyboard. */
int pgetch_macro(psrsalsaApplication *application, verbose_definition verbose);

/* Tries to find out the username of the person running the
   program. Memory will be allocated. Returns 1 on success, 0
   on error. */
int getUsername(char **username, verbose_definition verbose);

/*
  Put command line in string txt with maximum length length
 */
void constructCommandLineString(char *txt, int length, int argc, char **argv, verbose_definition verbose);

/* Set size to the size of the array. Returns 1 on success, 0 on
   error. */
int getMachinename(char *hostname, int size, verbose_definition verbose);

//START REGION DEVELOP

/* Try to find out what the IP address is of hostname. ip_host is
   assumed to be (at least) 128 bytes. Returns 0 if lookup failed. The
   verbose output is sent to verbose_fptr, which can for instance be
   stdio. */
int getIPofHost(char *hostname, char *ip_host, FILE *verbose_fptr, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE
/* Parse the input string for the nth word (counting from 1). Words
   are defined to be separated by ascii code separator (i.e. ' ' for a
   space). The return value is the ptr to the start of the nth word in
   the input string (or NULL if not found). The total number of words
   in the input string is returned as nrwords. The individual words
   cannot be larger than MaxPickWordFromString_WordLength bytes, or
   else the conversion to a string using sscanf will stop. If
   replacetabs is set, all tabs in the input string are replaced by a
   space. Any trailing spaces are ignored, as well as \n and \r. The
   input string is not altered by this function. The return pointer is
   just a pointer to the start of the requested word. The string is
   not null terminated after the word. */
char * pickWordFromString(char *string, int n, int *nrwords, int replacetabs, char separator, verbose_definition verbose);

//START REGION DEVELOP

/* This function finds the first occurrence of the substring needle in
   the string haystack. This function is similar to C function strstr,
   except that it works on a block of memory with a given size, rather
   than on an ascii null terminated string. It returns a pointer to
   the beginning of the substring, or NULL if the substring is not
   found.
 */
char *searchStringInMem(const char *haystack, int haystacksize, const char *needle, int needlesize, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

// Given string orig, replace every occurance of rep with with. The
// return value is the result (memory will be allocated, so you must free
// the result if result is non-NULL). The function returns NULL if memory
// cannot be allocated.
char *str_replace(char *orig, char *rep, char *with, verbose_definition verbose);

//START REGION DEVELOP

/* This function skips lines in a file. It returns 1 on success, 0
   when reaching EOF */
int skipLinesInFile(FILE *fptr, int skiplines, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Assuming the lines in the text file (fin is the file pointer and
   the file is already opened) are less than maxlinelength characters
   long, read in the next line in txt. The variable txt should be
   allocated with maxlinelength characters. If the line starts with
   character skipChar (maybe set to '#'), the line is ignored (set to
   zero to ignore this feature). Returns 0 if there is an error
   (reached EOF?), or the actual nr of lines being read in from the
   file (i.e. 1 + nr of lines which started with skipChar). */
int ascii_file_get_next_line(FILE *fin, char *txt, int maxlinelength, int skipChar, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Assuming the lines in the text file (fin is the file pointer and
   the file is already opened) are less than maxlinelength characters
   long, determine the nr of lines in the file. If the line starts
   with character skipChar (maybe set to '#'), the line is ignored
   (set to zero to ignore this feature). If autoNrColumns is set, the
   number of columns (nrColumns) is determined as well. An error is
   generated in the nr of columns appear to change. If autoNrColumns
   is set to zero, the nr of columns is compared with that in
   *nrColumns, unless that is set to the NULL pointer. If nrColumns is
   set to a negative number, each line is expected to have at least
   -nrColumns of columns. No rewind is done. The function returns 1 on
   success, 0 on error. */
int ascii_file_stats(FILE *fin, char skipChar, long *nrlines, int maxlinelength, int autoNrColumns, int *nrColumns, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* Takes inputname, cut extension and add new extension. inputname and
   extension are not changed. A check is done if length of outputname
   will not exceed outputnamelength. Returns 1 if successful or 0 on
   error. */
int change_filename_extension(char *inputname, char *outputname, char *extension, int outputnamelength, verbose_definition verbose);

//START REGION DEVELOP
//START REGION RELEASE

/* 
  This function opens the ascii file with name fname. The first
  skiplines number of lines are ignored. If the line starts with
  character skipChar (maybe set to '#'), the line is ignored as well
  (set to zero to ignore this feature). In total there are nrColumn
  columns in the file, but you can also set autoNrColumns to determine
  this number automatically. The data is in column number colnum
  (counting from 1). Enoung memory will be allocated to contain the
  nrdatapoints points which are read in. The data is multiplied with
  scale. If read_log is set, the base 10 logarithm is stored rather
  than the actual value. If mindata, maxdata and/or avdata are set to
  something else than NULL, these will be set to the minimum, maximum
  and average values being read in (after applying the read_log
  option). If verbose_stderr is set, the verbose output is sent to the
  stderr, which could be useful if it shouldn't interfere with other
  output of program which can be expected to be redirected. Verbose
  level determines nr of spaces before output. Lines cannot exceed 10k
  lenght. The individual words cannot be larger than 1000 bytes, or
  else there will be a memory overflow.

  Returns 0 on error, otherwise 1
*/
int read_ascii_column(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, int colnum, double scale, int read_log, float **data, float *mindata, float *maxdata, float *avdata, verbose_definition verbose, int verbose_stderr);
int read_ascii_column_double(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, int colnum, double scale, int read_log, double **data, double *mindata, double *maxdata, double *avdata, verbose_definition verbose, int verbose_stderr);
int read_ascii_column_int(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, int colnum, int **data, int *mindata, int *maxdata, double *avdata, verbose_definition verbose, int verbose_stderr);
int read_ascii_column_str(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, int colnum, char ***data, verbose_definition verbose, int verbose_stderr);

//START REGION DEVELOP

// Copy file from to destination. If dest already exist, it will be overwritten.
// Returns 0 on success
int cp(const char *from, const char *dest, verbose_definition verbose);

//START REGION DEVELOP

// Read in byte values represented as hex numbers as a long double variable.
//
// Return 0: Error
// Return 1: Success
int convert_hexstring_to_longdouble(char *string, long double *value, verbose_definition verbose);

// Write a long double variable as byte values represented as hex numbers to stream
void print_longdouble_as_hexstring(FILE *stream, long double value);

//START REGION DEVELOP

/* *********************************************
   Functions and definitions in linalg.c 
   ********************************************* */

//START REGION DEVELOP

/*
  Given the n*n matrix a and collection of answer vectors b, find solutions vertors x of the matrix eq: a*x=b. Each vector b thus has n elements. In total there are m solution vectors b provided in a matrix, which is just a list of numbers. The m solution vectors are ordered as:

b[] = [b11...b1n,b21...b2n,...,bm1...bmn]

So the elements of each vector b appears consecutive in memory.

The matrix a is ordered such that the first n elements correspond to a row that is used to get the first element of the solution vectors.

The output matrix x has the same format as b (memory is never allocated and is only filled when findsolutions is nonzero.

The output matrix a_inv has the same format as a, and will contain the inverse of the matrix a. Memory should already be allocated when findinverse is nonzero, otherwise this matrix is not filled.

If preserve_a is zero, it is not guaranteed that the input array a is perserved (will not be the case when gsl is used).

If checkcaninverted is nonzero, a check is done if the matrix can be inverted. The check is always done in NR, but gsl returns by default nan's/inf's when it cannot be inverted.

Return values:
0 = ok
1 = memory allocation error
2 = Matrix a cannot be inverted
3 = Method not implemented
 */

int linalg_solve_matrix_eq(double *a, int n, int preserve_a, double *b, int m, double *x, int findsolutions, double *a_inv, int findinverse, int checkcaninverted, verbose_definition verbose);

/* If returnerror is set, the singular matrix error which terminates
   program is captured and causes the function to return 0 rather than
   1. If not set, the program always returns 1, or terminates.

   Given the n*n matrix a and collection of answer vectors b, find solutions vertors x of the matrix eq: a*x=b. Each vector b thus has n elements. In total there are m solution vectors b provided in a matrix, which is just a list of numbers. The m solution vectors are ordered as:

   b[] = [b11...b1n,b21...b2n,...,bm1...bmn]

   So the elements of each vector b appears consecutive in memory.

   The matrix a is ordered such that the first n elements correspond to a row that is used to get the first element of the solution vectors.

   The input matrix matrixa (an n X n matrix) is replaced by its inverse.
   The input matrix matrixb (an n X m matrix) is replaced by the solution (same dimensions).

   Return 0 = Success
   Return 1 = Memory allocation error
   Return 2 = Singular matrix, no solution can be determined
 */
//START REGION RELEASE
int linalg_solve_matrix_eq_gauss_jordan(double *matrixa, double *matrixb, int n, int m, verbose_definition verbose);
//START REGION DEVELOP
int linalg_solve_matrix_eq_gauss_jordan_ld(long double *matrixa, long double *matrixb, int n, int m, verbose_definition verbose);

//START REGION DEVELOP

/* *********************************************
   Functions and definitions in binarymodels.c 
   ********************************************* */

//START REGION DEVELOP

double BTmodel(long double bbat, binary_def companion);
double T2model(long double bbat, int nCompanion, binary_def *companion);
long double solve_mass_function_m2(long double mf, long double sini, long double m1);
long double getEccentricAnomaly(long double phase, long double ecc);
double BTmodel_hyper(long double bbat, int nCompanion, binary_def *companion);
double T2model_hyper(long double bbat, int nCompanion, int calc_deriv, int select_companion, binary_def *companion);

//START REGION DEVELOP

/* *********************************************
   Functions and definitions in refractionlib.c 
   ********************************************* */

/*                                                                        
Initializes variables Raytrace_Magnetosphere, Raytrace_StartPosition and
Raytrace_EqnFlag

Input:
  StartPosition.r
  StartPosition.chi
Output:
  StartPosition  (all values set to start condition)
*/
void raytrace_Init(refraction_PulsarPhysicalParameters Pulsar,
                   refraction_PositionStruct *StartPosition, int EqnFlag, double f0);

/*
Calculates the plasma density of the Magnetosphere at position Pos were the
plasma density is 1 at (r=1, chi=chi_c).
Input:
  Pos.r
  Pos.chi
Output:
  Pos.N
  Pos.NormN
*/
void plasmaDensity(refraction_PositionStruct *Pos);

/*
Output:
  FRD.Succesful
  FRD.theta_f
Return values:
  0 - Ok
 -1 - Out of memory
*/
int calculateFRD(refraction_PulsarPhysicalParameters Pulsar, int Raytrace_EqnFlag,
		double r_max, double MaximalError, double MinimalStepSize,
		refraction_FinalRayDirection *FRD, int verbose);

/*
 0 = ok
-1 = Cannot open FileName
-2 = More than MaxNrConesRefraction cones defined
-3 = No number after Cone
-4 = Cannot set relative values first plasma cone
-5 = Cannor set chi_cn without chi_c
-6 = Unknown variable
*/
int refraction_ReadInfoFile(char *FileName, int *NrRays, int *EquationSet, 
                 double *LongitudeRange, int *FileNr, int Verbose, 
                 refraction_FinalRayDirection *FRD, refraction_PulsarPhysicalParameters *Pulsar, 
                 refraction_LineOfSight *LOS, refraction_InfoReadedStruct *IR);


double dlnNdchi(refraction_PositionStruct Pos);

/*
 1 = ok
 0 = Not ok
*/
int refraction_CheckEnoughInfoForPulseProfile(refraction_InfoReadedStruct IR);

/* 
Crates an equispaced beam pattern (in chi_s) from chi1 to chi2 using NrPoints points. 
Input:
  chi1 [chi_c cone 1]
  chi2 [chi_c cone 1]
  NrPoints
Output:
  FRD.NrPoints
  FRD.chis
Return:
  0 = Ok
 -1 = Out of memory
*/
int refraction_CreateEquiSpacedBeamPattern(int NrPoints, double chi1, double chi2,
                                refraction_FinalRayDirection *FRD);

/* 
Crates an equispaced f0 pattern. 
Return:
  0 = Ok
 -1 = Out of memory
*/
int refraction_CreateEquiSpacedf0Pattern(refraction_PulsarPhysicalParameters Pulsar, 
                              refraction_FinalRayDirection *FRD);

/*
Input:
  FRD.NrPoints
  FRD.EmissionHeight (if FRD.DisableVariableEmissionHeight != 0)
  Pulsar.NrCones
  Pulsar.Cone[].chi_c
  Pulsar.Cone[].chi_cn
Output:
  FRD.r0
  FRD.chi0
Return:
  0 = Ok
 -1 = Out of memory  
*/
int refraction_CalculateEmissionHeight(refraction_PulsarPhysicalParameters Pulsar, refraction_FinalRayDirection *FRD);

/*
Return values:
   0 = Ok
  -1 = Cannot open file
*/
int refraction_outputFRD(refraction_FinalRayDirection FRD, refraction_PulsarPhysicalParameters Pulsar, 
              int FileNr);

/*
Input:
  LOS.NrPoints
  LOS.alpha
  LOS.zeta
Output:
  LOS.pl
  LOS.phi
  LOS.theta
Return values:
  0 = ok
 -1 = out of memory
*/
int refraction_calculateLOS(refraction_LineOfSight *LOS);


/*
Return values
  0 = Ok
 -1 = Cannot Open File
*/
int refraction_outputLOS(refraction_LineOfSight LOS, int FileNr);

/*
Input:
  Stack.NrPulses
Output:
Return values
  0 = ok
 -1 = out of memory
*/
int refraction_calculatePulseStack(refraction_PulsarPhysicalParameters Pulsar, refraction_FinalRayDirection FRD, 
                       refraction_LineOfSight LOS, double LongitudeRange, refraction_PulseProfileStack *Stack);

void refraction_NormalizeStack(refraction_PulseProfileStack *Stack);


/*
Return Values:
   0 = Ok
  -1 = Cannot open file
*/
int refraction_outputStack(refraction_PulseProfileStack Stack, refraction_LineOfSight LOS, 
                refraction_PulsarPhysicalParameters Pulsar, int FileNr);

void freeStack(refraction_PulseProfileStack Stack);

void FreeFRD(refraction_FinalRayDirection FRD);
void freeLOS(refraction_LineOfSight LOS);


/*
 1 = ok
 0 = Not ok
*/
int refraction_CheckEnoughInfoForRaytrajectories(refraction_InfoReadedStruct IR);

/* 
Calculates trajectory of a single ray
Input:
  Ray.StartPoint.r
  Ray.StartPoint.chi
  Pulsar.Magnetosphere
  Raytrace_EqnFlag    - 0 = Weltevrede, 1 = Petrova, 2 = Petrova with bug
  r_max               - Integrate to this distance from the star
  MaximalError   - Maximal num. error per integration step (for example 1.0E-8)
  MinimalStepSize - Minamal allowed stepsize
  GetExtendedInformation - 0 = Don't calculate N, NormN, NrSolutions and eta 
                               for all points of ray trajectory
Output:
  Ray
Return values:
  0 - ok
 -1 - nummerical error
 -2 - cannot allocate memory
*/
int refraction_calculateSingleRay(refraction_RayParameters *Ray,
                      refraction_PulsarPhysicalParameters Pulsar, double r_max,
                      double MaximalError, double MinimalStepSize,
                      int Raytrace_EqnFlag, int GetExtendedInformation,
                      double f0, double gamma, int verbose);

/*
 1 = ok
 0 = Not ok
*/
int refraction_CheckEnoughInfoForSingleRay(refraction_InfoReadedStruct IR);

/* Output ray trajectory
 Input:
   Ray      - The ray
   Pulsar   - The physical parameters
   FileNr   - Output to file called ray99.dat if FileNr = 99
 Return values:
   0 - ok
  -1 - Cannot create file
*/
int refraction_outputSingleRay(refraction_RayParameters Ray, refraction_PulsarPhysicalParameters Pulsar, 
                    int FileNr);


/* 
Free all allocated memory of ray trajectory
*/
void FreeRay(refraction_RayParameters Ray);


/* ***********************************
   Functions in nr_f.c 
   *********************************** */

/*
Return values:
  0 - No errors
 -1 - "Step size too small in rkqs"
*/

int rkqs_adj(float y[], float dydx[], int n, float *x, float htry, float eps,
	float yscal[], float *hdid, float *hnext,
	     void (*derivs)(float, float [], float []));

/*
Return values:
  0 - No errors
 -1 - "Step size too small in odeint"
 -2 - "Too many steps in routine odeint"
*/
int odeint_adj(float ystart[], int nvar, float x1, float x2, float eps, float h1,
	float hmin, int *nok, int *nbad,
	void (*derivs)(float, float [], float []),
	int (*rkqs_adj)(float [], float [], int, float *, float, float, float [],
		    float *, float *, void (*)(float, float [], float [])));


/* ***********************************
   Functions in nr_d.c 
   *********************************** */
double *vector_d(long nl, long nh);
void free_vector_d(double *v, long nl, long nh);


void rkck_d(double y[], double dydx[], int n, double x, double h, double yout[],
	    double yerr[], void (*derivs)(double, double [], double []));

int rkqs_adj_d(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	       void (*derivs)(double, double [], double []));

int rkqs_d(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	   void (*derivs)(double, double [], double []));

int odeint_adj_d(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	int (*rkqs_d)(double [], double [], int, double *, double, double, double [],
		      double *, double *, void (*)(double, double [], double [])));

void ludcmp_d(double **a, int n, int *indx, double *d);
void lubksb_d(double **a, int n, int *indx, double b[]);

double **matrix_d(long nrl, long nrh, long ncl, long nch);
void free_matrix_d(double **m, long nrl, long nrh, long ncl, long nch);






/* ***********************************
   Functions in nr_ld.c 
   *********************************** */

long double *vector_ld(long nl, long nh);
long double **matrix_ld(long nrl, long nrh, long ncl, long nch);
void free_vector_ld(long double *v, long nl, long nh);
void free_matrix_ld(long double **m, long nrl, long nrh, long ncl, long nch);
long double amotry_ld(long double **p, long double y[], long double psum[], int ndim,
		      long double (*funk)(long double []), int ihi, long double fac);
void mrqmin_ld(long double x[], long double y[], long double sig[], int ndata, long double a[], int ia[],
	int ma, long double **covar, long double **alpha, long double *chisq,
	       void (*funcs)(long double, long double [], long double *, long double [], int), long double *alamda);
/* If returnerror is set, the singular matrix error which terminates
   program is captured and causes the function to return 0 rather than
   1. If not set, the program always returns 1, or terminates.
 */
int gaussj_ld(long double **a, int n, long double **b, int m, int returnerror);


