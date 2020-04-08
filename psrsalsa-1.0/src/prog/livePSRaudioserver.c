#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "cpgplot.h"
#include <netdb.h>          /* hostent struct, gethostbyname() */
#include <arpa/inet.h>      /* inet_ntoa() to format IP address */
#include <netinet/in.h>     /* in_addr structure */
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>
#include <inttypes.h>
#include "psrsalsa.h"

#ifdef __APPLE__
  #ifndef MSG_NOSIGNAL
// MSG_NOSIGNAL doesn't appears to be defined on Mac OSX. By defining
// it to be zero, this option of send() is effectively ignored.  I
// don't really understand what this flag is supposed to do, so it
// might cause problems on Macs....
     #define MSG_NOSIGNAL 0
  #endif
#endif

// One package is sent across every 50 ms by default. All numbers assume one freq channel. In reallity, the output is larger by a factor nrchan per 50 ms block. Frequency channel is the inner loop.
struct {
  double packetDuration;         // A packet of data will be sent repeatedly every packetDuration ms.
  long nrSamplesPerPackage;      // Each 50 ms package sent contains this number of floats
  //  long maxAllowedOversampling;   // The input data may have up to this factor extra floats per package
  //  long nrFreqChan;               // Default is 1 frequency channel being transmitted
}packet_def;

#define MaxPlotLength     80000        // In samples, this corresponds to 10 seconds
#define MaxNrTracesInPlot   100        // The maximum amount of pulses/lines in plot, not sure why really used
#define MaxAutopeakValue    100        // Maximum nr of cycles to use for determining peak values
#define MaxNrFreqChan      1000        // Maxumum number of frequency channels to support
#define MAXPENDING            8        // Max connection requests, not sure if really used as there is a simple 1-1 connection

/* Stuff TCP/IP server with PuMa II */
int serversock, clientsock;
unsigned int portnr_puma2;
struct sockaddr_in audioserver, puma2client;
FILE *input_file;

void Die(char *mess) { perror(mess); exit(1); }     /* Let program die, and say why */
void NotDie(char *mess) { perror(mess); }           /* Stay alive, but say what failed */
void connectPulsarSignalSender();                               /* (Re)connect server for PuMa II */
void Readparameters(char *filename, float *baseline, float *peak, int *resample, int *autobaseline, int *autopeak, float *period, int *plot_length_samples, float *tsamp_orig, char *title, int *NrPlottedPulses, float *delay_plot_seconds, long *pulsenumber, int *plotOnePulseAtTheTime, float *sonification_pwrlaw, int *reset_plot_after_new_connection, int *autopeak_nrinringbuffer, int *autopeak_curbufferpos);
void showParametersExplanation();
void plotPulsarSignal(float *Iserie, int pulse, int NrBins, int NrPulses, float baseline, float peak, float xviewportSize, float yviewportSize, char *title, int startbin, int endbin);
void plotSignalCheck(float *Iserie, float *Iserie2, int NrBins, float baseline, float peak, float xviewportSize, float yviewportSize, int mode, char *title, char *pulsarname);
void plotPulsarSignal_reset(int NrBins, int NrPulses, float xviewportSize, float yviewportSize, char *title, char *pulsarname);
float calc_median(float *buffer, int nrentries);


int main(int argc, char **argv)
{
  int streamsock, stopSendingData, plotOnePulseAtTheTime;
  struct sockaddr_in streamserver;
  char ip_host_stream[128], ip_host_localhost[128];
  unsigned int portnr_stream;
  int received, maxy_undefined, closeAtEndTransmission;
  long nrfreqchan;
  char hostname[100], pulsarname[100];

  char PlotDevice[100], PlotDevice2[100], title[100];
  float maxy, xviewportSize, yviewportSize, baseline, peak, average, period, delay_plot_seconds;
  double cur_plot_sample_time, orig_sample_time, stack_time;
  int deviceID, deviceID2, resample, autobaseline, autopeak, NrPlottedPulses, firsttime;
  long i, j; //, k;
  // float miny, 
  float *timeseries_package, *timeseries, *timeseries_sonified, tmpfloat, tsamp_orig, sonification_pwrlaw;
  long cur_plot_sample, cur_plot_sample_drawn, nosonic, InputNrSamplesPerTime, plotDevice2_defined, pulsenumber;
  int plot_length_samples, plot_after_this_number_of_samples; // , firsttimeconnected_flag
  int bytes, output_error_stream, output_error_puma2, output_stdout_only, output_stdout_aswell; //, pass
  int reset_plot_after_new_connection, streamserver_set;
  int16_t *audiobuf;
  float autopeak_ringbuffer[MaxAutopeakValue];
  int autopeak_nrinringbuffer, autopeak_curbufferpos;
  verbose_definition noverbose, verbose;

  /*  fprintf(stderr, "MSG_NOSIGNAL=%x\n", MSG_NOSIGNAL); */

  packet_def.packetDuration = 50;
  packet_def.nrSamplesPerPackage = 400;
  portnr_puma2 = 3000;
  portnr_stream = 3001;
  ip_host_stream[0] = 0;
  input_file = NULL;

  nrfreqchan = 1;
  strcpy(PlotDevice, "?");
  xviewportSize = 0.9;
  yviewportSize = 0.8;

  baseline = 0;
  peak = 1;
  resample = 1;
  pulsenumber = 0;
  autobaseline = 0;
  autopeak = 0;
  output_stdout_only = 0;
  nosonic = 0;
  plot_length_samples = 8000;
  tsamp_orig = 0.000125;
  period = 1;
  title[0] = 0;
  plotDevice2_defined = 0;
  delay_plot_seconds = 0;
  //  firsttimeconnected_flag = 1;
  plot_after_this_number_of_samples = packet_def.nrSamplesPerPackage;
  plotOnePulseAtTheTime = 1;
  sonification_pwrlaw = 2;
  pulsarname[0] = 0;
  reset_plot_after_new_connection = 1;
  autopeak_nrinringbuffer = 0;
  autopeak_curbufferpos = 0;
  cleanVerboseState(&noverbose);
  cleanVerboseState(&verbose);
  verbose.verbose = 1;
  output_stdout_aswell = 0;
  streamserver_set = 0;
  closeAtEndTransmission = 0;
  stopSendingData = 0;

  if(argc < 2) {
    fprintf(stderr, "Usage:livePSRaudioserver  parameterfile\n\n");
    fprintf(stderr, "This program listens on the port specified with the -p option for the incomming signal from the pulsar backend, or from the program livePSRsender. The input signal is assumed to be at a fixed rate of 8000 floating points per second. This signal is plotted (using parameters from parameterfile, which can be adjusted while the software is running), repacked to 16 bit signed integer values, and send out again to the stdout (so you can pipe it to an audio player), or to an other computer where the audio is generated (the stream server).\n\n");
    fprintf(stderr, "An example of the contents of the parameterfile (everything before : is ignored, information should be on the correct line number):\n");
    showParametersExplanation();
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "-h \"host name\" Set host name of stream server.\n");
    fprintf(stderr, "-i alternatively, set IP adress.\n");
    fprintf(stderr, "-p portnumber on which the signal is received. Default is %d.\n", portnr_puma2);
    fprintf(stderr, "-file Instead, read stream data from this binary file generated by livePSRsender.\n");
    fprintf(stderr, "-P portnumber to stream server to use. Default is %d.\n", portnr_stream);
    fprintf(stderr, "-D Set PGPLOT device, default is %s.\n", PlotDevice);    
    fprintf(stderr, "-D2 Set PGPLOT device for stack (and enables it).\n");    
    fprintf(stderr, "-s Send output to stdout instead of streamserver.\n");
    fprintf(stderr, "-S Send output to stdout as well as streamserver.\n");
    fprintf(stderr, "-nosonic  Do not display the sonification signal.\n");
    fprintf(stderr, "-x Set viewport x-range from 0.1 to this value (def is 0.9).\n");
    fprintf(stderr, "-y Set viewport y-range from 0.1 to this value (def is 0.9).\n");
    fprintf(stderr, "-close    Close after incomming signal stops.\n");
    return 0;
  }

  fprintf(stderr, "Try to find IP of computer where livePSRaudioserver is running\n");
  if(getMachinename(hostname, 100, noverbose) == 0)
    strcpy(hostname, "localhost");
  else
    fprintf(stderr, "  Machine identified as %s\n", hostname);
  if(getIPofHost(hostname, ip_host_localhost, stderr, verbose) == 0) {
    if(strcmp(hostname, "localhost") != 0) {
      strcpy(hostname, "localhost");
      fprintf(stderr, "Try out machine name localhost\n");
      getIPofHost(hostname, ip_host_localhost, stderr, verbose);     
    }
  }
  fprintf(stderr, "\n");

  if(argc > 2) {
    for(i = 1; i < argc-1; i++) {
      if(strcmp(argv[i], "-D") == 0) {
	strcpy(PlotDevice, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-D2") == 0) {
	strcpy(PlotDevice2, argv[i+1]);
	plotDevice2_defined = 1;
	i++;
      }else if(strcmp(argv[i], "-file") == 0) {
	input_file = fopen(argv[i+1], "r");
	if(input_file == NULL) {
	  fprintf(stderr, "Cannot open input file.");
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-x") == 0) {
	xviewportSize = atof(argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-y") == 0) {
	yviewportSize = atof(argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-nosonic") == 0) {
	nosonic = 1;
	yviewportSize = 0.85;
      }else if(strcmp(argv[i], "-close") == 0) {
	closeAtEndTransmission = 1;
      }else if(strcmp(argv[i], "-h") == 0) {
	strcpy(hostname, argv[i+1]);
	fprintf(stderr, "Try to find IP of the specified stream server: %s\n", hostname);
	if(getIPofHost(hostname, ip_host_stream, stderr, verbose) == 0) {
	  fprintf(stderr, "Cannot determine IP address of the stream server, specify this with the -i option.");
	}
	fprintf(stderr, "\n");
	streamserver_set = 1;
	i++;
      }else if(strcmp(argv[i], "-i") == 0) {
	strcpy(ip_host_stream, argv[i+1]);      
	streamserver_set = 1;
	i++;
      }else if(strcmp(argv[i], "-p") == 0) {
	portnr_puma2 = atoi(argv[i+1]);      
	i++;
      }else if(strcmp(argv[i], "-P") == 0) {
	portnr_stream = atoi(argv[i+1]);      
	i++;
      }else if(strcmp(argv[i], "-s") == 0) {
	output_stdout_only = 1;
      }else if(strcmp(argv[i], "-S") == 0) {
	output_stdout_aswell = 1;
      }else {
        fprintf(stderr, "Unknown option: %s\n", argv[i]);
	return 0;
      }
    }
  }
 
  if(output_stdout_aswell && output_stdout_only) {
    fprintf(stderr, "You cannot use -s and -S at the same time.\n");
    return 0;
  }

  if(streamserver_set && output_stdout_only) {
    fprintf(stderr, "If you want to use a stream server and stdout at the same time, use the -S option, not the -s option.\n");
    return 0;
  }

  /* Check if the stream server is specified, or if stdout is used */
  if(output_stdout_only == 0 || output_stdout_aswell) {
    if(ip_host_stream[0] == 0) {                                     
      fprintf(stderr, "Specify the computer on which audio output server is running (using -h or -i), or use -s.\n");
      return 0;
    }
  }

  // Open pgplot device
  deviceID = cpgopen(PlotDevice);
  cpgask(0);
  cpgpage();
  cpgsch(1);
  cpgslw(1);

  if(plotDevice2_defined) {
    deviceID2 = cpgopen(PlotDevice2);
    cpgask(0);
    cpgpage();
    cpgsch(1);
    cpgslw(1);
  }

  // Store the incomming 50ms packages
  timeseries_package = (float *)malloc(MaxNrFreqChan*packet_def.nrSamplesPerPackage*2*sizeof(float));
  // Store the whole specified period
  timeseries = (float *)malloc(MaxNrFreqChan*2*MaxPlotLength*sizeof(float));
  // The rescaled/baselined version, maybe also squared
  timeseries_sonified = (float *)malloc(2*MaxPlotLength*sizeof(float));
  // The created audio stream
  audiobuf = (int16_t *) calloc(packet_def.nrSamplesPerPackage, sizeof(int16_t));

  if(timeseries_package == NULL  || audiobuf == NULL || timeseries == NULL || timeseries_sonified == NULL) {
    fprintf(stderr, "Memory allocation error.\n");
    return 0;
  }
  for(i = 0; i < plot_length_samples; i++) {
    timeseries[i] = 0;
    timeseries_sonified[i] = 0;
  }

  /* Create socket for streaming server */

  if(output_stdout_only == 0) {
    fprintf(stderr, "Setup connection to stream server.\n");
    /*  if ((streamsock = socket(PF_INET, SOCK_DGRAM, IPPROTO_UDP)) < 0) {      */
    if ((streamsock = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) {       
      Die("  Failed to create socket");
    }else {
      fprintf(stderr, "  Created socked %d\n", streamsock);
    }
    memset(&streamserver, 0, sizeof(streamserver));       /* Clear struct */
    streamserver.sin_family = AF_INET;                   /* Internet/IP */
    streamserver.sin_addr.s_addr = inet_addr(ip_host_stream);   /* IP address */
    streamserver.sin_port = htons(portnr_stream);               /* server port */
    fprintf(stderr, "  Trying to establish connection to streaming server %s at port %d\n", ip_host_stream, portnr_stream);    
    if (connect(streamsock, (struct sockaddr *) &streamserver, sizeof(streamserver)) < 0) {
      NotDie("  Failed to connect with streamserver");
      fprintf(stderr, "  livePSRaudioserver will keep trying to connect to stream server\n");
      stopSendingData = 3;      /* Will not start with sending data, but keep trying to connect */
    }else {
      fprintf(stderr, "  Established connection with streamserver\n");    
      stopSendingData = 0;      /* Default is an established connection to audio server */
    }
  }

  NrPlottedPulses = 10;
  fprintf(stderr, "\nReading in parameter file for the first time\n");
  Readparameters(argv[argc-1], &baseline, &peak, &resample, &autobaseline, &autopeak, &period, &plot_length_samples, &tsamp_orig, title, &NrPlottedPulses, &delay_plot_seconds, &pulsenumber, &plotOnePulseAtTheTime, &sonification_pwrlaw, &reset_plot_after_new_connection, &autopeak_nrinringbuffer, &autopeak_curbufferpos);

  fprintf(stderr, "\n");
  if(input_file == NULL) {
    connectPulsarSignalSender();
  }

  //  k = 0;
  received = 0;
  //  pass = 1;
  output_error_puma2 = 1;
  output_error_stream = 1;
  cur_plot_sample = 0;
  cur_plot_sample_drawn = 0;
  cur_plot_sample_time = 0;
  orig_sample_time = 0;
  stack_time = 0;
  pulsenumber = 0;
  average = 0;
  maxy_undefined = 1;
  fprintf(stderr, "Codes are:\n  'r' read package ok\n  'w' write audio ok\n  '(number)' Received package size is wrong\n  'S' try open socket\n  'C' try reconnect to audio server\n");
  //Skip plotting first pulse, as baseline is not yet set.
  if(autobaseline)
    pulsenumber = -1;
  firsttime = 1;
  do {
    bytes = 0;
    fflush(stderr);
    if(input_file == NULL) {
      bytes = recv(clientsock, &timeseries_package[received/4], packet_def.nrSamplesPerPackage*4, 0);
      if(bytes < 1 && closeAtEndTransmission) {
	close(clientsock);
	close(serversock);
	close(streamsock);
	cpgend();
	fprintf(stderr, "\nTransmission ended\n");
	exit(0);
      }
    }else {
      bytes = fread(&timeseries_package[received/4], 1, packet_def.nrSamplesPerPackage*4, input_file);
      if(bytes == 0) {  // Reading from input file failed, end somewhat graciously
	//	close(clientsock);  // This is the incomming network connection
	//	close(serversock);  // This is the incomming network connection
	if(output_stdout_only == 0) {
	  close(streamsock);
	}
	cpgend();
	fprintf(stderr, "\nReached EOF\n");
	exit(0);
      }
    }
    if(bytes < 1) {
      if(output_error_puma2 < 5) {
	NotDie("\nConnection PuMa II failed"); 
	output_error_puma2++;
      }else {
	fprintf(stderr, "Going to restart server for PuMa II.\n");
	close(clientsock);
	close(serversock);
	connectPulsarSignalSender();
	output_error_puma2 = 1;
	received = 0;
	if(reset_plot_after_new_connection) {
	  pulsenumber = 0;
	  cur_plot_sample_drawn = 0;
	  cur_plot_sample = 0;
	}
	pulsarname[0] = 0;
      }
    }else {                               // Received package with pulsar signal
      output_error_puma2 = 1;
      if(bytes == packet_def.nrSamplesPerPackage*4)
	fprintf(stderr, "r");
      else
	fprintf(stderr, "(%d)", bytes);
      /*
      if(bytes > 0 && firsttimeconnected_flag) {
	int buffer_index;
	fprintf(stderr, "\nFirst floats to read from data stream are:");
	for(buffer_index = 0; buffer_index < packet_def.nrSamplesPerPackage*4/4; buffer_index++) {
	  fprintf(stderr, " %e", timeseries_package[buffer_index]);
	}
	fprintf(stderr, "\n");
	firsttimeconnected_flag = 0;
      }
      */
      received += bytes;
      if(received >= packet_def.nrSamplesPerPackage*4) {                 // Check if a complete 50 ms package has arrived
	int isheader;
	char txt[100];
	char *hdr_ptr;

	hdr_ptr = (char *)(&timeseries_package[0]);

	if(firsttime) {
	  fprintf(stderr, "\nFirst byte values received are: ");
	  for(i = 0; i < 20; i++) {
	    fprintf(stderr, " %d", hdr_ptr[i]);
	  }
	  fprintf(stderr, "\n\n");
	  firsttime = 0;
	}

	isheader = 0;
	strcpy(txt, "livePSRsender header");
	if(memcmp(txt, hdr_ptr, strlen(txt)) == 0) {
	  isheader = 1;
	}

	if(isheader) {
	  fprintf(stderr, "\nReading header package\n");
	  hdr_ptr = (char *)(&timeseries_package[0]);
	  hdr_ptr += strlen(txt)+1;
	  sscanf(hdr_ptr, "%s", pulsarname);
	  fprintf(stderr, "  Source name is: %s\n", pulsarname);
	  hdr_ptr += strlen(pulsarname)+1;
	  sscanf(hdr_ptr, "%s", txt);
	  sscanf(txt, "%ld", &nrfreqchan);
	  fprintf(stderr, "  Nr of frequency channels: %ld\n", nrfreqchan);
	  hdr_ptr += strlen(txt)+1;
	  if(nrfreqchan > MaxNrFreqChan) {
	    fprintf(stderr, "Maximum number of frequency channels exceeded\n");
	    exit(0);
	  }
	  
	  received = 0;                                     // Start writing at start again, otherwise header will remain in input
	}else {
	  if(resample > 0) {                                   // Resample data if requested (note, package is still packet_def.nrSamplesPerPackage samples, but samples are duplicated)
	    for(j = 0; j < packet_def.nrSamplesPerPackage/resample; j++) {
	      tmpfloat = 0;
	      for(i = 0; i < resample; i++)
		tmpfloat += timeseries_package[j*resample+i];
	      tmpfloat /= (float)resample;
	      for(i = 0; i < resample; i++)
		timeseries_package[j*resample+i] = tmpfloat;
	    }
	  }
	  for(j = 0; j < packet_def.nrSamplesPerPackage; j++) {           // Do the scaling and convert data to audio stream
	    average += timeseries_package[j];
	    //	    fprintf(stderr, "XXXX %f %f\n", timeseries_package[j], average);
	    if(maxy_undefined) {
	      maxy = timeseries_package[j];
	      maxy_undefined = 0;
	    }
	    if(timeseries_package[j] > maxy)
	      maxy = timeseries_package[j];
	    if(cur_plot_sample >= 2*MaxPlotLength-10) {
	      fprintf(stderr, "Buffer size exceeded\n");
	      exit(0);
	    }
	    if(cur_plot_sample < 0) {
	      fprintf(stderr, "Buffer size undeflow????\n");
	      exit(0);
	    }
	    timeseries[cur_plot_sample] = timeseries_package[j];
	    timeseries_sonified[cur_plot_sample] = timeseries_package[j] - baseline;
	    timeseries_sonified[cur_plot_sample] /= peak;
	    if(timeseries_sonified[cur_plot_sample] > 1)
	      timeseries_sonified[cur_plot_sample] = 1;
	    if(timeseries_sonified[cur_plot_sample] < -1)
	      timeseries_sonified[cur_plot_sample] = -1;
	    if(sonification_pwrlaw != 1)
	      timeseries_sonified[cur_plot_sample] = pow(timeseries_sonified[cur_plot_sample], sonification_pwrlaw);
	    /*	    audiobuf[j] = (int16_t) lrintf( 32768.0f*timeseries_package2[nrfloats+j] ); */
	    audiobuf[j] = (int16_t)  32768.0f*timeseries_sonified[cur_plot_sample];
	    cur_plot_sample++;
	    
	    /* Do pulse-shift correction because of difference data-rate actual data/sended data */
	    cur_plot_sample_time += 1.0/8000.0;           
	    InputNrSamplesPerTime = 0.001*packet_def.packetDuration/tsamp_orig;  /* Used in the sender software */
	    orig_sample_time += tsamp_orig*((double)InputNrSamplesPerTime/400.0);
	    if(8000.0*(orig_sample_time - cur_plot_sample_time)  > 1) {
	      timeseries[cur_plot_sample] = timeseries[cur_plot_sample-1];
	      timeseries_sonified[cur_plot_sample] = timeseries_sonified[cur_plot_sample-1];
	      cur_plot_sample++;
	      orig_sample_time -= 1.0/8000.0;      
	    }else if(8000.0*(orig_sample_time - cur_plot_sample_time)  < -1) {
	      cur_plot_sample--;  
	      orig_sample_time += 1.0/8000.0;      
	    }
	    
	    /* Correction pulse-shift because of non-integer number of samples in pulse period */
	    stack_time += period/(double)plot_length_samples;
	    if(8000.0*(cur_plot_sample_time - stack_time)  > 1) {
	      timeseries[cur_plot_sample] = timeseries[cur_plot_sample-1];
	      timeseries_sonified[cur_plot_sample] = timeseries_sonified[cur_plot_sample-1];
	      cur_plot_sample++;
	      stack_time += 1.0/8000.0;      
	    }else if(8000.0*(cur_plot_sample_time - stack_time)  < -1) {
	      cur_plot_sample--;  
	      stack_time -= 1.0/8000.0;      
	    }
	    
	    if(cur_plot_sample >= plot_length_samples +  delay_plot_seconds*packet_def.nrSamplesPerPackage/(0.001*packet_def.packetDuration)) {
	      average /= (float)plot_length_samples;
	      fprintf(stderr, "\n Av = %.3e (%.0f%%) ", average, 100.0*(average-baseline)/peak);
	      fprintf(stderr, " Pk = %.3e (%.0f%%)\n", maxy-baseline, 100.0*(maxy-baseline)/peak);
	      plotSignalCheck(timeseries, timeseries_sonified, plot_length_samples, baseline, peak, xviewportSize, yviewportSize, nosonic, title, pulsarname);
	      /* blaat */
	      if(plotDevice2_defined) {
		cpgslct(deviceID2);
		if(pulsenumber == 0 && plotOnePulseAtTheTime) {
		  plotPulsarSignal_reset(plot_length_samples, NrPlottedPulses, xviewportSize, yviewportSize, title, pulsarname);
		}
		if(pulsenumber >= 0) {
		  if(plotOnePulseAtTheTime)
		    plotPulsarSignal(timeseries, pulsenumber, plot_length_samples, NrPlottedPulses, baseline, peak, xviewportSize, yviewportSize, title, 0, plot_length_samples-1);
		  else
		    plotPulsarSignal(timeseries, pulsenumber, plot_length_samples, NrPlottedPulses, baseline, peak, xviewportSize, yviewportSize, title, cur_plot_sample_drawn, plot_length_samples-1);
		}
		pulsenumber++;
		if(pulsenumber >= NrPlottedPulses)
		  pulsenumber = 0;
		cpgslct(deviceID);
	      }
	      if(cur_plot_sample-plot_length_samples > 0) {
		memcpy(timeseries, &timeseries[plot_length_samples], 4*(cur_plot_sample-plot_length_samples));
	      }
	      cur_plot_sample -= plot_length_samples;
	      Readparameters(argv[argc-1], &baseline, &peak, &resample, &autobaseline, &autopeak, &period, &plot_length_samples, &tsamp_orig, title, &NrPlottedPulses, &delay_plot_seconds, &pulsenumber, &plotOnePulseAtTheTime, &sonification_pwrlaw, &reset_plot_after_new_connection, &autopeak_nrinringbuffer, &autopeak_curbufferpos);
	      if(autobaseline) {
		baseline = average;
		fprintf(stderr, "Change baseline value to %f\n", baseline);
	      }
	      if(autopeak) {
		if(autopeak_nrinringbuffer == autopeak) { //Ringbuffer full
		  //		  fprintf(stderr, "Write at ringbuffer pos %d: %f (%f %f)\n", autopeak_curbufferpos, maxy-baseline, maxy, baseline);
		  autopeak_ringbuffer[autopeak_curbufferpos++] = maxy-baseline;
		  if(autopeak_curbufferpos == autopeak)
		    autopeak_curbufferpos = 0;
		}else {
		  //		  fprintf(stderr, "Write at ringbuffer pos %d: %f (%f %f)\n", autopeak_nrinringbuffer, maxy-baseline, maxy, baseline);
		  autopeak_ringbuffer[autopeak_nrinringbuffer++] = maxy-baseline;
		  autopeak_curbufferpos = 0;
		}
		peak = calc_median(autopeak_ringbuffer, autopeak_nrinringbuffer);
		fprintf(stderr, "Change peak value to %f\n", peak);
	      }
	      average = 0;
	      maxy_undefined = 1;
	      cur_plot_sample_drawn = 0; 
	    }else if(plotOnePulseAtTheTime == 0 && (cur_plot_sample >= cur_plot_sample_drawn +  delay_plot_seconds*packet_def.nrSamplesPerPackage/(0.001*packet_def.packetDuration) + plot_after_this_number_of_samples)) {
	      
	      int max_sample_to_draw = cur_plot_sample- delay_plot_seconds*packet_def.nrSamplesPerPackage/(0.001*packet_def.packetDuration);
	      
	      if(plotDevice2_defined) {
		cpgslct(deviceID2);
		if(pulsenumber == 0 && cur_plot_sample_drawn == 0) {
		  plotPulsarSignal_reset(plot_length_samples, NrPlottedPulses, xviewportSize, yviewportSize, title, pulsarname);
		}
		if(pulsenumber >= 0)
		  plotPulsarSignal(timeseries, pulsenumber, plot_length_samples, NrPlottedPulses, baseline, peak, xviewportSize, yviewportSize, title, cur_plot_sample_drawn, max_sample_to_draw-1);
		cpgslct(deviceID);
		cur_plot_sample_drawn = max_sample_to_draw;
	      }
	    }
	  }
	  if(stopSendingData == 0 || output_stdout_only != 0 || output_stdout_aswell != 0) {
	    if(output_stdout_only == 0 || output_stdout_aswell != 0) {
	      if(send(streamsock, audiobuf, packet_def.nrSamplesPerPackage*4/2, MSG_NOSIGNAL) != packet_def.nrSamplesPerPackage*4/2) { 
		NotDie("\nFailed to send bytes to streamserver");
		stopSendingData = 1;
	      }else {
		fprintf(stderr, "w");
	      }
	    }
	    if(output_stdout_only != 0 || output_stdout_aswell != 0) {
	      fwrite(audiobuf, sizeof(int16_t), packet_def.nrSamplesPerPackage, stdout);
	      fprintf(stderr, "w");
	    }
	  }else if(stopSendingData == 1) {          /* Connection failed, so close socked */
	    close(streamsock);
	    stopSendingData++;
	  }else if(stopSendingData == 2) {        /* Try to open new socket */
	    fprintf(stderr, "S");
	    if ((streamsock = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) {
	      NotDie("\nFailed to create socket");
	    }else {
	      fprintf(stderr, "\nCreated socked %d\n", streamsock);
	      memset(&streamserver, 0, sizeof(streamserver));       /* Clear struct */
	      streamserver.sin_family = AF_INET;                   /* Internet/IP */
	      streamserver.sin_addr.s_addr = inet_addr(ip_host_stream);   /* IP address */
	      streamserver.sin_port = htons(portnr_stream);               /* server port */
	      stopSendingData++;
	    }
	  }else if(stopSendingData == 3) {       /* Try to reconnect to audio server */
	    fprintf(stderr, "C");
	    if (connect(streamsock, (struct sockaddr *) &streamserver, sizeof(streamserver)) < 0) {
	      if(output_error_stream == 1)
		NotDie("Failed to connect with server");
	      output_error_stream = 0;
	    }else {
	      fprintf(stderr, "\nConnected\n"); 
	      stopSendingData = 0;
	      output_error_stream = 1;
	      received = 0;	      
	    }
	  }
	  received = 0;
	}
      }   // End of if statement if full package is received
    }           /* End of if reveived package */
  }while(2 > 1);    /* Infinite loop */

  close(clientsock);
  close(serversock);
  close(streamsock);
  cpgend();
  /*  free(aubuffer); */
  /*  free(timeseries_package); */
  return 0;
}

void showParametersExplanation()
{
  fprintf(stderr, "  title (use '_' instead of space):               Live@JB_of_pulsar_%%P\n");
  fprintf(stderr, "  autobaseline_flag (0 = off, 1 = on):            1\n");
  fprintf(stderr, "  initial baseline value:                         0.0\n");
  fprintf(stderr, "  autopeak_flag (> 0 is nr cycles to use):       10\n");
  fprintf(stderr, "  initial peak amplitude:                         1.0\n");
  fprintf(stderr, "  resample factor (integer):                      1\n");
  fprintf(stderr, "  horizontal range in plot (seconds):             3.0\n");
  fprintf(stderr, "  tsamp_orig (live: 0.000125):                    0.000125\n");
  fprintf(stderr, "  nr of traces stacked in second plot device:     1\n");
  fprintf(stderr, "  only plot when trace is complete (1=yes, 0=no): 0\n");
  fprintf(stderr, "  reset plot after new connection:                1\n");
  fprintf(stderr, "  delay of plot seconds (to align with sound):    1.2\n");
  fprintf(stderr, "  powerlaw index for sonification:                2\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  In the title, %%P means pulsar name as read from header if sent by sender.\n");
  fprintf(stderr, "  If the autopeak_flag > 0, the median of the maxima found during the specified nr of cycles will be used to determine the peak value. Note that having a fixed amplitude gives less fluctuation in the produced sound.\n");
}


char *skipUntilAfterCharacter(char *txt, char delimiter)
{
  char *txtptr;
  int i, found;
  txtptr = txt;
  found = 0;
  for(i = 0; i < strlen(txt)-1; i++) {
    if(txtptr[0] == delimiter) {
      found = 1;
      break;
    }
    txtptr += 1;
  }
  if(found)
    return txtptr+1;
  else
    return NULL;
}

void Readparameters(char *filename, float *baseline, float *peak, int *resample, int *autobaseline, int *autopeak, float *period, int *plot_length_samples, float *tsamp_orig, char *title, int *NrPlottedPulses, float *delay_plot_seconds, long *pulsenumber, int *plotOnePulseAtTheTime, float *sonification_pwrlaw, int *reset_plot_after_new_connection, int *autopeak_nrinringbuffer, int *autopeak_curbufferpos)
{
  FILE *fin;
  char line[1000];
  float tmp_baseline_value, tmp_peak_value, tmp_period, tmp_tsamp, tmp_delay, tmp_sonification_pwrlaw;
  int tmp_resample, tmp_autobaseline, tmp_autopeak, tmp6, tmp_nrpulses, tmp_oneatthetime, tmp_reset_after_connection;
  char tmp_title[100];
  int linenr, ret;
  char *txt_ptr;
  fin = fopen(filename, "r");
  if(fin == NULL) {
    fprintf(stderr, "Cannot open %s\n\n", filename);
  }else {
    linenr = 1;

    if(fgets(line, 1000, fin) == NULL) {
      fprintf(stderr, "Cannot read in line %d from %s\n", linenr, filename);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    txt_ptr = skipUntilAfterCharacter(line, ':');
    if(txt_ptr == NULL) {
      fprintf(stderr, "Cannot find a : in line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    ret = sscanf(txt_ptr, "%s", tmp_title);
    if(ret != 1) {
      fprintf(stderr, "Cannot find a string after : on line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    linenr++;

    if(fgets(line, 1000, fin) == NULL) {
      fprintf(stderr, "Cannot read in line %d from %s\n", linenr, filename);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    txt_ptr = skipUntilAfterCharacter(line, ':');
    if(txt_ptr == NULL) {
      fprintf(stderr, "Cannot find a : in line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    ret = sscanf(txt_ptr, "%d", &tmp_autobaseline);
    if(ret != 1) {
      fprintf(stderr, "Cannot find an integer after : on line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    linenr++;

    if(fgets(line, 1000, fin) == NULL) {
      fprintf(stderr, "Cannot read in line %d from %s\n", linenr, filename);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    txt_ptr = skipUntilAfterCharacter(line, ':');
    if(txt_ptr == NULL) {
      fprintf(stderr, "Cannot find a : in line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    ret = sscanf(txt_ptr, "%f", &tmp_baseline_value);
    if(ret != 1) {
      fprintf(stderr, "Cannot find a float after : on line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    linenr++;

    if(fgets(line, 1000, fin) == NULL) {
      fprintf(stderr, "Cannot read in line %d from %s\n", linenr, filename);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    txt_ptr = skipUntilAfterCharacter(line, ':');
    if(txt_ptr == NULL) {
      fprintf(stderr, "Cannot find a : in line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    ret = sscanf(txt_ptr, "%d", &tmp_autopeak);
    if(ret != 1) {
      fprintf(stderr, "Cannot find an integer after : on line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    linenr++;

    if(fgets(line, 1000, fin) == NULL) {
      fprintf(stderr, "Cannot read in line %d from %s\n", linenr, filename);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    txt_ptr = skipUntilAfterCharacter(line, ':');
    if(txt_ptr == NULL) {
      fprintf(stderr, "Cannot find a : in line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    ret = sscanf(txt_ptr, "%f", &tmp_peak_value);
    if(ret != 1) {
      fprintf(stderr, "Cannot find a float after : on line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    linenr++;

    if(fgets(line, 1000, fin) == NULL) {
      fprintf(stderr, "Cannot read in line %d from %s\n", linenr, filename);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    txt_ptr = skipUntilAfterCharacter(line, ':');
    if(txt_ptr == NULL) {
      fprintf(stderr, "Cannot find a : in line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    ret = sscanf(txt_ptr, "%d", &tmp_resample);
    if(ret != 1) {
      fprintf(stderr, "Cannot find an integer after : on line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    linenr++;

    if(fgets(line, 1000, fin) == NULL) {
      fprintf(stderr, "Cannot read in line %d from %s\n", linenr, filename);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    txt_ptr = skipUntilAfterCharacter(line, ':');
    if(txt_ptr == NULL) {
      fprintf(stderr, "Cannot find a : in line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    ret = sscanf(txt_ptr, "%f", &tmp_period);
    if(ret != 1) {
      fprintf(stderr, "Cannot find a float after : on line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    linenr++;

    if(fgets(line, 1000, fin) == NULL) {
      fprintf(stderr, "Cannot read in line %d from %s\n", linenr, filename);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    txt_ptr = skipUntilAfterCharacter(line, ':');
    if(txt_ptr == NULL) {
      fprintf(stderr, "Cannot find a : in line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    ret = sscanf(txt_ptr, "%f", &tmp_tsamp);
    if(ret != 1) {
      fprintf(stderr, "Cannot find a float after : on line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    linenr++;

    if(fgets(line, 1000, fin) == NULL) {
      fprintf(stderr, "Cannot read in line %d from %s\n", linenr, filename);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    txt_ptr = skipUntilAfterCharacter(line, ':');
    if(txt_ptr == NULL) {
      fprintf(stderr, "Cannot find a : in line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    ret = sscanf(txt_ptr, "%d", &tmp_nrpulses);
    if(ret != 1) {
      fprintf(stderr, "Cannot find an integer after : on line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    linenr++;

    if(fgets(line, 1000, fin) == NULL) {
      fprintf(stderr, "Cannot read in line %d from %s\n", linenr, filename);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    txt_ptr = skipUntilAfterCharacter(line, ':');
    if(txt_ptr == NULL) {
      fprintf(stderr, "Cannot find a : in line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    ret = sscanf(txt_ptr, "%d", &tmp_oneatthetime);
    if(ret != 1) {
      fprintf(stderr, "Cannot find an integer after : on line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    linenr++;

    if(fgets(line, 1000, fin) == NULL) {
      fprintf(stderr, "Cannot read in line %d from %s\n", linenr, filename);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    txt_ptr = skipUntilAfterCharacter(line, ':');
    if(txt_ptr == NULL) {
      fprintf(stderr, "Cannot find a : in line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    ret = sscanf(txt_ptr, "%d", &tmp_reset_after_connection);
    if(ret != 1) {
      fprintf(stderr, "Cannot find an integer after : on line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    linenr++;

    if(fgets(line, 1000, fin) == NULL) {
      fprintf(stderr, "Cannot read in line %d from %s\n", linenr, filename);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    txt_ptr = skipUntilAfterCharacter(line, ':');
    if(txt_ptr == NULL) {
      fprintf(stderr, "Cannot find a : in line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    ret = sscanf(txt_ptr, "%f", &tmp_delay);
    if(ret != 1) {
      fprintf(stderr, "Cannot find a float after : on line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    linenr++;

    if(fgets(line, 1000, fin) == NULL) {
      fprintf(stderr, "Cannot read in line %d from %s\n", linenr, filename);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    txt_ptr = skipUntilAfterCharacter(line, ':');
    if(txt_ptr == NULL) {
      fprintf(stderr, "Cannot find a : in line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    ret = sscanf(txt_ptr, "%f", &tmp_sonification_pwrlaw);
    if(ret != 1) {
      fprintf(stderr, "Cannot find a float after : on line %d from %s: line=%s\n", linenr, filename, line);
      showParametersExplanation();
      fclose(fin);
      return;
    }
    linenr++;



    tmp6 = tmp_period*8000;
    if(tmp6 > MaxPlotLength) {
      fprintf(stderr,"Requested plot period is too large, it is set to 1 second\n");
      tmp6 = 8000;
    }else if(tmp6 < 1000) {
      fprintf(stderr,"Requested plot period is too small, it is set to 1 second\n");
      tmp6 = 8000;
    }
    if(tmp_nrpulses > MaxNrTracesInPlot || tmp_nrpulses < 1) {
      fprintf(stderr,"Requested number of plotted pulses is too large, it is set to 10\n");
      tmp6 = 10;
    }
    if((*baseline != tmp_baseline_value && tmp_autobaseline == 0)|| (*peak != tmp_peak_value && tmp_autopeak == 0) || *resample != tmp_resample || *autobaseline != tmp_autobaseline || *autopeak != tmp_autopeak || *period != tmp_period || *plot_length_samples != tmp6 || *tsamp_orig != tmp_tsamp || strcmp(tmp_title, title) != 0 || (*NrPlottedPulses != tmp_nrpulses) || (*delay_plot_seconds != tmp_delay) || (*plotOnePulseAtTheTime != tmp_oneatthetime) || (*sonification_pwrlaw != tmp_sonification_pwrlaw) || (*reset_plot_after_new_connection != tmp_reset_after_connection)) {
      *baseline = tmp_baseline_value;
      *peak = tmp_peak_value;
      *resample = tmp_resample;
      *autobaseline = tmp_autobaseline;
      if(tmp_autopeak < *autopeak) {
	*autopeak_nrinringbuffer = 0;
	*autopeak_curbufferpos = 0;
      }
      *autopeak = tmp_autopeak;
      *period = tmp_period;
      *plot_length_samples = tmp6;
      *tsamp_orig = tmp_tsamp;
      strcpy(title, tmp_title);
      if(*NrPlottedPulses != tmp_nrpulses)
	*pulsenumber = 0;
      *NrPlottedPulses = tmp_nrpulses;
      *delay_plot_seconds = tmp_delay;
      *plotOnePulseAtTheTime = tmp_oneatthetime;
      *sonification_pwrlaw = tmp_sonification_pwrlaw;
      *reset_plot_after_new_connection = tmp_reset_after_connection;
      fprintf(stderr, "New parameters: baseline = %f, peak = %f, resample factor = %d, autobaseline = %d, autopeak = %d, period = %f, tsamp_orig = %f, plot_length_samples = %d, title = \"%s\", NrPlottedPulses = %d, delay = %f, plot one pulse at the time = %d, sonification powerlaw = %f, reset plot after new connection = %d\n", *baseline, *peak, *resample, *autobaseline, *autopeak, *period, *tsamp_orig, *plot_length_samples, title, *NrPlottedPulses, *delay_plot_seconds, *plotOnePulseAtTheTime, *sonification_pwrlaw, *reset_plot_after_new_connection);
    }
    fclose(fin);
  }
}

// I think these are all related to the incomming signal
void connectPulsarSignalSender()
{
  /* Create socket for server for PuMa II */

  fprintf(stderr, "Setup connection to the sender of pulsar data.\n");
  if ((serversock = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) {
    Die("  Failed to create socket");
  }else {
    fprintf(stderr, "  Created socked %d\n", serversock);
  }
  /* Construct the server sockaddr_in structure */
  memset(&audioserver, 0, sizeof(audioserver));       /* Clear struct */
  audioserver.sin_family = AF_INET;                   /* Internet/IP */
  audioserver.sin_addr.s_addr = htonl(INADDR_ANY);    /* Incoming addr */
  audioserver.sin_port = htons(portnr_puma2);         /* server port */
  fprintf(stderr, "  Binding socket using port %d\n", portnr_puma2);    
  do {
    if (bind(serversock, (struct sockaddr *) &audioserver, sizeof(audioserver)) < 0) {
      NotDie("  Failed to bind the server socket");
      sleep(1);
    }else {
      break;
    }
  }while(portnr_puma2 > 0);                           /* Infinite loop */

  /* Listen on the server socket */
  fprintf(stderr, "  listen socket for sender to make contact\n");    
  if (listen(serversock, MAXPENDING) < 0) {
    Die("  Failed to listen on server socket");
  }
  unsigned int clientlen = sizeof(puma2client);
  /* Wait for client connection */
  if ((clientsock =
       accept(serversock, (struct sockaddr *) &puma2client,
	      &clientlen)) < 0) {
    Die("  Failed to accept client connection");
  }
  fprintf(stderr, "  pulsar signal connected (%s)\n", inet_ntoa(puma2client.sin_addr));
}

/* mode: 0 = Show Audio signal
         1 = No Audio
*/
void plotSignalCheck(float *timeseries, float *timeseries_sonified, int NrBins, float baseline, float peak, float xviewportSize, float yviewportSize, int nosonic, char *title, char *pulsarname)
{
  int i, j;
  char title2[200];

  memset(title2, 0, 200);
  i = 0;
  for(j = 0; j < 100; j++) {
    if(title[j] == 0) {
      break;
    }else if(title[j] == '_') {
      title2[i++] = ' '; 
    }else if(title[j] == '%') {
      if(title[j+1] == 'P') {
	strcat(title2, pulsarname);
	i = strlen(title2)+1;
	j++;
      }else {
	title2[i++] = title[j]; 
      }
    }else {
      title2[i++] = title[j]; 
    }
  }

  cpgpage();
  cpgbbuf();
  cpgsci(1);
  cpgsvp(0.15, xviewportSize, 0.2, yviewportSize);
  cpgswin(0,NrBins,baseline - 0.5*peak, baseline + 1.5*peak);
  /*  cpgbox("bcnst",0.0,0,"bcnst",0.0,0);  */
  cpgsch(2);
  cpgslw(8);
  cpgaxis("n", 0, baseline - 0.5*peak, NrBins, baseline - 0.5*peak, 0, NrBins/8000.0, 0, 0, 0.5, 0.0, 0.3, 1.0, 0); 
  if(nosonic == 1)
    cpgaxis("",  0, baseline + 1.5*peak, NrBins, baseline + 1.5*peak, 0, NrBins/8000.0, 0, 0, 0, 0.5, 0.3, -1.0, 0); 
  cpgsch(0.7);
  cpgaxis("", NrBins, baseline - 0.5*peak, NrBins, baseline + 1.5*peak, baseline - 0.5*peak, baseline - 1.5*peak, 0, 0, 0.5, 0, 0.5, -1.0, 0); 
  cpgaxis("", 0, baseline - 0.5*peak, 0, baseline + 1.5*peak, baseline - 0.5*peak, baseline + 1.5*peak, 0, 0, 0, 0.5, 0.5, -1.0, 0); 
  cpgslw(1);
  cpgaxis("n", 0, baseline - 0.5*peak, 0, baseline + 1.5*peak, baseline - 0.5*peak, baseline + 1.5*peak, 0, 0, 0, 0.5, 0.5, -1.0, 0); 
  cpgsch(2);
  cpgslw(8);
  cpglab("Time (sec)", "Radio Intensity", title2);
  cpgsch(1);
  cpgslw(1);

  cpgmove(0, timeseries[0]);
  for(j = 0; j < NrBins; j++) {
    cpgdraw(j, timeseries[j]);
  }


  cpgsci(2);
  cpgsls(4);
  cpgmove(0, baseline);
  cpgdraw(NrBins, baseline);
  cpgmove(0, baseline+peak);
  cpgdraw(NrBins, baseline+peak);
  cpgmove(0, baseline-peak);
  cpgdraw(NrBins, baseline-peak);
  cpgsls(1);

  if(nosonic == 0) {

    cpgsci(1);
    cpgsvp(0.15, xviewportSize, yviewportSize, 0.95);
    cpgswin(0,NrBins,0, 1);
    cpgbox("bcst",0.0,0,"bcnst",0.0,0);

    cpgsci(3);
    cpgmove(0, timeseries_sonified[0]); //*peak+baseline
    for(j = 0; j < NrBins; j++) {
      cpgdraw(j, timeseries_sonified[j]);
    }
  }
  cpgebuf();
}


void plotPulsarSignal_reset(int NrBins, int NrPulses, float xviewportSize, float yviewportSize, char *title, char *pulsarname)
{
  int i, j;
  char title2[200];

  memset(title2, 0, 200);
  i = 0;
  for(j = 0; j < 100; j++) {
    if(title[j] == 0) {
      break;
    }else if(title[j] == '_') {
      title2[i++] = ' '; 
    }else if(title[j] == '%') {
      if(title[j+1] == 'P') {
	strcat(title2, pulsarname);
	i = strlen(title2)+1;
	j++;
      }else {
	title2[i++] = title[j]; 
      }
    }else {
      title2[i++] = title[j]; 
    }
  }
  
  cpgpage(); 
  cpgsci(1);
  cpgsvp(0.15, xviewportSize, 0.2, yviewportSize);
  cpgswin(0,NrBins,-0.5, NrPulses + 0.5);
  cpgsch(2);
  cpgslw(8);
  cpgaxis("n", 0, -0.5, NrBins, -0.5, 0, NrBins/8000.0, 0, 0, 0.5, 0.0, 0.3, 1.0, 0); 
  cpgaxis("",  0, NrPulses + 0.5, NrBins, NrPulses + 0.5, 0, NrBins/8000.0, 0, 0, 0, 0.5, 0.3, -1.0, 0); 
  cpgsch(0.7);
  cpgaxis("", NrBins, -0.5, NrBins, NrPulses + 0.5, -0.5, NrPulses+0.5, 0, 0, 0.5, 0, 0.5, -1.0, 0); 
  cpgaxis("", 0, -0.5, 0, NrPulses+ 0.5, -0.5, NrPulses+ 0.5, 0, 0, 0, 0.5, 0.5, -1.0, 0); 
  cpgslw(1);
  cpgsch(2);
  cpgslw(8);
  if(NrPulses == 1)
    cpgaxis("n", 0, -0.5, 0, NrPulses +0.5, -0.5, NrPulses+ 0.5, 0, 0, 0, 0.5, 0.5, -1.0, 0); 
  else
    cpgaxis("n", 0, -0.5, 0, NrPulses +0.5, -0.5+1, NrPulses+ 0.5+1, 0, 0, 0, 0.5, 0.5, -1.0, 0); 
  if(NrPulses == 1)
    cpglab("Time (sec)", "Intensity", title2);
  else
    cpglab("Time (sec)", "Pulse number", title2);
  cpgsch(1);
  cpgslw(1);
}

void plotPulsarSignal(float *Iserie, int pulse, int NrBins, int NrPulses, float baseline, float peak, float xviewportSize, float yviewportSize, char *title, int startbin, int endbin)
{
  int j;
  
  //  fprintf(stderr, "XXXX: %d - %d\n", startbin, endbin);

  cpgmove(startbin, pulse + (Iserie[startbin]-baseline)/peak);
  for(j = startbin; j <= endbin; j++) {
    cpgdraw(j, pulse + (Iserie[j]-baseline)/peak);
  }
}

float calc_median(float *buffer, int nrentries)
{
  int i, j, value, best_i, best_value;
  int tmp_buffer[MaxAutopeakValue];

  if(nrentries <= 0)
    return 0;

  for(i = 0; i < nrentries; i++) {
    value = 0;
    for(j = 0; j < nrentries; j++) {
      if(buffer[i] < buffer[j])
	value++;
    }
    tmp_buffer[i] = value;
  }
  for(i = 0; i < nrentries; i++) {
    tmp_buffer[i] = abs(tmp_buffer[i] - 0.5*nrentries);
  }
  
  best_i = 0; 
  best_value = tmp_buffer[0];
  for(i = 0; i < nrentries; i++) {
    if(tmp_buffer[i] < best_value) {
      best_i = i;
      best_value = tmp_buffer[i];
    }
  }

  return buffer[best_i];
}
