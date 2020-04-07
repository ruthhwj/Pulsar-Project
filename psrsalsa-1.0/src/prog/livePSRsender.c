#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netdb.h>         
#include <arpa/inet.h>     
#include <netinet/in.h>   
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <math.h>
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
  long maxAllowedOversampling;   // The input data may have up to this factor extra floats per package
  long nrFreqChan;               // Default is 1 frequency channel being transmitted
}packet_def;

struct {
  int socketNr;                      // Socket for sending the data to audio server
  unsigned int portnr;               // port number of audio server
  struct sockaddr_in server;         // The server to connect with
  char hostname[1000];                // Name computer audio server is running on
  char ip_host[128];                 // IP of audio server
  int stopSendingData_flag;          // Flag used if connection fails
  FILE *output_file;                 // If sending to a file rather than an internet connection
}connection_def;

struct {
  int byteswap;                      // Byte swap needs to be applied to data being read in
  double start_skip_time;            // The requested amount of time to be skipped
  long nrInputSamplesPerPackage;     // Nr of samples than need to be read in from the input data that will form a package. This will be resampled to be nrSamplesPerPackage
  char psrname[1000];
  double samplingtime;
  int fileformat;
  long fileHeaderSize;
  FILE *fin;                         // File pointer (used for some formats)
  datafile_definition datafile;      // Header (used for other formats)
}input_data;


void Die(char *mess) { perror(mess); exit(1); }     /* Let program die, and say why */
void NotDie(char *mess) { perror(mess); }           /* Stay alive, but say what failed */
void openFile(char *filename);
void SentFile(int sentheader);
void setSamplingTime(double samplingtime);
void swapArray_sender(void *xptr, int nelem, int size);



/* defenition nanosleep function, sometimes needed, sometimes you have
   to leave it out. */
/* int nanosleep(const timespec *req, timespec *rem);   */

int main(int argc, char **argv)
{
  int i, ret, sentheader;
  verbose_definition verbose;

  packet_def.packetDuration = 50;
  packet_def.nrSamplesPerPackage = 400;
  packet_def.maxAllowedOversampling = 10;
  packet_def.nrFreqChan = 1;


  cleanVerboseState(&verbose);
  verbose.verbose = 1;
  connection_def.output_file = NULL;  // Default: sent over network

  connection_def.portnr = 3000;            /* Default port number of audio server */
  connection_def.ip_host[0] = 0;           /* No default ip of audio server */
  input_data.fileHeaderSize = 0;       // Default is headerless data
  input_data.fileformat = 0;           // Default is raw floats
  input_data.psrname[0] = 0;
  sentheader = 0;

  /* If there are no options given to program, show help */
  if(argc < 2) {
    printf("Usage: livePSRsender [options] file.\n\nThe input file will be read and the sample-rate (should be <= 8000 samples per second, if a non-integer nr of samples span %.0lf ms, the data rate will not be entirely accurate) will be sent over the network at a fixed rate of 8000 floating-point samples per second. The program livePSRaudioserver can receive this data stream.\n\n", packet_def.packetDuration);
    printf("-h         \"host_name\" Set the host name of computer where livePSRaudioserver is running. The signal will\n");
    printf("           be sent to this computer to the port specified with the -p option.\n");
    printf("-i         alternatively, set its IP adress.\n");
    printf("-p         portnumber to use to connect to audio server. Default is %d.\n", connection_def.portnr);
    printf("-file      Write out stream to this file instead\n");
    printf("-byteswap  The floating points of the input file are byte-swapped\n");
    printf("-bskip     Skip this amount of seconds in bytes file (can be used together with -tskip)\n");
    printf("-tskip     Skip this amount of seconds in input file (can be used together with -bskip)\n");
    printf("-tsamp     Set sampling time of the input data (in sec). Default is %f, unless -puma1 or -sigproc is specifies.\n", 0.001*packet_def.packetDuration/(float)packet_def.nrSamplesPerPackage);
    printf("-psrname   Set the source name of the pulsar\n");
    printf("-header    Sent a package with header information before sending the actual data.\n");
    printf("-puma1     Data file is puma1 format. Implies -bskip 4504 and -byteswap and the sampling time/source name will be determined from the file.\n");
    printf("-sigproc   Data file is sigproc format. Implies -bskip 316 and the sampling time/source name will be determined from the file.\n");
    return 0;
  }
  
  /* Parse command line */
  input_data.byteswap = 0;
  input_data.nrInputSamplesPerPackage = packet_def.nrSamplesPerPackage;          // Default. Gets adjusted when sampling time of data is known
  input_data.start_skip_time = 0;
  for(i = 1; i < argc-1; i++) {
    if(strcmp(argv[i], "-h") == 0) {
      strcpy(connection_def.hostname, argv[i+1]);
      fprintf(stderr, "Try to find IP of the specified audio server: %s\n", connection_def.hostname);
      if(getIPofHost(connection_def.hostname, connection_def.ip_host, stderr, verbose) == 0) {
	fprintf(stderr, "Cannot determine IP address of the stream server, specify this with the -i option.");
      }
      i++;
   }else if(strcmp(argv[i], "-header") == 0) {
      sentheader = 1;
   }else if(strcmp(argv[i], "-i") == 0) {
      strcpy(connection_def.ip_host, argv[i+1]);      
      i++;
    }else if(strcmp(argv[i], "-file") == 0) {
      connection_def.output_file = fopen(argv[i+1], "w");
      if(connection_def.output_file == NULL) {
	fprintf(stderr, "Cannot open output file.");
	return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-p") == 0) {
      ret = sscanf(argv[i+1], "%d", &connection_def.portnr);
      if(ret != 1) {
	fprintf(stderr, "Cannot read in command line option %d.", connection_def.portnr);
	return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-puma1") == 0) {
      input_data.fileHeaderSize = 4504;
      input_data.byteswap = 1;
      input_data.fileformat = 1;
    }else if(strcmp(argv[i], "-sigproc") == 0) {
      input_data.fileHeaderSize = 316;
      input_data.byteswap = 0;
      input_data.fileformat = 2;
    }else if(strcmp(argv[i], "-bskip") == 0) {
      input_data.fileHeaderSize = atoi(argv[i+1]);      
      i++;
    }else if(strcmp(argv[i], "-tskip") == 0) {
      input_data.start_skip_time = atof(argv[i+1]);
      i++;
    }else if(strcmp(argv[i], "-psrname") == 0) {
      ret = sscanf(argv[i+1], "%s", input_data.psrname);
      if(ret != 1) {
	fprintf(stderr, "Cannot parse option %s: need one string\n", argv[i]);
	return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-tsamp") == 0) {
      setSamplingTime(atof(argv[i+1]));
      i++;
    }else if(strcmp(argv[i], "-byteswap") == 0) {
      input_data.byteswap = 1;
    }else {
      printf("\nUnknown option: %s\n", argv[i]);
      return 0;
    }
  }

  /* Check if the server is specified */
  if(connection_def.ip_host[0] == 0 && connection_def.output_file == NULL) {                                     
    printf("Specify the computer on which livePSRaudioserver is running (using -h or -i).\n");
    return 0;
  }

  openFile(argv[argc-1]);
  if(connection_def.output_file == NULL) {
    printf("Create socket for connection with audio server.\n");
    if ((connection_def.socketNr = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) {
      Die("Failed to create socket");
    }else {
      fprintf(stdout, "  Created socked %d\n", connection_def.socketNr);
    }
    memset(&connection_def.server, 0, sizeof(connection_def.server));                 /* Clear struct */
    connection_def.server.sin_family = AF_INET;                             /* Internet/IP */
    connection_def.server.sin_addr.s_addr = inet_addr(connection_def.ip_host);             /* IP address */
    connection_def.server.sin_port = htons(connection_def.portnr);                         /* connection_def.server port */
    fprintf(stdout, "Trying to establish connection to %s at port %d\n", connection_def.ip_host, connection_def.portnr);    
    if (connect(connection_def.socketNr, (struct sockaddr *) &connection_def.server, sizeof(connection_def.server)) < 0) {
      NotDie("  Failed to connect with server, but will keep trying");
      connection_def.stopSendingData_flag = 3;      /* Will not start with sending data, but keep trying to connect */
    }else {
      fprintf(stdout, "  Established connection\n\n");    
      connection_def.stopSendingData_flag = 0;      /* Default is an established connection to audio server */
    }
  }else {
    connection_def.stopSendingData_flag = 0;
  }

  printf("Going to send file (%s)\n", argv[argc-1]);
  SentFile(sentheader);

  if(connection_def.output_file == NULL) {
    close(connection_def.socketNr);              /* Close connection */
  }
  return 0;
}


// When opening a sigproc file rather than an other format, datafile is used rather than fptr
void openFile(char *filename)
{
  verbose_definition noverbose, verbose;
  cleanVerboseState(&noverbose);
  cleanVerboseState(&verbose);
  verbose.verbose = 1;

  input_data.fin = fopen(filename, "rb");
  if((input_data.fin) == NULL) {
    printf("Cannot open %s\n\n", filename);
    exit(0);
  }
  if(input_data.fileformat == 1) {   // Got a puma1 file
    /*   
    //Test to find out where sampling time is stored in header
    Header_type puma_hdr;
    void *position, *position2;
    position = &(puma_hdr.redn.DeltaTime);
    position2 = &puma_hdr;
    printf("position of sampling time is %p = %d\n", position-position2, position-position2);
    */
    //Test to find out where sampling time is stored in header
    printf("Get sampling time from puma file:\n");
    fseek(input_data.fin, 4040, SEEK_SET);
    fread(&input_data.samplingtime, 8, 1, input_data.fin);
    swapArray_sender(&input_data.samplingtime, 1, 8);
    setSamplingTime(input_data.samplingtime);

    fseek(input_data.fin, 776, SEEK_SET);
    fread(input_data.psrname, 24, 1, input_data.fin);
    printf("Source name = '%s'\n", input_data.psrname);

    printf("\n");

  }else if(input_data.fileformat == 2) {   // Got a sigproc file

    
    if(openPSRData(&input_data.datafile, filename, 0, 0, 0, 0, verbose) == 0) {
      printf("Cannot open %s\n\n", filename);
      exit(0);
    }
    if(readHeaderPSRData(&input_data.datafile, 1, 0, verbose) != 1) {
      printf("Cannot read header from %s\n\n", filename);
      exit(0);
    }
    input_data.fileHeaderSize = input_data.datafile.datastart;

    /*
    fread(&header, 316, 1, *fptr);
    printf("Get sampling time from sigproc file:\n");
    char header[317], substring[6];
    char *hdr_ptr;
    double *samplingtime;
    strcpy(substring, "tsamp");
    hdr_ptr = searchStringInMem(header, 316, substring, strlen(substring), noverbose);
    if(hdr_ptr == NULL) {
      fprintf(stderr, "Cannot find sampling time in sigproc file. Do not use the -sigproc option and use the -tsamp option.\n");
      exit(0);
    }
    hdr_ptr += strlen(substring);
    samplingtime = (double *)hdr_ptr;
    setSamplingTime(*samplingtime, &nrInputFloatsPerPackage, packagesize);
    */
    input_data.samplingtime = get_tsamp(input_data.datafile, 0, verbose);
    setSamplingTime(input_data.samplingtime);
    
    /*
    strcpy(substring, "source_name");
    hdr_ptr = searchStringInMem(header, 316, substring, strlen(substring), noverbose);
    if(hdr_ptr == NULL) {
      fprintf(stderr, "Cannot find source name in sigproc file. Do not use the -sigproc option and use the -XXXXX option.\n");
      exit(0);
    }
    hdr_ptr += strlen(substring);
    for(i = 0; i < 13; i++) {  //Skip some 0 and backspace characters which are there for some reason
      if(hdr_ptr[0] < 32)
	hdr_ptr++;
    }
    for(i = 0; i < 13; i++) {  //Make the first non normal character a 0 (i.e. remove line feed etc)
      if(hdr_ptr[i] < 32 || hdr_ptr[i] > 126) {
	hdr_ptr[i] = 0;
	break;
      }
    }
    sscanf(hdr_ptr, "%s", psrname);
    printf("Source name = '%s'\n", psrname);
    */

    strcpy(input_data.psrname, input_data.datafile.psrname);
    printf("Source name = '%s'\n", input_data.psrname);
    printf("\n");

    packet_def.nrFreqChan = input_data.datafile.NrFreqChan;

    /*
    printf("%d     '%c'\n", hdr_ptr[0], hdr_ptr[0]);
    printf("%d     '%c'\n", hdr_ptr[1], hdr_ptr[1]);
    printf("%d     '%c'\n", hdr_ptr[2], hdr_ptr[2]);
    printf("%d     '%c'\n", hdr_ptr[3], hdr_ptr[3]);
    printf("%d     '%c'\n", hdr_ptr[4], hdr_ptr[4]);
    printf("%d     '%c'\n", hdr_ptr[5], hdr_ptr[5]);
    printf("%d     '%c'\n", hdr_ptr[6], hdr_ptr[6]);
    printf("%d     '%c'\n", hdr_ptr[7], hdr_ptr[7]);
    printf("%d     '%c'\n", hdr_ptr[8], hdr_ptr[8]);
    printf("%d     '%c'\n", hdr_ptr[9], hdr_ptr[9]);
    printf("%d     '%c'\n", hdr_ptr[10], hdr_ptr[10]);
    printf("%d     '%c'\n", hdr_ptr[11], hdr_ptr[11]);
    printf("%d     '%c'\n", hdr_ptr[12], hdr_ptr[12]);
    printf("%d     '%c'\n", hdr_ptr[13], hdr_ptr[13]);
    printf("%d     '%c'\n", hdr_ptr[14], hdr_ptr[14]);
    printf("%d     '%c'\n", hdr_ptr[15], hdr_ptr[15]);
    printf("%d     '%c'\n", hdr_ptr[16], hdr_ptr[16]);
    printf("%d     '%c'\n", hdr_ptr[17], hdr_ptr[17]);
    printf("%d     '%c'\n", hdr_ptr[18], hdr_ptr[18]);
    printf("%d     '%c'\n", hdr_ptr[19], hdr_ptr[19]);
    */  
  }
  fseek(input_data.fin, input_data.fileHeaderSize+input_data.start_skip_time*1000*4.0*input_data.nrInputSamplesPerPackage/packet_def.packetDuration, SEEK_SET);               /* Skip header of file*/
}

void printcodes()
{
  printf("Codes are:\n  '.' send ok\n  'S' try open socket\n  'C' try reconnect audio server\n");
}

void SentFile(int sentheader)
{
  float *bufferfloats,*bufferfloats_in;    // bufferfloats_in will be filled, and resampled in bufferfloats and sent
  int received;  // Nr of time samples (can have multiple freq channels) received
  int firsttime, nrdotsonscreen;
  int output_error;                  /* Flag used to stop repeated error messages */
  struct timespec req_delay;               /* Variable to set delay of nanosleep */
  int i, j;
  //  long k, readnumber;
  double expected_timedif, timediff;
  struct timeval cur_time, start_time;
  verbose_definition verbose;
  cleanVerboseState(&verbose);

  // Same memory as float buffer, but of different type
  char* buffer_in = malloc(packet_def.maxAllowedOversampling*sizeof(float)*packet_def.nrSamplesPerPackage*packet_def.nrFreqChan);
  char* buffer = malloc(sizeof(float)*packet_def.nrSamplesPerPackage*packet_def.nrFreqChan);
  bufferfloats = (float *)buffer;
  bufferfloats_in = (float *)buffer_in;
  if(buffer == NULL || buffer_in == NULL) {
    fprintf(stderr, "Cannot allocate buffer size\n");
    exit(0);
  }




  output_error = 1;                              /* Output errors */
  //  k = 0;
  expected_timedif = 0;
  firsttime = 1;
  nrdotsonscreen = 0;
  //  readnumber = 0;  // Used in sigproc format, the current read number.

  gettimeofday(&start_time, NULL);


  do {
    fflush(stdout);                             /* Make sure the printf's are showed up */

    /* Read a block of floats. Decide if it should be a header or data. */
    if(firsttime && sentheader) {   // Fill a complete package with the header
      char *hdr_ptr;
      char txt[100];
      hdr_ptr = buffer_in;
      received = input_data.nrInputSamplesPerPackage;  // The header size will be a full package
      memset(buffer_in, 0, received*sizeof(float)*packet_def.nrFreqChan);

      strcpy(txt, "livePSRsender header");
      strcpy(hdr_ptr, txt);
      hdr_ptr += strlen(txt)+1;

      if(strlen(input_data.psrname) == 0)        //Avoid empty name, as makes less easy to read header.
	strcpy(input_data.psrname, " ");
      strcpy(hdr_ptr, input_data.psrname);
      hdr_ptr += strlen(input_data.psrname)+1;

      sprintf(txt, "%ld", packet_def.nrFreqChan);
      strcpy(hdr_ptr, txt);
      hdr_ptr += strlen(txt)+1;

    }else {
      if(input_data.fileformat == 2) {
	/*
	if(packet_def.nrFreqChan == 1) {
	  if(readPulsePSRData(&input_data.datafile, 0, 0, 0, readnumber*input_data.nrInputSamplesPerPackage, input_data.nrInputSamplesPerPackage, bufferfloats_in, verbose) != 1) {
	    received = 0;
	  }else {
	    received = input_data.nrInputSamplesPerPackage;
	  }
	  readnumber ++;
	}else {
	  long freqnr, binnr;
	  received = 1;
	  for(binnr = readnumber*input_data.nrInputSamplesPerPackage; binnr < (readnumber+1)*input_data.nrInputSamplesPerPackage; binnr++) {
	    if(received == 1) {
	      for(freqnr = 0; freqnr < packet_def.nrFreqChan; freqnr++) {
		if(readPulsePSRData(&input_data.datafile, 0, 0, freqnr, binnr, 1, &(bufferfloats_in[(binnr-readnumber*input_data.nrInputSamplesPerPackage)*packet_def.nrFreqChan+freqnr]), verbose) != 1) {
		  received = 0;
		  break;
		}else {
		  received = 1;
		}
	      }
	    }
	  }
	  if(received == 1)
	    received = input_data.nrInputSamplesPerPackage;
	  readnumber ++;
	}
	*/
	// Read all frequency channel data for the given number of input samples
	received = fread(buffer_in, 4, input_data.nrInputSamplesPerPackage*packet_def.nrFreqChan, input_data.fin);
	received /= packet_def.nrFreqChan;  // Received is the number of time samples
      }else {
	received = fread(buffer_in, 4, input_data.nrInputSamplesPerPackage, input_data.fin);  
	if(input_data.byteswap) 
	  swapArray_sender(buffer_in, input_data.nrInputSamplesPerPackage, 4);
      }
    }

    expected_timedif += 0.001*packet_def.packetDuration;
    gettimeofday(&cur_time, NULL);
    timediff = expected_timedif - (cur_time.tv_sec + 0.000001*cur_time.tv_usec - (start_time.tv_sec + 0.000001*start_time.tv_usec)); 
    /*      printf("%lf\n", timediff);  */
    if(timediff > 0.1) {
      if(timediff > 0.2)
	printf("Glitch: %f seconds\n", timediff-0.1);
      req_delay.tv_sec = 0; 
      req_delay.tv_nsec = (timediff-0.1)*1000000000;
      /*	printf("\nSleep for %f sec to keep synchronized.\n", timediff-0.1); */
      nanosleep(&req_delay, NULL);               
    }
    if(received == input_data.nrInputSamplesPerPackage) {                          /* Check if read was succesfull */
      if(input_data.nrInputSamplesPerPackage <= packet_def.nrSamplesPerPackage) {
	if(firsttime && sentheader) {                            // Do not resample if it is a header package
	  if(packet_def.nrSamplesPerPackage >= input_data.nrInputSamplesPerPackage) {
	    memset(buffer, 0, sizeof(float)*packet_def.nrSamplesPerPackage*packet_def.nrFreqChan);  // Set the full packet to be sent to zero
	    memcpy(buffer, buffer_in, input_data.nrInputSamplesPerPackage*sizeof(float)*packet_def.nrFreqChan);  // The copy is smaller because of resampling that normally can oversamples the data
	  }else {
	    memcpy(buffer, buffer_in, packet_def.nrSamplesPerPackage*sizeof(float)*packet_def.nrFreqChan);  // Don't copy more than can fit in the package
	  }
      //	  memcpy(bufferfloats, bufferfloats_in, sizeof(float)*packet_def.nrSamplesPerPackage*packet_def.nrFreqChan);  
	}else {
	  for(i = 0; i < packet_def.nrSamplesPerPackage; i++) {
	    j = i*(input_data.nrInputSamplesPerPackage-1.0)/(float)(packet_def.nrSamplesPerPackage-1.0);
	    if(packet_def.nrFreqChan == 1) {
	      bufferfloats[i] = bufferfloats_in[j];
	    }else {
	      long freqnr;
	      for(freqnr = 0; freqnr < packet_def.nrFreqChan; freqnr++) {
		bufferfloats[i*packet_def.nrFreqChan+freqnr] = bufferfloats_in[j*packet_def.nrFreqChan+freqnr];		
	      }
	    }
	  }
	}
      }else {
	printf("Cannot handle data with higher time-resolution that the output resolution (%ld >= %ld)!\n", input_data.nrInputSamplesPerPackage, packet_def.nrSamplesPerPackage);
	exit(0);
      }
      received = sizeof(float)*packet_def.nrSamplesPerPackage;

      if(firsttime) {
	printf("First byte values sent are: ");
	for(i = 0; i < 20; i++) {
	  //	  printf("%d '%c'", buffer[i], buffer[i]);
	  printf("%d ", buffer_in[i]);
	}
	printf("\n\n");
	printcodes();
	firsttime = 0;
      }

      if(connection_def.stopSendingData_flag == 0 || connection_def.output_file != NULL) {                        /* Check if connected */
	if(connection_def.output_file == NULL) {
	  if(send(connection_def.socketNr, buffer, received*packet_def.nrFreqChan, MSG_NOSIGNAL) != received*packet_def.nrFreqChan) {
	    NotDie("\nFailed to send bytes to client");
	    received = 4*packet_def.nrSamplesPerPackage;
	    connection_def.stopSendingData_flag = 1;                       /* Connection failed, try to reconnect */
	  }else {
	    printf(".");
	    nrdotsonscreen++;
	    if(nrdotsonscreen == 2000) {
	      printf("\n");
	      printcodes();
	      nrdotsonscreen = 0;
	    }
	  }
	}else {
	  fwrite(buffer, 1, received*packet_def.nrFreqChan, connection_def.output_file);
	}
      }else if(connection_def.stopSendingData_flag == 1) {          /* Connection failed, so close socked */
	close(connection_def.socketNr);
	connection_def.stopSendingData_flag++;
      }else if(connection_def.stopSendingData_flag == 2) {        /* Try to open new socket */
	printf("S");
	if ((connection_def.socketNr = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) {
	  NotDie("\nFailed to create socket");
	}else {
	  fprintf(stdout, "\nCreated socked %d\n", connection_def.socketNr);
	  memset(&connection_def.server, 0, sizeof(connection_def.server));       /* Clear struct */
	  connection_def.server.sin_family = AF_INET;                  /* Internet/IP */
	  connection_def.server.sin_addr.s_addr = inet_addr(connection_def.ip_host);  /* IP address */
	  connection_def.server.sin_port = htons(connection_def.portnr);       /* server port */
	  connection_def.stopSendingData_flag++;
	}
      }else if(connection_def.stopSendingData_flag == 3) {       /* Try to reconnect to audio server */
	printf("C");
	if (connect(connection_def.socketNr, (struct sockaddr *) &connection_def.server, sizeof(connection_def.server)) < 0) {
	  if(output_error == 1)
	    NotDie("Failed to connect with server");
	  output_error = 0;
	}else {
	  fprintf(stdout, "\nConnected\n"); 
	  connection_def.stopSendingData_flag = 0;
	  output_error = 1;
	}
      }
    }
  }while(received == 4*packet_def.nrSamplesPerPackage);
  if(input_data.fileformat == 2) {
    closePSRData(&input_data.datafile, 0, verbose);
  }else {
    fclose(input_data.fin);
  }
  printf("\nEnd of file\n");

  free(buffer_in);
  free(buffer);
}


/* redwards func to swap an arbitrary word */
/* swaps from the out-in */
void swapWord_sender(void *xptr, int size)
{
  unsigned char *cptr = (unsigned char *)xptr;
  unsigned char tmp;
  int i;

  for (i=0; i < size/2; i++)
  {
    tmp = cptr[i];
    cptr[i] = cptr[size-1-i];
    cptr[size-1-i] = tmp;
  }
}

/* redwards func to swap an array */
void swapArray_sender(void *xptr, int nelem, int size)
{
  int i;
  unsigned char *cptr = (unsigned char *)xptr;
  for (i=0; i < nelem; i++)
  {
    swapWord_sender(cptr, size);
    cptr += size;
  }
}

void setSamplingTime(double samplingtime)
{
  input_data.nrInputSamplesPerPackage = 0.001*packet_def.packetDuration/samplingtime;
  printf("Sampling time is %lf samples per second.\n", samplingtime);
  printf("%ld samples in the input data corresponds to %.0lf ms.\n", input_data.nrInputSamplesPerPackage, packet_def.packetDuration);
  if(input_data.nrInputSamplesPerPackage > packet_def.maxAllowedOversampling*4*packet_def.nrSamplesPerPackage) {
    printf("Please reduce resolution of the input data!\n");
    exit(0);
  }
}
