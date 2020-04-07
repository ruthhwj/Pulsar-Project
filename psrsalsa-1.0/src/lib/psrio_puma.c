//START REGION RELEASE
#define _FILE_OFFSET_BITS 64
#define _USE_LARGEFILE 1
#define _LARGEFILE_SOURCE 1

#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
//START REGION DEVELOP
//#include "pumadata.h"
//#include "bigendian.h"
//#include "libpuma.h"
//#include "fpuma.h"
//#include "sla_wrap.h"
//START REGION RELEASE
#include "psrsalsa.h"


//#ifdef __cplusplus
//#define PUMA_STD std::
//#include <ctime>
//#include <cstdio>
//#else
#define PUMA_STD
//#include <time.h>
//#include <stdio.h>
//#endif

/* $Header: /home/bws/Repository/src/PuMa/libs/puma/puma.h,v 1.6 2004/01/05 13:41:17 redwards Exp $ */

#ifndef _PUMA_H  /* this whole header file */
#define _PUMA_H


/* redwards mods to make this file compliant with C++ standard */


/* puma.h - header file for TMS-PuMa shared memory areas */

/* pm 18 Sep 98 - JLLV */


/*                |========================================|               *
 * ===============|Shared memory header file per PuMa crate|============== *
 *                |========================================|               *
 *                                                                         *
 * This header file defines the shared memory via which each PuMa crate    *
 * can communicate with TMS.                                               *
 * The shared memory consists of a TMS area, which is considered read-only *
 * for PuMa, and similarly a PuMa-crate area which is read-only for TMS.   *
 * The variables of these areas are defined at the end of this file, the   *
 * rest contains the definitions of the various types and constants.       *
 *                                                                         *
 * ======================================================================= */


/*                                |=====|                                  *
 * ===============================| N.B.|================================= *
 *                                |=====|                                  *
 *                                                                         *
 * For the SHARC software, only 'MAXPACKS is required and "time.h" is non  *
 * existing, so the SHARC Makefile defines TIME to skip most of this file: *
 *                                                                         *
 * ======================================================================= */

#define MMAX(A, B) (A > B? A: B)  /* local (re-)definition of MAX */

#define MAXDMS          4 /* The number of different DMs that mode 2 will handle */
#define MAXSTOKESPARAMS 4 /* That's the way it is . . . */

#define MAXPACKS        (MMAX(MAXDMS, MAXSTOKESPARAMS)) /* for adjustment tokens */

/* 'packcntrl' indexing for offset and scale arrays. For 4 DMs use 0..3.
   Also used by SHARCs (dd) */

#define PACKX           0
#define PACKY           1

#define PACKI           0
#define PACKQ           1
#define PACKU           2
#define PACKV           3

/* Mode 1, 2 and 3 specific - expressed in # of points: */

#define MAXFFTSIZE   8192
#define MINFFTSIZE     32

/* To remember which modeX programme ran last time: */

#define NoMODE         -9
#define MODEM1         -1
#define MODE0           0
#define MODE1           1
#define MODE2           2
#define MODE3           3
#define MODE4           4


#ifndef TIME  /* ============ skip till the end of this include file ============= */



/* BWS change 2003-04-06 as not needed! #define HEADERFORMATID "DPC_1.0" Dutch Pulsar Community*/

#include <time.h>


/* -------------------------- C O N S T A N T S --------------------------- */


/* P A R A M E T E R S */

/* General parameters: */

#ifndef Boolean
#define Boolean int  /* defined in "utils.h" as 'boolean' - vuschil mot ur wese */
#endif
#ifndef FALSE          /* also defined in "utils.h" */
#define FALSE   0
#endif
#ifndef TRUE           /* also defined in "utils.h" */
#define TRUE    1
#endif

#define NOREQUESTRECEIVED 0
#define WORKINGONREQUEST  1
#define IGNORE            "_IgNoRe_"  /* for as yet unknown RCS headers */

#define TXTLEN      64 /* Standard length of most TMS text fields */
#define NAMELEN     16 /* For observatory and unique observation ID */
#define COMMENTLEN 256
#define FILELEN     64 /* Maximum length of a file name */
#define CMDLEN     256 /* Maximum length of a command line */
#define DIRLEN     256 /* Maximum length of a full (sub-)dir name */
#define PATHLEN    256 /* Maximum length of a full path name */

/* Overall PuMa parameters: */

#define NCRATES              2 /* Only 2? */
#define CRATE0               0
#define CRATE1               1

#define MODE0REP    "mode0rep" /* 'mode0.c'    - reporting */
#define MODE1REP    "mode1rep" /* 'mode1.c'    - reporting */
#define ADJUST      "adjust"   /* 'puma2obs.c' - adjustment scratch file */
#define ARCREP      "arcrep"   /* 'puma2obs.c' - reporting */
#define DELREP      "delrep"   /* 'deleteobs.c'- reporting */
#define BINDAT      "bindat"   /* 'modes.c' & 'puma2obs.c' - file with all obs params */
#define OBSFNS      "obsfns"   /* 'thread.c'   - sizes + names of all PuMa obs files */
#define ASTRO       "puma"     /* 'puma2obs.c' - the tape files, geen hoofdletters van MarcoK */

/* PuMa parameters per crate: */

#define NOFCLUSTERS          4 /* Number of clusters in this crate */
#define MAXDISKS             4 /* Number of disks serviced in this crate */
#define MAXPARTITIONS        5 /* 4 data partitions and 1 log & scratch partition */
#define MAXDSPSPERBOARD      6 /* Number of DSPs per board */
#define MAXBOARDSPERCLUSTER  4 /* Number of DSP boards per cluster */
#define MAXDSPBOARDS        16 /* Number of DSP boards per crate */
			       /* Number of DPSs in this crate */
#define MAXDSPS             (NOFCLUSTERS * MAXBOARDSPERCLUSTER * MAXDSPSPERBOARD)

/* WSRT parameters: */

#define MAXFREQBANDS   8 /* Number of WSRT frequency bands */
#define MAXTELESCOPES 14 /* Number of radio telescopes in the array */



/* R E Q U E S T S */

/* Requests from TMS to PuMa: */

#define NoRequest         0 /* PuMa supposed to be idle */
#define GoOffLine        -1 /* A message to PuMa daemon or operator to allow maintenance on PuMa */
#define Initialise        1 /* Prepare for an observation */
#define StartObservation  2 /* On next 10" tick start the observation */
#define AbortObservation  3 /* Unconditional orderly abort of observation */
#define GetNofFrames      4 /* Report number of frames taken so far */
#define GetDiskInfo       5 /* Give capacity and occupancy of all disks */
#define GetCrateNumber    6 /* Give the hardware number of this crate */
#define CopyToArchive     7 /* The filenames are in a specified TMS file */
#define AbortCopy         8 /* Unconditional orderly abort of CopyToArchive */
#define DeleteFiles       9 /* The filenames are in a specified TMS file */

/* Requests from PuMa to TMS: */

/*      NoRequest        0    Already defined */
#define GetGeneralInfo   1 /* Obtain the data for _all_ the fields in the TMS area */
#define GetTime          2 /* Obtain the current date & time (in steps of 1") */
#define WantOffLine      3 /* Operator requests PuMa off line for maintenance, confirmed by GoOffLine */
#define NewTape          4 /* PuMa requires a new DAT tape */


/* P u M a   E R R O R   M E S S A G E S */

/* These error messages are single bits and can thus be combined.
   Some of these error messages just indicate that something is wrong,
   the "(StatusType) itsPuMaStatus" will indicate what pieces of hardware
   are affected.                                                          */

#define NoError            0x00000000
#define AnyError           0x00000001 /* So the exit code of 'system' will show that there was any error */
#define DSPDown            0x00000002 /* At least one DSP is down */
#define DSPBoardDown       0x00000004 /* The master DSP is not responding */
#define ClusterDown        0x00000008 /* A whole, but required, cluster is down */
#define AnalogChannelDown  0x00000010 /* An A/D & Clock section is down of an XYchannel */
#define TemperatureWarning 0x00000020 /* Could be in the analog or the digital section */
#define PowerSupplyWarning 0x00000040 /* Could be in the analog or the digital section */
#define LastRequestFailed  0x00000080 /* For some reason PuMa couldn't carry out the last TMS request */
#define DiskFull           0x00000100 /* PuMa ran out of disk space during an observation */
#define OutOfSync          0x00000200 /* PuMa's own sense of time is off */
#define ClashesHappened    0x00000400 /* The data input PAL on a DSP board couldn't get the local bus */
#define HPMissedFrame      0x00000800 /* HP too slow or housekeeping error */
#define WrongScanObsNumber 0x00001000 /* Files of this scan already exist */

/* redwards commented out, breaks swin code */
/* XXXX #define FileNotFound       0x00002000*/ /* Like the PuMa obs data files not present */
#define DeleteFailed       0x00004000 /* Couldn't delete a file */
#define DataCorrupt        0x00008000 /* Wrong SHARC ids, etc. */
#define SHARC_Message      0x00010000 /* Any warning or error from any of the SHARCs */

#define KilledBySignal     0x00400000 /* Like Ctrl+C or run time error */
#define TapeError          0x00800000 /* Taping to the DAT drive failed */
#define LongTapeDelay      0x01000000 /* Expected tape time exceeded */
#define CopyAborted        0x02000000 /* for whatever reason */
#define ObsAborted         0x04000000 /* for whatever reason */
#define ForcedToAbort      0x08000000 /* TMS or 'operator' requested an abort */
#define RevisionError      0x10000000 /* Mismatch between TMS and actual RCS identifier */
#define SoftwareError      0x20000000 /* For debugging purposes only!... */
#define ExceptionReport    0x40000000 /* PuMa couldn't carry out the requested observation */
#define UnknownError       0x80000000 /* Something that shouldn't occur... */


/* P u M a   S T A T E S */

/* Most of these states are single bits and can thus be combined */

#define PuMaBootAndTest   -1 /* After power-up, testing all hardware */
#define PuMaIdle           0 /* What a boring life for a PuMa . . . */
#define PuMaInitialising   1 /* Preparing for an observation */
#define PuMaAwaitStart     2 /* Ready for the orchestrated "count down" from TMS */
#define PuMaObserving      4 /* Observation in progress */
#define PuMaAborting       8 /* Orderly aborting upon a TMS abort request */
#define PuMaArchiving     16 /* Copying certain disk files to archive */
#define PuMaDeleting      32 /* Deleting certain disk files */


/* ----------------------------- T Y P E S -------------------------------- */


/* G E N E R A L   T Y P E S */

typedef int      RequestType;           /* See list of requests above */
typedef int      RequestStatus;         /* To store NOREQUESTRECEIVED or WORKINGONREQUEST */
typedef long int MJD_t;                 /* Modified Julian Date, 5 digit number */
typedef char     ProgNameType[PATHLEN]; /* To store the prog name of the shared mem */

typedef enum {Observation, Test, Calibration} ObsReason;

typedef struct { int i0, i1;} ExtInt;


/* T M S   T Y P E S */

typedef struct { double itsDeclination;    /* [Radians] */
		 double itsRightAscension; /* [Radians] */
               } SkyPosType;

typedef struct { char       itsName[NAMELEN];  /* May be "Artificial Pulsar" */
                 ObsReason  itsReason;
                 SkyPosType itsSkyPosition;
                 char       itsEpoch[NAMELEN]; /* e.g. J2000, B1950, Apparent */
               } SourceType;

typedef struct { char       itsObsID   [NAMELEN]; /* A unique WSRT observation ID */
		 char       itsCalObsID[NAMELEN]; /* Its calibration counterpart (if any) */
                 int        itsModeNo;
		 MJD_t      itsStartDate;      /* Steps up when UTC --> 0 */
		 PUMA_STD time_t     itsStartTime;      /* In UTC */
		 long int   itsDuration;       /* [s] */
		 SourceType itsSource;
		 char       itsPI[TXTLEN];     /* Prime Investigator */
		 char       itsProposalNumber[TXTLEN];
	       } ObservationType;

typedef struct { Boolean isGeografical; /* FALSE means geocentral */
		 double  itsLongitude;  /* [Radians] */
		 double  itsLatitude;   /* [Radians] */
	       } CoordType;

typedef struct { char      itsName[TXTLEN];    /* Someting like "WSRT" */
		 CoordType itsCoord;
		 double    itsHeight;          /* [m] above mean sea level */
		 MJD_t     itsCurrentDate;     /* Steps up when UTC --> 0 */
		 PUMA_STD time_t    itsCurrentTime;     /* [s] UTC */
		 int       itsGPS_MASERoffset; /* [ns]  Is the +/- sign well defined?
							============================= */
		 Boolean   isUniverseClosed;   /* Beyond the scope of PuMa */
               } ObservatoryType;

typedef struct { Boolean isAvailable;               /* The RT is on-line */
                 char    itsFrontendID[TXTLEN];     /* Unique ID */
		 char    itsFrontendStatus[TXTLEN]; /* What sky freq band */
		 Boolean isNoiseSrcOff;             /* The noise source in teh front ends */
		 double  itsSysTemp;                /* Estimated system temperature in K */
               } TelescopeType;

typedef struct { double  itsMidFreq;    /* Something like 6.25 or 5 MHz */
                 double  itsWidth;      /* Something like 10, 5 or 2.5 MHz */
                 double  itsMidSkyFreq; /* Sky freq. corr. with itsMidFreq */
                 Boolean isNeeded;      /* Should this band be processed */
                 Boolean isFlipped;     /* Increasing or decreasing sky frequency */
	       } BandType;

typedef struct { char itsDir[DIRLEN];      /* Full (sub-)dir name where PuMa software is stored */
		 char itsFile[FILELEN];    /* Name of a file (data or executable) */
		 char itsRevision[TXTLEN]; /* The unique release number of the software in use */
	       } SWType;

typedef struct { Boolean         isOffsetDynamic;               /* Or fixed value */
		 Boolean         isScaleDynamic;                /* Or is it a fixed value */
		 long int        itsAdjustInterval;             /* [s] period */
	       } AdjustType;

typedef struct { char            itsComment[COMMENTLEN];        /* Any text */
		 char            itsCrateNumber;                /* CRATE0, CRATE1, BOTHCRATES */
		 ObservatoryType itsObservatory;
		 ObservationType itsReqObs;
		 TelescopeType   itsTelescopes[MAXTELESCOPES];
		 BandType        itsFreqBands    [MAXFREQBANDS];
		 int             itsBandToClusMap[MAXFREQBANDS];/* How the cables are connected */
		 AdjustType      itsAdjust;                     /* Dynamic updating of Offset and ScaleFactor */
		 SWType          itsDSPTestProg;                /* Initial test prog for SHARCs */
		 char            itsWarnDigiTemp;               /* Warning T for SHARC boards ºCelsius */
		 int             itsDATTapeSize;                /* in MB for archiving of astro files */
               } InitType;

typedef struct { Boolean isXUsed;
                 Boolean isYUsed;
               } PolUsedType;

typedef struct { SWType itsSW;
		 int    itsDecimation;
	       } FilterType;

typedef struct { InitType    itsInit;                              /* All modes */
                 SWType      itsMasterSoftware;                    /* All modes */
		 SWType      itsDSPSoftware [NOFCLUSTERS][MAXDSPSPERBOARD];
		 float       itsScaleFactor [NOFCLUSTERS][MAXPACKS];/* Mode 0..4 */
		 float       itsDCvalue     [NOFCLUSTERS][MAXPACKS];/* Mode 0..4 */
		 float       itsXPolScaleFac[NOFCLUSTERS];          /* Multiplication factor for X pol relative to Y */
		 double      itsNofSigmas;                         /* # of sigmas mapped onto outsample range */
                 PolUsedType itsBandPol;                           /* Mode 0 */
                 int         itsNofBitsPerSample;                  /* All modes */
                 FilterType  itsFilter;                            /* Mode 0, 4 */
                 Boolean     isStokesParamNeeded[MAXSTOKESPARAMS]; /* Mode 1..3 */
                 int         itsNofFFTFreqChannels;                /* Mode 1..3 */
                 int         itsNofSampsToAdd;                     /* Mode 1..3 */
                 int         itsNofSHARCsToAdd;                    /* Mode 1..3 */
                 double      itsDM[MAXDMS];                        /* Mode 2,
                                                              mode 3 and 4 use "itsDM[0]" only */
		 double      itsRM;                                /* Mode 3,4 */
               } ModeType;

typedef struct { char itsPath[PATHLEN]; /* Full path of file with file names for copy or del */
               } TaskType;


/* P u M a   T Y P E S */

typedef struct { long int itsMaxSpace [MAXPARTITIONS]; /* [MB] Total harddisk capacity */
                 long int itsFreeSpace[MAXPARTITIONS]; /* [MB] Total unused harddisk capacity */
               } DiskType;

typedef struct { Boolean isPos5VoltOK;
		 Boolean isNeg5VoltOK;
	       } AnaPwrType;

typedef struct { Boolean isPos5VoltOK;
	       } DigPwrType;

typedef struct { Boolean    isClusterOK;       /* The whole cluster is OK */
                 Boolean    isAnalogOK;        /* The analogue XYchannel is OK */
              /* Boolean    isFrequencyBandOK is not needed, because
                            isFrequencyBandOK = (isAnalogOK && isClusterOK) */
                 Boolean    isAnalogTempOK;    /* Temp < TempAlert of analog XYchannel */
                 AnaPwrType itsAnaPwrSupply;   /* Power supplies of analog XYchannel */
               } ClusterWarnType;

typedef struct { long int   itsNofClashes;    /* Data input local-bus clashes, should be 0 */
                 Boolean    isDSPBoardOK;     /* This board is OK */
                 Boolean    isDSPBoardTempOK; /* Temperature < TempAlert of DSP board */
                 DigPwrType itsDigPwrSupply;  /* Power supplies of the DSP board */
               } DSPBoardWarnType;

typedef struct { long int itsNofObsFiles; /* Number of files of last observation of this disk */
		 char     itsDir[DIRLEN]; /* The full dir on this disk where the data are stored,
					     if the disk is used (#ofFiles>0) then also the log
					     file and the exception file are stored there      */
               } StoreLocType;

typedef struct { char              itsCrateNumber;             /* Hardware crate number 0, 1 */
		 volatile long int itsState;                   /* See list of states above */
		 unsigned long int itsErrorMessage;            /* See list of error messages */
		 long int          itsNofCounter;              /* See requests from TMS to PuMa */
		 long int          itsTotCounter;              /* Total frames per SHARC to process */
		 char              itsUnit[TXTLEN];            /* The unit, like files or frames */
		 ClusterWarnType   itsCluster  [NOFCLUSTERS];
		 char              itsADCNumber[NOFCLUSTERS];  /* Unique ADC hardware number */
		 DSPBoardWarnType  itsDSPBoard   [MAXDSPBOARDS];
		 ExtInt            itsBoardNumber[MAXDSPBOARDS];
		 DiskType          itsDisks       [MAXDISKS];
		 PUMA_STD time_t            itsActualStartTime;         /* In seconds */
		 StoreLocType      itsLastObsFiles;
		 Boolean           isDSPOK[MAXDSPS];           /* The DSP is OK */
		 /* Feedback to TMS; can be used as start values for the next observation (ModeType) */
		 float             itsLastScaleFactor[NOFCLUSTERS][MAXPACKS];
		 float             itsLastDCValue    [NOFCLUSTERS][MAXPACKS];
		 float             itsLastXYRatio    [NOFCLUSTERS];
		 int               itsLastModeNumber;
	       } StatusType;


/* S H A R E D   M E M O R Y   A R E A S */

/* The TMS area: */
typedef struct { volatile RequestType   itsRequestTMStoPuMa;
			  RequestStatus itsRequestState;
			  ModeType      itsMode;
			  TaskType      itsOtherTask;
			  ProgNameType  itsProgName;
                 volatile Boolean       isActive;
               } TMSAreaType;

/* The PuMa area: */

typedef struct { volatile RequestType   itsRequestPuMatoTMS;
                          RequestStatus itsRequestState;
                          StatusType    itsPuMaStatus;
                          ProgNameType  itsProgName;
		 volatile Boolean       isActive;
               } PuMaAreaType;


/* ----------------------- D E C L A R A T I O N S ------------------------ */

#define  TMSSHMNAME "/shmemtms"
#define PuMaSHMNAME "/shmempuma"

//static char *puma_h_RCSid = "$Header: /home/bws/Repository/src/PuMa/libs/puma/puma.h,v 1.6 2004/01/05 13:41:17 redwards Exp $";




#endif  /* part that's excluded when TIME is defined, see above */

#endif  /* this whole header file */

/* ============ End of this PuMa sheared memory header file =============== */



/* $Header: /home/bws/Repository/src/PuMa/libs/puma/pumadata.h,v 1.4 2004/09/06 14:43:16 straten Exp $ */

#ifndef _PUMADATA_H
#define _PUMADATA_H

#define HEADERFORMATID "DPC_1.3" /* Dutch Pulsar Community */

//#include "puma.h"    /* for all the shared info like array sizes and exitcodes */

/* =========================================================== */
/*                                                             */
/* pumadata.h - header file for astronomical PuMa data files   */
/*   ----------- IMPORTANT NOTE ------------------------------ */
/*   ----------- IMPORTANT NOTE ------------------------------ */
/*   ----------- IMPORTANT NOTE ------------------------------ */
/*   ----------- IMPORTANT NOTE ------------------------------ */
/* NOTE!!! Changes here need also be made in the CVS, lynx and */
/*                        HPRT    versions!!!!!!!!!!!!         */
/*                                                             */
/*           ---------- DPC_1.3 ----------                     */
/* Change in the header so that there is a variable which      */
/* Describes the changes in the software versions for the      */
/* reduction software. Called RednSoftwareVer.                 */
/*     BWS 2003-06-04                                          */
/*           ---------- DPC_1.2 ----------                     */
/* Change in the adjustment method used                        */
/*            ---------- DPC_1.1 ----------                    */
/*                                                             */
/* Version: 24 Aug 99 - Marco, Lodie, Ben and Ramach.          */
/*                                                             */
/* BEWARE : At some places this file is 160 characters wide    */
/*                                                             */
/* =========================================================== */


/* =========================================================== */
/*                                                             */
/* General layout of a PuMa data file:                         */
/*                                                             */
/*     -------------------                                     */
/*    |                   |  Information about the observation */
/*    |      HEADER       |  fixed length of 4200 bytes        */
/*    |                   |  = sizeof(Header_type).            */
/*     -------------------                                     */
/*     -------------------                                     */
/*    |    ADJUSTMENT     |  Scaling and DC offset values, its */
/*    |    PARAMETER      |  length depends on the observation */
/*    |    BLOCK          |  duration and adjustment interval. */
/*     -------------------   Can also be polyco block          */
/*     -------------------                                     */
/*    |                   |                                    */
/*    |       DATA        |  Astronomical data, variable in    */
/*    |       BLOCK       |  length.                           */
/*    |                   |                                    */
/*    |                   |                                    */
/*    |                   |                                    */
/*    |                   |                                    */
/*    |                   |                                    */
/*    |                   |                                    */
/*    |                   |                                    */
/*    |                   |                                    */
/*    |                   |                                    */
/*    |                   |                                    */
/*    |                   |                                    */
/*     -------------------                                     */
/*                                                             */
/*                                                             */
/* The data will be divided in +/- 200 MB blocks (or smaller). */
/* Each block will be stored in a separate file with a         */
/* separate header.                                            */
/*                                                             */
/* Default file length will be +/-200 MB maximum, since if data*/
/* for each cluster is in a separate file, 8x200 MB = 1.6 GB   */
/* can be handled on a 2GB disk while still leaving some spare */
/* capacity.                                                   */
/*                                                             */
/* If the data of one observation is more than the maximum     */
/* capacity of one tape, it will be divided and spread over    */
/* more tapes, each data file will start with a header.        */
/*                                                             */
/* All data is taken in one mode. If clusters or crates are    */
/* in different modes, different output files will be          */
/* produced.                                                   */
/*                                                             */
/* Datafiles will be put on tape as TAR files.                 */
/*                                                             */
/* =========================================================== */


/* =========================================================== */
/*                                                             */
/* Log file contains the binary and ASCII version of the       */
/* header but contains also information about                  */
/* - what went wrong if ExitCodes are not 0                    */
/* - RCS numbers of sub-software                               */
/* - ID numbers of boards and ADC units                        */
/* - More details about crate and ADC unit status              */
/* - Planned starttime and reason for delay                    */
/*                                                             */
/* =========================================================== */


/* =========================================================== */
/*                                                             */
/* In this header the following sizes of variable types        */
/* are assumed:                                                */
/*                                                             */
/*        char           = 1 byte                              */
/*        Dummy          = 1 byte  (own definition)            */
/*        Boolean        = 4 bytes (defined as integer)        */
/*        int            = 4 bytes                             */
/*        float          = 4 bytes                             */
/*        double         = 8 bytes                             */
/*                                                             */
/*        N.B. the long int is no longer used                  */
/*                                                             */
/* Running a 'sizeof'-program gives :                          */
/*                                                             */
/*   Telescope_type   :  144  bytes                            */
/*   Band_type        :   32  bytes                            */
/*   General_type     :  488  bytes                            */
/*   Observatory_type :   40  bytes                            */
/*   Observation_type :  248  bytes                            */
/*   Target_type      :   48  bytes                            */
/*   Signalpath_type  : 2344  bytes                            */
/*   Mode_type        :  264  bytes                            */
/*   Software_type    :  512  bytes                            */
/*   Check_type       :   64  bytes                            */
/*   Reduction_type   :  496  bytes                            */
/*                      ----------- +                          */
/*   Header_type      : 4504  bytes    correct ...             */
/*                                                             */
/* =========================================================== */

/* --------  Software version control Integer     ------------ */
/*   Added by BWS in 4 June 2003 to take into account the possibility of major revisions
     in data reduction software that affect the output - reduced data in a major way.
     Every time that such a major change is introduced this integer should be increased
     and the reason for the change should be clearly mentioned and documented as well as
     the name of the software that was changed and also the date and author of the change.
     Software down the line from the software in which the change has been implemented will
     be able to then use this integer as a flag to see what if any modifications are necessary. */


# define defRednSoftwareVer 1

/* Explanations for  RednSoftwareVer setting */
/*    1         BWS        4 June 2003
      --- The impletation and use of this variable was started because of the incorrect frequency
      labelling of frequency channels which occured in both fpuma and fpumanr before this date.
      This results in data which is reduced before 4 June 2003 (and thus which have no variable
      called RednSoftwareVer in the header) having all frequency channels in observations which
      record Flipped bands being labelled with frequencies which are BW/(2*Nchan) too low. While
      those observations of Non-Flipped bands are labeeled with frequencies which are BW/(2*Nchan)
      too high. For example for a 10 MHz band with 64 frequency channels and a flipped band channel 0
      was being labelled with frequency 9.921875 MHz when in fact it should have been labelled 10.0 MHz.
      This results in the data being incorrectly dedispersed and in the case of Flipped data pulse(s) will appear
      to arrive earlier than it(they) should while in Non-Flipped data pulse(s) will be appear to arrive
      later than it(they) should. It will also result in some smearing of the pulse but this is (I believe)
      less than the sampling interval all cases I know of.
      This of course has important implications for both Timing and Alignment.
      Note that this only affects Mode 1 data and from now on is fixed in fpuma and fpumanr. */


/* --------- G E N E R A L   P A R A M E T E R S  ------------ */

#define ASTROFILESIZE (200 * 1024 * 1024) /* 200 MB - approximate size of tape files */

#define LONGNAMELEN 24


/* ----- G E N E R A L   V A R I A B L E   T Y P E S   ------- */


typedef char     Dummy;                 /* Dummy variable  */


/* -----------------   S U B - T Y P E S   ------------------- */


/*                 TYPE          NAME                              EXAMPLE            DET.BY  DESCRIPTION                                                 */


typedef struct {   Boolean       Active;                        /* TRUE               TMS     Is telescope operational and in adding box?                 */
		   Boolean       NoiseSrcOff;                   /* TRUE               TMS     Is internal 10s noise source off? TRUE is needed            */
		   char          FrontendID[TXTLEN];            /* "MFFE-dd"          TMS     dd is an index                                              */
		   char          FrontendStatus[TXTLEN];        /* "operational"      TMS     Operational/in repair/..                                    */
		   double        Tsys;                          /*                    TMS     Estimate/actual averaged system temp. (in K)                */
	       } Telescope_type;


typedef struct {   Boolean       Operational;                   /* TRUE               TMS     Is band o.k?                                                */
                   Boolean       NonFlip;                       /* TRUE               TMS     Is band flipped (reverse freq order within band)?           */
                   double        MidFreq;                       /* 6.25               TMS     Mid frequency of band (in MHz): 5.0 or 6.25                 */
                   double        Width;                         /* 10.0               TMS     Width of the band (in MHz): 2.5, 5.0 or 10.0                */
                   double        SkyMidFreq;                    /* 1405.0005          TMS     Mid sky frequency of band (in MHz)                          */
               } Band_type;




/* ---------------------   T Y P E S   ----------------------- */


/*                 TYPE          NAME                              EXAMPLE            DET.BY  DESCRIPTION                                    */


typedef struct {   char          HdrVer[NAMELEN];               /* "PuMa1.0"          HP      To distinguish future formats                  */
		   char          Platform[NAMELEN];             /* "PuMa"             HP      Name of machine: "PuMa" or "FFB"               */
		   char          ThisFileName[LONGNAMELEN];     /* ScanNum.<filenumber>.<clusternumber{0..MAXFREQBANDS-1}>                   */
		   char          ScanNum[NAMELEN];              /* "200100002                                                                */
								/*  YYYYWWWWW-FFC"               - FileSeriesNumber(FF)+Cluster(C).          */
								/*                               Change TMS year '98' into '1998'            */
		   char          Comment[COMMENTLEN];           /*                    TMS     Comment entered by operator                    */
		   int           NFiles;                        /*                    HP      Number of files over which this obs. is spread */
								/*                               starting with 0.                            */
								/*                               Beware: each cluster in separate file       */
		   int           FileNum;                       /*                    HP      Which file out of NFiles is this               */
		   char          TapeID[LONGNAMELEN];           /* "200100002A"       HP      Unique TapeID on which this data is saved      */
								/*                               Format: Unique WSRT obs id + letter(A,B,..) */
		   int           NTapes;                        /* 3                  HP      Number of tapes over which this obs. is spread */
		   int           TapeNum;                       /* 2                  HP      Which tape out of NumberOfTapes is this        */
								/*                               starting with 0                             */
		   int           ParBlkSize;                    /*                    HP      Number of bytes in second (parameter) block    */
		   int           DataBlkSize;                   /*                    HP      Number of bytes in data block                  */
		   Boolean       Cluster[MAXFREQBANDS];         /*                    HP      Is Data from this cluster in this file?        */
		   int           DataMJD;                       /*                    HP      MJD of first sample of this data block         */

                   Dummy AlignDouble[4];

		   double        DataTime;                      /*                    HP      Time of first sample (fraction of day)         */
		   Dummy         GeneralDummy[64];
	       } General_type;


typedef struct {   char          Name[NAMELEN];                 /* "WESTERBORK"       TMS     Name of the observatory                        */
                   double        Long;                          /* -0.11526451        TMS     Geocentral west longitude of observatory (rad) */
                   double        Lat;                           /* 0.92357448         TMS     Geocentral latitude of observatory (rad)       */
                   double        Height;                        /* 5                  TMS     Height of observatory above sea level (metres) */
               } Observatory_type;


typedef struct {   char          ObsName[TXTLEN];               /* "R/98A/01"         TMS+HP  Proposal name (source name from TMS extracted  */
                   char          ObsType[TXTLEN];               /*                    TMS     Research/test/calibration                      */
                   char          LastCalID[NAMELEN];            /* "200100001"        TMS     ScanNumber of last calibration observation     */
                   int           StMJD;                         /* 50432              TMS+HP  MJD at start of observation                    */
                   int           StTime;                        /* 43210              TMS+HP  Starttime (s after midnight, multiple of 10 s) */
                   double        StLST;                         /*                    TMS     Start Siderial time (fraction of 24 h)         */
                   double        MaserOffset;                   /* 0.0000432          TMS     (Averaged over obs.) offset of Maser (in s)    */
		   double        Dur;                           /* 3590.001           HP      Actual duration of observation (in s)          */
                                                                /*                            = multiple of sampling time                    */
                   double        IonRM;                         /*                    TMS     Future: ionospheric RM towards target          */
                   Dummy         ObservationDummy[64];
               } Observation_type;


typedef struct {   char          Pulsar[NAMELEN];               /* "PSR J0218+4232"   TMS     Using the Taylor (1993) convention             */
                   double        RA;                            /*                    TMS     RA of the target (in radians)                  */
                   double        Dec;                           /*                    TMS     Dec of the target (in radians)                 */
                   char          Epoch[NAMELEN];                /* "J2000.0"          TMS     Epoch of target coordinates                    */
               } Target_type;


typedef struct {   Telescope_type Tel[MAXTELESCOPES];           /*                    TMS     Active/FrontendID+Status/Syst.Temp             */
		   Band_type     Band[MAXFREQBANDS];            /*                    TMS     MidFreq/Width/MidSkyFreq/Operational/Flipped   */
                   char          Backend[NAMELEN];              /*                    TMS     DCB/DZB/DLB/DXB                                */
                   Boolean       ResFringeOff;                  /* TRUE               TMS     Is residual fringe stopping off? TRUE is needed*/
                   Boolean       AddBoxOn;                      /* TRUE               TMS     Is the adding box working? TRUE is needed      */
		   int           BandsToClusMap[MAXFREQBANDS];  /* 4,3,2,1,0,7,6,5    TMS     Band 0 is connected to cluster 4               */
                   Dummy         WSRTDummy[16];
	       } Signalpath_type;


typedef struct {   int           Nr;                            /*                    TMS     Mode of PuMa observation -1,0,1,2,4            */
                   Dummy         SpareDummy1[4];                /*                            To avoid problems with 64 bits machines        */
		   float         XPolScaleFac[MAXFREQBANDS];    /*                            Multiplication factor of X pol realtive to Y   */
                   Boolean       ActCluster[MAXFREQBANDS];      /*                    TMS+HP  BandOperational.and.AdcOk.and.ClusterOk        */
		   int           FIRFactor;                     /*                    TMS     Mode 0,4 ; 1,2,4,8,...     Mode -1,1,2 : 1     */
		   int           NSHARCsAdded;                  /*                    TMS     Mode 1,2 ; 1,2,3,6         Mode -1,0,4 : 1     */
		   int           NSampsAdded;                   /*                    TMS     Mode 1,2 ; 1,2,4,8,16,.. depends on FreqChans  */
		   int           NFreqInFile;                   /*                    TMS     Mode 0-4 ; 1,2,4,8,16,..FreqChans in this file */
								/*                               Normally one cluster per file is written    */
                   int           Tsamp;                         /*                    TMS     Mode 0-4 ; output sample interval in nano sec  */
                   Boolean       Iout;                          /*                    TMS     Mode 1,2 ; I pol output    Mode -1,0,4 : FALSE */
                   Boolean       Qout;                          /*                    TMS     Mode 1,2 ; Q pol output    Mode -1,0,4 : FALSE */
                   Boolean       Uout;                          /*                    TMS     Mode 1,2 ; U pol output    Mode -1,0,4 : FALSE */
                   Boolean       Vout;                          /*                    TMS     Mode 1,2 ; V pol output    Mode -1,0,4 : FALSE */
                   Boolean       Xout;                          /*                    TMS     Mode -1,0; X pol output    Mode 1,2,4  : FALSE */
                   Boolean       Yout;                          /*                    TMS     Mode -1,0; Y pol output    Mode 1,2,4  : FALSE */
                   Dummy         SpareDummy2[8];                /*                    HP      Future: Stokes P,Theta??                       */
                   int           NDMs;                          /*                    TMS     Mode 2 ; Num calc.dedisp.series  Mode: -1,4: 1 */
                   Dummy         SpareDummy3[4];                /*                                                             Mode: 0,1 : 0 */

                   Dummy AlignDouble[4];

		   double        DM[MAXDMS];                    /*                    TMS     Mode 2,4 ; Mode -1: DM[0] = 0                  */
                   double        RM;                            /*                    TMS     Mode 2                                         */
                   Boolean       DC_Dynamic;                    /*                    TMS     Is DC dynamic removed?                         */
                   Boolean       ScaleDynamic;                  /*                    TMS     Is scaling dynamic adjusted                    */
                   double        AdjustInterval;                /*                    TMS     Adjustment interval                            */
                   Boolean       FloatsOut;                     /*                    TMS     Is the data written in floats ?                */
                   int           BitsPerSamp;                   /*                    TMS     Number of output bits ; 1,2,4,8                */
		   double        SigmaRange;			/*                    TMS     Number of Sigmas mapped onto output range      */
                   Dummy         ModeDummy[64 - 8];             /* removed 8 bytes for SigmaRange */
	       } Mode_type;


typedef struct {   char          Master[TXTLEN];                /*                    TMS     Unique RCS number of Mastersoftware on HP      */
                   char          DSP[MAXDSPSPERBOARD][TXTLEN];  /*                    TMS     RCS number of DSP software                     */
		   char          Filter[TXTLEN];                /*                    TMS     RCS number decimation filter softw., Mode 0,4  */
	       } Software_type;


typedef struct {   int           ExitCodeObsy;                  /*                    TMS     Exit Code of observatory                       */
		   int           ExitCodePuMa;                  /*                    HP      Exit Code of PuMa                              */
		   int           ExitCodeClust[MAXFREQBANDS];   /*                    HP      Exit Code of clusters                          */
		   int           ExitCodeDataConsistency;       /*                    HP      Exit Code after checking the timeseries        */
		   Dummy         CheckDummy[20];
	       } Check_type;


typedef struct {   Boolean       Raw;                           /*                            In case of raw data                            */
		   int           MJDint;                        /*                            Integer part of MJD for reduced data           */
		   double        MJDfrac;                       /*                            Fractional part of MJD for reduced data        */
                                                                /*                              Midpoint of data in case of folded data      */ 
                                                                /*                              Starttime in all other cases                 */ 
		   double        DM;                            /*                            Used DM during reduction                       */
		   int           NTimeInts;                     /*                            Number of time subintegrations (samples for    */
								/*                              normal timeseries) in this file              */
		   int           NFreqs;                        /*                            Number of frequency bands in this file         */
		   double        DeltaTime;                     /*                            Delta time between subintegrations/samples     */
		   double        FreqCent;                      /*                            Central frequency of all these intervals       */
		   double        DeltaFreq;                     /*                            Delta freq. between bands (assume equidistant) */
		   Boolean       IsDedisp;                      /*                            In case of dedispersed data                    */
                   Boolean       IsCohDedisp;                   /*                            In case of coherently dedispersed data         */
                   int           CohFFTSize;                    /*                            Size of FFT during coherent dedispersion       */ 
		   Boolean       IsPwr;                         /*                            In case of powerspectrum                       */
		   char          Zapfile[NAMELEN];              /*                            Name of zap file or 'default'                  */
		   int           NBins;                         /*                            Number of bins in profile (or in spectrum)     */
		   Boolean       Folded;                        /*                            In case of folded data                         */
		   double        FoldPeriod;                    /*                            (Average ?) period at which is folded (s)      */
		   Boolean       Polyco;                        /*                            In case of polyco driven folding               */
		   int           NCoef;                         /*                            Number of polyco coefficients                  */
		   int           PolycoSpan;                    /*                            Span of polyco                                 */
                   Boolean       IsAdjusted;                    /*                            In case adjustments have been applied          */
                   Boolean       PolycoStored;                  /*                            Polyco stored in Parameter Block               */
                                                                /*                            ParBlkSize gives size of Polyco (in bytes)     */
		   Boolean       Bary;                          /*                            In case of barycentring                        */
		   Boolean       OI;                            /*                            Stokes I in output                             */
		   Boolean       OQ;                            /*                            Stokes Q in output                             */
		   Boolean       OU;                            /*                            Stokes U in output                             */
		   Boolean       OV;                            /*                            Stokes V in output                             */
		   Boolean       OP;                            /*                            Stokes P = Q^2 + U^2 in output                 */
		   Boolean       OTheta;                        /*                            Stokes Theta = arctan(U/Q) in output           */
		   Boolean       Op;                            /*                            Stokes p = P/I in output                       */
		   Boolean       Ov;                            /*                            Stokes v = V/I in output                       */
		   Boolean       Opoldeg;                       /*                            poldeg = (Q^2 + U^2 + V^2) / I in output       */
		   Boolean       OX;                            /*                            X in output                                    */
		   Boolean       OY;                            /*                            Y in output                                    */
		   int           TRedn;                         /*                            Time at moment of reduction (sec after 1970)   */
                   char          Command[CMDLEN];               /*                            Commandline of reduction                       */
		   int           RednSoftwareVer;           
                   int           GenType;                       /*                            Type of reduced data (other than:              */
                                                                /*                                 Raw, IsDedisp, IsCohDedisp, IsPwr         */
                                                                /*                                 Folded, Polyco, IsAdjusted, Bary)         */
                                                                /*                             1 = Standard Profile                          */
                                                                /*                             2 = Autocorrelation spectrum                  */
                                                                /*                             4 = Gated data                                */
                                                                /*                             8 = Dynamic spectrum                          */     
                  Dummy         ReductionDummy[56];
	       } Reduction_type;


/* GenType bit fields: */

#define PuMa_Standard_Profile          0x0001
#define PuMa_Autocorrelation_Spectrum  0x0002
#define PuMa_Gated_Data                0x0004
#define PuMa_Dynamic_Spectrum          0x0008
#define PuMa_Reference_Flux_Density    0x0010
#define PuMa_Absolute_Flux_Density     0x0020
#define PuMa_Feed_Corrected            0x0040
#define PuMa_Platform_Corrected        0x0080
#define PuMa_Polarization_Calibrated   0x0100
#define PuMa_RM_Corrected              0x0200


typedef struct {   General_type       gen;
		   Observatory_type   obsy;
		   Observation_type   obs;
		   Target_type        src;
		   Signalpath_type    WSRT;
		   Mode_type          mode;
		   Software_type      software;
		   Check_type         check;
		   Reduction_type     redn;
	       } Header_type;


/* Adjustment Parameter block                                                                     */
/*                                                                                                */
/* The DC offset and the Scale factor are adjusted for each cluster and for each Stokes parameter */
/* In mode 2 (in case of only Stokes I) it will be done for each DM as well                       */
/* It is not done for each frequency channel                                                      */

typedef struct { int   framenumber;
		 int   whatfor;   /* can be PACKX, PACKY, PACKI, PACKQ, PACKU, PACKV  */
		 float scale;     /* for conversion from the clipped int to the 'original' float */
		 float offset;    /* for conversion from the clipped int to the 'original' float */
	       } Adjustments;

/* Then the data follows in a bit stream.  */


/* -------------------------------------------------------------------------------------------------------------------------------------------- */
/*                                                                                                                                              */
/* MODE 0 : The incoming signal is baseband sampled. This mode will be used for off-line coherent dedispersion.                                 */
/* ------   In the reduction process the data will be handled per band, so per cluster. This is also how the data is written.                   */
/*          The data from different crates (set of 4 clusters) will be (almost) always on different tapes.                                      */
/*          In all cases data from different clusters will be in different files.                                                               */
/*                                                                                                                                              */
/*            Do i = 1, number_of_active_clusters                                                                                               */
/*              Do j = 1, number_of_time_samples                                                                                                */
/*                Do k = 1, number_of_active_polarisations                                                                                      */
/*                  Write, outputbits                                                                                                           */
/*                Enddo                                                                                                                         */
/*              Enddo                                                                                                                           */
/*            Enddo                                                                                                                             */
/*                                                                                                                                              */
/*                                                                                                                                              */
/*                             number_of_active_polarisations   =   Mode.XPolarisationUsed + Mode.YPolarisationsUsed                            */
/*                                                                                                                                              */
/*                                     number_of_time_samples   =   Observation.Duration / Mode.OutputSampleInterval                            */
/*                                                                                                                                              */
/*                                                                   7                                                                          */
/*                                  number_of_active_clusters   =   Sum   Mode.ActiveCluster[i]                                                 */
/*                                                                  i=0                                                                         */
/*                                                                                                                                              */
/*          So the data will look like:                                                                                                         */
/*                                                                                                                                              */
/*       ------------------------------------- ... ------------------           ---------------- ... ------------------                         */
/*      | Cluster 0 |         |          |             |  Cluster 0  |         | Cluster 1 |             |  Cluster 1  |                        */
/*      | Sample 0  |         | Sample 1 |             | Last Sample |         | Sample 0  |             | Last Sample |  Etc.                  */
/*      |    X      |    Y    |    X     |             |      Y      |         |    X      |             |      Y      |                        */
/*       ------------------------------------- ... ------------------           ---------------- ... ------------------                         */
/*                                                                                                                                              */
/*                                                                                                                                              */
/*                                                                                                                                              */
/*                                                                                                                                              */
/* MODE 1 : PuMa acts like a filterbank.                                                                                                        */
/* ------   Data from different clusters will be written to different files.                                                                    */
/*          The data from different crates (set of 4 clusters) will be (almost) always on different tapes.                                      */
/*          In all cases data from different clusters will be in different files.                                                               */
/*                                                                                                                                              */
/*            Do i = 1, number_of_active_clusters                                                                                               */
/*              Do j = 1, number_of_time_samples                                                                                                */
/*                Do k = 1, number_of_Stokes_parameters                                                                                         */
/*                  Do l = 1, number_of_frequency_channels                                                                                      */
/*                    Write, outputbits                                                                                                         */
/*                  Enddo                                                                                                                       */
/*                Enddo                                                                                                                         */
/*              Enddo                                                                                                                           */
/*            Enddo                                                                                                                             */
/*                                                                                                                                              */
/*                                                                   3                                                                          */
/*                                number_of_Stokes_parameters   =   Sum   Mode.StokesOutput[i]                                                  */
/*                                                                  i=0                                                                         */
/*                                                                                                                                              */
/*                               number_of_frequency_channels   =   Mode.NumFreqChannelsInFile                                                  */
/*                                                                                                                                              */
/*          So the data will look like (in case of 4 clusters in crate 0) :                                                                     */
/*                                                                                                                                              */
/*   ------------------ ... --------------------- ... ----------------------- ... -------------       ---------- ... ------------               */
/*  | Clus 0 |      |         |        |      |         |        |        |         | Clus 0   |     | Clus 1 |       |          |              */
/*  | Samp 0 |      |         |        |      |         |        | Samp 1 |         | LastSamp |     | Samp 0 |       | LastSamp |              */
/*  | Ch 0   | Ch 1 |         | LastCh | Ch 0 |         | LastCh | Ch 0   |         | LastChan |     | Ch 0   |       | LastCh   |              */
/*  |    I   |      |         |    I   |   Q  |         |    V   |    I   |         |    V     |     |    I   |       |   V      |              */
/*   ------------------ ... --------------------- ... ----------------------- ... -------------       ---------- ... ------------               */
/*                                                                                                                                              */
/*                                                                                                                                              */
/*                                                                                                                                              */
/*                                                                                                                                              */
/* MODE -1: Incoherent dedispersion with DM = 0. No FFT is done, incoming X and Y samples are squared.                                          */
/* -------  X, Y, X^2, Y^2 are saved.                                                                                                           */
/*          Fixed number of samples are added together to bring down the timeresolution. Final sampling time = 81.92 us                         */
/*          Per cluster only one frequency channel comes out.                                                                                   */
/*                                                                                                                                              */
/*                                                                                                                                              */
/*            Do i = 1, number_of_time_samples                                                                                                  */
/*                Write, outputbits_x                                                                                                           */
/*                Write, outputbits_y                                                                                                           */
/*                Write, outputbits_x2                                                                                                          */
/*                Write, outputbits_y2                                                                                                          */
/*            Enddo                                                                                                                             */
/*                                                                                                                                              */
/*          ---------------------------------------- ... --------------                     -------------- ...  --------------                  */
/*         | Crate 0   |    |     |     |                  |           |                   | Crate 0   |          |           |                 */
/*         | Cluster 0 |    |     |     |        |         |           |                   | Cluster 0 |          |           |                 */
/*         | Sample 0  |    |     |     | Samp 1 |         | LastSamp  |                   | Sample 0  |          | LastSamp  |                 */
/*         |     X     | Y  | X^2 | Y^2 |    X   |         |    Y^2    |                   |     X     |          |    Y^2    |                 */
/*          ---------------------------------------- ... --------------                     -------------- ...  --------------                  */
/*                                                                                                                                              */
/*                                                                                                                                              */
/*                                                                                                                                              */
/*                                                                                                                                              */
/* MODE 2, (3,) 4 are not described yet                                                                                                         */
/* -------------- It is proposed to make the loop over the DM's in Mode 2 (only Stokes I) the outerloop                                         */
/*                but it is not yet clear if this is possible                                                                                   */
/*                                                                                                                                              */
/*                                                                                                                                              */
/* -------------------------------------------------------------------------------------------------------------------------------------------- */


/* --------------------------------------------------------------------------------- */
/*                                                                                   */
/* How are the samples stored in the bytes?                                          */
/* ========================================                                          */
/*                                                                                   */
/* When data are read from a pumadata-file into RAM from a given address (char *)    */
/* onwards, then the samples are stored starting with the oldest sample on the       */
/* "right hand side", i.e. it starts at the least-significant-bit (LSB) of that      */
/* given address.                                                                    */
/*                                                                                   */
/* Example: 4 bits wide mode 0 data: 'x0' occupies the 4 lowest bits.                */
/*           In bits 4..7 the oldest 'y0' sample is stored,                          */
/*           in bits 8..11 the next 'x1' sample is stored and so on.                 */
/*                                                                                   */
/*                                                         MSB     LSB               */
/*            +----------+-----------+-----------+---//---+-----------+              */
/*     bits   |31 .... 0 | 31 .... 0 | 31 .... 0 |        | 31 .... 0 |              */
/*            +----------+-----------+-----------+---//---+-----------+              */
/* relative      addr N    addr N-1    addr N-2               addr 0                 */
/*                                                                                   */
/*                                       <----------- Time -----------               */
/*                                                                                   */
/*                             . . . . y x y x y x y x y x y x y x y x               */
/*                             . . . . 7 7 6 6 5 5 4 4 3 3 2 2 1 1 0 0               */
/*                                                                                   */
/* As samples (when in integer format) can be at most 8 bits wide, it is not         */
/* important how the four bytes of a (long) integer are mapped. In this context      */
/* it not relevant if we are dealing with big or little endian.                      */
/* Of course, when the on-board data reduction is large, the accuracy of 8 bits      */
/* could well insufficient and to this end PuMa will have the option to write        */
/* its data as 32-bits floats.                                                       */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */


/*************************************************************************************
 *                                                                                   *
 * A word or two on the various Scales and Offsets as used by PuMa:                  *
 * ================================================================                  *
 *                                                                                   *
 * There are three sets 1. for the user - should be 'intuitive' (see below),         *
 *                      2. for PuMa as used by the SHARCs and                        *
 *                      3. for the off-line reduction software.                      *
 *								                     *
 * Naming convention for the three sets are: UserXXX, PuMAXXX and ConvXXX.           *
 *                                                                                   *
 * The second set are used for mapping the 'float's that result from the             *
 * calculations done by the SHARCs to the clipped 'int's that are stored on PuMa's   *
 * harddisks. The clipping is done on n-bit integers (n in {1,2,4,8}) in the range   *
 * 0..maxval, where 0 is included and maxval is excluded (maxval = 2^n).             *
 *                                                                                   *
 * The third set gives the offset and scale for converting the n-bits clipped        *
 * integers back to the 'original' floats.                                           *
 * To facilitate the conversion by the off-line reduction software, the 'Conv'       *
 * values are stored in the adjustments in the "xxx.PuMa" astro files.               *
 * The conversion is done as in:                                                     *
 *                                                                                   *
 *      'original'_float = (clipped_int - ConvOffset) * ConvScale;                   *
 *                                                                                   *
 * For given PuMaScale and PuMaOffset the following expressions give the conversion  *
 * to the 'Conv' ones:                                                               *
 *                                                                                   *
 *      ConvScale(PuMaScale) = 1 / PuMaScale;                                        *
 * and                                                                               *
 *      ConvOffset(PuMaScale, PuMaOffset) = -PuMaScale * PuMaOffset;                 *
 *                                                                                   *
 *                                                                                   *
 * Some relevant quantities are:                                                     *
 *      a) NSIGMA - number of sigmas of the original floats to be mapped and         *
 *      b) MinP and MaxP - range of the relevant observational parameter.            *
 *      c) NBITS - number of bits for each clipped output sample,                    *
 *         {The quantitity maxval is calculated as 2^NBITS}                          *
 *                                                                                   *
 * With dynamic adjustustment of the PuMaScale and PuMaOffset the following          *
 * calculations are done.                                                            *
 * When the SHARCs prepare their output frames, they enter the clip-routine, where,  *
 * prior to the actual clipping for each observational parameter the floats and      *
 * their squares are accumulated and passed on to the HP743. The latter accumulates  *
 * these accumulated values even further, till either a thread has written the data  *
 * to a harddisk or a specified minimum 'adjustment period' is elapsed, whichever is *
 * the longest.                                                                      *
 * At that time from these two values (lin and sqr) the HP computes the averages and *
 * standard deviations of the noisy parameters and calculates new values for the     *
 * PuMaScale and PuMaOffset as used by the SHARCS in:                                *
 *                                                                                   *
 *      int-before-clipping = (float - PuMaOffset) * PuMaScale;                      *
 *                                                                                   *
 * The conversions from these averages and deviations is done as follows:            *
 *                                                                                   *
 *      PuMaScale(dev) = maxval / (NSIGMA * 2 * dev);                                *
 * and                                                                               *
 *      PuMaOffset(ave, dev) = ave - NSIGMA * dev;                                   *
 *                                                                                   *
 * From these two expressions it will be clear that the whole interval from          *
 *                                                                                   *
 *      ave - NSIGMA * dev .. ave + NSIGMA * dev                                     *
 *                                                                                   *
 * is mapped onto 0..maxval.                                                         *
 * It should be noted that this mapping is not ideal when the distribution of the    *
 * values of a given observational parameter is not symetric.                        *
 *                                                                                   *
 *                                                                                   *
 * The folowing is in case the dynamic scaling is not in use.                        *
 *                                                                                   *
 * For each observational parameter (i.e. for mode0 the X and Y polarisation and for *
 * mode 1 the 4 Stokes parameters) an estimate is made of the maximum range these    *
 * parameters can reach, as this describes the worst case mapping onto the clipped   *
 * integer. This range is denoted by MinP..MaxP.                                     *
 *                                                                                   *
 * The UserScale determines what range of MinP..MaxP is mapped onto 0..maxval,       *
 * so with UserScale = 0.5 only half that range is mapped.                           *
 *                                                                                   *
 * The UserOffset determines where that range is situated in MinP .. MaxP, so with   *
 * UserOffset = 0 that range is situated right in the middle of the MinP .. MaxP     *
 * range, its own bottom is found a quarter of the whole range above MinP.           *
 * Hence with UserOffset = -0.25 the mapped range is lowered by a quarter of the     *
 * whole range and MinP will be mapped to 0 and the mid of the range to maxval.      *
 *                                                                                   *
 * The User and PuMa scales and offsets are coverted as follows:                     *
 *                                                                                   *
 *      PuMaScale(UserScale) = maxval / (UserScale * (MaxP - MinP));                 *
 *                                                                                   *
 *      UserScale(PuMaScale) = maxval / (PuMaScale * (MaxP - MinP));                 *
 * and                                                                               *
 *      PuMaOffset(UserScale, UserOffset) =                                          *
 *        (UserOffset - 0.5 * UserScale) * (MaxP - MinP) +                           *
 *                                     0.5 * (MaxP + MinP);                          *
 *                                                                                   *
 *      UserOffset(PuMaScale, PuMaOffset) =                                          *
 *        (maxval / PuMaSCale + 2 * PuMaOffset - (MaxP + MinP)) /                    *
 *                                            (2 * (MaxP - MinP));                   *
 *                                                                                   *
 *************************************************************************************/


#endif

/* $Header: /home/bws/Repository/src/Soft/putils/bigendian.c,v 1.2 2002/02/26 12:45:36 redwards Exp $ */
int GetBEint(int *src)
{
  int dst;
  
  ((char*)&dst)[0] = ((char*)src)[3];
  ((char*)&dst)[1] = ((char*)src)[2];
  ((char*)&dst)[2] = ((char*)src)[1];
  ((char*)&dst)[3] = ((char*)src)[0];

  return dst;
}

short GetBEshort(short *src)
{
  short dst;
  
  ((char*)&dst)[0] = ((char*)src)[1];
  ((char*)&dst)[1] = ((char*)src)[0];

  return dst;
}

float GetBEfloat(float *src)
{
  float dst;
  
  ((char*)&dst)[0] = ((char*)src)[3];
  ((char*)&dst)[1] = ((char*)src)[2];
  ((char*)&dst)[2] = ((char*)src)[1];
  ((char*)&dst)[3] = ((char*)src)[0];

  return dst;
}

double GetBEdouble(double *src)
{
  double dst;
  
  ((char*)&dst)[0] = ((char*)src)[7];
  ((char*)&dst)[1] = ((char*)src)[6];
  ((char*)&dst)[2] = ((char*)src)[5];
  ((char*)&dst)[3] = ((char*)src)[4];
  ((char*)&dst)[4] = ((char*)src)[3];
  ((char*)&dst)[5] = ((char*)src)[2];
  ((char*)&dst)[6] = ((char*)src)[1];
  ((char*)&dst)[7] = ((char*)src)[0];

  return dst;
}

/*
Boolean GetBEbool(Boolean *src)
{
  if (*src != 0)
    return TRUE;
  return FALSE;
} */

void PutBEint(int *dst, int src)
{
  ((char*)dst)[0] = ((char*)&src)[3];
  ((char*)dst)[1] = ((char*)&src)[2];
  ((char*)dst)[2] = ((char*)&src)[1];
  ((char*)dst)[3] = ((char*)&src)[0];
}

void PutBEshort(short *dst, short src)
{
  ((char*)dst)[0] = ((char*)&src)[1];
  ((char*)dst)[1] = ((char*)&src)[0];
}

void PutBEfloat(float *dst, float src)
{
  ((char*)dst)[0] = ((char*)&src)[3];
  ((char*)dst)[1] = ((char*)&src)[2];
  ((char*)dst)[2] = ((char*)&src)[1];
  ((char*)dst)[3] = ((char*)&src)[0];
}

void PutBEdouble(double *dst, double src)
{
  ((char*)dst)[0] = ((char*)&src)[7];
  ((char*)dst)[1] = ((char*)&src)[6];
  ((char*)dst)[2] = ((char*)&src)[5];
  ((char*)dst)[3] = ((char*)&src)[4];
  ((char*)dst)[4] = ((char*)&src)[3];
  ((char*)dst)[5] = ((char*)&src)[2];
  ((char*)dst)[6] = ((char*)&src)[1];
  ((char*)dst)[7] = ((char*)&src)[0];
}

/*
void PutBEbool(Boolean *dst, Boolean src)
{
  if (src != FALSE)
    *dst = 0x01000000;
  else
    *dst = 0;
} */

/* redwards func to swap an arbitrary word */
/* swaps from the out-in */
void swapWord(void *xptr, int size)
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
void swapArray(void *xptr, int nelem, int size)
{
  int i;
  unsigned char *cptr = (unsigned char *)xptr;
  for (i=0; i < nelem; i++)
  {
    swapWord(cptr, size);
    cptr += size;
  }
}


/* Program to convert a header in big endian format to little endian
   format so that it can be read on a PC. It Makes us of the nice
   utilities in bigendian.c/.h written by RvdK. BWS 2001-03-02 */
void beheader_convert(Header_type *src, Header_type *dst)
{

  int i;

  /* Making them equal ensures that we only have to change the floats,
     doubles and booleans */
  *dst = *src;

  /* Let us do the conversion type by type. Note we treat Booleans as ints here. */
  
  /* Signal Path */
    for (i=0; i< MAXTELESCOPES; i++) {
      /* Ints */
      dst->WSRT.Tel[i].Active = GetBEint(&src->WSRT.Tel[i].Active);
      dst->WSRT.Tel[i].NoiseSrcOff = GetBEint(&src->WSRT.Tel[i].NoiseSrcOff);

      /* Doubles */
      dst->WSRT.Tel[i].Tsys = GetBEdouble(&src->WSRT.Tel[i].Tsys);
  }

    for (i=0; i< MAXFREQBANDS; i++) {
      /* Ints */
      dst->WSRT.Band[i].Operational = GetBEint(&src->WSRT.Band[i].Operational);
      dst->WSRT.Band[i].NonFlip = GetBEint(&src->WSRT.Band[i].NonFlip);
      dst->WSRT.BandsToClusMap[i] = GetBEint(&src->WSRT.BandsToClusMap[i]);
      dst->gen.Cluster[i] = GetBEint(&src->gen.Cluster[i]);


      /* Doubles */
      dst->WSRT.Band[i].MidFreq = GetBEdouble(&src->WSRT.Band[i].MidFreq);
      dst->WSRT.Band[i].Width = GetBEdouble(&src->WSRT.Band[i].Width);
      dst->WSRT.Band[i].SkyMidFreq = GetBEdouble(&src->WSRT.Band[i].SkyMidFreq);
    }

    /* Gen Type */

    dst->gen.NFiles = GetBEint(&src->gen.NFiles);
    dst->gen.FileNum = GetBEint(&src->gen.FileNum);
    dst->gen.NTapes = GetBEint(&src->gen.NTapes);
    dst->gen.TapeNum = GetBEint(&src->gen.TapeNum);
    dst->gen.ParBlkSize = GetBEint(&src->gen.ParBlkSize);
    dst->gen.DataBlkSize = GetBEint(&src->gen.DataBlkSize);
    dst->gen.DataMJD = GetBEint(&src->gen.DataMJD);

    dst->gen.DataTime = GetBEdouble(&src->gen.DataTime);


    /* Obsy */
    dst->obsy.Long = GetBEdouble(&src->obsy.Long);
    dst->obsy.Lat = GetBEdouble(&src->obsy.Lat);
    dst->obsy.Height = GetBEdouble(&src->obsy.Height);


    /* Obs */

    dst->obs.StMJD = GetBEint(&src->obs.StMJD);
    dst->obs.StTime = GetBEint(&src->obs.StTime);

    dst->obs.StLST = GetBEdouble(&src->obs.StLST);
    dst->obs.MaserOffset = GetBEdouble(&src->obs.MaserOffset);
    dst->obs.Dur = GetBEdouble(&src->obs.Dur);
    dst->obs.IonRM = GetBEdouble(&src->obs.IonRM);

    /* Target */

    dst->src.RA =  GetBEdouble(&src->src.RA);
    dst->src.Dec =  GetBEdouble(&src->src.Dec);

    /* Signal Path */
    
    dst->WSRT.ResFringeOff = GetBEint(&src->WSRT.ResFringeOff);
    dst->WSRT.AddBoxOn = GetBEint(&src->WSRT.AddBoxOn);

    /* Gen Type */

    /* INts */
    dst->mode.Nr = GetBEint(&src->mode.Nr);
    dst->mode.FIRFactor = GetBEint(&src->mode.FIRFactor);
    dst->mode. NSHARCsAdded= GetBEint(&src->mode.NSHARCsAdded);
    dst->mode.NSampsAdded = GetBEint(&src->mode.NSampsAdded);
    dst->mode.NFreqInFile = GetBEint(&src->mode.NFreqInFile);
    dst->mode.Tsamp = GetBEint(&src->mode.Tsamp);
    dst->mode.Iout = GetBEint(&src->mode.Iout);
    dst->mode.Qout = GetBEint(&src->mode.Qout);
    dst->mode.Uout = GetBEint(&src->mode.Uout);
    dst->mode.Vout = GetBEint(&src->mode.Vout);
    dst->mode.Xout = GetBEint(&src->mode.Xout);
    dst->mode.Yout = GetBEint(&src->mode.Yout);
    dst->mode.NDMs = GetBEint(&src->mode.NDMs);
    dst->mode.DC_Dynamic = GetBEint(&src->mode.DC_Dynamic);
    dst->mode.ScaleDynamic = GetBEint(&src->mode.ScaleDynamic);
    dst->mode.FloatsOut = GetBEint(&src->mode.FloatsOut);
    dst->mode.BitsPerSamp = GetBEint(&src->mode.BitsPerSamp);

    /* Doubles */

    for (i=0;i<MAXDMS;i++)
      dst->mode.DM[i] = GetBEdouble(&src->mode.DM[i]);

    dst->mode.RM = GetBEdouble(&src->mode.RM);
    dst->mode.AdjustInterval = GetBEdouble(&src->mode.AdjustInterval);
    dst->mode.SigmaRange = GetBEdouble(&src->mode.SigmaRange);

    for (i=0;i<MAXFREQBANDS;i++)
      {
	dst->mode.XPolScaleFac[i] = GetBEfloat(&src->mode.XPolScaleFac[i]);
	dst->mode.ActCluster[i] = GetBEint(&src->mode.ActCluster[i]);
	dst->check.ExitCodeClust[i] = GetBEint(&src->check.ExitCodeClust[i]);
      }

    /* Check */

    dst->check.ExitCodeObsy = GetBEint(&src->check.ExitCodeObsy);
    dst->check.ExitCodePuMa = GetBEint(&src->check.ExitCodePuMa);
    dst->check.ExitCodeDataConsistency = GetBEint(&src->check.ExitCodeDataConsistency);

    /* Reduction */

    dst->redn.Raw = GetBEint(&src->redn.Raw);
    dst->redn.MJDint = GetBEint(&src->redn.MJDint);
    dst->redn.NTimeInts = GetBEint(&src->redn.NTimeInts);
    dst->redn.NFreqs = GetBEint(&src->redn.NFreqs);
    dst->redn.IsDedisp = GetBEint(&src->redn.IsDedisp);
    dst->redn.IsCohDedisp = GetBEint(&src->redn.IsCohDedisp);
    dst->redn.CohFFTSize = GetBEint(&src->redn.CohFFTSize);
    dst->redn.IsPwr = GetBEint(&src->redn.IsPwr);
    dst->redn.NBins = GetBEint(&src->redn.NBins);
    dst->redn.Folded = GetBEint(&src->redn.Folded);
    dst->redn.Polyco = GetBEint(&src->redn.Polyco);
    dst->redn.NCoef = GetBEint(&src->redn.NCoef);
    dst->redn.PolycoSpan = GetBEint(&src->redn.PolycoSpan);
    dst->redn.IsAdjusted = GetBEint(&src->redn.IsAdjusted);
    dst->redn.PolycoStored = GetBEint(&src->redn.PolycoStored);
    dst->redn.Bary = GetBEint(&src->redn.Bary);
    dst->redn.OI = GetBEint(&src->redn.OI);
    dst->redn.OQ = GetBEint(&src->redn.OQ);
    dst->redn.OU = GetBEint(&src->redn.OU);
    dst->redn.OV = GetBEint(&src->redn.OV);
    dst->redn.OP = GetBEint(&src->redn.OP);
    dst->redn.OTheta = GetBEint(&src->redn.OTheta);
    dst->redn.Op = GetBEint(&src->redn.Op);
    dst->redn.Ov = GetBEint(&src->redn.Ov);
    dst->redn.Opoldeg = GetBEint(&src->redn.Opoldeg);
    dst->redn.OX = GetBEint(&src->redn.OX);
    dst->redn.OY = GetBEint(&src->redn.OY);
    dst->redn.TRedn = GetBEint(&src->redn.TRedn);

    dst->redn.MJDfrac = GetBEdouble(&src->redn.MJDfrac);
    dst->redn.DM = GetBEdouble(&src->redn.DM);
    dst->redn.DeltaTime = GetBEdouble(&src->redn.DeltaTime);
    dst->redn.FreqCent = GetBEdouble(&src->redn.FreqCent);
    dst->redn.DeltaFreq = GetBEdouble(&src->redn.DeltaFreq);
    dst->redn.FoldPeriod = GetBEdouble(&src->redn.FoldPeriod);
} 
  

/* Program to read in or write out the header and decide whether or not you
   are on a linux box and thus need to fiddle with the header structure or
   not.  BWS 2002-02-26 */
int prheader(Header_type *inphdr,FILE *srcfile)
{
  int              r;
#if defined ( __linux__) || defined (__alpha)
  Header_type      behdr;
  r = fread(&behdr,sizeof(Header_type),1,srcfile);
  beheader_convert(&behdr,inphdr);
#else
  r = fread(inphdr,sizeof(Header_type),1,srcfile);
#endif

  /* Check to see what version of the header and also what version of the data reduction software 
       we are using -- BWS 2003 - June - 04 */
  /* Checking whether or not we are not yet up to DPC version 1.3 */
  if (strcmp(inphdr->gen.HdrVer,"DPC_1.1")==0 || strcmp(inphdr->gen.HdrVer,"DPC_1.2")==0)
    {
         inphdr->redn.RednSoftwareVer = 0;
    }
   else
   {
       inphdr->redn.RednSoftwareVer = defRednSoftwareVer;
   }
  return r;
}

int pwheader(Header_type *outphdr,FILE *dstfile)
{
  int              r;
#if defined ( __linux__) || defined (__alpha)
  Header_type      behdr;
  beheader_convert(outphdr,&behdr);
  r = fwrite(&behdr,sizeof(Header_type),1,dstfile);
# else
  r = fwrite(outphdr,sizeof(Header_type),1,dstfile);
#endif

  return r;
}

/* Return the number of polarization channels in header */
int puma_nrpol(Header_type hdr)
{
  int NrPol;
  NrPol = 0;
  if((hdr.redn.OX) && (!hdr.redn.OY)) NrPol=1;
  if((!hdr.redn.OX) && (hdr.redn.OY)) NrPol=1;
  if((hdr.redn.OX) && (hdr.redn.OY))  NrPol=2;
  if(hdr.redn.OI) {
      NrPol = 1;
      if(hdr.redn.OQ && hdr.redn.OU && hdr.redn.OV) NrPol=4;
      //      if (hdr.redn.OTheta) NrPol=9;
      // 5 Sept 2016: Patrick changed this to be compatible with POLTYPE_ILVPAdPA, which only has 5 polarization channels
      if(hdr.redn.OP && hdr.redn.OV && hdr.redn.OTheta && hdr.redn.Op && hdr.redn.Ov == 0 && hdr.redn.Opoldeg == 0 && hdr.redn.OU == 0 && hdr.redn.OX == 0 && hdr.redn.OY == 0) NrPol=5;
      // POLTYPE_ILVPAdPATEldEl
      if(hdr.redn.OP && hdr.redn.OV && hdr.redn.OTheta && hdr.redn.Op && hdr.redn.Ov == 0 && hdr.redn.Opoldeg == 0 && hdr.redn.OU == 1 && hdr.redn.OX == 1 && hdr.redn.OY == 1) NrPol=8;
  }
  return NrPol;
}

#define ConvertArrayFromBE(xptr,n, size)   swapArray(xptr, n, size)

void pumaread(void *prptr, int size, int nelem,FILE *in)
{

#ifdef __alpha
  fread(prptr,size,nelem,in);
  ConvertArrayFromBE(prptr,nelem,size);
#endif

#ifdef __linux__
  fread(prptr,size,nelem,in);
  ConvertArrayFromBE(prptr,nelem,size);
#endif

#ifdef __hpux
  fread(prptr,size,nelem,in);
#endif
}

/* Does same as above except now writes out file in bigendian format*/


int pumawrite(void *prptr, int size, int nelem,FILE *out)
{
  int status;

#ifdef __linux__
  ConvertArrayFromBE(prptr,nelem,size);
  status = fwrite(prptr,size,nelem,out);
  ConvertArrayFromBE(prptr,nelem,size);
# else
  status = fwrite(prptr,size,nelem,out);
# endif
  return status;
}


int readWSRTHeader(datafile_definition *datafile, verbose_definition verbose)
{
  int i, version;
  char *dummy, bytevalue;
  int *idummy;
  float *fdummy;
  Header_type puma_hdr;

  prheader(&puma_hdr, datafile->fptr_hdr);
  datafile->NrSubints = puma_hdr.redn.NTimeInts;
  datafile->NrBins = puma_hdr.redn.NBins;
  datafile->NrBits = 8*sizeof(float);
  /* If the data is not folded, it is a dedispersed time series, so
     set nr pulses equal to 1. In header nrpulses=nrbins for some
     reason. Not completely sure if it is also in real data the case,
     but definitely when it is converted. */
  if(puma_hdr.redn.Folded == 0) {
    if(verbose.debug) {
      printf("DEBUG: redn.Folded flag is not set, assume data is searchmode data\n");
    }
    datafile->NrSubints = 1;
    //    datafile->dd_mode = 1;
    datafile->gentype = GENTYPE_SEARCHMODE;
    datafile->isFolded = 0;
    datafile->foldMode = FOLDMODE_UNKNOWN;
    datafile->fixedPeriod = -1;
  }else {
    if(verbose.debug) {
      printf("DEBUG: redn.Folded flag is set\n");
    }
    datafile->isFolded = 1;
    datafile->foldMode = FOLDMODE_FIXEDPERIOD;
    datafile->fixedPeriod = puma_hdr.redn.FoldPeriod;
  }
  datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
  datafile->fixedtsamp = puma_hdr.redn.DeltaTime;
  datafile->tsubMode = TSUBMODE_FIXEDTSUB;
  if(datafile->tsub_list != NULL)
    free(datafile->tsub_list);
  datafile->tsub_list = (double *)malloc(sizeof(double));
  if(datafile->tsub_list == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readWSRTHeader: Memory allocation error");
    return 0;
  }
  datafile->tsub_list[0] = puma_hdr.obs.Dur/(double)datafile->NrSubints;
  datafile->mjd_start = puma_hdr.redn.MJDint + puma_hdr.redn.MJDfrac;
  datafile->freqMode = FREQMODE_UNIFORM;
  if(datafile->freqlabel_list != NULL) {
    free(datafile->freqlabel_list);
    datafile->freqlabel_list = NULL;
  }
  //  datafile->freq_list = malloc(2*sizeof(double));
  //  if(datafile->freq_list == NULL) {
  //    fflush(stdout);
  //    printerror(verbose.debug, "ERROR readWSRTHeader: Memory allocation error.");
  //    return 0;
  //  }
  set_centre_frequency(datafile, puma_hdr.redn.FreqCent, verbose);
  /*  datafile->bw = puma_hdr.redn.DeltaFreq;
      datafile->channelbw = puma_hdr.WSRT.Band[0].Width/(float)puma_hdr.mode.NFreqInFile;*/
  double bw;
  bw = puma_hdr.WSRT.Band[0].Width;
  double channelbw = puma_hdr.redn.DeltaFreq;
  if(channelbw < 0 || bw < 0) {
    //    datafile->channelbw = -fabs(datafile->channelbw);
    bw = -fabs(bw);
  }
  if(set_bandwidth(datafile, bw, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readWSRTHeader: Bandwidth changing failed.");
    return 0;
  }
  datafile->NrFreqChan = puma_hdr.redn.NFreqs;
  datafile->ra = puma_hdr.src.RA;
  if(datafile->ra < 0)
    datafile->ra += 2*M_PI;
  datafile->dec = puma_hdr.src.Dec;
  datafile->dm = puma_hdr.redn.DM;
  datafile->rm = puma_hdr.mode.RM;
  datafile->isDeDisp = puma_hdr.redn.IsDedisp;
  if(set_psrname_PSRData(datafile, puma_hdr.src.Pulsar, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readWSRTHeader: Setting pulsar name failed.");
    return 0;
  }
  if(set_observatory_PSRData(datafile, puma_hdr.obsy.Name, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readWSRTHeader: Setting observatory name failed.");
    return 0;
  }
  //  datafile->telescope_long = puma_hdr.obsy.Long;
  //  datafile->telescope_lat = puma_hdr.obsy.Lat;
  if(verbose.debug) {
    printf("DEBUG: defined long=%lf deg lat=%lf deg height=%lf m\n", puma_hdr.obsy.Long*180.0/M_PI, puma_hdr.obsy.Lat*180.0/M_PI, puma_hdr.obsy.Height);
  }
  // Go from geodetic to geocentric coordinates using slalib
  /*
  double r, z, telc_lat, telc_long;
  csla_geoc(puma_hdr.obsy.Lat, 0*puma_hdr.obsy.Height, &r, &z);
  r *= 149597870700;  // From AU to meters
  z *= 149597870700;
  telc_lat = atan2(z, r);
  telc_long = puma_hdr.obsy.Long;
  if(verbose.debug) {
    printf("DEBUG: r=%lf m z=%lf m long=%lf deg lat=%lf deg\n", r, z, telc_long*180.0/M_PI, telc_lat*180.0/M_PI);
  }
  // Go from cylindrical distance to radial distance
  r = sqrt(r*r+z*z);
  datafile->telescope_X = r*cos(telc_lat)*cos(telc_long);
  datafile->telescope_Y = r*cos(telc_lat)*sin(telc_long);
  datafile->telescope_Z = r*sin(telc_lat);
  */

  // Go from geodetic to geocentric coordinates using tempo2 function
  tempo2_GRS80_to_ITRF(puma_hdr.obsy.Long, puma_hdr.obsy.Lat, puma_hdr.obsy.Height, &(datafile->telescope_X), &(datafile->telescope_Y), &(datafile->telescope_Z));
  if(verbose.debug) {
    printf("DEBUG: derived ITRF (X,Y,Z) m = (%lf, %lf, %lf) m\n", datafile->telescope_X, datafile->telescope_Y, datafile->telescope_Z);
  }
  double telc_lat, telc_long;
  telc_long = observatory_long_geodetic(*datafile);
  telc_lat = observatory_lat_geodetic(*datafile);
  if(verbose.debug) {
    printf("DEBUG: long=%lf deg lat=%lf (should be the header defined position)\n", telc_long*180.0/M_PI, telc_lat*180.0/M_PI);
  }
  // If it looks like this is the WSRT, update position to tempo2 values.
  // Especially the latitude calculation above does not appear to be entirely correct.
  if(telc_long*180.0/M_PI < 6.61 && telc_long*180.0/M_PI > 6.60) {
    if(telc_lat*180.0/M_PI < 53.0 && telc_lat*180.0/M_PI > 52.7) {
      datafile->telescope_X = 3828445.659;
      datafile->telescope_Y = 445223.600000;
      datafile->telescope_Z = 5064921.5677;
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readWSRTHeader: Updated derived ITRF position of telescope with accurate WSRT position.");
    }
  }

  if(set_institute_PSRData(datafile, puma_hdr.gen.Comment, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readWSRTHeader: Setting institute name failed.");
    return 0;
  }
  for(i = 0; i < strlen(datafile->psrname); i++) {
    if(datafile->psrname[i] == ' ')
      datafile->psrname[i] = '_';                    /* Make sure there are no spaces */
  }
  for(i = 0; i < strlen(datafile->institute); i++) {
    if(datafile->institute[i] == ' ')
      datafile->institute[i] = '_';                    /* Make sure there are no spaces */
  }
  if(set_institute_PSRData(datafile, "Converted", verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readWSRTHeader: Setting institute name failed.");
    return 0;
  }
  if(set_instrument_PSRData(datafile, "PuMa", verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readWSRTHeader: Setting instrument name failed.");
    return 0;
  }

  datafile->NrPols = puma_nrpol(puma_hdr);
  if(puma_hdr.redn.OI == TRUE && puma_hdr.redn.OP == TRUE && puma_hdr.redn.OV == TRUE && puma_hdr.redn.OQ == FALSE && puma_hdr.redn.OU == FALSE && puma_hdr.redn.OTheta == FALSE && puma_hdr.redn.Op == FALSE && puma_hdr.redn.Ov == FALSE && puma_hdr.redn.Opoldeg == FALSE && puma_hdr.redn.OX == FALSE && puma_hdr.redn.OY == FALSE) {
    datafile->NrPols = 3;
  }





  dummy = puma_hdr.redn.ReductionDummy;
  i = 1;
  if(dummy[0] != 'p')
    i = 0;
  dummy ++;
  if(dummy[0] != 's')
    i = 0;
  dummy ++;
  if(dummy[0] != 'o')
    i = 0;
  dummy ++;
  if(dummy[0] != 'f')
    i = 0;
  dummy ++;
  if(dummy[0] != 't')
    i = 0;
  dummy ++;

  if(i == 1) {
    idummy = (int *)dummy;
    dummy += sizeof(int);
    version = *idummy;

    if(verbose.debug)
      printf("DEBUG: Reading extended PuMa header version %d\n", version);

    if(version >= 1 && version <= 8) {
      idummy = (int *)dummy;
      dummy += sizeof(int);
      datafile->gentype = *idummy;
      if(verbose.debug) {
	printf("DEBUG: Extended header suggests that gentype=%d\n", datafile->gentype);
      }
      
      if(version >= 5) {   // Integer header parameter stored as a byte to save space
	bytevalue = *dummy;
	dummy += 1;
	datafile->xrangeset = bytevalue;
      }else {
	idummy = (int *)dummy;
	dummy += sizeof(int);
	datafile->xrangeset = *idummy;
      }

      fdummy = (float *)dummy;
      dummy += sizeof(float);
      datafile->xrange[0] = *fdummy;
      
      fdummy = (float *)dummy;
      dummy += sizeof(float);
      datafile->xrange[1] = *fdummy;
      
      if(version >= 5) {   // Integer header parameter stored as a byte to save space
	bytevalue = *dummy;
	dummy += 1;
	datafile->yrangeset = bytevalue;
      }else {
	idummy = (int *)dummy;
	dummy += sizeof(int);
	datafile->yrangeset = *idummy;
      }

      fdummy = (float *)dummy;
      dummy += sizeof(float);
      datafile->yrange[0] = *fdummy;
      
      fdummy = (float *)dummy;
      dummy += sizeof(float);
      datafile->yrange[1] = *fdummy;

      if(version >= 2) {
	if(version >= 5) {   // Integer header parameter stored as a byte to save space
	  bytevalue = *dummy;
	  dummy += 1;
	  datafile->poltype = bytevalue;
	}else {
	  idummy = (int *)dummy;
	  dummy += sizeof(int);
	  datafile->poltype = *idummy;
	}

	if(version >= 5) {   // Integer header parameter stored as a byte to save space
	  bytevalue = *dummy;
	  dummy += 1;
	  datafile->feedtype = bytevalue;
	}else {
	  idummy = (int *)dummy;
	  dummy += sizeof(int);
	  datafile->feedtype = *idummy;
	}
      }

      if(version >= 3) {
	if(version >= 5) {   // Integer header parameter stored as a byte to save space
	  bytevalue = *dummy;
	  dummy += 1;
	  datafile->isDeFarad = bytevalue;
	}else {
	  idummy = (int *)dummy;
	  dummy += sizeof(int);
	  datafile->isDeFarad = *idummy;
	}
      }

      if(version >= 4) {
	if(version >= 5) {   // Integer header parameter stored as a byte to save space
	  bytevalue = *dummy;
	  dummy += 1;
	  datafile->isDePar = bytevalue;
	}else {
	  idummy = (int *)dummy;
	  dummy += sizeof(int);
	  datafile->isDePar = *idummy;
	}
      }

      if(version >= 5) {
	bytevalue = *dummy;
	dummy += 1;
	datafile->cableSwap = bytevalue;

	bytevalue = *dummy;
	dummy += 1;
	datafile->cableSwapcor = bytevalue;

	//	printf("XXXXX %d %d total = %ld\n", datafile->cableSwap, datafile->cableSwapcor, dummy - puma_hdr.redn.ReductionDummy);
      }

      if(version >= 6) {
	fdummy = (float *)dummy;
	dummy += sizeof(float);
	datafile->telescope_X = *fdummy;

	fdummy = (float *)dummy;
	dummy += sizeof(float);
	datafile->telescope_Y = *fdummy;

	fdummy = (float *)dummy;
	dummy += sizeof(float);
	datafile->telescope_Z = *fdummy;
      }

      if(version >= 7) {
	fdummy = (float *)dummy;
	dummy += sizeof(float);
	datafile->freq_ref = *fdummy;
      }

      if(version >= 8) {
	bytevalue = *dummy;
	dummy += 1;
	datafile->isDebase = bytevalue;
      }
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readWSRTHeader: Undefined version of the extended header (%d).", version);
      return 0;
    }
  }

  /*
// 5 sept 2016 - Patrick: don't see how this can be true. Must have become obsolete?
  if(datafile->poltype == POLTYPE_ILVPAdPA) {
    if(datafile->NrPols != 3) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readWSRTHeader: For gentype ILVPA the nr of polarization is expected to be (wrongly) set to be 3 in the header, it was %ld. The number of polarizations channels will be set to 5. If incorrect, expect things to go wrong (badly).", datafile->NrPols);
      datafile->NrPols = 5;
    }
    datafile->NrPols = 5;
  }
  */

  return 1;
}


/* 
obsID                   e.g. 200605031
timefilenr              starting from 0
band                    0-7
NFiles                  Number of files over which this obs. is spread starting with 0.
FileNum                 Which file out of NFiles is this               
nrTimeSamples           In this file
DataMJD                 MJD of first sample of this data block  
DataTime                Time of first sample (fraction of day) 
ObsName                 e.g. S06B/001/J1808-0813
StartMJD                MJD at start of observation   
StTime                  Starttime (s after midnight, multiple of 10 s)
MaserOffset             (Averaged over obs.) offset of Maser (in s)  
Dur                     Actual duration of observation (in s)   
pulsar                  e.g. PSR_B1758-03
RA                      RA of the target (in radians)   
Dec                     Dec of the target (in radians)     
Epoch                   e.g. J2000
bands                   Information about eight frequency bands
                          Boolean       Operational;                   Is band o.k?                                      
                          Boolean       NonFlip;                       Is band flipped (reverse freq order within band)? 
                          double        MidFreq;                       Mid frequency of band (in MHz): 5.0 or 6.25       
                          double        Width;                         Width of the band (in MHz): 2.5, 5.0 or 10.0      
                          double        SkyMidFreq;                    Mid sky frequency of band (in MHz)                
BandsToClusMap          4,3,2,1,0,7,6,5    Band 0 is connected to cluster 4       
XPolScaleFac            Multiplication factor of X pol realtive to Y  for each band

FIRFactor;              Mode 0,4 ; 1,2,4,8,...     Mode -1,1,2 : 1     
NSHARCsAdded;           Mode 1,2 ; 1,2,3,6         Mode -1,0,4 : 1     
NSampsAdded;            Mode 1,2 ; 1,2,4,8,16,.. depends on FreqChans  
NFreqInFile;            Mode 0-4 ; 1,2,4,8,16,..FreqChans in this file 
			         Normally one cluster per file is written    
Tsamp;                  Mode 0-4 ; output sample interval in nano sec  
polmode                 1 = I
                        2 = XY
                        4 = IQUV
AdjustInterval          Adjustment interval       (e.g.  40 seconds)
BitsPerSamp             Number of output bits ; 1,2,4,8 
SigmaRange              Number of Sigmas mapped onto output range    
NrFrames                Nr of frames stored in puma file per polarization
*/

/*void FillPuMaHeader(Header_type *hdr, int obsID, int timefilenr, int freqband, int NFiles, int FileNum, long nrTimeSamples, int DataMJD, double DataTime, int StartMJD, char *ObsName, int StTime, double MaserOffset, double Dur, char *pulsar, double RA, double Dec, char *Epoch, Band_type *bands, int *BandsToClusMap, float *XPolScaleFac, int FIRFactor, int NSHARCsAdded, int NSampsAdded, int NFreqInFile, int Tsamp, int polmode, double AdjustInterval, int BitsPerSamp, double SigmaRange, long NrFrames, char *observatoryname) */

void FillPuMaHeader(Header_type *hdr, int obsID, int timefilenr, int freqband, int NFiles, int FileNum, int nrTimeSamples, int DataMJD, double DataTime, int StartMJD, char *ObsName, int StTime, double MaserOffset, double Dur, char *pulsar, double RA, double Dec, char *Epoch, Band_type *bands, int *BandsToClusMap, float *XPolScaleFac, int FIRFactor, int NSHARCsAdded, int NSampsAdded, int NFreqInFile, int Tsamp, int polmode, double AdjustInterval, int BitsPerSamp, double SigmaRange, int NrFrames, char *observatoryname, double longitude, double latitude, double height)
{
  int i;

  /* Clear memory */
  memset(hdr, 0, sizeof(Header_type));
  /*
  Header_type
  ===========
  General_type       gen;
  Observatory_type   obsy;
  Observation_type   obs;
  Target_type        src;
  Signalpath_type    WSRT;
  Mode_type          mode;
  Software_type      software;
  Check_type         check;
  Reduction_type     redn;

  */

  /* General_type */
  strcpy(hdr->gen.HdrVer, "DPC_1.3");
  strcpy(hdr->gen.Platform, "PuMaISim");              /* "PuMa Crate 0"  */
  sprintf(hdr->gen.ThisFileName, "%d.%05d.%d.puma", obsID, timefilenr, freqband);
  sprintf(hdr->gen.ScanNum, "%d", obsID);
  strcpy(hdr->gen.Comment, "PuMaISim");
  hdr->gen.NFiles=NFiles;
  hdr->gen.FileNum=FileNum;
  hdr->gen.TapeID[0] = 0;
  hdr->gen.NTapes = 0;
  hdr->gen.TapeNum = 0;
  hdr->gen.ParBlkSize = NrFrames*polmode*sizeof(Adjustments);                                   
  hdr->gen.DataBlkSize = nrTimeSamples*BitsPerSamp*polmode*NFreqInFile/8;
  hdr->gen.Cluster[0] = FALSE;
  hdr->gen.Cluster[1] = FALSE;
  hdr->gen.Cluster[2] = FALSE;
  hdr->gen.Cluster[3] = FALSE;
  hdr->gen.Cluster[4] = FALSE;
  hdr->gen.Cluster[5] = FALSE;
  hdr->gen.Cluster[6] = FALSE;
  hdr->gen.Cluster[7] = FALSE;
  hdr->gen.Cluster[freqband] = TRUE;
  hdr->gen.DataMJD = DataMJD;
  hdr->gen.DataTime = DataTime;  

  /* Observatory_type */
  strcpy(hdr->obsy.Name, observatoryname);
  hdr->obsy.Long = longitude;
  hdr->obsy.Lat = latitude;
  hdr->obsy.Height = height;

  /* Observation_type */
  strcpy(hdr->obs.ObsName, ObsName);
  strcpy(hdr->obs.ObsType, "observation");          
  strcpy(hdr->obs.LastCalID, "CalId");       
  hdr->obs.StMJD = StartMJD;                    
  hdr->obs.StTime = StTime;                   
  hdr->obs.StLST = 0;                                            /* Not set????? */
  hdr->obs.MaserOffset = MaserOffset;              
  hdr->obs.Dur = Dur;                      
  hdr->obs.IonRM = 0;                                            /* Not set????? */

  /* Target_type */
  strcpy(hdr->src.Pulsar, pulsar);
  i = -1;
  do{
    i++;
    if(i == NAMELEN)
      break;
    if(hdr->src.Pulsar[i] == ' ')
      hdr->src.Pulsar[i] = '_';
  }while(hdr->src.Pulsar[i] != 0);
  hdr->src.RA = RA;            
  hdr->src.Dec = Dec;           
  strcpy(hdr->src.Epoch, Epoch);  

  /* Signalpath_type */
  for(i = 0; i < MAXTELESCOPES; i++) {                          /* Don't set: Active/FrontendID+Status/Syst.Temp */
    hdr->WSRT.Tel[i].Active = FALSE;
    hdr->WSRT.Tel[i].NoiseSrcOff = FALSE;
    hdr->WSRT.Tel[i].FrontendID[0] = 0;
    hdr->WSRT.Tel[i].FrontendStatus[0] = 0;
    hdr->WSRT.Tel[i].Tsys = 0;
  }
  for(i = 0; i < MAXFREQBANDS; i++) {
    memcpy(&(hdr->WSRT.Band[i].Operational), &(bands[i]), sizeof(Band_type));
  }
  hdr->WSRT.Backend[0] = 0;                         /* DCB/DZB/DLB/DXB */
  hdr->WSRT.ResFringeOff = TRUE;                    /* Is residual fringe stopping off? TRUE is needed */
  hdr->WSRT.AddBoxOn = TRUE;                        /* Is the adding box working? TRUE is needed */
  for(i = 0; i < MAXFREQBANDS; i++) {
    hdr->WSRT.BandsToClusMap[i] = BandsToClusMap[i];
  }

  /* Mode_type  */
  hdr->mode.Nr = 1;                                       /* Mode of PuMa observation -1,0,1,2,4  */
  for(i = 0; i < MAXFREQBANDS; i++) {
    hdr->mode.XPolScaleFac[i] = XPolScaleFac[i];          /* Multiplication factor of X pol realtive to Y */
    hdr->mode.ActCluster[i] = bands[i].Operational;       /* BandOperational.and.AdcOk.and.ClusterOk  */
  }
  hdr->mode.FIRFactor = FIRFactor;
  hdr->mode.NSHARCsAdded = NSHARCsAdded;
  hdr->mode.NSampsAdded = NSampsAdded; 
  hdr->mode.NFreqInFile = NFreqInFile;
  hdr->mode.Tsamp = Tsamp;
  hdr->mode.Iout = FALSE;
  hdr->mode.Qout = FALSE;
  hdr->mode.Vout = FALSE;
  hdr->mode.Uout = FALSE;
  hdr->mode.Xout = FALSE;
  hdr->mode.Yout = FALSE;
  if(polmode == 1)
    hdr->mode.Iout = TRUE;
  else if(polmode == 2) {
    hdr->mode.Xout = TRUE;
    hdr->mode.Yout = TRUE;
  }else if(polmode == 4) {
    hdr->mode.Iout = TRUE;
    hdr->mode.Qout = TRUE;
    hdr->mode.Uout = TRUE;
    hdr->mode.Vout = TRUE;
  }
  hdr->mode.NDMs = 0;        
  for(i = 0; i < MAXDMS; i++)
    hdr->mode.DM[i] = 0.000000;
  hdr->mode.RM=0;     
  hdr->mode.DC_Dynamic = TRUE;                   /* Is DC dynamic removed?                         */
  hdr->mode.ScaleDynamic = TRUE;                 /* Is scaling dynamic adjusted                    */
  hdr->mode.AdjustInterval = AdjustInterval;       /* Adjustment interval                            */
  hdr->mode.FloatsOut=FALSE;                     /* Is the data written in floats ?                */
  /*  fprintf(stderr, "!!!!! Writing %d bits per sample (%d)\n", hdr->mode.BitsPerSamp, BitsPerSamp); */
  hdr->mode.BitsPerSamp = BitsPerSamp;           /* Number of output bits ; 1,2,4,8                */
  /*  fprintf(stderr, "!!!!! Writing %d bits per sample (%d)\n", hdr->mode.BitsPerSamp, BitsPerSamp); */
  
  /*  fprintf(stderr, "!!!!! Boolean = %d bytes\n", sizeof(Boolean)); 
  fprintf(stderr, "!!!!! int = %d bytes\n", sizeof(int)); 
  fprintf(stderr, "!!!!! short int = %d bytes\n", sizeof(short int)); 
  fprintf(stderr, "!!!!! long int = %d bytes\n", sizeof(long int)); 
  fprintf(stderr, "!!!!! long = %d bytes\n", sizeof(long)); 
  fprintf(stderr, "!!!!! double = %d bytes\n", sizeof(double)); 
  fprintf(stderr, "!!!!! float = %d bytes\n", sizeof(float)); 
  fprintf(stderr, "!!!!! Dummy = %d bytes\n", sizeof(Dummy)); 
  fprintf(stderr, "!!!!! Header_type = %d bytes\n", sizeof(Header_type)); */

  hdr->mode.SigmaRange = SigmaRange;             /* Number of Sigmas mapped onto output range      */


  /* Software type */
  strcpy(hdr->software.Master, "Created by PuMaISim");    /*   Unique RCS number of Mastersoftware on HP      */
  for(i = 0; i < MAXDSPSPERBOARD; i++)
    hdr->software.DSP[i][0] = 0;                            /* RCS number of DSP software                     */
  strcpy(hdr->software.Filter, "Created by PuMaISim");     /* RCS number decimation filter softw., Mode 0,4  */


  /* Check_type */
  hdr->check.ExitCodeObsy = 0;
  hdr->check.ExitCodePuMa = 0;
  for(i = 0; i < MAXFREQBANDS; i++)
    hdr->check.ExitCodeClust[i] = 0;
  hdr->check.ExitCodeDataConsistency = 0;

}


void beadj_convert_write(Adjustments outpadj, FILE *fout)
{ 
  //  int              r;
  //  Adjustments      beadj;

#ifdef __linux__
  Adjustments      beadj;
  PutBEint(&beadj.framenumber, outpadj.framenumber);
  PutBEfloat(&beadj.scale, outpadj.scale);
  PutBEfloat(&beadj.offset, outpadj.offset);
  //  r = 
  fwrite(&beadj,sizeof(Adjustments),1,fout);
#endif

#ifdef __alpha
  Adjustments      beadj;
  PutBEint(&beadj.framenumber, outpadj.framenumber);
  PutBEfloat(&beadj.scale, outpadj.scale);
  PutBEfloat(&beadj.offset, outpadj.offset);
  //  r = 
  fwrite(&beadj,sizeof(Adjustments),1,fout);
#endif

#ifdef __hpux
  //  r = 
  fwrite(&outpadj,sizeof(Adjustments),1,fout);
#endif
}


/* Assume J2000 
 Number of Sigmas mapped onto output range = 3 (no idea if it is used????)
*/
void FillPuMaHeaderSimple(Header_type *hdr, int obsID, int timefilenr, int freqband, int flipped, float BW, float SkyMidFreq, long nrTimeSamples, double DataMJD, double StartMJD, char *ObsName, double Dur, char *pulsar, double RA, double Dec, int NSHARCsAdded, int NSampsAdded, int NFreqInFile, int Tsamp, int polmode, double AdjustInterval, int BitsPerSamp, long NrFrames, char *observatoryname, double longitude, double latitude, double height)
{
  int DataMJDint, StartMJDint, i, BandsToClusMap[8];
  float XPolScaleFac[8];
  char Epoch[NAMELEN];
  Band_type bands[8];
  Boolean NonFlip;


  DataMJDint = DataMJD;
  DataMJD -= DataMJDint;
  StartMJDint = StartMJD;
  StartMJD -= StartMJDint;
  StartMJD *= 24*3600;               /* In sec after midnight */
  strcpy(Epoch, "J2006");
  if(flipped == 0)
    NonFlip = TRUE;
  else
    NonFlip = FALSE;
  for(i = 0; i < 8; i++) {
    bands[i].Operational = FALSE; 
    bands[i].NonFlip = NonFlip;
    bands[i].SkyMidFreq = SkyMidFreq;             
    bands[i].Width = BW;
    bands[i].MidFreq = 1.25;         /* No idea, probably not used */
  }
  bands[freqband].Operational = TRUE; 
  BandsToClusMap[0] = 0;
  BandsToClusMap[1] = 1;
  BandsToClusMap[2] = 2;
  BandsToClusMap[3] = 3;
  BandsToClusMap[4] = 4;
  BandsToClusMap[5] = 5;
  BandsToClusMap[6] = 6;
  BandsToClusMap[7] = 7;
  XPolScaleFac[0] = 1;
  XPolScaleFac[1] = 1;
  XPolScaleFac[2] = 1;
  XPolScaleFac[3] = 1;
  XPolScaleFac[4] = 1;
  XPolScaleFac[5] = 1;
  XPolScaleFac[6] = 1;
  XPolScaleFac[7] = 1;

  FillPuMaHeader(hdr,  obsID,  timefilenr,  freqband,  1,  1,  nrTimeSamples,  DataMJDint,  DataMJD,  StartMJDint, ObsName,  StartMJD,  0,  Dur, pulsar,  RA,  Dec, Epoch, bands, BandsToClusMap, XPolScaleFac,  0,  NSHARCsAdded,  NSampsAdded,  NFreqInFile,  Tsamp,  polmode,  AdjustInterval,  BitsPerSamp,  3,  NrFrames, observatoryname, longitude, latitude, height);
}

int writeWSRTHeader(datafile_definition datafile, verbose_definition verbose)
{
  Header_type puma_hdr;
  float *fdummy;
  int *idummy;
  char *dummy;

  if(datafile.freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeWSRTHeader: Writing this data format is only implemented when the frequency channels are uniformely separated.");
    return 0;
  }


  FillPuMaHeaderSimple(&puma_hdr, 200700000, 1, 0, 0, get_bandwidth(datafile, verbose), get_centre_frequency(datafile, verbose), datafile.NrBins, datafile.mjd_start, datafile.mjd_start, "", get_tobs(datafile, verbose), datafile.psrname, datafile.ra, datafile.dec, 1, 1, datafile.NrFreqChan, get_tsamp(datafile, 0, verbose), datafile.NrPols, 0, 8, 1, datafile.observatory, observatory_long_geodetic(datafile), observatory_lat_geodetic(datafile), observatory_height_geodetic(datafile));

  puma_hdr.gen.ParBlkSize = 0;
  strncpy(puma_hdr.gen.ScanNum, datafile.scanID, NAMELEN);
  puma_hdr.redn.NBins = datafile.NrBins;
  double period;
  int ret;
  if(datafile.isFolded) {
    ret = get_period(datafile, 0, &period, verbose);
    if(ret == 2) {
      printerror(verbose.debug, "ERROR writeWSRTHeader (%s): Cannot obtain period", datafile.filename);
      return 0;
    }
  }else {
    ret = 1;
    period = -1;
  }
  
  if(datafile.gentype == GENTYPE_SEARCHMODE || datafile.isFolded != 1 || ret == 1 || period < 0) {
    //Write out data as single subint with many bins
    puma_hdr.redn.NTimeInts = datafile.NrBins*datafile.NrSubints;
    puma_hdr.redn.NBins = puma_hdr.redn.NTimeInts;
    puma_hdr.redn.GenType = 0;
    puma_hdr.redn.Folded = 0;
  }else {
    puma_hdr.redn.NTimeInts = datafile.NrSubints;
    puma_hdr.redn.GenType = 4;
    puma_hdr.redn.Folded = 1;
  }
  puma_hdr.redn.FoldPeriod = datafile.fixedPeriod;
  puma_hdr.redn.DeltaTime = get_tsamp(datafile, 0, verbose);
  puma_hdr.redn.MJDint = (int)datafile.mjd_start;
  puma_hdr.redn.MJDfrac = datafile.mjd_start - puma_hdr.redn.MJDint;
  puma_hdr.redn.FreqCent =  get_centre_frequency(datafile, verbose);
  double chanbw;
  if(get_channelbandwidth(datafile, &chanbw, verbose) == 0) {
    printerror(verbose.debug, "ERROR writeWSRTHeader (%s): Cannot obtain channel bandwidth.", datafile.filename);
    return 0;
  }
  puma_hdr.redn.DeltaFreq = chanbw;
  puma_hdr.redn.DM = datafile.dm;
  puma_hdr.mode.RM = datafile.rm;
  puma_hdr.redn.NFreqs = datafile.NrFreqChan;
  puma_hdr.redn.Raw = FALSE;
  puma_hdr.redn.IsDedisp = datafile.isDeDisp;

  if(datafile.NrPols == 1) {
    puma_hdr.redn.OI = TRUE;
  }else if(datafile.NrPols == 2) {
    puma_hdr.redn.OX = TRUE;
    puma_hdr.redn.OY = TRUE;
  }else if(datafile.NrPols == 4) {
    puma_hdr.redn.OI = TRUE;
    puma_hdr.redn.OQ = TRUE;
    puma_hdr.redn.OU = TRUE;
    puma_hdr.redn.OV = TRUE;
  }else if(datafile.NrPols == 3) {
    puma_hdr.redn.OI = TRUE;
    puma_hdr.redn.OP = TRUE;
    puma_hdr.redn.OV = TRUE;
  }else if(datafile.NrPols == 5 && datafile.poltype == POLTYPE_ILVPAdPA) {
    puma_hdr.redn.OI = TRUE;
    puma_hdr.redn.OP = TRUE;
    puma_hdr.redn.OV = TRUE;
    puma_hdr.redn.OTheta = TRUE;
    puma_hdr.redn.Op = TRUE;  // Should be error on PA, but that is not defined in puma....
    puma_hdr.redn.OU = FALSE;  // Following should be false in order to distinguish this from POLTYPE_ILVPAdPATEldEl
    puma_hdr.redn.OX = FALSE;
    puma_hdr.redn.OY = FALSE;
  }else if(datafile.NrPols == 8 && datafile.poltype == POLTYPE_ILVPAdPATEldEl) {
    puma_hdr.redn.OI = TRUE;
    puma_hdr.redn.OP = TRUE;
    puma_hdr.redn.OV = TRUE;
    puma_hdr.redn.OTheta = TRUE;
    puma_hdr.redn.Op = TRUE;  // Should be error on PA, but that is not defined in puma....
    puma_hdr.redn.OU = TRUE;  // Following should be false in order to distinguish this from POLTYPE_ILVPAdPATEldEl
    puma_hdr.redn.OX = TRUE;
    puma_hdr.redn.OY = TRUE;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeWSRTHeader: Cannot write out %ld polarization channels", datafile.NrPols);
    return 0;
  }

  dummy = puma_hdr.redn.ReductionDummy;
  strcpy(dummy, "psoft");
  dummy += 5;

  idummy = (int *)dummy;
  dummy += sizeof(int);
  *idummy = 8;                /* The version of the extended header */

  idummy = (int *)dummy;
  dummy += sizeof(int);
  *idummy = datafile.gentype;

  *dummy = datafile.xrangeset;   // Write out as single byte to save space
  dummy += 1;

  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.xrange[0];

  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.xrange[1];

  *dummy = datafile.yrangeset;   // Write out as single byte to save space
  dummy += 1;

  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.yrange[0];

  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.yrange[1];

  *dummy = datafile.poltype;   // Write out as single byte to save space
  dummy += 1;

  *dummy = datafile.feedtype;   // Write out as single byte to save space
  dummy += 1;

  *dummy = datafile.isDeFarad;   // Write out as single byte to save space
  dummy += 1;

  *dummy = datafile.isDePar;   // Write out as single byte to save space
  dummy += 1;

  *dummy = datafile.cableSwap;   // Write out as single byte to save space
  dummy += 1;

  *dummy = datafile.cableSwapcor;   // Write out as single byte to save space
  dummy += 1;

  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.telescope_X;

  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.telescope_Y;

  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.telescope_Z;

  fdummy = (float *)dummy;
  dummy += sizeof(float);
  *fdummy = datafile.freq_ref;

  *dummy = datafile.isDebase;   // Write out as single byte to save space
  dummy += 1;

  //Total shouldn't exceed 56
  // Version 8 is at 54, so 2 bytes left...
  //  printf("XXXXX TOTAL SHOULD NOT EXCEED 56 - %d %d. tot=%ld\n", datafile.cableSwap, datafile.cableSwapcor, dummy-puma_hdr.redn.ReductionDummy);

  rewind(datafile.fptr_hdr);
  pwheader(&puma_hdr, datafile.fptr_hdr);

  return 1;
}



long long PuMaFilepos(long long PulseNumber, long long FreqChan, long long polarization, long long binnr, long long NrSubints, long long NrBins, long long NrFreqChan, long long data_start)
{
  long long puma_filepos;
  puma_filepos = polarization*NrFreqChan*NrSubints*NrBins;  /* Begin polarization block */
  puma_filepos += FreqChan*NrSubints*NrBins;                /* Begin frequency channel */
  puma_filepos += PulseNumber*NrBins;                      /* Begin first pulse */
  puma_filepos += binnr;
  puma_filepos *= sizeof(float);
  puma_filepos += data_start;
  return puma_filepos;
}

int readPulseWSRTData(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse)
{
  long long filepos;

  filepos = polarization*datafile.NrFreqChan*datafile.NrSubints*datafile.NrBins;  /* Begin polarization block */
  filepos += freq*datafile.NrSubints*datafile.NrBins;                /* Begin frequency channel */
  filepos += pulsenr*datafile.NrBins;                      /* Begin first pulse */
  filepos += binnr;
  filepos *= sizeof(float);
  filepos += datafile.datastart;
  fseeko(datafile.fptr, filepos, SEEK_SET);
  pumaread(pulse, sizeof(float), nrSamples, datafile.fptr);
  return 1;
}

int writePulseWSRTData(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse)
{
  long long filepos;

  filepos = polarization*datafile.NrFreqChan*datafile.NrSubints*datafile.NrBins;  /* Begin polarization block */
  filepos += freq*datafile.NrSubints*datafile.NrBins;                /* Begin frequency channel */
  filepos += pulsenr*datafile.NrBins;                      /* Begin first pulse */
  filepos += binnr;
  filepos *= sizeof(float);
  filepos += datafile.datastart;
  fseeko(datafile.fptr, filepos, SEEK_SET);
  pumawrite(pulse, sizeof(float), nrSamples, datafile.fptr);
  return 1;
}


int writePuMafile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  long n, f, p;
  fseeko(datafile.fptr, datafile.datastart, SEEK_SET);
  for(p = 0; p < datafile.NrPols; p++) {
    for(f = 0; f < datafile.NrFreqChan; f++) {
      for(n = 0; n < datafile.NrSubints; n++) {
	if(verbose.verbose && verbose.nocounters == 0) printf("writePuMafile: pulse %ld/%ld\r", n+1, datafile.NrSubints);
	pumawrite(&data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))], sizeof(float), datafile.NrBins, datafile.fptr);
      }
    }
  }
  if(verbose.verbose) printf("  Writing is done.              \n");
  return 1;
}

int readPuMafile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  long n, f, p;
  if(verbose.verbose) {
    printf("Start reading PuMa file\n");
  }
  fseeko(datafile.fptr, datafile.datastart, SEEK_SET);
  for(p = 0; p < datafile.NrPols; p++) {
    for(f = 0; f < datafile.NrFreqChan; f++) {
      for(n = 0; n < datafile.NrSubints; n++) {
	if(verbose.verbose && verbose.nocounters == 0) 
	  printf("  Progress reading PuMa file (%.1f%%)\r", 100.0*(n+(f+p*datafile.NrFreqChan)*datafile.NrSubints)/(float)(datafile.NrSubints*datafile.NrFreqChan*datafile.NrPols));
	pumaread(&data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))], sizeof(float), datafile.NrBins, datafile.fptr);
      }
    }
  }
  if(verbose.verbose) printf("  Reading is done.                           \n");
  return 1;
}



/* Write out the last history line to PuMa file. 
   Returns 1 on success, 0 on error */
int writeHistoryPuma(datafile_definition datafile, verbose_definition verbose)
{
  Header_type puma_hdr;
  char *txt, *txt_date, *user, *hostname, *txt_cmd;
  char questionmark[2];
  datafile_history_entry_definition *curHistoryEntry;
  questionmark[0] = '?';
  questionmark[1] = 0;
  txt = malloc(10000);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "Error writeHistoryPuma: Memory allocation error");
    return 0;
  }

  curHistoryEntry = &(datafile.history);
  // Find last entry
  while(curHistoryEntry->nextEntry != NULL) {
    curHistoryEntry = curHistoryEntry->nextEntry;
  }
  if(curHistoryEntry->timestamp != NULL)
    txt_date = curHistoryEntry->timestamp;
  else
    txt_date = questionmark;
  if(curHistoryEntry->user != NULL)
    user = curHistoryEntry->user;
  else
    user = questionmark;
  if(curHistoryEntry->hostname != NULL)
    hostname = curHistoryEntry->hostname;
  else
    hostname = questionmark;
  if(curHistoryEntry->cmd != NULL)
    txt_cmd = curHistoryEntry->cmd;
  else
    txt_cmd = questionmark;

  if(datafile.opened_flag) {
    /*
      fflush(datafile.fptr_hdr);
      fseeko(datafile.fptr_hdr, 0, SEEK_SET);
      fflush(datafile.fptr_hdr);
      prheader(&puma_hdr, datafile.fptr_hdr);
      fflush(datafile.fptr_hdr);
    */
    sprintf(txt, "%s %s@%s: %s", txt_date, user, hostname, txt_cmd);
    /*  fprintf(stderr, "Appending '%s' to history\n", txt); */
    memset(puma_hdr.redn.Command, 0, CMDLEN);
    strncpy(puma_hdr.redn.Command, txt, CMDLEN-1);
    
    /*  printf("XXXX offset %ld\n", puma_hdr.redn.Command-(char *)(&puma_hdr)); */
    fseeko(datafile.fptr_hdr, puma_hdr.redn.Command-(char *)(&puma_hdr), SEEK_SET);
    fwrite(puma_hdr.redn.Command, 1, CMDLEN, datafile.fptr_hdr);
    
    /*
      fseeko(datafile.fptr_hdr, 0, SEEK_SET);
      fflush(datafile.fptr_hdr);
      fprintf(stderr, "XXXX '%s'\n", puma_hdr.gen.HdrVer);
      pwheader(&puma_hdr, datafile.fptr_hdr); 
      fprintf(stderr, "XXXX '%s'\n", puma_hdr.gen.HdrVer);
      fflush(datafile.fptr_hdr);
      fprintf(stderr, "XXXX %Ld\n", datafile.datastart);
    */
    
    
    fseeko(datafile.fptr, datafile.datastart, SEEK_SET);
    fseeko(datafile.fptr_hdr, datafile.datastart, SEEK_SET);
  }else {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING writeHistoryPuma: File not opened, ignoring command");
  }
  free(txt);
  return 1;
}

/* Reads the history. 
   Returns 1 on success, 0 on error */
int readHistoryPuma(datafile_definition *datafile, verbose_definition verbose)
{
  Header_type puma_hdr;
  datafile_history_entry_definition *curHistoryEntry;

  if(datafile->opened_flag == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHistoryPuma: File already closed?");
    return 0;
  }
  fseeko(datafile->fptr_hdr, 0, SEEK_SET);
  prheader(&puma_hdr, datafile->fptr_hdr);

  curHistoryEntry = &(datafile->history);
  curHistoryEntry->timestamp = NULL;
  curHistoryEntry->cmd = NULL;
  curHistoryEntry->user = NULL;
  curHistoryEntry->hostname = NULL;
  curHistoryEntry->nextEntry = NULL;

  curHistoryEntry->cmd = malloc(strlen(puma_hdr.redn.Command)+1);
  if(curHistoryEntry->cmd == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHistoryPuma: Memory allocation error"); 
    return 0;
  }
  strcpy(curHistoryEntry->cmd, puma_hdr.redn.Command);

  //  printf("  %s\n", puma_hdr.redn.Command);

  return 1;
}
