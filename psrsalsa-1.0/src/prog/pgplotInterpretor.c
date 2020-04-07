#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <errno.h>
#include "cpgplot.h"
#include "psrsalsa.h"

int getfloats(FILE *fin, float *f1, int c1, float *f2, int c2, float *f3, int c3);
int getstrings(char *txt, int nrstrings, char *string1, char *string2, char *string3);

int main(int argc, char **argv)
{
  char txt[1000];
  FILE *fin;
  int i, continuerunning, debug, ret;
  fin = stdin;

  debug = 0;

  if(argc > 1) {
    for(i = 1; i < argc; i++) {
      if(strcmp(argv[i], "-i") == 0) {
	fin = fopen(argv[i+1], "r");
	if(fin == NULL) {
	  sprintf(txt, "Cannot open %s ", argv[i+1]);
	  perror(txt);
	  return 0;
	}
	i += 1;
      }else if(strcmp(argv[i], "-h") == 0) {
	printf("-i filename      open script\n");
	return 0;
      }else {
	fprintf(stderr, "Unknown option '%s'\n", argv[i]);
      }
    }
  }

  if(fin == stdin) {
    printf("Use pgplotInterpretor -i file to load a file with commands\nHit ctrl-D to exit and use 'help' for help.\n");
  }
  continuerunning = 1;
  do {
    if(fgets(txt, 1000, fin) == NULL) {
      printf("Reached EOF\n");
      continuerunning = 0;
    }else {
      ret = doPgplotInterprator(txt, debug);
      if(ret == 0 || ret == 100) {
	return 0;
      }
    }
  }while(continuerunning);
  
  fclose(fin);
  return 0;
}

