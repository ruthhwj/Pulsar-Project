//START REGION RELEASE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <pwd.h>
#include <unistd.h>
#include <netdb.h>
//#ifdef __APPLE__
//#include <curses.h>
//#else
//#include <termio.h>
#include <termios.h>
#include <sys/ioctl.h>
//#endif
#include "psrsalsa.h"
#include <fcntl.h>
#include <signal.h>
#include <stdarg.h>
#include <errno.h>
#include <math.h>
//START REGION DEVELOP
#include <arpa/inet.h>

//START REGION RELEASE

/* Define own getch function, which doesn't exist in all versions of
   C. */
int pgetch(void)
{
  //  int ch;
  //  ch = getchar();  // Waits for return after a line befor characters are read in....

  //  char ch;
  //  fread(&ch, 1, 1, stdin);  // Waits for return after a line befor characters are read in....

  //  char ch;
  //  fscanf(stdin, "%c", &ch);  // Waits for return after a line befor characters are read in....

  /*
    // The following requires the curses library
#ifdef __APPLE__
  int ch;
  initscr();     // Initialise the curses library
  cbreak();      // Don't wait for return
  while((ch = getch()) == ERR) {   // ERR is returned if no key is available
    usleep(1000);
  }
  nocbreak();
  endwin();      // Quit the curses library stuff. If not done, all the printf's seem to be garbled since \n goes to new line, but not to start of line. Not sure how to fix this. Still the initial text of terminal appears to disappear upon first call of pgetch(). 
#else
  */


  /*
//    This depends on termio.h, which doesn't seem to be supported on Mac OS X
  char ch;
  int fd = fileno(stdin);
  struct termio old_tty, new_tty;
  
  ioctl(fd, TCGETA, &old_tty);
  new_tty = old_tty;
  new_tty.c_lflag &= ~(ICANON | ECHO | ISIG);
  ioctl(fd, TCSETA, &new_tty);
  fread(&ch, 1, sizeof(ch), stdin);
  ioctl(fd, TCSETA, &old_tty);
  */
  // The following is the termios.h equivalent, which works both on linux and Max OS X
  struct termios oldattr, newattr;
  int ch;
  tcgetattr(STDIN_FILENO, &oldattr);
  newattr = oldattr;
  newattr.c_lflag &= ~( ICANON | ECHO | ISIG);    // ECHO means that the character appears on the terminal. ISIG means we deal with control-c and control-Z and stuff.
  tcsetattr(STDIN_FILENO, TCSANOW, &newattr);

  // If just want to wait for new character to be available, simply do:
  //  ch = getchar();

  // More complicated: check if character is available. If so, get it
  // from buffer. Can use this to wait for key from either terminal or
  // pgplot window for example.
  int gotkey = 0;
  int n;
  while(gotkey == 0) {
    // Check if key is available from terminal
    //    if(ioctl(fileno(stdin), I_NREAD, &n) == 0 && n > 0) {  // Internet said I_NREAD, I found FIONREAD in man pages.
    if(ioctl(fileno(stdin), FIONREAD, &n) == 0 && n > 0) {
      ch = getchar();
      gotkey = 1;
    }
    // Check if key is available from pgplot
    // ... [to be implemented]
    if(gotkey == 0) {
      usleep(10000);
    }
    //    printf("n = %d\n", n);
  }

  tcsetattr(STDIN_FILENO, TCSANOW, &oldattr);

  if(ch == 3) {
    fflush(stdout);
    fprintf(stderr, "pgetch: caught control-c\n");
    fprintf(stderr, "Terminating program\n");
    exit(0);
  }
  if(ch == 26) {
    fflush(stdout);
    fprintf(stderr, "pgetch: caught control-z\n");
    fprintf(stderr, "Sending suspend signal. Program is not terminated.\n");
    raise(SIGTSTP);
  }
  //#endif
  //  printf("XXXXX %d\n", ch);
  return ch;
}


/* Like function pgetch(), except that it uses the on the command line
   defined macro file if available instead of the keyboard. A ^ is
   interpreted as a ctrl-key. Line feeds/returns are ignored in
   macro's. If end of macro is reached, it will be closed and input
   will happen from keyboard. */
int pgetch_macro(psrsalsaApplication *application, verbose_definition verbose)
{
  int key;
  int ctrlstate;

  ctrlstate = 0;
  do {
    if(application->macro_ptr == NULL) { // If no macro defined, use the keyboard instead
      return pgetch();
    }else {
      key = fgetc(application->macro_ptr);
      if(key == EOF) {
	printf("\nReached end of macro, switching to keyboard input\n");
	fclose(application->macro_ptr);
	application->macro_ptr = NULL;
	ctrlstate = 0;
      }
    }
    if(key == '^') {
      if(ctrlstate == 0) {
	ctrlstate = 1;
      }else {
	fflush(stdout);
	printerror(verbose.debug, "ERROR pgetch_macro: key sequence ^^ is not allowed");
	exit(0);
      }
    }
  }while(key == '\n' || key == '\r' || key == '^' || key == EOF);

  if(ctrlstate) {
    if(key == 'c' || key == 'C') {
      fflush(stdout);
      fprintf(stderr, "pgetch_macro: caught control-c\n");
      exit(0);
    }
    if(key >= 'a' && key <= 'z') {   // ctrl-a to ctrl-z
      return key - 96;
    }else if(key >= 'A' && key <= 'Z') {  // ctrl-A to ctrl-Z
      return key - 64;
    }else if(key == 32) {   // ctrl-space
      return 0;
    }
  }
  
  return key;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Tries to find out the username of the person running the
   program. Memory will be allocated. Returns 1 on success, 0
   on error. */
int getUsername(char **username, verbose_definition verbose)
{
  /* Non-recurrent version with memory leak:
  struct passwd *pw;
  pw = getpwuid(geteuid());
  */
  // New code
  //  int ret;
  struct passwd pw, *result;
  char *tmpmem;
  tmpmem = malloc(MaxStringLength);  // This is a guess for a reasonable buffer size. If too small, an error will happen.
  if(tmpmem == NULL) {
    printwarning(verbose.debug, "ERROR getUsername: Memory allocation error.");
    *username = malloc(8);
    if(*username != NULL) {
      sprintf(*username, "Unknown");
    }
    return 0;
  }
  //  ret =
  getpwuid_r(geteuid(), &pw, tmpmem, MaxStringLength, &result);
  //  if(pw) {
  if(result != NULL) {
    //    user = pw->pw_name;
    *username = malloc(strlen(pw.pw_name)+1);
    if(*username == NULL) {
      printwarning(verbose.debug, "ERROR getUsername: Memory allocation error.");
      return 0;
    }
    strcpy(*username, pw.pw_name);
    if(verbose.debug) {
      printf("Username: %s\n", *username);
    }
    free(tmpmem);
    return 1;
  }
  free(tmpmem);

  // Try environment variable instead
  char *user;
  if((user = getenv ("USER")) == NULL) {
    /* Could consider to try to use the change_cmd_line function below instead */
    fflush(stdout);
    printerror(verbose.debug, "ERROR: getUsername failed to retrieve user name from USER environment variable.");
    *username = malloc(8);
    if(*username != NULL) {
      sprintf(*username, "Unknown");
    }
    return 0;
  }else {
    *username = malloc(strlen(user)+1);
    if(*username == NULL) {
      printwarning(verbose.debug, "ERROR getUsername: Memory allocation error.");
      return 0;
    }
    strcpy(*username, user);
    if(verbose.debug) {
      printf("Username: %s\n", *username);
    }
  }

  return 1;
}


/*
  Put command line in string txt with maximum length length
 */
void constructCommandLineString(char *txt, int length, int argc, char **argv, verbose_definition verbose)
{
  int i;
  txt[0] = 0;
  for(i = 0; i < argc; i++) {
    if(strlen(txt) + strlen(argv[i]) + 4 > length-1) {
      printwarning(verbose.debug, "WARNING constructCommandLineString: Truncating command line which is too long.");
      break;
    }
    if(strchr(argv[i], ' ') == NULL) {
      strcat(txt, argv[i]);
    }else {
      strcat(txt, "\"");
      strcat(txt, argv[i]);
      strcat(txt, "\"");
    }
    if(i != argc-1)
      strcat(txt, " ");
  }
}

//START REGION DEVELOP

/* Put the user and command line in the header. */
/*
void change_cmd_line(Header_type *hdr, int argc, char **argv)
{
  int i;
  FILE *tfile;
  char commline_s[256];
  sprintf(commline_s,"rm -f /tmp/who");
  system(commline_s);
  sprintf(commline_s,"whoami > /tmp/who");
  system(commline_s);
  
  tfile = fopen("/tmp/who","r");
  fscanf(tfile,"%[^\n]s",commline_s);
  fclose(tfile);

  strcat(commline_s," :");
  i = 0;
  while (i<argc)
    {
      strcat(commline_s,argv[i]);
      strcat(commline_s," ");
      i++;
    }
  strncpy(hdr->redn.Command,commline_s,256);
}
*/

//START REGION DEVELOP
//START REGION RELEASE

/* Set size to the size of the array. Returns 1 on success, 0 on
   error. */
int getMachinename(char *hostname, int size, verbose_definition verbose)
{
  hostname[size-1] = '\0';
  if(gethostname(hostname, size-1) == 0) {
    if(verbose.debug) {
      printf("Hostname: %s\n", hostname);
    }
    /*Obsolete, and it doesn't seem to do much on systems I know, so removed this lookup. It was causing a memory leak as well, since hostent should be freed on some systems, but it could be static as well......
Tried to replaced  the following code with the code below with the _r version to avoid memory leak, but that didn't resulted in solving the memory leak.
    struct hostent* h;
    h = gethostbyname(hostname);   // I guess this was to convert an IP address to a machine name, if required. Obsolete, use getaddrinfo() and getnameinfo() instead.
    if(h != NULL) {
      if(verbose.debug) {
	printf("h_name: %s\n", h->h_name); 
      }
      strncpy(hostname, h->h_name, size-1);
    }
    */

    /* Here is the _r version of the code, which still has the same memory leak....
    char *tmpmem;
    tmpmem = malloc(MaxStringLength);  // This is a guess for a reasonable buffer size. If too small, an error will happen.
    if(tmpmem == NULL) {
      printwarning(verbose.debug, "WARNING getMachinename: Memory allocation error.");
      sprintf(hostname, "Unknown");
      return 0;
    }
    struct hostent h, *result;
    int ret, h_errnop;
    ret = gethostbyname_r(hostname, &h, tmpmem, MaxStringLength, &result, &h_errnop);
    if(result == NULL) {
      printwarning(verbose.debug, "WARNING getMachinename: Call to gethostbyname_r() failed.");
      if(ret == ERANGE) {
	printwarning(verbose.debug, "WARNING getMachinename: Call to gethostbyname_r() failed because of a too small default temporary memory allocation.");
      }
      sprintf(hostname, "Unknown");
      free(tmpmem);
      return 0;
    }
    if(verbose.debug) {
      printf("h_name: %s\n", h.h_name); 
    }
    strncpy(hostname, h.h_name, size-1);
    free(tmpmem);
    */
    return 1;
  }else {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING getMachinename failed to find the machine name.");
    sprintf(hostname, "Unknown");
    return 0;
  }
}

//START REGION DEVELOP

/* Try to find out what the IP address is of hostname. ip_host is
   assumed to be (at least) 128 bytes. Returns 0 if lookup failed. The
   verbose output is sent to verbose_fptr, which can for instance be
   stdio. */
int getIPofHost(char *hostname, char *ip_host, FILE *verbose_fptr, verbose_definition verbose)
{
  struct hostent *host;     
  int j;

  ip_host[0] = 0;
  // Memory leak. I suspect this is unavoidable
  if ((host = gethostbyname(hostname)) == NULL) {     
    fflush(stdout);
    printwarning(verbose.debug, "getIPofHost:  (mini) nslookup failed on '%s'", hostname);
    return 0;
  }else {
    if(verbose.verbose) {
      fprintf(verbose_fptr, "Find out information about host %s\n",host->h_name); 
      fprintf(verbose_fptr, "  aliases:\n"); 
      for (j=0; host->h_aliases[j]; j++) { 
	fprintf(verbose_fptr, "    %s\n",host->h_aliases[j]); 
      } 
      fprintf(verbose_fptr, "  addresses:\n");
    }
    for (j=0; host->h_addr_list[j]; j++) { 
      char buffer[128]; 
      inet_ntop(PF_INET, host->h_addr_list[j],buffer,128); 
      if(verbose.verbose)
	fprintf(verbose_fptr, "    %s\n",buffer); 
      if(j == 0) {
	strcpy(ip_host, buffer);
      }
    }
  }
  if(verbose.verbose)
    fprintf(verbose_fptr, "  assume IP address of %s is %s\n", hostname, ip_host);
  return 1;
}

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
char * pickWordFromString(char *string, int n, int *nrwords, int replacetabs, char separator, verbose_definition verbose)
{
  char *ptr, *ret, *field, *string_mod, *format;
  int nrchars;

  if(n <= 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pickWordFromString: requested word number (%d) is not valid", n);
    if(verbose.debug) {
      printerror(verbose.debug, "ERROR pickWordFromString: string=%s", string);
    }
    exit(0);
  }

  //  if(separator != ' ' && separator != ',' && separator != ':') {
  //    fflush(stdout);
  //    fprintf(stderr, "ERROR pickWordFromString: This particular separator (ascii code %d) is not implemented", separator);
  //    exit(0);
  //  }

  // Allocate memory for string + 0 termination byte
  string_mod = malloc(strlen(string)+1);
  // Allocate memory to hold a word from the string
  field = malloc(MaxPickWordFromString_WordLength+1);
  format = malloc(20);
  if(string_mod == NULL || field == NULL || format == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR pickWordFromString: Memory allocation error");
    exit(0);
  }

  // Make a copy of the input string, so we can modify it
  strcpy(string_mod, string);
  //  printf("pickWordFromString: Got string: '%s'\n", string_mod);

  // Replace tabs with spaces
  if(replacetabs) {
    for(nrchars = 0; nrchars < strlen(string_mod); nrchars++) {
      if(string_mod[nrchars] == '\t')
	string_mod[nrchars] = ' ';
    }
  }
  //  printf("pickWordFromString: After tab replacement: '%s'\n", string_mod);

  // Ignore any trailing spaces and strip return
  for(nrchars = strlen(string_mod)-1; nrchars >= 0; nrchars--) {
    if(string_mod[nrchars] == '\n')
      string_mod[nrchars] = 0;
    else if(string_mod[nrchars] == '\r')
      string_mod[nrchars] = 0;
    else if(string_mod[nrchars] == ' ')
      string_mod[nrchars] = 0;
    else
      break;
  }
  //  printf("pickWordFromString: After removing return/trailing spaces: '%s'\n", string_mod);

  ptr = string_mod;
  *nrwords = 0;
  ret = NULL;
  do {
    while(*ptr == ' ')   /* Skip any extra spaces */
      ++ptr;
    if(*ptr == 0)   /* Reached end of string */
      break;
    nrchars = 0;
    // The scanf format "%1000[^ ]%n" means: convert at least 1000 characters in the output string, and the conversion accepts all characters except a space (indicated after the ^). After conversion stops, return the number of converted characters as well.
    sprintf(format, "%%%d[^%c]%%n", MaxPickWordFromString_WordLength, separator);
    //    printf("XXXXXX format=%s\n", format);
    sscanf(ptr, format, field, &nrchars);
    //    if(separator == ' ')
    //      sscanf(ptr, "%1000[^ ]%n", field, &nrchars);
    //    else if(separator == ',')
    //      sscanf(ptr, "%1000[^,]%n", field, &nrchars);
    //    else if(separator == ':')
    //      sscanf(ptr, "%1000[^:]%n", field, &nrchars);
    (*nrwords) ++;
    if(*nrwords == n)
      ret = ptr;
    /*    printf("nrchars=%d ('%s')\n", nrchars, field);
	  printf("field = \"%s\"\n", field); */
    ptr += nrchars; /* advance the pointer by the number of characters read */
    if ( *ptr != separator ) {
      /*      printf("Quit because of %d (%c)\n", *ptr, *ptr); */
      break; /* didn't find an expected delimiter, done? */
    }
    ++ptr; /* skip the delimiter */
  }while(1);
  //  printf("pickWordFromString: Identified '%s' as start word %d (out of %d)\n", ret, n, *nrwords);

  // Need to return pointer w.r.t. provided string rather than modified copy we worked on above
  if(ret != NULL)
    ret = ret-string_mod+string;  // Note this works because only made replacements and only cut at end
  free(string_mod);
  free(field);
  free(format);
  return ret;
}

//START REGION DEVELOP

/* This function finds the first occurrence of the substring needle in
   the string haystack. This function is similar to C function strstr,
   except that it works on a block of memory with a given size, rather
   than on an ascii null terminated string. It returns a pointer to
   the beginning of the substring, or NULL if the substring is not
   found.
 */
char *searchStringInMem(const char *haystack, int haystacksize, const char *needle, int needlesize, verbose_definition verbose)
{
  long index, found;
  char *hdr_ptr;

  if(needlesize > haystacksize) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR: searchStringInMem called with a needle larger than the haystack!");
    return NULL;
  }

  found = 0;
  hdr_ptr = (char *)haystack;
  for(index = 0; index <= haystacksize-needlesize; index++) {
    if(memcmp(needle, hdr_ptr, needlesize) == 0) {
      found = 1;
      break;
    }
    hdr_ptr++;
  }
  if(found)
    return hdr_ptr;
  return NULL;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Takes inputname, cut extension and add new extension. inputname and
   extension are not changed. A check is done if length of outputname
   will not exceed outputnamelength. Returns 1 if successful or 0 on
   error. */
int change_filename_extension(char *inputname, char *outputname, char *extension, int outputnamelength, verbose_definition verbose)
{
  int i;
  for(i = strlen(inputname); i >= 0; i--) {
    if(inputname[i] == '.') {
      break;
    }
  }
  if(i <= 0) {
    fflush(stdout);
    printerror(verbose.debug, "change_filename_extension: no extension in '%s'?", inputname);
    return 0;
  }
  if(strlen(extension) + i+1 >= outputnamelength) {
    fflush(stdout);
    printerror(verbose.debug, "change_filename_extension: outputnamelength too long");
    return 0;
  }
  memcpy(outputname, inputname, i+1);
  outputname[i+1] = 0;
  strcat(outputname, extension);
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* This function skips lines in a file. It returns 1 on success, 0
   when reaching EOF */
int skipLinesInFile(FILE *fptr, int skiplines, verbose_definition verbose)
{
  long i, j;
  if(skiplines > 0) {
    if(verbose.verbose)
      printf("Skipping %d lines\n", skiplines);
    for(i = 0; i < skiplines; i++) {
      do {
	j = fgetc(fptr);
	if(j == EOF) {
	  return 0;
	}
      }while(j !='\n');
    }
  }
  return 1;
}

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
int ascii_file_get_next_line(FILE *fin, char *txt, int maxlinelength, int skipChar, verbose_definition verbose)
{
  char *ret_ptr;
  int ret;
  ret = 0;
  do {
    ret_ptr = fgets(txt, maxlinelength, fin);
    if(ret_ptr == NULL) {
      return 0;
    }else {
      ret++;
      if(txt[0] != skipChar) {
	return ret;
      }
    }
  }while(txt[0] == skipChar);
  printerror(verbose.debug, "ascii_file_get_next_line: BUG!!!");
  exit(0);
}

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
int ascii_file_stats(FILE *fin, char skipChar, long *nrlines, int maxlinelength, int autoNrColumns, int *nrColumns, verbose_definition verbose)
{
  char *txt;
  int nrwords, ret;
  
  txt = malloc(maxlinelength);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ascii_file_stats: Cannot allocate temporary memory");
    return 0;    
  }
  if(verbose.debug) {
    printf("ascii_file_stats: skipChar=%c\n", skipChar);
  }

  *nrlines = 0;
  do {
    ret = ascii_file_get_next_line(fin, txt, maxlinelength, skipChar, verbose);
    if(verbose.debug) {
      printf("ascii_file_stats: line %ld=%s\n", *nrlines, txt);
    }
    if(ret != 0) {
      pickWordFromString(txt, 1, &nrwords, 1, ' ', verbose);
      //      printf("XXXXXX line %d: %d words\n", *nrlines, nrwords);
      if(autoNrColumns != 0 && *nrlines == 0) {
	*nrColumns = nrwords;
      }else {
	if(nrColumns != NULL) {
	  if(*nrColumns >= 0) {
	    if(*nrColumns != nrwords) {
	      fflush(stdout);
	      printerror(verbose.debug, "ascii_file_stats: Nr of columns on line %ld is not the expected %d. Possibly the number of columns is changing from line to line?", (*nrlines)+1, *nrColumns);
	      free(txt);
	      return 0;
	    }
	  }else {
	    if(nrwords < -(*nrColumns)) {
	      fflush(stdout);
	      printerror(verbose.debug, "ascii_file_stats: Nr of columns on line %ld is smaller than the expected %d. Possibly the number of columns is changing from line to line?", (*nrlines)+1, *nrColumns);
	      free(txt);
	      return 0;
	    }
	  }
	}
      }
      (*nrlines)++;
    }
  }while(ret != 0);

  free(txt);
  return 1;
}

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
  level detirmines nr of spaces before output. Lines cannot exceed 10k
  lenght. The individual words cannot be larger than 1000 bytes, or
  else there will be a memory overflow.

  Returns 0 on error
*/
// THERE IS ALSO A DOUBLE AND INT AND STR VERSION BELOW!!!!!!
int read_ascii_column(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, int colnum, double scale, int read_log, float **data, float *mindata, float *maxdata, float *avdata, verbose_definition verbose, int verbose_stderr)
{
  FILE *fin, *verbose_stream;
  long i, j, n, fpos, maxlinelength, linenr;
  int nrwords, ret;
  char *txt, *word_ptr;
  double minx, maxx, sumx;

  maxlinelength = 10240+1;
  

  if(verbose_stderr) {
    fflush(stdout);
    verbose_stream = stderr;
  }else {
    verbose_stream = stdout;
  }
  fin = fopen(fname, "r");
  if(fin == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column: Cannot open %s", fname);
    return 0;
  }else {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      fflush(stdout);
      fprintf(verbose_stream, "Opened file '%s'\n", fname);
    }
  }

  // Skip first lines
  if(skipLinesInFile(fin, skiplines, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column: Reached EOF while skipping first %d lines", skiplines);
    return 0;
  }
  fpos = ftell(fin);


  txt = malloc(maxlinelength);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column: Cannot allocate temporary memory");
    return 0;    
  }

  // Find out how many lines there are in the file, and determine the nr of columns
  if(ascii_file_stats(fin, skipChar, nrdatapoints, maxlinelength, autoNrColumns, &nrColumns, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column: Error in determining the nr of lines");
    free(txt);
    return 0;
  }

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  There are %ld datapoints\n", *nrdatapoints);
    if(autoNrColumns) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      fprintf(verbose_stream, "  There are %d columns\n", nrColumns);
    }
  }

  *data = (float *)malloc((*nrdatapoints)*sizeof(float));
  if(*data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column: Memory allocation error");
    free(txt);
    return 0;
  }

  fseek(fin, fpos, SEEK_SET);
  minx = maxx = NAN;
  sumx = 0;
  n = 0;
  linenr = 0;
  do {
    ret = ascii_file_get_next_line(fin, txt, maxlinelength, skipChar, verbose);
    if(ret != 0) {
      linenr += ret;
      word_ptr = pickWordFromString(txt, colnum, &nrwords, 1, ' ', verbose);
      if(n >= *nrdatapoints) {
	fflush(stdout);
	printerror(verbose.debug, "read_ascii_column: Nr of lines in file changed?????");
	free(txt);
	free(*data);
	return 0;
      }
      if(word_ptr == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "read_ascii_column: Cannot find column %d on line %ld", colnum, linenr);
	free(txt);
	free(*data);
	return 0;
      }
      j = sscanf(word_ptr, "%f", &((*data)[n]));
      if(j != 1) {
	fflush(stdout);
	printerror(verbose.debug, "read_ascii_column: Cannot interpret column %d on line %ld as a float", colnum, linenr);
	free(txt);
	free(*data);
	return 0;
      }

      (*data)[n] *= scale;
      if(read_log) {
	if((*data)[n] <= 0) {
	  printerror(verbose.debug, "read_ascii_column: Cannot take logarithm of a value <= 0");
	  return 0;
	}
	(*data)[n] = log10((*data)[n]);
      }
      if((*data)[n] < minx || n == 0) {
	minx = (*data)[n];
      }
      if((*data)[n] > maxx || n == 0) {
	maxx = (*data)[n];
      }
      sumx += (*data)[n];
      
      n++;
    }
  }while(ret != 0);

  free(txt);
  fclose(fin);

  if(verbose.verbose) {
    fflush(stdout);
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  %ld points loaded from %s with values between %lf and %lf\n", *nrdatapoints, fname, minx, maxx);
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    fprintf(verbose_stream, "  Average value = %lf\n", sumx/(double)(*nrdatapoints));
  }
  if(mindata != NULL)
    *mindata = minx;
  if(maxdata != NULL)
    *maxdata = maxx;
  if(avdata != NULL)
    *avdata = sumx/(double)(*nrdatapoints);

  return 1;
}

// SAME AS ABOVE, WITH FLOAT -> DOUBLE AND IN SCANF %f -> %lf and _column _column_double
int read_ascii_column_double(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, int colnum, double scale, int read_log, double **data, double *mindata, double *maxdata, double *avdata, verbose_definition verbose, int verbose_stderr)
{
  FILE *fin, *verbose_stream;
  long i, j, n, fpos, maxlinelength, linenr;
  int nrwords, ret;
  char *txt, *word_ptr;
  double minx, maxx, sumx;

  maxlinelength = 10240+1;
  

  if(verbose_stderr) {
    fflush(stdout);
    verbose_stream = stderr;
  }else {
    verbose_stream = stdout;
  }
  fin = fopen(fname, "r");
  if(fin == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_double: Cannot open %s", fname);
    return 0;
  }else {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      fflush(stdout);
      fprintf(verbose_stream, "Opened file '%s'\n", fname);
    }
  }

  // Skip first lines
  if(skipLinesInFile(fin, skiplines, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_double: Reached EOF while skipping first %d lines", skiplines);
    return 0;
  }
  fpos = ftell(fin);


  txt = malloc(maxlinelength);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_double: Cannot allocate temporary memory");
    return 0;    
  }

  // Find out how many lines there are in the file, and determine the nr of columns
  if(ascii_file_stats(fin, skipChar, nrdatapoints, maxlinelength, autoNrColumns, &nrColumns, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_double: Error in determining the nr of lines");
    free(txt);
    return 0;
  }

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  There are %ld datapoints\n", *nrdatapoints);
    if(autoNrColumns) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      fprintf(verbose_stream, "  There are %d columns\n", nrColumns);
    }
  }

  *data = (double *)malloc((*nrdatapoints)*sizeof(double));
  if(*data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_double: Memory allocation error");
    free(txt);
    return 0;
  }

  fseek(fin, fpos, SEEK_SET);
  minx = maxx = NAN;
  sumx = 0;
  n = 0;
  linenr = 0;
  do {
    ret = ascii_file_get_next_line(fin, txt, maxlinelength, skipChar, verbose);
    if(ret != 0) {
      linenr += ret;
      word_ptr = pickWordFromString(txt, colnum, &nrwords, 1, ' ', verbose);
      if(n >= *nrdatapoints) {
	fflush(stdout);
	printerror(verbose.debug, "read_ascii_column_double: Nr of lines in file changed?????");
	free(txt);
	free(*data);
	return 0;
      }
      if(word_ptr == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "read_ascii_column_double: Cannot find column %d on line %ld", colnum, linenr);
	free(txt);
	free(*data);
	return 0;
      }
      j = sscanf(word_ptr, "%lf", &((*data)[n]));
      if(j != 1) {
	fflush(stdout);
	printerror(verbose.debug, "read_ascii_column_double: Cannot interpret column %d on line %ld as a double", colnum, linenr);
	free(txt);
	free(*data);
	return 0;
      }

      (*data)[n] *= scale;
      if(read_log) {
	if((*data)[n] <= 0) {
	  printerror(verbose.debug, "read_ascii_column: Cannot take logarithm of a value <= 0");
	  return 0;
	}
	(*data)[n] = log10((*data)[n]);
      }
      if((*data)[n] < minx || n == 0) {
	minx = (*data)[n];
      }
      if((*data)[n] > maxx || n == 0) {
	maxx = (*data)[n];
      }
      sumx += (*data)[n];
      
      n++;
    }
  }while(ret != 0);

  free(txt);
  fclose(fin);

  if(verbose.verbose) {
    fflush(stdout);
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  %ld points loaded from %s with values between %lf and %lf\n", *nrdatapoints, fname, minx, maxx);
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    fprintf(verbose_stream, "  Average value = %lf\n", sumx/(double)(*nrdatapoints));
  }
  if(mindata != NULL)
    *mindata = minx;
  if(maxdata != NULL)
    *maxdata = maxx;
  if(avdata != NULL)
    *avdata = sumx/(double)(*nrdatapoints);

  return 1;
}

// SAME AS ABOVE, WITH DOUBLE -> int (except for avdata) AND IN SCANF %lf -> %d and _column_double _column_int
// Remove the scale option and read_log
int read_ascii_column_int(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, int colnum, int **data, int *mindata, int *maxdata, double *avdata, verbose_definition verbose, int verbose_stderr)
{
  FILE *fin, *verbose_stream;
  long i, j, n, fpos, maxlinelength, linenr;
  int nrwords, ret;
  char *txt, *word_ptr;
  long minx, maxx, sumx;

  maxlinelength = 10240+1;
  

  if(verbose_stderr) {
    fflush(stdout);
    verbose_stream = stderr;
  }else {
    verbose_stream = stdout;
  }
  fin = fopen(fname, "r");
  if(fin == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_int: Cannot open %s", fname);
    return 0;
  }else {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      fflush(stdout);
      fprintf(verbose_stream, "Opened file '%s'\n", fname);
    }
  }

  // Skip first lines
  if(skipLinesInFile(fin, skiplines, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_int: Reached EOF while skipping first %d lines", skiplines);
    return 0;
  }
  fpos = ftell(fin);


  txt = malloc(maxlinelength);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_int: Cannot allocate temporary memory");
    return 0;    
  }

  // Find out how many lines there are in the file, and determine the nr of columns
  if(ascii_file_stats(fin, skipChar, nrdatapoints, maxlinelength, autoNrColumns, &nrColumns, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_int: Error in determining the nr of lines");
    free(txt);
    return 0;
  }

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  There are %ld datapoints\n", *nrdatapoints);
    if(autoNrColumns) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      fprintf(verbose_stream, "  There are %d columns\n", nrColumns);
    }
  }

  *data = (int *)malloc((*nrdatapoints)*sizeof(int));
  if(*data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_int: Memory allocation error");
    free(txt);
    return 0;
  }

  fseek(fin, fpos, SEEK_SET);
  minx = maxx = 0;
  sumx = 0;
  n = 0;
  linenr = 0;
  do {
    ret = ascii_file_get_next_line(fin, txt, maxlinelength, skipChar, verbose);
    if(ret != 0) {
      linenr += ret;
      word_ptr = pickWordFromString(txt, colnum, &nrwords, 1, ' ', verbose);
      if(n >= *nrdatapoints) {
	fflush(stdout);
	printerror(verbose.debug, "read_ascii_column_int: Nr of lines in file changed?????");
	free(txt);
	free(*data);
	return 0;
      }
      if(word_ptr == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "read_ascii_column_int: Cannot find column %d on line %ld", colnum, linenr);
	free(txt);
	free(*data);
	return 0;
      }
      j = sscanf(word_ptr, "%d", &((*data)[n]));
      if(j != 1) {
	fflush(stdout);
	printerror(verbose.debug, "read_ascii_column_int: Cannot interpret column %d on line %ld as a int", colnum, linenr);
	free(txt);
	free(*data);
	return 0;
      }

      if((*data)[n] < minx || n == 0) {
	minx = (*data)[n];
      }
      if((*data)[n] > maxx || n == 0) {
	maxx = (*data)[n];
      }
      sumx += (*data)[n];
      
      n++;
    }
  }while(ret != 0);

  free(txt);
  fclose(fin);

  if(verbose.verbose) {
    fflush(stdout);
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  %ld points loaded from %s with values between %ld and %ld\n", *nrdatapoints, fname, minx, maxx);
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    fprintf(verbose_stream, "  Average value = %lf\n", sumx/(double)(*nrdatapoints));
  }
  if(mindata != NULL)
    *mindata = minx;
  if(maxdata != NULL)
    *maxdata = maxx;
  if(avdata != NULL)
    *avdata = sumx/(double)(*nrdatapoints);

  return 1;
}

// SAME AS ABOVE DOUBLE VERSION, WITH VARIOUS CHANGES TO MAKE IT WORK WITH STRINGS
// Remove the scale option and read_log and min/maxdata and avdata
// char **data-> char ***data
int read_ascii_column_str(char *fname, int skiplines, char skipChar, int nrColumns, int autoNrColumns, long *nrdatapoints, int colnum, char ***data, verbose_definition verbose, int verbose_stderr)
{
  FILE *fin, *verbose_stream;
  long i, n, fpos, maxlinelength, linenr;
  int nrwords, ret;
  char *txt, *word, *word_ptr;

  maxlinelength = 10240+1;
  

  if(verbose_stderr) {
    fflush(stdout);
    verbose_stream = stderr;
  }else {
    verbose_stream = stdout;
  }
  fin = fopen(fname, "r");
  if(fin == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_str: Cannot open %s", fname);
    return 0;
  }else {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      fflush(stdout);
      fprintf(verbose_stream, "Opened file '%s'\n", fname);
    }
  }

  // Skip first lines
  if(skipLinesInFile(fin, skiplines, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_str: Reached EOF while skipping first %d lines", skiplines);
    return 0;
  }
  fpos = ftell(fin);


  txt = malloc(maxlinelength);
  word = malloc(maxlinelength);
  if(txt == NULL || word == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_str: Cannot allocate temporary memory");
    return 0;    
  }

  // Find out how many lines there are in the file, and determine the nr of columns
  if(ascii_file_stats(fin, skipChar, nrdatapoints, maxlinelength, autoNrColumns, &nrColumns, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_str: Error in determining the nr of lines");
    free(txt);
    free(word);
    return 0;
  }

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  There are %ld datapoints\n", *nrdatapoints);
    if(autoNrColumns) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      fprintf(verbose_stream, "  There are %d columns\n", nrColumns);
    }
  }

  *data = (char **)malloc((*nrdatapoints)*sizeof(char *));
  if(*data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "read_ascii_column_str: Memory allocation error");
    free(word);
    free(txt);
    return 0;
  }

  fseek(fin, fpos, SEEK_SET);
  n = 0;
  linenr = 0;
  do {
    ret = ascii_file_get_next_line(fin, txt, maxlinelength, skipChar, verbose);
    if(ret != 0) {
      linenr += ret;
      word_ptr = pickWordFromString(txt, colnum, &nrwords, 1, ' ', verbose);
      if(n >= *nrdatapoints) {
	fflush(stdout);
	printerror(verbose.debug, "read_ascii_column_str: Nr of lines in file changed?????");
	free(word);
	free(txt);
	free(*data);
	return 0;
      }
      if(word_ptr == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "read_ascii_column_str: Cannot find column %d on line %ld", colnum, linenr);
	free(word);
	free(txt);
	free(*data);
	return 0;
      }
      sscanf(word_ptr, "%s", word);
      (*data)[n] = (char *)malloc(strlen(word)+1);
      if((*data)[n] == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "read_ascii_column_str: Memory allocation error");
	free(word);
	free(txt);
	return 0;
      }
      strcpy((*data)[n], word);
      n++;
    }
  }while(ret != 0);

  free(word);
  free(txt);
  fclose(fin);

  if(verbose.verbose) {
    fflush(stdout);
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    fflush(stdout);
    fprintf(verbose_stream, "  %ld words loaded from %s\n", *nrdatapoints, fname);
  }
  return 1;
}



//START REGION DEVELOP
//START REGION RELEASE

// Given string orig, replace every occurance of rep with with. The
// return value is the result (memory will be allocated, so you must free
// the result if result is non-NULL). The function returns NULL if memory
// cannot be allocated.
char *str_replace(char *orig, char *rep, char *with, verbose_definition verbose) 
{
  char *result; // the return string
  char *ins;    // the next insert point
  char *tmp;    // varies
  int len_rep;  // length of rep
  int len_with; // length of with
  int len_front; // distance between rep and end of last rep
  int count;    // number of replacements
  int new_string_length;

  if(orig == NULL)
    return NULL;
  if(rep == NULL)
    rep = "";
  len_rep = strlen(rep);
  if(with == NULL)
    with = "";
  len_with = strlen(with);

  ins = orig;
  for(count = 0; (tmp = strstr(ins, rep)) != NULL; ++count) {
    ins = tmp + len_rep;
  }

  // first time through the loop, all the variable are set correctly
  // from here on,
  //    tmp points to the end of the result string
  //    ins points to the next occurrence of rep in orig
  //    orig points to the remainder of orig after "end of rep"
  new_string_length = strlen(orig) + (len_with - len_rep) * count + 1;
  //  printf("XXXXX %d\n", new_string_length);
  result = malloc(new_string_length);
  tmp = result;

  if(result == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR str_replace: Cannot allocate memory");
    return NULL;    
  }

  while(count--) {
    ins = strstr(orig, rep);
    len_front = ins - orig;
    tmp = strncpy(tmp, orig, len_front) + len_front;
    tmp = strcpy(tmp, with) + len_with;
    orig += len_front + len_rep; // move to next "end of rep"
  }
  strcpy(tmp, orig);
  return result;
}

//START REGION DEVELOP
//START REGION RELEASE

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
void fprintf_color(FILE *destination, int color, const char *format, ...)
{
  va_list args;
  // If destination is a terminal (not redirected to a file), sent escape code to change colour
  if(isatty(fileno(destination))) {
    switch(color) {
    case 2: fprintf(destination, "\x1B[31m"); break;
    case 3: fprintf(destination, "\x1B[32m"); break;
    case 4: fprintf(destination, "\x1B[33m"); break;
    case 5: fprintf(destination, "\x1B[34m"); break;
    case 6: fprintf(destination, "\x1B[35m"); break;
    case 7: fprintf(destination, "\x1B[36m"); break;
    case 8: fprintf(destination, "\x1B[37m"); break;
    }
  }
  va_start(args, format);
  vfprintf(destination, format, args);
  // Reset colour to default
  if(isatty(fileno(destination))) {
    fprintf(destination, "\x1B[0m");
  }
}

//START REGION DEVELOP

void close_files_internal(int fd_to, int fd_from)
{
  int saved_errno;
  saved_errno = errno;

  close(fd_from);
  if(fd_to >= 0)
    close(fd_to);

  errno = saved_errno;
}


// Copy file from to destination. If dest already exist, it will be overwritten.
// Returns 0 on success
int cp(const char *from, const char *dest, verbose_definition verbose)
{
  int fd_dest, fd_from;
  char buf[4096];
  ssize_t nread;

  fd_dest = -1;
  fd_from = open(from, O_RDONLY);
  if (fd_from < 0) {
    fprintf(stderr, "ERROR: Opening file %s for reading failed.\n", from);
    close_files_internal(fd_dest, fd_from);
    return -1;
  }
  
  fd_dest = open(dest, O_WRONLY | O_CREAT | O_TRUNC, 0666);
  if (fd_dest < 0) {
    printerror(verbose.debug, "ERROR cp: Opening file %s for writing failed.\n", dest);
    close_files_internal(fd_dest, fd_from);
    return -1;
  }
  while (nread = read(fd_from, buf, sizeof buf), nread > 0)
    {
      char *out_ptr = buf;
      ssize_t nwritten;

      do {
	nwritten = write(fd_dest, out_ptr, nread);
	
	if (nwritten >= 0)
	  {
	    nread -= nwritten;
	    out_ptr += nwritten;
	  }
	else if (errno != EINTR)
	  {
	    printerror(verbose.debug, "ERROR cp: Writing data from file %s to %s failed.\n", from, dest);
	    close_files_internal(fd_dest, fd_from);
	    return -1;
	  }
      } while (nread > 0);
    }
  
    if (nread == 0)
    {
        if (close(fd_dest) < 0)
        {
	  printerror(verbose.debug, "ERROR cp: Writing data from file %s to %s failed.\n", from, dest);
	  close_files_internal(fd_dest, fd_from);
	  return -1;
        }
        close(fd_from);

        /* Success! */
        return 0;
    }

    printerror(verbose.debug, "ERROR cp: Writing data from file %s to %s failed.\n", from, dest);
    close_files_internal(fd_dest, fd_from);
    return -1;
}

//START REGION DEVELOP

// Read in byte values represented as hex numbers as a long double variable.
//
// Return 0: Error
// Return 1: Success
int convert_hexstring_to_longdouble(char *string, long double *value, verbose_definition verbose)
{
  unsigned char *byte;
  unsigned int intvalue;
  int bytenr;
  if(strlen(string) != 2*sizeof(long double)) {
    printerror(0, "Cannot parse %s as a long double (%d bytes)\n", string, sizeof(long double));
    return 0;
  }
  //  printf("Going to parse %s\n", string);
  for(bytenr = 0; bytenr < sizeof(long double); bytenr++) {
    sscanf(string+2*bytenr, "%02x", &intvalue);
    byte = (unsigned char *)value;
    byte += bytenr;
    *byte = intvalue;
  }
  return 1;
}

// Write a long double variable as byte values represented as hex numbers to stream
void print_longdouble_as_hexstring(FILE *stream, long double value)
{
  unsigned char *byte;
  int bytenr;
  for(bytenr = 0; bytenr < sizeof(long double); bytenr++) {
    byte = (unsigned char *)&value;
    byte += bytenr;
    fprintf(stream, "%02x", *byte);
  }
  fprintf(stream, " ");
}

//START REGION DEVELOP
