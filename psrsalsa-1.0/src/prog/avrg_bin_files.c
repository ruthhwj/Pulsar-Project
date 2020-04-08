// Can be compiled as stand-alone with gcc -Wall -o avrg_bin_files avrg_bin_files.c

#include <stdio.h>
#include <psrsalsa.h>

int main(int argc, char **argv)
{
  long filesize, i;
  unsigned char byte1, byte2, result;
  FILE *fin1, *fin2, *fout;
  if(argc != 4) {
    printf("This program reads in two files which should be of the same size. A byte at the time is read in from both input files (assumed to be unsigned) and the average of the two is written to the output file. This can be used in combination with ImageMagick to obtain transparent figures. See the documentation of ppolFit.\n\n");
    printf("Usage: avrg_bin_files file1 file2 output_file\n\n");
    printCitationInfo();
    return 0;
  }
  fin1 = fopen(argv[1], "r");
  if(fin1 == NULL) {
    printf("Cannot open first file=%s\n", argv[1]);
    return 0;
  }
  fin2 = fopen(argv[2], "r");
  if(fin2 == NULL) {
    printf("Cannot open second file=%s\n", argv[2]);
    return 0;
  }
  fout = fopen(argv[3], "w");
  if(fout == NULL) {
    printf("Cannot open output file=%s\n", argv[3]);
    return 0;
  }

  fseek(fin1, 0, SEEK_END);
  filesize = ftell(fin1);
  fseek(fin2, 0, SEEK_END);
  if(filesize != ftell(fin2)) {
    printf("Input files have different sizes!\n");
    return 0;
  }
  rewind(fin1);
  rewind(fin2);


  for(i = 0; i < filesize; i++) {
    byte1 = fgetc(fin1);
    byte2 = fgetc(fin2);
    result = (byte1+byte2)/2;
    fputc(result, fout);
  }

  fclose(fout);
  fclose(fin1);
  fclose(fin2);

  return 0;
}
