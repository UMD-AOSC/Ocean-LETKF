/*  runSWA computes the weekly dates of Altimeter data needed for a run
    and calls spltWkAlt to build each weekly data file
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Tm.h"

#define NWKS 7

void usage();

int nwks, nw2, idate;
int dow0, dow, yday, year0;
char prog[51];

char lstFile[81], grdFile[21], path[81], cmnd[181];
float Lts = -90.0, Ltn = 90.0;
int kmin = -1;

YMD sTm, wTm;

void main(argc,argv)
int argc;
char *argv[];
{
  FILE *fs;
  int n;
  char str[81];

  strcpy(prog,argv[0]);
  strcpy(lstFile,"EMPTY");
  strcpy(grdFile,"EMPTY");
  strcpy(path,"EMPTY");
  idate = -1;
  nwks = -1;
  n = 1;
  while (n < argc) {
    if (!strcmp(argv[n],"-d")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -d requires an integer as a parameter.\n");
        usage();
        exit(99);
      }
      sscanf(argv[n],"%d",&idate);
    } else if (!strcmp(argv[n],"-n")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -nrequires an integer as a parameter.\n");
        usage();
        exit(99);
      }
      sscanf(argv[n],"%d",&nwks);
    } else if (!strcmp(argv[n],"-f")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -f requires a file name as a parameter.\n");
        usage();
        exit(99);
      }
      strcpy(lstFile, argv[n]);
    } else if (!strcmp(argv[n],"-p")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -p requires a path as a parameter.\n");
        usage();
        exit(99);
      }
      strcpy(path, argv[n]);
    } else if (!strcmp(argv[n],"-m")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -m requires a file name as a parameter.\n");
        usage();
        exit(99);
      }
      strcpy(grdFile, argv[n]);
    } else if (!strcmp(argv[n],"-lt")) {
      n++;
      if (n >= argc || sscanf(argv[n],"%f",&Lts) != 1) {
        printf("Error: the option -lt requires 2 reals as parameters.\n");
        usage();
        exit(99);
      }
      n++;
      if (n >= argc || sscanf(argv[n],"%f",&Ltn) != 1) {
        printf("Error: the option -lt requires 2 reals as parameters.\n");
        usage();
        exit(99);
      }
    } else if (!strcmp(argv[n],"-k")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -k requires an integer as a parameter.\n");
        usage();
        exit(99);
      }
      sscanf(argv[n], "%d", &kmin);
    } else if (!strcmp(argv[n],"-h")) {
      usage();
      exit(99);
    } else {
      printf("Error: %s is not an option.\n", argv[n]);
      usage();
      exit(99);
    }
    n++;
  }
  if (idate < 0) {
    printf("Error: A start date for the model run must be given.\n");
    usage();
    exit(99);
  }
  if (nwks < 0) {
    printf("Error: The number of weeks of data to combine must be given.\n");
    usage();
    exit(99);
  }
  if (!strcmp(lstFile,"EMPTY")) {
    printf("Error: Input list file must be given.\n");
    usage();
    exit(99);
  }
  if (!strcmp(path,"EMPTY")) {
    printf("Error: A path must be given.\n");
    usage();
    exit(99);
  }
  if (!strcmp(grdFile,"EMPTY")) {
    printf("Error: A grid/mask file must be given.\n");
    usage();
    exit(99);
  }

  sTm.year = idate / 10000;
  sTm.month = 1;
  sTm.day = 1;
  year0 = sTm.year - 1;
  sTm.year0 = year0;
  dow0 = DayOfWeek(sTm);

  sTm.year = idate / 10000;
  sTm.month = (idate % 10000) / 100;
  sTm.day = idate % 100;
  year0 = sTm.year - 1;
  sTm.year0 = year0;
  sTm.yday = YearDay(sTm);

  dow = DayOfWeek(sTm);
  if (dow != 4) {
    yday = sTm.yday + 4 - dow;
    sTm.yday = yday;
    CalendarDay(&sTm);
  }
  sTm.year0 = year0;
  yday = YearDay(sTm);
  sTm.yday = yday;

  nw2 = nwks / 2;
  if (dow > 4) {
    wTm.yday = sTm.yday - 7*(nw2-1);
  } else {
    wTm.yday = sTm.yday - 7*nw2;
  }
  wTm.year0 = year0;

  fs = fopen("split.ksh","w");
  printf("Writing lines:");
  for (n = 0; n < nwks; n++) {
    CalendarDay(&wTm);
    idate = wTm.day + (wTm.month + wTm.year*100)*100;
    cmnd[0] = '\0';
    sprintf(cmnd,"./spltWkAlt4 -f %s -p %s -w %d -m %s -lt %g %g -k %d\n",
                            lstFile, path, idate, grdFile, Lts, Ltn, kmin);
//    fputs(cmnd,fs);
    printf("%s",cmnd);
    fprintf(fs,cmnd);
/*
    system(cmnd);
*/
    wTm.yday += 7;
  }
  fclose(fs);
}

/* ================================================================= */

void usage()
{
  printf("Usage:\n");
  printf(" %s -d date -n nWeeks ...\n", prog);
  printf("   -d date     - start date for model run (YYYYMMDD)\n");
  printf("   -n nWeeks   - number of weeks to combine\n");
  printf("  Options for call to spltWkAlt\n");
  printf("   -f lstFile  - file list of ascii altimetry cycles\n");
  printf("   -p path     - path for list and data files\n");
  printf("   -m grdFile  - \"t\" grid/mask file for MOM2 model\n");
  printf("   -lt LtS LtN - latitude boundaries\n");
  printf("   -k Kmin     - minimum depth index\n");

}
