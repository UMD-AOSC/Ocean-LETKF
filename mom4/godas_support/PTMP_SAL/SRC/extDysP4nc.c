/*  extDysP4nc extracts daily profile files        */

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Tm.h"

void usage();

int dt0, n1, n2, nDys, doT = 1, dbg = 0;
int year, month, day, year0, hour, minute;
YMD *tm;
char typ = 'O';

char path[121], prfFile[121];
char prog[51];

mode_t mode = 0644;
float spv = 999.999;

main(argc,argv)
int argc;
char *argv[];
{
  int i, k, n, nprf;
  char str[121];

  strcpy(prog,argv[0]);
  strcpy(path,"EMPTY");
  n = 1;
  while (n < argc) {
    if (!strcmp(argv[n],"-p")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -p requires a path as a parameter.\n");
        usage();
        exit(0);
      }
      strcpy(path, argv[n]);
    } else if (!strcmp(argv[n],"-d")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: for option -d a date must be given.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%d", &dt0);
    } else if (!strcmp(argv[n],"-b")) {
      n++;
      if (n >= argc || sscanf(argv[n],"%d", &n1) != 1) {
        printf("Error: for option -b 2 integers must be given.\n");
        usage();
        exit(0);
      }
      n++;
      if (n >= argc || sscanf(argv[n],"%d", &n2) != 1) {
        printf("Error: for option -b 2 integers must be given.\n");
        usage();
        exit(0);
      }
    } else if (!strcmp(argv[n],"-s")) {
      doT = 0;
    } else if (!strcmp(argv[n],"-t")) {
      n++;
      if (n >= argc || sscanf(argv[n],"%c", &typ) != 1) {
        printf("Error: for option -t a character must be given.\n");
        usage();
        exit(0);
      }
    } else if (!strcmp(argv[n],"-dbg")) {
      dbg = 1;
    } else if (!strcmp(argv[n],"-h")) {
      usage();
      exit(0);
    } else {
      printf("Error: %s is not an option.\n", argv[n]);
      usage();
      exit(0);
    }
    n++;
  }
  if (!strcmp(path,"EMPTY")) {
    printf("Error: A path must be given.\n");
    usage();
    exit(0);
  }
  if (dt0 == 0) {
    printf("Error: A date must be given.\n");
    usage();
    exit(0);
  }
  if (n1 == -99) {
    printf("Error: Bounds (n1, n2) on the date must be given.\n");
    usage();
    exit(0);
  }

  nDys = n2 - n1 + 1;
  tm = (YMD*)malloc(16*nDys);
  year = dt0 / 10000;
  tm->year = year;
  tm->month = (dt0 - year*10000) / 100;
  tm->day = dt0 % 100;
  year0 = year - 1;
  tm->year0 = year0;
  tm->yday = YearDay(*tm);
  tm->yday += n1;
  CalendarDay(tm);
  for (n = 1; n < nDys; n++) {
    (tm+n)->yday = (tm+n-1)->yday + 1;
    (tm+n)->year0 = year0;
    CalendarDay(tm+n);
  }

  nprf = 0;
  for (n = 0; n < nDys; n++) {
    if (doT) {
      sprintf(str,"tar -xvf %s/%4.4dtmp%c.tar %4.4d%2.2d%2.2dtmp.nc",path,
                     (tm+n)->year,typ,(tm+n)->year,(tm+n)->month,(tm+n)->day);
    } else {
      sprintf(str,"tar -xvf %s/%4.4dsal%c.tar %4.4d%2.2d%2.2dsal.nc",path,
                     (tm+n)->year,typ,(tm+n)->year,(tm+n)->month,(tm+n)->day);
    }
    system(str);
    nprf++;
  }

  printf("Profiles extracted: %d\n", nprf);

}

/* -------------------------------------------------------------- */

void usage()
{
  printf("Usage:\n");
  printf(" %s -p path gridFile [options]\n", prog);
  printf("   -p path      - path to daily profile files\n");
  printf("   -d date      - date (yyyymmdd)\n");
  printf("   -b n1 n2     - use data between date+n1 and date+n2\n");
  printf("   -s           - extract salinity files (df:temperature)\n");
  printf("   -t type      - O, M, S (df:O)\n");
}
