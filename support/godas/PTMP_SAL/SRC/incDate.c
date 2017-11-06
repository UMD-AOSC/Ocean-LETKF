#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#include "Tm.h"

#define YEAR0 1950

YMD tm;
int year, month, day, year0, inc;

main(argc,argv)
int argc;
char *argv[];
{
  int n;
  char str[41];

  if (argc < 3) {
    printf("Usage:\n");
    printf("%s date inc \n", argv[0]);
    printf("  date format: yyyymmdd\n");
    printf("  inc is increment in days\n");
    exit(0);
  }

  if (sscanf(argv[1],"%4d%2d%2d", &year,&month,&day) != 3 || sscanf(argv[2],"%d",&inc) != 1) {
    printf("Usage:\n");
    printf("%s date inc \n", argv[0]);
    printf("  date format: yyyymmdd\n");
    printf("  inc is increment in days\n");
    exit(0);
  }

  tm.year0 = YEAR0;
  tm.year = year;
  tm.month = month;
  tm.day = day;

  tm.yday = YearDay(tm);
  tm.yday += inc;

  CalendarDay(&tm);
  printf("%4d%2.2d%2.2d\n", tm.year, tm.month, tm.day);
}
