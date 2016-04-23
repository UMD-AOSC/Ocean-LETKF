/*  reads netCDF files and creates netCDF files ndays+2 in length of
    daily SBC records to be used by MOM4. 
    the files are assumed to be annual files on the gaussian grid.
    the files are assumed to be named according to the format
    prefix_YYYY_postfix.nc, where YYYY is a 4-digit year.
    The files are assumed to be flipped, scaled and combined as
      necessary when the annual nteCDF files were made                  */

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "netcdf.h"
#include "Tm.h"

#define LYEAR 2013

void initOutFile(), transField(), closeFiles(), chngInFile(), usage();

int icx, jcx;
double *xc, *yc, *tme;
float *c;

int ncidi, ncido;
int ndims, dimids[3], xdid, ydid, tdid;
int nvars, xvid, yvid, tvid, fvid, natts;
size_t len, ndx, start[3], count[3];
char str[41], aname[31], fname[21];

int day, month, year, daym, monthm, yearm, dayp, monthp, yearp;
YMD tm;
int yday, nreci, nreco;
int nmx = 0, nm2, nm, dt0 = 0, last_year = LYEAR;

char tunits[31];

char preFix[121], postFix[121], inFile[181];
char outFile[51];
char prog[81];

int DBG = 0;

main(argc,argv)
int argc;
char **argv;
{
  int i, j, jj, m, n, fd, fd2;
  float dx;

  strcpy(prog,argv[0]);
  strcpy(preFix,"EMPTY");
  strcpy(outFile,"EMPTY");
  strcpy(postFix,"EMPTY");
  n = 1;
  while (n < argc) {
    if (!strcmp(argv[n],"-f")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: for option -f input pre- and postfix must be given.\n");
        usage();
        exit(0);
      }
      strcpy(preFix, argv[n]);
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: for option -f input pre- and postfix must be given.\n");
        usage();
        exit(0);
      }
      strcpy(postFix, argv[n]);
    } else if (!strcmp(argv[n],"-o")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: for option -o an output file must be given.\n");
        usage();
        exit(0);
      }
      strcpy(outFile, argv[n]);
    } else if (!strcmp(argv[n],"-d")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: for option -d a date must be given.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%d", &dt0);
    } else if (!strcmp(argv[n],"-n")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: for option -n the number of days must be given.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%d", &nmx);
    } else if (!strcmp(argv[n],"-y")) {
      n++;
      if (n >= argc || sscanf(argv[n],"%d", &last_year) != 1) {
        printf("Error: for option -y an integer year must be given.\n");
        usage();
        exit(0);
      }
    } else if (!strcmp(argv[n],"-dbg")) {
      DBG = 1;
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
  if (!strcmp(preFix,"EMPTY")) {
    printf("Error: An input file prefix must be given.\n");
    usage();
    exit(0);
  }
  if (!strcmp(outFile,"EMPTY")) {
    printf("Error: An output file must be given.\n");
    usage();
    exit(0);
  }
  if (dt0 == 0) {
    printf("Error: A date must be given.\n");
    usage();
    exit(0);
  }
  if (nmx == 0) {
    printf("Error: The number of days must be given.\n");
    usage();
    exit(0);
  }

  year = dt0 / 10000;
  month = (dt0 - year*10000) / 100;
  day = dt0 % 100;
  tm.year = year;
  tm.month = month;
  tm.day = day;
  if (month == 1 && day == 1) {
    tm.year0 = year - 1;
  } else {
    tm.year0 = year;
  }
  yday = YearDay(tm);
  tm.yday = yday;

  tm.yday--;
  CalendarDay(&tm);
  yearm = tm.year;
  monthm = tm.month;
  daym = tm.day;

  nm2 = nmx+2;
  tme = (double*)malloc(8*nm2);

  nm = 0;
  sprintf(tunits,"days since %4.4d-%2.2d-%2.2d 12:00:00", tm.year, tm.month, tm.day);
  printf("rr-> %s\n", tunits);
  printf(" %2.2d/%2.2d/%2.2d r:%3d rr:%d\n", tm.year, tm.month, tm.day, tm.yday, nm);
  *(tme+nm) = (double)nm;

  sprintf(inFile,"%s_%4.4d_%s.nc", preFix, tm.year, postFix);
  printf("Open: %s\n", inFile);

  initOutFile();

  nreci = tm.yday - 1;
  nreco = nm;
  transField();

  tm.yday++;
  CalendarDay(&tm);
  if (tm.year > last_year) {
    printf("Error: last year of available data exceeded (%d > %d).\n", tm.year, last_year);
    exit(0);
  }

  for (m = 1; m < nm2; m++) {

    if (tm.year > yearm) {
      sprintf(inFile,"%s_%4.4d_%s.nc", preFix, tm.year, postFix);
      printf("Open: %s\n", inFile);
      chngInFile();
      tm.year0 = tm.year;
      tm.yday = 1;
      CalendarDay(&tm);
      yearm = tm.year;
    }

    nm++;
    printf(" %2.2d/%2.2d/%2.2d r:%3d rr:%d\n", tm.year, tm.month, tm.day, tm.yday, nm);
    *(tme+nm) = (double)nm;

    nreci = tm.yday - 1;
    nreco = nm;
    transField();

    if (m+1 < nm2) {
      tm.yday++;
      CalendarDay(&tm);
      if (tm.year > last_year) {
        printf("Error: last year of available data exceeded (%d > %d).\n", tm.year, last_year);
        exit(0);
      }
    }

  }

  closeFiles();
}

/* ================================================================= */

void initOutFile()
{
  int i, n, ncstat;

  ncstat = nc_open(inFile, NC_NOWRITE, &ncidi);
  if (ncstat != NC_NOERR) {
    printf("Cannot open %s\n", inFile);
    exit(0);
  }

  ncstat = nc_create(outFile, NC_NOCLOBBER, &ncido);
  if (ncstat != NC_NOERR) {
    printf("Cannot create %s\n", outFile);
    exit(0);
  }

/* definitions */
/* dimensions  */

  ncstat = nc_inq_ndims(ncidi, &ndims);
  if (ncstat != NC_NOERR) {
    printf("Cannot get number of dimensions\n");
    exit(0);
  }

  xdid = -1;
  ydid = -1;
  tdid = -1;
  for (n = 0; n < ndims; n++) {
    ncstat = nc_inq_dim(ncidi, n, str, &len);
    if (ncstat != NC_NOERR) {
      printf("Cannot get dimensions\n");
      exit(0);
    }
    if (!strncmp(str,"XAXIS_1",7)) {
      xdid = n;
      icx = len;
      ncstat = nc_def_dim(ncido,"XAXIS_1", len, &xdid);
      if (ncstat != NC_NOERR) {
        printf("Cannot define dimension %s\n", "XAXIS_1");
        exit(0);
      }
    } else if (!strncmp(str,"YAXIS_1",7)) {
      ydid = n;
      jcx = len;
      ncstat = nc_def_dim(ncido,"YAXIS_1", len, &ydid);
      if (ncstat != NC_NOERR) {
        printf("Cannot define dimension %s\n", "YAXIS_1");
        exit(0);
      }
    } else if (!strncmp(str,"TIME",4)) {
      tdid = n;
      ncstat = nc_def_dim(ncido,"TIME", NC_UNLIMITED, &tdid);
      if (ncstat != NC_NOERR) {
        printf("Cannot define dimension %s\n", "TIME");
        exit(0);
      }
    }
  }
  if (xdid < 0 || ydid < 0 || tdid < 0) {
    printf("Cannot get dimension ids: x:%d y:%d t:%d\n", xdid, ydid, tdid);
    exit(0);
  }

  xc = (double*)malloc(8*icx);
  yc = (double*)malloc(8*jcx);
  c = (float*)malloc(4*icx*jcx);

/* variables */

  ncstat = nc_inq_nvars(ncidi, &nvars);
  if (ncstat != NC_NOERR) {
    printf("Cannot get number of variables\n");
    exit(0);
  }

  xvid = -1;
  yvid = -1;
  tvid = -1;
  fvid = -1;
  for (n = 0; n < nvars; n++) {
    ncstat = nc_inq_varname(ncidi, n, str);
    if (ncstat != NC_NOERR) {
      printf("Cannot get variables\n");
      exit(0);
    }
    ncstat = nc_inq_varnatts(ncidi, n, &natts);
    if (ncstat != NC_NOERR) {
      printf("Cannot get number of atts for %s\n", str);
      exit(0);
    }
    if (!strncmp(str,"XAXIS_1",7)) {
      xvid = n;
      ndims = 1;
      ncstat = nc_def_var(ncido, "XAXIS_1", NC_DOUBLE, ndims, &xdid, &xvid);
      if (ncstat != NC_NOERR) {
        printf("Cannot define variable %s\n", "XAXIS_1");
        exit(0);
      }
      for (i = 0; i < natts; i++) {
        ncstat = nc_inq_attname(ncidi, xvid, i, aname);
        if (ncstat != NC_NOERR) {
          printf("Cannot get att name from %s\n", "XAXIS_1");
          exit(0);
        }
        ncstat = nc_copy_att(ncidi, xvid, aname, ncido, xvid);
        if (ncstat != NC_NOERR) {
          printf("Cannot copy att %s for %s\n", aname, "XAXIS_1");
          exit(0);
        }
      }
    } else if (!strncmp(str,"YAXIS_1",7)) {
      yvid = n;
      ndims = 1;
      ncstat = nc_def_var(ncido, "YAXIS_1", NC_DOUBLE, ndims, &ydid, &yvid);
      if (ncstat != NC_NOERR) {
        printf("Cannot define variable %s\n", "YAXIS_1");
        exit(0);
      }
      for (i = 0; i < natts; i++) {
        ncstat = nc_inq_attname(ncidi, yvid, i, aname);
        if (ncstat != NC_NOERR) {
          printf("Cannot get att name from %s\n", "YAXIS_1");
          exit(0);
        }
        ncstat = nc_copy_att(ncidi, yvid, aname, ncido, yvid);
        if (ncstat != NC_NOERR) {
          printf("Cannot copy att %s for %s\n", aname, "YAXIS_1");
          exit(0);
        }
      }
    } else if (!strncmp(str,"TIME",4)) {
      tvid = n;
      ndims = 1;
      ncstat = nc_def_var(ncido, "TIME", NC_DOUBLE, ndims, &tdid, &tvid);
      if (ncstat != NC_NOERR) {
        printf("Cannot define variable %s\n", "TIME");
        exit(0);
      }
      for (i = 0; i < natts; i++) {
        ncstat = nc_inq_attname(ncidi, tvid, i, aname);
        if (ncstat != NC_NOERR) {
          printf("Cannot get att name from %s\n", "TIME");
          exit(0);
        }
        if (!strncmp(aname,"units",5)) {
          len = strlen(tunits);
          ncstat = nc_put_att_text(ncido, tvid, "units", len, tunits);
          if (ncstat != NC_NOERR) {
            printf("Cannot put att %s for %s\n", "units", "TIME");
            exit(0);
          }
        } else {
          ncstat = nc_copy_att(ncidi, tvid, aname, ncido, tvid);
          if (ncstat != NC_NOERR) {
            printf("Cannot copy att %s for %s\n", aname, "TIME");
            exit(0);
          }
        }
      }
    } else {
      strcpy(fname,str);
      fvid = n;
      ndims = 3;
      *dimids = tdid;
      *(dimids+1) = ydid;
      *(dimids+2) = xdid;
      ncstat = nc_def_var(ncido, fname, NC_FLOAT, ndims, dimids, &fvid);
      if (ncstat != NC_NOERR) {
        printf("Cannot define variable %s\n", fname);
        exit(0);
      }
      for (i = 0; i < natts; i++) {
        ncstat = nc_inq_attname(ncidi, fvid, i, aname);
        if (ncstat != NC_NOERR) {
          printf("Cannot get att name from %s\n", fname);
          exit(0);
        }
        ncstat = nc_copy_att(ncidi, fvid, aname, ncido, fvid);
        if (ncstat != NC_NOERR) {
          printf("Cannot copy att %s for %s\n", aname, fname);
          exit(0);
        }
      }
    }
  }

  ncstat = nc_inq_natts(ncidi, &natts);
  if (ncstat != NC_NOERR) {
    printf("Cannot get number of global atts\n");
    exit(0);
  }
  for (i = 0; i < natts; i++) {
    ncstat = nc_inq_attname(ncidi, NC_GLOBAL, i, aname);
    if (ncstat != NC_NOERR) {
      printf("Cannot get global att\n");
      exit(0);
    }
    ncstat = nc_copy_att(ncidi, NC_GLOBAL, aname, ncido, NC_GLOBAL);
    if (ncstat != NC_NOERR) {
      printf("Cannot copy global att %s\n", aname);
      exit(0);
    }
  }

/* end definitions */

  ncstat = nc_enddef(ncido);
  if (ncstat != NC_NOERR) {
    printf("Cannot end def mode\n");
    exit(0);
  }

/* variables */

  ncstat = nc_get_var_double(ncidi, xvid, xc);
  if (ncstat != NC_NOERR) {
    printf("Cannot get variable %s\n", "XAXIS_1");
    exit(0);
  }

  ncstat = nc_put_var_double(ncido, xvid, xc);
  if (ncstat != NC_NOERR) {
    printf("Cannot put variable %s\n", "XAXIS_1");
    exit(0);
  }

  ncstat = nc_get_var_double(ncidi, yvid, yc);
  if (ncstat != NC_NOERR) {
    printf("Cannot get variable %s\n", "YAXIS_1");
    exit(0);
  }

  ncstat = nc_put_var_double(ncido, yvid, yc);
  if (ncstat != NC_NOERR) {
    printf("Cannot put variable %s\n", "YAXIS_1");
    exit(0);
  }

}

/* ================================================================= */

void transField()
{
  int ncstat;

  *start = nreco;
  *count = 1;
  ncstat = nc_put_vara_double(ncido, tvid, start, count, tme+nreco);
  if (ncstat != NC_NOERR) {
    printf("Cannot put variable in %s\n", "TIME");
    exit(0);
  }

  *start = nreci;
  *count = 1;
  *(start+1) = 0;
  *(count+1) = jcx;
  *(start+2) = 0;
  *(count+2) = icx;
  ncstat = nc_get_vara_float(ncidi, fvid, start, count, c);
  if (ncstat != NC_NOERR) {
    printf("Cannot put variable %s\n", fname);
    exit(0);
  }

  *start = nreco;
  *count = 1;
  *(start+1) = 0;
  *(count+1) = jcx;
  *(start+2) = 0;
  *(count+2) = icx;
  ncstat = nc_put_vara_float(ncido, fvid, start, count, c);
  if (ncstat != NC_NOERR) {
    printf("Cannot put variable %s\n", fname);
    exit(0);
  }
}

/* ================================================================= */

void chngInFile()
{
  int ncstat;

  ncstat = nc_close(ncidi);
  if (ncstat != NC_NOERR) {
    printf("Error closing %s\n", inFile);
  }

  ncstat = nc_open(inFile, NC_NOWRITE, &ncidi);
  if (ncstat != NC_NOERR) {
    printf("Cannot open %s\n", inFile);
    exit(0);
  }
}

/* ================================================================= */

void closeFiles()
{
  int ncstat;

  ncstat = nc_close(ncidi);
  if (ncstat != NC_NOERR) {
    printf("Error closing %s\n", inFile);
  }

  ncstat = nc_close(ncido);
  if (ncstat != NC_NOERR) {
    printf("Error closing %s\n", outFile);
  }

}

/* ================================================================= */

void usage()
{
  printf("Usage:\n");
  printf(" %s -f preFix [options]\n", prog);
  printf("   -f preFix postFix  - input path/prefix postfix\n");
  printf("   -d date            - date (yyyymmdd)\n");
  printf("   -n ndays           - number of days\n");
  printf("   -o outFile         - output file\n");
  printf("   -y last_year       - last year of available data\n");
}
