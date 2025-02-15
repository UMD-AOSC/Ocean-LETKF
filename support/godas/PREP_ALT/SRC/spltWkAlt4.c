//STEVE: This file creates the GODAS format altimetry observation data file
//-------------------------------------------------------------------------
/*  spltWkAlt4 reads the ascii altimeter cycle files and writes a file
    for a given week in the format used by GODAS  */

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

#include "netcdf.h"
#include "Tm.h"

/* if std. error = 2cm, then EVR = 0.25
                   3cm             0.1111
                   4cm             0.0625
                   5cm             0.04
                  10cm             0.01  */
#define EVR   0.0625

/* if latitude bounds are set (-lt), then increasing the std error will
   smoothly reduce the influence of the ssh data approaching the bounds.  */

#define BNDW  5.0
#define YEAR0 1985
#define MAPBREAK 80.0

void getM4Grid(), getDelLst(), writeFile(), usage();
int indx(), inbnds(), unDel();

int imx, jmx, kmx, npts, ngpts, np;
float *a, *xo, *yo, *xt, *yt, *zt, *zw;
double *nLev;
int ibuf[5];
float buf[4];
YMD tm0, tmw, tm1, tc0, tc1, tm;
float t0, tw, t1, tc, hgt, evr0 = EVR, *evr, *tr;
float ax, ay, mapBreak = MAPBREAK;
int did, yday, week0, week, dow0, dow, year0 = YEAR0;
int nfl;
char *inFile[2];
char lstFile[121], delFile[121], prog[121], grdFile[121];
char path[121], pfile[121];
int nob, nx, nxx, wid;
float Lts = -90.0, Ltn = 90.0, bndw = BNDW, Ltsb, Ltnb;
int ltbnd = 0, kmin = -1;
float *xx, *xy;
int *mask;

int cnt = 0;

mode_t mode = 0644;
float spv = 999.999, eps = 0.005, obnds = 99999.99;

main(argc,argv)
int argc;
char *argv[];
{
  FILE *fs, *fso;
  int n, yr0, mn0, dy0, yr1, mn1, dy1;
  char str[81], fl[31];

  strcpy(prog,argv[0]);
  strcpy(lstFile,"EMPTY");
  strcpy(grdFile,"EMPTY");
  strcpy(delFile,"EMPTY");
  strcpy(path,"EMPTY");
  nx = 0;
  week0 = -1;
  n = 1;
  while (n < argc) {
    if (!strcmp(argv[n],"-f")) {
     // file list of ascii altimetry cycles
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -f requires a file name as a parameter.\n");
        usage();
        exit(99);
      }
      strcpy(lstFile, argv[n]);
    } else if (!strcmp(argv[n],"-p")) {
      // path to list and data files
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -p requires a path as a parameter.\n");
        usage();
        exit(99);
      }
      strcpy(path, argv[n]);
    } else if (!strcmp(argv[n],"-m")) {
      // \"t\" grid/mask file for MOM4 model
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -m requires a file name as a parameter.\n");
        usage();
        exit(99);
      }
      strcpy(grdFile, argv[n]);
    } else if (!strcmp(argv[n],"-lt")) {
      // latitude boundaries
      n++;
      if (n >= argc || sscanf(argv[n],"%f",&Lts) != 1) {
        printf("Error: the option -lt requires 2 reals as parameters.\n");
        usage();
        exit(99);
      }
      Ltsb = Lts + bndw;
      n++;
      if (n >= argc || sscanf(argv[n],"%f",&Ltn) != 1) {
        printf("Error: the option -lt requires 2 reals as parameters.\n");
        usage();
        exit(99);
      }
      Ltnb = Ltn - bndw;
      ltbnd = 1;
    } else if (!strcmp(argv[n],"-k")) {
      // minimum depth index
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -k requires an integer as a parameter.\n");
        usage();
        exit(99);
      }
      sscanf(argv[n], "%d", &kmin);
    } else if (!strcmp(argv[n],"-d")) {
      // list of deletions
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -d requires a file name as a parameter.\n");
        usage();
        exit(99);
      }
      strcpy(delFile, argv[n]);
      nx = 1;
    } else if (!strcmp(argv[n],"-w")) {
      // date of desired week (yyyymmdd)
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -w requires an integer as a parameter.\n");
        usage();
        exit(99);
      }
      sscanf(argv[n],"%d", &did);
      tmw.year = did / 10000;
      tmw.month = 1;
      tmw.day = 1;
      tmw.year0 = tmw.year;
      dow0 = DayOfWeek(tmw);
      tmw.month = (did % 10000) / 100;
      tmw.day = did % 100;
      tmw.yday = YearDay(tmw);
      week0 = 1 + (tmw.yday + dow0 - 2) / 7;
      dow = DayOfWeek(tmw);
      if (dow != 4) {
        yday = tmw.yday + 4 - dow;
        if (yday < 1) {
          tmw.year0--;
          if ((tmw.year0)%4)
            yday += 365;
          else
            if (((tmw.year0)%100) || !((tmw.year0)%400))
              yday += 366;
            else
              yday += 365;
        }
        tmw.yday = yday;
        CalendarDay(&tmw);
      }
      tmw.year0 = year0;
      yday = YearDay(tmw);
      tmw.yday = yday;
      tm0.year0 = year0;
      tm0.yday = yday - 3;
      CalendarDay(&tm0);
      tm1.year0 = year0;
      tm1.yday = yday + 3;
      CalendarDay(&tm1);
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
  if (week0 < 1) {
    printf("Error: A date indicating the desired week must be given.\n");
    usage();
    exit(99);
  }

  t0 = (float)tm0.yday;
  t1 = (float)tm1.yday + 1.0;
  tw = 0.5*(t0 + t1);
  wid = tmw.day + 100*(tmw.month + 100*(tmw.year%100));

  printf("Heights for week %2d (%2.2d/%2.2d/%2.2d) between dates: ",
               week0, tmw.year%100, tmw.month, tmw.day);
  printf("%2.2d/%2.2d/%2.2d -> %2.2d/%2.2d/%2.2d\n",
               tm0.year%100, tm0.month, tm0.day,
               tm1.year%100, tm1.month, tm1.day);
  getM4Grid();

  if (nx) getDelLst();

  sprintf(pfile,"%s/%s", path, lstFile);

  //STEVE: read in the list of dates input from an external file
  nfl = 0;
  fs = fopen(pfile, "r");
  while (fgets(str, 80, fs)) {
    sscanf(str,"%4d%*c%2d%*c%2d %4d%*c%2d%*c%2d %d %s",
            &yr0, &mn0, &dy0, &yr1, &mn1, &dy1, &n, fl);
    tc0.year = yr0;
    tc0.month = mn0;
    tc0.day = dy0;
    tc0.year0 = year0;
    tc0.yday = YearDay(tc0);
    tc1.year = yr1;
    tc1.month = mn1;
    tc1.day = dy1;
    tc1.year0 = year0;
    tc1.yday = YearDay(tc1);
    if ((tm0.yday >= tc0.yday && tm0.yday <= tc1.yday) || 
           (tm1.yday >= tc0.yday && tm1.yday <= tc1.yday)) {
      inFile[nfl] = strdup(fl);
      nfl++;
    } else if (tc0.yday > tm1.yday) {
      break;
    }
  }
  fclose(fs);

  //STEVE: if the target date (tc0) is on the list...
  if (nfl) {
    //STEVE: this is where the raw altimetry data is read in:
    ngpts = 0;
    fso = fopen("scrtch", "w"); //STEVE: used to temporarily store the adjusted lat/lon version of the data
    for (n = 0; n < nfl; n++) {
      sprintf(pfile,"%s/%s", path, inFile[n]);
      fs = fopen(pfile, "r");
      while (fgets(str, 80, fs)) {
        // Read the lat, lon, <skip>, tc (day?), and surface height (anomaly? cm?)
        sscanf(str,"%f%f%*f%f%f", &ay, &ax, &tc, &hgt);
        yday = tc;
        //STEVE: if the day is in the range of the file:
        if (yday >= tm0.yday && yday <= tm1.yday) {
          //STEVE: put the lat between 0 and 180, for the godas-format output file
          ay += 90.0; 
          //STEVE: if the logitude is not on the model grid, then add 360 to get it back on, again for the godas-format output file
          if (ax < mapBreak) ax += 360.0;
          fprintf(fso,"%8.2f%8.2f%8.2f%8.2f\n", ay, ax, tc, hgt);
          ngpts++;
        }
      }
      fclose(fs);
    }
    fclose(fso);
      
    system("sort -o scrtch scrtch");

    xo = (float*)malloc(4*ngpts);
    yo = (float*)malloc(4*ngpts);
    a = (float*)malloc(4*ngpts);
    tr = (float*)malloc(4*ngpts);
    evr = (float*)malloc(4*ngpts);
    np = 0;
    nxx = 0;
    nob = 0;

    //STEVE: read the scratch file with the reformated and coallated obs data:
    fs = fopen("scrtch", "r");
    while (fgets(str, 80, fs)) {
      sscanf(str,"%f%f%f%f", &ay, &ax, &tc, &hgt);
      ay -= 90.0;
      ax -= 360.0;
      if (inbnds(ax,ay)) {
        if (unDel(ax,ay)) {
          *(xo+np) = ax;    //STEVE: This is the lon
          *(yo+np) = ay;    //STEVE: This is the lat
          *(tr+np) = tc;    //STEVE: This is the ??? day?
          *(a+np) = hgt;    //STEVE: This is the observed height (anomaly?)
          *(evr+np) = evr0; //STEVE: This is the inverse of the variance of the obs (obs error)

          //STEVE: taper the error with an exponential function if near the latitude limits
          if (ltbnd) {
            if ((Ltsb-ay) > 0.0) {
              *(evr+np) = evr0 * exp(2.0*(ay-Ltsb)/bndw);
            } else if ((ay-Ltnb) > 0.0) {
              *(evr+np) = evr0 * exp(2.0*(Ltnb-ay)/bndw);
            }
          }
          np++;
        } else {
          nxx++;
        }
      } else {
        nob++;
      }
    }
    npts = np;

    fclose(fs);

    system("rm scrtch");

  } else {

    printf("There is no data for (%d-%2.2d-%2.2d -> %d-%2.2d-%2.2d)\n",
            tm0.year, tm0.month, tm0.day, tm1.year, tm1.month, tm1.day);
    npts = 0;
  }

  writeFile();

  printf(" Points: %d    Deletions: %d    Out of Bounds: %d\n", npts,nxx,nob);

}

/* ================================================================= */

void getM4Grid()
{
  int i, ii, n, ncid, stat;
  int ndims, xdid, ydid;
  int nvars, xvid, yvid, nvid;
  size_t len;
  char str[41];

  stat = nc_open(grdFile, NC_NOWRITE, &ncid);
  if (stat != NC_NOERR) {
    printf("Cannot open %s\n", grdFile);
    exit(0);
  }

  stat = nc_inq_ndims(ncid, &ndims);
  if (stat != NC_NOERR) {
    printf("Cannot get number of dimensions\n");
    exit(0);
  }

  xdid = -1;
  ydid = -1;
  for (n = 0; n < ndims; n++) {
    stat = nc_inq_dim(ncid, n, str, &len);
    if (stat != NC_NOERR) {
      printf("Cannot get dimensions\n");
      exit(0);
    }
    if (!strncmp(str,"grid_x_T",8)) {
      xdid = n;
      imx = len;
    } else if (!strncmp(str,"grid_y_T",8)) {
      ydid = n;
      jmx = len;
    }
  }
  if (xdid < 0 || ydid < 0) {
    printf("Cannot get dimension ids: x:%d y:%d\n", xdid, ydid);
    exit(0);
  }

  xt = (float*)malloc(4*imx);
  yt = (float*)malloc(4*jmx);
  nLev = (double*)malloc(8*imx*jmx);
  mask = (int*)malloc(4*imx*jmx);

  stat = nc_inq_nvars(ncid, &nvars);
  if (stat != NC_NOERR) {
    printf("Cannot get number of variables\n");
    exit(0);
  }
  xvid = -1;
  yvid = -1;
  nvid = -1;
  for (n = 0; n < nvars; n++) {
    stat = nc_inq_varname(ncid, n, str);
    if (stat != NC_NOERR) {
      printf("Cannot get variables\n");
      exit(0);
    }
    if (!strncmp(str,"grid_x_T",8)) {
      xvid = n;
    } else if (!strncmp(str,"grid_y_T",8)) {
      yvid = n;
    } else if (!strncmp(str,"num_levels",10)) {
      nvid = n;
    }
  }
  if (xvid < 0 || yvid < 0 || nvid < 0) {
    printf("Cannot get variable ids: x:%d y:%d nl:%d\n", xvid, yvid, nvid);
    exit(0);
  }

  stat = nc_get_var_float(ncid, xvid, xt);
  if (stat != NC_NOERR) {
    printf("Cannot get variable %s\n", "grid_x_T");
    exit(0);
  }
  stat = nc_get_var_float(ncid, yvid, yt);
  if (stat != NC_NOERR) {
    printf("Cannot get variable %s\n", "grid_y_T");
    exit(0);
  }
  stat = nc_get_var_double(ncid, nvid, nLev);
  if (stat != NC_NOERR) {
    printf("Cannot get variable %s\n", "num_levels");
    exit(0);
  }

  stat = nc_close(ncid);
  if (stat != NC_NOERR) {
    printf("Error closing %s\n", grdFile);
  }

  for (n = 0; n < imx*jmx; n++) {
    *(mask+n) = (int)(*(nLev+n) + 0.001);
  }

}

/* ================================================================= */

void getDelLst()
{
  FILE *fx;
  int n;
  char str[121];

  if ((fx = fopen(delFile,"r")) == NULL) {
    printf("Cannot open deletion file:\n  %s\n", delFile);
    exit(99);
  }
  nx = 0;
  while (fgets(str,120,fx))
    nx++;
  fclose(fx);
  xx = (float*)malloc(4*nx);
  xy = (float*)malloc(4*nx);
  fx = fopen(delFile,"r");
  for (n = 0; n < nx; n++) {
    fgets(str,120,fx);
    sscanf(str,"%f%f", xx+n, xy+n);
  }
  fclose(fx);
}

/* ================================================================= */

int inbnds(x, y)
float x, y;
{
  int i, j, ii, jj;

  i = indx(x, xt, imx);
  j = indx(y, yt, jmx);
  ii = i + 1;
  jj = j + 1;

  if (i < 0 || j < 0) return(0);

  if (*(mask+i+j*imx) < kmin || *(mask+ii+j*imx) < kmin ||
     *(mask+i+jj*imx) < kmin || *(mask+ii+jj*imx) < kmin) return(0);

  if (y >= Lts && y <= Ltn) {
    return(1);
  } else {
    return(0);
  }

  return(0);
}

/* -------------------------------------------------------------- */

int indx(p, pt, nmx)
float p, *pt;
int nmx;
{
  int n, np;

  if (p >= *pt && p <= *(pt+nmx-1)) {
    n = 0;
    while (*(pt+n) < p)
      n++;
    np = (int)((p - *(pt+n-1)) / (*(pt+n) - *(pt+n-1)) + 0.5);
    n = n - 1 + np;
    return(n);
  }
  else
    return(-1);
}

/* ================================================================= */

int unDel(x, y)
float x, y;
{
  int n;

  if (!nx) return(1);

  for (n = 0; n < nx; n++) {
    if (fabs(*(xx+n)-x) < eps && fabs(*(xy+n)-y) < eps) return(0);
  }

  return(1);
}
/* ================================================================= */

void writeFile()
{
  int n, nib, nfb, fd;
  float hour;
  char filename[41];

  tm.year0 = year0;

  if (week0 == 1)
    sprintf(filename, "swssh.%4d%2.2d", tm1.year, week0);
  else
    sprintf(filename, "swssh.%4d%2.2d", tm0.year, week0);
  fd = creat(filename, mode);

  nib = 4;
  write(fd, &nib, 4);
  write(fd, &npts, 4);
  write(fd, &nib, 4);
  nib = 20;
  nfb = 16;
  for (n = 0; n < npts; n++) {
    tm.yday = *(tr+n);
    hour = (*(tr+n) - (float)tm.yday) * 24.0;
    CalendarDay(&tm);
    write(fd, &nib, 4);
    *ibuf = tm.year;
    *(ibuf+1) = tm.month;
    *(ibuf+2) = tm.day;
    *(ibuf+3) = hour;
    *(ibuf+4) = 0;
    write(fd, ibuf, 20);
    write(fd, &nib, 4);
    write(fd, &nfb, 4);
    *buf = *(xo+n);      //STEVE: lon
    *(buf+1) = *(yo+n);  //STEVE: lat
    *(buf+2) = *(a+n);   //STEVE: altimeter data
    *(buf+3) = *(evr+n); //STEVE: presumed inverse variance (for representing the obs error)
    write(fd, buf, 16);
    write(fd, &nfb, 4);
  }
  close(fd);

}

/* ============================================================ */

void usage()
{
  printf("Usage:\n");
  printf(" %s -f lstFile -w date -m grdFile [options]\n", prog);
  printf("   -f lstFile  - file list of ascii altimetry cycles\n");
  printf("   -p path     - path to list and data files\n");
  printf("   -w date     - date of desired week (yyyymmdd)\n");
  printf("   -m grdFile  - \"t\" grid/mask file for MOM4 model\n");
  printf("   -lt LtS LtN - latitude boundaries\n");
  printf("   -k Kmin     - minimum depth index\n");
  printf("   -d delFile  - list of deletions\n");
  printf("\n");
  printf(" The format of the file list:\n");
  printf(" startDate endDate nData fileName\n");
}
