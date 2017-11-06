/*  reads GRIB files of DAILY SST and creates files of specified length
    (N+2 days, N = length of run) daily SST records to be used by MOM4.
    the files are assumed to be on 1x1 grid.
    input naming is assumed to have the date embedded: *yyyymmdd*       */

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>

#include "netcdf.h"
#include "Tm.h"

void getM4Grid(), writeFile(), intrp2d(), intrpTri(), usage();

int icx, jcx, flip;
float *xc, *yc, *c, *c0;
float xc0, xc1, dxc, yc0, yc1, dyc;
int igx, jgx, ioff;
float *xg, *yg, *g, *g0, *wrk, xgmn, xgmx, ytr = -99.0;
double *xgd, *ygd, *tme;

int dt0 = 0, nDys = 0, nm2, jtr = 0;
int year, year0;
YMD *tm;
float tk0 = 273.15;

char *source = "Various_SST_Analyses";
char *tcal = "julian";
char tunits[31];
float mVal = 9999.0;

/* char inFile[81], path[51]; */
char inFile[181], path[121];
char outFile[51], gspFile[81];
char prog[81];

/* void main(argc,argv) */
int main(argc,argv)
int argc;
char **argv;
{
  FILE *fs;
  int i, ii, j, jj, n, len, fd;
  char *ps, str[242];

  strcpy(prog,argv[0]);
  strcpy(path,"EMPTY");
  strcpy(gspFile,"EMPTY");
  strcpy(outFile,"EMPTY");
  n = 1;
  while (n < argc) {
    if (!strcmp(argv[n],"-p")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: for option -p a path must be given.\n");
        usage();
        exit(0);
      }
      strcpy(path, argv[n]);
    } else if (!strcmp(argv[n],"-g")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: for option -g a grid/mask file must be given.\n");
        usage();
        exit(0);
      }
      strcpy(gspFile, argv[n]);
    } else if (!strcmp(argv[n],"-y")) {
      n++;
      if (n >= argc || sscanf(argv[n],"%f",&ytr) != 1) {
        printf("Error: for option -y a latitude must be given.\n");
        usage();
        exit(0);
      }
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
      sscanf(argv[n],"%d", &nDys);
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
  if (!strcmp(gspFile,"EMPTY")) {
    printf("Error: An grid/mask file must be given.\n");
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
  if (nDys == 0) {
    printf("Error: The number of days must be given.\n");
    usage();
    exit(0);
  }

  getM4Grid();

  while (*(yg+jtr) < ytr) {
    jtr++;
  }

  xgmn = *xg;
  xgmx = *(xg+igx-1);

  ioff = -1;
  for (i = 1; i < igx; i++) {
    if (*(xg+i-1) < 0.0 && *(xg+i) >= 0.0) {
      ioff = i;
      break;
    }
  }

  wrk = (float*)malloc(4*igx);
  for (i = 0; i < igx; i++) {
    ii = i+ioff < igx ? i+ioff : i+ioff-igx;
    *(wrk+i) = *(xg+ii);
  }
  for (i = 0; i < igx; i++) {
    *(xg+i) = *(wrk+i) >= 0.0 ? *(wrk+i) : *(wrk+i) + 360.0;
  }

  nm2 = nDys + 2;
  tm = (YMD*)malloc(16*nm2);
  year = dt0 / 10000;
  (tm+1)->year = year;
  (tm+1)->month = (dt0 - year*10000) / 100;
  (tm+1)->day = dt0 % 100;
  year0 = year - 1;
  (tm+1)->year0 = year0;
  (tm+1)->yday = YearDay(*(tm+1));

  tm->yday = (tm+1)->yday - 1;
  tm->year0 = year0;
  CalendarDay(tm);

  for (n = 2; n < nm2; n++) {
    (tm+n)->yday = tm->yday + n;
    (tm+n)->year0 = year0;
    CalendarDay(tm+n);
  }

  tme = (double*)malloc(8*nm2);
  g0 = (float*)malloc(4*igx*jgx);
  g = (float*)malloc(4*igx*jgx*nm2);

  sprintf(tunits,"days since %4.4d-%2.2d-%2.2d 00:00:00", tm->year, tm->month, tm->day);

  ioff = igx - ioff;
  for (n = 0; n < nm2; n++) {
    sprintf(str,"ls %s/*.%4.4d%2.2d%2.2d* > fname",path,
                              (tm+n)->year,(tm+n)->month,(tm+n)->day);
    system(str);
    fs = fopen("fname","r");
    fscanf(fs,"%s",inFile);
    fclose(fs);
    sprintf(str,"wgrib -V %s | head -7 | grep ' nx ' > finfo", inFile);
    system(str);
    fs = fopen("finfo","r");
    fgets(str,120,fs);
    fclose(fs);
    sscanf(strstr(str,"nx")+2,"%d%*s%d", &icx, &jcx);
    sprintf(str,"wgrib -V %s | head -7 | grep \" lat \" > finfo", inFile);
    system(str);
    fs = fopen("finfo","r");
    fgets(str,120,fs);
    fclose(fs);
    sscanf(strstr(str," lat ")+4,"%f%*s%f%*s%f", &yc0, &yc1, &dyc);
    sprintf(str,"wgrib -V %s | head -7 | grep \" long \" > finfo", inFile);
    system(str);
    fs = fopen("finfo","r");
    fgets(str,120,fs);
    fclose(fs);
    sscanf(strstr(str," long ")+5,"%f%*s%f%*s%f", &xc0, &xc1, &dxc);

    xc = (float*)realloc(xc,4*icx);
    yc = (float*)realloc(yc,4*jcx);
    c  = (float*)realloc(c,4*icx*jcx);
    c0 = (float*)realloc(c0,4*icx*jcx);
    for (i = 0; i < icx; i++) {
      *(xc+i) = xc0 + (float)i * dxc;
    }
    if (yc0 > yc1) {
      yc0 = yc1;
      flip = 1;
    } else {
      flip = 0;
    }
    for (j = 0; j < jcx; j++) {
      *(yc+j) = yc0 + (float)j * dyc;
    }

    sprintf(str,"wgrib -s %s | grep TMP | wgrib -i -s -ieee -nh %s", inFile, inFile);
    system(str);
    fd = open("dump", O_RDONLY);
    read(fd, c, 4*icx*jcx);
    close(fd);

    for (j = 0; j < jcx; j++) {
      if (flip) {
        jj = jcx - j - 1;
      } else {
        jj = j;
      }
      for (i = 0; i < icx; i++) {
        *(c0+i+jj*icx) = *(c+i+j*icx) - tk0;
      }
    }

    intrp2d();

    for (j = 0; j < jgx; j++) {
      for (i = 0; i < igx; i++) {
        ii = i+ioff < igx ? i+ioff : i+ioff-igx;
        *(g+i+(j+n*jgx)*igx) = *(g0+ii+j*igx);
      }
    }

    if (jtr) intrpTri(n);

    *(tme+n) = n;
  }
  sprintf(str,"rm fname finfo dump");
  system(str);

  for (i = 0; i < igx; i++) {
    ii = i+ioff < igx ? i+ioff : i+ioff-igx;
    *(wrk+i) = *(xg+ii);
  }
  for (i = 0; i < igx; i++) {
    *(xg+i) = *(wrk+i);
    if (*(xg+i) > xgmx) *(xg+i) -= 360.0;
  }

  writeFile();

}

/* ================================================================= */

void getM4Grid()
{
  int i, ii, n, ncid, stat;
  int ndims, xdid, ydid;
  int nvars, xvid, yvid, xdvid, ydvid;
  size_t len;
  char str[41];

  stat = nc_open(gspFile, NC_NOWRITE, &ncid);
  if (stat != NC_NOERR) {
    printf("Cannot open %s\n", gspFile);
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
      igx = len;
    } else if (!strncmp(str,"grid_y_T",8)) {
      ydid = n;
      jgx = len;
    }
  }
  if (xdid < 0 || ydid < 0) {
    printf("Cannot get dimension ids: x:%d y:%d\n", xdid, ydid);
    exit(0);
  }

  xg = (float*)malloc(4*igx);
  yg = (float*)malloc(4*jgx);
  xgd = (double*)malloc(8*igx*jgx);
  ygd = (double*)malloc(8*igx*jgx);

  stat = nc_inq_nvars(ncid, &nvars);
  if (stat != NC_NOERR) {
    printf("Cannot get number of variables\n");
    exit(0);
  }
  xvid = -1;
  yvid = -1;
  xdvid = -1;
  ydvid = -1;
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
    } else if (!strncmp(str,"x_T",3)) {
      xdvid = n;
    } else if (!strncmp(str,"y_T",3)) {
      ydvid = n;
    }
  }
  if (xvid < 0 || yvid < 0) {
    printf("Cannot get variable ids: x:%d y:%d\n", xvid, yvid);
    exit(0);
  }
  if (xdvid < 0 || ydvid < 0) {
    printf("Cannot get variable ids: x:%d y:%d\n", xvid, yvid);
    exit(0);
  }

  stat = nc_get_var_float(ncid, xvid, xg);
  if (stat != NC_NOERR) {
    printf("Cannot get variable %s\n", "grid_x_T");
    exit(0);
  }
  stat = nc_get_var_float(ncid, yvid, yg);
  if (stat != NC_NOERR) {
    printf("Cannot get variable %s\n", "grid_y_T");
    exit(0);
  }
  stat = nc_get_var_double(ncid, xdvid, xgd);
  if (stat != NC_NOERR) {
    printf("Cannot get variable %s\n", "x_T");
    exit(0);
  }
  stat = nc_get_var_double(ncid, ydvid, ygd);
  if (stat != NC_NOERR) {
    printf("Cannot get variable %s\n", "y_T");
    exit(0);
  }

  stat = nc_close(ncid);
  if (stat != NC_NOERR) {
    printf("Error closing %s\n", gspFile);
  }

}

/* ================================================================= */

void writeFile()
{
  int n, ncid, stat;
  int ndims, dimids[3], xdid, ydid, tdid;
  int xvid, yvid, tvid, fvid;
  size_t len, ndx;
  char str[41];

  stat = nc_create(outFile, NC_NOCLOBBER, &ncid);
  if (stat != NC_NOERR) {
    printf("Cannot create %s\n", outFile);
    exit(0);
  }

/* definitions */

  len = igx;
  stat = nc_def_dim(ncid,"grid_x_T", len, &xdid);
  if (stat != NC_NOERR) {
    printf("Cannot define dimension %s\n", "grid_x_T");
    exit(0);
  }

  len = jgx;
  stat = nc_def_dim(ncid,"grid_y_T", len, &ydid);
  if (stat != NC_NOERR) {
    printf("Cannot define dimension %s\n", "grid_y_T");
    exit(0);
  }

  len = nm2;
  stat = nc_def_dim(ncid,"time", NC_UNLIMITED, &tdid);
  if (stat != NC_NOERR) {
    printf("Cannot define dimension %s\n", "time");
    exit(0);
  }

  ndims = 1;
  stat = nc_def_var(ncid, "grid_x_T", NC_FLOAT, ndims, &xdid, &xvid);
  if (stat != NC_NOERR) {
    printf("Cannot define variable %s\n", "grid_x_T");
    exit(0);
  }

  ndims = 1;
  stat = nc_def_var(ncid, "grid_y_T", NC_FLOAT, ndims, &ydid, &yvid);
  if (stat != NC_NOERR) {
    printf("Cannot define variable %s\n", "grid_y_T");
    exit(0);
  }

  ndims = 1;
  stat = nc_def_var(ncid, "time", NC_DOUBLE, ndims, &tdid, &tvid);
  if (stat != NC_NOERR) {
    printf("Cannot define variable %s\n", "time");
    exit(0);
  }

  ndims = 3;
  *dimids = tdid;
  *(dimids+1) = ydid;
  *(dimids+2) = xdid;
  stat = nc_def_var(ncid, "temp", NC_FLOAT, ndims, dimids, &fvid);
  if (stat != NC_NOERR) {
    printf("Cannot define variable %s\n", "temp");
    exit(0);
  }

/* attributes */

  len = strlen("Nominal Longitude of T-cell center");
  stat = nc_put_att_text(ncid, xvid, "long_name", len, "Nominal Longitude of T-cell center");
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in %s\n", "long_name", "grid_x_T");
    exit(0);
  }
  len = strlen("degree_east");
  stat = nc_put_att_text(ncid, xvid, "units", len, "degree_east");
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in %s\n", "units", "grid_x_T");
    exit(0);
  }
  len = strlen("X");
  stat = nc_put_att_text(ncid, xvid, "cartesian_axis", len, "X");
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in %s\n", "cartesian_axis", "grid_x_T");
    exit(0);
  }

  len = strlen("Nominal Latitude of T-cell center");
  stat = nc_put_att_text(ncid, yvid, "long_name", len, "Nominal Latitude of T-cell center");
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in %s\n", "long_name", "grid_y_T");
    exit(0);
  }
  len = strlen("degree_north");
  stat = nc_put_att_text(ncid, yvid, "units", len, "degree_north");
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in %s\n", "units", "grid_y_T");
    exit(0);
  }
  len = strlen("Y");
  stat = nc_put_att_text(ncid, yvid, "cartesian_axis", len, "Y");
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in %s\n", "cartesian_axis", "grid_y_T");
    exit(0);
  }

  len = strlen("Time");
  stat = nc_put_att_text(ncid, tvid, "long_name", len, "Time");
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in %s\n", "long_name", "time");
    exit(0);
  }
  len = strlen(tunits);
  stat = nc_put_att_text(ncid, tvid, "units", len, tunits);
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in %s\n", "units", "time");
    exit(0);
  }
  len = strlen("T");
  stat = nc_put_att_text(ncid, tvid, "cartesian_axis", len, "T");
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in %s\n", "cartesian_axis", "time");
    exit(0);
  }
/*
  len = strlen("");
  stat = nc_put_att_text(ncid, tvid, "modulo", len, "");
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in %s\n", "modulo", "time");
    exit(0);
  }
*/
  len = strlen(tcal);
  stat = nc_put_att_text(ncid, tvid, "calendar", len, tcal);
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in %s\n", "calendar", "time");
    exit(0);
  }

  len = strlen("sea surface temperature");
  stat = nc_put_att_text(ncid, fvid, "long_name", len, "sea surface temperature");
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in %s\n", "long_name", "temp");
    exit(0);
  }
  len = strlen("deg C");
  stat = nc_put_att_text(ncid, fvid, "units", len, "deg C");
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in %s\n", "units", "temp");
    exit(0);
  }
  len = 1;
  stat = nc_put_att_float(ncid, fvid, "missing_value", NC_FLOAT, len, &mVal);
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in %s\n", "missing_value", "temp");
    exit(0);
  }
  len = 1;
  stat = nc_put_att_float(ncid, fvid, "FillValue", NC_FLOAT, len, &mVal);
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in %s\n", "FillValue", "temp");
    exit(0);
  }

  sprintf(str,"NCEP %s", source);
  len = strlen(str);
  stat = nc_put_att_text(ncid, NC_GLOBAL, "source", len, str);
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in GLOBAL\n", "source");
    exit(0);
  }
  sprintf(str,"Created %d", dt0);
  len = strlen(str);
  stat = nc_put_att_text(ncid, NC_GLOBAL, "history", len, str);
  if (stat != NC_NOERR) {
    printf("Cannot put attribute %s in GLOBAL\n", "history");
    exit(0);
  }

/* end definitions */

  stat = nc_enddef(ncid);
  if (stat != NC_NOERR) {
    printf("Cannot end def mode\n");
    exit(0);
  }

/* variables */

  stat = nc_put_var_float(ncid, xvid, xg);
  if (stat != NC_NOERR) {
    printf("Cannot put variable in %s\n", "grid_x_T");
    exit(0);
  }

  stat = nc_put_var_float(ncid, yvid, yg);
  if (stat != NC_NOERR) {
    printf("Cannot put variable in %s\n", "grid_y_T");
    exit(0);
  }

  for (n = 0; n < nm2; n++) {
    ndx = n;
    stat = nc_put_var1_double(ncid, tvid, &ndx, tme+n);
    if (stat != NC_NOERR) {
      printf("Cannot put variable in %s\n", "time");
      exit(0);
    }
  }

  stat = nc_put_var_float(ncid, fvid, g);
  if (stat != NC_NOERR) {
    printf("Cannot put variable in %s\n", "temp");
    exit(0);
  }

  stat = nc_close(ncid);
  if (stat != NC_NOERR) {
    printf("Error closing %s\n", outFile);
  }

}

/* ================================================================= */

void intrp2d()
{
  int i, j, iw, ie, js, jn, n;
  float dxw, dxe, dx, dys, dyn, dy;

  for (j = 0; j < jgx; j++) {
    js = 0;
    while (*(yc+js+1) < *(yg+j))
      js++;
    jn = js + 1;
    dyn = *(yc+jn) - *(yg+j);
    dys = *(yg+j) - *(yc+js);
    dy = *(yc+jn) - *(yc+js);

    for (i = 0; i < igx; i++) {
      if (*(xc+icx-1) < *(xg+i)) {
        iw = icx - 1;
        ie = 0;
        dxe = *(xc+ie) + 360.0 - *(xg+i);
        dxw = *(xg+i) - *(xc+iw);
        dx = *(xc+ie) + 360.0 - *(xc+iw);
      } else {
        iw = 0;
        while (*(xc+iw+1) < *(xg+i))
          iw++;
        ie = iw + 1;
        dxe = *(xc+ie) - *(xg+i);
        dxw = *(xg+i) - *(xc+iw);
        dx = *(xc+ie) - *(xc+iw);
      }
      *(g0+i+j*igx) = (dyn*(dxe * *(c0+iw+js*icx) + dxw * *(c0+ie+js*icx)) +
                                dys*(dxe * *(c0+iw+jn*icx) + dxw * *(c0+ie+jn*icx))) /
                                       (dx*dy);
    }
  }
}

/* ================================================================= */

void intrpTri(nt)
int nt;
{
  int i, j, iw, ie, js, jn, n;
  float dxw, dxe, dx, dys, dyn, dy, xgdn, cplr;

  j = jcx - 1;
  cplr = 0.0;
  for (i = 0; i < icx; i++) {
    cplr += *(c0+i+j*icx);
  }
  cplr /= (float)icx;

  for (j = jtr; j < jgx; j++) {
    for (i = 0; i < igx; i++) {
      if (*(ygd+i+j*igx) < *(yc+jcx-1)) {
        js = (int)((*(ygd+i+j*igx) - yc0)/dyc);
        jn = js + 1;
        dyn = *(yc+jn) - *(ygd+i+j*igx);
        dys = *(ygd+i+j*igx) - *(yc+js);
        dy = *(yc+jn) - *(yc+js);
        xgdn = *(xgd+i+j*igx) > 0.0 ? *(xgd+i+j*igx) : *(xgd+i+j*igx) + 360.0;
        if (xgdn > *xc) {
          iw = (int)((xgdn - xc0)/dxc);
          ie = iw + 1;
          dxe = *(xc+ie) - xgdn;
          dxw = xgdn - *(xc+iw);
          dx = *(xc+ie) - *(xc+iw);
        } else {
          iw = icx - 1;
          ie = 0;
          dxe = *(xc+ie) - xgdn;
          dxw = xgdn + 360.0 - *(xc+iw);
          dx = *(xc+ie) + 360.0 - *(xc+iw);
        }
        *(g+i+(j+nt*jgx)*igx) = (dyn*(dxe * *(c0+iw+js*icx) + dxw * *(c0+ie+js*icx)) +
                               dys*(dxe * *(c0+iw+jn*icx) + dxw * *(c0+ie+jn*icx))) /
                                         (dx*dy);
      } else {
        *(g+i+(j+nt*jgx)*igx) = cplr;
      }
    }
  }
}

/* ================================================================= */

void usage()
{
  printf("Usage:\n");
  printf(" %s -p path -g gspFile [options]\n", prog);
  printf("   -p path     - path to data files\n");
  printf("   -g gspFile  - MOM4 grid_spec netCDF file\n");
  printf("   -y ytr      - \"start\" latitude for tripolar grid\n");
  printf("   -d date     - date (yyyymmdd)\n");
  printf("   -n nDays    - number of days in run\n");
  printf("   -o outFile  - output file\n");
}
