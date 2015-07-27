/*  cmbDLstPs4nc writes daily GODAS input-format profile files.
    It assumes profiles have already been interpolated to model depths.
    A profile will be excluded if it lies outside the area specified 
    by the mask (inbnds)
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "netcdf.h"
#include "Tm.h"

#define NDX 150

/* the vertical temperature gradient is rescaled to vary between 0.0 and 1.0
   so to set the std errors ultimately between 1.0 and 2.5 use
     SE0 = 1.0 and SEF = 1.5
   then 1/e.v. = 1/(se*se)
   if std error = 0.5, then EVR = 4.0
                  1.0             1.0
                  1.5             0.4444
*/

#define SE0 0.05
#define SEF 0.15

#define STDE  0.1

#define SdZ  1
#define SqRt 1
#define KAV  5
#define TZ0F 1.0

float indx();
int inbnds();
void initPrfnc(), getPrfnc(), closePrfnc();
void getM4Grid(), cmpTz(), usage();

int q = 0;

int dt0, n1, n2, nDys = 0, first = 1, doT = 1;
int year, month, day, year0, hour, minute;
YMD *tm;
int imx, jmx, kmx, *tmsk;
float *xt, *yt, *zt, *ztw, *dz;
double *nlev;
float *a, *z, *t, *tz;
float *zw, *tw, *ew;
int nb, n4 = 4, n8 = 8, n10 = 10, n24 = 24, ibuf[6];
float tz0f = TZ0F, xmins, xplus, buf[2];
int ktx;
float tzk, tz0;

char lstFile[121], prfFile[121], grdFile[121], outFile[121], wklyFile[121];
char prog[51];
int kd, kamx = 0, ecnst = 0;
float se0 = SE0, seF = SEF, stde = STDE, kav = KAV;
float zeps = 2.0, teps = 0.00005;
int nout = 0, dbg = 0;

float xmx = -999.999, xmn = 999.999;
float ymx = -999.999, ymn = 999.999;

int np, pcnt, nc1, nc2;
int ncid, ndims, xdid, ydid, zdid;
int nvars, xvid, yvid, zvid, pfvid, ptvid, lvid, hvid, fvid;
int evid, mvid; //STEVE: for reading stde and model fields from synthetic observation netcdf file
size_t len, start[4], count[4];
float xhour, mVal;
char plat[9], sid[3], ptyp[3], pname[5];
float x, y, xi, yj, se;

mode_t mode = 0644;
float spv = 999.999;

//STEVE: make conversion optional:
int do_conv_pottmp=1;
int do_comp_se=1;

/******************************************************
  MAIN
*******************************************************/

main(argc,argv)
int argc;
char *argv[];
{
  FILE *fs;
  int i, k, n, nrd, nprf, fd, fw;
  char acru[9], sid[3], dtyp[3], qkey, str[151];
  int kw, kk;

  printf("In main: cmbDLstPs4nc_obsa.c\n");

  strcpy(prog,argv[0]);
  strcpy(lstFile,"EMPTY");
  strcpy(grdFile,"EMPTY");
  strcpy(outFile,"EMPTY");
  n = 1;
  while (n < argc) {
    if (!strcmp(argv[n],"-f")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -f requires a list file as a parameter.\n");
        usage();
        exit(0);
      }
      strcpy(lstFile, argv[n]);
    } else if (!strcmp(argv[n],"-o")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -o requires a filename as a parameter.\n");
        usage();
        exit(0);
      }
      strcpy(outFile, argv[n]);
    } else if (!strcmp(argv[n],"-g")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -g requires a filename as a parameter.\n");
        usage();
        exit(0);
      }
      strcpy(grdFile, argv[n]);
    } else if (!strcmp(argv[n],"-r")) {
      n++;
      if (n >= argc || sscanf(argv[n],"%f",&se0) != 1) {
        printf("Error: the option -r requires 2 real parameters.\n");
        usage();
        exit(0);
      }
      n++;
      if (n >= argc || sscanf(argv[n],"%f",&seF) != 1) {
        printf("Error: the option -r requires 2 real parameters.\n");
        usage();
        exit(0);
      }
    } else if (!strcmp(argv[n],"-c")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        n--;
      } else {
        sscanf(argv[n],"%f",&stde);
      }
      ecnst = 1;
    } else if (!strcmp(argv[n],"-k")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -k requires an integer as a parameter.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%d", &kamx);
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
  if (!strcmp(lstFile,"EMPTY")) {
    printf("Error: A list file must be given.\n");
    usage();
    exit(0);
  }
  if (!strcmp(grdFile,"EMPTY")) {
    printf("Error: A grid/mask file must be given.\n");
    usage();
    exit(0);
  }
  if (!strcmp(outFile,"EMPTY")) {
    printf("Error: An output file must be given.\n");
    usage();
    exit(0);
  }
  if (kamx < 1) {
    printf("Error: The maximum depth level for the data must be given.\n");
    usage();
    exit(0);
  }

  getM4Grid();

  a = (float*)malloc(8*kmx);
  z = (float*)malloc(4*kmx);
  t = (float*)malloc(4*kmx);
  tz = (float*)malloc(4*kmx);
  zw = (float*)malloc(4*kmx);
  tw = (float*)malloc(4*kmx);
  ew = (float*)malloc(4*kmx);

  // STEVE: open file for reading:
  fs = fopen(lstFile, "r");
  // Count the number of days in the file list
  nDys=0;
  while (fgets(prfFile,50,fs)) {
    nDys++;
  }
  // Reset pointer to the beginning of the file
  fseek(fs, 0, SEEK_SET);

  nprf = 0;
  nrd = 0;
  fw = creat("scrtch", mode);
  // Cycle through each one of the days in the file list to read in the dates:
  for (n = 0; n < nDys; n++) {
    fgets(prfFile,50,fs);
    len = strlen(prfFile) - 1;
    prfFile[len] = '\0';

    initPrfnc();
    //STEVE: edit 
    // year, month and day are read from 'Date' field in the netcdf file,
    // but the synthetic obs don't have the date stored in the file
    sscanf(prfFile,"%8d", &dt0);
    printf("Just read date: %d\n",dt0);
    year = dt0 / 10000;
    month = (dt0 - year*10000) / 100;
    day = dt0 % 100;
    //STEVE: end edit

    nrd += pcnt;

    for (np = 0; np < pcnt; np++) {
      getPrfnc(np);

      kd = 0;
      for (k = 0; k < kmx; k++) {
        if (*(tw+k) != mVal) {
          kd = k;
        }
      }
      kd++;

      hour = (int)xhour;
      minute = (int)((xhour - (float)hour) * 60.0);

      if (x < xmins) x += 360.0;
      if (x > xplus) x -= 360.0;

      xi = indx(x, xt, imx);
      yj = indx(y, yt, jmx);
      if (xi > 0.0 && yj > 0.0) {

        kw = 0;
        for (k = 0; k < kd; k++) {
          if (*(tw+k) != mVal) {
            *(t+k) = *(tw+k);
            *(z+k) = *(zw+k);
          } else {
            *(t+k) = spv;
            *(z+k) = *(zw+k);
          }
        }

        if (do_comp_se) {
          cmpTz();
        }

        kd = kd  < kamx ? kd : kamx;

        if (kd > 1 && inbnds(xi,yj)) {

          ktx = 0;
          for (k = 0; k < kd; k++) {
            if (*(tz+k) > *(tz+ktx) && *(tz+k) != spv) {
              ktx = k;
            }
          }
          tzk = *(tz+ktx);
          tz0 = tz0f * tzk;

          for (k = 0; k < ktx; k++) {
            *(tz+k) = tz0*(*(z+ktx) - *(z+k))/(*(z+ktx)) + tzk*(*(z+k))/(*(z+ktx));
          }

          for (k = 0; k < kd; k++) {
            *(a+2*k) = *(t+k);
            if (*(t+k) != spv) {
              // specify the standard deviation (obs error)
              if (do_comp_se) {
                if (ecnst) {
                  se = stde;
                } else {
                  se = se0 + seF * *(tz+k);
                }
              } else {
                // NOTE: for the synthetic observations, this should be read in:
                se = *(ew+k);
              }
              *(a+2*k+1) = 1.0 / (se*se);
            } else {
              se = spv;
              *(a+2*k+1) = 0.0;
            }
          }

          write(fw, plat, 8);
          write(fw, ptyp, 2);
          *ibuf = year;
          *(ibuf+1) = month;
          *(ibuf+2) = day;
          *(ibuf+3) = hour;
          *(ibuf+4) = minute;
          *(ibuf+5) = kd;
          write(fw, ibuf, n24);
          *buf = x;
          *(buf+1) = y;
          write(fw, buf, n8);
          nb = 8*kd;
          write(fw, a, nb);

          nprf++;
        }
      }
    }
    closePrfnc();
  }
  close(fw);

  fd = open("scrtch",O_RDONLY);
  fw = creat(outFile, mode);
  write(fw, &n4, 4);
  write(fw, &nprf, 4);
  write(fw, &n4, 4);
  for (n = 0; n < nprf; n++) {
    read(fd, plat, 8);
    read(fd, ptyp, 2);
    write(fw, &n10, 4);
    write(fw, plat, 8);
    write(fw, ptyp, 2);
    write(fw, &n10, 4);
    read(fd, ibuf, n24);
    write(fw, &n24, 4);
    write(fw, ibuf, n24);
    write(fw, &n24, 4);
    read(fd, buf, n8);
    write(fw, &n8, 4);
    write(fw, buf, n8);
    write(fw, &n8, 4);
    nb = 8 * *(ibuf+5);
    read(fd, a, nb);
    write(fw, &nb, 4);
    write(fw, a, nb);
    write(fw, &nb, 4);
  }
  close(fd);
  close(fw);

  system("rm scrtch");
/*
  system("rm *sal.nc");
*/

  printf("Profiles read:     %d\n", nrd);
  printf("Profiles retained: %d\n", nprf);

}

/* -------------------------------------------------------------- */

void getM4Grid()
{
  int i, ii, k, n, stat;
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
  zdid = -1;
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
    } else if (!strncmp(str,"zt",2)) {
      zdid = n;
      kmx = len;
    }
  }
  if (xdid < 0 || ydid < 0 || zdid < 0) {
    printf("Cannot get dimension ids: x:%d y:%d z:%d\n", xdid, ydid, zdid);
    exit(0);
  }

  xt = (float*)malloc(4*imx);
  yt = (float*)malloc(4*jmx);
  zt = (float*)malloc(4*kmx);
  nlev = (double*)malloc(8*imx*jmx);

  stat = nc_inq_nvars(ncid, &nvars);
  if (stat != NC_NOERR) {
    printf("Cannot get number of variables\n");
    exit(0);
  }
  xvid = -1;
  yvid = -1;
  zvid = -1;
  lvid = -1;
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
    } else if (!strncmp(str,"zt",2)) {
      zvid = n;
    } else if (!strncmp(str,"num_levels",10)) {
      lvid = n;
    }
  }
  if (xvid < 0 || yvid < 0 || zvid < 0) {
    printf("Cannot get variable ids: x:%d y:%d z:%d\n", xvid, yvid, zvid);
    exit(0);
  }
  if (lvid < 0) {
    printf("Cannot get levels id\n");
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
  stat = nc_get_var_float(ncid, zvid, zt);
  if (stat != NC_NOERR) {
    printf("Cannot get variable %s\n", "zt");
    exit(0);
  }
  stat = nc_get_var_double(ncid, lvid, nlev);
  if (stat != NC_NOERR) {
    printf("Cannot get variable %s\n", "num_levels");
    exit(0);
  }

  stat = nc_close(ncid);
  if (stat != NC_NOERR) {
    printf("Error closing %s\n", grdFile);
  }

  tmsk = (int*)malloc(4*imx*jmx);
  for (n = 0; n < imx*jmx; n++) {
    *(tmsk+n) = (int)(*(nlev+n)+0.01);
  }

  dz = (float*)malloc(4*kmx);
  *dz = 2.0* *zt;
  for (k = 1; k < kmx; k++)
    *(dz+k) = 2.0*(*(zt+k) - *(zt+k-1)) - *(dz+k-1);

  xmins = *xt;
  xplus = *(xt+imx-1);
}

/* -------------------------------------------------------------- */

void initPrfnc()
{
  int k, n, ncstat;
  char str[51];

  ncstat = nc_open(prfFile, NC_NOWRITE, &ncid);
  if (ncstat != NC_NOERR) {
    printf("Cannot open %s\n", prfFile);
    printf("%s\n",nc_strerror(ncstat));
    exit(-1);
  }

  ncstat = nc_inq_ndims(ncid, &ndims);
  if (ncstat != NC_NOERR) {
    printf("Cannot get number of dimensions\n");
    printf("%s\n",nc_strerror(ncstat));
    exit(-1);
  }

  pcnt = 0;
  for (n = 0; n < ndims; n++) {
    ncstat = nc_inq_dim(ncid, n, str, &len);
    if (ncstat != NC_NOERR) {
      printf("Cannot get dimensions\n");
      printf("%s\n",nc_strerror(ncstat));
      exit(-1);
    }
    if (!strncmp(str,"grid_z",6)) {
      k = len;
    } else if (!strncmp(str,"len1",4)) {
      nc1 = len;
    } else if (!strncmp(str,"len2",4)) {
      nc2 = len;
    } else if (!strncmp(str,"count",5)) {
      pcnt = len;
    }
  }

  if (k != kmx) {
    printf("Profile levels (%d) may not match grid levels (%d)\n", k, kmx);
    exit(0);
  }

  ncstat = nc_inq_nvars(ncid, &nvars);
  if (ncstat != NC_NOERR) {
    printf("Cannot get number of variables\n");
    printf("%s\n",nc_strerror(ncstat));
    exit(-1);
  }
  for (n = 0; n < nvars; n++) {
    ncstat = nc_inq_varname(ncid, n, str);
    if (ncstat != NC_NOERR) {
      printf("Cannot get variable name\n");
      printf("%s\n",nc_strerror(ncstat));
      exit(-1);
    }
    if (!strncmp(str,"grid_z",6)) {
      zvid = n;
    } else if (!strncmp(str,"xlon",4)) {
      xvid = n;
    } else if (!strncmp(str,"ylat",4)) {
      yvid = n;
    } else if (!strncmp(str,"hour",4)) {
      hvid = n;
    } else if (!strncmp(str,"plat",4)) {
      pfvid = n;
    } else if (!strncmp(str,"ptyp",4)) {
      ptvid = n;
/*
    } else if (!strncmp(str,"sid",3)) {
      sdvid = n;
    } else if (!strncmp(str,"qkey",4)) {
      qkvid = n;
*/
    } else if (!strncmp(str,"temp",4)) {
      fvid = n;
      doT = 1;
      strcpy(pname,"temp");
    } else if (!strncmp(str,"salt",4)) {
      fvid = n;
      doT = 0;
      strcpy(pname,"salt");
    } else if (!strncmp(str,"stde",4)) {
        evid = n;
    } else if (!strncmp(str,"model",5)) {
        mvid = n;
    }
  }
  ncstat = nc_get_att_float(ncid, fvid, "missing_value", &mVal);
  if (ncstat != NC_NOERR) {
    printf("Cannot get missing_value\n");
    printf("%s\n",nc_strerror(ncstat));
    exit(-1);
  }

  ncstat = nc_get_att_text(ncid, NC_GLOBAL, "Date", str);
  if (ncstat != NC_NOERR) {
    printf("Cannot get Global Date\n");
    printf("%s\n",nc_strerror(ncstat));
    exit(-1);
  }
  sscanf(str,"%d%*c%d%*c%d",&year,&month,&day);

  if (first) {
    *start = 0;
    *count = kmx;
    ncstat = nc_get_vara_float(ncid, zvid, start, count, zw);
    if (ncstat != NC_NOERR) {
      printf("Cannot get Z variable\n");
      printf("%s\n",nc_strerror(ncstat));
      exit(-1);
    }
    first = 0;
  }
}

/* -------------------------------------------------------------- */

void getPrfnc(np)
int np;
{
  int ncstat;

  *start = np;
  *count = 1;
  ncstat = nc_get_vara_float(ncid, xvid, start, count, &x);
  if (ncstat != NC_NOERR) {
    printf("Cannot get variable in %s\n", "xlon");
    printf("%s\n",nc_strerror(ncstat));
    exit(0);
  }

  ncstat = nc_get_vara_float(ncid, yvid, start, count, &y);
  if (ncstat != NC_NOERR) {
    printf("Cannot get variable in %s\n", "ylat");
    printf("%s\n",nc_strerror(ncstat));
    exit(0);
  }

  ncstat = nc_get_vara_float(ncid, hvid, start, count, &xhour);
  if (ncstat != NC_NOERR) {
    printf("Cannot get variable in %s\n", "hour");
    printf("%s\n",nc_strerror(ncstat));
    exit(0);
  }

  *start = np;
  *count = 1;
  *(start+1) = 0;
  *(count+1) = nc2;
  ncstat = nc_get_vara_text(ncid, pfvid, start, count, plat);
  if (ncstat != NC_NOERR) {
    printf("Cannot get variable in %s\n", "plat");
    printf("%s\n",nc_strerror(ncstat));
    exit(0);
  }

  *start = np;
  *count = 1;
  *(start+1) = 0;
  *(count+1) = nc1;
  ncstat = nc_get_vara_text(ncid, ptvid, start, count, ptyp);
  if (ncstat != NC_NOERR) {
    printf("Cannot get variable in %s\n", "ptyp");
    printf("%s\n",nc_strerror(ncstat));
    exit(0);
  }

/*

  ncstat = nc_get_vara_text(ncid, sdvid, start, count, sid);
  if (ncstat != NC_NOERR) {
    printf("Cannot get variable in %s\n", "sid");
    printf("%s\n",nc_strerror(ncstat));
    exit(0);
  }

  ncstat = nc_get_vara_text(ncid, qkvid, start, count, &qkey);
  if (ncstat != NC_NOERR) {
    printf("Cannot get variable in %s\n", "qkey");
    printf("%s\n",nc_strerror(ncstat));
    exit(0);
  }

*/

  *start = np;
  *count = 1;
  *(start+1) = 0;
  *(count+1) = kmx;
  ncstat = nc_get_vara_float(ncid, fvid, start, count, tw);
  if (ncstat != NC_NOERR) {
    printf("Cannot get variable in %s\n", pname);
    printf("%s\n",nc_strerror(ncstat));
    exit(0);
  }

  // adding the standard deviation from file:
/*
  *start = np;
  *count = 1;
  *(start+1) = 0;
  *(count+1) = kmx;
  ncstat = nc_get_vara_float(ncid, evid, start, count, ew);
  if (ncstat != NC_NOERR) {
    printf("Cannot get variable in stde\n");
    printf("%s\n",nc_strerror(ncstat));
    exit(0);
  } 
*/

}

/* -------------------------------------------------------------- */

void closePrfnc()
{
  int ncstat;

  ncstat = nc_close(ncid);
  if (ncstat != NC_NOERR) {
    printf("Cannot file %s\n", prfFile);
    printf("%s\n",nc_strerror(ncstat));
    exit(0);
  }
}

/* -------------------------------------------------------------- */

void cmpTz()
{
  int k, kk, km, kp, kv2, cnt;
  float tzmn, tzmx;

  for (k = 1; k < kd-1; k++) {
    if (*(t+k-1) != spv && *(t+k+1) != spv) {
      *(tz+k) = (*(t+k-1) - *(t+k+1)) / (*(z+k+1) - *(z+k-1));
      if (*(tz+k) < 0.0) *(tz+k) = 0.0;
    } else {
      *(tz+k) = spv;
    }
  }
  *tz = spv;
  *(tz+kd-1) = spv;

  for (k = 0; k < kd; k++) {
    if (*(tz+k) == spv && *(t+k) != spv) {
      km = -1;
      for (kk = k-1; kk >= 0; kk--) {
        if (*(tz+kk) != spv) {
          km = kk;
          break;
        }
      }
      kp = kd;
      for (kk = k+1; kk < kd; kk++) {
        if (*(tz+kk) != spv) {
          kp = kk;
          break;
        }
      }
      if (km >= 0 && kp < kd) {
        if ((k-km) <= (kp-k)) {
          *(tz+k) = *(tz+km);
        } else {
          *(tz+k) = *(tz+kp);
        }
      } else if (km >= 0) {
        *(tz+k) = *(tz+km);
      } else if (kp < kd) {
        *(tz+k) = *(tz+kp);
      } else {
        *(tz+k) = spv;
      }
    }
  }

  for (k = 0; k < kd; k++) {
    if (*(tz+k) == spv && *(t+k) != spv) {
      *(tz+k) = 0.0;
    }
  }

  if (SdZ) {
    for (k = 0; k < kd; k++) {
      if (*(tz+k) != spv) {
        *(tz+k) /= *(dz+k);
      }
    }
  }

  if (SqRt) {
    for (k = 0; k < kd; k++) {
      if (*(tz+k) != spv) {
        *(tz+k) = sqrt(fabs(*(tz+k)));
      }
    }
  }

  if (kav > 1) {
    kv2 = kav / 2;
    for (k = 0; k < kd; k++) {
      if (*(tz+k) != spv) {
        km = (k-kv2) > 0 ? (k-kv2) : 0;
        cnt = 1;
        *(tw+k) = *(tz+k);
        for (kk = k-1; kk >= km; kk--) {
          if (*(tz+kk) != spv) {
            *(tw+k) += *(tz+kk);
            cnt++;
          } else {
            break;
          }
        }

        kp = (k+kv2) < kd ? (k+kv2) : kd-1;
        for (kk = k+1; kk <= kp; kk++) {
          if (*(tz+kk) != spv) {
            *(tw+k) += *(tz+kk);
            cnt++;
          } else {
            break;
          }
        }
        if (cnt) *(tw+k) /= (float)cnt;
      } else {
        *(tw+k) = spv;
      }
    }
    for (k = 0; k < kd; k++) {
      *(tz+k) = *(tw+k);
    }
  }

  tzmn = spv;
  tzmx = -spv;
  for (k = 0; k < kd; k++) {
    if (*(tz+k) != spv) {
      tzmn = tzmn < *(tz+k) ? tzmn : *(tz+k);
      tzmx = tzmx > *(tz+k) ? tzmx : *(tz+k);
    }
  }
  tzmx -= tzmn;
  if (tzmx < teps) {
    for (k = 0; k < kd; k++) {
      if (*(tz+k) != spv) {
        *(tz+k) = 0.0;
      }
    }
  } else {
    for (k = 0; k < kd; k++) {
      if (*(tz+k) != spv) {
        *(tz+k) -= tzmn;
        *(tz+k) /= tzmx;
      }
    }
  }

}

/* -------------------------------------------------------------- */

float indx(p, pt, nmx)
float p, *pt;
int nmx;
{
  int n;

  if (p >= *pt && p <= *(pt+nmx-1)) {
    n = 0;
    while (*(pt+n) < p)
      n++;
    return((float)n + (p - *(pt+n-1)) / (*(pt+n) - *(pt+n-1)));
  } else {
    return(-1.0);
  }
}

/* -------------------------------------------------------------- */

int inbnds(xi, yj)
float xi, yj;
{
  int i, j;

  i = (int)xi - 1;
  i = i < 0 ? 0 : i;
  i = i > imx-2 ? imx-2 : i;
  j = (int)yj - 1;
  j = j < 0 ? 0 : j;
  j = j > jmx-2 ? jmx-2 : j;
  if (*(tmsk+i+j*imx) || *(tmsk+i+1+j*imx) ||
                  *(tmsk+i+(j+1)*imx) || *(tmsk+i+1+(j+1)*imx))
    return(1);

  nout++;
  return(0);
}

/* -------------------------------------------------------------- */

void usage()
{
  printf("Usage:\n");
  printf(" %s -f lstFile -g gridFile [options]\n", prog);
  printf("   -f lstFile   - list of daily profile files\n");
  printf("   -g gridFile  - netCDF grid_spec file for MOM4 model\n");
  printf("   -k LevMx     - maximum level for data assimilation\n");
  printf("   -r se0 seF   - set min/max range\n");
  printf("   -c [se]      - use constant std error\n");
  printf("   -o outFile   - output file \n");

  printf(" Defaults:\n");
  printf("  SdZ:  %d\n", SdZ);
  printf("  SqRt: %d\n", SqRt);
  printf("  KAV:  %d\n", KAV);
  printf(" Controling min/max:\n");
  printf("  SE0:  %4.2f\n", SE0);
  printf("  SEF:  %4.2f\n", SEF);
  printf(" Contstant std error:\n");
  printf("  SEC:  %4.2f\n", STDE);
}
