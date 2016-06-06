/*  mkEvNcr is an alternative version of mkEvNc, mkEvNcr reads an M4
    restart and mkEvNc reads a time_mean archive. Both write a single
    variance field for temperature for use by MOM3. The option -s causes
    it to also make a variance field for salinity.
    It is also possible to write a "geo3" file to aid in debugging.  */

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "netcdf.h"
#include "nc_loc.h"

void getTz(), getM4Grid(), writeFile(), writeGFile(), usage();
void nrmFile(), smthFile(), fltrFile(), cielFile(), sfcMnMx(), hFill();
void mkSEv();

mode_t mode = 0644;

float spv = 999.999, epsmx = 0.30;
float rpd = M_PI / 180.0;

int imx, jmx, kmx, kass = 30, iav = 11, jav = 5, kav = 5, ioff;
float *a, *b, *bs, *ev, *xt, *yt, *zt, *dz, *w;
double *nmlv;
int SdZ = 1, SqRt = 2, zTpr = 0, zTpr2 = 0, prdc = 0, Sfc = 1, geo3 = 0;
int Nrml = 1, ltS = -2, ltN = 2, lgW = 130, lgE = 250;
int *mask;
float VsMn = 0.30, VsMx = 0.80, gscl = 1.0, LsTpr = -90.0, LnTpr = 90.0;
float zTpf = 0.25;
int vbin[100];
float period, gciel = 1.0;
char prog[41], fileIn[91], fileOut[91], grdFile[91];
char bsn[3], ulTitle[91], urDate[91];

int dte, year, month, day, dbg = 0;

float *zm, *dzm, *wtk, *wtkp, *ab, *ap, *az;
double *buf;

/* Salinity 
    if Sm0 = 1, the Sal bkgd error will be zeroed in the mixed layer
    wherever the surface temperature is greater than isoT             */

int SVar = 0, Sm0 = 0;

/*  for salinity in MOM3:  S = (sal - 35.0) * 0.001
      sfctr = sqrt( 0.001 * 1000.0 * 1000.0 ) = 3.162e-5
    for salinity in MOM4, i.e. in psu units
      sfctr = sqrt( 0.001 ) = 3.162e-2           */
 
float Tzmn = 0.4, sfctr = 3.162e-2, isoT = 20.0, dTml = 0.5;
float *t0;
int *mlk;
char sfileOut[91];

/* netCDF */

int ncid;
int t2Vid, dtVid, xtVid, ytVid, ztVid, tVid, mlVid;
size_t start[4], count[4];
int ndims, nvars, ngatts, xdimid;
int nspdms;
float msV;
NCDIM dims;
NCVAR var;
NCATT att;

main(argc,argv)
int argc;
char *argv[];
{
  int j, n;

  strcpy(prog,argv[0]);
  strcpy(fileIn,"EMPTY");
  strcpy(fileOut,"EMPTY");
  strcpy(grdFile,"EMPTY");
  n = 1;
  while (n < argc) {
    if (!strcmp(argv[n],"-f")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -f requires a file name as a parameter.\n");
        usage();
        exit(0);
      }
      strcpy(fileIn, argv[n]);
    } else if (!strcmp(argv[n],"-o")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -o requires a file name as a parameter.\n");
        usage();
        exit(0);
      }
      strcpy(fileOut, argv[n]);
    } else if (!strcmp(argv[n],"-gr")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -gr requires a file name as a parameter.\n");
        usage();
        exit(0);
      }
      strcpy(grdFile, argv[n]);
    } else if (!strcmp(argv[n],"-d")) {
      n++;
      if (n >= argc || sscanf(argv[n],"%d",&dte) != 1) {
        printf("Error: the option -d requires an integer as parameter.\n");
        usage();
        exit(0);
      }
      year = dte/10000;
      month = (dte%10000)/100;
      day = dte%100;
    } else if (!strcmp(argv[n],"-p")) {
      prdc = 1;
    } else if (!strcmp(argv[n],"-g3")) {
      geo3 = 1;
    } else if (!strcmp(argv[n],"-k")) {
      n++;
      if (n >= argc || sscanf(argv[n],"%d",&kass) != 1) {
        printf("Error: the option -k requires an integer as parameter.\n");
        usage();
        exit(0);
      }
    } else if (!strcmp(argv[n],"-LtN")) {
      n++;
      if (n >= argc || sscanf(argv[n],"%d",&ltS) != 1) {
        printf("Error: the option -LtN requires 2 integers as parameters.\n");
        usage();
        exit(0);
      }
      n++;
      if (n >= argc || sscanf(argv[n],"%d",&ltN) != 1) {
        printf("Error: the option -LtN requires 2 integers as parameters.\n");
        usage();
        exit(0);
      }
      Nrml = 1;
    } else if (!strcmp(argv[n],"-LgN")) {
      n++;
      if (n >= argc || sscanf(argv[n],"%d",&lgW) != 1) {
        printf("Error: the option -LgN requires 2 integers as parameters.\n");
        usage();
        exit(0);
      }
      if (lgW < 0) lgW += 360;
      n++;
      if (n >= argc || sscanf(argv[n],"%d",&lgE) != 1) {
        printf("Error: the option -LgN requires 2 integers as parameters.\n");
        usage();
        exit(0);
      }
      if (lgE < 0) lgE += 360;
      Nrml = 1;
    } else if (!strcmp(argv[n],"-gS")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -gS requires a real as a parameter.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%f",&gscl);
    } else if (!strcmp(argv[n],"-gC")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -gC requires a real as a parameter.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%f",&gciel);
    } else if (!strcmp(argv[n],"-Av")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -Av requires 3 integers as parameters.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%d",&iav);
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -Av requires 3 integers as parameters.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%d",&jav);
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -Av requires 3 integers as parameters.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%d",&kav);
    } else if (!strcmp(argv[n],"-Vsfc")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -Vsfc requires 2 reals as parameters.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%f",&VsMn);
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -Vsfc requires 2 reals as parameters.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%f",&VsMx);
      Sfc = 1;
    } else if (!strcmp(argv[n],"-SdZ")) {
      SdZ = 1;
    } else if (!strcmp(argv[n],"-SqRt")) {
      n++;
      if (n >= argc || sscanf(argv[n],"%d",&SqRt) != 1) {
        printf("Error: the option -SqRt requires an integer as parameter.\n");
        usage();
        exit(0);
      }
    } else if (!strcmp(argv[n],"-zTpr")) {
      n++;
      if (n >= argc || sscanf(argv[n],"%f",&zTpf) != 1) {
        printf("Error: the option -zTpr requires a positive real <= 1.\n");
        usage();
        exit(0);
      }
      zTpr = 1;
      zTpr2 = 0;
    } else if (!strcmp(argv[n],"-zTpr2")) {
      n++;
      if (n >= argc || sscanf(argv[n],"%f",&zTpf) != 1) {
        printf("Error: the option -zTpr2 requires a positive real <= 1.\n");
        usage();
        exit(0);
      }
      zTpr2 = 1;
      zTpr = 0;
    } else if (!strcmp(argv[n],"-yTpr")) {
      n++;
      if (n >= argc || sscanf(argv[n],"%f",&LsTpr) != 1) {
        printf("Error: the option -yTpr requires 2 reals as parameters.\n");
        usage();
        exit(0);
      }
      n++;
      if (n >= argc || sscanf(argv[n],"%f",&LnTpr) != 1) {
        printf("Error: the option -yTpr requires 2 reals as parameters.\n");
        usage();
        exit(0);
      }
    } else if (!strcmp(argv[n],"-Sv")) {
      SVar = 1;
    } else if (!strcmp(argv[n],"-Sm0")) {
      Sm0 = 1;
    } else if (!strcmp(argv[n],"-so")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -so requires a file name as a parameter.\n");
        usage();
        exit(0);
      }
      strcpy(sfileOut, argv[n]);
    } else if (!strcmp(argv[n],"-st")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -st requires a real argument.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%f", &Tzmn);
    } else if (!strcmp(argv[n],"-sc")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -sc requires a positive real argument.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%f", &sfctr);
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
  if (!strcmp(fileIn,"EMPTY")) {
    printf("Error: An input file must be given.\n");
    usage();
    exit(0);
  }
  if (!strcmp(fileOut,"EMPTY")) {
    printf("Error: An output file must be given.\n");
    usage();
    exit(0);
  }
  if (!strcmp(grdFile,"EMPTY")) {
    printf("Error: An grid_spec file must be given.\n");
    usage();
    exit(0);
  }
  if (!Nrml) {
    gciel = -99.0;
  }
  if (SVar && !strcmp(sfileOut,"EMPTY")) {
    printf("Error: An output file for salinity must be given.\n");
    usage();
    exit(0);
  }
  if (!SVar) Sm0 = 0;

  getM4Grid();

  getTz();

  fltrFile();

  nrmFile();

  cielFile();

  smthFile(b);

  nrmFile();

  sfcMnMx();

  hFill();

  if (SVar) {
    mkSEv();
    smthFile(bs);
  }

  if (geo3) {
    writeGFile();
  } else {
    writeFile();
  }

  free(a);
  free(b);
  free(bs);
  free(xt);
  free(yt);
  free(zt);
  free(dz);
  free(mask);
}

/* ================================================================= */

void getM4Grid()
{
  int i, ii, j, k, n, ncid, stat;
  int ndims, xdid, ydid, zdid;
  int nvars, xvid, yvid, zvid, mvid;
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
  mask = (int*)malloc(4*imx*jmx);
  nmlv = (double*)malloc(8*imx*jmx);

  stat = nc_inq_nvars(ncid, &nvars);
  if (stat != NC_NOERR) {
    printf("Cannot get number of variables\n");
    exit(0);
  }
  xvid = -1;
  yvid = -1;
  zvid = -1;
  mvid = -1;
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
      mvid = n;
    }
  }
  if (xvid < 0 || yvid < 0 || zvid < 0 || mvid < 0) {
    printf("Cannot get variable ids: x:%d y:%d z:%d m:%d\n", xvid, yvid, zvid, mvid);
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
  stat = nc_get_var_double(ncid, mvid, nmlv);
  if (stat != NC_NOERR) {
    printf("Cannot get variable %s\n", "num_levels");
    exit(0);
  }

  stat = nc_close(ncid);
  if (stat != NC_NOERR) {
    printf("Error closing %s\n", grdFile);
  }

  w = (float*)malloc(4*imx);
  *w = *xt;
  ioff = 0;
  for (i = 1; i < imx; i++) {
    if (*(xt+i-1) < 0.0 && *(xt+i) >= 0.0) {
      ioff = i;
    }
    *(w+i) = *(xt+i);
  }
  for (i = 0; i < imx; i++) {
    ii = i + ioff;
    ii = ii < imx ? ii : ii - imx;
    *(xt+i) = *(w+ii) >= 0.0 ? *(w+ii) : *(w+ii) + 360.0;
  }

  for (j = 0; j < jmx; j++) {
    for (i = 0; i < imx; i++) {
      ii = i + ioff;
      ii = ii < imx ? ii : ii - imx;
      *(mask+i+j*imx) = (int) (*(nmlv+ii+j*imx) + 0.01);
    }
  }

  zm = (float*)malloc(4*kmx);
  dzm = (float*)malloc(4*kmx);
  wtk = (float*)malloc(4*kmx);
  wtkp = (float*)malloc(4*kmx);
  *zm = 0.0;
  for (k = 1; k < kmx; k++) {
    *(zm+k) = 0.5*(*(zt+k-1) + *(zt+k));
    *(dzm+k) = *(zt+k) - *(zt+k-1);
  }
  *dzm = *(dzm+1);
  for (k = 1; k < kmx-1; k++) {
    *(wtkp+k) = (*(zm+k+1) - *(zt+k)) / (*(zm+k+1) - *(zm+k));
    *(wtk+k) = (*(zt+k) - *(zm+k)) / (*(zm+k+1) - *(zm+k));
  }

}

/* =================================================================== */

void getTz()
{
  int i, ii, j, k, n, ka;
  float bmx, tzt;
  char str[121];
  int nc_stat;
  int dimid, varid, prd2;
  double dblTm;
  float period;

/*  printf("  Computing Tz from %s\n", fileIn); */

  nc_stat = nc_open(fileIn, NC_NOWRITE, &ncid);

  nc_stat = nc_inq(ncid, &ndims, &nvars, &ngatts, &xdimid);

  i = -1;
  j = -1;
  k = -1;
  for (dimid = 0; dimid < ndims; dimid++) {
    nc_stat = nc_inq_dim(ncid, dimid, dims.name, &dims.size);
    if (!strcasecmp(dims.name, "xaxis_1"))
      i = dims.size;
    else if (!strcasecmp(dims.name, "yaxis_1"))
      j = dims.size;
    else if (!strcasecmp(dims.name, "zaxis_1"))
      k = dims.size;
  }
  if (i != imx || j != jmx || k != kmx) {
    printf("Dimension mismatch: (%d,%d,%d) vs (%d,%d,%d)\n",i,j,k,imx,jmx,kmx);
    exit(0);
  }

  for (varid = 0; varid < nvars; varid++) {
    nc_stat = nc_inq_var(ncid, varid, var.name, &var.type, &var.ndims,
                                                      var.dims, &var.natts);
    if (!strcasecmp(var.name, "temp")) tVid = varid;
  }

  b = (float*)malloc(4*imx*jmx*kmx);
  a = (float*)malloc(4*imx*jmx*kmx);
  ab = (float*)malloc(4*kmx);
  ap = (float*)malloc(4*kmx);
  az = (float*)malloc(4*kmx);
  buf = (double*)malloc(8*imx*kmx);

  *start = 0;
  *count = 1;
  *(start+1) = 0;
  *(count+1) = kmx;
  *(start+2) = 0;
  *(count+2) = 1;
  *(start+3) = 0;
  *(count+3) = imx;
  for (j = 0; j < jmx; j++) {
    *(start+2) = j;
    nc_stat = nc_get_vara_double(ncid, tVid, start, count, buf);
    for (i = 0; i < imx; i++) {
      if (*(mask+i+j*imx)) {
        ii = i + ioff;
        ii = ii < imx ? ii : ii - imx;
        ka = *(mask+i+j*imx);
        for (k = 0; k < ka; k++)
          *(ab+k) = *(buf+ii+k*imx);

        for (k = 1; k < ka; k++) {
          *(ap+k) = (*(ab+k-1) - *(ab+k)) / *(dzm+k);
        }
        *ap = *(ap+1);
        for (k = 1; k < ka-1; k++) {
          *(az+k) = *(ap+k) * *(wtkp+k)  +  *(ap+k+1) * *(wtk+k);
        }
        *az = *(az+1);
        *(az+ka-1) = *(az+ka-2);

        for (k = 0; k < kmx; k++) {
          if (k < ka)
            *(b+i+(j+k*jmx)*imx) = *(az+k);
          else
            *(b+i+(j+k*jmx)*imx) = spv;
        }
      } else {
        for (k = 0; k < kmx; k++)
          *(b+i+(j+k*jmx)*imx) = spv;
      }
    }
  }

  if (Sm0) {

    t0 = (float*)malloc(4*imx*jmx);
    *start = 0;
    *count = 1;
    *(start+1) = 0;
    *(count+1) = 1;
    *(start+2) = 0;
    *(count+2) = jmx;
    *(start+3) = 0;
    *(count+3) = imx;
    nc_stat = nc_get_vara_float(ncid, tVid, start, count, t0);

    *start = 0;
    *count = 1;
    *(start+1) = 0;
    *(count+1) = kmx;
    *(start+2) = 0;
    *(count+2) = jmx;
    *(start+3) = 0;
    *(count+3) = imx;
    nc_stat = nc_get_vara_float(ncid, tVid, start, count, a);

    mlk = (int*)malloc(4*imx*jmx);
    for (j = 0; j < jmx; j++) {
      for (i = 0; i < imx; i++) {
        if (*(mask+i+j*imx) && *(t0+i+j*imx) > isoT) {
          tzt = *(t0+i+j*imx) - dTml;
          k = 1;
          while (*(a+i+(j+k*jmx)*imx) >= tzt && k < *(mask+i+j*imx)) {
            k++;
          }
          *(mlk+i+j*imx) = k;
        } else {
          *(mlk+i+j*imx) = 0;
        }
      }
    }
  }
  nc_stat = nc_close(ncid);

  bmx = 0.0;
  for (j = 0; j < jmx; j++) {
    for (i = 0; i < imx; i++) {
      for (k = 0; k < *(mask+i+j*imx); k++)
        bmx = *(b+i+(j+k*jmx)*imx) > bmx ? *(b+i+(j+k*jmx)*imx) : bmx;
    }
  }
  if (bmx > 0.0) {
    for (j = 0; j < jmx; j++) {
      for (i = 0; i < imx; i++) {
        for (k = 0; k < *(mask+i+j*imx); k++) {
          if (*(b+i+(j+k*jmx)*imx) < 0)
            *(b+i+(j+k*jmx)*imx) = 0.0;
          else
            *(b+i+(j+k*jmx)*imx) = *(b+i+(j+k*jmx)*imx) / bmx;
        }
      }
    }
  }

  dz = (float*)malloc(4*kmx);
  *dz = 2.0* *zt;
  for (k = 1; k < kmx; k++)
    *(dz+k) = 2.0*(*(zt+k) - *(zt+k-1)) - *(dz+k-1);

  for (j = 0; j < jmx; j++) {
    for (i = 0; i < imx; i++) {
      if (*(mask+i+j*imx) > 1)
        *(b+i+j*imx) = *(b+i+(j+jmx)*imx);
    }
  }
}

/* =================================================================== */

void nrmFile()
{
  int i, j, k;
  float lts, ltn, lgw, lge, bmx, fctr;

  bmx = 0.0;
  lts = (float)ltS;
  ltn = (float)ltN;
  lgw = (float)lgW;
  lge = (float)lgE;
  for (j = 0; j < jmx; j++) {
    if (*(yt+j) > ltn) break;
    if (*(yt+j) >= lts) {
      for (i = 0; i < imx; i++) {
        if (*(xt+i) > lge) break;
        if (*(xt+i) >= lgw) {
          for (k = 0; k < *(mask+i+j*imx); k++) {
            bmx = bmx > *(b+i+(j+k*jmx)*imx) ? bmx : *(b+i+(j+k*jmx)*imx);
          }
        }
      }
    }
  }
  fctr = gscl / bmx;
  for (j = 0; j < jmx; j++) {
    for (i = 0; i < imx; i++) {
      for (k = 0; k < *(mask+i+j*imx); k++) {
        *(b+i+(j+k*jmx)*imx) *= fctr;
      }
    }
  }
}

/* =================================================================== */

void cielFile()
{
  int i, j, k;

  if (gciel < 0.0) return;

  for (j = 0; j < jmx; j++) {
    for (i = 0; i < imx; i++) {
      for (k = 0; k < *(mask+i+j*imx); k++) {
        *(b+i+(j+k*jmx)*imx) = *(b+i+(j+k*jmx)*imx) < gciel ?
                                *(b+i+(j+k*jmx)*imx) : gciel;
      }
    }
  }
}

/* =================================================================== */

void smthFile(c)
float *c;
{
  int i, j, k, ii, i3, im, ip, jj, jm, jp, kk, km, kp, cnt;

  if (iav != 0 && jav != 0) {
    iav /= 2;
    jav /= 2;
    for (k = 0; k < kmx; k++) {
      for (j = 0; j < jmx; j++) {
        jm = (j-jav) > 0 ? (j-jav) : 0;
        jp = (j+jav) < jmx ? (j+jav) : jmx-1;
        for (i = 0; i < imx; i++) {
          if (k < *(mask+i+j*imx)) {
            if (prdc) {
              im = i-iav;
              ip = i+iav;
            } else {
              im = (i-iav) > 0 ? (i-iav) : 0;
              ip = (i+iav) < imx ? (i+iav) : imx-1;
            }
            cnt = 0;
            *(a+i+(j+k*jmx)*imx) = 0.0;
            for (jj = jm; jj <= jp; jj++) {
              for (i3 = im; i3 <= ip; i3++) {
                ii = i3;
                if (ii < 0) ii += imx;
                if (ii >= imx) ii -= imx;
                if (k < *(mask+ii+jj*imx)) {
                  *(a+i+(j+k*jmx)*imx) += *(c+ii+(jj+k*jmx)*imx);
                  cnt++;
                }
              }
            }
            if (cnt) *(a+i+(j+k*jmx)*imx) /= (float)cnt;
          } else {
            *(a+i+(j+k*jmx)*imx) = spv;
          }
        }
      }
    }
    iav = 2*iav + 1;
    jav = 2*jav + 1;
  } else {
    for (k = 0; k < kmx; k++) {
      for (j = 0; j < jmx; j++) {
        for (i = 0; i < imx; i++) {
          *(a+i+(j+k*jmx)*imx) = *(c+i+(j+k*jmx)*imx);
        }
      }
    }
  }

  if (kav > 1) {
    kav /= 2;
    for (j = 0; j < jmx; j++) {
      for (i = 0; i < imx; i++) {
        for (k = 0; k < kmx; k++) {
          if (k < *(mask+i+j*imx)) {
            km = (k-kav) > 0 ? (k-kav) : 0;
            kp = (k+kav) < kmx ? (k+kav) : kmx-1;
            cnt = 0;
            *(c+i+(j+k*jmx)*imx) = 0.0;
            for (kk = km; kk <= kp; kk++) {
              if (kk < *(mask+i+j*imx)) {
                *(c+i+(j+k*jmx)*imx) += *(a+i+(j+kk*jmx)*imx);
                cnt++;
              }
            }
            if (cnt) *(c+i+(j+k*jmx)*imx) /= (float)cnt;
          } else {
            *(c+i+(j+k*jmx)*imx) = spv;
          }
        }
      }
    }
    kav = 2*kav + 1;
  } else {
    for (k = 0; k < kmx; k++) {
      for (j = 0; j < jmx; j++) {
        for (i = 0; i < imx; i++) {
          *(c+i+(j+k*jmx)*imx) = *(a+i+(j+k*jmx)*imx);
        }
      }
    }
  }

}

/* =================================================================== */

void fltrFile()
{
  int i, j, k, kt;
  float tfz, bmx, zmx;

/* DBG the SdZ block first means variance is NOT sharply tapered, so it allows
   greater influence of observations DEEP DBG */
  if (SdZ) {
    for (j = 0; j < jmx; j++) {
      for (i = 0; i < imx; i++) {
        for (k = 0; k < *(mask+i+j*imx); k++) {
          *(b+i+(j+k*jmx)*imx) /= *(dz+k);
        }
      }
    }
  }

  if (SqRt) {
    for (j = 0; j < jmx; j++) {
      for (i = 0; i < imx; i++) {
        for (k = 0; k < *(mask+i+j*imx); k++) {
          if (*(b+i+(j+k*jmx)*imx) > 0.0) {
            *(b+i+(j+k*jmx)*imx) = sqrt(*(b+i+(j+k*jmx)*imx));
            if (SqRt >= 2) {
              *(b+i+(j+k*jmx)*imx) = sqrt(*(b+i+(j+k*jmx)*imx));
            }
          } else {
            *(b+i+(j+k*jmx)*imx) = 0.0;
          }
        }
      }
    }
  }
  if (zTpr || zTpr2) {
    zmx = *(zt+kmx-1);
    for (j = 0; j < jmx; j++) {
      for (i = 0; i < imx; i++) {
        kt = 0;
        bmx = 0.0;
        for (k = 0; k < *(mask+i+j*imx); k++) {
          if (*(b+i+(j+k*jmx)*imx) > bmx) {
            bmx = *(b+i+(j+k*jmx)*imx);
            kt = k;
          }
        }
        for (k = kt; k < *(mask+i+j*imx); k++) {
          tfz = (zmx - *(zt+k) + zTpf * (*(zt+k) - *(zt+kt))) /
                                                        (zmx - *(zt+kt));
          if (zTpr2) tfz *= tfz;
          *(b+i+(j+k*jmx)*imx) *= tfz;
        }
      }
    }
  }
/* DBG the SdZ block last means variance IS sharply tapered, so it allows
   less influence of observations deep and so more influence SHALLOW
  if (SdZ) {
    for (j = 0; j < jmx; j++) {
      for (i = 0; i < imx; i++) {
        for (k = 0; k < *(mask+i+j*imx); k++) {
          *(b+i+(j+k*jmx)*imx) /= *(dz+k);
        }
      }
    }
  }
 DBG */
}

/* =================================================================== */

void sfcMnMx()
{
  int i, j, k, klv;

  for (j = 0; j < jmx; j++) {
    for (i = 0; i < imx; i++) {
      if (*(b+i+j*imx) < VsMn) {
        klv = -1;
        for (k = 0; k < *(mask+i+j*imx); k++) {
/*  klv = k;  */
          if (*(b+i+(j+k*jmx)*imx) >= VsMn) {
            klv = k;
            break;
          }
        }
        for (k = 0; k < klv; k++) {
          *(b+i+(j+k*jmx)*imx) = VsMn;
        }
      } else if (*(b+i+j*imx) > VsMx) {
        klv = -1;
        for (k = 0; k < *(mask+i+j*imx); k++) {
          if (*(b+i+(j+k*jmx)*imx) <= VsMx) {
            klv = k;
            break;
          }
        }
        for (k = 0; k < klv; k++) {
          *(b+i+(j+k*jmx)*imx) = VsMx;
        }
      }
    }
  }
}

/* =================================================================== */

void hFill()
{
  int i, j, k, jm2, ii, jj, n;
  int ijd, im, ip, jm, jp, cn;
  float bf, eps;

  jm2 = jmx / 2;

  for (k = 0; k < kmx; k++) {

    for (n = 0; n < 100; n++) *(vbin+n) = 0;

    for (j = 0; j < jmx; j++) {
      for (i = 0; i < imx; i++) {
        if (k < *(mask+i+j*imx)) {
          n = *(b+i+(j+k*jmx)*imx) * 100.0;
          n = n < 100 ? n : 99;
          *(vbin+n) += 1;
        }
      }
    }

    i = 1;
    for (n = 1; n < 100; n++) {
      if (*(vbin+n) > *(vbin+i)) i = n;
    }
    eps = (float)i / 100.0;
    eps = eps < epsmx ? eps : epsmx;

    for (j = jm2; j >= 0; j--) {
      for (i = 0; i < imx; i++) {
        if (k < *(mask+i+j*imx) && *(b+i+(j+k*jmx)*imx) < eps) {
          cn = 0;
          bf = 0.0;
          ijd = 0;
          while (cn == 0) {
            ijd++;
            jp = j+ijd < jmx ? j+ijd : jmx-1;
            jm = j-ijd > 0 ? j-ijd : 0;
            ip = i+ijd;
            im = i-ijd;
            if (im < 0) {
              im += imx;
              for (jj = jm; jj <= jp; jj++) {
                for (ii = im; ii < imx; ii++) {
                  if (k < *(mask+ii+jj*imx) && *(b+ii+(jj+k*jmx)*imx) >= eps) {
                    bf += *(b+ii+(jj+k*jmx)*imx);
                    cn++;
                  }
                }
              }
              for (jj = jm; jj <= jp; jj++) {
                for (ii = 0; ii <= ip; ii++) {
                  if (k < *(mask+ii+jj*imx) && *(b+ii+(jj+k*jmx)*imx) >= eps) {
                    bf += *(b+ii+(jj+k*jmx)*imx);
                    cn++;
                  }
                }
              }
            } else if (ip >= imx) {
              for (jj = jm; jj <= jp; jj++) {
                for (ii = im; ii < imx; ii++) {
                  if (k < *(mask+ii+jj*imx) && *(b+ii+(jj+k*jmx)*imx) >= eps) {
                    bf += *(b+ii+(jj+k*jmx)*imx);
                    cn++;
                  }
                }
              }
              ip -= imx;
              for (jj = jm; jj <= jp; jj++) {
                for (ii = 0; ii <= ip; ii++) {
                  if (k < *(mask+ii+jj*imx) && *(b+ii+(jj+k*jmx)*imx) >= eps) {
                    bf += *(b+ii+(jj+k*jmx)*imx);
                    cn++;
                  }
                }
              }
            } else {
              for (jj = jm; jj <= jp; jj++) {
                for (ii = im; ii <= ip; ii++) {
                  if (k < *(mask+ii+jj*imx) && *(b+ii+(jj+k*jmx)*imx) >= eps) {
                    bf += *(b+ii+(jj+k*jmx)*imx);
                    cn++;
                  }
                }
              }
            }
          }
          *(b+i+(j+k*jmx)*imx) = bf / (float)cn;
        }
      }
    }

    for (j = jm2; j < jmx; j++) {
      for (i = 0; i < imx; i++) {
        if (k < *(mask+i+j*imx) && *(b+i+(j+k*jmx)*imx) < eps) {
          cn = 0;
          bf = 0.0;
          ijd = 0;
          while (cn == 0) {
            ijd++;
            jp = j+ijd < jmx ? j+ijd : jmx-1;
            jm = j-ijd > 0 ? j-ijd : 0;
            ip = i+ijd;
            im = i-ijd;
            if (im < 0) {
              im += imx;
              for (jj = jm; jj <= jp; jj++) {
                for (ii = im; ii < imx; ii++) {
                  if (k < *(mask+ii+jj*imx) && *(b+ii+(jj+k*jmx)*imx) >= eps) {
                    bf += *(b+ii+(jj+k*jmx)*imx);
                    cn++;
                  }
                }
              }
              for (jj = jm; jj <= jp; jj++) {
                for (ii = 0; ii <= ip; ii++) {
                  if (k < *(mask+ii+jj*imx) && *(b+ii+(jj+k*jmx)*imx) >= eps) {
                    bf += *(b+ii+(jj+k*jmx)*imx);
                    cn++;
                  }
                }
              }
            } else if (ip >= imx) {
              for (jj = jm; jj <= jp; jj++) {
                for (ii = im; ii < imx; ii++) {
                  if (k < *(mask+ii+jj*imx) && *(b+ii+(jj+k*jmx)*imx) >= eps) {
                    bf += *(b+ii+(jj+k*jmx)*imx);
                    cn++;
                  }
                }
              }
              ip -= imx;
              for (jj = jm; jj <= jp; jj++) {
                for (ii = 0; ii <= ip; ii++) {
                  if (k < *(mask+ii+jj*imx) && *(b+ii+(jj+k*jmx)*imx) >= eps) {
                    bf += *(b+ii+(jj+k*jmx)*imx);
                    cn++;
                  }
                }
              }
            } else {
              for (jj = jm; jj <= jp; jj++) {
                for (ii = im; ii <= ip; ii++) {
                  if (k < *(mask+ii+jj*imx) && *(b+ii+(jj+k*jmx)*imx) >= eps) {
                    bf += *(b+ii+(jj+k*jmx)*imx);
                    cn++;
                  }
                }
              }
            }
          }
          *(b+i+(j+k*jmx)*imx) = bf / (float)cn;
        }
      }
    }
  }
}

/* =================================================================== */

void mkSEv()
{
  int i, j, k, km;

  bs = (float*)malloc(4*imx*jmx*kmx);

  for (i = 0; i < imx*jmx*kmx; i++) {
    *(bs+i) = *(b+i);
  }

  for (j = 0; j < jmx; j++) {
    for (i = 0; i < imx; i++) {
      km = 0;
      for (k = 0; k < *(mask+i+j*imx); k++) {
        if (*(bs+i+(j+k*jmx)*imx) > Tzmn) km = k+1;
      }
      for (k = 0; k < km; k++) {
        *(bs+i+(j+k*jmx)*imx) = 1.0;
      }
      for (k = km; k < *(mask+i+j*imx); k++) {
        *(bs+i+(j+k*jmx)*imx) /= Tzmn;
      }

      if (Sm0) {
        if (*(mlk+i+j*imx)) {
          for (k = 0; k < *(mlk+i+j*imx); k++) {
            *(bs+i+(j+k*jmx)*imx) = 0.0;
          }
          k = *(mlk+i+j*imx) + 1;
          if (k < *(mask+i+j*imx)) {
            k--;
            *(bs+i+(j+k*jmx)*imx) = 0.5* *(bs+i+(j+(k+1)*jmx)*imx);
          }
        }
      }

      for (k = 0; k < *(mask+i+j*imx); k++) {
        *(bs+i+(j+k*jmx)*imx) *= sfctr;
      }

    }
  }

}

/* =================================================================== */

void writeFile()
{
  int i, ii, j, k, ime, jme, kme, nb, fd;

  ev = (float*)malloc(4*imx*jmx*kass);

  for (i = 0; i < imx; i++) {
    ii = i - ioff;
    ii = ii >= 0 ? ii : ii + imx;
    for (k = 0; k < kass; k++) {
      for (j = 0; j < jmx; j++) {
        *(ev+i+(j+k*jmx)*imx) = *(b+ii+(j+k*jmx)*imx);
      }
    }
  }

  fd = creat(fileOut, mode);

  nb = 3*4;
  write(fd, &nb, 4);
  write(fd, &year, 4);
  write(fd, &month, 4);
  write(fd, &day, 4);
  write(fd, &nb, 4);
  write(fd, &nb, 4);
  write(fd, &imx, 4);
  write(fd, &jmx, 4);
  write(fd, &kass, 4);
  write(fd, &nb, 4);
  nb = 4*imx*jmx*kass;
  write(fd, &nb, 4);
  write(fd, ev, nb);
  write(fd, &nb, 4);

  close(fd);

  if (SVar) {
    for (i = 0; i < imx; i++) {
      ii = i - ioff;
      ii = ii >= 0 ? ii : ii + imx;
      for (k = 0; k < kass; k++) {
        for (j = 0; j < jmx; j++) {
          *(ev+i+(j+k*jmx)*imx) = *(bs+ii+(j+k*jmx)*imx);
        }
      }
    }

    fd = creat(sfileOut, mode);

    nb = 3*4;
    write(fd, &nb, 4);
    write(fd, &year, 4);
    write(fd, &month, 4);
    write(fd, &day, 4);
    write(fd, &nb, 4);
    write(fd, &nb, 4);
    write(fd, &imx, 4);
    write(fd, &jmx, 4);
    write(fd, &kass, 4);
    write(fd, &nb, 4);
    nb = 4*imx*jmx*kass;
    write(fd, &nb, 4);
    write(fd, ev, nb);
    write(fd, &nb, 4);

    close(fd);
  }
}

/* =================================================================== */

void writeGFile()
{
  int i, j, k, ii, len, fd;
  char str[31], ch[3];

  strcpy(ulTitle, "rTv");
  sprintf(str," (A:%dx%dx%d) (", iav, jav, kav);
  strcat(ulTitle, str);
  if (SdZ || SqRt || zTpr || zTpr2) {
    if (SdZ) strcat(ulTitle, " SdZ");
    if (SqRt) strcat(ulTitle, " SqRt");
    if (zTpr) {
      sprintf(str, " zTpr:%g", zTpf);
      strcat(ulTitle, str);
    }
    if (zTpr2) {
      sprintf(str, " zTpr2:%g", zTpf);
      strcat(ulTitle, str);
    }
    if (Sfc) {
      sprintf(str," sV:%g %g", VsMn, VsMx);
      strcat(ulTitle, str);
    }
  } else {
    strcat(ulTitle, " No Filter");
  }
  strcat(ulTitle, " ) (");
  if (Nrml) {
    if (lgW > 180) lgW -= 360;
    if (lgE > 180) lgE -= 360;
    sprintf(str," N:%d,%d->%d,%d S:%g", ltS, lgW, ltN, lgE, gscl);
    if (gciel > 0.0) {
      strcat(ulTitle, str);
      sprintf(str," C:%g", gciel);
    }
  } else {
    sprintf(str," N: F:%g", gscl);
  }
  strcat(ulTitle, str);
  strcat(ulTitle, " )");

  sprintf(urDate,"%d/%2.2d/%2.2d",year,month,day);

  fd = creat(fileOut, mode);
  strcpy(bsn,"GL");
  write(fd, bsn, 2);
  strcpy(ch,"3D");
  write(fd, ch, 2);
  len = strlen(ulTitle);
  write(fd, &len, 4);
  write(fd, ulTitle, len);
  len = strlen(urDate);
  write(fd, &len, 4);
  write(fd, urDate, len);
  write(fd, &imx, 4);
  write(fd, &jmx, 4);
  write(fd, &kmx, 4);
  write(fd, xt, 4*imx);
  write(fd, yt, 4*jmx);
  write(fd, zt, 4*kmx);
  write(fd, b, 4*imx*jmx*kmx);
  close(fd);

  if (SVar) {
    *(ulTitle+1) = 'S';
    fd = creat(sfileOut, mode);
    strcpy(bsn,"GL");
    write(fd, bsn, 2);
    strcpy(ch,"3D");
    write(fd, ch, 2);
    len = strlen(ulTitle);
    write(fd, &len, 4);
    write(fd, ulTitle, len);
    len = strlen(urDate);
    write(fd, &len, 4);
    write(fd, urDate, len);
    write(fd, &imx, 4);
    write(fd, &jmx, 4);
    write(fd, &kmx, 4);
    write(fd, xt, 4*imx);
    write(fd, yt, 4*jmx);
    write(fd, zt, 4*kmx);
    write(fd, bs, 4*imx*jmx*kmx);
    close(fd);
  }
}

/* =================================================================== */

void usage()
{
  printf("Usage:\n");
  printf("  %s -f fileIn -o fileOut ...\n", prog);
  printf("    -f fileIn    - input restart netCDF file\n");
  printf("    -o fileOut   - output MOM3 file\n");
  printf("    -gr grdFile  - MOM4 grid_spec netCDF file\n");
  printf("    -d           - date (yyyymmdd) \n");
  printf("    -g3          - make output file in \'geo3\' format\n");
  printf("    -p           - input field is periodic\n");
  printf("    -k kass      - number of levels in assimilation\n");
  printf(" Filter options:\n");
  printf("    -SdZ         - scale by layer thickness\n");
  printf("    -SqRt n      - take square root of Tz (n= 1 or 2)\n");
  printf("    -Vsfc Mn Mx  - minimum/maximum values of variance at surface\n");
  printf("    -zTpr Tfctr  - taper portion of Tz below maximum (lin)\n");
  printf("    -zTpr2 Tfctr - taper portion of Tz below maximum (quad)\n");
  printf("    -yTpr Ls Ln  - taper Tz south (north) of Ls (Ln)\n");
  printf(" Smoothing options:\n");
  printf("    -Av i j k    - # grid points in averaging window\n");
  printf(" Normalization options:\n");
  printf("    -LtN ltS ltN - latitude limits for normalizing\n");
  printf("    -LgN lgW lgE - longitude limits for normalizing\n");
  printf("    -gS real     - global scale factor\n");
  printf("    -gC real     - global cieling (applied after gS)\n");
  printf(" Salinity options:\n");
  printf("    -Sv          - also do salinity variance\n");
  printf("    -so fileOut  - output MOM3 file for salinity v\n");
  printf("    -st Tzmn     - lower limit on input Tz field\n");
  printf("    -sc Sfctr    - multiply field by constant factor\n");
  printf("    -Sm0         - zero out bkgd error in mixed layer\n");
  printf("\n");
  printf("Defaults:\n");
  printf("    -k %d\n", kass);
  printf("    -SdZ\n");
  printf("    -SqRt %d\n", SqRt);
  printf("    -Vsfc %3.1f %3.1f\n", VsMn, VsMx);
  printf("    -LtN %d %d\n", ltS, ltN);
  printf("    -LgN %d %d\n", lgW, lgE);
  printf("    -gS %3.1f\n", gscl);
  printf("    -gC %3.1f\n", gciel);
  printf("    -Av %d %d %d\n", iav, jav, kav);
  printf("\n");
  printf("    -st %3.1f\n", Tzmn);
  printf("    -sc %8.3e\n", sfctr);
  printf("\n");
}
