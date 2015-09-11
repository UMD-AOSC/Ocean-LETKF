/*  STEVE: this version converted from synthetic profiles,
           performs pot temp conversion and computes stanard errors
 
    cmbDysTh4nc reads daily profile files, combines them and
    writes them in NEW M4 assimilation format.
    It reads T and S profiles (and does NOT convert T to Pot.T.)
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

#define NDX 1500

/* the vertical temperature gradient is rescaled to vary between 0.0 and 1.0
   so to set the std errors ultimately between 1.0 and 2.5 use
     SE0 = 1.0 and SEF = 1.5
   then EVR = 1/(se*se)
   if std error = 0.5, then EVR = 4.0
                  1.0             1.0
                  1.5             0.4444
*/

#define SE0 1.0
#define SEF 1.5

#define SdZ  1
#define SqRt 1
#define KAV  5

float indx(), atg(), press(), tempa(), theta(), zeta();
int inbnds();
void initPrfnc(), getHeadnc(), getPrfnc(), closePrfnc();
void getM4Grid(), Sfill(), cmpTz(), usage();

int q = 0;

int dt0, nDys, first = 1, doT = 1;
int year, month, day, year0, hour, minute;
YMD *tm;
int imx, jmx, kmx, *tmsk;
float *xt, *yt, *zt, *ztw, *dz;
double *nlev;
float *a, *z, *t, *tz;
float *zw, *tw, *sw, *ew;
int nb, n4 = 4, n8 = 8, n10 = 10, n24 = 24, ibuf[6];
float xmins, xplus, buf[2];

char prfTFile[121], prfSFile[121], grdFile[121], outFile[121];
char tFlst[31], sFlst[31], tlist[31], slist[31];
char prog[51];
int kd, kdt, kamx = 0;
float se0 = SE0, seF = SEF, kav = KAV;
float zeps = 2.0, teps = 0.00005;
int nout = 0, dbg = 0;

float xmx = -999.999, xmn = 999.999;
float ymx = -999.999, ymn = 999.999;

int np, npt, nps, pcnt, pcntt, nc1, nc2;
int ncid, ndims, xdid, ydid, zdid;
int ncidt, ncids;
int nvars, xvid, yvid, zvid, pfvid, ptvid, lvid, hvid, tvid, svid;
int evid, mvid; //STEVE: for reading stde and model fields from synthetic observation netcdf file
size_t len, start[4], count[4];
float xhour, mtVal, msVal; // mtVal and msVal are missing values for temp and salt, respectively
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
  FILE *fsft, *fsfs, *fslt, *fsls, *fslts;
  int done, ntwo, nsng, nt, ns;
  int i, k, n, sflg, nl, nlt, nrd, nprf, ndf, fd, fw, sln,  sln0;
  float p, p0;
  char acru[9], sid[3], dtyp[3], qkey, str[151];
  char tline[81], sline[81];
  int kw, kk;

  strcpy(prog,argv[0]);
  strcpy(grdFile,"EMPTY");
  strcpy(outFile,"EMPTY");

  //STEVE: Read all command line arguments:
  n = 1;
  while (n < argc) {
    // -f tFlst sFlst - temperature and salinity file lists
    if (!strcmp(argv[n],"-f")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -f requires 2 filenames as parameters.\n");
        usage();
        exit(0);
      }
      strcpy(tFlst, argv[n]);
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -f requires 2 filenames as parameters.\n");
        usage();
        exit(0);
      }
      strcpy(sFlst, argv[n]);
    // -o outFile   - output file
    } else if (!strcmp(argv[n],"-o")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -o requires a filename as a parameter.\n");
        usage();
        exit(0);
      }
      strcpy(outFile, argv[n]);
    //  -g gridFile  - netCDF grid_spec file for MOM4 model
    } else if (!strcmp(argv[n],"-g")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -g requires a filename as a parameter.\n");
        usage();
        exit(0);
      }
      strcpy(grdFile, argv[n]);
    // -r se0 seF   - set min/max range
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
    //  -k LevMx     - maximum level for data assimilation
    } else if (!strcmp(argv[n],"-k")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -k requires an integer as a parameter.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%d", &kamx);
    // set debug flag
    } else if (!strcmp(argv[n],"-dbg")) {
      dbg = 1;
    // list usage help
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
  if (!strcmp(tFlst,"EMPTY")) {
    printf("Error: 2 file lists (T and S) must be given.\n");
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
  
  // STEVE: open file for reading:
  fsft = fopen(tFlst, "r");
  nDys = 0;
  // Count the number of days in the file list
  while (fgets(prfTFile,120,fsft)) {
    nDys++;
  }
  // Reset pointer to the beginning of the file
  fseek(fsft, 0, SEEK_SET);
  tm = (YMD*)malloc(16*nDys);

  // Cycle through each one of the days in the file list to read in the dates:
  for (n = 0; n < nDys; n++) {
    fgets(prfTFile,120,fsft);
    sln = strlen(prfTFile) - 1;
    prfTFile[sln] = '\0';
//  sscanf(prfTFile,"%d", &dt0);
    sscanf(prfTFile,"%8d", &dt0);
    printf("Just read date: %d\n",dt0);
    year = dt0 / 10000;
    month = (dt0 - year*10000) / 100;
    day = dt0 % 100;
/*
    printf("n=%d\n",n);
    printf("year=%d\n",year);
    printf("month=%d\n",month);
    printf("day=%d\n",day);
*/
    if (!n) year0 = year;
    (tm+n)->year = year;
    (tm+n)->month = month;
    (tm+n)->day = day;
    (tm+n)->year0 = year0;
    (tm+n)->yday = YearDay(*(tm+n));
  }

  // Read in mom4 grid file  
  getM4Grid();

  // Allocate primary arrays for this routine:
  a = (float*)malloc(8*kmx);
  z = (float*)malloc(4*kmx);
  t = (float*)malloc(4*kmx);
  tz = (float*)malloc(4*kmx);
  zw = (float*)malloc(4*kmx);
  tw = (float*)malloc(4*kmx);
  sw = (float*)malloc(4*kmx);
  ew = (float*)malloc(4*kmx); //STEVE: added for error supplied via netcdf file

  nprf = 0;
  nrd = 0;
  // open/create a new scratch file
  fw = creat("scrtch", mode);

  printf("    Date     M/UnM\n");

  //STEVE: step through each day
  for (n = 0; n < nDys; n++) {

    //STEVE: get name of temperature profiles netcdf file
    sprintf(prfTFile,"%4.4d%2.2d%2.2dtmp.nc",(tm+n)->year,(tm+n)->month,(tm+n)->day);

    // Initialize netcdf file with temp profiles (e.g. open, get general info)
//    printf("DEBUG:: Calling :: initPrfnc('T'); ...\n");
    initPrfnc('T');
    //STEVE: edit 
    // year, month and day are read from 'Date' field in the netcdf file,
    // but the synthetic obs don't have the date stored in the file
    year=(tm+n)->year;
    month=(tm+n)->month;
    day=(tm+n)->day;
    //STEVE: end edit
    nrd += pcnt;
    pcntt = pcnt; // establish the observation count for this file

    // Create a more detailed list of observations
//    printf("DEBUG:: Calling :: fslt = fopen(\"tlist\",\"w\"); ...\n");
    fslt = fopen("tlist","w");
    for (np = 0; np < pcnt; np++) {
//      printf("DEBUG:: Calling :: getHeadnc('T', %5d); ...\n",np);
      getHeadnc('T', np);
      hour = (int)xhour;
      minute = (int)((xhour - (float)hour) * 60.0);
      if (x < 0.0) x += 360.0;
      y += 90.0;
      // For each observation, output the date, time, lat, lon, and index in a text file
      fprintf(fslt,"%4d%2.2d%2.2d %2.2d%2.2d %7.3f %8.3f 0 %5d\n", year, month, day, hour, minute, y, x, np);
      //STEVE: output to screen just to verify date is correct...
      //printf("%4d%2.2d%2.2d %2.2d%2.2d %7.3f %8.3f 0 %5d\n", year, month, day, hour, minute, y, x, np);
    }
    fclose(fslt);
    // Write sorted concatenation of all files (and output to overwrite original file)
//    printf("DEBUG:: Calling :: system(\"sort -o tlist tlist\"); ...\n");
    system("sort -o tlist tlist");

    sprintf(prfSFile,"%4.4d%2.2d%2.2dsal.nc",(tm+n)->year,(tm+n)->month,(tm+n)->day);

    // Initialize netcdf file with salt profiles (e.g. open, get general info)
//    printf("DEBUG:: Calling :: initPrfnc('S'); ...\n");
    initPrfnc('S');

//    printf("DEBUG:: Calling :: fsls = fopen(\"slist\",\"w\"); ...\n");
    fsls = fopen("slist","w");
    for (np = 0; np < pcnt; np++) {
//      printf("DEBUG:: Calling :: getHeadnc('S', %5d); ...\n",np);
      getHeadnc('S', np);
      hour = (int)xhour;
      minute = (int)((xhour - (float)hour) * 60.0);
      if (x < 0.0) x += 360.0;
      y += 90.0;
      fprintf(fsls,"%4d%2.2d%2.2d %2.2d%2.2d %7.3f %8.3f 1 %5d\n", year, month, day, hour, minute, y, x, np);
      //STEVE: output to screen just to verify date is correct...
//      printf("%4d%2.2d%2.2d %2.2d%2.2d %7.3f %8.3f 0 %5d\n", year, month, day, hour, minute, y, x, np);
    }
    fclose(fsls);
//    printf("DEBUG:: Calling :: system(\"sort -o slist slist\"); ...\n");
    system("sort -o slist slist");

    // Combine temp and salt observation lists into one file
//    printf("DEBUG:: Calling :: system(\"sort -o srtlist tlist slist\"); ...\n");
    system("sort -o srtlist tlist slist");

    if (pcnt != pcntt) {
      printf("Profile count difference: T=%d  S=%d\n", pcntt, pcnt);
    }

    pcnt = 0;
    ndf = 0;
    // Open list of all observations for reading
//    printf("DEBUG:: Calling :: fslt = fopen(\"srtlist\",\"r\"); ...\n");
    fslt = fopen("srtlist","r");
    // Write to new tslist
    fslts = fopen("tslist","w");
    done = 0;
    if (! fgets(tline,80,fslt)) done = 1;
    if (! fgets(sline,80,fslt)) done = 1;
    while (!done) {
//      printf("DEBUG:: pcnt = %d, ndf = , %d, done = %d \n",pcnt,ndf,done);
//      printf("DEBUG:: 0 tline = %s\n",tline);
      sln = 38; //strlen(tline) - 1;
      sln0=strlen(tline);
//      printf("DEBUG:: strlen(tline) = %d\n",sln0);
//      printf("DEBUG:: sln=%d\n",sln);
//      printf("DEBUG:: tline[sln] = %c\n",tline[sln]);
      if (sln < 38) {
        printf("DEBUG:: tline = %s\n",tline);
        printf("DEBUG:: strlen(tline) = %d\n",sln0);
        sln=38;
        int ii;
        for(ii=0;ii<sln;ii++) {
          printf("DEBUG:: tline[%d] = %c\n",ii,tline[ii]);
        }
      } //STEVE:: DEBUG
//      printf("DEBUG:: tline[sln] = %c\n",tline[sln]);
      tline[sln] = '\0';
//      printf("DEBUG:: 1a tline = %s\n",tline);
      sscanf(tline+30,"%d",&nt);
//      printf("DEBUG:: nt = %d \n", nt);
//      printf("DEBUG:: 1b tline = %s\n",tline);
//      printf("DEBUG:: 0 sline = %s\n",sline);
      sln = 38; //strlen(sline) - 1;
      sline[sln] = '\0';
//      printf("DEBUG:: 1a sline = %s\n",sline);
      sscanf(sline+30,"%d",&ns);
//      printf("DEBUG:: nt = %d \n", ns);
//      printf("DEBUG:: 1b sline = %s\n",sline);
      if (nt == 0 && ns == 1) {
//        printf("DEBUG:: 2a tline = %s\n",tline);
//        printf("DEBUG:: 2a sline = %s\n",sline);
        fprintf(fslts,"%s\n",tline);
        fprintf(fslts,"%s\n",sline);
        if (! fgets(tline,80,fslt)) done = 1;
        if (! fgets(sline,80,fslt)) done = 1;
        pcnt++;
      } else if (nt == 1 && ns == 0) {
//        printf("DEBUG:: 2b tline = %s\n",tline);
//        printf("DEBUG:: 2b sline = %s\n",sline);
        nt = ns;
        strcpy(tline, sline);
        if (! fgets(sline,80,fslt)) done = 1;
        ndf++;
      } else if (nt == 0 && ns == 0) {
//        printf("DEBUG:: 2c tline = %s\n",tline);
//        printf("DEBUG:: 2c sline = %s\n",sline);
        nt = ns;
        strcpy(tline, sline);
        if (! fgets(sline,80,fslt)) done = 1;
        ndf++;
      } else if (nt == 1 && ns == 1) {
//        printf("DEBUG:: 2d tline = %s\n",tline);
//        printf("DEBUG:: 2d sline = %s\n",sline);
        if (! fgets(tline,80,fslt)) done = 1;
        if (! fgets(sline,80,fslt)) done = 1;
        ndf++;
      }
    }
    fclose(fslt);
    fclose(fslts);
    
    printf("%4d-%2.2d-%2.2d : %d/%d\n", year, month, day, pcnt, ndf);

    // Open tslist just created for reading:
//    printf("DEBUG :: Opening tslist file...\n");
    fslts = fopen("tslist","r");

    for (np = 0; np < pcnt; np++) {
      fgets(tline,80,fslts);
      sln = 38; //strlen(tline) - 1;
      tline[sln] = '\0';
      sscanf(tline+33,"%d", &npt);

      // get temperature profile
//      printf("DEBUG:: Calling :: getPrfnc('T', %5d); ...\n",npt);
      getPrfnc('T', npt);
      fgets(sline,80,fslts);
      sln = 38; //strlen(sline) - 1;
      sline[sln] = '\0';
      sscanf(sline+33,"%d", &nps);

      // get salinity profile
//      printf("DEBUG:: Calling :: getPrfnc('S', %5d); ...\n",nps);
      getPrfnc('S', nps);

      // count the number of levels with data (non-missing values)
      kd = 0;
      sflg = 0;
      for (k = 0; k < kmx; k++) {
        if (*(tw+k) != mtVal) {
          kd = k;
          if (*(sw+k) == msVal) sflg++;
        }
      }
      kd++; // deepest level with observation data

      // if some temperature fields don't have salinity, then fill in
      if (sflg) {
        Sfill();
      }

      hour = (int)xhour;
      minute = (int)((xhour - (float)hour) * 60.0);

      if (x < xmins) x += 360.0;
      if (x > xplus) x -= 360.0;

      xi = indx(x, xt, imx);
      yj = indx(y, yt, jmx);
      if (xi > 0.0 && yj > 0.0) {

        kw = 0;
        for (k = 0; k < kd; k++) {
          if (*(tw+k) != mtVal) {
            if (do_conv_pottmp) { //STEVE: made this conditional on option
              // get pressure as a function of depth in meters and gravitaional acceleration
              // use for input to theta (pot temp computation)
              // NOTE: zw must be depth in meters, not model levels! (issue with synthetic obs)
              p = press(*(zw+k),980.0);
              p0 = 0.0;
              // computes the potential temperature conversion, given pressure, temperature, salinity, and reference pressure
              *(t+k) = theta(p,*(tw+k),*(sw+k),p0);
            } else{
              *(t+k) = *(tw+k);
            }
            // assign depth in meters
            *(z+k) = *(zw+k);
          } else {
            // Apply a spurious value identifier
            *(t+k) = spv;
            *(z+k) = *(zw+k);
          }
        }

        // Compute the gradient of t at the model levels. Also perform some other computations
        // tz is the scaling to apply to the se error profiles
        // a alternates between the temperature and the temperature error
        if (do_comp_se) {
          cmpTz();
        }
        // tz is used only for scaling the standard deviation (error) by depth, so it's not needed
        // if the observational errors are read in from the netcdf obs file

        kd = kd  < kamx ? kd : kamx;

        if (kd > 1 && inbnds(xi,yj)) {

          for (k = 0; k < kd; k++) {
            *(a+2*k) = *(t+k);
            if (*(t+k) != spv) {
              // specify the standard deviation (obs error)
              if (do_comp_se) {
                se = se0 + seF * *(tz+k);
              } else {
                // NOTE: for the synthetic observations, this should be read in:
                se = *(ew+k);
              }
              *(a+2*k+1) = 1.0 / (se*se);
            } else {
              *(a+2*k+1) = 0.0;
            }
          }

          // write out the observation data to the scratch file fw
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
    // Close the temp and salt files
    closePrfnc('T');
    closePrfnc('S');
    fclose(fslts);
  }
  close(fw);

  // Open the scratch file
  printf("DEBUG :: Opening scrtch file...\n");
  fd = open("scrtch",O_RDONLY);
  // Create a new file to write the final observation output to:
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
/* DBG
  system("rm *tmp.nc");
 DBG */

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

void initPrfnc(V)
char V;
{
  int k, n, ncstat;
  char str[51];

  if (V == 'T') {
    ncstat = nc_open(prfTFile, NC_NOWRITE, &ncidt);
    if (ncstat != NC_NOERR) {
      printf("Cannot open %s\n", prfTFile);
      printf("%s\n",nc_strerror(ncstat));
      exit(-1);
    }

    ncstat = nc_inq_ndims(ncidt, &ndims);
    if (ncstat != NC_NOERR) {
      printf("Cannot get number of dimensions\n");
      printf("%s\n",nc_strerror(ncstat));
      exit(-1);
    }

    pcnt = 0;
    for (n = 0; n < ndims; n++) {
      ncstat = nc_inq_dim(ncidt, n, str, &len);
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

    ncstat = nc_inq_nvars(ncidt, &nvars);
    if (ncstat != NC_NOERR) {
      printf("Cannot get number of variables\n");
      printf("%s\n",nc_strerror(ncstat));
      exit(-1);
    }
    for (n = 0; n < nvars; n++) {
      ncstat = nc_inq_varname(ncidt, n, str);
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
        tvid = n;
        doT = 1;
        strcpy(pname,"temp");
      } else if (!strncmp(str,"stde",4)) {
        evid = n;
      } else if (!strncmp(str,"model",5)) {
        mvid = n;
      }
    }
    ncstat = nc_get_att_float(ncidt, tvid, "missing_value", &mtVal);
    if (ncstat != NC_NOERR) {
      printf("Cannot get missing_value for T\n");
      printf("%s\n",nc_strerror(ncstat));
      exit(-1);
    }

    ncstat = nc_get_att_text(ncidt, NC_GLOBAL, "Date", str);
    if (ncstat != NC_NOERR) {
      printf("Cannot get Global Date\n");
      printf("%s\n",nc_strerror(ncstat));
      exit(-1);
    }
    sscanf(str,"%d%*c%d%*c%d",&year,&month,&day);

    if (first) {
      *start = 0;
      *count = kmx;
      ncstat = nc_get_vara_float(ncidt, zvid, start, count, zw);
      if (ncstat != NC_NOERR) {
        printf("Cannot get Z variable\n");
        printf("%s\n",nc_strerror(ncstat));
        exit(-1);
      }
      first = 0;
    }
  } else {
    ncstat = nc_open(prfSFile, NC_NOWRITE, &ncids);
    if (ncstat != NC_NOERR) {
      printf("Cannot open %s\n", prfSFile);
      printf("%s\n",nc_strerror(ncstat));
      exit(-1);
    }

    pcnt = 0;
    for (n = 0; n < ndims; n++) {
      ncstat = nc_inq_dim(ncids, n, str, &len);
      if (ncstat != NC_NOERR) {
        printf("Cannot get dimensions\n");
        printf("%s\n",nc_strerror(ncstat));
        exit(-1);
      }
      if (!strncmp(str,"count",5)) {
        pcnt = len;
      }
    }

    ncstat = nc_inq_nvars(ncids, &nvars);
    if (ncstat != NC_NOERR) {
      printf("Cannot get number of variables\n");
      printf("%s\n",nc_strerror(ncstat));
      exit(-1);
    }
    for (n = 0; n < nvars; n++) {
      ncstat = nc_inq_varname(ncids, n, str);
      if (ncstat != NC_NOERR) {
        printf("Cannot get variable name\n");
        printf("%s\n",nc_strerror(ncstat));
        exit(-1);
      }
      if (!strncmp(str,"salt",4)) {
        svid = n;
        strcpy(pname,"salt");
      }
    }
    ncstat = nc_get_att_float(ncids, svid, "missing_value", &msVal);
    if (ncstat != NC_NOERR) {
      printf("Cannot get missing_value for S\n");
      printf("%s\n",nc_strerror(ncstat));
      exit(-1);
    }
  }
}

/* -------------------------------------------------------------- */

void getHeadnc(V, np)
char V;
int np;
{
  int ncstat, ncidlc;

  if (V == 'T') {
    ncidlc = ncidt;
  } else {
    ncidlc = ncids;
  }

  *start = np;
  *count = 1;
  ncstat = nc_get_vara_float(ncidlc, xvid, start, count, &x);
  if (ncstat != NC_NOERR) {
    printf("Cannot get variable in %s\n", "xlon");
    printf("%s\n",nc_strerror(ncstat));
    exit(0);
  }

  ncstat = nc_get_vara_float(ncidlc, yvid, start, count, &y);
  if (ncstat != NC_NOERR) {
    printf("Cannot get variable in %s\n", "ylat");
    printf("%s\n",nc_strerror(ncstat));
    exit(0);
  }

  ncstat = nc_get_vara_float(ncidlc, hvid, start, count, &xhour);
  if (ncstat != NC_NOERR) {
    printf("Cannot get variable in %s\n", "hour");
    printf("%s\n",nc_strerror(ncstat));
    exit(0);
  }
}

/* -------------------------------------------------------------- */

void getPrfnc(V, np)
char V;
int np;
{
  int ncstat, ncidlc;

  if (V == 'T') {
    *start = np;
    *count = 1;
    ncstat = nc_get_vara_float(ncidt, xvid, start, count, &x);
    if (ncstat != NC_NOERR) {
      printf("Cannot get variable in %s\n", "xlon");
      printf("%s\n",nc_strerror(ncstat));
      exit(0);
    }

    ncstat = nc_get_vara_float(ncidt, yvid, start, count, &y);
    if (ncstat != NC_NOERR) {
      printf("Cannot get variable in %s\n", "ylat");
      printf("%s\n",nc_strerror(ncstat));
      exit(0);
    }

    ncstat = nc_get_vara_float(ncidt, hvid, start, count, &xhour);
    if (ncstat != NC_NOERR) {
      printf("Cannot get variable in %s\n", "hour");
      printf("%s\n",nc_strerror(ncstat));
      exit(0);
    }

    *start = np;
    *count = 1;
    *(start+1) = 0;
    *(count+1) = nc2;
    ncstat = nc_get_vara_text(ncidt, pfvid, start, count, plat);
    if (ncstat != NC_NOERR) {
      printf("Cannot get variable in %s\n", "plat");
      printf("%s\n",nc_strerror(ncstat));
      exit(0);
    }

    *start = np;
    *count = 1;
    *(start+1) = 0;
    *(count+1) = nc1;
    ncstat = nc_get_vara_text(ncidt, ptvid, start, count, ptyp);
    if (ncstat != NC_NOERR) {
      printf("Cannot get variable in %s\n", "ptyp");
      printf("%s\n",nc_strerror(ncstat));
      exit(0);
    }

/*

    ncstat = nc_get_vara_text(ncidt, sdvid, start, count, sid);
    if (ncstat != NC_NOERR) {
      printf("Cannot get variable in %s\n", "sid");
      printf("%s\n",nc_strerror(ncstat));
      exit(0);
    }

    ncstat = nc_get_vara_text(ncidt, qkvid, start, count, &qkey);
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
    ncstat = nc_get_vara_float(ncidt, tvid, start, count, tw);
    if (ncstat != NC_NOERR) {
      printf("Cannot get variable in %s\n", pname);
      printf("%s\n",nc_strerror(ncstat));
      exit(0);
    }

    // adding the standard deviation from file:
    *start = np;
    *count = 1;
    *(start+1) = 0;
    *(count+1) = kmx;
    ncstat = nc_get_vara_float(ncidt, evid, start, count, ew);
    if (ncstat != NC_NOERR) {
      printf("Cannot get variable in stde\n");
      printf("%s\n",nc_strerror(ncstat));
      exit(0);
    }
    
  } else {
    *start = np;
    *count = 1;
    *(start+1) = 0;
    *(count+1) = kmx;
    ncstat = nc_get_vara_float(ncids, svid, start, count, sw);
    if (ncstat != NC_NOERR) {
      printf("Cannot get variable in %s\n", pname);
      printf("%s\n",nc_strerror(ncstat));
      exit(0);
    }
  }

}

/* -------------------------------------------------------------- */

void closePrfnc(V)
char V;
{
  int ncstat;

  if (V == 'T') {
    ncstat = nc_close(ncidt);
    if (ncstat != NC_NOERR) {
      printf("Cannot file %s\n", prfTFile);
      printf("%s\n",nc_strerror(ncstat));
      exit(0);
    }
  } else {
    ncstat = nc_close(ncids);
    if (ncstat != NC_NOERR) {
      printf("Cannot file %s\n", prfSFile);
      printf("%s\n",nc_strerror(ncstat));
      exit(0);
    }
  }
}

/* -------------------------------------------------------------- */
//
// Interpolate to fill in missing values of salinity where there are values of temperature
//
void Sfill()
{
  int k, k1, k2, kk, km, kp;
  float dz, rm, rp;

  k1 = 0;
  if (*sw == msVal) {
    for (k1 = 1; k1 < kd; k1++) {
      if (*(sw+k1) != msVal) {
        for (k = 0; k < k1; k++) {
          *(sw+k) = *(sw+k1);
        }
        break;
      }
    }
  }
  k2 = kd-1;
  if (*(sw+kd-1) == msVal) {
    for (k2 = kd-2; k2 >=0; k2--) {
      if (*(sw+k2) != msVal) {
        for (k = kd-1; k > k2; k--) {
          *(sw+k) = *(sw+k2);
        }
        break;
      }
    }
  }
  for (k = k1; k < k2; k++) {
    if (*(sw+k) == msVal) {
      km = k-1;
      for (kp = k+1; kp <= k2; kp++) {
        if (*(sw+kp) != msVal) {
          dz = *(zw+kp) - *(zw+km);
          for (kk = k; kk <= kp; kk++) {
            rm = (*(zw+kp) - *(zw+kk)) / dz;
            rp = (*(zw+kk) - *(zw+km)) / dz;
            *(sw+kk) = *(sw+km)*rm + *(sw+kp)*rp;
          }
        }
      }
    }
  }
}
     

/* -------------------------------------------------------------- */

void cmpTz()
{
  int k, kk, km, kp, kv2, cnt;
  float tzmn, tzmx;

  // Estimate temperature gradients
  // i.e. every level is linearly fit from below and above
  for (k = 1; k < kd-1; k++) {
    if (*(t+k-1) != spv && *(t+k+1) != spv) {
      // Estimate temperature gradient
      *(tz+k) = (*(t+k-1) - *(t+k+1)) / (*(z+k+1) - *(z+k-1));
      // enforce non-negative temperature gradient
      if (*(tz+k) < 0.0) *(tz+k) = 0.0;
    } else {
      *(tz+k) = spv;
    }
  }
  *tz = spv;
  *(tz+kd-1) = spv;

  // For all levels with obs data,
  for (k = 0; k < kd; k++) {
    // if one of the values was computed with a missing value, then
    // we have to go back and find the nearest non-missing value to do the computation
    if (*(tz+k) == spv && *(t+k) != spv) {
      // k minus
      km = -1;
      // Go back to the first non-missing value
      for (kk = k-1; kk >= 0; kk--) {
        if (*(tz+kk) != spv) {
          km = kk;
          break;
        }
      }
      // k plus
      kp = kd;
      // Go forward to the first non-missing value
      for (kk = k+1; kk < kd; kk++) {
        if (*(tz+kk) != spv) {
          kp = kk;
          break;
        }
      }

      // if k-minus and k-plus are in the range of observed levels,
      if (km >= 0 && kp < kd) {
        // assign the value to level k of the nearest neighbor (km or kp)
        if ((k-km) <= (kp-k)) {
          *(tz+k) = *(tz+km);
        } else {
          *(tz+k) = *(tz+kp);
        }
      // otherwise, handle the boundary cases
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
    // divide tz by dz
    for (k = 0; k < kd; k++) {
      if (*(tz+k) != spv) {
        *(tz+k) /= *(dz+k);
      }
    }
  }

  if (SqRt) {
    // take sqaure root of tz
    for (k = 0; k < kd; k++) {
      if (*(tz+k) != spv) {
        if (*(tz+k) >  0.0)
          *(tz+k) = sqrt(*(tz+k));
        else
          *(tz+k) = 0.0;
      }
    }
  }

  //STEVE: I think this averages consecutive tz's to smooth out the observations
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

/*  Adiabatic temperature gradient  Bryden (DSR, v.20, p.404, 1973)
       p   = pressure (dbar)
       t   = temperature (C)
       s   = salinity (ppt)
       atg = adiabatic temperature gradient (C/dbar)
*/

float atg(float p, float t, float s)
{
  float ds, a;

  ds = s-35.0;
  a = (((-2.1687e-16*t+1.8676e-14)*t-4.6206e-13)*p
         +((2.7759e-12*t-1.1351e-10)*ds+((-5.4481e-14*t
         +8.733e-12)*t-6.7795e-10)*t+1.8741e-8))*p
         +(-4.2393e-8*t+1.8932e-6)*ds
         +((6.6228e-10*t-6.836e-8)*t+8.5258e-6)*t+3.5803e-5;
  return(a);
}

/* -------------------------------------------------------------- */

/*  Supply depth (zz) in meters and grav acc'l (g) in cm/sec**2 */

float press(float z, float g)
{
  int itr = 0;
  float p, a0;
  double e, ae, es;

  p = z*(1.0076+z*(2.3487e-6-z*1.2887e-11));
  e = zeta(p,g)-z;
  ae = fabs(e);
  es = ae*2.;
  while (ae > 0.01 && ae < es && itr++ < 20) {
    a0 = 0.972643+p*(1.32696e-5-p*(6.228e-12+p*1.885e-16));
    a0 = a0/(1.0+1.83e-5*p);
    p = p-((g+1.113e-4*p)/a0)*e*0.001;
    es = ae;
    e = zeta(p,g)-z;
    ae = fabs(e);
  }
  if (ae > es)
    printf("   PRESS DIVERGENCE");
  return(p);
}

/* -------------------------------------------------------------- */

/*   calculates the in situ temperature of sea water given a potential
     temperature and its reference pressure by integrating the
     adiabatic lapse rate using a fourth-order Runge-Kutta algorithm.

     p      =   pressure (dbar)
     th     =   potential temperature (C)
     s      =   salinity (ppt)
     pref   =   reference pressure (dbar)
     tempa  =   in situ temperature (C)
*/

float tempa(float p, float th, float s, float pref)
{
  int n;
  double sqrt2 = 0.7071067811865475;
  double a = 0.0000018;
  float del_p, del_t1, del_t2, del_t3, del_t4, tp, t;

  del_p = pref-p;
  del_t1 = del_p*atg(p,th,s);
  tp = th-0.5*del_t1;
  del_t2 = del_p*atg((p+0.5*del_p),tp,s);
  tp = th-(-0.5+sqrt2)*del_t1+(1.0-sqrt2)*del_t2;
  del_t3 = del_p*atg((p+0.5*del_p),tp,s);
  tp = th+sqrt2*del_t2+(1.0+sqrt2)*del_t3;
  del_t4 = del_p*atg(pref,tp,s);
  t = (th-(del_t1+(1.0-sqrt2)*del_t2*2.0 +
                   (1.0+sqrt2)*del_t3*2.0+del_t4)/6.0)+a*p;
  return(t);
}

/* -------------------------------------------------------------- */

/*   calculates the potential temperature of sea water at specified
     reference pressure by integrating the adiabatic lapse rate using
     a fourth-order Runge-Kutta algorithm.

     p      =   pressure (dbar)
     t      =   temperature (C)
     s      =   salinity (ppt)
     pref   =   reference pressure (dbar)
     theta  =   potential temperature (C)
*/

float theta(float p, float t, float s, float pref)
{
  double sqrt2 = 0.7071067811865475;
  float del_p, del_t1, del_t2, del_t3, del_t4, tp, th;

  del_p = pref-p;
  del_t1 = del_p*atg(p,t,s);
  tp = t+0.5*del_t1;
  del_t2 = del_p*atg((p+0.5*del_p),tp,s);
  tp = t+(-0.5+sqrt2)*del_t1+(1.0-sqrt2)*del_t2;
  del_t3 = del_p*atg((p+0.5*del_p),tp,s);
  tp = t-sqrt2*del_t2+(1.0+sqrt2)*del_t3;
  del_t4 = del_p*atg(pref,tp,s);
  th = (t+(del_t1+(1.0-sqrt2)*del_t2*2.0 +
                   (1.0+sqrt2)*del_t3*2.0+del_t4)/6.0);
  return(th);
}

/* -------------------------------------------------------------- */

/*
  Depth from pressure and gravity but ignoring dynamic height anomaly
  ref. Saunders and Fofonoff, DSR v. 23  1976
    p is pressure in dbar and glat is gravity
  rcm   Aug 1978
*/

float zeta(float p, float glat)
{
  float z;

  z = ((-3.434e-12*p+1.113e-7)*p+0.712953)*p+14190.7*log(1.0+1.83e-5*p);
  z = (z/(glat+1.113e-4*p))*1000.;

  return(z);
}

/* -------------------------------------------------------------- */

void usage()
{
  printf("Usage:\n");
  printf(" %s -f tFlst sFlst -g gridFile [options]\n", prog);
  printf("   -f tFlst sFlst - temperature and salinity file lists\n");
  printf("   -g gridFile  - netCDF grid_spec file for MOM4 model\n");
  printf("   -k LevMx     - maximum level for data assimilation\n");
  printf("   -r se0 seF   - set min/max range\n");
  printf("   -o outFile   - output file \n");

  printf(" Defaults:\n");
  printf("  SdZ:  %d\n", SdZ);
  printf("  SqRt: %d\n", SqRt);
  printf("  KAV:  %d\n", KAV);
  printf(" Controling min/max:\n");
  printf("  SE0:  %4.2f\n", SE0);
  printf("  SEF:  %4.2f\n", SEF);
}
