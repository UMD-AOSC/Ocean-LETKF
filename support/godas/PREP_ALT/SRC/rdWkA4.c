/*  rdWkA4 reads and lists a sample from a weeky altimeter file  */

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

void  usage();

int npts, nstrt, nend;
float buf[4];
int ibuf[5];
int year, month, day, hour;
float xo, yo, a, evr;
char wklyFile[51], prog[51];

main(argc,argv)
int argc;
char *argv[];
{
  int n, nb, fd;

  strcpy(prog,argv[0]);
  strcpy(wklyFile,"EMPTY");
  nstrt = 0;
  nend = 19;
  n = 1;
  while (n < argc) {
    if (!strcmp(argv[n],"-f")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -f requires a file name as a parameter.\n");
        usage();
        exit(0);
      }
      strcpy(wklyFile, argv[n]);
    } else if (!strcmp(argv[n],"-s")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -s requires an integer as a parameter.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%d", &nstrt);
      nstrt--;
    } else if (!strcmp(argv[n],"-e")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -e requires an integer as a parameter.\n");
        usage();
        exit(0);
      }
      sscanf(argv[n],"%d", &nend);
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
  if (!strcmp(wklyFile,"EMPTY")) {
    printf("Error: An input file must be given.\n");
    usage();
    exit(0);
  }

  fd = open(wklyFile, O_RDONLY);

  read(fd, &nb, 4);
  read(fd, &npts, 4);
  read(fd, &nb, 4);
  for (n = 0; n < npts; n++) {
    if (n >= nend) {
      break;
    }
    read(fd, &nb, 4);
    read(fd, ibuf, 20);
    read(fd, &nb, 4);
    read(fd, &nb, 4);
    read(fd, &buf, 16);
    read(fd, &nb, 4);
    if (n >= nstrt) {
      year = *(ibuf+0);
      month = *(ibuf+1);
      day = *(ibuf+2);
      hour = *(ibuf+3);
      xo = *(buf+0);
      yo = *(buf+1);
      a = *(buf+2);
      evr = *(buf+3);
      printf("%5d %4d-%2.2d-%2.2d : %2.2d %8.2f %8.2f %8.2f    %g\n",n,year,month,day,hour,xo,yo,a,evr);
    }
  }
  close (fd);
}

/* -------------------------------------------------------------- */

void usage()
{
  printf("Usage:\n");
  printf(" %s -f wklyFile [options]\n", prog);
  printf("   -f wklyFile - file of assm data\n");
  printf("   -s nstrt    - first record to print\n");
  printf("   -e nend     - last record to print\n");
}
