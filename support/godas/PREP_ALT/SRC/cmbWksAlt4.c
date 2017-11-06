/*  cmbWksAlt4 combines files produced by spltWkAlt4 into a multiweek
    file to be used by the assimilating model.  See also spltWkPrf4 */

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

#define MMX 20

void usage();

char listFile[81], wFiles[MMX][81], outFile[81], prog[51];
int nb, n20b = 20, n16b = 16, n4b = 4, nfls = 0;
float buf[4];
int ibuf[5];

mode_t mode = 0644;

main(argc,argv)
int argc;
char *argv[];
{
  FILE *fs;
  char str[81], cll[8];
  int fin, fout, n, k;
  int nobs;

  strcpy(prog,argv[0]);
  strcpy(listFile,"EMPTY");
  strcpy(outFile,"EMPTY");
  n = 1;
  while (n < argc) {
    if (!strcmp(argv[n],"-f")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -f requires a file name as a parameter.\n");
        usage();
        exit(99);
      }
      strcpy(listFile, argv[n]);
    } else if (!strcmp(argv[n],"-o")) {
      n++;
      if (n >= argc || !strncmp(argv[n],"-",1)) {
        printf("Error: the option -o requires a file name as a parameter.\n");
        usage();
        exit(99);
      }
      strcpy(outFile, argv[n]);
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
  if (!strcmp(listFile,"EMPTY")) {
    printf("Error: An input list file must be given.\n");
    usage();
    exit(99);
  }
  if (!strcmp(outFile,"EMPTY")) {
    printf("Error: An output file must be given.\n");
    usage();
    exit(99);
  }

  fs = fopen(listFile, "r");
  while (nfls < MMX && fgets(str,80,fs)) {
    sscanf(str, "%s", wFiles[nfls]);
    nfls++;
  }
  fclose(fs);


  n = 0;
  fout = creat(outFile, mode);
  while (n < nfls) {
    fin = open(wFiles[n], O_RDONLY);
    printf("   %s", wFiles[n]);
    read(fin, &nb, 4);
    read(fin, &nobs, 4);
    read(fin, &nb, 4);
    printf("   nobs: %d\n", nobs);
    write(fout, &n4b, 4);
    write(fout, &nobs, 4);
    write(fout, &n4b, 4);
    for (k = 0; k < nobs; k++) {
      read(fin, &nb, 4);
      read(fin, ibuf, 20);
      read(fin, &nb, 4);
      read(fin, &nb, 4);
      read(fin, buf, 16);
      read(fin, &nb, 4);
      write(fout, &n20b, 4);
      write(fout, ibuf, 20);
      write(fout, &n20b, 4);
      write(fout, &n16b, 4);
      write(fout, buf, 16);
      write(fout, &n16b, 4);
    }
    close(fin);
    n++;
  }
  close(fout);
}

/* ================================================================= */

void usage()
{
  printf("Usage:\n");
  printf(" %s -f listFile -o outFile\n", prog);
  printf("   -f listFile - list of files to be combined\n");
  printf("   -o outFile  - output file name\n");
}
