#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

typedef float fltarr[3];

int main  (int argc, char *argv[])
{
  int i;
  int num;
  char filename[256];
  FILE *fp;
  fltarr *pos;
  int verbose = 0;
  int shortheader = 0;
  int quiet = 0;
  i = 1;
  //  if (argc==1) usage();
  while (i<argc)
    {
      if (*argv[i]!='-')
        {
          strcpy(filename,argv[i]);
          i++;
        }
      else if (!strcmp(argv[i],"-v"))
	{
	  i++;
	  verbose = 1;
	} 
      else if (!strcmp(argv[i],"-q"))
        {
          i++;
          quiet = 1;
        }

    }

  fp = fopen(filename, "rb");
  fread(&num, sizeof(int), 1, fp);
  pos = calloc(num, sizeof(fltarr));
  fread(&pos[0], sizeof(fltarr), num, fp);
  fclose(fp);
  if (verbose)
    {
      printf("numIDs  | %d\n", num);
    }
  else
    for ( i=0; i<num; i++) 
      {
	printf("%f %f %f\n", pos[i][0], pos[i][1], pos[i][2]);
      }

  return 0;
}
