#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#ifdef LONGIDS
typedef unsigned long long IDtype;
#else 
typedef unsigned int IDtype;
#endif


int main  (int argc, char *argv[])
{
  int i;
  int num;
  IDtype *id;
  char filename[256];
  FILE *fp;
  float pos[3];
  float maxdist;
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
      else if (!strcmp(argv[i],"-s"))
        {
          i++;
          shortheader = 1;
        }
      else if (!strcmp(argv[i],"-q"))
        {
          i++;
          quiet = 1;
        }

    }

  if (!quiet)
  {
#ifdef LONGIDS
        fprintf(stderr, "compiled for LONGIDS (8-Bytes)\n");
#else 
        fprintf(stderr, "compiled for 4-Byte IDs\n");
#endif
  }  

  fp = fopen(filename, "rb");
  fread(&num, sizeof(int), 1, fp);
  if (!shortheader)
        { 
          fread(pos, sizeof(float), 3, fp);
          fread(&maxdist, sizeof(float), 1, fp);
        }       
  //  fread(&id[0], sizeof(IDtype), num, fp);
  id = calloc(num, sizeof(IDtype));
  //  printf("%d %d\n", num, sizeof(unsigned long));
  fread(&id[0], sizeof(IDtype), num, fp);
  fclose(fp);
  if (verbose)
    {
      printf("numIDs  | %d\n", num);
      if (!shortheader)
        {       
          printf("center  | %g %g %g\n", pos[0], pos[1], pos[2]);
          printf("maxdist | %f\n", maxdist);
        }    
    }
  else
    for ( i=0; i<num; i++) 
      {
	printf("%lu\n", id[i]);
      }

  return 0;
}
