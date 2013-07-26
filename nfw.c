/*
NFW-fit a Gadget-File

gcc -lm -o ../bin/altix/nfw nfw.c libgad.o lmfit-2.2/lmmin.o lmfit-2.2/lm_eval.o
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "libgad.h"

#define USE 63

void usage()
{
  	  fprintf(stderr," v0.01\n");
  	  fprintf(stderr," -i <input file name>\n");
  	  fprintf(stderr," -use <bitcode particle types to use (default 2¹)>\n\n");
	  exit(1);
}


int main (int argc, char *argv[])
{
  FILE *fp;
  char infile[256];
  int i,j,k, usepart;
  struct gadpart *part;
  gadpart_dist *wpart;
  struct header head;
  fltarr cm={0,0,0};
  double rv=0, sfl=0;
  int calccenter=1;
  
  i=1;
  usepart=USE;
  if (1==argc) usage();
  while (i<argc)
    {
      if (!strcmp(argv[i],"-i"))
	{
	  i++;
	  strcpy(infile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-cm")) {
	cm[0]=atof(argv[++i]);
	cm[1]=atof(argv[++i]);
	cm[2]=atof(argv[++i]);
	i++;
      }
      else if (!strcmp(argv[i],"-nocenter")) {
	calccenter=0;
	i++;
      }
      else if (!strcmp(argv[i],"-rv")) {
	rv=atof(argv[++i]);
	i++;
      }
      else if (!strcmp(argv[i],"-soft")) {
	sfl=atof(argv[++i]);
	i++;
      }
      else if (!strcmp(argv[i],"-use")) {
	i++;
	if (!strcmp(argv[i],"all")) usepart=63;
	else usepart=atoi(argv[i]);
	i++;
      } 
      else if (*argv[i]!='-')
	{
	  strcpy(infile,argv[i]);
	  i++;
	} else {
	usage();
      }
    }

  unsigned int numpart_all;
  if (!(numpart_all=readgadget_part(infile, &head, &part))) 
    {
      extern int libgaderr;
      printf("error reading file %s\nError Code %d\n",infile, libgaderr);
      exit(1);
    }

  unsigned int numpart=0;
  for (i=0; i<6; i++)
    {
      if (usepart&(1<<i)) {numpart+=head.npart[i];}
    }
  wpart=(struct gadpart_dist*) malloc (numpart*sizeof(struct gadpart_dist));
  j=0;
  for (i=0; i<numpart_all;i++)
    {
	if ((1<<(part[i].type))&usepart)
	{
	  cpygadpart(&(wpart[j].part),&(part[i]));
	  wpart[j].dist=distance(wpart[j].part.pos, cm);
	  //	  printf("%f  | ", wpart[j].dist);
	  j++;
	}
    }
  qsort(wpart, j, sizeof(gadpart_dist), cmp_dist);
  //  pcenter(wpart, j, rv, cm);

  /*********************************************************************

      Program code goes here

  *********************************************************************/

  double par[2]={0.005,20};
  double c;
  double mvir=0;
  double rcs;
  int vcnt;
  if (calccenter)
    {
      const double maxdist =2000.0;
      pcenter(wpart, numpart, maxdist, cm, 2);
    }
  if (rv==0) rv= r200(wpart, j, 200, head, &vcnt, &mvir);
  c = nfwfit(par, wpart, j, rv ,sfl, &rcs);
  printf("center %f %f %f\n", cm[0], cm[1], cm[2]);
  printf("Rvir %g | Vcnt %d | Mvir %g\n", rv, vcnt, mvir);
  printf("rs %g dc %g\n", par[1], par[0]);
  printf("reduced chi^2 %f\n", rcs);
  printf("c %f\n", c);
  return 0;
}
