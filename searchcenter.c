/*
Program to do something with gadget-files

gcc -std=c99 -lm -o a.out whatever.c libgad.o
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "libgad.h"
#include "OctTree.h"

#define USE 2

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
  char infile[256], outfile[256];
  int i,j,k, usepart;
  int verbose = 0;
  struct gadpart *part, *wpart;
  struct header head;
  double conv_dist = 1.;
  int exclude = 0;
  int write = 0;
  float ratio = .5;
  int maxdepth = 10;
  fltarr center = {36000, 36000, 36000};
  float srad = 15000;
  
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
      else if (*argv[i]!='-')
	{
	  strcpy(infile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-o"))
	{
	  i++;
	  strcpy(outfile,argv[i]);
	  i++;
	  write=1;
	}
      else if (!strcmp(argv[i],"-cm"))
	{
	  i++;
	  center[0]=atof(argv[i++]);
	  center[1]=atof(argv[i++]);
	  center[2]=atof(argv[i++]);
	}
      else if (!strcmp(argv[i],"-v"))
	{
	  i++;
	  verbose = 1;
	}
      else if (!strcmp(argv[i],"-cd"))
	{
	  i++;
	  conv_dist = atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-r"))
	{
	  i++;
	  ratio = atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-sr"))
	{
	  i++;
	  srad = atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-x"))
	{
	  i++;
	  exclude = atoi(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-d"))
	{
	  i++;
	  maxdepth = atoi(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-use")) {
	i++;
	if (!strcmp(argv[i],"all")) usepart=63;
	else usepart=atoi(argv[i]);
	i++;
      }
       else {
	usage();
      }
    }

  unsigned int numpart_all;
  extern int basic;
  basic = 1;
  if (!(numpart_all=readgadget_part(infile, &head, &part))) 
    {
      extern int libgaderr;
      printf("error reading file %s\nError Code %d\n",infile, libgaderr);
      exit(1);
    }
  extern float BOXSIZE;
  BOXSIZE = head.boxsize;

  if (verbose)
    {
      printf("building OctTree\n");
      fflush(stdout);
    }
  OctNode *onode;
  set_periodic_boundaries(head);
  extern double crit_dens;
  crit_dens *= pow(conv_dist,3);
  buildTreeBox(&onode, part, head, 200, maxdepth);
  if (exclude)
    {
      rejectNodes(onode, usepart, exclude, ratio);
    }

  /*********************************************************************

      Program code goes here

  *********************************************************************/
  
  gadpart **part_pnt;
  unsigned int npart=0;
  unsigned int bsize=0;

  findParticles( onode, center, srad, &part_pnt, &npart, &bsize);
  printf("npart: %d\n", npart);

  for (i=0; i<3; i++) center[i] = 0;
  gadpart_dist *dpart = (gadpart_dist*) malloc(npart* sizeof(gadpart_dist));
  for (i=0; i<npart; i++)
    {
      dpart[i].part = *part_pnt[i];
    }
  pcenter(dpart, npart, srad, center, usepart);
  free(dpart);
  printf("%f %f %f\n", center[0], center[1], center[2]);
  fp = fopen("cm.dat", "w");
  fprintf(fp,"%f %f %f\n", center[0], center[1], center[2]);
  fclose(fp);
  if (write)
    {
      gadpart *wpart = (gadpart*) malloc(npart* sizeof(gadpart));
      for ( i=0; i<npart; i++)
	wpart[i] = *(part_pnt[i]);
      struct header out = cphead(head, wpart, npart);
      writegadget_part(outfile, out, wpart);
      free(wpart);
    }
  return 0;
}
