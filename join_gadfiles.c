/*
Program to do join (numfiles>1) the snapshots of a gadget simulation

gcc -std=c99 -lm -o a.out whatever.c libgad.o
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "libgad.h"

#define USE 2

void usage()
{
  	  fprintf(stderr," v0.01\n");
  	  fprintf(stderr," -i <input file name>\n");
  	  fprintf(stderr," -o <outut file name>\n");
  	  fprintf(stderr," -use <bitcode particle types to use (default 2¹)>\n\n");
	  exit(1);
}


int main (int argc, char *argv[])
{
  FILE *fp;
  char infile[256];
  char outfile[256];
  int i,j,k, usepart;
  struct gadpart *part, *wpart;
  struct header head;
  double conv_dist=1.0;
  double conv_mass=1.0;
  int outputset=0;
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
      else if (!strcmp(argv[i],"-o"))
	{
	  i++;
	  strcpy(outfile,argv[i]);
	  i++;
	  outputset=1;
	}       
      else if (*argv[i]!='-')
	{
	  strcpy(infile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-dc"))
	{
	  i++;
	  conv_dist=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-mc"))
	{
	  i++;
	  conv_mass=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-use")) {
	i++;
	if (!strcmp(argv[i],"all")) usepart=63;
	else usepart=atoi(argv[i]);
	i++;
      } else {
	usage();
      }
    }

  if (!outputset) usage();

  unsigned int numpart_all;

  if (!(numpart_all=readgadget_part(infile, &head, &part))) 
    {
      extern int libgaderr;
      printf("error reading file %s\nError Code %d\n",infile, libgaderr);
      exit(1);
    }
  convertunits(&head, part, conv_mass, conv_dist);
  head.numfiles=1;
  /*********************************************************************

      Program code goes here

  *********************************************************************/

  writegadget_part(outfile, head, part);
  return 0;
}
