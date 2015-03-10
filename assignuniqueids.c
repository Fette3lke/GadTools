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

#define USE 2

void usage()
{
  	  fprintf(stderr," v0.01\n");
  	  fprintf(stderr," -i <input file name>\n");
  	  fprintf(stderr," -i <output file name>\n");
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

  sprintf(outfile, "gad.out");
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
	}       
      else if (*argv[i]!='-')
	{
	  strcpy(infile,argv[i]);
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

  if (!(numpart_all=readgadget_part(infile, &head, &part))) 
    {
      extern int libgaderr;
      printf("error reading file %s\nError Code %d\n",infile, libgaderr);
      exit(1);
    }
  extern float BOXSIZE;
  BOXSIZE = head.boxsize;

  /*********************************************************************

      Program code goes here

  *********************************************************************/

  for ( i = 0; i < numpart_all; i++ )
    {
      part[i].id = i;
    }

  writegadget_part(outfile, head, part);

  return 0;
}
