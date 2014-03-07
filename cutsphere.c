/*
Program to cut a sphere out of a gadget snapshot

gcc -std=c99 -lm -o a.out whatever.c libgad.o
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
          fprintf(stderr," -o <output file name>\n");          
          fprintf(stderr," -r <radius of the sphere>\n");
          fprintf(stderr," -cm <center of the sphere>\n");
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
  fltarr cm;
  double rad;
  int maxpart = 0;
  int basic_only = 0;

  
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
      else if (!strcmp(argv[i],"-cm"))
        {
          i++;
          cm[0]=atof(argv[i++]);
          cm[1]=atof(argv[i++]);
          cm[2]=atof(argv[i++]);
        }
      else if (!strcmp(argv[i],"-r"))
        {
          i++;
          rad = atof(argv[i++]);
        }  
      else if (!strcmp(argv[i],"-use")) {
	i++;
	if (!strcmp(argv[i],"all")) usepart=63;
	else usepart=atoi(argv[i]);
	i++;
      }
      else if (!strcmp(argv[i],"-b")) {
        basic_only=1;
	i++;
      }
      else if (!strcmp(argv[i],"-max")) 
      {
	i++;
	maxpart=atof(argv[i++]);
      } else {
	usage();
      }
    }

  unsigned int numpart_all;

  if (basic_only)
    {
      extern int basic;
      basic = 1;
    }
  

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
      struct header outhead;
      outhead = head;
      for ( j = 0; j < 6; j++)
      {
        outhead.npart[j] = 0;
        outhead.nall[j] = 0;
      }
      j = 0;
      if (maxpart==0)
        maxpart = numpart_all;
      wpart = (gadpart*) malloc (maxpart * sizeof(gadpart));
      for (i = 0; i < numpart_all; i++)
      {
        if (!(( 1 << part[i].type) & usepart))
                continue;
        if (distance(cm, part[i].pos) > rad)
                continue;

        wpart[j++] = part[i];
        outhead.npart[part[i].type]++;
        outhead.nall[part[i].type]++;
        if ( j>= maxpart)
                break;
      }

      writegadget_part(outfile, outhead, wpart);
      free(wpart);
  return 0;
}
