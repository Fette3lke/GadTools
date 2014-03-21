/*

gcc -fopenmp -lm  -lgsl -lgslcblas -lgad -L ./ mk_id_list.c -o ~/bin/mk_id_list OctTree.o

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "libgad.h"
#include "ompfuncs.h"

#define USE 63

// #ifdef LONGIDS
// typedef unsigned long long IDtype;
// #else 
// typedef unsigned int IDtype;
// #endif

int cmp_IDtype (const void *first, const void *second)
{
  IDtype *a = (IDtype  *)first;
  IDtype *b = (IDtype  *)second;
  if (*a > *b) return 1;
  else if (*a < *b) return -1;
  else return 0;
}


const int MAX_HALO_ID = 100000;
const float EXTEND = 500;
const float TRACE_FACTOR = 2.;
const float SEARCHDIST = 25;
const float MAXDIST = 3000.;
const float SOFTENING = 1.0;

void usage()
{
  fprintf(stderr," search positions of ID list - reads a list of IDs and creates position file of corresponding particles\n");
  fprintf(stderr,"\t-o  \t<ID list base file name>\n");
  fprintf(stderr,"\t-i  \t<snaphsot file name>\n");
  fprintf(stderr,"\t-max\t<max Halo ID>\n");
  fprintf(stderr,"\t-use\t<bitcode particle types to use (default 2¹)>\n\n");
  exit(1);
}


int main (int argc, char *argv[])
{
  FILE *fp;
  char infile[256];
  char outbase[256];
  char catname[256];
  char **output;
  int i,j,k, usepart;
  struct gadpart *part, *wpart;
  struct header head;
  int max_halo_id = MAX_HALO_ID;
  float extend = EXTEND;
  float trace_factor = TRACE_FACTOR;
  int verbose = 0;
  float searchdist = SEARCHDIST;
  float def_maxdist = MAXDIST;
  double conv_dist = 1.;
  int start_id = 0;
  int num_halos = 0;
  int write_catalogue = 0;
  int write_gad_file = 0;
  int outpos = 1;
  double soft = SOFTENING;

  strcpy(outbase,"idlist");
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
	  strcpy(outbase,argv[i]);
	  i++;
	} 
      else if (!strcmp(argv[i],"-c"))
	{
	  i++;
	  strcpy(catname,argv[i]);
	  write_catalogue = 1;
	  i++;
	}
      else if (!strcmp(argv[i],"-gad"))
        {
          i++;
          write_gad_file = 1;
        } 
      else if (!strcmp(argv[i],"-pos"))
        {
          i++;
          outpos = 1;
        } 
      else if (!strcmp(argv[i],"-v"))
	{
	  i++;
	  verbose = 1;
	} 
      else if (!strcmp(argv[i],"-s"))
	{
	  i++;
	  start_id = atoi(argv[i]);
	  i++;
	} 
      else if (!strcmp(argv[i],"-max"))
	{
	  i++;
	  max_halo_id = atoi(argv[i]);
	  i++;
	}       
      else if (!strcmp(argv[i],"-e"))
	{
	  i++;
	  extend = atof(argv[i]);
	  i++;
	}  
      else if (!strcmp(argv[i],"-sfl"))
	{
	  i++;
	  soft = atof(argv[i]);
	  i++;
	}       
      else if (!strcmp(argv[i],"-md"))
	{
	  i++;
	  def_maxdist = atof(argv[i]);
	  i++;
	}       
      else if (!strcmp(argv[i],"-cd"))
	{
	  i++;
	  conv_dist = atof(argv[i]);
	  i++;
	}  
      else if (!strcmp(argv[i],"-tf"))
	{
	  i++;
	  trace_factor = atof(argv[i]);
	  i++;
	}       
      else if (!strcmp(argv[i],"-sd"))
	{
	  i++;
	  searchdist = atof(argv[i]);
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


  if (verbose)
  {
    printf("reading snapshot\n");
    fflush(stdout);
  }

  unsigned int numpart_all;

  if (!(numpart_all=readgadget_part(infile, &head, &part))) 
    {
      extern int libgaderr;
      printf("error reading file %s\nError Code %d\n",infile, libgaderr);
      exit(1);
    }
  if (verbose)
  {
    printf("sorting snapshot\n");
    fflush(stdout);
  }
  myqsort(part, numpart_all, sizeof(gadpart), cmp_id);
  /*********************************************************************

      Program code goes here

  *********************************************************************/
  if (verbose)
  {
    printf("main loop...\n");
    fflush(stdout);
  }

  int haloid;
#pragma omp parallel for private (i,j,k) reduction (+ : num_halos)
  for ( haloid = start_id; haloid <= max_halo_id; haloid++ ) 
    {
      char idlistname[128];
      sprintf(idlistname, "%s_%d", outbase, haloid);
      char posfilename[128];
      sprintf(posfilename, "%s_positions_%d", outbase, haloid);
      FILE *fp = fopen(idlistname, "rb");
      if (fp == NULL) 
        {
          continue;
        }     
      num_halos++;

      int numids=0;
      fltarr center;
      IDtype *idlist;
      fltarr *pos = NULL;
      float maxdist = 0;
      fread(&numids, sizeof(int), 1, fp);
      fread(center, sizeof(float), 3, fp);
      if (verbose)
      {
	printf("haloid %d | numids %d | center %g %g %g\n", haloid, numids, center[0], center[1], center[2]);
	fflush(stdout);
      }
      if (numids)
	{
	  fread(&maxdist, sizeof(float), 1, fp);
	  idlist = calloc(numids, sizeof(IDtype));
	  fread(&idlist[0], sizeof(IDtype), numids, fp);
	  qsort(idlist, numids, sizeof(IDtype), cmp_IDtype);
	  pos = (fltarr*) calloc(numids, sizeof(fltarr));
	}
      fclose(fp);
      if (verbose)
	{
	  printf("haloid %d | center %g %g %g\n", haloid, center[0], center[1], center[2]);
	   fflush(stdout);  
	}

      if (numids)
	{
	  int numfnd = 0;
	  gadpart *start = part;
	  for ( i = 0; i < numids; i++ )
	    {
	      gadpart *fnd;
	      gadpart idpart;
	      idpart.id = idlist[i];
	      long int size = &part[numpart_all] - start;

	      fnd = bsearch( &idpart, start, size, sizeof(gadpart), cmp_id);
	      if (fnd != NULL)
		{
		  start = fnd;
		  for ( j = 0; j < 3; j++)
		    pos[numfnd][j] = fnd->pos[j] / head.boxsize;		
		  numfnd++;
		}
	      //	      if (numfnd >= numids) break;
	    }
	  if (verbose)
	    {
	      printf("haloid %d | numfnd %d\n", haloid, numfnd);
	      if (numfnd != numids)
		{
		  fprintf(stderr, "particle not found | halo %d\n", haloid);
		  exit(1);
		}
	      fflush(stdout);
	    }

	}


      int totnumids = numids;
      fp = fopen(posfilename, "w");
      fwrite(&totnumids, sizeof(int), 1, fp);
      if (numids)
	fwrite(&pos[0], sizeof(fltarr), numids, fp);      
      fclose(fp);

      if (numids)
	{
	  if (outpos)
	    free(pos);
	}
    }

  return 0;
}
