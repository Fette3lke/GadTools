/*

gcc -fopenmp -lm  -lgsl -lgslcblas -lgad -L ./ mk_id_list.c -o ~/bin/mk_id_list OctTree.o

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "libgad.h"
#include "OctTree.h"

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
  fprintf(stderr," Make ID list - reads a list of IDs, finds center of corresponding particles in snapshot, adds particles to list that are inside trace_factor * R200 and overwrites ID list\n");
  fprintf(stderr,"\t-o  \t<ID list base file name>\n");
  fprintf(stderr,"\t-i  \t<snaphsot file name>\n");
  fprintf(stderr,"\t-tf \t<trace_factor>\n");
  fprintf(stderr,"\t-max\t<max Halo ID>\n");
  fprintf(stderr,"\t-sd \t<search distance>\n");
  fprintf(stderr,"\t-md \t<max distance from cm>\n");
  fprintf(stderr,"\t-cd \t<convert distances (factor needed to get to kpc, for calculating the critical density>\n");
  fprintf(stderr,"\t-c  \t<write basic properties to catalogue file>\n");
  fprintf(stderr,"\t-pos\t<write position file (to use with MUSIC)>\n");
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
  int outpos = 0;
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

  output = (char**) malloc (sizeof(char*) * (max_halo_id+1));
  for (i = 0; i < (max_halo_id+1); i++ )
        output[i] = (char*) malloc(sizeof(char) * 256);

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
  buildTreeBox(&onode, part, head, 200, 15);

  if (verbose)
    printf("numpart %d \nchecktree %d\n", numpart_all, checkOctTree(onode));

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
          sprintf(output[haloid], "");
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
	  if (outpos)
	    pos = (fltarr*) calloc(numids, sizeof(fltarr));
	}
      if (maxdist == 0) maxdist = def_maxdist;

      fclose(fp);
      gadpart **part_pnt;
      unsigned int npart=0;
      unsigned int bsize=0;
      findParticles( onode, center, (maxdist + extend), &part_pnt, &npart, &bsize);
      if (verbose)
	{
	  printf("haloid %d | center %g %g %g | npart %d | maxdist + extend %g\n", haloid, center[0], center[1], center[2], npart, maxdist+extend);
//	  float mindist = distance(part_pnt[0]->pos, center);
//	  for ( i = 0; i < npart; i++ )
//	    {
//	      float dum = distance(part_pnt[i]->pos, center);
//	      if ( dum < mindist ) mindist = dum;
//	    }
//	  printf("mindist %g\n", mindist);
	   fflush(stdout);  
	}
//      if (verbose)
//	{
//	  for ( i = 0; i < 200000; i+=10000 )
//	    printf("%lu\n", part_pnt[i]->id);
//	}
      if (numids)
	{
	  gadpart_dist *wpart = (gadpart_dist *) malloc (sizeof(gadpart_dist) * numids);	 
        if (wpart==NULL)
        {
                fprintf(stderr, "unable to allocate memory\n");
                exit(1);
        } 
	  int numfnd = 0;
	  for ( i = 0; i < npart; i++ )
	    {
	      IDtype *fnd;
	      IDtype id_key = part_pnt[i] -> id;	      

	      fnd = bsearch( &id_key, idlist, numids, sizeof(IDtype), cmp_IDtype);
	      if (fnd != NULL)
		{
		  if (outpos)
		    {
		      for ( j = 0; j < 3; j++)
			pos[numfnd][j] = part_pnt[i]->pos[j] / head.boxsize;
		    }
		  wpart[numfnd++].part = *part_pnt[i];
		}
	      if (numfnd >= numids) break;
	    }
	  if (verbose)
                {
	       printf("haloid %d | numfnd %d\n", haloid, numfnd);
                fflush(stdout);
                }

	  pcenter(wpart, numfnd, maxdist, center, usepart);
	  free(wpart);
	}

      gadpart_dist *wpart = (gadpart_dist *) malloc (sizeof(gadpart_dist) * npart);
      gadpart *outpart;
      if (write_gad_file)
        outpart = (gadpart *) malloc (sizeof(gadpart) * npart);
      for ( i = 0; i < npart; i++ )
	{
	  wpart[i].part = *part_pnt[i];
          if (write_gad_file)
	         outpart[i] = *part_pnt[i];
	}
      free(part_pnt);
      pcenter(wpart, npart, searchdist, center, usepart);
      if (verbose)
	{
	  printf("haloid %d | center %g %g %g\n", haloid, center[0], center[1], center[2]);
          if (write_gad_file)
                {
        	  struct header outhead = head;
        	  outhead.npart[1]=npart;
        	  outhead.nall[1]=npart;
                  char gadfilename[128];
                  sprintf(gadfilename,"halo_%d.dat", haloid);
        	  writegadget_part(gadfilename, outhead, outpart);
                }
	}
      if (write_gad_file)
            free(outpart);
      int vcnt = 0;
      double mvir = 0;
      qsort(wpart, npart, sizeof(gadpart_dist), cmp_dist);
      double rvir = r200(wpart, npart, 200, head, &vcnt, &mvir);      
      if (verbose)
	{
	  printf("haloid %d | vcnt %d | rvir %g\n", haloid, vcnt, rvir);
	  printf("haloid %d | center %g %g %g\n", haloid, center[0], center[1], center[2]);
	}
      if ((write_catalogue) && (rvir > 0))
      {
      		double nfw_c;
      		double par[2];
      		par[0]=0.005;
                par[1]= 20. / conv_dist;
                double rcs;
      		nfw_c = nfwfit(par, wpart, npart, rvir, soft, &rcs);
      		double xoff = xoffset(wpart, npart, rvir, center);
		sprintf(output[haloid], "%8d\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\n", haloid, mvir, rvir, nfw_c, rcs, xoff, center[0], center[1], center[2]);
      }
      else
      {
                rvir=0;
                sprintf(output[haloid], "%8d\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\n", haloid, rvir, rvir, rvir, rvir, rvir, center[0], center[1], center[2]);
      }
      double dist = 0;
      i = 0;
      int num_new_ids = 0;
      fltarr *newpos = NULL;
      if (outpos)	
	newpos = (fltarr*) calloc (npart, sizeof(IDtype));	
      IDtype *newids = (IDtype*) calloc (npart, sizeof(IDtype));
      while (dist < (trace_factor * rvir))
	{
	  dist = wpart[i].dist;
	  IDtype *fnd;
	  IDtype id_key = wpart[i].part.id;
	  fnd = bsearch( &id_key, idlist, numids, sizeof(IDtype), cmp_IDtype);
	  if (fnd == NULL)
	    {
	      if (outpos)
		{
		  for (j = 0; j < 3; j++)
		    newpos[num_new_ids][j] = wpart[i].part.pos[j] / head.boxsize;
		}
	      newids[num_new_ids++] = id_key;
	    }
	  i++;
          if (i == npart)
                break;
	}
      free(wpart);
      maxdist = dist;
      if (verbose)
	printf("haloid %d |#%d particles added | maxdist: %g \n", haloid, num_new_ids, maxdist);

      fp = fopen(idlistname, "w");
      int totnumids = numids + num_new_ids;
      fwrite(&totnumids, sizeof(int), 1, fp);
      fwrite(center, sizeof(float), 3, fp);
      fwrite(&maxdist, sizeof(float), 1, fp);
      if (numids)
	fwrite(&idlist[0], sizeof(IDtype), numids, fp);      
      if (num_new_ids)
	fwrite(&newids[0], sizeof(IDtype), num_new_ids, fp);
      fclose(fp);

      if (outpos)
	{
	  fp = fopen(posfilename, "w");
	  fwrite(&totnumids, sizeof(int), 1, fp);
	  if (numids)
	    fwrite(&pos[0], sizeof(fltarr), numids, fp);      
	  if (num_new_ids)
	    fwrite(&newpos[0], sizeof(fltarr), num_new_ids, fp);
	  fclose(fp);
	}

      if (numids)
	{
	  free(idlist);
	  if (outpos)
	    free(pos);
	}
      if (npart)
	{
	  free(newids);
	  if (outpos)
	    free(newpos);
	}
    }

  if (write_catalogue)
  {	
	  fp =  fopen(catname, "w");
	  for (haloid = start_id; haloid <= max_halo_id; haloid++ )
	  {
	  	fprintf(fp, "%s", output[haloid]);
	  }
	  fclose(fp);
  }

  for (i = 0; i < (max_halo_id+1); i++ )
        free(output[i]);
  free(output);
  return 0;
}
