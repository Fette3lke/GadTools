/*
FINDCENTER v0.01

Program to find the center of mass of friends-of-friends groups found with mfof in a gadget simulation

icc -lm -openmp -o ../bin/altix/findcenter findcenter.c libgad.o 

gcc -std=c99 -lm -fopenmp -lgad -lgsl -lgslcblas -o ../bin/findcenter findcenter.c libgad.o
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <time.h>
#include "libgad.h"
#define	PI 3.14159265358979323846
#define GAP2BORDER 9000                                       //min distance to boarders to ignore periodic boundaries
#define OD_RAD 4000                                  //Radius for calculation of environmental density (kpc/h)
#define INCLUDE 0                                    //(friends + INCLUDE*kpc/h) are used to determine CM
#define NBOX 6                                       //Number of adjacent boxes that are used in every direction
#define GRIDSIZE 1000
#define h0 0.72
#define	MIN(a, b)   ((a)<(b)?(a):(b))
#define	MAX(a, b)   ((a)>(b)?(a):(b))
#define ABS(a) ((a) >= 0 ? (a) : -(a))
#define CMP(a,b) ((a)>(b)?(1):(-1))
#define PB(a,b) ((a)>(b)?(a-b):(a))
#define MOVE(a,b) PB(a+b/2,b)
#define MV(a,b) ((a)+(b)/2)%(b)
#define MOVEB(a) MOVE((a),boxsz)
//#define SQR(x) (x)*(x)
#define SOFTENING 5.00
#define G 6.6742e-11
#define Msun 1.989e30
#define kpc 3.08567758128e19

// #ifdef LONGIDS
// typedef unsigned long long IDtype;
// #else 
// typedef unsigned int IDtype;
// #endif

struct rank {
  int ind;
  int cnt;
};

void usage()
{
  fprintf(stderr, "Findcenter v0.01\n");
  fprintf(stderr, "-f <friendsfile>\n");
  fprintf(stderr, "-s <gadget snapshotfile>\n");
  fprintf(stderr, "-o <outputfile>\n");
  fprintf(stderr, "-v [verbose]\n");
  fprintf(stderr, "-gap [min-dist to boundaries for PB]\n");
  fprintf(stderr, "\n");
  exit(0);
}

int cmp_hcnt (const void *first, const void *second)
{
  struct rank *a = (struct rank *) first;
  struct rank *b = (struct rank *) second;
  if (a->cnt > b->cnt) return -1;
  else if (a->cnt < b->cnt) return +1;
  else return 0;
}


int main  (int argc, char *argv[])
{
  struct rank *halo;
  char friendsfile[256], gadgetfile[256], outfile[256];
  FILE *fp;
  struct header head;
  struct gadpart *part;
  int i,j,k,l,m,n,q, idum;
  unsigned int numpart, ndm, ngas, checknpart;
  int nhalo, minnum;
  int *iclus, *hcnt, **indlist;
  IDtype **idlist;
  float fdum;
  double ddum;
  int verbose=0, pnum=0, usegas=0;
  double boxsz;
  double GAP=GAP2BORDER;
  IDtype minID;


/***********************************************************************************
  START

***********************************************************************************/

  i=1;
  strcpy(outfile,"cm.txt");
  //  GAP=0;
  if (argc==1) usage();
  while (i<argc)
    {
      
      if (!strcmp(argv[i],"-f"))
	{
	  i++;
	  strcpy(friendsfile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-s"))
	{
	  i++;
	  strcpy(gadgetfile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-o"))
	{
	  i++;
	  strcpy(outfile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-gap"))
	{
	  i++;
	  GAP=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-v"))
	{
	  i++;
	  verbose=1;
	}
      else if (!strcmp(argv[i],"-gas"))
	{
	  i++;
	  usegas=1;
	}
      else if (!strcmp(argv[i],"-n"))
	{
	  i++;
	  pnum=1;
	}
      else usage();
    }

  

  //READ Gadget snapshot file

//#pragma omp parallel sections num_threads(2) private( i, j, k)
// {  
//#pragma omp section 
//   {

  if (!(numpart=readgadget_part(gadgetfile, &head, &part))) 
    {
      extern int libgaderr;
      printf("error reading gadgetfile, %d\n", libgaderr);
      exit(1);
    }
  ndm =head.npart[1];
  ngas=head.npart[0];
  boxsz=head.boxsize;
//   }
//#pragma omp section
//   {
     // READ friends-file
     fp=fopen(friendsfile,"r");                                         //read friendsfile
     fscanf(fp,"%d %lf %d %f %d", &checknpart, &ddum, &nhalo, &fdum, &minnum);
     if (ABS(ddum-head.time)>1e-4) printf("Possible Snapshot-Friendsfile Mismatch (Time, %g != %g)!\n", ddum, head.time);
     iclus=(int *)calloc(checknpart, sizeof(int));
     hcnt = (int *) calloc(nhalo, sizeof(int));
     indlist = (int **) malloc(nhalo * sizeof(int *));
     idlist = (IDtype **) malloc(nhalo * sizeof(IDtype *));
     halo= (struct rank *) malloc (nhalo * sizeof(struct rank));

     for (i=0; i < checknpart; i++)
       {
	 fscanf(fp,"%d %d", &idum, &iclus[i]);
	 if ((iclus[i]>nhalo) || (iclus[i]<0)) {printf("Error in friendsfile!\n"); exit(1);}
	 if (iclus[i]!=0) hcnt[iclus[i]-1]++;
       }
     fclose(fp);
     
     for (i=0; i< nhalo; i++)
       {
	 halo[i].ind=i;
	 halo[i].cnt=hcnt[i];
	 indlist[i] = (int *) malloc (hcnt[i]* sizeof(int));
	 if (indlist[i]==NULL) {printf("malloc failed (indlist)\n");exit(1);}
	 idlist[i] = (IDtype *) calloc (hcnt[i], sizeof(IDtype));
	 if (idlist[i]==NULL) {printf("malloc failed (idlist)\n");exit(1);}

       }
     qsort(halo, nhalo, sizeof(struct rank), cmp_hcnt);
//   }
// }

 {
     int *nind;
     nind = (int *) calloc(nhalo, sizeof(int));
     minID = part[0].id;
     for (i=0; i < checknpart; i++)
       {
	 if (iclus[i]!=0)
	   {
	     idlist[iclus[i]-1][nind[iclus[i]-1]] = part[i+ngas].id;
	     if (part[i+ngas].id < minID) minID = part[i+ngas].id;
	     if (part[i+ngas].id == 259725696)
	       {
		 printf("!#!#!#!#!#!#!\n");
	       }
	     indlist[iclus[i]-1][nind[iclus[i]-1]++]= i+ngas;	     
	   }
       }
     free(nind);
 }  
 if (!usegas)
   {
     if (checknpart!=ndm) printf("Possible Snapshot-Friendsfile Mismatch (particle number, ndm) %lu != %lu!\n", checknpart, ndm);
   }
 else
   {
     if (checknpart!=ngas) printf("Possible Snapshot-Friendsfile Mismatch (particle number, ngas) %lu != %lu!\n", checknpart, ngas);
   }

   if (verbose) {
#ifdef _OPENMP
#pragma omp parallel private(i)
     {
       i= omp_get_thread_num();
       //    printf("this is thread number %d \n", i);
       //    #pragma omp barrier
       if (i==0)
	 {
	   printf ("Number of threads: %d\n", omp_get_num_threads());
	 }
     } 
#endif
   }

   if (verbose) {printf("Number of fof-groups %d\n", nhalo);fflush(stdout);} 

/************************************************************************************/
/* Start of main loop                                                               */
/************************************************************************************/
#pragma omp parallel for private(i, j, k, l, m, n, q)
   //for (int tr_halo=0; tr_halo < nhalo; tr_halo++)
for ( i=0; i < nhalo; i++)
     {
       int tr_halo = halo[i].ind;
       if (hcnt ==NULL) {printf("malloc failed (hcnt)\n");exit(1);}
       if (halo == NULL) {printf("malloc failed (hdata)\n");exit(1);}
       
       l=indlist[tr_halo][0];
       char idlistname[128];
       sprintf(idlistname, "idlist_%d", i);
       FILE *fp = fopen(idlistname, "w");
       if (fp == NULL) continue;
       float maxdist = 0;
       fwrite(&hcnt[tr_halo], sizeof(int), 1, fp);
       fwrite(&(part[l].pos[0]), sizeof(float), 3, fp);
       fwrite(&maxdist, sizeof(float), 1, fp);
       fwrite(&(idlist[tr_halo][0]), sizeof(IDtype), hcnt[tr_halo], fp);
       fclose(fp);

     }

 printf("minID %lu\n", minID);

  return 0;
}
