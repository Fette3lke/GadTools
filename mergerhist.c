/*
Program to extract the stellar merger history from gadget files created with trace.c

gcc -std=c99 -lm -o ~/bin/altix/mergerhist mergerhist.c -lgad-altix KdTree.o -lgsl -lgslcblas

stan:
gcc -std=c99 -lm -o ~/bin/mergerhist mergerhist.c -lgad-stan -lgsl -lgslcblas stan/KdTree.o

gcc -std=c99 -lm -o ~/bin/mergerhist mergerhist.c -lgad-stan -lgsl -lgslcblas -fopenmp stan/KdTree-omp.o stan/ompfuncs.o

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "libgad.h"
#include "KdTree.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define NMAXSTARS 1000000
#define ACC_FRAC 0.15
#define FOF_LIMIT 20
#define NMAX_MERGER 1000
#define NMAX_SNAPSHOTS 100
#define LINKING_LENGTH 0.8                //can be overruled at runtime with -l

void usage()
{
  	  fprintf(stderr," v0.02\n");
  	  fprintf(stderr," mergerhist <halo name>\n");
	  exit(1);
}

struct star{
  int id;
  float dist;
  float gasdist;
  float frac;
  float idist;
  float ifrac;
  float ifrac_rhalf;
  int isnap;
  int snap;
  float a;
  int fnd;
  float a_acc;
  float pot;
  int index;
};

typedef struct Merger{
  int id;
  float start;
  float end;
  int n;
  float mass;
  float gmass;
  double gdisp;
  float ratio;
  int cid;
  double J[3];
  double Jint[3];
  double angle;
  double disp;
  double mean_sqr_vel;
} Merger;

int cmp_star_id(const void *a, const void *b)
{
  struct star *x= (struct star*)a;
  struct star *y= (struct star*)b;
  if (x->id > y->id) return 1;
  if (x->id < y->id) return -1;
  return 0;
}


int main (int argc, char *argv[])
{
  FILE *fp, *mergerfile, *partfile;
  int mergerind = 1;
  char infile[256], filename[256], haloname[256];
  int i,j,k;
  int startind=0;
  int endind=94;
  int snapind;
  struct gadpart *part, *wpart;
  struct header head;
  float fdum;
  int nacc=0;
  int idum=0;
  int debug = 0;
  double ll=LINKING_LENGTH;
  float fac =1;
  int imerger=0;
  int igm=0;
  float acc_frac = ACC_FRAC;
  double galmass[NMAX_SNAPSHOTS];
  double galrad[NMAX_SNAPSHOTS];
  double gal_disp[NMAX_SNAPSHOTS];
  double Jgal[NMAX_SNAPSHOTS][3];
  Merger *merger = (Merger*) malloc (sizeof(Merger) * NMAX_MERGER);
  float acclim = 0;

  i=1;
  if (1==argc) usage();
  while (i<argc)
    {
      if (!strcmp(argv[i],"-n"))
	{
	  i++;
	  strcpy(haloname,argv[i]);
	  i++;
	} 
      else if (*argv[i]!='-')
	{
	  strcpy(haloname,argv[i]);
	  i++;
      }
      else if (!strcmp(argv[i],"-l"))
	{
	  i++;
	  ll = atof(argv[i++]);
	}
      else if (!strcmp(argv[i],"-f"))
	{
	  i++;
	  fac = atof(argv[i++]);
	}
      else if (!strcmp(argv[i],"-af"))
	{
	  i++;
	  acc_frac = atof(argv[i++]);
	}
      else if (!strcmp(argv[i],"-al"))
	{
	  i++;
	  acclim = atof(argv[i++]);
	}

      else if (!strcmp(argv[i],"-test"))
	{
	  i++;
	  idum = atoi(argv[i++]);
	}
      else if (!strcmp(argv[i],"-debug"))
	{
	  i++;
	  debug=1;
	}
      else if (!strcmp(argv[i],"-i"))
	{
	  i++;
	  startind=atoi(argv[i]);
	  i++;
	  endind=atoi(argv[i]);
	  i++;
	  if (startind>endind)
	    {
	      int intdum=startind;
	      startind=endind;
	      endind=intdum;
	    }
	}
      else {
	usage();
      }
    }

  if (debug) printf("starting...\n");

  unsigned int numpart;
  float rvir = 0;
  sprintf(filename, "%s_094.idl", haloname);
  fp = fopen(filename, "r");
  if (fp == NULL)
    {
      fprintf(stderr, "%s not found\n", filename);
      //exit(1);
    }
  else
    {
      fscanf(fp, "%f %f %f %f %f %f", &fdum, &fdum, &rvir, &fdum, &fdum, &fdum);
      fclose(fp);         
      if (rvir)
	{
	  acclim = rvir * acc_frac;
	}
    }
  struct star dumstar;
  struct star *stars = (struct star*) malloc (sizeof(struct star) * NMAXSTARS );
  
  sprintf(filename,"%s.stars.insitu",haloname);
  fp = fopen(filename, "r");
  if (fp == NULL)
    {
      fprintf(stderr, "Insitu file not found\n");
      exit(1);
    }

if (debug) printf("reading star input file...\n");
  while (!feof(fp))
    {
      fscanf(fp,"%d %d %f %f %f %f %d %f %d %f %f %f %f\n", &stars[nacc].id, &stars[nacc].isnap, &stars[nacc].idist, &stars[nacc].dist, &stars[nacc].frac, &stars[nacc].gasdist, &stars[nacc].snap, &stars[nacc].a,  &stars[nacc].fnd, &stars[nacc].a_acc, &stars[nacc].pot, &stars[nacc].ifrac, &stars[nacc].ifrac_rhalf);
      //      if (debug) printf("nacc %6d\n", nacc);
      if (!acclim)
	{
	  if (stars[nacc].isnap >= 94)
	    {
	      acclim = (stars[nacc].idist / stars[nacc].ifrac)  * acc_frac;
	    }
	  else continue;
	}
      if ((stars[nacc].isnap >= 94) && (stars[nacc].idist <= (acclim)))
	{
	  if (( stars[nacc].frac > acc_frac ) || ( stars[nacc].fnd == 0)) nacc++;
	}
      
      if (nacc >= NMAXSTARS)
	{
	  fprintf(stderr, "Increase NMAXSTARS and recompile\n");
	  exit(1);
	}
    }
  fclose(fp);
  if (debug) printf("%d stars in input file | < %f..\n", nacc, acclim);

#ifdef _OPENMP
#pragma omp parallel private(i)
      {
	i= omp_get_thread_num();
	//    printf("this is thread number %d \n", i);
	//    #pragma omp barrier
	if (i==0)
	  {
	    if (debug) printf ("Number of threads: %d\n", omp_get_num_threads());
	  }
      } 
#endif



  if (debug) printf("sorting stars...\n");
  //  qsort(stars, nacc, sizeof(struct star), (__compar_fn_t)cmp_star_id);
  qsort(stars, nacc, sizeof(struct star), cmp_star_id);

  /*********************************************************************

      Program code goes here

  *********************************************************************/
  wpart = (gadpart *) malloc (sizeof(gadpart) * nacc);      
  double smooth=0;
  sprintf(filename,"%s.merger_part.dat", haloname);
  partfile = fopen(filename, "w");

  if ((partfile == NULL))
    {
      fprintf(stderr, "could not create file %s!\n", filename);
      exit(1);
    }

  if (endind - startind >=NMAX_SNAPSHOTS)
    {
      fprintf(stderr, "recompile with larger NMAX_SNAPSHOTS");
      exit(1);
    }

  if (debug) printf("comencing main loop\n");

  //    for (snapind=endind; snapind>=startind; snapind--)
  for (snapind=startind; snapind<=endind; snapind++)
    { 
      if (debug) printf("snapshot %3d\n", snapind);
      sprintf(infile, "%s_%03d.gad", haloname, snapind);
      //      printf("reading %s\n", infile);
      if (!(numpart=readgadget_part(infile, &head, &part))) 
	{
	  extern int libgaderr;
	  if (libgaderr==1)
	    {
	      fprintf(stderr, "not found, skipping: %s\n",infile);
	      continue;
	    }
	  else
	    {
	      fprintf(stderr, "error in file %s\nLibgad-Error Code: %d\n",infile, libgaderr);
	      exit(1);
	    }
	}

      int nstars = head.npart[4];
      int starind = head.npart[0] + head.npart[1] + head.npart[2] + head.npart[3];
      

      //convert to physical units:
      for ( i = 0; i < numpart; i++) 
	for ( j = 0; j < 3; j++)
	  {
	    part[i].vel[j] *= sqrt(head.time);
	    part[i].pos[j] *= head.time;
	  }



      //Calculate central galaxy mass
      igm= snapind-startind;
      if (debug) printf("buildung tree...\n");
      KdNode *all_root;
      initKdNode(&all_root, NULL);
      buildKdTree(all_root, &part[starind], nstars, 0);
      //      printf("checktree: %d\n", checkKdTree(all_root));
      //      buildKdTree(all_root, part, numpart, 0);
      //printf("before %d %d\n", checkKdTree(all_root), nstars);
      gadpart **galpart = NULL;
      int ngal = 0;      
      unsigned int buffer = 0;
      fltarr center ={0.0, 0.0, 0.0};

      if (debug) {printf("find central FOF group...\n"); fflush(stdout);}
      //      findGadparts(all_root, center, fac *  ll * head.time, &galpart, &ngal, &buffer);
      //      qsort(galpart, ngal, sizeof(gadpart *),(__compar_fn_t) cmp_pointer_id );
      if (debug) {printf("sorting...\n"); fflush(stdout);}
      //      qsort(galpart, ngal, sizeof(gadpart *), cmp_pointer_id );
      //      printf("FOF\n"); fflush(stdout);
      //      printf("before  %d %d\n", checkKdTree(all_root), nstars);fflush(stdout);
      if (debug) {printf("FOF search...\n"); fflush(stdout);}
      findFOF(all_root, center, ll * head.time, &galpart, &ngal, &buffer);
      //            printf("after  %d %d\n", checkKdTree(all_root), nstars);fflush(stdout);
      delKdNode(&all_root);

      if (debug) {printf("central FOF found...%d\n", ngal); fflush(stdout);}
   
      //qsort(galpart, ngal, sizeof(gadpart*), (__compar_fn_t) cmp_pointer_id );
      galmass[igm] = 0;
      for ( k=0; k<3; k++) Jgal[igm][k] = 0;
      double avg_vel[3] = {0.0, 0.0, 0.0};
      double avg_sqr_vel = 0;
      gadpart_dist *gal = (gadpart_dist *) malloc (sizeof (gadpart_dist) * ngal);
      for ( i = 0; i < ngal; i++ )
	{
	  galmass[igm] += galpart[i] -> mass;
	  gal[i].part = *(galpart[i]);
	  gal[i].dist = distance_nopb(galpart[i] -> pos, center);
	  //	  printf("%d\n", galpart[i] -> id);
	  double sqrdum = 0;
	  for ( j = 0; j < 3; j++)
	    {
	      avg_vel[j] += galpart[i] -> vel[j];
	      avg_sqr_vel += galpart[i] -> vel[j] * galpart[i] -> vel[j];
	    }

	  Jgal[igm][0] += galpart[i] -> mass * (galpart[i] -> pos[1] * galpart[i] -> vel[2] - galpart[i] -> pos[2] * galpart[i] -> vel[1]);
	  Jgal[igm][1] += galpart[i] -> mass * (galpart[i] -> pos[2] * galpart[i] -> vel[0] - galpart[i] -> pos[0] * galpart[i] -> vel[2]);
	  Jgal[igm][2] += galpart[i] -> mass * (galpart[i] -> pos[0] * galpart[i] -> vel[1] - galpart[i] -> pos[1] * galpart[i] -> vel[0]);
	}
      qsort(gal, ngal, sizeof(gadpart_dist), cmp_dist);
      double massdum = 0;
      i=0;
      while (massdum < (galmass[igm]/2.))
	{
	  galrad[igm] = gal[i].dist;
	  massdum += gal[i].part.mass;
	  i++;
	}
      free(gal);
      
      for ( j = 0; j < 3; j++) avg_vel[j] /= ngal;
      avg_sqr_vel /= ngal;
      gal_disp[igm] = SQR(avg_vel[0]) + SQR(avg_vel[1]) + SQR(avg_vel[2]);
      gal_disp[igm] = sqrt(gal_disp[igm]);
      gal_disp[igm] = avg_sqr_vel - gal_disp[igm] * gal_disp[igm];
      gal_disp[igm] = sqrt(gal_disp[igm]);

#pragma omp parallel for
      for ( i = 0; i < imerger; i++)
	{
	  if ((merger[i].cid == 0) || (merger[i].end != 0) ) continue;
	  gadpart dumpart;
	  dumpart.id = merger[i].cid;
	  gadpart *dumpointer = &dumpart;
	  gadpart **fnd;
	  //	  fnd = bsearch(&dumpointer ,galpart, ngal, sizeof(gadpart*), (__compar_fn_t) cmp_pointer_id );
	  fnd = bsearch(&dumpointer ,galpart, ngal, sizeof(gadpart*), cmp_pointer_id );
	  if ((fnd != NULL) && ( merger[i].end == 0 ))
	    {
	      merger[i].end   = head.time;
	      merger[i].ratio = merger[i].mass / galmass[igm - 1];
	    }
	}


      //Search for infalling structures

      //      qsort(&part[starind], nstars, sizeof(gadpart), (__compar_fn_t) cmp_id );
      //      qsort(&part[0], numpart, sizeof(gadpart), (__compar_fn_t) cmp_id );
      qsort(&part[0], numpart, sizeof(gadpart),  cmp_id );

      //      printf("search\n");fflush(stdout);
      j=0;
      if (debug) printf("looking up stars...\n");
      for ( i = 0; i < nacc; i++ )
	{	  
	  gadpart dumpart;
	  dumpart.id = stars[i].id;
	  gadpart *fnd;
//	  fnd = bsearch (&dumpart, &part[0], numpart, sizeof(gadpart), (__compar_fn_t) cmp_id );
	  fnd = bsearch (&dumpart, &part[0], numpart, sizeof(gadpart), cmp_id );
	  if (fnd != NULL)
	    {
	      if (fnd -> type == 4)
		{
		  wpart[j] = *fnd;
		  j++;
		}
	    }	  
	}
      int starsfnd = j;
      //qsort(&wpart[0], starsfnd, sizeof(gadpart), (__compar_fn_t) cmp_id );
      
      free(part);
 
      //            printf("fnd: %d\n", j);      
      if (starsfnd == 0) 
	{
	  if (debug) printf("no stars found, next snapshot..\n");
	  continue;
	}
      else
	{
	  if (debug) printf("found %d stars\n", starsfnd);
	}

      if (debug) printf("buildung stellar tree...\n");
      KdNode *root;
      initKdNode(&root, NULL);
      buildKdTree(root, wpart, starsfnd, 0);
      

      for ( i = 0; i < nacc; i++) stars[i].index = -1;

      if (debug) printf("searching substructures...\n");
#pragma omp parallel for 
      for ( i = 0; i < starsfnd; i++)
	{
	  struct star stardum, *fnd;
	  stardum.id = wpart[i].id;
	  //	  fnd = bsearch(&stardum, stars, nacc, sizeof(struct star), (__compar_fn_t)cmp_star_id);
	  fnd = bsearch(&stardum, stars, nacc, sizeof(struct star), cmp_star_id);
	  if (fnd != NULL)
	    {
	      fnd -> index = i;	      
	    }
	}

      double smoothadded = 0;
      int nadded = 0;

      if (debug) printf("aexpn %f\n", head.time);
      if (debug) { printf("check Tree...%d ?= %d\n", checkKdTree(root), starsfnd );fflush(stdout); }
      for ( i = 0; i < nacc; i++) 
	{
	  if (stars[i].index < 0)
	    {
	      //	      printf("star not found in snapshot. %d\n", i);
	      continue;
	    }

	  //if ( stars[i].a_acc <= (head.time+0.0001) ) continue;
	  if (( (stars[i].a_acc - 0.001) > (head.time) ) || ( stars[i].a_acc == 0 ) ) continue;
	  //	  if (( (stars[i].a_acc - 0.011) > (head.time) ) || ( stars[i].a_acc == 0 ) ) continue;
	  //	  if ( stars[i].a_acc > (head.time) )  continue;

	  gadpart **FOF = NULL;
	  int nFOF = 1;
	  unsigned int size = 10000;
	  FOF = (gadpart **) malloc (size * sizeof(gadpart*));
	  FOF[0] = &wpart[stars[i].index];

	  findFOF(root, FOF[0]->pos, ll * head.time, &FOF, &nFOF, &size);

	  double mergermass = 0;
	  double cm[3];
	  double cvel[3];
	  int cm_id;
	  int currentmerger_id = 0;
	  for (k=0; k<3; k++) 
	    {
	      cm[k]   = 0;
	      cvel[k] = 0;
	    }
	  for ( j = 0; j < nFOF; j++)
	    {
	      for (k=0; k<3; k++)
		{
		  cm[k]   += FOF[j] -> pos[k] * FOF[j] -> mass;
		  cvel[k] += FOF[j] -> vel[k] * FOF[j] -> mass;
		}
	      mergermass += FOF[j]->mass;
	    }
	  for (k=0; k<3; k++) cm[k]   /=  mergermass;
	  for (k=0; k<3; k++) cvel[k] /=  mergermass;

	  
	  if (nFOF >= FOF_LIMIT)
	    {
	      if (debug) { printf("nFOF > limit: check Tree...%d ?= %d\n", checkKdTree(root), starsfnd );fflush(stdout); }
	      if (debug) 
		{
		  printf("FOF group found...%d | # %d\n", mergerind, nFOF);
		  printf("central star: id %d | a_acc %f \n",stars[i].id, stars[i].a_acc);
		}

//	      for (k=0; k<3; k++) 
//		{
//		  printf("%f %f\n", cvel[k], cm[k]);
//		}
	      gadpart dumpart;
	      double avg_vel[3]={0.0, 0.0, 0.0};
	      double avg_sqr_vel=0;

	      for (k=0; k<3; k++) dumpart.pos[k] = cm[k];
	      cm_id = findNN(root, &dumpart) -> part -> id;	      
	      merger[imerger].id    = mergerind;
	      merger[imerger].start = head.time;
	      merger[imerger].n     = nFOF;
	      merger[imerger].mass  = mergermass;
	      merger[imerger].gmass = galmass[igm];
	      merger[imerger].gdisp = gal_disp[igm];
	      merger[imerger].ratio = merger[imerger].mass / galmass[igm-1];
	      merger[imerger].disp  = 0;
	      merger[imerger].cid   = cm_id;
	      for (k=0; k<3; k++) 
		{
		  merger[imerger].J[k]    =  0;
		  merger[imerger].Jint[k] =  0;
		}
	      for ( j = 0; j < nFOF; j++)
		{		  
		  merger[imerger].J[0] += FOF[j] -> mass * (FOF[j] -> pos[1] * FOF[j] -> vel[2] - FOF[j] -> pos[2] * FOF[j] -> vel[1]);
		  merger[imerger].J[1] += FOF[j] -> mass * (FOF[j] -> pos[2] * FOF[j] -> vel[0] - FOF[j] -> pos[0] * FOF[j] -> vel[2]);
		  merger[imerger].J[2] += FOF[j] -> mass * (FOF[j] -> pos[0] * FOF[j] -> vel[1] - FOF[j] -> pos[1] * FOF[j] -> vel[0]);
		  double veldum[3];
		  double posdum[3];
		  for (k=0; k<3; k++) 
		    {
		      veldum[k] = FOF[j] -> vel[k] - cvel[k];
		      posdum[k] = FOF[j] -> pos[k] - cm[k];

		      avg_vel[k]  += veldum[k];
		      avg_sqr_vel += veldum[k] * veldum[k];
		    }

		  merger[imerger].Jint[0] += FOF[j] -> mass * (posdum[1] * veldum[2] - posdum[2] * veldum[1]);
		  merger[imerger].Jint[1] += FOF[j] -> mass * (posdum[2] * veldum[0] - posdum[0] * veldum[2]);
		  merger[imerger].Jint[2] += FOF[j] -> mass * (posdum[0] * veldum[1] - posdum[1] * veldum[0]);
		  		  
		}
	      for (k=0; k<3; k++) avg_vel[k] /= nFOF;
	      avg_sqr_vel /= nFOF;
	      double ddum = SQR(avg_vel[0]) + SQR(avg_vel[1]) + SQR(avg_vel[2]);
	      ddum = sqrt(ddum);
	      merger[imerger].disp = sqrt(avg_sqr_vel - ddum * ddum);
	      
	      double sp = 0;
	      double abs[2] = {0, 0};
	      for (k=0; k<3; k++) 
		{
		  sp +=  merger[imerger].Jint[k] *  Jgal[igm][k];
		  abs[0] += merger[imerger].Jint[k] * merger[imerger].Jint[k];
		  abs[1] += Jgal[igm][k] * Jgal[igm][k];
		}

	      abs[0] = sqrt(abs[0]);
	      abs[1] = sqrt(abs[1]);

	      merger[imerger].angle = acos (sp / (abs[0] * abs[1]));
	      
	      currentmerger_id = mergerind++;
	      imerger++;	      
	      
	      if ((imerger +1) >= NMAX_MERGER)
		{
		  fprintf(stderr, "recompile with larger NMAX_MERGER");
		  exit(1);
		}

	    } else 
	    {
	      if ((debug) && (nFOF)) printf("below limit: %d | %f %f %f\n", nFOF, FOF[0]->pos[0], FOF[0]->pos[1], FOF[0]->pos[2]);
	      smoothadded += mergermass;
	      nadded += nFOF;
	    }
	  
	  int notfnd=0;
	  int last_fnd_id=-1;
	  
	  for ( j = 0; j < nFOF; j++)
	    {
	      struct star stardum, *fnd;
	      stardum.id = FOF[j] -> id;
	      //	      fnd = bsearch(&stardum, stars, nacc, sizeof(struct star), (__compar_fn_t)cmp_star_id);
	      fnd = bsearch(&stardum, stars, nacc, sizeof(struct star), cmp_star_id);
	      if (fnd != NULL)
		{
		  last_fnd_id= fnd -> id;
		  fprintf(partfile,"%7d %5d %3.2f\n", fnd -> id, currentmerger_id, head.time);		  
		  int index = (fnd - stars);
		  struct star stardum = stars[index];
		  for ( k = index; k < (nacc - 1); k++ )
		    {
		      stars[k] = stars[k + 1];
		    }
		  stars[nacc-1] = stardum;
		  nacc--;
		} else 
		{
		  notfnd++;
		}
	    }
	    
	  if ((debug) && (nFOF >= FOF_LIMIT)) 
	    {
//	      int ok = 0;
//	      for ( j = 0; j < nFOF; j++)
//		{
//		  if (FOF[j])
//		}
	      printf("FOF particles not found in list %d\n", notfnd); 
	      if ( (nFOF - notfnd) == 1)
		{
		  printf("last found id %d\n", last_fnd_id); 
		  for ( j = 0; j < nFOF; j++)
		    {
		      printf("%d %d | ", FOF[j] -> id, FOF[j] -> type);
		    }
		  printf("\n");
		}
	      fflush(stdout);
	    }
	  if (FOF != NULL)
	    free(FOF);
	}

      merger[imerger].id    = 0;
      merger[imerger].start = head.time;
      merger[imerger].end   = head.time;
      merger[imerger].n     = nadded;
      merger[imerger].mass  = smoothadded;
      merger[imerger].gmass = galmass[igm];
      merger[imerger].gdisp = gal_disp[igm];
      merger[imerger].ratio = merger[imerger].mass / galmass[igm-1];
      merger[imerger].cid   = 0;
      merger[imerger].angle = 0;
      merger[imerger].disp  = galrad[igm];
      int l;
      for ( l=0; l<3; l++) 
	{	 
	  merger[imerger].Jint[l] = 0;
	  //	  Jgal[igm][l] = 0;
	  merger[imerger].Jint[l] = Jgal[igm][l];
	}
      imerger++;
      smooth += smoothadded;

//      for (i = 0; i < nFOF; i++)
//	{
//	  printf("%6d %8.2f %8.2f %8.2f\n", FOF[i]->id, FOF[i]->pos[0],  FOF[i]->pos[1], FOF[i]->pos[2]);
//	}
      //      printf("nFOF %d %d\n", wpart[0].id, nFOF);
      
      delKdNode(&root);
    }

 free(wpart);      

 fclose(partfile);

 sprintf(filename,"%s.merger.dat", haloname);
 mergerfile = fopen(filename, "w");
  if ((mergerfile == NULL))
    {
      fprintf(stderr, "could not create file %s!\n", filename);
      exit(1);
    }
  for ( i = 0; i < imerger; i++)
    {
      fprintf(mergerfile, "%4.2f 0 0 %10.5g %10.5g %10.5g %4.2f %5d %7d %6d %10.5g %10.5g %10.5g %4.2f %4.2f %5.2f %5.2f %10.5g %10.5g %10.5g\n", merger[i].start , merger[i].gmass , merger[i].ratio , merger[i].mass,  merger[i].end , merger[i].cid, merger[i].n, merger[i].id, merger[i].Jint[0], merger[i].Jint[1], merger[i].Jint[2], merger[i].angle, (merger[i].angle * 360 / (2 * M_PI)), merger[i].gdisp, merger[i].disp, merger[i].J[0], merger[i].J[1], merger[i].J[2]);
    }
 // fprintf(mergerfile, "%4.2f %5d %10.5g %7d %6d\n", head.time, mergerind, mergermass, cm_id, nFOF);
 fclose(mergerfile);

 if (debug)
   {
     printf("%f %d %d\n", rvir, nacc, j);
     printf("%d %f\n", stars[nacc-1].id, stars[nacc-1].idist);
   }
 
 return 0;
}
