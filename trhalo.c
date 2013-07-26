/*
Program to trace back particles in a halo to the initial conditions
and/or to print distance of particles from the Halo center
requires gadget file at z=0 and initial conditions
requires center of mass (can be read in from a file created with virialmass.c)

gcc -std=c99 -lm -lgsl -lgslcblas -o trhalo trhalo.c libgad-gcc.o lmfit-2.2/lmmin.o lmfit-2.2/lm_eval.o KdTree.o

gcc -std=c99 -lm -lgsl -lgslcblas -o trhalo trhalo.c -lgad-altix KdTree.o
gcc -std=c99 -lm -lgsl -lgslcblas -o ../bin/altix/trhalo-nopot trhalo.c -lgad-altix-nopot KdTree.o -DNOPOT

stan:
gcc -lm -lgsl -lgslcblas -o ~/bin/trhalo trhalo.c -lgad-stan stan/KdTree.o 


*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include "libgad.h"
#include "KdTree.h"
#include "./lmmin.h"
#include "./lm_eval.h"
#define	PI 3.14159265358979323846
#define GAP2BORDER 16000                                     //min distance to boarders to ignore periodic boundaries
#define INCLUDE 0                                    //(friends + INCLUDE*kpc/h) are used to determine CM
#define TRACE_FACTOR 2.0                             //TRACE_FACTOR times virial radius will be traced back to IC
#define h0 0.72
#define VIROD 200 
#define DIM 512
#define ENVDENSRAD (4000 * h0)
#define	MIN(a, b)   ((a)<(b)?(a):(b))
#define	MAX(a, b)   ((a)>(b)?(a):(b))
#define ABS(a) ((a) >= 0 ? (a) : -(a))
#define CMP(a,b) ((a)>(b)?(1):(-1))
#define PB(a,b) ((a)>(b)?(a-b):(a))
#define MOVE(a,b) PB(a+b/2,b)
#define MOVEB(a) MOVE((a),boxsz)
#define MV(a,b) ((a)+(b)/2)%(b)
//#define SQR(x) (x)*(x)
#define SOFTENING 0
#define G 6.6742e-11
#define Msun 1.989e30
#define kpc 3.08567758128e19
#define DEBUG 1

#define VA_USE_PART_TYPE 2
#define CM_USE_PART_TYPE 2
#define INERTIA_USE_PART_TYPE 2
#define kNN 50

double cdens;

struct particle {int ind;float dist;};
clock_t t[2];


void usage()
{
  	  fprintf(stderr,"TrHalo v0.05\n");
	  fprintf(stderr," -f <snapshot file> -i <IC file> \n");
	  fprintf(stderr," -m <virialmass file> -n <Halo Number> or -cm <X Y Z>\n");
	  fprintf(stderr,"   if you want to search for CM in an area around the given coordinates\n");
	  fprintf(stderr,"   add -sr <radius> to include particles for the search with distance to cm < radius\n");
	  fprintf(stderr," -gr <gridfile> -tf <tracefactor default 2.0>\n");
	  fprintf(stderr," -jr <print j distribution to file> -cut <create gadget file of the halo>\n");
	  fprintf(stderr," -pb <force periodic boundary conditions (bitwise)>\n");
  	  fprintf(stderr," -n <name of the halo>\n");
  	  fprintf(stderr," -cid <center on particle with a given ID>\n");
  	  fprintf(stderr," -trid <file with particle IDs to trace>\n");
  	  fprintf(stderr," -r <fraction of virial radius for longest axis of Inertia Tensor>\n");
  	  fprintf(stderr," -ui <bitcode particle types for Inertia Tensor>\n");
  	  fprintf(stderr," -gap <min distance to Box-border for PB>\n");
  	  fprintf(stderr," -mcnt <min particle Count for CM Search>\n");
  	  fprintf(stderr," -sd <search distance>\n");

	  exit(1);
}


double fit_fct(double t, double *p)
{
  return log10((p[0]) / ((t/p[1]) * SQR(1+t/p[1])));
}

// void printout ( int n_par, double* par, int m_dat, double* fvec, 
//                        void *data, int iflag, int iter, int nfev )
//          {
//            // dummy function to catch fitting output
//          }


float step()
{
  float tm;
  t[1]=clock();
  tm=(float)(t[1]-t[0])/CLOCKS_PER_SEC;
  t[0]=clock();
  return tm;
}

/*
int cmp (struct particle *a, struct particle *b)
{
  if (a[0].dist > b[0].dist) return 1; else return -1;
}
*/
int cmp (const void *first, const void *second)
{
  struct particle *a = (struct particle*) first;
  struct particle *b = (struct particle*) second;
  if (a->dist > b->dist) return 1;
  else if (a->dist < b->dist) return -1;
  else return 0;
}

int cmpind (const void *first, const void *second)
{
  struct particle *a = (struct particle*) first;
  struct particle *b = (struct particle*) second;
  if (a->ind > b->ind) return 1;
  else if (a->ind < b->ind) return -1;
  else return 0;
}

int compare (float *a, float *b)
{
  if (*a > *b) return 1; else return -1;
}

int coord512 (int *result, long index)
{
  int c[3];
  int i;
  result[0]=floor(index%262144/512.0);
  result[1]=(index%512);
  result[2]=floor(index/262144.0);
  for (i=0; i<3; i++)
    if ((result[i]<0) || (result[i]>512)) return 1;
  return 0;
}


int main  (int argc, char *argv[])
{
  typedef float fltarr[3];
  struct particle *part, *vpart;
  gadpart *P;
  struct header ic,evolved, out;
  FILE *fp;
  char vmfile[256],icfile[256],evolvedfile[256],snfile[256],friendsfile[256],idlfile[256],gridfile[256],jrfile[256], gadfile[256], hname[256], parfile[256], idfile[256], datfile[256];
  unsigned int numpart=0, nummass=0;
  int blocksize,i,j,k,l,m,dum,dum2,tr_halo=0,tr_halo_cnt = 0,id,halo,checknpart,mindum,nhalo,count,notrace=1;
  int *tr_halo_id,*tr_halo_i,*tr_halo_i_dum, ca[3], size[3];
  //  fltarr *pos_ev,*pos_ic,*vel, *outpos, *outvel;
  fltarr *pos_ic;
  fltarr shift, diff;
  float *mass_ev, *mass_dum,*mass_ic, *outmass,dist,maxdist,*distance,halomass,halorad, halolmd;
  double posdum,max[3],min[3],srad=0, icmin[3], icmax[3], scm[3];
  double lmd,vcm[3]={0,0,0},jcm[3]={0,0,0},J[3]={0,0,0}, torq[3],rdum[3], vdum[3], massdum;
  int *id_ev, *outid,*iclus, halonpart;
  float tm,linkle,od, tr_factor;
  double mdens,boxsz,cm[3]={0,0,0}, expcm[3], masstot, gridsize, rlog[50], p[50], err[50], d, rad_ratio=0.3;
  short pb[3]={0,0,0};
  int minbox[3], maxbox[3], szbox[3], ***grid, boxind[3], alignf, mincnt=100;
  int nfun, npar, indmax[3], indmin[3], dojr=0, docut=0, vcnt, starsonly=0, dmonly=1, doidl=1, start, end, use_inertia=0, use_va, use_cm, cmset=0, nopb=0, cid=0, traceid=0, gal_all=0, addbh=0, printcm=0;
  double par[2],ddum, dx512, ddum1, totj, envdens, sqrenvdensrad, bgdens, sfl=SOFTENING, meff,  halfmrad, effrad;
  double GAP=GAP2BORDER, searchdist=0;
 
  t[0]=clock();
  dx512=72000.0/512.0;
  
  if (DEBUG) printf("DEBUG-Output enabled\n");fflush(stdout);

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //filenames are ignored if given as command line parameters in the same order

  strcpy(vmfile,"<none>");
  strcpy(gridfile,"<none>");
  strcpy(parfile,"halodata.txt");
  strcpy(evolvedfile,"");
  strcpy(icfile,"");
  strcpy(idlfile,"<none>");
  tr_factor=TRACE_FACTOR;
  use_inertia=INERTIA_USE_PART_TYPE;
  use_va=VA_USE_PART_TYPE;
  use_cm=CM_USE_PART_TYPE;
  //  tr_halo=123;
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  i=1;
  if (1==argc) usage();
  while (i<argc)
    {
      
      if (!strcmp(argv[i],"-f"))
	{
	  i++;
	  strcpy(evolvedfile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-i"))
	{
	  i++;
	  strcpy(icfile,argv[i]);
	  i++;
	  notrace=0;
	}
      else if (!strcmp(argv[i],"-m"))
	{
	  i++;
	  strcpy(vmfile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-idl"))
	{
	  i++;
	  strcpy(idlfile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-trid"))
	{
	  i++;
	  strcpy(idfile,argv[i]);
	  traceid=1;
	  i++;
	}
      else if (!strcmp(argv[i],"-n"))
	{
	  i++;
	  strcpy(hname,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-gr"))
	{
	  i++;
	  strcpy(gridfile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-n"))
	{
	  i++;
	  tr_halo=atoi(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-sr"))
	{
	  i++;
	  srad=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-cid"))
	{
	  i++;
	  cid=atoi(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-gap"))
	{
	  i++;
	  GAP=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-r"))
	{
	  i++;
	  rad_ratio=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-ui"))
	{
	  i++;
	  use_inertia=atoi(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-uva"))
	{
	  i++;
	  use_va=atoi(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-ucm"))
	{
	  i++;
	  use_cm=atoi(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-use"))
	{
	  i++;
	  use_cm=atoi(argv[i]);
	  use_va=use_cm;
	  use_inertia=use_cm;
	  i++;
	}
      else if (!strcmp(argv[i],"-soft"))
	{
	  i++;
	  sfl=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-sd"))
	{
	  i++;
	  searchdist=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-mcnt"))
	{
	  i++;
	  mincnt=atoi(argv[i]);
	  i++;
	}
       else if (!strcmp(argv[i],"-jr"))
	{
	  i++;
	  dojr=1;
	}
       else if (!strcmp(argv[i],"-bh"))
	{
	  i++;
	  addbh=1;
	}
       else if (!strcmp(argv[i],"-noidl"))
	{
	  i++;
	  doidl=0;
	}
      else if (!strcmp(argv[i],"-pcm"))
	{
	  i++;
	  printcm=1;
	}
      else if (!strcmp(argv[i],"-nopb"))
	{
	  i++;
	  nopb=1;
	}
      else if (!strcmp(argv[i],"-ga"))
	{
	  i++;
	  gal_all=1;
	}
      else if (!strcmp(argv[i],"-cut"))
	{
	  i++;
	  docut=1;
	}
      else if (!strcmp(argv[i],"-pb"))
	{
	  i++;
	  dum=atoi(argv[i]);
	  for (j=0; j<3; j++) if (dum&(1<<j)) pb[k]=1;
	  i++;
	}
      else if (!strcmp(argv[i],"-stars"))
	{
	  i++;
	  starsonly=1;
	}
      else if (!strcmp(argv[i],"-all"))
	{
	  i++;
	  dmonly=0;
	}
      else if (!strcmp(argv[i],"-tr"))
	{
	  i++;
          notrace=1;
	}
      else if (!strcmp(argv[i],"-tf"))
	{
	  i++;
	  tr_factor=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-cm"))
	{
	  i++;
	  cm[0]=atof(argv[i++]);
	  cm[1]=atof(argv[i++]);
	  cm[2]=atof(argv[i++]);
	  cmset=1;
	} else {
	  usage();
	}

    }
  
  if (strcmp(hname,""))
    {
      sprintf( idlfile, "%s.idl", hname);
      sprintf(gridfile, "%s.gr", hname);
      sprintf( datfile, "%s.dat", hname);
      sprintf(  jrfile, "%s.jr", hname);
      sprintf( gadfile, "%s.gad", hname);
    }

  if (!traceid)
    {

// auxiliary settings for fitting:

    lm_control_type control;
    lm_data_type data;
    lm_initialize_control(&control);

    data.user_func = fit_fct;
    data.user_t = rlog;
    data.user_y = p;
//---------------------------------------------------

    /*
   fp=fopen(evolvedfile,"r");

  fread(&blocksize,sizeof(int),1,fp);
  fread(&evolved,sizeof(struct header),1,fp);                        //read header of evolved gadget file
  fread(&blocksize,sizeof(int),1,fp);
  boxsz=evolved.boxsize;
  
  for (i=0; i<6; i++) 
    {
      numpart+=evolved.npart[i];
      if ((evolved.npart[i]!=0) || (evolved.massarr[i]==0)) nummass+=evolved.npart[i];
    }
  pos_ev =(fltarr *)malloc(sizeof(fltarr)*numpart+addbh);
  vel    =(fltarr *)malloc(sizeof(fltarr)*numpart+addbh);
  id_ev  =(int *)   malloc(sizeof(int)*numpart+addbh);
  mass_dum=(float *) malloc(sizeof(float)*nummass+addbh);
  mass_ev=(float *) malloc(sizeof(float)*numpart+addbh);

  fread(&blocksize,sizeof(int),1,fp);                                //read particle positions of evolved snapshot file
  fread(&P[0].pos,sizeof(fltarr),numpart,fp);
  fread(&blocksize,sizeof(int),1,fp);  
  
  fread(&blocksize,sizeof(int),1,fp);                                //skip velocities
  fread(&P[0].vel,sizeof(fltarr),numpart,fp);
  fread(&blocksize,sizeof(int),1,fp);

  fread(&blocksize,sizeof(int),1,fp);                                //read particle ids of evolved snapshot file
  fread(&P[0].id,sizeof(int),numpart,fp);
  fread(&blocksize,sizeof(int),1,fp);
  
  if (nummass>0)
    {
  fread(&blocksize,sizeof(int),1,fp);                                //read particle masses of evolved snapshot file
  fread(&mass_dum[0],sizeof(float),nummass,fp);
  fread(&blocksize,sizeof(int),1,fp);
    }

  fclose(fp);
*/
    numpart= readgadget_part(evolvedfile, &evolved, &P);
    if (numpart==0) 
	{
	  extern int libgaderr;
	  fprintf(stderr,"LibGad Error Code: %d\n", libgaderr);
	  exit(libgaderr);
	}
  masstot=0;
  j=0;
  for (i=0; i< numpart; i++) 
   {
     masstot+=P[i].mass;
     if ((cid) && P[i].id==cid )
       {
	 for (j=0; j<3; j++) cm[j]=P[i].pos[j];
       }
   } 
  mdens=masstot/pow(evolved.boxsize*evolved.time,3);

  cdens=(2.78e-8);
  cdens=cdens*(evolved.omegal+evolved.omega0*pow(evolved.time,-3));
  bgdens=cdens*evolved.omega0;
  dum=0;
  for (i=0; i<5; i++)
    {
      dum+=evolved.npart[i];
      if (DEBUG)
      if (evolved.npart[i]>0) printf("%d %f %f\n", i, P[dum-1].mass, P[dum].mass);
    }

  if (DEBUG) printf("readin evolved %.2f\n",step());fflush(stdout);
  printf("snapshot %s\n", evolvedfile);
  printf("ucm %d   ui %d   uva %d\n", use_cm, use_inertia, use_va);
  if (strcmp(vmfile,"<none>")!=0)
    {
      fp=fopen(vmfile,"r");                                         //read virialmassfile
      fscanf(fp,"%d %d %s %s", &checknpart, &nhalo, &snfile, &friendsfile);
      if (strcmp(snfile,evolvedfile)!=0) 
	{
	  printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	  printf("Warning: Evolved Snapshotfile and the one specified in %s are not identical!\n",vmfile);
	  printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	}
      if (tr_halo > nhalo)
	{
	  fprintf(stderr,"Illegal Halo index\n");
	  exit(2);
	}
      j=1;
      while (j < tr_halo)
	{
	  fscanf(fp,"%d %f %d %f %f %f %f %f", &j, &halomass, &halonpart, &halolmd, &halorad, &cm[0], &cm[1], &cm[2]);
	}
      fclose(fp);
    }

  if ( (!cm[0]) && (!cm[1]) && (!cm[2]) && !(cmset) )
    {
      fprintf(stderr,"Neither -m nor -cm option used\n");
      exit(3);
    }
  if ( (!cm[0]) && (!cm[1]) && (!cm[2])) 
    {
      nopb=1;
      //     tr_factor=1.0;
    }
  /*
  printf("n %d\n",j);
  printf("cm %f %f %f\n",cm[0], cm[1], cm[2]);
  printf("friendsfile %s\n",friendsfile);
  printf("snfile %s\n",snfile);
  printf("n %d\n",tr_halo);
  */

  tr_halo_cnt=0;
  masstot=0;

  if (GAP==GAP2BORDER) 
    {
      GAP = evolved.boxsize/4.;
    }

  for (k=0; k<3; k++)                                                         //check whether periodic boundaries are needed
    {
      if ( ((cm[k]<GAP) || (cm[k]>(evolved.boxsize-GAP) || pb[k])) && (!nopb))
	{
	  pb[k]=1;
	  cm[k]=MOVE(cm[k],boxsz);
	}
    }
  if (DEBUG) printf("search center \n");fflush(stdout);
  for (i=0; i<3; i++) scm[i]=cm[i];
if (srad>0)
  {
  tr_halo_id=(int *)malloc(numpart*sizeof(int));
  tr_halo_i =(int *)malloc(numpart*sizeof(int));
  maxdist=srad*srad;
  if (dmonly) {start=evolved.npart[0]; end=evolved.npart[0]+evolved.npart[1];}
  else if (starsonly) {start=evolved.npart[0]+evolved.npart[1]+evolved.npart[2]+evolved.npart[3]; end=start+evolved.npart[4];}
  else {start=0; end=numpart;}

  {start=0; end=numpart;}
  if ((end-start)<mincnt) exit(2);
  while (tr_halo_cnt<mincnt)
    {
      tr_halo_cnt=0;
      for (i=start; i < end; i++)
	{
	  dist=0;
	  for (j=0; j<3; j++) 
	    if (pb[j]) dist+=pow(MOVE(P[i].pos[j],boxsz)-cm[j],2);
	    else dist+=pow(P[i].pos[j]-cm[j],2);

	  if  ((dist<maxdist) && ( (1<<P[i].type) & use_cm) )
	    {
	      tr_halo_id[tr_halo_cnt]=P[i].id;
	      tr_halo_i[tr_halo_cnt++]=i;
	    }
	}
      maxdist*=4;
    }
  printf("Number of particles for CM calculation: %d\n",tr_halo_cnt );fflush(stdout);
  tr_halo_id=realloc(tr_halo_id,sizeof(int)*tr_halo_cnt);
  tr_halo_i =realloc(tr_halo_i ,sizeof(int)*tr_halo_cnt);
  
  }


  
  /*  
  for (i=0; i < checknpart; i++)
    {
      dum=1;
      for (k=0, k<3; k++)
	{
	  if (pb[k]) 
	    {
	      posdum=(MOVE(P[i].pos[k],boxsz));
	      if ( (posdum < (MOVE(cm[k],boxsz)-GAP)) || (posdum > (MOVE(cm[k],boxsz)+GAP )) ) dum=0;
	    } else {
	      posdum=P[i].pos[k];
	      if ( (posdum < (cm[k]-GAP)) || (posdum > (cm[k]+GAP)) ) dum=0;
	    }
	}

       if (dum) 
	{
	  tr_halo_id[tr_halo_cnt]=P[i].id;
	  tr_halo_i[tr_halo_cnt++]=i;
	}
    }

  tr_halo_id=realloc(tr_halo_id,sizeof(int)*tr_halo_cnt);
  tr_halo_i =realloc(tr_halo_i ,sizeof(int)*tr_halo_cnt);
  
  //  for (j=0; j<3; j++) printf("min %f max %f cm %f\n", min[j],max[j],cm[j]);
  //  pos_dum=(fltarr *)malloc(sizeof(fltarr)*numpart);
  tr_halo_i_dum =(int *)malloc(sizeof(int)*numpart);
  count=0; 


  for (i=0; i < numpart; i++)                                                //
    {
      dum=1;
      for (j=0; j<3; j++) 
	{
	  if (pb[j]) 
	    {
	      if (((MOVE(P[i].pos[j],boxsz)+INCLUDE) < min[j]) || ((MOVE(P[i].pos[j],boxsz)-INCLUDE) > max[j])) dum=0;
	    }
	    else {if (((P[i].pos[j]+INCLUDE) < min[j]) || ((P[i].pos[j]-INCLUDE) > max[j])) dum=0;}
	}
      if (dum) 
	{
	  //	  for (j=0; j<3; j++) pos_dum[count][j]=P[i].pos[j];
	  tr_halo_i_dum[count]=i;
 	  count++;
	}
    }
  maxdist=0;
  for (j=0; j<3; j++) 
    {
      //if (pb[j]) maxdist+=pow((MOVE(max[j],boxsz)-MOVE(min[j],boxsz))/2,2);else
       maxdist+=pow((max[j]-min[j])/2,2);
    }
  maxdist=sqrt(maxdist)*2;
  */

if (srad>0)
{
  double cvel[3]={0.,0.,0.};
  maxdist=srad;
  count=tr_halo_cnt;
do
  {
    // printf("Center of Mass for halo %d : %f %f %f particles: %d rad %f\n ", tr_halo, cm[0], cm[1], cm[2],count,maxdist);

    //pos_dum=realloc(pos_dum,sizeof(fltarr)*count);
    //tr_halo_i_dum=realloc(tr_halo_i_dum,sizeof(int)*count);

    //printf("ja\n");fflush(stdout);

  masstot=0;
  maxdist*=0.92;
  for (j=0; j<3; j++){rdum[j]=cm[j];}
  cm[0]=0;
  cm[1]=0;
  cm[2]=0;
  //  printf("%d\n",count);
  for (i=0; i < count; i++)
    {
      for (j=0; j<3; j++) 
	{
	  if (pb[j]) cm[j]+=MOVE(P[tr_halo_i[i]].pos[j],boxsz)*P[tr_halo_i[i]].mass;
	  else cm[j]+=P[tr_halo_i[i]].pos[j]*P[tr_halo_i[i]].mass;
	}
      masstot+=P[tr_halo_i[i]].mass;
      //      if (i<10)
      // printf("%f %f %f !\n", P[tr_halo_i[i]].pos[0]*P[tr_halo_i[i]].mass, cm[1], cm[2]);fflush(stdout);
      //printf("%d\n",i);fflush(stdout);
    }

  //  printf("%f %f %f !\n", cm[0], cm[1], cm[2]);fflush(stdout);

  for (j=0; j<3; j++) cm[j]=cm[j]/(masstot);
  
  dum=0;
  for (i=0; i < count; i++)                                                //
    {
      dist=0;
      for (j=0; j<3; j++) 
	{
	  if (pb[j]) dist += pow((MOVE(P[tr_halo_i[i]].pos[j],boxsz)-cm[j]),2);
	  else dist += pow((P[tr_halo_i[i]].pos[j]-cm[j]),2);
	}
      
      dist=sqrt(dist);
      if (dist < maxdist) 
	{
	  //	  for (j=0; j<3; j++) pos_dum[dum][j]=pos_dum[i][j];
	  tr_halo_i[dum]=tr_halo_i[i];
	  dum++;
	}
    }
  count=dum;
  if ((count<50) && (cvel[0]==0))
    {
      masstot=0;
      for (i=0; i< count; i++)
	{
	  for (j=0; j<3; j++) 
	    {
	      cvel[j]+=P[tr_halo_i[i]].vel[j]*P[tr_halo_i[i]].mass;
	    }
	  masstot+= P[tr_halo_i[i]].mass;
	}
      for (j=0; j<3; j++) cvel[j]=cvel[j]/(masstot);
      printf("%d particles for cvel-calculation\n",count);
      printf("Central Velocity %10.4f %10.4f %10.4f\n", cvel[0], cvel[1], cvel[2] );
    }
  // printf("count: %d  cm: %f %f %f change: %f% f% f rad: %f \n",count, cm[0],cm[1], cm[2], cm[0]-rdum[0], cm[1]-rdum[1], cm[2]-rdum[2], maxdist);
  } while (count>5);
 


 if (addbh)
   {
     if (use_cm!=16) printf("Black holed not centered on stars, are you sure? (think about using -ucm 16)\n");
     for (j=0; j<3; j++)
       {
	 P[numpart].pos[j]=cm[j];
	 P[numpart].vel[j]=cvel[j];
       }
     P[numpart].mass=P[0].mass;
     P[numpart].id=numpart+1;
     out=evolved;
     out.npart[5]++;
     out.nall[5]++;
     char outfilename[256];
     sprintf(outfilename,"%s-bh", evolvedfile);
     //     writegadget(outfilename, out, pos_ev, vel, id_ev, mass_ev);
     writegadget_part(outfilename, out, P);
   }

 }
/*
//printf("CM found %.2f now sort by distance and calculate virial radius+mass\n",step());
// free(pos_dum);
 free(tr_halo_i_dum);
   
 
part=(struct particle *)malloc(sizeof(struct particle)*tr_halo_cnt);
for (i=0; i< tr_halo_cnt; i++) 
       {
	 dist=0;
	 for (j=0; j<3; j++) 
	   {
	     if (pb[j]) dist += pow((MOVE(P[tr_halo_i[i]].pos[j],boxsz)-cm[j]),2);
	     else dist += pow((P[tr_halo_i[i]].pos[j]-cm[j]),2);
	   }
	 part[i].dist=sqrt(dist);
	 //part[i].mass=P[tr_halo_i[i]].mass;
	 //part[i].id=P[tr_halo_i[i]].id;
	 part[i].ind=tr_halo_i[i];
       } 
 qsort(&part[0],tr_halo_cnt,sizeof(struct particle),(void *)cmp);
 masstot=0;
     i=0;
     od=201;
     while (((od > 200) && (i<tr_halo_cnt)) || (i<11))                        //Add mass until overdensity drops below 200 or all friends are examined
       {
	 masstot+=P[part[i].ind].mass;
	 od=masstot/(pow(part[i].dist*evolved.time,3)*(4.0/3.0)*PI*mdens);
	 //printf("i %d od %f\n",i,od);
	 i++;
       }
     masstot-=P[part[i-1].ind].mass;
  */

  //  if (halorad!=0) maxdist=halorad; else maxdist=1000;

  if (DEBUG) printf("Determine Virial Radius \n");fflush(stdout);
  maxdist=GAP*GAP;
  effrad=30*evolved.hubparam;
  meff=0;
  sqrenvdensrad=SQR(ENVDENSRAD); 

     part=(struct particle *)malloc(sizeof(struct particle)*numpart);
     vpart=(struct particle *)malloc(sizeof(struct particle)*numpart);
     count=0;
     envdens=0;
     for (i=0; i< numpart; i++) 
       {
	 dist=0;
	 for (j=0; j<3; j++) 
	   {
	     if (pb[j]) dist += pow((MOVE(P[i].pos[j],boxsz)-cm[j]),2);
	     else dist += pow((P[i].pos[j]-cm[j]),2);
	     if (dist>maxdist) break;
	   }

	 if (dist < sqrenvdensrad)
	   {
	     envdens+=P[i].mass;
	   }

	 if (dist < maxdist) 
	   {
	     dist=sqrt(dist);
	     part[count].dist=dist;
	     part[count++].ind=i;
	     if (dist < effrad)
	       {
		 dum2=0;
		 for (m=0; m<6; m++)
		   {
		     dum2+=evolved.npart[m];
		     if (dum2>i) break;
		   }
		 if ((1<<m) & use_cm)
		   {
		     meff += P[i].mass;
		   }
	       }
	   }
	 	 
       }
     if (use_cm==1) printf("Gas mass inside 30 kpc %f\n", meff);
     if (use_cm==16) printf("Stellar mass inside 30 kpc %f\n", meff);
     if (use_cm==17) printf("Baryonic mass inside 30 kpc %f\n", meff);
     qsort(&part[0],count,sizeof(struct particle),(void *)cmp);                  //sort particles by distance to CM of halo
     
     masstot=0;
     i=0;
     od=201;
     double dm_mass=0;
     double b_mass=0;
     double gas_mass=0;
     double star_mass=0;
     while ((od > VIROD) || (i<5))                                                    //Add mass until overdensity drops below 200
       {
	 int ind = part[i].ind;
	 masstot+=P[part[i].ind].mass;
	 if ((P[ind].type == 0) ||(P[ind].type == 4)) 
	   {
	     b_mass += P[ind].mass;
	     if (P[ind].type == 0) gas_mass += P[ind].mass;
	     else star_mass += P[ind].mass;
	   } else 
	   {
	     dm_mass+= P[ind].mass;
	   }
	 od=masstot/(pow(part[i].dist*evolved.time,3)*(4.0/3.0)*PI*cdens);
	 //printf("i %d od %f\n",i,od);
	 vpart[i].dist=part[i].dist;
         vpart[i].ind=part[i].ind;
	 i++;
	 if (i > count)
	   {
	     fp=fopen("error_trhalo.dat","a");
	     fprintf(fp,"Halo %d is making trouble\n",tr_halo);
	     fclose(fp);
	     printf("halo %d is making trouble\n",tr_halo);
	     break;
	   }
	 if (i>numpart) {i--;break;}
	 if ((part[i].dist > searchdist) && (searchdist)) break;
       }
     if (DEBUG) printf("Virial Radius found count %d\n", count);fflush(stdout); 
     masstot-=P[part[i-1].ind].mass;
     if ((P[part[i-1].ind].type == 0) ||(P[part[i-1].ind].type == 4)) 
	   {
	     b_mass -= P[part[i-1].ind].mass;
	   } else 
	   {
	     dm_mass-= P[part[i-1].ind].mass;
	   }
     dum=i-1;  
     vcnt=i;
     halorad=part[vcnt-1].dist;
     envdens=envdens/((4.0/3.0)*PI*pow(ENVDENSRAD*evolved.time,3));
     vpart=realloc(vpart,sizeof(struct particle)*vcnt);
     qsort(&vpart[0],vcnt,sizeof(struct particle),(void *)cmpind);

     {
       gadpart_dist *tmppart= malloc (sizeof(gadpart_dist) * vcnt);
       for (i=0; i< vcnt; i++)
	 {
	   tmppart[i].dist=part[i].dist;
	   tmppart[i].part=P[part[i].ind];
	 }
       double par[2]={0.005,20};
       double rcs;
       double conc = nfwfit(par, tmppart, vcnt, halorad, sfl, &rcs);     
       printf("NFW: %f %f\nc %f\nRs %f\n", par[0], par[1], conc, rcs);
       free(tmppart);
     }

 for (j=0; j<3; j++) 
   {
     if (pb[j])
       {
	 cm[j]=MOVE(cm[j],boxsz);
	 printf("!!!periodic boundaries used in Dimension %d !!!\n",j);
       }
   }

 printf("\n halo %5d consists of %5d particles, virial Mass/h %5.2f\n rad %5.2f         od %4.2f\n",tr_halo,vcnt,masstot,part[i-1].dist,od);
 printf(" Baryonic Mass/h %5.2f\n DM Mass/h %5.2f \n\n",b_mass, dm_mass);
 printf("Center of Mass: %10.4f %10.4f %10.4f\n",cm[0]       ,cm[1]       ,cm[2]);
 printf("Shift:          %10.4f %10.4f %10.4f\n",cm[0]-scm[0],cm[1]-scm[1],cm[2]-scm[2]);
 if (printcm)
   {
     FILE *filep= fopen("cm.dat", "w");
     fprintf(filep, "%f %f %f\n",cm[0],cm[1],cm[2]);
     fclose(filep);
   }
 printf("particles in box: %d\n",count);
 printf("Environmental Density (R = %f ): %g\n", ENVDENSRAD,(envdens/bgdens)-1);
 if (DEBUG) printf("... \n");fflush(stdout); 
 tr_halo_id=(int *)malloc(count*sizeof(int));
 tr_halo_i =(int *)malloc(count*sizeof(int));
 tr_halo_cnt=0;

 if (docut)
   {
     gadpart *OUT= (gadpart * ) malloc (sizeof(gadpart)*vcnt);
     /*     outpos =(fltarr *)malloc(sizeof(fltarr)*vcnt);
     outvel =(fltarr *)malloc(sizeof(fltarr)*vcnt);
     outid =(int *)malloc(sizeof(int)*vcnt);
     outmass =(float *)malloc(sizeof(float)*vcnt);
     */
     out=evolved;
     for (i=0; i<6; i++)
       {
	 out.npart[i]=0;
	 out.nall[i]=0;
       }
     i=0;
     while (i<vcnt)
       {
	 dist=vpart[i].dist;
	 k=vpart[i].ind;
	 OUT[i]= P[k];
	 for (j=0; j<3; j++)
	   {
	     if (pb[j]) OUT[i].pos[j]=MOVEB(P[k].pos[j])-MOVEB(cm[j]);
	       else OUT[i].pos[j]=P[k].pos[j]-cm[j];
	     //	     OUT[i].vel[j]=P[k].vel[j];
	   }

	 /*
	 for (l=0; l<6; l++)
	   {
 	     dum+=evolved.npart[l];
	     if (k<dum) break;
	   }
	 */
	 l = P[k].type;
	 out.npart[l]++;
	 out.nall[l]++;
	 i++;
       }
     //     writegadget(gadfile, out, outpos, outvel, outid, outmass);
     if (DEBUG) printf("writing... \n");fflush(stdout); 
     writegadget_part(gadfile, out, OUT);
     free(OUT);
     /*
     free(outpos);
     free(outvel);
     free(outid);
     free(outmass);
     */
     if (DEBUG) printf("outputfile written... \n");fflush(stdout); 
   }

 maxdist=halorad*tr_factor;
 double hmeff=0;
 double galmass=0;
 double gasmass=0;
 double DMmass=0;
 double innertotm=0;
 double gsfr=0;
 double meanage=0;
 int nstars_gal=0;
 dist=0;
 i=0;
 for (j=0; j<50; j++)
   {
     p[j]=0;
     err[j]=0;
     rlog[j]=0;
   }

 masstot=0;
 while ((dist<maxdist) && (i<=numpart))
   { 
     dist=part[i].dist;
     dum2=0;
     for (m=0; m<6; m++)
       {
	 dum2+=evolved.npart[m];
	 if (dum2>part[i].ind) break;
       }
     if ((dist<halorad) && (m==1) && (dist>sfl)) 
       {
	 d=log10(dist);
	 j=floor(d/(log10(halorad)/50));
	 err[j]++;
	 p[j]+=P[part[i].ind].mass;
       }
     if  ( ((1<<m) & use_cm) && (hmeff < (meff/2.0)) )
       {
	 hmeff+=P[part[i].ind].mass;
	 effrad=dist;
       }
     if ((dist <= (halorad * 0.1)) )
       {
	 innertotm+= P[part[i].ind].mass;
	 if (m==4) 
	   {
	     galmass+= P[part[i].ind].mass;
	     double redsh = (1/P[part[i].ind].stellarage) - 1;
	     double sage = TIMEDIFF(redsh, (1/evolved.time)-1 );
	     meanage += sage;
	     nstars_gal++;
	   }
	 else if (m==0) 
	   {
	     gasmass+= P[part[i].ind].mass;
	     gsfr+= P[part[i].ind].sph->sfr;
	   }
	 else if (m==1) DMmass+= P[part[i].ind].mass;
       }
     i++;
   }
 meanage /= nstars_gal;
 double innerdens=innertotm/(pow((halorad * 0.1)*evolved.time,3)*(4.0/3.0)*PI);
 printf("hmeff: %f effrad: %f\n", hmeff, effrad);
 printf("Galaxy mass (M_star < 10 %% Rvir): %f \n", galmass);
 printf("mean age: %g \n", meanage);
 printf("SFR (< 10 %% Rvir): %g \n", gsfr);
 printf("specific SFR (< 10 %% Rvir): %g \n", gsfr / galmass);
 printf("Gas mass (M_gas < 10 %% Rvir): %f \n", gasmass);
 printf("darkmatter mass (M_gas < 10 %% Rvir): %f \n", DMmass);
 printf("Total mass < 10 %% Rvir: %f \n", innertotm);
 printf("overdensity < 10 %% Rvir: %g\n", innerdens/bgdens - 1);
 printf("density < 10 %% Rvir: %g\n", innerdens);
 if (DEBUG) printf("lambda calculation: \n");fflush(stdout); 
  nfun=0;
 k=0;

 /*fitting

 d=log10(halorad)/50;
 for (j=0; j<50; j++)
   {
     if (err[j]!=0)
       {
	 err[nfun]=(1/(sqrt(err[j])));
	 //	 err[nfun]=(1/((err[j])));
	 if (k==0) {
       rlog[nfun]=(pow(10, ((j+1)*d )))/2;
       p[nfun]=p[j]/((4.0/3.0)* PI * (pow(10, ((j+1)*d*3))));
       //       printf("j %d\n", j);

     }
       else {
	 ddum=p[j];
	 rlog[nfun]=(pow(10, ((j+1)*d )) + pow(10, (k*d)))/2;
	 p[nfun]=p[j]/((4.0/3.0)* PI * (pow(10, ((j+1)*d*3)) - pow(10, (k*d*3))));
	 if (p[nfun]<0) printf("%g %e %e \n", ddum, pow(10, ((j+1)*d*3)), pow(10, (k*d*3)) );
       }
 
//     if (rlog[nfun] < SOFTENING) 
//       {
//	 err[nfun]*=10;
//	 printf("S %d\n", nfun);
//       }
 
     k=j+1;
     nfun++;
       } else {
	 //	 p[j]=0;
	 //	 err[j]=1e10;
       }
   }
 
 for (j=0; j<nfun; j++)
   {
     //    printf("%2d: rlog %g  p %e err %e d %g d*j %g \n", j, rlog[j], p[j], err[j], pow(10,d*j) ,d*j );
          p[j]=log10(p[j]);

   }
 
 
 if (DEBUG) printf("fitting: ");fflush(stdout); 
 par[0]=1.00;
 par[1]=10;
 // lm_minimize(nfun, 2, par, err, lm_evaluate_default, printout, &data, &control);
 if (DEBUG) printf("done\n ");fflush(stdout); 
 
 printf("Delta_c %g   Scale Radius %g   Concentration factor c %g\n", par[0], par[1], halorad/par[1]);

  fitting*/

 // printf("%g %g\n", d,  pow(10,d*50));
 i--;
 if (doidl)
   {
     fp=fopen(idlfile,"w");
     fprintf(fp,"dist/mass %s %d %d %g %g %g %g %g\n",vmfile,tr_halo,i,masstot,halorad,cm[0],cm[1],cm[2]);
   }
 dist=0;
 i=0;
 while ((dist < maxdist) && (i<=numpart))
   {
     tr_halo_id[tr_halo_cnt]=P[part[i].ind].id;
     tr_halo_i[tr_halo_cnt++]=part[i].ind;
     dist=part[i].dist;
     if (doidl)
       fprintf(fp,"%g %g %d %d\n",dist, P[part[i].ind].mass, P[part[i].ind].id, P[part[i].ind].type );
     i++;
   }
 if (tr_halo_cnt> numpart) tr_halo_cnt=numpart;
 printf("number of particles within %g times virial radius %d\n trhalo_cnt %d\n", tr_factor ,i, tr_halo_cnt);
 if (doidl) fclose(fp);
 if (count!=tr_halo_cnt)
   {
     tr_halo_id=realloc(tr_halo_id,sizeof(int)*tr_halo_cnt);
     tr_halo_i =realloc(tr_halo_i ,sizeof(int)*tr_halo_cnt);
   }

 if (DEBUG) printf("... \n");fflush(stdout); 
 /**********************************************************************************************************************************************/

 // for (l=-5; l<=5; l++)
   
     count=vcnt; //+l*300;
 for (k=0; k<3; k++)
   {
     jcm[k]=0;
     vcm[k]=0;
     J[k]=0;
   }
 totj=0;
 massdum=0;

 for (j=0; j < count; j++)
   {
     for (k=0; k<3; k++)
       {
	 vcm[k]+=(P[part[j].ind].vel[k]*P[part[j].ind].mass);
	 if (pb[k]) jcm[k]+=(MOVE(P[part[j].ind].pos[k],boxsz)*P[part[j].ind].mass);
 	 else  jcm[k]+=(P[part[j].ind].pos[k]*P[part[j].ind].mass);
       }
      massdum+=P[part[j].ind].mass;
   }
 for (k=0; k<3; k++)
   {
     vcm[k]=vcm[k]/massdum;
     jcm[k]=jcm[k]/massdum;
     jcm[k]=jcm[k];
     //test
     //jcm[k]=cm[k];
   }
 if (dojr) fp=fopen(jrfile,"w");
 for (j=0; j < count; j++)
   {
     for (k=0; k<3; k++) 
       {
	 if (pb[k]) rdum[k]=MOVE(P[part[j].ind].pos[k],boxsz)-jcm[k];
	 else rdum[k]=P[part[j].ind].pos[k]-jcm[k];
	 vdum[k]=P[part[j].ind].vel[k]-vcm[k];
       }
     torq[0]=(rdum[1]*vdum[2]-rdum[2]*vdum[1]);
     torq[1]=(rdum[2]*vdum[0]-rdum[0]*vdum[2]);
     torq[2]=(rdum[0]*vdum[1]-rdum[1]*vdum[0]);

     ddum=0;
     for (k=0; k<3; k++) 
       {
	 J[k]+=torq[k]*P[part[j].ind].mass;
	 ddum+=torq[k]*torq[k];
       }
     totj+=sqrt(ddum)*P[part[j].ind].mass;
     
     ddum=0;
     ddum1=0;
     for (k=0; k<3; k++) {ddum+=SQR(rdum[k]); ddum1+=SQR(torq[k]);}
     if (dojr) fprintf(fp,"%g %g\n", sqrt(ddum), sqrt(ddum1)*P[part[j].ind].mass);
    }
 if (dojr) fclose(fp);


 lmd=0;
  for (k=0; k<3; k++) {lmd+=J[k]*J[k];}
  lmd=(sqrt(lmd)*1e10*(Msun/evolved.hubparam)*1e3*(kpc/evolved.hubparam))/(sqrt(2)*part[count-1].dist*(kpc/evolved.hubparam)*massdum*1e10*(Msun/evolved.hubparam));
 // lmd=(totj*1e10*(Msun/evolved.hubparam)*1e3*(kpc/evolved.hubparam))/(sqrt(2)*part[count-1].dist*(kpc/evolved.hubparam)*massdum*1e10*(Msun/evolved.hubparam));
 lmd=lmd/sqrt(G*massdum*1e10*(Msun/evolved.hubparam)/(part[count-1].dist*(kpc/evolved.hubparam)));

 printf("Spin: lambda %e  r %g\n",lmd, part[count-1].dist);
   

 //   fp=fopen(parfile,"w");
 //   fprintf(fp,"%g %g %g %g \n", lmd, par[0], par[1], halorad/par[1] );
 //   fclose(fp);
 /**********************************************************************************************************************************************/
   //Determine Halo-Shape

   if (use_inertia)
   {
     maxdist=rad_ratio*halorad;
     if (!(use_inertia&2)) 
       {
	 if (gal_all) maxdist=30*evolved.hubparam;
	 else maxdist=effrad;
       }
     double q=1, s=1, s_old;
     gsl_matrix *I = gsl_matrix_alloc (3, 3);
     gsl_vector *eval = gsl_vector_alloc (3);
     gsl_matrix *evec = gsl_matrix_alloc (3, 3);
     gsl_matrix *LU = gsl_matrix_alloc (3, 3);
     gsl_matrix *inv  = gsl_matrix_alloc (3, 3);     
     gsl_eigen_symmv_workspace * w =  gsl_eigen_symmv_alloc (3);
     gsl_matrix *rotation = gsl_matrix_alloc (3, 3);
     gsl_matrix_set_identity(rotation);
     gsl_matrix *resultmatrix = gsl_matrix_alloc (3, 3);
     gsl_matrix_set_zero(I);
     struct gadpart *wpart;
      
     wpart= (struct gadpart *)malloc(sizeof(struct gadpart)*vcnt);
     m=0;
     for (k=0; k < vcnt; k++)
       {
	 l=part[k].ind;
	 for (j=0; j < 3; j++)
	   {
	     if (pb[j]) wpart[k].pos[j]=MOVEB(P[l].pos[j])-MOVEB(cm[j]);
	       else wpart[k].pos[j]=P[l].pos[j]-cm[j];
	     wpart[k].vel[j]=P[l].vel[j];
	   }
	 wpart[k].mass=P[l].mass;
	 wpart[k].id=P[l].id;
	 dum=0;
	 for (m=0; m<6; m++)
	   {
 	     dum+=evolved.npart[m];
	     if (l<dum) break;
	   }
	 wpart[k].type=m;
       }    
     
     
     do
     {
     s_old=s;
     gsl_matrix_set_zero(I);
     for (k=0; k < vcnt; k++)
       {     
	 l=part[k].ind;
	 for (i=0; i < 3; i++)
	   for (j=i; j < 3; j++)
	     {
	       ddum =wpart[k].pos[i] * wpart[k].pos[j];
	       dist =SQR(wpart[k].pos[0])+SQR(wpart[k].pos[1]/q)+SQR(wpart[k].pos[2]/s);
	       ddum/= dist;
	       if ((sqrt(dist)<maxdist) && ((1<<wpart[k].type)&use_inertia)) 
		 {
		   gsl_matrix_set(I,i,j, gsl_matrix_get(I,i,j)+ddum);
		   if (i!=j) gsl_matrix_set(I,j,i, gsl_matrix_get(I,j,i)+ddum);
		 }
	     }

       }

     gsl_eigen_symmv (I, eval, evec, w);
     gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
     gsl_matrix_memcpy(LU, evec);
     /*
       for (i=0; i < 3; i++)
       {
       for (j=0; j < 3; j++)
       printf("%15.4g", gsl_matrix_get(I,i,j));
       printf("\n");
       }
     */
   
     
   for (i=0; i < 3; i++)
     {
       double eval_i  = gsl_vector_get (eval, i);
       //             gsl_vector_view evec_i  = gsl_matrix_column (evec, i);
       if (i==0) ddum = sqrt(eval_i);
       else if (i==1) q=sqrt(eval_i)/ddum;
       else           s=sqrt(eval_i)/ddum;
       //printf ("Ratio = %g\n", sqrt(eval_i)/ddum);
	     //             printf ("eigenvector = \n");
	     //         gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
     }
    
   gsl_permutation *perm = gsl_permutation_alloc (3);
   int sign;

   gsl_linalg_LU_decomp (LU, perm, &sign);
   gsl_linalg_LU_invert (LU, perm, inv);
  

   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
		   1.0, inv, rotation,
		   0.0, resultmatrix);
   gsl_matrix_memcpy (rotation, resultmatrix);

   //Rotate Particles 
  gsl_vector *oldpos = gsl_vector_alloc (3);
  gsl_vector *newpos = gsl_vector_alloc (3);  
  gsl_vector *oldvel = gsl_vector_alloc (3);
  gsl_vector *newvel = gsl_vector_alloc (3);
  ddum=0;
  for (k=0; k < vcnt; k++)
    {
      for (i=0; i<3; i++)
	{
	  gsl_vector_set(oldpos, i, wpart[k].pos[i]);
	  gsl_vector_set(oldvel, i, wpart[k].vel[i]);
	}
      gsl_blas_dgemv( CblasNoTrans, 1.0, inv, oldpos, 0.0, newpos);
      gsl_blas_dgemv( CblasNoTrans, 1.0, inv, oldvel, 0.0, newvel);
      //ddum += gsl_blas_dnrm2 (oldvel) - gsl_blas_dnrm2 (newvel);
      for (i=0; i<3; i++)
	{
	  wpart[k].pos[i]=gsl_vector_get(newpos, i);
	  wpart[k].vel[i]=gsl_vector_get(newvel, i);
	}
    }
  //printf("%g\n", ddum);
     } while ((ABS(s_old-s)/s) > 1e-2);

     if (docut)
       {
     sprintf( gadfile, "%s_rot", hname);
     FILE *matrixf=fopen(gadfile,"w");
     gsl_matrix_fwrite (matrixf, rotation);
     fclose(matrixf);
       }

     gsl_matrix_set_zero(I);
     for (k=0; k < vcnt; k++)
       {     
	 l=part[k].ind;
	 for (i=0; i < 3; i++)
	   for (j=i; j < 3; j++)
	     {
	       ddum =wpart[k].pos[i] * wpart[k].pos[j];
	       dist =SQR(wpart[k].pos[0])+SQR(wpart[k].pos[1]/q)+SQR(wpart[k].pos[2]/s);
	       ddum/= dist;
	       if ((sqrt(dist)<maxdist) && ((1<<wpart[k].type)&use_inertia)) 
		 {
		   gsl_matrix_set(I,i,j, gsl_matrix_get(I,i,j)+ddum);
		   if (i!=j) gsl_matrix_set(I,j,i, gsl_matrix_get(I,j,i)+ddum);
		 }
	     }
	 
       }
     printf("I: ");
     for (i=0; i < 3; i++) printf("%g ",gsl_matrix_get(I,i,i));
     printf("\n");

     printf("longest axis = %g r_ratio: %g \nq = %g\ns = %g\n", maxdist, rad_ratio, q, s);
     fflush(stdout);
/*************************************************************************************/     
/*Calculate Velocity anisotropy                                                      */

     gadpart      * vapart=malloc(vcnt* sizeof(gadpart));
     gadpart_dist * dpart =malloc(vcnt* sizeof(gadpart_dist));
     if ((vapart==NULL) || (dpart==NULL)) {fprintf(stderr, "Error allocating memory!\n");exit(1);}
     double Pi[3]={0.0, 0.0, 0.0};
     double T[3]={0.0, 0.0, 0.0};
     double Pi_PHI=0;
     double  T_PHI=0;
     double Pi_RR=0;
     double  T_RR=0;
     double delta;
     int num_va=0;
     for (i=0; i< vcnt; i++)
       {
	 if ((1<<(wpart[i].type))&use_va)
	   {
	     cpygadpart(&vapart[num_va]    , &wpart[i]);
	     cpygadpart(&dpart[num_va].part, &wpart[i]);
	     dpart[num_va].dist=sqrt(SQR(wpart[i].pos[0])+SQR(wpart[i].pos[1]/q)+SQR(wpart[i].pos[2]/s));
	     num_va++;
	   }
       }
     if (DEBUG) {printf("Begin VA Calculation %d %d\n", vcnt,  num_va);fflush(stdout);}
     if (num_va > 50)
       {

	 qsort(&dpart[0], num_va, sizeof(gadpart_dist), (void *) cmp_dist);
	 if (DEBUG) {printf("particles sorted \n");fflush(stdout);}
	 int distnum= gadsearch(dpart, maxdist, 0, num_va);
	 if (DEBUG) {printf("maxdist determined \n");fflush(stdout);}
	 //	 dpart= realloc(dpart, distnum*sizeof(gadpart_dist));

	 KdNode * root;
	 initKdNode(&root, NULL);
	 buildKdTree(root, vapart, num_va, 0);
	 if (DEBUG) {printf("KD-Tree completed \n");fflush(stdout);}

	 for (i=0; i<3; i++) Pi[j]=0;
	 int incstep;
	 if (distnum<120) incstep=1;
	 else if (distnum<300) incstep=5;
	 else if (distnum<500) incstep=10;
	 else if (distnum<5000) incstep=30;
	 else incstep=50;

	 for (i=0; i<distnum; i+=incstep)
	   {
	     masstot=0;
	     gadpart_dist * knn;
	     double knndist=findkNN(root, &dpart[i].part, 15, &knn, kNN);
	     double meanv[3]={0.0, 0.0, 0.0};
	     double sqrmv[3]={0.0, 0.0, 0.0};
	     double R   = sqrt(SQR(dpart[i].part.pos[0])+SQR(dpart[i].part.pos[1]));
	     double ang;
	     double vR,vR2, vRsum;
	     double vPHI, vPHI2, vPHIsum;
	     for (j=0; j<3; j++)
	       {
		 meanv[j]=dpart[i].part.vel[j]*dpart[i].part.mass;
		 sqrmv[j]=SQR(dpart[i].part.vel[j])*dpart[i].part.mass;
	       }
	     ang    = atan2(dpart[i].part.pos[1], dpart[i].part.pos[0]);
	     vR     = dpart[i].part.vel[0]* cos(ang)+ dpart[i].part.vel[1]* sin(ang);
	     vR2    = SQR(vR)*dpart[i].part.mass;
	     vRsum  = vR*dpart[i].part.mass;
	     vPHI   = -dpart[i].part.vel[0]*sin(ang)+dpart[i].part.vel[1]*cos(ang);
	     vPHI2  = SQR(vPHI)*dpart[i].part.mass;
	     vPHIsum= vPHI*dpart[i].part.mass;

	     masstot=dpart[i].part.mass;
	     for (k=0; k<kNN; k++)
	       {
		 for (j=0; j<3; j++)
		   {
		     meanv[j]+=knn[k].part.vel[j]*knn[k].part.mass;
		     sqrmv[j]+=SQR(knn[k].part.vel[j])*knn[k].part.mass;
		   }
		 //		 ang     = atan2(knn[k].part.pos[1], knn[k].part.pos[1]);
		 ang     = atan2(knn[k].part.pos[1], knn[k].part.pos[0]);
		 vR      = knn[k].part.vel[0]* cos(ang)+ knn[k].part.vel[1]* sin(ang);
		 vR2    += SQR(vR)*knn[k].part.mass;
		 vRsum  += vR*knn[k].part.mass;
		 vPHI    = -knn[k].part.vel[0]*sin(ang)+ knn[k].part.vel[1]* cos(ang);
		 vPHI2  += SQR(vPHI)*knn[k].part.mass;
		 vPHIsum+= vPHI*knn[k].part.mass;
		 masstot+=knn[k].part.mass;
	       }
	     for (j=0; j<3; j++)
	       {
		 meanv[j]=meanv[j]/masstot;
		 sqrmv[j]=sqrmv[j]/masstot;
		 double sigma=(sqrmv[j]-SQR(meanv[j]));
		 Pi[j]+=(sqrmv[j]-SQR(meanv[j]));
		 T[j]+=SQR(meanv[j]);
		 if (sigma < 0) printf("!sigmaerror!%g %g %g %g %g\n",sqrmv[j], meanv[j], SQR(meanv[j]), sigma, masstot);
	       }
	     vR2    = vR2/masstot;
	     vRsum  = vRsum/masstot;
	     vPHI2  = vPHI2/masstot;
	     vPHIsum= vPHIsum/masstot;
	     Pi_RR += vR2-SQR(vRsum);
	      T_RR += SQR(vRsum);
	     Pi_PHI+= vPHI2-SQR(vPHIsum);
	      T_PHI+= SQR(vPHIsum);
	   
	     //      printf("%5d %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n", i, dpart[i].dist, knndist, dpart[i].part.pos[0], dpart[i].part.pos[1], dpart[i].part.pos[2], Pi[i][0], Pi[i][1], Pi[i][2]);
	     free(knn);
	   }
	 
	 fprintf(stdout, "PI %g %g %g\n", Pi[0], Pi[1], Pi[2]);
	 double mPi=(Pi[0]+Pi[1])/2;
	 delta= (mPi-Pi[2])/mPi;
	 fprintf(stdout, "Pi_RR  %g\n", Pi_RR);
	 fprintf(stdout, "Pi_PHI %g\n", Pi_PHI);
	 fprintf(stdout, "delta  %g\n", delta);

	 fprintf(stdout, "T  %g %g %g\n",  T[0],  T[1],  T[2]);
	 fprintf(stdout, "T_RR  %g\n", T_RR);
	 fprintf(stdout, "T_PHI %g\n", T_PHI);
	 
	 free(root);
       }
     free(vapart);
     free(dpart);

     if (DEBUG) {printf("End of VA Calculation\n");fflush(stdout);} 
/*END VA Calculation************************************************************************************/ 

/*
  if (docut)
   {
     outpos =(fltarr *)malloc(sizeof(fltarr)*vcnt);
     outvel =(fltarr *)malloc(sizeof(fltarr)*vcnt);
     outid =(int *)malloc(sizeof(int)*vcnt);
     outmass =(float *)malloc(sizeof(float)*vcnt);
     out=evolved;
     qsort(&wpart[0],vcnt,sizeof(struct gadpart),cmp_type);
     for (i=0; i<6; i++)
       {
	 out.npart[i]=0;
	 out.nall[i]=0;
       }
     i=0;
     while (i<vcnt)
       {
	 for (j=0; j<3; j++)
	   {
	     outpos[i][j]=wpart[i].pos[j];
	     outvel[i][j]=wpart[i].P[j].vel;
	   }
	 outid[i]=wpart[i].id;
	 outmass[i]=wpart[i].mass;
	 l=wpart[i].type;
	 out.npart[l]++;
	 out.nall[l]++;
	 i++;
       }
     sprintf( gadfile, "%s_rot.gad", hname);
     writegadget(gadfile, out, outpos, outvel, outid, outmass);
     free(outpos);
     free(outvel);
     free(outid);
     free(outmass);
   }
*/
  
  gsl_matrix_free(I);
  gsl_matrix_free(evec);
  gsl_matrix_free(LU);
  gsl_matrix_free(inv);
  gsl_matrix_free(rotation);
  gsl_matrix_free(resultmatrix);
  gsl_vector_free(eval);
  gsl_eigen_symmv_free(w);
  free(wpart);

   }
   
    }else
    {
      fp=fopen(idfile, "r");
      fread(&tr_halo_cnt, sizeof(int), 1,fp);
      fread(cm, sizeof(double), 3,fp);
      tr_halo_id=(int*) malloc(tr_halo_cnt*sizeof(int));
      fread(tr_halo_id, sizeof(int), tr_halo_cnt,fp);
      fclose(fp);
      printf("TraceIDfile (%d particles) read, CM at %f %f %f\n", tr_halo_cnt, cm[0], cm[1], cm[2]);
      fflush(stdout);           
    }
  

 /**********************************************************************************************************************************************/
 if (notrace) 
   {
     if (DEBUG)
     printf("no IC-file specified, exitting...\n");
     exit(0);
   } 

 i=0;
 numpart=0;
                                                                                           //read in IC File
  fp=fopen(icfile,"r");

  fread(&blocksize,sizeof(int),1,fp);
  fread(&ic,sizeof(struct header),1,fp);                        //read header of ic gadget file
  fread(&blocksize,sizeof(int),1,fp);
  boxsz=ic.boxsize;
  
  for (i=0; i<6; i++) numpart+=ic.npart[i];
  pos_ic =(fltarr *)malloc(sizeof(fltarr)*numpart);
#ifdef LONGIDS
  long *id_ic;
  id_ic  =(long *)   malloc(sizeof(long)*numpart);
#else
  int *id_ic;
  id_ic  =(int *)   malloc(sizeof(int)*numpart);
#endif //LONGIDS

  //  mass_ic=(float *) malloc(sizeof(float)*numpart);

  fread(&blocksize,sizeof(int),1,fp);                                //read particle positions of ic snapshot file
  fread(&pos_ic[0],sizeof(fltarr),numpart,fp);
  fread(&blocksize,sizeof(int),1,fp);  
  
  fread(&blocksize,sizeof(int),1,fp);                                //skip velocities
  fseek(fp,blocksize,SEEK_CUR);
  fread(&blocksize,sizeof(int),1,fp);

  fread(&blocksize,sizeof(int),1,fp);                                //read particle ids of ic snapshot file
#ifdef LONGIDS
  fread(&id_ic[0],sizeof(long),numpart,fp);
#else
  fread(&id_ic[0],sizeof(int),numpart,fp);
#endif //LONGIDS
  fread(&blocksize,sizeof(int),1,fp);
  
//  fread(&blocksize,sizeof(int),1,fp);                                //read particle masses of ic snapshot file
//  fread(&mass_ic[0],sizeof(float),numpart,fp);
//  fread(&blocksize,sizeof(int),1,fp);
  fclose(fp);
  
  for (k=0; k<3; k++)                                                         //check whether periodic boundaries are needed
    {
      if ( ((cm[k]<GAP) || (cm[k]>(ic.boxsize-GAP) || pb[k])) && (!nopb))
	{
	  pb[k]=1;
	  //	      cm[k]=MOVE(cm[k],boxsz);
	}
    }

  //  j=P[part[0].ind].id-1;
  j=tr_halo_id[0]-1;
  coord512(ca,j);
  for (k=0; k<3; k++)
    {
      if (pb[k]) {icmin[k]=MV(ca[k], DIM)*dx512; icmax[k]=MV(ca[k], DIM)*dx512;}
      else {icmin[k]=ca[k]*dx512;icmax[k]=ca[k]*dx512;}
    }
  if (!traceid) fp=fopen("ind.dat","w");
  for (k=0; k<3; k++) {indmax[k]=j; indmin[k]=j;}
  for (i=1; i<tr_halo_cnt; i++)
    {
      //      j=P[part[i].ind].id-1;
      j=tr_halo_id[i]-1;
      coord512(ca, j);
      for (k=0; k<3; k++)
	{
	  if (pb[k])
	    {
	      ddum=MV(ca[k], DIM)*dx512;
	      if (ddum < icmin[k]) {icmin[k]= round(ddum); indmin[k]=j;}
	      if (ddum > icmax[k]) {icmax[k]= round(ddum); indmax[k]=j;}
	    } else
	    {
	      ddum=ca[k]*dx512;
	      if (ddum < icmin[k]) {icmin[k]= round(ddum); indmin[k]=j;}
	      if (ddum > icmax[k]) {icmax[k]= round(ddum); indmax[k]=j;}
	    }  
	}

      if (!traceid)
      fprintf(fp,"%g %g %g\n",P[tr_halo_i[i]].pos[0]-cm[0],P[tr_halo_i[i]].pos[1]-cm[1],P[tr_halo_i[i]].pos[2]-cm[2]);
      
    }
  if (!traceid) fclose(fp);

  
  alignf=2;
  gridsize=72000.0/(64.0*alignf);

  for (k=0; k<3; k++)
	{
	  min[k]=icmin[k];
	  max[k]=icmax[k];
	}

 if (DEBUG) 
   {
     printf("indmin: %d %d %d...\n",indmin[0], indmin[1], indmin[2]);fflush(stdout);
     printf("indmax: %d %d %d...\n",indmax[0], indmax[1], indmax[2]);fflush(stdout);
     printf("icmin: %f %f %f...\n",icmin[0], icmin[1], icmin[2]);fflush(stdout);
     printf("icmax: %f %f %f...\n",icmax[0], icmax[1], icmax[2]);fflush(stdout);

}

  if (strcmp(gridfile,"<none>"))              //build a grid and determine which cells contain particles of the halo
  {
  for (k=0; k<3; k++)
    {
      printf("%d\n", pb[k]);
      coord512(ca, indmin[k]);
      minbox[k]= (ca[k]>0) ? (ca[k]-1) : (ca[k]+DIM-1);
      coord512(ca, indmax[k]);
      maxbox[k]= (ca[k]<DIM) ? (ca[k]+1) : (ca[k]-DIM+1);
      if (pb[k]) szbox[k] =ceil((icmax[k]-MOVE(minbox[k]*dx512,boxsz))/(gridsize))+2;
	else szbox[k] =ceil((icmax[k]-(minbox[k]*dx512))/(gridsize))+2;
      size[k]  = (maxbox[k] > minbox[k]) ? (maxbox[k]-minbox[k]) : (maxbox[k]-minbox[k]+DIM);
      
      //if ((szbox[k]&1)==1) szbox[k]++;
      //      szbox[k]*=alignf;
      printf("%d %d %d \n",minbox[k],maxbox[k],szbox[k]);fflush(stdout);
    }

  grid=  (int ***) malloc (szbox[0]*sizeof(int **));                       //allocate 3-dimensional array

  for (i=0; i < szbox[0]; i++) 
    {
      grid[i]=  (int **) malloc (szbox[1] * sizeof(int *));
    }
  for (i=0; i < szbox[0]; i++) 
    for (j=0; j < szbox[1]; j++)
      {
	grid[i][j]=   (int *) calloc (szbox[2],sizeof(int ));
      }
  fp = fopen(datfile, "w");
  for (i=0; i<tr_halo_cnt; i++)
    {
      //      j=P[part[i].ind].id-1;
      j=tr_halo_id[i]-1;
      coord512(ca, j);
      for (k=0; k<3; k++)
	{
	  if (pb[k])
	    {
	      boxind[k]=floor(((MV(ca[k],512)*dx512)- MOVE(minbox[k]*dx512, boxsz))/gridsize );	
	    } else
	    {
	      boxind[k]=floor((ca[k]*dx512/gridsize)- (minbox[k]*alignf/8.0));	
	    }  
            fprintf(fp,"%g\t",(ca[k]*dx512/gridsize) );
	  if ((j==indmax[k]) ||(j==indmin[k]))  
	    {
	      printf("%d %g %g\n", k, boxind[k]*gridsize , minbox[k]*(gridsize*alignf/8.0));
	      printf("%g %d %d\n", pos_ic[j][k], boxind[k], minbox[k]);fflush(stdout);
	    }
	}
      //      if (DEBUG) {printf("%d %d %d...\n",boxind[0], boxind[1], boxind[2]);fflush(stdout);}
      grid[boxind[0]][boxind[1]][boxind[2]]++;
      fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\n", boxind[0], boxind[1], boxind[2] , ca[0], ca[1], ca[2]);
      //      if ((i<100) && (grid[boxind[0]][boxind[1]][boxind[2]]<100000)) grid[boxind[0]][boxind[1]][boxind[2]]+=100000;
    }
    fclose(fp);
    printf("%g\t%g\t%g\n!\n", dx512, gridsize, dx512/gridsize);
  if (DEBUG) {printf("grid built...\n");fflush(stdout);}
  fp=fopen(gridfile,"w");
  blocksize=28;
  fwrite(&blocksize,sizeof(int),1,fp);
  fwrite(&szbox[0],sizeof(int),3,fp);
  fwrite(&minbox[0],sizeof(int),3,fp);
  fwrite(&alignf,sizeof(int),1,fp);
  fwrite(&blocksize,sizeof(int),1,fp);

  blocksize=sizeof(int)*szbox[0]*szbox[1]*szbox[2];
  fwrite(&blocksize,sizeof(int),1,fp);
  for (k=0; k < szbox[2]; k++)
    for (j=0; j < szbox[1]; j++)
	for (i=0; i < szbox[0]; i++) 
	fwrite(&grid[i][j][k],sizeof(int),1,fp);
  fwrite(&blocksize,sizeof(int),1,fp);
  fclose(fp);
  }





  printf("Min/Max positions in IC file:\n");fflush(stdout);
  for (k=0; k<3; k++)
    {
      if (pb[k]) 
	{
	  min[k]=MOVE(min[k],boxsz);max[k]=MOVE(max[k],boxsz);
	  icmin[k]=MOVE(icmin[k], boxsz);icmax[k]=MOVE(icmax[k], boxsz);
	}

      printf("%d:code units: Min %5.2f Max %5.2f physical units: Min %5.2f Max %5.2f\n",k, min[k],max[k],min[k]/ic.hubparam,max[k]/ic.hubparam);
    }
  for (k=0; k<3; k++)
    {
      //printf("%d min %g max %g\n", k, min[k], max[k]);
    }
  printf("minx=%g;\n",icmin[0]);
  printf("maxx=%g;\n",icmax[0]);
  printf("miny=%g;\n",icmin[1]);
  printf("maxy=%g;\n",icmax[1]);
  printf("minz=%g;\n",icmin[2]);
  printf("maxz=%g;\n",icmax[2]); 
  printf("offx=%d;\n",minbox[0]);
  printf("offy=%d;\n",minbox[1]);  
  printf("offz=%d;\n",minbox[2]);
  printf("sizex=%d;\n",size[0]);
  printf("sizey=%d;\n",size[1]);
  printf("sizez=%d;\n",size[2]);
  
  printf("Expected CM in resimulation: ");
  for (k=0; k<3; k++)
    {
      diff[k]=icmax[k]-(minbox[k]*dx512);
      diff[k]=(diff[k]>0)? diff[k] : (diff[k]+boxsz);
      shift[k]=( (boxsz/2.0) - (icmin[k]+(diff[k]/2.0)));
      expcm[k]=(cm[k]+shift[k]) > 0 ? (cm[k]+shift[k]) : (cm[k]+shift[k])+boxsz;
      printf(" %g ", expcm[k]);
    }
  printf("\n");





  free(part);
  printf("Time needed: %.2f sec\n",((float)clock())/CLOCKS_PER_SEC);
  return 0;
}
