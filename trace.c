/*
Program to trace back a refined halo through a set of snapshots
(contains code overhead, since it is based on a copy of trhalomult.c)

icc -lm -lgsl -lgslcblas -o ../bin/altix/trace trace.c lmfit-2.2-icc/lmmin.o lmfit-2.2-icc/lm_eval.o libgad-icc.o KdTree-icc.o -I/usr/local/gsl-1.6

icc -lm -lgsl -lgslcblas -o ../bin/altix/trace trace.c -lgad-altix-icc KdTree-icc.o -I/usr/local/gsl-1.6
icc -lm -lgsl -lgslcblas -o ../bin/altix/trace_nopot trace.c -lgad-altix-icc-nopot KdTree-icc.o -I/usr/local/gsl-1.6 -DNOPOT

adhara:
icc -g -O2 -lm -lgsl -lgslcblas -o ../bin/other64/trace trace.c -lgad-adhara-icc adhara/KdTree-icc.o

stan:
gcc -lm -lgsl -lgslcblas -o ~/bin/trace trace.c -lgad-stan stan/KdTree.o -L/home/albireo/oser/libraries/stan/GSL/lib -L./lib/

gcc -lm -lgsl -lgslcblas -o ~/bin/trace-winds trace.c -lgad-stan-winds stan/KdTree-winds.o -L/home/albireo/oser/libraries/stan/GSL/lib -L./lib/ -DWINDS

gcc -fopenmp -lm -lgsl -lgslcblas -o ~/bin/trace trace.c -lgad-stan stan/KdTree.o -L/home/albireo/oser/libraries/stan/GSL/lib -L./lib/ stan/ompfuncs.o

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
//#include "./lmfit-2.2/lmmin.h"
//#include "./lmfit-2.2/lm_eval.h"
#include "./lmmin.h"
#include "./lm_eval.h"
#include "libgad.h"
#include "KdTree.h"

#define	PI 3.14159265358979323846
#define GAP2BORDER 10000                                     //min distance to boarders to ignore periodic boundaries
#define INCLUDE 0                                    //(friends + INCLUDE*kpc/h) are used to determine CM
#define NBOX 4 
#define TRACE_FACTOR 2.0                             //TRACE_FACTOR times virial radius will be traced back to IC
#define h0 0.72
#define DIM 512
#define ENVDENSRAD (4000 * h0)
#define INERTIA_START_RAD 0.4
#define INERTIA_USE_PART_TYPE 16
#define CM_USE_PART_TYPE 16
#define VA_USE_PART_TYPE 16
#define MIN_PARTICLE_COUNT 30
#define MAX_STAR_DIST (30*h0)
#define GRIDSIZE 1000
#define WORKSPACE 0.25                               //workspace for WORKSPACE*numpart particles is allocated, if glibc or segmentation fault occur, try to compile with a larger fraction
#define LARGE_WS 8
#define ADDSTAR_BUF 1000000
#define ADDGAS_BUF 1000000
#define kNN 50
#define	MIN(a, b)   ((a)<(b)?(a):(b))
#define	MAX(a, b)   ((a)>(b)?(a):(b))
#define ABS(a) ((a) >= 0 ? (a) : -(a))
#define CMP(a,b) ((a)>(b)?(1):(-1))
#define PB(a,b) ((a)>(b)?((a)-(b)):(a))
#define MOVE(a,b) PB((a)+((b)/2),(b))
#define MOVEB(a) MOVE((a),boxsz)
#define MV(a,b) ((a)+(b)/2)%(b)

#define SOFTENING 0.40
#define G 6.6742e-11
#define Msun 1.989e30
#define kpc 3.08567758128e19
#define MAXNHALOES 30000
#define MAXINT 10000000
#define NUMTIMESTEPS 94
#define NUMTRACE 50
#define ACC_FRAC 0.1
#define DENS_THRESH 1e-3                //Threshold for maxtemp of gas particles to be set

double cdens;

struct particle {int ind;float dist;};
clock_t t[2];

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
#ifdef POTENTIAL
  float pot;
#endif
};

struct gaspart{
  int id;
  float maxtemp;
  float Tvir;
  float Tvir_acc;
  float a_maxtemp;
  float frac_maxtemp;
  float a_acc;
  float T_acc;
  float a_star;
};

int cmp_star_id(const void *a, const void *b)
{
  struct star *x= (struct star*)a;
  struct star *y= (struct star*)b;
  if (x->id > y->id) return 1;
  if (x->id < y->id) return -1;
  return 0;
}

int cmp_gas_id(const void *a, const void *b)
{
  struct gaspart *x= (struct gaspart*)a;
  struct gaspart *y= (struct gaspart*)b;
  if (x->id > y->id) return 1;
  if (x->id < y->id) return -1;
  return 0;
}

double fitfct(double t, double *p)
{
  return log10((p[0]) / ((t/p[1]) * SQR(1+t/p[1])));
}

void dumprintout ( int n_par, double* par, int m_dat, double* fvec, 
                       void *data, int iflag, int iter, int nfev )
         {
           // dummy function to catch fitting output

         }


void usage()
{
  fprintf(stderr,"Trace v0.02\n");
  fprintf(stderr,"-f <snapshot basefilename>\n");
  fprintf(stderr,"-i <startindex> <endindex>\n");
  fprintf(stderr,"if you want to search for CM in an area around the given coordinates\n");
  fprintf(stderr,"add -sr <radius> to include particles for the search with distance to cm < radius\n");
  fprintf(stderr,"-n  <name of the halo>\n");
  fprintf(stderr,"-cm  <center of mass X Y Z>\n");
  fprintf(stderr,"-ui <bitcode for particle types to use for computation of Inertia Tensor>\n");
  fprintf(stderr,"-uva <bitcode for particle types to use for computation of Velocity Anisotropy>\n");
  fprintf(stderr,"-ucm <bitcode for particle types to use for computation of Center of Mass>\n");
  fprintf(stderr,"-cut <create Gadget-file for every halo>\n");
  fprintf(stderr,"-tf <trace factor>\n");
  fprintf(stderr,"-mincnt <minimum number of particles for CM-Search>\n");
  fprintf(stderr,"-ntr <Number of particles to trace back>\n");
  fprintf(stderr,"-cont <do not owerwrite existing files>\n");
  exit(1);
}

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
int cmp (struct particle *a, struct particle *b)
{
  return CMP(a[0].dist,b[0].dist);
}

int cmpind (struct particle *a, struct particle *b)
{
  return CMP(a[0].ind,b[0].ind);
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
  gadpart *part_ev;
  struct star *stars;
  struct gaspart *gas;
  struct header ic, evolved, out;
  FILE *fp, *outf;
  char cmfile[80],icfile[80],snapfile[256], basename[256],snfile[80],friendsfile[80],idlfile[80],asciifile[80],gridfile[80],jrfile[80], outfile[80], gadfile[80], starfile[80], prefix[20];
  unsigned int numpart=0, nummass=0;
  int blocksize,i,j,k,l,m,n,h,dum,tr_halo=0,tr_halo_cnt = 0,id,halo,checknpart,mindum,nhalo,count,notrace=0;
  int *tr_halo_id,*tr_halo_i, ca[3], size[3];
  fltarr *pos_ev,*pos_ic,*vel, *outpos, *outvel;
  float *mass_ev, *mass_dum,*mass_ic,dist,maxdist,*distance,halomass,halorad, halolmd, *outmass;
  double posdum,max[3],min[3],srad=25, icmin[3], icmax[3], envdens, bgdens;
  float lmd,vcm[3]={0,0,0},jcm[3]={0,0,0},J[3]={0,0,0}, torq[3],rdum[3], vdum[3], massdum;
  int *id_ev,*id_ic,*iclus, halonpart, *outid;
  float tm,linkle,od, tr_factor, maxstardist=0;
  double boxsz,cm[3]={0,0,0}, masstot, gridsize, rlog[50], p[50], err[50], d, rad_ratio;
  double cvel[3];
  short pb[3]={0,0,0};
  int minbox[3], maxbox[3], szbox[3], ****grid, ***gridsz, boxind[3], alignf;
  int nfun, npar, dojr=0, vcnt=0, use_inertia, use_va, use_cm;
  double par[2],ddum, dx512, ddum1, totj, rotvel, peakrotvel, peakr, vmass, sqrenvdensrad, sfl=SOFTENING;
  int  cutgad=0, idum, donotoverwrite=0, test=0, *cid, numalloc, numtrace;
  fltarr *hcm;
  int debug=0, maxhaloes, minpartcnt;
  int startind, endind, snapind, starcnt=0, gascnt=0, totalstarcnt=0, totalgascnt=0,  dostars=1;
  int addstar[ADDSTAR_BUF];
  int addgas[ADDGAS_BUF];
  double meff, effrad;
  double ll = 0.8;
  int doFOF = 1;
  int doascii = 0;
  int singlefile=0;
  double GAP=GAP2BORDER;
  double conv_dist=1.0;
  double conv_mass=1.0;
 
  //  float timestep[NUMTIMESTEPS];
 
  t[0]=clock();
  dx512=72000.0/512.0;
  
  rad_ratio=INERTIA_START_RAD;
  use_inertia=INERTIA_USE_PART_TYPE;
  use_va=VA_USE_PART_TYPE;
  use_cm=CM_USE_PART_TYPE;
  gridsize=GRIDSIZE;
  maxhaloes=MAXNHALOES;
  minpartcnt=MIN_PARTICLE_COUNT;
  numtrace=NUMTRACE;
  maxstardist=MAX_STAR_DIST;
  sprintf(prefix,"m");

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //filenames are ignored if given as command line parameters 

  strcpy(gridfile,"<none>");
  strcpy(cmfile,"<none>");
  strcpy(basename,"");
  strcpy(icfile,"");
  strcpy(idlfile,"idlinp.dat");
  strcpy(asciifile,"ascii.dat");
  tr_factor=TRACE_FACTOR;
  //  tr_halo=123;
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  i=1;
  if (argc==1) usage();
  while (i<argc)
    {
      
      if (!strcmp(argv[i],"-f"))
	{
	  i++;
	  strcpy(basename,argv[i]);
	  i++;
	}
      else if (*argv[i]!='-')
        {
          startind=endind=1;
          singlefile=1;
          strcpy(basename,argv[i]);
          i++;
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
      else if (!strcmp(argv[i],"-n"))
	{
	  i++;
	  strcpy(prefix,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-debug"))
	{
	  i++;
	  debug=1;
	}
      else if (!strcmp(argv[i],"-test"))
	{
	  i++;
	  test=1;
	}
      else if (!strcmp(argv[i],"-cont"))
	{
	  i++;
	  donotoverwrite=1;
	}
      else if (!strcmp(argv[i],"-gap"))
	{
	  i++;
	  GAP=atof(argv[i]);
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
      else if (!strcmp(argv[i],"-nofof"))
	{
	  i++;
	  doFOF=0;
	}
      else if (!strcmp(argv[i],"-pfof"))
	{
	  i++;
	  doFOF= doFOF | 2;
	}
      else if (!strcmp(argv[i],"-cutfof"))
	{
	  i++;
	  doFOF= doFOF | 4;   //not yet implemented
	}
      else if (!strcmp(argv[i],"-idl"))
	{
	  i++;
	  strcpy(idlfile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-ascii"))
	{
	  i++;
	  doascii = 1;
	  strcpy(asciifile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-da"))
	{
	  i++;
	  doascii = 1;
	}
      else if (!strcmp(argv[i],"-gr"))
	{
	  i++;
	  strcpy(gridfile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-mincnt"))
	{
	  i++;
	  minpartcnt=atoi(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-ntr"))
	{
	  i++;
	  numtrace=atoi(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-sr"))
	{
	  i++;
	  srad=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-gs"))
	{
	  i++;
	  gridsize=atof(argv[i]);
	  i++;
	}
       else if (!strcmp(argv[i],"-r"))
	{
	  i++;
	  rad_ratio=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-soft"))
	{
	  i++;
	  sfl=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-ll"))
	{
	  i++;
	  ll=atof(argv[i]);
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
      else if (!strcmp(argv[i],"-jr"))
	{
	  i++;
	  dojr=1;
	  strcpy(jrfile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-cut"))
	{
	  i++;
	  cutgad=1;
	}
      else if (!strcmp(argv[i],"-tr"))
	{
	  i++;
          notrace=1;
	}
      else if (!strcmp(argv[i],"-ds"))
	{
	  i++;
          dostars=0;
	}
      else if (!strcmp(argv[i],"-msd"))
	{
	  i++;
	  maxstardist=atof(argv[i]);
	  i++;
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
	  //	  cmset=1;
	}
       else {
	  usage();
	}

    }

  cid=(int *) calloc(numtrace, sizeof(int));


  /*
  fp=fopen("times.dat", "r");
  if (fp!=NULL)
    {
      for (i=0; i<NUMTIMESTEPS; i++) fscanf(fp,"%f", &timestep[i]);
      fclose(fp);
    } else printf("times.dat not found\n");
  */

  dum=0;
    

  sqrenvdensrad=SQR(ENVDENSRAD); 

  if (debug) printf("comencing main loop...\n");fflush(stdout);

/************************************************************************************/
/* Start of main loop                                                               */
/************************************************************************************/



  for (snapind=endind; snapind>=startind; snapind--)
    { 
      if (debug) {printf("Ind: %d %d %d\n", snapind, endind, startind);fflush(stdout);}
      sprintf(snapfile,"%s%03d", basename, snapind);
      if (singlefile) strcpy(snapfile, basename);
      /*
      if (access(outfile, 0)!=0) {
	printf("File %s does not exits!\n", snapfile);
	continue;
      }
      */

      //      numpart=readgadget(snapfile, &evolved, &pos_ev, &vel, &id_ev, &mass_ev);
      numpart=readgadget_part(snapfile, &evolved, &part_ev);   
      if (numpart==0) 
	{
	  extern int libgaderr;
	  fprintf(stderr,"LibGad Error Code: %d\n", libgaderr);
	  continue;
	}
      convertunits(&evolved, part_ev, conv_mass, conv_dist);

      if (evolved.npart[0]==0) {dostars=0;}
      if (debug) {printf("numpart %d\n", numpart);fflush(stdout);}
#ifndef NOGAS
      if ((debug) && (evolved.npart[0]) ) {printf("sph:  %g %g\n", part_ev[0].sph->u, part_ev[0].sph->rho);fflush(stdout);}
#endif
      printf("SnapInd: %d\nStart CM search at %f %f %f\n",snapind, cm[0], cm[1], cm[2]);
      if (cid[0]) 
	{
	  k=0;
	  for (j=0; j<3; j++) cm[j]=0;
	  masstot=0;
	  qsort(cid, numtrace, sizeof(int), cmp_int);

	  for (i=evolved.npart[0]; i< numpart; i++)
	    {
	      int *fnd=bsearch(&part_ev[i].id,cid, numtrace, sizeof(int), cmp_int);
	      if (fnd!=NULL)
		{
		  for (j=0; j<3; j++) 
		    {
		      if (pb[j]) cm[j]+=MOVE(part_ev[i].pos[j], boxsz)*part_ev[i].mass;
		      else cm[j]+=part_ev[i].pos[j]*part_ev[i].mass;
		    }
		  masstot+=part_ev[i].mass;
		  k++;
		}
	      if (k==numtrace) break;
	    }

	  //	  if (debug) printf("CID %d\n", cid);
	  //set start CM
	  for (j=0; j<3; j++) 
	    {
	      cm[j]=cm[j]/masstot;
	      if (pb[j]) cm[j]==MOVE(cm[j],boxsz);
	    }
	}
      cdens=(5.36e-8) * evolved.hubparam * evolved.hubparam;
      cdens=cdens*(evolved.omegal+evolved.omega0*pow(evolved.time,-3));
      bgdens=cdens*evolved.omega0;
      if (debug) {printf("...\n");fflush(stdout);}
      
      tr_halo_cnt=0;
      /*
      for (k=0; k<3; k++) {
	cm[k]=hcm[trhalo][k];
	x=floor(cm[0]/gridsize);
	y=floor(cm[1]/gridsize);
	z=floor(cm[2]/gridsize);
      }
      */
  sprintf(idlfile ,"%s_%03d.idl",prefix, snapind);
  sprintf(asciifile ,"%s_%03d.ascii",prefix, snapind);
  sprintf( jrfile ,"%s_%03d.jr", prefix, snapind);
  sprintf(outfile ,"%s_%03d.tr", prefix, snapind);
  sprintf(gadfile ,"%s_%03d.gad",prefix, snapind);
   if ((access(outfile, 0)==0) && (donotoverwrite)) 
    {
      continue;
    }
 
   if (singlefile)
    {
        sprintf(idlfile ,"%s.idl",prefix);
        sprintf(asciifile ,"%s.ascii",prefix);
	sprintf( jrfile ,"%s.jr", prefix);
	sprintf(outfile ,"%s.tr", prefix);
	sprintf(gadfile ,"%s.gad",prefix);
	outf=stdout;
    }
   else
     outf=fopen(outfile, "w");
  fprintf(outf, "ucm %d ui %d uva %d\n", use_cm, use_inertia, use_va);
  pb[0]=0;
  pb[1]=0;
  pb[2]=0;
  masstot=0;
  double icm[3];
  for (k=0; k<3; k++)                                                         //check whether periodic boundaries are needed
    {
      if ((cm[k]<GAP) || (cm[k]>(evolved.boxsize-GAP))) 
	{
	  pb[k]=1;
	  cm[k]=MOVE(cm[k],boxsz);
	}
      icm[k]=cm[k];
    }



  /*
  index=(int *) malloc((int)(numpart*(WORKSPACE*LARGE_WS))*sizeof(int));
  if (index==NULL)
    {
      fprintf(stderr, "memory allocation failed (index)\n");
      exit(1);
    }
  icnt=0;
  for (i=(x-NBOX); i <= (x+NBOX); i++)
  for (j=(y-NBOX); j <= (y+NBOX); j++)
  for (k=(z-NBOX); k <= (z+NBOX); k++)
    {
      if (icnt>=(WORKSPACE*LARGE_WS*numpart-1))
	{
	  fprintf(stderr, "out of memory (index), increase WORKSPACE and recompile\n");
	  exit(1);
	}
      dist=sqrt(pow(x-i,2)+pow(y-j,2)+pow(z-k,2));
      if (dist<(NBOX+1))
      {
      for (l=1; l <= grid[(i+n)%n][(j+n)%n][(k+n)%n][0]; l++)
      {
	index[icnt++]=grid[(i+n)%n][(j+n)%n][(k+n)%n][l];
      }
      }
    }
  //printf("halo %d %d x %d y %d z %d\n", tr_halo , icnt, x[tr_halo], y[tr_halo], z[tr_halo]);
  index=realloc(index, sizeof(int)*icnt);
  */

if (debug) {printf("search Center\n");fflush(stdout);}
  int iteration=0;
  double searchrad=srad/1.5;
  float allocsize=WORKSPACE/8;
  if (srad>0) do
    {
      tr_halo_cnt=0;
      for (k=0; k<3; k++) cm[k]=icm[k];
      searchrad*=1.5;
      allocsize=MIN(8*allocsize, WORKSPACE* LARGE_WS);
      if (srad>0)
	{
	  tr_halo_id=(int *)malloc((int)(numpart*allocsize)*sizeof(int));
	  tr_halo_i =(int *)malloc((int)(numpart*allocsize)*sizeof(int));
	  if ((tr_halo_id==NULL) || (tr_halo_i==NULL))
	    {
	      fprintf(stderr, "memory allocation failed\n");
	      exit(2);
	    }
	  maxdist=searchrad*searchrad;
	  if (debug) {printf("!\n");fflush(stdout);}
	  for (i=0; i < numpart; i++)
	    {
	      //	      i=index[l];
	      dist=0;
	      for (j=0; j<3; j++) 
		if (pb[j]) dist+=pow(MOVE(part_ev[i].pos[j],boxsz)-cm[j],2);
		else dist+=pow(part_ev[i].pos[j]-cm[j],2);

	      //	      if (debug) {printf("%%");fflush(stdout);}
	      int type= part_ev[i].type;
      
	      if ((dist<maxdist) && ( (1<<type) & use_cm) )
		{
		  tr_halo_id[tr_halo_cnt]=part_ev[i].id;
		  tr_halo_i[tr_halo_cnt++]=i;
		}
	      if (tr_halo_cnt>=(allocsize*numpart-1))
		{
		  fprintf(stderr, "out of memory (<srad), increase WORKSPACE and recompile!\n");
		  exit(1);
		}
	    }

	  //	  tr_halo_id=realloc(tr_halo_id,sizeof(int)*tr_halo_cnt);
	  //	  tr_halo_i =realloc(tr_halo_i ,sizeof(int)*tr_halo_cnt);
  
	}

      if (++iteration>3)
	{
	  tr_halo_cnt=0;
	}
      if (debug) {printf("iteration %d\n", iteration);fflush(stdout);}
      if (tr_halo_cnt < minpartcnt) 
	{

	  break;
	}
      if (srad>0)
	{
	  maxdist=searchrad;
	  count=tr_halo_cnt;
	  do
	    {
	      masstot=0;
	      maxdist*=0.92;
	      for (j=0; j<3; j++){rdum[j]=cm[j];}
	      cm[0]=0;
	      cm[1]=0;
	      cm[2]=0;
	      for (i=0; i < count; i++)
		{
		  for (j=0; j<3; j++) 
		    {
		      if (pb[j]) cm[j]+=MOVEB(part_ev[tr_halo_i[i]].pos[j]) * part_ev[tr_halo_i[i]].mass;
		      else cm[j]+=part_ev[tr_halo_i[i]].pos[j]*part_ev[tr_halo_i[i]].mass;
		    }
		  masstot+=part_ev[tr_halo_i[i]].mass;
		}


	      for (j=0; j<3; j++) cm[j]=cm[j]/(masstot);
  
	      dum=0;
	      for (i=0; i < count; i++)                                                //
		{
		  dist=0;
		  for (j=0; j<3; j++) 
		    {		      
		      if (pb[j]) dist += pow((MOVE(part_ev[tr_halo_i[i]].pos[j],boxsz)-cm[j]),2);
		      else dist += pow((part_ev[tr_halo_i[i]].pos[j]-cm[j]),2);
		    }
      
		  dist=sqrt(dist);
		  if (dist < maxdist) 
		    {
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
			  cvel[j] += part_ev[tr_halo_i[i]].vel[j] * part_ev[tr_halo_i[i]].mass;
			}
		      masstot+= part_ev[tr_halo_i[i]].mass;
		    }
		  for (j=0; j<3; j++) cvel[j]=cvel[j]/(masstot);
		}
	    } while (count>5);
	  free(tr_halo_id);
	  free(tr_halo_i);
	}
    } while  ((ABS(cm[0]-icm[0])>(searchrad/2)) ||
	      (ABS(cm[1]-icm[1])>(searchrad/2)) ||
	      (ABS(cm[2]-icm[2])>(searchrad/2)) );
  if (debug) {printf("Center Found\n");fflush(stdout);}

  if ((tr_halo_cnt < minpartcnt) && (srad>0))
    {

      fprintf(outf, "%g %g %g <No Halo found>\n" ,cm[0], cm[1], cm[2]);
      //      free(index);
      free(tr_halo_id);
      free(tr_halo_i);
      fclose(outf);
      continue;
    }


     maxdist=GAP*GAP;
     numalloc=(int)(numpart);
     part= (struct particle *)malloc(sizeof(struct particle)*numalloc);
     vpart=(struct particle *)malloc(sizeof(struct particle)*numalloc);
     if ((part==NULL) || (vpart==NULL))
       {
	 fprintf(stderr,"memory allocation failed!\n");
	 exit(1);
       }
     if (debug) {printf("memmory allocated\n");fflush(stdout);}

     count=0;
     envdens=0;

     for (i=0; i< numpart; i++) 
       {
	 //	 i=index[l];
	 dist=0;
	 for (j=0; j<3; j++) 
	   {
	     if (pb[j]) dist += pow((MOVE(part_ev[i].pos[j],boxsz)-cm[j]),2);
	     else dist += pow((part_ev[i].pos[j]-cm[j]),2);
	     if (dist>maxdist) break;
	   }

	 if (dist < sqrenvdensrad)
	   {
	     envdens+=part_ev[i].mass;
	   } 

	 if (dist < maxdist) 
	   {
	     dist=sqrt(dist);
	     part[count].dist=dist;
	     part[count++].ind=i;
	   }
	 
	 if (count>=numalloc) 
	   {
	     fprintf(stderr, "increase Workspace!\n");
	     exit(1);
	   }

       }

     if (debug) {printf("Close particles chosen\n");fflush(stdout);}
 

     //     qsort(&part[0],count,sizeof(struct particle), (__compar_fn_t)cmp);                  //sort particles by distance to CM of halo
     qsort(&part[0],count,sizeof(struct particle), cmp);                  //sort particles by distance to CM of halo
     //     part=realloc(part,sizeof(struct particle)*count);
     masstot=0;
     i=0;
     od=201;
     effrad=30*evolved.hubparam;
     meff=0;
     peakrotvel=0;
     double b_mass=0;
     double dm_mass=0;
     while ((od > 200) || (i<5))                                                    //Add mass until overdensity drops below 200
       {
	 int ind= part[i].ind;
	 masstot+=part_ev[ind].mass;
	 if ((part_ev[ind].type == 0) ||(part_ev[ind].type == 4)) 
	   {
	     b_mass += part_ev[ind].mass;
	   } else 
	   {
	     dm_mass+= part_ev[ind].mass;
	   }
	 if ((part[i].dist < effrad) && ( (1<<part_ev[ind].type) & use_cm ) )
	   {
	     meff +=part_ev[ind].mass;
	   }

	 od=masstot/(pow(part[i].dist*evolved.time,3)*(4.0/3.0)*PI*cdens);
	 //printf("i %d od %f\n",i,od);
	 rotvel=sqrt((masstot/(part[i].dist*evolved.time))*(G*Msun*1e10/kpc))*1e-3;
	 if (rotvel > peakrotvel) 
	   {
	     peakrotvel= rotvel;
	     peakr=part[i].dist;
	   }
	 vpart[i].dist=part[i].dist;
         vpart[i].ind=part[i].ind;
	 i++;
	 if (i > count)
	   {
	     fp=fopen("error_trhalo.dat","a");
	     fprintf(fp,"Halo %d is making trouble\n",snapind);
	     fclose(fp);
	     printf("halo %d is making trouble\n",snapind);
	     break;
	   }
       }

     if (debug) {fprintf(outf,"Virial radius calculated\n");fflush(stdout);}
     if (use_cm==1) fprintf(outf,"Gas mass inside 30 kpc %f\n", meff);
     if (use_cm==16) fprintf(outf,"Stellar mass inside 30 kpc %f\n", meff);
     if (use_cm==17) fprintf(outf,"Baryonic mass inside 30 kpc %f\n", meff);
 
     masstot-=part_ev[part[i-1].ind].mass;
     if ((part_ev[part[i-1].ind].type == 0) ||(part_ev[part[i-1].ind].type == 4)) 
	   {
	     b_mass -= part_ev[part[i-1].ind].mass;
	   } else 
	   {
	     dm_mass-= part_ev[part[i-1].ind].mass;
	   }
     dum=i-1;  
     vcnt=i;
     halorad=part[vcnt-1].dist;
     vpart=realloc(vpart,sizeof(struct particle)*vcnt);
     //     qsort(&vpart[0],vcnt,sizeof(struct particle),(__compar_fn_t)cmpind);
     qsort(&vpart[0],vcnt,sizeof(struct particle),cmpind);

     j=0;
     k=0;
     while (k<numtrace)
       {
	 if (part[j].ind>evolved.npart[0])
	   {
	     cid[k++]=part_ev[part[j].ind].id;
	   }
	 if ((++j) > vcnt) 
	   {
	     fprintf(stderr, "not enough DM particles inside the virial radius (%d).\n", k);
	     numtrace=k;
	   }
       }
     if (debug) {printf("Innermost DM particles determined\n");fflush(stdout);}

     envdens=envdens/((4.0/3.0)*PI*pow(ENVDENSRAD*evolved.time,3));


 for (j=0; j<3; j++) 
   {
     if (pb[j])
       {
	 cm[j]=MOVE(cm[j],boxsz);
	 fprintf(outf,"!!!periodic boundaries used in Dimension %d !!!\n",j);
       }
   }
 vmass=masstot;
 fprintf(outf,"redshift: %f\naexp: %f\n", 1./evolved.time - 1,evolved.time);
 fprintf(outf,"\n halo in snapshot %5d consists of %5d particles, virial Mass/h %5.2f\n rad %5.2f         od %4.2f\n",snapind,vcnt,masstot,halorad,od);
 fprintf(outf," Baryonic Mass/h %5.2f\n DM Mass/h %5.2f \n\n",b_mass, dm_mass);
 
 fprintf(outf,"Center of Mass  : %12.2f %12.2f %12.2f\n",cm[0],cm[1],cm[2]);
 fprintf(outf,"shift to cm-file: %+12.2f %+12.2f %+12.2f\n",cm[0]-icm[0],cm[1]-icm[1],cm[2]-icm[2]);
 fprintf(outf,"CM Vel          : %12.2f %12.2f %12.2f\n",cvel[0],cvel[1],cvel[2]);
 fprintf(outf,"particles in box: %d\n",count);
 fprintf(outf,"Environmental Radius: %f \n", ENVDENSRAD);
 fprintf(outf,"Environmental OD: %f %f\n", envdens/cdens, envdens/bgdens);
 fprintf(outf,"Peakrotvel: %f at %f kpc\n",peakrotvel, peakr);
 fprintf(outf,"    Rotvel: %f at Rvir\n",rotvel);


 tr_halo_id=(int *)malloc(count*sizeof(int));
 tr_halo_i =(int *)malloc(count*sizeof(int));
 if ((tr_halo_id==NULL) || (tr_halo_i==NULL))
	    {
	      fprintf(stderr, "memory allocation failed (tr_halo_?)\n");
	      exit(2);
	    }
 tr_halo_cnt=0;


 if (cutgad)
   {
     gadpart *outpart= (gadpart *) malloc (sizeof(gadpart)* vcnt);
     if ((outpart==NULL) )
	    {
	      fprintf(stderr, "memory allocation failed (outpart)\n");
	      exit(2);
	    }
//     outpos =(fltarr *)malloc(sizeof(fltarr)*i);
//     outvel =(fltarr *)malloc(sizeof(fltarr)*i);
//     outid =(int *)malloc(sizeof(int)*i);
//     outmass =(float *)malloc(sizeof(float)*i);
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
	 outpart[i]=part_ev[k];
	 for (j=0; j<3; j++)
	   {

	     if (pb[j]) outpart[i].pos[j]=MOVEB(outpart[i].pos[j]) - MOVEB(cm[j]);
	     else outpart[i].pos[j]=(outpart[i].pos[j]) - (cm[j]);
	     outpart[i].vel[j] -= cvel[j];
	   }
//	 dum=0;
//	 for (l=0; l<6; l++)
//	   {
// 	     dum+=evolved.npart[l];
//	     if (k<dum) break;
//	   }
	 out.npart[outpart[i].type]++;
	 out.nall[outpart[i].type]++;
	 i++;
       }
     if (debug) {printf("...writing output file %d %d\n", vcnt, i);fflush(stdout);}
     //     if (debug) {printf("...data-test %g %g %g\n", outpart[0].sph->rho, outpart[0].sph->u, outpart[0].pot);fflush(stdout);}
     writegadget_part(gadfile, out, outpart);
     if (debug) {printf("Gadget File written\n");fflush(stdout);}
     free(outpart);
//     free(outpos);
//     free(outvel);
//     free(outid);
//     free(outmass);

   }


 maxdist=halorad*tr_factor;
 double hmeff=0;
 double galmass=0;
 double gasmass=0;
 double DMmass=0;
 double innertotm=0;
 double galrad=0;
 double gsfr=0;
 double meanage=0;
 int nstars_gal=0;
 dist=0;
 i=0;
// for (j=0; j<50; j++)
//   {
//     p[j]=0;
//     err[j]=0;
//     rlog[j]=0;
//   }
 while (dist<maxdist) 
   { 
     dist=part[i].dist;
     m=part_ev[part[i].ind].type;
     if ((m==4) && (dostars) && (dist < halorad) && (snapind==endind)) totalstarcnt++;
     else if ((m==0) && (dostars) && (dist < (halorad * ACC_FRAC)) && (snapind==endind)) totalgascnt++;

//     if ((dist<halorad) && (m>0) && (m<4) && (dist>sfl))
//       {
//	 d=log10(dist);
//	 j=floor(d/(log10(halorad)/50));
//	 err[j]++;
//	 p[j]+=part_ev[part[i].ind].mass;
//       }
     if  ( ((1<<m) & use_cm) && (hmeff < (meff/2.0)) )
       {
	 hmeff+=part_ev[part[i].ind].mass;
	 effrad=dist;
       }
     if ((dist < (halorad * 0.1) ) )
       {
	 innertotm+=part_ev[part[i].ind].mass;
	 if ((m==4)) 
	   {
	     galmass+= part_ev[part[i].ind].mass;
	     double redsh = (1/part_ev[part[i].ind].stellarage) - 1;
	     double sage = TIMEDIFF(redsh, (1/evolved.time)-1 );
	     meanage += sage;
	     nstars_gal++;
	   }
	 else if ((m==0)) 
	   {
	     gasmass+= part_ev[part[i].ind].mass;
	     gsfr+= part_ev[part[i].ind].sph->sfr;
	   }
	 else if ((m==1)) DMmass+= part_ev[part[i].ind].mass;
       }
     i++;
   }
 meanage /= nstars_gal;
 if (debug) {printf("totalstarcnt %d\ntotalgascnt %d\n", totalstarcnt,totalgascnt);fflush(stdout);}
 if (use_cm&16)
   {
     i=0;
     double mdum=0;
     double total_mass=0;
     double dm_mass=0;
     dist=0;
     while (dist<maxdist) 
       { 
	 total_mass += part_ev[part[i].ind].mass;
	 dist=part[i].dist;
	 m=part_ev[part[i].ind].type;	 
	 if  (m==4) 
	   {
	     if (mdum < (galmass/2.0)) 
	       {
		 mdum+=part_ev[part[i].ind].mass;
		 galrad=dist;
	       } else break;
	   } else if ((m>0) && (m<4))
	   {
	     dm_mass += part_ev[part[i].ind].mass;
	   }

	 i++;
       }
      fprintf(outf,"3D half-mass-galaxy radius: %g\n", galrad);
      fprintf(outf,"total mass inside 3D half-mass-galaxy radius: %6g\n", total_mass);
      fprintf(outf,"stellar mass inside 3D half-mass-galaxy radius: %6g\n", mdum);
      fprintf(outf,"dark matter mass inside 3D half-mass-galaxy radius: %6g\n", dm_mass);
   }
 fprintf(outf,"hmeff: %f \neffrad: %f\n", hmeff, effrad);
 fprintf(outf,"Galaxy mass (M_star < 10 %% Rvir): %f \n", galmass);
 fprintf(outf,"mean age: %g \n", meanage);
 fprintf(outf,"SFR (< 10 %% Rvir): %g \n", gsfr);
 fprintf(outf,"specific SFR (< 10 %% Rvir): %g \n", gsfr / galmass);
 fprintf(outf,"Gas mass (M_gas < 10 %% Rvir): %f \n", gasmass);
 fprintf(outf,"darkmatter mass (M_dm < 10 %% Rvir): %f \n", DMmass);
 fprintf(outf,"Total mass < 10 %% Rvir: %f \n", innertotm);
 if (debug) {printf("Stars: %d\ni %d\n",totalstarcnt, i);fflush(stdout);}
 /****************************************************************************************/
 /*Find Star Particles inside Rvir********************************************************/
 if (dostars)
   {
     if (debug) {printf("%d %d %d\n", snapind, endind, (snapind==endind));fflush(stdout);}
     if (snapind==endind)
       {
	 starcnt=totalstarcnt;
	 gascnt =totalgascnt;
	 stars = (struct star*) calloc(starcnt, sizeof(struct star));
	 gas   = (struct gaspart*) calloc(gascnt, sizeof(struct gaspart));
	 if (debug) {printf("Memory for Stars allocated\n");fflush(stdout);}
	 j=0;
	 l=0;
	 k=0;
	 while ( (j<starcnt) || (k<gascnt) )
	   {
	     //	     idum=0;
//	     for (m=0; m<6; m++)
//	       {
//		 idum+=evolved.npart[m];
//		 if (idum>part[l].ind) break;
//	       }
	     if (j>starcnt) exit(1);
	     m = part_ev[part[l].ind].type;
	     if ((m==4) && (part[l].dist < halorad))
	       {
		 stars[j].id   =part_ev[part[l].ind].id;
		 stars[j].idist=part[l].dist;
		 stars[j].ifrac=part[l].dist/halorad;
		 stars[j].ifrac_rhalf=part[l].dist/galrad;
		 stars[j].isnap=snapind;
		 stars[j].dist=part[l].dist;
		 stars[j].frac=part[l].dist/halorad;
		 if (stars[j].frac < ACC_FRAC) stars[j].a_acc=evolved.time;
		 else stars[j].a_acc=0;
		 stars[j].snap=snapind;
		 stars[j].a=evolved.time;
#ifdef POTENTIAL
		 stars[j].pot=part_ev[part[l].ind].pot;
#endif
		 j++;
	       }
	     else if ((m==0) && (part[l].dist < halorad * ACC_FRAC))
	       {

		 gas[k].id = part_ev[part[l].ind].id;
		 double temp = temperature(part_ev[part[l].ind]);
		 float vcirc = sqrt((vmass/(halorad * evolved.time))*(G*Msun*1e10/kpc))*1e-3;
		 gas[k].Tvir_acc = SQR(vcirc / 167.) * 1e6;
		 double rho = part_ev[part[l].ind].sph->rho * SQR(HUB) * pow(evolved.time, -3);
		 if (rho < DENS_THRESH)
		   {
		     gas[k].maxtemp = temp;
		     gas[k].a_maxtemp = evolved.time;
		     gas[k].frac_maxtemp = part[l].dist / halorad;
		     gas[k].Tvir = SQR(vcirc / 167.) * 1e6;
		   }
		 else 
		   {
		     gas[k].maxtemp = 0;
		     gas[k].a_maxtemp = 0;
		     gas[k].frac_maxtemp = 0;
		     gas[k].Tvir = 0;
		   }
		 gas[k].a_acc = evolved.time;
		 gas[k].T_acc = temp;
		 gas[k].a_star = 0;
		 k++;
	       }
	     l++;
	   }
	 if (debug) {printf("sorting initial stars and gas...\n");fflush(stdout);}
	 //	 qsort(stars, starcnt, sizeof(struct star), (__compar_fn_t)cmp_star_id);
	 qsort(stars, starcnt, sizeof(struct star)   , cmp_star_id);
	 qsort(gas  , gascnt , sizeof(struct gaspart), cmp_gas_id);
	 if (debug) {printf("Stars sorted\n");fflush(stdout);}
       }
     else
       {
	 j=0;
	 l=0;
	 k=0;
	 n=0;
	 h=0;
	 struct star stardum;
	 struct gaspart gasdum;
	 if (debug) {printf("Search Stars..%d\n", starcnt);fflush(stdout);}
	 while ( ( (j<starcnt) || (k<gascnt) )  && (l < count))
	   {
//	     idum=0;
//	     for (m=0; m<6; m++)
//	       {
//		 idum+=evolved.npart[m];
//		 if (idum>part[l].ind) break;
//	       }
	     m= part_ev[part[l].ind].type;
	     //	     if (debug) {printf("#");fflush(stdout);}
	     if ((m==0) || (m==4))
	       {
		 stardum.id=part_ev[part[l].ind].id;
		 struct star *fnd=bsearch(&stardum, stars, starcnt, sizeof(struct star), cmp_star_id);
		 if (fnd!=NULL)
		   {
		     j++;
		     //		     if (debug) {printf("Found one!...\n");fflush(stdout);}
		     if (m==4)
		       {
			 fnd->dist=part[l].dist;
			 fnd->frac=part[l].dist/halorad;
			 if (fnd->frac < ACC_FRAC) fnd->a_acc=evolved.time;
			 fnd->a=evolved.time;
			 fnd->snap=snapind;
#ifdef POTENTIAL
			 fnd->pot=part_ev[part[l].ind].pot;
#endif			 
		       } else 
		       {
			 fnd->fnd=1; 
			 fnd->gasdist=part[l].dist;
			 if ( (fnd->ifrac <= ACC_FRAC)  )
			   {
			     addgas[h++]=l;
			     if (h>=ADDGAS_BUF) {fprintf(stderr, "Increase ADDGAS_BUF\n");exit(1);}
			   }
		       }
		   } else
		   {
		     if (m==4)
		       {
			 addstar[n++]=l;
			 if (n>=ADDSTAR_BUF) {fprintf(stderr, "Increase ADDSTAR_BUF\n");exit(1);}
		       }
		     else
		       {
			 gasdum.id=part_ev[part[l].ind].id;
			 struct gaspart *fndgas=bsearch(&gasdum, gas, gascnt, sizeof(struct gaspart), cmp_gas_id);
			 if (fndgas != NULL)
			   {
			     double temp = temperature(part_ev[part[l].ind]);
			     double rho = part_ev[part[l].ind].sph->rho * SQR(HUB) * pow(evolved.time, -3);
			     if ((fndgas->maxtemp < temp) && (rho < DENS_THRESH))
			       {
				 fndgas->maxtemp = temp;
				 float vcirc = sqrt((vmass/(halorad * evolved.time))*(G*Msun*1e10/kpc))*1e-3;
				 fndgas->Tvir = SQR(vcirc / 167.) * 1e6;
				 fndgas->a_maxtemp = evolved.time;
				 fndgas->frac_maxtemp = part[l].dist/halorad;
			       }
			     if (part[l].dist < (halorad * ACC_FRAC))
			       {
				 fndgas->a_acc = evolved.time;
				 fndgas->T_acc = temp;
				 float vcirc = sqrt((vmass/(halorad * evolved.time))*(G*Msun*1e10/kpc))*1e-3;
				 fndgas->Tvir_acc = SQR(vcirc / 167.) * 1e6;
			       }
			     k++;
			   }

		       }
		   }
		 //		 j++;
	       }
	     l++;
	   }
	 if (debug) {printf("Write Stars to list %d..\n",n);fflush(stdout);}
	 if (n)
	   {
	     totalstarcnt+=n;
	     stars=realloc(stars, totalstarcnt* sizeof(struct star));
	     k=totalstarcnt-1;
	     while (k>=n)
	       {
		 stars[k]=stars[k-n];
		 k--;
	       }
	     for (k=0; k<n; k++)
	       {
		 stars[k].id   =part_ev[part[addstar[k]].ind].id;
		 stars[k].idist=part[addstar[k]].dist;
		 stars[k].ifrac=part[addstar[k]].dist / halorad;
		 stars[k].ifrac_rhalf=part[addstar[k]].dist / galrad;
		 stars[k].isnap=snapind;
		 stars[k].dist=part[addstar[k]].dist;
		 stars[k].frac=part[addstar[k]].dist/halorad;
		 if (stars[k].frac < ACC_FRAC) stars[k].a_acc=evolved.time;
		 else stars[k].a_acc=0;
		 stars[k].snap=snapind;
		 stars[k].a=evolved.time;
#ifdef POTENTIAL
		 stars[k].pot=part_ev[part[addstar[k]].ind].pot;
#endif
	       }
	     starcnt+=n;
	   }
	 if (h)
	   {
	     totalgascnt+=h;
	     gas= realloc(gas, totalgascnt * sizeof(struct gaspart));
	     k=totalgascnt-1;
	     while (k>=h)
	       {
		 gas[k] = gas[k-h];
		 k--;
	       }
	     for ( k = 0; k < h; k++ )
	       {
		 idum = addgas[k];
		 gas[k].id = part_ev[part[idum].ind].id;
		 double temp = temperature(part_ev[part[idum].ind]);
		 float vcirc = sqrt((vmass/(halorad * evolved.time ))*(G*Msun*1e10/kpc))*1e-3;
		 double rho = part_ev[part[idum].ind].sph->rho * SQR(HUB) * pow(evolved.time, -3);
		 if (rho < DENS_THRESH)
		   {
		     gas[k].maxtemp = temp;
		     gas[k].Tvir = SQR(vcirc / 167.) * 1e6;
		     gas[k].a_maxtemp = evolved.time;
		     gas[k].frac_maxtemp = part[idum].dist / halorad;
		   }
		 else
		   {
		     gas[k].maxtemp = 0;
		     gas[k].Tvir = 0;
		     gas[k].a_maxtemp = 0;
		     gas[k].frac_maxtemp = 0;
		   }
		 if (part[idum].dist < (halorad * ACC_FRAC))
		   {
		     gas[k].a_acc = evolved.time;
		     gas[k].T_acc = temp;
		     gas[k].Tvir_acc = SQR(vcirc / 167.) * 1e6;
		   } 
		 else 
		   {
		     gas[k].a_acc=0;
		     gas[k].T_acc=0;
		     gas[k].Tvir_acc=0;
		   }
		 gas[k].a_star = evolved.time;
	       }
	     gascnt+=h;
	   }
	 
	 j=0;
	 k=starcnt-1;
	 while ((stars[k].fnd) && (k)) {k--;}
	 l=0;
	 if (debug) {printf("resorting Stars \n");fflush(stdout);}
	 while (j<k)
	   {
	     if (stars[j].fnd)
	       {
		 stardum=stars[j];
		 stars[j]=stars[k];
		 stars[k]=stardum;
		 while ((stars[k].fnd) && (k)) {k--;}
	       }
	     j++;
	   }
	 if (!stars[j].fnd) j++;
	 starcnt=j;
	 //	 qsort(stars, starcnt, sizeof(struct star), (__compar_fn_t)cmp_star_id);
	 qsort(stars, starcnt, sizeof(struct star), cmp_star_id);
	 qsort(gas, gascnt, sizeof(struct gaspart), cmp_gas_id);
       }
   }
 /****************************************************************************************/
 if (debug) printf("fitting...\n");fflush(stdout);
// nfun=0;
// k=0;
// d=log10(halorad)/50;
// for (j=0; j<50; j++)
//   {
//     if (err[j]!=0)
//       {
//     err[nfun]=1/sqrt(sqrt(err[j]));
//     if (k==0) {
//       rlog[nfun]=(pow(10, ((j+1)*d )))/2;
//       p[nfun]=p[j]/((4.0/3.0)* PI * (pow(10, ((j+1)*d*3))));
//       //       printf("j %d\n", j);
//
//     }
//       else {
//	 ddum=p[j];
//	 rlog[nfun]=(pow(10, ((j+1)*d )) + pow(10, (k*d)))/2;
//	 p[nfun]=p[j]/((4.0/3.0)* PI * (pow(10, ((j+1)*d*3)) - pow(10, (k*d*3))));
//	 if (p[nfun]<0) fprintf(outf,"%f %e %e \n", ddum, pow(10, ((j+1)*d*3)), pow(10, (k*d*3)) );
//       }
//     /*
//     if (rlog[nfun] < SOFTENING) 
//       {
//	 err[nfun]*=10;
//	 printf("S %d\n", nfun);
//       }
//     */
//     k=j+1;
//     nfun++;
//       } else {
//	 //	 p[j]=0;
//	 //	 err[j]=1e10;
//       }
//
//
//   }
// 
// for (j=0; j<nfun; j++)
//   {
//          p[j]=log10(p[j]);
//   }
// 
// // auxiliary settings for fitting:
//
//    lm_control_type control;
//    lm_data_type data;
//    lm_initialize_control(&control);
//
//    data.user_func = fitfct;
//    data.user_t = rlog;
//    data.user_y = p;
//---------------------------------------------------
//
// par[0]=1;
// par[1]=10;
// lm_minimize(nfun, 2, par, err, lm_evaluate_default, dumprintout, &data, &control);

 // printf("%f %f\n", d,  pow(10,d*50));
 i--;
 FILE *fp2;
 if (doascii)
   {
     fp2=fopen(asciifile,"w");
     
   }
 fp=fopen(idlfile,"w");
 fprintf(fp," %d %f %f %f %f %f\n",i,masstot,halorad,cm[0],cm[1],cm[2]);
 dist=0;
 i=0;
 while (dist < maxdist)
   {
     tr_halo_id[tr_halo_cnt]=part_ev[part[i].ind].id;
     tr_halo_i[tr_halo_cnt++]=part[i].ind;
     dist=part[i].dist;
     double temp;
     if (part_ev[part[i].ind].type == 0)
       {
	 temp = temperature(part_ev[part[i].ind]);
       }
     else
       {
	 temp=0;
       }
     fltarr veldum, raddum;
     for ( j = 0; j < 3; j++ )
       {
	 veldum[j] = part_ev[part[i].ind].vel[j] - cvel[j];
	 raddum[j] = part_ev[part[i].ind].pos[j] - cm[j];
       }
     double vrad = radvel(veldum, raddum);
     fprintf(fp,"%f %f %6d %3d %9.7g %8.2f\n",dist, part_ev[part[i].ind].mass, part_ev[part[i].ind].id, part_ev[part[i].ind].type, temp, vrad);
     if (doascii)
       {
	 for (j=0; j<3; j++) fprintf(fp2, "%12g ", part_ev[part[i].ind].pos[j]-cm[j]);
	 for (j=0; j<3; j++) fprintf(fp2, "%12g ", part_ev[part[i].ind].vel[j]-cvel[j]);
	 fprintf(fp2, " %6d ", part_ev[part[i].ind].id);
	 fprintf(fp2, " %6g ", part_ev[part[i].ind].mass);
	 fprintf(fp2, "\n");
       }
     i++;
   }
 fprintf(outf,"number of particles within %f times virial radius %d\n", tr_factor ,i);
 if (doascii)
   {
     fclose(fp2);
   }
 fclose(fp);
 if (debug) printf("realloc...\n");fflush(stdout);
 tr_halo_id=realloc(tr_halo_id,sizeof(int)*tr_halo_cnt);
 tr_halo_i =realloc(tr_halo_i ,sizeof(int)*tr_halo_cnt);


 /**********************************************************************************************************************************************/
 if (debug) printf("lambda...\n");fflush(stdout);
 // for (l=-5; l<=5; l++)
   {
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
	 vcm[k]+=(part_ev[part[j].ind].vel[k] * sqrt(evolved.time)*part_ev[part[j].ind].mass);
	 if (pb[k]) jcm[k]+=(MOVE(part_ev[part[j].ind].pos[k],boxsz)*part_ev[part[j].ind].mass);
	 else  jcm[k]+=(part_ev[part[j].ind].pos[k]*part_ev[part[j].ind].mass);
       }
      massdum+=part_ev[part[j].ind].mass;
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
	 if (pb[k]) rdum[k]=MOVE(part_ev[part[j].ind].pos[k],boxsz)-jcm[k];
	 else rdum[k]=part_ev[part[j].ind].pos[k]-jcm[k];
	 vdum[k]=part_ev[part[j].ind].vel[k]-vcm[k];
	 vdum[k] *= sqrt(evolved.time);
       }
     torq[0]=(rdum[1]*vdum[2]-rdum[2]*vdum[1]);
     torq[1]=(rdum[2]*vdum[0]-rdum[0]*vdum[2]);
     torq[2]=(rdum[0]*vdum[1]-rdum[1]*vdum[0]);

     ddum=0;
     for (k=0; k<3; k++) 
       {
	 J[k]+=torq[k]*part_ev[part[j].ind].mass;
	 ddum+=torq[k]*torq[k];
       }
     totj+=sqrt(ddum)*part_ev[part[j].ind].mass;
     
     ddum=0;
     ddum1=0;
     for (k=0; k<3; k++) {ddum+=SQR(rdum[k]); ddum1+=SQR(torq[k]);}
     if (dojr) fprintf(fp,"%g %g\n", sqrt(ddum), sqrt(ddum1)*part_ev[part[j].ind].mass);
    }
 if (dojr) fclose(fp);


 lmd=0;
  for (k=0; k<3; k++) {lmd+=J[k]*J[k];}
  lmd=(sqrt(lmd)*1e10*(Msun/evolved.hubparam)*1e3*(kpc/evolved.hubparam))/(sqrt(2)*part[count-1].dist*(kpc/evolved.hubparam)*massdum*1e10*(Msun/evolved.hubparam));
 // lmd=(totj*1e10*(Msun/evolved.hubparam)*1e3*(kpc/evolved.hubparam))/(sqrt(2)*part[count-1].dist*(kpc/evolved.hubparam)*massdum*1e10*(Msun/evolved.hubparam));
 lmd=lmd/sqrt(G*massdum*1e10*(Msun/evolved.hubparam)/(part[count-1].dist*(kpc/evolved.hubparam)));

 fprintf(outf,"Spin: lamda %e  r %g\n",lmd, part[count-1].dist);


 if (debug) printf("ui...\n");fflush(stdout);
 /*************************************************************************************************************************************************************/
 double q=1, s=1;
 double Pi[3]={0.0, 0.0, 0.0};
 double Pi_PHI=0;
 double Pi_RR=0;
 double delta;
 if (use_inertia)
   {
     maxdist=rad_ratio*halorad;
     if ((use_inertia&16)) 
       {
	 maxdist=galrad;
	 //	 maxdist=10 * HUB; //think about it (used in feldmann 2010)
       }
     double s_old;
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
	 wpart[k]=part_ev[l];
	 for (j=0; j < 3; j++)
	   {
	     if (pb[j]) wpart[k].pos[j]=MOVEB(part_ev[l].pos[j])-MOVEB(cm[j]);
	       else wpart[k].pos[j]=part_ev[l].pos[j]-cm[j];
	     wpart[k].vel[j] -= cvel[j];
	   }

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
  gsl_vector_free(oldpos);
  gsl_vector_free(newpos);
  gsl_vector_free(oldvel);
  gsl_vector_free(newvel);
  gsl_permutation_free(perm);
     } while ((ABS(s_old-s)/s) > 1e-2);
     
     char rotgadfile[128];
     if (singlefile)
       {
	 sprintf(rotgadfile ,"%s.rot",prefix);
       } else
       sprintf( rotgadfile, "%s_%03d.rot", prefix, snapind);

     FILE *matrixf=fopen(rotgadfile,"w");
     gsl_matrix_fwrite (matrixf, rotation);
     fclose(matrixf);

     fprintf(outf,"longest axis = %g * R_vir\nq = %g\ns = %g\n", rad_ratio, q, s);


     /***********************************************************************************************************************************/
     /*Density Profile Fit      *********************************************************************************************************/
     gadpart_dist * gpart =malloc(vcnt* sizeof(gadpart_dist));
     for (i=0; i< vcnt; i++) 
       {
	 gpart[i].part=wpart[i];
	 for (j=0; j<3; j++) 
	   {
	     gpart[i].part.vel[j] *= sqrt(evolved.time);
	   }
	 gpart[i].dist = sqrt(SQR(gpart[i].part.pos[0]) + SQR(gpart[i].part.pos[1]) + SQR(gpart[i].part.pos[2])) * evolved.time;
       }
     qsort(gpart, vcnt, sizeof(gadpart_dist), cmp_dist);

     double rcs;
     double concentration = nfwfit(par, gpart, vcnt, halorad * evolved.time, sfl*evolved.time, &rcs);

     fprintf(outf,"Delta_c %f   Scale Radius %f   Concentration factor c %f\n", par[0], par[1], halorad/par[1]);
     fprintf(outf,"Error NFW: %6g\n", rcs);

     double gamma = densproffit(par, gpart, vcnt, galrad*evolved.time, sfl*evolved.time, &rcs, 63);
     fprintf(outf,"logarithmic slope of density profile: %6g\n", gamma);
     fprintf(outf,"rho0 of fit: %6g\n", par[0]);
     fprintf(outf,"Error profile-fit: %6g\n", rcs);
     /***********************************************************************************************************************************/

     
     //***********************************************
     //Calculate projected half-mass-radii and VA
     if (use_cm == 16)
       {
	 
	 for (j=0; j<3; j++)
	   {
	     int dim[3];
	     dim[0]=j;
	     dim[1]=(j+1)%3;
	     dim[2]=(j+2)%3;
	     for (i=0; i< vcnt; i++) gpart[i].dist=sqrt(SQR(gpart[i].part.pos[dim[1]]) + SQR(gpart[i].part.pos[dim[2]]));
	     //	     fprintf(outf,"%g\n", gpart[10].dist);
	     //	     qsort(gpart, vcnt, sizeof(gadpart_dist),(__compar_fn_t)  cmp_dist);
	     qsort(gpart, vcnt, sizeof(gadpart_dist), cmp_dist);

	     /*
	       r_e apperture
	      */

	     int count=0;
	     double prad=0;
	     double mdum=0;
	     double veltot=0;
	     double sqrveltot=0;
	     k=0;
	     while (mdum < (galmass/2.0))
	       {
		 if (gpart[k].part.type == 4)
		   {
		     mdum+=gpart[k].part.mass;
		     prad =gpart[k].dist;
		     veltot   +=     gpart[k].part.vel[dim[0]] * gpart[k].part.mass;
		     sqrveltot+= SQR(gpart[k].part.vel[dim[0]]) *gpart[k].part.mass;
		     count++;
		   }
		 k++;
	       }
	     double halfmassrad = gpart[k-1].dist;
	     double velmean   =    veltot / mdum;
	     double sqrvelmean= sqrveltot / mdum;
	     double sigm = sqrvelmean - SQR(velmean);
	     fprintf(outf,"DIM %d (%d)| projected hm-radius: %6g | mean velocity variance %6g\n", j, count, prad, sqrt(sigm));

	     /*
	       1/2 * r_e aperture
	      */

	     count = 0;
	     prad  = 0;
	     mdum  = 0;
	     veltot= 0;
	     double total_mass= 0;
	     sqrveltot = 0;
	     k = 0;
	     while (gpart[k].dist <= (halfmassrad / 2.))
	       {
		 if (gpart[k].part.type == 4)
		   {
		     mdum+=gpart[k].part.mass;
		     prad =gpart[k].dist;
		     veltot   +=     gpart[k].part.vel[dim[0]] * gpart[k].part.mass;
		     sqrveltot+= SQR(gpart[k].part.vel[dim[0]]) *gpart[k].part.mass;
		     count++;
		   }
		 total_mass += gpart[k].part.mass;
		 k++;
	       }
	     velmean = veltot / mdum;
	     sqrvelmean = sqrveltot / mdum;
	     sigm = sqrvelmean - SQR(velmean);
	     fprintf(outf,"0.5 * R_1/2 aperture: DIM %d (%d)| projected hm-radius: %6g | mean velocity variance %6g | proj. stellar mass %6g | total proj. mass %6g\n", j, count, prad, sqrt(sigm), mdum, total_mass);


	     /*
	       1 kpc fixed radius
	      */

	     count = 0;
	     prad  = 0;
	     mdum  = 0;
	     veltot= 0;
	     sqrveltot = 0;
	     k = 0;
	     while (gpart[k].dist <= ( 1. * evolved.time ))
	       {
		 if (gpart[k].part.type == 4)
		   {
		     mdum+=gpart[k].part.mass;
		     prad =gpart[k].dist;
		     veltot   +=     gpart[k].part.vel[dim[0]] * gpart[k].part.mass;
		     sqrveltot+= SQR(gpart[k].part.vel[dim[0]]) *gpart[k].part.mass;
		     count++;
		   }
		 k++;
	       }
	     velmean = veltot / mdum;
	     sqrvelmean = sqrveltot / mdum;
	     sigm = sqrvelmean - SQR(velmean);
	     fprintf(outf,"1 kpc fixed: DIM %d (%d)| projected hm-radius: %6g | mean velocity variance %6g | proj. stellar mass %6g\n", j, count, prad, sqrt(sigm), mdum);


	     
	   }
	

       }
     free(gpart);


 if (debug) printf("uva...\n");fflush(stdout);
/*************************************************************************************/     
/*Calculate Velocity anisotropy                                                      */

     gadpart      * vapart=malloc(vcnt* sizeof(gadpart));
     gadpart_dist * dpart =malloc(vcnt* sizeof(gadpart_dist));
     int num_va=0;
     for (i=0; i< vcnt; i++)
       {
	 if ((1<<(wpart[i].type))&use_va)
	   {
	     //	     cpygadpart(&vapart[num_va]    , &wpart[i]);
	     //	     cpygadpart(&dpart[num_va].part, &wpart[i]);
	     vapart[num_va]=wpart[i];
	     dpart[num_va].part=wpart[i];
	     dpart[num_va].dist=sqrt(SQR(wpart[i].pos[0])+SQR(wpart[i].pos[1]/q)+SQR(wpart[i].pos[2]/s));
	     num_va++;
	   }
       }
     if (debug) printf("particles copied for VA-calculation...\n");fflush(stdout);
     if (num_va > 50)
       {

	 
	 //	 qsort(&dpart[0], num_va, sizeof(gadpart_dist), (__compar_fn_t) cmp_dist);
	 qsort(&dpart[0], num_va, sizeof(gadpart_dist),  cmp_dist);

	 int distnum= gadsearch(dpart, maxdist, 0, num_va);
	 //	 dpart= realloc(dpart, distnum*sizeof(gadpart_dist));

	 if ((test) && (distnum==-1) )
	   printf("!num_va: %d %f %f %d\n", num_va, q, s, snapind);
	 


	 if (debug) printf("Build Tree...\n");fflush(stdout);
	 KdNode * root;
	 initKdNode(&root, NULL);
	 buildKdTree(root, vapart, num_va, 0);
	 if (debug) printf("Tree built...\n");fflush(stdout);
	 if (debug) {
	   printf("check Tree...%d ?= %d\n", checkKdTree(root), num_va);
	   printf("min: %f %f %f %f %f\n", root->min, root->down->min, root->up->min, root->down->down->min, root->up->up->min);
	   printf("max: %f %f %f %f %f\n", root->max, root->down->max, root->up->max, root->down->down->max, root->up->up->max);
	   printf("dim: %d %d %d %d %d\n", root->dim, root->down->dim, root->up->dim, root->down->down->dim, root->up->up->dim);
	   fflush(stdout);
	 }
	 for (i=0; i<3; i++) Pi[j]=0;
	 int incstep;
	 if (distnum<120) incstep=1;
	 else if (distnum<300) incstep=5;
	 else if (distnum<500) incstep=10;
	 else if (distnum<5000) incstep=30;
	 else if (distnum<100000) incstep=50;
	 else incstep=250;



	 for (i=0; i<distnum; i+=incstep)
	   {
	     masstot=0;
	     gadpart_dist * knn;
	     double knndist=findkNN(root, &dpart[i].part, 0.5, &knn, kNN);
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
		 ang     = atan2(knn[k].part.pos[1], knn[k].part.pos[1]);
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
		 if (sigma < 0) printf("!sigmaerror!%g %g %g %g %g\n",sqrmv[j], meanv[j], SQR(meanv[j]), sigma, masstot);
	       }
	     vR2    = vR2/masstot;
	     vRsum  = vRsum/masstot;
	     vPHI2  = vPHI2/masstot;
	     vPHIsum= vPHIsum/masstot;
	     Pi_RR += vR2-SQR(vRsum);
	     Pi_PHI+= vPHI2-SQR(vPHIsum);
	   
	     //      printf("%5d %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n", i, dpart[i].dist, knndist, dpart[i].part.pos[0], dpart[i].part.pos[1], dpart[i].part.pos[2], Pi[i][0], Pi[i][1], Pi[i][2]);
	     free(knn);
	   }
	 
	 fprintf(outf, "PI %g %g %g\n", Pi[0], Pi[1], Pi[2]);
	 double mPi=(Pi[0]+Pi[1])/2;
	 delta= (mPi-Pi[2])/mPi;
	 fprintf(outf, "delta %g\n", delta);

	 if (doFOF)
	 {
	   if (debug) printf("look for central FOF-group...\n");fflush(stdout);
	   gadpart **FOF;
	   int nFOF = 0;      
	   int *done = NULL;
	   unsigned int buffer = 0;
	   fltarr center ={0.0, 0.0, 0.0};
	   
//	   findGadparts(root, center, 3 *  ll, &FOF, &nFOF, &buffer);
//	   //	   	   findGadparts(root, center,   ll, &FOF, &nFOF, &buffer);
//	   qsort(FOF, nFOF, sizeof(gadpart *), cmp_pointer_id );
//	   if (debug) printf("innermost particles determined and sorted...%d\n", nFOF);fflush(stdout);
//	   //testing
//	   
//	   if (debug) printf("sorted...\n");fflush(stdout);

	   findFOF(root, center, ll, &FOF, &nFOF, &buffer);

   	   if (debug) printf("central FOF-group determined...%d\n", nFOF);fflush(stdout);
//	   FILE* tmp;
//	   tmp = fopen("ids.tmp","w");
//	   for ( i = 0; i < nFOF; i++ )
//	     {
//	       int *fnd;
//	       fnd = bsearch( &(FOF[i] -> id), done, ndone, sizeof(int), cmp_int );
//	       if (fnd == NULL) fprintf(tmp,"%d\n", FOF[i] -> id);
//	     }
//	   fclose(tmp);
	 
	   double FOFmass = 0;
	   double FOFcenter[3] = { 0., 0., 0.};
	   for ( i = 0; i < nFOF; i++ )
	     {
	       FOFmass += FOF[i] -> mass;
	       for ( j = 0; j < 3; j++) FOFcenter[j] += FOF[i] -> pos[j] * FOF[i] -> mass;
	     }
	   double FOFdist = 0;
	   for ( j = 0; j < 3; j++) 
	     {
	       FOFcenter[j] /= FOFmass;
	       FOFdist += SQR(FOFcenter[j]);
	     }
	   FOFdist = sqrt(FOFdist);
	   fprintf(outf, "FOF-Mass: %f\n", FOFmass);
	   fprintf(outf, "FOF-center: %f %f %f\n", FOFcenter[0], FOFcenter[1], FOFcenter[2]);
	   fprintf(outf, "FOF-offset: %f\n", FOFdist);
	   
	   

	   if (use_cm&16)
	     {
	       i=0;
	       double mdum=0;
	       double fofrad=0;
	       dist=0;
	       while (dist<halorad) 
		 { 
		   dist=part[i].dist;
		   m=part_ev[part[i].ind].type;	 
		   if  (m==4) 
		     {
		       if (mdum < (FOFmass/2.0)) 
			 {
			   mdum+=part_ev[part[i].ind].mass;
			   fofrad=dist;
			 } else break;
		     }
		   i++;
		 }
	       fprintf(outf,"3D-half-mass-FOF radius: %g\n", fofrad);
	     }

	   if (doFOF & 2)
	     {
	       char FOFfile[128];
	       if (singlefile)
		 {
		   sprintf(FOFfile ,"%s.fof",prefix);
		 } else
		 {
		   sprintf(FOFfile ,"%s_%03d.fof",prefix, snapind);
		 }
	       FILE *FOFfp = fopen(FOFfile, "w");
	       for ( i = 0; i < nFOF; i++ )
		 {
		   fprintf(FOFfp, "%d\n", FOF[i] -> id);
		 }
	       fclose(FOFfp);
	     }

	   free(FOF);
	 }

	 free(root);
       }
     free(vapart);
     free(dpart);


/*************************************************************************************/     
     /*
  if (cutgad)
   {
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
	 l=wpart[i].type;
	 out.npart[l]++;
	 out.nall[l]++;
	 i++;
       }
     char rotgadfile[128];
     sprintf( rotgadfile, "%s_%03d_rot.gad", prefix, snapind);
     writegadget_part(rotgadfile, out, wpart);
   }
  */
  gsl_matrix_free(I);
  gsl_matrix_free(evec);
  gsl_matrix_free(LU);
  gsl_matrix_free(inv);
  gsl_matrix_free(resultmatrix);
  gsl_matrix_free(rotation);
  gsl_vector_free(eval);
  gsl_eigen_symmv_free(w);
  free(wpart);
   }

 /*************************************************************************************************************************************************************/
 //   endloop:
 fprintf(outf,"\n\n%g %g %g %g %g %g %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",cm[0],cm[1],cm[2],cm[0]-icm[0],cm[1]-icm[1],cm[2]-icm[2],vcnt,halorad, vmass ,par[0], par[1], halorad/par[1],(envdens/bgdens)-1,peakrotvel,peakr,rotvel,lmd, rad_ratio, q, s, Pi[0], Pi[1], Pi[2], Pi_RR, Pi_PHI);

   }

   fclose(outf);
   free(part);
   free(vpart);
   free(tr_halo_i);
   free(tr_halo_id);
   //   free(index);
#ifndef NOGAS
   free(part_ev[0].sph);
#endif
   free(part_ev);
//   free(pos_ev);
//   free(id_ev);
//   free(mass_ev);
//   free(vel);

 }
 /**********************************************************************************************************************************************/

  qsort(stars, totalstarcnt, sizeof(struct star), cmp_star_id);
  if (dostars)
    {
  sprintf(starfile ,"%s.stars.insitu",prefix);
  outf=fopen(starfile, "w");

  for (i = 0; i < totalstarcnt; i++)
    {
      fprintf(outf,"%8d %4d %8.2f %8.2f %8.7f %8.2f %4d %8.7f %1d %8.7f %8.7f %8.7f %8.7f\n", stars[i].id, stars[i].isnap, stars[i].idist, stars[i].dist, stars[i].frac, stars[i].gasdist, stars[i].snap, stars[i].a,  stars[i].fnd, stars[i].a_acc
#ifdef POTENTIAL
	      ,stars[i].pot
#else 
	      , 0.
#endif
	      , stars[i].ifrac, stars[i].ifrac_rhalf);
    }

  fclose(outf);

  sprintf(starfile ,"%s.gas.data",prefix);
  outf=fopen(starfile, "w");
  for (i = 0; i < totalgascnt; i++)
    {
      fprintf(outf,"%8d %8.4g %8.4g %8.4g %8.7f %8.7f %8.7f %8.4g %8.7f\n", gas[i].id, gas[i].maxtemp, gas[i].Tvir, gas[i].Tvir_acc, gas[i].a_maxtemp, gas[i].frac_maxtemp, gas[i].a_acc, gas[i].T_acc, gas[i].a_star);
    }
  fclose(outf);

    }
  printf("Time needed: %.2f sec\n",((float)clock())/CLOCKS_PER_SEC);
  return 0;
}
