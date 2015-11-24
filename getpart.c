/*

(contains code overhead, since it is based on a copy of trhalomult.c)

#icc -lm -lgsl -lgslcblas -o ../bin/altix/trace trace.c lmfit-2.2/lmmin.o lmfit-2.2/lm_eval.o libgad.o KdTree.o -I/usr/local/gsl-1.6

icc -lm -lgsl -lgslcblas -o ../bin/altix/getpart getpart.c -lgad-altix KdTree.o -I/usr/local/gsl-1.6


APPLE
gcc -lm -lgsl -lgslcblas -lgad-stan -o ~/bin/getpart getpart.c stan/KdTree.o 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "libgad.h"
#include "KdTree.h"

#define	PI 3.14159265358979323846
#define h0 0.72
#define DIM 512
#define ENVDENSRAD 2000
#define INERTIA_START_RAD 0.4
#define INERTIA_USE_PART_TYPE 2
#define CM_USE_PART_TYPE 2
#define VA_USE_PART_TYPE 2
#define MIN_PARTICLE_COUNT 50
#define GRIDSIZE 1000
#define WORKSPACE 0.1                               //workspace for WORKSPACE*numpart particles is allocated, if glibc or segmentation fault occur, try to compile with a larger fraction
#define LARGE_WS 8
#define kNN 50
#define	MIN(a, b)   ((a)<(b)?(a):(b))
#define	MAX(a, b)   ((a)>(b)?(a):(b))
#define ABS(a) ((a) >= 0 ? (a) : -(a))
#define CMP(a,b) ((a)>(b)?(1):(-1))
#define PB(a,b) ((a)>(b)?(a-b):(a))
#define MOVE(a,b) PB(a+b/2,b)
#define MOVEB(a) MOVE((a),boxsz)
#define MV(a,b) ((a)+(b)/2)%(b)
#define SOFTENING 5.00
#define G 6.6742e-11
#define Msun 1.989e30
#define kpc 3.08567758128e19
#define MAXNHALOES 30000
#define NUMTIMESTEPS 94
#define NUMTRACE 50
#define MAXPART 10000000

double cdens;

struct particle {int ind;float dist;};
clock_t t[2];

double fit_fct(double t, double *p)
{
  return log10((p[0]) / ((t/p[1]) * SQR(1+t/p[1])));
}

void usage()
{
  fprintf(stderr,"GetPart v0.01\n");
  fprintf(stderr,"getpart <gadgetfile> [options] \n");
  //  fprintf(stderr,"-i <startindex> <endindex>\n");
  //  fprintf(stderr,"-cut <create Gadget-file for every halo>\n");
  fprintf(stderr,"-trid <file with particle IDs to trace> or STDIN\n");
  fprintf(stderr,"-gd <trid is a gadget file>\n");
  fprintf(stderr,"");
  fprintf(stderr,"");
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
  struct header ic, evolved, out;
  FILE *fp, *outf;
  char cmfile[80],icfile[80],snapfile[256], basename[256],snfile[80],friendsfile[80],idlfile[80],gridfile[80],jrfile[80], outfile[80], gadfile[80], prefix[20], idfile[80];
  char output[256][5];
  unsigned int numpart=0, nummass=0;
  int blocksize,i,j,k,l,m,n,dum,tr_halo=0,tr_halo_cnt = 0,id,halo,checknpart,mindum,nhalo,count,notrace=0;
  int *tr_halo_id,*tr_halo_i, ca[3], size[3];
  fltarr *pos_ev,*pos_ic,*vel, *outpos, *outvel;
  float *mass_ev, *mass_dum,*mass_ic,dist,maxdist,*distance,halomass,halorad, halolmd, *outmass;
  double posdum,max[3],min[3],srad=0, icmin[3], icmax[3], envdens, bgdens;
  float lmd,vcm[3]={0,0,0},jcm[3]={0,0,0},J[3]={0,0,0}, torq[3],rdum[3], vdum[3], massdum;
  int *id_ev,*id_ic,*iclus, halonpart, *outid;
  float tm,linkle,od;
  double boxsz,cm[3]={0,0,0}, masstot, gridsize, rlog[50], p[50], err[50], d, rad_ratio;
  short pb[3]={0,0,0};
  int minbox[3], maxbox[3], szbox[3], ****grid, ***gridsz, boxind[3], alignf;
  int nfun, npar, dojr=0, vcnt=0, use_inertia, use_va, use_cm, use=0;
  double par[2],ddum, dx512, ddum1, totj, rotvel, peakrotvel, peakr, vmass, sqrenvdensrad;
  int  cutgad=0, idum, donotoverwrite=0, test=0, gadget=0, *cid, numalloc, numtrace;
  int central_galaxy = 0;
  fltarr *hcm;
  int debug=0, maxhaloes, minpartcnt;
  int startind, endind, snapind;
  int single=0;
  float ll = 0.8;
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
  sprintf(prefix,"tr");
  int pid=0, psa=0, psfr=0, ptemp=0, prho=0, pmass=0, ppos=0;
  int numfields=0;
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //filenames are ignored if given as command line parameters 

  strcpy(gridfile,"<none>");
  strcpy(cmfile,"<none>");
  strcpy(basename,"");
  strcpy(icfile,"");
  strcpy(idlfile,"idlinp.dat");
  //  tr_halo=123;
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  i=1;
  j=0;
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
          single=1;
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
      else if (!strcmp(argv[i],"-trid"))
	{
	  i++;
	  strcpy(idfile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-debug"))
	{
	  i++;
	  debug=1;
	}
      else if (!strcmp(argv[i],"-id"))
	{
	  i++;
	  pid=1;
	  strcpy(output[j++],"ID");
	}
      else if (!strcmp(argv[i],"-cg"))
	{
	  i++;
	  central_galaxy=1;
	}

      else if (!strcmp(argv[i],"-sa"))
	{
	  i++;
	  psa=1;
	  strcpy(output[j++],"SA");
	}
      else if (!strcmp(argv[i],"-sfr"))
	{
	  i++;
	  psfr=1;
	  strcpy(output[j++],"SFR");
	}
      else if (!strcmp(argv[i],"-pos"))
	{
	  i++;
	  ppos=1;
	  strcpy(output[j++],"POS");
	}
      else if (!strcmp(argv[i],"-pos_phys"))
	{
	  i++;
	  ppos=1;
	  strcpy(output[j++],"PPOS");
	}
      else if (!strcmp(argv[i],"-vel"))
	{
	  i++;
	  strcpy(output[j++],"VEL");
	}
      else if (!strcmp(argv[i],"-vel_phys"))
	{
	  i++;
	  strcpy(output[j++],"PVEL");
	}
      else if (!strcmp(argv[i],"-temp"))
	{
	  i++;
	  ptemp=1;
	  strcpy(output[j++],"TEMP");
	}
      else if (!strcmp(argv[i],"-rho"))
	{
	  i++;
	  prho=1;
	  strcpy(output[j++],"RHO");
	}
      else if (!strcmp(argv[i],"-mass"))
	{
	  i++;
	  pmass=1;
	  strcpy(output[j++],"MASS");
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
      else if (!strcmp(argv[i],"-gd"))
	{
	  i++;
	  gadget=1;
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
      else if (!strcmp(argv[i],"-ui"))
	{
	  i++;
	  use_inertia=atoi(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-use"))
	{
	  i++;
	  use=atoi(argv[i]);
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
      else if (!strcmp(argv[i],"-ll"))
	{
	  i++;
	  ll=atof(argv[i]);
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
  numfields=j;



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
/* Read in Particle Id-List                                                         */
/************************************************************************************/
  if (gadget)
    {
      fltarr *d1, *d2;
      float *d3;
      struct header rhead;
      numtrace=readgadget(idfile, &rhead, &d1, &d2, &cid, &d3);
      free(d1);
      free(d2);
      free(d3);
    }
  else
    {
      cid=(int *) calloc(MAXPART, sizeof(int));
      fp=fopen(idfile, "r");
      if (fp==NULL) fp=stdin;
      i=0;
      while (fscanf(fp,"%d", &cid[i++])==1)
	{
	  //      printf("%d\n", cid[i-1]);
	  if (feof(fp)) break;
	}
      numtrace=i;
      cid=realloc(cid, numtrace * sizeof(int));
      fclose(fp);
    }

  if (debug) printf("Sorting %d entries...\n", numtrace);
  qsort(cid, numtrace, sizeof(int), cmp_int);

/************************************************************************************/
/* Start of main loop                                                               */
/************************************************************************************/

//  struct gadpart *wpart=(struct gadpart*)malloc(numtrace * sizeof(struct gadpart));
  struct gadpart *part;
  gadpart * wpart;

  for (snapind=endind; snapind>=startind; snapind--)
    { 
      if (single)
	sprintf(snapfile,"%s", basename);
      else
	sprintf(snapfile,"%s%03d", basename, snapind);
      /*
      if (access(outfile, 0)!=0) {
	printf("File %s does not exits!\n", snapfile);
	continue;
      }
      */
     
      //      numpart=readgadget(snapfile, &evolved, &pos_ev, &vel, &id_ev, &mass_ev);
      numpart=readgadget_part(snapfile, &evolved, &part);
      if (numpart==0) continue;
      else if (debug) printf("number of particles in snapshot: %d\n",numpart);

      gadpart **galpart;
      int ngal = 0;      
      if (central_galaxy)
	{
	  int nstars = evolved.npart[4];
	  int starind = evolved.npart[0] + evolved.npart[1] + evolved.npart[2] + evolved.npart[3];

	  KdNode *all_root;
	  initKdNode(&all_root, NULL);
	  buildKdTree(all_root, &part[starind], nstars, 0);
	  unsigned int buffer = 0;
	  fltarr center ={0.0, 0.0, 0.0};
	  
	  findGadparts(all_root, center, 3 *  ll, &galpart, &ngal, &buffer);
	  //	  qsort(galpart, ngal, sizeof(gadpart *),(__compar_fn_t) cmp_pointer_id );
	  qsort(galpart, ngal, sizeof(gadpart *), cmp_pointer_id );

	  findFOF(all_root, center, ll, &galpart, &ngal, &buffer);
	  
	  delKdNode(&all_root);
	}

      if (cutgad) 
	{
	  wpart = (gadpart * ) malloc (sizeof(gadpart) * numpart);
	}
     out=evolved;
     for (i=0; i<6; i++)
       {
	 out.npart[i]=0;
	 out.nall[i]=0;
       }
     if (debug) printf("SnapInd: %d\n",snapind);
     k=0;
     //     char field[10][30];
        for (i=0; i< numpart; i++)
	    {
	      int *fnd=bsearch(&(part[i].id),cid, numtrace, sizeof(int), cmp_int);
	      if (fnd!=NULL)
		{
//		  cpygadpart(&wpart[k], &part[i]);
//		  out.npart[wpart[k].type]++;
//		  out.nall[wpart[k].type]++;
//		  if (pid)
//		    sprintf(field[0],"%d ", part[i].id);
//		  if (psa)
//		    sprintf(field[1]," %g ", part[i].stellarage);
//		  if (psfr)
//		    sprintf(field[2]," %g ", part[i].sph->sfr);
//		  if (ptemp)
//		    sprintf(field[3]," %g ", temperature(part[i]));
//		  if (prho)
//		    sprintf(field[4]," %g ", part[i].sph->rho);
//		  if (pmass)
//		    sprintf(field[5]," %g ", part[i].mass);


		  if (cutgad) 
		    {
		      out.npart[part[i].type]++;
		      out.nall [part[i].type]++;
		      wpart[k] = part[i];
		    }
		  k++;
		
		  for ( j = 0; j < numfields; j++)
		    {
		      if (!strcmp(output[j],"ID"))
			printf(" %d ", part[i].id);
		      else if (!strcmp(output[j],"SA"))
			printf(" %g ", part[i].stellarage);
		      else if (!strcmp(output[j],"SFR"))
			printf(" %g ", part[i].sph->sfr);
		      else if (!strcmp(output[j],"TEMP"))
			printf(" %g ", temperature(part[i]));
		      else if (!strcmp(output[j],"RHO"))
			printf(" %g ", part[i].sph->rho);
		      else if (!strcmp(output[j],"MASS"))
			printf(" %g ", part[i].mass);
		      else if (!strcmp(output[j],"POS"))
			printf(" %g %g %g ",part[i].pos[0], part[i].pos[1], part[i].pos[2]);
		      else if (!strcmp(output[j],"PPOS"))
			printf(" %g %g %g ",part[i].pos[0] * evolved.time / evolved.hubparam, part[i].pos[1] * evolved.time / evolved.hubparam, part[i].pos[2] * evolved.time / evolved.hubparam);
		      else if (!strcmp(output[j],"VEL"))
			printf(" %g %g %g ",part[i].vel[0], part[i].vel[1], part[i].vel[2]);
		      else if (!strcmp(output[j],"PVEL"))
			printf(" %g %g %g ",part[i].vel[0] * sqrt(evolved.time), part[i].vel[1] * sqrt(evolved.time), part[i].vel[2] * sqrt(evolved.time));
		 
		    }
		  printf("\n");
		}

	      if (k==numtrace) break;
	    }
	if (debug) printf("%d\n",k);
	  //	  if (debug) printf("CID %d\n", cid);
	  //set start CM
	
	  if (single)
	    sprintf(gadfile,"%s.gad",prefix);
	  else
	    sprintf(gadfile,"%s_%03d.gad",prefix, snapind);

      if (debug)
	{
	  for (i=0; i<6; i++)
	    {
	      printf("numpart %d: %d\n",i,out.npart[i]);
	      //	      out.nall[i]=0;
	    }
	}

      if (cutgad) 
	writegadget_part(gadfile, out, wpart);

      if (debug) {printf("Gadget File written\n");fflush(stdout);}



 /**********************************************************************************************************************************************/

 }
 /**********************************************************************************************************************************************/



  if (debug) printf("Time needed: %.2f sec\n",((float)clock())/CLOCKS_PER_SEC);
  return 0;
}
