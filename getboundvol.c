/*
Program to compute the bounding volume of a set of particles


*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "libgad.h"

#define	PI 3.14159265358979323846
#define h0 0.72
#define DIM 512
#define LARGE_WS 8
#define kNN 50
#define	MIN(a, b)   ((a)<(b)?(a):(b))
#define	MAX(a, b)   ((a)>(b)?(a):(b))
#define ABS(a) ((a) >= 0 ? (a) : -(a))
#define CMP(a,b) ((a)>(b)?(1):(-1))
#define SOFTENING 5.00
#define G 6.6742e-11
#define Msun 1.989e30
#define kpc 3.08567758128e19

double cdens;
const int MAXPART=10000000;

struct particle {int ind;float dist;};
clock_t t[2];

// double fit_fct(double t, double *p)
// {
//   return log10((p[0]) / ((t/p[1]) * SQR(1+t/p[1])));
// }

// void printout ( int n_par, double* par, int m_dat, double* fvec, 
//                        void *data, int iflag, int iter, int nfev )
//          {
//            // dummy function to catch fitting output

//          }


void usage()
{
  fprintf(stderr,"get Bounding Volume v0.01\n");
  fprintf(stderr,"-f <snapshot basefilename>\n");
  //  fprintf(stderr,"-cut <create Gadget-file for every halo>\n");
  fprintf(stderr,"-id <file with particle IDs to trace> otherwise stdin\n");
  fprintf(stderr,"-gd <id-file is a gadget file>\n");
  fprintf(stderr,"-x <increase by factor (shrink for < 1)>\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n");
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

// int cmp_int (const void *first, const void *second)
// {
//   int *a= (int *) first; 
//   int *b= (int *) second; 
//   if (*a>*b) return 1;
//   else if (*a<*b) return -1;
//   else return 0;
// }


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
  //  float timestep[NUMTIMESTEPS];


  char basename[256], prefix[256], idfile[256], snapfile[256], outfile[256];
  int debug=0, test=0, gadget=0, single=0;
  double dx512=72000.0/512.0;
  double cm[3]={0.,0.,0.};
  int startind, endind, snapind;
  int use = 63; //All particle types
  int i,j,k;

  int *cid, *insideID;
  FILE *fp;
  double inc_factor = 1.0;

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //filenames are ignored if given as command line parameters 
  strcpy(basename,"");
  sprintf(prefix,"boundIDs");
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
      else if (!strcmp(argv[i],"-id"))
	{
	  i++;
	  strcpy(idfile,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-use"))
	{
	  i++;
	  use = atoi(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-x"))
	{
	  i++;
	  inc_factor = atof(argv[i]);
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
      else if (!strcmp(argv[i],"-gd"))
	{
	  i++;
	  gadget=1;
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


  if (debug) printf("comencing main loop...\n");fflush(stdout);

/************************************************************************************/
/* Read in Particle Id-List                                                         */
/************************************************************************************/
  unsigned int numtrace =0;
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

  printf("Sorting %d entries...\n", numtrace);
  qsort(cid, numtrace, sizeof(int), cmp_int);

/************************************************************************************/
/* Start of main loop                                                               */
/************************************************************************************/

  struct gadpart *wpart=(struct gadpart*)malloc(numtrace * sizeof(struct gadpart));
  struct gadpart *part;
  unsigned int numpart;
  struct header evolved;

  insideID=(int *) calloc(MAXPART, sizeof(int));
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
      else printf("number of particles in snapshot: %d\n",numpart);

//     out=evolved;
//     for (i=0; i<6; i++)
//       {
//	 out.npart[i]=0;
//	 out.nall[i]=0;
//       }
      printf("SnapInd: %d\n",snapind);
      k=0;
      for (i=0; i< numpart; i++)
	{
	  int *fnd=bsearch(&(part[i].id),cid, numtrace, sizeof(int), cmp_int);
	  if (fnd!=NULL)
	    {
	      if ( (1<<part[i].type) & use)
		{
		  cpygadpart(&wpart[k], &part[i]);
		  //		      out.npart[wpart[k].type]++;
		  //		      out.nall[wpart[k].type]++;
		  k++;
		}
	    }
	  if (k==numtrace) break;
	}
      unsigned int numfound=k;
      if (debug) {printf("found %d particles in snapshot\n", numfound);fflush(stdout);}
      double rad=1.e6;
      double res[2];
      simplecenter(wpart, numfound, cm, use);
      for (i = 0; i < numfound; i++)
	{
	  for (j = 0; j < 3; j++)
	    wpart[i].pos[j] -= cm[j];
	}
      if (debug) {printf("CM %g %g %g\n", cm[0],cm[1],cm[2]);fflush(stdout);}
      gsl_matrix *rotation;
      rotategalaxy(wpart, numfound, rad, use, res, &rotation);
      if (debug) {printf("rotation %g %g\n", res[0], res[1]);fflush(stdout);}
      for (i = 0; i < numpart; i++)
	{
	  for (j = 0; j < 3; j++)
	    part[i].pos[j] -= cm[j];
	}
      rotatepart(part, numpart, rotation);
      fltarr min ={0.,0.,0.};
      fltarr max ={0.,0.,0.};
      double deform[3] = {1.,1.,1.};
      double correction[3] = {0., 0., 0.};

      for (i = 0; i < numfound; i++)
	{
	  for (j = 0; j < 3; j++)
	    {
	      if (wpart[i].pos[j] < min[j]) min[j] = wpart[i].pos[j];
	      if (wpart[i].pos[j] > max[j]) max[j] = wpart[i].pos[j];
	    }	    
	}
      
      for (j = 0; j < 3; j++)
	{	  
	  correction[j] = (max[j] + min[j]) / 2.;
	  deform[j] = 1./(max[j] - correction[j]);
//	  if (ABS(min[j]) > ABS(max[j])) deform[j] = 1./ABS(min[j]);
//	  else deform[j] = 1./ABS(max[j]);
//	  if (debug) printf("diff: %g\n", ABS(max[j]) - ABS(min[j]));
	}
      double maxdist = 0;
      for (i = 0; i < numfound; i++)
	{
	  double dist = 0;
	  for (j = 0; j < 3; j++)
	    {
	      dist += SQR((wpart[i].pos[j] - correction[j]) * deform[j]);
	    }
	  dist = sqrt(dist);
	  if (dist > maxdist) maxdist = dist;
	}
      if (debug) {printf("deform: %g | %g | %g | maxdist %g\n", deform[0], deform[1], deform[2], maxdist);fflush(stdout);}
      unsigned int numinside=0;
      for (i = 0; i < numpart; i++)
	{
	  if (! ( (1<<part[i].type) & use) )
	    continue;
	  double dist = 0;
	  for (j = 0; j < 3; j++)
	    {
	      dist += SQR( (part[i].pos[j] - correction[j]) * deform[j] / inc_factor);
	    }
	  dist = sqrt(dist);
	  if (dist <= maxdist) 
	    {
	      insideID[numinside++]= part[i].id;
	      //	      if (test) printf("%d\n", part[i].id);
	    } 

	}
      if (debug) {printf("found %d particles inside ellipsoid\n", numinside);fflush(stdout);}

      if (single)
	sprintf(outfile,"%s.dat",prefix);
      else
	sprintf(outfile,"%s_%03d.dat",prefix, snapind);
      
      FILE *outf = fopen(outfile, "w");
      fwrite(&numinside, sizeof(unsigned int), 1, outf);
      fwrite(&insideID[0], sizeof(int), numinside, outf);
      fclose(outf);
      if (debug) {printf("finished writing %d IDs\n", numinside);fflush(stdout);}
  
//      writegadget_part(gadfile, out, wpart);

//      if (debug) {printf("Gadget File written\n");fflush(stdout);}



 /**********************************************************************************************************************************************/
      free(part);
 }
 /**********************************************************************************************************************************************/



  return 0;
}
