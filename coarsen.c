/*
Program to combine particles in a gadget-file. 
Number of particles combined is dependent on the distance from a given center

icc -lm -openmp -o coarsen coarsen.c libgad.o

icc -lm -lgad-altix-nopot -openmp -o ../bin/altix/coarsen coarsen.c -DNOPOT

gcc -fopenmp -lm -lgad coarsen.c -o ~/bin/coarsen
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "libgad.h"
#define	PI 3.14159265358979323846
#define GAP 8000                                     //min distance to boarders to ignore periodic boundaries
#define h0 0.72
#define DIM 512
#define	MIN(a, b)   ((a)<(b)?(a):(b))
#define	MAX(a, b)   ((a)>(b)?(a):(b))
#define ABS(a) ((a) >= 0 ? (a) : -(a))
#define CMP(a,b) ((a)>(b)?(1):(-1))
#define PB(a,b) ((a)>(b)?(a-b):(a))
#define MOVE(a,b) PB(a+b/2,b)
#define MV(a,b) ((a)+(b)/2)%(b)
//#define SQR(x) (x)*(x)
#define G 6.6742e-11
#define Msun 1.989e30
#define kpc 3.08567758128e19

struct part
{
  fltarr pos;
  fltarr vel;
  int id;
  float mass;
};

int cmp_mass(struct part *a, struct part *b)
{
  if (a->mass > b->mass) return 1;
  else if (a->mass < b->mass) return -1;
  else return 0;
}

void usage()
{
  fprintf(stderr,"Coarsen v0.01\n");
  fprintf(stderr,"    -i <input file> -o <outputfile>\n");
  fprintf(stderr,"    -f <Coarse factor> \n");
  fprintf(stderr,"    -cm <centerX centerY centerZ>\n");
  fprintf(stderr,"    -d <min-distance to cm>\n");
  fprintf(stderr,"    -g <initial grid size>\n");
  fprintf(stderr,"    -im <initial mass (default=mass of particle 0)>\n");
  fprintf(stderr,"    -b <particles with mass greater imass*b will become bulge particles>\n");
  fprintf(stderr,"    -l <linear increase of mass with distance>\n");
  fprintf(stderr,"    -box <minx miny minz maxx maxy maxz>\n");
  fprintf(stderr,"    -all <combine all particles regardless of position>\n");
  fprintf(stderr,"    -max <maximum number of passes>\n");
  exit(1);
}

int main (int argc, char *argv[])
{
  fltarr *pos, *vel, *posnew, *velnew;
  float *mass, *massnew;
  int *id, *idnew;
  char infile[256], outfile[256];
  FILE *outf;
  struct part *partdata;
  struct header headin, headout;
  fltarr cm, pcm, vcm;
  float cf=100.0, dxgrid=100.0, fdum, keep;
  double ddum, dist, mindist, masscmb, imass=0.0, hbox, olddxgrid, omit, maxmindist;
  double a,b,normdist, maxdist, bulge=0.0, omitdist;
  fltarr min, max;
  int i,j,k,l, x,y,z, NBOX, iteration=0, go_on=1, cmball=0, maxiter=200, mincnt=0;
  int numpart, numpartnew, n,m, dum, icnt, idmin;
  int *friend, *index, cmb, pb[3], cnt, linear=0, box=0;
  int ****grid, ***gridsz;
  int debug=0;
  

  strcpy( infile,"<none>");
  strcpy(outfile,"<none>");
  i=1;
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
      else if (!strcmp(argv[i],"-f"))
	{
	  i++;
	  cf=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-d"))
	{
	  i++;
	  keep=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-max"))
	{
	  i++;
	  maxiter=atoi(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-cnt"))
	{
	  i++;
	  mincnt=atoi(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-all"))
	{
	  i++;
	  cmball=1;
	}
      else if (!strcmp(argv[i],"-debug"))
	{
	  i++;
	  debug=1;
	}
      else if (!strcmp(argv[i],"-im"))
	{
	  i++;
	  imass=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-b"))
	{
	  i++;
	  bulge=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-g"))
	{
	  i++;
	  dxgrid=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-l"))
	{
	  i++;
	  linear=1;
	}
      else if (!strcmp(argv[i],"-cm"))
	{
	  i++;
	  cm[0]=atof(argv[i++]);
	  cm[1]=atof(argv[i++]);
	  cm[2]=atof(argv[i++]);
	}
      else if (!strcmp(argv[i],"-box"))
	{
	  box=1;
	  i++;
	  min[0]=atof(argv[i++]);
	  min[1]=atof(argv[i++]);
	  min[2]=atof(argv[i++]);
	  max[0]=atof(argv[i++]);
	  max[1]=atof(argv[i++]);
	  max[2]=atof(argv[i++]);
	}
      else usage();
    }
  if (!strcmp(infile, outfile)) usage();



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
      fflush(stdout);

  numpart=readgadget(infile, &headin, &pos, &vel, &id, &mass);
  if (0==imass) imass=mass[0];
  hbox=headin.boxsize;
  double BOXSIZE=headin.boxsize;
  normdist=sqrt(3*SQR(hbox/2.0));
  omitdist=keep;
  if (box)
    {
      omitdist=hbox;
      normdist=0;
      for (i=0; i<3; i++)
	{
	  if (min[i]<max[i])
	    {
	      normdist+=SQR((BOXSIZE-(max[i]-min[i]))/2.0);
	      omitdist= (omitdist>(max[i]-min[i])) ? (max[i]-min[i]) : omitdist;
	    }
	  else
	    {
	      normdist+=SQR((min[i]-max[i])/2.0);
	      ddum=BOXSIZE-(min[i]-max[i]);
	      omitdist= (omitdist>ddum) ? (ddum) : omitdist;
	    }
	}
      normdist=sqrt(normdist);
      omitdist=omitdist/2.0;
      printf("normdist %g\n", normdist);fflush(stdout);
    }

  a=( (2*normdist - keep*cf) / (normdist*SQR(keep) - SQR(normdist)*keep) );
  b=( (2-a*keep*keep) / (keep) );
  if (linear)
    {
      a=(cf-2)/(normdist - keep);
      b=cf-a*normdist;
    }
  printf("Finished reading input file.\n");fflush(stdout);
  while (go_on)
    {
      iteration++;
      friend = (int *) calloc(numpart,sizeof(int));
      olddxgrid=dxgrid;
      dxgrid*=1.10;
      n=ceil(headin.boxsize/dxgrid);
      if (n<8) 
	{
	  n=8;
	  dxgrid=olddxgrid;
	}
      m=ceil(numpart/pow(n,3))*2;
      printf("Building Grid (dim: %d)\n", n);
      grid=  (int ****) malloc (n*sizeof(int ***));                       //allocate 4-dimensional array to store particle indexes
      gridsz=(int ***)  malloc (n*sizeof(int **));
      for (i=0; i < n; i++) 
	{
	  grid[i]=  (int ***) malloc (n * sizeof(int **));
	  gridsz[i]=(int **)  malloc (n * sizeof(int *));
	}
      for (i=0; i < n; i++) 
	for (j=0; j < n; j++)
	  {
	    grid[i][j]=   (int **) malloc (n * sizeof(int *));
	    gridsz[i][j]= (int *)  malloc (n * sizeof(int ));
	  }
      for (i=0; i < n; i++) 
	for (j=0; j < n; j++)
	  for (k=0; k < n; k++)
	    {
	      grid[i][j][k]= (int *) malloc (m * sizeof(int));
	      if (grid[i][j][k]==NULL)
		{
		  fprintf(stderr,"Grid memory allocation failed!\n");
		  exit(3);
		}
	      grid[i][j][k][0]=0;
	      gridsz[i][j][k]=m;
	    }

      j=0;
      //Populating the grid
      for (x=0; x< numpart; x++) 
	{
	  i=floor(pos[x][0]/dxgrid);
	  j=floor(pos[x][1]/dxgrid);
	  k=floor(pos[x][2]/dxgrid);
	  if (gridsz[i][j][k]<=(grid[i][j][k][0]+2)) 
	    {
	      grid[i][j][k]=realloc(grid[i][j][k],sizeof(int)*(gridsz[i][j][k]+m));
	      gridsz[i][j][k]+=m;
	    }
	  grid[i][j][k][0]++;
	  grid[i][j][k][grid[i][j][k][0]]=x;
	} 
      dum=0;
      for (i=0; i < n; i++) 
	for (j=0; j < n; j++)
	  for (k=0; k < n; k++) 
	    {
	      dum+=grid[i][j][k][0];
	      grid[i][j][k]=realloc(grid[i][j][k],sizeof(int)*(grid[i][j][k][0]+1));
	    }
      if (dum!=numpart)
	{
	  fprintf(stderr, "Something went wrong, Grid error!\n");
	  exit(2);
	}


      NBOX=1;
      printf("\nIteration: %d Gridsize: %g\n", iteration, dxgrid);
      cnt=0;
      maxmindist=0;
#pragma omp parallel for firstprivate(y,z, i,j,k,l, omit, index, dist, icnt, idmin, mindist, pcm, maxmindist) reduction (+ : cnt)
      for (x=0; x < n; x++) 
	{
	  index=(int *) malloc(numpart*sizeof(int));
	  for (y=0; y < n; y++)
	    for (z=0; z < n; z++) 
	      {
		if ((iteration > 3) && (!cmball))
		  {
		    pcm[0]=(x+0.5)*dxgrid;
		    pcm[1]=(y+0.5)*dxgrid;
		    pcm[2]=(z+0.5)*dxgrid;
		    dist=distance(pcm, cm);
		    omit=omitdist+((sqrt(iteration)/2.0)*dxgrid);
		    if (dist<omit)
		      {
			cnt++;
			continue;
		      }
		  }
		icnt=0;
		
		for (i=(x-NBOX); i <= (x+NBOX); i++)
		  for (j=(y-NBOX); j <= (y+NBOX); j++)
		    for (k=(z-NBOX); k <= (z+NBOX); k++)
		      {
			if ((iteration > 3) && (!cmball))
			  {
			    pcm[0]=(i+0.5)*dxgrid;
			    pcm[1]=(j+0.5)*dxgrid;
			    pcm[2]=(k+0.5)*dxgrid;
			    dist=distance(pcm, cm);
			    if (dist<omit)
			      {
				continue;
			      }
			  }
			for (l=1; l <= grid[(i+n)%n][(j+n)%n][(k+n)%n][0]; l++)
			  {
			    index[icnt++]=grid[(i+n)%n][(j+n)%n][(k+n)%n][l];
			  }
		      }

		for (i=1; i<= grid[x][y][z][0]; i++)
		  {
		    j=grid[x][y][z][i];
		    if (j!=index[0])
		      {
			idmin=index[0];
		      } else {
			idmin=index[1];
		      }
		    mindist=distance(pos[j], pos[idmin]);
		    for (l=0; l< icnt; l++)
		      {
			k=index[l];
			if (j!=k)
			  {
			    dist=distance(pos[j],pos[k]);
			    if (dist<mindist) 
			      {
				idmin=k;
				mindist=dist;
			      }
			  }
		      }
		    friend[j]=idmin;
		    if (mindist>maxmindist) 
		      {
			  maxmindist=mindist;
		      }
		  }
	      }
      	  if ((x)==0) printf("icnt %d maxmindist %g\n", icnt, maxmindist);fflush(stdout);
	  free(index);
	}
      // End of parallel region
      printf("gridcells omitted %d of %d\n", cnt, n*n*n);
      printf("Friends found\n");fflush(stdout);

      dum=0;
      numpartnew=0;
      posnew =(fltarr *)malloc(sizeof(fltarr)*numpart);
      velnew =(fltarr *)malloc(sizeof(fltarr)*numpart);
      idnew  =(int *)   malloc(sizeof(int)*numpart);
      massnew=(float *)malloc(sizeof(float)*numpart);
      if (massnew == NULL) 
	{
	  fprintf(stderr, "malloc failed!\n");
	  exit(1);
	}
      maxdist=0;
      for (i=0; i<numpart; i++) 
	{
	  cmb=0;
	  if ((friend[i] >= 0) && (i==friend[friend[i]]) && (i!=friend[i]))
	    {
	      for (j=0; j<3; j++) {pcm[j]=0;vcm[j]=0;pb[j]=0;}
	      for (j=0; j<3; j++) 
		{
		  if (ABS(pos[i][j]-pos[friend[i]][j]) < (hbox/2)) 
		    {
		      pcm[j]+=pos[i][j]*mass[i];
		      pcm[j]+=pos[friend[i]][j]*mass[friend[i]];
		    }
		  else
		    {
		      pb[j]=1;
		      pcm[j]+=MOVE(pos[i][j], hbox)*mass[i];
		      pcm[j]+=MOVE(pos[friend[i]][j], hbox)*mass[friend[i]];
		    }
		  vcm[j]+=vel[i][j]*mass[i];
		  vcm[j]+=vel[friend[i]][j]*mass[friend[i]];
		}
	      masscmb=mass[i]+mass[friend[i]];
	      for (j=0; j<3; j++) 
		{
		  if (pb[j])
		    pcm[j]=MOVE(pcm[j]/masscmb, hbox);
		  else pcm[j]=pcm[j]/masscmb;
		  vcm[j]=vcm[j]/masscmb;
		}
	      if (cmball) cmb=1;
	      else
		{
		  if (box) dist=distbox(pcm, min, max);
		  else dist=distance(pcm, cm);
		  if (dist>maxdist) maxdist=dist;
		  if (linear) {
		    if ((dist > keep) && ( (a*dist+b)  > (masscmb/imass))) cmb=1;
		  } else if ((dist > keep) && ( (a*SQR(dist)+b*dist)  > (masscmb/imass))) cmb=1;
		}
	    }
	  if (cmb)
	    {
	      for (j=0; j<3; j++) 
		{
		  posnew[numpartnew][j]=pcm[j];
		  velnew[numpartnew][j]=vcm[j];
		}
	      massnew[numpartnew]=masscmb;
	      idnew[numpartnew]=id[i];
	      friend[friend[i]]=-1;
	      friend[i]=-1;
	      numpartnew++;
	    } else if (friend[i]>=0) {
	      for (j=0; j<3; j++) 
		{
		  posnew[numpartnew][j]=pos[i][j];
		  velnew[numpartnew][j]=vel[i][j];
		}
	      massnew[numpartnew]=mass[i];
	      idnew[numpartnew]=id[i];
	      numpartnew++;
	    }
	}
      printf("# of particles: old %d new %d\n", numpart ,numpartnew);
      printf("maxdist: %g\n", maxdist);
      free(friend);
      for (i=0; i < n; i++) 
	for (j=0; j < n; j++)
	  for (k=0; k < n; k++)
	    {
	      free(grid[i][j][k]);
	    }
      for (i=0; i < n; i++) 
	for (j=0; j < n; j++)
	  {
	    free(grid[i][j]);
	    free(gridsz[i][j]);
	  }
      for (i=0; i < n; i++) 
	{
	  free(grid[i]);
	  free(gridsz[i]);
	}
      free(grid);
      free(gridsz);
      if (numpartnew==numpart) go_on=0;
      if (numpartnew< mincnt) go_on=0;
      if (iteration >= maxiter ) go_on=0;
      if (iteration==200) go_on=0;
      //	  printf("%10.2f != %g\n", pos[12][2], posnew[12][2]);
      //	  printf("%10.2f != %g\n", vel[13][1], velnew[13][1]);
      //	  printf("%10.2f != %g\n", mass[133], massnew[133]);
      memcpy(pos[0], posnew[0], sizeof(fltarr)*numpartnew);
      memcpy(vel[0], velnew[0], sizeof(fltarr)*numpartnew);
      memcpy(id, idnew, sizeof(int)*numpartnew);
      memcpy(mass, massnew, sizeof(float)*numpartnew);
      pos = realloc(pos, sizeof(fltarr)*numpartnew);
      vel = realloc(vel, sizeof(fltarr)*numpartnew);
      id  = realloc(id, sizeof(int)*numpartnew);
      mass= realloc(mass, sizeof(float)*numpartnew);
      //	  printf("%10.2f ?= %g\n", pos[12][2], posnew[12][2]);
      //	  printf("%10.2f ?= %g\n", vel[13][1], velnew[13][1]);
      //	  printf("%10.2f ?= %g\n", mass[133], massnew[133]);
      numpart=numpartnew;
      free(posnew);
      free(velnew);
      free(idnew);
      free(massnew);
    }

  posnew  =(fltarr *)malloc(sizeof(fltarr)*numpart);
  velnew  =(fltarr *)malloc(sizeof(fltarr)*numpart);
  idnew   =(int *)   malloc(sizeof(int)*numpart);
  massnew =(float *)malloc(sizeof(float)*numpart);
  partdata=(struct part *)malloc(sizeof(struct part)*numpart);

  headout=headin;
  for (j=0; j<6; j++)
    {
      headout.npart[j]=0;
      headout.nall[j]=0;
    }
  
  for (i=0; i<numpart; i++)
    {
      for (j=0; j<3; j++)
	{
	  partdata[i].pos[j]=pos[i][j];
	  partdata[i].vel[j]=vel[i][j];
	}
      partdata[i].id=id[i];
      partdata[i].mass=mass[i];
    }

  qsort(&partdata[0], numpart, sizeof(struct part),(void *)cmp_mass);

  for (i=0; i<numpart; i++)
    {
      for (j=0; j<3; j++)
	{
	  posnew[i][j] =partdata[i].pos[j];
	  velnew[i][j] =partdata[i].vel[j];
	}
      massnew[i]=partdata[i].mass;
      idnew[i]  =partdata[i].id;
//      if (massnew[i]==imass)
//	{
//	  headout.npart[1]++;
//	  headout.nall[1]++;	  
//	} else
      if ((0==bulge) || (massnew[i] <= (bulge*imass)))
	{
	  headout.npart[2]++;
	  headout.nall[2]++;
	} else {
	  headout.npart[3]++;
	  headout.nall[3]++;
	}
    }

  writegadget(outfile, headout, posnew, velnew, idnew, massnew);
  return 0;
}


