#define SKIP fread(&blocksize,sizeof(int),1,fp);
#define SKIPFORMAT2 if (snapformat == 2) fseek(fp, 16, SEEK_CUR);	\
  fread(&blocksize,sizeof(int),1,fp);
#define SKIP2 fread(&blocksize2,sizeof(int),1,fp);
#define BLOCK fwrite(&blocksize,sizeof(int),1,fp);
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "libgad.h"
//#include "./lmfit-2.2/lmmin.h"
//#include "./lmfit-2.2/lm_eval.h"
#include "./lmmin.h"
#include "./lm_eval.h"


float BOXSIZE;
int basic=0;
int libgaderr=0;
int search_dim=0;
double crit_dens = (2.78e-8);

int cmp_id (const void *first, const void *second)
{
  struct gadpart *a = (struct gadpart  *)first;
  struct gadpart *b = (struct gadpart  *)second;
  if (a->id > b->id) return 1;
  else if (a->id < b->id) return -1;
  else return 0;
}

int cmp_pointer_id(const void *a, const void *b)
{
  struct gadpart **x= (struct gadpart**)a;
  struct gadpart **y= (struct gadpart**)b;
  if ((*x)->id > (*y)->id) return 1;
  if ((*x)->id < (*y)->id) return -1;
  return 0;
}

int cmp_type (const void *first, const void *second)
{
  struct gadpart *a = (struct gadpart  *)first;
  struct gadpart *b = (struct gadpart  *)second;
  if (a->type > b->type) return 1;
  else if (a->type < b->type) return -1;
  else return 0;
}

int cmp_pos (const void *first, const void *second)
{
  struct gadpart *a = (struct gadpart  *)first;
  struct gadpart *b = (struct gadpart  *)second;
  if (a->pos[search_dim] > b->pos[search_dim]) return 1;
  else if (a->pos[search_dim] < b->pos[search_dim]) return -1;
  else return 0;
}

int cmp_x (const void *first, const void *second)
{
  struct gadpart *a = (struct gadpart  *)first;
  struct gadpart *b = (struct gadpart  *)second;
  if (a->pos[0] > b->pos[0]) return 1;
  else if (a->pos[0] < b->pos[0]) return -1;
  else return 0;
}
int cmp_y (const void *first, const void *second)
{
  struct gadpart *a = (struct gadpart  *)first;
  struct gadpart *b = (struct gadpart  *)second;
  if (a->pos[1] > b->pos[1]) return 1;
  else if (a->pos[1] < b->pos[1]) return -1;
  else return 0;
}
int cmp_z (const void *first, const void *second)
{
  struct gadpart *a = (struct gadpart  *)first;
  struct gadpart *b = (struct gadpart  *)second;
  if (a->pos[2] > b->pos[2]) return 1;
  else if (a->pos[2] < b->pos[2]) return -1;
  else return 0;
}

int cmp_dist (const void *first, const void *second)
{
  gadpart_dist *a = (gadpart_dist  *)first;
  gadpart_dist *b = (gadpart_dist  *)second;
  if (a->dist > b->dist) return 1;
  else if (a->dist < b->dist) return -1;
  else return 0;
}

int cmp_int (const void *first, const void *second)
{
  int *a= (int *) first; 
  int *b= (int *) second; 
  if (*a>*b) return 1;
  else if (*a<*b) return -1;
  else return 0;
}



unsigned int readgadget(char *filename, struct header *h, fltarr **p, fltarr **v, int **n, float **m)
{
  int blocksize, blocksize2, i,j, b=0, l, dum;
  unsigned int numpart=0, check=0;
  fltarr *pos, *vel;
  float *mass, *mass_dum;
  int *id;
  FILE *fp;
  fp=fopen(filename,"r");
  if (fp==NULL) return 0;
  SKIP;
  fread(h,sizeof(struct header),1,fp);            
  SKIP2;
  if (blocksize!=blocksize2) return 0;
  for (i=0; i<6; i++) {numpart+=h->npart[i]; if ((h->npart[i]!=0) && (h->massarr[i]==0)) b+=h->npart[i];}
  pos=(fltarr *)malloc(sizeof(fltarr)*numpart); 
  vel=(fltarr *)malloc(sizeof(fltarr)*numpart); 
  mass=(float *)malloc(sizeof(float)*numpart); 
  if (b) mass_dum=(float  *)malloc(sizeof(float)*b);
  id=  (int    *)malloc(sizeof(int)*numpart);
  SKIP;                                                
   if (!fread(&pos[0],sizeof(fltarr),numpart,fp)) return 0;
  SKIP2;
  if (blocksize!=blocksize2) return 0;
  SKIP;
   if (!fread(&vel[0],sizeof(fltarr),numpart,fp)) return 0;
  SKIP2;
  if (blocksize!=blocksize2) return 0;
  SKIP;
   if (!fread(&id[0],sizeof(int),numpart,fp)) return 0;
  SKIP2;
  if (blocksize!=blocksize2) return 0;
  if (b)
    {
  SKIP;
   if (!fread(&mass_dum[0],sizeof(float),b,fp)) return 0;
  SKIP2;
  if (blocksize!=blocksize2) return 0;
    }
  fclose(fp);
  *p=pos;
  *v=vel;
  *n=id;
  j=0;
  for (i=0; i< numpart; i++) 
   {
     dum=0;
     for (l=0; l<6; l++)
       {
	 dum+=h->npart[l];
	 if (dum>i) break;
       }
     if (h->massarr[l]!=0) mass[i]=h->massarr[l];
     else mass[i]=mass_dum[j++];
   } 
  *m=mass;
  for (i=0; i<6; i++) {check+=h->npart[i];}
  if (check==numpart) return numpart; else return 0;
}

int convertunits(struct header *head, struct gadpart *part, double convert_mass, double convert_distance)
{
  /*
    convert units (usually to kpc/h, 10^10 Msun/h)
   */

  unsigned int i;
  int j;
  unsigned int numpart = 0;
  for ( i = 0; i < 6; i++) 
    {
      numpart += head->nall[i];
      head->massarr[i] *= convert_mass;
    }
  if (convert_distance != 1.0)
    {
      head->boxsize *= convert_distance;
      for ( i = 0; i < numpart; i++)
	{
	  for ( j = 0; j < 3; j++)
	    part[i].pos[j] *= convert_distance;
	  if (part[i].type == 0)
	    {
	      part[i].sph->hsml *= convert_distance;
	    }
	}
    }
  if (convert_mass != 1.0)
    for ( i = 0; i < numpart; i++)
      {
	part[i].mass *= convert_mass;
      }

  return 0;
}

unsigned int readgadget_part(char *basefilename, struct header *h, struct gadpart **particle)
{
  char filename[256];
  int blocksize, blocksize2, i,j, b=0, l, dum, ngas, nstars, numfiles=1, fnr=0;
  // byte switch_endianess = 0;
  long pi=0;
  unsigned int numpart=0, numpart_all=0, check=0;
  fltarr *pos, *vel;
  float *mass, *mass_dum;
  sphdata *sphtmp;
  int snapformat = 1;
#ifdef POTENTIAL
  float *pot;
#endif
#ifdef LONGIDS
  long *id;
#else
  int *id;
#endif //LONGIDS
  int *type;
  FILE *fp;
  for ( fnr = 0; fnr < numfiles; fnr++ )
    {
      sprintf(filename,"%s", basefilename);
      fp=fopen(filename,"r");
      if (fp==NULL)
	{
	  sprintf(filename, "%s.%d", basefilename, fnr);
	  fp=fopen(filename,"r");
	}
      if (fp==NULL)
	{
	  libgaderr=1;
	  return 0;
	}
      SKIP;
      if (blocksize == 8)
	{
	  snapformat = 2;
	  fseek(fp, 12, SEEK_CUR);
	  SKIP;
	}
      //      if (blocksize == 65536) switch_endianess = 1;
      fread(h,sizeof(struct header),1,fp);            
      numfiles = h->numfiles;
      ngas=h->npart[0];
      nstars=h->npart[4];
      SKIP2;
      if (blocksize!=blocksize2) {libgaderr=256; return 0;}
      numpart = 0;
      b=0;
      for (i=0; i<6; i++) 
	{
	  numpart+=h->npart[i]; 
	  if (fnr==0) numpart_all += h->nall[i];
	  if ((h->npart[i]!=0) && (h->massarr[i]==0)) b+=h->npart[i];
	}
      pos=(fltarr *)malloc(sizeof(fltarr)*numpart); 
      if (fnr==0)
	{
	  *particle= (struct gadpart *) malloc(sizeof (struct gadpart)*numpart_all);	  
	  if (*particle==NULL) {libgaderr=99; return 0;}
	}
      SKIPFORMAT2;                                                
      if (!fread(&pos[0],sizeof(fltarr),numpart,fp)) {libgaderr=10;return 0;}
      SKIP2;
      if (blocksize!=blocksize2) {libgaderr=20;return 0;}
      for (i=0; i<numpart; i++)
	{
	  for (j=0; j<3; j++) (*particle)[pi + i].pos[j] =pos[i][j];
	}
      free(pos);

#ifndef NOVEL
      vel=(fltarr *)malloc(sizeof(fltarr)*numpart); 
      SKIPFORMAT2;
      if (!fread(&vel[0],sizeof(fltarr),numpart,fp)) {libgaderr=11;return 0;}
      SKIP2;
      if (blocksize!=blocksize2) {libgaderr=21;return 0;}
      for (i=0; i<numpart; i++)
	{
	  for (j=0; j<3; j++) (*particle)[pi + i].vel[j] =vel[i][j];
	}
      free(vel);
#else
      SKIPFORMAT2;
      fseek(fp,blocksize,SEEK_CUR);
      SKIP2;
      if (blocksize!=blocksize2) {libgaderr=21;return 0;}
#endif //NOVEL

      mass=(float *)malloc(sizeof(float)*numpart); 
      if (b) mass_dum=(float  *)malloc(sizeof(float)*b);
#ifdef LONGIDS
      id  =  (long    *)malloc(sizeof(long)*numpart);
      SKIPFORMAT2;
      if (!fread(&id[0],sizeof(long),numpart,fp)) {libgaderr=12;return 0;}
      SKIP2;
#else
      id  =  (int    *)malloc(sizeof(int)*numpart);
      SKIPFORMAT2;
      if (!fread(&id[0],sizeof(int),numpart,fp)) {libgaderr=12;return 0;}
      SKIP2;
#endif // LONGIDS
      type=  (int    *)malloc(sizeof(int)*numpart);
      if (type==NULL) {libgaderr=2; return 0;}
      if (blocksize!=blocksize2) {libgaderr=22;return 0;}
      if (b)
	{
	  SKIPFORMAT2;
	  if (!fread(&mass_dum[0],sizeof(float),b,fp)) {libgaderr=13;return 0;}
	  SKIP2;
	  if (blocksize!=blocksize2) {libgaderr=23;return 0;}
	}
      //  part->pos=pos;
      //  part->vel=vel;
      //  part->id=id;
      j=0;
      for (i=0; i< numpart; i++) 
	{
	  dum=0;
	  for (l=0; l<6; l++)
	    {
	      dum+=h->npart[l];
	      if (dum>i) break;
	    }
	  if (h->massarr[l]!=0) mass[i]=h->massarr[l];
	  else mass[i]=mass_dum[j++];
	  type[i]=l;
	} 

      for (i=0; i<numpart; i++)
	{
	  (*particle)[pi + i].id  =id[i];
	  (*particle)[pi + i].mass=mass[i];      
	  (*particle)[pi + i].type=type[i];
#ifndef NOGAS
	  (*particle)[pi + i].sph=NULL;
#if defined(WINDS) || defined(METALS)
	  (*particle)[pi + i].sd=NULL;
#endif //WINDS
#ifdef METALS
	  (*particle)[pi + i].metals=NULL;
#endif
	  (*particle)[pi + i].stellarage=0;
#endif
	}


#ifdef WINDS
      SKIPFORMAT2;
      pot=(float*) malloc (sizeof(float)* numpart);
      if (!fread(&pot[0],sizeof(float),numpart,fp)) {libgaderr=38;return 0;}
      SKIP2;
      if (blocksize!=blocksize2) {libgaderr=48;return 0;}
      for (i=0; i< (numpart); i++)
	{
	  (*particle)[pi + i].pot=pot[i];
	}
      free(pot);


#endif //WINDS

      /*********************************************************************************/
      // Read SPH Properties  
#ifndef NOGAS
      if (!basic)
	{
	  if (ngas)
	    {
	      sphtmp= (struct sphdata *) malloc(sizeof (struct sphdata)* ngas);
	      float *sph[6];
	      for (i=0; i<6; i++)
		{
		  sph[i]=(float  *)calloc(ngas, sizeof(float));
		}
	      for (i=0; i<6; i++)
		{
		  if ( (i==5) && (!(h->flg_sfr)) ) 
		    {
		      continue;
		    } 
		  SKIPFORMAT2;
		  if (!fread(&sph[i][0],sizeof(float),ngas,fp)) {libgaderr=30+i;return 0;}
		  SKIP2;
		  if (blocksize!=blocksize2) {libgaderr=40+i;return 0;}
		}
	      for (i=0; i< ngas; i++)
		{
		  sphtmp[i].u    =sph[0][i];
		  sphtmp[i].rho  =sph[1][i];
		  sphtmp[i].nelec=sph[2][i];
		  sphtmp[i].nh   =sph[3][i];
		  sphtmp[i].hsml =sph[4][i];
		  sphtmp[i].sfr  =sph[5][i];
		  //	      (*particle)[i].sph=&(sphtmp[i]);
		}
	      for (i=0; i<6; i++) free(sph[i]);

#ifdef WINDS

	      stardata *sdtmp = (struct stardata*) malloc (sizeof (struct stardata) * nstars);
	  
	      //Read delaytime
	      float *sphdum = (float  *)malloc(sizeof(float)*ngas);
	      SKIPFORMAT2;
	      if (!fread(&sphdum[0],sizeof(float),ngas,fp)) {libgaderr=36;return 0;}
	      SKIP2;
	      if (blocksize!=blocksize2) {libgaderr=46;return 0;}
	      for (i=0; i< ngas; i++) sphtmp[i].dtime = sphdum[i];
	      free(sphdum);

	      //Read Metals
	      float *metals= (float *) malloc( 4 * sizeof(float) * (ngas+nstars));
	      SKIPFORMAT2;
	      if (!fread(&metals[0],4 * sizeof(float),(ngas + nstars),fp)) {libgaderr=37;return 0;}
	      SKIP2;
	      if (blocksize!=blocksize2) {libgaderr=47;return 0;}
	      for (i=0; i< ngas; i++)
		{
		  for ( j = 0; j < 4; j++)
		    {
		      sphtmp[i].metals[j] = metals[ i*4 + j ];
		    }
		}
	      for (i=ngas; i< (ngas+nstars); i++)
		{
		  for ( j = 0; j < 4; j++)
		    {
		      sdtmp[i-ngas].metals[j] = metals[ i*4 + j ];
		    }
		}
	      free(metals);

	      //Read tmax
	      float *gsdum = (float  *)malloc(sizeof(float)*(ngas+nstars));
	      SKIPFORMAT2;
	      if (!fread(&gsdum[0],sizeof(float),(ngas+nstars),fp)) {libgaderr=38;return 0;}
	      SKIP2;
	      if (blocksize!=blocksize2) {libgaderr=48;return 0;}				  
	      for (i=0; i< ngas; i++)
		{
		  sphtmp[i].tmax = gsdum[ i ];	     
		}
	      for (i=ngas; i< (ngas+nstars); i++)
		{
		  sdtmp[i-ngas].tmax = gsdum[ i ];
		}
	  
	      //Read n_spawn
	      SKIPFORMAT2;
	      if (!fread(&gsdum[0],sizeof(float),(ngas+nstars),fp)) {libgaderr=39;return 0;}
	      SKIP2;
	      if (blocksize!=blocksize2) {libgaderr=49;return 0;}				  
	      for (i=0; i< ngas; i++)
		{
		  sphtmp[i].n_spawn = gsdum[ i ];	     
		}
	      for (i=ngas; i< (ngas+nstars); i++)
		{
		  sdtmp[i-ngas].n_spawn = gsdum[ i ];
		}
	      free(gsdum);

	      int stars_start = h->npart[0] + h->npart[1] + h->npart[2] + h->npart[3]; 
	      for (i=stars_start; i< (stars_start+nstars); i++)
		{
		  (*particle)[pi + i].sd=&(sdtmp[i-stars_start]);
		}

#endif //WINDS

	      for (i=0; i< ngas; i++)
		{
		  (*particle)[pi + i].sph=&(sphtmp[i]);
		}
	    }
	  /*********************************************************************************/
	  if ((nstars) || (h->npart[5]))
	    {
	      float *sa;
	      sa=(float*) malloc(sizeof(float)*(nstars + h->npart[5]));
	      SKIPFORMAT2;
	      if (!fread(&sa[0],sizeof(float),(nstars + h->npart[5]),fp)) {libgaderr=37;return 0;}
	      SKIP2;
	      if (blocksize!=blocksize2) {libgaderr=47;return 0;}
	      int start=h->npart[0] + h->npart[1] + h->npart[2] + h->npart[3]; 
	      for (i=0; i< (nstars + h->npart[5]); i++)
		{
		  (*particle)[pi + start + i].stellarage=sa[i];
		}
	      free(sa);
	    }
	  /*********************************************************************************/

#ifdef METALS
	  
	  if (nstars)
	    {
	      SKIPFORMAT2;
	      if ( blocksize == sizeof(int)*nstars )
		{
		  stardata *sdtmp = (struct stardata*) malloc (sizeof (struct stardata) * nstars);
		  int *let;
		  let=(int*) malloc(sizeof(int)*nstars);
	      
		  if (!fread(&let[0],sizeof(int),nstars,fp)) {libgaderr=38;return 0;}
		  SKIP2;
		  if (blocksize!=blocksize2) {libgaderr=48;return 0;}
		  int start=h->npart[0] + h->npart[1] + h->npart[2] + h->npart[3]; 
		  for (i=0; i< (nstars); i++)
		    {
		      sdtmp[i].let = let[i];

		    }
		  free(let);

		  float *initialmass;
		  initialmass=(float*) malloc(sizeof(float)*nstars);
		  SKIPFORMAT2;
		  if (!fread(&initialmass[0],sizeof(float),nstars,fp)) {libgaderr=39;return 0;}
		  SKIP2;
		  if (blocksize!=blocksize2) {libgaderr=49;return 0;}

		  for (i=0; i< (nstars); i++)
		    {
		      sdtmp[i].initialmass=initialmass[i];
		      (*particle)[pi + start + i].sd = &(sdtmp[i]);
		    }
		  free(initialmass);
		}
	      else
		{
		  fseek(fp, -sizeof(int), SEEK_CUR);
		}
	    }

	  if (((nstars) || (ngas)) && (h->flg_mtl))
	    {
	      float *metals;
	      metals = (float*) malloc(sizeof(float)*12*(nstars+ngas));
	      SKIPFORMAT2;
	      if (!fread(&metals[0],sizeof(float),12*(nstars+ngas),fp)) {libgaderr=110;return 0;}
	      SKIP2;
	      if (blocksize!=blocksize2) {libgaderr=120;return 0;}

	      for (i=0; i< (ngas); i++)
		{
		  (*particle)[pi + i].metals=&(metals[i*12]);
		}

	      int start=h->npart[0] + h->npart[1] + h->npart[2] + h->npart[3]; 
	      for (i=0; i< (nstars); i++)
		{
		  (*particle)[pi + start + i].metals = &(metals[(ngas+i)*12]);
		}
	      
	    }

#endif //METALS

	  /* BH data should be read in here, temporarily skipped */	
	  if (h->npart[5])
	    {
	      SKIPFORMAT2;
	      while ( blocksize == sizeof(float)*h->npart[5] )
		{
		  SKIPFORMAT2;
		}
	      fseek(fp, -sizeof(int), SEEK_CUR);
	    }

#ifdef POTENTIAL
#ifndef WINDS

	  SKIPFORMAT2;
	  if (!feof(fp))
	    {	      	      
	      pot=(float*) malloc (sizeof(float)* numpart);
	      if (!fread(&pot[0],sizeof(float),numpart,fp)) {libgaderr=68;return 0;}
	      SKIP2;
	      if (blocksize!=blocksize2) 
		{
		  libgaderr=69;
		  for (i=0; i< (numpart); i++)
		    {
		      (*particle)[pi + i].pot=0;
		    }
		  //		  return 0;
		}
	      else
		for (i=0; i< (numpart); i++)
		  {
		    (*particle)[pi + i].pot=pot[i];
		  }
	      free(pot);
	    } 
	  else 
	    {
	      for (i=0; i< (numpart); i++)
		{
		  (*particle)[pi + i].pot=0;
		}
	    }
#endif //WINDS
#endif //POTENTIAL

#ifdef METALS
	  if (ngas)
	    {
	      float *temp;
	      temp = (float*) malloc(sizeof(float)*ngas);
	      SKIPFORMAT2;
	      if (blocksize == sizeof(float) * ngas)
		{
		  if (!fread(&temp[0],sizeof(float),ngas,fp)) {libgaderr=111;return 0;}
		  SKIP2;
		  if (blocksize!=blocksize2) {libgaderr=121;return 0;}
		  
		  for (i=0; i< (ngas); i++)
		    {
		      (*particle)[pi + i].sph->temp=temp[i];
		    }
		}
	      else
		{
		  for (i=0; i< (ngas); i++)
		    {
		      (*particle)[pi + i].sph->temp=temperature((*particle)[pi + i]);
		    }
		  fseek(fp, -sizeof(int), SEEK_CUR);
		}
	      free(temp);	      
	    }
	  
#endif // METALS
	}
      fclose(fp);
#endif

      /*********************************************************************************/
      free(mass);
      free(id);
      free(type);
      if (b) free(mass_dum);
      pi += numpart;
    }

  for (i=0; i<6; i++) 
    {
      check+=h->nall[i];
      h->npart[i] = h->nall[i];
    }
  if (check==numpart_all) return numpart_all; else {libgaderr=30;return 0;}
}

#ifdef LONGIDS
unsigned int writegadget(char *filename, struct header h, fltarr *p, fltarr *v, long *n, float *m)
#else 
unsigned int writegadget(char *filename, struct header h, fltarr *p, fltarr *v, int *n, float *m)
#endif
{
  int blocksize, blocksize2, i,j, b=0, l, dum;
  unsigned int numpart=0, check=0, ngas=0;
  FILE *fp;
  for (i=0; i<6; i++)
    {
      numpart+=h.npart[i];
      h.massarr[i]=0;
    }
  ngas=h.npart[0];
  fp=fopen(filename,"w");
  blocksize=256;
  BLOCK
    fwrite(&h, sizeof(struct header),1,fp);
  BLOCK
    
  blocksize=numpart*12;
  BLOCK
    fwrite(&p[0], sizeof(fltarr), numpart, fp);
  BLOCK

  blocksize=numpart*12;
  BLOCK
    fwrite(&v[0], sizeof(fltarr), numpart, fp);
  BLOCK

#ifdef LONGIDS
  blocksize=numpart*sizeof(long);
  BLOCK
    fwrite(&n[0], sizeof(long), numpart, fp);
  BLOCK
#else  
  blocksize=numpart*4;
  BLOCK
    fwrite(&n[0], sizeof(int), numpart, fp);
  BLOCK
#endif

  blocksize=numpart*4;
  BLOCK
    fwrite(&m[0], sizeof(float), numpart, fp);
  BLOCK
    if (ngas)
      {
	dum=0;
	blocksize=ngas*4;
	BLOCK
	  for (i=0; i<ngas; i++) fwrite(&dum, sizeof(int), 1, fp);
	BLOCK
	BLOCK
	  for (i=0; i<ngas; i++) fwrite(&dum, sizeof(int), 1, fp);
	BLOCK
      }
  fclose(fp);
  return numpart;
}

unsigned int writegadget_part(char *filename, struct header h, struct gadpart *part)
{
#ifdef DEBUG
  printf("LIBGAD: entering writegadget_part\n");
#endif
  int blocksize, blocksize2, i,j, b=0, l, dum;
  unsigned int numpart=0, check=0, ngas=0;
  fltarr *p, *v;
#ifdef LONGIDS
  long *n;
#else
  int *n;
#endif
  float *m;
  FILE *fp;
  for (i=0; i<6; i++)
    {
      numpart+=h.npart[i];
      h.massarr[i]=0;
    }
  ngas=h.npart[0];
  int nstars= h.npart[4];
  qsort(part, numpart, sizeof(struct gadpart), cmp_type);
  p= (fltarr *) malloc (sizeof(fltarr)*numpart);
  v= (fltarr *) malloc (sizeof(fltarr)*numpart);
#ifdef LONGIDS
  n= (long *)    malloc (sizeof(long)*numpart);
#else
  n= (int *)    malloc (sizeof(int)*numpart);
#endif
  m= (float *)  malloc (sizeof(float)*numpart);
  if ( (p==NULL) || (v==NULL) || (n==NULL) || (m==NULL) )
    {
      libgaderr=-1;
      return 0;
    }
#ifdef POTENTIAL
 float *pot= (float *)  malloc (sizeof(float)*numpart);
#endif
  for (i=0; i<numpart; i++) {
    for (j=0; j<3; j++) p[i][j]=part[i].pos[j];
#ifndef NOVEL
    for (j=0; j<3; j++) v[i][j]=part[i].vel[j];
#else
    for (j=0; j<3; j++) v[i][j]=0;
#endif //NOVEL
    n[i]=part[i].id;
    m[i]=part[i].mass;
#ifdef POTENTIAL
    pot[i]=part[i].pot;
#endif
  }
  fp=fopen(filename,"w");
  if (fp==NULL) 
    {
      libgaderr=1;
      return 0;
    }
  blocksize=256;
  BLOCK
    fwrite(&h, sizeof(struct header),1,fp);
  BLOCK
  
  blocksize=numpart*12;
  BLOCK
    fwrite(&p[0], sizeof(fltarr), numpart, fp);
  BLOCK

  blocksize=numpart*12;
  BLOCK
    fwrite(&v[0], sizeof(fltarr), numpart, fp);
  BLOCK
  
#ifdef LONGIDS
  blocksize=numpart*sizeof(long);
  BLOCK
    fwrite(&n[0], sizeof(long), numpart, fp);
  BLOCK
#else
  blocksize=numpart*4;
  BLOCK
    fwrite(&n[0], sizeof(int), numpart, fp);
  BLOCK
#endif // LONGIDS

  blocksize=numpart*4;
  BLOCK
    fwrite(&m[0], sizeof(float), numpart, fp);
  BLOCK

#ifdef WINDS
  blocksize=numpart*4;
  BLOCK
    fwrite(pot, sizeof(float), numpart, fp);
  BLOCK
    free(pot);
#endif

#ifndef NOGAS
    if ((ngas) && (part[0].sph!=NULL))
      {
	dum=0;
	blocksize=ngas*4;
	float *sph[6];
	for (i=0; i<6; i++)
	  {
	    sph[i]=(float  *)malloc(sizeof(float)*ngas);
	  }
	for (i=0; i<ngas; i++)
	  {
	    sph[0][i]=part[i].sph->u;
	    sph[1][i]=part[i].sph->rho;
	    sph[2][i]=part[i].sph->nelec;
	    sph[3][i]=part[i].sph->nh;
	    sph[4][i]=part[i].sph->hsml;
	    sph[5][i]=part[i].sph->sfr;
	  }
	for (i=0; i<6; i++)
	  {
	    if ( (i==5) && (!(h.flg_sfr)) ) 
	      {
		continue;
	      } 
	    BLOCK
	      fwrite(&sph[i][0], sizeof(float), ngas, fp);
	    BLOCK
	  }
	for (i=0; i<6; i++) free(sph[i]);

#ifdef WINDS
	//dtime
	float *sphdum = (float  *)malloc(sizeof(float)*ngas);
	for (i=0; i<ngas; i++)
	  {
	    sphdum[i] = part[i].sph->dtime;
	  }
	blocksize=ngas*4;
	BLOCK
	  fwrite(&sphdum[0], sizeof(float), ngas, fp);
	BLOCK
	free(sphdum);


	//metals
	float *metals= (float *) malloc( 4 * sizeof(float) * (ngas+nstars));

	for (i=0; i< ngas; i++)
	  {
	    for ( j = 0; j < 4; j++)
	      {
		metals[ i*4 + j ] = part[i].sph->metals[j];
	      }
	  }
	int stars_start = h.npart[0] + h.npart[1] + h.npart[2] + h.npart[3];
	for (i=ngas; i< (ngas+nstars); i++)
	  {
	    for ( j = 0; j < 4; j++)
	      {
		if (part[stars_start + i - ngas ].sd == NULL) 
		  {
		    printf("NULL\n");fflush(stdout);
		    exit(0);
		  }
		metals[ i*4 + j ] = part[stars_start + i - ngas ].sd->metals[j];
	      }
	  }

	blocksize=(ngas+nstars) * sizeof(float) * 4;
	BLOCK
	  fwrite(&metals[0], 4 * sizeof(float), (ngas+nstars), fp);
	BLOCK
	  free (metals);

	//tmax
	float *gsdum = (float  *)malloc(sizeof(float)*(ngas+nstars));

	for (i=0; i< ngas; i++)
	  {
	    gsdum[ i ] = part[i].sph->tmax;
	  }

	for (i=ngas; i< (ngas+nstars); i++)
	  {
	    gsdum[ i ] = part[stars_start + i - ngas].sd->tmax;
	  }

	blocksize= (ngas + nstars) * sizeof(float);
	BLOCK
	  fwrite(&gsdum[0], sizeof(float), (ngas+nstars), fp);
	BLOCK

	  //n_spawn
	for (i=0; i< ngas; i++)
	  {
	    gsdum[ i ] = part[i].sph->n_spawn;
	  }
	for (i=ngas; i< (ngas+nstars); i++)
	  {
	    gsdum[ i ] = part[stars_start + i -ngas].sd->n_spawn;
	  }
	blocksize= (ngas + nstars) * sizeof(float);
	BLOCK
	  fwrite(&gsdum[0], sizeof(float), (ngas+nstars), fp);
	BLOCK 
	  free(gsdum);
	
#endif


      }
    if ((nstars) || h.npart[5] )
    {
      float *sa=(float  *)malloc(sizeof(float)*(nstars+h.npart[5]));
      int start= h.npart[0] + h.npart[1] + h.npart[2] + h.npart[3];
      for (i=0; i<nstars+h.npart[5]; i++)
	  {
	    sa[i]=part[i + start].stellarage;
	  }
      blocksize=(nstars+h.npart[5])*4;
      BLOCK
	fwrite(sa, sizeof(float), nstars+h.npart[5], fp);
      BLOCK
	free(sa);
    }
#endif
#ifdef METALS
  if (nstars)
    {
      int start= h.npart[0] + h.npart[1] + h.npart[2] + h.npart[3];
      if (part[start].sd !=NULL)
	{
	  int *let=(int  *)malloc(sizeof(int)*nstars);	  
	  for (i=0; i<nstars; i++)
	    {
	      let[i]=part[i + start].sd->let;
	    }
	  blocksize=nstars*4;
	  BLOCK
	    fwrite(let, sizeof(int), nstars, fp);
	  BLOCK
	    free(let);
	  
	  float *initialmass=(float  *)malloc(sizeof(float)*nstars);
	  for (i=0; i<nstars; i++)
	    {
	      initialmass[i]=part[i + start].sd->initialmass;
	    }
	  blocksize=nstars*4;
	  BLOCK
	    fwrite(initialmass, sizeof(float), nstars, fp);
	  BLOCK
	    free(initialmass);
	}
    }
  if (((nstars) || (ngas)) && (h.flg_mtl))
    {
      float *metals=(float  *)malloc(sizeof(float)*12*(nstars+ngas));

      for (i=0; i<ngas; i++)
	  {
	    for (j = 0; j < 12; j++ )
	      metals[i*12 + j]=part[i].metals[j];
	  }
      int start= h.npart[0] + h.npart[1] + h.npart[2] + h.npart[3];
      for (i=0; i<nstars; i++)
	  {
	    for (j = 0; j < 12; j++ )
	      metals[(i+ngas)*12 + j]=part[i + start].metals[j];
	  }
      blocksize=(nstars+ngas)*4*12;
      BLOCK
	fwrite(metals, sizeof(float), 12*(ngas+nstars), fp);
      BLOCK
	free(metals);
    }
#endif //METALS

  /* BH data should be written here, temporarily skipped */	
  if (h.npart[5])
    {
      for ( i = 0; i < 4; i++ )
	{
	  blocksize = sizeof(float) * h.npart[5];
	  float* fdum = (float*) calloc(h.npart[5], sizeof(float));
	  BLOCK
	    fwrite(fdum, sizeof(float), h.npart[5], fp);
	  BLOCK
	    free(fdum);
	}
    }

#ifdef POTENTIAL
#ifndef WINDS
  if (!basic)
    {
      blocksize=numpart*4;
      BLOCK
	fwrite(pot, sizeof(float), numpart, fp);
      BLOCK
	free(pot);
    }
#endif
#endif

#ifdef METALS
  if (ngas)
    {
      float *temp=(float  *)malloc(sizeof(float)*ngas);
      for (i=0; i<ngas; i++)
	  {
	    temp[i]=part[i].sph->temp;
	  }
      blocksize=ngas*4;
      BLOCK
	fwrite(temp, sizeof(float), ngas, fp);
      BLOCK
	free(temp);
    }
#endif
  fclose(fp);
  free(p);
  free(v);
  free(n);
  free(m);
  return numpart;
}

unsigned int readgadget_novel(char *filename, struct header *h, fltarr **p, int **n, float **m)
{
  int blocksize, blocksize2, i, b=0;
  unsigned int numpart=0, check=0;
  fltarr *pos;
  float *mass;
  int *id;
  FILE *fp;
  fp=fopen(filename,"r");
  if (fp==NULL) return 0;
  SKIP;
  fread(h,sizeof(struct header),1,fp);            
  SKIP2;
  if (blocksize!=blocksize2) return 0;
  for (i=0; i<6; i++) {numpart+=h->npart[i]; if ((h->npart[i]!=0) && (h->massarr[i]==0)) b+=h->npart[i];}
  pos=(fltarr *)malloc(sizeof(fltarr)*numpart); 
  if (b) mass=(float  *)malloc(sizeof(float)*b);
  id=  (int    *)malloc(sizeof(int)*numpart);
  SKIP;                                                
   if (!fread(&pos[0],sizeof(fltarr),numpart,fp)) return 0;
  SKIP2;
  if (blocksize!=blocksize2) return 0;
  SKIP;
  fseek(fp,blocksize,SEEK_CUR);
  SKIP2;
  if (blocksize!=blocksize2) return 0;
  SKIP;
   if (!fread(&id[0],sizeof(int),numpart,fp)) return 0;
  SKIP2;
  if (blocksize!=blocksize2) return 0;
  if (b)
    {
  SKIP;
   if (!fread(&mass[0],sizeof(float),b,fp)) return 0;
  SKIP2;
  if (blocksize!=blocksize2) return 0;
    }
  fclose(fp);
  *p=pos;
  *n=id;
  if (b) *m=mass;
  for (i=0; i<6; i++) {check+=h->npart[i];}
  if (check==numpart) return numpart; else return 0;
}

unsigned int readgadget_sph(char *filename, struct header *h, fltarr **p, fltarr **v, int **n, float **m, float **u, float **rho, float **d1, float **d2, float **d3, float **d4 , float **sa)
{
  int blocksize, blocksize2, i, noma=0;
  unsigned int numpart=0, check=0;
  fltarr *pos, *vel;
  float *mass, *stellarage, *sph[6];
  int *id;
  FILE *fp;
  fp=fopen(filename,"r");
  if (fp==NULL) return 0;
  SKIP;
  fread(h,sizeof(struct header),1,fp);            
  SKIP2;
  if (blocksize!=blocksize2) return 0;
  for (i=0; i<6; i++) {numpart+=h->npart[i]; if ((h->npart[i]!=0) && (h->massarr[i]==0)) noma+=h->npart[i];}
  pos=(fltarr *)malloc(sizeof(fltarr)*numpart); 
  vel=(fltarr *)malloc(sizeof(fltarr)*numpart); 
  if (noma) mass=(float  *)malloc(sizeof(float)*noma);
  id=  (int    *)malloc(sizeof(int)*numpart);
  for (i=0; i<6; i++)
    {
      sph[i]=(float  *)malloc(sizeof(float)*h->npart[0]);
    }

  stellarage =(float  *)malloc(sizeof(float)*h->npart[4]);

  SKIP;                                                
   if (!fread(&pos[0],sizeof(fltarr),numpart,fp)) return 0;
  SKIP2;
  if (blocksize!=blocksize2) return 0;
  SKIP;
   if (!fread(&vel[0],sizeof(fltarr),numpart,fp)) return 0;
  SKIP2;
  if (blocksize!=blocksize2) return 0;
  SKIP;
   if (!fread(&id[0],sizeof(int),numpart,fp)) return 0;
  SKIP2;
  if (blocksize!=blocksize2) return 0;
  if (noma)
    {
  SKIP;
   if (!fread(&mass[0],sizeof(float),noma,fp)) return 0;
  SKIP2;
  if (blocksize!=blocksize2) return 0;
    }
  for (i=0; i<6; i++)
    {
      SKIP;
      if (!fread(&sph[i][0],sizeof(float),h->npart[0],fp)) return 0;
      SKIP2;
    }
  if (h->npart[4]!=0)
    {
      SKIP;
      if (!fread(&stellarage[0],sizeof(float),h->npart[4],fp)) return 0;
      SKIP2;
    }
  fclose(fp);
  *p=pos;
  *v=vel;
  *n=id;
  if (noma) *m=mass;
  *u  =sph[0];
  *rho=sph[1];
  *d1 =sph[2];
  *d2 =sph[3];
  *d3 =sph[4];
  *d4 =sph[5];
  *sa =stellarage;
  for (i=0; i<6; i++) {check+=h->npart[i];}
  if (check==numpart) return numpart; else return 0;
}

double distance(fltarr a, fltarr b)
{
  int i;
  fltarr d;
  double dist=0;
  if (!BOXSIZE) BOXSIZE=72000;
  for (i=0; i<3; i++)
    {
      d[i]=ABS(a[i]-b[i]);
      d[i]=MIN(d[i], BOXSIZE-d[i]);
      dist+=SQR(d[i]);
    }
  return sqrt(dist);
}

double distance_nopb(fltarr a, fltarr b)
{
  int i;
  fltarr d;
  double dist=0;
  for (i=0; i<3; i++)
    {
      d[i]=ABS(a[i]-b[i]);
      dist+=SQR(d[i]);
    }
  return sqrt(dist);
}

double distbox(fltarr a, fltarr min, fltarr max)
{
  double dist=0, dmin, dmax;
  int i;
  if (!BOXSIZE) BOXSIZE=72000;

  for (i=0; i<3; i++)
    {
      if (min[i]<=max[i])
	{
	  if ((a[i]<min[i]) || (a[i]>max[i]))
	    {
	      dmin=ABS(a[i]-min[i]);
	      dmin=MIN(dmin, BOXSIZE-dmin);
	      dmax=ABS(a[i]-max[i]);
	      dmax=MIN(dmax, BOXSIZE-dmax);
	      dist+=SQR(MIN(dmin, dmax));
	    }
	}
      else
	{
	  if ((a[i]<min[i]) && (a[i]>max[i]))
	    {
	      dmin=ABS(a[i]-min[i]);
	      dmin=MIN(dmin, BOXSIZE-dmin);
	      dmax=ABS(a[i]-max[i]);
	      dmax=MIN(dmax, BOXSIZE-dmax);
	      dist+=SQR(MIN(dmin, dmax));
	    }
	}
    }
  dist=sqrt(dist);
  return dist;
}

double distbox_nopb(fltarr a, fltarr min, fltarr max)
{
  double dist=0, dmin, dmax;
  int i;

  for (i=0; i<3; i++)
    {
	  if ((a[i]<min[i]) || (a[i]>max[i]))
	    {
	      dmin=ABS(a[i]-min[i]);
	      dmax=ABS(a[i]-max[i]);
	      dist+=SQR(MIN(dmin, dmax));
	    }
    }
  dist=sqrt(dist);
  return dist;
}

void dummyfunction()
{
  extern float BOXSIZE;
  if (!BOXSIZE) BOXSIZE=72000;
  printf("%f\n", BOXSIZE);
}


void cpygadpart(gadpart * to, gadpart * from)
{
  memcpy(to, from, sizeof(gadpart));
//  int i;
//  for (i=0; i<3; i++)
//    {
//      to->pos[i]= from->pos[i];
//      to->vel[i]= from->vel[i];
//    }
//  to->id  = from->id;
//  to->mass= from->mass;
//  to->type= from->type;

}

int gadsearch(gadpart_dist *data, double toFind, int start, int end)
 {
    int mid=0;
    int location=-1;
    while ((start<=end) && (location==-1))
      {
	//	printf("TEST\n");fflush(stdout);
	if (data[mid].dist!=data[mid].dist) return -1;
	mid = (int) (start + (end - start)/2); 
	if (data[mid].dist== toFind) 
	  {
	    location=mid;
	    break;
	  }
	else if (data[mid].dist < toFind) start=mid+1;
	else if (data[mid].dist > toFind) end  =mid-1;
      }
    if (start>end) return start;
    return location;
 }

static double fit_fct(double t, double *p)
{
  return log10((p[0]) / ((t/p[1]) * SQR(1+t/p[1])));
}

double log_fit (double t, double *p)
{
  return ( p[0] + p[1] * log10(t) );
}

void printout ( int n_par, double* par, int m_dat, double* fvec, 
                       void *data, int iflag, int iter, int nfev )
         {
           // dummy function to catch fitting output

         }

double nfwfit(double *par, gadpart_dist *part, int cnt, double rv, double soft, double *rcs)
{
  const int nbins=50;
  //  return (double)2;
  double dist=0;
  int i=0,j,m;
  double err[nbins], p[nbins], rlog[nbins];
  //  qsort(part, cnt, sizeof(gadpart_dist), cmp_dist);
  for (j=0; j<nbins; j++)
   {
     p[j]=0;
     err[j]=0;
     rlog[j]=0;
   }
  double minr=(rv/pow(10,2.5));
  double lmr=log10(minr);
  double lrv=log10(rv);
  double range=(lrv-lmr);
  double dr=range/nbins;

  while ((dist<rv) && (i<cnt))
   { 
     dist=part[i].dist;
     //     int idum=0;
     m= part[i].part.type;
     if ((m>0) && (m<4))
       {
	 double d=log10(dist);
	 //	 j=floor(d/(log10(rv)/nbins));
	 j=floor((d-lmr)/dr);
	 if (j>=0)
	   {
      if (dist<soft)
         err[j] += 0.01;
      else
	       err[j]++;
	    p[j]+=part[i].part.mass;
	   }
       }
     i++;
   }

  int nfun=0;
  int k=0;
  for (j=0; j<nbins; j++)
   {
     if (err[j]!=0)
       {
	 err[nfun]=1/(sqrt(err[j]));
	 if (k==0) {
	   rlog[nfun]=(pow(10, ((j+1)*dr + lmr )))/2.;
	   p[nfun]=p[j]/((4.0/3.0)* M_PI * (pow(10, (((j+1)*dr+lmr)*3))));
	 }
	 else {
	   double ddum=p[j];
	   rlog[nfun]=(pow(10, ((j+1)*dr + lmr )) + pow(10, (k*dr) + lmr))/2.;
	   p[nfun]=p[j]/((4.0/3.0)* M_PI * (pow(10, (((j+1)*dr + lmr)*3)) - pow(10, ((k*dr + lmr)*3))));
	 }
	 k=j+1;
	 nfun++;
       } 
   }
  for (j=0; j<nfun; j++)
   {
     p[j]=log10(p[j]);
   }
  // auxiliary settings for fitting:

    lm_control_type control;
    lm_data_type data;
    lm_initialize_control(&control);

    data.user_func = fit_fct;
    data.user_t = rlog;
    data.user_y = p;
//---------------------------------------------------
    if (par[0]==0)
        par[0]=0.005;
    if (par[1]==0)
        par[1]=20;
    lm_minimize(nfun, 2, par, err, lm_evaluate_default, printout, &data, &control);
    *rcs=0;
    for (j=0; j<nfun; j++)
      {
	double dens=(p[j]);
	double fit=(fit_fct(rlog[j], par));
	//	(*rcs)+= SQR(pow(10,fit_fct(rlog[j], par))-p[j])/ABS(p[j]);
	(*rcs)+= SQR(fit-dens)/ABS(dens);
      }
    (*rcs)=(*rcs)/(nfun-2);
    return (rv/par[1]);
}

double densproffit(double *par, gadpart_dist *part, int cnt, double re, double soft, double *rcs, int type)
{
  const int nbins=32;
  //  return (double)2;
  double dist=0;
  int i=0,j,m;
  double err[nbins], p[nbins], rlog[nbins];
  //  qsort(part, cnt, sizeof(gadpart_dist), cmp_dist);
  for (j=0; j<nbins; j++)
   {
     p[j]=0;
     err[j]=0;
     rlog[j]=0;
   }
  double minr=MAX(soft, 0.1 * re);
  double maxr=2* re;
  double lminr=log10(minr);
  double lmaxr=log10(maxr);
  double range=(lmaxr-lminr);
  double dr=range/nbins;
  //  printf("%6g %6g %6g %6g %6g\n", minr, maxr, lminr, lmaxr, dr);
  while ((dist<maxr) && (i<cnt))
   { 
     dist=part[i].dist;
     //     int idum=0;
     m= part[i].part.type;
     if ( ( (1<<m) & type) && (dist>minr) )
       {
	 double d=log10(dist);
	 //	 j=floor(d/(log10(rv)/nbins));
	 j=floor((d-lminr)/dr);
	 if (j>=0)
	   {
	     err[j]++;
	     p[j]+=part[i].part.mass;
	   }
       }
     i++;
   }

  int nfun=0;
  int k=0;
  for (j=0; j<nbins; j++)
   {
     if (err[j]!=0)
       {
	 err[nfun]=1/(sqrt(err[j]));
	 if (k==0) {
	   rlog[nfun]=(pow(10, ((j+1)*dr + lminr )) + minr)/2.;
	   p[nfun]=p[j]/((4.0/3.0)* M_PI * (pow(10, (((j+1)*dr+lminr)*3)) - pow(10, (lminr*3)) ) );
	 }
	 else {
	   double ddum=p[j];
	   rlog[nfun]=(pow(10, ((j+1)*dr + lminr )) + pow(10, (k*dr) + lminr))/2.;
	   p[nfun]=p[j]/((4.0/3.0)* M_PI * (pow(10, (((j+1)*dr + lminr)*3)) - pow(10, ((k*dr + lminr)*3))));
	 }
	 k=j+1;
	 nfun++;
       } 
   }
  for (j=0; j<nfun; j++)
   {
     p[j]=log10(p[j]);
   }
  // auxiliary settings for fitting:

    lm_control_type control;
    lm_data_type data;
    lm_initialize_control(&control);

    data.user_func = log_fit;
    data.user_t = rlog;
    data.user_y = p;
//---------------------------------------------------
    par[0]=1.;
    par[1]=-1.5;
    lm_minimize(nfun, 2, par, err, lm_evaluate_default, printout, &data, &control);
    *rcs=0;
    for (j=0; j<nfun; j++)
      {
	double dens=pow(10,p[j]);
	//	double fit=pow(10,fit_fct(rlog[j], par));
	double fit=pow(10,log_fit(rlog[j], par));
	//	(*rcs)+= SQR(pow(10,fit_fct(rlog[j], par))-p[j])/ABS(p[j]);
	(*rcs)+= SQR(fit-dens)/ABS(dens);
	//	printf("%6g %6g %6g\n", rlog[j] / 0.72, log10(dens * 1.e10 * SQR(0.72) ), log10(fit * 1.e10 * SQR(0.72)) );
      }
    (*rcs)=(*rcs)/(nfun-2);
    return (par[1]);
}


void calcdist(gadpart_dist *gd, int cnt, float *center)
{
  int j,i;
  for (i= 0; i< cnt; i++)
    {
      gd[i].dist=distance(gd[i].part.pos, center);
    }
}

void simplecenter(gadpart *part, int cnt, double* cm, int use)
{
  int i,j;
  double totmass=0;
  for ( j = 0; j < 3; j++ )
    cm[j] = 0;
 
  for ( i = 0; i < cnt; i++ )
    {
      if (! (( 1 << part[i].type) & use)) continue;
      totmass += part[i].mass;
      for ( j = 0; j < 3; j++ )
	cm[j] += part[i].pos[j] * part[i].mass;
    }
  for ( j = 0; j < 3; j++ )
    cm[j] /= totmass; 
}

void pcenter(gadpart_dist *part, int cnt, double maxdist,  float* cm, int use)
{
  double tmpcm[3]={0,0,0};
  double dist, ddum, masstot=0;
  gadpart *tmppart= (gadpart*) malloc(cnt * sizeof(gadpart));
  int j,i;
  //  for (j=0; j<3; j++) cm[j]=0;
  if ( ( cm[0] == 0 ) &&  ( cm[1] == 0 ) && ( cm[2] == 0 ) )
    {
      for (i=0; i<cnt; i++)
	{
	  tmppart[i]=part[i].part;
	  if (! (( 1 << part[i].part.type) & use)) continue;
	  masstot+=part[i].part.mass;
	  for (j=0; j<3; j++)
	    cm[j]+=part[i].part.pos[j]*part[i].part.mass;      
	}
      for (j=0; j<3; j++)
	cm[j]/=masstot;
    }
  else
    {
      for (i=0; i<cnt; i++)
	{
	  tmppart[i]=part[i].part;
	}
    }
  int count=cnt;
  int idum;
  do
    {
      for (j=0; j<3; j++) tmpcm[j]=0;
      idum=0;
      masstot=0;
      for (i=0; i<count; i++)
	{
	  if (! (( 1 << tmppart[i].type) & use)) continue;
	  ddum= distance(cm, tmppart[i].pos);
	  if (ddum<maxdist)
	    {
	      for (j=0; j<3; j++) tmpcm[j]+=tmppart[i].pos[j] * tmppart[i].mass;
	      masstot+=tmppart[i].mass;
	      tmppart[idum++]=tmppart[i];
	    }
	}
      if (masstot)
	for (j=0; j<3; j++) 
	  {
	    tmpcm[j]/=masstot;
	    cm[j]=tmpcm[j];
	  }
      maxdist*=0.92;
      count=idum;
    } while (count> 15);
  calcdist(part, cnt, cm);
  free(tmppart);
  return;
}

#ifndef NOVEL
void findcenter(gadpart *part, int cnt, double maxdist, int use)
{
  double tmpcm[3];
  double dist, ddum, masstot=0;
  int *idlist= (int*) calloc(cnt, sizeof(int));
  int i,j,k;
  double cm[3];
  double cvel[3];
  if ( maxdist <= 0 ) maxdist = 20000.0;
  if (!BOXSIZE) BOXSIZE=72000.;
  for (j=0; j<3; j++) {cm[j]=0; cvel[j]=0;tmpcm[j]=0;}
  int count=0;
  for (i=0; i<cnt; i++)
    {
      if (! (( 1 << part[i].type) & use)) continue;
      masstot+=part[i].mass;
      idlist[count++]=i;
      for (j=0; j<3; j++)
	cm[j]+=part[i].pos[j]*part[i].mass;      
    }
  for (j=0; j<3; j++)
    cm[j]/=masstot;
  
  int idum;
  do
    {
      double extend = 0;
      for (j=0; j<3; j++) tmpcm[j]=0;
      idum=0;
      masstot=0;
      for (k=0; k<count; k++)
	{
	  i = idlist[k];
	  if (! (( 1 << part[i].type) & use)) continue;
	  ddum= sqrt( SQR(cm[0]-part[i].pos[0]) + SQR(cm[1]-part[i].pos[1]) + SQR(cm[2]-part[i].pos[2]));
	  if (ddum<maxdist)
	    {
	      for (j=0; j<3; j++) tmpcm[j]+=part[i].pos[j] * part[i].mass;
	      masstot+=part[i].mass;
	      idlist[idum++]=i;
	    }
	  if (ddum > extend) extend = ddum;
	}
      for (j=0; j<3; j++) 
	{
	  tmpcm[j]/=masstot;
	  cm[j]=tmpcm[j];
	}
#ifdef DEBUG
      printf("findcenter: %f %f %f | %f\n", cm[0], cm[1], cm[2], maxdist);
#endif
      maxdist= (maxdist > extend) ? extend : maxdist;
      maxdist*=0.9;
      count=idum;
      if ((count<50) && (cvel[0]==0))
	{
	  masstot=0;
	  for (k=0; k< count; k++)
	    {
	      i = idlist[k];
	      for (j=0; j<3; j++) 
		{
		  cvel[j] += part[i].vel[j] * part[i].mass;
		}
	      masstot+= part[i].mass;
	    }
	  for (j=0; j<3; j++) cvel[j]=cvel[j]/(masstot);
	}
    } while (count> 5);
  for (i=0; i<cnt; i++)
    {
      for (j=0; j<3; j++)
	{
	  part[i].pos[j] -= cm[j];
	  if ( part[i].pos[j] > (BOXSIZE/2.) )
	    {
	      part[i].pos[j] = BOXSIZE - part[i].pos[j];
	    }
	  else if ( part[i].pos[j] < -(BOXSIZE/2.) )
	    {
	      part[i].pos[j] = BOXSIZE + part[i].pos[j];
	    }
	  part[i].vel[j] -= cvel[j];
	}
    }
  free(idlist);
  return;
}


#ifdef GSL
void rotategalaxy(gadpart *part, int numpart, double rad, int use, double *res, gsl_matrix **rotmat )
{
  int i,j,k,l, cnt=0;

  int *idlist = (int*) calloc(numpart, sizeof(int)); 
  if ( rad <= 0 )
    {
      gadpart_dist *wpart =(gadpart_dist*)  malloc(numpart * sizeof(gadpart_dist));
      double masstot = 0;
      double massdum = 0;
      for ( i = 0; i < numpart; i++ )
	{
	  if (! (( 1 << part[i].type) & use)) 
	    {
	      continue;
	    }
	  idlist[cnt] = i;
	  wpart[cnt].part = part[i];
	  wpart[cnt].dist = sqrt( SQR(part[i].pos[0]) + SQR(part[i].pos[1]) + SQR(part[i].pos[2]) );
	  //	  if (wpart[cnt].dist<24) masstot += part[i].mass;
	  cnt++;
	}
      qsort(wpart, cnt, sizeof(gadpart_dist), cmp_dist);
      for ( i = 0; i < cnt; i++ )
	{
	  rad = wpart[i].dist;
	  massdum += wpart[i].part.mass;
	  if (massdum >= (masstot/2.) ) break;
	}
      free(wpart);
    }
  else
    for ( i = 0; i < numpart; i++ )
      {
	if (! (( 1 << part[i].type) & use)) continue;
	idlist[cnt++] = i;
      }



  double s_old, ddum;
  double q=1, s=1;
  gsl_matrix *I = gsl_matrix_alloc (3, 3);
  gsl_vector *eval = gsl_vector_alloc (3);
  gsl_matrix *evec = gsl_matrix_alloc (3, 3);
  gsl_matrix *LU = gsl_matrix_alloc (3, 3);
  gsl_matrix *inv  = gsl_matrix_alloc (3, 3);     
  gsl_eigen_symmv_workspace * w =  gsl_eigen_symmv_alloc (3);
  gsl_matrix *rotation = gsl_matrix_alloc (3, 3);
  *rotmat = rotation;
  gsl_matrix_set_identity(rotation);
  gsl_matrix *resultmatrix = gsl_matrix_alloc (3, 3);
  gsl_matrix_set_zero(I);

     
  do
    {
      double dist;
      s_old=s;
      gsl_matrix_set_zero(I);
      for (l=0; l < cnt; l++)
	{     
	  k=idlist[l];
	  for (i=0; i < 3; i++)
	    for (j=i; j < 3; j++)
	      {
		ddum =part[k].pos[i] * part[k].pos[j];
		dist =SQR(part[k].pos[0])+SQR(part[k].pos[1]/q)+SQR(part[k].pos[2]/s);
		ddum/= dist;
		if (sqrt(dist)<rad) 
		  {
		    gsl_matrix_set(I,i,j, gsl_matrix_get(I,i,j)+ddum);
		    if (i!=j) gsl_matrix_set(I,j,i, gsl_matrix_get(I,j,i)+ddum);
		  }
	      }
	}

      gsl_eigen_symmv (I, eval, evec, w);
      gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
      gsl_matrix_memcpy(LU, evec);
     
      for (i=0; i < 3; i++)
	{
	  double eval_i  = gsl_vector_get (eval, i);
	  if (i==0) ddum = sqrt(eval_i);
	  else if (i==1) q=sqrt(eval_i)/ddum;
	  else           s=sqrt(eval_i)/ddum;
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
      for (k=0; k < numpart; k++)
	{
	  for (i=0; i<3; i++)
	    {
	      gsl_vector_set(oldpos, i, part[k].pos[i]);
	      gsl_vector_set(oldvel, i, part[k].vel[i]);
	    }
	  gsl_blas_dgemv( CblasNoTrans, 1.0, inv, oldpos, 0.0, newpos);
	  gsl_blas_dgemv( CblasNoTrans, 1.0, inv, oldvel, 0.0, newvel);
	  for (i=0; i<3; i++)
	    {
	      part[k].pos[i]=gsl_vector_get(newpos, i);
	      part[k].vel[i]=gsl_vector_get(newvel, i);
	    }
	}
      gsl_vector_free(oldpos);
      gsl_vector_free(newpos);
      gsl_vector_free(oldvel);
      gsl_vector_free(newvel);
      gsl_permutation_free(perm);
    } while ((ABS(s_old-s)/s) > 1e-2);
  res[0] = q;
  res[1] = s;
  gsl_matrix_free(I);
  gsl_matrix_free(resultmatrix);
  //  gsl_matrix_free(rotation);
  gsl_matrix_free(inv);
  gsl_matrix_free(LU);
  gsl_matrix_free(evec);
  gsl_vector_free(eval);
  gsl_eigen_symmv_free(w);
  free(idlist);
  return;
}
#endif //GSL
#endif //NOVEL


void simplecm(gadpart ** part, int cnt, float * cm)
{
  double masstot=0;
  int i, j;
  cm[0]=0;
  cm[1]=0;
  cm[2]=0;
  
  for ( i = 0; i < cnt; i++)
    {
      for ( j = 0; j < 3; j++) cm[j] += part[i] -> pos[j];
      masstot += part[i] -> mass;
    }
  for ( j = 0; j < 3; j++) cm[j] /= masstot;
}

double r200(gadpart_dist* pd, int cnt, double denscontrast, struct header h, int *vcnt, double *mvir)
{
  //  qsort(pd, cnt, sizeof(gadpart_dist), cmp_dist);
  double cdens;
  cdens=crit_dens*(h.omegal+h.omega0*pow(h.time,-3));
  double masstot=0;
  int i=0;
  double od=denscontrast+1;
  while ((od > denscontrast) || (i<10))
    {
      masstot+=pd[i].part.mass;
      od=masstot/(pow(pd[i].dist*h.time,3)*(4.0/3.0)*M_PI*cdens);
      i++;
      if (i>cnt) {i--;break;}
    }
  if (i==10) return 0;
  *vcnt=i;
  *mvir=masstot-pd[i].part.mass;
  return pd[i-1].dist;
}

double xoffset(gadpart_dist* pd, int cnt, float rvir, float* center)
{
  int i,j;
  double totmass=0;
  double cm[3] = {0., 0., 0.};
  for ( j = 0; j < 3; j++ )
    cm[j] = 0;
 
  for ( i = 0; i < cnt; i++ )
    {
      if (pd[i].dist > rvir) continue;
      totmass += pd[i].part.mass;
      for ( j = 0; j < 3; j++ )
      cm[j] += pd[i].part.pos[j] * pd[i].part.mass;
    }
  for ( j = 0; j < 3; j++ )
    cm[j] /= totmass; 

  double dist;
  float fcm[3];
  for (j = 0; j < 3; j++)
    fcm[j] = (float) cm[j];
  dist = distance(center, fcm);
  return dist / rvir;
}

#ifndef NOGAS
double temperature(const gadpart part)
{
  if (part.sph==NULL) return 0;
  //  const double k = 1.3806504e-23;
  const double gamma_minus_1 = 5.0/3 -1;
  const double boltzmann = 1.3806e-16;
  const double protonm= 1.6726e-24;
  const double hydr_frac= 0.76;
  const double yhelium = ( 1 - hydr_frac ) / ( 4 * hydr_frac );
  const double mass_in_g= 1.989e43;
  const double length_cm=3.085678e21;
  const double vel_in_cm_per_s=1e5;
  const double time_in_s= length_cm / vel_in_cm_per_s;
  const double energy_cgs = mass_in_g * pow(length_cm, 2) / pow(time_in_s, 2);
  const double mu = (1 + 4 * yhelium) / (1 + yhelium + part.sph->nelec);    
  
  double temp = gamma_minus_1 / boltzmann * part.sph->u * protonm * mu;
  temp *= energy_cgs / mass_in_g;
  
  return temp;
}
#endif //NOGAS

double galage(double z, double omegam, double omegal, double h)
{
  double omegak = 1.0 - omegam - omegal;
  double h0 = 3.24078e-18 * h;
  double tsec;

  if ((omegam < 0) || (omegam > 1)) return 0;

  if (omegam == 1)
    tsec = 2.0/(3.0*h0)/pow((1.0+z),1.5);
  else if (omegak < 1.0e-10)
    {
      double tmp1 = 2.0/(3.0*sqrt(1.0-omegam)*h0);
      double tmp2 = sqrt((1.0 - omegam)/omegam)/pow((1.0+z), 1.5);
      double tmp3 = log(tmp2 + sqrt(1.0+tmp2*tmp2));
      tsec = tmp1 * tmp3;
    } else if (omegal < 1.0e-10)
    {
      double tmp1 = 1.0/h0 * omegam/(2.0*pow((1.0-omegam),1.5));
      double tmp2 = 2.0*(1.0-omegam)/(omegam*(1.0+z))+1.0;
      double tmp3 = 2.0 * 
	log(sqrt((tmp2+1.0)/2.0) 
	    + sqrt((tmp2-1.0)/2.0));

      tsec = tmp1 * ( sqrt(pow(tmp2,2)-1.0) - tmp3 );

    } else return 0;

  return tsec / 3.155815e7;
}

double timediff(double z1, double z2, double omegam, double omegal, double h)
{
  double time1 = galage (z1, omegam, omegal, h);
  double time2 = galage (z2, omegam, omegal, h);
  double diff = time1 - time2;
  if (diff < 0) return ((-1) * diff);
  else return diff;
}

double a2z(double a)
{
  return ((1 / a) - 1);
}

double z2a(double z)
{
  return (1 / (z + 1));
}
#ifdef GSL
#ifndef NOVEL
void rotatepart(gadpart *part, int numpart, const gsl_matrix *rotmat)
{
  int i,j,k; 
  gsl_vector *oldpos = gsl_vector_alloc (3);
  gsl_vector *newpos = gsl_vector_alloc (3);  
  gsl_vector *oldvel = gsl_vector_alloc (3);
  gsl_vector *newvel = gsl_vector_alloc (3);
  for (k=0; k < numpart; k++)	
    {
      for (i=0; i<3; i++)
	{
	  gsl_vector_set(oldpos, i, part[k].pos[i]);          
	  gsl_vector_set(oldvel, i, part[k].vel[i]);
	}
      gsl_blas_dgemv( CblasNoTrans, 1.0, rotmat, oldpos, 0.0, newpos);
      gsl_blas_dgemv( CblasNoTrans, 1.0, rotmat, oldvel, 0.0, newvel);
      //ddum += gsl_blas_dnrm2 (oldvel) - gsl_blas_dnrm2 (newvel);
      for (i=0; i<3; i++)
	{
	  part[k].pos[i]=gsl_vector_get(newpos, i);
	  part[k].vel[i]=gsl_vector_get(newvel, i);
	}
      
    }

}
#endif //NOVEL

void xrotate(double angle, gadpart *part, int numpart)
{
  double angle_rad = (angle / 180) * M_PI;
  gsl_matrix *rotmat = gsl_matrix_alloc (3, 3);
  gsl_matrix_set_zero(rotmat);
  gsl_matrix_set(rotmat,0,0,  1.0);
  gsl_matrix_set(rotmat,1,1,  cos(angle_rad));
  gsl_matrix_set(rotmat,1,2,  sin(angle_rad));
  gsl_matrix_set(rotmat,2,2,  cos(angle_rad));
  gsl_matrix_set(rotmat,2,1, -sin(angle_rad));
  rotatepart(part, numpart, rotmat);
  gsl_matrix_free(rotmat);
}

void yrotate(double angle, gadpart *part, int numpart)
{
  double angle_rad = (angle / 180) * M_PI;
  gsl_matrix *rotmat = gsl_matrix_alloc (3, 3);
  gsl_matrix_set_zero(rotmat);
  gsl_matrix_set(rotmat,1,1,  1.0);
  gsl_matrix_set(rotmat,0,0,  cos(angle_rad));
  gsl_matrix_set(rotmat,2,0,  sin(angle_rad));
  gsl_matrix_set(rotmat,2,2,  cos(angle_rad));
  gsl_matrix_set(rotmat,0,2, -sin(angle_rad));
  rotatepart(part, numpart, rotmat);
  gsl_matrix_free(rotmat);
}


void zrotate(double angle, gadpart *part, int numpart)
{
  double angle_rad = (angle / 180) * M_PI;
  gsl_matrix *rotmat = gsl_matrix_alloc (3, 3);
  gsl_matrix_set_zero(rotmat);
  gsl_matrix_set(rotmat,2,2,  1.0);
  gsl_matrix_set(rotmat,0,0,  cos(angle_rad));
  gsl_matrix_set(rotmat,1,1,  cos(angle_rad));
  gsl_matrix_set(rotmat,0,1,  sin(angle_rad));
  gsl_matrix_set(rotmat,1,0, -sin(angle_rad));
  rotatepart(part, numpart, rotmat);
  gsl_matrix_free(rotmat);
}
#endif //GSL

double angle(fltarr a, fltarr b)
{
  double anorm = sqrt( SQR(a[0]) + SQR(a[1]) + SQR(b[2]) );
  double bnorm = sqrt( SQR(b[0]) + SQR(b[1]) + SQR(b[2]) );
  double adotb = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  return acos(adotb / (anorm * bnorm));
}

double radvel(fltarr vel, fltarr rad)
{
  //  double velnorm = sqrt( SQR(vel[0]) + SQR(vel[1]) + SQR(vel[2]) );
  double radnorm = sqrt( SQR(rad[0]) + SQR(rad[1]) + SQR(rad[2]) );
  double vdotr = vel[0] * rad[0] + vel[1] * rad[1] + vel[2] * rad[2];
  //  double angle = acos(vdotr / (velnorm * radnorm));
  return (vdotr / radnorm);
  
}


struct header cphead(struct header head, gadpart* part, int cnt)
{
  struct header res = head;
  int i;
  for ( i=0; i<6; i++)
    {
      res.npart[i] = 0;
      res.nall[i] = 0;
    }

  for ( i=0; i<cnt; i++)
    {
      int t = part[i].type;
      res.npart[t]++;
      res.nall[t]++;
    }
  res.numfiles=1;
  return res;
}
