/*
Program to do retrieve fof-group-finder catalogue

gcc -lm -o get_group_catalogue  get_group_catalogue.c
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
//#include "libgad.h"

typedef struct
{
  int len;
  int file;
  int offset;
  int typecount[6];
  double typemass[6];
  float center[3];
  float vel[3];
  float mass;
  float sfr;
  float mdot;
  float mbh;
}
cat_data;

void usage()
{
  	  fprintf(stderr," Get_Group_Catalogue v0.01\n");
  	  fprintf(stderr," -i <input file name>\n");
  	  fprintf(stderr," -n <snapshot number>\n");
  	  fprintf(stderr," \n\n");
	  exit(1);
}

int id_sort_groups(const void *a, const void *b)
{
  if(((cat_data *) a)->len > ((cat_data *) b)->len)
    return -1;

  if(((cat_data *) a)->len < ((cat_data *) b)->len)
    return +1;

  return 0;
}

int id_sort_compare_key(const void *a, const void *b)
{
  if(*((unsigned long *) a) < *((unsigned long *) b))
    return -1;

  if(*((unsigned long *) a) > *((unsigned long *) b))
    return +1;

  return 0;
}


int get_group_catalogue(int Num, char *OutputDir, cat_data **ret_cat, int *ngroups)
{
  int Ngroups, Nids, TotNgroups, NTask, count, i, j, filenr;
  int *GroupLen, *GroupOff, *GroupFileNr, *GroupNr, *GroupTypeLen;
  long long TotNids;
  float *GroupTypeMass;
  float *GroupCenter, *GroupVel, *GroupSfr, *GroupMass;
  char buf[1000];
  cat_data *cat;
  FILE *fd;

  //  Num = *(int *) argv[1];


  sprintf(buf, "%s/groups_%03d/group_tab_%03d.0", OutputDir, Num, Num);
  if(!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/group_tab_%03d.0", OutputDir, Num);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  return -1;
	}
    }

  fread(&Ngroups, sizeof(int), 1, fd);
  fread(&TotNgroups, sizeof(int), 1, fd);
  *ngroups = TotNgroups;
  fread(&Nids, sizeof(int), 1, fd);
  fread(&TotNids, sizeof(long long), 1, fd);
  fread(&NTask, sizeof(int), 1, fd);
  fclose(fd);

  GroupLen = (int *) malloc(sizeof(int) * TotNgroups);
  GroupOff = (int *) malloc(sizeof(int) * TotNgroups);
  GroupFileNr = (int *) malloc(sizeof(int) * TotNgroups);
  GroupNr = (int *) malloc(sizeof(int) * TotNgroups);
  GroupTypeLen = (int *) malloc(sizeof(int) * TotNgroups * 6);
  GroupTypeMass = (float *) malloc(sizeof(float) * TotNgroups * 6);
  GroupCenter = (float *) malloc(sizeof(float) * TotNgroups * 3);
  GroupVel = (float *) malloc(sizeof(float) * TotNgroups * 3);
  GroupMass = (float *) malloc(sizeof(float) * TotNgroups);
  GroupSfr = (float *) malloc(sizeof(float) * TotNgroups);
  


  for(filenr = 0, count = 0; filenr < NTask; filenr++)
    {
      //      printf("filenr %d | NTask %d \n", filenr, NTask);
      sprintf(buf, "%s/groups_%03d/group_tab_%03d.%d", OutputDir, Num, Num, filenr);
      if(!(fd = fopen(buf, "r")))
	{
	  sprintf(buf, "%s/group_tab_%03d.%d", OutputDir, Num, filenr);
	  if(!(fd = fopen(buf, "r")))
	    {
	      printf("can't open file `%s'\n", buf);
	      return -1;
	    }
	}

      fread(&Ngroups, sizeof(int), 1, fd);
      fread(&TotNgroups, sizeof(int), 1, fd);
      fread(&Nids, sizeof(int), 1, fd);
      fread(&TotNids, sizeof(long long), 1, fd);
      fread(&NTask, sizeof(int), 1, fd);
      //      printf(" %d %d %d %d | %d\n", Ngroups, Nids, TotNgroups, NTask, TotNids);
      

      fread(&GroupLen[count], sizeof(int), Ngroups, fd);
      /* skip offset table */
      fread(&GroupOff[count], sizeof(int), Ngroups, fd);
      fread(&GroupMass[count], sizeof(float), Ngroups, fd);
      fread(&GroupCenter[3*count], 3*sizeof(float), Ngroups, fd);
      fread(&GroupVel[3*count], 3*sizeof(float), Ngroups, fd);
      fread(&GroupTypeLen[6*count], 6*sizeof(int), Ngroups, fd);
      fread(&GroupTypeMass[6*count], 6*sizeof(float), Ngroups, fd);
      fread(&GroupSfr[count], sizeof(float), Ngroups, fd);
      for(i = 0; i < Ngroups; i++)
	{
	  GroupFileNr[i + count] = filenr;
	  GroupNr[i + count] = i;
	}
      
      count += Ngroups;
      
      fclose(fd);
    }

  cat = calloc(TotNgroups, sizeof(cat_data));
  *ret_cat = cat;

//for(i = 0; i < 10; i++)
//  printf("%d\n", GroupLen[i]);

for(i = 0; i < TotNgroups; i++)
    {
      cat[i].len = GroupLen[i];
      cat[i].file = GroupFileNr[i];
      cat[i].offset = GroupOff[i];
      cat[i].mass = GroupMass[i];
      for(j=0; j<6; j++)
	{
	  cat[i].typecount[j] = GroupTypeLen[6*i+j];
	  cat[i].typemass[j] = GroupTypeMass[6*i+j];
	}
      for(j=0; j<3; j++)
	{
	  cat[i].center[j] = GroupCenter[3*i+j];
	  cat[i].vel[j] = GroupVel[3*i+j];
	}
      cat[i].sfr = GroupSfr[i];
    }
  qsort(cat, TotNgroups, sizeof(cat_data), id_sort_groups);

//  for(i = 0; i < TotNgroups; i++)
//    {
//      GroupLen[i] = cat[i].len;
//      GroupFileNr[i] = cat[i].file;
//      GroupNr[i] = cat[i].offset;
//      for(j=0; j<6; j++)
//	{
//	  GroupTypeLen[6*i+j] = cat[i].typecount[j];
//	  GroupTypeMass[6*i+j] = cat[i].typemass[j];
//	}
//
//      for(j=0; j<3; j++)
//	GroupCenter[3*i+j]= cat[i].center[j];
//      GroupSfr[i] = cat[i].sfr;
//    }
  free(GroupLen);
  free(GroupOff);
  free(GroupFileNr); 
  free(GroupNr);
  free(GroupTypeLen);
  free(GroupTypeMass);
  free(GroupCenter);
  free(GroupVel);
  free(GroupMass);
  free(GroupSfr);

  //  free(cat);

  return TotNgroups;
}



int get_group_indices(int Num, char *OutputDir, cat_data cat, unsigned long **IDS)
{
  int FileNr = cat.file;
  int offset = cat.offset;
  int GroupLen = cat.len;
  int NTask;
  FILE *fd;
  long long TotNids;
  int Ngroups, TotNgroups, Nids, idum;
  char buf[1024];
  unsigned long *ids;

  ids = malloc(GroupLen * sizeof(unsigned long));
  *IDS = ids;

  sprintf(buf, "%s/groups_%03d/group_ids_%03d.%d", OutputDir, Num, Num, FileNr);
  if(!(fd = fopen(buf, "r")))
    {
      sprintf(buf, "%s/group_ids_%03d.%d", OutputDir, Num, FileNr);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  return -1;
	}
    }
  fread(&Ngroups, sizeof(int), 1, fd);
  fread(&TotNgroups, sizeof(int), 1, fd);
  fread(&Nids, sizeof(int), 1, fd);
  fread(&TotNids, sizeof(long long), 1, fd);
  fread(&NTask, sizeof(int), 1, fd);
  fread(&idum, sizeof(int), 1, fd);
  fseek(fd, sizeof(unsigned long) * offset, SEEK_CUR);
  fread(ids, sizeof(unsigned long), GroupLen, fd);
  fclose(fd);
//  printf("%lu, %d, %d\n", ids[0], offset, sizeof(unsigned long long) );
  qsort(ids, GroupLen, sizeof(unsigned long), id_sort_compare_key);
//  printf("%lu\n", ids[0]);

  return 0;
}

void write_ids(char *filename, unsigned long *ids, cat_data cat)
{
  FILE *fp = fopen(filename, "w");
  unsigned int dum=0;
  fwrite(&dum, sizeof(unsigned int), 1, fp);
  fwrite(&(cat.center), sizeof(float), 3, fp);
  //  fwrite(ids, sizeof(unsigned long), len,fp);
  fclose(fp);
}

int main (int argc, char *argv[])
{
  int i;
  int ngroups=0;
  int num = 95;
  int write = 0;
  char infile[1024], outbase[1024];
  cat_data *cat;
  int min_len = 0;
  i=1;
  strcpy(infile,"./");
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
      else if (!strcmp(argv[i],"-n")) 
	{
	  i++;
	  num = atoi(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-o")) 
	{
	  i++;
	  strcpy(outbase, argv[i]);
	  write = 1;
	  i++;
	}
      else if (!strcmp(argv[i],"-m")) 
	{
	  i++;
	  min_len= atoi(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-use")) {
      } else {
	usage();
      }
    }

  get_group_catalogue(num, infile, &cat, &ngroups);
  unsigned int j=0;
    for ( i = 0; i < ngroups; i++) 
  //  for ( i = 0; i < 1; i++) 
  //  for ( i = 0; i < 3; i++) 
    {
//      printf("%d %10d %10d %10g\n", i, cat[i].len, cat[i].offset, cat[i].mass);
//      printf("%10d %d %d %d %d %d %d %d \n", i, cat[i].len, cat[i].typecount[0], cat[i].typecount[1], cat[i].typecount[2], cat[i].typecount[3], cat[i].typecount[4], cat[i].typecount[5]);
//      printf("%10d %10g %10g %10g %10g %10g %10g \n", i, cat[i].typemass[0], cat[i].typemass[1], cat[i].typemass[2], cat[i].typemass[3], cat[i].typemass[4], cat[i].typemass[5]);
//      printf("%10d %10g %10g %10g \n", i, cat[i].center[0], cat[i].center[1], cat[i].center[2]);
//      printf("%10d %10g %10g %10g \n", i, cat[i].vel[0], cat[i].vel[1], cat[i].vel[2]);
//      printf("%10d %10g %10g %10g \n", i, cat[i].sfr, cat[i].mdot, cat[i].mbh);

      printf("%6d %10d %10d %10d %10g %10g %10g %10g\n", i, cat[i].len, cat[i].offset, cat[i].file, cat[i].mass, cat[i].center[0], cat[i].center[1], cat[i].center[2]);
      if ( (write) && (cat[i].len >= min_len) )
	{      	  
	  unsigned long *ids;
	  get_group_indices(num, infile, cat[i], &ids);
	  char outfilename[1024];
	  sprintf(outfilename,"%s_%u", outbase, j++);
	  write_ids(outfilename, ids, cat[i]);
	  free(ids);
	}
    }

    
  return 0;
}
