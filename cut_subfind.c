/*

  stan: gcc -lm -lgad-stan -lgsl -lgslcblas cut_subfind.c -o ~/bin/cut_subfind

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "libgad.h"

#define FLAG_Group_VelDisp 0

int         Ngroups;
int         TotNgroups;
int         Nids;
long long   TotNids;
int         NTask;
int         Nsubgroups;
int         TotNsubgroups;
int         Offset;

int    *GroupLen, *GroupOffset, *GroupContaminationCount;
int    *GroupNsubs, *GroupFirstSub;
double *GroupMass, *GroupPosX, *GroupPosY, *GroupPosZ, *Group_M_Mean200, *Group_R_Mean200;
double *Group_M_Crit200, *Group_R_Crit200, *Group_M_TopHat200, *Group_R_TopHat200;
double *Group_VelDisp_Mean200, *Group_VelDisp_Crit200, *Group_VelDisp_TopHat200;
double *GroupContaminationMass;

int    *SubhaloLen, *SubhaloOffset, *SubhaloParent, *SubhaloIDMostbound, *SubhaloGrNr;
double *SubhaloMass, *SubhaloVelDisp, *SubhaloVmax, *SubhaloVmaxRad, *SubhaloHalfmassRad;
double *SubhaloPosX, *SubhaloPosY, *SubhaloPosZ;
double *SubhaloVelX, *SubhaloVelY, *SubhaloVelZ;
double *SubhaloCMX, *SubhaloCMY, *SubhaloCMZ;
double *SubhaloSpinX, *SubhaloSpinY, *SubhaloSpinZ;
double *SubhaloGasMass, *SubhaloHaloMass, *SubhaloDiscMass, *SubhaloBulgeMass, *SubhaloStarsMass, *SubhaloBndryMass;

int    *IDs;

int skip = 0;
int skip_sub = 0;

char filename[256];
unsigned int numpart;
struct header head;
gadpart *part;
float StarMassThreshold = 0.0;
float TotalMassThreshold = 0.0;
int use_cm = 16;
int use_ind = 0;

void usage()
{
  fprintf(stderr,"Cut SubFind Groups out of a GADGET Snapshot v0.2\n");
  fprintf(stderr,"-f <snapshot filename>\n");
  fprintf(stderr,"-s set threshold in stellar mass in code units\n");
  fprintf(stderr,"-m set threshold in total mass in code units\n");
  fprintf(stderr,"-use bitcode for particle types used for centering\n");
  fprintf(stderr,"-ind <use array index instead of particle ID>\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"\n");
  exit(1);
}


int main(int argc, char **argv)
{
	int  i, j, snapshot_number;	
	
	/* Read input and output files */
	i=1;
	if (argc==1) usage();
	while (i<argc)
	  {
      
	    if (!strcmp(argv[i],"-f"))
	      {
		i++;
		strcpy(filename,argv[i]);
		i++;
	      }
	    else if (*argv[i]!='-')
	      {
		strcpy(filename,argv[i]);
		i++;
	      }
	    else if (!strcmp(argv[i],"-i"))
	      {
		i++;
		strcpy(filename,argv[i]);
		i++;
	      }
	    else if (!strcmp(argv[i],"-use"))
	      {
		i++;
		use_cm=atoi(argv[i]);
		i++;
	      }
	    else if (!strcmp(argv[i],"-ind"))
	      {
		i++;
		use_ind = 1;
	      }
	    else if (!strcmp(argv[i],"-s"))
	      {		
		i++;
		StarMassThreshold = atof(argv[i]);
		i++;
	      }
	    else if (!strcmp(argv[i],"-m"))
	      {		
		i++;
		TotalMassThreshold = atof(argv[i]);
		i++;
	      }
	    else usage();
	  }
	char snap_number_dum[4];
	sprintf(snap_number_dum, "%.*s", 3, &filename[ strlen(filename)-3]);
	snapshot_number = atoi( snap_number_dum );
	printf("num %d\n", snapshot_number);

	load_groups(snapshot_number);
	
	load_group_ids(snapshot_number);

	load_snapshot(filename, snapshot_number);
	
	write_ascii(snapshot_number);
	
	return 0;
}


int write_ascii(int snapshot_number){
	
	FILE *ofpp;
	char property_fname[256];
	int i,j,k;
	
	sprintf(property_fname, "subhalos_%03d/SFproperties.%03d",snapshot_number, snapshot_number);
	ofpp = fopen(property_fname,"w");
	fprintf(ofpp,"# %d\n",TotNsubgroups);
	fprintf(ofpp,"ID ParentID GroupNr PosX PosY PosZ HalfMassRad TotalMass GasMass HaloMass DiscMass BulgeMass StarMass BndrMass IDmostbound\n");
	for (i=0; i<TotNsubgroups; i++) fprintf(ofpp,"%d %d %d %g %g %g %g %g %g %g %g %g %g %g %d\n",i, SubhaloParent[i], SubhaloGrNr[i],SubhaloPosX[i],SubhaloPosY[i],SubhaloPosZ[i], SubhaloHalfmassRad[i], SubhaloMass[i],SubhaloGasMass[i],SubhaloHaloMass[i],SubhaloDiscMass[i],SubhaloBulgeMass[i],SubhaloStarsMass[i],SubhaloBndryMass[i], SubhaloIDMostbound[i]);
	fclose(ofpp);
	
	
	return 0;
}

int load_snapshot(char* fname, int num)
{
  int i,j,k,l;
  char outfilename[256];
  printf("reading file: %s\n", fname);
  numpart = readgadget_part(fname, &head, &part);   
  if (numpart==0) 
    {
      extern int libgaderr;
      fprintf(stderr,"LibGad Error Code: %d\n", libgaderr);
      exit(1);
    }
  if (!use_ind)
    qsort(part, numpart, sizeof(gadpart),  cmp_id );
  i=0;

  char buf[256];
  sprintf(buf, "subhalos_%03d", num);
  mkdir(buf, 02755);

  for (j=0; j<TotNsubgroups; j++) 
    {
      if ((StarMassThreshold) && (SubhaloStarsMass[j] < StarMassThreshold )) continue;
      if ((TotalMassThreshold) && (SubhaloMass[j] < TotalMassThreshold )) continue;
      gadpart *wpart=(gadpart*) malloc (sizeof(gadpart) * SubhaloLen[j]);
      struct header outhead = head;
      for (i=0; i<6; i++)
	{
	  outhead.nall[i] = 0;
	  outhead.npart[i]= 0;
	}
      l=0;
      for (k=SubhaloOffset[j];k<SubhaloOffset[j]+SubhaloLen[j];k++)
	{   
	  gadpart dumpart;
	  dumpart.id = IDs[k];
	  gadpart *fnd;
	  if (!use_ind)
	    {
	      fnd = bsearch (&dumpart, &part[0], numpart, sizeof(gadpart), cmp_id );
	    }
	  else
	    {
	      fnd = &part[IDs[k]-1];
	    }
	  if (fnd != NULL)
	    {
	      wpart[l] = *fnd;
	      outhead.nall[ wpart[l].type ]++;
	      outhead.npart[ wpart[l].type ]++;
	      l++;
	    }	  
	}
      findcenter(wpart, SubhaloLen[j], -1, use_cm);
      sprintf(outfilename, "subhalos_%03d/subhalo.%05d.p%05d.g%05d.gad", num, j, SubhaloParent[j],SubhaloGrNr[j]);
      writegadget_part(outfilename, outhead, wpart);
      free (wpart);
    }
}



int load_groups(int snapshot_number)
{
	FILE *fd;
	char   buf[200];
	int    i,j,k,fnr;
	
	int   *locLen, *locOffset, *locContaminationCount;
	int   *locNsubs, *locFirstSub;
	float *locMass, *locPosX, *locPosY, *locPosZ, *loc_M_Mean200, *loc_R_Mean200;
	float *loc_M_Crit200, *loc_R_Crit200, *loc_M_TopHat200, *loc_R_TopHat200;
	float *loc_VelDisp_Mean200, *loc_VelDisp_Crit200, *loc_VelDisp_TopHat200;
	float *locContaminationMass;
	
	int   *locParent, *locIDMostbound, *locGrNr;
	float *locVelDisp, *locVmax, *locVmaxRad, *locHalfmassRad;
	float *locVelX, *locVelY, *locVelZ;
	float *locCMX, *locCMY, *locCMZ;
	float *locSpinX, *locSpinY, *locSpinZ;
	float *locGasMass, *locHaloMass, *locDiscMass, *locBulgeMass, *locStarsMass, *locBndryMass;
	
	fnr = 0;
		
	do {
		
		sprintf(buf, "groups_%03d/subhalo_tab_%03d.%d",snapshot_number,snapshot_number,fnr);		

		if(!(fd=fopen(buf,"r")))
		{
			printf("can't open file `%s`\n",buf);
			exit(0);
		}
		
		printf("reading `%s' ...\n",buf); fflush(stdout);
		
		//fread(&dummy, sizeof(dummy), 1, fd);
		//fread(&header, sizeof(header), 1, fd);
		fread(&Ngroups, sizeof(int), 1, fd);
		fread(&TotNgroups, sizeof(int), 1, fd);
		fread(&Nids, sizeof(int), 1, fd);
		fread(&TotNids, sizeof(long long), 1, fd);
		fread(&NTask, sizeof(int), 1, fd);
		fread(&Nsubgroups, sizeof(int), 1, fd);
		fread(&TotNsubgroups, sizeof(int), 1, fd);	
		
		/*
		 printf("%d\n",Ngroups);
		 printf("%d\n",TotNgroups);
		 printf("%d\n",Nids);
		 printf("%lld\n",TotNids);
		 printf("%d\n",NTask);
		 printf("%d\n",Nsubgroups);
		 printf("%d\n",TotNsubgroups);
		 */
		
		if (fnr==0){
			
			GroupLen = (int*) calloc (TotNgroups, sizeof(int));
			GroupOffset = (int*) calloc (TotNgroups, sizeof(int));
			GroupMass = (double*) calloc (TotNgroups, sizeof(double));
			GroupPosX = (double*) calloc (TotNgroups, sizeof(double));
			GroupPosY = (double*) calloc (TotNgroups, sizeof(double));
			GroupPosZ = (double*) calloc (TotNgroups, sizeof(double));
			Group_M_Mean200 = (double*) calloc (TotNgroups, sizeof(double));
			Group_R_Mean200 = (double*) calloc (TotNgroups, sizeof(double));
			Group_M_Crit200 = (double*) calloc (TotNgroups, sizeof(double));
			Group_R_Crit200 = (double*) calloc (TotNgroups, sizeof(double));
			Group_M_TopHat200 = (double*) calloc (TotNgroups, sizeof(double));
			Group_R_TopHat200 = (double*) calloc (TotNgroups, sizeof(double));
			Group_VelDisp_Mean200 = (double*) calloc (TotNgroups, sizeof(double));
			Group_VelDisp_Crit200 = (double*) calloc (TotNgroups, sizeof(double));
			Group_VelDisp_TopHat200 = (double*) calloc (TotNgroups, sizeof(double));
			GroupContaminationCount = (int*) calloc (TotNgroups, sizeof(int));
			GroupContaminationMass= (double*) calloc (TotNgroups, sizeof(double));
			GroupNsubs = (int*) calloc (TotNgroups, sizeof(int));
			GroupFirstSub = (int*) calloc (TotNgroups, sizeof(int));
			
			SubhaloLen = (int*) calloc (TotNsubgroups, sizeof(int));
			SubhaloOffset = (int*) calloc (TotNsubgroups, sizeof(int));
			SubhaloParent = (int*) calloc (TotNsubgroups, sizeof(int));
			SubhaloMass = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloPosX = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloPosY = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloPosZ = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloVelX = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloVelY = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloVelZ = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloCMX = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloCMY = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloCMZ = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloSpinX = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloSpinY = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloSpinZ = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloVelDisp = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloVmax = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloVmaxRad = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloHalfmassRad = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloIDMostbound = (int*) calloc (TotNsubgroups, sizeof(int));
			SubhaloGrNr = (int*) calloc (TotNsubgroups, sizeof(int));
			SubhaloGasMass = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloHaloMass = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloDiscMass = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloBulgeMass = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloStarsMass = (double*) calloc (TotNsubgroups, sizeof(double));
			SubhaloBndryMass = (double*) calloc (TotNsubgroups, sizeof(double));
		}
		
		if (Ngroups>0){
			
			locLen = (int*) calloc (Ngroups, sizeof(int));
			locOffset = (int*) calloc (Ngroups, sizeof(int));
			locMass = (float*) calloc (Ngroups, sizeof(float));
			locPosX = (float*) calloc (Ngroups, sizeof(float));
			locPosY = (float*) calloc (Ngroups, sizeof(float));
			locPosZ = (float*) calloc (Ngroups, sizeof(float));
			loc_M_Mean200 = (float*) calloc (Ngroups, sizeof(float));
			loc_R_Mean200 = (float*) calloc (Ngroups, sizeof(float));
			loc_M_Crit200 = (float*) calloc (Ngroups, sizeof(float));
			loc_R_Crit200 = (float*) calloc (Ngroups, sizeof(float));
			loc_M_TopHat200 = (float*) calloc (Ngroups, sizeof(float));
			loc_R_TopHat200 = (float*) calloc (Ngroups, sizeof(float));
			loc_VelDisp_Mean200 = (float*) calloc (Ngroups, sizeof(float));
			loc_VelDisp_Crit200 = (float*) calloc (Ngroups, sizeof(float));
			loc_VelDisp_TopHat200 = (float*) calloc (Ngroups, sizeof(float));
			locContaminationCount = (int*) calloc (Ngroups, sizeof(int));
			locContaminationMass= (float*) calloc (Ngroups, sizeof(float));
			locNsubs = (int*) calloc (Ngroups, sizeof(int));
			locFirstSub = (int*) calloc (Ngroups, sizeof(int));
			
			//printf("Ngroups = %d\n",Ngroups);
			
			for (i=0;i<Ngroups;i++) fread(&locLen[i], sizeof(int), 1, fd);
			for (i=0;i<Ngroups;i++) fread(&locOffset[i], sizeof(int), 1, fd);
			for (i=0;i<Ngroups;i++) fread(&locMass[i], sizeof(float), 1, fd);
			for (i=0;i<Ngroups;i++){
				fread(&locPosX[i], sizeof(float), 1, fd);
				fread(&locPosY[i], sizeof(float), 1, fd);
				fread(&locPosZ[i], sizeof(float), 1, fd);
			}
			for (i=0;i<Ngroups;i++) fread(&loc_M_Mean200[i], sizeof(float), 1, fd);
			for (i=0;i<Ngroups;i++) fread(&loc_R_Mean200[i], sizeof(float), 1, fd);
			for (i=0;i<Ngroups;i++) fread(&loc_M_Crit200[i], sizeof(float), 1, fd);
			for (i=0;i<Ngroups;i++) fread(&loc_R_Crit200[i], sizeof(float), 1, fd);
			for (i=0;i<Ngroups;i++) fread(&loc_M_TopHat200[i], sizeof(float), 1, fd);
			for (i=0;i<Ngroups;i++) fread(&loc_R_TopHat200[i], sizeof(float), 1, fd);
			if (FLAG_Group_VelDisp){
				for (i=0;i<Ngroups;i++) fread(&loc_VelDisp_Mean200[i], sizeof(float), 1, fd);
				for (i=0;i<Ngroups;i++) fread(&loc_VelDisp_Crit200[i], sizeof(float), 1, fd);
				for (i=0;i<Ngroups;i++) fread(&loc_VelDisp_TopHat200[i], sizeof(float), 1, fd);
			}
			for (i=0;i<Ngroups;i++) fread(&locContaminationCount[i], sizeof(int), 1, fd);
			for (i=0;i<Ngroups;i++) fread(&locContaminationMass[i], sizeof(float), 1, fd);
			for (i=0;i<Ngroups;i++) fread(&locNsubs[i], sizeof(int), 1, fd);
			for (i=0;i<Ngroups;i++) fread(&locFirstSub[i], sizeof(int), 1, fd);
			
			//for (i=0;i<Ngroups;i++) printf("%f %f %f\n",locPosX[i],locPosY[i],locPosZ[i]);
			//for (i=0;i<Ngroups;i++) printf("%d\n",locFirstSub[i]);
			
			for (i=0;i<Ngroups;i++){
				GroupLen[i+skip] = locLen[i];
				GroupOffset[i+skip] = locOffset[i];
				GroupMass[i+skip] = locMass[i];
				GroupPosX[i+skip] = locPosX[i];
				GroupPosY[i+skip] = locPosY[i];
				GroupPosZ[i+skip] = locPosZ[i];
				Group_M_Mean200[i+skip] = loc_M_Mean200[i];
				Group_R_Mean200[i+skip] = loc_R_Mean200[i];
				Group_M_Crit200[i+skip] = loc_M_Crit200[i];
				Group_R_Crit200[i+skip] = loc_R_Crit200[i];
				Group_M_TopHat200[i+skip] = loc_M_TopHat200[i];
				Group_R_TopHat200[i+skip] = loc_R_TopHat200[i];
				Group_VelDisp_Mean200[i+skip] = loc_VelDisp_Mean200[i];
				Group_VelDisp_Crit200[i+skip] = loc_VelDisp_Crit200[i];
				Group_VelDisp_TopHat200[i+skip] = loc_VelDisp_TopHat200[i];
				GroupContaminationCount[i+skip] = locContaminationCount[i];
				GroupContaminationMass[i+skip] = locContaminationMass[i];
				GroupNsubs[i+skip] = locNsubs[i];
				GroupFirstSub[i+skip] = locFirstSub[i];			
			}
			
			skip+=Ngroups;
			
			free(locLen);
			free(locOffset);
			free(locMass);
			free(locPosX);
			free(locPosY);
			free(locPosZ);
			free(loc_M_Mean200);
			free(loc_R_Mean200);
			free(loc_M_Crit200);
			free(loc_R_Crit200);
			free(loc_M_TopHat200);
			free(loc_R_TopHat200);
			free(loc_VelDisp_Mean200);
			free(loc_VelDisp_Crit200);
			free(loc_VelDisp_TopHat200);
			free(locContaminationCount);
			free(locContaminationMass);
			free(locNsubs);
			free(locFirstSub);
			
		}
		
		
		if (Nsubgroups>0){
			
			locLen = (int*) calloc (Nsubgroups, sizeof(int));
			locOffset = (int*) calloc (Nsubgroups, sizeof(int));
			locParent = (int*) calloc (Nsubgroups, sizeof(int));
			locMass = (float*) calloc (Nsubgroups, sizeof(float));
			locPosX = (float*) calloc (Nsubgroups, sizeof(float));
			locPosY = (float*) calloc (Nsubgroups, sizeof(float));
			locPosZ = (float*) calloc (Nsubgroups, sizeof(float));
			locVelX = (float*) calloc (Nsubgroups, sizeof(float));
			locVelY = (float*) calloc (Nsubgroups, sizeof(float));
			locVelZ = (float*) calloc (Nsubgroups, sizeof(float));
			locCMX = (float*) calloc (Nsubgroups, sizeof(float));
			locCMY = (float*) calloc (Nsubgroups, sizeof(float));
			locCMZ = (float*) calloc (Nsubgroups, sizeof(float));
			locSpinX = (float*) calloc (Nsubgroups, sizeof(float));
			locSpinY = (float*) calloc (Nsubgroups, sizeof(float));
			locSpinZ = (float*) calloc (Nsubgroups, sizeof(float));
			locVelDisp = (float*) calloc (Nsubgroups, sizeof(float));
			locVmax = (float*) calloc (Nsubgroups, sizeof(float));
			locVmaxRad = (float*) calloc (Nsubgroups, sizeof(float));
			locHalfmassRad = (float*) calloc (Nsubgroups, sizeof(float));
			locIDMostbound = (int*) calloc (Nsubgroups, sizeof(int));
			locGrNr = (int*) calloc (Nsubgroups, sizeof(int));
			locGasMass = (float*) calloc (Nsubgroups, sizeof(float));
			locHaloMass = (float*) calloc (Nsubgroups, sizeof(float));
			locDiscMass = (float*) calloc (Nsubgroups, sizeof(float));
			locBulgeMass = (float*) calloc (Nsubgroups, sizeof(float));
			locStarsMass = (float*) calloc (Nsubgroups, sizeof(float));
			locBndryMass = (float*) calloc (Nsubgroups, sizeof(float));
			
			for (i=0;i<Nsubgroups;i++) fread(&locLen[i], sizeof(int), 1, fd);
			for (i=0;i<Nsubgroups;i++) fread(&locOffset[i], sizeof(int), 1, fd);
			for (i=0;i<Nsubgroups;i++) fread(&locParent[i], sizeof(int), 1, fd);
			for (i=0;i<Nsubgroups;i++) fread(&locMass[i], sizeof(float), 1, fd);
			for (i=0;i<Nsubgroups;i++){
				fread(&locPosX[i], sizeof(float), 1, fd);
				fread(&locPosY[i], sizeof(float), 1, fd);
				fread(&locPosZ[i], sizeof(float), 1, fd);
			}
			for (i=0;i<Nsubgroups;i++){
				fread(&locVelX[i], sizeof(float), 1, fd);
				fread(&locVelY[i], sizeof(float), 1, fd);
				fread(&locVelZ[i], sizeof(float), 1, fd);
			}
			for (i=0;i<Nsubgroups;i++){
				fread(&locCMX[i], sizeof(float), 1, fd);
				fread(&locCMY[i], sizeof(float), 1, fd);
				fread(&locCMZ[i], sizeof(float), 1, fd);
			}
			for (i=0;i<Nsubgroups;i++){
				fread(&locSpinX[i], sizeof(float), 1, fd);
				fread(&locSpinY[i], sizeof(float), 1, fd);
				fread(&locSpinZ[i], sizeof(float), 1, fd);
			}
			for (i=0;i<Nsubgroups;i++) fread(&locVelDisp[i], sizeof(float), 1, fd);
			for (i=0;i<Nsubgroups;i++) fread(&locVmax[i], sizeof(float), 1, fd);
			for (i=0;i<Nsubgroups;i++) fread(&locVmaxRad[i], sizeof(float), 1, fd);
			for (i=0;i<Nsubgroups;i++) fread(&locHalfmassRad[i], sizeof(float), 1, fd);
			for (i=0;i<Nsubgroups;i++) fread(&locIDMostbound[i], sizeof(int), 1, fd);
			for (i=0;i<Nsubgroups;i++) fread(&locGrNr[i], sizeof(int), 1, fd);
			
			for (i=0;i<Nsubgroups;i++){
				fread(&locGasMass[i], sizeof(float), 1, fd);
				fread(&locHaloMass[i], sizeof(float), 1, fd);
				fread(&locDiscMass[i], sizeof(float), 1, fd);
				fread(&locBulgeMass[i], sizeof(float), 1, fd);
				fread(&locStarsMass[i], sizeof(float), 1, fd);
				fread(&locBndryMass[i], sizeof(float), 1, fd);
			}
			
			
			for (i=0;i<Nsubgroups;i++){
				SubhaloLen[i+skip_sub] = locLen[i];
				SubhaloOffset[i+skip_sub] = locOffset[i];
				SubhaloParent[i+skip_sub] = locParent[i];
				SubhaloMass[i+skip_sub] = locMass[i];
				SubhaloPosX[i+skip_sub] = locPosX[i];
				SubhaloPosY[i+skip_sub] = locPosY[i];
				SubhaloPosZ[i+skip_sub] = locPosZ[i];
				SubhaloVelX[i+skip_sub] = locVelX[i];
				SubhaloVelY[i+skip_sub] = locVelY[i];
				SubhaloVelZ[i+skip_sub] = locVelZ[i];
				SubhaloCMX[i+skip_sub] = locCMX[i];
				SubhaloCMY[i+skip_sub] = locCMY[i];
				SubhaloCMZ[i+skip_sub] = locCMZ[i];
				SubhaloSpinX[i+skip_sub] = locSpinX[i];
				SubhaloSpinY[i+skip_sub] = locSpinY[i];
				SubhaloSpinZ[i+skip_sub] = locSpinZ[i];
				SubhaloVelDisp[i+skip_sub] = locVelDisp[i];
				SubhaloVmax[i+skip_sub] = locVmax[i];
				SubhaloVmaxRad[i+skip_sub] = locVmaxRad[i];
				SubhaloHalfmassRad[i+skip_sub] = locHalfmassRad[i];
				SubhaloIDMostbound[i+skip_sub] = locIDMostbound[i];
				SubhaloGrNr[i+skip_sub] = locGrNr[i];
				SubhaloGasMass[i+skip_sub] = locGasMass[i];
				SubhaloHaloMass[i+skip_sub] = locHaloMass[i];
				SubhaloDiscMass[i+skip_sub] = locDiscMass[i];
				SubhaloBulgeMass[i+skip_sub] = locBulgeMass[i];
				SubhaloStarsMass[i+skip_sub] = locStarsMass[i];
				SubhaloBndryMass[i+skip_sub] = locBndryMass[i];
			}
			
			
			skip_sub+=Nsubgroups;
			
			free(locLen);
			free(locOffset);
			free(locParent);
			free(locMass);
			free(locPosX);
			free(locPosY);
			free(locPosZ);
			free(locVelX);
			free(locVelY);
			free(locVelZ);
			free(locCMX);
			free(locCMY);
			free(locCMZ);
			free(locSpinX);
			free(locSpinY);
			free(locSpinZ);
			free(locVelDisp);
			free(locVmax);
			free(locVmaxRad);
			free(locHalfmassRad);
			free(locIDMostbound);
			free(locGrNr);
			free(locGasMass);
			free(locHaloMass);
			free(locDiscMass);
			free(locBulgeMass);
			free(locStarsMass);
			free(locBndryMass);
			
		}
		
		fnr++;
		
	} while (fnr<NTask);	

	return 0;
	
}




int load_group_ids(int snapshot_number)
{

	FILE *fd;
	char   buf[200];
	int    i,j,k,fnr;
	
	int   *locIDs;
	
	fnr = 0;
	skip = 0;
	
	
	do {
		
		sprintf(buf, "groups_%03d/subhalo_ids_%03d.%d",snapshot_number,snapshot_number,fnr);		
		
		if(!(fd=fopen(buf,"r")))
		{
			printf("can't open file `%s`\n",buf);
			exit(0);
		}
		
		printf("reading `%s' ...\n",buf); fflush(stdout);
		
		fread(&Ngroups, sizeof(int), 1, fd);
		fread(&TotNgroups, sizeof(int), 1, fd);
		fread(&Nids, sizeof(int), 1, fd);
		fread(&TotNids, sizeof(long long), 1, fd);
		fread(&NTask, sizeof(int), 1, fd);
		fread(&Offset, sizeof(int), 1, fd);
		
		/*
		printf("%d\n",Ngroups);
		printf("%d\n",TotNgroups);
		printf("%d\n",Nids);
		printf("%lld\n",TotNids);
		printf("%d\n",NTask);
		printf("%d\n",Offset);
		*/
		 
		if (fnr==0) IDs = (int*) calloc (TotNids, sizeof(int));
		
		if (Nids>0){
			
			locIDs = (int*) calloc (Nids, sizeof(int));
			
			for (i=0;i<Nids;i++) fread(&locIDs[i], sizeof(int), 1, fd);
			for (i=0;i<Nids;i++) IDs[i+skip] = locIDs[i];
			
			skip+=Nids;
			
			free(locIDs);
		}
		
		fnr++;
	
	} while (fnr<NTask);	
	
	return 0;
	
}


