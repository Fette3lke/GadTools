/*
Program to check (and alter) the header of gadgetfiles
gcc -lm
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

int main  (int argc, char *argv[])
{
  typedef float fltarr[3];
  struct header
  {
    int npart[6];
    double massarr[6];
    double time;
    double redsh;
    int flg_sfr;
    int flg_fdbck;
    int nall[6];
    int flag_cool;
    int numfiles;
    double boxsize;
    double omega0;
    double omegal;
    double hubparam;
    int flg_age;
    int flg_mtl;
    //    float shift[3];
    int bytesleft[22];
  } head;
  FILE *fp, *out;
  char gadgetfile[180], outfile[180];
  unsigned int i=0,blocksize,dum, *id, blocks;
  unsigned int numpart=0, notgas=0;
  long long size=0;
  float *mass_ev, fdum, *mass;
  fltarr *pos0,*vel0;
  double mdens,masstot;
  int j=0,k=0,domass=0,verbose=0,check=0, test=0, unsplit=0, massarray=0, addgas=0, kill_massarray=0;
  int forcewrite=0;
 
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //filenames are ignored if given as command line parameters in the same order

  strcpy(gadgetfile,"<nofile>");
  strcpy(outfile,"<nofile>");
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  i=1;
  while (i<argc)
    {
      if (!strcmp(argv[i],"-o"))
	{
	  i++;
	  strcpy(outfile,argv[i]);
	} else
	  {
      if (!strcmp(argv[i],"-m")) domass=1;
      if (!strcmp(argv[i],"-v")) verbose=1;
      if (!strcmp(argv[i],"-c")) check=1;
      if (!strcmp(argv[i],"-t")) test=1;
      if (!strcmp(argv[i],"-u")) unsplit=1;
      if (!strcmp(argv[i],"-a")) addgas=1;
      if (!strcmp(argv[i],"-ma")) massarray=1;
      if (!strcmp(argv[i],"-x")) forcewrite=1;
      if (!strcmp(argv[i],"-kma")) {massarray=1;kill_massarray=1;}

      if ((!strcmp(argv[i],"-h")) || (!strcmp(argv[i],"--help")))  {i=1; break;}
      if (*argv[i]!='-') strcpy(gadgetfile,argv[i]);
	  }
      i++;
      
    }
  if ((i==1) || (!strcmp(gadgetfile,"<nofile>")))
    {
      
      fprintf(stderr,"\nphead <options> <gadgetfile>\n");
      fprintf(stderr,"[-c check blocksizes] [-u unsplit gadget file (for files created with grafic2gadget_split)] \n");
      fprintf(stderr,"[-ma create massarray] [-a add 0s for internal energy and density] [-o <outfile>]\n");
      fprintf(stderr,"[-kma remove massarray]\n");
      fprintf(stderr,"[-v verbose] [-m calculate total mass etc.] \n\n");
      exit(1);
    }
  //if (argc>1) strcpy(gadgetfile,argv[1]);


 

  if ((fp=fopen(gadgetfile,"r"))==NULL) 
    {
      fprintf(stderr,"\ncan't read %s\n\n",gadgetfile);
      exit(1);
    }

  if ((check) || (unsplit) || (addgas))
    {
      i=0;
      //      for (i=1; i<=5; i++)
      //{
	  while(fread(&blocksize,sizeof(int),1,fp))
	    {
	      if (check) printf("Blocksize %d: %ld\n",i ,blocksize);
	  if ((test) && (i==2))
	    {
	      dum=0;
	      while (fread(&fdum,sizeof(float),1,fp))
		{
	       if ((fdum<0)||(fdum>72000))
		 {
		   printf("position %u hex %8.X uint %10.u float %f\n", ftell(fp) ,fdum ,fdum, fdum);
		   //		   printf("position %u\n", ftell(fp));
		 }
	       
		}

	       //	       for (j=31; j>=0; j--) printf("%d", ((1<<j)&blocksize)>>j);
	       //	       printf("\n%u\n", ftell(fp));
	       exit(0);
	    }
	  fseek(fp,(blocksize),SEEK_CUR);
	  //printf("%u\n", ftell(fp));
	  fread(&blocksize,sizeof(int),1,fp);
	  if (check) printf("Blocksize %d: %ld\n",i ,blocksize);
	  size+=blocksize+8;
	  i++;
	  if (blocksize==0) break;
	    }
	  //}
      rewind(fp);
      printf("Size : %ld\n",size);
    }
  blocks=i;
  if (blocks!=7) unsplit=0;

  fread(&blocksize,sizeof(int),1,fp);
  if (blocksize==256) printf("Headersize OK\n"); 
  else printf("Bad Headersize %d\n",blocksize);
  fread(&head,sizeof(struct header),1,fp);                        //read header of evolved gadget file
//  printf("struct %d\n",sizeof(struct header));
  fread(&blocksize,sizeof(int),1,fp);


  for (i=0; i<6; i++) numpart+=head.npart[i];
  notgas=numpart-head.npart[0];
  mass_ev=(float *) malloc(sizeof(float)*numpart);

  if ((!unsplit) && (!massarray) && (!forcewrite))
    {
  dum=0;
  for (i=0; i<6; i++) if ((head.massarr[i]==0) && (head.npart[i]!=0)) dum=1;

  if (dum)
    {
  fread(&blocksize,sizeof(int),1,fp);                                //skip particle positions
  fseek(fp,numpart*12,SEEK_CUR);
  fread(&blocksize,sizeof(int),1,fp);  

  
  fread(&blocksize,sizeof(int),1,fp);                                //skip velocities
  fseek(fp,numpart*12,SEEK_CUR);
  fread(&blocksize,sizeof(int),1,fp);

  fread(&blocksize,sizeof(int),1,fp);                                //skip particle ids
  fseek(fp,numpart*4,SEEK_CUR);
  fread(&blocksize,sizeof(int),1,fp);

  
  fread(&blocksize,sizeof(int),1,fp);                                //read particle masses of evolved snapshot file
  if (check) printf("Blocksize Masses: %ld\n", blocksize);
  if (domass) fread(&mass_ev[0],sizeof(float),numpart,fp);
  else {fread(&mass_ev[0],sizeof(float),1,fp);   fseek(fp,(numpart-1)*4,SEEK_CUR);}
  

  fread(&blocksize,sizeof(int),1,fp);
    }

    } else if (outfile!="<nofile>") {
      out=fopen(outfile,"w");

      pos0=(fltarr *)malloc(sizeof(fltarr)*numpart); //allocate memory for particles
      vel0=(fltarr *)malloc(sizeof(fltarr)*numpart);
      
      fread(&blocksize,sizeof(int),1,fp);
      fread(&pos0[0],sizeof(fltarr),head.npart[0],fp);
      if (unsplit)
	{
      fread(&blocksize,sizeof(int),1,fp);  
      fread(&blocksize,sizeof(int),1,fp);
	}
      fread(&pos0[head.npart[0]],sizeof(fltarr),notgas,fp);
      fread(&blocksize,sizeof(int),1,fp);  
  
      fread(&blocksize,sizeof(int),1,fp);
      fread(&vel0[0],sizeof(fltarr),head.npart[0],fp);
      if (unsplit)
	{
      fread(&blocksize,sizeof(int),1,fp);
      fread(&blocksize,sizeof(int),1,fp);
	}
      fread(&vel0[head.npart[0]],sizeof(fltarr),notgas,fp);
      fread(&blocksize,sizeof(int),1,fp);  

      if ((massarray) || (forcewrite))
	  {
	    mass=(float *)malloc(sizeof(float)*numpart);
	    id=(int *)malloc(sizeof(int)*numpart);	  
	    fread(&blocksize,sizeof(int),1,fp);
	    fread(&id[0],sizeof(int),numpart,fp);
	    fread(&blocksize,sizeof(int),1,fp);
	    fread(&blocksize,sizeof(int),1,fp);
	    fread(&mass[0],sizeof(float),numpart,fp);
	    fread(&blocksize,sizeof(int),1,fp);
	    dum=0;
	    for (i=0; i<6; i++)
	      {
		if ((head.npart[i]>0) && (massarray)) head.massarr[i]=mass[dum];
		if (kill_massarray) head.massarr[i]=0;
		dum+=head.npart[i];
	      }
	    
	  } 
      
      blocksize=256;      
      fwrite(&blocksize,sizeof(int),1,out);          //write Header;
      fwrite(&head,sizeof(struct header),1,out);
      fwrite(&blocksize,sizeof(int),1,out);

      dum=numpart*12;
      fwrite(&dum,sizeof(int),1,out);
      fwrite(&pos0[0],sizeof(fltarr),numpart,out);
      fwrite(&dum,sizeof(int),1,out);

      dum=numpart*12;
      fwrite(&dum,sizeof(int),1,out);
      fwrite(&vel0[0],sizeof(fltarr),numpart,out);
      fwrite(&dum,sizeof(int),1,out);
      
      dum=numpart*4;
      fwrite(&dum,sizeof(int),1,out);
      fwrite(&id[0],sizeof(int),numpart,out);
      fwrite(&dum,sizeof(int),1,out);
      
      if ((!massarray) || (kill_massarray))
	{
	  blocksize=numpart*4;
	  fwrite(&blocksize,sizeof(int),1,out);
	  fwrite(&mass[0],sizeof(float),numpart,out);
	  fwrite(&blocksize,sizeof(int),1,out);
	  
	}
      if (addgas)
	{
      blocksize=head.npart[0]*4;
      if (blocksize>0)
	{
      dum=0;
       fwrite(&blocksize,sizeof(int),1,out); 
       for (i=0; i<head.npart[0]; i++) fwrite(&dum,sizeof(int),1,out);
       fwrite(&blocksize,sizeof(int),1,out); 

       fwrite(&blocksize,sizeof(int),1,out); 
       for (i=0; i<head.npart[0]; i++) fwrite(&dum,sizeof(int),1,out);
       fwrite(&blocksize,sizeof(int),1,out); 
	}
	} else {
      while (fread(&dum,sizeof(int),1,fp))
	{
	  fwrite(&dum,sizeof(int),1,out);
	}
	}
      
      fclose(out);
    }

  
  fclose(fp);
  if ((addgas) && (!unsplit) && (!massarray) && (outfile!="<nofile>") && (blocks<6)  )
    {
      fp=fopen(gadgetfile,"a");
      blocksize=head.npart[0]*4;
      if (blocksize>0)
	{
      dum=0;
       fwrite(&blocksize,sizeof(int),1,fp); 
       for (i=0; i<head.npart[0]; i++) fwrite(&dum,sizeof(int),1,fp);
       fwrite(&blocksize,sizeof(int),1,fp); 

       fwrite(&blocksize,sizeof(int),1,fp); 
       for (i=0; i<head.npart[0]; i++) fwrite(&dum,sizeof(int),1,fp);
       fwrite(&blocksize,sizeof(int),1,fp); 
	}
      fclose(fp);
    }
  
  masstot=0;
  if (domass)
   {
  for (i=0; i< numpart; i++) 
   {
       masstot+=mass_ev[i];
   } 
  mdens=masstot/pow(head.boxsize*head.time,3);
   }
 
  printf("%20s %12d\n","numpart",numpart);
  printf("%20s %1.10f\n","a",head.time);
  printf("%20s %10.2f\n","redshift",head.redsh);
  printf("%20s %10.2f\n","boxsize",head.boxsize);
  printf("%20s %10.2f\n","Omega_L",head.omegal);
  printf("%20s %10.2f\n","Omega_0",head.omega0);
  printf("%20s %10.2f\n","Hubparam",head.hubparam);
  printf("%20s %1.10g\n","mass of particle 0",mass_ev[0]);
  if (verbose) {
  printf("%20s %10d\n","Numfiles",head.numfiles);  
  printf("%20s %10d\n","Flag SFR",head.flg_sfr);
  printf("%20s %10d\n","Flag Feedback",head.flg_fdbck);
  printf("%20s %10d\n","Flag Cooling",head.flag_cool);
  printf("%20s %10d\n","Flag Age",head.flg_age);
  printf("%20s %10d\n","Flag Metals",head.flg_mtl);
  for (i=0; i<3; i++)
    //printf("%18s %d %10.2f\n","Shift",i,head.shift[i]); 
  printf("%20s %10d\n","unused",sizeof(head.bytesleft));
  for (i=0; i<6; i++)
  printf("%18s %d %10d\n","Npart",i,head.npart[i]);
  for (i=0; i<6; i++)
  printf("%18s %d %1.10f\n","MassArr",i,head.massarr[i]);
  for (i=0; i<6; i++)
  printf("%18s %d %10d\n","Nall",i,head.nall[i]);
  
  }
  if (domass) {
  printf("%20s %10.2f\n","total mass",masstot);
  printf("%20s %1.15e\n","mean density",mdens);
  }

  return 0;
}
