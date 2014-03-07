/*

Program to cut out a box of gadget file
*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void usage()
{
  fprintf(stderr,"Program to cut a box out of gadget file, creating a new gadget file\n");
  fprintf(stderr,"-i <input-file> -o <output-file> \n");
  fprintf(stderr,"-p <radius X Y Z> or -b <minx miny minz sizex sizey sizez> \n");
  fprintf(stderr,"-gr <gridsize minx miny minz sizex sizey sizez>\n");
  fprintf(stderr,"-gv <gridx gridy gridz minx miny minz sizex sizey sizez>\n");
  fprintf(stderr,"-rv [reverse: cut everyting but the box]\n");
  fprintf(stderr,"-ascii [output in ascii format]\n");
  fprintf(stderr,"-reset [new positions are relative to box offset]\n\n");
  exit(1);

}

int main (int argc, char *argv[])
{
  typedef float fltarr[3];
  FILE *fp;
  char in1[256],in2[256],out[256];
  int numpart[2],i,j = 0,k,l,blocksize,dummy,conv=-1,cut,a,b,c,x[3], ascii=0;
  double dum[3],boxsize[3],checkmass[3];
  double min[3],max[3],cm[3],rad=0,reset=0,dx=0;
  int grid[3]={0,0,0},gcoord[3],gsize[3],reverse=0;
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
    int bytesleft[22];
  } input0,output;

  fltarr *pos0,*vel0,*posnew,*velnew;
  float *mass0,*massnew;
  int *id,*idnew;

  int passes=1;

  /*
  struct particle 
  {
    int id;
    float pos[3];
    float vel[3];
    float mass;
  } *part1,*part2;
  */
 
#define SKIP fread(&dummy,sizeof(dummy),1,fp);    // to Skip the Block-Size-info at the beginning and end of each data-Block



  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //Edit this values, ignored if given as command line parameters in the same order

  strcpy(in1,"");                 //filename of first gadget-file, 80 characters at most (e.g. top-grid)
  strcpy(out,"");                 //filename of output file
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  cm[2]=0;
  max[2]=0;
  if (argc==1) usage();
  i=1;
  while (i<argc)
    {
      
      if (!strcmp(argv[i],"-h"))
	{
	  usage;
	}
      else if (!strcmp(argv[i],"-o"))
	{
	  i++;
	  strcpy(out,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-i"))
	{
	  i++;
	  strcpy(in1,argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-reset"))
	{
	  i++;
	  reset=1;
	}
      else if (!strcmp(argv[i],"-rv"))
	{
	  i++;
	  reverse=1;
	}
      else if (!strcmp(argv[i],"-ascii"))
	{
	  i++;
	  ascii=1;
	}
      else if (!strcmp(argv[i],"-gr"))
	{
	  i++;
	  grid[0]=atoi(argv[i++]);	  
	  grid[1]=grid[0];
	  grid[2]=grid[0];
	  gcoord[0]=atoi(argv[i++]);
	  gcoord[1]=atoi(argv[i++]);
	  gcoord[2]=atoi(argv[i++]);
	  gsize[0]=atoi(argv[i++]);
	  gsize[1]=atoi(argv[i++]);
	  gsize[2]=atoi(argv[i++]);
	}
      else if (!strcmp(argv[i],"-gv"))
	{
	  i++;
	  grid[0]=atoi(argv[i++]);	  
	  grid[1]=atoi(argv[i++]);	  
	  grid[2]=atoi(argv[i++]);	  
	  gcoord[0]=atoi(argv[i++]);
	  gcoord[1]=atoi(argv[i++]);
	  gcoord[2]=atoi(argv[i++]);
	  gsize[0]=atoi(argv[i++]);
	  gsize[1]=atoi(argv[i++]);
	  gsize[2]=atoi(argv[i++]);
	}
       else if (!strcmp(argv[i],"-b"))
	{
	  i++;
	  
	  min[0]=atof(argv[i++]);
	  min[1]=atof(argv[i++]);
	  min[2]=atof(argv[i++]);
	  max[0]=min[0]+atof(argv[i++]);
	  max[1]=min[1]+atof(argv[i++]);
	  max[2]=min[2]+atof(argv[i++]);	  
	}
      else if (!strcmp(argv[i],"-p"))
	{
	  i++;
	  rad  =atof(argv[i++]);
	  cm[0]=atof(argv[i++]);
	  cm[1]=atof(argv[i++]);
	  cm[2]=atof(argv[i++]);
	} else {
	  usage();
	}

    }

  if ((in1==out) || ((cm[2]==0) && (max[2]==0) && (grid[0]==0) )) usage();
  
  if (rad!=0)
    for (j=0; j<3; j++)
    {
      min[j]=cm[j]-rad;
      max[j]=cm[j]+rad;
    }

  
  //  for (i=0; i<6; i++) {input1.npart[i]=0;input2.npart[i];input1.massarr[i]=0;inp1.nall[i]=0;}
  numpart[0]=0;

  fp=fopen(in1,"r");
  SKIP;
  fread(&input0,sizeof(struct header),1,fp);            //read Header of input file 1
  SKIP;
  for (i=0; i<6; i++) {numpart[0]+=input0.npart[i];}
  pos0=(fltarr *)malloc(sizeof(fltarr)*numpart[0]); //allocate memory for particles in input file 1
  vel0=(fltarr *)malloc(sizeof(fltarr)*numpart[0]); 
  mass0=(float  *)malloc(sizeof(float)*numpart[0]);
  id=(int *)malloc(sizeof(int)*numpart[0]);

  SKIP;                                                 //read particle positions
   fread(&pos0[0],sizeof(fltarr),numpart[0],fp);
  SKIP;
  SKIP;
   fread(&vel0[0],sizeof(fltarr),numpart[0],fp);
  SKIP;
  SKIP;
   fread(&id[0],sizeof(int),numpart[0],fp);
  SKIP;
  SKIP;
   fread(&mass0[0],sizeof(float),numpart[0],fp);
  SKIP;

  if ( input0.npart[0] > 0) 
    {
      passes = 2;
    }

  // printf("%d\n",dummy);  

  fclose(fp);

  if (grid[0]!=0)
    for (j=0; j<3; j++)
    {
      dx=(input0.boxsize/grid[0]);
      min[j]=dx*gcoord[j];
      max[j]=min[j]+gsize[j]*dx;
    }


  posnew=(fltarr *)malloc(sizeof(fltarr)*numpart[0]); //allocate memory for particles in the new file
  velnew=(fltarr *)malloc(sizeof(fltarr)*numpart[0]); 
  massnew=(float  *)malloc(sizeof(float)*numpart[0]);
  idnew=(int *)malloc(sizeof(int)*numpart[0]);

  output=input0;
  for (i=0; i<6; i++) 
    {
      output.npart[i]=0;
      output.nall[i]=0;
    }
  numpart[1]=0;

 for (j=0; j<3; j++)
   {
     printf("min %f max %f\n", min[j],max[j]);
   }
  
 cut=0;
 if (grid[0]==0)
   {
  for (i=0; i< numpart[0]; i++)
    {
      k=1;
      if (!rad)
	{
	  for (j=0; j<3; j++) if ((pos0[i][j] < min[j]) || (pos0[i][j] > max[j])) k=0;
	}
      else
	{
	  double dist=0;
	  for (j=0; j<3; j++) dist += pow((pos0[i][j] - cm[j]),2);
	  dist = sqrt(dist);
	  if (dist > rad) k = 0;
				
	}
      
      
      if (((k) && (!reverse)) || ((!k) && (reverse)))
	{

	  for (j=0; j<3; j++)
	    {
	      if (reset) {posnew[numpart[1]][j] = (pos0[i][j]-min[j]);}
	      else posnew[numpart[1]][j] = pos0[i][j];
	      velnew[numpart[1]][j] = vel0[i][j];	      
	    }
	  idnew[numpart[1]]=id[i];
	  massnew[numpart[1]]=mass0[i];
	  dummy=0;
	  for (l=0; l<5; l++)
	    {
	      dummy+=input0.npart[l];
	      if (i < dummy) break;
	    }
	  output.npart[l]++;
	  output.nall[l]++;
	  numpart[1]++;
	} else cut++;

    }
   }
 else
   {
     if (!reverse)
     {
       for ( k = 0; k < passes; k++ )
	 for (c=gcoord[2]; c<(gsize[2]+gcoord[2]); c++)
	   for (a=gcoord[0]; a<(gsize[0]+gcoord[0]); a++)
	     for (b=gcoord[1]; b<(gsize[1]+gcoord[1]); b++)
	       {
		 x[0]=a;
		 x[1]=b;
		 x[2]=c;
		 for (j=0; j<3; j++)
		   {
		     if (x[j]<0) x[j]=x[j]+grid[j];		
		     if (x[j]>=grid[j]) x[j]=x[j]-grid[j];
		   }
		 int offset = 0;
		 if (k) offset = input0.npart[0];
		 i=x[1]+grid[1]*x[0]+grid[0]*grid[1]*x[2] + offset;		
		 for (j=0; j<3; j++)
		   {
		     if (reset)
		       {
			 posnew[numpart[1]][j] = (pos0[i][j]-min[j]);
			 if (posnew[numpart[1]][j]<(input0.boxsize/(-2))) posnew[numpart[1]][j]+=input0.boxsize;
		       }
		     else posnew[numpart[1]][j] = pos0[i][j];
		     velnew[numpart[1]][j] = vel0[i][j];	      
		   }
		 idnew[numpart[1]]=id[i];
		 massnew[numpart[1]]=mass0[i];
		 dummy=0;
		 for (l=0; l<5; l++)
		   {
		     dummy+=input0.npart[l];
		     if (i < dummy) break;
		   }
		 output.npart[l]++;
		 output.nall[l]++;
		 numpart[1]++;
	      
	}
     } else {
       for ( k = 0; k < passes; k++ )
	 for (c=0; c<grid[2]; c++)
	   for (a=0; a<grid[0]; a++)
	     for (b=0; b<grid[1]; b++)
	       {
		 //	     if ((b==gcoord[1]) && (a>=gcoord[0]) && (a<gcoord[0]+gsize[0]) && (c>=gcoord[2]) && (c<gcoord[2]+gsize[2])) b=b+gsize[1];
		 int skip=1;
		 if ((gcoord[0]+gsize[0])>grid[0])
		   {
		     if (!( ( a < (gcoord[0]+gsize[0])%grid[0]) || (a >= gcoord[0]))) skip=0;
		   }
		 else
		   {
		     if (!( ( a < (gcoord[0]+gsize[0])) && (a >= gcoord[0]))) skip=0;
		   }
		 if ((gcoord[1]+gsize[1])>grid[1])
		   {
		     if (!( ( b < (gcoord[1]+gsize[1])%grid[1]) || (b >= gcoord[1]))) skip=0;
		   }
		 else
		   {
		     if (!( ( b < (gcoord[1]+gsize[1])) && (b >= gcoord[1]))) skip=0;
		   }
		 if ((gcoord[2]+gsize[2])>grid[2])
		   {
		     if (!( ( c < (gcoord[2]+gsize[2])%grid[2]) || (c >= gcoord[2]))) skip=0;
		   }
		 else
		   {
		     if (!( ( c < (gcoord[2]+gsize[2])) && (c >= gcoord[2]))) skip=0;
		   }
		 //	     if (b>=grid[1]) continue;
		 if (skip) continue;
	      
		 int offset = 0;
		 if (k) offset = input0.npart[0];
		 
		 i=b+grid[1]*a+grid[0]*grid[1]*c;
		 for (j=0; j<3; j++)
		   {
		     if (reset) 
		       {
			 posnew[numpart[1]][j] = (pos0[i][j]-min[j]);
			 if (posnew[numpart[1]][j]<(input0.boxsize/(-2))) posnew[numpart[1]][j]+=input0.boxsize;
		       }
		     else posnew[numpart[1]][j] = pos0[i][j];
		     velnew[numpart[1]][j] = vel0[i][j];	      
		   }
		 idnew[numpart[1]]=id[i];
		 massnew[numpart[1]]=mass0[i];
		 dummy=0;
		 for (l=0; l<5; l++)
		   {
		     dummy+=input0.npart[l];
		     if (i < dummy) break;
		   }
		 output.npart[l]++;
		 output.nall[l]++;
		 numpart[1]++;	       
	}

     }
   }
  printf("new file: %d particles\n", numpart[1]);

  fp=fopen(out,"w");
  if (!ascii)
    {
  blocksize=256;
  fwrite(&blocksize,sizeof(int),1,fp);          //write Header;
  fwrite(&output,sizeof(struct header),1,fp);
  fwrite(&blocksize,sizeof(int),1,fp);

  blocksize=numpart[1]*4*3;                                   //write positions
  fwrite(&blocksize,sizeof(int),1,fp);
  fwrite(&posnew[0],sizeof(fltarr),numpart[1],fp);
  fwrite(&blocksize,sizeof(int),1,fp); 

  blocksize=numpart[1]*4*3;                                  //write velocities
  fwrite(&blocksize,sizeof(int),1,fp);
  fwrite(&velnew[0],sizeof(fltarr),numpart[1],fp);
  fwrite(&blocksize,sizeof(int),1,fp); 



  blocksize=numpart[1]*4;                                    //write particle IDs
  fwrite(&blocksize,sizeof(int),1,fp);
  fwrite(&idnew[0],sizeof(int),numpart[1],fp);
  fwrite(&blocksize,sizeof(int),1,fp); 

  blocksize=numpart[1]*4;                                    //write masses
  fwrite(&blocksize,sizeof(int),1,fp);
  fwrite(&massnew[0],sizeof(float),numpart[1],fp);
  fwrite(&blocksize,sizeof(int),1,fp); 
    }
  else 
    {
      for (i=0; i < 6; i++) fprintf(fp,"%d ", output.npart[i]);
      fprintf(fp,"\n");
      for (i=0; i<numpart[1]; i++)
	{
	  for (j=0; j <3; j++) fprintf(fp,"%12g ",posnew[i][j]);
	  for (j=0; j <3; j++) fprintf(fp,"%12g ",velnew[i][j]);
	  fprintf(fp, " %6d ", idnew[i]);
	  fprintf(fp, " %6g ", massnew[i]);
	  fprintf(fp,"\n");
	}
    }


  fclose(fp);
  
  return 0;
}
