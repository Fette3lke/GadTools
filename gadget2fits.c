/*
Program to produce fits files from gadget snapshots

gcc -fopenmp -lgad-stan -lgsl -lgslcblas gadget2fits.c -o ~/bin/gadget2fits -I$HOME/usr/include -L$HOME/usr/lib -lcfitsio
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "libgad.h"
#include "fitsio.h"

#define USE 16

void usage()
{
  	  fprintf(stderr," Create FITS file for gas surface density v0.01\n");
  	  fprintf(stderr," -i   <input file name>\n");
  	  fprintf(stderr," -t   <temperature threshold>\n");
  	  fprintf(stderr," -use <bitcode of particles to be used for inertia tensor>\n");
  	  fprintf(stderr," -r   <max-distance of particles to be considered for inertia tensor>\n");
  	  fprintf(stderr," -vb  <number of velocity bins>\n");
  	  fprintf(stderr," -vm  <sets range of velocity bins [-vm, vm]>\n");
  	  fprintf(stderr," -rf  <rotate file> (sets viewport)\n");
  	  fprintf(stderr," -srf <save rotate file>\n");
  	  fprintf(stderr," -s   <smoothing length for stellar particles>\n\n");
  	  fprintf(stderr," choose two:\n");
  	  fprintf(stderr,"    -b <boxsize>\n");
  	  fprintf(stderr,"    -g <gridsize>\n");
  	  fprintf(stderr,"    -p <pixelsize>\n\n");
  	  fprintf(stderr,"    -nophys   <do not convert pixelsize to physical units (instead use code units)>\n\n");
  	  fprintf(stderr," \n\n");
	  exit(1);
}

#ifndef _OPENMP
int omp_get_thread_num() {return 0;}
//int omp_get_num_threads() {return 1;}
#endif

int main (int argc, char *argv[])
{
  FILE *fp;
  char infile[256];
  char fitsfilename[256];
  char rotfile[256];
  int i,j,k,n, usepart;
  int gridsize = 0;
  int verbose = 0;
  long ii;
  int velbins = 64;
  double velmax = 320;
  double boxsize = 0;
  double binsize = 0;
  double rotate_dist = 0;
  double center[3] = {0., 0., 0.};
  double stellar_hsml = 0.4;
  struct gadpart *part, *wpart;
  struct header head;
  double tempthreshold=1.e5;
  int convert_phys = 1;
  int save_rotation = 0;
  int load_rotation = 0;

  double *dens;
  double **dens_tmp;

  double *sdens;
  double **sdens_tmp;

  int nproj = 7;
  int iproj;
  double proj_angle;
 
  int numthreads ;
#pragma omp parallel
  {
    int thread= omp_get_thread_num();
    if (thread == 0)
      numthreads = omp_get_num_threads();
  }

  dens_tmp = (double**) malloc(numthreads * sizeof(double *)); 
  sdens_tmp = (double**) malloc(numthreads * sizeof(double *)); 
  
  strcpy(fitsfilename,"default.fits");

  i=1;
  usepart=USE;
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
      else if (!strcmp(argv[i],"-o"))
	{
	  i++;
	  strcpy(fitsfilename,argv[i]);
	  i++;
	} 
      else if (!strcmp(argv[i],"-rf"))
	{
	  i++;
	  strcpy(rotfile,argv[i]);
	  load_rotation=1;
	  i++;
	} 
      else if (!strcmp(argv[i],"-srf"))
	{
	  i++;
	  strcpy(rotfile,argv[i]);
	  save_rotation=1;
	  i++;
	} 
      else if (!strcmp(argv[i],"-use")) {
	i++;
	if (!strcmp(argv[i],"all")) usepart=63;
	else usepart=atoi(argv[i]);
	i++;
      }
      else if (!strcmp(argv[i],"-g")) 
	{
	  i++;
	  gridsize=atoi(argv[i]);
	  i++;
	} 
      else if (!strcmp(argv[i],"-np")) 
	{
	  i++;
	  nproj=atoi(argv[i]);
	  i++;
	} 
      else if (!strcmp(argv[i],"-vb")) 
	{
	  i++;
	  velbins=atoi(argv[i]);
	  i++;
	} 
      else if (!strcmp(argv[i],"-vm")) 
	{
	  i++;
	  velmax=atof(argv[i]);
	  i++;
	} 
      else if (!strcmp(argv[i],"-b")) 
	{
	  i++;
	  boxsize=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-s")) 
	{
	  i++;
	  stellar_hsml=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-r")) 
	{
	  i++;
	  rotate_dist=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-p")) 
	{
	  i++;
	  binsize=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-nophys")) 
	{
	  i++;
	  convert_phys = 0;
	}
      else if (!strcmp(argv[i],"-t")) 
	{
	  i++;
	  tempthreshold=atof(argv[i]);
	  i++;
	}
      else if (!strcmp(argv[i],"-v")) 
	{
	  i++;
	  verbose=1;
	}
      else if (!strcmp(argv[i],"-c")) 
	{
	  i++;
	  center[0]=atof(argv[i]);
	  i++;
	  center[1]=atof(argv[i]);
	  i++;
	  center[2]=atof(argv[i]);
	  i++;
	} else {
	usage();
      }
    }

  unsigned int numpart_all;

  if (!(numpart_all=readgadget_part(infile, &head, &part))) 
    {
      extern int libgaderr;
      printf("error reading file %s\nError Code %d\n",infile, libgaderr);
      exit(1);
    }
  

  if (nproj > 1)
    {
      proj_angle = 90 / ((nproj-1) /2. );
    } else nproj = 1;

  if (convert_phys)
    {
      binsize = binsize * head.hubparam / head.time;
    }

  if ((binsize != 0) && (gridsize !=0) && (boxsize==0)) 
    {
      boxsize = gridsize * binsize;
    }
  else if ((binsize != 0) && (gridsize ==0) && (boxsize!=0)) 
    {
      gridsize = ceil(boxsize/binsize);
    }
  else if ((binsize == 0) && (gridsize !=0) && (boxsize!=0)) 
    {
      binsize = boxsize / gridsize;
    }
  else usage();

  if  (!strcmp(fitsfilename,"default.fits"))
    sprintf(fitsfilename,"!%s.fits",infile);
  else
    {
      char chdum[256];
      strcpy(chdum, fitsfilename);
      sprintf(fitsfilename,"!%s",chdum);
    }

  if (verbose) printf("building grid...\n");
  long g3 = gridsize * gridsize * velbins;
  long g2 = gridsize * gridsize;
  long g  = gridsize;

  dens = (double *) calloc (g3, sizeof( double ));
  sdens = (double *) calloc (g3, sizeof( double ));
  if (verbose) printf("allocating grid for %d processes...\n", numthreads);
#pragma omp parallel
    {
      int thread= omp_get_thread_num();
      dens_tmp[thread] = (double *) calloc (g3, sizeof(double));
      sdens_tmp[thread] = (double *) calloc (g3, sizeof(double));
    }
  if (verbose) printf("read snapshot...\n");
 
  if ( ( center[0] != 0 ) || ( center[1] != 0 ) || ( center[2] != 0 ) )
    {
      for ( n = 0; n < numpart_all; n++)
	{
	  for ( j = 0; j < 3; j++)
	    part[n].pos[j] -= center[j];
	}
    }
  if (verbose) printf("rotate system...\n");
  double ratios[2] = {0,0};
  gsl_matrix *rotation;
      //      rotategalaxy(wpart, num, -1, use_cm);//, &ratios[0], num);
  if (load_rotation)
    {
      FILE *matrixf=fopen(rotfile,"r");
      rotation = gsl_matrix_alloc(3,3);
      gsl_matrix_fread (matrixf, rotation);
      fclose(matrixf);
      rotatepart(part, numpart_all, rotation);
    }
  else if (usepart)
    {
      if (!rotate_dist) rotate_dist = boxsize/2.;
      rotategalaxy(part, numpart_all, rotate_dist, usepart, &ratios[0], &rotation);

      if (save_rotation)
	{
	  FILE *matrixf=fopen(rotfile,"w");
	  gsl_matrix_fwrite (matrixf, rotation);
	  fclose(matrixf);
	}
    }
  gsl_matrix_free(rotation);
  


  /*********************************************************************

      Program code goes here

  *********************************************************************/

  int d;
  const double boxhalf = boxsize / 2.;
  const double cellhalf = boxhalf / gridsize;
  const double cellvol = 8.0 * cellhalf * cellhalf * cellhalf;
  const double cellarea = 4.0 * cellhalf * cellhalf;
  const double binsizeinv = gridsize / boxsize;
  const starind = head.npart[0] + head.npart[1] + head.npart[2] + head.npart[3];
  double velbinsize = 2*velmax / velbins;
  const double convert = 10 * SQR(head.hubparam) * binsize * 1000 / SQR(head.time); //convert from 10^10Msun*h^2/kpc^3 comoving to Msun/pc^2 physical
  double total_proj_angle = 0;
  fitsfile *fptr;
  int status = 0;
  fits_create_file(&fptr, fitsfilename, &status);

  if (verbose) printf("binnning %d gas particles\n", head.npart[0], gridsize);
  if (verbose) printf("start binning...\n");

  for (iproj = 0; iproj < nproj; iproj++)
    {

      for ( i=0; i<g3; i++)
	{
	  dens[i] = 0;
	  sdens[i] = 0;
	}

#pragma omp parallel private(i)
      {
	int thread= omp_get_thread_num();
	for ( i=0; i<g3; i++)
	  {
	    dens_tmp[thread][i] = 0;
	    sdens_tmp[thread][i] = 0;
	  }
      }

      if (iproj)
	{
	  if (iproj == ((nproj-1)/2)+1)
	    {
	      xrotate(-total_proj_angle, part, numpart_all);
	      total_proj_angle = 0;
	    }
	  if (iproj <= (nproj-1)/2)
	    xrotate(proj_angle, part, numpart_all);
	  else
	    yrotate(proj_angle, part, numpart_all);
	  total_proj_angle += proj_angle;
	}
      if (verbose) printf("Projection Angle: %g\n", total_proj_angle);
//      char testname[256];
//      sprintf(testname,"test%d.gad", iproj);
//      writegadget_part(testname, head, part);
      
#pragma omp parallel for private(i, j, k)
      for ( n = 0; n < (head.npart[0]+head.npart[4]); n++ )
	{
	  long index;
	  double h;
	  int star = 0;
	  if (n < head.npart[0]) 
	    {
	      index = n;
	      if (temperature(part[index]) > tempthreshold) continue;
	      h = part[index].sph->hsml;
	    }
	  else
	    {
	      index = n + head.npart[1] + head.npart[2] + head.npart[3];
	      h = stellar_hsml;
	      star = 1;
	    }
	  double x = part[index].pos[0];
	  double y = part[index].pos[1];
	  double z = part[index].pos[2];
	  double vz =part[index].vel[2];
	  double pmass = part[index].mass;


	  if ( ((x+2*h) < -boxhalf) || ((x-2*h) > boxhalf) ) continue;
	  if ( ((y+2*h) < -boxhalf) || ((y-2*h) > boxhalf) ) continue;
	  if ( ((z+2*h) < -boxhalf) || ((z-2*h) > boxhalf) ) continue;
     
	  int ix = floor( (x + boxhalf) * binsizeinv );
	  int iy = floor( (y + boxhalf) * binsizeinv );
	  int iz = floor( (z + boxhalf) * binsizeinv );
	  int inc = 2.0 * h * binsizeinv;
      
	  double gx = ( ix / binsizeinv ) - boxhalf + cellhalf;
	  double gy = ( iy / binsizeinv ) - boxhalf + cellhalf;
	  double gz = ( iz / binsizeinv ) - boxhalf + cellhalf;

	  double dist2 = ( SQR( x - gx ) + SQR( y - gy ) + SQR( z - gz ) ) / SQR( h );
	  int thread = omp_get_thread_num();
	  if ( dist2 < 4.0 )
	    {
	      for ( i = (ix - inc); i <= (ix+inc); i++ )
		{
		  if ((i<0) || (i>= gridsize)) continue;
		  for ( j = (iy - inc); j <= (iy+inc); j++ )
		    {
		      if ((j<0) || (j>= gridsize)) continue;
		      for ( k = (iz - inc); k <= (iz+inc); k++ )
			{
			  if ((k<0) || (k>= gridsize)) continue;
			  gx = ( i / binsizeinv ) - boxhalf + cellhalf;
			  gy = ( j / binsizeinv ) - boxhalf + cellhalf;
			  gz = ( k / binsizeinv ) - boxhalf + cellhalf;
			  double val = 0.;
			  dist2 = ( SQR( x - gx ) + SQR( y - gy ) + SQR( z - gz ) ) / SQR( h );
			  if ( dist2 < 4.0 )
			    {
			      double weight = 0;
			      double dist = sqrt( dist2 );
			      if ( dist < 1.0 )
				{
				  weight = 1.0 -1.5 * dist2 + 0.75 * dist2 * dist;
				}
			      else
				{
				  double dif2 = 2.0 - dist;
				  weight = 0.25 * dif2 * dif2 * dif2;
				}
			      val = weight / (h*h*h) / M_PI * pmass;
			    }
			  if (val > 0)
			    {
			      long vind = floor((vz + velmax) / velbinsize);
//			      if (vind < 0) vind = 0;
//			      else if (vind >= velbins) vind = velbins - 1;
			      if (vind < 0) continue;
			      else if (vind >= velbins) continue;

			      long ind = i + g * j + g2 * vind;
			      {			
				if (star)
				  sdens_tmp[thread][ind] += val * convert;
				else
				  dens_tmp[thread][ind] += val * convert;
			      }
			    }
			}
		    }
		}

	    }
	  else
	    {
	      if ( (ix<0) || (ix >= gridsize) ) continue;
	      if ( (iy<0) || (iy >= gridsize) ) continue;
	      if ( (iz<0) || (iz >= gridsize) ) continue;

	      long vind = floor((vz + velmax) / velbinsize);
//	      if (vind < 0) vind = 0;
//	      else if (vind >= velbins) vind = velbins - 1;
	      if (vind < 0) continue;
	      else if (vind >= velbins) continue;

	      long ind = i + g * j + g2 * vind;
	      {
		if (star)
		  sdens_tmp[thread][ind] += pmass/cellvol * convert;
		else
		  dens_tmp[thread][ind] += pmass/cellvol * convert;
	      }
	    }
	}


      long fpixel = 1;
      long lpixel = g3;
      int naxis = 3;
      long npixels[3] = {gridsize, gridsize, velbins};
      double pixelsize = binsize / head.hubparam * head.time;
      long nelements = g3;
      double dum = 0;
      double aexpn = head.time;
      for ( ii = 0; ii < nelements; ii++)
	{
	  for ( i = 0; i < numthreads; i++ ) 
	    {
	      dens[ii] += dens_tmp[i][ii];
	      sdens[ii] += sdens_tmp[i][ii];
	    }
	}

      fits_create_img(fptr, DOUBLE_IMG, naxis, npixels, &status );
      fits_update_key(fptr, TDOUBLE, "TYPE", &dum , "Cold Gas Distribution", &status);
      fits_update_key(fptr, TDOUBLE, "PIXELSIZE", &pixelsize , "Size of Pixel in kpc", &status);
      fits_update_key(fptr, TDOUBLE, "VELBINSIZE", &velbinsize , "Size of Velocity bins in km/s", &status);
      fits_update_key(fptr, TDOUBLE, "VELBINMAX", &velmax , "Range of Velocities in km/s, [-velbinmax, velbinmax]", &status);
      fits_update_key(fptr, TDOUBLE, "PROJECTION ANGLE", &total_proj_angle, "Projection Angle in degrees (0 = face-on)", &status);
      fits_update_key(fptr, TDOUBLE, "AEXP", &aexpn, "Expansion factor a of the snapshot", &status);
      //  fits_write_img(fptr, TDOUBLE, fpixel, nelements, dens, &status);
      fits_write_3d_dbl(fptr, 1, g, g, g, g, velbins, dens, &status);

      fits_create_img(fptr, DOUBLE_IMG, naxis, npixels, &status );
      fits_update_key(fptr, TDOUBLE, "TYPE", &dum , "Stellar Distribution", &status);
      fits_update_key(fptr, TDOUBLE, "PIXELSIZE", &pixelsize , "Size of Pixel in kpc", &status);
      fits_update_key(fptr, TDOUBLE, "VELBINSIZE", &velbinsize , "Size of Velocity bins in km/s", &status);
      fits_update_key(fptr, TDOUBLE, "VELBINMAX", &velmax , "Range of Velocities in km/s, [-velbinmax, velbinmax]", &status);
      fits_update_key(fptr, TDOUBLE, "PROJECTION ANGLE", &total_proj_angle, "Projection Angle in degrees (0 = face-on)", &status);
      fits_update_key(fptr, TDOUBLE, "AEXP", &aexpn, "Expansion factor a of the snapshot", &status);
      //  fits_write_img(fptr, TDOUBLE, fpixel, nelements, dens, &status);
      fits_write_3d_dbl(fptr, 1, g, g, g, g, velbins, sdens, &status);

      //  fits_write_subset_dbl(fptr, 1, 3, npixels, &fpixel, &lpixel, dens, &status);
    }
  fits_close_file(fptr, &status);
  fits_report_error(stderr, status);
  //  free(surfdens);
  free(dens);
  return (status);
}
