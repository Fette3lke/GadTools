#ifndef LIBGAD
#define LIBGAD
#ifndef NOPOT
#define POTENTIAL
#endif
#include <stdio.h>
#include <math.h>

#ifndef NOGSL
#ifndef GSL
#define GSL
#endif
#endif

#ifdef GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#endif

#define	MIN(a, b)   ((a)<(b)?(a):(b))
#define	MAX(a, b)   ((a)>(b)?(a):(b))
#define ABS(a) ((a) >= 0 ? (a) : -(a))
#define SQR(x) ((x)*(x))
#define OMEGA_M 0.26
#define OMEGA_L 0.74
#define HUB 0.72
#define GRAV 6.6742e-11
#define MSUN 1.989e30
#define KPC 3.085678e19  //meters
#define sec_per_yr 3.155e7
#define GALAGE(a) galage(a, OMEGA_M, OMEGA_L, HUB)
#define TIMEDIFF(a,b) timediff(a, b, OMEGA_M, OMEGA_L, HUB)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef LONGIDS
typedef unsigned long long IDtype;
#else 
typedef unsigned int IDtype;
#endif


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
    int bytesleft[22];
};

typedef struct sphdata
{
  float u;
  float rho;
  float nelec;
  float nh;
  float hsml;
  float sfr;
#ifdef WINDS
  float dtime;
  float metals[4];
  float tmax;
  float n_spawn;
#endif //WINDS
#ifdef METALS
  float temp;
#endif
} sphdata;

#if defined(WINDS) || defined(METALS)
typedef struct stardata
{
#ifdef WINDS
  float metals[4];
  float tmax;
  float n_spawn;
#else
  int let;
  float initialmass;
#endif
} stardata;
#endif //WINDS || METALS

typedef struct gadpart
{
  fltarr pos;
#ifndef NOVEL
  fltarr vel;
#endif //NOVEL
  float mass;
  IDtype id;
  short type;
#ifndef NOGAS
  sphdata *sph;
  float stellarage;
#endif
#ifdef POTENTIAL
  float pot;
#endif

#if defined(WINDS) || defined(METALS)
  stardata *sd;
#endif //WINDS

#ifdef METALS
  float* metals;
#endif
} gadpart;

typedef struct gadpart_dist
{
  gadpart part;
  double dist;
} gadpart_dist;

int cmp_id (const void *first, const void *second);
int cmp_pointer_id(const void *a, const void *b);
int cmp_type (const void *first, const void *second);
int cmp_pos (const void *first, const void *second);
int cmp_x (const void *first, const void *second);
int cmp_y (const void *first, const void *second);
int cmp_z (const void *first, const void *second);
int cmp_dist (const void *first, const void *second);
int cmp_int (const void *first, const void *second);
int cmp_float (const void *first, const void *second);
unsigned int readgadget(char *filename, struct header *h, fltarr **p, fltarr **v, int **n, float **m);
unsigned int readgadget_part(char *filename, struct header *h,struct gadpart **part);
#ifdef LONGIDS
unsigned int writegadget(char *filename, struct header h, fltarr *p, fltarr *v, long *n, float *m);
#else
unsigned int writegadget(char *filename, struct header h, fltarr *p, fltarr *v, int *n, float *m);
#endif
unsigned int writegadget_part(char *filename, struct header h, struct gadpart *part);
unsigned int readgadget_novel(char *filename, struct header *h, fltarr **p, int **n, float **m);
unsigned int readgadget_sph(char *filename, struct header *h, fltarr **p, fltarr **v, int **n, float **m, float **u, float **rho, float **d1, float **d2, float **d3, float **d4 , float **sa);
int convertunits(struct header *head, struct gadpart *part, double convert_mass, double convert_distance);
double distance(fltarr a, fltarr b);
double distance_nopb(fltarr a, fltarr b);
double distbox(fltarr a, fltarr min, fltarr max);
double distbox_nopb(fltarr a, fltarr min, fltarr max);
void cpygadpart(gadpart * to, gadpart * from);
int gadsearch(gadpart_dist *data, double toFind, int start, int end);
double nfwfit(double *par, gadpart_dist *part, int cnt, double rv, double soft, double *rcs);
double densproffit(double *par, gadpart_dist *part, int cnt, double re, double soft, double *rcs, int type);
void calcdist(gadpart_dist *gd, int cnt, float *center);
void simplecenter(gadpart *part, int cnt, double* cm, int use);
void pcenter(gadpart_dist *part, int cnt, double, float*, int);
void findcenter(gadpart *part, int cnt, double maxdist, int use);
void simplecm(gadpart ** part, int cnt, float * cm);
double r200(gadpart_dist* pd, int cnt, double denscontrast, struct header h, int *vcnt, double *mvir);
double xoffset(gadpart_dist* pd, int cnt, float rvir, float* center);
#ifndef NOGAS
double temperature(const gadpart part);
#endif //NOGAS
void dummyfunction();
double galage(double z, double omegam, double omegal, double h);
double timediff(double z1, double z2, double omegam, double omegal, double h);
double a2z(double a);
double z2a(double z);
double angle(fltarr a, fltarr b);
double radvel(fltarr vel, fltarr rad);
struct header cphead(struct header head, gadpart* part, int cnt);

#ifdef GSL
void rotatepart(gadpart *part, int numpart, const gsl_matrix *rotmat);
void xrotate(double angle, gadpart *part, int numpart);
void yrotate(double angle, gadpart *part, int numpart);
void zrotate(double angle, gadpart *part, int numpart);
void rotategalaxy(gadpart *part, int numpart, double rad, int use, double *res, gsl_matrix **rotation);
#endif

#endif //LIBGAD
