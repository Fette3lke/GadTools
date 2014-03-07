#ifndef KDTREE
#define KDTREE

#include "libgad.h"
#include <math.h>

#define PBUFF 30000
#define MAX_ITERATIONS 100
#define MAXINT 10000000


typedef struct KdNode {
  struct KdNode * parent;
  struct KdNode * up;
  struct KdNode * down;
  gadpart * part;
  int id;
  float min;
  float max;
  short int dim;
} KdNode;

KdNode * initKdNode (KdNode ** p, KdNode * parent);
void delKdNode(KdNode **node);
void buildKdTree(KdNode * root, gadpart * partarr, unsigned int numpart, short int dim);
unsigned int checkKdTree(KdNode * root);
KdNode * findNN(KdNode * root, gadpart * part);
KdNode * KdNeighbor(KdNode * root, gadpart * part, double * dist); 
double findkNN(KdNode * root, gadpart * part, double dist, gadpart_dist ** result, int k);
void findparts(KdNode * root, fltarr pos, double dist, gadpart_dist ** result, unsigned int * numpart,unsigned int * size);
void findpartsPB(KdNode * root, fltarr pos, double dist, gadpart_dist ** result, unsigned int * numpart,unsigned int * size);
void findGadpartsPB(KdNode * root, fltarr pos, double dist, gadpart *** result, int * numpart, unsigned int * size);
void findGadparts(KdNode * root, fltarr pos, double dist, gadpart *** result, int * numpart, unsigned int * size);
void findNewGadparts(KdNode * root, fltarr pos, double dist, gadpart *** result, int * numpart, unsigned int * size);
void findFOF(KdNode * root, fltarr pos, double ll, gadpart *** result, int * numpart, unsigned int * size);
int checkTree(KdNode *root);
#endif
