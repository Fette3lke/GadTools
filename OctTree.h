#ifndef OCTTREE
#define OCTTREE

#include "libgad.h"
#include <math.h>

#define PBUFF 32768

typedef struct OctNode {
  struct OctNode *child[8];
  fltarr center;
  float radius;
  gadpart **part;
  unsigned int numpart;
} OctNode;

#ifdef _OPENMP
void set_num_threads();
#endif


OctNode * initNode();
void delOctNode (OctNode *onode);
void set_periodic_boundaries(struct header head);
int buildTreeBox(OctNode **onode, gadpart *part, struct header head, int maxparticles, int maxdepth);
int buildTree (OctNode *onode, gadpart **part, unsigned int count, fltarr center, float radius, unsigned int leafsize, unsigned int maxDepth, unsigned int currentDepth );
int checkOctTree(OctNode *onode);
float distOctNode(OctNode *onode, fltarr center);
float distPart(gadpart *part, fltarr pos);
unsigned long findParticles(OctNode *onode, fltarr pos, float dist, gadpart ***result, unsigned int *numpart, unsigned int *size);
OctNode * findLeaf(OctNode *onode, fltarr pos);
void rejectNodes(OctNode *onode, int atype, int rtype, float reject_ratio);
#endif
