#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#include "ompfuncs.h"
#endif

#include "OctTree.h"
#include "libgad.h"

static int num_threads=1;
static int periodic = 0;
static float boxsize = 0;


OctNode * initNode ()
{
  OctNode *onode;
  int i, j;
  onode = (OctNode*) malloc(sizeof(OctNode));
  onode->numpart = 0;
  for (j = 0; j < 3; j++)
    onode->center[j] = 0;
  onode->radius = 0;
  onode->part = NULL;
  for ( i=0; i<8; i++) 
    {
      onode->child[i] = NULL;
    }
  return onode;
}

void delOctNode (OctNode *onode)
{
  int i;
  for ( i = 0; i < 8; i++ )
    {
      if (onode->child[i]!=NULL)
	delOctNode(onode->child[i]);      
    }
  if (onode->numpart)
    {
      free(onode->part);
    }
  free(onode);
}

#ifdef _OPENMP
void set_num_threads()
{
#pragma omp parallel //private(num_threads)
  {
    int ithread= omp_get_thread_num();
    if (ithread == 0)
      {
	num_threads = omp_get_num_threads();
	//	printf("Number of Threads %d\n", num_threads);
      }
  }
}
#endif

void set_periodic_boundaries(struct header head)
{
  periodic =1;
  boxsize = head.boxsize;
}

int buildTreeBox(OctNode **onode, gadpart *part, struct header head, int maxparticles, int maxdepth)
{
  *onode = initNode();
  float radius = head.boxsize/2.;
  fltarr center={radius, radius, radius};
  unsigned int numpart_all = 0;
  int i;
  for ( i=0; i<6; i++) 
    {
      numpart_all += head.npart[i];
    }

  gadpart **pnt_part = (gadpart**) malloc (sizeof(gadpart*) * numpart_all);
#pragma omp parallel for
  for ( i = 0; i < numpart_all; i++) 
    {
      pnt_part[i] = &part[i];
    }

  set_num_threads();

  int ret_val;
  ret_val = buildTree(*onode, pnt_part, numpart_all, center, radius, maxparticles, maxdepth, 0);
//  free(pnt_part);
  return ret_val;
}

int buildTree(OctNode *onode, gadpart **part, unsigned int count, fltarr center, float radius, unsigned int leafsize, unsigned int maxDepth, unsigned int currentDepth )
{
#ifdef DEBUG
  printf("Octree.buildtree: count %d\t| center %f %f %f | radius %f | currentDepth %d \n", count, center[0], center[1], center[2], radius, currentDepth);
  fflush(stdout);
#endif
  int j;
  for (j=0; j<3; j++)
    onode -> center[j] = center[j];
  onode -> radius = radius;

  if (count <= leafsize || currentDepth >= maxDepth)
    {
      int i;
      onode -> numpart = count;
      onode -> part = (gadpart **) malloc (sizeof(gadpart*) * count);
      for ( i = 0; i < count; i++ )
	{
	  onode -> part[i] = part[i];
	}	
      free(part);
      return 1;
    }
  int i;
  unsigned int *childPointCounts[8];
  unsigned int *code = (int *) calloc (count, sizeof(int));
#pragma omp parallel for private(i)
  for ( i = 0; i < 8; i++ )
    {
      childPointCounts[i] = (int*) calloc (num_threads, sizeof(int));
    }

#pragma omp parallel for private(i)
  for ( i = 0; i < count; i++ )
    {
      int ithread = 0;
#ifdef _OPENMP
      ithread = omp_get_thread_num();
#endif
      code[i] = 0;

      if (part[i]->pos[0] > center[0]) code[i] |= 1;
      if (part[i]->pos[1] > center[1]) code[i] |= 2;
      if (part[i]->pos[2] > center[2]) code[i] |= 4; 

      childPointCounts[code[i]][ithread]++;
//      if (code[i]==7)
//	{
//	  printf("pos2 %g code %d | ithread %d\n", part[i]->pos[0], code[i], ithread);fflush(stdout);
//	  printf("cpc 7: %d\n",  childPointCounts[7][0]);fflush(stdout);
//	}
    }

  //  printf("num_treads %d\n", num_threads);fflush(stdout);
  for ( i = 0; i < 8; i++ )
    {
      for ( j = 1; j<num_threads; j++)
	{
	  childPointCounts[i][0] += childPointCounts[i][j];
	  //	  	  printf("childPointcount %d %d | %d\n", i, j, childPointCounts[i][j]);fflush(stdout);
	}
      //      printf("childPointcount0 %d\n", childPointCounts[i][0]);fflush(stdout);
    }

  gadpart **newlist[8];
  unsigned int newcount[8];

#pragma omp parallel for // firstprivate(radius)
  for ( i = 0; i < 8; i++ )
  {
    newcount[i] = 0;
    if (!childPointCounts[i][0]) continue;
    onode->child[i] = initNode();
    newlist[i] = (gadpart**) malloc (sizeof(gadpart*) * childPointCounts[i][0]);
    if (newlist==NULL)
          {
                  fprintf(stderr, "OctTree: fail to allocate memory in buildtree!\n");
                  exit(1);
          }
    int n;
    for ( n = 0; n < count; n++)
  	   {
  	   if (code[n]==i)
  	       {
  	           newlist[i][newcount[i]++] = part[n];
  	       }
  	   }
  }

  free(part);

#pragma omp parallel  for // firstprivate(radius)
  for ( i = 0; i < 8; i++ )
  {
    if (childPointCounts[i][0])
    {
      int n;
      fltarr newcenter = {-0.5,-0.5,-0.5};
      if (i & 1) newcenter[0] = 0.5;
      if (i & 2) newcenter[1] = 0.5;
      if (i & 4) newcenter[2] = 0.5;
      for ( n = 0; n < 3; n++ )
    	   newcenter[n] = center[n] + newcenter[n] * radius;
      float newradius = radius *  0.5;
      buildTree(onode->child[i], newlist[i], newcount[i], newcenter, newradius, leafsize, maxDepth, currentDepth+1);
 //     free(newlist[i]);
    }
    free(childPointCounts[i]);
  }

  free(code);
  return 1;
}

int checkOctTree(OctNode *onode)
{
  if (onode->numpart)
    return onode->numpart;
  int sum=0;
  int i;
  for ( i = 0; i < 8; i++ )
    {
      if (onode->child[i] != NULL)
	sum += checkOctTree(onode->child[i]);
    }
  return sum;
}

float distOctNode(OctNode *onode, fltarr center)
{
  float dist = 0;
  int dim;
  for ( dim=0; dim<3; dim++)
    {
      float dum=ABS(center[dim] - onode->center[dim]);
      if (periodic)
	dum = MIN(dum, boxsize-dum);
      if (dum < onode->radius) 
	continue;
//      if ((periodic) && (dum > (boxsize - onode->radius)) )
//	continue;
      dist += SQR(dum - onode->radius);
    }
  return sqrt(dist);  
}

float distPart(gadpart *part, fltarr pos)
{
  float dist = 0;
  int dim;
  for ( dim=0; dim<3; dim++)
    {
      float dum=ABS(pos[dim] - part->pos[dim]);
      if (periodic)
	dum = MIN(dum, boxsize-dum);
      dist += SQR(dum);
    }
  return sqrt(dist);  
}

unsigned long findParticles(OctNode *onode, fltarr pos, float dist, gadpart ***result, unsigned int *numpart, unsigned int *size)
{
  if ( distOctNode(onode, pos) > dist )
    return 0;  
  int i;
  int sum = 0;
  for ( i = 0; i < 8; i++ )
    {
      if (onode->child[i] == NULL)
	continue;
      sum += findParticles(onode->child[i], pos, dist, result, numpart, size);
    }
  int npart = 0;
  if (onode->numpart)
    {
      if (*size==0)
	{
	  //	  printf("allocating...\n");fflush(stdout);
	  *size = PBUFF;
	  *result = (gadpart**) calloc (*size, sizeof(gadpart*));

	} 
      if ((*numpart+onode->numpart) > *size)
	{
	  //	  printf("reallocating...\n");fflush(stdout);
	  while ( (*numpart+onode->numpart) > *size)
	    (*size)*=2;
	  *result = realloc(*result, (*size)*sizeof(gadpart*));
	}
      for ( i = 0; i < onode->numpart; i++ )
	{
	  if ( distPart(onode->part[i], pos) < dist )
	    {
	      (*result)[(*numpart)++] = (onode->part[i]);
	      npart++;
	    }
	}
    }
  return sum+npart;
}

OctNode * findLeaf(OctNode *onode, fltarr pos)
{
  OctNode *node = onode;
  while ( node->numpart == 0 ) 
    {
      int code = 0;
      if (pos[0] > node->center[0]) code |= 1;
      if (pos[1] > node->center[1]) code |= 2;
      if (pos[2] > node->center[2]) code |= 4; 
      if ( node->child[code] == NULL )
	{
	  return node;
	}
      node = node->child[code];
    }
  return node;
}

void rejectNodes(OctNode *onode, int atype ,int rtype, float reject_ratio)
{
  int i,j;
  if (onode->numpart)
    {
      float rejects = 0;
      float accepts = 0;
      for ( i = 0; i < onode->numpart; i++)
	{
	  if ((1<<(onode->part[i]->type)) & atype )
	    accepts++;
	  if ((1<<(onode->part[i]->type)) & rtype )
	    rejects++;
	}
      if ((accepts == 0) && (rejects == 0))
	{
	    return;
	}
      if ( rejects/(accepts+rejects) > reject_ratio )
	{
	  onode->numpart = 0;
	  free(onode->part);
	}
      return;
    }

#pragma omp parallel for // firstprivate(radius)
  for ( i = 0; i < 8; i++ )
  {
    if (onode->child[i] == NULL)
      continue;
    rejectNodes(onode->child[i], atype, rtype, reject_ratio);
  }  
  
  return;
}

