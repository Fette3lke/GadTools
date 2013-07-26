/*
  Actually it's not a KD-Tree but a 3D-Tree
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#include "ompfuncs.h"
#endif

#include "KdTree.h"
#include "libgad.h"


KdNode * initKdNode (KdNode ** p, KdNode * parent)
{
    *p=(KdNode *) malloc(sizeof(KdNode));
    if (*p==NULL) return NULL;
    (*p)->down=NULL;
    (*p)->up  =NULL;
    (*p)->parent =parent;
    (*p)->min=0;
    (*p)->max=0;
    (*p)->part=NULL;
    (*p)->id=-1;
    (*p)->dim =-1;
    return *p;
}

void delKdNode(KdNode **node)
{
  KdNode * tmpNode = *node;
  if ( (*node) -> id == 0)
    {
      fprintf(stderr, "STOP %d! delKdNode error\n", (*node) -> id);
      exit(1);
    }
    if ((*node)->up  !=NULL) delKdNode(&((*node)->up));
    if ((*node)->down!=NULL) delKdNode(&((*node)->down));
    if ((*node)->parent!=NULL) 
      {
	if ( (*node)->parent->up == (*node) ) 
	  {
	    (*node)->parent->up = NULL;
	  }
	else if ( (*node)->parent->down == (*node) ) 
	  {
	    (*node)->parent->down = NULL;
	  }
	else 
	  {
	    fprintf(stderr, "There is something wrong with your tree!\n");
//	    fprintf(stderr, "%d %d\n", (*node)->id, (*node)->parent->id);
//	    fprintf(stderr, "%f %d\n", (*node)->parent->min, (*node)->parent->dim);
//	    if ( (*node)->parent->part == NULL )
//	      	    fprintf(stderr, "kein parent part\n");
//	    fprintf(stderr, "%f %d\n", (*node)->part->pos[0], (*node)->part->id);
//	    fprintf(stderr, "%d\n", checkKdTree(*node));
	    exit(1);
	  }
      }
    free(tmpNode);
    tmpNode=NULL;
}

#ifndef _OPENMP
void buildKdTree(KdNode * root, gadpart * partarr, unsigned int numpart, short int dim)
{
  if (numpart>1)
    {
      extern int search_dim;
      search_dim=dim;
      qsort(partarr, numpart, sizeof(gadpart), cmp_pos);
      int median=numpart/2;
      root->part=&(partarr[median]);
      root->id=(partarr[median].id);
      root->min=partarr[0].pos[dim];
      root->max=partarr[numpart-1].pos[dim];
      root->dim=dim;
	if ((median)>0) 
	  {
	    if (initKdNode(&(root->down),root)==NULL) {fprintf(stderr,"failed to allocate memory\n");exit(1);}
	    buildKdTree(root->down, &(partarr[0]), (median), (dim+1)%3);
	  }
	if ((numpart-(median+1))>0) 
	  {
	    if (initKdNode(&(root->up),  root)==NULL) {fprintf(stderr,"failed to allocate memory\n");exit(1);}
	    buildKdTree(root->up  , &(partarr[median+1]), (numpart-(median+1)), (dim+1)%3);
	  }
    } 
  else 
    {
      root->part=&(partarr[0]);
      root->id=(partarr[0].id);
      root->min=partarr[0].pos[dim];
      root->max=partarr[0].pos[dim];
      root->dim=dim;
    }
}
#else
void buildKdTree(KdNode * root, gadpart * partarr, unsigned int numpart, short int dim)
{
  if (numpart>1)
    {
      void *a, *b;
      int median;
    switch (dim)
      {
      case 0:
	median=partarray(partarr, numpart, sizeof(gadpart), cmp_x, &a, &b);
	break;
      case 1:
	median=partarray(partarr, numpart, sizeof(gadpart), cmp_y, &a, &b);
	break;
      case 2:
	median=partarray(partarr, numpart, sizeof(gadpart), cmp_z, &a, &b);
	break;
      }
      gadpart *min=(gadpart *)a;
      gadpart *max=(gadpart *)b;
      root->part=&(partarr[median]);
      root->id=(partarr[median].id);
      root->min=min->pos[dim];
      root->max=max->pos[dim];
      if ((max->pos[dim] == partarr[median].pos[dim]) && (median>2))
	{
	  int j;
	  for (j=median; j<numpart; j++)
	    {
	      if (partarr[j].pos[dim] < partarr[median].pos[dim])
		{
		  printf("lower than median\n");
		  exit(1);
		}
	      if (partarr[j].pos[dim] > max->pos[dim])
		{
		  printf("larger than max\n");
		  exit(1);
		}
	    }

	}
      root->dim=dim;
#pragma omp parallel sections
      {
#pragma omp section      
	if ((median)>0) 
	  {
	    if (initKdNode(&(root->down),root)==NULL) {fprintf(stderr,"failed to allocate memory\n");exit(1);}
	    buildKdTree(root->down, partarr, (median), (dim+1)%3);
	  }
#pragma omp section
	if ((numpart-(median+1))>0) 
	  {
	    if (initKdNode(&(root->up),  root)==NULL) {fprintf(stderr,"failed to allocate memory\n");exit(1);}
	    buildKdTree(root->up  , &(partarr[median+1]), (numpart-(median+1)), (dim+1)%3);
	  }
      }
    } 
  else 
    {
      root->part=&(partarr[0]);
      root->id=(partarr[0].id);
      root->min=partarr[0].pos[dim];
      root->max=partarr[0].pos[dim];
      root->dim=dim;
    }
}
#endif //_OPENMP


//returns distance of fltarr to a Tree-node, 0 if inside

double distKdNode(fltarr pos, KdNode * node)
{
  fltarr min, max;
  short int dim = node->dim;
  min[dim]= node->min;
  max[dim]= node->max;
  if (node->parent!=NULL)
    {
      dim = node->parent->dim;
      min[dim]=node->parent->min;
      max[dim]=node->parent->max;

      if ( node->parent->up == (node))
	{
	  min[dim] = node->parent->part->pos[dim];
	}
      else 
	{
	  max[dim] = node->parent->part->pos[dim];
	}
      
      if (node->parent->parent!=NULL)
	{
	  dim = node->parent->parent->dim;
	  min[dim]=node->parent->parent->min;
	  max[dim]=node->parent->parent->max;
	  if ( node->parent->parent->up == (node->parent))
	    {
	      min[dim] = node->parent->parent->part->pos[dim];
	    }
	  else 
	    {
	      max[dim] = node->parent->parent->part->pos[dim];
	    }
	  return distbox_nopb(pos, min, max);	 
	} 
      else
	{
	  return 0;
//	  double dist;
//	  dist =SQR(ABS(pos[dim] - node->part->pos[dim]));
//	  double tmp = MIN(ABS(pos[node->parent->dim] - node->parent->part->pos[node->parent->dim]), ABS(pos[node->parent->dim] - node->part->pos[node->parent->dim]));
//	  dist+=SQR(tmp);
//	  return sqrt(dist);
	}
    } else return 0;

}

unsigned int checkKdTree(KdNode * root)
{
  float pos = root->part->pos[root->dim];
  float dist = distKdNode( root->part->pos, root);
  int up = 0;
  if (root->parent != NULL)
    {
      if (root->parent->up == root)
	{
	  up = 1;
	}

      int dim = root->parent->dim;
      pos = root->part->pos[dim];
      if ( ( (up) && ( pos < root-> parent->part-> pos[dim] ) ) || ( (!up) && ( pos > root-> parent->part-> pos[dim] ) ) || (dist > 0))
	{
	  KdNode *tmp=root;
	  int lvls=0;
	  while (tmp->parent != NULL)
	    {
	      tmp = tmp->parent;
	      lvls++;
	    }
	  fprintf(stderr, "particle out of bounds..%d...%f...%d!\n",root->part->id, dist, lvls);
	  fprintf(stderr, "parent..%f...pos... %f...up....%d!\n",root->parent->part->pos[dim], pos, up);
	  fprintf(stderr, "pos... %f %f %f\n",root->part->pos[0],root->part->pos[1],root->part->pos[2]);
	  fprintf(stderr, "min %f   max %f  dim %d\n",root->min,root->max,root->dim);
	  fprintf(stderr, "min %f   max %f  dim %d parent %f\n",root->parent->min,root->parent->max,root->parent->dim, root->parent->part->pos[root->parent->dim]);          
	  fprintf(stderr, "min %f   max %f  dim %d parent %f\n",root->parent->parent->min,root->parent->parent->max,root->parent->parent->dim, root->parent->parent->part->pos[root->parent->parent->dim]);

	  if (root->up != NULL)
	    fprintf(stderr, "up\n");
	  if (root->down != NULL)
	    fprintf(stderr, "down\n");
	  exit(1);
	}
    }
  if (((root->up)!=NULL) && ((root->down)!=NULL)) return (checkKdTree(root->up)+checkKdTree(root->down)+1);
  if ((root->up)!=NULL)   return (checkKdTree(root->up)+1);
  if ((root->down)!=NULL) return (checkKdTree(root->down)+1);
  return 1;
}

double distKdNodePB(fltarr pos, KdNode * node)
{
  fltarr min, max;
  short int dim = node->dim;
  min[dim]= node->min;
  max[dim]= node->max;
  if (node->parent!=NULL)
    {
      min[node->parent->dim]=node->parent->min;
      max[node->parent->dim]=node->parent->max;
      if (node->parent->parent!=NULL)
	{
	  min[node->parent->parent->dim]=node->parent->parent->min;
	  max[node->parent->parent->dim]=node->parent->parent->max;
	  return distbox(pos, min, max);
	} 
      else
	{
	  return 0;
//	  double dist;
//	  dist =SQR(ABS(pos[dim] - node->part->pos[dim]));
//	  double tmp = MIN(ABS(pos[node->parent->dim] - node->parent->part->pos[node->parent->dim]), ABS(pos[node->parent->dim] - node->part->pos[node->parent->dim]));
//	  dist+=SQR(tmp);
//	  return sqrt(dist);
	}
    } else return 0;

}

//Find Nearest Neighbour of a particle

KdNode * findNN(KdNode * root, gadpart * part)
{
  double dist=-1;
  return KdNeighbor(root, part, &dist);
}

KdNode * KdNeighbor(KdNode * root, gadpart * part, double * dist)
{
  double distr=distance_nopb(root->part->pos, part->pos);
  if (((distr < *dist) || (*dist<=0)) && (distr>0)) *dist=distr;
  short int dim = root->dim;
  fltarr min, max;
  if (distr==0) distr=root->max - root->min;
  double distup=0;
  double distdown=0;
  if (root->up  !=NULL) distup  =distKdNode(part->pos, root->up);
  if (root->down!=NULL) distdown=distKdNode(part->pos, root->down);

  if ((root->up!=NULL) && (root->down!=NULL))
    {
      //     printf("%f %f %f\n", distup, distdown, distr);
      if (((distup < *dist ) && (distdown < *dist )) || (*dist<=0))
	{
	  KdNode *a=KdNeighbor(root->up  , part, dist);
  	  KdNode *b=KdNeighbor(root->down, part, dist);
	  //	  (Vorsicht, glaub noch nicht, dass das stimmt) scheint zu funktionieren
	  double dista=distance_nopb(a->part->pos, part->pos);
	  double distb=distance_nopb(b->part->pos, part->pos);

	  //    printf("%f %f %f %f\n",distr, dista, distb, *dist);
	  
	  if (distr==0) distr=dista+distb;
	  if (dista==0) dista=distr+distb;
	  if (distb==0) distb=dista+distr;
	  
	  if ((dista<distb) && (dista<distr)) return a;
	  if ((distb<dista) && (distb<distr)) return b;
	  return root;
	}
      else if (distup < *dist )
	{
	  KdNode *a=KdNeighbor(root->up  , part, dist);
	  double dista=distance_nopb(a->part->pos, part->pos);
	  if (dista < distr) return a;
	  return root;
	}
      else if (distdown < *dist ) 
	{
	  KdNode *b=KdNeighbor(root->down  , part, dist);
	  double distb=distance_nopb(b->part->pos, part->pos);
	  if (distb < distr) return b;
	  return root;
	}
    }
  else if (root->up  !=NULL) 
    {
      KdNode *a=KdNeighbor(root->up  , part, dist);
      double dista=distance_nopb(a->part->pos, part->pos);
      if (dista < distr) return a;
      return root;
    }
  else if (root->down!=NULL) 
    {
      KdNode *b=KdNeighbor(root->down  , part, dist);
      double distb=distance_nopb(b->part->pos, part->pos);
      if (distb < distr) return b;
      return root;
    }
  
  return root;
}

//find k nearest Neighbours, returns distance of farest found neighbour
double findkNN(KdNode * root, gadpart * part, double dist, gadpart_dist ** result, int k)
{
  unsigned int num=0, size=0, iter=0;
  if (dist==0) return 0;
  while (num<k)
    {
      if (num!=0) free(*result);
      num=0;
      size=0;
      findparts(root, part->pos, dist, result, &num, &size);
      dist+=dist;
      iter++;
      if (iter > MAX_ITERATIONS)
	{
	  fprintf(stderr, "findkNN failed!!\n");
	  free(*result);
	  return 0;
	}
    }

  *result=realloc(*result, num*sizeof(gadpart_dist));

  qsort(*result, num, sizeof(gadpart_dist),cmp_dist);
  *result=realloc(*result, k*sizeof(gadpart_dist));
  return (*result)[k-1].dist;
}

//find particles closer than dist to pos -> stored in *result, numpart

void findparts(KdNode * root, fltarr pos, double dist, gadpart_dist ** result, unsigned int * numpart, unsigned int * size)
{
  double distup=0;
  double distdown=0;
  if (root->up  !=NULL) 
    {
      distup  =distKdNode(pos, root->up);
      if (distup   < dist) findparts(root->up  , pos, dist, result, numpart, size);
    }
  if (root->down!=NULL) 
    {
      distdown=distKdNode(pos, root->down);
      if (distdown < dist) findparts(root->down, pos, dist, result, numpart, size);
    }
  double distance=distance_nopb(pos, root->part->pos);
  if ((distance < dist) && (distance>0))
    {
      if (!(*numpart)) 
	{
	  *result = (gadpart_dist *) malloc (PBUFF*sizeof(gadpart_dist));
	  if (*result==NULL) {fprintf(stderr,"failed to allocate memory\n");exit(1);}
	  *size=PBUFF;
	  *numpart=1;
	}
      else
	{
	  (*numpart)++;
	  if (*numpart>=*size) 
	    {
	      *size+=PBUFF;
	      *result=realloc(*result, (*size)*sizeof(gadpart_dist));
	    }
	}
      //      memcpy(&((*result)[(*numpart)-1].part), root->part, sizeof(gadpart));
      cpygadpart(&((*result)[(*numpart)-1].part)  , root->part);
      (*result)[(*numpart)-1].dist=distance_nopb(root->part->pos, pos);
      //      (*result)[(*numpart)-1]=root->part;
    }
  return ;
}

void findpartsPB(KdNode * root, fltarr pos, double dist, gadpart_dist ** result, unsigned int * numpart, unsigned int * size)
{
  double distup=0;
  double distdown=0;
  if (root->up  !=NULL) 
    {
      distup  =distKdNodePB(pos, root->up);
      if (distup   < dist) findpartsPB(root->up  , pos, dist, result, numpart, size);
    }
  if (root->down!=NULL) 
    {
      distdown=distKdNodePB(pos, root->down);
      if (distdown < dist) findpartsPB(root->down, pos, dist, result, numpart, size);
    }
  double Dist=distance(pos, root->part->pos);
  if ((Dist < dist) && (Dist>0))
    {
      if (!(*numpart)) 
	{
	  *result = (gadpart_dist *) malloc (PBUFF*sizeof(gadpart_dist));
	  if (*result==NULL) {fprintf(stderr,"failed to allocate memory\n");exit(1);}
	  *size=PBUFF;
	  *numpart=1;
	}
      else
	{
	  (*numpart)++;
	  if (*numpart>=*size) 
	    {
	      *size+=PBUFF;
	      *result=realloc(*result, (*size)*sizeof(gadpart_dist));
	    }
	}
      //      memcpy(&((*result)[(*numpart)-1].part), root->part, sizeof(gadpart));
      cpygadpart(&((*result)[(*numpart)-1].part)  , root->part);
      (*result)[(*numpart)-1].dist=distance(root->part->pos, pos);
      //      (*result)[(*numpart)-1]=root->part;
    }
  return ;
}

void findGadpartsPB(KdNode * root, fltarr pos, double dist, gadpart *** result, int * numpart, unsigned int * size)
{
  //  printf("!");
  double distup=0;
  double distdown=0;
  if (root->up  !=NULL) 
    {
      distup  =distKdNodePB(pos, root->up);
      //      printf(" du%f ", distup);
      if (distup   < dist) findGadpartsPB(root->up  , pos, dist, result, numpart, size);
    }
  if (root->down!=NULL) 
    {
      distdown=distKdNodePB(pos, root->down);
      //      printf(" dd%f ", distdown);
      if (distdown < dist) findGadpartsPB(root->down, pos, dist, result, numpart, size);
    }
  double Dist=distance(pos, root->part->pos);
  if ((Dist < dist) && (Dist>0))
    {
      if (!(*numpart)) 
	{
	  *result = (gadpart **) malloc (PBUFF*sizeof(gadpart*));
	  if (*result==NULL) {fprintf(stderr,"failed to allocate memory\n");exit(1);}
	  *size=PBUFF;
	  *numpart=1;
	}
      else
	{
	  (*numpart)++;
	  if (*numpart>=*size) 
	    {
	      *size+=PBUFF;
	      *result=realloc(*result, (*size)*sizeof(gadpart*));
	    }
	}
      //      memcpy(&((*result)[(*numpart)-1].part), root->part, sizeof(gadpart));
      //      cpygadpart(&((*result)[(*numpart)-1].part)  , root->part);
      //      (*result)[(*numpart)-1].dist=distance(root->part->pos, pos);
      //      (*result)[(*numpart)-1]=root->part;
      (*result)[(*numpart)-1]=(root->part);
    }
  return ;
}

void findGadparts(KdNode * root, fltarr pos, double dist, gadpart *** result, int * numpart, unsigned int * size)
{
  //  printf("!");
  double distup=0;
  double distdown=0;
  if (root->up  !=NULL) 
    {
      distup  =distKdNode(pos, root->up);
      //      printf(" du%f ", distup);
      if (distup   < dist) findGadparts(root->up  , pos, dist, result, numpart, size);
    }
  if (root->down!=NULL) 
    {
      distdown=distKdNode(pos, root->down);
      //      printf(" dd%f ", distdown);
      if (distdown < dist) findGadparts(root->down, pos, dist, result, numpart, size);
    }
  double Dist=distance_nopb(pos, root->part->pos);
  if ((Dist < dist) && (Dist>0))
    {
      if (!(*numpart)) 
	{
	  *result = (gadpart **) malloc (PBUFF*sizeof(gadpart*));
	  if (*result==NULL) {fprintf(stderr,"failed to allocate memory\n");exit(1);}
	  *size=PBUFF;
	  *numpart=1;
	}
      else
	{
	  (*numpart)++;
	  if (*numpart>=*size) 
	    {
	      *size+=PBUFF;
	      *result=realloc(*result, (*size)*sizeof(gadpart*));
	    }
	}
      //      memcpy(&((*result)[(*numpart)-1].part), root->part, sizeof(gadpart));
      //      cpygadpart(&((*result)[(*numpart)-1].part)  , root->part);
      //      (*result)[(*numpart)-1].dist=distance(root->part->pos, pos);
      //      (*result)[(*numpart)-1]=root->part;
      (*result)[(*numpart)-1]=(root->part);
    }
  return ;
}

void findNewGadparts(KdNode * root, fltarr pos, double dist, gadpart *** result, int * numpart, unsigned int * size)
{
  //    printf("!");
  double distup=0;
  double distdown=0;
  if (root->up  !=NULL) 
    {
      distup  =distKdNode(pos, root->up);
      //      printf(" du%f ", distup);
      if (distup   < dist) findNewGadparts(root->up  , pos, dist, result, numpart, size);
    }
  if (root->down!=NULL) 
    {
      distdown=distKdNode(pos, root->down);
      //      printf(" dd%f ", distdown);
      if (distdown < dist) findNewGadparts(root->down, pos, dist, result, numpart, size);
    }
  double Dist=distance_nopb(pos, root->part->pos);
  if ((Dist < dist) && (Dist>0))
    {
      gadpart** fnd = NULL;
      fnd = bsearch(&(root -> part), (*result), *numpart, sizeof(gadpart*),  cmp_pointer_id);
      if ( fnd != NULL ) return;

      if (!(*numpart)) 
	{
	  *result = (gadpart **) malloc (PBUFF*sizeof(gadpart*));
	  if (*result==NULL) {fprintf(stderr,"failed to allocate memory\n");exit(1);}
	  *size=PBUFF;
	  *numpart=1;
	}
      else
	{
	  (*numpart)++;
	  //	  printf("%d ", *numpart);
	  if (*numpart>=*size) 
	    {
	      *size+=PBUFF;
	      *result=realloc(*result, (*size)*sizeof(gadpart*));
	    }
	}
      //      memcpy(&((*result)[(*numpart)-1].part), root->part, sizeof(gadpart));
      //      cpygadpart(&((*result)[(*numpart)-1].part)  , root->part);
      //      (*result)[(*numpart)-1].dist=distance(root->part->pos, pos);
      //      (*result)[(*numpart)-1]=root->part;
      int i = *numpart-2;
      while ( (i>=0) && ( root->part->id < (*result)[i]->id ) )
	{
	  (*result)[i+1] = (*result)[i];
	  i--;
	}
      //      printf("insert: %d | %d\n", *numpart-i, root->part->id);
      (*result)[i+1]=(root->part);
//      if (i<*numpart-3)
//	printf("%d %d %d\n", (*result)[i]->id, (*result)[i+1]->id, (*result)[i+2]->id);
//      else 
//	printf("%d %d\n", (*result)[i]->id, (*result)[i+1]->id);
      //      (*result)[(*numpart)-1]=(root->part);
    }
  return ;
}

//find FOF-group around pos, dist = linking length

static double maxdist = 0;

static void FindFOF(KdNode * root, fltarr pos, double dist, gadpart *** result, int * numpart, unsigned int * size)
{
  double distup=0;
  double distdown=0;
  int i;
  if (root->up  !=NULL) 
    {
      distup  =distKdNode(pos, root->up);
      //      printf(" du%f ", distup);
      if (distup   < ( dist + maxdist )) FindFOF(root->up  , pos, dist, result, numpart, size);
    }
  if (root->down!=NULL) 
    {
      distdown=distKdNode(pos, root->down);
      //      printf(" dd%f ", distdown);
      if (distdown < ( dist + maxdist )) FindFOF(root->down, pos, dist, result, numpart, size);
    }
  double Dist=distance_nopb(pos, root->part->pos);
  if ((Dist > maxdist + dist) || ( root -> id == -1) ) return;
  else 
    {
      gadpart **fnd;
      fnd = bsearch(&(root -> part), (*result), *numpart, sizeof(gadpart*),  cmp_pointer_id);
      if (( fnd != NULL) )// && (root -> parent !=NULL ) )//&& (root -> parent ->id))
	{
	  root -> id = -1;
	  //	  printf("fnd: %d %d %d %d\n", root -> id, root -> part -> id, (*fnd) ->  id, root -> parent ->id);
	  if (root -> parent !=NULL )
	  if( (root -> down == NULL) && (root -> up == NULL) && (root -> parent != NULL) ) 
	    {
	      //	    if (root -> parent -> down -> id == root -> id) delKdNode( &(root -> parent -> down ) );
	      //	    if (root -> parent -> up   -> id == root -> id) delKdNode( &(root -> parent -> up   ) );
	      //	      printf("*%d %d %d\n", root -> id, root ->parent -> id, root -> parent -> dim); fflush (stdout);

	      if ((root -> parent -> down != NULL) && (root -> parent -> down -> id == root -> id)) 
		{
		  delKdNode( &( root -> parent -> down ) );	  
		}
	      else if ( (root -> parent -> up != NULL) && (root -> parent -> up   -> id == root -> id)) 
		{
		  delKdNode( &(root -> parent -> up   ) );
		}
	    }
	  return;
      }
    }
  int cnt=0;
//   for ( i = 0; i < *numpart-1; i++ )
//	{
//	  if ((*result)[i]->id > (*result)[i+1]->id  )
//	    {
//	      printf("!#!#!\n");
//	      exit(1);
//	    }
//	  else cnt++;
//	}
//   printf("###%d %d### %d\n", cnt, *numpart, (*result)[*numpart-1] -> id);
  for ( i = 0; i < *numpart; i++ )
	{	 
	  double ddum = distance_nopb( (*result)[i]->pos, root->part->pos);
	  if ((ddum < Dist) ) Dist = ddum;
	  if (Dist < 1e-7) break;
	  if (Dist < dist) 
	    {
	      //printf("%d %d %f %f %d %d\n", *numpart, i, Dist, dist, root->id, (*result)[i]->id);
	      break;
	    }
	}
  if ((Dist < dist) && (Dist > 1e-7))
    {
      double ddum = distance_nopb(pos, root->part->pos);
      if (( ddum > maxdist )) maxdist = ddum;
      if (!(*numpart)) 
	{
	  *result = (gadpart **) malloc (PBUFF*sizeof(gadpart*));
	  if (*result==NULL) {fprintf(stderr,"failed to allocate memory\n");exit(1);}
	  *size=PBUFF;
	  *numpart=1;
	}
      else
	{
	  (*numpart)++;
	  if (*numpart>=*size) 
	    {
	      *size+=PBUFF;
	      *result=realloc(*result, (*size)*sizeof(gadpart*));
	    }
	}
      //      memcpy(&((*result)[(*numpart)-1].part), root->part, sizeof(gadpart));
      //      cpygadpart(&((*result)[(*numpart)-1].part)  , root->part);
      //      (*result)[(*numpart)-1].dist=distance(root->part->pos, pos);
      //      (*result)[(*numpart)-1]=root->part;

      //Insert particle in list

//      if ((*numpart == 1) || (root -> part -> id > (*result)[*numpart - 2] -> id) )
//	(*result)[(*numpart)-1]=(root->part);
//      else
//	for ( i = 0; i < (*numpart-1); i++)
//	  {
//	    if (root -> part -> id < (*result)[i] -> id)
//	      {
//		int j;
//		for ( j = (*numpart-1); j > i; j--)
//		  {
//		    (*result)[j] = (*result)[j-1];
//		  }		  
//		(*result)[i] = root -> part;
//		break;
//	      }
//	  }

      int j = *numpart - 2;
      while ((j>=0) && (root -> part -> id < (*result)[j] -> id))
	{
	  (*result)[j+1] = (*result)[j];
	  j--;
	}
      j++;
      (*result)[j] = root -> part;

      //      printf("new fof: %d %f \n", root -> id, Dist);

      if( (root -> down == NULL) && (root -> up == NULL) && (root -> parent != NULL) )
	{
	  //	  KdNode *tempNode = root -> parent ;
	  if ((root -> parent -> down != NULL) && (root -> parent -> down -> id == root -> id)) delKdNode( &( root -> parent -> down ) );	  
	  else if ( (root -> parent -> up != NULL) && (root -> parent -> up   -> id == root -> id)) delKdNode( &(root -> parent -> up   ) );
//	  if ((tempNode != NULL) && (tempNode -> parent != NULL)) 
//	    FindFOF(tempNode->parent, pos, dist, result, numpart, size);
	}
    }
  return ;
}

//void  findFOF(KdNode * root, fltarr pos, double dist, gadpart *** result, int * numpart, unsigned int * size)
//{
//  int nold = -1;
//  if (*numpart >= 1)
//    {
//      int i;
//      maxdist = 0;
//      for ( i = 0; i < *numpart; i++)
//	{
//	  double ddum = distance_nopb(pos, (*result)[i] -> pos);
//	  if (ddum > maxdist) maxdist = ddum;
//	}
//    } else maxdist = 0;
//  while (nold != *numpart)
//    {
//      nold = *numpart;
//      //printf("nold %d\n", nold); fflush(stdout);
//      FindFOF(root, pos, dist, result, numpart, size);
//    }
//  return;
//}


void  findFOF(KdNode * root, fltarr pos, double dist, gadpart *** result, int * numpart, unsigned int * size)
{
  int minfac = 0;
  int *done = NULL;

  if ((*numpart == 0) || (*result == NULL))
    {
      int fac = 3;
      findGadparts(root, pos, fac *  dist, result, numpart, size);
      qsort(*result, *numpart, sizeof(gadpart *), cmp_pointer_id );
      minfac = 0.92 * fac;
    }

  int nold = 0;
  int j = 0;
  int i = 0;
  done = (int *) malloc(sizeof(int) * MAXINT);
  int ndone = 0;
  while (nold < *numpart)
    {
      nold = *numpart;
      for ( i = 0; i < *numpart; i++ )
	{
	  int* fnd = NULL;
	  int insertid = (*result)[i]->id;

	  fnd = bsearch( &(insertid), done, ndone, sizeof(int), cmp_int );
		   
	  if ( fnd == NULL )
	    {
	      double Dist=distance_nopb(pos, (*result)[i]->pos);
	      if (Dist > minfac * dist) 
		{
		  findNewGadparts(root, (*result)[i]->pos, dist, result, numpart, size);
		}
	      if (ndone < MAXINT)
		{
		  j = ndone-1;
		  while ((j>=0) && (done[j] > insertid))
		    {
		      done[j+1] = done[j];
		      j--;
		    }
		  j++;
		  done[j] = insertid;
		  ndone++;
		}		       		       
	    }
	}
    }
  free(done);
}


int checkTree(KdNode *root)
  {
    int check=0;
    if (root->up  !=NULL) check+=checkTree(root->up);
    if (root->down!=NULL) check+=checkTree(root->down);
    if (root->part->id==root->id) return (check);
    else return (++check);
  }


