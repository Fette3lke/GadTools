#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "ompfuncs.h"

#define swapbytes(TYPE, i, j, n) {	     \
  register TYPE *a=(TYPE *)i;		     \
  register TYPE *b=(TYPE *)j;		     \
  long k=n;				     \
  do {					     \
    register TYPE t=*a;			     \
    *a++=*b;				     \
    *b++=t;				     \
  } while(--k>0);			     \
}


static __inline void swap(void* to, void* src, size_t element_size)
{
  //if (to!=src)
    {
      register char temp;
      register char *a=to;
      register char *b=src;
      int i;
      for (i=0; i<element_size; i++)
	{
	  temp=b[i];
	  b[i]=a[i];
	  a[i]=temp;
	}
    }
}

static __inline void bytecopy(void* to, void* src, size_t element_size)
{
  if (to!=src)
    {
      register char *a=to;
      register char *b=src;
      int i;
      for (i=0; i<element_size; i++)
	{
	  a[i]=b[i];
	}
    }
}

static __inline char *
med3(a, b, c, cmp)
	char *a, *b, *c;
	int (*cmp)(const void *, const void *);
{
	return cmp(a, b) < 0 ?
	       (cmp(b, c) < 0 ? b : (cmp(a, c) < 0 ? c : a ))
              :(cmp(b, c) > 0 ? b : (cmp(a, c) < 0 ? a : c ));
}


static __inline int partition(void* data, size_t num_elements, size_t element_size,
               int (*comparer)(const void *, const void *))
{
  int store=0, i;
  char* base=data;
  int pivotIndex=0;
  int d;
//  char *pm,*pl,*pn;
//  pm = (char *)data + (num_elements / 2) * element_size;
//  if (num_elements > 7) {
//    pl = data;
//    pn = (char *)data + (num_elements - 1) * element_size;
//    if (num_elements > 40) {
//      d = (num_elements / 8) * element_size;
//      pl = med3(pl, pl + d, pl + 2 * d, comparer);
//      pm = med3(pm - d, pm, pm + d, comparer);
//      pn = med3(pn - 2 * d, pn - d, pn, comparer);
//    }
//    pm = med3(pl, pm, pn, comparer);
//  }

  if (comparer(&base[0], &base[element_size*(num_elements-1)])==1)
    {
      if (comparer(&base[0], &base[element_size*((int)((num_elements-1)/2))])!=1)
	pivotIndex=0;
      else if (comparer(&base[element_size*(num_elements-1)], &base[element_size*((int)((num_elements-1)/2))])==1)
	pivotIndex=num_elements-1;
      else pivotIndex=(int)((num_elements-1)/2);
    } else 
    {
      if (comparer(&base[0], &base[element_size*((int)(num_elements-1)/2)])==1)
	pivotIndex=0;
      else if (comparer(&base[element_size*(num_elements-1)], &base[element_size*((int)((num_elements-1)/2))])!=1)
	pivotIndex=num_elements-1;
      else pivotIndex=(int)((num_elements-1)/2);
    }
  
  swapbytes(char , &base[element_size*(pivotIndex)], &base[element_size*(num_elements-1)], element_size);
  for (i=0; i<(num_elements-1); i++)
    {
      if (comparer(&base[element_size*i], &base[element_size*(num_elements-1)])!=1)
	{
	  if (i!=store) swapbytes(char ,&base[element_size*i], &base[element_size*store], element_size);
	  store++;
	}
    }
  swapbytes(char, &base[element_size*store], &base[element_size*(num_elements-1)], element_size);
  return store;
}

int partarray(void* data, size_t num_elements, size_t element_size,
		       int (*comparer)(const void *, const void *), void **min, void **max)
{
  int store=0, i;
  char* base=data;
  int pivotIndex=0;
  int d;
  char *pm,*pl,*pn;
  pm = (char *)data + (num_elements / 2) * element_size;
  if (num_elements > 7) {
    pl = data;
    pn = (char *)data + (num_elements - 1) * element_size;
    if (num_elements > 40) {
      d = (num_elements / 8) * element_size;
      pl = med3(pl, pl + d, pl + 2 * d, comparer);
      pm = med3(pm - d, pm, pm + d, comparer);
      pn = med3(pn - 2 * d, pn - d, pn, comparer);
    }
    pm = med3(pl, pm, pn, comparer);
  } 

//  if (comparer(&base[0], &base[element_size*(num_elements-1)])==1)
//    {
//      if (comparer(&base[0], &base[element_size*((int)((num_elements-1)/2))])!=1)
//	pivotIndex=0;
//      else if (comparer(&base[element_size*(num_elements-1)], &base[element_size*((int)((num_elements-1)/2))])==1)
//	pivotIndex=num_elements-1;
//      else pivotIndex=(int)((num_elements-1)/2);
//    } else 
//    {
//      if (comparer(&base[0], &base[element_size*((int)(num_elements-1)/2)])==1)
//	pivotIndex=0;
//      else if (comparer(&base[element_size*(num_elements-1)], &base[element_size*((int)((num_elements-1)/2))])!=1)
//	pivotIndex=num_elements-1;
//      else pivotIndex=(int)((num_elements-1)/2);
//    }
//  
//  swapbytes(char , &base[element_size*(pivotIndex)], &base[element_size*(num_elements-1)], element_size);
  *min=base;
  *max=base;
  swapbytes(char , pm, &base[element_size*(num_elements-1)], element_size);
  for (i=0; i<(num_elements-1); i++)
    {
      if (comparer(&base[element_size*i], &base[element_size*(num_elements-1)])!=1)
	{
	  if (i!=store) 
	    {
	      swapbytes(char ,&base[element_size*i], &base[element_size*store], element_size);
	    }
	  //if (comparer(*max, &base[element_size*store])==-1) *max=&base[element_size*store];
	  if (comparer(*min, &base[element_size*store])==1) *min=&base[element_size*store];
	  store++;
	} 
      //	if (comparer(*min, &base[element_size*i])==1) *min=&base[element_size*i];
      if (comparer(*max, &base[element_size*i])==-1) *max=&base[element_size*i];
      
    }
  if (store != i) swapbytes(char, &base[element_size*store], &base[element_size*(num_elements-1)], element_size);
  if (comparer(*max, &base[element_size*(num_elements-1)])==-1) *max=&base[element_size*(num_elements-1)];
  // else if (comparer(*min, &base[element_size*(num_elements-1)])==1) *min=&base[element_size*(num_elements-1)];
  return store;
}

static __inline void isort(void* data, size_t num_elements, size_t element_size,
               int (*comparer)(const void *, const void *))
{
  register int i, j=0;
  char *base=data;
  for (i=1; i< num_elements; i++)
    {
      j=i-1;
      while ((j>=0) && (comparer(&base[j*element_size], &base[(j+1)*element_size])==1))
	{
	  swapbytes(char ,&base[element_size*(j+1)],&base[element_size*j],element_size);
	  j--;
	}
    }
}

void myqsort(void* data, size_t num_elements, size_t element_size,
               int (*comparer)(const void *, const void *))
{
  if (num_elements<=1) return;
  if (num_elements<10) {isort(data, num_elements, element_size, comparer); return;}
  int index= partition(data, num_elements, element_size, comparer);
  char *base=data;
  if (index<(num_elements/2))
    {
#pragma omp parallel sections //firstprivate(base, p, q, r)
      {
#pragma omp section       
	myqsort(base, index, element_size, comparer);
#pragma omp section       
	myqsort(&base[element_size*(index+1)], num_elements-(index+1), element_size, comparer);
      }
    } else
    {
#pragma omp parallel sections //firstprivate(base, p, q, r)
      {
#pragma omp section 
	myqsort(&base[element_size*(index+1)], num_elements-(index+1), element_size, comparer);
#pragma omp section 
	myqsort(base, index, element_size, comparer);
      }
    }
}

void print_omp_info(void)
{


}
