#ifndef __OMP_FUNCS
#define __OMP_FUNCS
#include <omp.h>

void myqsort(void* data, size_t num_elements, size_t element_size,
	     int (*comparer)(const void *, const void *));

int partarray(void* data, size_t num_elements, size_t element_size,
	      int (*comparer)(const void *, const void *), void **min, void **max);

#endif
