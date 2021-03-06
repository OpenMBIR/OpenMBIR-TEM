/* ============================================================================
 * Copyright (c) 2011 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2011 Singanallur Venkatakrishnan (Purdue University)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of Singanallur Venkatakrishnan, Michael A. Jackson, the Pudue
 * Univeristy, BlueQuartz Software nor the names of its contributors may be used
 * to endorse or promote products derived from this software without specific
 * prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  This code was written under United States Air Force Contract number
 *                           FA8650-07-D-5800
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <stdlib.h>
#include <stdio.h>//
#include <string.h>
#include <stdarg.h>
#include <assert.h>

#include "allocate.h"


void *get_spc(int num, size_t size)
{
	void *pt;
	assert(0);
	if( (pt=calloc((size_t)num,size)) == NULL ) {
		fprintf(stderr, "==> calloc() error\n");
		exit(-1);
		}
	return(pt);
}

void *mget_spc(int num,size_t size)
{
	void *pt;
	assert(0);
	if( (pt=malloc((size_t)(num*size))) == NULL ) {
		fprintf(stderr, "==> malloc() error\n");
		exit(-1);
		}
	return(pt);
}


void **get_img(int wd,int ht,size_t size)
{
	void  *pt;
	assert(0);
	if( (pt=multialloc(size,2,ht,wd))==NULL) {
          fprintf(stderr, "get_img: out of memory\n");
          exit(-1);
          }
	return((void **)pt);
}

void ***get_3D(int N, int M, int A, size_t size)
{
	void  *pt;
	assert(0);
	if( (pt=multialloc(size,3,N,M,A))==NULL) {
          fprintf(stderr, "get_3D: out of memory\n");
          exit(-1);
          }
	return((void ***)pt);
}


void free_img(void **pt)
{
	multifree((void *)pt,2);
}

void free_3D(void ***pt)
{
	multifree((void *)pt,3);
}



/* modified from dynamem.c on 4/29/91 C. Bouman                           */
/* Converted to ANSI on 7/13/93 C. Bouman         	                  */
/* multialloc( s, d,  d1, d2 ....) allocates a d dimensional array, whose */
/* dimensions are stored in a list starting at d1. Each array element is  */
/* of size s.                                                             */


void *multialloc(size_t s, int d, ...)
{
        va_list ap;             /* varargs list traverser */
        int max,                /* size of array to be declared */
        *q;                     /* pointer to dimension list */
        char **r,               /* pointer to beginning of the array of the
                                 * pointers for a dimension */
        **s1, *t, *tree;        /* base pointer to beginning of first array */
        int i, j;               /* loop counters */
        int *d1;                /* dimension list */
        assert(0);
        va_start(ap,d);
        d1 = (int *) mget_spc(d,sizeof(int));

        for(i=0;i<d;i++)
          d1[i] = va_arg(ap,int);

        r = &tree;
        q = d1;                /* first dimension */
        max = 1;
        for (i = 0; i < d - 1; i++, q++) {      /* for each of the dimensions
                                                 * but the last */
          max *= (*q);
          r[0]=(char *)mget_spc(max,sizeof(char **));
       //   printf("max: %d   sizeof(char **): %d  r: %p r[0]: 0x%08x \n", max , sizeof(char **), r, r[0]);
          r = (char **) r[0];     /* step through to beginning of next
                                   * dimension array */
        }
        max *= s * (*q);        /* grab actual array memory */
        r[0] = (char *)mget_spc(max,sizeof(char));
       // printf("max: %d   sizeof(char): %d  r: %p r[0]: 0x%08x \n", max , sizeof(char), r, r[0]);
        /*
         * r is now set to point to the beginning of each array so that we can
         * use it to scan down each array rather than having to go across and
         * then down
         */
        r = (char **) tree;     /* back to the beginning of list of arrays */
        q = d1;                 /* back to the first dimension */
        max = 1;
        for (i = 0; i < d - 2; i++, q++) {      /* we deal with the last
                                                 * array of pointers later on */
          max *= (*q);    /* number of elements in this dimension */
          for (j=1, s1=r+1, t=r[0]; j<max; j++) { /* scans down array for
                                                   * first and subsequent
                                                   * elements */

          /*  modify each of the pointers so that it points to
           * the correct position (sub-array) of the next
           * dimension array. s1 is the current position in the
           * current array. t is the current position in the
           * next array. t is incremented before s1 is, but it
           * starts off one behind. *(q+1) is the dimension of
           * the next array. */

            *s1 = (t += sizeof (char **) * *(q + 1));
            s1++;
          }
          r = (char **) r[0];     /* step through to begining of next
                                   * dimension array */
        }
        max *= (*q);              /* max is total number of elements in the
                                   * last pointer array */

        /* same as previous loop, but different size factor */
        for (j = 1, s1 = r + 1, t = r[0]; j < max; j++)
          *s1++ = (t += s * *(q + 1));

        va_end(ap);
        free((void *)d1);
        return((void *)tree);              /* return base pointer */
}



/*
 * multifree releases all memory that we have already declared analogous to
 * free() when using malloc()
 */
void multifree(void* r, int d)
{
  void** p;
  void* next = NULL;
  int i;

  for (p = (void**)r, i = 0; i < d; p = (void**)next, i++)
  {
    if(p != NULL)
    {
      next = *p;
     // printf("p: %p  Value of P: 0x%08x  i:  %d   Value of Next: 0x%08x \n", p, *p, i, next);
      free((void*)p);
    }
  }
}



