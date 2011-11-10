#ifndef _allocate_h_
#define _allocate_h_



#include "EIMTomo/EIMTomo.h"

#ifdef EIMTOMO_HAVE_STDLIB_H
#include <stdlib.h>
#endif

#if EIMTOMO_HAVE_STDINT_H
#include <stdint.h>
#endif




#ifdef __cplusplus
extern "C"
{
#endif

  void *get_spc(int num, size_t size);
  void *mget_spc(int num, size_t size);
  void **get_img(int wd, int ht, size_t size);
  void free_img(void **pt);

  void ***get_3D(int wd, int ht, int dp, size_t size);
  void free_3D(void ***pt);

  void *multialloc(size_t s, int d, ...);
  void multifree(void *r,int d);

#ifdef __cplusplus
}
#endif



#endif /* _allocate_h_  */
