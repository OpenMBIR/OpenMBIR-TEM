#include <stdlib.h>

void *get_spc(int num, size_t size);
void *mget_spc(int num, size_t size);
void **get_img(int wd,int ht, size_t size);
void ***get_3D(int N, int M, int A, size_t size);
void free_img(void **pt);
void free_3D(void ***pt);
void *multialloc(size_t s, int d, ...);
void multifree(void *r,int d);



