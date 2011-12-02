


#if 1

#include "TomoEngine/mt/mt19937ar.h"

#else
#error DO NOT USE THIS FILE. Use the Mersenne Twister code instead.

double random2();
int random3();
void srandom2(unsigned long num);
void readseed();
void writeseed();
double normal();
double dexprand();
#endif

