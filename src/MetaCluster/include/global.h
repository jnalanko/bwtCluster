#ifndef __GLOBAL_H_

#define __GLOBAL_H_
#include "ULL4.h"
#include "Structs.h"
#include <omp.h>

const long long MAXLINE = 500;
const long long MAXLENGTH = 14000000;
extern ULL4 AndMask[129];//get last k bases i.e. 2k bits.
extern double confidence[129];

typedef unsigned NodeRef;
extern INodeRef* KmerMap;

extern int NodeCount;
extern IntNodeAloc INodePool;
extern SameTNAloc SNodePool;
/*
extern omp_lock_t Ipool_lock;
extern omp_lock_t Ipool_inc_lock;
*/
#endif
