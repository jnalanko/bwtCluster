/*
 * last fixed: 2012.01.05.
 * by Wang Yi.
 * */
#ifndef __TOM_BINNER_H_

#define __TOM_BINNER_H_

#include <string>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <omp.h>
#include "ULL4.h"
#include "USet.h"
#include "Structs2.h"

////////////////////////////////////////////////////////////////////////
//functions for KmerMap
void iniMap();
void countKmer(const int ReadNum, const int ReadLen, const ULL4* Reads, const ULL4* RevReads);
void iniArrayMap();
void mapReads(const int ReadNum, const int ReadLen, const ULL4* Reads, const ULL4* RevReads,bool* IsOrphan);
void testVector(const vector<pairlong>& ITR);
void testVCount();

void clearSameTypeNodes(const vector<pairlong>& ITR, USet* Type);
void rmSameTypeNodes(USet* Type);

void setPara(const vector<pairlong> &ITR, vector<int> &ITRPVEC,const int paranum);
//end 
////////////////////////////////////////////////////////////////////////

//is two ULL4 have exactly one base(2bits)difference;
bool isclosekmer(const ULL4 &a,const ULL4 &b);

//given two ULL4 that has one base differ, return the position and base of the different base
int kmerdiffer(const ULL4 &a,const ULL4 &b);
unsigned diffsearch(unsigned long long a);

string ull2str(unsigned long long a, int size);//turn the lower (size*2) bits into DNA string
string read2str(const ULL4& read,int ReadLen);

int ullDist(unsigned long long a, unsigned long long b);

ULL4 array2read(int a[],int readLen);
ULL4 array2rcread(int a[],int readLen);

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
void compa_read(int ReadLen, const int &position1,const int &position2,int &match,int &mismatch,const ULL4 &read1,const ULL4 &read2);

void USMerge(const vector<pairlong>& ITR, USet* Type, int LARGESIZE, int ReadLen,
	   	ULL4*  Reads, ULL4*  RevReads, int StrLen, int StrLenLowCover,bool* IsOrphan);
#endif
