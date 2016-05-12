#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "TomAlgorithm.h"

using namespace std;

void printRead(const ULL4 &read,int blank,int Len=75)
{
	for(int i=0;i<blank;++i)
		cout << ' ';
	for(int i=Len-1;i>=0;--i)
	{
		cout << ((read >> (i<<1)).low&3U);
	}
	cout << endl;
}

double* tomNormalize(double* distri, int size,bool newresult)
{
	if(newresult)
	{
		double sum=0.;
		for(int i=0;i<size;++i)
			sum += distri[i];
		if(sum==0)return NULL;

		double *result = new double[size];
		for(int i=0;i<size;++i)
			result[i] = distri[i]/sum;
		return result;
	}
	else
	{
		double sum=0.;
		for(int i=0;i<size;++i)
			sum += distri[i];
		if(sum==0)return NULL;
	
		for(int i=0;i<size;++i)
			distri[i] /= sum;
		return NULL;
	}
}

double* tomNormalize(int *distri, int size,bool remove)
{
	double sum=0.;
	for(int i=0;i<size;++i)
		sum += distri[i];
	if(sum==0)return NULL;

	double *result = new double[size];
	for(int i=0;i<size;++i)
		result[i] = distri[i]/sum;
	if(remove)
		delete[]distri;
	return result;
}
/*
inline int tomReverSize(int kmerLen)
{
	if(kmerLen & 1U)
		return ((1<<(kmerLen<<1))>>1);
	else
		return (((1<<(kmerLen<<1))+(1<<kmerLen))>>1);
}*/

//kmer is a number represents a kmer like AGGT or TTAC, return the reverse complement.
int tomReverComple(int kmer,int kmerLen)
{
	int result=0;
	unsigned orig=kmer;
	for(int i=0;i<kmerLen;++i)
	{
		unsigned digit=orig&3U;
		orig >>= 2;
		result <<= 2;
		result |= (~digit)&3U;
	}
	return result;
}

//return the position of each kmer(rever complement is considered to be the same) in the old kmer-vector
int* getReverPosition(const int KmerLen)
{
	size_t orisize = 1<<(KmerLen<<1);
	size_t newsize;
	if(KmerLen & 1)
		newsize = orisize>>1;
	else
		newsize = (orisize+(1<<KmerLen))>>1;
	int *result = new int[newsize];
	int count = 0;
	for(int i=0;i<orisize;++i)
	{
		if(tomReverComple(i,KmerLen)<i)continue;
		result[count++]=i;
	}
	return result;
}
