/*
 * last fixed: 2012.01.05.
 * by Wang Yi.
 * */
#ifndef __TOM_ALGO_H_

#define __TOM_ALGO_H_

#include <string>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include "ULL4.h"

// kmer is a number represents a kmer like AGGT or TTAC, return the reverse complement.
int tomReverComple(int kmer,int kmerLen=4);
// Consider the reverse complements, return the number of kmers
inline int tomReverSize(int kmerLen=4);
// return the position of each kmer(rever complement is considered to be the same) in the old kmer-vector
int* getReverPosition(const int KmerLen);

double* tomNormalize(int *distri, int size,bool remove=true);
double* tomNormalize(double* distri, int size,bool newresult=false);

void printRead(const ULL4 &read,int blank,int Len);

/*
template <typename T>
void swap(T a,T b);
*/

template <typename T>
double* tomNormalize_rever(T *distri, int KmerLen,bool remove);

template <typename T>
double tomAbsDis(T* distri1,T* distri2,int size);
template <typename T>
double tomEuclidDis(T* distri1,T* distri2,int size);
template <typename T>
double tomRootDis(T* distri1,T* distri2,int size);
template <typename T>
double tom3RootDis(T* distri1,T* distri2,int size);
template <typename T>
double tom4RootDis(T* distri1,T* distri2,int size);
template <typename T>
void tomQuicksort(T* a,int lo,int hi);
template <typename T>
void tomQuicksort(T* distri,int* order,int lo,int hi);
template <typename T>
int* toSpear(T* distri,int size);
template <typename T>
int* toNewSpear(T* distriori,int size);
template <typename T>
int spearDistance(T* distri1,T* distri2,int size);
template <typename T>
double newspearDistance(T* distri1,T* distri2,int size);


//// Normalize the distribution and remove the duplicate caused by reverse complement. Assume the reverse complement has already been computed.
template <typename T>
double* tomNormalize_rever(T *distri,int* Position,int KmerLen = 4,bool remove = true)
{
	size_t orisize = 1<<(KmerLen<<1);
	double* result = NULL;
	size_t newsize = 0;
	if(KmerLen & 1)
		newsize = orisize>>1;
	else
		newsize = (orisize+(1<<KmerLen))>>1;
	result = new double[newsize];
	for(int i=0;i<newsize;++i)
	{
		result[i] = distri[Position[i]];
	}
	tomNormalize(result,newsize);
	if(remove)delete[]distri;
	return result;
}
//Define the Template functions
////calculate the absolute distance
template <typename T>
double tomAbsDis(T* distri1,T* distri2,int size)
{
	double sum = 0.;
	for(int i=0;i<size;++i)
	{
		double substr=(distri1[i]-distri2[i]);
		sum += (substr>=0?substr:-substr);
		//sum += (distri1[i]-distri2[i]);
		//std::cout<<distri1[i]<<'\t'<<distri2[i]<<'\t'<<(distri1[i]-distri2[i])<<std::endl;
	}
	return sum;
}
////calculate the Euclidean distance
template <typename T>
double tomEuclidDis(T* distri1,T* distri2,int size)
{
	double sum = 0.;
	for(int i=0;i<size;++i)
		sum += (distri1[i]-distri2[i])*(distri1[i]-distri2[i]);
	return sum;
}
////calculate the 2-root distance
template <typename T>
double tomRootDis(T* distri1,T* distri2,int size)
{
	double sum = 0.;
	for(int i=0;i<size;++i)
	{
		double substr=(distri1[i]-distri2[i]);
		sum += sqrt(substr>=0?substr:-substr);
	}
	return sum;
}
////calculate the Euclidean distance
template <typename T>
double tom3RootDis(T* distri1,T* distri2,int size)
{
	double sum = 0.;
	for(int i=0;i<size;++i)
	{
		double substr=(distri1[i]-distri2[i]);
		sum += pow((substr>=0?substr:-substr),1.0/3);
	}
	return sum;
}
////calculate the Euclidean distance
template <typename T>
double tom4RootDis(T* distri1,T* distri2,int size)
{
	double sum = 0.;
	for(int i=0;i<size;++i)
	{
		double substr=(distri1[i]-distri2[i]);
		sum += pow((substr>=0?substr:-substr),0.25);
	}
	return sum;
}
////standard quicksort
template <typename T>
void tomQuicksort(T* a, int lo, int hi)
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
{
	int i=lo, j=hi;
	T x=a[(lo+hi)/2],h;
	do
	{	
		while (a[i]<x) i++; 
		while (a[j]>x) j--;
		if (i<=j)
		{
			h=a[i]; a[i]=a[j]; a[j]=h;
			i++; j--;
		}
	} while (i<=j);
	if (lo<j) tomQuicksort(a, lo, j);
	if (i<hi) tomQuicksort(a, i, hi);
}

/////get the order for spear man with quick sort
template <typename T>
void tomQuicksort(T* distri,int* order,int lo, int hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
	int i=lo, j=hi,inttemp;
	T x=distri[(lo+hi)/2],h;
	do
	{	
		while (distri[i]<x) i++; 
		while (distri[j]>x) j--;
		if (i<=j)
		{
			h=distri[i]; distri[i]=distri[j]; distri[j]=h;
			inttemp=order[i];order[i]=order[j];order[j]=inttemp;
			i++; j--;
		}
	} while (i<=j);
	if (lo<j) tomQuicksort(distri,order, lo, j);
	if (i<hi) tomQuicksort(distri,order, i, hi);
}
////get the spearman order for a double array
template <typename T>
int* toSpear(T* distriori,int size)
{
	int* order = new int[size];
	T* distri = new T[size];
	for(int i=0;i<size;++i)
		order[i]=i;
	for(int i=0;i<size;++i)
		distri[i] = distriori[i];
	T h;int inttemp;
	for(int i=0;i<size;++i)
	{
		int small=i;
		for(int j=i+1;j<size;++j)
		{
			if(distri[small]>distri[j])
				small = j;
		}
		if(small!=i)
		{
			h=distri[i]; distri[i]=distri[small]; distri[small]=h;
			inttemp=order[i];order[i]=order[small];order[small]=inttemp;
		}
	}
 //   tomQuicksort(distriori, order, 0, size-1);
	int* result = new int[size];
	for(int i=0;i<size;++i)
		result[order[i]] = i;
	delete[] order;
	delete[] distri;
	return result;
}
////get the new spearman order for a double array
template <typename T>
int* toNewSpear(T* distriori,int size)
{
	int* order = new int[size];
	T* distri = new T[size];
	for(int i=0;i<size;++i)
		order[i]=i;
	for(int i=0;i<size;++i)
		distri[i] = distriori[i];
	T h;int inttemp;
	for(int i=0;i<size;++i)
	{
		int small=i;
		for(int j=i+1;j<size;++j)
		{
			if(distri[small]>distri[j])
				small = j;
		}
		if(small!=i)
		{
			h=distri[i]; distri[i]=distri[small]; distri[small]=h;
			inttemp=order[i];order[i]=order[small];order[small]=inttemp;
		}
	}
	//tomQuicksort(distri, order, 0, size-1);
	int* result = new int[size];
	int* order2 = new int[size];
	order2[0] = 0;
	for(int i=1;i<size;++i)
	{
		if(distri[i]==distri[i-1])
			order2[i]=order2[i-1];
		else order2[i]=i;
	}
	for(int i=0;i<size;++i)
		result[order[i]] = order2[i];
	delete[] order;
	delete[] distri;
	delete[]order2;
	return result;
}
////calculate the spearman distance
template <typename T>
int spearDistance(T* distri1,T* distri2,int size)
{
	int* spear1 = toSpear(distri1,size);
	int* spear2 = toSpear(distri2,size);
	int sum = 0;
	for(int i=0;i<size;++i)
	{
		sum += abs(spear1[i]-spear2[i]);
		//std::cout<<spear1[i]<<'\t'<<spear2[i]<<std::endl;
	}
	delete[]spear1;
	delete[]spear2;
	return sum;
}
////calculate the new spearman distance
template <typename T>
double newspearDistance(T* distri1,T* distri2,int size)
{
	int* spear1 = toNewSpear(distri1,size);
	int* spear2 = toNewSpear(distri2,size);
	double sum = 0;
	for(int i=0;i<size;++i)
	{
		sum += abs(spear1[i]-spear2[i]);
		//std::cout<<spear1[i]<<'\t'<<spear2[i]<<std::endl;
	}
	delete[]spear1;
	delete[]spear2;
	return sum;
}

inline int tomReverSize(int kmerLen)
{
	if(kmerLen & 1U)
		return ((1<<(kmerLen<<1))>>1);
	else
		return (((1<<(kmerLen<<1))+(1<<kmerLen))>>1);
}

/*
template <typename T>
void swap(T a,T b)
{
	T c = a ;
	a = b;
	b = c;
}*/
#endif
