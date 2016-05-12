/*
 * used as neighbor. we need to record the #. of k-mers shared for two cluster, which is a primer step for further merge.
 * by WangYi.(ywang@cs.hku.hk)
 * 2012-02-10
 */

#ifndef _TOM_SPARSEMATRIX_

#define _TOM_SPARSEMATRIX_
#include <iostream>
#include <algorithm>
#include <map>

#include <string.h>
typedef unsigned int uint;
#define MCInner_RadixSort_32_16_OneRound(arr, item, auxArr, length, shift) { \
	memset(count, 0, 65536 * sizeof(uint)); \
	for(uint i = 0; i < length; i++) \
		++count[(arr[i].item >> shift) & 0xFFFF]; \
	position[0] = 0; \
	for(uint i = 1; i < 65536; i++) \
		position[i] = position[i-1] + count[i-1]; \
	for(uint i = 0; i < length; i++) \
		auxArr[position[(arr[i].item >> shift) & 0xFFFF]++] = arr[i]; \
}
#define MC_RadixSort_32_16(arr, item, auxArr, length) { \
	uint count[65536]; \
	uint position[65536]; \
	MCInner_RadixSort_32_16_OneRound(arr, item, auxArr, length, 0) \
	MCInner_RadixSort_32_16_OneRound(auxArr, item, arr, length, 16) \
}

using namespace std;

struct SpMtrxNode
{
	unsigned x,y;
	unsigned data;
	SpMtrxNode(unsigned x_, unsigned y_)
	{
		x = x_;
		y = y_;
		data = 1;
	}
	bool operator<(const SpMtrxNode& t)const
	{
		return ((x<t.x)||((x==t.x)&&(y<t.y)));
	}
	void set(unsigned x_, unsigned y_, unsigned data_)
	{
		x = x_;
		y = y_;
		data = data_;
	}
	SpMtrxNode(){}
	const SpMtrxNode &operator=(const SpMtrxNode &t)
	{
		x = t.x;
		y = t.y;
		data = t.data;
		return *this;
	}
//private:
};
union TomVCTNode
{
	unsigned xy[2];
	SpMtrxNode* smvect;
};
struct TomVCT
{
	unsigned size,data;
	TomVCTNode vect;
	void increase(unsigned x_, unsigned y_)
	{
		if(size == 0)
		{
			size = 1;
			data = 1;
			vect.xy[0] = x_;
			vect.xy[1] = y_;
		}
		else if(size==1)
		{
			if(vect.xy[0]==x_ && vect.xy[1]==y_)
			{
				++data;
			}
			else
			{
				size = 2;
				unsigned x1 = vect.xy[0];
				unsigned y1 = vect.xy[1];
				vect.smvect = new SpMtrxNode[2];
				vect.smvect[0].set(x1,y1,data);
				vect.smvect[1].set(x_,y_,1);
			}
		}
		else
		{
			for(int i=0;i<size;++i)
			{
				if(vect.smvect[i].x==x_ && vect.smvect[i].y==y_)
				{
					++vect.smvect[i].data;
					return;
				}
			}
			if(size & (size-1))
			{
				vect.smvect[size].set(x_,y_,1);
				++size;
			}
			else
			{
				//////get new vector
				SpMtrxNode* node = vect.smvect;
				vect.smvect = new SpMtrxNode[size<<1];
				std::copy(node,node+size,vect.smvect);
				vect.smvect[size].set(x_,y_,1);
				++size;
				delete node;
			}
		}
	}
};

class SparseMatrix
{
public:
	static const unsigned Size = 1U<<20;
	unsigned NumNodes;
	TomVCT *HashTable;
	SpMtrxNode* SortedNegbr;

	SparseMatrix()
	{
		HashTable = new TomVCT[Size];
		for(int i=0;i<Size;++i)
			HashTable[i].size = 0;
	}
	virtual ~SparseMatrix()
	{
		delete[] HashTable;
	}

	void insert(unsigned x_, unsigned y_)
	{
		unsigned id = (x_+y_)%1048573;
		HashTable[id].increase(x_, y_);
	}
	int toNeighbor(map<int, map<int,int> >& Neighbor)
	{
		process();
		for(int i=0;i<NumNodes;++i)
		{
			Neighbor[SortedNegbr[i].x][SortedNegbr[i].y] = SortedNegbr[i].data;
		}
		clear();
		return NumNodes;
	}

	void clear()
	{
		if(HashTable != NULL)
			delete[] HashTable;
		if(SortedNegbr != NULL)
			delete[] SortedNegbr;
		HashTable = NULL;
		SortedNegbr = NULL;
	}

private:
	void process()
	{
		NumNodes = 0;
		for(int i=0;i<Size;++i)
			NumNodes += HashTable[i].size;
		std::cerr<<"NumNodes:\t"<<NumNodes<<std::endl;

		int Aidx = 0;
		SortedNegbr = new SpMtrxNode[NumNodes];
		for(int i=0;i<Size;++i)
		{
			if(HashTable[i].size == 1)
			{
				SortedNegbr[Aidx].set(HashTable[i].vect.xy[0],HashTable[i].vect.xy[1],HashTable[i].data);
				++Aidx;
			}
			else
			{
				for(int j=0;j<HashTable[i].size;++j)
				{
					SortedNegbr[Aidx] = HashTable[i].vect.smvect[j];
					++Aidx;
				}
			}
		}
		/////////////////////////rm HashTable
		for(int i=0;i<Size;++i)
		{
			if(HashTable[i].size>1)
			{
				delete[]HashTable[i].vect.smvect;
			}
		}
		delete[] HashTable;
		HashTable = NULL;

		SpMtrxNode* auxItems = new SpMtrxNode[NumNodes];
		MC_RadixSort_32_16(SortedNegbr, y, auxItems, NumNodes);
		MC_RadixSort_32_16(SortedNegbr, x, auxItems, NumNodes);
		delete[] auxItems;
	}
	SparseMatrix(const SparseMatrix &uset){}
	const SparseMatrix &operator=(const SparseMatrix &uset){return *this;}
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//To replace the map<int, map<int, int> >neighbor;
struct ydataNode
{
	unsigned y,data;
	ydataNode(){}
	ydataNode(unsigned y_, unsigned data_)
	{
		y = y_;
		data = data_;
	}
	void set(unsigned y_, unsigned data_)
	{
		y = y_;
		data = data_;
	}
};
struct xvect
{
	unsigned Size,x;
	ydataNode* yVect;
	xvect(){}
	xvect(unsigned size_, unsigned x_)
	{
		Size = size_;
		x = x_;
		yVect = new ydataNode[Size];
	}
};
struct LinkHash
{
	static const int Size = 1U<<20;
	unsigned* XSizes;
	unsigned* Xcount;
	xvect** XVects;
	LinkHash(SparseMatrix& sm)
	{
		XSizes = new unsigned[Size];
		Xcount = new unsigned[Size];
		XVects = new xvect*[Size];
		for(int i=0;i<Size;++i)
		{
			XSizes[i] = 0;
			Xcount[i] = 0;
			XVects[i] = NULL;
		}
		////////////////////////////test
		for(int i=1;i<sm.NumNodes;++i)
		{
			if(sm.SortedNegbr[i]<sm.SortedNegbr[i-1])
				std::cerr<<"neighbor not sorted!"<<std::endl;
		}
		///////////////////////////////
		unsigned preX = sm.SortedNegbr[0].x;
		XSizes[preX%1048573] = 1;
		for(int i=1;i<sm.NumNodes;++i)
		{
			unsigned curX = sm.SortedNegbr[i].x;
			if(curX != preX)
			{
				++XSizes[curX%1048573];
			}
			preX = curX;
		}
		for(int i=0;i<Size;++i)
		{
			if(XSizes[i]>0)
				XVects[i] = new xvect[XSizes[i]];
		}
		for(int i=0;i<sm.NumNodes;++i)
		{
			unsigned id = sm.SortedNegbr[i].x%1048573;
		}
	}
};
//////////////////////end.
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
#endif
