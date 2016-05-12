#ifndef __MetaCluster_H_

#define __MetaCluster_H_
#include <iostream>
#include <string>
#include <vector>
#include "TomAlgorithm.h"

using namespace std;
class MetaCluster
{
public:
	explicit MetaCluster(int KmerLen_,int size_,int**spear_,int**KmerDistri_,int*SixTmer_,int GenoNum_,int** Component_,int kmeansize,int MaxSpecies_,int MinSpecies_)
	{
		KmerLen = KmerLen_;
		Size = size_;
		SpearRank = spear_;
		SixTMerSize = SixTmer_;
		KmerDistri = KmerDistri_;
		DisSize = 1<<(KmerLen<<1);

		MaxSpecies = MaxSpecies_;
		MinSpecies = MinSpecies_;

		GenoNum = GenoNum_;
		Component = new int*[Size];
		for(int i=0;i<Size;++i)
		{
			Component[i] = new int[GenoNum];
			for(int j=0;j<GenoNum;++j)
				Component[i][j] = Component_[i][j];
		}

		int sixall=0;for(int i=0;i<Size;++i)sixall+=SixTMerSize[i];
		//classes = Size/100;
		if(Size < MaxSpecies)
		{
			MaxSpecies = Size/2;
			cerr << "MaxSpecies is too large. We will start with half of the group number." << endl;
		}
			
		if(kmeansize==0)
			classes = MaxSpecies;//sixall/600000+1;
		else classes = kmeansize;

		cerr<<"Size:"<<Size<<",\t"<<sixall<<endl;
		cerr<<"classes:"<<classes<<",\t"<<sixall/1200000<<endl;
		cerr<<"KmerDistri last:"<<KmerDistri[0][0]<<"\t"<<KmerDistri[Size-1][255]<<",\t"<<SixTMerSize[Size-1]<<endl;
		distSum = new int[classes];
		distNum = new int[classes];
		distSumSquare = new int[classes];

		distMean = new int[classes];
		distSD = new int[classes];
//////////////////////////////////////////////////////
	//	isOutlier = new bool[Size];
		type = new int[Size];
		best = new int[Size];
		aux = new int[Size];
//////////////////////////////////////////////////////
		ReverPosition = getReverPosition(KmerLen);
		NewSize = tomReverSize(KmerLen);
//////////////////////////////////////////////////////
		centerDistri = new int*[classes];
		for(int i=0;i<classes;++i)
		{
			centerDistri[i] = new int[DisSize];
		}
		centerRank = new int*[classes];
		for(int i=0;i<classes;++i)
		{
			centerRank[i] = new int[NewSize];
		}
	}
	virtual ~MetaCluster(){}
////////////////////////////////////////////////////////
	double muiltkmeans(const int ROUND,int classes_);
	int MergeClusters(double shre_t);
	int iterMeta(const int ROUND,double shre_t);
///////////////////////////////////////////////////////
	int Size;
	int classes;
	int **SpearRank;
	int **KmerDistri;
	int *SixTMerSize;
///////////////////////////////////////////////////////
	int *type;
	int *best;
private:
	//static const int KmerLen = 4;
	//static const int DisSize = 1<<(KmerLen<<1);
	int KmerLen;
	int DisSize;
	int NewSize;

	int MaxSpecies;
	int MinSpecies;

	int *distSum;
	int *distNum;
	int *distSumSquare;

	int *distMean;
	int *distSD;

	int **centerRank;
	int **centerDistri;
//////////////////////////////////////////////////////
	int* ReverPosition;
//////////////////////////////////////////////////////
	//bool *isOutlier;
	int *aux;
//////////////////////////////////////////////////////
	int GenoNum;
	int** Component;
//////////////////////////////////////////////////////
	static const int Times=200;
//////////////////////////////////////////////////////
	void ComputeRank(int id);
	double Distance(const int* rank1,const int* rank2);
////////////////////////////////////////////////////////
	double DistributePoints();
	bool CalCenter();
	double kmeanscluster(int classes_);
	double DistanceAverage(vector<int> &kv1,vector<int> &kv2);
	double IntraDistance(vector<int> &kv1);

//////////////////////////////////////////////////////
	MetaCluster(){}
	/*MetaCluster(const MetaCluster &kmeans){}
	const MetaCluster &operator=(const MetaCluster &kmeans){return *this;}*/
};
#endif
