#include <iostream>
#include <omp.h>
#include <algorithm>
#include <vector>
#include <map>
#include <set>

#include "MetaCluster.h"
#include "TomAlgorithm.h"
#include "USet.h"

using namespace std;

double MetaCluster::Distance(const int*  rank1,const int* rank2)
{
	double result = 0;
	for(int i=0;i<NewSize;++i)
		result += fabs(rank1[i]-rank2[i]);
	return result;
}

void MetaCluster::ComputeRank(int centerid)
{
	delete[] centerRank[centerid];
	centerRank[centerid] = toSpear(tomNormalize_rever(centerDistri[centerid],ReverPosition,KmerLen,false),NewSize);
}

double MetaCluster::DistanceAverage(vector<int> &kv1, vector<int> &kv2)
{
    vector<double> sum(kv1.size());
#pragma omp parallel for
    for (int i = 0; i < (int)kv1.size(); ++i)
    {
        sum[i] = 0;
        for (int j = 0; j < (int)kv2.size(); ++j)
        {
            sum[i] += Distance(SpearRank[kv1[i]], SpearRank[kv2[j]]);
        }
    }

    double total = 0;
    for (unsigned i = 0; i < sum.size(); ++i)
        total += sum[i];

	double count = kv1.size()*((double) kv2.size());
    if(count==0)return 0;
    else
    return total / count;
}

double MetaCluster::IntraDistance(vector<int> &kv1)
{
    vector<double> sum(kv1.size());
#pragma omp parallel for
    for (int i = 0; i < (int)kv1.size(); ++i)
    {
        sum[i] = 0;
        for (int j = i+1; j < (int)kv1.size(); ++j)
        {
            sum[i] += Distance(SpearRank[kv1[i]], SpearRank[kv1[j]]);
        }
    }

    double total = 0;
    for (unsigned i = 0; i < sum.size(); ++i)
        total += sum[i];

	double count = kv1.size();count = count*(kv1.size()-1)/2;
    if(count==0)return 0;
    else
    return total / count;
}

double MetaCluster:: DistributePoints()
{
	double power = 0;
	for (int i = 0; i < classes; ++i)
	{
		distSum[i] = 0;
		distSumSquare[i] = 0;
		distNum[i] = 0;
	}

#pragma omp parallel for
	for (int i = 0; i < Size; ++i)
	{
	//	if (isOutlier[i])continue;

		double minimum = 1e100;
		int t = -1;
		for (int j = 0; j < classes; ++j)
		{
			double d = Distance(SpearRank[i], centerRank[j]);
			if (d < minimum)
			{
				minimum = d;
				t = j;
			}
		}

		if (minimum < 1e100)
		{
			type[i] = t;
#pragma omp atomic
			power += minimum;
#pragma omp atomic
			distSum[t] += minimum;
#pragma omp atomic
			distSumSquare[t] += minimum * minimum;
#pragma omp atomic
			distNum[t]++;
		}
		else
			type[i] = -1;
	}

	for (int i = 0; i < classes; ++i)
	{
		if(distNum[i]>1)
		{
			distMean[i] = distSum[i] / distNum[i];
			distSD[i] = sqrt((distSumSquare[i] - distSum[i]*distMean[i]) / (distNum[i] - 1));
		}
		else
		{
			distMean[i] = distSum[i] / (distNum[i]+2);
			distSD[i] = sqrt((distSumSquare[i] - distSum[i]*distMean[i]) / (distNum[i] + 1));
		}
	}
/*
	int count = 0;
	for (int i = 0; i < Size; ++i)
	{
		if (isOutlier[i] != true)
			++count;
	}*/

    if(Size==0)return 0;
    else
	return power / Size;
}

bool MetaCluster::CalCenter()
{
	bool isUpdate = false;
	for (int j = 0; j < classes; ++j)
	{
		int sum[DisSize];
		for (int l = 0; l < DisSize; ++l)
			sum[l]=0;
		int count = 0;
		for (int i = 0; i < Size; ++i)
		{
		//	if (!isOutlier[i] && type[i] == j)
			if ( type[i] == j)
			{
				for (int l = 0; l < DisSize; ++l)
					sum[l] += KmerDistri[i][l];
				++count;
			}
		}

	if(count>0)
		for (int l = 0; l < DisSize; ++l)
			sum[l] /= count;

		for (int l = 0; l < DisSize; ++l)
		{
			if (fabs(centerDistri[j][l] - sum[l]) > 1e-6)
				isUpdate = true;
			centerDistri[j][l] = sum[l];
		}

		ComputeRank(j);
	}

	return isUpdate;
}

double MetaCluster::kmeanscluster(int classes_)
{
	// Randomly select the cluster center.
	classes = classes_;
	for (int i = 0; i < Size; ++i)
	{
		aux[i] = i;
	//	isOutlier[i] = false;
		type[i] = -1;
	}
	set<int>candidate;
	for (int i = 0; i < classes && Size-i>0; ++i)
	{
		int j = i + rand()%(Size-i);
		while(candidate.count(aux[j])==1)
			j = i + rand()%(Size-i);
		candidate.insert(aux[j]);
		int tmp = aux[i];aux[i]=aux[j];aux[j]=tmp;
		for(int k=0;k<DisSize;++k)
			centerDistri[i][k] = KmerDistri[aux[i]][k];
		ComputeRank(i);
	}
	double power = 0;
	int upcount = 0;
	for (; upcount < Times; ++upcount)
	{
		power = DistributePoints();
		if (CalCenter() == false)
			break;
	}

	power = DistributePoints();

	int count = 0;
	for (int i = 0; i < Size; ++i)
	{
	//	if (isOutlier[i])continue;

		double d = Distance(SpearRank[i], centerRank[type[i]]);
		if (d > distMean[type[i]] + 2 * distSD[type[i]])
		{
			++count;
		}
	}
	power = 0;
	bool updated = true;
	while(updated && upcount<=Times)
	{
			power = DistributePoints();
			updated = CalCenter();
			++upcount;
	//		if(upcount%100==0)cerr<<"iterations: "<<upcount<<endl;
	}
	cerr<<"All iterations: "<<upcount<<endl;
	return power;
}
double MetaCluster::muiltkmeans(const int ROUND,int classes_)
{
	double optimal = 1e100;
	int iterid = 0;
	for(int i=0;i<ROUND;++i)
	{
		double power = kmeanscluster(classes_);
		if(power < optimal)
		{
			optimal = power;
			copy(type,type+Size,best);
			iterid = i;
		}
	}
	copy(best,best+Size,type);
	cerr << "iteration "<<iterid<<" is the best."<<endl;
}

int MetaCluster::MergeClusters(double cret_shre)
{
	USet uset(classes);
	vector<int> clusters[classes];
	double interdistance[classes][classes];
//////////////////////////////////////////////////////////////////////
	//test
	int typemax = -1000,typemin = 9999;
	for (int i = 0; i < Size; ++i)
	{
		if(best[i]>typemax)typemax=best[i];
		if(best[i]<typemin)typemin=best[i];
	}
	//cerr<<"typemin: "<<typemin<<"\ttypemax: "<<typemax<<endl;
//////////////////////////////////////////////////////////////////////
	for (int i = 0; i < Size; ++i)
	{
	//	if (!isOutlier[i])
		clusters[best[i]].push_back(i);
	}

	for (int i = 0; i < classes; ++i)
	{
		for (int j = i+1; j < classes; ++j)
		{
			interdistance[i][j] = DistanceAverage(clusters[i], clusters[j]);
			interdistance[j][i] = interdistance[i][j];
		}
	}

	int current_classes = classes;
	int oriclasses = classes;

	int* clusize = new int[classes];
	int* cluSTsize = new int[classes];
	double* intradistance = new double[classes];
	for(int i= 0;i<classes;i++)
	{
		clusize[i]=0;
		cluSTsize[i]=0;
	}
	for(int i= 0;i<Size;i++)
	{
		clusize[best[i]]++;
		cluSTsize[best[i]] += SixTMerSize[i];
	}
	 
	for (int i = 0; i < classes; ++i)
	{
		intradistance[i] = IntraDistance(clusters[i]);
	}
	
	bool hasabovesh = true;
	int edgesize = classes*(classes-1)/2;
	////////////////////////////////////
	int** lastcomp = new int*[classes];
	for(int i=0;i<classes;++i)
	{
		lastcomp[i] = new int[GenoNum];
		for(int j=0;j<GenoNum;++j)
			lastcomp[i][j] = 0;
	}
	for(int i=0;i<Size;++i)
	{
		for(int j=0;j<GenoNum;++j)
			lastcomp[best[i]][j] += Component[i][j];
	}
	////////////////////////////////////////////////////////////
	for(int i=0;i<edgesize;++i)
	{
		if (current_classes <= 1)
			break;
		hasabovesh = false;	
		double disttmp = 1e100;
		int from,to;
		for (unsigned j1 = 0; j1 < classes; ++j1)
		{
			if(uset.find(j1)!=j1)continue;
			for (unsigned j2 = j1+1; j2 < classes; ++j2)
			{
				if(uset.find(j2)!=j2)continue;
				if(interdistance[j1][j2]<disttmp && (intradistance[j1]+intradistance[j2])/interdistance[j1][j2]>cret_shre*2)
				{
					disttmp = interdistance[j1][j2];
					from = j1;to=j2;
					hasabovesh = true;
				}
			}
		}
		if(hasabovesh==false)break;
		--current_classes;
		////////////////////////////////////////////////output info
		int major1=0,major2=0,sum1=lastcomp[from][0],sum2=lastcomp[to][0];
		for(int j=1;j<GenoNum;++j)
		{
			sum1 += lastcomp[from][j];
			sum2 += lastcomp[to][j];
			if(lastcomp[from][j]>lastcomp[from][major1])
				major1 = j;
			if(lastcomp[to][j]>lastcomp[to][major2])
				major2 = j;
		}
		/*
		if(major1==major2)cout<<"1,";
		else cout<<"0,";
		cout<<from<<","<<to<<","<<major1<<","<<sum1<<","<<major2<<","<<sum2<<","<<current_classes<<","<<interdistance[from][to]<<","<<intradistance[from]<<","<<intradistance[to]<<","<<clusize[from]<<","<<clusize[to]<<","<<cluSTsize[from]<<","<<cluSTsize[to]<<endl;
		*/
		////////////////////////////////update intradistance
		double totaldis = interdistance[from][to]*clusize[from]*clusize[to];
		totaldis += intradistance[from]*clusize[from]*(clusize[from]-1)/2;
		totaldis += intradistance[to]*clusize[to]*(clusize[to]-1)/2;
		totaldis /= ((double)(clusize[from]+clusize[to])*(double)(clusize[from]+clusize[to]-1)/2);
		intradistance[from]=intradistance[to]=totaldis;
		////////////////////////////////update interdistance
		for(int j=0;j<classes;++j)
		{
			double interdis=(interdistance[j][from]*clusize[from]+interdistance[j][to]*clusize[to])/(clusize[from]+clusize[to]);
			interdistance[from][j] = interdis;
			interdistance[j][from] = interdis;
			interdistance[to][j] = interdis;
			interdistance[j][to] = interdis;
		}
		////////////////////////////////update cluster size
		clusize[from] += clusize[to];
		clusize[to] = clusize[from];
		cluSTsize[from] += cluSTsize[to];
		cluSTsize[to] = cluSTsize[from];
		////////////////////////////////update lastcomp
		for(int j=0;j<GenoNum;++j)
		{
			lastcomp[from][j] += lastcomp[to][j];
			lastcomp[to][j]  = lastcomp[from][j];
		}
		////////////////////////////////////
		uset.Union(from,to);
	}
	//cout<<"metacluster finished!"<<endl;
	///////////////////////////////////
		////////////////////////////////////////////////output info
	//cout<<"Left interdistance and intradistance!"<<endl;
	for (unsigned j1 = 0; j1 < classes; ++j1)
	{
		if(uset.find(j1)!=j1)continue;
		for (unsigned j2 = j1+1; j2 < classes; ++j2)
		{
			if(uset.find(j2)!=j2)continue;
				int major1=0,major2=0,sum1=lastcomp[j1][0],sum2=lastcomp[j2][0];
				for(int j=1;j<GenoNum;++j)
				{
					sum1 += lastcomp[j1][j];
					sum2 += lastcomp[j2][j];
					if(lastcomp[j1][j]>lastcomp[j1][major1])
						major1 = j;
					if(lastcomp[j2][j]>lastcomp[j2][major2])
						major2 = j;
				}
		/*	if(major1==major2)cout<<"1,";
			else cout<<"0,";
			cout<<j1<<","<<j2<<","<<major1<<","<<major2<<","<<sum1<<","<<sum2<<",";
			cout<<Distance(centerRank[j1], centerRank[j2])<<",";
			cout<<interdistance[j1][j2]<<","<<intradistance[j1]<<","<<intradistance[j2]<<","<<clusize[j1]<<","<<clusize[j2]<<","<<cluSTsize[j1]<<","<<cluSTsize[j2]<<endl;
			*/
		}
	}
	//cout<<"end of left interdistance and intradistance!"<<endl;
	///////////////////////////////////////////
	map<int, int> id;
	int index = 0;
	for (int i = 0; i < classes; ++i)
	{
		if (uset.find(i)==i && id.find(i) == id.end())
			id[i] = index++;
	}

	for (int i = 0; i < Size; ++i)
	{
	//	if (!isOutlier[i])
		best[i] = id[uset.find(best[i])];
	}

	/////////////////////////////////////////////
	set<int> testbest;
	for (int i = 0; i < Size; ++i)
		testbest.insert(best[i]);
	/*cerr<<"best :"<<endl;
	for(set<int>::const_iterator itr = testbest.begin();itr!=testbest.end();++itr)
		cerr<< *itr<<"\t";
	cerr<<endl;*/
	//classes = id.size();
	cerr << "index:" <<index-1<<",\tclasses:"<<classes<<"\tcurr_class:"<<current_classes<< endl;
	///////////////////////////////////////////
	return current_classes;
}

int MetaCluster::iterMeta(const int ROUND,double threshold)
{
	//int FirstType[Size];
	muiltkmeans(ROUND,classes);
	int lastsize = MergeClusters(threshold);
	//copy(best,best+Size,FirstType);
	cerr<<"before and after 0: "<<classes<<",\t"<<lastsize<<endl;
	bool converge = false;
	while((lastsize <= classes*0.90 && classes >=10) || (classes < 10 && lastsize < classes-2))//|| classes-lastsize>=3)
	{
		converge = true;
		muiltkmeans(ROUND,lastsize);
		lastsize = MergeClusters(threshold);
        if(lastsize >= classes -1)break;
		cerr<<"before and after 1: "<<classes<<",\t"<<lastsize<<endl;
		
		///////////////////////////////////////////////////////////
		//for less-merge cases
		
		if(lastsize > classes*0.90)
		{
			int oriclasses = classes;
			muiltkmeans(ROUND,oriclasses+1);
			lastsize = MergeClusters(threshold);
			if(lastsize > classes*0.90)
			{
				muiltkmeans(ROUND,oriclasses-1);
				lastsize = MergeClusters(threshold);
			}
		}
		/////////////////////////////////////////////////////
	}
	if(converge)
	{
		if(classes+1 > MinSpecies)
			muiltkmeans(ROUND,classes+1);
		else
			muiltkmeans(ROUND,MinSpecies);
	//	lastsize = MergeClusters(threshold);//2011.11.30
	}
	cerr<<"before and after 2: "<<classes<<",\t"<<lastsize<<endl;
	///////////////////////////////////////////
	///////////////////////////////////////////
	return lastsize;
}
