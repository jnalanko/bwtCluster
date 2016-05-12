/*
 * last fixed: 2012.02.04.
 * by Wang Yi.
 * */
#ifndef __TOM_STRUCTS_H_

#define __TOM_STRUCTS_H_
#include <iostream>
#include <algorithm>
#include <omp.h>
typedef unsigned INodeRef;
//typedef unsigned long long KmerType;
typedef unsigned  KmerType;
const unsigned INULL = ~0U;

//const int TRANLEN = 25;
//const unsigned long long TRANMASK = 0x3ffffffffffff;
//const unsigned HASHMASK = 0x3fffffff;
//const unsigned HASHSIZE = 1U<<30;

const int TRANLEN = 16;
const unsigned TRANMASK = 0xffffffff;
const unsigned HASHMASK = 0x3fffffff;
const unsigned HASHSIZE = 1U<<30;

////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////structure & function
#pragma pack(1)
struct NodeIDPosi
{
	unsigned id;
	unsigned char posi;
	bool operator<(const NodeIDPosi& t)const
	{
		return posi<t.posi;
	}
};

#pragma pack()
struct IntNode
{
	unsigned VSize;
	INodeRef next;
	KmerType kmer;
	NodeIDPosi* myvector;
	IntNode(KmerType kmer_,unsigned VSize_,INodeRef next_)
	{
		kmer = kmer_;
		VSize = VSize_;
		next = next_;
		myvector = NULL;
	}
	IntNode()
	{
		VSize = 0;
		myvector = NULL;
		next = INULL;
	}
	void push_back(unsigned ele,unsigned posi)
	{
		myvector[VSize].id=ele;
		myvector[VSize].posi=posi;
		++VSize;
	}
	void set(KmerType kmer_,unsigned VSize_)
	{
		kmer = kmer_;
		VSize = VSize_;
		myvector = NULL;
		next = INULL;
	}
	void newVector()
	{
		myvector = new NodeIDPosi[VSize];
	}

	void clear()
	{
		VSize = 0;
		if(myvector!=NULL)
		{
			delete[]myvector;
			myvector = NULL;
		}
	}
	void cpfrom(IntNode& node)
	{
		myvector = node.myvector;
		next = node.next;
		kmer = node.kmer;
		VSize = node.VSize;

		node.VSize = 0;
		node.myvector = NULL;
	}
private:
	IntNode(const IntNode &node){}
	const IntNode &operator=(const IntNode &node){return *this;}
};


//have to process memory myself.
class IntNodeAloc
{
public:
	const static int RowSize = 1U<<20;
	const static int RowNum = 1U<<12;
	const static int RowSizeBit = 20;
	const static unsigned MASK = 0xfffff;
	IntNode** AllNodes;
	omp_lock_t getnew_lock;

	IntNodeAloc()
	{
		NodeNum = 0;
		AllNodes = new IntNode*[RowNum];
		for(unsigned i=0;i<RowNum;++i)
		{
			AllNodes[i] = NULL;
		}
		omp_init_lock(&getnew_lock);
	}
	void clear()
	{
		for(unsigned i=0;i<NodeNum;++i)
		{
			AllNodes[i>>RowSizeBit][i&MASK].clear();
		}
		for(unsigned i=0;i<RowNum;++i)
		{
			if(AllNodes[i] != NULL)
				delete[] AllNodes[i];
			AllNodes[i] = NULL;
		}
		delete[] AllNodes;
		NodeNum = 0;
		AllNodes = NULL;
	}
	void reIni()
	{
		/*
		for(unsigned i=0;i<NodeNum;++i)
			AllNodes[i>>RowSizeBit][i&MASK].clear();
		*/
		for(unsigned i=0;i<RowNum;++i)
		{
			if(AllNodes[i] != NULL)
				delete[] AllNodes[i];
			AllNodes[i] = NULL;
		}
		NodeNum = 0;
	}
	/*
	virtual ~IntNodeAloc()
	{
		for(int i=0;i<NodeNum;++i)
			clrVect(i);
		for(int i=0;i<=(NodeNum-1)/RowSize;++i)
		{
			delete[] AllNodes[i];
		}
		delete[]AllNodes;
		AllNodes = NULL;
		omp_destroy_lock(&getnew_lock);
	}
	*/

	unsigned getNew(KmerType kmer_,unsigned VSize_)
	{
		omp_set_lock(&getnew_lock);
		unsigned rowNo = NodeNum>>RowSizeBit;
		if(AllNodes[rowNo]==NULL)
			AllNodes[rowNo] = new IntNode[RowSize];
		unsigned curNum = NodeNum++;
		omp_unset_lock(&getnew_lock);
		AllNodes[rowNo][curNum&MASK].set(kmer_,VSize_);
		return curNum;
	}

	//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////
	unsigned getNodeNum() const
	{
		return NodeNum;
	}
	unsigned getVSize(unsigned id)
	{
		return AllNodes[id>>RowSizeBit][id&MASK].VSize;
	}
	unsigned getNext(unsigned id)
	{
		return AllNodes[id>>RowSizeBit][id&MASK].next;
	}
	KmerType getKmer(unsigned id)
	{
		return AllNodes[id>>RowSizeBit][id&MASK].kmer;
	}
	NodeIDPosi* getVector(unsigned id)
	{
		return AllNodes[id>>RowSizeBit][id&MASK].myvector;
	}

	void setVSize(unsigned id,unsigned VSize_)
	{
		AllNodes[id>>RowSizeBit][id&MASK].VSize = VSize_;
	}
	void incVSize(unsigned id)
	{
		++AllNodes[id>>RowSizeBit][id&MASK].VSize;
	}
	unsigned setNext(unsigned id,unsigned next_)
	{
		AllNodes[id>>RowSizeBit][id&MASK].next = next_;
	}
	unsigned setKmer(unsigned id,KmerType kmer_)
	{
		AllNodes[id>>RowSizeBit][id&MASK].kmer = kmer_;
	}
	void newVector(unsigned id)
	{
		AllNodes[id>>RowSizeBit][id&MASK].newVector();
	}
	IntNode* getRef(unsigned id)
	{
		return &AllNodes[id>>RowSizeBit][id&MASK];
	}
	void push_back(unsigned id,unsigned ele,unsigned posi)
	{
		AllNodes[id>>RowSizeBit][id&MASK].push_back(ele,posi);
	}
	void clrVect(unsigned id)
	{
//		delete[] AllNodes[id>>RowSizeBit][id&MASK].myvector;
		AllNodes[id>>RowSizeBit][id&MASK].clear();
	}

	inline INodeRef updateNext(INodeRef* KmerMap,const INodeRef next_,unsigned hashid,unsigned nextnew)
	{
		if(KmerMap[hashid] == next_)
		{
			KmerMap[hashid] = nextnew;
			return hashid;
		}
		for(INodeRef curr = KmerMap[hashid];curr != INULL;curr = getNext(curr))
		{
			if(getNext(curr)==next_)
			{
				AllNodes[curr>>RowSizeBit][curr&MASK].next = nextnew;
				return curr;
			}
		}
		return INULL;
	}
	void shrinkSize1(INodeRef* KmerMap)
	{
		unsigned srcid = NodeNum-1, tgtid = 0;//move src node to tgt(target) node
		while(tgtid < srcid && AllNodes[tgtid>>RowSizeBit][tgtid&MASK].VSize >= 2)
			++tgtid;
		while(tgtid < srcid && AllNodes[srcid>>RowSizeBit][srcid&MASK].VSize < 2)
			--srcid;
		while(srcid > tgtid)
		{
		//	unsigned testhashid = getKmer(srcid)&HASHMASK;
			/*
			for(INodeRef curr = KmerMap[testhashid];curr != INULL;curr = getNext(curr))
			{
				if(getVSize(curr) < 2)
					std::cerr<<std::hex<<"in shrinksize 1:\t"<< testhashid<<"\t"<<srcid<<"\t"<<tgtid<<std::dec<<std::endl;
				if(testhashid==0x9e94951)
					std::cerr<<"9e94951:1\t"<<KmerMap[testhashid]<<"\t"<<curr<<"\t"<<getKmer(curr)<<"\t"<<getVSize(curr)<<"\t"<<getNext(curr)<<std::endl;
			}
			*/
			if(KmerMap != NULL)
			{
				unsigned result = updateNext(KmerMap, srcid, getKmer(srcid)&HASHMASK, tgtid);
				if(result==INULL)
				{
					std::cerr<<"ERROR in shrinksize 1:\t"<<tgtid<<"\t"<<std::hex<<getKmer(srcid)<<std::dec<<"\t"<<srcid<<std::endl;
				}
			}
			/*
			if(testhashid==0x9e94951)
				std::cerr<<"9e94951:2\t"<<KmerMap[testhashid]<<"\t"<<getVSize(KmerMap[testhashid])<<std::endl;
			if(result==INULL)
			{
				std::cerr<<tgtid<<"\t"<<std::hex<<getKmer(srcid)<<std::dec<<"\t"<<srcid<<std::endl;
				for(INodeRef curr = KmerMap[getKmer(srcid)&0xfffffff];curr != INULL;curr = getNext(curr))
					std::cerr<<"nodes:\t"<<std::hex<<getKmer(curr)<<"\t"<<(getKmer(curr)&0xfffffff)<<std::endl;
				std::cerr<<std::dec<<"updateNext result is INULL"<<std::endl;
			}
			if(testhashid==0x9e94951)
			{
				std::cerr<<"9e94951:3\t"<<getVSize(KmerMap[testhashid])<<std::endl;
				std::cerr<<"hashid:\t"<<tgtid<<"\t"<<srcid<<"\t"<<KmerMap[testhashid]<<std::endl;
			}*/
			AllNodes[tgtid>>RowSizeBit][tgtid&MASK].cpfrom(AllNodes[srcid>>RowSizeBit][srcid&MASK]);
			/*

			testhashid = getKmer(tgtid)&0xfffffff;
			if(testhashid==0x9e94951)
			{
				std::cerr<<"9e94951:4\t"<<getVSize(KmerMap[testhashid])<<"\t"<<getVSize(tgtid)<<std::endl;
				std::cerr<<"hashid:\t"<<tgtid<<"\t"<<srcid<<"\t"<<KmerMap[testhashid]<<std::endl;
			}
			for(INodeRef curr = KmerMap[testhashid];curr != INULL;curr = getNext(curr))
			{
				if(getVSize(curr) < 2)
				{
					std::cerr<<std::hex<<"in shrinksize 2:\t"<< testhashid<<"\t"<<getKmer(curr)<<"\t"<<std::dec<<srcid<<"\t"<<tgtid<<std::dec<<std::endl;
				}
				if(testhashid==0x9e94951)
					std::cerr<<"9e94951:\t"<<curr<<"\t"<<getKmer(curr)<<"\t"<<getVSize(curr)<<"\t"<<getNext(curr)<<std::endl;
			}*/
			--srcid;
			++tgtid;
			while(tgtid < srcid && AllNodes[tgtid>>RowSizeBit][tgtid&MASK].VSize >= 2)
				++tgtid;
			while(tgtid < srcid && AllNodes[srcid>>RowSizeBit][srcid&MASK].VSize < 2)
				--srcid;
		}
		const unsigned oldnum = (NodeNum-1)>>RowSizeBit;
		NodeNum=0;
		while(getVSize(NodeNum)>=2)
			++NodeNum;
		for(unsigned i=((NodeNum-1)>>RowSizeBit)+1;i<=oldnum;++i)
		{
			delete[]AllNodes[i];
			AllNodes[i] = NULL;
		}
		std::cout<<"deleted rows:\t"<<oldnum-((NodeNum-1)>>RowSizeBit)<<std::endl;
	}
private:
	unsigned NodeNum;
	IntNodeAloc(const IntNodeAloc &uset){}
	const IntNodeAloc &operator=(const IntNodeAloc &uset){return *this;}
};//INodePool;
////////////////////////////////////////////////////
////////for same NodeAloc
struct SameTNode
{
	unsigned count,type;
	KmerType kmer;
	SameTNode(){}
	void set(KmerType kmer_,unsigned count_,unsigned type_)
	{
		kmer = kmer_;
		count = count_;
		type = type_;
	}
	void getData(KmerType &kmer_, unsigned &count_, unsigned &type_)
	{
		kmer_ = kmer;
		count_ = count;
		type_ = type;
	}
};
class SameTNAloc
{
public:
	const static int RowSize = 1U<<20;
	const static int RowNum = 1U<<12;
	const static int RowSizeBit = 20;
	const static unsigned MASK = 0xfffff;
	SameTNode** NodeArray;
	omp_lock_t getnew_lock;

	SameTNAloc()
	{
		NodeNum = 0;
		NodeArray = new SameTNode*[RowNum];
		for(unsigned i=0;i<RowNum;++i)
		{
			NodeArray[i] = NULL;
		}
		omp_init_lock(&getnew_lock);
	}
	/*
	virtual ~SameTNAloc()
	{
		clear();
		omp_destroy_lock(&getnew_lock);
	}
	*/

	unsigned getNew(KmerType kmer_, unsigned count_, unsigned type_)
	{
		omp_set_lock(&getnew_lock);
		unsigned rowNo = NodeNum>>RowSizeBit;
		if(NodeArray[rowNo]==NULL)
			NodeArray[rowNo] = new SameTNode[RowSize];
		unsigned curNum = NodeNum++;
		omp_unset_lock(&getnew_lock);
		NodeArray[rowNo][curNum&MASK].set(kmer_, count_, type_);
		return curNum;
	}
	unsigned getNodeNum()
	{
		return NodeNum;
	}
	/*
	unsigned getTest()
	{
		return test;
	}
	*/

	unsigned getType(unsigned id)
	{
		return NodeArray[id>>RowSizeBit][id&MASK].type;
	}
	unsigned getCount(unsigned id)
	{
		return NodeArray[id>>RowSizeBit][id&MASK].count;
	}
	KmerType getKmer(unsigned id)
	{
		return NodeArray[id>>RowSizeBit][id&MASK].kmer;
	}
	void getData(unsigned id,KmerType &kmer_, unsigned &count_, unsigned &type_)
	{
		NodeArray[id>>RowSizeBit][id&MASK].getData(kmer_,count_,type_);
	}
	void clear()
	{
		for(unsigned i=0;i<RowNum;++i)
		{
			if(NodeArray[i]!=NULL)
				delete[]NodeArray[i];
			NodeArray[i] = NULL;
		}
		delete[]NodeArray;
	}
private:
	unsigned NodeNum;
//	unsigned test;
	SameTNAloc(const SameTNAloc &uset){}
	const SameTNAloc &operator=(const SameTNAloc &uset){return *this;}
};
#endif
