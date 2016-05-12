#ifndef __USET_H_

#define __USET_H_
class USet
{
public:
	explicit USet(int size_)
	{
		parent = new int[size_];
		rank = new int[size_];
		setsize = new int[size_];
		kmers = new int[size_];
		Size = size_;
		for(int i=0;i<size_;++i)
			parent[i]=i;
		for(int i=0;i<size_;++i)
			rank[i]=0;
		for(int i=0;i<size_;++i)
			setsize[i]=1;
		for(int i=0;i<size_;++i)
			kmers[i]=0;
	}
	virtual ~USet()
	{
		delete []parent;
		delete []rank;
		delete []setsize;
		delete []kmers;
	}

	void ReInitial();
	int length()const{return Size;}
	int size()const{return Size;}
	int getsize(const int&x);
	void make_set(const int &x);
	int find(const int &x);
	void Union(int x,int y);
	void inckmers(int i);
	inline int getkmer(int k){return kmers[k];}

private:
	int* parent,*rank,*setsize,*kmers;
	int Size;
	USet(){}
	USet(const USet &uset){}
	const USet &operator=(const USet &uset){return *this;}
};
#endif
