#include <iostream>
#include "USet.h"
void USet::inckmers( int i)
{
	++kmers[find(i)];
}
void USet::make_set(const int &x)
{
	parent[x]=x;
	setsize[x]=1;
	rank[x]=0;
}
int USet::getsize(const int &x)
{
	int y=find(x);
	return setsize[y];
}
int USet::find(const int &x)
{
	//std::cerr<<"Finding "<<x<<std::endl;
	int y,temp,root;
	root=x;
	while(parent[root]!=root)
		root=parent[root];
	y=x;
	while(parent[y]!=y)//compress path
	{
		temp=parent[y];
		parent[y]=root;
		y=temp;
	}
	return root;
}

void USet::Union(int x,int y)
{
	int u,v;
	//std::cerr<<"Union "<<x<<","<<y<<std::endl;
	u=find(x);v=find(y);
	if(u==v)return;
	//std::cerr<<"Unioning... "<<x<<':'<<u<<","<<y<<':'<<v<<std::endl;
	if(rank[u]<=rank[v])//Union operation
	{
		parent[u]=v;
		setsize[v] += setsize[u];
		kmers[v] += kmers[u];
		if(rank[u]==rank[v])
			rank[v]++;
	}
	else
	{
		parent[v]=u;
		setsize[u] += setsize[v];
		kmers[u] += kmers[v];
	}
}

void USet::ReInitial()
{
		for(int i=0;i<Size;++i)
			parent[i]=i;
		for(int i=0;i<Size;++i)
			rank[i]=0;
		for(int i=0;i<Size;++i)
			setsize[i]=1;
		for(int i=0;i<Size;++i)
			kmers[i]=0;
}
