/*
 * last fixed: 2012.02.04.
 * by Wang Yi.
 * */
#ifndef __TOM_STRUCTS2_H_

#define __TOM_STRUCTS2_H_
#include <iostream>
#include <algorithm>
#include "global.h"
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////
//structs
struct pairlong
{
	//INodeRef node;
	INodeRef node;
	pairlong(){node=INULL;}
	pairlong(INodeRef node_)
	{
		node = node_;
	}
	bool operator<(const pairlong& t)const
	{
	//	return ((node->VSize)<(t.node->VSize));
		return ((INodePool.getVSize(node))<(INodePool.getVSize(t.node)));
	}
	bool operator==(const pairlong& t)const
	{
	//	return (node->VSize==t.node->VSize);
		return ((INodePool.getVSize(node))==(INodePool.getVSize(t.node)));
	}
};
#endif
