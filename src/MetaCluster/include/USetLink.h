/*
 * Record the elements of each set with linked list.
 * NOTICE: When you iterate the links,for each type i, don't forget node i itself.
 * To save space, Node i is not in the list of firstele[i].
 *
 * by: Wang Yi: ywang@cs.hku.hk
 * @9/12/2011
 */
#ifndef __USETLINK_H_

#define __USETLINK_H_
#include <iostream>
#include "USet.h"
struct UsetNode
{
	int value;
	UsetNode* next;
};

class USetLink
{
public:
	explicit USetLink(USet* set);
	virtual ~USetLink()
	{
		for(int i=0;i<Size;++i)
		{
			while(firstele[i]!=NULL)
			{
				UsetNode* tmp = firstele[i];
				firstele[i] = firstele[i]->next;
				delete tmp;
			}
		}
		delete[]firstele;
		delete[]lastele;
		std::cerr << "USetLink deleted." << std::endl;
	}
	int length()const{return Size;}
	int size()const{return Size;}

	UsetNode** firstele;
	UsetNode** lastele;
private:
	int Size;
	USetLink(){}
	USetLink(const USetLink &uset){}
	const USetLink &operator=(const USetLink &uset){return *this;}
};
#endif
