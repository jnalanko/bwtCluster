/*
 * last fixed: 2012.01.05.
 * union set with relative position recorded.
 * by Wang Yi. ywang@cs.hku.hk
 * @12/09/2011
 */
#ifndef __RELAPOSI_H_

#define __RELAPOSI_H_
#include "USet.h"
class RelaPosi
{
public:
	explicit RelaPosi(int size_)
	{
		Size = size_;
		anchoridL = new int[size_];
		anchoridR = new int[size_];

		positionL = new int[size_];
		positionR = new int[size_];

		scoreL = new int[size_];
		scoreR = new int[size_];

		isRevL = new bool[size_];
		isRevR = new bool[size_];

		FreqL = new int[size_];
		FreqR = new int[size_];

		for(int i=0;i<size_;++i)
			positionL[i]=128;
		for(int i=0;i<size_;++i)
			positionR[i]=128;
		for(int i=0;i<size_;++i)
			anchoridL[i]=-1;
		for(int i=0;i<size_;++i)
			anchoridR[i]=-1;
		for(int i=0;i<size_;++i)
			isRevL[i]=false;
		for(int i=0;i<size_;++i)
			isRevR[i]=false;
		for(int i=0;i<size_;++i)
			scoreL[i]=0;
		for(int i=0;i<size_;++i)
			scoreR[i]=0;
		for(int i=0;i<size_;++i)
			FreqL[i]=0;
		for(int i=0;i<size_;++i)
			FreqR[i]=0;
	}
	virtual ~RelaPosi()
	{
		delete []positionL;
		delete []anchoridL;
		delete []isRevL;
		delete []scoreL;
		delete []positionR;
		delete []anchoridR;
		delete []isRevR;
		delete []scoreR;
		delete []FreqL;
		delete []FreqR;
	}
	int length()const{return Size;}
	int size()const{return Size;}

	//position(y) equals position(x)>>posi
	void setPosi(int x,int y,int deltaposi,int isRevx,int isRevy,int score);
	void getFreq();
	void printInfo(int id);

	int *anchoridL, *anchoridR;
	int *scoreL, *scoreR;
	int* positionL, *positionR;
	bool* isRevL, *isRevR;
private:
	int* FreqL, *FreqR;
	int Size;

	RelaPosi(){}
	RelaPosi(const RelaPosi &uset){}
	const RelaPosi &operator=(const RelaPosi &uset){return *this;}
};
#endif
