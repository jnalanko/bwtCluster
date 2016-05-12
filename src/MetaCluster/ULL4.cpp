#include <iostream>

#include "ULL4.h"
using namespace std;

ULL4::ULL4()
{
	high = 0ULL;
	mid_H =0ULL;
	mid_L =0ULL;
	low = 0ULL;
}
ULL4::ULL4(int i)
{
	high = 0ULL;
	mid_H = 0ULL;
	mid_L = 0ULL;
	low = i;
}
ULL4::ULL4(const ULL4& obj)
{
	high = obj.high;
	mid_H = obj.mid_H;
	mid_L = obj.mid_L;
	low = obj.low;
}
ULL4::ULL4(unsigned long long high_, unsigned long long mid_H_, unsigned long long mid_L_, unsigned long long low_)
{
	high = high_;
	mid_H = mid_H_;
	mid_L = mid_L_;
	low = low_;
}

void ULL4::setzero()
{
	high = 0ULL;
	mid_H = 0ULL;
	mid_L = 0ULL;
	low  = 0ULL;
}

bool ULL4::operator < (const ULL4& obj) const
{
	if(this->high < obj.high)return true;
	else if(this->high == obj.high && this->mid_H < obj.mid_H)return true;
	else if(this->high == obj.high && this->mid_H == obj.mid_H && this->mid_L < obj.mid_L)return true;
	else if(this->high == obj.high && this->mid_H == obj.mid_H && this->mid_L ==obj.mid_L && this->low < obj.low)return true;
	else return false;
}

inline bool ULL4::operator == (const ULL4& obj) const
{
	return (this->high == obj.high && this->mid_H == obj.mid_H && this->mid_L == obj.mid_L && this->low == obj.low);
}

inline bool ULL4::operator != (const ULL4& obj) const
{
	return !(this->high == obj.high && this->mid_H == obj.mid_H && this->mid_L == obj.mid_L && this->low == obj.low);
}

ULL4 ULL4::operator << (const int obj) const
{
	unsigned long long high_t=0ULL,mid_H_t=0ULL,mid_L_t=0ULL,low_t=0ULL;
	if(obj > 192)
	{
		high_t = low<<(obj-192);
	}
	else if(obj == 192)
	{
		high_t = low;
	}
	else if(obj > 128)
	{
		high_t  = mid_L << (obj-128);
   		high_t |= low >> (192-obj);
		mid_H_t = low << (obj-128);
	}
	else if(obj == 128)
	{
		high_t = mid_L;
		mid_H_t = low;
	}
	else if(obj > 64)
	{
		high_t   = (mid_H << (obj-64));
		high_t  |= (mid_L >> (128-obj));
		mid_H_t  = (mid_L << (obj-64));
		mid_H_t |= (low >> (128-obj));
		mid_L_t  = (low << (obj-64));
	}
	else if(obj == 64)
	{
		high_t   = mid_H;
		mid_H_t  = mid_L;
		mid_L_t  = low;
	}
	else if(obj>0)
	{
		high_t   = high << obj;
		high_t  |= mid_H >> (64-obj);
		mid_H_t  = mid_H << obj;
		mid_H_t |= mid_L >> (64-obj);
		mid_L_t  = mid_L << obj;
		mid_L_t |= low >> (64-obj);
		low_t  = low << obj;
	}
	else
	{
		high_t = high;
		mid_H_t  = mid_H;
		mid_L_t  = mid_L;
		low_t  = low;
	}
	return ULL4(high_t,mid_H_t,mid_L_t,low_t);
}

ULL4& ULL4::operator <<= (const int obj) 
{
	if(obj > 192)
	{
		high = low<<(obj-192);
		mid_H=0ULL;
		mid_L=0ULL;
		low = 0ULL;
	}
	else if(obj == 192)
	{
		high = low;
		mid_H=0ULL;
		mid_L=0ULL;
		low = 0ULL;
	}
	else if(obj > 128)
	{
		high    = mid_L << (obj-128);
   		high   |= low >> (192-obj);
		mid_H   = low << (obj-128);
		mid_L=0ULL;
		low = 0ULL;
	}
	else if(obj == 128)
	{
		high  = mid_L;
		mid_H = low;
		mid_L = 0ULL;
		low   = 0ULL;
	}
	else if(obj > 64)
	{
		high   = (mid_H << (obj-64));
		high  |= (mid_L >> (128-obj));
		mid_H  = (mid_L << (obj-64));
		mid_H |= (low >> (128-obj));
		mid_L  = (low << (obj-64));
		low   = 0ULL;
	}
	else if(obj == 64)
	{
		high  = mid_H;
		mid_H = mid_L;
		mid_L = low;
		low   = 0ULL;
	}
	else if(obj>0)
	{
		high <<= obj;
		high  |= mid_H >> (64-obj);
		mid_H <<= obj;
		mid_H |= mid_L >> (64-obj);
		mid_L <<= obj;
		mid_L |= low >> (64-obj);
		low <<= obj;
	}
	return (*this);
}

int ULL4::non0base()
{
	int result = 0;
	unsigned long long v=((high>>1)|(high))&0x5555555555555555;
	for(;v;++result)
	{
		v &= v-1;
	}
	v=((mid_H >>1)|(mid_H ))&0x5555555555555555;
	for(;v;++result)
	{
		v &= v-1;
	}
	v=((mid_L >>1)|(mid_L ))&0x5555555555555555;
	for(;v;++result)
	{
		v &= v-1;
	}
	v=((low >>1)|(low ))&0x5555555555555555;
	for(;v;++result)
	{
		v &= v-1;
	}
	return result;
}

ULL4 ULL4::operator >> (const int obj) const
{
	unsigned long long high_t=0ULL,mid_H_t=0ULL,mid_L_t=0ULL,low_t=0ULL;
	if(obj > 192)
	{
		low_t = high>>(obj-192);
	}
	else if(obj == 192)
	{
		low_t = high;
	}
	else if(obj >128)
	{
		low_t   = mid_H >> (obj-128);
   		low_t  |= high <<(128-obj);
		mid_L_t = high >> (obj-128);
	}
	else if(obj ==128)
	{
		low_t = mid_H;// >> (obj-64);
		mid_L_t  = high;// >> (obj-64);
	}
	else if(obj > 64)
	{
		low_t   =  mid_L >> obj-64;
		low_t   |= mid_H << (128-obj);
		mid_L_t =  mid_H >> obj-64;
		mid_L_t |= high << (128-obj);
		mid_H_t =  high >> obj-64;
	}
	else if(obj == 64)
	{
		mid_H_t = high;
		mid_L_t  = mid_H ;
		low_t  = mid_L;
	}
	else if(obj > 0)
	{
		low_t  = low >> obj;
		low_t |= mid_L << (64-obj);
		mid_L_t  = mid_L >> obj;
		mid_L_t |= mid_H << (64-obj);
		mid_H_t  = mid_H >> obj;
		mid_H_t |= high << (64-obj);
		high_t = high >> obj;
	}
	else
	{
		high_t = high;
		mid_H_t  = mid_H;
		mid_L_t  = mid_L;
		low_t  = low;
	}
	return ULL4(high_t,mid_H_t,mid_L_t,low_t);
}

ULL4& ULL4::operator >>= (const int obj) 
{
	unsigned long long high_t=0ULL,mid_H_t=0ULL,mid_L_t=0ULL,low_t=0ULL;
	if(obj > 192)
	{
		low = high>>(obj-192);
		mid_L = 0ULL;
		mid_H = 0ULL;
		high = 0ULL;
	}
	else if(obj == 192)
	{
		low_t = high;
		mid_L = 0ULL;
		mid_H = 0ULL;
		high = 0ULL;
	}
	else if(obj >128)
	{
		low   = mid_H >> (obj-128);
   		low  |= high <<(128-obj);
		mid_L = high >> (obj-128);
		mid_H = 0ULL;
		high = 0ULL;
	}
	else if(obj ==128)
	{
		low = mid_H;// >> (obj-64);
		mid_L  = high;// >> (obj-64);
		mid_H = 0ULL;
		high = 0ULL;
	}
	else if(obj > 64)
	{
		low   =  mid_L >> obj-64;
		low   |= mid_H << (128-obj);
		mid_L =  mid_H >> obj-64;
		mid_L |= high << (128-obj);
		mid_H =  high >> obj-64;
		high = 0ULL;
	}
	else if(obj == 64)
	{
		low   = mid_L;
		mid_L = mid_H ;
		mid_H = high;
		high = 0ULL;
	}
	else if(obj > 0)
	{
		low >>= obj;
		low |= mid_L << (64-obj);
		mid_L >>= obj;
		mid_L |= mid_H << (64-obj);
		mid_H >>= obj;
		mid_H |= high << (64-obj);
		high >>= obj;
	}
	return *this;
}

ULL4& ULL4::operator |= (const int obj)
{
	low |= obj;
	return *this;
}

ULL4& ULL4::operator |= (const ULL4& obj)
{
	low  |= obj.low;
	mid_L|= obj.mid_L;
	mid_H|= obj.mid_H;
	high |= obj.high;
	return *this;
}

ULL4 ULL4::operator | (const ULL4& obj)const
{
	ULL4 result;
	result.low = this->low | obj.low;
	result.mid_L = this->mid_L | obj.mid_L;
	result.mid_H = this->mid_H | obj.mid_H;
	result.high= this->high| obj.high;
	return result;
}

ULL4 ULL4::operator ^ (const ULL4& obj)const
{
	ULL4 result;
	result.low = this->low ^ obj.low;
	result.mid_L = this->mid_L ^ obj.mid_L;
	result.mid_H = this->mid_H ^ obj.mid_H;
	result.high= this->high^ obj.high;
	return result;
}

ULL4 ULL4::operator & (const ULL4& obj)const
{
	ULL4 result;
	result.low = this->low & obj.low;
	result.mid_L = this->mid_L & obj.mid_L;
	result.mid_H = this->mid_H & obj.mid_H;
	result.high= this->high& obj.high;
	return result;
}

ULL4& ULL4::operator ^= (const ULL4& obj)
{
	low   ^= obj.low;
	mid_L ^= obj.mid_L;
	mid_H ^= obj.mid_H;
	high  ^= obj.high;
	return *this;
}
ULL4& ULL4::operator &= (const ULL4& obj)
{
	low   &= obj.low;
	mid_L &= obj.mid_L;
	mid_H &= obj.mid_H;
	high  &= obj.high;
	return *this;
}

ULL4 ULL4::operator~() const
{
	return ULL4(~(this->high),~(this->mid_H),~(this->mid_L),~(this->low));
}

//position is [1,75];
/*void ULL4::setbase(int position,unsigned long long base)
{
	if(position > 97)
	{
		position -= 97;
		unsigned long long mask = 3ULL << (position*2);
		high &= (~mask);
		high |= (base << (position*2));
	}
	else if(position == 97)
	{
		high &= (~3ULL);
		high |= base;
	}
	else if(position > 65)
	{
		position -= 65;
		unsigned long long mask = 3ULL << (position*2);
		high &= (~mask);
		high |= (base << (position*2));
	}
	else if(position == 65)
	{
		high &= (~3ULL);
		high |= base;
	}
	else if(position > 33)
	{
		position -= 33;
		unsigned long long mask = 3ULL << (position*2);
		mid &= (~mask);
		mid |= (base << position*2);
	}
	else if(position == 33)
	{
		mid &= (~3ULL);
		mid |= base;
	}
	else if(position > 1)
	{
		--position;
		unsigned long long mask = 3ULL << (position*2);
		low &= (~mask);
		low |= (base << position*2);
	}
	else
	{
		low &= (~3ULL);
		low |= base;
	}
}*/

ostream& operator<<(ostream& os,const ULL4 &obj)
{
	os<<hex<<obj.high<<'-'<<obj.mid_H<<'_'<<obj.mid_L<<'-'<<obj.low;
	return os;
}
