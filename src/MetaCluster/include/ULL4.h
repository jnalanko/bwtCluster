/*
 * last fixed: 2012.01.05.
 * by Wang Yi.
 * */
#ifndef __ULL4_H_

#define __ULL4_H_
#include <iostream>
using namespace std;
class ULL4
{
public:
    unsigned long long high,mid_H,mid_L,low;

    ULL4();
    ULL4(int i);
    ULL4(const ULL4& obj);
    ULL4(unsigned long long high_,unsigned long long mid_H_,unsigned long long mid_L_, unsigned long long low_);
    void setzero();
    int non0base();//not non 0 bits

    bool operator < (const ULL4& obj) const;
    bool operator == (const ULL4& obj) const;
    bool operator != (const ULL4& obj) const;
    ULL4& operator <<= (const int obj);
    ULL4  operator << (const int obj) const;
    ULL4& operator >>= (const int obj);
    ULL4  operator >> (const int obj) const;
    ULL4& operator |= (const int obj);
    ULL4& operator |= (const ULL4& obj);
    ULL4 operator | (const ULL4& obj)const;
    ULL4& operator &= (const ULL4& obj);
    ULL4 operator & (const ULL4& obj)const;
    ULL4& operator ^= (const ULL4& obj);
    ULL4  operator ^ (const ULL4& obj) const;
    ULL4 operator~() const;
 //   void setbase(int position,unsigned long long base);
    friend ostream& operator<<(ostream& os,const ULL4 & obj);
};
#endif
