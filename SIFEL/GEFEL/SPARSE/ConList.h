#ifndef CONLIST_H__
#define CONLIST_H__

#include "Array.h"

DSS_NAMESPASE_BEGIN

class IntLinkArray
{
public:
	long* last;
	long* next;

	long first;
	long N;

	IntLinkArray()
	{
		last = next = NULL;
		first = 0;
		N = 0;
	}
				 	
	IntLinkArray(long* l,long* n,long N)
	{
		last = l;
		next = n;
		first = 0;
		this->N = N;
	}

	void Init(long* l,long* n,long N)
	{
		last = l;
		next = n;
		first = 0;
		this->N = N;
	}

	void Remove(long i)
	{
		//if(!(last[i]!=-2))	{i =i;}

		if (last[i]>=0)	next[last[i]] = next[i];
		if (next[i]>=0)	last[next[i]] = last[i];

		if (i == first)
			first = next[i];
		N--;

		//if(!((last[i]=-2)==-2))	{i =i;}
	}

};

DSS_NAMESPASE_END

#endif //CONLIST_H__
