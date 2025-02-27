// ColHash.h

#ifndef _COLHASH_H__
#define _COLHASH_H__

#include "Array.h"
#include "IntArrayList.h"

DSS_NAMESPASE_BEGIN

class ColHash
{
public:
	IntArrayList** buckets;
	IntArrayList occupied_buckets;
	int n;

	ColHash()
	{
		buckets = NULL;
	}

	ColHash(long n)
	{
		buckets = new IntArrayList*[n];
		memset(buckets,0,n*sizeof(IntArrayList*));
		this->n = n;
		//for (long i=0; i<n; i++)
		//	buckets[i] = new IntArrayList();
	}

	void Init(long n)
	{
		if (buckets) delete [] buckets;
		buckets = new IntArrayList*[n];
		memset(buckets,0,n*sizeof(IntArrayList*));
		this->n = n;
	}

	~ColHash()
	{
		for (int i=0; i<n; i++)
			if (buckets[i]) {delete buckets[i];buckets[i] = NULL;}
		if (buckets) delete [] buckets;
		buckets = NULL;
	}

	void AddValue(long bucket,long obj)
	{
		if (buckets[bucket]==NULL)
			buckets[bucket] = new IntArrayList();
		//buckets[bucket].Add(obj);

		long i = buckets[bucket]->Add(obj);
		if (i==0)
			occupied_buckets.Add(bucket);
	}

	void Clear()
	{
		long*pI = occupied_buckets.Items;
		for(long* p = pI+occupied_buckets.Count-1; p>=pI; p--)
			buckets[*p]->Clear();
		occupied_buckets.Clear();
	}
};

DSS_NAMESPASE_END

#endif //_COLHASH_H__
