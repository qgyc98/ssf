// BiSection.h 

#ifndef _BISECTION_H__
#define _BISECTION_H__

#include "SparseConectivityMtx.h"

DSS_NAMESPASE_BEGIN

class CMcKee
{
private:
	long n;

	SparseConectivityMtxII* mtx;
	long* p_node_level;
	long* p_order; //out
	// 0 - unsorted available
	//


public:
	long* nodes;  //in
	long size;
	long domA,domB;

	CMcKee();
	~CMcKee();

	void Init(SparseConectivityMtxII* mtx);

	BOOL IsAvailable(int v);

	void PrepareValid();

	long FindFirstNode();
	void ComputeLevels();
	void DivideByMidLevel();
};


class CBiSection
{
private:
	SparseConectivityMtxII* mtx;
	CMcKee mck;

public:
	CBiSection(SparseConectivityMtxII* mtx);
	void RecurBiSectOrder(IntArrayList* order);

private:
	void RecurBiSect(long* nodes,long size);
	void BiSect(long* nodes,long size, long& domA, long& domB);
};

DSS_NAMESPASE_END

#endif //_BISECTION_H__
