//#include "Ordering.h"

#ifndef _ORDERING_H__
#define _ORDERING_H__

#include "IntArrayList.h"
#include "BigMatrix.h"

DSS_NAMESPASE_BEGIN

struct IConectMatrix;

/// <summary>
/// Summary description for Ordering.
/// </summary>
class Ordering
{
public:
	enum Type
	{
		None = 0,
		ReverseCuthillMcKee = 1,
		CuthillMcKee = 2,
		MinimumDegree = 3,
		ApproxMinimumDegree = 4,
		ApproxMinimumDegreeIncomplete = 5,
		NestedGraphBisection = 6,
		MetisND = 7,
		ColAMD = 8,
		ApproxMinimumDegreeAA = 9
	};


public:
	IntArrayList* perm;
	IntArrayList* order;
	Type type;
	IConectMatrix* cm;

	Ordering(IntArrayList* perm,IntArrayList* order)
	{
		cm = NULL;
		type = None;
		this->perm  = perm;
		this->order = order;
	}

	Ordering(IntArrayList* order)
	{
		cm = NULL;
		this->order = order;
		perm = new IntArrayList(order->Count);
		perm->Alloc();
		for (long i=0; i<perm->Count; i++)
			perm->Items[order->Items[i]] = i;
	}

	virtual ~Ordering()
	{
		if (perm) {delete perm;perm = NULL;}
		if (order) {delete order;order = NULL;}
		if (cm) {delete cm;cm = NULL;}
	}
};

DSS_NAMESPASE_END

#endif //_ORDERING_H__
