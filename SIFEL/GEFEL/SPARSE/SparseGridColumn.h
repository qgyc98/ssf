#ifndef _GRIDCOLUMN_H__
#define _GRIDCOLUMN_H__

#include "IntArrayList.h"
#include "Array.h"

DSS_NAMESPASE_BEGIN

/// <summary>
/// Summary description for SparseGridColumn.
/// </summary>
class SparseGridColumn
{
	// used during sealed state
public:
	long column_start_idx;	// pointer into Columns_data array
	IntArrayList* IndexesUfa;
	long Entries;

	SparseGridColumn(IntArrayList* indexes)
	{
		Entries = 0;
		column_start_idx = -1;
		IndexesUfa = indexes;
		Entries = indexes->Count;
	}

	~SparseGridColumn()
	{
		if (IndexesUfa) {delete IndexesUfa;IndexesUfa=NULL;}
	}

	// Sealed procedures
	//---------------------------

	// Column must not be sealed
	long FindExistingBlockIndex(long bi)	// returns the pointer in the data
	{
		long myIndex=IndexesUfa->BinarySearch(bi);
		//System.Diagnostics.Debug.Assert(myIndex>=0,"Index not found !");
		//if (myIndex<0)
		//	return 0;
		return myIndex;
	}

	void SetValue(long bn,long bi, long si, long sj,double val,double* Columns_data,long &aux_idx)
	{
		if ((aux_idx>=IndexesUfa->Count) || (IndexesUfa->Items[aux_idx]!=bi))
			aux_idx = FindExistingBlockIndex(bi);
		if (aux_idx>=0)
			Columns_data[column_start_idx+bn*(bn*aux_idx+sj) + si] = val;
	}

	double GetValue(long bn,long bi, long si, long sj,double* Columns_data,long &aux_idx)
	{
		if ((aux_idx>=IndexesUfa->Count) || aux_idx<0 || (IndexesUfa->Items[aux_idx]!=bi))
			aux_idx = FindExistingBlockIndex(bi);
		if (aux_idx<0)
			return 0;
		return Columns_data[column_start_idx+bn*(bn*aux_idx+sj) + si];
	}

	void AddValue(long bn,long bi, long si, long sj,double val,double* Columns_data)
	{
		long newBlock = FindExistingBlockIndex(bi);
		Columns_data[column_start_idx+bn*bn*newBlock + si + sj*bn] += val;
	}
};

DSS_NAMESPASE_END

#endif //GRIDCOLUMN_H__
