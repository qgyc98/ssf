// SkyLineMtxLDL.h

#ifndef _SKYLINEMTXLDL_H__
#define _SKYLINEMTXLDL_H__

#include "SkyLineMtx.h"

DSS_NAMESPASE_BEGIN

/// <summary>
/// Summary description for SkyLineMtxLDL.
/// </summary>
class SkyLineMtxLDL : 
	public SkyLineMtx
{
public:
	SkyLineMtxLDL(SparseMatrixF& sm,Ordering* order,MathTracer* eMT);
	virtual ~SkyLineMtxLDL();
	
	void LoadMatrixData(SparseMatrixF& sm);	


public:
	virtual void Solve(double* b, double* x);

//ILargeMatrix 
	virtual double& ElementAt(int i, int j);
	virtual void LoadZeros();
	virtual void LoadMatrixNumbers(SparseMatrixF& sm);
	virtual void SolveLV(const LargeVector& b, LargeVector& x); 
	virtual void Factorize();
	virtual void MultiplyByVector(const LargeVectorAttach& x, LargeVectorAttach& y);

public:
	virtual void SchurComplementFactorization(int fixed_blocks);
	virtual void SolveA11(double* x,long fixed_blocks);
	virtual void Sub_A21_A11inv(double* x,long fixed_blocks);
	virtual void Sub_A11inv_A12(double* x,long fixed_blocks);
	virtual void WriteCondensedMatrixA22(double* a,Ordering* mcn,IntArrayList* lncn);

}; //class SkyLineMtxLDL 

DSS_NAMESPASE_END

#endif// _SKYLINEMTXLDL_H__
