// SparseGridMtxPD.h

#ifndef _SPARSEGRIDMTXPD_H__
#define _SPARSEGRIDMTXPD_H__

#include "SparseGridMtxLDL.h"

DSS_NAMESPASE_BEGIN

/// <summary>
/// Summary description for SparseGridMtx.
/// </summary>
class SparseGridMtxPD 
{

public:
	// Allocates new space according to bskl and reads old matrix with respect
	// to permutation blockP
	SparseGridMtxPD(SparseMatrixF& sm,BYTE block_size,Ordering* block_order,long fixed_blocks,MathTracer* MT);

	// Allocates new space according to bskl and reads old matrix with respect
	// to permutation blockP
	SparseGridMtxPD(SparseMatrixF& sm,BYTE block_size,Ordering* block_order,Ordering* node_order,long fixed_blocks,MathTracer* MT);

	SparseGridMtxPD(SparseGridMtx* pMatrix,long fixed_blocks);

	virtual ~SparseGridMtxPD();

private:
	SparseGridMtx* m_pMatrix;
	long m_lFixed_blocks;

public:
	// This functions computes LDL' decomposition of first (nblocks-fixed_bn) columns
	// The resulting submatrix will contain the schur complement A22 - A21 * A11^(-1) * A12
	void SchurComplementFactorization();
	void SolveA11(double* x);
	void Sub_A21_A11inv(double* x);
	void Sub_A11inv_A12(double* x);

	void WriteCondensedMatrixA22(double* a,Ordering* mcn,IntArrayList* lncn);

	SparseGridMtx* Matrix();

}; //class SparseGridMtx 

DSS_NAMESPASE_END

#endif// _SPARSEGRIDMTX_H__
