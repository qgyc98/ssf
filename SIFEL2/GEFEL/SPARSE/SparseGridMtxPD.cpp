// SparseGridMtxPD.cpp

#include "SparseGridMtxPD.h"

DSS_NAMESPASE_BEGIN

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SparseGridMtxPD::SparseGridMtxPD(SparseMatrixF& sm,BYTE block_size,Ordering* block_order,long fixed_blocks,MathTracer* MT)
{	
	this->m_pMatrix = new SparseGridMtxLDL(sm,block_size,block_order,MT);
	this->m_lFixed_blocks = fixed_blocks;
}

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SparseGridMtxPD::SparseGridMtxPD(SparseMatrixF& sm,BYTE block_size,Ordering* block_order,Ordering* node_order,long fixed_blocks,MathTracer* MT)
{	
	this->m_pMatrix = new SparseGridMtxLDL(sm,block_size,block_order,node_order,MT);
	this->m_lFixed_blocks = fixed_blocks;
}

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SparseGridMtxPD::SparseGridMtxPD(SparseGridMtx* pMatrix,long fixed_blocks)
{	
	this->m_pMatrix = pMatrix;
	this->m_lFixed_blocks = fixed_blocks;
}

SparseGridMtxPD::~SparseGridMtxPD()
{
	delete m_pMatrix;
}

SparseGridMtx* SparseGridMtxPD::Matrix()
{
	return m_pMatrix;
}

// This functions computes LDL' decomposition of first (nblocks-fixed_bn) columns
// The resulting submatrix will contain the schur complement A22 - A21 * A11^(-1) * A12
void SparseGridMtxPD::SchurComplementFactorization()
{
	m_pMatrix->SchurComplementFactorization(m_lFixed_blocks);
}

void SparseGridMtxPD::SolveA11(double* x)
{
	m_pMatrix->SolveA11(x,m_lFixed_blocks);
}

void SparseGridMtxPD::Sub_A21_A11inv(double* x)
{
	m_pMatrix->Sub_A21_A11inv(x,m_lFixed_blocks);
}

void SparseGridMtxPD::Sub_A11inv_A12(double* x)
{
	m_pMatrix->Sub_A11inv_A12(x,m_lFixed_blocks);
}

void SparseGridMtxPD::WriteCondensedMatrixA22(double* a,Ordering* mcn,IntArrayList* lncn)
{
	m_pMatrix->WriteCondensedMatrixA22(a,mcn,lncn);
}

DSS_NAMESPASE_END

