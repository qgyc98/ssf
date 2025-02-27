// SkyLineMtxLDL.cpp

#include "SkyLineMtxLDL.h"

DSS_NAMESPASE_BEGIN

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SkyLineMtxLDL::SkyLineMtxLDL(SparseMatrixF& sm,Ordering* order,MathTracer* eMT):SkyLineMtx(sm,order,eMT)
{	
	nonzeros = 0;
	this->eMT = eMT;
	this->order = order;
}

void SkyLineMtxLDL::LoadZeros()
{
}

void SkyLineMtxLDL::Solve(double* b, double* x)
{
	UNREFERENCED_PARAMETER(b);
	UNREFERENCED_PARAMETER(x);
}

double dummy = 0;
double& SkyLineMtxLDL::ElementAt(int i, int j)
{
	return dummy;
}

void SkyLineMtxLDL::LoadMatrixData(SparseMatrixF& sm)	
{
	long* nodeP = this->order?order->perm->Items:NULL;

	// diagonal
	for (ULONG j=0; j<sm.neq; j++)
	{
		long nj = (nodeP==NULL)?j:nodeP[j];
		int col_start = column_starts[nj]+nj-1;

		for (ULONG ad = sm.Adr(j); ad<sm.Adr(j+1); ad++)
		{
			long i = (long)sm.Ci(ad);
			long ni = (nodeP==NULL)?i:nodeP[i];
			double val = sm.a[ad];

			if (ni<nj)
				columndata[col_start-ni] = val;
			else if (ni==nj)
				D[ni] = val;
		}
	}
}


SkyLineMtxLDL::~SkyLineMtxLDL()
{
}

void SkyLineMtxLDL::LoadMatrixNumbers(SparseMatrixF& sm)
{
}
void SkyLineMtxLDL::SolveLV(const LargeVector& b, LargeVector& x)
{
} 
void SkyLineMtxLDL::Factorize()
{
}
void SkyLineMtxLDL::MultiplyByVector(const LargeVectorAttach& x, LargeVectorAttach& y) 
{
}

void SkyLineMtxLDL::SchurComplementFactorization(int fixed_blocks)
{
}

void SkyLineMtxLDL::SolveA11(double* x,long fixed_blocks) 
{
}

void SkyLineMtxLDL::Sub_A21_A11inv(double* x,long fixed_blocks) 
{
}

void SkyLineMtxLDL::Sub_A11inv_A12(double* x,long fixed_blocks) 
{
}
void SkyLineMtxLDL::WriteCondensedMatrixA22(double* a,Ordering* mcn,IntArrayList* lncn) 
{
}

DSS_NAMESPASE_END
