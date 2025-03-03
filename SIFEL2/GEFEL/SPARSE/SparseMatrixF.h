// SparseMatrixF.h

#ifndef _SPARSEMATRIXF_H__
#define _SPARSEMATRIXF_H__

#include "DSSAfx.h"

DSS_NAMESPASE_BEGIN

struct IConectMatrix;

struct SparseMatrixF
{
public:
	ULONG neq;
	double* a;

enum eOrientation
{
	eCompressedRows,      //ci - stores column indices [j]
	eCompressedColumns,   //ci - stores row indices [i]
	eSymmetric            // column indices <-> row indices, but only half is stored
};

private:
	ULONG* ci;
	ULONG* adr;
	ULONG Aoffset;
	ULONG Coffset;
	bool bJardaConvention;
	bool bLocalCopy;
	bool m_bIsSymmertric; // TRUE - only diagonal + one half is stored, FALSE - all entries are stored
	eOrientation m_eSparseOrientation;

public:
	SparseMatrixF();
	SparseMatrixF(unsigned long neq,double* a,unsigned long* ci,unsigned long *adr,ULONG Aofs = 0,ULONG Cofs = 0,bool JardaConvention = true,bool bIsSymetric=true,eOrientation sparseOri=eCompressedColumns);
	SparseMatrixF(IConectMatrix* pConMtx);
	~SparseMatrixF();

	void Delete();
	void Detach();
	void CreateLocalCopy();
	void AddNumbers (double alfa,double* a);
	void ScaleMatrix(double alfa);

	long Nonzeros();
	BOOL IsSymmetric() { return m_bIsSymmertric; }
	eOrientation GetOrientation() { return m_eSparseOrientation; }
	BOOL IsSamePattern(SparseMatrixF* sm);
	
	void GetA12block(double* pA12,long c);
	void LoadMatrix(FILE* stream);		
	void SaveMatrix(FILE *stream);

	inline ULONG Adr(int i)	const
	{
		if (bJardaConvention)
		return this->adr[i]-Aoffset;
		else
		{
			if (i==0) 
				return 0;
			else 
				return this->adr[i-1];
		}
	}

	inline ULONG Ci(int i) const
	{
		return (this->ci[i])-Coffset;
	}

	// c = A.b
	void MulMatrixByVector(double *b,double *c);

	// c = A.b
	void MulNonsymMatrixByVector(double *b,double *c);

	// c = A.b
	void MulSymMatrixByVector(double *b,double *c);
	

	/**
	   function multiplies %matrix by %vector

	   @param b - array containing %vector b
	   @param c - array containing resulting %vector c = A.b

	   JK
	*/
	void mxv_scr (double *b,double *c);

	void ReadDiagonal(double* dv);

};

DSS_NAMESPASE_END

#endif//_SPARSEMATRIXF_H__
