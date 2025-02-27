// DSSolver.h

#ifndef _DSSOLVER_H__
#define _DSSOLVER_H__

#include "DSSAfx.h"
#include "SparseConectivityMtx.h"
#include "IntArrayList.h"
//#include "SparseGridMtx.h"
#include "SparseGridMtxPD.h"

DSS_NAMESPASE_BEGIN

enum eDSMatrixType
{
	eDSSparseMatrix = 0,
	eDSSkylineMatrix = 1
};

enum eDSSolverType
{
	eDSSFactorizationLDLT = 0,
	eDSSFactorizationLLT = 1,
	eDSSFactorizationLU = 2,
	eDSSFactorizationLDLTIncomplete = 3,
	eDSSFactorizationLLTIncomplete = 4,
	eDSSDiagonalScaling = 5,
	eDSSFastCG = 6
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///
//		ISolver is the user interface for the DirectSparseSolver
//

struct ISolver 
{
	// Makes necessary intialization
	// 0 - produces sparse symmetrical LDL^T factorization (default)
	// 1 - produces sparse symmetrical LLT Cholesky factorization
	// 2 - produces sparse LU factorization on symmetrical pattern -> undirected graph
	// 3 - produces sparse incomplete LDL^T factorization 
	virtual long Initialize (unsigned char run_code ,eDSSolverType solverType = eDSSFactorizationLDLT, eDSMatrixType matrixType = eDSSparseMatrix) = 0;

	// Sets the type of the reordering of the sparse matrix to something else than default
	// returns true if possible 
	virtual BOOL SetOrderingType(Ordering::Type otype) = 0;

	// Registers the source to the DSSolver
	virtual BOOL LoadMatrix (ULONG neq,unsigned char block_size,double * a,ULONG * ci,ULONG * adr ) = 0;
	virtual BOOL LoadMatrix (SparseMatrixF* sm,unsigned char block_size) = 0;
	virtual BOOL SetMatrixPattern(SparseMatrixF* smt,unsigned char block_size) = 0;

	// Loads the "MatrixCodeNumbers" vector[n_blocks*block_size]
	// each entry in the vector[i] can be :
	// vector[i] >=  0  normal unknown			(means the corresponding row in SparseMatrixF [zero-based])
	// vector[i] == -1  this entry is not used	(no corresponding row in SparseMatrixF)
	// vector[i] <= -2  unknown to be left uncondensed	( -vector[i]-2 is the zero-based row index in SparseMatrixF)

	//  Ex:
		// 0,1,2,....... code numbers
		// -1   ........ eliminated row
		// -2, -3, -4 .. rows (0,1,2) ment to be left uncondensed 
	//  n_blocks		- number of blocks
	
	virtual BOOL LoadMCN (ULONG n_blocks,unsigned char block_size,long * mcn, BOOL bIsSchur ) = 0;


	// Creates the connectivity matrix 
	// Computes the MinimumDegree ordering
	// Allocates the memory for the SparseGrid matrix
	virtual BOOL PreFactorize ( ) = 0;

	// Loads zeros to the already allocated SparseGrid matrix 
	virtual void LoadZeros() = 0;

	// Loads or reloads the numeric values of the sparse matrix 'sm' to the already allocated SparseGrid matrix 
	virtual BOOL LoadNumbers (SparseMatrixF* sm) = 0;

	// Adds alfa multiple of matrix of the same size to already allocated matrix (only array a)
	virtual BOOL AddNumbers (double alfa,SparseMatrixF* smtx) = 0;

	// Muliplies the matrix by specified value
	virtual BOOL ScaleMatrix(double alfa) = 0;

	virtual SparseMatrixF* GetSparseMatrix() = 0;

	// Runs the LDL factorization on the already allocated and loaded matrix
	virtual BOOL ReFactorize ( ) = 0;

	// This function unites the three previous function into one step if necessary
	// PreFactorize ( ) + LoadNubers( ) + ReFactorize ( )
	// This is usefull for single-pass solutions
	virtual BOOL Factorize ( ) = 0;

	// When the matrix was factorized we can solve several A*r=f equations
	// if the pointers are equal (r==f) the result will overwrite the RHS vector (In that case the soultion is faster)
	// It is recomended that both vectors have size in a whole multiple of the 'block_size' or more
	virtual BOOL Solve (double * r,double * f ) = 0;

	// When the matrix was factorized we solve r from A*r=f equation
	// where f is vector with units. And we return (f-A*r)|(f-A*r)/N value.
	virtual double GetFactorizationError() = 0;

	// Disposes the solver internal data from the memory for us to be able to load another matrix.
	virtual void Dispose() = 0;

	// Contains calls Dispose function.
	virtual long Close ( ) = 0;

	virtual void SetMT(MathTracer* MT) = 0;
	virtual void SetMT() = 0;
	virtual void SetMTN() = 0;
    virtual void SetMTMtx() = 0;

	virtual BOOL IsFactorized() = 0;
	virtual BOOL IsAllocated() = 0;
	virtual BOOL IsInitialized() = 0;

	virtual ~ISolver(){};

	// |A11 A12| |r1| = |f1|
	// |A21 A22| |r2|   |f2|

	//  a - condensed matrix 
	//  lhs - left hand side 
	//  rhs - right hand side
	//  tc - type of computation
	//  tc=1 - return condensed matrix a=(A22-A21*A11inv*A12) and modify rhs part f2 = (f2-A21*A11inv*f1)
	//  tc=2 - return part of left hand side
	//  tc=3 - return (Kii lhs = rhs) solution
	virtual void condense(double *a,double *lhs,double *rhs,long tc) = 0;

	// A few words to different ordering systems used during sparse PD solution
	//
	// A-order [neq]  
	//		This is the order of data stored in user specified SparseMatrixF* sm.
	//		Usually, we get RHS and return LHS in this order system.
	//
	// B-order [neq] -> [n_blocks*block_size]
	//		You obtain B-order after applying MCN vector on A-order. 
	//		MCN re-introduces removed lines and DOFs from certain nodes.
	//
	//	ex:		for (i=0;i<neq;i++) lhsA[i] = tmpB[mcn->order[i]];
	//			for (i=0;i<neq;i++) tmpB[mcn->order[i]] = lhsA[i];
	//
	// C-order [neq]
	//		This is such a permutaion of A order where all noncondensed nodes are on the end of 
	//		the order.
	//					dom_order

	virtual void GetA12block(double *pA12) = 0;

	// c += A * b
	virtual void MulMatrixByVector(double *b, double *c) = 0;

	// Copute x by Preconditioned Conjugate gradient method
	virtual int PreCG(double* b, double* x,double epsilon,int max_iter) = 0;

	// Copute x by Conjugate gradient method
	virtual int CG(double* b, double* x,double epsilon,int max_iter) = 0;

	enum eState
	{
		None = 0,
		Initialized = 1, // basic construction
		Allocated = 2,   // PreFactorized
		Factorized = 3,  // Factorized
		ErrorInMCN = -1 
	};
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////

class DSSolver : public ISolver
{
	MathTracer MT;
    MathTracerNull MTN;
    MathTracerFile MTF;
	MathTracer* eMT;

	SparseMatrixF m_sm;
	SparseGridMtx* act_matrix; // points to matrix or matrixPD
	SparseGridMtx* matrix;
	SparseGridMtxPD* matrixPD;
	char str[512];
	ISolver::eState m_eState;

	eDSSolverType SolverType;
	eDSMatrixType MatrixType;
	unsigned char run_code;

	// PDD data
	double* tmpR;
	IntArrayList* dom_order;	

	// blocked matrix
	unsigned char blockSize;
	long n_blocks;
	long neq;
	Ordering* mcn;

	// Diagonal precoditioner
	LargeVector matrixD;
	SparseGridMtx* orig_matrix; // blockmatrix for CG multiplication

	// Dirichlets conditions - noncondensed DOFs
	IntArrayList* fixed;		// noncondensed blocks
	IntArrayList* lncn;			// noncondensed DOFs

	Ordering::Type	OrderingType;	//MinimumDegree,ApproxMinimumDegree,ApproxMinimumDegreeIncomplete,...
	BOOL m_bIsSchur; 

public:
	DSSolver(MathTracer* pMT = NULL);

	virtual ~DSSolver();
	virtual long Initialize (unsigned char run_code ,eDSSolverType solverType = eDSSFactorizationLDLT, eDSMatrixType matrixType  = eDSSparseMatrix);
	virtual void SetMT(MathTracer* pMT);
    virtual void SetMTN() {eMT = &MTN;};
    virtual void SetMT()  {eMT = &MT;};
    virtual void SetMTMtx() {if (matrix) matrix->SetMT(eMT);};
    virtual MathTracer* GivePMT(){return eMT;};
	virtual void Dispose();

	virtual BOOL SetOrderingType(Ordering::Type otype);
	virtual BOOL LoadMatrix(unsigned long neq,unsigned char block_size,double * a,unsigned long * ci,unsigned long * adr );
	virtual BOOL LoadMatrix(SparseMatrixF* smt,unsigned char block_size);
	virtual BOOL SetMatrixPattern(SparseMatrixF* smt,unsigned char block_size);

	virtual BOOL IsFactorized();
	virtual BOOL IsAllocated();
	virtual BOOL IsInitialized();
	virtual BOOL IsSchur(); // are some DOFs fixed for condensation?

	virtual BOOL Factorize ( );
	virtual BOOL PreFactorize();

	virtual void LoadZeros();
	virtual BOOL LoadNumbers (SparseMatrixF* sm);
	virtual BOOL AddNumbers (double alfa,SparseMatrixF* smtx);
	virtual BOOL ScaleMatrix(double alfa);
	virtual SparseMatrixF* GetSparseMatrix();

	virtual double& ElementAt(int i, int j);

	virtual BOOL ReFactorize( );

	virtual BOOL Solve (double * r,	double * f );
	virtual long Close ( );

	void StartSolverWriteInfo();
	void EndSolverWriteInfo();

	void SetSM(SparseMatrixF* sm);
	virtual double GetFactorizationError();

	// Loads the "MatrixCodeNumbers" vector[n_blocks*block_size]
	// each entry in the vector[i] can be :
	// vector[i] >=  0  normal unknown			(means the corresponding row in SparseMatrixF)
	// vector[i] == -1  this entry is not used	(no corresponding row in SparseMatrixF)
	virtual BOOL LoadMCN (ULONG n_blocks,unsigned char block_size,long * mcn, BOOL bIsSchur );

	virtual BOOL LoadMCN (IntArrayList& mcn);

	IntArrayList* GetFixedBlocks() { return this->fixed; }
	IntArrayList* GetFixedDOFs()   { return this->lncn; }

	virtual void condense(double *a,double *lhs,double *rhs,long tc);
	
	void GetA12block(double *pA12);

	// c = A * b
	void MulMatrixByVector(double *b, double *c);

	int PreCG(double* b, double* x,double epsilon,int max_iter);

	// Copute x by Conjugate gradient method
	int CG(double* b, double* x,double epsilon,int max_iter);

private:
	BOOL LoadMCN_int(IntArrayList* mcn_order);
	void ExpandMCN(IntArrayList& mcn);
	SparseGridMtx* CreateNewSparseGridMtx(IntArrayList* fixed = NULL);
	void WriteFactorizationInfo();
	BOOL CreateFixedArray(long no_noncondensed_DOFs);
	BOOL PreFactorizeSchur();
	void StoreFixedLastPermutation_dom_order();
};

DSS_NAMESPASE_END

#endif //_DSSOLVER_H__
