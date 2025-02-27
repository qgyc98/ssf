// DenseMatrix.h

#ifndef _DENSEMATRIX_H__
#define _DENSEMATRIX_H__

#include "DenseMatrixArithmeticsNN.h"

DSS_NAMESPASE_BEGIN

	/// <summary>
	/// Dense matrix stored by columns
	/// </summary>
class DenseMatrix
{
public:
	long n;
	double* data;		// Matrix stored by columns
	  
private:
	DenseMatrixArithmetics* dma;
	  
public:
	DenseMatrixArithmetics& DMA()
	{
		if (dma == NULL) 
			dma = DenseMatrixArithmetics::NewArithmetics(n);
		return *dma;			
	}
	  
	/*public double this[long i, long j]
	{
	get{return data[i+n*j];}
	set{data[i+n*j] = value;}
	}*/
	  
	void Add(long i,long j, double val)
	{
		data[i+n*j] += val;
	}

	DenseMatrix(long n)
	{
		dma = NULL;
		this->n = n;
		data = new double[n*n];
		memset(data,0,n*n*sizeof(double));
	}
	  
	  
	  DenseMatrix(long n,double d)
	{
	  dma = NULL;
	  this->n = n;
	  data = new double[n*n];
	  for (long i=0; i<n; i++) data[i+i*n] = d;
	}
	  
	  DenseMatrix(double* dataFrom,long start_idx,long n)
	{
	  dma = NULL;
	  this->n = n;
	  data = new double[n*n];
	  Array::Copy(dataFrom,start_idx,data,0,n*n);
	}
	  
	  ~DenseMatrix()
	{
	  if (dma)
		delete dma;
	  if (data)
		delete [] data;
	  }



	void CopyTo(DenseMatrix& blockA,long bn)
	{
		Array::Copy(this->data,blockA.data,bn*bn);
	}

	void CopyTo(double* dataTo,long start_idx)
	{
		Array::Copy(data,0,dataTo,start_idx,n*n);
	}

	void Clear()
	{
		Array::Clear(data,0,n*n);
	}

	static double InnerProduct(double *dp1,double *dp2,long len)
	{
		long i,len4;
		double	sum0, sum1, sum2, sum3;

		sum0 = sum1 = sum2 = sum3 = 0.0;

		len4 = len / 4;
		len  = len % 4;

		for ( i = 0; i < len4; i++ )
		{
			sum0 += dp1[4*i]*dp2[4*i];
			sum1 += dp1[4*i+1]*dp2[4*i+1];
			sum2 += dp1[4*i+2]*dp2[4*i+2];
			sum3 += dp1[4*i+3]*dp2[4*i+3];
		}
		sum0 += sum1 + sum2 + sum3;
		dp1 += 4*len4;	dp2 += 4*len4;

		for ( i = 0; i < len; i++ )
			sum0 += (*dp1++)*(*dp2++);

		return sum0;
	}

};

DSS_NAMESPASE_END

#endif //_DENSEMATRIX_H__
