#ifndef __SURFFIND_H__
#define __SURFFIND_H__

#include <iostream>
#include "spolecnaData.h"



class SurfFind
{
	void setgsurface(SpolecnaData &);

public:

	SurfFind(const long nn, const long ne, const long * nadjelel, const long **adjelel, 
		const long * nnod,const long ** nodes, const double **xyz);

	~SurfFind();

	///  vypis 
	friend std::ostream& operator<<(std::ostream &os, SurfFind &sf);	
	
	/// count of surfaces
	long surfaces;
   
	///  numbers of nodes on surfaces
	long * nn;

	///  numbers of adjacent elements on surfaces
	long * nae;

 	///  numbers of adjacent elements on surfaces
	long ** adjels;

	///  list of nodes on surfaces
	long ** nlists1;
	long ** nlists2;
	
	///  normal vector
	double *** nv;

	///  threshold
	//double threshold;

};


#endif


