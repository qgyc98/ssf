#include "surffind.h"
#include "oblasti.h"

using std::ostream;

SurfFind::SurfFind(const long nn, const long ne, 
				   const long *nadjelel, const long **adjelel, 
				   const long *nnod, const long **nodes, 
				   const double **xyz)

{
	SpolecnaData sd( nn, ne, nadjelel, adjelel,	nnod, nodes, xyz );
	setgsurface(sd);
}


SurfFind::~SurfFind()
{
	long i,j;
	FOR(i,surfaces)
	{
		FOR(j,nae[i]/2)delete [] nv[i][j];
		delete [] nv[i];
		delete [] nlists1[i];
		delete [] nlists2[i];
		delete [] adjels[i];
	}
	delete [] nn;
	delete [] nae;
	delete [] nv;
	delete [] nlists1;
	delete [] nlists2;
	delete [] adjels;
}


void SurfFind::setgsurface(SpolecnaData & sd)
{
	Oblasti o(sd);
	o.vyhledejOblasti();
	o.hledejStycnePlochy();
	size_t i;
	long j, pocetVektoru;
	surfaces = (long)o.nlists.size();
	nn = new long[surfaces];
	nae = new long[surfaces];
	nlists1 = new long*[surfaces];
	nlists2 = new long*[surfaces];
	adjels = new long*[surfaces];
	nv = new double**[surfaces];
	FOR(i,o.nlists.size()){
		nn[i] = (long)  o.nlists[i].size();
		nae[i] = (long) o.adjels[i].size();		
		nlists1[i] = new long [nn[i]];
		nlists2[i] = new long [nn[i]];
		pocetVektoru = nae[i] / 2;
		nv[i] = new double *[pocetVektoru];
		j = 0;
		for(set<MyPair, Less>::iterator it = o.nlists[i].begin(); it != o.nlists[i].end(); it++)
		{
			nlists1[i][j] = (*it).node1;
			nlists2[i][j] = (*it).node2;	
			j++;
		}
		adjels[i] = new long[nae[i]];
		FOR(j,nae[i])adjels[i][j] = o.adjels[i][j];
		FOR(j,pocetVektoru)nv[i][j] = o.nvectors[i][j];
	}	// FOR i

}


ostream& operator<<(ostream &os, SurfFind &sf)
{
	for (long i = 0; i < sf.surfaces; i++ )
	{
		os << "Stycna plocha " << i << "\n-------------\n\n";
		os << "Pocet prvku [nae] : " << sf.nae[i] << "\n\n";
		os << "Prvky [adjel] : ";
		for(long j = 0; j < sf.nae[i]; j++)
			os << sf.adjels[i][j]+1 << " ";
		os << "\n\n";
		os << "Pocet bodu [nn] : " << sf.nn[i] << "\n\n";
		os << "Dvojice stycnych bodu [nlist1 nlist2 |] : ";
		for(long j = 0; j < sf.nn[i]; j++)
			os << sf.nlists1[i][j]+1 << " " << sf.nlists2[i][j]+1 << " | ";
		os << "\n\n";
		os << "Vektory (x,y,z)| : ";
		for(long j = 0; j < sf.nae[i] / 2; j++)
		{	
			os << "("; 
			for (short k = 0; k < 3; k++) os << sf.nv[i][j][k] << ",";
			os << ")| ";
			
		}
		os << "\n\n\n\n";
	}
	return os;
}

