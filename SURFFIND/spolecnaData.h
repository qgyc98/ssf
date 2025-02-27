#ifndef __SPOLECNADATA_H
#define __SPOLECNADATA_H

// trida SpolecnaData slouzi jako rozhrani obsahujici vsechna dulezita data


class SpolecnaData 
{

public:
	// celkovy pocet bodu
	const long nn; 
	// ne = celkovy pocet prvku
	const long ne;
	// nadjelel = pocet sousednich prvku prvku
	const long *nadjelel;
	// adjelel = seznam sousednich prvku prvku
	const long **adjelel;
	// nnod = pocet bodu prvku
	const long * nnod;
	// nodes = seznam bodu prvku
	const long ** nodes;
	// pole souradnic. Radkovy index je cislo bodu. Sloupcovy soradnice x, y, z
	const double **xyz;
	

SpolecnaData (	const long nn, const long ne, const long *nadjelel, 
			    const long **adjelel, const long * nnod, 
				const long ** nodes, const double **xyz):
nn(nn),
ne(ne),	
nadjelel(nadjelel),
adjelel(adjelel), 
nnod(nnod),
nodes(nodes), 
xyz(xyz) 
{

}

static short * poleShort(const long pocetSloupcu,const short hodnota=0)
{
	short * pole = new short[pocetSloupcu];
	for (long j = 0; j < pocetSloupcu; j++)
	{
		pole[j] = hodnota;
	}
	return pole;
}

~SpolecnaData() {}

};

#endif

