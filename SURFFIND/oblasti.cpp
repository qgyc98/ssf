#include "oblasti.h"	
#include <iostream>
#include <vector>

// plocha jedne krychle nemuze sousedit s libovolnou plochou jine
// zde uvedene dvojice vyjadruji cisla ploch dvou ruznych krychli ktere spolu mohou sousedit
/*
// **** NEFUNGUJE OBECNE  ****
const short sousednostKrychle [KR_POCET_PLOCH][2] =	{
														{0,2},	// 1,3
														{1,3},	// 2,4
														{2,0},	// 3,1
														{3,1},	// 4,2
														{4,5},	// 5,6
														{5,4}	// 6,5
													};

*/
using namespace std;


Oblasti::Oblasti(SpolecnaData & sd):
sd (sd)
{
	short mez;
	poplpr = new short [sd.ne];
	seplpr = new Plocha**[sd.ne];
	for (long i = 0; i < sd.ne; i++)
	{
		
		sd.nnod[i] == KR_BODU ? mez = KR_POCET_PLOCH : mez = CT_POCET_PLOCH;		
		poplpr[i] = mez;
		seplpr[i] = new Plocha*[mez];
	}
	// zatim nemam o plochach prvku zadne informace
	znamePlochy = SpolecnaData::poleShort(sd.ne, 0);	
}

Oblasti::~Oblasti()
{
	for (long i = 0; i < sd.ne; i++)
	{	
		// dealokace ploch na i-tem radku
		for (long j = 0; j < this->poplpr[i]; j++)
			delete seplpr[i][j];
		// dealokace sloupcu i-teho radku
		delete [] seplpr[i];
	}	
	// dealokace radku pole
	delete [] seplpr;
	// dealokace pole poplpr
	delete [] this->poplpr;
	// dealokace oblastí
	for(size_t i = 0; i < oblasti.size();i++)
		delete oblasti[i];
	delete [] znamePlochy;
}

// pouzit algoritmus prohledavani do hloubky DFS
// na zacatku je stav vsech prvku oznacen jako FRESH
// do zasobniku se vlozi prvni prvek jeho stav se oznaci OPEN
// vstoupi se do while cyklu kde jsou nasledujici kroky:
// 1. vyjme se prvek z vrcholu zasobniku
// 2. prochazeji se sousedi prvku
//    2 a) prvni nalezeny soused se stavem FRESH je vlozen na zasobnik, jeho stav se zmeni na OPEN a preskoci se na krok 1 
//	  2 b) pokud zadny soused neni FRESH, ale existuje OPEN soused, vlozi se na zasobnik. Stav prvku se zmeni na CLOSED
//		   pokracuje se krokem 1
//	  2 c) pokud stav vsech sousedu prvku je CLOSED, bylo vyhledano oblast.Ulozi se do pole teles
//		   a algoritmus pokracuje krokem 3
// 3. jestlize existuje uzel,ktery je jeste FRESH, vlozi se na zasobnik a pokracuje se krokem 1.
//    Pokud je jiz stav vsech prvku CLOSED algoritmus konci...


void Oblasti::vyhledejOblasti()
{
	const short FRESH = -1;
	const short OPEN  =  0;
	const short CLOSE =  1;
	short * stav = SpolecnaData::poleShort(sd.ne, FRESH);
	stav[0] = OPEN;
	long openPrvek, topStack = 0;
    bool freshNalezen;
	Oblast * oblast = new Oblast(sd, poplpr, seplpr, znamePlochy);
	oblasti.push_back(oblast);		
	oblasti[0]->vlozPrvek(0);		
	while (true)
	{
		openPrvek = -1;
		freshNalezen = false;
		for (long i = 0; i < sd.nadjelel[topStack]; i++)
		{
			if (stav[ sd.adjelel[topStack][i] ] == FRESH)
			{
				stav[ sd.adjelel[topStack][i] ] = OPEN;
				topStack = sd.adjelel[topStack][i];
				oblasti[oblasti.size()-1]->vlozPrvek(sd.adjelel[topStack][i]);	
				freshNalezen = true;
				break;
			}
			if ( stav[ sd.adjelel[topStack][i] ] == OPEN ) 
				openPrvek = sd.adjelel[topStack][i]; 
		}	// for
		
		if (freshNalezen)continue;
		if (openPrvek == -1) 
		{
			
			if (najdiFreshStart(&topStack,stav,FRESH,OPEN))
			{
				oblast = new Oblast(sd,poplpr,seplpr, znamePlochy);
				oblasti.push_back(oblast);	
				oblasti[oblasti.size()-1]->vlozPrvek(topStack);	
				continue;
			}
			break;	// konec vsechno prolezeno
		}
		stav[topStack] = CLOSE;
		topStack = openPrvek;
	}	// while
	delete [] stav;
}



bool Oblasti::najdiFreshStart(long * topStack, short * stav, const short & FRESH, const short & OPEN)
{
	for(long i = 0; i < sd.ne; i++)
	{
		if (stav[i] == FRESH)
		{
			stav[i] = OPEN;		
			*topStack = i;		// ulozim cislo FRESH elementu
			return true;		// nalezen novy vychozi prvek
		}	
	}	
	return false;				// vsechny prvky jsou CLOSE
}
	


void Oblasti::hledejStycnePlochy()
{
	size_t i, k;
	set<MyPair, Less> nlist; 
	vector<long> adjel;
	deque <double *> nvector;
	set<Plocha *, Comparator>::iterator it1;
	set<Plocha *, Comparator>::iterator it2;
	//vector<long> stycnebody[5][20];
	Plocha * p1;
	Plocha * p2;
	FOR (k, oblasti.size()-1) // oblast ktera se bude porovnavat se vsemi ostatnimi 
	{
		for (it1 = oblasti[k]->mp.begin();it1 != oblasti[k]->mp.end(); it1++) // pro vsechny body oblasti
		{			
			p1 = *it1;
			if (p1->typ == STYCNA) continue;
			for (i = k+1; i < oblasti.size(); i++) // pro vsechny zbyla telesa
			{
				for (it2 = oblasti[i]->mp.begin();it2 != oblasti[i]->mp.end(); it2++) // pro vsechny body zbylych teles
				{
					p2 = *it2;
					if (p2->typ == STYCNA) continue;
					if (*p1 == *p2)
					{
						p1->typ = STYCNA;
						p2->typ = STYCNA;
						nlists.push_back(nlist);
						adjels.push_back(adjel);
						nvectors.push_back(nvector);
						projdiStycnouPlochu(p1,p2);
					}
				}	// pro vsechny body zbylych teles
			}	// pro vsechny zbyla telesa
			if (p1->typ == NEVI_SE)	
			{
				p1->typ = VNEJSI;
				znamePlochy[p1->prvek]++;
			}
		}	// oblast ktera se bude porovnavat se vsemi ostatnimi 
	}
}

void Oblasti::projdiStycnouPlochu(Plocha * p1, Plocha * p2)
{
	short i, j;
	Plocha * soused1, * soused2;
	Plocha * tmp1, * tmp2;
	bool majiStycnou;
	this->route.clear();
	routeAdd(p1,p2);
	while(!route.empty())
	{
		p2 = *(this->route.end()-1);
		p1 = *(this->route.end()-2);
		majiStycnou = false;
		FOR(i,poplpr[p1->prvek])
		{
			tmp1 = seplpr[p1->prvek][i];
			if ( tmp1->typ != VNITRNI ) continue;
			soused1 = tmp1->sousedni[0];
			if ( znamePlochy[soused1->prvek] == poplpr[soused1->prvek] )
			{
				tmp1->typ = CLOSE;
				continue;
			}
			FOR(j,poplpr[p2->prvek])
			{
				tmp2 = seplpr[p2->prvek][j];
				if (tmp2->typ != VNITRNI)continue;
				soused2 = tmp2->sousedni[0];
				if ( znamePlochy[soused2->prvek] == poplpr[soused2->prvek] )
				{
					tmp2->typ = CLOSE; 
					continue;
				}
				majiStycnou = stycnaPlocha(soused1, soused2); 
				if (majiStycnou)
				{
					tmp1->typ = CLOSE;
					tmp2->typ = CLOSE;
					soused1->typ = CLOSE;
					soused2->typ = CLOSE;
					break;
				}
			}	// FOR
			if (majiStycnou) break;
		}	// FOR
		if (majiStycnou) continue;
		route.pop_back();
		route.pop_back();
	}	//	while (!route.empty())
}



bool Oblasti::stycnaPlocha(Plocha * soused1, Plocha * soused2)
{
	short i, j;
	FOR(i,poplpr[soused1->prvek])
	{
		soused1 = seplpr[soused1->prvek][i]; 
		if ( soused1->typ != NEVI_SE ) continue;
		FOR(j,poplpr[soused1->prvek])
		{
			soused2 = seplpr[soused2->prvek][j];
			if ( soused2->typ != NEVI_SE ) continue;
			if ( *soused1 != *soused2 ) continue;
			soused1->typ = STYCNA;
			soused2->typ = STYCNA;
			routeAdd(soused1,soused2);
			if (soused1->pbp == soused2->pbp)return true;		// stejné prvky
			if (soused1->pbp < soused2->pbp) ctyrstenKrychle(soused1, soused2);
			ctyrstenKrychle(soused2, soused1);
			return true;
		}
	}	// FOR i
	return false;
}

void Oblasti::ctyrstenKrychle (Plocha * triangle, Plocha * rect)
{
	short i, j;
	Plocha * soused;
	FOR(i,poplpr[triangle->prvek])
	{
		if (seplpr[triangle->prvek][i]->typ != VNITRNI) continue;
		soused = triangle->sousedni[0];
		if (triangle->pbp != soused->pbp)continue;
		if ( znamePlochy[soused->prvek] == poplpr[soused->prvek] )
		{
			seplpr[triangle->prvek][i]->typ = CLOSE;
			continue;
		}
		FOR(j,poplpr[soused->prvek])
		{
			soused = seplpr[soused->prvek][j]; 
			if ( soused->typ != NEVI_SE ) continue;
			if ( *soused != *rect ) continue;
			soused->typ = STYCNA;
			routeAdd(rect,soused);
			return;
		}	// FOR j
	}	// FOR i
}


void Oblasti::ulozStycneBodyAPrvky(Plocha *p1, Plocha *p2)
{
	MyPair mp;
	for( short i = 0; i < KR_PLOCHA_BODU; i++)
		if (Plocha::shoda[i] > -1)
		{
			mp.node1 = p1->bodyp[i];
			mp.node2 = p2->bodyp[ Plocha::shoda[i] ];
			// kdyz tam ta dvojice jeste neni ulozit
			if (nlists[nlists.size()-1].find(mp) == nlists[nlists.size()-1].end())
				nlists[nlists.size()-1].insert(mp);
		}
	adjels[adjels.size()-1].push_back(p1->prvek);
	adjels[adjels.size()-1].push_back(p2->prvek);

}


void Oblasti::routeAdd(Plocha *p1, Plocha *p2)
{
	this->route.push_back(p1);
	this->route.push_back(p2);
	ulozStycneBodyAPrvky(p1,p2);
	normalovyVektor(p1);
}


void Oblasti::normalovyVektor(Plocha *p)
{
	short i;
	double u[3];
	double v[3];
	double * nv = new double[3];
	FOR(i,3)u[i]= sd.xyz[p->bodyp[2]][i] - sd.xyz[p->bodyp[0]][i];	
	FOR(i,3)v[i]= sd.xyz[p->bodyp[2]][i] - sd.xyz[p->bodyp[1]][i];	
	nv[0]= u[1] * v[2] - v[1] * u[2];
	nv[1]= u[2] * v[0] - v[2] * u[0];
	nv[2]= u[0] * v[1] - v[0] * u[1];
	if(!nvectors[nvectors.size()-1].size())
	{
		nvectors[nvectors.size()-1].push_back(nv);
		return;
	}
	// kontrola smeru pomoci skalarniho soucinu s poslednim vektorem
	double * nvlast = *(nvectors[nvectors.size()-1].end()-1); 
	u[0] = 0;
	FOR(i,3)u[0]+= nv[i] * nvlast[i]; 
	if (u[0] < 0)FOR(i,3)nv[i] *= -1;
	(*(nvectors.end()-1)).push_back(nv);	
}

