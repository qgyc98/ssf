#ifndef __OBLAST_H
#define __OBLAST_H
#include <list>
#include <set>
#include <vector>
#include "comparator.h"
#include "spolecnaData.h"



using std::set;
using std::list;
using std::vector;


 //  system cislovani sten (ploch) u krychli a ctyrstenu
 //  krychle 6 sten po 4 bodech
 const long systemKrychle [KR_POCET_PLOCH][KR_PLOCHA_BODU] = {
						 {0,3,7,4},	// 1,4,8,5
						 {1,0,4,5},	// 2,1,5,6 
						 {2,1,6,5},	// 3,2,7,6
						 {3,2,6,7},	// 4,3,7,8
						 {0,1,2,3},	// 1,2,3,4
						 {4,5,6,7}	// 5,6,7,8 
					 	 };
 
 //  ctyrsten 4 steny po 3 bodech
 const long systemCtyrsten [CT_POCET_PLOCH][CT_PLOCHA_BODU] = {
						 {0,1,2},	// 1,2,3
						 {1,2,3},	// 2,3,4
						 {0,2,3},	// 1,3,4
						 {0,1,3}	// 1,2,4
						 };

class Oblast
{
private:
	
	SpolecnaData & sd;
	
	//  pocty ploch jednotlivych prvku
	static short * poplpr;
	
	//  pole pointeru na plochy jednotlivych prvku 
	static Plocha *** seplpr;
	
	//  pocet ploch prvku u kterých znám typ. Tzn. typ neni NEVI_SE.
	static short * znamePlochy;
	
	void PlochyPrvkuSerazene(const long prvek, vector<list<long> > & plochyPrvku);
public:
	//  mnozina ploch obsahujici vnejsi plochy oblasti
	set <Plocha *, Comparator> mp;
	
	//  konstruktor
	Oblast(SpolecnaData & sd, short * poplpr, Plocha *** seplpr, short * znamePlochy);
	
	//  destruktor
	~Oblast();
	
	//  vlozi vsechny plochy prvku do mnoziny mnoziny ploch mp
	//  parametrem je cislo prvku
	void vlozPrvek(const long prvek);

};

#endif




