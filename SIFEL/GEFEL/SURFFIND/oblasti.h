#ifndef __OBLASTI_H
#define __OBLASTI_H

#include <deque>
#include <vector>
#include "oblast.h"
#include "mypair.h" 

using std::deque;
using std::vector;


class Oblasti
{
	bool najdiFreshStart(long * topStack, short * stav, const short & FRESH, const short & OPEN);
	
	void ctyrstenKrychle (Plocha * triangle,Plocha * rect);
	
	bool stycnaPlocha(Plocha * soused1, Plocha * soused2);
	
	void ulozStycneBodyAPrvky(Plocha * sp1, Plocha * sp2); 
	
	//  projde celou stycnou plochu, parametrem jsou dva ukazatele
	//  na plochy ze dvou oblasti, ktere tvori stycnou plochu
	void projdiStycnouPlochu(Plocha * ps1, Plocha * ps2 );
	
	// prida do route p1 a p2, zavola fci ulozStycneBodyAPrvky
	void routeAdd(Plocha *p1, Plocha *p2);

	//  spocita normalovy vektor na plochu a ulozi ho do nvectors
	void normalovyVektor(Plocha *p);
	
	// jednotlive oblasti(telesa)
	vector<Oblast*> oblasti;
	
	//  pocty alokovanych ploch prvku v poli seplpr
	short * poplpr;
	
	//  seznamy podezrelych (vnejsich, stycnych) ploch prvku
	Plocha *** seplpr;
	
	SpolecnaData & sd;
	
	// prosle stycne plochy pri jejim prohledavani
	deque <Plocha *> route;

	//  pocet ploch prvku u kterých znám typ. Tzn. typ neni NEVI_SE.
	short * znamePlochy;
public:
	//  konstruktor
	Oblasti(SpolecnaData & sd);	
	
	//  destruktor
	~Oblasti();
	
	// rozdeleni vstupnich dat na autonomni oblasti
	void vyhledejOblasti();    
	
	//  hledani shodnych ploch mezi oblastmi hrubou silou
	void hledejStycnePlochy();
	
	//  lists of nodes on surface
	vector< set <MyPair, Less> > nlists;	
	
	//  numbers of adjacent elements
	vector < vector<long> > adjels;		

	// normal vectors
	deque < deque<double *> > nvectors;
};



#endif




