#ifndef __PLOCHA_H
#define __PLOCHA_H

#include <iostream>

// celkovy pocet bodu krychle
#define KR_BODU 8
// pocet sten krychle
#define KR_POCET_PLOCH 6
// pocet bodu steny krychle
#define KR_PLOCHA_BODU 4
// celkovy pocet bodu ctyrstenu
#define CT_BODU 4
//pocet bodu sten ctyrstenu
#define CT_POCET_PLOCH 4
//pocet bodu steny ctyrstenu
#define CT_PLOCHA_BODU 3

#define FOR(i,MAX) for(i=0;i<MAX;i++)


enum typPlochy { VNEJSI, VNITRNI, STYCNA, NEVI_SE, CLOSE };

class Plocha
{
	// porovnává souøadnice x, y, z pro dva body. V pøípadì shody vrací true jinak false
	// vstupní argumenty jsou èísla bodù
	inline bool stejneXYZ(const long cb1,const long cb2)const;
	// porovnává dvì pole èísel bodù tvoøící plochy ( promenna bodyp) dvou rùzných ploch
	// porovnání je úspìšné pokud obì plochy tvoøí xyz se stejnými souøadnicemi
	bool cmp(const Plocha & p1,const Plocha & p2)const;
public:
	// pole udrzujici informace o poslednim porovnavani dvou ploch
	// sloupec pole shoda je index v p1 a hodnota je index v p2
	static short shoda[KR_PLOCHA_BODU];
	// ukazatel na dvojrozmìrné pole souøadnic. Øádkový index je èíslo bodu
	// sloupcový index jsou po øadì souøadnice x, y, z
	static double const** xyz;
	// cislo prvku
	long prvek;
	// cislo plochy
	long plocha;
	// poèet bodù plochy... krychle - 4, ctyruhelnik - 3
	long pbp;
	// pole cisel bodu tvorici plochu
	long * bodyp;
	// informace o tom zda je plocha vnejsi vnitrni stycna nebo se to zatim nevi
	typPlochy typ;
	// pocet ukazatelu ulozenych v poli sousedni
	short size_sousedni;
	// pole ukazatelu na plochy sousednich prvku. Platne pouze pro vnitrni plochy !
	Plocha ** sousedni;
	// plocha krychle sousedi se dvema ctyrstenama
	bool vyhodit;
	// konstruktor
	Plocha(const double **xyz, long * bodyp, long plochaBodu, long prvek, long plocha);
	// porovnává plochy, využívá k tomu funkci porovnej
	bool operator ==( Plocha &)const;
	bool operator !=( Plocha &)const; 
	~Plocha();
	void vlozSousedniPlochu(Plocha * p);
	// výpis èísel bodù plochy
	friend std::ostream& operator<<(std::ostream &os, Plocha &n);
};



#endif


