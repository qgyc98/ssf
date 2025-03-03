#include <iostream>
#include "oblast.h"

using namespace std;

//  pocty ploch jednotlivych prvku
short * Oblast::poplpr;

//  pole pointeru na plochy jednotlivych prvku 
Plocha *** Oblast::seplpr;

//  pocet ploch prvku u kterých znám typ. Tzn. typ neni NEVI_SE.
short * Oblast::znamePlochy;


Oblast::Oblast(SpolecnaData & sd, short * poplpr, Plocha *** seplpr, short * znamePlochy):
sd(sd)
 {
	 Oblast::poplpr = poplpr;
	 Oblast::seplpr = seplpr;
	 Oblast::znamePlochy = znamePlochy;
 }
 
 
 Oblast::~Oblast()
 {
	 // vsechny alokovane plochy budou dealokovany v destruktoru tridy Oblasti
	 // pole seplpr obsahuje pointery na všechny vytvorene instance tridy Plocha.
	 // (do mnoziny se ukladaji ukazatele, ne objekty)
	 // toto pole vytvari trida Oblasti. Je tedy jeho vlastnikem a stara se o dealokaci.
	 // tim padem se zde nic nedealokuje.
	 
 }
 
 // vytvari pole ploch prvku serazenych podle cisla bodu vzestupne
 // cislo plochy odpovida indexu rady v poli
 void Oblast::PlochyPrvkuSerazene(const long prvek, vector<list<long> > & plochyPrvku)
 {
	// pbp = pocet bodu plochy
	long pbp,rada,sloupec;
	// pocet ploch, cislo bodu
	long pp, cb;
	const long *srfsys;
	// pocet vrcholu (bodu) prvku
	if (sd.nnod[prvek] ==  KR_BODU){ 
		pp = KR_POCET_PLOCH;
		pbp = KR_PLOCHA_BODU;
		srfsys = &systemKrychle[0][0];
	}
	else { 
		pp = CT_POCET_PLOCH;
		pbp = CT_PLOCHA_BODU;
		srfsys = &systemCtyrsten[0][0];
	}
	// ili = iterator list
	vector<list<long> >::iterator ili;
	ili = plochyPrvku.begin();
	// ino = iterator node
	list<long>::iterator ino;
	rada = 0;
	while(ili != plochyPrvku.end()){
		sloupec = 0;
		cb = sd.nodes[prvek][*(srfsys + rada * pbp + sloupec)];
		ili->push_front(cb);
		ino = ili->begin();
		while(sloupec < pbp-1){
			sloupec++;
			cb = sd.nodes[prvek][*(srfsys + rada * pbp + sloupec)];
			while(ino != ili->end()){
				if (*ino > cb)
				{
					ili->insert(ino,cb);
					break;
				}
				ino++;
			}
			if (ino == ili->end())
				ili->push_back(cb);
			ino = ili->begin();
		}	
		ili++;
		rada++;
	}
 }

void Oblast::vlozPrvek(const long prvek)
{
	long * plocha;
	long * pint;
	Plocha * pinset;
	vector<list<long> > plochyPrvku(sd.nnod[prvek] == KR_BODU ? KR_POCET_PLOCH : CT_POCET_PLOCH);
	pair<set<Plocha *, Comparator >::iterator,bool> check;
	PlochyPrvkuSerazene(prvek,plochyPrvku);
	long cisloPlochy = 0;
	for (vector<list<long> >::iterator i= plochyPrvku.begin(); i!=plochyPrvku.end(); i++) 	
	{
		plocha = new long[i->size()];
		pint = &plocha[0];
		for (list<long>::iterator j= i->begin(); j!=i->end(); j++) 		
		{
			*pint = *j;
			pint++;
		}
		Plocha * p = new Plocha(sd.xyz, plocha, i->size(), prvek, cisloPlochy);
		// ulozim pointer na vnejsi plochu do pole seplpr ("SEZNAM PLOCH PRVKU")
		seplpr[p->prvek][p->plocha] = p;
		// fce insert vrati objekt pair<iterator, bool> kde false znamena ze se prvek nevlozil (duplicita)
		check = this->mp.insert(p);
		// kontrola zda byla plocha uspesne vlozena
		// pokud ne, jedna se o vnitrni plochy, promennou vnejsi nastavim na false a vlozim do seplpr pointer na druhou plochu
		if (!check.second)
		{
			// plocha ktera uz je v mnozine
			pinset = *(check.first);
			pinset->typ = VNITRNI; 
			p->typ = VNITRNI;
			pinset->vlozSousedniPlochu(p);
			p->vlozSousedniPlochu(pinset);
			// odstranim vnitrni plochu, pointer na platny objekt 
			// je ulozen v poli seplpr[pinset->prvek][pinset->plocha]
			// nalezena shoda mezi plochami se stejnym poctem uzlu -> vyhodit z mnoziny
			if (pinset->pbp == p->pbp){
				this->mp.erase(check.first);
				znamePlochy[p->prvek]++;
				znamePlochy[pinset->prvek]++;
			}
			// v mnozine je ulozena plocha steny krychle
			// vkladam plochu steny ctyrstenu, ktera s ni sousedi
			// je nutne nastavit priznak vyhodit u steny krychle na true
			// pokud uz je nastaven tak vyhodit z mnoziny
			else if (pinset->pbp > p->pbp)
			{
				znamePlochy[p->prvek]++;
				if (pinset->vyhodit){
					znamePlochy[pinset->prvek]++;				
					this->mp.erase(check.first);
				}
				else pinset->vyhodit = true;
			}
			// nejhorsi pripad. V mnozine je ulozena stena ctyrstenu (SC), ktera sousedi
			// s vkladanou stenou krychle(SK). SC se musi vyhodit z mnoziny.
			// U SK se nastavi promenna "vyhodit" na true a znovu se vlozi.
			// Mohou nastat dva pripady. Bud se to podari, tzn., ze
			// plocha druheho sousedniho ctyrstenu jeste v mnozine neni.A nebo se to opet
			// nepovede. Pak se musi vyhodit z mnoziny druha stena ctyrstenu, nastavit promennou
			// typ na VNITRNI a nastavit pointery sousednich ploch u SK a SC
			else 
			{
				// nastavim priznak vyhozeni plochy z mnoziny pri dalsi shode
				p->vyhodit = true;
				// vyhodim plochu ctyrstenu
				this->mp.erase(check.first);
				znamePlochy[pinset->prvek]++;
				// zkusim plochu krychle vlozit znovu
				check = this->mp.insert(p);
				if (!check.second)
				{	// znovu se to nepovedlo
			
					pinset = *(check.first);
					// nastavim typ plochy ctyrstenu na VNITRNI
					pinset->typ = VNITRNI;
					// nastavim promennou "sousedni" SC -> ukazuje na SK
					pinset->vlozSousedniPlochu(p);
					// nastavim promennou "sousedni" SK -> ukazuje na SC
					p->vlozSousedniPlochu(pinset);
					// vyhodim z mnoziny druhou stenu ctyrstenu
					this->mp.erase(check.first);
					znamePlochy[pinset->prvek]++;
					znamePlochy[p->prvek]++;
				}
			}
		}
		cisloPlochy++;
	}
	
}



