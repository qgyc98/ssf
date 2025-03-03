#include "plocha.h"

using namespace std;
// konstruktor
// parametry xyz = pole souradnic, bodyp = body plochy, pbp = pocet bodu plochy
// prvek = cislo prvku
// plocha = cislo plochy


short Plocha::shoda[KR_PLOCHA_BODU];
double const ** Plocha::xyz; 

Plocha::Plocha(const double **xyz, long * bodyp, long pbp, long prvek, long plocha):
//xyz(xyz),
prvek(prvek),
plocha(plocha),
pbp(pbp),
bodyp(bodyp),
typ(NEVI_SE),
size_sousedni(0),
vyhodit(false)
{
	Plocha::xyz = xyz;
}

// destruktor
Plocha::~Plocha()
{
	if (size_sousedni > 0) delete [] sousedni;
	delete [] bodyp; 
	//std::cout << "~Plocha()\n";
}


void Plocha::vlozSousedniPlochu(Plocha * sousedni_plocha)
{
	
	if (!size_sousedni)
	{
		sousedni = new Plocha *[++size_sousedni];
		sousedni[0] = sousedni_plocha;
		return;
	}
	
	Plocha ** p = new Plocha *[++size_sousedni];
	for (short i = 0; i < size_sousedni-1; i++)
		p[i] = sousedni[i];		
	
	p[size_sousedni - 1] = sousedni_plocha;
	delete [] sousedni;
	sousedni = p;	
}


bool Plocha::operator==( Plocha & p)const
{
	if(pbp < p.pbp)
	{
		return cmp(*this, p);
	}
	return cmp(p , *this); 
}

bool Plocha::operator!=( Plocha & p)const
{
	return !(*this == p);
}

// porovnava plochy na urovni souradnic bodu
// v pripade shody vsech souradnic bodu u obou ploch vrati true
bool Plocha::cmp(const Plocha & p1,const Plocha & p2)const
{
	short h;
	FOR(h,KR_PLOCHA_BODU) shoda[h]= -1;	// -1 = neporovnano; >=0  shodny bod pro dve plochy  
	bool equal;
	for(short i = 0; i < p1.pbp; i++)
	{	
		equal = false;
		for(short j = 0; j < p2.pbp; j++)
		{
			if (shoda[j] != -1)continue;
			if ( stejneXYZ(p1.bodyp[i], p2.bodyp[j]) )
			{
				equal = true;
				shoda[i] = j;
				break;			
			}
		}	
		if(!equal) return false;
	}	
	return true;
}


inline bool Plocha::stejneXYZ(const long cb1,const long cb2)const
{
	return ( xyz[cb1][0] == xyz[cb2][0] &&
		     xyz[cb1][1] == xyz[cb2][1] && 
			 xyz[cb1][2] == xyz[cb2][2] );
}


ostream& operator<<(ostream &os, Plocha &n)
{
	for(long i = 0; i < n.pbp; i++)
		os << n.bodyp[i]+1 << " ";
	return os;
}

