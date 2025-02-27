#ifndef __COMPARATOR_H
#define __COMPARATOR_H

// http://www.linuxselfhelp.com/HOWTO/C++Programming-HOWTO-17.html

#include "plocha.h"
#include <string.h>

#define _INC(i1,i2)i1++;i2++;

class Comparator
{
private:
	// funkce resi problem porovnani dvou poli, ktere maji ruzny pocet hodnot
	// a rozhoduje ktere znich je vetsi
	// vraci true pokud je p1 vetsi nez p2

        // The following memberes were commented out by TKo (29.8.2022) because comparator object for std::set must be const invokable.
        // Definition of mycmp function seems that these data meber can be replaced by local variables because they are set to zero at the beginning of each call of mycmp
	// short index1;
	// short index2;
	bool mycmp(const Plocha * p1, const Plocha * p2) const
	{
          // the following two lines were in the original code, they must be changed due to compiler restrictions on std::set comparator object
          // TKo 29.8.2022
          //index1 = 0;
          //index2 = 0; 
		short index1 = 0;
		short index2 = 0; 
		if (p1->pbp > p2->pbp)
		{
			while (index1 < p1->pbp && index2 < p2->pbp)
			{
				if (p1->bodyp[index1] > p2->bodyp[index2])  
					return (index1 == index2);
				else if ( p1->bodyp[index1] == p2->bodyp[index2] ){
					_INC(index1,index2);				
				}
				else index1++;
			}	
		}	
		else			// p1->pbp < p2->pbp
		{
			while (index1 < p1->pbp && index2 < p2->pbp)
			{
				if (p1->bodyp[index1] > p2->bodyp[index2])
					index2++;
				else if ( p1->bodyp[index1] == p2->bodyp[index2] ){
					_INC(index1,index2);				
				}
				else return (index1 < index2);
			}
		}
		return (index1 < p1->pbp);	// true -> p1 > p2
	}

public:
	// vrácená hodnota je true pokud je p2 > p1
	bool operator()(const Plocha * p1, const Plocha * p2) const
    { 		
		if (p1->pbp == p2->pbp)
			return (memcmp(p2->bodyp, p1->bodyp, p1->pbp * sizeof(long)) > 0);
		
		// pozor na poradi argumentu
		return mycmp(p2, p1);
	}
};
#endif


