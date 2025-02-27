#include "IO.h"
#include "surffind.h"
#include <iostream>
#include <time.h>

using std::cout;

int main (int argc, char * argv[])
{

	IO io("kvadr.top"); // vytvori jednotliva vstupni pole jako argumenty slouzici fci stycne plochy
	// *** TESTOVÁNÍ RYCHLOSTI PROGRAMU V ZÁVISLOTI NA ZMÌNÌ ALGORITMU***
/*	
	const long POCET_OPAKOVANI_SPUSTENI = 5;
	const long POCET_SPUSTENI_PROGRAMU  = 500;
	double prumer = 0;
	double t;
	for (long j = 0; j < POCET_OPAKOVANI_SPUSTENI; j++)
	{
		clock_t t1=clock();
		for(long i = 0; i < POCET_SPUSTENI_PROGRAMU;i++)
		{
	
*/	
	SurfFind sf( (const long)io.nn, (const long)io.ne, (const long*)io.nadjelel, 
				 (const long**)io.adjelel, (const long*)io.nnod, 
				 (const long**)io.nodes, (const double**)io.xyz );	
	
	cout << sf;
	
/*

	}
		clock_t t2=clock();
		t = (t2 - t1);
		prumer += t;
		cout << POCET_SPUSTENI_PROGRAMU << " vykonani trvalo : " << t / CLOCKS_PER_SEC << " s\n";
		
	}
	prumer /= ( CLOCKS_PER_SEC * POCET_OPAKOVANI_SPUSTENI );
	cout	<< '\n' << POCET_SPUSTENI_PROGRAMU  
			<< "x vykonat tento priklad zabere prumerne " 
			<< prumer << " s.\n";
	cout << "Prumerny cas jednoho vykonani je " << prumer / POCET_SPUSTENI_PROGRAMU << " s\n"; 
	
*/	
	return 0;	

}	
