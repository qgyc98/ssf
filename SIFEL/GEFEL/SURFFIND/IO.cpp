#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "IO.h"

using namespace std;


long * monoFieldAlloc(long kolikPrvku, ifstream & inFile);
long ** biFieldAlloc(long kolikPrvku,long *pprNaRadku, ifstream & inFile);
void preskocNazevPole(ifstream & inFile);
void deleteBiField(long pocetRadku,long ** pole);
void deleteBiField(long pocetRadku,double ** pole);
void vypisBiPole(long **poleSeznam, long * polePocet, long pocetRadku);



IO::IO (const char* fileName)
{
	ifstream inFile ( fileName , ifstream::in );
	// kontrola vstupniho souboru
	if(!inFile) 
	{
		cout << endl << "Nepodarilo se otevrit soubor " << fileName << " !\n";
		exit (1);
	}
	long nicDulezityho;
	// POLE xyz BOD SOURADNICE X Y Z
	inFile >> nn;				// ulozim pocet bodu
	xyz = new double *[nn];
	for (long i = 0; i < nn; i++)
	{
		xyz[i] = new double[3]; 
		inFile >> nicDulezityho; // preskocim cislo radku
		for (long j = 0; j < 3; j++)	// ulozim souradnice bodu
			inFile >> xyz[i][j];
		inFile >> nicDulezityho; // vlastnost bodu se zahazuje
	}
	// pole nodes SEZNAM CISEL BODU ELEMENTU
	inFile >> ne;	// ulozim pocet elementu
	nnod = new long [ne];
	nodes = new long *[ne]; 
	long tmp;
	for (long i = 0; i < ne; i++)
	{
		inFile >> nicDulezityho;	// preskocim cislo elementu
		inFile >> nicDulezityho;	// nactu pocet vrcholu
		nodes[i] = new long [nicDulezityho];	
		nnod[i] = nicDulezityho;
           
		for (long j = 0; j < nicDulezityho; j++)	// ulozim souradnice bodu
		{	
			inFile >> tmp;
			nodes[i][j] = tmp -1;	//  pozor ve vstupnim souboru sou prvky cislovane od 1 !!
		}
		inFile >> nicDulezityho;	// preskocim vlastnost
	}
	preskocNazevPole(inFile);
	// pole nadjnodnod POCET UZLU SOUSEDICICH S UZLEM
	
	
	
	nadjnodnod = new long [nn];
	for (long i = 0; i < nn; i++)
	{
		inFile >> nicDulezityho;	// preskocim cislo bodu
		inFile >> nadjnodnod[i];
	}
	
	preskocNazevPole(inFile);
				
	// pole adjnodnod;	CISLA SOUSEDICICH BODU S BODEM 
	adjnodnod = biFieldAlloc(nn, nadjnodnod, inFile); 
	preskocNazevPole(inFile);
		
	// pole nadjelnod CISLO_BODU|POCET ELEMENTU KTERE OBSAHUJI BOD CISLO_BODU 
	nadjelnod = new long [nn];
	for (long i = 0; i < nn; i++)
	{
		inFile >> nicDulezityho;	// preskocim cislo bodu
		inFile >> nadjelnod[i];
	}
	preskocNazevPole(inFile);
	// pole adjelnod CISLO_BODU|SEZNAM ELEMENTU KTERE OBSAHUJI BOD CISLO_BODU 
	adjelnod = biFieldAlloc(nn,nadjelnod, inFile); 
	preskocNazevPole(inFile);
		
	// pole nadjelel;
	nadjelel = new long [ne];
	for (long i = 0; i < ne; i++)
	{
		inFile >> nicDulezityho;	// preskocim cislo elementu
		inFile >> nadjelel[i];
	}
	preskocNazevPole(inFile);
	
	// pole adjelel;
	adjelel = biFieldAlloc(ne,nadjelel,inFile);
	inFile.close();
		
}


IO::~IO ()
{
	delete [] nnod;
	delete [] nadjnodnod;
	deleteBiField(nn,adjnodnod);		// cisla sousedicich bodu s bodem 
	delete [] nadjelnod;
	deleteBiField(nn, adjelnod);
	delete [] nadjelel;
	deleteBiField(ne,adjelel);	
	deleteBiField(ne,nodes);		// element = seznam cisel uzlu
	deleteBiField(nn,xyz);				// bod = souradnice x,y,z
}

void deleteBiField(long pocetRadku,long ** pole)
{
	for (long i = 0; i < pocetRadku; i++)
		delete [] pole[i]; 
	delete[] pole;
}

void deleteBiField(long pocetRadku,double ** pole)
{
	for (long i = 0; i < pocetRadku; i++)
		delete [] pole[i]; 
	delete[] pole;
}

void vypisBiPole(long **poleSeznam, long * polePocet, long pocetRadku)
{
	for (long i = 0; i < pocetRadku; i++)
	{
		cout << i << "    ";
		for (long j = 0; j < polePocet[i]; j++)
			cout << poleSeznam[i][j] << " ";
		cout << endl;
	}
	
}

long * monoFieldAlloc(long kolikPrvku, ifstream & inFile)
{
	long * pole = new long [kolikPrvku];
	long preskocit;
	inFile >> preskocit;	// preskocim cislo bodu nebo elementu	
	for (long i = 0; i < kolikPrvku; i++)
		inFile >> pole[i];	
	return pole;
}


long ** biFieldAlloc(long kolikPrvku,long *pprNaRadku, ifstream & inFile){
	long **pole = new long *[kolikPrvku];
	for (long i = 0; i < kolikPrvku; i++)
		pole[i] = monoFieldAlloc(pprNaRadku[i], inFile);	
	return pole;
}


void preskocNazevPole(ifstream & inFile){
	inFile.ignore(256, '\n');	
	inFile.ignore(256, '\n');
	inFile.ignore(256, '\n');
}
	
	

