#include<cstdio>
#include<iostream>
#include<fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
using namespace std;




enum Rozmery{
	_2D = 2, 
	_3D = 3
};


enum TypRezimu{
	NORMALNI = 1, // mapuji se spolecne uzly na globalni cislovani uzlu
	HRANOVY = 2, // mapuji se ocislovane spolecne hrany
	ROHOVE = 3
};
enum TypBodu{
	VNEJSI_HRANA =1,
	VNITRNI_HRANA = 5,
	VNEJSI_ROH = 6,
	VNITRNI_ROH = 7,
	VNITRNI_BOD = 8
};

enum Vlastnost{
	MAX_X = 3,
	MIN_X = 1,
	MAX_Y = 2,
	MIN_Y = 4,
	MIN_Z = 9,
	MAX_Z = 10,
	PLD_ROH = 5,
	PPD_ROH = 6,
	PLH_ROH = 7,
	PPH_ROH = 8,
	ZLD_ROH = 11,
	ZPD_ROH = 12,
	ZLH_ROH = 13,
	ZPH_ROH = 14,
	VNITRNI = 0,
	PPA_HRANA = 29,////////XPA = PRAVA
	LZ_HRANA = 27,//////
	LP_HRANA = 23,//////
	PAZ_HRANA = 25,////////XZ = ZADNI
	PAP_HRANA = 21,////////XP = PREDNI
	HL_HRANA = 30,//////// hrany horni steny
	HP_HRANA = 22,//////
	HZ_HRANA = 26,//////
	DL_HRANA = 31,//////// hrany spodni steny
	DP_HRANA = 20,//////
	DPA_HRANA = 28,////////XPA = PRAVA
	DZ_HRANA = 24//////








};


struct HranovyBod{
	TypBodu priznak;
	Vlastnost vlastnost;
	double* mapovani;
	double souradnice;
};


struct Domena{
	// realne rozmenry domeny [m]
	double rozmerX, rozmerY, rozmerZ;
	// pocet uzlu v obou osach, ktere tvori domenu
	int pocetUzluX, pocetUzluY, pocetUzluZ;
	// uzly na nich domena zacina v globalnim cislovanim
	int pocatecniUzelX, pocatecniUzelY, pocatecniUzelZ;
	// pocatecni realne souradnice uzlu
	double polohaPocatkuX,polohaPocatkuY,polohaPocatkuZ;
	// poradi domeny
	int poradi;
};

// funkce vypise informaci o domene
void vypisDomenu(Domena _domena){
	cout << "-----------------------------------------" << endl;
	cout << "Domena : " << _domena.poradi << endl;
	cout << "Rozmery : X =  " << _domena.rozmerX << ", Y =  " << _domena.rozmerY << endl;
	cout << "Pocatek realny : X =  " << _domena.polohaPocatkuX << ", Y =  " << _domena.polohaPocatkuY << endl;
	cout << "Pocet uzlu : X =  " << _domena.pocetUzluX << ", Y =  " << _domena.pocetUzluY << endl;
	cout << "Pocatecni uzly : X =  " << _domena.pocatecniUzelX << ", Y =  " << _domena.pocatecniUzelY << endl;
}


// funkce vypise mapovani
void vypisMapovani(long*** _mapovani, int _pocX, int _pocY, int _pocZ){
	cout << " Mapovani : " << _pocX << 'X' << _pocY << 'X' << _pocZ << endl;
	//getchar();
	cout << endl;
	for( int k = 0; k < _pocZ; k++){
		cout << " Zacatek roviny" << endl;
		for( int i = _pocY-1; i >= 0 ; i--){
			for( int j = 0; j< _pocX; j++){
				//cout << '.' ;k
				cout.width(4);
				cout << _mapovani[k][i][j] << " ";
			}
			cout << endl;
		}
		cout << " Konec roviny" << endl;
	}

}

// vypise mapovani pro konkretni domenu
void vypisMapovaniDomeny(Domena _domena, long** _mapovani){

	int uzlyXOd = _domena.pocatecniUzelX;
	int uzlyYOd = _domena.pocatecniUzelY;
	int uzlyXDo = _domena.pocetUzluX + uzlyXOd;
	int uzlyYDo = _domena.pocetUzluY + uzlyYOd;



	for( int i = uzlyXOd; i < uzlyXDo; i++){
		for( int j = uzlyYOd; j< uzlyYDo; j++){

			cout <<  _mapovani[j][i] << endl;
		}
	}
}



void vypisHranu(HranovyBod* _hrana, int _pocetUzlu){
	for(int i = 0; i < _pocetUzlu; i++){
		cout << _hrana[i].priznak << endl;
	}
}



string itos(int i)	// prevede int na string
{
	stringstream s;
	s << i;
	return s.str();
}

class ParGenQuad{
public:

	Domena ***poleDomen; 
	int pocetDomenX, pocetDomenY, pocetDomenZ;
	int rezim ;
	int rozmer;
	string jmenoSouboru;
	string jmenoVystSouboru;
	int pocetUzluX , pocetUzluY, pocetUzluZ ;
	long ***poleMapovani;
	HranovyBod *hranaX;
	HranovyBod *hranaY;
	HranovyBod *hranaZ;
	ParGenQuad(string _jmeno);
	~ParGenQuad();
	int ctiVstupniSoubor();
	int init();
	int zapisSoubory(); 
	void zapisUzly(ofstream &_vystupniSoubor, Domena _domena);
	void zapisElementy(ofstream &_vystupniSoubor, Domena _domena);
	void zapisMapovani(ofstream &_vystupniSoubor, Domena _domena);
	void vypisVlastnosti();
};



void ParGenQuad::vypisVlastnosti(){


	// vystup jednotlivych uzlu pro prvni domenu
	int uzlyX ;
	int uzlyY ;

	int vlastnost;
	int uzlyXOd ;
	int uzlyYOd;
	int uzlyXDo;
	int uzlyYDo ;
	Domena pomDom;	

	for( int k = 0; k < pocetUzluZ;k ++){
		for( int j = pocetUzluX-1; j >= 0 ; j--){
				//cout << "  uzel -1  " << (uzlyX -1) << endl;
			for( int i = 0; i< pocetUzluY ; i++){
				//cout << "(" << j << "," << i << "," << k << ")";
					//cout.width(5);
					vlastnost = VNITRNI;
					// vlastnosti na hranach
					if(hranaX[i].priznak == VNEJSI_HRANA){
						if(i == pocetUzluX-1)vlastnost = MAX_X;
						if(i == 0) vlastnost = MIN_X;
					}
					if(hranaY[j].priznak == VNEJSI_HRANA){
						if(j == 0) vlastnost = MIN_Y;
						if(j == pocetUzluY-1) vlastnost = MAX_Y;
					}
					if(hranaY[k].priznak == VNEJSI_HRANA){
						if(k == 0) vlastnost = MIN_Z;
						if(k == pocetUzluZ-1) vlastnost = MAX_Z;
					}
					if(((hranaY[j].priznak == VNEJSI_HRANA) && (hranaX[i].priznak == VNEJSI_HRANA))|
						((hranaY[j].priznak == VNEJSI_HRANA) && (hranaZ[k].priznak == VNEJSI_HRANA))|
						((hranaZ[k].priznak == VNEJSI_HRANA) && (hranaX[i].priznak == VNEJSI_HRANA))){		
						// vlastnosti v rohu
						if((i == 0) && (j == 0)&& (k == 0)) vlastnost = PPD_ROH;
						if((i == pocetUzluX-1) && (j == pocetUzluY-1)&& (k == 0)) vlastnost = PPH_ROH;
						if((i == 0) && (j == pocetUzluY-1)&& (k == 0)) vlastnost = PLH_ROH;
						if((i == pocetUzluX-1) && (j == 0)&& (k == 0)) vlastnost = PLD_ROH;
						if((i == 0) && (j == 0)&& (k == pocetUzluZ-1)) vlastnost = ZPD_ROH;
						if((i == pocetUzluX-1) && (j == pocetUzluY-1)&& (k == pocetUzluZ-1)) vlastnost = ZPH_ROH;
						if((i == 0) && (j == pocetUzluY-1)&& (k == pocetUzluZ-1)) vlastnost = ZLH_ROH;
						if((i == pocetUzluX-1) && (j == 0)&& (k == pocetUzluZ-1)) vlastnost = ZLD_ROH;
						if((vlastnost != PPD_ROH)&(vlastnost != PPH_ROH)&(vlastnost != PLH_ROH)&
							(vlastnost != PLD_ROH)&(vlastnost != ZPD_ROH)&(vlastnost != ZPH_ROH)&
							(vlastnost != ZLH_ROH)&(vlastnost != ZLD_ROH)){
						//		printf("%d, %d, %d\n",i,j,k);
						//PPA
								if((i == pocetDomenX-1)&& (k == 0)) vlastnost = PPA_HRANA;
						//PL
						//	if((i == 0)&& (k == 0)) vlastnost = PL_HRANA;
						//PH
							//if((j == 0)&& (k == 0)) vlastnost = PD_HRANA;
						//PD
							//if((j == pocetDomenY -1)&& (k == 0)) vlastnost = PH_HRANA;
						

						//ZPA
							//if((i == pocetDomenX-1)&& (k == pocetUzluZ-1)) vlastnost = ZPA_HRANA;
						//ZL
						//	if((i == 0)&& (k == pocetUzluZ-1)) vlastnost = ZL_HRANA;
						//ZH
							//if((j == 0)&& (k == pocetUzluZ-1)) vlastnost = ZD_HRANA;
						//ZD
							//if((j == pocetDomenY -1)&& (k == pocetUzluZ-1)) vlastnost = ZH_HRANA;



						//PAP
							if((k == 0)&& (i == pocetUzluX -1 )) vlastnost = PAP_HRANA;
						//PAZ
							if((k == pocetUzluZ - 1)&& (i == pocetUzluX -1 )) vlastnost = PAZ_HRANA;
						//PAH
							//if((j == 0)&& (i == pocetUzluX -1 )) vlastnost = PAD_HRANA;
						//PAD
						//if((j == pocetUzluY-1)&& (i == pocetUzluX -1 )) vlastnost = PAH_HRANA;
	

						//LAP
							if((k == 0)&& (i == 0)) vlastnost = LP_HRANA;
						//LZ
							if((k == pocetUzluZ - 1)&& (i == 0)) vlastnost = LZ_HRANA;
						//LH
						//	if((j == 0)&& (i == 0)) vlastnost = LD_HRANA;
						//LD
						//	if((j == pocetUzluY-1)&& (i == 0)) vlastnost = LH_HRANA;

						//HP
							if((k == 0)&& (j == pocetUzluY - 1)) vlastnost = HP_HRANA;
						//HPA
							if((k == pocetUzluZ - 1 )&& (j == pocetUzluY - 1)) vlastnost = HZ_HRANA;
						//HZ
							if((i == 0)&& (j == pocetUzluY - 1)) vlastnost = HL_HRANA;
						//HL
							if((i == pocetUzluX -1)&& (j == pocetUzluY - 1)) vlastnost = PPA_HRANA;


						//HP
							if((k == 0)&& (j == 0)) vlastnost = DP_HRANA;
						//HPA
							if((k == pocetUzluZ - 1 )&& (j == 0)) vlastnost = DZ_HRANA;
						//HZ
							if((i == 0)&& (j == 0)) vlastnost = DL_HRANA;
						//HL
							if((i == pocetUzluX -1)&& (j == 0)) vlastnost = DPA_HRANA;
					}
					}
					cout.width(4);
					
					cout << vlastnost ;
				}
			cout<< endl;
			}
		cout << "Dalsi rez" << endl << "---------------------------------"<< endl;
		}
}
// konstruktor
//vstupnim parametrem je jmeno souboru a typ rezimu
ParGenQuad::ParGenQuad(std::string _jmeno):jmenoSouboru(_jmeno){

}

ParGenQuad::~ParGenQuad(){
	delete[] hranaX;
	delete[] hranaY;

	for( int i = 0; i< pocetDomenY; i++){
		delete [] poleDomen[i];
	}
	delete[] poleDomen;


	for( int i = 0; i< pocetUzluY; i++){
		delete[] poleMapovani[i];
	}
	delete[] poleMapovani;
}
// nacte soubor a pripravi pracovni pole
int ParGenQuad::ctiVstupniSoubor(){
	ifstream vstupniSoubor;
	vstupniSoubor.open(jmenoSouboru.c_str());
	if(vstupniSoubor.fail()){
		cout << "Vstupni soubor nenalezen." << endl ;
		exit(10);
	}
	// nactu rezim 
	vstupniSoubor >> rezim;
	if (!((rezim != NORMALNI)||(rezim != HRANOVY)||(rezim != ROHOVE))){
		cout << "Chybne zadany rezim : " << rezim << endl;
		exit(20);
	}
	// nactu rozmer 
	vstupniSoubor >> rozmer;
	if (!((rozmer != _3D)||(rozmer != _2D))){
		cout << "Chybne zadany rozmer modelovani : " << rezim << endl;
		exit(20);
	}

	// nactu pocet domen
	vstupniSoubor >> pocetDomenX;
	vstupniSoubor >> pocetDomenY;
	if(rozmer == _3D){
		vstupniSoubor >> pocetDomenZ;
	}
	else{
		pocetDomenZ = 1;
	}
	cout << "Pocet domen v jednotlivych smerech(X, Y, Z) " << endl << pocetDomenX << " " << pocetDomenY << " " << pocetDomenZ << endl;


	poleDomen = new Domena**[pocetDomenZ];
	for( int j = 0; j< pocetDomenZ; j++){
		poleDomen[j] = new Domena*[pocetDomenY];
		for( int i = 0; i< pocetDomenY; i++){
			poleDomen[j][i] = new Domena[pocetDomenX];
		}
	}
	int pomPocet;
	pocetUzluX = 0;
	pocetUzluY = 0;
	pocetUzluZ = 0;
	int poradi = 1;
	double vzdalenostX = 0, vzdalenostY = 0, vzdalenostZ = 0;
	double pomRozmer;
	//cout << " Nacitam informace X ." << endl;
	//ctu X
	for( int k = 0; k< pocetDomenX; k++){
		vstupniSoubor >> pomRozmer;
		vstupniSoubor >> pomPocet;
		// nastavim hodnoty X pro vsechny domeny v rezu pro danne X
		for( int j = 0; j< pocetDomenY; j++){
			for( int i = 0; i< pocetDomenZ; i++){

				poleDomen[i][j][k].pocetUzluX = pomPocet+1;
				poleDomen[i][j][k].rozmerX = pomRozmer;
				poleDomen[i][j][k].pocatecniUzelX = pocetUzluX;
				poleDomen[i][j][k].poradi = poradi++; // cislujeme od jednicky
				poleDomen[i][j][k].polohaPocatkuX = vzdalenostX;
			}
		}
		//zvetsime pomocne promenne
		pocetUzluX += pomPocet ;
		vzdalenostX += pomRozmer;
	}
//	cout << "Informace X nacetny " << endl;

	//------------------------------------------------
	//cout << " Nacitam informace Y " << endl;
	//ctu Y
	for( int k = 0; k< pocetDomenY; k++){
		vstupniSoubor >> pomRozmer;
		vstupniSoubor >> pomPocet;
		// nastavim hodnoty Y pro vsechny domeny v rezu pro danne Y
		for( int j = 0; j< pocetDomenX; j++){
			for( int i = 0; i< pocetDomenZ; i++){

				poleDomen[i][k][j].pocetUzluY = pomPocet+1;
				poleDomen[i][k][j].rozmerY = pomRozmer;
				poleDomen[i][k][j].pocatecniUzelY = pocetUzluY;
				//poleDomen[i][k][j].poradi = (i*pocetDomenX + j)+1; // cislujeme od jednicky
				poleDomen[i][k][j].polohaPocatkuY = vzdalenostY;
			}
		}
		//zvetsime pomocne promenne
		pocetUzluY += pomPocet ;
		vzdalenostY += pomRozmer;
	}
//	cout << "Informace X nacetny " << endl;
	//------------------------------------------------
	if(rozmer == _3D){
	//	cout << " Nacitam informace Z " << endl;
	//	cout << " Pocet domen z = " << pocetDomenZ << endl;
		for( int k = 0; k< pocetDomenZ; k++){
			vstupniSoubor >> pomRozmer;
			vstupniSoubor >> pomPocet;
			// nastavim hodnoty Y pro vsechny domeny v rezu pro danne Y
			for( int j = 0; j< pocetDomenX; j++){
				for( int i = 0; i< pocetDomenY; i++){

					poleDomen[k][i][j].pocetUzluZ = pomPocet+1;
					poleDomen[k][i][j].rozmerZ = pomRozmer;
					poleDomen[k][i][j].pocatecniUzelZ = pocetUzluZ;
					//poleDomen[k][i][j].poradi = (i*pocetDomenX + j)+1; // cislujeme od jednicky
					poleDomen[k][i][j].polohaPocatkuZ = vzdalenostZ;
				}
			}
			//zvetsime pomocne promenne
			pocetUzluZ += pomPocet ;
			vzdalenostZ += pomRozmer;
		}
	//	cout << " Nacitam informace Z " << endl;
	}
	// pricteme posledni protoze ten jediny sme nepocitali
	pocetUzluX++;
	pocetUzluY++;
	pocetUzluZ++;
	return 1;
}







// provede potrebne operace a pripravi vse pro vystup
int ParGenQuad::init(){

	poleMapovani = new long**[pocetUzluZ];
	for( int i = 0; i< pocetUzluZ; i++){
		poleMapovani[i] = new long*[pocetUzluY];
		for( int j = 0; j< pocetUzluY; j++)	poleMapovani[i][j] = new long[pocetUzluX];
	}



	int pomocnaHodnota;
	for( int k = 0; k < pocetUzluZ; k++){
		for( int i = 0; i < pocetUzluY; i++){
			for( int j = 0; j< pocetUzluX; j++){
				// do pole nacpu hodnoty cisel podle poradi prvku v matici
				pomocnaHodnota = k*pocetUzluY*pocetUzluX + i * pocetUzluX + j;
				poleMapovani[k][i][j] = 0;//pomocnaHodnota ;
			}
		}
	}




	hranaX = new HranovyBod[pocetUzluX];
	hranaY = new HranovyBod[pocetUzluY];
	hranaZ = new HranovyBod[pocetUzluZ];
	double prirustek;
	double pomRozmer = 0;
	int pomPocet = 0;
	int pomSourad;

	// pro vsechny domeny zacneme pocitat mapovani X
	// smer Z neni nutne uvazovat, protoze pole priznaku je spolecne pro celou matici
	for( int i = 0; i< pocetDomenX; i++){
		// pocet uzlu ve smeru X
		pomPocet = poleDomen[0][0][i].pocetUzluX;
		// prirustek realnych souradnic mezi uzly
	
		prirustek = (double)poleDomen[0][0][i].rozmerX / (double)(pomPocet-1);
		//printf("pocet uzlu %d %f / %f%f\n",poleDomen[0][0][i].pocetUzluX ,(double)poleDomen[0][0][i].rozmerX , (double)(pomPocet-1), prirustek);
		// pro jednotlive domeny zaciname pocitat realne body uzlu
		for( int j = 0; j< pomPocet; j++){// pomPocet obsahuje pocet uzlu domeny
			pomSourad = j + poleDomen[0][0][i].pocatecniUzelX;
			if((j == 0) && (i != 0)) 	hranaX[pomSourad].priznak = VNITRNI_HRANA;
			else hranaX[pomSourad].priznak = VNITRNI_BOD;

			hranaX[pomSourad].souradnice = (prirustek*j) + poleDomen[0][0][i].polohaPocatkuX;
		}

		pomRozmer +=  poleDomen[0][0][i].rozmerX; // je nutne si pamatovat aktualni vzdalenost od pocatku
	}
	pomRozmer = 0;

	// pro vsechny domeny zacneme pocitat mapovani Y
	for( int i = 0; i< pocetDomenY; i++){
		// pocet uzlu ve smeru Y
		pomPocet = poleDomen[0][i][0].pocetUzluY;
		// prirustek realne vzdalenosti
		prirustek = poleDomen[0][i][0].rozmerY / (pomPocet-1);
			//printf("pocet uzlu %d %f / %f% f\n",poleDomen[0][i][0].pocetUzluY, (double)poleDomen[0][0][i].rozmerY , (double)(pomPocet-1), prirustek);
		// pro jednotlive domeny zaciname pocitat realne body uzlu
		for( int j = 0; j< pomPocet; j++){// pomPocet obsahuje pocet uzlu domeny
			pomSourad = j + poleDomen[0][i][0].pocatecniUzelY;
			if((j == 0) && (i != 0)) 	hranaY[pomSourad].priznak = VNITRNI_HRANA;
			else hranaY[pomSourad].priznak = VNITRNI_BOD;
			hranaY[pomSourad].souradnice = (prirustek*j) + poleDomen[0][i][0].polohaPocatkuY;
		}
		pomRozmer +=  poleDomen[0][i][0].rozmerY; // je nutne si pamatovat aktualni vzdalenost od pocatku
	}



	if(rozmer == _3D){

		// pro vsechny domeny zacneme pocitat mapovani Z
		for( int i = 0; i< pocetDomenZ; i++){
			// pocet uzlu ve smeru Z
			pomPocet = poleDomen[i][0][0].pocetUzluZ;
			// prirustek realne vzdalenosti
			prirustek = poleDomen[i][0][0].rozmerZ / (pomPocet-1);
			//	printf("pocet uzlu %d  %f / %f%f\n",poleDomen[i][0][0].pocetUzluZ, (double)poleDomen[0][0][i].rozmerZ , (double)(pomPocet-1), prirustek);
			// pro jednotlive domeny zaciname pocitat realne body uzlu
			for( int j = 0; j< pomPocet; j++){// pomPocet obsahuje pocet uzlu domeny
				pomSourad = j + poleDomen[i][0][0].pocatecniUzelZ;
				if((j == 0) && (i != 0)) 	hranaZ[pomSourad].priznak = VNITRNI_HRANA;
				else hranaZ[pomSourad].priznak = VNITRNI_BOD;
				hranaZ[pomSourad].souradnice = (prirustek*j) + poleDomen[i][0][0].polohaPocatkuZ;
			}
			pomRozmer +=  poleDomen[i][0][0].rozmerZ; // je nutne si pamatovat aktualni vzdalenost od pocatku
		}


	}






	hranaX[0].priznak = VNEJSI_HRANA;
	hranaX[pocetUzluX-1].priznak = VNEJSI_HRANA;
	hranaY[0].priznak = VNEJSI_HRANA;
	hranaY[pocetUzluY-1].priznak = VNEJSI_HRANA;
	hranaZ[0].priznak = VNEJSI_HRANA;
	hranaZ[pocetUzluZ-1].priznak = VNEJSI_HRANA;


	switch(rezim){
				case NORMALNI:
					cout << endl << "Cislovani spolecnych uzlu urcuje globalni cislovani" << endl;
					break;
				case HRANOVY:
					cout << endl << "Spolecne uzly jsou cislovany na hranach s rostoucim Y a pak s roustoucim X" << endl;
					break;
				case ROHOVE:
					cout << endl << "Spolecne uzly jsou cislovany na hranach s rostoucim Y a pak s roustoucim X" << endl
						<< ", rohove uzly se cisluji stejnym zpusobem v zapornych cislech" << endl;
					break;
				default: cout << "Spatne zadany rezim : " << rezim << endl;
	}

	int mapovaciUzel = 1;
	int mapovaciRoh = -1;
	for( int k = 0; k < pocetUzluZ; k++){
		for( int i = 0; i < pocetUzluX; i++){
			for( int j = 0; j< pocetUzluY; j++){
				//rezim = ROHOVE;
				switch(rezim){
				case NORMALNI:
					// odstraneno kvuli tomu, ze v tomto rezimu se maji vypisovat vsechny cisla uzlu v domene
					// ne jen ty co jsou na spojnici
					if(rozmer ==_3D){

						poleMapovani[k][j][i] = mapovaciUzel;//(i*pocetUzluY+ j);
						mapovaciUzel++;

					}
					else{


						poleMapovani[k][j][i] = mapovaciUzel;//(i*pocetUzluY+ j);
						mapovaciUzel++;

					}
					break;
				case HRANOVY:
					if(rozmer ==_3D){
						if((hranaX[i].priznak == VNITRNI_HRANA) || (hranaY[j].priznak == VNITRNI_HRANA)|| (hranaZ[k].priznak == VNITRNI_HRANA)){
							poleMapovani[k][j][i] = mapovaciUzel;
							mapovaciUzel++;
						}
					}
					else{
						if((hranaX[i].priznak == VNITRNI_HRANA) || (hranaY[j].priznak == VNITRNI_HRANA)){
							poleMapovani[k][j][i] = mapovaciUzel;
							mapovaciUzel++;
						}
					}


					break;
				case ROHOVE:
					if(rozmer ==_3D){

						if(((hranaX[i].priznak == VNITRNI_HRANA) && (hranaY[j].priznak == VNITRNI_HRANA)&& (hranaZ[k].priznak == VNITRNI_HRANA))||
							((hranaX[i].priznak == VNEJSI_HRANA) && (hranaY[j].priznak == VNITRNI_HRANA)&& (hranaZ[k].priznak == VNITRNI_HRANA))||
							((hranaX[i].priznak == VNITRNI_HRANA) && (hranaY[j].priznak == VNEJSI_HRANA)&& (hranaZ[k].priznak == VNITRNI_HRANA))||
							((hranaX[i].priznak == VNITRNI_HRANA) && (hranaY[j].priznak == VNITRNI_HRANA)&& (hranaZ[k].priznak == VNEJSI_HRANA))||
							((hranaX[i].priznak == VNEJSI_HRANA) && (hranaY[j].priznak == VNEJSI_HRANA)&& (hranaZ[k].priznak == VNITRNI_HRANA))||
							((hranaX[i].priznak == VNITRNI_HRANA) && (hranaY[j].priznak == VNEJSI_HRANA)&& (hranaZ[k].priznak == VNEJSI_HRANA))||
							((hranaX[i].priznak == VNEJSI_HRANA) && (hranaY[j].priznak == VNITRNI_HRANA)&& (hranaZ[k].priznak == VNEJSI_HRANA))||
							((hranaX[i].priznak == VNITRNI_HRANA) && (hranaY[j].priznak == VNEJSI_HRANA)))
						{
							poleMapovani[k][j][i] = mapovaciRoh;
							mapovaciRoh--;
						}
						else if((hranaX[i].priznak == VNITRNI_HRANA) || (hranaY[j].priznak == VNITRNI_HRANA)|| (hranaZ[k].priznak == VNITRNI_HRANA)){
							poleMapovani[k][j][i] = mapovaciUzel;
							mapovaciUzel++;
						}
					}
					else{
						if(((hranaX[i].priznak == VNITRNI_HRANA) && (hranaY[j].priznak == VNITRNI_HRANA))||
							((hranaX[i].priznak == VNEJSI_HRANA) && (hranaY[j].priznak == VNITRNI_HRANA))||
							((hranaX[i].priznak == VNITRNI_HRANA) && (hranaY[j].priznak == VNEJSI_HRANA)))
						{
							poleMapovani[k][j][i] = mapovaciRoh;
							mapovaciRoh--;
						}
						else{
							if((hranaX[i].priznak == VNITRNI_HRANA) || (hranaY[j].priznak == VNITRNI_HRANA))
							{
								poleMapovani[k][j][i] = mapovaciUzel;
								mapovaciUzel++;
							}
						}
					}

					break;
				default: cout << "Spatne zadany rezim : " << rezim << endl;

					}

				}
			}
		}
	

	////vypisMapovani(poleMapovani, pocetUzluX, pocetUzluY);
	////getchar();
	return 1;
}





int ParGenQuad::zapisSoubory(){
	ofstream vystup;
	if(rozmer ==_3D){
		for( int k = 0; k< pocetDomenZ; k++){
			for( int i = 0; i< pocetDomenY; i++){

				for( int j = 0; j< pocetDomenX; j++){
					Domena aktDomena = poleDomen[k][i][j];
					// otevreme soubor

					int pozicePripony = jmenoSouboru.find(".con");
					jmenoSouboru = jmenoSouboru.substr(0, pozicePripony);
					jmenoVystSouboru = jmenoSouboru + itos(aktDomena.poradi) + ".top";
					//cout << aktDomena.poradi<< " " << endl;

					vystup.open(jmenoVystSouboru.c_str());
					//zapiseme jednotlive casti vystupu
					zapisUzly(vystup, aktDomena);
					zapisElementy(vystup, aktDomena);
					zapisMapovani(vystup, aktDomena);
					vystup.close();
				}
			}
		}
	}
	else {

		for( int i = 0; i< pocetDomenY; i++){

			for( int j = 0; j< pocetDomenX; j++){
				Domena aktDomena = poleDomen[0][i][j];
				// otevreme soubor

				int pozicePripony = jmenoSouboru.find(".con");
				jmenoSouboru = jmenoSouboru.substr(0, pozicePripony);
				jmenoVystSouboru = jmenoSouboru + itos(aktDomena.poradi) + ".top";
				//cout << aktDomena.poradi<< " " << endl;

				vystup.open(jmenoVystSouboru.c_str());
				//zapiseme jednotlive casti vystupu
				zapisUzly(vystup, aktDomena);
				zapisElementy(vystup, aktDomena);
				zapisMapovani(vystup, aktDomena);
				vystup.close();
			}
		}
	}
	////vypisMapovani(poleMapovani, pocetUzluX, pocetUzluY);
	return 1;
}

void ParGenQuad::zapisUzly(ofstream &_vystupniSoubor, Domena _domena){
	int vlastnost;
	int uzlyX = _domena.pocetUzluX;
	int uzlyY = _domena.pocetUzluY;
	int uzlyXOd = _domena.pocatecniUzelX;
	int uzlyYOd = _domena.pocatecniUzelY;
	int uzlyXDo = _domena.pocetUzluX + uzlyXOd;
	int uzlyYDo = _domena.pocetUzluY + uzlyYOd;

	if(rozmer ==_3D){

		int uzlyZ = _domena.pocetUzluZ;
		int uzlyZOd = _domena.pocatecniUzelZ;
		int uzlyZDo = _domena.pocetUzluZ + uzlyZOd;
		_vystupniSoubor << uzlyX*uzlyY*uzlyZ << endl;
		for( int k = uzlyZOd; k < uzlyZDo;k ++){
			for( int i = uzlyXOd; i < uzlyXDo; i++){
				//cout << "  uzel -1  " << (uzlyX -1) << endl;
				for( int j = uzlyYOd; j< uzlyYDo; j++){

					//cout.width(5);
					vlastnost = VNITRNI;
					// vlastnosti na hranach
					if(hranaX[i].priznak == VNEJSI_HRANA){
						if(i == pocetUzluX-1)vlastnost = MAX_X;
						if(i == 0) vlastnost = MIN_X;
					}
					if(hranaY[j].priznak == VNEJSI_HRANA){
						if(j == 0) vlastnost = MIN_Y;
						if(j == pocetUzluY-1) vlastnost = MAX_Y;
					}
					if(hranaY[k].priznak == VNEJSI_HRANA){
						if(k == 0) vlastnost = MIN_Z;
						if(k == pocetUzluZ-1) vlastnost = MAX_Z;
					}
					if(((hranaY[j].priznak == VNEJSI_HRANA) && (hranaX[i].priznak == VNEJSI_HRANA))|
						((hranaY[j].priznak == VNEJSI_HRANA) && (hranaZ[k].priznak == VNEJSI_HRANA))|
						((hranaZ[k].priznak == VNEJSI_HRANA) && (hranaX[i].priznak == VNEJSI_HRANA))){		
						// vlastnosti v rohu
						if((i == 0) && (j == 0)&& (k == 0)) vlastnost = PPD_ROH;
						if((i == pocetUzluX-1) && (j == pocetUzluY-1)&& (k == 0)) vlastnost = PPH_ROH;
						if((i == 0) && (j == pocetUzluY-1)&& (k == 0)) vlastnost = PLH_ROH;
						if((i == pocetUzluX-1) && (j == 0)&& (k == 0)) vlastnost = PLD_ROH;
						if((i == 0) && (j == 0)&& (k == pocetUzluZ-1)) vlastnost = ZPD_ROH;
						if((i == pocetUzluX-1) && (j == pocetUzluY-1)&& (k == pocetUzluZ-1)) vlastnost = ZPH_ROH;
						if((i == 0) && (j == pocetUzluY-1)&& (k == pocetUzluZ-1)) vlastnost = ZLH_ROH;
						if((i == pocetUzluX-1) && (j == 0)&& (k == pocetUzluZ-1)) vlastnost = ZLD_ROH;
						if((vlastnost != PPD_ROH)&(vlastnost != PPH_ROH)&(vlastnost != PLH_ROH)&
							(vlastnost != PLD_ROH)&(vlastnost != ZPD_ROH)&(vlastnost != ZPH_ROH)&
							(vlastnost != ZLH_ROH)&(vlastnost != ZLD_ROH)){
								//printf("%d, %d, %d\n",i,j,k);
						//PPA
								if((i == pocetDomenX-1)&& (k == 0)) vlastnost = PPA_HRANA;
						//PL
							//if((i == 0)&& (k == 0)) vlastnost = PL_HRANA;
						//PH
						//	if((j == 0)&& (k == 0)) vlastnost = PD_HRANA;
						//PD
						//	if((j == pocetDomenY -1)&& (k == 0)) vlastnost = PH_HRANA;
						

						//ZPA
						//	if((i == pocetDomenX-1)&& (k == pocetUzluZ-1)) vlastnost = ZPA_HRANA;
						//ZL
						//	if((i == 0)&& (k == pocetUzluZ-1)) vlastnost = ZL_HRANA;
						//ZH
						//	if((j == 0)&& (k == pocetUzluZ-1)) vlastnost = ZD_HRANA;
						//ZD
						//	if((j == pocetDomenY -1)&& (k == pocetUzluZ-1)) vlastnost = ZH_HRANA;



						//PAP
							if((k == 0)&& (i == pocetUzluX -1 )) vlastnost = PAP_HRANA;
						//PAZ
							if((k == pocetUzluZ - 1)&& (i == pocetUzluX -1 )) vlastnost = PAZ_HRANA;
						//PAH
						//	if((j == 0)&& (i == pocetUzluX -1 )) vlastnost = PAD_HRANA;
						//PAD
						//	if((j == pocetUzluY-1)&& (i == pocetUzluX -1 )) vlastnost = PAH_HRANA;
	

						//LAP
							if((k == 0)&& (i == 0)) vlastnost = LP_HRANA;
						//LZ
							if((k == pocetUzluZ - 1)&& (i == 0)) vlastnost = LZ_HRANA;
						//LH
						//	if((j == 0)&& (i == 0)) vlastnost = LD_HRANA;
						//LD
						//	if((j == pocetUzluY-1)&& (i == 0)) vlastnost = LH_HRANA;

						//HP
							if((k == 0)&& (j == pocetUzluY - 1)) vlastnost = HP_HRANA;
						//HPA
							if((k == pocetUzluZ - 1 )&& (j == pocetUzluY - 1)) vlastnost = HZ_HRANA;
						//HZ
							if((i == 0)&& (j == pocetUzluY - 1)) vlastnost = HL_HRANA;
						//HL
							if((i == pocetUzluX -1)&& (j == pocetUzluY - 1)) vlastnost = PPA_HRANA;


						//HP
							if((k == 0)&& (j == 0)) vlastnost = DP_HRANA;
						//HPA
							if((k == pocetUzluZ - 1 )&& (j == 0)) vlastnost = DZ_HRANA;
						//HZ
							if((i == 0)&& (j == 0)) vlastnost = DL_HRANA;
						//HL
							if((i == pocetUzluX -1)&& (j == 0)) vlastnost = DPA_HRANA;
					}

					
					
					}

					_vystupniSoubor.precision(12);
					
					_vystupniSoubor <<fixed <<hranaX[i].souradnice << " " << hranaY[j].souradnice << " " << hranaZ[k].souradnice << " " << vlastnost << endl;
				}
			}
		}
	}
	else {
		_vystupniSoubor << uzlyX*uzlyY << endl;
		// vystup jednotlivych uzlu pro prvni domenu

		for( int i = uzlyXOd; i < uzlyXDo; i++){

			for( int j = uzlyYOd; j< uzlyYDo; j++){
				vlastnost = VNITRNI;
				// vlastnosti na hranach
				if(hranaX[i].priznak == VNEJSI_HRANA){
					if(i == pocetUzluX-1){
						vlastnost = MAX_X;
					}
					if(i == 0) vlastnost = MIN_X;
				}
				if(hranaY[j].priznak == VNEJSI_HRANA){
					if(j == 0) vlastnost = MIN_Y;
					if(j == pocetUzluY-1) vlastnost = MAX_Y;
				}
				if((hranaY[j].priznak == VNEJSI_HRANA) && (hranaX[i].priznak == VNEJSI_HRANA)){		
					// vlastnosti v rohu
					if((i == 0) && (j == 0)) vlastnost = PPD_ROH;
					if((i == pocetUzluX-1) && (j == pocetUzluY-1)) vlastnost = PPH_ROH;
					if((i == 0) && (j == pocetUzluY-1)) vlastnost = PLH_ROH;
					if((i == pocetUzluX-1) && (j == 0)) vlastnost = PLD_ROH;
				}
				_vystupniSoubor.precision(12);
				
				_vystupniSoubor << fixed << hranaX[i].souradnice << " " << hranaY[j].souradnice << " " << 0.0l << " " << vlastnost << endl;
			}
		}
	}
}
void ParGenQuad::zapisElementy(ofstream &_vystupniSoubor, Domena _domena){
	int uzlyX = _domena.pocetUzluX;
	int uzlyY = _domena.pocetUzluY;
	int uzlyXOd = _domena.pocatecniUzelX;
	int uzlyYOd = _domena.pocatecniUzelY;
	int uzlyXDo = _domena.pocetUzluX + uzlyXOd;
	int uzlyYDo = _domena.pocetUzluY + uzlyYOd;
	int pomPocetY = _domena.pocetUzluY;
	if(rozmer ==_3D){
		int pomPocetZ = _domena.pocetUzluZ;
		int uzlyZ = _domena.pocetUzluZ;
		int uzlyZOd = _domena.pocatecniUzelZ;
		int uzlyZDo = _domena.pocetUzluZ + uzlyZOd;
		for( int k = 0; k < uzlyZ-1; k++){
			for( int i = 0; i < uzlyX-1; i++){
				for( int j = 0; j< uzlyY-1; j++){
					_vystupniSoubor 
						<< "8 " << ((k*pomPocetY*pomPocetZ)+(i*pomPocetY)+(j))+1
						<< " "<< ((k*pomPocetY*pomPocetZ)+((i+1)*pomPocetY)+(j))+1
						<< " "<< ((k*pomPocetY*pomPocetZ)+((i+1)*pomPocetY)+(j+1))+1
						<< " "<< ((k*pomPocetY*pomPocetZ)+(i*pomPocetY)+(j+1))+1
						<< " "<< (((k+1)*pomPocetY*pomPocetZ)+(i*pomPocetY)+(j))+1
						<< " "<< (((k+1)*pomPocetY*pomPocetZ)+((i+1)*pomPocetY)+(j))+1
						<< " "<< (((k+1)*pomPocetY*pomPocetZ)+((i+1)*pomPocetY)+(j+1))+1
						<< " "<< (((k+1)*pomPocetY*pomPocetZ)+(i*pomPocetY)+(j+1))+1
						<< " " << endl;
				}
			}
		}

	}
	else {


		for( int i = 0; i < uzlyX-1; i++){
			for( int j = 0; j< uzlyY-1; j++){

				_vystupniSoubor << "4 " << (i*pomPocetY + j)+1 << " "<< ((i+1)*pomPocetY + j)+1 << " "<< ((i+1)*pomPocetY + (j+1))+1 << " "<< (i*pomPocetY + (j+1))+1
					<< " " << endl;
			}
		}
	}
}
void ParGenQuad::zapisMapovani(ofstream &_vystupniSoubor, Domena _domena)
{
	int uzlyXOd = _domena.pocatecniUzelX;
	int uzlyYOd = _domena.pocatecniUzelY;
	int uzlyXDo = _domena.pocetUzluX + uzlyXOd;
	int uzlyYDo = _domena.pocetUzluY + uzlyYOd;
	if( rozmer == _3D){
		int uzlyZOd = _domena.pocatecniUzelZ;
		int uzlyZDo = _domena.pocetUzluZ + uzlyZOd;
		for( int k = uzlyZOd; k < uzlyZDo; k++){
			for( int i = uzlyXOd; i < uzlyXDo; i++){
				for( int j = uzlyYOd; j< uzlyYDo; j++){

					_vystupniSoubor <<  poleMapovani[k][j][i] << endl;
				}
			}
		}
	} 
	else {
		for( int i = uzlyXOd; i < uzlyXDo; i++){
			for( int j = uzlyYOd; j< uzlyYDo; j++){
				_vystupniSoubor <<  poleMapovani[0][j][i] << endl;
			}
		}
	}
}
int main (int argc,char* argv[]){
	if( argc < 2) {
		cout << " maly pocet parametru " << endl;
		exit(30);
	}

	ParGenQuad *generator = new ParGenQuad( argv[1]);
	generator->ctiVstupniSoubor();
	generator->init();
	generator->zapisSoubory();
	//vypisMapovani(generator->poleMapovani, generator->pocetUzluX, generator->pocetUzluY, generator->pocetUzluZ);
	//generator->vypisVlastnosti();
	//getchar();
	return 1;
}
