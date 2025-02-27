#include <stdio.h>
#include <stdlib.h>
#include "siftop.h"
#include "meshgen2d.h"

/*
Funkce vytvori elementarni prvek v souradnicovem systemu ksi,eta a pote ho prevede do souradnicoveho systemu x,y podle zadanych krajnich bodu makroelementu.
Dale vytvori pole edgenod, kde jsou ulozeny krajni uzly jednotlivych elementu

neh    .... deleni v horizontalnim smeru (deleni pro hranu 0 a 2)
nev    .... deleni ve vertikalnim smeru (deleni pro hranu 1 a 3)
nodes  .... pole krajnich uzlu v makroelementu
non    .... indexy zadanych uzlu v makroelementu
x,y    .... pole vyslednych souradnic cele topologie
            u pole *x a *y predpokladam ze se do nich budou vkladat souradnice postupne ze vsech makroelementu a proto vstupuji do funkce jako parametr
j      .... udava cislo uzlu, od ktereho se bude generovat
edges  .... pole hran
edenod .... trojrozmerne pole ktere uklada uzly na j-té hrane k-tého elementu
k      .... parametr udava, jaky element se generuje (treti element -> k=2)
el     .... struktura s generovanym prvkem
*/
void gennodes(long neh, long nev, snode *nodes, long *non, snode *uzly, long &j, long **edges, long ***edgenod,long k, selement &el)
{
  long i,nnn,m,c=0,d=0, l, ll, nvp, nid;
  double dksi,deta;
  double a=0,b=0,n;
  double *ksi,*eta;
  double e,f,g,h;       //  Interpolacni funkce pro dvourozmerny prvek, e je pro pravy horni vrchol; f,g,h jsou vrcholy postupne proti smeru hod.

  
  
  if (edges[0][0]==1) c++; if (edges[2][0]==1) c++;
  if (edges[1][0]==1) d++; if (edges[3][0]==1) d++; 
  
  nnn = (neh+1)*(nev+1);
  nnn = nnn - c*(neh+1) - d*(nev+1) +  c*d;

  if (edges[0][0]==1) copy_similar_edgenod (edgenod,edges,neh,k,0);
  if (edges[1][0]==1) copy_similar_edgenod (edgenod,edges,nev,k,1);
  if (edges[2][0]==1) copy_similar_edgenod (edgenod,edges,neh,k,2);
  if (edges[3][0]==1) copy_similar_edgenod (edgenod,edges,nev,k,3);

  if (edges[1][0]==0)           // pokud se hrana neshoduje, vygeneruje uzly na hrane
  {
    if (edges[2][0]==1) {for (i=1;i<(nev+1);i++) edgenod[k][1][i]=j+i-1; edgenod[k][1][0]=edgenod[k][2][0];}
    else    for (i=0;i<(nev+1);i++) edgenod[k][1][i]=j+i;
    if (edges[0][0]==1) edgenod[k][1][nev]=edgenod[k][0][0];
  }                 
  if (edges[3][0]==0)           // pokud se hrana neshoduje, vygeneruje uzly na hrane
  {
    if (edges[0][0]==1) if (edges[2][0]==1) {for (i=1;i<(nev+1);i++) edgenod[k][3][i]=j+nnn-nev+i;   edgenod[k][3][0]=edgenod[k][2][neh];}
                        else   for (i=0;i<(nev+1);i++) edgenod[k][3][i]=j+nnn-nev+i;
    else                if (edges[2][0]==1) {for (i=1;i<(nev+1);i++) edgenod[k][3][i]=j+nnn-1-nev+i; edgenod[k][3][0]=edgenod[k][2][neh];}
                        else   for (i=0;i<(nev+1);i++) edgenod[k][3][i]=j+nnn-1-nev+i;
    if (edges[0][0]==1) edgenod[k][3][nev]=edgenod[k][0][neh];
  }

  ksi = new double [nnn];
  eta = new double [nnn];

  for (i=0;i<nnn;i++) ksi[i]=0.0;
  for (i=0;i<nnn;i++) eta[i]=0.0;

  dksi = (2.0)/neh;
  deta = (2.0)/nev;
  a=-1.0; b=-1.0; m=1;

  if (edges[1][0]==1) a = a + dksi;
  if (edges[2][0]==1) b = b + deta;
  
  for (i=0;i<nnn;i++)
  {
    ksi[i]=a;
    eta[i]=b;

    n = m*(nev+1)-1;   // pocitadlo posledniho prvku ve sloupci

    if (edges[0][0]==1) 
      n = n - m;
    else{                     // pokud se hrana neshoduje, vygeneruje uzly na hrane
      if (edges[2][0]==1)  
        if (edges[1][0]==1){
          edgenod[k][0][m]=j+n-m; 
          edgenod[k][0][0]=edgenod[k][1][nev];
        } 
        else 
          edgenod[k][0][m-1]=j+n-m;
      else   
        if (edges[1][0]==1){
          edgenod[k][0][m]=j+n;   
          edgenod[k][0][0]=edgenod[k][1][nev];
        } 
        else edgenod[k][0][m-1]=j+n;
      if (edges[3][0]==1)  
        edgenod[k][0][neh]=edgenod[k][3][nev];
    }
    if (edges[2][0]==1)
      n = n - m;
    else{                     // pokud se hrana neshoduje, vygeneruje uzly na hrane
      if (edges[0][0]==1)  
        if (edges[1][0]==1){
          edgenod[k][2][m]=j+n-nev+1; 
          edgenod[k][2][0]=edgenod[k][1][0];
        } 
        else 
          edgenod[k][2][m-1]=j+n-nev+1;
      else
        if (edges[1][0]==1){
          edgenod[k][2][m]=j+n-nev;   
          edgenod[k][2][0]=edgenod[k][1][0];
        } 
        else 
          edgenod[k][2][m-1]=j+n-nev;
      if (edges[3][0]==1)  edgenod[k][2][neh]=edgenod[k][3][0];
                   }

    if (i==n) {
      if (edges[2][0]==1) 
        b=-1.0+deta; 
      else 
        b=-1.0; 
      a=a+dksi; 
      m++;
    } 
    else 
      b=b+deta;
  }

  for (i=0;i<nnn;i++)
  {
    e = 0.25*(1+ksi[i])*(1+eta[i]);
    f = 0.25*(1-ksi[i])*(1+eta[i]);
    g = 0.25*(1-ksi[i])*(1-eta[i]);
    h = 0.25*(1+ksi[i])*(1-eta[i]);

    uzly[j+i].x = e*nodes[non[0]].x + f*nodes[non[1]].x + g*nodes[non[2]].x + h*nodes[non[3]].x;
    uzly[j+i].y = e*nodes[non[0]].y + f*nodes[non[1]].y + g*nodes[non[2]].y + h*nodes[non[3]].y;

    uzly[j+i].nprop = 2; // surface property + volume property

    // detect edge properties
    if (ksi[i] == 1.0)
      uzly[j+i].nprop++;
    if (ksi[i] == -1.0)
      uzly[j+i].nprop++;
    if (eta[i] == 1.0)
      uzly[j+i].nprop++;
    if (eta[i] == -1.0)
      uzly[j+i].nprop++;

    // detect vertex properties
    nvp = 0L;
    if (ksi[i] == 1.0 && eta[i] == 1.0){
      nvp = nodes[non[0]].searchpropent(evertex);
      nid = non[0];
      uzly[j+i].nprop += nvp;
    }
    if (ksi[i] == -1.0 && eta[i] == 1.0){
      nvp = nodes[non[1]].searchpropent(evertex);
      nid = non[1];
      uzly[j+i].nprop += nvp;
    }
    if (ksi[i] == -1.0 && eta[i] == -1.0){
      nvp = nodes[non[2]].searchpropent(evertex);
      nid = non[2];
      uzly[j+i].nprop += nvp;
    }
    if (ksi[i] == 1.0 && eta[i] == -1.0){
      nvp = nodes[non[3]].searchpropent(evertex);
      nid = non[3];
      uzly[j+i].nprop += nvp;
    }
    // allocate property arrays at node
    uzly[j+i].entid = new entityp[uzly[j+i].nprop];
    uzly[j+i].prop = new long[uzly[j+i].nprop];
    // vertex property assignment
    ll = 0L;
    if (nvp){
      for (l=0L; l<nodes[nid].nprop; l++){
        if (nodes[nid].entid[l] == evertex){
          uzly[j+i].entid[ll] = evertex;
          uzly[j+i].prop[ll]  = nodes[nid].prop[l];
          ll++;
        }
      }
    }
    // edge property assignment
    if (ksi[i] == 1.0){
      uzly[j+i].entid[ll] = ecurve;
      uzly[j+i].prop[ll] = el.propedg[3];
      ll++;
    }
    if (ksi[i] == -1.0){
      uzly[j+i].entid[ll] = ecurve;
      uzly[j+i].prop[ll] = el.propedg[1];
      ll++;
    }
    if (eta[i] == 1.0){
      uzly[j+i].entid[ll] = ecurve;
      uzly[j+i].prop[ll] = el.propedg[0];
      ll++;
    }
    if (eta[i] == -1.0){
      uzly[j+i].entid[ll] = ecurve;
      uzly[j+i].prop[ll] = el.propedg[2];
      ll++;
    }
    // surface property assignment
    if(el.propsurf){
      uzly[j+i].entid[ll] = esurface;
      uzly[j+i].prop[ll] = el.propsurf[0];
      ll++;
    }
    // volume property assignment
    uzly[j+i].entid[ll] = eregion;
    uzly[j+i].prop[ll] = el.prop;
  }
  j=j+nnn;  
}

/*
Funkce najde hrany prave generovaneho elementu, ktere uz existuji a vypise je do pole edges

elements ... pole elementu
k        ... parametr udava, jaky element se testuje (treti element -> k=2)
edges    ... pole hran, v kazde hrane ulozeny 3 hodnoty: prvni nabyva hodnot 0 a 1, pricemz 0 znamena,ze se hrana neshoduje, 1 znamena, ze se hrana shoduje
                                                         druha udava s jakym elementem hrana sousedi (pouze pokud prvni hodnota je 1)
                                                         treti udava s jakou hranou element sousedi
*/
void searchsimedg(selement *elements, long k, long **edges)
{
  long i,j,a=0,b=0,c=0;
  long jedna,dve,tri,ctyri;

  for (i=0;i<4;i++)
    for (j=0;j<3;j++) edges[i][j]=0;
  jedna=0,dve=0,tri=0,ctyri=0;
  
  for (i=0;i<k;i++)
  {
    for (j=0;j<4;j++)
      {
        if (elements[k].nodes[0]==elements[i].nodes[j]) {jedna = 1;if (a!=0) b=j+1; else a=j+1;}
        if (elements[k].nodes[1]==elements[i].nodes[j]) {dve   = 1;if (a!=0) b=j+1; else a=j+1;}
        if (elements[k].nodes[2]==elements[i].nodes[j]) {tri   = 1;if (a!=0) b=j+1; else a=j+1;}
        if (elements[k].nodes[3]==elements[i].nodes[j]) {ctyri = 1;if (a!=0) b=j+1; else a=j+1;}
      }

    if ((jedna==1)&&(dve==1))   {edges[0][0]=1;edges[0][1]=i;c=0;}
    if ((dve==1)  &&(tri==1))   {edges[1][0]=1;edges[1][1]=i;c=1;}
    if ((tri==1)  &&(ctyri==1)) {edges[2][0]=1;edges[2][1]=i;c=2;}
    if ((ctyri==1)&&(jedna==1)) {edges[3][0]=1;edges[3][1]=i;c=3;}

    if ((a==1)&&(b==2)) edges[c][2]=0;
    if ((a==2)&&(b==3)) edges[c][2]=1;
    if ((a==3)&&(b==4)) edges[c][2]=2;
    if ((a==1)&&(b==4)) edges[c][2]=3;

    jedna=0;dve=0;tri=0;ctyri=0;a=0;b=0;c=0;
  }
}
/*
Za pouziti pole, kde jsou poskladany uzly v hranach (edgenod), funkce posklada uzly do jednotlivych elementu

el       ... vstupni makroelement
elements ... pole elementu
edgenod  ... dvojrozmerne pole udavajici uzly na jednotlivych hranach daneho makroelementu
edges    ... pole hran, v kazde hrane ulozeny 2 hodnoty: prvni nabyva hodnot 0 a 1, pricemz 0 znamena,ze se hrana neshoduje, 1 znamena, ze se hrana shoduje
                                                         druha udava s jakym elementem hrana sousedi (pouze pokud prvni hodnota je 1)
l        ... udava prvni jiz (nove)vygenerovany uzel v danem makroelementu
k        ... udava cislo elementu, od ktereho se bude generovat
*/
void genelements (long neh, long nev, selement &el, selement *elements, long **edgenod, long **edges, long l, long &k)
{
  long nee;               // pocet elementu
  long i,n,m,f;           // m udava posledni element ve sloupci, n znaci sloupec
  long dif=0;             // rozdil po preskoceni sloupce
  long b=0;               // indikuje spodni radek

  if ((edges[0][0]==0) && (edges[1][0]==0) && (edges[2][0]==0)) l=l+nev+2;
  if ((edges[0][0]==1) && (edges[1][0]==0) && (edges[2][0]==0)) l=l+nev+1;
  if ((edges[0][0]==0) && (edges[1][0]==1) && (edges[2][0]==0)) l=l+1;
  if ((edges[0][0]==0) && (edges[1][0]==0) && (edges[2][0]==1)) l=l+nev;
  if ((edges[0][0]==1) && (edges[1][0]==0) && (edges[2][0]==1)) l=l+nev-1;
  if ((edges[0][0]==1) && (edges[1][0]==1) && (edges[2][0]==0)) l=l+1;

  nee = neh*nev; n=1;
  if (edges[0][0]==0) dif++;
  if (edges[2][0]==0) dif++;

  for (i=0;i<nee;i++)
  {
    elements[i+k].propsurf[0]=el.propsurf[0];
    elements[i+k].prop=el.prop;
    
    m=n*nev-1;
  
    if ((n==1) && (i==0) && (b==0)) f=1;                  // levy dolni roh
    if ((n==1) && (i!=0) && (i!=m) && (b==0)) f=2;        // prvni sloupec mimo kraje
    if ((n==1) && (i==m) && (b==0)) f=3;                  // levy horni roh
    if ((n!=neh) &&(b==1)) f=4;                           // spodni radek
    if ((n!=1) && (n!=neh) && (i!=m) && (b==0)) f=5;      // uprostred
    if ((n!=1) && (n!=neh) && (i==m) && (b==0))  f=6;     // horni radek
    if ((n==neh) && (b==1)) f=7;                          // pravy dolni roh
    if ((n==neh) && (i!=m) && (b==0)) f=8;                // posledni sloupec mimo kraje
    if ((n==neh) && (i==m) && (b==0)) f=9;                // pravy horni roh
  

    switch (f)
    {
      case 1:
        elements[i+k].nodes[0]=l;
        elements[i+k].nodes[1]=edgenod[1][1];
        elements[i+k].nodes[2]=edgenod[1][0];
        elements[i+k].nodes[3]=edgenod[2][1];
        elements[i+k].propedg[1]=el.propedg[1];
        elements[i+k].propedg[2]=el.propedg[2];
        l++;
        break;
      case 2:
        elements[i+k].nodes[0]=l;
        elements[i+k].nodes[1]=edgenod[1][i+1];
        elements[i+k].nodes[2]=edgenod[1][i];
        elements[i+k].nodes[3]=l-1;
        elements[i+k].propedg[1]=el.propedg[1];
        l++;
        break;
      case 3:
        elements[i+k].nodes[0]=edgenod[0][1];
        elements[i+k].nodes[1]=edgenod[0][0];
        elements[i+k].nodes[2]=edgenod[1][nev-1];
        elements[i+k].nodes[3]=l-1;
        elements[i+k].propedg[0]=el.propedg[0];
        elements[i+k].propedg[1]=el.propedg[1];
        b=1;
        l=l+dif;
        n++;
        break;
      case 4:
        elements[i+k].nodes[0]=l;
        elements[i+k].nodes[1]=l-(nev-1+dif);
        elements[i+k].nodes[2]=edgenod[2][n-1];
        elements[i+k].nodes[3]=edgenod[2][n];
        elements[i+k].propedg[2]=el.propedg[2];
        b=0;
        l++;
        break;
      case 5:
        elements[i+k].nodes[0]=l;
        elements[i+k].nodes[1]=l-(nev-1+dif);
        elements[i+k].nodes[2]=l-(nev+dif);
        elements[i+k].nodes[3]=l-1;
        l++;
        break;
      case 6:
        elements[i+k].nodes[0]=edgenod[0][n];
        elements[i+k].nodes[1]=edgenod[0][n-1];
        elements[i+k].nodes[2]=l-(nev+dif);
        elements[i+k].nodes[3]=l-1;
        elements[i+k].propedg[0]=el.propedg[0];
        l=l+dif;
        b=1;
        n++;
        break;
      case 7:
        elements[i+k].nodes[0]=edgenod[3][1];
        elements[i+k].nodes[1]=l-(nev-1+dif);
        elements[i+k].nodes[2]=edgenod[2][neh-1];
        elements[i+k].nodes[3]=edgenod[3][0];
        elements[i+k].propedg[2]=el.propedg[2];
        elements[i+k].propedg[3]=el.propedg[3];
        b=0;
        l++;
        break;
      case 8:
        elements[i+k].nodes[0]=edgenod[3][nev-(nee-i)+1];
        elements[i+k].nodes[1]=l-(nev-1+dif);
        elements[i+k].nodes[2]=l-(nev+dif);
        elements[i+k].nodes[3]=edgenod[3][nev-(nee-i)];
        elements[i+k].propedg[3]=el.propedg[3];
        l++;
        break;
      case 9:
        elements[i+k].nodes[0]=edgenod[0][neh];
        elements[i+k].nodes[1]=edgenod[0][neh-1];
        elements[i+k].nodes[2]=l-(nev+dif);
        elements[i+k].nodes[3]=edgenod[3][nev-1];
        elements[i+k].propedg[3]=el.propedg[3];
        elements[i+k].propedg[0]=el.propedg[0];
    }
  }
  k=k+nee;
}
/*
Funkce prekopiruje uzly na hranach z makroelementu, se kterym sousedi prave generovany makroelement

edgenod  ... dvojrozmerne pole udavajici uzly na jednotlivych hranach daneho makroelementu
edges    ... viz vyse
a        ... delka hrany (neh nebo nev +1)
k        ... index do pole edgenod (cislo makroelementu)
b        ... index hrany
*/
void copy_similar_edgenod (long ***edgenod, long **edges, long a, long k, long b)
{
  long i;
  for (i=0;i<a+1;i++)
  {
    if ((b==0)||(b==2))
    {
      if (edges[b][2]==b+1) 
        edgenod[k][b][i]=edgenod[edges[b][1]][b+1][a-i];
      else 
        edgenod[k][b][i]=edgenod[edges[b][1]][edges[b][2]][i];
    }
    if ((b==1)||(b==3))
    {
      if (edges[b][2]==b-1) 
        edgenod[k][b][i]=edgenod[edges[b][1]][b-1][a-i];
      else 
        edgenod[k][b][i]=edgenod[edges[b][1]][edges[b][2]][i];
    }
  }
}
/*
Funkce nacte ze souboru deleni na jednotlivych hranach zadane topologie

in      ... vstupni soubor
edgdiv  ... pole, do ktereho se uklada deleni hran (index je cislo hrany)
a       ... pocet zadanych hran
*/
void read_edgdiv (XFILE *in, long *&edgdiv, long &a)
{
  long i;
  
  in -> kwdmode = sequent_mode;

  
  xfscanf(in,"%k %ld","pocet_zadanych_hran",&a);
  edgdiv = new long [a];
  for (i=0;i<a;i++)
    edgdiv[i]=0;

  i=0;
  while (i<a)
    xfscanf(in,"%ld %k %ld",&i,"hrana",&edgdiv[i]); 
}

/*
Funkce pro zjisteni deleni na jednotlivych hranach k-teho makroelementu

edgdiv  ... pole, kde je ulozeno deleni hran prectene z topologie
elemdiv ... dvourozmerne pole udavajici deleni na i-te hrane k-teho elementu
propedg ... pole s vlastnostmi hrana k-teho makroelementu (ze vstupni topologie)
edges   ... viz vyse
k       ... udava index makroelementu
*/

void div_elements (long *edgdiv, long **elemdiv, long *propedg, long **edges, long k)
{
  long i;
  
  for (i=0;i<4;i++)
  {
    if (propedg[i]!=0)
      elemdiv[k][i]=edgdiv[propedg[i]-1];
    else
      if (edges[i][0]==1)
        elemdiv[k][i]=elemdiv[edges[i][1]][edges[i][2]];
  }

  if (elemdiv[k][0]==0)  elemdiv[k][0]=elemdiv[k][2];
  if (elemdiv[k][1]==0)  elemdiv[k][1]=elemdiv[k][3];
  if (elemdiv[k][2]==0)  elemdiv[k][2]=elemdiv[k][0];
  if (elemdiv[k][3]==0)  elemdiv[k][3]=elemdiv[k][1];
}




