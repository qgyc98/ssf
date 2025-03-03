#include <stdio.h>
#include <stdlib.h>
#include "siftop.h"
#include "bridgen.h"
#include "quadeq.h"



/**
  Funkce vytvori elementarni prvek (dle pozadovaneho deleni neh a nev) v souradnicovem systemu ksi,eta a pote ho prevede do souradnicoveho systemu x,y podle zadanych krajnich bodu makroelementu.
  Dale vytvori pole edgenod, kde jsou ulozeny krajni uzly jednotlivych elementu.
  Vysledkem funkce jsou dve pole souradnic uzlu x,y daneho makroelementu a naplneni pole edgenod.
  

  @param neh - deleni v horizontalnim smeru (deleni pro hranu 0 a 2)
  @param nev - deleni ve vertikalnim smeru (deleni pro hranu 1 a 3)
  @param nodes - pole krajnich uzlu v makroelementu
  @param non - indexy zadanych uzlu v makroelementu
  @param x,y - pole vyslednych souradnic cele topologie
             - u pole *x a *y predpokladam, ze se do nich budou vkladat souradnice postupne ze vsech makroelementu a proto vstupuji do funkce jako parametr
  @param j - udava cislo uzlu, od ktereho se bude generovat
  @param edges - pole hran (podrobneji viz nize)
  @param edenod - trojrozmerne pole ktere uklada uzly na j-t?hrane k-tého makroelementu
  @param elemdiv - deleni makroelemntu
  @param k - parametr udava, jaky element se generuje (treti element -> k=2)
*/
void gennodes(long neh, long nev, snode *nodes, long *non, snode *uzly, long &j, long **edges,long **corner, long ***edgenod,long **elemdiv, long k)
{
  long i,nnn,m,c=0,d=0,x=0,y=0;
  double dksi,deta;
  double a=0,b=0,n;
  double *ksi,*eta;
  double e,f,g,h;       //  Interpolacni funkce pro dvourozmerny prvek, e je pro pravy horni vrchol; f,g,h jsou vrcholy postupne proti smeru hod.

  
  if (edges[0][0]==1) c++;
  if (edges[2][0]==1) c++;
  if (edges[1][0]==1) d++;
  if (edges[3][0]==1) d++; 
  
  nnn = (neh+1)*(nev+1);
  nnn = nnn - c*(neh+1) - d*(nev+1) +  c*d;     // napocita pocet uzlu v makroelementu

  for (i=0;i<4;i++)
    if (corner[i][0]==1) nnn--;

/*  if (corner[0][0]==1) {edgenod[k][0][0]=edgenod[corner[0][1]][2][0];edgenod[k][3][nev]=edgenod[corner[0][1]][2][0];}   //
    if (corner[1][0]==1) {edgenod[k][1][0]=edgenod[corner[1][1]][3][0];edgenod[k][0][neh]=edgenod[corner[1][1]][3][0];}   //   prekopiruje shodujici se rohovy uzel do edgenod
    if (corner[2][0]==1) {edgenod[k][2][0]=edgenod[corner[2][1]][0][0];edgenod[k][1][nev]=edgenod[corner[2][1]][0][0];}   //
    if (corner[3][0]==1) {edgenod[k][3][0]=edgenod[corner[3][1]][1][0];edgenod[k][2][neh]=edgenod[corner[3][1]][1][0];}   // 
*/

  if (edges[0][0]==1) copy_similar_edgenod (edgenod,edges,neh,k,0); // prekopiruji se uzly na hrane, pokud se hrana shoduje
  if (edges[1][0]==1) copy_similar_edgenod (edgenod,edges,nev,k,1); // prekopiruji se uzly na hrane, pokud se hrana shoduje
  if (edges[2][0]==1) copy_similar_edgenod (edgenod,edges,neh,k,2); // prekopiruji se uzly na hrane, pokud se hrana shoduje
  if (edges[3][0]==1) copy_similar_edgenod (edgenod,edges,nev,k,3); // prekopiruji se uzly na hrane, pokud se hrana shoduje

  if (edges[1][0]==0)           // pokud se hrana neshoduje, vygeneruje uzly na hrane (hrana vlevo)  / naplni edgenod
  {
    if (edges[2][0]==1) {for (i=1;i<(nev+1);i++) edgenod[k][1][i]=j+i-1; edgenod[k][1][0]=edgenod[k][2][0];}
      else
      {
        if (corner[2][0]==1) {for (i=1;i<(nev+1);i++) edgenod[k][1][i]=j+i-1; edgenod[k][1][0]=edgenod[(corner[2][1])][0][(elemdiv[(corner[2][1])][3])];} // <- lahudka
          else
            for (i=0;i<(nev+1);i++) edgenod[k][1][i]=j+i;
      }
    if (edges[0][0]==1)  edgenod[k][1][nev]=edgenod[k][0][0];
    if (corner[1][0]==1) edgenod[k][1][nev]=edgenod[corner[1][1]][3][0];
  }
  //done

  if (corner[0][0]==1) y=1;
  
  if (edges[3][0]==0)           // pokud se hrana neshoduje, vygeneruje uzly na hrane  (hrana vpravo) / naplni edgenod
  {
    if (edges[0][0]==1) {if (edges[2][0]==1) {for (i=1;i<(nev+1);i++) edgenod[k][3][i]=j+nnn-nev+i;   edgenod[k][3][0]=edgenod[k][2][neh];}
                           else
                            if (corner[3][0]==1) {for (i=1;i<(nev+1);i++) edgenod[k][3][i]=j+nnn-nev+i;  edgenod[k][3][0]=edgenod[corner[3][1]][0][0];}
                              else
                                for (i=0;i<(nev+1);i++) edgenod[k][3][i]=j+nnn-nev+i;}


    else                if (edges[2][0]==1) {for (i=1;i<(nev+1);i++) edgenod[k][3][i]=j+nnn-1-nev+i+y; edgenod[k][3][0]=edgenod[k][2][neh];}
                          else
                            if (corner[3][0]==1) {for (i=1;i<(nev+1);i++) edgenod[k][3][i]=j+nnn-1-nev+i+y;  edgenod[k][3][0]=edgenod[corner[3][1]][0][0];}
                              else
                                for (i=0;i<(nev+1);i++) edgenod[k][3][i]=j+nnn-1-nev+i+y;
                                
    if (edges[0][0]==1) edgenod[k][3][nev]=edgenod[k][0][neh];
    if (corner[0][0]==1) edgenod[k][3][nev]=edgenod[corner[0][1]][1][0];
  }

  // done

  ksi = new double [nnn];
  eta = new double [nnn];

  for (i=0;i<nnn;i++) ksi[i]=0.0;
  for (i=0;i<nnn;i++) eta[i]=0.0;

  dksi = (2.0)/neh;
  deta = (2.0)/nev;
  a=-1.0; b=-1.0; m=1;

  if (edges[1][0]==1) a = a + dksi;
  if ((edges[2][0]==1)||(corner[2][0]==1)) b = b + deta;
  
  for (i=0;i<nnn;i++)
  {
    ksi[i]=a;
    eta[i]=b;

    n = m*(nev+1)-1-x;   // pocitadlo posledniho prvku ve sloupci

    if (corner[2][0]==1) n--;

    if (edges[0][0]==1) n = n - m;
                   else                     // pokud se hrana neshoduje, vygeneruje uzly na hrane
                   {
                     if (edges[2][0]==1)  if (edges[1][0]==1) {edgenod[k][0][m]=j+n-m; edgenod[k][0][0]=edgenod[k][1][nev];} else edgenod[k][0][m-1]=j+n-m;
                                   else   if (edges[1][0]==1) {edgenod[k][0][m]=j+n;   edgenod[k][0][0]=edgenod[k][1][nev];} else edgenod[k][0][m-1]=j+n;
                     if (edges[3][0]==1)  edgenod[k][0][neh]=edgenod[k][3][nev];
                   }
    if (edges[2][0]==1) n = n - m;
                   else                     // pokud se hrana neshoduje, vygeneruje uzly na hrane
                   {
                     if (edges[0][0]==1)  if (edges[1][0]==1) {edgenod[k][2][m]=j+n-nev+1; edgenod[k][2][0]=edgenod[k][1][0];} else edgenod[k][2][m-1]=j+n-nev+1;
                                   else   if (edges[1][0]==1) {edgenod[k][2][m]=j+n-nev;   edgenod[k][2][0]=edgenod[k][1][0];} else edgenod[k][2][m-1]=j+n-nev;
                     if (edges[3][0]==1)  edgenod[k][2][neh]=edgenod[k][3][0];
                   }
                   
    if ((corner[3][0]==1)&&(i==(nnn-nev))) {b+=deta;eta[i]=b;}
    if ((corner[3][0]==1)&&(edges[0][0]==1)&&(i==(nnn-nev+1))) {b+=deta;eta[i]=b;}
    
    if (i==n) {if (edges[2][0]==1) b=-1.0+deta; else b=-1.0; a=a+dksi; m++;} else b=b+deta;
    if ((corner[1][0]==1)&&(m==1)&&(i==(n-1))) {if (edges[2][0]==1) b=-1.0+deta; else b=-1.0; a=a+dksi; x=1; m++;}

  }
  if (corner[1][0]==1) edgenod[k][0][0]=edgenod[k][1][nev];
  if (corner[2][0]==1) edgenod[k][2][0]=edgenod[k][1][0];
  if (corner[3][0]==1) {edgenod[k][2][neh]=edgenod[k][3][0]; edgenod[k][0][neh]--;}
  if (corner[0][0]==1) edgenod[k][0][neh]=edgenod[k][3][nev];
  
  
  for (i=0;i<nnn;i++)
  {
    e = 0.25*(1+ksi[i])*(1+eta[i]);
    f = 0.25*(1-ksi[i])*(1+eta[i]);
    g = 0.25*(1-ksi[i])*(1-eta[i]);
    h = 0.25*(1+ksi[i])*(1-eta[i]);

    uzly[j+i].x = e*nodes[non[0]].x + f*nodes[non[1]].x + g*nodes[non[2]].x + h*nodes[non[3]].x;
    uzly[j+i].y = e*nodes[non[0]].y + f*nodes[non[1]].y + g*nodes[non[2]].y + h*nodes[non[3]].y;
    uzly[j+i].z = nodes[0].z;
  }
  j=j+nnn;
  
  delete [] ksi;
  delete [] eta;
}



/**
  Funkce zjisti, zda prave generovany makroelement sousedi s jinym jiz EXISTUJICI makroelementem, a tuto informaci zapise do pole edges zpusobem napsanym nize.
  Pro tuto funkci je obzvlast dulezite poradi generovani makroelementu. Pokud dojde k tomu, ze se zacne generovat makroelement, ktery sousedi s jinym pouze krajnim uzlem
  a ne hranou, funkce si toho nevsimne a dojde ke zdvojeni uzlu. Vysledkem funkce je tedy naplneni pole edges.


  @param elements - pole elementu
  @param k - parametr udava, jaky element se testuje (treti element -> k=2)
  @param edges - pole hran, v kazde hrane ulozeny 3 hodnoty: - prvni nabyva hodnot 0 a 1, pricemz 0 znamena,ze se hrana neshoduje, 1 znamena, ze se hrana shoduje
                                                             - druha udava s jakym elementem hrana sousedi (pouze pokud prvni hodnota je 1)
                                                             - treti udava s jakou hranou element sousedi
  @param corner - udava, pokud se prekryva makroel pouze rohem, index rika, o ktery roh se jedna, a hodnota 1 je pravda, 0 nepravda
*/
void searchsimedg(selement *&elements, long k, long **edges, long **corner)
{
  long i,j,l,a=0,b=0,c=0;
  long jedna,dve,tri,ctyri;

  for (i=0;i<4;i++)
  {
    for (j=0;j<3;j++) edges[i][j]=0;
    for (l=0;l<3;l++) corner[i][l]=0;
  }
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

   
    if ((jedna==1)&&(dve==0)&&(ctyri==0)) {corner[0][0]=1;corner[0][1]=i;}
    if ((jedna==0)&&(dve==1)&&(tri==0))   {corner[1][0]=1;corner[1][1]=i;}
    if ((dve==0)&&(tri==1)&&(ctyri==0))   {corner[2][0]=1;corner[2][1]=i;}
    if ((jedna==0)&&(tri==0)&&(ctyri==1)) {corner[3][0]=1;corner[3][1]=i;}
    
    if ((a==1)&&(b==2)) edges[c][2]=0;
    if ((a==2)&&(b==3)) edges[c][2]=1;
    if ((a==3)&&(b==4)) edges[c][2]=2;
    if ((a==1)&&(b==4)) edges[c][2]=3;

    jedna=0;dve=0;tri=0;ctyri=0;a=0;b=0;c=0;
  }

  // nasledujici serie prikazu opravuje chybu v nalezani "corneru", ktera muze nastat urcitym usporadanim zadani
  
  if (edges[0][0]==1) {corner[0][0]=0; corner[1][0]=0;}
  if (edges[1][0]==1) {corner[1][0]=0; corner[2][0]=0;}
  if (edges[2][0]==1) {corner[2][0]=0; corner[3][0]=0;}
  if (edges[3][0]==1) {corner[3][0]=0; corner[0][0]=0;}
}



/**
  Za pouziti pole, kde jsou poskladany uzly v hranach (edgenod), funkce zatridi uzly do jednotlivych elementu v danem makroelementu
  Vysledkem funkce je naplneni pole nodes, propedg a propsurf pro kazdy element v generovanem makroelementu.


  @param propedg - pole s vlastnosti hrana vstpujiciho makroelementu
  @param propsurf - pole s vlastnosti plocha vstupujicho makroelementu
  @param prope    - vlastnost vstupujicho makroelementu
  @param elements - pole elementu
  @param edgenod - dvojrozmerne pole udavajici uzly na jednotlivych hranach daneho makroelementu
  @param edges - pole hran, podrobneji viz vyse
  @param l - udava prvni jiz (nove)vygenerovany uzel v danem makroelementu
  @param k - udava cislo elementu, od ktereho se bude generovat
*/
void genelements (long neh, long nev, long prope, long *propedg, long *propsurf, selement *elements, long **edgenod, long **edges, long **corner, long l, long &k)
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
  if (corner[1][0]==1) l--;
  if (corner[2][0]==1) l--;

  nee = neh*nev; n=1;
  if (edges[0][0]==0) dif++;
  if (edges[2][0]==0) dif++;

  for (i=0;i<nee;i++)
  {
    elements[i+k].propsurf[0]=propsurf[0];
    elements[i+k].prop=prope;
    
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
        elements[i+k].propedg[1]=propedg[1];
        elements[i+k].propedg[2]=propedg[2];
        l++;
        break;
      case 2:
        elements[i+k].nodes[0]=l;
        elements[i+k].nodes[1]=edgenod[1][i+1];
        elements[i+k].nodes[2]=edgenod[1][i];
        elements[i+k].nodes[3]=l-1;
        elements[i+k].propedg[1]=propedg[1];
        l++;
        break;
      case 3:
        elements[i+k].nodes[0]=edgenod[0][1];
        elements[i+k].nodes[1]=edgenod[0][0];
        elements[i+k].nodes[2]=edgenod[1][nev-1];
        elements[i+k].nodes[3]=l-1;
        elements[i+k].propedg[0]=propedg[0];
        elements[i+k].propedg[1]=propedg[1];
        b=1;
        l=l+dif;
        n++;
        break;
      case 4:
        elements[i+k].nodes[0]=l;
        elements[i+k].nodes[1]=l-(nev-1+dif);
        elements[i+k].nodes[2]=edgenod[2][n-1];
        elements[i+k].nodes[3]=edgenod[2][n];
        elements[i+k].propedg[2]=propedg[2];
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
        elements[i+k].propedg[0]=propedg[0];
        l=l+dif;
        b=1;
        n++;
        break;
      case 7:
        elements[i+k].nodes[0]=edgenod[3][1];
        elements[i+k].nodes[1]=l-(nev-1+dif);
        elements[i+k].nodes[2]=edgenod[2][neh-1];
        elements[i+k].nodes[3]=edgenod[3][0];
        elements[i+k].propedg[2]=propedg[2];
        elements[i+k].propedg[3]=propedg[3];
        b=0;
        l++;
        break;
      case 8:
        elements[i+k].nodes[0]=edgenod[3][nev-(nee-i)+1];
        elements[i+k].nodes[1]=l-(nev-1+dif);
        elements[i+k].nodes[2]=l-(nev+dif);
        elements[i+k].nodes[3]=edgenod[3][nev-(nee-i)];
        elements[i+k].propedg[3]=propedg[3];
        l++;
        break;
      case 9:
        elements[i+k].nodes[0]=edgenod[0][neh];
        elements[i+k].nodes[1]=edgenod[0][neh-1];
        elements[i+k].nodes[2]=l-(nev+dif);
        elements[i+k].nodes[3]=edgenod[3][nev-1];
        elements[i+k].propedg[3]=propedg[3];
        elements[i+k].propedg[0]=propedg[0];
    }
  }
  k=k+nee;
}



/**
  Funkce prekopiruje uzly na hranach z makroelementu, se kterym sousedi prave generovany makroelement.
  Pomocna funkce, pro castecne zprehledneni kodu funkce gennodes.
  
  @param edgenod - dvojrozmerne pole udavajici uzly na jednotlivych hranach daneho makroelementu
  @param edges - pole hran, podrobneji viz vyse
  @param a - delka hrany (neh nebo nev +1)
  @param k - index do pole edgenod (cislo makroelementu)
  @param b - index hrany
*/
void copy_similar_edgenod (long ***edgenod, long **edges, long a, long k, long b)
{
  long i;
  for (i=0;i<a+1;i++)
  {
    if ((b==0)||(b==2))
    {
      if (edges[b][2]==b+1) edgenod[k][b][i]=edgenod[edges[b][1]][b+1][a-i];
      else edgenod[k][b][i]=edgenod[edges[b][1]][edges[b][2]][i];
    }
    if ((b==1)||(b==3))
    {
      if (edges[b][2]==b-1) edgenod[k][b][i]=edgenod[edges[b][1]][b-1][a-i];
      else edgenod[k][b][i]=edgenod[edges[b][1]][edges[b][2]][i];
    }
  }
}



/**
  Funkce nacte vstupni soubor a alokuje potrebne promenne

  @param in - vstupni soubor
  @param edgdiv - pole, do ktereho se uklada deleni hran (index je cislo hrany)
  @param arsec - pole useku
  @param cuts - pole zadanych char. rezu
  @param numcuts - pocet rezu
  @param topcuts - pole vygenerovanych char.rezu
*/
void read (XFILE *in, long *&edgdiv, section *&arsec, siftop *&cuts, long &numcuts, siftop *&topcuts, long &numsec)
{
  long a;           // pocet zadanych hran
  long i, j;
  char filename[20];
  XFILE *fn;

  printf("\n\n");
  printf("Cteni vstupu:\n");
  
  in -> kwdmode = sequent_mode;

  xfscanf(in,"%k %ld","pocet_charakteristickych_rezu",&numcuts);
  printf(" - Pocet charakteristickych rezu: %ld \n",numcuts);
  cuts  = new siftop [numcuts];
  topcuts = new siftop [numcuts]; for (i=0;i<numcuts;i++) {topcuts[i].ne=0;topcuts[i].nn=0;}
  i=0;
  printf(" - Cteni charakteristickych rezu... ");
  while (i<numcuts)
  {
    xfscanf(in,"%ld%19s",&i,filename);
    fn = xfopen(filename, "rt");
    cuts[i-1].read(fn,0,1);
    xfclose(fn);
  }
  printf("OK \n");
  
  xfscanf(in,"%k %ld","pocet_zadanych_hran",&a);
  printf(" - Pocet zadanych hran: %ld \n",a);
  edgdiv = new long [a];
  for (i=0;i<a;i++)
    edgdiv[i]=0;
  i=0;
  while (i<a)
    xfscanf(in,"%ld %k %ld",&i,"hrana",&edgdiv[i]);

  xfscanf(in,"%k %ld","pocet_zadanych_useku",&numsec);
  printf(" - Pocet zadanych useku: %ld \n",numsec);
  arsec=new section[numsec];
  i=0;
  printf(" - Cteni useku... ");
  for(j=0; j<numsec; j++)
  {
    xfscanf(in,"%ld", &i); 
    i--;
    xfscanf(in,"%k %ld %ld %k %le %k %ld %k %ld %+k %ld",
            "usek",&arsec[i].a,&arsec[i].b,"delka",&arsec[i].length,
            "deleni",&arsec[i].division,"prubeh",&arsec[i].quadlin, "propid", &arsec[i].prop);
    arsec[i].a--;
    arsec[i].b--;
/*    arsec[i-1].a+=-1;
    arsec[i-1].b+=-1;*/
  }
  printf("OK \n");
  
}



/**
  Funkce, ktera urci, jake deleni hran ma mit dany makroelement, podle sousedicich makroelemntu.
  Vysledkem je naplneni pole elemdiv.

  @param edgdiv - pole, kde je ulozeno deleni hran prectene z topologie
  @param elemdiv - dvourozmerne pole udavajici deleni na i-te hrane k-teho elementu
  @param propedg - pole s vlastnostmi hrana k-teho makroelementu (ze vstupni topologie)
  @param edges - pole hran, podrobneji viz vyse
  @param k - udava index makroelementu
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



/**
  Funkce vygeneruje 3d topologii (uzly i elementy) v jednom useku mezi dvema rezy.


  @param sec - prave generovanz usek, ve kterem jsou informace o krajnich rezech a delce
  @param topcuts - pole vygenerovanych char.rezu
  @param top - topologie, do ktere se uklada
  @param numel - celkove pocitadlo na elementy
  @param numnod - celkove pocitadlo na uzly
  @param zetko - udava zetovou souradnici, od ktere se generuje dany usek; po dokonceni generace se zvysi o delku useku,
                 udavajici zacatek dalsiho useku v poradi
               - zetova souradnice uvazuju ve smeru strednice mostu
                 (v mainu, pred vystupem do souboru, jsou osy upraveny: x ve smeru strednice, y a z v rovine prurezu)
*/
void gensection (section &sec, siftop *&topcuts, siftop &top, long &numel, long &numnod, double &zetko)
{
  long i,j,w;
  double krokx,kroky,krokz;
  long pocet=(sec.division+1);   // pocet rezu v useku (vcetne krajnich)
  long cn,ce,cp=0,cpa=0;         // parametry, ktere osetruji rozdil poctu elementu v jednotlivych rezech v danem useku
  double a=0.0,b=0.0,c=0.0;
  double pom;

  if (topcuts[sec.a].ne > topcuts[sec.b].ne)
  {
    cn = topcuts[sec.b].nn;
    ce = topcuts[sec.b].ne;
    cp =(topcuts[sec.a].nn)-(topcuts[sec.b].nn);
  }
  else
  {
    cn = topcuts[sec.a].nn;
    ce = topcuts[sec.a].ne;
  }
// cyklus vygeneruje vsechny souradnice v useku
  for (i=0;i<cn;i++)
  {
    if (sec.quadlin == 2)
      getparam(a,b,c,1.0,topcuts[sec.a].nodes[i].y,sec.length+1.0,topcuts[sec.b].nodes[i].y);
    krokx=((topcuts[sec.b].nodes[i].x)-(topcuts[sec.a].nodes[i].x))/(pocet-1);
    krokz=(sec.length)/(pocet-1);
    for (j=0;j<pocet;j++)
    {
      pom = j*krokz +1;
      if (sec.quadlin == 2)         // kvadraticky prubeh
      {
        kroky = a*pom*pom + b*pom + c ;
        top.nodes[numnod+i+cpa+j*cn].y=kroky;
      }
      if (sec.quadlin == 1)         // linearni prubeh
      {
        kroky=((topcuts[sec.b].nodes[i].y)-(topcuts[sec.a].nodes[i].y))/(pocet-1);
        kroky*=j;
        top.nodes[numnod+i+cpa+j*cn].y=topcuts[sec.a].nodes[i].y+kroky;
      }
      top.nodes[numnod+i+cpa+j*cn].x=topcuts[sec.a].nodes[i].x+j*krokx;
      top.nodes[numnod+i+cpa+j*cn].z=topcuts[sec.a].nodes[i].z+zetko+j*krokz;
      cpa=cp;      
    }
    cpa=0;
  }
  zetko+=sec.length;
  
// cyklus vygeneruje elementy
  for (w=0;w<(pocet-1);w++)
  {
    for (i=0;i<ce;i++)
    {
/*    // Puvodni cislovani
      top.elements[numel+i+w*ce].nodes[0]=numnod+cpa+topcuts[sec.a].elements[i].nodes[0]+w*cn;
      top.elements[numel+i+w*ce].nodes[1]=numnod+cpa+topcuts[sec.a].elements[i].nodes[1]+w*cn;
      top.elements[numel+i+w*ce].nodes[2]=numnod+cpa+topcuts[sec.a].elements[i].nodes[2]+w*cn;
      top.elements[numel+i+w*ce].nodes[3]=numnod+cpa+topcuts[sec.a].elements[i].nodes[3]+w*cn;

      top.elements[numel+i+w*ce].nodes[4]=numnod+cp+topcuts[sec.a].elements[i].nodes[0]+(w+1)*cn;
      top.elements[numel+i+w*ce].nodes[5]=numnod+cp+topcuts[sec.a].elements[i].nodes[1]+(w+1)*cn;
      top.elements[numel+i+w*ce].nodes[6]=numnod+cp+topcuts[sec.a].elements[i].nodes[2]+(w+1)*cn;
      top.elements[numel+i+w*ce].nodes[7]=numnod+cp+topcuts[sec.a].elements[i].nodes[3]+(w+1)*cn;
*/
      top.elements[numel+i+w*ce].nodes[0]=numnod+cpa+topcuts[sec.a].elements[i].nodes[0]+(w+1)*cn;
      top.elements[numel+i+w*ce].nodes[1]=numnod+cpa+topcuts[sec.a].elements[i].nodes[1]+(w+1)*cn;
      top.elements[numel+i+w*ce].nodes[2]=numnod+cpa+topcuts[sec.a].elements[i].nodes[2]+(w+1)*cn;
      top.elements[numel+i+w*ce].nodes[3]=numnod+cpa+topcuts[sec.a].elements[i].nodes[3]+(w+1)*cn;

      top.elements[numel+i+w*ce].nodes[4]=numnod+cp+topcuts[sec.a].elements[i].nodes[0]+w*cn;
      top.elements[numel+i+w*ce].nodes[5]=numnod+cp+topcuts[sec.a].elements[i].nodes[1]+w*cn;
      top.elements[numel+i+w*ce].nodes[6]=numnod+cp+topcuts[sec.a].elements[i].nodes[2]+w*cn;
      top.elements[numel+i+w*ce].nodes[7]=numnod+cp+topcuts[sec.a].elements[i].nodes[3]+w*cn;

      if (sec.prop < 0) 
      {
        // pokud se property nezadalo stejne v ramci celeho useku (tj. sec.prop<0), 
        // tak se bere jako objemova property z 1. char rezu (sec.a) daneho useku
        top.elements[numel+i+w*ce].prop=topcuts[sec.a].elements[i].prop;
      }
      else 
      {
        // objemova property zadana v ramci celeho useku stejna pro vsechny 3D prvky
        top.elements[numel+i+w*ce].prop=sec.prop;
      }
      top.elements[numel+i+w*ce].propsurf[1]=topcuts[sec.a].elements[i].propedg[0];
      top.elements[numel+i+w*ce].propsurf[2]=topcuts[sec.a].elements[i].propedg[1];
      top.elements[numel+i+w*ce].propsurf[3]=topcuts[sec.a].elements[i].propedg[2];
      top.elements[numel+i+w*ce].propsurf[0]=topcuts[sec.a].elements[i].propedg[3];
    }
    cpa=cp;
  }
  numnod+=(pocet-1)*cn+cp;
  numel +=(pocet-1)*ce;
}


/**
  Funkce vygeneruje 3d topologii (uzly i elementy) v jednom useku mezi dvema rezy.


  @param pcut_id - index prurezu z predchoziho useku nebo -1 pokud jde o uplne prvni usek
  @param sec - prave generovany usek, ve kterem jsou informace o krajnich rezech a delce
  @param topcuts - pole vygenerovanych char.rezu
  @param top - topologie, do ktere se uklada
  @param numel - celkove pocitadlo na elementy
  @param numnod - celkove pocitadlo na uzly
  @param zetko - udava zetovou souradnici, od ktere se generuje dany usek; po dokonceni generace se zvysi o delku useku,
                 udavajici zacatek dalsiho useku v poradi
               - zetova souradnice uvazuju ve smeru strednice mostu
                 (v mainu, pred vystupem do souboru, jsou osy upraveny: x ve smeru strednice, y a z v rovine prurezu)
*/
void gensection2 (long pcut_id, section &sec, siftop *&topcuts, siftop &top, long &numel, long &numnod, double &zetko)
{
  long i,j,w;
  double dx, dy, dz;
  long pocet = sec.division;   // pocet rezu v useku, ktere se budou generovat
  long cn, ce, pcn;
  double a=0.0, b=0.0, c=0.0;
  double pom;

  cn = topcuts[sec.a].nn;
  ce = topcuts[sec.a].ne;
  // cyklus vygeneruje vsechny souradnice v useku
  dz=(sec.length)/(pocet);
  for (i=0; i<cn; i++)
  {
    if (sec.quadlin == 2)
      getparam(a,b,c,1.0,topcuts[sec.a].nodes[i].y,sec.length+1.0,topcuts[sec.b].nodes[i].y);
    dx=((topcuts[sec.b].nodes[i].x)-(topcuts[sec.a].nodes[i].x))/(pocet);
    double az = zetko+dz;
    for (j=0; j<pocet; j++)
    {
      pom = j*dz +1;
      if (sec.quadlin == 2)         // kvadraticky prubeh
      {
        dy = a*pom*pom + b*pom + c ;
        top.nodes[numnod+i+j*cn].y=dy;
      }
      if (sec.quadlin == 1)         // linearni prubeh
      {
        dy=((topcuts[sec.b].nodes[i].y)-(topcuts[sec.a].nodes[i].y))/(pocet);
        dy*=j;
        top.nodes[numnod+i+j*cn].y=topcuts[sec.a].nodes[i].y+dy;
      }
      top.nodes[numnod+i+j*cn].x=topcuts[sec.a].nodes[i].x+j*dx;
      top.nodes[numnod+i+j*cn].z=topcuts[sec.a].nodes[i].z+az;
      az += dz;
    }    
  }
  zetko+=sec.length;
  
  // cyklus vygeneruje elementy
  if (pcut_id >= 0){
    pcn = topcuts[pcut_id].nn;
    if (pcn < cn)   pcn = cn;      
  }
  else
    pcn = cn;
  
  for (w=0; w<pocet; w++)
  {
    for (i=0;i<ce;i++)
    {
      long eid = numel+i+w*ce;
      top.elements[eid].nodes[0] = numnod + topcuts[sec.a].elements[i].nodes[0] + w*cn;
      top.elements[eid].nodes[1] = numnod + topcuts[sec.a].elements[i].nodes[1] + w*cn;
      top.elements[eid].nodes[2] = numnod + topcuts[sec.a].elements[i].nodes[2] + w*cn;
      top.elements[eid].nodes[3] = numnod + topcuts[sec.a].elements[i].nodes[3] + w*cn;

      top.elements[eid].nodes[4] = top.elements[eid].nodes[0] - pcn;
      top.elements[eid].nodes[5] = top.elements[eid].nodes[1] - pcn;
      top.elements[eid].nodes[6] = top.elements[eid].nodes[2] - pcn;
      top.elements[eid].nodes[7] = top.elements[eid].nodes[3] - pcn;

      if (sec.prop < 0) 
      {
        // pokud se property nezadalo stejne v ramci celeho useku (tj. sec.prop<0), 
        // tak se bere jako objemova property z 1. char rezu (sec.a) daneho useku
        top.elements[numel+i+w*ce].prop=topcuts[sec.a].elements[i].prop;
      }
      else 
      {
        // objemova property zadana v ramci celeho useku stejna pro vsechny 3D prvky
        top.elements[numel+i+w*ce].prop=sec.prop;
      }
      top.elements[eid].propsurf[1]=topcuts[sec.a].elements[i].propedg[0];
      top.elements[eid].propsurf[2]=topcuts[sec.a].elements[i].propedg[1];
      top.elements[eid].propsurf[3]=topcuts[sec.a].elements[i].propedg[2];
      top.elements[eid].propsurf[0]=topcuts[sec.a].elements[i].propedg[3];
    }
    pcn = cn;
  }
  numnod += pocet*cn;
  numel  += pocet*ce;
}


/**
  Pomocna funkce, ktera alokuje potrebne promenne.

  @param edges - pole hran
  @param elemdiv - dvourozmerne pole udavajici deleni na i-te hrane k-teho elementu
  @param numcuts - pocet rezu v useku
  @param cuts - pole zadanych charakteristickych rezu  
*/
void alloc (long **&edges, long ***&elemdiv, long numcuts, siftop *&cuts, long **&corner)
{
  long i,w,z;
  
  edges = new long*[4];
  for (i=0;i<4;i++)
  {
    edges[i]= new long[3];
    for (z=0;z<3;z++) edges [i][z]=0;
  }

  elemdiv = new long**[numcuts];
  for (w=0;w<numcuts;w++)
  {
    elemdiv[w] = new long*[cuts[w].ne];
    for (i=0;i<cuts[w].ne;i++)
    {
      elemdiv[w][i]= new long[4];
      for (z=0;z<4;z++) elemdiv[w][i][z]=0;
    }
  }

  corner = new long*[4];
  for (i=0;i<4;i++)
  {
    corner[i]= new long[3];
    for (z=0;z<3;z++) corner [i][z]=0;
  }
}


/**
  Pomocna funkce, ktera alokuje potrebne promenne.

  @param edges - pole hran
  @param elemdiv - dvourozmerne pole udavajici deleni na i-te hrane k-teho elementu
  @param numcuts - pocet rezu v useku
  @param cuts - pole zadanych charakteristickych rezu  
*/
void dealloc (long **edges, long ***elemdiv, long numcuts, siftop *cuts, long **corner)
{
  long i,w;
  
  for (i=0;i<4;i++)
    delete [] edges[i];
  delete [] edges;

  for (w=0;w<numcuts;w++)
  {
    for (i=0;i<cuts[w].ne;i++)
      delete [] elemdiv[w][i];
    delete [] elemdiv[w];
  }
  delete [] elemdiv;

  for (i=0;i<4;i++)
    delete [] corner[i];
  delete [] corner;
}


/**
  The function copies a set of nodes from the section topology sect to the bridge topology top. 
  The set of nodes to be copied is given by the id of the first node and the number of nodes  
  num_nodes to be copied. The nodes are copied to the array of nodes in the topology top starting 
  from the anid position. Previously allocated nodes in the topology top will be reallocated.

  @param sect[in] - topology of the section whose nodes are to be copied,
  @param ini_nodsec_id[in] - id of the initial node in the topology sect from which the copying 
                             will begin,
  @param num_nodes[in] - number of nodes to be copied from the section sect,
  @param anid[in/out] - current index in the array of nodes of the topology top from which 
                        to start placing copies of nodes,
  @param top[in/out] - the resulting topology where the nodes will be copied to,

  @return The function does not returns anything but the actual value of node index in the  
          argumnet anid. Copied nodes are stored in the array top.nodes.

  Created by Tomas Koudelka, 04.2023
*/
void copy_section_nodes2top(siftop &sect, long ini_nodesc_id, long num_nodes, long &anid, siftop &top)
{
  long id = ini_nodesc_id;
  for(long i=0; i<num_nodes; i++, id++, anid++){
    sect.nodes[id].copyto(top.nodes[anid]);
  }
}



/**
  The function copies a set of nodes from the section topology sect to the bridge topology top. 
  The set of nodes to be copied is given by the id of the first node and the number of nodes  
  num_nodes to be copied. The nodes are copied to the array of nodes in the topology top starting 
  from the anid position. Previously allocated nodes in the topology top will be reallocated. 
  Coordinates of the node copies are shifted by the prescribed coordinate increments dx, dy, dz.

  @param sect[in] - topology of the section whose nodes are to be copied,
  @param ini_nodsec_id[in] - id of the initial node in the topology sect from which the copying 
                             will begin,
  @param num_nodes[in] - number of nodes to be copied from the section sect,
  @param anid[in/out] - current index in the array of nodes of the topology top from which to 
                        start placing copies of nodes,
  @param dx[in] - shift of the x cooridnate of the copied nodes
  @param dy[in] - shift of the y cooridnate of the copied nodes
  @param dz[in] - shift of the z cooridnate of the copied nodes
  @param top[in/out] - the resulting topology where the nodes will be copied to,

  @return The function does not returns anything but the actual value of node index in the 
          argumnet anid. Copied nodes are stored in the array top.nodes.

  Created by Tomas Koudelka, 04.2023
*/
void copy_shift_set_section_nodes2top(siftop &sect, long ini_nodesc_id, long num_nodes, double dx,
                                      double dy, double dz, long &anid, siftop &top)
{
  long id = ini_nodesc_id;
  for(long i=0; i<num_nodes; i++, id++, anid++){
    sect.nodes[id].copyto(top.nodes[anid]);
    top.nodes[anid].x += dx;
    top.nodes[anid].y += dy;
    top.nodes[anid].z += dz;
  }
}
