#include "cebfipcontactmat.h"
#include "iotools.h"
#include "matrix.h"
#include "vector.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "intpoints.h"
#include <math.h>

cebfipcontactmat::cebfipcontactmat (void)
{
  fcc = 0.0; fct = 0.0; conf = 0; bond = 0; normal = 0.0;

  s = 0.0; ss = 0.0; sss = 0.0; taumax = 0.0; tauf = 0.0; alfa = 0.0;
}

cebfipcontactmat::~cebfipcontactmat (void)
{
  
}

void cebfipcontactmat::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %ld %ld %lf",&fcc,&fct,&conf,&bond,&normal);

  if (conf==1)  // neovinuty beton
  {
    s=ss=0.6e-3; alfa=0.4;
    if (bond==1) {sss=1.0e-3;taumax=2*sqrt(fcc/1e6)*1e6;tauf=0.3*sqrt(fcc/1e6)*1e6;}  // dobre podminky soudrznosti
    if (bond==2) {sss=2.5e-3;taumax=sqrt(fcc/1e6)*1e6;tauf=0.15*sqrt(fcc/1e6)*1e6;}   // spatne podminky soudrznosti
  }
  if (conf==2) // ovinuty beton
  {
    s=1.0e-3; ss=3.0e-3; alfa=0.4;
    if (bond==1) {taumax=2.5*sqrt(fcc/1e6)*1e6;tauf=1.0*sqrt(fcc/1e6)*1e6;}  // dobre podminky soudrznosti
    if (bond==2) {taumax=1.25*sqrt(fcc/1e6)*1e6;tauf=0.5*sqrt(fcc/1e6)*1e6;}   // spatne podminky soudrznosti
  }
}

void cebfipcontactmat::matstiff (matrix &d,long ipp)
{
  long i, n = Mm->ip[ipp].ncompstr;
  vector eps(n);

  for (i=0;i<n;i++)
    eps[i]=Mm->ip[ipp].strain[i];

  eps[0]=fabs(eps[0]);
  if (eps[0]<0.05*s)
    eps[0] = 0.05*s;
  if (eps[0]<0.8*s)  
    d[0][0]=taumax/pow(s,alfa)*alfa*pow(eps[0],(alfa-1));
  else
    d[0][0]=taumax/pow(s,alfa)*alfa*pow(0.8*s,(alfa-1));
    
  d[0][0] = d[0][0]*0.0314;        // prenasobeni obvodem prutu
  d[1][1] = normal;
}

void cebfipcontactmat::nlstresses (long ipp, long /*im*/, long /*ido*/)
{
  long i, n = Mm->ip[ipp].ncompstr;
  vector eps(n),sig(n),dpp(n),dpmax(n);

  //  initial values
  for (i=0;i<n;i++)
  {
    eps[i]=Mm->ip[ipp].strain[i];
    dpmax[i]=Mm->ip[ipp].eqother[i];
  }

  if (fabs(eps[0])<s)
    sig[0]=taumax*pow(fabs(eps[0])/s,alfa)*sgn(eps[0]);

  if ((fabs(eps[0])>s) && (fabs(eps[0])<ss))
    sig[0]=taumax*sgn(eps[0]);

  if ((fabs(eps[0])>ss) && (fabs(eps[0])<sss))
    sig[0]=taumax-(taumax-tauf)/(sss-ss)*(eps[0]-ss);

  if (fabs(eps[0])>sss)
    sig[0]=tauf*sgn(eps[0]);


  sig[1]=eps[1]*normal;
  if (sig[1]>fct)
    sig[1]=fct;
  if (sig[1]<-fcc)
    sig[1]=-fcc;

  for (i=0;i<n;i++)
    Mm->ip[ipp].stress[i]=sig[i];
  
}

void cebfipcontactmat::updateval (long ipp, long im, long ido)
{
  long i,n = Mm->givencompeqother(ipp, im);

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }
}

long cebfipcontactmat::sgn (double a)
{
  long b;
  if (a==0) b=0;
  if (a<0)  b=-1;
  if (a>0)  b=1;
  return (b);
}

