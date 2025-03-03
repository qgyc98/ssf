#include "winpast.h"
#include "stochdriver.h"
#include "global.h"
#include "mechmat.h"
#include "intpoints.h"

winpast::winpast (void)
{
  allocv (3,c1);  allocv (3,c2);
}
winpast::~winpast (void)
{
}

void winpast::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf",&c1[0],&c1[1],&c1[2]);
  xfscanf (in,"%lf %lf %lf",&c2[0],&c2[1],&c2[2]);
}



/**
  The function prints the material parameters to the opened text file given 
  by argument out.

  @param out - pointer to the opened text file

  @return The function does not return anything.

  Created by Tomas Koudelka, 10.2015
*/
void winpast::print (FILE *out)
{
  fprintf (out,"%le %le %le",c1[0],c1[1],c1[2]);
  fprintf (out," %le %le %le",c2[0],c2[1],c2[2]);
}

void winpast::matstiff (matrix &d,long ipp)
{
  switch (Mm->ip[ipp].ssst){
  case plbeam:{
    matstiff_soilplbeam (d);
    break;
  }
  case spacebeam:{
    matstiff_soilbeam (d);
    break;
  }
  case platek:{
    matstiff_soilplate (d);
    break;
  }
  case plates:{
    matstiff_soilplate (d);
    break;
  }
  case spacestress:{
//    matstiff_spacestr (d);
//    break;
  }
  default:{
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function winpast::matstiff (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}

void winpast::elmatstiff (matrix &d,long ipp)
{
  switch (Mm->ip[ipp].ssst){
  case plbeam:{
    matstiff_soilplbeam (d);
    break;
  }
  case spacebeam:{
    matstiff_soilbeam (d);
    break;
  }
  case platek:{
    matstiff_soilplate (d);
    break;
  }
  case plates:{
    matstiff_soilplate (d);
    break;
  }
  case spacestress:{
//    matstiff_spacestr (d);
//    break;
  }
  default:{
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function winpast::elmatstiff (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}

/**
   function creates stiffness matrix of the elastic
   isotropic material for plane beam elements
   
   @param d - stiffness matrix of the material
   
   11.9.2001
*/
void winpast::matstiff_soilplbeam (matrix &d)
{

  d[0][0] = c1[0];
  d[1][1] = c1[1];
  d[2][2] = c2[0];

}
void winpast::matstiff_soilbeam (matrix &d)
{

  d[0][0] = c1[0];
  d[1][1] = c1[1];
  d[2][2] = c1[2];
  d[3][3] = c2[0];
  d[4][4] = c2[1];
  d[5][5] = c2[2];

}


/**
   function creates stiffness matrix of the elastic
   isotropic material for plate elements
   
   @param d - stiffness matrix
   
   19.7.2001
*/
void winpast::matstiff_soilplate (matrix &d)
{

    
  d[0][0]=c1[0];    d[0][1]=0.0;      d[0][2]=0.0;
  d[1][0]=d[0][1];  d[1][1]=c2[0];    d[1][2]=0.0;
  d[2][0]=0.0;      d[2][1]=0.0;      d[2][2]=c2[0];
  
}

/**
   function creates stiffness matrix of the elastic
   isotropic material for 3D problems
   
   @param d - stiffness matrix of the material

   19.7.2001
*/
void winpast::matstiff_soilspacestr (matrix &/*d*/)
{
  /*
  double g,s;
  
  fillm(0.0,d);
  
  g = e/2.0/(1.0+nu);
  s = e/(1.0+nu)/(1.0-2.0*nu);
  
  d[0][0]=s*(1-nu);  d[0][1]=s*nu;     d[0][2]=s*nu;
  d[1][0]=d[0][1];   d[1][1]=d[0][0];  d[1][2]=d[0][1];
  d[2][0]=d[0][1];   d[2][1]=d[0][1];  d[2][2]=d[0][0];

  d[3][3]=g;         d[4][4]=g;        d[5][5]=g;
*/
}


void winpast::nlstresses (long ipp)
{
  long i, n = Mm->ip[ipp].ncompstr;
  vector eps(n),sig(n);
  strastrestate ssst = Mm->ip[ipp].ssst;
  matrix d(n,n);
  
  //  initial values
  for (i=0;i<n;i++){
    eps[i]=Mm->ip[ipp].strain[i];
  }
  
  matstiff (d,ssst);
  mxv (d,eps,sig);

  for (i=0;i<n;i++){
    Mm->ip[ipp].stress[i]=sig[i];
  }
  
}

