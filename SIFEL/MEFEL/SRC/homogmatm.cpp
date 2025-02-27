/*
  File:             homogmatm.cpp
  Author:           Tomas Krejci, 22/07/2022
  Purpose:          Calculates properties of homogenized elastic material
*/ 
#include "mechmat.h"
#include "probdesc.h"
#include "homogmatm.h"
#include "matrix.h"
#include "vector.h"
#include "stochdriver.h"
#include "global.h"
#include "intpoints.h"
#include "vecttens.h"



/**
  Constructor initializes data members to zero or default values.

  Created by TKr,
*/
homogmatm::homogmatm (void)
{
  allocm (6,6,dd);

  hom_mattypem = 0;
  hom_mattypem_number = 0;
}



/**
  Destructor is defined only for the formal purposes.

  Created by TKr,
*/
homogmatm::~homogmatm (void)
{
  destrm(dd);
}



/**
  Function reads material parameters from the opened text file.
   
  @param in - pointer to the opened XFILE

  @return The function does not return anything.

  Created by TKr,
*/
void homogmatm::read (XFILE *in)
{
  xfscanf (in,"%ld %ld",&hom_mattypem,&hom_mattypem_number);
  hom_mattypem_number = hom_mattypem_number-1;
}

/**
  Function prints material parameters into the opened text file.
   
  @param out - pointer to the opened FILE

  @return The function does not return anything.

  Created by TKr
*/
void homogmatm::print (FILE *out)
{
  fprintf (out,"  %ld %ld",hom_mattypem,hom_mattypem_number+1);
}


/**
  Function assembles stiffness %matrix of material.
   
  @param d - stiffness %matrix of material (output)
  @param ssst - strain/stress state

  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKr
*/
void homogmatm::matstiff (matrix &d,strastrestate ssst)
{
  switch (ssst){
  case planestress:{
    matstiff_plstress (d);
    break;
  }
  case planestrain:{
    matstiff_plstrain (d);
    break;
  }
  case axisymm:{
    matstiff_axi (d);
    break;
  }
  case spacestress:{
    matstiff_spacestr (d);
    break;
  }
  default:{
    print_err("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
  }
  }
}


/**
  Function creates stiffness matrix of the elastic
  isotropic material for 2D problems (plane stress).

  @param d - stiffness matrix of the material

  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKr
*/
void homogmatm::matstiff_plstress (matrix &d)
{
  d[0][0] = dd[0][0];  d[0][1] = dd[0][1];  d[0][2] = dd[0][2];
  d[1][0] = dd[1][0];  d[1][1] = dd[1][1];  d[1][2] = dd[1][2];
  d[2][0] = dd[2][0];  d[2][1] = dd[2][1];  d[2][2] = dd[2][2];
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for 2D problems (plane strain).

  @param d - stiffness %matrix of the material

  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKr
*/
void homogmatm::matstiff_plstrain (matrix &d)
{
  d[0][0] = dd[0][0];  d[0][1] = dd[0][1];  d[0][2] = dd[0][2];
  d[1][0] = dd[1][0];  d[1][1] = dd[1][1];  d[1][2] = dd[1][2];
  d[2][0] = dd[2][0];  d[2][1] = dd[2][1];  d[2][2] = dd[2][2];

  // zde rozmyslet
  if (d.m > 3)
    {
      d[0][3] = dd[0][1]; d[1][3] = dd[1][0];
      d[3][0] = dd[1][0]; d[3][1] = dd[1][0]; d[3][3] = dd[1][1];
    }
}



/**
  Function creates stiffness %matrix of the elastic
  isotropic material for 2D axisymmetric problems.

  @param d - stiffness %matrix of the material

  @return The function returns material stiffness %matrix in the parameter d.

  Created by TKr, 19.7.2001
*/
void homogmatm::matstiff_axi (matrix &d)
{
  d[0][0] = dd[0][0];  d[0][1] = dd[0][1];  d[0][2] = dd[0][2];  d[0][3] = dd[0][3];
  d[1][0] = dd[1][0];  d[1][1] = dd[1][1];  d[1][2] = dd[1][2];  d[1][3] = dd[1][3];
  d[2][0] = dd[2][0];  d[2][1] = dd[2][1];  d[2][2] = dd[2][2];  d[2][3] = dd[2][3];
  d[3][0] = dd[3][0];  d[3][1] = dd[3][1];  d[3][2] = dd[3][2];  d[3][3] = dd[3][3];
}


/**
  Function creates stiffness %matrix of the elastic
  isotropic material for 3D problems.
   
  @param d - stiffness %matrix of the material

  Created by TKr
*/
void homogmatm::matstiff_spacestr (matrix &d)
{
  d[0][0] = dd[0][0];  d[0][1] = dd[0][1];  d[0][2] = dd[0][2];  d[0][3] = dd[0][3];  d[0][4] = dd[0][4];  d[0][5] = dd[0][5];
  d[1][0] = dd[1][0];  d[1][1] = dd[1][1];  d[1][2] = dd[1][2];  d[1][3] = dd[1][3];  d[1][4] = dd[1][4];  d[1][5] = dd[1][5];
  d[2][0] = dd[2][0];  d[2][1] = dd[2][1];  d[2][2] = dd[2][2];  d[2][3] = dd[2][3];  d[2][4] = dd[2][4];  d[2][5] = dd[2][5];
  d[3][0] = dd[3][0];  d[3][1] = dd[3][1];  d[3][2] = dd[3][2];  d[3][3] = dd[3][3];  d[3][4] = dd[3][4];  d[3][5] = dd[3][5];
  d[4][0] = dd[4][0];  d[4][1] = dd[4][1];  d[4][2] = dd[4][2];  d[4][3] = dd[4][3];  d[4][4] = dd[4][4];  d[4][5] = dd[4][5];
  d[5][0] = dd[5][0];  d[5][1] = dd[5][1];  d[5][2] = dd[5][2];  d[5][3] = dd[5][3];  d[5][4] = dd[5][4];  d[5][5] = dd[5][5];
}




/**
  Function assembles compliance %matrix of material.
   
  @param c - complience %matrix of material (output)
  @param ssst - strain/stress state

  @return The function returns material complience %matrix in the parameter d.

  Created by TKr
*/
void homogmatm::matcompl (matrix &c,strastrestate ssst)
{
  matrix d;
  
  allocm (6,6,d); 
  
  matstiff (d,ssst);
  
  //c = invm(d);

  destrm(d);
}


/**
   function assembles the conductivity %matrix after homogenization
   
   @param d - array of %matrix entries
   
   Created by TKr
*/
void homogmatm::assemble_matrices (double *d,long ncomp,long dim)
{
  if (ncomp == 3){
    dd[0][0]=d[0];   dd[0][1]=d[1];   dd[0][2]=d[2];
    dd[1][0]=d[3];   dd[1][1]=d[4];   dd[1][2]=d[5];
    dd[2][0]=d[6];   dd[2][1]=d[7];   dd[2][2]=d[8];

    /* fprintf (Out,"\n");
       fprintf (Out,"\n\n d[0] =  %e",d[0]);
       fprintf (Out,"\n\n d[1] =  %e",d[1]);
       fprintf (Out,"\n\n d[2] =  %e",d[2]);
       fprintf (Out,"\n\n d[3] =  %e",d[3]);
       fprintf (Out,"\n\n d[4] =  %e",d[4]);
       fprintf (Out,"\n\n d[5] =  %e",d[5]);
       fprintf (Out,"\n\n d[6] =  %e",d[6]);
       fprintf (Out,"\n\n d[7] =  %e",d[7]);
       fprintf (Out,"\n\n d[8] =  %e",d[8]);
       fprintf (Out,"\n");
       fflush(Out);
    */
  }

  if (ncomp == 4){
    dd[0][0]=d[0];   dd[0][1]=d[1];   dd[0][2]=d[2];   dd[0][3]=d[3];
    dd[1][0]=d[4];   dd[1][1]=d[5];   dd[1][2]=d[6];   dd[1][3]=d[7];
    dd[2][0]=d[8];   dd[2][1]=d[9];   dd[2][2]=d[10];  dd[2][3]=d[11];
    dd[3][0]=d[12];  dd[3][1]=d[13];  dd[3][2]=d[14];  dd[3][3]=d[15];
  }

  if (ncomp == 6){
    dd[0][0]=d[0];   dd[0][1]=d[1];   dd[0][2]=d[2];   dd[0][3]=d[3];   dd[0][4]=d[4];   dd[0][5]=d[5];
    dd[1][0]=d[6];   dd[1][1]=d[7];   dd[1][2]=d[8];   dd[1][3]=d[9];   dd[1][4]=d[10];  dd[1][5]=d[11];
    dd[2][0]=d[12];  dd[2][1]=d[13];  dd[2][2]=d[14];  dd[2][3]=d[15];  dd[2][4]=d[16];  dd[2][5]=d[17];
    dd[3][0]=d[18];  dd[3][1]=d[19];  dd[3][2]=d[20];  dd[3][3]=d[21];  dd[3][4]=d[22];  dd[3][5]=d[23];
    dd[4][0]=d[24];  dd[4][1]=d[25];  dd[4][2]=d[26];  dd[4][3]=d[27];  dd[4][4]=d[28];  dd[4][5]=d[29];
    dd[5][0]=d[30];  dd[5][1]=d[31];  dd[5][2]=d[32];  dd[5][3]=d[33];  dd[5][4]=d[34];  dd[5][5]=d[35];

    /* fprintf (Out,"\n");
       fprintf (Out,"\n\n d[0] =  %e",d[0]);
       fprintf (Out,"\n\n d[1] =  %e",d[1]);
       fprintf (Out,"\n\n d[2] =  %e",d[2]);
       fprintf (Out,"\n\n d[3] =  %e",d[3]);
       fprintf (Out,"\n\n d[4] =  %e",d[4]);
       fprintf (Out,"\n\n d[5] =  %e",d[5]);
       fprintf (Out,"\n\n d[6] =  %e",d[6]);
       fprintf (Out,"\n\n d[7] =  %e",d[7]);
       fprintf (Out,"\n\n d[8] =  %e",d[8]);
       fprintf (Out,"\n\n d[9] =  %e",d[9]);
       fprintf (Out,"\n\n d[10] =  %e",d[10]);
       fprintf (Out,"\n\n d[11] =  %e",d[11]);
       fprintf (Out,"\n\n d[12] =  %e",d[12]);
       fprintf (Out,"\n\n d[13] =  %e",d[13]);
       fprintf (Out,"\n\n d[14] =  %e",d[14]);
       fprintf (Out,"\n\n d[15] =  %e",d[15]);
       fprintf (Out,"\n\n d[16] =  %e",d[16]);
       fprintf (Out,"\n\n d[17] =  %e",d[17]);
       fprintf (Out,"\n\n d[18] =  %e",d[18]);
       fprintf (Out,"\n\n d[19] =  %e",d[19]);
       fprintf (Out,"\n\n d[20] =  %e",d[20]);
       fprintf (Out,"\n\n d[21] =  %e",d[21]);
       fprintf (Out,"\n\n d[22] =  %e",d[22]);
       fprintf (Out,"\n\n d[23] =  %e",d[23]);
       fprintf (Out,"\n\n d[24] =  %e",d[24]);
       fprintf (Out,"\n\n d[25] =  %e",d[25]);
       fprintf (Out,"\n\n d[26] =  %e",d[26]);
       fprintf (Out,"\n\n d[27] =  %e",d[27]);
       fprintf (Out,"\n\n d[28] =  %e",d[28]);
       fprintf (Out,"\n\n d[29] =  %e",d[29]);
       fprintf (Out,"\n\n d[30] =  %e",d[30]);
       fprintf (Out,"\n\n d[31] =  %e",d[31]);
       fprintf (Out,"\n\n d[32] =  %e",d[32]);
       fprintf (Out,"\n\n d[33] =  %e",d[33]);
       fprintf (Out,"\n\n d[34] =  %e",d[34]);
       fprintf (Out,"\n\n d[35] =  %e",d[35]);
       fprintf (Out,"\n");
       fflush(Out);
    */
  }
}


/**
   This function initializes material model data
   with respect of consistency parameter gamma.

   @param ipp - integration point number
   @param ido - index of internal variables for given material in the ipp other array

   Created by TKr
*/
void homogmatm::initval(long ipp, long im, long ido)
{
  //Mm->ip[ipp].other[ido+0] = Mm->ip[ipp].eqother[ido+0] = 0.0;
  //Mm->ip[ipp].other[ido+1] = Mm->ip[ipp].eqother[ido+1] = 0.0;
}


/**
  Function computes true stresses.
   
  @param ipp - number of integration point
   
  @return The function does not return anything.

  Created by TKr
*/
void homogmatm::nlstresses (long ipp, long im, long ido)
{
  long i, n = Mm->ip[ipp].ncompstr;
  vector eps(n),sig(n),epst(ASTCKVEC(6));
  strastrestate ssst = Mm->ip[ipp].ssst;
  matrix d(n,n);
  
  //  initial values
  for (i=0;i<n;i++){
    eps[i]=Mm->ip[ipp].strain[i];
  }
  
  matstiff (d,ssst);
  mxv (d,eps,sig);

  /* //it is necessary to complete here:
    if (Mm->ip[ipp].ssst == planestress)
    Mm->ip[ipp].strain[3] = -nu / (1.0 - nu) * (eps[0]+eps[1]);
  */

  if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5))
  {
    for (i=0;i<n;i++)
      sig(i) += Mm->eigstresses[ipp][i];
  }
  for (i=0;i<n;i++){
    Mm->ip[ipp].stress[i]=sig(i);
  }
}

/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  Created by TKr
*/
void homogmatm::updateval (long ipp,long ido)
{
  //Mm->ip[ipp].eqother[ido+0] = Mm->ip[ipp].other[ido+0];     //
  //Mm->ip[ipp].eqother[ido+1] = Mm->ip[ipp].other[ido+1];     //
}



/**
  The function returns rate of the volumetric strain rate at the given integration point.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns rate of the volumetric strain.
  
  Created by TKr according to Tomas Koudelka
*/
double homogmatm::give_strain_vol(long ipp, long ido)
{
  //return Mm->ip[ipp].eqother[ido+1];
  return 0.0;
}



















