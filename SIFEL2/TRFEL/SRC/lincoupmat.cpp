/*
  File:             lincoupmat.cpp
  Author:           Jaroslav Kruis, 15.11.2008
  Purpose:          Calculates properties of general linear material model for coupled problems
*/ 

#include "lincoupmat.h"
#include "stochdrivert.h"

lincoupmat::lincoupmat (void)
{
  //  number of transported matters
  ntm=0;
  //  geometrical dimension of the problem
  dim=0;
  
  //  array of conductivity coefficients
  k=NULL;
  //  array of capacity coefficients
  c=NULL;
}
lincoupmat::~lincoupmat (void)
{
  delete [] k;
  delete [] c;
}



/**
   function reads input and material parameters
   
   @param in - input file
   
   JK, 15.11.2008
*/
void lincoupmat::read (XFILE *in)
{
  long i,n;
  
  //  number of transported matters
  xfscanf (in,"%k%ld","number_of_matters",&ntm);
  //  geometrical dimension of the problem
  xfscanf (in,"%k%ld","geom_dimension",&dim);
  
  if (ntm<1)
    print_err("number of transported matters is less than 1",__FILE__,__LINE__,__func__);
  if (dim<1)
    print_err("geometrical dimension is less than 1",__FILE__,__LINE__,__func__);
  if (dim<3)
    print_err("geometrical dimension is greater than 3",__FILE__,__LINE__,__func__);
  
  if (k!=NULL)
    delete [] k;
  if (c!=NULL)
    delete [] c;
  
  //  number of conductivity (capacaity) coefficients
  n = ntm*ntm*dim*dim;

  //  array of conductivity coefficients
  k = new double [n];
  //  array of capacity coefficients
  c = new double [n];
  
  //  reading of conductivity coefficients
  for (i=0;i<n;i++){
    xfscanf (in,"%le",k+i);
  }

  //  reading of capacity coefficients
  for (i=0;i<n;i++){
    xfscanf (in,"%le",c+i);
  }

}


/**
   function prints input and material parameters
   
   @param out - outut file
   
   JK, 16.11.2008
*/
void lincoupmat::print (FILE *out)
{
  long i,n;
  
  //  number of transported matters
  fprintf (out,"\n%ld",ntm);
  //  geometrical dimension of the problem
  fprintf (out," %ld\n",dim);

  //  number of conductivity (capacaity) coefficients
  n = ntm*ntm*dim*dim;

  //  printing of conductivity coefficients
  for (i=0;i<n;i++){
    fprintf (out,"%le ",k[i]);
  }
  fprintf (out,"\n");
  //  printing of capacity coefficients
  for (i=0;i<n;i++){
    fprintf (out,"%le ",c[i]);
  }
  fprintf (out,"\n");

}




/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   JK, 15.11.2008
*/
void lincoupmat::matcond (matrix &d,long ri,long ci,long ipp)
{
  switch (d.n){
  case 1:{
    matcond1d (d,ri,ci,ipp);//1D
    break;
  }
  case 2:{
    matcond2d (d,ri,ci,ipp);//2D
    break;
  }
  case 3:{
    matcond3d (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    print_err("unknown geometrical dimension of the problem is required",__FILE__,__LINE__,__func__);
  }
  }
}


/**
   function assembles conductivity %matrix of the material for 1D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index of block
   @param ci  - column index of block
   @param ipp - number of integration point
   
   JK, 16.11.2008
*/
void lincoupmat::matcond1d (matrix &d,long ri,long ci,long /*ipp*/)
{
  d[0][0] = give_k (ri,ci,0,0);
}

/**
   function assembles conductivity %matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index of block
   @param ci  - column index of block
   @param ipp - number of integration point
   
   JK, 16.11.2008
*/
void lincoupmat::matcond2d (matrix &d,long ri,long ci,long /*ipp*/)
{
  fillm(0.0,d);
  
  d[0][0] = give_k (ri,ci,0,0);
  d[0][1] = give_k (ri,ci,0,1);

  d[1][0] = give_k (ri,ci,1,0);
  d[1][1] = give_k (ri,ci,1,1);
}

/**
   function assembles conductivity %matrix of the material for 3D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index of block
   @param ci  - column index of block
   @param ipp - number of integration point
   
   JK, 16.11.2008
*/
void lincoupmat::matcond3d (matrix &d,long ri,long ci,long /*ipp*/)
{
  fillm(0.0,d);
  
  d[0][0] = give_k (ri,ci,0,0);
  d[0][1] = give_k (ri,ci,0,1);
  d[0][2] = give_k (ri,ci,0,2);

  d[1][0] = give_k (ri,ci,1,0);
  d[1][1] = give_k (ri,ci,1,1);
  d[1][2] = give_k (ri,ci,1,2);

  d[2][0] = give_k (ri,ci,2,0);
  d[2][1] = give_k (ri,ci,2,1);
  d[2][2] = give_k (ri,ci,2,2);
}


/**
   function assembles capacity %matrix of the material

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   JK, 16.11.2008
*/
void lincoupmat::matcap (double &cc,long ri,long ci,long /*ipp*/)
{
  cc = give_c (ri,ci,0,0);
}




/**
   function returns conductivity coefficient
   
   @param bri - row index of the block
   @param bci - column index of the block
   @param ri - row index of coefficient in the required block
   @param ci - column index of coefficient in the required block
   
   JK, 16.11.2008
*/
double lincoupmat::give_k (long bri,long bci, long ri,long ci)
{
  return (k[bri*ntm*dim*dim+bci*dim*dim+ri*dim+ci]);
}

/**
   function returns capacity coefficient
   
   @param bri - row index of the block
   @param bci - column index of the block
   @param ri - row index of coefficient in the required block
   @param ci - column index of coefficient in the required block
   
   JK, 16.11.2008
*/
double lincoupmat::give_c (long bri,long bci, long ri,long ci)
{
  return (c[bri*ntm*dim*dim+bci*dim*dim+ri*dim+ci]);
}

