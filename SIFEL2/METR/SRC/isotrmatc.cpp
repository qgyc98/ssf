#include "stochdrivert.h"
#include "constrel.h"
#include "coupmatu.h"
#include "globalc.h"
#include "global.h"
#include "globalt.h"
#include "isotrmatc.h"

isotrmatc::isotrmatc (void)
{
  k=0.0;  c=0.0;
}

isotrmatc::~isotrmatc (void)
{

}


/**
   function reads parameters
   
   @param in - input file

   TKr 17/07/2018
*/
void isotrmatc::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %lf %lf",&c,&k,&e,&nu,&alpha);
}


/**
   function prints parameters
   
   @param out - output file

   TKr 17/07/2018
*/
void isotrmatc::print(FILE *out)
{
  fprintf (out,"\n %lf %lf %lf %lf %lf",c,k,e,nu,alpha);
}


/**
   function computes conductivity %matrix of the material
   in the required integration point for upper block
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17/07/2018, TKr
*/
void isotrmatc::matcond_u (matrix &d,long ri,long ci,long ipp)
{
  long m;
  m = d.m;

  switch (m){
  case 1:{
    matcond1d_u (d,ri,ci,ipp);//1D
    break;
  }
  case 3:{
    matcond2d_u (d,ri,ci,ipp);//2D
    break;
  }
  case 4:{
    matcond2d_ax_u (d,ri,ci,ipp);//2D - axisymmetric
    break;
  }   
  case 6:{
    matcond3d_u (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    print_err("unknown number of components of stress tensor is required",__FILE__,__LINE__,__func__);
  }
  }
  
  if (ci == 0){
    matrix s(d.m,d.m);
    
    Cmu->matstiff (s,ipp);//only for temperature
    mxm(s,d,d);
    
    destrm (s);
  }
}


/**
   function creates conductivity %matrix of the material for 1D problems
   for upper block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17/07/2018, TKr
*/
void isotrmatc::matcond1d_u (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double k;
  k = 0.0;
  
  k = -get_alpha();
  
  fillm(0.0,d);
  
  d[0][0] = k;
}

/**
   function creates conductivity %matrix of the material for 2D problems
   for upper block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17/07/2018, TKr
*/
void isotrmatc::matcond2d_u (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double k;
  k = 0.0;
  
  k = -get_alpha();
  
  fillm(0.0,d);

  d[0][0] = k;
  d[1][0] = k;
  d[2][0] = 0.0;  
}

/**
   function creates conductivity %matrix of the material for 2D axisymmetric problems
   for upper block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void isotrmatc::matcond2d_ax_u (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double k;
  k = 0.0;
  
  k = -get_alpha();
  
  fillm(0.0,d);

  d[0][0] = k;
  d[1][0] = k;
  d[2][0] = k;  
  d[3][0] = 0.0; 
}

/**
   function creates conductivity %matrix of the material for 3D problems
   for upper block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17/07/2018, TKr
*/
void isotrmatc::matcond3d_u (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double k;
  k = 0.0;
  
  k = -get_alpha();
  
  fillm(0.0,d);
  
  d[0][0]=k;  
  d[1][0]=k;   
  d[2][0]=k;
  d[3][0]=0.0;   
  d[4][0]=0.0;  
  d[5][0]=0.0;
}


/**
   function computes capacity %matrix of the material
   in the required integration point for upper block
   
   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   17/07/2018, TKr
*/
void isotrmatc::matcap_u (matrix &d,long ri,long ci,long ipp)
{
  long m;

  m = d.m;

  switch (m){
  case 1:{
    matcap1d_u (d,ri,ci,ipp);//1D
    break;
  }
  case 3:{
    matcap2d_u (d,ri,ci,ipp);//2D
    break;
  }
  case 4:{
    matcap2d_ax_u (d,ri,ci,ipp);//2D - axisymmetric
    break;
  }   
  case 6:{
    matcap3d_u (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    print_err("unknown number of components of stress tensor is required",__FILE__,__LINE__,__func__);
  }
  }
}


/**
   function creates capacity %matrix of the material for 1D problems
   for upper block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17/07/2018, TKr
*/
void isotrmatc::matcap1d_u (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double c;
  
  c = 0.0;
  
  d[0][0] = c;
}


/**
   function creates capacity %matrix of the material for 2D problems
   for upper block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17/07/2018, TKr
*/
void isotrmatc::matcap2d_u (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double c;
  
  c = 0.0;
  
  fillm(0.0,d);

  d[0][0] = c;
  d[1][0] = c;
  d[2][0] = 0.0;
}



/**
   function creates capacity %matrix of the material for 2D axisymmetric problems
   for upper block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void isotrmatc::matcap2d_ax_u (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double c;
  
  c = 0.0;
  
  fillm(0.0,d);

  d[0][0] = c;
  d[1][0] = c;
  d[2][0] = c;
  d[3][0] = 0.0;
}


/**
   function creates capacity %matrix of the material for 3D problems
   for upper block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17/07/2018, TKr
*/
void isotrmatc::matcap3d_u (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double c;
  
  c = 0.0;
  
  fillm(0.0,d);

  d[0][0] = c;
  d[1][0] = c;
  d[2][0] = c;
  d[3][0] = 0.0;
  d[4][0] = 0.0;
  d[5][0] = 0.0;
}


/**
   function computes conductivity %matrix of the material
   in the required integration point for lower block
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17/07/2018, TKr
*/
void isotrmatc::matcond_l (matrix &d,long ri,long ci,long ipp)
{
  long m;
  m = d.n;

  switch (m){
  case 1:{
    matcond1d_l (d,ri,ci,ipp);//1D
    break;
  }
  case 3:{
    matcond2d_l (d,ri,ci,ipp);//2D
    break;
  }
  case 4:{
    matcond2d_ax_l (d,ri,ci,ipp);//2D - axisymmetric
    break;
  }   
  case 6:{
    matcond3d_l (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    print_err("unknown number of components of stress tensor is required",__FILE__,__LINE__,__func__);
  }
  }
}


/**
   function creates conductivity %matrix of the material for 1D problems
   for lower block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17/07/2018, TKr
*/
void isotrmatc::matcond1d_l (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double k;

  k = 0.0;
  
  d[0][0] = k;
}

/**
   function creates conductivity %matrix of the material for 2D problems
   for lower block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17/07/2018, TKr
*/
void isotrmatc::matcond2d_l (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double k;
  k = 0.0;
  
  fillm(0.0,d);
  
  d[0][0] = k;
  d[0][1] = k;
  d[0][2] = 0.0;  
}

/**
   function creates conductivity %matrix of the material for 2D axisymmetric problems
   for lower block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void isotrmatc::matcond2d_ax_l (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double k;
  k = 0.0;
  
  fillm(0.0,d);
  
  d[0][0] = k;
  d[0][1] = k;
  d[0][2] = k;  
  d[0][3] = 0.0;  
}

/**
   function creates conductivity %matrix of the material for 3D problems
   for lower block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17/07/2018, TKr
*/
void isotrmatc::matcond3d_l (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double k;
  k = 0.0;
  
  fillm(0.0,d);
  
  d[0][0]=k;  
  d[0][1]=k;   
  d[0][2]=k;
  d[0][3]=0.0;   
  d[0][4]=0.0;  
  d[0][5]=0.0;
}


/**
   function computes capacity %matrix of the material
   in the required integration point for lower block
   
   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   17/07/2018, TKr
*/
void isotrmatc::matcap_l (matrix &d,long ri,long ci,long ipp)
{
  long m;

  m = d.n;

  switch (m){
  case 1:{
    matcap1d_l (d,ri,ci,ipp);//1D
    break;
  }
  case 3:{
    matcap2d_l (d,ri,ci,ipp);//2D
    break;
  }
  case 4:{
    matcap2d_ax_l (d,ri,ci,ipp);//2D - axisymmetric
    break;
  }   
  case 6:{
    matcap3d_l (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    print_err("unknown number of components of stress tensor is required",__FILE__,__LINE__,__func__);
  }
  }
}


/**
   function creates capacity %matrix of the material for 1D problems
   for lower block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17/07/2018, TKr
*/
void isotrmatc::matcap1d_l (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double c;

  c = 0.0;
  
  d[0][0] = c;
}


/**
   function creates capacity %matrix of the material for 2D problems
   for lower block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17/07/2018, TKr
*/
void isotrmatc::matcap2d_l (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double c;

  c = 0.0;
  
  fillm(0.0,d);

  
  d[0][0] = c;
  d[0][1] = c;
  d[0][2] = 0.0;  
}



/**
   function creates capacity %matrix of the material for 2D axisymmetric problems
   for lower block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void isotrmatc::matcap2d_ax_l (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double c;

  c = 0.0;
  
  fillm(0.0,d);

  
  d[0][0] = c;
  d[0][1] = c;
  d[0][2] = c;  
  d[0][3] = 0.0;  
}


/**
   function creates capacity %matrix of the material for 3D problems
   for lower block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17/07/2018, TKr
*/
void isotrmatc::matcap3d_l (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double c;

  c = 0.0;
  
  fillm(0.0,d);
  
  d[0][0]=c;  
  d[0][1]=c;   
  d[0][2]=c;
  d[0][3]=0.0;   
  d[0][4]=0.0;  
  d[0][5]=0.0;
}


/**
   function computes volume part of right-hand side %matrix
   in the required integration point
   
   @param d - right-hand side %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void isotrmatc::rhs_u2 (matrix &d,long ri,long ci,long ipp)
{
  long m;
  m = d.m;
  
  switch (m){
  case 1:{
    rhs1d2 (d,ri,ci,ipp);//1D
    break;
  }
  case 3:{
    rhs2d2 (d,ri,ci,ipp);//2D
    break;
  }
  case 4:{
    rhs2d_ax_2 (d,ri,ci,ipp);//2D axisymmetric
    break;
  }
  case 6:{
    rhs3d2 (d,ri,ci,ipp);//3D      
    break;
  }
  default:{
    print_err("unknown number of components of stress tensor is required",__FILE__,__LINE__,__func__);
  }
  }
  
  if (ci == 0){
    matrix s(d.m,d.m);
    
    Cmu->matstiff (s,ipp);//only for temperature
    mxm(s,d,d);
    
    destrm (s);
  }
}




/**
   function creates volume right-hand side %matrix of the material for 1D problems

   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void isotrmatc::rhs1d2 (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double f;
  
  f = -get_alpha();
  
  fillm(0.0,d);

  d[0][0] = f;
}

/**
   function creates volume right-hand side %matrix of the material for 2D problems

   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void isotrmatc::rhs2d2 (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double f;
  
  f = -get_alpha();
  
  fillm(0.0,d);

  d[0][0] = f;
  d[1][0] = f;
  d[2][0] = 0.0;
}


/**
   function creates volume right-hand side %matrix of the material for 2D axisymmetric problems

   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void isotrmatc::rhs2d_ax_2 (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double f;
  
  f = -get_alpha();
  
  fillm(0.0,d);

  d[0][0] = f;
  d[1][0] = f;
  d[2][0] = f;
  d[3][0] = 0.0;
}

/**
   function creates volume right-hand side %matrix of the material for 3D problems

   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void isotrmatc::rhs3d2 (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  double f;
  
  f = -get_alpha();
  
  fillm(0.0,d);

  d[0][0]=f;  
  d[1][0]=f;   
  d[2][0]=f;
  d[3][0]=0.0;   
  d[4][0]=0.0;  
  d[5][0]=0.0;
}




/**
   function creates conductivity %matrix of the
   isotropic material

   @retval k - heat conductivity %matrix of the isotropic
*/
double isotrmatc::get_k()
{
  return(k);
}

/**
   function creates capacity %matrix of the
   isotropic material
   
   @retval c - heat capacity %matrix of the isotropic
*/
double isotrmatc::get_c()
{
  return(c);
}

/**
   function creates Young's modulus
   isotropic material
   
   @retval e - Young's modulus
*/
double isotrmatc::get_e()
{
  return(e);
}

/**
   function creates Young's modulus
   isotropic material
   
   @retval nu - Poisson's constant
*/
double isotrmatc::get_nu()
{
  return(nu);
}

/**
   function creates thermal dilatation coefficient
   isotropic material
   
   @retval alpha - thermal dilatation coefficient
*/
double isotrmatc::get_alpha()
{
  return(alpha);
}

void isotrmatc::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=0;i<atm.num;i++){
    switch (atm.atrib[i]){
    case 0:{
      c=val[0];
      break;
    }
    case 1:{
      k=val[1];
      break;
    }
    case 2:{
      alpha=val[2];
      break;
    }
    default:{
      print_err("wrong number of atribute",__FILE__,__LINE__,__func__);
    }
    }
  }
}
