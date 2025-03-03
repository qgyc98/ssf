#include "nlisotrmat.h"
#include "stochdrivert.h"
#include "globalt.h"

nlisotrmat::nlisotrmat (void)
{
  a=0.0;
  b=0.0;
  c=0.0;
/*
  omega=0.0;
  lambda=0.0;
  epsilon=0.0;
  sigma=0.0;
  w=0.0;
  */
}
nlisotrmat::~nlisotrmat (void)
{}


/**
   function reads material parameters
   
   @param in - input file
   
   JK, 18.10.2007
*/
void nlisotrmat::read (XFILE *in)
{
  //xfscanf (in,"%lf %lf %lf",&a,&b,&c);
  //xfscanf (in,"%lf %lf %lf %lf %lf",&omega,&lambda,&epsilon,&sigma,&w);
  xfscanf (in,"%lf",&c); //c=capacity
  k.read(in); 	    // k=conductivity
}

/**
   function prints material parameters
   
   @param out - outut file
   
   JK, 18.10.2007
*/
void nlisotrmat::print (FILE *out)
{
  fprintf (out,"  %le %le %le",a,b,c);
}

/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void nlisotrmat::matcond (matrix &d,long ri,long ci,long ipp)
{
  long n;
  n = d.n;
  
  switch (n){
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
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


/**
   function creates conductivity %matrix of the material for 1D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void nlisotrmat::matcond1d (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double v,kk;
  
  //  actual value
  v = Tm->ip[ipp].av[0];
  //if (v<1.0e-6)
  //v=280.0;
  //  actual coefficient of conductivity 
  kk = get_k (v);
  
  d[0][0] = kk;
}

/**
   function creates conductivity %matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void nlisotrmat::matcond2d (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double v,kk;

  //  actual value
  v = Tm->ip[ipp].av[0];
  //  actual coefficient of conductivity 
  kk = get_k (v);
  
  fillm(0.0,d);
  
  d[0][0] = kk;   d[0][1] = 0.0;
  d[1][0] = 0.0;  d[1][1] = kk;
}

/**
   function creates conductivity %matrix of the material for 3D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/

void nlisotrmat::matcond3d (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double v,kk;

  //  actual value
  v = Tm->ip[ipp].av[0];
  //  actual coefficient of conductivity 
  kk = get_k (v);
  
  fillm(0.0,d);
  
  d[0][0]=kk;   d[0][1]=0.0;  d[0][2]=0.0;
  d[1][0]=0.0;  d[1][1]=kk;   d[1][2]=0.0;
  d[2][0]=0.0;  d[2][1]=0.0;  d[2][2]=kk;
}


/**
  The function returns ordered dof names of primary unknowns 
  required by the model.

  @param dofname   - array of uknown name for particular nodal dofs (output)
                     dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
  @param ntm       - number of transported media = number of nodal dof = length of array dofname

  Created by Tomas Koudelka, 6.12.2013
*/

void nlisotrmat::give_dof_names(namevart *dofname, long ntm)
{
  if (ntm < 1)
  {
    print_err("the model defines %ld unknowns but number of transported media is %ld", __FILE__, __LINE__, __func__, 1, ntm);
    abort();
  }
  dofname[0] = trf_temperature;
}






/**
   function creates capacity %matrix of the material

   @param c   - capacity matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void nlisotrmat::matcap (double &cc,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  cc = 0.0;
  
  cc = get_c();
}




/**
   function creates conductivity coefficient of the isotropic material

   @param v - actual value of solved variable
   @retval k - conductivity coefficient of the isotropic material
   
   JK, 18.10.2007
*/

double nlisotrmat::get_k (double v)
{
  //return (a*(v-280.0)*(v-280.0)+b);
  //return (a+v/10.0);
  //return (a*v+b);
  //return a;
  //return ((1.0-omega)*lambda+4.0*omega*epsilon*sigma*w*v*v*v);
  
  double lambda;
  lambda= k.getval(v);
  //fprintf(Outt, " teplota %le  lambda %le", v, lambda);
  return lambda;
  
}

/**
   function creates specific heat of the isotropic material
   
   @retval c - specific heat of the isotropic material
*/
double nlisotrmat::get_c()
{
  return(c);
}


/**
   function changes parameters of conductivity and capacity from a table
   @param
*/
/*
void nlisotrmat::changeparam (atsel &atm,vector &val)
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
    default:{
      fprintf (stderr,"\n\n wrong number of atribute in function changeparam (%s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
}
*/

/**
  The function returns ordered dof names of primary unknowns 
  required by the model.

  @param dofname   - array of uknown name for particular nodal dofs (output)
                     dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
  @param ntm       - number of transported media = number of nodal dof = length of array dofname

  JK, 8. 9. 2014
*/
/*
void nlisotrmat::give_dof_names (namevart *dofname, long ntm)
{
  if (ntm < 1)
  {
    print_err("the model defines %ld unknowns but number of transported media is %ld", 
              __FILE__, __LINE__, __func__, 1, ntm);
    abort();
  }
  dofname[0] = trf_temperature;
}
*/
