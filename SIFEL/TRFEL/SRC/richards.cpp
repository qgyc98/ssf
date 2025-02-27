/*
    File:             richards.cpp
    Author:          
    Purpose:         
    Source:          
    Assumptions:


*/


#include "richards.h"
#include "globalt.h"


richards::richards()
{
  //  parameters of van Genuchten model
  alpha=0.0;
  n=0.0;
  m=0.0;
  //  dimension of problem solved
  dim=0;
  
  //  saturated hydraulic conductivity
  kksxx=0.0;
  kksxy=0.0;
  kksxz=0.0;
  kksyy=0.0;
  kksyz=0.0;
  kkszz=0.0;

  //  saturated water content
  thetas=0.0;
  //  residual water content
  thetar=0.0;
  
}

richards::~richards()
{

}

/**
   function reads material parameters
   
   @param in - input file
   
   1.12.2011
*/
void richards::read(XFILE *in)
{
    //  dimenison of problem solved
  //  the dimension is determined from the first element
  dim = Tt->give_dimension (0);
  if (dim<1 || dim>3){
    print_err("wrong dimension of problem solved (%ld)",__FILE__,__LINE__,__func__,dim);
  }
  
  switch (dim) {
    case 1:
      xfscanf (in,"%lf",&kksxx);
      break;
    case 2: 
       xfscanf (in,"%lf %lf %lf",&kksxx,&kksxy,&kksyy);
       break;
    case 3:
       xfscanf (in,"%lf %lf %lf %lf %lf %lf",&kksxx,&kksxy,&kksxz,&kksyy,&kksyz,&kkszz);
       break;    
  }
  //  parameters of van Genuchten model
  //  saturated hydraulic conductivity
  //  saturated water content
  //  residual water content
  //  specific storage
  xfscanf (in,"%lf %lf %lf %lf %lf %lf",&alpha,&n,&m,&thetas,&thetar,&storage);
}

/**
   function prints material parameters
   
   @param out - output file
   
   1.12.2011
*/
void richards::print(FILE *out)
{
  switch (dim) {
    case 1:
      fprintf (out,"%le",kksxx);
      break;
    case 2: 
      fprintf (out,"%le %le %le",kksxx,kksxy,kksyy);
       break;
    case 3:
      fprintf (out,"%le %le %le %le %le %le",kksxx,kksxy,kksxz,kksyy,kksyz,kkszz);
       break;    
  }

  //  parameters of van Genuchten model
  //  saturated hydraulic conductivity
  //  saturated water content
  //  residual water content
  fprintf (out," %le %le %le %le  %le\n",alpha,n,m,thetas,thetar);
}

/**
   function computes hydraulic conductivity
   
   @param h - hydraulic head
   
   1.12.2011
*/
double richards::kk_value (double h)
{
  if (h < 0){
    return (pow(1 - pow(-(alpha*h),m*n)/pow(1 + pow(-(alpha*h),n),m),2)/pow(1 + pow(-(alpha*h),n),m/2.0));
  }
  else{
    return 1.0;
  }

}

/**
   function computes the derivative of K_r (relative hydraulic conductivity) with respect to h (dK/dh van Genuchten).
   
   @param h - hydraulic head
   
   1.12.2011
   modified 24.3.2014
*/
double richards::dkkdh_value (double h)
{
  
  if (h < 0){
    return ((alpha*pow(-(alpha*h),-1 + n)*pow(1 + pow(-(alpha*h),n),-1 - m/2.)*pow(1 - pow(-(alpha*h),m*n)/pow(1 + pow(-(alpha*h),n),m),2)*m*n)/2. + (2*(1 - pow(-(alpha*h),m*n)/pow(1 + pow(-(alpha*h),n),m))*(-(alpha*pow(-(alpha*h),-1 + n + m*n)*pow(1 + pow(-(alpha*h),n),-1 - m)*m*n) + (alpha*pow(-(alpha*h),-1 + m*n)*m*n)/pow(1 + pow(-(alpha*h),n),m)))/pow(1 + pow(-(alpha*h),n),m/2.));
  }
  else{
    return 0.0;
  }
  
}


/**
   function computes second derivative of K with respect to h (d^2K/dh^2 van Genuchten).
   
   (possibly not needed, consider remove) 
   
   @param h - hydraulic head
   
   1.12.2011
*/
// double richards::ddkkddh_value (double h)
// {
//   if (h < 0){
//     return kkszz*( -(pow(alpha,2)*pow(-(alpha*h),-2 + n)*pow(1 + pow(-(alpha*h),n),-1 - m/2.)*
// 	  pow(1 - pow(-(alpha*h),m*n)/
// 	  pow(1 + pow(-(alpha*h),n),m),2)*m*(-1 + n)*n)/2.0
// 	    - (pow(alpha,2)*pow(-(alpha*h),-2 + 2*n)*
// 	    pow(1 + pow(-(alpha*h),n),-2 - m/2.)*
// 	  pow(1 - pow(-(alpha*h),m*n)/
// 	      pow(1 + pow(-(alpha*h),n),m),2)*(-1 - m/2.)*m*
// 	  pow(n,2))/2. + 2*alpha*pow(-(alpha*h),-1 + n)*
// 	    pow(1 + pow(-(alpha*h),n),-1 - m/2.)*
// 	  (1 - pow(-(alpha*h),m*n)/
// 	    pow(1 + pow(-(alpha*h),n),m))*m*n*
// 	  (-(alpha*pow(-(alpha*h),-1 + n + m*n)*
// 	    pow(1 + pow(-(alpha*h),n),-1 - m)*m*n) + 
// 	  (alpha*pow(-(alpha*h),-1 + m*n)*m*n)/
// 	  pow(1 + pow(-(alpha*h),n),m)) + 
// 	  (2*pow(-(alpha*pow(-(alpha*h),-1 + n + m*n)*
// 	      pow(1 + pow(-(alpha*h),n),-1 - m)*m*n) + 
// 	    (alpha*pow(-(alpha*h),-1 + m*n)*m*n)/
// 	    pow(1 + pow(-(alpha*h),n),m),2))/
// 	    pow(1 + pow(-(alpha*h),n),m/2.) + 
// 	  (2*(1 - pow(-(alpha*h),m*n)/
// 	    pow(1 + pow(-(alpha*h),n),m))*
// 	    (pow(alpha,2)*pow(-(alpha*h),-2 + 2*n + m*n)*
// 	    pow(1 + pow(-(alpha*h),n),-2 - m)*(-1 - m)*m*
// 	    pow(n,2) + pow(alpha,2)*
// 	    pow(-(alpha*h),-2 + n + m*n)*
// 	    pow(1 + pow(-(alpha*h),n),-1 - m)*pow(m,2)*
// 	    pow(n,2) - (pow(alpha,2)*
// 	      pow(-(alpha*h),-2 + m*n)*m*n*(-1 + m*n))/
// 	    pow(1 + pow(-(alpha*h),n),m) + 
// 	    pow(alpha,2)*pow(-(alpha*h),-2 + n + m*n)*
// 	    pow(1 + pow(-(alpha*h),n),-1 - m)*m*n*
// 	    (-1 + n + m*n)))/pow(1 + pow(-(alpha*h),n),m/2.)) ;
//   }
//   else{
//     return 0.0;
//   }
// 
// }

/**
   function computes water content (van Genuchten).
   
   @param h - hydraulic head
   
   27.3.2014
*/
double richards::theta_val (double h)
{
  if (h < 0){
    return (thetas - thetar)/pow(1 + pow(-(alpha*h),n),m) + thetar ; 
  }
  else{
    return thetas;
  }
}

/**
   function Darcian flux (Darcy-Buckingham law).
   
   @param h - hydraulic head
   
   27.3.2014
*/
vector richards::darcian_flux (long ipp)
{
  double h;
  vector gradH(dim);
  vector flux(dim);
  matrix K(dim,dim);
  int i, j;
  double Kr;
  
  h = Tm->ip[ipp].av[0] ;
  
  Kr = kk_value(h) ;
  
  switch  (dim) {
    case 1:
      K(0,0) = kksxx;
      break ; 
    case 2:
      K(0,0) = kksxx;
      K(0,1) = kksxy;
      K(1,0) = kksxy;
      K(1,1) = kksyy;
      break;
    case 3:
      K(0,0) = kksxx;
      K(0,1) = kksxy;
      K(0,2) = kksxz;
      K(1,0) = kksxy;
      K(1,1) = kksyy;
      K(1,2) = kksyz;
      K(2,0) = kksxz;
      K(2,1) = kksyz;
      K(2,2) = kkszz;
      break;
    }
  
  
  for (i = 1; i <= dim; i++) {
    for (j = 1; j <= dim; j++) {
      K(i-1,j-1) = K(i-1,j-1)*Kr*(-1.0) ;
    }
  }
  
  for (i = 1; i <= dim; i++) {
    gradH(i-1) = Tm->ip[ipp].grad[0][i-1];
  }
  
  mxv(K, gradH, flux) ; 
  
  return flux ;
  
  
}



/**
   function computes capacity (van Genuchten).
   
   @param h - hydraulic head
   
   1.12.2011
   modified 27.3.2014
*/
double richards::c_value (double h)
{
  if (h < 0){
    return alpha*pow(-(alpha*h),-1 + n)*pow(1 + pow(-(alpha*h),n),-1 - m)*m*n*(thetas - thetar) + 
      storage*theta_val(h)/thetas;
  }
  else{
    return storage;
  }
}



/**
   function computes derivative of capacity with respect to hydraulic head (van Genuchten).
   
   @param h - hydraulic head
   
   1.12.2011
*/
double richards::dcdh_value (double h)
{
  if (h < 0){
    return -(pow(alpha,2)*pow(-(alpha*h),-2 + n)*
	     pow(1 + pow(-(alpha*h),n),-1 - m)*m*(-1 + n)*n*
	     (thetas - thetar)) - pow(alpha,2)*pow(-(alpha*h),-2 + 2*n)*
      pow(1 + pow(-(alpha*h),n),-2 - m)*(-1 - m)*m*pow(n,2)*
      (thetas - thetar) + (alpha*pow(-(alpha*h),-1 + n)*
			     pow(1 + pow(-(alpha*h),n),-1 - m)*m*n*storage*
			     (thetas - thetar))/thetas;
  }
  else{
    return 0.0; 
  }
}


/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void richards::matcond (matrix &d,long ri,long ci,long ipp)
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
    print_err("unknown number of flux components is required",__FILE__,__LINE__,__func__);
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
void richards::matcond1d (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double h,k;
  
  //  actual hydraulic head
  h = Tm->ip[ipp].av[0];
  
  //  hydraulic conductivity
  k = kk_value (h);
  
  d[0][0] = k*kksxx;
}

/**
   function creates conductivity %matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void richards::matcond2d (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double h,k;
  
  //  actual hydraulic head
  h = Tm->ip[ipp].av[0];
  
  //  hydraulic conductivity
  k = kk_value (h);
  
  fillm(0.0,d);
  
  d[0][0] = k*kksxx;   d[0][1] = k*kksxy;
  d[1][0] = k*kksxy;   d[1][1] = k*kksyy;
  //fprintf (stdout,"\n vodivost  %le",d[0][0]);

}

/**
   function creates conductivity %matrix of the material for 3D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/

void richards::matcond3d (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double h,k;
  
  //  actual hydraulic head
  h = Tm->ip[ipp].av[0];
  
  //  hydraulic conductivity
  k = kk_value (h);
  
  fillm(0.0,d);
  
  d[0][0] = k*kksxx;   d[0][1] = k*kksxy;  d[0][2]=k*kksxz;
  d[1][0] = k*kksxy;   d[1][1] = k*kksyy;  d[1][2]=k*kksyz;
  d[2][0] = k*kksxz;   d[2][1] = k*kksyz;  d[2][2]=k*kkszz;
}

void richards::matcond2 (matrix &d,long ri,long ci,long ipp)
{
  long n;
  n = d.n;
  
  switch (n){
  case 1:{
    matcond1d_2 (d,ri,ci,ipp);//1D
    break;
  }
  case 2:{
    matcond2d_2 (d,ri,ci,ipp);//2D
    break;
  }
  case 3:{
    matcond3d_2 (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    print_err("unknown number of components is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function creates conductivity %matrix of the material for 1D problems (convective term)
   
   this %matrix is correct only for homogeneous material on finite elements
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void richards::matcond1d_2 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double h;
  
  //  actual hydraulic head
  h = Tm->ip[ipp].av[0];
  
  fillm(0.0,d);
  d[0][0]=0.0-kksxx*dkkdh_value (h);
}


/**
   function creates conductivity %matrix of the material for 1D problems (convective term)

   this %matrix is correct only for homogeneous material on finite elements
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
*/
void richards::matcond2d_2 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double h;
  
  //  actual hydraulic head
  h = Tm->ip[ipp].av[0];
  
  fillm(0.0,d);
  d[0][0]=0.0-kksxy*dkkdh_value (h);
  d[0][1]=0.0-kksyy*dkkdh_value (h);
  
  //fprintf (stdout,"\n konvekce %le",d[0][1]);
}

/**
   function creates conductivity %matrix of the material for 1D problems (convective term)

   this %matrix is correct only for homogeneous material on finite elements
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void richards::matcond3d_2 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double h;
  
  //  actual hydraulic head
  h = Tm->ip[ipp].av[0];
  
  fillm(0.0,d);
  d[0][0]=0.0 - kksxz*dkkdh_value (h);
  d[0][1]=0.0 - kksyz*dkkdh_value (h);
  d[0][2]=0.0 - kkszz*dkkdh_value (h);
}


/**
   function creates capacity %matrix of the material

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void richards::matcap (double &c,long /*ri*/,long /*ci*/,long ipp)
{
  double h;
  
  //  actual hydraulic head
  h = Tm->ip[ipp].av[0];
  
  //  capacity
  c = c_value (h);

}



/**
  The function returns ordered dof names of primary unknowns 
  required by the model.

  @param dofname   - array of uknown name for particular nodal dofs (output)
                     dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
  @param ntm       - number of transported media = number of nodal dof = length of array dofname

  10. 3. 2014, JK
*/
void richards::give_dof_names(namevart *dofname, long ntm)
{
  if (ntm < 1)
  {
    print_err("the model defines %ld unknowns but number of transported media is %ld", 
              __FILE__, __LINE__, __func__, 1, ntm);
    abort();
  }
  dofname[0] = trf_hydraulic_head;
}


/**
  This function computes first order Taylor serie error of the solution over time. 
  \f[ h_n \f] actual hydraulic head
  \f[ h_p \f] previous time level hydraulic head
  the error is defined as
  \f[ \epsilon(h_n) = \theta(h_p) +C(h_p)(h_p-h_n) - \theta(h_n) \f]
  (specific storage must be zero in order to use this function!!)
*/
double richards::taylor_error(long ipp)
{
  double hn;
  double hp;
  
  //  actual hydraulic head
  hn = Tm->ip[ipp].av[0];

  //previous hydraulic head
  hp = Tm->ip[ipp].pv[0];
  
  return theta_val(hp) + c_value(hp)*(hp-hn) - theta_val(hn) ; 
  
}
  
  

  

