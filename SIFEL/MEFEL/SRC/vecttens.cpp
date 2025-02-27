#include "vecttens.h"
#include "matrix.h"
#include "vector.h"
#include "alias.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>



/**
  Function creates tensor components from %vector components.
   
  @param v - %vector components
  @param t - tensor components (output)
  @param ssst - strain/stress state
  @param ss - strain/stress identifier

  @return The function returns required tensor in the parameter t.
   
  Created by JK, 25.3.2002
*/
void vector_tensor (vector &v,matrix &t,strastrestate ssst,strastre ss)
{
  nullm (t);
  
  if (ss==strain){
    switch (ssst){
    case bar:{
      t(0,0)=v(0);
      break;
    }
    case plbeam:{
      break;
    }
    case planecontact:{
      t(0,0) = v(0);
      t(1,1) = v(1);
      break;
    }
    case planestress:{
      t(0,0)=v(0);      t(0,1)=v(2)/2.0;  t(0,2)=0.0;
      t(1,0)=v(2)/2.0;  t(1,1)=v(1);      t(1,2)=0.0;
      t(2,0)=0.0;       t(2,1)=0.0;       t(2,2)=v(3);
      break;
    }
    case planestrain:{
      t(0,0)=v(0);      t(0,1)=v(2)/2.0;  t(0,2)=0.0;
      t(1,0)=v(2)/2.0;  t(1,1)=v(1);      t(1,2)=0.0;
      t(2,0)=0.0;       t(2,1)=0.0;       t(2,2)=v(3);
      break;
    }
    case axisymm:{
      t(0,0)=v(0);      t(0,1)=v(3)/2.0;  t(0,2)=0.0;
      t(1,0)=v(3)/2.0;  t(1,1)=v(1);      t(1,2)=0.0;
      t(2,0)=0.0;       t(2,1)=0.0;       t(2,2)=v(2);
      break;
    }
    case spacestress:{
      t(0,0)=v(0);      t(0,1)=v(5)/2.0;  t(0,2)=v(4)/2.0;
      t(1,0)=v(5)/2.0;  t(1,1)=v(1);      t(1,2)=v(3)/2.0;
      t(2,0)=v(4)/2.0;  t(2,1)=v(3)/2.0;  t(2,2)=v(2);
      break;
    }
    default:
      print_err("unknown strain state %d is required", __FILE__, __LINE__, __func__, ssst);
    }
  }
  
  if (ss==stress){
    switch (ssst){
    case bar:{
      t(0,0)=v(0);
      break;
    }
    case plbeam:{
      break;
    }
    case planecontact:{
      t(0,0) = v(0);
      t(1,1) = v(1);
      break;
    }
    case planestress:{
      t(0,0)=v(0);  t(0,1)=v(2);  t(0,2)=0.0;
      t(1,0)=v(2);  t(1,1)=v(1);  t(1,2)=0.0;
      t(2,0)=0.0;   t(2,1)=0.0;   t(2,2)=v(3);
      break;
    }
    case planestrain:{
      t(0,0)=v(0);  t(0,1)=v(2);  t(0,2)=0.0;
      t(1,0)=v(2);  t(1,1)=v(1);  t(1,2)=0.0;
      t(2,0)=0.0;   t(2,1)=0.0;   t(2,2)=v(3);
      break;
    }
    case axisymm:{
      t(0,0)=v(0);  t(0,1)=v(3);  t(0,2)=0.0;
      t(1,0)=v(3);  t(1,1)=v(1);  t(1,2)=0.0;
      t(2,0)=0.0;   t(2,1)=0.0;   t(2,2)=v(2);
      break;
    }
    case spacestress:{
      t(0,0)=v(0);  t(0,1)=v(5);  t(0,2)=v(4);
      t(1,0)=v(5);  t(1,1)=v(1);  t(1,2)=v(3);
      t(2,0)=v(4);  t(2,1)=v(3);  t(2,2)=v(2);
      break;
    }
    default:
      print_err("unknown stress state %d is required", __FILE__, __LINE__, __func__, ssst);
    }
  }

}


/**
  Function creates %vector components from tensor components.
   
  @param v - %vector components (output)
  @param t - tensor components
  @param ssst - strain/stress state
  @param ss - strain/stress identifier
   
  @return The function returns required %vector in the parameter v.

  Created by JK, 25.3.2002
*/
void tensor_vector (vector &v,matrix &t,strastrestate ssst,strastre ss)
{
  nullv (v);
  
  if (ss==strain){
    switch (ssst){
    case bar:{
      v(0)=t(0,0);
      break;
    }
    case plbeam:{
      break;
    }
    case planecontact:{
      v(0) = t(0,0);
      v(1) = t(1,1);
      break;
    }
    case planestress:{
      v(0)=t(0,0);  v(1)=t(1,1);  v(2)=t(0,1)*2.0;  v(3)=t(2,2);
      break;
    }
    case planestrain:{
      v(0)=t(0,0);  v(1)=t(1,1);  v(2)=t(0,1)*2.0;  v(3)=t(2,2);
      break;
    }
    case axisymm:{
      v(0)=t(0,0);  v(1)=t(1,1);  v(2)=t(2,2);  v(3)=t(0,1)*2.0;
      break;
    }
    case spacestress:{
      v(0)=t(0,0);      v(1)=t(1,1);      v(2)=t(2,2);
      v(3)=t(1,2)*2.0;  v(4)=t(0,2)*2.0;  v(5)=t(0,1)*2.0;
      break;
    }
    default:
      print_err("unknown strain state is required", __FILE__, __LINE__, __func__);
    }
  }
  
  if (ss==stress){
    switch (ssst){
    case bar:{
      v(0)=t(0,0);
      break;
    }
    case plbeam:{
      break;
    }
    case planecontact:{
      v(0) = t(0,0);
      v(1) = t(1,1);
      break;
    }
    case planestress:{
      v(0)=t(0,0);  v(1)=t(1,1);  v(2)=t(0,1);  v(3)=t(2,2);
      break;
    }
    case planestrain:{
      v(0)=t(0,0);  v(1)=t(1,1);  v(2)=t(0,1);  v(3)=t(2,2);
      break;
    }
    case axisymm:{
      v(0)=t(0,0);  v(1)=t(1,1);  v(2)=t(2,2);  v(3)=t(0,1);
      break;
    }
    case spacestress:{
      v(0)=t(0,0);  v(1)=t(1,1);  v(2)=t(2,2);
      v(3)=t(1,2);  v(4)=t(0,2);  v(5)=t(0,1);
      break;
    }
    default:
      print_err("unknown stress state is required", __FILE__, __LINE__, __func__);
    }
  }
}

/**
  Function creates %vector components from tensor components.
  The %vector is returned in the full form.
   
  @param v - %vector components (output)
  @param t - tensor components
  @param ssst - strain/stress state
  @param ss - strain/stress identifier
   
  @return The function returns required %vector in the parameter v.

  Created by TKo, 19.1.2015
*/
void tensor_vector_full (vector &v, matrix &t, strastre ss)
{
  nullv (v);
  
  if (ss==strain){
    v(0)=t(0,0);      v(1)=t(1,1);      v(2)=t(2,2);
    v(3)=t(1,2)*2.0;  v(4)=t(0,2)*2.0;  v(5)=t(0,1)*2.0;
  }
  
  if (ss==stress){
    v(0)=t(0,0);  v(1)=t(1,1);  v(2)=t(2,2);
    v(3)=t(1,2);  v(4)=t(0,2);  v(5)=t(0,1);
  }
}


/**
  The function returns reduced form (ncompstr components) of stress/strain %vector from the given
  full-length %vector (6 components) of stress/strain. 

  @param fv - full length (6 components) %vector of stress/strain - output parameter
  @param rv - reduced %vector (ncompstr components) %vector of stress/strain - input parameter
  @param ssst - stress/strain state indicator

  @return The function returns reduced %vector in the argument rv

  Created by Tomas Koudelka 11.2.2014
*/
void give_red_vector(const vector &fv, vector &rv, strastrestate ssst)
{
  switch (ssst){
    case bar:{
      rv(0)=fv(0);
      break;
    }
    case plbeam:{
      break;
    }
    case planecontact:{
      rv(0) = fv(0);
      rv(1) = fv(1);
      break;
    }
    case planestress:{
      rv(0)=fv(0);  rv(1)=fv(1);  rv(2)=fv(5);  rv(3)=fv(2);
      break;
    }
    case planestrain:{
      rv(0)=fv(0);  rv(1)=fv(1);  rv(2)=fv(5);  rv(3)=fv(2);
      break;
    }
    case axisymm:{
      rv(0)=fv(0);  rv(1)=fv(1);  rv(2)=fv(2);  rv(3)=fv(5);
      break;
    }
    case spacestress:{
      rv(0)=fv(0);  rv(1)=fv(1);  rv(2)=fv(2);
      rv(3)=fv(3);  rv(4)=fv(4);  rv(5)=fv(5);
      break;
    }
    default:
      print_err("unknown strain state is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function returns full-length stress/strain %vector (6 components) from the reduced 
  stress/strain %vector (ncompstr components).
  . 

  @param fv - full length (6 components) %vector of stress/strain - output parameter
  @param rv - reduced %vector (ncompstr components) %vector of stress/strain - input parameter
  @param ssst - stress/strain state indicator

  @return The function returns full-length %vector in the argument rv

  Created by Tomas Koudelka 11.2.2014
*/
void give_full_vector(vector &fv, const vector &rv, strastrestate ssst)
{
  switch (ssst){
    case bar:{
      fv(0)=rv(0);
      break;
    }
    case plbeam:{
      break;
    }
    case planecontact:{
      fv(0) = rv(0);
      fv(1) = rv(1);
      break;
    }
    case planestress:{
      fv(0)=rv(0);  fv(1)=rv(1);  fv(5)=rv(2);  fv(2)=rv(3);
      break;
    }
    case planestrain:{
      fv(0)=rv(0);  fv(1)=rv(1);  fv(5)=rv(2);  fv(2)=rv(3);
      break;
    }
    case axisymm:{
      fv(0)=rv(0);  fv(1)=rv(1);  fv(2)=rv(2);  fv(5)=rv(3);
      break;
    }
    case spacestress:{
      fv(0)=rv(0);  fv(1)=rv(1);  fv(2)=rv(2);
      fv(3)=rv(3);  fv(4)=rv(4);  fv(5)=rv(5);
      break;
    }
    default:
      print_err("unknown strain state is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function returns indices of the shear components for the reduced %vector notation of stress/strain vectors/matrices(tensors).
  It is intended for the doubling shear components due to conversion from tensor to vector/matrix format.

  @param[in] ssst - stress/strain state indicator
  @param[out] id - resulting vector of shear component indices, %vector should be allocated at the input
                   to the maximum number of stress/strain shear components, i.e. 3,
                   it will be reallocated according to the true number of shear components.


  @return The function returns indices in the argument id

  Created by Tomas Koudelka 5.2022
*/
void give_shear_indices(strastrestate ssst, ivector &id)
{
  switch (ssst){
    case bar:{
      // no shear components
      reallocv(0, id);
      break;
    }
    case plbeam:{
      // no shear components
      reallocv(0, id);
      break;
    }
    case planecontact:{
      // one shear component tau
      reallocv(1, id);
      id(0) = 1;
      break;
    }
    case planestress:{
      // one shear component tau_xy
      reallocv(1, id);
      id(0) = 2;
      break;
    }
    case planestrain:{
      // one shear component tau_xy
      reallocv(1, id);
      id(0) = 2;
      break;
    }
    case axisymm:{
      // one shear component tau_xz
      reallocv(1, id);
      id(0) = 3;
      break;
    }
    case spacestress:{
      // three shear components tau_yz, tau_xz, tau_xy
      reallocv(3, id);
      id(0) = 3;
      id(1) = 4;
      id(2) = 5;
      break;
    }
    default:
      print_err("unknown strain state is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function returns indices of all shear components for the reduced %vector notation of stress/strain vectors/matrices(tensors).
  Index of shear components not used in the given format is set to -1. It is intended for the conversion from full vector/matrix 
  format to the reduced one in the case that the vector/matrix has non standard dimensions 
  (vectors/matrices with block strutures - stacked vectors/matrices - containing also other components than stresses or strains).

  @param[in] ssst - stress/strain state indicator
  @param[out] id - resulting vector of shear component indices, %vector should be allocated at the input
                   to the maximum number of stress/strain shear components, i.e. 3,
                   id[0] represents index of tau_yz, 
                   id[1] represents index of tau_xz, 
                   id[2] represents index of tau_xy.
  Resulting values of component indices are being assumed as follows
                   id[i] = -1 for i-th shear components not used in the given reduced format, 
                   id[i] >= 0 for used shear components 

  @return The function returns indices in the argument id

  Created by Tomas Koudelka 5.2022
*/
void give_all_shear_indices(strastrestate ssst, ivector &id)
{
  fillv(-1, id);
  switch (ssst){
    case bar:{
      // no shear components
      break;
    }
    case plbeam:{
      // no shear components
      break;
    }
    case planecontact:{
      // one shear component tau_xy
      id(2) = 1;
      break;
    }
    case planestress:{
      // one shear component tau_xy
      id(2) = 2;
      break;
    }
    case planestrain:{
      // one shear component tau_xy
      id(2) = 2;
      break;
    }
    case axisymm:{
      // one shear component tau_xz
      id(1) = 3;
      break;
    }
    case spacestress:{
      // three shear components tau_yz, tau_xz, tau_xy
      id(0) = 3;
      id(1) = 4;
      id(2) = 5;
      break;
    }
    default:
      print_err("unknown strain state is required", __FILE__, __LINE__, __func__);
  }  
}



/**
  The function returns indices of all normal components for the reduced %vector notation of stress/strain vectors/matrices(tensors).
  Index of normal components not used in the given format is set to -1. It is intended for the conversion from full vector/matrix 
  format to the reduced one in the case that the vector/matrix has nonstandard dimensions 
  (vectors/matrices with block strutures - stacked vectors/matrices - containing also other components than stresses or strains).

  @param[in] ssst - stress/strain state indicator
  @param[out] id - resulting vector of normal component indices, %vector should be allocated at the input
                   to the maximum number of stress/strain normal components, i.e. 3,
                   id[0] represents index of sig_x, 
                   id[1] represents index of sig_y, 
                   id[2] represents index of sig_z.
  Resulting values of component indices are being assumed as follows
                   id[i] = -1 for i-th normal components not used in the given reduced format, 
                   id[i] >= 0 for used normal components 

  @return The function returns indices in the argument id

  Created by Tomas Koudelka 5.2022
*/
void give_all_normal_indices(strastrestate ssst, ivector &id)
{
  fillv(-1, id);
  switch (ssst){
    case bar:{
      // one normal component sig_x
      fillv(-1, id);
      id(0) = 0;
      break;
    }
    case plbeam:{
      // no normal components
      fillv(-1, id);
      break;
    }
    case planecontact:{
      // one normal component sig_n
      fillv(-1, id);
      id(0) = 0;
      break;
    }
    case planestress:{
      // three component tau_xy
      id(0) = 0;
      id(1) = 1;
      id(2) = 3;
      break;
    }
    case planestrain:{
      // one shear component tau_xy
      id(0) = 0;
      id(1) = 1;
      id(2) = 3;
      break;
    }
    case axisymm:{
      // one shear component tau_xz
      id(0) = 0;
      id(1) = 1;
      id(2) = 2;
      break;
    }
    case spacestress:{
      // three shear components tau_yz, tau_xz, tau_xy
      id(0) = 0;
      id(1) = 1;
      id(2) = 2;
      break;
    }
    default:
      print_err("unknown strain state is required", __FILE__, __LINE__, __func__);
  }  
}



/**
  The function returns reduced form (ncompstr components) of stress/strain %vector from the given
  full-length %vector (6 components) of stress/strain. 

  @param fv - array of full length %vector (6 components) - output parameter
  @param rv - array of reduced %vector (ncompstr components) - input parameter
  @param ssst - stress/strain state indicator

  @return The function returns reduced %vector components in the argument rv

  Created by Tomas Koudelka 11.2.2014
*/
void give_red_vector(double *fv, double *rv, strastrestate ssst)
{
  switch (ssst){
    case bar:{
      rv[0]=fv[0];
      break;
    }
    case plbeam:{
      break;
    }
    case planecontact:{
      rv[0] = fv[0];
      rv[1] = fv[1];
      break;
    }
    case planestress:{
      rv[0]=fv[0];  rv[1]=fv[1];  rv[2]=fv[5];  rv[3]=fv[2];
      break;
    }
    case planestrain:{
      rv[0]=fv[0];  rv[1]=fv[1];  rv[2]=fv[5];  rv[3]=fv[2];
      break;
    }
    case axisymm:{
      rv[0]=fv[0];  rv[1]=fv[1];  rv[2]=fv[2];  rv[3]=fv[5];
      break;
    }
    case spacestress:{
      rv[0]=fv[0];  rv[1]=fv[1];  rv[2]=fv[2];
      rv[3]=fv[3];  rv[4]=fv[4];  rv[5]=fv[5];
      break;
    }
    default:
      print_err("unknown strain state is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function returns full-length stress/strain %vector (6 components) from the reduced 
  stress/strain %vector (ncompstr components).
  . 

  @param fv - array of full length %vector (6 components) - output parameter
  @param rv - array of reduced %vector (ncompstr components) - input parameter
  @param ssst - stress/strain state indicator

  @return The function returns full-length %vector components in the argument rv

  Created by Tomas Koudelka 11.2.2014
*/
void give_full_vector(double *fv, double *rv, strastrestate ssst)
{
  switch (ssst){
    case bar:{
      fv[0]=rv[0];
      break;
    }
    case plbeam:{
      break;
    }
    case planecontact:{
      fv[0] = rv[0];
      fv[1] = rv[1];
      break;
    }
    case planestress:{
      fv[0]=rv[0];  fv[1]=rv[1];  fv[5]=rv[2];  fv[2]=rv[3];
      break;
    }
    case planestrain:{
      fv[0]=rv[0];  fv[1]=rv[1];  fv[5]=rv[2];  fv[2]=rv[3];
      break;
    }
    case axisymm:{
      fv[0]=rv[0];  fv[1]=rv[1];  fv[2]=rv[2];  fv[5]=rv[3];
      break;
    }
    case spacestress:{
      fv[0]=rv[0];  fv[1]=rv[1];  fv[2]=rv[2];
      fv[3]=rv[3];  fv[4]=rv[4];  fv[5]=rv[5];
      break;
    }
    default:
      print_err("unknown strain state is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function returns full-length stress/strain %vector (6 components) from the reduced 
  stress/strain %vector (ncompstr components).
  . 

  @param[out] fv - full length %vector (6 components)
  @param[in] rv - array of reduced %vector (ncompstr components)
  @param ssst - stress/strain state indicator

  @return The function returns full-length %vector components in the argument rv

  Created by Tomas Koudelka 11.2.2014
*/
void give_full_vector(vector &fv, double *rv, strastrestate ssst)
{
  switch (ssst){
    case bar:{
      fv[0]=rv[0];
      break;
    }
    case plbeam:{
      break;
    }
    case planecontact:{
      fv[0] = rv[0];
      fv[1] = rv[1];
      break;
    }
    case planestress:{
      fv[0]=rv[0];  fv[1]=rv[1];  fv[5]=rv[2];  fv[2]=rv[3];
      break;
    }
    case planestrain:{
      fv[0]=rv[0];  fv[1]=rv[1];  fv[5]=rv[2];  fv[2]=rv[3];
      break;
    }
    case axisymm:{
      fv[0]=rv[0];  fv[1]=rv[1];  fv[2]=rv[2];  fv[5]=rv[3];
      break;
    }
    case spacestress:{
      fv[0]=rv[0];  fv[1]=rv[1];  fv[2]=rv[2];
      fv[3]=rv[3];  fv[4]=rv[4];  fv[5]=rv[5];
      break;
    }
    default:
      print_err("unknown strain state is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function creates auxiliary %matrix m
  (useful for thermal strains, pore pressure, etc.)

  @param ssst - strain/stress state
  @param m - auxiliary %matrix (output)
   
  @returns The function returns auxiliary %matrix in the parameter m.

  Created by JK, 24.10.2004
*/
void tensor_vector_matrix (strastrestate ssst,matrix &m)
{
  nullm (m);
  
  switch (ssst){
  case bar:{
    m(0,0)=1.0;
    break;
  }
  case planestress:{
    m(0,0)=1.0;
    m(1,0)=1.0;
    m(2,0)=0.0;
    break;
  }
  case axisymm:{
    m(0,0)=1.0;
    m(1,0)=1.0;
    m(2,0)=1.0;
    m(3,0)=0.0;
    break;
  }
  case spacestress:{
    m(0,0)=1.0;
    m(1,0)=1.0;
    m(2,0)=1.0;
    m(3,0)=0.0;
    m(4,0)=0.0;
    m(5,0)=0.0;
    break;
  }
  default:
    print_err("unknown stress state is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function assembles reduced material %matrix corresponding to the fourth order tensor.
   
  @param[out] m    - %matrix form 
  @param[in]  t    - tensor form (%matrix 6x6)
  @param[in]  ssst - strain/stress state indicator

  @return The function returns required %matrix in the parameter m.

  Created by JK, 18.3.2005
*/
void tensor4_matrix (matrix &m, const matrix &t, strastrestate ssst)
{
  nullm (m);
  
  switch (ssst){
  case bar:{
    m(0,0)=t(0,0);
    break;
  }
  case plbeam:{
    break;
  }
  case planestress:{
    m(0,0)=t(0,0);  m(0,1)=t(0,1);  m(0,2)=t(0,5);
    m(1,0)=t(1,0);  m(1,1)=t(1,1);  m(1,2)=t(1,5);
    m(2,0)=t(5,0);  m(2,1)=t(5,1);  m(2,2)=t(5,5); 
    break;
  }
  case planestrain:{
    m(0,0)=t(0,0);  m(0,1)=t(0,1);  m(0,2)=t(0,5);  m(0,3)=t(0,2);
    m(1,0)=t(1,0);  m(1,1)=t(1,1);  m(1,2)=t(1,5);  m(1,3)=t(1,2);
    m(2,0)=t(5,0);  m(2,1)=t(5,1);  m(2,2)=t(5,5);  m(2,3)=t(5,2);
    m(3,0)=t(2,0);  m(3,1)=t(2,1);  m(3,2)=t(2,5);  m(3,3)=t(2,2);
    break;
  }
  case axisymm:{
    m(0,0)=t(0,0);  m(0,1)=t(0,1);  m(0,2)=t(0,2);  m(0,3)=t(0,5);
    m(1,0)=t(1,0);  m(1,1)=t(1,1);  m(1,2)=t(1,2);  m(1,3)=t(1,5);
    m(2,0)=t(2,0);  m(2,1)=t(2,1);  m(2,2)=t(2,2);  m(2,3)=t(2,5);
    m(3,0)=t(5,0);  m(3,1)=t(5,1);  m(3,2)=t(5,2);  m(3,3)=t(5,5);
    break;
  }
  case spacestress:{
    copym(t,m);
    break;
  }
  default:
    print_err("unknown strain/stress state is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function assembles reduced element stiffness %matrix corresponding to the full fourth order stiffness 
  tensor. In the case of plane-stress problem, the normal components x and y are modified in order to 
  take into account $f \sigma_z = 0 $f.
   
  @param[out] m    - reduced %matrix form
  @param[in]  t    - tensor form (%matrix 6x6)
  @param[in]  ssst - strain/stress state

  @return The function returns required %matrix in the parameter m.

  Created by TKo, 8.3.2016
*/
void tensor4_ematrix (matrix &m, const matrix &t, strastrestate ssst)
{
  nullm (m);
  
  switch (ssst){
  case bar:{
    m(0,0)=t(0,0);
    break;
  }
  case plbeam:{
    break;
  }
  case planestress:{
    // it is valid for isotropic or orthotropic material
    m(0,0)=t(0,0)-t(0,2)*t(2,0)/t(2,2);  m(0,1)=t(0,1)-t(0,2)*t(2,1)/t(2,2);  m(0,2)=t(0,5);
    m(1,0)=t(1,0)-t(1,2)*t(2,0)/t(2,2);  m(1,1)=t(1,1)-t(1,2)*t(2,1)/t(2,2);  m(1,2)=t(1,5);
    if ((fabs(t(5,0)/t(5,5)) < 1.0e-7) || (fabs(t(5,1)/t(5,5)) < 1.0e-7)){ // isotropic or orthotropic material
      m(2,0)=t(5,0);  m(2,1)=t(5,1);  m(2,2)=t(5,5);
    }
    else{
      // anistropic material
      print_err("conversion for the anisotropic material for plane-stress problem has not yet been implemened\n"
                "t_61=%le, t_62=%le\n", __FILE__, __LINE__, __func__, t(5,0), t(5,1));
      abort();
    }
    break;
  }
  case planestrain:{
    m(0,0)=t(0,0);  m(0,1)=t(0,1);  m(0,2)=t(0,5);
    m(1,0)=t(1,0);  m(1,1)=t(1,1);  m(1,2)=t(1,5);
    m(2,0)=t(5,0);  m(2,1)=t(5,1);  m(2,2)=t(5,5);
    break;
  }
  case axisymm:{
    m(0,0)=t(0,0);  m(0,1)=t(0,1);  m(0,2)=t(0,2);  m(0,3)=t(0,5);
    m(1,0)=t(1,0);  m(1,1)=t(1,1);  m(1,2)=t(1,2);  m(1,3)=t(1,5);
    m(2,0)=t(2,0);  m(2,1)=t(2,1);  m(2,2)=t(2,2);  m(2,3)=t(2,5);
    m(3,0)=t(5,0);  m(3,1)=t(5,1);  m(3,2)=t(5,2);  m(3,3)=t(5,5);
    break;
  }
  case spacestress:{
    copym(t,m);
    break;
  }
  default:
    print_err("unknown strain/stress state is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function assembles the fourth order tensor corresponding to the given %matrix 
  in reduced form.
   
  @param m[in] - %matrix in the reduced form 
  @param t[out] - tensor form (%matrix 6x6)
  @param ssst[in] - strain/stress state

  @return The function returns required %matrix in the parameter m.

  Created by Tomas Koudelka, 11.12.2014
*/
void matrix_tensor4 (const matrix &m, matrix &t, strastrestate ssst)
{
  nullm (t);

  switch (ssst){
  case bar:{
    t(0,0) = m(0,0);
    break;
  }
  case plbeam:{
    break;
  }
  case planestress:{
    t(0,0)=m(0,0);  t(0,1)=m(0,1);  t(0,5)=m(0,2);
    t(1,0)=m(1,0);  t(1,1)=m(1,1);  t(1,5)=m(1,2);
    t(5,0)=m(2,0);  t(5,1)=m(2,1);  t(5,5)=m(2,2);
    break;
  }
  case planestrain:{
    t(0,0)=m(0,0);  t(0,1)=m(0,1);  t(0,2)=m(0,2);  t(0,5)=m(0,3);
    t(1,0)=m(1,0);  t(1,1)=m(1,1);  t(1,2)=m(1,2);  t(1,5)=m(1,3);
    t(2,0)=m(2,0);  t(2,1)=m(2,1);  t(2,2)=m(2,2);  t(2,5)=m(2,3);
    t(5,0)=m(3,0);  t(5,1)=m(3,1);  t(5,2)=m(3,2);  t(5,5)=m(3,3);
    break;
  }
  case axisymm:{
    t(0,0)=m(0,0);  t(0,1)=m(0,1);  t(0,2)=m(0,2);  t(0,5)=m(0,3);
    t(1,0)=m(1,0);  t(1,1)=m(1,1);  t(1,2)=m(1,2);  t(1,5)=m(1,3);
    t(2,0)=m(2,0);  t(2,1)=m(2,1);  t(2,2)=m(2,2);  t(2,5)=m(2,3);
    t(5,0)=m(3,0);  t(5,1)=m(3,1);  t(5,2)=m(3,2);  t(5,5)=m(3,3);
    break;
  }
  case spacestress:{
    copym(m,t);
    break;
  }
  default:
    print_err("unknown strain/stress state is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function assembles %matrix corresponding to the fourth order tensor.
   
  @param m - %matrix form (output)
  @param t - tensor form
  @param ssst - strain/stress state
   
  @return The function returns required %matrix in the parameter m.

  Created by 18.3.2005, JK
*/
void red_rows_matrix (matrix &m, matrix &t, strastrestate ssst)
{
  long i,n;
  nullm (m);
  
  n=t.n;

  switch (ssst){
  case bar:{
    for (i=0;i<n;i++){
      m(0,i)=t(0,i);
    }
    break;
  }
  case plbeam:{
    break;
  }
  case planestress:{
    for (i=0;i<n;i++){
      m(0,i)=t(0,i);
      m(1,i)=t(1,i);
      m(2,i)=t(5,i);
    }
    break;
  }
  case planestrain:{
    for (i=0;i<n;i++){
      m(0,i)=t(0,i);
      m(1,i)=t(1,i);
      m(2,i)=t(5,i);
      m(3,i)=t(2,i);
    }
    break;
  }
  case axisymm:{
    for (i=0;i<n;i++){
      m(0,i)=t(0,i);
      m(1,i)=t(1,i);
      m(2,i)=t(2,i);
      m(3,i)=t(5,i);
    }
    break;
  }
  case spacestress:{
    copym(t,m);
    break;
  }
  default:
    print_err("unknown strain/stress state is required", __FILE__, __LINE__, __func__);
  }

}


/**
  Function assembles %matrix corresponding to the fourth order tensor.
   
  @param m - %matrix form (output)
  @param t - tensor form
  @param ssst - strain/stress state
   
  @return The function returns required %matrix in the parameter m.

  Created by 18.3.2005, JK
*/
void red_cols_matrix (matrix &m, matrix &t, strastrestate ssst)
{
  long i,n;
  nullm (m);
  
  n=t.m;

  switch (ssst){
  case bar:{
    for (i=0;i<n;i++){
      m(i,0)=t(i,0);
    }
    break;
  }
  case plbeam:{
    break;
  }
  case planestress:{
    for (i=0;i<n;i++){
      m(i,0)=t(i,0);
      m(i,1)=t(i,1);
      m(i,2)=t(i,5);
    }
    break;
  }
  case planestrain:{
    for (i=0;i<n;i++){
      m(i,0)=t(i,0);
      m(i,1)=t(i,1);
      m(i,2)=t(i,5);
      m(i,3)=t(i,2);
    }
    break;
  }
  case axisymm:{
    for (i=0;i<n;i++){
      m(i,0)=t(i,0);
      m(i,1)=t(i,1);
      m(i,2)=t(i,2);
      m(i,3)=t(i,5);
    }
    break;
  }
  case spacestress:{
    copym(t,m);
    break;
  }
  default:
    print_err("unknown strain/stress state is required", __FILE__, __LINE__, __func__);
  }

}


/**
  Function tries to guess stress/strain state from the number of
  components of some %vector (ncompstr at nodes). The result need not
  to be correct. It is only guess and it is useless for beam
  stress/strain state.

  @param ncomp - number of %vector components

  @retval bar for 1 component
  @retval plcontact for ncomp=2
  @retval planestrain for ncomp=3
  @retval planestrain for ncomp=4
  @retval spacestress for ncomp=6
  @retval -1 for any other ncomp value

  Created by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz, 4.10.2007 
*/
strastrestate guess_ssst(long ncomp)
{
  switch (ncomp)
  {
    case 1:
      return bar;
    case 2:
      return planecontact;
    case 3:
    case 4:
      return planestrain;
    case 6:
      return spacestress;
    default:
      return strastrestate(-1);
  }
  return strastrestate(-1);
}


/**
  The function returns number of stress/strain array components depending on the stress/strain state indicator.

  @param[in] ssst - required stress/strain state indicator

  @return The function returns number of stress/strain array components for the given stress strain indicator.
*/
long give_ncompstr(strastrestate ssst)
{
  switch (ssst){
    case bar:
      return 1;
    case plbeam:
      return 3;
    case spacebeam:
      return 6;
    case planestress:
    case planestrain:      
      return 4;
    case planecontact:
      return 2;
    case platek:
    case plates:
      return 3;
    case axisymm:
      return 4;
    case shell:
      return 6;
    case spacestress:
      return 6;
    default:
      print_err("unknown strain state is required", __FILE__, __LINE__, __func__);
      abort();
  }
  return 0;
}



/**
  The function transforms the second order stress/strain tensor l given in the local coordinate system to 
  tensor g in global coordinate system. Both tensors l and g are defined in the %vector engineering 
  notation but transformation %matrix T(3,3) must be given in the form x_g = T x_l where x_l(3) and x_g(3)
  are local and global vectors respectively.
  
  @param g - 2-nd order stress/strain tensor in the %vector engineering notation defined in 
             the global coordinate system (output)
  @param l - 2-nd order stress/strain tensor in the %vector engineering notation defined in 
             the local coordinate system (input)
  @param tmat - transformation %matrix T(3,3)
  @param ssst - stress/strain state indicator - defines the number of g and l components
  @param stra - stress/strain indicator - defines whether g and l represents stresses or strains
*/  
void lg_engvectortransf (vector &g, const vector &l, const matrix &tmat, strastrestate ssst, strastre stra)
{
  matrix t(ASTCKMAT(6,6));
  vector lf(ASTCKVEC(6)), gf(ASTCKVEC(6));

  give_full_vector(lf, l, ssst);

  switch (stra)
  {
    case strain:
    {
      // Transformation matrix for strains in vector notation i.e. eps_g = T_eps eps_l

      // the first row of transf. matrix
      t(0,0) = tmat(0,0)*tmat(0,0);
      t(0,1) = tmat(0,1)*tmat(0,1);
      t(0,2) = tmat(0,2)*tmat(0,2);
      t(0,3) = tmat(0,1)*tmat(0,2);
      t(0,4) = tmat(0,0)*tmat(0,2);
      t(0,5) = tmat(0,0)*tmat(0,1);
  
      // the second row of transf. matrix
      t(1,0) = tmat(1,0)*tmat(1,0);
      t(1,1) = tmat(1,1)*tmat(1,1);
      t(1,2) = tmat(1,2)*tmat(1,2);
      t(1,3) = tmat(1,1)*tmat(1,2);
      t(1,4) = tmat(1,0)*tmat(1,2);
      t(1,5) = tmat(1,0)*tmat(1,1);
      
      // the third row of transf. matrix
      t(2,0) = tmat(2,0)*tmat(2,0);
      t(2,1) = tmat(2,1)*tmat(2,1);
      t(2,2) = tmat(2,2)*tmat(2,2);
      t(2,3) = tmat(2,1)*tmat(2,2);
      t(2,4) = tmat(2,0)*tmat(2,2);
      t(2,5) = tmat(2,0)*tmat(2,1);
      
      // the fourth row of transf. matrix
      t(3,0) = 2.0*tmat(2,0)*tmat(1,0);
      t(3,1) = 2.0*tmat(2,1)*tmat(1,1);
      t(3,2) = 2.0*tmat(2,2)*tmat(1,2);
      t(3,3) = tmat(2,1)*tmat(1,2)+tmat(2,2)*tmat(1,1);
      t(3,4) = tmat(2,0)*tmat(1,2)+tmat(2,2)*tmat(1,0);
      t(3,5) = tmat(2,0)*tmat(1,1)+tmat(2,1)*tmat(1,0);
      
      // the fifth row of transf. matrix
      t(4,0) = 2.0*tmat(2,0)*tmat(0,0);
      t(4,1) = 2.0*tmat(2,1)*tmat(0,1);
      t(4,2) = 2.0*tmat(2,2)*tmat(0,2);
      t(4,3) = tmat(2,1)*tmat(0,2)+tmat(2,2)*tmat(0,1);
      t(4,4) = tmat(2,0)*tmat(0,2)+tmat(2,2)*tmat(0,0);
      t(4,5) = tmat(2,0)*tmat(0,1)+tmat(2,1)*tmat(0,0);
  
      // the sixth row of transf. matrix
      t(5,0) = 2.0*tmat(1,0)*tmat(0,0);
      t(5,1) = 2.0*tmat(1,1)*tmat(0,1);
      t(5,2) = 2.0*tmat(1,2)*tmat(0,2);
      t(5,3) = tmat(1,1)*tmat(0,2)+tmat(1,2)*tmat(0,1);
      t(5,4) = tmat(1,0)*tmat(0,2)+tmat(1,2)*tmat(0,0);
      t(5,5) = tmat(1,0)*tmat(0,1)+tmat(1,1)*tmat(0,0);

      break;
    }
    case stress:
    {
      // Transformation matrix for stresses in vector notation i.e. sig_g = T_sig sig_l

      // the first row of transf. matrix
      t(0,0) = tmat(0,0)*tmat(0,0);
      t(0,1) = tmat(0,1)*tmat(0,1);
      t(0,2) = tmat(0,2)*tmat(0,2);
      t(0,3) = 2.0*tmat(0,1)*tmat(0,2);
      t(0,4) = 2.0*tmat(0,0)*tmat(0,2);
      t(0,5) = 2.0*tmat(0,0)*tmat(0,1);
  
      // the second row of transf. matrix
      t(1,0) = tmat(1,0)*tmat(1,0);
      t(1,1) = tmat(1,1)*tmat(1,1);
      t(1,2) = tmat(1,2)*tmat(1,2);
      t(1,3) = 2.0*tmat(1,1)*tmat(1,2);
      t(1,4) = 2.0*tmat(1,0)*tmat(1,2);
      t(1,5) = 2.0*tmat(1,0)*tmat(1,1);
      
      // the third row of transf. matrix
      t(2,0) = tmat(2,0)*tmat(2,0);
      t(2,1) = tmat(2,1)*tmat(2,1);
      t(2,2) = tmat(2,2)*tmat(2,2);
      t(2,3) = 2.0*tmat(2,1)*tmat(2,2);
      t(2,4) = 2.0*tmat(2,0)*tmat(2,2);
      t(2,5) = 2.0*tmat(2,0)*tmat(2,1);
      
      // the fourth row of transf. matrix
      t(3,0) = tmat(2,0)*tmat(1,0);
      t(3,1) = tmat(2,1)*tmat(1,1);
      t(3,2) = tmat(2,2)*tmat(1,2);
      t(3,3) = tmat(2,1)*tmat(1,2)+tmat(2,2)*tmat(1,1);
      t(3,4) = tmat(2,0)*tmat(1,2)+tmat(2,2)*tmat(1,0);
      t(3,5) = tmat(2,0)*tmat(1,1)+tmat(2,1)*tmat(1,0);
      
      // the fifth row of transf. matrix
      t(4,0) = tmat(2,0)*tmat(0,0);
      t(4,1) = tmat(2,1)*tmat(0,1);
      t(4,2) = tmat(2,2)*tmat(0,2);
      t(4,3) = tmat(2,1)*tmat(0,2)+tmat(2,2)*tmat(0,1);
      t(4,4) = tmat(2,0)*tmat(0,2)+tmat(2,2)*tmat(0,0);
      t(4,5) = tmat(2,0)*tmat(0,1)+tmat(2,1)*tmat(0,0);
  
      // the sixth row of transf. matrix
      t(5,0) = tmat(1,0)*tmat(0,0);
      t(5,1) = tmat(1,1)*tmat(0,1);
      t(5,2) = tmat(1,2)*tmat(0,2);
      t(5,3) = tmat(1,1)*tmat(0,2)+tmat(1,2)*tmat(0,1);
      t(5,4) = tmat(1,0)*tmat(0,2)+tmat(1,2)*tmat(0,0);
      t(5,5) = tmat(1,0)*tmat(0,1)+tmat(1,1)*tmat(0,0);
     
      break;
    }
    default:
      print_err("unknown type of tensor quantity is required", __FILE__, __LINE__, __func__);
  }
  mxv(t,lf,gf);
  give_red_vector(gf, g, ssst);
}



/**
  The function transforms the second order stress/strain tensor g given in the global coordinate system to 
  tensor l in local coordinate system. Both tensors l and g are defined in the %vector engineering 
  notation but transformation %matrix T(3,3) must be given in the form x_g = T x_l where x_l(3) and x_g(3)
  are local and global vectors respectively.
  
  @param g - 2-nd order stress/strain tensor in the %vector engineering notation defined in 
             the global coordinate system (input)
  @param l - 2-nd order stress/strain tensor in the %vector engineering notation defined in 
             the local coordinate system (output)
  @param tmat - transformation %matrix T(3,3)
  @param ssst - stress/strain state indicator - defines the number of g and l components
  @param stra - stress/strain indicator - defines whether g and l represents stresses or strains
*/  
void gl_engvectortransf (const vector &g, vector &l, const matrix &tmat, strastrestate ssst, strastre stra)
{
  matrix t(6,6);
  vector lf(6), gf(6);

  give_full_vector(gf, g, ssst);

  switch (stra)
  {
    case strain:
    {
      // Transformation matrix for strains in vector notation i.e. eps_g = T_eps eps_l

      // the first row of transf. matrix
      t(0,0) = tmat(0,0)*tmat(0,0);
      t(0,1) = tmat(1,0)*tmat(1,0);
      t(0,2) = tmat(2,0)*tmat(2,0);
      t(0,3) = tmat(1,0)*tmat(2,0);
      t(0,4) = tmat(0,0)*tmat(2,0);
      t(0,5) = tmat(0,0)*tmat(1,0);
  
      // the second row of transf. matrix
      t(1,0) = tmat(0,1)*tmat(0,1);
      t(1,1) = tmat(1,1)*tmat(1,1);
      t(1,2) = tmat(2,1)*tmat(2,1);
      t(1,3) = tmat(1,1)*tmat(2,1);
      t(1,4) = tmat(0,1)*tmat(2,1);
      t(1,5) = tmat(0,1)*tmat(1,1);
      
      // the third row of transf. matrix
      t(2,0) = tmat(0,2)*tmat(0,2);
      t(2,1) = tmat(1,2)*tmat(1,2);
      t(2,2) = tmat(2,2)*tmat(2,2);
      t(2,3) = tmat(1,2)*tmat(2,2);
      t(2,4) = tmat(0,2)*tmat(2,2);
      t(2,5) = tmat(0,2)*tmat(1,2);
      
      // the fourth row of transf. matrix
      t(3,0) = 2.0*tmat(0,2)*tmat(0,1);
      t(3,1) = 2.0*tmat(1,2)*tmat(1,1);
      t(3,2) = 2.0*tmat(2,2)*tmat(2,1);
      t(3,3) = tmat(1,2)*tmat(2,1)+tmat(2,2)*tmat(1,1);
      t(3,4) = tmat(0,2)*tmat(2,1)+tmat(2,2)*tmat(0,1);
      t(3,5) = tmat(0,2)*tmat(1,1)+tmat(1,2)*tmat(0,1);
      
      // the fifth row of transf. matrix
      t(4,0) = 2.0*tmat(0,2)*tmat(0,0);
      t(4,1) = 2.0*tmat(1,2)*tmat(1,0);
      t(4,2) = 2.0*tmat(2,2)*tmat(2,0);
      t(4,3) = tmat(1,2)*tmat(2,0)+tmat(2,2)*tmat(1,0);
      t(4,4) = tmat(0,2)*tmat(2,0)+tmat(2,2)*tmat(0,0);
      t(4,5) = tmat(0,2)*tmat(1,0)+tmat(1,2)*tmat(0,0);
  
      // the sixth row of transf. matrix
      t(5,0) = 2.0*tmat(0,1)*tmat(0,0);
      t(5,1) = 2.0*tmat(1,1)*tmat(1,0);
      t(5,2) = 2.0*tmat(2,1)*tmat(2,0);
      t(5,3) = tmat(1,1)*tmat(2,0)+tmat(2,1)*tmat(1,0);
      t(5,4) = tmat(0,1)*tmat(2,0)+tmat(2,1)*tmat(0,0);
      t(5,5) = tmat(0,1)*tmat(1,0)+tmat(1,1)*tmat(0,0);

      break;
    }
    case stress:
    {
      // Transformation matrix for strains in vector notation i.e. eps_g = T_eps eps_l

      // the first row of transf. matrix
      t(0,0) = tmat(0,0)*tmat(0,0);
      t(0,1) = tmat(1,0)*tmat(1,0);
      t(0,2) = tmat(2,0)*tmat(2,0);
      t(0,3) = 2.0*tmat(1,0)*tmat(2,0);
      t(0,4) = 2.0*tmat(0,0)*tmat(2,0);
      t(0,5) = 2.0*tmat(0,0)*tmat(1,0);
  
      // the second row of transf. matrix
      t(1,0) = tmat(0,1)*tmat(0,1);
      t(1,1) = tmat(1,1)*tmat(1,1);
      t(1,2) = tmat(2,1)*tmat(2,1);
      t(1,3) = 2.0*tmat(1,1)*tmat(2,1);
      t(1,4) = 2.0*tmat(0,1)*tmat(2,1);
      t(1,5) = 2.0*tmat(0,1)*tmat(1,1);
    
      // the third row of transf. matrix
      t(2,0) = tmat(0,2)*tmat(0,2);
      t(2,1) = tmat(1,2)*tmat(1,2);
      t(2,2) = tmat(2,2)*tmat(2,2);
      t(2,3) = 2.0*tmat(1,2)*tmat(2,2);
      t(2,4) = 2.0*tmat(0,2)*tmat(2,2);
      t(2,5) = 2.0*tmat(0,2)*tmat(1,2);
      
      // the fourth row of transf. matrix
      t(3,0) = tmat(0,2)*tmat(0,1);
      t(3,1) = tmat(1,2)*tmat(1,1);
      t(3,2) = tmat(2,2)*tmat(2,1);
      t(3,3) = tmat(1,2)*tmat(2,1)+tmat(2,2)*tmat(1,1);
      t(3,4) = tmat(0,2)*tmat(2,1)+tmat(2,2)*tmat(0,1);
      t(3,5) = tmat(0,2)*tmat(1,1)+tmat(1,2)*tmat(0,1);
      
      // the fifth row of transf. matrix
      t(4,0) = tmat(0,2)*tmat(0,0);
      t(4,1) = tmat(1,2)*tmat(1,0);
      t(4,2) = tmat(2,2)*tmat(2,0);
      t(4,3) = tmat(1,2)*tmat(2,0)+tmat(2,2)*tmat(1,0);
      t(4,4) = tmat(0,2)*tmat(2,0)+tmat(2,2)*tmat(0,0);
      t(4,5) = tmat(0,2)*tmat(1,0)+tmat(1,2)*tmat(0,0);
  
      // the sixth row of transf. matrix
      t(5,0) = tmat(0,1)*tmat(0,0);
      t(5,1) = tmat(1,1)*tmat(1,0);
      t(5,2) = tmat(2,1)*tmat(2,0);
      t(5,3) = tmat(1,1)*tmat(2,0)+tmat(2,1)*tmat(1,0);
      t(5,4) = tmat(0,1)*tmat(2,0)+tmat(2,1)*tmat(0,0);
      t(5,5) = tmat(0,1)*tmat(1,0)+tmat(1,1)*tmat(0,0);
     
      break;
    }
    default:
      print_err("unknown type of tensor quantity is required", __FILE__, __LINE__, __func__);
  }
  mxv(t,gf,lf);
  give_red_vector(lf, l, ssst);
}



/**
  The function transforms the second order stress/strain tensor g given in the global coordinate system to 
  local coordinate system and calculates only required i-th local engineering component. Tensor g 
  is defined in the %vector engineering notation but transformation %matrix T(3,3) must be given in 
  the form x_g = T x_l where x_l(3) and x_g(3) are local and global vectors respectively.
  
  @param[in]  g - 2-nd order stress/strain tensor in the %vector engineering notation defined in 
                  the global coordinate system
  @param[out] l - i-th component 2-nd order stress/strain tensor in the %vector engineering notation defined in 
                  the local coordinate system (output)
  @param[in]  i - index of required stress/strain component in the egineering notation
  @param[in]  tmat - transformation %matrix T(3,3)
  @param[in]  ssst - stress/strain state indicator - defines the number of g and l components
  @param[in]  stra - stress/strain indicator - defines whether g and l represents stresses or strains
*/  
void gl_comp_engvectortransf (const vector &g, double &l, long i, const matrix &tmat, strastrestate ssst, strastre stra)
{
  matrix t(ASTCKMAT(6,6));
  vector gf(ASTCKVEC(6));

  give_full_vector(gf, g, ssst);

  switch (stra)
  {
    case strain:
    {
      // Transformation matrix for strains in vector notation i.e. eps_g = T_eps eps_l

      // the first row of transf. matrix
      t(0,0) = tmat(0,0)*tmat(0,0);
      t(0,1) = tmat(1,0)*tmat(1,0);
      t(0,2) = tmat(2,0)*tmat(2,0);
      t(0,3) = tmat(1,0)*tmat(2,0);
      t(0,4) = tmat(0,0)*tmat(2,0);
      t(0,5) = tmat(0,0)*tmat(1,0);
  
      // the second row of transf. matrix
      t(1,0) = tmat(0,1)*tmat(0,1);
      t(1,1) = tmat(1,1)*tmat(1,1);
      t(1,2) = tmat(2,1)*tmat(2,1);
      t(1,3) = tmat(1,1)*tmat(2,1);
      t(1,4) = tmat(0,1)*tmat(2,1);
      t(1,5) = tmat(0,1)*tmat(1,1);
      
      // the third row of transf. matrix
      t(2,0) = tmat(0,2)*tmat(0,2);
      t(2,1) = tmat(1,2)*tmat(1,2);
      t(2,2) = tmat(2,2)*tmat(2,2);
      t(2,3) = tmat(1,2)*tmat(2,2);
      t(2,4) = tmat(0,2)*tmat(2,2);
      t(2,5) = tmat(0,2)*tmat(1,2);
      
      // the fourth row of transf. matrix
      t(3,0) = 2.0*tmat(0,2)*tmat(0,1);
      t(3,1) = 2.0*tmat(1,2)*tmat(1,1);
      t(3,2) = 2.0*tmat(2,2)*tmat(2,1);
      t(3,3) = tmat(1,2)*tmat(2,1)+tmat(2,2)*tmat(1,1);
      t(3,4) = tmat(0,2)*tmat(2,1)+tmat(2,2)*tmat(0,1);
      t(3,5) = tmat(0,2)*tmat(1,1)+tmat(1,2)*tmat(0,1);
      
      // the fifth row of transf. matrix
      t(4,0) = 2.0*tmat(0,2)*tmat(0,0);
      t(4,1) = 2.0*tmat(1,2)*tmat(1,0);
      t(4,2) = 2.0*tmat(2,2)*tmat(2,0);
      t(4,3) = tmat(1,2)*tmat(2,0)+tmat(2,2)*tmat(1,0);
      t(4,4) = tmat(0,2)*tmat(2,0)+tmat(2,2)*tmat(0,0);
      t(4,5) = tmat(0,2)*tmat(1,0)+tmat(1,2)*tmat(0,0);
  
      // the sixth row of transf. matrix
      t(5,0) = 2.0*tmat(0,1)*tmat(0,0);
      t(5,1) = 2.0*tmat(1,1)*tmat(1,0);
      t(5,2) = 2.0*tmat(2,1)*tmat(2,0);
      t(5,3) = tmat(1,1)*tmat(2,0)+tmat(2,1)*tmat(1,0);
      t(5,4) = tmat(0,1)*tmat(2,0)+tmat(2,1)*tmat(0,0);
      t(5,5) = tmat(0,1)*tmat(1,0)+tmat(1,1)*tmat(0,0);

      break;
    }
    case stress:
    {
      // Transformation matrix for strains in vector notation i.e. eps_g = T_eps eps_l

      // the first row of transf. matrix
      t(0,0) = tmat(0,0)*tmat(0,0);
      t(0,1) = tmat(1,0)*tmat(1,0);
      t(0,2) = tmat(2,0)*tmat(2,0);
      t(0,3) = 2.0*tmat(1,0)*tmat(2,0);
      t(0,4) = 2.0*tmat(0,0)*tmat(2,0);
      t(0,5) = 2.0*tmat(0,0)*tmat(1,0);
  
      // the second row of transf. matrix
      t(1,0) = tmat(0,1)*tmat(0,1);
      t(1,1) = tmat(1,1)*tmat(1,1);
      t(1,2) = tmat(2,1)*tmat(2,1);
      t(1,3) = 2.0*tmat(1,1)*tmat(2,1);
      t(1,4) = 2.0*tmat(0,1)*tmat(2,1);
      t(1,5) = 2.0*tmat(0,1)*tmat(1,1);
    
      // the third row of transf. matrix
      t(2,0) = tmat(0,2)*tmat(0,2);
      t(2,1) = tmat(1,2)*tmat(1,2);
      t(2,2) = tmat(2,2)*tmat(2,2);
      t(2,3) = 2.0*tmat(1,2)*tmat(2,2);
      t(2,4) = 2.0*tmat(0,2)*tmat(2,2);
      t(2,5) = 2.0*tmat(0,2)*tmat(1,2);
      
      // the fourth row of transf. matrix
      t(3,0) = tmat(0,2)*tmat(0,1);
      t(3,1) = tmat(1,2)*tmat(1,1);
      t(3,2) = tmat(2,2)*tmat(2,1);
      t(3,3) = tmat(1,2)*tmat(2,1)+tmat(2,2)*tmat(1,1);
      t(3,4) = tmat(0,2)*tmat(2,1)+tmat(2,2)*tmat(0,1);
      t(3,5) = tmat(0,2)*tmat(1,1)+tmat(1,2)*tmat(0,1);
      
      // the fifth row of transf. matrix
      t(4,0) = tmat(0,2)*tmat(0,0);
      t(4,1) = tmat(1,2)*tmat(1,0);
      t(4,2) = tmat(2,2)*tmat(2,0);
      t(4,3) = tmat(1,2)*tmat(2,0)+tmat(2,2)*tmat(1,0);
      t(4,4) = tmat(0,2)*tmat(2,0)+tmat(2,2)*tmat(0,0);
      t(4,5) = tmat(0,2)*tmat(1,0)+tmat(1,2)*tmat(0,0);
  
      // the sixth row of transf. matrix
      t(5,0) = tmat(0,1)*tmat(0,0);
      t(5,1) = tmat(1,1)*tmat(1,0);
      t(5,2) = tmat(2,1)*tmat(2,0);
      t(5,3) = tmat(1,1)*tmat(2,0)+tmat(2,1)*tmat(1,0);
      t(5,4) = tmat(0,1)*tmat(2,0)+tmat(2,1)*tmat(0,0);
      t(5,5) = tmat(0,1)*tmat(1,0)+tmat(1,1)*tmat(0,0);
     
      break;
    }
    default:
      print_err("unknown type of tensor quantity is required", __FILE__, __LINE__, __func__);
  }
  mixv(t, gf, l, i);  
}



/**
  The function transforms the fourth order stiffness tensor L given in the local coordinate system to 
  tensor G in global coordinate system. Both tensors fourth order tensors L and G are defined in %matrix 
  engineering notation (G(6,6), L(6,6)  but transformation %matrix T(3,3) must be given in the form 
  x_g = T x_l where x_l(3) and x_g(3) are local and global vectors respectively.
  
  @param g[out]   - 4-th order stiffness tensor in the %matrix engineering notation G(6,6) defined in 
                    the global coordinate system
  @param l[in]    - 4-th order stiffness tensor in the %matrix engineering notation L(6,6) defined in 
                     the local coordinate system
  @param tmat[in] - transformation %matrix T(3,3)
  @param ssst[in] - stress/strain state indicator - defines the number of g and l components
*/  
void lg_tens4transf (matrix &g, const matrix &l, const matrix &tmat, strastrestate /*ssst*/)
{
  matrix t(ASTCKMAT(6,6)), aux(ASTCKMAT(6,6));

//  matrix_tensor4(l, lf, ssst);

  // Transformation matrix glob->loc (T_eps)^{-1} for strains in the vector notation, i.e., eps_l = (T_eps)^{-1} eps_g

  // the first row of transf. matrix
  t(0,0) = tmat(0,0)*tmat(0,0);
  t(0,1) = tmat(1,0)*tmat(1,0);
  t(0,2) = tmat(2,0)*tmat(2,0);
  t(0,3) = tmat(1,0)*tmat(2,0);
  t(0,4) = tmat(0,0)*tmat(2,0);
  t(0,5) = tmat(0,0)*tmat(1,0);
  
  // the second row of transf. matrix
  t(1,0) = tmat(0,1)*tmat(0,1);
  t(1,1) = tmat(1,1)*tmat(1,1);
  t(1,2) = tmat(2,1)*tmat(2,1);
  t(1,3) = tmat(1,1)*tmat(2,1);
  t(1,4) = tmat(0,1)*tmat(2,1);
  t(1,5) = tmat(0,1)*tmat(1,1);
  
  // the third row of transf. matrix
  t(2,0) = tmat(0,2)*tmat(0,2);
  t(2,1) = tmat(1,2)*tmat(1,2);
  t(2,2) = tmat(2,2)*tmat(2,2);
  t(2,3) = tmat(1,2)*tmat(2,2);
  t(2,4) = tmat(0,2)*tmat(2,2);
  t(2,5) = tmat(0,2)*tmat(1,2);
  
  // the fourth row of transf. matrix
  t(3,0) = 2.0*tmat(0,2)*tmat(0,1);
  t(3,1) = 2.0*tmat(1,2)*tmat(1,1);
  t(3,2) = 2.0*tmat(2,2)*tmat(2,1);
  t(3,3) = tmat(1,2)*tmat(2,1)+tmat(2,2)*tmat(1,1);
  t(3,4) = tmat(0,2)*tmat(2,1)+tmat(2,2)*tmat(0,1);
  t(3,5) = tmat(0,2)*tmat(1,1)+tmat(1,2)*tmat(0,1);
  
  // the fifth row of transf. matrix
  t(4,0) = 2.0*tmat(0,2)*tmat(0,0);
  t(4,1) = 2.0*tmat(1,2)*tmat(1,0);
  t(4,2) = 2.0*tmat(2,2)*tmat(2,0);
  t(4,3) = tmat(1,2)*tmat(2,0)+tmat(2,2)*tmat(1,0);
  t(4,4) = tmat(0,2)*tmat(2,0)+tmat(2,2)*tmat(0,0);
  t(4,5) = tmat(0,2)*tmat(1,0)+tmat(1,2)*tmat(0,0);
  
  // the sixth row of transf. matrix
  t(5,0) = 2.0*tmat(0,1)*tmat(0,0);
  t(5,1) = 2.0*tmat(1,1)*tmat(1,0);
  t(5,2) = 2.0*tmat(2,1)*tmat(2,0);
  t(5,3) = tmat(1,1)*tmat(2,0)+tmat(2,1)*tmat(1,0);
  t(5,4) = tmat(0,1)*tmat(2,0)+tmat(2,1)*tmat(0,0);
  t(5,5) = tmat(0,1)*tmat(1,0)+tmat(1,1)*tmat(0,0);

  // transformation of the stiffness tensor to the global coordinate system: G = ((T_{eps})^{-1})^T L (T_{eps}^{-1})
  // aux = ((T_{eps})^{-1})^T L
  mtxm(t, l, aux);
  // G = aux ((T_{eps})^{-1})
  mxm(aux, t, g);

  // convert G_f to reduced form G depending on ssst
//  tensor4_matrix(g, gf, ssst);
}



/**
  The function transforms the fourth order stiffness tensor G given in the global coordinate system to 
  tensor L in local coordinate system. Both tensors fourth order tensors L and G are defined in %matrix 
  engineering notation (G(6,6), L(6,6)  but transformation %matrix T(3,3) must be given in the form 
  x_g = T x_l where x_l(3) and x_g(3) are local and global vectors respectively.
  
  @param[in] g - 4-th order stiffness tensor in the %matrix engineering notation g(6,6) defined in 
                 the global coordinate system
  @param[out] l - 4-th order stiffness tensor in the %matrix engineering notation l(6,6) defined in 
                  the local coordinate system (input)
  @param[in] tmat - transformation %matrix T(3,3)
  @param[in] ssst - stress/strain state indicator - defines the number of g and l components
*/  
void gl_tens4transf (const matrix &g, matrix &l, const matrix &tmat, strastrestate /*ssst*/)
{
  matrix t(ASTCKMAT(6,6)), aux(ASTCKMAT(6,6));

//  matrix_tensor4(l, lf, ssst);

  // Transformation matrix loc->glob T_{eps} for strains in vector notation, i.e., eps_g = T_{eps} eps_l

  // the first row of transf. matrix
  t(0,0) = tmat(0,0)*tmat(0,0);
  t(0,1) = tmat(0,1)*tmat(0,1);
  t(0,2) = tmat(0,2)*tmat(0,2);
  t(0,3) = tmat(0,2)*tmat(0,1);
  t(0,4) = tmat(0,2)*tmat(0,0);
  t(0,5) = tmat(0,1)*tmat(0,0);
  
  // the second row of transf. matrix
  t(1,0) = tmat(1,0)*tmat(1,0);
  t(1,1) = tmat(1,1)*tmat(1,1);
  t(1,2) = tmat(1,2)*tmat(1,2);
  t(1,3) = tmat(1,2)*tmat(1,1);
  t(1,4) = tmat(1,2)*tmat(1,0);
  t(1,5) = tmat(1,1)*tmat(1,0);
  
  // the third row of transf. matrix
  t(2,0) = tmat(2,0)*tmat(2,0);
  t(2,1) = tmat(2,1)*tmat(2,1);
  t(2,2) = tmat(2,2)*tmat(2,2);
  t(2,3) = tmat(2,2)*tmat(2,1);
  t(2,4) = tmat(2,2)*tmat(2,0);
  t(2,5) = tmat(2,1)*tmat(2,0);
  
  // the fourth row of transf. matrix
  t(3,0) = 2.0*tmat(1,0)*tmat(2,0);
  t(3,1) = 2.0*tmat(1,1)*tmat(2,1);
  t(3,2) = 2.0*tmat(1,2)*tmat(2,2);
  t(3,3) = tmat(1,2)*tmat(2,1)+tmat(2,2)*tmat(1,1);
  t(3,4) = tmat(1,2)*tmat(2,0)+tmat(2,2)*tmat(1,0);
  t(3,5) = tmat(1,1)*tmat(2,0)+tmat(2,1)*tmat(1,0);
  
  // the fifth row of transf. matrix
  t(4,0) = 2.0*tmat(0,0)*tmat(2,0);
  t(4,1) = 2.0*tmat(0,1)*tmat(2,1);
  t(4,2) = 2.0*tmat(0,2)*tmat(2,2);
  t(4,3) = tmat(0,2)*tmat(2,1)+tmat(2,2)*tmat(0,1);
  t(4,4) = tmat(0,2)*tmat(2,0)+tmat(2,2)*tmat(0,0);
  t(4,5) = tmat(0,1)*tmat(2,0)+tmat(2,1)*tmat(0,0);
  
  // the sixth row of transf. matrix
  t(5,0) = 2.0*tmat(0,0)*tmat(1,0);
  t(5,1) = 2.0*tmat(0,1)*tmat(1,1);
  t(5,2) = 2.0*tmat(0,2)*tmat(1,2);
  t(5,3) = tmat(0,2)*tmat(1,1)+tmat(1,2)*tmat(0,1);
  t(5,4) = tmat(0,2)*tmat(1,0)+tmat(1,2)*tmat(0,0);
  t(5,5) = tmat(0,1)*tmat(1,0)+tmat(1,1)*tmat(0,0);

  // transformation of the stiffness tensor to the global coordinate system: L = (T_{eps})^T G T_{eps}
  // aux = (T_{eps})^T G
  mtxm (t, g, aux);
  // G = aux T_{eps}
  mxm (aux, t, l);

  // convert G_f to reduced form G depending on ssst
//  tensor4_matrix(g, gf, ssst);
}



/**
   Function converts long vectors (with 6 components) to short vectors 
   (with appropriate number of components)
   e.g. in the case of plane stress, vectors with 6 components are converted 
   to vectors with 3 components
   
   ordering of components is shown on the example of stresses
   lv(0) = sigma_x
   lv(1) = sigma_y
   lv(2) = sigma_z
   lv(3) = tau_yz
   lv(4) = tau_zx
   lv(5) = tau_xy
   
   @param[in]  lv - long %vector
   @param[out] sv - short %vector (output)
   @param[in]  ssst - strain/stress state
   
   @return The function returns required short %vector in the parameter sv.
   
   Created by JK, 17.2.2007
*/
void longvect_shortvect (vector &lv,vector &sv,strastrestate ssst)
{
  nullv (sv);
  
  switch (ssst){
  case bar:{
    sv(0)=lv(0);
    break;
  }
  case plbeam:{
    break;
  }
  case planestress:{
    sv(0)=lv(0);
    sv(1)=lv(1);
    sv(2)=lv(5);
    break;
  }
  case planestrain:{
    sv(0)=lv(0);
    sv(1)=lv(1);
    sv(2)=lv(5);
    break;
  }
  case axisymm:{
    sv(0)=lv(0);
    sv(1)=lv(1);
    sv(2)=lv(2);
    sv(3)=lv(5);
    break;
  }
  case spacestress:{
    sv(0)=lv(0);
    sv(1)=lv(1);
    sv(2)=lv(2);
    sv(3)=lv(3);
    sv(4)=lv(4);
    sv(5)=lv(5);
    break;
  }
  default:{
    print_err("unknown strain state is required", __FILE__, __LINE__, __func__);
  }
  }
}

// **************************************************************
// **************************************************************
// **************************************************************
// **************************************************************

/**
   function computes second invariant of a deviator of a stress tensor
   the input %vector contains the engineering components of the stress tensor
   components of the deviator are not needed

   ordering of components is shown on the example of stresses
   lv(0) = sigma_x
   lv(1) = sigma_y
   lv(2) = sigma_z
   lv(3) = tau_yz
   lv(4) = tau_zx
   lv(5) = tau_xy
   
   @param[in] lv - %vector of stress components
   
   17. 4. 2015, JK
*/
double j2_stress_invar (vector &lv)
{
  double j2;
  
  j2 = ((lv(0)-lv(1))*(lv(0)-lv(1)) + (lv(1)-lv(2))*(lv(1)-lv(2)) + (lv(2)-lv(0))*(lv(2)-lv(0)))/6.0;
  j2 += lv(3)*lv(3) + lv(4)*lv(4) + lv(5)*lv(5);
  
  return j2;
}



/**
  The function computes second invariant of a deviator of a strain tensor.
  The input %vector contains engineering components of the strain tensor.
  Components of the strain deviator are not needed on input.

  Ordering of components is shown on the example of strains
  lv(0) = eps_x
  lv(1) = eps_y
  lv(2) = eps_z
  lv(3) = gamma_yz
  lv(4) = gamma_zx
  lv(5) = gamma_xy
   
  @param[in] lv - %vector of strain components in Voigt's notation
   
  @return The function returns the resulting deviator in the argument dev.

  17. 4. 2015, JK
*/
double j2_strain_invar (vector &lv)
{
  double j2;
  
  j2 = ((lv(0)-lv(1))*(lv(0)-lv(1)) + (lv(1)-lv(2))*(lv(1)-lv(2)) + (lv(2)-lv(0))*(lv(2)-lv(0)))/6.0;
  j2 += (lv(3)*lv(3) + lv(4)*lv(4) + lv(5)*lv(5))/4.0;
  
  return j2;
}

/**
  The function computes deviatoric part of the given tensor in Voigt's notation.
   
  @param[in] tens - engineering components of the tensor, i.e. tensor in Voigt's notation
  @param[out] dev - engineering components of the deviator of the tensor, i.e. deviator in Voigt's notation
   
  Both vectors must contain 6 components, because this function cannot know the actual number of components,
  the actual components have to be determined outside of this function.
   
  @return The function returns the resulting deviator in the argument dev.

  28.10.2001, modified 21. 4. 2015
*/
void deviator (vector &tens,vector &dev)
{
  double m;
  
  //  isotropic part of the tensor
  m = (tens(0) + tens(1) + tens(2))/3.0;
  
  dev(0) = tens(0) - m;
  dev(1) = tens(1) - m;
  dev(2) = tens(2) - m;
  dev(3) = tens(3);
  dev(4) = tens(4);
  dev(5) = tens(5);
}



/**
  The function computes first invariant of a tensor.
  Input is a %vector with 6 engineering components of the tensor, i.e. tensor in Voigt's notation.
   
  @param[in] tens - engineering components of the tensor
   
  tens(0) = tens_11
  tens(1) = tens_22
  tens(2) = tens_33
  tens(3) = tens_23
  tens(4) = tens_32
  tens(5) = tens_13

  @return The function returns a value of the first invariant of the given tensor tens. 

  4.1.2002
*/
double first_invar (vector &tens)
{
  double inv;
  inv=tens(0)+tens(1)+tens(2);
  return inv;
}



/**
  The function computes second invariant of a stress tensor.
  Input is a %vector with six engineering components of the tensor, i.e. tensor in Voigt's notation.
   
  @param[in] tens - engineering components of the stress tensor
   
  @return The function returns a value of the first invariant of the given tensor tens. 

  4.1.2002
*/
double second_stress_invar (vector &tens)
{
  double inv;

  inv = tens(0)*tens(1) + tens(1)*tens(2) + tens(2)*tens(0);
  inv = inv - tens(3)*tens(3) - tens(4)*tens(4) - tens(5)*tens(5);
  return inv;
}



/**
  The function computes third invariant of a stress tensor.
  Input is a %vector with six engineering components of the tensor, i.e. tensor in Voigt's notation.
   
  @param[in] tens - engineering components of stress tensor
   
  @return The function returns a value of the first invariant of the given tensor tens. 

  4.1.2002
*/
double third_stress_invar (vector &tens)
{
  double inv;

  //  scheme of engineering components in tensor notation
  //  0 5 4
  //  5 1 3
  //  4 3 2
  
  inv = tens(0)*tens(1)*tens(2) + 2.0*tens(3)*tens(4)*tens(5);
  inv = inv - tens(0)*tens(3)*tens(3) - tens(1)*tens(4)*tens(4) - tens(2)*tens(5)*tens(5);
  
  return inv;
}



/**
   function computes second invariant of a strain tensor
   
   @param[in] tens - engineering components of the strain tensor
   
   4.1.2002
*/
double second_strain_invar (vector &tens)
{
  double inv;

  inv = tens(0)*tens(1) + tens(1)*tens(2) + tens(2)*tens(0);
  inv = inv - 1.0/4.0 *(tens(3)*tens(3) + tens(4)*tens(4) + tens(5)*tens(5));
  return inv;
}

/**
   function computes third invariant of a strain tensor
   
   @param[in] tens - engineering components of strain tensor
   
   4.1.2002
*/
double third_strain_invar (vector &tens)
{
  double inv;

  //  scheme of engineering components in tensor notation
  //  0 5 4
  //  5 1 3
  //  4 3 2
  
  inv = tens(0)*tens(1)*tens(2) + 1.0/4.0*(tens(3)*tens(4)*tens(5) - tens(0)*tens(3)*tens(3) - tens(1)*tens(4)*tens(4) - tens(2)*tens(5)*tens(5));
  
  return inv;
}



/**
  The function returns norm of stress tensor given in Voigt notation.

  @param[in] sig - stress components in Voigt notation
  
  @return The function returns norm of stress tensor.

  Created by Tomas Koudelka, 12.2015
*/
double tensor_stress_norm(vector &sig)
{
  long i;
  double norm = 0.0;

  for (i=0; i<3; i++)
    norm += sig(i) * sig(i);
  for (i=3; i<6; i++)
    norm += 2.0 * sig(i) * sig(i);

  norm = sqrt(norm);

  return norm;
}



/**
  The function returns norm of strain tensor given in Voigt notation.

  @param eps - strain components in Voigt notation
  
  @return The function returns norm of strain tensor.

  Created by Tomas Koudelka, 12.2015
*/
double tensor_strain_norm(vector &eps)
{
  long i;
  double norm = 0.0;

  for (i=0; i<3; i++)
    norm += eps(i) * eps(i);
  for (i=3; i<6; i++)
    norm += 0.5 * eps(i) * eps(i);

  norm = sqrt(norm);

  return norm;
}



/**
  The function computes normed stress tensor nsig from the given stress tensor sig.
  Both are given in Voigt notation.
  
  @param sig - stress components in Voigt notation
  @param nsig - resulting normed stress tensor in Voigt notation (output)

  @return The funtion returns the resulting normed stress tensor in the argument nsig.

  Created by Tomas Koudelka, 12.2015
*/
void normed_stress_tensor(vector &sig, vector &nsig)
{
  long i;
  double norm = 0.0;
  for (i=0; i<3; i++)
    norm += sig(i) * sig(i);
  for (i=3; i<6; i++)
    norm += 2.0 * sig(i) * sig(i);

  norm = sqrt(norm);

  for (i=0; i<6; i++)
    nsig(i) = sig(i) / norm;

  return;
}


/**
  The function computes normed strain tensor neps from the given strain tensor eps.
  Both are given in Voigt notation.
  
  @param eps - strain components in Voigt notation
  @param neps - resulting normed strain tensor in Voigt notation (output)

  @return The funtion returns the resulting normed strain tensor in the argument neps.

  Created by Tomas Koudelka, 12.2015
*/
void normed_strain_tensor(vector &eps, vector &neps)
{
  long i;
  double norm = 0.0;
  for (i=0; i<3; i++)
    norm += eps(i) * eps(i);
  for (i=3; i<6; i++)
    norm += 0.5 * eps(i) * eps(i);

  norm = sqrt(norm);

  for (i=0; i<6; i++)
    neps(i) = eps(i) / norm;

  return;
}


/**
  The function computes the first invariant of the second order tensor.

  @param[in] tens - tensor components stored as a %matrix 3x3

  @return The function returns the first invariant value.
*/
double first_invar (matrix &t)
{
  double inv = t(0,0) + t(1,1) + t(2,2);

  return inv;
}



/**
  The function computes the second invariant of the second order tensor t, which 
  is defined as T11*T22 + T11*T33 + T22*T33 - T12^2 - T13^2 - T23^2.

  @param[in] tens - tensor components stored as a %matrix 3x3

  @return The function returns the second invariant value.
*/
double second_invar (matrix &t)
{
  double inv = t(0,0)*t(1,1) + t(0,0)*t(2,2) + t(1,1)*t(2,2);
  inv -= t(0,1)*t(0,1) + t(0,2)*t(0,2) + t(1,2)*t(1,2);

  return inv;
}



/**
  The function computes the third invariant of the second order tensor, which is
  defined as det|T|.

  @param[in] tens - tensor components stored as a %matrix 3x3

  @return The function returns the third invariant value.
*/
double third_invar (matrix &t)
{
  double inv = t(0,0)*t(1,1)*t(2,2) + t(0,1)*t(1,2)*t(2,0) + t(0,2)*t(1,0)*t(2,1);
  inv -= t(0,2)*t(1,1)*t(2,0) + t(0,0)*t(1,2)*t(2,1) + t(0,1)*t(1,0)*t(2,2);

  return inv;
}
