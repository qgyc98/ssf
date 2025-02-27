#include "therisomat.h"
#include "global.h"
#include "mechmat.h"
#include "matrix.h"
#include "intpoints.h"



/**
  The constructor inializes attributes to zero values.
  
  Created by JK,
*/
therisomat::therisomat (void)
{
  alpha=0.0;
}



/**
  The destructor is defined only for the formal purposes.
  
  Created by JK,
*/
therisomat::~therisomat (void)
{

}



/**
  The function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened input text file
  
  @return The function does not return anything.

  Created by JK,
*/
void therisomat::read (XFILE *in)
{
  xfscanf (in,"%k%lf","alpha",&alpha);
}

/**
  The function prints material parameters into the opened text file given
  by the parameter out.

  @param out - pointer to the opened input text file
  
  @return The function does not return anything.

  Created by TKr, 02/01/2013
*/
void therisomat::print (FILE *out)
{
  fprintf (out,"%le",alpha);
}



/**
  The function assembles thermal dilatancy %matrix for the given stress/strain state.
 
  @param d - dilatancy %matrix (output)
  @param ssst - stress/strain state indicator

  @return The function returns required %matrix in the parameter d.

  Created by JK,
*/
void therisomat::matdilat (matrix &d,strastrestate ssst)
{
  switch (ssst){
  case bar:{
    matdilat_bar (d);
    break;
  }
  case plbeam:{
    matdilat_plbeam (d);
    break;
  }
  case planestress:{
    matdilat_plstress (d);
    break;
  }
  case planestrain:{
    matdilat_plstrain (d);
    break;
  }
  case platek:{
    matdilat_plate (d);
    break;
  }
  case plates:{
    matdilat_plate (d);
    break;
  }
  case axisymm:{
    matdilat_axi (d);
    break;
  }
  case spacestress:{
    matdilat_spacestr (d);
    break;
  }
  default:
    print_err("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function creates thermal %dilatancy %matrix of the elastic
  isotropic material for bar elements.
   
  @param d - thermal dilatancy %matrix of the material (output)
   
  @return The function returns required %matrix in the parameter d.

  Created by JK, 11.9.2001
*/
void therisomat::matdilat_bar (matrix &d)
{
  d[0][0] = alpha;
}



/**
  Function creates thermal dilatancy %matrix of the elastic
  isotropic material for plane beam elements.
   
  @param d - thermal dilatancy %matrix of the material
   
  @return The function returns required %matrix in the parameter d.

  Created by JK, 11.9.2001
*/
void therisomat::matdilat_plbeam (matrix &d)
{
  d[0][0] = alpha;
  d[1][1] = 0.0;
  d[2][2] = 0.0;//tady oprava
  //d[2][2] = alpha;
}



/**
  Function creates thermal dilatancy %matrix of the elastic
  isotropic material for 2D problems (plane stress).

  @param d - thermal dilatancy matrix of the material

  @return The function returns required %matrix in the parameter d.

  Created by JK, 19.7.2001
*/
void therisomat::matdilat_plstress (matrix &d)
{
  fillm(0.0,d);
  
  d[0][0] = alpha;  d[0][1] = 0.0;    d[0][2] = 0.0;
  d[1][0] = 0.0;    d[1][1] = alpha;  d[1][2] = 0.0;
  d[2][0] = 0.0;    d[2][1] = 0.0;    d[2][2] = 0.0;
}



/**
  Function creates thermal dilatancy %matrix of the elastic
  isotropic material for 2D problems (plane strain).

  @param d - thermal dilatancy %matrix of the material

  @return The function returns required %matrix in the parameter d.

  Created by JK, 19.7.2001
*/
void therisomat::matdilat_plstrain (matrix &d)
{
  fillm(0.0,d);
  
  d[0][0] = alpha;  d[0][1] = 0.0;     d[0][2] = 0.0;
  d[1][0] = 0.0;    d[1][1] = alpha;   d[1][2] = 0.0;
  d[2][0] = 0.0;    d[2][1] = 0.0;     d[2][2] = 0.0;
}



/**
  Function creates thermal dilatancy %matrix of the elastic
  isotropic material for 2D problems (axisymmetric problem).

  @param d - thermal dilatancy %matrix of the material

  @return The function returns required %matrix in the parameter d.

  Created by JK, 19.7.2001
*/
void therisomat::matdilat_axi (matrix &d)
{
  fillm(0.0,d);
  
  d[0][0]=alpha;  d[0][1]=0.0;    d[0][2]=0.0;
  d[1][0]=0.0;    d[1][1]=alpha;  d[1][2]=0.0;
  d[2][0]=0.0;    d[2][1]=0.0;    d[2][2]=alpha;

  d[3][3]=0.0;
}



/**
  Function creates thermal dilatancy %matrix of the elastic
  isotropic material for plate elements.
   
  @param d - thermal dilatancy %matrix
   
  @return The function returns required %matrix in the parameter d.

  Created by JK, 19.7.2001
*/
void therisomat::matdilat_plate (matrix &d)
{
  fillm(0.0,d);
  
  d[0][0]=alpha;  d[1][1]=alpha;
}

/**
  Function creates thermal dilatancy %matrix of the elastic
  isotropic material for 3D problems.
   
  @param d - thermal dilatancy %matrix of the material

  @return The function returns required %matrix in the parameter d.

  Created by JK, 19.7.2001
*/
void therisomat::matdilat_spacestr (matrix &d)
{
  fillm(0.0,d);
  
  d[0][0]=alpha;  d[1][1]=alpha;  d[2][2]=alpha;
}



/**
  Function computes strains caused by temperature changes at the given integration point.
  Function assumes temperature changes from the array Mm->nonmechq
  as difference between quantity 'temperature' and 'initial_temperature'.

  @param ipp - integration point id
   
  Created by JK, 4.6.2005
*/
void therisomat::temprstrains (long ipp)
{
  double dt;
  strastrestate sss;
  
  //  type of strain/stress state
  sss = Mm->ip[ipp].ssst;
  
  //  change of temperature
  dt = Mm->givenonmechq(temperature, ipp) - Mm->givenonmechq(initial_temperature, ipp);

  // thermal strains are stored, temperature change is calculated according to the actual value
  // stored in the nomechq arrays
  switch (sss){
  case bar:{
    Mm->tempstrains[ipp][0]=alpha*dt;
    break;
  }
  case plbeam:{
    Mm->tempstrains[ipp][0]=alpha*dt;
    Mm->tempstrains[ipp][1]=0.0;
    Mm->tempstrains[ipp][2]=0.0;//tady oprava
    //Mm->tempstrains[ipp][2]=alpha*dt;
    break;
  }
  case planestress:{
    Mm->tempstrains[ipp][0]=alpha*dt;
    Mm->tempstrains[ipp][1]=alpha*dt;
    Mm->tempstrains[ipp][3]=alpha*dt;
    break;
  }
  case planestrain:{
    Mm->tempstrains[ipp][0]=alpha*dt;
    Mm->tempstrains[ipp][1]=alpha*dt;
    Mm->tempstrains[ipp][3]=alpha*dt;
    break;
  }
  case axisymm:{
    Mm->tempstrains[ipp][0]=alpha*dt;
    Mm->tempstrains[ipp][1]=alpha*dt;
    Mm->tempstrains[ipp][2]=alpha*dt;
    break;
  }
  case spacestress:{
    Mm->tempstrains[ipp][0]=alpha*dt;
    Mm->tempstrains[ipp][1]=alpha*dt;
    Mm->tempstrains[ipp][2]=alpha*dt;
    Mm->tempstrains[ipp][3]=0.0;
    Mm->tempstrains[ipp][4]=0.0;
    Mm->tempstrains[ipp][5]=0.0;    
    break;
  }
  default:
    print_err("unknown strain/stress state is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function computes strains caused by temperature changes at the given integration point
  and accumulates them in the Mm->tempstrains array.
  Function assumes temperature changes from the array Mm->nonmechq
  as difference between quantity 'temperature' and 'initial_temperature'.

  @param ipp - integration point id
   
  Created by TKo, 4.10.2016
*/
void therisomat::cumultemprstrains (long ipp)
{
  double dt;
  strastrestate sss;
  
  //  type of strain/stress state
  sss = Mm->ip[ipp].ssst;
  
  //  change of temperature
  dt = Mm->givenonmechq(temperature, ipp) - Mm->givenonmechq(initial_temperature, ipp);

  // Thermal strain contributions are added because of different contributions at subloadcases 
  // in time dependent problems.
  // Mm->tempstrains array is cleaned at mefel_right_hand_side function before the 
  // given load case is being assembled.
  switch (sss){
  case bar:{
    Mm->tempstrains[ipp][0]+=alpha*dt;
    break;
  }
  case plbeam:{
    Mm->tempstrains[ipp][0]+=alpha*dt;
    Mm->tempstrains[ipp][1]+=0.0;
    Mm->tempstrains[ipp][2]+=0.0;//tady oprava
    //Mm->tempstrains[ipp][2]+=alpha*dt;
    break;
  }
  case planestress:{
    Mm->tempstrains[ipp][0]+=alpha*dt;
    Mm->tempstrains[ipp][1]+=alpha*dt;
    Mm->tempstrains[ipp][3]+=alpha*dt;
    break;
  }
  case planestrain:{
    Mm->tempstrains[ipp][0]+=alpha*dt;
    Mm->tempstrains[ipp][1]+=alpha*dt;
    Mm->tempstrains[ipp][3]+=alpha*dt;
    break;
  }
  case axisymm:{
    Mm->tempstrains[ipp][0]+=alpha*dt;
    Mm->tempstrains[ipp][1]+=alpha*dt;
    Mm->tempstrains[ipp][2]+=alpha*dt;
    break;
  }
  case spacestress:{
    Mm->tempstrains[ipp][0]+=alpha*dt;
    Mm->tempstrains[ipp][1]+=alpha*dt;
    Mm->tempstrains[ipp][2]+=alpha*dt;
    Mm->tempstrains[ipp][3]+=0.0;
    Mm->tempstrains[ipp][4]+=0.0;
    Mm->tempstrains[ipp][5]+=0.0;    
    break;
  }
  default:
    print_err("unknown strain/stress state is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The funtion marks required non-mechanical quantities in the array anmq.

  @param anmq - array with flags for used material types
                anmq[i] = 1 => qunatity type nonmechquant(i+1) is required
                anmq[i] = 0 => qunatity type nonmechquant(i+1) is not required

  @return The function does not return anything, but it may change content of anmq array.
*/
void therisomat::give_reqnmq(long *anmq)
{
  anmq[temperature-1] = 1;
  anmq[initial_temperature-1] = 1;
}
