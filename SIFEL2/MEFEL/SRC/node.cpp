#include "node.h"
#include "vector.h"
#include "iotools.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include <string.h>



/**
  Constructor initializes data members to zero or default values.

  Created by JK,
*/
node::node (void)
{
  //  type of cross section
  crst = (crsectype) 0;
  //  cross section id
  idcs=0;
  
  //  indicator of local coordinate system on node
  transf=0;
  //  base vectors of local coordinate system
  e1=NULL;  e2=NULL;  e3=NULL;
  
  react=0;  r=NULL;
  
  //  number of strain/stress components
  ncompstr=0;
  //  number of components in array other
  ncompother=0;
  //  array containing nodal strains
  strain=NULL;
  //  array containing nodal stresses
  stress=NULL;
  //  array containing nodal values of the array other
  other=NULL;
  //  number of contributions to the array strain
  ncontr_strain=NULL;
  //  number of contributions to the array stress
  ncontr_stress=NULL;
  //  number of contributions to the array other
  ncontr_other=0;
  //  volume used for strain contributions
  vol_strain=NULL;
  //  volume used for stress contributions
  vol_stress=NULL;
  //  volume used for other contributions
  vol_other=0;
  
  
  pstra = NULL;
  pstre = NULL;
  meaning = NULL;

  nodval = NULL;
}



/**
  Destructor releases allocated memory of the node object.

  Created by JK,
*/
node::~node (void)
{
  delete [] e1;  delete [] e2;  delete [] e3;

  delete [] strain;  delete [] stress;  delete [] other;
  delete [] ncontr_strain;  delete [] ncontr_stress;
  delete [] vol_strain;  delete [] vol_stress;

  delete [] r;  
  delete [] pstra;
  delete [] pstre;
  delete [] meaning;
  delete [] nodval;
}



/**
  Function reads nodal data from opened text file.
   
  @param in - pointer to teh opened text file
   
  @return The function does not return anything.

  Created by JK
*/
void node::read (XFILE *in)
{
  //  type of cross section
  xfscanf (in, "%m", &crsectype_kwdset, (int*)&crst);
  if (crst!=0){
    xfscanf (in,"%ld",&idcs);
    idcs--;
  }
  
  //  transformation to local system
  xfscanf (in,"%ld",&transf);
  if (transf!=0 && transf!=2 && transf!=3)
    print_err("wrong identification of local system", __FILE__, __LINE__, __func__);
  
  if (transf==2){
    e1 = new double [2];  e2 = new double [2];
    xfscanf (in,"%lf %lf",e1+0,e1+1);
    xfscanf (in,"%lf %lf",e2+0,e2+1);
  }
  if (transf==3){
    e1 = new double [3];  e2 = new double [3];  e3 = new double [3];
    xfscanf (in,"%lf %lf %lf",e1+0,e1+1,e1+2);
    xfscanf (in,"%lf %lf %lf",e2+0,e2+1,e2+2);
    xfscanf (in,"%lf %lf %lf",e3+0,e3+1,e3+2);
  }

}

/**
  Function prints nodal data into opened text file.
   
  @param out - pointer to teh opened text file
   
  @return The function does not return anything.

  Created by TKr, 02/01/2013
*/
void node::print (FILE *out)
{
  //  type of cross section
  fprintf (out,"  %d", int(crst));
  if (crst!=0){
    fprintf (out," %ld",idcs+1);
  }

  //  transformation to local system
  fprintf (out,"  %ld",transf);

  if (transf==2){
    fprintf (out," %lf %lf",e1[0],e1[1]);
    fprintf (out," %lf %lf",e2[0],e2[1]);
  }
  if (transf==3){
    fprintf (out," %lf %lf %lf",e1[0],e1[1],e1[2]);
    fprintf (out," %lf %lf %lf",e2[0],e2[1],e2[2]);
    fprintf (out," %lf %lf %lf",e3[0],e3[1],e3[2]);
  }
  fprintf (out,"\n");
}



/**
  Function allocates arrays strain, stress and other on nodes
   
  @param ncomp - number of strain/stress components
  @param ncompo - number of components of array other
  @param nlc - number of load cases
   
  @return The function does not return anything.

  Created by JK, 18.5.2002, 
  Revised by JK, 29.11.2006
*/
void node::alloc (long ncomp,long ncompo,long nlc)
{
  //  number of components of strain/stress
  ncompstr=ncomp;
  //  number of components of array other
  ncompother=ncompo;
  
  if (strain!=NULL)
    delete [] strain;
  if (stress!=NULL)
    delete [] stress;
  if (other!=NULL)
    delete [] other;
  
  strain = new double [nlc*ncompstr];
  stress = new double [nlc*ncompstr];
  other  = new double [nlc*ncompother];

  memset(strain, 0, sizeof(*strain)*nlc*ncompstr);
  memset(stress, 0, sizeof(*stress)*nlc*ncompstr);
  memset(other, 0, sizeof(*other)*nlc*ncompother);
  
  if (ncontr_strain!=NULL)
    delete [] ncontr_strain;
  if (ncontr_stress!=NULL)
    delete [] ncontr_stress;
  
  ncontr_strain = new long [nlc];
  ncontr_stress = new long [nlc];

  memset(ncontr_strain, 0, sizeof(*ncontr_strain)*nlc);
  memset(ncontr_stress, 0, sizeof(*ncontr_stress)*nlc);
}



/**
  Function allocates arrays strain, ncontr_strain and vol_strain.

  @param ncomp - number of stress/strain components
  @param nlc - number of load cases
   
  @return The function does not return anything.

  Created by JK, 16.7.2008
*/
void node::alloc_strain (long ncomp,long nlc)
{
  //  number of components of strain/stress
  ncompstr=ncomp;
  
  if (strain!=NULL)
    delete [] strain;
  
  strain = new double [nlc*ncompstr];
  memset(strain, 0, sizeof(*strain)*nlc*ncompstr);
  
  if (Mp->strainaver==1){
    if (ncontr_strain!=NULL)
      delete [] ncontr_strain;
    ncontr_strain = new long [nlc];
    memset(ncontr_strain, 0, sizeof(*ncontr_strain)*nlc);
  }
  if (Mp->strainaver==2){
    if (vol_strain!=NULL)
      delete [] vol_strain;
    vol_strain = new double [nlc];
    memset(vol_strain, 0, sizeof(*vol_strain)*nlc);
  }
}



/**
  Function allocates array stress, ncontr_stress and vol_stress.
   
  @param ncomp - number of stress/strain components
  @param nlc - number of load cases
   
  @return The function does not return anything.

  Created by JK, 16.7.2008
*/
void node::alloc_stress (long ncomp,long nlc)
{
  //  number of components of strain/stress
  ncompstr=ncomp;
  
  if (stress!=NULL)
    delete [] stress;
  
  stress = new double [nlc*ncompstr];
  memset(stress, 0, sizeof(*stress)*nlc*ncompstr);
  
  if (Mp->stressaver==1){
    if (ncontr_stress!=NULL)
      delete [] ncontr_stress;
    ncontr_stress = new long [nlc];
    memset(ncontr_stress, 0, sizeof(*ncontr_stress)*nlc);
  }
  if (Mp->stressaver==2){
    if (vol_stress!=NULL)
      delete [] vol_stress;
    vol_stress = new double [nlc];
    memset(vol_stress, 0, sizeof(*vol_stress)*nlc);
  }
}



/**
  Function allocates array other, ncontr_other and vol_other.
   
  @param ncomp - number of stress/strain components
   
  @return The function does not return anything.

  Created by JK, 16.7.2008
  Modified by TKo 14.11.2013
*/
void node::alloc_other (long ncompo)
{
  //  number of components of array other
  ncompother=ncompo;
  
  if (other!=NULL)
    delete [] other;
  
  other = new double [ncompother];
  memset(other, 0, sizeof(*other)*ncompother);
}



/**
  Function stores strain components.
   
  @param lcid - load case id
  @param fi - first index
  @param eps - %vector containing part of strain components

  @return The function does not return anything.

  Created by JK, 19.5.2002
*/
void node::storestrain (long lcid,long fi,vector &eps)
{
  long i,j,m,n;
  m=lcid*ncompstr;
  n=eps.n;
  
  if (fi==0)  ncontr_strain[lcid]++;
  
  j=0;
  for (i=m+fi;i<m+fi+n;i++){
    strain[i]+=eps[j];
    j++;
  }
}



/**
  Function stores strain components.
   
  @param lcid - load case id
  @param fi - first index
  @param vol - volume used for contribution
  @param eps - %vector containing part of strain components
   
  @return The function does not return anything.

  Created by JK, 19.5.2002
*/
void node::storestrain (long lcid,long fi,double vol,vector &eps)
{
  long i,j,m,n;
  m=lcid*ncompstr;
  n=eps.n;
  
  if (fi==0)  vol_strain[lcid]+=vol;
  
  j=0;
  for (i=m+fi;i<m+fi+n;i++){
    strain[i]+=eps[j];
    j++;
  }
}



/**
  Function stores strain components.
   
  @param lcid - load case id
  @param fi - first index
  @param ncomp - number of components
  @param eps - %vector containing part of strain components
   
  @return The function does not return anything.

  Created by JK, 19.5.2002
*/
void node::storestrain (long lcid,long fi,long ncomp,vector &eps)
{
  long i,j,m;
  m=lcid*ncompstr;

  if (fi==0)  ncontr_strain[lcid]++;
  
  j=fi;
  for (i=m+fi;i<m+fi+ncomp;i++){
    strain[i]+=eps[j];
    j++;
  }
}



/**
  Function stores stress components.
   
  @param lcid - load case id
  @param fi - first index
  @param sig - %vector containing part of stress components
   
  @return The function does not return anything.

  Created by JK, 25.9.2004
*/
void node::storestress (long lcid,long fi,vector &sig)
{
  long i,j,m,n;
  m=lcid*ncompstr;
  n=sig.n;
  
  if (fi==0)  ncontr_stress[lcid]++;
  
  j=0;
  for (i=m+fi;i<m+fi+n;i++){
    stress[i]+=sig[j];
    j++;
  }
}



/**
  Function stores stress components.
   
  @param lcid - load case id
  @param fi - first index
  @param vol - volume used for contribution
  @param sig - %vector containing part of stress components
   
  @return The function does not return anything.

  Created by JK, 25.9.2004
*/
void node::storestress (long lcid,long fi,double vol,vector &sig)
{
  long i,j,m,n;
  m=lcid*ncompstr;
  n=sig.n;
  
  if (fi==0)  vol_stress[lcid]+=vol;
  
  j=0;
  for (i=m+fi;i<m+fi+n;i++){
    stress[i]+=sig[j];
    j++;
  }
}



/**
  Function stores stress components.
   
  @param lcid - load case id
  @param fi - first index
  @param ncomp- number of components
  @param sig - %vector containing part of stress components
   
  @return The function does not return anything.

  Created by JK, 19.5.2002
*/
void node::storestress (long lcid,long fi,long ncomp,vector &sig)
{
  long i,j,m;
  m=lcid*ncompstr;
  
  if (fi==0)  ncontr_stress[lcid]++;
  
  j=fi;
  for (i=m+fi;i<m+fi+ncomp;i++){
    stress[i]+=sig[j];
    j++;
  }
}



/**
  Function stores other components.
   
  @param fi - first index
  @param ncomp- number of components
  @param otherv - array containing part of other components
   
  @return The function does not return anything.

  Created by JK, 19.5.2002
*/
void node::storeother (long fi,long ncomp,vector &otherv)
{
  long i;
  
  if (fi==0)  ncontr_other++;
  
  for (i=fi;i<fi+ncomp;i++){
    other[i]+=otherv[i];
  }
}



/**
  Function stores other components.
   
  @param fi - first index
  @param ncomp - number of components
  @param vol - volume used for contribution
  @param otherv - array containing part of other components
   
  @return The function does not return anything.

  Created by JK, 19.5.2002
*/
void node::storeother (long fi,long ncomp,double vol,vector &otherv)
{
  long i;
  
  if (fi==0)  vol_other+=vol;
  
  for (i=fi;i<fi+ncomp;i++){
    other[i]+=otherv[i];
  }
}



/**
  Function averages strain components.
   
  @param lcid - load case id
   
  @return The function does not return anything.

  Created by JK,
*/
void node::strain_averageval (long lcid)
{
  long i;
  
  if (Mp->strainaver==1){
    if (ncontr_strain[lcid]>0){
      for (i=0;i<ncompstr;i++){
	strain[lcid*ncompstr+i]/=ncontr_strain[lcid];
      }
    }
  }
  if (Mp->strainaver==2){
    if (vol_strain[lcid]>0){
      for (i=0;i<ncompstr;i++){
	strain[lcid*ncompstr+i]/=vol_strain[lcid];
      }
    }
  }
}



/**
  Function averages stress components.
   
  @param lcid - load case id
   
  @return The function does not return anything.

  Created by JK,
*/
void node::stress_averageval (long lcid)
{
  long i;

  if (Mp->stressaver==1){
    if (ncontr_stress[lcid]>0){
      for (i=0;i<ncompstr;i++){
	stress[lcid*ncompstr+i]/=ncontr_stress[lcid];
      }
    }
  }
  if (Mp->stressaver==2){
    if (vol_stress[lcid]>0){
      for (i=0;i<ncompstr;i++){
	stress[lcid*ncompstr+i]/=vol_stress[lcid];
      }
    }
  }
}



/**
  Function averages components of array other.
   
   
  @return The function does not return anything.

  Created by JK,
*/
void node::other_averageval ()
{
  long i;
  
  if (Mp->otheraver==1){
    if (ncontr_other>0){
      for (i=0;i<ncompother;i++){
	other[i]/=ncontr_other;
      }
    }
  }
  if (Mp->otheraver==2){
    if (vol_other>0){
      for (i=0;i<ncompother;i++){
	other[i]/=vol_other;
      }
    }
  }
  
}



/**
  Function sets strain components to zero
   
  @param lcid - load case id
   
  @return The function does not return anything.

  Created by JK,
*/
void node::nullstrain (long lcid)
{
  long i;
  for (i=0;i<ncompstr;i++){
    strain[lcid*ncompstr+i]=0.0;
  }
  if (Mp->strainaver==1)
    ncontr_strain[lcid]=0;
  if (Mp->strainaver==2)
    vol_strain[lcid]=0.0;
}



/**
  Function sets stress components to zero.
   
  @param lcid - load case id
   
  @return The function does not return anything.

  Created by JK,
*/
void node::nullstress (long lcid)
{
  long i;
  for (i=0;i<ncompstr;i++){
    stress[lcid*ncompstr+i]=0.0;
  }
  if (Mp->stressaver==1)
    ncontr_stress[lcid]=0;
  if (Mp->stressaver==2)
    vol_stress[lcid]=0.0;
}



/**
  Function sets components of array other to zero.
   
  @return The function does not return anything.

  Created by JK,
*/
void node::nullother ()
{
  long i;
  for (i=0;i<ncompother;i++){
    other[i]=0.0;
  }
  if (Mp->otheraver==1)
    ncontr_other=0;
  if (Mp->otheraver==2)
    vol_other=0.0;
}



/**
   The function returns compid-th component of other state variable array at the given node.

   @param compid - id of other array component [input]

   @return  The function returns compid-th component of other array.
   
   Created by Tomas Koudelka, 09.2018
*/
double node::giveother (long compid)
{
  if (compid < ncompother)
    return other[compid];
  else
    print_err("other component id %ld is out of range <1;%ld>", __FILE__, __LINE__, __func__, compid+1, ncompother);

  return 0.0;
}



/**
  The function allocates array of DOF meaning.

  @param nid - node id

  @return The function does not return anything.

  Created by JK,
*/
void node::alloc_meaning (long nid)
{
  long i,ndofn;
  
  ndofn = Mt->give_ndofn (nid);
  meaning = new long [ndofn];
  for (i=0;i<ndofn;i++){
    meaning[i]=-1;
  }
}



/**
  Function cleans all arrays defined at node.
   
  @param nlc - number of load cases
   
  @return The function does not return anything.

  Created by JK, 23.8.2005
*/
void node::clean (long nlc)
{
  long i;
  
  for (i=0;i<nlc;i++){
    nullstrain (i);
    nullstress (i);
    nullother ();    
  }
}



/**
  Function allocates arrays used for problems with
  changing nodes, elements and DOFs
   
  @param nid - node id
  
  Created by JK, 7.11.2006
*/
void node::alloc_growstr (long nid)
{
  long ndofn;
  
  ndofn = Mt->give_ndofn (nid);
  
  nodval = new double [ndofn];
  memset(nodval, 0, sizeof(*nodval)*ndofn);
}


void node::give_transfmat(long ndofn, matrix &tmat)
{
  tmat[0][0]=e1[0];
  tmat[1][0]=e1[1];
  
  tmat[0][1]=e2[0];
  tmat[1][1]=e2[1];
  if (transf == 3){
    tmat[0][2]=e3[0];
    tmat[1][2]=e3[1];

    tmat[2][0]=e1[2];
    tmat[2][1]=e2[2];
    tmat[2][2]=e3[2];      
  }
  else{
    if (ndofn == 3) // 2D beam
      tmat[2][2]=1.0;
  }
}
