#include <string.h>
#include "gnode.h"
#include "gtopology.h"
#include "gfunct.h"



gnode::gnode (void)
{
  //  number of degrees of freedom
  ndofn=0;
  //  array of code numbers
  cn=NULL;
  //  node coordinates
  x=0.0;  y=0.0;  z=0.0;
  
  //  it is defined with respect to application of sparse direct solver
  //  see function auxdatsparsesolver in gmatrix.cpp
  ai=-1;
  
  //  general time function
  tgf=NULL;
  
  //  array of the master nodes (if the node is hanging node)
  mnodes = NULL;
  //  natural coordinates of the hanging node
  natcoord = NULL;
  //  type of master entity
  masentity = noel;
}



gnode::~gnode (void)
{
  //  array containing code numbers
  delete [] cn;
  //  numbers (pointers) of general functions
  delete [] tgf;
  //  array of the master nodes (if the node is hanging node)
  delete [] mnodes;
  //  natural coordinates of the hanging node
  delete [] natcoord;
}



/**
   functions reads basic informations
   
   @param in - input stream
   
   JK
*/
void gnode::read (XFILE *in)
{
  long i;
  //  coordinates
  xfscanf (in,"%lf %lf %lf",&x,&y,&z);
  
  //  code numbers
  xfscanf (in,"%ld",&ndofn);
  if (ndofn>0){
    //  in this context, ndofn is the number of DOFs
    //  defined in the node
    cn = new long [ndofn];
    for (i=0;i<ndofn;i++){
      cn[i]=1;
    }
  }
  if (ndofn<0){
    //  this is a hanging node
    //  ndofn is the number of the master nodes
    mnodes = new long [0-ndofn];
    //  the master nodes are read
    for (i=0;i<0-ndofn;i++){
      xfscanf (in,"%ld",mnodes+i);
      mnodes[i]--;
    }
    natcoord = new double [3];
    //  the natrual coordinates are read
    xfscanf (in,"%le %le %le",natcoord+0,natcoord+1,natcoord+2);
    //  the type of master entity is read
    xfscanf (in,"%d",&masentity);
  }
}



/**
   functions prints basic informations
   
   @param in - input stream
   
   JK
*/
void gnode::print (FILE *out)
{
  long i;
  //  coordinates
  fprintf (out,"  %f %f %f  ",x,y,z);
  
  //  code numbers
  fprintf (out," %ld ",ndofn);
  if (ndofn<0){
    for (i=0;i<0-ndofn;i++){
      fprintf (out," %ld ",mnodes[i]+1);
    }
    //  the natrual coordinates are read
    fprintf (out,"  %e %e %e ",natcoord[0],natcoord[1],natcoord[2]);
    //  the type of master entity is read
    fprintf (out," %d ",masentity);
  }
}



/**
   function reads constraints from input file
   
   @param dofcontr - type of dof control 
                     0 = non-growing construction, 
                     1 = dofs controlled by time function used in growing construction
   @param in - input stream
   
   TKr
*/
void gnode::constr (long dofcontr,XFILE *in)
{
  long i;

  if (ndofn < 0)
    print_err("Hanging node cannot be constrained", __FILE__, __LINE__, __func__);

  //  classical approach, constraints are read
  if (dofcontr==0){
    for (i=0;i<ndofn;i++){
      xfscanf (in,"%ld",cn+i);
    }
  }
  
  //  type of functions for DOFs control
  if (dofcontr==1){
    tgf = new long [ndofn];
    for (i=0;i<ndofn;i++){
      xfscanf (in,"%ld",tgf+i);
      tgf[i]--;
    }
  }
}



/**
   function prints constraints from input file
   
   @param in - input stream
   
   TKr
*/
void gnode::print_constr (long dofcontr,FILE *out)
{
  long i;
  
  if (ndofn < 0)
    print_err("Hanging node cannot be constrained", __FILE__, __LINE__, __func__);

  //  classical approach, constraints are read
  if (dofcontr==0){
    for (i=0;i<ndofn;i++){
      if(cn[i] > 0)
	fprintf (out," 1 ");//this is only for printing
      else
	fprintf (out," %ld ",cn[i]);
    }
  }
  
  //  type of functions for DOFs control
  if (dofcontr==1){
    for (i=0;i<ndofn;i++){
      fprintf (out,"%ld",tgf[i]+1);
    }
  }
}



/**
   function returns number of degrees of freedom of the node
   
   JK
*/
long gnode::give_ndofn ()
{
  return ndofn;
}



/**
   function returns number of required degree of freedom
   
   @param m - required component
   
   JK
*/
long gnode::give_dof (long m)
{
  if (m>=ndofn || m<0){
    print_err("wrong number of component is required,\n %ld-th dof is out of interval <%ld;%ld>",
              __FILE__,__LINE__,__func__, m+1, 1, ndofn);
  }
  return cn[m];
}



/**
   function locates number of DOF to required position
   
   @param m - number of component
   @param num - number of DOF which will be saved
   
   JK
*/
void gnode::save_dof (long m,long num)
{
  if (m>=ndofn || m<0){
    print_err("wrong number of component is required,\n %ld-th dof is out of interval <%ld;%ld>",
              __FILE__,__LINE__,__func__, m+1, 1, ndofn);
  }
  cn[m]=num;
}



/**
   The function clears array of DOF numbers.
   
   
   Created by TKo, 7.6.2013
*/
void gnode::clear_dof ()
{
  if (ndofn > 0)
    memset(cn, 0, sizeof(*cn)*ndofn);
}



/**
   The function clears j-th DOF number.
  
   @param j - index of DOF whose DOF number should be cleared (set to zero)
   
   @return The function does not return anything but it changes the content of cn array.

   Created by TKo, 07.2016
*/
void gnode::clear_dof (long j)
{
  if (j < ndofn)
    cn[j] = 0L;
  else
    print_err("invalid DOF index at node is required, j=%ld is out of range <0, %ld>", __FILE__, __LINE__, __func__, j, ndofn-1);
}



/**
   function computes distance of first point (=this node) and second point defined by 'c'
   
   @param dim - dimension of problem
   @param c - cartesian coordinates of second point
   
   created  8.6.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
double gnode::distance2 (long dim,const double *c)
{
  if      (dim == 1)  return ((x-c[0])*(x-c[0]));
  else if (dim == 2)  return ((x-c[0])*(x-c[0]) + (y-c[1])*(y-c[1]));
  else if (dim == 3)  return ((x-c[0])*(x-c[0]) + (y-c[1])*(y-c[1]) + (z-c[2])*(z-c[2]));
  else { print_err("Wrong dimension",__FILE__,__LINE__,__func__);  return (-1); }
  return 0;
}



/**
   function changes DOFs
   
   @param gf - array of time functions
   @param time - actual time
   @param lnso - state indicator of node according to state of adjacent elements
   
   6.2.2006, JK
   Updated 21.5.2013, TKo
*/
void gnode::update_dofs (gfunct *gf,double time, long lnso)
{
  long i,j,k;
  
  for (i=0;i<ndofn;i++){
    //  time function id
    j=tgf[i];
    if (j < 0) // this dof is controlled by time functions of the adjacent elements
    {
      cn[i] = lnso;
    }
    else
    {
      //  time function at actual time
      k=gf[j].getval_long (time);
      cn[i]=k;
    }
  }
}



/**
   The function searches for changed DOFs.
   
   @param gf - array of time functions
   @param time - actual time
   @param prev_time - previous time
   @param lnso - actual state indicator of the node
   @param plnso - state indicator of the node in the previous time
   
   Created 6.2.2006, JK
   Updated 21.5.2013,TKo
*/
long gnode::search_changed_dofs (gfunct *gf,double time,double prev_time, long lnso, long plnso)
{
  long i,j,k,l,nch;
  
  //  number of changes
  nch=0;
    
  if (time>prev_time){
    for (i=0;i<ndofn;i++){
      //  time function id
      j=tgf[i];

      if (j < 0) // this dof is controlled by time functions of the adjacent elements
      {
        if (plnso != lnso) // node was switched on/off and now it is switched off/on
          nch++;
      }
      else
      {
        //  time function at actual time
        k=gf[j].getval_long (time);
        //  time function at previous time
        l=gf[j].getval_long (prev_time);

/*
        if ((l <= 0) &&  // there was prescribed displacement or rigid support
            (k <= 0))    // and actually, there is prescribed displacement or rigid support again =>
          continue;      // => no change of number of dofs
       // BUT in the above case, the code numbers must be actualized in the case of change from the 
       // rigid support to the prescribed displacement and vice versa
*/
        
        if (k == l)     // there are identical conditions compared with the previous time
          continue;
        
        // excluding the above two conditions there must be change
        nch++;
      }
    }
  }
  else{
    //  initial time step
    //  time = starting_time
    for (i=0;i<ndofn;i++){
      //  time function id
      j=tgf[i];
      if (j < 0) // this dof is controlled by time functions of the adjacent elements
      {
        if (plnso != lnso) // node was switched on/off and now it is switched off/on
          nch++;
      }
      else
      {
        //  time function at actual time
        k=gf[j].getval_long (time);
        
        cn[i]=k;
        nch++;
      }
    }
  }
  
  return nch;
}
