#include "interfacequadrilat.h"
#include "globalt.h"
#include "globmatt.h"

interfacequadrilat::interfacequadrilat (void)
{
  long i;
  
  //  number of nodes on element
  nne=4;
  //  geometric problem dimension (2D)
  ncomp=2;
  
  //  number of transported variables
  ntm=Tp->ntm;

  dofe = new long* [ntm];
  nip = new long* [ntm];
  ordering = new long* [ntm];
  for (i=0;i<ntm;i++){
    dofe[i] = new long [ntm];
    nip[i] = new long [ntm];
    ordering[i] = new long [nne];
  }
    
  switch (Tp->tmatt){
  case nomedium:{
    break;
  }
  case onemedium:{
    ordering[0][0]=1;  ordering[0][1]=2;  ordering[0][2]=3;  ordering[0][3]=4;
    dofe[0][0]=4;
    nip[0][0]=2;
    ndofe=4;
    break;
  }
  case twomediacoup:{
    ordering[0][0]=1;  ordering[0][1]=3;  ordering[0][2]=5;  ordering[0][3]=7;
    ordering[1][0]=2;  ordering[1][1]=4;  ordering[1][2]=6;  ordering[1][3]=8;
    
    
    if (Tp->savemode==0){
      nip[0][0]=2;       nip[0][1]=2;       nip[1][0]=2;       nip[1][1]=2;
    }
    if (Tp->savemode==1){
      nip[0][0]=2;       nip[0][1]=0;       nip[1][0]=0;       nip[1][1]=0;
    }
    
    
    dofe[0][0]=4;  dofe[0][1]=4;  dofe[1][0]=4;  dofe[1][1]=4;
    ndofe=8;
    break;
  }
  case threemediacoup:{
    ordering[0][0]=1;   ordering[0][1] =4;  ordering[0][2] =7;  ordering[0][3] =10;
    ordering[1][0]=2;   ordering[1][1] =5;  ordering[1][2] =8;  ordering[1][3] =11;
    ordering[2][0]=3;   ordering[2][1] =6;  ordering[2][2] =9;  ordering[2][3] =12;

    if (Tp->savemode==0){
      nip[0][0]=2;  nip[0][1]=2;  nip[0][2]=2;
      nip[1][0]=2;  nip[1][1]=2;  nip[1][2]=2;
      nip[2][0]=2;  nip[2][1]=2;  nip[2][2]=2;
    }
    if (Tp->savemode==1){
      nip[0][0]=2;  nip[0][1]=0;  nip[0][2]=0;
      nip[1][0]=0;  nip[1][1]=0;  nip[1][2]=0;
      nip[2][0]=0;  nip[2][1]=0;  nip[2][2]=0;
    }
    
    dofe[0][0]=4;  dofe[0][1]=4;  dofe[0][2]=4;
    dofe[1][0]=4;  dofe[1][1]=4;  dofe[1][2]=4;
    dofe[2][0]=4;  dofe[2][1]=4;  dofe[2][2]=4;


    ndofe=12;
    break;
  }
  case fourmediacoup:{
    ordering[0][0]=1;   ordering[0][1] =5;  ordering[0][2] =9;   ordering[0][3] =13;
    ordering[1][0]=2;   ordering[1][1] =6;  ordering[1][2] =10;  ordering[1][3] =14;
    ordering[2][0]=3;   ordering[2][1] =7;  ordering[2][2] =11;  ordering[2][3] =15;
    ordering[3][0]=4;   ordering[3][1] =8;  ordering[3][2] =12;  ordering[3][3] =16;

    if (Tp->savemode==0){
      nip[0][0]=2;  nip[0][1]=2;  nip[0][2]=2;  nip[0][3]=2;
      nip[1][0]=2;  nip[1][1]=2;  nip[1][2]=2;  nip[1][3]=2;
      nip[2][0]=2;  nip[2][1]=2;  nip[2][2]=2;  nip[2][3]=2;
      nip[3][0]=2;  nip[3][1]=2;  nip[3][2]=2;  nip[3][3]=2;
    }
    if (Tp->savemode==1){
      nip[0][0]=2;  nip[0][1]=0;  nip[0][2]=0;  nip[0][3]=0;
      nip[1][0]=0;  nip[1][1]=0;  nip[1][2]=0;  nip[1][3]=0;
      nip[2][0]=0;  nip[2][1]=0;  nip[2][2]=0;  nip[2][3]=0;
      nip[3][0]=0;  nip[3][1]=0;  nip[3][2]=0;  nip[3][3]=0;
    }
    
    dofe[0][0]=4;  dofe[0][1]=4;  dofe[0][2]=4;  dofe[0][3]=4;
    dofe[1][0]=4;  dofe[1][1]=4;  dofe[1][2]=4;  dofe[1][3]=4;
    dofe[2][0]=4;  dofe[2][1]=4;  dofe[2][2]=4;  dofe[2][3]=4;
    dofe[3][0]=4;  dofe[3][1]=4;  dofe[3][2]=4;  dofe[3][3]=4;

    ndofe=16;
    break;
  }
  default:{
    print_err("unknown number of transported matters is required",__FILE__,__LINE__,__func__);
  }
  }
}

interfacequadrilat::~interfacequadrilat (void)
{
  long i;

  for (i=0;i<ntm;i++){
    delete [] dofe[i];
    delete [] nip[i];
    delete [] ordering[i];
  }
  delete [] dofe;
  delete [] nip;
  delete [] ordering;
}

/**
   function assembles element code numbers
   they are used for localization between all element unknowns and unknowns related to one matter
   
   @param cn - code numbers
   @param ri - number of matter (usually row index)
*/
void interfacequadrilat::codnum (long *cn,long ri)
{
  long i;
  for (i=0;i<nne;i++){
    cn[i]=ordering[ri][i];
  }
}

/**
   function assembles gradient of %matrix of base functions

   @param gm - gradient %matrix
   @param x,y - array containing node coordinates
   @param l - the length of element (output)
   @param h - the width of element (input)
   
   JK, 4.1.2016
*/
void interfacequadrilat::grad_matrix (matrix &gm,vector &x,vector &y,double &l,double h)
{
  double sx,sy;
  
  //  direction vector
  sx = x[0]-x[1];
  sy = y[0]-y[1];
  
  //  lenght of the element
  l = sqrt(sx*sx+sy*sy);
  
  fillm (0.0,gm);
  
  gm[0][0]= 0.5/l;
  gm[0][1]=-0.5/l;
  gm[0][2]=-0.5/l;
  gm[0][3]= 0.5/l;

  gm[1][0]= 0.5/h;
  gm[1][1]= 0.5/h;
  gm[1][2]=-0.5/h;
  gm[1][3]=-0.5/h;
}

/**
   function computes conductivity %matrix of 2D problems for one transported matter
   
   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param km - conductivity %matrix

   JK, 4. 1.2016
*/
void interfacequadrilat::conductivity_matrix (long /*lcid*/,long eid,long ri,long ci,matrix &km)
{
  long ii;
  double thick,l,h,jac;
  ivector nodes(nne);
  vector x(nne),y(nne),t(nne);
  matrix gm(ncomp,dofe[ri][ci]),d;

  matrix n(1,dofe[ri][ci]);
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  thickness
  Tc->give_thickness (eid,nodes,t);
  //  node coordinates
  Tt->give_node_coord2d (x,y,eid);
  
  fillm (0.0,km);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];
  
  //  virtual width of the element
  h = Tm->ifacemat[Tm->ip[ii].idm].h;

  //  matrix of gradients
  grad_matrix (gm,x,y,l,h);
  
  //  matrix of conductivity of the material
  reallocm(ncomp,ncomp,d);
  Tm->matcond (d,ii,ri,ci);
  
  //  thickness
  thick = (t[0]+t[1]+t[2]+t[3])/4.0;
  
  if (l<0.0){
    print_err("\n negative lenght of interface quadrilateral element %ld", __FILE__, __LINE__, __func__,eid);
  }
  
  jac=thick*l*h;
  
  //  contribution to the conductivity matrix of the element
  bdbj (km.a,gm.a,d.a,jac,gm.m,gm.n);
  
}

/**
   function assembles resulting element conductivity %matrix

   @param eid - element id
   @param lcid - load case id
   @param km - resulting conductivity %matrix of one element

   JK, 4. 1. 2016
*/
void interfacequadrilat::res_conductivity_matrix (long eid,long /*lcid*/,matrix &km)
{
  long i,j,*rcn,*ccn;
  matrix lkm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      rcn = new long [dofe[i][j]];
      ccn = new long [dofe[i][j]];
      reallocm (dofe[i][j],dofe[i][j],lkm);
      conductivity_matrix (i,eid,i,j,lkm);
      codnum (rcn,i);
      codnum (ccn,j);
      mat_localize (km,lkm,rcn,ccn);
      delete [] rcn;  delete [] ccn;
    }
  }
}

/**
   function assembles resulting element capacity %matrix
   
   @param eid - element id
   @param cm - resulting capacity %matrix of one element

   JK, 4. 1. 2016
*/
void interfacequadrilat::res_capacity_matrix (long /*eid*/,matrix &/*cm*/)
{
  //  interface element does not have the capacity matrix
  //  there is no contribution to the resulting capacity matrix
}



/**
   Function returns transport (non-mechanical) quantities at nodes of element.
   The values of selected quantity is copied from the closest integration points 
   to element nodes.

   @param eid - element id
   @param nodval - %vector of nodal values
   @param ntq - type of non-mechanical quantity
   
   @return The function does not return anything.
   
   Created by TKo, 13.1.2016
*/
void interfacequadrilat::transq_nodval (long eid,vector &nodval,nonmechquant nmq)
{
  long i,ipid;
  ivector ipnum(nne);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipid=Tt->elements[eid].ipp[0][0];

  ipnum(0) = ipid;
  ipnum(1) = ipid;
  ipnum(2) = ipid+1;
  ipnum(3) = ipid+1;

  for (i=0;i<nne;i++)
  {
    //copy transport (non-mechanical) quantity from closest int. point
    nodval[i] = Tm->givetransq(nmq, ipnum[i]);
  }
}



/**
   Function returns initial values of transport (non-mechanical) quantities at nodes of element.
   The values of selected quantity is copied from the closest integration points 
   to element nodes.

   @param eid - element id
   @param nodval - %vector of nodal values
   @param ntq - type of non-mechanical quantity
   
   @return The function does not return anything.
   
   12/06/2012 TKr
   Modified by TKo, 10.10.2013
*/
void interfacequadrilat::transq_init_nodval (long eid,vector &nodval,nonmechquant nmq)
{
  long i,ipid;
  ivector ipnum(ASTCKIVEC(nne)), enod(ASTCKIVEC(nne));
  vector ipav;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipid=Tt->elements[eid].ipp[0][0];
  // element nodes
  Tt->give_elemnodes(eid, enod);

  ipnum(0) = ipid;
  ipnum(1) = ipid;
  ipnum(2) = ipid+1;
  ipnum(3) = ipid+1;

  for (i=0;i<nne;i++)
  {
    // create reference vector to the int. point av array
    makerefv(ipav, Tm->ip[ipnum[i]].av, Tp->ntm);
    // store initial nodal values to the integration point
    initnodval2 (enod[i], ipav);
    //copy transport (non-mechanical) quantity from closest int. point
    nodval[i] = Tm->givetransq(nmq, ipnum[i]);
  }
}



