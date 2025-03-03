#include "elementt.h"
#include "globalt.h"
#include <string.h>

long elementt::ntm = 0;

elementt::elementt (void)
{
  long i;
  
  //  type of element
  te = (elemtypet) 0;
  
  //  source connected with element
  //  this means no source
  source=0;
  
  //  array containing id of integration points
  ipp=NULL;
  //  type of cross section
  crst = (crsectypet) 0;
  //  cross section id
  idcs=0;
  
  //  number of transported matters
  ntm = Tp->ntm;
  
  //  types of material
  tm = NULL;
  //  id of materials
  idm = NULL;
  
  //  indicator of boundary conditions on element
  transi=new long [ntm];
  for (i=0;i<Tp->ntm;i++){
    //  this represents element with no boundary condition
    transi[i]=0;
  }
  
  //  array of initial nodal values
  //  it is used in problems with changing number of elements
  initnodval=NULL;

  // transformation matrix which transforms master node values to hanging node values
  tmat = NULL;
}

elementt::~elementt (void)
{
  long i;  
  if (ipp != NULL)
  {
    for (i=0; i<ntm; i++)
      delete [] ipp[i];
  }
  delete [] tm;
  delete [] idm;
  delete [] ipp;
  delete [] transi;
  delete [] initnodval;
  delete tmat;
}

/**
   function reads input data
   
   @param in - input file
   @param eid - element id
   
   JK, 25.11.2008
*/
void elementt::read (XFILE *in,long eid)
{
  long nne,ndofe;
  gelemtype get;
  
  //  elemet type
  xfscanf (in,"%d",(int*)&te);
  
  switch (te){
    case barlint:{
      get=linbar;
      if (Lbt==NULL)
        Lbt = new linbart;
      if (Tp->gdim < 1)
        Tp->gdim=1;
      break;
    }
    case barlint3d:{
      get=linbar;
      if (Lbt3d==NULL)
        Lbt3d = new linbart3d;
      if (Tp->gdim < 1)
        Tp->gdim=1;
      break;
    }
    case barlintax:{
      get=linbar;
      if (Lbt==NULL)
        Lbat = new linbartax;
      if (Tp->gdim < 1)
        Tp->gdim=1;
      break;
    }
    case barquadt:{
      get=quadbar;
      if (Qbt==NULL)
        Qbt = new quadbart;
      if (Tp->gdim < 1)
        Tp->gdim=1;
      break;
    }
    case barquadtax:{
      get=quadbar;
      if (Qbat==NULL)
        Qbat = new quadbartax;
      if (Tp->gdim < 1)
        Tp->gdim=1;
      break;
    }
    case trlint:{
      get=lintriag;
      if (Ltt==NULL)
        Ltt = new trlineart;
      if (Tp->gdim < 2)
        Tp->gdim=2;
      break;
    }
    case trlaxisym:{
      get=lintriag;
      if (Ltat==NULL)
        Ltat = new trlinaxisym;
      if (Tp->gdim < 2)
        Tp->gdim=2;
      break;
    }
    case quadlint:{
      get=linquad;
      if (Lqt==NULL)
        Lqt = new quadlineart;
      if (Tp->gdim < 2)
        Tp->gdim=2;
      break;
    }
    case quadquadt:{
      get=quadquad;
      if (Qqt==NULL)
        Qqt = new quadquadrilatt;
      if (Tp->gdim < 2)
        Tp->gdim=2;
      break;
    }
    case quadquadtax:{
      get=quadquad;
      if (Qqat==NULL)
        Qqat = new quadquadrilattax;
      if (Tp->gdim < 2)
        Tp->gdim=2;
      break;
    }
    case quadlaxisym:{
      get=linquad;
      if (Lqat==NULL)
        Lqat = new quadlinaxisym;
      if (Tp->gdim < 2)
        Tp->gdim=2;
      break;
    }
    case ifacequadel:{
      get=linquad;
      if (Ifcquadt==NULL)
        Ifcquadt = new interfacequadrilat;
      if (Tp->gdim < 2)
        Tp->gdim=2;
      break;
    }
    case lineartett:{
      get=lintetra;
      if (Ltett==NULL)
        Ltett = new lintett;
      if (Tp->gdim < 3)
        Tp->gdim=3;
      break;
    }
    case linearhext:{
      get=linhexa;
      if (Lht==NULL)
        Lht = new linhext;
      if (Tp->gdim < 3)
        Tp->gdim=3;
      break;
    }
    case quadratichext:{
      get=quadhexa;
      if (Qht==NULL)
        Qht = new quadhext;
      if (Tp->gdim < 3)
      Tp->gdim=3;
      break;
    }
    case gen2del:{
      get=gen2d;

      if (G2d==NULL)
        G2d = new gen2delem;
      if (Tp->gdim < 2)
        Tp->gdim=2;
      break;
    }
    case linearwedget:{
      get=linwed;
      
      if (Lwt==NULL)
        Lwt = new linwedget;
      if (Tp->gdim < 3)
        Tp->gdim=3;
      break;
    }
      
    default:{
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
    }
  }
  
  //  number of DOFs on element
  ndofe=Tt->give_ndofe (eid);
  //  number of nodes on element
  nne=Tt->give_nne(eid);
  
  if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
    Gtt->gelements[eid].read_gf (in,nne,ndofe,get,Gtt->nn);
    initnodval = new double [ndofe];
  }
  else{
    Gtt->gelements[eid].read (in,nne,ndofe,get,Gtt->nn);
  }
  
  //  type of cross section
  xfscanf (in,"%d",(int*)&crst);
  if (crst!=0){
    //  cross section id
    xfscanf (in,"%ld",&idcs);
    idcs--;
  }
  
  //  reading of material types and ids
  readmat (Tp->ntm*Tp->ntm,in);
}

/**
   function prints input data
   
   @param out - output file
   @param eid - element id
   
   JK, 25.11.2008
*/
void elementt::print (FILE *out,long eid)
{
  long nne,ndofe;
  gelemtype get;

  //  element type
  fprintf (out,"\n%d",(int)te);
  
  //  number of DOFs on element
  ndofe=Tt->give_ndofe (eid);
  //  number of nodes on element
  nne=Tt->give_nne(eid);
  //  type of element (see galias.h)
  get = Gtt->give_elem_type (eid);
  
  if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
    Gtt->gelements[eid].print_gf (out,nne,ndofe);
  }
  else{
    Gtt->gelements[eid].print (out,nne,ndofe,get);
  }

  fprintf (out," %d",crst);
  if (crst!=0){
    fprintf (out," %ld",idcs+1);
  }
  
  //  printing of material tpyes and ids
  printmat (Tp->ntm*Tp->ntm,out);  
}

/**
   function reads types of materials on elements
   
   @param m - number of types of material models
   @param in - input stream
   
*/
void elementt::readmat (long m,XFILE *in)
{
  long i;
  
  tm = new mattypet [m];
  idm = new long [m];
  
  for (i=0;i<m;i++){
    xfscanf (in,"%d %ld",(int*)&tm[i],&idm[i]);
    idm[i]--;
  }
}

/**
   function prints types of materials on elements
   
   @param m - number of types of material models
   @param in - input stream
   
*/
void elementt::printmat (long m,FILE *out)
{
  long i;
  
  for (i=0;i<m;i++){
    fprintf (out," %d %ld",tm[i],idm[i]+1);
  }
  fprintf (out,"\n");
}

/**
   function allocates array for initial nodal values
   initial nodal values are used in problems with growing number of elements
   
   @praram ndofe - number of DOFs on element
   
   JK, 3.3.2006
*/
void elementt::alloc_initnodval (long ndofe)
{
  initnodval = new double [ndofe];
}

/**
   function defines initial nodal values
   initial nodal values are used in problems with growing number of elements
   
   @param r[in] - %vector containing initial nodal values
   
   JK, 3.3.2006
*/
void elementt::initnodvalues (vector &r)
{
  long i;
  
  for (i=0;i<r.n;i++){
    initnodval[i]=r[i];
  }
  
}

/**
   function subtracts initial nodal values from total values
   function is used in problems with changing number of elements
   
   @param r - total nodal values
   @param ndofe - number of DOFs on element
   
   5.3.2006, JK
*/
void elementt::subtrinitnodval (double *r,long ndofe)
{
  long i;
  
  for (i=0;i<ndofe;i++){
    r[i]-=initnodval[i];
  }
}
