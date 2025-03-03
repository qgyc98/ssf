#include <stdlib.h>
#include "globmatt.h"
#include "globalt.h"
#include "lhsrhst.h"
#include "constrel.h"
#include "onemedium.h"
#include "twomedia.h"
#include "threemedia.h"
/**
   function assembles conductivity %matrix of the problem
   
   @param lcid - load case id
*/
void conductivity_matrix (long lcid)
{
  if (Kmat==NULL)  Kmat = new gmatrix;
  Kmat->ts=(storagetype)(*(int*)&Tp->tstorkm);
  Kmat->initiate (Gtt,Ndoft,Mesprt);
  //gmat_initialize (Tp,Kmat);
  Kmat->setval (Tp->ssle);

  long i,ndofe,*cn;
  matrix lm;
  
  for (i=0;i<Tt->ne;i++){
    ndofe = Gtt->give_ndofe (i);
    
    allocm (ndofe,ndofe,lm);
    conductmat (i,lcid,lm);

    cn = new long [ndofe];
    Gtt->give_code_numbers (i,cn);
    
    Kmat->localize (lm,cn,i);

    destrm (lm);
    delete [] cn;
  }
  
  Kmat->prepmat (Gtt,0.0,Mesprt);

}

/**
   function assembles capacity %matrix of the problem

   @param lcid - load case id
*/
void capacity_matrix (long lcid)
{
  if (Cmat==NULL)  Cmat = new gmatrix;
  Cmat->ts=(storagetype)(*(int*)&Tp->tstorcm);
  Cmat->initiate (Gtt,Ndoft,Mesprt);
  //gmat_initialize (Tp,Cmat);
  Cmat->setval (Tp->ssle);

  long i,ndofe,*cn;
  matrix lm;

  for (i=0;i<Tt->ne;i++){
    ndofe = Gtt->give_ndofe (i);

    allocm (ndofe,ndofe,lm);
    //charmat (i,lcid,m);
    capacmat (i,lcid,lm);

    cn = new long [ndofe];
    Gtt->give_code_numbers (i,cn);

    Cmat->localize (lm,cn,i);

    destrm (lm);
    delete [] cn;
  }

  Cmat->prepmat (Gtt,0.0,Mesprt);

}

/**
   @param r - %vector of residuum
   @param p - %vector d_{n+1} - d_{n}
   @param v - auxiliary %vector
   @param dt - actual time increment
   @param lcid - load case id
   
   JK, 2.2.2003
*/
void residuum (double *r,double *p,double *v,double dt,long n,long lcid)
{
  capacity_matrix (0);
  conductivity_matrix (0);
  trfel_right_hand_side (0,Lsrst->rhs,Ndoft);
  
  Cmat->gmxv (p,v);
  cmulv (1.0/dt,v,n);
  Kmat->gmxv (Lsrst->lhs,r);
  addv (r,v,n);
  subv (r,Lsrst->rhs,n);
}


/**
   function extracts values on one element

   @param lcid - number of load case
   @param r - vector of nodal values
   @param cn - array containing code numbers
   @param ndofe - number of DOFs on actual element

   9.7.2001
*/
void nodalvalues (long lcid,double *r,long *cn,long ndofe)
{
  long i,ii;
  
  if ((Tp->tprob==nonstationary_problem) || (Tp->tprob==nonlinear_nonstationary_problem)){
    for (i=0;i<ndofe;i++){
      ii=cn[i];
      if (ii<0)   r[i]=Tb->lc[0].pv[0-ii-1].getval();
      if (ii==0)  r[i]=0.0;
      if (ii>0){   
	r[i]=Lsrst->lhsi[ii-1]+Lsrst->lhs[ii-1];
	
	/*
	// ************************************************
	if (Tp->tmatt == threemediacoup)//bude predelano a presunuto jinam - pouze pro model z Padovy??!!
	  {
	    //tohle je jenom pro tri media - bude predelano
	    if (i == 0)//capillary pressure check //pokus??!!
	      if (r[0] <= 100.0){
		Lsrst->lhs[ii-1] = 100.0 - Lsrst->lhsi[ii-1];
		fprintf (Outt,"\n\n Uprava reseni strany pro kapilarni tlak");
	      }
	    
	    if (i == 1)//gas pressure check
	      if (r[1] <= 100.0){
		Lsrst->lhs[ii-1] = 100.0 - Lsrst->lhsi[ii-1];
		fprintf (Outt,"\n\n Uprava reseni strany pro tlak plynu");
	      }
	  }
	// ************************************************
	*/
      }
    }
  }
  else{
    for (i=0;i<ndofe;i++){
      ii=cn[i];
      if (ii<0)   r[i]=Tb->lc[lcid].pv[0-ii-1].getval();
      if (ii==0)  r[i]=0.0;
      if (ii>0)   r[i]=Lsrst->lhs[ii-1];
    }
  }
  /*
  // ************************************************
  if (Tp->tmatt == threemediacoup)//bude predelano a presunuto jinam - pouze pro model z Padovy??!!
    {
      double pgw;
      state_eq tt;
      
      //gas pressure check
      pgw = tt.get_pgw(r[0],r[2]);
      if (pgw > r[1]) 
	r[1] = pgw;
      
      //minimum gas pressure check
      if (r[1] <= 100.0)
	r[1] = 100.0;
      
      //minimum capillary pressure check
      if (r[0] <= 100.0)
	r[0] = 100.0;
    }
  // ************************************************
  */
}


/**
   function extracts values on one node
   
   @param lcid  - number of load case
   @param r     - allocated array for displacement
   @param idn   - node number

   9.7.2001
*/
void nodval (long lcid,double *r, long idn)
{
  long i,ii;
  
  if ((Tp->tprob==nonstationary_problem) || (Tp->tprob==nonlinear_nonstationary_problem)){
    for (i=0;i<Tt->give_ndofn (idn);i++){
      ii=Tt->give_dof(idn,i);
      if (ii<0)   r[i]=Tb->lc[0].pv[0-ii-1].getval();
      if (ii==0)  r[i]=0.0;
      if (ii>0)      if (ii>0){   
	r[i]=Lsrst->lhsi[ii-1]+Lsrst->lhs[ii-1];
	
	/*
	// ************************************************
	if (Tp->tmatt == threemediacoup)//bude predelano a presunuto jinam - pouze pro model z Padovy??!!
	  {
	    if (i == 0)//capillary pressure check //pokus??!!
	      if (r[0] <= 100.0){
		Lsrst->lhs[ii-1] = 100.0 - Lsrst->lhsi[ii-1];
		fprintf (Outt,"\n\n Uprava reseni strany pro kapilarni tlak");
	      }    
	    if (i == 1)//gas pressure check
	      if (r[1] <= 100.0){
		Lsrst->lhs[ii-1] = 100.0 - Lsrst->lhsi[ii-1];
		fprintf (Outt,"\n\n Uprava reseni strany pro tlak plynu");
	      }
	      }
	      // ************************************************
	      */
      }
    }
  }
  else{
    for (i=0;i<Tt->give_ndofn (idn);i++){
      ii=Tt->give_dof(idn,i);
      if (ii<0)   r[i]=Tb->lc[lcid].pv[0-ii-1].getval();
      if (ii==0)  r[i]=0.0;
      if (ii>0)   r[i]=Lsrst->lhs[ii-1];
    }
  }
  
  /*
  // ************************************************
  if (Tp->tmatt == threemediacoup)//bude predelano a presunuto jinam - pouze pro model z Padovy??!!
    {
      double pgw;
      state_eq tt;
      
      //gas pressure check
      pgw = tt.get_pgw(r[0],r[2]);
      if (pgw > r[1]) 
	r[1] = pgw;
      
      //minimum gas pressure check
      if (r[1] <= 100.0)
	r[1] = 100.0;
      
      //minimum capillary pressure check
      if (r[0] <= 100.0)
	r[0] = 100.0;
    }
  // ************************************************
  */
}


/**
   function extracts first time derivatives of nodal values on one element

   @param r - vector of derivatives of nodal values
   @param cn - array containing code numbers
   @param ndofe - number of DOFs on actual element

   4.3.2003
*/
void nodalderivatives (double *r,long *cn,long ndofe)
{
  long i,ii;
  
  if ((Tp->tprob==nonstationary_problem) || (Tp->tprob==nonlinear_nonstationary_problem)){
    for (i=0;i<ndofe;i++){
      ii=cn[i];
      if (ii<0)   r[i]=0.0;
      if (ii==0)  r[i]=0.0;
      if (ii>0)   r[i]=Lsrst->tdlhs[ii-1];
    }
  }
}

/**
   function extracts values on one element
   
   @param eid - element id
   @param r - vector of nodal values
   @param cn - array containing code numbers
   
   21.3.2004, JK
*/
void prescvalues (double *r,long *cn,long ndofe)
{
  long i,ii;
  
  for (i=0;i<ndofe;i++){
    ii=cn[i];
    if (ii<0)   r[i]=Tb->lc[0].pv[0-ii-1].getval();
    if (ii==0)  r[i]=0.0;
    if (ii>0)   r[i]=Lsrst->lhsi[ii-1];
  }
}

/**
   function extracts values on one element
   
   @param eid - element id
   @param r - vector of nodal values
   @param cn - array containing code numbers
   
   21.3.2004, JK
*/
void initialvalues (double *r,long *cn,long ndofe)
{
  long i,ii;
  
  for (i=0;i<ndofe;i++){
    ii=cn[i];
    if (ii<0)   r[i]=Tb->lc[0].pv[0-ii-1].getipv();
    if (ii==0)  r[i]=0.0;
    if (ii>0)   r[i]=Lsrst->lhsi[ii-1];
  }
}

/**
   function computes conductivity %matrix of one element
   
   @param eid - element id
   @param lcid - load case id
   @param km - conductivity %matrix of one element

*/
void conductmat (long eid,long lcid,matrix &km)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);
  
  switch (te){
  case barlint:{
    Lbt->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case barlintax:{
    Lbat->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case barquadt:{
    Qbt->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case barquadtax:{
    Qbat->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case trlint:{
    Ltt->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case trlaxisym:{
    Ltat->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case quadlint:{
    Lqt->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case quadquadt:{
    Qqt->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case quadquadtax:{
    Qqat->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case quadlaxisym:{
    Lqat->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case lineartett:{
    Ltett->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case linearhext:{
    Lht->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case quadratichext:{
    Qht->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown element type is required in function");
    fprintf (stderr,"\n conductmat (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
}

/**
   function computes capacity %matrix of one element

   @param eid - element id
   @param lcid - load case id
   @param km - capacity %matrix of one element

*/
void capacmat (long eid,long lcid,matrix &cm)
{
  elemtypet te;

  te=Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    Lbt->res_capacity_matrix (eid,cm);
    break;
  }
  case barlintax:{
    Lbat->res_capacity_matrix (eid,cm);
    break;
  }
  case barquadt:{
    Qbt->res_capacity_matrix (eid,cm);
    break;
  }
  case barquadtax:{
    Qbat->res_capacity_matrix (eid,cm);
    break;
  }
  case trlint:{
    Ltt->res_capacity_matrix (eid,cm);
    break;
  }
  case trlaxisym:{
    Ltat->res_capacity_matrix (eid,cm);
    break;
  }
  case quadlint:{
    Lqt->res_capacity_matrix (eid,cm);
    break;
  }
  case quadquadt:{
    Qqt->res_capacity_matrix (eid,cm);
    break;
  }
  case quadquadtax:{
    Qqat->res_capacity_matrix (eid,cm);
    break;
  }
  case quadlaxisym:{
    Lqat->res_capacity_matrix (eid,cm);
    break;
  }
  case lineartett:{
    Ltett->res_capacity_matrix (eid,cm);
    break;
  }
  case linearhext:{
    Lht->res_capacity_matrix (eid,cm);
    break;
  }
  case quadratichext:{
    Qht->res_capacity_matrix (eid,cm);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown element type is required in function");
    fprintf (stderr,"\n capacmat (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


/**
   function assembles %vector of nodal internal fluxes

   @param lcid - load case id
   @param intflux - array containing internal nodal fluxes
   
*/
void internal_fluxes (double *intflux,long n)
{
  long i,ne,ndofe,*cn;
  elemtypet te;
  vector lif;
  
  ne = Tt->ne;
  nullv (intflux,n);

  for (i=0;i<ne;i++){
    te = Tt->give_elem_type (i);
    ndofe = Gtt->give_ndofe (i);
    allocv (ndofe,lif);
    fillv (0.0,lif);

    cn = new long [ndofe];
    Gtt->give_code_numbers (i,cn);
    
    switch (te){
    case barlint:{
      Lbt->res_internal_fluxes (i,lif);
      break;
    }
    case barlintax:{
      Lbat->res_internal_fluxes (i,lif);
      break;
    }
   case barquadt:{
      Qbt->res_internal_fluxes (i,lif);
      break;
    }
   case barquadtax:{
      Qbat->res_internal_fluxes (i,lif);
      break;
    }
    case trlint:{
      Ltt->res_internal_fluxes (i,lif);
      break;
    }
    case trlaxisym:{
      Ltat->res_internal_fluxes (i,lif);
      break;
    }
    case quadlint:{
      Lqt->res_internal_fluxes (i,lif);
      break;
    }
    case quadquadt:{
      Qqt->res_internal_fluxes (i,lif);
      break;
    }
    case quadquadtax:{
      Qqat->res_internal_fluxes (i,lif);
      break;
    }
    case quadlaxisym:{
      Lqat->res_internal_fluxes (i,lif);
      break;
    }
    case lineartett:{
      Ltett->res_internal_fluxes (i,lif);
      break;
    }
    case linearhext:{
      Lht->res_internal_fluxes (i,lif);
      break;
    }
    case quadratichext:{
      Qht->res_internal_fluxes (i,lif);
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown element type is required in function internal_fluxes (file %s, line %d).\n",__FILE__,__LINE__);
      abort ();
    }
    }
    
    locglob (intflux,lif.a,cn,ndofe);
    destrv (lif);
    delete [] cn;
    
  }

  trfel_bound_flux (0,intflux,n);

}

/**
   function computes quantity values in integration points
   
   16.1.2002
*/
void approximation ()
{
  long i;
  
  switch (Tp->tmatt){
  case onemedium:{
    for (i=0;i<Tt->ne;i++){
      intpointvalues (i);
      intpointgradients (i);
    }
    break;
  }
  case twomediacoup:{
    for (i=0;i<Tt->ne;i++){
      intpointvalues (i);
      intpointgradients (i);
    }
    break;
  }
  case threemediacoup:{
    for (i=0;i<Tt->ne;i++){
      intpointvalues (i);
      intpointgradients (i);
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of element in function approximation (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }

}

void intpointvalues (long eid)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    Lbt->intpointval (eid);
    break;
  }
  case barlintax:{
    Lbat->intpointval (eid);
    break;
  }
  case barquadt:{
    Qbt->intpointval (eid);
    break;
  }
  case barquadtax:{
    Qbat->intpointval (eid);
    break;
  }
  case trlint:{
    Ltt->intpointval (eid);
    break;
  }
  case trlaxisym:{
    Ltat->intpointval (eid);
    break;
  }
  case quadlint:{
    Lqt->intpointval (eid);
    break;
  }
  case quadquadt:{
    Qqt->intpointval (eid);
    break;
  }
  case quadquadtax:{
    Qqat->intpointval (eid);
    break;
  }
  case quadlaxisym:{
    Lqat->intpointval (eid);
    break;
  }
  case lineartett:{
    Ltett->intpointval (eid);
    break;
  }
  case linearhext:{
    Lht->intpointval (eid);
    break;
  }
  case quadratichext:{
    Qht->intpointval (eid);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown element type is required in function intpointvalues (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
}

void intpointgradients (long eid)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    Lbt->intpointgrad (eid);
    break;
  }
  case barlintax:{
    Lbat->intpointgrad (eid);
    break;
  }
  case barquadt:{
    Qbt->intpointgrad (eid);
    break;
  }
  case barquadtax:{
    Qbat->intpointgrad (eid);
    break;
  }
  case trlint:{
    Ltt->intpointgrad (eid);
    break;
  }
  case trlaxisym:{
    Ltat->intpointgrad (eid);
    break;
  }
  case quadlint:{
    Lqt->intpointgrad (eid);
    break;
  }
  case quadquadt:{
    Qqt->intpointgrad (eid);
    break;
  }
  case quadquadtax:{
    Qqat->intpointgrad (eid);
    break;
  }
  case quadlaxisym:{
    Lqat->intpointgrad (eid);
    break;
  }
  case lineartett:{
    Ltett->intpointgrad (eid);
    break;
  }
  case linearhext:{
    Lht->intpointgrad (eid);
    break;
  }
  case quadratichext:{
    Qht->intpointgrad (eid);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown element type is required in function intpointvalues (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
}


/**
   function assembles part of the %vector of right hand side
   
   @param rhs - right hand side
   
   21.3.2004, JK
*/
void assemble_init (double *rhs)
{
  long i,ndofe,*cn,lcid=0,ii,kk;
  vector r,f;
  matrix km;
  long ne = Tt->ne;
  
  for (i=0;i<ne;i++){
    ndofe=Tt->give_ndofe(i);
    allocm (ndofe,ndofe,km);
    
    conductmat (i,lcid,km);
    
    allocv (ndofe,r);
    allocv (ndofe,f);
    cn = new long [ndofe];
    
    Tt->give_code_numbers (i,cn);
    
    prescvalues (r.a,cn,ndofe);
    //initialvalues (r.a,cn,ndofe);
    
    mxv (km,r,f);
    
    cmulv (-1.0,f);
    
    locglob (rhs,f.a,cn,ndofe);

    destrm (km);  destrv (f);  destrv (r);
    delete [] cn;
  }
}

/**
   function assembles right hand side of the problem
   
   @param lcid - load case id
   @param rhs - array containing right hand side
   
   12.3.2002
*/
void trfel_right_hand_side (long lcid,double *rhs,long n)
{
  double *av;
  av = new double [n];
  nullv (av,n);
  
  nullv (rhs,n);
  
  assemble_init (rhs);
  
  switch (Tp->tmatt){
  case onemedium:{
    Tb->lc[lcid].assemble (lcid,av+lcid*n);
    addv (rhs,av,n);
    break;
  }
  case twomediacoup:{
    Tb->lc[0].assemble (0,av);
    addv (rhs,av,n);
    nullv (av,n);
    Tb->lc[1].assemble (1,av);
    addv (rhs,av,n);
    break;
  }
  case threemediacoup:{
    Tb->lc[0].assemble (0,av);
    addv (rhs,av,n);
    nullv (av,n);
    Tb->lc[1].assemble (1,av);
    addv (rhs,av,n);
    nullv (av,n);
    Tb->lc[2].assemble (2,av);
    addv (rhs,av,n);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown transported matter is required in function");
    fprintf (stderr,"\n right_hand_side (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  delete [] av;
}


/**
   function assembles flux on the boundary
   
   @param lcid - load case id
   @param iflux - array containing flux (right hand side)
   
   12.3.2002
*/
void trfel_bound_flux (long lcid,double *iflux,long n)
{
  double *av;
  av = new double [n];
  
  switch (Tp->tmatt){
  case onemedium:{
    Tb->lc[lcid].assemble_flux (lcid,av+lcid*n,n);
    addv (iflux,av,n);
    break;
  }
  case twomediacoup:{
    Tb->lc[0].assemble_flux (0,av,n);
    addv (iflux,av,n);
    nullv (av,n);
    Tb->lc[1].assemble_flux (1,av,n);
    addv (iflux,av,n);
    break;
  }
  case threemediacoup:{
    Tb->lc[0].assemble_flux (0,av,n);
    addv (iflux,av,n);
    nullv (av,n);
    Tb->lc[1].assemble_flux (1,av,n);
    addv (iflux,av,n);
    nullv (av,n);
    Tb->lc[2].assemble_flux (2,av,n);
    addv (iflux,av,n);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown transported matter is required in function");
    fprintf (stderr,"\n right_hand_side (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  delete [] av;
}


/**
   function computes all required values
   this function is called e.g. before output print

   @param lcid - load case id

   1.4.2004, JK

   25.6.2004 changed by TKr
*/
void compute_req_valt (long lcid)
{
  //  gradients computation
  //if (Tp->gradcomp==1)
  //Tm->compute_nodegrads (lcid);
  //  fluxes computation
  //if (Tp->fluxcomp==1)
  //Tm->compute_nodefluxes (lcid);
  //  other values computation
  if (Tp->othercomp==1)
    Tm->compute_nodeothers (lcid);
}

/**
   function computes contributions from source of quantity
   
   @param lcid - load case id
   @param eid - element id
   @param nodval - array of source nodal values
   @param lv - local %vector
   
   JK, 25.6.2005
*/
void source_vector (long lcid,long eid,vector &nodval,vector &lv)
{
  elemtypet te;
  
  //  type of element
  te=Tt->elements[eid].te;
  
  switch (te){
  case barlint:{
    Lbt->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case barlintax:{
    Lbat->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case trlint:{
    Ltt->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case trlaxisym:{
    Ltat->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case quadlint:{
    Lqt->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case quadlaxisym:{
    Lqat->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case quadquadt:{
    Qqt->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case quadquadtax:{
    Qqat->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case lineartett:{
    Ltett->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case linearhext:{
    Lht->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case quadratichext:{
    Qht->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown element type is required in function source_vector (file %s, line %d",__FILE__,__LINE__);
  }
  }
  
}


void give_nodal_humid (double *gv,long mnt)
{
  long i,cn,compother;
  double *r;
  
  r = new double [Tp->ntm];
  
  switch (Tp->tmatt){
  case onemedium:{
    med1 m1;

    for (i=0;i<mnt;i++){
      cn=Tt->give_dof (i,0);
      if (cn>0)  gv[i] = Lsrst->lhsi[cn-1]+Lsrst->lhs[cn-1];
      if (cn<0)  gv[i] = Tb->lc[0].pv[0-cn-1].getval();
    }
    break;
    
    gv[i] = m1.compute_othervalues(1,0,r);
    if(gv[i] >= 1.0)
      gv[i] = 1.0;
  }
  case twomediacoup:{
    med2 m2;
    
    for (i=0;i<mnt;i++){
      cn=Tt->give_dof (i,0);
      if (cn>0)  r[0] = Lsrst->lhsi[cn-1]+Lsrst->lhs[cn-1];
      if (cn==0) r[0] = 0.0;
      if (cn<0)  r[0] = Tb->lc[0].pv[0-cn-1].getval();
      
      cn=Tt->give_dof (i,1);
      if (cn>0)  r[1] = Lsrst->lhsi[cn-1]+Lsrst->lhs[cn-1];
      if (cn==0) r[1] = 0.0;
      if (cn<0)  r[1] = Tb->lc[0].pv[0-cn-1].getval();
      
      compother = 1;
      if ((Tm->ip[0].tm) == kunzel)
	compother = 0;

      gv[i] = m2.compute_othervalues(compother,0,r);
      if(gv[i] >= 1.0)
	gv[i] = 1.0;
    }
    break;
  }
  case threemediacoup:{
    med3 m3;
    
    for (i=0;i<mnt;i++){
      cn=Tt->give_dof (i,0);
      if (cn>0)  r[0] = Lsrst->lhsi[cn-1]+Lsrst->lhs[cn-1];
      if (cn==0) r[0] = 0.0;
      if (cn<0)  r[0] = Tb->lc[0].pv[0-cn-1].getval();
      
      cn=Tt->give_dof (i,1);
      if (cn>0)  r[1] = Lsrst->lhsi[cn-1]+Lsrst->lhs[cn-1];
      if (cn==0) r[1] = 0.0;
      if (cn<0)  r[1] = Tb->lc[0].pv[0-cn-1].getval();
      
      cn=Tt->give_dof (i,2);
      if (cn>0)  r[2] = Lsrst->lhsi[cn-1]+Lsrst->lhs[cn-1];
      if (cn==0) r[2] = 0.0;
      if (cn<0)  r[2] = Tb->lc[0].pv[0-cn-1].getval();
      
      gv[i] = m3.compute_othervalues(3,0,r);
      if(gv[i] >= 1.0)
	gv[i] = 1.0;
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of transported matter in function (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  delete [] r;
}


/**
   function corrects solution of algebraic system of equations
   unacceptable trial values are replaced by the limit values
   
   14.7.2005
*/
void solution_correction ()
{
  long i,j,k,cn,ne,nne,ndofn,lcid,ipp;
  ivector nodes;
  vector nv;

  //  number of elements
  ne = Tt->ne;
  
  //  load case id, default value is 0
  lcid=0;
  
  //  loop over elements
  for (i=0;i<ne;i++){
    //  number nodes on element
    nne = Tt->give_nne (i);
    
    //  allocation
    allocv (nne,nodes);

    //  node numbers on element
    Tt->give_elemnodes (i,nodes);
    
    //  integration point id
    ipp=Tt->elements[i].ipp[0][0];
    
    //  extraction of nodal values    
    for (j=0;j<nne;j++){
      //  number of DOFs on node
      ndofn = Tt->give_ndofn (nodes[j]);
      //  allocation
      allocv (ndofn,nv);
      //  nodal values
      nodval (lcid,nv.a,nodes[j]);
      //  correction of values
      Tm->values_correction (nv,ipp);

      //  replacement of values in array containing the solution
      for (k=0;k<ndofn;k++){
	cn = Tt->give_dof (nodes[j],k);
	if (cn>0){
	  Lsrst->lhs[cn-1]=nv[k]-Lsrst->lhsi[cn-1];
	}
      }
      
      //
      destrv (nv);
    }
    
    destrv (nodes);
  }
  
}
