#include "globmatc.h"
#include "globmat.h"
#include "globmatt.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechbclc.h"
#include "globalg.h"
#include "globalt.h"
#include "globalc.h"
#include "loadcase.h"
#include "element.h"
#include "elemswitch.h"
#include "elemswitcht.h"
#include "elemswitchc.h"
#include "intpoints.h"
#include "node.h"
#include "element.h"
#include "dloadcase.h"



/**
   function assembles %matrix which contains stiffness and conductivity matrices
   
   @param lcid - load case id

   JK
*/
void zero_order_matrix (long lcid)
{
  long i,mndofe,tndofe;
  ivector rcn,ccn;
  matrix lm;
  
  if (D0mat==NULL)  D0mat = new gmatrix;
    
  //  provizorni gtopologie (kombinace Gtm a Gtt)
  
  //D0mat->ts=Cp->tstord0;
  D0mat->setval (Cp->ssle);
  D0mat->initiate (Gtu,Ndofc,Cp->tstord0,Mesprc,Outc);

  
  //  mazani provizorni topologie
    
  //D0mat->setval (Cp->ssle.tlinsol,Cp->ssle.prec.pt,Cp->ssle.ni,Cp->ssle.err);
  
  //fprintf (Outc,"\n\n\n na zacatku\nMATICE VODIVOSTI (TUHOSTI) \n\n\n");
  //D0mat->printmat (Outc);
  
  //***************************************
  //  contributions from stiffness matrices
  //***************************************
  //  contributions from MEFEL part
  for (i=0;i<Mt->ne;i++){
    mndofe = Gtm->give_ndofe (i);

    reallocm (RSTCKMAT(mndofe,mndofe,lm));
    stiffmat (lcid, i, lm);
    
    reallocv (RSTCKIVEC(mndofe,rcn));
    //rcn = new long [mndofe];
    Gtm->give_code_numbers (i,rcn.a);
    
    D0mat->localize (lm,rcn,i);
  }
  

  //fprintf (Outc,"\n\n\n mefel\nMATICE TUHOSTI \n\n\n");
  //D0mat->printmat (Outc);

  //******************************************
  //  contributions from conductivity matrices
  //******************************************
  //  contributions from TRFEL part
  for (i=0;i<Tt->ne;i++){
    tndofe = Gtt->give_ndofe (i);
    
    reallocm (RSTCKMAT(tndofe,tndofe,lm));
    conductmat (i,lcid,lm);
    
    reallocv (RSTCKIVEC(tndofe,rcn));
    Gtt->give_code_numbers (i,rcn.a);
    
    D0mat->localize (lm,rcn,i);
  }
    
  //fprintf (Outc,"\n\n\n trfel\nMATICE VODIVOSTI \n\n\n");
  //D0mat->printmat (Outc);

  //  contributions from METR upper part
  for (i=0;i<Ct->ne;i++){
    mndofe = Gtm->give_ndofe (i);
    tndofe = Gtt->give_ndofe (i);
    
    reallocm (RSTCKMAT(mndofe,tndofe,lm));
    upper_cond_coupl_mat (i,lcid,lm);
    
    reallocv (RSTCKIVEC(mndofe,rcn));
    reallocv (RSTCKIVEC(tndofe,ccn));

    Gtm->give_code_numbers (i,rcn.a);
    Gtt->give_code_numbers (i,ccn.a);
    
    D0mat->glocalize (lm,rcn,ccn);
  }
    
  //fprintf (Outc,"\n\n\n metru\nMATICE VODIVOSTI-TUHOSTI \n\n\n");
  //D0mat->printmat (Outc);
			  
  //  contributions from METR lower part
  for (i=0;i<Ct->ne;i++){
    mndofe = Gtm->give_ndofe (i);
    tndofe = Gtt->give_ndofe (i);
    
    reallocm (RSTCKMAT(tndofe,mndofe,lm));
    lower_cond_coupl_mat (i,lcid,lm);
    
    reallocv (RSTCKIVEC(tndofe,rcn));
    reallocv (RSTCKIVEC(mndofe,ccn));

    Gtm->give_code_numbers (i,ccn.a);
    Gtt->give_code_numbers (i,rcn.a);
    
    D0mat->glocalize (lm,rcn,ccn);
  }
  
  //fprintf (Outc,"\n\n\n metrl\nMATICE VODIVOSTI-TUHOSTI \n\n\n");
  //D0mat->printmat (Outc);
  
  //D0mat->prepmat (Gtu,Cp->zero,Mesprc);
  D0mat->prepmat (Cp->zero,Mesprc);
}


/**
   function assembles zero order %matrix
   the %matrix multiplies the %vector of nodal values
   
   JK, 11.4.2019
*/
void zero_order_matrix_fc ()
{
  long i,ndofe;
  ivector rcn,ccn;
  matrix km;
  
  if (D0mat==NULL)  D0mat = new gmatrix;
  
  D0mat->setval (Cp->ssle);
  D0mat->initiate (Gtu,Ndofc,Cp->tstord0,Mesprc,Outc);
  
  for (i=0;i<Ct->ne;i++){
    //  the number of all DOFs on element
    ndofe = Gtu->give_ndofe (i);

    //  zero order matrix of a single element
    reallocm (RSTCKMAT(ndofe,ndofe,km));
    zero_order_matrix (i,km);
    
    //  code numbers
    reallocv (RSTCKIVEC(ndofe,rcn));
    Gtu->give_code_numbers (i,rcn.a);
    
    //  localization of the local matrix
    D0mat->localize (km,rcn,i);
  }
  
  D0mat->prepmat (Cp->zero,Mesprc);
}


/**
   function assembles general capacity %matrix
   
   @param lcid - load case id

   JK
*/
void first_order_matrix (long lcid)
{
  long i,mndofe,tndofe;
  ivector rcn,ccn;
  matrix lm;

  if (D1mat==NULL)  D1mat = new gmatrix;

  //  sdruzena topologie
  //D1mat->ts=Cp->tstord1;
  D1mat->setval (Cp->ssle);
  D1mat->initiate (Gtu,Ndofc,Cp->tstord1,Mesprc,Outc);

  
  //D1mat->setval (Cp->ssle.tlinsol,Cp->ssle.prec.pt,Cp->ssle.ni,Cp->ssle.err);


  //  contributions from MEFEL part
  /* for (i=0;i<Mt->ne;i++){
     mndofe = Gtm->give_ndofe (i);
     
     reallocm (RSTCKMAT(mndofe,mndofe,lm));
     stiffmat (lcid, i, lm);
     
     reallocv (RSTCKIVEC(mndofe,rcn));
     //rcn = new long [mndofe];
     Gtm->give_code_numbers (i,rcn.a);
     
     D1mat->localize (lm,rcn,i);
     
     delete [] rcn;
     }
  */

  //  contributions from TRFEL part
  for (i=0;i<Tt->ne;i++){
    tndofe = Gtt->give_ndofe (i);
    
    reallocm (RSTCKMAT(tndofe,tndofe,lm));
    capacmat (i,lcid,lm);
    
    reallocv (RSTCKIVEC(tndofe,rcn));
    //rcn = new long [tndofe];
    Gtt->give_code_numbers (i,rcn.a);
    
    D1mat->localize (lm,rcn,i);

  }



  //  contributions from METR upper part
  for (i=0;i<Ct->ne;i++){
    mndofe = Gtm->give_ndofe (i);
    tndofe = Gtt->give_ndofe (i);
    
    reallocm (mndofe,tndofe,lm);
    upper_cap_coupl_mat (i,lcid,lm);
    
    reallocv (mndofe,rcn);
    //rcn = new long [mndofe];
    reallocv (tndofe,ccn);
    //ccn = new long [tndofe];
    Gtm->give_code_numbers (i,rcn.a);
    Gtt->give_code_numbers (i,ccn.a);
    
    D1mat->glocalize (lm,rcn,ccn);
    
  }

  //  contributions from METR lower part
  for (i=0;i<Ct->ne;i++){
    mndofe = Gtm->give_ndofe (i);
    tndofe = Gtt->give_ndofe (i);
    
    reallocm (tndofe,mndofe,lm);
    lower_cap_coupl_mat (i,lcid,lm);
    
    reallocv (tndofe,rcn);
    //rcn = new long [tndofe];
    reallocv (mndofe,ccn);
    //ccn = new long [mndofe];
    Gtm->give_code_numbers (i,ccn.a);
    Gtt->give_code_numbers (i,rcn.a);
    
    //pred lokalizaci
    //D1mat->printmat (Outc);
    //fflush(Outc);
    D1mat->glocalize (lm,rcn,ccn);
    //po lokalizaci
    //D1mat->printmat (Outc);
    //fflush(Outc);
  }
  


  //D1mat->prepmat (Gtu,Cp->zero,Mesprc);
  D1mat->prepmat (Cp->zero,Mesprc);
}


/**
   function assembles first order %matrix
   the %matrix multiplies the %vector of time derivatives of nodal values
   
   JK, 11.4.2019
*/
void first_order_matrix_fc ()
{
  long i,ndofe;
  ivector rcn,ccn;
  matrix cm;

  if (D1mat==NULL)  D1mat = new gmatrix;
  D1mat->setval (Cp->ssle);
  D1mat->initiate (Gtu,Ndofc,Cp->tstord1,Mesprc,Outc);
  
  for (i=0;i<Ct->ne;i++){
    //  the number of all DOFs on element
    ndofe = Gtu->give_ndofe (i);

    //  first order matrix of a single element
    reallocm (ndofe,ndofe,cm);
    first_order_matrix (i,cm);
    
    //  code numbers
    reallocv (ndofe,rcn);
    Gtu->give_code_numbers (i,rcn.a);

    //  localization of the local matrix
    D1mat->localize (cm,rcn,i);
  }

  D1mat->prepmat (Cp->zero,Mesprc);
}

/**
   function assembles residual %vector
   the %vector contains residuals in mechanical and transport parts
   mechanical contributions are obtained by the function compute_internal_forces
   transport part is obtained by %matrix-%vector multiplication
   
   @param res - residual %vector
   
   JK, 15.4.2019
*/
void residual_vector (double *res)
{
  long lcid=0;
  
  nullv (res,Ndofc);
  
  //  function assembles mechanical part of the residual vector
  internal_forces (lcid,res);
  
  
}



/**
   function computes %vector of generalized internal forces
   general internal force is notation of internal force or internal flux
   
   @param lcid   - load case id
   @param intfor - internal forces and fluxes

   JK, 28.7.2001
*/
void internal_gforces (long lcid,double *intfor)
{
  double *ifort,*iform;

  ifort= new double [Ndofc];
  nullv (ifort,Ndofc);
  iform= new double [Ndofc];
  nullv (iform,Ndofc);
  nullv (intfor,Ndofc);

  //  actual nodal forces MEFEL
  internal_forces (lcid,iform);
  addv (intfor,iform,Ndofc);

  //  actual nodal fluxes TRFEL
  internal_fluxes (ifort,Ndofc);
  addv (intfor,ifort,Ndofc);

  // actual nodal forces  and fluxes METR
  internal_coup_fluxes (lcid,intfor);
  internal_coup_forces (lcid,intfor);

  //not completed yet
  //add_internal_coup_forces_to_mefel();//correct forces for MEFEL printing
  //add_internal_coup_fluxes_to_trfel();//correct fluxes for TRFEL printing

  
  delete [] iform;
  delete [] ifort;
}



/**
   function computes %vector of increments of generalized internal forces
   general internal force is notation of internal force or internal flux
   
   @param lcid   - load case id
   @param intfor - increments of internal forces and fluxes

   12/9/2008, TKr
*/
void incr_internal_gforces (long lcid,double *intfor)
{
  nullv (intfor,Ndofc);
  
  //  actual increments of nodal forces MEFEL
  double *iform = new double [Ndofc];
  nullv (iform,Ndofc);
  incr_internal_forces (lcid,iform);
  addv (intfor,iform,Ndofc);
  delete [] iform;
  
  //  actual increments of nodal fluxes TRFEL
  //double *ifort = new double [Ndofc];
  //nullv (ifort,Ndofc);
  //incr_internal_fluxes (ifort,Ndofc);//not completed yet
  //addv (intfor,ifort,Ndofc);
  //delete [] ifort;

  // actual increments of nodal forces  and fluxes METR
  //incr_internal_coup_fluxes (av,Ndofc);//not completed yet
  //incr_internal_coup_forces (lcid,intfor);//not completed yet  
}



/**
   function extracts TRFEL values on one element

   @param lcid - number of load case
   @param r - vector of nodal values
   @param cn - array containing code numbers
   @param ndofe - number of DOFs on actual element

   JK, 9.7.2001
*/
void elemvaluesc (long lcid,double *r,long *cn,long ndofe)
{
  long i,ii;

  for (i=0;i<ndofe;i++){
    ii=cn[i];
    if (ii<0)   r[i]=Tb->lc[lcid].pv[0-ii-1].getval(Tp->time);
    if (ii==0)  r[i]=0.0;
    if (ii>0)   r[i]=Lsrst->lhsi[ii-1]+Lsrst->lhs[ii-1];
  }
}



/**
   function computes %vector of generalized right hand side
   
   @param lcid - load case id
   @param rhs - array of right hand side
   @param flv - array of load %vector caused by forces only (default value is NULL)
   
   @return The function returns assembled vectors of righ hand side and load vector caused by force.

   JK, 28.7.2001,
   Modified by Tomas Krejci 18/05/2012
*/
void metr_right_hand_side (long lcid,double *rhs,double *flv)
{
  double *av;

  av = new double [Ndofc];

  nullv (rhs,Ndofc);
  nullv (flv,Ndofc);
  nullv (av,Ndofc);

  mefel_right_hand_side (lcid,rhs,flv);//MEFEL part
  trfel_right_hand_side (lcid,av,Ndofc);//TRFEL part
  //mechanical influence of transport problem
  trfel_right_hand_side2 (0,av,Ndofc);

  addv (rhs,av,Ndofc);
  
  nullv (av,Ndofc);
  right_hand_side (lcid,av,Ndofc);//METR part
  addv (rhs,av,Ndofc);

  delete [] av;
}

/**
   function computes %vector of generalized right hand side
   it is used in problems with fully coupled material models
   
   @param rhs - array of right hand side
   
   JK, 10.4.2019
*/
void metr_right_hand_side_fc (double */*rhs*/)
{

}

/**
   function updates material data on element int. points
   
   14/05/2010, TKr
*/
void updateval()
{
  //updating of int. points material models
  Tm->updateipval ();
  Tm->update_aipval ();
  Mm->updateipval();
  Mm->update_aipval();
  Cmu->updateipval();
  //Cml->updateipval(); //not created now
}



/**
   function approximates data on element int. points
   
   JK, 28.7.2001
*/
void approximationc ()
{
  //  approximation of nodal values into integration points of TRFEL
  approximation ();  
  
  //  approximation of nodal values into integration points of METR
  approximationcoup ();
}



/**
   function approximates data on element int. points on coupling mesh
   
   JK, 28.7.2001
*/
void approximationcoup ()
{
  long i;
  
  switch (Cp->tmatt){
  case mech_onemedium:{
    for (i=0;i<Ct->ne;i++){
      intpointvaluesc (i);
      intpointstrainsc (i);
      intpointgradientsc (i);
    }
    break;
  }
  case mech_twomedia:{
    for (i=0;i<Ct->ne;i++){
      intpointvaluesc (i);
      intpointstrainsc (i);
      intpointgradientsc (i);
    }
    break;
  }
  case mech_threemedia:{
    for (i=0;i<Ct->ne;i++){
      intpointvaluesc (i);
      intpointstrainsc (i);
      intpointgradientsc (i);
    }
    break;
  }
  default:{
    print_err("unknown transported matter is required",__FILE__,__LINE__,__func__);
    abort();
  }
  }
  
}

/**
   function approximates nodal values to integration points
   
   JK, 11.4.2019
*/
void approximationc_fc ()
{
 
}



/**
   function assembles coupling right hand side of the problem
   
   @param lcid - load case id
   @param rhs - array containing right hand side
   
   09.03.2004, JK+TKr
*/
void right_hand_side (long lcid,double *rhs,long n)
{
  long i,mndofe,*cn;
  vector lv;
  double *av;
  av = new double [n];
  
  assemble_init_coup (rhs);

  switch (Cp->tmatt){
  case mech_onemedium:{
    nullv (av,n);
    assemble_coup (lcid,av,n);
    addv (rhs,av,n);
    break;
  }
  case mech_twomedia:{
    nullv (av,n);
    assemble_coup (lcid,av,n);
    addv (rhs,av,n);
    break;
  }
  case mech_threemedia:{
    nullv (av,n);
    assemble_coup (lcid,av,n);
    addv (rhs,av,n);
    break;
  }
  default:{
    print_err("unknown transported matter is required", __FILE__, __LINE__, __func__);
  }
  }
  
  //  contributions from volume integrals
  nullv (av,n);
  for (i=0;i<Ct->ne;i++){
    if (Gtm->leso[i]==1){
      //  only elements switched on are processed
      
      //  number of DOFs on element
      mndofe=Gtm->give_ndofe (i);
      reallocv (mndofe,lv);
      fillv (0.0,lv);
      
      //  function is defined in elemswitchc.cpp
      volume_rhs_vectorc (lv,lcid,i);
      
      
      //  code numbers
      cn = new long [mndofe];
      Gtm->give_code_numbers (i,cn);
      
      //  localization of element values to the global vector
      locglob (av,lv.a,cn,mndofe);
      
      delete [] cn;
    }  
  }  
  
  addv (rhs,av,n);
  delete [] av;
}



/**
   function assembles part of the %vector of coupling right hand side
   
   @param rhs - right hand side
   
   21.3.2004, JK+TKr
   09/11/2012 corrected by TKr,
   20/04/2023 another correction by TKr
*/
void assemble_init_coup (double *rhs)
{
  long i,j,mndofe,tndofe,lcid=0;
  double dt;
  ivector rcn,ccn;
  vector r,f,pr;
  matrix km,cm;
  
  for (i=0;i<Ct->ne;i++){
    mndofe = Gtm->give_ndofe (i);
    tndofe = Gtt->give_ndofe (i);
    
    // upper block for METR-MEFEL part
    reallocv (RSTCKVEC(tndofe,r));
    reallocv (RSTCKVEC(mndofe,f));
    reallocv (RSTCKVEC(tndofe,pr));
    reallocm (RSTCKMAT(mndofe,tndofe,km));
    reallocv (RSTCKIVEC(mndofe, rcn));
    reallocv (RSTCKIVEC(tndofe, ccn));

    upper_cond_coupl_mat (i,lcid,km);
    Gtm->give_code_numbers (i,rcn.a);
    Gtt->give_code_numbers (i,ccn.a);
    
    // prescribed values from TRFEL
    prescvalues (r.a,ccn.a,tndofe);

    mxv (km,r,f);
    cmulv (-1.0,f);
    
    locglob (rhs,f.a,rcn.a,mndofe);

    // contribution from time derivative of prescribed values on boundaries
    //  prescribed values from the previous time step
    prevprescvalues (pr.a,rcn.a,tndofe);
    //  the length of time step
    dt=Tp->timecont.forwarddt;
    //  computation of the time derivative of
    for (j=0;j<tndofe;j++){
      pr.a[j]=(r.a[j]-pr.a[j])/dt;
    }
    
    reallocm (RSTCKMAT(mndofe,tndofe,cm));
    upper_cap_coupl_mat (i,lcid,cm);
    check_math_errel(i);
    
    //  nodal values
    mxv (cm,pr,f);
    cmulv (-1.0,f);
    
    //  localization of nodal values
    locglob (rhs,f.a,rcn.a,mndofe);
    
 
    //lower block for METR-TRFEL part
    //conductivity - stiffness matrix:
    reallocv (RSTCKVEC(mndofe,r));
    reallocv (RSTCKVEC(tndofe,f));
    reallocv (RSTCKVEC(mndofe,pr));
    reallocm (RSTCKMAT(tndofe,mndofe,km));
    reallocv (RSTCKIVEC(tndofe, rcn));
    reallocv (RSTCKIVEC(mndofe, ccn));

    lower_cond_coupl_mat (i,lcid,km);
    
    Gtt->give_code_numbers (i,rcn.a);
    Gtm->give_code_numbers (i,ccn.a);
    
    // prescribed values (displacements) from MEFEL
    // prescdispl (r.a,ccn.a,mndofe); // not implemented in MEFEL
    elprdispl (lcid,i,r.a);

    mxv (km,r,f);
    cmulv (-1.0,f);
    
    locglob (rhs,f.a,rcn.a,tndofe);

    // contribution from time derivative of prescribed values on boundaries - displacements
    // prevprescdislp (pr.a,rcn.a,mndofe); // not implemented in MEFEL
    elprevprdispl (lcid,i,pr.a);
    //  the length of time step
    dt=Tp->timecont.forwarddt;
    //  computation of the time derivative of
    for (j=0;j<tndofe;j++){
      pr.a[j]=(r.a[j]-pr.a[j])/dt;
    }
    
    reallocm (RSTCKMAT(tndofe,mndofe,cm));
    lower_cap_coupl_mat (i,lcid,cm);
    check_math_errel(i);
    
    //  nodal values
    mxv (cm,pr,f);
    cmulv (-1.0,f);
    
    //  localization of nodal values
    locglob (rhs,f.a,ccn.a,tndofe);
  }
}



/**
  The function passes coupling data from TRFEL to MEFEL and vice 
  versa simulatneously.    

  @param lcid[in] - load case id

  14/05/2010, TKr
  Rewritten by Tomas Koudelka, 12.6.2018
*/
void pass_coup_data(long lcid)
{
  //  passing of required transport data to MEFEL
  trfel_mefel();
  //  passing of required mechanical data to TRFEL
  mefel_trfel(lcid);
}



/**
   function exports data from METR into MEFEL
   
   05/06/2012, TKr, not completed, tady pokracovat !!!
*/
void metr_mefel ()
{
}



/**
  The function exports data from MEFEL into TRFEL
  All required non-transport quantities (known in MEFEL) are recalculated to 
  TRFEL integration points and stored in Tm->nontransq.
   
  @param lcid[in] - load case id

  05/06/2012, TKr
  Modified by TKo 10.10.2013
  Rewritten by TKo 12.6.2018
*/
void mefel_trfel (long lcid)
{
  switch (Cp->dpt)
  {
    case pass_by_closest_ip:
      mefel_trfel_by_nodes();
      break;
    case pass_by_nodes_comp:
      mefel_trfel_by_nodes_comp();
      break;
    case pass_by_aux_ip:
      mefel_trfel_by_aip(lcid, Tm->tnip, TMipmap);
      break;
    case pass_by_copy_ip:
      mefel_trfel_copyip();
      break;     
    default:
      print_err("unknown type of coupled data passing %d is required", __FILE__, __LINE__, __func__);
      abort();
  }
}



/**
  The function passes coupling data between MEFEL and TRFEL by nodal values. Values are copied to nodes 
  from the closest int. points at the particular MEFEL elements and then passed to TRFEL 
  elements which approximate them to the TRFEL int. points.

  @return The function does not return anything.

  Created by Tomas Koudelka, 14.11.2013
*/
void mefel_trfel_by_nodes(void)
{
  long i, j, k, nnet, nnem, nnv, nip, ipid, nentq;
  ivector nodes;
  vector nodval, ipval;
  long antq[Tm->tnkntq];
  nontransquant *entqo; 

  for (i=0;i<Tt->ne;i++)
  {
    if (Gtt->leso[i]==1)
    {      
      memset(antq, 0, sizeof(*antq)*Tm->tnkntq);
      ipid=Tt->elements[i].ipp[0][0];
      for(j=0; j<Tt->give_tnip(i); j++)
      {
        Tm->give_reqntq(ipid, antq); // search for required non-transport quantities on the i-th integration point
        ipid++;
      }

      for (j=0, nentq=0; j<Tm->tnkntq; j++)
      {
        if (antq[j] == 1)   nentq++;  // compute number of required non-transport quantities on i-th element
      } 

      entqo = new nontransquant[nentq];      
      for (j=0, k=0; j<Tm->tnkntq; j++)
      {
        if (antq[j] == 1)
        {   
          entqo[k]= nontransquant(j+1);  // assemble arrray with ordering of required non-transport quantities on element
          k++;
        }
      } 
      
      nnet = Tt->give_nne (i); // number of nodes on transport element      
      nnem = Mt->give_nne (i); // number of nodes on mechanical element
      nnv = (nnet < nnem) ? nnet : nnem; // minimum number of nodes
      reallocv (RSTCKVEC(nnv,nodval));
      
      //  interpolation of individual non-transport (mechanical) quantities      
      nip = Tt->give_tnip (i); // number of int. points on trnasport element
      reallocv (nip,ipval);
      for (j=0; j<nentq; j++)
      { 
        if (Cp->bb==quad_quad || Cp->bb==lin_lin){
          //  nodal values of all required nontransport(=mechanical) quantities on mechanical element
          elem_mechq_nodval(i, nodval, entqo[j]);
        }
        if (Cp->bb==quad_lin){
          //  nodal values of all required nontransport(=mechanical) quantities on mechanical element
          elem_mechq_nodval2(i, nodval, entqo[j]);
        }

        if (nnet > nnem) // TRFEL element approximation is quadratic but MEFEL is linear
        {
          // it is assumed that order of corner nodes on a quadratic element 
          // is identical with order on linear one.

          // elem_intpointvalt2 (i, nodval, ipval);//tady toto jeste rozmyslet, jestli to ma vyznam by TKr??!!
          print_err("TRFEL element approximation order is higher than on MEFEL element", __FILE__, __LINE__, __func__);
          abort();
        }
        else
        {
          // MEFEL and TRFEL element approximation are identical =>
          // it is assumed that order of nodes on TRFEL and MEFEL elements are identical

          // MEFEL element approximation is quadratic but TRFEL is linear =>
          // nothing is necessary because it is assumed that order of corner nodes on a quadratic
          // element is identical with order on linear one.

          //  interpolation of nodal values to integration points
          elem_intpointvalt (i, nodval, ipval);
        }
      
        //  number of the first integration point
        ipid=Tt->elements[i].ipp[0][0];
      
        for (k=0;k<nip;k++)
        {
          Tm->storenontransq(entqo[j], ipid, ipval[k]);
          ipid++;
        }
      }
    }
  }
}



/**
  The function computes and passes nodal values used in the case of passing of coupling data
  between MEFEL and TRFEL. Nodal values are computed at particular MEFEL element nodes and then 
  passed to TRFEL elements which approximate them to the TRFEL int. points.

  @return The function does not return anything.

  Created by Tomas Koudelka, 14.11.2013
*/
void mefel_trfel_by_nodes_comp(void)
{
  long i, j, k, nnet, nnem, nnv, nip, ipid, nentq;
  ivector nodes;
  vector nodval, ipval, auxnv;
  long antq[Tm->tnkntq];
  nontransquant *entqo; 

  for (i=0;i<Tt->ne;i++)
  {
    if (Gtt->leso[i]==1)
    {      
      memset(antq, 0, sizeof(*antq)*Tm->tnkntq);
      ipid=Tt->elements[i].ipp[0][0];
      for(j=0; j<Tt->give_tnip(i); j++)
      {
        Tm->give_reqntq(ipid, antq); // search for required non-transport quantities on the i-th integration point
        ipid++;
      }

      for (j=0, nentq=0; j<Tm->tnkntq; j++)
      {
        if (antq[j] == 1)   nentq++;  // compute number of required non-transport quantities on i-th element
      } 

      entqo = new nontransquant[nentq];      
      for (j=0, k=0; j<Tm->tnkntq; j++)
      {
        if (antq[j] == 1)
        {   
          entqo[k]= nontransquant(j+1);  // assemble arrray with ordering of required non-transport quantities on element
          k++;
        }
      } 
      
      nnet = Tt->give_nne (i); // number of nodes on transport element      
      nnem = Mt->give_nne (i); // number of nodes on mechanical element
      nnv = (nnet < nnem) ? nnet : nnem; // minimum number of nodes
      reallocv (RSTCKVEC(nnv*Tm->nntq,nodval));
      
      if (Cp->bb==quad_quad || Cp->bb==lin_lin){
	//  nodal values of all required nontransport(=mechanical) quantities on mechanical element
	elem_mechq_nodval_comp(i, nodval, nnv, nentq, entqo);
      }
      if (Cp->bb==quad_lin){
	//  nodal values of all required nontransport(=mechanical) quantities on mechanical element
	elem_mechq_nodval_comp2(i, nodval, nnv, nentq, entqo);
      }

      //  interpolation of individual non-transport (mechanical) quantities      
      auxnv.n = nnet;
      nip = Tt->give_tnip (i); // number of int. points on trnasport element
      reallocv (RSTCKVEC(nip,ipval));
      for (j=0; j<nentq; j++)
      { 
        auxnv.a = nodval.a + j*nnv;

        if (nnet > nnem) // TRFEL element approximation is quadratic but MEFEL is linear
        {
          // it is assumed that order of corner nodes on a quadratic element 
          // is identical with order on linear one.

          // elem_intpointvalt2 (i, auxnv, ipval);//tady toto jeste rozmyslet, jestli to ma vyznam by TKr??!!
          print_err("TRFEL element approximation order is higher than on MEFEL element", __FILE__, __LINE__, __func__);
          abort();
        }
        else
        {
          // MEFEL and TRFEL element approximation are identical =>
          // it is assumed that order of nodes on TRFEL and MEFEL elements are identical

          // MEFEL element approximation is quadratic but TRFEL is linear =>
          // nothing is necessary because it is assumed that order of corner nodes on a quadratic
          // element is identical with order on linear one.

          //  interpolation of nodal values to integration points
          elem_intpointvalt (i, auxnv, ipval);
        }
      
        //  number of the first integration point
        ipid=Tt->elements[i].ipp[0][0];
      
        for (k=0;k<nip;k++)
        {
          Tm->storenontransq(entqo[j], ipid, ipval[k]);
          ipid++;
        }
      }
      delete [] entqo;
    }
  }
  auxnv.a = NULL; // avoid memory deallocation error
}



/**
  The function computes/passes coupling data between MEFEL and TRFEL. Data are taken from the auxiliary integration points 
  in MEFEL and stored in TRFEL nontransq arrays. 

  @param lcid[in] - load case id
  @param n[in]    - the number of auxiliary points in the mapping array ipm, in this case it should the 
                    number of regular int. points of the TRFEL module.
  @param ipm[in]  - integration point mapping array, 
                    ipm[i].ipp < 0 => auxiliary integration point must be used
                    ipm[i].ipp >= 0 => direct mapping to the regulal integration point, 
                                       no computation of values is performed, the values are assumed 
                                       to be computed at the main solution procedure of the problem

  @return The function does not return anything but it computes strains, stresses and state variables for the
          given load case at auxiliary integration points and passes them in TRFEL nontransq arrays.

  Created by Tomas Koudelka, 13.06.2018
*/
void mefel_trfel_by_aip(long lcid, long n, ipmap *ipm)
{
  long i, j, ipp, app;
  const long tnkntq = Tm->tnkntq;
  long antq[tnkntq];
  long nntq = Tm->nntq;
  nontransquant *ntqo = Tm->ntqo;
  intpoints *tmp_ip;
  long      *tmp_elip;
  double   **tmp_nonmechq;
  double   **tmp_ic;
  double   **tmp_eigstrains;
  double   **tmp_eigstresses;
  double   **tmp_tempstrains;

  // pass the coupling data to the MEFEL auxiliary int. points
  actualize_aip_nonmechq(n, ipm);
  // actualize values of eigenstrains/eigenstresses
  if (Mp->eigstrains && Mp->eigstrains != 3)
    Mb->aip_eigstrain_computation(n, ipm, Mp->time);
  // compute actual strains at MEFEL auxiliary int. points
  compute_aipstrains(lcid, n, ipm);
  // compute actual stresses and state variables at MEFEL auxiliary int. points
  compute_aipstresses(lcid);
  // update eqother array of auxiliary int. points at MEFEL
  Mm->update_aipval();

  //
  // processing of auxiliary int. points without direct mapping to regular int. points
  //
  
  // swap regular integration point auxiliary integration point arrays
  // from now, all functions will work with AUXILIARY integration points 
  tmp_ip          = Mm->ip;
  tmp_elip        = Mm->elip;
  tmp_nonmechq    = Mm->nonmechq;
  tmp_ic          = Mm->ic;
  tmp_eigstrains  = Mm->eigstrains;
  tmp_eigstresses = Mm->eigstresses;
  tmp_tempstrains = Mm->tempstrains;

  Mm->tnip         = Mm->tnaip;
  Mm->ip           = Mm->aip;
  Mm->elip         = Mm->elaip;
  Mm->nonmechq     = Mm->aip_nonmechq;
  Mm->ic           = Mm->aip_ic;
  Mm->eigstrains   = Mm->aip_eigstrains;
  Mm->eigstresses  = Mm->aip_eigstresses;
  Mm->tempstrains  = Mm->aip_tempstrains;

  for(i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app >= 0) // passing values from auxiliary int. point 
    {
      // obtain array of indicators of the required nontransport quantities
      // at the given MEFEL auxiliary int. point = TRFEL regular int. point, 
      // i.e. index i in the array of MEFEL aux. int points = ipp in TRFEL regular int. points
      memset(antq, 0, sizeof(*antq)*tnkntq);
      Tm->give_reqntq(i, antq);
      for (j=0; j<nntq; j++)
      {     
        if (antq[ntqo[j]-1])
        {
          // the quantity is required at the given TRFEL int. point ->
          // -> take value from the auxiliary int. point stored in MEFEL mechmat class
          // at this point arrays of regular int. points are swapped with the auxiliary one
          // and therefore all functions of mechmat work with aux. int. points
          Tm->nontransq[j][i] = Mm->givemechq(ntqo[j], app); 
        }
      }
    }
  }
  // swap regular integration point and auxiliary integration point arrays
  // from now, all functions will work with REGULAR integration points 
  Mm->tnip         = Mm->tnrip;
  Mm->ip           = tmp_ip;
  Mm->elip         = tmp_elip;
  Mm->nonmechq     = tmp_nonmechq;
  Mm->ic           = tmp_ic;
  Mm->eigstrains   = tmp_eigstrains;
  Mm->eigstresses  = tmp_eigstresses;
  Mm->tempstrains  = tmp_tempstrains;

  //
  // processing of auxiliary int. points mapped to the regular int. points directly
  //

  for(i=0; i<n; i++)
  {
    ipp = ipm[i].ipp;
    // check direct mapping between auxiliary and regular int. point
    if (ipp >= 0)
    {
      // obtain array of indicators of the required nontransport quantities
      // at the given MEFEL auxiliary int. point = TRFEL regular int. point, 
      // i.e. index i in the array of MEFEL aux. int points = ipp in TRFEL regular int. points
      memset(antq, 0, sizeof(*antq)*tnkntq);
      Tm->give_reqntq(i, antq);
      for (j=0; j<nntq; j++)
      {
        if (antq[ntqo[j]-1])
        {
          // the quantity is required at the given TRFEL int. point ->
          // -> take value from the regular int. point stored in MEFEL mechmat class
          Tm->nontransq[j][i] = Mm->givemechq(ntqo[j], ipp);
        }
      }
    }
  }
}



/**
  The function exports data of TRFEL into MEFEL.
  All required non-mechanical quantities (known in TRFEL) are recalculated to 
  MEFEL integration points and stored in Mm->nonmechq.
   
   
  05/06/2012, TKr
  Modified by TKo 10.10.2013
  Rewritten by TKo 12.6.2018
*/
void trfel_mefel ()
{
  switch (Cp->dpt)
  {
    case pass_by_closest_ip:
      trfel_mefel_by_nodes();
      break;
    case pass_by_nodes_comp:
      trfel_mefel_by_nodes_comp();
      break;
    case pass_by_aux_ip:
      trfel_mefel_by_aip(Mm->tnip, MTipmap);
      break;
    case pass_by_copy_ip:
      trfel_mefel_copyip();
      break;     
    default:
      print_err("unknown type of coupled data passing %d is required", __FILE__, __LINE__, __func__);
      abort();
  }
}



/**
  The function passes coupling data between TRFEL and MEFEL by nodal values. Values are copied to nodes 
  from the closest int. points at the particular TRFEL elements and then passed to MEFEL 
  elements which approximate them to the MEFEL int. points.

  @return The function does not return anything.

  Created by Tomas Koudelka, 5.2018
*/
void trfel_mefel_by_nodes(void)
{
  long i, j, k, nnet, nnem, nip, ipid, nenmq;
  ivector nodes;
  vector nodval, ipval, auxnv;
  long anmq[Mm->tnknmq];
  nonmechquant *enmqo = new nonmechquant[Mm->tnknmq];
  
  for (i=0;i<Mt->ne;i++)
  {
    if (Gtm->leso[i]==1)
    {
      nnem = Mt->give_nne (i); // number of nodes on mechanical element
      nip = Mt->give_tnip (i); // number of int. points on mechanical element
      reallocv (RSTCKIVEC(nnem,nodes));
      Mt->give_elemnodes (i, nodes); //  nodes on element
      
      memset(anmq, 0, sizeof(*anmq)*Mm->tnknmq);

      ipid=Mt->elements[i].ipp[0][0];
      for(j=0; j<Mt->give_tnip(i); j++)
      {
        for (k=0; k<Mm->ip[ipid].nm; k++)
          Mm->give_reqnmq(ipid, k, anmq); // search for required non-mechanical quantities on the i-th integration point
        ipid++;
      }

      for (j=0, nenmq=0; j<Mm->tnknmq; j++)
      {
        if (anmq[j] == initial_temperature) // initial temperature can be passed only once in init_trfel_mefel()
          continue;

        if (anmq[j] == 1)   nenmq++;  // compute number of required non-mechanical quantities on i-th element
      } 

      if (nenmq == 0)
        continue;

      reallocv (RSTCKVEC(nnem,nodval));
      reallocv (RSTCKVEC(nip,ipval));
      memset(enmqo, 0, sizeof(*enmqo)*Mm->tnknmq);      

      for (j=0, k=0; j<Mm->tnknmq; j++)
      {
        if (anmq[j] == initial_temperature) // initial temperature can be passed only once in init_trfel_mefel()
          continue;

        if (anmq[j] == 1)
        {
          enmqo[k]= nonmechquant(j+1);  // assemble array with ordering of required non-mechanical quantities on element
          k++;
        }
      } 


      nnet = Tt->give_nne (i); // number of nodes on transport element
      
      //  interpolation of individual non-mechanical(=transport) quantities
      auxnv.n = nnem;
      for (j=0; j<nenmq; j++)
      { 
        // nodal values of all required non-mechanical(=transport) quantities on transport element
        // because the order of elements in transport part must be the same or lower then in mechanical part
        // only one function elem_transq_nodval may be used
        elem_transq_nodval(i, nodval, enmqo[j]);

        if (nnem > nnet) // MEFEL element approximation is quadratic but TRFEL is linear
        {
          // it is assumed that order of corner nodes on a quadratic element 
          // is identical with order on linear one.

          elem_intpointval2 (i, nodval, ipval);
        }
        else
        {
          // MEFEL and TRFEL element approximation are identical =>
          // it is assumed that order of nodes on TRFEL and MEFEL elements are identical

          // TRFEL element approximation is quadratic but MEFEL is linear
          // nothing is necessary because it is assumed that order of corner nodes on a quadratic
          // element is identical with order on linear one.
          
          //  interpolation of nodal values to integration points
          elem_intpointval (i, nodval, ipval);
        }

      
        //  number of the first integration point
        ipid=Mt->elements[i].ipp[0][0];
      
        for (k=0;k<nip;k++)
        {
          Mm->storenonmechq(enmqo[j], ipid, ipval[k]);
          ipid++;
        }
      }
    }
  }
  delete [] enmqo;
}



/**
  The function computes and passes nodal values used in the case of passing data
  between TRFEL and MEFEL. Nodal values are computed at particular TRFEL element nodes and then 
  passed to MEFEL elements which approximate them to the MEFEL int. points.

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.12.2013
*/
void trfel_mefel_by_nodes_comp(void)
{
  long i, j, k, nnet, nnem, nnv, nip, ipid, nenmq;
  ivector nodes;
  vector nodval, ipval, auxnv;
  long anmq[Mm->tnknmq];
  nonmechquant *enmqo; 
  
  for (i=0;i<Mt->ne;i++)
  {
    if (Gtm->leso[i]==1)
    {
      memset(anmq, 0, sizeof(*anmq)*Mm->tnknmq);
      ipid=Mt->elements[i].ipp[0][0];
      for(j=0; j<Mt->give_tnip(i); j++)
      {
        for (k=0; k<Mm->ip[ipid].nm; k++)
          Mm->give_reqnmq(ipid, k, anmq); // search for required non-mechanical quantities on the i-th integration point
        ipid++;
      }

      for (j=0, nenmq=0; j<Mm->tnknmq; j++)
      {
        if (anmq[j] == initial_temperature) // initial temperature can be passed only once in init_trfel_mefel()
          continue;

        if (anmq[j] == 1)   nenmq++;  // compute number of required non-mechanical quantities on i-th element
      } 

      if (nenmq == 0)
        continue;

      enmqo = new nonmechquant[nenmq];      
      for (j=0, k=0; j<Mm->tnknmq; j++)
      {
        if (anmq[j] == initial_temperature) // initial temperature can be passed only once in init_trfel_mefel()
          continue;

        if (anmq[j] == 1)
        {   
          enmqo[k]= nonmechquant(j+1);  // assemble array with ordering of required non-mechanical quantities on element
          k++;
        }
      } 

      nnet = Tt->give_nne (i); // number of nodes on transport element      
      nnem = Mt->give_nne (i); // number of nodes on mechanical element
      nnv = (nnet < nnem) ? nnet : nnem; // minimum number of nodes
      reallocv (RSTCKVEC(nnv*nenmq,nodval));
      //      reallocv (RSTCKVEC(nnv*Mm->nnmq,nodval));
      
      // nodal values of all required non-mechanical(=transport) quantities on transport element
      // because the order of elements in transport part must be the same or lower then in mechanical part
      // only one function elem_transq_nodval_comp may be used
      elem_transq_nodval_comp(i, nodval, nnv, nenmq, enmqo);

      //  interpolation of individual non-mechanical(=transport) quantities
      auxnv.n = nnem;
      nip = Mt->give_tnip (i); // number of int. points on mechanical element
      reallocv (RSTCKVEC(nip,ipval));
      for (j=0; j<nenmq; j++)
      { 
        auxnv.a = nodval.a + j*nnv;

        if (nnem > nnet) // MEFEL element approximation is quadratic but TRFEL is linear
        {
          // it is assumed that order of corner nodes on a quadratic element 
          // is identical with order on linear one.

          elem_intpointval2 (i, auxnv, ipval);
        }
        else
        {
          // MEFEL and TRFEL element approximation are identical =>
          // it is assumed that order of nodes on TRFEL and MEFEL elements are identical

          // TRFEL element approximation is quadratic but MEFEL is linear
          // nothing is necessary because it is assumed that order of corner nodes on a quadratic
          // element is identical with order on linear one.
          
          //  interpolation of nodal values to integration points
          elem_intpointval (i, auxnv, ipval);
        }

      
        //  number of the first integration point
        ipid=Mt->elements[i].ipp[0][0];
      
        for (k=0;k<nip;k++)
        {
          Mm->storenonmechq(enmqo[j], ipid, ipval[k]);
          ipid++;
        }
      }
      delete [] enmqo;
    }
  }
  auxnv.a = NULL;
}



/**
  The function computes/passes coupling data from TRFEL to MEFEL. Data are taken from the auxiliary integration points 
  in TRFEL and stored in MEFEL nonmechq arrays. 

  @param n[in]   - the number of auxiliary points in the mapping array ipm, in this case it should the 
                   number of regular int. points of the MEFEL module.
  @param ipm[in] - integration point mapping array, 
                   ipm[i].ipp < 0 => auxiliary integration point must be used
                   ipm[i].ipp >= 0 => direct mapping to the regulal integration point, 
                                      no computation of values is performed, the values are assumed 
                                      to be computed at the main solution procedure of the problem

  @return The function does not return anything but it computes nodal values and state variables 
          at auxiliary integration points of TRFEL and passes them in MEFEL nonmechq arrays.

  Created by Tomas Koudelka, 13.06.2018
*/
void trfel_mefel_by_aip(long n, ipmap *ipm)
{
  long i, j, ipp, app;
  const long tnknmq = Mm->tnknmq;
  long anmq[tnknmq];
  long nnmq = Mm->nnmq;
  nonmechquant *nmqo = Mm->nmqo;
  intpointst *tmp_ip;
  long *tmp_elip;
  double *tmp_iv;

  // pass coupling data to the TRFEL auxiliary int. points
  actualize_aip_nontransq(n, ipm);
  // compute actual unknown values and their gradients at TRFEL auxiliary int. points
  aip_approximation(n, ipm);
  // compute other required values at TRFEL auxiliary int. points
  Tm->update_aipval();

  //
  // processing of auxiliary int. points given by natural coordinates
  //
  
  // swap regular integration point and auxiliary integration point arrays
  // from now, all functions will work with AUXILIARY integration points 
  tmp_ip = Tm->ip;
  tmp_elip = Tm->elip;
  tmp_iv = Tm->initval;
  Tm->ip = Tm->aip;
  Tm->elip = Tm->elaip;
  Tm->initval = Tm->aip_initval;
  for(i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app >= 0) // passing values from auxiliary int. point 
    {
      // obtain array of indicators of the required nonmechanical quantities
      // the given TRFEL auxiliary int. point = MEFEL regular int. point, 
      // i.e. index i in the array of TRFEL aux. int points = ipp in MEFEL regular int. points
      memset(anmq, 0, sizeof(*anmq)*tnknmq);
      for (j=0; j<Mm->ip[i].nm; j++)
        Mm->give_reqnmq(i, j, anmq);
      for (j=0; j<nnmq; j++)
      {     
        if (nmqo[j] == initial_temperature) // initial temperature can be passed only once in init_trfel_mefel()
          continue;

        if (anmq[nmqo[j]-1])
        {
          // the quantity is required at the given MEFEL int. point ->
          // -> take value from the auxiliary int. point stored in TRFEL transmat class
          // at this point arrays of regular int. points are swapped with the auxiliary one
          // and therefore all functions of transmat work with aux. int. points
          Mm->nonmechq[j][i] = Tm->givetransq(nmqo[j], app); 
        }
      }
    }
  }
  // swap regular integration point and auxiliary integration point arrays
  // from now, all functions will work with REGULAR integration points 
  Tm->ip = tmp_ip;
  Tm->elip = tmp_elip;
  Tm->initval = tmp_iv;

  //
  // processing of auxiliary int. points mapped to the regular int. points directly
  //
  for(i=0; i<n; i++)
  {
    ipp = ipm[i].ipp;
    // check direct mapping between auxiliary and regular int. point
    if (ipp >= 0)
    {
      // obtain array of indicators of the required nonmechanical quantities
      // the given TRFEL auxiliary int. point = MEFEL regular int. point, 
      // i.e. index i in the array of TRFEL aux. int points = ipp in MEFEL regular int. points
      memset(anmq, 0, sizeof(*anmq)*tnknmq);
      for (j=0; j<Mm->ip[i].nm; j++)
        Mm->give_reqnmq(i, j, anmq);
      for (j=0; j<nnmq; j++)
      {
        if (nmqo[j] == initial_temperature)// initial temperature can be passed only once in init_trfel_mefel()
          continue;

        if (anmq[nmqo[j]-1])
        {
          // the quantity is required at the given MEFEL int. point ->
          // -> take value from the auxiliary int. point stored in TRFEL transmat class
          // at this point arrays of regular int. points are swapped with the auxiliary one
          // and therefore all functions of transmat work with aux. int. points
          Mm->nonmechq[j][i] = Tm->givetransq(nmqo[j], ipp); 
        }
      }
    }
  }
}



/**
  The function passes initial coupling data of TRFEL into MEFEL.

  @return The function does not return anything but it changes data at arrays nonmechq and nontransq
          of mechmat and transmat respectively.

  Rewritten by Tomas Koudelka 13.6.2018
*/
void init_trfel_mefel()
{
  switch (Cp->dpt)
  {
    case pass_by_closest_ip:
      init_trfel_mefel_by_nodes();
      break;
    case pass_by_nodes_comp:
      init_trfel_mefel_by_nodes_comp();
      break;
    case pass_by_aux_ip:
      init_trfel_mefel_by_aip(Mm->tnip, MTipmap);
      break;
    case pass_by_copy_ip:
      init_trfel_mefel_copyip();
      break;
    default:
      print_err("unknown type of coupled data passing %d is required", __FILE__, __LINE__, __func__);
      abort();
  }
}



/**
  The function computes and passes initial nodal values used in the case of passing data
  between TRFEL and MEFEL. Initial nodal values are copied from the closest int. points 
  at particular TRFEL element nodes and then passed to MEFEL elements which approximate
  them to the MEFEL int. points.

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.12.2013
*/
void init_trfel_mefel_by_nodes(void)
{
  long i, j, k, nnet, nnem, nip, ipid, nenmq;
  ivector nodes;
  vector nodval, ipval, auxnv;
  long anmq[Mm->tnknmq];
  nonmechquant *enmqo = new nonmechquant[Mm->tnknmq];
  
  for (i=0;i<Mt->ne;i++)
    {
      
      if(i >= 7446)
      printf("Zde\n");

    if (Gtm->leso[i]==1)
    {
      nnem = Mt->give_nne (i); // number of nodes on mechanical element
      nip = Mt->give_tnip (i); // number of int. points on mechanical element
      reallocv (RSTCKIVEC(nnem,nodes));
      Mt->give_elemnodes (i, nodes); //  nodes on element
      
      memset(anmq, 0, sizeof(*anmq)*Mm->tnknmq);

      ipid=Mt->elements[i].ipp[0][0];
      for(j=0; j<Mt->give_tnip(i); j++)
      {
        for (k=0; k<Mm->ip[ipid].nm; k++)
          Mm->give_reqnmq(ipid, k, anmq); // search for required non-mechanical quantities on the i-th integration point
        ipid++;
      }

      for (j=0, nenmq=0; j<Mm->tnknmq; j++)
      {
        if (anmq[j] == 1)   nenmq++;  // compute number of required non-mechanical quantities on i-th element
      } 

      if (nenmq == 0)
        continue;

      reallocv (RSTCKVEC(nnem,nodval));
      reallocv (RSTCKVEC(nip,ipval));
      memset(enmqo, 0, sizeof(*enmqo)*Mm->tnknmq);      

      for (j=0, k=0; j<Mm->tnknmq; j++)
      {
        if (anmq[j] == 1)
        {
          enmqo[k]= nonmechquant(j+1);  // assemble array with ordering of required non-mechanical quantities on element
          k++;
        }
      } 


      nnet = Tt->give_nne (i); // number of nodes on transport element
      
      //  interpolation of individual non-mechanical(=transport) quantities
      auxnv.n = nnem;
      for (j=0; j<nenmq; j++)
      { 
        // initial nodal values of all required non-mechanical(=transport) quantities on transport element
        // because the order of elements in transport part must be the same or lower then in mechanical part
        // only one function elem_transq_init_nodval_comp may be used
        elem_transq_init_nodval(i, nodval, enmqo[j]);

        if (nnem > nnet) // MEFEL element approximation is quadratic but TRFEL is linear
        {
          // it is assumed that order of corner nodes on a quadratic element 
          // is identical with order on linear one.

          elem_intpointval2 (i, nodval, ipval);
        }
        else
        {
          // MEFEL and TRFEL element approximation are identical =>
          // it is assumed that order of nodes on TRFEL and MEFEL elements are identical

          // TRFEL element approximation is quadratic but MEFEL is linear
          // nothing is necessary because it is assumed that order of corner nodes on a quadratic
          // element is identical with order on linear one.
          
          //  interpolation of nodal values to integration points
          elem_intpointval (i, nodval, ipval);
        }

      
        //  number of the first integration point
        ipid=Mt->elements[i].ipp[0][0];
      
        for (k=0;k<nip;k++)
        {
          Mm->storenonmechq(enmqo[j], ipid, ipval[k]);
          ipid++;
        }
      }
    }
  }
  delete [] enmqo;
}



/**
  The function computes and passes initial nodal values used in the case of passing data
  between TRFEL and MEFEL. Initial nodal values are computed at particular TRFEL element nodes and then 
  passed to MEFEL elements which approximate them to the MEFEL int. points.

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.12.2013
*/
void init_trfel_mefel_by_nodes_comp(void)
{
  long i, j, k, nnet, nnem, nnv, nip, ipid, nenmq;
  ivector nodes;
  vector nodval, ipval, auxnv;
  long anmq[Mm->tnknmq];
  nonmechquant *enmqo; 
  
  for (i=0;i<Mt->ne;i++)
  {
    if (Gtm->leso[i]==1)
    {
      memset(anmq, 0, sizeof(*anmq)*Mm->tnknmq);
      ipid=Mt->elements[i].ipp[0][0];
      for(j=0; j<Mt->give_tnip(i); j++)
      {
        for (k=0; k<Mm->ip[ipid].nm; k++)
          Mm->give_reqnmq(ipid, k, anmq); // search for required non-mechanical quantities on the i-th integration point
        ipid++;
      }

      for (j=0, nenmq=0; j<Mm->tnknmq; j++)
      {
        if (anmq[j] == 1)   nenmq++;  // compute number of required non-mechanical quantities on i-th element
      } 

      if (nenmq == 0)
        continue;

      enmqo = new nonmechquant[nenmq];      
      for (j=0, k=0; j<Mm->tnknmq; j++)
      {
        if (anmq[j] == 1)
        {   
          enmqo[k]= nonmechquant(j+1);  // assemble array with ordering of required non-mechanical quantities on element
          k++;
        }
      } 

      nnet = Tt->give_nne (i); // number of nodes on transport element      
      nnem = Mt->give_nne (i); // number of nodes on mechanical element
      nnv = (nnet < nnem) ? nnet : nnem; // minimum number of nodes
      reallocv (RSTCKVEC(nnv*Mm->nnmq,nodval));
      
      // initial nodal values of all required non-mechanical(=transport) quantities on transport element
      // because the order of elements in transport part must be the same or lower then in mechanical part
      // only one function elem_transq_nodval_comp may be used
      elem_transq_init_nodval_comp(i, nodval, nnv, nenmq, enmqo);

      //  interpolation of individual non-mechanical(=transport) quantities
      auxnv.n = nnem;
      nip = Mt->give_tnip (i); // number of int. points on mechanical element
      reallocv (RSTCKVEC(nip,ipval));
      for (j=0; j<nenmq; j++)
      { 
        auxnv.a = nodval.a + j*nnv;

        if (nnem > nnet) // MEFEL element approximation is quadratic but TRFEL is linear
        {
          // it is assumed that order of corner nodes on a quadratic element 
          // is identical with order on linear one.

          elem_intpointval2 (i, auxnv, ipval);
        }
        else
        {
          // MEFEL and TRFEL element approximation are identical =>
          // it is assumed that order of nodes on TRFEL and MEFEL elements are identical

          // TRFEL element approximation is quadratic but MEFEL is linear
          // nothing is necessary because it is assumed that order of corner nodes on a quadratic
          // element is identical with order on linear one.
          
          //  interpolation of nodal values to integration points
          elem_intpointval (i, auxnv, ipval);
        }

      
        //  number of the first integration point
        ipid=Mt->elements[i].ipp[0][0];
      
        for (k=0;k<nip;k++)
        {
          Mm->storenonmechq(enmqo[j], ipid, ipval[k]);
          ipid++;
        }
      }
      delete [] enmqo;
    }
  }
  auxnv.a = NULL;
}



/**
  The function computes/passes intial coupling data from TRFEL to MEFEL. Data are taken from the auxiliary integration points 
  in TRFEL and stored in MEFEL nonmechq arrays. 

  @param n[in]   - the number of auxiliary points in the mapping array ipm, in this case it should the 
                   number of regular int. points of the MEFEL module.
  @param ipm[in] - integration point mapping array, 
                   ipm[i].ipp < 0 => auxiliary integration point must be used
                   ipm[i].ipp >= 0 => direct mapping to the regulal integration point, 
                                      no computation of initial values is performed, the values are assumed 
                                      to be computed at the main solution procedure of the problem

  @return The function does not return anything but it computes/takes initial nodal values and state variables 
          at auxiliary integration points of TRFEL and passes them in MEFEL nonmechq arrays.

  Created by Tomas Koudelka, 13.06.2018
*/
void init_trfel_mefel_by_aip(long n, ipmap *ipm)
{
  long i, j, ipp, app;
  const long tnknmq = Mm->tnknmq;
  long anmq[tnknmq];
  long nnmq = Mm->nnmq;
  nonmechquant *nmqo = Mm->nmqo;
  intpointst *tmp_ip;
  long *tmp_elip;
  double *tmp_iv;

  // compute actual unknown values and their gradients at TRFEL auxiliary int. points
  aip_initapproximation(n, ipm);
  // compute other required values at TRFEL auxiliary int. points
  Tm->update_aipval();

  //
  // processing of auxiliary int. points given by natural coordinates
  //
  
  // swap regular integration point and auxiliary integration point arrays
  // from now, all functions will work with AUXILIARY integration points 
  tmp_ip = Tm->ip;
  tmp_elip = Tm->elip;
  tmp_iv = Tm->initval;
  Tm->ip = Tm->aip;
  Tm->elip = Tm->elaip;
  Tm->initval = Tm->aip_initval;
  for(i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app >= 0) // passing values from auxiliary int. point 
    {
      // obtain array of indicators of the required nonmechanical quantities
      // the given TRFEL auxiliary int. point = MEFEL regular int. point, 
      // i.e. index i in the array of TRFEL aux. int points = ipp in MEFEL regular int. points
      memset(anmq, 0, sizeof(*anmq)*tnknmq);
      for (j=0; j<Mm->ip[i].nm; j++)
        Mm->give_reqnmq(i, j, anmq);
      for (j=0; j<nnmq; j++)
      {     
        if (anmq[nmqo[j]-1])
        {
          // the quantity is required at the given MEFEL int. point ->
          // -> take value from the auxiliary int. point stored in TRFEL transmat class
          // at this point arrays of regular int. points are swapped with the auxiliary one
          // and therefore all functions of transmat work with aux. int. points
          Mm->nonmechq[j][i] = Tm->givetransq(nmqo[j], app); 
        }
      }
    }
  }
  // swap regular integration point and auxiliary integration point arrays
  // from now, all functions will work with REGULAR integration points 
  Tm->ip = tmp_ip;
  Tm->elip = tmp_elip;
  Tm->initval = tmp_iv;

  //
  // processing of auxiliary int. points mapped to the regular int. points directly
  //
  for(i=0; i<n; i++)
  {
    ipp = ipm[i].ipp;
    // check direct mapping between auxiliary and regular int. point
    if (ipp >= 0)
    {
      // obtain array of indicators of the required nonmechanical quantities
      // the given TRFEL auxiliary int. point = MEFEL regular int. point, 
      // i.e. index i in the array of TRFEL aux. int points = ipp in MEFEL regular int. points
      memset(anmq, 0, sizeof(*anmq)*tnknmq);
      for (j=0; j<Mm->ip[i].nm; j++)
        Mm->give_reqnmq(i, j, anmq);
      for (j=0; j<nnmq; j++)
      {
        if (anmq[nmqo[j]-1])
        {
          // the quantity is required at the given MEFEL int. point ->
          // -> take value from the auxiliary int. point stored in TRFEL transmat class
          // at this point arrays of regular int. points are swapped with the auxiliary one
          // and therefore all functions of transmat work with aux. int. points
          Mm->nonmechq[j][i] = Tm->givetransq(nmqo[j], ipp); 
        }
      }
    }
  }
}



/**
  The function passes transport quantities to MEFEL auxiliary integration points.

  @param n[in]    - the number of auxiliary points in the mapping array ipm, in this case it should the 
                    number of regular int. points of the TRFEL module.
  @param ipm[in]  - integration point mapping array, 
                    ipm[i].ipp < 0 => auxiliary integration point must be used
                    ipm[i].ipp >= 0 => direct mapping to the regulal integration point, 
                                       no computation of values is performed, the values are assumed 
                                       to be computed at the main solution procedure of the problem

  @return The function does not return anything but it actualizes aip_nonmechq array in mechmat class (Mm->aip_nonmechq).

  Created by Tomas Koudelka, 20.06.2018
*/
void actualize_aip_nonmechq(long n, ipmap *ipm)
{
  long i, j, app;
  const long tnknmq = Mm->tnknmq;
  long anmq[tnknmq];
  long nnmq = Mm->nnmq;
  nonmechquant *nmqo = Mm->nmqo;
  intpoints *tmp_ip;
  long *tmp_elip;

  // swap regular integration point auxiliary integration point arrays
  // from now, all functions will work with AUXILIARY integration points 
  tmp_ip = Mm->ip;
  tmp_elip = Mm->elip;
  Mm->ip = Mm->aip;
  Mm->elip = Mm->elaip;

  for(i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app >= 0) // passing values from auxiliary int. point 
    {
      // obtain array of indicators of the required nonmechanical quantities
      // at the given MEFEL auxiliary int. point = TRFEL regular int. point, 
      // i.e. index i in the array of MEFEL auxiliar int. points = ipp in TRFEL regular int. points
      // at this point arrays of regular int. points are swapped with the auxiliary one
      // and therefore all functions of mechmat work with aux. int. points
      memset(anmq, 0, sizeof(*anmq)*tnknmq);
      for (j=0; j<Mm->ip[app].nm; j++)
        Mm->give_reqnmq(app, j, anmq);

      for (j=0; j<nnmq; j++)
      {
        if (anmq[nmqo[j]-1])
        {
          // the quantity is required at the given MEFEL auxiliary int. point ->
          // -> take value from the regular int. point stored in TRFEL transmat class
          Mm->aip_nonmechq[j][app] = Tm->givetransq(nmqo[j], i);
        }
      }
    }
  }

  // swap regular integration point and auxiliary integration point arrays
  // from now, all functions will work with REGULAR integration points 
  Mm->ip = tmp_ip;
  Mm->elip = tmp_elip;
}



/**
  The function passes mechanical quantities to TRFEL auxiliary integration points.

  @param n[in]    - the number of auxiliary points in the mapping array ipm, in this case it should the 
                    number of regular int. points of the MEFEL module.
  @param ipm[in]  - integration point mapping array, 
                    ipm[i].ipp < 0 => auxiliary integration point must be used
                    ipm[i].ipp >= 0 => direct mapping to the regulal integration point, 
                                       no computation of values is performed, the values are assumed 
                                       to be computed at the main solution procedure of the problem

  @return The function does not return anything but it actualizes aip_nontransq array in transmat class (Tm->aip_nontransq).

  Created by Tomas Koudelka, 20.06.2018
*/
void actualize_aip_nontransq(long n, ipmap *ipm)
{
  long i, j, app;
  const long tnkntq = Tm->tnkntq;
  long antq[tnkntq];
  long nntq = Tm->nntq;
  nontransquant *ntqo = Tm->ntqo;
  intpointst *tmp_ip;
  long *tmp_elip;
  double *tmp_iv;

  if (nntq == 0) // partially coupled problem, no nontransport quantities are required in the transport problem
    return;

  // swap regular integration point and auxiliary integration point arrays
  // from now, all functions will work with AUXILIARY integration points 
  tmp_ip = Tm->ip;
  tmp_elip = Tm->elip;
  tmp_iv = Tm->initval;
  Tm->ip = Tm->aip;
  Tm->elip = Tm->elaip;
  Tm->initval = Tm->aip_initval;

  for(i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app >= 0) // passing values from auxiliary int. point 
    {
      // obtain array of indicators of the required nontransport quantities
      // at the given TRFEL auxiliary int. point = MEFEL regular int. point, 
      // i.e. index i in the array of TRFEL auxiliar int. points = ipp in MEFEL regular int. points
      // at this point arrays of regular int. points are swapped with the auxiliary one
      // and therefore all functions of transmat work with aux. int. points
      memset(antq, 0, sizeof(*antq)*tnkntq);
      Tm->give_reqntq(app, antq); // query aux. point for the required nontransport quantities

      for (j=0; j<nntq; j++)
      {
        if (antq[ntqo[j]-1])
        {
          // the quantity is required at the given TRFEL auxiliary int. point ->
          // -> take value from the regular int. point stored in TRFEL transmat class
          Tm->aip_nontransq[j][app] = Mm->givemechq(ntqo[j], i);
        }
      }
    }
  }

  // swap regular integration point and auxiliary integration point arrays
  // from now, all functions will work with REGULAR integration points 
  Tm->ip = tmp_ip;
  Tm->elip = tmp_elip;
  Tm->initval = tmp_iv;
}



/**
  The function creates mapping between MEFEL and TRFEL regular integration points which results
  in initialization of arrays MTipmap and TMipmap connected with regular int. points of one FE mesh 
  with calculated natural coordinates of these int. points on the other FE mesh.
  In the case, that two regular int. points have the same position on both FE meshes,
  direct mapping is used with the help of indeces of regular int. points. 

  For each MEFEL regular int. point, array MTipmap contains corresponding TRFEL element id 
  and calculated natural coordinates or direct mapping to the TRFEL regular int. points.
  
  For each TRFEL regular int. point, array TMipmap contains corresponding MEFEL element id 
  and calculated natural coordinates or direct mapping to the MEFEL regular int. points.

  Arrays MTipmap and TMipmap are used in connection with the auxiliary integration points in MEFEL and TRFEL.

  @param ni[in] - the maximum number of iterations used in the natural coordinate computation
  @param err[in] - tolerance used in itartion solver for the natural coordinate computation
                   or maximum distance at which two points are assumed to be identical


  @return The function does not return anything but initializes arrays TMipmap and MTipmap.

  Created by Tomas Koudelka, 2017
*/
void metr_ip_mapping(long ni, double err)
{
  long tnipt = Tm->tnip;
  long tnipm = Mm->tnip;
  long nem = Mt->ne;
  long net = Tt->ne;
  long ii, jj, i, j, nip, ipp;
  vector coord(ASTCKVEC(3));

  MTipmap = new ipmap[tnipm];
  TMipmap = new ipmap[tnipt];

  for(i=0; i<nem; i++)
  {
    for(ii=0; ii<Mt->give_nb(i); ii++)
    {
      for(jj=0; jj<Mt->give_nb(i); jj++)
      {
        nip = Mt->give_nip(i, ii, jj);
        ipp = Mt->elements[i].ipp[ii][jj];
        for(j=0; j<nip; j++, ipp++)
        {
          ipcoord(i, ipp, ii, jj, coord);
          MTipmap[ipp].x = coord[0];
          MTipmap[ipp].y = coord[1];
          MTipmap[ipp].z = coord[2];
        }
      }
    }
  }
  for (ii=0; ii<Tp->ntm; ii++)
  {
    for (jj=0; jj<Tp->ntm; jj++)
    {
      for(i=0; i<net; i++)
      {
        nip = Tt->give_nip(i, ii, jj);
        ipp = Tt->elements[i].ipp[ii][jj];
        for(j=0; j<nip; j++, ipp++)
        {
          ipcoordt(i, ipp, ii, jj, coord);
          TMipmap[ipp].x = coord[0];
          TMipmap[ipp].y = coord[1];
          TMipmap[ipp].z = coord[2];
        }
      }
      if (Tp->savemode==1)  break;
    }
    if (Tp->savemode==1)  break;
  }

  mefel_trfel_ip_mapping(ni, err);
  trfel_mefel_ip_mapping(ni, err);

  for (i=0; i< tnipm; i++)
  {
    if (MTipmap[i].cerr > err)
    {
      print_err("natural coordinates are out of tolerance on MEFEL element %ld, ip=%ld\n"
                " xi=%le, eta=%le, zeta=%le, err=%le > %le, x=%le, y=%le, z=%le\n", 
                __FILE__, __LINE__, __func__, Mm->elip[i]+1, i+1, MTipmap[i].xi, MTipmap[i].eta, MTipmap[i].zeta,
                MTipmap[i].cerr, err, MTipmap[i].x, MTipmap[i].y, MTipmap[i].z);
      abort();
    }
    //fprintf(Out, "MT aux. point %4ld: eid=%4ld, xi=% le, eta=% le, zeta=% le, x=% le, y=% le, z=% le\n", 
    //        i+1, Mm->elip[i]+1, MTipmap[i].xi, MTipmap[i].eta, MTipmap[i].zeta, MTipmap[i].x, MTipmap[i].y, MTipmap[i].z);
  }  

  for (i=0; i< tnipt; i++)
  {
    if (TMipmap[i].cerr > err)
    {
      print_err("natural coordinates are out of tolerance on TRFEL element %ld, ip=%ld\n"
                " xi=%le, eta=%le, zeta=%le, err=%le > %le,\n x=%le, y=%le, z=%le\n", 
                __FILE__, __LINE__, __func__, Tm->elip[i]+1, i+1, TMipmap[i].xi, TMipmap[i].eta, TMipmap[i].zeta,
                TMipmap[i].cerr, err, TMipmap[i].x, TMipmap[i].y, TMipmap[i].z);
      abort();
    }
    //fprintf(Outt, "TM aux. point %4ld: eid=%4ld, xi=% le, eta=% le, zeta=% le, x=% le, y=% le, z=% le\n", 
    //        i+1, Tm->elip[i]+1, TMipmap[i].xi, TMipmap[i].eta, TMipmap[i].zeta, TMipmap[i].x, TMipmap[i].y, TMipmap[i].z);
  }  
}



/**
  The function creates mapping between TRFEL and MEFEL regular integration points which results
  in initialization of array TMipmap connected with regular int. points of TRFEL FE mesh 
  with calculated natural coordinates of these int. points on the MEFEL FE mesh.
  In the case, that two regular int. points have the same position on both FE meshes,
  direct mapping is used with the help of indeces of regular int. points. 

  For each TRFEL regular int. point, array TMipmap contains corresponding MEFEL element id 
  and calculated natural coordinates or direct mapping to the MEFEL regular int. points.

  Array TMipmap is used in connection with the auxiliary integration points in MEFEL.

  @param ni[in] - the maximum number of iterations used in the natural coordinate computation
  @param err[in] - tolerance used in itartion solver for the natural coordinate computation
                   or maximum distance at which two points are assumed to be identical


  @return The function does not return anything but initializes arrays TMipmap.

  Created by Tomas Koudelka, 2017
*/
void trfel_mefel_ip_mapping(long ni, double err)
{
  long nem = Mt->ne;
  long i, j, nne, eldim;
  vector cg(ASTCKVEC(3));
  vector x, y, z;
#ifndef  INC_OPENMP
  vector p(ASTCKVEC(3));
#endif  
  double d2, r2, xi, eta, zeta, cerr;
  long tnipt = Tm->tnip;
  long nrerr;
  char *ips = new char[tnipt];
  memset(ips, 0, sizeof(*ips)*tnipt);

  if (Mesprc == 1)
    fprintf(stdout, " Calculation of %ld auxiliary int. points TRFEL->MEFEL:\n", tnipt);
  if (Tm->nntq == 0){
    fprintf(stdout, " No mechanical quantities are required in TRFEL => no auxiliary points are calculated\n");
    for(j=0; j<tnipt; j++){
      TMipmap[j].cerr = err*err/10.0;
      TMipmap[j].ipp = 0;
      TMipmap[j].eid = 0;
    }
    return;
  }
  for(i=0; i<nem; i++)
  {
    if ((Mesprc == 1) && ((i%1000 == 0) || (i == nem-1))){
      fprintf(stdout, " Searching aux. points on %ld. element\r", i+1);
      fflush(stdout);
    }
    if(Gtm->leso[i])
    {
      centroid(i, cg);
      r2 = Gtm->max_sqrdist_nod_pt(i, cg);
#ifdef INC_OPENMP
    #pragma omp parallel        \
            num_threads(Numth), \
            private(j, nne, eldim, x, y, z, d2, xi, eta, zeta, cerr, nrerr)
     {
      vector p(ASTCKVEC(3));
      #pragma omp for 
#endif      
      for(j=0; j<tnipt; j++)
      {
        if (ips[j])
          continue;
        d2  = sqr(TMipmap[j].x-cg(0)) + sqr(TMipmap[j].y-cg(1)) + sqr(TMipmap[j].z-cg(2));
        if (d2 <= r2)
        {
          TMipmap[j].ipp = Mt->give_closest_ip_ncoord(i, TMipmap[j].x, TMipmap[j].y, TMipmap[j].z, MTipmap, err, xi, eta, zeta);
          if (TMipmap[j].ipp >= 0)  // the closest integration point on element matches the global coordinates of TMipmap[j] exactly
          {                      // => no further investigation is needed
            TMipmap[j].xi = xi;
            TMipmap[j].eta = eta;
            TMipmap[j].zeta = zeta;
            TMipmap[j].cerr = 0.0;
            ips[j] = 1;
            continue;
          }
          eldim = Mt->give_dimension(i);
          nne = Mt->give_nne(i);
          reallocv(RSTCKVEC(nne, x));
          reallocv(RSTCKVEC(nne, y));
          reallocv(RSTCKVEC(nne, z));
          switch (eldim)
          {
            case 1:
              Mt->give_node_coord3d(x, y, z, i);
              xi = eta  = zeta = 0.0;
              nrerr = point_natcoord_1d_3d(TMipmap[j].x, TMipmap[j].y, TMipmap[j].z, x, y, z, ni, err, xi, cerr);
              corr_nat_coord_bounds(Gtm->give_siftopelem_type(i), xi, eta, zeta);
              bf_1d_3d(p, x, y, z, xi);
              cerr = length(p, TMipmap[j].x, TMipmap[j].y, TMipmap[j].z);
              break;
            case 2:
              Mt->give_node_coord2d(x, y, i);
              xi = eta  = zeta = 0.0;
              nrerr = point_natcoord_2d(TMipmap[j].x, TMipmap[j].y, x, y, ni, err, xi, eta, cerr);
              corr_nat_coord_bounds(Gtm->give_siftopelem_type(i), xi, eta, zeta);
              bf_2d(p, x, y, xi, eta);
              cerr = length(p, TMipmap[j].x, TMipmap[j].y, TMipmap[j].z);
              break;
            case 3:
              Mt->give_node_coord3d(x, y, z, i);
              xi = eta  = zeta = 0.0;
              nrerr = point_natcoord_3d(TMipmap[j].x, TMipmap[j].y, TMipmap[j].z, x, y, z, ni, err, xi, eta, zeta, cerr);
              corr_nat_coord_bounds(Gtm->give_siftopelem_type(i), xi, eta, zeta);
              bf_3d(p, x, y, z, xi, eta, zeta);
              cerr = length(p, TMipmap[j].x, TMipmap[j].y, TMipmap[j].z);
              break;
            default:
              print_err("unknown dimension (dim=%ld) of the element %ld is required", __FILE__, __LINE__, __func__, eldim, i+1);
              abort();
          }

          if (nrerr > ni)
            fprintf(Outt,"solution of ip=%ld natural coordinates of the element %ld did not converge in %ld steps (file:%s, line:%d, %s)\n", j, Tm->elip[j], nrerr, __FILE__, __LINE__, __func__);

          TMipmap[j].xi = xi;
          TMipmap[j].eta = eta;
          TMipmap[j].zeta = zeta;
          TMipmap[j].cerr = cerr;
          TMipmap[j].eid = i;
          if (cerr <= err)
            ips[j] = 1;
        }
      }
#ifdef INC_OPENMP
     }// end of #pragma block
#endif      
    }
  }
  fprintf(stdout, "\n");
  delete [] ips;
}



/**
  The function creates mapping between MEFEL and TRFEL regular integration points which results
  in initialization of array MTipmap connected with regular int. points of MEFEL FE mesh 
  with calculated natural coordinates of these int. points on the TRFEL FE mesh.
  In the case, that two regular int. points have the same position on both FE meshes,
  direct mapping is used with the help of indeces of regular int. points. 

  For each MEFEL regular int. point, array MTipmap contains corresponding TRFEL element id 
  and calculated natural coordinates or direct mapping to the TRFEL regular int. points.

  Array MTipmap is used in connection with the auxiliary integration points in TRFEL.

  @param ni[in] - the maximum number of iterations used in the natural coordinate computation
  @param err[in] - tolerance used in itartion solver for the natural coordinate computation
                   or maximum distance at which two points are assumed to be identical


  @return The function does not return anything but initializes arrays MTipmap.

  Created by Tomas Koudelka, 2017
*/
void mefel_trfel_ip_mapping(long ni, double err)
{
  long net = Tt->ne;
  long i, j, eldim, nne;
  vector cg(ASTCKVEC(3));
  vector x, y, z;
#ifndef INC_OPENMP
  vector p(ASTCKVEC(3));
#endif
  double d2, r2, xi, eta, zeta, cerr;
  long tnipm = Mm->tnip;
  long nrerr;
  char *ips = new char[tnipm];
  memset(ips, 0, sizeof(*ips)*tnipm);

  if (Mesprc == 1)
    fprintf(stdout, "\n Calculation of %ld auxiliary int. points MEFEL->TRFEL:\n", tnipm);
  if (Mm->nnmq == 0){
    fprintf(stdout, " No transport quantities are required in MEFEL => no auxiliary points are calculated\n");
    for(j=0; j<tnipm; j++){
      MTipmap[j].cerr = err*err/10.0;
      MTipmap[j].ipp = 0;
      MTipmap[j].eid = 0;
    }
    return;
  }
  for(i=0; i<net; i++)
  {
    if ((Mesprc == 1) && ((i%1000 == 0) || (i == net-1))){
      fprintf(stdout, " Searching aux. points on %ld. element\r", i+1);
      fflush(stdout);
    }
    
    if(Gtt->leso[i])
    {
      centroidt(i, cg);
      r2 = Gtt->max_sqrdist_nod_pt(i, cg);
#ifdef INC_OPENMP
     #pragma omp parallel        \
             num_threads(Numth), \
             private(j, nne, eldim, x, y, z, d2, xi, eta, zeta, cerr, nrerr)
     {
       vector p(ASTCKVEC(3));
      #pragma omp for      
#endif      
      for(j=0; j<tnipm; j++)
      {
        if (ips[j])
          continue;
        d2  = sqr(MTipmap[j].x-cg(0)) + sqr(MTipmap[j].y-cg(1)) + sqr(MTipmap[j].z-cg(2));
        if (d2 <= r2)
        {
          eldim = Tt->give_dimension(i);
          nne = Tt->give_nne(i);
          reallocv(RSTCKVEC(nne, x));
          reallocv(RSTCKVEC(nne, y));
          reallocv(RSTCKVEC(nne, z));

          MTipmap[j].ipp = Tt->give_closest_ip_ncoord(i, MTipmap[j].x, MTipmap[j].y, MTipmap[j].z, TMipmap, err, xi, eta, zeta);
          if (MTipmap[j].ipp >= 0)  // the closest integration point on element matches th global coordinates of TMipmap[j] exactly
          {                      // => no further investigation is needed
            ips[j] = 1;
            MTipmap[j].xi = xi;
            MTipmap[j].eta = eta;
            MTipmap[j].zeta = zeta;
            MTipmap[j].cerr = 0.0;
            continue;
          }
          switch (eldim)
          {
            case 1:
              Tt->give_node_coord3d(x, y, z, i);
              xi = eta  = zeta = 0.0;
              nrerr = point_natcoord_1d_3d(MTipmap[j].x, MTipmap[j].y, MTipmap[j].z, x, y, z, ni, err, xi, cerr);
              corr_nat_coord_bounds(Gtt->give_siftopelem_type(i), xi, eta, zeta);
              bf_1d_3d(p, x, y, z, xi);
              cerr = length(p, MTipmap[j].x, MTipmap[j].y, MTipmap[j].z);
              break;
            case 2:
              Tt->give_node_coord2d(x, y, i);
              xi = eta  = zeta = 0.0;
              nrerr = point_natcoord_2d(MTipmap[j].x, MTipmap[j].y, x, y, ni, err, xi, eta, cerr);
              corr_nat_coord_bounds(Gtt->give_siftopelem_type(i), xi, eta, zeta);
              bf_2d(p, x, y, xi, eta);
              cerr = length(p, MTipmap[j].x, MTipmap[j].y, MTipmap[j].z);
              break;
            case 3:
              Tt->give_node_coord3d(x, y, z, i);
              xi = eta  = zeta = 0.0;
              nrerr = point_natcoord_3d(MTipmap[j].x, MTipmap[j].y, MTipmap[j].z, x, y, z, ni, err, xi, eta, zeta, cerr);
              corr_nat_coord_bounds(Gtt->give_siftopelem_type(i), xi, eta, zeta);
              bf_3d(p, x, y, z, xi, eta, zeta);
              cerr = length(p, MTipmap[j].x, MTipmap[j].y, MTipmap[j].z);
              break;
            default:
              print_err("unknown dimension (dim=%ld) of the element %ld is required", __FILE__, __LINE__, __func__, eldim, i+1);
              abort();
          }

          if (nrerr > ni)
            fprintf(Out,"solution of ip=%ld natural coordinates of the element %ld did not converge in %ld steps (file:%s, line:%d, %s)\n", j, Mm->elip[j], nrerr, __FILE__, __LINE__, __func__);

          MTipmap[j].xi = xi;
          MTipmap[j].eta = eta;
          MTipmap[j].zeta = zeta;
          MTipmap[j].cerr = cerr;
          MTipmap[j].eid = i;
          if (cerr <= err)
            ips[j] = 1;
        }
      }
#ifdef INC_OPENMP
     }// end of #pragma block
#endif      
    }
  }
  fprintf(stdout, "\n");
  delete [] ips;
}





/**
  The function approximates transport quantity from integration points to nodes.
   
  @param gv[out]    - allocated %vector of nodal nontransport quantity
  @param nodmap[in] - array of TRFEL->MEFEL nodal correspondence map, nodmap[i] = MEFEL node id of TRFEL i-th node
  @param nmq[in]    - type of transport (non-mechanical) quantity
 
  Created by TKo 04.2018
*/
void trfmef_give_transq_nodval (double *gv, long *nodmap, nonmechquant nmq)
{
  long i, j, nne;
  ivector nodes;
  vector nodval;
  
  for (i=0;i<Tt->ne;i++){
    if (Gtt->leso[i]==1){
      
      nne = Tt->give_nne (i);
      
      reallocv (RSTCKIVEC(nne,nodes));
      reallocv (RSTCKVEC(nne,nodval));
      
      //  nodes on element
      Tt->give_elemnodes (i,nodes);
      
      //  nodal values on element
      elem_transq_nodval (i, nodval, nmq);
            
      for (j=0;j<nne;j++){
	gv[nodmap[nodes[j]]] = nodval[j];
      }
    }
  }
}



/**
  The function creates map between TRFEL nodes and MEFEL nodes according to 
  TRFEL and MEFEL elements. The number of MEFEL and TRFEL elements must be the same or
  MEFEL mesh must have the first Tt->ne of elements same geometrical shape as in TRFEL
  The order of nodes on the first Tt->ne elements must be the same in TRFEL and MEFEL, but
  MEFEL elements may have greater number of nodes on elements if the positions of first 
  element nodes corresponds to positions on TRFEL elements.

  @param nodmap[out] - pointer to allocated array, which will be filled with TRFEL-->MEFEL nodal correspondence map
                       nodmap[i] = MEFEL node number of TRFEL i-th node

  @return The function does not return anything but it stores nodal correspondence map in the argument nodmap.

  Created by Tomas Koudelka, 04.2018
*/
void create_trfmef_nod_map(long *nodmap)
{
  long nnt = Tt->nn;
  long nnm = Mt->nn;
  long ne;
  long i, j, nnem, nnet;
  ivector enodm, enodt;

  if (Cp->tprob==fully_coupled_mech_trans)
    ne = Gtu->ne;
  else
  {
    ne = Gtm->ne;
    if (ne > Gtt->ne)
      ne = Gtt->ne;
  }

  // initialize TRFEL->MEFEL nodal map
  for (i=0; i<nnt; i++)
    nodmap[i] = -1;

  if (nnt == nnm) // the same number of nodes in TRFEL and MEFEL meshes => direct correspondence of nodes
  {
    for(i=0; i<nnt; i++)
      nodmap[i] = i;
  }
  else // the number of nodes in TRFEL and MEFEL meshes differs => node correspondence is established according to the minimum number of elements in meshes
  {
    for (i=0; i<ne; i++)
    {
      // get number of nodes on MEFEL and TRFEL i-th element
      nnem = Mt->give_nne(i);
      nnet = Tt->give_nne(i);
      if (nnet <= nnem)
      {
        reallocv(RSTCKIVEC(nnem, enodm));
        reallocv(RSTCKIVEC(nnet, enodt));

        // get elelment nodes on MEFEL and TRFEL i-th element
        Mt->give_elemnodes(i, enodm);
        Tt->give_elemnodes(i, enodt);

        // Elements in TRFEL and MEFEL meshes must be ordered in the same way.
        // Element nodes on TRFEL element and MEFEL element must be ordered in the same way,
        // but the number of nodes on MEFEL element may be higher than on TRFEL element
        for (j=0; j<nnet; j++)
          nodmap[enodt[j]] = enodm[j]; // create map between j-th TRFEL element node j-th MEFEL element node
      }
      else
      {
        print_err("incompatible number of nodes on element %ld (nne_trf=%ld > nne_mef=%ld)", __FILE__, __LINE__, __func__, i+1, nnet, nnem);
        abort();
      }
    }
    // check whether all TRFEL nodes were mapped on the MEFEL ones if the number of nodes differs
    for (i=0; i<nnt; i++)
    {
      if (nodmap[i] < 0)
      {
        print_err("TRFEL->MEFEL nodal correspondence cannot be established at node %ld", __FILE__, __LINE__, __func__, i+1);
        abort();
      }
    }
  }
}









/*******************************************************************
  The functions listed below are obsolete
 ********************************************************************/


/**
  The function exports initial data of TRFEL into MEFEL. Actualy, only initial temperature (known in TRFEL) 
  is recalculated to MEFEL integration points and stored in Mm->nonmechq if it is required.
  TRFEL initial nodal tempeartures are taken from the global %vector of initial nodal values directly.
  This function cannot handle cases of nodes at the interface between different material models and in the case of 
  retrieving data from closest int. point, the initial nodal values are rewrriten in course of element loop in trfmef_give_transq_nodval 
  and thus the same initial nodal value is not guaranteed for different element numbering!
   
  05/06/2012, TKr
  Modified by TKo 10.10.2013
  Renamed by TKo 12.6.2018, original name was init_trfel_mefel, now the function is substituted by the more general function init_trfel_mefel_by_nodes
*/
void init_trfel_mefel_orig ()
{
  long i, nnm;
  double *gv;
  
  //  number of nodes in MEFEL
  nnm = Mt->nn;
  
  gv = new double [nnm];
  for (i=0; i<Mm->nnmq; i++)
  {
    if (Mm->nmqo[i] == initial_temperature)
    {
      // initial temperature should be passed once at the computation begining
      nullv (gv,nnm);

      //  selection or interpolation of quantity
      nodal_nodal_values (gv, TM_nod_map, Mm->nmqo[i]); // ???!!! fci casem asi presunout do TRFELu (jestli to pujde)

      if (Cp->bb==quad_quad || Cp->bb==lin_lin){
        intpointval (gv,Mm->nmqo[i],1.0);
      }
      if (Cp->bb==quad_lin){
        intpointval2 (gv,Mm->nmqo[i]);
      }
    }
  }

  delete [] gv;
}



/**
   function exports data from MEFEL into TRFEL
   
   05/06/2012, TKr
   Modified by TKo 10.10.2013
   Renamed by TKo 12.6.2018, original name was mefel_trfel, now the function is substituted by the more general function mefel_trfel_by_nodes
*/
void mefel_trfel_orig ()
{
  long i, j, k, nne, nip, ipid;
  ivector nodes;
  vector nodval, ipval;
  
  //  interpolation of individual non-transport (mechanical) quantities
  for(i=0; i<Tm->nntq; i++)
  {
    for (j=0;j<Tt->ne;j++)
    {
      if (Gtt->leso[j]==1)
      {
        nne = Tt->give_nne (j);
        nip = Tt->give_tnip (j);
      
        reallocv (RSTCKIVEC(nne,nodes));
        reallocv (RSTCKVEC(nne,nodval));
        reallocv (RSTCKVEC(nip,ipval));
      
        //  nodes on element
        Tt->give_elemnodes (j, nodes);
      
	if (Cp->bb==quad_quad || Cp->bb==lin_lin){
	  //  nodal values on mechanical element
	  elem_mechq_nodval (j, nodval, Tm->ntqo[i]);
	}
	if (Cp->bb==quad_lin){
	  //  nodal values on mechanical element
	  elem_mechq_nodval2 (j, nodval, Tm->ntqo[i]);
	}

	//  interpolation of nodal values to integration points
	elem_intpointvalt (j, nodval, ipval);

        //  number of the first integration point
        ipid=Tt->elements[j].ipp[0][0];
      
        for (k=0;k<nip;k++)
        {
          Tm->storenontransq(Tm->ntqo[i], ipid, ipval[k]);
          ipid++;
        }
      }
    }
  }
}



/**
  The function exports data of TRFEL into MEFEL. All required non-mechanical quantities (known in TRFEL) are
  recalculated to MEFEL integration points and stored in Mm->nonmechq. TRFEL nodal values are taken 
  either from the global %vector of nodal values directly or from the closest int. point if they are known only at the
  material model level.
  This function cannot handle cases of nodes at the interface between different material models and in the case of 
  retrieving data from closest int. point, the nodal values are rewrriten in course of element loop in trfmef_give_transq_nodval 
  and thus the same nodal value is not guaranteed for different element numbering!
   
   
  05/06/2012, TKr
  Modified by TKo 10.10.2013
  Renamed by TKo 12.6.2018, original name was trfel_mefel, now the function is substituted by the more general function trfel_mefel_by_nodes
*/
void trfel_mefel_orig ()
{
  long i, nnm;
  double *gv;
 
  //  number of nodes in MEFEL
  nnm = Mt->nn;
  
  gv = new double [nnm];
  for (i=0; i<Mm->nnmq; i++)
  {
    if (Mm->nmqo[i] == initial_temperature)
    {
      // initial temperature can be passed only once in init_trfel_mefel()
      continue;
    }

    nullv (gv,nnm);
    //  selection or interpolation of quantity
    nodal_nodal_values (gv, TM_nod_map, Mm->nmqo[i]); // ???!!! fci casem asi presunout do TRFELu (jestli to pujde)

    if (Cp->bb==quad_quad || Cp->bb==lin_lin){
      intpointval(gv,Mm->nmqo[i],1.0);
    }
    if (Cp->bb==quad_lin){
      intpointval2(gv,Mm->nmqo[i]);
    }
  }

  delete [] gv;
}









/***********************************************************************
  The functions listed below are obsolete but used in some old solvers
*************************************************************************/



/**
   function approximates initial temperature computed in the TRFEL part into the MEFEL part
   
   21.11.2007, JK
*/
void approximation_inittemper ()
{
  long nnm;
  double *gv;
  
  //  number of nodes in MEFEL
  nnm = Mt->nn;
  
  gv = new double [nnm];
  nullv (gv,nnm);
  
  //  selection or interpolation of temperatures 
  nodal_nodal_values (gv, TM_nod_map, initial_temperature);
  
  if (Cp->bb==quad_quad || Cp->bb==lin_lin){
    intpointval (gv,initial_temperature,1.0);
  }
  if (Cp->bb==quad_lin){
    intpointval2 (gv,initial_temperature);
  }
  
  delete [] gv;
}



/**
   function approximates temperature computed in the TRFEL part into the MEFEL part
   
   21.6.2004, JK
*/
void approximation_temper ()
{
  long nnm;
  double *gv;
  
  //  number of nodes in MEFEL
  nnm = Mt->nn;
  
  gv = new double [nnm];
  nullv (gv,nnm);
  
  //  selection or interpolation of temperatures 
  nodal_nodal_values (gv, TM_nod_map, temperature);
  
  //  approximation of nodal data to the data at integration points
  if (Cp->bb==quad_lin){
    intpointval2 (gv,temperature);
  }
  if (Cp->bb==quad_quad || Cp->bb==lin_lin){
    intpointval (gv,temperature,1.0);
  }
  
  /* 
     fprintf (Outc,"\n\n kontrola hodnot na integracnich bodech \n");
     for (long i=0;i<Mm->tnip;i++){
     fprintf (Outc,"\n point %5ld   %le",i,Mm->tempr[i]);
     }
     fprintf (Outc,"\n");
  */

  delete [] gv;
}



/**
   function approximates humidity computed in the TRFEL part into the MEFEL part
   
   22.6.2004, JK
*/
void approximation_humid ()
{
  long nnm;
  double *gv;
  
  //  number of nodes in MEFEL
  nnm = Mt->nn;
  
  gv = new double [nnm];
  nullv (gv,nnm);
  
  //  selection or interpolation of humidities
  nodal_nodal_values (gv, TM_nod_map, rel_hum);
  
  //  approximation of nodal data to the data at integration points
  if (Cp->bb==quad_lin){
    intpointval2 (gv,rel_hum);
  }
  if (Cp->bb==quad_quad || Cp->bb==lin_lin){
    intpointval (gv,rel_hum,1.0);
  }
  
  delete [] gv;
}



/**
   function returns nodal values of selected quantity from the TRFEL part
   
   @param gv - array of nodal quantity, the number of components is equal to the number of nodes in MEFEL
   @param nodmap - array of TRFEL->MEFEL nodal correspondence map, nodmap[i] = MEFEL node id of TRFEL i-th node
   @param nmq - name of the value needed   

   24/8/2012, TKr
   Modified by TKo, 04.2018
*/
void nodal_nodal_values (double *gv, long *nodmap, nonmechquant nmq)
{
  long i,nnt,cn, idof;
  
  nnt = Tt->nn;

  switch (nmq){
    ////////////////////////////
  case temperature:{
    idof = Tp->give_var_dofid(trf_temperature);
    if (idof < 0)
    {
      print_err("temperature is not defined in the problem", __FILE__, __LINE__, __func__);
      abort();
    }
    for (i=0;i<nnt;i++)
    {
      cn=Tt->give_dof (i,idof);
      if (cn>0)  gv[nodmap[i]] = Lsrst->lhsi[cn-1]+Lsrst->lhs[cn-1];
      if (cn<0)  gv[nodmap[i]] = Tb->lc[0].pv[0-cn-1].getval(Tp->time);
    }
    break;
  }
    ////////////////////////////
  case initial_temperature:{
    idof = Tp->give_var_dofid(trf_temperature);
    if (idof < 0)
    {
      printf ("\n\n\n idof = %ld \n\n\n",idof);
      print_err("initial temperature is not defined in the problem", __FILE__, __LINE__, __func__);
      abort();
    }
    for (i=0;i<nnt;i++)
    {
      cn=Tt->give_dof (i,idof);
      if (cn>0)  gv[nodmap[i]] = Lsrst->lhsi[cn-1];
      if (cn<0)  gv[nodmap[i]] = Tb->lc[0].pv[0-cn-1].ipv;
    }
    break;
  }
    ////////////////////////////
  case water_pressure:{
    idof = Tp->give_var_dofid(trf_water_press);
    if (idof < 0)
    {
      print_err("water pressure is not defined in the problem", __FILE__, __LINE__, __func__);
      abort();
    }
    for (i=0;i<nnt;i++)
    {
      cn=Tt->give_dof (i,idof);
      if (cn>0)  gv[nodmap[i]] = Lsrst->lhsi[cn-1]+Lsrst->lhs[cn-1];
      if (cn<0)  gv[nodmap[i]] = Tb->lc[0].pv[0-cn-1].getval(Tp->time);
    }
    break;
  }
  case saturation_degree:{
    //approximation degree of saturation from int. points to nodes
    trfmef_give_transq_nodval (gv, nodmap, nmq);
    break;
  }
    ////////////////////////////
  case suction:{
    //approximation suction from int. points to nodes
    trfmef_give_transq_nodval (gv, nodmap, nmq);
    break;
  }
    ////////////////////////////
  case rel_hum:{
    give_nodal_humid(gv, nodmap);
    break;
  }
  case pore_pressure:{
    // water pressure + gass pressure + crystallization pressure from integration point
    // particular contributions should be scaled by degree of saturation depending 
    // on the used transport material model
    trfmef_give_transq_nodval (gv, nodmap, nmq);
    break;
  }
  case vol_moist_cont:{
    trfmef_give_transq_nodval(gv, nodmap, nmq);
    break;
  }
    ////////////////////////////
  default:{
    print_err("unknown type of nmq %d is required",__FILE__,__LINE__,__func__, nmq);
  }
  }  
}



/**
  The function assembles map of correspondece between 
  TRFEL and MEFEL integration points. The meshes must be identical and 
  integration order of stiffness matrix and conductivity matrix must be 
  the same. Resulting maps are stored in arrays tm_ip and mt_ip.

  @param tm_ip - trfel to mefel ip correspondence map - tm_ip[trf_ipp] = mef_ipp
  @param mt_ip - mefel to trfel ip correspondence map - mt_ip[mef_ipp] = trf_ipp

  @return The function returns maps in the tm_ip and mt_ip.

  Created by Tomas Koudelka, 11.11.2013
*/
void mefel_trfel_ip_mapping_old(long *tm_ip, long *mt_ip)
{
  elemtypet ett;
  long ippt, ippm;
  long iokm;
  long i, j, ii, jj, n;


  //temporarily commented for Temelin computation
  /* if (Mt->ne != Tt->ne)
     {
     print_err("number of elements of TRFEL(=%ld) and MEFEL(=%ld)\n meshes are not the same", 
     __FILE__, __LINE__, __func__, Tt->ne, Mt->ne);
     abort();
     }
  */

  // TRFEL -> MEFEL ip correspondence
  // tm_ip[trf_ipp] = mefel_ipp
  for(i=0; i<Tm->tnip; i++)
    tm_ip[i] = -1;
  // MEFEL -> TRFEL ip correspondence
  // mt_ip[mef_ipp] = trfel_ipp
  for(i=0; i<Mm->tnip; i++)
    mt_ip[i] = -1;

  for(i=0; i<Tt->ne; i++)
  {
    ett = Tt->give_elem_type(i);
    switch (ett)
    {
      case barlint:
      case barlint3d:
        for (ii=0;ii<Tp->ntm;ii++)
        {
          for (jj=0;jj<Tp->ntm;jj++)
          {
            ippt = Tt->elements[i].ipp[ii][jj];
            n  = Tt->give_nip(i, ii, jj);
            for(j=0; j < n; j++)
            {
              tm_ip[ippt] = Mt->elements[i].ipp[0][0];
              ippt++;
            }
          }
        }
        ippm = Mt->elements[i].ipp[0][0];
        for(j=0; j<Mt->give_tnip(i); j++)
        {
          mt_ip[ippm] = Tt->elements[i].ipp[0][0]+1;
          ippm++;
        }
        break;
      case trlint:
      case trlaxisym:
        for (ii=0;ii<Tp->ntm;ii++)
        {
          for (jj=0;jj<Tp->ntm;jj++)
          {
            ippt = Tt->elements[i].ipp[ii][jj];
            n  = Tt->give_nip(i, ii, jj);
            for(j=0; j < n; j++)
            {
              tm_ip[ippt] = Mt->elements[i].ipp[0][0];
              ippt++;
            }
          }
        }
        ippm = Mt->elements[i].ipp[0][0];
        for(j=0; j<Mt->give_tnip(i); j++)
        {
          mt_ip[ippm] = Tt->elements[i].ipp[0][0];
          ippm++;
        }
        break;
      case quadlint:
        for (ii=0;ii<Tp->ntm;ii++)
        {
          for (jj=0;jj<Tp->ntm;jj++)
          {
            ippt = Tt->elements[i].ipp[ii][jj];
            n  = Tt->give_nip(i, ii, jj);
            iokm = Tt->give_intordkm(i, ii, jj);
            iokm *= iokm;
            for(j=0; j < n; j++)
            {
              tm_ip[ippt] = Mt->elements[i].ipp[0][0]+(j%iokm);
              ippt++;
            }
          }
        }
        ippm = Mt->elements[i].ipp[0][0];
        n = Mt->give_nip(i, 0, 0);
        for(j=0; j<n; j++)
        {
          mt_ip[ippm] = Tt->elements[i].ipp[0][0]+j;
          ippm++;
        }
        break;
      case quadlaxisym:
        for (ii=0;ii<Tp->ntm;ii++)
        {
          for (jj=0;jj<Tp->ntm;jj++)
          {
            ippt = Tt->elements[i].ipp[ii][jj];
            n  = Tt->give_nip(i, ii, jj);
            iokm = Tt->give_intordkm(i, ii, jj);
            iokm *= iokm;
            for(j=0; j < n; j++)
            {
              tm_ip[ippt] = Mt->elements[i].ipp[0][0]+(j%iokm);
              ippt++;
            }
          }
        }
        ippm = Mt->elements[i].ipp[0][0];
        n = Mt->give_nip(i, 0, 0);
        for(j=0; j<n; j++)
        {
          mt_ip[ippm] = Tt->elements[i].ipp[0][0]+j;
          ippm++;
        }
        break;
      case ifacequadel:
        for (ii=0;ii<Tp->ntm;ii++)
        {
          for (jj=0;jj<Tp->ntm;jj++)
          {
            ippt = Tt->elements[i].ipp[ii][jj];
            n  = Tt->give_nip(i, ii, jj);
            for(j=0; j < n; j++)
            {
              tm_ip[ippt] = Mt->elements[i].ipp[0][0]+j;
              ippt++;
            }
          }
        }
        ippm = Mt->elements[i].ipp[0][0];
        for(j=0; j<Mt->give_tnip(i); j++)
        {
          mt_ip[ippm] = Tt->elements[i].ipp[0][0]+j;
          ippm++;
        }
        break;
      case barquadt:
        for (ii=0;ii<Tp->ntm;ii++)
        {
          for (jj=0;jj<Tp->ntm;jj++)
          {
            ippt = Tt->elements[i].ipp[ii][jj];
            n  = Tt->give_nip(i, ii, jj);
            iokm = Tt->give_intordkm(i, ii, jj);
            for(j=0; j < n; j++)
            {
              tm_ip[ippt] = Mt->elements[i].ipp[0][0]+(j%iokm);
              ippt++;
            }
          }
        }
        ippm = Mt->elements[i].ipp[0][0];
        for(j=0; j<Mt->give_tnip(i); j++)
        {
          mt_ip[ippm] = Tt->elements[i].ipp[0][0]+j;
          ippm++;
        }
        break;
      case quadquadt:
      case quadquadtax:
        for (ii=0;ii<Tp->ntm;ii++)
        {
          for (jj=0;jj<Tp->ntm;jj++)
          {
            ippt = Tt->elements[i].ipp[ii][jj];
            n  = Tt->give_nip(i, ii, jj);
            iokm = Tt->give_intordkm(i, ii, jj);
            iokm *= iokm;
            for(j=0; j < n; j++)
            {
              tm_ip[ippt] = Mt->elements[i].ipp[0][0]+(j%iokm);
              ippt++;
            }
          }
        }
        ippm = Mt->elements[i].ipp[0][0];
        for(j=0; j<Mt->give_tnip(i); j++)
        {
          mt_ip[ippm] = Tt->elements[i].ipp[0][0]+j;
          ippm++;
        }
        break;
      case lineartett:
      case linearhext:
      case quadratichext:
        for (ii=0;ii<Tp->ntm;ii++)
        {
          for (jj=0;jj<Tp->ntm;jj++)
          {
            ippt = Tt->elements[i].ipp[ii][jj];
            n  = Tt->give_nip(i, ii, jj);
            iokm = Tt->give_intordkm(i, ii, jj);
            iokm *= iokm*iokm;
            for(j=0; j < n; j++)
            {
              tm_ip[ippt] = Mt->elements[i].ipp[0][0]+(j%iokm);
              ippt++;
            }
          }
        }
        ippm = Mt->elements[i].ipp[0][0];
        for(j=0; j<Mt->give_tnip(i); j++)
        {
          mt_ip[ippm] = Tt->elements[i].ipp[0][0]+j;
          ippm++;
        }
        break;
      default:
        print_err("unknown element type is required", __FILE__, __LINE__, __func__);
    }
  }
}



/**
  The function transfers MEFEL quantities to TRFEL with respect to required nontransport 
  quantities in TRFEL. Meshes must be identical in TRFEL and MEFEL and
  individual quantities are copied between corresponding integration points.

  @return The function does not return anything but it changes values in the Tm->nontransq array.

  Created by Tomas Koudelka 14.11.2013
*/
void mefel_trfel_copyip(void)
{
  long i,j;
  double ipval;
  long antq[Tm->tnkntq];

  // transfer of MEFEL quantities to TRFEL
  for(i=0; i<Tm->tnip; i++)
  {
    memset(antq, 0, sizeof(*antq)*Tm->tnkntq);
    Tm->give_reqntq(i, antq); // search for required non-transport quantities on the i-th integration point

    for (j=0; j<Tm->nntq; j++)
    {
      if (antq[int(Tm->ntqo[j])-1] == 1) // the model requires j-th quantity of totally used quantities
      {  
        ipval = Mm->givemechq(Tm->ntqo[j], TM_ip[i]);
        Tm->storenontransq(Tm->ntqo[j], i, ipval);
      }
      else // the quantity is not required => store 0.0
        Tm->storenontransq(Tm->ntqo[j], i, 0.0);
    }
  }
}



/**
  The function transfers TRFEL quantities to MEFEL with respect to required nonmechanical
  quantities in MEFEL. Meshes must be identical in MEFEL and TRFEL and
  individual quantities are copied between corresponding integration points.

  @return The function does not return anything but it changes values in the Mm->nonmechq array.

  Created by Tomas Koudelka 14.11.2013
*/
void trfel_mefel_copyip(void)
{
  long i,j;
  double ipval;
  long anmq[Mm->tnknmq];

  // transfer of TRFEL quantities to MEFEL
  for(i=0; i<Mm->tnip; i++)
  {
    memset(anmq, 0, sizeof(*anmq)*Mm->tnknmq);
    for (j=0; j<Mm->ip[i].nm; j++)
      Mm->give_reqnmq(i, j, anmq); // search for required non-transport quantities on the i-th integration point

    for (j=0; j<Mm->nnmq; j++)
    {
      if (Mm->nmqo[j] == initial_temperature) // initial values may be obtained only for the initial time by call init_trfel_mefel_copyip
        continue;

      if (anmq[int(Mm->nmqo[j])-1] == 1) // the model requires j-th quantity of totally used quantities
      {  
        if (MT_ip[i] < 0) // the fifth integration point on plelemlq
        {
          ipval  = Tm->givetransq(Mm->nmqo[j], MT_ip[i-4]);
          ipval += Tm->givetransq(Mm->nmqo[j], MT_ip[i-3]);
          ipval += Tm->givetransq(Mm->nmqo[j], MT_ip[i-2]);
          ipval += Tm->givetransq(Mm->nmqo[j], MT_ip[i-1]);
          ipval /= 4.0;
        }
        else
          ipval = Tm->givetransq(Mm->nmqo[j], MT_ip[i]);
        Mm->storenonmechq(Mm->nmqo[j], i, ipval);
      }
      else // the quantity is not required on the i-th integartion point => store 0.0
        Mm->storenonmechq(Mm->nmqo[j], i, 0.0);
    }
  }
}



/**
  The function transfers TRFEL intial values of quantities to MEFEL with respect to required nonmechanical
  quantities in MEFEL. Meshes must be identical in MEFEL and TRFEL and
  individual quantities are copied between corresponding integration points.

  @return The function does not return anything but it changes values in the Mm->nonmechq array.

  Created by Tomas Koudelka 14.11.2013
*/
void init_trfel_mefel_copyip(void)
{
  long i,j;
  double ipval;
  long anmq[Mm->tnknmq];

  // transfer of TRFEL quantities to MEFEL
  for(i=0; i<Mm->tnip; i++)
  {
    memset(anmq, 0, sizeof(*anmq)*Mm->tnknmq);
    for (j=0; j<Mm->ip[i].nm; j++)
      Mm->give_reqnmq(i, j, anmq); // search for required non-transport quantities on the i-th integration point

    for (j=0; j<Mm->nnmq; j++)
    {
      if (anmq[int(Mm->nmqo[j])-1] == 1) // the model requires j-th quantity of totally used quantities
      {  
        if (MT_ip[i] < 0) // the fifth integration point on plelemlq
        {
          ipval  = Tm->givetransq(Mm->nmqo[j], MT_ip[i-4]);
          ipval += Tm->givetransq(Mm->nmqo[j], MT_ip[i-3]);
          ipval += Tm->givetransq(Mm->nmqo[j], MT_ip[i-2]);
          ipval += Tm->givetransq(Mm->nmqo[j], MT_ip[i-1]);
          ipval /= 4.0;
        }
        else
          ipval = Tm->givetransq(Mm->nmqo[j], MT_ip[i]);
        Mm->storenonmechq(Mm->nmqo[j], i, ipval);
      }
      else // the quantity is not required on the i-th integartion point => store 0.0
        Mm->storenonmechq(Mm->nmqo[j], i, 0.0);
    }
  }
}

/**
  The function extracts displacements on one element in fully coupled problems
   
  @param lcid - number of load case
  @param eid - element id
  @param r - allocated array for displacement, it is output parameter
   
  @return The function returns nodal displacements at the given element in the array r.

  JK, 15.4.2019
*/
void eldispl_fc (long lcid,long eid,double *r)
{
  long i,ii,ndofe,ndofem;
  ivector cn,nod,ncn;
  vector rr;
  
  /*
  //  the number of DOFs on element
  //  this number is equal to the number of DOFs without hanging nodes
  ndofe = Mt->give_ndofe (eid);
  //  the number of DOFs on element
  //  this number is equal to the number of DOFs with master nodes generated by hanging nodes
  ndofemn = Gtm->give_ndofe (eid);
  
  if (ndofe == ndofemn){
    //  there are no hanging nodes on element
    reallocv (RSTCKIVEC(ndofe,cn));
    reallocv (RSTCKVEC(ndofe,rr));
  }
  else{
    //  there are hanging nodes on element
    reallocv (RSTCKIVEC(ndofemn,cn));
    reallocv (RSTCKVEC(ndofemn,rr));
    ndofe=ndofemn;
  }
  */

  ndofem = ndofe = Gtm->give_ndofe (eid);
  reallocv (RSTCKIVEC(ndofe,cn));
  reallocv (RSTCKVEC(ndofe,rr));
  
  //  code numbers on element
  Gtu->give_code_numbers (eid,cn.a);
  
  switch (Cp->tprob){
  case fully_coupled_material:{
    //  lcid must be equal to zero
    for (i=0;i<ndofe;i++){
      ii=cn[i];
      if (ii<0)   rr[i]=Cb->dlc[lcid].get_pd (Cp->time, ii);
      if (ii==0)  rr[i]=0.0;
      if (ii>0)   rr[i]=Lsrsc->lhs[lcid*Ndofc+ii-1];
    }
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }//  end of the switch (Cp->tprob)
  
  /*
  //  the variables ndofe and ndofemn have to be obtained again
  //  because they are possibly rewritten
  //
  //  the number of DOFs on element
  //  this number is equal to the number of DOFs without hanging nodes
  ndofe = Mt->give_ndofe (eid);
  //  the number of DOFs on element
  //  this number is equal to the number of DOFs with master nodes generated by hanging nodes
  ndofemn = Gtm->give_ndofe (eid);
  
  if (ndofe != ndofemn){
    if ((Mp->homog!=3) || (Mp->homog!=5) || (Mp->homog!=7)){
      //  this case means hanging nodes
      
      //  there are hanging nodes on element
      mxv (Mt->elements[eid].tmat->a,rr.a,r,ndofe,ndofemn);
    }
    else{
      //  this case means homogenization
      //  the vector rr contains more components than the vector r
      //  only the "classical" components are copied
      copyv (rr.a,r,ndofe);
    }
  }
  else{
    //  there are no hanging nodes on element
    copyv (rr,r);
  }
  
  if (Mp->tprob == growing_mech_structure){
    //  subtraction of initial displacements
    Mt->elements[eid].subtrinitdispl (r, ndofe);
  }
  */
  
  copyv (rr,r);
  
  destrv (rr);
  destrv (cn);
}


