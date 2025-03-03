#include "quadrilatc.h"
#include "quadlineart.h"
#include "plelemlq.h"
#include "plelemqq.h"
#include "global.h"
#include "globalt.h"
#include "globalc.h"
#include "element.h"
#include "intpoints.h"
#include "globmatt.h"


quadrilatc::quadrilatc (void)
{
  long i;
  
  if (Cp->lbb==lin_lin){
    //  number of blocks of the mechanical element
    mnb=Pelq->nb;
    //  number of mechanical DOFs
    mndofe=Pelq->ndofe;
    //  number of strain/stress components
    tnmcomp=Pelq->tncomp;
    //  number of nodes in mechanical problems
    nnemp = Pelq->nne;
    //  number of nodes in transport problems
    nnetp = Lqt->nne;
    //  number of transport DOFs
    tndofe = Lqt->ndofe;
  }

  if (Cp->lbb==quad_lin){
    //  number of blocks of the mechanical element
    mnb=Peqq->nb;
    //  number of mechanical DOFs
    mndofe=Peqq->ndofe;
    //  number of strain/stress components
    tnmcomp=Peqq->tncomp;
    //  number of nodes in mechanical problems
    nnemp = Peqq->nne;
    //  number of nodes in transport problems
    nnetp = Lqt->nne;
    //  number of transport DOFs
    tndofe = Lqt->ndofe;
  }

  if (Cp->lbb==quad_quad){
    //  number of blocks of the mechanical element
    mnb=Peqq->nb;
    //  number of mechanical DOFs
    mndofe=Peqq->ndofe;
    //  number of strain/stress components
    tnmcomp=Peqq->tncomp;
    //  number of nodes in mechanical problems
    nnemp = Peqq->nne;
    //  number of nodes in transport problems
    nnetp = Qqt->nne;
    //  number of transport DOFs
    tndofe = Qqt->ndofe;
  }
  
  
  //  number of transported media
  ntm=Tp->ntm;
  
  
  
  intordvum=NULL;  intordvlm=NULL;
  nipu=NULL;  nipl=NULL;
  
  //  number of DOFs in transport part
  dofe = new long [ntm];
  //  number of integration points for upper part of matrices (mnb,ntm)
  nipu = new long* [mnb];
  //  number of integration points for lower part of matrices (ntm,mnb)
  nipl = new long* [ntm];
  //  order of numerical integration for upper part of matrices (mnb,ntm)
  intordvum = new long* [mnb];
  //  order of numerical integration for lower part of matrices (ntm,mnb)
  intordvlm = new long* [ntm];

  for (i=0;i<mnb;i++){
    nipu[i] = new long [ntm];
    intordvum[i] = new long [ntm];
  }
  for (i=0;i<ntm;i++){
    nipl[i] = new long [mnb];
    intordvlm[i] = new long [mnb];
  }
  
  //  ordering of unknowns from mechanical part
  mordering = new long [mndofe];
  
  for (i=0;i<mndofe;i++){
    mordering[i]=i+1;
  }
  
  
  if (Cp->lbb==lin_lin){
    switch (Tp->tmatt){
    case onemedium:{
      intordvum[0][0]=2;  intordvlm[0][0]=2;
      nipu[0][0]=4;  nipl[0][0]=4;
      dofe[0]=nnetp;  
      break;
    }
    case twomediacoup:{
      intordvum[0][0]=2;  intordvum[0][1]=2;
      intordvlm[0][0]=2;  intordvlm[1][0]=2;

      if (Cp->savemode==0){
	nipu[0][0]=4;       nipu[0][1]=4;
	nipl[0][0]=4;       nipl[1][0]=4;
      }
      if (Cp->savemode==1){
	nipu[0][0]=4;       nipu[0][1]=0;
	nipl[0][0]=4;       nipl[1][0]=0;
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;
      break;
    }
    case threemediacoup:{
      intordvum[0][0]=2;  intordvum[0][1]=2;  intordvum[0][2]=2;
      intordvlm[0][0]=2;  intordvlm[1][0]=2;  intordvlm[2][0]=2;
      
      if (Cp->savemode==0){
	nipu[0][0]=4;  nipu[0][1]=4;  nipu[0][2]=4;
	nipl[0][0]=4;  nipl[1][0]=4;  nipl[2][0]=4;
      }
      if (Cp->savemode==1){
	nipu[0][0]=4;  nipu[0][1]=0;  nipu[0][2]=0;
	nipl[0][0]=4;  nipl[1][0]=0;  nipl[2][0]=0;
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;  dofe[2]=nnetp;
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown type of transported matter is required");
      fprintf (stderr,"\n in function quadrilatc::quadrilatc (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }

  if (Cp->lbb==quad_lin){
    switch (Tp->tmatt){
    case onemedium:{
      intordvum[0][0]=2;  intordvlm[0][0]=2;
      nipu[0][0]=4;  nipl[0][0]=4;
      dofe[0]=nnetp;
      break;
    }
    case twomediacoup:{
      intordvum[0][0]=2;  intordvum[0][1]=2;
      intordvlm[0][0]=2;  intordvlm[1][0]=2;
      
      if (Cp->savemode==0){
	nipu[0][0]=4;       nipu[0][1]=4;
	nipl[0][0]=4;       nipl[1][0]=4;
      }
      if (Cp->savemode==1){
	nipu[0][0]=4;       nipu[0][1]=0;
	nipl[0][0]=4;       nipl[1][0]=0;
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;
      break;
    }
    case threemediacoup:{
      intordvum[0][0]=2;  intordvum[0][1]=2;  intordvum[0][2]=2;
      intordvlm[0][0]=2;  intordvlm[1][0]=2;  intordvlm[2][0]=2;
      
      if (Cp->savemode==0){
	nipu[0][0]=4;  nipu[0][1]=4;  nipu[0][2]=4;
	nipl[0][0]=4;  nipl[1][0]=4;  nipl[2][0]=4;
      }
      if (Cp->savemode==1){
	nipu[0][0]=4;  nipu[0][1]=0;  nipu[0][2]=0;
	nipl[0][0]=4;  nipl[1][0]=0;  nipl[2][0]=0;
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;  dofe[2]=nnetp;
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown type of transported matter is required");
      fprintf (stderr,"\n in function quadrilatc::quadrilatc (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
  
  if (Cp->lbb==quad_quad){
    switch (Tp->tmatt){
    case onemedium:{
      intordvum[0][0]=3;  intordvlm[0][0]=3;
      nipu[0][0]=9;  nipl[0][0]=9;
      dofe[0]=nnetp;
      break;
    }
    case twomediacoup:{
      intordvum[0][0]=3;  intordvum[0][1]=3;
      intordvlm[0][0]=3;  intordvlm[1][0]=3;
      
      if (Cp->savemode==0){
	nipu[0][0]=9;       nipu[0][1]=9;
	nipl[0][0]=9;       nipl[1][0]=9;
      }
      if (Cp->savemode==1){
	nipu[0][0]=9;       nipu[0][1]=0;
	nipl[0][0]=9;       nipl[1][0]=0;
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;
      break;
    }
    case threemediacoup:{
      intordvum[0][0]=3;  intordvum[0][1]=3;  intordvum[0][2]=3;
      intordvlm[0][0]=3;  intordvlm[1][0]=3;  intordvlm[2][0]=3;
      
      if (Cp->savemode==0){
	nipu[0][0]=9;  nipu[0][1]=9;  nipu[0][2]=9;
	nipl[0][0]=9;  nipl[1][0]=9;  nipl[2][0]=9;
      }
      if (Cp->savemode==1){
	nipu[0][0]=9;  nipu[0][1]=0;  nipu[0][2]=0;
	nipl[0][0]=9;  nipl[1][0]=0;  nipl[2][0]=0;
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;  dofe[2]=nnetp;
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown type of transported matter is required");
      fprintf (stderr,"\n in function quadrilatc::quadrilatc (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }


}

quadrilatc::~quadrilatc (void)
{
  long i;

  for (i=0;i<ntm;i++){
    delete [] nipl[i];
    delete [] intordvlm[i];
  }
  for (i=0;i<mnb;i++){
    delete [] nipu[i];
    delete [] intordvum[i];
  }
  delete [] dofe;
  delete [] nipu;
  delete [] nipl;
  delete [] intordvum;
  delete [] intordvlm;
}

void quadrilatc::eleminit (long eid)
{
  long ii,jj;

  //Ct->elements[eid].nb=nb;
  Ct->elements[eid].intordvum = new long* [mnb];
  Ct->elements[eid].intordvlm = new long* [ntm];
  Ct->elements[eid].nipu = new long* [mnb];
  Ct->elements[eid].nipl = new long* [ntm];
  
  for (ii=0;ii<mnb;ii++){
    Ct->elements[eid].intordvum[ii] = new long [ntm];
    Ct->elements[eid].nipu[ii] = new long [ntm];
    for (jj=0;jj<ntm;jj++){
      Ct->elements[eid].intordvum[ii][jj]=intordvum[ii][jj];
      Ct->elements[eid].nipu[ii][jj]=nipu[ii][jj];
    }
  }
  for (ii=0;ii<ntm;ii++){
    Ct->elements[eid].intordvlm[ii] = new long [mnb];
    Ct->elements[eid].nipl[ii] = new long [mnb];
    for (jj=0;jj<mnb;jj++){
      Ct->elements[eid].intordvlm[ii][jj]=intordvlm[ii][jj];
      Ct->elements[eid].nipl[ii][jj]=nipl[ii][jj];
    }
  }
}


/**
   function computes upper coupling conductivity %matrix of 2D problems
   
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting matrix
   @param vm - coupling %matrix
   
   JK, 9.1.2003
*/
void quadrilatc::upper_cond_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,j,ii;
  double thick,xi,eta,jac;
  ivector nodes(nnetp);
  vector x(nnemp),y(nnemp),w(intordvum[ri][ci]),gp(intordvum[ri][ci]),t(nnetp);
  matrix gm(tnmcomp,mndofe),d(tnmcomp,1),n(1,nnetp);

  Tt->give_elemnodes (eid,nodes);
  Tc->give_thickness (eid,nodes,t);

  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordvum[ri][ci]);

  fillm (0.0,vm);
  
  if (Cp->savemode==0){
    ii=Ct->elements[eid].ippu[ri][ci];
  }
  if (Cp->savemode==1){
    ii=Ct->elements[eid].ippu[0][0];
  }
  
  for (i=0;i<intordvum[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordvum[ri][ci];j++){
      eta=gp[j];
      
      if (Cp->lbb==lin_lin){
	Pelq->geom_matrix (gm,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->lbb==quad_lin){
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->lbb==quad_quad){
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Qqt->bf_matrix (n,xi,eta);
      }
      
      //  matrix of conductivity of the material
      Cmu->matcond (d,ii,ri,ci);
      
      if (Cp->lbb==lin_lin || Cp->lbb==quad_lin)
	thick = Lqt->approx (xi,eta,t);
      else
	thick = Qqt->approx (xi,eta,t);
      
      jac*=thick*w[i]*w[j];
      
      //  contribution to the conductivity matrix of the element
      bdbjac (vm,gm,d,n,jac);
      
      ii++;
    }
  }
  
}


/**
   function computes lower coupling conductivity %matrix of 2D problems
   
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param vm - coupling %matrix
   
   JK, 9.1.2003
*/
void quadrilatc::lower_cond_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,j,ii;
  double thick,xi,eta,ww1,ww2,jac;
  ivector nodes(nnetp);
  vector x(nnemp),y(nnemp),w(intordvlm[ri][ci]),gp(intordvlm[ri][ci]),t(nnetp);
  matrix gm(tnmcomp,mndofe),d(1,tnmcomp),n(1,nnetp);

  Tt->give_elemnodes (eid,nodes);
  Tc->give_thickness (eid,nodes,t);

  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordvum[ri][ci]);

  fillm (0.0,vm);
  
  if (Cp->savemode==0){
    ii=Ct->elements[eid].ippl[ri][ci];
  }
  if (Cp->savemode==1){
    ii=Ct->elements[eid].ippl[0][0];
  }

  for (i=0;i<intordvlm[ri][ci];i++){
    xi=gp[i];  ww1=w[i];
    for (j=0;j<intordvlm[ri][ci];j++){
      eta=gp[j];  ww2=w[j];
      
      if (Cp->lbb==lin_lin){
	Pelq->geom_matrix (gm,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->lbb==quad_lin){
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->lbb==quad_quad){
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Qqt->bf_matrix (n,xi,eta);
      }

      //  matrix of conductivity of the material
      Cml->matcond (d,ii,ri,ci);

      if (Cp->lbb==lin_lin || Cp->lbb==quad_lin)
	thick = Lqt->approx (xi,eta,t);
      else
	thick = Qqt->approx (xi,eta,t);
      
      jac*=thick*ww1*ww2;
      
      //  contribution to the conductivity matrix of the element
      bdbjac (vm,n,d,gm,jac);

      ii++;
    }
  }
  
}

/**
   function computes upper coupling capacity %matrix of 2D problems
   
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param vm - coupling %matrix
   
   JK, 9.1.2003
*/
void quadrilatc::upper_cap_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,j,ii;
  double thick,xi,eta,ww1,ww2,jac;
  ivector nodes(nnetp);
  vector x(nnemp),y(nnemp),w(intordvum[ri][ci]),gp(intordvum[ri][ci]),t(nnetp);
  matrix gm(tnmcomp,mndofe),d(tnmcomp,1),n(1,nnetp);

  Tt->give_elemnodes (eid,nodes);
  Tc->give_thickness (eid,nodes,t);

  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordvum[ri][ci]);

  fillm (0.0,vm);

  if (Cp->savemode==0){
    ii=Ct->elements[eid].ippu[ri][ci];
  }
  if (Cp->savemode==1){
    ii=Ct->elements[eid].ippu[0][0];
  }

  for (i=0;i<intordvum[ri][ci];i++){
    xi=gp[i];  ww1=w[i];
    for (j=0;j<intordvum[ri][ci];j++){
      eta=gp[j];  ww2=w[j];
      
      if (Cp->lbb==lin_lin){
	Pelq->geom_matrix (gm,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->lbb==quad_lin){
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->lbb==quad_quad){
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Qqt->bf_matrix (n,xi,eta);
      }

      
      //  matrix of conductivity of the material
      Cmu->matcap (d,ii,ri,ci);

      if (Cp->lbb==lin_lin || Cp->lbb==quad_lin)
	thick = Lqt->approx (xi,eta,t);
      else
	thick = Qqt->approx (xi,eta,t);

      jac*=thick*ww1*ww2;

      //  contribution to the conductivity matrix of the element
      bdbjac (vm,gm,d,n,jac);

      ii++;
    }
  }

}


/**
   function computes lower coupling capacity %matrix of 2D problems
   
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param vm - coupling %matrix
   
   JK, 9.1.2003
*/
void quadrilatc::lower_cap_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,j,ii;
  double thick,xi,eta,ww1,ww2,jac;
  ivector nodes(nnetp);
  vector x(nnemp),y(nnemp),w(intordvlm[ri][ci]),gp(intordvlm[ri][ci]),t(nnetp);
  matrix gm(tnmcomp,mndofe),d(1,tnmcomp),n(1,nnetp);

  Tt->give_elemnodes (eid,nodes);
  Tc->give_thickness (eid,nodes,t);

  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordvum[ri][ci]);

  fillm (0.0,vm);
  
  if (Cp->savemode==0){
    ii=Ct->elements[eid].ippl[ri][ci];
  }
  if (Cp->savemode==1){
    ii=Ct->elements[eid].ippl[0][0];
  }
  
  for (i=0;i<intordvlm[ri][ci];i++){
    xi=gp[i];  ww1=w[i];
    for (j=0;j<intordvlm[ri][ci];j++){
      eta=gp[j];  ww2=w[j];
      
      if (Cp->lbb==lin_lin){
	Pelq->geom_matrix (gm,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->lbb==quad_lin){
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->lbb==quad_quad){
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Qqt->bf_matrix (n,xi,eta);
      }
      
      //  matrix of conductivity of the material
      Cml->matcap (d,ii,ri,ci);
      
      if (Cp->lbb==lin_lin || Cp->lbb==quad_lin)
	thick = Lqt->approx (xi,eta,t);
      else
	thick = Qqt->approx (xi,eta,t);
      
      jac*=thick*ww1*ww2;
      
      //  contribution to the conductivity matrix of the element
      bdbjac (vm,n,d,gm,jac);

      ii++;
    }
  }

}




/**
   function assembles upper coupling conductivity %matrices into one element %matrix
   
   @param eid - element id
   @param vm - element upper coupling conductivity %matrix
   
   JK, 17.7.2005
 */
void quadrilatc::res_upper_cond_coup_matrix (long eid,matrix &vm)
{
  long i,j,*ccn;
  matrix lvm(mndofe,nnetp);
  
  ccn = new long [nnetp];

  fillm(0.0,vm);

  for (i=0;i<mnb;i++){
    for (j=0;j<ntm;j++){
      
      //  computation of submatrices
      upper_cond_coup_matrix (eid,i,j,lvm);
	
      if (Cp->lbb==lin_lin || Cp->lbb==quad_lin)
	Lqt->codnum (ccn,j);
      if (Cp->lbb==quad_quad)
	Qqt->codnum (ccn,j);
      
      mat_localize (vm,lvm,mordering,ccn);
    }
  }
  delete [] ccn;
}

/**
   function assembles lower coupling conductivity matrices into one element %matrix
   
   @param eid - element id
   @param vm - element lower coupling conductivity %matrix
   
   JK, 17.7.2005
 */
void quadrilatc::res_lower_cond_coup_matrix (long eid,matrix &vm)
{
  long i,j,*rcn;
  matrix lvm(nnetp,mndofe);
  
  rcn = new long [nnetp];

  fillm(0.0,vm);
  
  for (i=0;i<ntm;i++){
    for (j=0;j<mnb;j++){
      
      //  computation of submatrices
      lower_cond_coup_matrix (eid,i,j,lvm);
      
      if (Cp->lbb==lin_lin || Cp->lbb==quad_lin)
	Lqt->codnum (rcn,j);
      else
	Qqt->codnum (rcn,j);
      
      mat_localize (vm,lvm,rcn,mordering);
    }
  }
  delete [] rcn;
}

/**
   function assembles upper coupling capacity matrices into one element %matrix
   
   @param eid -element id
   @param vm - element upper coupling capacity %matrix
   
   JK, 17.7.2005
*/
void quadrilatc::res_upper_cap_coup_matrix (long eid,matrix &vm)
{
  long i,j,*ccn;
  matrix lvm(mndofe,nnetp);
  
  ccn = new long [nnetp];
  
  fillm(0.0,vm);
  
  for (i=0;i<mnb;i++){
    for (j=0;j<ntm;j++){
      
      //  computation of submatrices
      upper_cap_coup_matrix (eid,i,j,lvm);
      
      if (Cp->lbb==lin_lin || Cp->lbb==quad_lin)
	Lqt->codnum (ccn,j);
      else
	Qqt->codnum (ccn,j);
      
      mat_localize (vm,lvm,mordering,ccn);
    }
  }
  delete [] ccn;
}

/**
   function assembles lower coupling capacity matrices into one element %matrix
   
   @param eid -element id
   @param vm - element lower coupling capacity %matrix
   
   JK, 17.7.2005
*/
void quadrilatc::res_lower_cap_coup_matrix (long eid,matrix &vm)
{
  long i,j,*rcn;
  matrix lvm(nnetp,mndofe);
  
  rcn = new long [nnetp];

  fillm(0.0,vm);

  for (i=0;i<ntm;i++){
    for (j=0;j<mnb;j++){

      //  computation of submatrices
      lower_cap_coup_matrix (eid,i,j,lvm);
      
      if (Cp->lbb==lin_lin || Cp->lbb==quad_lin)
	Lqt->codnum (rcn,j);
      else
	Qqt->codnum (rcn,j);
      
      mat_localize (vm,lvm,rcn,mordering);
    }
  }
  delete [] rcn;
}

















/**
   function copies mechanical data from MEFEL into integration points
   of coupled part in METR
   
   @param eid - element id
   
   JK, 29.10.2004
*/
void quadrilatc::mefel_metr (long eid)
{
  long i,j,ippm,ippcu,ippcl;
  
  //  save mode - one set of integration points is used for all coupling blocks
  //  this strategy saves significantly memory but it is not the most general approach
  
  if (Cp->savemode==1){
    //  number of the first integration point on mechanical element
    //  one block formulation of mechanical problem is assumed
    ippm = Mt->elements[eid].ipp[0][0];
    
    //  number of the fisrt integration point on coupling element
    //  for upper part
    ippcu = Ct->elements[eid].ippu[0][0];
    
    //  number of the fisrt integration point on coupling element
    //  for upper part
    ippcl = Ct->elements[eid].ippl[0][0];
    
    for (j=0;j<nipu[0][0];j++){
      for (i=0;i<tnmcomp;i++){
	Cmu->ip[ippcu].strains[i]  = Mm->ip[ippm].strain[i];
	Cmu->ip[ippcu].stresses[i] = Mm->ip[ippm].stress[i];
      }
      ippm++;  ippcu++;
    }
    
    //  number of the first integration point on mechanical element
    //  one block formulation of mechanical problem is assumed
    ippm = Mt->elements[eid].ipp[0][0];
    for (j=0;j<nipl[0][0];j++){
      for (i=0;i<tnmcomp;i++){
	Cml->ip[ippcl].strains[i]  = Mm->ip[ippm].strain[i];
	Cml->ip[ippcl].stresses[i] = Mm->ip[ippm].stress[i];
      }
      ippm++;  ippcl++;
    }
  }
  
  
  //  each coupling block has its own set of integration points
  //  it is the most general approach but not implemented at this moment
  //  but the code is ready for this extension
  
  if (Cp->savemode==0){
    fprintf (stderr,"\n\n not implemented at this moment (file %s, line %d).\n",__FILE__,__LINE__);
  }
  
}

/**
   function copies transport data from TRFEL into integration points
   of coupled part in METR
   
   @param eid - element id
   
   JK, 29.10.2004
*/
void quadrilatc::trfel_metr (long eid)
{
  long i,j,ippt,ippcu,ippcl;
  
  //  save mode - one set of integration points is used for all coupling blocks
  //  this strategy saves significantly memory but it is not the most general approach
  
  if (Cp->savemode==1){
    //  number of the first integration point on transport element
    ippt = Tt->elements[eid].ipp[0][0];
    
    //  number of the fisrt integration point on coupling element
    //  for upper part
    ippcu = Ct->elements[eid].ippu[0][0];
    
    //  number of the fisrt integration point on coupling element
    //  for upper part
    ippcl = Ct->elements[eid].ippl[0][0];
    
    for (j=0;j<nipu[0][0];j++){
      for (i=0;i<ntm;i++){
	Cmu->ip[ippcu].av[i] = Tm->ip[ippt].av[i];
	//Cmu->ip[ippcu].pv[i] = Tm->ip[ippt].pv[i];
      }
      ippt++;  ippcu++;
    }
    
    //  number of the first integration point on transport element
    ippt = Mt->elements[eid].ipp[0][0];
    for (j=0;j<nipl[0][0];j++){
      for (i=0;i<ntm;i++){
	Cml->ip[ippcl].av[i] = Tm->ip[ippt].av[i];
	//Cml->ip[ippcl].pv[i] = Tm->ip[ippt].pv[i];
      }
      ippt++;  ippcl++;
    }
  }
  
  
  //  each coupling block has its own set of integration points
  //  it is the most general approach but not implemented at this moment
  //  but the code is ready for this extension
  
  if (Cp->savemode==0){
    fprintf (stderr,"\n\n not implemented at this moment (file %s, line %d).\n",__FILE__,__LINE__);
  }
  
}

/**
   function approximates temperature into mechanical integration points
   
   @param eid - element id
   
   JK, 31.10.2004
*/
void quadrilatc::trfel_mefel (long eid)
{
  long i,j,k,ipp,into;
  double xi,eta,val;
  ivector cn(tndofe);
  vector r(tndofe),t(nnetp),gp,w;
  
  //  code numbers of transport part
  Tt->give_code_numbers (eid,cn.a);
  //  nodal values of transport part
  nodalvalues (0,r.a,cn.a,tndofe);
  
  //  index of temperature
  k=Tp->ntm-1;
  
  //  selection of nodal temperatures from all nodal values
  for (i=0;i<dofe[k];i++){
    t[i]=r[Lqt->ordering[k][i]-1];
  }
  
  //  order of numerical integration of mechanical part
  into = Peqq->intordsm[0][0];
  
  allocv (into,gp);
  allocv (into,w);
  gauss_points (gp.a,w.a,into);
  
  //  number of the first integration point on element
  ipp=Mt->elements[eid].ipp[0][0];
  
  for (i=0;i<into;i++){
    xi=gp[i];
    for (j=0;j<into;j++){
      eta=gp[j];
      
      val = Lqt->approx (xi,eta,t);
      Mm->tempr[ipp]=val;
      ipp++;
    }
  }
  
  destrv (gp);  destrv (w);
  
}
