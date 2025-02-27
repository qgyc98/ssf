#include "hexahedc.h"
#include "global.h"
#include "mechtop.h"
#include "mechmat.h"
#include "globalt.h"
#include "globalc.h"
#include "linhex.h"
#include "quadhex.h"
#include "linhext.h"
#include "quadhext.h"
#include "genfile.h"
#include "globmatt.h"
#include "globmatc.h"
#include "globmat.h"
#include "node.h"
#include "element.h"
#include "intpoints.h"

hexahedc::hexahedc (void)
{
  long i,j;
  
  if (Cp->bb==lin_lin){
    //  number of blocks of the mechanical element
    mnb=Lhex->nb;
    //  number of mechanical DOFs
    mndofe=Lhex->ndofe;
    //  number of nodes in mechanical problems
    nnemp = Lhex->nne;
    //  number of strain/stress components
    tnmcomp=Lhex->tncomp;
    //  array of numbers of components in strain vectors
    mncomp=Lhex->ncomp;
    ssst = Lhex->ssst;
    //  number of nodes in transport problems
    nnetp = Lht->nne;
    //  number of transport DOFs
    tndofe = Lht->ndofe;
  }

  if (Cp->bb==quad_lin){
    //  number of blocks of the mechanical element
    mnb=Qhex->nb;
    //  number of mechanical DOFs
    mndofe=Qhex->ndofe;
    //  number of nodes in mechanical problems
    nnemp = Qhex->nne;
    //  number of strain/stress components
    tnmcomp=Qhex->tncomp;
    //  array of numbers of components in strain vectors
    mncomp=Qhex->ncomp;
    ssst = Qhex->ssst;
    //  number of nodes in transport problems
    nnetp = Lht->nne;
    //  number of transport DOFs
    tndofe = Lht->ndofe;
  }

  if (Cp->bb==quad_quad){
    //  number of blocks of the mechanical element
    mnb=Qhex->nb;
    //  number of mechanical DOFs
    mndofe=Qhex->ndofe;
    //  number of nodes in mechanical problems
    nnemp = Qhex->nne;
    //  number of strain/stress components
    tnmcomp=Qhex->tncomp;
    //  array of numbers of components in strain vectors
    mncomp=Qhex->ncomp;
    ssst = Qhex->ssst;
    //  number of nodes in transport problems
    nnetp = Qht->nne;
    //  number of transport DOFs
    tndofe = Qht->ndofe;
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
    
  //  ordering of unknowns from transport part
  tordering = new long* [ntm];
  for (i=0;i<ntm;i++){
    tordering[i] = new long [nnetp];
  }

  
  if (Cp->bb==lin_lin){
    switch (Tp->tmatt){
    case onemedium:{
      intordvum[0][0]=2;  intordvlm[0][0]=2;
      nipu[0][0]=8;  nipl[0][0]=8;
      dofe[0]=nnetp;
      break;
    }
    case twomediacoup:{
      intordvum[0][0]=2;  intordvum[0][1]=2;
      intordvlm[0][0]=2;  intordvlm[1][0]=2;

      if (Cp->savemode==0){
	nipu[0][0]=8;       nipu[0][1]=8;
	nipl[0][0]=8;       nipl[1][0]=8;
      }
      if (Cp->savemode==1){
	nipu[0][0]=8;       nipu[0][1]=0;
	nipl[0][0]=8;       nipl[1][0]=0;
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;
      break;
    }
    case threemediacoup:{
      intordvum[0][0]=2;  intordvum[0][1]=2;  intordvum[0][2]=2;
      intordvlm[0][0]=2;  intordvlm[1][0]=2;  intordvlm[2][0]=2;
      
      if (Cp->savemode==0){
	nipu[0][0]=8;  nipu[0][1]=8;  nipu[0][2]=8;
	nipl[0][0]=8;  nipl[1][0]=8;  nipl[2][0]=8;
      }
      if (Cp->savemode==1){
	nipu[0][0]=8;  nipu[0][1]=0;  nipu[0][2]=0;
	nipl[0][0]=8;  nipl[1][0]=0;  nipl[2][0]=0;
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;  dofe[2]=nnetp;
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown type of transported matter is required");
      fprintf (stderr,"\n in function hexahedc::hexahedc (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }

  if (Cp->bb==quad_lin){
    switch (Tp->tmatt){
    case onemedium:{
      intordvum[0][0]=2;  intordvlm[0][0]=2;
      nipu[0][0]=8;  nipl[0][0]=8;
      dofe[0]=nnetp;
      break;
    }
    case twomediacoup:{
      intordvum[0][0]=2;  intordvum[0][1]=2;
      intordvlm[0][0]=2;  intordvlm[1][0]=2;
      
      if (Cp->savemode==0){
	nipu[0][0]=8;       nipu[0][1]=8;
	nipl[0][0]=8;       nipl[1][0]=8;
      }
      if (Cp->savemode==1){
	nipu[0][0]=8;       nipu[0][1]=0;
	nipl[0][0]=8;       nipl[1][0]=0;
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;
      break;
    }
    case threemediacoup:{
      intordvum[0][0]=2;  intordvum[0][1]=2;  intordvum[0][2]=2;
      intordvlm[0][0]=2;  intordvlm[1][0]=2;  intordvlm[2][0]=2;
      
      if (Cp->savemode==0){
	nipu[0][0]=8;  nipu[0][1]=8;  nipu[0][2]=8;
	nipl[0][0]=8;  nipl[1][0]=8;  nipl[2][0]=8;
      }
      if (Cp->savemode==1){
	nipu[0][0]=8;  nipu[0][1]=0;  nipu[0][2]=0;
	nipl[0][0]=8;  nipl[1][0]=0;  nipl[2][0]=0;
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;  dofe[2]=nnetp;
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown type of transported matter is required");
      fprintf (stderr,"\n in function hexahedc::hexahedc (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
  
  if (Cp->bb==quad_quad){
    switch (Tp->tmatt){
    case onemedium:{
      intordvum[0][0]=3;  intordvlm[0][0]=3;
      nipu[0][0]=27;  nipl[0][0]=27;
      dofe[0]=nnetp;
      break;
    }
    case twomediacoup:{
      intordvum[0][0]=3;  intordvum[0][1]=3;
      intordvlm[0][0]=3;  intordvlm[1][0]=3;
      
      if (Cp->savemode==0){
	nipu[0][0]=27;       nipu[0][1]=27;
	nipl[0][0]=27;       nipl[1][0]=27;
      }
      if (Cp->savemode==1){
	nipu[0][0]=27;       nipu[0][1]=0;
	nipl[0][0]=27;       nipl[1][0]=0;
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;
      break;
    }
    case threemediacoup:{
      intordvum[0][0]=3;  intordvum[0][1]=3;  intordvum[0][2]=3;
      intordvlm[0][0]=3;  intordvlm[1][0]=3;  intordvlm[2][0]=3;
      
      if (Cp->savemode==0){
	nipu[0][0]=27;  nipu[0][1]=27;  nipu[0][2]=27;
	nipl[0][0]=27;  nipl[1][0]=27;  nipl[2][0]=27;
      }
      if (Cp->savemode==1){
	nipu[0][0]=27;  nipu[0][1]=0;  nipu[0][2]=0;
	nipl[0][0]=27;  nipl[1][0]=0;  nipl[2][0]=0;
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;  dofe[2]=nnetp;
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown type of transported matter is required");
      fprintf (stderr,"\n in function hexahedc::hexahedc (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }

  //  total number of integration points
  tnipu=0;
  for (i=0;i<mnb;i++){
    for (j=0;j<ntm;j++)
      tnipu+=nipu[i][j];
  }
  tnipl=0;
  for (i=0;i<mnb;i++){
    for (j=0;j<ntm;j++)
      tnipl+=nipl[i][j];
  }

}

hexahedc::~hexahedc (void)
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
  delete [] nipu;
  delete [] nipl;
  delete [] intordvum;
  delete [] intordvlm;
  delete [] dofe;
  delete [] mordering;
}

/*
  function initializes order of numerical integration, number of integration points,
  
  @param eid - element id
  
*/
/*
void hexahedc::eleminit (long eid)
{
  long ii,jj;

  Ct->elements[eid].nb = mnb;
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
*/


/**
   function computes values at integration points from nodal values
   
   @param eid - element id

   05/06/2018 TKr
*/

void hexahedc::intpointval (long eid)
{
  long i,j,k,kk,ii,jj,ipp;
  double xi,eta,zeta,val;
  vector r(ASTCKVEC(tndofe)),t(ASTCKVEC(nnetp)),l(ASTCKVEC(nnetp)),gp,w;
  
  elemvalues (eid, r);
  
  
  for (k=0;k<ntm;k++){
    
    for (i=0;i<dofe[k];i++){
      t[i]=r[tordering[k][i]-1];
    }
    
    for (i=0;i<dofe[k];i++){
      l[i]=r[tordering[k][i]-1];
    }
    
    ii=0;  jj=0;
    for (ii=0;ii<mnb;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
        
        //upper int. points
        reallocv (RSTCKVEC(intordvum[ii][jj],gp));
	reallocv (RSTCKVEC(intordvum[ii][jj],w));
        gauss_points (gp.a,w.a,intordvum[ii][jj]);
        
        ipp=Ct->elements[eid].ippu[ii][jj];
        for (i=0;i<intordvum[ii][jj];i++){
          xi=gp[i];
	  for (j=0;j<intordvum[ii][jj];j++){
	    eta=gp[j];
	    for(kk=0;kk<intordvum[ii][jj];kk++){
	      zeta=gp[kk];
	      
	      if (Cp->bb==lin_lin){
		val = Lht->approx (xi,eta,zeta,t);
	      }
	      if (Cp->bb==quad_lin){
		val = Lht->approx (xi,eta,zeta,t);
	      }
	      if (Cp->bb==quad_quad){
		val = Qht->approx (xi,eta,zeta,t);
	      }
	      
	      Cmu->ip[ipp].av[k]=val;
	      ipp++;
	    }
	  }
	}
        destrv (gp);  destrv (w);

        //lower int. points
        reallocv (RSTCKVEC(intordvlm[jj][ii],gp));
        reallocv (RSTCKVEC(intordvlm[jj][ii],w));
        gauss_points (gp.a,w.a,intordvlm[jj][ii]);
        
        ipp=Ct->elements[eid].ippl[jj][ii];
        for (i=0;i<intordvlm[jj][ii];i++){
          xi=gp[i];
          for (j=0;j<intordvlm[jj][ii];j++){
	    eta=gp[j];
	    for (kk=0;kk<intordvlm[jj][ii];kk++){
	      zeta=gp[kk];
	      
	      if (Cp->bb==lin_lin){
		val = Lht->approx (xi,eta,zeta,l);
	      }
	      if (Cp->bb==quad_lin){
		val = Lht->approx (xi,eta,zeta,l);
	      }
	      if (Cp->bb==quad_quad){
		val = Qht->approx (xi,eta,zeta,l);
	      }
	      
	      Cml->ip[ipp].av[k]=val;
	      ipp++;
	    }
	  }
        }
        if (Cp->savemode==1)  break;
      }
      if (Cp->savemode==1)  break;
    }
  }
}



/**
   function computes gradients in integration points from nodal values
   
   @param eid - element id
   
   05/06/2018 TKr

*/
void hexahedc::intpointgrad (long eid)
{

  long i,j,ii,jj,k,kk,ipp;
  double xi,eta,zeta,jac;
  vector x(nnetp),y(nnetp),z(nnetp),r(tndofe),t(nnetp),l(nnetp),gp,w,grad(3);
  matrix gm(3,nnetp);
  
  Ct->give_node_coord3d (x,y,z,eid);
  
  elemvalues (eid, r);
  
  for (k=0;k<Tp->ntm;k++){
    
    for (i=0;i<dofe[k];i++){
      t[i]=r[tordering[k][i]-1];
    }
    
    for (i=0;i<dofe[k];i++){
      l[i]=r[tordering[k][i]-1];
    }    
    
    ii=0;  jj=0;
    for (ii=0;ii<mnb;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
        
        //upper int. points
        allocv (intordvum[ii][jj],gp);
        allocv (intordvum[ii][jj],w);
        gauss_points (gp.a,w.a,intordvum[ii][jj]);
        
        ipp=Ct->elements[eid].ippu[ii][jj];
        for (i=0;i<intordvum[ii][jj];i++){
          xi=gp[i];
	  for (j=0;j<intordvum[ii][jj];j++){
	    eta=gp[j];
	    for (kk=0;kk<intordvum[ii][jj];kk++){
	      zeta=gp[kk];
	      
	      //  matrix of gradients
	      if (Cp->bb==lin_lin){
		Lht->grad_matrix (gm,x,y,z,xi,eta,zeta,jac);
	      }
	      if (Cp->bb==quad_lin){
		Lht->grad_matrix (gm,x,y,z,xi,eta,zeta,jac);
	      }
	      if (Cp->bb==quad_quad){
		Qht->grad_matrix (gm,x,y,z,xi,eta,zeta,jac);
	      }
	      
	      mxv (gm,t,grad);
	      Cmu->storegrad_cmu (k,ipp,grad);
	      ipp++;
	    }
	  }
	}
        destrv (gp);  destrv (w);
        
        //lower int. points
        allocv (intordvlm[jj][ii],gp);
        allocv (intordvlm[jj][ii],w);
        gauss_points (gp.a,w.a,intordvlm[jj][ii]);
        
        ipp=Ct->elements[eid].ippl[jj][ii];
        for (i=0;i<intordvlm[jj][ii];i++){
          xi=gp[i];
	  for (j=0;j<intordvlm[jj][ii];j++){
	    eta=gp[j];
	    for (kk=0;kk<intordvlm[jj][ii];kk++){
	      zeta=gp[kk];
	    
	    //  matrix of gradients
	    if (Cp->bb==lin_lin){
	      Lht->grad_matrix (gm,x,y,z,xi,eta,zeta,jac);
	    }
	    if (Cp->bb==quad_lin){
	      Lht->grad_matrix (gm,x,y,z,xi,eta,zeta,jac);
	    }
	    if (Cp->bb==quad_quad){
	      Qht->grad_matrix (gm,x,y,z,xi,eta,zeta,jac);
	    }
	    
	    mxv (gm,l,grad);
	    Cml->storegrad_cml (k,ipp,grad);
	    ipp++;
	    }
	  }
        }
        destrv (gp);  destrv (w);

	if (Cp->savemode==1)  break;
      }
      if (Cp->savemode==1)  break;
    } 
  }
}



/**
   function computes strains in integration points of element

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param x,y,z - arrays with node coordinates
   @param r - %vector of nodal displacements

   05/06/2018 TKr
*/
void hexahedc::mainip_strains (long eid,long ri,long ci,vector &x,vector &y,vector &z,vector &r)
{
  long lri,lci;
  long i,j,k,ii,ipp;
  double xi,eta,zeta,jac;
  vector gp,w,eps;
  matrix gm;

  if (Cp->savemode==0){
    lri=ri;  lci=ci;
  }
  if (Cp->savemode==1){
    lri=0;  lci=0;
  }
  
  for (ii=0;ii<mnb;ii++){
    if (intordvlm[ii][ii]==0)  continue;
    
    reallocv (RSTCKVEC(intordvlm[lri][lci],gp));
    reallocv (RSTCKVEC(intordvlm[lri][lci],w));
    reallocv (RSTCKVEC(mncomp[lri],eps));
    reallocm (RSTCKMAT(mncomp[lri],mndofe,gm));
    
    gauss_points (gp.a,w.a,intordvlm[lri][lci]);
    
    //lower int. points
    ipp=Ct->elements[eid].ippl[lri][lci];
    
    for (i=0;i<intordvlm[lri][lci];i++){
      xi=gp[i];
      for (j=0;j<intordvlm[lri][lci];j++){
	eta=gp[j];
	for (k=0;k<intordvlm[lri][lci];k++){
	  zeta=gp[k];
	  
	  if (Cp->bb==lin_lin){
	    Lhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  }
	  if (Cp->bb==quad_lin){
	    Qhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  }
	  if (Cp->bb==quad_quad){
	    Qhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  }
	  
	  mxv (gm,r,eps);
	  
	  Cml->storestrain_cml (ipp,0,eps);
	  ipp++;
	}
      }
    }

    //upper int. points - musi byt opacne ri,ci!
    ipp=Ct->elements[eid].ippu[lci][lri];
    for (i=0;i<intordvum[lci][lri];i++){
      xi=gp[i];
      for (j=0;j<intordvum[lci][lri];j++){
	eta=gp[j];
	for (k=0;k<intordvum[lci][lri];k++){
	  zeta=gp[k];
	  
	  if (Cp->bb==lin_lin){
	    Lhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  }
	  if (Cp->bb==quad_lin){
	    Qhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  }
	  if (Cp->bb==quad_quad){
	    Qhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  }
	  
	  mxv (gm,r,eps);
	  
	  Cmu->storestrain_cmu (ipp,0,eps);
	  ipp++;
	}
      }
    }
  }
}


/**
   function computes resulting strains in integration points of element

   @param eid - element id

   TKr, 22.3.2004
*/
void hexahedc::res_mainip_strains (long eid)
{
  long i;
  vector aux,x(ASTCKVEC(nnemp)),y(ASTCKVEC(nnemp)),z(ASTCKVEC(nnemp)),r(ASTCKVEC(mndofe));
  ivector nodes(ASTCKIVEC(nnemp));
  matrix tmat;
  
  Ct->give_node_coord3d (x,y,z,eid);
  Ct->give_elemnodes (eid,nodes);
  
  eldispl (0,eid,r.a);
    
  //  transformation of displacement vector
  long transf = Mt->locsystems (nodes);
  if (transf>0){
    reallocv (RSTCKVEC(mndofe,aux));
    reallocm (RSTCKMAT(mndofe,mndofe,tmat));

    if (Cp->bb==lin_lin)
      Lhex->transf_matrix (nodes,tmat);
    if (Cp->bb==quad_lin)
      Qhex->transf_matrix (nodes,tmat);
    if (Cp->bb==quad_quad)
      Qhex->transf_matrix (nodes,tmat);

    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  for (i=0;i<Tp->ntm;i++){
    mainip_strains (eid,i,0,x,y,z,r);
  }
}




/**
   function computes one block of the upper coupling stiffness-conductivity %matrix
   
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param vm - block of coupling %matrix
   
   JK, 10.4.2003
*/
void hexahedc::upper_cond_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,j,k,ipp;
  double xi,eta,zeta,jac;
  vector x(nnemp),y(nnemp),z(nnemp),w(intordvum[ri][ci]),gp(intordvum[ri][ci]);
  matrix gm(mncomp[ri],mndofe),d(mncomp[ri],1),n(1,dofe[ci]);
  
  //  integration points and weights
  gauss_points (gp.a,w.a,intordvum[ri][ci]);
  
  fillm (0.0,vm);
  
  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  
  //  number of integration point
  if (Cp->savemode==0){
    ipp=Ct->elements[eid].ippu[ri][ci];
  }
  if (Cp->savemode==1){
    ipp=Ct->elements[eid].ippu[0][0];
  }
  
  for (i=0;i<intordvum[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordvum[ri][ci];j++){
      eta=gp[j];
      for (k=0;k<intordvum[ri][ci];k++){
	zeta=gp[k];
	
	
	if (Cp->bb==lin_lin){
	  //  geometric matrix
	  Lhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  //Lhex->geom_matrix_block (gm,ri,x,y,z,xi,eta,zeta,jac);
	  //  matrix of approxiamtion functions
	  Lht->bf_matrix (n,xi,eta,zeta);
	}
	if (Cp->bb==quad_lin){
	  //  geometric matrix
	  Qhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  //Qhex->geom_matrix_block (gm,ri,x,y,z,xi,eta,zeta,jac);
	  //  matrix of approxiamtion functions
	  Lht->bf_matrix (n,xi,eta,zeta);
	}
	if (Cp->bb==quad_quad){
	  //  geometric matrix
	  Qhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  //Qhex->geom_matrix_block (gm,ri,x,y,z,xi,eta,zeta,jac);
	  //  matrix of approxiamtion functions
	  Qht->bf_matrix (n,xi,eta,zeta);
	}
	
	
	Cmu->matcond (d,ipp,ri,ci);
	
	jac*=w[i]*w[j]*w[k];

	//  contribution to the conductivity matrix of the element
	bdbjac (vm,gm,d,n,jac);
	
	ipp++;
      }
    }
  }

}


/**
   function computes one block of the lower coupling stiffness-conductivity %matrix
   
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param vm - coupling %matrix
   
   JK, 10.4.2003
*/
void hexahedc::lower_cond_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,j,k,ipp;
  double xi,eta,zeta,jac;
  vector x(nnemp),y(nnemp),z(nnemp),w(intordvlm[ri][ci]),gp(intordvlm[ri][ci]);
  matrix gm(mncomp[ci],mndofe),d(1,mncomp[ci]),n(1,dofe[ri]);
  
  //  integration points and weights
  gauss_points (gp.a,w.a,intordvlm[ri][ci]);
  
  fillm (0.0,vm);
  
  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  
  //  number of integration point
  if (Cp->savemode==0){
    ipp=Ct->elements[eid].ippl[ri][ci];
  }
  if (Cp->savemode==1){
    ipp=Ct->elements[eid].ippl[0][0];
  }
  
  for (i=0;i<intordvlm[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordvlm[ri][ci];j++){
      eta=gp[j];
      for (k=0;k<intordvlm[ri][ci];k++){
	zeta=gp[k];
	
	if (Cp->bb==lin_lin){
	  //  geometric matrix
	  Lhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  //Lhex->geom_matrix_block (gm,ci,x,y,z,xi,eta,zeta,jac);
	  //  matrix of approxiamtion functions
	  Lht->bf_matrix (n,xi,eta,zeta);
	}
	if (Cp->bb==quad_lin){
	  //  geometric matrix
	  Qhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  //Qhex->geom_matrix_block (gm,ci,x,y,z,xi,eta,zeta,jac);
	  //  matrix of approxiamtion functions
	  Lht->bf_matrix (n,xi,eta,zeta);
	}
	if (Cp->bb==quad_quad){
	  //  geometric matrix
	  Qhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  //Qhex->geom_matrix_block (gm,ci,x,y,z,xi,eta,zeta,jac);
	  //  matrix of approxiamtion functions
	  Qht->bf_matrix (n,xi,eta,zeta);
	}
	
	//  matrix of conductivity of the material
	Cml->matcond (d,ipp,ri,ci);
	
	jac*=w[i]*w[j]*w[k];

	//  contribution to the conductivity matrix of the element
	bdbjac (vm,n,d,gm,jac);
	
	ipp++;
      }
    }
  }
  
}

/**
   function computes one block of upper coupling capacity %matrix
   
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param vm - coupling %matrix
   
   JK, 9.1.2003
*/
void hexahedc::upper_cap_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,j,k,ipp;
  double xi,eta,zeta,jac;
  vector x(nnemp),y(nnemp),z(nnemp),w(intordvum[ri][ci]),gp(intordvum[ri][ci]);
  matrix gm(mncomp[ri],mndofe),d(mncomp[ri],1),n(1,dofe[ci]);
  
  //  integration points and weights
  gauss_points (gp.a,w.a,intordvum[ri][ci]);

  fillm (0.0,vm);

  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);

  //  number of integration point
  if (Cp->savemode==0){
    ipp=Ct->elements[eid].ippu[ri][ci];
  }
  if (Cp->savemode==1){
    ipp=Ct->elements[eid].ippu[0][0];
  }
    
  for (i=0;i<intordvum[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordvum[ri][ci];j++){
      eta=gp[j];
      for (k=0;k<intordvum[ri][ci];k++){
	zeta=gp[k];
    

	if (Cp->bb==lin_lin){
	  //  geometric matrix
	  Lhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  //Lhex->geom_matrix_block (gm,ri,x,y,z,xi,eta,zeta,jac);
	  //  matrix of approxiamtion functions
	  Lht->bf_matrix (n,xi,eta,zeta);
	}
	if (Cp->bb==quad_lin){
	  //  geometric matrix
	  Qhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  //Qhex->geom_matrix_block (gm,ri,x,y,z,xi,eta,zeta,jac);
	  //  matrix of approxiamtion functions
	  Lht->bf_matrix (n,xi,eta,zeta);
	}
	if (Cp->bb==quad_quad){
	  //  geometric matrix
	  Qhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  //Qhex->geom_matrix_block (gm,ri,x,y,z,xi,eta,zeta,jac);
	  //  matrix of approxiamtion functions
	  Qht->bf_matrix (n,xi,eta,zeta);
	}

	
	//  capacity matrix of the material
	Cmu->matcap (d,ipp,ri,ci);
	
	jac*=w[i]*w[j]*w[k];
	
	//  contribution to the conductivity matrix of the element
	bdbjac (vm,gm,d,n,jac);
	
	ipp++;
      }
    }
  }

}


/**
   function computes one block of coupling capacity %matrix
   
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param vm - coupling %matrix
   
   JK, 10.4.2003
*/
void hexahedc::lower_cap_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,j,k,ipp;
  double xi,eta,zeta,jac;
  vector x(nnemp),y(nnemp),z(nnemp),w(intordvlm[ri][ci]),gp(intordvlm[ri][ci]);
  matrix gm(mncomp[ci],mndofe),d(1,mncomp[ci]),n(1,dofe[ri]);
  
  //  integration points and weights
  gauss_points (gp.a,w.a,intordvlm[ri][ci]);

  fillm (0.0,vm);

  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  
  //  number of integration point
  if (Cp->savemode==0){
    ipp=Ct->elements[eid].ippl[ri][ci];
  }
  if (Cp->savemode==1){
    ipp=Ct->elements[eid].ippl[0][0];
  }
  
  
  for (i=0;i<intordvlm[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordvlm[ri][ci];j++){
      eta=gp[j];
      for (k=0;k<intordvlm[ri][ci];k++){
	zeta=gp[k];
    
	if (Cp->bb==lin_lin){
	  //  geometric matrix
	  Lhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  //Lhex->geom_matrix_block (gm,ci,x,y,z,xi,eta,zeta,jac);
	  //  matrix of approxiamtion functions
	  Lht->bf_matrix (n,xi,eta,zeta);
	}
	if (Cp->bb==quad_lin){
	  //  geometric matrix
	  Qhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  //Qhex->geom_matrix_block (gm,ci,x,y,z,xi,eta,zeta,jac);
	  //  matrix of approxiamtion functions
	  Lht->bf_matrix (n,xi,eta,zeta);
	}
	if (Cp->bb==quad_quad){
	  //  geometric matrix
	  Qhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  //Qhex->geom_matrix_block (gm,ci,x,y,z,xi,eta,zeta,jac);
	  //  matrix of approxiamtion functions
	  Qht->bf_matrix (n,xi,eta,zeta);
	}
	
	
	//  capacity matrix of the material
	Cml->matcap (d,ipp,ri,ci);
	
	jac*=w[i]*w[j]*w[k];
	
	//  contribution to the conductivity matrix of the element
	bdbjac (vm,n,d,gm,jac);
	
	ipp++;
      }
    }
  }
  
}





/**
   function assembles the upper coupling stiffness-conductivity %matrix
   
   @param eid - element id
   @param vm - coupling %matrix
   
   26.3.2004, JK
*/
void hexahedc::res_upper_cond_coup_matrix (long eid,matrix &vm)
{
  long ii,jj,*ccn;
  matrix lvm;
  
  fillm(0.0,vm);
  
  for (ii=0;ii<ntm;ii++){
    //  column code numbers
    ccn = new long [dofe[ii]];
    //  submatrix
    allocm (mndofe,dofe[ii],lvm);
    
    //  code numbers of required components
    if (Cp->bb==lin_lin || Cp->bb==quad_lin)
      Lht->codnum (ccn,ii);
    if (Cp->bb==quad_quad)
      Qht->codnum (ccn,ii);
    
    for (jj=0;jj<mnb;jj++){
      //  block of upper coupling stiffness-conductivity matrix
      upper_cond_coup_matrix (eid,jj,ii,lvm);
      //  localization of the block into the matrix
      mat_localize (vm,lvm,mordering,ccn);
    }

    delete [] ccn;
    destrm (lvm);
  }
}

/**
   function assembles the lower coupling stiffness-conductivity %matrix
   
   @param eid - element id
   @param vm - %matrix
   
   26.3.2004, JK
*/
void hexahedc::res_lower_cond_coup_matrix (long eid,matrix &vm)
{
  long ii,jj,*rcn;
  matrix lvm;
  
  
  fillm(0.0,vm);
  
  for (ii=0;ii<ntm;ii++){
    //  row code numbers
    rcn = new long [dofe[ii]];
    //  submatrix
    allocm (dofe[ii],mndofe,lvm);
    
    //  code numbers of required components
    if (Cp->bb==lin_lin || Cp->bb==quad_lin)
      Lht->codnum (rcn,ii);
    if (Cp->bb==quad_quad)
      Qht->codnum (rcn,ii);
    
    for (jj=0;jj<mnb;jj++){
      //  block of lower coupling stiffness-conductivity matrix
      lower_cond_coup_matrix (eid,ii,jj,lvm);
      //  localization of the block into the matrix
      mat_localize (vm,lvm,rcn,mordering);
    }

    destrm (lvm);
    delete [] rcn;
  }
}


/**
   function assembles the upper coupling capacity %matrix
   
   @param eid - element id
   @param vm - coupling %matrix
   
   26.3.2004, JK
*/
void hexahedc::res_upper_cap_coup_matrix (long eid,matrix &vm)
{
  long ii,jj,*ccn;
  matrix lvm;
  
  fillm(0.0,vm);
  
  for (ii=0;ii<ntm;ii++){
    //  allocation of one block of the matrix
    allocm (mndofe,dofe[ii],lvm);
    //  column code numbers
    ccn = new long [dofe[ii]];
    
    //  code numbers of required components
    if (Cp->bb==lin_lin || Cp->bb==quad_lin)
      Lht->codnum (ccn,ii);
    if (Cp->bb==quad_quad)
      Qht->codnum (ccn,ii);
    
    for (jj=0;jj<mnb;jj++){
      //  block of upper coupling capacity matrix
      upper_cap_coup_matrix (eid,jj,ii,lvm);
      //  localization of the block into the matrix
      mat_localize (vm,lvm,mordering,ccn);
    }

    destrm (lvm);
    delete [] ccn;
  }
}

/**
   function assembles the lower coupling capacity %matrix
   
   @param eid - element id
   @param vm - coupling %matrix
   
   26.3.2004, JK
*/
void hexahedc::res_lower_cap_coup_matrix (long eid,matrix &vm)
{
  long ii,jj,*rcn;
  matrix lvm;
  
  fillm(0.0,vm);
  
  for (ii=0;ii<ntm;ii++){
    //  allocation of one block of the matrix
    allocm (dofe[ii],mndofe,lvm);
    //  row code numbers
    rcn = new long [dofe[ii]];

    //  code numbers of required components
    if (Cp->bb==lin_lin || Cp->bb==quad_lin)
      Lht->codnum (rcn,ii);
    if (Cp->bb==quad_quad)
      Qht->codnum (rcn,ii);

    for (jj=0;jj<mnb;jj++){
      //  block of upper coupling stiffness-conductivity matrix
      lower_cap_coup_matrix (eid,ii,jj,lvm);
      //  localization of the block into the matrix
      mat_localize (vm,lvm,rcn,mordering);
    }
    
    destrm (lvm);
    delete [] rcn;
  }
}


/**
   function computes coupling %vector of 3D problems
   
   @param vm - coupling %vector
   @param nodval - %vector of nodal values
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   
   TKr, 26.4.2007
*/
void hexahedc::upper_cond_coup_vector (vector &tvm,vector &nodval,long eid,long ri,long ci)
{
  long i,j,k,ipp;
  double xi,eta,zeta,jac;
  vector x(nnemp),y(nnemp),z(nnemp),w(intordvum[ri][ci]),gp(intordvum[ri][ci]),av(mndofe);
  matrix gm(mncomp[ri],mndofe),d(mncomp[ri],1),n(1,dofe[ci]),vm(mndofe,dofe[ci]);
  
  //  integration points and weights
  gauss_points (gp.a,w.a,intordvum[ri][ci]);
  
  fillm (0.0,vm);
  
  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  
  //  number of integration point
  if (Cp->savemode==0){
    ipp=Ct->elements[eid].ippu[ri][ci];
  }
  if (Cp->savemode==1){
    ipp=Ct->elements[eid].ippu[0][0];
  }
  
  for (i=0;i<intordvum[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordvum[ri][ci];j++){
      eta=gp[j];
      for (k=0;k<intordvum[ri][ci];k++){
	zeta=gp[k];
	
	
	if (Cp->bb==lin_lin){
	  //  geometric matrix
	  Lhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  //Lhex->geom_matrix_block (gm,ri,x,y,z,xi,eta,zeta,jac);
	  //  matrix of approxiamtion functions
	  Lht->bf_matrix (n,xi,eta,zeta);
	}
	if (Cp->bb==quad_lin){
	  //  geometric matrix
	  Qhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  //Qhex->geom_matrix_block (gm,ri,x,y,z,xi,eta,zeta,jac);
	  //  matrix of approxiamtion functions
	  Lht->bf_matrix (n,xi,eta,zeta);
	}
	if (Cp->bb==quad_quad){
	  //  geometric matrix
	  Qhex->geom_matrix (gm,x,y,z,xi,eta,zeta,jac);
	  //Qhex->geom_matrix_block (gm,ri,x,y,z,xi,eta,zeta,jac);
	  //  matrix of approxiamtion functions
	  Qht->bf_matrix (n,xi,eta,zeta);
	}
	
	//  matrix of conductivity of the material
	Cmu->volume_rhs2 (d,ipp,ri,ci);
	
	jac*=w[i]*w[j]*w[k];
	
	//  contribution to the conductivity matrix of the element
	bdbjac (vm,gm,d,n,jac);
	
	ipp++;
      }
    }
  }
  
  mxv (vm,nodval,av);
  addv (av,tvm,tvm);
  
}




/**
   function computes
   
   @param f - 
   @param eid - element id

   TKr, 26.4.2007
*/
void hexahedc::res_upper_cond_coup_vector (vector &f,long eid)
{
  long ii,jj;
  long *lcn;
  vector r(ASTCKVEC(tndofe));
  vector lr;

  //  all initial transport values
  initialvalues (eid, r);
  
  nullv (f);
  
  for (ii=0;ii<ntm;ii++){
    //  code numbers for unknowns decribing one medium
    lcn = new long [dofe[ii]];
    //  initial values of one medium
    reallocv (dofe[ii],lr);
    
    //  ordering on element
    if (Cp->bb==lin_lin || Cp->bb==quad_lin)
      Lht->codnum (lcn,ii);
    else
      Qht->codnum (lcn,ii);
    
    //  extraction of initial values of one medium
    globloc (r.a,lr.a,lcn,dofe[ii]);
    
    for (jj=0;jj<mnb;jj++){
      //  contribution to vector
      upper_cond_coup_vector (f,lr,eid,jj,ii);
    }
    delete [] lcn;
  }
}



/*
   function computes internal forces

   @param lcid - number of load case
   @param eid - element id
   @param ifor - vector of internal forces
   
   4.7.2004, JK

void hexahedc::upper_internal_forces (long ri,long eid,vector &ifor)
{
  long i,j,k,l,ii,ipp;
  double xi,eta,zeta,jac;
  ivector cn(ndofe);
  vector x(nne),y(nne),z(nne),w(),gp,r(ndofe),eps(tncomp),sig,contr(ndofe);
  matrix gm(mncomp,mndofe);

  Mt->give_node_coord3d (x,y,z,eid);
  Mt->give_code_numbers (eid,cn.a);
  eldispl (0,r.a,cn.a,ndofe);
  
  fillv (0.0,ifor);

  for (ii=0;ii<nb;ii++){
    if (intordsm[ii][ii]==0)  continue;

    allocv (intordsm[ii][ii],gp);
    allocv (intordsm[ii][ii],w);
    allocv (ncomp[ii],sig);
    
    gauss_points (gp.a,w.a,intordsm[ii][ii]);
    ipp=Mt->elements[eid].ipp[ii][ii];
    
    for (i=0;i<intordvum[0][ii];i++){
      xi=gp[i];
      for (j=0;j<intordsm[ii][ii];j++){
	eta=gp[j];
	for (k=0;k<intordsm[ii][ii];k++){
	  zeta=gp[k];
	  
	  
	  appstrain (lcid,eid,xi,eta,zeta,0,tncomp,eps);
	  Mm->storestrain (lcid,ipp,eps);
	  
	  Mm->computenlstresses (ipp);
	  
	  Mm->givestress (lcid,ipp,cncomp[ii],ncomp[ii],sig);
	  
	  geom_matrix_block (gm,ii,x,y,z,xi,eta,zeta,jac);
	  mtxv (gm,sig,contr);
	  
	  cmulv (jac*w[i]*w[j]*w[k],contr);
	  
	  for (l=0;l<contr.n;l++){
	    ifor[l]+=contr[l];
	  }
	  
	  ipp++;
	}
      }
    }
    destrv (sig);  destrm (gm);  destrv (w);  destrv (gp);
  }
}

void quadhex::res_internal_forces (long lcid,long eid,vector &ifor)
{
  internal_forces (lcid,eid,ifor);
}
*/



/**
   function copies mechanical data from MEFEL into integration points
   of coupled part in METR
   
   @param eid - element id
   
   JK, 29.10.2004
*/
void hexahedc::mefel_metr (long eid)
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
	Cmu->ip[ippcu].strain[i] = Mm->ip[ippm].strain[i];
	Cmu->ip[ippcu].stress[i] = Mm->ip[ippm].stress[i];
      }
      ippm++;  ippcu++;
    }
    
    //  number of the first integration point on mechanical element
    //  one block formulation of mechanical problem is assumed
    ippm = Mt->elements[eid].ipp[0][0];
    for (j=0;j<nipl[0][0];j++){
      for (i=0;i<tnmcomp;i++){
	Cml->ip[ippcl].strain[i] = Mm->ip[ippm].strain[i];
	Cml->ip[ippcl].stress[i] = Mm->ip[ippm].stress[i];
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
void hexahedc::trfel_metr (long eid)
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
