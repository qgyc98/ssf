#include "axiquadc.h"
#include "axisymlq.h"
#include "axisymqq.h"
#include "globmat.h"
#include "global.h"
#include "mechtop.h"
#include "globalt.h"
#include "globalc.h"
#include "element.h"
#include "intpoints.h"
#include "globmatt.h"


axiquadc::axiquadc (void)
{
  long i,j;

  if (Cp->bb==lin_lin){
    //  number of blocks of the mechanical element
    mnb=Asymlq->nb;
    //  number of mechanical DOFs
    mndofe=Asymlq->ndofe;
    //  number of strain/stress components
    tnmcomp=Asymlq->tncomp;
    //  array of numbers of components in strain vectors
    mncomp=Asymlq->ncomp;
    //  number of nodes in mechanical problems
    nnemp = Asymlq->nne;
    ssst = Asymlq->ssst;
    //  number of nodes in transport problems
    nnetp = Lqat->nne;
    //  number of transport DOFs
    tndofe = Lqat->ndofe;
  }

  if (Cp->bb==quad_lin){
    //  number of blocks of the mechanical element
    mnb=Asymqq->nb;
    //  number of mechanical DOFs
    mndofe=Asymqq->ndofe;
    //  number of strain/stress components
    tnmcomp=Asymqq->tncomp;
    //  array of numbers of components in strain vectors
    mncomp=Asymqq->ncomp;
    //  number of nodes in mechanical problems
    nnemp = Asymqq->nne;
    ssst = Asymqq->ssst;
    //  number of nodes in transport problems
    nnetp = Lqat->nne;
    //  number of transport DOFs
    tndofe = Lqat->ndofe;
  }

  if (Cp->bb==quad_quad){
    //  number of blocks of the mechanical element
    mnb=Asymqq->nb;
    //  number of mechanical DOFs
    mndofe=Asymqq->ndofe;
    //  number of strain/stress components
    tnmcomp=Asymqq->tncomp;
    //  array of numbers of components in strain vectors
    mncomp=Asymqq->ncomp;
    //  number of nodes in mechanical problems
    nnemp = Asymqq->nne;
    ssst = Asymqq->ssst;
    //  number of nodes in transport problems
    nnetp = Qqat->nne;
    //  number of transport DOFs
    tndofe = Qqat->ndofe;
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
      nipu[0][0]=4;  nipl[0][0]=4;
      dofe[0]=nnetp;  

      tordering[0][0]=1;  tordering[0][1]=2; tordering[0][2]=3;  tordering[0][3]=4;
      break;
    }
    case twomediacoup:{
      intordvum[0][0]=2;  intordvum[0][1]=2;
      intordvlm[0][0]=2;  intordvlm[1][0]=2;

      tordering[0][0]=1;  tordering[0][1]=3; tordering[0][2]=5;  tordering[0][3]=7;
      tordering[1][0]=2;  tordering[1][1]=4; tordering[1][2]=6;  tordering[1][3]=8;

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
      
      tordering[0][0]=1;  tordering[0][1]=4; tordering[0][2]=7;  tordering[0][3]=10;
      tordering[1][0]=2;  tordering[1][1]=5; tordering[1][2]=8;  tordering[1][3]=11;
      tordering[2][0]=3;  tordering[2][1]=6; tordering[2][2]=9;  tordering[2][3]=12;
      
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
      fprintf (stderr,"\n in function axiquadc::axiquadc (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }

  if (Cp->bb==quad_lin){
    switch (Tp->tmatt){
    case onemedium:{
      intordvum[0][0]=2;  intordvlm[0][0]=2;
      nipu[0][0]=4;  nipl[0][0]=4;
      dofe[0]=nnetp;  

      tordering[0][0]=1;  tordering[0][1]=2; tordering[0][2]=3;  tordering[0][3]=4;
      break;
    }
    case twomediacoup:{
      intordvum[0][0]=2;  intordvum[0][1]=2;
      intordvlm[0][0]=2;  intordvlm[1][0]=2;

      tordering[0][0]=1;  tordering[0][1]=3; tordering[0][2]=5;  tordering[0][3]=7;
      tordering[1][0]=2;  tordering[1][1]=4; tordering[1][2]=6;  tordering[1][3]=8;
      
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
      
      tordering[0][0]=1;  tordering[0][1]=4; tordering[0][2]=7;  tordering[0][3]=10;
      tordering[1][0]=2;  tordering[1][1]=5; tordering[1][2]=8;  tordering[1][3]=11;
      tordering[2][0]=3;  tordering[2][1]=6; tordering[2][2]=9;  tordering[2][3]=12;
      
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
      fprintf (stderr,"\n in function axiquadc::axiquadc (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
  
  if (Cp->bb==quad_quad){
    switch (Tp->tmatt){
    case onemedium:{
      intordvum[0][0]=3;  intordvlm[0][0]=3;
      nipu[0][0]=9;  nipl[0][0]=9;
      dofe[0]=nnetp;  

      tordering[0][0]=1;  tordering[0][1]=2; tordering[0][2]=3;  tordering[0][3]=4; tordering[0][4]=5;  tordering[0][5]=6; tordering[0][6]=7;  tordering[0][7]=8;
      break;
    }
    case twomediacoup:{
      intordvum[0][0]=3;  intordvum[0][1]=3;
      intordvlm[0][0]=3;  intordvlm[1][0]=3;

      tordering[0][0]=1;  tordering[0][1]=3; tordering[0][2]=5;  tordering[0][3]=7; tordering[0][4]=9;  tordering[0][5]=11; tordering[0][6]=13;  tordering[0][7]=15;
      tordering[1][0]=2;  tordering[1][1]=4; tordering[1][2]=6;  tordering[1][3]=8; tordering[1][4]=10; tordering[1][5]=12; tordering[1][6]=14;  tordering[1][7]=16;
      
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
      
      tordering[0][0]=1;  tordering[0][1]=4; tordering[0][2]=7;  tordering[0][3]=10; tordering[0][4]=13;  tordering[0][5]=16; tordering[0][6]=19;  tordering[0][7]=22;
      tordering[1][0]=2;  tordering[1][1]=5; tordering[1][2]=8;  tordering[1][3]=11; tordering[1][4]=14;  tordering[1][5]=17; tordering[1][6]=20;  tordering[1][7]=23;
      tordering[2][0]=3;  tordering[2][1]=6; tordering[2][2]=9;  tordering[2][3]=12; tordering[2][4]=15;  tordering[2][5]=18; tordering[2][6]=21;  tordering[2][7]=24;    
      
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
      fprintf (stderr,"\n in function axiquadc::axiquadc (file %s, line %d).\n",__FILE__,__LINE__);
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
  for (i=0;i<ntm;i++){
    for (j=0;j<mnb;j++)
      tnipl+=nipl[i][j];
  }

}

/**
   destructor
*/
axiquadc::~axiquadc (void)
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
  delete [] mordering;
}

/**
   allocation and initiation of integration orders and numbers of integration points
   
   @param eid - element id
   
   JK, 24.10.2004
*/
/*
void axiquadc::eleminit (long eid)
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
   function assembles blocks of stiffness %matrix of material
   
   @param ri - row index
   @param ci - column index
   @param d - stiffness %matrix of material
   @param dd - required block of stiffness %matrix of material
*/
void axiquadc::dmatblockcol (long ri,long ci,matrix &d, matrix &dd)
{
  fillm (0.0,dd);

  if (mnb==1){
    dd[0][0]=d[0][0];
    dd[1][0]=d[1][0];
    dd[2][0]=d[2][0];
    dd[3][0]=d[3][0];
  }
  if (mnb==2){
    if (ri==0 && ci==0){
      dd[0][0]=d[0][0];
      dd[1][0]=d[1][0];
    }
    if (ri==1 && ci==0){
      dd[0][0]=d[2][0];
    }
    if (ri==2 && ci==0){
      dd[0][0]=d[3][0];
    }
  }
}

/**
   function assembles blocks of stiffness %matrix of material
   
   @param ri - row index
   @param ci - column index
   @param d - stiffness %matrix of material
   @param dd - required block of stiffness %matrix of material
*/
void axiquadc::dmatblockrow (long ri,long ci,matrix &d, matrix &dd)
{
  fillm (0.0,dd);
  
  if (mnb==1){
    dd[0][0]=d[0][0];
    dd[0][1]=d[0][1];
    dd[0][2]=d[0][2];
    dd[0][3]=d[0][3];
  }
  if (mnb==2){
    if (ri==0 && ci==0){
      dd[0][0]=d[0][0];
      dd[0][1]=d[0][1];
    }
    if (ri==0 && ci==1){
      dd[0][0]=d[0][2];
    }
    if (ri==0 && ci==2){
      dd[0][0]=d[0][3];
    }
  }
}




/**
   function computes values at integration points from nodal values
   
   @param eid - element id

   05/06/2018, TKr according to quadrilatc.cpp
   
*/

void axiquadc::intpointval (long eid)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,val;
  vector r(tndofe),t(nnetp),l(nnetp),gp,w;
  
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
        allocv (intordvum[ii][jj],gp);
        allocv (intordvum[ii][jj],w);
        gauss_points (gp.a,w.a,intordvum[ii][jj]);
        
        ipp=Ct->elements[eid].ippu[ii][jj];
        for (i=0;i<intordvum[ii][jj];i++){
          xi=gp[i];
	  for (j=0;j<intordvum[ii][jj];j++){
	    eta=gp[j];
	    
	    if (Cp->bb==lin_lin){
	      val = Lqat->approx (xi,eta,t);
	    }
	    if (Cp->bb==quad_lin){
	      val = Lqat->approx (xi,eta,t);
	    }
	    if (Cp->bb==quad_quad){
	      val = Qqat->approx (xi,eta,t);
	    }
	    
	    Cmu->ip[ipp].av[k]=val;
	    ipp++;
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
	    
	    if (Cp->bb==lin_lin){
	      val = Lqat->approx (xi,eta,l);
	    }
	    if (Cp->bb==quad_lin){
	      val = Lqat->approx (xi,eta,l);
	    }
	    if (Cp->bb==quad_quad){
	      val = Qqat->approx (xi,eta,l);
	    }
	    
	    Cml->ip[ipp].av[k]=val;
	    ipp++;
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
   function computes gradients in integration points from nodal values
   
   @param eid - element id
   
   05/06/2018, TKr according to quadrilatc.cpp

*/
void axiquadc::intpointgrad (long eid)
{

  long i,j,ii,jj,k,ipp;
  double xi,eta,jac;
  vector x(nnetp),y(nnetp),r(tndofe),t(nnetp),l(nnetp),gp,w,grad(2);
  matrix gm(2,nnetp);
  
  //Ct->give_node_coord2d (x,y,eid);
  Tt->give_node_coord2d (x,y,eid);  

  elemvalues (eid,r);
  
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
	    
	    //  matrix of gradients
	    if (Cp->bb==lin_lin){
	      Lqat->grad_matrix (gm,x,y,xi,eta,jac);
	    }
	    if (Cp->bb==quad_lin){
	      Lqat->grad_matrix (gm,x,y,xi,eta,jac);
	    }
	    if (Cp->bb==quad_quad){
	      Qqat->grad_matrix (gm,x,y,xi,eta,jac);
	    }
	    
	    mxv (gm,t,grad);
	    Cmu->storegrad_cmu (k,ipp,grad);
	    ipp++;
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
	    
	    //  matrix of gradients
	    if (Cp->bb==lin_lin){
	      Lqat->grad_matrix (gm,x,y,xi,eta,jac);
	    }
	    if (Cp->bb==quad_lin){
	      Lqat->grad_matrix (gm,x,y,xi,eta,jac);
	    }
	    if (Cp->bb==quad_quad){
	      Qqat->grad_matrix (gm,x,y,xi,eta,jac);
	    }
	    
	    mxv (gm,l,grad);
	    Cml->storegrad_cml (k,ipp,grad);
	    ipp++;
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
   @param x,y - arrays with node coordinates
   @param r - %vector of nodal displacements

   05/06/2018 TKr
*/
void axiquadc::mainip_strains (long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long lri,lci;
  long i,j,ii,ipp;
  double xi,eta,jac;
  vector gp,w,eps;
  matrix gm;

  if (Cp->savemode==0){
    lri=ri;  lci=ci;
  }
  if (Cp->savemode==1){
    lri=0;  lci=0;
  }
  
  for (ii=0;ii<mnb;ii++){
    //if (intordvlm[ii][ii]==0)  continue;
    
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
	
	if (Cp->bb==lin_lin){
	  Asymlq->geom_matrix (gm,x,y,xi,eta,jac);
	}
	if (Cp->bb==quad_lin){
	  Asymqq->geom_matrix (gm,x,y,xi,eta,jac);
	}
	if (Cp->bb==quad_quad){
	  Asymqq->geom_matrix (gm,x,y,xi,eta,jac);
	}
	
	mxv (gm,r,eps);
	
	Cml->storestrain_cml (ipp,0,eps);
	ipp++;
     }
    }

    //upper int. points - musi byt opacne ri,ci!
    ipp=Ct->elements[eid].ippu[lci][lri];
    for (i=0;i<intordvum[lci][lri];i++){
      xi=gp[i];
      for (j=0;j<intordvum[lci][lri];j++){
	eta=gp[j];
	
	if (Cp->bb==lin_lin){
	  Asymlq->geom_matrix (gm,x,y,xi,eta,jac);
	}
	if (Cp->bb==quad_lin){
	  Asymqq->geom_matrix (gm,x,y,xi,eta,jac);
	}
	if (Cp->bb==quad_quad){
	  Asymqq->geom_matrix (gm,x,y,xi,eta,jac);
	}
	
	mxv (gm,r,eps);
	
	Cmu->storestrain_cmu (ipp,0,eps);
	ipp++;
      }
    }
  }
}


/**
   function computes resulting strains in integration points of element

   @param eid - element id

   TKr, 22.3.2004
*/
void axiquadc::res_mainip_strains (long eid)
{
  long i;
  vector aux,x(ASTCKVEC(nnemp)),y(ASTCKVEC(nnemp)),r(ASTCKVEC(mndofe));
  ivector nodes(ASTCKIVEC(nnemp));
  matrix tmat;
  
  Ct->give_node_coord2d (x,y,eid);
  Ct->give_elemnodes (eid,nodes);
  
  eldispl (0,eid,r.a);
    
  //  transformation of displacement vector
  long transf = Mt->locsystems (nodes);
  if (transf>0){
    reallocv (RSTCKVEC(mndofe,aux));
    reallocm (RSTCKMAT(mndofe,mndofe,tmat));

    if (Cp->bb==lin_lin)
      Asymlq->transf_matrix (nodes,tmat);
    if (Cp->bb==quad_lin)
      Asymqq->transf_matrix (nodes,tmat);
    if (Cp->bb==quad_quad)
      Asymqq->transf_matrix (nodes,tmat);

    //locglobtransf (aux,r,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
  }
  
  for (i=0;i<Tp->ntm;i++){
    mainip_strains (eid,i,0,x,y,r);
  }
}





/**
   function computes upper coupling conductivity %matrix of 2D problems
   it is coupling submatrix between the stiffness and conductivity %matrices
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param vm - coupling %matrix
   
   JK, 24.10.2004
*/
void axiquadc::upper_cond_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,j,ipp;
  double xi,eta,jac,r;
  vector x(nnemp),y(nnemp),w(intordvum[ri][ci]),gp(intordvum[ri][ci]);
  matrix gm(mncomp[ri],mndofe),d(tnmcomp,1),dd(mncomp[ri],1),n(1,dofe[ci]);
  
  //  node coordinates of mechanical problem
  Mt->give_node_coord2d (x,y,eid);
  
  //  integration points and weights
  gauss_points (gp.a,w.a,intordvum[ri][ci]);

  fillm (0.0,vm);
  
  if (Cp->savemode==0)
    ipp=Ct->elements[eid].ippu[ri][ci];
  if (Cp->savemode==1)
    ipp=Ct->elements[eid].ippu[0][0];
  
  for (i=0;i<intordvum[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordvum[ri][ci];j++){
      eta=gp[j];
      
      if (Cp->bb==lin_lin){
	//  geometric matrix of axisymmetric quadrilateral mechanical element with linear approximation functions
	Asymlq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	//  matrix of approximation functions of axisymmetric quadrilateral transport element with linear approximation functions
	Lqat->bf_matrix (n,xi,eta);
	//  radius
	r = Asymlq->approx (xi,eta,x);
      }
      if (Cp->bb==quad_lin){
	//  geometric matrix of axisymmetric quadrilateral mechanical element with quadratic approximation functions
	Asymqq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	//  matrix of approximation functions of axisymmetric quadrilateral transport element with linear approximation functions
	Lqat->bf_matrix (n,xi,eta);
	//  radius
	r = Asymqq->approx (xi,eta,x);
      }
      if (Cp->bb==quad_quad){
	//  geometric matrix of axisymmetric quadrilateral mechanical element with quadratic approximation functions
	Asymqq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	//  matrix of approximation functions of axisymmetric quadrilateral transport element with quadratic approximation functions
	Qqat->bf_matrix (n,xi,eta);
	//  radius
	r = Asymqq->approx (xi,eta,x);
      }
      
      
      //  material matrix
      Cmu->matcond (d,ipp,ri,ci);
      
      jac*=w[i]*w[j]*r;
      
      bdbjac (vm,gm,d,n,jac);
      
      ipp++;
    }
  }
}

/**
   function computes lower coupling conductivity %matrix of 2D problems
   it is coupling submatrix between the stiffness and conductivity %matrices
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param vm - coupling %matrix
   
   JK, 24.10.2004
*/
void axiquadc::lower_cond_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,j,ipp;
  double xi,eta,jac,r;
  vector x(nnemp),y(nnemp),w(intordvlm[ri][ci]),gp(intordvlm[ri][ci]);
  matrix gm(tnmcomp,mndofe),d(1,tnmcomp),n(1,nnetp);
  
  //  node coordinates of mechanical problem
  Mt->give_node_coord2d (x,y,eid);
  
  //  integration points and weights
  gauss_points (gp.a,w.a,intordvlm[ri][ci]);

  fillm (0.0,vm);
  
  if (Cp->savemode==0)
    ipp=Ct->elements[eid].ippl[ri][ci];
  if (Cp->savemode==1)
    ipp=Ct->elements[eid].ippl[0][0];
  
  for (i=0;i<intordvlm[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordvlm[ri][ci];j++){
      eta=gp[j];
      
      if (Cp->bb==lin_lin){
	//  geometric matrix of axisymmetric quadrilateral mechanical element with linear approximation functions
	Asymlq->geom_matrix (gm,x,y,xi,eta,jac);
	//  matrix of approximation functions of axisymmetric quadrilateral transport element with linear approximation functions
	Lqat->bf_matrix (n,xi,eta);
	//  radius
	r = Asymlq->approx (xi,eta,x);
      }
      if (Cp->bb==quad_lin){
	//  geometric matrix of axisymmetric quadrilateral mechanical element with quadratic approximation functions
	Asymqq->geom_matrix (gm,x,y,xi,eta,jac);
	//  matrix of approximation functions of axisymmetric quadrilateral transport element with linear approximation functions
	Lqat->bf_matrix (n,xi,eta);
	//  radius
	r = Asymqq->approx (xi,eta,x);
      }
      if (Cp->bb==quad_quad){
	//  geometric matrix of axisymmetric quadrilateral mechanical element with quadratic approximation functions
	Asymqq->geom_matrix (gm,x,y,xi,eta,jac);
	//  matrix of approximation functions of axisymmetric quadrilateral transport element with quadratic approximation functions
	Qqat->bf_matrix (n,xi,eta);
	//  radius
	r = Asymqq->approx (xi,eta,x);
      }
      
      
      //  material matrix
      Cml->matcond (d,ipp,ri,ci);
      
      jac*=w[i]*w[j]*r;
      
      bdbjac (vm,n,d,gm,jac);
      
      ipp++;
    }
  }
}


/**
   function computes upper coupling capacity %matrix of 2D problems
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param vm - coupling %matrix
   
   JK, 24.10.2004
*/
void axiquadc::upper_cap_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,j,ipp;
  double xi,eta,jac,r;
  vector x(nnemp),y(nnemp),w(intordvum[ri][ci]),gp(intordvum[ri][ci]);
  matrix gm(tnmcomp,mndofe),d(tnmcomp,1),n(1,nnetp);
  
  //  node coordinates of mechanical problem
  Mt->give_node_coord2d (x,y,eid);

  fillm (0.0,vm);
  
  //  integration points and weights
  gauss_points (gp.a,w.a,intordvum[ri][ci]);
  
  if (Cp->savemode==0)
    ipp=Ct->elements[eid].ippu[ri][ci];
  if (Cp->savemode==1)
    ipp=Ct->elements[eid].ippu[0][0];
  
  for (i=0;i<intordvum[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordvum[ri][ci];j++){
      eta=gp[j];
      
      if (Cp->bb==lin_lin){
	//  geometric matrix of axisymmetric quadrilateral mechanical element with linear approximation functions
	Asymlq->geom_matrix (gm,x,y,xi,eta,jac);
	//  matrix of approximation functions of axisymmetric quadrilateral transport element with linear approximation functions
	Lqat->bf_matrix (n,xi,eta);
	//  radius
	r = Asymlq->approx (xi,eta,x);
      }
      if (Cp->bb==quad_lin){
	//  geometric matrix of axisymmetric quadrilateral mechanical element with quadratic approximation functions
	Asymqq->geom_matrix (gm,x,y,xi,eta,jac);
	//  matrix of approximation functions of axisymmetric quadrilateral transport element with linear approximation functions
	Lqat->bf_matrix (n,xi,eta);
	//  radius
	r = Asymqq->approx (xi,eta,x);
      }
      if (Cp->bb==quad_quad){
	//  geometric matrix of axisymmetric quadrilateral mechanical element with quadratic approximation functions
	Asymqq->geom_matrix (gm,x,y,xi,eta,jac);
	//  matrix of approximation functions of axisymmetric quadrilateral transport element with quadratic approximation functions
	Qqat->bf_matrix (n,xi,eta);
	//  radius
	r = Asymqq->approx (xi,eta,x);
      }
      
      
      //  material capacity matrix
      Cmu->matcap (d,ipp,ri,ci);
      
      jac*=w[i]*w[j]*r;
      
      bdbjac (vm,gm,d,n,jac);
      
      ipp++;
    }
  }
}

/**
   function computes lower coupling capacity %matrix of 2D problems
   
   @param eid - element id
   @param ri,ci - row and column indices
   @param vm - coupling %matrix
   
   JK, 24.10.2004
*/
void axiquadc::lower_cap_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,j,ipp;
  double xi,eta,jac,r;
  vector x(nnemp),y(nnemp),w(intordvlm[ri][ci]),gp(intordvlm[ri][ci]);
  matrix gm(tnmcomp,mndofe),d(1,tnmcomp),n(1,nnetp);
  
  //  node coordinates of mechanical problem
  Mt->give_node_coord2d (x,y,eid);

  fillm (0.0,vm);
  
  //  integration points and weights
  gauss_points (gp.a,w.a,intordvlm[ri][ci]);
  
  if (Cp->savemode==0)
    ipp=Ct->elements[eid].ippl[ri][ci];
  if (Cp->savemode==1)
    ipp=Ct->elements[eid].ippl[0][0];
  
  for (i=0;i<intordvlm[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordvlm[ri][ci];j++){
      eta=gp[j];
      
      if (Cp->bb==lin_lin){
	//  geometric matrix of axisymmetric quadrilateral mechanical element with linear approximation functions
	Asymlq->geom_matrix (gm,x,y,xi,eta,jac);
	//  matrix of approximation functions of axisymmetric quadrilateral transport element with linear approximation functions
	Lqat->bf_matrix (n,xi,eta);
	//  radius
	r = Asymlq->approx (xi,eta,x);
      }
      if (Cp->bb==quad_lin){
	//  geometric matrix of axisymmetric quadrilateral mechanical element with quadratic approximation functions
	Asymqq->geom_matrix (gm,x,y,xi,eta,jac);
	//  matrix of approximation functions of axisymmetric quadrilateral transport element with linear approximation functions
	Lqat->bf_matrix (n,xi,eta);
	//  radius
	r = Asymqq->approx (xi,eta,x);
      }
      if (Cp->bb==quad_quad){
	//  geometric matrix of axisymmetric quadrilateral mechanical element with quadratic approximation functions
	Asymqq->geom_matrix (gm,x,y,xi,eta,jac);
	//  matrix of approximation functions of axisymmetric quadrilateral transport element with quadratic approximation functions
	Qqat->bf_matrix (n,xi,eta);
	//  radius
	r = Asymqq->approx (xi,eta,x);
      }
      
      
      //  material matrix
      Cml->matcap (d,ipp,ri,ci);
      
      jac*=w[i]*w[j]*r;
      
      bdbjac (vm,n,d,gm,jac);
      
      ipp++;
    }
  }
}


/**
   function assembles upper coupling stiffness-conductivity matrix
   
   @param eid - element id
   @param vm - upper coupling stiffness-conductivity matrix
   
   JK, 24.10.2004
*/
void axiquadc::res_upper_cond_coup_matrix (long eid,matrix &vm)
{
  long i,j,*ccn;
  matrix lvm(mndofe,nnetp);
  
  ccn = new long [nnetp];
  
  fillm(0.0,vm);

  for (i=0;i<mnb;i++){
    for (j=0;j<ntm;j++){

      //  computation of submatrices
      upper_cond_coup_matrix (eid,i,j,lvm);
      
      if (Cp->bb==lin_lin || Cp->bb==quad_lin){
	Lqat->codnum (ccn,j);
      }
      if (Cp->bb==quad_quad){
	Qqat->codnum (ccn,j);
      }
      
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
void axiquadc::res_lower_cond_coup_matrix (long eid,matrix &vm)
{
  long i,j,*rcn;
  matrix lvm(nnetp,mndofe);
  
  rcn = new long [nnetp];

  fillm(0.0,vm);
  
  for (i=0;i<ntm;i++){
    for (j=0;j<mnb;j++){
      
      //  computation of submatrices
      lower_cond_coup_matrix (eid,i,j,lvm);
      
      if (Cp->bb==lin_lin || Cp->bb==quad_lin){
	Lqat->codnum (rcn,j);
      }
      if (Cp->bb==quad_quad){
	Qqat->codnum (rcn,j);
      }
      
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
void axiquadc::res_upper_cap_coup_matrix (long eid,matrix &vm)
{
  long i,j,*ccn;
  matrix lvm(mndofe,nnetp);
  
  ccn = new long [nnetp];
  
  fillm(0.0,vm);
  
  for (i=0;i<mnb;i++){
    for (j=0;j<ntm;j++){
      
      //  computation of submatrices
      upper_cap_coup_matrix (eid,i,j,lvm);
      
      if (Cp->bb==lin_lin || Cp->bb==quad_lin){
	Lqat->codnum (ccn,j);
      }
      if (Cp->bb==quad_quad){
	Qqat->codnum (ccn,j);
      }
      
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
void axiquadc::res_lower_cap_coup_matrix (long eid,matrix &vm)
{
  long i,j,*rcn;
  matrix lvm(nnetp,mndofe);
  
  rcn = new long [nnetp];

  fillm(0.0,vm);

  for (i=0;i<ntm;i++){
    for (j=0;j<mnb;j++){

      //  computation of submatrices
      lower_cap_coup_matrix (eid,i,j,lvm);
      
      if (Cp->bb==lin_lin || Cp->bb==quad_lin){
	Lqat->codnum (rcn,j);
      }
      if (Cp->bb==quad_quad){
	Qqat->codnum (rcn,j);
      }
      
      mat_localize (vm,lvm,rcn,mordering);
    }
  }
  delete [] rcn;
}





/**
   function computes coupling %vector of axisymmetric problems
   
   @param vm - coupling %vector
   @param nodval - %vector of nodal values
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   
   JK, 19.7.2005
*/
void axiquadc::upper_cond_coup_vector (vector &tvm,vector &nodval,long eid,long ri,long ci)
{
  long i,j,ipp;
  double xi,eta,jac,r;
  ivector nodes(nnetp);
  vector x(nnemp),y(nnemp),w(intordvum[ri][ci]),gp(intordvum[ri][ci]),av(mndofe),t(nnetp);
  matrix gm(mncomp[ri],mndofe),d(tnmcomp,1),dd(mncomp[ri],1),n(1,dofe[ci]),ucm(mndofe,dofe[ci]);
  
  //  nodes on elements
  Tt->give_elemnodes (eid,nodes);
  //  thicknesses
  Tc->give_thickness (eid,nodes,t);
  //  node coordinates
  Mt->give_node_coord2d (x,y,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordvum[ri][ci]);
  
  fillm (0.0,ucm);
  
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
      
      if (Cp->bb==lin_lin){
	Asymlq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	//  radius
	r = Asymlq->approx (xi,eta,x);
	Lqat->bf_matrix (n,xi,eta);
      }
      if (Cp->bb==quad_lin){
	Asymqq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	//  radius
	r = Asymqq->approx (xi,eta,x);
	Lqat->bf_matrix (n,xi,eta);
      }
      if (Cp->bb==quad_quad){
	Asymqq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	//  radius
	r = Asymqq->approx (xi,eta,x);
	Qqat->bf_matrix (n,xi,eta);
      }
      
      
      //  matrix of conductivity of the material
      Cmu->volume_rhs2 (d,ipp,ri,ci);
      
      //  extraction of appropriate block
      dmatblockcol (ri,ci,d,dd);
      
      jac*=r*w[i]*w[j];
      
      //  contribution to the conductivity matrix of the element
      bdbjac (ucm,gm,dd,n,jac);
      
      ipp++;
    }
  }
  
  mxv (ucm,nodval,av);
  addv (av,tvm,tvm);
  
}




/**
   function computes
   
   @param f - 
   @param eid - element id

   JK, 19.7.2005
*/
void axiquadc::res_upper_cond_coup_vector (vector &f,long eid)
{
  long i,j;
  long *lcn;
  vector r(ASTCKVEC(tndofe));
  vector lr;

  //  all initial transport values
  initialvalues (eid, r);
  
  fillv (0.0,f);
  
  for (i=0;i<ntm;i++){
    //  code numbers for unknowns decribing one medium
    lcn = new long [dofe[i]];
    //  initial values of one medium
    allocv (dofe[i],lr);
    
    //  ordering on element
    if (Cp->bb==lin_lin || Cp->bb==quad_lin)
      Lqat->codnum (lcn,i);
    else
      Qqat->codnum (lcn,i);
    
    //  extraction of initial values of one medium
    globloc (r.a,lr.a,lcn,dofe[i]);
    
    for (j=0;j<mnb;j++){
      //  contribution to vector
      upper_cond_coup_vector (f,lr,eid,j,i);
    }
    delete [] lcn;
    destrv (lr);
  }
}














/**
   function copies mechanical data from MEFEL into integration points
   of coupled part in METR
   
   @param eid - element id
   
   JK, 29.10.2004
*/
/* void axiquadc::mefel_metr (long eid)
   {
   long i,ippm,ippcu,ippcl;
   
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
   
   for (i=0;i<tnmcomp;i++){
   Cmu->ip[ippcu].strains[i]  = Mm->ip[ippm].strain[i];
   Cmu->ip[ippcu].stresses[i] = Mm->ip[ippm].stress[i];
   }
   ippcu++;  ippm+=2;
   
   for (i=0;i<tnmcomp;i++){
   Cmu->ip[ippcu].strains[i]  = Mm->ip[ippm].strain[i];
   Cmu->ip[ippcu].stresses[i] = Mm->ip[ippm].stress[i];
   }
   ippcu++;  ippm+=4;
   
   for (i=0;i<tnmcomp;i++){
   Cmu->ip[ippcu].strains[i]  = Mm->ip[ippm].strain[i];
   Cmu->ip[ippcu].stresses[i] = Mm->ip[ippm].stress[i];
   }
   ippcu++;  ippm+=2;
   
   for (i=0;i<tnmcomp;i++){
   Cmu->ip[ippcu].strains[i]  = Mm->ip[ippm].strain[i];
   Cmu->ip[ippcu].stresses[i] = Mm->ip[ippm].stress[i];
   }
   
   //  number of the first integration point on mechanical element
   //  one block formulation of mechanical problem is assumed
   ippm = Mt->elements[eid].ipp[0][0];
   
   for (i=0;i<tnmcomp;i++){
   Cml->ip[ippcl].strains[i]  = Mm->ip[ippm].strain[i];
   Cml->ip[ippcl].stresses[i] = Mm->ip[ippm].stress[i];
   }
   ippcl++;  ippm+=2;
   
   for (i=0;i<tnmcomp;i++){
   Cmu->ip[ippcu].strains[i]  = Mm->ip[ippm].strain[i];
   Cmu->ip[ippcu].stresses[i] = Mm->ip[ippm].stress[i];
   }
   ippcu++;  ippm+=4;
   
   for (i=0;i<tnmcomp;i++){
   Cmu->ip[ippcu].strains[i]  = Mm->ip[ippm].strain[i];
   Cmu->ip[ippcu].stresses[i] = Mm->ip[ippm].stress[i];
   }
   ippcu++;  ippm+=2;
   
   for (i=0;i<tnmcomp;i++){
   Cmu->ip[ippcu].strains[i]  = Mm->ip[ippm].strain[i];
   Cmu->ip[ippcu].stresses[i] = Mm->ip[ippm].stress[i];
   }
   
   }
  
  
   //  each coupling block has its own set of integration points
   //  it is the most general approach but not implemented at this moment
   //  but the code is ready for this extension
   
   if (Cp->savemode==0){
   fprintf (stderr,"\n\n not implemented at this moment (file %s, line %d).\n",__FILE__,__LINE__);
   }
   
   }
*/

/**
   function copies transport data from TRFEL into integration points
   of coupled part in METR
   
   @param eid - element id
   
   JK, 29.10.2004
*/
/* void axiquadc::trfel_metr (long eid)
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
*/

/**
   function approximates temperature into mechanical integration points
   
   @param eid - element id
   
   JK, 31.10.2004
*/
/* void axiquadc::trfel_mefel (long eid)
   {
   long i,j,k,ipp,into;
   double xi,eta,val;
   ivector cn(tndofe);
   vector r(tndofe),t(nnetp),gp,w;
   
   //  code numbers of transport part
   Tt->give_code_numbers (eid,cn.a);
   //  nodal values of transport part
   elemvalues (0,eid,r.a,cn.a,tndofe);
   
   //  index of temperature
   k=Tp->ntm-1;
   
   //  selection of nodal temperatures from all nodal values
   for (i=0;i<dofe[k];i++){
   t[i]=r[Lqat->ordering[k][i]-1];
   }
   
   //  order of numerical integration of mechanical part
   into = Asymqq->intordsm[0][0];
   
   allocv (into,gp);
   allocv (into,w);
   gauss_points (gp.a,w.a,into);
   
   //  number of the first integration point on element
   ipp=Mt->elements[eid].ipp[0][0];
   
   for (i=0;i<into;i++){
   xi=gp[i];
   for (j=0;j<into;j++){
   eta=gp[j];
   
   val = Lqat->approx (xi,eta,t);
   Mm->tempr[ipp]=val;
   ipp++;
   }
   }
   
   destrv (gp);  destrv (w);
   
   }
*/
