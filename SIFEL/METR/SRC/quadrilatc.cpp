#include "quadrilatc.h"
#include "quadlineart.h"
#include "plelemlq.h"
#include "plelemqq.h"
#include "globmat.h"
#include "global.h"
#include "mechtop.h"
#include "globalt.h"
#include "globalc.h"
#include "element.h"
#include "alias.h"
#include "intpoints.h"
#include "globmatt.h"


quadrilatc::quadrilatc (void)
{
  long i,j;

  // ???!!! debug
  // Actually, ssst cannot be simply determined because the state is stored on the general mechanical element class (element.cpp)
  // but not in the particular element type class (e.g. plelemqq.cpp).
  // At this moment, ssst is used to determine number of stress/strain components at integration points
  // which is the same both for planestrain as well as planestress problems, i.e. 4, so ssst is set to plainstrain
  // It must be solved later
  ssst = planestrain;
  if (Cp->bb==lin_lin){
    //  number of blocks of the mechanical element
    mnb=Pelq->nb;
    //  number of mechanical DOFs
    mndofe=Pelq->ndofe;
    //  number of strain/stress components
    tnmcomp=Pelq->tncomp;
    //  array of numbers of components in strain vectors
    mncomp=Pelq->ncomp;
    //  number of nodes in mechanical problems
    nnemp = Pelq->nne;
    //ssst = Pelq->ssst;
    //  number of nodes in transport problems
    nnetp = Lqt->nne;
    //  number of transport DOFs
    tndofe = Lqt->ndofe;
  }

  if (Cp->bb==quad_lin){
    //  number of blocks of the mechanical element
    mnb=Peqq->nb;
    //  number of mechanical DOFs
    mndofe=Peqq->ndofe;
    //  number of strain/stress components
    tnmcomp=Peqq->tncomp;
    //  array of numbers of components in strain vectors
    mncomp=Peqq->ncomp;
    //  number of nodes in mechanical problems
    nnemp = Peqq->nne;
    //ssst = Peqq->ssst;
    //  number of nodes in transport problems
    nnetp = Lqt->nne;
    //  number of transport DOFs
    tndofe = Lqt->ndofe;
  }

  if (Cp->bb==quad_quad){
    //  number of blocks of the mechanical element
    mnb=Peqq->nb;
    //  number of mechanical DOFs
    mndofe=Peqq->ndofe;
    //  number of strain/stress components
    tnmcomp=Peqq->tncomp;
    //  array of numbers of components in strain vectors
    mncomp=Peqq->ncomp;
    //  number of nodes in mechanical problems
    nnemp = Peqq->nne;
    //ssst = Peqq->ssst;
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
  
  //  ordering of unknowns from transport part
  tordering = new long* [ntm];
  for (i=0;i<ntm;i++){
    tordering[i] = new long [nnetp];
  }
  
  
  if (Cp->bb==lin_lin){
    switch (Tp->tmatt){
    case onemedium:{
      if (mnb==1){
	intordvum[0][0]=2;  intordvlm[0][0]=2;
      }
      if (mnb==2){
	intordvum[0][0]=2;  intordvlm[0][0]=2;
	intordvum[1][0]=2;  intordvlm[0][1]=2;
      }

      tordering[0][0]=1;  tordering[0][1]=2; tordering[0][2]=3;  tordering[0][3]=4;

      if (Cp->savemode==0){
	if (mnb==1){
	  nipu[0][0]=4;  nipl[0][0]=4;
	}
	if (mnb==2){
	  nipu[0][0]=4;  nipl[0][0]=4;
	  nipu[1][0]=4;  nipl[0][1]=4;
	}
      }
      if (Cp->savemode==1){
	if (mnb==1){
	  nipu[0][0]=4;  nipl[0][0]=4;
	}
	if (mnb==2){
	  nipu[0][0]=4;  nipl[0][0]=4;
	  nipu[1][0]=0;  nipl[0][1]=0;
	}
      }

      dofe[0]=nnetp;  
      break;
    }
    case twomediacoup:{
      if (mnb==1){
	intordvum[0][0]=2;  intordvum[0][1]=2;
	intordvlm[0][0]=2;  intordvlm[1][0]=2;
      }
      if (mnb==2){
	intordvum[0][0]=2;  intordvum[0][1]=2;  intordvum[1][0]=2;   intordvum[1][1]=2;
	intordvlm[0][0]=2;  intordvlm[0][1]=2;  intordvlm[1][0]=2;   intordvlm[1][1]=2;
      }
      
      tordering[0][0]=1;  tordering[0][1]=3; tordering[0][2]=5;  tordering[0][3]=7;
      tordering[1][0]=2;  tordering[1][1]=4; tordering[1][2]=6;  tordering[1][3]=8;
      
      if (Cp->savemode==0){
	if (mnb==1){
	  nipu[0][0]=4;       nipu[0][1]=4;
	  nipl[0][0]=4;       nipl[1][0]=4;
	}
	if (mnb==2){
	  nipu[0][0]=4;       nipu[0][1]=4; nipu[1][0]=4;       nipu[1][1]=4;
	  nipl[0][0]=4;       nipl[0][1]=4; nipl[1][0]=4;       nipl[1][1]=4;
	}
      }
      if (Cp->savemode==1){
	if (mnb==1){
	  nipu[0][0]=4;       nipu[0][1]=0;
	  nipl[0][0]=4;       nipl[1][0]=0;
	}
	if (mnb==2){
	  nipu[0][0]=4;       nipu[0][1]=0; nipu[1][0]=0;       nipu[1][1]=0;
	  nipl[0][0]=4;       nipl[0][1]=0; nipl[1][0]=0;       nipl[1][1]=0;
	}
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;
      break;
    }
    case threemediacoup:{
      intordvum[0][0]=2;  intordvum[0][1]=2;  intordvum[0][2]=2; intordvum[1][0]=2;  intordvum[1][1]=2;  intordvum[1][2]=2;
      intordvlm[0][0]=2;  intordvlm[1][0]=2;  intordvlm[2][0]=2; intordvlm[0][1]=2;  intordvlm[1][1]=2;  intordvlm[2][1]=2; 
      
      tordering[0][0]=1;  tordering[0][1]=4; tordering[0][2]=7;  tordering[0][3]=10;
      tordering[1][0]=2;  tordering[1][1]=5; tordering[1][2]=8;  tordering[1][3]=11;
      tordering[2][0]=3;  tordering[2][1]=6; tordering[2][2]=9;  tordering[2][3]=12;

      if (Cp->savemode==0){
	nipu[0][0]=4;  nipu[0][1]=4;  nipu[0][2]=4; nipu[1][0]=4;  nipu[1][1]=4;  nipu[1][2]=4;
	nipl[0][0]=4;  nipl[1][0]=4;  nipl[2][0]=4; nipl[0][1]=4;  nipl[1][1]=4;  nipl[2][1]=4;
      }
      if (Cp->savemode==1){
	nipu[0][0]=4;  nipu[0][1]=0;  nipu[0][2]=0; nipu[1][0]=0;  nipu[1][1]=0;  nipu[1][2]=0;
	nipl[0][0]=4;  nipl[1][0]=0;  nipl[2][0]=0; nipl[0][1]=0;  nipl[1][1]=0;  nipl[2][1]=0;
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
      fprintf (stderr,"\n in function quadrilatc::quadrilatc (file %s, line %d).\n",__FILE__,__LINE__);
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
      fprintf (stderr,"\n in function quadrilatc::quadrilatc (file %s, line %d).\n",__FILE__,__LINE__);
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
      tnipl+=nipl[j][i];
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

/*
void quadrilatc::eleminit (long eid)
{
  long ii,jj;

  Ct->elements[eid].nb=mnb;
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
void quadrilatc::dmatblockcol (long ri,long ci,matrix &d, matrix &dd)
{
  fillm (0.0,dd);

  if (mnb==1){
    dd[0][0]=d[0][0];
    dd[1][0]=d[1][0];
    dd[2][0]=d[2][0];
  }
  if (mnb==2){
    if (ri==0 && ci==0){
      dd[0][0]=d[0][0];
      dd[1][0]=d[1][0];
    }
    if (ri==1 && ci==0){
      dd[0][0]=d[2][0];
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
void quadrilatc::dmatblockrow (long ri,long ci,matrix &d, matrix &dd)
{
  fillm (0.0,dd);
  
  if (mnb==1){
    dd[0][0]=d[0][0];
    dd[0][1]=d[0][1];
    dd[0][2]=d[0][2];
  }
  if (mnb==2){
    if (ri==0 && ci==0){
      dd[0][0]=d[0][0];
      dd[0][1]=d[0][1];
    }
    if (ri==0 && ci==1){
      dd[0][0]=d[0][2];
    }
  }
}



/**
   function computes values at integration points from nodal values
   
   @param eid - element id

   TKr, 01/10/2012
   
*/

void quadrilatc::intpointval (long eid)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,val;
  vector r(tndofe),t(nnetp),l(nnetp),gp,w;
  
  elemvalues (eid,r);
  
  
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
	      val = Lqt->approx (xi,eta,t);
	    }
	    if (Cp->bb==quad_lin){
	      val = Lqt->approx (xi,eta,t);
	    }
	    if (Cp->bb==quad_quad){
	      val = Qqt->approx (xi,eta,t);
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
	      val = Lqt->approx (xi,eta,l);
	    }
	    if (Cp->bb==quad_lin){
	      val = Lqt->approx (xi,eta,l);
	    }
	    if (Cp->bb==quad_quad){
	      val = Qqt->approx (xi,eta,l);
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
   
   TKr, 01/10/2012

*/
void quadrilatc::intpointgrad (long eid)
{

  long i,j,ii,jj,k,ipp;
  double xi,eta,jac;
  vector x(nnetp),y(nnetp),r(tndofe),t(nnetp),l(nnetp),gp,w,grad(2);
  matrix gm(2,nnetp);
  
  //Ct->give_node_coord2d (x,y,eid);
  Tt->give_node_coord2d (x,y,eid);

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
	    
	    //  matrix of gradients
	    if (Cp->bb==lin_lin){
	      Lqt->grad_matrix (gm,x,y,xi,eta,jac);
	    }
	    if (Cp->bb==quad_lin){
	      Lqt->grad_matrix (gm,x,y,xi,eta,jac);
	    }
	    if (Cp->bb==quad_quad){
	      Qqt->grad_matrix (gm,x,y,xi,eta,jac);
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
	      Lqt->grad_matrix (gm,x,y,xi,eta,jac);
	    }
	    if (Cp->bb==quad_lin){
	      Lqt->grad_matrix (gm,x,y,xi,eta,jac);
	    }
	    if (Cp->bb==quad_quad){
	      Qqt->grad_matrix (gm,x,y,xi,eta,jac);
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
void quadrilatc::mainip_strains (long eid,long ri,long ci,vector &x,vector &y,vector &r)
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
	  Pelq->geom_matrix_block (gm,ii,x,y,xi,eta,jac);
	}
	if (Cp->bb==quad_lin){
	  Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	}
	if (Cp->bb==quad_quad){
	  Peqq->geom_matrix (gm,x,y,xi,eta,jac);
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
	  Pelq->geom_matrix_block (gm,ii,x,y,xi,eta,jac);
	}
	if (Cp->bb==quad_lin){
	  Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	}
	if (Cp->bb==quad_quad){
	  Peqq->geom_matrix (gm,x,y,xi,eta,jac);
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
void quadrilatc::res_mainip_strains (long eid)
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
      Pelq->transf_matrix (nodes,tmat);
    if (Cp->bb==quad_lin)
      Peqq->transf_matrix (nodes,tmat);
    if (Cp->bb==quad_quad)
      Peqq->transf_matrix (nodes,tmat);

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
  matrix gm(mncomp[ri],mndofe),d(tnmcomp,1),dd(mncomp[ri],1),n(1,dofe[ci]);

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
      
      if (Cp->bb==lin_lin){
	Pelq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->bb==quad_lin){
	//Peqq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->bb==quad_quad){
	//Peqq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Qqt->bf_matrix (n,xi,eta);
      }
      
      //  matrix of conductivity of the material
      Cmu->matcond (d,ii,ri,ci);
      
      dmatblockcol (ri,ci,d,dd);

      if (Cp->bb==lin_lin || Cp->bb==quad_lin)
	thick = Lqt->approx (xi,eta,t);
      else
	thick = Qqt->approx (xi,eta,t);
      
      jac*=thick*w[i]*w[j];
      
      //  contribution to the conductivity matrix of the element
      bdbjac (vm,gm,dd,n,jac);
      
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
  matrix gm(mncomp[ci],mndofe),d(1,tnmcomp),dd(1,mncomp[ci]),n(1,dofe[ri]);

  Tt->give_elemnodes (eid,nodes);
  Tc->give_thickness (eid,nodes,t);

  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordvlm[ri][ci]);

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
      
      if (Cp->bb==lin_lin){
	Pelq->geom_matrix_block (gm,ci,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->bb==quad_lin){
	//Peqq->geom_matrix_block (gm,ci,x,y,xi,eta,jac);
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->bb==quad_quad){
	//Peqq->geom_matrix_block (gm,ci,x,y,xi,eta,jac);
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Qqt->bf_matrix (n,xi,eta);
      }

      //  matrix of conductivity of the material
      Cml->matcond (d,ii,ri,ci);

      dmatblockrow (ri,ci,d,dd);

      if (Cp->bb==lin_lin || Cp->bb==quad_lin)
	thick = Lqt->approx (xi,eta,t);
      else
	thick = Qqt->approx (xi,eta,t);
      
      jac*=thick*ww1*ww2;
      
      //  contribution to the conductivity matrix of the element
      bdbjac (vm,n,dd,gm,jac);

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
  matrix gm(mncomp[ri],mndofe),d(tnmcomp,1),dd(mncomp[ri],1),n(1,dofe[ci]);

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
      
      if (Cp->bb==lin_lin){
	Pelq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->bb==quad_lin){
	//Peqq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->bb==quad_quad){
	//Peqq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Qqt->bf_matrix (n,xi,eta);
      }

      
      //  material capacity matrix
      Cmu->matcap (d,ii,ri,ci);
      
      dmatblockcol (ri,ci,d,dd);

      if (Cp->bb==lin_lin || Cp->bb==quad_lin)
	thick = Lqt->approx (xi,eta,t);
      else
	thick = Qqt->approx (xi,eta,t);

      jac*=thick*ww1*ww2;

      //  contribution to the conductivity matrix of the element
      bdbjac (vm,gm,dd,n,jac);

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
  matrix gm(mncomp[ci],mndofe),d(1,tnmcomp),dd(1,mncomp[ci]),n(1,dofe[ri]);

  Tt->give_elemnodes (eid,nodes);
  Tc->give_thickness (eid,nodes,t);

  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordvlm[ri][ci]);

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
      
      if (Cp->bb==lin_lin){
	Pelq->geom_matrix_block (gm,ci,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->bb==quad_lin){
	//Peqq->geom_matrix_block (gm,ci,x,y,xi,eta,jac);
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->bb==quad_quad){
	//Peqq->geom_matrix_block (gm,ci,x,y,xi,eta,jac);
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Qqt->bf_matrix (n,xi,eta);
      }
      
      //  matrix of conductivity of the material
      Cml->matcap (d,ii,ri,ci);

      dmatblockrow (ri,ci,d,dd);
      
      if (Cp->bb==lin_lin || Cp->bb==quad_lin)
	thick = Lqt->approx (xi,eta,t);
      else
	thick = Qqt->approx (xi,eta,t);
      
      jac*=thick*ww1*ww2;
      
      //  contribution to the conductivity matrix of the element
      bdbjac (vm,n,dd,gm,jac);

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
  matrix lvm;
  
  fillm(0.0,vm);
  
  for (i=0;i<ntm;i++){
    
    ccn = new long [dofe[i]];
    allocm (mndofe,dofe[i],lvm);
    
    //  code numbers of required components
    if (Cp->bb==lin_lin || Cp->bb==quad_lin)
      Lqt->codnum (ccn,i);
    if (Cp->bb==quad_quad)
      Qqt->codnum (ccn,i);
    
    for (j=0;j<mnb;j++){
      //  computation of submatrices
      upper_cond_coup_matrix (eid,j,i,lvm);
      //  localization of the block into the matrix
      mat_localize (vm,lvm,mordering,ccn);
    }

    delete [] ccn;
    destrm (lvm);
  }
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
  matrix lvm;
  
  fillm(0.0,vm);
  
  for (i=0;i<ntm;i++){
    rcn = new long [nnetp];
    allocm (dofe[i],mndofe,lvm);
    
    if (Cp->bb==lin_lin || Cp->bb==quad_lin)
      Lqt->codnum (rcn,i);
    else
      Qqt->codnum (rcn,i);
    
    for (j=0;j<mnb;j++){
      //  computation of submatrices
      lower_cond_coup_matrix (eid,i,j,lvm);
      //  localization of the block into the matrix      
      mat_localize (vm,lvm,rcn,mordering);
    }
    delete [] rcn;
    destrm (lvm);
  }
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
  matrix lvm;
  
  
  fillm(0.0,vm);
  
  for (i=0;i<ntm;i++){
    ccn = new long [nnetp];
    allocm (mndofe,dofe[i],lvm);
    
    if (Cp->bb==lin_lin || Cp->bb==quad_lin)
      Lqt->codnum (ccn,i);
    else
      Qqt->codnum (ccn,i);
    
    for (j=0;j<mnb;j++){
      //  computation of submatrices
      upper_cap_coup_matrix (eid,j,i,lvm);
      //  localization of the block into the matrix
      mat_localize (vm,lvm,mordering,ccn);
    }
    delete [] ccn;
    destrm (lvm);
  }
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
  matrix lvm;
  

  fillm(0.0,vm);

  for (i=0;i<ntm;i++){
    rcn = new long [nnetp];
    allocm (dofe[i],mndofe,lvm);
    
    if (Cp->bb==lin_lin || Cp->bb==quad_lin)
      Lqt->codnum (rcn,i);
    else
      Qqt->codnum (rcn,i);
    
    for (j=0;j<mnb;j++){
      //  computation of submatrices
      lower_cap_coup_matrix (eid,i,j,lvm);
      //  localization of the block into the matrix
      mat_localize (vm,lvm,rcn,mordering);
    }
    delete [] rcn;
    destrm (lvm);
  }
}


/**
   function computes coupling %vector of 2D problems
   
   @param vm - coupling %vector
   @param nodval - %vector of nodal values
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   
   JK, 19.7.2005
*/
void quadrilatc::upper_cond_coup_vector (vector &tvm,vector &nodval,long eid,long ri,long ci)
{
  long i,j,ipp;
  double xi,eta,jac,thick;
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
	Pelq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->bb==quad_lin){
	//Peqq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Lqt->bf_matrix (n,xi,eta);
      }
      if (Cp->bb==quad_quad){
	//Peqq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
	Qqt->bf_matrix (n,xi,eta);
      }
      
      
      //  matrix of conductivity of the material
      Cmu->volume_rhs2 (d,ipp,ri,ci);
      
      //  extraction of appropriate block
      dmatblockcol (ri,ci,d,dd);
      
      if (Cp->bb==lin_lin || Cp->bb==quad_lin)
	thick = Lqt->approx (xi,eta,t);
      else
	thick = Qqt->approx (xi,eta,t);
      
      jac*=thick*w[i]*w[j];
       
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
void quadrilatc::res_upper_cond_coup_vector (vector &f,long eid)
{
  long i,j;
  long *lcn;
  vector r(ASTCKVEC(tndofe));
  vector lr;
  
  //  all initial transport values
  initialvalues (eid, r);
  
  nullv (f);
  
  for (i=0;i<ntm;i++){
    //  code numbers for unknowns decribing one medium
    lcn = new long [dofe[i]];
    //  initial values of one medium
    reallocv (dofe[i],lr);
    
    //  ordering on element
    if (Cp->bb==lin_lin || Cp->bb==quad_lin)
      Lqt->codnum (lcn,i);
    else
      Qqt->codnum (lcn,i);
    
    //  extraction of initial values of one medium
    globloc (r.a,lr.a,lcn,dofe[i]);
    
    for (j=0;j<mnb;j++){
      //  contribution to vector
      upper_cond_coup_vector (f,lr,eid,j,i);
    }
    delete [] lcn;
  }
}






/**
   function computes internal forces caused by particular medium
   
   @param ri,ci - row and column indices
   @param eid - element id
   @param ifo - %vector of internal forces
   
   TKr 23/09/2021 according to JK
*/
void quadrilatc::upper_internal_forces (long ri,long ci,long eid,vector &ifo)
{
  long i,j,ipp;
  double thick,eta,xi,ww1,ww2,jac;
  ivector nodes(nnetp);
  vector x(nnemp),y(nnemp),w(intordvum[ri][ci]),gp(intordvum[ri][ci]),fl(tnmcomp),contr(mndofe),gx(nnemp),gy(nnemp),s(2),t(nnetp);
  matrix gm(tnmcomp,mndofe),d(tnmcomp,1);
  

  Tt->give_elemnodes (eid,nodes);
  Tc->give_thickness (eid,nodes,t);
  Mt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordvum[ri][ci]);
  
  fillv (0.0,ifo);
  
  //  first integration point on element
  if (Cp->savemode==0){
    ipp=Ct->elements[eid].ippu[ri][ci];
  }
  if (Cp->savemode==1){
    ipp=Ct->elements[eid].ippu[0][0];
  }
  
  
  for (i=0;i<intordvum[ri][ci];i++){
    xi=gp[i];  ww1=w[i];
    for (j=0;j<intordvum[ri][ci];j++){
      eta=gp[j];  ww2=w[j];
      
      //  real stress evaluation
      Cmu->computenlstresses (d,ri,ci,ipp);
      
      Cmu->givestresses_cmu (ipp,0,fl);
      
      //  geometric matrix of the element
      if (Cp->bb==lin_lin){
	Pelq->geom_matrix_block (gm,ri,x,y,xi,eta,jac);
	
      }
      if (Cp->bb==quad_lin){
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
      }
      if (Cp->bb==quad_quad){
	Peqq->geom_matrix (gm,x,y,xi,eta,jac);
      }
      
      //  contribution to the result
      mtxv (gm,fl,contr);
      
      
      if (Cp->bb==lin_lin || Cp->bb==quad_lin)
	thick = Lqt->approx (xi,eta,t);
      else
	thick = Qqt->approx (xi,eta,t);
      
      jac*=thick*ww1*ww2;
      
      cmulv (jac,contr);
      
      addv (contr,ifo,ifo);
      
      ipp++;
    }
  }
}

    
/**
   function computes internal forces caused by transport effects
   
   @param eid - element id
   @param ifor - %vector of internal forces
   
   TKr 23/09/2021 according to JK
*/
void quadrilatc::res_upper_internal_forces (long eid,vector &ifor)
{
  long i;
  ivector nodes (nnemp);
  vector v(mndofe);
  vector lif(mndofe);

  fillv(0.0,ifor);
  
  for (i=0;i<mnb;i++){
    fillv(0.0,lif);
    upper_internal_forces (i,0,eid,lif);
    addv (lif,ifor,ifor);
  }
  
  //  this part must be tested??!!
  //  transformation of nodal forces
  //  (in the case of nodal coordinate systems)
  /* Mt->give_elemnodes (eid,nodes);
     long transf = Mt->locsystems (nodes);
     if (transf>0){
     matrix tmat (mndofe,mndofe);
     transf_matrix (nodes,tmat);
     //globloctransf (ifor,v,tmat);
     glvectortransf (ifor,v,tmat);
     copyv (v,ifor);
     }
  */
}


/**
   function computes internal fluxes
   
   @param ri,ci - row and column indices
   @param eid - element id
   @param ifl - %vector of internal fluxes
   
   JK, 17.7.2005, corrected by TKr 25.4.2007
*/
/* void quadrilatc::lower_internal_fluxes (long ri,long ci,long eid,vector &ifl)
   {
   long i,ipp;
   double jac,area,xi;
   vector x(nnetp),gx(nnetp),gy(nnetp),w(intordvlm[ri][ci]),gp(intordvlm[ri][ci]),fl(mndofe),contr(tnmcomp);
   matrix gm(1,nnetp),d(1,tnmcomp);
   
   Tc->give_area (eid,area);
   Tt->give_node_coord2d (gx,gy,eid);
   gauss_points (gp.a,w.a,intordvlm[ri][ci]);
   giveloccoord (gx,gy,x);
   
   fillv (0.0,ifl);
   
   //  first integration point on element
   if (Cp->savemode==0){
   ipp=Ct->elements[eid].ippl[ri][ci];
   }
   if (Cp->savemode==1){
   ipp=Ct->elements[eid].ippl[0][0];
   }
   
   
   for (i=0;i<intordvlm[ri][ci];i++){
   xi=gp[i];
   
   //  real fluxes evaluation
   Cml->computenlfluxes (d,ci,ipp);
   
   Cml->givefluxes_cml (ci,ipp,fl);
   
   
   //  gradient matrix
   if (Cp->bb==lin_lin){
   Lbt->grad_matrix (gm,x,xi,jac);
   }
   if (Cp->bb==quad_lin){
   Lbt->grad_matrix (gm,x,xi,jac);
   }
   if (Cp->bb==quad_quad){
   Qbt->grad_matrix (gm,x,xi,jac);
   }
   
   
   //  contribution to the result
   mtxv (gm,fl,contr);
   cmulv (area*jac*w[i],contr);
   
   addv (contr,ifl,ifl);
   
   ipp++;
   }
   }
*/


/**
   function computes internal fluxes caused by mechanical effects
   
   @param eid - element id
   @param ifl - %vector of internal fluxes
   
   JK, 17.7.2005
*/
/* void quadrilatc::res_lower_internal_fluxes (long eid,vector &ifl)
   {
   long i,*cn;
   vector lif,tdnv;
   matrix cm;
   
   fillv (0.0,ifl);
   cn = new long [tnmcomp];
   allocv (tnmcomp,lif);
   
   for (i=0;i<ntm;i++){
   Qbt->codnum (cn,i);
   lower_internal_fluxes (i,0,eid,lif);
   locglob (ifl.a,lif.a,cn,tnmcomp);
   } 
   
   delete [] cn;
   destrv (lif);
   }
*/




/**
   function computes contributions to the right-hand side - volume integral
      
   \int_{Omega} B^T D dOmega
   
   @param tmv - transmission %vector of one matter
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)
   
   TKr, 23/09/2021 - must be corrected
*/
/* void quadrilatc::volume_rhs_vector (long lcid,long eid,long ri,long ci,vector &vrhs)
   {
   long i,ii;
   double area,xi,ww,jac;
   vector x(nnemp),w(intordvum[ri][ci]),gp(intordvum[ri][ci]),gx(nnemp),gy(nnemp),s(2);
   matrix gm(tnmcomp,mndofe),d;
   matrix km(mndofe,1);
   
   Mc->give_area (eid,area);
   Mt->give_node_coord1d (x,eid);
   gauss_points (gp.a,w.a,intordvum[ri][ci]);
   
   fillm (0.0,km);
   fillv (0.0,vrhs);
   
   if (Cp->savemode==0)
   ii=Ct->elements[eid].ippu[ri][ci];
   if (Cp->savemode==1)
   ii=Ct->elements[eid].ippu[0][0];
   
   for (i=0;i<intordvum[ri][ci];i++){
   xi=gp[i];  ww=w[i];
   
   //  matrix of gradients
   if (Cp->bb==lin_lin){
   Bar2d->geom_matrix (gm,x,jac);
   }
   if (Cp->bb==quad_lin){      
   //  computation of direction vector
   Barq2d->dirvect (s,gx,gy);
   
   Barq2d->geom_matrix (gm,x,s,xi,jac);
   }
   if (Cp->bb==quad_quad){
   //  computation of direction vector
   Barq2d->dirvect (s,gx,gy);
   
   Barq2d->geom_matrix (gm,x,s,xi,jac);
   }
   
   //  matrix of conductivity of the material
   allocm(tnmcomp,1,d);
   
   Cmu->volume_rhs2 (d,ii,ri,ci,tnmcomp);
   
   jac*=area*ww;
   
   //  contribution to the volume_rhs integral of the element
   nnjac (km,gm,d,jac);
   destrm(d);
   
   ii++;
   }
   
   for (i=0;i<vrhs.n;i++){
   vrhs[i] = km[i][0];
   }
   
   }
*/


/**
   function assembles resulting element volume right-hand side

   @param eid - element id
   @param lcid - load case id
   @param f - resulting volume right-hand side %vector of one element

   TKr, 23/09/20211 - must be corrected
*/
/* void quadrilatc::res_volume_rhs_vector (vector &f,long eid,long lcid)
{
  long i,*cn;
  vector lf;
  
  for (i=0;i<mnb;i++){
  cn = new long [mndofe];
  allocv (mndofe,lf);
  
  if (Cp->bb==lin_lin || Cp->bb==quad_lin)
  qt->codnum (ccn,i);
  else
  Qqt->codnum (ccn,i);
  
  volume_rhs_vector (lcid,eid,i,i,lf);
  locglob (f.a,lf.a,cn,mndofe);
  delete [] cn;
  destrv (lf);
  }
  }
*/


/**
   function copies mechanical data from MEFEL into integration points
   of coupled part in METR
   
   @param eid - element id
   
   JK, 29.10.2004
*/
/* void quadrilatc::mefel_metr (long eid)
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
*/

/**
   function copies transport data from TRFEL into integration points
   of coupled part in METR
   
   @param eid - element id
   
   JK, 29.10.2004
*/
/* void quadrilatc::trfel_metr (long eid)
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
/* void quadrilatc::trfel_mefel (long eid)
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
*/


