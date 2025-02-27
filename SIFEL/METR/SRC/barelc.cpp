#include "barelc.h"
#include "global.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "globalt.h"
#include "globalc.h"
#include "barel2d.h"
#include "barelq2d.h"
#include "linbart.h"
#include "quadbart.h"
#include "genfile.h"
#include "globmatt.h"
#include "globmatc.h"
#include "globmat.h"
#include "node.h"
#include "element.h"
#include "intpoints.h"

barelc::barelc (void)
{
  long i,j;
  
  if (Cp->bb==lin_lin){
    //  number of blocks of the mechanical element
    mnb=Bar2d->nb;
    //  number of mechanical DOFs
    mndofe=Bar2d->ndofe;
    //  number of strain/stress components
    tnmcomp=Bar2d->tncomp;
    //  number of nodes in mechanical problems
    nnemp = Bar2d->nne;
    ssst = Bar2d->ssst;
    //  number of nodes in transport problems
    nnetp = Lbt->nne;
    //  number of transport DOFs
    tndofe = Lbt->ndofe;
  }

  if (Cp->bb==quad_lin){
    //  number of blocks of the mechanical element
    mnb=Barq2d->nb;
    //  number of mechanical DOFs
    mndofe=Barq2d->ndofe;
    //  number of strain/stress components
    tnmcomp=Barq2d->tncomp;
    //  number of nodes in mechanical problems
    nnemp = Barq2d->nne;
    ssst = Barq2d->ssst;
    //  number of nodes in transport problems
    nnetp = Lbt->nne;
    //  number of transport DOFs
    tndofe = Lbt->ndofe;
  }

  if (Cp->bb==quad_quad){
    //  number of blocks of the mechanical element
    mnb=Barq2d->nb;
    //  number of mechanical DOFs
    mndofe=Barq2d->ndofe;
    //  number of strain/stress components
    tnmcomp=Barq2d->tncomp;
    //  number of nodes in mechanical problems
    nnemp = Barq2d->nne;
    ssst = Barq2d->ssst;
    //  number of nodes in transport problems
    nnetp = Qbt->nne;
    //  number of transport DOFs
    tndofe = Qbt->ndofe;
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
      intordvum[0][0]=1;  intordvlm[0][0]=1;
      nipu[0][0]=1;  nipl[0][0]=1;
      dofe[0]=nnetp;
      tordering[0][0]=1;  tordering[0][1]=2;
      break;
    }
    case twomediacoup:{
      tordering[0][0]=1;  tordering[0][1]=3;
      tordering[1][0]=2;  tordering[1][1]=4;

      intordvum[0][0]=1;  intordvum[0][1]=1;
      intordvlm[0][0]=1;  intordvlm[1][0]=1;

      if (Cp->savemode==0){
      nipu[0][0]=1;       nipu[0][1]=1;
      nipl[0][0]=1;       nipl[1][0]=1;
      }
      if (Cp->savemode==1){
      nipu[0][0]=1;       nipu[0][1]=0;
      nipl[0][0]=1;       nipl[1][0]=0;
      }

      dofe[0]=nnetp;  dofe[1]=nnetp;
      break;
    }
    case threemediacoup:{
      tordering[0][0]=1;   tordering[0][1]=4;
      tordering[1][0]=2;   tordering[1][1]=5;
      tordering[2][0]=3;   tordering[2][1]=6;

      intordvum[0][0]=1;  intordvum[0][1]=1;  intordvum[0][2]=1;
      intordvlm[0][0]=1;  intordvlm[1][0]=1;  intordvlm[2][0]=1;
      
      if (Cp->savemode==0){
        nipu[0][0]=1;  nipu[0][1]=1;  nipu[0][2]=1;
        nipl[0][0]=1;  nipl[1][0]=1;  nipl[2][0]=1;
      }
      if (Cp->savemode==1){
        nipu[0][0]=1;  nipu[0][1]=0;  nipu[0][2]=0;
        nipl[0][0]=1;  nipl[1][0]=0;  nipl[2][0]=0;
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;  dofe[2]=nnetp;

      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown type of transported matter is required");
      fprintf (stderr,"\n in function barelc (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }

  if (Cp->bb==quad_lin){
    switch (Tp->tmatt){
    case onemedium:{
      intordvum[0][0]=2;  intordvlm[0][0]=2;
      nipu[0][0]=2;  nipl[0][0]=2;
      dofe[0]=nnetp;
      tordering[0][0]=1;  tordering[0][1]=2;
      break;
    }
    case twomediacoup:{
      tordering[0][0]=1;  tordering[0][1]=3;
      tordering[1][0]=2;  tordering[1][1]=4;

      intordvum[0][0]=2;  intordvum[0][1]=2;
      intordvlm[0][0]=2;  intordvlm[1][0]=2;

      if (Cp->savemode==0){
      nipu[0][0]=2;       nipu[0][1]=2;
      nipl[0][0]=2;       nipl[1][0]=2;
      }
      if (Cp->savemode==1){
      nipu[0][0]=2;       nipu[0][1]=0;
      nipl[0][0]=2;       nipl[1][0]=0;
      }

      dofe[0]=nnetp;  dofe[1]=nnetp;
      break;
    }
    case threemediacoup:{
      tordering[0][0]=1;   tordering[0][1]=4;
      tordering[1][0]=2;   tordering[1][1]=5;
      tordering[2][0]=3;   tordering[2][1]=6;
    
      intordvum[0][0]=2;  intordvum[0][1]=2;  intordvum[0][2]=2;
      intordvlm[0][0]=2;  intordvlm[1][0]=2;  intordvlm[2][0]=2;
      
      if (Cp->savemode==0){
        nipu[0][0]=2;  nipu[0][1]=2;  nipu[0][2]=2;
        nipl[0][0]=2;  nipl[1][0]=2;  nipl[2][0]=2;
      }
      if (Cp->savemode==1){
        nipu[0][0]=2;  nipu[0][1]=0;  nipu[0][2]=0;
        nipl[0][0]=2;  nipl[1][0]=0;  nipl[2][0]=0;
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;  dofe[2]=nnetp;

      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown type of transported matter is required");
      fprintf (stderr,"\n in function barelc (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }

  if (Cp->bb==quad_quad){
    switch (Tp->tmatt){
    case onemedium:{
      intordvum[0][0]=2;  intordvlm[0][0]=2;
      nipu[0][0]=2;  nipl[0][0]=2;
      dofe[0]=nnetp;
      tordering[0][0]=1;  tordering[0][1]=2;  tordering[0][2]=3;
      break;
    }
    case twomediacoup:{
      tordering[0][0]=1;  tordering[0][1]=3;  tordering[0][2]=5;
      tordering[1][0]=2;  tordering[1][1]=4;  tordering[1][2]=6;

      intordvum[0][0]=2;  intordvum[0][1]=2;
      intordvlm[0][0]=2;  intordvlm[1][0]=2;

      if (Cp->savemode==0){
      nipu[0][0]=2;       nipu[0][1]=2;
      nipl[0][0]=2;       nipl[1][0]=2;
      }
      if (Cp->savemode==1){
      nipu[0][0]=2;       nipu[0][1]=0;
      nipl[0][0]=2;       nipl[1][0]=0;
      }

      dofe[0]=nnetp;  dofe[1]=nnetp;
      break;
    }
    case threemediacoup:{
      tordering[0][0]=1;   tordering[0][1]=4;   tordering[0][2]=7;
      tordering[1][0]=2;   tordering[1][1]=5;   tordering[1][2]=8;
      tordering[2][0]=3;   tordering[2][1]=6;   tordering[2][2]=9;

      intordvum[0][0]=2;  intordvum[0][1]=2;  intordvum[0][2]=2;
      intordvlm[0][0]=2;  intordvlm[1][0]=2;  intordvlm[2][0]=2;
      
      if (Cp->savemode==0){
        nipu[0][0]=2;  nipu[0][1]=2;  nipu[0][2]=2;
        nipl[0][0]=2;  nipl[1][0]=2;  nipl[2][0]=2;
      }
      if (Cp->savemode==1){
        nipu[0][0]=2;  nipu[0][1]=0;  nipu[0][2]=0;
        nipl[0][0]=2;  nipl[1][0]=0;  nipl[2][0]=0;
      }
      
      dofe[0]=nnetp;  dofe[1]=nnetp;  dofe[2]=nnetp;

      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown type of transported matter is required");
      fprintf (stderr,"\n in function barelc (file %s, line %d).\n",__FILE__,__LINE__);
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

  /*  for (i=0;i<mnb;i++){
      for (j=0;j<ntm;j++)
      tnipl+=nipl[i][j];
      }
  */
  
  for (i=0;i<ntm;i++){
    for (j=0;j<mnb;j++)
      tnipl+=nipl[i][j];
  }

}

barelc::~barelc (void)
{
  long i;

  for (i=0;i<ntm;i++){
    delete [] nipl[i];
    delete [] intordvlm[i];
    delete [] tordering[i];
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
  delete [] tordering;
}

/*
  function initializes order of numerical integration, number of integration points,
  
  @param eid - element id
  
*/
/*
void barelc::eleminit (long eid)
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

void barelc::giveloccoord(vector &x, vector &y,vector &lx)
{
  double l;
  
  if(x.n == 2){
    lx[0] = 0.0;
    l = sqr(x[1] - x[0]) + sqr(y[1] - y[0]);
    l = sqrt(l);
    lx[1] = l;
  }
  if(x.n == 3){
    lx[0] = 0.0;
    l = sqr(x[1] - x[0]) + sqr(y[1] - y[0]);
    l = sqrt(l);
    lx[1] = l;
    l = sqr(x[2] - x[0]) + sqr(y[2] - y[0]);
    l = sqrt(l);
    lx[2] = l;
  }
}


/**
   transformation %matrix x_g = T x_l
*/
void barelc::transf_matrix (ivector &nodes,matrix &tmat)
{
  long i;
  fillm (0.0,tmat);

  for (i=0;i<tmat.m;i++){
    tmat[i][i]=1.0;
  }

  for (i=0;i<nodes.n;i++){//temp.
    if (Mt->nodes[nodes[i]].transf>0){//temp.
      tmat[i*2][i*2]=Mt->nodes[nodes[i]].e1[0];    tmat[i*2][i*2+1]=Mt->nodes[nodes[i]].e2[0];
      tmat[i*2+1][i*2]=Mt->nodes[nodes[i]].e1[1];  tmat[i*2+1][i*2+1]=Mt->nodes[nodes[i]].e2[1];
    }
  }

}


/**
   function computes values at integration points from nodal values
   
   @param eid - element id

   TKr, 22.3.3004
   
*/

void barelc::intpointval (long eid)
{
  long i,k,ii,jj,ipp;
  double xi,val;

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
          
          if (Cp->bb==lin_lin){
            val = Lbt->approx (xi,t);
          }
          if (Cp->bb==quad_lin){
            val = Lbt->approx (xi,t);
          }
          if (Cp->bb==quad_quad){
            val = Qbt->approx (xi,t);
          }

          Cmu->ip[ipp].av[k]=val;
          ipp++;
        }
        destrv (gp);  destrv (w);

        //lower int. points
        allocv (intordvlm[jj][ii],gp);
        allocv (intordvlm[jj][ii],w);
        gauss_points (gp.a,w.a,intordvlm[jj][ii]);
        
        ipp=Ct->elements[eid].ippl[jj][ii];
        for (i=0;i<intordvlm[jj][ii];i++){
          xi=gp[i];
          
          if (Cp->bb==lin_lin){
            val = Lbt->approx (xi,l);
          }
          if (Cp->bb==quad_lin){
            val = Lbt->approx (xi,l);
          }
          if (Cp->bb==quad_quad){
            val = Qbt->approx (xi,l);
          }

          Cml->ip[ipp].av[k]=val;
          ipp++;
        }
        
        destrv (gp);  destrv (w);

        if (Cp->savemode==1)  break;//added 14/03/2011

      }
      if (Cp->savemode==1)  break;//added 14/03/2011
    }
  }
}



/**
   function computes gradients in integration points from nodal values
   
   @param eid - element id
   
   TKr, 22.3.2004

*/
void barelc::intpointgrad (long eid)
{

  long i,ii,jj,k,ipp;
  double xi,jac;
  vector x(nnetp),y(nnetp),r(tndofe),t(nnetp),l(nnetp),gp,w,grad(1);
  matrix gm(1,nnetp);
  
  //Ct->give_node_coord2d (x,y,eid)
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
          
          //  matrix of gradients
          if (Cp->bb==lin_lin){
            Lbt->grad_matrix (gm,x,xi,jac);
          }
          if (Cp->bb==quad_lin){
            Lbt->grad_matrix (gm,x,xi,jac);
          }
          if (Cp->bb==quad_quad){
            Qbt->grad_matrix (gm,x,xi,jac);
          }

          mxv (gm,t,grad);
          Cmu->storegrad_cmu (k,ipp,grad);
          ipp++;
        }
        destrv (gp);  destrv (w);
        
        //lower int. points
        allocv (intordvlm[jj][ii],gp);
        allocv (intordvlm[jj][ii],w);
        gauss_points (gp.a,w.a,intordvlm[jj][ii]);
        
        ipp=Ct->elements[eid].ippl[jj][ii];
        for (i=0;i<intordvlm[jj][ii];i++){
          xi=gp[i];
                  
          //  matrix of gradients
          if (Cp->bb==lin_lin){
            Lbt->grad_matrix (gm,x,xi,jac);
          }
          if (Cp->bb==quad_lin){
            Lbt->grad_matrix (gm,x,xi,jac);
          }
          if (Cp->bb==quad_quad){
            Qbt->grad_matrix (gm,x,xi,jac);
          }

          mxv (gm,l,grad);
          Cml->storegrad_cml (k,ipp,grad);
          ipp++;
        }
        
        destrv (gp);  destrv (w);

	if (Cp->savemode==1)  break;//added 14/03/2011
      }
      if (Cp->savemode==1)  break;//added 14/03/2011
    } 
  }
}


/**
   function computes strains in integration points of element

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param x - arrays with node coordinates
   @param s - direction %vector
   @param r - %vector of nodal displacements

   TKr, 22.3.2004
*/
void barelc::mainip_strains (long eid,long ri,long ci,vector &x,vector &y,vector &r)
{
  long lri,lci;
  long i,ii,ipp;
  double xi,jac;
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
    reallocv (RSTCKVEC(1,eps));
    reallocm (RSTCKMAT(1,mndofe,gm));
    
    gauss_points (gp.a,w.a,intordvlm[lri][lci]);
    
    //lower int. points
    ipp=Ct->elements[eid].ippl[lri][lci];
    
    for (i=0;i<intordvlm[lri][lci];i++){
      xi=gp[i];
      
      if (Cp->bb==lin_lin){
        Bar2d->geom_matrix (gm,x,y,jac);
      }
      if (Cp->bb==quad_lin){
        Barq2d->geom_matrix (gm,x,y,xi,jac);
      }
      if (Cp->bb==quad_quad){
        Barq2d->geom_matrix (gm,x,y,xi,jac);
      }
      
      mxv (gm,r,eps);
      
      Cml->storestrain_cml (ipp,0,eps);
      ipp++;
    }
    
    //upper int. points - musi byt opacne ri,ci!
    ipp=Ct->elements[eid].ippu[lci][lri];
    for (i=0;i<intordvum[lci][lri];i++){
      xi=gp[i];
      
      if (Cp->bb==lin_lin){
        Bar2d->geom_matrix (gm,x,y,jac);
      }
      if (Cp->bb==quad_lin){
        Barq2d->geom_matrix (gm,x,y,xi,jac);
      }
      if (Cp->bb==quad_quad){
        Barq2d->geom_matrix (gm,x,y,xi,jac);
      }
      
      mxv (gm,r,eps);
      
      Cmu->storestrain_cmu (ipp,0,eps);
      ipp++;
    }
  }
}


/**
   function computes resulting strains in integration points of element

   @param eid - element id

   TKr, 22.3.2004
*/
void barelc::res_mainip_strains (long eid)
{
  long i,transf;
  vector aux,gx(ASTCKVEC(nnemp)),gy(ASTCKVEC(nnemp)),x(ASTCKVEC(nnemp)),r(ASTCKVEC(mndofe)),s(ASTCKVEC(2));
  ivector nodes(ASTCKIVEC(nnemp));
  matrix tmat;
  
  Ct->give_node_coord2d (gx,gy,eid);
  Ct->give_elemnodes (eid,nodes);
  
  eldispl (0,eid,r.a);
  
  //  transformation of displacement vector
  //  (in the case of nodal coordinate systems)
  transf = Mt->locsystems (nodes);
  if (transf>0){
    allocv (mndofe,aux);
    allocm (mndofe,mndofe,tmat);
    transf_matrix (nodes,tmat);
    lgvectortransf (aux,r,tmat);
    copyv (aux,r);
    destrv (aux);
    destrm (tmat);
  }
  
  giveloccoord (gx,gy,x);
  
  if ((Cp->bb==quad_lin) || (Cp->bb==quad_quad))
    Barq2d->dirvect (s,gx,gy);
  
  for (i=0;i<Tp->ntm;i++){
    mainip_strains (eid,i,0,x,s,r);
  }
}

/**
   function computes upper coupling conductivity %matrix of 1D problems
   
   @param lcid - number of load case
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param vm - coupling %matrix
   
   JK, 10.4.2003, corrected by TKr 25.4.2007
*/
void barelc::upper_cond_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,ipp;
  double area,xi,jac;
  vector x(nnemp),w(intordvum[ri][ci]),gp(intordvum[ri][ci]),gx(nnemp),gy(nnemp),s(2);
  matrix gm(tnmcomp,mndofe),d(tnmcomp,1),n(1,nnetp);
  
  Mc->give_area (eid,area);
  Mt->give_node_coord2d (gx,gy,eid);
  gauss_points (gp.a,w.a,intordvum[ri][ci]);
  giveloccoord (gx,gy,x);
  
  fillm (0.0,vm);
  
  if (Cp->savemode==0){
    ipp=Ct->elements[eid].ippu[ri][ci];
  }
  if (Cp->savemode==1){
    ipp=Ct->elements[eid].ippu[0][0];
  }

  for (i=0;i<intordvum[ri][ci];i++){
    xi=gp[i];

    if (Cp->bb==lin_lin){
      Bar2d->geom_matrix (gm,gx,gy,jac);
      Lbt->bf_matrix (n,xi);
    }
    if (Cp->bb==quad_lin){      
      //  computation of direction vector
      //Barq2d->dirvect (s,gx,gy);

      Barq2d->geom_matrix (gm,gx,gy,xi,jac);
      Lbt->bf_matrix (n,xi);
    }
    if (Cp->bb==quad_quad){
      //  computation of direction vector
      Barq2d->dirvect (s,gx,gy);
      
      Barq2d->geom_matrix (gm,x,s,xi,jac);
      Qbt->bf_matrix (n,xi);
    }
    
    //  matrix of conductivity of the material
    Cmu->matcond (d,ipp,ri,ci);
    
    jac*=area*w[i];
    
    //  contribution to the conductivity matrix of the element
    bdbjac (vm,gm,d,n,jac);
    
    ipp++;
  }
}


/**
   function computes lower coupling conductivity %matrix of 1D problems
   
   @param lcid - number of load case
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting matrix
   @param vm - coupling %matrix
   
   JK, 10.4.2003, corrected by TKr 25.4.2007
*/
void barelc::lower_cond_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,ipp;
  double area,xi,jac;
  vector x(nnemp),w(intordvlm[ri][ci]),gp(intordvlm[ri][ci]),gx(nnemp),gy(nnemp),s(2);
  matrix gm(tnmcomp,mndofe),d(1,tnmcomp),n(1,nnetp);

  Mc->give_area (eid,area);
  Mt->give_node_coord2d (gx,gy,eid);
  gauss_points (gp.a,w.a,intordvlm[ri][ci]);
  giveloccoord (gx,gy,x);

  fillm (0.0,vm);

  if (Cp->savemode==0){
    ipp=Ct->elements[eid].ippl[ri][ci];
  }
  if (Cp->savemode==1){
    ipp=Ct->elements[eid].ippl[0][0];
  }
  
  for (i=0;i<intordvlm[ri][ci];i++){
    xi=gp[i];
    
    if (Cp->bb==lin_lin){
      Bar2d->geom_matrix (gm,gx,gy,jac);
      Lbt->bf_matrix (n,xi);
    }
    if (Cp->bb==quad_lin){
      //  computation of direction vector
      //Barq2d->dirvect (s,gx,gy);

      Barq2d->geom_matrix (gm,gx,gy,xi,jac);
      Lbt->bf_matrix (n,xi);
    }
    if (Cp->bb==quad_quad){
      //  computation of direction vector
      Barq2d->dirvect (s,gx,gy);

      Barq2d->geom_matrix (gm,x,s,xi,jac);
      Qbt->bf_matrix (n,xi);
    }
    
    //  matrix of conductivity of the material
    Cml->matcond (d,ipp,ri,ci);
    
    jac*=area*w[i];
    
    //  contribution to the conductivity matrix of the element
    bdbjac (vm,n,d,gm,jac);
    
    ipp++;
  }
  
}

/**
   function computes upper coupling capacity %matrix of 1D problems
   
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting matrix
   @param vm - coupling %matrix
   
   JK, 9.1.2003, corrected by TKr 25.4.2007
*/
void barelc::upper_cap_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,ipp;
  double area,xi,jac;
  vector x(nnemp),w(intordvum[ri][ci]),gp(intordvum[ri][ci]),gx(nnemp),gy(nnemp),s(2);
  matrix gm(tnmcomp,mndofe),d(tnmcomp,1),n(1,nnetp);

  Mc->give_area (eid,area);
  Mt->give_node_coord2d (gx,gy,eid);
  gauss_points (gp.a,w.a,intordvum[ri][ci]);
  giveloccoord (gx,gy,x);

  fillm (0.0,vm);

  if (Cp->savemode==0){
    ipp=Ct->elements[eid].ippu[ri][ci];
  }
  if (Cp->savemode==1){
    ipp=Ct->elements[eid].ippu[0][0];
  }

  for (i=0;i<intordvum[ri][ci];i++){
    xi=gp[i];
    
    if (Cp->bb==lin_lin){
      Bar2d->geom_matrix (gm,gx,gy,jac);
      Lbt->bf_matrix (n,xi);
    }
    if (Cp->bb==quad_lin){
      //  computation of direction vector
      //Barq2d->dirvect (s,gx,gy);
      
      Barq2d->geom_matrix (gm,gx,gy,xi,jac);
      Lbt->bf_matrix (n,xi);
    }
    if (Cp->bb==quad_quad){
      //  computation of direction vector
      Barq2d->dirvect (s,gx,gy);

      Barq2d->geom_matrix (gm,x,s,xi,jac);
      Qbt->bf_matrix (n,xi);
    }
    
    //  matrix of conductivity of the material
    Cmu->matcap (d,ipp,ri,ci);
    
    jac*=area*w[i];
    
    //  contribution to the conductivity matrix of the element
    bdbjac (vm,gm,d,n,jac);
    
    ipp++;
  }
}


/**
   function computes lower coupling capacity %matrix of 1D problems
   
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param vm - coupling %matrix
   
   JK, 10.4.2003, corrected by TKr 25.4.2007
*/
void barelc::lower_cap_coup_matrix (long eid,long ri,long ci,matrix &vm)
{
  long i,ipp;
  double area,xi,jac;
  vector x(nnemp),w(intordvlm[ri][ci]),gp(intordvlm[ri][ci]),gx(nnemp),gy(nnemp),s(2);
  matrix gm(tnmcomp,mndofe),d(1,tnmcomp),n(1,nnetp);
  
  Mc->give_area (eid,area);
  Mt->give_node_coord2d (gx,gy,eid);
  gauss_points (gp.a,w.a,intordvlm[ri][ci]);
  giveloccoord (gx,gy,x);

  fillm (0.0,vm);
  
  if (Cp->savemode==0){
    ipp=Ct->elements[eid].ippl[ri][ci];
  }
  if (Cp->savemode==1){
    ipp=Ct->elements[eid].ippl[0][0];
  }

  for (i=0;i<intordvlm[ri][ci];i++){
    xi=gp[i];
    
    if (Cp->bb==lin_lin){
      Bar2d->geom_matrix (gm,gx,gy,jac);
      Lbt->bf_matrix (n,xi);
    }
    if (Cp->bb==quad_lin){
      //  computation of direction vector
      //Barq2d->dirvect (s,gx,gy);
      
      Barq2d->geom_matrix (gm,gx,gy,xi,jac);
      Lbt->bf_matrix (n,xi);
    }
    if (Cp->bb==quad_quad){
      //  computation of direction vector
      Barq2d->dirvect (s,gx,gy);

      Barq2d->geom_matrix (gm,x,s,xi,jac);
      Qbt->bf_matrix (n,xi);
    }
    
    //  matrix of conductivity of the material
    Cml->matcap (d,ipp,ri,ci);
    
    jac*=area*w[i];
    
    //  contribution to the conductivity matrix of the element
    bdbjac (vm,n,d,gm,jac);
    
    ipp++;
  }
}



/**
   function assembles upper coupling conductivity %matrices into one element %matrix
   
   @param eid - element id
   @param vm - element upper coupling conductivity %matrix
   
   JK, 17.7.2005
*/
void barelc::res_upper_cond_coup_matrix (long eid,matrix &vm)
{
  long i,j,*ccn;
  matrix lvm(mndofe,nnetp);
  
  ccn = new long [nnetp];

  fillm(0.0,vm);
  
  for (i=0;i<mnb;i++){
    for (j=0;j<ntm;j++){
      
      //  computation of submatrices
      upper_cond_coup_matrix (eid,i,j,lvm);
      
      if (Cp->bb==lin_lin || Cp->bb==quad_lin)
        Lbt->codnum (ccn,j);
      else
        Qbt->codnum (ccn,j);
      
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
void barelc::res_lower_cond_coup_matrix (long eid,matrix &vm)
{
  long i,j,*rcn;
  matrix lvm(nnetp,mndofe);
  
  rcn = new long [nnetp];
  
  fillm(0.0,vm);
  
  for (i=0;i<ntm;i++){
    for (j=0;j<mnb;j++){
      
      //  computation of submatrices
      lower_cond_coup_matrix (eid,i,j,lvm);
      
      if (Cp->bb==lin_lin || Cp->bb==quad_lin)
        Lbt->codnum (rcn,j);
      else
        Qbt->codnum (rcn,j);
      
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
void barelc::res_upper_cap_coup_matrix (long eid,matrix &vm)
{
  long i,j,*ccn;
  matrix lvm(mndofe,nnetp);
  
  ccn = new long [nnetp];

  fillm(0.0,vm);

  for (i=0;i<mnb;i++){
    for (j=0;j<ntm;j++){
      
      //  computation of submatrices
      upper_cap_coup_matrix (eid,i,j,lvm);
      
      if (Cp->bb==lin_lin || Cp->bb==quad_lin)
        Lbt->codnum (ccn,j);
      else
        Qbt->codnum (ccn,j);

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
void barelc::res_lower_cap_coup_matrix (long eid,matrix &vm)
{
  long i,j,*rcn;
  matrix lvm(nnetp,mndofe);
  
  rcn = new long [nnetp];

  fillm(0.0,vm);

  for (i=0;i<ntm;i++){
    for (j=0;j<mnb;j++){
      
      //  computation of submatrices
      lower_cap_coup_matrix (eid,i,j,lvm);
      
      if (Cp->bb==lin_lin || Cp->bb==quad_lin)
        Lbt->codnum (rcn,j);
      else
        Qbt->codnum (rcn,j);

      mat_localize (vm,lvm,rcn,mordering);
    }
  }
  delete [] rcn;
}







/**
   function computes coupling %vector of 1D problems
   
   @param vm - coupling %vector
   @param nodval - %vector of nodal values
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   
   JK, 10.4.2003, corrected by TKr 25.4.2007
*/
void barelc::upper_cond_coup_vector (vector &tvm,vector &nodval,long eid,long ri,long ci)
{
  long i,ipp;
  double area,xi,jac;
  vector x(nnemp),w(intordvum[ri][ci]),gp(intordvum[ri][ci]),av(mndofe),gx(nnemp),gy(nnemp),s(2);
  matrix gm(tnmcomp,mndofe),d(tnmcomp,1),n(1,nnetp),ucm(mndofe,nnetp);
  
  Mc->give_area (eid,area);  
  Mt->give_node_coord2d (gx,gy,eid);
  gauss_points (gp.a,w.a,intordvum[ri][ci]);
  giveloccoord (gx,gy,x);

  fillm (0.0,ucm);
  
  if (Cp->savemode==0){
    ipp=Ct->elements[eid].ippu[ri][ci];
  }
  if (Cp->savemode==1){
    ipp=Ct->elements[eid].ippu[0][0];
  }
  
  for (i=0;i<intordvum[ri][ci];i++){
    xi=gp[i];

    if (Cp->bb==lin_lin){
      Bar2d->geom_matrix (gm,gx,gy,jac);
      Lbt->bf_matrix (n,xi);
    }
    if (Cp->bb==quad_lin){      
      //  computation of direction vector
      //Barq2d->dirvect (s,gx,gy);

      Barq2d->geom_matrix (gm,gx,gy,xi,jac);
      Lbt->bf_matrix (n,xi);
    }
    if (Cp->bb==quad_quad){
      //  computation of direction vector
      Barq2d->dirvect (s,gx,gy);

      Barq2d->geom_matrix (gm,x,s,xi,jac);
      Qbt->bf_matrix (n,xi);
    }
    
    //  matrix of conductivity of the material
    Cmu->volume_rhs2 (d,ipp,ri,ci);
    
    jac*=area*w[i];
    
    //  contribution to the conductivity matrix of the element
    bdbjac (ucm,gm,d,n,jac);
    
    ipp++;
  }
  
  mxv (ucm,nodval,av);
  addv (av,tvm,tvm);
}




/**
   function computes
   
   @param f - %vector of right-hand side
   @param eid - element id

   JK, 15.7.2005
*/
void barelc::res_upper_cond_coup_vector (vector &f,long eid)
{
  long i;
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
      Lbt->codnum (lcn,i);
    else
      Qbt->codnum (lcn,i);
    
    //  extraction of initial values of one medium
    globloc (r.a, lr.a, lcn, dofe[i]);
    
    upper_cond_coup_vector (f,lr,eid,0,i);
    
    delete [] lcn;
  }
}






/**
   function computes internal forces caused by particular medium
   
   @param ri,ci - row and column indices
   @param eid - element id
   @param ifo - %vector of internal forces
   
   JK, 17.7.2005, corrected by TKr 25.4.2007
*/
void barelc::upper_internal_forces (long ri,long ci,long eid,vector &ifo)
{
  long i,ipp;
  double jac,area,xi;
  vector x(nnemp),w(intordvum[ri][ci]),gp(intordvum[ri][ci]),fl(tnmcomp),contr(mndofe),gx(nnemp),gy(nnemp),s(2);
  matrix gm(tnmcomp,mndofe),d(tnmcomp,1);
  
  Mc->give_area (eid,area);
  Mt->give_node_coord2d (gx,gy,eid);
  gauss_points (gp.a,w.a,intordvum[ri][ci]);
  giveloccoord (gx,gy,x);
  
  fillv (0.0,ifo);
  
  //  first integration point on element
  if (Cp->savemode==0){
    ipp=Ct->elements[eid].ippu[ri][ci];
  }
  if (Cp->savemode==1){
    ipp=Ct->elements[eid].ippu[0][0];
  }
  
  
  for (i=0;i<intordvum[ri][ci];i++){
    xi=gp[i];
    
    //  real stress evaluation
    Cmu->computenlstresses (d,ri,ci,ipp);
    
    Cmu->givestresses_cmu (ipp,ri,fl);
    
    //  geometric matrix of the element
    if (Cp->bb==lin_lin){
      Bar2d->geom_matrix (gm,gx,gy,jac);
    }
    if (Cp->bb==quad_lin){
      //  computation of direction vector
      //Barq2d->dirvect (s,gx,gy);
      
      Barq2d->geom_matrix (gm,gx,gy,xi,jac);
    }
    if (Cp->bb==quad_quad){
      //  computation of direction vector
      Barq2d->dirvect (s,gx,gy);
      
      Barq2d->geom_matrix (gm,x,s,xi,jac);
    }
    
    //  contribution to the result
    mtxv (gm,fl,contr);
    cmulv (area*jac*w[i],contr);
    
    addv (contr,ifo,ifo);
    
    ipp++;
  }
}

    
/**
   function computes internal forces caused by transport effects
   
   @param eid - element id
   @param ifor - %vector of internal forces
   
   JK, 17.7.2005
*/
void barelc::res_upper_internal_forces (long eid,vector &ifor)
{
  long i;
  vector lif(mndofe);

  fillv(0.0,ifor);
  
  for (i=0;i<mnb;i++){
    fillv(0.0,lif);
    upper_internal_forces (i,0,eid,lif);
    addv (lif,ifor,ifor);
  }
}

/**
   function computes internal fluxes
   
   @param ri,ci - row and column indices
   @param eid - element id
   @param ifl - %vector of internal fluxes
   
   JK, 17.7.2005, corrected by TKr 25.4.2007
*/
void barelc::lower_internal_fluxes (long ri,long ci,long eid,vector &ifl)
{
  long i,ipp;
  double jac,area,xi;
  ivector nodes(nnetp);
  vector x(nnetp),gx(nnetp),gy(nnetp),w(intordvlm[ri][ci]),gp(intordvlm[ri][ci]),fl(mndofe),contr(tnmcomp),a(nnetp);
  matrix gm(1,nnetp),d(1,tnmcomp);
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  //  node coordinates
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

    //  area of cross section
    area = Lbt->approx (xi,a);
    
    //  contribution to the result
    mtxv (gm,fl,contr);
    cmulv (area*jac*w[i],contr);
    
    addv (contr,ifl,ifl);
    
    ipp++;
  }
}



/**
   function computes internal fluxes caused by mechanical effects
   
   @param eid - element id
   @param ifl - %vector of internal fluxes
   
   JK, 17.7.2005
*/
void barelc::res_lower_internal_fluxes (long eid,vector &ifl)
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





/**
   function computes contributions to the right-hand side - volume integral
      
   \int_{Omega} B^T D dOmega
   
   @param tmv - transmission %vector of one matter
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)
   
   TKr, 09/03/2011 - must be corrected
*/
/* void barelc::volume_rhs_vector (long lcid,long eid,long ri,long ci,vector &vrhs)
   {
   long i,ii;
   double area,xi,ww,jac;
   vector x(nnemp),w(intordvum[ri][ci]),gp(intordvum[ri][ci]),gx(nnemp),gy(nnemp),s(2);
   matrix gm(tnmcomp,mndofe),d;
   matrix km(mndofe,1);
   
   Mc->give_area (eid,area);
   Mt->give_node_coord2d (gx,gy,eid);
   gauss_points (gp.a,w.a,intordvum[ri][ci]);
   giveloccoord (gx,gy,x);
   
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
   Bar2d->geom_matrix (gm,gx,gy,jac);
   }
   if (Cp->bb==quad_lin){      
   //  computation of direction vector
   //Barq2d->dirvect (s,gx,gy);
   
   Barq2d->geom_matrix (gm,gx,gy,xi,jac);
   }
   if (Cp->bb==quad_quad){
   //  computation of direction vector
   //Barq2d->dirvect (s,gx,gy);
   
   Barq2d->geom_matrix (gm,gx,gy,xi,jac);
   }
   
   //  matrix of conductivity of the material
   allocm(tnmcomp,1,d);
   
   Cmu->volume_rhs1 (d,ii,ri,ci,tnmcomp);
   
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

   TKr, 09/03/2011 - must be corrected
*/
/* void barelc::res_volume_rhs_vector (vector &f,long eid,long lcid)
   {
   long i,*cn;
   vector lf;
   
   for (i=0;i<mnb;i++){
   
   cn = new long [mndofe];
   allocv (mndofe,lf);
   
   fillv (0.0,lf);
   
   if (Cp->bb==lin_lin)
   Lbt->codnum (cn,i);//mechanical code numbers
   else
   Qbt->codnum (cn,i);//mechanical code numbers 
   
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
/* void barelc::mefel_metr (long eid)
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
/* void barelc::trfel_metr (long eid)
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
