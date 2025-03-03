#include <stdlib.h>
#include "axisymc.h"
#include "vector.h"
#include "intp.h"
#include "global.h"
#include "globalc.h"
#include "axisymqq.h"
#include "globmatc.h"

axisymc::axisymc (void)
{
  //  the number of DOFs
  //  there are 8 nodes in mechanics and 4 nodes in transport process
  //  therefore 8*2+4=20
  ndofe = 20;
  
  //  the number of DOFs in mechanics
  mndofe=16;
  //  the number of DOFs in transport
  tndofe=4;
  
  //  the number of nodes for mechanics
  mnne=8;
  //  the number of nodes for transport
  tnne=4;
  
  //  the number of edges
  ned=4;
  //  the number of nodes on an edge for transport
  tnned=2;

  //  integration
  intordlin = 2;
  intordquad = 3;
  
  //  the number of integration points
  nip = 13;
  
  //  the number of displacement components
  ncompdispl = 2;
  //  the number of strain/stress components
  ncompstr = 4;
  //  strain/stress state
  ssst = Asymqq->ssst;
  //  the number of transported media
  ntm = 1;
  //  the number of gradient/flux components
  ncompgrad = 2;
  
  
  //  there are two contributions: mechanics and transport process
  ordering = new long* [2];
  ordering[0] = new long [mndofe];
  ordering[1] = new long [tndofe];
  
  ordering[0][0]  = 1;
  ordering[0][1]  = 2;
  ordering[0][2]  = 3;
  ordering[0][3]  = 4;
  ordering[0][4]  = 5;
  ordering[0][5]  = 6;
  ordering[0][6]  = 7;
  ordering[0][7]  = 8;
  ordering[0][8]  = 9;
  ordering[0][9]  = 10;
  ordering[0][10] = 11;
  ordering[0][11] = 12;
  ordering[0][12] = 13;
  ordering[0][13] = 14;
  ordering[0][14] = 15;
  ordering[0][15] = 16;

  ordering[1][0] = 17;
  ordering[1][1] = 18;
  ordering[1][2] = 19;
  ordering[1][3] = 20;

}

axisymc::~axisymc (void)
{
  delete [] ordering[0];
  delete [] ordering[1];
  delete [] ordering;
}

/**
   function assembles code numbers for particular blocks
   
   @param cn - code numbers
   @param ri - row index
   
   JK, 11.4.2019
*/
void axisymc::codnum (long *cn,long ri)
{
  long i,n;
  
  if (ri==0){
    n=mndofe;
  }
  if (ri==1){
    n=tndofe;
  }
  
  for (i=0;i<n;i++){
    cn[i]=ordering[ri][i];
  }
}

void axisymc::stiffness_matrix (long /*lcid*/,long eid,matrix &sm)
{
  long i,j,ipp;
  double xi,eta,jac,r;
  ivector nodes(ASTCKIVEC(mnne));
  vector x(ASTCKVEC(mnne)),y(ASTCKVEC(mnne)),w(ASTCKVEC(intordquad)),gp(ASTCKVEC(intordquad));
  matrix gm(ASTCKMAT(ncompstr,mndofe)),d(ASTCKMAT(ncompstr,ncompstr));
  
  //CMt->give_elemnodes (eid,nodes);
  //Mt->give_node_coord2d (x,y,eid);
  Gtm->give_node_coord2d (x,y,eid);

  fillm (0.0,sm);
  
  gauss_points (gp.a,w.a,intordquad);
  
  //  number of first integration point
  //  there are 3x3=9 integration points for quadratic approximation functions
  //  and 2x2=4 integration points for linear approximation functions
  ipp=Ct->elements[eid].ipp;
  
  for (i=0;i<intordquad;i++){
    xi=gp[i];
    for (j=0;j<intordquad;j++){
      eta=gp[j];
      
      //  geometric matrix
      Asymqq->geom_matrix (gm,x,y,xi,eta,jac);
      
      if(jac < 0.0){
	print_err("wrong numbering of nodes on element number %ld, negative volume! jac = %e", __FILE__, __LINE__, __func__, eid+1, jac);
	abort();
      }
      
      //  matrix of stiffness of the material
      Cm->matstiff (d,ipp);
      
      check_math_errel(eid);
     
      r = Asymqq->approx (xi,eta,x);
      jac*=w[i]*w[j]*r;
      
      //  contribution to the stiffness matrix of the element
      bdbjac (sm,gm,d,gm,jac);
      
      ipp++;
      
    }
  }
  
  /*
  //  transformation of stiffness matrix
  transf = Mt->locsystems (nodes);
  if (transf>0){
    matrix tmat (ASTCKMAT(ndofe,ndofe));
    transf_matrix (nodes,tmat);
    glmatrixtransf (sm,tmat);
  }
  */
}


void axisymc::conductivity_matrix (long /*lcid*/,long eid,matrix &km)
{
  long i,j,ipp;
  double xinp,xi,eta,ww1,ww2,jac;
  //ivector nodes(nne);
  vector x(ASTCKVEC(tnne)),y(ASTCKVEC(tnne)),xx(ASTCKVEC(mnne)),yy(ASTCKVEC(mnne)),w(ASTCKVEC(intordlin)),gp(ASTCKVEC(intordlin));
  matrix gm(ASTCKMAT(ncompgrad,tndofe)),d(ASTCKMAT(ncompgrad,ncompgrad)),n(ASTCKMAT(1,tndofe));
  
  //Tt->give_elemnodes (eid,nodes);
  
  //  this function returns all nodes on element
  Gtm->give_node_coord2d (xx,yy,eid);
  //  only the corner nodes are used in this part
  x[0]=xx[0];  y[0]=yy[0];
  x[1]=xx[1];  y[1]=yy[1];
  x[2]=xx[2];  y[2]=yy[2];
  x[3]=xx[3];  y[3]=yy[3];
  
  
  gauss_points (gp.a,w.a,intordlin);

  fillm (0.0,km);
  
  //  number of first integration point
  //  there are 3x3=9 integration points for quadratic approximation functions
  //  and 2x2=4 integration points for linear approximation functions
  ipp=Ct->elements[eid].ipp+9;
  
  for (i=0;i<intordlin;i++){
    xi=gp[i];  ww1=w[i];
    for (j=0;j<intordlin;j++){
      eta=gp[j];  ww2=w[j];
      
      //  matrix of gradients
      Lqat->grad_matrix (gm,x,y,xi,eta,jac);
      
      //  matrix of conductivity of the material
      Cm->matcond (d,ipp);
      check_math_errel(eid);

      xinp = Lqat->approx (xi,eta,x);

      jac*=xinp*ww1*ww2;

      //  contribution to the conductivity matrix of the element
      bdbj (km.a,gm.a,d.a,jac,gm.m,gm.n);
      
      //convective terms
      //reallocm(1,ncomp,d);
      
      //Tm->matcond2(d,ipp);
      //Lqat->bf_matrix (n, xi, eta);
      //bdbjac(km, n, d, gm, jac);	
      
      ipp++;  
    }
  }
  
  /*
  if ((Tt->elements[eid].transi[lcid]==3) || (Tt->elements[eid].transi[lcid]==4)){
    transmission_matrix (lcid,eid,ri,ci,km);
  }
  */
}

/**
   function computes capacity %matrix of 2D problems for one transported matter
   finite element with bilinear approximation functions

   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param cm - capacity %matrix

   JK, 3.4.2019
*/
void axisymc::c_pp_matrix (long eid,matrix &cm)
{
  long i,j,ipp;
  double jac,xi,eta,w1,w2,xinp,c;
  //ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(tnne)),y(ASTCKVEC(tnne)),xx(ASTCKVEC(mnne)),yy(ASTCKVEC(mnne)),w(ASTCKVEC(intordquad)),gp(ASTCKVEC(intordquad));
  matrix n(ASTCKMAT(1,tndofe));

  //  this function returns all nodes on element
  Gtm->give_node_coord2d (xx,yy,eid);
  //  only the corner nodes are used in this part
  x[0]=xx[0];  y[0]=yy[0];
  x[1]=xx[1];  y[1]=yy[1];
  x[2]=xx[2];  y[2]=yy[2];
  x[3]=xx[3];  y[3]=yy[3];
  
  gauss_points (gp.a,w.a,intordquad);
  
  fillm (0.0,cm);

  ipp=Ct->elements[eid].ipp;

  //Tt->give_elemnodes (eid,nodes);
  

  //  number of first integration point
  //  there are 3x3=9 integration points for quadratic approximation functions
  //  and 2x2=4 integration points for linear approximation functions
  ipp=Ct->elements[eid].ipp;
  
  for (i=0;i<intordquad;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordquad;j++){
      eta=gp[j];  w2=w[j];
      
      jac_2d (jac,x,y,xi,eta);
      Lqat->bf_matrix (n,xi,eta);
      
      c=Cm->c_pp_coeff (ipp);
      
      check_math_errel(eid);
      
      xinp = Lqat->approx (xi,eta,x);
      jac*=w1*w2*xinp*c;
      
      nnj (cm.a,n.a,jac,n.m,n.n);
      ipp++;
    }
  }

}

/**
   function computes capacity %matrix of 2D problems for one transported matter
   finite element with bilinear approximation functions

   @param eid - number of element
   @param ri,ci - row and column indeces of the computed block in the resulting %matrix
   @param cm - capacity %matrix

   JK, 3.4.2019
*/
void axisymc::c_up_matrix (long eid,matrix &cm)
{
  long i,j,ipp;
  double jac,xi,eta,w1,w2,xinp,c;
  //ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(tnne)),y(ASTCKVEC(tnne)),xx(ASTCKVEC(mnne)),yy(ASTCKVEC(mnne)),w(ASTCKVEC(intordquad)),gp(ASTCKVEC(intordquad));
  matrix n(ASTCKMAT(1,tndofe)),gm(ASTCKMAT(ncompstr,mndofe)),m(3,1);
  
  m[0][0]=1.0;
  m[1][0]=1.0;
  m[2][0]=0.0;
  
  //  this function returns all nodes on element
  Gtm->give_node_coord2d (xx,yy,eid);
  //  only the corner nodes are used in this part
  x[0]=xx[0];  y[0]=yy[0];
  x[1]=xx[1];  y[1]=yy[1];
  x[2]=xx[2];  y[2]=yy[2];
  x[3]=xx[3];  y[3]=yy[3];
  
  gauss_points (gp.a,w.a,intordquad);
  
  fillm (0.0,cm);

  ipp=Ct->elements[eid].ipp;

  //Tt->give_elemnodes (eid,nodes);
  

  //  number of first integration point
  //  there are 3x3=9 integration points for quadratic approximation functions
  //  and 2x2=4 integration points for linear approximation functions
  ipp=Ct->elements[eid].ipp;
  
  for (i=0;i<intordquad;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordquad;j++){
      eta=gp[j];  w2=w[j];
      
      jac_2d (jac,x,y,xi,eta);
      Lqat->bf_matrix (n,xi,eta);
      
      c=Cm->c_up_coeff (ipp);
      
      check_math_errel(eid);
      
      xinp = Lqat->approx (xi,eta,x);
      jac*=w1*w2*xinp*c;
      
      bdbjac (cm,gm,m,n,jac);
      ipp++;
    }
  }

}

/**
   function computes capacity %matrix of 2D problems for one transported matter
   finite element with bilinear approximation functions

   @param eid - number of element
   @param ri,ci - row and column indeces of the computed block in the resulting %matrix
   @param cm - capacity %matrix

   JK, 3.4.2019
*/
void axisymc::c_pu_matrix (long eid,matrix &cm)
{
  long i,j,ipp;
  double jac,xi,eta,w1,w2,xinp,c;
  //ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(tnne)),y(ASTCKVEC(tnne)),xx(ASTCKVEC(mnne)),yy(ASTCKVEC(mnne)),w(ASTCKVEC(intordquad)),gp(ASTCKVEC(intordquad));
  matrix n(ASTCKMAT(1,tndofe)),gm(ASTCKMAT(ncompstr,mndofe)),m(1,3);
  
  m[0][0]=1.0;
  m[0][1]=1.0;
  m[0][2]=0.0;

  //  this function returns all nodes on element
  Gtm->give_node_coord2d (xx,yy,eid);
  //  only the corner nodes are used in this part
  x[0]=xx[0];  y[0]=yy[0];
  x[1]=xx[1];  y[1]=yy[1];
  x[2]=xx[2];  y[2]=yy[2];
  x[3]=xx[3];  y[3]=yy[3];
  
  gauss_points (gp.a,w.a,intordquad);
  
  fillm (0.0,cm);

  ipp=Ct->elements[eid].ipp;

  //Tt->give_elemnodes (eid,nodes);
  

  //  number of first integration point
  //  there are 3x3=9 integration points for quadratic approximation functions
  //  and 2x2=4 integration points for linear approximation functions
  ipp=Ct->elements[eid].ipp;
  
  for (i=0;i<intordquad;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordquad;j++){
      eta=gp[j];  w2=w[j];
      
      jac_2d (jac,x,y,xi,eta);
      //  matrix of approximation functions for pressure
      Lqat->bf_matrix (n,xi,eta);
      //  strain-displacement matrix
      Asymqq->geom_matrix (gm,xx,yy,xi,eta,jac);
      
      c=Cm->c_pu_coeff (ipp);
      
      check_math_errel(eid);
      
      xinp = Lqat->approx (xi,eta,x);
      jac*=w1*w2*xinp*c;
      
      bdbjac (cm,n,m,gm,jac);
      ipp++;
    }
  }

}

/**
   function computes zero order %matrix of 2D problems for one transported matter
   zero order %matrix multiplies the vector of nodal values
   
   @param eid - number of element
   @param km - zero order %matrix

   JK, 11.4.2019
*/
void axisymc::zero_order_matrix (long eid,matrix &km)
{
  long lcid=0;
  ivector rcn;
  matrix m;
  
  //  conductivity matrix
  reallocm (tndofe,tndofe,m);
  conductivity_matrix (lcid,eid,m);
  reallocv (tndofe,rcn);
  codnum (rcn.a,1);
  mat_localize (km,m,rcn.a,rcn.a);
 
}

/**
   function computes first order %matrix of 2D problems for one transported matter
   first order %matrix multiplies the vector of time derivatives of nodal values
   
   @param eid - number of element
   @param cm - first order %matrix

   JK, 11.4.2019
*/
void axisymc::first_order_matrix (long eid,matrix &cm)
{
  long lcid=0;
  ivector rcn,ccn;
  matrix m;
  
  //  stiffness matrix
  reallocm (mndofe,mndofe,m);
  stiffness_matrix (lcid,eid,m);
  reallocv (mndofe,rcn);
  codnum (rcn.a,0);
  mat_localize (cm,m,rcn.a,rcn.a);
  
  //  permeability matrix
  reallocm (tndofe,tndofe,m);
  c_pp_matrix (eid,m);
  reallocv (tndofe,rcn);
  codnum (rcn.a,1);
  mat_localize (cm,m,rcn.a,rcn.a);
  
  //  coupling matrix C_up
  reallocm (mndofe,tndofe,m);
  c_up_matrix (eid,m);
  reallocv (mndofe,rcn);
  reallocv (tndofe,ccn);
  codnum (rcn.a,0);
  codnum (ccn.a,1);
  mat_localize (cm,m,rcn.a,ccn.a);
  

  //  coupling matrix C_pu
  reallocm (tndofe,mndofe,m);
  c_pu_matrix (eid,m);
  reallocv (tndofe,rcn);
  reallocv (mndofe,ccn);
  codnum (rcn.a,1);
  codnum (ccn.a,0);
  mat_localize (cm,m,rcn.a,ccn.a);

}


/**
   function computes strains in all integration points on element
   
   @param eid - element id
   
   JK, 11.4.2019
*/
void axisymc::ip_strains (long eid)
{
  long i,j,ipp,lcid=0;
  double xi,eta,jac;
  vector x(ASTCKVEC(mnne)),y(ASTCKVEC(mnne)),r(ASTCKVEC(mndofe)),gp,w,eps(ASTCKVEC(ncompstr));
  matrix gm(ASTCKMAT(ncompstr,mndofe));
  
  eldispl_fc (lcid,eid,r.a);
  
  // ************************************************************
  //  values for the integration points for quadratic functions
  // ************************************************************
  reallocv (RSTCKVEC(intordquad,gp));
  reallocv (RSTCKVEC(intordquad,w));
  
  gauss_points (gp.a,w.a,intordquad);
  
  //  id of the first integration point
  ipp=Ct->elements[eid].ipp;
  
  for (i=0;i<intordquad;i++){
    xi=gp[i];
    for (j=0;j<intordquad;j++){
      eta=gp[j];
      
      Asymqq->geom_matrix (gm,x,y,xi,eta,jac);
      mxv (gm,r,eps);
      
      Cm->storestrain (ipp,eps);
      ipp++;
    }
  }

  // ************************************************************
  //  values for integration points for linear functions
  // ************************************************************
  reallocv (RSTCKVEC(intordlin,gp));
  reallocv (RSTCKVEC(intordlin,w));
  
  gauss_points (gp.a,w.a,intordlin);
  for (i=0;i<intordlin;i++){
    xi=gp[i];
    for (j=0;j<intordlin;j++){
      eta=gp[j];
      
      Asymqq->geom_matrix (gm,x,y,xi,eta,jac);
      mxv (gm,r,eps);
      
      Cm->storestrain (ipp,eps);
      ipp++;
    }
  }
  
}
