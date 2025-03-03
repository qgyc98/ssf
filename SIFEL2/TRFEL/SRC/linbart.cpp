/*
  File:     linbart.cpp
  Author:   Jaroslav Kruis, 31.3.2003
  Purpose:  onedimensional element with linear approximation functions
*/

#include "globalt.h"
#include "linbart.h"
#include "genfile.h"
#include "globmatt.h"
#include "loadcaset.h"
#include "loadelt.h"
#include "elemswitcht.h"

linbart::linbart (void)
{
  long i;
  
  //  number of nodes on element
  nne=2;
  //  geometric problem dimension (1D)
  ncomp=1;
  //  number of end nodes
  nen=2;
  //  the number of edges
  ned=1;
  //  the number of nodes on one edge
  nned=0;

  //  number of transported variables
  ntm=Tp->ntm;


  dofe = new long* [ntm];
  nip = new long* [ntm];
  intordkm = new long* [ntm];
  intordcm = new long* [ntm];
  ordering = new long* [ntm];
  for (i=0;i<ntm;i++){
    dofe[i] = new long [ntm];
    nip[i] = new long [ntm];
    intordkm[i] = new long [ntm];
    intordcm[i] = new long [ntm];
    ordering[i] = new long [nne];
  }
  
  
  switch (Tp->tmatt){
  case onemedium:{
    ordering[0][0]=1;  ordering[0][1]=2;
    dofe[0][0]=2;  intordkm[0][0]=3;  intordcm[0][0]=3;  nip[0][0]=6;
    ndofe=2;  napfun=1;
    break;
  }
  case twomediacoup:{
    ordering[0][0]=1;  ordering[0][1]=3;
    ordering[1][0]=2;  ordering[1][1]=4;

    intordkm[0][0]=1;  intordkm[0][1]=1;  intordkm[1][0]=1;  intordkm[1][1]=1;
    intordcm[0][0]=2;  intordcm[0][1]=2;  intordcm[1][0]=2;  intordcm[1][1]=2;
    
    if (Tp->savemode==0){
      nip[0][0]=3;       nip[0][1]=3;       nip[1][0]=3;       nip[1][1]=3;
    }
    if (Tp->savemode==1){
      nip[0][0]=3;       nip[0][1]=0;       nip[1][0]=0;       nip[1][1]=0;
    }
    
    dofe[0][0]=2;  dofe[0][1]=2;  dofe[1][0]=2;  dofe[1][1]=2;
    ndofe=4;  napfun=2;
    break;
  }
  case threemediacoup:{
    ordering[0][0]=1;   ordering[0][1]=4;
    ordering[1][0]=2;   ordering[1][1]=5;
    ordering[2][0]=3;   ordering[2][1]=6;

    intordkm[0][0]=1;  intordkm[0][1]=1;  intordkm[0][2]=1;
    intordkm[1][0]=1;  intordkm[1][1]=1;  intordkm[1][2]=1;
    intordkm[2][0]=1;  intordkm[2][1]=1;  intordkm[2][2]=1;

    intordcm[0][0]=2;  intordcm[0][1]=2;  intordcm[0][2]=2;
    intordcm[1][0]=2;  intordcm[1][1]=2;  intordcm[1][2]=2;
    intordcm[2][0]=2;  intordcm[2][1]=2;  intordcm[2][2]=2;

    if (Tp->savemode==0){
      nip[0][0]=3;  nip[0][1]=3;  nip[0][2]=3;
      nip[1][0]=3;  nip[1][1]=3;  nip[1][2]=3;
      nip[2][0]=3;  nip[2][1]=3;  nip[2][2]=3;
    }
    if (Tp->savemode==1){
      nip[0][0]=3;  nip[0][1]=0;  nip[0][2]=0;
      nip[1][0]=0;  nip[1][1]=0;  nip[1][2]=0;
      nip[2][0]=0;  nip[2][1]=0;  nip[2][2]=0;
    }
    
    dofe[0][0]=2;  dofe[0][1]=2;  dofe[0][2]=2;
    dofe[1][0]=2;  dofe[1][1]=2;  dofe[1][2]=2;
    dofe[2][0]=2;  dofe[2][1]=2;  dofe[2][2]=2;

    ndofe=6;  napfun=3;
    break;
  }
  case fourmediacoup:{
    ordering[0][0]=1;   ordering[0][1]=5;
    ordering[1][0]=2;   ordering[1][1]=6;
    ordering[2][0]=3;   ordering[2][1]=7;
    ordering[3][0]=4;   ordering[3][1]=8;

    intordkm[0][0]=1;  intordkm[0][1]=1;  intordkm[0][2]=1;  intordkm[0][3]=1;
    intordkm[1][0]=1;  intordkm[1][1]=1;  intordkm[1][2]=1;  intordkm[1][3]=1;
    intordkm[2][0]=1;  intordkm[2][1]=1;  intordkm[2][2]=1;  intordkm[2][3]=1;
    intordkm[3][0]=1;  intordkm[3][1]=1;  intordkm[3][2]=1;  intordkm[3][3]=1;

    intordcm[0][0]=2;  intordcm[0][1]=2;  intordcm[0][2]=2;  intordcm[0][3]=2;
    intordcm[1][0]=2;  intordcm[1][1]=2;  intordcm[1][2]=2;  intordcm[1][3]=2;
    intordcm[2][0]=2;  intordcm[2][1]=2;  intordcm[2][2]=2;  intordcm[2][3]=2;
    intordcm[3][0]=2;  intordcm[3][1]=2;  intordcm[3][2]=2;  intordcm[3][3]=2;

    if (Tp->savemode==0){
      nip[0][0]=3;  nip[0][1]=3;  nip[0][2]=3;  nip[0][3]=3;
      nip[1][0]=3;  nip[1][1]=3;  nip[1][2]=3;  nip[1][3]=3;
      nip[2][0]=3;  nip[2][1]=3;  nip[2][2]=3;  nip[2][3]=3;
      nip[3][0]=3;  nip[3][1]=3;  nip[3][2]=3;  nip[3][3]=3;
    }
    if (Tp->savemode==1){
      nip[0][0]=3;  nip[0][1]=0;  nip[0][2]=0;  nip[0][3]=0;
      nip[1][0]=0;  nip[1][1]=0;  nip[1][2]=0;  nip[1][3]=0;
      nip[2][0]=0;  nip[2][1]=0;  nip[2][2]=0;  nip[2][3]=0;
      nip[3][0]=0;  nip[3][1]=0;  nip[3][2]=0;  nip[3][3]=0;
    }
    
    dofe[0][0]=2;  dofe[0][1]=2;  dofe[0][2]=2;  dofe[0][3]=2;
    dofe[1][0]=2;  dofe[1][1]=2;  dofe[1][2]=2;  dofe[1][3]=2;
    dofe[2][0]=2;  dofe[2][1]=2;  dofe[2][2]=2;  dofe[2][3]=2;
    dofe[3][0]=2;  dofe[3][1]=2;  dofe[3][2]=2;  dofe[3][3]=2;

    ndofe=8;  napfun=4;
    break;
  }
  default:{
    print_err("unknown type of transported matter is required",__FILE__,__LINE__,__func__);
  }
  }
}

linbart::~linbart (void)
{
  long i;

  for (i=0;i<ntm;i++){
    delete [] dofe[i];
    delete [] nip[i];
    delete [] intordkm[i];
    delete [] intordcm[i];
    delete [] ordering[i];
  }
  delete [] dofe;
  delete [] nip;
  delete [] intordkm;
  delete [] intordcm;
  delete [] ordering;
}

/**
   function assembles code numbers for one medium
   
   @param cn - code numbers
   @param ri - row index
   
   JK
*/
void linbart::codnum (long *cn,long ri)
{
  long i;
  for (i=0;i<nne;i++){
    cn[i]=ordering[ri][i];
  }
}

/**
   function approximates function defined by nodal values

   @param xi - natural coordinate on element
   @param nodval - %vector of nodal values
   
   JK, 31.3.2002
*/
double linbart::approx (double xi,vector &nodval)
{
  double f;
  vector bf(nne);
  
  bf_lin_1d (bf.a,xi);
  
  scprd (bf,nodval,f);
  
  return f;
}

/**
   function computes values in integration points from nodal values

   @param eid - element id

   JK, 31.3.2003
*/
void linbart::intpointval (long eid)
{
  long i,k,ii,jj,ipp;
  double xi,val;
  vector r(ASTCKVEC(ndofe)),t(ASTCKVEC(nne)),gp,w;
  
  elemvalues(eid, r);
  
  for (k=0;k<ntm;k++){
    
    //  nodal values of a single variable
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
	
	//  integration points for the conductivity matrix
	
	reallocv (RSTCKVEC(intordkm[ii][jj],gp));
	reallocv (RSTCKVEC(intordkm[ii][jj],w));
	gauss_points (gp.a,w.a,intordkm[ii][jj]);
	
	ipp=Tt->elements[eid].ipp[ii][jj];
	for (i=0;i<intordkm[ii][jj];i++){
	  xi=gp[i];
	  
	  val = approx (xi,t);
	  Tm->ip[ipp].av[k]=val;
	  ipp++;
	}
	
	//  integration points for the capacity matrix
	
	reallocv (RSTCKVEC(intordcm[ii][jj],gp));
	reallocv (RSTCKVEC(intordcm[ii][jj],w));
	gauss_points (gp.a,w.a,intordcm[ii][jj]);
	
	for (i=0;i<intordcm[ii][jj];i++){
	  xi=gp[i];
	  
	  val = approx (xi,t);
	  Tm->ip[ipp].av[k]=val;
	  ipp++;
	}
	
	
	if (Tp->savemode==1)  break;
      }
      if (Tp->savemode==1)  break;
    }
  }
}


/**
   function computes values in integration points from nodal values
   this function is used for arbitrary variable, not only for
   variables used as unknowns in the problem
   
   @param eid    - element id
   @param nodval - nodal values
   @param ipval  - value at integration points
   
   @retval Function returns approximated values at integration points in the %vector ipval.

   12/06/2012 TKr according to JK
*/
void linbart::intpointval (long eid,vector &nodval,vector &ipval)
{
  long i,kk,ii,jj,ipp;
  double xi;
  vector gp,w;
  
  kk = 0;
  
  for (ii=0;ii<Tp->ntm;ii++){
    for (jj=0;jj<Tp->ntm;jj++){
      
      //  interpolation to integration points of conductivity matrix
      reallocv (RSTCKVEC(intordkm[ii][jj],gp));
      reallocv (RSTCKVEC(intordkm[ii][jj],w));
      gauss_points (gp.a,w.a,intordkm[ii][jj]);
      
      ipp=Tt->elements[eid].ipp[ii][jj];
      for (i=0;i<intordkm[ii][jj];i++){
        xi=gp[i];
        //  value in integration point
        ipval[kk] = approx (xi,nodval);
	
        ipp++;
        kk++;
      }
      
      
      //  interpolation to integration points of capacity matrix
      reallocv (RSTCKVEC(intordcm[ii][jj],gp));
      reallocv (RSTCKVEC(intordcm[ii][jj],w));
      gauss_points (gp.a,w.a,intordcm[ii][jj]);
      
      for (i=0;i<intordcm[ii][jj];i++){
        xi=gp[i];
        //  value in integration point
        ipval[kk] = approx (xi,nodval);

        ipp++;
        kk++;
      }
	
	
      if (Tp->savemode==1)  break;
    }
    if (Tp->savemode==1)  break;
  }
}


/**
   function computes initial values in integration points from initial nodal values

   @param eid - element id

   TKo, 4.7.2018
*/
void linbart::initintpointval (long eid)
{
  long i,k,ii,jj,ipp,ndofn,cndofn;
  double xi,val;
  ivector enod(ASTCKIVEC(nne));
  vector rr, r(ASTCKVEC(ndofe)), t(ASTCKVEC(nne)), gp, w;
  
  Tt->give_elemnodes (eid,enod);

  for(i=0, cndofn=0; i<nne; i++)
  {
    ndofn = Tt->give_ndofn(enod[i]);
    // make reference rr to the vector of nodal values on element
    makerefv(rr, r.a+cndofn, ndofn);
    // get initial nodal values and store them in the vector r with the help of reference rr
    initnodval2(enod[i], rr);
    cndofn += ndofn;
  }
  
  for (k=0;k<ntm;k++){
    
    //  nodal values of a single variable
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
	
	//  integration points for the conductivity matrix
	
	reallocv (RSTCKVEC(intordkm[ii][jj],gp));
	reallocv (RSTCKVEC(intordkm[ii][jj],w));
	gauss_points (gp.a,w.a,intordkm[ii][jj]);
	
	ipp=Tt->elements[eid].ipp[ii][jj];
	for (i=0;i<intordkm[ii][jj];i++){
	  xi=gp[i];
	  
	  val = approx (xi,t);
	  Tm->ip[ipp].av[k]=val;
	  ipp++;
	}
	
	//  integration points for the capacity matrix
	
	reallocv (RSTCKVEC(intordcm[ii][jj],gp));
	reallocv (RSTCKVEC(intordcm[ii][jj],w));
	gauss_points (gp.a,w.a,intordcm[ii][jj]);
	
	for (i=0;i<intordcm[ii][jj];i++){
	  xi=gp[i];
	  
	  val = approx (xi,t);
	  Tm->ip[ipp].av[k]=val;
	  ipp++;
	}
	
	
	if (Tp->savemode==1)  break;
      }
      if (Tp->savemode==1)  break;
    }
  }
}


/**
   function computes values in integration points from nodal values

   @param eid - element id

   JK, 31.3.2003
*/
void linbart::intpointgrad (long eid)
{
  long i,ii,jj,k,ipp;
  double xi,jac;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),r(ASTCKVEC(ndofe)),t(ASTCKVEC(nne));
  vector gp,w,grad(ASTCKVEC(ncomp));
  matrix gm(ncomp,nne);
  
  Tt->give_node_coord2d (x,y,eid);
  elemvalues(eid, r);

  for (k=0;k<Tp->ntm;k++){
    
    //  nodal values of a single variable
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
	
	//  integration points for the conductivity matrix
	
	reallocv (RSTCKVEC(intordkm[ii][jj],gp));
	reallocv (RSTCKVEC(intordkm[ii][jj],w));
	gauss_points (gp.a,w.a,intordkm[ii][jj]);
	
	ipp=Tt->elements[eid].ipp[ii][jj];
	for (i=0;i<intordkm[ii][jj];i++){
	  xi=gp[i];
	  
	  //  matrix of gradients
	  grad_matrix (gm,x,xi,jac);
	  mxv (gm,t,grad);
	  Tm->storegrad (k,ipp,grad);
	  ipp++;
	}
	
	//  integration points for the capacity matrix
	
	reallocv (RSTCKVEC(intordcm[ii][jj],gp));
	reallocv (RSTCKVEC(intordcm[ii][jj],w));
	gauss_points (gp.a,w.a,intordcm[ii][jj]);
	
	for (i=0;i<intordcm[ii][jj];i++){
	  xi=gp[i];
	  
	  //  matrix of gradients
	  grad_matrix (gm,x,xi,jac);
	  mxv (gm,t,grad);
	  Tm->storegrad (k,ipp,grad);
	  ipp++;
	}
	
	if (Tp->savemode==1)  break;
      }
      if (Tp->savemode==1)  break;
    }
  }
}

/**
   function approximates nodal values of array other to integration points
   
   @param eid - element id
   
   JK, 17.9.2005
*/
void linbart::intpointother (long eid)
{
  long i, k, ii, jj, ipp, ncompo, nodid;
  double xi, val;
  ivector nodes(nne);
  vector t(nne), r, gp, w;
  
  //  nodes of required element
  Tt->give_elemnodes (eid,nodes);
  
  //  first node number
  nodid=nodes[0];

  //  number of components
  ncompo=Tt->nodes[nodid].ncompother;
  reallocv(RSTCKVEC(ncompo*nne, r));
  
  //  nodal values of array other
  nodalotherval (nodes, r);
  
  for (k=0;k<ncompo;k++){
    
    //  nodal values of a single variable
    for (i=0;i<nne;i++){
      t[i]=r[i*ncompo+k];
    }
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
	
	//  integration points for the conductivity matrix

	reallocv (intordkm[ii][jj],gp);
	reallocv (intordkm[ii][jj],w);
	gauss_points (gp.a,w.a,intordkm[ii][jj]);
	
	ipp=Tt->elements[eid].ipp[ii][jj];
	for (i=0;i<intordkm[ii][jj];i++){
	  xi=gp[i];
	  
	  val = approx (xi,t);
	  Tm->ip[ipp].other[k]=val;
	  ipp++;
	}
	
	//  integration points for the capacity matrix

	reallocv (intordcm[ii][jj],gp);
	reallocv (intordcm[ii][jj],w);
	gauss_points (gp.a,w.a,intordcm[ii][jj]);
	
	for (i=0;i<intordcm[ii][jj];i++){
	  xi=gp[i];
	  
	  val = approx (xi,t);
	  Tm->ip[ipp].other[k]=val;
	  ipp++;
	}
	
	if (Tp->savemode==1)  break;
      }
      if (Tp->savemode==1)  break;
    }
  }
}

/**
   function assembles approximation or test functions into a %matrix
   
   @param w - %matrix of test functions
   @param xi - natural coordinate
   @param v - velocity of advection
   
   2. 10. 2016, JK
*/
void linbart::btf_matrix (matrix &n,double xi,double v)
{
  if (Tp->advect==0)
    bf_matrix (n,xi);
  if (Tp->advect==1)
    tf_matrix (n,xi,v);
}


/**
   function assembles %matrix of base functions
   
   @param n - %matrix of base functions
   @param xi - natural coordinate
   
   JK, 31.3.2002
*/
void linbart::bf_matrix (matrix &n,double xi)
{
  nullm (n);
  bf_lin_1d (n.a,xi);
}

/**
   function assembles %matrix of test functions
   
   @param w - %matrix of test functions
   @param xi - natural coordinate
   @param v - velocity of advection
   
   JK, 9. 9. 2016
*/
void linbart::tf_matrix (matrix &w,double xi,double v)
{
  double alpha=1.0;
  matrix n(1,dofe[0][0]);
  nullm (n);
  bf_lin_1d (n.a,xi);

  double h=1.2500000000e-03;
  double nor=1.507500e-07;
  double pe=nor*h/2.0/9.52e-12;
  alpha=(exp(pe)+exp(0.0-pe))/(exp(pe)-exp(0.0-pe))-1.0/pe;
  //alpha=1.0;

  if (v<0.0){
    w[0][0]=n[0][0]+alpha/2.0;
    w[0][1]=n[0][1]-alpha/2.0;
  }else{
    w[0][0]=n[0][0]-alpha/2.0;
    w[0][1]=n[0][1]+alpha/2.0;
  }
}


/**
   function assembles gradient of %matrix of base functions
   
   @param gm - gradient %matrix
   @param x - array containing node coordinates
   @param xi - natural coordinate
   @param jac - jacobian
   
   JK, 31.3.2002
*/
void linbart::grad_matrix (matrix &gm,vector &x,double xi,double &jac)
{
  long i;
  vector dx(nne);
  dx_bf_lin_1d (dx.a);
  derivatives_1d (dx,jac,x,xi);
  
  for (i=0;i<nne;i++){
    gm[0][i]=dx[i];
  }
}

/**
   function computes conductivity %matrix of 1D problems for one transported matter
   finite element with linear approximation functions
   
   @param lcid - number of load case
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param km - conductivity %matrix
   
   JK, 31.3.2002
*/
void linbart::conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,ii;
  double area,xi,ww,jac;
  ivector nodes(nne);
  vector x(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]),a(nne);
  matrix gm(ncomp,dofe[ri][ci]),d;
  matrix n(1,dofe[ri][ci]);
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  nullm (km);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];
  
  for (i=0;i<intordkm[ri][ci];i++){
    xi=gp[i];  ww=w[i];
    
    //  matrix of gradients
    grad_matrix (gm,x,xi,jac);
    
    //  matrix of conductivity of the material
    reallocm(ncomp,ncomp,d);
    Tm->matcond (d,ii,ri,ci);

    //  area of cross section
    area = approx (xi,a);
    
    jac*=area*ww;
    
    //  contribution to the conductivity matrix of the element
    bdbj (km.a,gm.a,d.a,jac,gm.m,gm.n);
    
    //convective terms
    reallocm (1,ncomp,d);
    Tm->matcond2 (d,ii,ri,ci);
    btf_matrix (n, xi,d[0][0]);
    bdbjac (km, n, d, gm, jac);
    
    ii++;
  }
  
  if ((Tt->elements[eid].transi[lcid]==3) || (Tt->elements[eid].transi[lcid]==4)){
    //  additional matrix due to transmission
    transmission_matrix (lcid,eid,ri,ci,km);
  }
  
}

/**
   function computes capacity %matrix of 1D problems for one transported matter
   finite element with bilinear approximation functions
   
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param cm - capacity %matrix
   
   JK, 31.3.2002
*/
void linbart::capacity_matrix (long eid,long ri,long ci,matrix &cm)
{
  long i,ii;
  double jac,xi,ww,area,rho,c;
  ivector nodes(nne);
  vector x(nne),w(intordcm[ri][ci]),gp(intordcm[ri][ci]),a(nne);
  matrix n(1,dofe[ri][ci]);
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  Tc->give_densitye (eid,rho);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordcm[ri][ci]);
  
  nullm (cm);

  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci]+intordkm[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0]+intordkm[0][0];

  for (i=0;i<intordcm[ri][ci];i++){
    xi=gp[i];  ww=w[i];
    
    bf_matrix (n,xi);
    
    //  area of cross section
    area = approx (xi,a);
    c=Tm->capcoeff (ii,ri,ci);
    jac_1d (jac,x,xi);
    jac*=area*ww*rho*c;
    
    nnj (cm.a,n.a,jac,n.m,n.n);
    ii++;
  }
  
}


/**
   function computes source %vector of one matter on one element
   
   \int_{Omega} N^T N d Omega . s
   
   @param sv - source %vector of one matter
   @param nodval - array of nodal values
   @param eid - element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)
   
   JK, 31.3.2002
*/
void linbart::quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci)
{
  long i;
  double jac,xi,ww,area;
  ivector nodes(nne);
  vector x(nne),w(intordcm[ri][ci]),gp(intordcm[ri][ci]),v(dofe[ri][ci]),a(nne);
  matrix n(1,dofe[ri][ci]),nm(dofe[ri][ci],dofe[ri][ci]);
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordcm[ri][ci]);
  
  nullm (nm);
  
  for (i=0;i<intordcm[ri][ci];i++){
    xi=gp[i];  ww=w[i];
    
    bf_matrix (n,xi);
    
    jac_1d (jac,x,xi);
    
    //  area of cross section
    area = approx (xi,a);

    jac*=ww*area;
    
    nnj (nm.a,n.a,jac,n.m,n.n);
  }
  mxv (nm,nodval,v);  addv (sv,v,sv);
  
}



/**
   function computes transmission complement to the conductivity %matrix for one matter
   
   \int_{Gamma_3} N^T c_{tr} N dGamma
   
   @param lcid - number of load case
   @param eid - number of element
   @param ri,ci - row and column indices
   @param km - part of the conductivity %matrix
   
   JK, 31.3.2002
   TKr, 30.1.2004 - new added
*/
void linbart::transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,neleml,ii;
  bocontypet bc[2]; // array for boundary condition indicators

  double xi,area;
  double new_trc;//added
  long nn;//added 

  ivector nodes(nne);
  vector x(nne),trc(2),a(nne);
  matrix n(1,dofe[ri][ci]);
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  
  for (i=0;i<Tb->lc[lcid].neb;i++){
    if (Tb->lc[lcid].elemload[i].eid==eid){
      neleml=i;
      break;
    }
  }

  //  indicators of boundary conditions
  Tb->lc[lcid].elemload[neleml].give_bc (bc);
  //  transmission coefficients
  Tb->lc[lcid].elemload[neleml].give_trc (Tp->time,Tb->lc[lcid].nodval,lcid,ci,trc);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];

 
  if (bc[0]==5 || bc[0]>10){
    xi=-1.0;
    
    bf_matrix (n,xi);
    
    /***********************************************///added
    nn = nodes[0];
    
    Tm->transmission_transcoeff(new_trc,trc[0],ri,ci,nn,bc[0],ii);
    //printf("new_trc   = %e\n",new_trc);

    trc[0]=new_trc;
    /***********************************************///added
    //  area of cross section
    area = approx (xi,a);
    
    nnj (km.a,n.a,trc[0]*area,n.m,n.n);
  }

  if (bc[1]==5 || bc[1]>10){
    xi=1.0;
    
    bf_matrix (n,xi);

    /***********************************************///added
    nn = nodes[1];
    
    Tm->transmission_transcoeff(new_trc,trc[1],ri,ci,nn,bc[1],ii);
    //printf("new_trc   = %e\n",new_trc);
    
    trc[1]=new_trc;
    /***********************************************///added
    //  area of cross section
    area = approx (xi,a);

    nnj (km.a,n.a,trc[1]*area,n.m,n.n);
  }
}



/**
   function computes contributions to the transmission %vector
   
   \int_{Gamma_3} N^T c_{tr} N dGamma * nodal_external_value
   
   @param tmv - transmission %vector of one matter
   @param lcid - load case id (corresponds to the row index)
   @param eid - element id
   @param leid - loaded element id
   @param cid - component id (corresponds to the column index)
   
   JK, 5.10.2001
   TKr, 30.1.2002 - new added
*/
void linbart::transmission_vector (vector &tmv,long lcid,long eid,long leid,long cid)
{
  long ii,nid;
  bocontypet bc[2];  //  array for boundary condition indicators
  double xi,area,new_nodval,tr;
  ivector nodes(nne);
  vector x(nne),trc(2),trcn(nne),nodval(2),av(dofe[lcid][cid]),v(dofe[lcid][cid]),a(nne);
  vector trr(2);
  matrix n(1,dofe[lcid][cid]),km(dofe[lcid][cid],dofe[lcid][cid]);
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  
  //  indicators of boundary conditions
  Tb->lc[lcid].elemload[leid].give_bc (bc);
  //  nodal values of exterior variable
  Tb->lc[lcid].elemload[leid].give_external_nodval (lcid,cid,nodval);
  //  transmission coefficients
  Tb->lc[lcid].elemload[leid].give_trc (Tp->time,Tb->lc[lcid].nodval,lcid,cid,trc);
  //  transmission/radiation coefficients
  Tb->lc[lcid].elemload[leid].give_trr (Tp->time,Tb->lc[lcid].nodval,trr);

  //  integration point number
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[lcid][cid];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];

  if (bc[0]==5 || bc[0]>10){
    nullm (km);
    xi=-1.0;
    
    bf_matrix (n,xi);

    //  node id
    nid = nodes[0];
    tr = 0.0;
    if(bc[0] == 90)
      tr = trr[0]/trc[0];
    
    
    Tm->transmission_nodval (new_nodval,nodval[0],tr,lcid,cid,nid,bc[0],ii);
    
    av[0]=new_nodval;  av[1]=0.0;
    
    //  area of cross section
    area = approx (xi,a);
    
    //fprintf (stdout,"\n linbart  trc[0]   %le",trc[0]);
    //fprintf (stdout,"\n linbart  val[0]   %le",av[0]);

    nnj (km.a,n.a,trc[0]*area,n.m,n.n);
    
    mxv (km,av,v);  addv (tmv,v,tmv);
  }
  if (bc[1]==5 || bc[1]>10){
    nullm (km);
    xi=1.0;
    
    bf_matrix (n,xi);

    /***********************************************///added
    nid = nodes[1];
    tr = 0.0;
    if(bc[0] == 90)
      tr = trr[1]/trc[1];
    
    Tm->transmission_nodval(new_nodval,nodval[1],tr,lcid,cid,nid,bc[1],ii);
   
    av[0]=0.0;  av[1]=new_nodval;
    /***********************************************/
    //  area of cross section
    area = approx (xi,a);
          
    //fprintf (stdout,"\n linbart  trc[1]   %le",trc[1]);
    //fprintf (stdout,"\n linbart  val[1]   %le",av[1]);
    nnj (km.a,n.a,trc[1]*area,n.m,n.n);

    mxv (km,av,v);  addv (tmv,v,tmv);
  }
}



/**
   function computes contribution to the convection %vector

   \int_{Gamma_2} N^T N dGamma * nodal_flux_values
   
   @param f - convection %vector of one matter
   @param lcid - load case id
   @param eid - element id
   @param leid - loaded element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)

   JK, 8.10.2001
*/
void linbart::convection_vector (vector &f,long lcid,long eid,long leid,long ri,long ci)
{
  bocontypet bc[2]; // array for boundary condition indicators
  double xi,area;
  ivector nodes(nne);
  vector x(nne),nodval(2),av(dofe[ri][ci]),v(dofe[ri][ci]),a(nne);
  matrix n(1,dofe[ri][ci]),km(dofe[ri][ci],dofe[ri][ci]);

  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  
  //  indicators of boundary conditions
  Tb->lc[lcid].elemload[leid].give_bc (bc);
  //  nodal values of exterior variable
  Tb->lc[lcid].elemload[leid].give_nodval (lcid,nodval);


  if (bc[0]==2 || bc[0]==3 || bc[0]==5){
    nullm (km);
    xi=-1.0;

    bf_matrix (n,xi);
    //  area of cross section
    area = approx (xi,a);
    nnj (km.a,n.a,area,n.m,n.n);

    av[0]=nodval[0];  av[1]=0.0;
    mxv (km,av,v);  addv (f,v,f);
  }
  if (bc[1]==2 || bc[1]==3 || bc[1]==5){
    nullm (km);
    xi=1.0;
    bf_matrix (n,xi);
    //  area of cross section
    area = approx (xi,a);
    nnj (km.a,n.a,area,n.m,n.n);

    av[0]=0.0;  av[1]=nodval[1];
    mxv (km,av,v);  addv (f,v,f);
  }
}

/**
   function computes internal fluxes of 1D problems for one transported matter
   finite element with linear approximation functions
   
   @param lcid - number of load case
   @param eid - number of element
   @param ifl - %vector of internal fluxes
   
   JK, 31.3.2002
*/
void linbart::internal_fluxes (long lcid,long eid,vector &ifl)
{
  long i,ipp;
  double jac,area,xi;
  ivector nodes(nne);
  vector x(nne),w(intordkm[lcid][lcid]),gp(intordkm[lcid][lcid]),fl(ncomp),contr(dofe[lcid][lcid]),a(nne);
  matrix gm(ncomp,dofe[lcid][lcid]),d(ncomp,ncomp);
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[lcid][lcid]);
  
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[lcid][lcid];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];


  for (i=0;i<intordkm[lcid][lcid];i++){
    xi=gp[i];
    
    Tm->computenlfluxes (lcid,ipp);
    
    Tm->givefluxes (lcid,ipp,fl);

    grad_matrix (gm,x,gp[i],jac);
    
    mtxv (gm,fl,contr);
    
    //  area of cross section
    area = approx (xi,a);
    cmulv (area*jac*w[i],contr);
    
    addv (contr,ifl,ifl);
    
    ipp++;
  }
}

/**
   function computes energy accumulated in an element
   finite element with linear approximation functions
   
   @param lcid - number of load case
   @param eid - number of element
   @param ifl - %vector of internal fluxes
   
   JK, 9.4.2019
*/
void linbart::accumul_energy (long lcid,long eid,vector &acum)
{
  long i,ipp;
  double jac,area,xi,ae,v;
  ivector nodes(nne);
  vector x(nne),w(intordcm[lcid][lcid]),gp(intordcm[lcid][lcid]),contr(dofe[lcid][lcid]),a(nne);
  matrix n(1,dofe[lcid][lcid]);
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordcm[lcid][lcid]);
  
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[lcid][lcid]+intordkm[lcid][lcid];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0]+intordkm[0][0];


  for (i=0;i<intordcm[lcid][lcid];i++){
    xi=gp[i];
    
    Tm->computenlcapacities (lcid,ipp);
    
    //  velocity has to be obtained from somewhere
    v=0.0;

    btf_matrix (n,xi,v);
    
    contr[0]=n[0][0];
    contr[1]=n[0][1];

    ae = Tm->ip[ipp].cap[lcid];
    
    //  area of cross section
    area = approx (xi,a);
    jac_1d (jac,x,xi);
    cmulv (ae*area*jac*w[i],contr);
    
    addv (contr,acum,acum);
    
    ipp++;
  }
}


/**
   function computes advection %matrix of 1D problems for one transported matter
   finite element with linear approximation functions
   
   @param lcid - number of load case
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param hm - advection %matrix
   
   JK, 1. 5. 2016
*/
void linbart::advection_matrix (long /*lcid*/,long eid,long ri,long ci,matrix &hm)
{
  long i,ii;
  double area,xi,ww,jac;
  ivector nodes(nne);
  vector x(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]),a(nne);
  matrix v(1,ncomp);
  matrix gm(ncomp,dofe[ri][ci]), n(1,dofe[ri][ci]);
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  nullm (hm);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];
  
  for (i=0;i<intordkm[ri][ci];i++){
    xi=gp[i];  ww=w[i];
    
    //  matrix of approximation functions
    bf_matrix (n,xi);
    //  matrix of gradients
    grad_matrix (gm,x,xi,jac);
    
    //  vector of velocity
    Tm->matcond2(v,ii,ri,ci);
    
    //  area of cross section
    area = approx (xi,a);
    
    jac*=area*ww;
    
    //  contribution to the advection matrix of the element
    bdbjac (hm,n,v,gm,jac);
    
    ii++;
  }
  
}



/**
   function assembles resulting element conductivity %matrix

   @param eid - element id
   @param lcid - load case id
   @param km - resulting conductivity %matrix of one element

   JK, 6.1.2002
*/
void linbart::res_conductivity_matrix (long eid,long /*lcid*/,matrix &km)
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
   function computes contributions to the right-hand side - volume integral
      
   \int_{Omega} B^T D dOmega
   
   @param tmv - transmission %vector of one matter
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)
   
   TKr, 16.5.2007
*/
void linbart::volume_rhs_vector (long /*lcid*/,long eid,long ri,long ci,vector &vrhs)
{
  long i,ii;
  double area,xi,ww,jac;
  ivector nodes(nne);
  vector x(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]),a(nne);
  matrix gm(ncomp,dofe[ri][ci]),d;
  matrix km(dofe[ri][ci],1);

  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  nullm (km);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];
  
  for (i=0;i<intordkm[ri][ci];i++){
    xi=gp[i];  ww=w[i];
    
    //  matrix of gradients
    grad_matrix (gm,x,xi,jac);
    
    //  matrix of conductivity of the material
    reallocm(ncomp,1,d);

    Tm->volume_rhs (d,ii,ri,ci,ncomp);
    
    //  area of cross section
    area = approx (xi,a);

    jac*=area*ww;
    
    //  contribution to the volume_rhs integral of the element
    nnjac (km,gm,d,jac);
    
    ii++;
  }

  for (i=0;i<vrhs.n;i++){
    vrhs[i] = km[i][0];
  }

}



/**
   function computes contributions to the right-hand side - volume integral of the second type
      
   \int_{Omega} N^T const dOmega
   
   @param vrhs - volume right-hand side %vector of one matter
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)
   
   TKr, 08/06/2018
*/
void linbart::volume_rhs_vector2 (long /*lcid*/,long eid,long ri,long ci,vector &vrhs)
{
  long i,ii,ipp;
  double jac,xi,w1,area,c;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordcm[ri][ci]),gp(intordcm[ri][ci]),a(nne);
  matrix n(1,dofe[ri][ci]);
  vector nn(dofe[ri][ci]);
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordcm[ri][ci]);
  
  nullv (vrhs);
  
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  for (i=0;i<intordcm[ri][ci];i++){
    xi=gp[i];  w1=w[i];
    
    jac_1d (jac,x,xi);
    
    bf_matrix (n,xi);
    
    c=Tm->volume_rhs2 (ipp,ri,ci);
    
    area = approx (xi,a);
    
    jac*=w1*area*c;
    
    for (ii=0;ii<n.n;ii++){
      nn[ii] = n[0][ii];
    }
    
    addmultv(vrhs,nn,jac);
    ipp++;
  }
}




/**
   function assembles resulting element volume right-hand side

   @param eid - element id
   @param lcid - load case id
   @param f - resulting volume right-hand side %vector of one element

   TKr, 16.5.2007
*/
void linbart::res_volume_rhs_vector (vector &f,long eid,long /*lcid*/)
{
  long i,*cn;
  vector lf;

  for (i=0;i<ntm;i++){
    cn = new long [dofe[i][i]];
    reallocv (dofe[i][i],lf);
    codnum (cn,i);
    volume_rhs_vector (i,eid,i,i,lf);
    locglob (f.a,lf.a,cn,dofe[i][i]);
    delete [] cn;
  }
}


/**
   function assembles resulting element volume right-hand side of the second type

   @param eid - element id
   @param lcid - load case id
   @param f - resulting volume right-hand side %vector of one element

   TKr, 08/06/2018
*/
void linbart::res_volume_rhs_vector2 (vector &f,long eid,long /*lcid*/)
{
  long i,*cn;
  vector lf;

  for (i=0;i<ntm;i++){
    cn = new long [dofe[i][i]];
    reallocv (dofe[i][i],lf);
    codnum (cn,i);
    volume_rhs_vector2 (i,eid,i,i,lf);
    locglob (f.a,lf.a,cn,dofe[i][i]);
    delete [] cn;
  }
}




/**
   function assembles resulting element capacity %matrix
   
   @param eid - element id
   @param cm - resulting capacity %matrix of one element

   JK, 6.1.2002
*/
void linbart::res_capacity_matrix (long eid,matrix &cm)
{
  long i,j,*rcn,*ccn;
  matrix lcm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      rcn = new long [dofe[i][j]];
      ccn = new long [dofe[i][j]];
      reallocm (dofe[i][j],dofe[i][j],lcm);
      capacity_matrix (eid,i,j,lcm);

      // diagonalization of capacity matrix on one element
      if (Tp->diagcap == 1){
	diagonalization (lcm);
      }
      
      codnum (rcn,i);
      codnum (ccn,j);
      mat_localize (cm,lcm,rcn,ccn);
      delete [] rcn;  delete [] ccn;
    }
  }
}

/**
   function assembles resulting element convection %vector

   @param f - resulting convection %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   
   JK, 6.1.2002
*/
void linbart::res_convection_vector (vector &f,long lcid,long eid,long leid)
{
  long *cn;
  vector lf;

  //  transi[lcid]==2 - element contains boundary with prescribed flux
  //  transi[lcid]==4 - element contains boundary with prescribed transmission and prescribed flux
  if ((Tt->elements[eid].transi[lcid]==2)||(Tt->elements[eid].transi[lcid]==4)){
    //  array for code numbers
    cn = new long [dofe[lcid][lcid]];
    //  array for nodal values of one matter
    reallocv (dofe[lcid][lcid],lf);

    convection_vector (lf,lcid,eid,leid,lcid,lcid);

    //  code numbers
    codnum (cn,lcid);
    //  localization to element vector
    locglob (f.a,lf.a,cn,dofe[lcid][lcid]);
    delete [] cn;
    
    cmulv (-1.0,f,f);
  }
}

/**
   function assembles resulting element transmission %vector

   @param f - resulting transmission %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - loaded element id
   
   JK, 6.1.2002
*/
void linbart::res_transmission_vector (vector &f,long lcid,long eid,long leid)
{
  long *cn,i;
  vector lf;
  
  cn = new long [dofe[lcid][lcid]];
  //  code numbers
  codnum (cn,lcid);
  reallocv (dofe[lcid][lcid],lf);
  
  //  transi[i]==3 - element contains boundary with prescribed transmission
  //  transi[i]==4 - element contains boundary with prescribed transmission and prescribed flux
  if ((Tt->elements[eid].transi[lcid]==3)||(Tt->elements[eid].transi[lcid]==4)){
    
    for (i=0;i<ntm;i++){    
      nullv (lf.a,lf.n);

      transmission_vector (lf,lcid,eid,leid,i);
      
      locglob (f.a,lf.a,cn,dofe[lcid][lcid]);
    }
  }
  
  delete [] cn;
}

/**
   function assembles resulting element source %vector

   @param sv - resulting source %vector of one element
   @param lcid - load case id
   @param eid - element id
   
   JK, 6.1.2002
*/
void linbart::res_quantity_source_vector (vector &sv,vector &nodval,long lcid,long eid)
{
  long *cn;
  vector lsv;

  cn = new long [dofe[lcid][lcid]];
  reallocv (dofe[lcid][lcid],lsv);
  codnum (cn,lcid);
  
  quantity_source_vector (lsv,nodval,eid,lcid,lcid);
  locglob (sv.a,lsv.a,cn,dofe[lcid][lcid]);
  
  delete [] cn;
}



/**
   function assembles resulting element internal fluxes vector

   @param eid - element id
   @param elemif - resulting internal fluxes %vector of one element
   
   JK, 6.1.2002
*/
void linbart::res_internal_fluxes (long eid,vector &elemif)
{
  long i;
  vector lif, tdnv, nodval;
  ivector cn;
  matrix cm;
  

  for (i=0; i<ntm; i++){
    reallocv(RSTCKIVEC(dofe[i][i], cn));
    reallocv(RSTCKVEC(dofe[i][i], lif));
    codnum(cn.a, i);
    internal_fluxes(i, eid, lif);
    locglob(elemif.a, lif.a, cn.a, dofe[i][i]);
  } 

  reallocv(RSTCKVEC(ndofe, tdnv));
  reallocv(RSTCKVEC(ndofe, lif));
  reallocm(RSTCKMAT(ndofe, ndofe, cm));
  
  nodalderivatives(eid, tdnv);
  res_capacity_matrix(eid, cm);
  mxv(cm, tdnv, lif);
  
  subv(lif, elemif, elemif);
  
  if (Tt->elements[eid].source==1){
    reallocv(RSTCKVEC(nne, nodval));
    nullv(nodval);
    nullv(lif);
    
    Tb->lc[0].sourcenodalvalues(eid, nodval);
    
    source_vector(0, eid, nodval, lif);    
    addv(elemif, lif, elemif);
  }
}


/**
   function computes nodal energies for one transported matter
   finite element with linear approximation functions
   
   @param lcid - number of load case
   @param eid - number of element
   @param nener - %vector of nodal energy
   
   JK, 8.4.2019
*/
void linbart::nodal_energies (long lcid,long eid,vector &nener)
{
  long i,ipp;
  double jac,area,xi;
  ivector nodes(nne);
  vector x(nne),w(intordkm[lcid][lcid]),gp(intordkm[lcid][lcid]),fl(ncomp),contr(dofe[lcid][lcid]),a(nne);
  matrix gm(ncomp,dofe[lcid][lcid]),d(ncomp,ncomp);
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[lcid][lcid]);
  
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[lcid][lcid];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];


  for (i=0;i<intordkm[lcid][lcid];i++){
    xi=gp[i];
    
    Tm->computenlfluxes (lcid,ipp);
    
    Tm->givefluxes (lcid,ipp,fl);

    grad_matrix (gm,x,gp[i],jac);
    
    mtxv (gm,fl,contr);
    
    //  area of cross section
    area = approx (xi,a);
    cmulv (area*jac*w[i],contr);
    
    addv (contr,nener,nener);
    
    ipp++;
  }
}


/**
   function assembles resulting nodal energies in nodes on element

   @param eid - element id
   @param elemne - resulting %vector of nodal energies
   
   JK, 8.4.2019
*/
void linbart::res_nodal_energy (long eid,vector &elemne,double dt)
{
  long i;
  vector nener, tdnv, nodval;
  ivector cn;
  matrix cm;
  

  for (i=0; i<ntm; i++){
    reallocv(RSTCKIVEC(dofe[i][i], cn));
    reallocv(RSTCKVEC(dofe[i][i], nener));
    codnum(cn.a, i);
    internal_fluxes(i, eid, nener);
    locglob(elemne.a, nener.a, cn.a, dofe[i][i]);
  } 
  
  cmulv (dt,elemne);
  
  reallocv(RSTCKVEC(ndofe, tdnv));
  reallocv(RSTCKVEC(ndofe, nener));
  reallocm(RSTCKMAT(ndofe, ndofe, cm));
  
  //nodalderivatives(eid, tdnv);
  //res_capacity_matrix(eid, cm);
  //mxv(cm, tdnv, nener);

  for (i=0; i<ntm; i++){
    reallocv(RSTCKIVEC(dofe[i][i], cn));
    reallocv(RSTCKVEC(dofe[i][i], nener));
    codnum(cn.a, i);
    accumul_energy (i,eid,nener);
    cmulv(-1.0,nener);
    locglob(elemne.a, nener.a, cn.a, dofe[i][i]);
  } 
  
  if (Tt->elements[eid].source==1){
    reallocv(RSTCKVEC(nne, nodval));
    nullv(nodval);
    nullv(nener);
    
    Tb->lc[0].sourcenodalvalues(eid, nodval);
    
    source_vector(0, eid, nodval, nener);    
    addv(elemne, nener, elemne);
  }
}



/**
   function assembles resulting element advection %matrix

   @param eid - element id
   @param lcid - load case id
   @param hm - resulting advection %matrix of one element

   JK, 2. 5. 2016
*/
void linbart::res_advection_matrix (long eid,long /*lcid*/,matrix &hm)
{
  long i,j,*rcn,*ccn;
  matrix lhm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      rcn = new long [dofe[i][j]];
      ccn = new long [dofe[i][j]];
      reallocm (dofe[i][j],dofe[i][j],lhm);
      advection_matrix (i,eid,i,j,lhm);
      codnum (rcn,i);
      codnum (ccn,j);
      mat_localize (hm,lhm,rcn,ccn);
      delete [] rcn;  delete [] ccn;
    }
  }
}


/**
   function computes element quantity integral
   the quantity is stored in nodes

   @param eid - element id
   @param nodval - %vector of quantity nodal values

   @retval f - element quantity integral

   TKr, 30.1.2004
*/
double linbart::total_integral (long eid,vector &nodval)
{
  long i;
  double area,xi,ww,jac,value,f;
  ivector nodes(nne);
  vector x(nne),w(2),gp(2),a(nne);

  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  
  gauss_points (gp.a,w.a,2);
  
  f = 0.0;

  for (i=0;i<2;i++){
    xi=gp[i];  ww=w[i];
    
    jac_1d (jac,x,xi);

    value = approx (xi,nodval);
    
    //  area of cross section
    area = approx (xi,a);

    jac*=area*ww*value;
    
    f = f + jac;
  }
  return(f);
}

/**
   function computes element quantity integral
   the quantity is stored in integration points

   @param eid - element id
   @param varid - id of variable required
   
   @retval f - element quantity integral

   JK, 3. 10. 2013
*/
double linbart::total_integral_ip (long eid,long varid)
{
  long i,ii;
  double area,xi,ww,jac,value,f;
  ivector nodes(nne);
  vector x(nne),w(2),gp(2),a(nne);

  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  
  gauss_points (gp.a,w.a,2);
  
  f = 0.0;
  
  long ri=0;
  long ci=0;
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];

  for (i=0;i<2;i++){
    xi=gp[i];  ww=w[i];
    
    jac_1d (jac,x,xi);

    value = Tm->ip[ii].eqother[varid];
    
    //  area of cross section
    area = approx (xi,a);

    jac*=area*ww*value;
    
    f = f + jac;
  }
  return f;
}

/**************************************************************///added
/**
   function computes contributions to the boundary flux from transmission %vector
   
   \int_{Gamma_3} N^T c_{tr} N dGamma * nodal_external_value

   @param tmv - %vector of boundary fluxes
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   
   TKr, 28.2.2004
*/
void linbart::boundary_flux (vector &tmv,long lcid,long eid,long leid,long ri,long ci)
{
  long ii;
  double xi,area;
  double new_nodval,tr;
  long nn;
  bocontypet bc[2];
  
  ivector nodes(nne);
  vector x(nne),trc(2),trcn(nne),nodval(2),av(dofe[ri][ci]),v(dofe[ri][ci]),a(nne);
  vector trr(2);
  matrix n(1,dofe[ri][ci]),km(dofe[ri][ci],dofe[ri][ci]);
  
  //  element nodes
  Tt->give_elemnodes (eid,nodes);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  
  //  indicators of boundary conditions
  Tb->lc[lcid].elemload[leid].give_bc (bc);
  //  nodal values of exterior variable
  Tb->lc[lcid].elemload[leid].give_nodval (lcid,nodval);
  //  transmission coefficients
  Tb->lc[lcid].elemload[leid].give_trc (Tp->time,Tb->lc[lcid].nodval,lcid,ci,trc);
  //  transmission/radiation coefficients
  Tb->lc[lcid].elemload[leid].give_trr (Tp->time,Tb->lc[lcid].nodval,trr);

  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];

  if (bc[0]==5 || bc[0]>10){
    nullm (km);
    xi=-1.0;
    
    bf_matrix (n,xi);

    nn = nodes[0];
    tr = 0.0;
    if(bc[0] == 90)
      tr = trr[0]/trc[0];
    
    Tm->transmission_flux(new_nodval,nodval[0],tr,ri,ci,nn,bc[0],ii);
   
    av[0]=new_nodval;  av[1]=0.0;
    
    //  area of cross section
    area = approx (xi,a);
    nnj (km.a,n.a,trc[0]*area,n.m,n.n);
    
    mxv (km,av,v);  addv (tmv,v,tmv);
  }
  if (bc[1]==5 || bc[1]>10){
    nullm (km);
    xi=1.0;
    
    bf_matrix (n,xi);

    nn = nodes[1];
    tr = 0.0;
    if(bc[0] == 90)
      tr = trr[1]/trc[1];
    
    Tm->transmission_flux(new_nodval,nodval[1],tr,ri,ci,nn,bc[1],ii);
   
    av[0]=0.0;  av[1]=new_nodval;
          
    //  area of cross section
    area = approx (xi,a);
    nnj (km.a,n.a,trc[1]*area,n.m,n.n);

    mxv (km,av,v);  addv (tmv,v,tmv);
  }
}



/**
   function assembles resulting element boundary flux %vector

   @param f - resulting boundary flux %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   
   TKr, 28.2.2004
*/
void linbart::res_boundary_flux (vector &f,long lcid,long eid,long leid)
{
  long *cn,i;
  vector lf;

  cn = new long [dofe[lcid][lcid]];
  codnum (cn,lcid);
  
  //coupled b.c.
  for (i=0;i<ntm;i++){    
    reallocv (dofe[lcid][lcid],lf);
    
    if ((Tt->elements[eid].transi[i]==3)||(Tt->elements[eid].transi[i]==4))
      boundary_flux (lf,i,eid,leid,lcid,i);
    
    locglob (f.a,lf.a,cn,dofe[lcid][lcid]);
  }
  
  delete [] cn;
}


/**
   function computes correct fluxes at integration points on element

   @param eid - element id
   
   TKr, 01/02/2010
*/
void linbart::intpointflux (long eid)
{
  long i,ii,jj,k,ipp;
  
  for (k=0;k<Tp->ntm;k++){
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
	
	ipp=Tt->elements[eid].ipp[ii][jj];
	
	for (i=0;i<intordkm[ii][jj];i++){	  
	  //  computation of correct fluxes
	  if (Tp->fluxcomp==1)
	    Tm->computenlfluxes (k,ipp);
	  
	  ipp++;
	}
	
	for (i=0;i<intordcm[ii][jj];i++){	  
	  //  computation of correct fluxes
	  if (Tp->fluxcomp==1)
	    Tm->computenlfluxes (k,ipp);
	  
	  ipp++;
	}
	if (Tp->savemode==1)  break;
      }
      if (Tp->savemode==1)  break;
    }
  }
}

/**
   function computes gradients in nodes of element

   @param eid - element id
   
   TKr, 01/02/2010
*/
void linbart::nod_grads_ip (long eid)
{
  long i,j,k,ipp;
  ivector ipnum(nne),nod(nne);
  vector grad(ncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  nodip_bar (ipp,intordkm[0][0],ipnum);

  //  node numbers of the element
  Tt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    for (k=0;k<Tp->ntm;k++){
      //  strains at the closest integration point
      Tm->givegrad (k,ipnum[i],grad);
      
      //  storage of strains to the node
      j=nod[i];
      Tt->nodes[j].storegrad (k,grad);
    }
  }
  
}


/**
   function computes fluxes in nodes of element

   @param eid - element id
   
   TKr, 01/02/2010
*/
void linbart::nod_fluxes_ip (long eid)
{
  long i,j,k,ipp;
  ivector ipnum(nne),nod(nne);
  vector flux(ncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  nodip_bar (ipp,intordkm[0][0],ipnum);
  
  //  node numbers of the element
  Tt->give_elemnodes (eid,nod);
  
  intpointflux (eid);

  for (i=0;i<nne;i++){
    for (k=0;k<Tp->ntm;k++){
      //  strains at the closest integration point
      Tm->givefluxes (k,ipnum[i],flux);
      
      //  storage of strains to the node
      j=nod[i];
      Tt->nodes[j].storeflux (k,flux);
    }
  }
  
}



/**
   function computes other values directly in nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index

   TKr, 03/02/2010
*/
void linbart::nod_others_comp (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp,ncompother,ii;
  vector h,r(ASTCKVEC(ndofe));
  ivector nod(ASTCKIVEC(nne));
  double other;
  
  Tt->give_elemnodes (eid,nod);  
  elemvalues(eid, r);
  ncompother = Tm->givencompother();
  reallocv (RSTCKVEC(ntm,h));

  ii = 0;
  for (i=0;i<nne;i++){
    if (Tp->savemode==0)
      ipp=Tt->elements[eid].ipp[ri][ci];
    if (Tp->savemode==1)
      ipp=Tt->elements[eid].ipp[0][0];
    
    for (j=0;j<ntm;j++){
      h[j] = r[ii];
      ii++;
    }
    
    for(k=0;k<ncompother;k++){
      other = Tm->givecompother(k,ipp,h.a);
      Tt->nodes[nod[i]].storeother (k,other);
    }
  }
}



/**
   Function returns transport (non-mechanical) quantities at nodes of element.
   The values of selected quantity is copied from the closest integration points 
   to element nodes.

   @param eid - element id
   @param nodval - %vector of nodal values
   @param ntq - type of non-mechanical quantity
   
   @return The function does not return anything.
   
   12/06/2012 TKr
   Modified by TKo, 10.10.2013
*/
void linbart::transq_nodval (long eid,vector &nodval,nonmechquant nmq)
{
  long i,ipid;
  ivector ipnum(nne);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipid=Tt->elements[eid].ipp[0][0];
  nodip_bar (ipid,intordkm[0][0],ipnum);

  for (i=0;i<nne;i++)
    {
      //copy transport (non-mechanical) quantity from closest int. point
      nodval[i] = Tm->givetransq(nmq, ipnum[i]);
    }
}



/**
   Function computes transport (non-mechanical) quantities in nodes of element.

   @param eid - element id
   @param nodval - %vector of nodal values of all required quantities, i.e., 
                   nodal value of i-th quantity in j-th node is given by nodval[i*ncnv+j] where ncnv 
                   is the number of calculated nodes on eid-th element.
   @param ncne - number of computed nodes on element (only first ncne of nodes is calculated)
   @param nq - number of required transport quantities
   @param qt - array of types of required non-mechanical, i.e. transport quantities
   
   @return The function does not return anything.
   
   Modified by TKo, 10.12.2013
*/
void linbart::transq_nodval_comp (long eid,vector &nodval, long ncne, long nq, nonmechquant *qt)
{
  long i,j,ipid;
  long ncompo;
  ivector ipnum(nne), enod(nne);
  intpointst ipb; // backup of the first integration point of element
  vector ipav;
    
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipid=Tt->elements[eid].ipp[0][0];

  // element nodes
  Tt->give_elemnodes(eid, enod);

  //  numbers of integration points closest to element nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodip_bar (ipid,intordkm[0][0],ipnum);

  // calculate nodal gradients - actually, only one rarely used material model depends on gradients
  // nod_grads_comp(eid);

  // store original content of the first integration point on element because 
  // it becomes working int. point for nodal values calculations on the given element
  ipb.copy(Tm->ip[ipid], ntm, 1);

  // The first integration point will be used for computation of nodal values temporarily
  // then the original content of the first integration point will be restored
  for (i=0;i<ncne;i++)
  {
    makerefv(ipav, Tm->ip[ipid].av, Tp->ntm);
    // store nodal values to the first (working) integration point
    nodalval (enod[i], ipav);
    // store gradients to the first (working) integration point
    // for(j=0; j<ntm; j++)   
    //   Tm->ip[ipid].storegrad(j, Tt->nodes[enod[i]].gradient[j]);
    
    ncompo = Tm->ip[ipnum[i]].ncompeqother;
    if (ncompo){
      // take eqother values for the given node from the closest integration point 
      Tm->storeeqother(ipid, 0, ncompo, Tm->ip[ipnum[i]].eqother);
    }
    ncompo = Tm->ip[ipnum[i]].ncompother;
    if (ncompo){
      // take eqother values for the given node from the closest integration point 
      Tm->storeother(ipid, 0, ncompo, Tm->ip[ipnum[i]].other);
    }
    
    Tm->mat_aux_values(ipid);

    //give calculated transport quantity from the first int. point
    for (j=0; j<nq; j++)
      nodval[j*ncne+i] = Tm->givetransq(qt[j], ipid);
  }
  // restore original integration point content of other/eqother arrays
  Tm->ip[ipid].copy(ipb, ntm, 0);
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
void linbart::transq_init_nodval (long eid,vector &nodval,nonmechquant nmq)
{
  long i,ipid;
  ivector ipnum(ASTCKIVEC(nne)), enod(ASTCKIVEC(nne));
  vector ipav;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipid=Tt->elements[eid].ipp[0][0];
  nodip_bar (ipid,intordkm[0][0],ipnum);
  // element nodes
  Tt->give_elemnodes(eid, enod);

  for (i=0;i<nne;i++)
  {
    // create reference vector to the int. point av array
    makerefv(ipav, Tm->ip[ipnum[i]].av, Tp->ntm);
    // store initial nodal values to the integration point
    initnodval2(enod[i], ipav);
    //copy transport (non-mechanical) quantity from closest int. point
    nodval[i] = Tm->givetransq(nmq, ipnum[i]);
  }
}



/**
   Function computes initial transport (non-mechanical) quantities in nodes of element.

   @param eid - element id
   @param nodval - %vector of nodal values of all required quantities, i.e., 
                   nodal value of i-th quantity in j-th node is given by nodval[i*ncnv+j] where ncnv 
                   is the number of calculated nodes on eid-th element.
   @param ncne - number of computed nodes on element (only first ncne of nodes is calculated)
   @param nq - number of required transport quantities
   @param qt - array of types of required non-mechanical, i.e. transport quantities
   
   @return The function does not return anything.
   
   Created by TKo, 5.6.2018
*/
void linbart::transq_init_nodval_comp (long eid,vector &nodval, long ncne, long nq, nonmechquant *qt)
{
  long i,j,ipid;
  long ncompo;
  ivector ipnum(ASTCKIVEC(nne)), enod(ASTCKIVEC(nne));
  intpointst ipb; // backup of the first integration point of element
  vector ipav;
    
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipid=Tt->elements[eid].ipp[0][0];

  // element nodes
  Tt->give_elemnodes(eid, enod);

  //  numbers of integration points closest to element nodes
  //  (function is from the file GEFEL/ordering.cpp)
  nodip_bar (ipid,intordkm[0][0],ipnum);

  // calculate nodal gradients - actually, only one rarely used material model depends on gradients
  // nod_grads_comp(eid);

  // store original content of the first integration point on element because 
  // it becomes working int. point for nodal values calculations on the given element
  ipb.copy(Tm->ip[ipid], ntm, 1);

  // The first integration point will be used for computation of nodal values temporarily
  // then the original content of the first integration point will be restored
  for (i=0;i<ncne;i++)
  {
    // create reference vector to the int. point av array
    makerefv(ipav, Tm->ip[ipid].av, Tp->ntm);
    // store initial nodal values to the first (working) integration point
    initnodval2 (enod[i], ipav);
    // store gradients to the first (working) integration point
    // for(j=0; j<ntm; j++)   
    //   Tm->ip[ipid].storegrad(j, Tt->nodes[enod[i]].gradient[j]);
    
    ncompo = Tm->ip[ipnum[i]].ncompeqother;
    if (ncompo){
      // take eqother values for the given node from the closest integration point 
      Tm->storeeqother(ipid, 0, ncompo, Tm->ip[ipnum[i]].eqother);
    }
    ncompo = Tm->ip[ipnum[i]].ncompother;
    if (ncompo){
      // take eqother values for the given node from the closest integration point 
      Tm->storeother(ipid, 0, ncompo, Tm->ip[ipnum[i]].other);
    }
    
    Tm->mat_aux_values(ipid);

    //give calculated transport quantity from the first int. point
    for (j=0; j<nq; j++)
      nodval[j*ncne+i] = Tm->givetransq(qt[j], ipid);
  }
  // restore original integration point content of other/eqother arrays
  Tm->ip[ipid].copy(ipb, ntm, 0);
}



/**
   Function assembles global coordinates of integration points.
   
   @param eid - element id
   @param ipp - integration point pointer
   @param ri - row index of the block of integration points /according to transported media/ (input)
   @param ci - column index of the block of integration points /according to transported media/ (input)
   @param coord - array containing coordinates of integration points (output)
   
   @retval Function returns coordinates of integration points in the %vector coord.
   @retval 0

   TKo, 12.2016
*/
long linbart::ipcoord (long eid, long ipp, long ri, long ci, vector &coord)
{
  long i, ii;
  double xi;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), w, gp;

  Tt->give_node_coord2d(x, y, eid);
  
  //  integration points for the conductivity matrix
  reallocv(RSTCKVEC(intordkm[ri][ci], gp));
  reallocv(RSTCKVEC(intordkm[ri][ci], w));
  gauss_points(gp.a, w.a, intordkm[ri][ci]);
	
  ii = Tt->elements[eid].ipp[ri][ci];
  for (i=0; i<intordkm[ri][ci]; i++)
  {
    if (ii == ipp)
    {
      xi=gp[i];
      coord[0] = approx(xi, x);
      coord[1] = approx(xi, y);
      coord[2] = 0.0;
      return 0;
    }
    ii++;
  }
	
  //  integration points for the capacity matrix
  reallocv(RSTCKVEC(intordcm[ri][ci], gp));
  reallocv(RSTCKVEC(intordcm[ri][ci], w));
  gauss_points (gp.a, w.a, intordcm[ri][ci]);
  
  for (i=0; i<intordcm[ri][ci]; i++)
  {
    if (ii == ipp)
    {
      xi = gp[i];
      coord[0] = approx(xi, x);
      coord[1] = approx(xi, y);
      coord[2] = 0.0;
      return 0;
    }
    ii++;
  }

  return 1;
}



/**
   Function computes natural coordinates of integration point.
   
   @param eid - element id (input)
   @param ipp - integration point pointer (input)
   @param ncoord - array containing coordinates of integration points (output)
   
   @retval coord - function returns coordinates of integration points in the %vector ncoord.
   @retval 0 - on success
   @retval 1 - integration point with ipp was not found

   TKo, 12.2016
*/
long linbart::ipncoord (long eid, long ipp, vector &ncoord)
{
  long i, ri, ci, ii;
  double xi;
  vector w, gp;

  for (ri=0;ri<ntm;ri++)
  {
    for (ci=0;ci<ntm;ci++)
    {
      //  integration points for the conductivity matrix
      reallocv(RSTCKVEC(intordkm[ri][ci], gp));
      reallocv(RSTCKVEC(intordkm[ri][ci], w));
      gauss_points(gp.a, w.a, intordkm[ri][ci]);
	
      ii = Tt->elements[eid].ipp[ri][ci];
      for (i=0; i<intordkm[ri][ci]; i++)
      {
        if (ii == ipp)
        {
          xi=gp[i];
          ncoord[0] = xi;
          ncoord[1] = 0.0;
          ncoord[2] = 0.0;
          return 0;
        }
        ii++;
      }
	
      //  integration points for the capacity matrix
      reallocv(RSTCKVEC(intordcm[ri][ci], gp));
      reallocv(RSTCKVEC(intordcm[ri][ci], w));
      gauss_points (gp.a, w.a, intordcm[ri][ci]);
  
      for (i=0; i<intordcm[ri][ci]; i++)
      {
        if (ii == ipp)
        {
          xi = gp[i];
          ncoord[0] = xi;
          ncoord[1] = 0.0;
          ncoord[2] = 0.0;
          return 0;
        }
        ii++;
      }
      if (Tp->savemode==1)  break;
    }
    if (Tp->savemode==1)  break;
  }

  return 1;
}



/**
   Function computes the total flux through a surface of an element
   
   @param lcid - load case id
   @param eid - element id (id in the list of all elements)
   @param beid - id in the list of selected elements
   @param flux - array of fluxes computed
   
   JK, 8. 3. 2018
*/
void linbart::surface_flux (long lcid,long eid,long beid,double *fluxes)
{
  long i,ipp,fid;
  double q,qn,area;
  ivector nodes(nne);
  bocontypet *bc;
  vector flux(ncomp),a(nne),n(ncomp);
  
  //  id of integration point for determination of the D matrix
  ipp=Tt->elements[eid].ipp[0][0];

  Tm->computenlfluxes (lcid,ipp);
  Tm->givefluxes (lcid,ipp,flux);

  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);


  //  indicators of boundary conditions
  bc = new bocontypet [nen];
  Tb->bf[lcid].elemload[beid].give_bc (bc);
  
  //  loop over end nodes
  for (i=0;i<nen;i++){
    
    if (bc[i]==1){
      //  function constructs the outer normal vector
      bar_normal_vectors (i,n);
      //  area of the i-th element surface
      area = a[i];
      //  normal density flux
      scprd (flux,n,qn);
      //  normal flux
      q = qn*area;
      //  id of flux
      fid = Tb->bf[lcid].elemload[beid].nvid[i];
      
      //  flux is stored to the appropriate position
      fluxes[fid]+=q;
    }
  }
  
}
