/*
  File:     linbartax.cpp
  Author:   Jaroslav Kruis, 31.3.2003
  Purpose:  onedimensional element with linear approximation functions
*/

#include "globalt.h"
#include "linbartax.h"
#include "genfile.h"
#include "globmatt.h"

linbartax::linbartax (void)
{
  long i;
  
  //  number of nodes on element
  nne=2;
  //  geometric problem dimension (1D)
  ncomp=1;
  //  number of end nodes
  nen=2;
  
  ned=2;
  nned=1;

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
    dofe[0][0]=2;  intordkm[0][0]=2;  intordcm[0][0]=2;  nip[0][0]=4;
    ndofe=2;  napfun=1;
    break;
  }
  case twomediacoup:{
    ordering[0][0]=1;  ordering[0][1]=3;
    ordering[1][0]=2;  ordering[1][1]=4;

    intordkm[0][0]=2;  intordkm[0][1]=2;  intordkm[1][0]=2;  intordkm[1][1]=2;
    intordcm[0][0]=2;  intordcm[0][1]=2;  intordcm[1][0]=2;  intordcm[1][1]=2;
    
    if (Tp->savemode==0){
      nip[0][0]=4;       nip[0][1]=4;       nip[1][0]=4;       nip[1][1]=4;
    }
    if (Tp->savemode==1){
      nip[0][0]=4;       nip[0][1]=0;       nip[1][0]=0;       nip[1][1]=0;
    }
    
    dofe[0][0]=2;  dofe[0][1]=2;  dofe[1][0]=2;  dofe[1][1]=2;
    ndofe=4;  napfun=2;
    break;
  }
  case threemediacoup:{
    ordering[0][0]=1;   ordering[0][1]=4;
    ordering[1][0]=2;   ordering[1][1]=5;
    ordering[2][0]=3;   ordering[2][1]=6;

    intordkm[0][0]=2;  intordkm[0][1]=2;  intordkm[0][2]=2;
    intordkm[1][0]=2;  intordkm[1][1]=2;  intordkm[1][2]=2;
    intordkm[2][0]=2;  intordkm[2][1]=2;  intordkm[2][2]=2;

    intordcm[0][0]=2;  intordcm[0][1]=2;  intordcm[0][2]=2;
    intordcm[1][0]=2;  intordcm[1][1]=2;  intordcm[1][2]=2;
    intordcm[2][0]=2;  intordcm[2][1]=2;  intordcm[2][2]=2;

    if (Tp->savemode==0){
      nip[0][0]=4;  nip[0][1]=4;  nip[0][2]=4;
      nip[1][0]=4;  nip[1][1]=4;  nip[1][2]=4;
      nip[2][0]=4;  nip[2][1]=4;  nip[2][2]=4;
    }
    if (Tp->savemode==1){
      nip[0][0]=4;  nip[0][1]=0;  nip[0][2]=0;
      nip[1][0]=0;  nip[1][1]=0;  nip[1][2]=0;
      nip[2][0]=0;  nip[2][1]=0;  nip[2][2]=0;
    }
    
    dofe[0][0]=2;  dofe[0][1]=2;  dofe[0][2]=2;
    dofe[1][0]=2;  dofe[1][1]=2;  dofe[1][2]=2;
    dofe[2][0]=2;  dofe[2][1]=2;  dofe[2][2]=2;

    ndofe=6;  napfun=3;
    break;
  }
  default:{
    print_err("unknown type of transported matter is required",__FILE__,__LINE__,__func__);
  }
  }
}

linbartax::~linbartax (void)
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

void linbartax::codnum (long *cn,long ri)
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
double linbartax::approx (double xi,vector &nodval)
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
void linbartax::intpointval (long eid)
{
  long i,k,ii,jj,ipp;
  double xi,val;
  vector r(ASTCKVEC(ndofe)),t(ASTCKVEC(nne)),gp,w;
  
  elemvalues(eid, r);
  
  
  for (k=0;k<ntm;k++){
    
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
	
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
void linbartax::intpointgrad (long eid)
{
  long i,ii,jj,k,ipp;
  double xi,jac;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), r(ASTCKVEC(ndofe)), t(ASTCKVEC(nne));
  vector gp, w, grad(ASTCKVEC(ncomp));
  matrix gm(ASTCKMAT(ncomp,nne));
  
  Tt->give_node_coord2d (x,y,eid);
  elemvalues(eid, r);

  for (k=0;k<Tp->ntm;k++){

    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }

    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){

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
void linbartax::intpointother (long eid)
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
  nodalotherval(nodes, r);
  
  for (k=0;k<ncompo;k++){
    
    for (i=0;i<nne;i++){
      t[i]=r[i*ncompo+k];
    }
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
	
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
   function assembles matrix of base functions
   
   @param n - %matrix of base functions
   @param xi - natural coordinate
   
   JK, 31.3.2002
*/
void linbartax::bf_matrix (matrix &n,double xi)
{
  nullm (n);
  bf_lin_1d (n.a,xi);
}

/**
   function assembles gradient of matrix of base functions
   
   @param gm - gradient %matrix
   @param x - array containing node coordinates
   @param xi - natural coordinate
   
   JK, 31.3.2002
*/
void linbartax::grad_matrix (matrix &gm,vector &x,double xi,double &jac)
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
   function computes conductivity matrix of 1D problems for one transported matter
   finite element with linear approximation functions
   
   @param lcid - number of load case
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting matrix
   @param km - conductivity %matrix
   
   JK, 31.3.2002
*/
void linbartax::conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,ii;
  double area,xinp,xi,ww,jac;
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
    
    xinp = approx (xi,x);
    //  area of cross section
    area = approx (xi,a);

    jac*=area*ww*xinp;
    //jac*=area*ww;//pro cylindr??!!
    
    //  contribution to the conductivity matrix of the element
    bdbj (km.a,gm.a,d.a,jac,gm.m,gm.n);
    
    //convective terms
    reallocm(1,ncomp,d);
    
    Tm->matcond2(d,ii,ri,ci);
    
    bf_matrix (n, xi);
    bdbjac(km, n, d, gm, jac);
    //bdbjac(km, n, d, gm, jac/xinp);//pro cylindr??!!
    
    ii++;
  }
  
  if ((Tt->elements[eid].transi[lcid]==3) || (Tt->elements[eid].transi[lcid]==4)){
    transmission_matrix (lcid,eid,ri,ci,km);
  }
}

/**
   function computes capacity matrix of 1D problems for one transported matter
   finite element with bilinear approximation functions
   
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting matrix
   @param cm - capacity %matrix
   
   JK, 31.3.2002
*/
void linbartax::capacity_matrix (long eid,long ri,long ci,matrix &cm)
{
  long i,ii;
  double jac,xi,ww,area,xinp,rho,c;
  ivector nodes(nne);
  vector x(nne),w(intordcm[ri][ci]),gp(intordcm[ri][ci]),a(nne);
  matrix n(1,dofe[ri][ci]);
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  array of cross section areas in nodes
  Tc->give_area (eid,nodes,a);
  //  node coordinates
  Tt->give_node_coord1d (x,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  //  density
  Tc->give_densitye (eid,rho);

  
  nullm (cm);

  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci]+intordkm[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0]+intordkm[0][0];

  for (i=0;i<intordcm[ri][ci];i++){
    xi=gp[i];  ww=w[i];
    
    bf_matrix (n,xi);
    
    c=Tm->capcoeff (ii,ri,ci);
    jac_1d (jac,x,xi);
    xinp = approx (xi,x);
    //  area of cross section
    area = approx (xi,a);
    jac*=area*ww*rho*c*xinp;;
    
    nnj (cm.a,n.a,jac,n.m,n.n);
    ii++;
  }
  
}


/**
   function computes source vector of one matter on one element
   
   \int_{Omega} N^T N d Omega . s
   
   @param sv - source %vector of one matter
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)
   
   JK, 31.3.2002
*/
void linbartax::quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci)
{
  long i;
  double jac,xi,ww,area,xinp;
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
    
    xinp = approx (xi,x);

    //  area of cross section
    area = approx (xi,a);
    jac*=ww*area*xinp;
    
    nnj (nm.a,n.a,jac,n.m,n.n);
  }
  mxv (nm,nodval,v);  addv (sv,v,sv);
  
}



/**
   function computes transmission complement to the conductivity matrix for one matter
   
   \int_{Gamma_3} N^T c_{tr} N dGamma
   
   @param lcid - number of load case
   @param eid - number of element
   @param ri,ci - row and column indices
   @param km - part of the conductivity %matrix
   
   JK, 31.3.2002
   TKr, 30.1.2004 - new added
*/
void linbartax::transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,leid,ii;
  double xi,area,xinp;
  double new_trc;//added
  long nn;//added 
  bocontypet bc[2];

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
      leid=i;
      break;
    }
  }
  
  //  indicators of boundary conditions
  Tb->lc[lcid].elemload[leid].give_bc (bc);
  //  transmission coefficients
  Tb->lc[lcid].elemload[leid].give_trc (Tp->time,Tb->lc[lcid].nodval,lcid,ci,trc);

  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];

 
  if (bc[0]==4 || bc[0]>10){
    xi=-1.0;
    
    bf_matrix (n,xi);
    
    /***********************************************///added
    nn = nodes[0];
    
    Tm->transmission_transcoeff(new_trc,trc[0],ri,ci,nn,bc[0],ii);
    
    trc[0]=new_trc;
    /***********************************************///added
    
    xinp = approx (xi,x);
    //  area of cross section
    area = approx (xi,a);
    
    nnj (km.a,n.a,trc[0]*area*xinp,n.m,n.n);
  }

  if (bc[1]==4 || bc[1]>10){
    xi=1.0;
    
    bf_matrix (n,xi);

    /***********************************************///added
    nn = nodes[1];
    
    Tm->transmission_transcoeff(new_trc,trc[1],ri,ci,nn,bc[1],ii);
    
    trc[1]=new_trc;
    /***********************************************///added

    xinp = approx (xi,x);
    //  area of cross section
    area = approx (xi,a);

    nnj (km.a,n.a,trc[1]*area*xinp,n.m,n.n);
  }
}

/**
   function computes contributions to the transmission vector
   
   \int_{Gamma_3} N^T c_{tr} N dGamma * nodal_external_value
   
   @param tmv - transmission %vector of one matter
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)
   
   JK, 5.10.2001
   TKr, 30.1.2002 - new added
*/
void linbartax::transmission_vector (vector &tmv,long lcid,long eid,long leid,long ri,long ci)
{
  long ii;
  double xi,area,xinp;
  double new_nodval,tr;//added
  long nn;//added
  bocontypet bc[2];
  
  ivector nodes(nne);
  vector x(nne),trc(2),trcn(nne),nodval(2),av(dofe[ri][ci]),v(dofe[ri][ci]),a(nne);
  vector trr(2);//added
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
  //  transmission coefficients
  Tb->lc[lcid].elemload[leid].give_trc (Tp->time,Tb->lc[lcid].nodval,lcid,ci,trc);
  //  transmission/radiation coefficients
  Tb->lc[lcid].elemload[leid].give_trr (Tp->time,Tb->lc[lcid].nodval,trr);


  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];

  if (bc[0]==4 || bc[0]>10){
    nullm (km);
    xi=-1.0;
    
    bf_matrix (n,xi);

    /***********************************************///added
    nn = nodes[0];
    tr = 0.0;
    if(bc[0] == 90)
      tr = trr[0]/trc[0];
    
    Tm->transmission_nodval(new_nodval,nodval[0],tr,ri,ci,nn,bc[0],ii);
   
    av[0]=new_nodval;  av[1]=0.0;
    /***********************************************/

    xinp = approx (xi,x);
    //  area of cross section
    area = approx (xi,a);

    nnj (km.a,n.a,trc[0]*area*xinp,n.m,n.n);
    
    mxv (km,av,v);  addv (tmv,v,tmv);
  }
  if (bc[1]==4 || bc[1]>10){
    nullm (km);
    xi=1.0;
    
    bf_matrix (n,xi);

    /***********************************************///added
    nn = nodes[1];
    tr = 0.0;
    if(bc[0] == 90)
      tr = trr[1]/trc[1];
    
    Tm->transmission_nodval(new_nodval,nodval[1],tr,ri,ci,nn,bc[1],ii);
   
    av[0]=0.0;  av[1]=new_nodval;
    /***********************************************/

    xinp = approx (xi,x);
    //  area of cross section
    area = approx (xi,a);
          
    nnj (km.a,n.a,trc[1]*area*xinp,n.m,n.n);

    mxv (km,av,v);  addv (tmv,v,tmv);
  }
}



/**
   function computes contribution to the convection vector

   \int_{Gamma_2} N^T N dGamma * nodal_flux_values
   
   @param f - convection %vector of one matter
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)

   JK, 8.10.2001
*/
void linbartax::convection_vector (vector &f,long lcid,long eid,long leid,long ri,long ci)
{
  bocontypet bc[2];
  double xi,area,xinp;
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


  if (bc[0]==2 || bc[0]==3 || bc[0]==4){
    nullm (km);
    xi=-1.0;

    bf_matrix (n,xi);

    xinp = approx (xi,x);
    //  area of cross section
    area = approx (xi,a);

    nnj (km.a,n.a,area*xinp,n.m,n.n);

    av[0]=nodval[0];  av[1]=0.0;
    mxv (km,av,v);  addv (f,v,f);
  }
  if (bc[1]==2 || bc[1]==3 || bc[1]==4){
    nullm (km);
    xi=1.0;
    bf_matrix (n,xi);

    xinp = approx (xi,x);
    //  area of cross section
    area = approx (xi,a);

    nnj (km.a,n.a,area*xinp,n.m,n.n);

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
void linbartax::internal_fluxes (long lcid,long eid,vector &ifl)
{
  long i,ipp;
  double jac,area,xinp,xi;
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
    xinp = approx (xi,x);

    Tm->computenlfluxes (lcid,ipp);
    
    Tm->givefluxes (lcid,ipp,fl);

    grad_matrix (gm,x,gp[i],jac);
    
    mtxv (gm,fl,contr);
    
    //  area of cross section
    area = approx (xi,a);
    cmulv (area*jac*w[i]*xinp,contr);
    
    addv (contr,ifl,ifl);
    
    ipp++;
  }
}

/**
   function assembles resulting element conductivity matrix

   @param eid - element id
   @param lcid - load case id
   @param km - resulting conductivity %matrix of one element

   JK, 6.1.2002
*/
void linbartax::res_conductivity_matrix (long eid,long /*lcid*/,matrix &km)
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
   function assembles resulting element capacity matrix
   
   @param eid - element id
   @param cm - resulting capacity %matrix of one element

   JK, 6.1.2002
*/
void linbartax::res_capacity_matrix (long eid,matrix &cm)
{
  long i,j,*rcn,*ccn;
  matrix lcm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      rcn = new long [dofe[i][j]];
      ccn = new long [dofe[i][j]];
      reallocm (dofe[i][j],dofe[i][j],lcm);
      capacity_matrix (eid,i,j,lcm);

      ////////////////////////////////////////////////////
      // diagonalization of capacity matrix on one element
      if (Tp->diagcap == 1){
	diagonalization (lcm);
      }
      ////////////////////////////////////////////////////

      codnum (rcn,i);
      codnum (ccn,j);
      mat_localize (cm,lcm,rcn,ccn);
      delete [] rcn;  delete [] ccn;
    }
  }
}

/**
   function assembles resulting element convection vector

   @param f - resulting convection %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   
   JK, 6.1.2002
*/
void linbartax::res_convection_vector (vector &f,long lcid,long eid,long leid)
{
  long *cn;
  vector lf;

  cn = new long [dofe[lcid][lcid]];
  reallocv (dofe[lcid][lcid],lf);
  convection_vector (lf,lcid,eid,leid,lcid,lcid);
  codnum (cn,lcid);
  locglob (f.a,lf.a,cn,dofe[lcid][lcid]);
  delete [] cn;

  cmulv (-1.0,f,f);
  
}

/**
   function assembles resulting element transmission vector

   @param f - resulting transmission %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   
   JK, 6.1.2002
*/
void linbartax::res_transmission_vector (vector &f,long lcid,long eid,long leid)
{
  long *cn,i;
  vector lf;

  cn = new long [dofe[lcid][lcid]];
  codnum (cn,lcid);

  //coupled b.c.
  for (i=0;i<ntm;i++){    
    reallocv (dofe[lcid][lcid],lf);

    if ((Tt->elements[eid].transi[i]==3) || (Tt->elements[eid].transi[i]==4))
      transmission_vector (lf,i,eid,leid,lcid,i);
    
    locglob (f.a,lf.a,cn,dofe[lcid][lcid]);
  }
  
  delete [] cn;
 
}

/**
   function assembles resulting element source vector

   @param sv - resulting source %vector of one element
   @param lcid - load case id
   @param eid - element id
   
   JK, 6.1.2002
*/
void linbartax::res_quantity_source_vector (vector &sv,vector &nodval,long lcid,long eid)
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
void linbartax::res_internal_fluxes (long eid, vector &elemif)
{
  long i;
  ivector cn;
  vector lif,tdnv;
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
  
  subv(elemif, lif, elemif);
}


/**
   function computes element quantity integral
   
   @param eid - element id
   @param nodval - %vector of quantity nodal values

   @retval f - element quantity integral

   TKr, 30.1.2004
*/
double linbartax::total_integral(long eid,vector &nodval)
{
  long i;
  double area,xinp,xi,ww,jac,value,f;
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
    
    xinp = approx (xi,x);
    value = approx (xi,nodval);
    //  area of cross section
    area = approx (xi,a);
    
    jac*=area*ww*value*xinp;
    
    f = f + jac;
  }
  return(f);
}

/**************************************************************///added
/**
   function computes contributions to the boundary flux from transmission vector
   
   \int_{Gamma_3} N^T c_{tr} N dGamma * nodal_external_value

   @param tmv - %vector of boundary fluxes
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   @param ri,ci - row and column indices of the computed block in the resulting matrix
   
   TKr, 28.2.2004
*/
void linbartax::boundary_flux (vector &tmv,long lcid,long eid,long leid,long ri,long ci)
{
  long ii;
  double xi,area,xinp;
  double new_nodval,tr;
  long nn;
  bocontypet bc[2];
  
  ivector nodes(nne);
  vector x(nne),trc(2),trcn(nne),nodval(2),av(dofe[ri][ci]),v(dofe[ri][ci]),a(nne);
  vector trr(2);
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
  //  transmission coefficients
  Tb->lc[lcid].elemload[leid].give_trc (Tp->time,Tb->lc[lcid].nodval,lcid,ci,trc);
  //  transmission/radiation coefficients
  Tb->lc[lcid].elemload[leid].give_trr (Tp->time,Tb->lc[lcid].nodval,trr);

  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];

  if (bc[0]==4 || bc[0]>10){
    nullm (km);
    xi=-1.0;
    
    bf_matrix (n,xi);

    nn = nodes[0];
    tr = 0.0;
    if(bc[0] == 90)
      tr = trr[0]/trc[0];
    
    Tm->transmission_flux(new_nodval,nodval[0],tr,ri,ci,nn,bc[0],ii);
   
    av[0]=new_nodval;  av[1]=0.0;
    
    xinp = approx (xi,x);
    //  area of cross section
    area = approx (xi,a);

    nnj (km.a,n.a,trc[0]*area*xinp,n.m,n.n);
    
    mxv (km,av,v);  addv (tmv,v,tmv);
  }
  if (bc[1]==4 || bc[1]>10){
    nullm (km);
    xi=1.0;
    
    bf_matrix (n,xi);

    nn = nodes[1];
    tr = 0.0;
    if(bc[0] == 90)
      tr = trr[1]/trc[1];
    
    Tm->transmission_flux(new_nodval,nodval[1],tr,ri,ci,nn,bc[1],ii);
   
    av[0]=0.0;  av[1]=new_nodval;

    xinp = approx (xi,x);
    //  area of cross section
    area = approx (xi,a);

    nnj (km.a,n.a,trc[1]*area*xinp,n.m,n.n);

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
void linbartax::res_boundary_flux (vector &f,long lcid,long eid,long leid)
{
  long *cn,i;
  vector lf;

  cn = new long [dofe[lcid][lcid]];
  codnum (cn,lcid);
  
  //coupled b.c.
  for (i=0;i<ntm;i++){    
    reallocv (RSTCKVEC(dofe[lcid][lcid],lf));
    
    if ((Tt->elements[eid].transi[i]==3) || (Tt->elements[eid].transi[i]==4))
      boundary_flux (lf,i,eid,leid,lcid,i);
    
    locglob (f.a,lf.a,cn,dofe[lcid][lcid]);
  }
  
  delete [] cn;
}

/**
   function computes other values in nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index

   10.5.2002
*/
void linbartax::nod_others (long /*lcid*/,long eid,long ri,long ci)
{
  long i, j, k, ipp, ncompother, ii;
  vector other, h, r(ASTCKVEC(ndofe));
  ivector nod(ASTCKIVEC(nne));
  
  Tt->give_elemnodes (eid,nod);  
  elemvalues(eid, r);
  ncompother = Tm->givencompother();
  reallocv (RSTCKVEC(ncompother,other));
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
    
    for(k=0;k<ncompother;k++)
      other[k] = Tm->givecompother(k,ipp,h.a);
  }
}



/**
   Function computes global coordinates of integration point.
   
   @param eid - element id (input)
   @param ipp - integration point pointer (input)
   @param ri - row index of the block of integration points /according to transported media/ (input)
   @param ci - column index of the block of integration points /according to transported media/ (input)
   @param ncoord - %vector containing coordinates of integration points (output)
   
   @retval coord - function returns coordinates of integration points in the %vector coord.
   @retval 0 - on success
   @retval 1 - integration point with ipp was not found

   TKo, 12.2016
*/
long linbartax::ipcoord (long eid, long ipp, long ri, long ci, vector &coord)
{
  long i, ii;
  double xi;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), w, gp;

  Tt->give_node_coord2d(x, y, eid);
  
  //  integration points for the conductivity matrix
  reallocv(RSTCKVEC(intordkm[ri][ci], gp));
  reallocv(RSTCKVEC(intordkm[ri][ci], w));
  gauss_points(gp.a, w.a, intordkm[ri][ci]);
	
  ii=Tt->elements[eid].ipp[ri][ci];
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
   
   @retval ncoord - function returns coordinates of integration points in the %vector ncoord.
   @retval 0 - on success
   @retval 1 - integration point with ipp was not found

   TKo, 12.2016
*/
long linbartax::ipncoord (long eid, long ipp, vector &ncoord)
{
  long i, ii, ri, ci;
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
	
      ii=Tt->elements[eid].ipp[ri][ci];
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
   function computes contributions to the right-hand side - volume integral
      
   \int_{Omega} B^T D dOmega
   
   @param tmv - transmission %vector of one matter
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)
   
   TKr, 19/10/2023
*/
void linbartax::volume_rhs_vector (long /*lcid*/,long eid,long ri,long ci,vector &vrhs)
{
  long i,ii;
  double area,xi,xinp,ww,jac;
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
    
    xinp = approx (xi,x);
    //  area of cross section
    area = approx (xi,a);

    jac*=area*ww*xinp;
    
    //  contribution to the volume_rhs integral of the element
    nnjac (km,gm,d,jac);
    
    ii++;
  }

  for (i=0;i<vrhs.n;i++){
    vrhs[i] = km[i][0];
  }

}



/**
   function assembles resulting element volume right-hand side

   @param eid - element id
   @param lcid - load case id
   @param f - resulting volume right-hand side %vector of one element

   TKr, 19/10/2023
*/
void linbartax::res_volume_rhs_vector (vector &f,long eid,long /*lcid*/)
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
   function computes contributions to the right-hand side - volume integral of the second type
      
   \int_{Omega} N^T const dOmega
   
   @param vrhs - volume right-hand side %vector of one matter
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)
   
   TKr, 19/10/2023
*/
void linbartax::volume_rhs_vector2 (long /*lcid*/,long eid,long ri,long ci,vector &vrhs)
{
  long i,ii,ipp;
  double jac,xi,xinp,w1,area,c;
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

    xinp = approx (xi,x);
    area = approx (xi,a);
    
    jac*=w1*area*c*xinp;
    
    for (ii=0;ii<n.n;ii++){
      nn[ii] = n[0][ii];
    }
    
    addmultv(vrhs,nn,jac);
    ipp++;
  }
}



/**
   function assembles resulting element volume right-hand side of the second type

   @param eid - element id
   @param lcid - load case id
   @param f - resulting volume right-hand side %vector of one element

   TKr, 19/10/2023
*/
void linbartax::res_volume_rhs_vector2 (vector &f,long eid,long /*lcid*/)
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

