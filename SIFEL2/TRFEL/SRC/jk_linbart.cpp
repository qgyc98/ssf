/*
  File:     linbart.cpp
  Author:   Jaroslav Kruis, 31.3.2003
  Purpose:  onedimensional element with linear approximation functions
*/

#include "globalt.h"
#include "linbart.h"
#include "genfile.h"
#include "globmatt.h"

linbart::linbart (void)
{
  long i;
  
  nne=2;  ned=2;  ncomp=1; nned=1;
  
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
    fprintf (stderr,"\n\n unknown type of transported matter is required");
    fprintf (stderr,"\n in function linbart::linbart (file %s, line %d).\n",__FILE__,__LINE__);
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
  ivector cn(ndofe);
  vector r(ndofe),t(nne),gp,w;
  
  Tt->give_code_numbers (eid,cn.a);
  nodalvalues (0,r.a,cn.a,ndofe);
  
  
  for (k=0;k<ntm;k++){
    
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
	
	allocv (intordkm[ii][jj],gp);
	allocv (intordkm[ii][jj],w);
	gauss_points (gp.a,w.a,intordkm[ii][jj]);
	
	ipp=Tt->elements[eid].ipp[ii][jj];
	for (i=0;i<intordkm[ii][jj];i++){
	  xi=gp[i];
	  
	  val = approx (xi,t);
	  Tm->ip[ipp].av[k]=val;
	  ipp++;
	}
	
	destrv (gp);  destrv (w);
	
	allocv (intordcm[ii][jj],gp);
	allocv (intordcm[ii][jj],w);
	gauss_points (gp.a,w.a,intordcm[ii][jj]);
	
	for (i=0;i<intordcm[ii][jj];i++){
	  xi=gp[i];
	  
	  val = approx (xi,t);
	  Tm->ip[ipp].av[k]=val;
	  ipp++;
	}
	
	destrv (gp);  destrv (w);
	
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
  ivector cn(ndofe);
  vector x(nne),y(nne),r(ndofe),t(nne),gp,w,grad(ncomp);
  matrix gm(ncomp,nne);
  
  Tt->give_node_coord2d (x,y,eid);
  Tt->give_code_numbers (eid,cn.a);
  nodalvalues (0,r.a,cn.a,ndofe);

  for (k=0;k<Tp->ntm;k++){

    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }

    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){

	allocv (intordkm[ii][jj],gp);
	allocv (intordkm[ii][jj],w);
	gauss_points (gp.a,w.a,intordkm[ii][jj]);

	ipp=Tt->elements[eid].ipp[ii][jj];
	for (i=0;i<intordkm[ii][jj];i++){
	  xi=gp[i];
	  
	  //  matrix of gradients
	  grad_matrix (gm,x,xi,jac);
	  mxv (gm,t,grad);
	  Tm->storegrad (k,ipp,0,grad);
	  ipp++;
	}
	destrv (gp);  destrv (w);
	
	allocv (intordcm[ii][jj],gp);
	allocv (intordcm[ii][jj],w);
	gauss_points (gp.a,w.a,intordcm[ii][jj]);
	
	for (i=0;i<intordcm[ii][jj];i++){
	  xi=gp[i];
	  
	  //  matrix of gradients
	  grad_matrix (gm,x,xi,jac);
	  mxv (gm,t,grad);
	  Tm->storegrad (k,ipp,0,grad);
	  ipp++;
	}
	
	destrv (gp);  destrv (w);
	
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
void linbart::bf_matrix (matrix &n,double xi)
{
  fillm (0.0,n);
  bf_lin_1d (n.a,xi);
}

/**
   function assembles gradient of matrix of base functions
   
   @param gm - gradient %matrix
   @param x - array containing node coordinates
   @param xi - natural coordinate
   
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
   function computes conductivity matrix of 1D problems for one transported matter
   finite element with linear approximation functions
   
   @param lcid - number of load case
   @param eid - number of element
   @param ri,ci - row and column indeces of the computed block in the resulting matrix
   @param km - conductivity %matrix
   
   JK, 31.3.2002
*/
void linbart::conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,ii;
  double area,xi,ww,jac;
  ivector nodes(nne);
  vector x(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]);
  matrix gm(ncomp,dofe[ri][ci]),d;
  matrix n(1,dofe[ri][ci]);

  Tt->give_elemnodes (eid,nodes);
  Tc->give_area (eid,area);
  Tt->give_node_coord1d (x,eid);
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  fillm (0.0,km);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];
  
  for (i=0;i<intordkm[ri][ci];i++){
    xi=gp[i];  ww=w[i];
    
    //  matrix of gradients
    grad_matrix (gm,x,xi,jac);
    
    //  matrix of conductivity of the material
    allocm(ncomp,ncomp,d);
    Tm->matcond (d,ii,ri,ci,ncomp);
    
    jac*=area*ww;
    
    //  contribution to the conductivity matrix of the element
    bdbj (km.a,gm.a,d.a,jac,gm.m,gm.n);
    destrm(d);

    //convective terms
    if (Tp->tmatt == threemediacoup)
      {
	allocm(1,ncomp,d);

	Tm->matcond2(d,ii,ri,ci,ncomp);
	
	bf_matrix (n, xi);
	bdbjac(km, n, d, gm, jac);
	destrm(d);
      }
    
    ii++;
  }
  
  /*  //old
      if (Tt->elements[eid].transi[lcid]==1){
      if (ri == ci)
      }
  */
  
  //coupled b.c.//new
  if (Tt->elements[eid].transi[ci]==1){
    transmission_matrix (ci,eid,ri,ci,km);
  }
}

/**
   function computes capacity matrix of 1D problems for one transported matter
   finite element with bilinear approximation functions
   
   @param eid - number of element
   @param ri,ci - row and column indeces of the computed block in the resulting matrix
   @param cm - capacity %matrix
   
   JK, 31.3.2002
*/
void linbart::capacity_matrix (long eid,long ri,long ci,matrix &cm)
{
  long i,ii;
  double jac,xi,ww,area,rho,c;
  ivector nodes(nne);
  vector x(nne),w(intordcm[ri][ci]),gp(intordcm[ri][ci]);
  matrix n(1,dofe[ri][ci]);
  
  Tt->give_elemnodes (eid,nodes);
  Tc->give_area (eid,area);
  Tc->give_densitye (eid,rho);
  Tt->give_node_coord1d (x,eid);
  gauss_points (gp.a,w.a,intordcm[ri][ci]);
  
  fillm (0.0,cm);

  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci]+intordkm[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0]+intordkm[0][0];

  for (i=0;i<intordcm[ri][ci];i++){
    xi=gp[i];  ww=w[i];
    
    bf_matrix (n,xi);
    
    c=Tm->capcoeff (ii,ri,ci);
    jac_1d (jac,x,xi);
    jac*=area*ww*rho*c;
    
    nnj (cm.a,n.a,jac,n.m,n.n);
    ii++;
  }
  
}


/**
   function computes source vector of one matter on one element
   
   \int_{Omega} N^T N d Omega . s
   
   @param sv - source %vector of one matter
   @param nodval - array of nodal values
   @param eid - element id
   @param sourl - type of source definition (node/element)
   @param ri,ci - row and column indices of the block (ri must be equal to ci)
   
   JK, 31.3.2002
*/
void linbart::quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci)
{
  long i,j,k,ii;
  double s,jac,xi,ww,area;
  vector x(nne),w(intordcm[ri][ci]),gp(intordcm[ri][ci]),v(dofe[ri][ci]);
  matrix n(1,dofe[ri][ci]),nm(dofe[ri][ci],dofe[ri][ci]);
  
  Tc->give_area (eid,area);
  Tt->give_node_coord1d (x,eid);
  gauss_points (gp.a,w.a,intordcm[ri][ci]);
  
  fillm (0.0,nm);
  
  for (i=0;i<intordcm[ri][ci];i++){
    xi=gp[i];  ww=w[i];
    
    bf_matrix (n,xi);
    
    jac_1d (jac,x,xi);
    
    jac*=ww*area;
    
    nnj (nm.a,n.a,jac,n.m,n.n);
  }
  mxv (nm,nodval,v);  addv (sv,v,sv);
  
}



/**
   function computes transmission complement to the conductivity matrix for one matter
   
   \int_{Gamma_3} N^T c_{tr} N dGamma
   
   @param lcid - number of load case
   @param eid - number of element
   @param ri,ci - row and column indeces
   @param km - part of the conductivity %matrix
   
   JK, 31.3.2002
   TKr, 30.1.2004 - new added
*/
void linbart::transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,neleml,*bc,ii;
  double xi,area;
  double new_trc;//added
  long nn;//added 

  ivector nodes(nne);
  vector x(nne),trc(2);
  matrix n(1,dofe[ri][ci]);
  
  bc = new long [2];
  
  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord1d (x,eid);
  Tc->give_area (eid,area);
  
  for (i=0;i<Tb->lc[lcid].neb;i++){
    if (Tb->lc[lcid].elemload[i].eid==eid){
      neleml=i;
      break;
    }
  }
  
  for (i=0;i<2;i++){
    bc[i]=Tb->lc[lcid].elemload[neleml].bc[i];
  }
  for (i=0;i<2;i++){
    trc[i]=Tb->lc[lcid].elemload[neleml].trc[i];
  }
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];

 
  if (bc[0]>=30){
    xi=-1.0;
    
    bf_matrix (n,xi);
    
    /***********************************************///added
    nn = nodes[0];
    
    Tm->transmission_transcoeff(new_trc,trc[0],ri,ci,nn,bc[0],ii);
    //printf("new_trc   = %e\n",new_trc);

    trc[0]=new_trc;
    /***********************************************///added
    
    nnj (km.a,n.a,trc[0]*area,n.m,n.n);
  }

  if (bc[1]>=30){
    xi=1.0;
    
    bf_matrix (n,xi);

    /***********************************************///added
    nn = nodes[1];
    
    Tm->transmission_transcoeff(new_trc,trc[1],ri,ci,nn,bc[1],ii);
    //printf("new_trc   = %e\n",new_trc);
    
    trc[1]=new_trc;
    /***********************************************///added

    nnj (km.a,n.a,trc[1]*area,n.m,n.n);
  }
  
  delete [] bc;
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
void linbart::transmission_vector (vector &tmv,long lcid,long eid,long leid,long ri,long ci)
{
  long i,*bc,ii;
  double xi,area;
  double new_nodval,tr;//added
  long nn;//added
  
  ivector nodes(nne);
  vector x(nne),trc(2),trcn(nne),nodval(2),av(dofe[ri][ci]),v(dofe[ri][ci]);
  vector trr(2);//added
  matrix n(1,dofe[ri][ci]),km(dofe[ri][ci],dofe[ri][ci]);
  
  bc = new long [2];

  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord1d (x,eid);
  Tc->give_area (eid,area);
  
  for (i=0;i<2;i++){
    bc[i]=Tb->lc[lcid].elemload[leid].bc[i];
  }
  for (i=0;i<2;i++){
    trc[i]=Tb->lc[lcid].elemload[leid].trc[i];
  }
  /*  for (i=0;i<2;i++){
     //nodval[i]=Tb->lc[lcid].elemload[leid].nodval[i];
     }
  */
  Tb->lc[lcid].elemload[leid].getval(lcid,nodval);

  for (i=0;i<2;i++){//added
    trr[i]=Tb->lc[lcid].elemload[leid].trr[i];
  }  

  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];

  if (bc[0]>=30){
    fillm (0.0,km);
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
    
    nnj (km.a,n.a,trc[0]*area,n.m,n.n);
    
    mxv (km,av,v);  addv (tmv,v,tmv);
  }
  if (bc[1]>=30){
    fillm (0.0,km);
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
          
    nnj (km.a,n.a,trc[1]*area,n.m,n.n);

    mxv (km,av,v);  addv (tmv,v,tmv);
  }
  
  delete [] bc;
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
void linbart::convection_vector (vector &f,long lcid,long eid,long leid,long ri,long ci)
{
  long i,*bc;
  double xi,area;
  ivector nodes(nne);
  vector x(nne),nodval(2),av(dofe[ri][ci]),v(dofe[ri][ci]);
  matrix n(1,dofe[ri][ci]),km(dofe[ri][ci],dofe[ri][ci]);

  bc = new long [2];

  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord1d (x,eid);
  Tc->give_area (eid,area);

  for (i=0;i<2;i++){
    bc[i]=Tb->lc[lcid].elemload[leid].bc[i];
  }
  /*   for (i=0;i<2;i++){
       nodval[i]=Tb->lc[lcid].elemload[leid].nodval[i];
       }
  */
  Tb->lc[lcid].elemload[leid].getval(lcid,nodval);

  if (bc[0]==2){
    fillm (0.0,km);
    xi=-1.0;

    bf_matrix (n,xi);
    nnj (km.a,n.a,area,n.m,n.n);

    av[0]=nodval[0];  av[1]=0.0;
    mxv (km,av,v);  addv (f,v,f);
  }
  if (bc[1]==2){
    fillm (0.0,km);
    xi=1.0;
    bf_matrix (n,xi);
    nnj (km.a,n.a,area,n.m,n.n);
    av[0]=0.0;  av[1]=nodval[1];
    mxv (km,av,v);  addv (f,v,f);
  }
  
  delete [] bc;
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
  double jac,area;
  vector x(nne),w(intordkm[lcid][lcid]),gp(intordkm[lcid][lcid]),fl(ncomp),contr(dofe[lcid][lcid]);
  matrix gm(ncomp,dofe[lcid][lcid]),d(ncomp,ncomp);
  
  Tt->give_node_coord1d (x,eid);
  Tc->give_area (eid,area);
  gauss_points (gp.a,w.a,intordkm[lcid][lcid]);
  
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[lcid][lcid];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];


  for (i=0;i<intordkm[lcid][lcid];i++){
    
    Tm->computenlfluxes (lcid,ipp);
    
    Tm->givefluxes (lcid,ipp,0,fl);

    grad_matrix (gm,x,gp[i],jac);
    
    mtxv (gm,fl,contr);
    
    cmulv (area*jac*w[i],contr);
    
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
void linbart::res_conductivity_matrix (long eid,long lcid,matrix &km)
{
  long i,j,*rcn,*ccn;
  matrix lkm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      rcn = new long [dofe[i][j]];
      ccn = new long [dofe[i][j]];
      allocm (dofe[i][j],dofe[i][j],lkm);
      conductivity_matrix (i,eid,i,j,lkm);
      codnum (rcn,i);
      codnum (ccn,j);
      mat_localize (km,lkm,rcn,ccn);
      destrm (lkm);
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
void linbart::res_capacity_matrix (long eid,matrix &cm)
{
  long i,j,*rcn,*ccn,k,l;
  double s;
  matrix lcm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      rcn = new long [dofe[i][j]];
      ccn = new long [dofe[i][j]];
      allocm (dofe[i][j],dofe[i][j],lcm);
      capacity_matrix (eid,i,j,lcm);

      ////////////////////////////////////////////////////
      // diagonalization of capacity matrix on one element
      if (Tp->diagcap == 1){
	for (k=0;k<dofe[i][j];k++){
	  s=0.0;
	  for (l=0;l<dofe[i][j];l++){
	    s+=lcm[k][l];
	    lcm[k][l]=0.0;
	  }
	  lcm[k][k]=s;
	}
      }
      ////////////////////////////////////////////////////

      codnum (rcn,i);
      codnum (ccn,j);
      mat_localize (cm,lcm,rcn,ccn);
      destrm (lcm);
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
void linbart::res_convection_vector (vector &f,long lcid,long eid,long leid)
{
  long *cn;
  vector lf;

  cn = new long [dofe[lcid][lcid]];
  allocv (dofe[lcid][lcid],lf);
  convection_vector (lf,lcid,eid,leid,lcid,lcid);
  codnum (cn,lcid);
  locglob (f.a,lf.a,cn,dofe[lcid][lcid]);
  destrv (lf);
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
void linbart::res_transmission_vector (vector &f,long lcid,long eid,long leid)
{
  long *cn,i;
  vector lf;

  cn = new long [dofe[lcid][lcid]];
  codnum (cn,lcid);

  //coupled b.c.
  for (i=0;i<ntm;i++){    
    allocv (dofe[lcid][lcid],lf);
    
    if (Tt->elements[eid].transi[i]==1)
      //if ((Tt->elements[eid].transi[i]==1) && (Tb->lc[i].elemload[leid].bc != NULL))
      //if ((Tt->elements[eid].transi[i]==1) && (Tt->elements[eid].transi[lcid]==1))
      transmission_vector (lf,i,eid,leid,lcid,i);
    
    locglob (f.a,lf.a,cn,dofe[lcid][lcid]);
    destrv (lf);
       
  }
  
  delete [] cn;
  
  /*  //old 
      cn = new long [dofe[lcid][lcid]];
      allocv (dofe[lcid][lcid],lf);
      transmission_vector (lf,lcid,eid,leid,lcid,lcid);
      codnum (cn,lcid);
      locglob (f.a,lf.a,cn,dofe[lcid][lcid]);
      destrv (lf);
      delete [] cn;
  */
}

/**
   function assembles resulting element source vector

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
  allocv (dofe[lcid][lcid],lsv);
  codnum (cn,lcid);
  
  quantity_source_vector (lsv,nodval,eid,lcid,lcid);
  locglob (sv.a,lsv.a,cn,dofe[lcid][lcid]);
  
  destrv (lsv);
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
  long i,*cn;
  vector lif,tdnv;
  matrix cm;
  
  for (i=0;i<ntm;i++){
    cn = new long [dofe[i][i]];
    allocv (dofe[i][i],lif);
    codnum (cn,i);
    internal_fluxes (i,eid,lif);
    locglob (elemif.a,lif.a,cn,dofe[i][i]);
     delete [] cn;
     destrv (lif);
  } 
  
  cn = new long [ndofe];
  Gtt->give_code_numbers (eid,cn);
  allocv (ndofe,tdnv);
  allocv (ndofe,lif);
  allocm (ndofe,ndofe,cm);
  
  nodalderivatives (tdnv.a,cn,ndofe);
  res_capacity_matrix (eid,cm);
  mxv (cm,tdnv,lif);
  
  subv (elemif,lif,elemif);
  
  destrm (cm);
  destrv (lif);  destrv (tdnv);
}


/**
   function computes element quantity integral
   
   @param eid - element id
   @param nodval - %vector of quantity nodal values

   @retval f - element quantity integral

   TKr, 30.1.2004
*/
double linbart::total_integral(long eid,vector &nodval)
{
  long i;
  double area,xi,ww,jac,value,f;
  ivector nodes(nne);
  vector x(nne),w(2),gp(2);

  Tt->give_elemnodes (eid,nodes);
  Tc->give_area (eid,area);
  Tt->give_node_coord1d (x,eid);

  gauss_points (gp.a,w.a,2);
  
  f = 0.0;

  for (i=0;i<2;i++){
    xi=gp[i];  ww=w[i];
    
    jac_1d (jac,x,xi);

    value = approx (xi,nodval);
    
    jac*=area*ww*value;
    
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
   @param ri,ci - row and column indeces of the computed block in the resulting matrix
   
   TKr, 28.2.2004
*/
void linbart::boundary_flux (vector &tmv,long lcid,long eid,long leid,long ri,long ci)
{
  long i,*bc,ii;
  double xi,area;
  double new_nodval,tr;
  long nn;
  
  ivector nodes(nne);
  vector x(nne),trc(2),trcn(nne),nodval(2),av(dofe[ri][ci]),v(dofe[ri][ci]);
  vector trr(2);
  matrix n(1,dofe[ri][ci]),km(dofe[ri][ci],dofe[ri][ci]);
  
  bc = new long [2];

  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord1d (x,eid);
  Tc->give_area (eid,area);
  
  for (i=0;i<2;i++){
    bc[i]=Tb->lc[lcid].elemload[leid].bc[i];
  }
  for (i=0;i<2;i++){
    trc[i]=Tb->lc[lcid].elemload[leid].trc[i];
  }
  /*  for (i=0;i<2;i++){
      nodval[i]=Tb->lc[lcid].elemload[leid].nodval[i];
      }
  */
  Tb->lc[lcid].elemload[leid].getval(lcid,nodval);

  for (i=0;i<2;i++){
    trr[i]=Tb->lc[lcid].elemload[leid].trr[i];
  }  

  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];

  if (bc[0]>=30){
    fillm (0.0,km);
    xi=-1.0;
    
    bf_matrix (n,xi);

    nn = nodes[0];
    tr = 0.0;
    if(bc[0] == 90)
      tr = trr[0]/trc[0];
    
    Tm->transmission_flux(new_nodval,nodval[0],tr,ri,ci,nn,bc[0],ii);
   
    av[0]=new_nodval;  av[1]=0.0;
    
    nnj (km.a,n.a,trc[0]*area,n.m,n.n);
    
    mxv (km,av,v);  addv (tmv,v,tmv);
  }
  if (bc[1]>=30){
    fillm (0.0,km);
    xi=1.0;
    
    bf_matrix (n,xi);

    nn = nodes[1];
    tr = 0.0;
    if(bc[0] == 90)
      tr = trr[1]/trc[1];
    
    Tm->transmission_flux(new_nodval,nodval[1],tr,ri,ci,nn,bc[1],ii);
   
    av[0]=0.0;  av[1]=new_nodval;
          
    nnj (km.a,n.a,trc[1]*area,n.m,n.n);

    mxv (km,av,v);  addv (tmv,v,tmv);
  }
  
  delete [] bc;
}

/**
   function assembles resulting element boundary flux vector

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
    allocv (dofe[lcid][lcid],lf);
    
    if (Tt->elements[eid].transi[i]==1)
    //if ((Tt->elements[eid].transi[i]==1) && (Tb->lc[i].elemload[leid].bc != NULL))
      boundary_flux (lf,i,eid,leid,lcid,i);
    
    locglob (f.a,lf.a,cn,dofe[lcid][lcid]);
    destrv (lf);
    
  }
  
  delete [] cn;
   
/*  //old 
    cn = new long [dofe[lcid][lcid]];
    allocv (dofe[lcid][lcid],lf);
    boundary_flux (lf,lcid,eid,leid,lcid,lcid);
    codnum (cn,lcid);
    locglob (f.a,lf.a,cn,dofe[lcid][lcid]);
    destrv (lf);
    delete [] cn;
*/
}

/**
   function computes other values in nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index

   10.5.2002
*/
void linbart::nod_others (long lcid,long eid,long ri,long ci)
{
  long i,j,k,ipp,ncompother,ii;
  vector other,h,r(ndofe);
  ivector nod(nne),cn(ndofe);
  
  Tt->give_elemnodes (eid,nod);  
  Tt->give_code_numbers (eid,cn.a);
  nodalvalues (lcid,r.a,cn.a,ndofe);
  ncompother = Tm->givencompother();
  allocv (ncompother,other);
  allocv (ntm,h);

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
    
    Tt->nodes[nod[i]].storeother (ncompother,other);
  }

  destrv (other);  destrv (cn);  destrv (nod);  destrv (r);
  destrv (h);
}
