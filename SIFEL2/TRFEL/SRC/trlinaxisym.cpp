/*
  File:     trlinaxisym.cpp
  Author:   Tomas Krejci, 31.3.2002
  Purpose:  twodimensional triangular element with linear approximation functions for axisymmetric problem
*/

#include "globalt.h"
#include "trlinaxisym.h"
#include "genfile.h"
#include "globmatt.h"

trlinaxisym::trlinaxisym (void)
{
  long i;

  //  number of nodes on element
  nne=3;
  //  number of edges
  ned=3;
  //  number of nodes on one edge
  nned=2;
  //  geometric problem dimension (2D)
  ncomp=2;

  intordb=2;
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
  case nomedium:{  break; }
  case onemedium:{
    ordering[0][0]=1;  ordering[0][1]=2;  ordering[0][2]=3;
    dofe[0][0]=3;  intordkm[0][0]=1;  intordcm[0][0]=3;  nip[0][0]=4;
    ndofe=3;  napfun=1;
    break;
  }
  case twomediacoup:{
    ordering[0][0]=1;  ordering[0][1]=3;  ordering[0][2]=5;
    ordering[1][0]=2;  ordering[1][1]=4;  ordering[1][2]=6;
    
    intordkm[0][0]=1;  intordkm[0][1]=1;  intordkm[1][0]=1;  intordkm[1][1]=1;
    intordcm[0][0]=3;  intordcm[0][1]=3;  intordcm[1][0]=3;  intordcm[1][1]=3;

    if (Tp->savemode==0){
      nip[0][0]=4;       nip[0][1]=4;       nip[1][0]=4;       nip[1][1]=4;
    }
    if (Tp->savemode==1){
      nip[0][0]=4;       nip[0][1]=0;       nip[1][0]=0;       nip[1][1]=0;
    }
    
    dofe[0][0]=3;  dofe[0][1]=3;  dofe[1][0]=3;  dofe[1][1]=3;
    ndofe=6;  napfun=2;
    break;
  }
  case threemediacoup:{
    ordering[0][0]=1;   ordering[0][1] =4;  ordering[0][2] =7;
    ordering[1][0]=2;   ordering[1][1] =5;  ordering[1][2] =8;
    ordering[2][0]=3;   ordering[2][1] =6;  ordering[2][2] =9;

    intordkm[0][0]=1;  intordkm[0][1]=1;  intordkm[0][2]=1;
    intordkm[1][0]=1;  intordkm[1][1]=1;  intordkm[1][2]=1;
    intordkm[2][0]=1;  intordkm[2][1]=1;  intordkm[2][2]=1;

    intordcm[0][0]=3;  intordcm[0][1]=3;  intordcm[0][2]=3;
    intordcm[1][0]=3;  intordcm[1][1]=3;  intordcm[1][2]=3;
    intordcm[2][0]=3;  intordcm[2][1]=3;  intordcm[2][2]=3;
    
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
    
    dofe[0][0]=3;  dofe[0][1]=3;  dofe[0][2]=3;
    dofe[1][0]=3;  dofe[1][1]=3;  dofe[1][2]=3;
    dofe[2][0]=3;  dofe[2][1]=3;  dofe[2][2]=3;

    ndofe=9;  napfun=3;
    break;
  }
  default:{
    print_err("unknown type of transported matter is required",__FILE__,__LINE__,__func__);
  }
  }

}

trlinaxisym::~trlinaxisym (void)
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

void trlinaxisym::codnum (long *cn,long ri)
{
  long i;
  for (i=0;i<nne;i++){
    cn[i]=ordering[ri][i];
  }
}

/**
   function returns coordinates of integration points

   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param coord - array of coordinates

   19.1.2002
*/
void trlinaxisym::ipcoordblock (long eid,long ri,long ci,double **coord)
{
  vector x(nne),y(nne),w(intordkm[ri][ci]),gp1(intordkm[ri][ci]),gp2(intordkm[ri][ci]);
  
  gauss_points_tr (gp1.a,gp2.a,w.a,intordkm[ri][ci]);
  Tt->give_node_coord2d (x,y,eid);
  
  for (long i=0;i<intordkm[ri][ci];i++){
    coord[i][0]=approx_nat (gp1[i],gp2[i],x);
    coord[i][1]=approx_nat (gp1[i],gp2[i],y);
    coord[i][2]=0.0;
  }
}

/**
   function computes area of the element
   
   @param eid - element id
   
   05/11/2020, according to JK trlineart.cpp
*/
double trlinaxisym::element_area (long eid)
{
  double area;
  vector x(nne),y(nne);
  
  Tt->give_node_coord2d (x,y,eid);
  
  //  det is equal to double area of the element
  area = ((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]))/2.0;
  
  return area;
}

/**
   function approximates function defined by nodal values

   @param areacoord - vector containing area coordinates
   @param nodval - nodal values

   JK, 25.9.2001
*/
double trlinaxisym::approx (vector &areacoord,vector &nodval)
{
  double f;
  scprd (areacoord,nodval,f);
  return f;
}

/**
   function approximates function defined by nodal values

   @param xi,eta - natural coordinates
   @param nodval - nodal values

*/
double trlinaxisym::approx_nat (double xi,double eta,vector &nodval)
{
  double f;
  vector areacoord(3);
  areacoord[0]=xi;
  areacoord[1]=eta;
  areacoord[2]=1.0-areacoord[0]-areacoord[1];
  scprd (areacoord,nodval,f);
  return f;
}

/**
   function computes values in integration points from nodal values

   eid - element id

   JK, 25.9.2001
*/
void trlinaxisym::intpointval (long eid)
{
  long i,k,ii,jj,ipp;
  double val;
  vector areacoord(ASTCKVEC(3)), r(ASTCKVEC(ndofe)), t(ASTCKVEC(nne)),gp1,gp2,w;

  elemvalues(eid, r);

  for (k=0;k<Tp->ntm;k++){

    //  nodal values of a single variable
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }

    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){

	reallocv (RSTCKVEC(intordkm[ii][jj],gp1));
	reallocv (RSTCKVEC(intordkm[ii][jj],gp2));
	reallocv (RSTCKVEC(intordkm[ii][jj],w));
	gauss_points_tr (gp1.a,gp2.a,w.a,intordkm[ii][jj]);

	ipp=Tt->elements[eid].ipp[ii][jj];
	for (i=0;i<intordkm[ii][jj];i++){
	  areacoord[0]=gp1[i];  areacoord[1]=gp2[i];
	  areacoord[2]=1.0-areacoord[0]-areacoord[1];

	  val = approx (areacoord,t);
	  Tm->ip[ipp].av[k]=val;
	  ipp++;
	}
	
	reallocv (RSTCKVEC(intordcm[ii][jj],gp1));
	reallocv (RSTCKVEC(intordcm[ii][jj],gp2));
	reallocv (RSTCKVEC(intordcm[ii][jj],w));
	gauss_points_tr (gp1.a,gp2.a,w.a,intordcm[ii][jj]);
	
	for (i=0;i<intordcm[ii][jj];i++){
	  areacoord[0]=gp1[i];  areacoord[1]=gp2[i];
	  areacoord[2]=1.0-areacoord[0]-areacoord[1];
	  
	  val = approx (areacoord,t);
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

   eid - element id

   05/11/2020 TKr
*/
void trlinaxisym::initintpointval (long eid)
{
  long i,k,ii,jj,ipp,ndofn,cndofn;
  double val;
  vector areacoord(ASTCKVEC(3)), r(ASTCKVEC(ndofe)), t(ASTCKVEC(nne)),gp1,gp2,w;
  ivector enod(ASTCKIVEC(nne));
  vector rr;

  //  initial nodal values on element
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

  for (k=0;k<Tp->ntm;k++){

    //  nodal values of a single variable
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }

    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){

	reallocv (RSTCKVEC(intordkm[ii][jj],gp1));
	reallocv (RSTCKVEC(intordkm[ii][jj],gp2));
	reallocv (RSTCKVEC(intordkm[ii][jj],w));
	gauss_points_tr (gp1.a,gp2.a,w.a,intordkm[ii][jj]);

	ipp=Tt->elements[eid].ipp[ii][jj];
	for (i=0;i<intordkm[ii][jj];i++){
	  areacoord[0]=gp1[i];  areacoord[1]=gp2[i];
	  areacoord[2]=1.0-areacoord[0]-areacoord[1];

	  val = approx (areacoord,t);
	  Tm->ip[ipp].av[k]=val;
	  ipp++;
	}
	
	reallocv (RSTCKVEC(intordcm[ii][jj],gp1));
	reallocv (RSTCKVEC(intordcm[ii][jj],gp2));
	reallocv (RSTCKVEC(intordcm[ii][jj],w));
	gauss_points_tr (gp1.a,gp2.a,w.a,intordcm[ii][jj]);
	
	for (i=0;i<intordcm[ii][jj];i++){
	  areacoord[0]=gp1[i];  areacoord[1]=gp2[i];
	  areacoord[2]=1.0-areacoord[0]-areacoord[1];
	  
	  val = approx (areacoord,t);
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
   
   eid - element id

   JK, 25.9.2001
*/
void trlinaxisym::intpointgrad (long eid)
{
  long i,k,ii,jj,ipp;
  double det;
  vector areacoord(ASTCKVEC(3)), r(ASTCKVEC(ndofe)), t(ASTCKVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector gp1, gp2, w, b(ASTCKVEC(3)), c(ASTCKVEC(3)), grad(ASTCKVEC(2));
  matrix gm(ASTCKMAT(ncomp,nne));
  
  elemvalues(eid, r);  
  Tt->give_node_coord2d (x,y,eid);
  
  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  plsb (b.a,y.a,det);
  plsc (c.a,x.a,det);
  
  for (k=0;k<Tp->ntm;k++){
    
    //  nodal values of a single variable
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
	
	reallocv (RSTCKVEC(intordkm[ii][jj],gp1));
	reallocv (RSTCKVEC(intordkm[ii][jj],gp2));
	reallocv (RSTCKVEC(intordkm[ii][jj],w));
	gauss_points_tr (gp1.a,gp2.a,w.a,intordkm[ii][jj]);

	ipp=Tt->elements[eid].ipp[ii][jj];
	for (i=0;i<intordkm[ii][jj];i++){
	  areacoord[0]=gp1[i];  areacoord[1]=gp2[i];
	  areacoord[2]=1.0-areacoord[0]-areacoord[1];

	  // gradient matrix
	  grad_matrix (gm,b,c);

	  mxv (gm,t,grad);

	  Tm->storegrad (k,ipp,grad);

	  ipp++;
	}

	reallocv (RSTCKVEC(intordcm[ii][jj],gp1));
	reallocv (RSTCKVEC(intordcm[ii][jj],gp2));
	reallocv (RSTCKVEC(intordcm[ii][jj],w));
	gauss_points_tr (gp1.a,gp2.a,w.a,intordcm[ii][jj]);

	for (i=0;i<intordcm[ii][jj];i++){
	  areacoord[0]=gp1[i];  areacoord[1]=gp2[i];
	  areacoord[2]=1.0-areacoord[0]-areacoord[1];

	  // gradient matrix
	  grad_matrix (gm,b,c);

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
void trlinaxisym::intpointother (long eid)
{
  long i, k, ii, jj, ipp, ncompo, nodid;
  double val;
  ivector nodes(nne);
  vector r, areacoord(3), t(nne), gp1, gp2, w;

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

	reallocv (intordkm[ii][jj],gp1);
	reallocv (intordkm[ii][jj],gp2);
	reallocv (intordkm[ii][jj],w);
	gauss_points_tr (gp1.a,gp2.a,w.a,intordkm[ii][jj]);

	ipp=Tt->elements[eid].ipp[ii][jj];
	for (i=0;i<intordkm[ii][jj];i++){
	  areacoord[0]=gp1[i];  areacoord[1]=gp2[i];
	  areacoord[2]=1.0-areacoord[0]-areacoord[1];

	  val = approx (areacoord,t);
	  Tm->ip[ipp].other[k]=val;
	  ipp++;

	}
	
	reallocv (intordcm[ii][jj],gp1);
	reallocv (intordcm[ii][jj],gp2);
	reallocv (intordcm[ii][jj],w);
	gauss_points_tr (gp1.a,gp2.a,w.a,intordcm[ii][jj]);
	
	for (i=0;i<intordcm[ii][jj];i++){
	  areacoord[0]=gp1[i];  areacoord[1]=gp2[i];
	  areacoord[2]=1.0-areacoord[0]-areacoord[1];
	  
	  val = approx (areacoord,t);
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
   function returns %matrix of base function

   @param n - array of approximation functions
   @param areacoord - area coordinates

   JK, 5.11.2001
*/
void trlinaxisym::bf_matrix (matrix &n,vector &areacoord)
{
  nullm (n);

  n[0][0]=areacoord[0];
  n[0][1]=areacoord[1];
  n[0][2]=areacoord[2];
}

/**
   function returns one approximation function evaluated in required point

   @param f - value of approximation function
   @param areacoord - area coordinates
   @param i - number of approximation function

   JK, 1.2.2003
*/
void trlinaxisym::give_approx_fun (double &f,vector &areacoord,long i)
{
  f=areacoord[i];
}

/**
   function computes capacity %matrix for one medium

   @param eid - element id
   @param ri,ci - row and column indices
   @param cm - capacity %matrix

   JK, 5.11.2001
*/
void trlinaxisym::capacity_matrix (long eid,long ri,long ci,matrix &cm)
{
  long i,ii;
  double jac,det,xinp,rho,c;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),w(intordcm[ri][ci]),gp1(intordcm[ri][ci]),gp2(intordcm[ri][ci]),t(nne),dens(nne);
  matrix n(1,dofe[ri][ci]);

  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord2d (x,y,eid);
  
  Tc->give_density (eid,nodes,dens);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  gauss_points_tr (gp1.a,gp2.a,w.a,intordcm[ri][ci]);

  nullm (cm);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci]+intordkm[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0]+intordkm[0][0];

  for (i=0;i<intordcm[ri][ci];i++){
    areacoord[0]=gp1[i];
    areacoord[1]=gp2[i];
    areacoord[2]=1.0-gp1[i]-gp2[i];

    bf_matrix (n,areacoord);

    xinp = approx (areacoord,x);
    rho = approx (areacoord,dens);
    c=Tm->capcoeff (ii,ri,ci);

    jac=w[i]*xinp*rho*det*c;

    nnj (cm.a,n.a,jac,n.m,n.n);

    ii++;
  }

}

/**
   function assembles gradient of %matrix of base functions

   @param gm - gradient %matrix
   @param b,c - array containing node coordinates

   JK, 25.9.2001
*/
void trlinaxisym::grad_matrix (matrix &gm,vector &b,vector &c)
{
  long i;

  nullm (gm);

  for (i=0;i<nne;i++){
    gm[0][i]=c[i];
    gm[1][i]=b[i];
  }
}


/**
   function computes conductivity %matrix of one medium

   @param lcid - load case id
   @param eid - element id
   @param km - conductivity %matrix

   JK, 5.11.2001
*/
void trlinaxisym::conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long ii;
  double jac,det,xinp;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),b(3),c(3),t(nne);
  matrix gm(ncomp,dofe[ri][ci]),d(ncomp,ncomp);

  matrix n(1,dofe[ri][ci]);

  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord2d (x,y,eid);
  

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  areacoord[0]=1.0/3.0;
  areacoord[1]=1.0/3.0;
  areacoord[2]=1.0/3.0;

  plsb (b.a,y.a,det);
  plsc (c.a,x.a,det);

  nullm (km);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];

  // gradient matrix
  grad_matrix (gm,b,c);

  //  matrix of conductivity of material
  Tm->matcond (d,ii,ri,ci);

  xinp = approx (areacoord,x);

  //  det is equal to double area of the element
  jac=xinp*det/2.0;

  //  contribution to the stiffness matrix of the element
  bdbj (km.a,gm.a,d.a,jac,gm.m,gm.n);
  
  //convective terms
  reallocm(1,ncomp,d);
  
  Tm->matcond2(d,ii,ri,ci);
  bf_matrix (n, areacoord);
  bdbjac(km, n, d, gm, jac);	
  
  if ((Tt->elements[eid].transi[lcid]==3)||(Tt->elements[eid].transi[lcid]==4)){
    transmission_matrix (lcid,eid,ri,ci,km);
  }
}

/**
   function computes source %vector of one matter on one element

   \int_{Omega} N^T N d Omega . s

   @param sv - source %vector of one matter
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)

   JK, 4.10.2001
*/
void trlinaxisym::quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci)
{
  long i;
  double jac,det,xinp;
  vector x(nne),y(nne),areacoord(3),w(intordcm[ri][ci]),gp1(intordcm[ri][ci]),gp2(intordcm[ri][ci]),t(nne),v(dofe[ri][ci]);
  matrix n(1,dofe[ri][ci]),nm(dofe[ri][ci],dofe[ri][ci]);
  
  Tt->give_node_coord2d (x,y,eid);
  gauss_points_tr (gp1.a,gp2.a,w.a,intordcm[ri][ci]);
  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
    
  nullm (nm);

  for (i=0;i<intordcm[ri][ci];i++){    
    areacoord[0]=gp1[i];
    areacoord[1]=gp2[i];
    areacoord[2]=1.0-gp1[i]-gp2[i];
    
    bf_matrix (n,areacoord);
    
    xinp = approx (areacoord,x);
    
    jac=w[i]*xinp*det;
    
    nnj (nm.a,n.a,jac,n.m,n.n); 
  }   
  mxv (nm,nodval,v);  addv (sv,v,sv);

}



/**
   function computes internal fluxes of 1D problems for one transported matter
   finite element with linear approximation functions
   
   @param lcid - number of load case
   @param eid - number of element
   @param ifl - %vector of internal fluxes
   
   JK, 31.3.2002
*/
void trlinaxisym::internal_fluxes (long lcid,long eid,vector &ifl)
{
  long ipp;
  double det,xinp;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),b(3),c(3),t(nne),fl(ncomp),contr(dofe[lcid][lcid]);
  matrix gm(ncomp,dofe[lcid][lcid]);
  
  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord2d (x,y,eid);
  
  
  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  
  areacoord[0]=1.0/3.0;
  areacoord[1]=1.0/3.0;
  areacoord[2]=1.0/3.0;

  plsb (b.a,y.a,det);
  plsc (c.a,x.a,det);
  
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[lcid][lcid];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  xinp = approx (areacoord,x);
  
  Tm->computenlfluxes (lcid,ipp);
  
  Tm->givefluxes (lcid,ipp,fl);

  // gradient matrix
  grad_matrix (gm,b,c);
  
  mtxv (gm,fl,contr);
  
  cmulv (det/2.0*xinp,contr);
  
  addv (contr,ifl,ifl);
}

/**
   function assembles resulting element conductivity %matrix

   @param eid - element id
   @param lcid - load case id
   @param km - resulting conductivity %matrix of one element

   JK, 6.1.2002
*/
void trlinaxisym::res_conductivity_matrix (long eid,long /*lcid*/,matrix &km)
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
   function assembles resulting element capacity %matrix
   
   @param eid - element id
   @param cm - resulting capacity %matrix of one element

   JK, 6.1.2002
*/
void trlinaxisym::res_capacity_matrix (long eid,matrix &cm)
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
   function assembles resulting element convection %vector

   @param f - resulting convection %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   
   JK, 6.1.2002
*/
void trlinaxisym::res_convection_vector (vector &f,long lcid,long eid,long leid)
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
   function assembles resulting element transmission %vector

   @param f - resulting transmission %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   
   JK, 6.1.2002
*/
void trlinaxisym::res_transmission_vector (vector &f,long lcid,long eid,long leid)
{
  long *cn,i;
  vector lf;

  cn = new long [dofe[lcid][lcid]];
  codnum (cn,lcid);

  //coupled b.c.
  for (i=0;i<ntm;i++){    
    reallocv (dofe[lcid][lcid],lf);
    
    if ((Tt->elements[eid].transi[i]==3)||(Tt->elements[eid].transi[i]==4)){
      transmission_vector (lf,i,eid,leid,lcid,i);
    }
    
    locglob (f.a,lf.a,cn,dofe[lcid][lcid]);
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
void trlinaxisym::res_quantity_source_vector (vector &sv,vector &nodval,long lcid,long eid)
{
  long *cn;
  vector lsv;

  cn = new long [dofe[lcid][lcid]];
  reallocv (dofe[lcid][lcid],lsv);
  quantity_source_vector (lsv,nodval,eid,lcid,lcid);
  codnum (cn,lcid);
  locglob (sv.a,lsv.a,cn,dofe[lcid][lcid]);
  delete [] cn;
}

/**
   function assembles resulting element internal fluxes %vector

   @param eid - element id
   @param elemif - resulting internal fluxes %vector of one element
   
   JK, 6.1.2002
*/
void trlinaxisym::res_internal_fluxes (long eid,vector &elemif)
{
  long i,*cn;
  vector lif,tdnv;
  matrix cm;
  
  for (i=0;i<ntm;i++){
    cn = new long [dofe[i][i]];
    reallocv (dofe[i][i],lif);
    internal_fluxes (i,eid,lif);
    codnum (cn,i);
    locglob (elemif.a,lif.a,cn,dofe[i][i]);
    delete [] cn;
  }
  
  cn = new long [ndofe];
  Gtt->give_code_numbers (eid,cn);
  reallocv (ndofe,tdnv);
  reallocv (ndofe,lif);
  reallocm (ndofe,ndofe,cm);
  nodalderivatives(eid, tdnv);
  res_capacity_matrix (eid,cm);
  mxv (cm,tdnv,lif);
  subv (elemif,lif,elemif);
}

/**
   function computes element quantity integral
   
   @param eid - element id
   @param nodval - %vector of quantity nodal values

   @retval f - element quantity integral

   TKr, 30.1.2004
*/
double trlinaxisym::total_integral (long eid,vector &nodval)
{
  long i;
  double jac,det,xinp,value,f;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),w(3),gp1(3),gp2(3),t(nne);

  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord2d (x,y,eid);
  

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  gauss_points_tr (gp1.a,gp2.a,w.a,3);

  f = 0.0;

  for (i=0;i<3;i++){
    areacoord[0]=gp1[i];
    areacoord[1]=gp2[i];
    areacoord[2]=1.0-gp1[i]-gp2[i];

    xinp = approx (areacoord,x);
    value = approx (areacoord,nodval);

    jac=w[i]*xinp*value*det;

    f = f + jac;
  }
  return(f);
}


/**
   function assembles resulting element boundary flux %vector

   @param f - resulting boundary flux %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   
   TKr, 28.2.2004
*/
void trlinaxisym::res_boundary_flux (vector &f,long lcid,long eid,long leid)
{
  long *cn,i;
  vector lf;

  cn = new long [dofe[lcid][lcid]];
  codnum (cn,lcid);
  
  //coupled b.c.
  for (i=0;i<ntm;i++){    
    reallocv (dofe[lcid][lcid],lf);    
    
    if ((Tt->elements[eid].transi[i]==3)||(Tt->elements[eid].transi[i]==4)){
      boundary_flux (lf,i,eid,leid,lcid,i);
    }
    
    locglob (f.a,lf.a,cn,dofe[lcid][lcid]);
  }
  
  delete [] cn;
}



/**
   function computes contributions to the right-hand side - volume integral
      
   \int_{Omega} B^T D dOmega
   
   @param vrhs - volume right-hand side %vector of one matter
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)
   
   TKr, 06/11/2020
*/
void trlinaxisym::volume_rhs_vector (long /*lcid*/,long eid,long ri,long ci,vector &vrhs)
{
  long i,ii;
  double jac,det,xinp;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),b(3),c(3),t(nne);
  matrix gm(ncomp,dofe[ri][ci]),d;
  matrix km(ASTCKMAT(dofe[ri][ci], 2));
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  node coordinates
  Tt->give_node_coord2d (x,y,eid);
  //  thickness
  Tc->give_thickness (eid,nodes,t);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  if (det<0.0){
    fprintf (stdout,"\n negative determinant of element %ld",eid);
    fprintf (stdout,"\n det = %e",det);
    
    det=fabs(det);
    abort();
  }

  areacoord[0]=1.0/3.0;
  areacoord[1]=1.0/3.0;
  areacoord[2]=1.0/3.0;

  plsb (b.a,y.a,det);
  plsc (c.a,x.a,det);

  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];
  
  
  // gradient matrix
  grad_matrix (gm,b,c);

  //  matrix of conductivity of material
  reallocm (ncomp,ncomp,d);

  Tm->volume_rhs (d,ii,ri,ci,ncomp);

  xinp = approx (areacoord,x);

  //  det is equal to double area of the element
  jac=xinp*det/2.0;
  
  //  contribution to the volume_rhs integral of the element
  nnjac (km,gm,d,jac);
      
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
   
   TKr, 06/11/2020
*/
void trlinaxisym::volume_rhs_vector2 (long /*lcid*/,long eid,long ri,long ci,vector &vrhs)
{
  long i,k,ii;
  double jac,det,xinp,c;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),w(intordcm[ri][ci]),gp1(intordcm[ri][ci]),gp2(intordcm[ri][ci]),t(nne),dens(nne);
  matrix n(1,dofe[ri][ci]);
  vector nn(ASTCKVEC(dofe[ri][ci]));

  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord2d (x,y,eid);
  Tc->give_thickness (eid,nodes,t);
  Tc->give_density (eid,nodes,dens);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  gauss_points_tr (gp1.a,gp2.a,w.a,intordcm[ri][ci]);

  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci]+intordkm[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0]+intordkm[0][0];
  
  
  for (i=0;i<intordcm[ri][ci];i++){
    areacoord[0]=gp1[i];
    areacoord[1]=gp2[i];
    areacoord[2]=1.0-gp1[i]-gp2[i];
    
    bf_matrix (n,areacoord);
    
    xinp = approx (areacoord,x);
    
    c=Tm->volume_rhs2 (ii,ri,ci);
    
    jac=w[i]*xinp*det*c;
    
    for (k=0;k<n.n;k++){
      nn[k] = n[0][k];
    }
    
    addmultv(vrhs,nn,jac);
    
    ii++;
  }
}



/**
   function assembles resulting element volume right-hand side

   @param eid - element id
   @param lcid - load case id
   @param f - resulting volume right-hand side %vector of one element

   TKr, 05/11/2020
*/
void trlinaxisym::res_volume_rhs_vector (vector &f,long eid,long /*lcid*/)
{
  long i;
  ivector cn;
  vector lf;

  for (i=0;i<ntm;i++){
    reallocv(RSTCKIVEC(dofe[i][i], cn));
    reallocv(RSTCKVEC(dofe[i][i], lf));
    codnum (cn.a,i);
    volume_rhs_vector (i,eid,i,i,lf);
    locglob (f.a,lf.a,cn.a,dofe[i][i]);
  }
}



/**
   function assembles resulting element volume right-hand side of the second type

   @param eid - element id
   @param lcid - load case id
   @param f - resulting volume right-hand side %vector of one element

   TKr, 05/11/2020
*/
void trlinaxisym::res_volume_rhs_vector2 (vector &f,long eid,long /*lcid*/)
{
  long i;
  ivector cn;
  vector lf;

  for (i=0;i<ntm;i++){
    reallocv(RSTCKIVEC(dofe[i][i], cn));
    reallocv(RSTCKVEC(dofe[i][i], lf));
    codnum (cn.a,i);
    volume_rhs_vector2 (i,eid,i,i,lf);
    locglob (f.a,lf.a,cn.a,dofe[i][i]);
  }
}




/**
   function computes correct fluxes at integration points on element

   @param eid - element id
   
   TKr, 01/02/2010
*/
void trlinaxisym::intpointflux (long eid)
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
   function computes other values in nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index

   10.5.2002
*/
void trlinaxisym::nod_others (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp,ncompother,ii;
  vector other,h,r(ASTCKVEC(ndofe));
  ivector nod(ASTCKIVEC(nne));
  
  Tt->give_elemnodes (eid,nod);  
  elemvalues(eid, r);
  ncompother = Tm->givencompother();
  reallocv(RSTCKVEC(ncompother,other));
  reallocv(RSTCKVEC(ntm,h));

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
   function computes other values in nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index

   05/11/2020 TKr according to trlineart.cpp
*/
void trlinaxisym::nod_others_comp (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp,ncompother,ii;
  vector h,r(ASTCKVEC(ndofe));
  ivector nod(ASTCKIVEC(nne));
  double area,other;
  
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
    
    //  storage of eqother components to the node
    j=nod[i];
    for (k=0; k<ncompother; k++){
      other = Tm->givecompother(k,ipp,h.a);
      if (Tp->otheraver==0 || Tp->otheraver==1)
	Tt->nodes[j].storeother (k,other);
      if (Tp->otheraver==2){
	area = element_area (eid);
	Tt->nodes[j].storeother (k,area,other);
      }
    }
  }
}










/**
   function computes contribution to the convection %vector

   \int_{Gamma_2} N^T N dGamma * nodal_flux_values
   
   @param v - array of nodal fluxes
   @param lcid - load case id
   @param eid - element id
   @param leid - id of loaded element
   @param ri,ci - row and column indices
   
   JK, 19.8.2004
*/
void trlinaxisym::convection_vector (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]),t(nne);
  vector list(nned*ned),trc(nned*ned),nodval(nne),av(nne),coef(nne);
  matrix km(nne,nne);
  
  nullv (v);
  
  //  element nodes
  Tt->give_elemnodes (eid,nodes);
  //  node coordinates on element
  Tt->give_node_coord2d (x,y,eid);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  
  //  indicators of boundary conditions
  bc = new bocontypet[ned];
  Tb->lc[lcid].elemload[leid].give_bc (bc);
  //  nodal values of exterior variable
  Tb->lc[lcid].elemload[leid].give_nodval (lcid,list);
  
  //  auxiliary coefficients, necessary for function edge_integral
  for (i=0;i<nned*ned;i++){
    trc[i]=1.0;
  }

  //  integration point number
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  //  loop over edges
  for (i=0;i<ned;i++){
    
    if (bc[i]==2 || bc[i]==3){
      //  nodal values on actual edge
      edgenodeval (i,nodval,list);
      
      edgenodeval (i,coef,trc);
  
      nullm (km);

      //  matrix obtained from integration over edge
      edge_integral (i,x,y,intordkm[ri][ci],gp,w,coef,km);
      
      mxv (km,nodval,av);
      
      addv (v,av,v);
    }
  }
  
  delete [] bc;
}

/**
   function computes transmission complement to the conductivity %matrix for one matter
   
   \int_{Gamma_3} N^T c_{tr} N dGamma
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param km - conductivity %matrix
   
   JK, 19.8.2004
*/
void trlinaxisym::transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,ipp,leid;
  bocontypet *bc;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]);
  vector trc(nned*ned),coeff(nne),t(nne);
  
  //  element nodes
  Tt->give_elemnodes (eid,nodes);
  //  node coordinates on element
  Tt->give_node_coord2d (x,y,eid);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);

  //  correspondence between element ordering and boundary element ordering
  for (i=0;i<Tb->lc[lcid].neb;i++){
    if (Tb->lc[lcid].elemload[i].eid==eid){
      leid=i;
      break;
    }
  }

  //  indicators of boundary conditions
  bc = new bocontypet[ned];
  Tb->lc[lcid].elemload[leid].give_bc (bc);
  //  transmission coefficients
  Tb->lc[lcid].elemload[leid].give_trc (Tp->time,Tb->lc[lcid].nodval,lcid,ci,trc);
  
  //  integration point number
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];

  //  loop over edges
  for (i=0;i<ned;i++){
    
    if (bc[i]>10){
      //  transformation of coefficients
      transf_coeff (i,coeff,trc,eid,ri,ci,ipp,bc);
      
      //  matrix obtained from integration over edges
      edge_integral (i,x,y,intordkm[ri][ci],gp,w,coeff,km);
    }
  }
  
  delete [] bc;
}

/**
   function computes contributions to the transmission %vector
   
   \int_{Gamma_3} N^T c_{tr} N dGamma * nodal_external_value
   
   @param tmv - transmission %vector of one matter
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)
   
   JK, 5.10.2001
   TKr, 30.1.2002 - new added
*/
void trlinaxisym::transmission_vector (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]),t(nne);
  vector list(nned*ned),trc(nned*ned),trr(nned*ned),nodval(nne),av(nne),coef(nne);
  matrix km(nne,nne);
  
  nullv (v);
  
  //  element nodes
  Tt->give_elemnodes (eid,nodes);
  //  node coordinates on element
  Tt->give_node_coord2d (x,y,eid);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  
  //  indicators of boundary conditions
  bc = new bocontypet[ned];
  Tb->lc[lcid].elemload[leid].give_bc (bc);
  //  nodal values of exterior variable
  Tb->lc[lcid].elemload[leid].give_nodval (lcid,list);
  //  transmission coefficients
  Tb->lc[lcid].elemload[leid].give_trc (Tp->time,Tb->lc[lcid].nodval,lcid,ci,trc);
  //  transmission/radiation coefficients
  Tb->lc[lcid].elemload[leid].give_trr (Tp->time,Tb->lc[lcid].nodval,trr);
  
  //  integration point number
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  for (i=0;i<ned;i++){
    
    if (bc[i]>10){
      //  transformation of nodal values
      transf_val (i,nodval,list,trc,trr,eid,ri,ci,ipp,bc);
      
      edgenodeval (i,coef,trc);
  
      nullm (km);
      //
      edge_integral (i,x,y,intordkm[ri][ci],gp,w,coef,km);
      
      mxv (km,nodval,av);
      
      addv (v,av,v);
    }
  }
  
  delete [] bc;
}


/**
   function computes contributions to the boundary flux from transmission %vector
   
   \int_{Gamma_3} N^T c_{tr} N dGamma * nodal_external_value

   @param v - %vector of boundary fluxes
   @param lcid - load case id
   @param eid - element id
   @param leid - loaded element id
   @param ri,ci - row and column indices
   
   TKr, 28.2.2004
*/
void trlinaxisym::boundary_flux (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]),t(nne);
  vector list(nned*ned),trc(nned*ned),trr(nned*ned),nodval(nne),av(nne),coef(nne);
  matrix km(nne,nne);
  
  nullv (v);
  
  //  element nodes
  Tt->give_elemnodes (eid,nodes);
  //  node coordinates on element
  Tt->give_node_coord2d (x,y,eid);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  
  //  indicators of boundary conditions
  bc = new bocontypet[ned];
  Tb->lc[lcid].elemload[leid].give_bc (bc);
  //  nodal values of exterior variable
  Tb->lc[lcid].elemload[leid].give_nodval (lcid,list);
  //  transmission coefficients
  Tb->lc[lcid].elemload[leid].give_trc (Tp->time,Tb->lc[lcid].nodval,lcid,ci,trc);
  //  transmission/radiation coefficients
  Tb->lc[lcid].elemload[leid].give_trr (Tp->time,Tb->lc[lcid].nodval,trr);
  
  //  integration point number
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  for (i=0;i<ned;i++){
    
    if (bc[i]>10){
      //  transformation of nodal values
      transf_flux (i,nodval,list,trc,trr,eid,ri,ci,ipp,bc);
      
      edgenodeval (i,coef,trc);
  
      nullm (km);
      //
      edge_integral (i,x,y,intordkm[ri][ci],gp,w,coef,km);
      
      mxv (km,nodval,av);
      
      addv (v,av,v);
    }
  }
  
  delete [] bc;
}


/**
   function integrates N^T c N over edges
   
   @param edg - number of edge
   @param x,y - arrays of nodal coordinates
   @param intord - order of numerical integration
   @param gp - array of coordinates of integration points
   @param w - array of weights
   @param coef - array of nodal values of coefficient
   @param km - resulting %matrix
   
   JK
*/
void trlinaxisym::edge_integral (long edg,vector &x,vector &y,long /*intord*/,vector &gp,vector &w,
			       vector &coef,matrix &km)
{
  long i;
  double jac,ipval,xinp;
  vector av(nne),areacoord(3);
  matrix n(1,nne);
  
  
  if (edg==0){
    for (i=0;i<intordb;i++){
      areacoord[0]=(1.0-gp[i])/2.0;  areacoord[1]=(1.0+gp[i])/2.0;
      areacoord[2]=1.0-areacoord[0]-areacoord[1];
      
      bf_matrix (n,areacoord);
      
      jac1d_2d (jac,x,y,gp[i],edg);
      xinp = approx (areacoord,x);
      ipval = approx (areacoord,coef);
      
      jac*=w[i]*xinp*ipval;
      
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

  if (edg==1){
    for (i=0;i<intordb;i++){
      areacoord[0]=0.0;  areacoord[1]=(1.0-gp[i])/2.0;
      areacoord[2]=1.0-areacoord[0]-areacoord[1];
      
      bf_matrix (n,areacoord);
      
      jac1d_2d (jac,x,y,gp[i],edg);
      xinp = approx (areacoord,x);
      ipval = approx (areacoord,coef);

      jac*=w[i]*xinp*ipval;
      
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

  if (edg==2){
    for (i=0;i<intordb;i++){
      areacoord[0]=(1.0+gp[i])/2.0;  areacoord[1]=0.0;
      areacoord[2]=1.0-areacoord[0]-areacoord[1];
      
      bf_matrix (n,areacoord);
      
      jac1d_2d (jac,x,y,gp[i],2);
      xinp = approx (areacoord,x);
      ipval = approx (areacoord,coef);

      jac*=w[i]*xinp*ipval;
      
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

}

void trlinaxisym::transf_flux (long edg,vector &coeff,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double tr,new_flux;
  ivector nodes(nne),edgenod(nned);

  //  zeroing of array of coefficients
  nullv (coeff);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required edge
  lintriangle_edgnod (edgenod.a,edg);

  //  actual position in the array list
  j=edg*nned;

  for (i=0;i<nned;i++){
    if (bc[edg]==90){
      tr = trr[j+i]/trc[j+i];
    }

    //  node number
    k=edgenod[i];
    
    Tm->transmission_flux(new_flux,list[j+i],tr,ri,ci,nodes[k],bc[edg],ipp);
    
    coeff[k]=new_flux;
  }
  
}


void trlinaxisym::transf_coeff (long edg,vector &coeff,vector &list,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double new_coeff;
  ivector nodes(nne),edgenod(nned);

  //  zeroing of array of coefficients
  nullv (coeff);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required edge
  lintriangle_edgnod (edgenod.a,edg);
  
  //  actual position in the array list
  j=edg*nned;

  for (i=0;i<nned;i++){
    //  node number
    k=edgenod[i];
    
    Tm->transmission_transcoeff(new_coeff,list[j+i],ri,ci,nodes[k],bc[edg],ipp);
    
    coeff[k]=new_coeff;
  }
  
}


/**
   @param edg - number of required edge
   @param nodval - array of transformed nodal values
   @param list - array of nodal values defined on all edges
   @param trc -
   @param trr - 
   @param ri,ci - row and column indices
   @param ipp - integration point number
   @param bc - array defining boundary conditions
   
   JK, 19.8.2004
*/
void trlinaxisym::transf_val (long edg,vector &nodval,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double tr,new_nodval;
  ivector nodes(nne),edgenod(nned);
  
  //  zeroing of array of nodal values
  nullv (nodval);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required edge
  lintriangle_edgnod (edgenod.a,edg);
  
  //  actual position in the array list
  j=edg*nned;
  
  //  loop over number of nodes on one edge
  for (i=0;i<nned;i++){
    if (bc[edg]==90){
      tr = trr[j+i]/trc[j+i];
    }
    
    //  node number
    k=edgenod[i];
    
    Tm->transmission_nodval(new_nodval,list[j+i],tr,ri,ci,nodes[k],bc[edg],ipp);
    
    //  transformed nodal values
    nodval[k]=new_nodval;
  }
  
}

/**
   function picks up nodal values on required edges
   
   @param edg - number of required edge
   @param nodval - array of nodal values
   @param list - array of nodal values defined on all edges
   
   JK, 19.8.2004
*/
void trlinaxisym::edgenodeval (long edg,vector &nodval,vector &list)
{
  long i,j;
  ivector edgenod(nned);
  
  nullv (nodval);
  lintriangle_edgnod (edgenod.a,edg);
  
  for (i=0;i<nned;i++){
    j=edgenod[i];
    nodval[j]=list[edg*nned+i];
  }
}



/**
   Function assembles global coordinates of integration points.
   
   @param eid - element id (input)
   @param ipp - integration point pointer (input)
   @param ri - row index of the block of integration points /according to transported media/ (input)
   @param ci - column index of the block of integration points /according to transported media/ (input)
   @param coord - array containing coordinates of integration points (output)
   
   @retval coord - function returns coordinates of integration points in the %vector coord.
   @retval 0 - on success
   @retval 1 - integration point with ipp was not found

   TKo, 12.2016
*/
long trlinaxisym::ipcoord (long eid, long ipp, long ri, long ci, vector &coord)
{
  long i, ii;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), w, gp1, gp2, areacoord(3);

  Tt->give_node_coord2d(x, y, eid);
  
  //  integration points for the conductivity matrix
  reallocv(RSTCKVEC(intordkm[ri][ci], gp1));
  reallocv(RSTCKVEC(intordkm[ri][ci], gp2));
  reallocv(RSTCKVEC(intordkm[ri][ci], w));
  gauss_points_tr(gp1.a, gp2.a, w.a, intordkm[ri][ci]);
	
  ii = Tt->elements[eid].ipp[ri][ci];
  for (i=0; i<intordkm[ri][ci]; i++)
  {
    if (ii == ipp)
    {
      areacoord[0]=gp1[i];  areacoord[1]=gp2[i];
      areacoord[2]=1.0-areacoord[0]-areacoord[1];
      coord[0] = approx(areacoord, x);
      coord[1] = approx(areacoord, y);
      coord[2] = 0.0;
      return 0;
    }
    ii++;
  }
	
  //  integration points for the capacity matrix
  reallocv(RSTCKVEC(intordcm[ri][ci], gp1));
  reallocv(RSTCKVEC(intordcm[ri][ci], gp2));
  reallocv(RSTCKVEC(intordcm[ri][ci], w));
  gauss_points_tr(gp1.a, gp2.a, w.a, intordcm[ri][ci]);
  
  for (i=0; i<intordcm[ri][ci]; i++)
  {
    if (ii == ipp)
    {
      areacoord[0]=gp1[i];  areacoord[1]=gp2[i];
      areacoord[2]=1.0-areacoord[0]-areacoord[1];
      coord[0] = approx(areacoord, x);
      coord[1] = approx(areacoord, y);
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
long trlinaxisym::ipncoord (long eid, long ipp, vector &ncoord)
{
  long i, ii, ri, ci;
  vector w, gp1, gp2;

  for (ri=0;ri<ntm;ri++)
  {
    for (ci=0;ci<ntm;ci++)
    {
      //  integration points for the conductivity matrix
      reallocv(RSTCKVEC(intordkm[ri][ci], gp1));
      reallocv(RSTCKVEC(intordkm[ri][ci], gp2));
      reallocv(RSTCKVEC(intordkm[ri][ci], w));
      gauss_points_tr(gp1.a, gp2.a, w.a, intordkm[ri][ci]);
	
      ii = Tt->elements[eid].ipp[ri][ci];
      for (i=0; i<intordkm[ri][ci]; i++)
      {
        if (ii == ipp)
        {
          ncoord[0] = gp1[i];
          ncoord[1] = gp2[i];
          ncoord[2] = 0.0;
          return 0;
        }
        ii++;
      }
	
      //  integration points for the capacity matrix
      reallocv(RSTCKVEC(intordcm[ri][ci], gp1));
      reallocv(RSTCKVEC(intordcm[ri][ci], gp2));
      reallocv(RSTCKVEC(intordcm[ri][ci], w));
      gauss_points_tr(gp1.a, gp2.a, w.a, intordcm[ri][ci]);
  
      for (i=0; i<intordcm[ri][ci]; i++)
      {
        if (ii == ipp)
        {
          ncoord[0] = gp1[i];
          ncoord[1] = gp2[i];
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
