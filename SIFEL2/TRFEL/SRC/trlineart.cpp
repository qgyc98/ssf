/*
  File:     trlineart.cpp
  Author:   Jaroslav Kruis, 31.3.2001
  Purpose:  twodimensional triangular element with linear approximation functions
*/

#include "globalt.h"
#include "trlineart.h"
#include "genfile.h"
#include "globmatt.h"

trlineart::trlineart (void)
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
    dofe[0][0]=3;  intordkm[0][0]=1;  intordcm[0][0]=3;
    nip[0][0]=4;
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
    print_err ("unknown type of transported matter is required",__FILE__,__LINE__,__func__);
  }
  }

}

trlineart::~trlineart (void)
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

void trlineart::codnum (long *cn,long ri)
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
void trlineart::ipcoordblock (long eid,long ri,long ci,double **coord)
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
   
   23. 3. 2017, JK
*/
double trlineart::element_area (long eid)
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

   @param areacoord - %vector containing area coordinates
   @param nodval - nodal values

   JK, 25.9.2001
*/
double trlineart::approx (vector &areacoord,vector &nodval)
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
double trlineart::approx_nat (double xi,double eta,vector &nodval)
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
void trlineart::intpointval (long eid)
{
  long i,k,ii,jj,ipp;
  double val;
  vector areacoord(ASTCKVEC(3)), r(ASTCKVEC(ndofe)), t(ASTCKVEC(nne)), gp1, gp2, w;

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
   function computes values in integration points from nodal values for PUC

   eid - element id

   TKr, 08/09/2010
*/
void trlineart::intpointval_puc (long eid)
{
  long i,k,ii,jj,ipp;
  double val;
  vector areacoord(3),r(ndofe),t(nne),gp1,gp2,w;

  elemvalues_puc(eid, r);
  
  /*
  for (long ijk=0;ijk<ndofe;ijk++){
    fprintf (stdout,"\n nodval ra %le", r.a[ijk]);
  }
  */

  for (k=0;k<Tp->ntm;k++){

    //  nodal values of a single variable
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
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
	  Tm->ip[ipp].av[k]=val;
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
void trlineart::intpointgrad (long eid)
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
void trlineart::intpointother (long eid)
{
  long i, k, ii, jj, ipp, ncompo, nodid;
  double val;
  ivector nodes(nne);
  vector areacoord(3), t(nne), r, gp1, gp2, w;

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
void trlineart::bf_matrix (matrix &n,vector &areacoord)
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
void trlineart::give_approx_fun (double &f,vector &areacoord,long i)
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
void trlineart::capacity_matrix (long eid,long ri,long ci,matrix &cm)
{
  long i,ii;
  double jac,det,thick,rho,c;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),w(intordcm[ri][ci]),gp1(intordcm[ri][ci]),gp2(intordcm[ri][ci]),t(nne),dens(nne);
  matrix n(1,dofe[ri][ci]);

  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord2d (x,y,eid);
  Tc->give_thickness (eid,nodes,t);
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

    thick = approx (areacoord,t);
    rho = approx (areacoord,dens);
    c=Tm->capcoeff (ii,ri,ci);

    jac=w[i]*thick*rho*det*c;

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
void trlineart::grad_matrix (matrix &gm,vector &b,vector &c)
{
  long i;

  nullm (gm);

  for (i=0;i<nne;i++){
    gm[0][i]=b[i];
    gm[1][i]=c[i];
  }
}



/**
   function computes conductivity %matrix of one medium

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param km - conductivity %matrix

   JK, 5.11.2001
*/
void trlineart::conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long ii;
  double jac,det,thick;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),b(3),c(3),t(nne);
  matrix gm(ncomp,dofe[ri][ci]),d;

  matrix n(1,dofe[ri][ci]);
  
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
  
  nullm (km);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];
  
  
  // gradient matrix
  grad_matrix (gm,b,c);

  //  matrix of conductivity of material
  reallocm (ncomp,ncomp,d);
  Tm->matcond (d,ii,ri,ci);

  thick = approx (areacoord,t);

  //  det is equal to double area of the element
  jac=thick*det/2.0;

  //  contribution to the stiffness matrix of the element
  bdbj (km.a,gm.a,d.a,jac,gm.m,gm.n);
  
  reallocm (1,ncomp,d);
  
  Tm->matcond2(d,ii,ri,ci);
  bf_matrix (n, areacoord);
  bdbjac(km, n, d, gm, jac);
  
  if ((Tt->elements[eid].transi[lcid]==3)||(Tt->elements[eid].transi[lcid]==4)){
    //  additional matrix due to transmission
    transmission_matrix (lcid,eid,ri,ci,km);
  }
}


/**
   function computes L %matrix

   L = \int_{\Omega} D B {\rm d} \Omega

   @param lcid - load case id
   @param eid - element id
   @param lm - L %matrix

   TKr, 16/08/2010
*/
void trlineart::l_matrix (long /*lcid*/,long eid,long ri,long ci,matrix &lm)
{
  long ii;
  double jac,det,thick;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),b(3),c(3),t(nne);
  matrix gm(ncomp,dofe[ri][ci]),d(ncomp,ncomp);

  matrix n(1,dofe[ri][ci]);
  
  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord2d (x,y,eid);
  Tc->give_thickness (eid,nodes,t);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  areacoord[0]=1.0/3.0;
  areacoord[1]=1.0/3.0;
  areacoord[2]=1.0/3.0;

  plsb (b.a,y.a,det);
  plsc (c.a,x.a,det);

  nullm (lm);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];
  
  // gradient matrix
  grad_matrix (gm,b,c);

  //  matrix of conductivity of material
  Tm->matcond (d,ii,ri,ci);

  thick = approx (areacoord,t);

  //  det is equal to double area of the element
  jac=thick*det/2.0;

  //  contribution to the L matrix of the element
  mxm(d,gm,lm);
  cmulm(jac,lm);
}


/**
   function computes L^T (L transposed) %matrix

   L^T = \int_{\Omega} B^T D {\rm d} \Omega

   @param lcid - load case id
   @param eid - element id
   @param lm - L^T %matrix

   TKr, 16/08/2010
*/
void trlineart::l_t_matrix (long /*lcid*/,long eid,long ri,long ci,matrix &lm)
{
  long ii;
  double jac,det,thick;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),b(3),c(3),t(nne);
  matrix gm(ncomp,dofe[ri][ci]),d(ncomp,ncomp);

  matrix n(1,dofe[ri][ci]);
  
  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord2d (x,y,eid);
  Tc->give_thickness (eid,nodes,t);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  areacoord[0]=1.0/3.0;
  areacoord[1]=1.0/3.0;
  areacoord[2]=1.0/3.0;

  plsb (b.a,y.a,det);
  plsc (c.a,x.a,det);

  nullm (lm);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];
  
  // gradient matrix
  grad_matrix (gm,b,c);

  //  matrix of conductivity of material
  Tm->matcond (d,ii,ri,ci);

  thick = approx (areacoord,t);

  //  det is equal to double area of the element
  jac=thick*det/2.0;

  //  contribution to the L matrix of the element
  mtxm(gm,d,lm);
  cmulm(jac,lm);
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
void trlineart::quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci)
{
  long i;
  double jac,det,thick;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),w(intordcm[ri][ci]),gp1(intordcm[ri][ci]),gp2(intordcm[ri][ci]),t(nne),v(dofe[ri][ci]);
  matrix n(1,dofe[ri][ci]),nm(dofe[ri][ci],dofe[ri][ci]);
  
  Tt->give_elemnodes (eid,nodes);

  Tc->give_thickness (eid,nodes,t);
  
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
    
    thick = approx (areacoord,t);
    
    jac=w[i]*thick*det;
    
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
void trlineart::internal_fluxes (long lcid,long eid,vector &ifl)
{
  long ipp;
  double det,thick;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),b(3),c(3),t(nne),fl(ncomp),contr(dofe[lcid][lcid]);
  matrix gm(ncomp,dofe[lcid][lcid]);
  
  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord2d (x,y,eid);
  Tc->give_thickness (eid,nodes,t);
  
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

  thick = approx (areacoord,t);
  
  Tm->computenlfluxes (lcid,ipp);
  
  Tm->givefluxes (lcid,ipp,fl);

  // gradient matrix
  grad_matrix (gm,b,c);
  
  mtxv (gm,fl,contr);
  
  cmulv (det/2.0*thick,contr);
  
  addv (contr,ifl,ifl);
}

/**
   function computes advection %matrix of one medium

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param hm - advection %matrix

   JK, 7.9.2016
*/
void trlineart::advection_matrix (long /*lcid*/,long eid,long ri,long ci,matrix &hm)
{
  long ii;
  double jac,det,thick;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),b(3),c(3),t(nne);
  matrix gm(ncomp,dofe[ri][ci]),n(1,dofe[ri][ci]),v(1,ncomp);
  
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

  nullm (hm);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];
  
  
  //  matrix of approximation functions
  bf_matrix (n,areacoord);

  // gradient matrix
  grad_matrix (gm,b,c);

  //  vector of velocity
  Tm->matcond2 (v,ii,ri,ci);

  //  thickness
  thick = approx (areacoord,t);

  //  det is equal to double area of the element
  jac=thick*det/2.0;

  //  contribution to the advection matrix of the element
  bdbjac (hm,n,v,gm,jac);
  
}


/**
   function assembles resulting element conductivity %matrix

   @param eid - element id
   @param lcid - load case id
   @param km - resulting conductivity %matrix of one element

   JK, 6.1.2002
*/
void trlineart::res_conductivity_matrix (long eid,long /*lcid*/,matrix &km)
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
   function assembles resulting element L %matrix

   L^T = \int_{\Omega} B^T D {\rm d} \Omega

   @param eid - element id
   @param lcid - load case id
   @param lm - resulting L %matrix of one element

   TKr, 16/08/2010
*/
void trlineart::res_l_matrix (long eid,long /*lcid*/,matrix &lm)
{
  long i,j,*rcn,*ccn;
  matrix lkm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      rcn = new long [2];
      ccn = new long [dofe[i][j]];
      if (i==0){
	rcn[0]=1;
	rcn[1]=2;
      }
      if (i==1){
	rcn[0]=3;
	rcn[1]=4;
      }
      reallocm (ncomp,dofe[i][j],lkm);
      l_matrix (i,eid,i,j,lkm);
      codnum (ccn,j);
      
      mat_localize (lm,lkm,rcn,ccn);
      delete [] rcn;  delete [] ccn;
    }
  }
}



/**
   function assembles resulting element L^T (L transposed) %matrix

   L^T = \int_{\Omega} B^T D {\rm d} \Omega

   @param eid - element id
   @param lcid - load case id
   @param lm - resulting L^T %matrix of one element

   TKr, 16/08/2010
*/
void trlineart::res_l_t_matrix (long eid,long /*lcid*/,matrix &lm)
{
  long i,j,*rcn,*ccn;
  matrix lkm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      rcn = new long [dofe[i][j]];
      ccn = new long [2];
      if (j==0){
	ccn[0]=1;
	ccn[1]=2;
      }
      if (j==1){
	ccn[0]=3;
	ccn[1]=4;
      }
      reallocm (dofe[i][j],ncomp,lkm);
      l_t_matrix (i,eid,i,j,lkm);
      codnum (rcn,i);
      
      mat_localize (lm,lkm,rcn,ccn);
      delete [] rcn;  delete [] ccn;
    }
  }
}




/**
   function assembles average D %matrix

   @param eid - element id
   @param lm - resulting D %matrix of one element

   TKr, 16/08/2010
*/
void trlineart::averd_matrix (long eid,matrix &lm)
{
  long i,j,ii,*rcn,*ccn;
  matrix d(ncomp,ncomp);
  
  rcn = new long [ncomp];
  ccn = new long [ncomp];

  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      
      if (Tp->savemode==0)
	ii=Tt->elements[eid].ipp[i][j];
      if (Tp->savemode==1)
	ii=Tt->elements[eid].ipp[0][0];
      
      //  matrix of conductivity of material
      Tm->matcond (d,ii,i,j);
      
      if (i==0){
	rcn[0]=1;
	rcn[1]=2;
      }
      if (i==1){
	rcn[0]=3;
	rcn[1]=4;
      }

      if (j==0){
	ccn[0]=1;
	ccn[1]=2;
      }
      if (j==1){
	ccn[0]=3;
	ccn[1]=4;
      }

      mat_localize (lm,d,rcn,ccn);
    }
  }
  
  delete [] rcn;  delete [] ccn;
}





/**
   function assembles average C %matrix

   @param eid - element id
   @param lm - resulting C %matrix of one element

   TKr, 16/08/2010
*/
void trlineart::averc_matrix (long eid,matrix &lm)
{
  long i,j,ii;
  double c,rho;
  ivector nodes(nne);
  vector dens(nne),areacoord(3);
  
  Tt->give_elemnodes (eid,nodes);
  Tc->give_density (eid,nodes,dens);
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      
      if (Tp->savemode==0)
	ii=Tt->elements[eid].ipp[i][j];
      if (Tp->savemode==1)
	ii=Tt->elements[eid].ipp[0][0];
      
      areacoord[0]=0.0;
      areacoord[1]=0.0;
      areacoord[2]=1.0-areacoord[0]-areacoord[1];

      //  coefficient of capacity of material
      c = Tm->capcoeff (ii,i,j);
      rho = approx (areacoord,dens);
      
      lm[i][j] = rho*c*Tm->ip[ii].av[j];
    }
  }
}


/**
   function assembles area of one element

   @param eid - element id

   TKr, 07/09/2010
*/
double trlineart::elem_area (long eid)
{
  double det,area;
  vector x(nne),y(nne);
  
  Tt->give_node_coord2d (x,y,eid);
  
  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  area = det/2.0;

  return area;
}





/**
   function assembles resulting element capacity %matrix
   
   @param eid - element id
   @param cm - resulting capacity %matrix of one element

   JK, 6.1.2002
*/
void trlineart::res_capacity_matrix (long eid,matrix &cm)
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
void trlineart::res_convection_vector (vector &f,long lcid,long eid,long leid)
{
  long *cn;
  vector lf;

  //  transi[lcid]==2 - element contains boundary with prescribed flux
  //  transi[lcid]==4 - element contains boundary with prescribed transmission and prescribed flux
  if ((Tt->elements[eid].transi[lcid]==2)||(Tt->elements[eid].transi[lcid]==4)){
    
    //  array for code numbers
    cn = new long [dofe[lcid][lcid]];
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
   @param leid - element id
   
   JK, 6.1.2002
*/
void trlineart::res_transmission_vector (vector &f,long lcid,long eid,long leid)
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
      nullv (lf);
      
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
void trlineart::res_quantity_source_vector (vector &sv,vector &nodval,long lcid,long eid)
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
void trlineart::res_internal_fluxes (long eid,vector &elemif)
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
  delete [] cn;
}

/**
   function assembles resulting element advection %matrix

   @param eid - element id
   @param lcid - load case id
   @param km - resulting advection %matrix of one element

   JK, 7.9.2016
*/
void trlineart::res_advection_matrix (long eid,long /*lcid*/,matrix &km)
{
  long i,j,*rcn,*ccn;
  matrix lkm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      rcn = new long [dofe[i][j]];
      ccn = new long [dofe[i][j]];
      reallocm (dofe[i][j],dofe[i][j],lkm);
      
      advection_matrix (i,eid,i,j,lkm);

      codnum (rcn,i);
      codnum (ccn,j);
      mat_localize (km,lkm,rcn,ccn);
      delete [] rcn;  delete [] ccn;
    }
  }
}

/**
   function computes element quantity integral
   
   @param eid - element id
   @param nodval - %vector of quantity nodal values

   @retval f - element quantity integral

   TKr, 30.1.2004
*/
double trlineart::total_integral (long eid,vector &nodval)
{
  long i;
  double jac,det,thick,value,f;
  ivector nodes(nne);
  vector x(nne),y(nne),areacoord(3),w(3),gp1(3),gp2(3),t(nne);

  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord2d (x,y,eid);
  Tc->give_thickness (eid,nodes,t);

  //  det is equal to double area of the element
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

  gauss_points_tr (gp1.a,gp2.a,w.a,3);

  f = 0.0;

  for (i=0;i<3;i++){
    areacoord[0]=gp1[i];
    areacoord[1]=gp2[i];
    areacoord[2]=1.0-gp1[i]-gp2[i];

    thick = approx (areacoord,t);
    value = approx (areacoord,nodval);

    jac=w[i]*thick*value*det;

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
void trlineart::res_boundary_flux (vector &f,long lcid,long eid,long leid)
{
  long *cn,i;
  vector lf;

  cn = new long [dofe[lcid][lcid]];
  codnum (cn,lcid);
  
  //coupled b.c.
  for (i=0;i<ntm;i++){    
    reallocv (dofe[lcid][lcid],lf);
    
    if ((Tt->elements[eid].transi[i]==1) || (Tt->elements[eid].transi[i]==2)){
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
void trlineart::volume_rhs_vector (long /*lcid*/,long eid,long ri,long ci,vector &vrhs)
{
  long i,ii;
  double jac,det,thick;
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

  thick = approx (areacoord,t);

  //  det is equal to double area of the element
  jac=thick*det/2.0;
  
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
void trlineart::volume_rhs_vector2 (long /*lcid*/,long eid,long ri,long ci,vector &vrhs)
{
  long i,k,ii;
  double jac,det,thick,c;
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
    
    thick = approx (areacoord,t);
    
    c=Tm->volume_rhs2 (ii,ri,ci);
    
    jac=w[i]*thick*det*c;
    
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
void trlineart::res_volume_rhs_vector (vector &f,long eid,long /*lcid*/)
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
void trlineart::res_volume_rhs_vector2 (vector &f,long eid,long /*lcid*/)
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
void trlineart::intpointflux (long eid)
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
   
   23. 3. 2017, JK
*/
void trlineart::nod_grads_ip (long eid)
{
  long i,j,k,ipp;
  double area;
  ivector nod(nne);
  vector grad(ncomp);

  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  
  //  node numbers of the element
  Tt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    for (k=0;k<Tp->ntm;k++){
      //  gradients at the closest integration point
      Tm->givegrad (k,ipp,grad);
      
      //  storage of gradients to the node
      if (Tp->gradaver==0 || Tp->gradaver==1){
	//  gradients will be averaged by the number of contributions
	j=nod[i];
	Tt->nodes[j].storegrad (k,grad);
      }
      if (Tp->gradaver==2){
	//  gradients will be averaged by the volume
	area = element_area (eid);
	j=nod[i];
	Tt->nodes[j].storegrad (k,area,grad);
      }
    }
  }
  
}


/**
   function computes fluxes in nodes of element

   @param eid - element id
   
   23. 3. 2017, JK
*/
void trlineart::nod_fluxes_ip (long eid)
{
  long i,j,k,ipp;
  double area;
  ivector nod(nne);
  vector flux(ncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  
  //  node numbers of the element
  Tt->give_elemnodes (eid,nod);
  
  intpointflux (eid);

  for (i=0;i<nne;i++){
    for (k=0;k<Tp->ntm;k++){
      //  fluxes at the closest integration point
      Tm->givefluxes (k,ipp,flux);
      
      //  storage of fluxes to the node
      if (Tp->fluxaver==0 || Tp->fluxaver==1){
	//  fluxes will not be averaged or will be averaged by the number of contributions
	j=nod[i];
	Tt->nodes[j].storeflux (k,flux);
      }
      if (Tp->fluxaver==2){
	//  fluxes will be averaged by the volume
	area = element_area (eid);
	j=nod[i];
	Tt->nodes[j].storeflux (k,area,flux);
      }
    }
  }
  
}


/**
   function computes other values in nodes of element

   @param eid - element id
   
   23. 3. 2017, JK
*/
void trlineart::nod_eqother_ip (long eid)
{
  long i,j,k,ipp,ncompo;
  double area;
  ivector nod(nne);
  vector eqother;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  
  //  node numbers of the element
  Tt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    ncompo = Tm->givencompeqother (ipp,0);
    reallocv (ncompo,eqother);
    Tm->giveeqother (ipp,0,ncompo,eqother.a);
    
    //  storage of eqother components to the node
    j=nod[i];
    for (k=0; k<ncompo; k++){
      if (Tp->eqotheraver==0 || Tp->eqotheraver==1)
	Tt->nodes[j].storeeqother (k,eqother(k));
      if (Tp->eqotheraver==2){
	area = element_area (eid);
	Tt->nodes[j].storeeqother (k,area,eqother(k));
      }
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
void trlineart::nod_others_comp (long /*lcid*/,long eid,long ri,long ci)
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
   function computes nodal fluxes from boundary values
   
   \int_{Gamma_2} N^T N dGamma * nodal_flux_values

   @param v - array of nodal fluxes
   @param lcid - load case id
   @param eid - element id
   @param leid - id of loaded element
   @param ri,ci - row and column indices
   
   JK, 19.8.2004
*/
void trlineart::convection_vector (vector &v,long lcid,long eid,long leid,long ri,long ci)
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
  //  thicknesses at nodes
  Tc->give_thickness (eid,nodes,t);
  //  node coordinates on element
  Tt->give_node_coord2d (x,y,eid);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  
  //  indicators of boundary conditions
  bc = new bocontypet [ned];
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
    
    if (bc[i]==2 || bc[i]==3 || bc[i]==4 || bc[i]==5){
      //  nodal values on actual edge
      edgenodeval (i,nodval,list);
      
      edgenodeval (i,coef,trc);
      
      nullm (km);

      //  matrix obtained from integration over edge
      edge_integral (i,x,y,intordkm[ri][ci],gp,w,t,coef,km);
      
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
void trlineart::transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,ipp,leid;
  bocontypet *bc;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]);
  vector trc(nned*ned),coeff(nne),t(nne);
  
  //  element nodes
  Tt->give_elemnodes (eid,nodes);
  //  thicknesses at nodes
  Tc->give_thickness (eid,nodes,t);
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
  bc = new bocontypet [ned];
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
    
    if (bc[i]==4 || bc[i]==5 || bc[i]>10){
      //  transformation of coefficients
      transf_coeff (i,coeff,trc,eid,ri,ci,ipp,bc);
      
      //  matrix obtained from integration over edges
      edge_integral (i,x,y,intordkm[ri][ci],gp,w,t,coeff,km);
    }
  }
  
  delete [] bc;
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
void trlineart::transmission_vector (vector &v,long lcid,long eid,long leid,long cid)
{
  long i,ipp;
  bocontypet *bc;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),w(intordkm[lcid][cid]),gp(intordkm[lcid][cid]),t(nne);
  vector list(nned*ned),trc(nned*ned),trr(nned*ned),nodval(nne),av(nne),coef(nne);
  matrix km(nne,nne);

  nullv (v);
  
  //  element nodes
  Tt->give_elemnodes (eid,nodes);
  //  thicknesses at nodes
  Tc->give_thickness (eid,nodes,t);
  //  node coordinates on element
  Tt->give_node_coord2d (x,y,eid);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[lcid][cid]);
  
  //  indicators of boundary conditions
  bc = new bocontypet [ned];
  Tb->lc[lcid].elemload[leid].give_bc (bc);
  //  nodal values of exterior variable
  Tb->lc[lcid].elemload[leid].give_external_nodval (lcid,cid,list);
  //  transmission coefficients
  Tb->lc[lcid].elemload[leid].give_trc (Tp->time,Tb->lc[lcid].nodval,lcid,cid,trc);
  //  transmission/radiation coefficients
  Tb->lc[lcid].elemload[leid].give_trr (Tp->time,Tb->lc[lcid].nodval,trr);
  
  //  integration point number
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[lcid][cid];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  for (i=0;i<ned;i++){
    
    if (bc[i]==4 || bc[i]==5 || bc[i]>10){
      //  transformation of nodal values
      transf_val (i,nodval,list,trc,trr,eid,lcid,cid,ipp,bc);
      
      edgenodeval (i,coef,trc);
      
      nullm (km);
      //
      edge_integral (i,x,y,intordkm[lcid][cid],gp,w,t,coef,km);
      
      mxv (km,nodval,av);

      addv (v,av,v);
    }
  }
  
  delete [] bc;
}


/**
   function computes contributions to the boundary flux from transmission %vector
   
   \int_{Gamma_3} N^T c_{tr} N dGamma * nodal_external_value

   @param tmv - %vector of boundary fluxes
   @param lcid - load case id
   @param eid - element id
   @param leid - loaded element id
   @param ri,ci - row and column indices
   
   TKr, 28.2.2004
*/
void trlineart::boundary_flux (vector &v,long lcid,long eid,long leid,long ri,long ci)
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
  //  thicknesses at nodes
  Tc->give_thickness (eid,nodes,t);
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
    
    if (bc[i]==4 || bc[i]==5 || bc[i]>10){
      //  transformation of nodal values
      transf_flux (i,nodval,list,trc,trr,eid,ri,ci,ipp,bc);
      
      edgenodeval (i,coef,trc);
      
      nullm (km);
      //
      edge_integral (i,x,y,intordkm[ri][ci],gp,w,t,coef,km);
      
      mxv (km,nodval,av);
      
      addv (v,av,v);
    }
  }
  
  delete [] bc;
}


/**
   function integrates N^T c N over edges
   
   @param edg - edge id (number of edge)
   @param x, y - coordinates of element nodes
   @param intord - order of numerical integration
   @param gp, w  - coordinates and weights of integration points
   @param t - nodal thicknesses
   @param coef - array of nodal values of coefficient
   @param km - output %matrix
   
   JK
*/
void trlineart::edge_integral (long edg,vector &x,vector &y,long intord,vector &gp,vector &w,
			       vector &t,vector &coef,matrix &km)
{
  long i;
  double jac,ipval,thick;
  vector av(nne),areacoord(3);
  matrix n(1,nne);
  
  if (edg==0){
    for (i=0;i<intord;i++){
      areacoord[0]=(1.0-gp[i])/2.0;  areacoord[1]=(1.0+gp[i])/2.0;
      areacoord[2]=1.0-areacoord[0]-areacoord[1];
      
      bf_matrix (n,areacoord);
      
      jac1d_2d (jac,x,y,gp[i],edg);
      thick = approx (areacoord,t);
      ipval = approx (areacoord,coef);
      
      jac*=w[i]*thick*ipval;
      
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

  if (edg==1){
    for (i=0;i<intord;i++){
      areacoord[0]=0.0;  areacoord[1]=(1.0-gp[i])/2.0;
      areacoord[2]=1.0-areacoord[0]-areacoord[1];
      
      bf_matrix (n,areacoord);
      
      jac1d_2d (jac,x,y,gp[i],edg);
      thick = approx (areacoord,t);
      ipval = approx (areacoord,coef);

      jac*=w[i]*thick*ipval;
      
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

  if (edg==2){
    for (i=0;i<intord;i++){
      areacoord[0]=(1.0+gp[i])/2.0;  areacoord[1]=0.0;
      areacoord[2]=1.0-areacoord[0]-areacoord[1];
      
      bf_matrix (n,areacoord);
      
      jac1d_2d (jac,x,y,gp[i],2);
      thick = approx (areacoord,t);
      ipval = approx (areacoord,coef);

      jac*=w[i]*thick*ipval;
      
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

}

void trlineart::transf_flux (long edg,vector &coeff,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
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


void trlineart::transf_coeff (long edg,vector &coeff,vector &list,long eid,long ri,long ci,long ipp,bocontypet *bc)
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
void trlineart::transf_val (long edg,vector &nodval,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
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
void trlineart::edgenodeval (long edg,vector &nodval,vector &list)
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
   function selects values and gradients from the global level
   
   @param eid - element id
   @param counter - actual position in the array buff
   @param buff - array containing selected components
   
   JK, 25.3.2011
*/
void trlineart::higher_to_lower_level (long eid,long *counter,double *buff)
{
  long i,ipp;
  vector gr(2);
  
  //  id of integration point
  ipp=Tt->elements[eid].ipp[0][0];
  
  //  loop over the number of transported materials
  for (i=0;i<ntm;i++){
    
    //  value
    buff[counter[0]]=Tm->ip[ipp].av[i];
    counter[0]++;
    
    //  components of gradient
    Tm->givegrad (i,ipp,gr);
    buff[counter[0]]=gr[0];
    counter[0]++;
    buff[counter[0]]=gr[1];
    counter[0]++;
    
  }
  
}


/**
   Function
   1. computes "e2" - square of "energy" norm of error of solution at element
      e2 = integral{ (der_fine - der)T * D * (der_fine - der) }
      der_fine = recovered gradients, values in nodes are defined by "rderfull" and interpolated by base functions over element
      der = gradiensts obtained by FEM
      D_inv = inverse stiffness matrix of material
   2. computes "u2" - square of energy norm of gradient at element
      u2 = integral{ gradT * D * grad }
      grad = gradient obtained by FEM
      D = matrix of conductivity
   3. computes area of element and returns
   4. computes "sizel" (characteristic size of element)
   
   @param eid - element id
   @param e2 - empty(returned) value; 
   @param u2 - empty(returned) value; 
   @param sizel - empty(returned) value; 
   @param rderfull - 1Darray of gradients in nodes; dimmension is gt->nn x ncomp
   
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
double trlineart :: compute_error (long eid, vector *rderfull, int mattid, double &e2, double &u2, double &sizel)
{
  //   ????????? bude to pro vice medii? ????
  //   pak bude chyba pro kazde medium asi jina, ne?
  //   pak bych mel mit ve volani fce argument ri, ci nebo nejak jinak rozpoznat medium ...
  
  long i, ipp;
  double det, contr;
  ivector nodes(nne);
  vector x(nne), y(nne), areacoord(3), t(nne);
  vector der(ncomp), der_fine(ncomp), der_err(ncomp);
  matrix d(ncomp,ncomp);
  
  Tt->give_elemnodes (eid, nodes);
  Tt->give_node_coord2d (x, y, eid);
  Tc->give_thickness (eid, nodes, t);
  
  // kdyz to nebude savemode, tak misto 0,0 tam bude ri ci podle media
  if (Tp->savemode==1)  ipp = Tt->elements[eid].ipp[0][0];
  else { print_err("save mode is required at element %ld", __FILE__, __LINE__, __func__, eid+1);  exit (1); }
  
  det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
  // matrix of conductivity of material
  Tm->matcond (d, ipp, mattid, mattid);
  
  // compute u2
  fillv (1.0/3.0, areacoord);
  Tm->givegrad (mattid, ipp, der);
  vxmxv(der, d, contr);
  u2 = contr * 0.5*det*approx(areacoord,t); // contr * jacobian
  
  
  // compute e2
  long intord = 3; // =intordcm
  vector gp1(intord), gp2(intord), w(intord);
  
  gauss_points_tr (gp1.a, gp2.a, w.a, intord);
  
  e2 = 0;
  for (i=0; i<intord; i++) {
    areacoord[0]=gp1[i];
    areacoord[1]=gp2[i];
    areacoord[2]=1.0-gp1[i]-gp2[i];
    
    // der above
    give_der_star (areacoord, rderfull, nodes, der_fine, Tt->nn);
    subv (der_fine, der, der_err);
    
    vxmxv (der_err, d, contr);
    
    e2 += contr * w[i]*det*approx(areacoord,t);
  }
  
  sizel = sqrt(1.1547005384*det);
  return  det/2.0;
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
long trlineart::ipcoord (long eid, long ipp, long ri, long ci, vector &coord)
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
long trlineart::ipncoord (long eid, long ipp, vector &ncoord)
{
  long i, ri, ci, ii;
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
