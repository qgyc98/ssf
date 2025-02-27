/*
  File:     quadlinaxisym.cpp
  Author:   Tomas Krejci, 31.3.2002
  Purpose:  twodimensional quadrilateral element with linear approximation functions for axisymmetric problem
*/

#include "globalt.h"
#include "quadlinaxisym.h"
#include "genfile.h"
#include "globmatt.h"

quadlinaxisym::quadlinaxisym (void)
{
  long i;

  //  number of nodes on element
  nne=4;
  //  number of edges
  ned=4;
  //  number of nodes on one edge
  nned=2;
  //  geometric problem dimension (2D)
  ncomp=2;

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
    ordering[0][0]=1;  ordering[0][1]=2;  ordering[0][2]=3;  ordering[0][3]=4;
    dofe[0][0]=4;  intordkm[0][0]=2;  intordcm[0][0]=2;  nip[0][0]=8;
    ndofe=4;  napfun=1;
    break;
  }
  case twomediacoup:{
    ordering[0][0]=1;  ordering[0][1]=3;  ordering[0][2]=5;  ordering[0][3]=7;
    ordering[1][0]=2;  ordering[1][1]=4;  ordering[1][2]=6;  ordering[1][3]=8;

    intordkm[0][0]=2;  intordkm[0][1]=2;  intordkm[1][0]=2;  intordkm[1][1]=2;
    intordcm[0][0]=2;  intordcm[0][1]=2;  intordcm[1][0]=2;  intordcm[1][1]=2;
    
    if (Tp->savemode==0){
      nip[0][0]=8;       nip[0][1]=8;       nip[1][0]=8;       nip[1][1]=8;
    }
    if (Tp->savemode==1){
      nip[0][0]=8;       nip[0][1]=0;       nip[1][0]=0;       nip[1][1]=0;
    }

    dofe[0][0]=4;  dofe[0][1]=4;  dofe[1][0]=4;  dofe[1][1]=4;
    ndofe=8;  napfun=2;
    break;
  }
  case threemediacoup:{
    ordering[0][0]=1;   ordering[0][1] =4;  ordering[0][2] =7;  ordering[0][3] =10;
    ordering[1][0]=2;   ordering[1][1] =5;  ordering[1][2] =8;  ordering[1][3] =11;
    ordering[2][0]=3;   ordering[2][1] =6;  ordering[2][2] =9;  ordering[2][3] =12;

    intordkm[0][0]=2;  intordkm[0][1]=2;  intordkm[0][2]=2;
    intordkm[1][0]=2;  intordkm[1][1]=2;  intordkm[1][2]=2;
    intordkm[2][0]=2;  intordkm[2][1]=2;  intordkm[2][2]=2;

    intordcm[0][0]=2;  intordcm[0][1]=2;  intordcm[0][2]=2;
    intordcm[1][0]=2;  intordcm[1][1]=2;  intordcm[1][2]=2;
    intordcm[2][0]=2;  intordcm[2][1]=2;  intordcm[2][2]=2;
    
    if (Tp->savemode==0){
      nip[0][0]=8;  nip[0][1]=8;  nip[0][2]=8;
      nip[1][0]=8;  nip[1][1]=8;  nip[1][2]=8;
      nip[2][0]=8;  nip[2][1]=8;  nip[2][2]=8;
    }
    if (Tp->savemode==1){
      nip[0][0]=8;  nip[0][1]=0;  nip[0][2]=0;
      nip[1][0]=0;  nip[1][1]=0;  nip[1][2]=0;
      nip[2][0]=0;  nip[2][1]=0;  nip[2][2]=0;
    }
    
    dofe[0][0]=4;  dofe[0][1]=4;  dofe[0][2]=4;
    dofe[1][0]=4;  dofe[1][1]=4;  dofe[1][2]=4;
    dofe[2][0]=4;  dofe[2][1]=4;  dofe[2][2]=4;

    ndofe=12;  napfun=3;
    break;
  }
  default:{
    print_err("unknown type of transported matter is required",__FILE__,__LINE__,__func__);
  }
  }
}

quadlinaxisym::~quadlinaxisym (void)
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

void quadlinaxisym::codnum (long *cn,long ri)
{
  long i;
  for (i=0;i<nne;i++){
    cn[i]=ordering[ri][i];
  }
}

/**
   function computes area of the element
   
   @param eid - element id
   
   23. 3. 2017, JK
*/
double quadlinaxisym::element_area (long eid)
{
  double area;
  vector x(nne),y(nne);

  //  node coordinates
  Tt->give_node_coord2d (x,y,eid);

  area = (x[1]-x[0])*(y[3]-y[0])-(x[3]-x[0])*(y[1]-y[0]);
  area += (x[3]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[3]-y[2]);
  area /=2.0;
  
  return area;
}

/**
   function approximates function defined by nodal values
   
   @param xi,eta - coordinates on element
   @param nodval - %vector of nodal values

   JK, 25.9.2001
*/
double quadlinaxisym::approx (double xi,double eta,vector &nodval)
{
  double f;
  vector bf(ASTCKVEC(nne));
  
  bf_lin_4_2d (bf.a,xi,eta);
  
  scprd (bf,nodval,f);

  return f;
}



/**
   function assembles coordinates of integration points in block [ri][ci]
   
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param ipcoord - array containing coordinates of integration points
   
   JK, 8.5.2002 + TKr 14.11.2003
*/
void quadlinaxisym::ipcoordblock (long eid,long ri,long ci,double **coord)
{
  long i,j,k;
  double xi,eta;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w(ASTCKVEC(intordkm[ri][ci])),gp(ASTCKVEC(intordkm[ri][ci]));
  
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  Tt->give_node_coord2d (x,y,eid);
  
  k=0;
  for (i=0;i<intordkm[ri][ci];i++){
    xi=gp[i];
    for (j=0;j<intordkm[ri][ci];j++){
      eta=gp[j];
      
      coord[k][0]=approx (xi,eta,x);
      coord[k][1]=approx (xi,eta,y);
      coord[k++][2]=0.0;
    }
  }
}

/**
   function computes values in integration points from nodal values

   eid - element id

   JK, 25.9.2001
*/
void quadlinaxisym::intpointval (long eid)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,val;
  vector r(ASTCKVEC(ndofe)),t(ASTCKVEC(nne)),gp,w;

  elemvalues(eid, r);
  
  for (k=0;k<Tp->ntm;k++){

    //  nodal values of a single variable
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
	  for (j=0;j<intordkm[ii][jj];j++){
	    eta=gp[j];
	    
	    val = approx (xi,eta,t);
	    Tm->ip[ipp].av[k]=val;
	    ipp++;
	  }
	}

	reallocv (RSTCKVEC(intordcm[ii][jj],gp));
	reallocv (RSTCKVEC(intordcm[ii][jj],w));
	gauss_points (gp.a,w.a,intordcm[ii][jj]);
	
	for (i=0;i<intordcm[ii][jj];i++){
	  xi=gp[i];
	  for (j=0;j<intordcm[ii][jj];j++){
	    eta=gp[j];
	    
	    val = approx (xi,eta,t);
	    Tm->ip[ipp].av[k]=val;
	    ipp++;
	  }
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

   25/05/2016 TKr according to JK
*/
void quadlinaxisym::intpointval (long eid,vector &nodval,vector &ipval)
{
  long i,j,kk,ii,jj,ipp;
  double xi,eta;
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
        for (j=0;j<intordkm[ii][jj];j++){
          eta=gp[j];
	    
          ipval[kk] = approx (xi,eta,nodval);

          ipp++;
          kk++;
        }
      }
      
      
      //  interpolation to integration points of capacity matrix
      reallocv (RSTCKVEC(intordcm[ii][jj],gp));
      reallocv (RSTCKVEC(intordcm[ii][jj],w));
      gauss_points (gp.a,w.a,intordcm[ii][jj]);
	
      for (i=0;i<intordcm[ii][jj];i++){
        xi=gp[i];
        for (j=0;j<intordcm[ii][jj];j++){
          eta=gp[j];
	    
          ipval[kk] = approx (xi,eta,nodval);

          ipp++;
          kk++; // TKo ???!!!
        }
      }
	
      if (Tp->savemode==1)  break;
    }
    if (Tp->savemode==1)  break;
  }
}



/**
  The function computes initial values in integration points from initial nodal values.
  (function interpolates the initial nodal values to integration points)

  @param eid - element id
 
  @return The function does not return anything but stores computed values at int. points of the given element.

  Created by TKo, 4.7.2018
*/
void quadlinaxisym::initintpointval (long eid)
{
  long i,j,k,ii,jj,ipp,ndofn,cndofn;
  double xi,eta,val;
  ivector enod(ASTCKIVEC(nne));
  vector rr, r(ASTCKVEC(ndofe)), t(ASTCKVEC(nne)), gp, w;

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
	
	reallocv (RSTCKVEC(intordkm[ii][jj],gp));
	reallocv (RSTCKVEC(intordkm[ii][jj],w));
	gauss_points (gp.a,w.a,intordkm[ii][jj]);
	
	ipp=Tt->elements[eid].ipp[ii][jj];
	for (i=0;i<intordkm[ii][jj];i++){
	  xi=gp[i];
	  for (j=0;j<intordkm[ii][jj];j++){
	    eta=gp[j];
	    
	    val = approx (xi,eta,t);
	    Tm->ip[ipp].av[k]=val;
	    ipp++;
	  }
	}

	reallocv (RSTCKVEC(intordcm[ii][jj],gp));
	reallocv (RSTCKVEC(intordcm[ii][jj],w));
	gauss_points (gp.a,w.a,intordcm[ii][jj]);
	
	for (i=0;i<intordcm[ii][jj];i++){
	  xi=gp[i];
	  for (j=0;j<intordcm[ii][jj];j++){
	    eta=gp[j];
	    
	    val = approx (xi,eta,t);
	    Tm->ip[ipp].av[k]=val;
	    ipp++;
	  }
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
void quadlinaxisym::intpointgrad (long eid)
{
  long i,j,ii,jj,k,ipp;
  double xi,eta,jac;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),r(ASTCKVEC(ndofe));
  vector gp, w, t(ASTCKVEC(nne)),grad(ASTCKVEC(ncomp));
  matrix gm(ASTCKMAT(ncomp,nne));
  
  Tt->give_node_coord2d (x,y,eid);
  elemvalues(eid, r);

  for (k=0;k<Tp->ntm;k++){

    //  nodal values of a single variable
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
	  for (j=0;j<intordkm[ii][jj];j++){
	    eta=gp[j];

	    //  matrix of gradients
	    grad_matrix (gm,x,y,xi,eta,jac);
	    mxv (gm,t,grad);
	    Tm->storegrad (k,ipp,grad);
	    ipp++;
	  }
	}
	
	reallocv (RSTCKVEC(intordcm[ii][jj],gp));
	reallocv (RSTCKVEC(intordcm[ii][jj],w));
	gauss_points (gp.a,w.a,intordcm[ii][jj]);

	for (i=0;i<intordcm[ii][jj];i++){
	  xi=gp[i];
	  for (j=0;j<intordcm[ii][jj];j++){
	    eta=gp[j];

	    //  matrix of gradients
	    grad_matrix (gm,x,y,xi,eta,jac);
	    mxv (gm,t,grad);
	    Tm->storegrad (k,ipp,grad);
	    ipp++;
	  }
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
void quadlinaxisym::intpointother (long eid)
{
  long i, j, k, ii, jj, ipp, ncompo, nodid;
  double xi, eta, val;
  ivector nodes(ASTCKIVEC(nne));
  vector r, t(ASTCKVEC(nne)), gp, w;
  
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
	
	reallocv (RSTCKVEC(intordkm[ii][jj],gp));
	reallocv (RSTCKVEC(intordkm[ii][jj],w));
	gauss_points (gp.a,w.a,intordkm[ii][jj]);
	
	ipp=Tt->elements[eid].ipp[ii][jj];
	for (i=0;i<intordkm[ii][jj];i++){
	  xi=gp[i];
	  for (j=0;j<intordkm[ii][jj];j++){
	    eta=gp[j];
	    
	    val = approx (xi,eta,t);
	    Tm->ip[ipp].other[k]=val;
	    ipp++;
	  }
	}
	
	reallocv (RSTCKVEC(intordcm[ii][jj],gp));
	reallocv (RSTCKVEC(intordcm[ii][jj],w));
	gauss_points (gp.a,w.a,intordcm[ii][jj]);
	
	for (i=0;i<intordcm[ii][jj];i++){
	  xi=gp[i];
	  for (j=0;j<intordcm[ii][jj];j++){
	    eta=gp[j];
	    
	    val = approx (xi,eta,t);
	    Tm->ip[ipp].other[k]=val;
	    ipp++;
	  }
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
   @param xi,eta - natural coordinates

   JK, 25.9.2001
*/
void quadlinaxisym::bf_matrix (matrix &n,double xi,double eta)
{
  nullm (n);
  bf_lin_4_2d (n.a,xi,eta);
}

/**
   function returns one approximation function evaluated in required point

   @param f - value of approximation function
   @param xi,eta - natural coordinates
   @param i - number of approximation function

   JK, 1.2.2003
*/
void quadlinaxisym::give_approx_fun (double &f,double xi,double eta,long i)
{
  double n[4];
  bf_lin_4_2d (n,xi,eta);
  f=n[i];
}


/**
   function assembles gradient of matrix of base functions

   @param gm - gradient %matrix
   @param x,y - array containing node coordinates
   @param xi,eta - natural coorodinates
   @param jac - Jacobian

   JK, 25.9.2001
*/
void quadlinaxisym::grad_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac)
{
  long i;
  vector dx(ASTCKVEC(nne)),dy(ASTCKVEC(nne));

  dx_bf_lin_4_2d (dx.a,eta);
  dy_bf_lin_4_2d (dy.a,xi);

  derivatives_2d (dx,dy,jac,x,y,xi,eta);
 
  nullm (gm);

  for (i=0;i<nne;i++){
    gm[0][i]=dx[i];
    gm[1][i]=dy[i];
  }

}



/**
   function computes conductivity %matrix of 2D problems for one transported matter
   finite element with bilinear approximation functions

   @param lcid - number of load case
   @param eid - number of element
   @param ri,ci - row and column indeces of the computed block in the resulting %matrix
   @param km - conductivity %matrix

   JK, 25.9.2001
*/
void quadlinaxisym::conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,j,ii;
  double xinp,xi,eta,ww1,ww2,jac;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w(ASTCKVEC(intordkm[ri][ci])),gp(ASTCKVEC(intordkm[ri][ci]));
  matrix gm(ASTCKMAT(ncomp,dofe[ri][ci])),d;

  matrix n(ASTCKMAT(1,dofe[ri][ci]));
 
  Tt->give_elemnodes (eid,nodes);
  

  Tt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordkm[ri][ci]);

  nullm (km);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];

  for (i=0;i<intordkm[ri][ci];i++){
    xi=gp[i];  ww1=w[i];
    for (j=0;j<intordkm[ri][ci];j++){
      eta=gp[j];  ww2=w[j];

      //  matrix of gradients
      grad_matrix (gm,x,y,xi,eta,jac);

         
      if(jac < 0.0){
	//jac = fabs (jac);
	print_err("wrong numbering of nodes on element number %ld, negative volume! jac = %e", __FILE__, __LINE__, __func__, eid+1, jac);
	abort();
      }

      //  matrix of conductivity of the material
      reallocm(RSTCKMAT(ncomp,ncomp,d));
      Tm->matcond (d,ii,ri,ci);

      check_math_errel(eid);

      xinp = approx (xi,eta,x);

      jac*=xinp*ww1*ww2;

      //  contribution to the conductivity matrix of the element
      bdbj (km.a,gm.a,d.a,jac,gm.m,gm.n);
      
      reallocm(RSTCKMAT(1,ncomp,d));
      
      Tm->matcond2(d,ii,ri,ci);
      bf_matrix (n, xi, eta);
      bdbjac(km, n, d, gm, jac);

      check_math_errel(eid);
      
      ii++;  
    }
  }
  
  if ((Tt->elements[eid].transi[lcid]==3)||(Tt->elements[eid].transi[lcid]==4)){
    transmission_matrix (lcid,eid,ri,ci,km);
  }
  
}



/**
   function computes capacity %matrix of 2D problems for one transported matter
   finite element with bilinear approximation functions

   @param eid - number of element
   @param ri,ci - row and column indeces of the computed block in the resulting %matrix
   @param cm - capacity %matrix

   JK, 4.10.2001
*/
void quadlinaxisym::capacity_matrix (long eid,long ri,long ci,matrix &cm)
{
  long i,j,ii;
  double jac,xi,eta,w1,w2,xinp,rho,c;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w(ASTCKVEC(intordcm[ri][ci])),gp(ASTCKVEC(intordcm[ri][ci])),dens(ASTCKVEC(nne));
  matrix n(ASTCKMAT(1,dofe[ri][ci]));

  Tt->give_elemnodes (eid,nodes);
  
  Tc->give_density (eid,nodes,dens);

  Tt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordcm[ri][ci]);

  nullm (cm);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci]+intordkm[ri][ci]*intordkm[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0]+intordkm[0][0]*intordkm[0][0];

  for (i=0;i<intordcm[ri][ci];i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordcm[ri][ci];j++){
      eta=gp[j];  w2=w[j];
      jac_2d (jac,x,y,xi,eta);
      bf_matrix (n,xi,eta);

      check_math_errel(eid);

      c=Tm->capcoeff (ii,ri,ci);

      check_math_errel(eid);

      xinp = approx (xi,eta,x);
      rho = approx (xi,eta,dens);
      jac*=w1*w2*xinp*rho*c;

      nnj (cm.a,n.a,jac,n.m,n.n);
      ii++;
    }
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
void quadlinaxisym::quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci)
{
  long i,j;
  double jac,xi,eta,w1,w2,xinp;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w(ASTCKVEC(intordcm[ri][ci])),gp(ASTCKVEC(intordcm[ri][ci])),t(ASTCKVEC(nne)),dens(ASTCKVEC(nne)),v(ASTCKVEC(dofe[ri][ci]));
  matrix n(ASTCKMAT(1,dofe[ri][ci])),nm(ASTCKMAT(dofe[ri][ci],dofe[ri][ci]));

  Tt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordcm[ri][ci]);
  
  nullm (nm);

  for (i=0;i<intordcm[ri][ci];i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordcm[ri][ci];j++){
      eta=gp[j];  w2=w[j];
      jac_2d (jac,x,y,xi,eta);
      bf_matrix (n,xi,eta);

      xinp = approx (xi,eta,x);
      jac*=w1*w2*xinp;
      
      nnj (nm.a,n.a,jac,n.m,n.n);
    }
  }
  mxv (nm,nodval,v);  addv (sv,v,sv);

}






/**
   function computes internal fluxes for one transported matter
   finite element with linear approximation functions
   
   @param lcid - number of load case
   @param eid - number of element
   @param ifl - %vector of internal fluxes
   
   JK, 31.3.2002
*/
void quadlinaxisym::internal_fluxes (long lcid,long eid,vector &ifl)
{
  long i,j,ipp;
  double xi,eta,jac,xinp;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),t(ASTCKVEC(nne)),fl(ASTCKVEC(ncomp));
  vector w,gp,contr(ASTCKVEC(dofe[lcid][lcid]));
  matrix gm(ASTCKMAT(ncomp,dofe[lcid][lcid]));
  

  Tt->give_elemnodes (eid,nodes);
  
  Tt->give_node_coord2d (x,y,eid);
  
  reallocv (RSTCKVEC(intordkm[lcid][lcid],gp));
  reallocv (RSTCKVEC(intordkm[lcid][lcid],w));
  
  gauss_points (gp.a,w.a,intordkm[lcid][lcid]);

  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[lcid][lcid];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  for (i=0;i<intordkm[lcid][lcid];i++){
    xi=gp[i];
    for (j=0;j<intordkm[lcid][lcid];j++){
      eta=gp[j];
      xinp = approx (xi,eta,x);
      
      Tm->computenlfluxes (lcid,ipp);
      
      Tm->givefluxes (lcid,ipp,fl);

      grad_matrix (gm,x,y,xi,eta,jac);

      mtxv (gm,fl,contr);

      cmulv (jac*w[i]*w[j]*xinp,contr);
      
      addv (contr,ifl,ifl);
      
      ipp++;
    }
  }
  
}

/**
   function assembles resulting element conductivity matrix

   @param eid - element id
   @param lcid - load case id
   @param km - resulting conductivity %matrix of one element

   JK, 6.1.2002
*/
void quadlinaxisym::res_conductivity_matrix (long eid,long /*lcid*/,matrix &km)
{
  long i,j,*rcn,*ccn;
  matrix lkm;

  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      rcn = new long [dofe[i][j]];
      ccn = new long [dofe[i][j]];
      reallocm (RSTCKMAT(dofe[i][j],dofe[i][j],lkm));
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
void quadlinaxisym::res_capacity_matrix (long eid,matrix &cm)
{
  long i,j,*rcn,*ccn;
  matrix lcm;

  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      rcn = new long [dofe[i][j]];
      ccn = new long [dofe[i][j]];
      reallocm (RSTCKMAT(dofe[i][j],dofe[i][j],lcm));
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
   
   TKr 13/09/2017, JK, 6.1.2002
*/
void quadlinaxisym::res_convection_vector (vector &f,long lcid,long eid,long leid)
{
  long *cn;
  vector lf;
  
  //  transi[lcid]==2 - element contains boundary with prescribed flux
  //  transi[lcid]==4 - element contains boundary with prescribed transmission and prescribed flux
  if ((Tt->elements[eid].transi[lcid]==2)||(Tt->elements[eid].transi[lcid]==4)){
    //  array for code numbers
    cn = new long [dofe[lcid][lcid]];
    //  array for nodal values of one matter
    reallocv (RSTCKVEC(dofe[lcid][lcid],lf));
    
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
   
   TKr 13/09/2017, JK, 6.1.2002
*/
void quadlinaxisym::res_transmission_vector (vector &f,long lcid,long eid,long leid)
{
  long *cn,i;
  vector lf;

  cn = new long [dofe[lcid][lcid]];
  //  code numbers
  codnum (cn,lcid);
  reallocv (RSTCKVEC(dofe[lcid][lcid],lf));
  
  //  transi[i]==3 - element contains boundary with prescribed transmission
  //  transi[i]==4 - element contains boundary with prescribed transmission and prescribed flux
  if ((Tt->elements[eid].transi[lcid]==3)||(Tt->elements[eid].transi[lcid]==4)){
    
    //coupled b.c.
    for (i=0;i<ntm;i++){    
      nullv (lf);
      
      transmission_vector (lf,lcid,eid,leid,i);
      
      locglob (f.a,lf.a,cn,dofe[lcid][lcid]);
    }
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
void quadlinaxisym::res_quantity_source_vector (vector &sv,vector &nodval,long lcid,long eid)
{
  long *cn;
  vector lsv;

  cn = new long [dofe[lcid][lcid]];
  reallocv (RSTCKVEC(dofe[lcid][lcid],lsv));
  quantity_source_vector (lsv,nodval,eid,lcid,lcid);
  codnum (cn,lcid);
  locglob (sv.a,lsv.a,cn,dofe[lcid][lcid]);
  delete [] cn;
}

/**
   function assembles resulting element internal fluxes vector

   @param eid - element id
   @param elemif - resulting internal fluxes %vector of one element
   
   JK, 6.1.2002
*/
void quadlinaxisym::res_internal_fluxes (long eid,vector &elemif)
{
  long i;
  vector lif, tdnv;
  ivector cn;
  matrix cm;

  for (i=0;i<ntm;i++){
    reallocv (RSTCKIVEC(dofe[i][i], cn));
    reallocv (RSTCKVEC(dofe[i][i], lif));
    internal_fluxes (i, eid, lif);
    codnum (cn.a, i);
    locglob (elemif.a, lif.a, cn.a, dofe[i][i]);
  }

  reallocv (RSTCKVEC(ndofe, tdnv));
  reallocv (RSTCKVEC(ndofe, lif));
  reallocm (RSTCKMAT(ndofe, ndofe, cm));
  nodalderivatives(eid, tdnv);
  res_capacity_matrix (eid, cm);
  mxv (cm, tdnv, lif);

  subv (elemif, lif, elemif);
}

/**
   function computes element quantity integral
   
   @param eid - element id
   @param nodval - %vector of quantity nodal values

   @retval f - element quantity integral

   TKr, 30.1.2004
*/
double quadlinaxisym::total_integral(long eid,vector &nodval)
{
  long i,j;
  double xinp,xi,eta,ww1,ww2,jac,value,f;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w(ASTCKVEC(2)),gp(ASTCKVEC(2)),t(ASTCKVEC(nne));

  Tt->give_elemnodes (eid,nodes);
  

  Tt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,2);

  f = 0.0;

  for (i=0;i<2;i++){
    xi=gp[i];  ww1=w[i];
    for (j=0;j<2;j++){
      eta=gp[j];  ww2=w[j];

      jac_2d (jac,x,y,xi,eta);
      xinp = approx (xi,eta,x);
      value = approx (xi,eta,nodval);

      jac*=xinp*ww1*ww2*value;
      f = f + jac;
    }
  }
  return(f);
}


/**
   function assembles resulting element boundary flux vector

   @param f - resulting boundary flux %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   
   TKr, 28.2.2004
*/
void quadlinaxisym::res_boundary_flux (vector &f,long lcid,long eid,long leid)
{
  long *cn,i;
  vector lf;

  cn = new long [dofe[lcid][lcid]];
  codnum (cn,lcid);
  
  //coupled b.c.
  for (i=0;i<ntm;i++){    
    reallocv (RSTCKVEC(dofe[lcid][lcid],lf));
        
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
   
   TKr, 16/07/2013
*/
void quadlinaxisym::volume_rhs_vector (long /*lcid*/,long eid,long ri,long ci,vector &vrhs)
{
  long i,j,ii;
  double xinp,xi,eta,ww1,ww2,jac;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w(ASTCKVEC(intordkm[ri][ci])),gp(ASTCKVEC(intordkm[ri][ci])),a(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(ncomp,dofe[ri][ci])),d;
  matrix km(ASTCKMAT(dofe[ri][ci],2));
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  node coordinates
  Tt->give_node_coord2d (x,y,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  
  nullm (km);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];
  
  for (i=0;i<intordkm[ri][ci];i++){
    xi=gp[i];  ww1=w[i];
    for (j=0;j<intordkm[ri][ci];j++){
      eta=gp[j];  ww2=w[j];
      //  matrix of gradients
      grad_matrix (gm,x,y,xi,eta,jac);
      
      //  matrix of conductivity of the material
      reallocm(RSTCKMAT(ncomp,ncomp,d));
      
      Tm->volume_rhs (d,ii,ri,ci,ncomp);
      
      xinp = approx (xi,eta,x);

      jac*=ww1*ww2*xinp;
      
      //  contribution to the volume_rhs integral of the element
      nnjac (km,gm,d,jac);
      
      ii++;
    }
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
   
   TKr, 16/05/2018
*/
void quadlinaxisym::volume_rhs_vector2 (long /*lcid*/,long eid,long ri,long ci,vector &vrhs)
{
  long i,j,k,ipp;
  double jac,xi,eta,w1,w2,xinp,c;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordcm[ri][ci]),gp(intordcm[ri][ci]),t(nne);
  matrix n(1,dofe[ri][ci]);
  vector nn(dofe[ri][ci]);
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  node coordinates
  Tt->give_node_coord2d (x,y,eid);
  //  thickness
  Tc->give_thickness (eid,nodes,t);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordcm[ri][ci]);
  
  nullv (vrhs);
  
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  for (i=0;i<intordcm[ri][ci];i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordcm[ri][ci];j++){
      eta=gp[j];  w2=w[j];
      jac_2d (jac,x,y,xi,eta);
      bf_matrix (n,xi,eta);
      
      c=Tm->volume_rhs2 (ipp,ri,ci);

      xinp = approx (xi,eta,x);

      jac*=w1*w2*xinp*c;
   
      for (k=0;k<n.n;k++){
	nn[k] = n[0][k];
      }

      addmultv(vrhs,nn,jac);
      ipp++;
    }
  }
}



/**
   function assembles resulting element volume right-hand side

   @param eid - element id
   @param lcid - load case id
   @param f - resulting volume right-hand side %vector of one element

   TKr, 16/07/2013
*/
void quadlinaxisym::res_volume_rhs_vector (vector &f,long eid,long /*lcid*/)
{
  long i,*cn;
  vector lf;

  for (i=0;i<ntm;i++){
    cn = new long [dofe[i][i]];
    reallocv (RSTCKVEC(dofe[i][i],lf));
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

   TKr, 16/05/2018
*/
void quadlinaxisym::res_volume_rhs_vector2 (vector &f,long eid,long /*lcid*/)
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
   function computes correct fluxes at integration points on element

   @param eid - element id
   
   TKo according TKr, 04.01.2019
*/
void quadlinaxisym::intpointflux (long eid)
{
  long i,j,ii,jj,k,ipp;
  
  for (k=0;k<Tp->ntm;k++){
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
	
	ipp=Tt->elements[eid].ipp[ii][jj];
	
	for (i=0;i<intordkm[ii][jj];i++){
	  for (j=0;j<intordkm[ii][jj];j++){
	    //  computation of correct fluxes
	    if (Tp->fluxcomp==1)
	      Tm->computenlfluxes (k,ipp);
	    
	    ipp++;
	  }
	}
	for (i=0;i<intordcm[ii][jj];i++){
	  for (j=0;j<intordcm[ii][jj];j++){
	    //  computation of correct fluxes
	    if (Tp->fluxcomp==1)
	      Tm->computenlfluxes (k,ipp);
	    
	    ipp++;
	  }
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
   
   TKo according TKr, 04/01/2019
*/
void quadlinaxisym::nod_grads_ip (long eid)
{
  long i,j,k,ipp;
  double area;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  vector grad(ASTCKVEC(ncomp));

  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  nodip_planelq (ipp,intordkm[0][0],ipnum);
  
  //  node numbers of the element
  Tt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    for (k=0;k<Tp->ntm;k++){
      //  gradients at the closest integration point
      Tm->givegrad (k,ipnum[i],grad);
      
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
   
   TKo according TKr, 04/01/2019
*/
void quadlinaxisym::nod_fluxes_ip (long eid)
{
  long i,j,k,ipp;
  double area;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  vector flux(ASTCKVEC(ncomp));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  nodip_planelq (ipp,intordkm[0][0],ipnum);
  
  //  node numbers of the element
  Tt->give_elemnodes (eid,nod);
  
  intpointflux (eid);

  for (i=0;i<nne;i++){
    for (k=0;k<Tp->ntm;k++){
      //  fluxes at the closest integration point
      Tm->givefluxes (k,ipnum[i],flux);
      
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
   
   Tomas Koudelka, 04.01.2019
*/
void quadlinaxisym::nod_eqother_ip (long eid)
{
  long i,j,k,ipp,ncompo;
  double area;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  vector eqother;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  nodip_planelq (ipp,intordkm[0][0],ipnum);

  
  //  node numbers of the element
  Tt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    ncompo = Tm->givencompeqother (ipnum[i],0);
    reallocv (RSTCKVEC(ncompo,eqother));
    Tm->giveeqother (ipnum[i],0,ncompo,eqother.a);
    
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
void quadlinaxisym::nod_others_comp (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp,ncompother,ii;
  vector h,r(ASTCKVEC(ndofe));
  ivector nod(ASTCKIVEC(nne));
  double area, other;
  
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
   
   TKr 13/09/2017, JK, 19.8.2004
*/
void quadlinaxisym::convection_vector (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w(ASTCKVEC(intordkm[ri][ci])),gp(ASTCKVEC(intordkm[ri][ci]));
  vector list(ASTCKVEC(nned*ned)),trc(ASTCKVEC(nned*ned)),nodval(ASTCKVEC(nne)),av(ASTCKVEC(nne)),coef(ASTCKVEC(nne));
  matrix km(ASTCKMAT(nne,nne));
  
  nullv (v);
  
  //  element nodes
  Tt->give_elemnodes (eid,nodes);
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
    
    if (bc[i]==2 || bc[i]==3 || bc[i]==5){
      //  nodal values on actual edge
      edgenodeval (i,nodval,list);
      //  nodal values of auxiliary coefficients
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
   
   TKr 13/09/2017, JK, 19.8.2004
*/
void quadlinaxisym::transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,j,ipp,leid;
  bocontypet *bc;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w(ASTCKVEC(intordkm[ri][ci])),gp(ASTCKVEC(intordkm[ri][ci]));
  vector trc(ASTCKVEC(nned*ned)),coeff(ASTCKVEC(nne));
  vector trr(ASTCKVEC(nned*ned));

  //  element nodes
  Tt->give_elemnodes (eid,nodes);
  //  node coordinates on element
  Tt->give_node_coord2d (x,y,eid);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);

  //  indicators of boundary conditions
  bc = new bocontypet [ned];

  //  correspondence between element ordering and boundary element ordering
  for (i=0;i<Tb->lc[lcid].neb;i++){
    if (Tb->lc[lcid].elemload[i].eid==eid){
      leid=i;
      
      Tb->lc[lcid].elemload[leid].give_bc (bc);
      //  transmission coefficients
      Tb->lc[lcid].elemload[leid].give_trc (Tp->time,Tb->lc[lcid].nodval,lcid,ci,trc);
      //  transmission/radiation coefficients
      Tb->lc[lcid].elemload[leid].give_trr (Tp->time,Tb->lc[lcid].nodval,trr);
      
      //  integration point number
      if (Tp->savemode==0)
	ipp=Tt->elements[eid].ipp[ri][ci];
      if (Tp->savemode==1)
	ipp=Tt->elements[eid].ipp[0][0];
      
      //  loop over edges
      for (j=0;j<ned;j++){
	
	if (bc[j]==5 || bc[j]>10){
	  // combined b. c.
	  if (bc[j]==100){
	    //  transformation of coefficients
	    transf_coeff (j,coeff,trr,eid,ri,ci,ipp,bc);
	  }
	  else{
	    //  transformation of coefficients
	    transf_coeff (j,coeff,trc,eid,ri,ci,ipp,bc);
	  }
	  //  matrix obtained from integration over edges
	  edge_integral (j,x,y,intordkm[ri][ci],gp,w,coeff,km);
	}
      }
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
   TKr 13/09/2017
*/
void quadlinaxisym::transmission_vector (vector &v,long lcid,long eid,long leid,long cid)
{
  long i,ipp;
  bocontypet *bc;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne)),w(ASTCKVEC(intordkm[lcid][cid])),gp(ASTCKVEC(intordkm[lcid][cid]));
  vector list(ASTCKVEC(nned*ned)),trc(ASTCKVEC(nned*ned)),trr(ASTCKVEC(nned*ned)),nodval(ASTCKVEC(nne)),av(ASTCKVEC(nne)),coef(ASTCKVEC(nne));
  matrix km(ASTCKMAT(nne,nne));
  
  nullv (v);
  
  //  element nodes
  Tt->give_elemnodes (eid,nodes);
  //  node coordinates on element
  Tt->give_node_coord2d (x,y,eid);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[lcid][cid]);
  

  //  indicators of boundary conditions
  bc = new bocontypet[ned];
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
    
    if (bc[i]==5 || bc[i]>10){
      //  transformation of nodal values
      transf_val (i,nodval,list,trc,trr,eid,lcid,cid,ipp,bc);
      //  determination of transmission coefficients on edge
      transf_coeff (i,coef,trc,eid,lcid,cid,ipp,bc);
  
      nullm (km);
      //  matrix obtained from integration over edges
      edge_integral (i,x,y,intordkm[lcid][cid],gp,w,coef,km);
      //  matrix-vector multiplication
      mxv (km,nodval,av);
      //  contribution added to element vector
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
   
   TKr 13/09/2017, TKr, 28.2.2004
*/
void quadlinaxisym::boundary_flux (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w(ASTCKVEC(intordkm[ri][ci])),gp(ASTCKVEC(intordkm[ri][ci]));
  vector list(ASTCKVEC(nned*ned)),trc(ASTCKVEC(nned*ned)),trr(ASTCKVEC(nned*ned)),nodval(ASTCKVEC(nne)),av(ASTCKVEC(nne)),coef(ASTCKVEC(nne));
  matrix km(ASTCKMAT(nne,nne));
  
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
    
    if (bc[i]==5 || bc[i]>10){
      //  transformation of nodal values
      transf_flux (i,nodval,list,trc,trr,eid,ri,ci,ipp,bc);
      //  determination of transmission coefficients on edge
      //transf_coeff (i,coef,trc,eid,ri,ci,ipp,bc);
      //new:
      edgenodeval (i,coef,trc);
 
      nullm (km);
      //  matrix obtained from integration over edges
      edge_integral (i,x,y,intordkm[ri][ci],gp,w,coef,km);
      //  matrix-vector multiplication
      mxv (km,nodval,av);
      //  contribution added to element vector
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
void quadlinaxisym::edge_integral (long edg,vector &x,vector &y,long intord,vector &gp,vector &w,
				 vector &coef,matrix &km)
{
  long i;
  double xi,eta,jac,ipval,xinp;
  vector av(ASTCKVEC(nne));
  matrix n(ASTCKMAT(1,nne));
  
  if (edg==0){
    eta=1.0;
    for (i=0;i<intord;i++){
      xi=gp[i];
      
      bf_matrix (n,xi,eta);
      jac1d_2d (jac,x,y,xi,edg);
      
      ipval=approx (xi,eta,coef);
      xinp=approx (xi,eta,x);
      
      jac*=w[i]*ipval*xinp;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

  if (edg==1){
    xi=-1.0;
    for (i=0;i<intord;i++){
      eta=gp[i];
      
      bf_matrix (n,xi,eta);
      jac1d_2d (jac,x,y,eta,edg);
      
      ipval=approx (xi,eta,coef);
      xinp=approx (xi,eta,x);
      
      jac*=w[i]*ipval*xinp;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

  if (edg==2){
    eta=-1.0;
    for (i=0;i<intord;i++){
      xi=gp[i];
      
      bf_matrix (n,xi,eta);
      jac1d_2d (jac,x,y,xi,edg);
      
      ipval=approx (xi,eta,coef);
      xinp=approx (xi,eta,x);
      
      jac*=w[i]*ipval*xinp;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

  if (edg==3){
    xi=1.0;
    for (i=0;i<intord;i++){
      eta=gp[i];
      
      bf_matrix (n,xi,eta);
      jac1d_2d (jac,x,y,eta,edg);
      
      ipval=approx (xi,eta,coef);
      xinp=approx (xi,eta,x);
      
      jac*=w[i]*ipval*xinp;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

}

void quadlinaxisym::transf_flux (long edg,vector &coeff,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double tr,new_flux;
  ivector nodes(ASTCKIVEC(nne)),edgenod(ASTCKIVEC(nned));

  //  zeroing of array of coefficients
  nullv (coeff);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required edge
  linquadrilat_edgnod (edgenod.a,edg);

  //  actual position in the array list
  j=edg*nned;

  for (i=0;i<nned;i++){
    tr = trc[j+i];
    if (bc[edg]==90){
      tr = trr[j+i]/trc[j+i];
    }
    
    //  node number
    k=edgenod[i];
    
    Tm->transmission_flux(new_flux,list[j+i],tr,ri,ci,nodes[k],bc[edg],ipp);
    
    coeff[k]=new_flux;
  }
  
}


void quadlinaxisym::transf_coeff (long edg,vector &coeff,vector &list,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double new_coeff;
  ivector nodes(ASTCKIVEC(nne)),edgenod(ASTCKIVEC(nned));

  //  zeroing of array of coefficients
  nullv (coeff);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required edge
  linquadrilat_edgnod (edgenod.a,edg);
  
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
void quadlinaxisym::transf_val (long edg,vector &nodval,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double tr,new_nodval;
  ivector nodes(ASTCKIVEC(nne)),edgenod(ASTCKIVEC(nned));
  
  //  zeroing of array of nodal values
  nullv (nodval);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required edge
  linquadrilat_edgnod (edgenod.a,edg);
  
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
void quadlinaxisym::edgenodeval (long edg,vector &nodval,vector &list)
{
  long i,j;
  ivector edgenod(ASTCKIVEC(nned));
  
  nullv (nodval);
  
  linquadrilat_edgnod (edgenod.a,edg);
  
  for (i=0;i<nned;i++){
    j=edgenod[i];
    nodval[j]=list[edg*nned+i];
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
   
   Created by 25/05/2016 TKr according to TKo
*/
void quadlinaxisym::transq_nodval (long eid,vector &nodval,nonmechquant nmq)
{
  long i,ipid;
  ivector ipnum(ASTCKIVEC(nne));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipid=Tt->elements[eid].ipp[0][0];
  nodip_planelq (ipid,intordkm[0][0], ipnum);

  for (i=0;i<nne;i++)
    {
      //copy transport (non-mechanical) quantity from closest int. point
      nodval[i] = Tm->givetransq(nmq, ipnum[i]);
    }
}



/**
   Function computes transport (non-mechanical) quantities at nodes of element.

   @param eid - element id
   @param nodval - %vector of nodal values of all required quantities, i.e., 
                   nodal value of i-th quantity in j-th node is given by nodval[i*ncnv+j] where ncnv 
                   is the number of calculated nodes on eid-th element.
   @param ncne - number of computed nodes on element (only first ncne of nodes is calculated)
   @param nq - number of required transport quantities
   @param qt - array of types of required non-mechanical, i.e. transport quantities
   
   @return The function does not return anything.
   
   Created by 25/05/2016 TKr according to TKo and TKr
*/
void quadlinaxisym::transq_nodval_comp (long eid,vector &nodval, long ncne, long nq, nonmechquant *qt)
{
  long i,j,ipid,ncompo;
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
  nodip_planelq (ipid,intordkm[0][0], ipnum);

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
      if (ipnum[i] == ipid) // values of the first int. point may be rewritten -> take values from bakup
	Tm->storeeqother(ipid, 0, ncompo, ipb.eqother);
      else
	Tm->storeeqother(ipid, 0, ncompo, Tm->ip[ipnum[i]].eqother);
    }
    ncompo = Tm->ip[ipnum[i]].ncompother;
    if (ncompo){
      // take eqother values for the given node from the closest integration point 
      if (ipnum[i] == ipid) // values of the first int. point may be rewritten -> take values from bakup
	Tm->storeother(ipid, 0, ncompo, ipb.other);
      else
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
void quadlinaxisym::transq_init_nodval (long eid,vector &nodval,nonmechquant nmq)
{
  long i,ipid;
  ivector ipnum(ASTCKIVEC(nne)), enod(ASTCKIVEC(nne));
  vector ipav;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipid=Tt->elements[eid].ipp[0][0];
  nodip_planelq (ipid,intordkm[0][0], ipnum);
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
   Function computes initial transport (non-mechanical) quantities at nodes of element.

   @param eid - element id
   @param nodval - %vector of intial nodal values of all required quantities, i.e., 
                   nodal value of i-th quantity in j-th node is given by nodval[i*ncnv+j] where ncnv 
                   is the number of calculated nodes on eid-th element.
   @param ncne - number of computed nodes on element (only first ncne of nodes is calculated)
   @param nq - number of required transport quantities
   @param qt - array of types of required non-mechanical, i.e. transport quantities
   
   @return The function does not return anything.
   
   Created by TKo 5.6.2018
*/
void quadlinaxisym::transq_init_nodval_comp (long eid,vector &nodval, long ncne, long nq, nonmechquant *qt)
{
  long i,j,ipid,ncompo;
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
  nodip_planelq (ipid,intordkm[0][0], ipnum);

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
      if (ipnum[i] == ipid) // values of the first int. point may be rewritten -> take values from bakup
	Tm->storeeqother(ipid, 0, ncompo, ipb.eqother);
      else
	Tm->storeeqother(ipid, 0, ncompo, Tm->ip[ipnum[i]].eqother);
    }
    ncompo = Tm->ip[ipnum[i]].ncompother;
    if (ncompo){
      // take eqother values for the given node from the closest integration point 
      if (ipnum[i] == ipid) // values of the first int. point may be rewritten -> take values from bakup
	Tm->storeother(ipid, 0, ncompo, ipb.other);
      else
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
long quadlinaxisym::ipcoord (long eid, long ipp, long ri, long ci, vector &coord)
{
  long i, j, ii;
  double xi, eta;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), w, gp;

  Tt->give_node_coord2d(x, y, eid);
  
  //  integration points for the conductivity matrix
  reallocv(RSTCKVEC(intordkm[ri][ci], gp));
  reallocv(RSTCKVEC(intordkm[ri][ci], w));
  gauss_points(gp.a, w.a, intordkm[ri][ci]);
	
  ii = Tt->elements[eid].ipp[ri][ci];
  for (i=0; i<intordkm[ri][ci]; i++)
  {
    for (j=0; j<intordkm[ri][ci]; j++)
    {    
      if (ii == ipp)
      {
        xi = gp[i];
        eta = gp[j];
        coord[0] = approx(xi, eta, x);
        coord[1] = approx(xi, eta, y);
        coord[2] = 0.0;
        return 0;
      }
      ii++;
    }
  }
	
  //  integration points for the capacity matrix
  reallocv(RSTCKVEC(intordcm[ri][ci], gp));
  reallocv(RSTCKVEC(intordcm[ri][ci], w));
  gauss_points (gp.a, w.a, intordcm[ri][ci]);
  
  for (i=0; i<intordcm[ri][ci]; i++)
  {
    for (j=0; j<intordcm[ri][ci]; j++)
    {    
      if (ii == ipp)
      {
        xi = gp[i];
        eta = gp[j];
        coord[0] = approx(xi, eta, x);
        coord[1] = approx(xi, eta, y);
        coord[2] = 0.0;
        return 0;
      }
      ii++;
    }
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
long quadlinaxisym::ipncoord (long eid, long ipp, vector &ncoord)
{
  long i, j, ii, ri, ci;
  double xi, eta;
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
        for (j=0; j<intordkm[ri][ci]; j++)
        {    
          if (ii == ipp)
          {
            xi = gp[i];
            eta = gp[j];
            ncoord[0] = xi;
            ncoord[1] = eta;
            ncoord[2] = 0.0;
            return 0;
          }
          ii++;
        }
      }
	
      //  integration points for the capacity matrix
      reallocv(RSTCKVEC(intordcm[ri][ci], gp));
      reallocv(RSTCKVEC(intordcm[ri][ci], w));
      gauss_points (gp.a, w.a, intordcm[ri][ci]);
  
      for (i=0; i<intordcm[ri][ci]; i++)
      {
        for (j=0; j<intordcm[ri][ci]; j++)
        {    
          if (ii == ipp)
          {
            xi = gp[i];
            eta = gp[j];
            ncoord[0] = xi;
            ncoord[1] = eta;
            ncoord[2] = 0.0;
            return 0;
          }
          ii++;
        }
      }
      if (Tp->savemode==1)  break;
    }
    if (Tp->savemode==1)  break;
  }

  return 1;
}



/**
   Function computes the total flux through an edge of an element
   
   @param lcid - load case id
   @param eid - element id (id in the list of all elements)
   @param beid - id in the list of selected elements
   @param flux - array of fluxes computed
   
   TKr 29/11/2019 according to JK
*/
void quadlinaxisym::surface_flux (long lcid,long eid,long beid,double *fluxes)
{

  long i,ipp,fid;
  double q,qn,area;
  bocontypet *bc;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  vector w, gp;
  vector flux(ASTCKVEC(ncomp)),fl(ASTCKVEC(ncomp)),n(ASTCKVEC(ncomp));
  matrix gm(ASTCKMAT(ncomp,nne));
  
  nullv (flux);
  
  Tt->give_node_coord2d (x,y,eid);

  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[lcid][lcid];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  //the first int. point
  Tm->computenlfluxes (lcid,ipp);    
  Tm->givefluxes (lcid,ipp,flux);

  //or volume average
  /* 
     reallocv (RSTCKVEC(intordkm[lcid][lcid],gp));
     reallocv (RSTCKVEC(intordkm[lcid][lcid],w));
     
     gauss_points (gp.a,w.a,intordkm[lcid][lcid]);
     
     for (i=0;i<intordkm[lcid][lcid];i++){
     xi=gp[i];
     for (j=0;j<intordkm[lcid][lcid];j++){
     eta=gp[j];
     xinp = approx (xi,eta,x);
     
     Tm->computenlfluxes (lcid,ipp);
     
     Tm->givefluxes (lcid,ipp,fl);
     
     grad_matrix (gm,x,y,xi,eta,jac);
     
     //mtxv (gm,fl,contr);//this was for internal_fluxes = \int B^T D B r d\Omega
     
     //cmulv (jac*w[i]*w[j]*xinp,fl);
     cmulv (jac*w[i]*w[j],fl);
     
     //  the total flux over the element volume
     addv (fl,flux,flux);
     
     ipp++;
     }
     }
     
     //average flux on the element
     area = element_area (eid);
     cmulv (1.0/area,flux); 
  */
  
  //  indicators of boundary conditions
  bc = new bocontypet [ned];
  Tb->bf[lcid].elemload[beid].give_bc (bc);
  
  //  loop over surfaces
  for (i=0;i<ned;i++){
    
    if (bc[i]==1){
      //  function constructs the outer normal vector
      quadlin_normal_vectors (i,x,y,n);
      //  area of the i-th element surface
      area = quadlinaxisym_surface_area (i,x,y);
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

  delete [] bc;
}
