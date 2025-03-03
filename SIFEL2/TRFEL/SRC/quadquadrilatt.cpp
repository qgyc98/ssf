/*
  File:     quadquadrilatt.cpp
  Author:   Jaroslav Kruis, 31.3.2001
  Purpose:  twodimensional quadrilateral element with quadratic approximation functions
*/

#include "globalt.h"
#include "quadquadrilatt.h"
#include "genfile.h"
#include "globmatt.h"

quadquadrilatt::quadquadrilatt (void)
{
  long i;

  //  number of nodes on element
  nne=8;
  //  number of edges
  ned=4;
  //  number of nodes on one edge
  nned=3;
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
    ordering[0][4]=5;  ordering[0][5]=6;  ordering[0][6]=7;  ordering[0][7]=8;
    dofe[0][0]=8;  intordkm[0][0]=3;  intordcm[0][0]=3;
    nip[0][0]=18;
    ndofe=8;  napfun=1;
    break;
  }
  case twomediacoup:{
    ordering[0][0]=1;   ordering[0][1]=3;   ordering[0][2]=5;   ordering[0][3]=7;
    ordering[0][4]=9;   ordering[0][5]=11;  ordering[0][6]=13;  ordering[0][7]=15;

    ordering[1][0]=2;   ordering[1][1]=4;   ordering[1][2]=6;   ordering[1][3]=8;
    ordering[1][4]=10;  ordering[1][5]=12;  ordering[1][6]=14;  ordering[1][7]=16;
    
    intordkm[0][0]=3;  intordkm[0][1]=3;  intordkm[1][0]=3;  intordkm[1][1]=3;
    intordcm[0][0]=3;  intordcm[0][1]=3;  intordcm[1][0]=3;  intordcm[1][1]=3;
    
    if (Tp->savemode==0){
      nip[0][0]=18;      nip[0][1]=18;      nip[1][0]=18;      nip[1][1]=18;
    }
    if (Tp->savemode==1){
      nip[0][0]=18;      nip[0][1]=0;      nip[1][0]=0;      nip[1][1]=0;
    }
    
    dofe[0][0]=8;      dofe[0][1]=8;      dofe[1][0]=8;      dofe[1][1]=8;

    ndofe=16;  napfun=2;
    break;
  }
  case threemediacoup:{
    ordering[0][0]=1;   ordering[0][1] =4;  ordering[0][2] =7;  ordering[0][3] =10;
    ordering[0][4]=13;  ordering[0][5] =16; ordering[0][6] =19; ordering[0][7] =22;

    ordering[1][0]=2;   ordering[1][1] =5;  ordering[1][2] =8;  ordering[1][3] =11;
    ordering[1][4]=14;  ordering[1][5] =17; ordering[1][6] =20; ordering[1][7] =23;

    ordering[2][0]=3;   ordering[2][1] =6;  ordering[2][2] =9;  ordering[2][3] =12;
    ordering[2][4]=15;  ordering[2][5] =18; ordering[2][6] =21; ordering[2][7] =24;

    intordkm[0][0]=3;  intordkm[0][1]=3;  intordkm[0][2]=3;
    intordkm[1][0]=3;  intordkm[1][1]=3;  intordkm[1][2]=3;
    intordkm[2][0]=3;  intordkm[2][1]=3;  intordkm[2][2]=3;

    intordcm[0][0]=3;  intordcm[0][1]=3;  intordcm[0][2]=3;
    intordcm[1][0]=3;  intordcm[1][1]=3;  intordcm[1][2]=3;
    intordcm[2][0]=3;  intordcm[2][1]=3;  intordcm[2][2]=3;
    
    if (Tp->savemode==0){
      nip[0][0]=18;  nip[0][1]=18;  nip[0][2]=18;
      nip[1][0]=18;  nip[1][1]=18;  nip[1][2]=18;
      nip[2][0]=18;  nip[2][1]=18;  nip[2][2]=18;
    }
    if (Tp->savemode==1){
      nip[0][0]=18;  nip[0][1]=0;   nip[0][2]=0;
      nip[1][0]=0;   nip[1][1]=0;   nip[1][2]=0;
      nip[2][0]=0;   nip[2][1]=0;   nip[2][2]=0;
    }
    
    dofe[0][0]=8;  dofe[0][1]=8;  dofe[0][2]=8;
    dofe[1][0]=8;  dofe[1][1]=8;  dofe[1][2]=8;
    dofe[2][0]=8;  dofe[2][1]=8;  dofe[2][2]=8;

    ndofe=24;  napfun=3;
    break;
  }
  default:{
    print_err("unknown type of transported matter is required",__FILE__,__LINE__,__func__);
  }
  }
}

quadquadrilatt::~quadquadrilatt (void)
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

void quadquadrilatt::codnum (long *cn,long ri)
{
  long i;
  for (i=0;i<nne;i++){
    cn[i]=ordering[ri][i];
  }
}

/**
   function approximates function defined by nodal values
   
   @param xi,eta - coordinates on element
   @param nodval - %vector of nodal values

   JK, 8.5.2002
*/
double quadquadrilatt::approx (double xi,double eta,vector &nodval)
{
  double f;
  vector bf(nne);
  
  bf_quad_4_2d (bf.a,xi,eta);
  
  scprd (bf,nodval,f);

  return f;
}

/**
   function assembles coordinates of integration points in block [ri][ci]
   
   @param eid - element id
   @param ri - row index
   @param ci - column index
   @param ipcoord - array containing coordinates of integration points
   
   8.5.2002
*/
void quadquadrilatt::ipcoordblock (long eid,long ri,long ci,double **coord)
{
  long i,j,k;
  double xi,eta;
  vector x(nne),y(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]);
  
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

   @param eid - element id

  JK, 8.5.2002 
*/
void quadquadrilatt::intpointval (long eid)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,val;
  vector r(ASTCKVEC(ndofe)), t(ASTCKVEC(nne)), gp, w;

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
	  for (j=0;j<intordkm[ii][jj];j++){
	    eta=gp[j];
	    
	    val = approx (xi,eta,t);
	    Tm->ip[ipp].av[k]=val;
	    ipp++;
	  }
	}
	
	//  integration points for the capacity matrix
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
   
   @param eid - element id
   
   JK, 8.5.2002
*/
void quadquadrilatt::intpointgrad (long eid)
{
  long i,j,ii,jj,k,ipp;
  double xi,eta,jac;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), r(ASTCKVEC(ndofe)), t(ASTCKVEC(nne));
  vector gp,w,grad(ASTCKVEC(ncomp));
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

	//  integration points for the conductivity matrix
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
	

	//  integration points for the capacity matrix
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
void quadquadrilatt::intpointother (long eid)
{
  long i, j, k, ii, jj, ipp, ncompo, nodid;
  double xi, eta, val;
  ivector nodes(nne);
  vector r, t(nne), gp, w;
  
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
	  for (j=0;j<intordkm[ii][jj];j++){
	    eta=gp[j];
	    
	    val = approx (xi,eta,t);
	    Tm->ip[ipp].other[k]=val;
	    ipp++;
	  }
	}
	
	
	//  integration points for the capacity matrix
	reallocv (intordcm[ii][jj],gp);
	reallocv (intordcm[ii][jj],w);
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
   function assembles %matrix of base functions

   @param n - %matrix of base functions
   @param xi,eta - natural coordinates

   JK, 25.9.2001
*/
void quadquadrilatt::bf_matrix (matrix &n,double xi,double eta)
{
  fillm (0.0,n);
  bf_quad_4_2d (n.a,xi,eta);
}

/**
   function returns one approximation function evaluated in required point

   @param f - value of approximation function
   @param xi,eta - natural coordinates
   @param i - number of approximation function

   JK, 1.2.2003
*/
void quadquadrilatt::give_approx_fun (double &f,double xi,double eta,long i)
{
  double n[4];
  bf_quad_4_2d (n,xi,eta);
  f=n[i];
}


/**
   function assembles gradient of %matrix of base functions

   @param gm - gradient %matrix
   @param x,y - array containing node coordinates
   @param xi,eta - natural coorodinates
   @param jac - Jacobian

   JK, 25.9.2001
*/
void quadquadrilatt::grad_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac)
{
  long i;
  vector dx(nne),dy(nne);

  dx_bf_quad_4_2d (dx.a,xi,eta);
  dy_bf_quad_4_2d (dy.a,xi,eta);

  derivatives_2d (dx,dy,jac,x,y,xi,eta);

  fillm (0.0,gm);

  for (i=0;i<nne;i++){
    gm[0][i]=dx[i];
    gm[1][i]=dy[i];
  }

}



/**
   function computes conductivity %matrix of 2D problems for one transported matter
   finite element with quadratic approximation functions

   @param lcid - number of load case
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param km - conductivity %matrix

   JK, 25.9.2001
*/
void quadquadrilatt::conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,j,ii;
  double thick,xi,eta,ww1,ww2,jac;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]),t(nne);
  matrix gm(ncomp,dofe[ri][ci]),d;

  matrix n(1,dofe[ri][ci]);
 
  Tt->give_elemnodes (eid,nodes);
  Tc->give_thickness (eid,nodes,t);

  Tt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordkm[ri][ci]);

  fillm (0.0,km);
  
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
      reallocm(ncomp,ncomp,d);
      Tm->matcond (d,ii,ri,ci);

      thick = approx (xi,eta,t);

      if (jac<0.0){
	print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	jac=fabs(jac);
      }

      jac*=thick*ww1*ww2;

      //  contribution to the conductivity matrix of the element
      bdbj (km.a,gm.a,d.a,jac,gm.m,gm.n);
      
      //convective terms
      reallocm(1,ncomp,d);
      Tm->matcond2(d,ii,ri,ci);
      bf_matrix (n, xi, eta);
      bdbjac(km, n, d, gm, jac);	
      
      ii++;  
    }
  }
  
  if ((Tt->elements[eid].transi[lcid]==3)||(Tt->elements[eid].transi[lcid]==4)){
    //  additional matrix due to transmission
    transmission_matrix (lcid,eid,ri,ci,km);
  }
}



/**
   function computes capacity %matrix of 2D problems for one transported matter
   finite element with quadratic approximation functions

   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param cm - capacity %matrix

   JK, 4.10.2001
*/
void quadquadrilatt::capacity_matrix (long eid,long ri,long ci,matrix &cm)
{
  long i,j,ii;
  double jac,xi,eta,w1,w2,thick,rho,c;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordcm[ri][ci]),gp(intordcm[ri][ci]),t(nne),dens(nne);
  matrix n(1,dofe[ri][ci]);

  Tt->give_elemnodes (eid,nodes);
  Tc->give_thickness (eid,nodes,t);
  Tc->give_density (eid,nodes,dens);

  Tt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordcm[ri][ci]);

  fillm (0.0,cm);
  
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

      c=Tm->capcoeff (ii,ri,ci);
      thick = approx (xi,eta,t);
      rho = approx (xi,eta,dens);

      if (jac<0.0){
	print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	jac=fabs(jac);
      }

      jac*=w1*w2*thick*rho*c;

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
void quadquadrilatt::quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci)
{
  long i,j;
  double jac,xi,eta,w1,w2,thick;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordcm[ri][ci]),gp(intordcm[ri][ci]),t(nne),dens(nne),v(dofe[ri][ci]);
  matrix n(1,dofe[ri][ci]),nm(dofe[ri][ci],dofe[ri][ci]);

  Tt->give_elemnodes (eid,nodes);

  Tc->give_thickness (eid,nodes,t);
  
  Tt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordcm[ri][ci]);
  
  fillm (0.0,nm);

  for (i=0;i<intordcm[ri][ci];i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordcm[ri][ci];j++){
      eta=gp[j];  w2=w[j];
      jac_2d (jac,x,y,xi,eta);
      bf_matrix (n,xi,eta);

      thick = approx (xi,eta,t);

      if (jac<0.0){
	print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	jac=fabs(jac);
      }

      jac*=w1*w2*thick;
      
      nnj (nm.a,n.a,jac,n.m,n.n);
    }
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
void quadquadrilatt::internal_fluxes (long lcid,long eid,vector &ifl)
{
  long i,j,ipp;
  double xi,eta,jac,thick;
  ivector nodes(nne);
  vector x(nne),y(nne),w,gp,t(nne),fl(ncomp),contr(dofe[lcid][lcid]);
  matrix gm(ncomp,dofe[lcid][lcid]);
  

  Tt->give_elemnodes (eid,nodes);
  Tc->give_thickness (eid,nodes,t);
  Tt->give_node_coord2d (x,y,eid);
  
  reallocv (intordkm[lcid][lcid],gp);
  reallocv (intordkm[lcid][lcid],w);
  
  gauss_points (gp.a,w.a,intordkm[lcid][lcid]);
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[lcid][lcid];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  for (i=0;i<intordkm[lcid][lcid];i++){
    xi=gp[i];
    for (j=0;j<intordkm[lcid][lcid];j++){
      eta=gp[j];
      thick = approx (xi,eta,t);
      
      Tm->computenlfluxes (lcid,ipp);
      
      Tm->givefluxes (lcid,ipp,fl);

      grad_matrix (gm,x,y,xi,eta,jac);
      if (jac<0.0){
	print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	jac=fabs(jac);
      }

      mtxv (gm,fl,contr);

      cmulv (jac*w[i]*w[j]*thick,contr);
      
      addv (contr,ifl,ifl);
      
      ipp++;
    }
  }
  
}

/**
   function assembles resulting element conductivity %matrix

   @param eid - element id
   @param lcid - load case id
   @param km - resulting conductivity %matrix of one element

   JK, 6.1.2002
*/
void quadquadrilatt::res_conductivity_matrix (long eid,long /*lcid*/,matrix &km)
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
void quadquadrilatt::res_capacity_matrix (long eid,matrix &cm)
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
   function assembles resulting element convection vector

   @param f - resulting convection %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   
   JK, 6.1.2002
*/
void quadquadrilatt::res_convection_vector (vector &f,long lcid,long eid,long leid)
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
void quadquadrilatt::res_transmission_vector (vector &f,long lcid,long eid,long leid)
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
   function assembles resulting element source vector

   @param sv - resulting source %vector of one element
   @param lcid - load case id
   @param eid - element id
   
   JK, 6.1.2002
*/
void quadquadrilatt::res_quantity_source_vector (vector &sv,vector &nodval,long lcid,long eid)
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
   function assembles resulting element internal fluxes vector

   @param eid - element id
   @param elemif - resulting internal fluxes %vector of one element
   
   JK, 6.1.2002
*/
void quadquadrilatt::res_internal_fluxes (long eid,vector &elemif)
{
  long i;
  vector lif, tdnv;
  ivector cn;
  matrix cm;

  for (i=0;i<ntm;i++){
    reallocv (RSTCKIVEC(dofe[i][i], cn));
    reallocv (RSTCKVEC(dofe[i][i],lif));
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
double quadquadrilatt::total_integral(long eid,vector &nodval)
{
  long i,j;
  double thick,xi,eta,ww1,ww2,jac,value,f;
  ivector nodes(nne);
  vector x(nne),y(nne),w(3),gp(3),t(nne);

  Tt->give_elemnodes (eid,nodes);
  Tc->give_thickness (eid,nodes,t);

  Tt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,3);

  f = 0.0;

  for (i=0;i<3;i++){
    xi=gp[i];  ww1=w[i];
    for (j=0;j<3;j++){
      eta=gp[j];  ww2=w[j];

      jac_2d (jac,x,y,xi,eta);
      thick = approx (xi,eta,t);
      value = approx (xi,eta,nodval);

      jac*=thick*ww1*ww2*value;
      f = f + jac;
    }
  }
  return(f);
}



/**
   function assembles resulting element boundary flux %vector

   @param f - resulting boundary flux %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - loaded element id
   
   TKr, 28.2.2004
*/
void quadquadrilatt::res_boundary_flux (vector &f,long lcid,long eid,long leid)
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
   function computes other values in nodes of element

   @param lcid - load case id
   @param eid - element id
   @param ri - row index
   @param ci - column index

   10.5.2002
*/
void quadquadrilatt::nod_others (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp,ncompother,ii;
  vector other,h,r(ASTCKVEC(ndofe));
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
    
    //Tt->nodes[nod[i]].storeother (ncompother,other);
  }
}

















/**
   function computes contribution to the convection vector

   \int_{Gamma_2} N^T N dGamma * nodal_flux_values
   
   @param v - array of nodal fluxes
   @param lcid - 
   @param eid - element id
   @param leid - id of loaded element
   @param ri,ci - row and column indices
   
   JK, 19.8.2004
*/
void quadquadrilatt::convection_vector (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]),t(nne);
  vector list(nned*ned),trc(nned*ned),nodval(nne),av(nne),coef(nne);
  matrix km(nne,nne);
  
  fillv (0.0,v);
  
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
    
    if (bc[i]==2 || bc[i]==3){
      //  nodal values on actual edge
      edgenodeval (i,nodval,list);
      
      edgenodeval (i,coef,trc);
  
      fillm (0.0,km);

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
void quadquadrilatt::transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km)
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
    
    if (bc[i]>10){
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
   
   @param v - transmission %vector of one matter
   @param lcid - load case id
   @param eid - element id
   @param leid - loaded element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)
   
   JK, 5.10.2001
   TKr, 30.1.2002 - new added
*/
void quadquadrilatt::transmission_vector (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]),t(nne);
  vector list(nned*ned),trc(nned*ned),trr(nned*ned),nodval(nne),av(nne),coef(nne);
  matrix km(nne,nne);
  
  fillv (0.0,v);
  
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
  
      fillm (0.0,km);
      //
      edge_integral (i,x,y,intordkm[ri][ci],gp,w,t,coef,km);
      
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
   @param leid - element id
   @param ri,ci - row and column indices
   
   TKr, 28.2.2004
*/
void quadquadrilatt::boundary_flux (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  ivector nodes(nne);
  vector x(nne),y(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]),t(nne);
  vector list(nned*ned),trc(nned*ned),trr(nned*ned),nodval(nne),av(nne),coef(nne);
  matrix km(nne,nne);
  
  fillv (0.0,v);
  
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
  
      fillm (0.0,km);
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
void quadquadrilatt::edge_integral (long edg,vector &x,vector &y,long intord,vector &gp,vector &w,
				    vector &t,vector &coef,matrix &km)
{
  long i;
  double xi,eta,jac,ipval,thick;
  matrix n(1,nne);
  
  if (edg==0){
    eta=1.0;
    for (i=0;i<intord;i++){
      xi=gp[i];
      
      bf_matrix (n,xi,eta);
      jac1d_2d (jac,x,y,xi,edg);
      
      ipval=approx (xi,eta,coef);
      thick=approx (xi,eta,t);
      
      jac*=w[i]*ipval*thick;
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
      thick=approx (xi,eta,t);
      
      jac*=w[i]*ipval*thick;
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
      thick=approx (xi,eta,t);
      
      jac*=w[i]*ipval*thick;
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
      thick=approx (xi,eta,t);
      
      jac*=w[i]*ipval*thick;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }

}

/**
   function
   
   @param edg - edge id
   @param coeff - array of 
   @param list - array of nodal values prescribed on boundary
   @param trc - array of transmission coefficients
   @param trr - array of transmission/radiation coefficients
   @param eid - element id
   @param ri,ci - row and column indices
   @param ipp - integration point id
   @param bc - array describing type of boundary condition
*/
void quadquadrilatt::transf_flux (long edg,vector &coeff,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double tr,new_flux;
  ivector nodes(nne),edgenod(nned);

  //  zeroing of array of coefficients
  fillv (0.0,coeff);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required edge
  quadquadrilat_edgnod (edgenod.a,edg);

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


void quadquadrilatt::transf_coeff (long edg,vector &coeff,vector &list,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double new_coeff;
  ivector nodes(nne),edgenod(nned);

  //  zeroing of array of coefficients
  fillv (0.0,coeff);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required edge
  quadquadrilat_edgnod (edgenod.a,edg);
  
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
void quadquadrilatt::transf_val (long edg,vector &nodval,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double tr,new_nodval;
  ivector nodes(nne),edgenod(nned);
  
  //  zeroing of array of nodal values
  fillv (0.0,nodval);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required edge
  quadquadrilat_edgnod (edgenod.a,edg);
  
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
void quadquadrilatt::edgenodeval (long edg,vector &nodval,vector &list)
{
  long i,j;
  ivector edgenod(nned);
  
  fillv (0.0,nodval);
  quadquadrilat_edgnod (edgenod.a,edg);
  
  for (i=0;i<nned;i++){
    j=edgenod[i];
    nodval[j]=list[edg*nned+i];
  }
}



/**
   function computes correct fluxes at integration points on element

   @param eid - element id
   
   TKr, 04/12/2010
*/
void quadquadrilatt::intpointflux (long eid)
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
   
   TKr, 04/12/2010
*/
void quadquadrilatt::nod_grads_ip (long eid)
{
  long i,j,k,ipp;
  ivector ipnum(nne),nod(nne);
  vector grad(ncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  nodip_planelq (ipp,intordkm[0][0],ipnum);
  
  //  node numbers of the element
  Tt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    for (k=0;k<Tp->ntm;k++){
      //  strains at the closest integration point
      Tm->givegrad (k,ipnum[i],grad);
      
      //  storage of strains to the node
      j=nod[i];
      Tt->nodes[j].storegrad (k,ncomp,grad);
    }
  }
  
}


/**
   function computes fluxes in nodes of element

   @param eid - element id
   
   TKr, 04/12/2010
*/
void quadquadrilatt::nod_fluxes_ip (long eid)
{
  long i,j,k,ipp;
  ivector ipnum(nne),nod(nne);
  vector flux(ncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  nodip_planelq (ipp,intordkm[0][0],ipnum);
  
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

   TKr, 04/12/2010
*/
void quadquadrilatt::nod_others_comp (long /*lcid*/,long eid,long ri,long ci)
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
long quadquadrilatt::ipcoord (long eid, long ipp, long ri, long ci, vector &coord)
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
long quadquadrilatt::ipncoord (long eid, long ipp, vector &ncoord)
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

