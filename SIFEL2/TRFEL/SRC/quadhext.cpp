#include "globalt.h"
#include "quadhext.h"
#include "genfile.h"
#include "globmatt.h"

quadhext::quadhext (void)
{
  long i;

  //  number of nodes on element
  nne=20;
  //  number of edges
  ned=12;
  //  number of nodes on one edge
  nned=3;
  //  number of surfaces
  nsurf=6;
  //  number of nodes on one surface
  nnsurf=8;
  //  geometric problem dimension (3D)
  ncomp=3;
  
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
    ordering[0][0]=1;    ordering[0][1]=2;    ordering[0][2]=3;    ordering[0][3]=4;
    ordering[0][4]=5;    ordering[0][5]=6;    ordering[0][6]=7;    ordering[0][7]=8;
    ordering[0][8]=9;    ordering[0][9]=10;   ordering[0][10]=11;  ordering[0][11]=12;
    ordering[0][12]=13;  ordering[0][13]=14;  ordering[0][14]=15;  ordering[0][15]=16;
    ordering[0][16]=17;  ordering[0][17]=18;  ordering[0][18]=19;  ordering[0][19]=20;
    
    dofe[0][0]=20;  intordkm[0][0]=3;  intordcm[0][0]=3;  nip[0][0]=54;
    ndofe=20;  napfun=1;
    break;
  }
    
  case twomediacoup:{
    ordering[0][0]=1;    ordering[0][1]=3;    ordering[0][2]=5;    ordering[0][3]=7;
    ordering[0][4]=9;    ordering[0][5]=11;   ordering[0][6]=13;   ordering[0][7]=15;
    ordering[0][8]=17;   ordering[0][9]=19;   ordering[0][10]=21;  ordering[0][11]=23;
    ordering[0][12]=25;  ordering[0][13]=27;  ordering[0][14]=29;  ordering[0][15]=31;
    ordering[0][16]=33;  ordering[0][17]=35;  ordering[0][18]=37;  ordering[0][19]=39;

    ordering[1][0]=2;    ordering[1][1]=4;    ordering[1][2]=6;    ordering[1][3]=8;
    ordering[1][4]=10;   ordering[1][5]=12;   ordering[1][6]=14;   ordering[1][7]=16;
    ordering[1][8]=18;   ordering[1][9]=20;   ordering[1][10]=22;  ordering[1][11]=24;
    ordering[1][12]=26;  ordering[1][13]=28;  ordering[1][14]=30;  ordering[1][15]=32;
    ordering[1][16]=34;  ordering[1][17]=36;  ordering[1][18]=38;  ordering[1][19]=40;


    intordkm[0][0]=3;  intordkm[0][1]=3;  intordkm[1][0]=3;  intordkm[1][1]=3;
    intordcm[0][0]=3;  intordcm[0][1]=3;  intordcm[1][0]=3;  intordcm[1][1]=3;
    
    if (Tp->savemode==0){
      nip[0][0]=54;      nip[0][1]=54;      nip[1][0]=54;      nip[1][1]=54;
    }
    if (Tp->savemode==1){
      nip[0][0]=54;      nip[0][1]=0;      nip[1][0]=0;      nip[1][1]=0;
    }
    
    dofe[0][0]=20;     dofe[0][1]=20;     dofe[1][0]=20;     dofe[1][1]=20;
    ndofe=40;  napfun=2;
    break;
  }

  case threemediacoup:{
    ordering[0][0]=1;    ordering[0][1]=4;    ordering[0][2]=7;    ordering[0][3]=10;
    ordering[0][4]=13;   ordering[0][5]=16;   ordering[0][6]=19;   ordering[0][7]=22;
    ordering[0][8]=25;   ordering[0][9]=28;   ordering[0][10]=31;  ordering[0][11]=34;
    ordering[0][12]=37;  ordering[0][13]=40;  ordering[0][14]=43;  ordering[0][15]=46;
    ordering[0][16]=49;  ordering[0][17]=52;  ordering[0][18]=55;  ordering[0][19]=58;

    ordering[1][0]=2;    ordering[1][1]=5;    ordering[1][2]=8;    ordering[1][3]=11;
    ordering[1][4]=14;   ordering[1][5]=17;   ordering[1][6]=20;   ordering[1][7]=23;
    ordering[1][8]=26;   ordering[1][9]=29;   ordering[1][10]=32;  ordering[1][11]=35;
    ordering[1][12]=38;  ordering[1][13]=41;  ordering[1][14]=44;  ordering[1][15]=47;
    ordering[1][16]=50;  ordering[1][17]=53;  ordering[1][18]=56;  ordering[1][19]=59;

    ordering[2][0]=3;    ordering[2][1]=6;    ordering[2][2]=9;    ordering[2][3]=12;
    ordering[2][4]=15;   ordering[2][5]=18;   ordering[2][6]=21;   ordering[2][7]=24;
    ordering[2][8]=27;   ordering[2][9]=30;   ordering[2][10]=33;  ordering[2][11]=36;
    ordering[2][12]=39;  ordering[2][13]=42;  ordering[2][14]=45;  ordering[2][15]=48;
    ordering[2][16]=51;  ordering[2][17]=54;  ordering[2][18]=57;  ordering[2][19]=60;

    intordkm[0][0]=3;  intordkm[0][1]=3;  intordkm[0][2]=3;
    intordkm[1][0]=3;  intordkm[1][1]=3;  intordkm[1][2]=3;
    intordkm[2][0]=3;  intordkm[2][1]=3;  intordkm[2][2]=3;

    intordcm[0][0]=3;  intordcm[0][1]=3;  intordcm[0][2]=3;
    intordcm[1][0]=3;  intordcm[1][1]=3;  intordcm[1][2]=3;
    intordcm[2][0]=3;  intordcm[2][1]=3;  intordcm[2][2]=3;
    
    if (Tp->savemode==0){
      nip[0][0]=54;      nip[0][1]=54;      nip[0][2]=54;
      nip[1][0]=54;      nip[1][1]=54;      nip[1][2]=54;
      nip[2][0]=54;      nip[2][1]=54;      nip[2][2]=54;
    }
    if (Tp->savemode==1){
      nip[0][0]=54;      nip[0][1]=0;       nip[0][2]=0; 
      nip[1][0]=0;       nip[1][1]=0;       nip[1][2]=0; 
      nip[2][0]=0;       nip[2][1]=0;       nip[2][2]=0; 
    }

    dofe[0][0]=20;     dofe[0][1]=20;     dofe[0][2]=20;
    dofe[1][0]=20;     dofe[1][1]=20;     dofe[1][2]=20;
    dofe[2][0]=20;     dofe[2][1]=20;     dofe[2][2]=20;

    ndofe=60;  napfun=3;
    break;
  }
  default:{
    print_err("unknown type of transported matter is required",__FILE__,__LINE__,__func__);
  }
  }
}


quadhext::~quadhext (void)
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
   function assembles code numbers
   the code numbers are used only on the element level
   
   @param cn - array containing code numbers
   @param ri - id of the transported media
   
   16.3.2004, JK
*/
void quadhext::codnum (long *cn,long ri)
{
  long i;
  for (i=0;i<nne;i++){
    cn[i]=ordering[ri][i];
  }
}

/**
   function approximates function defined by nodal values
   
   @param xi,eta,zeta - coordinates on element
   @param nodval - nodal values
   
   16.3.2004, JK
*/
double quadhext::approx (double xi,double eta,double zeta, vector &nodval)
{
  double f;
  vector bf(nne);
  bf_quad_hex_3d (bf.a,xi,eta,zeta);
  
  scprd (bf,nodval,f);

  return f;
}

/**
   function computes values at integration points from nodal values
   
   @param eid - element id
   
   16.3.2004, JK
*/
void quadhext::intpointval (long eid)
{
  long i,j,k,l,ii,jj,ipp;
  double xi,eta,zeta,val;
  vector r(ASTCKVEC(ndofe)),t(ASTCKVEC(nne)),gp,w;
  
  elemvalues(eid, r);
  
  for (l=0;l<Tp->ntm;l++){
    
    //  nodal values of a single variable
    for (i=0;i<dofe[l][l];i++){
      t[i]=r[ordering[l][i]-1];
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
	    for(k=0;k<intordkm[ii][jj];k++){
	      zeta=gp[k];
	      val = approx (xi,eta,zeta,t);
	      Tm->ip[ipp].av[l]=val;
	      ipp++;
	    }
	  }
	}
	
	reallocv (RSTCKVEC(intordcm[ii][jj],gp));
	reallocv (RSTCKVEC(intordcm[ii][jj],w));
	gauss_points (gp.a,w.a,intordcm[ii][jj]);
	
	for (i=0;i<intordcm[ii][jj];i++){
	  xi=gp[i];
	  for (j=0;j<intordcm[ii][jj];j++){
	    eta=gp[j];
	    for(k=0;k<intordcm[ii][jj];k++){
	      zeta=gp[k];
	      val = approx (xi,eta,zeta,t);
	      Tm->ip[ipp].av[l]=val;
	      ipp++;
	    }
	  }
	}
	
	if (Tp->savemode==1)  break;
      }
      if (Tp->savemode==1)  break;
    }
  }
  
}

/**
   function computes values at integration points from nodal values
   
   @param eid - element id
   
   16.3.2004, JK
*/
void quadhext::intpointgrad (long eid)
{
  long i,j,k,l,ii,jj,ipp;
  double xi,eta,zeta,jac;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne)), r(ASTCKVEC(ndofe));
  vector gp, w, grad(ASTCKVEC(ncomp)), t(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(ncomp,nne));
  
  Tt->give_node_coord3d (x,y,z,eid);
  elemvalues(eid, r);
  
  for (l=0;l<Tp->ntm;l++){
    
    //  nodal values of a single variable
    for (i=0;i<dofe[l][l];i++){
      t[i]=r[ordering[l][i]-1];
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
	    for(k=0;k<intordkm[ii][jj];k++){
	      zeta=gp[k];
	      
	      grad_matrix (gm,x,y,z,xi,eta,zeta,jac);
	      mxv (gm,t,grad);
	      Tm->storegrad (l,ipp,grad);
	      ipp++;
	    }
	  }
	}
	
	reallocv (RSTCKVEC(intordcm[ii][jj],gp));
	reallocv (RSTCKVEC(intordcm[ii][jj],w));
	gauss_points (gp.a,w.a,intordcm[ii][jj]);
	
	for (i=0;i<intordcm[ii][jj];i++){
	  xi=gp[i];
	  for (j=0;j<intordcm[ii][jj];j++){
	    eta=gp[j];
	    for(k=0;k<intordcm[ii][jj];k++){
	      zeta=gp[k];
	      
	      grad_matrix (gm,x,y,z,xi,eta,zeta,jac);
	      mxv (gm,t,grad);
	      Tm->storegrad (l,ipp,grad);
	      ipp++;
	    }
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
void quadhext::intpointother (long eid)
{
  long i, j, k, l, ii, jj, ipp, ncompo, nodid;
  double xi, eta, zeta, val;
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
	  for (j=0;j<intordkm[ii][jj];j++){
	    eta=gp[j];
	    for (l=0;l<intordkm[ii][jj];l++){
	      zeta=gp[l];
	      
	      val = approx (xi,eta,zeta,t);
	      Tm->ip[ipp].other[k]=val;
	      ipp++;
	    }
	  }
	}
	
	reallocv (intordcm[ii][jj],gp);
	reallocv (intordcm[ii][jj],w);
	gauss_points (gp.a,w.a,intordcm[ii][jj]);
	
	for (i=0;i<intordcm[ii][jj];i++){
	  xi=gp[i];
	  for (j=0;j<intordcm[ii][jj];j++){
	    eta=gp[j];
	    for (l=0;l<intordkm[ii][jj];l++){
	      zeta=gp[l];
	      
	      val = approx (xi,eta,zeta,t);
	      Tm->ip[ipp].other[k]=val;
	      ipp++;
	    }
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
   @param xi,eta,zeta - natural coordinates
   
   JK, 16.3.2004
*/
void quadhext::bf_matrix (matrix &n,double xi,double eta,double zeta)
{
  nullm (n);
  bf_quad_hex_3d (n.a,xi,eta,zeta);
}

/**
   function assembles %matrix of gradients of base functions
   
   @param gm - gradient %matrix
   @param x,y,z - array containing node coordinates
   @param xi,eta,zeta - natural coorodinates
   @param jac - Jacobian
   
   JK, 16.3.2004
*/
void quadhext::grad_matrix (matrix &gm,vector &x,vector &y,vector &z,double xi,double eta,double zeta,double &jac)
{
  long i;
  vector dx(nne),dy(nne),dz(nne);

  dx_bf_quad_hex_3d (dx.a,xi,eta,zeta);
  dy_bf_quad_hex_3d (dy.a,xi,eta,zeta);
  dz_bf_quad_hex_3d (dz.a,xi,eta,zeta);
  
  derivatives_3d (dx,dy,dz,jac,x,y,z,xi,eta,zeta);
  
  nullm (gm);
  
  for (i=0;i<nne;i++){
    gm[0][i]=dx[i];
    gm[1][i]=dy[i];
    gm[2][i]=dz[i];
  }
  
}

/**
   function computes conductivity %matrix of one transported medium
   
   @param lcid - load case id
   @param eid - element id
   @param ri, ci - row and column index
   @param km - conductivity %matrix
   
   JK, 16.3.2004
*/
void quadhext::conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,j,k,ii;
  double xi,eta,zeta,ww1,ww2,ww3,jac;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]),t(nne);
  matrix gm(ncomp,dofe[ri][ci]),d;
  matrix n(1,dofe[ri][ci]);

  Tt->give_elemnodes (eid,nodes);
    
  Tt->give_node_coord3d (x,y,z,eid);
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
      for (k=0;k<intordkm[ri][ci];k++){
	zeta=gp[k]; ww3=w[k];
	//  matrix of gradients
	grad_matrix (gm,x,y,z,xi,eta,zeta,jac);
	
	//  matrix of conductivity of the material
	reallocm(ncomp,ncomp,d);
	Tm->matcond (d,ii,ri,ci);
	
	jac*=ww1*ww2*ww3;
	
	//  contribution to the conductivity matrix of the element
	bdbj (km.a,gm.a,d.a,jac,gm.m,gm.n);
	
	//convective terms
	reallocm(1,ncomp,d);
	
	Tm->matcond2(d,ii,ri,ci);
	bf_matrix (n, xi, eta,zeta);
	bdbjac(km, n, d, gm, jac);	
	
	ii++;
      }
    }
  }
  
  if ((Tt->elements[eid].transi[lcid]==3)||(Tt->elements[eid].transi[lcid]==4)){
    transmission_matrix (lcid,eid,ri,ci,km);
  }
  
}


/**
   function computes capacity %matrix of one transported medium
   
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param cm - capacity %matrix
   
   JK, 16.3.2004
*/
void quadhext::capacity_matrix (long eid,long ri,long ci,matrix &cm)
{
  long i,j,k,ii;
  double jac,xi,eta,zeta,w1,w2,w3,rho,c;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),w(intordcm[ri][ci]),gp(intordcm[ri][ci]),t(nne),dens(nne);
  matrix n(1,dofe[ri][ci]);
  
  Tt->give_elemnodes (eid,nodes);
  Tc->give_density (eid,nodes,dens);

  Tt->give_node_coord3d (x,y,z,eid);
  gauss_points (gp.a,w.a,intordcm[ri][ci]);
  
  nullm (cm);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci]+intordkm[ri][ci]*intordkm[ri][ci]*intordkm[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0]+intordkm[0][0]*intordkm[0][0]*intordkm[0][0];
  
  
  for (i=0;i<intordcm[ri][ci];i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordcm[ri][ci];j++){
      eta=gp[j];  w2=w[j];
      for (k=0;k<intordcm[ri][ci];k++){
	zeta=gp[k]; w3=w[k];
	jac_3d (jac,x,y,z,xi,eta,zeta);
	bf_matrix (n,xi,eta,zeta);
	
	c=Tm->capcoeff (ii,ri,ci);
	rho = approx (xi,eta,zeta,dens);
	jac*=w1*w2*w3*rho*c;
	
	nnj (cm.a,n.a,jac,n.m,n.n);
	ii++;
      }
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
   
   JK, 16.3.2004
*/
void quadhext::quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci)
{
  long i,j,k;
  double jac,xi,eta,zeta,w1,w2,w3;
  vector x(nne),y(nne),z(nne),w(intordcm[ri][ci]),gp(intordcm[ri][ci]),dens(nne),v(dofe[ri][ci]);
  matrix n(1,dofe[ri][ci]),nm(dofe[ri][ci],dofe[ri][ci]);
  
  Tt->give_node_coord3d (x,y,z,eid);
  gauss_points (gp.a,w.a,intordcm[ri][ci]);
  
  nullm (nm);
  
  for (i=0;i<intordcm[ri][ci];i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordcm[ri][ci];j++){
      eta=gp[j];  w2=w[j];
      for (k=0;k<intordcm[ri][ci];k++){
	zeta=gp[k];  w3=w[k];

	jac_3d (jac,x,y,z,xi,eta,zeta);
	bf_matrix (n,xi,eta,zeta);
	
	jac*=w1*w2*w3;
	
	nnj (nm.a,n.a,jac,n.m,n.n);
      }
    }
  }
  mxv (nm,nodval,v);  addv (sv,v,sv);
  
}






/**
   function computes internal fluxes
   they are used in nonlinear analysis where they are compared with the prescribed fluxes
   
   @param lcid - load case id
   @param eid - elemet id
   @param ifl - %vector containing fluxes
   
   JK, 16.3.2004
*/
void quadhext::internal_fluxes (long lcid,long eid,vector &ifl)
{
  long i,j,k,ipp;
  double xi,eta,zeta,jac;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),w,gp,fl(ncomp),contr(dofe[lcid][lcid]);
  matrix gm(ncomp,dofe[lcid][lcid]);
  
  
  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord3d (x,y,z,eid);
    
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
      for (k=0;k<intordkm[lcid][lcid];k++){
	zeta=gp[k];
	
	Tm->computenlfluxes (lcid,ipp);
	
	Tm->givefluxes (lcid,ipp,fl);
	
	grad_matrix (gm,x,y,z,xi,eta,zeta,jac);
	
	mtxv (gm,fl,contr);
	
	cmulv (jac*w[i]*w[j]*w[k],contr);
	
	addv (contr,ifl,ifl);
	
	ipp++;
      }
    }
  }
  
}

/**
   function assembles resulting element conductivity %matrix

   @param eid - element id
   @param lcid - load case id
   @param km - resulting conductivity %matrix of one element
   
   JK, 16.3.2004
*/
void quadhext::res_conductivity_matrix (long eid,long /*lcid*/,matrix &km)
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
   function assembles resultant element capacity %matrix

   @param eid - element id
   @param cm - resulting capacity %matrix of one element
   
   JK, 16.3.2004
*/
void quadhext::res_capacity_matrix (long eid,matrix &cm)
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
   function computes resulting element convection %vector
   
   @param f - resulting convection %vector of one element
   @param lcid - load case id
   @param eid - element id
   
   JK, 16.3.2004
*/
void quadhext::res_convection_vector (vector &f,long lcid,long eid,long leid)
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
   function computes resulting element transmission %vector
   
   @param f - resulting transmission %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   
   JK, 16.3.2004
*/
void quadhext::res_transmission_vector (vector &f,long lcid,long eid,long leid)
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
   function computes resulting element source %vector
   
   @param sv - resulting source %vector of one element
   @param lcid - load case id
   @param eid - element id
   
   JK, 16.3.2004
*/
void quadhext::res_quantity_source_vector (vector &sv,vector &nodval,long lcid,long eid)
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
   function computes resulting element internal fluxes %vector
   
   @param eid - element id
   @param elemif - resulting %vector of internal fluxes on element
   
   JK, 16.3.2004
*/
void quadhext::res_internal_fluxes (long eid,vector &elemif)
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
double quadhext::total_integral(long eid,vector &nodval)
{
  long i,j,k;
  double jac,xi,eta,zeta,w1,w2,w3,value,f;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),w(2),gp(2);
  
  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord3d (x,y,z,eid);
  gauss_points (gp.a,w.a,2);

  f = 0.0;

  for (i=0;i<2;i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<2;j++){
      eta=gp[j];  w2=w[j];
      for (k=0;k<2;k++){
	zeta=gp[k]; w3=w[k];

	jac_3d (jac,x,y,z,xi,eta,zeta);
	value = approx (xi,eta,zeta,nodval);

	jac*=w1*w2*w3*value;
	f = f + jac;
      }
    }
  }
  return(f);
}





/**
   function assembles resulting element boundary flux %vector

   @param f - resulting boundary flux %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   
   TKr, 15.4.2004
*/
void quadhext::res_boundary_flux (vector &f,long lcid,long eid,long leid)
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
void quadhext::nod_others (long /*lcid*/,long eid,long ri,long ci)
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
    
  }

}











/**
   function computes contribution to the flux vector

   \int_{Gamma_2} N^T N dGamma * nodal_heat_flux_values
   
   @param v - array of nodal fluxes
   @param lcid - 
   @param eid - element id
   @param leid - id of loaded element
   @param ri,ci - row and column indices
   
   JK, 19.8.2004
*/
void quadhext::convection_vector (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  vector x(nne),y(nne),z(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]);
  vector list(nnsurf*nsurf),trc(nnsurf*nsurf),nodval(nne),av(nne),coef(nne);
  matrix km(nne,nne);
  
  nullv (v);
  
  //  node coordinates on element
  Tt->give_node_coord3d (x,y,z,eid);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  
  //  indicators of boundary conditions
  bc = new bocontypet [nsurf];
  Tb->lc[lcid].elemload[leid].give_bc (bc);
  //  nodal values of exterior variable
  Tb->lc[lcid].elemload[leid].give_nodval (lcid,list);
  
  //  auxiliary coefficients, necessary for function surface_integral
  for (i=0;i<nnsurf*nsurf;i++){
    trc[i]=1.0;
  }

  //  integration point number
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  //  loop over surfaces
  for (i=0;i<nsurf;i++){
    
    if (bc[i]==2 || bc[i]==3){
      //  nodal values on actual surface
      surfnodeval (i,nodval,list);
      
      surfnodeval (i,coef,trc);

      nullm (km);

      //  matrix obtained from integration over surface
      surface_integral (i,x,y,z,intordkm[ri][ci],gp,w,coef,km);
      
      mxv (km,nodval,av);
      
      addv (v,av,v);
    }
  }
  
  delete [] bc;
}

/**
   function computes transmission complement to the conductivity %matrix for one matter
   
   \int_{Gamma_3} N^T c_{tr} N dGamma
   
   @param lcid -
   @param eid - element id
   @param ri,ci - row and column indices
   @param km - 
   
   JK, 19.8.2004
*/
void quadhext::transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,ipp,leid;
  bocontypet *bc;
  vector x(nne),y(nne),z(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]);
  vector trc(nnsurf*nsurf),coeff(nne);
  
  //  node coordinates on element
  Tt->give_node_coord3d (x,y,z,eid);
  
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
  bc = new bocontypet [nsurf];
  Tb->lc[lcid].elemload[leid].give_bc (bc);
  //  transmission coefficients
  Tb->lc[lcid].elemload[leid].give_trc (Tp->time,Tb->lc[lcid].nodval,lcid,ci,trc);
  
  //  integration point number
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];

  //  loop over surfaces
  for (i=0;i<nsurf;i++){
    
    if (bc[i]>10){
      //  transformation of coefficients
      transf_coeff (i,coeff,trc,eid,ri,ci,ipp,bc);
      
      //  matrix obtained from integration over surface
      surface_integral (i,x,y,z,intordkm[ri][ci],gp,w,coeff,km);
    }
  }
  
  delete [] bc;
}

/**
   function computes contributions to the transmission %vector

   \int_{Gamma_3} N^T c_{tr} N dGamma * nodal_external_temperature
   
   @param v - transmission %vector
   @param lcid - load case id
   @param eid - element id
   @param leid - loaded element id
   @param ri,ci - row and column indices
   
   JK+TKr, 15.4.2004
*/
void quadhext::transmission_vector (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  vector x(nne),y(nne),z(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]);
  vector list(nnsurf*nsurf),trc(nnsurf*nsurf),trr(nnsurf*nsurf),nodval(nne),av(nne),coef(nne);
  matrix km(nne,nne);
  
  nullv (v);
  
  //  node coordinates on element
  Tt->give_node_coord3d (x,y,z,eid);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  
  //  indicators of boundary conditions
  bc = new bocontypet [nsurf];
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
  
  for (i=0;i<nsurf;i++){
    
    if (bc[i]>10){
      //  transformation of nodal values
      transf_val (i,nodval,list,trc,trr,eid,ri,ci,ipp,bc);
      //  determination of transmission coefficients on edge
      transf_coeff (i,coef,trc,eid,ri,ci,ipp,bc);

      nullm (km);
      //  matrix obtained from integration over surfaces
      surface_integral (i,x,y,z,intordkm[ri][ci],gp,w,coef,km);
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
   
   \int_{Gamma_3} N^T c_{tr} N dGamma * nodal_external_temperature
   
   @param v - %vector of nodal values
   @param lcid - load case id
   @param eid - element id
   @param leid - id of loaded element
   @param ri,ci - row and column indices
   
   JK
*/
void quadhext::boundary_flux (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  vector x(nne),y(nne),z(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]);
  vector list(nnsurf*nsurf),trc(nnsurf*nsurf),trr(nnsurf*nsurf),nodval(nne),av(nne),coef(nne);
  matrix km(nne,nne);
  
  nullv (v);
  
  //  node coordinates on element
  Tt->give_node_coord3d (x,y,z,eid);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  
  //  indicators of boundary conditions
  bc = new bocontypet [nsurf];
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
  
  for (i=0;i<nsurf;i++){
    
    if (bc[i]>10){
      //  transformation of nodal values
      transf_flux (i,nodval,list,trc,trr,eid,ri,ci,ipp,bc);
      //  determination of transmission coefficients on edge
      transf_coeff (i,coef,trc,eid,ri,ci,ipp,bc);

      nullm (km);
      //  matrix obtained from integration over surfaces
      surface_integral (i,x,y,z,intordkm[ri][ci],gp,w,coef,km);
      //  matrix-vector multiplication
      mxv (km,nodval,av);
      //  contribution added to element vector
      addv (v,av,v);
    }
  }
  
  delete [] bc;
}


/**
   function integrates N^T c N over surface
   
   @param surf - surface id (number of surface)
   @param x,y,z - coordinates of element nodes
   @param intord - order of numerical integration
   @param gp, w  - coordinates and weights of integration points
   @param coef - array of nodal values of coefficient
   @param km - output %matrix
   
   JK
*/
void quadhext::surface_integral (long surf,vector &x,vector &y,vector &z,long intord,vector &gp,vector &w,
				 vector &coef,matrix &km)
{
  long i,j;
  double xi,eta,zeta,jac,ipval;
  matrix n(1,nne);
  
  if (surf==0){
    xi=1.0;
    for (i=0;i<intord;i++){
      eta=gp[i];
      for (j=0;j<intord;j++){
	zeta=gp[j];
	
	bf_matrix (n,xi,eta,zeta);
	jac2d_3d (jac,x,y,z,eta,zeta,0);
	
	ipval=approx (xi,eta,zeta,coef);
	
	jac*=w[i]*w[j]*ipval;
	nnj (km.a,n.a,jac,n.m,n.n);
      }
    }
  }

  if (surf==1){
    eta=1.0;
    for (i=0;i<intord;i++){
      xi=gp[i];
      for (j=0;j<intord;j++){
	zeta=gp[j];
	
	bf_matrix (n,xi,eta,zeta);
	jac2d_3d (jac,x,y,z,xi,zeta,1);
	
	ipval=approx (xi,eta,zeta,coef);
	
	jac*=w[i]*w[j]*ipval;
	nnj (km.a,n.a,jac,n.m,n.n);
      }
    }
  }

  if (surf==2){
    xi=-1.0;
    for (i=0;i<intord;i++){
      eta=gp[i];
      for (j=0;j<intord;j++){
	zeta=gp[j];
	
	bf_matrix (n,xi,eta,zeta);
	jac2d_3d (jac,x,y,z,eta,zeta,2);
	
	ipval=approx (xi,eta,zeta,coef);
	
	jac*=w[i]*w[j]*ipval;
	nnj (km.a,n.a,jac,n.m,n.n);
      }
    }
  }

  if (surf==3){
    eta=-1.0;
    for (i=0;i<intord;i++){
      xi=gp[i];
      for (j=0;j<intord;j++){
	zeta=gp[j];
	
	bf_matrix (n,xi,eta,zeta);
	jac2d_3d (jac,x,y,z,xi,zeta,3);
	
	ipval=approx (xi,eta,zeta,coef);
	
	jac*=w[i]*w[j]*ipval;
	nnj (km.a,n.a,jac,n.m,n.n);
      }
    }
  }

  if (surf==4){
    zeta=1.0;
    for (i=0;i<intord;i++){
      xi=gp[i];
      for (j=0;j<intord;j++){
	eta=gp[j];
	
	bf_matrix (n,xi,eta,zeta);
	jac2d_3d (jac,x,y,z,xi,eta,4);
	
	ipval=approx (xi,eta,zeta,coef);
	
	jac*=w[i]*w[j]*ipval;
	nnj (km.a,n.a,jac,n.m,n.n);
      }
    }
  }

  if (surf==5){
    zeta=-1.0;
    for (i=0;i<intord;i++){
      xi=gp[i];
      for (j=0;j<intord;j++){
	eta=gp[j];
	
	bf_matrix (n,xi,eta,zeta);
	jac2d_3d (jac,x,y,z,xi,eta,5);
	
	ipval=approx (xi,eta,zeta,coef);
	
	jac*=w[i]*w[j]*ipval;
	nnj (km.a,n.a,jac,n.m,n.n);
      }
    }
  }

}

void quadhext::transf_flux (long surf,vector &coeff,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double tr,new_flux;
  ivector nodes(nne),surfnod(nnsurf);

  //  zeroing of array of coefficients
  nullv (coeff);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required surface
  quadhexahedral_surfnod (surfnod.a,surf);

  //  actual position in the array list
  j=surf*nnsurf;

  for (i=0;i<nnsurf;i++){
    if (bc[surf]==90){
      tr = trr[j+i]/trc[j+i];
    }

    //  node number
    k=surfnod[i];
    
    Tm->transmission_flux(new_flux,list[j+i],tr,ri,ci,nodes[k],bc[surf],ipp);
    
    coeff[k]=new_flux;
  }
  
}


void quadhext::transf_coeff (long surf,vector &coeff,vector &list,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double new_coeff;
  ivector nodes(nne),surfnod(nnsurf);

  //  zeroing of array of coefficients
  nullv (coeff);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required surface
  quadhexahedral_surfnod (surfnod.a,surf);

  //  actual position in the array list
  j=surf*nnsurf;

  for (i=0;i<nnsurf;i++){
    //  node number
    k=surfnod[i];
    
    Tm->transmission_transcoeff(new_coeff,list[j+i],ri,ci,nodes[k],bc[surf],ipp);
    
    coeff[k]=new_coeff;
  }
  
}


/**
   @param surf - number of required surface
   @param nodval - array of transformed nodal values
   @param list - array of nodal values defined on all surfaces
   @param trc -
   @param trr - 
   @param ri,ci - row and column indices
   @param ipp - integration point number
   @param bc - array defining boundary conditions
   
   JK, 19.8.2004
*/
void quadhext::transf_val (long surf,vector &nodval,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double tr,new_nodval;
  ivector nodes(nne),surfnod(nnsurf);
  
  //  zeroing of array of nodal values
  nullv (nodval);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required surface
  quadhexahedral_surfnod (surfnod.a,surf);
  
  //  actual position in the array list
  j=surf*nnsurf;
  
  //  loop over number of nodes on one surface
  for (i=0;i<nnsurf;i++){
    if (bc[surf]==90){
      tr = trr[j+i]/trc[j+i];
    }
    
    //  node number
    k=surfnod[i];
    
    Tm->transmission_nodval(new_nodval,list[j+i],tr,ri,ci,nodes[k],bc[surf],ipp);
    
    //  transformed nodal values
    nodval[k]=new_nodval;
  }
  
}

/**
   function picks up nodal values on required surface
   
   @param surf - number of required surface
   @param nodval - array of nodal values
   @param list - array of nodal values defined on all surfaces
   
   JK, 19.8.2004
*/
void quadhext::surfnodeval (long surf,vector &nodval,vector &list)
{
  long i,j;
  ivector surfnod(nnsurf);
  
  nullv (nodval);
  quadhexahedral_surfnod (surfnod.a,surf);
  
  for (i=0;i<nnsurf;i++){
    j=surfnod[i];
    nodval[j]=list[surf*nnsurf+i];
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
long quadhext::ipcoord (long eid, long ipp, long ri, long ci, vector &coord)
{
  long i, j, k, ii;
  double xi, eta, zeta;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne)), w, gp;

  Tt->give_node_coord3d(x, y, z, eid);
  
  //  integration points for the conductivity matrix
  reallocv(RSTCKVEC(intordkm[ri][ci], gp));
  reallocv(RSTCKVEC(intordkm[ri][ci], w));
  gauss_points(gp.a, w.a, intordkm[ri][ci]);
	
  ii = Tt->elements[eid].ipp[ri][ci];
  for (i=0; i<intordkm[ri][ci]; i++)
  {
    for (j=0; j<intordkm[ri][ci]; j++)
    {    
      for (k=0; k<intordkm[ri][ci]; k++)
      {    
        if (ii == ipp)
        {
          xi = gp[i];
          eta = gp[j];
          zeta = gp[k];
          coord[0] = approx(xi, eta, zeta, x);
          coord[1] = approx(xi, eta, zeta, y);
          coord[2] = approx(xi, eta, zeta, z);
          return 0;
        }
        ii++;
      }
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
      for (k=0; k<intordcm[ri][ci]; k++)
      {    
        if (ii == ipp)
        {
          xi = gp[i];
          eta = gp[j];
          zeta = gp[k];
          coord[0] = approx(xi, eta, zeta, x);
          coord[1] = approx(xi, eta, zeta, y);
          coord[2] = approx(xi, eta, zeta, z);
          return 0;
        }
        ii++;
      }
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
long quadhext::ipncoord (long eid, long ipp, vector &ncoord)
{
  long i, j, k, ii, ri, ci;
  double xi, eta, zeta;
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
          for (k=0; k<intordkm[ri][ci]; k++)
          {    
            if (ii == ipp)
            {
              xi = gp[i];
              eta = gp[j];
              zeta = gp[k];
              ncoord[0] = xi;
              ncoord[1] = eta;
              ncoord[2] = zeta;
              return 0;
            }
            ii++;
          }
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
          for (k=0; k<intordcm[ri][ci]; k++)
          {    
            if (ii == ipp)
            {
              xi = gp[i];
              eta = gp[j];
              zeta = gp[k];
              ncoord[0] = xi;
              ncoord[1] = eta;
              ncoord[2] = zeta;
              return 0;
            }
            ii++;
          }
        }
      }
      if (Tp->savemode==1)  break;
    }
    if (Tp->savemode==1)  break;
  }

  return 1;
}
