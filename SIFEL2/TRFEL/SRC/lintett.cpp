/*
  File:     lintett.cpp
  Author:   Jaroslav Kruis, 31.3.2003
  Purpose:  tetrahedral element with linear approximation functions
*/

#include "globalt.h"
#include "lintett.h"
#include "genfile.h"
#include "globmatt.h"
#include <math.h>

lintett::lintett (void)
{
  long i;

  //  number of nodes on element
  nne=4;
  //  number of edges
  ned=6;
  //  number of nodes on one edge
  nned=2;
  //  number of surfaces
  nsurf=4;
  //  number of nodes on one surface
  nnsurf=3;
  //  geometric problem dimension (3D)
  ncomp=3;

  intordb = 3;
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
    
    dofe[0][0]=4;  intordkm[0][0]=1;  intordcm[0][0]=4;  nip[0][0]=5;
    ndofe=4;  napfun=1;
    break;
  }

  case twomediacoup:{
    ordering[0][0]=1;  ordering[0][1]=3;  ordering[0][2]=5;  ordering[0][3]=7;
    ordering[1][0]=2;  ordering[1][1]=4;  ordering[1][2]=6;  ordering[1][3]=8;

    intordkm[0][0]=1;  intordkm[0][1]=1;  intordkm[1][0]=1;  intordkm[1][1]=1;
    intordcm[0][0]=4;  intordcm[0][1]=4;  intordcm[1][0]=4;  intordcm[1][1]=4;

    if (Tp->savemode==0){
      nip[0][0]=5;       nip[0][1]=5;       nip[1][0]=5;       nip[1][1]=5;
    }
    if (Tp->savemode==1){
      nip[0][0]=5;       nip[0][1]=0;       nip[1][0]=0;       nip[1][1]=0;
    }


    dofe[0][0]=4;  dofe[0][1]=4;  dofe[1][0]=4;  dofe[1][1]=4;
    ndofe=8;  napfun=2;
    break;
  }

  case threemediacoup:{
    ordering[0][0]=1;  ordering[0][1]=4;  ordering[0][2]=7;  ordering[0][3]=10;
    ordering[1][0]=2;  ordering[1][1]=5;  ordering[1][2]=8;  ordering[1][3]=11;
    ordering[2][0]=3;  ordering[2][1]=6;  ordering[2][2]=9;  ordering[2][3]=12;

    intordkm[0][0]=1;  intordkm[0][1]=1;  intordkm[0][2]=1;
    intordkm[1][0]=1;  intordkm[1][1]=1;  intordkm[1][2]=1;
    intordkm[2][0]=1;  intordkm[2][1]=1;  intordkm[2][2]=1;

    intordcm[0][0]=4;  intordcm[0][1]=4;  intordcm[0][2]=4;
    intordcm[1][0]=4;  intordcm[1][1]=4;  intordcm[1][2]=4;
    intordcm[2][0]=4;  intordcm[2][1]=4;  intordcm[2][2]=4;

    if (Tp->savemode==0){
      nip[0][0]=5;       nip[0][1]=5;       nip[0][2]=5;
      nip[1][0]=5;       nip[1][1]=5;       nip[1][2]=5;
      nip[2][0]=5;       nip[2][1]=5;       nip[2][2]=5;
    }

    if (Tp->savemode==1){
      nip[0][0]=5;       nip[0][1]=0;       nip[0][2]=0;
      nip[1][0]=0;       nip[1][1]=0;       nip[1][2]=0;
      nip[2][0]=0;       nip[2][1]=0;       nip[2][2]=0;
    }

    dofe[0][0]=4;      dofe[0][1]=4;      dofe[0][2]=4;
    dofe[1][0]=4;      dofe[1][1]=4;      dofe[1][2]=4;
    dofe[2][0]=4;      dofe[2][1]=4;      dofe[2][2]=4;

    ndofe=12;  napfun=3;
    break;
  }

  default:{
    print_err("unknown type of transported matter is required",__FILE__,__LINE__,__func__);
  }
  }
}


lintett::~lintett (void)
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

void lintett::codnum (ivector &cn,long ri)
{
  long i;
  for (i=0;i<nne;i++){
    cn[i]=ordering[ri][i];
  }
}

/**
   function computes volume of element
   
   @param eid - element id
   
   JK, 24. 3. 2017
*/
double lintett::element_volume (long eid)
{
  double jac,vol;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  Tt->give_node_coord3d (x,y,z,eid);
  
  //jac_3d (jac,x,y,z,0.25,0.25,0.25);
  jac = det3d(x.a, y.a, z.a);
  
  if(jac < 0.0){
    //jac = fabs (jac);
    print_err("wrong numbering of nodes on 3D element number %ld, negative volume! det = %e", __FILE__, __LINE__, __func__, eid+1, jac);
    abort();
  }
  
  vol = jac/6.0;
  
  return vol;
}

/**
   function approximates function defined by nodal values
   
   @param volcoord - volume coordinates
   @param nodval - %vector nodal values
   
   JK, 24.3.2002
*/
double lintett::approx (vector &volcoord,vector &nodval)
{
  double f;

  scprd (volcoord,nodval,f);
  
  return f;
}



/**
   The function computes approximated function value with the help of nodal values.
   
   @param xi, eta, zeta - natural coordinates
   @param nodval - nodal values
   
   20.8.2001
*/
double lintett::approx_nat(double xi, double eta, double zeta, vector &nodval)
{
  double f;
  vector volcoord(ASTCKVEC(4));
  
  volcoord[0]=xi;
  volcoord[1]=eta;
  volcoord[2]=zeta;
  volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
  
  scprd (volcoord,nodval,f);
  
  return f;
}



/**
   function computes values in integration points from nodal values
   
   @param eid - element id
   
   JK, 24.3.2002
*/
void lintett::intpointval (long eid)
{
  long i,k,ii,jj,ipp;
  double val;
  vector volcoord(ASTCKVEC(4)),r(ASTCKVEC(ndofe)),t(ASTCKVEC(nne));
  vector gp1,gp2,gp3,w;
  
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
        reallocv (RSTCKVEC(intordkm[ii][jj],gp3));
        reallocv (RSTCKVEC(intordkm[ii][jj],w));
        gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordkm[ii][jj]);
        
        ipp=Tt->elements[eid].ipp[ii][jj];
        for (i=0;i<intordkm[ii][jj];i++){
          volcoord[0]=gp1[i];
          volcoord[1]=gp2[i];
          volcoord[2]=gp3[i];
          volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
          
          val = approx (volcoord,t);
          Tm->ip[ipp].av[k]=val;
          ipp++;
        }
	
        reallocv (RSTCKVEC(intordcm[ii][jj],gp1));
        reallocv (RSTCKVEC(intordcm[ii][jj],gp2));
        reallocv (RSTCKVEC(intordcm[ii][jj],gp3));
        reallocv (RSTCKVEC(intordcm[ii][jj],w));
        gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordcm[ii][jj]);
        for (i=0;i<intordcm[ii][jj];i++){
          volcoord[0]=gp1[i];
          volcoord[1]=gp2[i];
          volcoord[2]=gp3[i];
          volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
          
          val = approx (volcoord,t);
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
   
   @param eid - element id
   
   TKr, 05/04/2011
*/
void lintett::intpointval_puc (long eid)
{
  long i, k, ii, jj, ipp;
  double val;
  vector volcoord(ASTCKVEC(4)), r(ASTCKVEC(ndofe)), t(ASTCKVEC(nne));
  vector gp1, gp2, gp3, w;
  
  elemvalues_puc(eid, r);
  
  for (k=0;k<Tp->ntm;k++){
    
    //  nodal values of a single variable
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
        
        reallocv (RSTCKVEC(intordkm[ii][jj], gp1));
        reallocv (RSTCKVEC(intordkm[ii][jj], gp2));
        reallocv (RSTCKVEC(intordkm[ii][jj], gp3));
        reallocv (RSTCKVEC(intordkm[ii][jj], w));
        gauss_points_tet (gp1.a, gp2.a, gp3.a, w.a, intordkm[ii][jj]);
        
        ipp=Tt->elements[eid].ipp[ii][jj];
        for (i=0;i<intordkm[ii][jj];i++){
          volcoord[0]=gp1[i];
          volcoord[1]=gp2[i];
          volcoord[2]=gp3[i];
          volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
          
          val = approx (volcoord,t);
          Tm->ip[ipp].av[k]=val;
          ipp++;
        }
	
        reallocv (RSTCKVEC(intordcm[ii][jj],gp1));
        reallocv (RSTCKVEC(intordcm[ii][jj],gp2));
        reallocv (RSTCKVEC(intordcm[ii][jj],gp3));
        reallocv (RSTCKVEC(intordcm[ii][jj],w));
        gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordcm[ii][jj]);
        for (i=0;i<intordcm[ii][jj];i++){
          volcoord[0]=gp1[i];
          volcoord[1]=gp2[i];
          volcoord[2]=gp3[i];
          volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
          
          val = approx (volcoord,t);
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
  The function computes initial values in integration points from initial nodal values.
  (function interpolates the initial nodal values to integration points)

  @param eid - element id
 
  @return The function does not return anything but stores computed values at int. points of the given element.

  Created by TKo, 4.7.2018
*/
void lintett::initintpointval (long eid)
{
  long i,k,ii,jj,ipp,ndofn,cndofn;
  double val;
  ivector enod(ASTCKIVEC(nne));
  vector volcoord(ASTCKVEC(4)), r(ASTCKVEC(ndofe)), t(ASTCKVEC(nne));
  vector rr, gp1, gp2, gp3, w;
  
  //  initial nodal values on element
  Tt->give_elemnodes (eid,enod);

  for(i=0, cndofn=0; i<nne; i++)
  {
    ndofn = Tt->give_ndofn (enod[i]);
    // make reference rr to the vector of nodal values on element
    makerefv (rr, r.a+cndofn, ndofn);
    // get initial nodal values and store them in the vector r with the help of reference rr
    initnodval2 (enod[i], rr);
    cndofn += ndofn;
  }
  
  for (k=0;k<Tp->ntm;k++){
    
    //  nodal values of a single variable
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
	
        reallocv (RSTCKVEC(intordkm[ii][jj], gp1));
        reallocv (RSTCKVEC(intordkm[ii][jj], gp2));
        reallocv (RSTCKVEC(intordkm[ii][jj], gp3));
        reallocv (RSTCKVEC(intordkm[ii][jj], w));
        gauss_points_tet (gp1.a, gp2.a, gp3.a, w.a, intordkm[ii][jj]);
        
        ipp=Tt->elements[eid].ipp[ii][jj];
        for (i=0;i<intordkm[ii][jj];i++){
          volcoord[0]=gp1[i];
          volcoord[1]=gp2[i];
          volcoord[2]=gp3[i];
          volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
          
          val = approx (volcoord,t);
          Tm->ip[ipp].av[k]=val;
          ipp++;
        }
	
        reallocv (RSTCKVEC(intordcm[ii][jj], gp1));
        reallocv (RSTCKVEC(intordcm[ii][jj], gp2));
        reallocv (RSTCKVEC(intordcm[ii][jj], gp3));
        reallocv (RSTCKVEC(intordcm[ii][jj], w));
        gauss_points_tet (gp1.a, gp2.a, gp3.a, w.a, intordcm[ii][jj]);
        for (i=0;i<intordcm[ii][jj];i++){
          volcoord[0]=gp1[i];
          volcoord[1]=gp2[i];
          volcoord[2]=gp3[i];
          volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
          
          val = approx (volcoord, t);
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
   
   JK, 24.3.2002
*/
void lintett::intpointgrad (long eid)
{
  long i,k,ii,jj,ipp;
  double det;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne)),volcoord(ASTCKVEC(4));
  vector bb(ASTCKVEC(4)),cc(ASTCKVEC(4)),dd(ASTCKVEC(4)),r(ASTCKVEC(ndofe)),t(ASTCKVEC(nne));
  vector gp1,gp2,gp3,w;
  vector grad(ASTCKVEC(ncomp));
  matrix gm(ASTCKMAT(ncomp,nne));

  Tt->give_node_coord3d (x,y,z,eid);
  elemvalues(eid, r);
  
  det = det3d (x.a,y.a,z.a);
  volb_3d (bb.a,y.a,z.a,det);
  volc_3d (cc.a,x.a,z.a,det);
  vold_3d (dd.a,x.a,y.a,det);
  
  //  matrix of gradients
  grad_matrix (gm,bb,cc,dd);

  for (k=0;k<Tp->ntm;k++){
    
    //  nodal values of a single variable
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
        
        reallocv (RSTCKVEC(intordkm[ii][jj],gp1));
        reallocv (RSTCKVEC(intordkm[ii][jj],gp2));
        reallocv (RSTCKVEC(intordkm[ii][jj],gp3));
        reallocv (RSTCKVEC(intordkm[ii][jj],w));
        gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordkm[ii][jj]);
        
        ipp=Tt->elements[eid].ipp[ii][jj];
        for (i=0;i<intordkm[ii][jj];i++){
          volcoord[0]=gp1[i];
          volcoord[1]=gp2[i];
          volcoord[2]=gp3[i];
          volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
          
          mxv (gm,t,grad);
          //Tm->storegrad (k,ipp,0,grad);
          Tm->storegrad (k,ipp,grad);
          ipp++;
        }
	
        reallocv (RSTCKVEC(intordcm[ii][jj],gp1));
        reallocv (RSTCKVEC(intordcm[ii][jj],gp2));
        reallocv (RSTCKVEC(intordcm[ii][jj],gp3));
        reallocv (RSTCKVEC(intordcm[ii][jj],w));
        gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordcm[ii][jj]);
        for (i=0;i<intordcm[ii][jj];i++){
          volcoord[0]=gp1[i];
          volcoord[1]=gp2[i];
          volcoord[2]=gp3[i];
          volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
          
          mxv (gm,t,grad);
          //Tm->storegrad (k,ipp,0,grad);
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
void lintett::intpointother (long eid)
{
  long i, k, ii, jj, ipp, ncompo, nodid;
  double val;
  ivector nodes(ASTCKIVEC(nne));
  vector volcoord(ASTCKVEC(4)), t(ASTCKVEC(nne));
  vector gp1, gp2, gp3, w, r;
  
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
        
        reallocv (RSTCKVEC(intordkm[ii][jj],gp1));
        reallocv (RSTCKVEC(intordkm[ii][jj],gp2));
        reallocv (RSTCKVEC(intordkm[ii][jj],gp3));
        reallocv (RSTCKVEC(intordkm[ii][jj],w));
        gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordkm[ii][jj]);
        
        ipp=Tt->elements[eid].ipp[ii][jj];
        for (i=0;i<intordkm[ii][jj];i++){
          volcoord[0]=gp1[i];
          volcoord[1]=gp2[i];
          volcoord[2]=gp3[i];
          volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
          
          val = approx (volcoord,t);
          Tm->ip[ipp].other[k]=val;
          ipp++;
        }
	
        reallocv (RSTCKVEC(intordcm[ii][jj],gp1));
        reallocv (RSTCKVEC(intordcm[ii][jj],gp2));
        reallocv (RSTCKVEC(intordcm[ii][jj],gp3));
        reallocv (RSTCKVEC(intordcm[ii][jj],w));
        gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordcm[ii][jj]);
        for (i=0;i<intordcm[ii][jj];i++){
          volcoord[0]=gp1[i];
          volcoord[1]=gp2[i];
          volcoord[2]=gp3[i];
          volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
          
          val = approx (volcoord,t);
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
   function assembles %matrix of base functions
   
   @param n - %matrix of base functions
   @param volcoord - %vector of volume coordinates
   
   JK, 24.3.2002
*/
void lintett::bf_matrix (matrix &n,vector &volcoord)
{
  nullm (n);
  
  n[0][0]=volcoord[0];
  n[0][1]=volcoord[1];
  n[0][2]=volcoord[2];
  n[0][3]=volcoord[3];
}

/**
   function assembles gradient of %matrix of base functions
   
   @param gm - gradient %matrix
   @param b,c,d - array containing coefficients of volume coordinates
   
   JK, 24.3.2002
*/
void lintett::grad_matrix (matrix &gm,vector &b,vector &c,vector &d)
{
  long i;
  nullm (gm);
  
  for (i=0;i<nne;i++){
    gm[0][i]=b[i];
    gm[1][i]=c[i];
    gm[2][i]=d[i];
  }
  
}

/**
   function computes conductivity %matrix of one transported medium
   
   @param lcid - load case id
   @param eid - element id
   @param ri, ci - row and column index
   @param km - conductivity %matrix
   
   JK, 24.3.2002
*/
void lintett::conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long ipp;
  double jac,det;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne)),bb(ASTCKVEC(4));
  vector cc(ASTCKVEC(4)),dd(ASTCKVEC(4)),volcoord(ASTCKVEC(4));
  vector gp1,gp2,gp3,w;
  matrix gm(ASTCKMAT(ncomp,dofe[ri][ci])),d(ASTCKMAT(ncomp,ncomp));
  
  matrix n(ASTCKMAT(1,dofe[ri][ci]));

  Tt->give_node_coord3d (x,y,z,eid);
  
  nullm (km);

  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];

  det = det3d (x.a,y.a,z.a);
      
  if(det < 0.0){
    //det = fabs(det);
    print_err("wrong numbering of nodes on 3D element number %ld, negative volume! det = %e",__FILE__,__LINE__,__func__,eid+1,det);
    abort();
  }

  volb_3d (bb.a,y.a,z.a,det);
  volc_3d (cc.a,x.a,z.a,det);
  vold_3d (dd.a,x.a,y.a,det);
  
  reallocv (RSTCKVEC(intordkm[ri][ci],gp1));
  reallocv (RSTCKVEC(intordkm[ri][ci],gp2));
  reallocv (RSTCKVEC(intordkm[ri][ci],gp3));
  reallocv (RSTCKVEC(intordkm[ri][ci],w));

  gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordkm[ri][ci]);

  volcoord[0]=gp1[0];
  volcoord[1]=gp2[0];
  volcoord[2]=gp3[0];

  volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];

  //  matrix of gradients
  grad_matrix (gm,bb,cc,dd);
  
  //  matrix of conductivity of the material
  Tm->matcond (d,ipp,ri,ci);
  
  jac=fabs(det)/6.0;
  
  //  contribution to the conductivity matrix of the element
  bdbj (km.a,gm.a,d.a,jac,gm.m,gm.n);
  
  //convective terms
  if (Tp->tmatt == threemediacoup)
  {
    reallocm(RSTCKMAT(1,ncomp,d));
      
    Tm->matcond2(d,ipp,ri,ci);
    bf_matrix (n,volcoord);
    bdbjac(km, n, d, gm, jac);	
  }
  
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

   TKr, 05/04/2011
*/
void lintett::l_matrix (long /*lcid*/,long eid,long ri,long ci,matrix &lm)
{
  long ipp;
  double jac,det;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne)),bb(ASTCKVEC(4)),cc(ASTCKVEC(4));
  vector dd(ASTCKVEC(4)),volcoord(ASTCKVEC(4));
  vector gp1,gp2,gp3,w;
  matrix gm(ASTCKMAT(ncomp,dofe[ri][ci])),d(ASTCKMAT(ncomp,ncomp));
  
  matrix n(ASTCKMAT(1,dofe[ri][ci]));

  Tt->give_node_coord3d (x,y,z,eid);
  
  nullm (lm);

  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];

  det = det3d (x.a,y.a,z.a);
  volb_3d (bb.a,y.a,z.a,det);
  volc_3d (cc.a,x.a,z.a,det);
  vold_3d (dd.a,x.a,y.a,det);
  
  reallocv (RSTCKVEC(intordkm[ri][ci],gp1));
  reallocv (RSTCKVEC(intordkm[ri][ci],gp2));
  reallocv (RSTCKVEC(intordkm[ri][ci],gp3));
  reallocv (RSTCKVEC(intordkm[ri][ci],w));

  gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordkm[ri][ci]);

  volcoord[0]=gp1[0];
  volcoord[1]=gp2[0];
  volcoord[2]=gp3[0];

  volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];

  //  matrix of gradients
  grad_matrix (gm,bb,cc,dd);
  
  //  matrix of conductivity of the material
  Tm->matcond (d,ipp,ri,ci);
  
  jac=fabs(det)/6.0;
  
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

   TKr, 05/04/2011
*/
void lintett::l_t_matrix (long /*lcid*/,long eid,long ri,long ci,matrix &lm)
{
  long ipp;
  double jac,det;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne)),bb(ASTCKVEC(4)),cc(ASTCKVEC(4)),dd(ASTCKVEC(4)),volcoord(ASTCKVEC(4));
  vector gp1,gp2,gp3,w;
  matrix gm(ASTCKMAT(ncomp,dofe[ri][ci])),d(ASTCKMAT(ncomp,ncomp));
  
  matrix n(ASTCKMAT(1,dofe[ri][ci]));

  Tt->give_node_coord3d (x,y,z,eid);
  
  nullm (lm);

  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];

  det = det3d (x.a,y.a,z.a);
  volb_3d (bb.a,y.a,z.a,det);
  volc_3d (cc.a,x.a,z.a,det);
  vold_3d (dd.a,x.a,y.a,det);
  
  reallocv (RSTCKVEC(intordkm[ri][ci],gp1));
  reallocv (RSTCKVEC(intordkm[ri][ci],gp2));
  reallocv (RSTCKVEC(intordkm[ri][ci],gp3));
  reallocv (RSTCKVEC(intordkm[ri][ci],w));

  gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordkm[ri][ci]);

  volcoord[0]=gp1[0];
  volcoord[1]=gp2[0];
  volcoord[2]=gp3[0];

  volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];

  //  matrix of gradients
  grad_matrix (gm,bb,cc,dd);
  
  //  matrix of conductivity of the material
  Tm->matcond (d,ipp,ri,ci);
  
  jac=fabs(det)/6.0;
  
  //  contribution to the L matrix of the element
  mtxm(gm,d,lm);
  cmulm(jac,lm);

}




/**
   function computes capacity %matrix of one transported matter
   finite element with tri-linear approximation functions
   
   @param eid - number of element
   @param ri,ci - row and column indeces of the computed block in the resulting %matrix
   @param cm - capacity %matrix
   
   JK, 24.3.2002
*/
void lintett::capacity_matrix (long eid,long ri,long ci,matrix &cm)
{
  long i,ii;
  double jac,det,rho,c;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne)),w(ASTCKVEC(intordcm[ri][ci])),gp1(ASTCKVEC(intordcm[ri][ci]));
  vector gp2(ASTCKVEC(intordcm[ri][ci])),gp3(ASTCKVEC(intordcm[ri][ci]));
  vector dens(ASTCKVEC(nne)),volcoord(ASTCKVEC(4));
  matrix n(ASTCKMAT(1,dofe[ri][ci]));
  
  Tt->give_elemnodes (eid,nodes);
  Tc->give_density (eid,nodes,dens);
  
  Tt->give_node_coord3d (x,y,z,eid);
  gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordcm[ri][ci]);
  det = det3d (x.a,y.a,z.a);
  
  nullm (cm);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci]+intordkm[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0]+intordkm[0][0];

  for (i=0;i<intordcm[ri][ci];i++){
    volcoord[0]=gp1[i];
    volcoord[1]=gp2[i];
    volcoord[2]=gp3[i];
    volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
    
    bf_matrix (n,volcoord);
    
    c=Tm->capcoeff (ii,ri,ci);
    rho = approx (volcoord,dens);
    jac=fabs(det)*w[i]*rho*c;
    
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
   @param ri,ci - row and column indeces of the block (ri must be equal to ci)
   
   JK, 24.3.2002
*/
void lintett::quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci)
{
  long i;
  double jac,det;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne)),w(ASTCKVEC(intordcm[ri][ci]));
  vector gp1(ASTCKVEC(intordcm[ri][ci])), gp2(ASTCKVEC(intordcm[ri][ci])),gp3(ASTCKVEC(intordcm[ri][ci])),dens(ASTCKVEC(nne));
  vector volcoord(ASTCKVEC(4)),v(ASTCKVEC(dofe[ri][ci]));
  matrix n(ASTCKMAT(1,dofe[ri][ci])),nm(ASTCKMAT(dofe[ri][ci],dofe[ri][ci]));
  
  Tt->give_node_coord3d (x,y,z,eid);
  gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordcm[ri][ci]);
  det = det3d (x.a,y.a,z.a);
  
  nullm (nm);
  
  for (i=0;i<intordcm[ri][ci];i++){
    volcoord[0]=gp1[i];
    volcoord[1]=gp2[i];
    volcoord[2]=gp3[i];
    volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
    
    bf_matrix (n,volcoord);
    
    jac=w[i]*fabs(det);
    
    nnj (nm.a,n.a,jac,n.m,n.n);
    
  }
  mxv (nm,nodval,v);  addv (sv,v,sv);
  
}


/**
   function computes internal fluxes
   
   @param lcid - load case id
   @param eid - element id
   @param ifl - %vector of internal fluxes
   
   JK, 3.12.2002
*/
void lintett::internal_fluxes (long lcid,long eid,vector &ifl)
{
  long ipp;
  double det,jac;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  vector bb(ASTCKVEC(4)),cc(ASTCKVEC(4)),dd(ASTCKVEC(4)),fl(ASTCKVEC(ncomp)),contr(ASTCKVEC(dofe[lcid][lcid]));
  matrix gm(ASTCKMAT(ncomp,dofe[lcid][lcid]));
  
  Tt->give_elemnodes (eid,nodes);
  Tt->give_node_coord3d (x,y,z,eid);
  
  det = det3d (x.a,y.a,z.a);
  volb_3d (bb.a,y.a,z.a,det);
  volc_3d (cc.a,x.a,z.a,det);
  vold_3d (dd.a,x.a,y.a,det);

  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[lcid][lcid];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];

  Tm->computenlfluxes (lcid,ipp);
  
  //Tm->givefluxes (lcid,ipp,0,fl);
  Tm->givefluxes (lcid,ipp,fl);

  //  matrix of gradients
  grad_matrix (gm,bb,cc,dd);
  
  mtxv (gm,fl,contr);
  
  jac=fabs(det)/6.0;
  cmulv (jac,contr);

  addv (contr,ifl,ifl);
  
}

/**
   function assembles resulting element conductivity %matrix

   @param eid - element id
   @param lcid - load case id
   @param km - resulting conductivity %matrix of one element

   JK, 6.1.2002
*/
void lintett::res_conductivity_matrix (long eid,long /*lcid*/,matrix &km)
{
  long i,j;
  ivector rcn, ccn;
  matrix lkm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      reallocv(RSTCKIVEC(dofe[i][j], rcn));
      reallocv(RSTCKIVEC(dofe[i][j], ccn));
      //rcn = new long [dofe[i][j]];
      //ccn = new long [dofe[i][j]];
      reallocm (RSTCKMAT(dofe[i][j],dofe[i][j],lkm));
      conductivity_matrix (i,eid,i,j,lkm);
      codnum (rcn,i);
      codnum (ccn,j);
      mat_localize (km,lkm,rcn.a,ccn.a);
      //delete [] rcn;  delete [] ccn;
    }
  }
}


/**
   function assembles resulting element L %matrix

   L^T = \int_{\Omega} B^T D {\rm d} \Omega

   @param eid - element id
   @param lcid - load case id
   @param lm - resulting L %matrix of one element

   TKr, 05/04/2011
*/
void lintett::res_l_matrix (long eid,long /*lcid*/,matrix &lm)
{
  long i,j;
  ivector rcn, ccn;
  matrix lkm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      reallocv(RSTCKIVEC(3, rcn));
      reallocv(RSTCKIVEC(dofe[i][j], ccn));
      //rcn = new long [3];
      //ccn = new long [dofe[i][j]];
      if (i==0){
	rcn[0]=1;
	rcn[1]=2;
	rcn[2]=3;
      }
      if (i==1){
	rcn[0]=4;
	rcn[1]=5;
	rcn[2]=6;
      }
      reallocm (RSTCKMAT(ncomp,dofe[i][j],lkm));
      l_matrix (i,eid,i,j,lkm);
      codnum (ccn,j);
      
      mat_localize (lm,lkm,rcn.a,ccn.a);
      //delete [] rcn;  delete [] ccn;
    }
  }
}



/**
   function assembles resulting element L^T (L transposed) %matrix

   L^T = \int_{\Omega} B^T D {\rm d} \Omega

   @param eid - element id
   @param lcid - load case id
   @param lm - resulting L^T %matrix of one element

   TKr, 05/04/2011
*/
void lintett::res_l_t_matrix (long eid,long /*lcid*/,matrix &lm)
{
  long i,j;
  ivector rcn, ccn;
  matrix lkm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      reallocv(RSTCKIVEC(dofe[i][j], rcn));
      reallocv(RSTCKIVEC(3, ccn));
      //rcn = new long [dofe[i][j]];
      //ccn = new long [3];
      if (j==0){
	ccn[0]=1;
	ccn[1]=2;
	ccn[2]=3;
      }
      if (j==1){
	ccn[0]=4;
	ccn[1]=5;
	ccn[2]=6;
      }
      reallocm (RSTCKMAT(dofe[i][j],ncomp,lkm));
      l_t_matrix (i,eid,i,j,lkm);
      codnum (rcn,i);
      
      mat_localize (lm,lkm,rcn.a,ccn.a);
      //delete [] rcn;  delete [] ccn;
    }
  }
}





/**
   function assembles average D %matrix

   @param eid - element id
   @param lm - resulting D %matrix of one element

   TKr, 05/04/2011
*/
void lintett::averd_matrix (long eid,matrix &lm)
{
  long i,j,ii;
  ivector rcn(ASTCKIVEC(ncomp)),ccn(ASTCKIVEC(ncomp));
  matrix d(ASTCKMAT(ncomp,ncomp));
  
  //rcn = new long [ncomp];
  //ccn = new long [ncomp];

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
	rcn[2]=3;
      }
      if (i==1){
	rcn[0]=4;
	rcn[1]=5;
	rcn[2]=6;
      }

      if (j==0){
	ccn[0]=1;
	ccn[1]=2;
	ccn[2]=3;
      }
      if (j==1){
	ccn[0]=4;
	ccn[1]=5;
	ccn[2]=6;
      }

      mat_localize (lm,d,rcn.a,ccn.a);
    }
  }
  
  //delete [] rcn;  delete [] ccn;
}





/**
   function assembles average C %matrix

   @param eid - element id
   @param lm - resulting C %matrix of one element

   TKr, 05/04/2011
*/
void lintett::averc_matrix (long eid,matrix &lm)
{
  long i,j,ii;
  double c,rho;
  ivector nodes(ASTCKIVEC(nne));
  vector dens(ASTCKVEC(nne)),volcoord(ASTCKVEC(4));
  
  Tt->give_elemnodes (eid,nodes);
  Tc->give_density (eid,nodes,dens);
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      
      if (Tp->savemode==0)
	ii=Tt->elements[eid].ipp[i][j];
      if (Tp->savemode==1)
	ii=Tt->elements[eid].ipp[0][0];
      
      //temporarily
      volcoord[0] = 0.0;
      volcoord[1] = 0.0;
      volcoord[2] = 0.0;
      volcoord[3] = 1.0-volcoord[0]-volcoord[1]-volcoord[2];

      //  coefficient of capacity of material
      c = Tm->capcoeff (ii,i,j);
      rho = approx (volcoord,dens);

      lm[i][j] = rho*c*Tm->ip[ii].av[j];
    }
  }
}




/**
   function assembles volume of one element

   @param eid - element id

   TKr, 05/04/2011
*/
double lintett::elem_volume (long eid)
{
  double det,volume;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  
  Tt->give_node_coord3d (x,y,z,eid);
  
  //  det is equal to six times volume of the element
  det = det3d (x.a,y.a,z.a);

  if (det<0.0){
    //print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
    det=fabs(det);
  }

  volume = det/6.0;

  return volume;
}


/**
   function assembles resulting element capacity %matrix
   
   @param eid - element id
   @param cm - resulting capacity %matrix of one element

   JK, 6.1.2002
*/
void lintett::res_capacity_matrix (long eid,matrix &cm)
{
  long i,j;
  ivector rcn,ccn;
  matrix lcm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      reallocv(RSTCKIVEC(dofe[i][j], rcn));
      reallocv(RSTCKIVEC(dofe[i][j], ccn));
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
      mat_localize (cm,lcm,rcn.a,ccn.a);
    }
  }
}

/**
   function assembles resulting element convection %vector

   @param f - resulting convection %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   
   JK, 6.1.2002, 18. 4. 2014
*/
void lintett::res_convection_vector (vector &f,long lcid,long eid,long leid)
{
  ivector cn;
  vector lf;

  //  transi[lcid]==2 - element contains boundary with prescribed flux
  //  transi[lcid]==4 - element contains boundary with prescribed transmission and prescribed flux
  if ((Tt->elements[eid].transi[lcid]==2)||(Tt->elements[eid].transi[lcid]==4)){
    //  array for code numbers
    reallocv(RSTCKIVEC(dofe[lcid][lcid], cn));
    reallocv(RSTCKVEC(dofe[lcid][lcid],lf));
    
    convection_vector (lf,lcid,eid,leid,lcid,lcid);
    
    //  code numbers
    codnum (cn,lcid);
    //  localization to element vector
    locglob (f.a,lf.a,cn.a,dofe[lcid][lcid]);
    
    cmulv (-1.0,f,f);
  }
}

/**
   function assembles resulting element transmission %vector

   @param f - resulting transmission %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   
   JK, 6.1.2002, 18. 4. 2014
*/
void lintett::res_transmission_vector (vector &f,long lcid,long eid,long leid)
{
  long i;
  ivector cn;
  vector lf;
  
  reallocv(RSTCKIVEC(dofe[lcid][lcid], cn));
  //  code numbers
  codnum (cn,lcid);
  reallocv (RSTCKVEC(dofe[lcid][lcid],lf));
  
  //  transi[i]==3 - element contains boundary with prescribed transmission
  //  transi[i]==4 - element contains boundary with prescribed transmission and prescribed flux
  if ((Tt->elements[eid].transi[lcid]==3)||(Tt->elements[eid].transi[lcid]==4)){
    
    for (i=0;i<ntm;i++){    
      nullv (lf);
      
      transmission_vector (lf,lcid,eid,leid,i);
      
      locglob (f.a,lf.a,cn.a,dofe[lcid][lcid]);
    }
  }
}

/**
   function assembles resulting element source vector

   @param sv - resulting source %vector of one element
   @param lcid - load case id
   @param eid - element id
   
   JK, 6.1.2002
*/
void lintett::res_quantity_source_vector (vector &sv,vector &nodval,long lcid,long eid)
{
  ivector cn;
  vector lsv;

  reallocv(RSTCKIVEC(dofe[lcid][lcid], cn));
  reallocv (RSTCKVEC(dofe[lcid][lcid],lsv));
  quantity_source_vector (lsv,nodval,eid,lcid,lcid);
  codnum (cn,lcid);
  locglob (sv.a,lsv.a,cn.a,dofe[lcid][lcid]);
}

/**
   function assembles resulting element internal fluxes vector

   @param eid - element id
   @param elemif - resulting internal fluxes %vector of one element
   
   JK, 6.1.2002
*/
void lintett::res_internal_fluxes (long eid,vector &elemif)
{
  long i;
  ivector cn;
  vector lif,tdnv;
  matrix cm;

  for (i=0;i<ntm;i++){
    reallocv (RSTCKIVEC(dofe[i][i], cn));
    reallocv (RSTCKVEC(dofe[i][i],lif));
    internal_fluxes (i,eid,lif);
    codnum (cn, i);
    locglob (elemif.a, lif.a, cn.a, dofe[i][i]);
  }

  reallocv(RSTCKIVEC(ndofe, cn));
  Gtt->give_code_numbers (eid, cn.a);
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
double lintett::total_integral(long eid,vector &nodval)
{
  long i;
  double jac,det,value,f;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  vector w(ASTCKVEC(4)),gp1(ASTCKVEC(4)),gp2(ASTCKVEC(4)),gp3(ASTCKVEC(4));
  vector dens(ASTCKVEC(nne)),volcoord(ASTCKVEC(4));
  
  Tt->give_elemnodes (eid,nodes);
  
  Tt->give_node_coord3d (x,y,z,eid);
  gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,4);
  det = det3d (x.a,y.a,z.a);
  
  f = 0.0;
  
  for (i=0;i<4;i++){
    volcoord[0]=gp1[i];
    volcoord[1]=gp2[i];
    volcoord[2]=gp3[i];
    volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
    
    value = approx (volcoord,nodval);
    jac=fabs(det)*w[i]*value;
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
void lintett::res_boundary_flux (vector &f,long lcid,long eid,long leid)
{
  long i;
  ivector cn;
  vector lf;

  reallocv(RSTCKIVEC(dofe[lcid][lcid], cn));
  codnum (cn,lcid);
  
  //coupled b.c.
  for (i=0;i<ntm;i++){    
    reallocv (RSTCKVEC(dofe[lcid][lcid],lf));
    
    if ((Tt->elements[eid].transi[i]==3)||(Tt->elements[eid].transi[i]==4)){
      boundary_flux (lf,i,eid,leid,lcid,i);
    }
    
    locglob (f.a,lf.a,cn.a,dofe[lcid][lcid]);
  }
}



/**
   function computes fluxes in integration points

   @param eid - element id

   JK, 12.8.2014, revised by TKr 29/09/2022
*/
void lintett::intpointflux (long eid)
{
  long i,l,ii,jj,ipp;
  
  for (l=0;l<Tp->ntm;l++){
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
	
	ipp=Tt->elements[eid].ipp[ii][jj];
	
	for (i=0;i<intordkm[ii][jj];i++){
	  //  computation of correct fluxes
	  if (Tp->fluxcomp==1)
	    Tm->computenlfluxes (l,ipp);
	  
	  ipp++;
	}
	for (i=0;i<intordcm[ii][jj];i++){
	  //  computation of correct fluxes
	  if (Tp->fluxcomp==1)
	    Tm->computenlfluxes (l,ipp);
	  
	  ipp++;
	}
	if (Tp->savemode==1)  break;
      }
      if (Tp->savemode==1)  break;
    }
  }
  
  /*  
      long i,j,k,l,ll,ii,jj,ipp;
      vector grad(ASTCKVEC(ncomp)),fl(ASTCKVEC(ncomp)),cfl(ASTCKVEC(ncomp));
      matrix d(ASTCKMAT(ncomp,ncomp));
      
      //  loop over the number of media in the problem
      for (l=0;l<Tp->ntm;l++){
      
      //  loop over all integration points on element
      for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
      
      //  id of the actual integration point
      ipp=Tt->elements[eid].ipp[ii][jj];
      for (i=0;i<intordkm[ii][jj];i++){
      for (j=0;j<intordkm[ii][jj];j++){
      for(k=0;k<intordkm[ii][jj];k++){
      
      //  loop over the number of media
      //  flux in every point is the sum of contribution from each medium
      for (ll=0;ll<Tp->ntm;ll++){
      //  gradients
      Tm->givegrad (ll,ipp,grad);
      //  block of conductivity matrix of material
      Tm->matcond (d,ipp,l,ll);
      //  flux caused by the actual gradient
      mxv (d,grad,fl);
      //  summation of flux contributions
      addv (fl,cfl,cfl);
      }
      Tm->storeflux (l,ipp,cfl);
      ipp++;
      }
      }
      }
      }
      if (Tp->savemode==1)
      break;
      }
      }
  */
}


/**
   function computes gradients in nodes of element

   @param eid - element id
   
   24. 3. 2017, JK
*/
void lintett::nod_grads_ip (long eid)
{
  long i,j,k,ipp;
  double vol;
  ivector ipnum(ASTCKIVEC(nne)),nod(ASTCKIVEC(nne));
  vector grad(ASTCKVEC(ncomp));

  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  nodip_lintet (ipp,ipnum);
  
  //  node numbers of the element
  Tt->give_elemnodes (eid,nod);
  // element volume
  if (Tp->gradaver == 2)
    vol = element_volume(eid);
  
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
	j=nod[i];
	Tt->nodes[j].storegrad (k,vol,grad);
      }
    }
  }
  
}


/**
   function computes fluxes in nodes of element

   @param eid - element id
   
   24. 3. 2017, JK
*/
void lintett::nod_fluxes_ip (long eid)
{
  long i,j,k,ipp;
  double vol;
  ivector ipnum(ASTCKIVEC(nne)),nod(ASTCKIVEC(nne));
  vector flux(ASTCKVEC(ncomp));
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  nodip_lintet (ipp,ipnum);
  
  //  node numbers of the element
  Tt->give_elemnodes (eid,nod);
  // element volume
  if (Tp->fluxaver == 2)
    vol = element_volume(eid);
  
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
	j=nod[i];
	Tt->nodes[j].storeflux (k,vol,flux);
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
void lintett::nod_others_comp (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp,ncompother,ii;
  vector h;
  vector r(ASTCKVEC(ndofe));
  ivector nod(ASTCKIVEC(nne));
  double vol,other;
  
  Tt->give_elemnodes (eid,nod);  
  elemvalues(eid, r);
  ncompother = Tm->givencompother();
  reallocv (RSTCKVEC(ntm,h));
  // element volume
  if (Tp->otheraver == 2)
    vol = element_volume(eid);
  
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
	Tt->nodes[j].storeother (k,vol,other);
      }
    }
  }
  
}


/**
   function computes other values in nodes of element

   @param eid - element id
   
   24. 3. 2017, JK
*/
void lintett::nod_eqother_ip (long eid)
{
  long i,j,k,ipp,ncompo;
  double vol;
  ivector ipnum(ASTCKIVEC(nne)),nod(ASTCKIVEC(nne));
  vector eqother;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  nodip_lintet (ipp,ipnum);
  
  //  node numbers of the element
  Tt->give_elemnodes (eid,nod);
  // element volume
  if (Tp->eqotheraver == 2)
    vol = element_volume(eid);
  
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
	Tt->nodes[j].storeeqother (k,vol,eqother(k));
      }
    }
  }
}





// This bellow added by TKr 18/03/2008
/**
   function computes contribution to the flux %vector
   
   \int_{Gamma_2} N^T N dGamma * nodal_flux_values
   
   @param v - %vector of node fluxes
   @param lcid - load case id
   @param eid - element id
   @param leid - loaded element id
   @param ri,ci - row and column indices
   
   JK, 8.10.2001, TKr, 18.3.2008
*/
void lintett::convection_vector (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordb)),gp1(ASTCKVEC(intordb)),gp2(ASTCKVEC(intordb));
  vector list(ASTCKVEC(3*nsurf)),trc(ASTCKVEC(3*nsurf)),nodval(ASTCKVEC(nne));
  vector av(ASTCKVEC(nne)),coef(ASTCKVEC(nne));
  matrix km(ASTCKMAT(nne,nne));
  
  nullv (v);
  
  //  node coordinates on element
  Tt->give_node_coord3d (x,y,z,eid);
  
  //  coordinates and weights of integration points
  gauss_points_tr (gp1.a,gp2.a,w.a,intordb);
  
  //  indicators of boundary conditions
  bc = new bocontypet [nsurf];
  Tb->lc[lcid].elemload[leid].give_bc (bc);
  
  //  auxiliary coefficients, necessary for function surface_integral
  for (i=0;i<3*nsurf;i++){
    trc[i]=1.0;
  }

  //  nodal values of exterior variable
  Tb->lc[lcid].elemload[leid].give_nodval (lcid,list);
  
  //  integration point number
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  //  loop over surfaces
  for (i=0;i<nsurf;i++){
    
    if (bc[i]==2 || bc[i]==3 || bc[i]==5){
      //  nodal values on actual surface
      surfnodeval (i,nodval,list);
      
      nullm (km);

      surfnodeval (i,coef,trc);

      //  matrix obtained from integration over surface
      surface_integral (i,x,y,z,intordb,gp1,gp2,w,coef,km);
      
      mxv (km,nodval,av);
      
      addv (v,av,v);
    }
  }
  
  delete [] bc;
}

/**
   function computes transmission complement to the conductivity %matrix for one matter
   
   \int_{Gamma_3} N^T c_{tr} N dGamma
   
   @param lcid - number of load case
   @param eid - number of element
   @param ri,ci - row and column indeces
   @param km - part of the conductivity %matrix
   
   JK, 24.3.2002
   TKr, 30.1.2004 - new added
*/
void lintett::transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,j,ipp,leid;
  bocontypet *bc;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordb)),gp1(ASTCKVEC(intordb)),gp2(ASTCKVEC(intordb));
  vector trc(ASTCKVEC(3*nsurf)),coeff(ASTCKVEC(nne));
  
  //  node coordinates on element
  Tt->give_node_coord3d (x,y,z,eid);
  
  //  coordinates and weights of integration points
  gauss_points_tr (gp1.a,gp2.a,w.a,intordb);

  //  correspondence between element ordering and boundary element ordering
  for (i=0;i<Tb->lc[lcid].neb;i++){
    if (Tb->lc[lcid].elemload[i].eid==eid){
      leid=i;
      //break;//debug for St. Ann's church
      //}
      //}

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
      for (j=0;j<nsurf;j++){
	
	if (bc[j]==5 || bc[j]>10){
	  //  transformation of coefficients
	  transf_coeff (j,coeff,trc,eid,ri,ci,ipp,bc);
	  
	  //  matrix obtained from integration over surface
	  surface_integral (j,x,y,z,intordb,gp1,gp2,w,coeff,km);
	}
      }
      
      delete [] bc;
    }////debug for St. Ann's church
  }////debug for St. Ann's church
}

/**
   function computes contributions to the transmission vector
   
   \int_{Gamma_3} N^T c_{tr} N dGamma * nodal_external_temperature
   
   @param tmv - transmission %vector of one matter
   @param lcid - load case id (corresponds to the row index)
   @param eid - element id
   @param leid - loaded element id
   @param cid - component id (corresponds to the column index)
      
   JK, 24.3.2002
   TKr, 30.1.2004 - new added
*/
void lintett::transmission_vector (vector &v,long lcid,long eid,long leid,long cid)
{
  long i,ipp;
  bocontypet *bc;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordb)),gp1(ASTCKVEC(intordb)),gp2(ASTCKVEC(intordb));
  vector list(ASTCKVEC(3*nsurf)),trc(ASTCKVEC(3*nsurf)),trr(ASTCKVEC(3*nsurf));
  vector nodval(ASTCKVEC(nne)),av(ASTCKVEC(nne)),coef(ASTCKVEC(nne));
  matrix km(ASTCKMAT(nne,nne));
  
  nullv (v);
  
  //  node coordinates on element
  Tt->give_node_coord3d (x,y,z,eid);
  
  //  coordinates and weights of integration points
  gauss_points_tr (gp1.a,gp2.a,w.a,intordb);
  
  //  indicators of boundary conditions
  bc = new bocontypet [nsurf];
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
  
  for (i=0;i<nsurf;i++){
    
    if (bc[i]==5 || bc[i]>10){
      //  transformation of nodal values
      transf_val (i,nodval,list,trc,trr,eid,lcid,cid,ipp,bc);
      //  determination of transmission coefficients on edge
      transf_coeff (i,coef,trc,eid,lcid,cid,ipp,bc);

      nullm (km);
      //  matrix obtained from integration over surfaces
      surface_integral (i,x,y,z,intordb,gp1,gp2,w,coef,km);
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

   @param v - %vector of boundary fluxes
   @param lcid - load case id
   @param eid - element id
   @param leid - loaded element id
   @param ri,ci - row and column indices
   
   TKr, 28.2.2004
*/
void lintett::boundary_flux (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordb)),gp1(ASTCKVEC(intordb)),gp2(ASTCKVEC(intordb));
  vector list(ASTCKVEC(3*nsurf)),trc(ASTCKVEC(3*nsurf)),trr(ASTCKVEC(3*nsurf));
  vector nodval(ASTCKVEC(nne)),av(ASTCKVEC(nne)),coef(ASTCKVEC(nne));
  matrix km(ASTCKMAT(nne,nne));
  
  nullv (v);
  
  //  node coordinates on element
  Tt->give_node_coord3d (x,y,z,eid);
  
  //  coordinates and weights of integration points
  gauss_points_tr (gp1.a,gp2.a,w.a,intordb);
  
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
    
    if (bc[i]==5 || bc[i]>10){
      //  transformation of nodal values
      transf_flux (i,nodval,list,trc,trr,eid,ri,ci,ipp,bc);
      //  determination of transmission coefficients on edge
      transf_coeff (i,coef,trc,eid,ri,ci,ipp,bc);

      nullm (km);
      //  matrix obtained from integration over surfaces
      surface_integral (i,x,y,z,intordb,gp1,gp2,w,coef,km);
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
   @param gp1,gp2, w  - coordinates and weights of integration points
   @param coef - array of nodal values of coefficient
   @param km - output %matrix
   
   JK
*/
void lintett::surface_integral (long surf,vector &x,vector &y,vector &z,long intord,vector &gp1,vector &gp2,vector &w,
                                vector &coef,matrix &km)
{
  long i;
  double jac,ipval;
  vector volcoord(ASTCKVEC(4));
  matrix n(ASTCKMAT(1,nne));
  
  if (surf==0){
    for (i=0;i<intord;i++){
      volcoord[0]=0.0;
      volcoord[1]=gp1[i];
      volcoord[2]=gp2[i];
      volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
      
      bf_matrix (n,volcoord);
      
      jac2d_3d (jac,x,y,z,gp1[i],gp2[i],0);
        
      ipval=approx (volcoord,coef);
      
      jac*=w[i]*ipval;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }
  
  if (surf==1){
    for (i=0;i<intord;i++){
      volcoord[0]=gp1[i];
      volcoord[1]=0.0;
      volcoord[2]=gp2[i];
      volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
      
      bf_matrix (n,volcoord);
      
      jac2d_3d (jac,x,y,z,gp1[i],gp2[i],1);
        
      ipval=approx (volcoord,coef);
      
      jac*=w[i]*ipval;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }
  
  if (surf==2){
    for (i=0;i<intord;i++){
      volcoord[0]=gp1[i];
      volcoord[1]=gp2[i];
      volcoord[2]=0.0;
      volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
      
      bf_matrix (n,volcoord);
      
      jac2d_3d (jac,x,y,z,gp1[i],gp2[i],2);
        
      ipval=approx (volcoord,coef);
      
      jac*=w[i]*ipval;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }
  
  if (surf==3){
    for (i=0;i<intord;i++){
      volcoord[0]=gp1[i];
      volcoord[1]=gp2[i];
      volcoord[3]=0.0;
      volcoord[2]=1.0-volcoord[0]-volcoord[1]-volcoord[3];
      
      bf_matrix (n,volcoord);
      
      jac2d_3d (jac,x,y,z,gp1[i],gp2[i],3);
        
      ipval=approx (volcoord,coef);
      
      jac*=w[i]*ipval;
      nnj (km.a,n.a,jac,n.m,n.n);
    }
  }
  

}

void lintett::transf_flux (long surf,vector &coeff,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double tr,new_flux;
  ivector nodes(ASTCKIVEC(nne)),surfnod(ASTCKIVEC(nnsurf));

  //  zeroing of array of coefficients
  nullv (coeff);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required surface
  lintetrahedral_surfnod (surfnod.a,surf);

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


void lintett::transf_coeff (long surf,vector &coeff,vector &list,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double new_coeff;
  ivector nodes(ASTCKIVEC(nne)),surfnod(ASTCKIVEC(nnsurf));

  //  zeroing of array of coefficients
  nullv (coeff);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required surface
  lintetrahedral_surfnod (surfnod.a,surf);

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
void lintett::transf_val (long surf,vector &nodval,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double tr,new_nodval;
  ivector nodes(ASTCKIVEC(nne)),surfnod(ASTCKIVEC(nnsurf));
  
  //  zeroing of array of nodal values
  nullv (nodval);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required surface
  lintetrahedral_surfnod (surfnod.a,surf);
  
  //  actual position in the array list
  j=surf*nnsurf;
  
  //  loop over number of nodes on one surface
  for (i=0;i<nnsurf;i++){
    tr = 0.0;
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
void lintett::surfnodeval (long surf,vector &nodval,vector &list)
{
  long i,j;
  ivector surfnod(ASTCKIVEC(nnsurf));
  
  nullv (nodval);
  lintetrahedral_surfnod (surfnod.a,surf);
  
  for (i=0;i<nnsurf;i++){
    j=surfnod[i];
    nodval[j]=list[surf*nnsurf+i];
  }
}


/**
   function selects values and gradients from the global level
   
   @param eid - element id
   @param counter - actual position in the array buff
   @param buff - array containing selected components
   
   JK, 23.8.2011
*/
void lintett::higher_to_lower_level (long eid,long *counter,double *buff)
{
  long i,ipp;
  vector gr(ASTCKVEC(3));
  
  //  id of integration point
  ipp=Tt->elements[eid].ipp[0][0];
  

  //fprintf (Outt,"\n tisk na prvku c.= %ld pred: \n",eid);
  //for (long ijk=0;ijk<4;ijk++){
  //fprintf (Outt,"\n buff1 %ld   %le",ijk,buff[ijk]);
  //}

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
    buff[counter[0]]=gr[2];
    counter[0]++;
    
  }

  //fprintf (Outt,"\n counter na konci kazdeho prvku = %ld po: \n",counter[0]);
  //fprintf (Outt,"\n tisk na prvku c.= %ld po: \n",eid);
  //for (long ijk=0;ijk<4;ijk++){
  //fprintf (Outt,"\n buff1 %ld   %le",ijk,buff[ijk]);
  //}

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
long lintett::ipcoord (long eid, long ipp, long ri, long ci, vector &coord)
{
  long i, ii;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne));
  vector volcoord(ASTCKVEC(4));
  vector w, gp1, gp2, gp3;

  Tt->give_node_coord3d(x, y, z, eid);
  
  //  integration points for the conductivity matrix
  reallocv(RSTCKVEC(intordkm[ri][ci], gp1));
  reallocv(RSTCKVEC(intordkm[ri][ci], gp2));
  reallocv(RSTCKVEC(intordkm[ri][ci], gp3));
  reallocv(RSTCKVEC(intordkm[ri][ci], w));
  gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordkm[ri][ci]);
	
  ii = Tt->elements[eid].ipp[ri][ci];
  for (i=0; i<intordkm[ri][ci]; i++)
  {
    if (ii == ipp)
    {
      volcoord[0]=gp1[i];
      volcoord[1]=gp2[i];
      volcoord[2]=gp3[i];
      volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
      coord[0] = approx(volcoord, x);
      coord[1] = approx(volcoord, y);
      coord[2] = approx(volcoord, z);
      return 0;
    }
    ii++;
  }
	
  //  integration points for the capacity matrix
  reallocv(RSTCKVEC(intordcm[ri][ci], gp1));
  reallocv(RSTCKVEC(intordcm[ri][ci], gp2));
  reallocv(RSTCKVEC(intordcm[ri][ci], gp3));
  reallocv(RSTCKVEC(intordcm[ri][ci], w));
  gauss_points_tet (gp1.a, gp2.a, gp3.a, w.a, intordcm[ri][ci]);
  
  for (i=0; i<intordcm[ri][ci]; i++)
  {
    if (ii == ipp)
    {
      volcoord[0]=gp1[i];
      volcoord[1]=gp2[i];
      volcoord[2]=gp3[i];
      volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
      coord[0] = approx(volcoord, x);
      coord[1] = approx(volcoord, y);
      coord[2] = approx(volcoord, z);
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
long lintett::ipncoord (long eid, long ipp, vector &ncoord)
{
  long i, ii, ri, ci;
  vector w, gp1, gp2, gp3;

  for (ri=0;ri<ntm;ri++)
  {
    for (ci=0;ci<ntm;ci++)
    {
      //  integration points for the conductivity matrix
      reallocv(RSTCKVEC(intordkm[ri][ci], gp1));
      reallocv(RSTCKVEC(intordkm[ri][ci], gp2));
      reallocv(RSTCKVEC(intordkm[ri][ci], gp3));
      reallocv(RSTCKVEC(intordkm[ri][ci], w));
      gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,intordkm[ri][ci]);
	
      ii = Tt->elements[eid].ipp[ri][ci];
      for (i=0; i<intordkm[ri][ci]; i++)
      {
        if (ii == ipp)
        {
          ncoord[0] = gp1[i];
          ncoord[1] = gp2[i];
          ncoord[2] = gp3[i];
          return 0;
        }
        ii++;
      }
	
      //  integration points for the capacity matrix
      reallocv(RSTCKVEC(intordcm[ri][ci], gp1));
      reallocv(RSTCKVEC(intordcm[ri][ci], gp2));
      reallocv(RSTCKVEC(intordcm[ri][ci], gp3));
      reallocv(RSTCKVEC(intordcm[ri][ci], w));
      gauss_points_tet (gp1.a, gp2.a, gp3.a, w.a, intordcm[ri][ci]);
  
      for (i=0; i<intordcm[ri][ci]; i++)
      {
        if (ii == ipp)
        {
          ncoord[0] = gp1[i];
          ncoord[1] = gp2[i];
          ncoord[2] = gp3[i];
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
   
   JK, 2. 3. 2018
*/
void lintett::surface_flux (long lcid,long eid,long beid,double *fluxes)
{

  long i,k,ipp,fid;
  double det,q,qn,area;
  bocontypet *bc;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  vector bb(ASTCKVEC(4)),cc(ASTCKVEC(4)),dd(ASTCKVEC(4)),volcoord(ASTCKVEC(4));
  vector r(ASTCKVEC(ndofe)),t(ASTCKVEC(nne)),grad(ASTCKVEC(ncomp)),flux(ASTCKVEC(ncomp)),fl(ASTCKVEC(ncomp)),n(ASTCKVEC(ncomp));
  matrix gm(ASTCKMAT(ncomp,nne)),d(ASTCKMAT(ncomp,ncomp));
  
  
  //  node coordinates
  Tt->give_node_coord3d (x,y,z,eid);
  //  all nodal values on element
  elemvalues(eid, r);
  
  //  id of integration point for determination of the D matrix
  ipp=Tt->elements[eid].ipp[0][0];
  
  //  determinant for the volume coordinates
  det = det3d (x.a,y.a,z.a);
  
  if(det < 0.0){
    //det = fabs(det);
    print_err("wrong numbering of nodes on 3D element number %ld, negative volume! det = %e",__FILE__,__LINE__,__func__,eid+1,det);
    abort();
  }
  
  //  volume coordinates
  volb_3d (bb.a,y.a,z.a,det);
  volc_3d (cc.a,x.a,z.a,det);
  vold_3d (dd.a,x.a,y.a,det);
  
  //  matrix of gradients
  //  it is a constant matrix because the approximation functions are linear
  //  therefore, the gradients and fluxes are also constant everywhere in the element
  grad_matrix (gm,bb,cc,dd);
  
  nullv (flux);
  
  for (k=0;k<Tp->ntm;k++){
    //  loop over the number of variables in the problem
    
    //  nodal values of the k-th variable
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }
    
    //  gradient of the k-th variable
    mxv (gm,t,grad);
    //  conductivity matrix of the material
    Tm->matcond (d,ipp,lcid,k);
    //  flux caused by the k-th gradient
    mxv (d,grad,fl);
    //  the total flux
    addv (fl,flux,flux);
  }//  end of the loop over the number of variables in the problem
  
  //  indicators of boundary conditions
  bc = new bocontypet [nsurf];
  Tb->bf[lcid].elemload[beid].give_bc (bc);
  
  //  loop over surfaces
  for (i=0;i<nsurf;i++){
    
    if (bc[i]==1){
      //  function constructs the outer normal vector
      tetrahedra_normal_vectors (i,x,y,z,n);
      //  area of the i-th element surface
      area = tetrahedra_surface_areas (i,x,y,z);
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



/**
   Function returns transport (non-mechanical) quantities at nodes of element.
   The values of selected quantity is copied from the closest integration points 
   to element nodes.

   @param eid - element id
   @param nodval - %vector of nodal values
   @param ntq - type of non-mechanical quantity
   
   @return The function does not return anything.
   
   Created by TKo, 3.2.2014
*/
void lintett::transq_nodval (long eid,vector &nodval,nonmechquant nmq)
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
   
   Modified by TKo, 31.1.2014
   
   Created by 26/09/2012 TKr
*/
void lintett::transq_nodval_comp (long eid,vector &nodval, long ncne, long nq, nonmechquant *qt)
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
void lintett::transq_init_nodval (long eid,vector &nodval,nonmechquant nmq)
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
   @param nodval - %vector of initial nodal values of all required quantities, i.e., 
                   nodal value of i-th quantity in j-th node is given by nodval[i*ncnv+j] where ncnv 
                   is the number of calculated nodes on eid-th element.
   @param ncne - number of computed nodes on element (only first ncne of nodes is calculated)
   @param nq - number of required transport quantities
   @param qt - array of types of required non-mechanical, i.e. transport quantities
   
   @return The function does not return anything.
   
   Modified by TKo, 31.1.2014
   
   Created by 26/09/2012 TKr
*/
void lintett::transq_init_nodval_comp (long eid,vector &nodval, long ncne, long nq, nonmechquant *qt)
{
  long i,j,ipid,ncompo;
  ivector ipnum(ASTCKIVEC(nne)), enod(ASTCKIVEC(nne));
  vector ipav;
  intpointst ipb; // backup of the first integration point of element
    
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

