/*
  File:     quadlineart.cpp
  Author:   Jaroslav Kruis, 31.3.2001
  Purpose:  twodimensional quadrilateral element with linear approximation functions
*/

#include "globalt.h"
#include "quadlineart.h"
#include "genfile.h"
#include "globmatt.h"

quadlineart::quadlineart (void)
{
  long i;
  
  //  number of nodes on element
  nne=4;
  //  number of edges
  ned=4;
  //  number of nodes on one edge
  nned=2;
  //  geometrical dimension
  ncomp=2;
  
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
  case nomedium:{  break; }
  case onemedium:{
    ordering[0][0]=1;  ordering[0][1]=2;  ordering[0][2]=3;  ordering[0][3]=4;
    dofe[0][0]=4;  intordkm[0][0]=2;  intordcm[0][0]=2;
    nip[0][0]=8;
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
  case fourmediacoup:{
    ordering[0][0]=1;   ordering[0][1] =5;  ordering[0][2] =9;   ordering[0][3] =13;
    ordering[1][0]=2;   ordering[1][1] =6;  ordering[1][2] =10;  ordering[1][3] =14;
    ordering[2][0]=3;   ordering[2][1] =7;  ordering[2][2] =11;  ordering[2][3] =15;
    ordering[3][0]=4;   ordering[3][1] =8;  ordering[3][2] =12;  ordering[3][3] =16;

    intordkm[0][0]=2;  intordkm[0][1]=2;  intordkm[0][2]=2;  intordkm[0][3]=2;
    intordkm[1][0]=2;  intordkm[1][1]=2;  intordkm[1][2]=2;  intordkm[1][3]=2;
    intordkm[2][0]=2;  intordkm[2][1]=2;  intordkm[2][2]=2;  intordkm[2][3]=2;
    intordkm[3][0]=2;  intordkm[3][1]=2;  intordkm[3][2]=2;  intordkm[3][3]=2;

    intordcm[0][0]=2;  intordcm[0][1]=2;  intordcm[0][2]=2;  intordcm[0][3]=2;
    intordcm[1][0]=2;  intordcm[1][1]=2;  intordcm[1][2]=2;  intordcm[1][3]=2;
    intordcm[2][0]=2;  intordcm[2][1]=2;  intordcm[2][2]=2;  intordcm[2][3]=2;
    intordcm[3][0]=2;  intordcm[3][1]=2;  intordcm[3][2]=2;  intordcm[3][3]=2;
    
    if (Tp->savemode==0){
      nip[0][0]=8;  nip[0][1]=8;  nip[0][2]=8;  nip[0][3]=8;
      nip[1][0]=8;  nip[1][1]=8;  nip[1][2]=8;  nip[1][3]=8;
      nip[2][0]=8;  nip[2][1]=8;  nip[2][2]=8;  nip[2][3]=8;
      nip[3][0]=8;  nip[3][1]=8;  nip[3][2]=8;  nip[3][3]=8;
    }
    if (Tp->savemode==1){
      nip[0][0]=8;  nip[0][1]=0;  nip[0][2]=0;  nip[0][3]=0;
      nip[1][0]=0;  nip[1][1]=0;  nip[1][2]=0;  nip[1][3]=0;
      nip[2][0]=0;  nip[2][1]=0;  nip[2][2]=0;  nip[2][3]=0;
      nip[3][0]=0;  nip[3][1]=0;  nip[3][2]=0;  nip[3][3]=0;
    }
    
    dofe[0][0]=4;  dofe[0][1]=4;  dofe[0][2]=4;  dofe[0][3]=4;
    dofe[1][0]=4;  dofe[1][1]=4;  dofe[1][2]=4;  dofe[1][3]=4;
    dofe[2][0]=4;  dofe[2][1]=4;  dofe[2][2]=4;  dofe[2][3]=4;
    dofe[3][0]=4;  dofe[3][1]=4;  dofe[3][2]=4;  dofe[3][3]=4;

    ndofe=16;  napfun=4;
    break;
  }
  default:{
    print_err("unknown number of transported matters is required",__FILE__,__LINE__,__func__);
  }
  }
}

quadlineart::~quadlineart (void)
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
   function assembles element code numbers
   they are used for localization between all element unknowns and unknowns related to one matter
   
   @param cn - code numbers
   @param ri - number of matter (usually row index)
*/
void quadlineart::codnum (long *cn,long ri)
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
double quadlineart::element_area (long eid)
{
  double area;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));

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

   JK, 8.5.2002
*/
double quadlineart::approx (double xi,double eta,vector &nodval)
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
void quadlineart::ipcoordblock (long eid,long ri,long ci,double **coord)
{
  long i,j,k;
  double xi,eta;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordkm[ri][ci])), gp(ASTCKVEC(intordkm[ri][ci]));
  
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
   (function interpolates the nodal values to integration points)

   eid - element id

   JK, 25.9.2001
*/
void quadlineart::intpointval (long eid)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,val;
  vector r(ASTCKVEC(ndofe)), t(ASTCKVEC(nne)), gp, w;
  
  //  nodal values on element
  elemvalues(eid, r);
  
  for (k=0;k<Tp->ntm;k++){
    
    //  nodal values of a single variable
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1]; 
    }
    
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
	    
	    val = approx (xi,eta,t);
	    Tm->ip[ipp].av[k]=val;
	    ipp++;
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

   12/06/2012 TKr according to JK
*/
void quadlineart::intpointval (long eid,vector &nodval,vector &ipval)
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
	  kk++;
        }
      }
	
      if (Tp->savemode==1)  break;
    }
    if (Tp->savemode==1)  break;
  }
}


/**
   function computes values in integration points from nodal values for PUC
   (function interpolates the nodal values to integration points)

   eid - element id

   TKr, 05/04/2010
*/
void quadlineart::intpointval_puc (long eid)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,val;
  vector r(ASTCKVEC(ndofe)),t(ASTCKVEC(nne)),gp,w;
  
  //  nodal values on element
  elemvalues_puc(eid, r);
  
  for (k=0;k<Tp->ntm;k++){
    
    //  nodal values of a single variable
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }
    
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
	    
	    val = approx (xi,eta,t);
	    Tm->ip[ipp].av[k]=val;
	    ipp++;
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
  The function computes initial values in integration points from initial nodal values.
  (function interpolates the initial nodal values to integration points)

  @param eid - element id
 
  @return The function does not return anything but stores computed values at int. points of the given element.

  Created by TKo, 4.7.2018
*/
void quadlineart::initintpointval (long eid)
{
  long i,j,k,ii,jj,ipp, ndofn, cndofn;
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
	
	//  interpolation to integration points of conductivity matrix
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
	
	//  interpolation to integration points of capacity matrix
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
   function computes gradients in integration points from nodal values
   
   eid - element id
   
   JK, 25.9.2001   
*/
void quadlineart::intpointgrad (long eid)
{
  long i,j,ii,jj,k,ipp;
  double xi,eta,jac;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),r(ASTCKVEC(ndofe)),t(ASTCKVEC(nne));
  vector gp,w,grad(ASTCKVEC(ncomp));
  matrix gm(ASTCKMAT(ncomp,nne));
  
  //  nodal coordinates
  Tt->give_node_coord2d (x,y,eid);
  //  nodal values on element
  elemvalues(eid, r);
  
  for (k=0;k<Tp->ntm;k++){
    
    //  nodal values of a single variable
    for (i=0;i<dofe[k][k];i++){
      t[i]=r[ordering[k][i]-1];
    }
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
	
	//  computation of gradients in integration points of conductivity matrix
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
	
	//  computation of gradients in integration points of capacity matrix
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
   The function computes values of array other in the integration points
   
   @param eid - element id
   
   TKo, 6.2022
*/
void quadlineart::intpointother (long eid)
{
  long i, j, ipp, ncompo, tnip;
  double otherval;
  vector ipval(ASTCKVEC(ntm));
  
  ipp = Tt->elements[eid].ipp[0][0];
  tnip = Tt->give_tnip(eid);

  for (i=0; i < tnip; i++){
    // number of components in the other array
    ncompo = Tm->ip[ipp].ncompother;

    // assamble vector of values of primary unknowns in the given integration point
    for (j=0; j < ntm; j++)
      ipval(j) = Tm->ip[ipp].av[j];

    for (j=0; j<ncompo; j++){
      // compute j-th other value at the given integration point ipp
      Tm->ip[ipp].other[j] = otherval = Tm->givecompother(j, ipp, ipval.a);      
    }
    ipp++;
  }
}



/**
   function assembles approximation or test functions into a %matrix

   @param n - %matrix of test functions
   @param xi,eta - natural coordinates
   @param v - %vector of velocity
   @param eid - element id
   
   JK, 2. 10. 2016
*/
void quadlineart::btf_matrix (matrix &n,double xi,double eta,matrix &d,long eid)
{
  if (Tp->advect==0)
    bf_matrix (n,xi,eta);
  if (Tp->advect==1)
    tf_matrix (n,xi,eta,d,eid);
}

/**
   function assembles %matrix of base functions

   @param n - %matrix of base functions
   @param xi,eta - natural coordinates

   JK, 25.9.2001
*/
void quadlineart::bf_matrix (matrix &n,double xi,double eta)
{
  nullm (n);
  bf_lin_4_2d (n.a,xi,eta);
}

/**
   function assembles %matrix of test functions

   @param n - %matrix of test functions
   @param xi,eta - natural coordinates
   @param v - %vector of velocity
   @param eid - element id
   
   JK, 23. 9. 2016
*/
void quadlineart::tf_matrix (matrix &w,double xi,double eta,matrix &d,long eid)
{
  long i;
  double alpha,jac,h,nor,pe;
  vector dx(ASTCKVEC(nne)), dy(ASTCKVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), v(ASTCKVEC(2));
  matrix n(ASTCKMAT(1, dofe[0][0]));
  
  nullm (n);
  
  //  node coordinates
  Tt->give_node_coord2d (x,y,eid);

  v[0]=d[0][0];
  v[1]=d[0][1];
	
  //  approximation functions
  bf_lin_4_2d (n.a,xi,eta);

  //  derivatives of the approximation functions with respect to the natural coordinates
  dx_bf_lin_4_2d (dx.a,eta);
  dy_bf_lin_4_2d (dy.a,xi);
  //  derivatives of approximation functions with respect to the x,y coordinates
  derivatives_2d (dx,dy,jac,x,y,xi,eta);
  
  //  norm of the velocity vector
  scprd (v,v,nor);
  nor = sqrt(nor);
  
  if (nor<Tp->zero){
    print_err("nonpositive norm of the advection velocity vector",__FILE__, __LINE__, __func__);
  }
  
  //  characteristic size of the element
  h = Gtt->give_characteristic_size (eid);
  //h=1.2500000000e-03;

  pe=nor*h/2.0/9.52e-12;
  alpha=(exp(pe)+exp(0.0-pe))/(exp(pe)-exp(0.0-pe))-1.0/pe;
  alpha =1.0;

  for (i=0;i<nne;i++){
    w[0][i]=n[0][i]+(dx[i]*v[0]+dy[i]*v[1])*alpha*h/nor/2.0;
  }
}

/**
   function assembles gradient of %matrix of test or base functions

   @param gm - gradient %matrix
   @param x,y - array containing node coordinates
   @param xi,eta - natural coordinates
   @param jac - jacobian
   @param eid - element id
   @param ipp - id of integration point
   @param ri,ci - row and column indices
   
   JK, 8. 10. 2016
*/
void quadlineart::grad_btf_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac,long eid,long ipp,long ri,long ci)
{
  if (Tp->advect==0){
    //  Galerkin-Bubnov method
    //  test functions are equal to the approximation functions
    grad_matrix (gm,x,y,xi,eta,jac);
  }
  if (Tp->advect==1){
    //  Galerkin-Petrov method
    //  test functions are different from the approximation functions
    matrix d(ASTCKMAT(1,ncomp));
    //  velocity vector
    Tm->matcond2 (d,ipp,ri,ci);
    //  gradients of test functions
    grad_tf_matrix (gm,x,y,xi,eta,d,eid);
  }
}

/**
   function assembles gradient of %matrix of base functions

   @param gm - gradient %matrix
   @param x,y - array containing node coordinates
   @param xi,eta - natural coordinates
   @param jac - Jacobian

   JK, 25.9.2001
*/
void quadlineart::grad_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac)
{
  long i;
  vector dx(ASTCKVEC(nne)), dy(ASTCKVEC(nne));

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
   function assembles gradient of %matrix of test functions

   @param gm - gradient %matrix
   @param x,y - array containing node coordinates
   @param xi,eta - natural coordinates
   @param d - %column matrix of velocity
   @param eid - element id

   JK, 8. 10. 2016
*/
void quadlineart::grad_tf_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,matrix &d,long eid)
{
  long i;
  double h,nor,jac,alpha=1.0;
  vector dx(ASTCKVEC(nne)),dy(ASTCKVEC(nne));
  vector dxdx(ASTCKVEC(nne)),dxdy(ASTCKVEC(nne)),dydy(ASTCKVEC(nne));
  vector v(ASTCKVEC(nne));
  
  dx_bf_lin_4_2d (dx.a,eta);
  dy_bf_lin_4_2d (dy.a,xi);
  
  //  first derivatives
  derivatives_2d (dx,dy,jac,x,y,xi,eta);
  
  //  second derivatives
  sec_derivatives_2d (quadrilateral_bilinear,dxdx,dxdy,dydy,x,y,xi,eta);
  
  //  velocity vector
  v[0]=d[0][0];
  v[1]=d[0][1];
  
  //  norm of the velocity vector
  scprd (v,v,nor);
  nor = sqrt(nor);
  
  if (nor<Tp->zero){
    print_err("nonpositive norm of the advection velocity vector",__FILE__, __LINE__, __func__);
  }
  
  //  characteristic size of the element
  h = Gtt->give_characteristic_size (eid);

  nullm (gm);
  
  for (i=0;i<nne;i++){
    gm[0][i]=dx[i] + alpha*h/nor*(dxdx[i]*v[0]+dxdy[i]*v[1]);
    gm[1][i]=dy[i] + alpha*h/nor*(dxdy[i]*v[0]+dydy[i]*v[1]);
  }

}



/**
   function computes conductivity %matrix of 2D problems for one transported matter
   finite element with bilinear approximation functions

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param km - conductivity %matrix

   JK, 25.9.2001
*/
void quadlineart::conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,j,ii;
  double thick,xi,eta,ww1,ww2,jac;
  ivector nodes(nne);
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordkm[ri][ci])),gp(ASTCKVEC(intordkm[ri][ci])),t(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(ncomp,dofe[ri][ci])),gmt(ASTCKMAT(ncomp,dofe[ri][ci])),d;
  matrix n(ASTCKMAT(1,dofe[ri][ci]));
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  thickness
  Tc->give_thickness (eid,nodes,t);
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
      
      //  matrix of gradients of approximation functions
      grad_matrix (gm,x,y,xi,eta,jac);
      //  matrix of gradients of test functions
      //grad_btf_matrix (gmt,x,y,xi,eta,jac,eid,ii,ri,ci);
      
      //  matrix of conductivity of the material
      reallocm(RSTCKMAT(ncomp,ncomp,d));
      Tm->matcond (d,ii,ri,ci);
      
      //  thickness
      thick = approx (xi,eta,t);
      
      if (jac<0.0){
	print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	//jac=fabs(jac);
	abort();
      }
      
      jac*=thick*ww1*ww2;
      
      //  contribution to the conductivity matrix of the element
      //bdbjac(km, gmt, d, gm, jac);	
      bdbjac(km, gm, d, gm, jac);	
      
      //  convective terms
      reallocm(RSTCKMAT(1,ncomp,d));
      Tm->matcond2(d,ii,ri,ci);
      btf_matrix (n,xi,eta,d,eid);
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
   function computes L %matrix

   L = \int_{\Omega} D B {\rm d} \Omega

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param lm - L %matrix

   TKr, 05/04/2011
*/
void quadlineart::l_matrix (long /*lcid*/,long eid,long ri,long ci,matrix &lm)
{
  long i,j,ii;
  double thick,xi,eta,ww1,ww2,jac;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), w(ASTCKVEC(intordkm[ri][ci]));
  vector gp(ASTCKVEC(intordkm[ri][ci])), t(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(ncomp,dofe[ri][ci])),d;
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  thickness
  Tc->give_thickness (eid,nodes,t);
  //  node coordinates
  Tt->give_node_coord2d (x,y,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  
  nullm (lm);
  
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
      Tm->matcond (d,ii,ri,ci);
      
      //  thickness
      thick = approx (xi,eta,t);
      
      if (jac<0.0){
	print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	//jac=fabs(jac);
	abort();
      }
      
      jac*=thick*ww1*ww2;
      
      //  contribution to the L matrix of the element
      mxm(d,gm,lm);
      cmulm(jac,lm);
      
      ii++;  
    }
  }
}


/**
   function computes L^T (L transposed) %matrix

   L^T = \int_{\Omega} B^T D {\rm d} \Omega

   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param lm - L^T %matrix

   TKr, 05/04/2011
*/
void quadlineart::l_t_matrix (long /*lcid*/,long eid,long ri,long ci,matrix &lm)
{
  long i,j,ii;
  double thick,xi,eta,ww1,ww2,jac;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), w(ASTCKVEC(intordkm[ri][ci]));
  vector gp(ASTCKVEC(intordkm[ri][ci])), t(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(ncomp,dofe[ri][ci])), d;
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  thickness
  Tc->give_thickness (eid,nodes,t);
  //  node coordinates
  Tt->give_node_coord2d (x,y,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);

  nullm (lm);
  
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
      Tm->matcond (d,ii,ri,ci);
      
      //  thickness
      thick = approx (xi,eta,t);
      
      if (jac<0.0){
	print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	//jac=fabs(jac);
	abort();
      }
      
      jac*=thick*ww1*ww2;
      
      //  contribution to the L matrix of the element
      mtxm(gm,d,lm);
      cmulm(jac,lm);
      
      ii++;  
    }
  }
}



/**
   function computes capacity %matrix of 2D problems for one transported matter
   finite element with bilinear approximation functions

   @param eid - element id
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param cm - capacity %matrix

   JK, 4.10.2001
*/
void quadlineart::capacity_matrix (long eid,long ri,long ci,matrix &cm)
{
  long i,j,ii;
  double jac,xi,eta,w1,w2,thick,rho,c;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordcm[ri][ci])), gp(ASTCKVEC(intordcm[ri][ci])), t(ASTCKVEC(nne));
  vector dens(ASTCKVEC(nne));
  matrix n(ASTCKMAT(1,dofe[ri][ci]));

  Tt->give_elemnodes (eid,nodes);
  Tc->give_thickness (eid,nodes,t);
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
      
      if (Tp->react_capac==capacity)
	c=Tm->capcoeff (ii,ri,ci);
      if (Tp->react_capac==reaction)
	c=Tm->reactcoeff (ii,ri,ci);

      thick = approx (xi,eta,t);
      rho = approx (xi,eta,dens);
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
void quadlineart::quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci)
{
  long i,j;
  double jac,xi,eta,w1,w2,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordcm[ri][ci])), gp(ASTCKVEC(intordcm[ri][ci]));
  vector t(ASTCKVEC(nne)), dens(ASTCKVEC(nne)), v(ASTCKVEC(dofe[ri][ci]));
  matrix n(ASTCKMAT(1,dofe[ri][ci])), nm(ASTCKMAT(dofe[ri][ci],dofe[ri][ci]));
  
  Tt->give_elemnodes (eid,nodes);
  
  Tc->give_thickness (eid,nodes,t);
  
  Tt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,intordcm[ri][ci]);
  
  nullm (nm);
  
  for (i=0;i<intordcm[ri][ci];i++){
    xi=gp[i];  w1=w[i];
    for (j=0;j<intordcm[ri][ci];j++){
      eta=gp[j];  w2=w[j];
      jac_2d (jac,x,y,xi,eta);
      //  matrix of shape functions
      bf_matrix (n,xi,eta);
      //  thickness
      thick = approx (xi,eta,t);

      if (jac<0.0){
	print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	abort();
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
void quadlineart::internal_fluxes (long lcid,long eid,vector &ifl)
{
  long i,j,ipp;
  double xi,eta,jac,thick;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), w, gp;
  vector t(ASTCKVEC(nne)), fl(ASTCKVEC(ncomp)), contr(ASTCKVEC(dofe[lcid][lcid]));
  matrix gm(ASTCKMAT(ncomp,dofe[lcid][lcid]));
  

  Tt->give_elemnodes (eid,nodes);
  Tc->give_thickness (eid,nodes,t);
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
      
      //  thickness
      thick = approx (xi,eta,t);
      //  computation of fluxes
      Tm->computenlfluxes (lcid,ipp);
      
      Tm->givefluxes (lcid,ipp,fl);
      //  matrix of gradients
      grad_matrix (gm,x,y,xi,eta,jac);
      
      mtxv (gm,fl,contr);
      
      if (jac<0.0){
	print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	jac=fabs(jac);
      }
      
      cmulv (jac*w[i]*w[j]*thick,contr);
      
      addv (contr,ifl,ifl);
      
      ipp++;
    }
  }
  
}

/**
   function computes advection %matrix of 2D problems for one transported matter
   finite element with bilinear approximation functions

   @param lcid - number of load case
   @param eid - element id
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param hm - advection %matrix

   JK, 2. 5. 2016
*/
void quadlineart::advection_matrix (long /*lcid*/,long eid,long ri,long ci,matrix &hm)
{
  long i,j,ii;
  double thick,xi,eta,ww1,ww2,jac;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordkm[ri][ci])), gp(ASTCKVEC(intordkm[ri][ci])), t(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(ncomp,dofe[ri][ci])), n(ASTCKMAT(1,dofe[ri][ci])), v(ASTCKMAT(1,ncomp));
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  thickness
  Tc->give_thickness (eid,nodes,t);
  //  node coordinates
  Tt->give_node_coord2d (x,y,eid);
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  
  nullm (hm);
  
  if (Tp->savemode==0)
    ii=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ii=Tt->elements[eid].ipp[0][0];
  
  for (i=0;i<intordkm[ri][ci];i++){
    xi=gp[i];  ww1=w[i];
    for (j=0;j<intordkm[ri][ci];j++){
      eta=gp[j];  ww2=w[j];
      
      //  matrix of approximation functions
      bf_matrix (n,xi,eta);
      //  matrix of gradients
      grad_matrix (gm,x,y,xi,eta,jac);
      
      //  vector of velocity
      Tm->matcond2(v,ii,ri,ci);
      
      //  thickness
      thick = approx (xi,eta,t);
      
      if (jac<0.0){
	print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	jac=fabs(jac);
      }
      
      jac*=thick*ww1*ww2;
      
      //  contribution to the advection matrix of the element
      bdbjac (hm,n,v,gm,jac);
      
      ii++;  
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
void quadlineart::res_conductivity_matrix (long eid,long /*lcid*/,matrix &km)
{
  long i,j;
  ivector rcn, ccn;
  matrix lkm;

  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      reallocv(RSTCKIVEC(dofe[i][j], rcn));
      reallocv(RSTCKIVEC(dofe[i][j], ccn));
      reallocm(RSTCKMAT(dofe[i][j], dofe[i][j], lkm));

      conductivity_matrix (i,eid,i,j,lkm);
      
      codnum (rcn.a,i);
      codnum (ccn.a,j);
      mat_localize (km,lkm,rcn.a,ccn.a);
    }
  }
}



/**
   function assembles resulting element L %matrix

   L^T = \int_{\Omega} B^T D {\rm d} \Omega

   @param eid - element id
   @param lcid - load case id
   @param lm - resulting L %matrix of one element

   TKr, 05/04/2010
*/
void quadlineart::res_l_matrix (long eid,long /*lcid*/,matrix &lm)
{
  long i, j;
  ivector rcn, ccn;
  matrix lkm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      reallocv(RSTCKIVEC(2, rcn));
      reallocv(RSTCKIVEC(dofe[i][j], ccn));
      if (i==0){
	rcn[0]=1;
	rcn[1]=2;
      }
      if (i==1){
	rcn[0]=3;
	rcn[1]=4;
      }
      reallocm (RSTCKMAT(ncomp,dofe[i][j],lkm));
      l_matrix (i,eid,i,j,lkm);
      codnum (ccn.a,j);
      
      mat_localize (lm,lkm,rcn.a,ccn.a);
    }
  }
}



/**
   function assembles resulting element L^T (L transposed) %matrix

   L^T = \int_{\Omega} B^T D {\rm d} \Omega

   @param eid - element id
   @param lcid - load case id
   @param lm - resulting L^T %matrix of one element

   TKr, 05/04/2010
*/
void quadlineart::res_l_t_matrix (long eid,long /*lcid*/,matrix &lm)
{
  long i, j;
  ivector rcn, ccn;
  matrix lkm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      reallocv(RSTCKIVEC(dofe[i][j], rcn));
      reallocv(RSTCKIVEC(2, ccn));

      if (j==0){
	ccn[0]=1;
	ccn[1]=2;
      }
      if (j==1){
	ccn[0]=3;
	ccn[1]=4;
      }
      reallocm (RSTCKMAT(dofe[i][j],ncomp,lkm));
      l_t_matrix (i,eid,i,j,lkm);
      codnum (rcn.a,i);
      
      mat_localize (lm,lkm,rcn.a,ccn.a);
    }
  }
}




/**
   function assembles average D %matrix

   @param eid - element id
   @param lm - resulting D %matrix of one element

   TKr, 05/04/2010
*/
void quadlineart::averd_matrix (long eid,matrix &lm)
{
  long i, j, ii;
  ivector rcn(ASTCKIVEC(ncomp)), ccn(ASTCKIVEC(ncomp));
  matrix d(ASTCKMAT(ncomp,ncomp));
  matrix dd(ASTCKMAT(ncomp,ncomp));
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      
      if (Tp->savemode==0)
	ii=Tt->elements[eid].ipp[i][j];
      if (Tp->savemode==1)
	ii=Tt->elements[eid].ipp[0][0];
      
      ///////////////////////////////////////
      //this is new: it should be tested
      //averaging over int. points
      /* fillm (0.0,d);//null matrix
	 nii = 0;
	 for (k=0;k<intordkm[i][j];k++){
	 for (l=0;l<intordkm[i][j];l++){
	 Tm->matcond (dd,ii,i,j);
	 d[i][j] += dd[i][j];
	 nii++;
	 ii++; 
	 }
	 }
	 d[i][j] = d[i][j]/nii;//averaging on int. points
      */
      ///////////////////////////////////////
      //this is old:
      //  matrix of conductivity of material
      Tm->matcond (d,ii,i,j);
      ///////////////////////////////////////
      
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

      mat_localize (lm,d,rcn.a,ccn.a);
    }
  }
}





/**
   function assembles average C %matrix

   @param eid - element id
   @param lm - resulting C %matrix of one element

   TKr, 05/04/2011
*/
void quadlineart::averc_matrix (long eid,matrix &lm)
{
  long i,j,ii;
  double c,rho;
  ivector nodes(ASTCKIVEC(nne));
  vector dens(ASTCKVEC(nne));

  Tt->give_elemnodes (eid,nodes);
  Tc->give_density (eid,nodes,dens);
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      
      if (Tp->savemode==0)
	ii=Tt->elements[eid].ipp[i][j]+intordkm[i][j]*intordkm[i][j];
      if (Tp->savemode==1)
	ii=Tt->elements[eid].ipp[0][0]+intordkm[0][0]*intordkm[0][0];
      
            
      ///////////////////////////////////////
      //this is new: it should be tested
      //averaging over int. points
      /* c = 0.0;//null matrix
	 nii = 0;
	 for (k=0;k<intordcm[i][j];k++){
	 for (l=0;l<intordcm[i][j];l++){
	 cc = Tm->capcoeff (ii,i,j);
	 c += cc;
	 nii++;
	 ii++; 
	 }
	 }
	 c = c/nii;//averaging on int. points
      */
      ///////////////////////////////////////
      //this is old:
      //  coefficient of capacity of material
      c = Tm->capcoeff (ii,i,j);
      ///////////////////////////////////////


      //temporarily
      rho = approx (0.0,0.0,dens);
      
      lm[i][j] = rho*c*Tm->ip[ii].av[j];
      //lm[i][j] = rho*c;//here, it should be disscussed
    }
  }
}


/**
   function assembles area of one element

   @param eid - element id

   TKr, 05/04/2011
*/
double quadlineart::elem_area (long eid)
{
  double area,jac;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(ncomp,dofe[0][0]));

  Tt->give_node_coord2d (x,y,eid);

  //  matrix of gradients - jacobian
  grad_matrix (gm,x,y,0.0,1.0,jac);

  //  area is equal to four times jac of the element
  area = jac*4.0;

  return area;
}




/**
   function assembles resulting element capacity %matrix
   
   @param eid - element id
   @param cm - resulting capacity %matrix of one element

   JK, 6.1.2002
*/
void quadlineart::res_capacity_matrix (long eid, matrix &cm)
{
  long i,j;
  ivector rcn, ccn;
  matrix lcm;

  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      reallocv(RSTCKIVEC(dofe[i][j], rcn));
      reallocv(RSTCKIVEC(dofe[i][j], ccn));
      reallocm(RSTCKMAT(dofe[i][j],dofe[i][j],lcm));
      //  capacity matrix for one matter
      capacity_matrix (eid,i,j,lcm);
      
      // diagonalization of capacity matrix on one element
      if (Tp->diagcap == 1){
	diagonalization (lcm);
      }
      
      codnum (rcn.a,i);
      codnum (ccn.a,j);
      //  localization to element capacity matrix
      mat_localize (cm,lcm,rcn.a,ccn.a);
    }
  }
}

/**
   function assembles resulting element reaction %matrix

   the reaction %matrix is generated by reaction term in
   diffusion-reaction equation
   c du/dt = k d^2u/dt^2 + r u
   where c is the capacity, k is the permeability,
   r is reactivity (it is a material parameter) and
   u is the unknown function

   @param eid - element id
   @param cm - resulting capacity %matrix of one element

   JK, 12. 6. 2019
*/
void quadlineart::res_reaction_matrix (long eid, matrix &rm)
{
  long i, j, k, l;
  ivector rcn, ccn;
  double s;
  matrix lrm;

  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      reallocv(RSTCKIVEC(dofe[i][j], rcn));
      reallocv(RSTCKIVEC(dofe[i][j], ccn));
      reallocm(RSTCKMAT(dofe[i][j],dofe[i][j],lrm));
      //  reaction matrix for one matter
      capacity_matrix (eid,i,j,lrm);
      
      // diagonalization of reaction matrix on one element
      if (Tp->diagreact == 1){
	for (k=0;k<dofe[i][j];k++){
	  s=0.0;
	  for (l=0;l<dofe[i][j];l++){
	    s+=lrm[k][l];
	    lrm[k][l]=0.0;
	  }
	  lrm[k][k]=s;
	}
      }
      
      codnum (rcn.a,i);
      codnum (ccn.a,j);
      //  localization to element reaction matrix
      mat_localize (rm,lrm,rcn.a,ccn.a);
    }
  }
}


/**
   function assembles resulting element convection %vector

   @param f - resulting convection %vector of one element
   @param lcid - load case id
   @param eid - element id
   @param leid - id of element with boundary condition
   
   JK, 6.1.2002
*/
void quadlineart::res_convection_vector (vector &f, long lcid, long eid, long leid)
{
  ivector cn;
  vector lf;
  
  //  transi[lcid]==2 - element contains boundary with prescribed flux
  //  transi[lcid]==4 - element contains boundary with prescribed transmission and prescribed flux
  if ((Tt->elements[eid].transi[lcid]==2)||(Tt->elements[eid].transi[lcid]==4)){
    //  array for code numbers
    reallocv(RSTCKIVEC(dofe[lcid][lcid], cn));;
    //  array for nodal values of one matter
    reallocv (RSTCKVEC(dofe[lcid][lcid],lf));
    
    convection_vector (lf,lcid,eid,leid,lcid);
    
    //  code numbers
    codnum (cn.a,lcid);
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
   @param leid - id of element with boundary condition
   
   JK, 6.1.2002, 18. 4. 2014
*/
void quadlineart::res_transmission_vector (vector &f, long lcid, long eid, long leid)
{
  long i;
  vector lf;
  ivector cn(ASTCKIVEC(dofe[lcid][lcid]));
  //  code numbers
  codnum (cn.a,lcid);
  reallocv (RSTCKVEC(dofe[lcid][lcid],lf));
  
  //  transi[i]==3 - element contains boundary with prescribed transmission
  //  transi[i]==4 - element contains boundary with prescribed transmission and prescribed flux
  if ((Tt->elements[eid].transi[lcid]==3)||(Tt->elements[eid].transi[lcid]==4)){
    
    for (i=0;i<ntm;i++){    
      nullv (lf);
      
      transmission_vector (lf,lcid,eid,leid,i);
      
      locglob (f.a, lf.a, cn.a, dofe[lcid][lcid]);
    }
  }
}



/**
   function assembles resulting element source %vector

   @param sv - resulting source %vector of one element
   @param nodval - array of nodal values
   @param lcid - load case id
   @param eid - element id
   
   JK, 6.1.2002
*/
void quadlineart::res_quantity_source_vector (vector &sv,vector &nodval,long lcid,long eid)
{
  vector lsv(ASTCKVEC(dofe[lcid][lcid]));
  ivector cn(ASTCKIVEC(dofe[lcid][lcid]));

  quantity_source_vector (lsv,nodval,eid,lcid,lcid);
  codnum (cn.a,lcid);
  locglob (sv.a,lsv.a,cn.a,dofe[lcid][lcid]);
}



/**
   function assembles resulting element internal fluxes %vector

   @param eid - element id
   @param elemif - resulting internal fluxes %vector of one element
   
   JK, 6.1.2002
*/
void quadlineart::res_internal_fluxes (long eid,vector &elemif)
{
  long i;
  ivector cn;
  vector lif, tdnv;
  matrix cm;

  for (i=0;i<ntm;i++){
    reallocv (RSTCKIVEC(dofe[i][i], cn));
    reallocv (RSTCKVEC(dofe[i][i],lif));
    internal_fluxes (i,eid,lif);
    codnum (cn.a,i);
    locglob (elemif.a,lif.a,cn.a,dofe[i][i]);
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
   function assembles resulting element advection %matrix

   @param eid - element id
   @param lcid - load case id
   @param hm - resulting advection %matrix of one element

   JK, 2. 5. 2016
*/
void quadlineart::res_advection_matrix (long eid,long /*lcid*/,matrix &hm)
{
  long i, j;
  ivector rcn, ccn;
  matrix lhm;

  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      reallocv(RSTCKIVEC(dofe[i][j], rcn));
      reallocv(RSTCKIVEC(dofe[i][j], ccn));
      reallocm(RSTCKMAT(dofe[i][j],dofe[i][j],lhm));

      advection_matrix (i,eid,i,j,lhm);
      
      codnum (rcn.a,i);
      codnum (ccn.a,j);
      mat_localize (hm,lhm,rcn.a,ccn.a);
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
double quadlineart::total_integral(long eid,vector &nodval)
{
  long i,j;
  double thick,xi,eta,ww1,ww2,jac,value,f;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),w(ASTCKVEC(2)),gp(ASTCKVEC(2)),t(ASTCKVEC(nne));

  Tt->give_elemnodes (eid,nodes);
  Tc->give_thickness (eid,nodes,t);

  Tt->give_node_coord2d (x,y,eid);
  gauss_points (gp.a,w.a,2);

  f = 0.0;

  for (i=0;i<2;i++){
    xi=gp[i];  ww1=w[i];
    for (j=0;j<2;j++){
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
   @param leid - element id
   
   TKr, 28.2.2004
*/
void quadlineart::res_boundary_flux (vector &f,long lcid,long eid,long leid)
{
  long i;
  ivector cn(ASTCKIVEC(dofe[lcid][lcid]));
  vector lf;

  codnum (cn.a,lcid);
  
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
   function computes contributions to the right-hand side - volume integral
      
   \int_{Omega} B^T D dOmega
   
   @param vrhs - volume right-hand side %vector of one matter
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices of the block (ri must be equal to ci)
   
   TKr, 10/12/2013
*/
void quadlineart::volume_rhs_vector (long /*lcid*/,long eid,long ri,long ci,vector &vrhs)
{
  long i,j,ii;
  double thick,xi,eta,ww1,ww2,jac;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), w(ASTCKVEC(intordkm[ri][ci]));
  vector gp(ASTCKVEC(intordkm[ri][ci])), t(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(ncomp, dofe[ri][ci])), d;
  matrix km(ASTCKMAT(dofe[ri][ci], 2));
  
  //  nodes on required elements
  Tt->give_elemnodes (eid,nodes);
  //  node coordinates
  Tt->give_node_coord2d (x,y,eid);
  //  thickness
  Tc->give_thickness (eid,nodes,t);
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
      
      //  thickness
      thick = approx (xi,eta,t);
      
      jac*=ww1*ww2*thick;
      
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
void quadlineart::volume_rhs_vector2 (long /*lcid*/,long eid,long ri,long ci,vector &vrhs)
{
  long i,j,k,ipp;
  double jac,xi,eta,w1,w2,thick,c;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordcm[ri][ci])), gp(ASTCKVEC(intordcm[ri][ci])), t(ASTCKVEC(nne));
  matrix n(ASTCKMAT(1,dofe[ri][ci]));
  vector nn(ASTCKVEC(dofe[ri][ci]));
  
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

      thick = approx (xi,eta,t);

      jac*=w1*w2*thick*c;
   
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

   TKr, 10/12/2013
*/
void quadlineart::res_volume_rhs_vector (vector &f,long eid,long /*lcid*/)
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

   TKr, 16/05/2018
*/
void quadlineart::res_volume_rhs_vector2 (vector &f,long eid,long /*lcid*/)
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
   function computes nodal fluxes caused by prescribed fluxes
   (Neumann boundary condition)
   
   @param v - array of nodal fluxes
   @param lcid - load case id (corresponds to the row index)
   @param eid - element id
   @param leid - loaded element id
   @param cid - component id (corresponds to the column index)
   
   JK, 19.8.2004
*/
void quadlineart::convection_vector (vector &v,long lcid,long eid,long leid,long cid)
{
  long i,ipp;
  bocontypet *bc;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordkm[lcid][cid])), gp(ASTCKVEC(intordkm[lcid][cid])), t(ASTCKVEC(nne));
  vector list(ASTCKVEC(nned*ned)), trc(ASTCKVEC(nned*ned)), nodval(ASTCKVEC(nne));
  vector av(ASTCKVEC(nne)),coef(ASTCKVEC(nne));
  matrix km(ASTCKMAT(nne,nne));
  
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
  Tb->lc[lcid].elemload[leid].give_nodval (lcid,list);
  
  //  auxiliary coefficients, necessary for function edge_integral
  for (i=0;i<nned*ned;i++){
    trc[i]=1.0;
  }
  
  //  integration point number
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[lcid][cid];
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
      edge_integral (i,x,y,intordkm[lcid][cid],gp,w,t,coef,km);
      
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
   @param eid - element id
   @param ri,ci - row and column indices
   @param km - part of the conductivity %matrix
   
   JK, 19.8.2004
*/
void quadlineart::transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,j,ipp,leid;
  bocontypet *bc;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordkm[ri][ci])), gp(ASTCKVEC(intordkm[ri][ci]));
  vector trc(ASTCKVEC(nned*ned)), coeff(ASTCKVEC(nne)), t(ASTCKVEC(nne));
  
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

  //  correspondence between element ordering and boundary element ordering
  for (i=0;i<Tb->lc[lcid].neb;i++){
    if (Tb->lc[lcid].elemload[i].eid==eid){
      leid=i;
      
      Tb->lc[lcid].elemload[leid].give_bc (bc);
      //  transmission coefficients
      Tb->lc[lcid].elemload[leid].give_trc (Tp->time,Tb->lc[lcid].nodval,lcid,ci,trc);
      
      //  integration point number
      if (Tp->savemode==0)
	ipp=Tt->elements[eid].ipp[ri][ci];
      if (Tp->savemode==1)
	ipp=Tt->elements[eid].ipp[0][0];
      
      //  loop over edges
      for (j=0;j<ned;j++){
	
	if (bc[j]==5 || bc[j]>10){
	  //  transformation of coefficients
	  transf_coeff (j,coeff,trc,eid,ri,ci,ipp,bc,0);
	  //  matrix obtained from integration over edges
	  edge_integral (j,x,y,intordkm[ri][ci],gp,w,t,coeff,km);
	}
      }
    }
  }
  
  delete [] bc;
}

/**
   function computes nodal values on element caused by
   prescribed transmission (Newton boundary condition)
   
   function computes \int \kappa T_{ext}
   
   @param v - %vector of nodal values
   @param lcid - load case id (corresponds to the row index)
   @param eid - element id
   @param leid - loaded element id
   @param cid - component id (corresponds to the column index)
   
   JK, 18. 4. 2014
*/
void quadlineart::transmission_vector (vector &v,long lcid,long eid,long leid,long cid)
{
  long i,ipp;
  bocontypet *bc;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordkm[lcid][cid])), gp(ASTCKVEC(intordkm[lcid][cid])), t(ASTCKVEC(nne));
  vector list(ASTCKVEC(nned*ned)), trc(ASTCKVEC(nned*ned)), trr(ASTCKVEC(nned*ned));
  vector nodval(ASTCKVEC(nne)), av(ASTCKVEC(nne)), coef(ASTCKVEC(nne));
  matrix km(ASTCKMAT(nne,nne));
  
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
  
  
  //  loop over the nmber of edges
  for (i=0;i<ned;i++){
    
    if (bc[i]==5 || bc[i]>10){
      //  transformation of prescribed nodal values
      transf_val (i,nodval,list,trc,trr,eid,lcid,cid,ipp,bc);
      //  determination of transmission coefficients on edge
      transf_coeff (i,coef,trc,eid,lcid,cid,ipp,bc,1);

      /* if (bc[i]==50){
	 fprintf(Outt,"nodval[0] = %e\n",nodval[0]);
	 fprintf(Outt,"nodval[1] = %e\n",nodval[1]);
	 fprintf(Outt,"nodval[2] = %e\n",nodval[2]);
	 fprintf(Outt,"nodval[3] = %e\n",nodval[3]);
	 fprintf(Outt,"coef[0] = %e\n",coef[0]);
	 fprintf(Outt,"coef[1] = %e\n",coef[1]);
	 fprintf(Outt,"coef[2] = %e\n",coef[2]);
	 fprintf(Outt,"coef[3] = %e\n",coef[3]);
	 fprintf(Outt,"\n");
	 fflush(Outt);
	 }
      */

      nullm (km);
      //  matrix obtained from integration over edges
      edge_integral (i,x,y,intordkm[lcid][cid],gp,w,t,coef,km);
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

   @param tmv - %vector of boundary fluxes
   @param lcid - load case id
   @param eid - element id
   @param leid - element id
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   
   TKr, 28.2.2004
*/
void quadlineart::boundary_flux (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  ivector nodes(ASTCKIVEC(nne));
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne));
  vector w(ASTCKVEC(intordkm[ri][ci])), gp(ASTCKVEC(intordkm[ri][ci])), t(ASTCKVEC(nne));
  vector list(ASTCKVEC(nned*ned)), trc(ASTCKVEC(nned*ned)), trr(ASTCKVEC(nned*ned));
  vector nodval(ASTCKVEC(nne)), av(ASTCKVEC(nne)), coef(ASTCKVEC(nne));
  matrix km(ASTCKMAT(nne,nne));
  
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
void quadlineart::edge_integral (long edg,vector &x,vector &y,long intord,vector &gp,vector &w,
				 vector &t,vector &coef,matrix &km)
{
  long i;
  double xi,eta,jac,ipval,thick;
  matrix n(ASTCKMAT(1,nne));
  
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
void quadlineart::transf_flux (long edg,vector &coeff,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double tr,new_flux;
  ivector nodes(ASTCKIVEC(nne)), edgenod(ASTCKIVEC(nned));

  //  zeroing of array of coefficients
  nullv (coeff);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required edge
  linquadrilat_edgnod (edgenod.a,edg);

  //  actual position in the array list
  j=edg*nned;

  for (i=0;i<nned;i++){
    tr = 0.0;
    if (bc[edg]==90){
      tr = trr[j+i]/trc[j+i];
    }

    //  node number
    k=edgenod[i];
    
    Tm->transmission_flux(new_flux,list[j+i],tr,ri,ci,nodes[k],bc[edg],ipp);
    
    coeff[k]=new_flux;
  }
  
}


/**
   function
   
   @param edg - edge id
   @param coeff - array of 
   @param list - array of nodal values prescribed on boundary
   @param eid - element id
   @param ri,ci - row and column indices
   @param ipp - integration point id
   @param bc - array describing type of boundary condition
   @param flag - coefficient is computing for what 0=matrix,1=loading vector
*/
void quadlineart::transf_coeff (long edg,vector &coeff,vector &list,long eid,long ri,long ci,long ipp,bocontypet *bc,int flag)
{
  long i,j,k;
  double new_coeff;
  ivector nodes(ASTCKIVEC(nne)), edgenod(ASTCKIVEC(nned));

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
    
    Tm->transmission_transcoeff(new_coeff,list[j+i],ri,ci,nodes[k],bc[edg],ipp,flag);
    
    coeff[k]=new_coeff;
  }
  
}


/**
   function prepares nodal values of transported quantity, transmission
   coefficents and transmission/radiation coefficients
   
   @param edg - number of required edge (edge id)
   @param nodval - array of transformed nodal values
   @param list - array of nodal values defined on all edges
   @param trc - array of transmission coefficients
   @param trr - array of transmission/radiation coefficients
   @param ri,ci - row and column indices
   @param ipp - integration point number
   @param bc - array defining boundary conditions
   
   JK, 19.8.2004
*/
void quadlineart::transf_val (long edg,vector &nodval,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double tr,new_nodval;
  ivector nodes(ASTCKIVEC(nne)), edgenod(ASTCKIVEC(nned));
  
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
    if (bc[edg] == 90) {
      tr = trr[j+i]/trc[j+i];
    }
    else
     tr = 0.0;
    
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
void quadlineart::edgenodeval (long edg,vector &nodval,vector &list)
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
   function computes correct fluxes at integration points on element

   @param eid - element id
   
   TKr, 01/02/2010
*/
void quadlineart::intpointflux (long eid)
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
   
   TKr, 01/02/2010
*/
void quadlineart::nod_grads_ip (long eid)
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
   
   TKr, 01/02/2010
*/
void quadlineart::nod_fluxes_ip (long eid)
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
   
   Tomas Koudelka, 6.5.2014
*/
void quadlineart::nod_eqother_ip (long eid)
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

   @param eid - element id
   
   TKr, 05/05/2022 according to TKo nod_eqother_ip (long eid)
*/
void quadlineart::nod_other_ip (long eid)
{
  long i,j,k,ipp,ncompo;
  double area;
  ivector ipnum(ASTCKIVEC(nne)), nod(ASTCKIVEC(nne));
  vector other;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  nodip_planelq (ipp,intordkm[0][0],ipnum);

  
  //  node numbers of the element
  Tt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    ncompo = Tm->givencompother ();
    reallocv (RSTCKVEC(ncompo,other));
    Tm->giveother (ipnum[i],0,ncompo,other.a);
    
    //  storage of other components to the node
    j=nod[i];
    for (k=0; k<ncompo; k++){
      if (Tp->otheraver==0 || Tp->otheraver==1)
	Tt->nodes[j].storeother (k,other(k));
      if (Tp->otheraver==2){
	area = element_area (eid);
	Tt->nodes[j].storeother (k,area,other(k));
      }
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
void quadlineart::nod_others_comp (long /*lcid*/,long eid,long ri,long ci)
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
   function selects values and gradients from the global level
   
   @param eid - element id
   @param counter - actual position in the array buff
   @param buff - array containing selected components
   
   JK, 19.4.2011
*/
void quadlineart::higher_to_lower_level (long eid,long *counter,double *buff)
{
  long i,ipp,k,nip;
  double av,gr0,gr1;
  vector gr(ASTCKVEC(2));
  
  nip=Tt->give_tnip (eid);
  
  //  loop over the number of transported materials
  for (i=0;i<ntm;i++){
    ///////////////////////////////////////
    //this is new: it should be tested
    //averaging over int. points
    av = 0.0;
    gr0 = 0.0;
    gr1 = 0.0;
    //over all integration points of this element
    
    //  id of the first integration point
    ipp=Tt->elements[eid].ipp[0][0];
    for (k=0;k<nip;k++){
      av += Tm->ip[ipp].av[i];
      Tm->givegrad (i,ipp,gr);
      //fprintf(Outt,"ipp = %ld, av[%ld] = %e, gr[0] =  %e, gr[1] =  %e\n",ipp,i,Tm->ip[ipp].av[i],gr[0],gr[1]);
      gr0 += gr[0];
      gr1 += gr[1];
      ipp++;
    }

    av = av/nip;
    gr0 = gr0/nip;
    gr1 = gr1/nip;
    ///////////////////////////////////////
    
    //  values
    //old:
    //buff[counter[0]]=Tm->ip[ipp].av[i];
    //new:
    buff[counter[0]] = av;
    counter[0]++;
    
    //  components of gradient
    //old:
    //Tm->givegrad (i,ipp,gr);
    //buff[counter[0]]=gr[0];
    //new:
    buff[counter[0]] = gr0;
    counter[0]++;
    //old:
    //buff[counter[0]]=gr[1];
    //new:
    buff[counter[0]] = gr1;
    counter[0]++;

    //fprintf(Outt,"eid = %ld, av = %e, gr[0] =  %e, gr[1] =  %e\n",eid,av,gr0,gr1);
    //fflush(Outt);
  }
}


/**
   Function 1.
   1. computes "e2" - square of energy norm of error of solution on element
      e2 = integral{ (sig_star - sig_roof)T * D_inv * (sig_star - sig_roof) }
      sig_star = recovered stress, values in nodes are defined by "rsigfull" and over element are interpolated by base functions
      sig_roof = stress obtained by FEM
      D_inv = inverse stiffness matrix of material
   2. computes "u2" - square of energy norm of strain on element
      u2 = integral{ epsT * D * eps }
      eps = strain obtained by FEM
      D = stiffness matrix of material
   3. computes area of element and adds to "volume" (sum of areas of previous elements)
   4. computes "sizel" (characteristic size of element)
   
   @param eid - element id
   @param volume - sum of areas of previous elements
   @param e2 - empty(returned) value; 
   @param u2 - empty(returned) value; 
   @param sizel - empty(returned) value; 
   @param rsigfull - 1Darray of stress in nodes; dimmension is ncomp*gt->nn after this manner {(val[0]; ... ;val[gt->nn])[0] ; ... ; (val[0]; ... ;val[gt->nn])[ncomp]}
   
   created  3.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
double quadlineart :: compute_error (long eid, vector *rderfull, int mattid, double &e2, double &u2, double &sizel)
{
  long intord = 2;
  long i, j, ipp;
  double xi, eta, contr, jac, area;
  ivector nodes(ASTCKIVEC(nne)),cn;
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), t(ASTCKVEC(nne));
  vector gp(ASTCKVEC(intord)), w(ASTCKVEC(intord)), bf(ASTCKVEC(nne)), r;
  vector der(ASTCKVEC(ncomp)), der_fine(ASTCKVEC(ncomp)), der_err(ASTCKVEC(ncomp)), eps(ASTCKVEC(ncomp));
  matrix d(ASTCKMAT(ncomp,ncomp)), gm;
  
  Tt->give_elemnodes (eid, nodes);
  Tt->give_node_coord2d (x, y, eid);
  Tc->give_thickness (eid, nodes, t);
  
  // kdyz to nebude savemode, tak misto 0,0 tam bude ri ci podle media
  if (Tp->savemode==1)  ipp = Tt->elements[eid].ipp[0][0];
  else { print_err("save mode is not 1", __FILE__, __LINE__, __func__);  exit (1); }
  
  Tm->matcond (d, ipp, mattid, mattid);
  
  gauss_points (gp.a, w.a, intord);
  
  e2 = u2 = 0;
  for (i=0; i<intord; i++) {
    xi=gp[i];
    for (j=0; j<intord; j++) {
      eta=gp[j];
      
      Tm->givegrad (mattid, ipp, der);
      ipp++;
      jac_2d (jac, x, y, xi, eta);
      
      bf_lin_4_2d (bf.a, xi, eta);
      give_der_star (bf, rderfull, nodes, der_fine, Tt->nn);
      subv (der_fine, der, der_err);
      
      jac *= w[i]*w[j] * approx (xi,eta,t);
      
      vxmxv (der,     d, contr);    u2 += contr * jac;
      vxmxv (der_err, d, contr);    e2 += contr * jac;
    }
  }
  
  area = ( (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]) + (x[2]-x[0])*(y[3]-y[0])-(x[3]-x[0])*(y[2]-y[0]) ) / 2.0;
  sizel = sqrt(area);
  return area;
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
void quadlineart::transq_nodval (long eid,vector &nodval,nonmechquant nmq)
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
void quadlineart::transq_nodval_comp (long eid,vector &nodval, long ncne, long nq, nonmechquant *qt)
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
void quadlineart::transq_init_nodval (long eid,vector &nodval,nonmechquant nmq)
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
void quadlineart::transq_init_nodval_comp (long eid,vector &nodval, long ncne, long nq, nonmechquant *qt)
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



/**
   Function assembles global coordinates of integration point.
   
   @param eid - element id (input)
   @param ipp - integration point pointer (input)
   @param ri - row index of the block of integration points /according to transported media/ (input)
   @param ci - column index of the block of integration points /according to transported media/ (input)
   @param coord - array containing coordinates of integration point (output)
   
   @retval coord - function returns coordinates of integration point in the %vector coord.
   @retval 0 - on success
   @retval 1 - integration point with ipp was not found

   TKo, 12.2016
*/
long quadlineart::ipcoord (long eid, long ipp, long ri, long ci, vector &coord)
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
long quadlineart::ipncoord (long eid, long ipp, vector &ncoord)
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
void quadlineart::surface_flux (long lcid,long eid,long beid,double *fluxes)
{

  long i,ipp,fid;
  double q,qn,area;
  bocontypet *bc;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne));
  ivector nodes(nne);
  vector w, gp;
  vector flux(ASTCKVEC(ncomp)),fl(ASTCKVEC(ncomp)),n(ASTCKVEC(ncomp)),t(ASTCKVEC(nne));
  matrix gm(ASTCKMAT(ncomp,nne));
  
  nullv (flux);
  //  node coordinates
  Tt->give_node_coord2d (x,y,eid);

  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[lcid][lcid];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  //the first int. point
  Tm->computenlfluxes (lcid,ipp);    
  Tm->givefluxes (lcid,ipp,flux);
  
  //or volume average  
  /*  //  nodes on required elements
      Tt->give_elemnodes (eid,nodes);
      //  thickness
      Tc->give_thickness (eid,nodes,t);
      
      reallocv (RSTCKVEC(intordkm[lcid][lcid],gp));
      reallocv (RSTCKVEC(intordkm[lcid][lcid],w));
      
      //  coordinates and weights of integration points  
      gauss_points (gp.a,w.a,intordkm[lcid][lcid]);
      
      for (i=0;i<intordkm[lcid][lcid];i++){
      xi=gp[i];
      for (j=0;j<intordkm[lcid][lcid];j++){
      eta=gp[j];
      //  thickness
      thick = approx (xi,eta,t);
      
      Tm->computenlfluxes (lcid,ipp);
      
      Tm->givefluxes (lcid,ipp,fl);
      
      grad_matrix (gm,x,y,xi,eta,jac);
      
      //mtxv (gm,fl,contr);//this was for internal_fluxes = \int B^T D B r d\Omega
      
      //cmulv (jac*w[i]*w[j]*thick,fl);
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
      area = quadlin_edge_length (i,x,y);
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
