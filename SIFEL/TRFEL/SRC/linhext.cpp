/*
  File:     linhext.cpp
  Author:   Jaroslav Kruis, 31.3.2003
  Purpose:  3D hexahedral element with linear approximation functions
*/

#include "globalt.h"
#include "linhext.h"
#include "genfile.h"
#include "globmatt.h"

linhext::linhext (void)
{
  long i;

  //  number of nodes on element
  nne=8;
  //  number of edges
  ned=12;
  //  number of nodes on one edge
  nned=2;
  //  number of surfaces
  nsurf=6;
  //  number of nodes on one surface
  nnsurf=4;
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
    ordering[0][0]=1;  ordering[0][1]=2;  ordering[0][2]=3;  ordering[0][3]=4; 
    ordering[0][4]=5;  ordering[0][5]=6;  ordering[0][6]=7;  ordering[0][7]=8;
    
    dofe[0][0]=8;  intordkm[0][0]=2;  intordcm[0][0]=2;
    nip[0][0]=16;
    ndofe=8;  napfun=1;
    break;
  }

  case twomediacoup:{
    ordering[0][0]=1;  ordering[0][1]=3;  ordering[0][2]=5;  ordering[0][3]=7;
    ordering[0][4]=9;  ordering[0][5]=11; ordering[0][6]=13; ordering[0][7]=15;
    ordering[1][0]=2;  ordering[1][1]=4;  ordering[1][2]=6;  ordering[1][3]=8;
    ordering[1][4]=10; ordering[1][5]=12; ordering[1][6]=14; ordering[1][7]=16;

    intordkm[0][0]=2;  intordkm[0][1]=2;  intordkm[1][0]=2;  intordkm[1][1]=2;
    intordcm[0][0]=2;  intordcm[0][1]=2;  intordcm[1][0]=2;  intordcm[1][1]=2;
    
    if (Tp->savemode==0){
      nip[0][0]=16;      nip[0][1]=16;      nip[1][0]=16;      nip[1][1]=16;
    }
    if (Tp->savemode==1){
      nip[0][0]=16;      nip[0][1]=0;      nip[1][0]=0;      nip[1][1]=0;
    }
    
    dofe[0][0]=8;  dofe[0][1]=8;  dofe[1][0]=8;  dofe[1][1]=8;
    ndofe=16;  napfun=2;
    break;
  }

  case threemediacoup:{
    ordering[0][0]=1;  ordering[0][1]=4;  ordering[0][2]=7;  ordering[0][3]=10;
    ordering[1][0]=2;  ordering[1][1]=5;  ordering[1][2]=8;  ordering[1][3]=11;
    ordering[2][0]=3;  ordering[2][1]=6;  ordering[2][2]=9;  ordering[2][3]=12;

    ordering[0][4]=13; ordering[0][5]=16; ordering[0][6]=19; ordering[0][7]=22;
    ordering[1][4]=14; ordering[1][5]=17; ordering[1][6]=20; ordering[1][7]=23;
    ordering[2][4]=15; ordering[2][5]=18; ordering[2][6]=21; ordering[2][7]=24;

    intordkm[0][0]=2;  intordkm[0][1]=2;  intordkm[0][2]=2;
    intordkm[1][0]=2;  intordkm[1][1]=2;  intordkm[1][2]=2;
    intordkm[2][0]=2;  intordkm[2][1]=2;  intordkm[2][2]=2;

    intordcm[0][0]=2;  intordcm[0][1]=2;  intordcm[0][2]=2;
    intordcm[1][0]=2;  intordcm[1][1]=2;  intordcm[1][2]=2;
    intordcm[2][0]=2;  intordcm[2][1]=2;  intordcm[2][2]=2;

    if (Tp->savemode==0){
      nip[0][0]=16;      nip[0][1]=16;      nip[0][2]=16;
      nip[1][0]=16;      nip[1][1]=16;      nip[1][2]=16;
      nip[2][0]=16;      nip[2][1]=16;      nip[2][2]=16;
    }
    if (Tp->savemode==1){
      nip[0][0]=16;      nip[0][1]=0;      nip[0][2]=0;
      nip[1][0]=0;       nip[1][1]=0;      nip[1][2]=0;
      nip[2][0]=0;       nip[2][1]=0;      nip[2][2]=0;
    }
    
    dofe[0][0]=8;      dofe[0][1]=8;      dofe[0][2]=8;
    dofe[1][0]=8;      dofe[1][1]=8;      dofe[1][2]=8;
    dofe[2][0]=8;      dofe[2][1]=8;      dofe[2][2]=8;

    ndofe=24;  napfun=3;
    break;
  }
  default:{
    print_err("unknown type of transported matter is required",__FILE__,__LINE__,__func__);
  }
  }
}


linhext::~linhext (void)
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

void linhext::codnum (long *cn,long ri)
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
double linhext::element_volume (long eid)
{
  long i,j,k,ri,ci;
  ri=0;  ci=0;
  double xi,eta,zeta,ww1,ww2,ww3,jac,vol;
  vector x(nne),y(nne),z(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]),t(nne);
  
  Tt->give_node_coord3d (x,y,z,eid);
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  
  vol=0.0;
  for (i=0;i<intordkm[ri][ci];i++){
    xi=gp[i];  ww1=w[i];
    for (j=0;j<intordkm[ri][ci];j++){
      eta=gp[j];  ww2=w[j];
      for (k=0;k<intordkm[ri][ci];k++){
	zeta=gp[k]; ww3=w[k];
	
	jac_3d (jac,x,y,z,xi,eta,zeta);
	
	if (jac<0.0){
	  print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	  jac=fabs(jac);
	}
	
	jac*=ww1*ww2*ww3;
	vol+=jac;
      }
    }
  }
  
  return vol;
}

/**
   function approximates function defined by nodal values
   
   @param xi,eta,zeta - coordinates on element
   @param nodval - %vector of nodal values

   JK, 31.3.2002
*/
double linhext::approx (double xi,double eta,double zeta, vector &nodval)
{
  double f;
  vector bf(nne);
  bf_lin_hex_3d (bf.a,xi,eta,zeta);
    
  scprd (bf,nodval,f);

  return f;
}

/**
   function computes values in integration points from nodal values
   @param eid - element id

   JK, 31.3.2002
*/
void linhext::intpointval (long eid)
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
   function computes values in integration points from nodal values for PUC
   @param eid - element id

   TKr, 05/04/2011
*/
void linhext::intpointval_puc (long eid)
{
  long i, j, k, l, ii, jj, ipp;
  double xi, eta, zeta, val;
  vector r(ASTCKVEC(ndofe)), t(ASTCKVEC(nne)), gp, w;
  
  elemvalues_puc(eid, r);
  
  for (l=0;l<Tp->ntm;l++){
    
    //  nodal values of a single variable
    for (i=0;i<dofe[l][l];i++){
      t[i]=r[ordering[l][i]-1];
    }
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
	
	reallocv(RSTCKVEC(intordkm[ii][jj], gp));
	reallocv(RSTCKVEC(intordkm[ii][jj], w));
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
	
	reallocv(RSTCKVEC(intordcm[ii][jj], gp));
	reallocv(RSTCKVEC(intordcm[ii][jj], w));
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
  The function computes initial values in integration points from initial nodal values.
  (function interpolates the initial nodal values to integration points)

  @param eid - element id
 
  @return The function does not return anything but stores computed values at int. points of the given element.

  Created by TKo, 4.7.2018
*/
void linhext::initintpointval (long eid)
{
  long i,j,k,l,ii,jj,ipp,ndofn,cndofn;
  double xi,eta,zeta,val;
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
  
  for (l=0;l<Tp->ntm;l++){
    
    //  nodal values of a single variable
    for (i=0;i<dofe[l][l];i++){
      t[i]=r[ordering[l][i]-1];
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
	    for(k=0;k<intordkm[ii][jj];k++){
	      zeta=gp[k];
	      val = approx (xi,eta,zeta,t);
	      Tm->ip[ipp].av[l]=val;
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
   function computes gradients in integration points
   
   @param eid - element id

   JK, 31.3.2002
*/
void linhext::intpointgrad (long eid)
{
  long i,j,k,l,ii,jj,ipp;
  double xi,eta,zeta,jac;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne)),r(ASTCKVEC(ndofe));
  vector gp,w,grad(ASTCKVEC(ncomp)),t(ASTCKVEC(nne));
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
   function computes fluxes in integration points

   @param eid - element id

   JK, 12.8.2014, revised by TKr 29/09/2022
*/
void linhext::intpointflux (long eid)
{
  long i,j,k,l,ii,jj,ipp;
  
  for (l=0;l<Tp->ntm;l++){
    
    for (ii=0;ii<Tp->ntm;ii++){
      for (jj=0;jj<Tp->ntm;jj++){
	
	ipp=Tt->elements[eid].ipp[ii][jj];
	
	for (i=0;i<intordkm[ii][jj];i++){
	  for (j=0;j<intordkm[ii][jj];j++){
            for(k=0;k<intordkm[ii][jj];k++){
	      //  computation of correct fluxes
	      if (Tp->fluxcomp==1)
		Tm->computenlfluxes (l,ipp);
	      
	      ipp++;
	    }
	  }
	}
	for (i=0;i<intordcm[ii][jj];i++){
	  for (j=0;j<intordcm[ii][jj];j++){
            for(k=0;k<intordcm[ii][jj];k++){
	      //  computation of correct fluxes
	      if (Tp->fluxcomp==1)
		Tm->computenlfluxes (l,ipp);
	      
	      ipp++;
	    }
	  }
	}
	if (Tp->savemode==1)  break;
      }
      if (Tp->savemode==1)  break;
    }
  }
  
  /*  
      long i,j,k,l,ll,ii,jj,ipp;
      vector grad(ncomp),fl(ncomp),cfl(ncomp);
      matrix d(ncomp,ncomp);
      
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
   function computes average of flux of selected medium
   
   @param lcid - load case id, it is id of the medium
   @param eid - element id
   @param avfl - %vector of average flux
   
   JK, 12. 8. 2014
*/
void linhext::average_flux (long lcid,long eid,vector &avfl)
{
  long i,j,k,ii,jj,ipp;
  double xi,eta,zeta,jac,vol;
  vector x(nne),y(nne),z(nne),gp,w,grad(ncomp),fl(ncomp);
  matrix d(ncomp,ncomp);
  
  Tt->give_node_coord3d (x,y,z,eid);
  vol = elem_volume (eid);
  
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
	  for(k=0;k<intordkm[ii][jj];k++){
	    zeta=gp[k];
	    
	    jac_3d (jac,x,y,z,xi,eta,zeta);
	    
	    Tm->givefluxes (lcid,ipp,fl);
	    cmulv(jac*w[i]*w[j]*w[k]/vol,fl);
	    addv(fl,avfl,avfl);
	    
	    ipp++;
	  }
	}
      }
      if (Tp->savemode==1)
	break;
    }
  }
  
}

/**
   function approximates nodal values of array other to integration points
   
   @param eid - element id
   
   JK, 17.9.2005
*/
void linhext::intpointother (long eid)
{
  long i, j, k, l, ii, jj, ipp, ncompo, nodid;
  double xi, eta, zeta, val;
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
   
   MK, 26.2.2002
*/
void linhext::bf_matrix (matrix &n,double xi,double eta,double zeta)
{
  nullm (n);
  bf_lin_hex_3d (n.a,xi,eta,zeta);
}

/**
   function assembles gradient of %matrix of base functions
   
   @param gm - gradient %matrix
   @param x,y,z - array containing node coordinates
   @param xi,eta,zeta - natural coorodinates
   @param jac - Jacobian
   
   MK, 27.2.2002
*/
void linhext::grad_matrix (matrix &gm,vector &x,vector &y,vector &z,double xi,double eta,double zeta,double &jac)
{
  long i;
  vector dx(nne),dy(nne),dz(nne);

  dx_bf_lin_hex_3d (dx.a,eta,zeta);
  dy_bf_lin_hex_3d (dy.a,xi,zeta);
  dz_bf_lin_hex_3d (dz.a,xi,eta);

  derivatives_3d (dx,dy,dz,jac,x,y,z,xi,eta,zeta);
  
  nullm (gm);
  
  for (i=0;i<nne;i++){
    gm[0][i]=dx[i];
    gm[1][i]=dy[i];
    gm[2][i]=dz[i];
  }
  
}

/**
   function computes conductivity %matrix of one transported meduim
   
   @param lcid - load case id
   @param eid - element id
   @param ri, ci - row and column index
   @param km - conductivity %matrix
   
   JK, 9.3.2002
*/
void linhext::conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,j,k,ipp;
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
    ipp=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  
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
	Tm->matcond (d,ipp,ri,ci);
	
	if (jac<0.0){
	  print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	  jac=fabs(jac);
	}
	
	jac*=ww1*ww2*ww3;
	
	//  contribution to the conductivity matrix of the element
	bdbj (km.a,gm.a,d.a,jac,gm.m,gm.n);
	
	//convective terms
	reallocm(1,ncomp,d);
	
	Tm->matcond2(d,ipp,ri,ci);
	bf_matrix (n, xi, eta,zeta);
	bdbjac(km, n, d, gm, jac);	
	
	ipp++;
      }
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
   @param lm - L %matrix

   TKr, 05/04/2011
*/
void linhext::l_matrix (long /*lcid*/,long eid,long ri,long ci,matrix &lm)
{
  long i,j,k,ii;
  double xi,eta,zeta,ww1,ww2,ww3,jac;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]),t(nne);
  matrix gm(ncomp,dofe[ri][ci]),d;

  Tt->give_elemnodes (eid,nodes);
    
  Tt->give_node_coord3d (x,y,z,eid);
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
      for (k=0;k<intordkm[ri][ci];k++){
	zeta=gp[k]; ww3=w[k];
	//  matrix of gradients
	grad_matrix (gm,x,y,z,xi,eta,zeta,jac);
	
	//  matrix of conductivity of the material
	reallocm(ncomp,ncomp,d);
	Tm->matcond (d,ii,ri,ci);
	
	if (jac<0.0){
	  print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	  jac=fabs(jac);
	}
	
	jac*=ww1*ww2*ww3;
	
	//  contribution to the L matrix of the element
	mxm(d,gm,lm);
	cmulm(jac,lm);
	
	ii++;
      }
    }
  }
}



/**
   function computes L^T (L transposed) %matrix

   L^T = \int_{\Omega} B^T D {\rm d} \Omega

   @param lcid - load case id
   @param eid - element id
   @param lm - L^T %matrix

   TKr, 05/04/2011
*/
void linhext::l_t_matrix (long /*lcid*/,long eid,long ri,long ci,matrix &lm)
{
  long i,j,k,ii;
  double xi,eta,zeta,ww1,ww2,ww3,jac;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]),t(nne);
  matrix gm(ncomp,dofe[ri][ci]),d;

  Tt->give_elemnodes (eid,nodes);
    
  Tt->give_node_coord3d (x,y,z,eid);
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
      for (k=0;k<intordkm[ri][ci];k++){
	zeta=gp[k]; ww3=w[k];
	//  matrix of gradients
	grad_matrix (gm,x,y,z,xi,eta,zeta,jac);
	
	//  matrix of conductivity of the material
	reallocm(ncomp,ncomp,d);
	Tm->matcond (d,ii,ri,ci);
	
	if (jac<0.0){
	  print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	  jac=fabs(jac);
	}
	
	jac*=ww1*ww2*ww3;
	
	//  contribution to the L matrix of the element
	mtxm(gm,d,lm);
	cmulm(jac,lm);
	
	ii++;
      }
    }
  }
}



/**
   function computes capacity %matrix of one transported matter
   finite element with tri-linear approximation functions
   
   @param eid - number of element
   @param ri,ci - row and column indices of the computed block in the resulting %matrix
   @param cm - capacity %matrix
   
   MK, 26.2.2002
*/
void linhext::capacity_matrix (long eid,long ri,long ci,matrix &cm)
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
	
	if (jac<0.0){
	  print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	  jac=fabs(jac);
	}
	
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
   
   JK, 9.3.2002
*/
void linhext::quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci)
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
	
	if (jac<0.0){
	  print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	  jac=fabs(jac);
	}
	
	jac*=w1*w2*w3;
	
	nnj (nm.a,n.a,jac,n.m,n.n);
      }
    }
  }
  mxv (nm,nodval,v);  addv (sv,v,sv);
  
}


/**
   function computes internal fluxes of 3D problems for one transported matter
   
   @param lcid - number of load case
   @param eid - number of element
   @param ifl - %vector of internal fluxes
   
   JK, 31.3.2002
*/
void linhext::internal_fluxes (long lcid,long eid,vector &ifl)
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
	
	if (jac<0.0){
	  print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	  jac=fabs(jac);
	}
	
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
   
   JK, 9.3.2002
*/
void linhext::res_conductivity_matrix (long eid,long /*lcid*/,matrix &km)
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
   
   TKr, 12/5/2011
*/
void linhext::volume_rhs_vector (long /*lcid*/,long eid,long ri,long ci,vector &vrhs)
{
  long i,j,k,ii;
  double xi,eta,zeta,ww1,ww2,ww3,jac;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]);
  matrix gm(ncomp,dofe[ri][ci]),d;
  matrix km(dofe[ri][ci],3);
  
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
	
	Tm->volume_rhs (d,ii,ri,ci,ncomp);
	
	if (jac<0.0){
	  print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	  jac=fabs(jac);
	}
	
	jac*=ww1*ww2*ww3;
	
	//  contribution to the volume_rhs integral of the element
	nnjac (km,gm,d,jac);
	
	ii++;
      }
    }
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

   TKr, 12/5/2011
*/
void linhext::res_volume_rhs_vector (vector &f,long eid,long /*lcid*/)
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
   function assembles resulting element L %matrix

   L^T = \int_{\Omega} B^T D {\rm d} \Omega

   @param eid - element id
   @param lcid - load case id
   @param lm - resulting L %matrix of one element

   TKr, 05/04/2011
*/
void linhext::res_l_matrix (long eid,long /*lcid*/,matrix &lm)
{
  long i,j,*rcn,*ccn;
  matrix lkm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      rcn = new long [3];
      ccn = new long [dofe[i][j]];
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

   TKr, 05/04/2011
*/
void linhext::res_l_t_matrix (long eid,long /*lcid*/,matrix &lm)
{
  long i,j,*rcn,*ccn;
  matrix lkm;
  
  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      rcn = new long [dofe[i][j]];
      ccn = new long [3];
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

   TKr, 05/04/2011
*/
void linhext::averd_matrix (long eid,matrix &lm)
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

      mat_localize (lm,d,rcn,ccn);
    }
  }
  
  delete [] rcn;  delete [] ccn;
}





/**
   function assembles average C %matrix

   @param eid - element id
   @param lm - resulting C %matrix of one element

   TKr, 05/04/2011
*/
void linhext::averc_matrix (long eid,matrix &lm)
{
  long i,j,ii;
  double c,rho;
  ivector nodes(nne);
  vector dens(nne);

  Tt->give_elemnodes (eid,nodes);
  Tc->give_density (eid,nodes,dens);

  for (i=0;i<ntm;i++){
    for (j=0;j<ntm;j++){
      
      if (Tp->savemode==0)
	ii=Tt->elements[eid].ipp[i][j];
      if (Tp->savemode==1)
	ii=Tt->elements[eid].ipp[0][0];
      
      //  coefficient of capacity of material
      c = Tm->capcoeff (ii,i,j);
      //temporarily
      rho = approx (0.0,0.0,0.0,dens);
      
      lm[i][j] = rho*c*Tm->ip[ii].av[j];
    }
  }
}




/**
   function assembles volume of one element

   @param eid - element id

   TKr, 05/04/2011
*/
double linhext::elem_volume (long eid)
{
  double jac,volume;
  vector x(nne),y(nne),z(nne);
  matrix gm(ncomp,dofe[0][0]);

  Tt->give_node_coord3d (x,y,z,eid);
  
  //  matrix of gradients - jacobian
  grad_matrix (gm,x,y,z,1.0,1.0,1.0,jac);
  
  if (jac<0.0){
    print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
    jac=fabs(jac);
  }
  
  //  volume is equal to eight times jac of the element
  volume = jac*8.0;

  return volume;
}




/**
   function assembles resultant element capacity %matrix

   @param eid - element id
   @param cm - resulting capacity %matrix of one element
   
   JK, 9.3.2002
*/
void linhext::res_capacity_matrix (long eid,matrix &cm)
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
   @param leid - loaded element id
   
   JK, 6.1.2002
*/
void linhext::res_convection_vector (vector &f,long lcid,long eid,long leid)
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
void linhext::res_transmission_vector (vector &f,long lcid,long eid,long leid)
{
  long *cn,i;
  vector lf;
  
  cn = new long [dofe[lcid][lcid]];
  //  code numbers
  codnum (cn,lcid);
  reallocv (dofe[lcid][lcid],lf);
  
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
void linhext::res_quantity_source_vector (vector &sv,vector &nodval,long lcid,long eid)
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
   function assembles resulting element internal fluxes %vector

   @param eid - element id
   @param elemif - resulting internal fluxes %vector of one element
   
   JK, 6.1.2002
*/
void linhext::res_internal_fluxes (long eid, vector &elemif)
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
  nodalderivatives (eid, tdnv);
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
double linhext::total_integral(long eid,vector &nodval)
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

	if (jac<0.0){
	  print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	  jac=fabs(jac);
	}
	
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
   
   TKr, 28.2.2004
*/
void linhext::res_boundary_flux (vector &f,long lcid,long eid,long leid)
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
   function computes gradients in nodes of element

   @param eid - element id
   
   24. 3. 2017, JK
*/
void linhext::nod_grads_ip (long eid)
{
  long i,j,k,ipp;
  double vol;
  ivector ipnum(nne),nod(nne);
  vector grad(ncomp);

  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  nodip_linhex (ipp,intordkm[0][0],ipnum);
  
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
	vol = element_volume (eid);
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
void linhext::nod_fluxes_ip (long eid)
{
  long i,j,k,ipp;
  double vol;
  ivector ipnum(nne),nod(nne);
  vector flux(ncomp);
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  nodip_linhex (ipp,intordkm[0][0],ipnum);
  
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
	vol = element_volume (eid);
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
void linhext::nod_others_comp (long /*lcid*/,long eid,long ri,long ci)
{
  long i,j,k,ipp,ncompother,ii;
  vector h,r(ASTCKVEC(ndofe));
  ivector nod(ASTCKIVEC(nne));
  double vol,other;
  
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
	vol = element_volume (eid);
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
void linhext::nod_eqother_ip (long eid)
{
  long i,j,k,ipp,ncompo;
  double vol;
  ivector ipnum(nne),nod(nne);
  vector eqother;
  
  //  numbers of integration points closest to nodes
  //  (function is from the file GEFEL/ordering.cpp)
  ipp=Tt->elements[eid].ipp[0][0];
  nodip_linhex (ipp,intordkm[0][0],ipnum);
  
  //  node numbers of the element
  Tt->give_elemnodes (eid,nod);
  
  for (i=0;i<nne;i++){
    ncompo = Tm->givencompeqother (ipnum[i],0);
    reallocv (ncompo,eqother);
    Tm->giveeqother (ipnum[i],0,ncompo,eqother.a);
    
    //  storage of eqother components to the node
    j=nod[i];
    for (k=0; k<ncompo; k++){
      if (Tp->eqotheraver==0 || Tp->eqotheraver==1)
	Tt->nodes[j].storeeqother (k,eqother(k));
      if (Tp->eqotheraver==2){
	vol = element_volume (eid);
	Tt->nodes[j].storeeqother (k,vol,eqother(k));
      }
    }
  }
}
















/**
   function computes nodal fluxes from boundary values
   
   \int_{Gamma_2} N^T N dGamma * nodal_flux_values

   @param v - array of nodal fluxes
   @param lcid - 
   @param eid - element id
   @param leid - id of loaded element
   @param ri,ci - row and column indices
   
   JK, 19.8.2004
*/
void linhext::convection_vector (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  vector x(nne),y(nne),z(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]);
  vector list(4*nsurf),trc(4*nsurf),nodval(nne),av(nne),coef(nne);
  matrix km(nne,nne);
  
  nullv (v);
  
  //  node coordinates on element
  Tt->give_node_coord3d (x,y,z,eid);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[ri][ci]);
  
  //  boundary conditions
  bc = new bocontypet [nsurf];
  //  indicators of boundary conditions
  Tb->lc[lcid].elemload[leid].give_bc (bc);
  
  //  auxiliary coefficients, necessary for function surface_integral
  for (i=0;i<4*nsurf;i++){
    trc[i]=1.0;
  }

  //  nodal values on all surfaces
  Tb->lc[lcid].elemload[leid].give_nodval (lcid,list);
  
  //  integration point number
  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[ri][ci];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  //  loop over surfaces
  for (i=0;i<nsurf;i++){
    
    if (bc[i]==2 || bc[i]==3 || bc[i]==4 || bc[i]==5){
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
   
   @param lcid - load case id
   @param eid - element id
   @param ri,ci - row and column indices
   @param km - transmission part of conductivity %matrix
   
   JK, 19.8.2004
*/
void linhext::transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km)
{
  long i,ipp,leid;
  bocontypet *bc;
  vector x(nne),y(nne),z(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]);
  vector trc(4*nsurf),coeff(nne);
  
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
  bc = new bocontypet[nsurf];
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
    
    if (bc[i]==4 || bc[i]==5 || bc[i]>10){
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

   \int_{Gamma_3} N^T c_{tr} N dGamma * nodal_external_value
   
   @param tmv - transmission %vector
   @param lcid - load case id
   @param eid - element id
   @param leid - loaded element id
   @param ri,ci - row and column indices
   
   JK, 5.10.2001
   TKr, 30.1.2004 - new added
*/
void linhext::transmission_vector (vector &v,long lcid,long eid,long leid,long cid)
{
  long i,ipp;
  bocontypet *bc;
  vector x(nne),y(nne),z(nne),w(intordkm[lcid][cid]),gp(intordkm[lcid][cid]);
  vector list(4*nsurf),trc(4*nsurf),trr(4*nsurf),nodval(nne),av(nne),coef(nne);
  matrix km(nne,nne);
  
  nullv (v);
  
  //  node coordinates on element
  Tt->give_node_coord3d (x,y,z,eid);
  
  //  coordinates and weights of integration points
  gauss_points (gp.a,w.a,intordkm[lcid][cid]);
  
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
    
    if (bc[i]==4 || bc[i]==5 || bc[i]>10){
      //  transformation of nodal values
      transf_val (i,nodval,list,trc,trr,eid,lcid,cid,ipp,bc);
      //  determination of transmission coefficients on edge
      transf_coeff (i,coef,trc,eid,lcid,cid,ipp,bc);
      
      nullm (km);
      //  matrix obtained from integration over surfaces
      surface_integral (i,x,y,z,intordkm[lcid][cid],gp,w,coef,km);
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
   
   TKr, 28.2.2004
*/
void linhext::boundary_flux (vector &v,long lcid,long eid,long leid,long ri,long ci)
{
  long i,ipp;
  bocontypet *bc;
  vector x(nne),y(nne),z(nne),w(intordkm[ri][ci]),gp(intordkm[ri][ci]);
  vector list(4*nsurf),trc(4*nsurf),trr(4*nsurf),nodval(nne),av(nne),coef(nne);
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
    
    if (bc[i]==4 || bc[i]==5 || bc[i]>10){
      //  transformation of nodal values
      transf_flux (i,nodval,list,trc,trr,eid,ri,ci,ipp,bc);
      //  determination of transmission coefficients on edge
      transf_coeff (i,coef,trc,eid,ri,ci,ipp,bc);

      nullm (km);
      //  matrix obtained from integration over surface
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
void linhext::surface_integral (long surf,vector &x,vector &y,vector &z,long intord,vector &gp,vector &w,
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
void linhext::transf_flux (long surf,vector &coeff,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double tr,new_flux;
  ivector nodes(nne),surfnod(nnsurf);

  //  zeroing of array of coefficients
  nullv (coeff);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required surface
  linhexahedral_surfnod (surfnod.a,surf);

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


void linhext::transf_coeff (long surf,vector &coeff,vector &list,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double new_coeff;
  ivector nodes(nne),surfnod(nnsurf);

  //  zeroing of array of coefficients
  nullv (coeff);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required surface
  linhexahedral_surfnod (surfnod.a,surf);

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
void linhext::transf_val (long surf,vector &nodval,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc)
{
  long i,j,k;
  double tr,new_nodval;
  ivector nodes(nne),surfnod(nnsurf);
  
  //  zeroing of array of nodal values
  nullv (nodval);
  
  //  global numbers of nodes of element
  Tt->give_elemnodes (eid,nodes);

  //  local numbers of nodes on required surface
  linhexahedral_surfnod (surfnod.a,surf);
  
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
void linhext::surfnodeval (long surf,vector &nodval,vector &list)
{
  long i,j;
  ivector surfnod(nnsurf);
  
  nullv (nodval);
  linhexahedral_surfnod (surfnod.a,surf);
  
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
void linhext::higher_to_lower_level (long eid,long *counter,double *buff)
{
  long i,ipp;
  vector gr(3);
  
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
    buff[counter[0]]=gr[2];
    counter[0]++;
    
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
long linhext::ipcoord (long eid, long ipp, long ri, long ci, vector &coord)
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
   
   @retval ncoord - function returns natural coordinates of integration points in the %vector ncoord.
   @retval 0 - on success
   @retval 1 - integration point with ipp was not found

   TKo, 12.2016
*/
long linhext::ipncoord (long eid, long ipp, vector &ncoord)
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



/**
   Function computes the total flux through a surface of an element
   
   @param lcid - load case id
   @param eid - element id (id in the list of all elements)
   @param beid - id in the list of selected elements
   @param flux - array of fluxes computed
   
   TKr 26/09/2022 according to JK
*/
void linhext::surface_flux (long lcid,long eid,long beid,double *fluxes)
{

  long i,ipp,fid;
  double q,qn,area;
  bocontypet *bc;
  vector x(ASTCKVEC(nne)),y(ASTCKVEC(nne)),z(ASTCKVEC(nne));
  vector bb(ASTCKVEC(4)),cc(ASTCKVEC(4)),dd(ASTCKVEC(4)),volcoord(ASTCKVEC(4));
  vector r(ASTCKVEC(ndofe)),t(ASTCKVEC(nne)),grad(ASTCKVEC(ncomp)),flux(ASTCKVEC(ncomp)),fl(ASTCKVEC(ncomp)),n(ASTCKVEC(ncomp));
  matrix gm(ASTCKMAT(ncomp,nne)),d(ASTCKMAT(ncomp,ncomp));
  
  
  //  node coordinates
  Tt->give_node_coord3d (x,y,z,eid);

  if (Tp->savemode==0)
    ipp=Tt->elements[eid].ipp[lcid][lcid];
  if (Tp->savemode==1)
    ipp=Tt->elements[eid].ipp[0][0];
  
  //the first int. point
  Tm->computenlfluxes (lcid,ipp);    
  Tm->givefluxes (lcid,ipp,flux);

  //  indicators of boundary conditions
  bc = new bocontypet [nsurf];
  Tb->bf[lcid].elemload[beid].give_bc (bc);
  
  //  loop over surfaces
  for (i=0;i<nsurf;i++){
    
    if (bc[i]==1){
      //  function constructs the outer normal vector
      hexahedra_normal_vectors (i,x,y,z,n);
      //  area of the i-th element surface
      area = hexahedra_surface_areas (i,x,y,z);
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
