#include "coupmatu.h"
#include "globalc.h"
#include "mechtop.h"
#include "constrelcu.h"
#include "onemediumc.h"
#include "twomediac.h"
#include "threemediac.h"

coupmatu::coupmatu (void)
{
  tnip=0;
  ip=NULL;
  
  itrmc=NULL;
  sejtkrmc=NULL;
  consol_awf1c = NULL;
  consol_wf1c = NULL;
  consol_wf2c = NULL;
  consol_awf2c = NULL;
  consol_hwf2c = NULL;
  consol_hawf3c = NULL;

  itrmc=NULL;
  concretec=NULL;
  baroghelc = NULL;
  C60baroghelc = NULL;
  C30baroghelc = NULL;
  o30bazantc = NULL;
  C60bazantc = NULL;

  tenchc = NULL;

}

coupmatu::~coupmatu (void)
{
  delete [] ip; 

  delete [] itrmc;
  delete [] sejtkrmc;
  delete [] consol_awf1c;
  delete [] consol_wf1c;
  delete [] consol_wf2c;
  delete [] consol_awf2c;
  delete [] consol_hwf2c;
  delete [] consol_hawf3c;
  delete [] concretec;
  delete [] baroghelc;
  delete [] C60baroghelc;
  delete [] C30baroghelc;
  delete [] o30bazantc;
  delete [] C60bazantc;

  delete [] tenchc;  
}

void coupmatu::ipalloc (void)
{
  long i,j,ii,jj,ne,ntm,gdim,ippu,mnb,nip,ncompstr;
  strastrestate ssst;
  
  //  number of transported media
  ntm=Tp->ntm;
  //  number of geometric dimensions (1D, 2D or 3D)
  gdim = Tp->gdim;
  //  number of elements
  ne=Ct->ne;
  
  for (i=0;i<ne;i++){
    //  number of mechanical blocks
    mnb=Ct->give_mnb (i);
    //  number of strain/stress components
    ncompstr=Ct->give_ncompstr (i);
    
    for (ii=0;ii<mnb;ii++){
      //  type of strain/stress state
      ssst=Ct->give_ssst (i,ii);

      for (jj=0;jj<ntm;jj++){
	ippu=Ct->elements[i].ippu[ii][jj];
	nip=Ct->give_upper_nip (i,ii,jj);

	for (j=0;j<nip;j++){
	  ip[ippu].ncompstr=ncompstr;
	  ip[ippu].ssst=ssst;
	  ip[ippu].stress = new double [ncompstr];
	  ip[ippu].strain = new double [ncompstr];
	  
	  if ((ssst == planestrain) || (ssst == planestress))
	    ip[ippu].ncompstr=4;
	  
	  ippu++;
	}
      }
    }
  }
  
  for (i=0;i<tnip;i++){
    ip[i].av = new double [ntm];
    ip[i].pv = new double [ntm];
    ip[i].grad = new double* [ntm];
    ip[i].fluxes = new double* [ntm];
    
    for (j=0;j<ntm;j++){
      ip[i].grad[j] = new double [gdim];
      ip[i].fluxes[j] = new double [gdim];
    }
  }
  
}


/**
   function computes pointers on integration points
   
   offdiagonal submatrices are divided into blocks
   each block has row and column index
   number of rows is equal to number of blocks in mechanics
   number of columns is equal to number of transported matters
   
   25.2.2003
*/
long coupmatu::intpnum (void)
{
  long i,j,k,ne,n,nb;
  elemtypec te;
  ne=Ct->ne;
  
  n=0;
  for (i=0;i<ne;i++){
    te = Ct->give_elem_type (i);
    nb=Mt->give_nb (i);
    Ct->elements[i].ippu = new long* [nb];
    Ct->elements[i].nb = nb;
    for (j=0;j<nb;j++){
      Ct->elements[i].ippu[j] = new long [Tp->ntm];
      for (k=0;k<Tp->ntm;k++){
	Ct->elements[i].ippu[j][k]=n;
	n+=Ct->give_upper_nip (i,j,k);
      }
    }
  } 
  return n;
}

/**
   function allocates integration points
*/
void coupmatu::intpointalloc ()
{
  ip = new intpointsc [tnip];
  ipalloc ();
}

/**
   function initializes basic data in integration points
 */
void coupmatu::intpointinit ()
{
  long i,j,k,ii,jj,nb,ntm,nip,ipp;
  
  for (i=0;i<Tt->ne;i++){
    ntm=Tp->ntm;  k=0;
    nb=Mt->give_nb (i);
    for (ii=0;ii<nb;ii++){
      for (jj=0;jj<ntm;jj++){
	nip=Ct->give_upper_nip(i,ii,jj);
	ipp=Ct->elements[i].ippu[ii][jj];
	for (j=0;j<nip;j++){
	  ip[ipp].tm   = Ct->elements[i].tmu[k];
	  ip[ipp].idm  = Ct->elements[i].idmu[k];
	  ipp++;
	}
	k++;
      }
    }
  }
}


/**
   function initiate values at integration points
*/
void coupmatu::initmaterialmodels (void)
{
  long i,j,ipp,nip;
  
  for (i=0;i<Ct->ne;i++){
    if (Gtt->leso[i]==1){
      ipp=Ct->elements[i].ippu[0][0];
      nip=Ct->give_upper_tnip (i);
      for (j=0;j<nip;j++){
	initvalues (ipp,0,0);
	ipp++;
      }
    }
  }
}


/**
   function initiate material models and other values
   
   @param ipp - integration point pointer
   @param im  - index of material type for given ip
   @param ido - index in array eq_other   

   07/05/2010, TKr
*/
void coupmatu::initvalues (long /*ipp*/,long /*im*/,long /*ido*/)
{
  switch (Cp->tmatt){
  case mech_onemedium:{
    medc1 mc1;
    
    //mc1.initvalues (ipp);
    break;
  }
  case mech_twomedia:{
    medc2 mc2;
    
    //mc2.initvalues (ipp);
    break;
  }
  case mech_threemedia:{
    medc3 mc3;
    
    //mc3.initvalues (ipp);
    break;
  }
  default:{
    print_err("unknown media type is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
  Function updates values of internal variables at integration points,
  it is used for nonlinear computations.

  @return The function does not return anything.
  
  11/05/2010 TKr 
*/
void coupmatu::updateipval (void)
{
  long i,j,ipp,nip;
  
  //updates only used elements(int. points)
  for (i=0;i<Ct->ne;i++){
    if (Gtu->leso[i]==1){
      ipp=Ct->elements[i].ippu[0][0];
      nip=Ct->give_upper_tnip (i);
      for (j=0;j<nip;j++){
	updateipvalmat (ipp,0,0);
	ipp++;
      }
    }
  }
}



/**
  Function updates values of internal variables (eqother array) at integration points,
  it is used for nonlinear computations.
   
  @param ipp - number of integration point
  @param im - index of material type
  @param ido - index in array other
  
  11/05/2010 TKr
*/
void coupmatu::updateipvalmat (long /*ipp*/,long /*im*/,long /*ido*/)
{
  switch (Cp->tmatt){
  case mech_onemedium:{
    medc1 mc1;
    
    //mc1.updatevalues (ipp);
    break;
  }
  case mech_twomedia:{
    medc2 mc2;
    
    //mc2.updatevalues (ipp);
    break;
  }
  case mech_threemedia:{
    medc3 mc3;
    
    //mc3.updatevalues (ipp);
    break;
  }
  default:{
    print_err("unknown media type is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
   function computes actual stiffness %matrix of the material
   in the required integration point
   
   @param d - actual stiffness %matrix
   @param ipp - number of integration point
   @param im  - index of the material in the ipp tm and idm arrays. The standard value is zero.
   @param ido - index of internal variables in the ip's ipp other array
   
   TKr, 17.7.2001
*/
void coupmatu::matstiff (matrix &d,long ipp)
{
  switch (Cp->tmatt){
  case mech_onemedium:
  case mech_twomedia:
  case mech_threemedia:{
    state_eqcu c;
    c.matstiff (d,ip[ipp].ssst,ipp);//this is elastic stiffness matrix - has to be changed to actual (tangent or secant)
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function matstiff (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }  
}



/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ipp - number of integration point
   @param ri,ci - row and column indices
   
   17.07.2005
*/
void coupmatu::matcond (matrix &d,long ipp,long ri,long ci)
{
  switch (Cp->tmatt){//  transported matter
  case mech_onemedium:{
    medc1 mc1;
    
    mc1.matcond_u(d,ri,ci,ipp);
    break;    
  }
  case mech_twomedia:{
    medc2 mc2;
    
    mc2.matcond_u(d,ri,ci,ipp);    
    break;    
  }
  case mech_threemedia:{
    medc3 mc3;
    
    mc3.matcond_u(d,ri,ci,ipp);    
    break;    
  }
  default:{
    fprintf (stderr,"\n unknown number of transported media is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}

/**
   function computes capacity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ipp - number of integration point
   @param ri,ci - row and column indices
   
   17.07.2005
*/
void coupmatu::matcap (matrix &d,long ipp,long ri,long ci)
{   
  switch (Cp->tmatt){//  transported matter
  case mech_onemedium:{
    medc1 mc1;

    mc1.matcap_u(d,ri,ci,ipp);
    break;    
  }
  case mech_twomedia:{
    medc2 mc2;

    mc2.matcap_u(d,ri,ci,ipp);
    break;    
  }
  case mech_threemedia:{
    medc3 mc3;

    mc3.matcap_u(d,ri,ci,ipp);
    break;    
  }
  default:{
    fprintf (stderr,"\n unknown number of transported media is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }  
}




/**
   function computes right-hand side matrix of the material - contribution into right-hand side from volume integrals
   in the required integration point
   
   @param d - right-hand side %matrix of material
   @param ipp - number of integration point
   @param ri,ci - row and column indices
   @param ncomp - number of components
   
   25.9.2001
*/
void coupmatu::volume_rhs1 (matrix &d,long ipp,long ri,long ci,long /*ncomp*/)
{
  switch (Cp->tmatt){//  transported matter
  case mech_onemedium:{
    medc1 mc1;
    
    mc1.rhs_u1(d,ri,ci,ipp);    
    break;    
  }
  case mech_twomedia:{
    medc2 mc2;
    
    mc2.rhs_u1(d,ri,ci,ipp);    
    break;    
  }
  case mech_threemedia:{
    medc3 mc3;
    
    mc3.rhs_u1(d,ri,ci,ipp);    
    break;    
  }
  default:{
    fprintf (stderr,"\n unknown number of transported media is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}



/**
   function computes right-hand side matrix of the material - contribution into right-hand side from volume integrals
   in the required integration point
   
   @param d - right-hand side %matrix of material
   @param ipp - number of integration point
   @param ri,ci - row and column indices
   @param ncomp - number of components
   
   25.9.2001
*/
void coupmatu::volume_rhs2 (matrix &d,long ipp,long ri,long ci)
{
  fillm(0.0,d);

  switch (Cp->tmatt){//  transported matter
  case mech_onemedium:{
    medc1 mc1;

    mc1.rhs_u2(d,ri,ci,ipp);
    break;    
  }
  case mech_twomedia:{
    medc2 mc2;

    mc2.rhs_u2(d,ri,ci,ipp);
    break;    
  }
  case mech_threemedia:{
    medc3 mc3;
    
    mc3.rhs_u2(d,ri,ci,ipp);
    break;    
  }
  default:{
    fprintf (stderr,"\n unknown number of transported media is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}




/**
   function computes correct stresses

   @param ri,ci - row and column indices
   @param ipp - number of integration point
   
   22/06/2018, TKr
*/
void coupmatu::computenlstresses (matrix &d,long /*ri*/, long /*ci*/, long ipp)
{
  long j,k,m,n;
  
  m = d.m;     //number of strain/stress components
  n = Tp->ntm; //number of transported media

  for (j=0;j<m;j++)//null flux
    ip[ipp].stress[j] = 0.0;

  //matcond contribution = pore pressures contrib.
  
  for (j=0;j<n;j++){//number of transported media
    matcond (d,ipp,0,j);//ri = 0; ci = j
    for (k=0;k<m;k++){//number of strain/stress components
      ip[ipp].stress[k] += d[k][0]*ip[ipp].av[j];
    }
  } 
  
  //matcap contribution = pore pressures rate contrib. - not included in this approach
  


  //special case for matcond contribution from temperature is not implemented yet, see below the source commented:

  /*     
	 long j,k,n,m;
	 
	 m = d.m;
	 n = d.n;
	 
	 switch (Cp->tmatt){//  transported matter
	 case mech_onemedium:{
	 
	 switch (Cp->mednam){//  names of METR transported media
	 case mech_heat:{
	 
	 for (j=0;j<m;j++)//null flux
	 ip[ipp].stresses[j] = 0.0;
	 
	 //for temperature
	 matrix s(d.m,d.m),sd(d.m,d.n);
	 
	 ut.matcond(sd,ipp);
	 matstiff (s,ipp);
	 mxm(s,sd,d);
	 for (j=0;j<m;j++)
	 for (k=0;k<n;k++)
	 ip[ipp].stresses[j] -= d[j][k]*ip[ipp].av[0];//t//sign+
	 
	 break;
	 }
	 default:{
	 fprintf (stderr,"\n unknown name of transported media is required in function (%s, line %d).\n",__FILE__,__LINE__);
	 }
	 }
	 
	 break;    
	 }
	 case mech_threemedia:{
	 
	 switch (Cp->mednam){//  names of METR transported media
	 case mech_heat_moisture:{
	 
	 for (j=0;j<m;j++)//null flux
	 ip[ipp].stresses[j] = 0.0;
	 
	 //for capillary pressure pc
	 uc.matcond(d,ipp);
	 for (j=0;j<m;j++)
	 for (k=0;k<n;k++)
	 ip[ipp].stresses[j] -= d[j][k]*ip[ipp].av[0];//pc//sign+
	 
	 //for gas pressure pg
	 ug.matcond(d,ipp);
	 for (j=0;j<m;j++)
	 for (k=0;k<n;k++)
	 ip[ipp].stresses[j] -= d[j][k]*ip[ipp].av[1];//pg//sign+
	 
	 //for temperature
	 matrix s(d.m,d.m),sd(d.m,d.n);
	 
	 ut.matcond(sd,ipp);
	 matstiff (s,ipp);
	 mxm(s,sd,d);
	 for (j=0;j<m;j++)
	 for (k=0;k<n;k++)
	 ip[ipp].stresses[j] -= d[j][k]*ip[ipp].av[2];//t//sign+
	 
	 break;
	 }
	 default:{
	 fprintf (stderr,"\n unknown name of transported media is required in function (%s, line %d).\n",__FILE__,__LINE__);
	 }
	 }
	 
	 break;    
	 }
	 default:{
	 fprintf (stderr,"\n\n unknown material is required in function computenlstresses (%s, line %d).\n",__FILE__,__LINE__);
	 abort ();
	 }
	 }
  */
}


/**
   function stores components of stresses to integration points
   
   @param ipp - number of integration point
   @param fi - first read index
   @param gr - array containing components of stresses
   
   4.12.2002
*/
void coupmatu::storestresses_cmu (long ipp,long fi,vector &fl)
{
  long i,j,ncomp;
  ncomp=fl.n+fi;
  
  j=0;
  for (i=fi;i<ncomp;i++){
    ip[ipp].stress[i]=fl[j];
    j++;
  }
  
}

/**
   function returns components of flux (stress)
   
   @param ipp - number of integration point
   @param fi - first read index
   @param fl - array containing components of flux

   22.12.2002
*/
void coupmatu::givestresses_cmu (long ipp,long fi,vector &fl)
{
  long i,j,ncomp;
  ncomp=fl.n+fi;
  
  j=0;
  for (i=fi;i<ncomp;i++){
    fl[j] = ip[ipp].stress[i];
    j++;
  }  
}

/**
   function stores computed strains
   
   @param ipp - integration point pointer
   @param fi - first index
   @param eps - %vector containing strain components
   
   17.7.2001
*/
void coupmatu::storestrain_cmu (long ipp,long fi,vector &eps)
{
  long i,j,ncomp=eps.n;
  j=0;
  for (i=fi;i<fi+ncomp;i++){
    ip[ipp].strain[i]=eps[j];
    j++;
  }
}

/**
   function restores strains to %vector eps
   
   @param ipp - integration point pointer
   @param eps - %vector containing strain components
   
   5.8.2001
*/
void coupmatu::givestrain_cmu (long ipp,long fi,vector &eps)
{
  long i,j,ncomp=eps.n;
  j=0;
  for (i=fi;i<fi+ncomp;i++){
    eps[j]=ip[ipp].strain[i];
    j++;
  }
}


/**
   function stores components of gradients to integration points
   
   @param lcid - load case id = number of unknown
   @param ipp - number of integration point
   @param gr - array containing components of gradient
   
   4.12.2002
*/
void coupmatu::storegrad_cmu (long lcid,long ipp,vector &gr)
{
  long i,j,ncomp;
  ncomp=gr.n;
  
  j=0;
  for (i=0;i<ncomp;i++){
    ip[ipp].grad[lcid][i]=gr[j];
    j++;
  }
  
}


/**
   function returns components of gradients
   
   @param lcid - load case id = number ot unknown
   @param ipp - number of integration point
   @param gr - array containing components of flux

   14/05/2010, TKr
*/
void coupmatu::givegrad_cmu (long lcid,long ipp,vector &gr)
{
  long i,ncomp;
  ncomp=gr.n;
  
  for (i=0;i<ncomp;i++){
    gr[i] = ip[ipp].grad[lcid][i];
  }
}

/**
   function stores components of flux to integration points
   
   @param lcid - load case id = number of unknown
   @param ipp - number of integration point
   @param fl - array containing components of flux
   
   14/05/2010, TKr
*/
void coupmatu::storeflux_cmu (long lcid,long ipp,vector &fl)
{
  long i,ncomp;
  ncomp=fl.n;
  
  for (i=0;i<ncomp;i++){
    ip[ipp].fluxes[lcid][i]=fl[i];
  }
  
}

/**
   function returns components of flux
   
   @param lcid - load case id = number ot unknown
   @param ipp - number of integration point
   @param fl - array containing components of flux

   14/05/2010, TKr
*/
void coupmatu::givefluxes_cmu (long lcid,long ipp,vector &fl)
{
  long i,ncomp;
  ncomp=fl.n;
  
  for (i=0;i<ncomp;i++){
    fl[i] = ip[ipp].fluxes[lcid][i];
  }  
}


/**
   function reads integration points, material characteristics
   
   @param in - input stream
   
   21.7.2001
*/
void coupmatu::read (XFILE *in)
{
  //  computation of number of all integration points
  tnip = intpnum ();
  if (Mesprc==1){
    fprintf (stdout,"\n number of integration points in coupmatu %ld",tnip);
  }
  
  //  reading of integration points
  //readip (in);
  intpointalloc ();
  intpointinit ();
  
  //  reading of material characteristics
  readmatchar (in);
}

void coupmatu::readmatchar (XFILE *in)
{
  long i,j,k;
  long numtype;
  mattypec mattype;
  
  xfscanf (in,"%ld",&nmt);
  if (nmt<1){
    fprintf (stderr,"\n\n wrong number of material types in function coupmatu::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
  }
  if (Mesprc==1)  fprintf (stdout,"\n number of different types of materials  %ld",nmt);

  for (i=0;i<nmt;i++){
    xfscanf (in,"%d %ld",(int*)&mattype,&numtype);
    if (numtype<1){
      fprintf (stderr,"\n\n wrong number of material characteristics in function coupmatu::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
    }
    
    switch (mattype){//METR material type
     
    case isotransmatc:{
      if (Mesprc==1)  fprintf (stdout,"\n number of material models for concrete at high temperature %ld",numtype);
      itrmc = new isotrmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of material models for concrete at high temperature");
	  fprintf (stderr,"\n in function coupmatu::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	itrmc[k-1].read (in);
      }
      break;
    }
     
    case sejtkrc:{
      if (Mesprc==1)  fprintf (stdout,"\n number of material models for one-phase flow in deforming medium %ld",numtype);
      sejtkrmc = new sejtkrmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of material models for one-phase flow in deforming medium is required");
	  fprintf (stderr,"\n in function coupmatu::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	sejtkrmc[k-1].read (in);
      }
      break;
    }

    case consolawf1c:{
      if (Mesprc==1)  fprintf (stdout,"\n number of material models for one-phase flow in deforming medium %ld",numtype);
      consol_awf1c = new con_awf1matc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of material models for one-phase flow in deforming medium is required");
	  fprintf (stderr,"\n in function coupmatu::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	consol_awf1c[k-1].read (in);
      }
      break;
    }

    case consolwf1c:{
      if (Mesprc==1)  fprintf (stdout,"\n number of material models for one-phase flow in deforming medium %ld",numtype);
      consol_wf1c = new con_wf1matc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of material models for one-phase flow in deforming medium is required");
	  fprintf (stderr,"\n in function coupmatu::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	consol_wf1c[k-1].read (in);
      }
      break;
    }

    case concreteBc:{
      if (Mesprc==1)  fprintf (stdout,"\n number of material models for concrete at high temperature %ld",numtype);
      concretec = new concreteBmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of material models for concrete at high temperature");
	  fprintf (stderr,"\n in function coupmatu::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	concretec[k-1].read (in);
      }
      break;
    }

    case baroghelBc:{
      if (Mesprc==1)  fprintf (stdout,"\n number of baroghel material models for concrete %ld",numtype);
      baroghelc = new baroghelmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of baroghel material models for concrete");
	  fprintf (stderr,"\n in function coupmatu::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	baroghelc[k-1].read (in);
      }
      break;
    }
      

    case C60baroghelBc:{
      if (Mesprc==1)  fprintf (stdout,"\n number of baroghel material models for concrete %ld",numtype);
      C60baroghelc = new C60barmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of baroghel material models for concrete");
	  fprintf (stderr,"\n in function coupmatu::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	C60baroghelc[k-1].read (in);
      }
      break;
    }

    case C30baroghelBc:{
      if (Mesprc==1)  fprintf (stdout,"\n number of baroghel material models for concrete %ld",numtype);
      C30baroghelc = new C30barmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of baroghel material models for concrete");
	  fprintf (stderr,"\n in function coupmatu::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	C30baroghelc[k-1].read (in);
      }
      break;
    }


    case C60bazantBc:{
      if (Mesprc==1)  fprintf (stdout,"\n number of baroghel material models for concrete %ld",numtype);
      C60bazantc = new C60bazmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of baroghel material models for concrete");
	  fprintf (stderr,"\n in function coupmatu::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	C60bazantc[k-1].read (in);
      }
      break;
    }
      

    case o30bazantBc:{
      if (Mesprc==1)  fprintf (stdout,"\n number of baroghel material models for concrete %ld",numtype);
      o30bazantc = new o30bazmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of baroghel material models for concrete");
	  fprintf (stderr,"\n in function coupmatu::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	o30bazantc[k-1].read (in);
      }
      break;
    }
      
    case glasgowc:{
      if (Mesprc==1)  fprintf (stdout,"\n number of Glasgow material models %ld",numtype);
      tenchc = new glasgowmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of Glasgow material models for concrete");
	  fprintf (stderr,"\n in function coupmatu::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	tenchc[k-1].read (in);
      }
      break;
    }
      

    case consolawf2c:{
      if (Mesprc==1)  fprintf (stdout,"\n number of material models for two-phase flow in deforming medium %ld",numtype);
      consol_awf2c = new con_awf2matc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  print_err("wrong number of material models for two-phase flow in deforming medium is required",__FILE__,__LINE__,__func__);
	}
	consol_awf2c[k-1].read (in);
      }
      break;
    }
     

    case consolwf2c:{
      if (Mesprc==1)  fprintf (stdout,"\n number of material models for two-phase flow in deforming medium %ld",numtype);
      consol_wf2c = new con_wf2matc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  print_err("wrong number of material models for two-phase flow in deforming medium is required",__FILE__,__LINE__,__func__);
	}
	consol_wf2c[k-1].read (in);
      }
      break;
    }
     

    case consolhwf2c:{
      if (Mesprc==1)  fprintf (stdout,"\n number of material models for two-phase flow in deforming medium %ld",numtype);
      consol_hwf2c = new con_hwf2matc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  print_err("wrong number of material models for two-phase flow in deforming medium is required",__FILE__,__LINE__,__func__);
	}
	consol_hwf2c[k-1].read (in);
      }
      break;
    }
     
    case consolhawf3c:{
      if (Mesprc==1)  fprintf (stdout,"\n number of material models for three-phase flow in deforming medium %ld",numtype);
      consol_hawf3c = new con_hawf3matc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  print_err("wrong number of material models for two-phase flow in deforming medium is required",__FILE__,__LINE__,__func__);
	}
	consol_hawf3c[k-1].read (in);
      }
      break;
    }
     
    case glascoup:{
      if (Mesprc==1)  fprintf (stdout,"\n number of Glasgow material models %ld",numtype);
      gcm = new glasgowcoup [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of Glasgow material models for concrete");
	  fprintf (stderr,"\n in function coupmatu::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	gcm[k-1].read (in);
      }
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown material type is required in function coupmatu::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
    }
    }
    
  } 
}
