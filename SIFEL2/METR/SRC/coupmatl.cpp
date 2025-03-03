#include "coupmatl.h"
#include "globalc.h"
#include "mechtop.h"
#include "constrelcl.h"
#include "onemediumc.h"
#include "twomediac.h"
#include "threemediac.h"

coupmatl::coupmatl (void)
{
  tnip=0;
  ip=NULL;
  
  itrmc=NULL;
  sejtkrmc=NULL;

  concretec=NULL;
  baroghelc = NULL;
  C60baroghelc = NULL;
  C30baroghelc = NULL;
  o30bazantc = NULL;
  C60bazantc = NULL;

  consol_awf1c = NULL;
  consol_wf1c = NULL;
  consol_wf2c = NULL;
  consol_awf2c = NULL;
  consol_hwf2c = NULL;
  consol_hawf3c = NULL;

  tenchc = NULL;
}

coupmatl::~coupmatl (void)
{
  delete [] ip;

  delete [] itrmc;
  delete [] sejtkrmc;

  delete [] concretec;
  delete [] baroghelc;
  delete [] C60baroghelc;
  delete [] C30baroghelc;
  delete [] o30bazantc;
  delete [] C60bazantc;

  delete [] consol_awf1c;
  delete [] consol_wf1c;
  delete [] consol_wf2c;
  delete [] consol_awf2c;
  delete [] consol_hwf2c;
  delete [] consol_hawf3c;

  delete [] tenchc;  
}

/**
   function allocates auxiliary arrays on integration points
*/
/*
void coupmatl::ipalloc (void)
{
  long i,j;
  
  switch (Cp->tmatt){
  case mech_onemedium:{
    long ne,ii,jj,ippu,mnb,ntm,nip,ncompstr;
    strastrestate ssst;
    
    ne=Ct->ne;
    
    for (i=0;i<ne;i++){
      mnb=Ct->give_mnb (i);
      ntm = Tp->ntm;
      ncompstr=Ct->give_ncompstr (i);
      for (ii=0;ii<mnb;ii++){
	ssst=Ct->give_ssst (i,ii);
	for (jj=0;jj<ntm;jj++){
	  ippu=Ct->elements[i].ippu[ii][jj];
	  nip=Ct->give_upper_nip (i,ii,jj);
	  for (j=0;j<nip;j++){
	    ip[ippu].ncompstrstr=ncompstr;
	    ip[ippu].ssst=ssst;
	    ip[ippu].stresses = new double [ncompstr];
	    ip[ippu].strains = new double [ncompstr];
	    
	    if ((ssst == planestrain) || (ssst == planestress))
	      ip[ippu].ncompstrstr=4;
	    
	    ippu++;
	  }
	}
      }
    }

    for (i=0;i<tnip;i++){
      ip[i].av = new double [1];
      ip[i].grad = new double* [1];
      ip[i].fluxes = new double* [1];
      for (j=0;j<1;j++){
	ip[i].grad[j] = new double [Tp->gdim];
	ip[i].fluxes[j] = new double [Tp->gdim];
      }
    }
    break;
  }
   case mech_threemedia:{
    long ne,ii,jj,ippu,mnb,ntm,nip,ncompstr;
    strastrestate mssst;
    
    ne=Ct->ne;
    
    for (i=0;i<ne;i++){
      mnb=Ct->give_mnb (i);
      ntm = Tp->ntm;
      ncompstr=Ct->give_ncompstr (i);
      for (ii=0;ii<mnb;ii++){
	mssst=Ct->give_mssst (i,ii);
	for (jj=0;jj<ntm;jj++){
	  ippu=Ct->elements[i].ippu[ii][jj];
	  nip=Ct->give_upper_nip (i,ii,jj);
	  for (j=0;j<nip;j++){
	    ip[ippu].ncompstrstr=ncompstr;
	    ip[ippu].mssst=mssst;
	    ip[ippu].stresses = new double [ncompstr];
	    ip[ippu].strains = new double [ncompstr];
	    
	    if ((mssst == planestrain) || (mssst == planestress))
	      ip[ippu].ncompstrstr=4;
	    
	    ippu++;
	  }
	}
      }
    }

    for (i=0;i<tnip;i++){
      ip[i].av = new double [3];
      ip[i].grad = new double* [3];
      ip[i].fluxes = new double* [3];
      for (j=0;j<3;j++){
	ip[i].grad[j] = new double [Tp->gdim];
	ip[i].fluxes[j] = new double [Tp->gdim];
      }
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown problem type is required in function coupmatl::ipaloc (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
}
*/

void coupmatl::ipalloc (void)
{
  long i,j,ii,jj,ne,ntm,gdim,ippl,mnb,nip,ncompstr;
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
	ippl=Ct->elements[i].ippl[jj][ii];
	nip=Ct->give_lower_nip (i,jj,ii);

	for (j=0;j<nip;j++){
	  ip[ippl].ncompstr=ncompstr;
	  ip[ippl].ssst=ssst;
	  ip[ippl].stress = new double [ncompstr];
	  ip[ippl].strain = new double [ncompstr];
	  
	  if ((ssst == planestrain) || (ssst == planestress))
	    ip[ippl].ncompstr=4;
	  
	  ippl++;
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
   number of rows is equal to number of transported matters
   number of columns is equal to number of blocks in mechanics

   25.2.2003
*/
long coupmatl::intpnum (void)
{
  long i,j,k,ne,n,nb;
  elemtypec te;
  ne=Ct->ne;
  
  n=0;
  for (i=0;i<ne;i++){
    te = Ct->give_elem_type (i);
    nb=Mt->give_nb (i);
    Ct->elements[i].ippl = new long* [Tp->ntm];
    for (j=0;j<Tp->ntm;j++){
      Ct->elements[i].ippl[j] = new long [nb];
      for (k=0;k<nb;k++){
	Ct->elements[i].ippl[j][k]=n;
	n+=Ct->give_lower_nip (i,j,k);
      }
    }
  } 
  return n;
}

/**
   function allocates integration points
*/
void coupmatl::intpointalloc ()
{
  ip = new intpointsc [tnip];
  ipalloc ();
}

/**
   function initializes basic data in integration points
 */
void coupmatl::intpointinit ()
{
  long i,j,k,ii,jj,nb,ntm,nip,ipp;
  
  for (i=0;i<Tt->ne;i++){
    ntm=Tp->ntm;  k=0;
    nb=Mt->give_nb (i);
    for (ii=0;ii<ntm;ii++){
      for (jj=0;jj<nb;jj++){
	nip=Ct->give_lower_nip(i,ii,jj);
	ipp=Ct->elements[i].ippl[ii][jj];
	for (j=0;j<nip;j++){
	  ip[ipp].tm   = Ct->elements[i].tml[k];
	  ip[ipp].idm  = Ct->elements[i].idml[k];
	  ipp++;
	}
	k++;
      }
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
void coupmatl::matcond (matrix &d,long ipp,long ri,long ci)
{
  switch (Cp->tmatt){//  transported matter
  case mech_onemedium:{
    medc1 mc1;
    
    mc1.matcond_l(d,ri,ci,ipp);    
    break;    
  }
  case mech_twomedia:{
    medc2 mc2;
    
    mc2.matcond_l(d,ri,ci,ipp);    
    break;    
  }
  case mech_threemedia:{
    medc3 mc3;
    
    mc3.matcond_l(d,ri,ci,ipp);    
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
void coupmatl::matcap (matrix &d,long ipp,long ri,long ci)
{   
  switch (Cp->tmatt){//  transported matter
  case mech_onemedium:{
    medc1 mc1;

    mc1.matcap_l(d,ri,ci,ipp);
    break;    
  }
  case mech_twomedia:{
    medc2 mc2;

    mc2.matcap_l(d,ri,ci,ipp);
    break;    
  }
  case mech_threemedia:{
    medc3 mc3;

    mc3.matcap_l(d,ri,ci,ipp);
    break;    
  }
  default:{
    fprintf (stderr,"\n unknown number of transported media is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }  
}


/**
   function computes correct fluxes from actual gradients
   @param lcid - load case id
   @param ipp - number of integration point
   
   4.12.2002
*/
void coupmatl::computenlfluxes (matrix &d,long /*lcid*/,long ipp)//zkontrolovat slozky pro strainy
{
  long j,k,n,m;
  
  m = Tp->ntm; //number of transported media
  n = d.n;     //number of strain/stress components
  
  //null flux
  for (j=0;j<m;j++)
    for (k=0;k<Tp->gdim;k++) //geometric dimensions (1D, 2D or 3D)
      ip[ipp].fluxes[j][k] = 0.0;
  
  for (j=0;j<n;j++){//number of strain/stress components
    matcap (d,ipp,j,0);//ri = j; ci = 0
    for (k=0;k<m;k++){//number of transported media
      ip[ipp].fluxes[j][k] += d[0][j]*ip[ipp].strain[k];//the strain rate must be here
    }
  } 

  /*     
	 long j,k,m,n;
	 
	 m = d.m;
	 n = d.n;
	 
	 switch (Cp->tmatt){//  transported matter
	 case mech_onemedium:{
	 
	 switch (Cp->mednam){//  names of METR transported media
	 case mech_heat:{
	 
	 switch (lcid){//transported media
	 case 0:{
	 for (j=0;j<m;j++)//null flux
	 ip[ipp].fluxes[lcid][j] = 0.0;
	 
	 tu.matcap(d,ipp);
	 for (j=0;j<m;j++)
	 for (k=0;k<n;k++)
	 ip[ipp].fluxes[0][j] -= d[j][k]*ip[ipp].strains[k];//strain
	 
	 break;
	 }
	 default:{
	 fprintf (stderr,"\n unknown lcid is required in function (%s, line %d).\n",__FILE__,__LINE__);
	 }
	 }	
	 
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
	 
	 switch (lcid){//transported media
	 case 0:{
	 for (j=0;j<m;j++)//null flux
	 ip[ipp].fluxes[lcid][j] = 0.0;
	 
	 cu.matcap(d,ipp);
	 for (j=0;j<m;j++)
	 for (k=0;k<n;k++)
	 ip[ipp].fluxes[0][j] -= d[j][k]*ip[ipp].strains[k];//strain
	 
	 break;
	 }
	 case 1:{
	 for (j=0;j<m;j++)//null flux
	 ip[ipp].fluxes[lcid][j] = 0.0;
	 
	 gu.matcap(d,ipp);
	 for (j=0;j<m;j++)
	 for (k=0;k<n;k++)
	 ip[ipp].fluxes[1][j] -= d[j][k]*ip[ipp].strains[k];//strain
	 
	 break;
	 }
	 case 2:{
	 for (j=0;j<m;j++)//null flux
	 ip[ipp].fluxes[lcid][j] = 0.0;
	 
	 tu.matcap(d,ipp);
	 for (j=0;j<m;j++)
	 for (k=0;k<n;k++)
	 ip[ipp].fluxes[2][j] -= d[j][k]*ip[ipp].strains[k];//strain
	 break;
	 }
	 default:{
	 fprintf (stderr,"\n unknown lcid is required in function (%s, line %d).\n",__FILE__,__LINE__);
	 }
	 }	
	 
	 break;
	 }
	 default:{
	 fprintf (stderr,"\n unknown name of transported media is required in function (%s, line %d).\n",__FILE__,__LINE__);
	 }
	 }
	 
	 break;    
	 }
	 default:{
	 fprintf (stderr,"\n\n unknown material is required in function computenlfluxes (%s, line %d).\n",__FILE__,__LINE__);
	 abort ();
	 }
	 }
  */
}

/**
   function stores computed strains
   
   @param ipp - integration point pointer
   @param fi - first index
   @param eps - %vector containing strain components
   
   17.7.2001
*/
void coupmatl::storestrain_cml (long ipp,long fi,vector &eps)
{
  long i,j,ncomp=eps.n;
  j=0;
  for (i=fi;i<fi+ncomp;i++){
    ip[ipp].strain[i]=eps[j];
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
void coupmatl::storegrad_cml (long lcid,long ipp,vector &gr)
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
   
   @param lcid - load case id = number of unknown
   @param ipp - number of integration point
   @param gr - array containing components of flux

   14/05/2010, TKr
*/
void coupmatl::givegrad_cml (long lcid,long ipp,vector &gr)
{
  long i,ncomp;
  ncomp=gr.n;
  
  for (i=0;i<ncomp;i++){
    gr[i] = ip[ipp].grad[lcid][i];
  }
}

/**
   function stores components of flux to integration points
   
   @param lcid - load case id
   @param ipp - number of integration point
   @param fl - array containing components of flux
   
   14/05/2010, TKr
*/
void coupmatl::storeflux_cml (long lcid,long ipp,vector &fl)
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
void coupmatl::givefluxes_cml (long lcid,long ipp,vector &fl)
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
void coupmatl::read (XFILE *in)
{
  //  computation of number of all integration points
  tnip = intpnum ();
  if (Mesprc==1){
    fprintf (stdout,"\n number of integration points coupmatl %ld",tnip);
  }
  
  //  reading of integration points
  //readip (in);
  intpointalloc ();
  intpointinit ();
  
  //  reading of material characteristics
  readmatchar (in);
}

void coupmatl::readmatchar (XFILE *in)
{
  long i,j,k;
  long numtype;
  mattypec mattype;
  
  xfscanf (in,"%ld",&nmt);
  if (nmt<1){
    fprintf (stderr,"\n\n wrong number of material types in function coupmatl::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
  }
  if (Mesprc==1)  fprintf (stdout,"\n number of different types of materials  %ld",nmt);
  
  for (i=0;i<nmt;i++){
    xfscanf (in,"%d %ld",(int*)&mattype,&numtype);
    if (numtype<1){
      fprintf (stderr,"\n\n wrong number of material characteristics in function coupmatl::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
    }
    
    switch (mattype){//METR material type     
      
    case isotransmatc:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for concrete at high temperature %ld",numtype);
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
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for one-phase flow in deforming medium %ld",numtype);
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
      if (Mesprc==1)  fprintf (stdout,"\n number of material models for two-phase flow in deforming medium %ld",numtype);
      consol_awf1c = new con_awf1matc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  print_err("wrong number of material models for two-phase flow in deforming medium is required",__FILE__,__LINE__,__func__);
	}
	consol_awf1c[k-1].read (in);
      }
      break;
    }

    case consolwf1c:{
      if (Mesprc==1)  fprintf (stdout,"\n number of material models for two-phase flow in deforming medium %ld",numtype);
      consol_wf1c = new con_wf1matc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  print_err("wrong number of material models for one-phase flow in deforming medium is required",__FILE__,__LINE__,__func__);
	}
	consol_wf1c[k-1].read (in);
      }
      break;
    }
      
   case concreteBc:{
      if (Mesprt==1)  fprintf (stdout,"\n number of material models for concrete at high temperature %ld",numtype);
      concretec = new concreteBmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of material models for concrete at high temperature");
	  fprintf (stderr,"\n in function coupmatl::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	concretec[k-1].read (in);
      }
      break;
    }

    case baroghelBc:{
      if (Mesprt==1)  fprintf (stdout,"\n number of baroghel material models for concrete %ld",numtype);
      baroghelc = new baroghelmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of baroghel material models for concrete");
	  fprintf (stderr,"\n in function coupmatl::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	baroghelc[k-1].read (in);
      }
      break;
    }
      

    case C60baroghelBc:{
      if (Mesprt==1)  fprintf (stdout,"\n number of baroghel material models for concrete %ld",numtype);
      C60baroghelc = new C60barmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of baroghel material models for concrete");
	  fprintf (stderr,"\n in function coupmatl::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	C60baroghelc[k-1].read (in);
      }
      break;
    }

    case C30baroghelBc:{
      if (Mesprt==1)  fprintf (stdout,"\n number of baroghel material models for concrete %ld",numtype);
      C30baroghelc = new C30barmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of baroghel material models for concrete");
	  fprintf (stderr,"\n in function coupmatl::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	C30baroghelc[k-1].read (in);
      }
      break;
    }


    case o30bazantBc:{
      if (Mesprt==1)  fprintf (stdout,"\n number of baroghel material models for concrete %ld",numtype);
      o30bazantc = new o30bazmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of baroghel material models for concrete");
	  fprintf (stderr,"\n in function coupmatl::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	o30bazantc[k-1].read (in);
      }
      break;
    }
      

    case C60bazantBc:{
      if (Mesprt==1)  fprintf (stdout,"\n number of baroghel material models for concrete %ld",numtype);
      C60bazantc = new C60bazmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of baroghel material models for concrete");
	  fprintf (stderr,"\n in function coupmatl::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	C60bazantc[k-1].read (in);
      }
      break;
    }
      
    case glasgowc:{
      if (Mesprt==1)  fprintf (stdout,"\n number of Glasgow material models %ld",numtype);
      tenchc = new glasgowmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of Glasgow material models for concrete");
	  fprintf (stderr,"\n in function coupmatl::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
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
	  print_err("wrong number of material models for three-phase flow in deforming medium is required",__FILE__,__LINE__,__func__);
	}
	consol_hawf3c[k-1].read (in);
      }
      break;
    }
      
    case glascoup:{
      if (Mesprt==1)  fprintf (stdout,"\n number of Glasgow material models %ld",numtype);
      gcm = new glasgowcoup [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
	  fprintf (stderr,"\n\n wrong number of Glasgow material models for concrete");
	  fprintf (stderr,"\n in function coupmatl::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
	}
	gcm[k-1].read (in);
      }
      break;
    }

    default:{
      fprintf (stderr,"\n\n unknown material type is required in function coupmatl::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
}
