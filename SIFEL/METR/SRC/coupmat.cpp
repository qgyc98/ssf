#include "coupmat.h"
#include "globalc.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "constrelcl.h"
#include "onemediumc.h"
#include "twomediac.h"
#include "threemediac.h"

coupmat::coupmat (void)
{
  //  the number of integration points
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
  consol_awf2c = NULL;
  consol_hwf2c = NULL;
  consol_hawf3c = NULL;

  tenchc = NULL;
  
  //  fully coupled material model for bentnites
  bento = NULL;
}

coupmat::~coupmat (void)
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
  delete [] consol_awf2c;
  delete [] consol_hwf2c;
  delete [] consol_hawf3c;

  delete [] tenchc;

  delete [] bento;
}

/**
   The function returns the number of components of ipp's other array.
   Material type is given by the im parameter, which means index of material
   in the array tm of given integration point ipp.
   
   @param ipp - integration point pointer
   @param im - index of material
   
   @return The function returns total number of components of other array.  
   
   JK, 10.4.2019
*/
long coupmat::givencompother (long ipp,long /*im*/)
{
  long ncompo = 0;
  
  switch (ip[ipp].tm){
  case bentonite:{
    ncompo = bento[ip[ipp].idm].givencompeqother (ipp);
    break;
  }
  default:{
    print_err("unknown material type is required", __FILE__, __LINE__, __func__);
  }
  }
  
  return ncompo;
}

/**
   The function returns the number of components of ipp's eqother array.
   Material type is given by the im parameter, which means index of material
   in the array tm of given integration point ipp.
   
   @param ipp - integration point pointer
   @param im - index of material
   
   @return The function returns total number of components of other array.  
   
   JK, 10.4.2019
*/
long coupmat::givencompeqother (long ipp,long /*im*/)
{
  long ncompeqo = 0;
  
  switch (ip[ipp].tm){
  case bentonite:{
    ncompeqo = bento[ip[ipp].idm].givencompeqother (ipp);
    break;
  }
  default:{
    print_err("unknown material type is required", __FILE__, __LINE__, __func__);
  }
  }
  
  return ncompeqo;
}


/*
void coupmat::ipalloc (void)
{
  long i,j,ii,jj,ne,ntm,gdim,ippl,mnb,nip,mncomp;
  strastrestate mssst;
  
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
    mncomp=Ct->give_mncomp (i);

    for (ii=0;ii<mnb;ii++){
      //  type of strain/stress state
      mssst=Ct->give_mssst (i,ii);

      for (jj=0;jj<ntm;jj++){
	ippl=Ct->elements[i].ippl[jj][ii];
	nip=Ct->give_lower_nip (i,jj,ii);

	for (j=0;j<nip;j++){
	  ip[ippl].nmcomp=mncomp;
	  ip[ippl].mssst=mssst;
	  ip[ippl].stresses = new double [mncomp];
	  ip[ippl].strains = new double [mncomp];
	  
	  if ((mssst == planestrain) || (mssst == planestress))
	    ip[ippl].nmcomp=4;
	  
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
*/

/**
   function computes pointers on integration points
   
   offdiagonal submatrices are divided into blocks
   each block has row and column index
   number of rows is equal to number of transported matters
   number of columns is equal to number of blocks in mechanics

   25.2.2003
*/
 /*
long coupmat::intpnum (void)
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
 */

/**
   function allocates integration points
*/
void coupmat::intpointalloc ()
{
  ip = new intpointsc [tnip];
}

/**
   function initializes basic data in integration points
   JK, 10.4.2019
*/
void coupmat::intpointinit ()
{
  long i,j,nip,ipp,ncompdispl,ncompstr,ntm,ncompgrad;
  strastrestate ssst;
  
  //  element - integration point mapping
  elip = new long [tnip];

  for (i=0;i<Ct->ne;i++){
    //  number of components in array displacement
    ncompdispl = Ct->give_ncompdispl (i);
    //  number of components in arrays strains/stresses
    ncompstr = Ct->give_ncompstr (i);
    //  strain/stress state
    ssst = Ct->give_ssst (i,0);
    //  number of transported media
    ntm = Ct->give_ntm (i);
    //  number of components in arrays gradients/fluxes
    ncompgrad = Ct->give_ncompgrad (i);
    //  number of integration points on element
    nip = Ct->give_nip (i);
    //  number of the first integration point
    ipp = Ct->elements[i].ipp;
    
    for (j=0;j<nip;j++){
      if (elip != NULL)
	elip[ipp] = i;
      
      ip[ipp].ncompdispl = ncompdispl;
      ip[ipp].ncompstr = ncompstr;
      ip[ipp].ntm = ntm;
      ip[ipp].ncompgrad = ncompgrad;
      ip[ipp].ssst = ssst;
      ip[ipp].tm = Ct->elements[i].tm;
      ip[ipp].idm = Ct->elements[i].idm;
      
      ipp++;
      
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
    /*
void coupmat::matcond (matrix &d,long ipp,long ri,long ci)
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
    */

/**
   function computes capacity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ipp - number of integration point
   @param ri,ci - row and column indices
   
   17.07.2005
*/
     /*
void coupmat::matcap (matrix &d,long ipp,long ri,long ci)
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
     */


/**
   function computes correct fluxes from actual gradients
   @param lcid - load case id
   @param ipp - number of integration point
   
   4.12.2002
*/
      /*
void coupmat::computenlfluxes (matrix &d,long lcid,long ipp)//zkontrolovat slozky pro strainy
{
   
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

}
  */

/**
   function stores computed strains
   
   @param ipp - integration point pointer
   @param fi - first index
   @param eps - %vector containing strain components
   
   17.7.2001
*/
       /*
void coupmat::storestrain_cml (long ipp,long fi,vector &eps)
{
  long i,j,ncomp=eps.n;
  j=0;
  for (i=fi;i<fi+ncomp;i++){
    ip[ipp].strains[i]=eps[j];
    j++;
  }
}
       */

/**
   function stores components of gradients to integration points
   
   @param lcid - load case id = number of unknown
   @param ipp - number of integration point
   @param gr - array containing components of gradient
   
   4.12.2002
*/
	/*
void coupmat::storegrad_cml (long lcid,long ipp,vector &gr)
{
  long i,j,ncomp;
  ncomp=gr.n;
  
  j=0;
  for (i=0;i<ncomp;i++){
    ip[ipp].grad[lcid][i]=gr[j];
    j++;
  }
  
}
	*/


/**
   function returns components of gradients
   
   @param lcid - load case id = number of unknown
   @param ipp - number of integration point
   @param gr - array containing components of flux

   14/05/2010, TKr
*/
	 /*
void coupmat::givegrad_cml (long lcid,long ipp,vector &gr)
{
  long i,ncomp;
  ncomp=gr.n;
  
  for (i=0;i<ncomp;i++){
    gr[i] = ip[ipp].grad[lcid][i];
  }
}
	 */

/**
   function stores components of flux to integration points
   
   @param lcid - load case id
   @param ipp - number of integration point
   @param fl - array containing components of flux
   
   14/05/2010, TKr
*/
/*
void coupmat::storeflux_cml (long lcid,long ipp,vector &fl)
{
  long i,ncomp;
  ncomp=fl.n;
  
  for (i=0;i<ncomp;i++){
    ip[ipp].fluxes[lcid][i]=fl[i];
  }
  
}
*/

/**
   function returns components of flux
   
   @param lcid - load case id = number ot unknown
   @param ipp - number of integration point
   @param fl - array containing components of flux

   14/05/2010, TKr
*/
 /*
void coupmat::givefluxes_cml (long lcid,long ipp,vector &fl)
{
  long i,ncomp;
  ncomp=fl.n;
  
  for (i=0;i<ncomp;i++){
    fl[i] = ip[ipp].fluxes[lcid][i];
  }  
}
 */

/**
   function reads integration points, material characteristics
   
   @param in - input stream
   
   21.7.2001
*/
void coupmat::read (XFILE *in)
{
  //  computation of the number of all integration points
  tnip = intpnum_fc ();
  if (Mesprc==1){
    fprintf (stdout,"\n the number of integration points coupmat %ld",tnip);
  }
  
  //  reading of integration points
  //readip (in);
  intpointalloc ();
  intpointinit ();
  
  //  reading of material characteristics
  readmatchar (in);
}

void coupmat::readmatchar (XFILE *in)
{
  long i,j,k;
  long numtype;
  mattypec mattype;
  
  xfscanf (in,"%ld",&nmt);
  if (nmt<1){
    fprintf (stderr,"\n\n wrong number of material types in function coupmat::readmatchar (%s, line %d).\n",__FILE__,__LINE__);
  }
  if (Mesprc==1)  fprintf (stdout,"\n number of different types of materials  %ld",nmt);
  
  for (i=0;i<nmt;i++){
    xfscanf (in,"%d %ld",(int*)&mattype,&numtype);
    if (numtype<1){
      print_err("wrong number of material characteristics",__FILE__,__LINE__,__func__);
    }
    
    switch (mattype){//METR material type     
      
    case isotransmatc:{
      if (Mesprc==1)  fprintf (stdout,"\n number of material models for concrete at high temperature %ld",numtype);
      itrmc = new isotrmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
          print_err("wrong number of isotransmat materials",__FILE__,__LINE__,__func__);
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
          print_err("wrong number of material models for one-phase flow in deforming medium",__FILE__,__LINE__,__func__);
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
      
    case concreteBc:{
      if (Mesprc==1)  fprintf (stdout,"\n number of material models for concrete at high temperature %ld",numtype);
      concretec = new concreteBmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
          print_err("wrong number of material models for concrete at high temperature",__FILE__,__LINE__,__func__);
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
          print_err("wrong number of baroghel material models for concrete",__FILE__,__LINE__,__func__);
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
          print_err("wrong number of baroghel material models for concrete",__FILE__,__LINE__,__func__);
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
          print_err("wrong number of baroghel material models for concrete",__FILE__,__LINE__,__func__);
	}
	C30baroghelc[k-1].read (in);
      }
      break;
    }


    case o30bazantBc:{
      if (Mesprc==1)  fprintf (stdout,"\n number of bazant material models for concrete %ld",numtype);
      o30bazantc = new o30bazmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
          print_err("wrong number of bazant material model",__FILE__,__LINE__,__func__);
	}
	o30bazantc[k-1].read (in);
      }
      break;
    }
      

    case C60bazantBc:{
      if (Mesprc==1)  fprintf (stdout,"\n number of bazant material models for concrete %ld",numtype);
      C60bazantc = new C60bazmatc [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
          print_err("wrong number of bazant material model",__FILE__,__LINE__,__func__);
	}
	C60bazantc[k-1].read (in);
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
          print_err("wrong number of Glasgow material models for concrete",__FILE__,__LINE__,__func__);
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
      if (Mesprc==1)  fprintf (stdout,"\n number of Glasgow material models %ld",numtype);
      gcm = new glasgowcoup [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
          print_err("wrong number of Glasgow material models for concrete",__FILE__,__LINE__,__func__);
	}
	gcm[k-1].read (in);
      }
      break;
    }
      
    case bentonite:{
      if (Mesprc==1)  fprintf (stdout,"\n number of Masin-Schrefler material models %ld",numtype);
      bento = new bentonitemat [numtype];
      for (j=0;j<numtype;j++){
	k=numtype+1;
	xfscanf (in,"%ld",&k);
	if (k>numtype || k<1){
          print_err("wrong number of Masin-Schrefler material models for bentonite",__FILE__,__LINE__,__func__);
	}
	bento[k-1].read (in);
      }
      break;

    }

    default:{
      print_err("unknown material type is required",__FILE__,__LINE__,__func__);
    }
    }
  }
}


/**
  Function computes stiffness %matrix of the material
  in the required integration point.

  @param d - stiffness %matrix (output)
  @param ipp - number of integration point
   
  @return The function returns stiffness %matrix in the parametr d.

  31.3.2019, JK
*/
void coupmat::matstiff (matrix &d,long ipp)
{
  long i;
  
  switch (ip[ipp].tm){
  case bentonite:{
    i=ip[ipp].idm;
    bento[i].matstiff (d,ipp);
    break;
  }
    
  default:
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
    
  }
}

/**
   function assembles conductivity/permeability %matrix
   
   @param d - conductivity/permeability %matrix
   @param ipp - number of integration point
   
   2.4.2019, JK
 */
void coupmat::matcond (matrix &d,long ipp)
{
  long i;
  
  switch (ip[ipp].tm){
  case bentonite:{
    i=ip[ipp].idm;
    bento[i].matcond (d);
    break;
  }
  default:
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
    
  }
  
}

/**
   function computes capacity coefficient
   
   @param ipp - number of integration point
   
   3.4.2019, JK
 */
double coupmat::c_pp_coeff (long ipp)
{
  long i;
  double c;
  
  switch (ip[ipp].tm){
  case bentonite:{
    i=ip[ipp].idm;
    c = bento[i].c_pp_coeff (ip[ipp]);
    break;
  }
  default:
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  
  return c;
}

/**
   function computes capacity coefficient
   
   @param ipp - number of integration point
   
   3.4.2019, JK
 */
double coupmat::c_up_coeff (long ipp)
{
  long i;
  double c;
  
  switch (ip[ipp].tm){
  case bentonite:{
    i=ip[ipp].idm;
    c = bento[i].c_up_coeff (ip[ipp]);
    break;
  }
  default:
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  
  return c;
}

/**
   function computes capacity coefficient
   
   @param ipp - number of integration point
   
   3.4.2019, JK
 */
double coupmat::c_pu_coeff (long ipp)
{
  long i;
  double c;
  
  switch (ip[ipp].tm){
  case bentonite:{
    i=ip[ipp].idm;
    c = bento[i].c_pu_coeff (ip[ipp]);
    break;
  }
  default:
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  
  return c;
}

/**
   function computes the number of integration points
   in the whole problem
   
   2.4.2019, JK
*/
long coupmat::intpnum_fc (void)
{
  long i,n;
  
  n=0;
  for (i=0;i<Ct->ne;i++){
    Ct->elements[i].ipp=n;
    n+=Ct->give_nip (i);
  }
  return n;
}

/**
   Function initializes values of internal material variables
   at integration point ipp for given material im.
   
   @param lcid - load case id
   @param ipp  - integration point pointer
   @param im   - index of material type for given ip
   @param ido  - index of internal variables for given material in the ipp eqother array
   
   @return The function deos not return anything.
   
   Created by Tomas Koudelka, 11.4.2019
*/
void coupmat::initvalues (long ipp)
{
  switch (ip[ipp].tm){
  case bentonite:{
    /*
    if (Mp->strainstate == 0)
      {
        compute_ipstrains(lcid);
        Mp->strainstate = 1;
      }
    */
    bento[ip[ipp].idm].initval (ipp);
    break;
  }
  default:
    print_err("unknown material type is required", __FILE__,__LINE__,__func__);
  }
}



/**
  Function stores components of array other.
   
  @param ipp - integration point pointer
  @param fi - first index
  @param ncomp - number of required components
  @param comp - array containing components
   
  @return The function does not return anything.

  Created by TKo, 8.4.2019
*/
void coupmat::storeother (long ipp,long fi,long ncomp,double *comp)
{
  long i,j;
  j=0;
  for (i=fi;i<fi+ncomp;i++){
    ip[ipp].other[i]=comp[j];
    j++;
  }
}

/**
  The function stores all strains.

  @param ipp - integration point pointer
  @param eps - %vector containing strain components
  
  @return The function does not return anything.
  
  Created by JK, 11.4.2019
*/
void coupmat::storestrain (long ipp,vector &eps)
{
  long i,ncompstr;
  
  ncompstr = ip[ipp].ncompstr;
  if (ncompstr != eps.n)
    print_err("different number of strain components (in intpoint %ld, in the vector %ld", __FILE__, __LINE__, __func__,ncompstr,eps.n);
  
  for (i=0;i<ncompstr;i++){
    ip[ipp].strain[i]=eps[i];
  }
}
