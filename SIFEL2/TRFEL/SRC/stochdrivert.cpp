#include "stochdrivert.h"
#include "matrix.h"
#include "vector.h"
#include "globalt.h"
#include "isotrmat.h"
#include "crsection1d.h"
#include "crsection2d.h"
#include "crsection3d.h"



stochdrivert::stochdrivert (void)
{
  nsmt=0;  nscs=0;
  nsampl=0;  nstochvar=0;  nprunknowns=0;
  mt=NULL;  idm=NULL;  atm=NULL;
  cst=NULL;  idcs=NULL;  atcs=NULL;
  
  nna=NULL;  pnv=NULL;
}



stochdrivert::~stochdrivert (void)
{
  delete [] mt;  delete [] idm;  delete [] atm;
  delete [] cst;  delete [] idcs;  delete [] atcs;
  delete [] nna;  delete [] pnv;
}

/**
   The function reads input data form the opened text file.
   
   @param in - pointer to the opened XFILE
*/
void stochdrivert::read (XFILE *in)
{
  long i;
  XFILE *datin;
  
  xfscanf(in, " %a", auxfile);
  datin = xfopen (auxfile,"r");
  
  //  the number of stochastic materials
  xfscanf (in,"%ld",&nsmt);
  if (Mesprt==1)  fprintf (stdout,"\n number of stochastic materials  %ld",nsmt);
  if (nsmt>0){
    mt = new mattypet [nsmt];
    idm = new long [nsmt];
    atm = new atsel [nsmt];
    for (i=0;i<nsmt;i++){
      xfscanf (in,"%d %ld",(int*)mt+i,idm+i);
      idm[i]--;
      atm[i].read(in);
    }
  }
  
  //  the number of stochastic cross sections
  xfscanf (in,"%ld",&nscs);
  if (Mesprt==1)  fprintf (stdout,"\n number of stochastic cross sections  %ld",nscs);
  if (nscs>0){
    cst = new crsectypet [nscs];
    idcs = new long [nscs];
    atcs = new atsel [nscs];
    for (i=0;i<nscs;i++){
      xfscanf (in,"%d %ld",(int*)cst+i,idcs+i);
      idcs[i]--;
      atcs[i].read(in);
    }
  }
  
  //  the number of printed nodal values
  xfscanf (in,"%ld",&npnv);
  if (Mesprt==1)  fprintf (stdout,"\n number of printed nodal values  %ld",npnv);
  if (npnv>0){
    nna = new long [npnv];
    pnv = new atsel [npnv];
    for (i=0;i<npnv;i++){
      xfscanf (in,"%ld",nna+i);
      nna[i]--;
      pnv[i].read (in);
    }
  }
  
  readtable (datin);
  
  compute_nprunknowns ();
  
  xfclose (datin);
}



void stochdrivert::compute_nprunknowns ()
{
  long i;

  nprunknowns=0;
  for (i=0;i<npnv;i++){
    nprunknowns+=pnv[i].num;
  }
  
  allocm (nsampl,nprunknowns,stochtabout);
}



/**
  The function reads table of stochastic values.
  number of rows = number of samples
  number of columns = number of stochastic variables in the problem
   
  13.3.2003
*/
void stochdrivert::readtable (XFILE *in)
{
  xfscanf (in,"%ld %ld",&nsampl,&nstochvar);
  allocm (nsampl,nstochvar,stochtabin);
  readm(in, stochtabin);
}



/**
  The function writes table of stochastic values.
  number of rows = number of samples
  number of columns = number of stochastic variables in the problem
   
  13.3.2003
*/
void stochdrivert::writetable ()
{
  FILE *out;
  out = fopen (auxfile,"w");

  fprintf (out,"%ld\n",nprunknowns);
  printm(stochtabout,out,10,20);
  
  fclose (out);
}



/**
  The function changes stochastic values.
   
  13.3.2003
*/
void stochdrivert::changevalues (long sampleid)
{
  long i,j,k,l,m;
  vector val;
  
  j=0;
  for (i=0;i<nsmt;i++){
    k=atm[i].num;  m=0;
    allocv (k,val);
    for (l=j;l<j+k;l++){
      val[m]=stochtabin[sampleid][l];
      m++;
    }
    j+=k;
    changematerials (i,val);
    destrv (val);
  }

  for (i=0;i<nscs;i++){
    k=atcs[i].num;  m=0;
    allocv (k,val);
    for (l=j;l<j+k;l++){
      val[m]=stochtabin[sampleid][l];
      m++;
    }
    j+=k;
    changecrsections (i,val);
    destrv (val);
  }
}



/**
  The function imports stochastic values of the selected property 
  from the external source. It is used for interface between TRFEL and CTL library.

  @param dat - stochastic values of the parameters for the given property (material/cross-section)

  @return The function does not return anything, it changes parameters of the selected property stored 
          out of the stochdriver.

  12.2012
  
  Created by TKo+JK
*/
void stochdrivert::importvalues (double *dat)
{
  vector val;
  
  allocv (atm[0].num, val);
  copyv(dat, val);
  changematerials (0,val);
  destrv (val);
/*
  void importvalues (proptypet valtype, long id, double *dat);
  long i,j,k,l,m;
  vector val;
  switch (valtype)
  {
    case matelt:
      allocv (atm[id].num, val);
      copyv(dat, val);
      changematerials (id,val);
      destrv (val);
      break;
    case crosssect:
      allocv (atcs[id].num,val);
      copyv(dat, val);
      changecrsections (id,val);
      destrv (val);
      break;
    default:
      print_err("Unknown type of imported values is required", __FILE__, __LINE__, __func__);
  }
*/
}



void stochdrivert::changematerials (long id,vector &val)
{
  switch (mt[id]){
  case isotransmat:{
    Tm->itrm[idm[id]].changeparam (atm[id],val);
    break;
  }
  case damisotransmat:{
    Tm->damitrm[idm[id]].changeparam (atm[id],val);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
}



void stochdrivert::changecrsections (long id,vector &val)
{
  switch (cst[id]){
  case nocrosssectiont:{
    break;
  }
  case crsec1dt:{
    Tc->cs1d[idcs[id]].changeparam (atcs[id],val);
    break;
  }
  case crsec2dt:{
    Tc->cs2d[idcs[id]].changeparam (atcs[id],val);
    break;
  }
  case crsec3dt:{
    Tc->cs3d[idcs[id]].changeparam (atcs[id],val);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
}



void stochdrivert::extractor (long sampleid)
{
  long i,j,k,ci,nid;
  
  ci=0;
  for (i=0;i<npnv;i++){
    nid=nna[i];
    for (j=0;j<pnv[i].num;j++){
      k=Tt->give_dof (nid,pnv[i].atrib[j])-1;
      if (k<0){
	print_err("wrong number of DOF in function extractor",__FILE__,__LINE__,__func__);
      }
      stochtabout[sampleid][ci]=Lsrst->lhs[k];
      ci++;
    }
  }
}



/**
  The function exports required nodal values to the array val. It is used for
  interface between TRFEL and CTL library.

  @param - pointer to the allocated array where the required values will be exported

  @return The resulting values are stored in the array val.

  12.2012

  Created by TKo+JK  
*/
void stochdrivert::exportvalues (double *val)
{
  long i,j,k,ci,nid;
  
  ci=0;
  for (i=0;i<npnv;i++){
    nid=nna[i];
    for (j=0;j<pnv[i].num;j++){
      k=Tt->give_dof (nid,pnv[i].atrib[j])-1;
      if (k<0){
	print_err("Wrong number of DOF", __FILE__, __LINE__, __func__);
      }
      val[ci]=Lsrst->lhs[k];
      ci++;
    }
  }
}
