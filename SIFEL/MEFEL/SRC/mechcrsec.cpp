#include "mechcrsec.h"
#include "global.h"
#include "mechmat.h"
#include "mechtop.h"
#include "matrix.h"
#include "vector.h"
#include "crsec2dbar.h"
#include "crsecplstr.h"
#include "crsec2dbeam.h"
#include "crsec3dbeam.h"
#include "crsec3d.h"
#include "crsecnod.h"
#include "crseclayer.h"
#include "element.h"
#include "node.h"



/**
  Constructor initializes data members to zero or default values.

  Created by JK,
*/
mechcrsec::mechcrsec ()
{
  //  number of cross section types
  ncst=0;
  //  type of cross sections
  cstype = NULL;
  //  number of instances of particular cross sections
  numtype = NULL;

  cs2dbar=NULL;
  cs2dbeam=NULL;
  cs3dbeam=NULL;
  csplstr=NULL;
  cs3d=NULL;
  csnod=NULL;
  cslay=NULL;
}



/**
  Destructor releases allocated memory of the mechcrsec object.

  Created by JK,
*/
mechcrsec::~mechcrsec ()
{
  delete [] cstype;
  delete [] numtype;
  delete [] cs2dbar;
  delete [] cs2dbeam;
  delete [] cs3dbeam;
  delete [] csplstr;
  delete [] cs3d;
  delete [] csnod;
  delete [] cslay;
}



/**
  Function reads cross section characteristics.

  @param in - pointer to the opened text XFILE

  @return The function does not return anything.

  Created by JK, 22.7.2001
  Modified by TKo 26.6.2014
*/
void mechcrsec::read (XFILE *in)
{
  long i;

  xfscanf (in, "%k%ld", "num_crsec_types", &ncst);
  if (ncst<0)
    print_err("negative number of cross-sections", __FILE__, __LINE__, __func__);
  
  //  types of cross sections
  cstype = new crsectype [ncst];
  //  number of instances of particular cross section types
  numtype = new long [ncst];


  for (i=0;i<ncst;i++){
    xfscanf (in, "%k%m %k%ld", "crstype", &crsectype_kwdset, &cstype[i], "num_inst", &numtype[i]);
    if (numtype[i]<1)
      print_err("negative number of cross-section", __FILE__, __LINE__, __func__);

    readcrsectype(in, cstype[i], numtype[i]);
  }
}



/**
  The function reads numtype of cross section parameters for the given type of cross 
  section from the opened text file.

  @param in - pointer to the opened text XFILE
  @param ct - type of cross section read
  @param numt - number of parameter sets that will be read

  @return The function does not return anything.

  Created by Tomas Koudelka according to old version of function read 26.6.2014
*/
void mechcrsec::readcrsectype(XFILE *in, crsectype ct, long numt)
{
  long j, k;
  switch (ct)
  {
    case nocrosssection:{
      break;
    }
    case csbar2d:{
      cs2dbar = new crsec2dbar [numt];
      for (j=0;j<numt;j++){
	xfscanf (in,"%ld",&k);
	if ((k<1) || (k>numt))
        {
          print_err("index %ld of cross section parameter set is out of range <1,%ld>", 
                    __FILE__, __LINE__, __func__, k, numt);
        }
	cs2dbar[k-1].read (in);
      }
      break;
    }
    case csbeam2d:{
      cs2dbeam = new crsec2dbeam [numt];
      for (j=0;j<numt;j++){
	xfscanf (in,"%ld",&k);
	if ((k<1) || (k>numt))
        {
          print_err("index %ld of cross section parameter set is out of range <1,%ld>", 
                    __FILE__, __LINE__, __func__, k, numt);
        }
	cs2dbeam[k-1].read (in);
      }
      break;
    }
    case csbeam3d:{
      cs3dbeam = new crsec3dbeam [numt];
      for (j=0;j<numt;j++){
	xfscanf (in,"%ld",&k);
	if ((k<1) || (k>numt))
        {
          print_err("index %ld of cross section parameter set is out of range <1,%ld>", 
                    __FILE__, __LINE__, __func__, k, numt);
        }
	cs3dbeam[k-1].read (in);
      }
      break;
    }
    case csplanestr:{
      csplstr = new crsecplstr [numt];
      for (j=0;j<numt;j++){
	xfscanf (in,"%ld",&k);       
	if ((k<1) || (k>numt))
        {
          print_err("index %ld of cross section parameter set is out of range <1,%ld>", 
                    __FILE__, __LINE__, __func__, k, numt);
        }
	csplstr[k-1].read (in);
      }
      break;
    }
    case cs3dprob:{
      cs3d = new crsec3d [numt];
      for (j=0;j<numt;j++){
	xfscanf (in,"%ld",&k);
	if ((k<1) || (k>numt))
        {
          print_err("index %ld of cross section parameter set is out of range <1,%ld>", 
                    __FILE__, __LINE__, __func__, k, numt);
        }
	cs3d[k-1].read (in);
      }
      break;
    }
    case csnodal:{
      csnod = new crsecnod [numt];
      for (j=0;j<numt;j++){
	xfscanf (in,"%ld",&k);
	if ((k<1) || (k>numt))
        {
          print_err("index %ld of cross section parameter set is out of range <1,%ld>", 
                    __FILE__, __LINE__, __func__, k, numt);
        }
	csnod[k-1].read (in);
      }
      break;
    }
    case cslayer:
      cslay = new crseclayer [numt];
      for (j=0;j<numt;j++){
        xfscanf (in,"%ld",&k);
        if (k<1){}
        if (k>numt){}
        cslay[k-1].read (in);
      }
      break;
    default:
      print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
  }
}



/**
  Function prints cross section characteristics.

  @param out - pointer to the opened text FILE

  @return The function does not return anything.

  Created by TKr, 02/03/2013
  Modified by TKo, 1.7.2014
*/
void mechcrsec::print (FILE *out)
{
  long i,j;

  fprintf (out,"\n\n ## cross-sections:");
  fprintf (out,"\n%ld\n",ncst);
  
  for (i=0;i<ncst;i++){
    fprintf (out,"\n%d %ld",(int)cstype[i],numtype[i]);

    for (j=0;j<numtype[i];j++)
    {
      fprintf (out,"\n%ld ",j+1);
      printcrschar(out, cstype[i], j);
    }
  }
}



/**
  Function prints cross section characteristics for the given cross section type and index.

  @param out - pointer to the opened text FILE
  @param ct  - required cross section type
  @param numinst - index of required cross section pramater set

  @return The function does not return anything.

  According to original function print from TKr -  TKo, 1.7.2014
*/
void mechcrsec::printcrschar (FILE *out, crsectype ct, long numinst)
{
  switch (ct)
  {
    case nocrosssection:
      break;

    case csbar2d:{
      cs2dbar[numinst].print (out);
      break;
    }

    case csbeam2d:{
      cs2dbeam[numinst].print (out);
      break;
    }

    case csbeam3d:{
      cs3dbeam[numinst].print (out);
      break;
    }

    case csplanestr:{
      csplstr[numinst].print (out);
      break;
    }

    case cs3dprob:{
      cs3d[numinst].print (out);
      break;
    }

    case csnodal:{
      csnod[numinst].print (out);
      break;
    }

    default:
      print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
  }
}



/**
  Function returns thicknesses in nodes.
   
  @param nod - vector containing numbers of nodes defining element
  @param t - array containing thicknesses (output)
   
  @return The function returns thickness at the given node in the parameter t

  Created by JK, 22.7.2001
*/
void mechcrsec::give_thickn (ivector &nod,vector &t)
{
  long i,n;
  n=nod.n;
  for (i=0;i<n;i++){
    switch (Mt->nodes[nod[i]].crst){
    case nocrosssection:{
      break;
    }
    case csbar2d:{
      break;
    }
    case csbeam2d:{
      break;
    }
    case csbeam3d:{
      break;
    }
    case csplanestr:{
      t[i]=csplstr[Mt->nodes[nod[i]].idcs].t;
      break;
    }
    case cslayer:{
      //      t[i]=cslay[Mt->nodes[nod[i]].idcs].th;
      t[i]=1.0; // artificial thickness - thickness integration on plate elements is performed inside the layer model
      break;
    }
    default:{
      print_err("unknown cross section type is required", __FILE__, __LINE__, __func__);
    }
    }
  }
}



/**
  Function returns thickness on element.
   
  @param eid - element id
  @param t - thickness (output)
  
  @return The function returns thickness of the given element in the parameter t.
  
  Created by JK, 22.7.2001
*/
void mechcrsec::give_thicke (long eid,double &t)
{
  switch (Mt->elements[eid].crst){
  case nocrosssection:{
    break;
  }
  case csbar2d:{
    break;
  }
  case csbeam2d:{
    break;
  }
  case csbeam3d:{
    break;
  }
  case csplanestr:{
    t=csplstr[Mt->elements[eid].idcs].t;
    break;
  }
  case cslayer:{
    //    t=cslay[Mt->elements[eid].idcs].th;
    t=1.0; // artificial thickness - thickness integration on plate elements is performed inside the layer model
    break;
  }
  default:
    print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
    abort();
  }
}



/**
  Function returns thickness on element.
   
  @param eid - element id
  @param nodes - %vector of node numbers of the given element
  @param th - %vector of thicknesses at nodes (output)
  
  @return The function returns thickness at nodes of the given element in the parameter th.
  
  Created by JK, 22.7.2001
*/
void mechcrsec::give_thickness (long eid,ivector &nodes,vector &th)
{
  long i,nne;
  double t;
  
  if (Mt->elements[eid].crst==0){
    give_thickn (nodes,th);
  }
  else{
    give_thicke (eid,t);
    nne = Mt->give_nne (eid);
    for (i=0;i<nne;i++){
      th[i]=t;
    }
  }
}



/**
  Function returns thickness of the given crossection type and id.
   
  @param crst - cross-section type
  @param idcst - cross-sectio type id
  
  @return The function returns thickness of the given cross-section type and id.
  
  Created by JK, 22.7.2001
*/
double mechcrsec::give_onethickness (crsectype crst,long idcs)
{
  double t=0.0;
  
  switch (crst){
  case nocrosssection:{
    break;
  }
  case csbar2d:{
    break;
  }
  case csbeam2d:{
    t=cs2dbeam[idcs].t;
    break;
  }
  case csbeam3d:{
    break;
  }
  case csplanestr:{
    t=csplstr[idcs].t;
    break;
  }
  case csnodal:{
    t=csnod[idcs].t;
    break;
  }
  case cslayer:{
    t=cslay[idcs].th;
    break;
  }
  default:{
    print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
    abort();
  }
  }
  
  return t;
}



/**
  Function returns area of cross-section.
   
  @param eid - element id
  @param a - cross-sectio area of the element (output)

  @return The function returns cross-section area of the given element in the parameter a.

  Created by JK, 10.8.2001
*/
void mechcrsec::give_area (long eid,double &a)
{
  switch (Mt->elements[eid].crst){
  case nocrosssection:{
    break;
  }
  case csbar2d:{
    a=cs2dbar[Mt->elements[eid].idcs].a;
    break;
  }
  case csbeam2d:{
    a=cs2dbeam[Mt->elements[eid].idcs].a;
    break;
  }
  case csbeam3d:{
    a=cs3dbeam[Mt->elements[eid].idcs].a;
    break;
  }
  default:{
    print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
  Function returns moment(s) of inertia of cross-section.
   
  @param eid - element id
  @param i - array containing moments of inertia (output)
   
  @return The function returns cross-section moments of inertia in the parameter i.

  Created by JK, 10.8.2001
*/
void mechcrsec::give_mominer (long eid,double *i)
{
  switch (Mt->elements[eid].crst){
  case nocrosssection:{
    break;
  }
  case csbeam2d:{
    //*i=cs2dbeam[Mt->elements[eid].idcs].iy;
    cs2dbeam[Mt->elements[eid].idcs].give_moments (i);
    break;
  }
  case csbeam3d:{
    cs3dbeam[Mt->elements[eid].idcs].give_moments (i);
    break;
  }
  default:{
    print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
  Function returns shear coefficient(s) of cross-section
  The shear coefficient for rectangle is 5/6.
   
  @param eid - element id
  @param shearcoeff - array containing shear coefficients (output)
   
  @return The function returns shear coefficients in the parameter shearcoeff.

  Created by JK, 10.8.2001
*/
void mechcrsec::give_shearcoeff (long eid,double *shearcoeff)
{
  switch (Mt->elements[eid].crst){
  case nocrosssection:{
    break;
  }
  case csbeam2d:{
    //*shearcoeff=cs2dbeam[Mt->elements[eid].idcs].shearcoeff;
    cs2dbeam[Mt->elements[eid].idcs].give_shearcoeff (shearcoeff);
    break;
  }
  case csbeam3d:{
    cs3dbeam[Mt->elements[eid].idcs].give_shearcoeff (shearcoeff);
    break;
  }
  default:{
    print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
  Function returns shear %vector axis [z] cross-section.
  Shear coefficient for rectangle is 5/6.
   
  @param eid - element id
  @param vec - shear %vector (output)
   
  @return The function returns shear %vector axis in the parameter vec.

  Created by JK, 10.8.2001
*/
void mechcrsec::give_vectorlcs (long eid,vector &vec)
{
  switch (Mt->elements[eid].crst){
  case csbeam3d:{
    copyv(cs3dbeam[Mt->elements[eid].idcs].lcs,vec);
    break;
  }
  default:{
    print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
  Fnction returns densities in nodes.
   
  @param nod - vector containing numbers of nodes defining element
  @param rho - array containing densities (output)
   
  @return The function returns nodal densities in the parameter rho.

  Created by JK, 22.7.2001
*/
void mechcrsec::give_densityn (ivector &nod,vector &rho)
{
  long i,n;
  n=nod.n;
  for (i=0;i<n;i++){
    switch (Mt->nodes[nod[i]].crst){
    case nocrosssection:{
      break;
    }
    case csbar2d:{
      break;
    }
    case csbeam2d:{
      break;
    }
    case csbeam3d:{
      break;
    }
    case csplanestr:{
      rho[i]=csplstr[Mt->nodes[nod[i]].idcs].rho;
      break;
    }
    case cs3dprob:{
      rho[i]=cs3d[Mt->nodes[nod[i]].idcs].rho;
      break;
    }
    case cslayer:{
      rho[i]=cslay[Mt->nodes[nod[i]].idcs].rho;
      break;
    }
    default:{
      print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
    }
    }
  }
}



/**
  Function returns density on element.
   
  @param eid - element id
  @param rho - density (output)
   
  @return The function returns density in the parameter rho.

  Created by JK 10.8.2001
*/
void mechcrsec::give_densitye (long eid,double &rho)
{
  switch (Mt->elements[eid].crst){
  case nocrosssection:{    
    break;
  }
  case csbar2d:{
    rho=cs2dbar[Mt->elements[eid].idcs].rho;
    break;
  }
  case csbeam2d:{
    rho=cs2dbeam[Mt->elements[eid].idcs].rho;
    break;
  }
  case csbeam3d:{
    rho=cs3dbeam[Mt->elements[eid].idcs].rho;
    break;
  }
  case csplanestr:{
    rho=csplstr[Mt->elements[eid].idcs].rho;
    break;
  }
  case cs3dprob:{
    rho=cs3d[Mt->elements[eid].idcs].rho;
    break;
  }
  default:{
    print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
    abort();
  }
  }
}



/**
  Function returns density on element.
   
  @param eid - element id
  @param nodes - %vector containing numbers of nodes defining element
  @param dens - %vector of densities in the given nodes (output)
   
  @return The function returns nodal densities in the parameter dens.

  Created by JK 10.8.2001
*/
void mechcrsec::give_density (long eid,ivector &nodes,vector &dens)
{
  long i,nne;
  double d = 0.0;
  
  if (Mt->elements[eid].crst==0){
    give_densityn (nodes,dens);
  }
  else{
    give_densitye (eid,d);
    nne = Mt->give_nne (eid);
    for (i=0;i<nne;i++){
      dens[i]=d;
    }
  }
}



/**
  Function returns weight of concentrated mass
   
  @param cr - type of cross section
  @param idcs - number of cross section
   
  @return The function returns weight of concentrated mass for the given 
          cross-section type and id.

  Created by JK, 25.7.2005
*/
double mechcrsec::give_weight (crsectype cr,long idcs)
{
  double m = 0.0;
  
  switch (cr){
  case nocrosssection:{
    break;
  }
  case csbar2d:{
    m=cs2dbar[idcs].m;
    break;
  }
  case csbeam2d:{
    break;
  }
  case csbeam3d:{
    break;
  }
  case csplanestr:{
    m=csplstr[idcs].m;
    break;
  }
  case csnodal:{
    m=csnod[idcs].m;
    break;
  }
  case cslayer:{
    m=cslay[idcs].m;
    break;
  }
  default:{
    print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return m;
}


/**
  Function returns array of layer thicknesses on element.
   
  @param eid - element id
  @param *t - thickness (output)
  
  @return
  
  J.Fiedler, 10/2014
*/
double* mechcrsec::give_layer_thicke (long eid)
{
  double *t=NULL;
  
  switch (Mt->elements[eid].crst){
    case cslayer:{
      t = cslay[Mt->elements[eid].idcs].layth;
      break;
    }
    default:
      print_err("unknown type of layered cross section is required", __FILE__, __LINE__, __func__);
  }
  return t;
}


/**
  Function returns array of layer z-coordinates on element.
   
  @param eid - element id
  @param *z - z-coordinates (output)
  
  @return
  
  J.Fiedler, 10/2014
*/
double* mechcrsec::give_layer_zcoord (long eid)
{
  double *z=NULL;
  
  switch (Mt->elements[eid].crst){
    case cslayer:{
      z = cslay[Mt->elements[eid].idcs].layz;
      break;
    }
    default:
      print_err("unknown type of layered cross section is required", __FILE__, __LINE__, __func__);
  }
  return z;
}


/**
  Function returns number of layers on element.
   
  @param eid - element id
  
  @return
  
  J.Fiedler, 10/2014
*/
long mechcrsec::give_num_lay (long eid)
{
  long nl = 0.0;
  
  switch (Mt->elements[eid].crst){
    case cslayer:{
      nl = cslay[Mt->elements[eid].idcs].nl;
      break;
    }
    default:
      print_err("unknown type of layered cross section is required", __FILE__, __LINE__, __func__);
      abort();
  }
  return nl;
}
