#include "transcrsec.h"
#include "globalt.h"

transcrsec::transcrsec (void)
{
  //  number of cross section types
  ncst=0;
  //  type of cross sections
  cstype = NULL;
  //  number of instances of particular cross sections
  numtype = NULL;

  cs1d = NULL;
  cs2d = NULL;
  cs3d = NULL;
}
transcrsec::~transcrsec (void)
{
  delete [] cstype;
  delete [] numtype;
  delete [] cs1d;
  delete [] cs2d;
  delete [] cs3d;
}

/**
   function reads cross section characteristics

   @param in - input stream

   22.7.2001
*/
void transcrsec::read (XFILE *in)
{
  long i;

  xfscanf (in,"%k%ld","number_csect_types",&ncst);
  if (ncst<0){
    print_err("number of cross section types is less than 0 (line=%ld)",__FILE__,__LINE__,__func__, in->line);
    abort();
  }
  if (Mesprt==1){
    fprintf (stdout,"\n number of cross section types  %ld",ncst);
  }
  
  //  types of cross sections
  cstype = new crsectypet [ncst];
  //  number of instances of particular cross section types
  numtype = new long [ncst];
  
  for (i=0;i<ncst;i++){
    xfscanf (in, "%k%m %k%ld", "crstype", &crsectypet_kwdset, &cstype[i], "num_inst", &numtype[i]);
    if (numtype[i]<1){
      print_err("number of instances of cross section is less than 1",__FILE__,__LINE__,__func__);
    }
    
    readcrsectype(in, cstype[i], numtype[i]);
  }
}



/**
  The function reads numtype of cross section parameters for the given type of cross 
  section from the opened text file.

  @param in - pointer to the opened text XFILE
  @param cstype - type of cross section read
  @param numt - number of parameter sets that will be read

  @return The function does not return anything.

  Created by Tomas Koudelka according to old version of function read 26.6.2014
*/
void transcrsec::readcrsectype(XFILE *in, crsectypet ct, long numt)
{
  long j,k;

  switch (ct)
  {
    case nocrosssectiont:{
      break;
    }
    case crsec1dt:{
      if (Mesprt==1)  fprintf (stdout,"\n number of 1D cross sections  %ld",numt);
      cs1d = new crsection1d [numt];
      for (j=0;j<numt;j++){
	xfscanf (in,"%ld",&k);
	if ((k<1) || (k>numt))
        {
          print_err("index %ld of cross section parameter set is out of range <1,%ld>", 
                    __FILE__, __LINE__, __func__, k, numt);
        }
	cs1d[k-1].read (in);
      }
      break;
    }
    case crsec2dt:{
      if (Mesprt==1)  fprintf (stdout,"\n number of 2D cross sections  %ld",numt);
      cs2d = new crsection2d [numt];
      for (j=0;j<numt;j++){
	xfscanf (in,"%ld",&k);
	if ((k<1) || (k>numt))
        {
          print_err("index %ld of cross section parameter set is out of range <1,%ld>", 
                    __FILE__, __LINE__, __func__, k, numt);
        }
	cs2d[k-1].read (in);
      }
      break;
    }
    case crsec3dt:{
      if (Mesprt==1)  fprintf (stdout,"\n number of 3D cross sections  %ld",numt);
      cs3d = new crsection3d [numt];
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
    default:{
      print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
    }
  }
}



/**
   function prints cross section characteristics

   @param out - output stream

   TKr 3.1.2006
*/
void transcrsec::print (FILE *out)
{
  long i,j;

  fprintf(out,"\n\n## cross-sections:");
  fprintf (out,"\n%ld\n",ncst);

  for (i=0;i<ncst;i++){
    fprintf (out,"\n%d %ld",(int)cstype[i],numtype[i]);
    
    for (j=0;j<numtype[i];j++){
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

  According to original function print from TKr -  TKo, 4.7.2014
*/
void transcrsec::printcrschar (FILE *out, crsectypet ct, long numinst)
{
  switch (ct)
  {
    case nocrosssectiont:{
      break;
    }
    case crsec1dt:{
      cs1d[numinst].print (out);
      break;
    }
    case crsec2dt:{
      cs2d[numinst].print (out);
      break;
    }
    case crsec3dt:{
      cs3d[numinst].print (out);
      break;
    }
    default:{
      print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
    }
  }
}



/**
   function returns thicknesses in nodes
   
   @param nod - vector containing numbers of nodes defining element
   @param t - array containing thicknesses
   
   22.7.2001
*/
void transcrsec::give_thickn (ivector &nod,vector &t)
{
  long i,n;
  n=nod.n;
  for (i=0;i<n;i++){
    switch (Tt->nodes[nod[i]].crst){
    case nocrosssectiont:{
      break;
    }
    case crsec2dt:{
      t[i]=cs2d[Tt->nodes[nod[i]].idcs].t;
      break;
    }
    default:{
      print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
    }
    }
  }
}

/**
   function returns thickness on element
   
   @param eid - element id
   @param t - thickness
   
   22.7.2001
*/
void transcrsec::give_thicke (long eid,double &t)
{
  switch (Tt->elements[eid].crst){
  case nocrosssectiont:{
    break;
  }
  case crsec2dt:{
    t=cs2d[Tt->elements[eid].idcs].t;
    break;
  }
  default:{
    print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function assembles thicknesses in nodes
   
   @param eid - element id
   @param nodes - array of node numbers
   @param th - array of thicknesses in nodes
   
   JK
*/
void transcrsec::give_thickness (long eid,ivector &nodes,vector &th)
{
  long i,nne;
  double t;
  
  if (Tt->elements[eid].crst==0){
    give_thickn (nodes,th);
  }
  else{
    give_thicke (eid,t);
    nne = Tt->give_nne (eid);
    for (i=0;i<nne;i++){
      th[i]=t;
    }
  }
}

/**
   function returns area of cross-section on elements

   10.8.2001
*/
void transcrsec::give_areae (long eid,double &a)
{
  switch (Tt->elements[eid].crst){
  case nocrosssectiont:{
    break;
  }
  case crsec1dt:{
    a=cs1d[Tt->elements[eid].idcs].a;
    break;
  }
  default:{
    print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function returns area of cross-section in nodes
   
   @param nod - vector containing numbers of nodes defining element
   @param t - array containing thicknesses
   
   JK
   10.8.2001, 21.2.2013
*/
void transcrsec::give_arean (ivector &nod,vector &a)
{
  long i,n;
  n=nod.n;
  for (i=0;i<n;i++){
    switch (Tt->nodes[nod[i]].crst){
    case nocrosssectiont:{
      break;
    }
    case crsec1dt:{
      a[i]=cs1d[Tt->nodes[nod[i]].idcs].a;
      break;
    }
    default:{
      print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
    }
    }
  }
}

/**
   function returns area of cross-section
   
   @param eid - element id
   @param nodes - array of node numbers
   @param a - array of areas of cross sections in nodes

   JK
   10.8.2001, 21.2.2013
*/
void transcrsec::give_area (long eid,ivector &nodes,vector &a)
{
  long i,nne;
  double aa;
  
  if (Tt->elements[eid].crst==0){
    give_arean (nodes,a);
  }
  else{
    give_areae (eid,aa);
    //  the number of nodes on element
    nne = Tt->give_nne (eid);
    for (i=0;i<nne;i++){
      a[i]=aa;
    }
  }
}

/**
   function returns densities in nodes
   
   @param nod - vector containing numbers of nodes defining element
   @param rho - array containing densities

   22.7.2001
*/
void transcrsec::give_densityn (ivector &nod,vector &rho)
{
  long i,n;
  n=nod.n;
  for (i=0;i<n;i++){
    switch (Tt->nodes[nod[i]].crst){
    case nocrosssectiont:{
      break;
    }
    case crsec1dt:{
      rho[i]=cs1d[Tt->nodes[nod[i]].idcs].rho;
      break;
    }
    case crsec2dt:{
      rho[i]=cs2d[Tt->nodes[nod[i]].idcs].rho;
      break;
    }
    case crsec3dt:{
      rho[i]=cs3d[Tt->nodes[nod[i]].idcs].rho;
      break;
    }
    default:{
      print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
    }
    }
  }
}

/**
   function returns density on element
   
   @param eid - element id
   @param rho - density

   22.7.2001
*/
void transcrsec::give_densitye (long eid,double &rho)
{
  switch (Tt->elements[eid].crst){
  case nocrosssectiont:{
    break;
  }
  case crsec1dt:{
    rho=cs1d[Tt->elements[eid].idcs].rho;
    break;
  }
  case crsec2dt:{
    rho=cs2d[Tt->elements[eid].idcs].rho;
    break;
  }
  case crsec3dt:{
    rho=cs3d[Tt->elements[eid].idcs].rho;
    break;
  }
  default:{
    print_err("unknown cross section type is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function assembles densities in nodes
   
   @param eid - element id
   @param nodes - array of node numbers
   @param dens - array of densities in nodes
   
   JK
*/
void transcrsec::give_density (long eid,ivector &nodes,vector &dens)
{
  long i,nne;
  double d;
  
  if (Tt->elements[eid].crst==0){
    give_densityn (nodes,dens);
  }
  else{
    give_densitye (eid,d);
    nne = Tt->give_nne (eid);
    for (i=0;i<nne;i++){
      dens[i]=d;
    }
  }
}
