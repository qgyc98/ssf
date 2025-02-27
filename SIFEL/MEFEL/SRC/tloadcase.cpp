#include "tloadcase.h"
#include "matrix.h"
#include "vector.h"
#include "gfunct.h"
#include "loadcase.h"
#include "global.h"
#include "gtopology.h"
#include "globmat.h"
#include "dloadn.h"
#include "dloadpd.h"
#include "element.h"
#include "node.h"
#include <math.h>



tloadcase::tloadcase (void)
{
  nslc=0;
  nln = 0;  nle = 0;  npd = 0;
  lon = NULL;  pd = NULL;
  slc = NULL;  gf=NULL;
  direction = NULL;
  seism=NULL;
}

tloadcase::~tloadcase (void)
{
  delete [] lon;  delete [] pd;
  delete [] slc;  delete [] gf;
  delete [] direction;
  delete [] seism;
}

/**
   function reads load case characteristics
   
   @param in - input stream
   
   JK, TKo, 3.6.2005
*/
void tloadcase::read (FILE *in)
{
  long i;
  
  //  type of time dependent load
  fscanf (in,"%d",(int*)&tdl);
  
  switch (tdl)
  {
    case forcedload: 
    {
      //  number of subload cases
      fscanf (in,"%ld",&nslc);
    
      slc = new loadcase [nslc];
      gf = new gfunct [nslc];
    
      for (i=0;i<nslc;i++)
      {
        slc[i].read (in);
        gf[i].read (in);
      }
    
      break;
    }
    case seismicload:  
    {
      //  number of seismic acceleration components
      fscanf (in,"%ld",&nslc);
      direction = new dirdynload [nslc];
      gf = new gfunct [nslc];

      for (i=0;i<nslc;i++)
        fscanf (in,"%d",(int*)&direction[i]);

      for (i=0;i<nslc;i++)
        gf[i].read (in);

      seism = new double [Ndofm*nslc];
      nullv (seism,Ndofm*nslc);
      seisminit (seism);
      break;
    }
    case forceload_dyn: 
    {
      //  loaded nodes
      fscanf (in,"%ld",&nln);
      lon = new dloadn [nln];
      for (i=0;i<nln;i++)
        lon[i].read (in);

      //  loaded elements
//    fscanf (in,"%ld",&nle);
//    loe = new dloadel [nle];
//    for (i=0;i<nle;i++)
//      loe[i].read (in);
    
      //  prescribed displacements
      fscanf (in,"%ld",&npd);
      pd = new dloadpd [npd];
      for (i=0;i<npd;i++){
        pd[i].read (in);
      }
      break;
    }
  
    default:
      fprintf (stderr,"\n\n unknown type of load case is required in function read (file %s, line %d).\n",__FILE__,__LINE__);
  }
}

