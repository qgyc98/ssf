#include "intpointsc.h"
#include "globalc.h"
#include <string.h>

intpointsc::intpointsc ()
{
  //  material type
  tm = (mattypec) 0;
  //  material id
  idm=0;
  
  //  the number of transported media
  ntm=0;
  //  the number of components in gradient/flux arrays
  ncompgrad=0;
  
  //  actual values in transport processes
  av=NULL;
  //  values from the previous step in transport processes
  pv=NULL;
  //  gradients of transport values
  grad=NULL;
  //  fluxes of transport values
  fluxes=NULL;
  
  //  number of components of displacement
  ncompdispl=0;
  //  number of components of stress/strain array
  ncompstr=0;
  //  number of components of eqother array
  ncompeqother=0;
  //  number of components of other array
  ncompother=0;
  
  //  array of displacements
  displ=NULL;
  //  array of strains
  strain=NULL;
  //  array of stresses
  stress=NULL;
  //  other components
  other = NULL;
  //  equilibriated components of other array
  eqother = NULL; 
  
  //  stress-strain state
  ssst = (strastrestate) 0;
  
}

intpointsc::~intpointsc (void)
{
  
  delete [] displ;
  delete [] strain;
  delete [] stress;
  delete [] other;
  delete [] eqother;

  delete [] av;
  delete [] pv;
  for (long i=0; i<Tp->ntm; i++){
    delete [] grad[i];
    delete [] fluxes[i];
  }
  delete [] grad;
  delete [] fluxes;

}

/**
   function allocates arrays on integration points
   
   @param ipp - integration point id
   
   JK, 10.4.2019
*/
void intpointsc::alloc (long ipp)
{
  long i,j;
  
  displ = new double [ncompdispl];
  memset(displ, 0, sizeof(*displ)*ncompdispl);

  strain = new double [ncompstr];
  memset(strain, 0, sizeof(*strain)*ncompstr);

  stress = new double [ncompstr];
  memset(stress, 0, sizeof(*stress)*ncompstr);

  //  number of components of array other
  ncompother = Cm->givencompother (ipp,0);
  if (ncompother > 0){
    other = new double[ncompother];
    memset(other, 0, sizeof(*other)*ncompother);
  }
  
  //  number of components of array eqother
  ncompeqother = Cm->givencompeqother (ipp,0);
  if (ncompeqother > 0){
    eqother = new double[ncompeqother];
    memset(eqother, 0, sizeof(*eqother)*ncompeqother);
  }
  
  
  av = new double [ntm];
  memset(strain, 0, sizeof(*av)*ntm);
  
  pv = new double [ntm];
  memset(strain, 0, sizeof(*av)*ntm);
  
  grad = new double* [ntm];
  fluxes = new double* [ntm];
  for (i=0;i<ntm;i++){
    grad[i] = new double [ncompgrad];
    fluxes[i] = new double [ncompgrad];
  }
  
  for (i=0;i<ntm;i++){
    av[i]=0.0;
    pv[i]=0.0;
    for (j=0;j<ncompgrad;j++){
      grad[i][j]=0.0;
      fluxes[i][j]=0.0;
    }
  }

}

