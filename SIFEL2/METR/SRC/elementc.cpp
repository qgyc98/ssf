#include <string.h>
#include "elementc.h"
#include "globalt.h"
#include "globalc.h"
#include "gtopology.h"
#include "intpoints.h"


elementc::elementc (void)
{
  te = (elemtypec) 0;
  nb=0;
  ippu=NULL;  ippl=NULL;  intordvum=NULL;  intordvlm=NULL;
  tmu=NULL;  tml=NULL;  idmu=NULL;  idml=NULL;
}

elementc::~elementc (void)
{
  long i;

  delete [] tmu;
  delete [] idmu;
  delete [] tml;
  delete [] idml;

  for(i=0; i<Tp->ntm; i++)
  {
    delete [] ippl[i];
    //    delete [] intorvlm[i];
  }
  delete [] ippl;
  //  delete [] intordvlm;

  for(i=0; i<nb; i++)
  {
    //    delete [] intorvum[i];
    delete [] ippu[i];
  }
  delete [] ippu;
  //  delete [] intordvlm;
}

void elementc::read (XFILE *in,long eid)
{
  long ndofe,nne;
  
  xfscanf (in,"%d",(int*)&te);
  
  switch (te){
  case coupbar:{       if (Cbar==NULL)    Cbar = new barelc ();        break;  }
  case coupquad:{      if (Cquad==NULL)   Cquad = new quadrilatc ();   break;  }
  case coupaxiquad:{   if (Caxiq==NULL)   Caxiq = new axiquadc ();     break;  }
  case axisymfc:{      if (Caxifc==NULL)  Caxifc = new axisymc ();     break;  }
  case couphex:{       if (Chex==NULL)    Chex = new hexahedc ();      break;  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
  
  if (Cp->tprob==fully_coupled_material){
    gelemtype get=quadquad;
    
    ndofe=Ct->give_ndofe (eid);
    nne = Ct->give_nne (eid);
    
    Gtu->gelements[eid].read (in,nne,ndofe,get,Gtu->nn);
    
    readmatfc (in);
    
  }else{
    readmat (1,in);
  }
  
}

void elementc::readmat (long m,XFILE *in)
{
  long i;
  
  tmu = new mattypec [m];
  idmu = new long [m];
  tml = new mattypec [m];
  idml = new long [m];

  for (i=0;i<m;i++){
    xfscanf (in,"%d %ld",(int*)&tmu[i],&idmu[i]);
    idmu[i]--;
  }
  for (i=0;i<m;i++){
    xfscanf (in,"%d %ld",(int*)&tml[i],&idml[i]);
    idml[i]--;
  }
}

void elementc::readmatfc (XFILE *in)
{
  xfscanf (in,"%d %ld",(int*)&tm,&idm);
  idm--;
}

