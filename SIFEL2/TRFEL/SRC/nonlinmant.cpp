#include <string.h>
#include "nonlinmant.h"

nonlinmant::nonlinmant (void)
{
  /*
  tnlinsol = (nonlinsolvertype) 0;
  displnorm = (displacementnorm) 0;

  nial = 0;  niilal = 0;  erral = 0.0;  dlal = 0.0;  dlmaxal = 0.0;  psial = 0.0;
  nsnal=0;  selnodal=NULL;  hdbackupal=0;  nsdofal=0;  seldofal=NULL;

  ninr = 0;  niilnr = 0;  errnr = 0.0;  incrnr = 0.0;  minincrnr = 0.0;
  */
  /*
  errl = 0.0;
  nienr = 0;
  hdbr = 0;
  hdbid = 0;
  backupfname[0] = 0;
  
  selnodal=NULL;
  seldofal=NULL;
  */
}

nonlinmant::~nonlinmant (void)
{
  /*
  delete [] selnodal;
  delete [] seldofal;
*/
}

void nonlinmant::read (XFILE *in,long mespr)
{
  //  type of nonlinear solver
  xfscanf (in,"%k%d","type_of_nonlin_solver",(int*)&tnlinsol);
  
  switch (tnlinsol){
  case newtont:{
    //  maximum number of increments
    //  maximum number of iterations in inner loop
    //  required norm of residual
    //  magnitude of increment of loading
    //  minimum magnitude of increment of loading
    //  maximum magnitude of increment of loading
    if (mespr==1)  fprintf (stdout,"\n system of nonlinear equations will be solved by Newton-Raphson method");
    xfscanf (in,"%ld %ld %lf %lf %lf %lf",&ninr,&niilnr,&errnr,&incrnr,&minincrnr,&maxincrnr);
    break;
  }
  default:{
    print_err("unknown solver of nonlinear system of equations",__FILE__,__LINE__,__func__);
  }
  }
  
}

void nonlinmant::print (FILE *out)
{
  fprintf (out,"%d\n",(int)tnlinsol);
  
  switch (tnlinsol){
  case newtont:{
    fprintf (out, "%ld %ld %e %e %e %e\n", ninr, niilnr, errnr, incrnr, minincrnr, maxincrnr);
    break;      
  }
  default:{
    print_err("unknown solver of nonlinear system of equations",__FILE__,__LINE__,__func__);
  }
  }
  
}

/*
void nonlinmant::initiate (nonlinmant &nm)
{
  long i;

  tnlinsol = nm.tnlinsol;

  //  ARC-LENGTH METHOD
  displnorm = nm.displnorm;
  hdbackupal = nm.hdbackupal;
  nial = nm.nial;
  niilal = nm.niilal;
  erral = nm.erral;
  dlal = nm.dlal;
  dlminal = nm.dlminal;
  dlmaxal = nm.dlmaxal;
  psial = nm.psial;
  
  nsnal = nm.nsnal;
  nsdofal = nm.nsdofal;
  probdimal = nm.probdimal;

  nxal = nm.nxal;
  nyal = nm.nyal;
  nzal = nm.nzal;
  
  
  switch (displnorm){
  case alldofs:{  break; }
  case seldofs:{
    selnodal = new long [nsdofal];
    seldofal = new long [nsdofal];
    
    for (i=0;i<nsdofal;i++){
      selnodal[i]=nm.selnodal[i];
      seldofal[i]=nm.seldofal[i];
    }
    
    break;
  }
  case seldofscoord:{
    seldofal = new long [nsdofal];
    selnodcoord = new double [3*nsdofal];
    
    for (i=0;i<nsdofal;i++){
      selnodcoord[3*i]=nm.selnodcoord[3*i];
      selnodcoord[3*i+1]=nm.selnodcoord[3*i+1];
      selnodcoord[3*i+2]=nm.selnodcoord[3*i+2];
      seldofal[i]=nm.seldofal[i];
    }
    
    break;
  }
  case selecnodes:{
    selnodal = new long [nsnal];
    
    for (i=0;i<nsnal;i++){
      selnodal[i]=nm.selnodal[i];
    }
    
    break;
  }
  case nodesdistincr:{
    nsnal=2;
    selnodal = new long [nsnal];
    
    for (i=0;i<nsnal;i++){
      selnodal[i]=nm.selnodal[i];
    }
    
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown norm measurement of displacement increment (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  
  //  NEWTON-RAPHSON METHOD
  ninr = nm.ninr;
  niilnr = nm.niilnr;
  errnr = nm.errnr;
  incrnr = nm.incrnr;
  minincrnr = nm.minincrnr;
  rvlam = nm.rvlam;
  nienr = nm.nienr;
  hdbr  = nm.hdbr;
  if (hdbr)
  {
    hdbid = nm.hdbid;
    for (i=0; i < long(strlen(nm.backupfname)); i++)
      backupfname[i] = nm.backupfname[i];
  }
}
*/
