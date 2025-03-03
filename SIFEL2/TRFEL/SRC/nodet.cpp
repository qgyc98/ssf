#include "nodet.h"
#include "globalt.h"

nodet::nodet (void)
{
  //  number of degrees of freedom
  ndofn=0;
  //  number of components in array of gradients
  ncompgrad=0;
  //  number of components in array other
  ncompother = 0;
  //  number of components in array eqother
  ncompeqother=0;
  
  //  type of cross section
  crst = (crsectypet) 0;
  //  cross section id
  idcs=0;
  
  //  array of contributions to the gradients
  ncontr_grad=NULL;
  vol_grad=NULL;
  //  array of contributions to the fluxes
  ncontr_flux=NULL;
  vol_flux=NULL;
  //  number of contributions to the array other
  ncontr_other=NULL;
  vol_other=NULL;
  //  number of contributions to the array eqother
  ncontr_eqother=NULL;
  vol_eqother=NULL;
  
  //  array of gradients
  gradient=NULL;
  //  array of fluxes
  flux=NULL;
  //  array of auxiliary components
  other=NULL;
  //  array of equilibriated auxiliary components
  eqother=NULL;
  
  
  //  actual nodal values
  nodval=NULL;
  //  nodal values from previous step
  nodvalp=NULL;
  //  initial nodal values
  nodvali=NULL;
  //  time derivatives of nodal values
  nodvalt=NULL;
  // vector of weights for master node contributions at hanging node
  mnw = NULL;
}

nodet::~nodet (void)
{
  long i;

  for (i=0; i<Tp->ntm; i++){
    if (gradient != NULL)
      delete [] gradient[i];
    if (flux != NULL)
      delete [] flux[i];
  }
  delete [] gradient;
  delete [] flux;
  
  if (ncontr_grad!=NULL)
    delete [] ncontr_grad;
  if (ncontr_flux!=NULL)
    delete [] ncontr_flux;
  if (ncontr_other!=NULL)
    delete [] ncontr_other;
  if (ncontr_eqother!=NULL)
    delete [] ncontr_eqother;
  if (other!=NULL)
    delete [] other;
  if (eqother!=NULL)
    delete [] eqother;
  
  if (nodval!=NULL)
    delete [] nodval;
  if (nodvalp!=NULL)
    delete [] nodvalp; 
  if (nodvali!=NULL)
    delete [] nodvali; 
  if (nodvalt!=NULL)
    delete [] nodvalt;
 
  delete mnw;
}

/**
   function reads data stored on nodes
   
   @param in - input file
   @param ndof - number of DOFs on node
   
   JK, 10.7.2008 - revision
*/
void nodet::read (XFILE *in,long ndof)
{
  //  number of DOFs on node
  ndofn=ndof;
  
  //  type of cross section
  xfscanf (in,"%d",(int*)&crst);
  if (crst!=0){
    xfscanf (in,"%ld",&idcs);
    idcs--;
  }
}

/**
   function prints data
   
   @param out - output file
*/
void nodet::print (FILE *out)
{
  //  type of cross section
  fprintf (out," %d ",(int)crst);
  if (crst!=0){
    fprintf (out," %ld\n",idcs+1);
  }
}




/**
   function allocates gradient arrays on nodes
   
   @param ncompo - number components

   18.5.2002
   modified by TKr 21/01/2010
*/
void nodet::alloc_grad (long ncomp)
{
  long i;
  
  ncompgrad=ncomp;
  gradient = new double* [Tp->ntm];
  if (Tp->gradaver==1)
    ncontr_grad = new long [Tp->ntm];
  if (Tp->gradaver==2)
    vol_grad = new double [Tp->ntm];
  
  for (i=0;i<Tp->ntm;i++){
    gradient[i] = new double [ncomp];
  }
  nullgrad();
}

/**
   function allocates fluxes arrays on nodes
   
   @param ncompo - number components
   
   18.5.2002
   modified by TKr 21/01/2010
*/
void nodet::alloc_flux (long ncomp)
{
  long i;
  
  ncompgrad=ncomp;
  flux = new double* [Tp->ntm];
  if (Tp->fluxaver==1)
    ncontr_flux = new long [Tp->ntm];
  if (Tp->fluxaver==2)
    vol_flux = new double [Tp->ntm];
  
  for (i=0;i<Tp->ntm;i++){
    flux[i] = new double [ncomp];
  }
  nullflux();
}


/**
   function allocates other arrays on nodes
   
   @param ncompo - number components
   
   18.5.2002
*/
void nodet::alloc_other (long ncompo)
{
  ncompother=ncompo;
  if (Tp->otheraver==1)
    ncontr_other = new long [ncompo];
  if (Tp->otheraver==2)
    vol_other = new double [ncompo];

  other  = new double [ncompo];
  nullother ();
}


/**
   function allocates eq_other arrays on nodes
   
   @param ncompo - number components
   
   18.5.2002
*/
void nodet::alloc_eqother (long ncompeqo)
{
  ncompeqother=ncompeqo;
  if (Tp->eqotheraver==1)
    ncontr_eqother = new long [ncompeqo];
  if (Tp->eqotheraver==2)
    vol_eqother = new double [ncompeqo];

  eqother  = new double [ncompeqo];
  nulleqother ();
}

/**
   function allocates array nodval for actual nodal values on nodes
   
   JK, 10.7.2008
*/
void nodet::alloc_nodval ()
{
  nodval = new double [ndofn];
  nullnodval ();
}

/**
   function allocates array nodvali for initial nodal values on nodes
   
   JK, 10.7.2008
*/
void nodet::alloc_nodvali ()
{
  nodvali = new double [ndofn];
  nullnodvali ();
}

/**
   function allocates array nodvalp for previous nodal values on nodes
   
   JK, 10.7.2008
*/
void nodet::alloc_nodvalp ()
{
  nodvalp = new double [ndofn];
  nullnodvalp ();
}

/**
   function allocates array nodvalt for time derivatives of actual nodal values on nodes
   
   JK, 10.7.2008
*/
void nodet::alloc_nodvalt ()
{
  nodvalt = new double [ndofn];
  nullnodvalt ();
}


/**
   function sets up all components in array gradient to zero
   
   JK, 10.7.2008
   modified by TKr 21/01/2010
*/
void nodet::nullgrad ()
{
  long i,j;

  for (i=0;i<Tp->ntm;i++){
    if (Tp->gradaver==1)
      ncontr_grad[i] = 0;
    if (Tp->gradaver==2)
      vol_grad[i] = 0.0;
    for (j=0;j<ncompgrad;j++)
      gradient[i][j]=0.0;
  }
}

/**
   function sets up all components in array flux to zero
   
   JK, 10.7.2008
   modified by TKr 21/01/2010
*/
void nodet::nullflux ()
{
  long i,j;
  
  for (i=0;i<Tp->ntm;i++){
    if (Tp->fluxaver==1)
      ncontr_flux[i] = 0;
    if (Tp->fluxaver==2)
      vol_flux[i] = 0.0;
    for (j=0;j<ncompgrad;j++)
      flux[i][j]=0.0;
  }
}

/**
   function sets up all components in array other to zero
   
*/
void nodet::nullother ()
{
  long i;
  
  for (i=0;i<ncompother;i++){
    other[i]=0.0;
    if (Tp->otheraver==1)
      ncontr_other[i]=0;
    if (Tp->otheraver==2)
      vol_other[i]=0;
  }
  
}

/**
   function sets up all components in array eqother to zero
   
   @param ntm - number of transported media
*/
void nodet::nulleqother ()
{
  long i;
  
  for (i=0;i<ncompeqother;i++){
    eqother[i]=0.0;
    if (Tp->eqotheraver==1)
      ncontr_eqother[i]=0;
    if (Tp->eqotheraver==2)
      vol_eqother[i]=0;
  }

}

/**
   function cleans array nodval
   
   JK, 10.7.2008
*/
void nodet::nullnodval ()
{
  long i;
  for (i=0;i<ndofn;i++){
    nodval[i]=0.0;
  }
}

/**
   function cleans array nodvalp
   
   JK, 10.7.2008
*/
void nodet::nullnodvalp ()
{
  long i;
  for (i=0;i<ndofn;i++){
    nodvalp[i]=0.0;
  }
}

/**
   function cleans array nodvali
   
   JK, 10.7.2008
*/
void nodet::nullnodvali ()
{
  long i;
  for (i=0;i<ndofn;i++){
    nodvali[i]=0.0;
  }
}

/**
   function cleans array nodvalt
   
   JK, 10.7.2008
*/
void nodet::nullnodvalt ()
{
  long i;
  for (i=0;i<ndofn;i++){
    nodvalt[i]=0.0;
  }
}



/**
   function stores gradient components
   
   @param lcid - loadcase id = number of unknown
   @param gradv - array containing part of gradient components
   
   19.5.2002
   modified by TKr 21/01/2010
*/
void nodet::storegrad (long lcid,vector &gradv)
{
  long i;
  
  ncontr_grad[lcid]++;
  
  for (i=0;i<ncompgrad;i++){
    gradient[lcid][i]+=gradv[i];
  }
}


/**
   function stores gradient components
   
   @param lcid - loadcase id = number of unknown
   @param vol - volume used for contribution
   @param gradv - array containing part of gradient components
   
   19.5.2002
   modified by TKr 21/01/2010
*/
void nodet::storegrad (long lcid,double vol,vector &gradv)
{
  long i;
  
  vol_grad[lcid]+=vol;
  
  for (i=0;i<ncompgrad;i++){
    gradient[lcid][i]+=vol*gradv[i];
  }
}


/**
   function stores stress components
   
   @param lcid - loadcase id = number of unknown
   @param fluxv - array containing part of flux components
   
   19.5.2002
   modified by TKr 21/01/2010
*/
void nodet::storeflux (long lcid,vector &fluxv)
{
  long i;
  
  ncontr_flux[lcid]++;
  
  for (i=0;i<ncompgrad;i++){
    flux[lcid][i]+=fluxv[i];
  }
}

/**
   function stores gradient components
   
   @param lcid - loadcase id = number of unknown
   @param vol - volume used for contribution
   @param fluxv - array containing part of flux components
   
   19.5.2002
   modified by TKr 21/01/2010
*/
void nodet::storeflux (long lcid,double vol,vector &fluxv)
{
  long i;
  
  vol_flux[lcid]+=vol;
  
  for (i=0;i<ncompgrad;i++){
    flux[lcid][i]+=vol*fluxv[i];
  }
}


/**
   function stores other components

   @param ncomp - number of component
   @param otherv - other value
   
   19.5.2002
*/
void nodet::storeother (long ncomp,double otherv)
{
  ncontr_other[ncomp]++;
  other[ncomp]+=otherv;
}

/**
   function stores other components

   @param ncomp - number of component
   @param vol - volume
   @param otherv - other value
   
   23. 3. 2017, JK
*/
void nodet::storeother (long ncomp,double vol,double otherv)
{
  vol_other[ncomp]+=vol;
  other[ncomp]+=otherv*vol;
}

/**
   function stores eq_other components

   @param ncompeq - number of component
   @param eqotherv - array containing part of eq_other components
   
   19.5.2002
*/
void nodet::storeeqother (long ncompeq,double eqotherv)
{
  ncontr_eqother[ncompeq]++;
  eqother[ncompeq]+=eqotherv;
}

/**
   function stores eq_other components

   @param ncompeq - number of component
   @param vol - volume
   @param eqotherv - array containing part of eq_other components
   
   23. 3. 2017, JK
*/
void nodet::storeeqother (long ncompeq,double vol,double eqotherv)
{
  vol_eqother[ncompeq]+=vol;
  eqother[ncompeq]+=eqotherv*vol;
}



/**
   The function returns nodal gradient components of the given lcid
   
   @param lcid - loadcase id = number of unknown [in]
   @param gv - %vector of gradient components [out]
   
   @return The function does not return anything but it stores gradient components in argument gv.

   Created by Tomas Koudelka, 02.2018
*/
void nodet::givegrad (long lcid, vector &gv)
{
  if (lcid < ndofn)
    copyv(gradient[lcid], gv);
  else
    print_err("gradient load case id %ld is out of allowed range <1;%ld>", __FILE__, __LINE__, __func__, lcid+1, ndofn); 
}



/**
   The function returns compid-th nodal gradient components of the given lcid
   
   @param lcid - loadcase id = number of unknown [in]
   @param compid - gradient component id [in]
  
   @return The function returns required gradient component.

   Created by Tomas Koudelka, 02.2018
*/
double nodet::givegrad (long lcid, long compid)
{
  if ((lcid < ndofn) && (compid < ncompgrad))
    return gradient[lcid][compid];

  if (compid >= ndofn)
    print_err("gradient load case id %ld is out of range <1;%ld>", __FILE__, __LINE__, __func__, lcid+1, ndofn);

  if (compid >= ncompgrad)
    print_err("gradient component id %ld is out of range <1;%ld>", __FILE__, __LINE__, __func__, compid+1, ncompgrad);

  return 0.0;
}


/**
   The function returns nodal flux components of the given lcid
   
   @param lcid - loadcase id = number of unknown [in]
   @param fv - %vector of flux components [out]
   
   @return The function does not return anything but it stores flux components in argument fv.

   Created by Tomas Koudelka, 02.2018
*/
void nodet::giveflux (long lcid, vector &fv)
{
  if (lcid < ndofn)
    copyv(flux[lcid], fv);
  else
    print_err("flux load case id %ld is out of allowed range <1;%ld>", __FILE__, __LINE__, __func__, lcid+1, ndofn); 
}



/**
   The function returns compid-th nodal flux components of the given lcid
   
   @param lcid - loadcase id = number of unknown [in]
   @param compid - flux component id [in]
  
   @return The function returns required flux component.

   Created by Tomas Koudelka, 02.2018
*/
double nodet::giveflux (long lcid, long compid)
{
  if ((lcid < ndofn) && (compid < ncompgrad))
    return flux[lcid][compid];

  if (compid >= ndofn)
    print_err("flux load case id %ld is out of range <1;%ld>", __FILE__, __LINE__, __func__, lcid+1, ndofn);

  if (compid >= ncompgrad)
    print_err("flux component id %ld is out of range <1;%ld>", __FILE__, __LINE__, __func__, compid+1, ncompgrad);

  return 0.0;
}



/**
   The function returns compid-th component of other state variable array at the given node.

   @param compid - id of other array component [input]

   @return  The function returns compid-th component of other array.
   
   Created by Tomas Koudelka, 02.2018
*/
double nodet::giveother (long compid)
{
  if (compid < ncompother)
    return other[compid];
  else
    print_err("other component id %ld is out of range <1;%ld>", __FILE__, __LINE__, __func__, compid+1, ncompother);

  return 0.0;
}



/**
   The function copies nodal other state variable array at the given node to the vector othv.

   @param othv - %vector of state variables [output]

   @return  The function does not return anything but stores the nodal state variable vector eqother in othv argument.
   
   Created by Tomas Koudelka, 02.2018
*/
void nodet::giveother (vector &othv)
{
  copyv(other, othv);
}



/**
   The function returns compid-th component of eqother state variable array at the given node.

   @param compid - id of eqother array component [input]

   @return  The function returns compid-th component of eqother array.
   
   Created by Tomas Koudelka, 02.2018
*/
double nodet::giveeqother (long compid)
{
  if (compid < ncompeqother)
    return eqother[compid];
  else
    print_err("eqother component id %ld is out of range <1;%ld>", __FILE__, __LINE__, __func__, compid+1, ncompeqother);

  return 0.0;
}



/**
   The function copies nodal eqother state variable array at the given node to the vector eqothv.

   @param eqothv - %vector of state variables [output]

   @return  The function does not return anything but stores the nodal state variable vector eqother in eqothv argument.
   
   Created by Tomas Koudelka, 02.2018
*/
void nodet::giveeqother (vector &eqothv)
{
  copyv(eqother, eqothv);
}



/**
   function computes average gradient at node
   modified by TKr 21/01/2010
*/
void nodet::grad_averageval ()
{
  long i,j;
  
  if (Tp->gradaver==1){
    for (i=0;i<Tp->ntm;i++){
      if (ncontr_grad[i]>0){
	for (j=0;j<ncompgrad;j++){
          if (ncontr_grad[i] > 0)
            gradient[i][j]/=ncontr_grad[i];
	}
      }
    }
  }

  if (Tp->gradaver==2){
    for (i=0;i<Tp->ntm;i++){
      if (vol_grad[i]>0){
	for (j=0;j<ncompgrad;j++){
          if (vol_grad[i] > 0.0)
            gradient[i][j]/=vol_grad[i];
	}
      }
    }
  }

}

/**
   function computes average flux at node
   modified by TKr 21/01/2010
*/
void nodet::flux_averageval ()
{
  long i,j;
  
  if (Tp->fluxaver==1){
    for (i=0;i<Tp->ntm;i++){
      if (ncontr_flux[i]>0){
	for (j=0;j<ncompgrad;j++){
          if (ncontr_flux[i] > 0)
            flux[i][j]/=ncontr_flux[i];
	}
      }
    }
  }

  if (Tp->fluxaver==2){
    for (i=0;i<Tp->ntm;i++){
      if (vol_flux[i]>0){
	for (j=0;j<ncompgrad;j++){
          if (vol_flux[i] > 0.0)
            flux[i][j]/=vol_flux[i];
	}
      }
    }
  }

}


/**
   function computes average values from array other
   modified by TKr 21/01/2010
*/
void nodet::other_averageval ()
{
  long i;
  
  if (Tp->otheraver==1){
    if (ncontr_other){
      for (i=0;i<ncompother;i++){
        if (ncontr_other[i] > 0)
          other[i]/=ncontr_other[i];
      }
    }
  }
  
  if (Tp->otheraver==2){
    if (vol_other){
      for (i=0;i<ncompother;i++){
        if (vol_other[i] > 0.0)
          other[i]/=vol_other[i];
      }
    }
  }

}


/**
   function computes average values from array eq_other
   modified by TKr 21/01/2010
*/
void nodet::eqother_averageval ()
{
  long i;

  if (Tp->eqotheraver==1){
    if (ncontr_eqother){
      for (i=0;i<ncompeqother;i++){
        if (ncontr_eqother[i] > 0)
          eqother[i]/=ncontr_eqother[i];
      }
    }
  }
  
  if (Tp->eqotheraver==2){
    if (vol_eqother){
      for (i=0;i<ncompeqother;i++){
        if (vol_eqother[i] > 0.0)
          eqother[i]/=vol_eqother[i];
      }
    }
  }

}


/**
   function is used in function transmat.cpp::compute_jump

   @param in - array of input data
   
   JK, 7.1.2008
*/
void nodet::give_values (double *in, double *inp, double *ineq)
{
  long i;
  
  //  actual nodal values
  for (i=0;i<ndofn;i++){
    in[i]=nodval[i];
  }
  
  //  nodal values from the previous time step
  for (i=0;i<ndofn;i++){
    inp[i]=nodvalp[i];
  }
  
  //for (i=0;i<ncompeqother-5;i++){
  //ineq[i]=eqother[i];
  //}
  for (i=0;i<5;i++){
    ineq[i]=eqother[i];
  }
}

/**
   function is used in function transmat.cpp::compute_jump

   @param out - array of saved data
   
   JK, 7.1.2008
*/
void nodet::save_values (double *out)
{
  long i;
  
  for (i=0;i<5;i++){
    eqother[i]=out[i];
  }
}

/**
   function is used in function transmat.cpp::compute_jump

   function saves nodal values to array nodval
   
   @param nv - array of nodal values
   
   JK, 8.1.2008
*/
void nodet::save_nodval (double *nv)
{
  long i;
  
  for (i=0;i<ndofn;i++){
    nodval[i]=nv[i];
  }
}

/**
   function moves array of actual values to array of previous values
   
   JK, 29.5.2007
*/
void nodet::actual_previous_change ()
{
  long i;
    
  for (i=0;i<ndofn;i++){
    nodvalp[i]=nodval[i];
  }
}

/**
   function returns %vector of flux in node
   
   @param lcid - load case id
   @param nv - array containing nodal flux
   
   JK, 15. 2. 2011
*/
void nodet::give_flux (long lcid,double *nv)
{
  long i;
  
  for (i=0;i<ndofn;i++){
    nv[i]=flux[lcid][i];
  }
}

