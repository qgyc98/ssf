#include "gsurf.h"
#include "gnode.h"
#include "gtopology.h"
#include <math.h>

gsurf::gsurf (void)
{
  //  number of nodes on surface
  nn=0;
  //  number of DOFs defined on nodes
  ndofn=NULL;
  
  //  number of reference element
  re=-1;
  //  list of adjacent elements
  adjel=NULL;
  
  //  list of node numbers
  nlistfn = NULL;
  nlistln = NULL;
  //  array of code numbers (for Lagrange multipliers)
  cnmult = NULL;
  
  //  normal vector
  nv = NULL;

  //  threshold
  threshold=1.0e-6;
}


gsurf::~gsurf (void)
{
  long i;
  
  delete [] ndofn;
  delete [] adjel;
  delete [] nlistfn;
  delete [] nlistln;
  delete [] nv;
  
  for (i=0;i<nn;i++){
    delete [] cnmult[i];
  }
  delete [] cnmult;
}

/**
   function computes normal %vector
   
   @param top - pointer to topology
   
   JK, 9.7.2007
*/
void gsurf::normal_vector (gnode */*gnodes*/)
{
  double norm;
  
  nv = new double [3];
  
  print_err("nv array is used uninitialized!", __FILE__, __LINE__, __func__);
  norm = nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2];
  
  if (fabs(norm)<threshold){
    print_err("zero length of the normal vector", __FILE__, __LINE__, __func__);
  }
  
  norm = sqrt(norm);

  nv[0]/=norm;
  nv[1]/=norm;
  nv[2]/=norm;
}

/**
   function prints
   
   @param out - output stream
   
   JK, 11.7.2007
*/
void gsurf::print (FILE *out)
{
  long i;
  
  fprintf (out,"\n");
  fprintf (out,"\n number of nodes    %ld",nn);
  fprintf (out,"\n node multiplicity  %ld",nm);
  
  if (nlistfn!=NULL){
    fprintf (out,"\n surface nodes\n");
    for (i=0;i<nn;i++){
      fprintf (out,"   %ld",nlistfn[i]);
    }
    fprintf (out,"\n");
    for (i=0;i<nn;i++){
      fprintf (out,"   %ld",nlistln[i]);
    }
  }
  if (nv!=NULL){
    fprintf (out,"\n normal vector");
    fprintf (out,"\n %f  %f  %f",nv[0],nv[1],nv[2]);
  }

}

/**
   function returns normal %vector
   
   JK, 28.8.2007
*/
void gsurf::give_norvect (double *v)
{
  v[0]=nv[0];
  v[1]=nv[1];
  v[2]=nv[2];
}

/**
   function checks normal %vector
   
   JK, 21.10.2007
*/
void gsurf::check_normal (vector &/*x*/,vector &/*y*/,ivector &/*nod*/)
{
  /*
  long i,nne;
  double ss,xc,yc,vx,vy,xf,yf;
  
  //  number of nodes on element
  nne = x.n;
  
  //  center of gravity of element
  xc=0.0;  yc=0.0;
  for (i=0;i<nne;i++){
    xc+=x[i];
    yc+=y[i];
  }
  xc/=nne;
  yc/=nne;
  
  //  coordinates of the first node
  for (i=0;i<nne;i++){
    if (fn==nod[i]){
      xf=x[i];
      yf=y[i];
      break;
    }
  }
  
  //  auxiliary vector (center of gravity -> first node)
  vx=xf-xc;
  vy=yf-yc;
  
  //  dot product of auxiliary vector and normal vector
  ss = vx*nv[0] + vy*nv[1];
  
  if (ss<0.0){
    //  normal vector is inner, not outer vector
    //  orientation has to be changed
    nv[0]*=-1.0;
    nv[1]*=-1.0;
    nv[2]*=-1.0;
  }
  
  */
}

/**
   function allocates array for code numbers of Lagrange multipliers
   
   @param nccnfn - number of components in array cn of the first nodes
   @param nccnln - number of components in array cn of the last nodes
   
   JK, 6.1.2008
*/
void gsurf::alloc_cn (long nccnfn,long nccnln)
{
  long i,ndofnf,ndofnl;
  
  //  number of DOFs defined on first nodes
  ndofnf=nccnfn;
  //  number of DOFs defined on last nodes
  ndofnl=nccnln;
  
  if (ndofnf!=ndofnl){
    
  }
  else{}
    
  if (cnmult!=NULL){
    for (i=0;i<nn;i++){
      delete [] cnmult[i];
    }
    delete [] cnmult;
  }
  
  cnmult = new long* [nn];
  for (i=0;i<nn;i++){
    cnmult[i] = new long [ndofnf];
  }
}

/**
   function assembles node numbers of first nodes
   
   @param fnn - array containing node numbers of first nodes
   
   JK, 6.8.2008
*/
void gsurf::give_first_node_numbers (long *fnn)
{
  long i;
  
  for (i=0;i<nn;i++){
    fnn[i]=nlistfn[i];
  }
}

/**
   function assembles node numbers of last nodes
   
   @param lnn - array containing node numbers of last nodes
   
   JK, 6.8.2008
*/
void gsurf::give_last_node_numbers (long *lnn)
{
  long i;
  
  for (i=0;i<nn;i++){
    lnn[i]=nlistln[i+nn];
  }
}


/**
   function assembles code numbers of Lagrange multipliers
   defined between selected nodes
   
   @param fln - first or last node indicator
          fln=1 - first nodes are assumed
	  fln=2 - last nodes are assumed
   @param mcn - code numbers of selected multipliers
   
   JK, 6.8.2008
*/
void gsurf::give_mult_code_numbers (long /*fln*/,long */*mcn*/)
{
/*
  long i;
  
  if (fln==1){
    for (i=0;i<ndofnf;i++){
      mcn[i]=cnfn[i];
    }
  }
  if (fln==2){
    for (i=0;i<ndofnl;i++){
      mcn[i]=cnln[i];
    }
  }
*/
}

/*
ostream& operator<<(ostream &os, gsurf &gs)
{
	os << "Stycna plocha\n-------------\n\n";
	os << "Pocet prvku [nae] : " << gs.nae << "\n\n";
	os << "Prvky [adjel] : ";
	for(long i = 0; i < gs.nae; i++)
		os << gs.adjel[i]+1 << " ";
	os << "\n\n";
	os << "Pocet bodu [nn] : " << gs.nn << "\n\n";
	os << "Dvojice stycnych bodu [nlist1 nlist 2 |] : ";
	for(long i = 0; i < gs.nn; i++)
		os << gs.nlist1[i]+1 << " " << gs.nlist2[i]+1 << " | ";
	cout << "\n\n\n\n";
	
	return os;
}

*/
