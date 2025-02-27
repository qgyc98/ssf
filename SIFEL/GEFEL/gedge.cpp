#include "gedge.h"
#include "gnode.h"
#include "gtopology.h"
#include <math.h>

gedge::gedge (void)
{
  //  number of nodes on edge
  nn=0;
  //  number of first node
  fn=-1;
  //  number of last node
  ln=-1;
  //  node multiplicity
  nm=0;
  //  number of DOFs defined on first nodes
  ndofnf=0;
  //  number of DOFs defined on last nodes
  ndofnl=0;


  //  previous edge
  prev=-1;
  //  next edge
  next=-1;
  
  //  number of reference element
  re=-1;
  //  list of adjacent elements
  adjel=NULL;

  //  list of node numbers
  nlist = NULL;
  //  array of code numbers of the first nodes
  cnfn = NULL;
  //  array of code numbers of the last nodes
  cnln = NULL;
  
  //  edge length
  l=0.0;
  //  direction vector
  dv = NULL;
  //  normal vector
  nv = NULL;

  //  threshold
  threshold=1.0e-6;
}


gedge::~gedge (void)
{
  delete [] adjel;
  delete [] nlist;
  delete [] cnfn;
  delete [] cnln;
  delete [] dv;
  delete [] nv;
}

/**
   function computes direction %vector
   
   @param top - pointer to topology
   
   JK, 9.7.2007
*/
void gedge::direction_vector (gnode *gnodes)
{
  if (fn<0){
    print_err("first node on edge is not defined", __FILE__, __LINE__, __func__);
  }
  if (ln<0){
    print_err("last node on edge is not defined", __FILE__, __LINE__, __func__);
  }
  
  dv = new double [3];
  
  dv[0] = gnodes[ln].x - gnodes[fn].x;
  dv[1] = gnodes[ln].y - gnodes[fn].y;
  dv[2] = gnodes[ln].z - gnodes[fn].z;
  
  //  edge length
  l = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
  
  if (fabs(l)<threshold){
    print_err("zero length of the direction vector", __FILE__, __LINE__, __func__);
  }
  
  l = sqrt(l);

  dv[0]/=l;
  dv[1]/=l;
  dv[2]/=l;
}

/**
   function computes normal %vector
   
   @param top - pointer to topology
   
   JK, 9.7.2007
*/
void gedge::normal_vector (gnode *gnodes)
{
  double norm;
  
  if (dv==NULL){
    direction_vector (gnodes);
  }
  
  nv = new double [3];
  
  if (fabs(dv[2])<threshold){
    //  problem is in plane x-y
    nv[2]=0.0;
    
    if (fabs(dv[0])<threshold){
      nv[0]=1.0;
      nv[1]=0.0;
    }
    else{
      nv[1]=1.0;
      nv[0]=0.0-dv[1]*nv[1]/dv[0];
    }
  }
  else{
    if (fabs(dv[1])<threshold){
      nv[1]=0.0;
      
      if (fabs(dv[0])<threshold){
	nv[0]=1.0;
	nv[2]=0.0;
      }
      else{
	nv[2]=1.0;
	nv[0]=0.0-dv[2]*nv[2]/dv[0];
      }
    }
    else{
      nv[2]=0.0;
      nv[0]=1.0;
      nv[1]=0.0-dv[0]*nv[0]/dv[1];
    }
  }
  
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
void gedge::print (FILE *out)
{
  long i;
  
  fprintf (out,"\n");
  fprintf (out,"\n number of nodes    %ld",nn);
  fprintf (out,"\n node multiplicity  %ld",nm);
  fprintf (out,"\n number of previous edge  %ld",prev);
  fprintf (out,"\n number of next edge      %ld",next);

  if (nlist!=NULL){
    fprintf (out,"\n egde nodes\n");
    for (i=0;i<nn;i++){
      fprintf (out,"   %ld",nlist[i]);
    }
  }
  if (dv!=NULL){
    fprintf (out,"\n direction vector");
    fprintf (out,"\n %f  %f  %f",dv[0],dv[1],dv[2]);
  }
  if (nv!=NULL){
    fprintf (out,"\n normal vector");
    fprintf (out,"\n %f  %f  %f",nv[0],nv[1],nv[2]);
  }
  
}

/**
   function returns directional %vector
   
   JK, 28.8.2007
*/
void gedge::give_dirvect (double *v)
{
  v[0]=dv[0];
  v[1]=dv[1];
}

/**
   function returns normal %vector
   
   JK, 28.8.2007
*/
void gedge::give_norvect (double *v)
{
  v[0]=nv[0];
  v[1]=nv[1];
}

/**
   function checks normal %vector
   
   JK, 21.10.2007
*/
void gedge::check_normal (vector &x,vector &y,ivector &nod)
{
  long i,nne;
  double ss,xc,yc,vx,vy,xf=0.0,yf=0.0;
  
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
  
  
}

/**
   function allocates array for code numbers of Lagrange multipliers
   
   @param nccnfn - number of components in array cn of the first nodes
   @param nccnln - number of components in array cn of the last nodes
   
   JK, 6.1.2008
*/
void gedge::alloc_cn (long nccnfn,long nccnln)
{
  //  number of DOFs defined on first nodes
  ndofnf=nccnfn;
  //  number of DOFs defined on last nodes
  ndofnl=nccnln;
  
  if (cnfn!=NULL)
    delete [] cnfn;
  cnfn = new long [nccnfn];
  if (cnln!=NULL)
    delete [] cnln;
  cnln = new long [nccnln];
}

/**
   function assembles node numbers of first nodes
   
   @param fnn - array containing node numbers of first nodes
   
   JK, 6.8.2008
*/
void gedge::give_first_node_numbers (long *fnn)
{
  fnn[0]=nlist[0];
  fnn[1]=nlist[1];
}

/**
   function assembles node numbers of last nodes
   
   @param lnn - array containing node numbers of last nodes
   
   JK, 6.8.2008
*/
void gedge::give_last_node_numbers (long *lnn)
{
  lnn[0]=nlist[2];
  lnn[1]=nlist[3];
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
void gedge::give_mult_code_numbers (long fln,long *mcn)
{
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
}
