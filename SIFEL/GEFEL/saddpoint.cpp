#include "saddpoint.h"

saddpoint::saddpoint ()
{
  //  number of primary unknowns
  nrdof=0;
  //  number of dual unknowns (Lagrange multipliers)
  nm=0;
  
  //  allocation of dense matrix
  dm = new densemat ();
  
  c=NULL;
}

saddpoint::~saddpoint ()
{
  delete [] c;
  delete dm;
}




void saddpoint::initiate (seqselnodes */*selnodfeti*/,gtopology *top,FILE *out)
{
  long i,j,k,l,ndofn,dof,nid;
  long *red;
  red = new long [ns+1];
  
  //  number of DOFs (unknowns) in coarse problem = total number of boundary DOFs
  //ndofcp =selnodfeti->tndof;
  
  ncdofd = new long [ns];
  for (i=0;i<ns;i++){
    //ncdofd[i] = selnodfeti->ndofdom[i];
  }
  
  
  //  array containing code numbers contributing to the coarse problem
  //edofs = selnodfeti->ldof;
  edofs = new long* [ns];
  for (i=0;i<ns;i++){
    edofs[i] = new long [ncdofd[i]];
  }
  
  for (i=0;i<ns;i++){
    for (j=0;j<ncdofd[i];j++){
      //edofs[i][j] = selnodfeti->ldof[i][j];
    }
  }

  red[0]=0;
  nid=0;
  for (i=0;i<ns;i++){
    red[i+1]=0;
    for (j=0;j<top->stop->nnsd[i];j++){
      ndofn = top->give_ndofn (nid);
      for (k=0;k<ndofn;k++){
	dof = top->give_dof (nid,k);
	if (dof>red[i+1]){
	  red[i+1]=dof;
	}
      }
      nid++;
    }
  }
  
  for (i=0;i<ns;i++){
    for (j=0;j<ncdofd[i];j++){
      if (edofs[i][j]>0)
	edofs[i][j] -= red[i];
      else
	edofs[i][j] += red[i];
    }
  }

  fprintf (out,"\n\n\n osel \n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n");
    for (j=0;j<ncdofd[i];j++){
      fprintf (out,"\n edofs %ld",edofs[i][j]);
    }
  }
  fprintf (out,"\n\n\n osel \n");








  ccn = new long* [ns];
  for (i=0;i<ns;i++){
    ccn[i] = new long [ncdofd[i]];
  }
  
  for (i=0;i<ns;i++){
    for (j=0;j<ncdofd[i];j++){
      //ccn[i][j] = selnodfeti->cnmas[i][j];
    }
  }
  
  delete [] red;
  
  




  nsid = new long [top->nn];
  l=0;
  for (i=0;i<ns;i++){
    for (j=0;j<top->stop->nnsd[i];j++){
      nsid[l]=i;
      l++;
    }
  }

  //  number of DOFs on subdomains
  ndofdom = new long [ns];
  for (i=0;i<ns;i++){
    ndofdom[i]=0;
  }
  
  for (i=0;i<top->nn;i++){
    ndofn=top->give_ndofn (i);
    l=nsid[i];
    for (j=0;j<ndofn;j++){
      k=top->give_dof (i,j);
      if (k>0)
	ndofdom[l]++;
    }
  }
  
  cndom = new long* [ns];
  for (i=0;i<ns;i++){
    cndom[i] = new long [ndofdom[i]];
    ndofdom[i]=0;
  }
  
  for (i=0;i<top->nn;i++){
    ndofn=top->give_ndofn (i);
    l=nsid[i];
    for (j=0;j<ndofn;j++){
      k=top->give_dof (i,j);
      if (k>0){
	cndom[l][ndofdom[l]]=k;
	ndofdom[l]++;
      }
    }
  }

  
}

/*
void saddpoint::matrix_e ()
{
  long i,j,k,l,m;
  
  e = new double [nrdof*ndofcp];
  for (i=0;i<nrdof*ndofcp;i++){
    e[i]=0.0;
  }
  
  for (l=0;l<ns;l++){
    for (i=0;i<ncdofd[l];i++){
      j=edofs[l][i];
      m=cndom[l][j];
      k=ccn[l][i]-1;
      if (j<0){
	e[k*nrdof+m-1]=-1.0;
      }
      else{
	e[k*nrdof+m-1]=1.0;
      }
    }
  }
}
*/

void saddpoint::matrix_e ()
{

}


void saddpoint::get_jumps (long nc,double *jumps)
{
  long i;
  
  if (c==NULL)
    c = new double [nc];
  
  for (i=0;i<nc;i++){
    c[i]=jumps[i];
    //c[i]=-5.0;
    //c[i]=0.0;
    //c[i]=-2.0e-3;
  }
}

/**
   function solves system of equations
   
   @param top - general topology
   @param gm - general %matrix
   @param lhs - array of solution
   @param rhs - %vector of the right hand side
   @param out - output file
   
   JK, 28.11.2007
*/
void saddpoint::solve_system (gtopology *top,gmatrix *gm,double *lhs,double *rhs,FILE *out)
{
  long i,j;
  double limit=1.0e-7;
  double *condmat,*condvect,*x,*y;
  
  //  number of reduced/boundary/interface DOFs
  nrdof=top->nbdof;
  //  number of internal DOFs
  nidof=top->nidof;

  //  array for condensed matrix
  condmat = new double [nrdof*nrdof];
  //  array for condensed vector
  condvect = new double [nrdof];
  
  //  condensation of internal DOFs
  gm->condense (top,condmat,condvect,lhs,rhs,nrdof,1,out);
  
  
  //  number of unknowns in the final system
  n=nrdof+nrdof/2;
  
  //  number of multipliers
  nm=nrdof/2;
  
  //  allocation of internal arrays of the object dm
  dm->alloc (n);
  x = new double [n];
  y = new double [n];
  
  for (i=0;i<nrdof;i++){
    for (j=0;j<nrdof;j++){
      dm->a[i*n+j]=condmat[i*nrdof+j];
    }
  }

  
  for (i=0;i<nm;i++){
    dm->a[i*n+nrdof+i]=1.0;
  }
  j=nm;
  for (i=0;i<nm;i++){
    dm->a[j*n+nrdof+i]=-1.0;
    j++;
  }
  
  j=nrdof;
  for (i=0;i<nm;i++){
    dm->a[j*n+i]=1.0;
    j++;
  }
  
  j=nrdof;
  for (i=0;i<nm;i++){
    dm->a[j*n+nm+i]=-1.0;
    j++;
  }
  
  j=nidof;
  for (i=0;i<nrdof;i++){
    y[i]=rhs[j];
    j++;
  }
  
  j=nrdof;
  for (i=0;i<nm;i++){
    y[j]=c[i];
    j++;
  }

  //  solution of the final system of equations
  dm->gemp (x,y,1,limit,1);
  
  gm->condense (top,condmat,x,lhs,rhs,nrdof,2,out);
  
  //  deallocation of the matrix
  dm->dealloc();
  delete [] x;
  delete [] y;
  delete [] condmat;
  delete [] condvect;
}
