#include "psolver.h"
#include <string.h>
#include <math.h>

psolver::psolver (int np,int mr,int nd)
{
  nproc=np;  myrank=mr;  ndom=nd;

  tdd = (ldomdectype) 0;
  rssol = (lredsystsolver) 0;
  
  ngdof=0;  ndof=0;  indof=0;
  nrbm=0;  enrbm=0;  lithr=0.0;
  
  nicg=0;  anicg=0;  errcg=0.0;  aerrcg=0.0;
  
  hsize=0;  limit=0.0;  zero=0.0;
}

psolver::~psolver()
{

}

void psolver::movedata (paral *plg)
{
  ndof=plg->ndof;
  indof=plg->indof;
}


/**
   function solves reduced system of equations from primal
   domain decomposition by LDL decomposition
   matrix of reduced system is assembled on master processor
   in skyline storage
   
   condmat - array containing reduced matrix from one subdomain
   condvect - array containing reduced vector of right hand side
   
   2.12.2001
*/
void psolver::redsys_parldl (paral *plg,double *condmat,double *condvect,FILE *out)
{
  int ind;
  long i,j,n,maxnrdof,sizebuff;
  skyline rsm_sky;
  double *rhs,*lhs;
  double *pole;
  //  char *buff;
  MPI_Status stat;
  
  n=ndof-indof;
  maxnrdof = plg->maxnrdof;
  if (myrank==0)  ngdof=plg->ngdof;
  
  sizebuff = (maxnrdof*maxnrdof+maxnrdof)*sizeof(double);
  //buff = new char [sizebuff];

  pole = new double[maxnrdof*maxnrdof+maxnrdof];


  if (myrank==0){
    rsm_sky.allocadr (ngdof);
    lhs = new double [ngdof];
    memset (lhs,0,ngdof*sizeof(double));
    rhs = new double [ngdof];
    memset (rhs,0,ngdof*sizeof(double));


    for (i=0;i<nproc;i++){
      rsm_sky.column_lengths_elem (plg->masgcn[i],plg->nrdofdom[i]);
    }
    rsm_sky.addresses ();
    rsm_sky.neglobmat ();
    rsm_sky.allocglomat ();
    
    rsm_sky.localized (condmat,plg->masgcn[ndom],plg->nrdofdom[ndom]);
    locglob (rhs,condvect,plg->masgcn[ndom],plg->nrdofdom[ndom]);
    

    for (i=1;i<nproc;i++){
      //MPI_Recv(buff,sizebuff,MPI_BYTE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      MPI_Recv(pole,maxnrdof*maxnrdof+maxnrdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

      
      j=stat.MPI_TAG;  ind=0;
      //MPI_Unpack(buff,sizebuff,&ind,condmat,plg->nrdofdom[j]*plg->nrdofdom[j],MPI_DOUBLE,MPI_COMM_WORLD);
      //MPI_Unpack(buff,sizebuff,&ind,condvect,plg->nrdofdom[j],MPI_DOUBLE,MPI_COMM_WORLD);


      long ijk=0;
      for (long ij=0;ij<plg->nrdofdom[j];ij++){
	for (long ji=0;ji<plg->nrdofdom[j];ji++){
	  condmat[ij*plg->nrdofdom[j]+ji]=pole[ijk];
	  ijk++;
	}
      }
      for (long ij=0;ij<plg->nrdofdom[j];ij++){
	condvect[ij]=pole[ijk];
	ijk++;
      }
      
      
      rsm_sky.localized (condmat,plg->masgcn[j],plg->nrdofdom[j]);
      locglob (rhs,condvect,plg->masgcn[j],plg->nrdofdom[j]);
      
    }
    

    rsm_sky.ldl_sky (lhs,rhs,zero,1);
    
    
    fprintf (out,"\n\n\n SOLUTION OF REDUCED SYSTEM OF EQUATIONS \n");
    for (i=0;i<ngdof;i++){
      fprintf (out,"\n lhs %4ld   %lf",i,lhs[i]);
    }
    fprintf (out,"\n\n\n");
    

    for (i=0;i<nproc;i++){
      j=plg->domproc[i];

      if (j==myrank)  continue;
      nullvr (condvect,plg->nrdofdom[i]);
      globloc (lhs,condvect,plg->masgcn[i],plg->nrdofdom[i]);
      ind=0;
      //MPI_Pack(condvect,plg->nrdofdom[i],MPI_DOUBLE,buff,sizebuff,&ind,MPI_COMM_WORLD);
      //MPI_Send(buff,sizebuff,MPI_BYTE,j,ndom,MPI_COMM_WORLD);
      
      for (long ij=0;ij<plg->nrdofdom[i];ij++){
	pole[ij]=condvect[ij];
      }
      
      MPI_Send(pole,maxnrdof*maxnrdof+maxnrdof,MPI_DOUBLE,j,ndom,MPI_COMM_WORLD);

    }
    
    nullvr (condvect,plg->nrdofdom[ndom]);
    globloc (lhs,condvect,plg->masgcn[ndom],plg->nrdofdom[ndom]);
    
    delete [] lhs;  delete [] rhs;
  }
  else{
    
    long ijk=0;
    for (long ij=0;ij<n;ij++){
      for (long ji=0;ji<n;ji++){
	pole[ijk]=condmat[ij*n+ji];
	ijk++;
      }
    }
    for (long ij=0;ij<n;ij++){
      pole[ijk]=condvect[ij];
      ijk++;
    }

    
    //ind=0;
    //MPI_Pack(condmat,n*n,MPI_DOUBLE,buff,sizebuff,&ind,MPI_COMM_WORLD);
    //MPI_Pack(condvect,n,MPI_DOUBLE,buff,sizebuff,&ind,MPI_COMM_WORLD);
    //MPI_Send(buff,sizebuff,MPI_BYTE,0,ndom,MPI_COMM_WORLD);
    MPI_Send(pole,maxnrdof*maxnrdof+maxnrdof,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);

    //MPI_Recv(buff,sizebuff,MPI_BYTE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    MPI_Recv(pole,maxnrdof*maxnrdof+maxnrdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

    for (long ij=0;ij<n;ij++){
      condvect[ij]=pole[ij];
    }
    
    //ind=0;
    //MPI_Unpack(buff,sizebuff,&ind,condvect,n,MPI_DOUBLE,MPI_COMM_WORLD);
    
  }


  //delete [] buff;
  
}





/**
   function solves reduced system of equations from primal
   domain decomposition by LDL decomposition
   matrix of reduced system is assembled on master processor
   in skyline storage
   
   condmat - array containing reduced matrix from one subdomain
   condvect - array containing reduced vector of right hand side
   
   2.12.2001
*/
void psolver::redsys_parlu (paral *plg,double *condmat,double *condvect,FILE *out)
{
  int ind;
  long i,j,n,maxnrdof,sizebuff;
  dskyline rsm_dsky;
  double *rhs,*lhs;
  double *pole;
  //char *buff;
  MPI_Status stat;
  
  n=ndof-indof;
  maxnrdof = plg->maxnrdof;
  if (myrank==0)  ngdof=plg->ngdof;
  
  sizebuff = (maxnrdof*maxnrdof+maxnrdof)*sizeof(double);
  //buff = new char [sizebuff];
  
  pole = new double[maxnrdof*maxnrdof+maxnrdof];

  if (myrank==0){
    rsm_dsky.allocadr (ngdof);
    lhs = new double [ngdof];
    memset (lhs,0,ngdof*sizeof(double));
    rhs = new double [ngdof];
    memset (rhs,0,ngdof*sizeof(double));


    for (i=0;i<nproc;i++){
      rsm_dsky.column_lengths_elem (plg->masgcn[i],plg->nrdofdom[i]);
    }
    rsm_dsky.addresses ();
    rsm_dsky.neglobmat ();
    rsm_dsky.allocglomat ();
    
    rsm_dsky.localized (condmat,plg->masgcn[ndom],plg->nrdofdom[ndom]);
    locglob (rhs,condvect,plg->masgcn[ndom],plg->nrdofdom[ndom]);
    

    for (i=1;i<nproc;i++){
      //MPI_Recv(buff,sizebuff,MPI_BYTE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      MPI_Recv(pole,maxnrdof*maxnrdof+maxnrdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      j=stat.MPI_TAG;  ind=0;
      //MPI_Unpack(buff,sizebuff,&ind,condmat,plg->nrdofdom[j]*plg->nrdofdom[j],MPI_DOUBLE,MPI_COMM_WORLD);
      //MPI_Unpack(buff,sizebuff,&ind,condvect,plg->nrdofdom[j],MPI_DOUBLE,MPI_COMM_WORLD);
      
     long ijk=0;
      for (long ij=0;ij<plg->nrdofdom[j];ij++){
	for (long ji=0;ji<plg->nrdofdom[j];ji++){
	  condmat[ij*plg->nrdofdom[j]+ji]=pole[ijk];
	  ijk++;
	}
      }
      for (long ij=0;ij<plg->nrdofdom[j];ij++){
	condvect[ij]=pole[ijk];
	ijk++;
      }
 
      rsm_dsky.localized (condmat,plg->masgcn[j],plg->nrdofdom[j]);
      locglob (rhs,condvect,plg->masgcn[j],plg->nrdofdom[j]);
      
    }
    

    rsm_dsky.lu_dsky (lhs,rhs,zero,1);
    
    
    fprintf (out,"\n\n\n SOLUTION OF REDUCED SYSTEM OF EQUATIONS \n");
    for (i=0;i<ngdof;i++){
      fprintf (out,"\n lhs   %le",lhs[i]);
    }
    fprintf (out,"\n\n\n");
    

    for (i=0;i<nproc;i++){
      j=plg->domproc[i];

      if (j==myrank)  continue;
      nullvr (condvect,plg->nrdofdom[i]);
      globloc (lhs,condvect,plg->masgcn[i],plg->nrdofdom[i]);
      ind=0;
      //MPI_Pack(condvect,plg->nrdofdom[i],MPI_DOUBLE,buff,sizebuff,&ind,MPI_COMM_WORLD);
      //MPI_Send(buff,sizebuff,MPI_BYTE,j,ndom,MPI_COMM_WORLD);

      for (long ij=0;ij<plg->nrdofdom[i];ij++){
	pole[ij]=condvect[ij];
      }

      MPI_Send(pole,maxnrdof*maxnrdof+maxnrdof,MPI_DOUBLE,j,ndom,MPI_COMM_WORLD);

    }
    
    nullvr (condvect,plg->nrdofdom[ndom]);
    globloc (lhs,condvect,plg->masgcn[ndom],plg->nrdofdom[ndom]);
    
    delete [] lhs;  delete [] rhs;
  }
  else{

    long ijk=0;
    for (long ij=0;ij<n;ij++){
      for (long ji=0;ji<n;ji++){
	pole[ijk]=condmat[ij*n+ji];
	ijk++;
      }
    }
    for (long ij=0;ij<n;ij++){
      pole[ijk]=condvect[ij];
      ijk++;
    }

    //ind=0;
    //MPI_Pack(condmat,n*n,MPI_DOUBLE,buff,sizebuff,&ind,MPI_COMM_WORLD);
    //MPI_Pack(condvect,n,MPI_DOUBLE,buff,sizebuff,&ind,MPI_COMM_WORLD);
    //MPI_Send(buff,sizebuff,MPI_BYTE,0,ndom,MPI_COMM_WORLD);

    MPI_Send(pole,maxnrdof*maxnrdof+maxnrdof,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);

    //MPI_Recv(buff,sizebuff,MPI_BYTE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    MPI_Recv(pole,maxnrdof*maxnrdof+maxnrdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

    for (long ij=0;ij<n;ij++){
      condvect[ij]=pole[ij];
    }

    //ind=0;
    //MPI_Unpack(buff,sizebuff,&ind,condvect,n,MPI_DOUBLE,MPI_COMM_WORLD);

  }


  //delete [] buff;
  
}

/**
   function solves reduced system of equations from primal
   domain decomposition by Gauss elimination
   matrix of reduced system is assembled on master processor
   in skyline storage
   
   condmat - array containing reduced matrix from one subdomain
   condvect - array containing reduced vector of right hand side
   
   2.5.2002
*/
void psolver::redsys_pargemp (paral *plg,double *condmat,double *condvect,FILE *out)
{
  int ind;
  long i,j,n,maxnrdof,sizebuff;
  densemat rsm_dm;
  double *rhs,*lhs;
  double *pole;
  //char *buff;
  MPI_Status stat;
  
  n=ndof-indof;
  maxnrdof = plg->maxnrdof;
  if (myrank==0)  ngdof=plg->ngdof;
  
  sizebuff = (maxnrdof*maxnrdof+maxnrdof)*sizeof(double);
  //buff = new char [sizebuff];
  
  pole = new double[maxnrdof*maxnrdof+maxnrdof];

  if (myrank==0){
    rsm_dm.alloc (ngdof);
    lhs = new double [ngdof];
    memset (lhs,0,ngdof*sizeof(double));
    rhs = new double [ngdof];
    memset (rhs,0,ngdof*sizeof(double));


    rsm_dm.localized (condmat,plg->masgcn[ndom],plg->nrdofdom[ndom]);
    locglob (rhs,condvect,plg->masgcn[ndom],plg->nrdofdom[ndom]);
    

    for (i=1;i<nproc;i++){
      //MPI_Recv(buff,sizebuff,MPI_BYTE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      MPI_Recv(pole,maxnrdof*maxnrdof+maxnrdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      //j=stat.MPI_TAG;  ind=0;
      //MPI_Unpack(buff,sizebuff,&ind,condmat,plg->nrdofdom[j]*plg->nrdofdom[j],MPI_DOUBLE,MPI_COMM_WORLD);
      //MPI_Unpack(buff,sizebuff,&ind,condvect,plg->nrdofdom[j],MPI_DOUBLE,MPI_COMM_WORLD);
      
      long ijk=0;
      for (long ij=0;ij<plg->nrdofdom[j];ij++){
	for (long ji=0;ji<plg->nrdofdom[j];ji++){
	  condmat[ij*plg->nrdofdom[j]+ji]=pole[ijk];
	  ijk++;
	}
      }
      for (long ij=0;ij<plg->nrdofdom[j];ij++){
	condvect[ij]=pole[ijk];
	ijk++;
      }

      rsm_dm.localized (condmat,plg->masgcn[j],plg->nrdofdom[j]);
      locglob (rhs,condvect,plg->masgcn[j],plg->nrdofdom[j]);
      
    }
    

    rsm_dm.gemp (lhs,rhs,1,zero,1);
    
    
    fprintf (out,"\n\n\n SOLUTION OF REDUCED SYSTEM OF EQUATIONS \n");
    for (i=0;i<ngdof;i++){
      fprintf (out,"\n lhs   %le",lhs[i]);
    }
    fprintf (out,"\n\n\n");
    

    for (i=0;i<nproc;i++){
      j=plg->domproc[i];

      if (j==myrank)  continue;
      nullvr (condvect,plg->nrdofdom[i]);
      globloc (lhs,condvect,plg->masgcn[i],plg->nrdofdom[i]);
      ind=0;
      //MPI_Pack(condvect,plg->nrdofdom[i],MPI_DOUBLE,buff,sizebuff,&ind,MPI_COMM_WORLD);
      //MPI_Send(buff,sizebuff,MPI_BYTE,j,ndom,MPI_COMM_WORLD);

      for (long ij=0;ij<plg->nrdofdom[i];ij++){
	pole[ij]=condvect[ij];
      }
      
      MPI_Send(pole,maxnrdof*maxnrdof+maxnrdof,MPI_DOUBLE,j,ndom,MPI_COMM_WORLD);

    }
    
    nullvr (condvect,plg->nrdofdom[ndom]);
    globloc (lhs,condvect,plg->masgcn[ndom],plg->nrdofdom[ndom]);
    
    delete [] lhs;  delete [] rhs;
  }
  else{

    long ijk=0;
    for (long ij=0;ij<n;ij++){
      for (long ji=0;ji<n;ji++){
	pole[ijk]=condmat[ij*n+ji];
	ijk++;
      }
    }
    for (long ij=0;ij<n;ij++){
      pole[ijk]=condvect[ij];
      ijk++;
    }
    
    //ind=0;
    //MPI_Pack(condmat,n*n,MPI_DOUBLE,buff,sizebuff,&ind,MPI_COMM_WORLD);
    //MPI_Pack(condvect,n,MPI_DOUBLE,buff,sizebuff,&ind,MPI_COMM_WORLD);
    //MPI_Send(buff,sizebuff,MPI_BYTE,0,ndom,MPI_COMM_WORLD);
    MPI_Send(pole,maxnrdof*maxnrdof+maxnrdof,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);

    //MPI_Recv(buff,sizebuff,MPI_BYTE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

    MPI_Recv(pole,maxnrdof*maxnrdof+maxnrdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

    for (long ij=0;ij<n;ij++){
      condvect[ij]=pole[ij];
    }

    //ind=0;
    //MPI_Unpack(buff,sizebuff,&ind,condvect,n,MPI_DOUBLE,MPI_COMM_WORLD);

  }


  //delete [] buff;
  
}


/**
   function solves reduced system of equations from primal
   domain decomposition by conjugate gradient method
   matrix of reduced system is not assembled
   iteration process is directed by master processor
   
   condmat - array containing reduced matrix from one subdomain
   condvect - array containing reduced vector of right hand side
   
   2.12.2001
*/
void psolver::redsys_parcg (paral *plg,double *condmat,double *condvect,FILE *out)
{
  int ind;
  long i,j,k,n,maxnrdof,sizebuff;
  double alpha,beta,nom,denom,norrhs;
  double *lhs,*rhs,*buff,*p,*d,*r;
  MPI_Status stat;
  
  n=ndof-indof;
  maxnrdof = plg->maxnrdof;
  if (myrank==0)  ngdof=plg->ngdof;

    sizebuff = maxnrdof+1;
    buff = new double [sizebuff];
    memset (buff,0,sizebuff*sizeof(double));
    
    //  right hand side assembling
    if (myrank==0){
      lhs = new double [ngdof];
      memset (lhs,0,ngdof*sizeof(double));
      rhs = new double [ngdof];
      memset (rhs,0,ngdof*sizeof(double));
      
      locglob (rhs,condvect,plg->masgcn[ndom],n);
      for (i=1;i<nproc;i++){
	MPI_Recv(buff,sizebuff,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	j=stat.MPI_TAG;  ind=0;
	locglob (rhs,buff,plg->masgcn[j],plg->nrdofdom[j]);
      }
    }
    else{
      copydarr (condvect,buff,n);
      MPI_Send(buff,sizebuff,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);
    }
    
    //  iteration loop block
    if (myrank==0){
      r = new double [ngdof];
      d = new double [ngdof];
      p = new double [ngdof];

      //  iteration initialization
      norrhs=0.0;  nom=0.0;
      for (i=0;i<ngdof;i++){
	lhs[i]=0.0;
	r[i]=rhs[i];
	d[i]=rhs[i];
	norrhs+=rhs[i]*rhs[i];
	nom+=r[i]*r[i];
      }
      
      //  iteration loop
      for (i=0;i<nicg;i++){
	
	//  direction vector scattering
	for (j=0;j<nproc;j++){
	  k=plg->domproc[j];
	  if (k==myrank)  continue;
	  memset (buff,0,sizebuff*sizeof(double));
	  globloc (d,buff,plg->masgcn[j],plg->nrdofdom[j]);
	  MPI_Send(buff,sizebuff,MPI_DOUBLE,k,ndom,MPI_COMM_WORLD);
	}
	memset (buff,0,sizebuff*sizeof(double));
	globloc (d,buff,plg->masgcn[ndom],plg->nrdofdom[ndom]);
	
	//  matrix-vector multiplication
	mv (condmat,buff,condvect,n,n);
	
	//  matrix-vector multiplication collection
	memset (p,0,ngdof*sizeof(double));
	locglob (p,condvect,plg->masgcn[ndom],n);
	for (j=1;j<nproc;j++){
	  MPI_Recv(buff,sizebuff,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	  k=stat.MPI_TAG;
	  locglob (p,buff,plg->masgcn[k],plg->nrdofdom[k]);
	}
	
	//  denominator
	denom=ss(d,p,ngdof);
	
	if (fabs(denom)<zero){
	  if (i!=nicg-1){
	    buff[maxnrdof]=1.0;
	    for (j=1;j<nproc;j++){
	      MPI_Send(buff,sizebuff,MPI_DOUBLE,j,ndom,MPI_COMM_WORLD);
	    }
	  }
	  fprintf (stderr,"\n\n denominator in alpha expression in parallel conjugate gradient method is equal to zero.\n");
	  break;
	}
	
	//  coefficient alpha
	alpha=nom/denom;
	
	//  new global vectors
	for (j=0;j<ngdof;j++){
	  r[j]-=alpha*p[j];
	  lhs[j]+=alpha*d[j];
	}
	
	
	denom=nom;
	nom=ss(r,r,ngdof);
	
	fprintf (stdout,"\n iteration %ld    norres/norrhs %le",i,nom/norrhs);
	fprintf (out,"\n iteration %ld    norres/norrhs %le",i,nom/norrhs);
	if (nom/norrhs<errcg){
	  if (i!=nicg-1){
	    buff[maxnrdof]=1.0;
	    for (j=1;j<nproc;j++){
	      MPI_Send(buff,sizebuff,MPI_DOUBLE,j,ndom,MPI_COMM_WORLD);
	    }
	  }
	  break;
	}
	
	//  beta coefficient
	beta=nom/denom;
	
	//  new direction vector
	for (j=0;j<ngdof;j++){
	  d[j]=r[j]+beta*d[j];
	}

      }
      anicg=i;
      aerrcg=nom/norrhs;
      
      fprintf (out,"\n\n\n REDUCED SYSTEM SOLUTION (on master)\n");
      for (i=0;i<ngdof;i++){
	fprintf (out,"\n lhs %ld     %le",i,lhs[i]);
      }

    }
    else{
      //  iteration loop
      for (i=0;i<nicg;i++){
	//  direction vector receipt
	MPI_Recv(buff,sizebuff,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	
	if (buff[maxnrdof]>0.5)  break;
	
	//  matrix-vector multiplication
	mv (condmat,buff,condvect,n,n);
	
	//  matrix-vector multiplication collection
	copydarr (condvect,buff,n);
	MPI_Send(buff,sizebuff,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);
	
      }
      anicg=i;
    }
    
    //  global solution distribution
    if (myrank==0){
      for (i=0;i<nproc;i++){
	j=plg->domproc[i];
	if (j==myrank)  continue;
	memset (buff,0,sizebuff*sizeof(double));
	globloc (lhs,buff,plg->masgcn[i],plg->nrdofdom[i]);
	MPI_Send(buff,sizebuff,MPI_DOUBLE,j,ndom,MPI_COMM_WORLD);
      }
      nullvr (condvect,plg->nrdofdom[ndom]);
      globloc (lhs,condvect,plg->masgcn[ndom],plg->nrdofdom[ndom]);

    }
    else{
      MPI_Recv(buff,sizebuff,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      copydarr (buff,condvect,n);
    }
    
    delete [] buff;

}
















//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************

/**
   function computes size of matrix H
   
   @param rbm - array containing rigid body motions
   @param maxnrbm - maximum number of rigid body motions
   @param ngdof - number of global DOFs = number of Lagrange multipliers
   
   2.12.2001
*/
void psolver::hmatrixsize (paral *plg,double *rbm,long &maxnrbm)
{
  long i;
  long *ibuff;
  MPI_Status stat;
  
  ibuff = new long [3];
  ibuff[0]=nrbm;
  ibuff[1]=0;
  ibuff[2]=0;
  
  if (myrank==0){
    plg->nrbmdom = new long [nproc];
    plg->rbmadr = new long [nproc+1];
    for (i=0;i<nproc+1;i++){
      plg->rbmadr[i]=0;
    }
    maxnrbm=nrbm;  plg->nrbmdom[ndom]=nrbm;
    for (i=1;i<nproc;i++){
      MPI_Recv(ibuff,3,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      if (maxnrbm<ibuff[0])  maxnrbm=ibuff[0];
      plg->nrbmdom[stat.MPI_TAG]=ibuff[0];
    }
    
    plg->rbmadr[0]=0;
    for (i=1;i<nproc;i++){
      plg->rbmadr[i+1]=plg->rbmadr[i]+plg->nrbmdom[i];
      
      //fprintf (stderr,"\n i %ld    nrbmdom %ld      rbmadr %ld",i,plg->nrbmdom[i],plg->rbmadr[i+1]);


    }
    
    ibuff[0]=maxnrbm;
    ibuff[1]=plg->ngdof;
    ibuff[2]=plg->rbmadr[nproc];
    for (i=1;i<nproc;i++){
      MPI_Send(ibuff,3,MPI_LONG,i,ndom,MPI_COMM_WORLD);
    }
    
    ngdof=plg->ngdof;
    hsize=plg->rbmadr[nproc];
  }
  else{
    MPI_Send(ibuff,3,MPI_LONG,0,ndom,MPI_COMM_WORLD);
    MPI_Recv(ibuff,3,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    maxnrbm=ibuff[0];  ngdof=ibuff[1];  hsize=ibuff[2];
  }
  delete [] ibuff;
  
}


/**
   function assembles matrix H
   columns are RBM after localization
   
   @param h - matrix H
   @param rbm - array containing rigid body motions
   @param maxnrbm - maximum number of rigid body motions
   @param ngdof - number of global DOFs = number of Lagrange multipliers
   
   2.12.2001
*/
void psolver::hmatrix (gtopology *top,paral *plg,double *h,double *rbm,long maxnrbm)
{
  long i,j,k,l,m;
  double *buff;
  MPI_Status stat;
  
  buff = new double [maxnrbm*ngdof];
  for (i=0;i<maxnrbm*ngdof;i++){
    buff[i]=0.0;
  }


  /*
  fprintf (out,"\n\n\n kontrola RBM PRED LOCGLOBFETI");
  for (i=0;i<ndof;i++){
    fprintf (out,"\n");
    for (j=0;j<ppd->nrbm;j++){
      fprintf (out," %lf",rbm[j*ndof+i]);
    }
  }
  */


  for (i=0;i<nrbm;i++){
    plg->locglobfeti (top,buff+i*ngdof,rbm+i*ndof);
  }
  
  
  /*
  fprintf (out,"\n\n\n KONTROLA RBM PO LOCGLOBFETI");
  for (i=0;i<ngdof;i++){
    fprintf (out,"\n");
    for (j=0;j<ppd->nrbm;j++){
      fprintf (out," %lf",buff[j*ngdof+i]);
    }
  }
  */


  
  if (myrank==0){
    hsize = plg->rbmadr[nproc];
    

    m=0;
    for (j=0;j<plg->nrbmdom[ndom];j++){
      l=plg->rbmadr[ndom]+j;
      for (k=0;k<ngdof;k++){
	h[l]=0.0-buff[m];
	l+=hsize;  m++;
      }

      //fprintf (stderr,"\n konecne l je rovno %ld",l);

    }



    for (i=1;i<nproc;i++){
      MPI_Recv(buff,maxnrbm*ngdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      

      m=0;
      for (j=0;j<plg->nrbmdom[stat.MPI_TAG];j++){
	l=plg->rbmadr[stat.MPI_TAG]+j;
	for (k=0;k<ngdof;k++){
	  h[l]=0.0-buff[m];
	  l+=hsize;  m++;
	}

	//fprintf (stderr,"\n j %ld    konecne l je rovno %ld",j,l);

      }

    }
  }
  else{
    MPI_Send(buff,maxnrbm*ngdof,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);
  }

  MPI_Barrier (MPI_COMM_WORLD);

}



/**
   function assembles vector q
   
   @param q - vector q
   @param rbm - array containing rigid body motions
   @param f - right hand side
   @param maxnrbm - maximum number of rigid body motions
   
   6.3.2002
*/
void psolver::qvector (paral *plg,double *q,double *rbm,double *f,long maxnrbm)
{
  long i,j,l,m;
  double *buff;
  MPI_Status stat;
  

  buff = new double [maxnrbm];
  
  for (i=0;i<nrbm;i++){
    buff[i]=0.0-ss(rbm+i*ndof,f,ndof);
  }
  

  if (myrank==0){
    

    l=plg->rbmadr[ndom];  m=0;
    for (j=0;j<plg->nrbmdom[ndom];j++){
      q[l]=buff[m];  l++;  m++;
    }

    
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,maxnrbm,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

      

      l=plg->rbmadr[stat.MPI_TAG];  m=0;
      for (j=0;j<plg->nrbmdom[stat.MPI_TAG];j++){
	q[l]=buff[m];  l++;  m++;
      }

      
    }
    
  }
  else{
    MPI_Send(buff,maxnrbm,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);
  }
  

  delete [] buff;

}

/**
   function computes projection in FETI method
   v_new = v_old - H . (H^T . H)^{-1} . H^T . v_old
   function overwrites vector v by projected vector
   
   @param v - %vector
   @param h - matrix H
   @param h1 - matrix (H^T.H)^{-1}
   @param ngdof,hsize - number of rows and columns of matrix H
   
   ngdof - number of Lagrange multipliers
   hsize - total number of all rigid body motions
   
   JK, 7.3.2002
*/
void psolver::feti_projection (double *v,double *h,double *h1)
{
  long i;
  double *p,*q,*w;
  
  p = new double [hsize];
  q = new double [hsize];
  w = new double [ngdof];

  //  p = H^T.v
  //mtvc (h,v,p,ngdof,hsize);
  mtv (h,v,p,ngdof,hsize);
  
  
  /*
  fprintf (out,"\n\n\n kontrola H^T.v\n");
  for (i=0;i<hsize;i++){
    fprintf (out," %lf",p[i]);
  }
  */

  
  //  q = (H^T.H)^{-1}.p
  mv (h1,p,q,hsize,hsize);


  /*
  fprintf (out,"\n\n\n kontrola (H^T.H)^{-1}.H^T.v\n");
  for (i=0;i<hsize;i++){
    fprintf (out," %lf",q[i]);
  }
  */


  //  w = H.q
  mv (h,q,w,ngdof,hsize);
  
  
  /*
  fprintf (out,"\n\n\n kontrola H.(H^T.H)^{-1}.H^T.v\n");
  for (i=0;i<ngdof;i++){
    fprintf (out," %lf",w[i]);
  }
  */
  
  //  v_new = v_old - w
  for (i=0;i<ngdof;i++){
    v[i]-=w[i];
  }
  
  delete [] w;  delete [] q;  delete [] p;
}




/**
   function performs modified conjugate gradient method
   
   @param w - vector of Lagrange multipliers
   @param h - matrix H
   @param h1 - matrix (H^T.H)^{-1}
   @param ngdof,hsize - number of rows and columns of matrix H
   
   ngdof - number of Lagrange multipliers
   hsize - total number of all rigid body motions
   
   JK, 7.3.2002
*/
void psolver::mpcg (gtopology *top,paral *plg,gmatrix *gm,
		    double *w,double *rhs,double *q,double *h,double *h1,long *rbmi,FILE *out)
{
  long i,j;
  double nom,denom,alpha,beta;
  double *d,*dd,*g,*p,*pp,*buff;
  MPI_Status stat;
  
  //  direction vector allocation
  d = new double [ngdof+1];

  d[ngdof]=0.0;

  dd = new double [ndof];
  //  residuum vector allocation
  g = new double [ngdof];
  //  auxiliary vector allocation
  p = new double [ngdof];
  pp = new double [ndof];
  
  
  /*
  fprintf (stdout,"\n\n pocet iteraci %ld",ni);
  fprintf (out,"\n\n\n\n kontrola RHS\n");
  for (i=0;i<ndof;i++){
    fprintf (out," %lf",rhs[i]);
  }

  fprintf (out,"\n\n\n kontrola indexu zavislych rovnic");
  for (i=0;i<ppd->enrbm;i++){
    fprintf (out,"   %ld",rbmi[i]);
  }
  */
  
  /*
  fprintf (out,"\n\n\n kontrola matice tuhosti po eliminaci");
  for (i=0;i<ndof;i++){
    fprintf (out,"\n");
    for (j=Sm_sky->adr[i];j<Sm_sky->adr[i+1];j++){
      fprintf (out," %lf",Sm_sky->a[j]);
    }
  }
  */
  
  
  if (myrank==0){
    
    buff = new double [hsize];
    mv (h1,q,buff,hsize,hsize);
    mv (h,buff,w,ngdof,hsize);
    delete [] buff;
    
    /*
    fprintf (out,"\n\n\n kontrola vektoru q \n");
    for (i=0;i<hsize;i++){
      fprintf (out," %le",q[i]);
    }

    fprintf (out,"\n\n\n lagrangeovy multiplikatory pred ldl_feti_sky \n");
    for (i=0;i<ngdof;i++){
      fprintf (out," %le",w[i]);
    }
    */
    
    
    buff = new double [ngdof];
    
    for (j=1;j<nproc;j++){
      MPI_Send (w,ngdof,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
    }
    nullvr (dd,ndof);
    plg->globlocfeti (top,w,dd);
    subv (dd,rhs,ndof);
    
    
    /*
    fprintf (out,"\n\n\n inicializace --- kontrola vektoru prave strany dd pred funkci funkci ldl_feti_sky \n");
    for (i=0;i<ndof;i++){
      fprintf (out," %le",dd[i]);
    }
    */
    
    
    //Sm_sky->ldl_feti_sky (pp,dd,ppd->nrbm,rbmi,ppd->zero);
    gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
    
    
    /*
    fprintf (out,"\n\n\n inicializace --- kontrola vektoru reseni pp funkci ldl_feti_sky \n");
    for (i=0;i<ndof;i++){
      fprintf (out," %le",pp[i]);
    }
    */
    
    
    nullvr (g,ngdof);
    plg->locglobfeti (top,g,pp);
    for (j=1;j<nproc;j++){
      MPI_Recv(buff,ngdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      addv (g,buff,ngdof);
    }
    
    
    /*
    fprintf (out,"\n\n\n  inicializace    kontrola vektoru reziduii v inicializacni fazi pred projekci \n");
    for (i=0;i<ngdof;i++){
      fprintf (out," %le",g[i]);
    }
    */
    
    
    feti_projection (g,h,h1);
    
    
    /*
    fprintf (out,"\n\n\n  inicializace    kontrola vektoru reziduii v inicializacni fazi po projekci \n");
    for (i=0;i<ngdof;i++){
      fprintf (out," %le",g[i]);
    }
    */
    
    
    //  initialization of direction vector
    for (i=0;i<ngdof;i++){
      d[i]=0.0-g[i];
    }
    
    //  nominator evaluation
    nom = ss (g,g,ngdof);
    

    //fprintf (stdout,"\n\n nominator %lf",nom);

    
    /*************************/
    /*  main iteration loop  */
    /*************************/
    for (i=0;i<nicg;i++){
      
      //fprintf (stdout,"\n\n JEDEME GRADIENTY, proc %d, iterace %ld",myrank,i);

      /******************************************/
      /*  auxiliary vector computation K.d = p  */
      /******************************************/
      
      
      for (j=1;j<nproc;j++){
	MPI_Send (d,ngdof+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
      }
      

      //  computation of K^+ d on master
      nullvr (dd,ndof);
      plg->globlocfeti (top,d,dd);
      
      /*
      fprintf (out,"\n\n\n krok %ld        kontrola vektoru prave strany pred funkci ldl_feti_sky \n",i);
      for (j=0;j<ndof;j++){
	fprintf (out," %le",dd[j]);
      }
      
      fprintf (out,"\n\n\n krok %ld         kontrola vektoru reseni pp PRED  PRED funkci ldl_feti_sky \n",i);
      for (j=0;j<ndof;j++){
	fprintf (out," %le",pp[j]);
      }
      */
      
      
      //Sm_sky->ldl_feti_sky (pp,dd,ppd->nrbm,rbmi,ppd->zero);
      gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
      
      
      /*
      fprintf (out,"\n\n\n krok %ld         kontrola vektoru reseni pp funkci ldl_feti_sky \n",i);
      for (j=0;j<ndof;j++){
	fprintf (out," %le",pp[j]);
      }
      */
      
      
      nullvr (p,ngdof);
      plg->locglobfeti (top,p,pp);
      
      
      for (j=1;j<nproc;j++){
	MPI_Recv(buff,ngdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	addv (p,buff,ngdof);
      }
      
      
      /*
      fprintf (out,"\n\n\n krok %ld        kontrola pomocneho vektoru K.d\n",i);
      for (j=0;j<ngdof;j++){
	fprintf (out," %le",p[j]);
      }
      */
      /*
      fprintf (out,"\n\n\n krok %ld        kontrola pomocneho vektoru d\n",i);
      for (j=0;j<ngdof;j++){
	fprintf (out," %lf",d[j]);
      }
      */
      
      
      //  denominator of alpha
      denom = ss (d,p,ngdof);
      
      //fprintf (stdout,"\n denominator u alpha   %le",denom);

      fprintf (out,"\n\n kontrola citatele a jmenovatele pred alpha %lf / %lf",nom,denom);
      //fprintf (stderr,"\n\n kontrola citatele a jmenovatele pred alpha %lf / %lf",nom,denom);

      if (fabs(denom)<zero){
	fprintf (stderr,"\n\n zero denominator in modified conjugate gradient method (%s, line %d).\n",__FILE__,__LINE__);
	
	fprintf (stderr,"\n\n POKUD JSME TADY, TAK");
	
	if (i!=nicg-1){
	  d[ngdof]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send (d,ngdof+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	  }
	}
	break;
      }
      
      //fprintf (stderr," SEM SE NELZE DOSTAT\n");

      /*******************/
      /*  vypocet alpha  */
      /*******************/
      alpha = nom/denom;
      
      /**************************************************************/
      /*  vypocet noveho gradientu g a nove aproximace neznamych x  */
      /**************************************************************/
      for (j=0;j<ngdof;j++){
	w[j]+=alpha*d[j];
	g[j]+=alpha*p[j];
      }
      
      feti_projection (g,h,h1);
      
      
      //*******************************************************
      //*******************************************************
      //*******************************************************
      //  vlozka s predpodminenim
      //*******************************************************
      //*******************************************************
      //*******************************************************
      
      
      for (j=1;j<nproc;j++){
	MPI_Send (g,ngdof,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
      }

      nullvr (dd,ndof);
      plg->globlocfeti (top,g,dd);
      nullvr (pp,ndof);
      //lumpedprec (dd,pp);
      for (j=0;j<ndof;j++){
	pp[j]=dd[j];
      }
      nullvr (p,ngdof);
      plg->locglobfeti (top,p,pp);

      for (j=1;j<nproc;j++){
	MPI_Recv(buff,ngdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	addv (p,buff,ngdof);
      }

      feti_projection (p,h,h1);

      //*******************************************************
      //*******************************************************
      //*******************************************************
      //  konec vlozky s predpodminenim
      //*******************************************************
      //*******************************************************
      //*******************************************************


      denom = nom;
      if (fabs(denom)<zero){
	fprintf (stderr,"\n\n zero denominator of beta in function (%s, line %d).\n",__FILE__,__LINE__);
	
	if (i!=nicg-1){
	  d[ngdof]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send (d,ngdof+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	  }
	}
	break;
      }
      
      //nom = ss (g,g,ngdof);
      nom=ss(p,g,ngdof);
      
      fprintf (stdout,"\n iteration  %4ld        norm g   %le",i,nom);
      fprintf (out,"\n iteration  %4ld        norm g   %le",i,nom);
      
      fprintf (out,"\n kontrola citatele a jmenovatele pred betou  %lf / %lf",nom,denom);
      //fprintf (stderr,"\n kontrola citatele a jmenovatele pred betou  %lf / %lf",nom,denom);
      
      if (nom<errcg){
	if (i!=nicg-1){
	  d[ngdof]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send (d,ngdof+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	  }
	}
	break;
      }
      
      
      //  computation of beta coefficient
      beta = nom/denom;
      
      //  new direction vector
      for (j=0;j<ngdof;j++){
	//d[j]=beta*d[j]-g[j];
	d[j]=beta*d[j]-p[j];
      }
      
    }
    
    anicg=i;  aerrcg=nom;
    
    
    fprintf (out,"\n\n\n\n kontrola Lagrangeovych multiplikatoru \n");
    for (i=0;i<ngdof;i++){
      fprintf (out,"\n lagr. mult %4ld    %le",i,w[i]);
    }
    
  }
  else{
    
    MPI_Recv (d,ngdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    nullvr (dd,ndof);
    plg->globlocfeti (top,d,dd);
    subv (dd,rhs,ndof);
    
    
    /*
      fprintf (out,"\n\n\n inicializace --- kontrola vektoru prave strany dd pred funkci funkci ldl_feti_sky \n");
      for (i=0;i<ndof;i++){
      fprintf (out," %le",dd[i]);
      }
    */
	
	
    //Sm_sky->ldl_feti_sky (pp,dd,ppd->nrbm,rbmi,ppd->zero);
    gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
    
    
    /*
    fprintf (out,"\n\n\n inicializace --- kontrola vektoru reseni pp funkci ldl_feti_sky \n");
    for (i=0;i<ndof;i++){
      fprintf (out," %le",pp[i]);
    }
    */
    

    nullvr (p,ngdof);
    plg->locglobfeti (top,p,pp);
    MPI_Send (p,ngdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
    

    for (i=0;i<nicg;i++){

      //fprintf (stdout,"\n\n JEDEME GRADIENTY, proc %d, iterace %ld",myrank,i);
      
      MPI_Recv (d,ngdof+1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      if (d[ngdof]>0.5){  break;  }
      
      nullvr (dd,ndof);
      plg->globlocfeti (top,d,dd);
      
      
      /*
      fprintf (out,"\n\n\n krok %ld         kontrola vektoru prave strany pred funkci ldl_feti_sky \n",i);
      for (j=0;j<ndof;j++){
	fprintf (out," %le",dd[j]);
      }
      */
      
      
      //Sm_sky->ldl_feti_sky (pp,dd,ppd->nrbm,rbmi,ppd->zero);
      gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
      
      
      /*
      fprintf (out,"\n\n\n krok %ld          kontrola vektoru reseni pp funkci ldl_feti_sky \n",i);
      for (j=0;j<ndof;j++){
	fprintf (out," %le",pp[j]);
      }
      */
      
      
      nullvr (p,ngdof);
      plg->locglobfeti (top,p,pp);
      
      
      MPI_Send (p,ngdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
      
      //*******************************************************
      //*******************************************************
      //*******************************************************
      //  vlozka s predpodminenim
      //*******************************************************
      //*******************************************************
      //*******************************************************

      MPI_Recv (p,ngdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      nullvr (dd,ndof);
      plg->globlocfeti (top,p,dd);
      nullvr (pp,ndof);
      //lumpedprec (dd,pp);
      for (j=0;j<ndof;j++){
	pp[j]=dd[j];
      }
      nullvr (p,ngdof);
      plg->locglobfeti (top,p,pp);
      MPI_Send (p,ngdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
      
      //*******************************************************
      //*******************************************************
      //*******************************************************
      //  konec vlozky s predpodminenim
      //*******************************************************
      //*******************************************************
      //*******************************************************
      
    }
  }
  
}


/**
   function computes displacements from Lagrange multipliers
   function is used in FETI method
   
   @param w - array containing Lagrange multipliers
   @param d - array containing displacements
   @param f - array containing load forces
   @param rbm - array containing rigid body motions
   @param h - matrix H (columns are rigid body motions)
   @param h1 - matrix (H^T . H)^{-1}
   @param ngdof - number of Lagrange multipliers
   @param hsize - number of columns of matrix H
   
   13.3.2002
*/
void psolver::lagrmultdispl (gtopology *top,paral *plg,gmatrix *gm,
			     double *w,double *d,double *f,double *rbm,long *rbmi,
			     double *h,double *h1)
{
  long i,j,k;
  double *g,*p,*pp,*av,*gamma;
  MPI_Status stat;
  
  g = new double [ngdof];
  p = new double [ngdof];
  pp = new double [ndof];
  av = new double [hsize];
  

  if (myrank==0){
    
    
    /*
    fprintf (out,"\n\n\n kontrola Lagrangeovych multiplikatoru \n");
    for (j=0;j<ngdof;j++){
      fprintf (out," %lf",w[j]);
    }
    */
    
    
    for (j=1;j<nproc;j++){
      MPI_Send (w,ngdof,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
    }
    nullvr (pp,ndof);
    plg->globlocfeti (top,w,pp);

    
    /*
    fprintf (out,"\n\n\n kontrola vektoru pp \n");
    for (j=0;j<ndof;j++){
      fprintf (out," %lf",pp[j]);
    }
    */
    

    subv (pp,f,ndof);
    
    /*
    fprintf (out,"\n\n\n kontrola vektoru pp \n");
    for (j=0;j<ndof;j++){
      fprintf (out," %lf",pp[j]);
    }
    */
    
    /*
    fprintf (out,"\n\n\n kontrola matice tuhosti po eliminaci");
    for (i=0;i<ndof;i++){
      fprintf (out,"\n");
      for (j=Sm_sky->adr[i];j<Sm_sky->adr[i+1];j++){
	fprintf (out," %lf",Sm_sky->a[j]);
      }
    }
    */
    
    //Sm_sky->ldl_feti_sky (d,pp,pdd->nrbm,rbmi,ppd->zero);
    gm->ldl_feti (d,pp,nrbm,rbmi,zero);
    
    /*
    fprintf (out,"\n\n\n kontrola vektoru pp \n");
    for (j=0;j<ndof;j++){
      fprintf (out," %lf",pp[j]);
    }

    fprintf (out,"\n\n\n kontrola vektoru d \n");
    for (j=0;j<ndof;j++){
      fprintf (out," %lf",pp[j]);
    }
    */
    

    nullvr (g,ngdof);
    plg->locglobfeti (top,g,d);
    for (j=1;j<nproc;j++){
      MPI_Recv(p,ngdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      addv (g,p,ngdof);
    }
    
    gamma = new double [hsize];
    
    mtv (h,g,av,ngdof,hsize);
    mv (h1,av,gamma,hsize,hsize);
    scalarray (gamma,-1.0,hsize);
    
    for (j=0;j<nproc;j++){
      k=plg->rbmadr[j];
      for (i=0;i<plg->nrbmdom[j];i++){
	av[i]=gamma[k];  k++;
      }
      if (plg->domproc[j]!=0)
	MPI_Send (av,hsize,MPI_DOUBLE,plg->domproc[j],myrank,MPI_COMM_WORLD);
    }

  }

  else{
    MPI_Recv (p,ngdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    nullvr (pp,ndof);
    plg->globlocfeti (top,p,pp);

    /*
    fprintf (out,"\n\n\n kontrola vektoru pp \n");
    for (j=0;j<ndof;j++){
      fprintf (out," %lf",pp[j]);
    }
    */
    
    subv (pp,f,ndof);
    
    /*
    fprintf (out,"\n\n\n kontrola vektoru pp \n");
    for (j=0;j<ndof;j++){
      fprintf (out," %lf",pp[j]);
    }
    */


    //Sm_sky->ldl_feti_sky (d,pp,ppd->nrbm,rbmi,ppd->zero);
    gm->ldl_feti (d,pp,nrbm,rbmi,zero);
    nullvr (g,ngdof);
    plg->locglobfeti (top,g,d);
    MPI_Send (g,ngdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
    
    
    MPI_Recv (av,hsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    
  }
  
  
  /*
  fprintf (out,"\n\n kontrola posunu bez RBM\n");
  for (i=0;i<ndof;i++){
    fprintf (out," %lf",d[i]);
  }
  */
  
  scalarray (d,-1.0,ndof);
  
  for (i=0;i<ndof;i++){
    for (j=0;j<nrbm;j++){
      d[i]+=rbm[j*ndof+i]*av[j];
    }
  }
  
  
  /*
  fprintf (out,"\n\n kontrola posunu s RBM\n");
  for (i=0;i<ndof;i++){
    fprintf (out," %lf",d[i]);
  }
  */
  
  
  delete [] av;  delete [] pp;  delete [] p;  delete [] g;
}

/**
   function makes preconditioning by lumped matrix
   
   @param v - input vector
   @param pv - preconditioned vector
   
   8.4.2002
*/
void psolver::lumpedprec (gtopology *top,double *v,double *pv)
{
/*
  long i,j,k,ndofe;
  ivector cn;
  vector av,lv;
  matrix sm;
  
  for (i=0;i<top->ne;i++){
    ndofe=top->get_ndofe (i);
    allocm (ndofe,ndofe,sm);
    allocv (ndofe,cn);
    allocv (ndofe,av);
    allocv (ndofe,lv);
    stiffmat (i,0,sm);
    top->give_code_numbers (i,cn.a);
    for (j=0;j<ndofe;j++){
      av[j]=0.0;
      k=cn[j]-1;
      if (k>-1){
	av[j]=v[k];
      }
    }
    mxv (sm,av,lv);
    for (j=0;j<ndofe;j++){
      k=cn[j]-1;
      if (k>-1){
	pv[k]+=lv[j];
      }
    }
    destrm (sm);
    destrv (lv);
    destrv (av);
    destrv (cn);
  }
*/
}




//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************
//********************************************************************

void psolver::par_linear_solver (gtopology *top,paral *plg,gmatrix *gm,
				 double *lhs,double *rhs,FILE *out)
{
  long i,j,maxnrbm;
  long *rbmi;
  double *condmat,*condvect,*rbm,*h,*q,*hh,*ih,*ee,*lm;
  
  switch (tdd){
  case lprimaldd:{
    //  allocation of reduced matrix
    switch (rssol){
    case lpgemp:{
      if (myrank==0){
	condmat = new double [plg->maxnrdof*plg->maxnrdof];
	memset (condmat,0,plg->maxnrdof*plg->maxnrdof*sizeof(double));
      }
      else{
	condmat = new double [(ndof-indof)*(ndof-indof)];
	memset (condmat,0,(ndof-indof)*(ndof-indof)*sizeof(double));
      }
      break;
    }
    case lplu:{
      if (myrank==0){
	condmat = new double [plg->maxnrdof*plg->maxnrdof];
	memset (condmat,0,plg->maxnrdof*plg->maxnrdof*sizeof(double));
      }
      else{
	condmat = new double [(ndof-indof)*(ndof-indof)];
	memset (condmat,0,(ndof-indof)*(ndof-indof)*sizeof(double));
      }
      break;
    }
    case lpldl:{
      if (myrank==0){
	condmat = new double [plg->maxnrdof*plg->maxnrdof];
	memset (condmat,0,plg->maxnrdof*plg->maxnrdof*sizeof(double));
      }
      else{
	condmat = new double [(ndof-indof)*(ndof-indof)];
	memset (condmat,0,(ndof-indof)*(ndof-indof)*sizeof(double));
      }
      break;
    }
    case lpcg:{
      condmat = new double [(ndof-indof)*(ndof-indof)];
      memset (condmat,0,(ndof-indof)*(ndof-indof)*sizeof(double));
      break;
    }
    default:{
      fprintf (stderr,"\n\n (%d) unknown solver of reduced system is required",myrank);
      fprintf (stderr,"\n in function par_sol_red_sys (%s, line %d).\n",__FILE__,__LINE__);
    }
    }
    condvect = new double [plg->maxnrdof];
    memset (condvect,0,plg->maxnrdof*sizeof(double));
    
    //  elimination of inner DOFs = static condensation of inner DOFs
    //Sm_sky->ldlkon_sky (condmat,lhs,rhs,ndof-Indof,1);
    gm->condense (condmat,lhs,rhs,ndof-indof,1);
    
    j=0;
    for (i=indof;i<ndof;i++){
      condvect[j]=rhs[i];  j++;
    }
    
    /*
    fprintf (out,"\n\n\n kontrola zkondenzovanych matic %d \n",myrank);
    for (i=0;i<ndof-Indof;i++){
      fprintf (out,"\n");
      for (j=0;j<ndof-Indof;j++){
	fprintf (out," %lf",condmat[i*(ndof-Indof)+j]);
      }
    }
    fprintf (out,"\n\n kontrola zkondenzovaneho vektoru\n");
    for (i=0;i<ndof-Indof;i++){
      fprintf (out,"  %lf",condvect[i]);
    }
    */
    

    //  reduced system solution
    switch (rssol){
    case lpgemp:{
      redsys_pargemp (plg,condmat,condvect,out);
      break;
    }
    case lplu:{
      redsys_parlu (plg,condmat,condvect,out);
      break;
    }
    case lpldl:{
      redsys_parldl (plg,condmat,condvect,out);
      break;
    }
    case lpcg:{
      redsys_parcg (plg,condmat,condvect,out);
      break;
    }
    default:{
      fprintf (stderr,"\n\n (%d) unknown solver of reduced system is required",myrank);
      fprintf (stderr,"\n in function par_sol_red_sys (%s, line %d).\n",__FILE__,__LINE__);
    }
    }
    
    
    j=0;
    for (i=indof;i<ndof;i++){
      lhs[i]=condvect[j];  j++;
    }
    
    //  back substitution on subdomains
    //Sm_sky->ldlkon_sky (condmat,lhs,rhs,ndof-Indof,2);
    gm->condense (condmat,lhs,rhs,ndof-indof,2);

    delete [] condmat;
    delete [] condvect;

    break;
  }
    //*******************************************************************
    //*******************************************************************
    //*******************************************************************
    //*******************************************************************
  case lfetidd:{
    
    
    //  rigid body motion array allocation
    rbm = new double [enrbm*ndof];
    memset (rbm,0,enrbm*ndof*sizeof(double));
    rbmi = new long [enrbm];
    for (i=0;i<enrbm;i++){
      rbmi[i]=-1;
    }

    fprintf (out,"\n enrbm   %ld",enrbm);
    fprintf (out,"\n ndof    %ld",ndof);
    

    //  computation of kernel of system matrix
    //Sm_sky->ker (rbm,nse,ppd->se,ppd->ense,ppd->lithr,2);
    //Sm_sky->ker (rbm,ppd->nrbm,rbmi,ppd->enrbm,ppd->lithr,3);
    gm->kernel (rbm,nrbm,rbmi,enrbm,lithr,3);

    fprintf (stdout,"\n proc %d  nrbm   %ld",myrank,nrbm);
    fprintf (stdout,"\n proc %d  enrbm  %ld",myrank,enrbm);
    fprintf (stdout,"\n proc %d  lithr  %e",myrank,lithr);

    fprintf (out,"\n proc %d  nrbm   %ld",myrank,nrbm);
    fprintf (out,"\n proc %d  enrbm  %ld",myrank,enrbm);
    fprintf (out,"\n proc %d  lithr  %e",myrank,lithr);
    
    //  investigation of H matrix sizes
    hmatrixsize (plg,rbm,maxnrbm);
    

    
    fprintf (out,"\n\n\n kontrola RBM");
    for (i=0;i<ndof;i++){
      fprintf (out,"\n");
      for (j=0;j<nrbm;j++){
	fprintf (out," %lf",rbm[j*ndof+i]);
      }
    }
    
    
    
    if (myrank==0){
      hsize=plg->rbmadr[nproc];
      h = new double [ngdof*hsize];
      q = new double [hsize];
      lm = new double [ngdof];
      
      fprintf (out,"\n hsize   %ld",hsize);
      
      for (i=0;i<=nproc;i++){
	fprintf (out,"\n rbmadr %4ld   %4ld",i,plg->rbmadr[i]);
      }

    }
    
    
    //  assembling of H matrix
    hmatrix (top,plg,h,rbm,maxnrbm);
    
    
    

    if (myrank==0){
      fprintf (out,"\n\n\n\n kontrola matice H");
      for (i=0;i<ngdof;i++){
	fprintf (out,"\n");
	for (j=0;j<hsize;j++){
	  fprintf (out," %lf",h[i*hsize+j]);
	}
      }
    }




    //  assembling of q vector
    qvector (plg,q,rbm,rhs,maxnrbm);
    
    //  H^T . H and (H^T . H)^{-1.0} evaluation
    if (myrank==0){
      ih = new double [hsize*hsize];
      hh = new double [hsize*hsize];
      ee = new double [hsize*hsize];
      for (i=0;i<hsize;i++){
	for (j=0;j<hsize;j++){
	  ee[i*hsize+j]=0.0;
	}
	ee[i*hsize+i]=1.0;
      }
      
      //  old
      //mtmccr (h,h,hh,ngdof,hsize,hsize);

      mtm (h,h,hh,ngdof,hsize,hsize);


      if (myrank==0){
	fprintf (out,"\n\n\n\n kontrola matice (H^T.H)");
	for (i=0;i<hsize;i++){
	  fprintf (out,"\n");
	  for (j=0;j<hsize;j++){
	    fprintf (out," %lf",hh[i*hsize+j]);
	  }
	}
      }


      gemp (hh,ih,ee,hsize,hsize,zero,1);
      
      
      

      if (myrank==0){
	fprintf (out,"\n\n\n\n kontrola matice (H^T.H)^-1");
	for (i=0;i<hsize;i++){
	  fprintf (out,"\n");
	  for (j=0;j<hsize;j++){
	    fprintf (out," %lf",ih[i*hsize+j]);
	  }
	}
      }


      
      mtm (h,h,hh,ngdof,hsize,hsize);
      mm (hh,ih,ee,hsize,hsize,hsize);
      

      if (myrank==0){
	fprintf (out,"\n\n\n\n kontrola matice jednotkove matice");
	for (i=0;i<hsize;i++){
	  fprintf (out,"\n");
	  for (j=0;j<hsize;j++){
	    fprintf (out," %lf",ee[i*hsize+j]);
	  }
	}
      }

      
      
      delete [] ee;
    }

    //  modified conjugate gradient method
    

    mpcg (top,plg,gm,lm,rhs,q,h,ih,rbmi,out);
    
    fprintf (out,"\n\n hsize %ld",hsize);

    lagrmultdispl (top,plg,gm,lm,lhs,rhs,rbm,rbmi,h,ih);
    
    break;
  }
  default:{
    fprintf (stderr,"\n\n (%d) unknown domain decomposition solver is required",myrank);
    fprintf (stderr,"\n in function linear_solver (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

}

