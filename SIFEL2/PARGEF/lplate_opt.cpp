#include "lplate.h"

/**
   constructor
   
   @param np - number of processors
   @param mr - rank of processor
   @param nd - number of domain
   
   JK
*/
lplate::lplate (int np,int mr,int nd,gtopology *top)
{
  long i,j,k;
  
  zero=1.0e-12;

  nproc=np;  myrank=mr;  ndom=nd;
  
  nl = nproc;
  nn = top->nn;
  ndofn = top->gnodes[0].ndofn;
  nnmult = 4;
  
  if (myrank==0){
    lgnodes = new lgnode [nn];
    for (i=0;i<nn;i++){
      lgnodes[i].cn = new long* [nl+1];
      for (j=0;j<=nl;j++){
	lgnodes[i].cn[j] = new long [nnmult];
	for (k=0;k<nnmult;k++){
	  lgnodes[i].cn[j][k]=0;
	}
      }
    }
    
    cn = new long** [nl];
    for (i=0;i<nl;i++){
      cn[i] = new long* [nn];
      for (j=0;j<nn;j++){
	cn[i][j] = new long [ndofn];
	for (k=0;k<ndofn;k++){
	  cn[i][j][k]=0;
	}
      }
    }
    
    thick = new double* [nl];
    for (i=0;i<nl;i++){
      thick[i] = new double [nn];
    }
    
    constrmat = new double* [(nl-1)*nn*nnmult];
    for (i=0;i<(nl-1)*nn*nnmult;i++){
      constrmat[i] = new double [nl*2];
      for (j=0;j<nl*2;j++){
	constrmat[i][j]=0.0;
      }
    }
    
  }
  
}

lplate::~lplate ()
{
  long i,j;
  
  if (myrank==0){
    
    for (i=0;i<nl;i++){
      delete [] thick[i];
    }
    delete [] thick;
    
    for (i=0;i<nl;i++){
      for (j=0;j<nn;j++){
	delete [] cn[i][j];
      }
      delete [] cn[i];
    }
    delete [] cn;
    
    for (i=0;i<(nl-1)*nn*nnmult;i++){
      delete [] constrmat[i];
    }
    delete [] constrmat;
  }
  
}


/**
   function generates global code numbers

   function assembles local code numbers of subdomains (code
   numbers of nodal displacements)
   then generates global code numbers (code
   number of multipliers)
   
   @param top - general topology
   @param domproc - array with domain-processor correspondation
   @param out - outpu stream
   
   JK
*/
void lplate::globcnnum_lpp (gtopology *top,long *domproc,FILE *out)
{
  long i,j,k,buffsize;
  long *buff;
  MPI_Status stat;

  //  code numbers assembling
  buffsize = nn*ndofn;
  buff = new long [buffsize];
  
  k=0;
  for (i=0;i<nn;i++){
    for (j=0;j<ndofn;j++){
      buff[k]=top->gnodes[i].cn[j];
      k++;
    }
  }
  
  if (myrank==0){
    assemble_cn (buff,ndom);
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      assemble_cn (buff,domproc[stat.MPI_TAG]);
    }
    
    //  searching for maximum number of DOFs on subdomain
    maxndofdom ();
    
    for (i=1;i<nproc;i++){
      MPI_Send (&maxndof,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
    /*
    //  beginning of auxiliary print
    fprintf (out,"\n\n\n kontrola sebranych kodovych cisel \n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n\n DOMAIN number %ld",i);
      for (j=0;j<nn;j++){
	fprintf (out,"\n node %4ld   ",j);
	for (k=0;k<ndofn;k++){
	  fprintf (out,"  %ld",cn[i][j][k]);
	}
      }
    }
    //  end of auxiliary print
    */

    //  generation of code numbers of multipliers
    codenum_mult (nmult);
    
    /*
    //  beginning of auxiliary print
    fprintf (out,"\n\n\n kontrola nagenerovanych kodovych cisel \n");
    for (i=0;i<nn;i++){
      fprintf (out,"\n layered node number %4ld  ",i);
      for (j=0;j<=nl;j++){
	fprintf (out,"   |%ld|  ",j);
	for (k=0;k<nnmult;k++){
	  fprintf (out,"  %ld",lgnodes[i].cn[j][k]);
	}
      }
    }
    //  end of auxiliary print
    */


    globcnn_auxiliary ();
    
    fprintf (out,"\n\n kontrola adr \n");
    for (i=0;i<nproc+1;i++){
      fprintf (out,"\n %4ld    %ld",i,adr[i]);
    }
    
    //  beginning of auxiliary print
    fprintf (out,"\n\n\n kontrola globalnich kodovych cisel \n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n\n DOMAIN number %ld",i);
      for (j=0;j<nn;j++){
	fprintf (out,"\n node %4ld   ",j);
	for (k=0;k<ndofn;k++){
	  fprintf (out,"  %ld",gcn[i][j][k]);
	}
      }
    }
    //  end of auxiliary print

  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (&maxndof,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;

}

void lplate::globcnn_auxiliary ()
{
  long i,j,k;
  
  gcn = new long** [nl];
  for (i=0;i<nl;i++){
    gcn[i] = new long* [nn];
    for (j=0;j<nn;j++){
      gcn[i][j] = new long [ndofn];
      for (k=0;k<ndofn;k++){
	gcn[i][j][k]=0;
      }
    }
  }
  
  adr = new long [nproc+1];
  adr[0]=0;
  
  ngdof=1;
  for (i=0;i<nl;i++){
    for (j=0;j<nn;j++){
      for (k=0;k<ndofn;k++){
	if (cn[i][j][k]>0){
	  gcn[i][j][k]=ngdof;
	  ngdof++;
	}
      }
    }
    adr[i+1]=ngdof-1;
  }
  ngdof--;
  
}


/**
   function collects thicknesses at nodes from subdomains
   they are stored at the master processor
   function assembles the constraint matrix on the master processor
   
   @param th - array containing thicknesses in nodes
   @param domproc - array containing domain-processor correspondation
   @param out - output stream

   JK, 23.2.2003
*/
void lplate::constraintmat (double *th,long *domproc,FILE *out)
{
  long i;
  double *rbuff;
  MPI_Status stat;
  
  //   thicknesses assembling
  rbuff = new double [nn];
  for (i=0;i<nn;i++){
    rbuff[i]=th[i];
  }
  
  if (myrank==0){
    assemble_thicknesses (rbuff,ndom);
    
    for (i=1;i<nproc;i++){
      MPI_Recv (rbuff,nn,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      assemble_thicknesses (rbuff,domproc[stat.MPI_TAG]);
    }
  }
  else{
    MPI_Send (rbuff,nn,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] rbuff;
  
  
  
  if (myrank==0){
    assemble_constr ();

    //  for debugging only
    //assemble_constr_dm ();
  }
}

/**
   function investigates maximum number of degrees of freedom on domain
   function can be used only on master processor
   
   JK, 12.8.2003
*/
void lplate::maxndofdom ()
{
  long i,j,k;
  
  maxndof=0;
  for (i=0;i<nproc;i++){
    for (j=0;j<nn;j++){
      for (k=0;k<ndofn;k++){
	if (cn[i][j][k]>maxndof)  maxndof=cn[i][j][k];
      }
    }
  }
}

/**
   function assembles array of nodal code numbers on master processor
   can be used only on master processor
   
   @param buff - array with one layer data
   @param ndom - number of layer
   
   JK, 23.2.2003
*/
void lplate::assemble_cn (long *buff,long ndom)
{
  long i,j,k;
  
  k=0;
  for (i=0;i<nn;i++){
    for (j=0;j<ndofn;j++){
      cn[ndom][i][j]=buff[k];
      k++;
    }
  }
}

/**
   function assembles array of thicknesses on master processor
   can be used only on master processor
   
   @param buff - array with one layer data
   @param ndom - number of layer
   
   JK, 23.2.2003
*/
void lplate::assemble_thicknesses (double *buff,long ndom)
{
  long i;
  
  for (i=0;i<nn;i++){
    thick[ndom][i]=buff[i];
  }
}

/**
   function defines code numbers for multipliers
   can be used only on master processor
   
   @param nm - number of Lagrange multipliers
   
   JK, 23.2.2003
*/
void lplate::codenum_mult (long &nm)
{
  long i,j;
  
  nm=1;
  for (i=1;i<nl;i++){
    for (j=0;j<nn;j++){
      if (cn[i-1][j][0]>0 || cn[i-1][j][4]>0 || cn[i][j][0]>0 || cn[i][j][4]>0){
	lgnodes[j].cn[i][0]=nm;  nm++;
      }
      if (cn[i-1][j][1]>0 || cn[i-1][j][3]>0 || cn[i][j][1]>0 || cn[i][j][3]>0){
	lgnodes[j].cn[i][1]=nm;  nm++;
      }
      if (cn[i-1][j][2]>0 || cn[i][j][2]>0){
	lgnodes[j].cn[i][2]=nm;  nm++;
      }
      if (cn[i-1][j][5]>0 || cn[i][j][5]>0){
	lgnodes[j].cn[i][3]=nm;  nm++;
      }
    }
  }
  nm--;
}

/**
   function assembles constraint matrix on the master processor
   matrix has special pattern which lead to the special storage
   there are at most two nonzero entries in one row before orthonormalization
   and two.(number of subdomains) nonzero entries after orthonormalization
   
   Lagrange multipliers are defined gradually on interfaces (one interface after another)
   it is different system than in sequential description where multipliers are
   generated graduelly on nodes (one node after another, across interfaces)
   
   this function can be used only on master processor
   
   JK, 24.2.2003
*/
void lplate::assemble_constr ()
{
  long i,j,ii;
  double lthickness,uthickness;
  
  ii=0;
  for (i=1;i<nl;i++){
    for (j=0;j<nn;j++){
      lthickness = thick[i-1][j];
      uthickness = thick[i][j];
      
      if (cn[i-1][j][0]>0)  constrmat[ii][(i-1)*2+0]=-1.0;
      if (cn[i-1][j][4]>0)  constrmat[ii][(i-1)*2+1]=lthickness/(-2.0);
      if (cn[i][j][0]>0)    constrmat[ii][i*2+0]=1.0;
      if (cn[i][j][4]>0)    constrmat[ii][i*2+1]=uthickness/(-2.0);
      ii++;
      
      if (cn[i-1][j][1]>0)  constrmat[ii][(i-1)*2+0]=-1.0;
      if (cn[i-1][j][3]>0)  constrmat[ii][(i-1)*2+1]=lthickness/2.0;
      if (cn[i][j][1]>0)    constrmat[ii][i*2+0]=1.0;
      if (cn[i][j][3]>0)    constrmat[ii][i*2+1]=uthickness/2.0;
      ii++;
      
      if (cn[i-1][j][2]>0)  constrmat[ii][(i-1)*2+0]=-1.0;
      if (cn[i][j][2]>0)    constrmat[ii][i*2+0]=1.0;
      ii++;
      
      if (cn[i-1][j][5]>0)  constrmat[ii][(i-1)*2+0]=-1.0;
      if (cn[i][j][5]>0)    constrmat[ii][i*2+0]=1.0;
      ii++;
    }
  }

}

/**
   function assembles constarint matrix stored in dense matrix storage
   
   JK, 12.8.2003
*/
void lplate::assemble_constr_dm ()
{
  long i,j,k,l,*cnd,*cnm;
  double tl,tu;
  matrix lcm;
  
  cnd = new long [2*ndofn];
  cnm = new long [nnmult];
  allocm (nnmult,2*ndofn,lcm);
  
  fprintf (stdout,"\n pocet vsech multiplikatoru %ld",nmult);
  fprintf (stdout,"\n pocet vsech posunu  %ld",ndof);
  fprintf (stdout,"\n pocet vsech posunu  %ld",ndof*nproc);
  

  allocm (nmult,ngdof,dcm);

  for (i=1;i<nl;i++){
    for (j=0;j<nn;j++){

      tl = thick[i-1][j];
      tu = thick[i][j];
      
      for (k=0;k<ndofn;k++){
	cnd[k]=gcn[i-1][j][k];
      }
      for (k=0;k<ndofn;k++){
	cnd[k+ndofn]=gcn[i][j][k];
      }
      for (k=0;k<nnmult;k++){
	cnm[k]=lgnodes[j].cn[i][k];
      }
      
      /*
      fprintf (stdout,"\n\n");
      for (k=0;k<2*ndofn;k++){
	fprintf (stdout,"  %ld",cnd[k]);
      }
      fprintf (stdout,"X");
      for (k=0;k<nnmult;k++){
	fprintf (stdout,"  %ld",cnm[k]);
      }
      */
      

      fillm (0.0,lcm);
      
      lcm[0][0] = -1.0;
      lcm[1][1] = -1.0;
      lcm[2][2] = -1.0;
      lcm[1][3] = tl/2.0;
      lcm[0][4] = tl/(-2.0);
      lcm[3][5] = -1.0;
      
      lcm[0][6]  = 1.0;
      lcm[1][7]  = 1.0;
      lcm[2][8]  = 1.0;
      lcm[1][9]  = tu/2.0;
      lcm[0][10] = tu/(-2.0);
      lcm[3][11] = 1.0;
      

      for (k=0;k<nnmult;k++){
	if (cnm[k]==0)  continue;
	for (l=0;l<2*ndofn;l++){
	  if (cnd[l]==0)  continue;
	  else  dcm[cnm[k]-1][cnd[l]-1]=lcm[k][l];
	}
      }
      
    }
  }
  
  delete [] cnd;
  delete [] cnm;
  destrm (lcm);

}


/**
   optimalizovat!

   function orthonormalizes the constraint matrix
   constraint matrix is stored in a special way because it is extremely sparse
   
   function can be used only on master processor
   
   JK, 26.2.2003
*/
void lplate::orthonormalization ()
{
  long i,j,k,ii,jj;
  double alpha,nor;
  
  //  normalization of the first block
  for (i=0;i<nn*nnmult;i++){
    nor = constrmat[i][0]*constrmat[i][0] + constrmat[i][1]*constrmat[i][1];
    nor += constrmat[i][2]*constrmat[i][2] + constrmat[i][3]*constrmat[i][3];
    if (nor>zero){
      nor=sqrt(nor);
      constrmat[i][0]/=nor;
      constrmat[i][1]/=nor;
      constrmat[i][2]/=nor;
      constrmat[i][3]/=nor;
    }
  }
  
  ii=nn*nnmult;
  for (i=1;i<nl-1;i++){
    jj=ii-nn*nnmult;
    for (j=0;j<nn*nnmult;j++){
      alpha=constrmat[ii][i*2+0]*constrmat[jj][(i+1)*2+0] + constrmat[ii][i*2+1]*constrmat[jj][(i+1)*2+1];
      nor=0.0;
      for (k=0;k<2*(i+2);k++){
	constrmat[ii][k]-=alpha*constrmat[jj][k];
	nor+=constrmat[ii][k]*constrmat[ii][k];
      }
      if (nor>zero){
	nor=sqrt(nor);
	for (k=0;k<2*(i+2);k++){
	  constrmat[ii][k]/=nor;
	}
      }
      ii++;  jj++;
    }
  }
}

/**
   function orthonormalizes rows of matrix stored in dense matrix format
   
   JK, 12.8.2003
*/
void lplate::orthonormalization_dm ()
{
  long i,j,k;
  double nor,alpha;
  
  nor=0.0;
  for (i=0;i<ngdof;i++){
    nor+=dcm[0][i]*dcm[0][i];
  }
  nor=sqrt(nor);
  if (nor>zero){
    for (i=0;i<ngdof;i++){
      dcm[0][i]/=nor;
    }
  }
  
  for (i=0;i<nmult;i++){
    for (j=i+1;j<nmult;j++){
      alpha=0.0;
      for (k=0;k<ngdof;k++){
	alpha+=dcm[j][k]*dcm[i][k];
      }
      nor=0.0;
      for (k=0;k<ngdof;k++){
	dcm[j][k]-=alpha*dcm[i][k];
	nor+=dcm[j][k]*dcm[j][k];
      }
      nor=sqrt(nor);
      if (nor>zero){
	for (k=0;k<ngdof;k++){
	  dcm[j][k]/=nor;
	}
      }
    }
  }

}

/**
   function localizes components of local %vector into global %vector
   can be used only on master processor
   
   G lv = gv
   
   @param gv - global %vector
   @param lv - local %vector
   @param ns - number of subdomain
   
   JK, 4.3.2003
*/
void lplate::locglob (double *gv,double *lv,long ns)
{
  long i,j,li,ii,cnl,cn1,cn2;
  double cm1,cm2,lv1,lv2;
  
  li=ns-1;
  if (li<0)  li=0;
  ii=li*nn*nnmult;
  for (i=li;i<nl-1;i++){
    for (j=0;j<nn;j++){

      cnl=lgnodes[j].cn[i+1][0];
      if (cnl==0){  ii++; }
      else{
	cm1=constrmat[ii][ns*2+0];
	cm2=constrmat[ii][ns*2+1];
	cn1=cn[ns][j][0];
	if (cn1>0)  lv1=lv[cn1-1];
	else        lv1=0.0;
	cn2=cn[ns][j][4];
	if (cn2>0)  lv2=lv[cn2-1];
	else        lv2=0.0;
	gv[cnl-1]+=cm1*lv1+cm2*lv2;
	ii++;
      }

      cnl=lgnodes[j].cn[i+1][1];
      if (cnl==0){  ii++; }
      else{
	cm1=constrmat[ii][ns*2+0];
	cm2=constrmat[ii][ns*2+1];
	cn1=cn[ns][j][1];
	if (cn1>0)  lv1=lv[cn1-1];
	else        lv1=0.0;
	cn2=cn[ns][j][3];
	if (cn2>0)  lv2=lv[cn2-1];
	else        lv2=0.0;
	gv[cnl-1]+=cm1*lv1+cm2*lv2;
	ii++;
      }


      cnl=lgnodes[j].cn[i+1][2];
      if (cnl==0){  ii++; }
      else{
	cm1=constrmat[ii][ns*2+0];
	cn1=cn[ns][j][2];
	if (cn1>0)  lv1=lv[cn1-1];
	else        lv1=0.0;
	gv[cnl-1]+=cm1*lv1;
	ii++;
      }

      cnl=lgnodes[j].cn[i+1][3];
      if (cnl==0){  ii++; }
      else{
	cm1=constrmat[ii][ns*2+0];
	cn1=cn[ns][j][5];
	if (cn1>0)  lv1=lv[cn1-1];
	else        lv1=0.0;
	gv[cnl-1]+=cm1*lv1;
	ii++;
      }
    }
  }
}


void lplate::locglob_auxiliary (double *gv,double *lv,long did,FILE *out)
{
  long i,j;
  double s,*av;
  
  av = new double [ngdof];
  
  /*
  nullvr (blb,ngdof);
  for (i=0;i<nn;i++){
    for (j=0;j<ndofn;j++){
      if (gcn[did][i][j]>0){
	blb[gcn[did][i][j]-1]=lv[cn[did][i][j]-1];
      }
    }
  }
  */

  nullvr (av,ngdof);
  j=0;
  for (i=adr[did];i<adr[did+1];i++){
    av[i]=lv[j];  j++;
  }
  
  /*
  fprintf (out,"\n\n kontrola locglob, domena %ld",did);
  for (i=0;i<ngdof;i++){
    fprintf (out,"\n %4ld      %e     %e",i,blb[i],av[i]);
  }
  */
  
  for (i=0;i<nmult;i++){
    s=0.0;
    for (j=0;j<ngdof;j++){
      s+=dcm[i][j]*av[j];
    }
    gv[i]+=s;
  }
  
  delete [] av;
}
void lplate::globloc_auxiliary (double *gv,double *lv,long did)
{
  long i,j;
  double s,*av;
  
  av = new double [ngdof];

  nullvr (av,ngdof);
  for (i=0;i<ngdof;i++){
    s=0.0;
    for (j=0;j<nmult;j++){
      s+=dcm[j][i]*gv[j];
    }
    av[i]=s;
  }
  
  /*
  for (i=0;i<nn;i++){
    for (j=0;j<ndofn;j++){
      if (gcn[did][i][j]>0){
	lv[cn[did][i][j]-1]=av[gcn[did][i][j]-1];
      }
    }
  }
  */

  j=0;
  for (i=adr[did];i<adr[did+1];i++){
    lv[j]=av[i];  j++;
  }
  
  delete [] av;
}


/**
   function localizes components of global %vector into local %vector
   can be used only on master processor
   
   G^T gv = lv
   
   @param gv - global %vector
   @param lv - local %vector
   @param ns - number of subdomain
   
   JK, 4.3.2003
*/
void lplate::globloc (double *gv,double *lv,long ns)
{
  long i,j,li,ii,cnl,cn1,cn2;
  double cm1,cm2,cgv;
  
  li=ns-1;
  if (li<0)  li=0;
  ii=li*nn*nnmult;
  for (i=li;i<nl-1;i++){
    for (j=0;j<nn;j++){

      cnl=lgnodes[j].cn[i+1][0];
      if (cnl==0){  ii++; }
      else{
	cm1=constrmat[ii][ns*2+0];
	cm2=constrmat[ii][ns*2+1];
	cgv=gv[cnl-1];
	
	cn1=cn[ns][j][0];
	if (cn1>0)  lv[cn1-1]+=cgv*cm1;
	cn2=cn[ns][j][4];
	if (cn2>0)  lv[cn2-1]+=cgv*cm2;
	ii++;
      }

      cnl=lgnodes[j].cn[i+1][1];
      if (cnl==0){  ii++; }
      else{
	cm1=constrmat[ii][ns*2+0];
	cm2=constrmat[ii][ns*2+1];
	cgv=gv[cnl-1];
	
	cn1=cn[ns][j][1];
	if (cn1>0)  lv[cn1-1]+=cgv*cm1;
	cn2=cn[ns][j][3];
	if (cn2>0)  lv[cn2-1]+=cgv*cm2;
	ii++;
      }

      cnl=lgnodes[j].cn[i+1][2];
      if (cnl==0){  ii++; }
      else{
	cm1=constrmat[ii][ns*2+0];
	cgv=gv[cnl-1];
	
	cn1=cn[ns][j][2];
	if (cn1>0)  lv[cn1-1]+=cgv*cm1;
	ii++;
      }
      
      cnl=lgnodes[j].cn[i+1][3];
      if (cnl==0){  ii++; }
      else{
	cm1=constrmat[ii][ns*2+0];
	cgv=gv[cnl-1];

	cn1=cn[ns][j][5];
	if (cn1>0)  lv[cn1-1]+=cgv*cm1;
	ii++;
      }
    }
  }

}

/**
   function prints constraint matrix stored in efficient storage scheme
   
   @param out - output stream
   
   JK
*/
void lplate::printconstrmat (FILE *out)
{
  long i,j;
  
  fprintf (out,"\n\n\n Constraint matrix \n");
  for (i=0;i<(nl-1)*nn*4;i++){
    fprintf (out,"\n row %4ld ",i);
    for (j=0;j<nl*2;j++){
      fprintf (out,"  %6.3e",constrmat[i][j]);
      
    }
  }
}

/**
   function prints constraint matrix stored in dense matrix format
   
   @param out - output stream
   
   JK
*/
void lplate::printconstrmat_dm (FILE *out)
{
  long i,j;
  
  fprintf (out,"\n\n\n Constraint matrix \n");
  for (i=0;i<nmult;i++){
    fprintf (out,"\n row %4ld ",i);
    for (j=0;j<nproc*ndof;j++){
      fprintf (out,"  %6.3e",dcm[i][j]);
      
    }
  }
}


/**
   function solves system of equations by conjugate gradient method
   
   @param gm - general %matrix, %matrix of one subdomain
   @param w - %vector of unknown multipliers
   @param rhs - %vector of right hand side, %vector of nodal load
   @param domproc - array containing domain-processor correspondation
   @param out - output stream
   
   JK, 12.8.2003
*/
void lplate::cg (gmatrix *gm,double *w,double *rhs,long *domproc,FILE *out)
{
  long i,j;
  double nom,denom,alpha,beta,norrhs,pr1;
  double *d,*dd,*r,*p,*pp;
  MPI_Status stat;
  
  fprintf (stdout,"\n ndof     %ld",ndof);
  fprintf (stdout,"\n maxndof  %ld",maxndof);
  fprintf (stdout,"\n nmult    %ld",nmult);

  fprintf (out,"\n ndof     %ld",ndof);
  fprintf (out,"\n maxndof  %ld",maxndof);
  fprintf (out,"\n nmult    %ld",nmult);


  //  direction vector allocation
  dd = new double [maxndof+1];
  //  auxiliary vector allocation
  pp = new double [maxndof+1];
  
  nom=0.0;
  for (i=0;i<ndof;i++){
    nom+=rhs[i]*rhs[i];
  }
  
  fprintf (out,"\n\n kontrola normy prave strany  %e",nom);
  
  if (myrank==0){
    r = new double [nmult];
    d = new double [nmult];
    p = new double [nmult];
    
    nullvr (w,nmult);
    
    //  G^T w
    for (j=1;j<nproc;j++){
      nullvr (dd,maxndof+1);
      globloc (w,dd,domproc[j]);
      //globloc_auxiliary (w,dd,domproc[j]);
      MPI_Send (dd,maxndof,MPI_DOUBLE,j,ndom,MPI_COMM_WORLD);
    }
    nullvr (dd,maxndof+1);
    globloc (w,dd,domproc[0]);
    //globloc_auxiliary (w,dd,domproc[0]);
    

    fprintf (out,"\n\n kontrola vektoru dd \n");
    for (i=0;i<maxndof;i++){
      fprintf (out,"\n   %e",dd[i]);
    }
    
    //  f-G^T w
    for (j=0;j<ndof;j++){
      dd[j]=rhs[j]-dd[j];
    }
    

    //  K^+ (f - G^T w)
    gm->back_substitution (pp,dd);
    
    // G K^+ (f - G^T w)
    nullvr (r,nmult);
    locglob (r,pp,domproc[0]);
    //locglob_auxiliary (r,pp,domproc[0],out);
    //  norm of right hand side
    norrhs=nom;
    
    for (j=1;j<nproc;j++){
      MPI_Recv(pp,maxndof+1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      locglob (r,pp,domproc[stat.MPI_TAG]);
      
      /*
      fprintf (out,"\n\n kontrola zasilky, domena %ld \n",domproc[stat.MPI_TAG]);
      for (i=0;i<maxndof;i++){
	fprintf (out,"\n %4ld     %e",i,pp[i]);
      }
      */

      //locglob_auxiliary (r,pp,domproc[stat.MPI_TAG],out);
      norrhs+=pp[maxndof];
    }
    
    fprintf (out,"\n norm of right hand side vector  %e,",norrhs);
    fprintf (stdout,"\n norm of right hand side vector  %e,",norrhs);
    

    fprintf (out,"\n\n kontrola vektor reziduii\n");
    for (i=0;i<nmult;i++){
      fprintf (out,"\n %4ld   %e",i,r[i]);
    }


    //  initialization of direction vector
    for (i=0;i<nmult;i++){
      d[i]=r[i];
    }
    
    //  nominator evaluation
    nom = ss (r,r,nmult);
    

    // *************************
    //  main iteration loop
    // *************************
    for (i=0;i<nicg;i++){
      
      // ***********************************************
      //  auxiliary vector computation G K^+ G^T d = p
      // ***********************************************
      
      //  G^T d
      for (j=1;j<nproc;j++){
	nullvr (dd,maxndof+1);
	globloc (d,dd,domproc[j]);
	//globloc_auxiliary (d,dd,domproc[j]);
	MPI_Send (dd,maxndof+1,MPI_DOUBLE,j,ndom,MPI_COMM_WORLD);
      }
      
      nullvr (dd,maxndof+1);
      globloc (d,dd,domproc[0]);
      //globloc_auxiliary (d,dd,domproc[0]);
      
      //  K^+ G^T d
      gm->back_substitution (pp,dd);
      
      
      //  G K^+ G^T d
      nullvr (p,nmult);
      locglob (p,pp,domproc[0]);
      //locglob_auxiliary (p,pp,domproc[0],out);
      
      for (j=1;j<nproc;j++){
	MPI_Recv(pp,maxndof+1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	locglob (p,pp,domproc[stat.MPI_TAG]);
	//locglob_auxiliary (p,pp,domproc[stat.MPI_TAG],out);
      }
      
      //  denominator of alpha
      denom = ss (d,p,nmult);
      //  for print only
      pr1=denom;

      if (fabs(denom)<zero){
	fprintf (stderr,"\n\n zero denominator in modified conjugate gradient method (file %s, line %d).\n",__FILE__,__LINE__);
	
	if (i!=nicg-1){
	  dd[maxndof]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send (dd,maxndof+1,MPI_DOUBLE,j,0,MPI_COMM_WORLD);
	  }
	}
	break;
      }
      
      // *******************
      //  vypocet alpha
      // *******************
      alpha = nom/denom;
      
      // **************************************************************
      //  vypocet noveho gradientu g a nove aproximace neznamych x
      // **************************************************************
      for (j=0;j<nmult;j++){
	w[j]+=alpha*d[j];
	r[j]-=alpha*p[j];
      }
      
      denom = nom;
      if (fabs(denom)<zero){
	fprintf (stderr,"\n\n zero denominator of beta in function (file %s, line %d).\n",__FILE__,__LINE__);
	
	if (i!=nicg-1){
	  dd[maxndof]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send (dd,maxndof+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	  }
	}
	break;
      }
      
      nom = ss (r,r,nmult);
      
      fprintf (stdout,"\n iteration  %4ld        norm r   %e   mod. error   %e    alpha denom  %e",i,sqrt(nom),sqrt(nom)/sqrt(norrhs),pr1);
      fprintf (out,"\n iteration  %4ld        norm r   %e   mod. error   %e    alpha denom  %e",i,sqrt(nom),sqrt(nom)/sqrt(norrhs),pr1);

      
      if (sqrt(nom)/sqrt(norrhs)<errcg){
	if (i!=nicg-1){
	  dd[maxndof]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send (dd,maxndof+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	  }
	}
	break;
      }
      
      
      //  computation of beta coefficient
      beta = nom/denom;
      
      //  new direction vector
      for (j=0;j<nmult;j++){
	d[j]=beta*d[j]+r[j];
      }
      
    }

    
    anicg=i;  aerrcg=nom;
    
    fprintf (out,"\n\n\n\n kontrola Lagrangeovych multiplikatoru \n");
    for (i=0;i<nmult;i++){
      fprintf (out,"\n lagr. mult %4ld    %e",i,w[i]);
    }

  }
  else{
    //  G^T w
    MPI_Recv (dd,maxndof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

    fprintf (out,"\n\n kontrola vektoru dd \n");
    for (i=0;i<maxndof;i++){
      fprintf (out,"\n   %e",dd[i]);
    }
    
    //  f-G^T w
    for (j=0;j<ndof;j++){
      dd[j]=rhs[j]-dd[j];
    }
    

    // K^+ (f - G^T w)
    gm->back_substitution (pp,dd);
    
    //  norm of right hand side vector
    pp[maxndof]=nom;
    
    MPI_Send (pp,maxndof+1,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);
    

    for (i=0;i<nicg;i++){
      
      //  G^T d
      MPI_Recv (dd,maxndof+1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      if (dd[maxndof]>0.5){  break;  }
      
      //  K^+ G^T d
      gm->back_substitution (pp,dd);
      
      // G K^+ G^T d
      MPI_Send (pp,maxndof+1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
      
    }

  }
  
  if (myrank==0){
    delete [] r;
    delete [] d;
    delete [] p;
  }
  



  delete [] dd;
  delete [] pp;

}

/**
   function computes nodal displacements from known Lagrange multipliers
   
   @param gm - general %matrix, %matrix of one subdomain
   @param lhs - %vector of nodal displacements
   @param rhs - %vector of nodal forces
   @param w - %vector of Lagrange multipliers
   @param domproc - array containing domain-processor correspondation
   
   JK, 12.8.2003
*/
void lplate::nodaldisplacements (gmatrix *gm,double *lhs,double *rhs,double *w,long *domproc)
{
  long j;
  double *dd;
  MPI_Status stat;
  
  //  direction vector allocation
  dd = new double [maxndof+1];

  if (myrank==0){
    //  G^T w
    for (j=1;j<nproc;j++){
      nullvr (dd,maxndof+1);
      globloc (w,dd,domproc[j]);
      //globloc_auxiliary (w,dd,domproc[j]);
      MPI_Send (dd,maxndof,MPI_DOUBLE,j,ndom,MPI_COMM_WORLD);
    }
    
    nullvr (dd,maxndof+1);
    globloc (w,dd,domproc[0]);
    //globloc_auxiliary (w,dd,domproc[0]);

    //  f-G^T w
    for (j=0;j<ndof;j++){
      dd[j]=rhs[j]-dd[j];
    }

    //  K^+ (f - G^T w)
    gm->back_substitution (lhs,dd);

  }
  else{
    //  G^T w
    MPI_Recv (dd,maxndof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    
    //  f-G^T w
    for (j=0;j<ndof;j++){
      dd[j]=rhs[j]-dd[j];
    }

    //  K^+ (f - G^T w)
    gm->back_substitution (lhs,dd);

  }
  
  delete [] dd;
}

/**
   function solves system of equations
   
   @param gm - matrix of subdomain
   @param lhs - %vector of unknowns
   @param rhs - %vector of right hand side
   
   JK
*/
void lplate::solve_system (gtopology *top,gmatrix *gm,
			   long *domproc,double *lhs,double *rhs,FILE *out)
{
  double *w;
  
  /*
  fprintf (out,"\n\n kontrola prave strany\n");
  long i;
  for (i=0;i<ndof;i++){
    fprintf (out,"  %e",rhs[i]);
  }
  */
  
  //  array for Lagrange multipliers
  if (myrank==0)
    w = new double [nmult];
  
  //  decomposition of matrices of subdomains
  gm->decompose_matrix ();
  
  //if (myrank==0)  printconstrmat (out);
  //if (myrank==0)  printconstrmat_dm (out);
  
  //  orthonormalization of constraint matrix
  //if (myrank==0)  orthonormalization ();
  //if (myrank==0)  orthonormalization_dm ();
  
  //if (myrank==0)  printconstrmat (out);
  //if (myrank==0)  printconstrmat_dm (out);
  
  //  solution of reduced system by conjugate gradient method
  cg (gm,w,rhs,domproc,out);
  
  //  computation of nodal displacements
  nodaldisplacements (gm,lhs,rhs,w,domproc);
  
  if (myrank==0)
    delete [] w;
  
}
