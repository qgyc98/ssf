#include "mpi.h"
#include "plssolver.h"
#include "pglobal.h"
#include "seqfile.h"
#include "genfile.h"
#include <string.h>

void par_solve_nonlinear_statics ()
{
  long i,j,k,n,ni,ini,stop,modif,li;
  double a0,a1,a2,d,l1,l2,dl,dlmax,dlmin,psi,lambda,blambda,dlambda,ddlambda,norf,norfp,norv,zero,ierr;
  double lambdao,ss1,ss2,ss3,ss4,ss5;
  double *r,*ra,*ddr,*u,*v,*f,*fa,*fp,*fc,*fi;
  
  FILE *gr;
  gr = fopen (Mp->auxfile,"wt");
  
  //  number of rows of the matrix
  n = Ndof;
  //  maximum number of increments
  ni = Mp->nial;
  //  maximum number of iterations in one increment
  ini = Mp->niilal;
  //  computer zero
  zero = Mp->zero;
  //  required error in inner loop
  ierr = Mp->erral;
  //  length of arc
  dl = Mp->dlal;
  //  maximum length of arc
  dlmax = Mp->dlmaxal;
  //  minimum length of arc
  dlmin = Mp->dlminal;
  //  displacement-loading driving switch
  psi = Mp->psial;
  
  //  allocation of auxiliary arrays
  r = new double [n];
  ddr = new double [n];
  u = new double [n];
  v = new double [n];
  f = new double [n];
  fa = new double [n];
  fi = new double [n];
  
  //  allocation of backup arrays for state variables in integration points
  Mm->alloc_backup ();
  
  //  initialization phase
  ra = Lsrs->give_lhs (lcid);
  fp = Lsrs->give_rhs (lcid*2);
  fc = Lsrs->give_rhs (lcid*2+1);
  
  //  array containing selected DOFs
  seldofinit ();
  
  if (Mp->hdbackupal==2){
    arclopen (li,n,lambda,dl,ra,fp);    
  }
  else{
    lambda=0.0;
    lambdao=0.0;
    li=0;
  }
  
  //  norm of proportionality vector
  //norfp = ss (fp,fp,n);
  //  norfp = par_norfp ();
  
  modif=0;
  
  //***************************
  //  main iteration loop  ****
  //***************************
  for (i=li;i<ni;i++){
    
    //  backup of left hand side vector
    for (j=0;j<n;j++){
      r[j]=ra[j];
    }
    //  backup of reached lambda parameter
    blambda=lambda;
    //  backup of state variables
    Mm->backup_state_var ();
    
    
    aux_nonlin_print (gr,ra,lambda);
    //fprintf (gr,"%e %e\n",ra[0],lambda);
    
    
    
    fprintf (stdout,"\n arc-length: increment %ld   lambda %e  dl %e",i,lambda,dl);
    if (Mespr==1){
      
      //fprintf (stderr,"\n\n arc-length: increment %ld  lambda  %e    dl %e",i,lambda,dl);
      //fprintf (Out,"\n\n *******************************************************************");
      //fprintf (Out,"\n arc-length: increment %ld,   lambda=%e   dl=%e",i,lambda,dl);
      
      //print_displacements (Out,0);
      //print_strains (Out,0);
      //print_stresses (Out,0);
      //print_other (Out,0);
      
      if (Mp->adp==1)
	Mp->ado->print (Out,i,lambda,0);
      
    }
    
    //  assembling of tangent stiffness matrix
    stiffness_matrix (lcid);
    
    //  backup of the fp, in ldl_sky the right hand side will be destroyed
    for (j=0;j<n;j++){
      f[j]=fp[j];
    }
    //  solution of K(r).v=F
    //Smat->solve_system (v,f);
    Psol->par_linear_solver (Gt,Smat,v,f,Out,Mespr);
    
    
    //  generalized norm of displacement increments
    //norv = par_displincr (v,n);
    
    
    //if (fabs(lambdao)>fabs(lambda)){
    //dlambda = -1.0*dl*dl/sqrt(norv+psi*psi*norfp);
    //
    //else{
    //dlambda = dl*dl/sqrt(norv+psi*psi*norfp);
    //}
    //lambdao=lambda;
    
    
    dlambda = dl/sqrt(norv+psi*psi*norfp);
    
    //fprintf (stderr,"\n dlambda  %e",dlambda);
    
    
    for (j=0;j<n;j++){
      ddr[j]=dlambda*v[j];
      ra[j]+=ddr[j];
      fa[j]=fc[j]+(lambda+dlambda)*fp[j];
    }
    
    
    
    ddlambda=dlambda;
    
    //  computation of internal forces
    internal_forces (lcid,fi);
    
    for (k=0;k<n;k++){
      f[k]=fa[k]-fi[k];
    }
    
    norf=ss(f,f,n);
    norf=sqrt(norf);
    
    if (Mespr==1)  fprintf (stdout,"\n %e %e norf %e",lambda,dl,norf);
    
    if (norf<ierr){
      //******************************************
      //  no inner iteration loop is required  ***
      //******************************************
      modif++;

      if (modif>1){
	//  arc length modification
	dl*=2.0;
	if (dl>dlmax)  dl=dlmax;
	modif=0;
	
	if (Mespr==1){
	  fprintf (stdout,"\n arc length must be modified (dl increase) dl=%e    norf=%e",dl,norf);
	  fprintf (Out,"\n arc length must be modified (dl increase) dl=%e    norf=%e",dl,norf);
	}
      }
      
      lambda+=dlambda;
      continue;
    }
    
    //****************************
    //  inner iteration loop  ****
    //****************************
    stop=0;
    for (j=0;j<ini;j++){
      
      //  back substitution
      //Smat->solve_system (u,f);
      Psol->par_linear_solver (Gt,Smat,u,f,Out,Mespr);
      
      
      for (k=0;k<n;k++){
	f[k]=ddr[k];
	ddr[k]+=u[k];
      }
      
      //  coefficient of quadratic equation
      a2=norv+psi*psi*norfp;
      a1=2.0*ss(ddr,v,n)+2.0*ddlambda*psi*psi*norfp;
      a0=ss(ddr,ddr,n)+ddlambda*ddlambda*psi*psi*norfp-dl*dl;
      
      /*
      fprintf (Out,"\n\n\n Kontrola kvadraticke rovnice");
      fprintf (Out,"\n norfp     %15.10le",norfp);
      fprintf (Out,"\n norv      %15.10le",norv);
      fprintf (Out,"\n (ddr,v)   %15.10le",ss(ddr,v,n));
      fprintf (Out,"\n (ddr,ddr) %15.10le",ss(ddr,ddr,n));
      fprintf (Out,"\n dl        %15.10le",dl);
      fprintf (Out,"\n ddlambda  %15.10le",ddlambda);
      fprintf (Out,"\n psi       %15.10le",psi);
      fprintf (Out,"\n a2        %15.10le",a2);
      fprintf (Out,"\n a1        %15.10le",a1);
      fprintf (Out,"\n a0        %15.10le",a0);
      fprintf (Out,"\n discrim   %15.10le",a1*a1-4.0*a2*a0);
      fprintf (Out,"\n\n");
      */
      
      //  solution of quadratic equation
      if (fabs(a2)<zero){
	if (fabs(a1)<zero){
	  if (fabs(a0)>zero)  fprintf (stderr,"\n\n nonsolvable constrained condition in function arclength");
	  else  fprintf (stderr,"\n\n infinite number of solution of constrained condition in function arclength");
	}
	else{
	  dlambda=(0.0-a0)/a1;
	}
      }
      else{
	d=a1*a1-4.0*a2*a0;
	
	if (d<0.0){
	  fprintf (stderr,"\n negative discriminant in function arclength");
	  break;
	}
	
	
	l1=(0.0-a1+sqrt(d))/2.0/a2;
	l2=(0.0-a1-sqrt(d))/2.0/a2;
	
	//if (fabs(l1)<fabs(l2))  dlambda=l1;
	//else  dlambda=l2;
      }
      
      //  zacatek novinky
      
      ss1=0.0;  ss2=0.0;
      ss3=0.0;  ss4=0.0;  ss5=0.0;
      for (k=0;k<n;k++){
	ss1+=(ddr[k]+l1*v[k])*f[k];
	ss2+=(ddr[k]+l2*v[k])*f[k];
	ss3+=(ddr[k]+l1*v[k])*(ddr[k]+l1*v[k]);
	ss4+=(ddr[k]+l2*v[k])*(ddr[k]+l2*v[k]);
	ss5+=f[k]*f[k];
      }
      
      if (ss1/ss3/ss5>ss2/ss4/ss5)  dlambda=l1;
      else          dlambda=l2;
      //  konec novinky
      
      for (k=0;k<n;k++){
	ddr[k]+=dlambda*v[k];
	ra[k]+=u[k]+dlambda*v[k];
	fa[k]+=dlambda*fp[k];
      }
      ddlambda+=dlambda;
      
      //  computation of internal forces
      internal_forces (lcid,fi);
      
      for (k=0;k<n;k++){
	f[k]=fa[k]-fi[k];
      }
      
      norf=ss(f,f,n);
      norf=sqrt(norf);
      
      if (Mespr==1){
	fprintf (stdout,"\n increment  %ld     inner loop  %ld    norf=%e",i,j,norf);
	fprintf (Out,"\n increment  %ld     inner loop  %ld    norf=%e",i,j,norf);
      }
      
      if (norf<ierr){
	lambda+=ddlambda;
	stop=1;  break;
      }
    }
    
    modif=0;
    
    if (stop==0){
      //  modification of the arc length
      dl/=2.0;
      if (dl<dlmin){  dl=dlmin;  break; }
      
      if (Mespr==1){
	fprintf (stdout,"\n arc length must be modified (dl decrease) dl=%e    norf=%e",dl,norf);
	fprintf (Out,"\n arc length must be modified (dl decrease) dl=%e    norf=%e",dl,norf);
      }
      
      //  restoring of left hand side vector
      
      for (j=0;j<n;j++){
	ra[j]=r[j];
      }
      
      //  restoring of lambda parameter
      lambda=blambda;
      //  restoring of state variables
      Mm->restore_state_var ();
      
      
    }
  }
  
  
  if (Mp->hdbackupal==1){
    arclsave (i,n,lambda,dl,ra,fp);
  }
  
  
  delete [] r;		    
  delete [] fi;  
  delete [] fa;  
  delete [] f;
  delete [] v;  
  delete [] u;  
  delete [] ddr;
  
  fclose (gr); 
  
}

/**
   function assembles code numbers of selected
   
   JK, 23.6.2004
*/
void parseldofinit ()
{
  long i,j,k,l,ndofn;
  double x1,x2,y1,y2,z1,z2,length;
  
  switch (Mp->nlman->displnorm){
  case alldofs:{  break; }
  case seldofs:
  case seldofscoord:{
    for (i=0;i<Mp->nlman->nsdofal;i++){
      j=Mp->nlman->seldofal[i];
      Mp->nlman->seldofal[i]=Mt->give_dof (Mp->nlman->selnodal[i],j)-1;
    }
    break;
  }
  case selnodes:{
    Mp->nlman->nsdofal=0;
    for (i=0;i<Mp->nlman->nsnal;i++){
      Mp->nlman->nsdofal+=Mt->give_ndofn (Mp->nlman->selnodal[i]);
    }
    Mp->nlman->seldofal = new long [Mp->nlman->nsdofal];
    k=0;
    for (i=0;i<Mp->nlman->nsnal;i++){
      ndofn=Mt->give_ndofn (Mp->nlman->selnodal[i]);
      for (j=0;j<ndofn;j++){
	l=Mt->give_dof (Mp->nlman->selnodal[i],j);
	if (l>0){
	  Mp->nlman->seldofal[k]=l-1;
	  k++;
	}
      }
    }
    Mp->nlman->nsdofal=k;
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown displacement norm is required in function seldofinit (%s, line %d)",__FILE__,__LINE__);
  }
  }

}
