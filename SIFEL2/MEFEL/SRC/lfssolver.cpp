#include "lfssolver.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "lhsrhs.h"
#include "globmat.h"
#include "loadcase.h"
#include "gmatrix.h"
#include "mechprint.h"
#include "flsubdom.h"
#include "elemswitch.h"


/**
   function solves linear problem with floating subdomains
   it may be used for pull-out tests
   
   JK, 11. 10. 2007, modified 3. 8. 2012
*/
void solve_linear_floating_subdomains ()
{
  long i;
  double *lhs,*rhs;

  //  stiffness matrix assembling
  stiffness_matrix (0);
  
  //  array is allocated in Lsrs
  lhs = Lsrs->give_lhs (0);
  rhs = Lsrs->give_rhs (0);
  
  //  loop over load cases
  for (i=0;i<Lsrs->nlc;i++){
    //  right hand side of system
    mefel_right_hand_side (i,rhs);
  }
  
  //  solution of equation system
  for (i=0;i<Lsrs->nlc;i++){
    Mp->ssle->solve_system (Gtm,Smat,Lsrs->give_lhs(i),Lsrs->give_rhs(i),Out);
  }
  
  print_init(-1, "wt");    
  for (i=0;i<Lsrs->nlc;i++){
    //  computes and prints required quantities
    //compute_ipstrains (i);
    //compute_ipstresses (i);
    compute_req_val (i);
    print_step(i, 0, 0.0, NULL);    
  }
  print_close();

}





/**
   function solves linear problem with floating subdomains
   it may be used for pull-out tests
   
   JK, 11.10.2007
*/
void solve_linear_floating_subdomains2 ()
{


  //Out=fopen("dupnode","w");

  Gtm->edges (Out);
  Gtm->edgenode_sorting ();
  Gtm->prev_next_edges ();
  Gtm->edge_series (Out);
  Gtm->edge_dirvect ();
  Gtm->edge_normvect ();
  
  Gtm->create_ltg (Out);
  
  long i;
  
  /*
  for (i=0;i<Gtm->nged;i++){
    fprintf (Out,"\n %6ld %6ld       %6ld %6ld",Gtm->gedges[i].nlist[0],Gtm->gedges[i].nlist[1],Gtm->gedges[i].nlist[2],Gtm->gedges[i].nlist[3]);
  }
  
  fprintf (Out,"\n");
  for (i=0;i<Gtm->nged;i++){
    fprintf (Out,"\n %6ld       %6ld",Gtm->gedges[i].nlist[0]+1,Gtm->gedges[i].nlist[1]+1);
  }
  i=Gtm->nged-1;
  fprintf (Out,"\n %6ld       %6ld",Gtm->gedges[i].nlist[2]+1,Gtm->gedges[i].nlist[3]+1);
  
  







  //fprintf (Out,"\n\n kontrola prev a next \n");
  //for (i=0;i<Gtm->nged;i++){
  //fprintf (Out,"\n hrana %6ld   prev %6ld  next %6ld",i,Gtm->gedges[i].prev,Gtm->gedges[i].next);
  //}
  //fprintf (Out,"\n\n");

  
  fprintf (Out,"\n\n kontrola souslednosti\n");
  long *av;
  av = new long [Gtm->nged];

  for (i=0;i<Gtm->nged;i++){
    av[i]=-1;
  }
  for (i=0;i<Gtm->nged;i++){
    if (Gtm->gedges[i].prev==-1){
      a=i;
      n=Gtm->gedges[i].next;
      
      fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
      fprintf (Out,"\n %6ld %6ld",Gtm->gedges[a].nlist[0],Gtm->gedges[a].nlist[2]);
      //fprintf (Out,"\n %6ld",Gtm->gedges[a].nlist[0]);
      
      av[i]=1;
      break;
    }
  }
  long ui=Gtm->nged;
  
  for (i=1;i<ui;i++){
    a=n;
    n=Gtm->gedges[n].next;
    av[a]=1;
    fprintf (Out,"\n %6ld",Gtm->gedges[a].nlist[0]);
    //fprintf (Out,"\n %6ld %6ld",Gtm->gedges[a].nlist[0],Gtm->gedges[a].nlist[2]);
    fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
    if (n<0){
      for (j=0;j<Gtm->nged;j++){
	if (Gtm->gedges[j].prev==-1 && av[j]==-1){
	  a=j;
	  n=Gtm->gedges[j].next;
	  av[j]=1;
	  ui--;
	  fprintf (Out,"\n %6ld",Gtm->gedges[a].nlist[0]);
	  /fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
	  break;
	}
      }
    }
  }
  fprintf (Out,"\n %6ld",Gtm->gedges[a].nlist[2]);
  

  fprintf (Out,"\n");


  for (i=0;i<Gtm->nged;i++){
    av[i]=-1;
  }
  for (i=0;i<Gtm->nged;i++){
    if (Gtm->gedges[i].prev==-1){
      a=i;
      n=Gtm->gedges[i].next;

      //fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
      //fprintf (Out,"\n %6ld %6ld ",Gtm->gedges[a].nlist[1]-205,Gtm->gedges[a].nlist[3]-205);
      fprintf (Out,"\n %6ld",Gtm->gedges[a].nlist[1]-205);

      av[i]=1;
      break;
    }
  }
  ui=Gtm->nged;
  
  for (i=1;i<ui;i++){
    a=n;
    n=Gtm->gedges[n].next;
    av[a]=1;
    fprintf (Out,"\n %6ld",Gtm->gedges[a].nlist[1]-205);
    //fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
    if (n<0){
      for (j=0;j<Gtm->nged;j++){
	if (Gtm->gedges[j].prev==-1 && av[j]==-1){
	  a=j;
	  n=Gtm->gedges[j].next;
	  av[j]=1;
	  ui--;
	  fprintf (Out,"\n %6ld",Gtm->gedges[a].nlist[1]-205);
	  //fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
	  break;
	}
      }
    }
  }
  fprintf (Out,"\n %6ld",Gtm->gedges[a].nlist[3]-205);





  fprintf (Out,"\n");


  for (i=0;i<Gtm->nged;i++){
    av[i]=-1;
  }
  for (i=0;i<Gtm->nged;i++){
    if (Gtm->gedges[i].prev==-1){
      a=i;
      n=Gtm->gedges[i].next;

      fprintf (Out,"\n %lf %lf   %lf %lf ",Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1]);
      //fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
      //fprintf (Out,"\n %6ld %6ld ",Gtm->gedges[a].nlist[1]-205,Gtm->gedges[a].nlist[3]-205);

      av[i]=1;
      break;
    }
  }
  ui=Gtm->nged;
  
  for (i=1;i<ui;i++){
    a=n;
    n=Gtm->gedges[n].next;
    av[a]=1;
    //fprintf (Out,"\n %6ld %6ld",Gtm->gedges[a].nlist[1]-205,Gtm->gedges[a].nlist[3]-205);
    //fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
    fprintf (Out,"\n %lf %lf   %lf %lf",Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1]);
    if (n<0){
      for (j=0;j<Gtm->nged;j++){
	if (Gtm->gedges[j].prev==-1 && av[j]==-1){
	  a=j;
	  n=Gtm->gedges[j].next;
	  av[j]=1;
	  ui--;
	  fprintf (Out,"\n %6ld %6ld",Gtm->gedges[a].nlist[1]-205,Gtm->gedges[a].nlist[3]-205);
	  //fprintf (Out,"\n %6ld -> %6ld   %lf %lf %lf      %lf %lf %lf",a,n,Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
	  fprintf (Out,"\n %lf %lf   %lf %lf",Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1]);
	  break;
	}
      }
    }
  }















  //a=0; n=Gtm->gedges[a].next;
  //fprintf (Out,"\n\n kontrola serie \n");
  //for (i=0;i<Gtm->nged;i++){
  //fprintf (Out,"\n%ld -> %ld",a,n);
  //a=n;
  //n=Gtm->gedges[a].next;
  //}
  


  fprintf (Out,"\n\n kontrola souslednosti\n");
  long *av;
  av = new long [Gtm->nged];
  
  for (i=0;i<Gtm->nged;i++){
    av[i]=-1;
  }
  for (i=0;i<Gtm->nged;i++){
    if (Gtm->gedges[i].prev==-1){
      a=i;
      n=Gtm->gedges[i].next;
      
      fprintf (Out,"\n %6ld %6ld     %6ld %6ld",Gtm->gedges[a].nlist[0],Gtm->gedges[a].nlist[2],Gtm->gedges[a].nlist[1],Gtm->gedges[a].nlist[3]);
      fprintf (Out,"         %lf %lf %lf      %lf %lf %lf",Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
      
      av[i]=1;
      break;
    }
  }
  long ui=Gtm->nged;
  
  for (i=1;i<ui;i++){
    a=n;
    n=Gtm->gedges[n].next;
    av[a]=1;
    fprintf (Out,"\n %6ld %6ld     %6ld %6ld",Gtm->gedges[a].nlist[0],Gtm->gedges[a].nlist[2],Gtm->gedges[a].nlist[1],Gtm->gedges[a].nlist[3]);
    fprintf (Out,"         %lf %lf %lf      %lf %lf %lf",Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
    if (n<0){
      for (j=0;j<Gtm->nged;j++){
	if (Gtm->gedges[j].prev==-1 && av[j]==-1){
	  a=j;
	  n=Gtm->gedges[j].next;
	  av[j]=1;
	  ui--;
	  fprintf (Out,"\n %6ld %6ld     %6ld %6ld",Gtm->gedges[a].nlist[0],Gtm->gedges[a].nlist[2],Gtm->gedges[a].nlist[1],Gtm->gedges[a].nlist[3]);
	  fprintf (Out,"         %lf %lf %lf      %lf %lf %lf",Gtm->gedges[a].dv[0],Gtm->gedges[a].dv[1],Gtm->gedges[a].dv[2],Gtm->gedges[a].nv[0],Gtm->gedges[a].nv[1],Gtm->gedges[a].nv[2]);
	  break;
	}
      }
    }
  }
  

  fprintf (Out,"\n\n\n");






  for (i=0;i<Gtm->nged;i++){
    fprintf (Out,"\n\n hrana cislo  %ld",i);
    Gtm->gedges[i].print (Out);
  }

  
  fprintf (Out,"\n\n\n kontrola serii \n");
  for (i=0;i<Gtm->nged;i++){
    fprintf (Out,"\n edge %6ld je v serii  %6ld",i,Gtm->edgeser[i]);
  }

  //fclose (Out);

  */



  double *lhs,*rhs;
  
  //  stiffness matrix assembling
  stiffness_matrix (0);
  
  //  array is allocated in Lsrs
  lhs = Lsrs->give_lhs (0);
  rhs = Lsrs->give_rhs (0);
  
  //  loop over load cases
  for (i=0;i<Lsrs->nlc;i++){
    //  right hand side of system
    mefel_right_hand_side (i,rhs);
  }
  
  
  
  /*
  Smat->printmat(Out);
  FILE *p;
  p = fopen ("vektor.txt","w");
  fprintf (p,"%ld\n",Ndofm);
  for (i=0;i<Ndofm;i++){
    fprintf (p,"%20.16le\n",rhs[i]);
  }
  fclose (p);
  */
  
  
  
  
  //Smat->printdiag(Out);
  
  //  solution of equation system
  for (i=0;i<Lsrs->nlc;i++){
    Mp->ssle->solve_system (Gtm,Smat,Lsrs->give_lhs(i),Lsrs->give_rhs(i),Out);
    //Smat->solve_system (Lsrs->give_lhs(i),Lsrs->give_rhs(i));
  }
 
  print_init(-1, "wt");    
  for (i=0;i<Lsrs->nlc;i++){
    //  computes and prints required quantities
    //compute_ipstrains (i);
    //compute_ipstresses (i);
    compute_req_val (i);
    print_step(i, 0, 0.0, NULL);    
  }
  print_close();

 




}

/**
   function solves linear problem with floating subdomains
   it may be used for pull-out tests
   
   JK, 5.11.2005
*/
void solve_linear_floating_subdomains_old ()
{
  /*
  long i,lcid;
  double *lhs,*rhs;
  //flsubdom *fsd;
  
  lcid=0;
  
  Fsd = new flsubdom;
  
  lhs=Lsrs->give_lhs (lcid);
  rhs=Lsrs->give_rhs (lcid);
  
  mefel_right_hand_side (lcid,rhs);
  
  stiffness_matrix (lcid);
  
  Fsd->initialization (Ndofm,Mp->ense,rhs);
  
  Fsd->solve_lin_alg_system (lhs,rhs);
  
  print_init(-1, "wt");    
  for (i=0;i<Lsrs->nlc;i++){
    compute_ipstrains (i);
    compute_ipstresses (i);
    compute_req_val (i);
    print_step(i, 0, 0.0, NULL);    
  }
  print_close();
  */





  long i;
  double *lhs,*rhs;
  
  //  stiffness matrix assembling
  stiffness_matrix (0);
  
  //  array is allocated in Lsrs
  lhs = Lsrs->give_lhs (0);
  rhs = Lsrs->give_rhs (0);
  
  //  loop over load cases
  for (i=0;i<Lsrs->nlc;i++){
    //  right hand side of system
    mefel_right_hand_side (i,rhs);
  }
  


  
  
  /*
  Smat->printmat(Out);
  FILE *p;
  p = fopen ("vektor.txt","w");
  fprintf (p,"%ld\n",Ndofm);
  for (i=0;i<Ndofm;i++){
    fprintf (p,"%20.16le\n",rhs[i]);
  }
  fclose (p);
  */
  
  
  
  
  //Smat->printdiag(Out);
  
  double *condmat,*condvect;
  long nrdof=Ndofm-Gtm->nbdof;
  condmat = new double [nrdof*nrdof];
  condvect = new double [nrdof];
  
  //  solution of equation system
  for (i=0;i<Lsrs->nlc;i++){
    //Mp->ssle.solve_system (Gtm,Smat,Lsrs->give_lhs(i),Lsrs->give_rhs(i),Out);
    //Smat->solve_system (Lsrs->give_lhs(i),Lsrs->give_rhs(i));
    Smat->condense (Gtm,condmat,condvect,Lsrs->give_lhs(i),Lsrs->give_rhs(i),nrdof,1,Out);
  }
  
  long j;
  FILE *out;
  out = fopen ("matice.txt","w");
  for (i=0;i<nrdof;i++){
    for (j=0;j<nrdof;j++){
      fprintf (out,"%le ",condmat[i*nrdof+j]);
    }
    fprintf (out,"\n");
  }
  fclose (out);
  out = fopen ("vektor.txt","w");
  for (i=Gtm->nbdof;i<Ndofm;i++){
    fprintf (out,"%le\n",rhs[i]);
  }
  fclose (out);
  
  
  out = fopen ("nezname.txt","w");
  for (i=0;i<Mt->nn;i++){
    fprintf (out,"%6ld %6ld\n",Gtm->gnodes[i].cn[0],Gtm->gnodes[i].cn[1]);
  }
  fclose (out);
  


  print_init(-1, "wt");    
  for (i=0;i<Lsrs->nlc;i++){
    //  computes and prints required quantities
    //compute_ipstrains (i);
    //compute_ipstresses (i);
    compute_req_val (i);
    print_step(i, 0, 0.0, NULL);    
  }
  print_close();

  delete [] condmat;
  delete [] condvect;
}

