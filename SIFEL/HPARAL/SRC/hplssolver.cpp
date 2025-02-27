#include <mpi.h>
#include "hplssolver.h"
#include "global.h"
#include "globmat.h"
#include "hpglobal.h"
#include "seqfilesm.h"

#include <string.h>

/**
   function computes linear static problem on a macro-level
   material parameters are not defined in material models, they are
   obtained from a meso or micro level by a homogenization method
   
   JK, 4. 6. 2013
*/
void parallel_homogenization_linear_statics_tiles ()
{
  long i,j,combmin,combmax,ncomb,tid;
  double *lhs,*rhs,**condm,**condvlhs,**condvrhs;
//  double *fl,*fi;
  precond prec;  
//  char name[1100];
  char nameout[1100],namegid[1100],namedat[1100];
  
  strcpy (nameout,Outdm->outfn);
  strcpy (namegid,Outdm->outgrfn);
  strcpy (namedat,Outdm->outdiagfn);

  // *******************************************************************************************
  // *******************************************************************************************
  // **  all processors assemble, factorize and store matrices of each tile type
  // **  no communication is needed and it has to save a time
  // **  in future, a variant where matrices will be distributed via MPI could be implemented
  // *******************************************************************************************
  // *******************************************************************************************

  condm = new double* [Ntiletypes];
  condvlhs = new double* [Ntiletypes];
  condvrhs = new double* [Ntiletypes];

  // cycle over all types of tiles
  for (i=0;i<Ntiletypes;i++){
    //  array for the condensed matrices
    condm[i] = new double [Nbdof[i]*Nbdof[i]];
    memset (condm[i],0,(Nbdof[i]*Nbdof[i])*sizeof(double));
    
    //  array for the condensed vectors
    condvlhs[i] = new double [Nbdof[i]];
    memset (condvlhs[i],0,Nbdof[i]*sizeof(double));
    
    //  array for the condensed vectors
    condvrhs[i] = new double [Nbdof[i]];
    memset (condvrhs[i],0,Nbdof[i]*sizeof(double));
    
    //  condensation - elimination of all internal DOFs
    SSmat[i].condense (Gtm,condm[i],condvrhs[i],Slhs[i],Srhs[i],Nbdof[i],1,Out);
  }


  Ndofm=0;
  for (i=0;i<Ntiles;i++){
    for (j=0;j<Nbdoftiling[i];j++){
      if (Ndofm<Cnbn[i][j])
        Ndofm=Cnbn[i][j];
    }
  }
  
  gtopology *top;
  top = new gtopology;
  top->initiate (Cnbn,Nbdoftiling,Ntiles);


  if (Smat!=NULL)
    delete Smat;
  
  if (Mp!=NULL)
    delete Mp;
  Mp = new probdesc;
  

  Mp->ssle->tlinsol=ldl;
  Mp->tstorsm=skyline_matrix;
  Mp->tprob=linear_statics;
  
  Mp->straincomp = 1;
  Mp->strainpos  = 2;
  Mp->strainaver = 1;
  
  Mp->stresscomp = 1;
  Mp->stresspos  = 2;
  Mp->stressaver = 1;
  
  Mp->otherstate=0;
  
  Smat = new gmatrix;
  Smat->setval (Mp->ssle);
  Smat->initiate (top,Ndofm,Mp->tstorsm,1,Out);
  
  fprintf (stderr,"\n\n\n\n\n pocet DOFs %ld",Ndofm);
  fprintf (stderr,"\n ukladani %d",Smat->ts);
  fprintf (stderr,"\n solver Smat %d",Smat->tlinsol);
  fprintf (stderr,"\n solver Mp   %d",Mp->ssle->tlinsol);

  lhs = new double [Ndofm];
  rhs = new double [Ndofm];
  
  ncomb=Ncomb/Nproc;
  combmin=Myrank*ncomb;
  combmax=(Myrank+1)*ncomb;

  fprintf (stderr,"\n myrank %d   ncomb %ld  combmin %ld   combmax %ld",Myrank,ncomb,combmin,combmax);

  for (i=combmin;i<combmax;i++){
    Smat->sky->nullsky ();
    nullv (rhs,Ndofm);
    
    for (j=0;j<Ntiles;j++){
      //  tile id
      tid = Tiling[i][j];
      //  the last parameter is useless
      Smat->localized (condm[tid],Cnbn[j],Nbdoftiling[tid],0);
      //  localization of the right hand side
      locglob (rhs,condvrhs[tid],Cnbn[j],Nbdoftiling[tid]);
    }
    
    Smat->solve_system (top,prec,lhs,rhs,Out);
    //fprintf (stderr,"\n rhs %le     Smat %le  Ndofm %ld",rhs[Ndofm-1],Smat->sky->a[0],Ndofm);
    for (j=0;j<Ndofm;j++){
      fprintf (stdout,"glob lhs %ld   %le\n",j,lhs[j]);
    }
    
    for (long ii=0;ii<3;ii++){
      
      //  tile id
      tid = Tiling[i][ii];
      fprintf (stderr,"\n myrank %d   i %ld    tid %ld",Myrank,i,tid);
      globloc (lhs,condvlhs[tid],Cnbn[ii],Nbdoftiling[tid]);
      //  back substitution
      SSmat[tid].condense (top,condm[tid],condvlhs[tid],Slhs[tid],Srhs[tid],Nbdoftiling[tid],2,Out);
      
      for (j=0;j<Nidof[tid]+Nbdoftiling[tid];j++){
        Lsrs->lhs[j]=Slhs[tid][j];
        fprintf (stdout,"%ld lhs %ld   %le\n",ii,j,Slhs[tid][j]);
      }
    }
    
    /*
    //  names of output files
    strcpy (name,nameout);
    sprintf (&name[strlen (nameout)],"%d.out",i);
    strcpy (Outdm->outfn, name);
    if (Outdm->gf != grfmt_no)
      {
        strcpy (name,namegid);
        sprintf (&name[strlen (namegid)],"%d",Myrank+1);
        strcpy (Outdm->outgrfn, name);
      }
    if (Outdm->ndiag > 0)
      {
        strcpy (name,namedat);
        sprintf (&name[strlen (namedat)],"%d.dat",Myrank+1);
        strcpy (Outdm->outdiagfn, name);
      }  
    
    print_init(-1, "wt");
    
    Mm = Mmm;
    Mt = Mtt;
    //  computes required quantities
    compute_req_val (0);
    

    if (Outdm->nog.selnforce.st != sel_no)  // nodal forces are required
      {
        // rhs was rewritten by the solver call -> rebuild load vector
        mefel_right_hand_side (0,Lsrs->give_rhs(0),fl);
      }
    
    // in artificial time 0.0 the nodal forces contains only the load vector fl
    print_step(0, 0, 0.0, fl);
    
    // if the nodal forces are required for graphical output
    // reaction + load vector are printed in the artificial time 1.0
    if (Outdm->nog.selnforce.st != sel_no){
      //  indicator of strain computation
      //  it means, strains at integration points have not been computed
      Mp->strainstate=0;
      //  indicator of stress computation
      //  it means, stresses at integration points have not been computed
      Mp->stressstate=0;
      //  indicator of computation of other array
      //  it means, stresses at integration points have not been computed
      Mp->otherstate=0;
      // calculate nodal forces and reactions
      internal_forces(0, fi);
      
      // in artificial time 1.0 the nodal forces contains the load vector + reactions
      print_step(0, 1, 1.0, fi);
    }


    print_close();
    */
  }

}



/**
   function computes linear static problem on a macro-level
   material parameters are obtained by homogenization
   on microscale, stationary problems are solved
   these microproblems can be defined on every finite element or
   on aggregates of finite elements, the aggregates can send
   data from all aggregated elements or averaged data which
   has the same form as data from a single element
   

   material parameters are not defined in material models, they are
   obtained from a meso or micro level by a homogenization method
   
   TKr, 14/07/2022
*/
void parallel_homogenization_linear_statics ()
{

  long i,j,k,ii,jj,kk,n,aux,nelidom,buffsize1,buffsize2,ncomp,dim;
  long *neldom,*buff;
  double end_time,*lhs,*rhs,*fl,*fi;
  double *buff1,*buff2;
  MPI_Status stat;
  
  //  allocation of buffer
  //  buff[0] - the number of connections between the macroproblem and microproblem
  //  buff[1] - the maximum number of connections between the macroproblem and microproblem
  buff = new long [2];
  
  if (Myrank==0){  
    switch (Mp->mami){
    case elem_domain:{
      //  each element is connected with its own microproblem
      buff[0]=1;
      buff[1]=1;
      
      break;
    }
    case aggregate_domain:{
      //  generation of auxiliary data needed in parallel computation
      //  elements of macroproblem are aggregated, each finite element
      //  sends its data to microproblem and obtains its matrices
      
      //  map between finite elements of macroproblem and microproblems
      //  this map is read in the subroutine transtop::read_seq_top
      //  eldom[i]=j - the i-th finite element of macroproblem is connected to the j-th microproblem
      
      //  the number of elements on subdomains
      //  it is the number of macroelements connected to one microproblem
      //  the number of subdomains--microproblems must be Nproc-1
      //  neldom[i]=j - j finite elements of the macroproblem are connected to the i-th microproblem
      //  neldom[0]=0 - the master processor (myrank=0) contains no microproblem
      neldom = new long [Nproc];
      for (i=0;i<Nproc;i++){
	neldom[i]=0;
      }
      for (i=0;i<Mt->ne;i++){
	neldom[Gtm->stop->eldom[i]]++;
      }
      buffsize1=0;
      for (i=0;i<Nproc;i++){
	if (buffsize1<neldom[i])
	  buffsize1=neldom[i];
      }
      
      buff[1]=buffsize1;
      
      break;
    }
    case material_aggregate_domain:{
      //one domain represents one material = one processor
      //the aggregate covers elements of one materyal type
      neldom = new long [Nproc];
      for (i=0;i<Nproc;i++){
	neldom[i]=0;
      }
      
      for (i=0;i<Mt->ne;i++)
	neldom[Gtm->stop->eldom[i]] = 1;
      
      buff[0]=1;
      buff[1]=1;
      
      break;
    }
    case eff_aggregate_domain:{

      neldom = new long [Nproc];
      for (i=0;i<Nproc;i++){
	neldom[i]=0;
      }

      for (i=0;i<Mt->ne;i++){
	if(Gtm->stop->eldom[i] != Mt->elements[i].idm[0])
	  neldom[Gtt->stop->eldom[i]]++;
      }
      
      //  the aggregate sends averaged data to microproblem
      //  only aggregate sends and receives data, the same
      //  matrices are used on all elements in the aggregate
      buff[0]=1;
      buff[1]=1;
      
      break;
    }
    default:{
      par_print_err(Myrank,"unknown type of aggregation is required",__FILE__,__LINE__,__func__);
    }
    }    
  
    for (i=1;i<Nproc;i++){
      if (elem_domain == aggregate_domain)
	buff[0]=neldom[i];

      MPI_Send (buff,2,MPI_LONG,i,Myrank,MPI_COMM_WORLD);
      //MPI_Send_mine (buff,2,MPI_LONG,i,Myrank,MPI_COMM_WORLD);
    }
    if (elem_domain == aggregate_domain)
      buff[0]=neldom[0];
  }
  else{
    MPI_Recv (buff,2,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    //MPI_Recv_mine (buff,2,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  

  //  the number of connections between the macroproblem and microproblem
  //  the number of entities (elements or subdomains) of macroproblem connected with microproblem
  nelidom=buff[0];
  //  the maximum number of connections between the macroproblem and microproblem
  //  the maximum number of entities (elements or subdomains) of macroproblem connected with microproblem
  buffsize1=buff[1];
  
  //  dimension of the problem, one, two or three dimensional case
  //  the dimension is defined with respect to the first finite element
  dim = Mt->give_dimension (0);
  //  the number of stress/strain components
  ncomp = Mt->give_ncomp (0,0);
  
  buffsize2=buffsize1;
  //  buffer for displacments and strain components
  buffsize1*=(dim+ncomp+1);
  //  buffer for the stiffness matrix
  buffsize2*=(ncomp*ncomp);
  
  fprintf (stdout,"\n procesor %d  nelidom   %ld",Myrank,nelidom);
  fprintf (stdout,"\n procesor %d  buffsize1 %ld",Myrank,buffsize1);
  fprintf (stdout,"\n procesor %d  buffsize2 %ld",Myrank,buffsize2);
  
  
  delete [] buff;
  buff1 = new double [buffsize1];
  nullv (buff1,buffsize1);
  buff2 = new double [buffsize2];
  nullv (buff2,buffsize2);
  
  
  /**************************************************************************/

  if (Myrank==0){
    //  the number of unknonws
    //  it is the number of unknowns in macroproblem in the case of the master processor
    //  otherwise, it is the number of unknowns in microproblems
    //  number of mechanical degrees of freedom
    n=Ndofm;
    //  vector of prescribed force loads, it does not contain forces caused by temperature, etc.
    fl   = new double [n];
    nullv (fl,n);
    //  vector of resulting internal forces
    fi   = new double [n];
    nullv (fi,n);
    
    rhs = Lsrs->give_rhs (0);
  }

  // fictious time
  end_time = 0.0;
  
  if (Myrank==0){
    for (ii=1;ii<Nproc;ii++){
      MPI_Send (&end_time,1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
      //MPI_Send_mine (&end_time,1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (&end_time,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    //MPI_Recv_mine (&end_time,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  fprintf (stdout,"\n procesor %d",Myrank);
  
  if (Myrank==0){
    
    // **********************************
    //  solution of microscale problems
    // **********************************
    for (ii=1;ii<Nproc;ii++){
      
      //  function assembles data from the macroproblem for the ii-th microproblem to the array buff1
      // tady pole buff1 musi mit rozmer podle nelidom??!! nebo podle poctu elementu??!!
      higher_to_lower_levelm (ii,buffsize1,buff1);
      
      MPI_Send (buff1,buffsize1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
      //MPI_Send_mine (buff1,buffsize1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (buff1,buffsize1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    //MPI_Recv_mine(buff1,buffsize1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    
    for (ii=0;ii<nelidom;ii++){
      
      fprintf (Out,"\n\n *********************************************");
      fprintf (Out,"\n Myrank  = %ld",Myrank);
      fprintf (Out,"\n nelidom = %ld",ii);
      fprintf (Out,"\n\n pred vypoctem paral_mechanical_homogenization:");
      for (long ijk=0;ijk<buffsize1;ijk++){
	//for (long ijk=0;ijk<1;ijk++){
	  fprintf (Out,"\n buff1 %ld   %le",ijk,buff1[ijk]);
      }
     
      fprintf (Out,"\n");
      for (long ijk=0;ijk<buffsize2;ijk++){
	//for (long ijk=0;ijk<1;ijk++){
	fprintf (Out,"\n buff2 %ld   %le",ijk,buff2[ijk]);
      }
      fflush(Out);
      
      paral_mechanical_homogenization (buff1+ii*(dim+ncomp+1),buff2+ii*(ncomp*ncomp));
      
      fprintf (Out,"\n\n po vypoctu paral_mechanical_homogenization:");
      for (long ijk=0;ijk<buffsize1;ijk++){
	//for (long ijk=0;ijk<1;ijk++){
	fprintf (Out,"\n buff1 %ld   %le",ijk,buff1[ijk]);
      }
      fprintf (Out,"\n");
      for (long ijk=0;ijk<buffsize2;ijk++){
	//for (long ijk=0;ijk<1;ijk++){
	fprintf (Out,"\n buff2 %ld   %le",ijk,buff2[ijk]);
      }
      fflush(Out);
    }
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  /**************************************************************************/
  
  if (Myrank==0){
    for (ii=1;ii<Nproc;ii++){
      MPI_Recv (buff2,buffsize2,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      //MPI_Recv_mine(buff2,buffsize2,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      k=stat.MPI_TAG;
      jj=0;
      
      for (j=0;j<Mt->ne;j++){
	
	//fprintf (Out,"\n Myrank  = %ld",Myrank);
	//fprintf (Out,"\n element no. = %ld",j);
	//fprintf (Out,"\n Gtm->stop->eldom[j] = %ld",Gtm->stop->eldom[j]);
	//fprintf (Out,"\n ii = %ld",ii);
	//fprintf (Out,"\n k  =stat.MPI_TAG = %ld",k);
	
	if (Gtm->stop->eldom[j]==k){
	  //  only elements connected with the k-th microscale are assumed
	  
	  switch (Mp->mami){
	  case elem_domain:{
	    //  each element is connected with its own microproblem
	    Mm->hommatm[j].assemble_matrices (buff2,ncomp,dim);
	    break;
	  }
	    
	  case aggregate_domain:{
	    //  generation of auxiliary data needed in parallel computation
	    Mm->hommatm[j].assemble_matrices (buff2+jj*(ncomp*ncomp),ncomp,dim);
	    jj++;
	    break;
	  }
	  case material_aggregate_domain:{
	    //  aggregate sends averaged data to microproblem
	    //  index jj is not incremented in this case because all elements in aggregate??!!
	    //  obtain the same conductivity and capacity matrices
	    //
	    //  index of subdomain-aggregate = index of material (processor)
	    
	    //fprintf (Out,"\n Myrank  = %ld",Myrank);
	    //fprintf (Out,"\n from Myrank  = %ld",ii);
	    //fprintf (Out,"\n element  = %ld",j);
	    //fprintf (Out,"\n\n pred ulozenim matice na materialu c. %ld:",j);
	    //for (long ijk=0;ijk<buffsize2;ijk++){
	    //fprintf (Out,"\n buff2 %ld   %le",ijk,buff2[ijk]);
	    //}
	    
	    Mm->hommatm[j].assemble_matrices (buff2,ncomp,dim);
	    
	    //fflush(Out);
	    break;
	  }
	  case eff_aggregate_domain:{
	    //  aggregate sends averaged data to microproblem
	    //  index jj is not incremented in this case because all elements in aggregate
	    //  obtain the same conductivity and capacity matrices
	    //
	    //  index of subdomain-aggregate = index of element material
	    kk = Mt->elements[i].idm[0];
	    Mm->hommatm[kk].assemble_matrices (buff2+kk*(ncomp*ncomp),ncomp,dim);
	    break;
	  }
	    
	  default:{
	    par_print_err(Myrank,"unknown type of aggregation is required",__FILE__,__LINE__,__func__);
	  }
	  }
	  
	}
      }//  end loop over the number of elements
    }//  end loop over the number of processors
  }//  end of the master processor part
  else{
    MPI_Send (buff2,buffsize2,MPI_DOUBLE,0,Myrank,MPI_COMM_WORLD);
    //MPI_Send_mine (buff2,buffsize2,MPI_DOUBLE,0,Myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  fflush (stdout);
    
  if (Myrank==0){
    
    //  stiffness matrix assembling
    stiffness_matrix (0);
    
    //  loop over load cases
    for (i=0;i<Lsrs->nlc;i++){
      //  right hand side of system
      mefel_right_hand_side (i,rhs);
    }
    
    Smat->printmat (Out);
    
    //  solution of equation system
    for (i=0;i<Lsrs->nlc;i++){
      Mp->ssle->solve_system (Gtm,Smat,Lsrs->give_lhs(i),Lsrs->give_rhs(i),Out);
    }
    
    print_init(-1, "wt");    
    // clean temperature strains because they are stored in one array used for all load cases
    Mm->nulltempstrains();
    rhs = Lsrs->give_rhs(0);
    for (i=0;i<Lsrs->nlc;i++){
      //  indicator of strain computation
      //  it means, strains at integration points have not been computed
      Mp->strainstate=0;
      //  indicator of stress computation
      //  it means, stresses at integration points have not been computed
      Mp->stressstate=0;
      //  indicator of computation of other array
      //  it means, stresses at integration points have not been computed
      Mp->otherstate=0;
      // left hand side vector for debugging
      lhs =  Lsrs->give_lhs(i);
      // rhs was rewritten by the solver call -> rebuild load vector
      aux = Mp->homog;
      Mp->homog = 0;
      mefel_right_hand_side(i,rhs,fl);
      Mp->homog = aux; 
      //  computes required quantities
      compute_req_val (i); 
      // in artificial time 0.0 the nodal forces contains only the load vector fl
      print_step(i, 0, 0.0, fl);
      
      // if the nodal forces are required for graphical output
      // reaction + load vector are printed in the artificial time 1.0
      if (Outdm->nog.selnforce.st != sel_no)
	{
	  //  indicator of strain computation
	  //  it means, strains at integration points have not been computed
	  Mp->strainstate=0;
	  //  indicator of stress computation
	  //  it means, stresses at integration points have not been computed
	  Mp->stressstate=0;
	  //  indicator of computation of other array
	  //  it means, stresses at integration points have not been computed
	  Mp->otherstate=0;
	  // calculate nodal forces and reactions
	  internal_forces(i, fi);
	  //  computes required quantities
	  compute_req_val (i);
	  //  prints required quantities    
	  // in artificial time 1.0 the nodal forces contains the load vector + reactions
	  Outdm->print_graphics(Outdm->outgr, i, 1.0, 1, fi);
	}
      // clean temperature strains because they are stored in one array used for all load cases
      Mm->nulltempstrains();
    }
    
  }

  fprintf (stdout,"\n procesor %d",Myrank);
  fflush (stdout);
  
  if (Myrank==0){
    print_close();
    
    delete [] fl;  
    delete [] fi;
  }
  
  delete [] buff1;
  delete [] buff2;
}





/**
   function selects appropriate components of higher level of a multiscale
   problem which are sent to lower level
   
   @param llid - lower level id
   @param ncbuff - the number of components in the array buff
   @param buff - array for selected components
   
   TKr 14/07/2022 according to JK
*/
void higher_to_lower_levelm (long llid,long ncbuff,double *buff)
{
  long i,j,nc,*counter,ncomp,dim;
  double area,ar,*aux;
  strastrestate ssst;

  counter = new long [1];
  counter[0]=0;
  
  switch (Mp->mami){
  case elem_domain:{
    //  each element is connected with a microproblem
    
    //  loop over the number of elements
    for (i=0;i<Mt->ne;i++){
      
      if (Gtm->stop->eldom[i]==llid){
	//  only elements connected with the microid-th microscale are assumed
	//  the following function is located in the file MEFEL/SRC/elemswitch
	higher_to_lower_level_elemm (i,counter,buff);
      }
    }//  end of loop over the number of elements
    break;
  }
  case aggregate_domain:{
    //  each element is connected with a microproblem
    //  elements in one aggregate are connected with the same microproblem
    
    //  loop over the number of elements
    for (i=0;i<Mt->ne;i++){
      
      if (Gtt->stop->eldom[i]==llid){
	//  the following function is located in the file MEFEL/SRC/elemswitch
	higher_to_lower_level_elemm (i,counter,buff);
      }
    }//  end of loop over the number of elements
    break;
  }
  case material_aggregate_domain:{
    //  each aggregate (elements of one material) is connected with a microproblem
    //  averaged data is sent to the microproblem
      
    //  dimension of the problem, one, two or three dimensional case
    //  the dimension is defined with respect to the first finite element
    dim = Mt->give_dimension (0);
    //  the number of transported materials
    ncomp = Mt->give_ncomp (0,0);
    aux = new double [ncbuff];

    for (i=0;i<ncbuff;i++){
      buff[i]=0.0;//erase actual position of array buff
    }

    // total area of the aggregate
    ar = 0.0;

    //  loop over the number of elements
    for (i=0;i<Mt->ne;i++){
      
      //  generally, the averaging should take into account the size
      //  of particular finite elements, not only their number
      //  if needed, the function higher_to_lower_level_elem
      //  will return also area or volume of finite elements
      
      
      if (Gtm->stop->eldom[i]==llid){
	//  the following function is located in the file TRFEL/SRC/elemswitcht
	counter[0]=0;
	for (long ijk=0;ijk<ncbuff;ijk++){
	  aux[ijk]=0.0;//erase array
	}
	higher_to_lower_level_elemm (i,counter,aux);//array aux is filled with displacements and their strains 

	//averaging according to element area/volume:
	ssst=Mt->give_ssst (i,0);
	if(ssst == spacestress)//for 3D elements
	  area = Mt->give_volume(i);
	else
	  area = Mt->give_area(i);
	area = fabs(area);
	//fprintf (Out,"\n elem %ld area = %lf",i,area);

	for (j=0;j<ncbuff;j++){
	  buff[j]+=aux[j]*area;
	}
	ar += area;
      }
    }//  end of loop over the number of elements
    
    //fprintf (Out,"\n\n pred prumerovanim buff1:");
    //fprintf (Out,"\n Proc No.= llid = %ld po: \n",llid);
    //for (long ijk=0;ijk<ncbuff;ijk++){
    //  fprintf (Out,"\n buff1 %ld   %le",ijk,buff[ijk]);
    //}
  
    //  averaging of the values
    fprintf (stdout,"\n\n Averaging of values of material aggregate on processor No. %ld",llid);
    for (i=0;i<ncbuff;i++){
      buff[i]=buff[i]/ar;
    }
    
    //fprintf (Out,"\n\n zprumerovane buff1:");
    //fprintf (Out,"\n Proc No.= llid = %ld po: \n",llid);
    //for (long ijk=0;ijk<ncbuff;ijk++){
    //fprintf (Out,"\n buff1 %ld   %le",ijk,buff[ijk]);
    //}
    
    delete [] aux;
    break;
  }
  case eff_aggregate_domain:{
    //  each aggregate is connected with a microproblem
    //  averaged data is sent to the microproblem
    
    aux = new double [ncbuff];
    for (i=0;i<ncbuff;i++){
      aux[i]=0.0;
      buff[i]=0.0;
    }
    
    //  the number of contributions = the number of elements in the aggregate
    nc=0;
    
    //  loop over the number of elements
    for (i=0;i<Mt->ne;i++){
      
      //  generally, the averaging should take into account the size
      //  of particular finite elements, not only their number
      //  if needed, the function higher_to_lower_level_elem
      //  will return also area or volume of finite elements
      
      
      if (Gtm->stop->eldom[i]==llid){
	//  the following function is located in the file MEFEL/SRC/elemswitch
	higher_to_lower_level_elemm (i,counter,buff);
	counter[0]=0;
	//nc++;

	//averaging according to element area/volume:
	ssst=Mt->give_ssst (i,0);
	if(ssst == spacestress)//for 3D elements
	  area = Mt->give_volume(i);
	else
	  area = Mt->give_area(i);
	area = fabs(area);
	for (j=0;j<ncbuff;j++){
	  aux[j]+=buff[j]/area;
	}
      }
    }//  end of loop over the number of elements
    
    //  averaging of the values
    //for (i=0;i<nc;i++){
    //buff[i]=aux[i]/nc;
    //}
    
    delete [] aux;
    break;
  }
  default:{
    par_print_err (Myrank,"unknown type of aggregation is required",__FILE__,__LINE__,__func__);
  }
  }
  
  delete [] counter;
  
}
