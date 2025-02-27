#include <mpi.h>
#include "hptrftiles.h"
#include "hpglobal.h"
#include "globalt.h"
#include "globmatt.h"
#include "seqfilest.h"

#include <string.h>

/**
   
   JK, 8. 7. 2014
*/
void parallel_trfel_tiles ()
{
  long i,j,ii,combmin,combmax,ncomb,tid;
  double *lhs,*rhs,**condm,**condvlhs,**condvrhs;
  precond prec;  
  char name[1100],nameout[1100],namegid[1100],namedat[1100];
  
  strcpy (nameout,Outdt->outfn);
  strcpy (namegid,Outdt->outgrfn);
  strcpy (namedat,Outdt->outdiagfn);

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

  //  stiffness matrices will be stored there
  gmatrix *smatbackup;
  smatbackup = new gmatrix [Ntiletypes];
  
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
    
    //  backup of the matrices
    SSmat[i].copygm (smatbackup[i]);

    //  condensation - elimination of all internal DOFs
    SSmat[i].condense (Gtt,condm[i],condvrhs[i],Slhs[i],Srhs[i],Nbdof[i],1,Outt);
  }


  Ndoft=0;
  for (i=0;i<Ntiles;i++){
    for (j=0;j<Nbdoftiling[i];j++){
      if (Ndoft<Cnbn[i][j])
        Ndoft=Cnbn[i][j];
    }
  }
  
  gtopology *top;
  top = new gtopology;
  top->initiate (Cnbn,Nbdoftiling,Ntiles);


  if (Kmat!=NULL)
    delete Kmat;
  
  if (Tp!=NULL)
    delete Tp;
  Tp = new probdesct;
  
  Tp->gradcomp = 1;
  Tp->gradpos  = 2;
  Tp->gradaver = 1;
  Tp->fluxcomp = 1;
  Tp->fluxpos  = 2;
  Tp->fluxaver = 1;
  Tp->othercomp = 1;
  Tp->otherpos  = 2;
  Tp->otheraver = 1;
  Tp->eqothercomp = 1;
  Tp->eqotherpos  = 2;
  Tp->eqotheraver = 1;
  
  //Tp->ssle->tlinsol=ldl;
  Tp->ssle->tlinsol=cg;
  Tp->ssle->ni=1000;
  Tp->ssle->res=1.0e-8;
  Tp->ssle->iv=0;
  Tp->ssle->prec.pt=noprecond;
  //Tp->tstorkm=skyline_matrix;
  Tp->tstorkm=symm_comp_rows;
  Tp->tprob=stationary_problem;
  
  Kmat = new gmatrix;
  Kmat->setval (Tp->ssle);
  Kmat->initiate (top,Ndoft,Tp->tstorkm,1,Outt);


  fprintf (stderr,"\n\n\n\n\n pocet DOFs %ld",Ndoft);
  fprintf (stderr,"\n ukladani %d",Kmat->ts);
  fprintf (stderr,"\n solver Kmat %d",Kmat->tlinsol);
  fprintf (stderr,"\n solver Mp   %d",Tp->ssle->tlinsol);

  lhs = new double [Ndoft];
  rhs = new double [Ndoft];
  
  ncomb=Ncomb/Nproc;
  combmin=Myrank*ncomb;
  combmax=(Myrank+1)*ncomb;

  fprintf (stderr,"\n myrank %d   ncomb %ld  combmin %ld   combmax %ld",Myrank,ncomb,combmin,combmax);


  for (i=combmin;i<combmax;i++){

    //Kmat->sky->nullsky ();
    nullv (rhs,Ndoft);

    for (j=0;j<Ntiles;j++){
      //  tile id
      tid = Tiling[i][j];
      fprintf (stderr,"\n myrank %d   i %ld    tid %ld",Myrank,i,tid);
     
      //  the last parameter is useless
      Kmat->localized (condm[tid],Cnbn[j],Nbdoftiling[j],0);


      
      //  localization of the right hand side
      //for (ii=0;ii<Nbdoftiling[j];ii++){
      //if (Cnbn[j][ii]<0){
      //condvlhs[tid][ii]=Pdofval[0-Cnbn[j][ii]-1];
      //}else{
      //condvlhs[tid][ii]=0.0;
      //}
      //}
      
      //  new localization of the right hand side
      for (ii=0;ii<Nbdoftiling[j];ii++){
        if (Cnbn[j][ii]<0){
          Slhs[tid][ii+Nidof[tid]]=Pdofval[0-Cnbn[j][ii]-1];
        }else{
          Slhs[tid][ii+Nidof[tid]]=0.0;
        }
      }
      for (ii=0;ii<Nidof[tid];ii++){
	Slhs[tid][ii]=0.0;
      }
      

      //for (ii=0;ii<Nidof[tid]+Nbdoftiling[j];ii++){
      //fprintf (stderr,"tile lhs  %6ld   %le\n",tid,Slhs[j][ii]);
      //}
      
      //  toto vyhodit
      //for (ii=0;ii<Nbdoftiling[j];ii++){
      //s=0.0;
      //for (jj=0;jj<Nbdoftiling[j];jj++){
      //s+=condm[tid][ii*Nbdoftiling[j]+jj]*condvlhs[tid][jj];
      //}
      //condvrhs[tid][ii]=s;
      //}
      
      //  zde musi byt nasobeni matice s vektorem predepsanych hodnot
      smatbackup[tid].gmxv (Slhs[tid],Srhs[tid]);
      cmulv(-1.0,Srhs[tid],Nidof[tid]+Nbdoftiling[j]);
      
      //for (ii=0;ii<Nidof[tid]+Nbdoftiling[j];ii++){
      //fprintf (stderr,"tile rhs  %6ld   %le\n",tid,Srhs[j][ii]);
      //}

      //  modification of the right hand side
      SSmat[tid].condense (Gtt,condm[tid],condvrhs[tid],Slhs[tid],Srhs[tid],Nbdof[tid],4,Outt);

      //for (ii=0;ii<Nbdoftiling[j];ii++){
      //fprintf (stderr,"tile condvrhs  %6ld   %le\n",tid,condvrhs[j][ii]);
      //}

      locglob (rhs,condvrhs[tid],Cnbn[j],Nbdoftiling[j]);
      
    }//  end of for (j=0;j<Ntiles;j++){
    

    //long kij=0;
    //for (long ijk=0;ijk<Kmat->scr->n;ijk++){
    //for (long jik=Kmat->scr->adr[ijk];jik<Kmat->scr->adr[ijk+1];jik++){
    //  fprintf (stdout,"\n Kmat %ld  %ld %le",ijk,Kmat->scr->ci[jik],Kmat->scr->a[kij]);
    //  kij++;
    //}
    //}
    //long ijk=0;
    //for (ijk=0;ijk<Ndoft;ijk++){
    //fprintf (stdout,"rhs  %6ld   %le\n",ijk,rhs[ijk]);
    //}
    

    Tp->ssle->solve_system (top,Kmat,lhs,rhs,Outt);
    //Kmat->solve_system (top,*Tp->ssle->prec.pt,lhs,rhs,Outt);
    //fprintf (stderr,"\n rhs %le     Kmat %le  Ndoft %ld",rhs[Ndoft-1],Kmat->sky->a[0],Ndoft);
    //for (ii=0;ii<Ndoft;ii++){
    //fprintf (stderr,"myrank %ld master lhs  %6ld   %15.12le\n",Myrank,ii,lhs[ii]);
    //}
    
    //  id in tiling
    //long did=1;
    long did;
    
    for (did=0;did<Ntiles;did++){
      
      //  tile id
      tid = Tiling[i][did];
      //fprintf (stderr,"\n myrank %ld   i %ld    tid %ld\n",Myrank,i,tid);
      //globloc (lhs,condvlhs[tid],Cnbn[0],Nbdoftiling[0]);
      
      for (ii=0;ii<Nbdoftiling[did];ii++){
        if (Cnbn[did][ii]>0)
          condvlhs[tid][ii]=lhs[Cnbn[did][ii]-1];
        if (Cnbn[did][ii]<0)
          condvlhs[tid][ii]=Pdofval[0-Cnbn[did][ii]-1];
        //fprintf (stdout,"condvlhs %ld  %le\n",ii,condvlhs[tid][ii]);
      }
      //fprintf (stdout,"\n condm %le %le %le %le\n",condm[tid][0],condm[tid][1],condm[tid][2],condm[tid][3]);
      
      nullv (Srhs[tid],Nidof[tid]+Nbdoftiling[did]);
      //  back substitution
      SSmat[tid].condense (top,condm[tid],condvlhs[tid],Slhs[tid],Srhs[tid],Nbdoftiling[did],2,Outt);
      
      //fprintf (stderr,"\n myrank %ld   Nidof[tid]+Nbdoftiling[4]  %ld",Myrank,Nidof[tid]+Nbdoftiling[did]);
      for (j=0;j<Nidof[tid]+Nbdoftiling[did];j++){
        //Lsrst->lhs[j]=Slhs[0][j];
      }
      
      for (ii=0;ii<Ntiles;ii++){
        //fprintf (stdout,"Nbdoftiling %2ld  %ld\n",ii,Nbdoftiling[ii]);
      }
      for (ii=0;ii<Ntiletypes;ii++){
        //fprintf (stdout,"Nbdof       %2ld  %ld\n",ii,Nbdof[ii]);
      }
      
      //  names of output files
      strcpy (name,nameout);
      sprintf (&name[strlen (nameout)],"%ld.out",i);
      strcpy (Outdt->outfn, name);
      
      if (Outdt->gf != grftt_no){
        strcpy (name,namegid);
        sprintf (&name[strlen (namegid)],"%d",Myrank+1);
        strcpy (Outdt->outgrfn, name);
      }
      if (Outdt->ndiag > 0){
        strcpy (name,namedat);
        sprintf (&name[strlen (namedat)],"%d.dat",Myrank+1);
        strcpy (Outdt->outdiagfn, name);
      }  
      
      FILE *out;
      char name[100];
      sprintf (name,"dlazdice%d.out",Myrank+1);
      
      out = fopen (name,"w");
      for (j=0;j<Nidof[tid]+Nbdoftiling[did];j++){
        fprintf (out,"%7ld   % 20.15le\n",j,Lsrst->lhs[j]);
        //fprintf (out,"%7ld   % 20.15le\n",j,Slhs[did][j]);
      }
      //SSmat[tid].printmat (out);
      fclose (out);
      
      //if (Tp!=NULL)
      //delete Tp;
      //Tp = new probdesct;
      //Tp->gradcomp = 1;
      //Tp->fluxcomp = 1;
      //  fluxes are in integration points
      //Tp->fluxpos  = 1;
      
      
      //compute_req_valt (0);
      
      //print_initt(-1, "wt");    
      //print_stept(i, 0, 0.0, rhs);
      //print_closet();
      
      
    }//  end of for (did=0;did<Ntiles;did++){
  }//  end of for (i=combmin;i<combmax;i++){


}


/**
   function solves problems with tiles
   tile types are distributed on processors
   solution of reduced system of equations is done by the conjugate gradient method
   %matrix of the reduced system is not assembled
   
   JK, 9. 9. 2014
*/
void parallel_trfel_tiles_parsolver (XFILE *in)
{
  long i,j,k,m,ii,ijk,tndoft,ndofcp,sizebuff,gndofe,ntil;
  long npdof;
  long nicgsch,anicgsch;
  long *cn,*tileincidence,*tiling,*pdof;
  long **tileposition;
  long sizelbuff,*lbuff;
  double norrhs,nom,denom,alpha,beta,errcgsch,aerrcgsch,zero;
  double *lhs,*rhs,*condmat,*condvlhs,*condvrhs,*r,*d,*p,*buff;
  char nameout[1100],namegid[1100],namedat[1100];
  gtopology *gtop;
  MPI_Status stat;
  FILE *out1,*out2,*out3,*out4,*out5;
  time_t t1,t2,t3,t4,t5,t6,t7;
  FILE *casy;
  
  if (Myrank==0){
    casy = fopen ("casy.txt","w");
  }
  
  sizelbuff=6;
  lbuff = new long [sizelbuff];
  
  t1 = time (NULL);

  gtop = NULL;
 
  /*
  out1 = fopen ("out1","w");
  out2 = fopen ("out2","w");
  out3 = fopen ("out3","w");
  out4 = fopen ("out4","w");
  out5 = fopen ("out5","w");
  */

  nicgsch = 100;
  errcgsch = 1.0e-6;
  zero = 1.0e-15;
  
  strcpy (nameout,Outdt->outfn);
  strcpy (namegid,Outdt->outgrfn);
  strcpy (namedat,Outdt->outdiagfn);

  // *******************************************************************************************
  // *******************************************************************************************
  //   each processor assembles, factorizes and stores matrices of a single tile type
  // *******************************************************************************************
  // *******************************************************************************************
  
  //  backup of the conductivity matrix will be needed
  gmatrix *kmatbackup;
  kmatbackup = new gmatrix [1];
  
  //  array for the condensed matrices
  condmat = new double [Nbdof[Myrank]*Nbdof[Myrank]];
  memset (condmat,0,(Nbdof[Myrank]*Nbdof[Myrank])*sizeof(double));
  
  //  array for the condensed vectors of the solution
  condvlhs = new double [Nbdof[Myrank]];
  memset (condvlhs,0,Nbdof[Myrank]*sizeof(double));
  
  //  array for the condensed vectors of the right hand side
  condvrhs = new double [Nbdof[Myrank]];
  memset (condvrhs,0,Nbdof[Myrank]*sizeof(double));
  
  t2 = time (NULL);

  //  assembling of conductivity matrix of a tile
  conductivity_matrix (0);
  
  //  backup of the matrices
  Kmat->copygm (*kmatbackup);
  
  t3 = time (NULL);

  //  condensation - elimination of all internal DOFs
  Kmat->condense (Gtt,condmat,condvrhs,Lsrst->lhs,Lsrst->rhs,Nbdof[Myrank],1,Outt);

  t4 = time (NULL);

  // ************************************
  // ************************************
  //  solution of a sequence of tilings
  //  the number of tiles in a tiling may be variable
  // ************************************
  // ************************************

  fprintf (stdout,"\n\n pocet kombinaci %ld\n\n",Ncomb);
  for (ijk=0;ijk<Ncomb;ijk++){
    fprintf (stdout,"\n\n kombinace %ld\n\n",ijk);
    
    //  the number of tiles in a tiling
    xfscanf (in,"%ld",&Ntiles);
    
    //  the numbers of boundary DOFs on tiles in tiling
    Nbdoftiling = new long [Ntiles];
    for (i=0;i<Ntiles;i++){
      xfscanf (in,"%ld",Nbdoftiling+i);
    }
    
    // *************************************************
    //  the code numbers of Schur complements on tiles
    // *************************************************
    //  the total number of DOFs in the reduced problem
    tndoft=0;
    Cnbn = new long* [Ntiles];
    for (i=0;i<Ntiles;i++){
      Cnbn[i] = new long [Nbdoftiling[i]];
      
      for (j=0;j<Nbdoftiling[i];j++){
	xfscanf (in,"%ld",&Cnbn[i][j]);
	if (tndoft<Cnbn[i][j])
	  tndoft = Cnbn[i][j];
      }
    }
      
    gtop = new gtopology;
    
    //  numbers of unknowns (code numbers) are copied to the instant
    //  of the class gtopology
    gtop->initiate (Cnbn,Nbdoftiling,Ntiles);
    
    //  id of tiles used in a tiling
    tiling = new long [Ntiles];
    for (i=0;i<Ntiles;i++){
      xfscanf (in,"%ld",tiling+i);
    }
    
    //  the array tiling could contain a tile number several times
    //  for example, tiling = {3,4,3,4,6,2,0,6}	
    //  in the array tileincidence, the following numbers are stored in the case of 12 types of tiles
    //  tileincidence = {1,0,1,2,2,0,2,0,0,0,0,0}
    //  tileposition = {{6};{x};{5};{0,2};{1,3};{x};{4,7};{x};{x};{x};{x};{x}}
    tileincidence = new long [Ntiletypes];
    for (i=0;i<Ntiletypes;i++){
      tileincidence[i]=0;
    }
    for (i=0;i<Ntiles;i++){
      tileincidence[tiling[i]]++;
    }
    
    tileposition = new long* [Ntiletypes];
    for (i=0;i<Ntiletypes;i++){
      tileposition[i] = new long [tileincidence[i]];
    }
    
    for (i=0;i<Ntiletypes;i++){
      j=0;
      for (k=0;k<Ntiles;k++){
	if (i==tiling[k]){
	  tileposition[i][j]=k;
	  j++;
	}
      }
    }
    
    //  size of buffer
    sizebuff=0;
    for (i=0;i<Ntiletypes;i++){
      if (sizebuff < Nbdof[i]*tileincidence[i])
	sizebuff = Nbdof[i]*tileincidence[i];
    }
    
    //  size of buffer
    //  the buffer contains vectors of particular tiles
    //  at the end, there must be a double number for iteration managing
    sizebuff++;
    //  the number of contributions from tile assigned to the processor
    ntil=tileincidence[Myrank];
    
    if (ntil>0){
      Lsrst->lhs = new double [Ndoft*ntil];
      Lsrst->rhs = new double [Ndoft*ntil];
      
      memset (Lsrst->lhs,0,(Ndoft*ntil)*sizeof(double));
      memset (Lsrst->rhs,0,(Ndoft*ntil)*sizeof(double));
    }
    
    buff = new double [sizebuff];
    memset (buff,0,sizebuff*sizeof(double));
    
    // *********************************
    //  reading of boundary conditions
    // *********************************
    //  the number of prescribed DOFs
    xfscanf (in,"%ld",&npdof);
    pdof = new long [npdof];
    Pdofval = new double [npdof];
    
    //  prescribed values
    for (i=0;i<npdof;i++){
      xfscanf (in,"%ld %le",pdof+i,Pdofval+i);
    }
    
    
    //  loop over the number of contributions from the tile assigned to the processor
    for (i=0;i<ntil;i++){
      k=tileposition[Myrank][i];
      
      //  new localization of the right hand side
      for (j=0;j<Nbdoftiling[Myrank];j++){
        if (Cnbn[k][j]<0){
          Lsrst->lhs[i*Ndoft+j+Nidof[Myrank]]=Pdofval[0-Cnbn[k][j]-1];
        }else{
          Lsrst->lhs[i*Ndoft+j+Nidof[Myrank]]=0.0;
        }
      }
      for (j=0;j<Nidof[Myrank];j++){
        Lsrst->lhs[i*Ndoft+j]=0.0;
      }
      
      /*
      if (Myrank==0){
        fprintf (out1,"\n lhs");
        for (j=0;j<Ndoft;j++){
          fprintf (out1,"%le ",Lsrst->lhs[i*Ndoft+j]);
        }
        fprintf (out1,"\n negm %ld ",kmatbackup->sky->negm);
        for (j=0;j<kmatbackup->sky->negm;j++){
          fprintf (out1,"%le ",kmatbackup->sky->a[j]);
        }
      }
      if (Myrank==1){
        fprintf (out2,"\n lhs");
        for (j=0;j<Ndoft;j++){
          fprintf (out2,"%le ",Lsrst->lhs[i*Ndoft+j]);
        }
        fprintf (out2,"\n negm %ld ",kmatbackup->sky->negm);
        for (j=0;j<kmatbackup->sky->negm;j++){
          fprintf (out2,"%le ",kmatbackup->sky->a[j]);
        }
      }
      if (Myrank==2){
        fprintf (out3,"\n lhs");
        for (j=0;j<Ndoft;j++){
          fprintf (out3,"%le ",Lsrst->lhs[i*Ndoft+j]);
        }
        fprintf (out3,"\n negm %ld ",kmatbackup->sky->negm);
        for (j=0;j<kmatbackup->sky->negm;j++){
          fprintf (out3,"%le ",kmatbackup->sky->a[j]);
        }
      }
      */
      
      kmatbackup->gmxv (Lsrst->lhs+i*Ndoft,Lsrst->rhs+i*Ndoft);
      cmulv(-1.0,Lsrst->rhs+i*Ndoft,Ndoft);
      
      /*
      if (Myrank==0){
        fprintf (out1,"\n rhs");
        for (j=0;j<Ndoft;j++){
          fprintf (out1,"%le ",Lsrst->rhs[i*Ndoft+j]);
        }
      }
      if (Myrank==1){
        fprintf (out2,"\n rhs");
        for (j=0;j<Ndoft;j++){
          fprintf (out2,"%le ",Lsrst->rhs[i*Ndoft+j]);
        }
      }
      if (Myrank==2){
        fprintf (out3,"\n rhs");
        for (j=0;j<Ndoft;j++){
          fprintf (out3,"%le ",Lsrst->rhs[i*Ndoft+j]);
        }
      }
      */

      Kmat->condense (Gtt,condmat,condvrhs,Lsrst->lhs+i*Ndoft,Lsrst->rhs+i*Ndoft,Nbdof[Myrank],4,Outt);
      
      /*
      if (Myrank==0){
        fprintf (out1,"\n rhs");
        for (j=0;j<Ndoft;j++){
          fprintf (out1,"%le ",Lsrst->rhs[i*Ndoft+j]);
        }
      }
      if (Myrank==1){
        fprintf (out2,"\n rhs");
        for (j=0;j<Ndoft;j++){
          fprintf (out2,"%le ",Lsrst->rhs[i*Ndoft+j]);
        }
      }
      if (Myrank==2){
        fprintf (out3,"\n rhs");
        for (j=0;j<Ndoft;j++){
          fprintf (out3,"%le ",Lsrst->rhs[i*Ndoft+j]);
        }
      }
      */
      
      for (j=0;j<Nbdof[Myrank];j++){
        buff[i*Nbdof[Myrank]+j]=Lsrst->rhs[i*Ndoft+j+Nidof[Myrank]];
      }
    }//  end of the loop for (i=0;i<ntil;i++){
    
    /*
    if (Myrank==0){
      fprintf (out1,"\n buff");
      for (j=0;j<sizebuff;j++){
        fprintf (out1,"%le ",buff[j]);
      }
    }
    if (Myrank==1){
      fprintf (out2,"\n buff");
      for (j=0;j<sizebuff;j++){
        fprintf (out2,"%le ",buff[j]);
      }
    }
    if (Myrank==2){
      fprintf (out3,"\n buff");
      for (j=0;j<sizebuff;j++){
        fprintf (out3,"%le ",buff[j]);
      }
    }
    if (Myrank==3){
      fprintf (out4,"\n buff, sizebuff %ld\n",sizebuff);
      for (j=0;j<sizebuff;j++){
        fprintf (out4,"%le ",buff[j]);
      }
    }
    */
      
    MPI_Barrier (MPI_COMM_WORLD);
    
    if (Myrank==0){
      
      lhs = new double [tndoft];
      rhs = new double [tndoft];
      for (i=0;i<tndoft;i++){
        lhs[i]=0.0;
        rhs[i]=0.0;
      }
      
      //  contribution from the master processor
      m=0;
      //  loop over the number of contributions from the m-th processor/tile
      for (j=0;j<tileincidence[m];j++){
	//  position of the m-th tile in the tiling
        k=tileposition[m][j];
	//  the number of generalized code numbers of the k-th tile location in the tiling
        gndofe = gtop->give_ndofe (k);
        cn = new long [gndofe];
	//  generalized code numbers of the k-th tile location in the tiling
        gtop->give_code_numbers (k,cn);
	//  localization of the contribution to the rhs array
        locglob (rhs,buff+j*gndofe,cn,gndofe);
        delete [] cn;
      }
      
      //  contributions from slave processors
      for (i=1;i<Ntiletypes;i++){
        MPI_Recv(buff,sizebuff,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	//  processor/tile id where the data is from
        m=stat.MPI_TAG;
	//  loop over the number of contributions from the m-th processor/tile
        for (j=0;j<tileincidence[m];j++){
	  //  position of the m-th tile in the tiling
          k=tileposition[m][j];
	  //  the number of generalized code numbers of the k-th tile location in the tiling
          gndofe = gtop->give_ndofe (k);
          cn = new long [gndofe];
	  //  generalized code numbers of the k-th tile location in the tiling
          gtop->give_code_numbers (k,cn);
	  //  localization of the contribution to the rhs array
          locglob (rhs,buff+j*gndofe,cn,gndofe);
          delete [] cn;
        }
      }
      
      
    }//  end of the statement if (Myrank==0){
    else{
      MPI_Send(buff,sizebuff,MPI_DOUBLE,0,Myrank,MPI_COMM_WORLD);
    }//  end of the statement else
    MPI_Barrier (MPI_COMM_WORLD);
    
    
    // *******************************************************************
    // *******************************************************************
    //  solution of the reduced problem by the conjugate gradient method
    // *******************************************************************
    // *******************************************************************
    
    t5 = time (NULL);

    ndofcp = tndoft;
      
    if (Myrank==0){
      
      /*
      if (Myrank==0){
        
        fprintf (out1,"\n ndofcp %ld",ndofcp);
        fprintf (out1,"\n rhs\n");
        for (j=0;j<ndofcp;j++){
          fprintf (out1,"%le ",rhs[j]);
        }
      }
      
      fprintf (stdout,"\n\n myrank %ld  ntiles %ld\n",Myrank,Ntiles);
      for (i=0;i<Ntiles;i++){
        fprintf (stdout,"%ld ",Nbdoftiling[i]);
      }
      
      fprintf (stdout,"\n cnbn");
      for (i=0;i<Ntiles;i++){
        fprintf (stdout,"\n");
        for (j=0;j<Nbdoftiling[i];j++){
          fprintf (stdout,"%ld ",Cnbn[i][j]);
        }
      }
      
      fprintf (stdout,"\n tileincidence\n");
      for (i=0;i<Ntiletypes;i++){
        fprintf (stdout,"%ld ",tileincidence[i]);
      }
      
      fprintf (stdout,"\n tileposition\n");
      for (i=0;i<Ntiletypes;i++){
        for (j=0;j<tileincidence[i];j++){
          fprintf (stdout,"%ld ",tileposition[i][j]);
        }
      }
      */
      

      fprintf (stdout,"\n\n JSME TADY ndofcp %ld\n\n",ndofcp);

      
      r = new double [ndofcp];
      d = new double [ndofcp];
      p = new double [ndofcp];
      
      //  iteration initialization
      norrhs=0.0;  nom=0.0;
      for (i=0;i<ndofcp;i++){
        lhs[i]=0.0;
        r[i]=rhs[i];
        d[i]=rhs[i];
        norrhs+=rhs[i]*rhs[i];
        nom+=r[i]*r[i];
      }
      fprintf (stdout,"\n norrhs  %e",norrhs);
      
      //  iteration loop
      for (ii=0;ii<nicgsch;ii++){
        
        //  direction vector scattering
        for (i=1;i<Ntiletypes;i++){
          memset (buff,0,sizebuff*sizeof(double));
          for (j=0;j<tileincidence[i];j++){
            k=tileposition[i][j];
            gndofe = gtop->give_ndofe (k);
            cn = new long [gndofe];
            gtop->give_code_numbers (k,cn);
            globloc (d,buff+j*gndofe,cn,gndofe);
            delete [] cn;
          }
          MPI_Send(buff,sizebuff,MPI_DOUBLE,i,Myrank,MPI_COMM_WORLD);
        }
        
        i=0;
        memset (buff,0,sizebuff*sizeof(double));
        for (j=0;j<tileincidence[i];j++){
          k=tileposition[i][j];
          gndofe = gtop->give_ndofe (k);
          cn = new long [gndofe];
          gtop->give_code_numbers (k,cn);
          globloc (d,buff+j*gndofe,cn,gndofe);
          delete [] cn;
        }
        
        //  matrix-vector multiplication collection
        memset (p,0,ndofcp*sizeof(double));
        
        for (i=0;i<ntil;i++){
          //  matrix-vector multiplication
          mxv (condmat,buff+i*Nbdof[Myrank],condvlhs,Nbdof[Myrank],Nbdof[Myrank]);
          copyv (condvlhs,buff+i*Nbdof[Myrank],Nbdof[Myrank]);
        }
        
        i=0;
        for (j=0;j<tileincidence[i];j++){
          k=tileposition[i][j];
          gndofe = gtop->give_ndofe (k);
          cn = new long [gndofe];
          gtop->give_code_numbers (k,cn);
          locglob (p,buff+j*gndofe,cn,gndofe);
          delete [] cn;
        }
        

        for (i=1;i<Ntiletypes;i++){
          MPI_Recv(buff,sizebuff,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
          m=stat.MPI_TAG;
          for (j=0;j<tileincidence[m];j++){
            k=tileposition[m][j];
            gndofe = gtop->give_ndofe (k);
            cn = new long [gndofe];
            gtop->give_code_numbers (k,cn);
            locglob (p,buff+j*gndofe,cn,gndofe);
            delete [] cn;
          }
        }
	
        //  denominator
        denom=ss(d,p,ndofcp);
        
        if (fabs(denom)<zero){
          if (ii!=nicgsch-1){
            buff[sizebuff-1]=1.0;
            for (j=1;j<Nproc;j++){
              MPI_Send(buff,sizebuff,MPI_DOUBLE,j,Myrank,MPI_COMM_WORLD);
            }
          }
          fprintf (stderr,"\n\n denominator in alpha expression in parallel conjugate gradient method is equal to zero.\n");
          break;
        }
        
        //  coefficient alpha
        alpha=nom/denom;
        
        //  new global vectors
        for (j=0;j<ndofcp;j++){
          r[j]-=alpha*p[j];
          lhs[j]+=alpha*d[j];
        }
        
        
        denom=nom;
        nom=ss(r,r,ndofcp);
        
        //if (i%100==0)  fprintf (stdout,"\n iteration %ld    norres/norrhs %e",i,nom/norrhs);
        fprintf (stdout,"\n iteration %ld    norres %e   norrhs %e   norres/norrhs %e",ii,nom,norrhs,nom/norrhs);
        if (nom/norrhs<errcgsch){
          fprintf (stdout,"\n last iteration %ld    norres/norrhs %e",ii,nom/norrhs);
          if (ii!=nicgsch-1){
            buff[sizebuff-1]=1.0;
            for (j=1;j<Nproc;j++){
              MPI_Send(buff,sizebuff,MPI_DOUBLE,j,Myrank,MPI_COMM_WORLD);
            }
          }
          break;
        }
        
        //  beta coefficient
        beta=nom/denom;
        
        //  new direction vector
        for (j=0;j<ndofcp;j++){
          d[j]=r[j]+beta*d[j];
        }
        
      }//  end of the statement for (ii=0;ii<nicgsch;ii++){
      
      anicgsch=i;
      aerrcgsch=nom/norrhs;
      
      delete [] r;
      delete [] d;
      delete [] p;
      
      delete [] rhs;
    }//  end of the statement if (Myrank==0)
    else{
      //  iteration loop
      for (i=0;i<nicgsch;i++){
        //  direction vector receipt
        MPI_Recv(buff,sizebuff,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        
        if (buff[sizebuff-1]>0.5)  break;
        
        for (j=0;j<ntil;j++){
          //  matrix-vector multiplication
          mxv (condmat,buff+j*Nbdof[Myrank],condvlhs,Nbdof[Myrank],Nbdof[Myrank]);
          copyv (condvlhs,buff+j*Nbdof[Myrank],Nbdof[Myrank]);
        }
        
        MPI_Send(buff,sizebuff,MPI_DOUBLE,0,Myrank,MPI_COMM_WORLD);
        
      }
      anicgsch=i;
    }//  end of the statement else
    MPI_Barrier (MPI_COMM_WORLD);
  
    
    // *******************************************************************
    // *******************************************************************
    //  back substitution on tiles
    // *******************************************************************
    // *******************************************************************
    
    
    if (Myrank==0){
      
      for (i=1;i<Ntiletypes;i++){
        memset (buff,0,sizebuff*sizeof(double));
        for (j=0;j<tileincidence[i];j++){
          k=tileposition[i][j];
          gndofe = gtop->give_ndofe (k);
          cn = new long [gndofe];
          gtop->give_code_numbers (k,cn);
          globloc (lhs,buff+j*gndofe,cn,gndofe);
          delete [] cn;
        }
        MPI_Send(buff,sizebuff,MPI_DOUBLE,i,Myrank,MPI_COMM_WORLD);
      }
      
      i=0;
      memset (buff,0,sizebuff*sizeof(double));
      for (j=0;j<tileincidence[i];j++){
        k=tileposition[i][j];
        gndofe = gtop->give_ndofe (k);
        cn = new long [gndofe];
        gtop->give_code_numbers (k,cn);
        globloc (lhs,buff+j*gndofe,cn,gndofe);
        delete [] cn;
      }
      
      delete [] lhs;
    }
    else{
      MPI_Recv(buff,sizebuff,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }//  end of the statement else
    MPI_Barrier (MPI_COMM_WORLD);
    
    t6 = time (NULL);

    for (i=0;i<ntil;i++){
      for (j=0;j<Nbdof[Myrank];j++){
        condvlhs[j]=buff[i*Nbdof[Myrank]+j];
      }
      
      //  back substitution
      Kmat->condense (Gtt,condmat,condvlhs,Lsrst->lhs+i*Ndoft,Lsrst->rhs+i*Ndoft,Nbdof[Myrank],2,Outt);
      
      /*
      fprintf (Outt,"\n\n\n Vysledky \n\n");
      for (k=0;k<Ndoft;k++){
	fprintf (Outt,"\n Lsrst %6ld   % 20.14le",k,Lsrst->lhs[i*Ndoft+k]);
      }
      fprintf (Outt,"\n");
      */

    }
    
    t7 = time (NULL);
      
    print_initt(-1, "wt");
    compute_req_valt (0);
    print_stept (0, 0, 0.0, NULL);
    

    delete [] pdof;
    delete [] Pdofval;

    delete [] buff;

    if (ntil>1){
      delete [] Lsrst->lhs;
      delete [] Lsrst->rhs;
    }

    for (i=0;i<Ntiletypes;i++){
      delete [] tileposition[i];
    }
    delete [] tileposition;
    
    delete [] tileincidence;
    
    delete [] tiling;
    
    delete gtop;

    for (i=0;i<Ntiles;i++)
      delete [] Cnbn[i];
    delete [] Cnbn;
    
    delete [] Nbdoftiling;

    

    //  array allocation
    lbuff[0] = t2-t1;
    //  matrix assembling
    lbuff[1] = t3-t2;
    //  condensation
    lbuff[2] = t4-t3;
    //  reading and assembling of new tiling
    lbuff[3] = t5-t4;
    //  conjugate gradient on master
    lbuff[4] = t6-t5;
    //  back substitution on tiling
    lbuff[5] = t7-t6;
   
    if (Myrank==0){
      
      fprintf (casy,"procesor %3ld %3ld\n",0,ijk);
      fprintf (casy,"array allocation    %ld\n",lbuff[0]);
      fprintf (casy,"matrix assembling   %ld\n",lbuff[1]);
      fprintf (casy,"condensation        %ld\n",lbuff[2]);
      fprintf (casy,"reading and assembling of new tiling   %ld\n",lbuff[3]);
      fprintf (casy,"conjugate gradient on master           %ld\n",lbuff[4]);
      fprintf (casy,"back substitution on tiling            %ld\n",lbuff[5]);
      fprintf (casy,"\n\n\n");
      
      fprintf (stdout,"procesor %3ld %3ld\n",0,ijk);
      fprintf (stdout,"array allocation    %ld\n",lbuff[0]);
      fprintf (stdout,"matrix assembling   %ld\n",lbuff[1]);
      fprintf (stdout,"condensation        %ld\n",lbuff[2]);
      fprintf (stdout,"reading and assembling of new tiling   %ld\n",lbuff[3]);
      fprintf (stdout,"conjugate gradient on master           %ld\n",lbuff[4]);
      fprintf (stdout,"back substitution on tiling            %ld\n",lbuff[5]);
      fprintf (stdout,"\n\n\n");

      for (i=1;i<Nproc;i++){
        MPI_Recv(lbuff,sizelbuff,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	//  processor/tile id where the data is from
        m=stat.MPI_TAG;
        
        fprintf (casy,"procesor %3ld %3ld\n",m,ijk);
        fprintf (casy,"array allocation    %ld\n",lbuff[0]);
        fprintf (casy,"matrix assembling   %ld\n",lbuff[1]);
        fprintf (casy,"condensation        %ld\n",lbuff[2]);
        fprintf (casy,"reading and assembling of new tiling   %ld\n",lbuff[3]);
        fprintf (casy,"conjugate gradient on master           %ld\n",lbuff[4]);
        fprintf (casy,"back substitution on tiling            %ld\n",lbuff[5]);
        fprintf (casy,"\n\n\n");
        
        fprintf (stdout,"procesor %3ld %3ld\n",m,ijk);
        fprintf (stdout,"array allocation    %ld\n",lbuff[0]);
        fprintf (stdout,"matrix assembling   %ld\n",lbuff[1]);
        fprintf (stdout,"condensation        %ld\n",lbuff[2]);
        fprintf (stdout,"reading and assembling of new tiling   %ld\n",lbuff[3]);
        fprintf (stdout,"conjugate gradient on master           %ld\n",lbuff[4]);
        fprintf (stdout,"back substitution on tiling            %ld\n",lbuff[5]);
        fprintf (stdout,"\n\n\n");
        
      }
      
    }
    else{
      MPI_Send(lbuff,sizelbuff,MPI_LONG,0,Myrank,MPI_COMM_WORLD);
    }
    MPI_Barrier (MPI_COMM_WORLD);

    
  }//  end of the statement for (ijk=0;ijk<Ncomb;ijk++)
  
  xfclose (in);
  if (Myrank==0){
    fclose (casy);
  }
  
  /*
  fclose (out1);
  fclose (out2);
  fclose (out3);
  fclose (out4);
  fclose (out5);
  */
}


