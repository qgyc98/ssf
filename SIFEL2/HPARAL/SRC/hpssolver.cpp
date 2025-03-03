#include "hpssolver.h"
#include "pglobal.h"
#include "solverm.h"
#include "mpi.h"

void par_solve_stochastic_problem (stochdriver *stochd)
{
  long i,j,l,ncs,nrm,buffsize,stop;
  double *buff;
  MPI_Status stat;
  
  stop=0;

  //  computation and distribution of buffsize
  if (Myrank==0){
    if (stochd->nstochvar > stochd->nprunknowns){
      //  number of input parameters is greater than number of output parameters
      buffsize = stochd->nstochvar + 1;
    }
    else{
      //  number of output parameters is greater than number of input parameters
      buffsize = stochd->nprunknowns + 1;
    }
    for (i=1;i<Nproc;i++){
      MPI_Send (&buffsize,1,MPI_LONG,i,Myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (&buffsize,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  buff = new double [buffsize];
  
  //  evaluation of samples
  if (Myrank==0){
    //  number of computed sample
    ncs=0;
    //  number of received messages
    nrm=0;
    if (stochd->nsampl>Nproc){
      for (i=1;i<Nproc;i++){
	//  generation of new input data
	stochd->assemble_new_values (ncs);
	ncs++;
	
	//  copying of input data into buffer
	stochd->give_new_invalues (buff);
	
	buff[buffsize-1]=0.0;
	
	
	//  scatering of data to appropriate processors
	MPI_Send (buff,buffsize,MPI_DOUBLE,i,Myrank,MPI_COMM_WORLD);
      }
      do{
	//  receiving of output
	MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	l=stat.MPI_TAG;
	nrm++;
	
	fprintf (stdout,"\n myrank %ld   nrm %ld     stat %ld",Myrank,nrm,l);
	
	//  import of output data from buffer to stochdriver
	stochd->save_new_outvalues (buff);
	//  save results in object stochdriver
	//  ta nula neni dobre, provizorium
	stochd->save_results (nrm-1);
	
	//  generation of new input data
	stochd->assemble_new_values (ncs);
	ncs++;

	fprintf (stdout,"\n myrank %ld   ncs %ld",Myrank,ncs);

	//  copying of input data into buffer
	stochd->give_new_invalues (buff);

	buff[buffsize-1]=0.0;

	
	//  scatering of data to appropriate processors
	MPI_Send (buff,buffsize,MPI_DOUBLE,l,Myrank,MPI_COMM_WORLD);
	
	if (ncs==stochd->nsampl){
	  stop=1;
	}
	
      }while(stop==0);
      
      for (i=nrm;i<stochd->nsampl;i++){
	MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      }
      for (i=1;i<Nproc;i++){
	buff[buffsize-1]=1.0;
	MPI_Send (buff,buffsize,MPI_DOUBLE,i,Myrank,MPI_COMM_WORLD);
      }
    }
    else{
      fprintf (stderr,"\n\n number of samples is less than number of used processors,");
      fprintf (stderr,"\n use number of processors equal to number of samples,\n\n");
      abort ();
    }
  }
  else{
    //  slave processors
    do{
      MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      if (buff[buffsize-1]>0.5)
	stop=1;
      
      if (stop==0){
	//  import of new values to object stochdriver
	stochd->save_new_invalues (buff);
	//  update of stochastic values
	stochd->replace_values ();
	
	Lsrs->clean_lhs ();
	Mm->clean_ip ();
	Mt->clean_nodes ();

	//  deterministic computation
	solve_mefel_deterministic_problem ();
	
	//  extraction of output data
	stochd->extractor ();
	//  import of computed values from stochdriver to buffer
	stochd->give_new_outvalues (buff);
	
	MPI_Send (buff,buffsize,MPI_DOUBLE,0,Myrank,MPI_COMM_WORLD);
      }
    }while(stop==0);
    
  }

  if (Myrank==0){
    if (Mp->stochasticcalc==3){
      //for (i=0;i<stochd->nprunknowns;i++){
      //fprintf (Out,"\n\n fuzzy number %ld",i);
      //stochd->fn[i].print (Out);
      //}
      
      j=0;
      for (i=0;i<stochd->npnd;i++){
	//fprintf (Out,"\n\n fuzzy posun %ld",i);
	stochd->fn[j].print (Out);
	j++;
	stochd->fn[j].print (Out);
	j++;
	fprintf (Out,"\n");
      }
      
      /*
	long k;
	for (i=0;i<stochd->npev;i++){
	fprintf (Out,"\n\n");
	
	for (k=0;k<6;k++){
	fprintf (Out,"\n");
	stochd->fn[j].print (Out);
	j++;
	}
	}
      */
      
      
      fprintf (Out,"\n\n\n\n");
      for (i=0;i<stochd->npev*6;i+=3){
	fprintf (Out,"\n");
	stochd->fn[j+i].print (Out);
      }
      fprintf (Out,"\n\n\n\n");
      for (i=1;i<stochd->npev*6;i+=3){
	fprintf (Out,"\n");
	stochd->fn[j+i].print (Out);
      }
      fprintf (Out,"\n\n\n\n");
      for (i=2;i<stochd->npev*6;i+=3){
	fprintf (Out,"\n");
	stochd->fn[j+i].print (Out);
      }
      
    }
    
  }
  
}

