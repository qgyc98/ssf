#include "mpi.h"
#include "hpssolver.h"
#include "genfile.h"
#include "pglobal.h"
  /*
void par_solve_stochastic_problem (stochdriver *stochd)
{

  long i,l,buffsize;
  double buff;
  MPI_Status stat;
  
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
      MPI_Send (buffsize,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (buffsize,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  
  
  //  evaluation of samples
  if (Myrank==0){
    //  number of computed sample
    ncs=0;
    if (stochd->nsampl>Nproc){
      for (i=1;i<Nproc;i++){
	//  generation of new input data
	stochd->assemble_new_values (ncs);
	ncs++;
	//  copying of input data into buffer
	stochd->give_new_values (buff);
	//  scatering of data to appropriate processors
	MPI_Send (buff,buffsize,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
      }
      do{
	//  receiving of output
	MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	l=stat.MPI_TAG;
	
	//  import of output data from buffer to stochdriver
	stochd->save_new_outvalues (buff);
	//  save results in object stochdriver
	//  ta nula neni dobre, provizorium
	stochd->save_results (0);
	
	//  generation of new input data
	stochd->assemble_new_values (ncs);
	ncs++;
	//  copying of input data into buffer
	stochd->give_new_values (buff);
	//  scatering of data to appropriate processors
	MPI_Send (buff,buffsize,MPI_DOUBLE,l,myrank,MPI_COMM_WORLD);

	if (ncs==stochd->nsampl){
	  for (i=1;i<Nproc;i++){
	    buff[buffsize-1]=1.0;
	    MPI_Send (buff,buffsize,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
	  }
	  break;
	}

      }while();
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
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      if (buff[buffsize-1]>0.5)  break;
      
      //  import of new values to object stochdriver
      stochd->save_new_values (buff);
      //  update of stochastic values
      stochd->replace_values ();
      
      //  deterministic computation
      solve_deterministic_problem ();
      
      //  extraction of output data
      stochd->extractor (0);
      //  import of computed values from stochdriver to buffer
      stochd->save_new_outvalues (buff);
      
      MPI_Send (buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
    }while();
    
  }
  

}

*/  
