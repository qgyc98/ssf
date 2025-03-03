void par_solve_stochastic_problem (stochdriver *stochd)
{
  long i,j,buffsize;
  double *buff;
  MPI_Status stat;
  
  if (Myrank==0){
    if (St->nsampl<Nproc-1){
      fprintf (stderr,"\n\n the number of samples is less than the number of processors");
      fprintf (stderr,"\n the program finishes, reduce the number of processors (file %s, line %d)\n",__FILE__,__LINE__);
      buffsize=0;
    }
    else{
      buffsize=St->nstochvar+1;
    }
    
    for (i=1;i<Nproc;i++){
      MPI_Send (&buffsize,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
  }
  else{
    MPI_Recv (&buffsize,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  
  if (buffsize>0){
    buff = new double [buffsize];
    
    if (Myrank==0){
      sampleid=0;
      for (i=1;i<Nproc;i++){
	St->assemble_new_values (sampleid);
	sampleid++;
	St->give_new_values (buff);
	MPI_Send (buff,buffsize,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
      }
    }
    else{
      sampleid=0;
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      St->save_new_values (buff);
      St->replace_values ();
      
      Lsrs->clean_lhs ();
      Mm->clean_ip ();
      Mt->clean_nodes ();
      solve_mefel_deterministic_problem ();
      St->extractor (sampleid);
      sampleid++;
    }
    
    
    if (Myrank==0){
      do{
	MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	St->assemble_new_values (sampleid);
	sampleid++;
	St->give_new_values (buff);
	buff[buffsize]=1.0;
	MPI_Send (buff,buffsize,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
	
      }
      while(sampleid<St->nsampl);
      
    }
    else{
      do{
	MPI_Send (buff,buffsize,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
	MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	St->save_new_values (buff);
	St->replace_values ();
	
	Lsrs->clean_lhs ();
	Mm->clean_ip ();
	Mt->clean_nodes ();
	solve_mefel_deterministic_problem ();
	St->extractor (sampleid);
	sampleid++;

	stop=buff[buffsize];
      }
      while(stop>0);
    }
    
  }



    
      for (i=1;i<St->nsampl+1;i++){
	MPI_Send (buff,buffsize,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
	MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	
      }
    }
  }
  
  
  
  if (Myrank==0){
    ncs=0;
    if (stochd->nsampl>Nproc){
      for (i=1;i<Nproc;i++){
	stochd->input_data (i,buff);
	MPI_Send (buff,buffsize,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
      }
      do{
	MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	l=stat.MPI_TAG;
	stochd->output_data (ncs,buff);
	ncs++;
	if (ncs==stochd->nsampl){
	  for (i=1;i<Nproc;i++){
	    buff[buffsize-1]=1.0;
	    MPI_Send (buff,buffsize,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
	  }
	  break;
	}
	stochd->input_data (ncs,buff);
	MPI_Send (buff,buffsize,MPI_DOUBLE,l,myrank,MPI_COMM_WORLD);
      }while();
    }
    else{
      abort ();
    }
  }
  else{
    //  slave processors
    do{
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      if (buff[buffsize-1]>0.5)  break;
      
      stochd->changevalues (buff);
      solve_deterministic_problem ();
      stochd->extractor (buff);
      
      MPI_Send (buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
    }while();
    
  }
  
  
}

