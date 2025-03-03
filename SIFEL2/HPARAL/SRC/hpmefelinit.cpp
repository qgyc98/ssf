#include "hpmefelinit.h"
#include "hpglobal.h"
#include "xfile.h"
#include "stacktrace.h"

/**
   function reads all data from the input file
   
   @param argv - pointer to input file names
   
   TKr, 26/07/2022
*/
void hpmefel_init (int argc,const char *argv[])
{
  long i,mp;
  char name[1100];
  MPI_Status stat;
  const char *aux;
  
  //  names of input files
  strcpy (name,argv[1]);
  sprintf (&name[strlen (argv[1])],"%d.in",Myrank+1);
  aux = argv[1];
  argv[1] = name;

  //  reading of the classical MEFEL input file
  mefel_init (argc,(const char **)(argv));
  argv[1] = aux;
  
  
  if (Myrank==0){
    switch (Mp->tprob){
    case linear_statics:{
      mp=1;
      break;
    }
    case mech_timedependent_prob:{
      mp=15;
      break;
    }
    default:{
      print_err("unknown type of mechanical problem is required",__FILE__,__LINE__,__func__);
    }
    }

    for (i=1;i<Nproc;i++){
      MPI_Send (&mp,1,MPI_LONG,i,Myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (&mp,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  if (mp==1){
    Mtprobm=linear_statics;
  }
  if (mp==15){
    Mtprobm=mech_timedependent_prob;
  }
}


/**
   function reads all data from the input file
   
   @param argv - pointer to input file names
   
   JK, 4.6.2013
   TKr, actualized 12/07/2022
*/
void hpmefel_init_tiles (int argc,const char *argv[])
{
  long i,j,lcid;
  char name[1100],namebackup[1100];
  XFILE *in;
  
  //  in these problems lcid must be always equal to 0
  lcid=0;
  
  // *****************************************************
  //  data about tiles, tilings and combinations is read
  // *****************************************************
  in = xfopen(argv[1],"r");
  
  

  //  the number of types of tiles
  xfscanf (in,"%ld",&Ntiletypes);
  //  the number of tiles in a tiling
  xfscanf (in,"%ld",&Ntiles);
  fprintf (stderr,"\n myrank %d  Ntiletypes %ld    Ntiles %ld",Myrank,Ntiletypes, Ntiles);

  //  the numbers of internal DOFs on tiles
  Nidof = new long [Ntiletypes];
  for (i=0;i<Ntiletypes;i++){
    xfscanf (in,"%ld",Nidof+i);
  }
  
  //  the numbers of internal DOFs on tiles
  Nbdof = new long [Ntiletypes];
  
  //  the numbers of boundary DOFs on tiles in tiling
  fprintf (stderr,"\n Myrank %d  Nbdoftiling ", Myrank);
  Nbdoftiling = new long [Ntiles];
  for (i=0;i<Ntiles;i++){
    xfscanf (in,"%ld",Nbdoftiling+i);
    fprintf (stderr," %ld",Nbdoftiling[i]);
  }
  
  //  the code numbers of Schur complements on tiles
  Cnbn = new long* [Ntiles];
  for (i=0;i<Ntiles;i++){
    Cnbn[i] = new long [Nbdoftiling[i]];
    for (j=0;j<Nbdoftiling[i];j++){
      xfscanf (in,"%ld",&Cnbn[i][j]);
    }
    fprintf (stderr,"\n Myrank %d  i %ld j %ld",Myrank,i,j);
  }

  FILE *out;
  out = fopen ("jouda","w");
  for (i=0;i<Ntiles;i++){
    fprintf (out,"\n Myrank %d\n",Myrank);
    for (j=0;j<Nbdoftiling[i];j++){
      fprintf (out," %4ld",Cnbn[i][j]);
    }
  }
  fprintf (out,"\n");
  fclose (out);
  
  //  the number of tilings
  xfscanf (in,"%ld",&Ncomb);
  fprintf (stderr,"\n myrank %d  Ncomb %ld",Myrank,Ncomb);

  //  the list of tiles in tilings
  Tiling = new long* [Ncomb];
  for (i=0;i<Ncomb;i++){
    Tiling[i] = new long [Ntiles];
    for (j=0;j<Ntiles;j++){
      xfscanf (in,"%ld",&Tiling[i][j]);
    }
  }
  xfclose (in);



  // **********************************************************
  //  reading input data of tiles
  //  all tiles will be solved on each processor in order to
  //  save time because no communication will be needed
  // **********************************************************
  
  //  allocation of instances of gmatrix
  //  stiffness matrices will be stored there
  SSmat = new gmatrix [Ntiletypes];
  //  allocation of arrays for right hand sides
  Srhs = new double* [Ntiletypes];
  //  allocation of arrays for right hand sides
  Slhs = new double* [Ntiletypes];


  
  //  backup of the name of input files
  strcpy (namebackup,argv[1]);
  //namebackup=argv[1];
 
  //  loop over the number of tiles
  for (i=0;i<Ntiletypes;i++){
    
    //  names of input files
    //  name stored as a char in argv[1] is copied to the char variable denoted name
    strcpy (name,namebackup);
    //  strlen returns the length of char stored in argv[1]
    //  adds text from format (%d.in) to the end of char stored in name
    //sprintf (&name[strlen (argv[1])],"%d.in",Myrank+1);
    sprintf (&name[strlen (namebackup)],"%ld.in",i+1);
    //  the name is copied back to argv[1] because it is used in the function mefelinit
    argv[1] = name;
    
    //  reading of the classical MEFEL input file
    mefel_init (argc,(const char **)(argv));
    
    Nbdof[i]=Ndofm-Nidof[i];
    
    //  assembling of stiffness matrix of a tile
    stiffness_matrix (lcid);
    
    for (j=0;j<4;j++){
      fprintf (stdout,"%ld nodes %ld %ld\n",i,Gtm->gnodes[j].cn[0],Gtm->gnodes[j].cn[1]);
    }
    for (j=0;j<Smat->sky->n+1;j++){
      fprintf (stdout," %ld adr %ld \n",i,Smat->sky->adr[j]);
    }
    for (j=0;j<Smat->sky->n;j++){
      for (long k=Smat->sky->adr[j];k<Smat->sky->adr[j+1];k++){
        fprintf (stdout," %ld Smat %ld %ld  %le\n",i,j,k,Smat->sky->a[k]);
      }
    }

    //  allocation of array for right hand side
    Srhs[i] = new double [Ndofm];
    //  allocation of array for left hand side
    Slhs[i] = new double [Ndofm];
    fprintf(stderr, "Na domene %d je %ld stupnu volnosti\n", Myrank+1, Ndofm);
    
    //  assembling of right hand side
    mefel_right_hand_side (lcid,Srhs[i]);
    
    for (j=0;j<Ndofm;j++){
      fprintf (stdout," %ld rhs  %ld   %le\n",i,j,Srhs[i][j]);
    }

    //  copy of the stiffness matrix to arrays
    //  the stiffness matrix will be deleted
    Smat->copygm (SSmat[i]);
    
    switch (Mp->tprob){
    case linear_statics:{
      Mtprobm=linear_statics;
      break;
    }
    case mat_nonlinear_statics:{
      Mtprobm=mat_nonlinear_statics;
      break;
    }
    case mech_timedependent_prob:{
      Mtprobm=mech_timedependent_prob;
      break;
    }
    default:{
      
    }
    }
    
    if (i!=4){
      //  destruction of all objects
      delete Smat;
      Smat=NULL;
      delete Mp;
      Mp=NULL;
      delete Gtm;
      delete Mt;
      delete Mm;
      delete Mc;
      delete Mb;
      //delete Outdm;
    }
    else{
      delete Smat;
      Smat=NULL;
      delete Mp;
      Mp=NULL;
      delete Gtm;
      delete Mc;
      delete Mb;
      
      Mmm = Mm;
      Mtt = Mt;
    }
    
  }//  end of the loop over the number of tiles
  
  argv[1]=namebackup;
  


}
