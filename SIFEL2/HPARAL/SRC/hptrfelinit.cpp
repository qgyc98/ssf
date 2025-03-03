#include "hptrfelinit.h"
#include "hpglobal.h"
#include "xfile.h"
#include "stacktrace.h"

/**
   function reads all data from the input file
   
   @param argv - pointer to input file names
   
   JK, 22.9.2010
   TKr, actualized 12/07/2022
*/
void hptrfel_init (int argc,const char *argv[])
{
  long i,mpt;
  char name[1100];
  MPI_Status stat;
  const char *aux;
  
  //  names of input files
  strcpy (name,argv[1]);
  sprintf (&name[strlen (argv[1])],"%d.in",Myrank+1);
  aux = argv[1];
  argv[1] = name;

  //  reading of the classical TRFEL input file
  trfel_init (argc,(const char **)(argv));
  argv[1] = aux;
  
  
  if (Myrank==0){
    switch (Tp->tprob){
    case stationary_problem:{
      mpt=50;
      break;
    }
    case nonstationary_problem:{
      mpt=60;
      break;
    }
    case nonlinear_nonstationary_problem:{
      mpt=61;
      break;
    }
    default:{
      print_err("unknown type of transport problem is required ",__FILE__, __LINE__, __func__);
    }
    }

    for (i=1;i<Nproc;i++){
      MPI_Send (&mpt,1,MPI_LONG,i,Myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (&mpt,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  if (mpt==50){
    Mtprob=stationary_problem;
  }
  if (mpt==60){
    Mtprob=nonstationary_problem;
  }
  if (mpt==61){
    Mtprob=nonlinear_nonstationary_problem;
  }
}


/**
   function reads all data from the input file
   
   @param argv - pointer to input file names
   
   JK, 4.6.2013
   TKr, actualized 12/07/2022
*/
void hptrfel_init_tiles (int argc,const char *argv[])
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
  
  //  the numbers of internal nodes on tiles
  Nidof = new long [Ntiletypes];
  for (i=0;i<Ntiletypes;i++){
    xfscanf (in,"%ld",Nidof+i);
  }
  
  //  the numbers of internal DOFs on tiles
  Nbdof = new long [Ntiletypes];
  
  //  the numbers of boundary DOFs on tiles in tiling
  Nbdoftiling = new long [Ntiles];
  for (i=0;i<Ntiles;i++){
    xfscanf (in,"%ld",Nbdoftiling+i);
  }
  
  //  the code numbers of Schur complements on tiles
  long tndof=0;
  Cnbn = new long* [Ntiles];
  for (i=0;i<Ntiles;i++){
    Cnbn[i] = new long [Nbdoftiling[i]];
    for (j=0;j<Nbdoftiling[i];j++){
      xfscanf (in,"%ld",&Cnbn[i][j]);
      if (tndof<Cnbn[i][j])
        tndof = Cnbn[i][j];
    }
  }
  
  //  the number of tilings
  xfscanf (in,"%ld",&Ncomb);

  fprintf (stderr,"procesor %2d  pocet typu dlazdic %ld  pocet dlazdic v dlazdeni %ld   tndof %ld\n",Myrank,Ntiletypes,Ntiles,tndof);

  //  the list of tiles in tilings
  Tiling = new long* [Ncomb];
  for (i=0;i<Ncomb;i++){
    Tiling[i] = new long [Ntiles];
    for (j=0;j<Ntiles;j++){
      xfscanf (in,"%ld",&Tiling[i][j]);
    }
  }

  long *pdof;
  long *aux;
  aux = new long [tndof];
  for (i=0;i<tndof;i++){
    aux[i]=0;
  }
  
  //  the number of prescribed DOFs
  long npdof;
  xfscanf (in,"%ld",&npdof);
  pdof = new long [npdof];
  Pdofval = new double [npdof];
  
  for (i=0;i<npdof;i++){
    xfscanf (in,"%ld %le",pdof+i,Pdofval+i);
    //aux[pdof[i]-1]=0-i-1;
  }
  
  xfclose (in);
  /*
  j=1;
  for (i=0;i<tndof;i++){
    if (aux[i]==0){
      aux[i]=j;
      j++;
    }
  }
  
  for (j=0;j<Ntiles;j++){
    for (k=0;k<Nbdoftiling[j];k++){
      Cnbn[j][k]=aux[Cnbn[j][k]-1];
    }
  }
  */
  //fprintf (stdout,"\n\n\n NA PROSTREDNIM UZLU JE NEZNAMA CISLO %ld \n\n\n\n",aux[139]);
  
  delete [] aux;
  delete [] pdof;
  

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
  const char *bargv = argv[1];
    
  FILE *out;
  out = fopen ("kontrola","w");
    
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
    trfel_init (argc,(const char **)(argv));
    
    Nbdof[i]=Ndoft-Nidof[i];
    
    //  assembling of conductivity matrix of a tile
    conductivity_matrix (0);
    
    //  allocation of array for right hand side
    Srhs[i] = new double [Ndoft];
    //  allocation of array for left hand side
    Slhs[i] = new double [Ndoft];
    fprintf(stderr, "Na domene %d je %ld stupnu volnosti\n", Myrank+1, Ndoft);
    
    //  assembling of right hand side
    trfel_right_hand_side (0,Srhs[i],Ndoft);
    //Kmat->printmat (out);
    //  copy of the stiffness matrix to arrays
    //  the stiffness matrix will be deleted
    Kmat->copygm (SSmat[i]);
    
    //  destruction of all objects

    delete Kmat;
    Kmat=NULL;
    delete Gtt;
    delete Tt;
    delete Tm;
    delete Tc;
    delete Tb;
    //delete Outdt;
    delete Tp;
    Tp=NULL;

  }//  end of the loop over the number of tiles
  fclose (out);
  argv[1]=bargv;
  
  Mtprob=trfel_tiles;

}



/**
   this function is used in the general implementation of tilings
   
   the number of processors has to be equal to the number of types of tiles
   each processor computes the Schur complement of the assigned tile
   the master processor solves tilings
   
   function reads all data from the input file
   
   @param argv - pointer to input file names
   
   JK, 9. 9. 2014
   TKr, actualized 12/07/2022
*/
void hptrfel_init_tiles_parsolver (int argc,const char *argv[], XFILE* &in)
{
  long i,lcid;
  char name[1100];
  
  //  in these problems lcid must be always equal to 0
  lcid=0;
  
  // **********************************************************
  //  reading input data of tiles
  //  each tile is assigned to a single processor
  // **********************************************************
  
  //  backup of the name of input files
  const char *bargv = argv[1];
  
  FILE *out;
  out = fopen ("kontrola","w");
  
  //  names of input files
  //  name stored as a char in argv[1] is copied to the char variable denoted name
  strcpy (name,argv[1]);
  //  strlen returns the length of char stored in argv[1]
  //  adds text from format (%d.in) to the end of char stored in name
  sprintf (&name[strlen (bargv)],"%d.in",Myrank+1);
  //  the name is copied back to argv[1] because it is used in the function trfelinit
  argv[1] = name;
  
  //  reading of the classical TRFEL input file
  trfel_init (argc,(const char **)(argv));
  
  //  restoration of the original name in argv[1]
  argv[1]=bargv;
  
  
  
  // *********************************************************************
  //  file with details about tiles and tiling is read by all processors
  // *********************************************************************
  in = xfopen(argv[1],"r");
  
  
  //  the number of types of tiles
  xfscanf (in,"%ld",&Ntiletypes);
  
  //  the numbers of internal DOFs on tiles
  Nidof = new long [Ntiletypes];
  //  the numbers of boundary/interface DOFs on tiles
  Nbdof = new long [Ntiletypes];
  for (i=0;i<Ntiletypes;i++){
    xfscanf (in,"%ld %ld",Nidof+i,Nbdof+i);
  }
  
  //  the number of all tilings solved in this execution
  //  different numbers of tiles in tilings are allowed
  xfscanf (in,"%ld",&Ncomb);
  
  fclose (out);
  
  Mtprob=trfel_tiles;
  
}
