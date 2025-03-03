#include "mpi.h"
#include "pprobdesc.h"
#include "pglobal.h"
#include "genfile.h"
#include <string.h>

pprobdesc::pprobdesc (void)
{

}

pprobdesc::~pprobdesc (void)
{
}

/**
   function shifts indices of nodes
   
   this function is used if results are observed in GiD
*/
void pprobdesc::shift_indices ()
{
  long i,j,k,*nnarray,*nearray,*buff;
  MPI_Status stat;

  buff = new long [2];
  
  buff[0]=Gtm->nn;
  buff[1]=Gtm->ne;
  
  if (Myrank==0){
    nnarray = new long [Nproc];
    nearray = new long [Nproc];
    memset (nnarray,0,sizeof(*nnarray)*Nproc);
    memset (nearray,0,sizeof(*nearray)*Nproc);

    //  master contribution
    j=Psolm->domproc[0];
    nnarray[j]=buff[0];
    nearray[j]=buff[1];
    
    //  slaves contributions
    for (i=1;i<Nproc;i++){
      MPI_Recv (buff,2,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=Psolm->domproc[stat.MPI_TAG];
      nnarray[j]=buff[0];
      nearray[j]=buff[1];
    }
    
    j=nnarray[0];
    nnarray[0]=1;
    for (i=1;i<Nproc;i++){
      k=nnarray[i];
      nnarray[i]=nnarray[i-1]+j;
      j=k;
    }
    j=nearray[0];
    nearray[0]=1;
    for (i=1;i<Nproc;i++){
      k=nearray[i];
      nearray[i]=nearray[i-1]+j;
      j=k;
    }
    
    for (i=1;i<Nproc;i++){
      buff[0]=nnarray[i];
      buff[1]=nearray[i];
      MPI_Send (buff,2,MPI_LONG,i,Myrank,MPI_COMM_WORLD);
    }
    
    buff[0]=nnarray[0];
    buff[1]=nearray[0];
    
    delete [] nnarray;  delete [] nearray;
  }
  else{
    MPI_Send (buff,2,MPI_LONG,0,Myrank,MPI_COMM_WORLD);
    MPI_Recv (buff,2,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  
  //  first node index on subdomain (in global ordering; used for GiD)
  fni = buff[0];
  //  first element index on subdomain (in global ordering; used for GiD)
  fei = buff[1];
  
  delete [] buff;
}

