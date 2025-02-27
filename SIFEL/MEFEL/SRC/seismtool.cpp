#include "seismtool.h"
#include "global.h"
#include "mechtop.h"
#include "node.h"

seismtool::seismtool (void)
{
  nsac=0;
  
  direction = NULL;
  seism = NULL;
  gf = NULL;
  name = NULL;
}

seismtool::~seismtool (void)
{
  delete [] direction;
  delete [] seism;
  delete [] gf;
  delete [] name;
}

/**
   function reads seismic load from the opened text file
   
   @param in - pointer to input file
   
   JK
*/
void seismtool::read (XFILE *in)
{
  long i,j;
  
  //  number of seismic acceleration components
  xfscanf (in,"%ld",&nsac);
  //  directions of seismic loads
  direction = new dirdynload [nsac];
  
  //  amplitude functions or response spectra
  gf = new gfunct [nsac];
  for (i=0;i<nsac;i++){
    xfscanf (in,"%m", &dirdynload_kwdset, (int*)&direction[i]);
  }
  
  xfscanf (in,"%ld",&j);
  if (j==0){
    for (i=0;i<nsac;i++){
      gf[i].read (in);
    }
  }
  if (j==1){
    name = new char[1001];
    XFILE *accelin;
    xfscanf(in, " %1000a", name);

    accelin = xfopen (name,"r");
    accelin->warning = 1;
    accelin->kwdmode = ignore_kwd;
    //in->kwdmode = sequent_mode;
    accelin->ignorecase = 1;
    
    for (i=0;i<nsac;i++){
      gf[i].read (accelin);
    }
    xfclose (accelin);
  }
  
  //  array containing components of right hand side
  seism = new double [Ndofm*nsac];
  nullv (seism,Ndofm*nsac);
  //  initiation of auxiliary vectors used in seismic loading
  seisminit (seism);
}



/**
   function prints seismic load to the opened text file
   
   @param in - pointer to input file
   
   JK
*/
void seismtool::print (FILE *out)
{
  long i;
  
  //  number of seismic acceleration components
  fprintf (out,"%ld\n", nsac);
  
  //  amplitude functions or response spectra
  for (i=0; i<nsac; i++){
    fprintf (out,"%d ", int(direction[i]));
  }
  fprintf(out, "\n");

  if (name == NULL){
    fprintf(out, "0\n");
    for (i=0; i<nsac; i++){
      gf[i].print (out);
    }
  }
  else
    fprintf(out, "1\n%1000s\n", name);    
}



/**
   function assembles auxiliary vector used in the right hand side
   in seismic loading
   
   @param seis - auxiliary vector
   
   1.2.2005, JK
*/
void seismtool::seisminit (double *seism)
{
  long i,j,k,n,ndofn,m,*cn;
  
  n=Ndofm;
  //  initiation of meaning of DOFs
  Mt->alloc_meaning ();
  //  definition of meaning of DOFs
  Mt->define_meaning ();
  
  //  loop over nodes
  for (i=0;i<Mt->nn;i++){
    //  number of DOFs on node
    ndofn = Mt->give_ndofn (i);
    cn = new long [ndofn];
    //  code numbers of node DOFs
    Mt->give_node_code_numbers (i,cn);
    
    //  loop over node DOFs
    for (j=0;j<ndofn;j++){
      //  meaning of DOF
      m = Mt->nodes[i].meaning[j];
      if (cn[j]>0){
	for (k=0;k<nsac;k++){
	  if (m == direction[k])
	    seism[n*k+cn[j]-1]=1.0;
	}
      }
    }
    
    delete [] cn;
  }
  
}

/**
   function computes contributions to right hand side from seismic load
   
   @param rhs - pointer to right hand side
   @param time - actual time
   
   JK, 21.8.2005
*/
void seismtool::assemble (double *rhs,double time)
{
  long i,j,n;
  double a,*aux;
  
  n=Ndofm;
  //  auxiliary array
  aux = new double [n];
  
  for (i=0;i<nsac;i++){
    nullv (aux,n);
    
    copyv (seism+i*n,aux,n);
    
    //  scale factor
    a=gf[i].getval(time);
    //a=1.0;

    //  updating of rhs array
    for (j=0;j<n;j++){
      rhs[j]-=a*aux[j];
    }
  }
  
  delete [] aux;
}

