#include "fuzzygen.h"
#include <math.h>


fuzzygen::fuzzygen (void)
{
  //  number of fuzzy variables
  nfv=0;
  //  number of alpha-cuts
  nalph=0;
  //
  nval=0;
  
  //  number of combinations in one alpha-cut
  ncomb=0;
  //  total number of combinations
  tncomb=0;
  
  //  number of actual alpha-cut
  actalph=0;
  //  number of actual combination in actual alpha-cut
  actcomb=0;
  
  
  alpha=NULL;
  aux=NULL;
  fn = NULL;
}



fuzzygen::~fuzzygen ( void )
{
  delete [] alpha;
  delete [] aux;
  delete [] fn;
}

/**
   function reads input data
   
   @param in - input file
   
   JK, 18.8.2005
*/
void fuzzygen::read (XFILE *in)
{
  long i;

  //  number of fuzzy variables
  xfscanf (in,"%ld",&nfv);
  //  number of alpha-cuts
  xfscanf (in,"%ld",&nalph);
  
  //  number of values in one array
  nval=nalph*2-1;
  
  if (nfv<1){
    fprintf (stderr,"\n\n number of fuzzy variables is less than 1 in function fuzzygen::read (file %s, line %d)\n",
	     __FILE__,__LINE__);
  }
  if (nalph<0){
    fprintf (stderr,"\n\n number of alpha-cuts is negative in function fuzzygen::read (file %s, line %d)\n",
	     __FILE__,__LINE__);
  }

  if (alpha==NULL){
    alpha = new double [nalph];
  }
  if (aux==NULL){
    aux = new long [nfv];
  }
  
  for (i=0;i<nalph;i++){
    xfscanf (in,"%lf",alpha+i);
  }
  
  //  allocation of fuzzy variables
  fn = new fuzzynum [nfv];
  
  for (i=0;i<nfv;i++){
    fn[i].read (in);
    //fn[i].onearray ();
  }
  
  //  total number of combinations
  totcombnumber ();
}

/**
   function computes total number of combinations for given number of variables
   and number of alpha-cuts
   
   JK, 18.8.2005
*/
void fuzzygen::totcombnumber ()
{
  //  number of combinations for one alpha-cut
  combnumber ();
  //  total number of combinations
  tncomb=nalph*ncomb;
  
  fprintf (stdout,"\n total number of combinations  %ld",tncomb);
}

/**
   function computes number of combinations for given number of variables
   and for given alpha-cut
   
   JK, 18.8.2005, revised 17.4.2007
*/
void fuzzygen::combnumber ()
{
  long i;
  
  ncomb=1;
  for (i=0;i<nfv;i++){
    ncomb *= 2;
  }
  
}

/**
   function generates combinations based on alpha-cuts
   
   there are minimum and maximum values for given alpha
   aux[i]=0 means minimum value on the i-th variable
   aux[i]=1 means maximum value on the i-th variable
   the latest component of the array aux is changing fastest
   
   @param alphid - alpha-cut id
   @param avi - array of input values
   
   JK, 8.10.2005
*/
void fuzzygen::gener_alphacuts (long alphid,double *avi)
{
  long i,change;
  
  //  actual values assembling
  for (i=0;i<nfv;i++){
    if (aux[i]==0){
      avi[i]=fn[i].give_min (alpha[alphid]);
    }
    else{
      avi[i]=fn[i].give_max (alpha[alphid]);
    }
  }
  
  
  if (actcomb!=ncomb-1){
    //  update of auxiliary array
    change=0;
    aux[nfv-1]++;
    
    if (aux[nfv-1]>1){
      change=1;
      aux[nfv-1]=0;
    }
    
    for (i=nfv-2;i>-1;i--){
      if (change==1){
	aux[i]++;
	if (aux[i]>1){
	  aux[i]=0;
	  change=1;
	}
	else{
	  change=0;
	}
      }
    }
  }
}

/**
   function generates combinations based on alpha cuts
   auxiliary array aux is used, aux contains nfv (number
   of fuzzy variables) components
   
   JK, 9.10.2005
*/
void fuzzygen::gener_allcomb (double *avi)
{
  long i,j;
  
  i=nfv-1;
  if (aux[i]<nval){
    avi[i]=fn[i].give_val (aux[i]);
    aux[i]++;
    j=0;
  }
  else{
    aux[i]=0;
    avi[i]=fn[i].give_val (aux[i]);
    aux[i]++;
    j=1;
  }

  for (i=nfv-2;i>-1;i--){
    if (j==1)
      aux[i]++;
    if (aux[i]<nval){
      avi[i]=fn[i].give_val (aux[i]);
      j=0;
    }
    else{
      aux[i]=0;
      avi[i]=fn[i].give_val (aux[i]);
      j=1;
    }
  }
  
}


/**
   function defines new input values in fuzzy computations
   
   @param sampleid - sample id
   @param avi - array of input values
   @param out - output stream
   
   JK, 17.4.2007
*/
void fuzzygen::give_new_values (long sampleid,double *avi,FILE *out)
{
  long i;
  
  if (actalph==0 && actcomb==0)
    fprintf (stdout,"\n\n alpha-cut number %ld is solved (alpha is %f)\n",actalph,alpha[actalph]);
  
  if (actcomb==ncomb){
    //  new alpha-cut is prepared
    actalph++;
    if (actalph>nalph){
      fprintf (stderr,"\n\n number of actual alpha-cut is greater than the number of all alpha-cuts");
      fprintf (stderr,"\n in function give_new_values (file %s, line %d)\n",__FILE__,__LINE__);
    }
    actcomb=0;
    
    fprintf (stdout,"\n\n alpha-cut number %ld is solved (alpha is %f)\n",actalph,alpha[actalph]);
  }
  
  if (actcomb==0){
    for (i=0;i<nfv;i++){
      aux[i]=0;
    }
  }
  
  //  generation of data in actual combination
  gener_alphacuts (actalph,avi);
  

  fprintf (out,"\n sampleid %ld,  actalph %ld,  actcomb %ld   ",sampleid,actalph,actcomb);
  for (i=0;i<nfv;i++){
    fprintf (out," %e",avi[i]);
  }


  actcomb++;

}

/**
   function saves results to fuzzy numbers
   
   @param fnum - array of fuzzy numbers for output
   @param nprunknowns - number of printed output variables
   @param avo - array of output values
   
   JK, 17.4.2007
*/
void fuzzygen::save_values (fuzzynum *fnum,long nprunknowns,double *avo)
{
  long i;
  double min,max,al;
  
  al=alpha[actalph];

  for (i=0;i<nprunknowns;i++){
    min=fnum[i].give_min (al);
    max=fnum[i].give_max (al);
    if (avo[i]<min)
      fnum[i].save_alp_min (actalph,al,avo[i]);
    if (avo[i]>max)
      fnum[i].save_alp_max (actalph,al,avo[i]);
  }
}
