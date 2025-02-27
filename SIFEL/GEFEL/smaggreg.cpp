#include "smaggreg.h"

smaggreg::smaggreg ()
{
  //  number of unknowns (in the original problem)
  nu=0;
  //  number of aggregates
  na=0;
  //  dimension of the kernel of aggregates
  dimker=0;
  
  //  numbers of nodes in aggregates
  //nnagr=NULL;
  //  list of node numbers in aggregates
  //lnnagr=NULL;
  
  //  numbers of unknowns in aggregates
  nuagr=NULL;
  //  list of unknowns id in aggregates
  uagr=NULL;

  //  localization matrices
  lm=NULL;
}

smaggreg::~smaggreg ()
{
  
}

/**
   function creates aggregates
   
   function assmebles list of nodes for each aggregates (lnnagr)
   function computes the numbers of nodes in each aggregate
   
   JK, 23.12.2006
*/
/*
void smaggreg::aggregates ()
{
  long i;
  
  if (nnagr!=NULL){
    delete [] nnagr;
  }

  nnagr=new long [na];
  
  //dodelat
  //nnagr[i]=?
  
  for (i=0;i<na;i++){
    if (lnnagr[i]!=NULL)
      delete [] lnnagr[i];
    lnnagr[i] = new long [nnagr[i]];
  }
  
  //  dodelat
  //lnnagr[][]=?
}
*/

/**
   function assembles tentative prolongator
   
   in Brezina thesis, page 68, tentative prolongator is assemled
   from vectors located in columns in the %matrix \hat{P}
   GEFEL contains class locmatrix which is based on rows
   Therefore, transposed prolongator is constructed here and
   algorithm is modified with respect to this fact

   JK, 23.12.2006
*/
 /*
void smaggreg::tentative_prolongator ()
{
  long i;
  
  //  number of aggregates
  //na=
  //  kernel dimension for each aggregate
  //dimker=
  //  number of unknowns
  //nu=

  tp.initiate_var (na*dimker,nu);

}
 */

/**
   function assembles localization matrices
   
   JK, 23.12.2006
*/
void smaggreg::localization_matrices (gtopology &gt)
{
  long i,j,k,l,ndofn,nid,*aux;
  
  //  allocation of localization matrices
  if (lm!=NULL){
    delete [] lm;
  }
  lm = new locmatrix [na];
  
  
  //  function assembles list of number of nodes in each aggregate
  //  funkce sestavuje pole poctu uzlu na jednotlivych agregatech
  aggreg.assemble_lnnagr ();
  //  function assembles list of node numbers (node id) in each aggregate
  //  funkce sestavuje cisla uzlu na jednotlivych agregatech
  aggreg.assemble_lnagr ();
  
  // ****************************************************
  //  determination of number of unknowns on aggregates
  //  vypocet poctu neznamych na agregatech
  // ****************************************************
  if (nuagr!=NULL)
    delete [] nuagr;
  nuagr = new long [na];
  
  //  loop over aggregates
  for (i=0;i<na;i++){
    nuagr[i]=0;
    //  loop over all nodes in the i-th aggregate
    for (j=0;j<aggreg.lnnagr[i];j++){
      nid = aggreg.lnagr[i][j];
      ndofn = gt.give_ndofn (nid);
      for (k=0;k<ndofn;k++){
	l=gt.give_dof (nid,k);
	if (l>0)
	  nuagr[i]++;
      }
    }
  }
  
  // ***********************************
  //  list of unknown id on aggregates
  // ***********************************
  if (uagr!=NULL){
    for (i=0;i<na;i++){
      delete [] uagr[i];
    }
    delete [] uagr;
  }
  
  uagr = new long* [na];
  for (i=0;i<na;i++){
    uagr[i] = new long [nuagr[i]];
  }
  
  //  loop over aggregates
  for (i=0;i<na;i++){
    nuagr[i]=0;
    //  loop over all nodes in the i-th aggregate
    for (j=0;j<aggreg.lnnagr[i];j++){
      nid = aggreg.lnagr[i][j];
      ndofn = gt.give_ndofn (nid);
      for (k=0;k<ndofn;k++){
	l=gt.give_dof (nid,k);
	if (l>0){
	  uagr[i][nuagr[i]]=l-1;
	  nuagr[i]++;
	}
      }
    }
  }
  
  // **************************************
  //  assembling of localization matrices
  // **************************************
  //  loop over aggregates
  for (i=0;i<na;i++){
    //  function defines number of rows and columns of localization matrix
    //  nuagr[i] is number of rows, nu is number of columns
    lm[i].initiate_var (nuagr[i],nu);
  }
  
  //  auxiliary array
  aux = new long [nu];
  for (i=0;i<nu;i++){
    aux[i]=1;
  }
  for (i=0;i<na;i++){
    //  number of nonzero matrix entries in rows
    lm[i].initiate_nncr (aux);
    //  addresses of first entries in rows
    lm[i].addresses ();
    //  column indices
    lm[i].initiate_ci (uagr[i]);
  }
  delete [] aux;
  
}

/**
   function assembles local aggregated matrices
   
   @param sm - stiffness %matrix of the whole system
   
   JK, 23.12.2006
*/
void smaggreg::local_aggregated_matrices (gmatrix &sm)
{
  long i;
  
  gm = new gmatrix [na];
  
  for (i=0;i<na;i++){
    lm[i].lmxmxlmt01 (sm,gm[i]);
  }
  
}

/**
   function decomposes (factorizes) %matrix of one aggregate
   
   @param sm - %matrix of one aggregate
   JK, 6.1.2007
*/
void smaggreg::decompose_aggreg_matrix (gmatrix &sm)
{
  sm.decompose_matrix ();
}

/**
   function assembles smoothed prolongator
   
   JK, 6.1.2007
*/
void smaggreg::smoothed_prolong ()
{
  long i,j,k,*aux,*cind;
  double *val;
  
  //  
  sp.initiate_var (na*dimker,nu);
  
  aux = new long [na*dimker];
  aggreg.assemble_lnuagr_sp (aux);
  
  //  dodelat vyhledani cisel neznamych podle seznamu dotcenych uzlu
  
  
  k=0;
  for (i=0;i<na;i++){
    for (j=0;j<dimker;j++){
      cind=new long [aux[k]];
      val=new double [aux[k]];
      
      //  vygenerovat prislusny RBM a ulozit do rbmv
      //  ziskam val o delce dlouheho vektoru
      //aggreg.assemble_smoothed_prol (i,rbmv,cind,val);
      
      sp.initiate_ci (k,cind);

      sp.initiate_lm (k,val);

      delete [] cind;
      delete [] val;
      k++;
    }
  }
}
