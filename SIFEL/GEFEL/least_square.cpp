#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "least_square.h"
#include "gadaptivity.h"
#include "vector.h"
#include "matrix.h"
#include "mathem.h"
#include "gtopology.h"
#include "basefun.h"

/*
  // max number of nodes on element is 8
  
  do compute_patches dodelat volbu ncoef podle poctu sp. pointu
  zjistit zda nekde nejde pouzit mathem:shaker
  
  vydeletovat atributy

  konstruktor trochu vykosit, je tam nejakej bordel

*/


least_square::least_square (long dim,long nvals,gtopology *gt,long flag)
{
  least_square::gt = gt;
  least_square::dim = dim;
  least_square::flag = flag;
  least_square::nvals = nvals;
  nn = gt->nn;
  ne = gt->ne;
  
  nsp = NULL;
  spcoord = NULL;
  spvalue = NULL;
  nodvalue = NULL;
  
  insidelem = new long [ne];
  insidenod = new long* [nn];
  maxcoord = new double* [ne];
  
  nadjel = 0;
  adjel = NULL;
  
  long i;
  for (i=0;i<nn;i++){
    insidenod[i] = new long [1];
    insidenod[i][0] = 1;
  }
  for (i=0;i<ne;i++)
    maxcoord[i] = new double[2*dim];
}

least_square::~least_square (void)
{
  if (insidenod != NULL)  for (long i=0;i<nn;i++)  delete [] insidenod[i];
  if (maxcoord  != NULL)  for (long i=0;i<ne;i++)  delete [] maxcoord[i];
  if (spcoord   != NULL)  for (long i=0;i<ne;i++)  delete [] spcoord[i];
  
  delete [] nsp;
  delete [] insidenod;
  delete [] maxcoord;
  delete [] spcoord;
  //  delete [] nodvalue;
}


/**
   SPR Superconvergent patch recovery
   This is main function of class least_square.
   It smooths discontinuous 'values' which are known in 'sampling points'.
   Answer is array of smoothed values in nodes.
   Sampling points ARE IDENTICAL to main integration points in first block on elements.
   
   @param spvalue - 2D array of rough values in sampling points, dimension is (gt->ne) x (nip*tncomp)
   @param nodvalue - answer, 1D array, size is (least_square::nn * least_square::nvals)
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void least_square::spr_default (FILE *out,double **spvalue,double *nodvalue)
{
  fprintf (stdout,"\n VVV===========  SPR_smoothing  ===========VVV");
  
  if (gt->nadjelnod == NULL)
    gt->adjacelem (out);
  
  least_square::spvalue  = spvalue;
  least_square::nodvalue = nodvalue;
  
  nsp = new long [ne];
  spcoord = new double* [ne];
  
  nsp_spcoord_maxcoord_assembling ('a');
  insidenod_assembling ();
  compute_patches_spr (0);
  
  least_square::spvalue = NULL;
  least_square::nodvalue = NULL;
  
  fprintf (stdout,"\n AAA===========  SPR_smoothing  ===========AAA");
}


/**
   SPR Superconvergent patch recovery
   This is main function of class least_square.
   It smooths discontinuous 'values' which are known in 'sampling points'.
   Answer is array of smoothed values in nodes.
   Sampling points DO NOT HAVE TO BE IDENTICAL to main integration points in first block on elements,
   they are defined by arrays nsp,spcoord and spvalue.
   
   @param nsp - number of sampling points on element, dimension is least_square::ne
   @param spcoord - cartesian coorinates of sampling points, dimension is (least_square::ne x nsp[i]*least_square::dim)
   @param spvalue - 'values' in sampling points, dimension is (least_square::ne x nsp[i]*least_square::nvals)
   @param nodvalue - answer, dimension id is (least_square::nn * least_square::nvals)
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void least_square::spr_optional (FILE *out,long *nsp,double **spcoord,double **spvalue,double *nodvalue,long cut_flag)
{
  fprintf (stdout,"\n VVV===========  SPR_smoothing  ===========VVV");
  
  if (gt->nadjelnod == NULL)
    gt->adjacelem (out);
  
  least_square::nsp      = nsp;
  least_square::spcoord  = spcoord;
  least_square::spvalue  = spvalue;
  least_square::nodvalue = nodvalue;
  
  nsp_spcoord_maxcoord_assembling ('m');
  insidenod_assembling ();
  compute_patches_spr (cut_flag);
  
  least_square::nsp = NULL;
  least_square::spcoord = NULL;
  least_square::spvalue = NULL;
  least_square::nodvalue = NULL;
  
  fprintf (stdout,"\n AAA===========  SPR_smoothing  ===========AAA");
}


/**
   ONLY FOR LINEAR TRIANGLE (312) !!!!!!!!!!
   
   This is main function of class least_square.
   It aproximates values from nodes to 'sampling points'.
   Answer is array of aproximated values in 'sampling points'.
   Sampling points are defined by arrays nsp and spcoord.
   
   @param nsp - number of sampling points on element, dimension is least_square::ne
   @param spcoord - cartesian coorinates of sampling points, dimension is (least_square::ne x nsp[i]*least_square::dim)
   @param nodvalue - known values in nodes, dimension is (least_square::nn * least_square::nvals)
   @param spvalue - answer - smoothed 'values' in sampling points, dimension is (least_square::ne x nsp[i]*least_square::nvals)
   
   created  3.2.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void least_square::L2_nod2sp (FILE *out,long *nsp,double **spcoord,double *nodvalue,double **spvalue,char typ_patch)
{
  fprintf (stdout,"\n VVV===========  L2_nod2sp  ===========VVV");
  
  if (gt->nadjelnod == NULL)
    gt->adjacelem (out);
  
  least_square::nsp = nsp;
  least_square::spcoord = spcoord;
  least_square::spvalue = spvalue;
  least_square::nodvalue = nodvalue;
  
  nsp_spcoord_maxcoord_assembling ('m');
  insidenod_assembling ();
  compute_patches_nod2sp (typ_patch);
  
  least_square::nsp = NULL;
  least_square::spcoord = NULL;
  least_square::spvalue = NULL;
  least_square::nodvalue = NULL;
  
  fprintf (stdout,"\n AAA===========  L2_nod2sp  ===========AAA");
}


/**
   ONLY FOR LINEAR TRIANGLE (312) !!!!!!!!!!
   
   This is main function of class least_square.
   It aproximates values from known 'sampling points' to another 'sampling points' = 'answer points'.
   Answer is array of aproximated values in another 'sampling points' = 'answer points'.
   ????????
   Sampling and answer points ARE IDENTICAL to main integration points on elements,
   ????????
   they are defined by arrays nsp, spcoord ... .
   
   @param nsp - number of sampling points on element, dimension is least_square::ne
   @param spcoord - cartesian coorinates of sampling points, dimension is (least_square::ne x nsp[i]*least_square::dim)
   @param spvalue - answer - smoothed 'values' in sampling points, dimension is (least_square::ne x nsp[i]*least_square::nvals)
   @param nodvalue - known values in nodes, dimension is (least_square::nn * least_square::nvals)
   
   created  4.5.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void least_square::L2_sp2sp (FILE *out,double **spvalue,long nap,const long *parentel,const double **apcoord,double **apvalue)
{

  // dodelat reakci na flag
  
  fprintf (stdout,"\n ==============  SPR_2_smoothing  ==============");
  
  if (gt->nadjelnod == NULL)
    gt->adjacelem (out);
  
  least_square::spvalue = spvalue;
  
  long *nadjap,**adjap;
  
  nadjap = new long [nn];
  adjap = new long* [nn];
  
  nsp = new long [ne];
  spcoord = new double* [ne];
  
  
  nsp_spcoord_maxcoord_assembling ('a');
  insidenod_assembling ();
  adjap_assembling (nap,parentel,apcoord,nadjap,adjap);
  compute_patches_sp2sp (nap,apcoord,nadjap,(const long**)adjap,apvalue);
  
  
  delete [] nadjap;
  for (long i=0;i<nn;i++) delete [] adjap[i];
  delete [] adjap;
  
  least_square::spvalue = NULL;
}


/**
   Function fills up arrays:
    'nsp' (number of sampling points) by default values,
    'spcoord' by natural coordinates of default sampling points,
    'maxcoord' by extreme element coordinates.
   
   @param nsma_flag - determines which array will be actually filled: 'a' = all arrays, 'm' = only array maxcoord
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void least_square::nsp_spcoord_maxcoord_assembling (char nsma_flag)
{
  long i,j,k,nne;
  double xi,eta;
  vector x,y,z,bf;
  
  
  for (i=0;i<ne;i++){
    nne = gt->gelements[i].nne;
    if (dim==2) { reallocv (nne,x); reallocv (nne,y);                  gt->give_node_coord2d (x,y,i);   }
    if (dim==3) { reallocv (nne,x); reallocv (nne,y); reallocv (nne,z);  gt->give_node_coord3d (x,y,z,i); }
    
    if (nsma_flag == 'a'){
      
      switch (gt->gelements[i].auxinf){
      case 312:{
	nsp[i] = 1;
	spcoord[i] = new double[2];
	
	spcoord[i][0] = (x[0] + x[1] + x[2])/3.0;
	spcoord[i][1] = (y[0] + y[1] + y[2])/3.0;
	break;
      }
      case 412:{
	nsp[i] = 1;
	spcoord[i] = new double[2];
	
	spcoord[i][0] = (x[0] + x[1] + x[2]+ x[3])/4.0;
	spcoord[i][1] = (y[0] + y[1] + y[2]+ y[3])/4.0;
	break;
      }
      case 622:{
	nsp[i] = 3;
	spcoord[i] = new double[6];
	
	reallocv (nne,bf);
	xi = eta = 1.0/6.0;
	for (j=0;j<3;j++){
	  if (j==1)  xi=4.0/6.0;
	  if (j==2) {xi=1.0/6.0; eta=4.0/6.0;}
	  
	  bf_quad_3_2d (bf.a,xi,eta);
	  scprd (bf,x,spcoord[i][2*j  ]);
	  scprd (bf,y,spcoord[i][2*j+1]);
	}
	break;
      }
      case 822:{
	nsp[i] = 4;
	spcoord[i] = new double[8];
	
	reallocv (nne,bf);
	xi = eta = 0.577350269189626;
	for (j=0;j<2;j++){
	  xi = -xi;
	  for (k=0;k<2;k++){
	    eta = -eta;
	    
	    bf_quad_4_2d (bf.a,xi,eta);
	    scprd (bf,x,spcoord[i][4*j+2*k  ]);
	    scprd (bf,y,spcoord[i][4*j+2*k+1]);
	  }
	}
	break;
      }
      case 413:{
	nsp[i] = 1;
	spcoord[i] = new double[3];
	
	spcoord[i][0] = (x[0] + x[1] + x[2] + x[3])/4.0;
	spcoord[i][1] = (y[0] + y[1] + y[2] + y[3])/4.0;
	spcoord[i][2] = (z[0] + z[1] + z[2] + z[3])/4.0;
	
	break;
      }
      default:{
	fprintf (stderr,"\n\n wrong dimdegnne in function nsp_spcoord_maxcoord_assembling (%s, line %d)",__FILE__,__LINE__);
	break;
      }}
    }
    
    if (nsma_flag == 'a' || nsma_flag == 'm'){
      
      switch (gt->gelements[i].auxinf){
      case 312:
      case 622:{
	maxmin_3 (x.a,maxcoord[i][0],maxcoord[i][1]);
	maxmin_3 (y.a,maxcoord[i][2],maxcoord[i][3]);
	break;
      }
      case 412:
      case 822:{
	maxmin_4 (x.a,maxcoord[i][0],maxcoord[i][1]);
	maxmin_4 (y.a,maxcoord[i][2],maxcoord[i][3]);
	break;
      }
      case 413:{
	maxmin_4 (x.a,maxcoord[i][0],maxcoord[i][1]);
	maxmin_4 (y.a,maxcoord[i][2],maxcoord[i][3]);
	maxmin_4 (z.a,maxcoord[i][4],maxcoord[i][5]);
	break;
      }
      default:{
	fprintf (stderr,"\n\n wrong dimdegnne in function nsp_spcoord_maxcoord_assembling (%s, line %d)",__FILE__,__LINE__);
	break;
      }}
    }
  }
}


/**
   loop 1 :
   assembles vector of node position in domain
   2D: insidenod[i]== 1 -> inside node     of domain
       insidenod[i]==-1 -> node on edge    of domain
       insidenod[i]==-2 -> node at vertex  of domain
   3D: insidenod[i]== 1 -> inside node     of domain
       insidenod[i]==-1 -> node on boundary area  of domain
       insidenod[i]==-2 -> node on edge    of domain
       insidenod[i]==-3 -> node at vertex  of domain
       
   loop 2 :
   for every "node at vertex of domain"=vn is found "oposite node"=on  =>
   =>  inside[on][inside[on][0]++ - 1] = nv
   where 'on' is "oposite node" already "inside[on][0]" times.
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void least_square::insidenod_assembling (void)
{
  long i,j,k,l,nne,*n,in[8];            // max number of nodes on element is 8
  long jnod1,jnod2,sepnod,on,sum;
  long jnod,joinel,sepnodes[3];
  long aux3[3],n4,aux4[4];
  
  
  // loop 1
  for (i=0;i<ne;i++){
    nne = gt->give_nne (i);
    n = gt->gelements[i].nodes;
    for (j=0;j<nne;j++) in[j] = 0;
    
    switch (gt->gelements[i].auxinf){
    case 312:{
      if (edge_position (gt,n[0],n[1]) == 1) {in[0]--; in[1]--;}
      if (edge_position (gt,n[1],n[2]) == 1) {in[1]--; in[2]--;}
      if (edge_position (gt,n[2],n[0]) == 1) {in[2]--; in[0]--;}
      break;
    }
    case 412:{
      if (edge_position (gt,n[0],n[1]) == 1) {in[0]--; in[1]--;}
      if (edge_position (gt,n[1],n[2]) == 1) {in[1]--; in[2]--;}
      if (edge_position (gt,n[2],n[3]) == 1) {in[2]--; in[3]--;}
      if (edge_position (gt,n[3],n[0]) == 1) {in[3]--; in[0]--;}
      break;
    }
    case 622:{
      if (gt->nadjelnod[n[3]] == 1) {in[0]--; in[1]--; in[3]--;}
      if (gt->nadjelnod[n[4]] == 1) {in[1]--; in[2]--; in[4]--;}
      if (gt->nadjelnod[n[5]] == 1) {in[2]--; in[0]--; in[5]--;}
      break;
    }
    case 822:{
      if (gt->nadjelnod[n[4]] == 1) {in[0]--; in[1]--; in[4]--;}
      if (gt->nadjelnod[n[5]] == 1) {in[1]--; in[2]--; in[5]--;}
      if (gt->nadjelnod[n[6]] == 1) {in[2]--; in[3]--; in[6]--;}
      if (gt->nadjelnod[n[7]] == 1) {in[3]--; in[0]--; in[7]--;}
      break;
    }
    case 413:{
      if (surface_position (gt,n[0],n[1],n[2]) == 1) {in[0]--; in[1]--; in[2]--;}
      if (surface_position (gt,n[0],n[2],n[3]) == 1) {in[0]--; in[2]--; in[3]--;}
      if (surface_position (gt,n[0],n[3],n[1]) == 1) {in[0]--; in[3]--; in[1]--;}
      if (surface_position (gt,n[1],n[2],n[3]) == 1) {in[1]--; in[2]--; in[3]--;}
      break;
    }
    default:{
      fprintf (stderr,"\n\n wrong dimdegnne in function insidenod_assembling (%s, line %d)",__FILE__,__LINE__);
      break;
    }
    }
    
    for (j=0;j<nne;j++)  if (in[j] < insidenod[n[j]][0] && in[j] != 0)  insidenod[n[j]][0] = in[j];
  }
  
  
  // loop 2
  for (i=0;i<ne;i++){
    nne = gt->gelements[i].nne;
    n = gt->gelements[i].nodes;
    
    insidelem[i] = 1;
    for (j=0;j<nne;j++)
      if (insidenod[n[j]][0] < 0){
	insidelem[i] = 0;
	break;
      }
    
    switch (gt->gelements[i].auxinf){
    case 312:{
      if (insidenod[n[0]][0] == -1 && insidenod[n[1]][0] == -1 && insidenod[n[2]][0] == -1)
	fprintf (stderr,"\n\n !!!!!!!  wrong mesh - element %ld (%s, line %d)  !!!!!!!!! \n\n",i+1,__FILE__,__LINE__);
      
      if (insidenod[n[0]][0] == -2 || insidenod[n[1]][0] == -2 || insidenod[n[2]][0] == -2){
	if       (insidenod[n[0]][0] == -2) {jnod1 = n[1]; jnod2 = n[2]; sepnod = n[0];}
	else if  (insidenod[n[1]][0] == -2) {jnod1 = n[0]; jnod2 = n[2]; sepnod = n[1];}
	     else                           {jnod1 = n[0]; jnod2 = n[1]; sepnod = n[2];}
	
	for (j=0;j<gt->nadjelnod[jnod1];j++)
	  for (k=0;k<gt->nadjelnod[jnod2];k++)
	    if ( gt->adjelnod[jnod1][j] == gt->adjelnod[jnod2][k] && gt->adjelnod[jnod2][k] != i){
	      on = gt->adjelnod[jnod1][j];
	      j = k = 9999999;
	    }
	
	if      (gt->gelements[on].nodes[0] != jnod1 && gt->gelements[on].nodes[0] != jnod2) on = gt->gelements[on].nodes[0];
	else if (gt->gelements[on].nodes[1] != jnod1 && gt->gelements[on].nodes[1] != jnod2) on = gt->gelements[on].nodes[1];
	     else                                                                            on = gt->gelements[on].nodes[2];
	
	if (insidenod[on][0] < 0)
	  fprintf (stderr,"\n\n !!!!!!!  wrong mesh - element %ld (%s, line %d)  !!!!!!!!! \n\n",i+1,__FILE__,__LINE__);
	else
	  if (insidenod[on][0] == 1){
	    delete [] insidenod[on];
	    insidenod[on] = new long [10];
	    insidenod[on][0] = 2;
	  }
	  else
	    insidenod[on][0]++;
	
	insidenod[on][insidenod[on][0]-1] = sepnod;
      }
      break;
    }
    case 412:{
      sum = 0;
      for (j=0;j<nne;j++)
	if (insidenod[n[j]][0] < 0)
	  sum += insidenod[n[j]][0];
      
      if (sum < -4 || sum == -3 || (sum == -4 && insidenod[n[0]][0] == insidenod[n[1]][0]))
	fprintf (stderr,"\n\n !!!!!!!  wrong mesh - element %ld (%s, line %d)  !!!!!!!!! \n\n",i+1,__FILE__,__LINE__);
      
      break;
    }
    case 622:{
      for (j=3;j<6;j++)
	if      (insidenod[n[j]][0] ==  1)   insidenod[n[j]][0] =  0;
	else if (insidenod[n[j]][0] == -1)   insidenod[n[j]][0] = -7;
      
      if (insidenod[n[0]][0] == -1 && insidenod[n[1]][0] == -1 && insidenod[n[2]][0] == -1)
	fprintf (stderr,"\n\n !!!!!!!  wrong mesh - element %ld (%s, line %d)  !!!!!!!!! \n\n",i+1,__FILE__,__LINE__);
      
      if (insidenod[n[0]][0] == -2 || insidenod[n[1]][0] == -2 || insidenod[n[2]][0] == -2){
	if      (insidenod[n[4]][0] == 0) {jnod = n[4]; insidenod[n[4]][0] = -7; sepnodes[0] = n[0]; sepnodes[1] = n[3]; sepnodes[2] = n[5];}
	else if (insidenod[n[5]][0] == 0) {jnod = n[5]; insidenod[n[5]][0] = -7; sepnodes[0] = n[1]; sepnodes[1] = n[3]; sepnodes[2] = n[4];}
	     else                         {jnod = n[3]; insidenod[n[3]][0] = -7; sepnodes[0] = n[2]; sepnodes[1] = n[4]; sepnodes[2] = n[5];}
	
	joinel = (gt->adjelnod[jnod][0] == i) ? gt->adjelnod[jnod][1] : gt->adjelnod[jnod][0];
	
	if       (gt->gelements[joinel].nodes[4] == jnod) on = gt->gelements[joinel].nodes[0];
	else if  (gt->gelements[joinel].nodes[5] == jnod) on = gt->gelements[joinel].nodes[1];
	     else                                         on = gt->gelements[joinel].nodes[2];
	
	if (insidenod[on][0] < 0)
	  fprintf (stderr,"\n\n !!!!!!!  wrong mesh - element %ld (%s, line %d)  !!!!!!!!! \n\n",i+1,__FILE__,__LINE__);
	else
	  if (insidenod[on][0] == 1){
	    delete [] insidenod[on];
	    insidenod[on] = new long [30];
	    insidenod[on][0] = 4;
	  }
	  else
	    insidenod[on][0] += 3;
	
	for (k=0,j=insidenod[on][0]-3;j<insidenod[on][0];k++,j++)
	  insidenod[on][j] = sepnodes[k];
      }
      break;
    }
    case 822:{
      for (j=4;j<8;j++)
	if      (insidenod[n[j]][0] ==  1)   insidenod[n[j]][0] =  0;
	else if (insidenod[n[j]][0] == -1)   insidenod[n[j]][0] = -7;

      sum = 0;
      for (j=0;j<4;j++)
	if (insidenod[n[j]][0] < 0)
	  sum += insidenod[n[j]][0];
      
      if (sum < -4 || sum == -3 || (sum == -4 && insidenod[n[0]][0] == insidenod[n[1]][0]))
	fprintf (stderr,"\n\n !!!!!!!  wrong mesh - element %ld (%s, line %d)  !!!!!!!!! \n\n",i+1,__FILE__,__LINE__);

      break;
    }
    case 413:{
      if (insidenod[n[0]][0] == -3 || insidenod[n[1]][0] == -3 || insidenod[n[2]][0] == -3 || insidenod[n[3]][0] == -3){
	for (j=0,k=0;j<4;j++){
	  if (insidenod[n[j]][0] != -3) {
	    aux3[k] = n[j];
	    k++;
	  }
	  else n4 = n[j];
	}
	on = opposite_node (gt,aux3,i);

	if (insidenod[on][0] < 0)
	  fprintf (stderr,"\n\n !!!!!!!  wrong mesh - element %ld (%s, line %d)  !!!!!!!!! \n\n",i+1,__FILE__,__LINE__);
	else
	  if (insidenod[on][0] == 1){
	    delete [] insidenod[on];
	    insidenod[on] = new long [40];
	    insidenod[on][0] = 2;
	    insidenod[on][1] = n4;
	  }
	  else{
	    insidenod[on][insidenod[on][0]] = n4;
	    insidenod[on][0]++;
	  }
      }
      
      else if (insidenod[n[0]][0] != 0 && insidenod[n[1]][0] != 0 && insidenod[n[2]][0] != 0 && insidenod[n[3]][0] != 0){
	for (j=0;j<nne;j++)
	  aux4[j] = -2;
	
	for (j=0;j<nne;j++){
	  for (l=0,k=0;k<nne;k++)
          {
	    if (k!=j)
	      aux3[l++] = n[k];
          }
	  
	  on = opposite_node (gt,aux3,i);
	  
	  if (on>=0)
          {
	    if (insidenod[on][0]>0){
	      for (l=0;l<nne;l++)
		if (l!=j) aux4[l] = -1;
	      if (aux4[j] == -2)
		aux4[j] = on;
	    }
          }
	}
	
	for (j=0;j<nne;j++)
        {
	  if (aux4[j] >= 0)
          {
	    if (insidenod[aux4[j]][0] == 1){
	      delete [] insidenod[aux4[j]];
	      insidenod[aux4[j]] = new long [40];
	      insidenod[aux4[j]][0] = 2;
	      insidenod[aux4[j]][1] = n[j];
	    }
	    else {
	      insidenod[aux4[j]][insidenod[aux4[j]][0]] = n[j];
	      insidenod[aux4[j]][0]++;
	    }
          }
        }
      }
      break;
    }
    default:{
      fprintf (stderr,"\n\n wrong dimdegnne in function insidenod_assembling (%s, line %d)",__FILE__,__LINE__);
      break;
    }
    }
  }
  n = NULL;
}


/**
   
 */
void least_square::adjap_assembling (long nap,const long *parentel,const double **apcoord,long *nadjap,long **adjap)
{
  long i,j,k,nne,*n,*aux;
  double dist,ldist;
  
  aux = new long [nap];
  
  for (i=0;i<nap;i++){
    nne = gt->gelements[parentel[i]].nne;
    n = gt->gelements[parentel[i]].nodes;
    
    ldist = 1.0e100;
    for (j=0;j<nne;j++){
      if (insidenod[n[j]][0] != -1){
	dist = gt->gnodes[n[j]].distance2 (dim,apcoord[i]);
	if (dist < ldist){
	  ldist = dist;
	  aux[i] = n[j];
	}
      }
    }
  }
  
  for (i=0;i<nn;i++)
    if (insidenod[i][0] > 1)
      for (j=1;j<insidenod[i][0];j++)
	for (k=0;k<nap;k++)
	  if (aux[k] == insidenod[i][j])
	    aux[k] = i;
  
  // kontrola
  for (i=0;i<nap;i++)
    if (insidenod[aux[i]][0] < 1)
      fprintf (stderr,"\n\n something is wrong in function adjap_assembling\n\n\n");
  
  memset (nadjap,0,nn*sizeof(long));
  for (i=0;i<nap;i++)
    nadjap[aux[i]]++;
  
  for (i=0;i<nn;i++){
    adjap[i] = new long [nadjap[i]];
    nadjap[i] = 0;
  }
  
  for (i=0;i<nap;i++)
    adjap[aux[i]][nadjap[aux[i]]++] = i;
  
}



/**
   Function creates element patch around every vertex node and from sampling
   points inside of patch interpolates values to nodes.
   
   @param cut_flag - flag, whether to null values on the "null border" = border of region where known values are zero
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void least_square::compute_patches_spr (long cut_flag)
{
  long i,j,k,eid,auxl,ncoef;
  double auxd;
  vector coef_normcoord(2*dim),*a;
  ivector magnitude(nn);
  
  nullv (nodvalue,nn*nvals);
  
  switch (gt->gelements[0].auxinf%100){
  case 12:{  ncoef = 3; break; }
  case 22:{  ncoef = 6; break; }
  case 13:{  ncoef = 4; break; }
  default:{
    fprintf (stderr,"\n\n unknown element type is required in function compute_patches_spr (%s, line %d)",__FILE__,__LINE__);
    break;
  }}
  
  a = new vector[nvals];
  for (i=0;i<nvals;i++)
    reallocv (ncoef,a[i]);
  
  for (i=0;i<nn;i++){
    if (insidenod[i][0] < 1) continue;
    
    nadjel = gt->nadjelnod[i];
    adjel = gt->adjelnod[i];
    
    normal_coordinates_ae   (coef_normcoord);
    polynom_coefficients_ae (coef_normcoord,a);
    nodvalue_assembling_ae  (coef_normcoord,a,i,magnitude);
  }
  
  for (i=0;i<nn;i++){
    if (magnitude[i] > 1)
      for (j=0;j<nvals;j++)
	nodvalue[i+j*nn] /= (double)magnitude[i];
    else if (magnitude[i] == 0)
      fprintf (stderr,"\n !!!!!!!  stress of node %ld was not computed  !!!!!!!!! \n",i+1);
  }
  
  // border cutting
  if (cut_flag & 1){
    for (i=0;i<nn;i++){
      auxl=0;
      for (j=0;j<gt->nadjelnod[i];j++){
	eid=gt->adjelnod[i][j];
	auxd=0.0;
	for (k=0;k<nsp[eid]*nvals;k++)
	  auxd+=spvalue[eid][k];
	
	if (auxd==0.0) auxl++;
      }
      if ((gt->nadjelnod[i]-auxl)<=auxl)
	for (j=0;j<nvals;j++)
	  nodvalue[i+j*nn]=0.0;
      
    }
  }
  delete [] a;
}

/**
   Function creates element patch around every vertex node (typ_patch=='n') or inside element (typ_patch=='e')
   and from nodes on patch interpolates values to sampling points.
   
   @param typ_patch - flag for definition of type of patch
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void least_square::compute_patches_nod2sp (char typ_patch)
{
  long i,j,k;
  double **magnitude;
  vector coef_normcoord(2*dim),*a;
  
  // allocation
  magnitude = new double* [ne];
  for (i=0;i<ne;i++){
    magnitude[i] = new double [nsp[i]];
    memset (magnitude[i],0,nsp[i]*sizeof(double));
    memset (spvalue[i],0,nsp[i]*nvals*sizeof(double));
  }
  
  a = new vector[nvals];
  
  // computation
  if (typ_patch == 'n'){
    for (i=0;i<nn;i++){
      if (insidenod[i][0] < 1) continue;
      
      nadjel = gt->nadjelnod[i];
      adjel = gt->adjelnod[i];
      
      normal_coordinates_ae       (coef_normcoord);
      polynom_coefficients_inv_ae (coef_normcoord,a,'n');
      spvalue_assembling_patch_ae (coef_normcoord,a,-1,magnitude);
    }
  }
  else if (typ_patch == 'e'){
    for (i=0;i<ne;i++){
      if (insidelem[i]==0) continue;
      
      nadjel = gt->nadjelel[i];
      adjel = gt->adjelel[i];
      
      normal_coordinates_ae       (coef_normcoord);
      polynom_coefficients_inv_ae (coef_normcoord,a,'e');
      spvalue_assembling_patch_ae (coef_normcoord,a,i,magnitude);
    }
  }
  else fprintf (stderr,"\n\n typ_patch is wrong in function compute_patches_nod2sp (%s, line %d).\n",__FILE__,__LINE__);
  
  for (i=0;i<ne;i++)
    for (j=0;j<nsp[i];j++){
      if (magnitude[i][j] == 0.0)      fprintf (stderr,"\n !!!!!!!  stress of element %ld was not computed  !!!!!!!!! \n",i+1);
      for (k=0;k<nvals;k++)
	spvalue[i][j*nvals+k] /= magnitude[i][j];
    }
  
  // dealocation
  destrm (magnitude,ne);
  delete [] a;
}

/**
   Function creates element patch around every vertex node and from sampling
   points inside of patch interpolates values to answer points.
   
   created  4.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void least_square::compute_patches_sp2sp (long nap,const double **apcoord,const long *nadjap,const long **adjap,double **apvalue)
{
  long i,ncoef;
  vector coef_normcoord(2*dim),*a;
  
  a = new vector[nvals];
  for (i=0;i<nvals;i++)
    reallocv (3,a[i]);
  
  for (i=0;i<nn;i++){
    if (insidenod[i][0] < 1) continue;
    
    nadjel = gt->nadjelnod[i];
    adjel = gt->adjelnod[i];
    
    //if (nadjel>5)   ncoef = 6;   // dodelat pro libovolny prvek, toto je pro 312
    //else              ncoef = 3;
    
    ncoef = 3;
    
    normal_coordinates_ae   (coef_normcoord);
    polynom_coefficients_ae (coef_normcoord,a);
    apvalue_assembling      (coef_normcoord,a,nadjap[i],adjap[i],apcoord,apvalue);
  }
  
  delete [] a;
}


/**
   Function computes coeficients for computing of normed coordinates in patch.
   This function is called for every patch.
   
   @param coef_normcoord - answer = array of coeficients for computing of normed coordinates
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void least_square::normal_coordinates_ae (vector &coef_normcoord)
{
  long i,j;
  ivector mm(2*dim);
  
  fillv (adjel[0],mm);
  
  for (i=1;i<nadjel;i++)
    for (j=0;j<2*dim;j+=2){
      if (maxcoord[adjel[i]][j  ] > maxcoord[mm[j  ]][j  ]) mm[j  ] = adjel[i];
      if (maxcoord[adjel[i]][j+1] < maxcoord[mm[j+1]][j+1]) mm[j+1] = adjel[i];
    }
  
  for (i=0;i<2*dim;i+=2){
    coef_normcoord[i  ] = (maxcoord[mm[i]][i] + maxcoord[mm[i+1]][i+1])/(maxcoord[mm[i+1]][i+1] - maxcoord[mm[i]][i]);
    coef_normcoord[i+1] = 2.0/(maxcoord[mm[i]][i] - maxcoord[mm[i+1]][i+1]);
  }
}

/**
   Function computes 
   
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void least_square::polynom_coefficients_ae (const vector &coef_normcoord,vector *a)
{
  long i,j,k,l,m,ncoef=a[0].n;
  double zero=1.0e-20;
  vector normcoord(dim),pol(ncoef),contrv(ncoef);
  matrix ptp(ncoef,ncoef),contrm(ncoef,ncoef);
  
  for (i=0;i<nvals;i++)
    fillv (0.0,a[i]);
  
  for (i=0;i<nadjel;i++){
    for (j=0,l=0,m=0;j<nsp[adjel[i]];j++){
      for (k=0;k<dim;k++){
	normcoord[k] = coef_normcoord[2*k] + coef_normcoord[2*k+1] * spcoord[adjel[i]][l];
	l++;
      }
      
      polynom (ncoef,normcoord.a,pol.a);
      
      for (k=0;k<nvals;k++){
	cmulv (spvalue[adjel[i]][m],pol,contrv);
	addv (a[k],contrv,a[k]);
	m++;
      }
      
      vxv (pol,pol,contrm);
      addm (ptp,contrm,ptp);
    }
  }

  for (i=0;i<nvals;i++)
    gause (ptp,contrm,a[i],zero);
}

/**
   Function computes 
   
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void least_square::polynom_coefficients_inv_ae (const vector &coef_normcoord,vector *a,char typ_patch)
{
  long i,j,k,nid,nne,ncoef,nadjn,adjanod[51],*nodes;
  double zero=1.0e-20;

  nadjn = 0;
  for (i=0;i<nadjel;i++){
    nne = gt->gelements[adjel[i]].nne;
    nodes = gt->gelements[adjel[i]].nodes;
    for (j=0;j<nne;j++){
      for (k=0;k<nadjn;k++)
	if (nodes[j]==adjanod[k])   break;
      
      if (k==nadjn)
	adjanod[nadjn++]=nodes[j];
      
    }
  }
  nodes=NULL;
  
  // ncoef, reallocation of a[i]
  if (typ_patch=='n'){
    if (nadjel>4)   ncoef = 6;   // only for 312
    else            ncoef = 3;
  }
  else if (typ_patch=='e'){
    if (nadjel>10)  ncoef = 10;  // only for 312
    else            ncoef = 6;
  }
  else fprintf (stderr,"\n\n typ_patch is wrong in function polynom_coefficients_inv_ae (%s, line %d).\n",__FILE__,__LINE__);
  
  if (ncoef!=a[0].n)
    for (i=0;i<nvals;i++) { reallocv (ncoef,a[i]); }
  else
    for (i=0;i<nvals;i++)   nullv (a[i].a,a[i].n);
  
  // computation
  vector normcoord(dim),pol(ncoef),contrv(ncoef);
  matrix ptp(ncoef,ncoef),contrm(ncoef,ncoef);
  
  for (i=0;i<nadjn;i++){
    nid = adjanod[i];
    
    normcoord[0] = coef_normcoord[0] + coef_normcoord[1] * gt->gnodes[nid].x;
    normcoord[1] = coef_normcoord[2] + coef_normcoord[3] * gt->gnodes[nid].y;
    if (dim==3)
      normcoord[2] = coef_normcoord[4] + coef_normcoord[5] * gt->gnodes[nid].z;
    
    polynom (ncoef,normcoord.a,pol.a);
    
    for (j=0;j<nvals;j++){
      cmulv (nodvalue[nid+j*nn],pol,contrv);
      addv (a[j],contrv,a[j]);
    }
    
    vxv (pol,pol,contrm);
    addm (ptp,contrm,ptp);
  }
  
  for (i=0;i<nvals;i++)
    gause (ptp,contrm,a[i],zero);
}

/**
   Function computes values of `polynoms` in nodes over all patch.
   
   @param coef_normcoord - array of coeficients for computing of normed coordinates
   @param a              - array of vector of coefficients of `polynoms`
   @param id             - id of middle node
   @param magnitude      - 2D array of sums of magnitudes of values in sampling points
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void least_square::nodvalue_assembling_ae (const vector &coef_normcoord,const vector *a,long nid,ivector &magnitude)
{
  long i,j,nne,*nodes,ncoef=a[0].n;
  ivector aux(nn);
  
  sigma (ncoef,nid,coef_normcoord,a,magnitude);
  
  for (i=0;i<nadjel;i++){
    nne = gt->give_nne (adjel[i]);
    nodes = gt->gelements[adjel[i]].nodes;
    for (j=0;j<nne;j++)
      aux[nodes[j]]++;
  }
  nodes=NULL;
  
  for (i=1;i<insidenod[nid][0];i++)
    aux[insidenod[nid][i]] += 2;
  
  for (i=0;i<nn;i++)
    if ((aux[i] > 1 && insidenod[i][0] < 1) || (aux[i] > 0 && insidenod[i][0] < -1))
      sigma (ncoef,i,coef_normcoord,a,magnitude);
}

/**
   Function computes values of `polynoms` in sampling points over all patch.
   
   @param coef_normcoord - array of coeficients for computing of normed coordinates
   @param a              - array of vector of coefficients of `polynoms`
   @param id             - id of middle element; if patch is round node, id == -1
   @param magnitude      - 2D array of sums of magnitudes of values in sampling points
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void least_square::spvalue_assembling_patch_ae (const vector &coef_normcoord,const vector *a,long id,double **magnitude)
{
  long i,j,k,eid;
  double contr,mag;
  vector pol(a[0].n),normcoord(dim);
  
  for (i=0;i<nadjel;i++){
    eid = adjel[i];
    for (j=0;j<nsp[eid];j++){
      
      if (eid == id) mag = 10;
      else           mag = 1.0;
      
      magnitude[eid][j] += mag;
      
      for (k=0;k<dim;k++)
	normcoord[k] = coef_normcoord[2*k] + coef_normcoord[2*k+1] * spcoord[eid][k+j*dim];
      
      polynom (a[0].n,normcoord.a,pol.a);
      
      for (k=0;k<nvals;k++){
	scprd (pol,a[k],contr);
	spvalue[eid][j*nvals+k] += contr*mag;
      }
    }
  }
}


/**
   Function computes values of `polynoms` in answer points (over all patch ?????).
   
   @param coef_normcoord - array of coeficients for computing of normed coordinates
   @param a              - array of vector of coefficients of `polynoms`
   @param nadjap         -
   @param adjap          -
   @param 
   @param 
   
   created  4.5.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void least_square::apvalue_assembling (const vector &coef_normcoord,const vector *a,long nadjap,const long *adjap,const double **apcoord,double **apvalue)
{
  long i,j,ncoef=a[0].n;
  vector pol(ncoef),normcoord(dim);
  
  for (i=0;i<nadjap;i++){
    for (j=0;j<dim;j++)
      normcoord[j] = coef_normcoord[2*j] + coef_normcoord[2*j+1] * apcoord[adjap[i]][j];
    
    polynom (ncoef,normcoord.a,pol.a);
    
    for (j=0;j<nvals;j++)
      scprd (pol,a[j],apvalue[adjap[i]][j]);
    
  }
}

/**
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void least_square::sigma (long ncoef,long nid,const vector &coef_normcoord,const vector *a,ivector &magnitude)
{
  long i;
  double contr;
  vector pol(ncoef),normcoord(dim);
  
  magnitude[nid]++;
  
  normcoord[0] = coef_normcoord[0] + coef_normcoord[1] * gt->gnodes[nid].x;
  normcoord[1] = coef_normcoord[2] + coef_normcoord[3] * gt->gnodes[nid].y;
  if (dim == 3)
    normcoord[2] = coef_normcoord[4] + coef_normcoord[5] * gt->gnodes[nid].z;
  
  polynom (ncoef,normcoord.a,pol.a);
  
  for (i=0;i<nvals;i++){
    scprd (pol,a[i],contr);
    nodvalue[nid+i*nn] += contr;
  }
}

/**
   Function computes 
   
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void least_square::polynom (long ncoef,const double *normcoord,double *p)
{
  if (dim==2) {
    if (ncoef<= 0) { fprintf (stderr,"\n\n wrong ncoef in function polynom (%s, line %d).\n",__FILE__,__LINE__); exit (1); }
    if (ncoef>= 1)   p[ 0] = 1.0;
    if (ncoef>= 2)   p[ 1] = normcoord[0];
    if (ncoef>= 3)   p[ 2] = normcoord[1];
    if (ncoef>= 4)   p[ 3] = normcoord[0] * normcoord[1];
    if (ncoef>= 5)   p[ 4] = normcoord[0] * normcoord[0];
    if (ncoef>= 6)   p[ 5] = normcoord[1] * normcoord[1];
    if (ncoef>= 7)   p[ 6] = normcoord[0] * normcoord[0] * normcoord[1];
    if (ncoef>= 8)   p[ 7] = normcoord[0] * normcoord[1] * normcoord[1];
    if (ncoef>= 9)   p[ 8] = normcoord[0] * normcoord[0] * normcoord[0];
    if (ncoef>=10)   p[ 9] = normcoord[1] * normcoord[1] * normcoord[1];
    if (ncoef>=11)   p[10] = normcoord[0] * normcoord[0] * normcoord[0] * normcoord[1];
    if (ncoef>=12)   p[11] = normcoord[0] * normcoord[0] * normcoord[1] * normcoord[1];
    if (ncoef>=13)   p[12] = normcoord[0] * normcoord[1] * normcoord[1] * normcoord[1];
    if (ncoef>=14)   p[13] = normcoord[0] * normcoord[0] * normcoord[0] * normcoord[0];
    if (ncoef>=15)   p[14] = normcoord[1] * normcoord[1] * normcoord[1] * normcoord[1];
    if (ncoef >15) { fprintf (stderr,"\n\n wrong ncoef in function polynom (%s, line %d).\n",__FILE__,__LINE__); exit (1); }
  }
  else if (dim==3) {
    if (ncoef<= 0) { fprintf (stderr,"\n\n wrong ncoef in function polynom (%s, line %d).\n",__FILE__,__LINE__); exit (1); }
    if (ncoef>= 1)   p[0] = 1.0;
    if (ncoef>= 2)   p[1] = normcoord[0];
    if (ncoef>= 3)   p[2] = normcoord[1];
    if (ncoef>= 4)   p[3] = normcoord[2];
    if (ncoef > 4) { fprintf (stderr,"\n\n wrong ncoef in function polynom (%s, line %d).\n",__FILE__,__LINE__); exit (1); }
  }
  else { fprintf (stderr,"\n\n wrong dim in function polynom (%s, line %d).\n",__FILE__,__LINE__); exit (1); }
}


/*
  // *** PRINT ***
  
  //fprintf_l_2D (stdout,spr_est.insidenod,nn,1,"spr_est.insidenod",3);
  //for (i=0;i<ne;i++){
    //fprintf_d_1D (stdout,spr_est.maxcoord[i],4,"spr_est.insidenod",4);
    //fprintf_d_1D (stdout,spr_est.spcoord[i],2,"spr_est.insidenod",3);
    //fprintf_d_1D (stdout,spr_est.spvalue[i],3,"spr_est.insidenod",3);
  //}
  
  if (0){
    fprintf (stdout,"\n spr_est.insidenod:");
    for (i=0;i<nn;i++){
      fprintf (stdout,"\n %5ld",spr_est.insidenod[i][0]);
      if (spr_est.insidenod[i][0]>0)
	for (j=1;j<spr_est.insidenod[i][0];j++)
	  fprintf (stdout,"%5ld",spr_est.insidenod[i][j]);
    }
  }

   created  16.10.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
double least_square::give_magnitude (long eid,long nid,double *nc)
{
  double mag;
  
  switch (gt->gelements[eid].auxinf){
  case 312:{
    
    if      (nid == gt->gelements[eid].nodes[0])  mag = nc[0];
    else if (nid == gt->gelements[eid].nodes[1])  mag = nc[1];
    else if (nid == gt->gelements[eid].nodes[2])  mag = 1-nc[0]-nc[1];
    else fprintf (stderr,"\n\n something is wrong in function give_magnitude (%s, line %d).\n",__FILE__,__LINE__);
    
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown nnedegdim is required in function give_magnitude (%s, line %d).\n",__FILE__,__LINE__);
    return 0;
  }}
  
  if (mag < 0.05)  mag = 0.05;
  
  return mag;
}

   created  16.10.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
long least_square::give_nearnode (long eid,double *nc)
{
  long nid;
  
  switch (gt->gelements[eid].auxinf){
  case 312:{
    
    if      (0.5 < nc[0]        ) nid = gt->gelements[eid].nodes[0];
    else if (0.5 < nc[1]        ) nid = gt->gelements[eid].nodes[1];
    else if (0.5 < 1-nc[0]-nc[1]) nid = gt->gelements[eid].nodes[2];
    else                          nid = -1;
    
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown nnedegdim is required in function give_nearnode (%s, line %d).\n",__FILE__,__LINE__);
    return -2;
  }}
  
  return nid;
}

   Functions assembles lists of adjacent elements to each element of problem.
   Adjacent means elements share one edge.
   
   @ param nadjelel - 
   @ param adjelel  - 
   
   created  21.10.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
void least_square::adjelel_assembling (long *nadjelel,long **adjelel)
{
  long i,*n;
  
  for (i=0;i<ne;i++){
    n = gt->gelements[i].nodes;
    
    switch (gt->gelements[i].auxinf){
    case 312:{
      adjelel[i] = new long [nadjelel[i] = 3];
      adjelel[i][0] = adjelem2edge (gt,n[0],n[1],i);
      adjelel[i][1] = adjelem2edge (gt,n[1],n[2],i);
      adjelel[i][2] = adjelem2edge (gt,n[2],n[0],i);
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown nnedegdim is required in function adjelel_assembling (%s, line %d).\n",__FILE__,__LINE__);
    }}
  }
}

void least_square::nat_coord_assembling (double **nc)
{
  long i,j,nne,oldnne;
  vector x,y;
  
  oldnne = 3;
  reallocv (oldnne,x);
  reallocv (oldnne,y);
  
  for (i=0;i<ne;i++){
    nne = gt->give_nne (i);
    if (nne!=oldnne) { oldnne=nne; reallocv(nne,x); reallocv(nne,y); } 
    gt->give_node_coord2d (x,y,i);
    
    switch (gt->gelements[i].auxinf){
    case 312:{
      for (j=0;j<nsp[i];j++)
	nc_lin_3_2d (spcoord[i][j*2],spcoord[i][j*2+1],x.a,y.a,nc[i][j*2],nc[i][j*2+1]);
      
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown nnedegdim is required in function give_magnitude (%s, line %d).\n",__FILE__,__LINE__);
    }}
  }
}

   Function computes values of `polynoms` in sampling points over all patch.
   
   @param coef_normcoord - array of coeficients for computing of normed coordinates
   @param a              - array of vector of coefficients of `polynoms`
   @param id             - id of middle node or element 
   @param magnitude      - 2D array of sums of magnitudes of values in sampling points
   
   created  21.10.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz

void least_square::spvalue_assembling_patch2_ae (const vector &coef_normcoord,const vector *a,long id,double **magnitude)
{
  long i,j,k,eid;
  double contr,mag;
  vector pol(a[0].n),normcoord(dim);
  
  for (i=0;i<nadjel;i++){
    eid = adjel[i];
    if (eid != id && insidelem[eid]!=0)
      continue;
    
    for (j=0;j<nsp[eid];j++){
      
      mag = 1.0;
      magnitude[eid][j] += mag;
      
      for (k=0;k<dim;k++)
	normcoord[k] = coef_normcoord[2*k] + coef_normcoord[2*k+1] * spcoord[eid][k+j*dim];
      
      polynom (a[0].n,normcoord.a,pol.a);
      
      for (k=0;k<nvals;k++){
	scprd (pol,a[k],contr);
	spvalue[eid][j*nvals+k] += contr*mag;
      }
    }
  }
}
*/
