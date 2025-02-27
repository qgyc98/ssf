#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "gadaptivity.h"
#include "gtopology.h"
#include "vector.h"
#include "matrix.h"
#include "basefun.h"
#include "mathem.h"



/*
// max number of nodes on element is 8
*/
/// CONSTRUCTOR
patch_averaging :: patch_averaging (gtopology *gt, long dim, long nvals, long flag)
{
  patch_averaging::gt = gt;
  
  patch_averaging::dim = dim;
  patch_averaging::nvals = nvals;
  nn = gt->nn;
  ne = gt->ne;
  
  spcoord = maxcoord = NULL;
  insidenod = NULL;
  
  patch_averaging::flag = flag;
}

patch_averaging :: ~patch_averaging (void)
{
  if (insidenod != NULL)  for (long i=0;i<nn;i++)  delete [] insidenod[i];
  if (maxcoord  != NULL)  for (long i=0;i<ne;i++)  delete [] maxcoord[i];
  if (spcoord   != NULL)  for (long i=0;i<ne;i++)  delete [] spcoord[i];
  // 
  // delete [] nsp;
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
void patch_averaging :: solve (FILE *out, const matrix *spvalue, vector *nodvalue)
{
  fprintf (stdout,"\n VVV===========  SPR_smoothing  ===========VVV");
  
  if (gt->nadjelnod == NULL)
    gt->adjacelem (out);
  
  spcoord = new double* [ne];
  
  nsp_spcoord_maxcoord_assembling (spvalue, 'a');
  insidenod_assembling            ();
  compute_patches_spr             (spvalue, nodvalue, 0);
  
  fprintf (stdout,"\n AAA===========  SPR_smoothing  ===========AAA");
}


/**
   Function fills up arrays:
    'nsp' (number of sampling points) by default values,
    'spcoord' by natural coordinates of default sampling points,
    'maxcoord' by extreme element coordinates.
   
   @param nsma_flag - determines which array will be actually filled: 'a' = all arrays, 'm' = only array maxcoord
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void patch_averaging :: nsp_spcoord_maxcoord_assembling (const matrix *spvalue, char nsma_flag)
{
  long i,j,k,nne, nsp;
  double xi,eta;
  vector xyz[3], bf;
  allocv (20, xyz[0]);
  allocv (20, xyz[1]);
  allocv (20, xyz[2]);
  
  maxcoord = new double* [ne];
  
  ///
  for (i=0; i<ne; i++) {
    xyz[0].n = xyz[1].n = xyz[2].n = nne = gt->gelements[i].nne;
    gt->give_node_coord3d (xyz[0], xyz[1], xyz[2], i);
    
    nsp = spvalue[i].m;
    spcoord[i] = new double[nsp*dim];
    
    if (nsma_flag == 'a'){
      
      if (nsp==1) {
	for (j=0; j<dim; j++) {
	  spcoord[i][j] = 0.0;
	  for (k=0; k<nne; k++)  spcoord[i][j] += xyz[j][k];
	  spcoord[i][j] /= nne;
	}
      }
      else {
	allocv (nne, bf);
	
	switch (gt->gelements[i].auxinf) {
	case 622:{
	  if (nsp != 3)  { print_err("wrong mesh - element %ld", __FILE__, __LINE__, __func__, i+1);  exit (1); }
	  
	  xi = eta = 1.0/6.0;
	  for (j=0; j<3; j++) {
	    if (j==1)  xi=4.0/6.0;
	    if (j==2) {xi=1.0/6.0; eta=4.0/6.0;}
	    
	    bf_quad_3_2d (bf.a,xi,eta);
	    scprd (bf, xyz[0], spcoord[i][2*j  ]);
	    scprd (bf, xyz[1], spcoord[i][2*j+1]);
	  }
	  break;
	}
	case 412:
	case 822: {
	  if (nsp != 4)  { print_err("wrong mesh - element %ld", __FILE__, __LINE__, __func__, i+1);  exit (1); }
	  
	  xi = eta = 0.577350269189626;
	  for (j=0; j<2; j++) {
	    xi = -xi;
	    for (k=0; k<2; k++) {
	      eta = -eta;
	      
	      if (gt->gelements[i].auxinf == 412)  bf_lin_4_2d  (bf.a, xi, eta);
	      if (gt->gelements[i].auxinf == 822)  bf_quad_4_2d (bf.a, xi, eta);
	      
	      scprd (bf, xyz[0], spcoord[i][4*j + 2*k  ]);
	      scprd (bf, xyz[1], spcoord[i][4*j + 2*k+1]);
	    }
	  }
	  break;
	}
	default: { print_err("wrong dimdegnne in function nsp_spcoord_maxcoord_assembling", __FILE__, __LINE__, __func__);  exit (1); }
	}
	
	destrv (bf);
      }
    }
    
    maxcoord[i] = new double[2*dim];    
    
    if (nsma_flag == 'a' || nsma_flag == 'm'){
      
      switch (gt->gelements[i].auxinf){
      case 312:
      case 622:{
	maxmin_3 (xyz[0].a,maxcoord[i][0],maxcoord[i][1]);
	maxmin_3 (xyz[1].a,maxcoord[i][2],maxcoord[i][3]);
	break;
      }
      case 412:
      case 822:
      case 413:{
	;           maxmin_4 (xyz[0].a, maxcoord[i][0], maxcoord[i][1]);
	;           maxmin_4 (xyz[1].a, maxcoord[i][2], maxcoord[i][3]);
	if (dim==3) maxmin_4 (xyz[2].a, maxcoord[i][4], maxcoord[i][5]);
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
void patch_averaging :: insidenod_assembling (void)
{
  long i,j,k,l,nne,*n,in[8];            // max number of nodes on element is 8
  long jnod1,jnod2,sepnod,on,sum;
  long jnod,joinel,sepnodes[3];
  long aux3[3],n4,aux4[4];
  
  //
  insidenod = new long* [nn];
  for (i=0; i<nn; i++) {
    insidenod[i] = new long [1];
    insidenod[i][0] = 1;
  }
  
  ///
  /// loop 1
  ///
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
  
  
  ///
  /// loop 2
  ///
  //insidelem = new long [ne];
  
  for (i=0;i<ne;i++){
    nne = gt->gelements[i].nne;
    n = gt->gelements[i].nodes;
    
    //insidelem[i] = 1;
    //for (j=0;j<nne;j++)
    //  if (insidenod[n[j]][0] < 0){
    //	insidelem[i] = 0;
    //	break;
    //  }
    
    switch (gt->gelements[i].auxinf){
    case 312:{
      if (insidenod[n[0]][0] == -1 && insidenod[n[1]][0] == -1 && insidenod[n[2]][0] == -1)
	{ print_err("wrong mesh - element %ld", __FILE__, __LINE__, __func__, i+1);  exit (1); }
      
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
	
	if (insidenod[on][0] < 0)  { print_err("wrong mesh - element %ld", __FILE__, __LINE__, __func__, i+1);  exit (1); }
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
	{ print_err("wrong mesh - element %ld", __FILE__, __LINE__, __func__, i+1);  exit (1); }
      
      break;
    }
    case 622:{
      for (j=3;j<6;j++)
	if      (insidenod[n[j]][0] ==  1)   insidenod[n[j]][0] =  0;
	else if (insidenod[n[j]][0] == -1)   insidenod[n[j]][0] = -7;
      
      if (insidenod[n[0]][0] == -1 && insidenod[n[1]][0] == -1 && insidenod[n[2]][0] == -1)
	{ print_err("wrong mesh - element %ld", __FILE__, __LINE__, __func__, i+1);  exit (1); }
      
      if (insidenod[n[0]][0] == -2 || insidenod[n[1]][0] == -2 || insidenod[n[2]][0] == -2){
	if      (insidenod[n[4]][0] == 0) {jnod = n[4]; insidenod[n[4]][0] = -7; sepnodes[0] = n[0]; sepnodes[1] = n[3]; sepnodes[2] = n[5];}
	else if (insidenod[n[5]][0] == 0) {jnod = n[5]; insidenod[n[5]][0] = -7; sepnodes[0] = n[1]; sepnodes[1] = n[3]; sepnodes[2] = n[4];}
	     else                         {jnod = n[3]; insidenod[n[3]][0] = -7; sepnodes[0] = n[2]; sepnodes[1] = n[4]; sepnodes[2] = n[5];}
	
	joinel = (gt->adjelnod[jnod][0] == i) ? gt->adjelnod[jnod][1] : gt->adjelnod[jnod][0];
	
	if       (gt->gelements[joinel].nodes[4] == jnod) on = gt->gelements[joinel].nodes[0];
	else if  (gt->gelements[joinel].nodes[5] == jnod) on = gt->gelements[joinel].nodes[1];
	     else                                         on = gt->gelements[joinel].nodes[2];
	
	if (insidenod[on][0] < 0)  { print_err("wrong mesh - element %ld", __FILE__, __LINE__, __func__, i+1);  exit (1); }
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
	{ print_err("wrong mesh - element %ld", __FILE__, __LINE__, __func__, i+1);  exit (1); }
      
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
	
	if (insidenod[on][0] < 0)  { print_err("wrong mesh - element %ld", __FILE__, __LINE__, __func__, i+1);  exit (1); }
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
   Function creates element patch around every vertex node and from sampling
   points inside of patch interpolates values to nodes.
   
   @param cut_flag - flag, whether to null values on the "null border" = border of region where known values are zero
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void patch_averaging :: compute_patches_spr (const matrix *spvalue, vector *nodvalue, long cut_flag)
{
  long i,j,k,eid,auxl,ncoef, nadjel;
  const long *adjel;
  double auxd;
  vector coef_normcoord(2*dim),*a;
  ivector magnitude(nn);
  
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
    allocv (ncoef,a[i]);
  
  for (i=0;i<nn;i++){
    if (insidenod[i][0] < 1) continue;
    
    nadjel = gt->nadjelnod[i];
    adjel = gt->adjelnod[i];
    
    normal_coordinates_ae   (nadjel, adjel, coef_normcoord);
    polynom_coefficients_ae (nadjel, adjel, coef_normcoord, a, spvalue);
    nodvalue_assembling_ae  (nadjel, adjel, coef_normcoord, a, nodvalue, i, magnitude);
  }
  
  for (i=0; i<nn; i++)
    if (magnitude[i] > 1)
      for (j=0; j<nvals; j++)
	nodvalue[i].a[j] /= (double)magnitude[i];
    else
      if (magnitude[i] == 0)
	{ print_err("stress of node %ld was not computed", __FILE__, __LINE__, __func__, i+1);  exit (1); }
  
  // border cutting
  if (cut_flag & 1){
    for (i=0;i<nn;i++){
      auxl=0;
      for (j=0;j<gt->nadjelnod[i];j++){
	eid=gt->adjelnod[i][j];
	auxd=0.0;
	for (k=0; k<spvalue[eid].m*spvalue[eid].n; k++)
	  auxd += spvalue[eid].a[k];
	
	if (auxd==0.0) auxl++;
      }
      if ((gt->nadjelnod[i]-auxl)<=auxl)
	for (j=0;j<nvals;j++)
	  nodvalue[i].a[j]=0.0;
      
    }
  }
  
  for (i=0;i<nvals;i++)  destrv (a[i]);
  delete [] a;
}

/**
   Function computes coeficients for computing of normed coordinates in patch.
   This function is called for every patch.
   
   @param coef_normcoord - answer = array of coeficients for computing of normed coordinates
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void patch_averaging :: normal_coordinates_ae (long nadjel, const long *adjel, vector &coef_normcoord)
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
void patch_averaging :: polynom_coefficients_ae (long nadjel, const long *adjel, const vector &coef_normcoord, vector *a, const matrix *spvalue)
{
  long i, j, k, l, ncoef=a[0].n;
  double zero=1.0e-20;
  vector normcoord(dim),pol(ncoef),contrv(ncoef);
  matrix ptp(ncoef,ncoef),contrm(ncoef,ncoef);
  
  for (i=0;i<nvals;i++)
    fillv (0.0,a[i]);
  
  for (i=0; i<nadjel; i++) {
    for (j=0,l=0; j<spvalue[adjel[i]].m; j++) {
      for (k=0;k<dim;k++){
	normcoord[k] = coef_normcoord[2*k] + coef_normcoord[2*k+1] * spcoord[adjel[i]][l];
	l++;
      }
      
      polynom (ncoef,normcoord.a,pol.a);
      
      for (k=0; k<nvals; k++) {
	cmulv (spvalue[adjel[i]](j,k), pol, contrv);
	addv (a[k],contrv,a[k]);
      }
      
      vxv (pol,pol,contrm);
      addm (ptp,contrm,ptp);
    }
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
void patch_averaging :: nodvalue_assembling_ae (long nadjel, const long *adjel, const vector &coef_normcoord, const vector *a, vector *nodvalue, long nid, ivector &magnitude)
{
  long i,j,nne,*nodes,ncoef=a[0].n;
  ivector aux(nn);
  
  sigma (nodvalue, ncoef,nid,coef_normcoord,a,magnitude);
  
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
      sigma (nodvalue, ncoef,i,coef_normcoord,a,magnitude);
}


/**
   Function computes 
   
   
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void patch_averaging :: polynom (long ncoef,const double *normcoord,double *p)
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



/**
   created  3.3.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void patch_averaging :: sigma (vector *nodvalue, long ncoef, long nid, const vector &coef_normcoord, const vector *a, ivector &magnitude)
{
  long i;
  double contr;
  vector pol(ncoef),normcoord(dim);
  
  magnitude[nid]++;
  
  ;              normcoord[0] = coef_normcoord[0] + coef_normcoord[1] * gt->gnodes[nid].x;
  ;              normcoord[1] = coef_normcoord[2] + coef_normcoord[3] * gt->gnodes[nid].y;
  if (dim == 3)  normcoord[2] = coef_normcoord[4] + coef_normcoord[5] * gt->gnodes[nid].z;
  
  polynom (ncoef, normcoord.a, pol.a);
  
  for (i=0; i<nvals; i++) {
    scprd (pol, a[i], contr);
    nodvalue[nid].a[i] += contr;
  }
}































/**
   Function interpolates values (stress or strain = derivatives) from nodes of element into selected point on element.
   
   @param bf - approximation functions on element
   @param rderfull - 1Darray of values in nodes; dimmension is ncomp*gt->nn after this manner {(val[0]; ... ;val[gt->nn])[0] ; ... ; (val[0]; ... ;val[gt->nn])[ncomp]}
   @param nodes - 1Darray of node numbers of element
   @param nn - total number of nodes on domain
   @param der_star - empty(returned) vector, dimension is ncomp
   
   created 19.8.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
      //   votoc "vector &der_star" za "long nn"
void give_der_star (vector &bf, vector *rderfull, ivector &nodes, vector &der_star, long /*nn*/)
{
  long i,j,nne,ncomp;
  
  nne = nodes.n;
  ncomp = der_star.n;
  
  fillv (0.0,der_star);
  
  for (i=0; i<ncomp; i++)
    for (j=0; j<nne; j++)
      der_star[i] += bf[j] * rderfull[nodes[j]].a[i];
  
}

void ntnmtov (matrix &ntnm,vector &ntnv)
{
  long i,j,nne;
  double mij,mii;

  nne = ntnv.n;

  mij = 0;
  for (i=0;i<nne;i++)
    for (j=0;j<nne;j++)
      mij += ntnm[i][j];

  mii = 0;
  for (i=0;i<nne;i++)
    mii += ntnm[i][i];

  for (i=0;i<nne;i++)
    ntnv[i] = ntnm[i][i] * mij / mii;
}


//  vytiskne do souboru typfile isolinie hodnot na prvku(valnod)
void print_contoures (gtopology *gt,char *file,double **midpoints,double *valel)
{
  long i,j,nmidp,nn,ne;
  long *nadjelnod,**adjelnod,**auxnod,*nodes;
  double **auxxyzv;
  vector valnod;
  
  nn = gt->nn;
  ne = gt->ne;
  
  auxxyzv = new double* [nn+ne];
  auxnod = new long* [ne];
  
  allocv (nn,valnod);
  
  // computing of valnod
  nadjelnod = gt->nadjelnod;
  adjelnod = gt->adjelnod;
  
  for (i=0;i<nn;i++){
    for (j=0;j<nadjelnod[i];j++){
      valnod[i] += valel[adjelnod[i][j]]; 
    }
    valnod[i] /= (double)nadjelnod[i]; 
  }
  
  nadjelnod = NULL;
  adjelnod = NULL;
  
  // assemble aux*
  for (i=0;i<nn;i++){
    auxxyzv[i] = new double [4];
    auxxyzv[i][0] = gt->gnodes[i].x;
    auxxyzv[i][1] = gt->gnodes[i].y;
    auxxyzv[i][2] = gt->gnodes[i].z;
    auxxyzv[i][3] = valnod[i];
  }
  for (i=0;i<ne;i++){
    if (midpoints[i] != NULL){
      auxxyzv[i+nn] = new double [4];
      auxxyzv[i+nn][0] = midpoints[i][0];
      auxxyzv[i+nn][1] = midpoints[i][1];
      auxxyzv[i+nn][2] = 0.0;
      auxxyzv[i+nn][3] = valel[i];
    }
    else
      auxxyzv[i+nn] = NULL;
  }
  
  nmidp = nn;
  for (i=0;i<ne;i++){
    nodes = gt->gelements[i].nodes;
    nmidp++;
    
    switch (gt->gelements[i].auxinf){
    case 312:{
      auxnod[i] = new long [4];
      auxnod[i][0] = 1;   auxnod[i][1] = nodes[0]+1;   auxnod[i][2] = nodes[1]+1;   auxnod[i][3] = nodes[2]+1;
      break;
    }
    case 412:{
      auxnod[i] = new long [13];
      auxnod[i][0] = 4;   auxnod[i][ 1] = nodes[0]+1;   auxnod[i][ 2] = nodes[1]+1;   auxnod[i][ 3] = nmidp;
                          auxnod[i][ 4] = nodes[1]+1;   auxnod[i][ 5] = nodes[2]+1;   auxnod[i][ 6] = nmidp;
                          auxnod[i][ 7] = nodes[2]+1;   auxnod[i][ 8] = nodes[3]+1;   auxnod[i][ 9] = nmidp;
                          auxnod[i][10] = nodes[3]+1;   auxnod[i][11] = nodes[0]+1;   auxnod[i][12] = nmidp;
      break;
    }
    case 622:{
      auxnod[i] = new long [13];
      auxnod[i][0] = 4;   auxnod[i][ 1] = nodes[0]+1;   auxnod[i][ 2] = nodes[3]+1;   auxnod[i][ 3] = nodes[5]+1;
                          auxnod[i][ 4] = nodes[3]+1;   auxnod[i][ 5] = nodes[1]+1;   auxnod[i][ 6] = nodes[4]+1;
                          auxnod[i][ 7] = nodes[4]+1;   auxnod[i][ 8] = nodes[2]+1;   auxnod[i][ 9] = nodes[5]+1;
                          auxnod[i][10] = nodes[3]+1;   auxnod[i][11] = nodes[4]+1;   auxnod[i][12] = nodes[5]+1;
      break;
    }
    case 822:{
      auxnod[i] = new long [25];
      auxnod[i][0] = 8;   auxnod[i][ 1] = nodes[0]+1;   auxnod[i][ 2] = nodes[4]+1;   auxnod[i][ 3] = nmidp;
                          auxnod[i][ 4] = nodes[4]+1;   auxnod[i][ 5] = nodes[1]+1;   auxnod[i][ 6] = nmidp;
                          auxnod[i][ 7] = nodes[1]+1;   auxnod[i][ 8] = nodes[5]+1;   auxnod[i][ 9] = nmidp;
                          auxnod[i][10] = nodes[5]+1;   auxnod[i][11] = nodes[2]+1;   auxnod[i][12] = nmidp;
                          auxnod[i][13] = nodes[2]+1;   auxnod[i][14] = nodes[6]+1;   auxnod[i][15] = nmidp;
                          auxnod[i][16] = nodes[6]+1;   auxnod[i][17] = nodes[3]+1;   auxnod[i][18] = nmidp;
                          auxnod[i][19] = nodes[3]+1;   auxnod[i][20] = nodes[7]+1;   auxnod[i][21] = nmidp;
                          auxnod[i][22] = nodes[7]+1;   auxnod[i][23] = nodes[0]+1;   auxnod[i][24] = nmidp;
      break;
    }
    case 413:{
      auxnod[i] = new long [5];
      auxnod[i][0] = 1;   auxnod[i][1] = nodes[0]+1;   auxnod[i][2] = nodes[1]+1;   auxnod[i][3] = nodes[2]+1;   auxnod[i][4] = nodes[3]+1;
      break;
    }
    default:{
      fprintf (stderr,"\n\n wrong dimdegnne in function");
      fprintf (stderr," print_contoures (%s, line %d)",__FILE__,__LINE__);
      break;
    }
    }
  }
  
  print_confile (file,nn,ne,auxxyzv,auxnod,gt->gelements[0].auxinf % 10);
  
  nodes = NULL;
  for (i=0;i<nn+ne;i++)
    delete [] auxxyzv[i];
  for (i=0;i<ne;i++)
    delete [] auxnod[i];
  delete [] auxxyzv;  delete [] auxnod;
}

void print_confile (const char *file,long nnod,long nel,double **auxxyzv,long **auxnod,long dim)
{
  long i,j,nbgn,nbgtr,nbgtet;
  ivector newnum(nnod+nel);
  FILE *con;

  con = fopen (file,"w");                    if (con==NULL) fprintf (stderr,"\n .con file has not been specified.\n Try it again!\n");

  fprintf (con,"00000000010000000002\n");
  fprintf (con,"0000000001000000000200000000030000000004\n\n");

  nbgn = 0;
  for (i=0;i<nnod+nel;i++){
    if (auxxyzv[i] != NULL){
      fprintf (con,"%6ld%15.4f%15.4f%15.4f%15.4f\n",++nbgn,auxxyzv[i][0],auxxyzv[i][1],auxxyzv[i][2],auxxyzv[i][3]);
      newnum[i] = nbgn; 
    }
  }

  fprintf (con,"\n");

  nbgtr = nbgtet = 0;
  if (dim == 2)
    for (i=0;i<nel;i++)
      for (j=0;j<auxnod[i][0];j++)
	fprintf (con,"%6ld%6ld%6ld%6ld\n",++nbgtr,newnum[auxnod[i][j*3+1]-1],newnum[auxnod[i][j*3+2]-1],newnum[auxnod[i][j*3+3]-1]);
  if (dim == 3)
    for (i=0;i<nel;i++)
      for (j=0;j<auxnod[i][0];j++)
	fprintf (con,"%6ld%6ld%6ld%6ld%6ld\n",++nbgtet,newnum[auxnod[i][j*3+1]-1],newnum[auxnod[i][j*3+2]-1],newnum[auxnod[i][j*3+3]-1],newnum[auxnod[i][j*3+4]-1]);

  rewind (con);
  fprintf (con,"%10d%10d\n",3,1);
  fprintf (con,"%10ld%10d%10ld%10ld\n",nbgn,0,nbgtr,nbgtet);

  fclose (con);
}


// ************************************
//          OTHER FUNCTIONS
// ************************************
void fprintf_matrix (FILE *stream, matrix &mx, char s[])
{
  long m = mx.m;
  long n = mx.n;   fprintf(stream,"\n %s(%ld,%ld):",s,m,n);
  long i,j;
  for (i=0;i<m;i++){
    for(j=0;j<n;j++){
      if (j%18==0)    fprintf(stream," \n");
      fprintf(stream,"%12.4f",mx[i][j]);
    }
    fprintf(stream,"  **");
  }
}

void fprintf_vector (FILE *stream,vector &v,char s[],long c)
{
  long n = v.n - 1;
  if (c==0) c = 1000000;
  fprintf (stream,"\n %s(%ld):\n [ 1]: ",s,v.n);
  for (long i=0; i<n; i++){
    fprintf (stream,"%15.8f",v[i]);
    if ((i+1)%c==0) fprintf(stream,"\n [%2ld]: ",i+2);
  }
  fprintf (stream,"%15.8f",v[n]);
}

void fprintf_ivector (FILE *stream,ivector &v,char s[],long c)
{
  long n = v.n - 1;
  if (c==0) c = 1000000;
  fprintf (stream,"\n %s(%ld):\n [ 1]: ",s,v.n);
  for (long i=0; i<n; i++){
    fprintf (stream,"%5ld",v[i]);
    if ((i+1)%c==0) fprintf(stream,"\n [%2ld]: ",i+2);
  }
  fprintf (stream,"%5ld",v[n]);
}

void fprintf_d_1D (FILE *stream,double *p,long p_n,char s[],long c)
{
  long n = p_n - 1;
  if (c==0) c = 1000000;
  fprintf (stream,"\n %s(%ld):\n [ 1]: ",s,p_n);
  for (long i=0; i<n; i++){
    fprintf (stream,"%14.6f",p[i]);
    if ((i+1)%c==0) fprintf(stream,"\n [%2ld]: ",i+2);
  }
  fprintf (stream,"%14.6f",p[n]);
}

void fprintf_l_1D (FILE *stream,long *p,long p_n,char s[],long c)
{
  long n = p_n - 1;
  if (c==0) c = 1000000;
  fprintf (stream,"\n %s(%ld):\n [ 1]: ",s,p_n);
  for (long i=0; i<n; i++){
    fprintf (stream,"%5ld",p[i]);
    if ((i+1)%c==0) fprintf(stream,"\n [%2ld]: ",i+2);
  }
  fprintf (stream,"%5ld",p[n]);
}

void fprintf_l_2D (FILE *stream,long **p,long p_n,long p_m,char s[],long c)
{
  long i,j;

  for (i=0;i<p_n;i++){
    if (c==0) c = 1000000;
    fprintf (stream,"\n %s[%2ld]:",s,i);
    for (j=0;j<p_m;j++){
      fprintf (stream,"%5ld",p[i][j]);
      if ((j+1)%c==0 && (j+1)!=p_m) fprintf(stream,"\n   ");
    }
  }
}


/**
   Function detects count of elements containing linear edge defined by two nodes.
   Array gt->adjelnod has to be in ascending order.
   count == 1 => the edge is on the boundary of domain
   count == 2 => the edge is inside of domain
   
   @param gt - gtopology of domain
   @param node1,node2 - nodes of the edge
   
   @return - count of elements containing the edge
   
   created 19.8.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
long edge_position (gtopology *gt,long node1,long node2)
{
  long c,a1,a2,n1,n2,*p1,*p2;
  
  a1=a2=c=0;
  n1 = gt->nadjelnod[node1];      p1 = gt->adjelnod[node1];
  n2 = gt->nadjelnod[node2];	 p2 = gt->adjelnod[node2];
  
  for (;;)
    if (p1[a1] > p2[a2])        { a2++;             if (a2 == n2)                        return (c); }
    else  if (p1[a1] == p2[a2]) { a1++; a2++; c++;  if (c  == 2 || a2 == n2 || a1 == n1) return (c); }
          else                  { a1++;             if (a1 == n1)                        return (c); }
}

/**
   Function detects count of elements containing linear surface defined by three nodes.
   Array gt->adjelnod has to be in ascending order.
   count == 1 => the surface is on the boundary of domain
   count == 2 => the surface is inside of domain
   
   @param gt - gtopology of domain
   @param node1,node2,node3 - nodes of the surface
   
   @return - count of elements containing the surface
   
   created 19.8.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
long surface_position (gtopology *gt,long node1,long node2,long node3)
{
  long c,a1,a2,a3,n1,n2,n3,*p1,*p2,*p3;
  
  a1=a2=a3=c=0;
  n1 = gt->nadjelnod[node1];      p1 = gt->adjelnod[node1];
  n2 = gt->nadjelnod[node2];	 p2 = gt->adjelnod[node2];
  n3 = gt->nadjelnod[node3];	 p3 = gt->adjelnod[node3];
  
  for (;;)
    if (p1[a1] > p2[a2])                    { a2++;                   if (a2 == n2)                                    return (c); }
    else  if (p1[a1] < p2[a2])              { a1++;                   if (a1 == n1)                                    return (c); }
          else  if (p1[a1] > p3[a3])        { a3++;                   if (a3 == n3)                                    return (c); }
                else  if (p1[a1] == p3[a3]) { a1++; a2++; a3++; c++;  if (c  == 2 || a1 == n1 || a2 == n2 || a3 == n3) return (c); }
                      else                  { a1++; a2++;             if (a1 == n1 || a2 == n2)                        return (c); }
}

/**
   created  21.10.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
long adjelem2edge (gtopology *gt,long node1,long node2,long eid)
{
  long a1,a2,n1,n2,*p1,*p2;
  
  a1=a2=0;
  n1 = gt->nadjelnod[node1];      p1 = gt->adjelnod[node1];
  n2 = gt->nadjelnod[node2];	 p2 = gt->adjelnod[node2];
  
  for (;;)
    if      (p1[a1] > p2[a2]) { a2++;       if (a2 == n2)             return (-1); }
    else if (p1[a1] < p2[a2]) { a1++;       if (a1 == n1)             return (-1); }
    else if (p1[a1] == eid)   { a1++; a2++; if (a2 == n2 || a1 == n1) return (-1); }
    else return (p1[a1]);
}

long opposite_node (gtopology *gt,long *nod,long eid)
{
  long i,j,k;
  long *nodes;
  
  nodes = NULL;
  
  for (i=0;i<gt->nadjelnod[nod[0]];i++)
    for (j=0;j<gt->nadjelnod[nod[1]];j++)
      for (k=0;k<gt->nadjelnod[nod[2]];k++)
	if ( gt->adjelnod[nod[0]][i] == gt->adjelnod[nod[1]][j] && gt->adjelnod[nod[1]][j] == gt->adjelnod[nod[2]][k] && gt->adjelnod[nod[0]][i] != eid){
	  nodes = gt->gelements[gt->adjelnod[nod[0]][i]].nodes;
	  i = j = k = 999999;
	}
  
  if (nodes == NULL)                                                     return (-1);
  if (nodes[0] != nod[0] && nodes[0] != nod[1] && nodes[0] != nod[2])    return (nodes[0]);
  if (nodes[1] != nod[0] && nodes[1] != nod[1] && nodes[1] != nod[2])    return (nodes[1]);
  if (nodes[2] != nod[0] && nodes[2] != nod[1] && nodes[2] != nod[2])    return (nodes[2]);
  return (nodes[3]);
}

/**
   
   created 19.8.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
   
   long common_numbers (long n1,long *p1,long n2,long *p2)
   {
   long c=0;
   register int i,j;
   
   for (i=0;i<n1;i++)
   for (j=0;j<n2;j++)
   if (p1[i]==p2[j]){
   c++;
   break;
   }
   
   return c;
   }
*/


/**
   function compiles data(from array 'valel') for printing file.dx(.ex)
   
   @param gt - gtopology
   @param path - path
   @param file - name of "dx"("ex") file; if file=NULL file is generated automatically from 'caption'
   @param caption - name=caption of 'value'
   @param valel - array of 'values' on elements, one element = one value, => dimension is gt->ne
   @param flag - type of file - "d" = file.dx(for OpenDx) , "e" = file.ex(for Elixir)
   
   created  20.11.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void print_valel (FILE *stream, gtopology *gt, const char */*path*/, const char *file, const char */*caption*/, double *valel, char flag)
{
  long i,j;
  double *valnod;
  // char *tcaption, *fff=NULL;
  
  //  // due to literals passed in caption
  //  tcaption = new char[strlen(caption)];
  //  i = 0;
  //  do  tcaption[i] = caption[i];
  //  while (caption[i++] != '\0');
  
  valnod = new double [gt->nn];
  
  //  //
  //  if (file == NULL){
  //    fff = new char [255];
  //    if (flag == 'e')  sprintf (fff, "%s%s.ex", path, tcaption);
  //    if (flag == 'd')  sprintf (fff, "%s%s.dx", path, tcaption);
  //    file = fff;
  //  }
  
  if (gt->nadjelnod == NULL)
    gt->adjacelem (stream);
  
  for (i=0; i<gt->nn; i++) {
    valnod[i] = 0.0;
    for (j=0; j<gt->nadjelnod[i]; j++)  valnod[i] += valel[gt->adjelnod[i][j]];
    valnod[i] /= (double)gt->nadjelnod[i];
  }
  
  if (flag == 'e')  print_ex (gt,file, valnod, valel);
  //if (flag == 'd')  print_dx (gt,file,&valnod,&valel,'n',&(i=1),&tcaption,1);
  
  //if (fff != NULL){ file = NULL; delete [] fff; }
  delete [] valnod;
  //delete [] tcaption;
}
