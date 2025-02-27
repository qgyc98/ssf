#include "elemparticle.h"
#include "global.h"
#include "mechmat.h"
#include "mechtop.h"
#include "mechcrsec.h"
#include "globmat.h"
#include "genfile.h"
#include "element.h"
#include "node.h"
#include "intpoints.h"
#include "loadcase.h"
#include <math.h>


elemparticle::elemparticle (long cnne,long cdim)
{
  long i;
  
  //  dimension of solved problem
  dim = cdim;

  //  number nodes (particles) on element
  //  it is read from input file
  nne=cnne;
  //nne=5;
  
  //  number of DOFs
  ndofe = nne*cdim;
  //ndofe=5;
  
  //  number of blocks
  //  there are no blocks, number of blocks is set to one
  //  with respect of element philosophy of the code
  nb=1;
  
  //  number of integration points
  //  one interatomic potential is assumed in one element/cell
  //  more general case can be modelled after definition of several
  //  integration points and several potentials
  nip = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
  }
  
  nip[0][0]=1;
  
  tnip=0;
  for (i=0;i<nb;i++){
    tnip+=nip[i][i];
  }
  
}

elemparticle::~elemparticle (void)
{
  long i;
  
  for (i=0;i<nb;i++){
    delete [] nip[i];
  }
  delete [] nip;

}





/**
   function computes direction vectors

   @param eid - element id
   @param i,j - particle ids
   @param s - direction %vector
   @param x - %vector of coordinates
   @param u - %vector of displacements
   
   JK, 19.6.2005
*/
void elemparticle::direction_vector_1d (long /*eid*/,long i,long j,vector &s,vector &x,vector &u)
{
  s[0]=x[j]+u[j] - x[i]-u[i];
}

/**
   function computes direction vectors

   @param eid - element id
   @param i,j - particle ids
   @param s - direction %vector
   @param x,y - vectors of coordinates
   @param u,v - vectors of displacements
   
   JK, 19.6.2005
*/
void elemparticle::direction_vector_2d (long /*eid*/,long i,long j,vector &s,
					vector &x,vector &y,vector &u,vector &v)
{
  s[0]=x[j]+u[j] - x[i]-u[i];
  s[1]=y[j]+v[j] - y[i]-v[i];
}

/**
   function computes direction vectors

   @param eid - element id
   @param i,j - particle ids
   @param s - direction %vector
   @param x,y,z - vectors of coordinates
   @param u,v,w - vectors of displacements
   
   JK, 19.6.2005
*/
void elemparticle::direction_vector_3d (long /*eid*/,long i,long j,vector &s,
					vector &x,vector &y,vector &z,
					vector &u,vector &v,vector &w)
{
  s[0]=x[j]+u[j] - x[i]-u[i];
  s[1]=y[j]+v[j] - y[i]-v[i];
  s[2]=z[j]+w[j] - z[i]-w[i];
}

/**
   function computes stiffness in 1D
   
   @param ipp - integration point id
   @param k - stiffness %matrix (diagonal submatrices)
   @param s - direction %vector
   
   JK, 23.6.2005
*/
void elemparticle::stiffmat_1d_kii (long ipp,matrix &k,vector &s)
{
  double f,g,r;
  
  //  norm of the direction vector
  r = normv(s);
  
  //  first derivative of particle potential with respect to particle distance
  f = Mm->give_first_derivative (ipp,r);
  
  //  second derivative of particle potential with respect to particle distance
  g = Mm->give_second_derivative (ipp,r);
  
  //  part of the stiffness matrix
  k[0][0] = g*s[0]*s[0]/r/r + f*(1.0/r-s[0]*s[0]/r/r/r);
}

/**
   function computes stiffness in 1D
   
   @param ipp - integration point id
   @param k - stiffness %matrix (offdiagonal submatrices)
   @param s - direction %vector
   
   JK, 23.6.2005
*/
void elemparticle::stiffmat_1d_kij (long ipp,matrix &k,vector &s)
{
  double f,g,r;
  
  //  norm of the direction vector
  r = normv(s);
  
  //  first derivative of particle potential with respect to particle distance
  f = Mm->give_first_derivative (ipp,r);
  
  //  second derivative of particle potential with respect to particle distance
  g = Mm->give_second_derivative (ipp,r);
  
  //  part of the stiffness matrix
  k[0][0] = -1.0*g*s[0]*s[0]/r/r + f*(-1.0/r+s[0]*s[0]/r/r/r);
}

/**
   function computes stiffness %matrix of cell of particles
   
   @param eid - element id
   @param sm - stiffness %matrix
   
   JK, 23.6.2005
*/
void elemparticle::stiffness_matrix_1d (long eid,matrix &sm)
{
  long i,j,ipp;
  ivector nodes(nne);
  vector x(nne),u(nne),s(1);
  matrix k(1,1);
  
  //  node numbers
  Mt->give_elemnodes (eid,nodes);
  //  node coordinates
  Mt->give_node_coord1d (x,eid);
  //  node displacements
  //if (Mp->phase==1){
    Mt->give_noddispl_1d (nodes,u);
    //}
  
  ipp = Mt->elements[eid].ipp[0][0];
  
  fillm (0.0, sm);
  
  for (i=0;i<nne;i++){
    for (j=i+1;j<nne;j++){
      direction_vector_1d (eid,i,j,s,x,u);
      
      stiffmat_1d_kii (ipp,k,s);
      
      //if (k[0][0]<1.0e2)  k[0][0]=1.0e2;
      
      sm[i][i]+=k[0][0];
      sm[j][j]+=k[0][0];

      stiffmat_1d_kij (ipp,k,s);
      sm[i][j]+=k[0][0];
      sm[j][i]+=k[0][0];
    }
  }


  //fprintf (Out,"\nsm [0][0]  %le\n",sm[0][0]);
  
  /*
  if (Mp->phase==1){
    for (i=0;i<nne;i++){
      for (j=0;j<nne;j++){
	sm[i][j]*=-1.0;
      }
      if (fabs(sm[i][i])<1.0e-5){
	if (sm[i][i]<0.0)
	  sm[i][i]=-1.0e-5;
	if (sm[i][i]>0.0)
	  sm[i][i]=1.0e-5;
      }
    }
    
  }
  */

  fprintf (Out,"\n\n STIFFNESS MATRIX");
  for (i=0;i<nne;i++){
    fprintf (Out,"\n");
    for (j=0;j<nne;j++){
      fprintf (Out,"  %le",sm[i][j]);
    }
  }
  fprintf (Out,"\n");
  
  if (sm[0][0]<1.0e-2){
    sm[0][0]+=1.0;
    sm[0][1]-=1.0;
    sm[1][0]-=1.0;
    sm[1][1]+=1.0;
  }

}






/**
   function computes stiffness in 2D
   
   @param ipp - integration point id
   @param k - stiffness %matrix (diagonal submatrices)
   @param s - direction %vector
   
   JK, 26.9.2005
*/
void elemparticle::stiffmat_2d_kii (long ipp,matrix &k,vector &s)
{
  double f,g,r;
  
  //  norm of the direction vector
  r = normv(s);
  
  //  first derivative of particle potential with respect to particle distance
  f = Mm->give_first_derivative (ipp,r);
  
  //  second derivative of particle potential with respect to particle distance
  g = Mm->give_second_derivative (ipp,r);
  
  //  part of the stiffness matrix
  k[0][0] = g*s[0]*s[0]/r/r + f*(1.0/r-s[0]*s[0]/r/r/r);
  k[0][1] = g*s[0]*s[1]/r/r - f*s[0]*s[1]/r/r/r;

  k[1][0] = k[0][1];
  k[1][1] = g*s[1]*s[1]/r/r + f*(1.0/r-s[1]*s[1]/r/r/r);
}

/**
   function computes stiffness in 2D
   
   @param ipp - integration point id
   @param k - stiffness %matrix (offdiagonal submatrices)
   @param s - direction %vector
   
   JK, 26.9.2005
*/
void elemparticle::stiffmat_2d_kij (long ipp,matrix &k,vector &s)
{
  double f,g,r;
  
  //  norm of the direction vector
  r = normv(s);
  
  //  first derivative of particle potential with respect to particle distance
  f = Mm->give_first_derivative (ipp,r);
  
  //  second derivative of particle potential with respect to particle distance
  g = Mm->give_second_derivative (ipp,r);
  
  //  part of the stiffness matrix
  k[0][0] = -1.0*g*s[0]*s[0]/r/r + f*(-1.0/r+s[0]*s[0]/r/r/r);
  k[0][1] = -1.0*g*s[0]*s[1]/r/r + f*s[0]*s[1]/r/r/r;

  k[1][0] = k[0][1];
  k[1][1] = -1.0*g*s[1]*s[1]/r/r + f*(-1.0/r+s[1]*s[1]/r/r/r);
}

/**
   function computes stiffness %matrix of cell of particles
   
   @param eid - element id
   @param sm - stiffness %matrix
   
   JK, 26.9.2005
*/
void elemparticle::stiffness_matrix_2d (long eid,matrix &sm)
{
  long i,j,ipp;
  ivector nodes(nne);
  vector x(nne),y(nne),u(nne),v(nne),s(2);
  matrix k(2,2);
  
  //  node numbers
  Mt->give_elemnodes (eid,nodes);
  //  node coordinates
  Mt->give_node_coord2d (x,y,eid);
  //  node displacements
  //if (Mp->phase==1){
  Mt->give_noddispl_2d (nodes,u,v);
  //}
  
  ipp = Mt->elements[eid].ipp[0][0];
  
  fillm (0.0,sm);
  
  for (i=0;i<nne;i++){
    for (j=i+1;j<nne;j++){
      direction_vector_2d (eid,i,j,s,x,y,u,v);
      
      stiffmat_2d_kii (ipp,k,s);
      
      //if (k[0][0]<1.0e2)  k[0][0]=1.0e2;
      
      sm[i*2+0][i*2+0]+=k[0][0];
      sm[i*2+0][i*2+1]+=k[0][1];
      sm[i*2+1][i*2+0]+=k[1][0];
      sm[i*2+1][i*2+1]+=k[1][1];

      sm[j*2+0][j*2+0]+=k[0][0];
      sm[j*2+0][j*2+1]+=k[0][1];
      sm[j*2+1][j*2+0]+=k[1][0];
      sm[j*2+1][j*2+1]+=k[1][1];

      stiffmat_2d_kij (ipp,k,s);

      sm[i*2+0][j*2+0]+=k[0][0];
      sm[i*2+0][j*2+1]+=k[0][1];
      sm[i*2+1][j*2+0]+=k[1][0];
      sm[i*2+1][j*2+1]+=k[1][1];

      sm[j*2+0][i*2+0]+=k[0][0];
      sm[j*2+0][i*2+1]+=k[0][1];
      sm[j*2+1][i*2+0]+=k[1][0];
      sm[j*2+1][i*2+1]+=k[1][1];
    }
  }


  //fprintf (Out,"\nsm [0][0]  %le\n",sm[0][0]);
  
  /*
  if (Mp->phase==1){
    for (i=0;i<nne;i++){
      for (j=0;j<nne;j++){
	sm[i][j]*=-1.0;
      }
      if (fabs(sm[i][i])<1.0e-5){
	if (sm[i][i]<0.0)
	  sm[i][i]=-1.0e-5;
	if (sm[i][i]>0.0)
	  sm[i][i]=1.0e-5;
      }
    }
    
  }
  */

  fprintf (Out,"\n\n STIFFNESS MATRIX");
  for (i=0;i<nne;i++){
    fprintf (Out,"\n");
    for (j=0;j<nne;j++){
      fprintf (Out,"  %le",sm[i][j]);
    }
  }
  fprintf (Out,"\n");
    
}




/**
   function computes diagonal stiffness matrix block in 3D
   
   @param ipp - integration point id
   @param k - stiffness %matrix (diagonal submatrices)
   @param s - direction %vector
   
   JK, 26.9.2005
*/
void elemparticle::stiffmat_3d_kii (long ipp,matrix &k,vector &s)
{
  double f,g,r;
  
  //  norm of the direction vector
  r = normv(s);
  
  //  first derivative of particle potential with respect to particle distance
  f = Mm->give_first_derivative (ipp,r);
  
  //  second derivative of particle potential with respect to particle distance
  g = Mm->give_second_derivative (ipp,r);
  
  //  part of the stiffness matrix
  k[0][0] = g*s[0]*s[0]/r/r + f*(1.0/r-s[0]*s[0]/r/r/r);
  k[0][1] = g*s[0]*s[1]/r/r - f*s[0]*s[1]/r/r/r;
  k[0][2] = g*s[0]*s[2]/r/r - f*s[0]*s[2]/r/r/r;

  k[1][0] = k[0][1];
  k[1][1] = g*s[1]*s[1]/r/r + f*(1.0/r-s[1]*s[1]/r/r/r);
  k[1][2] = g*s[1]*s[2]/r/r - f*s[1]*s[2]/r/r/r;

  k[2][0] = k[0][2];
  k[2][1] = k[1][2];
  k[2][2] = g*s[2]*s[2]/r/r + f*(1.0/r-s[2]*s[2]/r/r/r);
}

/**
   function computes off-diagonal stiffness matrix block in 3D
   
   @param ipp - integration point id
   @param k - stiffness %matrix (offdiagonal submatrices)
   @param s - direction %vector
   
   JK, 26.9.2005
*/
void elemparticle::stiffmat_3d_kij (long ipp,matrix &k,vector &s)
{
  double f,g,r;
  
  //  norm of the direction vector
  r = normv (s);
  
  //  first derivative of particle potential with respect to particle distance
  f = Mm->give_first_derivative (ipp,r);
  
  //  second derivative of particle potential with respect to particle distance
  g = Mm->give_second_derivative (ipp,r);
  
  //  part of the stiffness matrix
  k[0][0] = -1.0*g*s[0]*s[0]/r/r + f*(-1.0/r+s[0]*s[0]/r/r/r);
  k[0][1] = -1.0*g*s[0]*s[1]/r/r + f*s[0]*s[1]/r/r/r;
  k[0][2] = -1.0*g*s[0]*s[2]/r/r + f*s[0]*s[2]/r/r/r;

  k[1][0] = k[0][1];
  k[1][1] = -1.0*g*s[1]*s[1]/r/r + f*(-1.0/r+s[1]*s[1]/r/r/r);
  k[1][2] = -1.0*g*s[1]*s[2]/r/r + f*s[1]*s[2]/r/r/r;

  k[2][0] = k[0][2];
  k[2][1] = k[1][2];
  k[2][2] = -1.0*g*s[2]*s[2]/r/r + f*(-1.0/r+s[2]*s[2]/r/r/r);
}

/**
   function computes stiffness %matrix of cell of particles in 3D
   
   @param eid - element id
   @param sm - stiffness %matrix
   
   JK, 26.9.2005
*/
void elemparticle::stiffness_matrix_3d (long eid,matrix &sm)
{
  long i,j,ipp,ii,jj;
  ivector nodes(nne);
  vector x(nne),y(nne),u(nne),v(nne),s(2);
  matrix k(3,3);
  
  //  node numbers
  Mt->give_elemnodes (eid,nodes);
  //  node coordinates
  Mt->give_node_coord2d (x,y,eid);
  //  node displacements
  //if (Mp->phase==1){
  Mt->give_noddispl_2d (nodes,u,v);
  //}
  
  ipp = Mt->elements[eid].ipp[0][0];
  
  fillm (0.0,sm);
  
  for (i=0;i<nne;i++){
    for (j=i+1;j<nne;j++){
      direction_vector_2d (eid,i,j,s,x,y,u,v);
      
      stiffmat_2d_kii (ipp,k,s);
      for (ii=0; ii<3; ii++)
      {
        for(jj=0; jj<3; jj++)
	{
	  sm[i*3+ii][i*3+jj] += k[ii][jj];
	  sm[j*3+ii][j*3+jj] += k[ii][jj];
	}
      }

      stiffmat_2d_kij (ipp,k,s);
      for (ii=0; ii<3; ii++)
      {
        for(jj=0; jj<3; jj++)
	{
	  sm[i*3+ii][j*3+jj] += k[ii][jj];
	  sm[j*3+ii][i*3+jj] += k[ii][jj];
	}
      }
    }
  }
  fprintf (Out,"\n\n STIFFNESS MATRIX");
  for (i=0;i<nne;i++){
    fprintf (Out,"\n");
    for (j=0;j<nne;j++){
      fprintf (Out,"  %le",sm[i][j]);
    }
  }
  fprintf (Out,"\n");
    
}




/**
   

*/
void elemparticle::res_stiffness_matrix (long eid,matrix &sm)
{
  switch (dim){
  case 1:{
    stiffness_matrix_1d (eid,sm);
    break;
  }
  case 2:{
    stiffness_matrix_2d (eid,sm);
    break;
  }
  case 3:{
    stiffness_matrix_3d (eid,sm);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown dimension of problem is required in function res_stiffness_matrix (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


void elemparticle::forces_1d (long ipp,vector &fij,vector &s)
{
  double r,d;
  
  //  norm of the direction vector
  r = normv(s);
  
  //  first derivative of particle potential with respect to particle distance
  d = Mm->give_first_derivative (ipp,r);
  
  fij[0] = -s[0]*d/r;
}

void elemparticle::inter_forces_1d (long eid,vector &f)
{
  long i,j,ipp;
  ivector nodes(nne);
  vector x(nne),u(nne),s(1),fij(1);
  
  //  node numbers
  Mt->give_elemnodes (eid,nodes);
  //  node coordinates
  Mt->give_node_coord1d (x,eid);
  //  node displacements
  Mt->give_noddispl_1d (nodes,u);

  ipp = Mt->elements[eid].ipp[0][0];
  
  fprintf (Out,"\n u   %le   %le\n",u[0],u[1]);

  fillv (0.0,f);

  for (i=0;i<nne;i++){
    for (j=i;j<nne;j++){
      if (i != j){
	direction_vector_1d (eid,i,j,s,x,u);
	forces_1d (ipp,fij,s);
	
	fprintf (Out,"\n s  %le   fij  %le\n",s[0],fij[0]);

	f[i]+=fij[0];
	f[j]-=fij[0];
      }
    }
  }
  copyv(f.a, Mm->ip[ipp].stress, f.n);
}


void elemparticle::forces_2d (long ipp,vector &fij,vector &s)
{
  double r,d;
  
  //  norm of the direction vector
  r = normv(s);
  //  first derivative of particle potential with respect to particle distance
  d = Mm->give_first_derivative (ipp,r);
  
  fij[0] = -s[0]*d/r;
  fij[1] = -s[1]*d/r;
}

void elemparticle::inter_forces_2d (long eid,vector &f)
{
  long i,j,ipp,ii;
  ivector nodes(nne);
  vector x(nne),y(nne),u(nne),v(nne),s(2),fij(2);
  
  //  node numbers
  Mt->give_elemnodes (eid,nodes);
  //  node coordinates
  Mt->give_node_coord2d (x,y,eid);
  //  node displacements
  Mt->give_noddispl_2d (nodes,u, v);
  // integration point
  ipp = Mt->elements[eid].ipp[0][0];
  
  fillv (0.0,f);
  for (i=0;i<nne;i++)
  {
    for (j=i+1;j<nne;j++)
    {
      direction_vector_2d (eid,i,j,s,x,y,u,v);
      forces_2d (ipp,fij,s);
      for(ii=0; ii<2; ii++) // loop over directions of forces
      {
        f[i*2+ii]+=fij[ii];
        f[j*2+ii]-=fij[ii];
      }
    }
  }
}


void elemparticle::forces_3d (long ipp,vector &fij,vector &s)
{
  double r,d;
  
  //  norm of the direction vector
  r = normv(s);
  //  first derivative of particle potential with respect to particle distance
  d = Mm->give_first_derivative (ipp,r);
  
  fij[0] = -s[0]*d/r;
  fij[1] = -s[1]*d/r;
  fij[2] = -s[2]*d/r;
}

void elemparticle::inter_forces_3d (long eid,vector &f)
{
  long i,j,ii,ipp;
  ivector nodes(nne);
  vector x(nne),y(nne),z(nne),u(nne),v(nne),w(nne),s(3),fij(3);
  
  //  node numbers
  Mt->give_elemnodes (eid,nodes);
  //  node coordinates
  Mt->give_node_coord3d (x,y,z,eid);
  //  node displacements
  Mt->give_noddispl_3d (nodes,u,v,w);
  // integration point
  ipp = Mt->elements[eid].ipp[0][0];
  
  fillv (0.0,f);
  for (i=0;i<nne;i++)
  {
    for (j=i+1;j<nne;j++)
    {
      direction_vector_3d (eid,i,j,s,x,y,z,u,v,w);
      forces_3d (ipp,fij,s);
      for(ii=0; ii<3; ii++) // loop over directions of forces
      {
        f[i*3+ii]+=fij[ii];
        f[j*3+ii]-=fij[ii];
      }
    }
  }
}


void elemparticle::res_internal_forces (long eid,vector &ifor)
{
  switch (dim){
  case 1:{
    inter_forces_1d (eid,ifor);
    break;
  }
  case 2:{
    inter_forces_2d (eid,ifor);
    break;
  }
  case 3:{
    inter_forces_3d (eid,ifor);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown dimension of problem is required in function res_stiffness_matrix (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }
}
