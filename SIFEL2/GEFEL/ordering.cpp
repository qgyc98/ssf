#include <stdlib.h>
#include "ordering.h"

/**
   function returns numbers of end points on 1D elements
   
   @param nodes - array of end point numbers
   
   JK, 12.10.2008
*/
void linbar_endpoints (long *nodes)
{
  nodes[0]=0;
  nodes[1]=1;
}

/**
   function returns numbers of end points on 1D elements
   
   @param nodes - array of end point numbers
   
   JK, 12.10.2008
*/
void quadbar_endpoints (long *nodes)
{
  nodes[0]=0;
  nodes[1]=1;
}

/**
   function returns numbers of nodes on edges
   of triangular elements with three nodes
   
   @param edgenod - array of edge node numbers
   @param edg - number of required edge
   
   JK, 24.8.2004
*/
void lintriangle_edgnod (long *edgenod,long edg)
{
  if (edg==0){
    edgenod[0]=0;
    edgenod[1]=1;
  }
  if (edg==1){
    edgenod[0]=1;
    edgenod[1]=2;
  }
  if (edg==2){
    edgenod[0]=2;
    edgenod[1]=0;
  }
}

/**
   function returns numbers of nodes on edges
   of triangular elements with six nodes
   
   @param edgenod - array of edge node numbers
   @param edg - number of required edge
   


*/ 
void quadtriangle_edgnod (long *edgenod,long edg)
{
  if (edg==0){
    edgenod[0]=0;  edgenod[1]=1;  edgenod[2]=3; 
  }
  if (edg==1){
    edgenod[0]=1;  edgenod[1]=2;  edgenod[2]=4; 
  }
  if (edg==2){ 
    edgenod[0]=2;  edgenod[1]=0; edgenod[2]=5; 
  }
}


/**
   function returns numbers of nodes on edges
   of quadrilateral elements with four nodes
   
   @param edgenod - array of edge node numbers
   @param edg - number of required edge
   
   JK, 24.8.2004
*/
void linquadrilat_edgnod (long *edgenod,long edg)
{
  if (edg==0){
    edgenod[0]=0;
    edgenod[1]=1;
  }
  if (edg==1){
    edgenod[0]=1;
    edgenod[1]=2;
  }
  if (edg==2){
    edgenod[0]=2;
    edgenod[1]=3;
  }
  if (edg==3){
    edgenod[0]=3;
    edgenod[1]=0;
  }
}



/**
   function returns numbers of nodes on edges
   of quadrilateral elements with eight nodes
   
   @param edgenod - array of edge node numbers
   @param edg - number of required edge
   
   JK, 24.8.2004
*/
void quadquadrilat_edgnod (long *edgenod,long edg)
{
  if (edg==0){
    edgenod[0]=0;
    edgenod[1]=1;
    edgenod[2]=4;
  }
  if (edg==1){
    edgenod[0]=1;
    edgenod[1]=2;
    edgenod[2]=5;
  }
  if (edg==2){
    edgenod[0]=2;
    edgenod[1]=3;
    edgenod[2]=6;
  }
  if (edg==3){
    edgenod[0]=3;
    edgenod[1]=0;
    edgenod[2]=7;
  }
}



/**
  The function returns numbers of nodes on edges
  of quadrilateral elements with twelve nodes.
   
  @param edgenod[out] - array of edge node numbers
  @param edg[in] - number of required edge
   
  @return The function returns edge node numbers in the argument edgenod.

  Created by Tomas Koudelka, 13.3.2019
*/
void cubicquadrilat_edgnod (long *edgenod, long edg)
{
  if (edg==0){
    edgenod[0]=0;
    edgenod[1]=1;
    edgenod[2]=4;
    edgenod[3]=5;
  }
  if (edg==1){
    edgenod[0]=1;
    edgenod[1]=2;
    edgenod[2]=6;
    edgenod[3]=7;
  }
  if (edg==2){
    edgenod[0]=2;
    edgenod[1]=3;
    edgenod[2]=8;
    edgenod[3]=9;
  }
  if (edg==3){
    edgenod[0]=3;
    edgenod[1]=0;
    edgenod[2]=10;
    edgenod[3]=11;
  }
}



/**
   function returns numbers of nodes on edges
   of quadrilateral elements with eight nodes
   
   @param edgenod - array of edge node numbers
   @param edg - number of required edge

*/
void lintetrahedral_edgnod(long *edgenod,long edg)
{
  if (edg==0){
    edgenod[0]=0;  edgenod[1]=1;
  }
  if (edg==1){
    edgenod[0]=1;  edgenod[1]=2; 
  }
  if (edg==2){ 
    edgenod[0]=2;  edgenod[1]=0; 
  }
  if (edg==3){
    edgenod[0]=1;  edgenod[1]=3; 
  }
  if (edg==4){
    edgenod[0]=3;  edgenod[1]=2; 
  }
  if (edg==5){
    edgenod[0]=3;  edgenod[1]=0; 
  }
  
}



/**
   Function returns numbers of nodes on edges
   of tetrahedron elements with ten nodes (quadratic approximation).
   
   @param edgenod - array of edge node numbers
   @param edg - number of required edge

   TKo 1.2009

*/
void quadtetrahedral_edgnod(long *edgenod,long edg)
{
  
  if (edg==0){
    edgenod[0]=0;  edgenod[1]=1;  edgenod[2]=4;
  }
  if (edg==1){
    edgenod[0]=1;  edgenod[1]=2;   edgenod[2]=5;
  }
  if (edg==2){ 
    edgenod[0]=2;  edgenod[1]=0;   edgenod[2]=6;
  }
  if (edg==3){
    edgenod[0]=1;  edgenod[1]=3;   edgenod[2]=8;
  }
  if (edg==4){
    edgenod[0]=3;  edgenod[1]=2;   edgenod[2]=9;
  }
  if (edg==5){
    edgenod[0]=3;  edgenod[1]=0;   edgenod[2]=7;
  }
  
}



/**
   Function returns numbers of nodes on edges
   of hexahedron elements with eight nodes (linear approximation).
   
   @param edgenod - array of edge node numbers
   @param edg - number of required edge
   
   JK
*/
void linhexahedral_edgnod(long *edgenod,long edg)
{
  
  if (edg==0){
    edgenod[0]=0;  edgenod[1]=1; 
  }
  if (edg==1){
    edgenod[0]=1;  edgenod[1]=2; 
  }
  if (edg==2){
    edgenod[0]=2;  edgenod[1]=3; 
  }
  if (edg==3){
    edgenod[0]=3;  edgenod[1]=0; 
  }  
  if (edg==4){
    edgenod[0]=0;  edgenod[1]=4; 
  }
  if (edg==5){
    edgenod[0]=1;  edgenod[1]=5; 
  }
  if (edg==6){
    edgenod[0]=2;  edgenod[1]=6; 
  }
  if (edg==7){
    edgenod[0]=3;  edgenod[1]=7; 
  }
  if (edg==8){
    edgenod[0]=4;  edgenod[1]=5; 
  }
  if (edg==9){
    edgenod[0]=5;  edgenod[1]=6; 
  }
  if (edg==10){ 
    edgenod[0]=6;  edgenod[1]=7; 
  }
  if (edg==11){
    edgenod[0]=7;  edgenod[1]=4; 
  }
  
}



/**
   Function returns numbers of nodes on edges
   of hexahedron elements with twenty nodes (quadratic approximation).
   
   @param edgenod - array of edge node numbers
   @param edg - number of required edge

   JK
*/
void quadhexahedral_edgnod (long *edgenod,long edg)
{
  if (edg==0){
    edgenod[0]=0;  edgenod[1]=1;  edgenod[2]= 8; 
  }
  if (edg==1){
    edgenod[0]=1;  edgenod[1]=2;  edgenod[2]= 9;
  }
  if (edg==2){
    edgenod[0]=2;  edgenod[1]=3;  edgenod[2]= 10;
  }
  if (edg==3){
    edgenod[0]=3;  edgenod[1]=0;  edgenod[2]= 11;
  }  
  if (edg==4){
    edgenod[0]=0;  edgenod[1]=4;  edgenod[2]= 12;
  }
  if (edg==5){
    edgenod[0]=1;  edgenod[1]=5;  edgenod[2]= 13;
  }
  if (edg==6){
    edgenod[0]=2;  edgenod[1]=6;  edgenod[2]= 14;
  }
  if (edg==7){
    edgenod[0]=3;  edgenod[1]=7;  edgenod[2]= 15;
  }
  if (edg==8){
    edgenod[0]=4;  edgenod[1]=5;  edgenod[2]= 16;
  }
  if (edg==9){
    edgenod[0]=5;  edgenod[1]=6;  edgenod[2]= 17;
  }
  if (edg==10){ 
    edgenod[0]=6;  edgenod[1]=7;  edgenod[2]= 18;
  }
  if (edg==11){
    edgenod[0]=7;  edgenod[1]=4;  edgenod[2]= 19;
  }
}

/**
   function returns numbers of nodes on surfaces
   of hexahedral elements with eight nodes
   
   @param surfnod - array of edge node numbers
   @param surf - number of required edge
   
   JK, 24.8.2004
*/
void linhexahedral_surfnod (long *surfnod,long surf)
{
  if (surf==0){
    surfnod[0]=0;
    surfnod[1]=3;
    surfnod[2]=7;
    surfnod[3]=4;
  }
  if (surf==1){
    surfnod[0]=1;
    surfnod[1]=0;
    surfnod[2]=4;
    surfnod[3]=5;
  }
  if (surf==2){
    surfnod[0]=2;
    surfnod[1]=1;
    surfnod[2]=5;
    surfnod[3]=6;
  }
  if (surf==3){
    surfnod[0]=3;
    surfnod[1]=2;
    surfnod[2]=6;
    surfnod[3]=7;
  }
  if (surf==4){
    surfnod[0]=0;
    surfnod[1]=1;
    surfnod[2]=2;
    surfnod[3]=3;
  }
  if (surf==5){
    surfnod[0]=4;
    surfnod[1]=5;
    surfnod[2]=6;
    surfnod[3]=7;
  }

}

/**
   function returns numbers of nodes on surfaces
   of hexahedral elements with twenty nodes
   
   @param surfnod - array of edge node numbers
   @param surf - number of required edge
   
   JK, 24.8.2004
*/
void quadhexahedral_surfnod (long *surfnod,long surf)
{
  if (surf==0){
    surfnod[0]=0;
    surfnod[1]=3;
    surfnod[2]=7;
    surfnod[3]=4;
    surfnod[4]=11;
    surfnod[5]=15;
    surfnod[6]=19;
    surfnod[7]=12;
  }
  if (surf==1){
    surfnod[0]=1;
    surfnod[1]=0;
    surfnod[2]=4;
    surfnod[3]=5;
    surfnod[4]=8;
    surfnod[5]=12;
    surfnod[6]=16;
    surfnod[7]=13;
  }
  if (surf==2){
    surfnod[0]=2;
    surfnod[1]=1;
    surfnod[2]=5;
    surfnod[3]=6;
    surfnod[4]=9;
    surfnod[5]=13;
    surfnod[6]=17;
    surfnod[7]=14;
  }
  if (surf==3){
    surfnod[0]=3;
    surfnod[1]=2;
    surfnod[2]=6;
    surfnod[3]=7;
    surfnod[4]=10;
    surfnod[5]=14;
    surfnod[6]=18;
    surfnod[7]=15;
  }
  if (surf==4){
    surfnod[0]=0;
    surfnod[1]=1;
    surfnod[2]=2;
    surfnod[3]=3;
    surfnod[4]=8;
    surfnod[5]=9;
    surfnod[6]=10;
    surfnod[7]=11;
  }
  if (surf==5){
    surfnod[0]=4;
    surfnod[1]=5;
    surfnod[2]=6;
    surfnod[3]=7;
    surfnod[4]=16;
    surfnod[5]=17;
    surfnod[6]=18;
    surfnod[7]=19;
  }

}



/**
   Function returns numbers of nodes on surfaces
   of tetrahedral elements with four nodes.
   
   @param surfnod - array of edge node numbers
   @param surf - number of required edge
   
   JK, 24.8.2004
*/
void lintetrahedral_surfnod (long *surfnod,long surf)
{
  if (surf==0){
    surfnod[0]=2;  surfnod[1]=1; surfnod[2]=3; 
  }
  if (surf==1){
    surfnod[0]=0;  surfnod[1]=2; surfnod[2]=3; 
  }
  if (surf==2){
    surfnod[0]=1;  surfnod[1]=0; surfnod[2]=3; 
  }
  if (surf==3){
    surfnod[0]=0;  surfnod[1]=1; surfnod[2]=2; 
  }
}



/**
   Function returns numbers of nodes on surfaces
   of tetrahedral elements with ten nodes.
   
   @param surfnod - array of edge node numbers
   @param surf - id of required surface
   
   TKo, 1.2009
*/
void quadtetrahedral_surfnod (long *surfnod,long surf)
{
  switch (surf)
  {
    case 0:
      surfnod[0]=2;  surfnod[1]=1; surfnod[2]=3;
      surfnod[3]=5;  surfnod[4]=8; surfnod[5]=9;  
      break;
    case 1:
      surfnod[0]=0;  surfnod[1]=2; surfnod[2]=3;
      surfnod[3]=6;  surfnod[4]=9; surfnod[5]=7;  
      break;
    case 2:
      surfnod[0]=1;  surfnod[1]=0; surfnod[2]=3; 
      surfnod[3]=4;  surfnod[4]=7; surfnod[5]=8;   
      break;
    case 3:
      surfnod[0]=0;  surfnod[1]=1; surfnod[2]=2; 
      surfnod[3]=4;  surfnod[4]=5; surfnod[5]=6;  
      break;
    default:
      print_err("unknown index of surface is required", __FILE__, __LINE__, __func__);
      abort();
  }
}


/**
   Function returns numbers of nodes on edges
   of wedge elements with six nodes (linear approximation).
   
   @param edgenod - array of edge node numbers
   @param edg - number of required edge
   
   JK, 03.2024
*/
void linwedge_edgnod(long *edgenod,long edg)
{
  
  if (edg==0){
    edgenod[0]=0;  edgenod[1]=1; 
  }
  if (edg==1){
    edgenod[0]=1;  edgenod[1]=2; 
  }
  if (edg==2){
    edgenod[0]=2;  edgenod[1]=0; 
  }
  if (edg==3){
    edgenod[0]=0;  edgenod[1]=3; 
  }  
  if (edg==4){
    edgenod[0]=1;  edgenod[1]=4; 
  }
  if (edg==5){
    edgenod[0]=2;  edgenod[1]=5; 
  }
  if (edg==6){
    edgenod[0]=3;  edgenod[1]=4; 
  }
  if (edg==7){
    edgenod[0]=4;  edgenod[1]=5; 
  }
  if (edg==8){
    edgenod[0]=5;  edgenod[1]=3; 
  }
}

/**
   function returns numbers of nodes on surfaces
   of linear wedge elements with six nodes
   
   @param surfnod - array of edge node numbers
   @param surf - number of required edge
   
   JK, 03.2024
*/
void linwedge_surfnod (long *surfnod,long surf)
{
  if (surf==0){
    surfnod[0]=2;
    surfnod[1]=1;
    surfnod[2]=4;
    surfnod[3]=5;
  }
  if (surf==1){
    surfnod[0]=0;
    surfnod[1]=2;
    surfnod[2]=5;
    surfnod[3]=3;
  }
  if (surf==2){
    surfnod[0]=1;
    surfnod[1]=0;
    surfnod[2]=3;
    surfnod[3]=4;
  }
  if (surf==3){
    surfnod[0]=0;
    surfnod[1]=1;
    surfnod[2]=2;
  }
  if (surf==4){
    surfnod[0]=4;
    surfnod[1]=3;
    surfnod[2]=5;
  }
}



/**
   function assembles natural coordinates of nodes
   on bar or beam with 2 nodes
   
   @param xi - natural coordinate of nodes
   
   JK, 29.11.2006
*/
void nodcoord_bar (vector &xi)
{
  xi[0] = -1.0;
  xi[1] =  1.0;
}

/**
   function assembles natural coordinates of nodes
   on bar or beam with 3 nodes
   
   @param xi - natural coordinate of nodes
   
   TKo, 23.5.2018
*/
void nodcoord_barq (vector &xi)
{
  xi[0] = -1.0;
  xi[1] =  1.0;
  xi[2] =  0.0;
}

/**
   function assembles natural coordinates of nodes
   on plane triangular element with 3 nodes
   
   @param xi,eta - natural coordinates of nodes
   
   JK, 27.11.2006
*/
void nodcoord_planelt (vector &xi,vector &eta)
{
  xi[0] =  1.0;  eta[0] =  0.0;
  xi[1] =  0.0;  eta[1] =  1.0;
  xi[2] =  0.0;  eta[2] =  0.0;
}


/**
   function assembles natural coordinates of nodes
   on plane quadrilateral element with 4 nodes (bilinear approximation functions)
   
   @param xi,eta - natural coordinates of nodes
   
   JK, 23.11.2006
*/
void nodcoord_planelq (vector &xi,vector &eta)
{
  xi[0] =  1.0;  eta[0] =  1.0;
  xi[1] = -1.0;  eta[1] =  1.0;
  xi[2] = -1.0;  eta[2] = -1.0;
  xi[3] =  1.0;  eta[3] = -1.0;
}



/**
   function assembles natural coordinates of nodes
   on plane quadrilateral element with 8 nodes (biquadratic approximation functions)
   
   @param xi,eta - natural coordinates of nodes

   JK, 23.11.2006
*/
void nodcoord_planeqq (vector &xi,vector &eta)
{
  xi[0] =  1.0;  eta[0] =  1.0;
  xi[1] = -1.0;  eta[1] =  1.0;
  xi[2] = -1.0;  eta[2] = -1.0;
  xi[3] =  1.0;  eta[3] = -1.0;
  xi[4] =  0.0;  eta[4] =  1.0;
  xi[5] = -1.0;  eta[5] =  0.0;
  xi[6] =  0.0;  eta[6] = -1.0;
  xi[7] =  1.0;  eta[7] =  0.0;
}



/**
  The function assembles natural coordinates of nodes
  on plane quadrilateral element with 12 nodes (bicubic approximation functions).
   
  @param xi[out],eta[out] - vector of natural coordinates of nodes
  
  @return The function returns coordinates in the arguments ksi and eta.

  Created by Tomas Koudelka, 13.3.2019
*/
void nodcoord_planecq (vector &xi,vector &eta)
{
  const double third = 1.0/3.0;
  xi[0] =  1.0;  eta[0] =  1.0;
  xi[1] = -1.0;  eta[1] =  1.0;
  xi[2] = -1.0;  eta[2] = -1.0;
  xi[3] =  1.0;  eta[3] = -1.0;

  xi[4] =   third;  eta[4] =  1.0;
  xi[5] =  -third;  eta[5] =  1.0;

  xi[6] =  -1.0;  eta[6] =  third;
  xi[7] =  -1.0;  eta[7] = -third;

  xi[8] =  -third;  eta[8] = -1.0;
  xi[9] =   third;  eta[9] = -1.0;

  xi[10] =  1.0;  eta[10] = -third;
  xi[11] =  1.0;  eta[11] =  third;
}



/**
   function assembles natural coordinates of nodes
   on tetrahedral element with 4 nodes
   
   @param xi - array containing natural coordinates xi
   @param eta - array containing natrual coordinates eta
   @param zeta - array containing natrual coordinates zeta
   
   JK, 2.10.2008
*/
void nodcoord_lintet (vector &xi,vector &eta,vector &zeta)
{
  xi[0] = 1.0;  eta[0] = 0.0;  zeta[0]=0.0;
  xi[1] = 0.0;  eta[1] = 1.0;  zeta[1]=0.0;
  xi[2] = 0.0;  eta[2] = 0.0;  zeta[2]=1.0;
  xi[3] = 0.0;  eta[3] = 0.0;  zeta[3]=0.0;
}

/**
   function assembles natural coordinates of nodes
   on tetrahedral element with 10 nodes
   
   @param xi - array containing natural coordinates xi
   @param eta - array containing natrual coordinates eta
   @param zeta - array containing natrual coordinates zeta
   
   TKo, 1.2009
*/
void nodcoord_quadtet (vector &xi,vector &eta,vector &zeta)
{
  xi[0] = 1.0;  eta[0] = 0.0;  zeta[0] = 0.0;
  xi[1] = 0.0;  eta[1] = 1.0;  zeta[1] = 0.0;
  xi[2] = 0.0;  eta[2] = 0.0;  zeta[2] = 1.0;
  xi[3] = 0.0;  eta[3] = 0.0;  zeta[3] = 0.0;
  xi[4] = 0.5;  eta[4] = 0.5;  zeta[4] = 0.0;
  xi[5] = 0.0;  eta[5] = 0.5;  zeta[5] = 0.5;
  xi[6] = 0.5;  eta[6] = 0.0;  zeta[6] = 0.0;
  xi[7] = 0.5;  eta[7] = 0.0;  zeta[7] = 0.0;
  xi[8] = 0.0;  eta[8] = 0.5;  zeta[8] = 0.0;
  xi[9] = 0.0;  eta[9] = 0.0;  zeta[9] = 0.5;
}

/**
   function assembles natural coordinates of nodes
   on hexahedral element with 8 nodes (trilinear approximation functions)
   
   @param xi,eta,zeta - natural coordinates of nodes

   JK, 23.11.2006
*/
void nodcoord_linhex (vector &xi,vector &eta,vector &zeta)
{
  xi[0] =  1.0;  eta[0] =  1.0;  zeta[0] =  1.0;
  xi[1] = -1.0;  eta[1] =  1.0;  zeta[1] =  1.0;
  xi[2] = -1.0;  eta[2] = -1.0;  zeta[2] =  1.0;
  xi[3] =  1.0;  eta[3] = -1.0;  zeta[3] =  1.0;
  xi[4] =  1.0;  eta[4] =  1.0;  zeta[4] = -1.0;
  xi[5] = -1.0;  eta[5] =  1.0;  zeta[5] = -1.0;
  xi[6] = -1.0;  eta[6] = -1.0;  zeta[6] = -1.0;
  xi[7] =  1.0;  eta[7] = -1.0;  zeta[7] = -1.0;
}

/**
   function assembles natural coordinates of nodes
   on hexahedral element with 20 nodes (triquadratic approximation functions)
   
   @param xi,eta,zeta - natural coordinates of nodes

   JK, 23.11.2006
*/
void nodcoord_quadhex (vector &xi,vector &eta,vector &zeta)
{
  xi[0]  =  1.0;  eta[0]  =  1.0;  zeta[0]  =  1.0;
  xi[1]  = -1.0;  eta[1]  =  1.0;  zeta[1]  =  1.0;
  xi[2]  = -1.0;  eta[2]  = -1.0;  zeta[2]  =  1.0;
  xi[3]  =  1.0;  eta[3]  = -1.0;  zeta[3]  =  1.0;
  xi[4]  =  1.0;  eta[4]  =  1.0;  zeta[4]  = -1.0;
  xi[5]  = -1.0;  eta[5]  =  1.0;  zeta[5]  = -1.0;
  xi[6]  = -1.0;  eta[6]  = -1.0;  zeta[6]  = -1.0;
  xi[7]  =  1.0;  eta[7]  = -1.0;  zeta[7]  = -1.0;
  xi[8]  =  1.0;  eta[8]  =  1.0;  zeta[8]  =  1.0;
  xi[9]  = -1.0;  eta[9]  =  1.0;  zeta[9]  =  1.0;
  xi[10] = -1.0;  eta[10] = -1.0;  zeta[10] =  1.0;
  xi[11] =  1.0;  eta[11] = -1.0;  zeta[11] =  1.0;
  xi[12] =  1.0;  eta[12] =  1.0;  zeta[12] = -1.0;
  xi[13] = -1.0;  eta[13] =  1.0;  zeta[13] = -1.0;
  xi[14] = -1.0;  eta[14] = -1.0;  zeta[14] = -1.0;
  xi[15] =  1.0;  eta[15] = -1.0;  zeta[15] = -1.0;
  xi[16] =  1.0;  eta[16] =  1.0;  zeta[16] =  1.0;
  xi[17] = -1.0;  eta[17] =  1.0;  zeta[17] =  1.0;
  xi[18] = -1.0;  eta[18] = -1.0;  zeta[18] =  1.0;
  xi[19] =  1.0;  eta[19] = -1.0;  zeta[19] =  1.0;
}

/**
   function assembles natural coordinates of nodes
   on linear wedge with 6 nodes
   
   @param xi,eta,zeta - natural coordinates of nodes

   JK, 03.2024
*/
void nodcoord_linwedge (vector &xi,vector &eta,vector &zeta)
{
  xi[0] =  1.0;  eta[0] =  0.0;  zeta[0] =  1.0;
  xi[1] =  0.0;  eta[1] =  1.0;  zeta[1] =  1.0;
  xi[2] =  0.0;  eta[2] =  0.0;  zeta[2] =  1.0;
  xi[3] =  1.0;  eta[3] =  0.0;  zeta[3] = -1.0;
  xi[4] =  0.0;  eta[4] =  1.0;  zeta[4] = -1.0;
  xi[5] =  0.0;  eta[5] =  0.0;  zeta[5] = -1.0;
}

/**
   function assembles map between nodes and the closest integration points
   on bar or beam elements, elements have two nodes
   
   @param n - order of numerical integration
   @param ipnum - array of map between nodes and integration points
   
   ipnum[i]=j - the closest integration point to the i-th node is the j-th integration point
   
   JK, 29.11.2006
*/
void nodip_bar (long i,long n,ivector &ipnum)
{
  ipnum[0]=i;
  ipnum[1]=i+n-1;
}

/**
   function assembles map between nodes and the closest integration points
   on bar or beam elements, elements have two nodes
   
   @param n - order of numerical integration
   @param ipnum - array of map between nodes and integration points
   
   ipnum[i]=j - the closest integration point to the i-th node is the j-th integration point
   
   TKo, 14.5.2018
*/
void nodip_barq (long i,long n,ivector &ipnum)
{
  ipnum[0]=i;
  ipnum[1]=i+n-1;
  ipnum[2]=i+1;
}

/**
   function assembles map between nodes and the closest integration points
   on plane quadrilateral element with 4 nodes
   
   @param n - order of numerical integration
   @param ipnum - array of map between nodes and integration points
   
   ipnum[i]=j - the closest integration point to the i-th node is the j-th integration point
   
   JK, 23.11.2006
*/
void nodip_planelq (long i,long n,ivector &ipnum)
{
  ipnum[0]=i+n*(n-1)+n-1;
  ipnum[1]=i+n-1;
  ipnum[2]=i+0;
  ipnum[3]=i+n*(n-1);
}

/**
   function assembles map between nodes and the closest integration points
   on plane quadrilateral element with 4 nodes
   
   @param i - integration point id
   @param n - order of numerical integration
   @param ipnum - array of map between nodes and integration points
   
   ipnum[i]=j - the closest integration point to the i-th node is the j-th integration point
   
   JK, 23.11.2006
*/
void nodip_planelt (long i,long n,ivector &ipnum)
{
  switch (n){
  case 1:{
    ipnum[0]=i;
    ipnum[1]=i;
    ipnum[2]=i;
    break;
  }
  case 3:{
    ipnum[0]=i;
    ipnum[1]=i+1;
    ipnum[2]=i+2;
    break;
  }
  default:{
    print_err("unknown order of numerical integration is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function assembles map between nodes and the closest integration points
   on space tetrahedral element with 4 nodes
   
   @param i - integration point id
   @param ipnum - array of map between nodes and integration points
   
   ipnum[i]=j - the closest integration point to the i-th node is the j-th integration point
   
   TKo, 11.2008
*/
void nodip_lintet (long i,ivector &ipnum)
{
  ipnum[0]=i;
  ipnum[1]=i;
  ipnum[2]=i;
  ipnum[3]=i;
}

/**
   function assembles map between nodes and the closest integration points
   on space tetrahedral element with 10 nodes
   
   @param i - integration point id
   @param ipnum - array of map between nodes and integration points
   
   ipnum[i]=j - the closest integration point to the i-th node is the j-th integration point
   
   TKo, 11.2008
*/
void nodip_quadtet (long i,ivector &ipnum)
{
  ipnum[0]=i+0;
  ipnum[1]=i+1;
  ipnum[2]=i+2;
  ipnum[3]=i+3;
  ipnum[4]=i+0;
  ipnum[5]=i+1;
  ipnum[6]=i+2;
  ipnum[7]=i+0;
  ipnum[8]=i+1;
  ipnum[9]=i+2;
}


/**
   function assembles map between nodes and the closest integration points
   on space hexahedral element with 8 nodes
   
   @param n - order of numerical integration
   @param ipnum - array of map between nodes and integration points
   
   ipnum[i]=j - the closest integration point to the i-th node is the j-th integration point
   
   JK, 23.11.2006
*/
void nodip_linhex (long i,long n,ivector &ipnum)
{
  ipnum[0]=i+n*n*n-1;
  ipnum[1]=i+n*n-1;
  ipnum[2]=i+n-1;
  ipnum[3]=i+n*n*(n-1)+n-1;
  ipnum[4]=i+n*n*(n-1)+n*(n-1);
  ipnum[5]=i+n*(n-1);
  ipnum[6]=i+0;
  ipnum[7]=i+n*n*(n-1);
}

/**
   function assembles map between nodes and the closest integration points
   on space hexahedral element with 20 nodes
   
   @param n - order of numerical integration
   @param ipnum - array of map between nodes and integration points
   
   ipnum[i]=j - the closest integration point to the i-th node is the j-th integration point
   
   JK, 27.11.2006
*/
void nodip_quadhex (long i,long n,ivector &ipnum)
{
  ipnum[0]=i+n*n*n-1;
  ipnum[1]=i+n*n-1;
  ipnum[2]=i+n-1;
  ipnum[3]=i+n*n*(n-1)+n-1;

  ipnum[4]=i+n*n*(n-1)+n*(n-1);
  ipnum[5]=i+n*(n-1);
  ipnum[6]=i;
  ipnum[7]=i+n*n*(n-1);

  ipnum[8]=i+n*n+n*n-1;
  ipnum[9]=i+n+n-1;
  ipnum[10]=i+n*n+n-1;
  ipnum[11]=i+n*n*(n-1)+n+n-1;

  ipnum[12]=i+n*n*(n-1)+n*(n-1)+1;
  ipnum[13]=i+n*(n-1)+1;
  ipnum[14]=i+1;
  ipnum[15]=i+n*n*(n-1)+1;

  ipnum[16]=i+n*n+n*(n-1);
  ipnum[17]=i+n;
  ipnum[18]=i+n*n;
  ipnum[19]=i+n*n*(n-1)+n;
}
