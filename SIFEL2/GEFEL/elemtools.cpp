#include <stdlib.h>
#include <math.h>
#include "genfile.h"
#include "elemtools.h"
#include "matrix.h"

/**
   function generates outer normal %vector to required element end
   
   @param enid - id of end node
   @param n - %vector of the outer normal
   
   JK, 8.3.2018
*/
void bar_normal_vectors (long enid,vector &n)
{
  if (enid==0){
    n[0]=-1.0;
  }
  if (enid==1){
    n[0]=1.0;
  }
}

/**
   function generates outer normal %vector to required element edge
   
   @param edgeid - id of edge
   @param x,y - arrays of node coordinates
   @param n - %vector of the outer normal
   
   JK, 3.3.2018
*/
void triangle_normal_vectors (long edgeid,vector &x,vector &y,vector &n)
{
  double norm,ss,alpha;
  vector t1(2),t2(2);
  
  if (edgeid==0){
    // direction vector of the required edge
    t1[0]=x[1]-x[0];
    t1[1]=y[1]-y[0];
    //  auxiliary vector
    t2[0]=x[0]-x[2];
    t2[1]=y[0]-y[2];
  }
  if (edgeid==1){
    // direction vector of the required edge
    t1[0]=x[2]-x[1];
    t1[1]=y[2]-y[1];
    //  auxiliary vector
    t2[0]=x[1]-x[0];
    t2[1]=y[1]-y[0];
  }
  if (edgeid==2){
    // direction vector of the required edge
    t1[0]=x[0]-x[2];
    t1[1]=y[0]-y[2];
    //  auxiliary vector
    t2[0]=x[2]-x[1];
    t2[1]=y[2]-y[1];
  }
  
  //  length of the direction vector
  norm = t1[0]*t1[0]+t1[1]*t1[1];
  if (norm<1.0e-10){
    print_err("nonpositive length of the direction vector",__FILE__,__LINE__,__func__);
    abort ();
  }

  scprd (t1,t2,ss);
  alpha = (0.0-ss)/norm;
  
  n[0]=t2[0]-alpha*t1[0];
  n[1]=t2[1]-alpha*t1[1];

  norm = sqrt(n[0]*n[0]+n[1]*n[1]);
  n[0]/=norm;
  n[1]/=norm;

}


/**
   function generates outer normal %vector to required linear quadrilateral element edge
   
   @param edgeid - id of edge
   @param x,y - arrays of node coordinates
   @param n - %vector of the outer normal
   
   TKr, 29/11/2019
*/
void quadlin_normal_vectors (long edgeid,vector &x,vector &y,vector &n)
{
  double norm;
  vector t1(3);
  
  if (edgeid==0){
    // direction vector of the required edge
    t1[0]=x[1]-x[0];
    t1[1]=y[1]-y[0];
  }
  if (edgeid==1){
    // direction vector of the required edge
    t1[0]=x[2]-x[1];
    t1[1]=y[2]-y[1];
  }
  if (edgeid==2){
    // direction vector of the required edge
    t1[0]=x[3]-x[2];
    t1[1]=y[3]-y[2];
  }
  if (edgeid==3){
    // direction vector of the required edge
    t1[0]=x[0]-x[3];
    t1[1]=y[0]-y[3];
  }
  
  //  length of the direction vector
  norm = t1[0]*t1[0]+t1[1]*t1[1];
  if (norm<1.0e-10){
    print_err("nonpositive length of the direction vector",__FILE__,__LINE__,__func__);
    abort ();
  }

  norm = sqrt(norm);

  n[0] = t1[1];
  n[1] = -t1[0];
  
  n[0]/=norm;
  n[1]/=norm;
}

/**
   function computes length of linear quadrilateral element egdes
   
   @param surfid - id of surface
   @param x,y - arrays of node coordinates
   
   TKr, 29/11/2019
*/
double quadlin_edge_length (long edgeid,vector &x,vector &y)
{
  double norm;
  vector t1(3);

  if (edgeid==0){
    // direction vector of the required edge
    t1[0]=x[1]-x[0];
    t1[1]=y[1]-y[0];
  }
  if (edgeid==1){
    // direction vector of the required edge
    t1[0]=x[2]-x[1];
    t1[1]=y[2]-y[1];
  }
  if (edgeid==2){
    // direction vector of the required edge
    t1[0]=x[3]-x[2];
    t1[1]=y[3]-y[2];
  }
  if (edgeid==3){
    // direction vector of the required edge
    t1[0]=x[0]-x[3];
    t1[1]=y[0]-y[3];
  }
  
  //  length of the direction vector
  norm = t1[0]*t1[0]+t1[1]*t1[1];
  if (norm<1.0e-10){
    print_err("nonpositive length of the direction vector",__FILE__,__LINE__,__func__);
    abort ();
  }

  norm = sqrt(norm);
  
  return norm;
}


/**
   function computes area of linear axisymmetric quadrilateral element surfaces
   
   @param surfid - id of surface
   @param x,y - arrays of node coordinates
   
   TKr, 29/11/2019
*/
double quadlinaxisym_surface_area (long surfid,vector &x,vector &y)
{
  double norm,r;
  vector t1(3);

  if (surfid==0){
    // direction vector of the required edge
    t1[0]=x[1]-x[0];
    t1[1]=y[1]-y[0];
    r = (x[1]+x[0])/2.0;
  }
  if (surfid==1){
    // direction vector of the required edge
    t1[0]=x[2]-x[1];
    t1[1]=y[2]-y[1];
    r = (x[2]+x[1])/2.0;
  }
  if (surfid==2){
    // direction vector of the required edge
    t1[0]=x[3]-x[2];
    t1[1]=y[3]-y[2];
    r = (x[3]+x[2])/2.0;
  }
  if (surfid==3){
    // direction vector of the required edge
    t1[0]=x[0]-x[3];
    t1[1]=y[0]-y[3];
    r = (x[0]+x[3])/2.0;
  }
  
  //  length of the direction vector
  norm = t1[0]*t1[0]+t1[1]*t1[1];
  if (norm<1.0e-10){
    print_err("nonpositive length of the direction vector",__FILE__,__LINE__,__func__);
    abort ();
  }

  norm = sqrt(norm)*r;
  
  return norm;
}


/**
   function computes volume of an tetrahedral finite element
   
   @param eid - element id (it is used only in the case of nonpositive volume)
   @param x,y,z - arrays of node coordinates
   
   JK, 3. 3. 2018
*/
double tetrahedra_volume (long eid,vector &x,vector &y,vector &z)
{
  double vol,det;
  
  det = det3d (x.a, y.a, z.a);
  
  if(det <= 0.0){
    print_err("nonpositive volume of an element %d (the volume is %e)", __FILE__, __LINE__, __func__, eid+1, det);
    abort();
  }
  
  vol = det/6.0;
  
  return vol;
}


/**
   function generates outer normal %vector to required element surface
   
   @param surfid - id of surface
   @param x,y,z - arrays of node coordinates
   @param n - %vector of the outer normal
   
   JK, 3.3.2018
*/
void tetrahedra_normal_vectors (long surfid,vector &x,vector &y,vector &z,vector &n)
{
  double norm;
  vector t1(3),t2(3);
  
  if (surfid==0){
    // first in-plane vector
    t1[0]=x[2]-x[3];
    t1[1]=y[2]-y[3];
    t1[2]=z[2]-z[3];
    //  second in-plane vector
    t2[0]=x[1]-x[3];
    t2[1]=y[1]-y[3];
    t2[2]=z[1]-z[3];
  }
  if (surfid==1){
    // first in-plane vector
    t1[0]=x[0]-x[3];
    t1[1]=y[0]-y[3];
    t1[2]=z[0]-z[3];
    //  second in-plane vector
    t2[0]=x[2]-x[3];
    t2[1]=y[2]-y[3];
    t2[2]=z[2]-z[3];
  }
  if (surfid==2){
    // first in-plane vector
    t1[0]=x[1]-x[3];
    t1[1]=y[1]-y[3];
    t1[2]=z[1]-z[3];
    //  second in-plane vector
    t2[0]=x[0]-x[3];
    t2[1]=y[0]-y[3];
    t2[2]=z[0]-z[3];
  }
  if (surfid==3){
    // first in-plane vector
    t1[0]=x[1]-x[0];
    t1[1]=y[1]-y[0];
    t1[2]=z[1]-z[0];
    //  second in-plane vector
    t2[0]=x[2]-x[0];
    t2[1]=y[2]-y[0];
    t2[2]=z[2]-z[0];
  }
  
  //  outer normal vector to the first element surface
  n[0]=t1[1]*t2[2]-t1[2]*t2[1];
  n[1]=t1[2]*t2[0]-t1[0]*t2[2];
  n[2]=t1[0]*t2[1]-t1[1]*t2[0];
  
  //  length of the normal vector
  norm = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  if (norm<1.0e-10){
    print_err("nonpositive length of the normal vector",__FILE__,__LINE__,__func__);
    abort ();
  }
  
  n[0]/=norm;
  n[1]/=norm;
  n[2]/=norm;
}

/**
   function computes areas of the tetrahedron surfaces
   
   @param surfid - id of surface
   @param x,y,z - arrays of node coordinates
   
   JK, 4.3.2018
*/
double tetrahedra_surface_areas (long surfid,vector &x,vector &y,vector &z)
{
  double area,l,l1,l2,l3;
  
  if (surfid==0){
    l1=sqrt((x[2]-x[1])*(x[2]-x[1]) + (y[2]-y[1])*(y[2]-y[1]) + (z[2]-z[1])*(z[2]-z[1]));
    l2=sqrt((x[3]-x[2])*(x[3]-x[2]) + (y[3]-y[2])*(y[3]-y[2]) + (z[3]-z[2])*(z[3]-z[2]));
    l3=sqrt((x[1]-x[3])*(x[1]-x[3]) + (y[1]-y[3])*(y[1]-y[3]) + (z[1]-z[3])*(z[1]-z[3]));
  }
  if (surfid==1){
    l1=sqrt((x[0]-x[2])*(x[0]-x[2]) + (y[0]-y[2])*(y[0]-y[2]) + (z[0]-z[2])*(z[0]-z[2]));
    l2=sqrt((x[0]-x[3])*(x[0]-x[3]) + (y[0]-y[3])*(y[0]-y[3]) + (z[0]-z[3])*(z[0]-z[3]));
    l3=sqrt((x[3]-x[2])*(x[3]-x[2]) + (y[3]-y[2])*(y[3]-y[2]) + (z[3]-z[2])*(z[3]-z[2]));
  }
  if (surfid==2){
    l1=sqrt((x[1]-x[0])*(x[1]-x[0]) + (y[1]-y[0])*(y[1]-y[0]) + (z[1]-z[0])*(z[1]-z[0]));
    l2=sqrt((x[3]-x[1])*(x[3]-x[1]) + (y[3]-y[1])*(y[3]-y[1]) + (z[3]-z[1])*(z[3]-z[1]));
    l3=sqrt((x[3]-x[0])*(x[3]-x[0]) + (y[3]-y[0])*(y[3]-y[0]) + (z[3]-z[0])*(z[3]-z[0]));
  }
  if (surfid==3){
    l1=sqrt((x[1]-x[0])*(x[1]-x[0]) + (y[1]-y[0])*(y[1]-y[0]) + (z[1]-z[0])*(z[1]-z[0]));
    l2=sqrt((x[2]-x[1])*(x[2]-x[1]) + (y[2]-y[1])*(y[2]-y[1]) + (z[2]-z[1])*(z[2]-z[1]));
    l3=sqrt((x[2]-x[0])*(x[2]-x[0]) + (y[2]-y[0])*(y[2]-y[0]) + (z[2]-z[0])*(z[2]-z[0]));
  }
  
  // Heron's formula for the triangle area
  l=(l1+l2+l3)/2.0;
  area=sqrt(l*(l-l1)*(l-l2)*(l-l3));
  
  return area;
}



/**
   function computes volume of an hexahedral finite element
      
   @param eid - element id (it is used only in the case of nonpositive volume)
   @param x,y,z - arrays of node coordinates

   
   TKr 26/09/2022 according to JK
*/
double hexahedra_volume (long eid,vector &x,vector &y,vector &z)
{
  long i,j,k;
  double xi,eta,zeta,ww1,ww2,ww3,jac,vol;
  vector w(2),gp(2);
  
  gauss_points (gp.a,w.a,2);
  
  vol=0.0;
  for (i=0;i<2;i++){
    xi=gp[i];  ww1=w[i];
    for (j=0;j<2;j++){
      eta=gp[j];  ww2=w[j];
      for (k=0;k<2;k++){
	zeta=gp[k]; ww3=w[k];
	
	jac_3d (jac,x,y,z,xi,eta,zeta);
	
	if (jac<0.0){
	  print_err("\n negative Jacobian on element %ld", __FILE__, __LINE__, __func__,eid);
	  jac=fabs(jac);
	}
	
	jac*=ww1*ww2*ww3;
	vol+=jac;
      }
    }
  }
  
  return vol;
}


/**
   function generates outer normal %vector to required element surface
   
   @param surfid - id of surface
   @param x,y,z - arrays of node coordinates
   @param n - %vector of the outer normal
   
   TKr 26/09/2022 according to JK
*/
void hexahedra_normal_vectors (long surfid,vector &x,vector &y,vector &z,vector &n)
{
  double norm;
  double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
  vector t1(3),t2(3);
  
  //two vectors are selected - one edge and one diagonal
  
  if (surfid==0){
    x1 = x[0];
    x2 = x[3];
    x3 = x[7];
    x4 = x[4];
    y1 = y[0];
    y2 = y[3];
    y3 = y[7];
    y4 = y[4];
    z1 = z[0];
    z2 = z[3];
    z3 = z[7];
    z4 = z[4];
  }
  if (surfid==1){
    x1 = x[1];
    x2 = x[0];
    x3 = x[4];
    x4 = x[5];
    y1 = y[1];
    y2 = y[0];
    y3 = y[4];
    y4 = y[5];
    z1 = z[1];
    z2 = z[0];
    z3 = z[4];
    z4 = z[5];
  }
  if (surfid==2){
    x1 = x[2];
    x2 = x[1];
    x3 = x[5];
    x4 = x[6];
    y1 = y[2];
    y2 = y[1];
    y3 = y[5];
    y4 = y[6];
    z1 = z[2];
    z2 = z[1];
    z3 = z[5];
    z4 = z[6];
  }
  if (surfid==3){
    x1 = x[3];
    x2 = x[2];
    x3 = x[6];
    x4 = x[7];
    y1 = y[3];
    y2 = y[2];
    y3 = y[6];
    y4 = y[7];
    z1 = z[3];
    z2 = z[2];
    z3 = z[6];
    z4 = z[7];
  }
  if (surfid==4){
    x1 = x[0];
    x2 = x[1];
    x3 = x[2];
    x4 = x[3];
    y1 = y[0];
    y2 = y[1];
    y3 = y[2];
    y4 = y[3];
    z1 = z[0];
    z2 = z[1];
    z3 = z[2];
    z4 = z[3];
  }
  if (surfid==5){
    x1 = x[4];
    x2 = x[5];
    x3 = x[6];
    x4 = x[7];
    y1 = y[4];
    y2 = y[5];
    y3 = y[6];
    y4 = y[7];
    z1 = z[4];
    z2 = z[5];
    z3 = z[6];
    z4 = z[7];
  }
  
  // first in-plane vector
  t1[0]=x2-x1;
  t1[1]=y2-y1;
  t1[2]=z2-z1;
  //  second in-plane vector
  t2[0]=x3-x1;
  t2[1]=y3-y1;
  t2[2]=z3-z1;
  
  //  outer normal vector to the first element surface
  n[0]=t1[1]*t2[2]-t1[2]*t2[1];
  n[1]=t1[2]*t2[0]-t1[0]*t2[2];
  n[2]=t1[0]*t2[1]-t1[1]*t2[0];
  
  //  length of the normal vector
  norm = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  if (norm<1.0e-10){
    print_err("nonpositive length of the normal vector",__FILE__,__LINE__,__func__);
    abort ();
  }
  
  n[0]/=norm;
  n[1]/=norm;
  n[2]/=norm;
}

/**
   function computes areas of the tetrahedron surfaces
   
   @param surfid - id of surface
   @param x,y,z - arrays of node coordinates
   
   TKr 26/09/2022 according to JK
*/
double hexahedra_surface_areas (long surfid,vector &x,vector &y,vector &z)
{
  double area,l,l1,l2,l3,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
  
  if (surfid==0){
    x1 = x[0];
    x2 = x[3];
    x3 = x[7];
    x4 = x[4];
    y1 = y[0];
    y2 = y[3];
    y3 = y[7];
    y4 = y[4];
    z1 = z[0];
    z2 = z[3];
    z3 = z[7];
    z4 = z[4];
  }
  if (surfid==1){
    x1 = x[1];
    x2 = x[0];
    x3 = x[4];
    x4 = x[5];
    y1 = y[1];
    y2 = y[0];
    y3 = y[4];
    y4 = y[5];
    z1 = z[1];
    z2 = z[0];
    z3 = z[4];
    z4 = z[5];
  }
  if (surfid==2){
    x1 = x[2];
    x2 = x[1];
    x3 = x[5];
    x4 = x[6];
    y1 = y[2];
    y2 = y[1];
    y3 = y[5];
    y4 = y[6];
    z1 = z[2];
    z2 = z[1];
    z3 = z[5];
    z4 = z[6];
  }
  if (surfid==3){
    x1 = x[3];
    x2 = x[2];
    x3 = x[6];
    x4 = x[7];
    y1 = y[3];
    y2 = y[2];
    y3 = y[6];
    y4 = y[7];
    z1 = z[3];
    z2 = z[2];
    z3 = z[6];
    z4 = z[7];
  }
  if (surfid==4){
    x1 = x[0];
    x2 = x[1];
    x3 = x[2];
    x4 = x[3];
    y1 = y[0];
    y2 = y[1];
    y3 = y[2];
    y4 = y[3];
    z1 = z[0];
    z2 = z[1];
    z3 = z[2];
    z4 = z[3];
  }
  if (surfid==5){
    x1 = x[4];
    x2 = x[5];
    x3 = x[6];
    x4 = x[7];
    y1 = y[4];
    y2 = y[5];
    y3 = y[6];
    y4 = y[7];
    z1 = z[4];
    z2 = z[5];
    z3 = z[6];
    z4 = z[7];
  }
  
  //the quad (1,2,3,4) is split into two triangles
  //the first triangle with nodes 1,2,4
  l1=sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
  l2=sqrt((x4-x2)*(x4-x2) + (y4-y2)*(y4-y2) + (z4-z2)*(z4-z2));
  l3=sqrt((x1-x4)*(x1-x4) + (y1-y4)*(y1-y4) + (z1-z4)*(z1-z4));
  
  // Heron's formula for the first triangle area
  l=(l1+l2+l3)/2.0;
  area = sqrt(l*(l-l1)*(l-l2)*(l-l3));
  
  //second triangle with nodes 2,3,4
  l1=sqrt((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2));
  l2=sqrt((x4-x3)*(x4-x3) + (y4-y3)*(y4-y3) + (z4-z3)*(z4-z3));
  l3=sqrt((x2-x4)*(x2-x4) + (y2-y4)*(y2-y4) + (z2-z4)*(z2-z4));
  
  // Heron's formula for the triangle area
  l=(l1+l2+l3)/2.0;
  area = area + sqrt(l*(l-l1)*(l-l2)*(l-l3));

  return area;
}
