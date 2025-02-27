#include "crsec3dbeam.h"
#include "vector.h"
#include "global.h"
#include "probdesc.h"
#include "stochdriver.h"

crsec3dbeam::crsec3dbeam (void)
{
  a=0.0;  ix=0.0;  iy=0.0;  iz=0.0;  shearcoeffy=0.0;  shearcoeffz=0.0;  rho=0.0;
  allocv(3, lcs);
  fillv(0.0, lcs);
}
crsec3dbeam::~crsec3dbeam (void)
{
  
}

/**
   function read parameters of cross-section of 3D beams
   
   @param in - input stream
*/
void crsec3dbeam::read (XFILE *in)
{
  switch (Mp->tprob){
  case linear_statics:{
    xfscanf (in,"%k%lf  %k%lf %k%lf %k%lf  %k%lf %k%lf", "a", &a, "ix", &ix, "iy", &iy, "iz", &iz,
             "kappa_y", &shearcoeffy, "kappa_z", &shearcoeffz);
    xfscanf(in, "%k", "loc_z");
    readv(in, lcs);//  function from the file vector.cpp
    break;
  }
  case eigen_dynamics:{
    xfscanf (in,"%k%lf  %k%lf %k%lf %k%lf  %k%lf %k%lf", "a", &a, "ix", &ix, "iy", &iy, "iz", &iz,
             "kappa_y", &shearcoeffy, "kappa_z", &shearcoeffz, "rho", &rho);
    xfscanf(in, "%k", "loc_z");
    readv(in, lcs);//  function from the file vector.cpp
    break;
  }
  case forced_dynamics:{
    xfscanf (in,"%k%lf  %k%lf %k%lf %k%lf  %k%lf %k%lf", "a", &a, "ix", &ix, "iy", &iy, "iz", &iz,
             "kappa_y", &shearcoeffy, "kappa_z", &shearcoeffz, "rho", &rho);
    xfscanf(in, "%k", "loc_z");
    readv(in, lcs);//  function from the file vector.cpp
    break;
  }
  case growing_mech_structure:
  case mech_timedependent_prob:{
    xfscanf (in,"%k%lf  %k%lf %k%lf %k%lf  %k%lf %k%lf", "a", &a, "ix", &ix, "iy", &iy, "iz", &iz,
             "kappa_y", &shearcoeffy, "kappa_z", &shearcoeffz);
    xfscanf(in, "%k", "loc_z");
    readv(in, lcs);//  function from the file vector.cpp
    break;
  }
  case earth_pressure:{
    xfscanf (in,"%k%lf  %k%lf %k%lf %k%lf  %k%lf %k%lf", "a", &a, "ix", &ix, "iy", &iy, "iz", &iz,
             "kappa_y", &shearcoeffy, "kappa_z", &shearcoeffz);
    xfscanf(in, "%k", "loc_z");
    readv(in, lcs);//  function from the file vector.cpp
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}


/**
   function prints parameters of cross-section of 3D beams
   
   @param out - output stream

   TKr, 11/02/2013, according to read (XFILE *in)
*/
void crsec3dbeam::print (FILE *out)
{
  switch (Mp->tprob){
  case linear_statics:{
    fprintf (out,"\n %le  %le %le %le  %le %le",a,ix,iy,iz,shearcoeffy,shearcoeffz);
    printv(out, lcs);//  function from the file vector.cpp
    break;
  }
  case eigen_dynamics:{
    fprintf (out,"\n %le  %le %le %le  %le %le  %le",a,ix,iy,iz,shearcoeffy,shearcoeffz,rho);
    printv(out, lcs);//  function from the file vector.cpp
    break;
  }
  case forced_dynamics:{
    fprintf (out,"\n %le  %le %le %le  %le %le  %le",a,ix,iy,iz,shearcoeffy,shearcoeffz,rho);
    printv(out, lcs);//  function from the file vector.cpp
    break;
  }
  case growing_mech_structure:
  case mech_timedependent_prob:{
    fprintf (out,"\n %le  %le %le %le  %le %le",a,ix,iy,iz,shearcoeffy,shearcoeffz);
    printv(out, lcs);//  function from the file vector.cpp
    break;
  }
  case earth_pressure:{
    fprintf (out,"\n %le  %le %le %le  %le %le",a,ix,iy,iz,shearcoeffy,shearcoeffz);
    printv(out, lcs);//  function from the file vector.cpp
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function returns area of cross section
*/
double crsec3dbeam::give_area ()
{
  return a;
}

/**
   function returns quadratic moment with respect to x axis
*/
double crsec3dbeam::give_ix ()
{
  return ix;
}

/**
   function returns moment of inertia with respect to y axis
*/
double crsec3dbeam::give_iy ()
{
  return iy;
}

/**
   function returns moment of inertia with respect to z axis
*/
double crsec3dbeam::give_iz ()
{
  return iz;
}

/**
   function returns all quadratic moments
*/
void crsec3dbeam::give_moments (double *i)
{
  i[0]=ix;
  i[1]=iy;
  i[2]=iz;
}

/**
   function returns shear coefficients
   the coefficients are used in Mindlin theory
*/
void crsec3dbeam::give_shearcoeff (double *sc)
{
  sc[0]=shearcoeffy;
  sc[1]=shearcoeffz;
}

/**
   function returns density
*/
double crsec3dbeam::give_density ()
{
  return rho;
}

void crsec3dbeam::changeparam (atsel &atcs,vector &val)
{
  long i;
  
  for (i=0;i<atcs.num;i++){
    switch (atcs.atrib[i]){
    case 0:{
      a=val[i];
      break;
    }
    case 1:{
      ix=val[i];
      break;
    }
    case 2:{
      iy=val[i];
      break;
    }
    case 3:{
      iz=val[i];
      break;
    }
    case 4:{
      shearcoeffy=val[i];
      break;
    }
    case 5:{
      shearcoeffz=val[i];
      break;
    }
    case 6:{
      rho=val[i];
      break;
    }
    default:{
      print_err("wrong atribute number is required",__FILE__,__LINE__,__func__);
    }
    }
  }
}
