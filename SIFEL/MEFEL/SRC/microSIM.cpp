#include "microSIM.h"
#include "global.h"
#include "mechmat.h"
#include <math.h>
#include <stdlib.h>



microSIM::microSIM (void)
{
  c1 = 6.20e-1;
  c2 = 2.76;
  c4 = 70.;
  c5 = 2.50;
  c6 = 1.30;
  c7 = 50.;
  c8 = 8.0;
  c9 = 1.30;
  c10 = 0.73;
  c11 = 0.2;
  c12 = 7000.;
  c13 = 0.23;
  c14 = 0.8;
  c15 = 1.0;
  c16 = 0.02;
  c17 = 0.01;
  c18 = 1.0;
  c19 = 0.4; 
  mu = 1.0;
  
  kronecker.a = new double [6];
  kronecker.n = 6;

  // kronecker delta
  kronecker [0]= 1.;
  kronecker [1]= 1.;
  kronecker [2]= 1.;
  kronecker [3]= 0.;
  kronecker [4]= 0.;
  kronecker [5]= 0.;
}
microSIM::~microSIM (void)
{
}

void microSIM::read(XFILE *in)
{
  xfscanf(in,"%ld %lf %lf %lf %lf %lf %lf %lf %lf",
	 &numberOfMicroplanes,&k1,&k2,&k3,&k4,&c3,&c20,&e,&nu);
  
  ev = e/(1-2*nu);
  ed = 5*e/(2+3*mu)/(1+nu);
  et = mu*ed;
  //compute microplaneweights and projection tensors
  initializeData (numberOfMicroplanes);
}

void microSIM::initializeData (long numberOfMicroplanes )
{
  //computes microplane normals, weights and projection tensors
  long i,i4,mPlane;
  long ij[6][2] = {{1,1},{2,2},{3,3},{2,3},{3,1},{1,2}};
  matrix microplaneNormals(numberOfMicroplanes,3);
  
  projN.a = new double [numberOfMicroplanes*6];
  projN.m = numberOfMicroplanes;
  projN.n = 6;
  
  microplaneWeights.a = new double [numberOfMicroplanes];
  microplaneWeights.n = numberOfMicroplanes;
  

  // definition of intergation points and weights
  
  if (numberOfMicroplanes == 28) {
    for (i=0; i<4; i++) {
      microplaneWeights[i]    = 0.0160714276E0;
      microplaneWeights[i+4]  = 0.0204744730E0;
      microplaneWeights[i+8]  = 0.0204744730E0;
      microplaneWeights[i+12] = 0.0204744730E0;
      microplaneWeights[i+16] = 0.0158350505E0;
      microplaneWeights[i+20] = 0.0158350505E0;
      microplaneWeights[i+24] = 0.0158350505E0;
    }
    double help[7][3]={{0.577350259E0, 0.577350259E0, 0.577350259E0},
		       {0.935113132E0, 0.250562787E0, 0.250562787E0},
		       {0.250562787E0, 0.935113132E0, 0.250562787E0},
		       {0.250562787E0, 0.250562787E0, 0.935113132E0},
		       {0.186156720E0, 0.694746614E0, 0.694746614E0},
		       {0.694746614E0, 0.186156720E0, 0.694746614E0},
		       {0.694746614E0, 0.694746614E0, 0.186156720E0}};
    for (i=0; i< 7; i++) {
      i4 = i * 4;
      microplaneNormals[i4  ][0] = help[i][0];
      microplaneNormals[i4  ][1] = help[i][1];
      microplaneNormals[i4  ][2] = help[i][2];
      microplaneNormals[i4+1][0] = help[i][0];
      microplaneNormals[i4+1][1] = help[i][1];
      microplaneNormals[i4+1][2] =-help[i][2];
      microplaneNormals[i4+2][0] = help[i][0];
      microplaneNormals[i4+2][1] =-help[i][1];
      microplaneNormals[i4+2][2] = help[i][2];
      microplaneNormals[i4+3][0] = help[i][0];
      microplaneNormals[i4+3][1] =-help[i][1];
      microplaneNormals[i4+3][2] =-help[i][2];
    }
  } 
  else if (numberOfMicroplanes == 21){
    for (i=0; i<3; i++) microplaneWeights[i]    = 0.02652141274;
    for (i=3; i<9; i++) microplaneWeights[i]    = 0.01993014153;
    for (i=9; i<21; i++) microplaneWeights[i]    = 0.02507124272;
    
    microplaneNormals[0][0] =  microplaneNormals[1][1] =  microplaneNormals[2][2] = 1. ;
    microplaneNormals[0][1] =  microplaneNormals[0][2] =  microplaneNormals[1][0] = 
      microplaneNormals[1][2] = microplaneNormals[2][0] = microplaneNormals[2][1] = 0.;
    microplaneNormals[3][0] =  microplaneNormals[3][1] =  microplaneNormals[4][0] = 
      microplaneNormals[5][0] = microplaneNormals[5][2] = microplaneNormals[6][0] =
      microplaneNormals[7][1] = microplaneNormals[7][2] = microplaneNormals[8][1] =
      0.7071067812;
    microplaneNormals[4][1] =  microplaneNormals[6][2] =  microplaneNormals[8][2] = 
      -0.7071067812;
    microplaneNormals[3][2] =  microplaneNormals[4][2] =  microplaneNormals[5][1] =
      microplaneNormals[6][1] =  microplaneNormals[7][0] =  microplaneNormals[8][0] = 0.;
    
    
    double help [3][3] = {{0.3879072746, 0.3879072746, 0.8360956240},
			    {0.3879072746, 0.8360956240, 0.3879072746},
			    {0.8360956240, 0.3879072746, 0.3879072746}};
    
    for(i=0; i<3;i++){
      i4 = i * 4;
      microplaneNormals[9+i4][0] = help [i][0];
      microplaneNormals[9+i4][1] = help [i][1];
      microplaneNormals[9+i4][2] = help [i][2];
      
      microplaneNormals[10+i4][0] = help [i][0];
      microplaneNormals[10+i4][1] = help [i][1];
      microplaneNormals[10+i4][2] = -help [i][2];
      
      microplaneNormals[11+i4][0] = help [i][0];
      microplaneNormals[11+i4][1] = -help [i][1];
      microplaneNormals[11+i4][2] = help [i][2];
      
      microplaneNormals[12+i4][0] = help [i][0];
      microplaneNormals[12+i4][1] = -help [i][1];
      microplaneNormals[12+i4][2] = -help [i][2];
    }
  }
  else  if(numberOfMicroplanes == 61){
    
    double help [61][4]={
      {1.000000000000, 0.000000000000, 0.000000000000, 0.00795844204678},
      {0.745355992500, 0.0           , 0.666666666667, 0.00795844204678},
      {0.745355992500,-0.577350269190,-0.333333333333, 0.00795844204678},
      {0.745355992500, 0.577350269190,-0.333333333333, 0.00795844204678},
      {0.333333333333, 0.577350269190, 0.745355992500, 0.00795844204678},
      {0.333333333333,-0.577350269190, 0.745355992500, 0.00795844204678},
      {0.333333333333,-0.934172358963, 0.127322003750, 0.00795844204678},
      {0.333333333333,-0.356822089773,-0.872677996250, 0.00795844204678},
      {0.333333333333, 0.356822089773,-0.872677996250, 0.00795844204678},
      {0.333333333333, 0.934172358963, 0.127322003750, 0.00795844204678},
      {0.794654472292,-0.525731112119, 0.303530999103, 0.0105155242892},
      {0.794654472292, 0.0           ,-0.607061998207, 0.0105155242892},
      {0.794654472292, 0.525731112119, 0.303530999103, 0.0105155242892},
      {0.187592474085, 0.0           , 0.982246946377, 0.0105155242892},
      {0.187592474085,-0.850650808352,-0.491123473188, 0.0105155242892},
      {0.187592474085, 0.850650808352,-0.491123473188, 0.0105155242892},
      {0.934172358963, 0.0           , 0.356822089773, 0.0100119364272},
      {0.934172358963,-0.309016994375,-0.178411044887, 0.0100119364272},
      {0.934172358963, 0.309016994375,-0.178411044887, 0.0100119364272},
      {0.577350269190, 0.309016994375, 0.755761314076, 0.0100119364272},
      {0.577350269190,-0.309016994375, 0.755761314076, 0.0100119364272},
      {0.577350269190,-0.809016994375,-0.110264089708, 0.0100119364272},
      {0.577350269190,-0.5           ,-0.645497224368, 0.0100119364272},
      {0.577350269190, 0.5           ,-0.645497224368, 0.0100119364272},
      {0.577350269190, 0.809016994375,-0.110264089708, 0.0100119364272},
      {0.356822089773,-0.809016994375, 0.467086179481, 0.0100119364272},
      {0.356822089773, 0.0           ,-0.934172358963, 0.0100119364272},
      {0.356822089773, 0.809016994375, 0.467086179481, 0.0100119364272},
      {0.0           , 0.5           , 0.866025403784, 0.0100119364272},
      {0.0           ,-1.0           , 0.0           , 0.0100119364272},
      {0.0           , 0.5           ,-0.866025403784, 0.0100119364272},
      {0.947273580412,-0.277496978165, 0.160212955043, 0.00690477957966},
      {0.812864676392,-0.277496978165, 0.512100034157, 0.00690477957966},
      {0.595386501297,-0.582240127941, 0.553634669695, 0.00690477957966},
      {0.595386501297,-0.770581752342, 0.227417407053, 0.00690477957966},
      {0.812864676392,-0.582240127941,-0.015730584514, 0.00690477957966},
      {0.492438766306,-0.753742692223,-0.435173546254, 0.00690477957966},
      {0.274960591212,-0.942084316623,-0.192025554687, 0.00690477957966},
      {-0.076926487903,-0.942084316623,-0.326434458707, 0.00690477957966},
      {-0.076926487903,-0.753742692223,-0.652651721349, 0.00690477957966},
      {0.274960591212,-0.637341166847,-0.719856173359, 0.00690477957966},
      {0.947273580412, 0.0           ,-0.320425910085, 0.00690477957966},
      {0.812864676392,-0.304743149777,-0.496369449643, 0.00690477957966},
      {0.595386501297,-0.188341624401,-0.781052076747, 0.00690477957966},
      {0.595386501297, 0.188341624401,-0.781052076747, 0.00690477957966},
      {0.812864676392, 0.304743149777,-0.496369449643, 0.00690477957966},
      {0.492438766306, 0.753742692223,-0.435173546254, 0.00690477957966},
      {0.274960591212, 0.637341166847,-0.719856173359, 0.00690477957966},
      {-0.076926487903, 0.753742692223,-0.652651721349, 0.00690477957966},
      {-0.076926487903, 0.942084316623,-0.326434458707, 0.00690477957966},
      {0.274960591212, 0.942084316623,-0.192025554687, 0.00690477957966},
      {0.947273580412, 0.277496978165, 0.160212955043, 0.00690477957966},
      {0.812864676392, 0.582240127941,-0.015730584514, 0.00690477957966},
      {0.595386501297, 0.770581752342, 0.227417407053, 0.00690477957966},
      {0.595386501297, 0.582240127941, 0.553634669695, 0.00690477957966},
      {0.812864676392, 0.277496978165, 0.512100034157, 0.00690477957966},
      {0.492438766306, 0.0           , 0.870347092509, 0.00690477957966},
      {0.274960591212, 0.304743149777, 0.911881728046, 0.00690477957966},
      {-0.076926487903, 0.188341624401, 0.979086180056, 0.00690477957966},
      {-0.076926487903,-0.188341624401, 0.979086180056, 0.00690477957966},
      {0.274960591212,-0.304743149777, 0.911881728046, 0.00690477957966}};
    
    for(i=0; i<numberOfMicroplanes; i++){
      microplaneNormals[i][0]= help[i][0];
      microplaneNormals[i][1]= help[i][1];
      microplaneNormals[i][2]= help[i][2];
      microplaneWeights[i]   = help[i][3];
    }
  }else{
    fprintf (stderr,"microSIM::initializeData   Unsupported number of microplanes");
  }

// compute projection tensors for each microplane
  long ii,jj;
  vector n(3), m(3), l(3);

  for (mPlane=0; mPlane < numberOfMicroplanes; mPlane++) {
    n[0] = microplaneNormals[mPlane][0];
    n[1] = microplaneNormals[mPlane][1];
    n[2] = microplaneNormals[mPlane][2];

    for (i=0; i<6; i++) {
      ii = ij[i][0]-1;
      jj = ij[i][1]-1;

      projN[mPlane][i] = n(ii)*n(jj);
    }
  }
  return;
}

void microSIM::matstiff (matrix &d,long ipp,long /*ido*/)
  //  function returns elastic stiffness matrix
  //  d - elastic stiffness matrix
{
  Mm->elmatstiff (d,ipp);
}

void microSIM::nlstresses (long ipp,long ido)
{
  long i, i_other,mPlaneIndex;
  double previousSigmaV,previousSigmaN,previousSigmaD;
  double sev,sed;
  double cv,cd,meanN=0.;
  double sigmaN,sigmaV,sigmaD;
  double epsV,depsV, epsD, epsN, depsD,depsN;
  vector macroStrain(6), macroStrainIncrement(6), macroSigma(6);
  
  //  initial values of strain vector and strain increment
  for (i=0;i<6;i++){
    macroStrain[i] = Mm->ip[ipp].strain[i];
    macroStrainIncrement[i] = Mm->ip[ipp].strain[i] - Mm->ip[ipp].eqother[ido+i+numberOfMicroplanes+1];

    //UPDATE history variables= store total macroStrain
    Mm->ip[ipp].other[ido+i+numberOfMicroplanes+1]=Mm->ip[ipp].strain[i];
  }

  //volumetric microstrain and microstrain increment
  epsV = (macroStrain[0]+macroStrain[1]+macroStrain[2])/3.0;
  depsV= (macroStrainIncrement[0]+macroStrainIncrement[1]+macroStrainIncrement[2])/3.0;

  //ask for previous volumetric stress
  previousSigmaV=Mm->ip[ipp].eqother[ido+0];
  
  //loop over all microplanes
  for (mPlaneIndex=0; mPlaneIndex < numberOfMicroplanes; mPlaneIndex++) {
    
    //compute total microstrains on microplane and its increments
    epsN=0., depsN=0.;
    for (i=0; i<6; i++) {
      epsN += projN [mPlaneIndex][i]*macroStrain[i];
      depsN += projN [mPlaneIndex][i]*macroStrainIncrement[i];
    }
    epsD = epsN - epsV;
    depsD = depsN - depsV;
    
    //ask for previous stress components on microplane
    i_other=ido+mPlaneIndex+1;//poradi komponenty v poli history parametru L,M,N
    previousSigmaN=Mm->ip[ipp].eqother[i_other];
    previousSigmaD=previousSigmaN-previousSigmaV;
    
    //evaluation of new microplane stresses

    //secant moduli
    cv = ev;
    cd = ed;
    
    //volumetric microstress
    sev=previousSigmaV + cv*depsV;
    sigmaV=minim( maxim(sev,FVminus(epsV)), FVplus(epsV));
    //deviatoric microstress
    sed=previousSigmaD + cd*depsD; 
    sigmaD=minim(maxim(sed,FDminus(epsD)),FDplus(epsD));
    //normal microstress
    sigmaN=minim( (sigmaV+sigmaD) ,FN(epsN,previousSigmaV));
     
    //mean normal microstress (computed after sweeping throgh all microplanes)
    meanN += sigmaN* microplaneWeights[mPlaneIndex]*6.0;

    //UPDATE history variables, i.e. store microstresses
    Mm->ip[ipp].other[i_other]=sigmaN;
    
  }//end 1st loop over nmp
    
  //volumetric stress is the same for all  mplanes
  //and does not need to be homogenized .
  //Only updating accordinging to mean normal stress must be done.
     if (sigmaV > (meanN/3.))  {sigmaV = meanN/3.;}

  //UPDATE history variables= store posibly new sigmaV
  Mm->ip[ipp].other[ido+0]=sigmaV;
  
  fillv(0,macroSigma);//zero macroSigma

  //compute macrostresstensor
  for (mPlaneIndex=0; mPlaneIndex < numberOfMicroplanes; mPlaneIndex++) {
    i_other=ido+mPlaneIndex+1;//poradi komponenty v poli history parametru L,M,N
    sigmaN=Mm->ip[ipp].other[i_other];
    //recalculation of sigmaD according to possibly new sigmaV 
    sigmaD = sigmaN - sigmaV;
    
    
    for (i=0; i<6; i++) {
      macroSigma[i]+=(projN[mPlaneIndex][i]-kronecker[i]/3.)*sigmaD*microplaneWeights[mPlaneIndex]*6.0;  
    }
    
  }//end 2nd loop over nmp
  
  //2nd constraint, addition of volumetric part
  macroSigma[0]+=sigmaV;
  macroSigma[1]+=sigmaV;
  macroSigma[2]+=sigmaV;
  
  //send new macrostress vector to the solver
  Mm->storestress(0,ipp,macroSigma);

}

void microSIM::updateval(long ipp,long im,long ido)
{
  long i,n = Mm->givencompother(ipp, im);

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }
}

inline double
microSIM::macbra(double x) /* Macauley bracket = positive part of x */
{
  return(maxim(x,0.0));
}

inline double
microSIM::FVplus (double epsV)
/*positive volumetric boundary */
{
  return(ev*k1*c13/(1+(c14/k1)*macbra(epsV-c13*c15*k1)));
}

inline double
microSIM::FVminus (double epsV)
/*negative volumetric boundary */
{
return(-e*k1*k3*exp(-epsV/(k1*k4)));
}

inline double
microSIM::FDminus(double epsD)
 /*negative deviatoric boundary */
 {
   double a;
   a=macbra(-epsD-c8*c9*k1)/(k1*c7);
   return(-e*k1*c8/(1+a*a));
 }

inline double
microSIM::FDplus(double epsD)
/*positive deviatoric bondary */
{
  double a;
  a=(macbra(epsD-c5*c6*k1)/(k1*c20*c7));
  return (e*k1*c5/(1+a*a));
}

inline double
microSIM::FN(double epsN,double sigmaV)
 /*normal boundary */
{
  return(e*k1*c1*exp(-macbra(epsN-c1*c2*k1)/(k1*c3+macbra(-c4*(sigmaV/ev)))));
}

inline double
microSIM::maxim (double a,double b)
{
  return(a>b ? a:b);
}

inline double
microSIM::minim (double a,double b)
{
  return(a<b ? a:b);
}
