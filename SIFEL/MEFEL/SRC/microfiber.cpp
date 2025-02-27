#include "microfiber.h"
#include "global.h"
#include "mechmat.h"
#include "intpoints.h"
#include <math.h>
#include <stdlib.h>





microfiber::microfiber (void)
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
microfiber::~microfiber (void)
{
}

void microfiber::read(XFILE *in)
{
  xfscanf(in,"%ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
	 &numberOfMicroplanes,&k1,&k2,&k3,&k4,&c3,&c20,&e,&nu,&k10,&k11,&k12,&k13);
  
  ev = e/(1-2*nu);
  ed = 5*e/(2+3*mu)/(1+nu);
  et = mu*ed;
  //compute microplaneweights and projection tensors
  initializeData (numberOfMicroplanes);
}

void microfiber::initializeData (long numberOfMicroplanes )
{
  //computes microplane normals, weights and projection tensors
  long i,i4,mPlane;
  long ij[6][2] = {{1,1},{2,2},{3,3},{2,3},{3,1},{1,2}};
  matrix microplaneNormals(numberOfMicroplanes,3);
  
  projL.a = new double [numberOfMicroplanes*6];
  projL.m = numberOfMicroplanes;
  projL.n = 6;

  projM.a = new double [numberOfMicroplanes*6];
  projM.m = numberOfMicroplanes;
  projM.n = 6;

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

    /* -------- */
    FILE * norm;
    norm= fopen("normals.dat","w");
    for(i=0;i<28;i++){
      fprintf(norm,"%ld %f %f %f\n",i,microplaneNormals[i][0],microplaneNormals[i][1],microplaneNormals[i][2]);
    }
    /* ------- */

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
    fprintf (stderr,"microfiber::initializeData   Unsupported number of microplanes");
  }
  
// compute projection tensors for each microplane
  long ii,jj;
  vector n(3), m(3), l(3);
  double aux;
  
  for (mPlane=0; mPlane < numberOfMicroplanes; mPlane++) {
    n[0] = microplaneNormals[mPlane][0];
    n[1] = microplaneNormals[mPlane][1];
    n[2] = microplaneNormals[mPlane][2];
    
    if ((mPlane+1)%3 == 0) {
      aux = sqrt (n[0]*n[0] + n[1]*n[1]);
      if(fabs(aux) > 1.0e-12){
	m[0] = n[1]/aux;
	m[1] = -n[0]/aux;
	m[2] = 0.0;
      }
      else {
	m[0] = 1.0;
	m[1] = 0.0;
	m[2] = 0.0;
      }
    } else if ((mPlane+1)%3 == 1) {
      aux = sqrt (n[1]*n[1] + n[2]*n[2]);
      if(fabs(aux) > 1.0e-12){
	m[0] = 0.0;
	m[1] = n[2]/aux;
	m[2] = -n[1]/aux;
      }
      else {
	m[0] = 0.0;
	m[1] = 1.0;
	m[2] = 0.0;
      }
    } else  {
      aux = sqrt (n[0]*n[0] + n[2]*n[2]);
      if(fabs(aux) > 1.0e-12){
	m[0] = n[2]/aux;
	m[1] = 0.0;
	m[2] = -n[0]/aux;
      }
      else {
	m[0] = 0.0;
	m[1] = 0.0;
	m[2] = 1.0;
      }
    }
    
    //l.beVectorProductOf (m,n);
    crprd(m,n,l);
    
    for (i=0; i<6; i++) {
      ii = ij[i][0]-1;
      jj = ij[i][1]-1;
      
      projN[mPlane][i] = n(ii)*n(jj);
      projM[mPlane][i] = 0.5*(m(ii)*n(jj) + m(jj)*n(ii)) ;
      projL[mPlane][i] = 0.5*(l(ii)*n(jj) + l(jj)*n(ii)) ;
    }
  }
  return;
}

void microfiber::matstiff (matrix &d,long ipp,long /*ido*/)
  //  function returns elastic stiffness matrix
  //  d - elastic stiffness matrix
{
  Mm->elmatstiff (d,ipp);
}

void microfiber::nlstresses (long ipp,long ido)
{
  long i, i_other,mPlaneIndex;
  double previousSigmaV,previousSigmaL,previousSigmaM,previousSigmaN,previousSigmaD,
    previousSigmaN_fiber,previousEpsN;
  double sel,sem,sev,sed,f,senf;
  double cv,cd,meanN=0., fiber_module=0;
  //double sigmaN,sigmaL,sigmaM,sigmaV,sigmaD;
  double sigmaL,sigmaM,sigmaV,sigmaD;
  double sigmaN_concrete, sigmaD_concrete, sigmaN_fiber;
  double epsV,depsV, epsD, epsN, depsD, depsM, depsL,depsN;
  vector macroStrain(6), macroStrainIncrement(6), macroSigma(6);
  
  //  initial values of strain vector and strain increment
  for (i=0;i<6;i++){
    macroStrain[i] = Mm->ip[ipp].strain[i];
    macroStrainIncrement[i] = Mm->ip[ipp].strain[i] - Mm->ip[ipp].eqother[ido+i+numberOfMicroplanes*5+1];

    //UPDATE history variables= store total macroStrain
    Mm->ip[ipp].other[ido+i+numberOfMicroplanes*5+1]=Mm->ip[ipp].strain[i];
  }
  
  //volumetric microstrain and microstrain increment
  epsV = (macroStrain[0]+macroStrain[1]+macroStrain[2])/3.0;
  depsV= (macroStrainIncrement[0]+macroStrainIncrement[1]+macroStrainIncrement[2])/3.0;

  //ask for previous volumetric stress
  previousSigmaV=Mm->ip[ipp].eqother[ido+0];
  
  //loop over all microplanes
  for (mPlaneIndex=0; mPlaneIndex < numberOfMicroplanes; mPlaneIndex++) {
    
    //compute total microstrains on microplane and its increments
    epsN=0., depsL=0., depsM=0.; depsN=0.;
    for (i=0; i<6; i++) {
      epsN += projN [mPlaneIndex][i]*macroStrain[i];
      depsL += projL [mPlaneIndex][i]*macroStrainIncrement[i];
      depsM += projM [mPlaneIndex][i]*macroStrainIncrement[i];
      depsN += projN [mPlaneIndex][i]*macroStrainIncrement[i];
    }
    epsD = epsN - epsV;
    depsD = depsN - depsV;
    
    //ask for previous stress components on microplane
    i_other=ido+(mPlaneIndex+1)*5-4;//poradi komponenty v poli history parametru L,M,N,sigN_fib,EpsN_fib
    previousSigmaL=Mm->ip[ipp].eqother[i_other];
    previousSigmaM=Mm->ip[ipp].eqother[i_other+1];
    previousSigmaN=Mm->ip[ipp].eqother[i_other+2];
    previousSigmaD=previousSigmaN-previousSigmaV;

    previousSigmaN_fiber=Mm->ip[ipp].eqother[i_other+3];
    previousEpsN=Mm->ip[ipp].eqother[i_other+4];

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
    sigmaN_concrete=minim( (sigmaV+sigmaD) ,FN(epsN,previousSigmaV));

    //napeti ve vlakne:
    fiber_module=(previousEpsN!=0.0 ? (previousSigmaN_fiber/previousEpsN) : 0.0 );
    senf=previousSigmaN_fiber+fiber_module*depsN;

   
    /*if(depsN<0) {
if (ipp==1) printf("mp=%d : depsN=%ef : senf=%lf : FF=%lf\n",mPlaneIndex,depsN,senf,FFiber(macbra(epsN)));
}*/

    sigmaN_fiber=(depsN<0.0 ? senf : FFiber(macbra(epsN)));
    
    //// sigmaN_fiber=macbra( minim(senf, FFiber(macbra(epsN))));
    //sigmaN_fiber= FFiber(macbra(epsN));
    // sigmaN=sigmaN_concrete+ sigmaN_fiber;

    //shear microstresses
    sel=previousSigmaL + et*depsL;
    sem=previousSigmaM + et*depsM;    
    f=FT(sigmaN_concrete,epsV);
    if (sel>f) sigmaL=f; else if (sel<-f) sigmaL=-f; else sigmaL=sel;  
    if (sem>f) sigmaM=f; else if (sem<-f) sigmaM=-f; else sigmaM=sem;
  
    //mean normal microstress (computed after sweeping throgh all microplanes)
    meanN += sigmaN_concrete* microplaneWeights[mPlaneIndex]*6.0;
      
    //UPDATE history variables, i.e. store microstresses
    Mm->ip[ipp].other[i_other]=sigmaL;
    Mm->ip[ipp].other[i_other+1]=sigmaM;
    Mm->ip[ipp].other[i_other+2]=sigmaN_concrete;
    Mm->ip[ipp].other[i_other+3]=sigmaN_fiber;
    Mm->ip[ipp].other[i_other+4]=epsN;
    
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
    i_other=ido+(mPlaneIndex+1)*5-4;//poradi komponenty v poli history parametru L,M,N
    sigmaL=Mm->ip[ipp].other[i_other];
    sigmaM=Mm->ip[ipp].other[i_other+1];
    sigmaN_concrete=Mm->ip[ipp].other[i_other+2];
    //recalculation of sigmaD according to possibly new sigmaV 
    sigmaD_concrete = sigmaN_concrete - sigmaV;
    sigmaN_fiber=Mm->ip[ipp].other[i_other+3];
    
    
    
      for (i=0; i<6; i++) {
    
	  macroSigma[i] +=((projN[mPlaneIndex][i]-kronecker[i]/3.)*sigmaD_concrete + 
		       projN[mPlaneIndex][i]*sigmaN_fiber +
		       projL[mPlaneIndex][i]*sigmaL +
		       projM[mPlaneIndex][i]*sigmaM)* microplaneWeights[mPlaneIndex]*6.0;
		       }
	

      /*  macroSigma[i] += ((projN[mPlaneIndex][i]-kronecker[i]/3.)*sigmaD_concrete + 
	  projN[mPlaneIndex][i]*sigmaN_fiber)* microplaneWeights[mPlaneIndex]*6.0; */
      


    
  }//end 2nd loop over nmp
  
  //2nd constraint, addition of volumetric part
  macroSigma[0]+=sigmaV;
  macroSigma[1]+=sigmaV;
  macroSigma[2]+=sigmaV;
  
  //send new macrostress vector to the solver
  Mm->storestress(0,ipp,macroSigma);

}

void microfiber::updateval(long ipp,long im,long ido)
{
  long i,n = Mm->givencompother(ipp, im);

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }
}

inline double
microfiber::macbra(double x) /* Macauley bracket = positive part of x */
{
  return(maxim(x,0.0));
}

inline double
microfiber::FVplus (double epsV)
/*positive volumetric boundary */
{
  return(ev*k1*c13/(1+(c14/k1)*macbra(epsV-c13*c15*k1)));
}

inline double
microfiber::FVminus (double epsV)
/*negative volumetric boundary */
{
return(-e*k1*k3*exp(-epsV/(k1*k4)));
}

inline double
microfiber::FDminus(double epsD)
 /*negative deviatoric boundary */
 {
   double a;
   a=macbra(-epsD-c8*c9*k1)/(k1*c7);
   return(-e*k1*c8/(1+a*a));
 }

inline double
microfiber::FDplus(double epsD)
/*positive deviatoric bondary */
{
  double a;
  a=(macbra(epsD-c5*c6*k1)/(k1*c20*c7));
  return (e*k1*c5/(1+a*a));
}

inline double
microfiber::FN(double epsN,double sigmaV)
 /*normal boundary */
{
  return(e*k1*c1*exp(-macbra(epsN-c1*c2*k1)/(k1*c3+macbra(-c4*(sigmaV/ev)))));
}

inline double
microfiber::FT (double sigmaN,double epsV)
/*shear boundary */
{
  double a, sn0;
  sn0= et*k1*c11/(1+c12*fabs(epsV));
  a=macbra(-sigmaN+sn0);
  return(et*k1*k2*c10*a/(et*k1*k2+c10*a));
}

inline double
microfiber::FFiber (double epsN)
/*fiber boundary */
{
  double f1, f2, a, x;
  x=epsN-k10;
  f1=e*k13*(2*sqrt(x/k11) - x/k11);
  a=1-(epsN-k10-k11)/k12;
  f2=e*k13*a*a;

  if ((epsN<k10) || (epsN>(k10+k11+k12)) ) return(0.0);
  else { 
    if( epsN<=(k10+k11) ) return (f1);
    else return (f2);
  }
}

inline double
microfiber::maxim (double a,double b)
{
  return(a>b ? a:b);
}

inline double
microfiber::minim (double a,double b)
{
  return(a<b ? a:b);
}
