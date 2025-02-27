#include "stdlib.h"
#include "math.h"
#include "homogt.h"
#include "matrix.h"
#include "vector.h"


//all homogenization methods provide a set of data E-nu or k-mu

/*symmetric means to be invariant if MTRX is swapped with AGGR*/
HomogData::HomogData(const double iE_m, const double inu_m, const double if_m,
		     const double iE_i, const double inu_i, const double iHirsch)

{
  E_m = iE_m; nu_m = inu_m; f_m = if_m;
  E_i = iE_i; nu_i = inu_i; f_i = (1-if_m);
  H_chi=iHirsch;
       
  ENuToKMu(E_m, nu_m, k_m, mu_m);
  ENuToKMu(E_i, nu_i, k_i, mu_i);
}

/*constructor for two phases of PhaseMatrix*/
HomogData::HomogData(double *PhaseMatrix){
  f_m=*(PhaseMatrix+0);
  E_m=*(PhaseMatrix+1);
  nu_m=*(PhaseMatrix+2);
  f_i=*(PhaseMatrix+3);
  E_i=*(PhaseMatrix+4);
  nu_i=*(PhaseMatrix+5);

  if((f_m+f_i)<0.99 || (f_m+f_i)>1.01){
    printf("volumetric fraction of phases 0,1 is not unity, terminating\n");
    exit(1);
  }

  ENuToKMu(E_m, nu_m, k_m, mu_m);
  ENuToKMu(E_i, nu_i, k_i, mu_i);
}



HomogData::HomogData()
{
}

/*parallel model by Voigt, matrix and inclusions values are symmetric*/
void HomogData::Voigt(void) {

  k_hmg=k_m*f_m+k_i*f_i;
  mu_hmg=mu_m*f_m+mu_i*f_i;

  KMuToENu(k_hmg, mu_hmg, E_hmg, nu_hmg);
}

/*serial model by Reuss, matrix and inclusions values are symmetric*/
void HomogData::Reuss(void){

  k_hmg=1/(f_m/k_m+f_i/k_i); 
  mu_hmg=1/(f_m/mu_m+f_i/mu_i); 
  KMuToENu(k_hmg, mu_hmg, E_hmg, nu_hmg);
}


/*Hirsh model as a combination of Voigt and Reuss (coef. chi)*/
/*standard value of 0.5 for concrete is suggested*/
void HomogData::Hirsch(void){
  E_hmg=1/((1.-H_chi)*(f_m/E_m+f_i/E_i)+H_chi/(f_m*E_m+f_i*E_i));
  /*nu_hmg unknown*/
}


/*Hansen model spherical aggregate in spherical paste, for Poisson's number=0.2*/
/*model is not symmetric*/
void HomogData::Hansen(void){
  E_hmg=E_i*((f_i*E_i+(1+f_m)*E_m)/((1.+f_m)*E_i+f_i*E_m));
  /*nu_hmg unknown*/
}


/*Counto model a prismatic aggregate in a prismatic paste*/
/*model is not symmetric*/
void HomogData::Counto(void){

  E_hmg=1./((1.-sqrt(f_m))/E_i+1./(E_i*((1.-sqrt(f_m))/sqrt(f_m))+E_m));
  /*nu_hmg unknown*/
}


/*Mori-Tanaka homogenization method*/
/*model is not symmetric*/
/*for derivation see Bernard et al.:A multiscale micromechanics-hydration model... CCR,pp. 1293-1309, 2003*/
void HomogData::MoriTanaka(void){
  double t, k_MT, mu_MT;

  alpha_m = 3. * k_m / ( 3. * k_m + 4 * mu_m );
  beta_m  = 6. * ( k_m + 2. * mu_m ) / ( 5. * ( 3. * k_m + 4. * mu_m ) );

  t = 1. + alpha_m * ( k_i / k_m - 1. );

  k_MT  = f_m * k_m + f_i * k_i / t;
  k_MT /= f_m + f_i / t + ( 1. - f_m - f_i ) / ( 1. - alpha_m );

  t = 1. + beta_m * ( mu_i / mu_m - 1. );

  mu_MT  = f_m * mu_m + f_i * mu_i / t;
  mu_MT /= f_m + f_i / t + ( 1. - f_m - f_i ) / ( 1. - beta_m );

  KMuToENu( k_MT, mu_MT, E_hmg, nu_hmg );

  k_hmg=k_MT;
  mu_hmg=mu_MT;
}

/*
Mori-Tanaka homogenization method
model is not symmetric
matrix phase is dominant and inclusions are in spherical form
scheme applies for at least two phases - information are in PhaseMatrix
Phasematrix[][0] is vol. fraction of a given phase
Phasematrix[][1] is intrinsic E modulus
Phasematrix[][2] is intrinsic Poisson's ratio
NumRows specifies, how many rows should be read from the matrix
MatrixRow is the reference phase of this scheme (0 is the first phase)
for derivation see Bernard et al.:A multiscale micromechanics-hydration model... CCR,pp. 1293-1309, 2003
*/

void HomogData::MT_mtrx(double *PhaseMatrix, int NumRows, int MatrixRow){
  double f_r, E_r, nu_r, k_r, mu_r;
  double fr_tot = 0.;
  double nom_k_MT=0. , denom_k_MT=0. , nom_mu_MT=0., denom_mu_MT=0.;
  //double t, k_MT, mu_MT;

  //determine matrix parameters 
  f_m=*(PhaseMatrix+3*MatrixRow+0);
  E_m=*(PhaseMatrix+3*MatrixRow+1); 
  nu_m=*(PhaseMatrix+3*MatrixRow+2);

  ENuToKMu(E_m, nu_m, k_m, mu_m);

  //auxiliary parameters
  alpha_m = 3. * k_m / ( 3. * k_m + 4 * mu_m );
  beta_m  = 6. * ( k_m + 2. * mu_m ) / ( 5. * ( 3. * k_m + 4. * mu_m ) );

  //cycle through all phases, including matrix phase
  for (int r=0; r<NumRows; r++){
    f_r=*(PhaseMatrix+3*r+0); 
    E_r=*(PhaseMatrix+3*r+1);
    nu_r=*(PhaseMatrix+3*r+2);
    
    fr_tot +=f_r;
    ENuToKMu(E_r, nu_r, k_r, mu_r);

    nom_k_MT += f_r*k_r/(1.+alpha_m*(k_r/k_m-1.));
    denom_k_MT += f_r/(1.+alpha_m*(k_r/k_m-1.));

    nom_mu_MT += f_r*mu_r/(1.+beta_m*(mu_r/mu_m-1.));
    denom_mu_MT += f_r/(1.+beta_m*(mu_r/mu_m-1.));
  }

  if(fr_tot<0.99 || fr_tot>1.01){
    printf("volumetric fraction of phases 0-%d is not unity (= %lf), terminating (MT_mtrx)\n",  NumRows, fr_tot);
    exit(1);
  }
  
  k_hmg=nom_k_MT/denom_k_MT;
  mu_hmg=nom_mu_MT/denom_mu_MT;
  
  KMuToENu(k_hmg, mu_hmg, E_hmg, nu_hmg);
  
}



/*
Self-consistent scheme of Hill and Budiansky
model is not symmetric, only for two phases
for derivation see Bernard et al.:A multiscale micromechanics-hydration model... CCR,pp. 1293-1309, 2003*/
void HomogData::SelfConsistent(void){
  double k_SCS, mu_SCS, t_m, t_i;
 
  /*estimation for the beginning of iteration*/  
  k_SCS=10;
  mu_SCS=.3;

  /*iteration for nonlinear equations*/
  for(int i=1; i<100; i++){
  
    alpha_m=3.*k_SCS/(3.*k_SCS+4*mu_SCS);
    beta_m=6.*(k_SCS+2.*mu_SCS)/(5.*(3.*k_SCS+4.*mu_SCS));
  
    t_m=1.+alpha_m*(k_m/k_SCS-1.);
    t_i=1.+alpha_m*(k_i/k_SCS-1.);

    k_SCS=f_m*k_m/t_m+f_i*k_i/t_i;
    k_SCS/=(f_m/t_m+f_i/t_i);

    /*k_SCS is updated*/
    t_m=1.+beta_m*(mu_m/mu_SCS-1.);
    t_i=1.+beta_m*(mu_i/mu_SCS-1.);

    mu_SCS=f_m*mu_m/t_m+f_i*mu_i/t_i;
    mu_SCS/=(f_m/t_m+f_i/t_i);
  }

  KMuToENu( k_SCS, mu_SCS, E_hmg, nu_hmg ); 
  
  k_hmg=k_SCS;
  mu_hmg=mu_SCS;

  temp1=alpha_m;
  temp2=beta_m;
  temp3=k_SCS;
  temp4=mu_SCS;
}


/*
Self-consistent scheme for many phases of Hill and Budiansky
model is not symmetric
for derivation see Bernard et al.:A multiscale micromechanics-hydration model... CCR,pp. 1293-1309, 2003
the meaning of PhaseMatrix and NumRows applies as in Mori-Tanaka scheme
NumRows is >=1
*/
void HomogData::SCS(double *PhaseMatrix, int NumRows){
  double f_r, E_r, nu_r, k_r, mu_r;
  double k_SCS, mu_SCS, nom_k_SCS, denom_k_SCS, nom_mu_SCS, denom_mu_SCS;
  double fr_tot;
  
  /*first estimation*/
  k_SCS=10.;
  mu_SCS=.3;

  /*iteration for nonlinear equations*/
  for(int i=1; i<100; i++){
    fr_tot = 0.; 
    nom_k_SCS=0;
    denom_k_SCS=0;
    nom_mu_SCS=0;
    denom_mu_SCS=0;
    for (int r=0; r<NumRows; r++){
       f_r=*(PhaseMatrix+3*r+0); 
      E_r=*(PhaseMatrix+3*r+1); 
      nu_r=*(PhaseMatrix+3*r+2);
 
      fr_tot += f_r;
     
      ENuToKMu(E_r, nu_r, k_r, mu_r);

      alpha_m=3.*k_SCS/(3.*k_SCS+4*mu_SCS);
      beta_m=6.*(k_SCS+2.*mu_SCS)/(5.*(3.*k_SCS+4.*mu_SCS));
  
      nom_k_SCS+=f_r*k_r/(1+alpha_m*(k_r/k_SCS-1));
      denom_k_SCS+=f_r/(1+alpha_m*(k_r/k_SCS-1));

      nom_mu_SCS+=f_r*mu_r/(1+beta_m*(mu_r/mu_SCS-1));
      denom_mu_SCS+=f_r/(1+beta_m*(mu_r/mu_SCS-1));
    }
    k_SCS=nom_k_SCS/denom_k_SCS;
    mu_SCS=nom_mu_SCS/denom_mu_SCS;
  }

  temp3=k_SCS;
  temp4=mu_SCS;

  if(fr_tot<0.99 || fr_tot>1.01){
    printf("volumetric fraction of phases 0-%d is not unity (= %lf), terminating (SCS)\n", NumRows, fr_tot);
    exit(1);
  }

  KMuToENu(k_SCS, mu_SCS, E_hmg, nu_hmg ); 

  k_hmg=k_SCS;
  mu_hmg=mu_SCS;
}


//parallel configuration
void HomogData::Voigt(double *PhaseMatrix, int NumRows){

  E_hmg = 0.;
  nu_hmg = 0.;
  for (int r=0; r<NumRows; r++){
    E_hmg+=(*(PhaseMatrix+3*r+0))*(*(PhaseMatrix+3*r+1));
    nu_hmg+=(*(PhaseMatrix+3*r+0))*(*(PhaseMatrix+3*r+2));
  }
  ENuToKMu(E_hmg, nu_hmg, k_hmg, mu_hmg); 
}

//serial configuration
void HomogData::Reuss(double *PhaseMatrix, int NumRows){

  E_hmg = 0.;
  nu_hmg = 0.;
  for (int r=0; r<NumRows; r++){
    E_hmg+=(*(PhaseMatrix+3*r+0))/(*(PhaseMatrix+3*r+1));
    nu_hmg+=(*(PhaseMatrix+3*r+0))/(*(PhaseMatrix+3*r+2));
  }
  E_hmg=1/E_hmg;
  nu_hmg=1/nu_hmg;
  ENuToKMu(E_hmg, nu_hmg, k_hmg, mu_hmg); 
}


/*Hashin-Shtrikman lower and upper bounds, only 2 phases
model is symmetric - then upper and lower bounds swap
k_m<k_i and mu_m<mu_i or
k_m>k_i and mu_m>mu_i
if not satisfied, e.g. Walpole bounds must be used
*/
void HomogData::HashinShtrikman(void){
  double k_HS, mu_HS;

  if (((k_m < k_i) && (mu_m > mu_i)) || ((k_m > k_i) && (mu_m < mu_i))){
    printf("Non well ordered two-phase material, use Walpole bounds\n");
    printf("k[0]=%f k[1]=%f mu[0]=%f mu[1]=%f\n", k_m, k_i, mu_m, mu_i);
    exit(0);
  }

  /*lower bounds if k_m<k_i && mu_m<mu_i*/
  k_HS=k_m+f_i/(1./(k_i-k_m)+3.*f_m/(3.*k_m+4.*mu_m));
  mu_HS=mu_m+f_i/(1./(mu_i-mu_m)+6.*f_m*(k_m+2.*mu_m)/(5.*mu_m*(3.*k_m+4.*mu_m)));


  //experimental-is different - why?? Sigmund:A new class of extremal composites, JMechSolPhy 397-428, 2000
  //  mu_HS=mu_m+f_i/(1./(mu_i-mu_m)+(1-f_i)*(k_m+2.*mu_m)/(2.*mu_m*(k_m+mu_m)));


  KMuToENu(k_HS, mu_HS, E_hmg, nu_hmg);

  k_hmg=k_HS;
  mu_hmg=mu_HS;
  
  /*upper bounds if k_m<k_i && mu_m<mu_i*/
  k_HS=k_i+f_m/(1./(k_m-k_i)+3.*f_i/(3.*k_i+4.*mu_i));

  mu_HS=mu_i+f_m/(1./(mu_m-mu_i)+6.*f_i*(k_i+2.*mu_i)/(5.*mu_i*(3.*k_i+4.*mu_i)));

  KMuToENu(k_HS, mu_HS, E_hmg_2, nu_hmg_2);

  k_hmg_2=k_HS;
  mu_hmg_2=mu_HS;

 
  /*exchange lower and upper bounds*/
 
  if (E_hmg>E_hmg_2){
    Swap(E_hmg,E_hmg_2);
  }
  if (nu_hmg>nu_hmg_2){
    Swap(nu_hmg,nu_hmg_2);
  }

  if (k_hmg>k_hmg_2){
    Swap(k_hmg,k_hmg_2);
  }
  if (mu_hmg>mu_hmg_2){
    Swap(mu_hmg,mu_hmg_2);
  }
}


/*
Hashin-Shtrikman-Walpole lower and upper bounds, arbitrary # of phases
phases do not have to satisfy k_0<k_1<...k_n together with mu_0<mu_1<...mu_n
model is symmetric - then upper and lower bound swap
see Berryman: Mixture theories for rock properties - www.agu.org/reference/rock/16_berryman.pdf
the meaning of PhaseMatrix and NumRows applies as in Mori-Tanaka scheme
*/

void HomogData::WalpoleMulti(double *PhaseMatrix, int NumRows){
  double k_low, mu_low, k_high, mu_high;
  double k_min, k_max, mu_min, mu_max;
  double **PhaseMatrixKMu;
  PhaseMatrixKMu = new double*[NumRows];
  long i;
  for(i=0; i<NumRows; i++)
	 PhaseMatrixKMu[i]= new double[3];
  
 int r;
  for(r=0; r<NumRows; r++){
    //material fraction is copied
    PhaseMatrixKMu[r][0]=*(PhaseMatrix+3*r+0);
    //K and nu is recalculated from E, nu
    ENuToKMu(*(PhaseMatrix+3*r+1), *(PhaseMatrix+3*r+2), PhaseMatrixKMu[r][1], PhaseMatrixKMu[r][2]);
  }

 
  //find maximum and minimum values
  k_min=PhaseMatrixKMu[0][1];
  k_max=0;
  mu_min=PhaseMatrixKMu[0][2];
  mu_max=0;
  for(r=0; r<NumRows; r++){
    if(PhaseMatrixKMu[r][1]<k_min)
      k_min=PhaseMatrixKMu[r][1];
    if(PhaseMatrixKMu[r][2]<mu_min)
      mu_min=PhaseMatrixKMu[r][2];
    if(PhaseMatrixKMu[r][1]>k_max)
      k_max=PhaseMatrixKMu[r][1];
    if(PhaseMatrixKMu[r][2]>mu_max)
      mu_max=PhaseMatrixKMu[r][2];
  }

  //calculate bounds
  k_low=Lambda(*PhaseMatrixKMu, NumRows, k_min);
  k_high=Lambda(*PhaseMatrixKMu, NumRows, k_max);
  mu_low=Gamma(*PhaseMatrixKMu, NumRows, Zeta(k_min, mu_min));
  mu_high=Gamma(*PhaseMatrixKMu, NumRows, Zeta(k_max, mu_max));

  //convert to E, nu
  KMuToENu(k_low, mu_low, E_hmg, nu_hmg);
  KMuToENu(k_high, mu_high, E_hmg_2, nu_hmg_2);

  k_hmg=k_low;
  mu_hmg=mu_low;
  k_hmg_2=k_high;
  mu_hmg_2=mu_high;
}


//help function for bulk modulus, mu refers to min or max value
double HomogData::Lambda(double *PhaseMatrixKMu, int NumRows, double mu){
  double lambda=0;
  for (int r=0; r<NumRows; r++){
    lambda+=*(PhaseMatrixKMu+3*r+0)/(*(PhaseMatrixKMu+3*r+1)+4./3.*mu);
  }
  lambda=1./lambda;
  lambda-=4./3.*mu;
  return lambda;
}

//help function for shear modulus, k refers to min or max value
double HomogData::Gamma(double *PhaseMatrixKMu, int NumRows, double k){
  double gamma=0;
  for (int r=0; r<NumRows; r++){
    gamma+=*(PhaseMatrixKMu+3*r+0)/(*(PhaseMatrixKMu+3*r+2)+k);
  }
  gamma=1./gamma;
  gamma-=k;
  return gamma;
}

//help function
double HomogData::Zeta(double k, double mu){
  return mu/6.*(9.*k+8.*mu)/(k+2.*mu);
}


/*
Homogenization scheme of Herve and Zaoui for n-spherical isotropic domains
see Herve E., Zaoui A. n-layered Inclusion-based Micromechanical Modelling. 
Int. J. Engng Sci. 31, pp 1-10, 1993.

The sum of n phases must yield a unity, (n+1) may have arbitrary radius

The scheme is derived from stress and displacement equivalence at boundaries among phases
homogenized bulk and shear response is calculated
NumPhases is the number of phases including homogeneous medium (n+1)
In homogenized values, (n+1) phase may have arbitrary values - if it has the right 
homogenized values of the assembly, the equations (49) and (51) are satisfied exactly without solution

In case of n=2, the Hashin-Shtrikman bounds are calculated if the well-ordered moduli are used 
Homogenized bulk modulus is 100 % OK (correspond to HS bounds for 2 phases)
Homogenized shear modulus is probably not 100 % OK (does not correspond to HS bounds for 2 phases)
The problem is not in numerics, neither solver, neither matrices M

Also, when one phase is formally replaced for two phases with the same vol.
content and properties as the replaced phase, the result are the same

Shear modulus is with small difference when compared to HS bound for 2 phases
typical usage for mortars and concrete is
phase 0-aggregate
phase 1-ITZ
phase 2 (n) -cement paste
phase 3 (n+1) -equivalent conrete media
*/

void HomogData::HerveZaoui(double *PhaseMatrix, int NumPhases){
  vector r(NumPhases-1),k(NumPhases),mu(NumPhases);
  //  int phase;
  //variables for HYDROSTATIC pressure
  vector Q11(NumPhases-1), Q21(NumPhases-1), F(NumPhases), G(NumPhases);
  matrix J(2,2), Jinv(2,2), N(2,2), Nhelp(2,2), Q(2,2);
  //set lambda_0 for strain control - arbitrary value for homogenization
  double lambda_0=0.735; 
  //variables for simple SHEAR
  vector A(NumPhases), B(NumPhases), C(NumPhases), D(NumPhases);
  vector P_v(4), Phelp(4);
  matrix L(4,4), Linv(4,4), M(4,4), Mhelp(4,4), Z(4,4), M_test(4,4);
  matrix *P_arr;//array to store individual P matrices
  //set gamma for displacement approach - arbitrary value for homogenization
  double gamma=2.;
  //variables for homogenized values
  long double a,b,c,nu_n,sqr;

  long double ak, bk, ck, dk, ek, fk, alphak, nuk, nukp1, koef;
  int i, phi;

  if(NumPhases<3){
    printf("Number of considered phases must be at least 3 (including hom. medium)\n");
    exit(1);
  }
  
  F[NumPhases-1]=lambda_0/3.;
  A[NumPhases-1]=gamma;

  //radii of spheres, M_PI defined in cmath, assume 1 as the maximal radius
  temp1=0.;
  for(i=0; i<NumPhases-1;i++){
    temp1+=*(PhaseMatrix+3*i+0);
    r[i]=1.*pow(temp1,1./3.);
    //    printf("r[%d]=%f\n", i, r[i]);
  }

  if(r[NumPhases-2]<0.99 || r[NumPhases-2]>1.01){
    printf("volumetric fraction of phases 0-%d is not unity, terminating\n", NumPhases-2);
    exit(1);
  }

  
  for(i=0; i<NumPhases;i++){
    ENuToKMu(*(PhaseMatrix+3*i+1), *(PhaseMatrix+3*i+2), k[i], mu[i]);
    //printf("k[%d]=%f mu[%d]=%f\n", i, k[i], i, mu[i]);
  }
  
  /*a)hydrostatic pressure - conditions of displacement and 
    radial stress compatibility between adjacent phases
  */

  //create ones matrix for the first multiplication
  Nhelp[0][0]=1.;
  Nhelp[1][1]=1.;
  
  //calculate Q
  for(phi=0;phi<=NumPhases-2;phi++){
    //calculate inverse of J_{phi+1}
    fillJ(J,r[phi],mu,k,phi+1);
    invm(J, Jinv, 1.e-10);
    //calculate J_{phi}
    fillJ(J,r[phi],mu,k,phi);
    //calulate N^{phi}=J_{phi+1}^{-1}*J_{phi}      
    mxm(Jinv,J, N);
    //calculate a part of Q^{phi}
    mxm(N,Nhelp,Q);
    //store part of matrix Q in vector Q11 and Q21
    Q11[phi]=Q[0][0];
    Q21[phi]=Q[1][0];
    //copy Q to Nhelp and use in the next cycle
    copym(Q, Nhelp);
  }

  //calculate V (contain vector components F,G), F[0] and G[0] make no sense
  for(phi=1;phi<=NumPhases-1;phi++){
    F[phi]=F[NumPhases-1]*Q11[phi-1]/Q11[NumPhases-2];
    G[phi]=F[NumPhases-1]*Q21[phi-1]/Q11[NumPhases-2];
    //printf("F[%d]=%f G[%d]=%f\n", phi, F[phi], phi, G[phi]); 
  }
  
  //check displacement and stress continuity, phi is the more distant phase
  for(phi=1;phi<=NumPhases-1;phi++){
    //printf("r=%f displ(in)=%f displ(out)=%f\n", r[phi-1], F[phi-1]*r[phi-1]+G[phi-1]/(r[phi-1]*r[phi-1]), F[phi]*r[phi-1]+G[phi]/(r[phi-1]*r[phi-1])); 
    //printf("r=%f stress(in)=%f stress(out)=%f\n", r[phi-1], 3.*k[phi-1]*F[phi-1]-4.*mu[phi-1]*G[phi-1]/(r[phi-1]*r[phi-1]*r[phi-1]), 3.*k[phi]*F[phi]-4.*mu[phi]*G[phi]/(r[phi-1]*r[phi-1]*r[phi-1]));
  }
  
  /*replace n+1 layer with equivalent homogeneous media
  calculate homogenized k
  */
  k_hmg=(3.*k[NumPhases-2]*pow(r[NumPhases-2],3)*Q11[NumPhases-3]-
	4.*mu[NumPhases-2]*Q21[NumPhases-3])/
    (3.*(pow(r[NumPhases-2],3)*Q11[NumPhases-3]+Q21[NumPhases-3]));
  //printf("Homogenized k=%f\n", k_hmg);

  destrv(Q11);
  destrv(Q21);
  destrv(F);
  destrv(G);

  destrm(J);
  destrm(Jinv);
  destrm(N);
  destrm(Nhelp);
  destrm(Q);



  /*b)simple shear - continuity of radial and phi displacement and
    two shear stresses
  */
  //create ones matrix for the first multiplication
  Mhelp[0][0]=1.;
  Mhelp[1][1]=1.;
  Mhelp[2][2]=1.;
  Mhelp[3][3]=1.;

  //allocate memory for P matrices, store there P11, P12, P21, P22
  P_arr=new matrix[NumPhases-1];

  //calculate P
  for(phi=0;phi<=NumPhases-2;phi++){
    //calculate inverse of L_{phi+1}
    fillL(L,r[phi],mu,k,phi+1);
    invm(L, Linv, 1.e-10);
    //calculate L_{phi}
    fillL(L,r[phi],mu,k,phi);
    //calulate M^{phi}=L_{phi+1}^{-1}*L_{phi}      
    mxm(Linv,L,M);//multiplication order OK
    //test routine directly for M
    nuk=*(PhaseMatrix+3*phi+2);
    nukp1=*(PhaseMatrix+3*(phi+1)+2);

    //this part calculate inverse of M directly
    //M=M_test so the L, Linv and M assembly are OK
    ak=mu[phi]/mu[phi+1]*(7.+5.*nuk)*(7.-10.*nukp1)-(7.-10*nuk)*(7.+5.*nukp1);
    bk=4.*(7.-10.*nuk)+mu[phi]/mu[phi+1]*(7.+5.*nuk);
    ck=(7.-5.*nukp1)+2.*mu[phi]/mu[phi+1]*(4.-5.*nukp1);
    dk=(7.+5.*nukp1)+4.*mu[phi]/mu[phi+1]*(7.-10.*nukp1);
    ek=2*(4.-5.*nuk)+mu[phi]/mu[phi+1]*(7.-5.*nuk);
    fk=(4.-5.*nuk)*(7.-5.*nukp1)-mu[phi]/mu[phi+1]*(4.-5.*nukp1)*(7.-5.*nuk);
    alphak=mu[phi]/mu[phi+1]-1.;
    koef=1./(5.*(1-nukp1));

    M_test[0][0]=koef*ck/3.;
    M_test[0][1]=koef*pow(r[phi],2)*(3*bk-7*ck)/(5.*(1.-2*nuk));
    M_test[0][2]=koef*(-12.)*alphak/pow(r[phi],5);
    M_test[0][3]=koef*4.*(fk-27.*alphak)/(15.*(1.-2*nuk)*pow(r[phi],3));

    M_test[1][0]=koef*0.;
    M_test[1][1]=koef*(1.-2*nukp1)*bk/(7.*(1.-2.*nuk));
    M_test[1][2]=koef*(-20.)*(1.-2.*nukp1)*alphak/(7*pow(r[phi],7));
    M_test[1][3]=koef*(-12.)*alphak*(1.-2.*nukp1)/(7.*(1.-2.*nuk)*(pow(r[phi],5)));

    M_test[2][0]=koef*pow(r[phi],5)*alphak/2.;
    M_test[2][1]=koef*(-pow(r[phi],7))*(2*ak+147*alphak)/(70.*(1.-2.*nuk));
    M_test[2][2]=koef*dk/7.;
    M_test[2][3]=koef*pow(r[phi],2)*(105*(1.-nukp1)+12.*alphak*(7.-10.*nukp1)-7.*ek)/(35.*(1.-2*nuk));

    M_test[3][0]=koef*(-5.)/6.*(1.-2*nukp1)*alphak*pow(r[phi],3);
    M_test[3][1]=koef*7.*(1.-2*nukp1)*alphak*pow(r[phi],5)/(2.*(1.-2.*nuk));
    M_test[3][2]=koef*0.;
    M_test[3][3]=koef*ek*(1.-2.*nukp1)/(3.*(1.-2*nuk));
 
    /*
    for(int ii=0;ii<4;ii++){
      for(int jj=0;jj<4;jj++){
	M[ii][jj]=M_test[ii][jj];
      }
    }
    */


   //allocate matrix P[]
    allocm(4,4,P_arr[phi]);
    //calculate a part of P[phi]^{phi}
    mxm(M,Mhelp,P_arr[phi]);//M*M_help=M_help*M, but why?????
    //copy P to Mhelp and use in the next cycle
    copym(P_arr[phi], Mhelp);
  }

  /*
  P11, P12, P21, P22 is contained in matrix array P_arr[].a[]
  P11[phi]=P_arr[phi].a[0+0];
  P12[phi]=P_arr[phi].a[0+1];
  P21[phi]=P_arr[phi].a[4+0];
  P22[phi]=P_arr[phi].a[4+1];
  */
  //calculate vector W (contain scalar components A,B,C,D for each step)

  P_v[0]=P_arr[NumPhases-2].a[4+1];
  P_v[1]=-P_arr[NumPhases-2].a[4+0];
  P_v[2]=P_v[3]=0;
  for(phi=1;phi<=NumPhases-1;phi++){
    mxv(P_arr[phi-1], P_v, Phelp);
    cmulv(A[NumPhases-1]/(P_arr[NumPhases-2].a[4+1]*P_arr[NumPhases-2].a[0+0]-P_arr[NumPhases-2].a[0+1]*P_arr[NumPhases-2].a[4+0]),Phelp, Phelp);
    A[phi]=Phelp[0];
    B[phi]=Phelp[1];
    C[phi]=Phelp[2];
    D[phi]=Phelp[3];
    //printf("A[%d]=%f B[%d]=%f C[%d]=%f D[%d]=%f\n", phi, A[phi], phi, B[phi], phi, C[phi], phi, D[phi]); 
  }

  //check displacement and stress continuity, phi is the less distant phase
  for(phi=0;phi<=NumPhases-2;phi++){
    //phase closer to the center
    fillL(L,r[phi],mu,k,phi);
    Phelp[0]=A[phi];
    Phelp[1]=B[phi];
    Phelp[2]=C[phi];
    Phelp[3]=D[phi];
    mxv(L, Phelp, P_v);
    //printf("(in) n=%d r=%f ur=%f uphi=%f srp=%f srf=%f\n",phi, r[phi], P_v[0],P_v[1], P_v[2], P_v[3]);

    fillL(L,r[phi],mu,k,phi+1);
    Phelp[0]=A[phi+1];
    Phelp[1]=B[phi+1];
    Phelp[2]=C[phi+1];
    Phelp[3]=D[phi+1];
    mxv(L, Phelp, P_v);
    //printf("(out) n=%d r=%f ur=%f uphi=%f srp=%f srf=%f\n", phi,r[phi], P_v[0],P_v[1], P_v[2], P_v[3]);
  }

  /*replace n+1 layer with equivalent homogeneous media
  calculate homogenized mu
  */
  for(i=0;i<=3;i++){
    for(int j=0;j<=3;j++){
      Z[i][j]=P_arr[NumPhases-3].a[4*i+0]*P_arr[NumPhases-3].a[4*j+1]-
	P_arr[NumPhases-3].a[4*j+0]*P_arr[NumPhases-3].a[4*i+1];
      //      printf("Z=%f ", Z[i][j]);
    }
  }
  //from this point the equilibrium equations are OK, error is somewhere above for the shear modulus
  //equilibrium equation
  //printf("ZERO_P %f\n", P_arr[NumPhases-2].a[12+0]*P_arr[NumPhases-2].a[4+1]-P_arr[NumPhases-2].a[12+1]*P_arr[NumPhases-2].a[4+0]);

  //printf("ZERO_MZ %f\n", M[3][0]*M[1][1]*Z[0][1]+M[3][0]*M[1][2]*Z[0][2]+M[3][0]*M[1][3]*Z[0][3]+M[3][1]*M[1][2]*Z[1][2]+M[3][3]*M[1][2]*Z[3][2]+(M[3][1]*M[1][3]-M[3][3]*M[1][1])*Z[1][3]);


  nu_n=(3.*k[NumPhases-2]-2.*mu[NumPhases-2])/(6.*k[NumPhases-2]+2.*mu[NumPhases-2]);
  
  a=4.*pow(r[NumPhases-2],10)*(1.-2.*nu_n)*(7.-10.*nu_n)*Z[0][1]+
    20.*pow(r[NumPhases-2],7)*(7.-12.*nu_n+8.*nu_n*nu_n)*Z[3][1]+
    12.*pow(r[NumPhases-2],5)*(1.-2.*nu_n)*(Z[0][3]-7.*Z[1][2])+
    20.*pow(r[NumPhases-2],3)*(1.-2.*nu_n)*(1.-2.*nu_n)*Z[0][2]+
    16.*(4.-5.*nu_n)*(1.-2.*nu_n)*Z[3][2];
  
  b=3.*pow(r[NumPhases-2],10)*(1.-2.*nu_n)*(15.*nu_n-7.)*Z[0][1]+
    60.*pow(r[NumPhases-2],7)*(nu_n-3.)*nu_n*Z[3][1]-
    24.*pow(r[NumPhases-2],5)*(1.-2.*nu_n)*(Z[0][3]-7.*Z[1][2])-
    40.*pow(r[NumPhases-2],3)*(1.-2.*nu_n)*(1.-2.*nu_n)*Z[0][2]-
    8.*(1.-5.*nu_n)*(1.-2.*nu_n)*Z[3][2];

  c=-pow(r[NumPhases-2],10)*(1.-2.*nu_n)*(7.+5*nu_n)*Z[0][1]+
    10.*pow(r[NumPhases-2],7)*(7.-nu_n*nu_n)*Z[3][1]+
    12.*pow(r[NumPhases-2],5)*(1.-2.*nu_n)*(Z[0][3]-7.*Z[1][2])+
    20.*pow(r[NumPhases-2],3)*(1.-2.*nu_n)*(1.-2.*nu_n)*Z[0][2]-
    8.*(7.-5.*nu_n)*(1.-2.*nu_n)*Z[3][2];



  //solve quadratic equation, report higher result
  if((sqr=b*b-4.*a*c)<0){
    printf("shear modulus does not yield real number, terminating\n");
    exit(1);
  }

  if((-b-pow(sqr, 0.5))>=0){
    printf("Two solutions for shear modulus were found, continuing with the higher value\n");
  }

  //usually this higher value is reported
  mu_hmg=mu[NumPhases-2]*(-b+pow(sqr, 0.5))/(2.*a);
  //uncomment if lower value of mu should be used, only when positive
  //mu_hmg=mu[NumPhases-2]*(-b+pow(sqr, 0.5))/(2.*a);
  
  //printf("Homogenized mu=%f\n", mu_hmg);
  
  KMuToENu(k_hmg, mu_hmg, E_hmg, nu_hmg);   

  destrv(r);
  destrv(k);
  destrv(mu);

  destrv(A);
  destrv(B);
  destrv(C);
  destrv(D);
  destrv(P_v);
  destrv(Phelp);

  destrm(L);
  destrm(Linv);
  destrm(M);
  destrm(Mhelp);
  destrm(Z);
  delete [] P_arr;
}

//help function for filling J matrix
void HomogData::fillJ(matrix &J, double r, const vector &mu, const vector &k, int phase){
  J[0][0]=r;
  J[0][1]=1./(pow(r,2));
  J[1][0]=3*k[phase];
  J[1][1]=-4.*mu[phase]/(pow(r,3));
}

//help function for filling L matrix
void HomogData::fillL(matrix &L, double r, const vector &mu, const vector &k, int phase){
  double mu_l=mu[phase];
  double k_l=k[phase];
  double nu_l=(3.*k_l-2.*mu_l)/(6.*k_l+2.*mu_l);

  L[0][0]=r;
  L[0][1]=-6.*nu_l*pow(r,3)/(1.-2.*nu_l);
  L[0][2]=3./pow(r,4);
  L[0][3]=(5.-4.*nu_l)/((1.-2.*nu_l)*pow(r,2));

  L[1][0]=r;
  L[1][1]=-(7.-4.*nu_l)*pow(r,3)/(1.-2.*nu_l);
  L[1][2]=-2./pow(r,4);
  L[1][3]=2./pow(r,2);
  
  L[2][0]=mu_l;
  L[2][1]=3.*nu_l*mu_l*pow(r,2)/(1.-2.*nu_l);
  L[2][2]=-12.*mu_l/pow(r,5);
  L[2][3]=2.*(nu_l-5.)*mu_l/((1.-2.*nu_l)*pow(r,3));

  L[3][0]=mu_l;
  L[3][1]=-(7.+2*nu_l)*mu_l*pow(r,2)/(1.-2.*nu_l);
  L[3][2]=8.*mu_l/pow(r,5);
  L[3][3]=2.*(1+nu_l)*mu_l/((1.-2.*nu_l)*pow(r,3));
}



/*
Kuster-Toksoz model from Monteiro:Elastic Moduli Lightweight Aggregate Model,
Cement and Concrete Research, 1995, No.2, pp. 276-280
assuming spheres as inclusions (aggregates), derived from
wave propagation velocity, 
model is not symmetric
inclusions are aggregates (mu_i, k_i)
*/
void HomogData::KusterToksoz(void){
  double k_KT, mu_KT, nom, denom;
  nom=1+(4*mu_m*(k_i-k_m)/((3*k_i*+4*mu_m)*k_m))*f_i;
  denom=1-(3*(k_i-k_m)/(3*k_i+4*mu_m))*f_i;
  k_KT=k_m*nom/denom;
  
  nom=(6*k_m+12*mu_m)*mu_i+(9*k_m+8*mu_m)*(f_m*mu_m+f_i*mu_i);
  denom=(9*k_m+8*mu_m)*mu_m+(6*k_m+12*mu_m)*(f_m*mu_i+f_i*mu_m);
  mu_KT=mu_m*nom/denom;

  KMuToENu( k_KT, mu_KT, E_hmg, nu_hmg );
  
  k_hmg=k_KT;
  mu_hmg=mu_KT;
}

void HomogData::ENuToKMu(const double oE, const double onu, double &ok, double &omu )
{
  omu = oE / ( 2. * ( 1. + onu ) );
  ok = oE / ( 3. * ( 1. - 2. * onu ) );
}

void HomogData::KMuToENu(const double ok, const double omu, double &oE, double &onu )
{
  oE  = 9. * ok * omu / ( 3. * ok + omu );
  onu = ( 3. * ok - 2. * omu ) / ( 6. * ok + 2 * omu );
}

void HomogData::Swap(double &A, double &B){
  double C;
  C=A;
  A=B;
  B=C;
}
