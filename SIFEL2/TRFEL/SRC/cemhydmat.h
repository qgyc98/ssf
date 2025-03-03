//responsible: smilauer@cml.fsv.cvut.cz

#ifndef CEMHYDMAT_H
#define CEMHYDMAT_H

/*CEMHYD3D has been developed at NIST, released under GNU public licence*/
/*programmed by D.P.Bentz*/
/*modified by smilauer@cml.fsv.cvut.cz*/

#include <stdio.h>

#define TRFEL_CEM

#ifdef TRFEL_CEM
#include "xfile.h"
#include"cmlfile.h"
struct matrix;
#else
#include"cmlfile.h"
#endif

typedef struct FCOMPLEX {float r,i;} fcomplex_cem;

class cemhydmat
{
 public:
  cemhydmat (void);
  cemhydmat (cemhydmat *CemP);
  ~cemhydmat (void);
  FILE *in;
  void read (char *inp);
  double GetHeat(double GiveTemp, double TargTime);
  double MoveCycles(double GiveTemp, int cycles);
  int MoveToDoH(double GiveTemp, double DesiredDoH, int maxcyc);
  int MoveToTime(double GiveTemp, double TargTime);
  double GetTotCemHeat(void);
  double GetTotHeat(void);
  double GetCp(void);
  double GetDoHLastCyc(void);
  double GetDoHActual(void);
  int GetCycNum(void);
  double GetCycTime(void);
  void CreateHDCSH(void);
  void PercolateForOutput(void);
  double GetWcr(void);
  void GetInputParams(char *my_string);
  void constructor_init(void);
  void AnalyticHomogenizationPaste(double &E, double &nu, int perc_unperc_flag);
  void AnalyticHomogenizationConcrete(double E_paste_inp, double nu_paste_inp, double *E_paste, double *nu_paste, double *E_mortar, double *nu_mortar, double &E_concrete, double &nu_concrete);
  void GetInitClinkerPhases(double &c3s, double &c2s, double &c3a, double &c4af, double &gypsum, double &hemi, double & anh);

#ifdef TRFEL_CEM
  //needed from trfel
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);
  void matcap (double &cc,long ri,long ci,long ipp);
  double get_k();
  double get_c();
  void read (XFILE *in);
  void print (FILE *out);
  void updateval(long ipp,long im,long ido);
  void initvalues(long ipp,long im,long ido);
  double get_source(long ipp,long im,long ido);

  //average temperature
  double average_temp;
  //number of int. points with cemhyd meterial
  long n_ipp;
  //inital material time for growing problems
  double init_material_time;
  //  coefficient of conductivity
  double conductivity;
  //  coefficient of capacity
  double capacity;
  char *files;
#endif


  int readcem (char *inp);
  int SYSSIZE;
  int SYSIZE;
  private:
  
  double TargTime_glob;
  /*define dimension size for image reconstruction and hydration*/
  /*Following parameters may be changed if you know what they are for*/
  int NEIGHBORS;  /* number of neighbors to consider (6, 18, or 26) in dissolution */
 
  int BoxSize;    /*int describing vicinity of CSH*/
  int SolidLimit; /*how many solid phase voxels must be in a box (max. <=(2*BoxSize+1)^3)*/
  long MAXTRIES; /* maximum number of random tries for sphere placement */
  long int MAXCYC_SEAL; /* Maximum number of cycles of sealed hydration (originally MAXCYC in disrealnew.c */
  
  /*Following parameters should not be changed*/
  long SYSIZE_POW3;

  /* Note that each particle must have a separate ID to allow for flocculation */
  int CEM;        /* and greater */
  int CEMID;	        /* phase identifier for cement */
  int C2SID;		/* phase identified for C2S cement */
  int GYPID;	        /* phase identifier for gypsum */	
  int HEMIHYDRATE;   /* phase identifier for hemihydrate */
  int POZZID;        /* phase identifier for pozzolanic material */
  int INERTID;       /* phase identifier for inert material */
  int SLAGID;	/* phase identifier for slag */
  int AGG;          /* phase identifier for flat aggregate */
  int FLYASH;       /* phase identifier for all fly ash components */

  long NPARTC;
  long BURNTG;     /* this value must be at least 100 > NPARTC */
  int NUMSIZES;     /* maximum number of different particle sizes */

  //+distrib3d
  long MAXSPH;   /* maximum number of elements in a spherical template */
  double PI;
  
  /*define heat capacities for all components in J/g/C*/
  /*including free and bound water*/
  double Cp_pozz; 
  double Cp_CH;
  double Cp_h2o;     /* Cp for free water */
  double Cp_bh2o;    /* Cp for bound water */
  double WN;         /* water bound per gram of cement during hydration */
  double WCHSH;      /* water imbibed per gram of cement during chemical shrinkage (estimate) */

  int CUBEMAX;
  int CUBEMIN;      /* Minimum cube size for checking pore size */
  long SYSIZEM1; /* System size -1 */

  double DISBIAS;  /* Dissolution bias- to change all dissolution rates */
  double DISMIN;  /* Minimum dissolution for C3S dissolution */
  double DISMIN2;  /* Minimum dissolution for C2S dissolution */
  double DISMINSLAG;  /* Minimum dissolution for SLAG dissolution */
  double DISMINASG;  /* Minimum dissolution for ASG dissolution */
  double DISMINCAS2;  /* Minimum dissolution for CAS2 dissolution */
  double DISMIN_C3A_0;  /* Minimum dissolution for C3A dissolution */
  double DISMIN_C4AF_0;  /* Minimum dissolution for C4AF dissolution */
  double DETTRMAX;
  double DGYPMAX;
  double DCACO3MAX;
  double DCACL2MAX;
  double DCAS2MAX;
  double CHCRIT;
  double C3AH6CRIT;
  double C3AH6GROW;  /* Probability for C3AH6 growth */
  double CHGROW;      /* Probability for CH growth */
  double CHGROWAGG;    /* Probability for CH growth on aggregate surface */
  double ETTRGROW;   /* Probability for ettringite growth */
  double C3AETTR;    /* Probability for reaction of diffusing C3A with ettringite  */
  double C3AGYP;      /*  Probability for diffusing C3A to react with diffusing gypsum */
  /* diffusing anhydrite, and diffusing hemihydrate */
  double SOLIDC3AGYP;    /* Probability of solid C3A to react with diffusing sulfate */
  double SOLIDC4AFGYP;  /* Probability of solid C4AF to react with diffusing sulfate */
  double PPOZZ;      /* base probability for pozzolanic reaction */
  double PCSH2CSH;  /* probability for CSH dissolution */
  /* for conversion of C-S-H to pozz. C-S-H */
  double A0_CHSOL;  /* Parameters for variation of CH solubility with */
  double A1_CHSOL;   /* temperature (data from Taylor- Cement Chemistry) */
  /* changed CSHSCALE to 70000 6/15/01  to better model induction CS */
  double CSHSCALE; /*scale factor for CSH controlling induction */
  double C3AH6_SCALE;/*scale factor for C3AH6 controlling induction of aluminates */

  int BURNT;      /* label for a burnt pixel <255 (in char type arrays) */
  long SIZE2D; /* size of matrices for holding burning locations */
  /* functions defining coordinates for burning in any of three directions */
  // definition of general equation of plane xy, yz, xz through cube
  long cx(int x,int y,int z,int a,int b,int c);
  long cy(int x,int y,int z,int a,int b,int c);
  long cz(int x,int y,int z,int a,int b,int c);

  long SIZESET;
  double AGRATE;        /* Probability of gypsum absorption by CSH */
  double VOLFACTOR; /* dm per pixel Note- dm*dm*dm = Liters */
  double MASSFACTOR;  /* cm per pixel - specific gravities in g/cm^3 */
  double MMNa;
  double MMK;
  double MMNa2O;
  double MMK2O;
  double BNa;     /* From Taylor paper in liters (31 mL/1000/ 100 g) */
  double BK;     /* From Taylor paper in liters (20 mL/1000/ 100 g) */
  double BprimeNa; /* From Taylor paper in liters (3 mL/1000/ 1 g POZZ) */
  double BprimeK; /* From Taylor paper in liters (3.3 mL/1000/ 1 g POZZ) */
  double KspCH25C;
  double KspGypsum;
  double KspSyngenite;
  double SpecgravSyngenite;   /* Source Taylor, H.F.W., Cement Chemistry */
  double KperSyn;   /* moles of K+ per mole of syngenite */

  double activeA0;     /* A at 295 K (from Ken) */
  double activeB0;   /* B at 295 K (from Ken) */
  /* z are the absolute charges (valences) per ion */
  double zCa;
  double zSO4;
  double zOH;
  double zNa;
  double zK;
  /* a is similar to an ionic radius (in Angstroms) */
  double aK;
  double aCa;
  double aOH;
  double aNa;
  double aSO4;     /* Estimate as S ionic radii + O ionic diameter */
  /* Ionic conductivities (From Snyder, Feng, Keen, and Mason) */
  /* and from CRC Hanbook of Chemistry and Physics (1983) pp. D-175 */
  /* pore solution conductivity = sum (zi * [i]*lambdai) */
  /* lambdai = (lambdai_0/(1.+Gi*(Istrength^0.5))) */
  /* where Istrength is in units of M (mol/L) */
  double lambdaOH_0;  /* Units: S cm-cm eq.^(-1) */
  double lambdaNa_0;
  double lambdaK_0;
  double lambdaSO4_0;
  double lambdaCa_0;     /* Note that CRC has 60./2 for this */
  double GOH;   /* Units: (eq.^2 mol/L)^(-0.5) */
  double GK;
  double GNa;
  double GCa;
  double GSO4;
  double cm2perL2m;  /* Conversion from cm2/Liter to 1/m */
  double EPSS;
  double MAXIT;
  double EPSP;
  int MAXM;
 
  int xoff[27];
  int yoff[27];
  int zoff[27];
 
  //random generator
  long IA;
  long IM;
  long IQ;
  int IR;
  int NTAB;
  double EPS;
  double NDIV; //= 1.0/(1.0+(IM-1.0)/NTAB);
  double RNMX;// = (1.0-EPS);
  double AM; //= (1.0/IM);
  int iy;
  int *iv;

  //phases
  int POROSITY;
  int C3S;
  int C2S;
  int C3A;
  int C4AF;
  int GYPSUM;
  int HEMIHYD;
  int ANHYDRITE;
  int POZZ;
  int INERT;
  int SLAG;
  int ASG;		  /* aluminosilicate glass */
  int CAS2;
  int CH;
  int CSH;
  int C3AH6;
  int ETTR;
  int ETTRC4AF;     /* Iron-rich stable ettringite phase */
  int AFM;
  int FH3;
  int POZZCSH;
  int SLAGCSH;      /* Slag gel-hydration product */
  int CACL2;
  int FREIDEL;	  /* Freidel's salt */
  int STRAT;	  /* stratlingite (C2ASH8) */
  int GYPSUMS;	/* Gypsum formed from hemihydrate and anhydrite */
  int CACO3;
  int AFMC;
  int INERTAGG;
  int ABSGYP;
  int DIFFCSH;
  int DIFFCH;
  int DIFFGYP;
  int DIFFC3A;
  int DIFFC4A;
  int DIFFFH3;
  int DIFFETTR;
  int DIFFCACO3;
  int DIFFAS;
  int DIFFANH;
  int DIFFHEM;
  int DIFFCAS2;
  int DIFFCACL2;
  int EMPTYP;         /*Empty porosity due to self desiccation*/
  int HDCSH;
  int OFFSET;           /*Offset for highlighted potentially soluble pixel*/
 
  //genpartnew
/* Note that each particle must have a separate ID to allow for flocculation */
//+distrib3d
  char ***mic; //char mic [SYSIZE] [SYSIZE] [SYSIZE];
 private:
  cmlfile *F;
  double ran1(int *idum);
  void addagg(void);
  int chksph(int xin,int yin,int zin,int radd,int wflg,int phasein,int phase2);
  int gsphere(int numgen,long int *numeach,int *sizeeach, int *pheach);
  int create(void);
  void drawfloc(int xin,int yin,int zin,int radd,int phasein,int phase2);
  int chkfloc(int xin,int yin,int zin,int radd);
  void makefloc(void);
  void measure(void);
  void measagg(void);
  void connect(void);
  void outmic(void);
  int genpartnew(void);
  void alloc_char_3D(char***(&mic), long SYSIZE);
  void dealloc_char_3D(char***(&mic), long SYSIZE);
  void alloc_long_3D(long***(&mic), long SYSIZE);
  void dealloc_long_3D(long***(&mic), long SYSIZE);
  void alloc_int_3D(int ***(&mask), long SYSIZE); 
  void dealloc_int_3D(int ***(&mask), long SYSIZE); 
  void alloc_shortint_3D(short int ***(&mic), long SYSIZE);
  void dealloc_shortint_3D(short int ***(&mic), long SYSIZE);
  void alloc_double_3D(double ***(&mic), long SYSIZE);
  void dealloc_double_3D(double ***(&mic), long SYSIZE);
 
  char ***micorig;//char micorig [SYSIZE] [SYSIZE] [SYSIZE];
  long int ***micpart;//long int micpart [SYSIZE] [SYSIZE] [SYSIZE];

  //genpartnew
  /* data structure for clusters to be used in flocculation */
  struct cluster{
    int partid;	/* index for particle */
    int clustid;  	/* ID for cluster to which this particle belongs */
    int partphase;  /* phase identifier for this particle */
    int x,y,z,r;	/* particle centroid and radius in pixels */
    struct cluster *nextpart;  /* pointer to next particle in cluster */
  };

  /* 3-D particle structure (each particle has own ID) stored in array cement */
  /* 3-D microstructure is stored in 3-D array cemreal */
  
  //define long int cement [SYSSIZE+1] [SYSSIZE+1] [SYSSIZE+1];
  long int ***cement; 
  //define long int cemreal [SYSSIZE+1] [SYSSIZE+1] [SYSSIZE+1];
  long int ***cemreal;

  int npart,aggsize;     /* global number of particles and size of aggregate */
  int iseed, nseed, *seed;      /* random number seed- global */
  int dispdist;   /* dispersion distance in pixels */
  int clusleft;   /* number of clusters in system */
  /* parameters to aid in obtaining correct sulfate content */
  long int n_sulfate,target_sulfate,n_total,target_total,volpart[47];
  long int n_anhydrite,target_anhydrite,n_hemi,target_hemi;
  double probgyp,probhem,probanh; /* probability of gypsum particle instead of cement */
  /* and probabilities of anhydrite and hemihydrate */
  /* relative to total sulfate */
  //  struct cluster *clust[NPARTC];	/* limit of NPARTC particles/clusters */
  struct cluster **clust;

  //distrib3d
  int maketemp(int size);   
  void phcount(void);
  int surfpix(int xin,int yin,int zin);
  float rhcalc(int phin);
  int countem(int xp,int yp,int zp,int phin);
  void sysinit(int ph1,int ph2);
  void sysscan(int ph1,int ph2);
  int procsol(int nsearch);
  int procair(int nsearch);
  int movepix(int ntomove,int ph1,int ph2);
  void sinter3d(int ph1id,int ph2id,float rhtarget);
  void stat3d(void);
  void rand3d(int phasein,int phaseout, float xpt);
  void distrib3d(void);

  //int mask[SYSIZE+1][SYSIZE+1][SYSIZE+1];
  int ***mask;
  //unsigned short int curvature [SYSSIZE+1] [SYSSIZE+1] [SYSSIZE+1];
  int ***curvature;
  long int volume[50],surface[50];
  int nsph;
  int *xsph,*ysph,*zsph;
  long int nsolid[1500],nair[1500];

  //disrealnew_30, burn3d, burnset, hydreal, burn_phases, nrutils, complex
 public:
  int ***mic_CSH;
  int ***ArrPerc;
  int ***ConnNumbers;
  double *PhaseFrac;
  //double E_CSH_hmg,nu_CSH_hmg;
  //double E_CSH_hmg;
  //double SH_hmg_1;
 private:
  void init(void);
  int chckedge(int xck,int yck,int zck);
  void passone(int low,int high,int cycid,int cshexflag);
  int loccsh(int xcur,int ycur,int zcur,int extent);
  int countbox(int boxsize,int qx,int qy,int qz);
  int countboxc(int boxsize,int qx,int qy,int qz);
  void makeinert(long int ndesire);
  void extslagcsh(int xpres,int ypres,int zpres);
  void dissolve(int cycle);
  void addrand(int randid,long int nneed);
  void measuresurf(void);
  void resaturate(void);
  void outputUnpercImageFile(void);
  void readhydrparam(void);
  void disrealnew_init(void);
  void disrealnew(double GiveTemp, double TargTime, int flag);
  int burn3d(int npix, int d1, int d2, int d3);
  int burnset(int d1, int d2, int d3);
  void parthyd(void);
  int moveone(int *xloc, int *yloc,int *zloc,int *act, int sumold);
  int edgecnt(int xck, int yck, int zck, int ph1, int ph2, int ph3);
  void extcsh(void);
  int movecsh(int xcur, int ycur, int zcur, int finalstep, int cycorig);
  void extfh3(int xpres, int ypres, int zpres);
  int extettr(int xpres, int ypres, int zpres, int etype);
  void extch(void);
  void extgyps(int xpres, int ypres, int zpres);
  int moveanh(int xcur, int ycur, int zcur, int finalstep, float nucprgyp);
  int movehem(int xcur, int ycur, int zcur, int finalstep, float nucprgyp);
  int extfreidel(int xpres, int ypres, int zpres);
  int extstrat(int xpres, int ypres, int zpres);
  int movegyp(int xcur, int ycur, int zcur, int finalstep);
  int movecacl2(int xcur, int ycur, int zcur, int finalstep);  
  int movecas2(int xcur, int ycur, int zcur, int finalstep);
  int moveas(int xcur, int ycur, int zcur, int finalstep);
  int movecaco3(int xcur,int ycur,int zcur,int finalstep);
  void extafm(int xpres, int ypres, int zpres);
  int moveettr(int xcur, int ycur, int zcur, int finalstep);
  void extpozz(int xpres, int ypres, int zpres);
  int movefh3(int xcur, int ycur, int zcur, int finalstep, float nucprob);
  int movech(int xcur, int ycur, int zcur, int finalstep, float nucprob);
  void extc3ah6(int xpres, int ypres, int zpres);
  int movec3a(int xcur, int ycur, int zcur, int finalstep, float nucprob);
  int movec4a(int xcur, int ycur, int zcur, int finalstep, float nucprob);
  void hydrate(int fincyc, int stepmax, float chpar1, float chpar2, float hgpar1, float hgpar2, float fhpar1, float fhpar2, float gypar1, float gypar2);
  void laguer(fcomplex_cem a[],int m,fcomplex_cem *x,float eps,int polish);
  void zroots(fcomplex_cem a[],int m,fcomplex_cem roots[],int polish);
  void pHpred(void);

  int IsSolidPhase(int phase);
  void burn_phases(int d1, int d2, int d3);
  int IsConnected(int cx, int cy, int cz, int dx, int dy, int dz);
  void GenerateConnNumbers(void);
  void outputImageFilePerc(void);
  void WriteUnsortedList(int px, int py, int pz);
  void CountPercolation(int &tot_perc, int &tot_unperc);
  inline int AdjCoord(int coord);
  int NumSol(int cx, int cy, int cz);
  void CSHbox(unsigned int *CSH_vicinity);
  void nrerror(const char *error_text);
  float* vector(int nl,int nh);
  int* ivector(int nl,int nh);
  double* dvector(int nl,int nh);
  //float** matrix(int nrl,int nrh,int ncl,int nch);
  float** matrix_cem(int nrl,int nrh,int ncl,int nch);
  double** dmatrix(int nrl,int nrh,int ncl,int nch);
  int** imatrix(int nrl,int nrh,int ncl,int nch);
  float** submatrix(float **a,int oldrl,int oldrh,int oldcl,int oldch,int newrl,int newcl);
  void free_vector(float *v, int nl);
  void free_ivector(int *v,int nl);
  void free_dvector(double *v,int nl);
  void free_matrix(float **m,int nrl,int nrh, int ncl);
  void free_dmatrix(double **m,int nrl,int nrh, int ncl);
  void free_imatrix(int **m,int nrl,int nrh,int ncl);
  void free_submatrix(float *b,int nrl);
  float** convert_matrix(float *a,int nrl,int nrh,int ncl,int nch);
  void free_convert_matrix(float **b,int nrl);
  double *last_values;
  long int *phase;

  fcomplex_cem Cadd(fcomplex_cem a,fcomplex_cem b);
  fcomplex_cem Csub(fcomplex_cem a,fcomplex_cem b);
  fcomplex_cem Cmul(fcomplex_cem a,fcomplex_cem b);
  fcomplex_cem Complex(float re,float im);
  fcomplex_cem Conjg(fcomplex_cem z);
  fcomplex_cem Cdiv(fcomplex_cem a,fcomplex_cem b);
  float Cabs(fcomplex_cem z);
  fcomplex_cem Csqrt(fcomplex_cem z);
  fcomplex_cem RCmul(float x,fcomplex_cem a);

  /* data structure for diffusing species - to be dynamically allocated */
  /* Use of a doubly linked list to allow for easy maintenance */
  /* (i.e. insertion and deletion) */
  /* Added 11/94 */
  /* Note that if SYSIZE exceeds 256, need to change x, y, and z to */
  /* int variables */
  struct ants{
    unsigned char x,y,z,id;
    int cycbirth;
    struct ants *nextant;
    struct ants *prevant;
  };
  
  /* data structure for elements to remove to simulate self-desiccation */
  /* once again a doubly linked list */
  struct togo{
    int x,y,z,npore;
    struct togo *nexttogo;
    struct togo *prevtogo;
  };
  //short int cshage [SYSIZE] [SYSIZE] [SYSIZE];
  short int ***cshage;
  //short int faces [SYSIZE] [SYSIZE] [SYSIZE];
  short int ***faces;
  unsigned int *CSH_vicinity;//[(2*BoxSize+1)*(2*BoxSize+1)*(2*BoxSize+1)+1];
  /* counts for dissolved and solid species */
  long int *discount, *count;
  long int ncshplategrow,ncshplateinit;
  /* Counts for pozzolan reacted, initial pozzolan, gypsum, ettringite,
     initial porosity, and aluminosilicate reacted */
  long int npr,nfill,ncsbar,netbar,porinit,nasr,nslagr,slagemptyp;
  /* Initial clinker phase counts */
  long int c3sinit,c2sinit,c3ainit,c4afinit,anhinit,heminit,chold,chnew;
  long int nmade,ngoing,gypready,poregone,poretodo,countpore;
  long int countkeep,water_left,water_off,pore_off;
  int ncyc,cyccnt,cubesize,sealed,outfreq,ImgOut;
  int icyc,burnfreq,setfreq,setflag,sf1,sf2,sf3,porefl1,porefl2,porefl3;
  /* Parameters for kinetic modelling ---- maturity approach */
  double ind_time,temp_0,temp_cur,time_step,time_cur,E_act,beta;
  /*define heat conversion factor for the cement mixture, includes all solids as
    count[INERT]+count[SLAG]+count[POZZ]+count[CACL2]+count[ASG]+count[CAS2]*/
  double heat_cf;
  float w_to_c,s_to_c,krate,totfract,tfractw04,fractwithfill;
  float tfractw05,surffract,pfract,pfractw05,sulf_conc;
  long int scntcement,scnttotal;
  float U_coeff,T_ambient;
  double alpha_cur,alpha_last,heat_old,heat_new,cemmass,mass_agg,mass_water,mass_fill,Cp_now, Cp_agg,Cp_cement, cem_m3;
  double alpha,CH_mass,mass_CH,mass_fill_pozz,E_act_pozz,chs_new,cemmasswgyp;
  float flyashmass,alpha_fa_cur;
  double E_act_slag;
  double TargDoHelas;
  /* Arrays for variable CSH molar volume and water consumption */
  float *molarvcsh,*watercsh;
  float heatsum,molesh2o,saturation;
  /* Arrays for dissolution probabilities for each phase */
  float *disprob,*disbase;
  float gypabsprob,ppozz;
  /* Arrays for specific gravities, molar volumes, heats of formation, and */
  /* molar water consumption for each phase */
  float *specgrav,*molarv,*heatf,*waterc;
  /* Solubility flags and diffusing species created for each phase */
  /* Also flag for C1.7SH4.0 to C1.1SH3.9 conversion */
  int *soluble,*creates;
  int csh2flag,adiaflag,chflag,nummovsl;
  float cs_acc;     /* increases disprob[C3S] and disprob[C2S] if gypsum is present */
  float ca_acc;    /* increases disprob[C3A] and disprob[C4AF] if gypsum is present */
  float dismin_c3a;
  float dismin_c4af;
  float gsratio2,onepixelbias;
  /* Slag probabilities */
  float p1slag;   /*  probability SLAG is converted to SLAGCSH */
  float p2slag;   /*  probability SLAG is converted to POROSITY or EMPTYP */
  float p3slag;   /*  probability adjoining pixel is converted to SLAGCSH */
  float p4slag;   /*  probability CH is consumed during SLAG reaction */
  float p5slag;   /*  probability a C3A diffusing species is created */
  double slagcasi,slaghydcasi;   /* Ca/Si ratios for SLAG and SLAGCSH */
  float slagh2osi;     /* H/S ratio of SLAGCSH */
  double slagc3a;       /* C3A/slag molar ratio */
  double siperslag;     /* S ratio of SLAG (per mole) */
  double slagreact;     /* Base dissolution reactivity factor for SLAG */
  long int DIFFCHdeficit,slaginit;   /* Deficit in CH due to SLAG reaction */
  long int slagcum,chgone;
  long int nch_slag;    /* number of CH consumed by SLAG reaction */
  long int sulf_cur;
  long int sulf_solid;
  char heatname[80],adianame[80],phasname[80],ppsname[80],ptsaname[80],phrname[80];
  char chshrname[80],micname[80];
  char cmdnew[120],pHname[80],fileroot[80];
  struct ants *headant,*tailant;
  FILE *heatfile,*chsfile,*ptmpfile,*movfile,*pHfile,*micfile, *fileperc,*percfile, *disprobfile, *phasfile, *perc_phases, *CSHfile;
  //fileperc used in burn3d.cpp, percfile in burnset.cpp
  /* Variables for alkali predictions */
  double pH_cur,totsodium,totpotassium,rssodium,rspotassium;
  /* Array for whether pH influences phase solubility  -- added 2/12/02 */
  float *pHeffect;
  float pHfactor;
  int pHactive,resatcyc,cshgeom;
  /* Make conccaplus global to speed up execution and moles_syn_precip */
  /* global to accumulate properly */
  double conccaplus,moles_syn_precip,concsulfate; 
  int primevalues[6];
  int cshboxsize;		/* Box size for addition of extra diffusing C-S-H */
  //int newmat [SYSIZE][SYSIZE][SYSIZE];  moved to burn_phases and burnset
  
  struct percolatedpath{
    int x, y, z;
    struct percolatedpath *next;
    struct percolatedpath *prev;
  };
  
  struct percolatedpath *last, *current;

  //disrealnew
  long int adiabatic_curing;
  int ntimes;
  int cycflag;
  int phydfreq;
  double InitTime;
  double pnucch,pscalech,pnuchg,pscalehg,pnucfh3,pscalefh3;
  double pnucgyp,pscalegyp;
  float thtimelo,thtimehi,thtemplo,thtemphi;
  double mass_cement,mass_cem_now,mass_cur,kpozz,kslag;
  FILE *adiafile,*thfile,*elasfile;
  double LastHydrTime, LastCallTime,PrevHydrTime;
  double LastCycHeat, LastTotHeat, PrevCycHeat;
  long int LastCycCnt, Calculate_elastic_homogenization;
	double Vol_cement_clinker_gypsum, Vol_cement_SCM, Vol_water, Vol_FA, Vol_CA, Vol_inert_filler, Vol_entrained_entrapped_air, Grain_average_FA, Grain_average_CA, ITZ_thickness, ITZ_Young_red, Young_SCM, Poisson_SCM, Young_FA, Poisson_FA, Young_CA, Poisson_CA, Young_inert, Poisson_inert;

};

#endif
