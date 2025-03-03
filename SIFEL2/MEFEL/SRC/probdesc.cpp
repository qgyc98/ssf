#include <string.h>
#include "galias.h"
#include "probdesc.h"
#include "global.h"
#include "gmatrix.h"
#include "intools.h"
#include "iotools.h"
#include "aggregator.h"
#include "adaptivity.h"



/**
  Constructor initializes data members to zero or default values.

  Created by JK
  Modified by TKo, TKr
*/
probdesc::probdesc ()
{
  kwdsw = nokwd;
  niilnr=0;  errnr=0.0;

  //  type of stiffness matrix storage
  tstorsm = (storagetype) 0;
  //  type of mass matrix storage
  tstormm = (storagetype) 0;
  
  // diagonalization of mass matrix
  //  mass matrix is not diagonalized
  diagmass=0;
  
  
  tepsol = (epsolvertype) 0;

  //  initialization of damping description
  damp = nodamp;
  dmass=0.0;  dstiff=0.0;
  
  //  compute stresses during evaluation of internal forces
  strcomp = 1;

  //  compute eigenstresses during evaluation of right han side
  eigstrcomp = 1;

  //  for open dx
  detnodstrain = 0;
  
  stochasticcalc=0;
  eigstrains=0;
  temperature=0;
  pore_press=0;
  
  adaptivityflag=0;
  
  limit=0.0;  zero=1.0e-10;

  //neigv=0;  nev=0;  nies=0;  anies=0;  erres=0.0;  aerres=0.0;
  //nijmr=0;  njacthr=0;  jacthr=NULL;

  alphafvn=0.0;  deltafvn=0.0;

  time = 0.0;
  dtime = 0.0;
  dtmin=zero;

  istep = 0;
  jstep = 0;
  //phase = 0;
  nonlocphase = 0;
  //  internal forces are computed at nodes
  nodeintfor = 1;
  lambda = 0.0;
  dlambda = 0.0;
  
  niep =0;  stepep = 0.0; errep = 0.0;  rlerrep = 0.0;  nselnodep = 0;  selnodep=NULL; loadcoef = NULL;

  matmodel = (materialmodel) 0;
  
  //  strains are not computed
  straincomp=0;
  //  strains are not averaged
  strainaver=0;
  //  strain position
  strainpos=0;
  //  stresses are not computed
  stresscomp=0;
  //  stresses are not averaged
  stressaver=0;
  //  stress position
  stresspos=0;
  //  components of array other are not computed
  othercomp = 0;
  //  components of array other are not averaged
  otheraver = 0;
  //  position of component of array other
  otherpos=0;
  //  reactions are not computed
  reactcomp=0;
  
  strainstate=0;
  stressstate=0;
  otherstate=0;   
  
  memset(name, 0, sizeof(*name)*1001);

  tlm=0;
  
  //  detection of hanging nodes
  hangnodesdetect=no;
  
  //  indicator of homogenization
  homog=0;
  //  indicator stress-strain state for homogenization
  hstrastre = (strastrestate)0;
  //  direction of unit load for homogenization for homog=7 and homog=8
  homogdir = 0;
  //  type of macro-micro problem correspondence in homogenization
  mami = (macromicrotype) 0;


  //  the number of Wang's tiles
  ntiletypes=0;

  mcnne=0;  dim=0;
  
  ense=0;
  
  nivsm=0;
  path = filename = suffix = NULL;
  selnodep = NULL;
  loadcoef = NULL;
  nlman = NULL;

  //  type of time printing
  //tpr = days;
  tpr = seconds;
  //tpr = minutes;
  //tpr = hours;

  // initial displacements of elements are taken only from interface nodes
  // in case of growing problems
  comp_inidispl = no;
  // prescribed rotation of selected nodes at the initial displacement computation
  // in case of growing problem
  rot_inidispl = no;

  // smooth element removal procedure
  cpsmooth = no;     // cpsmooth=yes procedure will be applied on tractions in growing mechanical problem 
  cpincrnr = 0.0;    // initial load coeficient increment in the procedure
  cpminincrnr = 0.0; // minimum load coeficient increment in the procedure

  //  type of load balancing algorithm
  //  this choice means that the code does not stop after negm is known
  lbtype=0;
  
  //
  ssle = new slesolv();

  // backup on harddisk
  hdbcont.hdbtype = nohdb;  // no backup
  hdbcont.hdbfmtr = text;
  hdbcont.hdbfmts = text;
  hdbcont.hdbnames[0] = 0;
  hdbcont.hdbnamer[0] = 0;
  hdbcont.prec=15;
  hdbcont.selother_r  = NULL;
  hdbcont.selother_s  = NULL;
  hdbcont.selother_id = NULL;
}



/**
  Destructor releases allocated memory of the probdesc object.

  Created by JK
  Modified by TKo, TKr
*/
probdesc::~probdesc ()
{
  delete [] selnodep;
  delete [] loadcoef;
  delete [] nlman;

  delete [] path; 
  delete [] filename; 
  delete [] suffix;
  delete ssle;
}



/**
  The function reads basic data about solved problem.
   
  @param in - input stream
   
  @return The function does not return anything.

  Created by JK, 20.7.2001
  Modified by TKo, TKr
*/
void probdesc::read (XFILE *in)
{
  long i;
  
  //  problem name
  //inputln (in,name,1000);
  //  cteme1
  xfscanf (in,"% a",name);
  
  //  detail messages
  Mespr=0;
  //  cteme2
  xfscanf (in,"%k%ld","mespr",&Mespr);
  if (Mespr==1)  fprintf (stdout,"\n detail information will be printed");
  else fprintf (stdout,"\n only important messages will be printed");

  //  problem type
  xfscanf (in,"%k%m","problemtype",&problemtype_kwdset,(int*)&tprob);
  

  // *************************
  //  computation of strains
  // *************************
  xfscanf (in,"%k%ld","straincomp", &straincomp);
  if (straincomp!=0 && straincomp!=1)
    print_err("wrong definition of strain computation", __FILE__, __LINE__, __func__);
  if (straincomp==1){
    xfscanf (in,"%k%ld%k%ld","strainpos", &strainpos, "strainaver", &strainaver);
    
    if (strainpos!=1 && strainpos!=2 && strainpos!=3){
      print_err("wrong definition of position where strains are computed", __FILE__, __LINE__, __func__);
    }
    if (strainaver!=0 && strainaver!=1 && strainaver!=2)
      print_err("wrong definition of strain averaging", __FILE__, __LINE__, __func__);
  }
  
  if (Mespr==1){
    if (straincomp==0)  fprintf (stdout,"\n strains will not be computed and stored");
    else{
      fprintf (stdout,"\n strains will be computed and stored");
      if (strainpos==1)   fprintf (stdout,"\n strains will be computed at integration points");
      if (strainpos==2)   fprintf (stdout,"\n strains will be copied to nodes from the closest integration point");
      if (strainpos==3)   fprintf (stdout,"\n strains will be computed at nodes");
      
      if (strainaver==0)  fprintf (stdout,"\n strains will not be averaged");
      if (strainaver==1)  fprintf (stdout,"\n strains will be averaged at nodes");
    }
  }
  
  
  // **************************
  //  computation of stresses
  // **************************
  xfscanf (in,"%k%ld","stresscomp", &stresscomp);
  if (stresscomp!=0 && stresscomp!=1)
    print_err("wrong definition of stress computation", __FILE__, __LINE__, __func__);
  if (stresscomp==1){
    xfscanf (in,"%k%ld%k%ld","stresspos", &stresspos, "stressaver", &stressaver);

    if (stresspos!=1 && stresspos!=2){
      print_err("wrong definition of position where stresses are computed", __FILE__, __LINE__, __func__);
    }
    if (stressaver!=0 && stressaver!=1 && stressaver!=2)
      print_err("wrong definition of stress averaging", __FILE__, __LINE__, __func__);
  }
  
  if (Mespr==1){
    if (stresscomp==0)  fprintf (stdout,"\n stresses will not be computed and stored");
    else{
      fprintf (stdout,"\n stresses will be computed and stored");
      if (stresspos==1)   fprintf (stdout,"\n stresses will be computed at integration points");
      if (stresspos==2)   fprintf (stdout,"\n stresses will be copied to nodes from the closest integration point");
      
      if (stressaver==0)  fprintf (stdout,"\n stresses will not be averaged");
      if (stressaver==1)  fprintf (stdout,"\n stresses will be averaged at nodes");
    }
  }
  
  // ****************************
  //  computation of other values
  // ****************************
  xfscanf (in,"%k%ld","othercomp", &othercomp);
  if (othercomp!=0 && othercomp!=1)
    print_err("wrong definition of other values computation", __FILE__, __LINE__, __func__);
  if (othercomp==1){
    xfscanf (in,"%k%ld%k%ld","otherpos", &otherpos, "otheraver", &otheraver);

    if (otherpos!=1 && otherpos!=2){
      print_err("wrong definition of position where other are computed", __FILE__, __LINE__, __func__);
    }
    if (otheraver!=0 && otheraver!=1 && otheraver!=2)
      print_err("wrong definition of other averaging", __FILE__, __LINE__, __func__);
  }
  
  if (Mespr==1){
    if (othercomp==1){
      fprintf (stdout,"\n other will be computed and stored");
      if (otherpos==1)   fprintf (stdout,"\n other will be computed at integration points");
      if (otherpos==2)   fprintf (stdout,"\n other will be copied to nodes from the closest integration point");
      
      if (otheraver==0)  fprintf (stdout,"\n other will not be averaged");
      if (otheraver==1)  fprintf (stdout,"\n other will be averaged at nodes");
    }
    else
      fprintf (stdout,"\n other will not be computed and stored");
  }
  
  // ***************************
  //  computation of reactions
  // ***************************
  xfscanf (in,"%k%ld","reactcomp", &reactcomp);
  if (reactcomp!=0 && reactcomp!=1)
    print_err("wrong definition of reaction computation", __FILE__, __LINE__, __func__);
  if (Mespr==1){
    if (reactcomp==0)  fprintf (stdout,"\n reactions will not be computed and stored");
    else  fprintf (stdout,"\n reactions will be computed and stored");
  }
  

  // **************************
  //   adaptivity computation
  // **************************
  xfscanf (in,"%k%ld","adaptivity",&adaptivityflag);
  if (adaptivityflag) {
    Ada = new adaptivity();
    Ada->readinit(in);
  }
  
  // **********************************
  //  stochastic or fuzzy computation
  // **********************************
  //  stochasticcalc=0 - deterministic computation
  //  stochasticcalc=1 - stochastic/fuzzy computation, data are read all at once
  //  stochasticcalc=2 - stochastic/fuzzy computation, data are read sequentially
  //  stochasticcalc=3 - stochastic/fuzzy computation, data are generated in the code
  xfscanf (in,"%k%ld","stochasticcalc",&stochasticcalc);
  if (Mespr==1){
    if (stochasticcalc==0)
      fprintf (stdout,"\n deterministic computation will be computed");
    if (stochasticcalc>0)
      fprintf (stdout,"\n non-deterministic (stochastic/fuzzy) computation will be computed");
    
  }
  
  // ****************************
  //  homogenization procedures
  // ****************************
  xfscanf (in,"%k%ld","homogenization",&homog);
  if (Mespr==1){
    if (homog==0)
      fprintf (stdout,"\n homogenization will not be computed");
    if (homog!=0)//changed 1.2.2006, smilauer@cml.fsv.cvut.cz
      fprintf (stdout,"\n homogenization will be computed");
    if (homog==2){
      //  Wang's tiles are used
      xfscanf (in,"%k%ld","numtiles",&ntiletypes);

      //  array ltg will be read in mechtop::read
      //Gtm->rst=1;
      
      //  code numbers will be generated with respect to the Schur complement method
      //  internal DOFs are ordered first, interface DOFs are ordered afterwards
      //  mozna to nebude treba
      //Gtm->cngen=2;
    }

    if (homog==3){
      fprintf (stdout,"\n general homogenization with prescribed macro-stresses");
    }
    if (homog==4){
      fprintf (stdout,"\n general homogenization with prescribed macro-strains");
    }

    if (homog==5 || homog==7){
      fprintf (stdout,"\n compliance matrix computed from homogenization with prescribed macro-stresses");
      xfscanf (in,"%m", &strastrestate_kwdset, (int*)&hstrastre);
    }
    if (homog==6 || homog==8){
      fprintf (stdout,"\n stiffness matrix computed from homogenization with prescribed macro-strains");
      xfscanf (in,"%m", &strastrestate_kwdset, (int*)&hstrastre);
    }
    if (homog==7 || homog==8){
      xfscanf (in,"%k%ld","homogdir",&homogdir);
      fprintf (stdout,"\n unifom load is prescribed in direction %ld",homogdir);
    }
    if (homog==9){
      fprintf (stdout,"\n general homogenization with mixed prescribed macro-strains/stresses");
    }
    if (homog==10){
      //  type of macro-micro problem correspondence
      //  mami=1 - elements are connected with microproblems
      //  mami=2 - aggregates of elements are connected with microproblems
      xfscanf (in,"%k%ld","macmiccorres",&mami);
      Gtm->rst=1;
    }
  }
  
  // ***************************
  //  renumbering of the nodes
  // ***************************
  xfscanf (in,"%k%m","noderenumber",&noderenumb_kwdset,(int*)&Gtm->nodren);
  if (Mespr==1){
    switch (Gtm->nodren){
    case no_renumbering:{
      fprintf (stdout,"\n nodes will not be renumbered");
      break;
    }
    case cuthill_mckee:{
      fprintf (stdout,"\n nodes will be renumbered by Cuthill-McKee algorithm");
      break;
    }
    case rev_cuthill_mckee:{
      fprintf (stdout,"\n nodes will be renumbered by reverse Cuthill-McKee algorithm");
      break;
    }
    case sloan:{
      fprintf (stdout,"\n nodes will be renumbered by Sloan algorithm");
      break;
    }
    default:{
      print_err("unknown type of node reordering is required",__FILE__,__LINE__,__func__);
    }
    }
  }

  // *****************************
  //  detection of hanging nodes
  // *****************************
  /*
  xfscanf (in,"%k%m","hangnodesdetect",&answertype_kwdset,(int*)&hangnodesdetect);
  if (hangnodesdetect==1){
    barlist.read (in);
  }
  */
  
  switch (tprob){
    // *****************
    //  LINEAR STATICS
    // *****************
  case linear_statics:{
    if (Mespr==1)  fprintf (stdout,"\n linear static problem");
    
    //  tstorsm - type of storage of stiffness matrix
    xfscanf (in,"%k%m","stiffmatstor",&storagetype_kwdset,(int*)&tstorsm);
    
    //  data about solver of system of linear equations
    ssle->read (Gtm,in,Mespr);
    
    break;
  }
    // *********************************************************************************
    //  EIGENDYNAMICS - EIGENVALUES AND EIGENVECTORS (EIGENFREQUENCIES AND EIGENMODES)
    // *********************************************************************************
  case eigen_dynamics:{
    if (Mespr==1)  fprintf (stdout,"\n eigenvalue vibration problem");
    
    //  type of storage of stiffness matrix
    //  type of storage of mass matrix
    xfscanf (in,"%k%m%k%m","stiffmatstor",&storagetype_kwdset,(int*)&tstorsm,"massmatstor",&storagetype_kwdset,(int*)&tstormm);
    
    //  data about solver of eigenvalues and eigenvectors
    eigsol.read (in);
    
    //  data about solver of system of linear equations
    ssle->read (Gtm,in,Mespr);

    // diagonalization of mass matrix
    xfscanf (in,"%k%ld","massdiagonalization",&diagmass);

    break;
  }
    // *************************************************************************************
    //  LINEAR STABILITY - EIGENVALUES AND EIGENVECTORS (CRITICAL LOAD AND BUCKLING MODES)
    // *************************************************************************************
  case linear_stability:{
    if (Mespr==1)  fprintf (stdout,"\n linear stability problem");
    
    //  type of storage of stiffness matrix
    //  type of storage of initial stiffness matrix
    xfscanf (in,"%k%m%k%m","stiffmatstor",&storagetype_kwdset,(int*)&tstorsm,"initstressmatstor",&storagetype_kwdset,(int*)&tstorgm);
    
    //  data about solver of eigenvalues and eigenvectors
    eigsol.read (in);
    
    //  data about solver of system of linear equations
    ssle->read (Gtm,in,Mespr);

    break;
  }
    // ******************
    //  FORCED DYNAMICS
    // ******************
  case forced_dynamics:{
    if (Mespr==1)  fprintf (stdout,"\n forced vibration problem");
    
    //  type of solver of forced vibration problem
    xfscanf (in,"%k%m","forceddyn",&forcedsolver_kwdset,(int*)&tforvib);
    
    switch (tforvib){
    case newmark:{
      //  optimal values
      //  alphafvn = 0.25
      //  deltafvn = 0.5
      xfscanf (in,"%lf %lf",&alphafvn,&deltafvn);
      break;
    }
    case findiff:{
      break;
    }
    case explfindiff:{
      break;
    }
    case modal_analysis:{
      //  data about eigenmode/eigenvectors
      eigsol.read (in);
      break;
    }
    default:{
      print_err("unknown solver of forced vibration is required",__FILE__,__LINE__,__func__);
    }
    }
    
    //  type of damping
    xfscanf (in,"%k%m","damping",&damping_kwdset,(int*)&damp);
    switch (damp){
    case nodamp:{
      if (Mespr==1)  fprintf (stdout,"\n no damping is prescribed");
      break;
    }
    case massdamp:{
      xfscanf (in,"%lf",&dmass);
      if (Mespr==1)  fprintf (stdout,"\n damping proportional to the mass matrix is prescribed");
      break;
    }
    case stiffdamp:{
      xfscanf (in,"%lf",&dstiff);
      if (Mespr==1)  fprintf (stdout,"\n damping proportional to the stiffness matrix is prescribed");
      break;
    }
    case rayleighdamp:{
      xfscanf (in,"%lf %lf",&dstiff,&dmass);
      if (Mespr==1)  fprintf (stdout,"\n Rayleigh damping (proportional to the mass and stiffness matrices) is prescribed");
      break;
    }
    case elemdamp:{
      if (Mespr==1)  fprintf (stdout,"\n damping is defined on finite elements");
      break;
    }
    default:{
      print_err("unknown type of damping is required",__FILE__,__LINE__,__func__);
    }
    }

    //  type of storage of stiffness matrix
    //  type of storage of mass matrix
    xfscanf (in,"%k%m%k%m","stiffmatstor",&storagetype_kwdset,(int*)&tstorsm,"massmatstor",&storagetype_kwdset,(int*)&tstormm);
    //  type of storage of damping matrix is defined by type of damping
    if (damp==elemdamp){
      xfscanf (in,"%k%m","dampmatstor",&storagetype_kwdset,(int*)&tstordm);
    }

    // diagonalization of mass matrix
    //  diagmass=1 - mass matrix is diagonalized
    xfscanf (in,"%k%ld","massdiagonalization",&diagmass);
    
    //  time controller
    timecon.read (in);

    //  data about solver of system of linear equations
    ssle->read (Gtm,in,Mespr);

    break;
  }
    // *******************************
    //  MATERIALLY NONLINEAR STATICS
    // *******************************
  case mat_nonlinear_statics:{
    if (Mespr==1)  fprintf (stdout,"\n nonlinear static problem");
    
    //  data about nonlinear solver
    nlman = new nonlinman[1];
    nlman->read (in,Mespr);
    
    // backup controler
    hdbcont.read(in);
    if ((Mespr==1) && hdbcont.hdbtype)
      fprintf(stdout, "\n backup will be performed %d", hdbcont.hdbtype);

    //  type of storage of stiffness matrix
    xfscanf (in,"%k%m","stiffmatstor", &storagetype_kwdset, (int*)&tstorsm);
    
    //  data about solver of system of linear equations
    ssle->read (Gtm,in,Mespr);
    
    break;
  }
    // *********************************************************************
    //  TIME-DEPENDENT MECHANICAL PROBLEMS WITH NEGLIGIBLE INERTIAL FORCES
    // *********************************************************************
  case mech_timedependent_prob:{
    if (Mespr==1)  fprintf (stdout,"\n system of nonlinear equations will be solved by visco-solver");
    
    //  time controller
    timecon.read (in);

    //  type of time printing
    xfscanf (in,"%k%m","timetypeprin", &timetypeprin_kwdset, &tpr);

    // backup controler
    hdbcont.read(in);
    if ((Mespr==1) && hdbcont.hdbtype)
      fprintf(stdout, "\n backup will be performed %d", hdbcont.hdbtype);

    nlman = new nonlinman[1];
    //  maximum number of iterations in inner loop
    //  required norm of residuals
    xfscanf (in,"%k%ld%k%lf","nr_num_iter",&nlman->niilnr,"nr_error",&nlman->errnr);
    xfscanf (in, "%k%m", "resid_norm_type", &resnormt_kwdset, (int*)&nlman->rnormtnr);
    
    // checking of divergency of inner iteration loop
    xfscanf(in, "%k%m", "check_div", &flagsw_kwdset, (int *)&nlman->check_div);
    if (Mp->nlman->check_div == on)
    {
      xfscanf(in, "%k%ld", "div_min_steps", &nlman->div_min_steps);
      if (Mp->nlman->div_min_steps < 3)
      {
        print_warning("the minimum number of performed steps in the divergency check is < 3.\n"
                      "It will be set to 3 automatically", __FILE__, __LINE__, __func__);
       Mp->nlman->div_min_steps = 3;
      }
    }

    //  type of storage of stiffness matrix
    //  type of stiffness matrix (initial, tangent)
    xfscanf (in,"%k%m%k%m","stiffmatstor", &storagetype_kwdset, (int*)&tstorsm, "stiffmat_type", &stiffmatrix_kwdset, (int*)&nlman->stmat);
    if(nlman->stmat == ijth_tangent_stiff)
      xfscanf (in,"%ld %ld", &nlman->ithstiffchange, &nlman->jthstiffchange);

    //  data about solver of system of linear equations
    ssle->read (Gtm,in,Mespr);

    break;
  }
    // **********************************************************************************************
    //  TIME-DEPENDENT MECHANICAL PROBLEMS WITH NEGLIGIBLE INERTIAL FORCES AND GRADUAL CONSTRUCTION
    // **********************************************************************************************
  case growing_mech_structure:{
    if (Mespr==1)  fprintf (stdout,"\n system with gradual construction will be solved by visco-solver");
    
    //  time controller
    timecon.read (in);

    //  type of time printing
    xfscanf (in,"%k%m","timetypeprin", &timetypeprin_kwdset, &tpr);

    // calculation of initial displacements due to interface nodes
    xfscanf (in, "%k%m", "comp_inidispl", &answertype_kwdset, &comp_inidispl);

    // prescribed initial displacements by rotation of selected nodes
    if (comp_inidispl == yes)
      xfscanf (in, "%k%m", "rot_inidispl", &answertype_kwdset, &rot_inidispl);

    // smooth element removal procedure (gradual application of traction due to element removal at interface nodes)
    xfscanf (in, "%k%m", "smooth_removal", &answertype_kwdset, &cpsmooth);
    if (cpsmooth == yes)
      // initial and minimum load coefficient increment in smoothed element removal procedure
      xfscanf (in, "%k%le %k%le", "cpincrnr", &cpincrnr, "cpminincrnr", &cpminincrnr);

    // backup controler
    hdbcont.read(in);

    if ((Mespr==1) && hdbcont.hdbtype)
      fprintf(stdout, "\n backup will be performed");

    nlman = new nonlinman[1];

    //  maximum number of iterations in inner loop
    //  required norm of residuals
    xfscanf (in,"%k%ld%k%lf","nr_num_iter",&nlman->niilnr,"nr_error",&nlman->errnr);
    
    //  type of nonlinear solver
    xfscanf (in,"%k%m","nr_solvtype", &nonlintimesolvertype_kwdset, (int*)&nrsolv);

    //  type of storage of stiffness matrix
    //  type of stiffness matrix (initial, tangent)
    xfscanf (in,"%k%m%k%m","stiffmatstor", &storagetype_kwdset, (int*)&tstorsm, "stiffmat_type", &stiffmatrix_kwdset, (int*)&nlman->stmat);
    if(nlman->stmat == ijth_tangent_stiff)
      xfscanf (in,"%ld %ld", &nlman->ithstiffchange, &nlman->jthstiffchange);
    
    //  data about solver of system of linear equations
    ssle->read (Gtm,in,Mespr);
    
    //  DOFs will be switched on or off with respect to time function
    Gtm->dofcontr=1;
    
    break;
  }
    // *********************************************************************
    //  GENERAL TIME-DEPENDENT MECHANICAL PROBLEMS WITH NEGLIGIBLE INERTIAL FORCES,
    //  GROWING STRUCTURES AND VARIABLE SETUP OF NONLINEAR SOLVER
    // *********************************************************************
  /*
  case mech_general_timedep_prob:{
    if (Mespr==1)  fprintf (stdout,"\n system of nonlinear equations will be solved by visco-solver");
    
    //  time controller
    timecon.read (in);

    // backup controler
    hdbcont.read(in);
    if ((Mespr==1) && hdbcont.hdbtype)
      fprintf(stdout, "\n backup will be performed %d", hdbcont.hdbtype);

    // reading setup of hdbackup, nonlinear managers and output in the important times
    if (timecon.nit == 0)   m = 1;
    else                    m = timecon.nit;
    nlman       = new nonlinman [m];
    hdbtime     = new answertype[m];
    outdrv_time = new answertype[m];
   

    //  maximum number of iterations in inner loop
    //  required norm of residuals
    xfscanf (in,"%k%ld%k%lf","nr_num_iter",&niilnr,"nr_error",&errnr);
    
    //  type of nonlinear solver
    xfscanf (in,"%k%m","nr_solvtype", &nonlintimesolvertype_kwdset, (int*)&nrsolv);

    //  type of storage of stiffness matrix
    //  type of stiffness matrix (initial, tangent)
    xfscanf (in,"%k%m%k%m","stiffmatstor", &storagetype_kwdset, (int*)&tstorsm, "stiffmat_type", &stiffmatrix_kwdset, (int*)&nlman->stmat);

    //  data about solver of system of linear equations
    ssle->read (Gtm,in,Mespr);

    break;
  }*/
    // **************************
    //  EARTH PRESSURE PROBLEMS
    // **************************
  case earth_pressure:{
    if (Mespr==1)  fprintf (stdout,"\n system of nonlinear equations will be solved by earth pressure solver");
    
    xfscanf(in, "%d", (int*)&tepsol);
    switch (tepsol){
    case gep_sol:{
      xfscanf (in,"%ld %lf %lf %lf",&niep,&errep,&stepep,&rlerrep);
      
      //  manager of solver of nonlinear equations
      nlman = new nonlinman[1];
      nlman->read (in,Mespr);
      
      //  type of storage of stiffness matrix
      xfscanf (in,"%k%m","stiffmatstor", &storagetype_kwdset, (int*)&tstorsm);
      //  data about solver of system of linear equations
      ssle->read (Gtm,in,Mespr);
      break;  
    }
    case gepvarsup_sol:{
      xfscanf (in,"%ld %lf %lf",&niep,&errep, &rlerrep);
      
      //  manager of solver of nonlinear equations
      nlman = new nonlinman[1];
      nlman->read (in,Mespr);
      
      xfscanf (in, "%ld", &nselnodep);
      if (niep < nselnodep)
        {
          fprintf(stderr, "\n\n Number of iteration steps is less then number of selected nodes");
          fprintf(stderr, "\n in file %s, line %d", __FILE__, __LINE__);
        }
      selnodep = new long [nselnodep];
      memset(selnodep, 0, sizeof(*selnodep)*nselnodep);
      for (i = 0; i < nselnodep; i++)
        {
          xfscanf (in, "%ld", selnodep+i);
          selnodep[i]--;
        }
      //  type of storage of stiffness matrix
      xfscanf (in,"%k%m","stiffmatstor", &storagetype_kwdset, (int*)&tstorsm);
      //  data about solver of system of linear equations
      ssle->read (Gtm,in,Mespr);
      break;
    }
    case epplast_sol:{
      xfscanf (in, "%ld", &niep);
      loadcoef = new double[niep];
      memset(loadcoef, 0, niep*sizeof(*loadcoef));
      for (i = 0; i < niep; i++)
	xfscanf(in, "%le", loadcoef+i);
      
      //  manager of solver of nonlinear equations
      nlman = new nonlinman[1];
      nlman->read (in,Mespr);
      
      //  type of storage of stiffness matrix
      //  type of stiffness matrix (initial, tangent)
      xfscanf (in,"%k%m%k%m","stiffmatstor", &storagetype_kwdset, (int*)&tstorsm, "stiffmat_type", &stiffmatrix_kwdset, (int*)&nlman->stmat);
      //  data about solver of system of linear equations
      ssle->read (Gtm,in,Mespr);
      break;
    }
    default:{
      print_err("uknown type of earth pressure is required",__FILE__,__LINE__,__func__);
    }
    }
    
    break;
  }
  // **********************************
  //  LAYERED LINEAR STATICS PROBLEMS
  // **********************************
  case layered_linear_statics:{
    if (Mespr==1)  fprintf (stdout,"\n linear layered static problem");
    
    //  tstorsm - type of storage of stiffness matrix
    xfscanf (in,"%d",(int*)&tstorsm);
    
    //  data about solver of system of linear equations
    ssle->read (Gtm,in,Mespr);
    
    break;
  }
    
  // *************************************
  //  LINEAR FLOATING SUBDOMAIN PROBLEMS
  // *************************************
  case lin_floating_subdomain:{
    if (Mespr==1)  fprintf (stdout,"\n linear floating subdomain problem");
    
    //  tstorsm - type of storage of stiffness matrix
    xfscanf (in,"%k%m","stiffmatstor",&storagetype_kwdset,(int*)&tstorsm);
    
    //  estimated number of rigid body motions (dimension of kernel of stiffness matrix)
    //xfscanf (in,"%d %le",&ense,&limit);
    
    //  data about solver of system of linear equations
    ssle->read (Gtm,in,Mespr);

    //  sequential topology will be read in function mechtop::read
    //Gtm->rst=1;

    break;
  }

  // *************************************
  //  NONLINEAR FLOATING SUBDOMAIN PROBLEMS
  // *************************************
  case nonlin_floating_subdomain:{
    if (Mespr==1)  fprintf (stdout,"\n nonlinear floating subdomain problem");
    
    //  data about nonlinear solver
    //nlman = new nonlinman[1];
    //nlman->read (in,Mespr);

    //  type of storage of stiffness matrix
    xfscanf (in,"%k%m","stiffmatstor", &storagetype_kwdset, (int*)&tstorsm);
    
    //  threshold for zero pivot detection
    //  the number of increments
    xfscanf (in,"%le %ld",&limit,&nincr);
    
    //  data about conjugate gradient method
    ssle->read (Gtm,in,Mespr);

    break;
  }

    // ****************************
    //  NONLINEAR FORCED DYNAMICS
    // ****************************
  case nonlinear_dynamics:{
    if (Mespr==1)  fprintf (stdout,"\n problem of nonlinear dynamics");
    
    //  type of solver of forced vibration problem
    xfscanf (in,"%k%m","forceddyn",&forcedsolver_kwdset,(int*)&tforvib);
    
    //  data about nonlinear solver
    nlman = new nonlinman[1];
    nlman->read (in,Mespr);

    switch (tforvib){
    case newmark:{
      xfscanf (in,"%lf %lf",&alphafvn,&deltafvn);
      break;
    }
      
    default:{
      print_err("unknown solver of forced vibration is required",__FILE__,__LINE__,__func__);
    }
    }
    
    //  type of damping
    xfscanf (in,"%k%m","damping",&damping_kwdset,(int*)&damp);
    switch (damp){
    case nodamp:{
      if (Mespr==1)  fprintf (stdout,"\n no damping is prescribed");
      break;
    }
    case massdamp:{
      xfscanf (in,"%lf",&dmass);
      if (Mespr==1)  fprintf (stdout,"\n damping proportional to the mass matrix is prescribed");
      break;
    }
    case stiffdamp:{
      xfscanf (in,"%lf",&dstiff);
      if (Mespr==1)  fprintf (stdout,"\n damping proportional to the stiffness matrix is prescribed");
      break;
    }
    case rayleighdamp:{
      xfscanf (in,"%lf %lf",&dstiff,&dmass);
      if (Mespr==1)  fprintf (stdout,"\n Rayleigh damping (proportional to the mass and stiffness matrices) is prescribed");
      break;
    }
    default:{
      print_err("unknown type of damping is required",__FILE__,__LINE__,__func__);
    }
    }
    
    //  type of storage of stiffness matrix
    //  type of storage of mass matrix
    //  type of storage of damping matrix is defined by type of damping
    xfscanf (in,"%k%m%k%m","stiffmatstor",&storagetype_kwdset,(int*)&tstorsm,"massmatstor",&storagetype_kwdset,(int*)&tstormm);
    
    //  time controller
    timecon.read (in);
    
    //  data about solver of system of linear equations
    ssle->read (Gtm,in,Mespr);
    
    // diagonalization of mass matrix
    xfscanf (in,"%k%ld","massdiagonalization",&diagmass);

    break;
  }
    
    // *******************************************
    //  PROBLEMS OF HEMIVARIATIONAL INEQUALITIES
    // *******************************************
  case hemivar_inequal:{
    if (Mespr==1)  fprintf (stdout,"\n problem of hemivariational inequalities");
    
    //  number of subdomains
    //xfscanf (in,"%ld",&nsub);

    //  tstorsm - type of storage of stiffness matrix
    xfscanf (in,"%k%m","stiffmatstor",&storagetype_kwdset,(int*)&tstorsm);
    
    //  data about solver of system of linear equations
    ssle->read (Gtm,in,Mespr);
    
    break;
  }
    // *****************
    //  LOAD BALANCING
    // *****************
  case load_balancing:{
    if (Mespr==1)  fprintf (stdout,"\n problem of load balancing");
    
    //  type of load balancing algorithm
    //  lbtype=1 - program stops after negm is known
    //  lbtype=0 - program goes on after negm is known
    xfscanf (in,"%ld",&lbtype);

    //  tstorsm - type of storage of stiffness matrix
    xfscanf (in,"%k%m","stiffmatstor",&storagetype_kwdset,(int*)&tstorsm);
    
    //  data about solver of system of linear equations
    ssle->read (Gtm,in,Mespr);
    
    //  code numbers will be generated with respect to the Schur complement method
    //  internal DOFs are ordered first, interface DOFs are ordered afterwards
    Gtm->cngen=2;

    //  sequential topology will be read in function mechtop::read
    Gtm->rst=1;

    break;
  }
    
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}




/**
  The function prints data for pre-processor.

  @param out - pointer to opened text file

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void probdesc::print (FILE *out)
{
  long i;

  //  problem name
  fprintf (out, "%s\n", name);

  //  detail messages
  fprintf (out, "%ld\n", Mespr);

  //  problem type
  fprintf (out, "%d\n", (int)tprob);
  
  //  computation of strains
  fprintf (out,"%ld ",straincomp);
  if (straincomp==1)  fprintf (out,"%ld %ld ",strainpos,strainaver);
  
  //  computation of stresses
  fprintf (out,"%ld ",stresscomp);
  if (stresscomp==1)  fprintf (out,"%ld %ld ",stresspos,stressaver);
  
  //  computation of values stored in array other
  fprintf (out,"%ld ",othercomp);
  if (othercomp==1)  fprintf (out,"%ld %ld ",otherpos,otheraver);

  //  computation of reactions
  fprintf (out,"%ld\n",reactcomp);
  
  // adaptivity calculation
  if (adaptivityflag)
    Ada->printinit(out);
  else
    fprintf (out,"0\n");
  
  //  stochastic or fuzzy computation
  fprintf (out,"%ld\n",stochasticcalc);
  
  //  homogenization procedures
  fprintf (out,"%ld\n",homog);
  
  if (homog==5 || homog==7){
    fprintf (out," %d\n", hstrastre);
  }
  if (homog==6 || homog==8){
    fprintf (out," %d\n", hstrastre);
  }
  if (homog==7 || homog==8){
    fprintf (out," %ld\n", homogdir);
  }

  //  renumbering of the nodes
  fprintf (out,"%d\n",Gtm->nodren);
  

  
  switch (tprob){
    // *****************
    //  LINEAR STATICS
    // *****************
  case linear_statics:{
    
    //  tstorsm - type of storage of stiffness matrix
    //  tlinsol - type of solver of linear equation system
    //fprintf (out, "%d %d\n", (int)tstorsm, (int)tlinsol);
    fprintf (out, "%d\n", (int)tstorsm);
    
    ssle->print (out);
    //print_linsoltype (out);
    
    fprintf(out, "\n");
    break;
  }
    
    // *********************************************************************************
    //  EIGENDYNAMICS - EIGENVALUES AND EIGENVECTORS (EIGENFREQUENCIES AND EIGENMODES)
    // *********************************************************************************
  case eigen_dynamics:{
    fprintf (out, "%d %d\n", (int)tstorsm, (int)tstormm);
    
    //  data about solver of eigenvalues and eigenvectors
    eigsol.print (out);

    //  data about solver of system of linear equations
    ssle->print (out);
    
    fprintf (out,"%ld\n",diagmass);
    
    break;
  }
    // *************************************************************************************
    //  LINEAR STABILITY - EIGENVALUES AND EIGENVECTORS (CRITICAL LOAD AND BUCKLING MODES)
    // *************************************************************************************
  case linear_stability:{
    fprintf (out, "%d %d\n", (int)tstorsm, (int)tstorgm);
    
    //  data about solver of eigenvalues and eigenvectors
    eigsol.print (out);

    //  data about solver of system of linear equations
    ssle->print (out);
    
    break;
  }


    // ******************
    //  FORCED DYNAMICS
    // ******************
  case forced_dynamics:
    {
      fprintf (out,"%d\n",(int)tforvib);
      
      switch (tforvib){
      case newmark:{
        fprintf (out, "%e %e\n",alphafvn, deltafvn);
        break;
      }
      case findiff:{
        break;
      }
      default:{
	print_err("unknown solver of forced vibration is required",__FILE__,__LINE__,__func__);
      }
      }
      
      
      fprintf (out,"%d\n",(int)damp);
      
      switch (damp){
      case nodamp:{
	break;
      }
      case massdamp:{
	fprintf (out,"%e\n",dmass);
	break;
      }
      case stiffdamp:{
	fprintf (out,"%e\n",dstiff);
	break;
      }
      case rayleighdamp:{
	fprintf (out,"%e %e\n",dstiff,dmass);
	break;
      }
      case elemdamp:{
	break;
      }
      default:{
	print_err("unknown type of damping is required",__FILE__,__LINE__,__func__);
      }
      }
      
      fprintf (out,"%d %d\n",(int)tstorsm,(int)tstormm);
      if (damp==elemdamp){
	fprintf (out,"%d\n",(int)tstordm);
      }
      
      timecon.print (out);

      ssle->print (out);

      fprintf (out,"%ld\n",diagmass);
      
     break;
    }
    // ********************
    //  NONLINEAR STATICS
    // ********************
  case mat_nonlinear_statics:
    {

      //  data about nonlinear solver
      nlman->print (out);
      
      // backup controler
      hdbcont.print (out);

      fprintf (out, "\n%d \n", (int)tstorsm);
      
      ssle->print (out);
      
      break;
    }
    // *********************************************************************
    //  TIME-DEPENDENT MECHANICAL PROBLEMS WITH NEGLIGIBLE INERTIAL FORCES
    // *********************************************************************
  case mech_timedependent_prob:{
    //  time controller
    timecon.print (out);

    //  type of time printing
    fprintf (out,"%d\n", tpr);
    
    // backup controler
    hdbcont.print (out);

    //  maximum number of iterations in inner loop
    //  required norm of residuals
    fprintf (out,"%ld %le %d\n",nlman->niilnr, nlman->errnr, nlman->rnormtnr);
    
    // checking of divergency of inner iteration loop
    fprintf (out,"%d\n",nlman->check_div);
    if (Mp->nlman->check_div == on){
      fprintf (out,"%ld",nlman->div_min_steps);
    }

    //  type of nonlinear solver
    //fprintf (out," %d\n",nrsolv);

    //  storage type of stiffness matrix
    //  type of stiffness matrix
    fprintf (out,"%d %d\n",(int)tstorsm,(int)nlman->stmat);
    if(nlman->stmat == ijth_tangent_stiff)
      fprintf (out,"%ld %ld\n", nlman->ithstiffchange, nlman->jthstiffchange);
    
    ssle->print (out);    
    break;
  }
    
    
    // **********************************************************************************************
    //  TIME-DEPENDENT MECHANICAL PROBLEMS WITH NEGLIGIBLE INERTIAL FORCES AND GRADUAL CONSTRUCTION
    // **********************************************************************************************
  case growing_mech_structure:{
    
    //  time controller
    timecon.print (out);

    //  type of time printing
    fprintf (out,"%d\n", tpr);
        
    // calculation of initial displacements due to interface nodes
    fprintf (out, "%d\n", comp_inidispl);

    // prescribed initial displacements by rotation of selected nodes
    if (comp_inidispl == yes)
      fprintf (out, "%d\n", rot_inidispl);

    // smooth element removal procedure (gradual application of traction due to element removal at interface nodes)
    fprintf (out, "%d\n", cpsmooth);
    if (cpsmooth == yes)
      // initial and minimum load coefficient increment in smoothed element removal procedure
      fprintf (out, "%le %le\n", cpincrnr, cpminincrnr);

    // backup controler
    hdbcont.print (out);

    //  maximum number of iterations in inner loop
    //  required norm of residuals
    fprintf (out,"%ld %le",nlman->niilnr,nlman->errnr);

    //  type of nonlinear solver
    fprintf (out," %d\n",nrsolv);

    fprintf (out,"%d %d\n",(int)tstorsm,(int)nlman->stmat);
    if(nlman->stmat == ijth_tangent_stiff)
      fprintf (out,"%ld %ld\n", nlman->ithstiffchange, nlman->jthstiffchange);
    
    ssle->print (out);
    
    break;
  }
    
    
    // **************************
    //  EARTH PRESSURE PROBLEMS
    // **************************
  case earth_pressure:{
    fprintf(out, "%d", (int)tepsol);
    switch (tepsol){
    case gep_sol:{
      fprintf (out,"%ld %le %le %le",niep,errep,stepep,rlerrep);
      
      nlman->print (out);
      
      //  type of storage of stiffness matrix
      fprintf (out,"%d ",(int)tstorsm);
      //  data about solver of system of linear equations
      ssle->print (out);
      break;  
    }
    case gepvarsup_sol:{
      fprintf (out,"%ld %e %e\n",niep,errep,rlerrep);
      
      nlman->print (out);
      
      fprintf (out, "%ld\n", nselnodep);
      if (niep < nselnodep)
	{
	  fprintf(stderr, "\n\n Number of iteration steps is less then number of selected nodes");
	  fprintf(stderr, "\n in file %s, line %d", __FILE__, __LINE__);
	}
      for (i = 0; i < nselnodep; i++)
	fprintf (out, "%ld ", selnodep[i]);
      //  type of storage of stiffness matrix
      fprintf (out,"\n%d ",(int)tstorsm);
      //  data about solver of system of linear equations
      ssle->print (out);
      break;  
    }
    case epplast_sol:{
      fprintf (out, " %ld\n", niep);
      for (i = 0; i < niep; i++)
	fprintf (out, "%f ", loadcoef[i]);
      
      nlman->print (out);
      
      //  type of storage of stiffness matrix
      //  type of stiffness matrix (initial, tangent)
      fprintf (out,"%d %d\n",(int)tstorsm,(int)nlman->stmat);
      //  data about solver of system of linear equations
      ssle->print (out);
      break;
    }
    default:{
      print_err("unknown type of earth pressure is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }
    // **********************************
  //  LAYERED LINEAR STATICS PROBLEMS
  // **********************************
  case layered_linear_statics:{

    //  tstorsm - type of storage of stiffness matrix
    //  tlinsol - type of solver of linear equation system
    //fprintf (out, "%d %d\n", (int)tstorsm, (int)tlinsol);
    fprintf (out, "%d\n", (int)tstorsm);
    
    ssle->print (out);
    //print_linsoltype (out);
    
    fprintf(out, "\n");
    
    break;
  }
    
    // *************************************
    //  LINEAR FLOATING SUBDOMAIN PROBLEMS
    // *************************************
  case lin_floating_subdomain:{
    
    fprintf (out, "%d\n", (int)tstorsm);
    
    //fprintf (out,"%ld %le\n",ense,limit);
    
    ssle->print (out);
    
    break;
  }
    // *************************************
    //  NONLINEAR FLOATING SUBDOMAIN PROBLEMS
    // *************************************
  case nonlin_floating_subdomain:{
    nlman->print (out);
    
    //fprintf (out, "%d %d\n", (int)tstorsm,(int)stmat);
    
    fprintf (out,"%ld %le\n",ense,limit);
    
    ssle->print (out);
    
    break;
  }
    
  case nonlinear_dynamics:{
    if (Mespr==1)  fprintf (stdout,"\n problem of nonlinear dynamics");
    
    
    break;
  }
    
    // *******************************************
    //  PROBLEMS OF HEMIVARIATIONAL INEQUALITIES
    // *******************************************
  case hemivar_inequal:{
    
    //  number of subdomains
    //fprintf (out,"%ld\n",nsub);

    //  tstorsm - type of storage of stiffness matrix
    fprintf (out, "%d\n", (int)tstorsm);
    
    //  data about solver of system of linear equations
    ssle->print (out);
    
    break;
  }
    // *****************
    //  LOAD BALANCING
    // *****************
  case load_balancing:{
    if (Mespr==1)  fprintf (stdout,"\n problem of load balancing");

    
    //  type of load balancing algorithm
    fprintf (out,"%ld\n",lbtype);

    //  tstorsm - type of storage of stiffness matrix
    fprintf (out,"%d\n",(int)tstorsm);
    
    //  data about solver of system of linear equations
    ssle->print (out);
    
    break;
  }
    
  default:{
    print_err("unknown type of problem is required",__FILE__,__LINE__,__func__);
    break;
  }
  }
  
  fprintf (out,"\n");
}
