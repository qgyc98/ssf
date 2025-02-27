#include "probdesct.h"
#include "globalt.h"
#include "aggregator.h"
#include <string.h>
#include "intools.h"
#include "iotools.h"
#include "gmatrix.h"
#include "adaptivityt.h"


probdesct::probdesct (void)
{
  long i;

  kwdsw = nokwd;
  
  //  number of transported media
  ntm=0;
  
  tmatt=nomedium;
  //  name of transported medium
  mednam=heat;

  //  array of dof indices of individual implemented primary variables(= nodal unknowns)
  for (i=0; i<tnkv; i++)
    var_dofid[i] = -1; // -1 = i-th variable is not defined

  // medium names associated with particualr DOFs
  dof_mednam = NULL;

  // string with particular medium names separated by \0
  dof_mednam_str = NULL;

  // array with ordering of used primary variables(= nodal unknowns)
  dofname = NULL;
  
  //  smaller number of integration points
  savemode=1;
  
  //  type of solver of nonstationary problems
  tnpsolver = nonstatsolver(0);
  
  //  type of transport solver
  trsolv=fullnewtont;
  //trsolv=modnewtont;

  //  type of residuum computation
  trestype=lrhst;
  //trestype=fluxest;

  //  convergence control for fully newton-raphson
  convergcontrolt=no;
  
  tstorkm = (storagetype) 0; 
  tstorcm = (storagetype) 0; 
  
  adaptivityflag=0;
  
  limit=0.0;  zero=1.0e-20;
  
  time=0.0;  alpha=1.0;
  istep = 0;
  jstep = 0;
  
  err=0.0;  nii=0;
  errarr = NULL;
  threshrhs = NULL;

  stochasticcalc=0;
  
  //  capacity matrix is not diagonalized
  diagcap=0;
  //  reaction matrix is not diagonalized
  diagreact=1;
  //  reaction/capacity computation
  //  reaction and capacity matrices are computed by the same function
  //  default value leads to computation of the capacity matrix
  react_capac=capacity;
  
  //  presence of advection
  advect = 0L;
  //  presence of reaction term
  react = 0L;
  
  gdim=0;  
  scale=NULL;

  tgravity = gr_no;
  memset(gr, 0, sizeof(*gr)*3);

  // integration point values have not been computed
  ipvcomp = 0;
  //  gradients are not computed and stored
  gradcomp = 0;
  //  fluxes are not computed and stored
  fluxcomp = 0;
  //  components of array other are not computed and stored
  othercomp = 0;
  //  components of array eqother are computed and stored
  eqothercomp = 0;
  
  //  gradients are not averaged
  gradaver = 0;
  //  fluxes are not averaged
  fluxaver = 0;
  //  array other is not averaged
  otheraver = 0;
  //  array eqother is not averaged
  eqotheraver = 0;

  //  gradients position
  gradpos = 0;
  //  fluxes position
  fluxpos = 0;
  //  array other position
  otherpos = 0;
  //  array eqother position
  eqotherpos = 0;
  
  //  actual nodal values are not stored on node
  nvs=0;
  //  previous nodal values are not stored on node
  pnvs=0;
  //  initial nodal values are not stored on node
  invs=0;
  //  time derivatives of actual nodal values are not stored on node
  tdnvs=0;
  
  //  homogenization
  homogt=0;
  //  type of macro-micro problem correspondence in homogenization
  mami = (macromicrotype) 0;
  
  memset (name,0,sizeof(*name)*1001);
  path = filename = suffix = NULL;

  //  type of time printing
  //tprt = dayst;
  tprt = secondst;
  //tprt = minutest;

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

probdesct::~probdesct (void)
{
  delete [] errarr;
  delete [] threshrhs;
  delete [] scale;

  delete [] path; 
  delete [] filename; 
  delete [] suffix;
  delete ssle;
  delete [] dofname;
  delete [] dof_mednam;
  delete [] dof_mednam_str;
}


/**  function reads basic data about solved problem
     JK, 25.9.2001
     
     @param in - input file
*/
void probdesct::read (XFILE *in)
{
  long i, l;
  char *amnstr;

  //  problem name
  xfscanf (in,"% a",name);
  
  //  detail messages
  Mesprt=0;
  xfscanf (in,"%k%ld","mesprt",&Mesprt);
  if (Mesprt==1)  fprintf (stdout,"\n detail information will be printed");
  else fprintf (stdout,"\n only important messages will be printed");
  
  //  problem type
  xfscanf (in,"%k%m","problemtype",&problemtypet_kwdset,(int*)&tprob);
  
  //  number of transported matters
  xfscanf (in,"%k%m","transmatter",&transmattert_kwdset,(int*)&tmatt);
  switch (tmatt){
  case onemedium:{        ntm=1;  break; }
  case twomediacoup:{     ntm=2;  break; }
  case threemediacoup:{   ntm=3;  break; }
  case fourmediacoup:{    ntm=4;  break; }
  default:{
    print_err("unknown type of transported matter is required",__FILE__,__LINE__,__func__);
  }
  }
  if (Mesprt==1)
    fprintf(stdout, "\n the number of transported media is %ld", ntm);

  //  names of transported media
  xfscanf (in,"%k%m","mednames",&mednamest_kwdset,(int*)&mednam);
  
  // DOF vs. media name map
  dof_mednam = new char*[ntm];
  memset(dof_mednam, 0, sizeof(*dof_mednam)*ntm);
  l = long(strlen(mednamest_kwdset.get_str(mednam)));
  // one string with media names separated by \0 
  dof_mednam_str = new char[l+1];
  strncpy(dof_mednam_str, mednamest_kwdset.get_str(mednam), l+1);
  // initial media name string parsing 
  // amnstr is the pointer to the actual medium name string
  amnstr = strtok(dof_mednam_str, "_");
  i = 0;

  /* do
     {
     if (amnstr == NULL) // no more medium name string
     {
     print_err("medium name '%s' violates the naming convention.\n"
     " The number of '_' characters in the medium name must be %ld (ntm-1)",
     __FILE__, __LINE__, __func__, mednamest_kwdset.get_str(mednam), ntm-1);
     }
     dof_mednam[i] = amnstr;    
     amnstr = strtok(NULL, "_");
     i++;
     } while (i<ntm);
  */

  // scales of transported media
  scale = new double[ntm];
  memset(scale, 0, sizeof(*scale)*ntm);

  // array of ordered dof names
  dofname = new namevart[ntm];
  memset(dofname, 0, sizeof(*dofname)*ntm);

  xfscanf(in,"%k","scale");
  for (i=0; i<ntm; i++)
    xfscanf (in,"%lf",scale+i);

  // *************************
  //  computation of gradients
  // *************************
  xfscanf (in,"%k%ld","gradcomp",&gradcomp);
  if (gradcomp!=0 && gradcomp!=1)
    print_err("\n\n wrong identification of gradient computation", __FILE__,__LINE__, __func__);
  if (gradcomp==1) {
    xfscanf (in,"%k%m%k%ld","gradpos", &elempositiont_kwdset, &gradpos, "gradaver", &gradaver);
    
    if (gradpos!=1 && gradpos!=2 && gradpos!=3){
      print_err("wrong definition of position where gradients are computed", __FILE__, __LINE__,__func__);
    }
    if (gradaver!=0 && gradaver!=1 && gradaver!=2)
      print_err("\n\n wrong identification of gradient averaging",__FILE__,__LINE__,__func__);
  }
  
  if (Mesprt==1){
    if (gradcomp==0)  fprintf (stdout,"\n gradients will not be computed and stored");
    else{
      fprintf (stdout,"\n gradients will be computed and stored");
      if (gradpos==1)   fprintf (stdout,"\n gradients will be computed at integration points");
      if (gradpos==2)   fprintf (stdout,"\n gradients will be copied to nodes from the closest integration point");
      if (gradpos==3)   fprintf (stdout,"\n gradients will be computed at nodes");
      
      if (gradaver==0)  fprintf (stdout,"\n gradients will not be averaged");
      if (gradaver==1)  fprintf (stdout,"\n gradients will be averaged at nodes by the number of contributions");
      if (gradaver==2)  fprintf (stdout,"\n gradients will be averaged at nodes by volume");
    }
  }
  

  // *************************
  //  computation of fluxes
  // *************************
  xfscanf (in,"%k%ld","fluxcomp",&fluxcomp);
  if (fluxcomp!=0 && fluxcomp!=1)
    print_err("\n\n wrong identification of flux computation", __FILE__,__LINE__, __func__);
  if (fluxcomp==1) {
    xfscanf (in,"%k%m%k%ld","fluxpos", &elempositiont_kwdset, &fluxpos, "fluxaver", &fluxaver);
    
    if (fluxpos!=1 && fluxpos!=2 && fluxpos!=3){
      print_err("wrong definition of position where fluxes are computed", __FILE__, __LINE__,__func__);
    }
    if (fluxaver!=0 && fluxaver!=1 && fluxaver!=2)
      print_err("\n\n wrong identification of flux averaging",__FILE__,__LINE__,__func__);
  }
  
  if (Mesprt==1){
    if (fluxcomp==0)  fprintf (stdout,"\n fluxes will not be computed and stored");
    else{
      fprintf (stdout,"\n fluxes will be computed and stored");
      if (fluxpos==1)   fprintf (stdout,"\n fluxes will be computed at integration points");
      if (fluxpos==2)   fprintf (stdout,"\n fluxes will be copied to nodes from the closest integration point");
      if (fluxpos==3)   fprintf (stdout,"\n fluxes will be computed at nodes");
      
      if (fluxaver==0)  fprintf (stdout,"\n fluxes will not be averaged");
      if (fluxaver==1)  fprintf (stdout,"\n fluxes will be averaged at nodes by the number of contributions");
      if (fluxaver==2)  fprintf (stdout,"\n fluxes will be averaged at nodes by volume");
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

    if (otherpos!=1 && otherpos!=2 && otherpos!=3){
      print_err("wrong definition of position where other are computed", __FILE__, __LINE__, __func__);
    }
    if (otheraver!=0 && otheraver!=1 && otheraver!=2)
      print_err("wrong definition of other averaging", __FILE__, __LINE__, __func__);
  }
  
  if (Mesprt==1){
    if (othercomp==1){
      fprintf (stdout,"\n other will be computed and stored");
      if (otherpos==1)   fprintf (stdout,"\n other will be computed at integration points");
      if (otherpos==2)   fprintf (stdout,"\n other will be copied to nodes from the closest integration point");
      if (otherpos==3)   fprintf (stdout,"\n other will be computed at nodes");
      
      if (otheraver==0)  fprintf (stdout,"\n other will not be averaged");
      if (otheraver==1)  fprintf (stdout,"\n other will be averaged at nodes by the number of contributions");
      if (otheraver==2)  fprintf (stdout,"\n other will be averaged at nodes by volume");
    }
    else
      fprintf (stdout,"\n other will not be computed and stored");
  }

    
  // ****************************
  //  computation of eqother values
  // ****************************
  xfscanf (in,"%k%ld","eqothercomp", &eqothercomp);
  if (eqothercomp!=0 && eqothercomp!=1)
    print_err("wrong definition of eqother values computation", __FILE__, __LINE__, __func__);
  if (eqothercomp==1){
    xfscanf (in,"%k%ld%k%ld","eqotherpos", &eqotherpos, "eqotheraver", &eqotheraver);

    if (eqotherpos!=1 && eqotherpos!=2 && eqotherpos!=3){
      print_err("wrong definition of position where eqother are computed", __FILE__, __LINE__, __func__);
    }
    if (eqotheraver!=0 && eqotheraver!=1 && eqotheraver!=2)
      print_err("wrong definition of eqother averaging", __FILE__, __LINE__, __func__);
  }
  
  if (Mesprt==1){
    if (eqothercomp==1){
      fprintf (stdout,"\n eqother will be computed and stored");
      if (eqotherpos==1)   fprintf (stdout,"\n eqother will be computed at integration points");
      if (eqotherpos==2)   fprintf (stdout,"\n eqother will be copied to nodes from the closest integration point");
      if (eqotherpos==3)   fprintf (stdout,"\n eqother will be computed at nodes");
      
      if (eqotheraver==0)  fprintf (stdout,"\n eqother will not be averaged");
      if (eqotheraver==1)  fprintf (stdout,"\n eqother will be averaged at nodes by the number of contributions");
      if (eqotheraver==1)  fprintf (stdout,"\n eqother will be averaged at nodes by volume");
    }
    else
      fprintf (stdout,"\n eqother will not be computed and stored");
  }
  
  //  gravity acceleration is taken into account
  xfscanf (in,"%k%m","gravityacceleration",&gravityaccelerationt_kwdset,(int*)&tgravity);
  
  if(tgravity == gr_yes){
    xfscanf (in,"%lf %lf %lf",gr+0,gr+1,gr+2);
  }
  
  // **************************
  //   adaptivity computation
  // **************************
  xfscanf (in,"%k%ld","adaptivity",&adaptivityflag);
  if (adaptivityflag) {
    Adat = new adaptivityt();
    Adat->readinit(in);
  }
  
  // **********************************
  //  stochastic or fuzzy computation
  // **********************************
  //  stochasticcalc=0 - deterministic computation
  //  stochasticcalc=1 - stochastic/fuzzy computation, data are read all at once
  //  stochasticcalc=2 - stochastic/fuzzy computation, data are read sequentially
  //  stochasticcalc=3 - stochastic/fuzzy computation, data are generated in the code
  xfscanf (in,"%k%ld","stochasticcalc",&stochasticcalc);
  
  // ****************************
  //  homogenization procedures
  // ****************************
  //  homogt=0 - no homogenization
  //  homogt=1 - homogenization on a single processor
  //  homogt=2 - homogenization on a parallel computer (many processors are used)
  xfscanf (in,"%k%ld","homogenization",&homogt);
  if (Mesprt==1){
    if (homogt==0)
      fprintf (stdout,"\n homogenization will not be computed");
    if (homogt!=0)
      fprintf (stdout,"\n homogenization will be computed");
  }
  if (homogt==2){
    //  type of macro-micro problem correspondence
    //  mami=1 - elements are connected with microproblems
    //  mami=2 - aggregates of elements are connected with microproblems
    xfscanf (in,"%k%ld","macmiccorres",&mami);
    Gtt->rst=1;
  }
  
  
  // ***************************
  //  renumbering of the nodes
  // ***************************
  xfscanf (in,"%k%m","noderenumber",&noderenumb_kwdset,(int*)&Gtt->nodren);
  if (Mesprt==1){
    switch (Gtt->nodren){
    case no_renumbering:{
      fprintf (stdout,"\n nodes will not be renumbered");
      break;
    }
    case rev_cuthill_mckee:{
      fprintf (stdout,"\n nodes will be renumbered by reverse Cuthill-McKee algorithm");
      break;
    }
    default:{
      print_err("unknown type of node reordering is required",__FILE__,__LINE__,__func__);
    }
    }
  }
  
  // ****************************
  //  presence of advection terms
  // ****************************
  //  advect=0 - no advection
  //  advect=1 - there is advection described by velocity vector
  xfscanf (in,"%k%ld","advectioncomp",&advect);
  if (Mesprt==1){
    switch (advect){
    case 0:{
      fprintf (stdout,"\n advection terms will not be assumed");
      break;
    }
    case 1:{
      fprintf (stdout,"\n advection terms  will be assumed");
      break;
    }
    default:{
      print_err("unknown type of advection computation is required",__FILE__,__LINE__,__func__);
    }
    }
  }

  // ****************************
  //  presence of reaction term
  // ****************************
  //  react=0 - no reaction
  //  react=1 - there is reaction term
  xfscanf (in,"%k%ld","reactioncomp",&react);
  if (Mesprt==1){
    switch (react){
    case 0:{
      fprintf (stdout,"\n reaction terms will not be assumed");
      break;
    }
    case 1:{
      fprintf (stdout,"\n reaction terms will be assumed");
      break;
    }
    default:{
      print_err("unknown type of reaction computation is required",__FILE__,__LINE__,__func__);
    }
    }
  }
  
  switch (tprob){
    // ****************************
    //  LINEAR STATIONARY PROBLEM
    // ****************************
  case stationary_problem:{
    if (Mesprt==1)  fprintf (stdout,"\n linear stationary problem");
    
    //  tstorsm - type of storage of conductivity matrix
    xfscanf (in,"%k%m","conductmatstor",&storagetype_kwdset,(int*)&tstorkm);
    
    //  data about solver of system of linear equations
    ssle->read (Gtt,in,Mesprt);
    
    break;
  }

    // ****************************
    //  NONLINEAR STATIONARY PROBLEM
    // ****************************
  case nonlinear_stationary_problem:{
    if (Mesprt==1)  fprintf (stdout,"\n nonlinear stationary problem");
    
    //  data about nonlinear solver
    nlman.read (in,Mesprt);
    
    //  type of residuum computation
    xfscanf (in,"%k%m","residuumcomputationtype",&transpresiduumtype_kwdset,(int*)&trestype);

    //  tstorsm - type of storage of conductivity matrix
    xfscanf (in,"%k%m","conductmatstor",&storagetype_kwdset,(int*)&tstorkm);
    
    //  data about solver of system of linear equations
    ssle->read (Gtt,in,Mesprt);
    
    break;
  }

    // *******************************
    //  LINEAR NONSTATIONARY PROBLEM
    // *******************************
  case nonstationary_problem:{
    if (Mesprt==1)  fprintf (stdout,"\n linear nonstationary problem");
    
    //  type of solver of nonstationary problem
    //  1 - trapezoidal method, d-version
    //  2 - trapezoidal method, v-version
    //  11 - forward finite differences
    xfscanf (in,"%k%m","nonstatsolver",&nonstatsolver_kwdset,(int*)&tnpsolver);
    
    switch (tnpsolver){
    case trapezoidd:{
      xfscanf (in,"%k%lf","alpha_integration", &alpha);
      break;
    }
    case trapezoidv:{
      xfscanf (in,"%k%lf","alpha_integration", &alpha);
      break;
    }
    case optim_trapezoidd:{
      xfscanf (in,"%k%lf","alpha_integration", &alpha);
      break;
    }
    case forwardfindif:{
      break;
    }
    default:{
      print_err("unknown type of nonstationary solver is required", __FILE__, __LINE__, __func__);      
    }
    }

    //  time controller
    timecont.read (in);

    xfscanf (in,"%k%m","timetypeprint", &timetypeprint_kwdset, &tprt);

    // backup controler
    hdbcont.read(in);
    if ((Mesprt==1) && hdbcont.hdbtype)
      fprintf(stdout, "\n backup will be performed %d", hdbcont.hdbtype);


    
    xfscanf (in,"%k%m%k%m","conductmatstor",&storagetype_kwdset,(int*)&tstorkm,"capacmatstor",&storagetype_kwdset,(int*)&tstorcm);
    
    //  data about solver of system of linear equations
    ssle->read (Gtt,in,Mesprt);
    
    // diagonalization of capacity matrix
    xfscanf (in,"%k%ld","capdiagonalization",&diagcap);
    if (react>0){
      // diagonalization of reaction matrix
      xfscanf (in,"%k%ld","reactdiagonalization",&diagreact);
    }
    
    break;
  }
    
    // ***************************************
    //  LINEAR NONSTATIONARY GROWING PROBLEM
    // ***************************************
  case growing_np_problem:{
    if (Mesprt==1)  fprintf (stdout,"\n nonstationary problem with growing number of elements");
    
    //  time controller
    timecont.read (in);

    xfscanf (in,"%k%m","timetypeprint", &timetypeprint_kwdset, &tprt);

    // backup controler
    hdbcont.read(in);
    if ((Mesprt==1) && hdbcont.hdbtype)
      fprintf(stdout, "\n backup will be performed %d", hdbcont.hdbtype);

    xfscanf (in,"%k%lf","alpha_integration", &alpha);

    xfscanf (in,"%k%m%k%m","conductmatstor",&storagetype_kwdset,(int*)&tstorkm,"capacmatstor",&storagetype_kwdset,(int*)&tstorcm);
    
    //  data about solver of system of linear equations
    ssle->read (Gtt,in,Mesprt);
  
    // diagonalization of capacity matrix
    xfscanf (in,"%k%ld","capdiagonalization",&diagcap);
    if (react>0){
      // diagonalization of reaction matrix
      xfscanf (in,"%k%ld","reactdiagonalization",&diagreact);
    }
    
    //  degrees of freedom are controlled by a time function
    Gtt->dofcontr=1;
    
    //  this problem requires storage of actual nodal values on nodes
    nvs=1;
    //  this problem requires storage of initial nodal values on nodes
    invs=1;
    //  this problem requires storage of time derivatives of nodal values on nodes
    tdnvs=1;
    
    break;
  }
    
    // **********************************
    //  NONLINEAR NONSTATIONARY PROBLEM
    // **********************************
  case nonlinear_nonstationary_problem:{
    if (Mesprt==1)  fprintf (stdout,"\n nonlinear nonstationary problem");
    
    //  time controller
    timecont.read (in);

    xfscanf (in,"%k%m","timetypeprint", &timetypeprint_kwdset, &tprt);
    fprintf(stdout, "\n time printing will be in %d", tprt);
    
    // backup controler
    hdbcont.read(in);
    if ((Mesprt==1) && hdbcont.hdbtype)
      fprintf(stdout, "\n backup will be performed %d", hdbcont.hdbtype);
    else
      fprintf(stdout, "\n backup will not be performed = %d", hdbcont.hdbtype);

    errarr = new double [ntm];
    threshrhs = new double [ntm];

    //xfscanf (in,"%k%lf %lf %lf %ld","alpha_integration", &alpha,&err,&threshrhs,&nii);
    xfscanf (in,"%k%lf %ld","alpha_integration", &alpha,&nii);

    for (i=0;i<ntm;i++){
      xfscanf (in,"%lf %lf",errarr+i,threshrhs+i);
    }

    
    //  type of transport solver
    xfscanf (in,"%k%m","transportsolvertype",&transpsolver_kwdset,(int*)&trsolv);
    
    //  type of residuum computation
    xfscanf (in,"%k%m","residuumcomputationtype",&transpresiduumtype_kwdset,(int*)&trestype);
    
    //  convergence control for fully newton-raphson
    xfscanf (in,"%k%m","convergencecontrol",&answertype_kwdset,(int*)&convergcontrolt);


    xfscanf (in,"%k%m%k%m","conductmatstor",&storagetype_kwdset,(int*)&tstorkm,"capacmatstor",&storagetype_kwdset,(int*)&tstorcm);
    
    //  data about solver of system of linear equations
    ssle->read (Gtt,in,Mesprt);
    
    // diagonalization of capacity matrix
    xfscanf (in,"%k%ld","capdiagonalization",&diagcap);
    if (react>0){
      // diagonalization of reaction matrix
      xfscanf (in,"%k%ld","reactdiagonalization",&diagreact);
    }
    
    break;
  }
    
    // ***************************************
    //  NONLINEAR NONSTATIONARY GROWING PROBLEM
    // ***************************************
  case growing_np_problem_nonlin:{
    if (Mesprt==1)  fprintf (stdout,"\n nonlinear nonstationary problem with growing number of elements");
    
    //  time controller
    timecont.read (in);

    xfscanf (in,"%k%m","timetypeprint", &timetypeprint_kwdset, &tprt);

    // backup controler
    hdbcont.read(in);
    if ((Mesprt==1) && hdbcont.hdbtype)
      fprintf(stdout, "\n backup will be performed %d", hdbcont.hdbtype);

    //xfscanf (in,"%k%lf %lf %lf %ld","alpha_integration",&alpha,&err,&threshrhs,&nii);
    xfscanf (in,"%k%lf %ld","alpha_integration", &alpha,&nii);

    errarr = new double [ntm];
    threshrhs = new double [ntm];

    for (i=0;i<ntm;i++){
      xfscanf (in,"%lf %lf",errarr+i,threshrhs+i);
    }

    //  type of transport solver
    xfscanf (in,"%k%m","transportsolvertype",&transpsolver_kwdset,(int*)&trsolv);
    
    //  type of residuum computation
    xfscanf (in,"%k%m","residuumcomputationtype",&transpresiduumtype_kwdset,(int*)&trestype);
    
    //  convergence control for fully newton-raphson
    xfscanf (in,"%k%m","convergencecontrol",&answertype_kwdset,(int*)&convergcontrolt);

    xfscanf (in,"%k%m%k%m","conductmatstor",&storagetype_kwdset,(int*)&tstorkm,"capacmatstor",&storagetype_kwdset,(int*)&tstorcm);
    
    //  data about solver of system of linear equations
    ssle->read (Gtt,in,Mesprt);
  
    // diagonalization of capacity matrix
    xfscanf (in,"%k%ld","capdiagonalization",&diagcap);
    if (react>0){
      // diagonalization of reaction matrix
      xfscanf (in,"%k%ld","reactdiagonalization",&diagreact);
    }
    
    //  degrees of freedom are controlled by a time function
    Gtt->dofcontr=1;
    
    //  this problem requires storage of actual nodal values on nodes
    nvs=1;
    //  this problem requires storage of initial nodal values on nodes
    invs=1;
    //  this problem requires storage of time derivatives of nodal values on nodes
    tdnvs=1;
 
    break;
  }
    
    // *******************************************
    //  NONSTATIONARY PROBLEM WITH DISCONTINUITY
    // *******************************************
  case discont_nonstat_problem:{
    if (Mesprt==1)  fprintf (stdout,"\n discontinuous linear_nonstationary problem");
    
    //  time controller
    timecont.read (in);
    
    xfscanf (in,"%k%m","timetypeprint", &timetypeprint_kwdset, &tprt);
 
    // backup controler
    hdbcont.read(in);
    if ((Mesprt==1) && hdbcont.hdbtype)
      fprintf(stdout, "\n backup will be performed %d", hdbcont.hdbtype);

    xfscanf (in,"%k%lf","alpha_integration",&alpha);
    
    xfscanf (in,"%k%m%k%m","conductmatstor",&storagetype_kwdset,(int*)&tstorkm,"capacmatstor",&storagetype_kwdset,(int*)&tstorcm);
    
    //  data about solver of system of linear equations
    ssle->read (Gtt,in,Mesprt);
  
    // diagonalization of capacity matrix
    xfscanf (in,"%k%ld","capdiagonalization",&diagcap);
    if (react>0){
      // diagonalization of reaction matrix
      xfscanf (in,"%k%ld","reactdiagonalization",&diagreact);
    }
    
    Gtt->edtype=jumps;

    break;
  }
  case discont_nonlin_nonstat_problem:{
    if (Mesprt==1)  fprintf (stdout,"\n discontinuous nonlinear_nonstationary problem");
    
    //  time controller
    timecont.read (in);

    xfscanf (in,"%k%m","timetypeprint", &timetypeprint_kwdset, &tprt);
 
    // backup controler
    hdbcont.read(in);
    if ((Mesprt==1) && hdbcont.hdbtype)
      fprintf(stdout, "\n backup will be performed %d", hdbcont.hdbtype);

    errarr = new double [ntm];
    threshrhs = new double [ntm];
    
    xfscanf (in,"%k%lf %ld","alpha_integration",&alpha,&nii);
    
    for (i=0;i<ntm;i++){
      xfscanf (in,"%lf %lf",errarr+i,threshrhs+i);
    }
    
    xfscanf (in,"%k%m%k%m","conductmatstor",&storagetype_kwdset,(int*)&tstorkm,"capacmatstor",&storagetype_kwdset,(int*)&tstorcm);
    
    //  data about solver of system of linear equations
    ssle->read (Gtt,in,Mesprt);
  
    // diagonalization of capacity matrix
    xfscanf (in,"%k%ld","capdiagonalization",&diagcap);
    if (react>0){
      // diagonalization of reaction matrix
      xfscanf (in,"%k%ld","reactdiagonalization",&diagreact);
    }
    
    Gtt->edtype=jumps;

    break;
  }
    
    // ***************************************
    //  provisional connection with Hermes2D
    // ***************************************
  case hermes:{
    
    break;
  }

  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  // **********************************
  //  special requirements of solvers
  // **********************************
  //  JK, oprava 10.8.2010
  //if (ssle->tlinsol==sschur){
  //  sequential topology will be read in function mechtop::read
  //Gtt->rst=1;
  //Gtt->cngen=0;
  // }
  //if (ssle->tlinsol==sfeti){
  //  sequential topology will be read in function transtop::read
  //Gtt->rst=1;
  //
  //if (ssle->feti.prec==feti_dirichlet){
  //  code numbers will be generated with respect to the Schur complement method
  //  internal DOFs are ordered first, interface DOFs are ordered afterwards
  //Gtt->cngen=2;
  //}
  //}
  //if (ssle->tlinsol==sfetidef){
  //  sequential topology will be read in function transtop::read
  //Gtt->rst=2;
  //}
  //
  //if (ssle->prec.pt==boss){
  //if (ssle->prec.agg->impl==2){
  //  sequential topology will be read in function transtop::read
  //  Gtt->rst=1;
  //}
  //}
  //
  //  JK, konec opravy 10.8.2010
  
}

/**  function prints basic data about solved problem
     JK, 25.9.2001
     
     @param in - output file
*/
void probdesct::print (FILE *out)
{
  long i;

  //  problem name
  fprintf (out, "%s\n", name);

  //  detail messages
  fprintf (out,"%ld\n", Mesprt);

  //  problem type
  fprintf (out,"%d", tprob);

  //  transported matter
  fprintf (out," %d", tmatt);
  
  //  names of transported media
  fprintf (out, " %d", mednam);
  
  // scale of transported media
  for (i=0; i<ntm; i++)
    fprintf(out, " %e", scale[i]);
  fprintf(out, "\n");

  // *************************
  //  computation of gradients
  // *************************
  fprintf (out,"%ld # gradients\n",gradcomp);
  if (gradcomp==1) {
    fprintf (out,"%ld %ld\n",gradpos,gradaver);
  }
  // *************************
  //  computation of fluxes
  // *************************
  fprintf (out,"%ld # fluxes\n",fluxcomp);
  if (fluxcomp==1) {
    fprintf (out,"%ld %ld\n",fluxpos,fluxaver);
  }
  // ****************************
  //  computation of other values
  // ****************************
  fprintf (out,"%ld # other values\n",othercomp);
  if (othercomp==1){
    fprintf (out,"%ld %ld\n",otherpos,otheraver);
  }
  // ****************************
  //  computation of eqother values
  // ****************************
  fprintf (out,"%ld # eqother values\n",eqothercomp);
  if (eqothercomp==1){
    fprintf (out,"%ld %ld\n",eqotherpos,eqotheraver);
  }
  
  //  gravity acceleration
  fprintf (out,"%d # gravity\n",tgravity);
  if(tgravity == gr_yes){
    fprintf (out,"%le %le %le\n",gr[0],gr[1],gr[2]);
  }
  
  // adaptivity calculation
  if (adaptivityflag)
    Adat->printinit(out);
  else
    fprintf (out,"0 # adaptivity\n");
  
  
  // stochastic calculation
  fprintf (out,"%ld # deterministic/stochastic/fuzzy\n",stochasticcalc);

  
  //  homogenization procedures
  fprintf (out,"%ld # homogenization\n",homogt);


  // ***************************
  //  renumbering of the nodes
  // ***************************
  fprintf (out,"%d # node renumbering\n",Gtt->nodren);

  // ************************
  //  presence of advection
  // ************************
  fprintf (out,"%ld # advection\n",advect);

  // ************************
  //  presence of advection
  // ************************
  fprintf (out,"%ld # reaction\n",react);

  switch (tprob){
    // ****************************
    //  LINEAR STATIONARY PROBLEM
    // ****************************
  case stationary_problem:{
    fprintf (out,"%d\n", tstorkm);
    
    ssle->print (out);

    break;
  }
    
    // ****************************
    //  NONLINEAR STATIONARY PROBLEM
    // ****************************
  case nonlinear_stationary_problem:{
    
    nlman.print (out);
    
    fprintf (out,"%d ", trestype);

    fprintf (out,"%d\n", tstorkm);
    
    ssle->print (out);

    break;
  }
    
    // *******************************
    //  LINEAR NONSTATIONARY PROBLEM
    // *******************************
  case nonstationary_problem:{
    
    fprintf (out,"%d\n", tnpsolver);

    switch (tnpsolver){
    case trapezoidd:{
      fprintf (out, "%le\n",alpha);
      break;
    }
    case trapezoidv:{
      fprintf (out, "%le\n",alpha);
      break;
    }
    case optim_trapezoidd:{
      fprintf (out, "%le\n",alpha);
      break;
    }
    case forwardfindif:{
      break;
    }
    default:{
      print_err("unknown type of nonstationary solver is required", __FILE__, __LINE__, __func__);      
    }
    }
    
    timecont.print (out);

    fprintf (out,"%d\n", tprt);

    hdbcont.print(out);
    
    fprintf (out,"%d %d ",tstorkm,tstorcm);
    
    ssle->print (out);

    fprintf (out,"%ld\n",diagcap);
    if (react>0)
      fprintf (out,"%ld\n",diagreact);
    
    fprintf (out,"\n");    
    
    break;
  }

    // ***************************************
    //  LINEAR NONSTATIONARY GROWING PROBLEM
    // ***************************************
  case growing_np_problem:{
    
    timecont.print (out);

    fprintf (out,"%d\n", tprt);

    hdbcont.print(out);

    fprintf (out, "%le\n",alpha);
    
    fprintf (out,"%d %d ",tstorkm,tstorcm);
    
    ssle->print (out);

    fprintf (out,"%ld\n",diagcap);
    if (react>0)
      fprintf (out,"%ld\n",diagreact);

    fprintf (out,"\n");    
    
    break;
  }

    // **********************************
    //  NONLINEAR NONSTATIONARY PROBLEM
    // **********************************
  case nonlinear_nonstationary_problem:{
    
    timecont.print (out);

    fprintf (out,"%d\n", tprt);

    hdbcont.print(out);

    fprintf (out, "%le %ld\n", alpha,nii);
    
    for (i=0;i<ntm;i++){
      fprintf (out, "%le %le\n", errarr[i], threshrhs[i]);
    }
    
    fprintf (out,"%d %d %d\n",trsolv,trestype,convergcontrolt);

    fprintf (out,"%d %d \n",tstorkm,tstorcm);

    ssle->print (out);

    fprintf (out,"%ld\n",diagcap);
    if (react>0)
      fprintf (out,"%ld\n",diagreact);

    fprintf (out,"\n");    

    break;
  }

    // ***************************************
    //  NONLINEAR NONSTATIONARY GROWING PROBLEM
    // ***************************************
  case growing_np_problem_nonlin:{
    
    timecont.print (out);

    fprintf (out,"%d\n", tprt);

    hdbcont.print(out);

    fprintf (out, "%le %ld\n",alpha,nii);
    for (i=0;i<ntm;i++){
      fprintf (out, "%le %le\n",errarr[i],threshrhs[i]);
    }
    
    fprintf (out,"%d %d %d\n",trsolv,trestype,convergcontrolt);
    
    fprintf (out,"%d %d ",tstorkm,tstorcm);
    
    ssle->print (out);

    fprintf (out,"%ld\n",diagcap);
    if (react>0)
      fprintf (out,"%ld\n",diagreact);
    
    fprintf (out,"\n");    
    
    break;
  }
    
    // *******************************************
    //  NONSTATIONARY PROBLEM WITH DISCONTINUITY
    // *******************************************
  case discont_nonstat_problem:
  case discont_nonlin_nonstat_problem:{
    
    //  time controller
    timecont.print (out);

    
    fprintf (out,"%d\n", tprt);

    hdbcont.print(out);

    fprintf (out,"%le %ld\n",alpha,nii);
    for (i=0;i<ntm;i++){
      fprintf (out,"%le %le\n",errarr[i],threshrhs[i]);
    }

    fprintf (out,"%d %d\n",tstorkm,tstorcm);
    
    //  data about solver of system of linear equations
    ssle->print (out);
  
    // diagonalization of capacity matrix
    fprintf (out,"%ld\n",diagcap);
    if (react>0)
      fprintf (out,"%ld\n",diagreact);
    
    break;
  }

  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }  
}



long probdesct::give_var_dofid(namevart var)
{
  return var_dofid[int(var)-1];
}
