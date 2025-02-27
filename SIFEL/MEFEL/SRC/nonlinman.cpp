#include "nonlinman.h"
#include <string.h>
#include <stdlib.h>


nonlinman::nonlinman (void)
{

  //  type of solver of nonlinear algebraic equation system
  tnlinsol = (nonlinsolvertype) 0;
  //  norm measure of displacement increments
  displnorm = (displacementnorm) 0;

  nial = 0;
  niilal = 0;
  erral = 0.0;
  dlal = 0.0;
  dlmaxal = 0.0;
  psial = 0.0;
  nsnal=0;
  selnodal=NULL;
  hdbackupal=0;
  nsdofal=0;
  seldofal=NULL;
  nmstrcomp=0;
  mstrastre=NULL;
  mstrid=NULL;
  //  determination of increment of load parameter
  dlam=nodetermination;

  check_tot_al=off;

  ninr = 0;
  niilnr = 0;
  errnr = 0.0;
  incrnr = 0.0;
  minincrnr = 0.0;
  rnormtnr = rel_load_norm;

  nienr = 0;
  hdbr = 0;
  hdbid = 0;
  backupfname[0] = 0;

  //  type of stiffness matrix
  stmat = initial_stiff;
  //  number of steps with constant stiffness %matrix
  ithstiffchange=1;
  jthstiffchange=1;

  check_div = off;
  div_min_steps = 3;
  divc_step = NULL;
  divc_err  = NULL;
  check_lambdar = off;
  lambdar = 0.0;
  errl = 0.0;

  // dissipation increment
  tau_ini = tau_lim = tau = 0.0;
}



nonlinman::~nonlinman (void)
{
  delete [] selnodal;
  delete [] seldofal;
  delete [] mstrastre;
  delete [] mstrid;
  delete [] divc_step;
  delete [] divc_err;
}



/**
  The function reads parameters for individual types of solvers.
  
  @param in - pointer to the opened text file
  @param mespr - flag for message printing control (0 = no messages are printed, 1 = messages are sentto stdout)

  @return The function does not return anything.

  Created by Jaroslav Kruis, 2008
  Modified by Tomas Koudelka, 08.2011
*/
void nonlinman::read (XFILE *in,long mespr)
{
  //  type of nonlinear solver
  //  1- arc-length
  //  2 - newton method
  //  3 - arclrv1
  //  4 - newtonrv1
  //  5 - newtonrestart
  //  6 - displctrl
  //  7 - displctrlrv
  //  8 - arclrv
  //  10 - adaptram
  //  20 - arclg
  //  30 - newtong
  xfscanf (in, "%k%m", "type_of_nonlin_solver", &nonlinsolvertype_kwdset, (int*)&tnlinsol);

  //  type of stiffness matrix (initial, tangent)
  //  1 - initial_stiff
  //  2 - tangent_stiff
  //  3 - secant_stiff
  //  5 - incr_tangent_stiff
  //  6 - ijth_tangent_stiff
  xfscanf (in,"%k%m","stiffmat_type", &stiffmatrix_kwdset, (int*)&stmat);
  switch (stmat){
  case initial_stiff:
  case tangent_stiff:
  case secant_stiff:
  case incr_tangent_stiff:{
    break;
  }
  case ijth_tangent_stiff:{
    xfscanf (in,"%ld %ld",&ithstiffchange,&jthstiffchange);
    break;
  }
  default:{
    print_err("unknown type of stiffness matrix is required",__FILE__,__LINE__,__func__);
  }
  }

  switch (tnlinsol){
  case arcl:
  case adaptram:
  case arclrv1:
  case arclrv:
  case arclg:{
    if (mespr==1){
      fprintf (stdout,"\n system of nonlinear equations will be solved by arc-length method (variant %s)",
               nonlinsolvertype_kwdset.get_str(tnlinsol));
    }
    read_arclength(in, tnlinsol, mespr);
    break;
  }
  case newton:
  case newtonrv1:
  case displctrl:
  case displctrlrv:
  case newtong:
    if (mespr==1){
      fprintf (stdout,"\n system of nonlinear equations will be solved by Newton-Raphson method (variant %s)",
               nonlinsolvertype_kwdset.get_str(tnlinsol));
    }
    read_newtonraphson(in, tnlinsol, mespr);
    break;
  /*
  case newtonrestart:{
    if (mespr==1)  fprintf (stdout,"\n system of nonlinear equations will be solved by Newton-Raphson method with ");
    xfscanf (in,"%ld %ld %lf %lf %lf",&ninr,&niilnr,&errnr,&incrnr,&minincrnr);
    xfscanf (in, "%ld %ld", &nienr, &hdbr);
    if (hdbr)
      xfscanf(in, "%s %ld", backupfname, &hdbid);
    if ((mespr==1) && (hdbr))  fprintf(stdout, "\n restart will be performed from file: %s\n record number : %ld\n", backupfname, hdbid);
    break;
  }
  */
  case dissip_incr:
    if (mespr==1){
      fprintf (stdout,"\n system of nonlinear equations will be solved by dissipation increment method");
    }
    read_dissipincr(in, mespr);
    break;
  default:
    print_err("unknown solver of nonlinear system of equations is required", __FILE__, __LINE__, __func__);
  }
  if (check_div == on)
  {
    divc_step = new double[div_min_steps];
    divc_err  = new double[div_min_steps];
  }
}



/**
  The function prints parameters for individual types of solvers to the opened text file.
  
  @param out - pointer to the opened text file for output

  @return The function does not return anything.

  Created by Jaroslav Kruis, 2008
  Modified by Tomas Koudelka, 08.2011
*/
void nonlinman::print (FILE *out)
{
  //  type of nonlinear solver
  fprintf (out,"%d\n",(int)tnlinsol);
  
  fprintf (out, "\n%d \n", (int)stmat);
  switch (stmat){
  case initial_stiff:
  case tangent_stiff:
  case secant_stiff:
  case incr_tangent_stiff:{
    break;
  }
  case ijth_tangent_stiff:{
    fprintf (out,"%ld %ld\n",ithstiffchange,jthstiffchange);
    break;
  }
  default:{
    print_err("unknown type of stiffness matrix is required",__FILE__,__LINE__,__func__);
  }
  }

  switch (tnlinsol){
  case arcl:
  case arclrv1:
  case arclrv:
  case arclg:
    print_arclength(out, tnlinsol);
    break;
  case newton:
  case newtonrv1:
  case displctrl:
  case displctrlrv:
  case newtong:
    print_newtonraphson(out, tnlinsol);
    break;
  /*
  case newtonrestart:{
    fprintf (out,"%ld %ld %le %le %le %d\n",ninr,niilnr,errnr,incrnr,minincrnr, rnormtnr);
    fprintf (out, "%ld %ld\n", nienr, hdbr);
    if (hdbr)
      fprintf(out, "%s\n %ld\n", backupfname, hdbid);
    break;
  }*/
  default:{}
  }
  
}



/**
  The function reads setup data for the arclength method.

  @param[in] in - pointer to the opened text file where the data will be read from
  @param[in] nlst - the given type of the nonlinear solver, i.e. variant of the arclength method
  @param[in] mespr - message printing flag

  @return The function does not return anything.

  Created by TKo, 03.2024
*/
void nonlinman::read_arclength(XFILE *in, nonlinsolvertype nlst, long mespr)
{
  xfscanf (in, "%k%m", "lambda_determination", &detlambda_kwdset, &dlam);
  xfscanf (in, "%k%ld%k%ld", "al_num_steps", &nial, "al_num_iter", &niilal);
  xfscanf (in, "%k%lf%k%lf", "al_error", &erral, "al_init_length", &dlal);
  xfscanf (in, "%k%lf%k%lf", "al_min_length", &dlminal, "al_max_length", &dlmaxal);      
  xfscanf (in, "%k%lf", "al_psi", &psial);
  
  if ((nlst == arclrv1) || (nlst == arclrv))
    xfscanf(in, "%k%lf", "reqval_lambda",&lambdar);
  if (nlst == arclrv)
    xfscanf(in, "%k%lf", "reqerr_lambda",&errl);
  
  xfscanf (in, "%k%m", "resid_norm_type", &resnormt_kwdset, (int*)&rnormtnr);
    
  if (nlst == arclg){
    // checking of divergency of inner iteration loop
    xfscanf(in, "%k%m", "check_div", &flagsw_kwdset, &check_div);
    if (check_div == on){
      xfscanf(in, "%k%ld", "div_min_steps", &div_min_steps);
      if (div_min_steps < 3){
        print_warning("the minimum number of performed steps in the divergency check is < 3.\n"
                      "It will be set to 3 automatically", __FILE__, __LINE__, __func__);
        div_min_steps = 3;
      }
    }
    
    // checking of maximum value of the total (cumulative) arc length
    xfscanf(in, "%k%m", "check_total_arcl", &flagsw_kwdset, &check_tot_al);
    if (check_tot_al == on)
      xfscanf(in, "%k%lf", "max_total_arc_length", &max_tot_al);
    // checking of  required value of lambda coefficent for proportional load vector
    xfscanf(in, "%k%m", "check_req_lambda", &flagsw_kwdset, &check_lambdar);
    if (check_lambdar == on){
      xfscanf(in, "%k%lf", "reqval_lambda",&lambdar);
      // tolerance for lambda tresholds
      xfscanf(in, "%k%lf", "reqerr_lambda",&errl);
    }
  }
  //fscanff (in,"%ld %d\n",hdbackupal,displnorm);
  //  1 - alldofs
  //  2 - seldofs
  //  3 - seldofscoord
  //  6 - selecnodes
  //  8 - nodesdistincr
  xfscanf (in,"%k%m","al_displ_contr_type", &displacementnorm_kwdset, (int*)&displnorm);
  read_displnorm(in, displnorm);
}



/**
  The function reads setup data for the Newton-Raphson method.

  @param[in] in - pointer to the opened text file where the data will be read from
  @param[in] nlst - the given type of the nonlinear solver, i.e. variant of the Newton-Raphson method
  @param[in] mespr - message printing flag

  @return The function does not return anything.

  Created by TKo, 03.2024
*/
void nonlinman::read_newtonraphson(XFILE *in, nonlinsolvertype nlst, long mespr)
{
  xfscanf (in, "%k%ld%k%ld", "nr_num_steps", &ninr, "nr_num_iter", &niilnr);
  xfscanf (in, "%k%lf%k%lf", "nr_error", &errnr, "nr_init_incr", &incrnr);
  xfscanf (in, "%k%lf%k%lf", "nr_minincr", &minincrnr, "nr_maxincr", &maxincrnr);
  xfscanf (in, "%k%m", "resid_norm_type", &resnormt_kwdset, (int*)&rnormtnr);
  if (nlst == newtonrv1){
    xfscanf (in, "%k%lf", "req_lambda", &lambdar);
    if (mespr==1)  fprintf(stdout, "\n requried value of lambda=% e", lambdar);
  }
  if (nlst == displctrlrv){
    xfscanf (in, "%k%lf%k%lf", "req_lambda", &lambdar, "req_lam_err", &errl);
    if (mespr==1)  fprintf(stdout, "\n requried value of lambda=% e, required error of lambda=% e", lambdar, errl);    
  }
  if (nlst == newtong){
    // checking of divergency of inner iteration loop
    xfscanf(in, "%k%m", "check_div", &flagsw_kwdset, &check_div);
    if (check_div == on)
    {
      xfscanf(in, "%k%ld", "div_min_steps", &div_min_steps);
      if (div_min_steps < 3)
      {
        print_warning("the minimum number of performed steps in the divergency check is < 3.\n"
                      "It will be set to 3 automatically", __FILE__, __LINE__, __func__);
        div_min_steps = 3;
      }
    }
  
    // checking of  required value of lambda coefficent for proportional load vector
    xfscanf(in, "%k%m", "check_req_lambda", &flagsw_kwdset, &check_lambdar);
    if (check_lambdar == on)
    {
      xfscanf(in, "%k%lf", "reqval_lambda",&lambdar);
      // tolerance for lambda tresholds
      xfscanf(in, "%k%lf", "reqerr_lambda",&errl);
    }
  }
}


    
/**
  The function reads setup data for the nonlinear solver with the constrained increment of dissipation.

  @param[in] in - pointer to the opened text file where the data will be read from
  @param[in] nlst - the given type of the nonlinear solver, i.e. variant of the Newton-Raphson method
  @param[in] mespr - message printing flag

  @return The function does not return anything.

  Created by TKo, 03.2024
*/
void nonlinman::read_dissipincr(XFILE *in, long mespr)
{
  xfscanf(in, "%k%le", "tau_ini", &tau_ini);
  xfscanf(in, "%k%le", "tau_max", &tau_max);
  xfscanf(in, "%k%le", "tau_lim", &tau_lim);
  xfscanf(in, "%k", "alt_arclength_setup");
  xfscanf (in, "%k%m", "type_of_alt_nonlin_solver", &nonlinsolvertype_kwdset, (int*)&altnlst);
  if ((altnlst == arcl) || (altnlst == arclrv1) || (altnlst == arclrv) || (altnlst == arclg))
    read_arclength(in, altnlst, mespr);
  else{
    print_err("wrong type of alternative nonlinear solver (%d) is required",
              __FILE__, __LINE__, __func__, int(altnlst));
    abort();
  }
}



/**
  The function reads setup for particular types of displacement norm computation.

  @param[in] in - pointer to the opened XFILE which the data will be read from
  @param[in] dn - type of required displacement norm

  @return The function does not return anything.

  Created by TKo, 03.2024
*/
void nonlinman::read_displnorm(XFILE *in, displacementnorm dn)
{
  switch (dn){
    case alldofs:{  break; }
    case seldofs:{
      read_seldof (in);
      break;
    }
    case seldofscoord:{
      read_seldofcoord (in);
      break;
    }
    case selmstr:{
      read_selmstr(in);
      break;
    }
    case selecnodes:{
      read_selnod (in);
      break;
    }
    case nodesdistincr:{
      read_selnoddistincr(in);
      break;
    }
    default:
      print_err("unknown norm measurement of displacement increment", __FILE__, __LINE__, __func__);
  }
}



/**
  The function reads selected dofs for arclength control.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by Jaroslav Kruis, 2008
  Modified by TKo, 03.2024
*/
void nonlinman::read_seldof (XFILE *in)
{
  long i;
  xfscanf (in,"%k%ld","num_sel_dofs",&nsdofal);
  selnodal = new long [nsdofal];
  seldofal = new long [nsdofal];
  for (i=0;i<nsdofal;i++){
    xfscanf (in,"%ld %ld",selnodal+i,seldofal+i);
    selnodal[i]--;  seldofal[i]--;
  }
}



/**
  The function reads coordinates of selected dofs for arclength control.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by Jaroslav Kruis, 2008
*/
void nonlinman::read_seldofcoord (XFILE *in)
{
  long i;
  xfscanf (in,"%k%ld", "num_sel_dof_coords", &nsdofal);
  selnodal = new long [nsdofal];
  seldofal = new long [nsdofal];
  selnodcoord = new double [3*nsdofal];
  for (i=0; i<nsdofal; i++){
    xfscanf (in,"%lf %lf %lf %ld",selnodcoord+3*i,selnodcoord+3*i+1,selnodcoord+3*i+2,seldofal+i);
    seldofal[i]--;
  }
}



/**
  The function reads selected macro-strain or macro-stress components for arclength control.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by Tomas Koudelka, 02.2020
*/
void nonlinman::read_selmstr(XFILE *in)
{
  long i;
  
  xfscanf(in, "%k%ld", "numcomp", &nmstrcomp);
  mstrastre = new strastre[nmstrcomp];
  mstrid = new long[nmstrcomp];
  for (i=0;i<nmstrcomp;i++){
    xfscanf (in,"%k%m%k%ld", "mstrcomptype", &strastre_kwdset, mstrastre+i, "mstrid", mstrid+i);
    mstrid[i]--;
  }  
}


/**
  The function reads selected nodes for arclength control.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by Jaroslav Kruis, 2008
*/
void nonlinman::read_selnod (XFILE *in)
{
  long i;
  xfscanf (in,"%k%ld", "num_sel_nodes", &nsnal);
  selnodal = new long [nsnal];
  for (i=0; i<nsnal; i++){
    xfscanf (in,"%ld",selnodal+i);
    selnodal[i]--;
  }
}



/**
  The function reads two selected nodes for arclength control by the distance increment of two nodes.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by TKo, 03.2024
*/
void nonlinman::read_selnoddistincr (XFILE *in)
{
  long i;
  
  xfscanf (in,"%k%ld","probdimal",&probdimal);
  nsnal=2;
  selnodal = new long [nsnal];
  for (i=0; i<nsnal; i++){
    xfscanf (in,"%ld",selnodal+i);
    selnodal[i]--;
  }
}



/**
  The function prints setup data for the arclength method to the text file.

  @param[in,out] out - pointer to the opened text file where the data will be printed to
  @param[in] nlst - the given type of the nonlinear solver, i.e. variant of the arclength method

  @return The function does not return anything.

  Created by TKo, 03.2024
*/
void nonlinman::print_arclength(FILE *out, nonlinsolvertype nlst)
{
  fprintf (out,"%d\n",(int)dlam);
  fprintf (out, "%ld %ld %e %e %e %e %e\n", nial, niilal, erral, dlal,
           dlminal, dlmaxal, psial);
  if ((nlst == arclrv1) || (nlst == arclrv))
    fprintf(out, "%f\n", lambdar);
  if (tnlinsol == arclrv)
    fprintf(out, "%f\n", errl);
  
  fprintf (out, "%d\n", (int)rnormtnr);
  if (nlst == arclg){
    // checking of divergency of inner iteration loop
    fprintf(out, "%d", check_div);
    if (check_div == on)
      fprintf(out, " %ld\n", div_min_steps+1);
    else
      fprintf(out, "\n");
    
    // checking of maximum value of the total (cumulative) arc length
    fprintf(out, "%d", check_tot_al);
    if (check_tot_al == on)
      fprintf(out, " %le\n", max_tot_al);
    else
      fprintf(out, "\n");
    
    // checking of  required value of lambda coefficent for proportional load vector
    fprintf(out, "%d", check_lambdar);
    if (check_lambdar == on)
      fprintf(out, " %le %le\n", lambdar, errl);
    else
      fprintf(out, "\n");
  }
  //fprintf (out,"%ld %d\n",hdbackupal,displnorm);
  fprintf (out,"%d\n",displnorm);
  print_displnorm(out, displnorm);
}



/**
  The function reads setup data for the Newton-Raphson method.

  @param[in,out] out - pointer to the opened text file where the data will be printed to
  @param[in] nlst - the given type of the nonlinear solver, i.e. variant of the Newton-Raphson method

  @return The function does not return anything.

  Created by TKo, 03.2024
*/
void nonlinman::print_newtonraphson(FILE *out, nonlinsolvertype nlst)
{
  fprintf (out, "%ld %ld %le %le %le %le %d", ninr, niilnr, errnr, incrnr, minincrnr, maxincrnr, rnormtnr);
  if (nlst == newtonrv1){
    fprintf (out, " %le", lambdar);
  }
  if (nlst == displctrlrv){
    fprintf (out, " %le %le", lambdar, errl);
  }
  fprintf(out, "\n");
  /*
  case newtonrestart:{
    fprintf (out,"%ld %ld %le %le %le %d\n",ninr,niilnr,errnr,incrnr,minincrnr, rnormtnr);
    fprintf (out, "%ld %ld\n", nienr, hdbr);
    if (hdbr)
      fprintf(out, "%s\n %ld\n", backupfname, hdbid);
    break;
  }*/
  if (nlst == newtong){
    // checking of divergency of inner iteration loop
    fprintf(out, "%d", check_div);
    if (check_div == on)
      fprintf(out, " %ld\n", div_min_steps+1);
    else
      fprintf(out, "\n");
    // checking of  required value of lambda coefficent for proportional load vector
    fprintf(out, "%d", check_lambdar);
    if (check_lambdar == on)
      fprintf(out, " %le %le\n", lambdar, errl);
    else
      fprintf(out, "\n");
  }
}



/**
  The function prints setup data for the nonlinear solver with the constrained increment of dissipation.

  @param[in, out] out - pointer to the opened text file where the data will be printed to

  @return The function does not return anything.

  Created by TKo, 03.2024
*/
void nonlinman::print_dissipincr(FILE *out)
{
  fprintf(out, "%le", tau_ini);
  fprintf(out, "%le", tau_max);
  fprintf(out, "%le", tau_lim);
  fprintf(out, "%d", int(altnlst));
  if ((altnlst == arcl) || (altnlst == arclrv1) || (altnlst == arclrv) || (altnlst == arclg))
    print_arclength(out, altnlst);
  else{
    print_err("wrong type of alternative nonlinear solver (%d) is required",
              __FILE__, __LINE__, __func__, int(altnlst));
    abort();
  }
}



/**
  The function printss setup of the actual type of the displacement norm computation.

  @param[in,out] out - pointer to the opened FILE where the data will be printed to
  @param[in] dn - type of required displacement norm

  @return The function does not return anything.

  Created by TKo, 03.2024
*/
void nonlinman::print_displnorm(FILE *out, displacementnorm dn)
{
  switch (dn){
    case alldofs:{  break; }
    case seldofs:{
      print_seldof (out);
      break;
    }
    case seldofscoord:{
      print_seldofcoord (out);
      break;
    }
    case selmstr:{
      print_selmstr(out);
      break;
    }
    case selecnodes:{
      print_selnod (out);
      break;
    }
    case nodesdistincr:{
      print_selnoddistincr(out);
      break;
    }
    default:
      print_err("unknown norm measurement of displacement increment", __FILE__, __LINE__, __func__);
  }
}



/**
  The function prints selected dofs for arclength control.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by Jaroslav Kruis, 2008
*/
void nonlinman::print_seldof (FILE *out)
{
  long i;
  
  fprintf (out,"%ld\n",nsdofal);
  for (i=0;i<nsdofal;i++){
    fprintf (out,"%ld %ld\n",selnodal[i]+1,seldofal[i]+1);
  }
  
}



/**
  The function prints coordinates of selected dofs for arclength control.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by Jaroslav Kruis, 2008
*/
void nonlinman::print_seldofcoord (FILE *out)
{
  long i;
  
  fprintf (out,"%ld\n",nsdofal);
  for (i=0;i<nsdofal;i++){
    fprintf (out,"%f %f %f %ld\n",selnodcoord[3*i],selnodcoord[3*i+1],selnodcoord[3*i+2],seldofal[i]+1);
  }
  
}



/**
  The function prints selected macro-strain or macro-stress components for arclength control.

  @param out - pointer to the opened text file

  @return The function does not return anything.

  Created by Tomas Koudelka, 02.2020
*/
void nonlinman::print_selmstr(FILE *out)
{
  long i;
  
  fprintf(out, "%ld\n", nmstrcomp);
  for (i=0;i<nmstrcomp;i++){    
    fprintf(out,"%d%ld", mstrastre[i], mstrid[i]+1);
  }  
}



/**
  The function prints selected nodes for arclength control.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by Jaroslav Kruis, 2008
*/
void nonlinman::print_selnod (FILE *out)
{
  long i;
  
  fprintf (out,"%ld\n",nsnal);
  for (i=0;i<nsnal;i++){
    fprintf (out,"%ld\n",selnodal[i]+1);
  }
  
}



/**
  The function prints selected nodes for arclength control.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by Jaroslav Kruis, 2008
*/
void nonlinman::print_selnoddistincr (FILE *out)
{
  fprintf (out,"%ld\n",probdimal);
  print_selnod(out);  
}



/**
  The function initializes data from the another nonlinman structure.

  @param nm - nonlinman structure used fro initialization

  @raturn The function does not return anything.

  Created by Jaroslav Kruis, 2008
*/
void nonlinman::initiate (nonlinman &nm)
{
  long i;

  tnlinsol = nm.tnlinsol;

  stmat = nm.stmat;
  ithstiffchange = nm.ithstiffchange;
  jthstiffchange = nm.jthstiffchange;

  //  ARC-LENGTH METHOD
  displnorm = nm.displnorm;
  dlam = nm.dlam;
  hdbackupal = nm.hdbackupal;
  nial = nm.nial;
  niilal = nm.niilal;
  erral = nm.erral;
  dlal = nm.dlal;
  dlminal = nm.dlminal;
  dlmaxal = nm.dlmaxal;
  psial = nm.psial;
  
  nsnal = nm.nsnal;
  nsdofal = nm.nsdofal;
  probdimal = nm.probdimal;
  check_tot_al = nm.check_tot_al;
  max_tot_al = nm.max_tot_al;

  nxal = nm.nxal;
  nyal = nm.nyal;
  nzal = nm.nzal;
  
  
  switch (displnorm){
  case alldofs:{  break; }
  case seldofs:{
    selnodal = new long [nsdofal];
    seldofal = new long [nsdofal];
    
    for (i=0;i<nsdofal;i++){
      selnodal[i]=nm.selnodal[i];
      seldofal[i]=nm.seldofal[i];
    }
    
    break;
  }
  case seldofscoord:{
    seldofal = new long [nsdofal];
    selnodcoord = new double [3*nsdofal];
    
    for (i=0;i<nsdofal;i++){
      selnodcoord[3*i]=nm.selnodcoord[3*i];
      selnodcoord[3*i+1]=nm.selnodcoord[3*i+1];
      selnodcoord[3*i+2]=nm.selnodcoord[3*i+2];
      seldofal[i]=nm.seldofal[i];
    }
    
    break;
  }
  case selecnodes:{
    selnodal = new long [nsnal];
    
    for (i=0;i<nsnal;i++){
      selnodal[i]=nm.selnodal[i];
    }
    
    break;
  }
  case nodesdistincr:{
    nsnal=2;
    selnodal = new long [nsnal];
    
    for (i=0;i<nsnal;i++){
      selnodal[i]=nm.selnodal[i];
    }
    
    break;
  }
  default:
    print_err("unknown norm measurement of displacement increment", __FILE__, __LINE__, __func__);
  }
  
  // COMMON SETUP for both Newton-Raphson and arclength methods
  check_lambdar = nm.check_lambdar;
  lambdar = nm.check_lambdar;
  errl = nm.errl;

  check_div = nm.check_div;
  div_min_steps = nm.div_min_steps;
  if (check_div == on)
  {
    divc_step = new double[div_min_steps];
    divc_err = new double[div_min_steps]; 
    for(i=0; i<div_min_steps; i++)
    {
      divc_step[i] = nm.divc_step[i];
      divc_err[i] = nm.divc_err[i];
    }
  }
  
  //  NEWTON-RAPHSON METHOD
  ninr = nm.ninr;
  niilnr = nm.niilnr;
  errnr = nm.errnr;
  incrnr = nm.incrnr;
  minincrnr = nm.minincrnr;
  maxincrnr = nm.maxincrnr;
  lambdar = nm.lambdar;
  nienr = nm.nienr;
  hdbr  = nm.hdbr;
  if (hdbr)
  {
    hdbid = nm.hdbid;
    for (i=0; i < long(strlen(nm.backupfname)); i++)
      backupfname[i] = nm.backupfname[i];
  }
}
