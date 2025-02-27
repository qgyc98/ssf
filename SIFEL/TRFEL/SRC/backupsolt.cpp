#include <stdio.h>
#include <stdlib.h>
#include "globalt.h"
#include "backupsolt.h"


/**
  Function saves all necessary data to backup file in the given time step. The process is driven by 
  backup controler Tp->hdbcont. The data can be used for restarting analysis from the given step.
  
  @param r - %vector of unknowns
  @param dr - %vector of increments of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - time controller used at solver
  @param n    - number of unknowns (size of %vectors r and fp)

  9.2008 TKo
  12.2008 TKr
*/
void solvert_save (double *r,double *dr,double *fp,long ni,double time,double dt,timecontr &tc,long n)
{
  switch (Tp->hdbcont.hdbfmts)
  {
    case text:
      solvert_save_text(r, dr, fp, ni, time, dt, tc, n);
      break;
    case binary:
      solvert_save_binary(r, dr, fp, ni, time, dt, tc, n);
      break;
    default: 
      print_err("unknown type of backup file format is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function saves all necessary data to the text backup file in the given time step. The process is driven by 
  backup controler Tp->hdbcont. The data can be used for restarting analysis from the given step.
  
  @param r - %vector of unknowns
  @param dr - %vector of increments of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - time controller used at solver
  @param n    - number of unknowns (size of %vectors r and fp)

  9.2008 TKo  
  12.2008 TKr
*/
void solvert_save_text (double *r,double *dr,double *fp,long ni,double time,double dt,timecontr &tc,long n)
{
  switch (Tp->hdbcont.hdbtype)
  {
    case hdbs_single:
    case hdbrs_single:
      solvert_save_text_single(r, dr, fp, ni, time, dt, tc, n);
      break;
    case hdbs_multiple:
    case hdbrs_multiple:
      solvert_save_text_multiple(r, dr, fp, ni, time, dt, tc, n);
      break;
    default:
      print_err("unknown type of backup is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function saves all necessary data to the binary backup file in the given time step. The process is driven by 
  backup controler Tp->hdbcont. The data can be used for restarting analysis from the given step.
  
  @param r - %vector of unknowns
  @param dr - %vector of increments of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - time controller used at solver
  @param n    - number of unknowns (size of %vectors r and fp)

  9.2008 TKo  
  12.2008 TKr
*/
void solvert_save_binary (double *r,double *dr,double *fp,long ni,double time,double dt,timecontr &tc,long n)
{
  switch (Tp->hdbcont.hdbtype)
  {
    case hdbs_single:
    case hdbrs_single:
      solvert_save_binary_single(r, dr, fp, ni, time, dt, tc, n);
      break;
    case hdbs_multiple:
    case hdbrs_multiple:
      solvert_save_binary_multiple(r, dr, fp, ni, time, dt, tc, n);
      break;
    default:
      print_err("unknown type of backup is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function saves all necessary data to the single text backup file in the given time step. The process is driven by 
  backup controler Tp->hdbcont. The data can be used for restarting analysis from the given step.
  
  @param r - %vector of unknowns
  @param dr - %vector of increments of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - time controller used at solver
  @param n    - number of unknowns (size of %vectors r and fp)

  9.2008 TKo  
  12.2008 TKr
*/
void solvert_save_text_single (double *r,double *dr,double *fp,long ni,double time,double dt,timecontr &tc,long n)
{
  long i;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];
  int prec = (int)Tp->hdbcont.prec;

  sprintf(name, "%s.%ld.bac",Tp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  fprintf (aux,"\n\n\n\n");

  // storage of time controller data
  tc.save_txt(aux, prec);
  //  attained time of solver
  fprintf (aux,"%.*le\n",prec,time);
  //  attained time increment of solver
  fprintf (aux,"%.*le\n",prec,dt);
  //  attained number of steps
  fprintf (aux,"%ld\n",ni);
  //  number of degrees of freedom
  fprintf (aux,"%ld\n",n);

  //  left hand side vector
  for (i=0;i<n;i++){
    fprintf (aux,"%.*le\n",prec,r[i]);
  }
  
  fprintf (aux,"\n");
  
  //  increment of left hand side vector
  for (i=0;i<n;i++){
    fprintf (aux,"%.*le\n",prec,dr[i]);
  }
  
  fprintf (aux,"\n");

  //  vector of the right hand side from previous step
  for (i=0;i<n;i++){
    fprintf (aux,"%.*le\n",prec,fp[i]);
  }
  
  //  number of integration points in the problem
  fprintf (aux,"\n\n %ld\n",Tm->tnip);

  //  data on integration points
  Tm->save_intpointst_txt (aux, Tp->hdbcont.selelems, Tp->hdbcont.selother_s);
  
  fclose (aux);

  if (Tp->hdbcont.rmold)
  {
    if (Tp->hdbcont.rmold_id>=0)
    {
      sprintf(name, "%s.%ld.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
      remove(name);
    }
    Tp->hdbcont.rmold_id = ni;
  }
}



/**
  Function saves all necessary data to the several text backup files in the given time step.
  Individual quantities are saved into separated files, each file contains given quantity for all
  integration points or all data needed by solver. The file names consist of name specified 
  in the Tp->hdbcont, time step id, description of the saved quantity and .bac suffix. Description
  of the quantity can be :
  solver - file contains date needed by solver
  strain - array of strains at all integration points
  stress - array of stresses at all integration points
  other  - selected components of eqother array at all integration points
  nonloc - adjacent integration points and their distances at all integration points (used only for nonlocal material models)
           These points are saved only once in the adjacip() function.
  The process is driven by the backup controler Tp->hdbcont. The data can be used for restarting 
  analysis from the given step.
  
  @param r - %vector of unknowns
  @param dr - %vector of increments of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - time controller used at solver
  @param n    - number of unknowns (size of %vectors r and fp)

  9.2008 TKo  
  12.2008 TKr
*/
void solvert_save_text_multiple (double *r,double *dr,double *fp,long ni,double time,double dt,timecontr &tc,long n)
{
  long i;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];
  int prec = (int)Tp->hdbcont.prec;

  sprintf(name, "%s.%ld.solver.bac",Tp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }

  // storage of time controller data
  tc.save_txt(aux, prec);
  //  attained time of solver
  fprintf (aux,"%.*le\n",prec,time);
  //  attained time increment of solver
  fprintf (aux,"%.*le\n",prec,dt);
  //  attained number of steps
  fprintf (aux,"%ld\n",ni);
  //  number of degrees of freedom
  fprintf (aux,"%ld\n",n);

  //  left hand side vector
  for (i=0;i<n;i++){
    fprintf (aux,"%.*le\n",prec,r[i]);
  }
  
  fprintf (aux,"\n");
    
  //  increment of left hand side vector
  for (i=0;i<n;i++){
    fprintf (aux,"%.*le\n",prec,dr[i]);
  }
  
  fprintf (aux,"\n");

  //  vector of the right hand side from previous step
  for (i=0;i<n;i++){
    fprintf (aux,"%.*le\n",prec,fp[i]);
  }
  
  //  number of integration points in the problem
  fprintf (aux,"\n\n %ld\n",Tm->tnip);

  fclose (aux);

  //  data on integration points
  Tm->save_intpointst_txt (ni, Tp->hdbcont.selelems, Tp->hdbcont.selother_s);

  if (Tp->hdbcont.rmold)
  {
    if (Tp->hdbcont.rmold_id>=0)
    {
      sprintf(name, "%s.%ld.solver.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
      remove(name);
      sprintf(name, "%s.%ld.grad.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
      remove(name);
      sprintf(name, "%s.%ld.flux.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
      remove(name);
      sprintf(name, "%s.%ld.other.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
      remove(name);
      sprintf(name, "%s.%ld.mefquant.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
      remove(name);
      if (Tm->tnaip > 0)
      {
        sprintf(name, "%s.%ld.agrad.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
        remove(name);
        sprintf(name, "%s.%ld.aflux.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
        remove(name);
        sprintf(name, "%s.%ld.aother.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
        remove(name);
      }
    }
    Tp->hdbcont.rmold_id = ni;
  }  
}



/**
  Function saves all necessary data to the single binary backup file in the given time step. The process is driven by 
  backup controler Tp->hdbcont. The data can be used for restarting analysis from the given step.
  
  @param r - %vector of unknowns
  @param dr - %vector of increments of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - time controller used at solver
  @param n    - number of unknowns (size of %vectors r and fp)

  9.2008 TKo  
  12.2008 TKr
*/
void solvert_save_binary_single (double *r,double *dr,double *fp,long ni,double time,double dt,timecontr &tc,long n)
{
  FILE *aux;
  char name[FNAMELEN+30];
  char emsg[FNAMELEN+100];

  sprintf(name, "%s.%ld.bac",Tp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }

  // storage of time controller data
  tc.save_bin(aux);
  //  attained time of solver
  fwrite (&time, sizeof(time), 1, aux);
  //  attained time increment of solver
  fwrite (&dt, sizeof(dt), 1, aux);
  //  attained number of steps
  fwrite (&ni, sizeof(ni), 1, aux);
  //  number of degrees of freedom
  fwrite (&n, sizeof(n), 1, aux);

  //  left hand side vector
  fwrite(r, sizeof(*r), n, aux);
  
  //  increment of left hand side vector
  fwrite(dr, sizeof(*dr), n, aux);
  
  //  vector of the right hand side from previous step
  fwrite(fp, sizeof(*fp), n, aux);
  
  //  number of integration points in the problem
  fwrite (&Tm->tnip, sizeof(Tm->tnip), 1, aux);

  //  data on integration points
  Tm->save_intpointst_bin (aux, Tp->hdbcont.selelems, Tp->hdbcont.selother_s); 
  fclose (aux);

  if (Tp->hdbcont.rmold)
  {
    if (Tp->hdbcont.rmold_id>=0)
    {
      sprintf(name, "%s.%ld.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
      remove(name);
    }
    Tp->hdbcont.rmold_id = ni;
  }
}



/**
  Function saves all necessary data to the several binary backup files in the given time step.
  Individual quantities are saved into separated files, each file contains given quantity for all
  integration points or all data needed by solver. The file names consist of name specified 
  in the Tp->hdbcont, time step id, description of the saved quantity and .bac suffix. Description
  of the quantity can be :
  solver - file contains date needed by solver
  strain - array of strains at all integration points
  stress - array of stresses at all integration points
  other  - selected components of eqother array at all integration points
  nonloc - adjacent integration points and their distances at all integration points (used only for nonlocal material models)
           These points are saved only once in the adjacip() function.
  The process is driven by the backup controler Tp->hdbcont. The data can be used for restarting 
  analysis from the given step.
  
  @param r - %vector of unknowns
  @param dr - %vector of increments of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - time controller used at solver
  @param n    - number of unknowns (size of %vectors r and fp)

  9.2008 TKo  
  12.2008 TKr
*/
void solvert_save_binary_multiple (double *r,double *dr,double *fp,long ni,double time,double dt,timecontr &tc,long n)
{
  FILE *aux;
  char name[FNAMELEN+30];
  char emsg[FNAMELEN+100];

  sprintf(name, "%s.%ld.solver.bac",Tp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }

  // storage of time controller data
  tc.save_bin(aux);
  //  attained time of solver
  fwrite (&time, sizeof(time), 1, aux);
  //  attained time increment of solver
  fwrite (&dt, sizeof(dt), 1, aux);
  //  attained number of steps
  fwrite (&ni, sizeof(ni), 1, aux);
  //  number of degrees of freedom
  fwrite (&n, sizeof(n), 1, aux);

  //  left hand side vector
  fwrite(r, sizeof(*r), n, aux);

  //  increment of left hand side vector
  fwrite(dr, sizeof(*dr), n, aux);
  
  //  vector of the right hand side from previous step
  fwrite(fp, sizeof(*fp), n, aux);
  
  //  number of integration points in the problem
  fwrite (&Tm->tnip, sizeof(Tm->tnip), 1, aux);

  fclose (aux);

  //  data on integration points
  Tm->save_intpointst_bin (ni, Tp->hdbcont.selelems, Tp->hdbcont.selother_s); 

  if (Tp->hdbcont.rmold)
  {
    if (Tp->hdbcont.rmold_id>=0)
    {
      sprintf(name, "%s.%ld.solver.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
      remove(name);
      sprintf(name, "%s.%ld.grad.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
      remove(name);
      sprintf(name, "%s.%ld.flux.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
      remove(name);
      sprintf(name, "%s.%ld.other.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
      remove(name);
      sprintf(name, "%s.%ld.mefquant.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
      remove(name);
      if (Tm->tnaip > 0)
      {
        sprintf(name, "%s.%ld.agrad.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
        remove(name);
        sprintf(name, "%s.%ld.aflux.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
        remove(name);
        sprintf(name, "%s.%ld.aother.bac",Tp->hdbcont.hdbnames,Tp->hdbcont.rmold_id);
        remove(name);
      }
    }
    Tp->hdbcont.rmold_id = ni;
  }  
}



/**
  Function restores all necessary data from the backup file. The process is driven by 
  backup controler Tp->hdbcont. Data from solver are passed through parameters of the function. 
  
  @param r - %vector of unknowns
  @param dr - %vector of increments of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - time controller used at solver
  @param n    - number of unknowns (size of %vectors r and fp)

  9.2008 TKo  
  12.2008 TKr
*/
void solvert_restore (double *r,double *dr,double *fp,long &ni,double &time,double &dt,timecontr &tc,long &n)
{
  switch (Tp->hdbcont.hdbfmtr)
  {
    case text:
      solvert_restore_text(r, dr, fp, ni, time, dt, tc, n);
      break;
    case binary:
      solvert_restore_binary(r, dr, fp, ni, time, dt, tc, n);
      break;
    default: 
      print_err("unknown type of backup file format is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function restores all necessary data from the text backup file. The process is driven by 
  backup controler Tp->hdbcont. Data from solver are passed through parameters of the function. 
  
  @param r - %vector of unknowns
  @param dr - %vector of increments of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - time controller used at solver
  @param n    - number of unknowns (size of %vectors r and fp)

  9.2008 TKo  
  12.2008 TKr
*/
void solvert_restore_text (double *r,double *dr,double *fp,long &ni,double &time,double &dt,timecontr &tc,long &n)
{
  switch (Tp->hdbcont.hdbtype)
  {
    case hdbr_single:
    case hdbrs_single:
      solvert_restore_text_single(r, dr, fp, ni, time, dt, tc, n);
      break;
    case hdbr_multiple:
    case hdbrs_multiple:
      solvert_restore_text_multiple(r, dr, fp, ni, time, dt, tc, n);
      break;
    default:
      print_err("unknown type of backup is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function restores all necessary data from the binary backup file. The process is driven by 
  backup controler Tp->hdbcont. Data from solver are passed through parameters of the function. 
  
  @param r - %vector of unknowns
  @param dr - %vector of increments of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - time controller used at solver
  @param n    - number of unknowns (size of %vectors r and fp)

  9.2008 TKo  
  12.2008 TKr
*/
void solvert_restore_binary (double *r,double *dr,double *fp,long &ni,double &time,double &dt,timecontr &tc,long &n)
{
  switch (Tp->hdbcont.hdbtype)
  {
    case hdbr_single:
    case hdbrs_single:
      solvert_restore_binary_single(r, dr, fp, ni, time, dt, tc, n);
      break;
    case hdbr_multiple:
    case hdbrs_multiple:
      solvert_restore_binary_multiple(r, dr, fp, ni, time, dt, tc, n);
      break;
    default:
      print_err("unknown type of backup is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function restores all necessary data from the one huge text backup file. The process is driven by 
  backup controler Tp->hdbcont. Data from solver are passed through parameters of the function. 
  
  @param r - %vector of unknowns
  @param dr - %vector of increments of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - time controller used at solver
  @param n    - number of unknowns (size of %vectors r and fp)

  9.2008 TKo  
  12.2008 TKr
*/
void solvert_restore_text_single (double *r,double *dr,double *fp,long &ni,double &time,double &dt,timecontr &tc,long &n)
{
  long i;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];

  aux = fopen(Tp->hdbcont.hdbnamer,"rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }

  // restorage of time controller data  
  tc.restore_txt(aux);
  //  attained time of solver
  fscanf (aux,"%le",&time);
  //  attained time increment of solver
  fscanf (aux,"%le",&dt);
  //  attained number of steps
  fscanf (aux,"%ld\n",&ni);
  //  number of degrees of freedom
  fscanf (aux,"%ld\n",&n);

  if (n!=Ndoft){
    print_err("number of DOFs in backup file and in the actual problem is not same\n", __FILE__, __LINE__, __func__);
    abort();
  }

  //  left hand side vector
  for (i=0;i<n;i++){
    fscanf (aux,"%le",r+i);
  }
  //  increment of left hand side vector
  for (i=0;i<n;i++){
    fscanf (aux,"%le",dr+i);
  }
  //  vector of the right hand side from previous step
  for (i=0;i<n;i++){
    fscanf (aux,"%le",fp+i);
  }
  
  //  number of integration points in the problem
  fscanf (aux,"%ld",&i);
  
  if (i!=Tm->tnip){
    print_err("number of integration points in backup file and in the actual problem is not same\n", __FILE__, __LINE__, __func__);
    abort();
  }

  //  data on integration points
  Tm->restore_intpointst_txt (aux, Tp->hdbcont.selelemr, Tp->hdbcont.selother_r, Tp->hdbcont.selother_id);

  fclose(aux);
}



/**
  Function restores all necessary data from the separated text backup file. The process is driven by 
  backup controler Tp->hdbcont. Data from solver are passed through parameters of the function. 
  Individual quantities have to be saved into separated files, each file contains given quantity for all
  integration points or all data needed by solver. The file names consist of name specified 
  in the Tp->hdbcont, time step id, description of the saved quantity and .bac suffix. Description
  of the quantity can be :
  solver - file contains date needed by solver
  grad   - array of gradients at all integration points
  fluxes - array of fluxes at all integration points
  other  - selected components of eqother array at all integration points
           These points are saved only once in the adjacip() function.
  The process is driven by the backup controler Tp->hdbcont. The data can be used for restarting 
  analysis from the given step.
  
  @param r - %vector of unknowns
  @param dr - %vector of increments of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - time controller used at solver
  @param n    - number of unknowns (size of %vectors r and fp)

  9.2008 TKo  
  12.2008 TKr
*/
void solvert_restore_text_multiple (double *r,double *dr,double *fp,long &ni,double &time,double &dt,timecontr &tc,long &n)
{
  long i;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];

  sprintf(name, "%s.solver.bac",Tp->hdbcont.hdbnamer);
  aux = fopen(name,"rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  
  // restorage of time controller data  
  tc.restore_txt(aux);
  //  attained time of solver
  fscanf (aux,"%le",&time);
  //  attained time increment of solver
  fscanf (aux,"%le",&dt);
  //  attained number of steps
  fscanf (aux,"%ld",&ni);
  //  number of degrees of freedom
  fscanf (aux,"%ld",&n);

  if (n!=Ndoft){
    print_err("number of DOFs in backup file and in the actual problem is not same\n", __FILE__, __LINE__, __func__);
    abort();
  }

  //  left hand side vector
  for (i=0;i<n;i++){
    fscanf (aux,"%le",r+i);
  }

  //  increment of left hand side vector
  for (i=0;i<n;i++){
    fscanf (aux,"%le",dr+i);
  }
  //  vector of the right hand side from previous step
  for (i=0;i<n;i++){
    fscanf (aux,"%le",fp+i);
  }
  
  //  number of integration points in the problem
  fscanf (aux,"%ld",&i);
  
  if (i!=Tm->tnip){
    print_err("number of integration points in backup file and in the actual problem is not same\n", __FILE__, __LINE__, __func__);
    abort();
  }

  fclose(aux);

  //  data on integration points
  Tm->restore_intpointst_txt (Tp->hdbcont.selelemr, Tp->hdbcont.selother_r, Tp->hdbcont.selother_id);
}



/**
  Function restores all necessary data from the one huge binary backup file. The process is driven by 
  backup controler Tp->hdbcont. Data from solver are passed through parameters of the function. 
  
  @param r - %vector of unknowns
  @param dr - %vector of increments of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - time controller used at solver
  @param n    - number of unknowns (size of %vectors r and fp)

  9.2008 TKo  
  12.2008 TKr
*/
void solvert_restore_binary_single (double *r,double *dr,double *fp,long &ni,double &time,double &dt,timecontr &tc,long &n)
{
  FILE *aux;
  char name[FNAMELEN+30];
  char emsg[FNAMELEN+100];

  aux = fopen(Tp->hdbcont.hdbnamer,"rb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }

  // restorage of time controller data
  tc.restore_bin(aux);
  //  attained time of solver
  fread (&time, sizeof(time), 1, aux);
  //  attained time increment of solver
  fread (&dt, sizeof(dt), 1, aux);
  //  attained number of steps
  fread (&ni, sizeof(ni), 1, aux);
  //  number of degrees of freedom
  fread (&n, sizeof(n), 1, aux);

  //  left hand side vector
  fread(r, sizeof(*r), n, aux);

  //  increment of left hand side vector
  fread(dr, sizeof(*dr), n, aux);
  
  //  vector of the right hand side from previous step
  fread(fp, sizeof(*fp), n, aux);
  
  //  number of integration points in the problem
  fread (&Tm->tnip, sizeof(Tm->tnip), 1, aux);

  //  data on integration points
  Tm->restore_intpointst_bin (aux, Tp->hdbcont.selelemr, Tp->hdbcont.selother_r, Tp->hdbcont.selother_id);
  fclose (aux);
}



/**
  Function restores all necessary data from the separated text backup file. The process is driven by 
  backup controler Tp->hdbcont. Data from solver are passed through parameters of the function. 
  Individual quantities have to be saved into separated files, each file contains given quantity for all
  integration points or all data needed by solver. The file names consist of name specified 
  in the Tp->hdbcont, time step id, description of the saved quantity and .bac suffix. Description
  of the quantity can be :
  solver - file contains date needed by solver
  grad   - array of gradients at all integration points
  fluxes - array of fluxes at all integration points
  other  - selected components of eqother array at all integration points

           These points are saved only once in the adjacip() function.
  The process is driven by the backup controler Tp->hdbcont. The data can be used for restarting 
  analysis from the given step.
  
  @param r - %vector of unknowns
  @param dr - %vector of increments of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - time controller used at solver
  @param n    - number of unknowns (size of %vectors r and fp)

  9.2008 TKo  
  12.2008 TKr
*/
void solvert_restore_binary_multiple (double *r,double *dr,double *fp,long &ni,double &time,double &dt,timecontr &tc,long &n)
{
  FILE *aux;
  char name[FNAMELEN+30];
  char emsg[FNAMELEN+100];

  sprintf(name, "%s.solver.bac",Tp->hdbcont.hdbnamer);
  aux = fopen(name,"rb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }

  // restorage of time controller data
  tc.restore_bin(aux);
  //  attained time of solver
  fread (&time, sizeof(time), 1, aux);
  //  attained time increment of solver
  fread (&dt, sizeof(dt), 1, aux);
  //  attained number of steps
  fread (&ni, sizeof(ni), 1, aux);
  //  number of degrees of freedom
  fread (&n, sizeof(n), 1, aux);

  //  left hand side vector
  fread(r, sizeof(*r), n, aux);

  //  increment of left hand side vector
  fread(dr, sizeof(*dr), n, aux);
  
  //  vector of the right hand side from previous step
  fread(fp, sizeof(*fp), n, aux);
  
  //  number of integration points in the problem
  fread (&Tm->tnip, sizeof(Tm->tnip), 1, aux);

  fclose (aux);

  //  data on integration points
  Tm->restore_intpointst_bin (Tp->hdbcont.selelemr, Tp->hdbcont.selother_r, Tp->hdbcont.selother_id);
}



