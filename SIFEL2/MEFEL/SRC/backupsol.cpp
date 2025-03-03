#include "backupsol.h"
#include "hdbcontr.h"
#include "galias.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "iotools.h"
#include "node.h"
#include "element.h"
#include "mathem.h"

#include <stdio.h>
#include <stdlib.h>


/**
  Function saves all necessary data to backup file in the given time step. The process is driven by 
  backup controler Mp->hdbcont. The data can be used for restarting analysis from the given step.
  
  @param r - %vector of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - pointer to the time controller used at solver (if NULL, no time controler data are saved)
  @param n    - number of unknowns (size of %vectors r and fp)

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void solver_save (double *r,double *fp,long ni,double time, double dt, timecontr *tc, long n)
{
  switch (Mp->hdbcont.hdbfmts)
  {
    case text:
      solver_save_text(r, fp, ni, time, dt, tc, n);
      break;
    case binary:
      solver_save_binary(r, fp, ni, time, dt, tc, n);
      break;
    default: 
      print_err("unknown type of backup file format is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function saves all necessary data to the text backup file in the given time step. The process is driven by 
  backup controler Mp->hdbcont. The data can be used for restarting analysis from the given step.
  
  @param r - %vector of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - pointer to the time controller used at solver (if NULL, no time controler data are saved)
  @param n    - number of unknowns (size of %vectors r and fp)

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void solver_save_text (double *r,double *fp,long ni,double time, double dt, timecontr *tc, long n)
{
  switch (Mp->hdbcont.hdbtype)
  {
    case hdbs_single:
    case hdbrs_single:
      solver_save_text_single(r, fp, ni, time, dt, tc, n);
      break;
    case hdbs_multiple:
    case hdbrs_multiple:
      solver_save_text_multiple(r, fp, ni, time, dt, tc, n);
      break;
    default:
      print_err("unknown type of backup is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function saves all necessary data to the binary backup file in the given time step. The process is driven by 
  backup controler Mp->hdbcont. The data can be used for restarting analysis from the given step.
  
  @param r - %vector of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - pointer to the time controller used at solver (if NULL, no time controler data are saved)
  @param n    - number of unknowns (size of %vectors r and fp)

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void solver_save_binary (double *r,double *fp,long ni,double time, double dt, timecontr *tc, long n)
{
  switch (Mp->hdbcont.hdbtype)
  {
    case hdbs_single:
    case hdbrs_single:
      solver_save_binary_single(r, fp, ni, time, dt, tc, n);
      break;
    case hdbs_multiple:
    case hdbrs_multiple:
      solver_save_binary_multiple(r, fp, ni, time, dt, tc, n);
      break;
    default:
      print_err("unknown type of backup is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function saves all necessary data to the single text backup file in the given time step. The process is driven by 
  backup controler Mp->hdbcont. The data can be used for restarting analysis from the given step.
  
  @param r - %vector of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - pointer to the time controller used at solver (if NULL, no time controler data are saved)
  @param n    - number of unknowns (size of %vectors r and fp)

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void solver_save_text_single (double *r,double *fp,long ni,double time, double dt, timecontr *tc, long n)
{
  long i, j, ndofn, ndofe;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];
  int prec = (int)Mp->hdbcont.prec;

  sprintf(name, "%s.%ld.bac",Mp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  fprintf (aux,"\n\n\n\n");

  //  storing time controller data
  if (tc)
    tc->save_txt(aux, prec);
  //  attained solver time
  fprintf (aux,"%.*le\n",prec,time);
  //  attained solver time increment
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
  
  //  vector of the right hand side from previous step
  for (i=0;i<n;i++){
    fprintf (aux,"%.*le\n",prec,fp[i]);
  }
  
  //  number of integration points in the problem
  fprintf (aux,"\n\n %ld\n",Mm->tnip);

  //  data on integration points
  Mm->save_intpoints_txt (aux, Mp->hdbcont.selelems, Mp->hdbcont.selother_s);
  
  // save initial displacements at nodes and elements for growing structures
  if (Mp->tprob == growing_mech_structure)
  {
    for (i=0;i<Mt->nn;i++)
    {
      if (Mt->nodedispl[i] == NULL)
        continue;
      ndofn=Mt->give_ndofn(i); 
      for (j=0; j<ndofn; j++)
        fprintf(aux,"%.*le\n",prec,Mt->nodedispl[i][j]);
    }
    for (i=0;i<Mt->ne;i++)
    {
      ndofe=Mt->give_ndofe(i); 
      for (j=0; j<ndofe; j++)
        fprintf(aux,"%.*le\n",prec,Mt->elements[i].initdispl[j]);
    }
  }

  // save auxiliary integration points if they are defined in the problem
  if (Mm->tnaip>0)
    Mm->save_auxintpoints_txt(aux, Mp->hdbcont.selelems, Mp->hdbcont.selother_s);

  fclose (aux);
  if (Mp->hdbcont.rmold)
  {
    if (Mp->hdbcont.rmold_id>=0)
    {
      sprintf(name, "%s.%ld.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
      remove(name);
    }
    Mp->hdbcont.rmold_id = ni;
  }
}



/**
  Function saves all necessary data to the several text backup files in the given time step.
  Individual quantities are saved into separated files, each file contains given quantity for all
  integration points or all data needed by solver. The file names consist of name specified 
  in the Mp->hdbcont, time step id, description of the saved quantity and .bac suffix. Description
  of the quantity can be:
  solver - file contains date needed by solver
  initdispl - file contains initial displacements at nodes and elements (created only in case of  growing structures)
  strain - array of strains at all integration points
  stress - array of stresses at all integration points
  other  - selected components of eqother array at all integration points
  nonloc - adjacent integration points and their distances at all integration points (used only for nonlocal material models)
           These points are saved only once in the adjacip() function.
  The process is driven by the backup controler Mp->hdbcont. The data can be used for restarting 
  analysis from the given step.
  
  @param r - %vector of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - pointer to the time controller used at solver (if NULL, no time controler data are saved)
  @param n    - number of unknowns (size of %vectors r and fp)

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void solver_save_text_multiple (double *r,double *fp,long ni,double time, double dt, timecontr *tc, long n)
{
  long i, j, ndofn, ndofe;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];
  int prec = (int)Mp->hdbcont.prec;

  sprintf(name, "%s.%ld.solver.bac",Mp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }

  //  storing time controller data
  if (tc)
    tc->save_txt(aux, prec);
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
  
  //  vector of the right hand side from previous step
  for (i=0;i<n;i++){
    fprintf (aux,"%.*le\n",prec,fp[i]);
  }
  
  //  number of integration points in the problem
  fprintf (aux,"\n\n %ld\n",Mm->tnip);

  fclose (aux);

  ////tady je to zakomentovano, protoze se to maze znovu dole na konci teto funkce s ostatnimi soubory
  //if ((Mp->hdbcont.rmold) && (Mp->hdbcont.rmold_id>=0))
  //{
  //sprintf(name, "%s.%ld.solver.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
  //remove(name);
  //}

  //  data on integration points
  Mm->save_intpoints_txt (ni, Mp->hdbcont.selelems, Mp->hdbcont.selother_s);

  // save initial displacements at nodes and elements for growing structures
  if (Mp->tprob == growing_mech_structure)
  {
    sprintf(name, "%s.%ld.initdispl.bac",Mp->hdbcont.hdbnames, ni);
    aux = fopen(name,"wt");
    if (aux==NULL)
    {
      sprintf(emsg, "cannot open backup file %s", name);
      print_err(emsg, __FILE__, __LINE__, __func__);
      abort();
    }
    for (i=0;i<Mt->nn;i++)
    {
      if (Mt->nodedispl[i] == NULL)
        continue;
      ndofn=Mt->give_ndofn(i); 
      for (j=0; j<ndofn; j++)
        fprintf(aux,"%.*le\n",prec,Mt->nodedispl[i][j]);
    }
    for (i=0;i<Mt->ne;i++)
    {
      ndofe=Mt->give_ndofe(i); 
      for (j=0; j<ndofe; j++)
        fprintf(aux,"%.*le\n",prec,Mt->elements[i].initdispl[j]);
    }
  }
 
  //  save data on auxiliary integration points if they are defined in the problem
  if (Mm->tnaip > 0)
    Mm->save_auxintpoints_txt (ni, Mp->hdbcont.selelems, Mp->hdbcont.selother_s);

  if (Mp->hdbcont.rmold)
  {
    if (Mp->hdbcont.rmold_id>=0)
    {
      sprintf(name, "%s.%ld.solver.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
      remove(name);
      if (Mp->tprob == growing_mech_structure){
	sprintf(name, "%s.%ld.initdispl.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
	remove(name);
      }
      sprintf(name, "%s.%ld.strain.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
      remove(name);
      sprintf(name, "%s.%ld.stress.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
      remove(name);
      sprintf(name, "%s.%ld.other.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
      remove(name);
      sprintf(name, "%s.%ld.trfquant.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
      remove(name);
      if (Mm->tnaip > 0)
      {
        sprintf(name, "%s.%ld.astrain.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
        remove(name);
        sprintf(name, "%s.%ld.astress.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
        remove(name);
        sprintf(name, "%s.%ld.aother.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
        remove(name);
      }
    }
    Mp->hdbcont.rmold_id = ni;
  }
}



/**
  Function saves all necessary data to the single binary backup file in the given time step. The process is driven by 
  backup controler Mp->hdbcont. The data can be used for restarting analysis from the given step.
  
  @param r - %vector of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - pointer to the time controller used at solver (if NULL, no time controler data are saved)
  @param n    - number of unknowns (size of %vectors r and fp)

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void solver_save_binary_single (double *r,double *fp,long ni,double time, double dt, timecontr *tc, long n)
{
  long i, ndofn, ndofe;
  FILE *aux;
  char name[FNAMELEN+30];
  char emsg[FNAMELEN+100];

  sprintf(name, "%s.%ld.bac",Mp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }

  //  storing time controller data
  if (tc)
    tc->save_bin(aux);
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
  
  //  vector of the right hand side from previous step
  fwrite(fp, sizeof(*fp), n, aux);
  
  //  number of integration points in the problem
  fwrite (&Mm->tnip, sizeof(Mm->tnip), 1, aux);

  //  data on integration points
  Mm->save_intpoints_bin (aux, Mp->hdbcont.selelems, Mp->hdbcont.selother_s); 

  // save initial displacements at nodes and elements for growing structures
  if (Mp->tprob == growing_mech_structure)
  {
    for (i=0;i<Mt->nn;i++)
    {
      if (Mt->nodedispl[i] == NULL)
        continue;
      ndofn=Mt->give_ndofn(i); 
      fwrite(Mt->nodedispl[i], sizeof(*Mt->nodedispl[i]), ndofn, aux);
    }
    for (i=0;i<Mt->ne;i++)
    {
      ndofe=Mt->give_ndofe(i); 
      fwrite(Mt->elements[i].initdispl, sizeof(*Mt->elements[i].initdispl), ndofe, aux);
    }
  }

  //  save data on auxiliary integration points if they are defined in the problem
  if (Mm->tnaip > 0)
    Mm->save_auxintpoints_bin (aux, Mp->hdbcont.selelems, Mp->hdbcont.selother_s); 

  fclose (aux);
  if (Mp->hdbcont.rmold)
  {
    if (Mp->hdbcont.rmold_id>=0)
    {
      sprintf(name, "%s.%ld.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
      remove(name);
    }
    Mp->hdbcont.rmold_id = ni;
  }
}



/**
  Function saves all necessary data to the several binary backup files in the given time step.
  Individual quantities are saved into separated files, each file contains given quantity for all
  integration points or all data needed by solver. The file names consist of name specified 
  in the Mp->hdbcont, time step id, description of the saved quantity and .bac suffix. Description
  of the quantity can be:
  solver - file contains date needed by solver
  strain - array of strains at all integration points
  stress - array of stresses at all integration points
  other  - selected components of eqother array at all integration points
  nonloc - adjacent integration points and their distances at all integration points (used only for nonlocal material models)
           These points are saved only once in the adjacip() function.
  The process is driven by the backup controler Mp->hdbcont. The data can be used for restarting 
  analysis from the given step.
  
  @param r - %vector of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - pointer to the time controller used at solver (if NULL, no time controler data are saved)
  @param n    - number of unknowns (size of %vectors r and fp)

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void solver_save_binary_multiple (double *r,double *fp,long ni,double time, double dt, timecontr *tc, long n)
{
  long i, ndofn, ndofe;
  FILE *aux;
  char name[FNAMELEN+30];
  char emsg[FNAMELEN+100];
  
  sprintf(name, "%s.%ld.solver.bac",Mp->hdbcont.hdbnames, ni);
  aux = fopen(name,"wb");
  if (aux==NULL)
    {
      sprintf(emsg, "cannot open backup file %s", name);
      print_err(emsg, __FILE__, __LINE__, __func__);
      abort();
    }
  
  //  storing time controller data
  if (tc)
    tc->save_bin(aux);
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
  
  //  vector of the right hand side from previous step
  fwrite(fp, sizeof(*fp), n, aux);
  
  //  number of integration points in the problem
  fwrite (&Mm->tnip, sizeof(Mm->tnip), 1, aux);
  
  fclose (aux);
  
  //tady je to zakomentovano, protoze se to maze znovu dole na konci teto funkce s ostatnimi soubory
  //if ((Mp->hdbcont.rmold) && (Mp->hdbcont.rmold_id>=0))
  //{
  //sprintf(name, "%s.%ld.solver.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
  //remove(name);
  //}
  
  //  data on integration points
  Mm->save_intpoints_bin (ni, Mp->hdbcont.selelems, Mp->hdbcont.selother_s);
  
  // save initial displacements at nodes and elements for growing structures
  if (Mp->tprob == growing_mech_structure)
    {
      sprintf(name, "%s.%ld.initdispl.bac",Mp->hdbcont.hdbnames, ni);
      aux = fopen(name,"wb");
      if (aux==NULL)
	{
	  sprintf(emsg, "cannot open backup file %s", name);
	  print_err(emsg, __FILE__, __LINE__, __func__);
	  abort();
	}
      for (i=0;i<Mt->nn;i++)
	{
	  if (Mt->nodedispl[i] == NULL)
	    continue;
	  ndofn=Mt->give_ndofn(i); 
	  fwrite(Mt->nodedispl[i], sizeof(*Mt->nodedispl[i]), ndofn, aux);
	}
      for (i=0;i<Mt->ne;i++)
	{
	  ndofe=Mt->give_ndofe(i); 
	  fwrite(Mt->elements[i].initdispl, sizeof(*Mt->elements[i].initdispl), ndofe, aux);
	}
    }
  
  //  save data on auxiliary integration points if they are defined in the problem
  if (Mm->tnaip > 0)
    Mm->save_auxintpoints_bin (ni, Mp->hdbcont.selelems, Mp->hdbcont.selother_s);

  if (Mp->hdbcont.rmold)
    {
      if (Mp->hdbcont.rmold_id>=0)
	{
	  sprintf(name, "%s.%ld.solver.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
	  remove(name);
	  if (Mp->tprob == growing_mech_structure){
	    sprintf(name, "%s.%ld.initdispl.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
	    remove(name);
	  }
	  sprintf(name, "%s.%ld.strain.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
	  remove(name);
	  sprintf(name, "%s.%ld.stress.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
	  remove(name);
	  sprintf(name, "%s.%ld.other.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
	  remove(name);
	  sprintf(name, "%s.%ld.trfquant.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
	  remove(name);
          if (Mm->tnaip > 0)
          {
            sprintf(name, "%s.%ld.astrain.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
            remove(name);
            sprintf(name, "%s.%ld.astress.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
            remove(name);
            sprintf(name, "%s.%ld.aother.bac",Mp->hdbcont.hdbnames,Mp->hdbcont.rmold_id);
            remove(name);
          }
	}
      Mp->hdbcont.rmold_id = ni;
    }
}



/**
  Function restores all necessary data from the backup file. The process is driven by 
  backup controler Mp->hdbcont. Data from solver are passed through parameters of the function. 
  
  @param r - %vector of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - pointer to the time controller used at solver (if NULL, no time controler data are saved)
  @param n    - number of unknowns (size of %vectors r and fp)

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void solver_restore (double *r,double *fp,long &ni,double &time, double &dt, timecontr *tc, long &n)
{
  switch (Mp->hdbcont.hdbfmtr)
  {
    case text:
      solver_restore_text(r, fp, ni, time, dt, tc, n);
      break;
    case binary:
      solver_restore_binary(r, fp, ni, time, dt, tc, n);
      break;
    default: 
      print_err("unknown type of backup file format is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function restores all necessary data from the text backup file. The process is driven by 
  backup controler Mp->hdbcont. Data from solver are passed through parameters of the function. 
  
  @param r - %vector of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - pointer to the time controller used at solver (if NULL, no time controler data are saved)
  @param n    - number of unknowns (size of %vectors r and fp)

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void solver_restore_text (double *r,double *fp,long &ni,double &time, double &dt, timecontr *tc, long &n)
{
  switch (Mp->hdbcont.hdbtype)
  {
    case hdbr_single:
    case hdbrs_single:
      solver_restore_text_single(r, fp, ni, time, dt, tc, n);
      break;
    case hdbr_multiple:
    case hdbrs_multiple:
      solver_restore_text_multiple(r, fp, ni, time, dt, tc, n);
      break;
    default:
      print_err("unknown type of backup is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function restores all necessary data from the binary backup file. The process is driven by 
  backup controler Mp->hdbcont. Data from solver are passed through parameters of the function. 
  
  @param r - %vector of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - pointer to the time controller used at solver (if NULL, no time controler data are saved)
  @param n    - number of unknowns (size of %vectors r and fp)

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void solver_restore_binary (double *r,double *fp,long &ni,double &time, double &dt, timecontr *tc, long &n)
{
  switch (Mp->hdbcont.hdbtype)
  {
    case hdbr_single:
    case hdbrs_single:
      solver_restore_binary_single(r, fp, ni, time, dt, tc, n);
      break;
    case hdbr_multiple:
    case hdbrs_multiple:
      solver_restore_binary_multiple(r, fp, ni, time, dt, tc, n);
      break;
    default:
      print_err("unknown type of backup is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function restores all necessary data from the one huge text backup file. The process is driven by 
  backup controler Mp->hdbcont. Data from solver are passed through parameters of the function. 
  
  @param r - %vector of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - pointer to the time controller used at solver (if NULL, no time controler data are saved)
  @param n    - number of unknowns (size of %vectors r and fp)

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void solver_restore_text_single (double *r, double *fp, long &ni, double &time, double &dt, timecontr *tc, long &n)
{
  long i, j, ndofn, ndofe;
  FILE *aux;
  char emsg[FNAMELEN+100];

  aux = fopen(Mp->hdbcont.hdbnamer,"rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", Mp->hdbcont.hdbnamer);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }

  // restorage of time controller data
  if (tc)
    tc->restore_txt(aux);  
  //  attained time of solver
  fscanf (aux,"%le",&time);
  //  attained time increment of solver
  fscanf (aux,"%le",&dt);
  //  attained number of steps
  fscanf (aux,"%ld",&ni);
  //  number of degrees of freedom
  fscanf (aux,"%ld",&n);

  if (n!=Ndofm){
    print_err("number of DOFs in backup file and in the actual problem is not same\n", __FILE__, __LINE__, __func__);
    abort();
  }

  //  left hand side vector
  for (i=0;i<n;i++){
    fscanf (aux,"%le",r+i);
  }
  //  vector of the right hand side from previous step
  for (i=0;i<n;i++){
    fscanf (aux,"%le",fp+i);
  }
  
  //  number of integration points in the problem
  fscanf (aux,"%ld",&i);
  
  if (i!=Mm->tnip){
    print_err("number of integration points in backup file and in the actual problem is not same\n", __FILE__, __LINE__, __func__);
    abort();
  }

  //  data on integration points
  Mm->restore_intpoints_txt (aux, Mp->hdbcont.selelemr, Mp->hdbcont.selother_r, Mp->hdbcont.selother_id);

  // restore initial displacements at nodes and elements for growing structures
  if (Mp->tprob == growing_mech_structure)
  {
    for (i=0;i<Mt->nn;i++)
    {
      if (Mt->nodedispl[i] == NULL)
        continue;
      ndofn = Mt->give_ndofn(i); 
      for (j=0; j<ndofn; j++)
        fscanf(aux,"%le", Mt->nodedispl[i]+j);
    }
    for (i=0;i<Mt->ne;i++)
    {
      ndofe = Mt->give_ndofe(i); 
      for (j=0; j<ndofe; j++)
        fscanf(aux,"%le", Mt->elements[i].initdispl+j);
    }
  }

  // restore data on auxiliary integration points if they are defined in the problem
  if (Mm->tnaip > 0)
    Mm->restore_auxintpoints_txt (aux, Mp->hdbcont.selelemr, Mp->hdbcont.selother_r, Mp->hdbcont.selother_id);

  fclose(aux);
}



/**
  Function restores all necessary data from the separated text backup file. The process is driven by 
  backup controler Mp->hdbcont. Data from solver are passed through parameters of the function. 
  Individual quantities have to be saved into separated files, each file contains given quantity for all
  integration points or all data needed by solver. The file names consist of name specified 
  in the Mp->hdbcont, time step id, description of the saved quantity and .bac suffix. Description
  of the quantity can be:
  solver - file contains date needed by solver
  strain - array of strains at all integration points
  stress - array of stresses at all integration points
  other  - selected components of eqother array at all integration points
  nonloc - adjacent integration points and their distances at all integration points (used only for nonlocal material models)
           These points are saved only once in the adjacip() function.
  The process is driven by the backup controler Mp->hdbcont. The data can be used for restarting 
  analysis from the given step.
  
  @param r - %vector of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - pointer to the time controller used at solver (if NULL, no time controler data are saved)
  @param n    - number of unknowns (size of %vectors r and fp)

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void solver_restore_text_multiple (double *r,double *fp,long &ni,double &time, double &dt, timecontr *tc, long &n)
{
  long i, j, ndofn, ndofe;
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];

  sprintf(name, "%s.solver.bac",Mp->hdbcont.hdbnamer);
  aux = fopen(name,"rt");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }
  
  // restorage of time controller data
  if (tc)
    tc->restore_txt(aux);  
  //  attained time of solver
  fscanf (aux,"%le",&time);
  //  attained time increment of solver
  fscanf (aux,"%le",&dt);
  //  attained number of steps
  fscanf (aux,"%ld",&ni);
  //  number of degrees of freedom
  fscanf (aux,"%ld",&n);

  if (n!=Ndofm){
    print_err("number of DOFs in backup file and in the actual problem is not same\n", __FILE__, __LINE__, __func__);
    abort();
  }

  //  left hand side vector
  for (i=0;i<n;i++){
    fscanf (aux,"%le",r+i);
  }
  //  vector of the right hand side from previous step
  for (i=0;i<n;i++){
    fscanf (aux,"%le",fp+i);
  }
  
  //  number of integration points in the problem
  fscanf (aux,"%ld",&i);
  
  if (i!=Mm->tnip){
    print_err("number of integration points in backup file and in the actual problem is not same\n", __FILE__, __LINE__, __func__);
    abort();
  }

  fclose(aux);

  //  data on integration points
  Mm->restore_intpoints_txt (Mp->hdbcont.selelemr, Mp->hdbcont.selother_r, Mp->hdbcont.selother_id);

  // restore initial displacements at nodes and elements for growing structures
  if (Mp->tprob == growing_mech_structure)
  {
    sprintf(name, "%s.%ld.initdispl.bac",Mp->hdbcont.hdbnames, ni);
    aux = fopen(name,"rt");
    if (aux==NULL)
    {
      sprintf(emsg, "cannot open backup file %s", name);
      print_err(emsg, __FILE__, __LINE__, __func__);
      abort();
    }
    for (i=0;i<Mt->nn;i++)
    {
      if (Mt->nodedispl[i] == NULL)
        continue;
      ndofn = Mt->give_ndofn(i); 
      for (j=0; j<ndofn; j++)
        fscanf(aux,"%le", Mt->nodedispl[i]+j);
    }
    for (i=0;i<Mt->ne;i++)
    {
      ndofe = Mt->give_ndofe(i); 
      for (j=0; j<ndofe; j++)
        fscanf(aux,"%le", Mt->elements[i].initdispl+j);
    }
  }
  // restore data on auxiliary integration points if they are defined in the problem
  if (Mm->tnaip > 0)
    Mm->restore_auxintpoints_txt (Mp->hdbcont.selelemr, Mp->hdbcont.selother_r, Mp->hdbcont.selother_id);
}



/**
  Function restores all necessary data from the one huge binary backup file. The process is driven by 
  backup controler Mp->hdbcont. Data from solver are passed through parameters of the function. 
  
  @param r - %vector of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - pointer to the time controller used at solver (if NULL, no time controler data are saved)
  @param n    - number of unknowns (size of %vectors r and fp)

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void solver_restore_binary_single (double *r,double *fp,long &ni,double &time, double &dt, timecontr *tc, long &n)
{
  long i, ndofn, ndofe;
  FILE *aux;
  char emsg[FNAMELEN+100];

  aux = fopen(Mp->hdbcont.hdbnamer,"rb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", Mp->hdbcont.hdbnamer);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }

  // restorage of time controller data
  if (tc)
    tc->restore_bin(aux);  
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
  
  //  vector of the right hand side from previous step
  fread(fp, sizeof(*fp), n, aux);
  
  //  number of integration points in the problem
  //  fread (&Mm->tnip, sizeof(Mm->tnip), 1, aux);
  fread (&i, sizeof(i), 1, aux);

  if (i!=Mm->tnip){
    print_err("number of integration points in backup file and in the actual problem is not same\n", __FILE__, __LINE__, __func__);
    abort();
  }

  //  data on integration points
  Mm->restore_intpoints_bin (aux, Mp->hdbcont.selelemr, Mp->hdbcont.selother_r, Mp->hdbcont.selother_id);

  // restore initial displacements at nodes and elements for growing structures
  if (Mp->tprob == growing_mech_structure)
  {
    for (i=0;i<Mt->nn;i++)
    {
      if (Mt->nodedispl[i] == NULL)
        continue;
      ndofn = Mt->give_ndofn(i); 
      fread(Mt->nodedispl[i], sizeof(*Mt->nodedispl[i]), ndofn, aux);
    }
    for (i=0;i<Mt->ne;i++)
    {
      ndofe = Mt->give_ndofe(i); 
      fread(Mt->elements[i].initdispl, sizeof(*Mt->elements[i].initdispl), ndofe, aux);
    }
  }

  // restore data on auxiliary integration points if they are defined in the problem
  if (Mm->tnaip > 0)
    Mm->restore_auxintpoints_bin (aux, Mp->hdbcont.selelemr, Mp->hdbcont.selother_r, Mp->hdbcont.selother_id);

  fclose (aux);
}



/**
  Function restores all necessary data from the separated text backup file. The process is driven by 
  backup controler Mp->hdbcont. Data from solver are passed through parameters of the function. 
  Individual quantities have to be saved into separated files, each file contains given quantity for all
  integration points or all data needed by solver. The file names consist of name specified 
  in the Mp->hdbcont, time step id, description of the saved quantity and .bac suffix. Description
  of the quantity can be:
  solver - file contains date needed by solver
  strain - array of strains at all integration points
  stress - array of stresses at all integration points
  other  - selected components of eqother array at all integration points
  nonloc - adjacent integration points and their distances at all integration points (used only for nonlocal material models)
           These points are saved only once in the adjacip() function.
  The process is driven by the backup controler Mp->hdbcont. The data can be used for restarting 
  analysis from the given step.
  
  @param r - %vector of unknowns
  @param fp - %vector of actual right hand side
  @param ni - time step id
  @param time - actual time of saved step of solver
  @param dt   - actual time increment of saved step of solver
  @param tc   - pointer to the time controller used at solver (if NULL, no time controler data are saved)
  @param n    - number of unknowns (size of %vectors r and fp)

  @return The function does not return anything.

  Created by Tomas Koudelka, 9.2008
*/
void solver_restore_binary_multiple (double *r,double *fp,long &ni,double &time, double &dt, timecontr *tc, long &n)
{
  long i, ndofn, ndofe;
  FILE *aux;
  char name[FNAMELEN+30];
  char emsg[FNAMELEN+100];

  sprintf(name, "%s.solver.bac",Mp->hdbcont.hdbnamer);
  aux = fopen(name,"rb");
  if (aux==NULL)
  {
    sprintf(emsg, "cannot open backup file %s", name);
    print_err(emsg, __FILE__, __LINE__, __func__);
    abort();
  }

  // restorage of time controller data
  if (tc)
    tc->restore_bin(aux);  
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
  
  //  vector of the right hand side from previous step
  fread(fp, sizeof(*fp), n, aux);
  
  //  number of integration points in the problem
  //  fread (&Mm->tnip, sizeof(Mm->tnip), 1, aux);
  fread (&i, sizeof(i), 1, aux);

  if (i!=Mm->tnip){
    print_err("number of integration points in backup file and in the actual problem is not same\n", __FILE__, __LINE__, __func__);
    abort();
  }

  fclose (aux);

  //  data on integration points
  Mm->restore_intpoints_bin (Mp->hdbcont.selelemr, Mp->hdbcont.selother_r, Mp->hdbcont.selother_id);

  // restore initial displacements at nodes and elements for growing structures
  if (Mp->tprob == growing_mech_structure)
  {
    sprintf(name, "%s.%ld.initdispl.bac",Mp->hdbcont.hdbnames, ni);
    aux = fopen(name,"rb");
    if (aux==NULL)
    {
      sprintf(emsg, "cannot open backup file %s", name);
      print_err(emsg, __FILE__, __LINE__, __func__);
      abort();
    }
    for (i=0;i<Mt->nn;i++)
    {
      if (Mt->nodedispl[i] == NULL)
        continue;
      ndofn = Mt->give_ndofn(i); 
      fread(Mt->nodedispl[i], sizeof(*Mt->nodedispl[i]), ndofn, aux);
    }
    for (i=0;i<Mt->ne;i++)
    {
      ndofe = Mt->give_ndofe(i); 
      fread(Mt->elements[i].initdispl, sizeof(*Mt->elements[i].initdispl), ndofe, aux);
    }
  }

  // restore data on auxiliary integration points if they are defined in the problem
  if (Mm->tnaip > 0)
    Mm->restore_auxintpoints_bin (Mp->hdbcont.selelemr, Mp->hdbcont.selother_r, Mp->hdbcont.selother_id);
}
