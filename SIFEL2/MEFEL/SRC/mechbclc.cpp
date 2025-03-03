#include <stdlib.h>
#include <string.h>

#include "mechbclc.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "alias.h"
#include "loadcase.h"
#include "dloadcase.h"
#include "inicd.h"
#include "axisrotrec.h"
#include "vector.h"
#include "matrix.h"
#include "elemhead.h"
#include "gtopology.h"
#include "element.h"
#include "node.h"
#include "elemswitch.h"
#include "intpoints.h"
#include "sequent.h"
#include "vector.h"


/**
  Constructor initializes data members to zero or default values.

  Created by JK,
  Modified by TKo, TKr,
*/
mechbclc::mechbclc ()
{
  //  number of load cases
  nlc=0L;
  //  number of initial conditions
  nico = 0L;
  //  load case
  lc=NULL;
  //  time dependent load case
  dlc = NULL;
  // initial conditions
  ico = NULL;
  // prescribed initial displacements by rotations
  arotrec = NULL;   

  ncsum = 0L;
  sumcomp         = NULL;
  reactsumcomp    = NULL;
  pd_reactsumcomp = NULL;

  //  the number of general functions describing the eigenstrains
  ngfes=0L;
  //  array of general functions describing the eigenstrains
  eigstrfun = NULL;
  ssst = (strastrestate)0;
}



/**
  Destructor releases allocated memory of the mechbclc object.

  Created by JK,
  Modified by TKo, TKr,
*/
mechbclc::~mechbclc ()
{
  if (lc!=NULL)
    delete [] lc;
  if (dlc!=NULL)
    delete [] dlc;
  if (ico!=NULL)
    delete [] ico;
  if (sumcomp!=NULL)
    delete [] sumcomp;
  if (reactsumcomp!=NULL)
    delete [] reactsumcomp;
  if (pd_reactsumcomp!=NULL)
    delete [] pd_reactsumcomp;
  if (eigstrfun!=NULL)
    delete [] eigstrfun;
  if (arotrec)
    delete [] arotrec;
}



/**
  Function reads data about boundary conditions and load cases.
   
  @param in - pointer to the opened XFILE
   
  @return The function does not return anything.

  Created by JK,
  Modified by Tomas Koudelka,
*/
void mechbclc::read (XFILE *in)
{
  long i;
  vector nodval,ipval;
  
  switch (Mp->tprob){
  case linear_statics:{
    xfscanf (in,"%k%ld","number_of_load_cases",&nlc);
    if (nlc<0){
      print_err("wrong number of load cases", __FILE__, __LINE__, __func__);
      abort ();
    }
    lc = new loadcase [nlc];
    for (i=0;i<nlc;i++){
      lc[i].read (in);
    }
    break;
  }
  case eigen_dynamics:{
    if (Mp->straincomp==1 || Mp->stresscomp==1)
      nlc=Mp->eigsol.neigv;
    break;
  }
  case linear_stability:{
    xfscanf (in,"%k%ld","number_of_load_cases",&nlc);
    if (nlc<0){
      print_err("wrong number of load cases", __FILE__, __LINE__, __func__);
      abort ();
    }
    lc = new loadcase [nlc];
    for (i=0;i<nlc;i++){
      lc[i].read (in);
    }
    break;
  }
  case mat_nonlinear_statics:{
    xfscanf (in,"%k%ld","number_of_load_cases",&nlc);
    if (nlc!=1){
      print_err("wrong number of load cases", __FILE__, __LINE__, __func__);
      abort ();
    }
    lc = new loadcase [2*nlc];
    //  vector of proportional load is read first
    //  vector of constant (non-proportional) load is read second
    for (i=0;i<2*nlc;i++){
      lc[i].read (in);
    }
    readinic(in);
    break;
  }
  case geom_nonlinear_statics:{
    xfscanf (in,"%k%ld","number_of_load_cases",&nlc);
    if (nlc!=1){
      print_err("wrong number of load cases", __FILE__, __LINE__, __func__);
      abort ();
    }
    lc = new loadcase [2*nlc];
    for (i=0;i<2*nlc;i++){
      lc[i].read (in);
    }
    readinic(in);
    break;
  }
  case forced_dynamics:{
    xfscanf (in,"%k%ld","number_of_load_cases",&nlc);
    if (nlc<0){
      print_err("wrong number of load cases", __FILE__, __LINE__, __func__);
      abort ();
    }
    lc = new loadcase [nlc];
    dlc = new dloadcase [nlc];
    for (i=0;i<nlc;i++){
      lc[i].read (in);
      dlc[i].read (in);
    }
    readinic(in);
    break;
  }
  case mech_timedependent_prob:{
    xfscanf (in,"%k%ld","number_of_load_cases",&nlc);
    if (nlc!=1){
      print_err("wrong number of load cases", __FILE__, __LINE__, __func__);
      abort ();
    }
    dlc = new dloadcase [nlc];
    for (i=0;i<nlc;i++){
      dlc[i].read (in);
    }
    //  initial values for material models
    readinic(in);
    break;
  }
  case growing_mech_structure:{
    xfscanf (in,"%k%ld","number_of_load_cases",&nlc);
    if (nlc!=1){
      print_err("wrong number of load cases", __FILE__, __LINE__, __func__);
      abort ();
    }
    dlc = new dloadcase [nlc];
    for (i=0;i<nlc;i++){
      dlc[i].read (in);
    }
    readinic(in);
    if (Mp->rot_inidispl == yes)
      read_rotinidispl(in);
    break;
  }
  case earth_pressure:{
    xfscanf (in,"%k%ld","number_of_load_cases",&nlc);
    if (Mp->tepsol < epplast_sol)
    {
      lc = new loadcase [nlc];
      dlc = new dloadcase [nlc];
      for (i=0;i<nlc;i++){
        lc[i].read (in);
        dlc[i].read (in);
      }
      readinic(in);
    }
    else
    {
      lc = new loadcase [2*nlc];
      for (i=0;i<2*nlc;i++){
        lc[i].read (in);
      }
      readinic(in);
    }
    break;
  }
  case layered_linear_statics:{
    xfscanf (in,"%k%ld","number_of_load_cases",&nlc);
    if (nlc<0){
      print_err("wrong number of load cases", __FILE__, __LINE__, __func__);
      abort ();
    }
    lc = new loadcase [nlc];
    for (i=0;i<nlc;i++){
      lc[i].read (in);
    }
    break;
  }
    
  case lin_floating_subdomain:{
    xfscanf (in,"%k%ld","number_of_load_cases",&nlc);
    if (nlc<0){
      print_err("wrong number of load cases", __FILE__, __LINE__, __func__);
      abort ();
    }
    lc = new loadcase [nlc];
    for (i=0;i<nlc;i++){
      lc[i].read (in);
    }
    break;
  }

  case nonlin_floating_subdomain:{
    xfscanf (in,"%k%ld","number_of_load_cases",&nlc);
    if (nlc!=1){
      print_err("wrong number of load cases", __FILE__, __LINE__, __func__);
      abort ();
    }
    lc = new loadcase [2*nlc];
    for (i=0;i<2*nlc;i++){
      lc[i].read (in);
    }
    break;
  }

  case nonlinear_dynamics:{
    xfscanf (in,"%k%ld","number_of_load_cases",&nlc);
    if (nlc!=1){
      print_err("wrong number of load cases", __FILE__, __LINE__, __func__);
      abort ();
    }
    lc = new loadcase [nlc];
    dlc = new dloadcase [nlc];
    for (i=0;i<nlc;i++){
      lc[i].read (in);
      dlc[i].read (in);
    }
    break;
  }
    
  case hemivar_inequal:{
    xfscanf (in,"%k%ld","number_of_load_cases",&nlc);
    if (nlc<0){
      print_err("wrong number of load cases", __FILE__, __LINE__, __func__);
      abort ();
    }
    lc = new loadcase [nlc];
    for (i=0;i<nlc;i++){
      lc[i].read (in);
    }
    break;
  }
    
  case load_balancing:{
    xfscanf (in,"%k%ld","number_of_load_cases",&nlc);
    if (nlc<0){
      print_err("wrong number of load cases", __FILE__, __LINE__, __func__);
      abort ();
    }
    lc = new loadcase [nlc];
    for (i=0;i<nlc;i++){
      lc[i].read (in);
    }
    break;
  }

  default:{
    print_err("unknown problem type", __FILE__, __LINE__, __func__);
  }
  }
  
  // **********************
  //  eigenstrain reading
  // **********************
  read_eigenstrains (in);
  
}


/**
  Function prints data about boundary conditions and load cases.
   
  @param out - pointer to the opened FILE
   
  @return The function does not return anything.

  TKr, 07/02/2013 according to read (XFILE *in)
*/
void mechbclc::print (FILE *out)
{
  long i;
  
  fprintf (out,"\n number of load cases:");
  switch (Mp->tprob){
  case linear_statics:{
    fprintf (out,"\n  %ld",nlc);
    for (i=0;i<nlc;i++){
      lc[i].print (out);
    }
    break;
  }
  case eigen_dynamics:{
    break;
  }
  case linear_stability:{
    fprintf (out,"\n  %ld",nlc);
    for (i=0;i<nlc;i++){
      lc[i].print (out);
    }
    break;
  }
  case mat_nonlinear_statics:{
    fprintf (out,"\n  %ld",nlc);
    for (i=0;i<2*nlc;i++){
      lc[i].print (out);
    }
    printinic(out);
    break;
  }
  case geom_nonlinear_statics:{
    fprintf (out,"\n  %ld",nlc);
    for (i=0;i<2*nlc;i++){
      lc[i].print (out);
    }
    printinic(out);
    break;
  }
  case forced_dynamics:{
    fprintf (out,"\n  %ld",nlc);
    for (i=0;i<nlc;i++){
      lc[i].print (out);
      dlc[i].print (out);
    }
    break;
  }
  case mech_timedependent_prob:{
    fprintf (out,"\n  %ld",nlc);
    for (i=0;i<nlc;i++){
      dlc[i].print (out);
    }
    printinic(out);
    break;
  }
  case growing_mech_structure:{
    fprintf (out,"\n  %ld",nlc);
    for (i=0;i<nlc;i++){
      dlc[i].print (out);
    }
    printinic(out);
    break;
  }
  case earth_pressure:{
    fprintf (out,"\n  %ld",nlc);
    if (Mp->tepsol < epplast_sol)
    {
      for (i=0;i<nlc;i++){
        lc[i].print (out);
        dlc[i].print (out);
      }
      printinic(out);
    }
    else
    {
      for (i=0;i<2*nlc;i++){
        lc[i].print (out);
      }
      printinic(out);
    }
    break;
  }
  case layered_linear_statics:{
    fprintf (out,"\n  %ld",nlc);
    for (i=0;i<nlc;i++){
      lc[i].print (out);
    }
    break;
  }
    
  case lin_floating_subdomain:{
    fprintf (out,"\n  %ld",nlc);
    for (i=0;i<nlc;i++){
      lc[i].print (out);
    }
    break;
  }

  case nonlin_floating_subdomain:{
    fprintf (out,"\n  %ld",nlc);
    for (i=0;i<2*nlc;i++){
      lc[i].print (out);
    }
    break;
  }

  case nonlinear_dynamics:{
    fprintf (out,"\n  %ld",nlc);
    for (i=0;i<nlc;i++){
      lc[i].print (out);
      dlc[i].print (out);
    }
    break;
  }
    
  case hemivar_inequal:{
    fprintf (out,"\n  %ld",nlc);
    for (i=0;i<nlc;i++){
      lc[i].print (out);
    }
    break;
  }
    
  case load_balancing:{
    fprintf (out,"\n  %ld",nlc);
    for (i=0;i<nlc;i++){
      lc[i].print (out);
    }
    break;
  }

  default:{
    print_err("unknown problem type", __FILE__, __LINE__, __func__);
  }
  }
  
  // **********************
  //  eigenstrain printing
  // **********************
  print_eigenstrains (out);
  
}


/**
   function reads eigenstrains
   
   @param in - input file
   
   JK, 8.11.2012
*/
void mechbclc::read_eigenstrains (XFILE *in)
{
  long i,j,k;
  long ncomp,nip,ipid;
  long *id=NULL;
  
  //  Mp->eigstrains = 0 - eigenstrains are not defined
  //  Mp->eigstrains = 1 - eigenstrains are defined for the whole problem
  //  Mp->eigstrains = 2 - eigenstrains are defined on each element
  //  Mp->eigstrains = 3 - eigenstrains are defined outside of MEFEL (for the coupled problems)
  //  Mp->eigstrains = 4 - eigenstresses are defined for the whole problem
  //  Mp->eigstrains = 5 - eigenstresses are defined on each element
  xfscanf (in,"%k%ld","eigenstrains",&(Mp->eigstrains));
  if (Mp->eigstrains<0 || Mp->eigstrains>5){
    print_err("wrong definition of eigenstrains in file %s (line %ld)", __FILE__, __LINE__, __func__, in->fname, in->line);
    abort ();
  }
  
  if (Mespr == 1){
    if (Mp->eigstrains==0)
      fprintf(stdout, "\n eigenstrains are not defined");
    if (Mp->eigstrains==1)
      fprintf(stdout, "\n eigenstrains are defined for the whole problem by a single set of values");
    if (Mp->eigstrains==2)
      fprintf(stdout, "\n eigenstrains are defined independently for each element");
    if (Mp->eigstrains==3)
      fprintf(stdout, "\n eigenstrains are defined independently in TRFEL");
    if (Mp->eigstrains==4)
      fprintf(stdout, "\n eigenstresses are defined for the whole problem by a single set of values");
    if (Mp->eigstrains==5)
      fprintf(stdout, "\n eigenstresses are defined independently for each element");
  }
  
  
  if (Mp->eigstrains>0){
    //  eigenstrains are defined
    
    //  array of eigenstrain id
    if (Mm->eigstrid==NULL)
      Mm->eigstrid = new long* [Mm->tnip];
    //  array of eigenstrains in integration points
    if (Mm->eigstrains==NULL)
      Mm->eigstrains = new double* [Mm->tnip];
    //  array of eigenstresses in integration points
    if (Mm->eigstresses==NULL)
      Mm->eigstresses = new double* [Mm->tnip];
    
    id = new long [6];
  }
  
  
  // *************************************************
  //  eigenstrains are defined for the whole problem
  // *************************************************
  if (Mp->eigstrains==1){
    //  type of strain/stress state
    //  strastre = 1 - bar
    //  strastre = 10 - plane stress
    //  strastre = 30 - space strain/stress
    xfscanf (in,"%m", &strastrestate_kwdset, (int *)&ssst);
    
    switch (ssst){
    case bar:{
      xfscanf (in,"%ld",id);
      ncomp = 1;
      break;
    }
    case planestress:{
      xfscanf (in,"%ld %ld %ld %ld",id,id+1,id+2,id+3);
      ncomp = 4;
      break;
    }
    case planestrain:{
      xfscanf (in,"%ld %ld %ld %ld",id,id+1,id+2,id+3);
      ncomp = 4;
      break;
    }
    case axisymm:{
      xfscanf(in,"%ld %ld %ld %ld", id, id+1, id+2, id+3);
      ncomp = 4;
      break;
    }
    case spacestress:{
      xfscanf (in,"%ld %ld %ld %ld %ld %ld",id,id+1,id+2,id+3,id+4,id+5);
      ncomp = 6;
      break;
    }
    default:{
      print_err("unknown type of strain-stress state", __FILE__, __LINE__, __func__);
      abort();
    }
    }//  end of the switch (strastre){
    
    for (i=0;i<Mm->tnip;i++){
      //  array of eigenstrains in integration points
      Mm->eigstrains[i] = new double [ncomp];
      memset(Mm->eigstrains[i], 0, sizeof(*Mm->eigstrains[i])*ncomp);
      //  array of eigenstresses in integration points
      Mm->eigstresses[i] = new double [ncomp];
      memset(Mm->eigstresses[i], 0, sizeof(*Mm->eigstresses[i])*ncomp);
      //  array of eigenstrain id
      Mm->eigstrid[i] = new long [ncomp];
      
      for (j=0;j<ncomp;j++){
        Mm->eigstrid[i][j]=id[j]-1;
      }
    }
    
  }//  end of the statement if (Mp->eigstrains==1){
  
  
  // ***************************************
  //  eigenstrains are defined on elements
  // ***************************************
  if (Mp->eigstrains==2){
    
    //  loop over the number of elements
    for (i=0;i<Mt->ne;i++){
      //  strain/stress state on element
      ssst = Mt->give_ssst (i,0);
      
      switch (ssst){
      case bar:{
        xfscanf (in,"%ld",id);
        ncomp = 1;
        break;
      }
      case planestress:{
	xfscanf (in,"%ld %ld %ld %ld",id,id+1,id+2,id+3);
	ncomp = 4;
	break;
      }
      case planestrain:{
        xfscanf (in,"%ld %ld %ld %ld",id,id+1,id+2,id+3);
        ncomp = 4;
        break;
      }
      case axisymm:{
        xfscanf(in,"%ld %ld %ld %ld", id, id+1, id+2, id+3);
        ncomp = 4;
        break;
      }
      case spacestress:{
	xfscanf (in,"%ld %ld %ld %ld %ld %ld",id,id+1,id+2,id+3,id+4,id+5);
	ncomp = 6;
	break;
      }
      default:{
	print_err("unknown type of strain-stress state", __FILE__, __LINE__, __func__);
        abort();
      }
      }//  end of the switch (strastre){
      
      //  number of integration points on element
      nip = Mt->give_tnip (i);
      //  number of the first integration point
      ipid=Mt->elements[i].ipp[0][0];
      
      //  loop over the number of integration points on element
      for (j=0;j<nip;j++){
	Mm->eigstrains[ipid] = new double [ncomp];
        memset(Mm->eigstrains[ipid], 0, sizeof(*Mm->eigstrains[ipid])*ncomp);
	Mm->eigstresses[ipid] = new double [ncomp];
        memset(Mm->eigstresses[ipid], 0, sizeof(*Mm->eigstresses[ipid])*ncomp);
	Mm->eigstrid[ipid] = new long [ncomp];
	
	for (k=0;k<ncomp;k++){
	  Mm->eigstrid[ipid][k] = id[k]-1;
	}
	
	ipid++;
      }//  end of the loop over the number of integration points on element
      
    }//  end of the loop over the number of elements 
  }//  end of the if (Mp->eigstrains==2)
  
  
  if (Mp->eigstrains==3){
    //  loop over the number of integration points in problem
    for (j=0;j<Mm->tnip;j++){
      ncomp = Mm->ip[j].ncompstr;
      Mm->eigstrains[j] = new double[ncomp];
      memset(Mm->eigstrains[j], 0, sizeof(*Mm->eigstrains[j])*ncomp);
      Mm->eigstresses[j] = new double [ncomp];
      memset(Mm->eigstresses[j], 0, sizeof(*Mm->eigstresses[j])*ncomp);
    }//  end of the loop over the number of integration points on element
    
  }//  end of the statement if (Mp->eigstrains==3){
  
  
  // ***************************************
  //  eigenstresses are defined on elements
  // ***************************************
  if (Mp->eigstrains==4){
    //  type of strain/stress state
    //  strastre = 1 - bar
    //  strastre = 10 - plane stress
    //  strastre = 30 - space strain/stress
    xfscanf (in,"%m", &strastrestate_kwdset, (int *)&ssst);
    
    switch (ssst){
    case bar:{
      xfscanf (in,"%ld",id);
      ncomp=1;
      break;
    }
    case planestress:{
      xfscanf (in,"%ld %ld %ld %ld",id,id+1,id+2,id+3);
      ncomp = 4;
      break;
    }
    case planestrain:{
      xfscanf (in,"%ld %ld %ld %ld",id, id+1, id+2, id+3);
      ncomp = 4;
      break;
    }
    case axisymm:{
      xfscanf(in,"%ld %ld %ld %ld", id, id+1, id+2, id+3);
      ncomp = 4;
      break;
    }
    case spacestress:{
      xfscanf (in,"%ld %ld %ld %ld %ld %ld",id,id+1,id+2,id+3,id+4,id+5);
      ncomp = 6;
      break;
    }
    default:{
      print_err("unknown type of strain-stress state", __FILE__, __LINE__, __func__);
      abort();
    }
    }//  end of the switch (strastre){
    
    for (i=0;i<Mm->tnip;i++){
      //  array of eigenstrains in integration points
      Mm->eigstrains[i] = new double [ncomp];
      memset(Mm->eigstrains[i], 0, sizeof(*Mm->eigstrains[i])*ncomp);
      //  array of eigenstresses in integration points
      Mm->eigstresses[i] = new double [ncomp];
      memset(Mm->eigstresses[i], 0, sizeof(*Mm->eigstresses[i])*ncomp);
      //  array of eigenstrain id
      Mm->eigstrid[i] = new long [ncomp];
      
      for (j=0;j<ncomp;j++){
	Mm->eigstrid[i][j]=id[j]-1;
      }
    }
  }//  end of the statement if (Mp->eigstrains==1){
  
  
  // ***************************************
  //  eigenstresses are defined on elements
  // ***************************************
  if (Mp->eigstrains==5){
    
    //  loop over the number of elements
    for (i=0;i<Mt->ne;i++){
      //  strain/stress state on element
      ssst = Mt->give_ssst (i,0);
      
      switch (ssst){
      case bar:{
	xfscanf (in,"%ld",id);
	ncomp = 1;
	break;
      }
      case planestress:{
	xfscanf (in,"%ld %ld %ld %ld", id, id+1, id+2, id+3);
	ncomp = 4;
	break;
      }
      case planestrain:{
	xfscanf (in,"%ld %ld %ld %ld", id, id+1, id+2, id+3);
	ncomp = 4;
	break;
      }
      case axisymm:{
        xfscanf(in,"%ld %ld %ld %ld", id, id+1, id+2, id+3);
        ncomp = 4;
        break;
      }
      case spacestress:{
	xfscanf (in,"%ld %ld %ld %ld %ld %ld",id,id+1,id+2,id+3,id+4,id+5);
	ncomp = 6;
	break;
      }
      default:{
	print_err("unknown type of strain-stress state", __FILE__, __LINE__, __func__);
        abort();
      }
      }//  end of the switch (strastre){
      
      //  number of integration points on element
      nip = Mt->give_tnip (i);
      //  number of the first integration point
      ipid=Mt->elements[i].ipp[0][0];
      
      //  loop over the number of integration points on element
      for (j=0;j<nip;j++){
      //  array of eigenstrains in integration points
        Mm->eigstrains[ipid] = new double [ncomp];
        memset(Mm->eigstrains[ipid], 0, sizeof(*Mm->eigstrains[ipid])*ncomp);
	Mm->eigstresses[ipid] = new double [ncomp];
        memset(Mm->eigstresses[ipid], 0, sizeof(*Mm->eigstresses[ipid])*ncomp);
	Mm->eigstrid[ipid] = new long [ncomp];
	
	for (k=0;k<ncomp;k++){
	  Mm->eigstrid[ipid][k] = id[k]-1;
	}
	
	ipid++;
      }//  end of the loop over the number of integration points on element
      
    }//  end of the loop over the number of elements 
  }//  end of the if (Mp->eigstrains==2)



  if (Mp->eigstrains>0){
    delete [] id;
  }
  
  // ************************************************
  //  general functions describing the eigenstrains
  // ************************************************
  if (Mp->eigstrains==1 || Mp->eigstrains==2 || Mp->eigstrains==4 || Mp->eigstrains==5){
    //  the number of general functions needed for eigenstrain description
    xfscanf (in,"%k%ld","number_of_gfunct_eigenstr",&ngfes);
    
    eigstrfun = new gfunct [ngfes];
    
    //  loop over the number of general functions
    for (i=0;i<ngfes;i++){
      eigstrfun[i].read (in);
    }//  end of the loop over the number of general functions
    
    //  computation of eigstrains
    eigstrain_computation (Mp->time);
    
  }//  end of the statement if (Mp->eigstrains==1 || Mp->eigstrains==2){

}



/**
   function prints eigenstrains, according to read_eigenstrains (XFILE *in)
   
   @param out - output file
   
   TKr, 07/02/2013
*/
void mechbclc::print_eigenstrains (FILE *out)
{
  long i;
  long ipid;
  
  //  Mp->eigstrains = 0 - eigenstrains are not defined
  //  Mp->eigstrains = 1 - eigenstrains are defined for the whole problem
  //  Mp->eigstrains = 2 - eigenstrains are defined on each element
  //  Mp->eigstrains = 3 - eigenstrains are defined outside of MEFEL (for the coupled problems)
  fprintf (out,"\n\n%ld",Mp->eigstrains);
  
  // *************************************************
  //  eigenstrains are defined for the whole problem
  // *************************************************
  if ((Mp->eigstrains==1) || (Mp->eigstrains==4)){
    //  type of strain/stress state
    //  strastre = 1 - bar
    //  strastre = 10 - plane stress
    //  strastre = 30 - space strain/stress
    fprintf (out,"\n\n%d",(int)ssst);
    
    switch (ssst){
    case bar:{
      fprintf (out,"\n  %ld",Mm->eigstrid[0][0]+1);
      break;
    }
    case planestress:{
      fprintf (out,"\n  %ld %ld %ld %ld",Mm->eigstrid[0][0]+1,Mm->eigstrid[0][1]+1,Mm->eigstrid[0][2]+1,Mm->eigstrid[0][3]+1);
      break;
    }
    case planestrain:{
      fprintf (out,"\n  %ld %ld %ld %ld",Mm->eigstrid[0][0]+1,Mm->eigstrid[0][1]+1,Mm->eigstrid[0][2]+1,Mm->eigstrid[0][3]+1);
      break;
    }
    case axisymm:{
      fprintf (out,"\n  %ld %ld %ld %ld",Mm->eigstrid[0][0]+1,Mm->eigstrid[0][1]+1,Mm->eigstrid[0][2]+1,Mm->eigstrid[0][3]+1);
      break;
    }
    case spacestress:{
      fprintf (out,"\n  %ld %ld %ld %ld %ld %ld",Mm->eigstrid[0][0]+1,Mm->eigstrid[0][1]+1,Mm->eigstrid[0][2]+1
	       ,Mm->eigstrid[0][3]+1,Mm->eigstrid[0][4]+1,Mm->eigstrid[0][5]+1);
      break;
    }
    default:{
      print_err("unknown type of strain-stress state", __FILE__, __LINE__, __func__);
    }
    }//  end of the switch (strastre){
  }//  end of the statement if (Mp->eigstrains==1){
  
  
  // ***************************************
  //  eigenstrains are defined on elements
  // ***************************************
  if ((Mp->eigstrains==2) || (Mp->eigstrains==5)){
    
    //  loop over the number of elements
    for (i=0;i<Mt->ne;i++){
      
      //  strain/stress state on element
      ssst = Mt->give_ssst (i,0);
      
      //  number of the first integration point
      ipid=Mt->elements[i].ipp[0][0];
      
      switch (ssst){
      case bar:{
	fprintf (out,"\n  %ld",Mm->eigstrid[ipid][0]+1);
	break;
      }
      case planestress:{
	fprintf (out,"\n  %ld %ld %ld %ld",Mm->eigstrid[ipid][0]+1,Mm->eigstrid[ipid][1]+1,Mm->eigstrid[ipid][2]+1,Mm->eigstrid[ipid][3]+1);
	break;
      }
      case planestrain:{
	fprintf (out,"\n  %ld %ld %ld %ld",Mm->eigstrid[ipid][0]+1,Mm->eigstrid[ipid][1]+1,Mm->eigstrid[ipid][2]+1,Mm->eigstrid[ipid][3]+1);
	break;
      }
      case axisymm:{
	fprintf (out,"\n  %ld %ld %ld %ld",Mm->eigstrid[ipid][0]+1,Mm->eigstrid[ipid][1]+1,Mm->eigstrid[ipid][2]+1,Mm->eigstrid[ipid][3]+1);
	break;
      }
      case spacestress:{
	fprintf (out,"\n  %ld %ld %ld %ld %ld %ld",Mm->eigstrid[ipid][0]+1,Mm->eigstrid[ipid][1]+1,Mm->eigstrid[ipid][2]+1
		 ,Mm->eigstrid[ipid][3]+1,Mm->eigstrid[ipid][4]+1,Mm->eigstrid[ipid][5]+1);
	break;
      }
      default:{
	print_err("unknown type of strain-stress state", __FILE__, __LINE__, __func__);
      }
      }//  end of the switch (strastre){
    }//  end of the loop over the number of elements 
  }//  end of the if (Mp->eigstrains==2)


  // ************************************************
  //  general functions describing the eigenstrains
  // ************************************************
  if (Mp->eigstrains==1 || Mp->eigstrains==2 || Mp->eigstrains==4 || Mp->eigstrains==5){
    //  the number of general functions needed for eigenstrain description
    fprintf (out,"\n\n  %ld",ngfes);
    
    //  loop over the number of general functions
    for (i=0;i<ngfes;i++){
      eigstrfun[i].print (out);
    }//  end of the loop over the number of general functions
    
  }//  end of the statement if (Mp->eigstrains==1 || Mp->eigstrains==2){

}



/**
  The function computes eigenstrains/eigenstresses at regular integration points.
 
  @param[in] time - actual time
   
  8. 11. 2012
*/
void mechbclc::eigstrain_computation (double time)
{
  long i, j, k, ii, jj, ipp, nip, ncomp, nb, gfid;
  vector coord(4); // three spatial coordinates  + time
  const char *namev[4];
  namev[0] = "x";
  namev[1] = "y";
  namev[2] = "z";
  namev[3] = "t";

  for (i=0;i<Gtm->ne;i++){
    nb = Mt->give_nb (i);
    for (ii=0;ii<nb;ii++){
      for (jj=0;jj<nb;jj++){
        ipp=Mt->elements[i].ipp[ii][jj];
        nip = Mt->give_nip (i,ii,jj);
        for (j=0;j<nip;j++){
          //  coordinates of integration points
          ipcoord (i,ipp,ii,jj,coord);
          //  the true number of strain components in the integration point
          ncomp = Mm->ip[ipp].ncompstr;
          coord(3) = time;
          //  loop over the number of eigenstrain components
          for (k=0;k<ncomp;k++){
            //  id of the general function
            gfid = Mm->eigstrid[ipp][k];
         
            if (Mp->eigstrains < 3)
              Mm->eigstrains[ipp][k] = eigstrfun[gfid].getval(coord, namev);
            if (Mp->eigstrains > 3)
              Mm->eigstresses[ipp][k] = eigstrfun[gfid].getval(coord, namev);
          }//  end of the loop over the number of eigenstrain components
          ipp++;
        } // end of loop over the number of ip in the given block
      }
    }
  }
}



/**
  The function computes eigenstrains/eigenstresses at auxiliary int. points.
   
  @param[in] n - number of auxiliary integration points (dimension of array ipm)
  @param[in] ipm - the array of integration point mapping objects
  @param[in] time - actual time for the evaulation of eigenstrain/eigenstress component functions

  The function does not return anything but it changes content of arrays aip_eigstrid, aip_eigstrains and aip_eigstresses.
   
  Created by TKo 10. 2021
*/
void mechbclc::aip_eigstrain_computation (long n, ipmap *ipm, double time)
{
  long i, k, app, ncomp, gfid;
  vector coord(4); // three spatial coordinates  + time
  const char *namev[4];
  namev[0] = "x";
  namev[1] = "y";
  namev[2] = "z";
  namev[3] = "t";

  if (Mm->aip == NULL) // there are no auxiliary int. points
    return;
  
  for (i=0; i<n; i++){
    if (ipm[i].ipp >= 0)  // direct mapping to regular integration point
      continue;
    //  coordinates of auxiliary integration points
    coord(0) = ipm[i].x;
    coord(1) = ipm[i].y;
    coord(2) = ipm[i].z;
    coord(3) = time;
    app = ipm[i].app;
    //  the true number of strain components in the integration point
    ncomp = Mm->aip[app].ncompstr;
    //  loop over the number of eigenstrain components
    for (k=0; k<ncomp; k++){
      //  id of the general function
      gfid = Mm->aip_eigstrid[app][k];
      if (Mp->eigstrains < 3)
        Mm->aip_eigstrains[app][k] = eigstrfun[gfid].getval(coord, namev);
      if (Mp->eigstrains > 3)
        Mm->aip_eigstresses[app][k] = eigstrfun[gfid].getval(coord, namev);
    }
  }
}



/**
   This method reads section with data about initial conditions in the nodes
   and if it is necessarry it allocates memory for mechmat ic array.
   
   @param in is pointer to the structure with opened file.
   
   @retval 0 - on succes
   @retval 1 - error reading number of nodes with initial condition
   @retval 2 - error reading node number
   
   Created by Tomas Koudelka,
*/
long mechbclc::readinic (XFILE *in)
{
  long i, n;
  
  //  the number of nodes where an initial condition or initial value of a material model is read
  xfscanf(in, "%k%ld", "num_ini_cond", &nico);

  if ((nico < 0) && (nico > Mt->nn))
  {
    print_err("invalid number of inital conditions", __FILE__, __LINE__, __func__);
    return 1;
  }
  if (nico == 0)
    return 0;
  
  //  if there is a nonzero number of nodes with initial condition or value,
  //  nn instances of the class inicd are created
  ico = new inicd[Mt->nn];
  for (i = 0L; i < nico; i++)
  {
    xfscanf(in, "%ld", &n);
    if ((n < 0) || (n > Mt->nn))
    {
      print_err("invalid node number", __FILE__, __LINE__, __func__);
      return 2;
    }

    //  reading of initial conditions or initial values
    ico[n-1].read(in);

    if ((ico[n-1].type & inicond) && (Mm->ic == NULL))
    {
      Mm->ic = new double *[Mm->tnip];
      memset(Mm->ic, 0, sizeof(*Mm->ic)*Mm->tnip);
    }
  }
  return 0;
}



/**
  This function prints section with data about initial conditions in the nodes
  and if it is necessarry it allocates memory for mechmat ic array.
  
  @param out is pointer to the structure with opened file.

  @retval 0 - on succes
  @retval 1 - error reading number of nodes with initial condition
  @retval 2 - error reading node number

  TKr, 11/02/2013 according to readinic (XFILE *in)
  Modified by TKo 7.2016
*/
long mechbclc::printinic (FILE *out)
{
  long i, nn=Mt->nn, aux=0L;

  fprintf(out, "\n ## initial conditions:");
  fprintf(out, "\n %ld\n", nico);

  for (i = 0L; i < nn; i++)
  {
    if ((ico[i].type == none) && (ico[i].nval == 0L))
    {
      fprintf(out, "\n %ld ",i+1);
      ico[aux].print(out);
      aux++;
    }
  }
  if (nico != aux)
    print_err("number of detected initial conditions at nodes (%ld) differs from the input (nico=%ld)", __FILE__, __LINE__, __func__, nico, aux);

  return 0;
}



/**
  Function computes initial values at integration points
  from initial nodal values.
   
  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void mechbclc::inicipval(void)
{
  long i, j, k, l, nne, nval, eid, aux;
  elemtype et;
  ivector enod;
  matrix  nodval;
  inictype *ictn;   // type of initial condition in the nodes


  for (i = 0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==1)
    {      
      nne = Mt->give_nne(i);
      reallocv(RSTCKIVEC(nne, enod));
      // find maximum number of initial values prescribed at nodes
      // initial displacement values are not taken into account
      Mt->give_elemnodes (i, enod);
      nval = ico[enod[0]].nval;
      if (ico[enod[0]].type & inidisp)
        nval -= Gtm->give_ndofn(enod[0]);
      if (ico[enod[0]].type & iniderdisp)
        nval -= Gtm->give_ndofn(enod[0]);
      for (j = 0; j < nne; j++)
      {
        aux = ico[enod[j]].nval;
        if (ico[enod[j]].type & inidisp)
          aux -= Gtm->give_ndofn(enod[j]);
        if (ico[enod[j]].type & iniderdisp)
          aux -= Gtm->give_ndofn(enod[j]);
        if (nval < aux)
          nval = aux;
      }
      reallocm (RSTCKMAT(nne, nval, nodval));
      fillm(0.0, nodval);
      ictn = new inictype[nne];
      for (j = 0; j < nne; j++)
      {
        ictn[j] = ico[enod[j]].type;
        l = 0L;
        if (ictn[j] & inidisp)  // skip initial nodal displacements
          l += Gtm->give_ndofn(enod[j]);
        if (ictn[j] & iniderdisp)  // skip initial time derivatives of nodal displacements
          l += Gtm->give_ndofn(enod[j]);
        for (k = 0; k < nval; k++)
        {
          nodval[j][k] = ico[enod[j]].val[k+l];
        }
      }
      eid = i;
      et = Mt->give_elem_type(eid);
      switch (et)
      {
	case bar2d:{              Bar2d->inicipval(eid, 0, 0, nodval, ictn);   break; }
	case barq2d:{             Barq2d->inicipval(eid, 0, 0, nodval, ictn);  break; }
	  //      case beam2d:{             Beam2d->inicipval(eid, 0, 0, nodval, ictn);  break; }
	  //      case beam3d:{             Beam3d->inicipval(eid, 0, 0, nodval, ictn);  break; }
	  //      case subsoilbeam:{        Sbeam->inicipval(eid, 0, 0, nodval, ictn);   break; }
	case spring_1:{           Spring->inicipval(eid, 0, 0, nodval, ictn);  break; }
	case spring_2:{           Spring->inicipval(eid, 0, 0, nodval, ictn);  break; }
	case spring_3:{           Spring->inicipval(eid, 0, 0, nodval, ictn);  break; }
	case spring_4:{           Spring->inicipval(eid, 0, 0, nodval, ictn);  break; }
	case spring_5:{           Spring->inicipval(eid, 0, 0, nodval, ictn);  break; }
	case spring_6:{           Spring->inicipval(eid, 0, 0, nodval, ictn);  break; }
	case planeelementlt:{     Pelt->inicipval(eid, 0, 0, nodval, ictn);    break; }
	case planeelementqt:{     Peqt->inicipval(eid, 0, 0, nodval, ictn);    break; }
	case planeelementrotlt:{  Perlt->inicipval(eid, 0, 0, nodval, ictn);   break; }
	case planeelementlq:{     Pelq->inicipval(eid, 0, 0, nodval, ictn);    break; }
	case planeelementqq:{     Peqq->inicipval(eid, 0, 0, nodval, ictn);    break; }
	case planeelementrotlq:{  Perlq->inicipval(eid, 0, 0, nodval, ictn);   break; }
	case planeelementsubqt:{  Pesqt->inicipval(eid, 0, 0, nodval, ictn);   break; }
	case cctel:{              Cct->inicipval(eid, 0, 0, nodval, ictn);     break; }
	case dktel:{              Dkt->inicipval(eid, 0, 0, nodval, ictn);     break; }
	case dstel:{              Dst->inicipval(eid, 0, 0, nodval, ictn);     break; }
	case q4plateel:{          Q4pl->inicipval(eid, 0, 0, nodval, ictn);    break; }
	  //      case subsoilplatetr:{     Spltr->inicipval(eid, 0, 0, nodval, ictn);   break; }
	  //      case subsoilplateq:{      Splq->inicipval(eid, 0, 0, nodval, ictn);    break; }
	case axisymmlq:{          Asymlq->inicipval(eid, 0, 0, nodval, ictn);  break; }
	case axisymmqq:{          Asymqq->inicipval(eid, 0, 0, nodval, ictn);  break; }
	case axisymmcq:{          Asymcq->inicipval(eid, 0, 0, nodval, ictn);  break; }
	case axisymmlt:{          Asymlt->inicipval(eid, 0, 0, nodval, ictn);  break; }
	case axisymmqt:{          Asymqt->inicipval(eid, 0, 0, nodval, ictn);  break; }
	case shelltrelem:{        Shtr->inicipval(eid, 0, 0, nodval, ictn);    break; }
	  //case shellqelem:{         Shq->inicipval(eid, 0, 0, nodval, ictn);     break; }
	case lineartet:{          Ltet->inicipval(eid, 0, 0, nodval, ictn);    break; }
	case quadrtet:{           Qtet->inicipval(eid, 0, 0, nodval, ictn);    break; }
	case linearhex:{          Lhex->inicipval(eid, 0, 0, nodval, ictn);    break; }
	case quadrhex:{           Qhex->inicipval(eid, 0, 0, nodval, ictn);    break; }
        case tetralatt:{                                                       break; }
	default:
        {
          print_err("unknown element type is required", __FILE__, __LINE__, __func__);
        }
      }
      delete [] ictn;
    }
  }
}



/**
  The function reads parameters for the calculation of prescribed initial displacement of selecetd nodes.
  The prescribed initial displacements are calculated by the rotation of selected nodes about an axis
  given by two points and angle given by directional cosine. The function is used in the growing mechanical problems
  in the case that the initial displacements of new active parts of structure should be calculated. Fo each new independent 
  part of structure that is being activated, one record consisting of time, axis definition, selected nodes, and directional 
  cosines may be specified. The selected nodes are expected to be at the surface of the new active part of structure and 
  the vector of their prescribed initial displacements will obtained by the difference between the original position and 
  position obtained by rotation of selected nodes about the given axis.

  @param in - pointer to the opened input XFILE structure

  @return The function does not return anything but stores paramters of the rotation in the mechbclc structure.

  Created by Tomas Koudelka 07.2016
*/
void mechbclc::read_rotinidispl(XFILE *in)
{
  long i;

  xfscanf(in, "%k %ld", "num_axis", &naxis);
  if ((naxis < 0) || (naxis > Gtm->nn))
    print_err("invalid number of axis, %ld is out of range <0; %ld>", __FILE__, __LINE__, __func__, naxis, Gtm->nn);

  arotrec = new axisrotrec[naxis];

  for (i=0; i<naxis; i++)
    arotrec[i].read(in);
}



/**
  The function clears nodal DOF numbers at nodes with prescribed initial displacements.
  The function is used at simulations of gradual construction process.

  @param ifn - array of flags of interface nodes between new active part and old part that have 
               to be excluded even though that there are defined initial displacements 
               (these initial displacement values are intended for the usage in another time step
                when these active intefrace nodes are just activated and they are not on the interface).
  @param time - actual time

  @return The function does not return anything but changes the content of DOF number array of particular nodes.
               
  Created by Tomas Koudelka, 07.2016
*/
void mechbclc::clear_rotinidispl(long *ifn, double time)
{
  long i, j, k, min_ndofn;

  if (dlc == NULL)
  {
    print_err("no time dependent load case has been defined,\n"
              " the function is intended for the growing time dependent problem only",
               __FILE__, __LINE__, __func__);
    abort();
  }
  // Detect the number of prescribed initial displacements at the current time
  for (i=0; i<naxis; i++)
  {
    if ((time < arotrec[i].stime) || (time > arotrec[i].etime))
      continue; // skip axis defined for another time steps

    for(j=0; j<Gtm->nn; j++)
    {
      if ((Gtm->lnso[j] == 0L) || ifn[j])
        continue;  // skip inactive or interface nodes

      // check for selected nodes
      if (arotrec[i].seln.presence_id(j))
      {
        min_ndofn = min2(Gtm->give_ndofn(j), 3);
        for (k=0; k<min_ndofn; k++)
        {
          if (arotrec[i].dap[k] == yes)
            Gtm->gnodes[j].clear_dof(k);
        }
      }
    }
  }
}



/**
  The function stores prescribed initial displacements by rotations at selected nodes.
  The function is used at simulations of gradual construction process for the calculation of stresses 
  due to prescribed initial displacements. 

  @param ifn - array of flags of interface nodes between new active part and old part that have 
               to be excluded even though that there are defined initial displacements 
               (these initial displacement values are intended for the usage in another time step
                when these active intefrace nodes are just activated and they are not on the interface).
  @param time - actual time
  @param   r - displacement %vector where the values of the initial displacements are localized.

  @return The function does not return anything but changes the content of DOF number array of particular nodes.
               
  Created by Tomas Koudelka, 07.2016
*/
void mechbclc::store_rotinidispl(long *ifn, double time, double *r)
{
  long i, j, k, cn, min_ndofn;
  vector coord(ASTCKVEC(3)), ur(ASTCKVEC(3));

  if (dlc == NULL)
  {
    print_err("no time dependent load case has been defined,\n"
              " the function is intended for the growing time dependent problem only",
               __FILE__, __LINE__, __func__);
    abort();
  }
  // Detect the number of prescribed initial displacements at the current time
  for (i=0; i<naxis; i++)
  {
    if ((time < arotrec[i].stime) || (time > arotrec[i].etime))
      continue; // skip axis defined for another time steps

    for(j=0; j<Gtm->nn; j++)
    {
      if ((Gtm->lnso[j] == 0L) || ifn[j])
        continue;  // skip inactive or interface nodes

      // check for selected nodes
      if (arotrec[i].seln.presence_id(j))
      {
        min_ndofn = min2(Gtm->give_ndofn(j), 3);
        Gtm->give_node_coord(j, coord);
        arotrec[i].compute_rotcoord(coord, ur);
        for (k=0; k<min_ndofn; k++)
        {
          cn = Gtm->give_dof(j, k);
          if ((cn > 0) && (arotrec[i].dap[k] == yes))
            r[cn-1] = ur[k];
        }
      }
    }
  }
}



/**
  The function computes the number of different DOFs with the prescribed initial displacements by rotation.
  These prescribed initial displacements are applied in course of initial displacement calculation of new 
  activated part of structure at simulations of gradual construction process.  

  @param ifn - array of flags of interface nodes between new active part and old part that have 
               to be excluded even though that there are defined initial displacements 
               (these initial displacement values are intended for the usage in another time step
                when these active intefrace nodes are just activated and they are not on the interface).
  @param time - actual time

  @return The function returns the number of DOFs with prescribed initial displacements by rotation about the given axis.
 
  Created by Tomas Koudelka, 08.2016
*/
long mechbclc::num_dofs_rotinidispl(long *ifn, double time)
{
  long i, j, k, min_ndofn, nrpid = 0L;

  // Detect the number of prescribed initial displacements at the current time
  for (i=0; i<naxis; i++)
  {
    if ((time < arotrec[i].stime) || (time > arotrec[i].etime))
      continue; // skip axis defined for another time steps

    for(j=0; j<Gtm->nn; j++)
    {
      if ((Gtm->lnso[j] == 0L) || ifn[j])
        continue;  // skip inactive or interface nodes

      // check for selected nodes
      if (arotrec[i].seln.presence_id(j))
      {
        min_ndofn = min2(Gtm->give_ndofn(j), 3);
        for (k=0; k<min_ndofn; k++)
        {
          if (arotrec[i].dap[k] == yes)
            nrpid++;
        }
      }
    }
  }
  return nrpid;
}



/**
  The function actualizes vector of displacements of axis definition nodes
  for all active axes. It is used in the prescribed initial displacements 
  by rotations about the given axis.

  @param lcid - load case id
  @param time - actual time

  @return The function does not return anything but it may change the content of
          particular components of arotrec array.

  Created by Tomas Koudelka, 08.2016
*/
void mechbclc::actualize_displ_def_node(long lcid, double time)
{
  long i;

  for (i=0; i<naxis; i++)
    arotrec[i].store_displ_def_node(lcid, time);

  return;
}



/**
  The function applies prescribed initial displacements at selected nodes of the new activated part of structure
  obtained by the rotation of the selected nodes about the given axis. These prescribed initial displacements are 
  applied in course of initial displacement calculation of new activated part of structure at simulations of 
  gradual construction process.

  @param lcid - time dependent load case (dloadcase) id
  @param ifn - array of flags of interface nodes between new active part and old part that have 
               to be excluded even though that there are defined initial displacements 
               (these initial displacement values are intended for the usage in another time step
                when these active intefrace nodes are just activated and they are not on the interface).
  @param time - actual time

  @return The function does not return anything but changes the content of DOF number array of particular nodes.
               
  Created by Tomas Koudelka, 07.2016
*/
void mechbclc::apply_rotinidispl(long lcid, long *ifn, double time)
{
  long i, j, k, max_npd, ipid;
  vector coord(ASTCKVEC(3)), ur(ASTCKVEC(3));

  if (dlc == NULL)
  {
    print_err("no time dependent load case has been defined,\n"
              " the function is intended for the growing time dependent problem only",
               __FILE__, __LINE__, __func__);
    abort();
  }

  if (dlc[lcid].nrpid == 0L)  // there are no prescribed initial displacements defined in the problem
    return;

  // maximum of prescribed displacements in the given dloadcase
  max_npd = dlc[lcid].max_npd;

  
  // Setup of array pid with prescribed intial displacements of the new active part of structure
  ipid = dlc[lcid].npid;
  for (i=0; i<naxis; i++)
  {
    if ((time < arotrec[i].stime) || (time > arotrec[i].etime))
      continue; // skip axis defined for another time steps

    for(j=0; j<Gtm->nn; j++)
    {
      if ((Gtm->lnso[j] == 0L) || ifn[j])
        continue;  // skip inactive or interface nodes

      // check for selected nodes
      if (arotrec[i].seln.presence_id(j))
      {
        Gtm->give_node_coord(j, coord);
        arotrec[i].compute_rotcoord(coord, ur);
        for (k=0; k<min2(Gtm->give_ndofn(j), 3); k++)
        {
          if (arotrec[i].dap[k] == yes)
          {
            dlc[lcid].pid[ipid] = ur[k];
            ipid++;
            Gtm->gnodes[j].save_dof(k, -(max_npd+ipid));
          }
        }
      }
    }
  }
}



/**
  The function clears nodal DOF numbers at nodes with prescribed initial displacements.
  The function is used at simulations of gradual construction process.

  @param ifn - array of flags of interface nodes between new active part and old part that have 
               to be excluded even though that there are defined initial displacements 
               (these initial displacement values are intended for the usage in another time step
                when these active intefrace nodes are just activated and they are not on the interface).

  @return The function does not return anything but changes the content of DOF number array of particular nodes.
               
  Created by Tomas Koudelka, 07.2016
*/
void mechbclc::clear_inidispl(long *ifn)
{
  long i;

  if (nico == 0)  // there are no initial conditions
    return;

  for(i=0; i<Gtm->nn; i++)
  {
    if ((Gtm->lnso[i] == 0) || (ifn[i] == 1))
      continue;  // skip inactive or interface nodes

    if ((ico[i].type & inidisp)) // initial displacement for all DOFs are specified
    {
      if ((ico[i].type & inidisp_x) || (ico[i].type & inidisp_y) || (ico[i].type & inidisp_z))
      {
        if ((ico[i].type & inidisp_x))
          Gtm->gnodes[i].clear_dof(0);
        if ((ico[i].type & inidisp_y))
          Gtm->gnodes[i].clear_dof(1);
        if ((ico[i].type & inidisp_z))
          Gtm->gnodes[i].clear_dof(2);
      }
      else
        Gtm->gnodes[i].clear_dof();
    }
  }
}



/**
  The function stores prescribed initial displacements at selected nodes to the displacement %vector.
  The function is used at simulations of gradual construction process for the calculation of stresses 
  due to prescribed initial displacements.

  @param ifn - array of flags of interface nodes between new active part and old part that have 
               to be excluded even though that there are defined initial displacements 
               (these initial displacement values are intended for the usage in another time step
                when these active intefrace nodes are just activated and they are not on the interface).
  @param   r - displacement %vector where the values of the initial displacements are localized.

  @return The function does not return anything but changes the content of array r for the particular nodes.
               
  Created by Tomas Koudelka, 08.2016
*/
void mechbclc::store_inidispl(long *ifn, double *r)
{
  long i, j, cn, ndofn;

  if (nico == 0)  // there are no initial conditions
    return;

  for(i=0; i<Gtm->nn; i++)
  {
    if ((Gtm->lnso[i] == 0) || (ifn[i] == 1))
      continue;  // skip inactive or interface nodes

    if ((ico[i].type & inidisp)) // initial displacement for all DOFs are specified
    {
      if ((ico[i].type & inidisp_x) || (ico[i].type & inidisp_y) || (ico[i].type & inidisp_z))
      {
        if ((ico[i].type & inidisp_x))
          r[Gtm->give_dof(i, 0)] = ico[i].val[0];
        if ((ico[i].type & inidisp_y))
          r[Gtm->give_dof(i, 1)] = ico[i].val[1];
        if ((ico[i].type & inidisp_z))
          r[Gtm->give_dof(i, 2)] = ico[i].val[2];
      }
      else
      {
        ndofn = Gtm->give_ndofn(i);
        for(j=0; j<ndofn; j++)
        {
          cn = Gtm->give_dof(i, j);
          if (cn > 0)
            r[cn-1] = ico[i].val[j];
        }
      }
    }
  }
}



/**
  The function computes the number of different DOFs with the prescribed initial displacements given at initial conditions.
  These prescribed initial displacements are applied in course of initial displacement calculation of new 
  activated part of structure at simulations of gradual construction process.  

  @param ifn - array of flags of interface nodes between new active part and old part that have 
               to be excluded even though that there are defined initial displacements 
               (these initial displacement values are intended for the usage in another time step
                when these active intefrace nodes are just activated and they are not on the interface).

  @return The function returns the number of DOFs with prescribed initial displacements given at initial conditions.
 
  Created by Tomas Koudelka, 08.2016
*/
long mechbclc::num_dofs_inidispl(long *ifn)
{
  long i, npid = 0L;

  if (nico == 0)  // there are no initial conditions
    return npid;

  // detect the number of prescribed initial displacements
  for(i=0; i<Gtm->nn; i++)
  {
    if ((Gtm->lnso[i] == 0L) || ifn[i])
      continue;  // skip inactive or interface nodes

    if ((ico[i].type & inidisp)) // initial displacement for all DOFs are specified
    {
      if ((ico[i].type & inidisp_x) || (ico[i].type & inidisp_y) || (ico[i].type & inidisp_z))
      {
        if ((ico[i].type & inidisp_x))
          npid++;
        if ((ico[i].type & inidisp_y))
          npid++;
        if ((ico[i].type & inidisp_z))
          npid++;
      }
      else
        npid += Gtm->give_ndofn(i);
    }
  }
  return npid;
}



/**
  The function applies prescribed initial displacements at nodes of the new activated part of structure.
  These prescribed initial displacements are applied in course of initial displacement calculation of 
  new activated part of structure at simulations of gradual construction process.

  @param lcid - time dependent load case (dloadcase) id
  @param ifn - array of flags of interface nodes between new active part and old part that have 
               to be excluded even though that there are defined initial displacements 
               (these initial displacement values are intended for the usage in another time step
                when these active intefrace nodes are just activated and they are not on the interface).

  @return The function does not return anything but changes the content of DOF number array of particular nodes.
               
  Created by Tomas Koudelka, 07.2016
*/
void mechbclc::apply_inidispl(long lcid, long *ifn)
{
  long i, j, ipid, max_npd, ndofn;

  if (dlc == NULL)
  {
    print_err("no time dependent load case has been defined,\n"
              " the function is intended for the growing time dependent problem only",
               __FILE__, __LINE__, __func__);
    abort();
  }

  if (dlc[lcid].npid == 0L)  // there are no initial displacements defined in the problem but some other initial values do
    return;

  // maximum of prescribed displacements in the given dloadcase
  max_npd = dlc[lcid].max_npd;


  ipid = 0L;
  for(i=0; i<Gtm->nn; i++)
  {
    if ((Gtm->lnso[i] == 0L) || ifn[i])
      continue;  // skip inactive or interface nodes

    if ((ico[i].type & inidisp)) // initial displacement for all DOFs are specified
    {
      if ((ico[i].type & inidisp_x) || (ico[i].type & inidisp_y) || (ico[i].type & inidisp_z))
      { // only selected directions of intial displacements specified are taken into account
        if (ico[i].type & inidisp_x)
        {
          dlc[lcid].pid[ipid] = ico[i].val[0];    
          ipid++;
          Gtm->gnodes[i].save_dof(0, -(max_npd+ipid));
        }
        if ((ico[i].type & inidisp_y))
        {
          dlc[lcid].pid[ipid] = ico[i].val[1];    
          ipid++;
          Gtm->gnodes[i].save_dof(1, -(max_npd+ipid));
        }
        if ((ico[i].type & inidisp_z))
        {
          dlc[lcid].pid[ipid] = ico[i].val[2];    
          ipid++;
          Gtm->gnodes[i].save_dof(2, -(max_npd+ipid));
        }
      }
      else // all initial displacements specified are taken into account
      {
        ndofn = Gtm->give_ndofn(i);
        for (j=0L; j<ndofn; j++)
        {
          dlc[lcid].pid[ipid] = ico[i].val[j];
          ipid++;
          Gtm->gnodes[i].save_dof(j, -(max_npd+ipid));
        }    
      }
    }
  }
}



/**
  Function allocates array containing sums of components 
  in particular directions of the nodal force %vector.
   
  @return The function does not return anything.

  Created by JK, 6.10.2003
*/
void mechbclc::alloc_sumcomp ()
{
  ncsum = Mt->give_ndofn (0);
  sumcomp = new double [ncsum];
  reactsumcomp = new double [ncsum];
  pd_reactsumcomp = new double [ncsum];
  memset(sumcomp, 0, sizeof(*sumcomp)*ncsum);
  memset(reactsumcomp, 0, sizeof(*reactsumcomp)*ncsum);
  memset(pd_reactsumcomp, 0, sizeof(*pd_reactsumcomp)*ncsum);
}



/**
  Function computes sums of components in particular directions 
  of the nodal force %vector.
   
  @param rhs - %vector of right hand side (nodal forces)
  
  Created by JK, 6.10.2003
*/
void mechbclc::comp_sum (double *rhs)
{
  long i,j,ndofn,cn;
  long *aux;
  
  aux=new long [Ndofm];
  for (i=0;i<Ndofm;i++){
    aux[i]=0;
  }
  
  nullv (sumcomp,ncsum);
  
  for (i=0;i<Mt->nn;i++){
    ndofn=Gtm->gnodes[i].ndofn;
    for (j=0;j<ndofn;j++){
      cn=Mt->give_dof (i,j);
      if (cn>0){
	if (aux[cn-1]==0){
	  sumcomp[j]+=rhs[cn-1];
	  aux[cn-1]=1;
	}
      }
    }
  }

  delete [] aux;
}



/**
  The function computes sums of components in particular directions of the reactions at supports.
  Particular sums of components are stored in reactsumcomp array.

  @return The function does not return anything.

  Created by Tomas Koudelka,  12.2009
*/
void mechbclc::comp_sum_react()
{
  long i,j,ndofn,cn;
  
  nullv (reactsumcomp,ncsum);
  
  for (i=0;i<Mt->nn;i++){
    ndofn=Gtm->gnodes[i].ndofn;
    for (j=0;j<ndofn;j++){
      cn=Mt->give_dof(i,j);
      if (cn==0)
       reactsumcomp[j]+=Mt->nodes[i].r[j];
    }
  }
}



/**
  The function computes sums of components in particular directions of the reactions caused by prescribed displacements.
  Particular sums of components are stored in pd_reactsumcomp array.
   
  @return The function does not return anything.

  Created by Tomas Koudelka,  12.2009
*/
void mechbclc::comp_sum_pdreact()
{
  long i,j,ndofn,cn;
  
  nullv (sumcomp,ncsum);
  
  for (i=0;i<Mt->nn;i++){
    ndofn=Gtm->gnodes[i].ndofn;
    for (j=0;j<ndofn;j++){
      cn=Mt->give_dof (i,j);
      if (cn<0)
        pd_reactsumcomp[j]+=Mt->nodes[i].r[j];
    }
  }
}



/**
  Function returns sums of components in particular directions of the nodal force %vector.
   
  Created by JK, 6.10.2003
*/
void mechbclc::give_comp_sum (double *sum)
{
  long i;
  for (i=0;i<ncsum;i++)
    sum[i]=sumcomp[i];
}



/**
  The function computes temperature strains at auxiliary integration points
  for the given load case lcid. It is intended for the use in coupled problems esspecially.
  The load case must be type of time dependent load case, dloadcase instance must be defined.

  @param lcid[in] - load case id
  @param n[in]    - the number of required auxiliary points in the mapping array ipm
  @param ipm[in]  - integration point mapping array, 
                    ipm[i].ipp < 0 => auxiliary integration point must be used
                    ipm[i].ipp >= 0 => direct mapping to the regulal integration point, 
                                       no computation of strains is performed, the strains are assumed 
                                       to be computed at the main solution procedure of the problem

  @return The function does not return anything but it updates content of mechmat::aip_tempstains array.

  Created by Tomas Koudelka, 18.6.2018
*/
void mechbclc::aip_temperstrains(long lcid, long n, ipmap *ipm)
{
  if (dlc)
    // compute cumulative temperature strains from all subloadcases
    dlc[lcid].aip_temperstrains(lcid, n, ipm);
  else
  {
    print_err("invalid load case type is required for the computation"
              "  of temperature strains at auxiliary int. points", __FILE__, __LINE__, __func__);
    abort();
  }
}



/**
  The function returns the number of prescribed macro-stress components in the given load case. 
  It is itended for the homogenization problems.

  @param[in] lcid - load case id

  @return The number of prescribed macro-stress components.
*/
long mechbclc::give_num_mstress_comp(long lcid)
{
  long ret = 0;
  if (lc){
    ret = lc[lcid].give_num_mstress_comp();
    return ret;
  }
  ret = dlc[lcid].give_num_mstress_comp();
  return ret;
}



/**
  The function returns the number of prescribed macro-strain components in the given load case. 
  It is itended for the homogenization problems.

  @param[in] lcid - load case id

  @return The number of prescribed macro-strain components.
*/
long mechbclc::give_num_mstrain_comp(long lcid)
{
  long ret = 0;
  if (lc){
    ret = lc[lcid].give_num_mstrain_comp();
    return ret;
  }
  if (dlc){
    ret = dlc[lcid].give_num_mstrain_comp();
    return ret;
  }
  return ret;
}



/**
  The function returns pointer to the array of types of prescribed macro-value
  components. It is itended for the homogenization problems (Mp->homog > 2).

  @param[in] lcid - load case id

  @return The pointer to the array of prescribed macro-value types.
*/
strastre* mechbclc::give_mstrastre(long lcid)
{
  strastre *ret = NULL;
  if (lc){
    ret = lc[lcid].give_mstrastre();
    return ret;
  }
  if (dlc){
    ret = dlc[lcid].give_mstrastre();
    return ret;
  }

  return ret;
}



/**
  The function returns pointer to the array of code numbers of macro-stress 
  components. It is itended for the homogenization problems  (Mp->homog > 2).

  @param[in] lcid - load case id

  @return The pointer to the array of macro-stress component code (DOF) numbers,
          dimension of the array is the maximum total number of stress/strain components of used elements.
          Particular components of the array contains either positive nonzero value which represents the code (DOF) number
          or zero value for components where no macro-stress value has been prescribed.
*/
long* mechbclc::give_mstress_cn(long lcid)
{
  long *ret = NULL;
  if (lc){
    ret = lc[lcid].give_mstress_cn();
    return ret;
  }
  ret = dlc[lcid].give_mstress_cn();
  return ret;
}



/**
  The function %vector of actual prescribed macro-strain components in the given load case. 
  It is itended for the homogenization problems.

  @param[in] lcid - load case id
  @param[in] time - actual time or load coefficient
  @param[out] mstra - the resulting %vector of macro-strain components, it must be allocated to hold ncompmv components

  @return The function returns actual macro-strain components in the argument mstra.
*/
void mechbclc::give_mstrains(long lcid, double time, vector &mstra)
{
  nullv(mstra);
  
  if (lc){
    lc[lcid].give_mstrains(mstra);
    if (Mp->tprob == mat_nonlinear_statics){
      vector aux(ASTCKVEC(mstra.n));
      lc[lcid+1].give_mstrains(aux);
      addmultv(aux, mstra, time, mstra);
    }
    return;
  }
  dlc[lcid].give_mstrains(time, mstra);
}



/**
  The function %vector of actual prescribed macro-stress components in the given load case. 
  It is itended for the homogenization problems.

  @param[in] lcid - load case id
  @param[in] time - actual time or load coefficient
  @param[out] mstre - the resulting %vector of macro-stress components, it must be allocated to hold ncompmv components

  @return The function returns actual macro-stress components in the argument mstre.
*/
void mechbclc::give_mstresses(long lcid, double time, vector &mstre)
{
  nullv(mstre);
  
  if (lc){
    lc[lcid].give_mstresses(mstre);
    if (Mp->tprob == mat_nonlinear_statics){
      vector aux(ASTCKVEC(mstre.n));
      lc[lcid+1].give_mstresses(aux);
      addmultv(aux, mstre, time, mstre);
    }
    return;
  }
  dlc[lcid].give_mstresses(time, mstre);
}



