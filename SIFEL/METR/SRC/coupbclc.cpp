#include "coupbclc.h"
#include "globalc.h"
#include "aliasc.h"
//#include "loadcase.h"
//#include "dloadcase.h"
#include "inicd.h"
#include "gfunct.h"

coupbclc::coupbclc (void)
{
  nlc=2;

  //  the number of nodes where inital values of material models are defined
  nico = 0L;
  //  nodal inital values of material models
  ico = NULL;

  //  the number of general functions describing the eigenstrains
  ngfes=0L;
  //  array of general functions describing the eigenstrains
  eigstrfun = NULL;
  strastre = (strastrestate)0;
}

coupbclc::~coupbclc (void)
{
  if (ico!=NULL)
    delete [] ico;
  if (eigstrfun!=NULL)
    delete [] eigstrfun;
}

void coupbclc::read (FILE */*in*/)
{
}

/**
  This method reads section with data about initial conditions in the nodes
  and if it is necessarry it allocates memory for mechmat ic array.
  
  @param in is pointer to the structure with opened file.

  @retval 0 - on succes
  @retval 1 - error reading number of nodes with initial condition
  @retval 2 - error reading node number

  Created by Tomas Koudelka, 11.4.2019
*/
long coupbclc::readinic (XFILE *in)
{
  long i, n;
  
  //  the number of nodes where inital values of material models are defined
  xfscanf(in, "%k%ld", "num_ini_cond", &nico);
  
  if ((nico < 0) && (nico > Ct->nn)){
    print_err("invalid number of inital conditions", __FILE__, __LINE__, __func__);
    return 1;
  }
  if (nico == 0)
    return 0;

  ico = new inicd[Ct->nn];
  for (i = 0L; i < nico; i++){
    xfscanf(in, "%ld", &n);
    if ((n < 0) || (n > Ct->nn)){
      print_err("invalid node number", __FILE__, __LINE__, __func__);
      return 2;
    }

    ico[n-1].read(in);

    if ((ico[n-1].type & inicond) && (Cm->ic == NULL)){
      Cm->ic = new double *[Cm->tnip];
      memset(Cm->ic, 0, sizeof(*Cm->ic)*Cm->tnip);
    }
  }
  return 0;
}


/**
   function reads eigenstrains
   
   @param in - input file
   
   JK, 11.4.2019
*/
void coupbclc::read_eigenstrains (XFILE *in)
{
  long i,j,k;
  long ncomp,nip,ipid;
  long *id = NULL;
  
  //  Cp->eigstrains = 0 - eigenstrains are not defined
  //  Cp->eigstrains = 1 - eigenstrains are defined for the whole problem
  //  Cp->eigstrains = 2 - eigenstrains are defined on each element
  //  Cp->eigstrains = 3 - eigenstrains are defined outside of MEFEL (for the coupled problems)
  //  Cp->eigstrains = 4 - eigenstresses are defined for the whole problem
  //  Cp->eigstrains = 5 - eigenstresses are defined on each element
  xfscanf (in,"%k%ld","eigenstrains",&(Cp->eigstrains));
  if (Cp->eigstrains<0 || Cp->eigstrains>5){
    print_err("wrong definition of eigenstrains in file %s (line %ld)", __FILE__, __LINE__, __func__, in->fname, in->line);
    abort ();
  }
  
  if (Mesprc == 1){
    if (Cp->eigstrains==0)
      fprintf(stdout, "\n eigenstrains are not defined");
    if (Cp->eigstrains==1)
      fprintf(stdout, "\n eigenstrains are defined for the whole problem by a single set of values");
    if (Cp->eigstrains==2)
      fprintf(stdout, "\n eigenstrains are defined independently for each element");
    if (Cp->eigstrains==3)
      fprintf(stdout, "\n eigenstrains are defined independently in TRFEL");
    if (Cp->eigstrains==4)
      fprintf(stdout, "\n eigenstresses are defined for the whole problem by a single set of values");
    if (Cp->eigstrains==5)
      fprintf(stdout, "\n eigenstresses are defined independently for each element");
  }
  
  
  if (Cp->eigstrains>0){
    //  eigenstrains are defined
    
    //  array of eigenstrain id
    if (Cm->eigstrid==NULL)
      Cm->eigstrid = new long* [Cm->tnip];
    //  array of eigenstrains in integration points
    if (Cm->eigstrains==NULL)
      Cm->eigstrains = new double* [Cm->tnip];
    //  array of eigenstresses in integration points
    if (Cm->eigstresses==NULL)
      Cm->eigstresses = new double* [Cm->tnip];
    
    id = new long [6];
  }
  
  
  // *************************************************
  //  eigenstrains are defined for the whole problem
  // *************************************************
  if (Cp->eigstrains==1){
    //  type of strain/stress state
    //  strastre = 1 - bar
    //  strastre = 10 - plane stress
    //  strastre = 30 - space strain/stress
    xfscanf (in,"%m", &strastrestate_kwdset, (int *)&strastre);
    
    switch (strastre){
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
    }
    }//  end of the switch (strastre){
    
    for (i=0;i<Cm->tnip;i++){
      //  array of eigenstrains in integration points
      Cm->eigstrains[i] = new double [ncomp];
      memset(Cm->eigstrains[i], 0, sizeof(*Cm->eigstrains[i])*ncomp);
      //  array of eigenstresses in integration points
      Cm->eigstresses[i] = new double [ncomp];
      memset(Cm->eigstresses[i], 0, sizeof(*Cm->eigstresses[i])*ncomp);
      //  array of eigenstrain id
      Cm->eigstrid[i] = new long [ncomp];
      
      for (j=0;j<ncomp;j++){
        Cm->eigstrid[i][j]=id[j]-1;
      }
    }
    
  }//  end of the statement if (Mp->eigstrains==1){
  
  
  // ***************************************
  //  eigenstrains are defined on elements
  // ***************************************
  if (Cp->eigstrains==2){
    
    //  loop over the number of elements
    for (i=0;i<Ct->ne;i++){
      //  strain/stress state on element
      strastre = Ct->give_ssst (i,0);
      
      switch (strastre){
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
      }
      }//  end of the switch (strastre){
      
      //  number of integration points on element
      nip = Ct->give_nip (i);
      //  number of the first integration point
      ipid=Ct->elements[i].ipp;
      
      //  loop over the number of integration points on element
      for (j=0;j<nip;j++){
	Cm->eigstrains[ipid] = new double [ncomp];
        memset(Cm->eigstrains[ipid], 0, sizeof(*Cm->eigstrains[ipid])*ncomp);
	Cm->eigstresses[ipid] = new double [ncomp];
        memset(Cm->eigstresses[ipid], 0, sizeof(*Cm->eigstresses[ipid])*ncomp);
	Cm->eigstrid[ipid] = new long [ncomp];
	
	for (k=0;k<ncomp;k++){
	  Cm->eigstrid[ipid][k] = id[k]-1;
	}
	
	ipid++;
      }//  end of the loop over the number of integration points on element
      
    }//  end of the loop over the number of elements 
  }//  end of the if (Mp->eigstrains==2)
  
  
  if (Cp->eigstrains==3){
    //  loop over the number of integration points in problem
    for (j=0;j<Cm->tnip;j++){
      ncomp = Cm->ip[j].ncompstr;
      Cm->eigstrains[j] = new double[ncomp];
      memset(Cm->eigstrains[j], 0, sizeof(*Cm->eigstrains[j])*ncomp);
      Cm->eigstresses[j] = new double [ncomp];
      memset(Cm->eigstresses[j], 0, sizeof(*Cm->eigstresses[j])*ncomp);
    }//  end of the loop over the number of integration points on element
    
  }//  end of the statement if (Mp->eigstrains==3){
  
  
  // ***************************************
  //  eigenstresses are defined on elements
  // ***************************************
  if (Cp->eigstrains==4){
    //  type of strain/stress state
    //  strastre = 1 - bar
    //  strastre = 10 - plane stress
    //  strastre = 30 - space strain/stress
    xfscanf (in,"%m", &strastrestate_kwdset, (int *)&strastre);
    
    switch (strastre){
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
    }
    }//  end of the switch (strastre){
    
    for (i=0;i<Cm->tnip;i++){
      //  array of eigenstrains in integration points
      Cm->eigstrains[i] = new double [ncomp];
      memset(Cm->eigstrains[i], 0, sizeof(*Cm->eigstrains[i])*ncomp);
      //  array of eigenstresses in integration points
      Cm->eigstresses[i] = new double [ncomp];
      memset(Cm->eigstresses[i], 0, sizeof(*Cm->eigstresses[i])*ncomp);
      //  array of eigenstrain id
      Cm->eigstrid[i] = new long [ncomp];
      
      for (j=0;j<ncomp;j++){
	Cm->eigstrid[i][j]=id[j]-1;
      }
    }
    
  }//  end of the statement if (Mp->eigstrains==1){
  
  
  // ***************************************
  //  eigenstresses are defined on elements
  // ***************************************
  if (Cp->eigstrains==5){
    
    //  loop over the number of elements
    for (i=0;i<Ct->ne;i++){
      //  strain/stress state on element
      strastre = Ct->give_ssst (i,0);
      
      switch (strastre){
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
      }
      }//  end of the switch (strastre){
      
      //  number of integration points on element
      nip = Ct->give_nip (i);
      //  number of the first integration point
      ipid=Ct->elements[i].ipp;
      
      //  loop over the number of integration points on element
      for (j=0;j<nip;j++){
      //  array of eigenstrains in integration points
        Cm->eigstrains[ipid] = new double [ncomp];
        memset(Cm->eigstrains[ipid], 0, sizeof(*Cm->eigstrains[ipid])*ncomp);
	Cm->eigstresses[ipid] = new double [ncomp];
        memset(Cm->eigstresses[ipid], 0, sizeof(*Cm->eigstresses[ipid])*ncomp);
	Cm->eigstrid[ipid] = new long [ncomp];
	
	for (k=0;k<ncomp;k++){
	  Cm->eigstrid[ipid][k] = id[k]-1;
	}
	
	ipid++;
      }//  end of the loop over the number of integration points on element
      
    }//  end of the loop over the number of elements 
  }//  end of the if (Mp->eigstrains==2)



  if (Cp->eigstrains>0){
    delete [] id;
  }
  
  // ************************************************
  //  general functions describing the eigenstrains
  // ************************************************
  if (Cp->eigstrains==1 || Cp->eigstrains==2 || Cp->eigstrains==4 || Cp->eigstrains==5){
    //  the number of general functions needed for eigenstrain description
    xfscanf (in,"%k%ld","number_of_gfunct_eigenstr",&ngfes);
    
    eigstrfun = new gfunct [ngfes];
    
    //  loop over the number of general functions
    for (i=0;i<ngfes;i++){
      eigstrfun[i].read (in);
    }//  end of the loop over the number of general functions
    
    //  computation of eigstrains
    eigstrain_computation (Cp->time);
    
  }//  end of the statement if (Mp->eigstrains==1 || Mp->eigstrains==2){

}

/**
   function computes eigenstrains
   
   @param time - actual time
   
   11.4.2019
*/
void coupbclc::eigstrain_computation (double time)
{
  long i, j, k, ipp, nip, ncomp, gfid;
  vector coord(4); // three spatial coordinates  + time
  const char* namev[4] = { "x", "y", "z", "t" };
  
  for (i=0;i<Ct->ne;i++){
    ipp=Ct->elements[i].ipp;
    nip = Ct->give_nip (i);
    for (j=0;j<nip;j++){

      //  coordinates of integration points
      //Asymlq->ipcoord(eid, ipp, ri, ci, ipcoord);
      //Asymqq->ipcoord(eid, ipp, ri, ci, ipcoord);
      //ipcoord (i,ipp,ii,jj,coord);

      //  the true number of strain components in the integration point
      ncomp = Cm->ip[ipp].ncompstr;
      coord(3) = time;
      //  loop over the number of eigenstrain components
      for (k=0;k<ncomp;k++){
	//  id of the general function
	gfid = Cm->eigstrid[ipp][k];
	
	if (Cp->eigstrains < 3)
	  Cm->eigstrains[ipp][k] = eigstrfun[gfid].getval(coord, namev);
	if (Cp->eigstrains > 3)
	  Cm->eigstresses[ipp][k] = eigstrfun[gfid].getval(coord, namev);
      }//  end of the loop over the number of eigenstrain components
      ipp++;
    } // end of loop over the number of ip in the given block
  }
}
