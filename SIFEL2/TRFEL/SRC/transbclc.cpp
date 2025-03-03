#include "transbclc.h"
#include "globalt.h"

transbclc::transbclc (void)
{
  // number of load cases
  nlc=0;
  //  load case pointers
  lc=NULL;

  //  the number of boundaries with output fluxes
  nbf=0;
  //  output boundary fluxes
  bf=NULL;
  //  fluxes
  fluxes=NULL;
}

transbclc::~transbclc (void)
{
  delete [] lc;
  delete [] bf;
  delete [] fluxes;
}

/**
   function reads data about load cases
   
   @param in - pointer to input file
   
   JK
*/
void transbclc::read (XFILE *in)
{
  long i,j,k,eid,nbo,nnbo;
  
  xfscanf (in,"%k%ld","number_of_load_cases",&nlc);
  lc = new loadcaset [nlc];
  for (i=0;i<nlc;i++){
    lc[i].read (in,i);
  }
  
  //  the number of instances bf is equal to the number of load cases
  bf = new boundfluxes [nlc];
  for (i=0;i<nlc;i++){
    bf[i].read (in,i);
  }
  
  //  number of required fluxes
  nbf=0;
  for (i=0;i<nlc;i++){
    for (j=0;j<bf[i].neb;j++){
      //  element id
      eid=bf[i].elemload[j].eid;
      //  nbo - the number of boundary object (edges, surfaces)
      Tt->give_nbobjects (eid,nbo,nnbo);
      for (k=0;k<nbo;k++){
	if (nbf<bf[i].elemload[j].nvid[k]+1)
	  nbf=bf[i].elemload[j].nvid[k]+1;
      }
    }
  }
  
  fprintf (stdout,"\n the number of fluxes required  %ld",nbf);

  fluxes = new double [nbf];
  
}


/**
   function prints data about load cases
   
   @param out - pointer to output file
   
   TKr
*/
void transbclc::print (FILE *out)
{
  long i;
  
  fprintf(out,"\n\n## boundary conditions:\n");
  fprintf (out,"\n");
  fprintf (out,"\n %ld #number of loadcases",nlc);
  for (i=0;i<nlc;i++){
    lc[i].print (out,i);
  }

  fprintf (out,"\n");
  for (i=0;i<nlc;i++){
    bf[i].print (out,i);
  }
}

/**
   function searches for elements with nodes where source is defined
   
   JK, 25.6.2005
*/
void transbclc::elemsource ()
{
  long i;
  
  for (i=0;i<nlc;i++){
    lc[i].elemsource ();
  }
  
}

