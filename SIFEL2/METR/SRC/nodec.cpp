#include <string.h>
#include "nodec.h"
#include "vector.h"
#include "global.h"



/**
   Constructor initializes data members to zero or default values.
   
   Created by JK,
*/
nodec::nodec (void)
{
  //  indicator of local coordinate system on node
  transf=0;
  //  base vectors of local coordinate system
  e1=NULL;  e2=NULL;  e3=NULL;
  
  react=0;  r=NULL;
  
  //  number of strain/stress components
  ncompstr=0;
  //  number of components in array other
  ncompother=0;
  //  array containing nodal strains
  strain=NULL;
  //  array containing nodal stresses
  stress=NULL;
  //  array containing nodal values of the array other
  other=NULL;
  //  number of contributions to the array strain
  ncontr_strain=NULL;
  //  number of contributions to the array stress
  ncontr_stress=NULL;
  //  number of contributions to the array other
  ncontr_other=0;
  //  volume used for strain contributions
  vol_strain=NULL;
  //  volume used for stress contributions
  vol_stress=NULL;
  //  volume used for other contributions
  vol_other=0;
  
  
  pstra = NULL;
  pstre = NULL;
  meaning = NULL;

  nodval = NULL;
}



/**
  Destructor releases allocated memory of the node object.

  Created by JK,
*/
nodec::~nodec (void)
{
  delete [] e1;  delete [] e2;  delete [] e3;

  delete [] strain;  delete [] stress;  delete [] other;
  delete [] ncontr_strain;  delete [] ncontr_stress;
  delete [] vol_strain;  delete [] vol_stress;

  delete [] r;  
  delete [] pstra;
  delete [] pstre;
  delete [] meaning;
  delete [] nodval;
}



/**
  Function reads nodal data from opened text file.
   
  @param in - pointer to teh opened text file
   
  @return The function does not return anything.

  Created by JK
*/
void nodec::read (XFILE *in)
{
  //  transformation to local system
  xfscanf (in,"%ld",&transf);
  if (transf!=0 && transf!=2 && transf!=3)
    print_err("wrong identification of local system", __FILE__, __LINE__, __func__);
  
  if (transf==2){
    e1 = new double [2];  e2 = new double [2];
    xfscanf (in,"%lf %lf",e1+0,e1+1);
    xfscanf (in,"%lf %lf",e2+0,e2+1);
  }
  if (transf==3){
    e1 = new double [3];  e2 = new double [3];  e3 = new double [3];
    xfscanf (in,"%lf %lf %lf",e1+0,e1+1,e1+2);
    xfscanf (in,"%lf %lf %lf",e2+0,e2+1,e2+2);
    xfscanf (in,"%lf %lf %lf",e3+0,e3+1,e3+2);
  }

}

/**
  Function prints nodal data into opened text file.
   
  @param out - pointer to teh opened text file
   
  @return The function does not return anything.

  Created by TKr, 02/01/2013
*/
void nodec::print (FILE *out)
{
  //  type of cross section
  fprintf (out,"  %d", int(crst));
  if (crst!=0){
    fprintf (out," %ld",idcs+1);
  }

  //  transformation to local system
  fprintf (out,"  %ld",transf);

  if (transf==2){
    fprintf (out," %lf %lf",e1[0],e1[1]);
    fprintf (out," %lf %lf",e2[0],e2[1]);
  }
  if (transf==3){
    fprintf (out," %lf %lf %lf",e1[0],e1[1],e1[2]);
    fprintf (out," %lf %lf %lf",e2[0],e2[1],e2[2]);
    fprintf (out," %lf %lf %lf",e3[0],e3[1],e3[2]);
  }
  fprintf (out,"\n");
}
