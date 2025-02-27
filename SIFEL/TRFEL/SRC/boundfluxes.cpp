#include "boundfluxes.h"
#include "globalt.h"
#include "aliast.h"
#include "globmatt.h"
#include "genfile.h"
#include "elemswitcht.h"
#include <math.h>
#include <limits.h>

boundfluxes::boundfluxes (void)
{
  //  load case id
  lcid=0;
  
  //  number of element with boundary conditions
  neb = 0;
  
  elemload = NULL;
}

boundfluxes::~boundfluxes (void)
{
  delete [] elemload;
}

/**
   function reads load case characteristics
   
   @param in - pointer to input stream
   @param lcid - load case id
   
   JK, 4. 3. 2018
*/
void boundfluxes::read (XFILE *in,long lcid)
{
  long i,eid,nbo,nnbo;
  
  // **********************************************
  //  elements with boundary conditions
  // **********************************************
  //  number of elements with boundary conditions
  xfscanf (in,"%ld",&neb);
  if (neb<0)
    print_err("negative number of elements with boundary conditions",__FILE__,__LINE__,__func__);

  elemload = new loadelt [neb];
  if (Mesprt==1)
    fprintf (stdout,"\n number of elements with boundary edges %ld",neb);
  
  for (i=0;i<neb;i++){
    //  element id
    xfscanf (in,"%ld",&eid);
    eid--;
    //  check of element id
    if (eid > Tt->ne-1){
      print_err("Element number in boundary conditions is greater than total number of elements", __FILE__, __LINE__, __func__);
      abort();
    }
    if (eid < 0){
      print_err("Element number in boundary conditions is less than 0", __FILE__, __LINE__, __func__);
      abort();
    }
    
    //  number of boundary objects (end nodes, edges, surfaces)
    //  number of nodes on each boundary object
    Tt->give_nbobjects (eid,nbo,nnbo);
    
    elemload[i].read_comp (in,lcid,eid,nbo,nnbo);
  }  
}


/**
   function prints load case characteristics
   
   @param out - pointer to output stream
   @param lcid - load case id
   
   JK, 4. 3. 2018
*/
void boundfluxes::print (FILE *out,long lcid)
{
  long i;
  
  //fprintf (out,"\n%ld\n",lcid);

  // **********************************************
  //  number of elements with boundary conditions
  // **********************************************
  fprintf(out,"\n## loaded elements:");
  fprintf (out,"\n%ld\n\n",neb);
  for (i=0;i<neb;i++){
    elemload[i].print (out,lcid);
  }
  
}


