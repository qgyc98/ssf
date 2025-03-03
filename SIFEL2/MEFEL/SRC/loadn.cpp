#include "loadn.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "gtopology.h"
#include "intools.h"
#include "stochdriver.h"



/**
  The constructor inializes attributes to zero values.

  Created by JK,
*/
loadn::loadn()
{
  nid = nlc = 0L;
  f = NULL;
}



/**
  The destructor deallocates used memory.

  Created by JK,
*/
loadn::~loadn()
{
  delete [] f;
}



/**
  The function reads element load from the opened text file
  given by the parameter in.

  @param in - pointer to the opened text file

  @retval 0 - on success

  Created by JK,
  Modified by Tomas Koudelka, 06.2009
*/
long loadn::read(XFILE *in)
{
  long ndofn;

  xfscanf(in, "%ld",&nid);
  if (nid<1)
    print_err("number of loaded node is nonpositive", __FILE__, __LINE__, __func__);
  if (nid>Mt->nn)
    print_err("number of loaded node (%ld) > total number of nodes (%ld)", __FILE__, __LINE__, __func__, nid, Mt->nn);
  nid--;
  ndofn = Mt->give_ndofn (nid);
  f = new double [ndofn];
  for (long i = 0; i < ndofn; i++)
    xfscanf(in, "%le", &f[i]);
  return(0);
}



/**
  The function prints element load into the opened text file
  given by the parameter out.

  @param out - pointer to the opened text file

  @retval 0 - on success

  TKr, 07/02/2013 according to read(XFILE *in)
*/
long loadn::print(FILE *out)
{
  long ndofn;

  fprintf (out, "\n  %ld",nid+1);
  ndofn = Mt->give_ndofn (nid);

  for (long i = 0; i < ndofn; i++)
    fprintf(out, "  %le", f[i]);
  return(0);
}



/**
  The function reads nodal load data from the opened text file
  given by the parameter in. It is used in the old version of mechprep preprocessor.

  @param in - pointer to the opened text file
  @param ndof - maximum nuber of dofs in the nodes 
  @param lc - total number of load cases

  @retval 0 - on success
  @retval 2 - wrong load case number

  Created by Tomas Koudelka
*/
long loadn::read_prop(FILE *in, long ndof, long lc)
{
  getlong(in, nlc);
  if ((nlc < 1) || (nlc > lc))
    return(2);
  f = new double [ndof];
  for (long i = 0; i < ndof; i++)
    getdouble(in, f[i]);
  return(0);
}



/**
  The function changes material parameters for stochastic analysis.
   
  @param atln - selected load components (components which are changed)
  @param val - array containing new values of parameters
   
  @return The function does not return anything.

  Created by JK,
*/
void loadn::changeparam (atsel &atln,vector &val)
{
  long i;
  
  for (i=0;i<atln.num;i++){
    f[atln.atrib[i]]=val[i];
  }
}
