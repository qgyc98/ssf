#include <string.h>
#include <stdlib.h>
#include "bocont.h"
#include "parser.h"
#include "globalt.h"



/**
  Constructor initializes data members
  Parameters :
    none
  Returns :
    nothing

  Created by TKo, 09.2010
*/
bocont::bocont()
{
  lcid = -1;
  con  =  0.0;
  iv   =  0.0;
  cgf  =  NULL;
}



/**
  Destructor releases allocated memory
  Parameters :
    none
  Returns :
    nothing

  Created by TKo, 09.2010
*/
bocont::~bocont()
{
  if (cgf)
    delete cgf;
}



/**
  Method reads data from the text file specified by parameter in, data are
    in preprocessor format

    @param in    - pointer to the opened text file, where the data will be read
    @param ndofn - number of dofs in given node with prescribed condition

    @retval 0 - on success
    @retval 1 - load case id is not in range <1, ndofn>
    @retval 2 - parser error

  Created by TKo, 09.2010
*/
long bocont::read(XFILE *in, long ndofn)
{
  xfscanf(in, "%k%ld", "lc_id", &lcid); // reading direction
  if ((lcid < 1) || (lcid > ndofn))
  // load case id should be in range <1, ndof>
    return(1);
  lcid--;
  xfscanf(in, "%k%le", "ini_cd", &iv);
  switch (Tp->tprob)
  {
    case stationary_problem :
    case nonlinear_stationary_problem:
      xfscanf(in, "%k%le", "cond", &con); // reading condition
      break;
    case nonstationary_problem :
    case nonlinear_nonstationary_problem :
    case growing_np_problem:
    case growing_np_problem_nonlin:
    case discont_nonstat_problem:
    case discont_nonlin_nonstat_problem:
      cgf = new gfunct;
      if (cgf->read(in))
      // error parsing gfunct
        return(2);
      break;
    default :
      print_err("Unknown type of problem required", __FILE__, __LINE__, __func__);
      abort();
  }
  return(0);
}



/**
  Method prints data to the text file specified by parameter out.

  @param in    - pointer to the opened text file, where the data will be read
  @param ndofn - number of dofs in given node with prescribed condition

  @retval 0 - on success

  Created by TKo, 09.2010
*/
long bocont::print(FILE *out)
{
  fprintf(out, "%le ", iv);
  if (cgf)
    cgf->print(out);
  else
    fprintf(out, "%le\n", con);

  return(0);
}



/**
  Method compares the actual object with the instance tbc.

  @param tbc  - object of boundary condition which will be compared to

  @retval 0 - in the case of identical boundary conditions
  @retval 1 - in the case of difference between objects

  Created by TKo, 09.2010
*/
long bocont::compare(bocont &tbc)
{
  // check identical load case id
  if (lcid != tbc.lcid)
    return 1;

  // check identical initial conditions
  if (iv != tbc.iv)
    return 1;

  // check identical conditons

  if (cgf) // nonstationary problems
    return cgf->compare(*tbc.cgf);

  if (con != tbc.con) // stationary problems
    return 1;

  return 0;
}
