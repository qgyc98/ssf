#include "advectvel.h"
#include "iotools.h"
#include "mathem.h"
#include <math.h>
#include <string.h>
#include <ctype.h>



/**
  This constructor initializes data to the zero values
*/
advectvel::advectvel()
{
  ncomp = 0L;
  lcid  = 0L;
  v = NULL;
}



/**
  This destructor deallocates used memory
*/
advectvel::~advectvel()
{
  delete [] v;
}



/**
  This function reads data about entity load from the text file given by the parameter in

  @param in - pointer to the opned text file
  @param ndof - total number of nodal DOFs i.e. total number of load cases

  @retval 0 - on succes
  @retval 1 - unsupported function type is required
  @retval 2 - invalid type of load coordinate system
  @retval 3 - function has insufficient number of parameters, wrong names of parameters or invalid function type

  Created by Tomas Koudelka, 5.2016
*/
long advectvel::read(XFILE *in, long /*ndof*/)
{
  long i;
  long var_x, var_y, var_z;

  // read load case id, i.e. medium id or DOF id
  xfscanf(in, "%k%ld", "lc_id", &lcid);
  // if ((lcid < 1) || (lcid > ndof))
  // {
  //   print_err("Invalid load case id %ld, it must be in range <1;%ld>\n", __FILE__ ,__LINE__, __func__, lcid, ndof);
  //   return 1;
  // }      
  // the following two condition statements can be removed after implementation of full support for more media
  if (lcid > 1)
  {
    print_err("Advection velocities for more media (lcid > 1) have not been implemented yet\n", __FILE__ ,__LINE__, __func__);
    return 1;
  }      
  if (lcid < 1)
  {
    print_err("Invalid load case id must be > 0 (lcid=%ld)\n", __FILE__ ,__LINE__, __func__, lcid);
    return 1;
  }      

  // read the number of velocity components which is equaled to the space dimension of the problem solved
  xfscanf(in, "%k%ld", "ncomp", &ncomp);
  if ((ncomp < 1) && (ncomp > 3))
  {
    print_err("invalid number of components of velocity vector - ncomp=%ld is out of range <1;3>", 
              __FILE__, __LINE__,__func__, ncomp);
    return 2;
  }

  // read particluar velocity components
  xfscanf(in, "%k", "velocity_comp");
  v = new gfunct[ncomp];
  for(i=0; i<ncomp; i++)
  {
    v[i].read(in);
    switch(v[i].tfunc) // check the type of supported function and tehir parameters
    {
      case constant:
        break;
      case pars:
        if (v[i].eq[0]->Variables.count() > 3)
        {
          print_err("Line %ld, col %ld, file %s: number of parameters %ld of velocity distribution function > 3 for component %ld", 
                    __FILE__, __LINE__, __func__, in->line, in->col, in->fname, v[i].eq[0]->Variables.count(), i+1);
          return 3;
        }
        if (var2coord(v[i].eq[0], var_x, var_y, var_z))
        {
          print_err("Line %ld, col %ld, file %s: wrong names of parameters of velocity distribution function for component %ld."
                    "\n Valid names are x, y or z.", 
                    __FILE__, __LINE__, __func__, in->line, in->col, in->fname, i+1);
          return 3;
        }
        break;
      default:
        print_err("Line %ld, col %ld, file %s: unsupported type %d of velocity distribution function\n"
                  " for velocity component %ld", __FILE__, __LINE__, __func__, in->line, in->col, int(v[i].tfunc), i+1);
        return 3;
    }
  }
  return 0;
}



/**
  The function returns indeces of 'x' parameter, 'y' parameter and 'z' parameter
  of the function defined in eq. Performed comparison is case insesitive.

  @param eq - pointer to the parsed function
  @param var_x - index of x-coordinate (used by eq->Variables.at function)
  @param var_y - index of y-coordinate (used by eq->Variables.at function)
  @param var_z - index of z-coordinate (used by eq->Variables.at function)
  
  @retval 0 - on success
  @retval 1 - wrong variable names used
  @return Indeces are returned via parameters var_x, var_y, var_z.
    Negative value of the index is returned in case that it is not found.
*/
long advectvel::var2coord(Equation *eq, long &var_x, long &var_y, long &var_z)
{
  variable *var;
  long i, var_det = 0L;
  char vn;

  var_x = var_y = var_z = -1;

  for (i=0; i<eq->Variables.count(); i++)
  {
    var = eq->Variables.at(i);  // get i-th function variable
    vn = tolower(var->Name[0]); // take the first character of the function variable name
    if (vn == 'x') // x-coordinate detected -> store its id
    {
      var_x = i;
      var_det++;
    }
    if (vn == 'y') // y-coordinate detected -> store its id
    {
      var_y = i;
      var_det++;
    }
    if (vn == 'z') // z-coordinate detected -> store its id
    {
      var_z = i;
      var_det++;
    }
  }
  if (eq->Variables.count() != var_det)
    return 1;
  return 0;
}   



/**
  This function returns value of compid-th advection velocity component at
  the given node

  @param tn - node where the load should be approximated
  @param compid - velocity component id
  @param vc - velocity component value (output)

  @return 0 - on succes
  @retval 1 - in case of invalid component id
  @retval 2 - in case of unknown type of load function
*/
long advectvel::getval(snode &tn, long compid, double &vc)
{
  long var_x, var_y, var_z;
  variable *var;

  if ((compid < 0) || (compid >= ncomp))
  {
    print_err("invalid component id %ld is required. Component must be from range <1;%ld>",
              __FILE__, __LINE__, __func__, compid+1, ncomp);
    return 1;
  }
  switch (v[compid].tfunc)
  {
    case constant:
      vc = v[compid].f;
      break;
    case pars:
      var2coord(v[compid].eq[0], var_x, var_y, var_z);
      // setting of x coordinate
      if (var_x >= 0)
      {
        var = v[compid].eq[0]->Variables.at(var_x);
        var->Value = tn.x;
      }
      // setting of y coordinate
      if (var_y >= 0)
      {
        var = v[compid].eq[0]->Variables.at(var_y);
        var->Value = tn.y;
      }
      // setting of z coordinate
      if (var_z >= 0)
      {
        var = v[compid].eq[0]->Variables.at(var_z);
        var->Value = tn.z;
      }
      vc = v[compid].eq[0]->Evaluate();
      break;
    default:
      return 2;
  }
  return(0);
}
