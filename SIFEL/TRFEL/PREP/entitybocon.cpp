#include "iotools.h"
#include "mathem.h"
#include "globalt.h"
#include "entitybocon.h"
#include <math.h>
#include <string.h>
#include <ctype.h>



/**
  This constructor initializes data to the zero values

  Created by TKo, 09.2010    
*/
entitybocon::entitybocon()
{
  it = interpoltype(0);
  dc = 0;
  dt = 0;
}



/**
  This destructor deallocates used memory

  Created by TKo, 09.2010    
*/
entitybocon::~entitybocon()
{
}



/**
  This function reads data about entity load from the text file given by the parameter in

  @param in - pointer to the opned text file

  @retval 0 - on succes
  @retval 1 - no coordinates nor time variables has been detected

  Created by TKo, 09.2010    
*/
long entitybocon::read(XFILE *in)
{
  long res;

  gf.read(in);
  switch (gf.tfunc)
  {
    case constant:
    case tab:
    case tab_file:
      return 0;
    case pars:
      res = checkvar2coord();
      if(res < 1)
        return 1;        
/*      if (dt == 0)
      { // no time dependency was detected in parser -> function depends on coordinates only
        // the function will be converted to the constant type -> 'it' flag is set to arbitrary nonzero value
        // in order to signal, that the expression does not contain time.
        // It will not be used for specification of any kind of interpolation.
        it = interpoltype(20);
      }*/
      return 0;
    case pars_set:
      res = checkvar2coord();
      if(res < 1)
        return 1;
      if (dt == 0)
      { // no time dependency was detected in parser -> function depends on coordinates only
        // the function will be converted to the table -> interpolation type is required
        xfscanf(in, "%k%m", "approx_type", &interpoltype_kwdset, (int*)&it);
      }
      return 0;
    default:
      print_err("unsupported type of gfunct - 'constant', 'tab', 'pars' or 'pars_set' can be specified", __FILE__, __LINE__, __func__);
      abort();
  }
  return 0;
  
}



/**
  The function detects 'x' parameter, 'y' parameter ,'z' parameter and time 't'
  in the gf. Performed comparison is case insesitive. In the case of time parameter 't'
  detection, the dc is set to nonzero value. In the case of coordinate detection, the dt is 
  set to nonzero value.

  @return  The function returns the number of detected parameters x,y,z,t.
 
  Created by TKo, 09.2010    
*/
long entitybocon::checkvar2coord()
{
  long i, j;
  long ret = 0;
  long n = 0;

  var_x = var_y = var_z = -1;
  dt = 0;
  dc = 0;

  switch(gf.tfunc)
  {
    case pars:
    case pars_set:
      for(i=0; i<gf.neqs; i++)
      {
        for(j=0; j<gf.eq[i]->Variables.count(); j++)
        {
          if (tolower(gf.var[i][j]->Name[0]) == 'x')
          {
            var_x = j;
            dc = 1;
            n++;
          }
          if (tolower(gf.var[i][j]->Name[0]) == 'y')
          {
            var_y = j;
            dc = 1;
            n++;
          }
          if (tolower(gf.var[i][j]->Name[0]) == 'z')
          {
            var_z = j;
            dc = 1;
            n++;
          }
          if (tolower(gf.var[i][j]->Name[0]) == 't')
          {
            dt = 1;
            n++;
          }
        }
        if (n != gf.eq[0]->Variables.count())
        {
          print_err("unknown parameters of function detected,\n"
                    " parameters must be denoted by x,y,z and t only.", __FILE__, __LINE__, __func__);
          abort();
        }
        if (n > ret)
          ret = n;
        n = 0;
      }
      break;
    default:
      print_err("unknown type of gfunct is required"
                " only 'pars' and 'pars_set' can be used", __FILE__, __LINE__, __func__);
      abort();
  }
  return ret;
}   



/**
  The function detects 'x' parameter, 'y' parameter ,'z' parameter and time 't'
  in the gf. Performed comparison is case insesitive.

  @param i - index of parser expression whose variables should be detected. It can be in range <0;f.neqs)
  
  @return  Class attributes var_x, var_y, var_z are set to non-negative values in the case of their 
           detection in parsed string. The function returns number of detected parameters x,y,z,t.
 
  Created by TKo, 09.2010    
*/
long entitybocon::var2coord(long i)
{
  long j;
  long ret = 0;

  var_x = var_y = var_z = -1;

  if ((i < 0) && (i>=gf.neqs))
  {
    print_err("index of parsed function is out of range", __FILE__, __LINE__, __func__);
    abort();
  }

  switch(gf.tfunc)
  {
    case pars:
    case pars_set:
      for(j=0; j<gf.eq[i]->Variables.count(); j++)
      {
        if (tolower(gf.var[i][j]->Name[0]) == 'x')
        {
          var_x = j;
          ret++;
        }
        if (tolower(gf.var[i][j]->Name[0]) == 'y')
        {
          var_y = j;
          ret++;
        }
        if (tolower(gf.var[i][j]->Name[0]) == 'z')
        {
          var_z = j;
          ret++;
        }
        if (tolower(gf.var[i][j]->Name[0]) == 't')
        {
          ret++;
        }
      }
      if (ret != gf.eq[i]->Variables.count())
      {
        print_err("unknown parameters of function detected,\n"
                  " parameters must be denoted by x,y,z and t only.", __FILE__, __LINE__, __func__);
        abort();
      }
      return ret;
    default:
      print_err("unknown type of gfunct is required"
                " only 'pars' and 'pars_set' can be used", __FILE__, __LINE__, __func__);
      abort();
  }
  return ret;
}   



/**
  The function approximates values of BC on the given entity to the given element node.
  The function requires gf.tfunc to be pars or pars_set. Corresponding
  expressions have to depend on one of space coordinate at least (x,y,z) or
  one of space coordinate(x,y,z) at least and time (t)
  
  @param tn - node where the BC should be approximated
  @param tgf - object of gfunct that is initialized by the BC function for the given node

  @return The function sets new function in tgf according to coordinates of the given node.
  @retval 0 - on success
  @retval 1 - in the case of wrong type of func (it can be pars or pars_set)
  @retval 2 - in the case of an error in parsing new function in tgf
*/
long entitybocon::getval(snode &tn, gfunct &tgf)
{
  long i, j;
  Parser par;
    
  switch (gf.tfunc)
  {
    case constant:  
    case tab:
    case tab_file:  
      tgf.copy(gf);
      break;
    case pars:
      if (dc)
      {
        if (dt)
        {
          // expression depends on x,y,z and time -> conversion to other pars with substituted coordinates
          // which depends only on time
          tgf.tfunc  = pars;
          tgf.neqs   = 1;
          tgf.eq     = new Equation* [1];
          tgf.var    = new variable**[1];
          tgf.func   = new char*     [1];
          tgf.func[0] = substbcstr(gf.func[0], tn);
          tgf.eq[0]   = par.TextToTree(tgf.func[0]);
          if (tgf.eq[0] == NULL)
          {
            print_err("parsing expression", __FILE__, __LINE__, __func__);
            return (2);
          }
          tgf.var[0]  = new variable* [tgf.eq[0]->Variables.count()];
          for(j=0; j < tgf.eq[0]->Variables.count(); j++)
            gf.var[0][j] = tgf.eq[0]->Variables.at(j);
        }
        else
        {
          // expression depends on x,y and z only -> conversion to constant function
          tgf.tfunc = constant;
          var2coord(0);
          setvars(tn, 0);
          tgf.f = gf.eq[0]->Evaluate();
        }
      }
      else
      {
        // expression does not depend on coordinates -> copy plain parser
        tgf.copy(gf);
      }
      break;
    case pars_set:
      if (dc)
      {
        if (dt)
        {
          // particular expressions depend on x,y,z and t -> conversion to pars_set with substituted coordinates
          // which depends only on time
          tgf.tfunc  = pars_set;
          tgf.neqs   = gf.neqs;
          tgf.eq     = new Equation* [gf.neqs];
          tgf.var    = new variable**[gf.neqs];
          tgf.func   = new char*     [gf.neqs];
          tgf.limval = new double    [gf.neqs];
          for(i=0; i<gf.neqs; i++)
          {
            tgf.limval[i] = gf.limval[i];
            tgf.func[i]   = substbcstr(gf.func[i], tn);
            tgf.eq[i]     = par.TextToTree(tgf.func[i]);
            if (tgf.eq[i] == NULL)
            {
              print_err("parsing expression", __FILE__, __LINE__, __func__);
              return (2);
            }
            tgf.var[i] = new variable* [tgf.eq[i]->Variables.count()];
            for(j=0; j < gf.eq[i]->Variables.count(); j++)
              tgf.var[i][j] = tgf.eq[i]->Variables.at(j);
          }
        }
        else
        {
          // particular expressions depend on x,y and z only -> conversion to tablefunc
          tgf.tfunc = tab;          
          // create new table one item larger than original pars_set due to lower limit
          tgf.tabf = new tablefunct(gf.neqs+1);       
          // set type of interpolation 
          tgf.tabf->itype = it;
          for(i=0; i<gf.neqs; i++)
          {
            tgf.tabf->x[i+1] = gf.limval[i];
            var2coord(i);
            setvars(tn, i);
            tgf.tabf->y[i+1] = gf.eq[i]->Evaluate();
          }
          // add lower limit for correct evaluation of table value
          // in the first interval (it is not included in the original
          // definition of pars_set
          tgf.tabf->x[0] = tgf.tabf->x[1]-1.0;
          tgf.tabf->y[0] = tgf.tabf->y[1];
        }
      }
      else
      {
        // particular expressions do not depend on coordinates -> copy plain pars_set
        tgf.copy(gf);
      }
      break;
    default:
      print_err("unknown type of function is required", __FILE__, __LINE__, __func__);
      return 1;
  }
  return(0);
}



/**
  The function sets x,y and z variables of func to the coordinates of 
  node tn. Attributes var_x, var_y and var_z have to contain indeces of coordinate variables
  for i-th expression of gf.

  @param tn - structure of node whose coordinates will be used
  @param i  - index of expression of func whose variables will be set
    
  @return The function sets variables of gf.eq[i] to the coordinates of given node.

  Created by Tomas Koudelka, 09.2010
*/
void entitybocon::setvars(snode &tn, long i)
{
  if (var_x > -1)
    gf.var[i][var_x]->Value = tn.x;

  if (var_y > -1)
    gf.var[i][var_y]->Value = tn.y;

  if (var_z > -1)
    gf.var[i][var_z]->Value = tn.z;
}



/**
  The function substitutes occurence of x, y and z variables in the expr by
  the real values of coordinates of the given node tn. A new string with substituted 
  values is returned.
  Class attributes var_x, var_y and var_z have to be initialized and
  they have negative values in the case that the given variable is not present in the string.

  @param expr - string with the expression
  @param tn   - node structure with coordinates

  @return The fuction returns string with substituted occurences of coordinates by the real values.
  
  Created by TKo, 09.2010
 */
char *entitybocon::substbcstr(const char *expr, snode &tn)
{
  char *src, *ret;
  char *loc, *ploc, *aptr;
  long l = long(strlen(expr));
  long ln;
  long i, nx, ny, nz;
  int  n;
  long lalpha, ralpha;

  // compute number of occurences of the coordinates
  src = new char[l];
  memcpy(src, expr, l);
  for(i=0; i<l; i++)
    src[i] = tolower(src[i]);
  loc = src;
  nx = ny = nz = 0L;
  do
  {  
    loc = strpbrk(loc, "xyz");
    if (loc)
    {
      if (loc > src)
        lalpha = isalpha(*(loc-1)); // test for left neighbour char of located coordinate
      else
        lalpha = 0;                 // located coordinate is at the beginning of the src

      ralpha = isalpha(*(loc+1));   // test for right neighbour char of located coordinate
                                    // loc can be incremented because the last char has to be \0x0
      if (lalpha+ralpha == 0)
      {
        // both neighbour chars are not alpha -> detected occurence of coordinate 
        // variable(not name of a function)
        if (*loc == 'x')
          nx++;
        if (*loc == 'y')
          ny++;
        if (*loc == 'z')
          nz++;
      }
    }
  }while(loc);

  // create new string with substituted coordinates
  ln = l + (nx+ny+nz)*(1+1+1+1+15+1+1+3+1)+1;  // calculate necessary size of the new string 
                 // ( - x . XX e + X ) \0
  ret = new char[ln];
  aptr = ret;
  loc = ploc = src;
  ploc--; // for the first copying of the precedings chars the ploc has to 
          //refer before the first char of src
  do
  {  
    loc = strpbrk(loc, "xyz");
    if (loc)
    {
      if (loc > src)
        lalpha = isalpha(*(loc-1)); // test for left neighbour char of located coordinate
      else
        lalpha = 0;                 // located coordinate is at the beginning of the src

      ralpha = isalpha(*(loc+1));   // test for right neighbour char of located coordinate
                                    // loc can be incremented because the last char has to be \0x0
      if (lalpha+ralpha == 0)
      {
        // both neighbour chars are not alpha -> detected occurence of coordinate 
        // variable(not name of a function)

        // copy preceding chars before the actual occurence of the coordinate
        memcpy(ret, ploc+1, loc-ploc-1);
        aptr += loc-ploc-1; // change actual position in the returned string

        // substitute the coordinate variables by real values from node
        if (*loc == 'x')
        {
          sprintf(aptr, "(%.15le)%n", tn.x, &n);
          aptr += n;
        }
        if (*loc == 'y')
        {
          sprintf(aptr, "(%.15le)%n", tn.y, &n);
          aptr += n;
        }
        if (*loc == 'z')
        {
          sprintf(aptr, "(%.15le)%n", tn.z, &n);
          aptr += n;
        }
        ploc = loc; // ploc has to contain pointer to the coordinate located in the previous step
      }
    }
  }while(loc);
  
  // copy rest of original string including terminating \0x0
  memcpy(aptr, ploc+1, (src+l+1)-(ploc+1));
                     // (ptr to end of src[including term \0) - (ptr to successive char after last substituted char)
  aptr += src+l+1-ploc-1;

  delete [] src;

  return ret;
}
