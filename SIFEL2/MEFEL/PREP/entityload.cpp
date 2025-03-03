#include "entityload.h"
#include "iotools.h"
#include "gfunct.h"
#include "mathem.h"
#include <math.h>
#include <string.h>
#include <ctype.h>



/**
  This constructor initializes data to the zero values
*/
entityload::entityload()
{
  ncomp = 0L;
  nlc   = 0L;
  nslc  = 0L;
  ft    = constant;
  lgcs  = 0L;
  val   = NULL;
  func  = NULL;
  tfunc = NULL;
}



/**
  This destructor deallocates used memory
*/
entityload::~entityload()
{
  long i;
  if (func)
  {
    for (i=0; i <ncomp; i++)
    {
      delete func[i];
      delete [] tfunc[i];
    }
    delete [] func;
    delete [] tfunc;
  }
  delete [] val;
}



/**
  This function reads data about entity load from the text file given by the parameter in

  @param in - pointer to the opned text file
  @param lc - total number of load cases
  @param slc - pointer to array with number of subloadcases for particular load cases, default value is NULL

  @retval 0 - on succes
  @retval 1 - unsupported function type is required
  @retval 2 - invalid type of load coordinate system
  @retval 3 - function has insufficient number of parameters
  @retval 4 - wrong name of function parameters
*/
long entityload::read(XFILE *in, long lc, long *slc=NULL)
{
  Parser par;
  long var_x, var_y, var_z;

  xfscanf(in, "%k%ld", "lc_id", &nlc);
  if ((nlc < 1) || (nlc > lc))
  {
    print_err("Line %ld, col %ld, file %s: invalid load case index %ld.\n Valid range is <1,%ld>", 
              __FILE__, __LINE__, __func__, in->line, in->col, in->fname, nlc, lc);
    return(1);
  }
  if (slc)
  {
    xfscanf(in, "%k%ld", "slc_id", &nslc);
    if ((nslc < 1) || (nslc > slc[nlc-1]))
    {
      print_err("Line %ld, col %ld, file %s: invalid subload case index %ld.\n Valid range is <1,%ld>", 
                __FILE__, __LINE__, __func__, in->line, in->col, in->fname, nslc, slc[nlc-1]);
      return(1);
    }
  }
  long i;
  xfscanf(in, "%k%ld%k%m", "ncomp", &ncomp, "func_type", &generalfunct_kwdset, &ft);
  xfscanf(in, "%k%ld", "coord_sys", &lgcs);
  if ((lgcs < 1) || (lgcs > 2))
  {
    print_err("Line %ld, col %ld, file %s: unknown type of coordinate system %ld.\n Valid values are 1 or 2", 
              __FILE__, __LINE__, __func__, in->line, in->col, in->fname, lgcs);
    return 2;
  }
  xfscanf(in, "%k", "load_comp");
  switch (ft)
  {
    case constant:
      val = new double [ncomp];
      memset (val, 0, sizeof(*val)*ncomp);
      for (i=0; i<ncomp; i++)
        xfscanf(in, "%le", val+i);
      break;
    case pars:
      func  = new Equation*[ncomp];      
      tfunc = new char*[ncomp];
      memset (func, 0, sizeof(*func)*ncomp);
      memset (tfunc, 0, sizeof(*tfunc)*ncomp);
      for (i=0; i<ncomp; i++)
      {
        tfunc[i] = new char[1001];
        xfscanf(in, "%1000s", tfunc[i]);
        func[i] = par.TextToTree(tfunc[i]);
        // check number of variables
        if (func[i]->Variables.count() > 3)
        {
          print_err("Line %ld, col %ld, file %s: number of parameters %ld of load distribution function > 3", 
                    __FILE__, __LINE__, __func__, in->line, in->col, in->fname, func[i]->Variables.count());
          return 3;
        }
        if (var2coord(func[i], var_x, var_y, var_z))
        {
          print_err("Line %ld, col %ld, file %s: wrong names of parameters of load distribution function."
                    "\n Valid names are x, y or z.", 
                    __FILE__, __LINE__, __func__, in->line, in->col, in->fname);
          return 4;
        }
      }
      break;
    default:
      print_err("Line %ld, col %ld, file %s: unknown type %d of load distribution function", 
                __FILE__, __LINE__, __func__, in->line, in->col, int(ft));
      return 5;
  }
  return(0);
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
long entityload::var2coord(Equation *eq, long &var_x, long &var_y, long &var_z)
{
  variable *var;
  long i, var_det = 0L;
  char vn;

  var_x = var_y = var_z = -1;



  for (i=0; i<eq->Variables.count(); i++)
  {
    var = eq->Variables.at(i);
    vn = tolower(var->Name[0]);
    if (vn == 'x')
    {
      var_x = i;
      var_det++;
    }
    if (vn == 'y')
    {
      var_y = i;
      var_det++;
    }
    if (vn == 'z')
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
  This function approximates values on the loaded entity to the element node which
  lay on it. It is used linear approximation.

  @param tn - node where the load should be approximated
  @param nv - array with approximated nodal values

  @retval 0 - on succes
  @retval 1 - in case of unknown type of load function
*/
long entityload::getval(snode &tn, double *nv)
{
  long i, var_x, var_y, var_z;
  variable *var;
    
  switch (ft)
  {
    case constant:
      for (i=0; i<ncomp; i++)
        nv[i] = val[i];
      break;
    case pars:
      for (i=0; i<ncomp; i++)
      {
        var2coord(func[i], var_x, var_y, var_z);
        // setting of x coordinate
        if (var_x >= 0)
        {
          var = func[i]->Variables.at(var_x);
          var->Value = tn.x;
        }
        // setting of y coordinate
        if (var_y >= 0)
        {
          var = func[i]->Variables.at(var_y);
          var->Value = tn.y;
        }
        // setting of z coordinate
        if (var_z >= 0)
        {
          var = func[i]->Variables.at(var_z);
          var->Value = tn.z;
        }
        nv[i] = func[i]->Evaluate();
      }
      break;
    default:
      return 1;
  }
  return(0);
}



/**
  The function creates new object of the loadel type 
  which contains edge loads generated for given element and 
  its edges with given property.

  @param el - element whose edge load will be generated
  @param prop - edge property id of the element el where the edge load will be applied
  @param nodes - array with node structures
  @param edgid   - array with edge indeces for the given property number prop detected on additional edges
  @param nedgid - the length of array edgid

  @return The function returns pointer to the new allocated object of the loadel type.

  Created by Tomas Koudelka, 06.2009
  Modified by Tomas Koudelka, 08.2011
*/
loadel *entityload::edge2loadel(selement &el, long prop, snode *nodes, long *edgid, long nedgid)
{
  long i, j, k, l, id, tnne;
  loadel *ret=new loadel;
  double *nv;

  for (i=0, tnne=0; i<el.ned; i++)
    tnne += el.nned[i];
  ret->tel  = edge;
  ret->ned  = el.ned;
  ret->nnve = tnne*ncomp;
  ret->nlc  = nlc;
  ret->nslc = nslc;
  ret->le = new long[ret->ned];
  memset(ret->le, 0, sizeof(*ret->le)*ret->ned);
  ret->nodvale = new double[ret->nnve];
  memset(ret->nodvale, 0, sizeof(*ret->nodvale)*ret->nnve);
  nv = new double[ncomp];  

  if (el.propedg)
  {
    id = 0;
    for(i=0; i<el.ned; i++)
    {
      if (el.propedg[i] == prop)
      {
        ret->le[i] = lgcs;
        for(j=0; j<el.nned[i]; j++)
        {
          if (getval(nodes[el.nodes[el.edgenod[i][j]]], nv))
          {
            delete ret;
            delete [] nv;
            return NULL;
          }
          
          for(k=0; k<ncomp; k++)
          {
            ret->nodvale[id] += nv[k];
            id++;
          }
        }
      }
      else
        id += el.nned[i]*ncomp;
    }
  }
  
  if (edgid && nedgid)
  {
    // there are found edge indeces comming from additional edges
    id = 0;
    l = 0;  
    for(i=0; i<el.ned; i++)
    {
      if (i == edgid[l])
      {
        ret->le[i] = lgcs;
        for(j=0; j<el.nned[i]; j++)
        {
          if (getval(nodes[el.nodes[el.edgenod[i][j]]], nv))
          {
            delete ret;
            delete [] nv;
            return NULL;
          }
        
          for(k=0; k<ncomp; k++)
          {
            ret->nodvale[id] += nv[k];
            id++;
          }
        }
        if (l < nedgid-1)
          l++;
        else
          break;
      }
      else
        id += el.nned[i]*ncomp;
    }
  }

  delete [] nv;
  return ret;
}



/**
  The function creates new object of the loadel type 
  which contains surface loads generated for given element and 
  its surfaces with given property.

  @param el - element whose surface load will be generated
  @param prop - surface property id of the element el where the surface load will be applied
  @param nodes - array with node structures
  @param surfid   - array with surface indeces for the given property number prop detected on additional surfaces
  @param nsurfid - the length of array surf

  @return The function returns pointer to the new allocated object of the loadel type.

  Created by Tomas Koudelka, 06.2009
  Modified by Tomas Koudelka, 08.2011
*/
loadel *entityload::surface2loadel(selement &el, long prop, snode *nodes, long *surfid, long nsurfid)
{
  long i, j, k, l, id, tnns;
  loadel *ret=new loadel;
  double *nv;

  for (i=0, tnns=0; i<el.nsurf; i++)
    tnns += el.nnsurf[i];
  ret->tel    = surface;
  ret->nsurf  = el.nsurf;
  ret->nnvs   = tnns*ncomp;
  ret->nlc    = nlc;
  ret->nslc   = nslc;
  ret->ls = new long[ret->nsurf];
  memset(ret->ls, 0, sizeof(*ret->ls)*ret->nsurf);
  ret->nodvals = new double[ret->nnvs];
  memset(ret->nodvals, 0, sizeof(*ret->nodvals)*ret->nnvs);
  nv = new double[ncomp];  

  if (el.propsurf)
  {
    id = 0;
    for(i=0; i<el.nsurf; i++)
    {
      if (el.propsurf[i] == prop)
      {
        ret->ls[i] = lgcs;
        for(j=0; j<el.nnsurf[i]; j++)
        {
          if (getval(nodes[el.nodes[el.surfnod[i][j]]], nv))
          {
            delete ret;
            delete [] nv;
            return NULL;
          }
        
          for(k=0; k<ncomp; k++)
          {
            ret->nodvals[id] += nv[k];
            id++;
          }
        }
      }
      else
        id += el.nnsurf[i]*ncomp;
    }
  }

  if (surfid && nsurfid)
  {
    // there are found surface indeces comming from additional surfaces
    id = 0;
    l = 0;  
    for(i=0; i<el.nsurf; i++)
    {
      if (i == surfid[l])
      {
        l++;
        ret->ls[i] = lgcs;
        for(j=0; j<el.nnsurf[i]; j++)
        {
          if (getval(nodes[el.nodes[el.surfnod[i][j]]], nv))
          {
            delete ret;
            delete [] nv;
            return NULL;
          }
        
          for(k=0; k<ncomp; k++)
          {
            ret->nodvals[id] += nv[k];
            id++;
          }
        }
        if (l < nsurfid-1)
          l++;
        else
          break;
      }
      else
        id += el.nnsurf[i]*ncomp;
    }
  }

  delete [] nv;
  return ret;
}



/**
  The function creates new object of the loadel type 
  which contains volume loads generated for given element and 
  its volume with given property.

  @param el - element whose volume load will be generated
  @param prop - surface property id of the element el where the volume load will be applied
  @param nodes - array with node structures

  @return The function returns pointer to the new allocated object of the loadel type.

  Created by Tomas Koudelka, 06.2009
*/
loadel *entityload::vol2loadel(selement &el, long prop, snode *nodes)
{
  long i, k, id;
  loadel *ret=new loadel;
  double *nv;

  nv = new double[ncomp];  
  if (el.prop == prop)
  {
    ret->nnvv = el.nne*ncomp;
    ret->tel    = volume;
    ret->nlc    = nlc;
    ret->nslc   = nslc;
    ret->nodvalv = new double[ret->nnvv];
    memset(ret->nodvalv, 0, sizeof(*ret->nodvalv)*ret->nnvv);
     
    id = 0;
    for (i=0; i<el.nne; i++)       
    {
      if (getval(nodes[el.nodes[i]], nv))
      {
        delete ret;
        delete [] nv;
        return NULL;
      }
      for(k=0; k<ncomp; k++)
      {
        ret->nodvalv[id] += nv[k];
        id++;
      }
    }
  }
  delete [] nv;
  return ret;
}
