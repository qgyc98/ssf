#include "entitytdload.h"
#include "iotools.h"
#include "mathem.h"
#include "gfunct.h"
#include <math.h>
#include <string.h>
#include <ctype.h>



/**
  This constructor initializes data to the zero values
*/
entitytdload::entitytdload()
{
  ncomp = 0L;
  nlc   = 0L;
  lgcs  = 0L;
  func  = NULL;
}



/**
  This destructor deallocates used memory
*/
entitytdload::~entitytdload()
{
  delete [] func;
}



/**
  This function reads data about entity load from the text file given by the parameter in

  @param in - pointer to the opned text file
  @param lc - total number of load cases

  @retval 0 - on succes
  @retval 1 - unsupported function type is required
  @retval 2 - invalid type of load coordinate system
  @retval 3 - function has too much number of parameters
  @retval 4 - wrong name of function parameters
*/
long entitytdload::read(XFILE *in, long lc)
{
  long nspatc;

  xfscanf(in, "%k%ld", "lc_id", &nlc);
  if ((nlc < 1) || (nlc > lc))
  {
    print_err("Line %ld, col %ld, file %s: invalid load case index %ld.\n Valid range is <1,%ld>", 
              __FILE__, __LINE__, __func__, in->line, in->col, in->fname, nlc, lc);
    return(1);
  }
  long i;
  xfscanf(in, "%k%ld", "ncomp", &ncomp);
  xfscanf(in, "%k%ld", "coord_sys", &lgcs);
  if ((lgcs < 1) || (lgcs > 2))
  {
    print_err("Line %ld, col %ld, file %s: unknown type of coordinate system %ld.\n Valid values are 1 or 2", 
              __FILE__, __LINE__, __func__, in->line, in->col, in->fname, lgcs);
    return 2;
  }
  xfscanf(in, "%k", "load_comp");
  func = new gfunct[ncomp];
  for(i=0; i<ncomp; i++)
  {
    func[i].read(in);
    switch(func[i].tfunc)
    {
      case constant:
      case tab:
      case tab_file:
        break;
      case pars:
        if (func[i].eq[0]->Variables.count() > 4)
        {
          print_err("Line %ld, col %ld, file %s: number of parameters %ld of load distribution function > 4", 
                    __FILE__, __LINE__, __func__, in->line, in->col, in->fname, func[i].eq[0]->Variables.count());
          return 3;
        }
        if (var2coord(func[i].eq[0], nspatc))
        {
          print_err("Line %ld, col %ld, file %s: wrong names of parameters of load distribution function."
                    "\n Valid names are x, y, z or t.", 
                    __FILE__, __LINE__, __func__, in->line, in->col, in->fname);
          return 4;
        }
        break;
      default:
        print_err("Line %ld, col %ld, file %s: unknown type %d of load distribution function", 
                  __FILE__, __LINE__, __func__, in->line, in->col, in->fname, int(func[i].tfunc));
        return 5;
    }
  }
  return(0);
}



/**
  The function checks the function defined eq for proper variable names
  that may be one of the following characters <x,y,z,t> (case does not matter).
  The number of different spatial coordinates used in the definition of function eq is stored in the
  parameter nspatc (if coordinates x and z are used in the function definition then nspatc = 2).

  @param eq - pointer to the parsed function
  @param nspatc - number of different spatial coordinates used in the function definition (output)
  
  @retval 0 - on success
  @retval 1 - wrong variable names used
*/
long entitytdload::var2coord(Equation *eq, long &nspatc)
{
  variable *var;
  long i, var_det = 0L;
  char vn;
  nspatc = 0L;

  for (i=0; i<eq->Variables.count(); i++)
  {
    var = eq->Variables.at(i);
    vn = tolower(var->Name[0]);
    if (vn == 'x')
    {
      nspatc++;
      var_det++;
    }
    if (vn == 'y')
    {
      nspatc++;
      var_det++;
    }
    if (vn == 'z')
    {
      nspatc++;
      var_det++;
    }
    if (vn == 't')
      var_det++;
  }
  if (eq->Variables.count() != var_det)
    return 1;
  return 0;
}   



/**
  This function approximates values on the loaded entity to the element node which
  lay on it. It is used linear approximation.

  @param tn - node where the load should be approximated
  @param nv - array with time functions created form approximated nodal values

  @retval 0 - on succes
  @retval 1 - in case of unknown type of load function
*/
long entitytdload::getval(snode &tn, gfunct *nv)
{
  long i, j, noc, nspatc = 0;
  char *cp;
  char *psrc, *pdest;
  double coord;
  Parser par;
    

  for (i=0; i<ncomp; i++)
  {
    switch (nv[i].tfunc)
    {
      case constant:
      case tab:
        nv[i].copy(func[i]);
        break;
      case pars:
      {
        var2coord(func[i].eq[0], nspatc);
        // setting of spatial coordinates in the parser string
        if (nspatc)
        {
          // calculate the total number of occurences of spatial coordinates
          noc = 0;
          while ((cp = strpbrk(func[i].func[0], "xXyYzZ")))
            noc++;
          // create new parser instance in the i-th load component
          nv[i].neqs = 1;
          nv[i].eq   = new Equation* [nv[i].neqs];
          nv[i].var  = new variable**[nv[i].neqs];
          nv[i].func = new char*     [nv[i].neqs];
          nv[i].func[0] = new char[256+noc*25];         
          // replace occurences in the parser string by real spatial coordinates of node tn
          cp = strpbrk(func[i].func[0], "xXyYzZ");
          psrc = func[i].func[0];   // set pointer to the parser string part preceding the coordinate occurence
          pdest = nv[i].func[0]; // set pointer to the resulting parser string where the output should be printed out
          while (cp)
          {
            // set value of the actual coordinate occurence
            if (tolower(*cp) == 'x')   coord = tn.x;
            if (tolower(*cp) == 'y')   coord = tn.y;
            if (tolower(*cp) == 'z')   coord = tn.z;
            // copy parser string before the actual coordinate occurence
            memcpy(pdest, psrc, sizeof(*psrc)*(cp-psrc));
            // set position of the resulting parser string for the printing the actual cooridnate 
            pdest += cp-psrc;
            // print the value of the nodal coordinate for the given occurence
            pdest += sprintf(pdest, "(%.15le)", coord);
            // actualize pointer to the beginning of searched part of the original parser string
            psrc = cp + 1;
            cp = strpbrk(psrc, "xXyYzZ");
          }
          // convert string to equation tree
          nv[i].eq[0] = par.TextToTree(nv[i].func[0]);
          if (nv[i].eq[0] == NULL)
	  {
	    print_err("Error parsing expression of %ld. load component",__FILE__,__LINE__,__func__, i+1);
	    return (1);
	  }
          // Compress parser string with the help of TreeToText conversion performed in the string parsing.
          // It is supposed that equation tree evaluates operations with real numbers that replaced spatial coordinates
          // and thus the resulting string will be compressed.
          strcpy(nv[i].func[0], nv[i].eq[0]->EquationText);
          // prepare variable pointers
          nv[i].var[0] = NULL;
          if (nv[i].eq[0]->Variables.count())
            nv[i].var[0] = new variable* [nv[i].eq[0]->Variables.count()];
          for(j=0; j < nv[i].eq[0]->Variables.count(); j++)
            nv[i].var[0][j] = nv[i].eq[0]->Variables.at(j);
        }
        else
          nv[i].copy(func[i]);
        break;
      }
      default:
        return 1;
    }
  }
  return(0);
}



/**
  The function creates new object of the dloadel type 
  which contains edge loads generated for given element and 
  its edges with given property.

  @param el - element whose edge load will be generated
  @param prop - edge property id of the element el where the edge load will be applied
  @param nodes - array with node structures
  @param edgid   - array with edge indeces for the given property number prop detected on additional edges
  @param nedgid - the length of array edgid

  @return The function returns pointer to the new allocated object of the dloadel type.

  Created by Tomas Koudelka, 06.2009
  Modified by Tomas Koudelka, 08.2011
*/
dloadel *entitytdload::tdedge2dloadel(selement &el, long prop, snode *nodes, long *edgid, long nedgid)
{
  long i, j, k, l, id, tnne;
  long *nodvale_set; // array of indicators of set ret->nodvale values (functions gfunct)
  dloadel *ret=new dloadel;
  gfunct *nv;

  for (i=0, tnne=0; i<el.ned; i++)
    tnne += el.nned[i];
  ret->tel  = edge;
  ret->ned  = el.ned;
  ret->nnve = tnne*ncomp;
  ret->nlc  = nlc;
  ret->le = new long[ret->ned];
  memset(ret->le, 0, sizeof(*ret->le)*ret->ned);
  ret->nodvale = new gfunct[ret->nnve];
  nodvale_set = new long[ret->nnve];
  memset(nodvale_set, 0, sizeof(*nodvale_set)*ret->nnve);
  nv = new gfunct[ncomp];  

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
            delete [] nodvale_set;
            return NULL;
          }
          
          for(k=0; k<ncomp; k++)
          {
            ret->nodvale[id].copy(nv[k]);
            nodvale_set[id] = 1;
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
            delete [] nodvale_set;
            return NULL;
          }
        
          for(k=0; k<ncomp; k++)
          {
            if (nodvale_set[id])
              ret->nodvale[id].merge(nv[k]);
            else
            {
              ret->nodvale[id].copy(nv[k]);
              nodvale_set[id] = 1;
            }
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
  delete [] nodvale_set;
  return ret;
}



/**
  The function creates new object of the dloadel type 
  which contains surface loads generated for given element and 
  its surfaces with given property.

  @param el - element whose surface load will be generated
  @param prop - surface property id of the element el where the surface load will be applied
  @param nodes - array with node structures
  @param surfid   - array with surface indeces for the given property number prop detected on additional surfaces
  @param nsurfid - the length of array surf

  @return The function returns pointer to the new allocated object of the dloadel type.

  Created by Tomas Koudelka, 06.2009
  Modified by Tomas Koudelka, 08.2011
*/
dloadel *entitytdload::tdsurface2dloadel(selement &el, long prop, snode *nodes, long *surfid, long nsurfid)
{
  long i, j, k, l, id, tnns;
  long *nodvals_set; // array of indicators of set ret->nodvals values (functions gfunct)
  dloadel *ret=new dloadel;
  gfunct *nv;

  for (i=0, tnns=0; i<el.nsurf; i++)
    tnns += el.nnsurf[i];
  ret->tel    = surface;
  ret->nsurf  = el.nsurf;
  ret->nnvs   = tnns*ncomp;
  ret->nlc    = nlc;
  ret->ls = new long[ret->nsurf];
  memset(ret->ls, 0, sizeof(*ret->ls)*ret->nsurf);
  ret->nodvals = new gfunct[ret->nnvs];
  nodvals_set = new long[ret->nnvs];
  memset(nodvals_set, 0, sizeof(*nodvals_set)*ret->nnvs);
  nv = new gfunct[ncomp];  

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
            delete [] nodvals_set;
            return NULL;
          }
        
          for(k=0; k<ncomp; k++)
          {
            ret->nodvals[id].copy(nv[k]);
            nodvals_set[id] = 1;
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
            delete [] nodvals_set;
            return NULL;
          }
        
          for(k=0; k<ncomp; k++)
          {
            if (nodvals_set[id])
              ret->nodvals[id].merge(nv[k]);
            else
            {
              ret->nodvals[id].copy(nv[k]);
              nodvals_set[id] = 1;
            }
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
  The function creates new object of the dloadel type 
  which contains volume loads generated for given element and 
  its volume with given property.

  @param el - element whose volume load will be generated
  @param prop - surface property id of the element el where the volume load will be applied
  @param nodes - array with node structures

  @return The function returns pointer to the new allocated object of the dloadel type.

  Created by Tomas Koudelka, 06.2009
*/
dloadel *entitytdload::tdvol2dloadel(selement &el, long prop, snode *nodes)
{
  long i, k, id;
  dloadel *ret=new dloadel;
  gfunct *nv;

  nv = new gfunct[ncomp];  
  if (el.prop == prop)
  {
    ret->nnvv = el.nne*ncomp;
    ret->tel    = volume;
    ret->nlc    = nlc;
    ret->nodvalv = new gfunct[ret->nnvv];
     
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
        ret->nodvalv[id].copy(nv[k]);
        id++;
      }
    }
  }
  delete [] nv;
  return ret;
}
