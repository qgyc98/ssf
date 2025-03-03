#include <string.h>
#include "graphmat.h"
#include "global.h"
#include "mechmat.h"
#include "intpoints.h"
#include "tablefunct.h"
#include "matrix.h"
#include "parser.h"



/**
  The constructor inializes attributes to zero or default values.

  Created by Tomas Koudelka
*/
graphmat::graphmat(void)
{
   gt = glinear;
   numf = 0L;
   k = 0.0;
   func = NULL;
   deq = NULL;
   eq = NULL;
   deq = NULL;
   tab = NULL;
}



/**
  The destructor deallocates used memory.

  Created by Tomas Koudelka
*/
graphmat::~graphmat(void)
{
   for (long i = 0; i < numf; i++)
   {
     delete [] func[i];
     delete deq[i];
     delete eq[i];
   }
   delete [] func;
   delete [] deq;
   delete [] eq;
   delete [] limval;
   
   delete tab;
}



/**
  The function reads material parameters from the opened text file given
  by the parameter in. First it reads the graph type attribute which
  denotes type of used input of the diagram. Then the appropriate
  input function is called.

  @param in - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - unable to read graph type
  @retval 2 - error parsing single expression
  @retval 3 - error parsing multiple expression
  @retval 4 - unknown graph type
  
  Created by Tomas Koudelka
*/
long graphmat::read (XFILE *in)
{
  long i;
  Parser par;

  xfscanf(in, "%m", &graphtype_kwdset, (int *)&gt);
  if ((gt < glinear) || (gt > gfunc_ser))
  {
    print_err("unknown type of graph", __FILE__, __LINE__, __func__);
    return(1);
  }

  switch (gt)
  {
    case glinear:
    {
      // stiffness parameter is required
      xfscanf (in,"%lf", &k);
      break;
    }
    case gtable:
    {
      // stress-strain diagram is described by the table
      tab = new tablefunct;
      tab->read(in);
      break;
    }
    case gfunc:
    {
      // stress-strain diagram is described by the expression
      eq = new Equation* [1];
      func = new char* [1];
      func[0] = new char[256];
      memset(func[0], 0, sizeof(*func[0])*256);
      xfscanf(in, "%255s", func[0]);
      eq[0] = par.TextToTree(func[0]);
      if (eq[0] == NULL)
      {
        print_err("parsing expression in function", __FILE__, __LINE__, __func__);
        return(2);
      }
      deq = new Equation* [1];
      deq[0] = new Equation;
      eq[0]->Differentiate(deq[0], 1, &par);
      break;
    }
    case gfunc_ser:
    {
      // stress-strain diagram is described by several different expressions
      // according to attained strains
      xfscanf(in, "%k%ld", "numf", &numf);
      func = new char* [numf];
      memset(func, 0, sizeof(*func)*numf);
      eq  = new Equation* [numf];
      deq = new Equation* [numf];
      limval = new double [numf];
      for (i = 0; i < numf; i++)
      {
        func[i] = new char[256];
        memset(func[i], 0, sizeof(*func[i])*256);
        xfscanf(in, "%le", limval+i); // strain limit value
        xfscanf(in, "%255s", func[i]);   // expression
        eq[i]  = par.TextToTree(func[i]);
        if (eq[i] == NULL)
        {
          print_err("parsing expression in function", __FILE__, __LINE__, __func__);
          return(3);
        }
        deq[i] = new Equation;
        eq[i]->Differentiate(deq[i], 1, &par);
      }
      break;
    }
    default:
    {
      print_err("unknown type of graph is required", __FILE__, __LINE__, __func__);
      return (4);
    }
  }
  return(0);
}



/**
  The function prints material parameters to the opened text file given
  by the parameter out.

  @param out - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - unknown graph type
  
  Created by Tomas Koudelka
*/
long graphmat::print (FILE *out)
{
  long i;

  fprintf(out, "%d ", (int)gt);

  switch (gt)
  {
    case glinear:
    {
      fprintf (out,"%le", k);
      break;
    }
    case gtable:
    {
      tab->print(out);
      break;
    }
    case gfunc:
    {
      fprintf(out, "%s", func[0]);
      break;
    }
    case gfunc_ser:
    {
      fprintf(out, "%ld", numf);
      for (i = 0; i < numf; i++)
      {
        fprintf(out, "\n%le ", limval[i]);
        fprintf(out, "%s", func[i]);
      }
      break;
    }
    default:
    {
      print_err("unknown type of graph is required", __FILE__, __LINE__, __func__);
      return (1);
    }
  }
  return(0);
}



/**
  The function computes material stiffnes matrix.

  @param d - allocated matrix structure for material stiffness %matrix (output)
  @param ipp - integration point number

  @return The function returns assembled stiffness %matrix in the parameter d

  Created by Tomas Koudelka
*/
void graphmat::matstiff (matrix &d, long ipp)
{
  long i, id;
  double val;

  switch (gt)
  {
    case glinear:
      for (i = 0; i < d.m; i++)
        d[i][i] = k;
      break;
    case gtable:
      for (i = 0; i < d.m; i++)
      {
        if (i >= Mm->ip[ipp].ncompstr)
        {
          d[i][i] = 0.0;
          break;
        }
        // d[i][i] = tab->getval2(Mm->ip[ipp].strain[i],k); // TKr version - stiffness is supposed to be returned by tab->getval but stresses cannot be calculated simply in this case (must be integrated)
        tab->getval2(Mm->ip[ipp].strain[i], d[i][i]);  // original version - stress is supposed to be returned by tab->getval
      }
      break;
    case gfunc:
      for (i = 0; i < d.m; i++)
      {
        if (i >= Mm->ip[ipp].ncompstr)
        {
          d[i][i] = 0.0;
          break;
        }
        if (deq[0]->Variables.at(0))
          deq[0]->Variables.at(0)->Value = Mm->ip[ipp].strain[i];
        d[i][i] = deq[0]->Evaluate();
      }
      break;
    case gfunc_ser:
    {
      for (i = 0; i < d.m; i++)
      {
        for (id = 0; id < numf; id++)
        {
          val = Mm->ip[ipp].strain[i];
          if (limval[id] >= val)
          {
            if (i >= Mm->ip[ipp].ncompstr)
            {
              d[i][i] = 0.0;
              break;
            }
            if (deq[id]->Variables.at(0))
              deq[id]->Variables.at(0)->Value = Mm->ip[ipp].strain[i];
            d[i][i] = deq[id]->Evaluate();
            break;
          }
          if (id == numf-1)
          {
            if (deq[id]->Variables.at(0))
              deq[id]->Variables.at(0)->Value = Mm->ip[ipp].strain[i];
            d[i][i] = deq[id]->Evaluate();
          }
        }
      }
      break;
    }
    default:
      print_err("unknown graph type is required", __FILE__, __LINE__, __func__);
  }
  return;
}



/**
  The function computes stresses at given integration point ipp,
  depending on the reached strains.
  The stress of the given integration point is actualized.

  @param ipp - integration point number in the mechmat ip array.
  
  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void graphmat::nlstresses (long ipp)
  //
{
  long i, id, nc=Mm->ip[ipp].ncompstr;

  switch (gt)
  {
    case glinear:
      for (i = 0; i < nc; i++)
        Mm->ip[ipp].stress[i] = k * Mm->ip[ipp].strain[i];
      break;
    case gtable:
      for (i = 0; i < nc; i++)
      {
        // TKr version - stiffness is supposed to be returned by tab->getval
        // stress cannot be calculated in this way, it must be integrated !
        /*
        k = tab->getval2(Mm->ip[ipp].strain[i],val);
        Mm->ip[ipp].stress[i] = k * Mm->ip[ipp].strain[i];
        */
        Mm->ip[ipp].stress[i] = tab->getval(Mm->ip[ipp].strain[i]);
      }
      break;
    case gfunc:
      for (i = 0; i < nc; i++)
      {
        if (eq[0]->Variables.at(0))
          eq[0]->Variables.at(0)->Value = Mm->ip[ipp].strain[i];
        Mm->ip[ipp].stress[i] = eq[0]->Evaluate();
      }
      break;
    case gfunc_ser:
    {
      for (i = 0; i < nc; i++)
      {
        for (id = 0; id < numf; id++)
        {
          if (limval[id] >= Mm->ip[ipp].strain[i])
          {
            if (eq[id]->Variables.at(0))
              eq[id]->Variables.at(0)->Value = Mm->ip[ipp].strain[i];
            Mm->ip[ipp].stress[i] = eq[id]->Evaluate();
            break;
          }
          if (id == numf-1)
          {
            if (eq[id]->Variables.at(0))
              eq[id]->Variables.at(0)->Value = Mm->ip[ipp].strain[i];
            Mm->ip[ipp].stress[i] = eq[id]->Evaluate();
          }
        }
      }
      break;
    }
    default:
      print_err("unknown graph type is required", __FILE__, __LINE__, __func__);
  }
}
