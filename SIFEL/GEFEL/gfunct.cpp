#include "gfunct.h"
#include "tablefunct.h"
#include "itablefunct.h"
#include "intools.h"
#include "vector.h"
#include <string.h>



/**
  This constructor initializes attributes to zero values.
*/
gfunct::gfunct()
{
  //  number of expression sets
  neqs = 0L;
  //  number of gfunction sets
  ngf = 0L;
  //  constant value
  f = 0.0;
  tfunc=constant;
  eq = NULL;
  deq = NULL;
  var = NULL;
  func = NULL;
  limval = NULL;
  tabf = NULL;
  itabf = NULL;
  gfs = NULL;
}

/**
   second constructor
   
   @param tf - type of general function
   @param nr - the number of rows in the table, otherwise it is ignored
*/
gfunct::gfunct(generalfunct tf,long nr)
{
  //  number of expression sets
  neqs = 0L;
  //  number of gfunction sets
  ngf = 0L;
  //  constant value
  f = 0.0;
  
  tfunc=tf;
  
  eq = NULL;
  deq = NULL;
  var = NULL;
  func = NULL;
  limval = NULL;
  itabf = NULL;
  gfs = NULL;
  
  switch (tfunc){
  case constant:{
    
    break;
  }
  case tab_file:
  case tab:{
    tabf = new tablefunct (nr);
    break;
  }
  default:{
    print_err("unknown type of general function is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  Destructor deallocates memory used by the expression
*/
gfunct::~gfunct()
{
  long i;
  for (i = 0; i < neqs; i++){
    delete [] func[i];
    delete [] var[i];
    delete eq[i];
    delete deq[i];
  }
  delete [] func;
  delete [] eq;
  delete [] deq;
  delete [] var;
  delete [] gfs;
  delete tabf;
  delete itabf;
}



/**
  This function reads data from the text file given by the parameter in. The first it reads
  type of function see galias.h - timefunct. Then it depending on the function type it
  reads data about function.

  Parameters:
  @param in - pointer to opened text file from which data will be read.

  Returns:
  @retval 0 - on succes
  @retval 1 - error parsing expression in case pars function type
  @retval 2 - number of expressions in case pars_set function type is <= 0
  @retval 3 - error parsing expression in case pars_set function type
  @retval 4 - number of functions in case gfunc_set function type is <= 0
  @retval 5 - unknown type of function
*/
long gfunct::read(XFILE *in)
{
  long i, j;
  Parser par;

  //  constant=0
  //  pars=1
  //  tab=2
  //  pars_set=3
  //  tab_file=4
  //  multpv=10
  //  itab=20
  //  gfunc_set=30
  xfscanf (in,"%k%m","funct_type", &generalfunct_kwdset, (int*)&tfunc);

  switch (tfunc){
    case constant:{
      xfscanf (in,"%k%lf","const_val", &f);
      break;
    }
    case pars:{
      neqs = 1;
      eq   = new Equation* [neqs];
      deq = new Equation* [neqs];
      var  = new variable**[neqs];
      func = new char*     [neqs];
      memset(deq, 0, sizeof(*deq)*neqs);
      for (i = 0; i < neqs; i++)
      {
        func[i]=new char [256];
        xfscanf(in, "%k%255s", "func_formula", func[i]);
        //fgets(func[i], 255, in->file);
        eq[i] = par.TextToTree(func[i]);
        if (eq[i] == NULL)
	  {
	    print_err("Error parsing expression",__FILE__,__LINE__,__func__);
	    return (1);
	  }
        var[i] = NULL;
        if (eq[i]->Variables.count())
          var[i] = new variable* [eq[i]->Variables.count()];
        for(j=0; j < eq[i]->Variables.count(); j++)
          var[i][j] = eq[i]->Variables.at(j);
        if (eq[0]->Variables.count() == 1)
        {
          deq[0] = new Equation;
          eq[0]->Differentiate(deq[0], 1, &par);
        }
      }
      break;
    }
    case tab:{
      tabf = new tablefunct;
      tabf->read(in);
      break;
    }
    case tab_file:{
      tabf = new tablefunct;
      tabf->read_data_file(in);
      break;
    }
    case pars_set:{
      xfscanf (in, "%k%ld", "num_funct", &neqs);
      if (neqs <= 0)
      {
        print_err("Number of expressions in set <= 0",__FILE__,__LINE__,__func__);
        return (2);
      }
      eq     = new Equation* [neqs];
      deq    = new Equation* [neqs];
      var    = new variable**[neqs];
      func   = new char*     [neqs];
      limval = new double    [neqs];
      memset(deq, 0, sizeof(*deq)*neqs);
      for (i = 0; i < neqs; i++)
      {
        xfscanf(in, "%k%le", "limval", limval+i);
        func[i]=new char [256];
        xfscanf(in, "%k%255s", "func_formula", func[i]);
        eq[i] = par.TextToTree(func[i]);
        if (eq[i] == NULL)
        {
          print_err("parsing expression", __FILE__, __LINE__, __func__);
          return (3);
        }
        var[i] = NULL;
        if (eq[i]->Variables.count())
          var[i] = new variable* [eq[i]->Variables.count()];
        for(j=0; j < eq[i]->Variables.count(); j++)
          var[i][j] = eq[i]->Variables.at(j);
        if (eq[i]->Variables.count() == 1)
        {
          deq[i] = new Equation;
          eq[i]->Differentiate(deq[i], 1, &par);
        }
      }
      break;
    }
    case itab:{
      itabf = new itablefunct;
      itabf->read(in);
      break;
    }
    case gfunc_set:{
      xfscanf (in, "%k%ld", "num_funct", &ngf);    
      if (ngf <= 0)
      {
        print_err("Number of functions in set <= 0",__FILE__,__LINE__,__func__);
        return (4);
      }
      gfs = new gfunct[ngf];
      for (i = 0; i < ngf; i++)
      {
        gfs[i].read(in);
      }
      break;
    }
    default:{
      print_err("unknown type of general function is required", __FILE__, __LINE__, __func__);
      return (5);
    }
  }
  return(0);
}



/**
  This function computes value y of given funtion expression for variable x whose
  value is given by the t parameter.
  
  Parameters:
  @param t - given value of function variable

  Returns:
  Function value at parameter t.
*/
double gfunct::getval(double t)
{
  long i;
  double ret = 0.0;

  switch (tfunc)
  {
    case constant:
      ret=f;
      break;
    case pars:
      if (var[0])
        var[0][0]->Value = t;
      ret = eq[0]->Evaluate();
      break;
    case tab_file:
    case tab:
      ret=tabf[0].getval(t);
      break;
    case pars_set:
    {
      for (i = 0; i < neqs; i++)
      {
        if (limval[i] >= t)
        {
          if (var[i])
            var[i][0]->Value = t;
          ret = eq[i]->Evaluate();
          break;
        }
        if (i == neqs-1)
        {
          if (var[i])
            var[i][0]->Value = t;
          ret = eq[i]->Evaluate();
        }
      }
      break;
    }
    case gfunc_set:
    {
      for (i = 0; i < ngf; i++)
      {
        ret += gfs[i].getval(t);
      }
      break;
    }
    default:
    {
      print_err("unknown type of general function is required - %d, this=%p", __FILE__, __LINE__, __func__, tfunc, this);
    }
  }
  
  return(ret);
}

/**
  This function computes value x of given funtion expression for variable y whose
  value is given by the t parameter.
  
  Parameters:
  @param t - given value of function variable

  Returns:
  Function value at parameter t.
*/
double gfunct::getinvval(double t)
{
  double ret = 0.0;
  
  switch (tfunc){
  case tab:{
    ret=tabf[0].getinvval(t);
    break;
  }
  default:{
    print_err("unknown type of general function is required - %d, this=%p", __FILE__, __LINE__, __func__, tfunc, this);
  }
  }
  
  return(ret);
}


/**
  This function computes value of given funtion expression for variable
  %vector p. If number p components is greater than 1 only constant and pars types can 
  be handled properly.

  Parameters:
  @param p - %vector of variable values (the order is x,y,z,t usually)
  @param namevar - array of strings with variable names for the each p component 
                   namevar[i] = pointer to the string with the name of the i-th variable

  Returns:
  Function value at parameter p.
*/
double gfunct::getval(vector &p, const char *namevar[])
{
  long i, j;
  double ret = 0.0;

  if (p.n == 1)
  {
    ret = getval(p(0)); 
    return ret;
  }
  switch (tfunc)
  {
    case constant:
      ret = f;
      break;
    case tab: // time is supposed to be the last component of p
      ret=tabf->getval(p(p.n-1));
      break;
    case pars:
      for(i=0; i<p.n; i++)
      {
        for (j=0; j < eq[0]->Variables.count(); j++)
        {
          if (strcmp(namevar[i], var[0][j]->Name) == 0)
          {
            var[0][j]->Value = p(i);
            break;
          }
        }
      }

      ret = eq[0]->Evaluate();
      break;
    default:
      print_err("unknown type of general function is required", __FILE__, __LINE__, __func__);
  }

  return ret;
}



/**
  This function computes value of given funtion expression for variables whose
  values are given by the t and t0 parameter.

  Parameters:
  @param t - given value of the first variable of function
  @param t0  - given value of the second auxiliar variable of function

  Returns:
  Function value at parameter t.
*/
double gfunct::getval(double t, double t0)
{
  long i;
  double ret = 0.0;

 
  switch (tfunc)
  {
    case constant:
      ret=f;
      break;
    case pars:
      if (eq[0]->Variables.count() < 2)
      {
        getval(t);
        break;
      }
      if (var[0])
      {
        var[0][0]->Value = t;
        if (var[0][1])
          var[0][1]->Value = t0;
      }
      ret = eq[0]->Evaluate();
      break;
    case tab_file:
    case tab:
      ret=tabf[0].getval(t);
      break;
    case pars_set:
    {
      for (i = 0; i < neqs; i++)
      {
        if (limval[i] >= t)
        {
          if (eq[i]->Variables.count() < 2)
          {
            getval(t);
            break;
          }
          if (var[i])
          {
            var[i][0]->Value = t;
            if (var[i][1])
              var[i][1]->Value = t0;
          }
          ret = eq[i]->Evaluate();
          break;
        }
        if (i == neqs-1)
        {
          if (eq[i]->Variables.count() < 2)
          {
            getval(t);
            break;
          }
          if (var[i])
          {
            var[i][0]->Value = t;
            if (var[i][1])
              var[i][1]->Value = t;
          }
          ret = eq[i]->Evaluate();
        }
      }
      break;
    }
    case gfunc_set:
    {
      for (i = 0; i < ngf; i++)
      {
        ret += gfs[i].getval(t,t0);
      }
      break;
    }
    default:
    {
      print_err("unknown type of general function is required", __FILE__, __LINE__, __func__);
    }
  }

  return(ret);
}



/**
  This function computes value of given funtion expression for variable whose
  value is given by the t parameter.

  Parameters:
  @param t - given value of function variable

  Returns:
  Function integer value at parameter t.
*/
long gfunct::getval_long(double t)
{
  long i, ret=0;
  
  switch (tfunc)
  {
    case itab:{
      ret=itabf[0].getval(t);
      break;
    }
    case gfunc_set:
    {
      ret = 0;
      for (i = 0; i < ngf; i++)
      {
        ret += gfs[i].getval_long(t);
      }
      break;
    }
    default:{
      print_err("unknown type of general function is required", __FILE__, __LINE__, __func__);
    }
  }
  
  return(ret);
}

/**
  This function computes derivative of given funtion expression for variable whose
  value is given by the t parameter.

  Parameters:
  @param t - given value of function variable

  Returns:
  derivative of the function at parameter t.
  
  JK, 7.10.2011
*/
double gfunct::getderiv (double t)
{
  double ret = 0.0;

  switch (tfunc){
  case constant:{
    ret=0.0;
    break;
  }
  case tab_file:
  case tab:{
    ret=tabf[0].getderiv (t);
    break;
  }
  case pars:
    if (deq[0])
    {
      if (var[0][0])
        var[0][0]->Value = t;
      ret = deq[0]->Evaluate();
    }
    break;
  case pars_set:
    for (long i = 0; i < neqs; i++) 
    {
      if (limval[i] >= t)
      {
        if (deq[i])
        {
          var[i][0]->Value = t;
          ret = deq[i]->Evaluate();
          break;
        }
      }
      if (i == neqs-1)
      {
        if (deq[i])
        {
          var[i][0]->Value = t;
          ret = deq[i]->Evaluate();
        }
      }
    }
    break;
  default:{
    print_err("unknown type of general function is required", __FILE__, __LINE__, __func__);
  }
  }
  
  return(ret);
}

/**
  This function computes derivative of given funtion expression for variable whose
  value is given by the t parameter.

  Parameters:
  @param t - given value of function variable

  Returns:
  derivative of the function at parameter t.
  
  JK, 7.10.2011
*/
double gfunct::getinvderiv (double t)
{
  double ret = 0.0;

  switch (tfunc){
  case constant:{
    ret=0.0;
    break;
  }
  case tab_file:
  case tab:{
    ret=tabf[0].inv_derivative (t);
    break;
  }
  default:{
    print_err("unknown type of general function is required", __FILE__, __LINE__, __func__);
  }
  }
  
  return(ret);
}



/**
  This function reads data from the text file given by the parameter in. It is
  used in the mechprep preprocessor.

  Parameters:
  @param in - pointer to opened text file from which data will be read.

  Returns:
  @retval 0 - on succes
  @retval 1 - error parsing expression in case pars function type
  @retval 2 - number of expression in case pars_set function type is <= 0
  @retval 3 - error parsing expression in case pars_set function type
  @retval 4 - number of functions in case gfunc_set function type is <= 0
  @retval 5 - unknown type of function
*/
long gfunct::read_prop(FILE *in)
{
  long i;
  Parser par;

  getint(in, (int &)tfunc);
  switch (tfunc)
  {
    case constant:
      getdouble(in, f);
      break;
    case pars:
      neqs = 1;
      eq   = new Equation* [neqs];
      func = new char*     [neqs];
      for (i = 0; i < neqs; i++)
      {
        func[i]=new char [256];
        getstring(in, func[i], 255);
        eq[i] = par.TextToTree(func[i]);
        if (eq[i] == NULL)
	  {
	    print_err("Error parsing expression", __FILE__, __LINE__, __func__);
	    return(1);
	  }
      }
      break;
    case tab_file:
    case tab:
      tabf = new tablefunct;
      tabf->read_prop(in);
      break;
    case pars_set:
    {
      getlong(in, neqs);
      if (neqs <= 0)
      {
	print_err("Number of expressions in set <= 0", __FILE__, __LINE__, __func__);
        return (2);
      }
      eq     = new Equation* [neqs];
      func   = new char*     [neqs];
      limval = new double    [neqs];
      for (i = 0; i < neqs; i++)
      {
        getdouble(in, limval[i]);
        func[i]=new char [256];
        getstring(in, func[i], 255);
        eq[i] = par.TextToTree(func[i]);
        if (eq[i] == NULL)
        {
	  print_err("Error parsing expression", __FILE__, __LINE__, __func__);
          return(3);
        }
      }
      break;
    }
    case itab:
      itabf = new itablefunct;
      itabf->read_prop(in);
      break;
    case gfunc_set:
    {
      getlong(in, ngf);
      if (ngf <= 0)
      {
	print_err("Number of functions in set <= 0", __FILE__, __LINE__, __func__);
        return (4);
      }
      for (i = 0; i < ngf; i++)
      {
        gfs[i].read_prop(in);
      }
      break;
    }
    default:
      print_err("unknown type of general function is required", __FILE__, __LINE__, __func__);
      return(5);
  }
  return(0);
}



/**
  This function writes data to the text file given by the parameter in. It is
  used in the mechprep preprocessor.

  Parameters:
  @param out - pointer to opened text file from which data will be written.

  Returns:
  @retval 0 - on succes
  @retval 1 - unknown type of function
*/
long gfunct::print(FILE *out)
{
  Parser par;
  long i;

  fprintf (out,"\n%d ",tfunc);

  switch (tfunc){
  case constant:{
    fprintf (out,"%e\n",f);
    break;
  }
  case pars:{
    fprintf(out, " %s\n", func[0]);
    break;
  }
  case tab:{
    tabf->print(out);
    break;
  }
  case tab_file:{
    tabf->print_data_file(out);
    break;
  }
  case pars_set:{
    fprintf (out, "%ld\n", neqs);
    for (i=0; i<neqs; i++)
      fprintf(out, " %s\n", func[i]);
    break;
  }
  case itab:{
    itabf->print(out);
    break;
  }
  case gfunc_set:{
    fprintf (out, "%ld\n", ngf);
    for (i = 0; i < ngf; i++)
    {
      gfs[i].print(out);
    }
    break;
  }
  default:{
    print_err("unknown type of general function is required", __FILE__, __LINE__, __func__);
    return (1);
  }
  }
  fprintf (out,"\n");

  return(0);
}



/**
   function initiates general function from gf
   
   JK
*/
void gfunct::initiate (gfunct &gf)
{
  long i;

  tfunc = gf.tfunc;
  
  switch (tfunc){
  case constant:{
    f=gf.f;
    break;
  }
  case tab_file:
  case tab:{
    tabf = new tablefunct;
    tabf->itype = gf.tabf->itype;
    tabf->asize = gf.tabf->asize;
    tabf->x=new double[tabf->asize];
    tabf->y=new double[tabf->asize];
    for (i=0;i<tabf->asize;i++){
      tabf->x[i]=gf.tabf->x[i];
      tabf->y[i]=gf.tabf->y[i];
    }
    break;
  }
  default:{
    print_err("unknown type of general function is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
   The function initiates general function of table type.

   Parameters:   
   @param nr - number of rows

   Returns:
   @retval The function does not return anything.
   
   JK, 13.6.2005
*/
void gfunct::init_tab (long nr)
{
  tfunc = tab;
  tabf = new tablefunct;
  tabf->itype = piecewiselin;
  tabf->asize = nr;
  tabf->x = new double [nr];
  tabf->y = new double [nr];
}



/**
   The function copies all data from subroutine argument to the actual object.

   Parameters:
   @param gf - function whose data will be copied to the actual object
   
   Returns:
   @retval The function does not return anything.

   JK, 15.11.2008
   TKo, 6.2009
*/
void gfunct::copy (gfunct &gf)
{
  long i, j;
  Parser par;

  switch (gf.tfunc)
  {
    case constant:{
      //  type of general function
      tfunc = gf.tfunc;
      // number of expression sets
      neqs = gf.neqs;
      // number of general functions
      ngf = gf.ngf;
      //  constant value
      f = gf.f;
      break;
    }
    case pars:
    case pars_set:{
      if (neqs)
      {
        for (i=0; i<neqs; i++)
        {
          delete eq[i];
          delete deq[i];
          delete [] var[i];
          delete [] func[i];
        }
        delete [] eq;
        delete [] deq;
        delete [] var;
        delete [] func;
      }
      //  type of general function
      tfunc = gf.tfunc;
      // number of expression sets
      neqs = gf.neqs;
      // number of general functions
      ngf = gf.ngf;
      eq   = new Equation* [gf.neqs];
      deq  = new Equation* [gf.neqs];
      memset(deq, 0, sizeof(*deq)*gf.neqs);
      var  = new variable**[gf.neqs];
      func = new char*     [gf.neqs];
      for (i = 0; i < gf.neqs; i++)
      {
        func[i]=new char [256];
        strcpy(func[i], gf.func[i]);
        eq[i] = par.TextToTree(func[i]);
        var[i] = NULL;
        if (eq[i]->Variables.count())
          var[i] = new variable* [eq[i]->Variables.count()];
        for(j=0; j < eq[i]->Variables.count(); j++)
          var[i][j] = eq[i]->Variables.at(j);
        if (eq[i]->Variables.count() == 1)
        {
          deq[i] = new Equation;
          eq[i]->Differentiate(deq[i], 1, &par);
        }
      }
      break;
    }
    case tab_file:
    case tab:{
      //  type of general function
      tfunc = gf.tfunc;
      // number of expression sets
      neqs = gf.neqs;
      // number of general functions
      ngf = gf.ngf;
      if (tabf)
        delete tabf;
      tabf = new tablefunct;
      tabf->copy (gf.tabf[0]);
      break;
    }
    case itab:{
      //  type of general function
      tfunc = gf.tfunc;
      // number of expression sets
      neqs = gf.neqs;
      // number of general functions
      ngf = gf.ngf;
      if (itabf)
        delete itabf;
      itabf = new itablefunct;
      itabf->copy (gf.itabf[0]);
      break;
    }
    case gfunc_set:
    {
      //  type of general function
      tfunc = gf.tfunc;
      // number of expression sets
      neqs = gf.neqs;
      // number of general functions
      ngf = gf.ngf;     
      if (ngf)
        delete [] gfs;
      gfs = new gfunct[ngf];
      for (i = 0; i < gf.ngf; i++)
      {
        gfs[i].copy(gf.gfs[i]);
      }
      break;
    }
    default:{
      print_err("unknown type of general function is required",__FILE__,__LINE__,__func__);
    }  
  }  
  return;
}



/**
   Function merges all data from subroutine argument gf with the actual object.

   Parameters:
   @param gf - function whose data will be merged with the actual object

   Returns:
   @retval The function does not return anything.
   
   TKo, 6.2009
*/
void gfunct::merge (gfunct &gf)
{
  long i;
  gfunct *tmp;

  if ((tfunc == gfunc_set) && (gf.tfunc != gfunc_set))
  { // actual object is set of general functions and merged object is basic type (i.e. parser, table, ...)
    tmp = new gfunct[ngf+1];
    for(i=0; i<ngf; i++)
      tmp[i].copy(gfs[i]);
    tmp[ngf].copy(gf);
    ngf++;
    delete [] gfs;
    // number of general functions
    gfs = tmp;
    return;
  }
  if ((gf.tfunc == gfunc_set) && (tfunc != gfunc_set))
  {// actual object is is basic type (i.e. parser, table, ...) and merged object is set of general functions
    for (i = 0; i < neqs; i++){
      delete [] func[i];
      delete [] var[i];
      delete eq[i];
    }
    delete [] func;
    delete [] eq;
    delete [] var;
    delete [] gfs;
    delete tabf;
    delete itabf;
    gfs = new gfunct[gf.ngf+1];
    gfs[0].copy(*this);
    for(i=0; i<gf.ngf; i++)
      gfs[i+1].copy(gf.gfs[i]);
    //  type of general function
    tfunc = gfunc_set;    
    // number of general functions
    ngf = gf.ngf+1;
    return;
  }
  if ((gf.tfunc == gfunc_set) && (tfunc == gfunc_set))
  {// actual object and merged object are sets of general functions
    tmp = new gfunct[ngf+gf.ngf];
    for(i=0; i<ngf; i++)
      tmp[i].copy(gfs[i]);
    for(i=0; i<gf.ngf; i++)
      tmp[i+ngf].copy(gf.gfs[i]);
    // number of general functions
    ngf+=gf.ngf;
    delete [] gfs;
    gfs = tmp;
    return;
  }
  switch (tfunc){// actual object and merged object are basic types (i.e. parser, table, ...)
    case constant:
    case pars:
    case pars_set:
    case tab_file:
    case tab:
    case itab:
      gfs = new gfunct[2];
      gfs[0].copy(*this); 
      gfs[1].copy(gf); 
      //  type of general function
      tfunc = gfunc_set;
      // number of general functions
      ngf = 2;
      break;
    default:{
      print_err("unknown type of general function is required",__FILE__,__LINE__,__func__);
    }
  }
  return;
}



/**
   The function compares all data from subroutine argument with the actual object.

   Parameters:
   @param gf - function whose data will be compared with the actual object
   
   Returns:
   @retval 0 - in case that objects represents the same general function
   @retval 1 - in case that objects are different or unknown type of general function

   TKo, 09.2010
*/
long gfunct::compare (gfunct &gf)
{
  long i;
  Parser par;


  // check for the same type of gfunct
  if (tfunc != gf.tfunc)
    return 1;

  switch (gf.tfunc)
  {
    case constant:{
      //  constant value
      if (f == gf.f)
        return 0;
      else
        return 1;
    }
    case pars:
      // check for number of expression sets
      if (neqs != gf.neqs)
        return 1;
      if (strcmp(func[0], gf.func[0]))
        return 1;
      // all function strings are mutually identical
      return 0;
    case pars_set:{
      // check for number of expression sets
      if (neqs != gf.neqs)
        return 1;
      for (i = 0; i < gf.neqs; i++)
      {
        if (strcmp(func[i], gf.func[i]))
          return 1;
        if (limval[i] != gf.limval[i])
          return 1;
      }
      // all function strings are mutually identical
      return 0;
    }
    case tab_file:
    case tab:{
      return tabf->compare(gf.tabf[0]);
    }
    case itab:{
      return itabf->compare(gf.itabf[0]);
    }
    case gfunc_set:
    {
      // check for number of general functions
      if (ngf != gf.ngf)
        return 1;
      for (i = 0; i < gf.ngf; i++)
      {
        if (gfs[i].compare(gf.gfs[i]))
         return 1;
      }
      return 0;
    }
    default:{
      print_err("unknown type of general function is required",__FILE__,__LINE__,__func__);
    }  
  }  
  return 1;
}


