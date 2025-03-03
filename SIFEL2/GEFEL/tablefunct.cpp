/*
    File:                    tablefunct.cpp
    Author:                  Martin Kulhavy
    Purpose:                 Class computes function value at exact point, where the function is given by the table x,y
                             Implemented interpolations: dirac <), lagrange, piecewiselinear 
    Revision history:        Martin Kulhavy 20.11.2001 - file was created
*/

#include <stdlib.h>
#include <string.h>
#include "tablefunct.h"
#include "intools.h"



/**
   The constructor initializes data members to zero values

   created  20.11.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
*/  
tablefunct::tablefunct ()
{
  itype = piecewiselin;
  asize=0;
  x=NULL;
  y=NULL;
  file1 = new char[1000];
  memset(file1, 0, sizeof(*file1)*1000);
  dealloc = true;
}



/**
   The constructor allocates table with n items.

   @param n - the number of table items

   Created  by TKo, 09.2010
*/  
tablefunct::tablefunct (long n)
{
  itype = piecewiselin;
  asize = n; 
  x = new double [n]; 
  y = new double [n];
  file1 = new char[1000];
  memset(file1, 0, sizeof(*file1)*1000);
  dealloc = true;
}


/**
   The constructor allocates table with n items.

   @param[in] n - the number of table items
   @param[in] it - type of interpolation
   
   Created by JK, 11. 10. 2013
*/  
tablefunct::tablefunct (long n,long it)
{
  if (interpoltype_kwdset.check_int(it)){
    print_err ("unknown type of interpolation is required", __FILE__, __LINE__, __func__);
    abort();
  }

  asize = n; 
  x = new double [n]; 
  y = new double [n];
  file1 = new char[1000];
  memset(file1, 0, sizeof(*file1)*1000);
  dealloc = true;
}


/**
   The constructor makes  table with n items.

   @param[in] xptr - pointer ot the array of x values
   @param[in] yptr - pointer ot the array of y values
   @param[in] n - the number of table items
   @param[in] it - type of interpolation, default value is piecewiselin
   
   Created by TKo, 02. 2024
*/  
tablefunct::tablefunct (double *xptr, double *yptr, long n, interpoltype it)
{
  if (interpoltype_kwdset.check_int(it)){
    print_err ("unknown type of interpolation is required", __FILE__, __LINE__, __func__);
    abort();
  }

  asize = n; 
  x = xptr;
  y = yptr;
  file1 = new char[1000];
  memset(file1, 0, sizeof(*file1)*1000);
  dealloc = false;
}


/**
   The destructor deallocates the memory where the vectors x, y are stored
   Parameters:

   created  20.11.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
*/
tablefunct::~tablefunct()
{
  if (dealloc){
    delete [] x;
    delete [] y;
  }
  delete [] file1;
}



/**
   The function allocates and reads (stores) input parameters (asize, x, y) for interpolation from the file 
   Example: 1 4 0.0 100.0 200.0 300.0 2.4   1.6   8.9  10.0 = itype asize x[0] y[0] x[1] y[1] x[2] y[2] x[3] y[3]

   Parameters: 
   @param in - input file
   
   created  20.11.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
*/
void tablefunct::read(XFILE *in)
{
  //  itype - type of interpolation
  //  itype=1 - piecewise linear interpolation
  //  itype=2 - piecewise constant interpolation

  //  asize - the number of components of arrays x and y
  xfscanf (in,"%k%m%k%ld","approx_type", &interpoltype_kwdset, (int*)&itype, "ntab_items", &asize);
  
  x=new double[asize];
  y=new double[asize];
  
  readval (in);
}

/**
   The function allocates and reads (stores) input parameters (asize, x, y) for interpolation from the data file 
   Example: 1 4 0.0 100.0 200.0 300.0 2.4   1.6   8.9  10.0 = itype asize x[0] y[0] x[1] y[1] x[2] y[2] x[3] y[3]

   Parameters: 
   @param in is input file
   
   16/11/2013, TKr
*/
void tablefunct::read_data_file(XFILE *in)
{
  XFILE *in_data;

  xfscanf(in, " %a", file1);

  in_data = xfopen(file1,"r");

  //  itype - type of interpolation
  //  itype=1 - piecewise linear interpolation

  //  asize - the number of components of arrays x and y
  xfscanf (in_data,"%k%m%k%ld","approx_type", &interpoltype_kwdset, (int*)&itype, "ntab_items", &asize);
  
  x=new double[asize];
  y=new double[asize];
  
  readval (in_data);

  xfclose(in_data);
}


/**
   The function allocates and reads (stores) input parameters (asize, x, y) for interpolation from the data file 
   Example: 1 4 0.0 100.0 200.0 300.0 2.4   1.6   8.9  10.0 = itype asize x[0] y[0] x[1] y[1] x[2] y[2] x[3] y[3]

   Parameters: 
   @param in is input file
   
   26/11/2018, TKr
*/
void tablefunct::read_data_file2(XFILE *in)
{
  XFILE *in_data;

  //  itype - type of interpolation
  //  itype=1 - piecewise linear interpolation

  //  asize - the number of components of arrays x and y
  xfscanf (in,"%k%m%k%ld","approx_type", &interpoltype_kwdset, (int*)&itype, "ntab_items", &asize);

  xfscanf(in, " %a", file1);

  in_data = xfopen(file1,"r");
  
  x=new double[asize];
  y=new double[asize];
  
  readval (in_data);

  xfclose(in_data);
}


/**
   function read values of the table
   
   @param in - input file
   
   created  22.9.2011 by TKo, JK
*/
void tablefunct::readval (XFILE *in)
{
  long i;
  for (i=0;i<asize;i++){
    xfscanf (in,"%lf",&x[i]);
    xfscanf (in,"%lf",&y[i]);
  }

  //  x components are checked
  //  they have to satisfy condition x[0] < x[1] < x[2] < ... <x[asize-1]
  //datacheck ();
}



/**
   The function writes input parameters (asize, x, y) for interpolation

   Parameters: 
   @param out - output file
   
   created  22.11.2002 by JK
*/
void tablefunct::print(FILE *out)
{
  long i;
  fprintf (out,"\n %d %ld",itype,asize);
  for (i=0;i<asize;i++){
    fprintf (out,"\n %le ",x[i]);
    fprintf (out,"%le ",y[i]);
  }
}

/**
   The function writes input parameters for interpolation from file

   Parameters: 
   @param out - output file
   
   created  24/06/2013 by TKr
*/
void tablefunct::print_data_file(FILE *out)
{
  fprintf (out,"\n %s",file1);
}


/**
   The function writes input parameters for interpolation from file

   Parameters: 
   @param out - output file
   
   created  26/11/2018 by TKr
*/
void tablefunct::print_data_file2(FILE *out)
{
  fprintf (out,"\n %d %ld",itype,asize);
  fprintf (out,"\n %s",file1);
}


/**
   The function allocates and reads (stores) input parameters (asize, x, y) for interpolation from the file
   Example: 1 4 0.0 100.0 200.0 300.0 2.4   1.6   8.9  10.0 = itype asize x[0] y[0] x[1] y[1] x[2] y[2] x[3] y[3].
   It was used in the old version of preprocessor and it is obsolate now.
   Parameters:
   @param in is input file

   created  20.11.2002 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
void tablefunct::read_prop(FILE *in)
{
  long i;
  getint(in, (int&)itype);
  getlong(in, asize);
  x=new double[asize];
  y=new double[asize];
  for (i=0;i<asize;i++)
  {
    getdouble (in, x[i]);
    getdouble (in, y[i]);
  }
}



/**
   function checks whether the x components are in strictly increasing order
   there cannot be two identical x components and the following sequence must be
   satisfied
   x[0] < x[1] < x[2] < ... <x[asize-1]
   
   JK, 23.9.2011
*/
void tablefunct::datacheck ()
{
  long i;
  
  for (i=0;i<asize-1;i++){
    if (x[i]<x[i+1]){
      continue;
    }
    else{
      print_err("the %ld-th component (%lf, %lf) does not satisfy condition x[0] < x[1] < x[2] < ... <x[asize-1]", __FILE__, __LINE__, __func__,i,x[i],x[i+1]);
    }
  }
}



/**
   function computes the interpolated value with the help of piecewise constant interpolation
   for x which satisfies x[i] < x <= x[i+1], the function values has the form f(x)=f(x[i+1])
   
   Parameters:
   @param temp is x-coordinate, where the function value will be found
   
   created  20.11.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
   JK, 23.9.2011

   TKo+JK 4.7.2013 changed evaluation of piecewise const function 
                   from (x_i; x_{i+1}> => y[i+1]
                   to   <x_i; x_{i+1}) => y[i]
*/
double tablefunct::piecewise_const_interpol (double temp)
{
  long i;
  double f=0.0;
  
  if (temp<x[0] || temp>=x[asize-1]){
    print_err("required value %le is out of range <%le;%le>", __FILE__, __LINE__, __func__, temp, x[0], x[asize-1]);
    exit(2);
  }
  for (i=0;i<asize-1;i++){
    if (x[i]<=temp && temp<x[i+1]){
      f=y[i];
      break;
    }
  }
  return f;
}



/**
   The function returns the interpolated value at exact point (temp) by Dirac interpolation,
   if the point lies out of the range, extrapolation from the limit intervals takes effect.
   Parameters:
   @param temp is x-coordinate, where the function value will be found
   
   created  8.8.2005 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
double tablefunct::diracinterpol2(double temp)
{
  long i;
  double iv=0.0;
  for (i=0;i<asize;i++){
    if (i==asize-1){
      iv=y[i];
      break;
    }
    if ((i==0) && (x[i] >= temp))
    {
      iv=y[i];
      break;
    }
    if (temp>x[i] && temp<=x[i+1]){
      iv=y[i+1];
      break;
    }
  }
  return iv;
}



/**
   The function returns the interpolated value at exact point (temp) by Dirac interpolation,
   if the point lies out of the range, extrapolation from the limit intervals takes effect.
   Parameters:
   @param temp is x-coordinate, where the function value will be found
   @param k - unused parametr is set to 0.0
   
   created  8.8.2005 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
double tablefunct::diracinterpol2(double temp, double &k)
{
  long i;
  double iv=0.0;
  for (i=0;i<asize;i++){
    if (i==asize-1){
      iv=y[i];
      break;
    }
    if ((i==0) && (x[i] >= temp))
    {
      iv=y[i];
      break;
    }
    if (temp>x[i] && temp<=x[i+1]){
      iv=y[i+1];
      break;
    }
  }
  k = 0.0;
  return iv;
}



/**
   function computes the interpolated value with the help of piecewise constant interpolation
   for x which satisfies x[i] <= x <= x[i+1], the function values has the form
   f(x)=(f(x[i+1])-f(x[i]))/(x[i+1]-x[i])*(x-x[i]) + f(x[i])
   
   Parameters:
   @param temp is x-coordinate, where the function value will be found
   
   created  20.11.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
   JK, 23.9.2011
*/
double tablefunct::piecewise_linear_interpol (double temp)
{
  long i;
  double f=0.0,k;
  
  if (temp<x[0] || temp>x[asize-1]){
    print_err("required value %le is out of range <%le;%le>", __FILE__, __LINE__, __func__, temp, x[0], x[asize-1]);
    exit(2);
  }
  
  for (i=0;i<asize-1;i++){
    if (x[i]<=temp && temp<=x[i+1]){
      k =(y[i+1]-y[i])/(x[i+1]-x[i]);
      f=y[i]+(temp-x[i])*k;
      break;
    }
  }
  
  return f;
}


/**
   function computes the interpolated value with the help of piecewise constant interpolation
   for x which satisfies x[i] <= x <= x[i+1], the function values has the form
   f(x)=(f(x[i+1])-f(x[i]))/(x[i+1]-x[i])*(x-x[i]) + f(x[i])
   
   Parameters:
   @param temp is x-coordinate, where the function value will be found
   
   created  20.11.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
   JK, 23.9.2011
*/
double tablefunct::inverse_piecewise_linear_interpol (double temp)
{
  long i;
  double f=0.0,k;
  
  if (temp<y[0] || temp>y[asize-1]){
    print_err("required value %le is out of range <%le;%le>", __FILE__, __LINE__, __func__, temp, y[0], y[asize-1]);
    exit(2);
  }
  
  for (i=0;i<asize-1;i++){
    if (y[i]<=temp && temp<=y[i+1]){
      k =(x[i+1]-x[i])/(y[i+1]-y[i]);
      f=x[i]+(temp-y[i])*k;
      break;
    }
  }
  
  return f;
}

/**
   The function returns the interpolated value at exact point (temp) by linear interpolation,
   if the point lies out of the range, extrapolation from the limit intervals takes effect.
   Parameters:
   @param temp is x-coordinate, where the function value will be found
   @param k is derivative of linear function
   
   created  8.8.2005 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
double tablefunct::lininterpol2 (double temp, double &k)
{
  //this function is priceless for instance for assembling load value at current time etc.
  
  double q,iv;

  int i,h=0;
  
  for (i=0;i<asize;i++){
    if ((x[i] >= temp) || (i == asize-1)){
      if (i == 0)
      {
        k =(y[i+1]-y[i])/(x[i+1]-x[i]);
        q =y[i]-k*x[i];
        iv=k*temp+q;
        h=1;
        break;
      }     
      k =(y[i]-y[i-1])/(x[i]-x[i-1]);
      q =y[i]-k*x[i];
      iv=k*temp+q;
      h=1;
      break;
    }
  }
  if(h==0){
    print_err("out of range value is required", __FILE__, __LINE__, __func__);
    exit(2);
  }
  
  return iv;
}



/**
   The function returns the interpolated value at exact point (temp) by linear interpolation of inverse function!!!,
   if the point lies out of the range, extrapolation from the limit intervals takes effect.
   Parameters:
   @param temp is x-coordinate, where the function value will be found
   @param k is derivative of linear function
   
   modified 11.2.2007 by Tomas Krejci, krejci@cml.fsv.cvut.cz
*/
double tablefunct::lininterpol3 (double temp, double &k)
{
  //this function is priceless for instance for assembling load value at current time etc.
  
  double q,iv;

  int i,h=0;
  
  for (i=0;i<asize;i++){
    if ((y[i] >= temp) || (i == asize-1)){
      if (i == 0)
      {
        k =(x[i+1]-x[i])/(y[i+1]-y[i]);
        q =x[i]-k*y[i];
        iv=k*temp+q;
        h=1;
        break;
      }     
      k =(x[i]-x[i-1])/(y[i]-y[i-1]);
      q =x[i]-k*y[i];
      iv=k*temp+q;
      h=1;
      break;
    }
  }
  if(h==0){
    print_err("out of range value is required", __FILE__, __LINE__, __func__);
    exit(2);
  }
  
  return iv;
}



/**
   The function returns the interpolated value at exact point (temp) by Lagrangian n-1 degree polynom
   Parameters:
   @param temp is x-coordinate, where the function value will be found
   
   created  20.11.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
*/
double tablefunct::laginterpol(double temp)
{
  double y1=0.,nom,denom;
  long i,j;
 
  for(i=0; i<asize; i++){
    nom=y[i];
    denom=1.;
    for(j=0; j<asize; j++){
      if(i==j) continue;
      nom=nom*(temp-x[j]);
      denom=denom*(x[i]-x[j]);
      if(denom<1.0e-8) continue; // skip subsequent multipile data for identical xx
    }
    y1=y1+nom/denom;
  }
 
  return y1;
}
      

/**
   function computes first derivative dy/dx

   Parameters:
   @param temp is x-coordinate, where the function value will be found
   
   JK, 23.9.2011
*/
double tablefunct::derivative (double temp)
{
  long i;
  double d=0.0;

  if (temp<x[0] || temp>x[asize-1]){
    print_err("required value %le is out of range <%le, %le>", __FILE__, __LINE__, __func__, temp, x[0], x[asize-1]);
    exit(2);
  }

  for (i=0;i<asize-1;i++){
    if (x[i]<=temp && temp<=x[i+1]){
      d=(y[i+1]-y[i])/(x[i+1]-x[i]);
      break;
    }
  }
  
  return d;
}

/**
   function computes first derivative dx/dy

   Parameters:
   @param temp is y-coordinate, where the function value will be found
   
   JK, 23.9.2011
*/
double tablefunct::inv_derivative (double temp)
{
  long i;
  double d=0.0;

  if (temp<y[0] || temp>y[asize-1]){
    print_err("required value is out of range", __FILE__, __LINE__, __func__);
    exit(2);
  }

  for (i=0;i<asize-1;i++){
    if (y[i]<=temp && temp<=y[i+1]){
      d=(x[i+1]-x[i])/(y[i+1]-y[i]);
      break;
    }
  }
  
  return d;
}



/**
   The function returns the interpolated value at exact point (temp), interpolation is chosen according to 
   variable itype. Dirac, linear and Lagrangian n-1 interpolations are supported, see tablefunct.h  
   Parameters:
   @param temp is x-coordinate, where the function value will be found
   
   created  20.11.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
*/
double tablefunct::getval(double temp)
{
    
  double iv;
  if (asize<=1 || x==NULL || y==NULL){
    print_err("wrong size of arrays or wrong number of components", __FILE__, __LINE__, __func__);
    exit(2);
  }
  
  switch(itype){
  case piecewiselin:{
    iv=piecewise_linear_interpol (temp);
    break;
  }
  case piecewiseconst:{
    iv=piecewise_const_interpol (temp);
    break;
  }
  case dirac2:{
    double dummy;
    iv=diracinterpol2(temp, dummy);
    break;
  }
  case lagrange:{
    iv=laginterpol(temp);
    break;
  }
  case piecewiselin2:{//added by TKr 2.11.2007
    iv=piecewise_linear_interpol (temp);
    break;
  }
  default:{
    print_err("unknown type of interpolation", __FILE__, __LINE__, __func__);
    exit(2);
  }
  }   
  
  return iv;
}


/**
   The function returns the interpolated value at exact point (temp), interpolation is chosen according to 
   variable itype. Dirac, linear and Lagrangian n-1 interpolations are supported, see tablefunct.h  
   Parameters:
   @param temp is x-coordinate, where the function value will be found
   
   created  2. 10. 2013 by JK
*/
double tablefunct::getinvval(double temp)
{
    
  double iv;
  if (asize<=1 || x==NULL || y==NULL){
    print_err("wrong size of arrays or wrong number of components", __FILE__, __LINE__, __func__);
    exit(2);
  }
  
  switch(itype){
  case piecewiselin:{
    iv=inverse_piecewise_linear_interpol (temp);
    break;
  }
  default:{
    print_err("unknown type of interpolation", __FILE__, __LINE__, __func__);
    exit(2);
  }
  }   
  
  return iv;
}


/**
   The function returns the interpolated value at exact point (temp), interpolation is chosen according to 
   variable itype. Dirac, linear and Lagrangian n-1 interpolations are supported, see tablefunct.h  
   Parameters:
   @param temp is x-coordinate, where the function value will be found
   @param k is derivative of the function at given interval
   
   created  05.2007 by Tomas Koudelka koudelka@cml.fsv.cvut.cz
*/
double tablefunct::getval2(double temp, double &k)
{
    
  double iv;
  if (asize<=1 || x==NULL || y==NULL){
    print_err("wrong size of the array %ld is required or table values are not initialized",__FILE__,__LINE__,__func__, asize);
    exit(2);
  }
  
  switch(itype){
  case piecewiselin:
  case piecewiselin2:{
    iv=lininterpol2(temp, k);
    break;
  }
  case dirac2:{
    iv=diracinterpol2(temp, k);
    break;
  }
  default:{
    print_err("unknown type of interpolation %d is required",__FILE__,__LINE__,__func__, int(itype));
    exit(2);
  }
  }   
  
  return iv;
}



/**
   The function returns the interpolated value at exact point (temp) of inverse function!!!, interpolation is chosen according to 
   variable itype. Dirac, linear and Lagrangian n-1 interpolations are supported, see tablefunct.h  
   Parameters:
   @param temp is x-coordinate, where the function value will be found
   @param k is derivative of the function at given interval
   
   modified  11.2007 by Tomas Krejci krejci@cml.fsv.cvut.cz
*/
double tablefunct::getval3(double temp, double &k)
{
    
  double iv;
  if (asize<=1 || x==NULL || y==NULL){
    print_err("unknown definition of material parameter is required",__FILE__,__LINE__,__func__);
    exit(2);
  }
  
  switch(itype){
  case piecewiselin:
  case piecewiselin2:{
    iv=lininterpol3(temp, k);
    break;
  }
  default:{
    print_err("unknown definition of material parameter is required",__FILE__,__LINE__,__func__);
    exit(2);
  }
  }   
  
  return iv;
}



/**
   function computes first derivative
   
   @param temp is x-coordinate, where the function value will be found
   
   JK, 23.9.2011
*/
double tablefunct::getderiv (double temp)
{
  double d;

  if (asize<=1 || x==NULL || y==NULL){
    print_err("wrong size of arrays or wrong number of components", __FILE__, __LINE__, __func__);
    exit(2);
  }
  
  switch(itype){
  case piecewiselin:{
    d = derivative (temp);
    break;
  }
  default:{
    print_err("unknown type of interpolation", __FILE__, __LINE__, __func__);
    exit(2);
  }
  }   
  
  return d;
}



/**
   The function copies all data from subroutine argument to the actual object.
   
   @param tf - object of tablefunct for the copying

   @return The function copies the actual object in the parameter tf to the given instance.  

   JK, 15.11.2008
*/
void tablefunct::copy (tablefunct &tf)
{
  long i;
  
  //  type of interpolation
  itype = tf.itype;
  //  number of function values
  asize = tf.asize;
  
  if (x!=NULL)
    delete [] x;
  x = new double [asize];
  
  if (y!=NULL)
    delete [] y;
  y = new double [asize];

  for (i=0;i<asize;i++){
    x[i]=tf.x[i];
    y[i]=tf.y[i];
  }
  size_t l = strlen(tf.file1);
  if (l > 0){
    if (file1)   delete [] file1;
    file1 = new char[1000];
    memset(file1, 0, sizeof(*file1)*1000);
    memcpy(file1, tf.file1, l+1);
  }
  else{
    if (file1)   delete [] file1;
    file1 = new char[1000];
    memset(file1, 0, sizeof(*file1)*1000);
  }
}



/**
   The function compares all data from subroutine argument with the actual object.

   @param tf - object of tablefunct for comparison

   @retval 0 - objects are identical
   @retval 1 - objects differs  

   TKo, 09.2010
*/
long tablefunct::compare(tablefunct &tf)
{
  long i;
  
  //  check for type of interpolation
  if (itype != tf.itype)
    return 1;
    
  //  check for number of function values
  if (asize != tf.asize)
    return 1; 
  
  // check for table values
  for (i=0;i<asize;i++){
    if (x[i] != tf.x[i])
      return 1;
    if (y[i] != tf.y[i])
      return 1;
  }
  return 0;
}
