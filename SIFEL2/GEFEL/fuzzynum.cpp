#include "fuzzynum.h"
#include <math.h>


fuzzynum::fuzzynum (void)
{
  //  number of alpha-cuts
  nalph=0;
  //  required deviation/error
  err=1.0e-8;
  
  //  array of alpha values
  alph=NULL;
  //  array of minimum values
  xmin=NULL;
  //  array of maximum values
  xmax=NULL;
  
  //  number of values on axis
  nval=0;
  //  array of values on axis
  x = NULL;
}



fuzzynum::~fuzzynum (void)
{
  delete [] alph;
  delete [] xmin;
  delete [] xmax;
  delete [] x;
}

/**
   function reads fuzzy number from input file
   
   @param in - input file
   
   JK, 18.8.2005
*/
void fuzzynum::read (XFILE *in)
{
  long i;
  
  //  number of alpha-cuts
  xfscanf (in,"%ld",&nalph);
  
  if (nalph<0){
    fprintf (stderr,"\n\n number of alpha-cuts is negative in function fuzzynum::read (file %s, line %d)\n",
	     __FILE__,__LINE__);
  }
  
  if (alph==NULL){
    alph = new double [nalph];
  }
  if (xmin==NULL){
    xmin = new double [nalph];
  }
  if (xmax==NULL){
    xmax = new double [nalph];
  }
  
  //  minimum and maximum values
  for (i=0;i<nalph;i++){
    xfscanf(in,"%lf %lf %lf",alph+i,xmin+i,xmax+i);
  }
  
}

/**
   function initiates fuzzy number
   
   @param n - number of alpha cuts
   @param al - array of alpha cuts
   
   JK
*/
void fuzzynum::initiate (long n,double *al)
{
  long i;
  
  //  number of alpha-cuts
  nalph=n;
  
  if (nalph<0){
    fprintf (stderr,"\n\n number of alpha-cuts is negative in function fuzzynum::initiate (file %s, line %d)\n",
	     __FILE__,__LINE__);
  }
  
  if (alph==NULL){
    alph = new double [nalph];
  }
  if (xmin==NULL){
    xmin = new double [nalph];
  }
  if (xmax==NULL){
    xmax = new double [nalph];
  }
  
  //  minimum and maximum values
  for (i=0;i<nalph;i++){
    alph[i] =  al[i];
    xmin[i] =  1.0e40;
    xmax[i] = -1.0e40;
  }
  
}




/**
   function returns minimum value of alpha-cut
   alpha from segment <0;1> is given
   
   @param alpha - alpha-cut value
   
   JK, 18.8.2005
*/
double fuzzynum::give_min (double alpha)
{
  long i,j;
  double min=0.0;
  
  j=0;
  for (i=0;i<nalph;i++){
    if (fabs(alpha-alph[i])<err){
      min=xmin[i];
      j=1;
      break;
    }
  }
  
  if (j==0){
    print_err(" required value of alpha-cut is not found for alpha = %f.",__FILE__,__LINE__, __func__, alpha);
  }
  
  return min;
}

/**
   function returns minimum value of alpha-cut
   
   @param acutid - number of alpha-cut
   
   JK, 18.8.2005
*/
double fuzzynum::give_min (long acutid)
{
  if (acutid>nalph){
    fprintf (stderr,"\n\n number of required alpha-cut is greater than number of all alpha-cuts");
    fprintf (stderr,"\n in function give_min (file %s, line %d)\n",__FILE__,__LINE__);
  }
  
  return xmin[acutid];
}

/**
   function returns maximum value of alpha-cut
   alpha from segment <0;1> is given
   
   @param alpha - alpha-cut value
   
   JK, 18.8.2005
*/
double fuzzynum::give_max (double alpha)
{
  long i,j;
  double max=0.0;
  
  j=0;
  for (i=0;i<nalph;i++){
    if (fabs(alpha-alph[i])<err){
      max=xmax[i];
      j=1;
      break;
    }
  }
  
  if (j==0){
    print_err(" required value of alpha-cut is not found for alpha = %f.", __FILE__, __LINE__, __func__, alpha);
  }
  
  return max;
}

/**
   function returns maximum value of alpha-cut
   
   @param acutid - number of alpha-cut
   
   JK, 18.8.2005
*/
double fuzzynum::give_max (long acutid)
{
  if (acutid>nalph){
    fprintf (stderr,"\n\n number of required alpha-cut is greater than number of all alpha-cuts");
    fprintf (stderr,"\n in function give_max (file %s, line %d)\n",__FILE__,__LINE__);
  }
  
  return xmax[acutid];
}




















/**
   function saves minimum value of alpha-cut
   
   @param alpha - alpha-cut value
   @param min - minimum value alpha-cut
   
   JK, 18.8.2005
*/
void fuzzynum::save_min (double alpha,double min)
{
  long i,j;
  
  j=0;
  for (i=0;i<nalph;i++){
    if (fabs(alpha-alph[i])<err){
      xmin[i]=min;
      j=1;
      break;
    }
  }
  
  if (j==0){
    fprintf (stderr,"\n\n required value of alpha-cut is not found in function save_min");
    fprintf (stderr,"\n (file %s, line %d)\n",__FILE__,__LINE__);
  }
  
}

/**
   function saves minimum value of alpha-cut
   
   @param acutid - number of alpha-cut
   @param min - minimum value of aplha-cut
   
   JK, 18.8.2005
*/
void fuzzynum::save_min (long acutid,double min)
{
  if (acutid>nalph){
    fprintf (stderr,"\n\n number of required alpha-cut is greater than number of all alpha-cuts");
    fprintf (stderr,"\n in function save_min (file %s, line %d)\n",__FILE__,__LINE__);
  }
  
  xmin[acutid]=min;
}

/**
   function saves maximum value and alpha value
   
   @param acutid - number of alpha-cut
   @param al - value of alpha
   @param min - minimum value of aplha-cut
   
   JK, 22.8.2005
*/
void fuzzynum::save_alp_min (long acutid,double al,double min)
{
  if (acutid>nalph){
    fprintf (stderr,"\n\n number of required alpha-cut is greater than number of all alpha-cuts");
    fprintf (stderr,"\n in function save_alp_min (file %s, line %d)\n",__FILE__,__LINE__);
  }
  
  alph[acutid]=al;
  xmin[acutid]=min;
}







/**
   function saves maximum value of alpha-cut
   
   @param alpha - alpha-cut value
   @param max - maximum value alpha-cut
   
   JK, 18.8.2005
*/
void fuzzynum::save_max (double alpha,double max)
{
  long i,j;
  
  j=0;
  for (i=0;i<nalph;i++){
    if (fabs(alpha-alph[i])<err){
      xmax[i]=max;
      j=1;
      break;
    }
  }
  
  if (j==0){
    fprintf (stderr,"\n\n required value of alpha-cut is not found in function save_max");
    fprintf (stderr,"\n (file %s, line %d)\n",__FILE__,__LINE__);
  }
  
}

/**
   function saves maximum value of alpha-cut
   
   @param acutid - number of alpha-cut
   @param max - maximum value of aplha-cut
   
   JK, 18.8.2005
*/
void fuzzynum::save_max (long acutid,double max)
{
  if (acutid>nalph){
    fprintf (stderr,"\n\n number of required alpha-cut is greater than number of all alpha-cuts");
    fprintf (stderr,"\n in function save_max (file %s, line %d)\n",__FILE__,__LINE__);
  }
  
  xmax[acutid]=max;
}

/**
   function saves maximum value and alpha value
   
   @param acutid - number of alpha-cut
   @param al - value of alpha
   @param max - maximum value of aplha-cut
   
   JK, 22.8.2005
*/
void fuzzynum::save_alp_max (long acutid,double al,double max)
{
  if (acutid>nalph){
    fprintf (stderr,"\n\n number of required alpha-cut is greater than number of all alpha-cuts");
    fprintf (stderr,"\n in function save_alp_max (file %s, line %d)\n",__FILE__,__LINE__);
  }
  
  alph[acutid]=al;
  xmax[acutid]=max;
}

/**
   function prints detail informations about fuzzy number
   
   @param out - output stream
   
   JK, 17.4.2007
*/
void fuzzynum::print (FILE *out)
{
  long i;
  
  fprintf (out,"\n number of alpha-cuts %ld",nalph);
  for (i=0;i<nalph;i++){
    fprintf (out,"\n alpha %15.12e   min %15.12e   max %15.12e",alph[i],xmin[i],xmax[i]);
  }
}

/**
   function prints only the minimum and maximum values
   
   @param out - output stream
   
   JK, 17.4.2007
*/
void fuzzynum::print_minmax (FILE *out)
{
  fprintf (out,"\n  %15.12e %15.12e",xmin[0],xmax[0]);
}

/**
   function creates one array of alpha-cut values instead of two arrays
   (one array for minimum and one array for maximum values)
   
   JK, 9.10.2005
*/
/*
void fuzzynum::onearray ()
{
  long i,j;

  nval=nalph*2-1;
  if (x==NULL){
    x=new double [nval];
  }
  else{
    delete [] x;
    x=new double [nval];
  }
  
  j=0;
  for (i=0;i<nalph;i++){
    x[j]=xmin[i];
    j++;
  }
  for (i=nalph-2;i>-1;i--){
    x[j]=xmax[i];
    j++;
  }
}
*/
/**
   function returns value from list of alpha-cuts
   
   JK, 9.10.2005
*/
double fuzzynum::give_val (long i)
{
  return x[i];
}


/**
   function saves value to the fuzzy number
   
   it is intended for problems where minimum and maximum values from a set of values is stored
   nalph should be 1

   if val is less than minimum value, the minimum value is overwritten by the value
   if val is greater than maximum value, the maximum value is overwritten by the value
   
   @param val - stored value
   
   JK, 4.12.2007
*/
void fuzzynum::save_value (double val)
{
  if (xmin[0]>val)
    xmin[0]=val;
  if (xmax[0]<val)
    xmax[0]=val;
}

/**
   function initiates fuzzy number for storage of minimum and maximum values
   
   JK, 4.12.2007
*/
void fuzzynum::minmax_init ()
{
  //  number of alpha-cuts must be 1
  nalph=1;
  
  if (xmin!=NULL)
    delete [] xmin;
  xmin = new double [1];
  xmin[0] =  1.0e40;
  
  if (xmax!=NULL)
    delete [] xmax;
  xmax = new double [1];
  xmax[0] = -1.0e40;
}

