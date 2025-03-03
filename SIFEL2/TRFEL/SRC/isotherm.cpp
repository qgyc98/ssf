#include "isotherm.h"
#include "globalt.h"


isotherm::isotherm ()
{
  //  type of isotherm
  //  the type has to be read from the input file
  isothermtype = noisotherm;
  
}

isotherm::~isotherm ()
{
}

/**
   function reads isotherm
   
   @param in - input file
   
   JK, 2. 10. 2013
*/
void isotherm::read(XFILE *in)
{
  //  type of isotherm
  xfscanf(in,"%k%m", "isothermtype", &isotypet_kwdset,(int*)&isothermtype);
  
  switch (isothermtype){
  case ithdata:{
    //  the isotherm is given by a table
    
    //  reading of the table
    data.read (in);
    
    break;
  }
  case grunewaldroot:{
    //  Grunewald-Root
    rooti.read (in);
    break;
  }
  case hansen:{
    //  Hansen
    hanseni.read (in);
    break;
  }
  default:{
    print_err("unknown type of isotherm is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (isothermtype){

}



/**
   function prints isotherm
   
   @param out - output file
   
   JK, 6. 10. 2013
*/
void isotherm::print (FILE *out)
{
  fprintf(out,"%d\n", isothermtype);

  switch (isothermtype){
  case ithdata:{
    //  isotherm is given by data
    
    data.print (out);
    fprintf (out,"\n");
    
    break;
  }
  case grunewaldroot:{
    //  Grunewald-Root
    rooti.print (out);
    break;
  }
  case hansen:{
    //  Hansen
    hanseni.print (out);
    break;
  }
  default:{
    print_err("unknown type of isotherm is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (isothermtype){
  
}

/**
   function evaluates the value of sorption isotherm
   
   @param in[in] - input value
   @param ipp[in] - integration point pointer (optional)
   @param eid[in] - element id (optional)
   
   9. 10. 2013
*/
double isotherm::isotherm_value (double in, long ipp, long eid)
{
  //  isother value
  double iv;
  
  switch (isothermtype){
  case ithdata:{
    if (data.tfunc == tab)
    {
      if ((in < data.tabf->x[0]) || (in > data.tabf->x[data.tabf->asize-1]))
      {
        print_err("required value %le is out of table range <%le;%le> on ip=%ld, eid=%ld\n", 
                  __FILE__, __LINE__, __func__, in, data.tabf->x[0], data.tabf->x[data.tabf->asize-1], ipp, eid);
        abort();
      }
    }
    iv = data.getval (in);
    break;
  }
  case grunewaldroot:{
    //  Grunewald-Root
    //  in - relative humidity
    //  iv - volumetric moisture content
    iv = rooti.sorption_isotherm (in);
    break;
  }
  case hansen:{
    //  Hansen
    //  in - relative humidity
    //  iv - volumetric moisture content
    hanseni.hansen_sorption_isotherm (in);
    break;
  }
  default:{
    print_err("unknown type of isotherm is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (isothermtype){

  return iv;
}

/**
   @param in - input value
   
   9. 10. 2013
*/
double isotherm::inverse_isotherm_value (double in)
{
  //  inverse isotherm value
  double iiv;
  
  switch (isothermtype){
  case ithdata:{
    iiv = data.getinvval (in);
    break;
  }
  case grunewaldroot:{
    //  Grunewald-Root
    //  in - volumetric moisture content
    //  iiv - relative humidity
    iiv = rooti.inverse_sorption_isotherm (in);
    break;
  }
  case hansen:{
    //  Hansen
    //  in - volumetric moisture content
    //  iiv - relative humidity
    iiv = hanseni.hansen_inverse_sorption_isotherm (in);
    break;
  }
  default:{
    print_err("unknown type of isotherm is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (isothermtype){

  return iiv;
}

/**
   @param in - input value
   
   9. 10. 2013
*/
double isotherm::derivative_isotherm_value (double in)
{
  //  derivative of the isotherm
  double di;
  
  switch (isothermtype){
  case ithdata:{
    di = data.getderiv (in);
    break;
  }
  case grunewaldroot:{
    //  Grunewald-Root
    //  in - relative humidity
    //  di - derivative of the sorption isotherm with respect to relative humidity
    di = rooti.derivative_sorption_isotherm (in);
    break;
  }
  case hansen:{
    //  Hansen
    di = hanseni.derivative_relhum (in);
    break;
  }
  default:{
    print_err("unknown type of isotherm is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (isothermtype){

  return di;
}

/**
   @param in - input value
   
   9. 10. 2013
*/
double isotherm::derivative_inverse_isotherm_value (double in)
{
  //  derivative of the inverse isotherm
  double dii;
  
  switch (isothermtype){
  case ithdata:{
    dii = data.getinvderiv (in);
    break;
  }
  case grunewaldroot:{
    //  Grunewald-Root
    //  in - volumetric moisture content
    //  dii - derivative of the sorption isotherm with respect to volumetric moisture content
    dii = rooti.derivative_inverse_sorption_isotherm (in);
    break;
  }
  case hansen:{
    //  Hansen
    break;
  }
  default:{
    print_err("unknown type of isotherm is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (isothermtype){

  return dii;
}
