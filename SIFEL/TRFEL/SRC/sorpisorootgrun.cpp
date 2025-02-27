#include "sorpisorootgrun.h"

sorpisorootgrun::sorpisorootgrun (void)
{
  //  the hygroscopic moisture content
  whyg=0.0;
  //  the maximum hygroscopic relative humidity
  rhhyg=0.0;
  //  the saturated volumetric moisture content
  wsat=0.0;

}

sorpisorootgrun::~sorpisorootgrun (void)
{
}

/**
   function reads material parameters
   
   @param in - input file
   
   19. 11. 2012, JM
*/
void sorpisorootgrun::read (XFILE *in)
{
  //  the hygroscopic moisture content
  //  the maximum hygroscopic relative humidity
  //  the saturated volumetric moisture content
  xfscanf(in,"%le %le %le",&whyg,&rhhyg,&wsat);
}

/**
   function prints material parameters

   @param out - output file

   19. 11. 2012, JM
*/
void sorpisorootgrun::print (FILE *out)
{
  fprintf(out, "%le %le %le",whyg,rhhyg,wsat);

}

/**
   function computes value of the Root (Grunewald) sorption isotherm

   @param rh - the relative humidity
   
   JM, 19. 11. 2012
*/
double sorpisorootgrun::sorption_isotherm (double rh)
{
  double moist;
  
  if (rh < rhhyg){
    moist = (1.0-sqrt(1.0-rh)) * whyg/(1.0-sqrt(1.0-rhhyg));
  }else{
    if (rh > 1.0){
      moist = wsat;
    }else{
      moist = whyg + (rh - rhhyg)/(1.0-rhhyg) * (wsat - whyg);
    }
  }
  
  return moist;
}

/**
   function computes value of the Root (Grunewald) inverse sorption isotherm

   @param w - volumetric moisture content
   
   JM, 19. 11. 2012
*/
double sorpisorootgrun::inverse_sorption_isotherm (double w)
{
  double phi,whygs;
  
  if (w < whyg){
    //  w hyd star
    whygs = whyg/(1.0-sqrt(1.0-rhhyg));
    phi = 1.0 - (1.0-w/whygs)*(1.0-w/whygs);
  }else{
    if (w > wsat){
      phi = 1.0;
    }else{
      phi = rhhyg + (1.0-rhhyg)*(w-whyg)/(wsat-whyg);
    }
  }
  
  return phi;
}




/**
   function computes derivative of the Root (Grunewald) sorption isotherm
   with respect to the relative humidity

   @param rh - the relative humidity
   
   JM, 19. 11. 2012
*/
double sorpisorootgrun::derivative_sorption_isotherm (double rh)
{
  double dersi,whygs;
  
  if (rh > 1.0)
    rh = 1.0;
  if (rh < 0.0)
    rh = 0.0;
  
  if (rh < rhhyg){
    //  w hyd star
    whygs = whyg/(1.0-sqrt(1.0-rhhyg));
	//dersi = whygs/2.0/(1-sqrt(rh));
	dersi = whygs/2.0/(sqrt (1-rh));
  }else{
    dersi = (wsat - whyg)/(1.0 - rhhyg);
  }
  
  return dersi;
}

/**
   function computes derivative of the inverse Root (Grunewald) sorption isotherm
   with respect to the volumetric moisture content 

   @param w - volumetric moisture content
   
   JM, 19. 11. 2012
*/
double sorpisorootgrun::derivative_inverse_sorption_isotherm (double w)
{
  double dersi,whygs;
  
  if (w > wsat){
    w = wsat;}
  if (w < 0.0){
    w = 0.0;
  }
  
  if (w < whyg){
    //  w hyd star
    whygs = whyg/(1.0-sqrt(1.0-rhhyg));
    dersi = 2.0*(whygs-w)/whygs/whygs;
  }else{
    dersi = (1.0 - rhhyg)/(wsat - whyg);
  }
  
  return dersi;
}
