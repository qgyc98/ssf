#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "climatcond.h"
#include "globalt.h"
#include "aliast.h"
#include "iotools.h"

climatcond::climatcond (void)
{
  //  Reference temperature for enthalpy (0 C)
  T0 = 273.15;
  //  Specific heat capacity of pure liquid water
  Cwat = 4180.0;
  //  Specific enthalpy of water vapour at 0 C
  hvap = 2256000.0;
  //  Specific heat capacity of water vapour
  cvap = 1617.0;    
  
  //  last time when values were computed
  last_time = -1.0;
  //  computer zero
  zero = 1.0e-10;
  
  
  //  type of material model (see alias.h)
  mattyp = 0;
  //  wall orientation in radians (rotation around vertical axis)
  basicdirect = 0.0;
  //  wall inclination in radians (rotation around horizontal axis)
  basicinclin = 0.0;
  
  
  KIND1 = KIND2 = KIND3 = KIND4 = KIND5 = KIND6 = WCorPH = RHorVP = 0;

  cfmax = WindFact = WindFactVD = coefAlf = AbsCoef = Albedo = Latitude = EmisCoeff = ExchCoeff =  ExchCoeffVD = RainCoeff = MinRainTemp = MinRainFlux = u = a = n = alfaK = Ax = Bx = Cx = 0.0;

  SI =  numberSI = d = 0;
  SIdata = NULL;
  
  aa=0.0;  b=0.0;
  eps=0.0;
  m=0.0;
  alpha=0.0;
  sig=5.67e-8;
  tn = 0;


  //  temperature
  TepVal=0.0;
  //  relative humidity
  RHVal=0.0;
  //  diffuse radiation flux on a horizontal plane
  DifSVal=0.0;
  //  direct radiation flux on a horizontal plane
  DirSVal=0.0;
  //  wind direction
  WDVal=0.0;
  //  wind velocity
  WSVal=0.0;
  //  vertical rain flow density, normal to the ground
  VRVal=0.0;
  //  long wave radiation emitted by a horizontal surface
  LWVal=0.0;
  //long wave radiation reflected on a horizontal surface by the clouds and particles in the atmosfere
  SLWVal=0.0;
  // rain flow density normal to the surface area of the wall
  NorRainVal=0.0;
  // imposed gaseous mass flux (dry air+ water vapor)
  GasFluxVal=0.0;
  // imposed heat flux
  HeatFluxVal=0.0;
  // imposed short wave radiation flux normal to the wall surface
  ShWRadVal=0.0;
  // imposed long wave radiation flux normal to the wall surface
  LoWRadVal=0.0;
  //  array containing radiation (JK+Tomaschko)
  RadVal=0.0;
  //  array containing temperature Tb (JK+Klier)
  TempbVal=0.0;

  files0 = new char[1000];
  files1 = new char[1000];
  files2 = new char[1000];
  files3 = new char[1000];
  files4 = new char[1000];
  files5 = new char[1000];
  files6 = new char[1000];
  files7 = new char[1000];
  files8 = new char[1000];
  files14 = new char[1000];
  files15 = new char[1000];

}

climatcond::~climatcond (void)
{
  if(SIdata != NULL)
    delete [] SIdata ;
  
  
  delete [] files0;
  delete [] files1;
  delete [] files2;
  delete [] files3;
  delete [] files4;
  delete [] files5;
  delete [] files6;
  delete [] files7;
  delete [] files8;
  delete [] files14;
  delete [] files15;
  
}

//---------------------------------------------------------------------------

/**
   The function reads data about a climatic condition from the opened text file.

   @param in - pointer to the opened text file
   @param readclimfile - flag for reading of climatic data file (used in the preprocessor).

   Default value of the parameter is 1.

   4. 10. 2013
*/
void climatcond::read (XFILE *in, long readclimfile)
{
  long klier=0;

  //  type of material model (see alias.h)
  xfscanf (in,"%d",(long*)&mattyp);
  //  wall orientation in radians (rotation around vertical axis)
  xfscanf (in,"%lf",&basicdirect);
  //  wall inclination in radians (rotation around horizontal axis)
  xfscanf (in,"%lf",&basicinclin);

  

  //  type of heat conduction
  xfscanf (in,"%ld",&KIND1);

  switch (KIND1){
  case 0:{
	//  no heat conduction
	break;
  }
  case 12:{
	//  heat exchange coefficent
	xfscanf (in,"%lf",&coefAlf);
	break;
  }
  case 13:{
	//  extended heat conduction
	//  wind factor
	xfscanf (in,"%lf %lf",&coefAlf, &WindFact);
	break;
  }
  case 14:{
	//  JK+Tomaschko
	//  aa - coefficient at velocity
	//  b - absolute coefficient
	//  eps - coefficient of emissivity
	//  alpha - absorption coefficient
	//  m - parameter of cloudiness
	xfscanf (in,"%lf %lf   %lf %lf %lf",&aa,&b,&eps,&alpha,&m);
	break;
  }
  case 15:{
	//  JK+Klier
	//  b - transmission coefficient for wind influence
	//  aa -
	//  eps - radiation coefficient
	xfscanf (in,"%lf %lf %lf",&b,&aa,&eps);

	klier=1;
	break;
  }
  case 16:{
	//  tint - temperature
	//  coefAlf - heat exchange coefficent
	xfscanf (in,"%lf %lf",&tint,&coefAlf);
	break;
  }
  default:{
	  print_err("unknown definition of Read clim.data-KIND1 is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (KIND1){



  //  type of short wave radiation
  xfscanf (in,"%ld",&KIND2);

  switch (KIND2){
  case 0:{
	//  no short wave radiation
	break;
  }
  case 21:{
	//  absorption coefficient for short wave radiation
	xfscanf (in,"%lf",&AbsCoef);
	break;
  }
  case 22:{
	//  absorption coefficient for short wave radiation
	//  ground reflextion coefficient
	//  geographical latitude of location in radians
	xfscanf (in,"%lf %lf %lf",&AbsCoef,&Albedo,&Latitude);
	break;
  }
  default:{
	print_err("unknown definition of Read clim.data-KIND2 is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (KIND2){



  //  type of long wave radiation heat flux
  xfscanf (in,"%ld",&KIND3);

  switch (KIND3){
  case 0:{
	//  no long wave radiation heat flux
	break;
  }
  case 31:{
	//  emission coefficient for long wave radiation (-) (it is between 0 and 1)
	xfscanf (in,"%lf",&EmisCoeff);
	break;
  }
  case 32:{
	//  emission coefficient for long wave radiation (-) (it is between 0 and 1)
	xfscanf (in,"%lf",&EmisCoeff);
	break;
  }
  default:{
	print_err("unknown definition of Read clim.data-KIND3 is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (KIND3)


  //   type of water contact
  xfscanf (in,"%ld",&KIND4);

  switch (KIND4){
  case 0:{
	//  no water contact
	break;
  }
  default:{
	print_err("unknown definition of Read clim.data-KIND4 is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (KIND4){




  //  type of mass of convective liquid water flux
  xfscanf (in,"%ld",&KIND5);

  switch (KIND5){
  case 0:{
	//  no rain
	break;
  }
  case 51:{
	//  rain exposure coefficient kg/m^2/s
	//  rain coefficient (it is between 0 and 1)
	//  minimum rain temperature from -40 C to 0 C (K)
	//  minimum normal rain intensity (l/m^2/s)
	xfscanf (in,"%lf %lf %lf %lf",&ExchCoeff,&RainCoeff,&MinRainTemp,&MinRainFlux);
	break;
  }
  case 52:{
	//  rain exposure coefficient kg/m^2/s
	//  rain coefficient (it is between 0 and 1)
	//  minimum rain temperature from -40 C to 0 C (K)
	//  minimum normal rain intensity (l/m^2/s)
	xfscanf (in,"%lf %lf %lf %lf",&ExchCoeff,&RainCoeff,&MinRainTemp,&MinRainFlux);
	break;
  }
  default:{
	print_err("unknown definition of Read clim.data-KIND5 is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (KIND5){




  ///  type of mass of diffusive water vapour flux
  xfscanf (in,"%ld",&KIND6);

  switch (KIND6){
  case 0:{
	//  no mass of diffusive water vapour flux
	break;
  }
  case 61:{
	// water vapor exchange coefficient
	// type - relative humidity or environmental water vapour pressure
	// type = 0 ... relative humidity
	// type = 1 ... nvironmental water vapour pressure
	xfscanf (in,"%lf %ld",&ExchCoeffVD,&RHorVP);
	break;
  }
  case 62:{
	// water vapor exchange coefficient
	// type - relative humidity or environmental water vapour pressure
	// type = 0 ... relative humidity
	// type = 1 ... nvironmental water vapour pressure
	//  wind factor (-) (greater or equal to zero)
	xfscanf (in,"%lf %ld %lf",&ExchCoeffVD,&RHorVP,&WindFactVD);
	break;
  }
  case 63:{
	// water vapor exchange coefficient
	// type - relative humidity or environmental water vapour pressure
	// type = 0 ... relative humidity
	// type = 1 ... nvironmental water vapour pressure

	xfscanf (in,"%lf %ld %lf",&ExchCoeffVD,&RHorVP,&alfaK);
	break;
  }
  case 64:{
	// alfa K - coefficient water vapor exchange coefficient
	// type - relative humidity or environmental water vapour pressure
	// type = 0 ... relative humidity
	// type = 1 ... nvironmental water vapour pressure

	xfscanf (in,"%lf %ld %lf",&ExchCoeffVD,&RHorVP,&alfaK);
	break;
  }
  case 65:{
	// Ax, Bx, Cx - coefficient of ...... water vapor exchange coefficient
	// type - relative humidity or environmental water vapour pressure
	// type = 0 ... relative humidity
	// type = 1 ... nvironmental water vapour pressure
	xfscanf (in,"%lf %ld %lf %lf",&Ax,&RHorVP,&Bx,&Cx);
	break;
  }
  case 66:{
	//  d - number of  .......water vapor exchange coefficient
	// type - relative humidity or environmental water vapour pressure
	// type = 0 ... relative humidity
	// type = 1 ... nvironmental water vapour pressure
	xfscanf (in,"%ld %d",&RHorVP,&d);

	exchangecoeff.tfunc=tab;
	exchangecoeff.tabf = new tablefunct();
	exchangecoeff.tabf->asize = d;
	exchangecoeff.read (in);
	
	break;
  }

  case 67:{
	// constant relative humidity of ambient air
	// ExchCoeffVD - water vapor exchange coefficient
	// constant temperature of ambiente air
	xfscanf (in,"%lf %lf %lf",&rhint,&ExchCoeffVD,&rhtemp);
	break;
  }
  case 68:{
	// constant relative humidity of ambient air
	// ExchCoeffVD - water vapor exchange coefficient
	// constant temperature of ambiente air
	xfscanf (in,"%lf %lf %lf",&rhint,&ExchCoeffVD,&rhtemp);
	break;
  }
  case 70:{
	xfscanf (in,"%lf %lf %lf",&rhint, &rhalfa, &rhtemp);
    break;
  }
  default:{
    print_err ("unknown definition of Read clim.data-KIND6 is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (KIND6){


  //  file with temperature
  xfscanf(in, " %a", files0);
  if (readclimfile)
    read_climatic_files(files0,0);

  //  file with the relative humidity
  xfscanf(in, " %a", files1);
  if (readclimfile)
	read_climatic_files(files1,1);

  //  file with diffuse sun radiation
  xfscanf(in, " %a", files2);
  if (readclimfile)
	read_climatic_files(files2,2);

  //  file with direct sun radiation
  xfscanf(in, " %a", files3);
  if (readclimfile)
	read_climatic_files(files3,3);

  //  file with wind direction
  xfscanf(in, " %a", files4);
  if (readclimfile)
    read_climatic_files(files4,4);
  
  //  file with wind velocity
  xfscanf(in, " %a", files5);
  if (readclimfile)
    read_climatic_files(files5,5);
  
  //  file with rain flux density
  xfscanf(in, " %a", files6);
  if (readclimfile)
    read_climatic_files(files6,6);
  
  //  file with long wave radiation
  xfscanf(in, " %a", files7);
  if (readclimfile)
    read_climatic_files(files7,7);
  
  //  file with heat emission of surrounding ground
  xfscanf(in, " %a", files8);
  if (readclimfile)
    read_climatic_files(files8,8);
  
  //  file with radiation (JK + Tomaschko)
  xfscanf(in, " %a", files14);
  if (readclimfile)
	read_climatic_files(files14,14);
  
  if (klier==1){
    //  file with radiation (JK + Klier)
    xfscanf(in, " %a", files15);
    if (readclimfile)
      read_climatic_files(files15,15);
  }

}

/**
   function prints data to output file

   @param out - output file

   4. 10. 2013
*/
void climatcond::print(FILE *out)
{
  fprintf (out,"\n %ld ",mattyp);
  fprintf (out," %le ",basicdirect);
  fprintf (out," %le ",basicinclin);


  fprintf (out,"\n %ld ",KIND1);

  switch (KIND1){
  case 0:{
	//  KIND1 =  0 - no heat conduction
	break;
  }
  case 11:{
	//  KIND1 = 11 - prescribed heat flux (on boundary)
	fprintf (out," %e ",coefAlf);
	break;
  }
  case 12:{
	//  KIND1 = 12 - heat transmission
	fprintf (out," %e ",coefAlf);
	break;
  }
  case 13:{
	//  KIND1 = 13 - heat transmission extended by the wind factor
	fprintf (out," %e %e ",coefAlf,WindFact);
	break;
  }
  case 14:{
	//  aa - coefficient at velocity
	//  b - absolute coefficient
	//  eps - coefficient of emissivity
	//  alpha - absorption coefficient
	//  m - parameter of cloudiness
	//  KIND1 = 14 - JK+Tomaschko
	fprintf (out," %e %e  %e %e %e ",aa,b,eps,alpha,m);
  break;
  }
  case 15:{
	//  KIND1 = 15 - JK+Klier
	//  b - transmission coefficient for wind influence
	//  aa -
	//  eps - radiation coefficient
	fprintf (out," %ld %le %le",d,aa,eps);
	break;
  }
  case 16:{
	//  KIND1 = 16 - constant heat transmission
	fprintf (out," %e %e",tint,coefAlf);
	break;
  }
  default:{
	print_err("unknown definition of Read clim.data-KIND1 is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (KIND1){



  fprintf (out,"\n %ld ",KIND2 );

  switch (KIND2){
  case 0:{
	//  KIND2 = 0 - no short wave radiation heat flux
	break;
  }
  case 21:{
	//  KIND2 = 21 - short wave radiation heat flux
	fprintf (out," %e ",AbsCoef  );
	break;
  }
  case 22:{
	//  KIND2 = 22 - short wave radiation heat flux extended by the wind
	fprintf (out," %e %e %e ",AbsCoef , Albedo , Latitude );
	break;
  }
  default:{
	print_err ("unknown definition of clim.data-KIND2 is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (KIND2){


  fprintf (out,"\n %ld ",KIND3);

  switch (KIND3){
  case 0:{
	//  KIND3 = 0 - no long wave radiation heat flux
	break;
  }
  case 31:{
	//  KIND3 = 31 - long wave radiation heat flux
	fprintf (out," %e ",EmisCoeff );
	break;
  }
  case 32:{
	//  KIND3 = 32 - long wave radiation heat flux extended by sky radiation
	fprintf (out," %e ",EmisCoeff );
	break;
  }
  default:{
	print_err ("unknown definition of clim.data-KIND3 is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (KIND3){


  fprintf (out,"\n %ld ",KIND4);

  switch (KIND4){
  case 0:{
	//  KIND4 = 0 - no water contact
	break;
  }
  case 41:{
	//  KIND4 = 41 - water contact
	fprintf (out," %ld ",WCorPH );
	break;
  }
  default:{
	print_err ("unknown definition of clim.data-KIND4 is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (KIND4){


  fprintf (out,"\n %ld ",KIND5);


  switch (KIND5){
  case 0:{
	//  KIND5 =  0 - no rain
	break;
  }
  case 51:{
	//  KIND5 = 51 - prescribed water flux caused by rain
	//  rain exposure coefficient kg/m^2/s
	//  rain coefficient (it is between 0 and 1)
	//  minimum rain temperature from -40 C to 0 C (K)
	//  minimum normal rain intensity (l/m^2/s)
	fprintf (out," %e %e %e %e ",ExchCoeff , RainCoeff , MinRainTemp, MinRainFlux);
	break;
  }
  case 52:{
	//  KIND5 = 52 - prescribed water flux caused by rain extended by wind effect
	//  rain exposure coefficient kg/m^2/s
	//  rain coefficient (it is between 0 and 1)
	//  minimum rain temperature from -40 C to 0 C (K)
	//  minimum normal rain intensity (l/m^2/s)
	fprintf (out," %e %e %e %e ",ExchCoeff , RainCoeff , MinRainTemp, MinRainFlux);
	break;
  }
  default:{
	print_err ("unknown definition of clim.data-KIND5 is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (KIND5){


  fprintf (out,"\n %ld ",KIND6);

   // mass of diffusive water vapor flux
  switch (KIND6){
  case 0:{
	break;
  }
  case 61:{
	// water vapor exchange coefficient
	// type - relative humidity or environmental water vapour pressure
	fprintf (out," %e %ld ",ExchCoeffVD, RHorVP);
	break;
  }
  case 62:{
	// water vapor exchange coefficient
	// type - relative humidity or environmental water vapour pressure
	//  wind factor (-) (greater or equal to zero)
	fprintf (out," %e %ld %e ", ExchCoeffVD, RHorVP, WindFactVD);
	break;
  }
  case 63:{
	// water vapor exchange coefficient
	// type - relative humidity or environmental water vapour pressure
	fprintf (out," %e %ld %e ", ExchCoeffVD, RHorVP, alfaK);
	break;
  }
  case 64:{
	// water vapor exchange coefficient
	// type - relative humidity or environmental water vapour pressure
	fprintf (out," %e %ld %e ", ExchCoeffVD, RHorVP, alfaK);
	break;
  }
  case 65:{
	fprintf (out," %e %ld %e %e ", Ax, RHorVP, Bx, Cx);
	break;
  }
  case 66:{
	fprintf (out," %ld ",RHorVP);
	exchangecoeff.print (out);
	break;
  }
  case 67:{
    // ExchCoeffVD - constant relative humidity of ambient air
    // water vapor exchange coefficient
    // constant temperature of ambiente air
    fprintf (out," %e %e %e",rhint, ExchCoeffVD, rhtemp);
    break;
  }
  default:{
    print_err ("unknown definition of clim.data-KIND6 is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (KIND6){


  fprintf (out,"\n%s",files0);
  fprintf (out,"\n%s",files1);
  fprintf (out,"\n%s",files2);
  fprintf (out,"\n%s",files3);
  fprintf (out,"\n%s",files4);
  fprintf (out,"\n%s",files5);
  fprintf (out,"\n%s",files6);
  fprintf (out,"\n%s",files7);
  fprintf (out,"\n%s",files8);
  fprintf (out,"\n%s",files14);
  if (KIND1==15)
    fprintf (out,"\n%s",files15);

  fprintf (out,"\n");
}

//---------------------------------------------------------------------------

/**
   function reads files with climatic data

   @param - file name
   
   4. 10. 2013
*/
void climatcond::read_climatic_files (char *Cfile, long cod)
{
  long time;
  long i;
  double Val;

  XFILE *ffile;
 /* char *zahlavi;



  if (strcmp(Cfile,"0") == 0)
	return;


  ffile = xfopen(Cfile,"r");
  if (ffile==NULL){
	print_err ("Input file %s for climatic condition is not found. \n",__FILE__,__LINE__,__func__,Cfile);
	abort();
  }
  zahlavi = new char[ffile->give_maxlnsize()];

  do{
	xfscanf (ffile," %a",zahlavi);
  }while((strstr(zahlavi,"Data=")) == NULL);
   */

  //xfscanf (ffile,"%ld %ld",&tn,&cc);
  if (strcmp(Cfile,"0") == 0)
	return;

  ffile = xfopen(Cfile,"r");

  if (ffile==NULL){
	print_err ("Input file %s for climatic condition is not found. \n",__FILE__,__LINE__,__func__,Cfile);
	abort();
  }

  xfscanf (ffile,"%ld",&tn);

  switch (cod){
  case 0:{
	// temperature
	TepValue.tfunc=tab;
	TepValue.tabf = new tablefunct (tn,1);
	break;
  }
  case 1:{
	// relative humidity
	RHValue.tfunc=tab;
	RHValue.tabf = new tablefunct (tn,1);
	break;
  }
  case 2:{
	//  diffuse radiation flux on a horizontal plane
	DifSValue.tfunc=tab;
	DifSValue.tabf = new tablefunct (tn,1);
	break;
  }
  case 3:{
	//  direct radiation flux on a horizontal plane
	DirSValue.tfunc=tab;
	DirSValue.tabf = new tablefunct (tn,1);
	break;
  }
  case 4:{
	//  wind direction
	WDValue.tfunc=tab;
	WDValue.tabf = new tablefunct (tn,1);
	break;
  }
  case 5:{
	//  wind velocity
	WSValue.tfunc=tab;
	WSValue.tabf = new tablefunct (tn,1);
	break;
  }
  case 6:{
    //  vertical rain flow density, normal to the ground
	VRValue.tfunc=tab;
	VRValue.tabf = new tablefunct (tn,1);
	break;
  }
  case 7:{
	//  long wave radiation emitted by a horizontal surface
	LWValue.tfunc=tab;
	LWValue.tabf = new tablefunct (tn,1);
    break;
  }
  case 8:{
    //long wave radiation reflected on a horizontal surface by the clouds and particles in the atmosfere
	SLWValue.tfunc=tab;
	SLWValue.tabf = new tablefunct (tn,1);
    break;
  }
  case 9:{
    // rain flow density normal to the surface area of the wall
	NorRainValue.tfunc=tab;
	NorRainValue.tabf = new tablefunct (tn,1);
    break;
  }
  case 10:{
    // imposed gaseous mass flux (dry air+ water vapor)
	GasFluxValue.tfunc=tab;
	GasFluxValue.tabf = new tablefunct (tn,1);
    break;
  }
  case 11:{
    // imposed heat flux
	HeatFluxValue.tfunc=tab;
	HeatFluxValue.tabf = new tablefunct (tn,1);
    break;
  }
  case 12:{
    // imposed short wave radiation flux normal to the wall surface
	ShWRadValue.tfunc=tab;
	ShWRadValue.tabf = new tablefunct (tn,1);
    break;
  }
  case 13:{
    // imposed long wave radiation flux normal to the wall surface
	LoWRadValue.tfunc=tab;
	LoWRadValue.tabf = new tablefunct (tn,1);
	break;
  }
  case 14:{
	//  array containing radiation (JK+Tomaschko)
	RadValue.tfunc=tab;
	RadValue.tabf = new tablefunct (tn,1);
	break;
  }
  case 15:{
	//  array containing temperature Tb (JK+Klier)
	TempbValue.tfunc=tab;
	TempbValue.tabf = new tablefunct (tn,1);
	break;
  }
  default:{
	print_err ("unknown definition of ReadCllimaticFiles is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (cod){





  for (i=0;i<tn;i++){
	//i = 0;
	//while ((xfscanf (ffile,"%ld %ld %ld %ld %ld %lf", &years, &days, &hours, &minutes, &sec, &Val)) != EOF){

  //  xfscanf (ffile,"%ld %ld %ld %ld %ld %lf", &years, &days, &hours, &minutes, &sec, &Val);

	//  time in seconds
   // time = (years * 31536000) + (days * 86400) + (hours * 3600) + (minutes * 60) + sec;
	xfscanf (ffile,"%ld %lf", &time, &Val);
	switch (cod){
	case 0:{
	  // temperature
	  TepValue.tabf->x[i]=time;
	  TepValue.tabf->y[i]=Val;
      break;
    }
    case 1:{
      // relative humidity
      RHValue.tabf->x[i]=time;
      RHValue.tabf->y[i]=Val;
      break;
    }
    case 2:{
      //  diffuse radiation flux on a horizontal plane
      DifSValue.tabf->x[i]=time;
      DifSValue.tabf->y[i]=Val;
      break;
    }
    case 3:{
      //  direct radiation flux on a horizontal plane
      DirSValue.tabf->x[i]=time;
      DirSValue.tabf->y[i]=Val;
      break;
    }
    case 4:{
      //  wind direction
      WDValue.tabf->x[i]=time;
      WDValue.tabf->y[i]=Val;
      break;
    }
    case 5:{
      //  wind velocity
      WSValue.tabf->x[i]=time;
      WSValue.tabf->y[i]=Val;
      break;
    }
    case 6:{
      //  vertical rain flow density, normal to the ground
      VRValue.tabf->x[i]=time;
      VRValue.tabf->y[i]=Val;
      break;
    }
    case 7:{
      //  long wave radiation emitted by a horizontal surface
      LWValue.tabf->x[i]=time;
	  LWValue.tabf->y[i]=Val;
      break;
    }
    case 8:{
      //long wave radiation reflected on a horizontal surface by the clouds and particles in the atmosfere
      SLWValue.tabf->x[i]=time;
      SLWValue.tabf->y[i]=Val;
      break;
    }
    case 9:{
      // rain flow density normal to the surface area of the wall
      NorRainValue.tabf->x[i]=time;
      NorRainValue.tabf->y[i]=Val;
      break;
    }
    case 10:{
      // imposed gaseous mass flux (dry air+ water vapor)
      GasFluxValue.tabf->x[i]=time;
      GasFluxValue.tabf->y[i]=Val;
      break;
    }
    case 11:{
      // imposed heat flux
      HeatFluxValue.tabf->x[i]=time;
      HeatFluxValue.tabf->y[i]=Val;
      break;
    }
    case 12:{
      // imposed short wave radiation flux normal to the wall surface
      ShWRadValue.tabf->x[i]=time;
      ShWRadValue.tabf->y[i]=Val;
      break;
    }
    case 13:{
      // imposed long wave radiation flux normal to the wall surface
      LoWRadValue.tabf->x[i]=time;
      LoWRadValue.tabf->y[i]=Val;
      break;
    }
    case 14:{
      //  array containing radiation (JK+Tomaschko)
      RadValue.tabf->x[i]=time;
      RadValue.tabf->y[i]=Val;
      break;
    }
    case 15:{
      //  array containing temperature Tb (JK+Klier)
	  TempbValue.tabf->x[i]=time;
      TempbValue.tabf->y[i]=Val;
      break;
    }
    default:{
      print_err ("unknown definition of ReadCllimaticFiles is required",__FILE__,__LINE__,__func__);
    }
    }//  end of switch (cod){
    
    //i++;
  }//  end of the statement while ((fscanf (ffile,"%ld %ld %ld %ld %ld %lf", &years, &days, &hours, &minutes, &sec, &Tep)) != EOF){
  
  
  xfclose (ffile);
 // delete [] zahlavi;
}

//---------------------------------------------------------------------------

/**
   function computes required fluxes

   @param flid - flux id
   @param nval - nodal values of density fluxes
   @param ipp - integration point id
   @param nid - node id

   JM 2011, revised JK 4. 10. 2013
*/
double climatcond::give_flux (long flid, double *nval,long ipp,long nid)
{
  double flux,time;
  mednamest mn;
  
  //  the actual time
  time = Tp->time;
  
  //  all actual values are checked and they are eventually computed
  prepare_values (time);

  //  name/names of transported media
  mn = Tp->mednam;
  
  switch (mn){
  case 1:{
    //  heat transport
    flux = give_heat_flux (0,nval[0],ipp, nid);
	break;
  }
  case 2:{
    //  moisture transport
    flux = give_moisture_flux (nval[0], -1000,ipp, nid);
    break;
  }
  case 10:{
    //  heat and moisture transport
    switch (flid){
    case 0:{
      //  moisture flux
      flux = give_moisture_flux (nval[0], nval[1],ipp, nid);
      break;
    }
    case 1:{
      //  heat flux
      flux = give_heat_flux (nval[0],nval[1],ipp, nid);
	  break;
    }
    default:{
      print_err ("unknown flux is required",__FILE__,__LINE__,__func__);
    }
    }//  end of switch (flid){
    break;
  }
  case 25:{
    //  moisture and salt transport with salt crystallization
    switch (flid){
    case 0:{
      //  moisture transport
      flux = give_moisture_flux (nval[0], 293.15,ipp, nid);
      break;
    }
    case 1:{
      //  salt transport
      flux = give_salt_flux (nval[0],nval[1], 293.15, nid);
      break;
    }
    default:{
      print_err ("unknown flux is required",__FILE__,__LINE__,__func__);
    }
    }//  end of switch (flid){
    break;
  }
  case 30:{
    //  heat, moisture and salt transport with crystallization
    switch (flid){
	case 0:{
      //  moisture flux
	  flux = give_moisture_flux (nval[0],nval[3],ipp,nid);
	  break;
	}
    case 1:{
      //  salt flux
      flux = give_salt_flux (nval[0],nval[1],nval[3],nid);
      break;
    }
    case 3:{
      //  heat flux
      flux = give_heat_flux (nval[0],nval[3],ipp,nid);
	  break;
    }
    default:{
      print_err ("unknown flux is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }
  default:{
    print_err ("unknown medium/media flux is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (mn)
  
  return flux;
}




// *****************************************************************************
// *****************************************************************************
//
//  HEAT FLUXES
//
// *****************************************************************************
// *****************************************************************************




/**
   function computes the heat flux
   
   @param nval_0 - moisture
   @param temp - temperature of structure
   
   4. 10. 2013
*/
double climatcond::give_heat_flux (double nval_0,double temp,long ipp, long nid)
{
  double JH1, JH2, JH3, JH7, JH8;
  double time;
  mednamest mn;

  //  meduim/media name
  mn = Tp->mednam;

  //  actual time
  time = Tp->time;

  switch (KIND1){
  case 0:{
	JH1 = 0;
	break;
  }
  case 11:{
	//KIND1 = 11 - prescribed heat flux (on boundary)
	JH1 = condit_heat_flux (11,temp);
	break;
  }
  case 12:{
	//  KIND1 = 12 - heat transmission
	JH1 = condit_heat_flux (12,temp);
	break;
  }
  case 13:{
	//  KIND1 = 13 - heat transmission extended by the wind factor
	JH1 = condit_heat_flux (13,temp);
	break;
  }
  case 14:{
	//  KIND1 = 14 - JK+Tomaschko
	JH1 = condit_heat_flux (14,temp);
	break;
  }
  case 15:{
	//  KIND1 = 15 - JK+Klier
	JH1 = condit_heat_flux (15,temp);
	break;
  }
  case 16:{
	//  KIND1 = 16 - constant heat transmission
	JH1 = condit_heat_flux (16,temp);
	break;
  }
  default:{
	print_err ("unknown definition of Get_Heat_Flux-KIND1 is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (KIND1)



  //  type of short wave radiation heat flux
  switch (KIND2){
  case 0:{
   //  KIND2 = 0 - no short wave radiation heat flux
   JH2 = 0;
	break;
  }
  case 21:{
	//  KIND2 = 21 - short wave radiation heat flux
	JH2 =-1* condit_short_wave_radiat_flux(21,time);
	break;
  }
  case 22:{
	//  KIND2 = 22 - short wave radiation heat flux extended by the wind
	JH2 =-1* condit_short_wave_radiat_flux(22,time);
	break;
  }
  default:{
	print_err ("unknown definition of Get_Heat_Flux-KIND2 is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (KIND2)


  //  type of long wave radiation heat flux
  switch(KIND3){
  case 0:{
	//  KIND3 = 0 - no long wave radiation heat flux
	JH3 = 0;
	break;
  }
  case 31:{
	//  KIND3 = 31 - long wave radiation heat flux
	JH3 =-1* condit_long_wave_radiat_flux (31);
	break;
  }
  case 32:{
	//  KIND3 = 32 - long wave radiation heat flux extended by sky radiation
	JH3 =-1* condit_long_wave_radiat_flux (32);
	break;
  }
  default:{
	print_err ("unknown definition of Get_Heat_Flux-KIND3 is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (KIND3)


  switch(mn){
  case 30:
  case 10:{
	switch (KIND5){
	case 0:{
	  JH7 = 0;
	  break;
	}
	case 51:{
	  //internal energy of convective liquid water flux (J/m2/s)
	  JH7 = condit_rain_heat_flux(51,nval_0,nid);
	  break;
	}
	case 52:{
	  // internal energy of convective liquid water flux (J/m2/s)
	  JH7 = condit_rain_heat_flux(52,nval_0,nid);
	  break;
	}
	default:{
	  print_err ("unknown definition of Get_Heat_Flux-KIND5 is required",__FILE__,__LINE__,__func__);
	}
	}//end of switch (KIND5)
	break;
  }
  default:{
	JH7 = 0;
	break;
  }
  }//  end of switch (mn)-meduim/media name


  switch(mn){
  case 30:
  case 10:{
	switch (KIND6){
	case 0:{
	  JH8 =0;
	  break;
	}
	case 61:{
	   // enthalpy of diffusive water vapor flux (J/m2/s)
	   JH8 = condit_vapour_diffusion_heat_flux(61,nval_0,temp,ipp, nid);
	   break;
	}
	case 62:{
	  // enthalpy of diffusive water vapor flux (J/m2/s) extended by the wind factor
	// The extended model uses an exchange coefficient that depends on the wall direction
	// and on the wind direction and wind speed
	  JH8 = condit_vapour_diffusion_heat_flux(62,nval_0,temp,ipp, nid);
	  break;
	}
	case 63:{
	  // enthalpy of diffusive water vapor flux (J/m2/s)
	  JH8 = condit_vapour_diffusion_heat_flux(63,nval_0,temp,ipp, nid);
	  break;
	}
	case 64:{
	  // enthalpy of diffusive water vapor flux (J/m2/s)
	  JH8 = condit_vapour_diffusion_heat_flux(64,nval_0,temp,ipp, nid);
	  break;
	}
	case 65:{
	  // enthalpy of diffusive water vapor flux (J/m2/s)
	  JH8 = condit_vapour_diffusion_heat_flux(65,nval_0,temp,ipp, nid);
	  break;
	}
	case 66:{
	   //enthalpy of diffusive water vapor flux (J/m2/s)
	   JH8 = condit_vapour_diffusion_heat_flux(66,nval_0,temp,ipp, nid) ;
	   break;
	}
	case 67:{
	  // enthalpy of diffusive water vapor flux (J/m2/s)
	  JH8 = condit_vapour_diffusion_heat_flux(67,nval_0,temp,ipp, nid) ;
	  break;
	}
	default:{
	  print_err ("unknown definition of Get_Heat_Flux-KIND6 is required",__FILE__,__LINE__,__func__);
	}
	}//end of switch (KIND6)
	break;
  }
  default:{
	JH8 = 0;
	break;
  }
  }//  end of switch (mn)-meduim/media name


  return(JH1 +JH2 +JH3 + JH7 + JH8);
}



/**
   function computes heat flux

   @param tkind - type of flux computation
   @param temperstr - temperature of the structure (domain)
   @param time - actual time

   JM, revised 4. 10. 2013
*/
double climatcond::condit_heat_flux (long tkind, double temperstr)
{
  double betaWind,coefAlf0;


  switch (tkind){
  case 11:{
	//  flux defined from climatic file
	//  KIND1 = 11 - prescribed heat flux (on boundary)
	Jkq = HeatFluxVal;
	break;
  }
  case 12:{
	//  KIND1 = 12 - heat transmission
	// coefAlf - heat exchange coefficient
	Jkq = -1.0 * coefAlf * (TepVal - temperstr);
	break;
  }
  case 13:{
	//  KIND1 = 13 - heat transmission extended by the wind factor
	// The extended model uses an exchange coefficient that depends on the wall direction
	// and on the wind direction and wind speed
	// WDVal - wind direction
	// WSVal - wind velocity
	// WindFact - wind factor (-) (greater or equal to zero)
	if ((fabs(basicdirect - WDVal)) <= M_PI){
	  betaWind = fabs (basicdirect - WDVal);
	}else{
	  betaWind = 2.0*M_PI - (fabs (basicdirect - WDVal));
	}

	if (betaWind < (M_PI/2.0)){
	  coefAlf0 = coefAlf + WindFact*sqrt(WSVal);
	}else{
	  coefAlf0 = coefAlf;
	}

	Jkq = -1.0 *(coefAlf0 * (TepVal - temperstr));

	break;
  }
  case 14:{
	//  KIND1 = 14 - JK+Tomaschko
	//  aa - coefficient at velocity
	//  b - absolute coefficient
	//  eps - coefficient of emissivity
	//  alpha - absorption coefficient
	//  m - parameter of cloudiness
	coefAlf0 = b + aa * WSVal;
	Jkq = -1.0*(coefAlf0 *(TepVal - temperstr));
	Jkq += -1.0*eps*sig*(TepVal*TepVal*TepVal*TepVal - temperstr*temperstr*temperstr*temperstr);
	Jkq += -1.0*alpha*m*RadVal;

	break;
  }
  case 15:{
	//  KIND1 = 15 - JK+Klier
	//  aa,b - coefficients for computation of transmission coefficent
	//  eps,sig - coefficients for Stefan-Boltzmann law
	Jkq = b *(temperstr - TepVal);
	Jkq += aa * (temperstr - TempbVal);
	Jkq += eps*sig*(temperstr*temperstr*temperstr*temperstr - TepVal*TepVal*TepVal*TepVal);
	Jkq += RadVal;

	break;
  }
  case 16:{
	//  KIND1 = 16 - constant heat transmission
	//  coefAlf - heat exchange coefficient
	//  tint - constant temperature of ambient air
	Jkq = -1.0 * coefAlf *((tint+273.15) - temperstr);

	break;
  }

  default:{
	print_err("unknown definition of ConditHeatFlux is required",__FILE__,__LINE__,__func__);
  }
  }//  end f switch (kind)


  return Jkq;
}


/**
   function computes short wave radiation flux

   @param swkind - type of flux computation
   @param time - actual time

   JM, revised 5. 10. 2013
*/
double climatcond::condit_short_wave_radiat_flux (long swkind, double time)
{
  double MAXSUNDECL = 0.4089306437;
  double HOURFREQU = 0.26179939;
  double YEARFREQU =0.017214206;
  double SPRINGBEGIN =80;
  double th, dt,day, sdt,cdt, hold, hnew, sold, snew, alf, r;
  double JkQdif, JkQdir,Jq;

  switch (swkind){
  case 21:{
	//  KIND2 = 21 - short wave radiation heat flux
	// imposed short wave radiation flux normal to the wall surface
	//  absorption coefficient for short wave radiation (-) (it is between 0 and 1)
	Jq=AbsCoef*ShWRadVal;
	break;
  }
  case 22:{
	//  KIND2 = 22 - short wave radiation heat flux extended by the wind factor
	//  The extended model uses an exchange coefficient that depends on the wall direction
	//  and on the wind direction and wind speed
	//  DifSVal - diffuse radiation flux on a horizontal plane
	//  DirSVal - direct radiation flux on a horizontal plane
	//  basicdirect - wall orientation in radians (rotation around vertical axis)
	//  basicinclin - wall inclination in radians (rotation around horizontal axis)
	//  Albedo - ground reflextion coefficient

	th = time/3600;
	day = time/86400;
	dt = MAXSUNDECL * sin((YEARFREQU * (day-SPRINGBEGIN )));

	sdt = sin(Latitude)*sin(dt);
	cdt = cos(Latitude)*cos(dt);

	hold = asin(sdt - cdt*cos(HOURFREQU * (th - 2)));
	hnew = asin(sdt - cdt*cos(HOURFREQU * (th - 1)));

	sold = cos(dt)*(sin(HOURFREQU * (th-2)))/(cos(hold));
	snew = cos(dt)*(sin(HOURFREQU * (th-1)))/(cos(hnew));

	if((th>12) && (snew > sold)){
		alf = 2*M_PI + asin(snew);
	}
	if ((th <= 12) || (snew <= sold)){
		alf = M_PI - asin(snew);
	}

	r = cos(basicinclin) + sin(basicinclin)*(cos(alf -basicdirect))/(tan(hnew));

	if (r<0.0){
		JkQdir = 0.0;
	}else{
		if (r>0.5){
			JkQdir = 0.5 *DirSVal;
		}else{
			JkQdir = r * DirSVal;
		}
	}

	JkQdif  = DifSVal*cos(basicinclin/2.0)*cos(basicinclin/2.0) + Albedo *(DifSVal +DirSVal)*sin(basicinclin/2.0)*sin(basicinclin/2.0);

	Jq=AbsCoef*(JkQdif +JkQdir);

	break;
  }
  default:{
	print_err("unknown definition of ConditHeatFlux is required",__FILE__,__LINE__,__func__);
  }
  }//  end f switch (kind)

  return (Jq);

}

/**
   function computes long wave radiation flux
   @param kwkind - type of flux computation

   JM, revised 5. 10. 2013
*/
double climatcond::condit_long_wave_radiat_flux (long lwkind)
{
  double flux;

  switch (lwkind){
  case 31:{
	//  KIND3 = 31 - long wave radiation heat flux
	//  LoWRadVal - imposed long wave radiation flux normal to the wall surface
	//  EmisCoeff - emission coefficient for long wave radiation (-) (it is between 0 and 1)
	flux = EmisCoeff*LoWRadVal;
	break;
  }
  case 32:{
	//  KIND3 = 32 - long wave radiation heat flux extended by sky radiation
	//  LWVal - long wave radiation emitted by a horizontal surface
	//  SLWVal - long wave radiation reflected on a horizontal surface by the clouds and particles in the atmosfere
	//  EmisCoeff - emission coefficient for long wave radiation (-) (it is between 0 and 1)

	flux = (EmisCoeff * (LWVal + SLWVal));
	break;
  }
  default:{
	print_err("unknown definition of ConditLongWaveRadiatFlux is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (lwkind){

  return (flux);
}

/**
   function computes internal energy of convective liquid water flux (J/m2/s)

   5. 10. 2013
*/
double climatcond::condit_rain_heat_flux (long rainkind, double nval_0,long nid)
{
  double RainMoistureFlux, TemperRain;
  double uwT;  // Specific internal energy of the liquid water uw(T)

  //  rainkind = KIND5 = 51 - prescribed water flux caused by rain
  //  rainkind = KIND5 = 52 - prescribed water flux caused by rain extended by wind effect
  // RainMoistureFlux - mass of convective liquid water flux

  RainMoistureFlux =-1.0 * condit_rain_moisture_flux (rainkind, nval_0,nid);

  TemperRain = pow (RHVal,0.1247) * (109.8 + TepVal) - 109.8;
  // Cwat - Specific heat capacity of pure liquid water
  uwT = Cwat *(TemperRain - T0);

  return uwT*RainMoistureFlux;
}

/**
   function computes enthalpy of diffusive water vapor flux (J/m2/s)

   5. 10. 2013
*/
double climatcond::condit_vapour_diffusion_heat_flux (long vdkind,double nval_0, double nval_1,long ipp, long nid)
{
  double JkmvDiff, Tvap, hv;

  //  type of mass of diffusive water vapour flux
  //  vdkind = KIND6 = 0 - no mass of diffusive water vapour flux
  //  vdkind = KIND6 = 61 - mass of diffusive water vapor flux
  //  vdkind = KIND6 = 62 - mass of diffusive water vapor flux extended by the wind factor
  //  The extended model uses an exchange coefficient that depends on the wall direction
  //  and on the wind direction and wind speed
  //  vdkind = KIND6 = 63 -mass of diffusive water vapour flux - water vapor exchange coefficient is given by function
  //  vdkind = KIND6 = 64 - mass of diffusive water vapour flux - water vapor exchange coefficient is given by function
  //  vdkind = KIND6 = 65 - mass of diffusive water vapour flux - water vapor exchange coefficient is given by function
  //  vdkind = KIND6 = 66 - mass of diffusive water vapour flux - water vapor exchange coefficient is given by a table
  //  vdkind = KIND6 = 67 - mass of diffusive water vapour flux - constant temperature, relative humidity, water vapor exchange coefficient
  // JkmvDiff - mass of diffusive water vapor flux

  JkmvDiff = condit_vapour_diffusion_moisture_flux (vdkind,nval_0,nval_1,ipp, nid);

  Tvap = (TepVal + nval_1 )/2;
  // hv - specific enthalpy of water vapor
  // cvap - Specific heat capacity of water vapour

  hv = cvap * (Tvap -T0) + hvap ;

  return (-1.0 * hv * JkmvDiff);
}


// *****************************************************************************
// *****************************************************************************
//
//  MOISTURE FLUXES
//
// *****************************************************************************
// *****************************************************************************


/**
   function computes the moisture flux density

   @param time - actual time

   4. 10. 2013
*/
double climatcond::give_moisture_flux (double nval_0, double nval_1,long ipp, long nid)
{
  /*
   nval_0 - rh for Kunzel type[1]
	  - pv for Kunzel type[2-4]
	  - w for Grunewaldmat, salt1-4, millymat,

   nval_1 - temperature
  */
  double jh4, jh5, jh6;


  switch (KIND4){
  case 0:{
	//  KIND4 = 0 - no water contact
	jh4 = 0;
	break;
  }
  case 41:{
	jh4 = condit_water_contact ();
	break;
  }
  default:{
	print_err("unknown definition of Get_moisture_flux-KIND4 is required",__FILE__,__LINE__,__func__);
  }
  }// end of switch(KIND4)

  switch (KIND5){
  case 0:{
	//  KIND5 =  0 - no rain
	jh5 = 0;
	break;
  }
  case 51:{
	//  KIND5 = 51 - prescribed water flux caused by rain
	jh5 = -1.0 * condit_rain_moisture_flux (51,nval_0,nid);
	break;
  }
  case 52:{
	//  KIND5 = 52 - prescribed water flux caused by rain extended by wind effect
	jh5 = -1.0 * condit_rain_moisture_flux (52,nval_0,nid);
	break;
  }
  default:{
	print_err("unknown definition of Get_moisture_flux-KIND5 is required",__FILE__,__LINE__,__func__);
  }
  }// end of switch(KIND5)

  if (nval_1 != -1000){
	switch(KIND6){
	case 0:{
	  jh6 = 0;
	  break;
	}
	case 61:{
	  jh6 = condit_vapour_diffusion_moisture_flux (61,nval_0,nval_1,ipp, nid);
	  break;
	}
	case 62:{
	// KIND6 = 62 - mass of diffusive water vapor flux extended by the wind factor
	//  The extended model uses an exchange coefficient that depends on the wall direction
	//  and on the wind direction and wind speed
	  jh6 = condit_vapour_diffusion_moisture_flux (62,nval_0,nval_1,ipp, nid);
	  break;
	}
	case 63:{
	  jh6 = condit_vapour_diffusion_moisture_flux (63,nval_0,nval_1,ipp, nid);
	  break;
	}
	case 64:{
	  jh6 = condit_vapour_diffusion_moisture_flux (64,nval_0,nval_1,ipp, nid);
	  break;
	}
	case 65:{
	  jh6 = condit_vapour_diffusion_moisture_flux (65,nval_0,nval_1,ipp, nid);
	  break;
	}
	case 66:{
	  jh6 = condit_vapour_diffusion_moisture_flux (66,nval_0,nval_1,ipp, nid);
	  break;
	}
	case 67:{
	  jh6 = condit_vapour_diffusion_moisture_flux (67,nval_0,nval_1,ipp, nid);
	  break;
	}
	case 68:{
	  double betaExch, pv,VapPres;
	  betaExch = ExchCoeffVD;
	  pv = nval_0;
	  VapPres = rhint;
	  jh6 = (betaExch *(pv - VapPres));  // pv je hodnota v konstrukci, VapPress predepsanan na exterieru
	  fprintf(Outt, " tok je  %le\n",jh6);
	break;
	}
	default:{
	  print_err("unknown definition of Get_moisture_flux-KIND6 is required",__FILE__,__LINE__,__func__);
	}
	} // end of switch(KIND6)
  }

  return(jh4 +jh5 +jh6);
}







/**
   function
   
   5. 10. 2013
*/
double climatcond::condit_water_contact ()
{
  return 0.0;
}


/**
   function computes mass of convective liqiud water flux

   @param rainkind - type of flux computation
   @param nid - node id

   JM, revised 5. 10. 2013
*/
double climatcond::condit_rain_moisture_flux (long rainkind, double nval_0,long nid)
{
  double NorRain, MoistureAkt, Smc;
  double kwind, betawind, krainEff, jkNorRain, jkmaxWat;
  double c1, c2, c3, c4, c5;

  //VRVal = VRVal/3600.0;

  if (WCorPH == 1)
	abort();


  switch (mattyp){
  case 150:{
    MoistureAkt = nval_0; //moisture content u(kg/kg)
    Smc = 10000; // jeste nutno dodelat;
    break;
  }
  case 154:{
    MoistureAkt = Tt->nodes[nid].eqother[1];
    Smc = Tt->nodes[nid].eqother[3];
    break;
  }
  case 155:{
    MoistureAkt = Tt->nodes[nid].eqother[0];
    Smc = Tt->nodes[nid].eqother[2];
    
    //Tm->kun[]->sorption_isotherm_value ();
    break;
  }
  case 158:{
    MoistureAkt = Tt->nodes[nid].eqother[0];
    Smc = Tt->nodes[nid].eqother[3];
    break;
  }
  case 203:{
    MoistureAkt = nval_0;
    Smc = Tt->nodes[nid].eqother[3];
    break;
  }
  case 202:{
    break;
  }
  case 157:{
    MoistureAkt = nval_0;
    break;
  }
  case 159:{
    MoistureAkt = nval_0;
    Smc = Tt->nodes[nid].eqother[3];
    break;
  }
  case 156:{
    MoistureAkt = nval_0;
    Smc = Tt->nodes[nid].eqother[2];
    break;
  }
  default:{
    print_err("unknown definition of Condit_Rain_moistureFlux is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (mattyp)


  // --------------------------------------

  //  toto nema cenu, je to double
  if (Smc == 0.0)
    return 0.0;
  
  // neni uz to davno spocitane
  //WDVal = WDVal*M_PI/180.0;
  
  //  rainkind = KIND5 = 51 - prescribed water flux caused by rain
  //  rainkind = KIND5 = 52 - prescribed water flux caused by rain extended by wind effect
  switch (rainkind){
  case 51:{
    kwind = 1;
    NorRain = VRVal;
    break;
  }
  case 52:{
    // WDVal - wind direction
    // WSVal - wind velocity
    // WindFact - wind factor (-) (greater or equal to zero)
    if ((fabs(basicdirect - WDVal)) <= M_PI){
      betawind = fabs(basicdirect - WDVal);
    }else{
      betawind = 2*M_PI - fabs(basicdirect - WDVal);
    }
    
    if ((betawind <(M_PI/2)) && (WSVal >0) && (VRVal >0)){
      
      c1 = cos (betawind);
      c2 = 1141*sqrt(3600*VRVal/(WSVal *WSVal*WSVal*WSVal));
      c3 = sqrt (1 + c2);
      c5 = sqrt(3600*VRVal);
      c5 = sqrt(c5);
      c4 = exp(-12/(5*c5));
      kwind = c1/c3 * c4;
      
    }else{
      kwind = 0;
	}

	NorRain = VRVal ;

	break;
  }
  default:{
	print_err("unknown definition of Condit_Rain_moistureFlux is required",__FILE__,__LINE__,__func__);
  }
  }


  if (NorRain <= MinRainFlux){
	krainEff = 1;
  }else{
	krainEff = RainCoeff ;
  }

  if (TepVal <= MinRainTemp){
	jkNorRain = 0;
  }else{
	jkNorRain = kwind * krainEff * NorRain;
  }

  jkmaxWat = ExchCoeff * (Smc - MoistureAkt);

  if (jkNorRain <jkmaxWat){
	return jkNorRain;
  }else{
	return jkmaxWat;
  }

}



/**
   function computes mass of diffusive water vapor flux
   
   @param VDKind - type of mass of diffusive water vapour flux
   @param nval_0 -
   @param nval_1 -
   @param ipp - integration point id
   @param nid - node id

   5. 10. 2013
*/
double climatcond::condit_vapour_diffusion_moisture_flux (long VDKind,double nval_0, double nval_1,long ipp, long nid)
{

  long idm;
  double VapPres=0.0, VapPresF=0.0, MoistureAkt, RHAkt=0.0;
  double betaWind, psate, psati, betaExch, pv, smc;
  double rov,M,R;

  //  water molecular weight
  M = 0.01801528;
  // universal gas constant
  R = 8.314472;


  //pro kunzel je w(m3/m3); pro bazped je u(kg/kg)
  smc = MatConstC[1];


  //  switch with respect to material type
  switch (mattyp){
  case 150:{
	//  Bazant-Pedersen material model

	//pro bazped vstupuje u(kg/kg)=RHAkt, t(K)=TemperAkt
	//moisture content u(kg/kg)
	MoistureAkt = nval_0;
	//  material id
	idm = Tm->ip[ipp].idm;
	RHAkt = Tm->bazped[idm].inverse_sorption_isotherm (MoistureAkt);

	//water vapor pressure in the structure
	psati = exp(23.5771 - 4042.9/(nval_1   - 37.58));
	pv = psati * RHAkt;
	break;
  }

  case 154:{
	//  Milly material model

	MoistureAkt = Tt->nodes[nid].eqother[1];
	smc = Tt->nodes[nid].eqother[3];
	RHAkt = Tt->nodes[nid].eqother[2];

	//water vapour pressure in the structure
	psati = exp(23.5771 - 4042.9/(nval_1   - 37.58));
	pv = psati * RHAkt ;
	break;
  }

  case 155:{
	//  Kunzel material model

	MoistureAkt = Tt->nodes[nid].eqother[0];
	smc = Tt->nodes[nid].eqother[2];
	switch (int(Tt->nodes[nid].eqother[4])){
	case 1:{
	  RHAkt = nval_0;

	  //water vapour pressure in the structure
	  psati = exp(23.5771 - 4042.9/(nval_1   - 37.58));
	  pv = psati * RHAkt ;
	  break;
	}

	case 2:
	case 3:
	case 4:{
	  pv = nval_0;
	  break;
	}
	case 5:
	case 6:{
	  //RHAkt = Tt->nodes[nid].eqother[1];
	  pv = nval_0;
	  MoistureAkt = 0.0;
	  break;
	}
	case 7:{
	  //RHAkt = Tt->nodes[nid].eqother[1];

	  rov = nval_0;
	  pv = rov*R*nval_1/M;
	  MoistureAkt = 0.0001;
	  break;
	}
	case 8:{
	  //RHAkt = Tt->nodes[nid].eqother[1];

	  //rov = nval_0;
	  pv = nval_0;
	  MoistureAkt = Tt->nodes[nid].eqother[0];
	  TepVal = 298.15;
	  break;
	}
	default:{

	}
	}
	break;
  }

  case 158:{
	//  Kunzel2 material model

	smc = Tt->nodes[nid].eqother[3];
	MoistureAkt = Tt->nodes[nid].eqother[0];

	//water vapor pressure in the structure
	psati = exp(23.5771 - 4042.9/(nval_1   - 37.58));
	pv = psati * RHAkt ;
	break;
  }

  case 159:{
	//  simplediscmat material model

	MoistureAkt = nval_0;
	RHAkt = Tt->nodes[nid].eqother[1];
	smc = Tt->nodes[nid].eqother[3];

	//water vapor pressure in the structure
	psati = exp(23.5771 - 4042.9/(nval_1   - 37.58));
	pv = psati * RHAkt ;
	break;
  }

  case 203:{
	//  saltmat4 material model

	MoistureAkt = nval_0;
	RHAkt = Tt->nodes[nid].eqother[1];
	smc = Tt->nodes[nid].eqother[3];
    
    //water vapor pressure in the structure
    psati = exp(23.5771 - 4042.9/(nval_1   - 37.58));
    pv = psati * RHAkt ;
    break;
  }

  case 202:
  case 201:
  case 200:{
	//  saltmat1, saltmat2, saltmat3 material models
    
    MoistureAkt = nval_0;
    RHAkt = Tt->nodes[nid].eqother[0];
	smc = Tt->nodes[nid].eqother[2];
    //water vapor pressure in the structure
	psati = exp(23.5771 - 4042.9/(273.15   - 37.58));
    pv = psati * RHAkt ;
    break;
  }

  case 157:
  case 156:{
    //  Grunewald and Devries material models
    
    MoistureAkt = nval_0;
    RHAkt = Tt->nodes[nid].eqother[0];
    smc = Tt->nodes[nid].eqother[2];
	//water vapor pressure in the structure
    psati = exp(23.5771 - 4042.9/(nval_1   - 37.58));
    pv = psati * RHAkt ;
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (mattyp)
  
  
  
  
  if (smc == 0)
    return(0);
  
  
  // ----------------------------------------------------------------------------------
  // nutne doplnit zavislost pro model bazped, kde je neznama moisture content u(kg/kg)
  //  vdkind - type of mass of diffusive water vapour flux
  //  vdkind = KIND6 = 0 - no mass of diffusive water vapour flux
  //  vdkind = KIND6 = 61 - mass of diffusive water vapor flux
  //  vdkind = KIND6 = 62 - mass of diffusive water vapor flux extended by the wind factor
  //  The extended model uses an exchange coefficient that depends on the wall direction
  //  and on the wind direction and wind speed
  //  vdkind = KIND6 = 63 -mass of diffusive water vapour flux - water vapor exchange coefficient is given by function
  //  vdkind = KIND6 = 64 - mass of diffusive water vapour flux - water vapor exchange coefficient is given by function
  //  vdkind = KIND6 = 65 - mass of diffusive water vapour flux - water vapor exchange coefficient is given by function
  //  vdkind = KIND6 = 66 - mass of diffusive water vapour flux - water vapor exchange coefficient is given by a table
  //  vdkind = KIND6 = 67 - mass of diffusive water vapour flux - constant temperature, relative humidity, water vapor exchange coefficient


  switch (VDKind){
  case 61:{
	// ExchCoeffVD - water vapor exchange coefficient
	betaExch = ExchCoeffVD;
	break;
  }
  case 62:{
	// WDVal - wind direction
	// WSVal - wind velocity
	// WindFact - wind factor (-) (greater or equal to zero)
	// ExchCoeffVD - water vapor exchange coefficient
	WDVal = (M_PI/180.0)*WDVal;
	if (fabs (basicdirect - WDVal) <= M_PI){
	  betaWind = fabs (basicdirect-WDVal);
	}
	else{
	  betaWind = 2*M_PI - fabs(basicdirect-WDVal);
	}

	if (betaWind < (M_PI/2)){
	  betaExch = ExchCoeffVD +WindFactVD * sqrt(WSVal);
	}
	else{
	  betaExch = ExchCoeffVD ;
	}
	break;
  }
  case 63:{
	// water vapor exchange coefficient is given by function
	betaExch = ExchCoeffVD*MoistureAkt + alfaK;
	break;
  }
  case 64:{
	// water vapor exchange coefficient is given by function
	betaExch = ExchCoeffVD*exp(alfaK*MoistureAkt);
	break;
  }
  case 65:{
	// water vapor exchange coefficient is given by function
	betaExch = Ax*MoistureAkt*MoistureAkt + Bx*MoistureAkt + Cx;
	break;
  }
  case 66:{
    // water vapor exchange coefficient is given by table
    betaExch = exchangecoeff.getval (MoistureAkt);
    break;
  }
  case 67:{
	// rhint - constant relative humidity of ambient air
	// ExchCoeffVD - water vapor exchange coefficient
	// rhtemp - constant temperature of ambiente air
	betaExch = ExchCoeffVD;
	//  toto je divne
	TepVal = rhtemp + 273.15;
	//  toto je divne
	RHVal=rhint/100;
	break;
  }
  default:{
	print_err("unknown type of VDKind is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (VDKind)



  // musim doplnit i moznost TRY tlaku misto RH
  //external water vapor pressure
  psate = exp(23.5771 - 4042.9/(TepVal - 37.58));


  if (RHorVP == 0){
	VapPres = psate * RHVal;
  }else{
	VapPres = VapPresF;
  }

  //VapPres = 1585.02;
  double tokp;

  if ((VapPres  - pv) <0)
	{
	  // return (betaExch *(VapPres -pv ));

	  tokp =  (betaExch *(pv - VapPres));
	}
  else
	{
	  // return(((Smc -MoistureAkt)/Smc)*betaExch *(VapPres -pv ));
	 tokp = (((smc-MoistureAkt)/smc)*betaExch *(pv -VapPres));
      //tokp =  (betaExch *(pv - VapPres));

	}
 //fprintf(Outt2, "%ld  %lf %lf  %lf %lf %lf\n",nid,  t, smc, RHAkt,MoistureAkt, RelHum);// kapa a decko
  if(mattyp==155)
  {
	if (int(Tt->nodes[nid].eqother[4])==7) {
	double rove;
	rove = VapPres*M/(R*TepVal);

	}
  }
  
  
  //fprintf (Outt,"\n BAF  ipp %5ld   tokp %le  beta %le  Tep %le   RH %le",ipp,tokp,betaExch,TepVal,RHVal);
  return(tokp);
}















/**
   function

   JM 2011
*/
double climatcond::give_salt_flux(double w, double cf, double temp, long nid)
{
  double jh6salt;

  //  type of mass of diffusive water vapour flux
  //  KIND6 = 0 - no mass of diffusive water vapour flux
  //  KIND6 = 61 - mass of diffusive water vapor flux
  //  KIND6 = 62 - mass of diffusive water vapor flux extended by the wind factor
  //  The extended model uses an exchange coefficient that depends on the wall direction
  //  and on the wind direction and wind speed
  //  KIND6 = 63 -mass of diffusive water vapour flux - water vapor exchange coefficient is given by function
  //  KIND6 = 64 - mass of diffusive water vapour flux - water vapor exchange coefficient is given by function
  //  KIND6 = 65 - mass of diffusive water vapour flux - water vapor exchange coefficient is given by function
  //  KIND6 = 66 - mass of diffusive water vapour flux - water vapor exchange coefficient is given by a table
  //  KIND6 = 67 - mass of diffusive water vapour flux - constant temperature, relative humidity, water vapor exchange coefficien


  if (temp != -1000){
	switch(KIND6){
	case 0:{
	  jh6salt = 0;
	  break;
	}
	case 61:{
	  jh6salt = condit_vapour_diffusion_moisture_flux_salt (61,w,cf,temp, nid);
	  break;
	}
	case 62:{
	  jh6salt = condit_vapour_diffusion_moisture_flux_salt (62,w,cf,temp, nid);
	  break;
	}
	case 63:{
	  jh6salt = condit_vapour_diffusion_moisture_flux_salt (63,w,cf,temp, nid);
	  break;
	}
	case 64:{
	  jh6salt = condit_vapour_diffusion_moisture_flux_salt (64,w,cf,temp, nid);
	  break;
	}
	case 65:{
	  jh6salt = condit_vapour_diffusion_moisture_flux_salt (65,w,cf,temp, nid);
	  break;
	}
	case 66:{
	  jh6salt = condit_vapour_diffusion_moisture_flux_salt (66,w,cf,temp, nid);
      break;
    }
    default:{
      print_err("unknown definition of Get_salt_flux-KIND6 is required",__FILE__,__LINE__,__func__);
    }
    }
  }//  end of the if (temp != -1000){

  
  return jh6salt;
}


/**
   function computes values of all external variables

   @param time - the actual time

   JK, 5. 10. 2013
*/
void climatcond::compute_values (double time)
{
  //  temperature
  if (TepValue.tabf != NULL){
    if (TepValue.tabf->asize>0){
      TepVal = TepValue.getval (time);
      TepVal = TepVal + 273.15;
    }
  }
  //  relative humidity
  if (RHValue.tabf != NULL){
    if (RHValue.tabf->asize>0){
      RHVal = RHValue.getval (time);
      RHVal = RHVal/100.0;
    }
  }
  //  diffuse radiation flux on a horizontal plane
  if (DifSValue.tabf != NULL){
    if (DifSValue.tabf->asize>0){
      DifSVal = DifSValue.getval (time);
    }
  }
  //  direct radiation flux on a horizontal plane
  if (DirSValue.tabf != NULL){
    if (DirSValue.tabf->asize>0){
      DirSVal = DirSValue.getval (time);
    }
  }
  //  wind direction
  if (WDValue.tabf != NULL){
    if (WDValue.tabf->asize>0){
      WDVal = WDValue.getval (time);
      WDVal = WDVal * M_PI / 180.0;
    }
  }
  //  wind velocity
  if (WSValue.tabf != NULL){
    if (WSValue.tabf->asize>0){
      WSVal = WSValue.getval (time);
    }
  }
  //  vertical rain flow density, normal to the ground
  if (VRValue.tabf != NULL){
    if (VRValue.tabf->asize>0){
      VRVal = VRValue.getval (time);
      VRVal = VRVal/3600.0;
    }
  }
  //  long wave radiation emitted by a horizontal surface
  if (LWValue.tabf != NULL){
    if (LWValue.tabf->asize>0){
      LWVal = LWValue.getval (time);
    }
  }
  //long wave radiation reflected on a horizontal surface by the clouds and particles in the atmosfere
  if (SLWValue.tabf != NULL){
    if (SLWValue.tabf->asize>0){
      SLWVal = SLWValue.getval (time);
    }
  }
  // rain flow density normal to the surface area of the wall
  if (NorRainValue.tabf != NULL){
    if (NorRainValue.tabf->asize>0){
      NorRainVal = NorRainValue.getval (time);
    }
  }
  // imposed gaseous mass flux (dry air+ water vapor)
  if (GasFluxValue.tabf != NULL){
    if (GasFluxValue.tabf->asize>0){
      GasFluxVal = GasFluxValue.getval (time);
    }
  }
  // imposed heat flux
  if (HeatFluxValue.tabf != NULL){
    if (HeatFluxValue.tabf->asize>0){
      HeatFluxVal = HeatFluxValue.getval (time);
    }
  }
  // imposed short wave radiation flux normal to the wall surface
  if (ShWRadValue.tabf != NULL){
    if (ShWRadValue.tabf->asize>0){
      ShWRadVal = ShWRadValue.getval (time);
    }
  }
  // imposed long wave radiation flux normal to the wall surface
  if (LoWRadValue.tabf != NULL){
    if (LoWRadValue.tabf->asize>0){
      LoWRadVal = LoWRadValue.getval (time);
    }
  }
  //
  if (RadValue.tabf != NULL){
    if (RadValue.tabf->asize>0){
      RadVal = RadValue.getval (time);
    }
  }
  //
  if (TempbValue.tabf != NULL){
    if (TempbValue.tabf->asize>0){
      TempbVal = TempbValue.getval (time);
    }
  }
  
}

/**
   function prepares the actual values

   @param time - the actual time

   5. 10. 2013
*/
void climatcond::prepare_values (double time)
{
  //  time reduction to segment <0; 31536000>
  while (time >= 31536000){
    time = time - 31536000;
  }
  
  if (fabs(time-last_time) < zero){
    //  all values were computed for the time required
    
  }else{
    compute_values (time);
  }
  
}



/**
   function
   
   5. 10. 2013
*/
double climatcond::condit_vapour_diffusion_moisture_flux_salt (long VDKind,double w, double cf, double temp, long nid)
{
  double VapPres=0.0, WindFactVD=0.0, VapPresF=0.0;
  double rh,betaWind, psat, betaExch, pv, smc;
  

  switch (mattyp){
  case 203:{
    rh = Tt->nodes[nid].eqother[1];
    smc = Tt->nodes[nid].eqother[3];
    cfmax = Tt->nodes[nid].eqother[4];
    break;
  }
  case 202:
  case 201:
  case 200:{
    //  toto je divne, temp sem leze zvenci a tady se prepisuje
    temp = 283.15; // temperature specimentu
    switch (MatCharC[0]){
    case 30:{
      //hvezda = (MatFunceC[0][0])/(1-sqrt(1-(MatFunceC[0][1])));
      //grunw1->give_data_si_root_fi(w,temp,0,MatFunceC[0][0],MatConstC[1],MatFunceC[0][1],hvezda,rh );
      break;}
    case 2:{
      //sorption_izoterm_data(w,rh);
      break;
    }
    }
    
    break;
  }
  default:{
    print_err("unknown definition of Condit_Vapour_Diffusion_moisture_flux_salt is required",__FILE__,__LINE__,__func__);
  }
  }
  
  
  switch (VDKind){
  case 0:{
	// ExchCoeffVD - water vapor exchange coefficient
	betaExch = ExchCoeffVD ;
	break;
  }
  case 61:{
    // WDVal - wind direction
	// WSVal - wind velocity
	// WindFact - wind factor (-) (greater or equal to zero)
	// ExchCoeffVD - water vapor exchange coefficient

	if (fabs (basicdirect-WDVal) <= M_PI){
	  betaWind = fabs (basicdirect-WDVal);
	}else{
	  betaWind = 2*M_PI - fabs(basicdirect-WDVal);
	}

	if (betaWind < (M_PI/2)){
	  betaExch = ExchCoeffVD +WindFactVD * sqrt(WSVal);
	}else{
	  betaExch = ExchCoeffVD ;
	}
	break;
  }
  case 62:{
	// water vapor exchange coefficient is given by function
	betaExch = ExchCoeffVD*w + alfaK;
	break;
  }
  case 63:{
	// water vapor exchange coefficient is given by function
	betaExch = ExchCoeffVD*exp(alfaK*w);
	break;
  }
  case 64:{
	// water vapor exchange coefficient is given by function
	betaExch = Ax*w*w + Bx*w + Cx;
	break;
  }
  case 65:{
    // water vapor exchange coefficient is given by table
    betaExch = exchangecoeff.getval (w);
    break;
  }
  case 67:{
    // rhint - constant relative humidity of ambient air
    // ExchCoeffVD - water vapor exchange coefficient
    // rhtemp - constant temperature of ambiente air
    betaExch = ExchCoeffVD;
    //  toto je divne
    TepVal = rhtemp;
    //  toto je divne
    RHVal=rhint/100;
    break;
  }
  default:{
    print_err("unknown definition of Condit_Vapour_Diffusion_moisture_flux_salt-VDKind is required",__FILE__,__LINE__,__func__);
  }
  }
  
  
  // musim doplnit i moznost TRY tlaku misto RH
  
  psat = exp(23.5771 - 4042.9/(TepVal  - 37.58));
  
  if (RHorVP == 0){
    VapPres = psat * RHVal;
  }else{
    VapPres = VapPresF;
  }

  psat = exp(23.5771 - 4042.9/(temp - 37.58));
  pv = psat * rh ;
  
  // tady je to podle popisu od Roberta - modelovani vzniku vykvetu
  // urcite nutna revize, spise preformulovani, preprogramovani...
  // nevim, zda  je zachovana Ct v celem vzorku
  double tok;
  if (cf >= cfmax){
    
    if ((VapPres  - pv) <0){
      if(Tt->nodes[nid].eqother[11] == 1){
	tok = betaExch *(pv - VapPres)*cf/1000;
	Tt->nodes[nid].eqother[11]= Tt->nodes[nid].eqother[11] + 1;
	Tt->nodes[nid].eqother[12]= Tt->nodes[nid].eqother[12] + tok;
	Tt->nodes[nid].eqother[13]= tok;
	return(tok);
      }else{
	Tt->nodes[nid].eqother[11]= Tt->nodes[nid].eqother[11]+1;
	if(Tt->nodes[nid].eqother[11] ==5)
	  Tt->nodes[nid].eqother[11] = 1;
	return (Tt->nodes[nid].eqother[13]);
      }
    }else{
	  if(Tt->nodes[nid].eqother[11] == 1){
	tok = ((smc -w)/smc)*betaExch *(pv -VapPres)*cf/1000;
	Tt->nodes[nid].eqother[11]= Tt->nodes[nid].eqother[11]+1;
	Tt->nodes[nid].eqother[12]= Tt->nodes[nid].eqother[12]+ tok;
	Tt->nodes[nid].eqother[13]= tok;
	return(tok);
      }else{
	Tt->nodes[nid].eqother[11]= Tt->nodes[nid].eqother[11]+1;
	if(Tt->nodes[nid].eqother[11] ==5)
	  Tt->nodes[nid].eqother[11] = 1;
	return(Tt->nodes[nid].eqother[13]);
      }
    }
  }else{
    if( Tt->nodes[nid].eqother[12] > 1.0e-20){
      if ((VapPres  - pv) <0){
	if(Tt->nodes[nid].eqother[11] == 1){
	  tok =betaExch *(pv - VapPres)*cf/1000;
	  if (fabs(tok) > Tt->nodes[nid].eqother[12])
	    tok = Tt->nodes[nid].eqother[12];
	  Tt->nodes[nid].eqother[11]= Tt->nodes[nid].eqother[11]+1;
	  Tt->nodes[nid].eqother[12]= Tt->nodes[nid].eqother[12]-tok;
	  Tt->nodes[nid].eqother[13]= tok;
	  return(tok);
	}else{
	  Tt->nodes[nid].eqother[11]= Tt->nodes[nid].eqother[11]+1;
	  if(Tt->nodes[nid].eqother[11] ==5)
	    Tt->nodes[nid].eqother[11] = 1;
	  return(Tt->nodes[nid].eqother[13]);
	}
      }else{
	if(Tt->nodes[nid].eqother[11] == 1){
	  tok = ((smc -w)/smc)*betaExch *(pv -VapPres)*cf/1000;
	  if (fabs(tok) > Tt->nodes[nid].eqother[12])
	    tok = Tt->nodes[nid].eqother[12];
	  Tt->nodes[nid].eqother[11]= Tt->nodes[nid].eqother[11]+1;
	  Tt->nodes[nid].eqother[12]= Tt->nodes[nid].eqother[12]-tok;
	  Tt->nodes[nid].eqother[13]= tok;
	  return(tok);
	}else{
	  Tt->nodes[nid].eqother[11]= Tt->nodes[nid].eqother[11]+1;
	  if(Tt->nodes[nid].eqother[11] ==5)
	    Tt->nodes[nid].eqother[11] = 1;
	  return(Tt->nodes[nid].eqother[13]);
	}
      }
    }else{
	  Tt->nodes[nid].eqother[11] = 1;
      Tt->nodes[nid].eqother[12] = 0.0;
      Tt->nodes[nid].eqother[13] = 0.0;
      return(0);
    }
  }
}


