#include "climatcond2.h"
#include "globalt.h"

climatcond2::climatcond2 (void)
{
  mattyp = 0;
  //  the Ludolf number
  //M_PI = 3.14159265358979323846;
  //  Reference temperature for enthalpy (0 C)
  T0 = 273.15;
  //  water molecular weight
  M = 0.01801528;
  // universal gas constant
  R = 8.314472;
  //  Specific heat capacity of pure liquid water
  Cwat = 4180.0;
  //  Specific heat capacity of water vapour
  cvap = 1617.0;    
  //  Specific enthalpy of water vapour at 0 C
  hvap = 2256000.0;
 
  //  last time when values were computed
  last_time = -1.0;
  //  computer zero
  zero = 1.0e-10;

  KIND1 = KIND2 = KIND3 = KIND4 = KIND5 = KIND6 = 0;

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

  files0 = new char[1000];
  files1 = new char[1000];
  files2 = new char[1000];
  files3 = new char[1000];
  files4 = new char[1000];
  files5 = new char[1000];
  files6 = new char[1000];
  files7 = new char[1000];
  files8 = new char[1000];
}

climatcond2::~climatcond2 (void)
{
  delete [] files0;
  delete [] files1;
  delete [] files2;
  delete [] files3;
  delete [] files4;
  delete [] files5;
  delete [] files6;
  delete [] files7;
  delete [] files8;
}

/**
   The function reads data about a climatic condition from the opened text file.

   @param in - pointer to the opened text file
   @param readclimfile - flag for reading of climatic data file (used in the preprocessor).

   Default value of the parameter is 1.

   4. 10. 2013
*/
void climatcond2::read (XFILE *in, long readclimfile)
{
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
  
  
  
  
  //  type of mass of diffusive water vapour flux
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
    xfscanf (in,"%lf",&ExchCoeffVD);
    break;
  }
  case 62:{
    // water vapor exchange coefficient
    // type - relative humidity or environmental water vapour pressure
    // type = 0 ... relative humidity
    // type = 1 ... nvironmental water vapour pressure
    //  wind factor (-) (greater or equal to zero)
    xfscanf (in,"%lf %lf",&ExchCoeffVD,&WindFactVD);
    break;
  }
  case 63:{
    // water vapor exchange coefficient
    // type - relative humidity or environmental water vapour pressure
    // type = 0 ... relative humidity
    // type = 1 ... nvironmental water vapour pressure
    
    xfscanf (in,"%lf %lf",&ExchCoeffVD,&alfaK);
    break;
  }
  case 64:{
    // alfa K - coefficient water vapor exchange coefficient
    // type - relative humidity or environmental water vapour pressure
    // type = 0 ... relative humidity
    // type = 1 ... nvironmental water vapour pressure
    
    xfscanf (in,"%lf %lf",&ExchCoeffVD,&alfaK);
    break;
  }
  case 65:{
    // Ax, Bx, Cx - coefficient of ...... water vapor exchange coefficient
    // type - relative humidity or environmental water vapour pressure
    // type = 0 ... relative humidity
    // type = 1 ... nvironmental water vapour pressure
    xfscanf (in,"%lf %lf %lf",&Ax,&Bx,&Cx);
    break;
  }
  case 66:{
    //  d - number of  .......water vapor exchange coefficient
    // type - relative humidity or environmental water vapour pressure
    // type = 0 ... relative humidity
    // type = 1 ... nvironmental water vapour pressure
    xfscanf (in,"%d",&d);
    
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
  
}

/**
   function prints data to output file

   @param out - output file

   4. 10. 2013
*/
void climatcond2::print(FILE *out)
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
	fprintf (out," %e ",ExchCoeffVD);
	break;
  }
  case 62:{
	// water vapor exchange coefficient
	// type - relative humidity or environmental water vapour pressure
	//  wind factor (-) (greater or equal to zero)
	fprintf (out," %e %e ", ExchCoeffVD, WindFactVD);
	break;
  }
  case 63:{
	// water vapor exchange coefficient
	// type - relative humidity or environmental water vapour pressure
	fprintf (out," %e %e ", ExchCoeffVD, alfaK);
	break;
  }
  case 64:{
	// water vapor exchange coefficient
	// type - relative humidity or environmental water vapour pressure
	fprintf (out," %e %e ", ExchCoeffVD, alfaK);
	break;
  }
  case 65:{
	fprintf (out," %e %e %e ", Ax, Bx, Cx);
	break;
  }
  case 66:{
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
  case 68:{
    // constant relative humidity of ambient air
    // ExchCoeffVD - water vapor exchange coefficient
    // constant temperature of ambiente air
    fprintf (out," %e %e %e",rhint,ExchCoeffVD,rhtemp);
    break;
  }
  case 70:{
    fprintf (out," %e %e %e",rhint, rhalfa, rhtemp);
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

  fprintf (out,"\n");
}


/**
   function reads files with climatic data

   @param - file name
   
   11. 4. 2014
*/
void climatcond2::read_climatic_files (char *Cfile, long cod)
{
  long time;
  long i;
  double Val;

  XFILE *ffile;
  
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
  default:{
    print_err ("unknown definition of ReadCllimaticFiles is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (cod){
  
  
  
  for (i=0;i<tn;i++){
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
    default:{
      print_err ("unknown definition of ReadCllimaticFiles is required",__FILE__,__LINE__,__func__);
    }
    }//  end of switch (cod){
    
  }//  end of the statement while ((fscanf (ffile,"%ld %ld %ld %ld %ld %lf", &years, &days, &hours, &minutes, &sec, &Tep)) != EOF){
  
  
  xfclose (ffile);
}


/**
   function prepares the actual values

   @param time - the actual time

   5. 10. 2013
*/
void climatcond2::prepare_values (double time)
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
   function computes values of all external variables

   @param time - the actual time

   JK, 5. 10. 2013
*/
void climatcond2::compute_values (double time)
{
  //  temperature
  if (TepValue.tabf != NULL){
    if (TepValue.tabf->asize>0){
      TepVal = TepValue.getval (time);
      TepVal = TepVal + 273.15;
//	fprintf (stdout,"\n tepval   %le",TepVal);
    }
  }
  //  relative humidity
  if (RHValue.tabf != NULL){
    if (RHValue.tabf->asize>0){
      RHVal = RHValue.getval (time);
	  RHVal = RHVal/100.0;
	  if(RHVal <0.0) RHVal=0.0;
    }
//    fprintf (stdout,"\n RHVal   %le",RHVal);
  }
  //  diffuse radiation flux on a horizontal plane
  if (DifSValue.tabf != NULL){
    if (DifSValue.tabf->asize>0){
      DifSVal = DifSValue.getval (time);
    }
//    fprintf (stdout,"\n DifSVal   %le",DifSVal);
  }
  //  direct radiation flux on a horizontal plane
  if (DirSValue.tabf != NULL){
    if (DirSValue.tabf->asize>0){
      DirSVal = DirSValue.getval (time);
    }
//    fprintf (stdout,"\n DirSVal   %le",DirSVal);
  }
  //  wind direction
  if (WDValue.tabf != NULL){
    if (WDValue.tabf->asize>0){
      WDVal = WDValue.getval (time);
	  WDVal = WDVal * M_PI / 180.0;
	  if(WDVal <0.0) WDVal=0.0;
    }
//    fprintf (stdout,"\n WDVal   %le",WDVal);
  }
  //  wind velocity
  if (WSValue.tabf != NULL){
    if (WSValue.tabf->asize>0){
	  WSVal = WSValue.getval (time);
	  if(WSVal <0.0) WSVal=0.0;
    }
//    fprintf (stdout,"\n WSVal   %le",WSVal);
  }
  //  vertical rain flow density, normal to the ground
  if (VRValue.tabf != NULL){
    if (VRValue.tabf->asize>0){
      VRVal = VRValue.getval (time);
	  VRVal = VRVal/3600.0;
	 // if(VRVal <0.0) VRVal=0.0;
	  if(VRVal <1e-8) VRVal = 0.0;
    }
//    fprintf (stdout,"\n VRVal   %le",VRVal);
  }
  //  long wave radiation emitted by a horizontal surface
  if (LWValue.tabf != NULL){
    if (LWValue.tabf->asize>0){
      LWVal = LWValue.getval (time);
    }
//    fprintf (stdout,"\n LWVal   %le",LWVal);
  }
  //long wave radiation reflected on a horizontal surface by the clouds and particles in the atmosfere
  if (SLWValue.tabf != NULL){
    if (SLWValue.tabf->asize>0){
      SLWVal = SLWValue.getval (time);
    }
//    fprintf (stdout,"\n SLWVal   %le",SLWVal);
  }
  // rain flow density normal to the surface area of the wall
  if (NorRainValue.tabf != NULL){
    if (NorRainValue.tabf->asize>0){
	  NorRainVal = NorRainValue.getval (time);
	  if(NorRainVal <0.0) NorRainVal=0.0;
    }
//    fprintf (stdout,"\n NorRainVal   %le",NorRainVal);
  }
  // imposed gaseous mass flux (dry air+ water vapor)
  if (GasFluxValue.tabf != NULL){
    if (GasFluxValue.tabf->asize>0){
      GasFluxVal = GasFluxValue.getval (time);
    }
//    fprintf (stdout,"\n GasFluxVal   %le",GasFluxVal);
  }
  // imposed heat flux
  if (HeatFluxValue.tabf != NULL){
    if (HeatFluxValue.tabf->asize>0){
      HeatFluxVal = HeatFluxValue.getval (time);
    }
//    fprintf (stdout,"\n HeatFluxVal   %le",HeatFluxVal);
  }
  // imposed short wave radiation flux normal to the wall surface
  if (ShWRadValue.tabf != NULL){
    if (ShWRadValue.tabf->asize>0){
      ShWRadVal = ShWRadValue.getval (time);
    }
//    fprintf (stdout,"\n ShWRadVal   %le",ShWRadVal);
  }
  // imposed long wave radiation flux normal to the wall surface
  if (LoWRadValue.tabf != NULL){
    if (LoWRadValue.tabf->asize>0){
      LoWRadVal = LoWRadValue.getval (time);
    }
//    fprintf (stdout,"\n LoWRadVal   %le",LoWRadVal);
  }
}

/**
   function computes required fluxes

   @param flid - flux id
   @param nval - nodal values
   @param ipp - integration point id
   @param nid - node id

   JM 2011, revised JK 4. 10. 2013
*/
double climatcond2::give_flux (long flid, double *nval,long ipp,long nid)
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
  default:{
    print_err ("unknown medium/media flux is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (mn)
  
  return flux;
}

// *****************************************************************************
// *****************************************************************************
//
//  HEAT FLUXES, VALUES AND COEFFICIENTS
//
// *****************************************************************************
// *****************************************************************************

/**
   function computes the heat flux
   
   @param nval_0 - moisture
   @param temp - temperature of structure
   
   4. 10. 2013
*/
double climatcond2::give_heat_flux (double /*nval_0*/,double /*temp*/,long /*ipp*/, long nid)
{
  double JH1, JH2, JH3, JH7;
  double time;
  mednamest mn;

  //  meduim/media name
  mn = Tp->mednam;

  //  actual time
  time = Tp->time;

  switch (KIND1){
  case 0:{
    JH1 = 0.0;
    break;
  }
  case 12:{
    JH1 = 0.0;
    break;
  }
  case 13:{
    JH1 = 0.0;
    break;
  }
  case 16:{
    JH1=0.0;
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
  case 10:{
    switch (KIND5){
    case 0:{
      JH7 = 0;
      break;
    }
    case 51:{
      //internal energy of convective liquid water flux (J/m2/s)
      JH7 = condit_rain_heat_flux(51,nid);
      break;
    }
    case 52:{
      // internal energy of convective liquid water flux (J/m2/s)
      JH7 = condit_rain_heat_flux(52,nid);
      break;
    }
    default:{
      print_err ("unknown definition of Get_Heat_Flux-KIND5 is required",__FILE__,__LINE__,__func__);
    }
    }//end of switch (KIND5)
    break;
  }

  default:{
    JH7 = 0.0;
    break;
  }
  }//  end of switch (mn)-meduim/media name
  

  return(JH1 +JH2 +JH3 + JH7);
}


/**
   function computes short wave radiation flux

   @param swkind - type of flux computation
   @param time - actual time

   JM, revised 5. 10. 2013, 11. 4. 2014
*/
double climatcond2::condit_short_wave_radiat_flux (long swkind, double time)
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
    print_err("unknown definition is required",__FILE__,__LINE__,__func__);
  }
  }//  end f switch (kind)

  return (Jq);

}

/**
   function computes long wave radiation flux
   @param kwkind - type of flux computation

   JM, revised 5. 10. 2013
*/
double climatcond2::condit_long_wave_radiat_flux (long lwkind)
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
double climatcond2::condit_rain_heat_flux (long rainkind,long nid)
{
  double RainMoistureFlux, TemperRain;
  double uwT;  // Specific internal energy of the liquid water uw(T)

  //  rainkind = KIND5 = 51 - prescribed water flux caused by rain
  //  rainkind = KIND5 = 52 - prescribed water flux caused by rain extended by wind effect
  // RainMoistureFlux - mass of convective liquid water flux

  RainMoistureFlux =-1.0 * condit_rain_moisture_flux (rainkind,nid);

  TemperRain = pow (RHVal,0.1247) * (109.8 + TepVal) - 109.8;
  // Cwat - Specific heat capacity of pure liquid water
  uwT = Cwat *(TemperRain - T0);

  return uwT*RainMoistureFlux;
}


/**
   function computes transmission coefficient for heat flux
   
   @param VDKind - type of mass of diffusive water vapour flux
   
   JK 10. 4. 2014
*/
double climatcond2::transmission_coefficient_tt (long tkind)
{
  double excoef,betaWind,coefAlf0;
  
  switch (tkind){
    case 0:{
    excoef = 0.0;
    break;    
    }
    case 12:{
    //  KIND1 = 12 - heat transmission
    // coefAlf - heat exchange coefficient
    excoef = coefAlf;
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
    
    excoef = coefAlf0;
    break;
  }
  case 16:{
    //  KIND1 = 16 - constant heat transmission
    //  coefAlf - heat exchange coefficient
    excoef = coefAlf;
    
    break;
  }
    
  default:{
    print_err("unknown type is required",__FILE__,__LINE__,__func__);
  }
  }//  end f switch (kind)
  
  return excoef;
}

/**
   function computes temperature of exterior
   it is used in the Newton boundary condition
   
   @param VDKind - type of mass of diffusive water vapour flux
   
   JK 11. 4. 2014
*/
double climatcond2::transmission_exterior_value_tt (long tkind)
{
  double tval;
  
  switch (tkind){
  case 0:{
    //  KIND1 = 0 - no condition for heat
    tval = -100000.0;
    break;
  }
  case 12:{
    //  KIND1 = 12 - heat transmission
    tval = TepVal;
    break;
  }
  case 13:{
    //  KIND1 = 13 - heat transmission extended by the wind factor
    // The extended model uses an exchange coefficient that depends on the wall direction
    // and on the wind direction and wind speed
    tval = TepVal;
    break;
  }
  case 16:{
    //  KIND1 = 16 - constant heat transmission
    //  tint - constant temperature of ambient air
    tval = tint+273.15;
    
    break;
  }
    
  default:{
    print_err("unknown type is required",__FILE__,__LINE__,__func__);
  }
  }//  end f switch (kind)
  
  return tval;
}

/**
   function computes transmission coefficient tw
   
   @param VDKind - type of mass of diffusive water vapour flux
   @param temp - temperature of structure
   
   JK 11. 4. 2014
*/
double climatcond2::transmission_coefficient_tw (long VDKind,double temp)
{

  double c,betaWind,betaExch, tval;


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
  case 0:{
    betaExch = 0.0;
    tval=temp;
    break;
  }
  case 61:{
	// ExchCoeffVD - water vapor exchange coefficient
	betaExch = ExchCoeffVD;
	tval = TepVal;
	break;
  }
  case 62:{
	// WDVal - wind direction
	// WSVal - wind velocity
	// WindFact - wind factor (-) (greater or equal to zero)
	// ExchCoeffVD - water vapor exchange coefficient
	tval = TepVal;
	WDVal = (M_PI/180.0)*WDVal;
	if (fabs (basicdirect - WDVal) <= M_PI){
	  betaWind = fabs (basicdirect-WDVal);
	}
	else{
	  betaWind = 2*M_PI - fabs(basicdirect-WDVal);
	}

	if (betaWind < (M_PI/2)){
	  betaExch = ExchCoeffVD*(1+WindFactVD * sqrt(WSVal));
	}
	else{
	  betaExch = ExchCoeffVD;
	}
	break;
  }
    //case 66:{
    // water vapor exchange coefficient is given by table
    //betaExch = exchangecoeff.getval (MoistureAkt);
    //break;
    //}
  case 67:{
    // ExchCoeffVD - water vapor exchange coefficient
    betaExch = ExchCoeffVD;
    tval = rhtemp + 273.15;
    break;
  }
  default:{
    print_err("unknown type of VDKind is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (VDKind)
  
  c = cvap*((tval-temp)/2.0-T0) + hvap;

//  fprintf (stdout,"\n tw   %le",betaExch*c);
  
  return betaExch*c;
}

// *****************************************************************************
// *****************************************************************************
//
//  MOISTURE FLUXES, VALUES AND COEFFICIENTS
//
// *****************************************************************************
// *****************************************************************************


/**
   function computes the moisture flux density

   @param time - actual time

   10. 4. 2014
*/
double climatcond2::give_moisture_flux (double /*nval_0*/, double /*nval_1*/,long /*ipp*/, long nid)
{
  /*
   nval_0 - rh for Kunzel type[1]
	  - pv for Kunzel type[2-4]
	  - w for Grunewaldmat, salt1-4, millymat,

   nval_1 - temperature
  */
  double jh4, jh5;


  switch (KIND4){
  case 0:{
	//  KIND4 = 0 - no water contact
	jh4 = 0.0;
	break;
  }
  case 41:{
    jh4 = 0.0;
	break;
  }
  default:{
	print_err("unknown definition of Get_moisture_flux-KIND4 is required",__FILE__,__LINE__,__func__);
  }
  }// end of switch(KIND4)

  switch (KIND5){
  case 0:{
	//  KIND5 =  0 - no rain
	jh5 = 0.0;
	break;
  }
  case 51:{
    //  KIND5 = 51 - prescribed water flux caused by rain
    jh5 = -1.0 * condit_rain_moisture_flux (51,nid);
	break;
  }
  case 52:{
	//  KIND5 = 52 - prescribed water flux caused by rain extended by wind effect
	jh5 = -1.0 * condit_rain_moisture_flux (52,nid);
	break;
  }
  default:{
	print_err("unknown definition of Get_moisture_flux-KIND5 is required",__FILE__,__LINE__,__func__);
  }
  }// end of switch(KIND5)

  return (jh4 +jh5);
}


/**
   function computes mass of convective liqiud water flux

   @param rainkind - type of flux computation
   @param nid - node id
   

   JM, revised 5. 10. 2013
*/
double climatcond2::condit_rain_moisture_flux (long rainkind,long nid)
{
  double NorRain, MoistureAkt, Smc;
  double kwind, betawind, krainEff, jkNorRain, jkmaxWat;
  double c1, c2, c3, c4, c5;
 
  //  returns volumetric moisture content stored in nodes
  MoistureAkt = Tm->give_nodal_vol_moist_cont (nid,mattyp);
 
  //  returns saturated volumetric moisture content stored in nodes
  Smc = Tm->give_nodal_sat_vol_moist_cont (nid,mattyp);
  
  //  toto nema cenu, je to double
  if (Smc == 0.0)
    return 0.0;
  
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
//      fprintf (stdout,"\n wsval   %le",WSVal);

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
   
   JK 10. 4. 2014
*/
double climatcond2::transmission_coefficient_ww (long VDKind)
{

  double betaWind,betaExch;


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
  case 0:{
	betaExch = 0.0;
	break;
  }
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
//       fprintf (stdout,"\n WSVALs   %le",WSVal);

	if (betaWind < (M_PI/2)){
	  betaExch = ExchCoeffVD*(1 +WindFactVD * sqrt(WSVal));
	}
	else{
	  betaExch = ExchCoeffVD;
	}
	break;
  }
    //case 66:{
    // water vapor exchange coefficient is given by table
    //betaExch = exchangecoeff.getval (MoistureAkt);
    //break;
    //}
  case 67:{
    // ExchCoeffVD - water vapor exchange coefficient
    betaExch = ExchCoeffVD;
    break;
  }
  default:{
    print_err("unknown type of VDKind is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (VDKind)


//  fprintf (stdout,"\n betaExch   %le",betaExch);

  return betaExch;
}

/**
   function computes water vapour pressure in exterior
   it is used in the Newton boundary condition
   
   @param VDKind - type of mass of diffusive water vapour flux
   
   JK 11. 4. 2014
*/
double climatcond2::transmission_exterior_value_ww (long VDKind)
{
  double vappres,psate, tval;
  tval = 298.15;
  
  switch (VDKind){
  case 0:{
    break;
  }
  case 61:{
//    tval = TepVal;
//    tval = 298.15;
    break;
  }
  case 62:{
//    tval = TepVal;
//    tval = 298.15;
    break;
  }
  case 67:{
    // rhtemp - constant temperature of ambiente air
    tval = rhtemp + 273.15;
    RHVal=rhint/100;
    break;
  }
  default:{
    print_err("unknown type of VDKind is required",__FILE__,__LINE__,__func__);
  }
  }//  end of switch (VDKind)
  
  
  //external water vapor pressure
  psate = exp(23.5771 - 4042.9/(tval - 37.58));
  
  vappres = psate * RHVal;

//  fprintf (stdout,"\n vappres   %le",vappres);

  return vappres;
}

/**
   function returns external value of the cid-th variable
   it is used in the Newton boundary condition
   
   @param cid - component id
   
   JK 11. 4. 2014
*/
double climatcond2::external_nodval (long cid)
{
  double eval;
  
  //heat transfer
  switch (Tp->mednam){//names of transported media
  case heat:{
    eval = transmission_exterior_value_tt (KIND1);
    break;
  }
  case moisture:
  case water:
  case heat_moisture:{
    //other transfer analyses
    if (cid==0){
      //  moisture or pressure
      eval = transmission_exterior_value_ww (KIND6);    
    }
    if (cid==1){
      //  temperature
      eval = transmission_exterior_value_tt (KIND1);
    }
    break;
  }
  default:{
    print_err("\n Unknown media name is required in ", __FILE__, __LINE__, __func__);
  }
  } 
  
  //  fprintf (stdout,"\n eval   %le",eval);
  
  return eval;
}

/**
   function returns transmission coefficient for the lcid-th flux caused by the cid-th gradient
   it is used in the Newton boundary condition
   
   @param lcid - load case id
   @param cid - component id
   
   JK 11. 4. 2014
*/
double climatcond2::transmission_coeff (long lcid,long cid,double temp)
{
  double trcoef;

  trcoef=0.0;
  
  //heat transfer
  switch (Tp->mednam){//names of transported media
  case heat:{
    trcoef = transmission_coefficient_tt (KIND1);
    break;
  }
  case moisture:
  case water:
  case heat_moisture:{
    //other transfer analyses
    if (lcid==0){
      if (cid==0){
	//  moisture flux caused by the moisture gradient
	trcoef = transmission_coefficient_ww (KIND6);
      }
      if (cid==1){
	//  moisture flux caused by the temperature gradient
	//  there is no coupling
      }
    }
    if (lcid==1){
      if (cid==0){
	//  heat flux cause by the moisture gradient
	trcoef = transmission_coefficient_tw (KIND6,temp);
      }
      if (cid==1){
	//  heat flux cause by the temperature gradient
	trcoef = transmission_coefficient_tt (KIND1);
      }
    }
    break;
  }
  default:{
    print_err("\n Unknown media name is required in ", __FILE__, __LINE__, __func__);
  }
  } 
  
  //  fprintf (stdout,"\n trcoef   %le",trcoef);

  return trcoef;
}
