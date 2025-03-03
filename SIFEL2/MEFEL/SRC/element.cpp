#include <string.h>
#include "element.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "elemhead.h"
#include "gtopology.h"
#include "intpoints.h"



/**
  Constructor initializes data members to zero or default values.

  Created by JK,
*/
element::element (void)
{
  //  material type
  tm = NULL;
  //  material id
  idm = NULL;
  //  number of material types
  nm = 0;
  //  type of element
  te = (elemtype) 0;
  
  
  //  indicator of prescribed displacements on element
  prescdispl=0;
  //  indicator of temperature changes on element
  presctemp=0;
  //  computation of reactions
  react=0;
  //  array of integration point ids
  ipp=NULL;
  //  type of cross section
  crst = (crsectype) 0;
  //  cross section id
  idcs=0;
  //  stress/strain state
  ssst = (strastrestate)0;
  //  number of strain/stress components
  ncomp=0;
  //  array of initial displacements
  initdispl=NULL;
  //  transformation matrix due to hanging nodes
  tmat=NULL;
}



/**
  Destructor releases allocated memory of the element object.

  Created by JK,
*/
element::~element (void)
{
  long i,nb;
  
  delete [] tm;
  delete [] idm;
  delete [] initdispl;
  
  //  the number of blocks
  nb = Mt->give_nb_te (te);
  
  if (ipp != NULL){
    for (i=0; i < nb; i++){
      delete [] ipp[i];
    }
    delete [] ipp;
  }
  delete tmat;
}



/**
  Function reads input data about elements from the opened text file.
   
  @param in - pointer to the opened XFILE
  @param eid - element id

  @return The function does not return anything.

  Created by JK,
*/
void element::read (XFILE *in,long eid)
{
  long nne,ndofe;
  gelemtype get=noelem;
  
  //  element type
  xfscanf (in,"%m", &elemtype_kwdset, (int*)&te);
  
  switch (te){
  case bar2d:{                if (Bar2d==NULL)     Bar2d = new barel2d;                break;  }
  case bar3d:{                if (Bar3d==NULL)     Bar3d = new barel3d;                break;  }
  case barq2d:{               if (Barq2d==NULL)    Barq2d = new barelq2d;              break;  }
  case barq3d:{               if (Barq3d==NULL)    Barq3d = new barelq3d;              break;  }
  case beam2d:{               if (Beam2d==NULL)    Beam2d = new beamel2d;              break;  }
  case beam3d:{               if (Beam3d==NULL)    Beam3d = new beamel3d;              break;  }
  case beamg3d:{              if (Beam3dg==NULL)   Beam3dg = new beamgen3d;            break;  }
  case subsoilbeam:{          if (Sbeam==NULL)     Sbeam = new soilbeam;               break;  }
  case spring_1:{             if (Spring==NULL)    Spring = new springel;              break;  }
  case spring_2:{             if (Spring==NULL)    Spring = new springel;              break;  }
  case spring_3:{             if (Spring==NULL)    Spring = new springel;              break;  }
  case spring_4:{             if (Spring==NULL)    Spring = new springel;              break;  }
  case spring_5:{             if (Spring==NULL)    Spring = new springel;              break;  }
  case spring_6:{             if (Spring==NULL)    Spring = new springel;              break;  }
  case planeelementlt:{       if (Pelt==NULL)      Pelt = new planeelemlt;             break;  }
  case planeelementqt:{       if (Peqt==NULL)      Peqt = new planeelemqt;             break;  }
  case planeelementrotlt:{    if (Perlt==NULL)     Perlt = new planeelemrotlt;         break;  }
  case planeelementlq:{       if (Pelq==NULL)      Pelq = new planeelemlq;             break;  }
  case planeelementqq:{       if (Peqq==NULL)      Peqq = new planeelemqq;             break;  }
  case planeelementrotlq:{    if (Perlq==NULL)     Perlq = new planeelemrotlq;         break;  }
  case planeelementsubqt:{    if (Pesqt==NULL)     Pesqt = new planeelemsubqt;         break;  }
  case planequadinterface:{   if (Pqifc==NULL)     Pqifc = new plquadinterface;        break;  }
  case cctel:{                if (Cct==NULL)       Cct = new cctelem;                  break;  }
  case dktel:{                if (Dkt==NULL)       Dkt = new dktelem;                  break;  }
  case dstel:{                if (Dst==NULL)       Dst = new dstelem;                  break;  }
  case q4plateel:{            if (Q4pl==NULL)      Q4pl = new q4plate;                 break;  }
    //case argyristr:{            if (Argtr==NULL)     Argtr = new ArgyrisTriangle;   break;  }
  case argyristr:{            if (Argtrpl==NULL)   Argtrpl = new argyrisplate;         break;  }
  case quadkirch:{            if (Qkirch==NULL)    Qkirch = new quadrilatkirch;        break;  }
  case dkqel:{                if (Dkqelem==NULL)   Dkqelem = new dkq;                  break;  }
  case subsoilplatetr:{       if (Spltr==NULL)     Spltr = new soilplatetr;            break;  }
  case subsoilplateq:{        if (Splq==NULL)      Splq = new soilplateq;              break;  }
  case axisymmlt:{            if (Asymlt==NULL)    Asymlt = new axisymlt;              break;  }
  case axisymmqt:{            if (Asymqt==NULL)    Asymqt = new axisymqt;              break;  }
  case axisymmlq:{            if (Asymlq==NULL)    Asymlq = new axisymlq;              break;  }
  case axisymmqq:{            if (Asymqq==NULL)    Asymqq = new axisymqq;              break;  }
  case axisymmcq:{            if (Asymcq==NULL)    Asymcq = new axisymcq;              break;  }
  case axisymmlqintface:{     if (Asymlqifc==NULL) Asymlqifc = new axisymlqinterface;  break;  }
  case shelltrelem:{          if (Shtr==NULL)      Shtr = new shelltr;                 break;  }
  case shelltrmelem:{         if (Shtrm==NULL)     Shtrm = new shelltrm;               break;  }
  case shellqelem:{           if (Shq==NULL)       Shq = new shellq;                   break;  }
  case lineartet:{            if (Ltet==NULL)      Ltet = new lintet;                  break;  }
  case quadrtet:{             if (Qtet==NULL)      Qtet = new quadtet;                 break;  }
  case linearhex:{            if (Lhex==NULL)      Lhex = new linhex;                  break;  }
  case quadrhex:{             if (Qhex==NULL)      Qhex = new quadhex;                 break;  }
  case lineartetrot:{         if (Ltetrot==NULL)   Ltetrot = new lintetrot;            break;  }
  case linearhexrot:{         if (Lhexrot==NULL)   Lhexrot = new linhexrot;            break;  }
  case linearwed:{            if (Lwed==NULL)      Lwed = new linwedge;                break;  }
  case quadrwed:{             if (Qwed==NULL)      Qwed = new quadwedge;               break;  }
  case hexintface:{           if (Hexifc==NULL)    Hexifc = new hexinterface;          break;  }
  case particleelem:{
    if (Pelem==NULL){
      long nne,dim;
      xfscanf (in,"%ld %ld",&nne,&dim);
      Pelem = new elemparticle (nne,dim);
    }
    break;
  }
  case tetralatt:{            if (Tlatt==NULL)     Tlatt = new tetralattice;      break;  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
  
  switch (te){
  case planeelementlt:{     xfscanf (in,"%m", &strastrestate_kwdset, (int*)&ssst);  break;  }
  case planeelementqt:{     xfscanf (in,"%m", &strastrestate_kwdset, (int*)&ssst);  break;  }
  case planeelementrotlt:{  xfscanf (in,"%m", &strastrestate_kwdset, (int*)&ssst);  break;  }
  case planeelementlq:{     xfscanf (in,"%m", &strastrestate_kwdset, (int*)&ssst);  break;  }
  case planeelementqq:{     xfscanf (in,"%m", &strastrestate_kwdset, (int*)&ssst);  break;  }
  case planeelementrotlq:{  xfscanf (in,"%m", &strastrestate_kwdset, (int*)&ssst);  break;  }
  case planeelementsubqt:{  xfscanf (in,"%m", &strastrestate_kwdset, (int*)&ssst);  break;  }
  default:{
  }
  }

  switch (te){
  case bar2d:{                get=linbar;      break;  }
  case bar3d:{                get=linbar;      break;  }
  case barq2d:{               get=linbar;      break;  }
  case barq3d:{               get=linbar;      break;  }
  case beam2d:{                                break;  }
  case beam3d:{                                break;  }
  case beamg3d:{                               break;  }
  case subsoilbeam:{                           break;  }
  case spring_1:{                              break;  }
  case spring_2:{                              break;  }
  case spring_3:{                              break;  }
  case spring_4:{                              break;  }
  case spring_5:{                              break;  }
  case spring_6:{                              break;  }
  case planeelementlt:{       get=lintriag;    break;  }
  case planeelementqt:{       get=quadtriag;   break;  }
  case planeelementrotlt:{    get=lintriag;    break;  }
  case planeelementlq:{       get=linquad;     break;  }
  case planeelementqq:{       get=quadquad;    break;  }
  case planeelementrotlq:{    get=linquad;     break;  }
  case planeelementsubqt:{                     break;  }
  case planequadinterface:{                    break;  }
  case cctel:{                get=lintriag;    break;  }
  case dktel:{                get=linquad;     break;  }
  case dstel:{                get=lintriag;    break;  }
  case q4plateel:{                             break;  }
  case argyristr:{                             break;  }
  case quadkirch:{                             break;  }
  case dkqel:{                                 break;  }
  case subsoilplatetr:{                        break;  }
  case subsoilplateq:{                         break;  }
  case axisymmlt:{            get=lintriag;    break;  }
  case axisymmqt:{            get=quadtriag;   break;  }
  case axisymmlq:{            get=linquad;     break;  }
  case axisymmqq:{            get=quadquad;    break;  }
  case axisymmcq:{            get=cubicquad;   break;  }
  case axisymmlqintface:{                      break;  }
  case shelltrelem:{                           break;  }
  case shelltrmelem:{                          break;  }
  case shellqelem:{                            break;  }
  case lineartet:{            get=lintetra;    break;  }
  case quadrtet:{             get=quadtetra;   break;  }
  case linearhex:{            get=linhexa;     break;  }
  case quadrhex:{             get=quadhexa;    break;  }
  case lineartetrot:{         get=lintetra;    break;  }
  case linearhexrot:{         get=linhexa;     break;  }
  case linearwed:{                             break;  }
  case quadrwed:{                              break;  }
  case hexintface:{                            break;  }
  case particleelem:{                          break;  }
  case tetralatt:{                             break;  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
  
  ndofe=Mt->give_ndofe (eid);
  nne = Mt->give_nne (eid);

  /* // this block was moved to mefelinit
  
  if ((Mp->homog == 3) || (Mp->homog == 5) || (Mp->homog == 7) || (Mp->homog == 9)){
    //  function modifies the number of degrees of freedom in the case of
    //  homogenization based on prescribed stresses
    ndofe+= Mt->give_tncomp(eid);
  }*/
  
  if (Mp->tprob != growing_mech_structure){
    Gtm->gelements[eid].read(in, nne, ndofe, get, Mt->nn);
  }
  else{
    Gtm->gelements[eid].read_gf(in, nne, ndofe, get, Mt->nn);
  }

  //////////////////////////////////
  //only for temelin quadtet
  switch (te){
  case quadrtet:{
    double vol;
    vol = Mt->give_volume(eid);
    if(vol < 0.0){
      printf("element volume for quadtet = %ld is negative\n\n",eid+1);
      print_err("wrong numbering in mesh element", __FILE__, __LINE__, __func__);
      exit(0);
    }
    
    break;}
  default:{}
  }
  //////////////////////////////////

  xfscanf (in,"%m", &crsectype_kwdset, (int*)&crst);
  if (crst!=0){
    xfscanf (in,"%ld",&idcs);
    idcs--;
  }
  
  readmat (in);
  
}




/**
  Function prints input data about elements into the opened text file.
   
  @param out - pointer to the opened FILE
  @param eid - element id

  @return The function does not return anything.

  Created by TKr, 02/01/2013
*/
void element::print (FILE *out,long eid)
{
  long nne,ndofe;
  gelemtype get;
  
  //  element type
  fprintf (out,"\n %d",te);
  
  //  the number of DOFs on element
  ndofe=Mt->give_ndofe (eid);
  //  the number of nodes on element
  nne = Mt->give_nne (eid);
  //  type of element (see galias.h)
  get = Gtm->give_elem_type (eid);
  
  if (Mp->tprob != growing_mech_structure){
    Gtm->gelements[eid].print (out,nne,ndofe,get);
  }
  else{
    Gtm->gelements[eid].print_gf (out,nne,ndofe);
  }
  
  fprintf (out,"  %d",crst);
  if (crst!=0){
    fprintf (out," %ld",idcs+1);
  }
  
  printmat (out);  
}



/**
  The function reads materials assigned to the particular integration 
  points of elements.

  Parameters:
  @param in - pointer to the opened text XFILE

  @return The function does not return anything.

  Created by JK, TKo,
*/
void element::readmat (XFILE *in)
{
  long i, idx, idxt;

  xfscanf (in, "%ld", &nm);
  tm  = new mattype[nm];
  idm = new long[nm];
  memset (tm,  0, sizeof(*tm)*nm);
  memset (idm, 0, sizeof(*idm)*nm);
  idx = 0;
  idxt = nm-1;
  for (i = 0; i < nm; i++)
  {
    xfscanf(in, "%k%m%k%ld", "mattype", &mattype_kwdset, (int*)(tm+idx), "num_inst", idm+idx);
    //xfscanf (in, "%d %ld", (int*)(tm+idx), idm+idx);
    idm[idx]--;
    /*    if (tm[idx] == therisodilat){
      tm[idxt] = tm[idx];
      idm[idxt] = idm[idx];
      idxt--;
      idx--;
      }*/
    idx++;
  }
}


/**
  The function prints materials assigned to the particular integration 
  points of elements.

  Parameters:
  @param out - pointer to the opened text FILE

  @return The function does not return anything.

  Created by TKr, 02/01/2013
*/
void element::printmat (FILE *out)
{
  long i;

  fprintf (out, " %ld", nm);
  for (i = 0; i < nm; i++)
    {
      fprintf(out, " %d %ld", tm[i],idm[i]+1);
    }
  fprintf (out,"\n");
}



/**
  Function allocates array for initial displacements.
  Initial displacements are used in problems with growing number of elements.
   
  @praram ndofe - number of DOFs on element
   
  Created by JK, 3.3.2006
*/
void element::alloc_initdispl (long ndofe)
{
  initdispl = new double [ndofe];
}



/**
  Function defines initial displacements
  initial displacements are used in problems with growing number of elements
   
  @param r - array containing initial displacements
  @param ndofe - number of DOFs on element
   
  Created by JK, 3.3.2006
*/
void element::initdisplacement (double *r,long ndofe)
{
  long i;
  
  for (i=0;i<ndofe;i++){
    initdispl[i]=r[i];
  }
  
}



/**
  Function subtracts initial displacements from total displacements.
  Function is used in problems with changing number of elements
   
  @param r - total displacement
  @param ndofe - number of DOFs on element
   
  Created by JK, 5.3.2006
*/
void element::subtrinitdispl (double *r,long ndofe)
{
  long i;
  
  for (i=0;i<ndofe;i++){
    r[i]-=initdispl[i];
  }
}



/**
  Function allocates array used for problems with
  changing nodes, elements and DOFs.
   
  @param eid - element id
   
  Created by JK, 7.11.2006
*/
void element::alloc_growstr (long eid)
{
  long i,ndofe;
  
  ndofe=Mt->give_ndofe (eid);

  initdispl = new double [ndofe];

  for (i=0;i<ndofe;i++){
    initdispl[i]=0.0;
  }
}

