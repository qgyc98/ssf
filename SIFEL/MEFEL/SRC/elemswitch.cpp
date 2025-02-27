#include "elemswitch.h"
#include "ipmap.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "mechcrsec.h"
#include "mechbclc.h"
#include "lhsrhs.h"
#include "elemhead.h"
#include "node.h"
#include "element.h"
#include "intpoints.h"
#include "loadcase.h"
#include "dloadcase.h"
#include "globmat.h"
#include "globalg.h"
//#undef INC_OPENMP
#ifdef INC_OPENMP
 #include <omp.h>
#endif

// *******************
// *******************
// ****  STRAINS  ****
// *******************
// *******************



/**
  The function computes strains at integration points.   
  All components at all integration points are evaluated.
   
  @param lcid - load case id
   
  @return The function does not return anything.

  Created by JK, 28.11.2006
*/
void compute_ipstrains (long lcid)
{
  long i,j,ne,te;
  long ipp, ncomp;
  // long nmatpoints;

  //  strain array cleaning
  Mm->cleanstrain ();
  
  //  number of elements
  ne=Mt->ne;
#ifdef INC_OPENMP  
#pragma omp parallel \
            num_threads(Numth),\
            private(i, te)

  {    
#pragma omp for
#endif
  for (i=0;i<ne;i++){
    if (Gtm->leso[i]==1){
      //  only elements switched on are taken into account
      
      //  type of element
      te = Mt->give_elem_type (i);
      
      switch (te){
	
      case bar2d:{
	Bar2d->res_ip_strains (lcid,i);
	break;
      }
      case bar3d:{
	Bar3d->res_ip_strains (lcid,i);
	break;
      }
      case barq2d:{
	Barq2d->res_mainip_strains (lcid,i);
	break;
      }
      case barq3d:{
	Barq3d->res_mainip_strains (lcid,i);
	break;
      }
	
      case beam2d:{
	Beam2d->nodal_displ (lcid,i);
	break;
      }
      case beam3d:{
	Beam3d->nodal_displ (lcid,i);
	break;
      }
      case beamg3d:{
	Beam3dg->nodal_displ (i,lcid);
	break;
      }
      case subsoilbeam:{
	Sbeam->strains (lcid, i, 0, 0);
	break;
      }
	
      case spring_1:{
	Spring->strains(i, lcid);
	break;
      }
      case spring_2:{
	Spring->strains(i, lcid);
	break;
      }
      case spring_3:{
	Spring->strains(i, lcid);
	break;
      }
      case spring_4:{
	Spring->strains(i, lcid);
	break;
      }
      case spring_5:{
	Spring->strains(i, lcid);
	break;
      }
      case spring_6:{
	Spring->strains(i, lcid);
	break;
      }
	
	
      case planeelementlt:{
	Pelt->res_ip_strains (lcid,i);
	break;
      }
      case planeelementqt:{
	Peqt->res_ip_strains (lcid,i);
	break;
      }
	
      case planeelementrotlt:{
	Perlt->res_ip_strains (lcid,i);
	break;
      }
	
      case planeelementlq:{
	Pelq->res_ip_strains (lcid,i);
	break;
      }
      case planeelementqq:{
	Peqq->res_ip_strains (lcid,i);
	break;
      }
      case planeelementrotlq:{
	Perlq->res_ip_strains (lcid,i);
	break;
      }
	
      case planeelementsubqt:{
	Pesqt->res_mainip_strains (lcid,i);
	break;
      }
      case planequadinterface:{
	Pqifc->res_mainip_strains (lcid,i);
	break;
      }
	
      case cctel:{
	Cct->res_ip_strains (lcid,i);
	break;
      }
      case dktel:{
	//Dkt->res_mainip_strains (lcid,i);
	Dkt->res_ip_strains (lcid,i);
	break;
      }
      case dstel:{
	Dst->res_ip_strains (lcid,i);
	break;
      }
      case q4plateel:{
	Q4pl->res_ip_strains (lcid,i);
	break;
      }
      case dkqel:{
	Dkqelem->res_ip_curvatures (lcid,i);
	break;
      }
	
      case subsoilplatetr:{
	Spltr->res_mainip_strains (lcid,i);
	break;
      }
	
      case subsoilplateq:{
	Splq->res_mainip_strains (lcid,i);
	break;
      }
	
	
      case shelltrelem:{
	Shtr->res_ip_strains (lcid,i);
	break;
      }
      case shelltrmelem:{
	Shtrm->res_ip_strains (lcid,i);
	break;
      }
      case shellqelem:{
	Shq->res_ip_strains (lcid,i);
	break;
      }
	
	
      case axisymmlt:{
	Asymlt->res_mainip_strains (lcid,i);
	break;
      }
      case axisymmqt:{
	Asymqt->res_mainip_strains (lcid,i);
	break;
      }
      case axisymmlq:{
	Asymlq->res_allip_strains (lcid,i);
	break;
      }
      case axisymmqq:{
	Asymqq->res_allip_strains (lcid,i);
	break;
      }
      case axisymmcq:{
	Asymcq->res_allip_strains (lcid,i);
	break;
      }
      case axisymmlqintface:{
        Asymlqifc->res_mainip_strains(lcid, i);
        break;
      }	
	
      case lineartet:{
        
	Ltet->res_ip_strains (lcid,i);
	break;
      }
      case quadrtet:{
	Qtet->res_ip_strains (lcid,i);
	break;
      }
      case linearhex:{
	Lhex->res_ip_strains (lcid,i);
	break;
      }
      case quadrhex:{
	Qhex->res_ip_strains (lcid,i);
	break;
      }
      case lineartetrot:{
	Ltetrot->res_ip_strains (lcid,i);
	break;
      }
      case linearhexrot:{
	Lhexrot->res_ip_strains (lcid,i);
	break;
      }
	
      case hexintface:{
        Hexifc->res_mainip_strains(lcid, i);
        break;
      }
	
      case particleelem:{
	break;
      }
      case tetralatt:{
	Tlatt->ip_strains (lcid,i,0,0);
	break;
      }
	
      default:{
	print_err("unknown element type is required", __FILE__, __LINE__, __func__);
      }
      }
      
    }
  }
#ifdef INC_OPENMP
  }
#endif
  //  function modifies strains
  //  eigenstrains are subtracted from the total strains
  //  it works only for the first load case
  Mm->eigstrmod ();
  
  if (Mb->give_num_mstress_comp(lcid)){
    //  homogenization, stress approach
    //  element functions compute fluctuation strains only
    //  the following function adds macroscopic strains computed in the system of equations
    
    ncomp = Mt->max_ncompstr;
    double *lhs = Lsrs->give_lhs (0);
    vector macrostrains (ncomp);
    strastre *mstrastre = Mb->give_mstrastre(lcid);
    j=Ndofm-Mb->give_num_mstress_comp(lcid);
    for (i=0;i<ncomp;i++){
      if (mstrastre[i] == stress){
        macrostrains[i]=lhs[j];
        j++;
      }
      else
        macrostrains[i] = 0.0;
    }
    Mm->add_macro_strains (lcid, macrostrains);
  }
  
  if (Mb->give_num_mstrain_comp(lcid)){
    
    //  homogenization, strain approach
    //  element functions compute fluctuation strains only
    //  the following function adds macroscopic strains computed in the system of equations
    ncomp = Mt->max_ncompstr;
    vector macrostrains(ASTCKVEC(ncomp));
    switch(Mp->tprob)
    {
      case linear_statics:
        copyv(Mb->lc[lcid].mstrain, macrostrains);
        Mm->add_macro_strains (lcid, macrostrains);
        break;
      case mat_nonlinear_statics:
        addmultv(Mb->lc[lcid+1].mstrain, Mb->lc[lcid].mstrain, Mp->lambda, macrostrains.a, ncomp);
        Mm->add_macro_strains (lcid, macrostrains);
        break;
      case mech_timedependent_prob:
        for(i=0; i<ncomp; i++)
          macrostrains(i) = Mb->dlc[lcid].get_macrostra(Mp->time, i);
        Mm->add_macro_strains (lcid, macrostrains);
        break;
      default:
        print_err("unknown type of problem is required", __FILE__, __LINE__, __func__);
    }
  }
  if (lcid != 0)
  {
    for (ipp=0; ipp<Mm->tnip; ipp++)
    {
      ncomp = Mm->ip[ipp].ncompstr;
      if (Gtm->leso[Mm->elip[ipp]]==1){
        memcpy(Mm->ip[ipp].strain, Mm->ip[ipp].strain+ncomp*lcid, sizeof(Mm->ip[ipp].strain)*ncomp);
      }
    }
  }
}

/**
  The function computes strains at material points
  the term material points is used for any point, where the strains are computed
  it can be integration point in MEFEL, integration point in TRFEL, a point defined by an user, etc.
  
  there can be several heaps of material points in the code,
  in such a case, this function is called for every heap of points independently
  this function cannot compute strains for all heaps within one call
  

  poznamka pro vyvojare: materialove body jsou setrideny podle te casti programu, ve ktere vznikly
  jsou-li body vygenerovany v TRFELu, budou obecne setrideny tak, ze pri jejich prochazeni
  v MEFELu se bude skakat chaoticky po mechanickych prvcich. Proto je treba mit pomocne pole,
  s jehoz pomoci bude mozne probirat body podle mechanickych prvku. Motivaci je, aby se
  vyridily vsechny body na jednotlivych mechanickych prvcich naraz.

  
  @param lcid - load case id
  @param matpoints - array of material points
  @param mpadr - cumulative number of material points on element
                 madr[i]=j - there are j material points on elements before the i-th one
  @param matpelemmap - map between material points and elements
                       it contains the number of components equal to the number of all points in a heap
         matpelemmap[i]=j 

  @return The function does not return anything.

  Created by TKo, JK, 25. 3. 2021
*/
void compute_strains (long lcid,intpoints *matpoints,long *mpadr,long *matpelemmap)
{
  long i,j,ne,te,ipp;
  long ncomp,nallmatpoints;

  //  number of elements
  ne=Mt->ne;
  
  //  number of all amterial points
  nallmatpoints = mpadr[ne];
  
  //  strain array cleaning
  for (i=0;i<nallmatpoints;i++){
    //matpoints[i].clean_strain (lcid);
  }
  
  for (i=0;i<ne;i++){
    if (Gtm->leso[i]==1){
      //  only elements switched on are taken into account
      
      
      
      
      //  type of element
      te = Mt->give_elem_type (i);
      
      switch (te){
	
      case bar2d:{
	Bar2d->res_ip_strains (lcid,i);
	break;
      }
      case bar3d:{
	Bar3d->res_ip_strains (lcid,i);
	break;
      }
      case barq2d:{
	Barq2d->res_mainip_strains (lcid,i);
	break;
      }
      case barq3d:{
	Barq3d->res_mainip_strains (lcid,i);
	break;
      }
	
      case beam2d:{
	Beam2d->nodal_displ (lcid,i);
	break;
      }
      case beam3d:{
	Beam3d->nodal_displ (lcid,i);
	break;
      }
      case beamg3d:{
	Beam3dg->nodal_displ (i,lcid);
	break;
      }
      case subsoilbeam:{
	Sbeam->strains (lcid, i, 0, 0);
	break;
      }
	
      case spring_1:{
	Spring->strains(i, lcid);
	break;
      }
      case spring_2:{
	Spring->strains(i, lcid);
	break;
      }
      case spring_3:{
	Spring->strains(i, lcid);
	break;
      }
      case spring_4:{
	Spring->strains(i, lcid);
	break;
      }
      case spring_5:{
	Spring->strains(i, lcid);
	break;
      }
      case spring_6:{
	Spring->strains(i, lcid);
	break;
      }
	
	
      case planeelementlt:{
	Pelt->res_ip_strains (lcid,i);
	break;
      }
      case planeelementqt:{
	Peqt->res_ip_strains (lcid,i);
	break;
      }
	
      case planeelementrotlt:{
	Perlt->res_ip_strains (lcid,i);
	break;
      }
	
      case planeelementlq:{
	Pelq->res_ip_strains (lcid,i);
	break;
      }
      case planeelementqq:{
	Peqq->res_ip_strains (lcid,i);
	break;
      }
      case planeelementrotlq:{
	Perlq->res_ip_strains (lcid,i);
	break;
      }
	
      case planeelementsubqt:{
	Pesqt->res_mainip_strains (lcid,i);
	break;
      }
      case planequadinterface:{
	Pqifc->res_mainip_strains (lcid,i);
	break;
      }
	
      case cctel:{
	Cct->res_ip_strains (lcid,i);
	break;
      }
      case dktel:{
	//Dkt->res_mainip_strains (lcid,i);
	Dkt->res_ip_strains (lcid,i);
	break;
      }
      case dstel:{
	Dst->res_ip_strains (lcid,i);
	break;
      }
      case q4plateel:{
	Q4pl->res_ip_strains (lcid,i);
	break;
      }
      case dkqel:{
	Dkqelem->res_ip_curvatures (lcid,i);
	break;
      }
	
      case subsoilplatetr:{
	Spltr->res_mainip_strains (lcid,i);
	break;
      }
	
      case subsoilplateq:{
	Splq->res_mainip_strains (lcid,i);
	break;
      }
	
	
      case shelltrelem:{
	Shtr->res_ip_strains (lcid,i);
	break;
      }
      case shelltrmelem:{
	Shtrm->res_ip_strains (lcid,i);
	break;
      }
      case shellqelem:{
	Shq->res_ip_strains (lcid,i);
	break;
      }
	
	
      case axisymmlt:{
	Asymlt->res_mainip_strains (lcid,i);
	break;
      }
      case axisymmqt:{
	Asymqt->res_mainip_strains (lcid,i);
	break;
      }
      case axisymmlq:{
	Asymlq->res_allip_strains (lcid,i);
	break;
      }
      case axisymmqq:{
	Asymqq->res_allip_strains (lcid,i);
	break;
      }
      case axisymmcq:{
	Asymcq->res_allip_strains (lcid,i);
	break;
      }
      case axisymmlqintface:{
        Asymlqifc->res_mainip_strains(lcid, i);
        break;
      }	
	
	
      case lineartet:{
	Ltet->res_ip_strains (lcid,i);
	break;
      }
      case quadrtet:{
	Qtet->res_ip_strains (lcid,i);
	break;
      }
      case linearhex:{
	Lhex->res_ip_strains2 (lcid,i,matpoints,mpadr,matpelemmap);
	
	break;
      }
      case quadrhex:{
	Qhex->res_ip_strains (lcid,i);
	break;
      }
      case lineartetrot:{
	Ltetrot->res_ip_strains (lcid,i);
	break;
      }
      case linearhexrot:{
	Lhexrot->res_ip_strains (lcid,i);
	break;
      }
      case hexintface:{
        Hexifc->res_mainip_strains(lcid, i);
        break;
      }
      case particleelem:{
	break;
      }
      case tetralatt:{
	Tlatt->ip_strains (lcid,i,0,0);
	break;
      }
	
      default:{
	print_err("unknown element type is required", __FILE__, __LINE__, __func__);
      }
      }
      
      
    }
  }

  //  function modifies strains
  //  eigenstrains are subtracted from the total strains
  //  it works only for the first load case
  Mm->eigstrmod ();
  
  if (Mb->give_num_mstress_comp(lcid)){
    //  homogenization, stress approach
    //  element functions compute fluctuation strains only
    //  the following function adds macroscopic strains computed in the system of equations
    
    ncomp = Mt->max_ncompstr;
    double *lhs = Lsrs->give_lhs (0);
    vector macrostrains (ncomp);
    strastre *mstrastre = Mb->give_mstrastre(lcid);
    j=Ndofm-Mb->give_num_mstress_comp(lcid);
    for (i=0;i<ncomp;i++){
      if (mstrastre[i] == stress){
        macrostrains[i]=lhs[j];
        j++;
      }
      else
        macrostrains[i] = 0.0;
    }
    Mm->add_macro_strains (lcid, macrostrains);
  }
  
  if (Mb->give_num_mstrain_comp(lcid)){
    
    //  homogenization, strain approach
    //  element functions compute fluctuation strains only
    //  the following function adds macroscopic strains computed in the system of equations
    ncomp = Mt->max_ncompstr;
    vector macrostrains(ASTCKVEC(ncomp));
    switch(Mp->tprob)
    {
      case linear_statics:
        copyv(Mb->lc[lcid].mstrain, macrostrains);
        Mm->add_macro_strains (lcid, macrostrains);
        break;
      case mat_nonlinear_statics:
        addmultv(Mb->lc[lcid+1].mstrain, Mb->lc[lcid].mstrain, Mp->lambda, macrostrains.a, ncomp);
        Mm->add_macro_strains (lcid, macrostrains);
        break;
      case mech_timedependent_prob:
        for(i=0; i<ncomp; i++)
          macrostrains(i) = Mb->dlc[lcid].get_macrostra(Mp->time, i);
        Mm->add_macro_strains (lcid, macrostrains);
        break;
      default:
        print_err("unknown type of problem is required", __FILE__, __LINE__, __func__);
    }
  }
  if (lcid != 0)
  {
    for (ipp=0; ipp<Mm->tnip; ipp++)
    {
      ncomp = Mm->ip[ipp].ncompstr;
      if (Gtm->leso[Mm->elip[ipp]]==1){
        memcpy(Mm->ip[ipp].strain, Mm->ip[ipp].strain+ncomp*lcid, sizeof(Mm->ip[ipp].strain)*ncomp);
      }
    }
  }
}



/**
  Function computes strains at element nodes.
  Nodal strains are averaged. In case of bar and beam elements 
  be carefull if averaging is reasonable.

  @param lcid - load case id
  
  @return The function does not return anything.
 
  Created by JK, 28.11.2006
*/
void compute_nodestrains (long lcid)
{
  long i;
  elemtype te;
  
  //  for open dx
  Mp->detnodstrain = 1;
  
  for (i=0;i<Mt->nn;i++){
    Mt->nodes[i].nullstrain (lcid);
  }
  
  for (i=0;i<Mt->ne;i++){
    if (Gtm->leso[i]==1){
      
      te = Mt->give_elem_type (i);
      
      switch (te){
	
      case bar2d:{
	Bar2d->nod_strains_ip (lcid, i, 0, 0);
	break;
      }
      case barq2d:{
	Barq2d->nod_strains_ip (lcid,i,0,0);
	break;
      }
      case bar3d:{//zakomentovano pro temelin
	//Bar3d->nod_strains_ip (lcid, i, 0, 0);
	break;
      }
      case barq3d:{//zakomentovano pro temelin
	//Barq3d->nod_strains_ip (lcid,i,0,0);
	break;
      }
	
      case beam2d:{
	Beam2d->nodal_displ (lcid,i);
	break;
      }
      case beam3d:{
	Beam3d->nodal_displ (lcid,i);
	break;
      }

      case spring_1:
      case spring_2:
      case spring_3:
      case spring_4:
      case spring_5:
      case spring_6:
	break;

      case planeelementlt:{
	Pelt->nod_strains_ip (lcid,i,0,0);
	break;
      }
      case planeelementqt:{
	Peqt->nod_strains (lcid,i,0,0);
	break;
      }
      case planeelementrotlt:{
	Perlt->nod_strains_ip (lcid,i,0,0);
	break;
      }
	
      case planeelementlq:{
	Pelq->nod_strains_ip (lcid,i,0,0);
	break;
      }
      case planeelementqq:{
	Peqq->nod_strains_ip (lcid,i,0,0);
	break;
      }
      case planeelementrotlq:{
	//Perlq->nod_strains (lcid,i,0,0);
	break;
      }
	
      case planeelementsubqt:{
	Pesqt->nod_strains (lcid,i,0,0);
	break;
      }
	
      case cctel:{
	Cct->nod_strains_ip (lcid,i,0,0);
	break;
      }
      case dktel:{
	//Dkt->nod_strains (lcid,i,0,0);
	break;
      }
      case dstel:{
	//Dst->nod_strains (lcid,i,0,0);
	break;
      }
      case q4plateel:{
	//Q4pl->nod_strains (lcid,i,0,0);
	break;
      }
	
      case subsoilplatetr:{
	break;
      }
	
      case subsoilplateq:{
	//Splq->nod_strains (lcid,i,0,0);
	break;
      }
	
      case axisymmlt:{
	Asymlt->nod_strains_ip (lcid,i);
	break;
      }
      case axisymmqt:{
	Asymqt->nod_strains_ip (lcid,i);
	break;
      }
      case axisymmlq:{
	Asymlq->nod_strains_ip (lcid,i,0,0);
	break;
      }
      case axisymmqq:{
	Asymqq->nod_strains_ip (lcid,i);
	break;
      }
      case axisymmcq:{
	Asymcq->nod_strains_ip (lcid,i);
	break;
      }
      case axisymmlqintface:{
        Asymlqifc->nod_strains_ip(lcid, i, 0, 0);
        break;
      }	
	
      case shelltrelem:{
	Shtr->nod_strains_ip (lcid,i);
	break;
      }
      case shelltrmelem:{
	Shtrm->nod_strains_ip (lcid,i);
	break;
      }
	//case shellqelem:{
	//Shq->nod_strains_ip (lcid,i);
	//break;
	//}
	
      case lineartet:{
	Ltet->nod_strains_ip (lcid,i,0,0);
	break;
      }
      case quadrtet:{
	Qtet->nod_strains_ip (lcid,i,0,0);
	break;
      }
      case linearhex:{
	Lhex->nod_strains_ip (lcid,i,0,0);
	break;
      }
      case quadrhex:{
	Qhex->nod_strains_ip (lcid,i,0,0);
	break;
      }
      case lineartetrot:{
	Ltetrot->nod_strains_ip (lcid,i,0,0);
	break;
      }
      case linearhexrot:{
	Lhexrot->nod_strains_ip (lcid,i,0,0);
	break;
      }
      case hexintface:{
        Hexifc->nod_strains_ip(lcid, i);
        break;
      }
      case particleelem:{
	break;
      }
      case tetralatt:{
	break;
      }
      default:{
	print_err("unknown element type is required", __FILE__, __LINE__, __func__);
      }
      }
    }
  }
  
  //  averaging of nodal strains
  for (i=0;i<Mt->nn;i++){
    //if (Gtm->lnso[i]==1){
    Mt->nodes[i].strain_averageval (lcid);
    //}
  }
}



/**
  Function computes strains at element nodes directly and nodal strains are averaged. 
  In case of bar and beam elements, be carefull if averaging is reasonable. The function 
  is intended for transfer of mechanical quantities in coupled problems.

  @param lcid - load case id
  
  @return The function does not return anything.
 
  Created by Tomas Koudelka, 14.11.2013
*/
void nodestrains_comp (long lcid)
{
  long i;
  elemtype te;
  
  for (i=0;i<Mt->nn;i++)
    Mt->nodes[i].nullstrain(lcid);
  
  for (i=0;i<Mt->ne;i++)
  {
    if (Gtm->leso[i]==1)
    {
      te = Mt->give_elem_type (i);
      
      switch (te)
      {
        case bar2d:
          Bar2d->nod_strains_comp(lcid, i);
          break;
        case barq2d:
          Barq2d->nod_strains_comp(lcid,i);
          break;
        case bar3d:
//          Bar3d->nod_strains_comp(lcid, i, 0, 0);
          break;
        case barq3d:
//          Barq3d->nod_strains_comp(lcid,i,0,0);
          break;
        case spring_1:
        case spring_2:
        case spring_3:
        case spring_4:
        case spring_5:
        case spring_6:
          break;

        case planeelementlt:
//          Pelt->nod_strains_comp(lcid,i,0,0);
          break;
        case planeelementqt:
//          Peqt->nod_strains_comp(lcid,i,0,0);
          break;
	
        case planeelementlq:
          Pelq->nod_strains_comp(lcid,i);
          break;
        case planeelementqq:
          Peqq->nod_strains_comp(lcid,i);
          break;
	
        case axisymmlt:
//          Asymlt->nod_strains_comp(lcid,i);
          break;
        case axisymmlq:
          Asymlq->nod_strains_comp(lcid,i);
          break;
        case axisymmqq:
          Asymqq->nod_strains_comp(lcid,i);
          break;
        case axisymmcq:
          Asymcq->nod_strains_comp(lcid,i);
          break;
        case axisymmlqintface:
          Asymlqifc->nod_strains_comp(lcid, i);
          break;
        case lineartet:
//          Ltet->nod_strains_comp(lcid,i,0,0);
          break;
        case quadrtet:
//          Qtet->nod_strains_comp(lcid,i,0,0);
          break;

        case linearhex:
//          Lhex->nod_strains_comp(lcid,i,0,0);
          break;
     
        case quadrhex:
//          Qhex->nod_strains_comp(lcid,i,0,0);
          break;
        default:
          print_err("unknown element type is required", __FILE__, __LINE__, __func__);
      }
    }
  }
  
  //  averaging of nodal strains
  for (i=0;i<Mt->nn;i++)
  {
    if (Gtm->lnso[i]==1)
      Mt->nodes[i].strain_averageval (lcid);
  }
}



/**
  The function computes strains at strain points.
   
  @param lcid - load case id
   
  @return The function does not return anything.
 
  Created by JK,
*/
void computestrains (long lcid)
{
  
  //  function computes strains at integration points
  compute_ipstrains (lcid);
  
  /*
  long j,nip,ipp,ncomp,k;
  fprintf (Out,"\n\n\n\n kontrola napeti a sil \n\n");
  for (i=0;i<Mt->ne;i++){
    fprintf (Out,"\n prvek %ld  ",i+1);
    ipp=Mt->elements[i].ipp[0][0];
    fprintf (Out,"ipp  %ld  ",ipp);
    nip=Mt->give_tnip (i);
    fprintf (Out,"nip  %ld  ",nip);
    for (j=0;j<nip;j++){
      fprintf (Out,"\n int. bod %ld  ",j+1);
      ncomp=Mm->ip[ipp].ncompstr;
      for (k=0;k<ncomp;k++){
  fprintf (Out,"  %lf",Mm->ip[ipp].strain[k]);
      }
      ipp++;
    }
  }
  */
  

  /*
  //  number of elements
  ne=Mt->ne;
  for (i=0;i<ne;i++){
    if (Gtm->leso[i]==1){
      
      //  type of element
      te = Mt->give_elem_type (i);
      
      switch (te){
      case bar2d:{
  Bar2d->strains (lcid,i,0,0);
  break;
      }
      case bar3d:{
  Bar3d->strains (lcid,i,0,0);
  break;
      }
      case barq2d:{
  Barq2d->strains (lcid,i,0,0);
  break;
      }
      case barq3d:{
  Barq3d->strains (lcid,i,0,0);
  break;
      }
      case beam2d:{
  Beam2d->strains (i,lcid);
  break;
      }
      case beam3d:{
  Beam3d->nodal_displ (i,lcid);
  break;
      }
      case subsoilbeam:{
  Sbeam->strains (lcid, i, 0, 0);
  break;
      }
      case spring_1:
      case spring_2:
      case spring_3:
      case spring_4:
      case spring_5:
      case spring_6:
  {
    Spring->strains (i,lcid);
    break;
  }
  
      case planeelementlt:{
  Pelt->strains (lcid,i,0,0);
  break;
      }
      case planeelementqt:{
  Peqt->strains (lcid,i,0,0);
  break;
      }
      case planeelementrotlt:{
  Perlt->strains (lcid,i,0,0);
  break;
      }
  
      case planeelementlq:{
  Pelq->strains (lcid,i,0,0);
  break;
      }
      case planeelementqq:{
  Peqq->strains (lcid,i,0,0);
  break;
      }
      case planeelementrotlq:{
  Perlq->strains (lcid,i,0,0);
  break;
      }
  
      case planeelementsubqt:{
  Pesqt->strains (lcid,i,0,0);
  break;
      }
  
      case cctel:{
  Cct->strains (lcid,i,0,0);
  break;
      }
      case dktel:{
  Dkt->strains (lcid,i,0,0);
  break;
      }
      case dstel:{
  Dst->strains (lcid,i,0,0);
  break;
      }
      case q4plateel:{
  Q4pl->strains (lcid,i,0,0);
  break;
      }
      case subsoilplateq:{
  Splq->strains (lcid,i,0,0);
  break;
      }
  
  
      case axisymmlt:{
  Asymlt->strains (lcid,i,0,0);
  break;
      }
      case axisymmlq:{
  Asymlq->strains (lcid,i,0,0);
  break;
      }
      case axisymmqq:{
  Asymqq->strains (lcid,i,0,0);
  break;
      }
      case axisymmcq:{
  Asymcq->strains (lcid,i,0,0);
  break;
      }
  
      case shelltrelem:{
  Shtr->strains (lcid,i);
  break;
      }
  
      case lineartet:{
  Ltet->strains (lcid,i,0,0);
  break;
      }
      case quadrtet:{
  Qtet->strains (lcid,i,0,0);
  break;
      }
      case linearhex:{
  Lhex->strains (lcid,i,0,0);
  break;
      }
      case quadrhex:{
  Qhex->strains (lcid,i,0,0);
  break;
      }
  
      case particleelem:{
  break;
      }
  
      default:{
  print_err("unknown element type is required", __FILE__, __LINE__, __func__);
      }
      }
      
    }
  }
  */
}



// ********************
// ********************
// ****  STRESSES  ****
// ********************
// ********************




/**
  The function computes all stresses at integration points.
   
  @param lcid - load case id
   
  @return The function does not return anything.
 
  Created by JK, 28.11.2006
*/
void compute_ipstresses (long lcid)
{
  long i, j, ne, tnip, ncomp, ipp=0L;
  elemtype te;
  
  //  number of elements
  ne=Mt->ne;
  for (i=0;i<ne;i++){
    tnip = Mt->give_tnip(i);
    if (Gtm->leso[i]==1){
      
      //  type of element
      te = Mt->give_elem_type (i);
      
      switch (te){
	
      case bar2d:{
	Bar2d->res_ip_stresses (lcid,i);
	break;
      }
      case bar3d:{
	Bar3d->res_ip_stresses (lcid,i);
	break;
      }
      case barq2d:{
	Barq2d->res_mainip_stresses (lcid,i);
	break;
      }
      case barq3d:{
	Barq3d->res_mainip_stresses (lcid,i);
	break;
      }
	
      case beam2d:{
	Beam2d->nodal_forces (lcid,i);
	break;
      }
	
	
      case beam3d:{
	Beam3d->nodal_forces (lcid,i);
	break;
      }
      case beamg3d:{
	Beam3dg->nodal_forces (lcid,i);
	break;
      }
	
  /*
    case subsoilbeam:{
    Sbeam->strains (lcid, i, 0, 0);
    break;
    }
  */
  
      case spring_1:{
	Spring->stresses(i, lcid);
	break;
      }
      case spring_2:{
	Spring->stresses(i, lcid);
	break;
      }
      case spring_3:{
	Spring->stresses(i, lcid);
	break;
      }
      case spring_4:{
	Spring->stresses(i, lcid);
	break;
      }
      case spring_5:{
	Spring->stresses(i, lcid);
	break;
      }
      case spring_6:{
	Spring->stresses(i, lcid);
	break;
      }
	
      case planeelementlt:{
	Pelt->res_ip_stresses (lcid,i);
	break;
      }
	
      case planeelementqt:{
	Peqt->res_ip_stresses (lcid,i);
	break;
      }
	
      case planeelementrotlt:{
	Perlt->res_ip_stresses (lcid,i);
	break;
      }
	
      case planeelementlq:{
	Pelq->res_ip_stresses (lcid,i);
	break;
      }
	
      case planeelementqq:{
	Peqq->res_ip_stresses (lcid,i);
	break;
      }
	
      case planeelementrotlq:{
	Perlq->res_ip_stresses (lcid,i);
	break;
      }
	/*	  
	  case planeelementsubqt:{
	  Pesqt->res_mainip_strains (lcid,i);
	  break;
	  }
	*/
	
      case planequadinterface:
        Pqifc->res_mainip_stresses(lcid,i);
        break;


      case cctel:{
	Cct->res_ip_stresses (lcid,i);
	break;
      }
      case dktel:{
	Dkt->res_ip_stresses (lcid,i);
	break;
      }
      case dstel:{
	Dst->res_ip_stresses (lcid,i);
	break;
      }
      case q4plateel:{
	Q4pl->res_ip_stresses (lcid,i);
	break;
      }
      case dkqel:{
	Dkqelem->res_moments (lcid,i);
	Dkqelem->res_forces (lcid,i);
	break;
      }
      case subsoilplatetr:{
	break;
      }
      case subsoilplateq:{
	//Splq->res_mainip_strains (lcid,i);
	break;
      }
	
	
      case axisymmlt:{
	Asymlt->res_mainip_stresses (lcid,i);
	break;
      }
      case axisymmqt:{
	Asymqt->res_mainip_stresses (lcid,i);
	break;
      }
      case axisymmlq:{
	Asymlq->res_allip_stresses (lcid,i);
	break;
      }
      case axisymmqq:{
	Asymqq->res_allip_stresses (lcid,i);
	break;
      }
      case axisymmcq:{
	Asymcq->res_allip_stresses (lcid,i);
	break;
      }
      case axisymmlqintface:{
        Asymlqifc->res_mainip_stresses(lcid, i);
        break;
      }	


      case shelltrelem:{
	Shtr->res_ip_stresses (lcid,i);
	break;
      }
      case shelltrmelem:{
	Shtrm->res_ip_stresses (lcid,i);
	break;
      }
      case shellqelem:{
	Shq->res_ip_stresses (lcid,i);
	break;
      }
	

	
      case lineartet:{
	Ltet->res_ip_stresses (lcid,i);
	break;
      }
      case quadrtet:{
	Qtet->res_ip_stresses (lcid,i);
	break;
      }
      case linearhex:{
	Lhex->res_ip_stresses (lcid,i);
	break;
      }
      case quadrhex:{
	Qhex->res_ip_stresses (lcid,i);
	break;
      }
      case lineartetrot:{
	Ltetrot->res_ip_stresses (lcid,i);
	break;
      }
      case linearhexrot:{
	Lhexrot->res_ip_stresses (lcid,i);
	break;
      }
      case hexintface:{
        Hexifc->res_mainip_stresses(lcid, i);
        break;
      }
	
      case particleelem:{
	break;
      }
	
      default:{
	print_err("unknown element type is required", __FILE__, __LINE__, __func__);
      }
      }
      // stresses are stored at the beginning of ip stress array
      // they have to be moved to the proper components of the stress array according to lcid
      if (lcid != 0)
      {
        for (j=0; j<tnip; j++)
        {
          ncomp = Mm->ip[ipp+j].ncompstr;
          memcpy(Mm->ip[ipp+j].stress+ncomp*lcid, Mm->ip[ipp+j].stress, sizeof(Mm->ip[ipp+j].stress)*ncomp);
        }
      }
    }
    ipp += tnip;
  }

  //find_extreme_ipstresses (lcid);
}

/**
   The function searches the minimum and maximum stresses in integration points.
   
   @param lcid - load case id
   
   @return The function does not return anything.
   
   Created by JK, 16. 11. 2012
*/
void find_extreme_ipstresses (long lcid)
{
  long i,ne;
  long nebar2d,nebar3d,nelinhex;
  vector minbar2d,minbar3d,minlinhex;
  vector maxbar2d,maxbar3d,maxlinhex;
  
  if (Bar2d != NULL){
    allocv (Bar2d->tncomp,minbar2d);
    fillv (1.0e80,minbar2d);
    vector maxbar2d(Bar2d->tncomp);
    fillv (-1.0e80,maxbar2d);
  }
  if (Bar3d != NULL){
    allocv  (Bar3d->tncomp,minbar3d);
    fillv (1.0e80,minbar3d);
    allocv (Bar3d->tncomp,maxbar3d);
    fillv (-1.0e80,maxbar3d);
  }
  if (Lhex != NULL){
    allocv (Lhex->tncomp,minlinhex);
    fillv (1.0e80,minlinhex);
    allocv (Lhex->tncomp,maxlinhex);
    fillv (-1.0e80,maxlinhex);
  }
  
  elemtype te;
  
  nebar2d=0;
  nebar3d=0;
  nelinhex=0;
  
  
  //  the number of elements
  ne=Mt->ne;
  
  //  loop over the number of elements
  for (i=0;i<ne;i++){
    if (Gtm->leso[i]==1){
      
      //  type of element
      te = Mt->give_elem_type (i);
      
      switch (te){
	
	/*
      case bar2d:{
	Bar2d->res_ip_stresses (lcid,i);
	break;
      }
	*/

      case bar3d:{
	Bar3d->find_extreme_stresses (minbar3d,maxbar3d,lcid,i);
	nebar3d++;
	break;
      }
	/*
      case barq2d:{
	Barq2d->res_mainip_stresses (lcid,i);
	break;
      }
      case barq3d:{
	Barq3d->res_mainip_stresses (lcid,i);
	break;
      }
	
      case beam2d:{
	Beam2d->nodal_forces (lcid,i);
	break;
      }
	
	
      case beam3d:{
	Beam3d->nodal_forces (lcid,i);
	break;
      }
      case beamg3d:{
	Beam3dg->nodal_forces (lcid,i);
	break;
      }
	
      
    case subsoilbeam:{
    Sbeam->strains (lcid, i, 0, 0);
    break;
    }

  
      case spring_1:{
	Spring->stresses(i, lcid);
	break;
      }
      case spring_2:{
	Spring->stresses(i, lcid);
	break;
      }
      case spring_3:{
	Spring->stresses(i, lcid);
	break;
      }
      case spring_4:{
	Spring->stresses(i, lcid);
	break;
      }
      case spring_5:{
	Spring->stresses(i, lcid);
	break;
      }
      case spring_6:{
	Spring->stresses(i, lcid);
	break;
      }
	
      case planeelementlt:{
	Pelt->res_ip_stresses (lcid,i);
	break;
      }

    case planeelementqt:{
    Peqt->res_mainip_strains (lcid,i);
    break;
    }
    


      case planeelementrotlt:{
	Perlt->res_ip_stresses (lcid,i);
	break;
      }
	
      case planeelementlq:{
	Pelq->res_ip_stresses (lcid,i);
	break;
      }
	
      case planeelementqq:{
	Peqq->res_ip_stresses (lcid,i);
	break;
      }

	  case planeelementrotlq:{
	  Perlq->res_mainip_strains (lcid,i);
	  break;
	  }
	  
	  case planeelementsubqt:{
	  Pesqt->res_mainip_strains (lcid,i);
	  break;
	  }

	
      case planequadinterface:
        Pqifc->res_mainip_stresses(lcid,i);
        break;
      case cctel:{
	Cct->res_ip_stresses (lcid,i);
	break;
      }
	
      case dktel:{
	Dkt->res_ip_stresses (lcid,i);
	break;
      }
      case dstel:{
	Dst->res_ip_stresses (lcid,i);
	break;
      }
      case q4plateel:{
	Q4pl->res_ip_stresses (lcid,i);
	break;
      }
	
      case subsoilplatetr:{
	break;
      }
	
      case subsoilplateq:{
	//Splq->res_mainip_strains (lcid,i);
	break;
      }
	
	
      case shelltrelem:{
	Shtr->res_ip_stresses (lcid,i);
	break;
      }
      case shellqelem:{
	Shq->res_ip_stresses (lcid,i);
	break;
      }
	
      case axisymmlt:{
	Asymlt->res_mainip_stresses (lcid,i);
	break;
      }
      case axisymmlq:{
	Asymlq->res_allip_stresses (lcid,i);
	break;
      }
      case axisymmqq:{
	Asymqq->res_allip_stresses (lcid,i);
	break;
      }
      case axisymmcq:{
	Asymcq->res_allip_stresses (lcid,i);
	break;
      }
	
      case lineartet:{
	Ltet->res_ip_stresses (lcid,i);
	break;
      }
      case quadrtet:{
	Qtet->res_ip_stresses (lcid,i);
	break;
      }
	*/
      case linearhex:{
	Lhex->find_extreme_stresses (minlinhex,maxlinhex,lcid,i);
	nelinhex++;
	break;
      }
	/*
      case quadrhex:{
	Qhex->res_ip_stresses (lcid,i);
	break;
      }
      case lineartetrot:{
	Ltetrot->res_ip_stresses (lcid,i);
	break;
      }
      case linearhexrot:{
	Lhexrot->res_ip_stresses (lcid,i);
	break;
      }
	
      case particleelem:{
	break;
      }
	*/
      default:{
	print_err("unknown element type is required", __FILE__, __LINE__, __func__);
      }
      }
    }
  }
  
  fprintf (stdout,"\n\n\n");
  if (Bar3d != NULL)
    fprintf (stdout,"\n linbar3d  %6ld        % le % le",nebar3d,minbar3d[0],maxbar3d[0]);
  if (Lhex != NULL){
    fprintf (stdout,"\n linhex    %6ld        % le % le % le % le % le % le\n",nelinhex,minlinhex[0],minlinhex[1],minlinhex[2],minlinhex[3],minlinhex[4],minlinhex[5]);
    fprintf (stdout,"                         % le % le % le % le % le % le",maxlinhex[0],maxlinhex[1],maxlinhex[2],maxlinhex[3],maxlinhex[4],maxlinhex[5]);
  }
  
  fprintf (stdout,"\n\n\n");
  
  destrv (minbar2d);
  destrv (minbar3d);
  destrv (minlinhex);

  destrv (maxbar2d);
  destrv (maxbar3d);
  destrv (maxlinhex);
}


/**
  The function computes stresses at element nodes.
  Nodal stresses are averaged. In case of bar and 
  beam elements, be carefull if averaging is reasonable.
  
  @param lcid - load case id
   
  @return The function does not return anything.
 
  Created by JK, 28.11.2006
*/
void compute_nodestresses (long lcid)
{
  long i;
  elemtype te;
  
  for (i=0;i<Mt->nn;i++){
    Mt->nodes[i].nullstress (lcid);
  }
  
  for (i=0;i<Mt->ne;i++){
    if (Gtm->leso[i]==1){
      
      te = Mt->give_elem_type (i);
      
      switch (te){
	
      case bar2d:{
	Bar2d->nod_stresses_ip (lcid, i, 0, 0);
	break;
      }
      case barq2d:{
	Barq2d->nod_stresses_ip (lcid,i,0,0);
	break;
      }
      case bar3d:{//zakomentovano pro temelin
	//Bar3d->nod_stresses_ip (lcid, i, 0, 0);
	break;
      }
      case barq3d:{//zakomentovano pro temelin
	//Barq3d->nod_stresses_ip (lcid,i,0,0);
	break;
      }
	
      case beam2d:{
	Beam2d->nodal_forces (lcid,i);
	break;
      }
      case beam3d:{
	Beam3d->nodal_forces (lcid,i);
	break;
      }
      case beamg3d:{
	Beam3dg->nodal_forces (lcid,i);
	break;
      }
	
	
	
      case spring_1:
      case spring_2:
      case spring_3:
      case spring_4:
      case spring_5:
      case spring_6:
	break;
      case planeelementlt:{
	Pelt->nod_stresses_ip (lcid,i,0,0);
	break;
      }
      case planeelementqt:{
	Peqt->nod_stresses (lcid,i,0,0);
	break;
      }
      case planeelementrotlt:{
	Perlt->nod_stresses_ip (lcid,i,0,0);
	break;
      }
	
      case planeelementlq:{
	Pelq->nod_stresses_ip (lcid,i,0,0);
	break;
      }
      case planeelementqq:{
	Peqq->nod_stresses_ip (lcid,i,0,0);
	break;
      }
      case planeelementrotlq:{
	//Perlq->nod_stresses (lcid,i,0,0);
	break;
      }
	
      case planeelementsubqt:{
	Pesqt->nod_stresses (lcid,i,0,0);
	break;
      }
	
      case cctel:{
	Cct->nod_stresses_ip (lcid,i,0,0);
	break;
      }
      case dktel:{
	//Dkt->nod_stresses (lcid,i,0,0);
	break;
      }
      case dstel:{
	//Dst->nod_stresses (lcid,i,0,0);
	break;
      }
      case q4plateel:{
	//Q4pl->nod_stresses (lcid,i,0,0);
	break;
      }
	
      case subsoilplatetr:{
	break;
      }
      case subsoilplateq:{
	break;
      }
	
      case axisymmlt:{
	Asymlt->nod_stresses_ip (lcid,i);
	break;
      }
      case axisymmqt:{
	Asymqt->nod_stresses_ip (lcid,i);
	break;
      }
      case axisymmlq:{
	Asymlq->nod_stresses_ip (lcid,i,0,0);
	break;
      }
      case axisymmqq:{
	Asymqq->nod_stresses_ip (lcid,i);
	break;
      }
      case axisymmcq:{
	Asymcq->nod_stresses_ip (lcid,i);
	break;
      }
      case axisymmlqintface:{
        Asymlqifc->nod_stresses_ip(lcid, i, 0, 0);
        break;
      }	
	
      case shelltrelem:{
	Shtr->nod_stresses_ip (lcid,i);
	break;
      }
      case shelltrmelem:{
	Shtrm->nod_stresses_ip (lcid,i);
	break;
      }
	//case shellqelem:{
	//Shq->nod_stresses_ip (lcid,i);
	//break;
	//}
	
      case lineartet:{
	Ltet->nod_stresses_ip (lcid,i,0,0);
	break;
      }
      case quadrtet:{
	Qtet->nod_stresses_ip (lcid,i,0,0);
	break;
      }
      case linearhex:{
	Lhex->nod_stresses_ip (lcid,i,0,0);
	break;
      }
      case quadrhex:{
	Qhex->nod_stresses_ip (lcid,i,0,0);
	break;
      }
      case lineartetrot:{
	Ltetrot->nod_stresses_ip (lcid,i,0,0);
	break;
      }
      case linearhexrot:{
	Lhexrot->nod_stresses_ip (lcid,i,0,0);
	break;
      }
      case hexintface:{
        Hexifc->nod_stresses_ip(lcid, i);
        break;
      }
      case particleelem:{
	break;
      }
      case tetralatt:{
	break;
      }
	
      default:{
	print_err("unknown element type is required", __FILE__, __LINE__, __func__);
      }
      }
      
    }
  }
  
  for (i=0;i<Mt->nn;i++){
    //Mt->nodes[i].strain_averageval (lcid);
    Mt->nodes[i].stress_averageval (lcid);
  }
  
}



/**
  The function computes stresses at stress points.
   
  @param lcid - load case id
   
  @return The function does not return anything.
 
  Created by JK,
*/
void computestresses (long lcid)
{
  long i,ne,te;
  
  //  all stress components at all integration points
  compute_ipstresses (lcid);
  
  //  number of elements
  ne=Mt->ne;
  for (i=0;i<ne;i++){
    if (Gtm->leso[i]==1){
      
      te = Mt->give_elem_type (i);
      
      switch (te){
	
      case bar2d:{
	Bar2d->stresses (lcid,i,0,0);
	break;
      }
      case bar3d:{
	Bar3d->stresses (lcid,i,0,0);
	break;
      }
      case barq2d:{
	Barq2d->stresses (lcid,i,0,0);
	break;
      }
      case barq3d:{
	Barq3d->stresses (lcid,i,0,0);
	break;
      }
	
      case beam2d:{
	Beam2d->nodal_forces (lcid,i);
	break;
      }
	
      case beam3d:{
	Beam3d->nodal_forces (lcid,i);
	break;
      }
	
      case spring_1:
      case spring_2:
      case spring_3:
      case spring_4:
      case spring_5:
      case spring_6:
	Spring->stresses(i, lcid);
	break;
      case planeelementlt:{
	Pelt->stresses (lcid,i,0,0);
	break;
      }
      case planeelementqt:{
	Peqt->stresses (lcid,i,0,0);
	break;
      }
      case planeelementrotlt:{
	Perlt->stresses (lcid,i,0,0);
	break;
      }
	
      case planeelementlq:{
	Pelq->stresses (lcid,i,0,0);
	break;
      }
      case planeelementqq:{
	Peqq->stresses (lcid,i,0,0);
	break;
      }
      case planeelementrotlq:{
	Perlq->stresses (lcid,i,0,0);
	break;
      }
	
      case planeelementsubqt:{
	Pesqt->stresses (lcid,i,0,0);
	break;
      }
	
      case cctel:{
	Cct->stresses (lcid,i,0,0);
	break;
      }
      case dktel:{
	Dkt->stresses (lcid,i,0,0);
	break;
      }
	
      case dstel:{
	Dst->stresses (lcid,i,0,0);
	break;
      }
	
      case q4plateel:{
	Q4pl->stresses (lcid,i,0,0);
	break;
      }
	
      case axisymmlt:{
	Asymlt->stresses (lcid,i,0,0);
	break;
      }
      case axisymmlq:{
	Asymlq->stresses (lcid,i,0,0);
	break;
      }
      case axisymmqq:{
	Asymqq->stresses (lcid,i,0,0);
	break;
      }
      case axisymmcq:{
	Asymcq->stresses (lcid,i,0,0);
	break;
      }
	
      case shelltrelem:{
	Shtr->stresses (lcid,i);
	break;
      }
      case shelltrmelem:{
	Shtrm->stresses (lcid,i);
	break;
      }
	//case shellqelem:{
	//Shq->stresses (lcid,i);
	//break;
	//}
	
      case lineartet:{
	Ltet->stresses (lcid,i,0,0);
	break;
      }
      case quadrtet:{
	Qtet->stresses (lcid,i,0,0);
	break;
      }
      case linearhex:{
	Lhex->stresses (lcid,i,0,0);
	break;
      }
      case quadrhex:{
	Qhex->stresses (lcid,i,0,0);
	break;
      }
      case lineartetrot:{
	Ltetrot->stresses (lcid,i,0,0);
	break;
      }
      case linearhexrot:{
	Lhexrot->stresses (lcid,i,0,0);
	break;
      }
      case particleelem:{
	break;
      }
      default:{
	print_err("unknown element type is required", __FILE__, __LINE__, __func__);
      }
      }
    }
  }
}



/**
  The function computes other values in nodes.
   
  @param lcid - load case id

  @return The function does not return anything.

  Created by JK, 
*/
void compute_nodeothers ()
{
  long i;
  elemtype te;
  
  for (i=0;i<Mt->nn;i++){
    Mt->nodes[i].nullother ();
  }
  
  for (i=0;i<Mt->ne;i++){
    if (Gtm->leso[i]==1){
      
      te = Mt->give_elem_type (i);
      
      switch (te){
	
      case bar2d:{
	Bar2d->nod_other_ip (i, 0, 0);
	break;
      }
      case bar3d:{//zakomentovano pro temelin
	//Bar3d->nod_eqother_ip (i, 0, 0);
	break;
      }
      case barq2d:{
	Barq2d->nod_other_ip (i, 0, 0);
	break;
      }
      case barq3d:{//zakomentovano pro temelin
	//Barq3d->nod_eqother_ip (i, 0, 0);
	break;
      }
	/*    case beam2d:{
	      Beam2d->stresses (i);
	      break;
	      }
	      case beam3d:{
	      Beam3d->stresses (i);
	      break;
	      }
	      case beam2d:{
	      break;
	      }
	*/
	
      case spring_1:
      case spring_2:
      case spring_3:
      case spring_4:
      case spring_5:
      case spring_6:
	break;
	
      case planeelementlt:{
	Pelt->nod_others (i,0,0);
	break;
      }
	
	/*
	  case planeelementrotlt:{
	  Perlt->nod_others (lcid,i,0,0);
	  break;
	  }
	*/
	
      case planeelementlq:{
	Pelq->nod_other_ip (i,0,0);
	break;
      }
	
      case planeelementqq:{
	Peqq->nod_other_ip (i,0,0);
	break;
      }
	/*
	  case planeelementrotlq:{
	  Perlq->nod_stresses (lcid,i,0,0);
	  break;
	  }
	  
	  case planeelementsubqt:{
	  Pesqt->nod_stresses (lcid,i,0,0);
	  break;
	  }
	  
	  case cctel:{
	  Cct->nod_stresses (lcid,i,0,0);
	  break;
	  }
	  case dktel:{
	  Dkt->nod_stresses (lcid,i,0,0);
	  break;
	  }
	  case dstel:{
	  Dst->nod_stresses (lcid,i,0,0);
	  break;
	  }
	  case q4plateel:{
	  Q4pl->nod_stresses (lcid,i,0,0);
	  break;
	  }*/
	
      case axisymmlt:{
	Asymlt->nod_other_ip (i);
	break;
      }
      case axisymmqt:{
	Asymqt->nod_other_ip (i);
	break;
      }
      case axisymmlq:{
	Asymlq->nod_other_ip (i);
	break;
      }
      case axisymmqq:{
	Asymqq->nod_other_ip (i);
	break;
      }
      case axisymmcq:{
	Asymcq->nod_other_ip (i);
	break;
      }
      case axisymmlqintface:{
        Asymlqifc->nod_other_ip(i);
        break;
      }
        
	
	/*
	  case shelltrelem:{
	  Shtr->nod_stresses (lcid,i);
	  break;
	  }
	*/
	
      case lineartet:{
	Ltet->nod_other_ip (i);
	break;
      }
      case quadrtet:{
	Qtet->nod_other_ip (i,0,0);
	break;
      }
      case linearhex:{
	Lhex->nod_other_ip (i,0,0);
	break;
      }
      case quadrhex:{
	Qhex->nod_other_ip (i,0,0);
	break;
      }
	
      default:{
	print_err("unknown element type is required", __FILE__, __LINE__, __func__);
      }
      }
    }
  }
  
  for (i=0;i<Mt->nn;i++){
    Mt->nodes[i].other_averageval ();
  }
}



/**
  The function interpolates nodal values to integration points on one element.

  @param eid - element id
  @param nodval - %vector of nodal values
  @param ipval - %vector of values at integration points (output)
   
  @return The function returns interpolated values in the parameter ipval.
 
  Created by JK, 29.11.2006
*/
void elem_intpointval (long eid,vector &nodval,vector &ipval)
{
  elemtype te;

  te = Mt->give_elem_type (eid);

  switch (te){
    case bar2d:{
      Bar2d->intpointval (eid,nodval,ipval);
      break;
    }
    case bar3d:{
      Bar3d->intpointval (eid,nodval,ipval);
      break;
    }
    case barq2d:{
      Barq2d->intpointval (eid,nodval,ipval);
      break;
    }
    case barq3d:{
      Barq3d->intpointval (eid,nodval,ipval);
      break;
    }
      
    case spring_1:{
      Spring->intpointval (eid,nodval,ipval);
      break;
    }
    case spring_2:{
      Spring->intpointval (eid,nodval,ipval);
      break;
    }
    case spring_3:{
      Spring->intpointval (eid,nodval,ipval);
      break;
    }
    case spring_4:{
      Spring->intpointval (eid,nodval,ipval);
      break;
    }
    case spring_5:{
      Spring->intpointval (eid,nodval,ipval);
      break;
    }
    case spring_6:{
      Spring->intpointval (eid,nodval,ipval);
      break;
    }
    
    case planeelementlt:{
      Pelt->intpointval (eid,nodval,ipval);
      break;
    }

    case planeelementlq:{
      Pelq->intpointval (eid,nodval,ipval);
      break;
    }
    case planeelementqq:{
      Peqq->intpointval (eid,nodval,ipval);
      break;
    }
    case planequadinterface:{
      Pqifc->intpointval (eid,nodval,ipval);
      break;
    }
    
    case axisymmlq:{
      Asymlq->intpointval (eid,nodval,ipval);
      break;
    }
    case axisymmqq:{
      Asymqq->intpointval (eid,nodval,ipval);
      break;
    }
    case axisymmcq:{
      Asymcq->intpointval (eid,nodval,ipval);
      break;
    }
    case axisymmlqintface:{
      Asymlqifc->intpointval(eid, nodval, ipval);
      break;
    }
        
    
    case lineartet:{
      Ltet->intpointval (eid,nodval,ipval);
      break;
    }
    case quadrtet:{
      Qtet->intpointval (eid,nodval,ipval);
      break;
    }
    case linearhex:{
      Lhex->intpointval (eid,nodval,ipval);
      break;
    }
    case quadrhex:{
      Qhex->intpointval (eid,nodval,ipval);
      break;
    }
    case lineartetrot:{
      Ltetrot->intpointval (eid,nodval,ipval);
      break;
    }
    case linearhexrot:{
      Lhexrot->intpointval (eid,nodval,ipval);
      break;
    }
    case hexintface:{
      Hexifc->intpointval(eid, nodval, ipval);
      break;
    }

    default:{
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
    }
  }
  
}



/**
  The function interpolates nodal values into integration points.

  @param gv - array containing values at all nodes of the mesh
  @param nmq - type of non-mechanical quantity
  @param scale - scale coefficient
   
  @return The function does not return anything.
 
  Created by JK, 21.6.2004, extended by TKr 31/08/2012
*/
void intpointval (double *gv,nonmechquant nmq,double scale)
{
  long i,j,nne,nip,ipid;
  ivector nodes;
  vector nodval,ipval;
  
  for (i=0;i<Mt->ne;i++){
    if (Gtm->leso[i]==1){
      
      nne = Mt->give_nne (i);
      nip = Mt->give_tnip (i);
      
      reallocv (RSTCKIVEC(nne,nodes));
      reallocv (RSTCKVEC(nne,nodval));
      reallocv (RSTCKVEC(nip,ipval));
      //  nodes on element
      Mt->give_elemnodes (i,nodes);
      //  nodal values on element
      for (j=0;j<nne;j++){
	nodval[j]=gv[nodes[j]]*scale;
      }
      //  interpolation of nodal values to integration points
      elem_intpointval (i,nodval,ipval);
      
      //  number of the first integration point
      ipid=Mt->elements[i].ipp[0][0];
      
      for (j=0;j<nip;j++)
      {
        Mm->storenonmechq(nmq, ipid, ipval[j]);
        ipid++;
      }
    }
  }
}



/**
  Function interpolates nodal values into integration points.
  Order of the used approximation functions is one degree less then on the given element.
   
  @param gv - array containing values at all nodes of the mesh
  @param nmq - type of non-mechanical quantity
   
  @return The function does not return anything.
 
  Created by 21.6.2004, JK, extended by TKr 31/08/2012
*/
void intpointval2 (double *gv,nonmechquant nmq)
{
  long i,j,nne,nip,ipp;
  ivector nodes;
  vector nodval,ipval;
  
  for (i=0;i<Mt->ne;i++){
    if (Gtm->leso[i]==1){
      nne = Mt->give_nne (i);
      nip = Mt->give_tnip (i);
      
      reallocv (nne,nodes);
      reallocv (nne,nodval);
      reallocv (nip,ipval);
      
      Mt->give_elemnodes (i,nodes);
      
      for (j=0;j<nne;j++){
	nodval[j]=gv[nodes[j]];
      }
      
      elem_intpointval2(i, nodval, ipval);
      
      //  number of the first integration point
      ipp=Mt->elements[i].ipp[0][0];
      
      for (j=0;j<nip;j++)
      {
        Mm->storenonmechq(nmq, ipp, ipval[j]);
        ipp++;
      }
    }
  }
}



/**
  Function interpolates nodal values into integration points.
  Order of the used approximation functions is one degree less then on the given element.
   
  @param eid - element id
  @param nodval - %vector containing values at corner nodes
  @param ipval - %vector of integration point values (output)
   
  @return The function returns resulting integration point values in the %vector ipval.
 
  Created by Tomas Koudelka, 10.12.2013
*/
void elem_intpointval2 (long eid, vector &nodval, vector &ipval)
{
  elemtype te;
  
  te = Mt->give_elem_type (eid);
    
  switch (te)
  {
    case bar2d:
      Bar2d->intpointval (eid,nodval,ipval);
      break;

    case barq2d:
      Barq2d->intpointval2 (eid,nodval,ipval);
      break;

    case barq3d:
      Barq3d->intpointval2 (eid,nodval,ipval);
      break;

    case spring_1:
      Spring->intpointval (eid,nodval,ipval);
      break;

    case spring_2:
      Spring->intpointval (eid,nodval,ipval);
      break;

    case spring_3:
      Spring->intpointval (eid,nodval,ipval);
      break;

    case spring_4:
      Spring->intpointval (eid,nodval,ipval);
      break;

    case spring_5:
      Spring->intpointval (eid,nodval,ipval);
      break;

    case spring_6:
      Spring->intpointval (eid,nodval,ipval);
      break;
	
    case planeelementqq:
      Peqq->intpointval2 (eid,nodval,ipval);
      break;

    case axisymmqq:
      Asymqq->intpointval2 (eid,nodval,ipval);
      break;
    case axisymmcq:
      Asymcq->intpointval2 (eid,nodval,ipval);
      break;
    case quadrtet:
      Qtet->intpointval2 (eid,nodval,ipval);
      break;
      
    case quadrhex:
      Qhex->intpointval2 (eid,nodval,ipval);
      break;
      
    default:
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function interpolates arbitrary quantity defined by its
  nodal values to the required inner point given by natural coordinates.
 
  @param eid    - element id
  @param nodval - %vector containing nodal values
  @param coord  - %vector of natural coordinates of the required point
   
  @return The function returns required interpolated value.
 
  Created by  JK, 23.2.2002
  Modified by TKo, 15.6.2018
*/
double interpolate (long eid, vector &nodval, vector &coord)
{
  double val;
  elemtype te;
  vector areacoord(3), volcoord(4);
  
  te = Mt->give_elem_type (eid);

  switch (te){
    case axisymmlt:{
      areacoord[0]=coord[0];  areacoord[1]=coord[1];  areacoord[2]=1.0-areacoord[0]-areacoord[1];
      val=Asymlt->approx (areacoord,nodval);
      break;
    }
    case axisymmqt:{
      areacoord[0]=coord[0];  areacoord[1]=coord[1];  areacoord[2]=1.0-areacoord[0]-areacoord[1];
      val=Asymqt->approx (areacoord,nodval);
      break;
    }
    case axisymmlq:{
      val=Asymlq->approx (coord[0],coord[1],nodval);
      break;
    }
    case axisymmqq:{
      val=Asymqq->approx (coord[0],coord[1],nodval);
      break;
    }
    case axisymmcq:{
      val=Asymcq->approx (coord[0],coord[1],nodval);
      break;
    }
    case axisymmlqintface:{
      Asymlqifc->approx(coord[0], nodval);
      break;
    }
        
    case planeelementlt:{
      areacoord[0]=coord[0];  areacoord[1]=coord[1];  areacoord[2]=1.0-areacoord[0]-areacoord[1];
      val=Pelt->approx (areacoord,nodval);
      break;
    }
    case planeelementqt:{
      areacoord[0]=coord[0];  areacoord[1]=coord[1];  areacoord[2]=1.0-areacoord[0]-areacoord[1];
      val=Peqt->approx (areacoord[0],areacoord[1],nodval);
      break;
    }
    case planeelementrotlt:{
      areacoord[0]=coord[0];  areacoord[1]=coord[1];  areacoord[2]=1.0-areacoord[0]-areacoord[1];
      val=Perlt->approx (areacoord,nodval);
      break;
    }
    case planeelementlq:{
      val=Pelq->approx (coord[0],coord[1],nodval);
      break;
    }
    case planeelementqq:{
      val=Peqq->approx (coord[0],coord[1],nodval);
      break;
    }
    case planeelementrotlq:{
      val=Perlq->approx (coord[0],coord[1],nodval);
      break;
    }
    case lineartet:{
      volcoord[0]=coord[0];  volcoord[1]=coord[1];  volcoord[2]=coord[2];  volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
      val=Ltet->approx (volcoord,nodval);
      break;
    }
    case quadrtet:{
      val=Qtet->approx (coord[0],coord[1],coord[2],nodval);
      break;
    }
    case linearhex:{
      val=Lhex->approx (coord[0],coord[1],coord[2],nodval);
      break;
    }
    case quadrhex:{
      val=Qhex->approx (coord[0],coord[1],coord[2],nodval);
      break;
    }
    default:{
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
      abort();
    }
  }
  
  return val;
}



/**
  The function computes global coordinates of center of gravity for the given element eid.
 
  @param eid - element id (input)
  @param coord - %vector containing resulting coordinates (output)
   
  @return The function returns required coordinates of centroid.
 
  Created by  TKo, 1.12.2016
*/
void centroid(long eid, vector &coord)
{
  elemtype te = Mt->give_elem_type (eid);
  long nne = Mt->give_nne (eid);
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne));
  vector areacoord(ASTCKVEC(3)),volcoord(ASTCKVEC(4));

  Mt->give_node_coord3d(x, y, z, eid);

  switch (te)
  {
    case bar2d:
      coord(0)=Bar2d->approx(0.0, x);
      coord(1)=Bar2d->approx(0.0, y);
      coord(2)=Bar2d->approx(0.0, z);
      break;
    case beam2d:
      coord(0)=Beam2d->approx(0.0, x);
      coord(1)=Beam2d->approx(0.0, y);
      coord(2)=Beam2d->approx(0.0, z);
      break;
    case bar3d:
      coord(0)=Bar3d->approx(0.0, x);
      coord(1)=Bar3d->approx(0.0, y);
      coord(2)=Bar3d->approx(0.0, z);
      break;
    case beam3d:
      coord(0)=Beam3d->approx(0.0, x);
      coord(1)=Beam3d->approx(0.0, y);
      coord(2)=Beam3d->approx(0.0, z);
      break;
    case barq2d:
      coord(0)=Barq2d->approx(0.0, x);
      coord(1)=Barq2d->approx(0.0, y);
      coord(2)=Barq2d->approx(0.0, z);
      break;
    case barq3d:
      coord(0)=Barq3d->approx(0.0, x);
      coord(1)=Barq3d->approx(0.0, y);
      coord(2)=Barq3d->approx(0.0, z);
      break;
    case axisymmlt:
      areacoord[0]=0.333333333333333;  areacoord[1]=0.333333333333333;  areacoord[2]=1.0-areacoord[0]-areacoord[1];
      coord(0)=Asymlt->approx(areacoord, x);
      coord(1)=Asymlt->approx(areacoord, y);
      coord(2)=Asymlt->approx(areacoord, z);
      break;
    case axisymmqt:
      areacoord[0]=0.333333333333333;  areacoord[1]=0.333333333333333;  areacoord[2]=1.0-areacoord[0]-areacoord[1];
      coord(0)=Asymqt->approx(areacoord, x);
      coord(1)=Asymqt->approx(areacoord, y);
      coord(2)=Asymqt->approx(areacoord, z);
      break;
    case axisymmlqintface:{
      coord(0)=Asymlqifc->approx(0.0, x);
      coord(1)=Asymlqifc->approx(0.0, y);
      coord(2)=Asymlqifc->approx(0.0, z);
      break;
    }
        
    case planeelementlt:
      areacoord[0]=0.333333333333333;  areacoord[1]=0.333333333333333;  areacoord[2]=1.0-areacoord[0]-areacoord[1];
      coord(0)=Pelt->approx(areacoord, x);
      coord(1)=Pelt->approx(areacoord, y);
      coord(2)=Pelt->approx(areacoord, z);
      break;
    case planeelementqt:
      areacoord[0]=0.333333333333333;  areacoord[1]=0.333333333333333;  areacoord[2]=1.0-areacoord[0]-areacoord[1];
      coord(0)=Peqt->approx(areacoord[0], areacoord[1], x);
      coord(1)=Peqt->approx(areacoord[0], areacoord[1], y);
      coord(2)=Peqt->approx(areacoord[0], areacoord[1], z);
      break;
    case planeelementrotlt:
      areacoord[0]=0.333333333333333;  areacoord[1]=0.333333333333333;  areacoord[2]=1.0-areacoord[0]-areacoord[1];
      coord(0)=Perlt->approx(areacoord, x);
      coord(1)=Perlt->approx(areacoord, y);
      coord(2)=Perlt->approx(areacoord, z);
      break;
    case axisymmlq:
      coord(0)=Asymlq->approx(0.0, 0.0, x);
      coord(1)=Asymlq->approx(0.0, 0.0, y);
      coord(2)=Asymlq->approx(0.0, 0.0, z);
      break;
    case planeelementlq:
      coord(0)=Pelq->approx(0.0, 0.0, x);
      coord(1)=Pelq->approx(0.0, 0.0, y);
      coord(2)=Pelq->approx(0.0, 0.0, z);
      break;
    case axisymmqq:
      coord(0)=Asymqq->approx(0.0, 0.0, x);
      coord(1)=Asymqq->approx(0.0, 0.0, y);
      coord(2)=Asymqq->approx(0.0, 0.0, z);
      break;
    case planeelementqq:
      coord(0)=Peqq->approx(0.0, 0.0, x);
      coord(1)=Peqq->approx(0.0, 0.0, y);
      coord(2)=Peqq->approx(0.0, 0.0, z);
      break;
    case axisymmcq:
      coord(0)=Asymcq->approx(0.0, 0.0, x);
      coord(1)=Asymcq->approx(0.0, 0.0, y);
      coord(2)=Asymcq->approx(0.0, 0.0, z);
      break;
    case planeelementrotlq:
      coord(0)=Perlq->approx(0.0, 0.0, x);
      coord(1)=Perlq->approx(0.0, 0.0, y);
      coord(2)=Perlq->approx(0.0, 0.0, z);
      break;
    case lineartet:
      volcoord[0]=0.25;  volcoord[1]=0.25;  volcoord[2]=0.25;  volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
      coord(0)=Ltet->approx(volcoord, x);
      coord(1)=Ltet->approx(volcoord, y);
      coord(2)=Ltet->approx(volcoord, z);
      break;
    case quadrtet:
      coord(0)=Qtet->approx(0.25,0.25,0.25, x);
      coord(1)=Qtet->approx(0.25,0.25,0.25, y);
      coord(2)=Qtet->approx(0.25,0.25,0.25, z);
      break;
    case linearhex:
      coord(0)=Lhex->approx(0.0, 0.0, 0.0, x);
      coord(1)=Lhex->approx(0.0, 0.0, 0.0, y);
      coord(2)=Lhex->approx(0.0, 0.0, 0.0, z);
      break;
    case quadrhex:
      coord(0)=Qhex->approx(0.0, 0.0, 0.0, x);
      coord(1)=Qhex->approx(0.0, 0.0, 0.0, y);
      coord(2)=Qhex->approx(0.0, 0.0, 0.0, z);
      break;
    default:
      print_err("unknown element type is required (eid=%ld)", __FILE__, __LINE__, __func__, eid);
  }
  return;
}



/**
  The function computes volume appropriate to integration points.
   
  @return The function does not return anything.

  Created by JK, 2.3.2004
*/
void ipvolume ()
{
  long i;
  elemtype te;
  
  if (Mm->ipv != NULL)
    delete [] Mm->ipv;
  
  Mm->ipv = new double [Mm->tnip];
  
  for (i=0;i<Mt->ne;i++){
    if (Gtm->leso[i]==1){
      
      te = Mt->give_elem_type (i);
      
      switch (te){
        case bar2d:{
          //Bar2d->ipvolume(i,0,0);
          break;
        }
        case bar3d:{
          //Bar3d->ipvolume(i,0,0);
          break;
        }
        case planeelementlt:{
          Pelt->ipvolume (i,0,0);
          break;
        }
        case planeelementqt:{
          Peqt->ipvolume (i,0,0);
          break;
        }
        case planeelementrotlt:{
          Perlt->ipvolume (i,0,0);
          break;
        }
        case planeelementlq:{
          Pelq->ipvolume (i,0,0);
          break;
        }
        case planeelementqq:{
          Peqq->ipvolume (i,0,0);
          break;
        }
        case planeelementrotlq:{
          Perlq->ipvolume (i,0,0);
          break;
        }
        case axisymmlqintface:{
          break;
        }
        case planequadinterface:{
          break;
        }
        case lineartet:{
          Ltet->ipvolume (i,0,0);
          break;
        }
        case quadrtet:{
          Qtet->ipvolume (i,0,0);
          break;
        }
        case linearhex:{
          Lhex->ipvolume (i,0,0);
          break;
        }
        case quadrhex:{
          Qhex->ipvolume (i,0,0);
          break;
        }
        case lineartetrot:{
          Ltetrot->ipvolume (i,0,0);
          break;
        }
        case linearhexrot:{
          Lhexrot->ipvolume (i,0,0);
          break;
        }
        case hexintface:{
          break;
        }

        default:
          print_err("unknown element type is required", __FILE__, __LINE__, __func__);
      }
    }
  }
}



/**
  The function computes contributions to internal forces from one finite element, i.e. 
  \f$ \int_{\Omega} \mathbf{B}^T \mathbf{\sigma} \mathrm{dV} \f$.
   
  @param i - element id
  @param lcid - load case id
  @param ifor - %vector of internal forces on one element (output)

  @return The function returns contributions in the parameter ifor.

  Created by JK, 3.11.2006
*/
void elem_internal_forces (long i,long lcid,vector &ifor)
{
  elemtype te;
  te = Mt->give_elem_type (i);
  
  switch (te){
    
    case bar2d:{
      Bar2d->res_internal_forces (lcid,i,ifor);
      break;
    }
    case bar3d:{
      Bar3d->res_internal_forces (lcid,i,ifor);
      break;
    }
    case barq2d:{
      Barq2d->res_internal_forces (lcid,i,ifor);
      break;
    }
    case barq3d:{
      Barq3d->res_internal_forces (lcid,i,ifor);
      break;
    }
    case beam2d:{
      Beam2d->res_internal_forces (lcid,i,ifor);
      break;
    }
    case beam3d:{
      Beam3d->res_internal_forces (lcid,i,ifor);
      break;
    }
    case subsoilbeam:{
      Sbeam->res_internal_forces (lcid,i,ifor);
      break;
    }
    
    case spring_1:{
      Spring->res_internal_forces (lcid,i,ifor);
      break;
    }
    case spring_2:{
      Spring->res_internal_forces (lcid,i,ifor);
      break;
    }
    case spring_3:{
      Spring->res_internal_forces (lcid,i,ifor);
      break;
    }
    case spring_4:{
      Spring->res_internal_forces (lcid,i,ifor);
      break;
    }
    case spring_5:{
      Spring->res_internal_forces (lcid,i,ifor);
      break;
    }
    case spring_6:{
      Spring->res_internal_forces (lcid,i,ifor);
      break;
    }
    
    case planeelementlt:{
      Pelt->res_internal_forces (lcid,i,ifor);
      break;
    }
    case planeelementqt:{
      Peqt->res_internal_forces (lcid,i,ifor);
      break;
    }
    case planeelementrotlt:{
      Perlt->res_internal_forces (lcid,i,ifor);
      break;
    }
    
    case planeelementlq:{
      Pelq->res_internal_forces (lcid,i,ifor);
      break;
    }
    case planeelementqq:{
      Peqq->res_internal_forces (lcid,i,ifor);
      break;
    }
    case planeelementrotlq:{
      Perlq->res_internal_forces (lcid,i,ifor);
      break;
    }
    
    case planeelementsubqt:{
      Pesqt->res_internal_forces (lcid,i,ifor);
      break;
    }
    
    case planequadinterface:{
      Pqifc->res_internal_forces (lcid,i,ifor);
      break;
    }
    
    
    case cctel:{
      Cct->res_internal_forces (lcid,i,ifor);
      break;
    }
    case dktel:{
      Dkt->res_internal_forces (lcid,i,ifor);
      break;
    }
    case dstel:{
      Dst->res_internal_forces (lcid,i,ifor);
      break;
    }
    
    case q4plateel:{
      Q4pl->res_internal_forces (lcid,i,ifor);
      break;
    }
    case dkqel:{
      Dkqelem->res_internal_forces (lcid,i,ifor);
      break;
    }
    case subsoilplateq:{
      Splq->res_internal_forces (lcid,i,ifor);
      break;
    }
    
    case subsoilplatetr:{
      Spltr->res_internal_forces (lcid,i,ifor);
      break;
    }

    case axisymmlt:{
      Asymlt->res_internal_forces (lcid,i,ifor);
      break;
    }
    case axisymmqt:{
      Asymqt->res_internal_forces (lcid,i,ifor);
      break;
    }
    case axisymmlq:{
      Asymlq->res_internal_forces (lcid,i,ifor);
      break;
    }
    case axisymmqq:{
      Asymqq->res_internal_forces (lcid,i,ifor);
      break;
    }
    case axisymmcq:{
      Asymcq->res_internal_forces (lcid,i,ifor);
      break;
    }
    case axisymmlqintface:{
      Asymlqifc->res_internal_forces(lcid, i, ifor);
      break;
    }
    
    case shelltrelem:{
      Shtr->res_internal_forces (lcid,i,ifor);
      break;
    }

    case shelltrmelem:{
      Shtrm->res_internal_forces (lcid,i,ifor);
      break;
    }

    case shellqelem:{
      Shq->res_internal_forces (lcid,i,ifor);
      break;
    }

    case lineartet:{
      Ltet->res_internal_forces (lcid,i,ifor);
      break;
    }
    case quadrtet:{
      Qtet->res_internal_forces (lcid,i,ifor);
      break;
    }
    case linearhex:{
      Lhex->res_internal_forces (lcid,i,ifor);
      break;
    }
    case quadrhex:{
      Qhex->res_internal_forces (lcid,i,ifor);
      break;
    }
    case lineartetrot:{
      Ltetrot->res_internal_forces (lcid,i,ifor);
      break;
    }
    case linearhexrot:{
      Lhexrot->res_internal_forces (lcid,i,ifor);
      break;
    }
    case hexintface:{
      Hexifc->res_internal_forces(lcid, i, ifor);
      break;
    }
    
    
    case particleelem:{
      Pelem->res_internal_forces (i,ifor);
      break;
    }
    case tetralatt:{
      Tlatt->internal_forces (lcid,i,0,0,ifor);
      break;
    }

    default:{
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
    }
  }  
}



/**
   The function integrates arbitrary selected quantity over finite element volume, e.g.
   it performs \int_{\Omega} \mbf{\sigma} d\Omega which results in integrated values that can 
   be used in the homogenization problems.

   
   @param eid - element id
   @param iq - type of integrated quantity (see alias.h)
   @param lcid - number of load case
   @param iv - integrated values (output)
   
   @return The function returns integrated values in the %vector iv

   Created by TKo, 23.1.2015
*/
void elem_volintegration_quant(long eid, integratedquant iq, long lcid, vector &iv)
{
  elemtype te;
  te = Mt->give_elem_type (eid);
  
  switch (te)
    {
    case planeelementlt:{
      Pelt->elem_volintegration_quant(eid, iq, lcid, iv);
      break;
    }
    case planeelementlq:{
      Pelq->elem_volintegration_quant(eid, iq, lcid, iv);
      break;
    }
    case axisymmqq:
      Asymqq->elem_volintegration_quant(eid, iq, lcid, iv);
      break;
    case axisymmcq:
      Asymcq->elem_volintegration_quant(eid, iq, lcid, iv);
      break;
    case lineartet:
      Ltet->elem_volintegration_quant(eid, iq, lcid, iv);
      break;
    default:
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
    }
}



/**
  The function computes stresses at one finite element.
  
  @param i - element id
  @param lcid - load case id
   
  @return The function does not return anything but it changes 
          stress array of integration points on the given element.

  Created by Tomas Koudelka, 10.6.2013
*/
void elem_nlstresses(long i,long lcid)
{
  elemtype te;
  te = Mt->give_elem_type (i);
  
  switch (te){
    
    case bar2d:{
      Bar2d->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case bar3d:{
      Bar3d->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case barq2d:{
      Barq2d->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case barq3d:{
      Barq3d->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case beam2d:
    case beam3d:
    case subsoilbeam:
    case spring_1:
    case spring_2:
    case spring_3:
    case spring_4:
    case spring_5:
    case spring_6:
      break;
    case planeelementlt:{
      Pelt->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case planeelementqt:{
      Peqt->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case planeelementrotlt:{
      Perlt->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case planeelementlq:{
      Pelq->compute_nlstress (lcid,i,0,0);
      break;
    }
    case planeelementqq:{
      Peqq->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case planeelementrotlq:{
      Perlq->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case planeelementsubqt:{
      Pesqt->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case planequadinterface:{
      Pqifc->compute_nlstress(lcid, i, 0, 0);    
      break;
    }
    case cctel:{
      Cct->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case q4plateel:{
      Q4pl->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case subsoilplateq:{
      Splq->compute_nlstress(lcid, i, 0, 0);
      break;
    }
      /*
        case subsoilplatetr:{
        Spltr->compute_nlstress(lcid, i, 0, 0);
        break;
        }
        case shelltrelem:{
        Shtr->compute_nlstress(lcid, i, 0, 0);
        break;
        }
        case shellqelem:{
        Shq->compute_nlstress(lcid, i, 0, 0);
        break;
        }
      */
    case axisymmlt:{
      Asymlt->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case axisymmlq:{
      Asymlq->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case axisymmqq:{
      Asymqq->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case axisymmcq:{
      Asymcq->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case axisymmlqintface:{
      Asymlqifc->compute_nlstress(lcid, i, 0, 0);
      break;
    }

    case lineartet:{
      Ltet->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case quadrtet:{
      Qtet->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case linearhex:{
      Lhex->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case quadrhex:{
      Qhex->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case lineartetrot:{
      Ltetrot->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case linearhexrot:{
      Lhexrot->compute_nlstress(lcid, i, 0, 0);
      break;
    }
    case hexintface:{
      Hexifc->compute_nlstress(lcid, i, 0, 0);
      break;
    }

    case particleelem:
      break;
    default:
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function computes local values which will be averaged for nonlocal material models
  at one finite element.
  
  @param i - element id
  @param lcid - load case id
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 7.2008
*/
void elem_local_values (long i,long lcid)
{
  elemtype te;
  te = Mt->give_elem_type (i);
  
  switch (te){
    
    case bar2d:{
      Bar2d->local_values(lcid, i, 0, 0);
      break;
    }
    case bar3d:{
      Bar3d->local_values(lcid, i, 0, 0);
      break;
    }
    case barq2d:{
      Barq2d->local_values(lcid, i, 0, 0);
      break;
    }
    case barq3d:{
      Barq3d->local_values(lcid, i, 0, 0);
      break;
    }
    case beam2d:
    case beam3d:
    case subsoilbeam:
    case spring_1:
    case spring_2:
    case spring_3:
    case spring_4:
    case spring_5:
    case spring_6:
      break;
    case planeelementlt:{
      Pelt->local_values(lcid, i, 0, 0);
      break;
    }
    case planeelementqt:{
      Peqt->local_values(lcid, i, 0, 0);
      break;
    }
    case planeelementrotlt:{
      Perlt->local_values(lcid, i, 0, 0);
      break;
    }
    case planeelementlq:{
      Pelq->local_values (lcid,i,0,0);
      break;
    }
    case planeelementqq:{
      Peqq->local_values(lcid, i, 0, 0);
      break;
    }
    case planeelementrotlq:{
      Perlq->local_values(lcid, i, 0, 0);
      break;
    }
    case planeelementsubqt:{
      Pesqt->local_values(lcid, i, 0, 0);
      break;
    }
    case axisymmlqintface:{
      break;
    }
    case planequadinterface:{
      break;
    }
    case cctel:{
      Cct->local_values(lcid, i, 0, 0);
      break;
    }
    case q4plateel:{
      Q4pl->local_values(lcid, i, 0, 0);
      break;
    }
    case subsoilplateq:{
      Splq->local_values(lcid, i, 0, 0);
      break;
    }
      /*
        case subsoilplatetr:{
        Spltr->local_values(lcid, i, 0, 0);
        break;
        }
        case shelltrelem:{
        Shtr->local_values(lcid, i, 0, 0);
        break;
        }
        case shellqelem:{
        Shq->local_values(lcid, i, 0, 0);
        break;
        }
      */
    case axisymmlt:{
      Asymlt->local_values(lcid, i, 0, 0);
      break;
    }
    case axisymmlq:{
      Asymlq->local_values(lcid, i, 0, 0);
      break;
    }
    case axisymmqq:{
      Asymqq->local_values(lcid, i, 0, 0);
      break;
    }
    case axisymmcq:{
      Asymcq->local_values(lcid, i, 0, 0);
      break;
    }

    case lineartet:{
      Ltet->local_values(lcid, i, 0, 0);
      break;
    }
    case quadrtet:{
      Qtet->local_values(lcid, i, 0, 0);
      break;
    }
    case linearhex:{
      Lhex->local_values(lcid, i, 0, 0);
      break;
    }
    case quadrhex:{
      Qhex->local_values(lcid, i, 0, 0);
      break;
    }
    case lineartetrot:{
      Ltetrot->local_values(lcid, i, 0, 0);
      break;
    }
    case linearhexrot:{
      Lhexrot->local_values(lcid, i, 0, 0);
      break;
    }
    case hexintface:{
      break;
    }

    case particleelem:
      break;
    default:
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function computes contributions to internal forces from one finite element.
  Contributions are evaluated nonlocally.
   
  @param i - element id
  @param lcid - load case id
  @param ifor - %vector of internal forces on one element (output)

  @return The function returns contributions in the parameter ifor.

  Created by JK, 3.11.2006
*/
void elem_nonloc_internal_forces (long i,long lcid,vector &ifor)
{
  elemtype te;
  te = Mt->give_elem_type (i);
  
  switch (te){
    
    case bar2d:{
      Bar2d->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case bar3d:{
      Bar3d->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case barq2d:{
      Barq2d->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case barq3d:{
      Barq3d->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case beam2d:{
      Beam2d->res_internal_forces (lcid,i,ifor);
      break;
    }
    case beam3d:{
      Beam3d->res_internal_forces (lcid,i,ifor);
      break;
    }
    case subsoilbeam:{
      Sbeam->res_internal_forces (lcid,i,ifor);
      break;
    }
    
    case spring_1:{
      Spring->res_internal_forces (lcid,i,ifor);
      break;
    }
    case spring_2:{
      Spring->res_internal_forces (lcid,i,ifor);
      break;
    }
    case spring_3:{
      Spring->res_internal_forces (lcid,i,ifor);
      break;
    }
    case spring_4:{
      Spring->res_internal_forces (lcid,i,ifor);
      break;
    }
    case spring_5:{
      Spring->res_internal_forces (lcid,i,ifor);
      break;
    }
    case spring_6:{
      Spring->res_internal_forces (lcid,i,ifor);
      break;
    }
    
    case planeelementlt:{
      Pelt->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case planeelementqt:{
      Peqt->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case planeelementrotlt:{
      Perlt->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    
    case planeelementlq:{
      Pelq->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case planeelementqq:{
      Peqq->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case planeelementrotlq:{
      Perlq->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    
    case planeelementsubqt:{
      Pesqt->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    
    case planequadinterface:{
      Pqifc->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    
    case cctel:{
      Cct->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    
    case q4plateel:{
      Q4pl->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case subsoilplateq:{
      Splq->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
      /*
        case subsoilplatetr:{
        Spltr->lres_nonloc_internal_forces (lcid,i,ifor);
        break;
        }
        case shelltrelem:{
        Shtr->res_nonloc_internal_forces (lcid,i,ifor);
        break;
        }
        case shellqelem:{
        Shq->res_nonloc_internal_forces (lcid,i,ifor);
        break;
        }
      */    

    case axisymmlt:{
      Asymlt->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case axisymmlq:{
      Asymlq->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case axisymmqq:{
      Asymqq->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case axisymmcq:{
      Asymcq->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case axisymmlqintface:{
      Asymlqifc->res_nonloc_internal_forces(lcid, i, ifor);
      break;
    }

    case lineartet:{
      Ltet->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case quadrtet:{
      Qtet->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case linearhex:{
      Lhex->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case quadrhex:{
      Qhex->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case lineartetrot:{
      Ltetrot->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case linearhexrot:{
      Lhexrot->res_nonloc_internal_forces (lcid,i,ifor);
      break;
    }
    case hexintface:{
      Hexifc->res_nonloc_internal_forces(lcid, i, ifor);
      break;
    }
    
    case particleelem:{
      break;
    }
    
    default:
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function computes contributions to increments of internal forces
  from one finite element.
 
  @param i - element id
  @param lcid - load case id
  @param ifor - %vector of increments of internal forces on one element (output)
  
  @return The function returns contributions in the parameter ifor.

  Created by JK, 29.4.2008
*/
void elem_incr_internal_forces (long i,long lcid,vector &ifor)
{
  elemtype te;
  te = Mt->give_elem_type (i);
  
  switch (te){
    
    case bar2d:{
      Bar2d->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case bar3d:{
      Bar3d->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case barq2d:{
      Barq2d->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case barq3d:{
      Barq3d->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
      /*
        case beam2d:{
        Beam2d->res_internal_forces (lcid,i,ifor);
        break;
        }
        case beam3d:{
        Beam3d->res_internal_forces (lcid,i,ifor);
        break;
        }
        case subsoilbeam:{
        Sbeam->res_internal_forces (lcid,i,ifor);
        break;
        }

      */  
    
    case spring_1:
    case spring_2:
    case spring_3:
    case spring_4:
    case spring_5:
    case spring_6:
      break;
    
    case planeelementlt:{
      Pelt->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case planeelementqt:{
      Peqt->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case planeelementrotlt:{
      Perlt->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    
    case planeelementlq:{
      Pelq->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case planeelementqq:{
      Peqq->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case planeelementrotlq:{
      Perlq->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    
    case planeelementsubqt:{
      Pesqt->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    
    case planequadinterface:{
      break;
    }
    
    case cctel:{
      Cct->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    
    case q4plateel:{
      Q4pl->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case subsoilplateq:{
      Splq->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    
    case axisymmlt:{
      Asymlt->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case axisymmqt:{
      Asymqt->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case axisymmlq:{
      Asymlq->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case axisymmqq:{
      Asymqq->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case axisymmcq:{
      Asymcq->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case axisymmlqintface:{
      Asymlqifc->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    
      /*
        case shelltrelem:{
        Shtr->res_incr_internal_forces (lcid,i,ifor);
        break;
        }
      */

    case lineartet:{
      Ltet->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case quadrtet:{
      Qtet->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case linearhex:{
      Lhex->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case quadrhex:{
      Qhex->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case lineartetrot:{
      Ltetrot->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case linearhexrot:{
      Lhexrot->res_incr_internal_forces (lcid,i,ifor);
      break;
    }
    case hexintface:{
      break;
    }
      /*    
            case particleelem:{
            Pelem->res_internal_forces (i,ifor);
            break;
            }
      */
    default:
      print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
}



/**
  The function computes nodal forces caused by eigenstrains on one element.
   
  @param lcid - load case id
  @param eid - element id
  @param nfor - %vector of nodal values (output)
  
  @return The function returns contributions in the parameter nfor.
   
  Created by JK,
*/
void elem_eigstrain_forces (long lcid,long eid,vector &nfor)
{
  switch (Mt->give_elem_type(eid)){
    
    case bar2d:{
      Bar2d->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case bar3d:{
      Bar3d->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case barq2d:{
      Barq2d->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case barq3d:{
      Barq3d->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }

    case planeelementlt:{
      Pelt->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case planeelementqt:{
      Peqt->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case planeelementrotlt:{
      Perlt->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }

    case planeelementlq:{
      Pelq->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case planeelementqq:{
      Peqq->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case planeelementrotlq:{
      Perlq->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case planeelementsubqt:{
      Pesqt->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case planequadinterface:{
      break;
    }

    case cctel:{
      Cct->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case q4plateel:{
      Q4pl->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case subsoilplateq:{
      Splq->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }


    case axisymmlt:{
      Asymlt->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case axisymmqt:{
      Asymqt->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case axisymmlq:{
      Asymlq->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case axisymmqq:{
      Asymqq->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case axisymmcq:{
      Asymcq->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case axisymmlqintface:{
      Asymlqifc->res_eigstrain_forces(lcid, eid, nfor);
      break;
    }
    
    case lineartet:{
      Ltet->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case quadrtet:{
      Qtet->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case linearhex:{
      Lhex->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case quadrhex:{
      Qhex->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case lineartetrot:{
      Ltetrot->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case linearhexrot:{
      Lhexrot->res_eigstrain_forces (lcid,eid,nfor);
      break;
    }
    case hexintface:{
      break;
    }
    case spring_1:
    case spring_2:
    case spring_3:
    case spring_4:
    case spring_5:
    case spring_6:
      break;

    default:
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function computes stiffness %matrix of required element.

  @param lcid - load case id
  @param eid - element id
  @param sm - stiffness %matrix of the element (output)

  @return The function returns assembled %matrix in the parameter sm.

  Created by JK,
*/
void stiffmat (long lcid,long eid,matrix &sm)
{
  elemtype te;
  
  te = Mt->give_elem_type (eid);

  switch (te){

    case bar2d:{
      Bar2d->res_stiffness_matrix (eid,sm);
      break;
    }
    case bar3d:{
      Bar3d->res_stiffness_matrix (eid,sm);
      break;
    }
    case barq2d:{
      Barq2d->res_stiffness_matrix (eid,sm);
      break;
    }
    case barq3d:{
      Barq3d->res_stiffness_matrix (eid,sm);
      break;
    }
    case beam2d:{
      Beam2d->res_stiffness_matrix (eid,sm);
      break;
    }
    case beam3d:{
      Beam3d->res_stiffness_matrix (eid,sm);
      break;
    }
    case beamg3d:{
      Beam3dg->res_stiffness_matrix (eid,sm);
      break;
    }
    case subsoilbeam:{
      Sbeam->res_stiffness_matrix (eid,sm);
      break;
    }

    case spring_1:{
      Spring->res_stiffness_matrix (eid,sm);
      break;
    }
    case spring_2:{
      Spring->res_stiffness_matrix (eid,sm);
      break;
    }
    case spring_3:{
      Spring->res_stiffness_matrix (eid,sm);
      break;
    }
    case spring_4:{
      Spring->res_stiffness_matrix (eid,sm);
      break;
    }
    case spring_5:{
      Spring->res_stiffness_matrix (eid,sm);
      break;
    }
    case spring_6:{
      Spring->res_stiffness_matrix (eid,sm);
      break;
    }

    case planeelementlt:{
      Pelt->res_stiffness_matrix (eid,sm);
      break;
    }
    case planeelementqt:{
      Peqt->res_stiffness_matrix (eid,sm);
      break;
    }
    case planeelementrotlt:{
      Perlt->res_stiffness_matrix (eid,sm);
      break;
    }

    case planeelementlq:{
      Pelq->res_stiffness_matrix (lcid,eid,sm);
      break;
    }
    case planeelementqq:{
      Peqq->res_stiffness_matrix (eid,sm);
      break;
    }
    case planeelementrotlq:{
      Perlq->res_stiffness_matrix (eid,sm);
      break;
    }
    
    case planeelementsubqt:{
      Pesqt->res_stiffness_matrix (eid,sm);
      break;
    }

    case planequadinterface:{
      Pqifc->res_stiffness_matrix (eid,sm);
      break;
    }
    
    case cctel:{
      Cct->res_stiffness_matrix (eid,sm);
      break;
    }
    case dktel:{
      Dkt->res_stiffness_matrix (eid,sm);
      break;
    }
    case dstel:{
      Dst->res_stiffness_matrix (eid,sm);
      break;
    }
    case q4plateel:{
      Q4pl->res_stiffness_matrix (eid,sm);
      break;
    }
    case argyristr:{
      //Argtr->res_stiffness_matrix (eid,sm);
      Argtrpl->res_stiffness_matrix (eid,sm);
      break;
    }
    case quadkirch:{
      Qkirch->res_stiffness_matrix (eid,sm);
      break;
    }
    case dkqel:{
      Dkqelem->res_stiffness_matrix (eid,sm);
      break;
    }

    case subsoilplatetr:{
      Spltr->res_stiffness_matrix (eid,sm);
      break;
    }
    case subsoilplateq:{
      Splq->res_stiffness_matrix (eid,sm);
      break;
    }
    
    case axisymmlt:{
      Asymlt->res_stiffness_matrix (eid,sm);
      break;
    }
    case axisymmqt:{
      Asymqt->res_stiffness_matrix (eid,sm);
      break;
    }
    case axisymmlq:{
      Asymlq->res_stiffness_matrix (eid,sm);
      break;
    }
    case axisymmqq:{
      Asymqq->res_stiffness_matrix (eid,sm);
      break;
    }
    case axisymmcq:{
      Asymcq->res_stiffness_matrix (eid,sm);
      break;
    }
    case axisymmlqintface:{
      Asymlqifc->res_stiffness_matrix(eid, sm);
      break;
    }
    
    case shelltrelem:{
      Shtr->res_stiffness_matrix (eid,sm);
      break;
    }
    case shelltrmelem:{
      Shtrm->res_stiffness_matrix (eid,sm);
      break;
    }
    case shellqelem:{
      Shq->res_stiffness_matrix (eid,sm);
      break;
    }

    case lineartet:{
      Ltet->res_stiffness_matrix (eid,sm);
      break;
    }
    case quadrtet:{
      Qtet->res_stiffness_matrix (eid,sm);
      break;
    }
    case linearhex:{
      Lhex->res_stiffness_matrix (lcid,eid,sm);
      break;
    }
    case quadrhex:{
      Qhex->res_stiffness_matrix (eid,sm);
      break;
    }
    case lineartetrot:{
      Ltetrot->res_stiffness_matrix (eid,sm);
      break;
    }
    case linearhexrot:{
      Lhexrot->res_stiffness_matrix (lcid,eid,sm);
      break;
    }
    case linearwed:{
      Lwed->res_stiffness_matrix (eid,sm);
      break;
    }
    case quadrwed:{
      Qwed->res_stiffness_matrix (eid,sm);
      break;
    }
    case hexintface:{
      Hexifc->res_stiffness_matrix (eid, sm);
      break;
    }
    
    case particleelem:{
      Pelem->res_stiffness_matrix (eid,sm);
      break;
    }
    default:
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}


/**
  The function computes bd %matrix of required element.

  @param lcid - load case id
  @param eid - element id
  @param bd - bd %matrix of the element (output)
  
  @return The function returns assembled %matrix in the parameter bd.
  
  Created by TKr, 23/10/2019
*/
void bdmatrix (long /*lcid*/,long eid,matrix &bd)
{
  elemtype te;
  
  te = Mt->give_elem_type (eid);
  
  switch (te){
    
  case planeelementlt:{
    Pelt->bd_matrix (eid,0,0,bd);
    break;
  }
  case planeelementlq:{
    Pelq->bd_matrix (eid,0,0,bd);
    break;
  }
  case lineartet:{
    Ltet->bd_matrix (eid,0,0,bd);
    break;
  }
  case linearhex:{
    Lhex->bd_matrix (eid,0,0,bd);
    break;
  }
  default:
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}


/**
  The function computes dd %matrix of required element.

  @param lcid - load case id
  @param eid - element id
  @param dd - dd %matrix of the element (output)

  @return The function returns assembled %matrix in the parameter dd.

  Created by TKr, 23/10/2019
*/
void ddmatrix (long /*lcid*/,long eid,matrix &dd)
{
  elemtype te;
  
  te = Mt->give_elem_type (eid);

  switch (te){
    
  case planeelementlt:{
    Pelt->dd_matrix (eid,0,0,dd);
    break;
  }
  case planeelementlq:{
    Pelq->dd_matrix (eid,0,0,dd);
    break;
  }
  case lineartet:{
    Ltet->dd_matrix (eid,0,0,dd);
    break;
  }
  case linearhex:{
    Lhex->dd_matrix (eid,0,0,dd);
    break;
  }
  default:
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}

/**
  The function computes bd %matrix of required element.

  @param lcid - load case id
  @param eid - element id
  @param db - db %matrix of the element (output)
  
  @return The function returns assembled %matrix in the parameter db.
  
  Created by TKr, 23/10/2019
*/
void dbmatrix (long /*lcid*/,long eid,matrix &/*db*/)
{
  elemtype te;
  
  te = Mt->give_elem_type (eid);
  
  switch (te){
    
  case planeelementlt:{
    //Pelt->db_matrix (eid,0,0,db);
    break;
  }
  case planeelementlq:{
    //Pelq->db_matrix (eid,0,0,db);
    break;
  }
  case lineartet:{
    //Ltet->db_matrix (eid,0,0,db);
    break;
  }
  case linearhex:{
    //Lhex->db_matrix (eid,0,0,db);
    break;
  }
  default:
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function computes mass %matrix of required element.

  @param lcid - load case id
  @param eid - element id
  @param mm - mass %matrix of the element (output)

  @return The function returns assembled %matrix in the parameter mm.

  Created by JK,
*/
void massmat (long /*lcid*/,long eid,matrix &mm)
{
  elemtype te;
  
  te = Mt->give_elem_type (eid);

  switch (te){

    case bar2d:{
      Bar2d->res_mass_matrix (eid,mm);
      break;
    }
    case barq2d:{
      Barq2d->res_mass_matrix (eid,mm);
      break;
    }
    case beam2d:{
      Beam2d->res_mass_matrix (eid,mm);
      break;
    }
    case bar3d:{
      Bar3d->res_mass_matrix (eid,mm);
      break;
    }
    case beam3d:{
      Beam3d->res_mass_matrix (eid,mm);
      break;
    }
    case barq3d:{
      Barq3d->res_mass_matrix (eid,mm);
      break;
    }
    case spring_1:
    case spring_2:
    case spring_3:
    case spring_4:
    case spring_5:
    case spring_6:
      break;
  
    case planeelementlt:{
      Pelt->res_mass_matrix (eid,mm);
      break;
    }
    case planeelementqt:{
      Peqt->res_mass_matrix (eid,mm);
      break;
    }
    case planeelementrotlt:{
      Perlt->res_mass_matrix (eid,mm);
      break;
    }
    case planeelementlq:{
      Pelq->res_mass_matrix (eid,mm);
      break;
    }
    case planeelementqq:{
      Peqq->res_mass_matrix (eid,mm);
      break;
    }
    case planeelementrotlq:{
      Perlq->res_mass_matrix (eid,mm);
      break;
    }

    case planeelementsubqt:{
      Pesqt->res_mass_matrix (eid,mm);
      break;
    }

    case cctel:{
      Cct->res_mass_matrix (eid,mm);
      break;
    }

    case axisymmlt:{
      Asymlt->res_mass_matrix (eid,mm);
      break;
    }
    case axisymmqt:{
      Asymqt->res_mass_matrix (eid,mm);
      break;
    }
    case axisymmlq:{
      Asymlq->res_mass_matrix (eid,mm);
      break;
    }
    case axisymmqq:{
      Asymqq->mass_matrix (eid,mm);
      break;
    }
    case axisymmcq:{
      Asymcq->mass_matrix (eid,mm);
      break;
    }
    case axisymmlqintface:{
      break;
    }

    case lineartet:{
      Ltet->res_mass_matrix (eid,mm);
      break;
    }
    case quadrtet:{
      Qtet->res_mass_matrix (eid,mm);
      break;
    }
    case linearhex:{
      Lhex->res_mass_matrix (eid,mm);
      break;
    }
    case quadrhex:{
      Qhex->res_mass_matrix (eid,mm);
      break;
    }
    case lineartetrot:{
      Ltetrot->mass_matrix (eid,mm);
      break;
    }
    case linearhexrot:{
      Lhexrot->res_mass_matrix (eid,mm);
      break;
    }
    
    case tetralatt:{
      Tlatt->create_facets (eid,mm);
      break;
    }

    default:
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function computes initial stiffness %matrix of required element.
   
  @param lcid - load case id
  @param eid - element id
  @param sm - initial stiffness %matrix of the element (output)

  @return The function returns assembled %matrix in the parameter sm.

  Created by JK,
*/
void initstiffmat (long lcid,long eid,matrix &sm)
{
  elemtype te;
  
  te = Mt->give_elem_type (eid);

  switch (te){

  case beam2d:{
    Beam2d->initstr_matrix_expl (lcid,eid,0,0,sm);
    break;
  }

  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }

}

/**
  The function computes damping %matrix of required element.

  @param lcid - load case id
  @param eid - element id
  @param dm - damping %matrix of the element (output)

  @return The function returns assembled %matrix in the parameter mm.

  Created by JK, 20.11.2017
*/
void dampmat (long /*lcid*/,long eid,matrix &dm)
{
  elemtype te;
  
  te = Mt->give_elem_type (eid);
  
  switch (te){
    
    /*
  case bar2d:{
    Bar2d->mass_matrix (eid,mm);
    break;
  }
  case barq2d:{
    Barq2d->res_mass_matrix (eid,mm);
    break;
  }
  case beam2d:{
    Beam2d->res_mass_matrix (eid,mm);
    break;
  }
  case bar3d:{
    Bar3d->res_mass_matrix (eid,mm);
    break;
  }
  case beam3d:{
    Beam3d->res_mass_matrix (eid,mm);
    break;
  }

  case spring_1:
  case spring_2:
  case spring_3:
  case spring_4:
  case spring_5:
  case spring_6:
    break;
  
  case planeelementlt:{
    Pelt->res_mass_matrix (eid,mm);
    break;
  }
  case planeelementqt:{
    Peqt->mass_matrix (eid,mm);
    break;
  }
  case planeelementrotlt:{
    Perlt->res_mass_matrix (eid,mm);
    break;
  }
  */
  case planeelementlq:{
    Pelq->res_damping_matrix (eid,dm);
    break;
  }
    /*
  case planeelementqq:{
    Peqq->res_mass_matrix (eid,mm);
    break;
  }
  case planeelementrotlq:{
    Perlq->res_mass_matrix (eid,mm);
    break;
  }

  case planeelementsubqt:{
    Pesqt->res_mass_matrix (eid,mm);
    break;
  }

  case cctel:{
    Cct->res_mass_matrix (eid,mm);
    break;
  }

  case axisymmlt:{
    Asymlt->mass_matrix (eid,mm);
    break;
  }
  case axisymmlq:{
    Asymlq->mass_matrix (eid,mm);
    break;
  }
  case axisymmqq:{
    Asymqq->mass_matrix (eid,mm);
    break;
  }
  case axisymmcq:{
    Asymcq->mass_matrix (eid,mm);
    break;
  }

  case lineartet:{
    Ltet->mass_matrix (eid,mm);
    break;
  }
  case quadrtet:{
    Qtet->mass_matrix (eid,mm);
    break;
  }
  case linearhex:{
    Lhex->mass_matrix (eid,mm);
    break;
  }
  case quadrhex:{
    Qhex->mass_matrix (eid,mm);
    break;
  }
  case lineartetrot:{
    Ltetrot->mass_matrix (eid,mm);
    break;
  }
  case linearhexrot:{
    Lhexrot->mass_matrix (eid,mm);
    break;
  }
  */
  default:
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}


/**
  The function computes load %matrix of required element.
   
  @param eid - element id
  @param sm - load %matrix of the element
   
  @return The function returns assembled %matrix in the parameter lm.

  Created by JK, 12.10.2008
*/
void loadmat(long eid, matrix &lm)
{
  elemtype et = Mt->give_elem_type (eid);
  switch (et){
    case bar2d:{
      Bar2d->res_load_matrix (eid,lm);
      break;
    }
    case barq2d:{
      Barq2d->res_load_matrix (eid,lm);
      break;
    }
    case beam3d:{
      Beam3d->load_matrix (eid,lm);
      break;
    }
    case bar3d:{
      Bar3d->res_load_matrix (eid,lm);
      break;
    }
    
    case planeelementlt:{
      Pelt->res_load_matrix (eid,lm);
      break;
    }
    case planeelementqt:{
      Peqt->load_matrix (eid,lm);
      break;
    }
    case planeelementlq:{
      Pelq->res_load_matrix (eid,lm);
      break;
    }
    case planeelementqq:{
      Peqq->res_load_matrix (eid,lm);
      break;
    }
    case planeelementrotlt:{
      Perlt->load_matrix (eid,lm);
      break;
    }
    case planeelementrotlq:{
      Perlq->res_load_matrix (eid,lm);
      break;
    }
    
    case cctel:{
      Cct->load_matrix (eid,lm);
      break;
    }
    case q4plateel:{
      Q4pl->load_matrix (eid,lm);
      break;
    }
    
    case axisymmlt:{
      Asymlt->load_matrix (eid,lm);
      break;
    }
    case axisymmlq:{
      Asymlq->load_matrix (eid,lm);
      break;
    }
    case axisymmqq:{
      Asymqq->load_matrix (eid,lm);
      break;
    }
    case axisymmcq:{
      Asymcq->load_matrix (eid,lm);
      break;
    }
    case axisymmlqintface:{
      break;
    }
    
    case lineartet:{
      Ltet->load_matrix (eid,lm);
      break;
    }
    case quadrtet:{
      Qtet->load_matrix (eid,lm);
      break;
    }
    case linearhex:{
      Lhex->load_matrix (eid,lm);
      break;
    }
    case quadrhex:{
      Qhex->load_matrix (eid,lm);
      break;
    }
    default:
      print_err("unknown element type %d is required", __FILE__, __LINE__, __func__, et);
  }
}



/**
  The function assembles transformation %matrix from local nodal coordinate system to the global coordinate system x_g = T x_l
  of the required element eid.

  @param eid - element id [in]
  @param enodes - %vector of nodal ids of the given element [in]
  @param tmat - transformation %matrix of the element [out]

  @return The function returns assembled %matrix in the parameter tmat.

  Created by Tomas Koudelka, 24.11.2017
*/
void elem_transf_matrix (long eid, ivector &enodes, matrix &tmat)
{
  elemtype te;
  
  te = Mt->give_elem_type (eid);

  switch (te){

    case bar2d:{
      Bar2d->transf_matrix(enodes, tmat);
      break;
    }
    case bar3d:{
      Bar3d->transf_matrix(enodes, tmat);
      break;
    }
    case barq2d:{
      Barq2d->transf_matrix(enodes, tmat);
      break;
    }
    case barq3d:{
      Barq3d->transf_matrix(enodes, tmat);
      break;
    }
    case beam2d:{
      Beam2d->transf_matrix(enodes, tmat);
      break;
    }
    case beam3d:{
      Beam3d->transf_matrix(enodes, tmat);
      break;
    }
    case beamg3d:{
      Beam3dg->transf_matrix(enodes, tmat);
      break;
    }
    case subsoilbeam:{
      Sbeam->transf_matrix(enodes, tmat);
      break;
    }

    case planeelementlt:{
      Pelt->transf_matrix(enodes, tmat);
      break;
    }
    case planeelementqt:{
      Peqt->transf_matrix(enodes, tmat);
      break;
    }
    case planeelementrotlt:{
      Perlt->transf_matrix(enodes, tmat);
      break;
    }

    case planeelementlq:{
      Pelq->transf_matrix(enodes, tmat);
      break;
    }
    case planeelementqq:{
      Peqq->transf_matrix(enodes, tmat);
      break;
    }
    case planeelementrotlq:{
      Perlq->transf_matrix(enodes, tmat);
      break;
    }
    
    case planeelementsubqt:{
      Pesqt->transf_matrix(enodes, tmat);
      break;
    }

    case planequadinterface:{
      Pqifc->transf_matrix(enodes, tmat);
      break;
    }
    
    case cctel:{
      Cct->transf_matrix(enodes, tmat);
      break;
    }
    case dktel:{
      Dkt->transf_matrix(enodes, tmat);
      break;
    }
    case dstel:{
      Dst->transf_matrix(enodes, tmat);
      break;
    }
    
    case subsoilplatetr:{
      Spltr->transf_matrix(enodes, tmat);
      break;
    }
    
    case axisymmlt:{
      Asymlt->transf_matrix(enodes, tmat);
      break;
    }
    case axisymmlq:{
      Asymlq->transf_matrix(enodes, tmat);
      break;
    }
    case axisymmqq:{
      Asymqq->transf_matrix(enodes, tmat);
      break;
    }
    case axisymmcq:{
      Asymcq->transf_matrix(enodes, tmat);
      break;
    }
    case axisymmlqintface:{
      Asymlqifc->transf_matrix(enodes, tmat);
      break;
    }

    case shelltrelem:{
      Shtr->node_transf_matrix(enodes, tmat);
      break;
    }
    case shelltrmelem:{
      Shtrm->node_transf_matrix(enodes, tmat);
      break;
    }
    case shellqelem:{
      Shq->node_transf_matrix (enodes, tmat);
      break;
    }

    case lineartet:{
      Ltet->transf_matrix(enodes, tmat);
      break;
    }
    case quadrtet:{
      Qtet->transf_matrix(enodes, tmat);
      break;
    }
    case linearhex:{
      Lhex->transf_matrix(enodes, tmat);
      break;
    }
    case quadrhex:{
      Qhex->transf_matrix(enodes, tmat);
      break;
    }
    case lineartetrot:{
      Ltetrot->transf_matrix(enodes, tmat);
      break;
    }
    case linearhexrot:{
      Lhexrot->transf_matrix(enodes, tmat);
      break;
    }
    case linearwed:{
      Lwed->transf_matrix(enodes, tmat);
      break;
    }
    case quadrwed:{
      Qwed->transf_matrix(enodes, tmat);
      break;
    }
    case hexintface:{
      Hexifc->transf_matrix(enodes, tmat);
      break;
    }
    default:
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function returns non-transport=mechanical quantities at element nodes. Nodal
  values are copied from the closest integration point on element.
  It is used for passing required mechanical quantities by TRFEL in coupled problems.

  @param eid - element id
  @param nodval - %vector of nodal values
  @param qt - type of mechanical quantity
   
  @return The function does not return anything but it stores nodal values in parameter nodval.
 
  Created 12/06/2012 TKr
  Modified by Tomas Koudelka 9.10.2013
*/
void elem_mechq_nodval(long eid, vector &nodval, nontransquant qt)
{
  elemtype te;

  te = Mt->give_elem_type (eid);

  switch (te){
  case bar2d:{
    Bar2d->mechq_nodval (eid,nodval,qt);
    break;
  }
  case planeelementlq:{
    Pelq->mechq_nodval (eid,nodval,qt);
    break;
  }
  case planeelementqq:{
    Peqq->mechq_nodval (eid,nodval,qt);
    break;
  }
  case axisymmlq:{
    Asymlq->mechq_nodval (eid,nodval,qt);
    break;
  }
  case axisymmqq:{
    Asymqq->mechq_nodval (eid,nodval,qt);
    break;
  }
  case axisymmcq:{
    Asymcq->mechq_nodval (eid,nodval,qt);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
   The function returns non-transport=mechanical quantities at element nodes. Nodal
   values are copied from the closest integration point on element.
   It is used for passing required mechanical quantities by TRFEL in coupled problems.
   
   @param eid - element id
   @param nodval - %vector of nodal values
   @param qt - type of mechanical quantity
   
   @return The function does not return anything but it stores nodal values in parameter nodval.
   
   01/11/2016 by TKr according to intpointval2
*/
void elem_mechq_nodval2(long eid, vector &nodval, nontransquant qt)
{
  elemtype te;
  
  te = Mt->give_elem_type (eid);
  
  switch (te){
  case bar2d:{
    Bar2d->mechq_nodval (eid,nodval,qt);
    break;
  }
  case barq2d:{
    Barq2d->mechq_nodval (eid,nodval,qt);
    break;
  }
  case planeelementlq:{
    Pelq->mechq_nodval (eid,nodval,qt);
    break;
  }
  case planeelementqq:{
    Peqq->mechq_nodval2 (eid,nodval,qt);
    break;
  }
  case axisymmlq:{
    Asymlq->mechq_nodval (eid,nodval,qt);
    break;
  }
  case axisymmqq:{
    Asymqq->mechq_nodval2 (eid,nodval,qt);
    break;
  }
  case axisymmcq:{
    Asymcq->mechq_nodval2 (eid,nodval,qt);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  The function returns non-transport=mechanical quantities from element nodes.
  Nodal values are computed directly at nodes according to strains.
  It is used for passing required mechanical quantities by TRFEL in coupled problems.

  @param eid - element id
  @param nodval - %vector of nodal values of all required quantities, i.e., 
                  nodal value of i-th quantity in j-th node is given by nodval[i*nne+j] where nne 
                  is the number of nodes on eid-th element.
  @param ncne - number of computed nodes on element (only first ncne of nodes is calculated)
  @param nq - number of required mechanical quantities
  @param qt - array of types of required mechanical quantities
   
  @return The function does not return anything but it stores nodal values in parameter nodval.
 
  Created 12/06/2012 TKr
  Modified by Tomas Koudelka 9.10.2013
*/
void elem_mechq_nodval_comp(long eid, vector &nodval, long ncne, long nq, nontransquant *qt)
{
  elemtype te;

  te = Mt->give_elem_type (eid);

  switch (te){
  case bar2d:{
    Bar2d->mechq_nodval_comp (eid,nodval,ncne,nq,qt);
    break;
  }
  case barq2d:{
    Barq2d->mechq_nodval_comp (eid,nodval,ncne,nq,qt);
    break;
  }
  case planeelementlq:{
    Pelq->mechq_nodval_comp (eid,nodval,ncne,nq,qt);
    break;
  }
  case planeelementqq:{
    Peqq->mechq_nodval_comp (eid,nodval,ncne,nq,qt);
    break;
  }
  case axisymmlq:{
    Asymlq->mechq_nodval_comp (eid,nodval,ncne,nq,qt);
    break;
  }
  case axisymmqq:{
    Asymqq->mechq_nodval_comp (eid,nodval,ncne,nq,qt);
    break;
  }
  case axisymmcq:{
    Asymcq->mechq_nodval_comp (eid,nodval,ncne,nq,qt);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}

/**
  The function returns non-transport=mechanical quantities from element nodes.
  Nodal values are computed directly at nodes according to strains.
  It is used for passing required mechanical quantities by TRFEL in coupled problems.

  @param eid - element id
  @param nodval - %vector of nodal values of all required quantities, i.e., 
                  nodal value of i-th quantity in j-th node is given by nodval[i*nne+j] where nne 
                  is the number of nodes on eid-th element.
  @param ncne - number of computed nodes on element (only first ncne of nodes is calculated)
  @param nq - number of required mechanical quantities
  @param qt - array of types of required mechanical quantities
   
  @return The function does not return anything but it stores nodal values in parameter nodval.
 
  01/11/2016 by TKr according to elem_mechq_nodval2
*/
void elem_mechq_nodval_comp2(long eid, vector &nodval, long ncne, long nq, nontransquant *qt)
{
  elemtype te;

  te = Mt->give_elem_type (eid);

  switch (te){
  case bar2d:{
    Bar2d->mechq_nodval_comp (eid,nodval,ncne,nq,qt);
    break;
  }
  case barq2d:{
    Barq2d->mechq_nodval_comp (eid,nodval,ncne,nq,qt);
    break;
  }
  case planeelementlq:{
    Pelq->mechq_nodval_comp (eid,nodval,ncne,nq,qt);
    break;
  }
  case planeelementqq:{
    Peqq->mechq_nodval_comp (eid,nodval,ncne,nq,qt);
    break;
  }
  case axisymmlq:{
    Asymlq->mechq_nodval_comp (eid,nodval,ncne,nq,qt);
    break;
  }
  case axisymmqq:{
    Asymqq->mechq_nodval_comp (eid,nodval,ncne,nq,qt);
    break;
  }
  case axisymmcq:{
    Asymcq->mechq_nodval_comp (eid,nodval,ncne,nq,qt);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}


/**
  Function computes global coordinates of integration point (ipp) on
  element (eid).
   
  Function is used in adjacip and in mesh transer values procedures.
   
  @param eid - element id
  @param ipp - integration point pointer
  @param ri,ci - row and column indices
  @param ipcoord - %vector containing coordinates of integration point (output)
   
  @return The function returns coordinates in the parameter ipcoord.

  Created by JK,
*/
void ipcoord (long eid,long ipp,long ri,long ci,vector &ipcoord)
{
  switch (Mt->give_elem_type(eid)){
    case bar2d:
      Bar2d->ipcoord (eid,ipp,ri,ci,ipcoord);
      break;
    case bar3d:
      Bar3d->ipcoord (eid,ipp,ri,ci,ipcoord);
      break;
    case barq2d:
      Barq2d->ipcoord (eid,ipp,ri,ci,ipcoord);
      break;
    case spring_1:
    case spring_2:
    case spring_3:
    case spring_4:
    case spring_5:
    case spring_6:
      Spring->ipcoord(eid, ipp, ipcoord);
      break;
    case axisymmlt:
      Asymlt->ipcoord(eid, ipp, ri, ci, ipcoord);
      break;
    case axisymmqt:
      Asymqt->ipcoord(eid, ipp, ri, ci, ipcoord);
      break;
    case axisymmlq:
      Asymlq->ipcoord(eid, ipp, ri, ci, ipcoord);
      break;
    case axisymmqq: 
      Asymqq->ipcoord(eid, ipp, ri, ci, ipcoord);
      break;
    case axisymmcq: 
      Asymcq->ipcoord(eid, ipp, ri, ci, ipcoord);
      break;
    case axisymmlqintface:{
      Asymlqifc->ipcoord(eid, ipp, ri, ci, ipcoord);
      break;
    }
    case planeelementlt:{
      Pelt->ipcoord (eid,ipp,ri,ci,ipcoord);
      break;
    }
    case planeelementqt:{
      Peqt->ipcoord (eid,ipp,ri,ci,ipcoord);
      break;
    }
    case planeelementrotlt:{
      Perlt->ipcoord (eid,ipp,ri,ci,ipcoord);
      break;
    }
    case planeelementlq:{
      Pelq->ipcoord (eid,ipp,ri,ci,ipcoord);
      break;
    }
    case planeelementqq:{
      Peqq->ipcoord (eid,ipp,ri,ci,ipcoord);
      break;
    }
    case planeelementrotlq:{
      Perlq->ipcoord (eid,ipp,ri,ci,ipcoord);
      break;
    }

    case lineartet:{
      Ltet->ipcoord (eid,ipp,ri,ci,ipcoord);
      break;
    }
    case quadrtet:{
      Qtet->ipcoord (eid,ipp,ri,ci,ipcoord);
      break;
    }
    case linearhex:{
      Lhex->ipcoord (eid,ipp,ri,ci,ipcoord);
      break;
    }
    case quadrhex:{
      Qhex->ipcoord (eid,ipp,ri,ci,ipcoord);
      break;
    }
    case hexintface:{
      Hexifc->ipcoord(eid,ipp, ri, ci, ipcoord);
      break;
    }
	
    default:{
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
    }
  }
}



/**
  Function computes natural coordinates of integration point (ipp) on
  element (eid).
   
  Function is used in mesh transer values procedures.
   
  @param eid - element id
  @param ipp - integration point pointer
  @param ipncoord - %vector containing natural coordinates of integration point (output)
   
  @return The function returns coordinates in the parameter ipcoord.

  Created by JK,
*/
void ipncoord (long eid,long ipp,vector &ipncoord)
{
  switch (Mt->give_elem_type(eid)){
    case bar2d:
      Bar2d->ipncoord (eid,ipp,ipncoord);
      break;
    case bar3d:
      Bar3d->ipncoord (eid,ipp,ipncoord);
      break;
    case barq2d:
      Barq2d->ipncoord (eid,ipp,ipncoord);
      break;
    case barq3d:
      Barq3d->ipncoord (eid,ipp,ipncoord);
      break;
    case axisymmlt:
      Asymlt->ipncoord(eid, ipp, ipncoord);
      break;
    case axisymmqt:
      Asymqt->ipncoord(eid, ipp, ipncoord);
      break;
    case axisymmlq:
      Asymlq->ipncoord(eid, ipp, ipncoord);
      break;
    case axisymmqq: 
      Asymqq->ipncoord(eid, ipp, ipncoord);
      break;
    case axisymmcq: 
      Asymcq->ipncoord(eid, ipp, ipncoord);
      break;
    case axisymmlqintface:
      Asymlqifc->ipncoord(eid, ipp, ipncoord);
      break;
    case planeelementlt:
      Pelt->ipncoord(eid, ipp, ipncoord);
      break;
    case planeelementqt:
      Peqt->ipncoord(eid, ipp, ipncoord);
      break;
    case planeelementrotlt:
      Perlt->ipncoord(eid, ipp, ipncoord);
      break;
    case planeelementlq:
      Pelq->ipncoord(eid, ipp, ipncoord);
      break;
    case planeelementqq:
      Peqq->ipncoord(eid, ipp, ipncoord);
      break;
    case planeelementrotlq:
      Perlq->ipncoord(eid, ipp, ipncoord);
      break;

    case lineartet:
      Ltet->ipncoord(eid,ipp,ipncoord);
      break;
    case quadrtet:
      Qtet->ipncoord(eid, ipp, ipncoord);
      break;
    case linearhex:
      Lhex->ipncoord(eid, ipp, ipncoord);
      break;
    case quadrhex:
      Qhex->ipncoord(eid, ipp, ipncoord);
      break;
    default:
      print_err("unknown element type %d is required", __FILE__, __LINE__, __func__, Mt->give_elem_type(eid));
  }
}



/**
  Function computes strain-displacement (geometric) matrix at point given by natural coordinates on
  element (eid).
   
  Function is used in mesh transer values procedures.
   
  @param eid - element id [in]
  @param xi, eta, zeta - natural coordinates of the given point [in]
  @param gm - resulting geometric %matrix [out]
   
  @return The function returns geometric matrix at point given by natrual coordinates in the argument gm.

  Created by Tomas Koudelka, 22.11.2017
*/
void geom_matrix(long eid, double xi, double eta, double zeta, matrix &gm)
{
  long nne = Mt->give_nne(eid);
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne));
  vector areacoord, l, nx, ny;
  double jac;

  switch (Mt->give_elem_type(eid)){
  case bar2d:{
    //  node cooridnates
    Mt->give_node_coord2d (x,y,eid);
    
    //  geometric matrix
    Bar2d->geom_matrix (gm,x,y,jac);
    break;
  }
    case bar3d:
      Bar3d->geom_matrix(gm,x,y,z,jac);
      break;
    case axisymmlt:
      reallocv(RSTCKVEC(3, areacoord));
      Mt->give_node_coord2d(x, y, eid);
      areacoord[0] = xi;
      areacoord[1] = eta;
      areacoord[2] = 1.0-xi-eta;
      Asymlt->geom_matrix(gm, areacoord, x, y);
      break;
    case axisymmqt:
      Mt->give_node_coord2d(x, y, eid);
      Asymqt->geom_matrix (gm,xi,eta,x,y,jac);
      break;
    case axisymmlq:
      Mt->give_node_coord2d(x, y, eid);
      Asymlq->geom_matrix(gm, x, y, xi, eta, jac);
      break;
    case axisymmqq: 
      Mt->give_node_coord2d(x, y, eid);
      Asymqq->geom_matrix(gm, x, y, xi, eta, jac);
      break;
    case axisymmcq: 
      Mt->give_node_coord2d(x, y, eid);
      Asymcq->geom_matrix(gm, x, y, xi, eta, jac);
      break;
    case axisymmlqintface:
      Mt->give_node_coord2d(x, y, eid);
      Asymlqifc->geom_matrix(gm, x, y, xi, jac);
      break;
    case planeelementlt:
      Mt->give_node_coord2d(x, y, eid);
      Pelt->geom_matrix(gm, x, y);
      break;
    case planeelementqt:
      Mt->give_node_coord2d(x, y, eid);
      Peqt->geom_matrix(gm, x, y, xi, eta, jac);
      break;
    case planeelementrotlt:
      reallocv(RSTCKVEC(3, areacoord));
      Mt->give_node_coord2d(x, y, eid);
      areacoord[0] = xi;
      areacoord[1] = eta;
      areacoord[2] = 1.0-xi-eta;
      Perlt->geom_matrix(gm, areacoord, x, y);
      break;
    case planeelementlq:
      Mt->give_node_coord2d(x, y, eid);
      Pelq->geom_matrix(gm, x, y, xi, eta, jac);
      break;
    case planeelementqq:
      Mt->give_node_coord2d(x, y, eid);
      Peqq->geom_matrix(gm, x, y, xi, eta, jac);
      break;
    case planeelementrotlq:
      reallocv(RSTCKVEC(4, l));
      reallocv(RSTCKVEC(4, nx));
      reallocv(RSTCKVEC(4, ny));
      Perlq->auxdata(x, y, l, nx, ny);
      Perlq->geom_matrix(gm, x, y, xi, eta, l, nx, ny, jac);
      break;

    case lineartet:
      Mt->give_node_coord3d(x, y, z, eid);
      Ltet->geom_matrix(gm, x, y, z, jac);
      break;
    case quadrtet:
      Mt->give_node_coord3d(x, y, z, eid);
      Qtet->geom_matrix(gm, x, y, z, xi, eta, zeta, jac);
      break;
    case linearhex:
      Mt->give_node_coord3d(x, y, z, eid);
      Lhex->geom_matrix(gm, x, y, z, xi, eta, zeta, jac);
      break;
    case quadrhex:
      Mt->give_node_coord3d(x, y, z, eid);
      Qhex->geom_matrix(gm, x, y, z, xi, eta, zeta, jac);
      break;
    default:
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function computes strains at auxiliary integration points.
  The resulting strains are stored in particular auxiliary integration points in mechmat 
  array Mm->aip. The function is intended for the transfer of values among meshes
  in coupled problems especially.

  @param lcid[in] - load case id
  @param n[in]    - the number of required auxiliary points in the mapping array ipm
  @param ipm[in]  - integration point mapping array, 
                    ipm[i].ipp < 0 => auxiliary integration point must be used
                    ipm[i].ipp >= 0 => direct mapping to the regulal integration point, 
                                       no computation of strains is performed, the strains are assumed 
                                       to be computed at the main solution procedure of the problem

  @return The function does not return anything but it stores computed strains at auxiliary 
          integration points at the array Mm->aip.

  Created by Tomas Koudelka, 22.11.2017
*/
void compute_aipstrains(long lcid, long n, ipmap *ipm)
{
  long i, j, nne, ndofe, eid, transf, ncompstr, ipp, app;
  vector eps, aux, r;
  ivector nodes;
  matrix tmat, gm;
  intpoints *tmp_ip;
  long      *tmp_elip;
  double   **tmp_nonmechq;
  double   **tmp_ic;
  double   **tmp_eigstrains;
  double   **tmp_eigstresses;
  double   **tmp_tempstrains;
  
  for (i=0; i<n; i++)
  {
    if (ipm[i].ipp >= 0)  // direct mapping to the regular integration point => no computation
      continue;

    eid = ipm[i].eid;
    ndofe = Mt->give_ndofe(eid);
    reallocv(RSTCKVEC(ndofe, r));
    eldispl(lcid, eid, r.a);
  
    //  transformation of displacement vector to the local nodal coord. systems
    nne = Mt->give_nne(eid);
    reallocv(RSTCKIVEC(nne, nodes));
    Mt->give_elemnodes(eid, nodes);
    transf = Mt->locsystems(nodes);
    if (transf>0)
    {
      reallocv(RSTCKVEC(ndofe, aux));
      reallocm(RSTCKMAT(ndofe, ndofe, tmat));
      elem_transf_matrix (eid, nodes, tmat);
      lgvectortransf (aux, r, tmat);
      copyv (aux,r);
    }

    // compute strain-displacement (geometric) matrix
    ncompstr = Mt->give_tncomp(eid); // number of strain components used on element level
    reallocm(RSTCKMAT(ncompstr, ndofe, gm));
    geom_matrix(eid, ipm[i].xi, ipm[i].eta, ipm[i].zeta, gm);

    // compute strains, eps = B.r    
    app = ipm[i].app;
    // the length of Mm->aip[app].strain array, i.e. strain components used on material level,
    // is/must be greater or equaled to the number of strain components used on element level
    makerefv(eps, Mm->aip[app].strain, ncompstr);
    mxv(gm, r, eps);
  }

  // compute actual values of eigenstrains and eigenstresses
  if (Mp->eigstrains)
    compute_aipeigstr(Mp->time, n, ipm);

  // compute temperature strains
  Mb->aip_temperstrains(lcid, n, ipm);

  // swap regular integration and point auxiliary integration point arrays
  // from now, all functions will work with AUXILIARY integration points 
  tmp_ip          = Mm->ip;
  tmp_elip        = Mm->elip;
  tmp_nonmechq    = Mm->nonmechq;
  tmp_ic          = Mm->ic;
  tmp_eigstrains  = Mm->eigstrains;
  tmp_eigstresses = Mm->eigstresses;
  tmp_tempstrains = Mm->tempstrains;

  Mm->tnip         = Mm->tnaip;
  Mm->ip           = Mm->aip;
  Mm->elip         = Mm->elaip;
  Mm->nonmechq     = Mm->aip_nonmechq;
  Mm->ic           = Mm->aip_ic;
  Mm->eigstrains   = Mm->aip_eigstrains;
  Mm->eigstresses  = Mm->aip_eigstresses;
  Mm->tempstrains  = Mm->aip_tempstrains;

  //  function modifies strains
  //  eigenstrains are subtracted from the total strains
  //  it works only for the first load case
  Mm->aip_eigstrmod();
  
  if (Mb->give_num_mstress_comp(lcid)){
    //  homogenization, stress approach
    //  element functions compute fluctuation strains only
    //  the following function adds macroscopic strains computed in the system of equations
    
    double *lhs = Lsrs->give_lhs(lcid);
    ncompstr = Mt->max_ncompstr;
    vector macrostrains(ncompstr);
    strastre *mstrastre = Mb->give_mstrastre(lcid);
    j = Ndofm - Mb->give_num_mstress_comp(lcid);
    for (i=0;i<ncompstr;i++){
      if (mstrastre[i] == stress){
        macrostrains[i] = lhs[j];
        j++;
      }
      else
        macrostrains[i] = 0.0;
    }
    Mm->add_macro_strains(lcid, macrostrains);
  }
  if (Mb->give_num_mstrain_comp(lcid)){
    //  homogenization, strain approach
    //  element functions compute fluctuation strains only
    //  the following function adds macroscopic strains computed in the system of equations
    ncompstr = Mt->max_ncompstr;
    vector macrostrains(ASTCKVEC(ncompstr));
    switch(Mp->tprob)
    {
      case linear_statics:
        copyv(Mb->lc[lcid].mstrain, macrostrains);
        Mm->add_macro_strains(lcid, macrostrains);
        break;
      case mat_nonlinear_statics:
        addmultv(Mb->lc[lcid+1].mstrain, Mb->lc[lcid].mstrain, Mp->lambda, macrostrains.a, ncompstr);
        Mm->add_macro_strains(lcid, macrostrains);
        break;
      case mech_timedependent_prob:
        for(i=0; i<ncompstr; i++)
          macrostrains(i) = Mb->dlc[lcid].get_macrostra(Mp->time, i);
        Mm->add_macro_strains(lcid, macrostrains);
        break;
      default:
        print_err("unknown type of problem is required", __FILE__, __LINE__, __func__);
    }
  }
  if (lcid != 0)
  {
    for (ipp=0; ipp<Mm->tnip; ipp++)
    {
      ncompstr = Mm->ip[ipp].ncompstr;
      if (Gtm->leso[Mm->elip[ipp]]==1){
        memcpy(Mm->ip[ipp].strain, Mm->ip[ipp].strain+ncompstr*lcid, sizeof(Mm->ip[ipp].strain)*ncompstr);
      }
    }
  }

  // swap regular integration and point auxiliary integration point arrays
  // from now, all functions will work with REGULAR integration points 
  Mm->tnip         = Mm->tnrip;
  Mm->ip           = tmp_ip;
  Mm->elip         = tmp_elip;
  Mm->nonmechq     = tmp_nonmechq;
  Mm->ic           = tmp_ic;
  Mm->eigstrains   = tmp_eigstrains;
  Mm->eigstresses  = tmp_eigstresses;
  Mm->tempstrains  = tmp_tempstrains;
}



/**
  The function computes stresses together with eventual state variables and other values at auxiliary integration 
  points on elements. The resulting stresses and state variables are stored in particular auxiliary integration 
  points in mechmat array Mm->aip. The function is intended for the transfer of values among meshes
  in coupled problems especially.

  @param lcid - load case id

  @return The function does not return anything but it stores computed stresses at auxiliary 
          integration points at the array Mm->aip.

  Created by Tomas Koudelka, 22.11.2017
*/
void compute_aipstresses(long /*lcid*/)
{
  long err, app, eid;
  intpoints *tmp_ip;
  long      *tmp_elip;
  double   **tmp_nonmechq;
  double   **tmp_ic;
  double   **tmp_eigstrains;
  double   **tmp_eigstresses;
  double   **tmp_tempstrains;

  // swap regular integration point and auxiliary integration point arrays
  // from now, all functions will work with AUXILIARY integration points 
  tmp_ip          = Mm->ip;
  tmp_elip        = Mm->elip;
  tmp_nonmechq    = Mm->nonmechq;
  tmp_ic          = Mm->ic;
  tmp_eigstrains  = Mm->eigstrains;
  tmp_eigstresses = Mm->eigstresses;
  tmp_tempstrains = Mm->tempstrains;

  Mm->tnip         = Mm->tnaip;
  Mm->ip           = Mm->aip;
  Mm->elip         = Mm->elaip;
  Mm->nonmechq     = Mm->aip_nonmechq;
  Mm->ic           = Mm->aip_ic;
  Mm->eigstrains   = Mm->aip_eigstrains;
  Mm->eigstresses  = Mm->aip_eigstresses;
  Mm->tempstrains  = Mm->aip_tempstrains;

#ifdef INC_OPENMP
  #pragma omp parallel num_threads(Numth)
  {
   #pragma omp for
#endif

    // just call computenlstresses procedure on all auxiliary integration points
    // and check for math errors
    for (app=0; app<Mm->tnaip; app++){
      eid = Mm->elip[app];
      if (Gtm->leso[eid] == 1){
        Mm->computenlstresses(app,Mm->ip[app]);
        //err = check_math_errel(Mm->elip[app]);
        //if (err)
	//abort();/debug??!!
      }
    }

#ifdef INC_OPENMP
  }
#endif
  // swap regular integration point and auxiliary integration point arrays
  // from now, all functions will work with REGULAR integration points 
  Mm->tnip         = Mm->tnrip;
  Mm->ip           = tmp_ip;
  Mm->elip         = tmp_elip;
  Mm->nonmechq     = tmp_nonmechq;
  Mm->ic           = tmp_ic;
  Mm->eigstrains   = tmp_eigstrains;
  Mm->eigstresses  = tmp_eigstresses;
  Mm->tempstrains  = tmp_tempstrains;
}



/**
  The function computes actual values of eigenstrains and eigenstresses at auxiliary integration points.
  The resulting strains are stored in particular auxiliary integration points in mechmat 
  array Mm->aip. The function is intended for the transfer of values amon meshes
  in coupled problems especially.

  @param time[in] - actual time
  @param n[in]    - the number of required auxiliary points in the mapping array ipm
  @param ipm[in]   - integration point mapping array, 
                     ipm[i].ipp < 0 => auxiliary integration point must be used
                     ipm[i].ipp >= 0 => direct mapping to the regulal integration point, 
                                        no computation of strains is performed, the strains are assumed 
                                        to be computed at the main solution procedure of the problem

  @return The function does not return anything but it stores computed eigentstrains at auxiliary 
          integration points at the arrays Mm->aip_eigstrains and Mm->aip_eigstresses.

  Created by Tomas Koudelka, 15.06.2018
*/
void compute_aipeigstr(double time, long n, ipmap *ipm)
{
  long i, j, ipp, app, ncomp, gfid;
  long eid, nne;
  vector x, y, z;
  vector natcoord(ASTCKVEC(3));
  vector coord(ASTCKVEC(4)); // three spatial coordinates  + time of the aux. int. point
  const char *namev[4];
  namev[0] = "x";
  namev[1] = "y";
  namev[2] = "z";
  namev[3] = "t";
  
  //  loop over the number of all auxiliary integration points
  for (i=0;i<n;i++)
  {
    app = ipm[i].app;
    if (app >=0)
    {
      //  the true number of strain components in the integration point
      ncomp = Mm->aip[app].ncompstr;
      // element id of the aux. int. point
      eid = Mm->elaip[app];
      // give the first regular int. point on the element due to the index in eigstrid array
      ipp = Mt->elements[eid].ipp[0][0];
      //  loop over the number of eigenstrain components
      for (j=0;j<ncomp;j++)
      {
        // id of the general function is obtained from the first regular int. point of the element
        // which the auxiliary int. point is connected with
        gfid = Mm->eigstrid[ipp][j];
        // prepare coordinates of the auxiliary int. points
        nne = Mt->give_nne(eid);
        reallocv(RSTCKVEC(nne, x));
        reallocv(RSTCKVEC(nne, y));
        reallocv(RSTCKVEC(nne, z));
        // nodal coordinates on eid-th element
        Mt->give_node_coord3d(x, y, z, eid);
        // natural coordinates of the aux. int. point
        natcoord(0) = ipm[app].xi;
        natcoord(1) = ipm[app].eta;
        natcoord(2) = ipm[app].zeta;
        // compute/interpolate global coordinates of the aux. int. point 
        coord(0) = interpolate(eid, x, natcoord);
        coord(1) = interpolate(eid, y, natcoord);
        coord(2) = interpolate(eid, z, natcoord);
        // prepare time coordinate
        coord(3) = time;
        if (Mp->eigstrains < 3)
          Mm->aip_eigstrains[app][j] = Mb->eigstrfun[gfid].getval(coord, namev);
        if (Mp->eigstrains > 3)
          Mm->aip_eigstresses[app][j] = Mb->eigstrfun[gfid].getval(coord, namev);
      }
    }
  }
}



/**
   function selects components from the macro level to meso level for homogenization
   
   @param eid - element id
   @param counter - actual position in the array buff
   @param buff - array containing components
   
   TKr 22/07/2022 accroding to JK
*/
void higher_to_lower_level_elemm (long eid,long */*counter*/,double */*buff*/)
{
  elemtype te;
  
  //  type of element
  te=Mt->give_elem_type (eid);
  
  switch (te){
    case planeelementqt:
      //Peqt->higher_to_lower_level (eid,counter,buff);
      break;
    case planeelementrotlt:
      //Perlt->higher_to_lower_level (eid,counter,buff);
      break;
    case planeelementlq:
      //Pelq->higher_to_lower_level (eid,counter,buff);
      break;
    case planeelementqq:
      //Peqq->higher_to_lower_level (eid,counter,buff);
      break;
    case axisymmlq:
      //Asymlq->higher_to_lower_level (eid,counter,buff);
      break;
    case axisymmqq: 
      //Asymqq->higher_to_lower_level (eid,counter,buff);
      break;
    case lineartet:
      //Ltet->higher_to_lower_level (eid,counter,buff);
      break;
    case quadrtet:
      //Qtet->higher_to_lower_level (eid,counter,buff);
      break;
    case linearhex:
      //Lhex->higher_to_lower_level (eid,counter,buff);
      break;
    case quadrhex:
      //Qhex->higher_to_lower_level (eid,counter,buff);
      break;
    default:
      print_err("unknown element type %d is required", __FILE__, __LINE__, __func__, Mt->give_elem_type(eid));
  }
}

