#include "homogmech.h"
#include "homog.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "mechbclc.h"
#include "element.h"
#include "intpoints.h"
#include "globmat.h"
#include "gtopology.h"
#include "math.h"
#include "solverm.h"
#include "lssolver.h"
#include "mechbclc.h"
#include "loadcase.h"
#include "elemswitch.h"

/**
   function solves homogenization for mechanical problems on PUC

   TKr, 07/01/2015
*/
void mechanical_homogenization ()
{
  long i,j,k,ii,jj,hncomp,tnipel,ipp;

  //  initiation of mechanical material models
  //Mm->initmaterialmodels();
  //  nodes - integration points interpolation
  //approximation_puc (); //if it is needed
  
  
  //solving of fluctuation parts
  // zatim pouze pro jeden zatezovaci stav a elasticky material:
  
  //tady doplnit cyklus pres vsechny jednotkove stavy podle typu homog. ulohy:
  Mp->stresscomp = 1; //for macrostresses computation
  Mp->straincomp = 1; //for macrostrains computation
  switch (Mp->hstrastre){
  case planestress:{
    hncomp = 3;
    break;
  }
  case planestrain:{
    hncomp = 3;
    break;
  }
  case axisymm:{
    hncomp = 4;
    break;
  }
  case spacestress:{
    hncomp = 6;
    break;
  }
  default:{
    print_err("unknown type of strain-stress state for mechanical homogeniaztion is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }//  end of the switch (strastre){

  //matrix allocation
  matrix d(ASTCKMAT(hncomp,hncomp));
  vector mstress(ASTCKVEC(hncomp));

  if(Mp->homog > 6){//options 7 or 8 for testing
    ii = Mp->homogdir-1; //unifom load direction
    
    //macro strain-stress loading:
    if (Mp->homog == 7)
      {
	for(jj=0;jj<hncomp;jj++)//null mstress          
	  Mb->lc[0].mstress[jj] = 0.0;
	
	Mb->lc[0].mstress[ii] = 1.0;//unit mstress for ii components
      }
    if (Mp->homog == 8)
      {
	for(jj=0;jj<hncomp;jj++)//null mstrain
	  Mb->lc[0].mstrain[jj] = 0.0;
	
	Mb->lc[0].mstrain[ii] = 1.0;//unit mstrain for ii components
      }

    switch (Mp->tprob){
    case linear_statics:{
      solve_linear_statics ();
      break;
    }
    default:{
      print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
    }
    }
    
    switch (Mp->homog){
    case 7:{ // macro-stress approach
      //compliance matrix computation from averaged strains
      for(i=0;i<hncomp;i++)
	d[i][ii] = macrostrains(0,i);
      break;
    }
    case 8:{ // macro-strain approach
      //stiffness matrix computation from averaged stresses (macrostress)
      for(i=0;i<hncomp;i++)
	d[i][ii] = macrostresses(0,i);
      break;
    }
    default:{
      print_err("unknown type of homogenization is required", __FILE__, __LINE__, __func__);
      abort();
    }
    }
  }
  else{
    
    //loop over strain-stress components; only for one (first) loadcase
    for(ii=0;ii<hncomp;ii++){
      
      //macro strain-stress loading:
      if (Mp->homog == 5)
	{
	  for(jj=0;jj<hncomp;jj++)//null mstress
	    Mb->lc[0].mstress[jj] = 0.0;
	  
	  Mb->lc[0].mstress[ii] = 1.0;//unit mstress at ii-th component
	}
      if (Mp->homog == 6)
	{
	  for(jj=0;jj<hncomp;jj++)//null mstrain
	    Mb->lc[0].mstrain[jj] = 0.0;
	  
	  Mb->lc[0].mstrain[ii] = 1.0;//unit mstrain at ii-th component
	}
      
      switch (Mp->tprob){
      case linear_statics:{
	solve_linear_statics ();
	break;
      }
      default:{
	print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
      }
      }
      
      switch (Mp->homog){
      case 5:{ // macro-stress approach
	//compliance matrix computation from averaged strains
	for(i=0;i<hncomp;i++)
	  d[i][ii] = macrostrains(0,i);
	  
	fprintf(Out,"\n Average stress control (%ld):\n",ii+1);

	//ip volumes needed
	if(Mm->ipv == NULL)
	  ipvolume ();
	
	for (jj=0; jj<Mm->max_ncompstre; jj++)
	  mstress[jj] = 0.0;
	
	for (k=0;k<Mt->ne;k++){
	  if (Gtm->leso[k]==1){
	    if(Mt->give_elem_type (k) == planeelementlq)
	      tnipel = 4;//only the first block for quadrilateral elelemnts 
	    else
	      tnipel = Mt->give_tnip(k);
	    ipp = Mt->elements[k].ipp[0][0];
	    for(j=0;j<tnipel;j++){
	      
	      for (jj=0; jj<Mm->max_ncompstre; jj++)
		mstress[jj]+=Mm->ip[ipp].stress[jj]*Mm->ipv[ipp];
	      
	      ipp++;
	    }
	  }
	}
	
	for (j=0; j<Mm->max_ncompstre; j++)
	  mstress[j]/=Mt->domvol;
	
	for (jj=0; jj<Mm->max_ncompstre; jj++)
	  fprintf(Out,"\n Macrostress(%ld) = %e\n",jj+1,mstress[jj]);

	break;
      }
      case 6:{ // macro-strain approach
	//stiffness matrix computation from averaged stresses (macrostress)
	for(i=0;i<hncomp;i++)
	  d[i][ii] = macrostresses(0,i);
	
	fprintf(Out,"\n Average strain control (%ld):\n",ii+1);

	for (jj=0; jj<Mm->max_ncompstre; jj++)
	  mstress[jj] = 0.0;

	for (k=0;k<Mt->ne;k++){
	  if (Gtm->leso[k]==1){
	    if(Mt->give_elem_type (k) == planeelementlq)
	      tnipel = 4;//only the first block for quadrilateral elelemnts 
	    else
	      tnipel = Mt->give_tnip(k);
	    ipp = Mt->elements[k].ipp[0][0];
	    for(j=0;j<tnipel;j++){
	      
	      for (jj=0; jj<Mm->max_ncompstre; jj++)
		mstress[jj]+=Mm->ip[ipp].strain[jj]*Mm->ipv[ipp];
	      //mstress[jj]+=Mm->ipv[ipp];
	      
	      ipp++;
	    }
	  }
	}
	
	for (j=0; j<Mm->max_ncompstre; j++)
	  mstress[j]/=Mt->domvol;
	
	for (jj=0; jj<Mm->max_ncompstre; jj++)
	  fprintf(Out,"\n Macrostrain(%ld) = %e\n",jj+1,mstress[jj]);

	break;
      }
      default:{
	print_err("unknown type of homogenization is required", __FILE__, __LINE__, __func__);
	abort();
      }
      }
    }//end of loop over strain-stress components
  }

  if (Mespr != 0){//printing of matrix into .log file
    switch (Mp->homog){
    case 5:
    case 7:{ // macro-stress approach
      fprintf(Out,"\n\n Resulting Compliance matrix C\n");
      break;
    }
    case 6:
    case 8:{ // macro-strain approach
      fprintf(Out,"\n\n Resulting Stiffness matrix D\n");
      break;
    }
    default:{
      print_err("unknown type of homogenization is required", __FILE__, __LINE__, __func__);
      abort();
    }
    }
    
    for (i=0;i<hncomp;i++){
      for (j=0;j<hncomp;j++)
	fprintf(Out,"%e  ",d[i][j]);
      fprintf(Out,"\n");
    }
    
    //printing for general anisotropic elastic material:
    switch (Mp->hstrastre){
    case planestress:
    case planestrain:{
      fprintf(Out,"\n\n Print for elastgmat2d:\n");
      fprintf(Out,"%le %le %le\n",d[0][0],d[0][1],d[0][2]);
      fprintf(Out,"%le %le\n",d[1][1],d[1][2]);
      fprintf(Out,"%le\n",d[2][2]);
      break;
    }
    case spacestress:{
      fprintf(Out,"\n\n Print for elastgmat3d:\n");
      fprintf(Out,"%le %le %le %le %le %le\n",d[0][0],d[0][1],d[0][2],d[0][3],d[0][4],d[0][5]);
      fprintf(Out,"%le %le %le %le %le\n",d[1][1],d[1][2],d[1][3],d[1][4],d[1][5]);
      fprintf(Out,"%le %le %le %le\n",d[2][2],d[2][3],d[2][4],d[2][5]);
      fprintf(Out,"%le %le %le\n",d[3][3],d[3][4],d[3][5]);
      fprintf(Out,"%le %le\n",d[4][4],d[4][5]);
      fprintf(Out,"%le\n",d[5][5]);
      break;
    }
    default:{
      print_err("unknown type of strain-stress state for mechanical homogeniaztion is required", __FILE__, __LINE__, __func__);
    }
    }//  end of the switch (strastre){
  
    
    // volume fraction printing
    double vol=0.0,*volf;
    strastrestate ssst;

    for (i=0; i <Mt->ne; i++){
      ssst=Mt->give_ssst (i,0);
      if(ssst == spacestress)//for 3D elements
	{
	  vol += Mt->give_volume(i);
	}
      else
	{
	  vol += Mt->give_area(i);
	}
    }
    fprintf(Out,"Domain volume = %le\n",vol);


    volf = new double [Mm->numtype[0]];

    for(ii=0;ii<Mm->numtype[0];ii++){
      volf[ii] = 0.0;
      for (i=0; i <Mt->ne; i++){
	if(ii == Mt->elements[i].idm[0]){
	  ssst=Mt->give_ssst (i,0);
	  if(ssst == spacestress)//for 3D elements
	    {
	      volf[ii] += Mt->give_volume(i);
	    }
	  else
	    {
	      volf[ii] += Mt->give_area(i);
	  }
	}
      }
      fprintf(Out,"Volume fraction (%ld)= %lf\n",ii+1,volf[ii]/vol);
    }

    delete [] volf;
  }
}

/**
   function solves homogenization for mechanical problems on PUC and returns overall (homogenized) stiffness matrix 

   @param unkn_r   - displacements and strains and stresses in integration point 
                     (in master node on PUC) stored as an array (one dimensional)
   @param matrix_k - homogenized stiffness matrix of PUC stored as an array (one dimensional)

   TKr, 12/07/2022
*/
void paral_mechanical_homogenization (double */*unkn_r*/,double *matrix_k)
{
  long i,j,k,ii,jj,hncomp,tnipel,ipp;
  
  //  initiation of mechanical material models
  //Mm->initmaterialmodels();
  //  nodes - integration points interpolation
  //approximation_puc (); //if it is needed
  
  //solving of fluctuation parts
  //only for the first loadcase and elastic material:
  
  Mp->stresscomp = 1; //for macrostresses computation
  Mp->straincomp = 1; //for macrostrains computation
  switch (Mp->hstrastre){
  case planestress:{
    hncomp = 3;
    break;
  }
  case planestrain:{
    hncomp = 3;
    break;
  }
  case axisymm:{
    hncomp = 4;
    break;
  }
  case spacestress:{
    hncomp = 6;
    break;
  }
  default:{
    print_err("unknown type of strain-stress state for mechanical homogeniaztion is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }//  end of the switch (strastre){

  //matrix allocation
  matrix d(ASTCKMAT(hncomp,hncomp));
  vector mstress(ASTCKVEC(hncomp));

  //loop over strain-stress components; only for one (first) loadcase
  for(ii=0;ii<hncomp;ii++){
    
    //only macro strain loading:
    if (Mp->homog == 6)
      {
	for(jj=0;jj<hncomp;jj++)//null mstrain
	  Mb->lc[0].mstrain[jj] = 0.0;
	
	Mb->lc[0].mstrain[ii] = 1.0;//unit mstrain at ii-th component
      }
    else{
      print_err("only macro strain loading is required",__FILE__,__LINE__,__func__);
      exit(0);
    }
    
    switch (Mp->tprob){
    case linear_statics:{
      solve_linear_statics ();
      break;
    }
    default:{
      print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
    }
    }
    
    switch (Mp->homog){
    case 6:{ // macro-strain approach
      //stiffness matrix computation from averaged stresses (macrostress)
      for(i=0;i<hncomp;i++){
	d[i][ii] = macrostresses(0,i);
	//fprintf (Out,"macrostresses[%ld][%ld] %16.12le\n",i,ii,d[i][ii]);
      }

      //fprintf(Out,"\n Average strain control (%ld):\n",ii+1);
      
      for (jj=0; jj<Mm->max_ncompstre; jj++)
	mstress[jj] = 0.0;
      
      for (k=0;k<Mt->ne;k++){
	if (Gtm->leso[k]==1){
	  if(Mt->give_elem_type (k) == planeelementlq)
	    tnipel = 4;//only the first block for quadrilateral elelemnts 
	  else
	    tnipel = Mt->give_tnip(k);
	  ipp = Mt->elements[k].ipp[0][0];
	  for(j=0;j<tnipel;j++){
	    
	    for (jj=0; jj<Mm->max_ncompstre; jj++)
	      mstress[jj]+=Mm->ip[ipp].strain[jj]*Mm->ipv[ipp];
	    //mstress[jj]+=Mm->ipv[ipp];
	    
	    ipp++;
	  }
	}
      }
      
      for (j=0; j<Mm->max_ncompstre; j++)
	mstress[j]/=Mt->domvol;

      //fprintf (Out,"Mt->domvol %16.12le\n",Mt->domvol);

      //for (jj=0; jj<Mm->max_ncompstre; jj++)
      //	fprintf(Out,"\n Macrostrain(%ld) = %e\n",jj+1,mstress[jj]);
      
      break;
    }
    default:{
      print_err("unknown type of homogenization is required", __FILE__, __LINE__, __func__);
      abort();
    }
    }
  }//end of loop over strain-stress components
  
  
  //storing of overall properties
  //Stiffness matrix D
  k = 0;
  for (i=0;i<hncomp;i++){
    for (j=0;j<hncomp;j++){
      matrix_k[k] = d[i][j];
      k++;
    }
  }

  ////////////////////////////////////////////////////////////
  // this printing below is only for debug:
  /* fprintf(Out,"\n\n vysledna matice D^M\n");
     for (i=0;i<hncomp;i++){
     for (j=0;j<hncomp;j++)
     fprintf(Out,"%e  \n",d[i][j]);
     fprintf(Out,"\n");
     }
  */
  /*
    if (Mespr != 0){//printing of matrix into .log file
    switch (Mp->homog){
    case 6:{// macro-strain approach
    fprintf(Out,"\n\n Resulting Stiffness matrix D\n");
    break;
    }
    default:{
    print_err("unknown type of homogenization is required", __FILE__, __LINE__, __func__);
    abort();
    }
    }
    
    for (i=0;i<hncomp;i++){
    for (j=0;j<hncomp;j++)
    fprintf(Out,"%e  ",d[i][j]);
    fprintf(Out,"\n");
    }
    
    //printing for general anisotropic elastic material:
    switch (Mp->hstrastre){
    case planestress:
    case planestrain:{
    fprintf(Out,"\n\n Print for elastgmat2d:\n");
    fprintf(Out,"%le %le %le\n",d[0][0],d[0][1],d[0][2]);
    fprintf(Out,"%le %le\n",d[1][1],d[1][2]);
    fprintf(Out,"%le\n",d[2][2]);
    break;
    }
    case spacestress:{
    fprintf(Out,"\n\n Print for elastgmat3d:\n");
    fprintf(Out,"%le %le %le %le %le %le\n",d[0][0],d[0][1],d[0][2],d[0][3],d[0][4],d[0][5]);
    fprintf(Out,"%le %le %le %le %le\n",d[1][1],d[1][2],d[1][3],d[1][4],d[1][5]);
    fprintf(Out,"%le %le %le %le\n",d[2][2],d[2][3],d[2][4],d[2][5]);
    fprintf(Out,"%le %le %le\n",d[3][3],d[3][4],d[3][5]);
    fprintf(Out,"%le %le\n",d[4][4],d[4][5]);
    fprintf(Out,"%le\n",d[5][5]);
    break;
    }
    default:{
    print_err("unknown type of strain-stress state for mechanical homogeniaztion is required", __FILE__, __LINE__, __func__);
    }
    }//  end of the switch (strastre){
    }
  */
  //////////////////////////////////////////////////////////
}
