#include "loadelt.h"
#include "globalt.h"
#include "aliast.h"
#include "globmatt.h"
#include <string.h>

loadelt::loadelt()
{
  //  loaded element id (number of loaded element)
  eid=-1;
  //  number of boundary objects (end nodes, edges, surfaces)
  nbo=0;
  //  number of nodes on boundary object
  nnbo=0;
  
  //  boundary conditions (for nongrowing problems)
  bc=NULL;
  //  boundary conditions (for growing problems)
  bcf=NULL;
  
  //  id of nodal values or climatic conditions (for nongrowing problems)
  nvid=NULL;
  //  id of nodal values or climatic conditions (for growing problems)
  nvidf=NULL;
  //  id of nodal values describing transmission coefficients (for nongrowing problems)
  trcid=NULL;
  //  id of nodal values describing transmission coefficients (for growing problems)
  trcidf=NULL;
  //  id of nodal values describing transmission/radiation coefficients (for nongrowing problems)
  trrid=NULL;
  //  id of nodal values describing transmission/radiation coefficients (for growing problems)
  trridf=NULL;
}



/**
   Constructor allocates arrays of id of bnodval objects.
   It is used in the preprocessor.
   
   @param id - element id
   @param numbo - number of boundary objects
   
   @return The constructor does not return anything.
   
   Created by TKo, 10.2010
*/
loadelt::loadelt(long id, long numbo)
{
  eid =id;
  nbo = numbo;
  if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
    //  indicator of boundary conditions
    bcf = new long [nbo];
    //  id of nodal values or climatic conditions
    nvidf = new long [nbo];
    //  id of nodal values describing transmission coefficients
    trcidf = new long [nbo];
    //  id of nodal values describing transmission/radiation coefficients
    trridf = new long [nbo];
  }

  //  indicators of boundary conditions
  bc = new bocontypet [nbo];
  memset(bc, 0, sizeof(*bc)*nbo);
  //  id of nodal values or climatic conditions
  nvid = new long [nbo];
  memset(nvid, 0, sizeof(*nvid)*nbo);
  //  id of nodal values describing transmission coefficients
  trcid = new long [nbo];
  memset(trcid, 0, sizeof(*trcid)*nbo);
  //  id of nodal values describing transmission/radiation coefficients
  trrid = new long [nbo];
  memset(trrid, 0, sizeof(*trrid)*nbo);
}


loadelt::~loadelt()
{
  if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin)
  {
    delete [] bcf;
    delete [] nvidf;
    delete [] trcidf;
    delete [] trridf;
  }
  
  delete [] bc;
  delete [] nvid;
  delete [] trcid;
  delete [] trrid;
}

/**
   function reads boundary conditions on elements
   
   @param in[in]     - pointer to the opened input file
   @param lcid[in]   - load case id
   @param elemid[in] - element id
   @param inbo[in]   - number of boundary objects on the given element, i.e. number of edges, surfaces, ...
   @param innbo[in]  - number of nodes which defines one boundary object on the given element, 
                       i.e. number of nodes on edge, surface, ...
*/
void loadelt::read (XFILE *in,long /*lcid*/,long elemid,long inbo,long innbo)
{
  long i;
  
  //  element id
  eid=elemid;
  //  the number of boundary objects (edges, surfaces)
  nbo=inbo;
  //  the number of nodes on the boundary objects
  nnbo = innbo;
  
  if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
    //  indicator of boundary conditions
    bcf = new long [nbo];
    //  id of nodal values or climatic conditions
    nvidf = new long [nbo];
    //  id of nodal values describing transmission coefficients
    trcidf = new long [nbo];
    //  id of nodal values describing transmission/radiation coefficients
    trridf = new long [nbo];
  }

  //  indicators of boundary conditions
  bc = new bocontypet [nbo];
  //  id of nodal values or climatic conditions
  nvid = new long [nbo];
  //  id of nodal values describing transmission coefficients
  trcid = new long [nbo];
  //  id of nodal values describing transmission/radiation coefficients
  trrid = new long [nbo];
  
  //  loop over boundary objects (end nodes, edges or surfaces)
  for (i=0;i<nbo;i++){
    
    if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){	
      //  type of boundary condition
      xfscanf (in,"%ld",bcf+i);
      bcf[i]--;
      //  id of nodal values
      xfscanf (in,"%ld",nvidf+i);
      nvidf[i]--;
      //  id of transmission coefficients
      xfscanf (in,"%ld",trcidf+i);
      trcidf[i]--;
      //  id of transmission/radiation coefficients
      xfscanf (in,"%ld",trridf+i);
      trridf[i]--;
    }
    else{
      //  type of boundary condition
      xfscanf (in,"%m", &bocontype_kwdset, bc+i);
      
      nvid[i]=-1;
      trcid[i]=-1;
      trrid[i]=-1;
      
      if (bc[i]>=presc_trmiss){
	//  prescribed transmission (Newton boundary conditon)
	
	//  id of exterior prescribed values
	xfscanf (in,"%ld",nvid+i);
	nvid[i]--;
	//  id of transmission coefficients
	xfscanf (in,"%ld",trcid+i);
	trcid[i]--;
	//  id of transmission/radiation coefficients
	xfscanf (in,"%ld",trrid+i);
	trrid[i]--;
      }
      if (bc[i]==presc_flux){
	//  prescribed fluxes

	//  id of prescribed fluxes (nodal values)
	xfscanf (in,"%ld",nvid+i);
	nvid[i]--;
	
      }
      if (bc[i]==det_climcond){
	// prescribed detail climatic conditions
	
	//  id of detail climatic conditions
	xfscanf (in,"%ld",nvid+i);
	nvid[i]--;
      }
      if (bc[i]==gen_climcond){
	// prescribed detail climatic conditions
	
	//  id of detail climatic conditions
	xfscanf (in,"%ld",nvid+i);
	nvid[i]--;
      }
    }
  }
}

/**
   function reads boundary conditions on elements needed for flux computation
   
   @param in - pointer to input file
   @param lcid - load case id
   
   TKr 23/02/2021 according to JK
*/
void loadelt::read_comp (XFILE *in,long /*lcid*/,long elemid,long inbo,long innbo)
{
  long i;
  
  //  element id
  eid=elemid;
  //  the number of boundary objects (edges, surfaces)
  nbo=inbo;
  //  the number of nodes on the boundary objects
  nnbo = innbo;
  
  //  indicators of boundary conditions
  bc = new bocontypet [nbo];
  //  id of nodal values or climatic conditions
  nvid = new long [nbo];
  
  //  loop over boundary objects (end nodes, edges or surfaces)
  for (i=0;i<nbo;i++){
    
    //  type of boundary condition
    xfscanf (in,"%m", &bocontype_kwdset, bc+i);
    
    nvid[i]=-1;
    
    if (bc[i]==presc_comp){
      //  prescribed fluxes
      
      //  id of prescribed fluxes (nodal values)
      xfscanf (in,"%ld",nvid+i);
      nvid[i]--;
      
    }
  }
}

/**
   determination of auxiliary indicators on elements
   they serve for efficient computation of fluxes
   only elements with appropriate transi indicators are performed
   
*/
void loadelt::bc_indicators (long &transi)
{
  long i;
  
  for (i=0;i<nbo;i++){
    if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
    }
    else{
      if (bc[i]==presc_flux || bc[i]==det_climcond){
	if (transi<2){
	  //  the element has not been modified yet
	  transi=2;
	}
	if (transi==2){
	  //  prescribed flux was determined earlier
	}
	if (transi==3){
	  //  transmission was detected earlier
	  transi=4;
	}
	if (transi==4){
	  //  prescribed flux and transmission were detected earlier
	}
      }
      
      if (bc[i]==gen_climcond){
	if (transi<4){
	  transi=4;
	}
	if (transi==4){
	  //  the element has been modified earlier
	}
      }
      
      if (bc[i]>=presc_trmiss){
	if (transi<2){
	  //  the element has not been modified yet
	  transi=3;
	}
	if (transi==2){
	  //  prescribed flux was determined earlier
	  transi=4;	
	}
	if (transi==3){
	  //  transmission was detected earlier
	}
	if (transi==4){
	  //  prescribed flux and transmission were detected earlier
	}
      }
    }
  }
}

/**
   function updates auxiliary indicators on elements for growing structures
   they serve for efficient computation of fluxes
   only elements with appropriate transi indicators are performed
   
   TKr 23/02/2021
*/
void loadelt::updatevalues (long /*lcid*/, long &transi)
{
  long i,j;
  
  for (i=0;i<nbo;i++){

    j = bcf[i];
    bc[i] = (bocontypet)Gtt->gf[j].getval_long (Tp->time);
    j = nvidf[i];
    nvid[i] = Gtt->gf[j].getval_long(Tp->time);
    nvid[i]--;
    j = trcidf[i];
    trcid[i] = Gtt->gf[j].getval_long(Tp->time);
    trcid[i]--;
    j = trridf[i];
    trrid[i] = Gtt->gf[j].getval_long(Tp->time);
    trrid[i]--;
    
    if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin)
      {
	if (bc[i]==presc_flux || bc[i]==det_climcond){
	  if (transi<2){
	    //  the element has not been modified yet
	    transi=2;
	  }
	  if (transi==2){
	    //  prescribed flux was determined earlier
	  }
	  if (transi==3){
	    //  transmission was detected earlier
	    transi=4;
	  }
	  if (transi==4){
	    //  prescribed flux and transmission were detected earlier
	  }
	}
	
	if (bc[i]==gen_climcond){
	  if (transi<4){
	    transi=4;
	  }
	  if (transi==4){
	    //  the element has been modified earlier
	  }
	}
	
	if (bc[i]>=presc_trmiss){
	  if (transi<2){
	    //  the element has not been modified yet
	    transi=3;
	  }
	  if (transi==2){
	    //  prescribed flux was determined earlier
	    transi=4;	
	  }
	  if (transi==3){
	    //  transmission was detected earlier
	  }
	  if (transi==4){
	    //  prescribed flux and transmission were detected earlier
	  }
	}
      }
  }
}

/**
   function prints boundary conditions on elements
   
   @param out - pointer to output file
   @param lcid - load case id
*/
void loadelt::print (FILE *out,long /*lcid*/)
{
  long i;
  
  //  element id
  fprintf (out,"\n\n %ld ",eid+1);
  
  //  loop over boundary objects (end nodes, edges or surfaces)
  for (i=0;i<nbo;i++){
    if (Tp->tprob == growing_np_problem || Tp->tprob == growing_np_problem_nonlin){
      //  type of boundary condition
      fprintf (out," %ld",bcf[i]+1);
      //  id of nodal values
      fprintf (out," %ld",nvidf[i]+1);
      //  id of transmission coefficients
      fprintf (out," %ld",trcidf[i]+1);
      //  id of transmission/radiation coefficients
      fprintf (out," %ld",trridf[i]+1);
    }
    else{
      fprintf (out,"\n %d", int(bc[i]));

      if (bc[i]==presc_flux){
	//  prescribed flux (Neumann boundary conditon)
	//  id of exterior prescribed values
	fprintf (out," %ld",nvid[i]+1);
      }

      if (bc[i]>=presc_trmiss){
	//  prescribed transmission (Newton boundary conditon)
	
	//  id of exterior prescribed values
	fprintf (out," %ld",nvid[i]+1);
	//  id of transmission coefficients
	fprintf (out," %ld",trcid[i]+1);
	//  id of transmission/radiation coefficients
	fprintf (out," %ld",trrid[i]+1);
      }
      if (bc[i]==det_climcond) {
	// prescribed detail climatic conditions
	
	//  id of detail climatic conditions
	fprintf (out," %ld",nvid[i]+1);
      }
      if (bc[i]==gen_climcond){
	// prescribed detail climatic conditions
	
	//  id of detail climatic conditions
	fprintf (out," %ld",nvid[i]+1);
      }
      
    }
  }
  
}




/**
  The function merges BC defined by the parameter ebc to the given object.

  @param ebc - merged element BC

  @retval 0 - on success
  @retval 1 - element has already assigned BC with different values
  @retval 2 - element has already assigned different type of BC

  Created by Tomas Koudelka 10.2010
*/
long loadelt::merge(loadelt &ebc)
{
  long i;
  for (i=0; i<nbo; i++)
  {
    if (ebc.bc[i] == no_bc)  // no BC -> no merging
      continue;
    if (bc[i])
    {
      if (bc[i] != ebc.bc[i])  // different types of BC cannot be merged
        return 2;

      if (nvid[i] != ebc.nvid[i]) // different id of bnodval -> different evolution of nodal values
        return 1;
      if (bc[i] >= presc_trmiss)
      {
        if (trcid[i] != ebc.trcid[i]) // different id of bnodval -> different evolution of transmission coefficient
          return 1;
        if (trrid[i] != ebc.trrid[i]) // different id of bnodval -> different evolution of radiation coefficient
          return 1;
      }
      // if no of the above differencies are found -> the condition are identical -> go to next boundary object
      continue;
    }
    else // no BC prscribed for i-th boundary object -> merge with ebc
    {
      bc[i] = ebc.bc[i];
      nvid[i] = ebc.nvid[i];
      trcid[i] = ebc.trcid[i];
      trrid[i] = ebc.trrid[i];
    }
  }
  return 0;
}



/**
  The function shifts indices of bnodvalt objects of transmission coefficients by
  nnv value and radiation coefficients by ntrc value. The shift is caused by
  a separate storage of bnodvalt objects for nodal values, transmission coefficients
  and radiation coefficients.

  @param nnv  - total number of bnodvalt objects for nodal values
  @param ntrc - total number of bnodvalt objects for transmission coefficients

  @return The function does not return anything, only indices are rewritten.

  Created by Tomas Koudelka, 10.2010
*/
void loadelt::renumber_id(long nnv, long ntrc)
{
  long i;
  for (i=0; i<nbo; i++)
  {
    if (bc[i] == no_bc)  // no BC -> no shift of indices
      continue;
    if (bc[i] >= presc_trmiss)
    {
      // Newton BC -> indices have to be shifted 
      // indices of bnodvalt objects of transmission coefficients 
      // have to be shifted by number of objects of nodal values
      trcid[i] += nnv;
      // indices of bnodvalt objects of radiation coefficients 
      // have to be shifted by number of objects of nodal values and
      // transmission coefficients
      trrid[i] += ntrc+nnv;
    }
  }  
}



/**
   function assembles indicators of boundary conditions
   
   @param bc - array of indicators of boundary conditions
   
   JK, 24.11.2008
*/
void loadelt::give_bc (bocontypet *av)
{
  long i;
  
  //  loop over boundary objects - end nodes, edges or surfaces
  for (i=0;i<nbo;i++){
    av[i]=bc[i];
  }
}

/**
   function assembles nodal values defined on boundaries
   there are values for all boundary objects on elements
   
   @param lcid - load case id
   @param fe - array of actual nodal values
   
   JK, 24.11.2008
*/
void loadelt::give_nodval (long lcid,vector &fe)
{
  long i,j,k,l,ndofn,ipp;
  vector av, nval;
  ivector bnodes(nnbo); // array of node numbers on boundary object
  
  //  all components are set to zero
  fillv (0.0,fe);
  
  //  first element material point
  ipp=Tt->elements[eid].ipp[0][0];
  
  k=0;
  //  loop over boundary objects - end nodes, edges or surfaces
  for (i=0;i<nbo;i++){
    if (bc[i]==no_bc){
      k+=nnbo;
    }

    if (bc[i]==presc_flux || bc[i]>=presc_trmiss){
      //  prescribed boundary fluxes (Neumann boundary condition)
      //  or prescribed external vaules (Newton boundary condition)
      
      //  id of nodal values
      l=nvid[i];
      
      if (l<0)
	print_err("negative id of nodal values",__FILE__,__LINE__,__func__);
      
      reallocv(RSTCKVEC(nnbo, av));
      
      //  nodal values for required boundary object
      Tb->lc[lcid].nodval[l].give_val(Tp->time, av);
      
      for (j=0;j<nnbo;j++){
	fe[k]=av[j];
	k++;
      }
    }
    
    
    
    
    if (bc[i]==det_climcond){
      /* *************************************************************** */
      /*  computation of boundary fluxes from climate condition records  */
      /* *************************************************************** */
      
      //  id of nodal values
      l=nvid[i];
      
      if (l<0)
	print_err("negative id of nodal values",__FILE__,__LINE__,__func__);
      
      if(l > Tb->lc[lcid].ncc-1){
	print_err("id of climatic conditions is greater than the number of clim. conditions", __FILE__, __LINE__, __func__);
	abort ();
      }
      
      //  determination of nodes on required boundary object (end node, edge, surface)
      Tt->give_bonodes (eid,i,bnodes);
      
      //  loop over boundary nodes
      for (j=0;j<nnbo;j++){
	//  number of DOFs of node
	ndofn = Tt->give_ndofn (bnodes[j]);
	//  nodal values of computed quantities
        reallocv(RSTCKVEC(ndofn, nval));
	//  nodal values on boundary
	nodalval(bnodes[j], nval);
	
	if(nvid[i] < 0)
	  fe[k]=0.0;//for growing_np_problem
	else
	  fe[k]=Tb->lc[lcid].climcond[l].give_flux(lcid, nval.a, ipp, bnodes[j]);
	
	k++;
      }
    }//  end of the statement if (bc[i]==det_climcond){

    if (bc[i]==gen_climcond){
      /* *************************************************************** */
      /*  computation of boundary fluxes from climate condition records  */
      /* *************************************************************** */
      
      //  id of nodal values
      l=nvid[i];
      
      if (l<0)
	print_err("negative id of nodal values",__FILE__,__LINE__,__func__);
      
      if(l > Tb->lc[lcid].ncc2-1){
	print_err("id of climatic conditions is greater than the number of clim. conditions", __FILE__, __LINE__, __func__);
	abort ();
      }
      
      //  determination of nodes on required boundary object (end node, edge, surface)
      Tt->give_bonodes (eid,i,bnodes);
      
      //  loop over boundary nodes
      for (j=0;j<nnbo;j++){
	//  number of DOFs of node
	ndofn = Tt->give_ndofn (bnodes[j]);
	//  nodal values of computed quantities
        reallocv(RSTCKVEC(ndofn, nval));
	//  nodal values on boundary
	nodalval(bnodes[j], nval);
	
	if(nvid[i] < 0)
	  fe[k]=0.0;//for growing_np_problem
	else
	  fe[k]=Tb->lc[lcid].climcond2[l].give_flux(lcid, nval.a, ipp, bnodes[j]);
	
	k++;
      }
    }//  end of the statement if (bc[i]==gen_climcond){


  }
}


/**
   function assembles nodal values defined on boundaries
   there are values for all boundary objects on elements
   
   @param lcid - load case id
   @param cid - contribution id
   @param fe - array of actual nodal values
   
   JK, 24.11.2008
*/
void loadelt::give_external_nodval (long lcid,long cid,vector &fe)
{
  long i,j,k,l,ndofn,ipp;
  vector nval;
  vector av;
  ivector bnodes(ASTCKIVEC(nnbo));  //  array of node numbers on boundary object
  
  //  all components are set to zero
  fillv (0.0,fe);
  
  //  first element material point
  ipp=Tt->elements[eid].ipp[0][0];
  
  k=0;
  //  loop over boundary objects - end nodes, edges or surfaces
  for (i=0;i<nbo;i++){
    if (bc[i]==no_bc){
      k+=nnbo;
    }

    if (bc[i]==presc_flux || bc[i]>=presc_trmiss){
      //  prescribed boundary fluxes (Neumann boundary condition)
      //  or prescribed external vaules (Newton boundary condition)
      
      //  id of nodal values
      l=nvid[i];
      
      if (l<0)
	print_err("negative id of nodal values",__FILE__,__LINE__,__func__);
      
      reallocv(RSTCKVEC(nnbo, av));
      
      //  nodal values for required boundary object
      Tb->lc[lcid].nodval[l].give_val(Tp->time, av);
      
      for (j=0;j<nnbo;j++){
	fe[k]=av[j];
	k++;
      }
    }
    
    
    if (bc[i]==gen_climcond){
      /* *************************************************************** */
      /*  computation of boundary fluxes from climate condition records  */
      /* *************************************************************** */
      
      //  id of nodal values
      l=nvid[i];
      
      if (l<0)
	print_err("negative id of nodal values",__FILE__,__LINE__,__func__);
      
      if(l > Tb->lc[lcid].ncc2-1){
	print_err("id of climatic conditions is greater than the number of clim. conditions", __FILE__, __LINE__, __func__);
	abort ();
      }
      
      //  determination of nodes on required boundary object (end node, edge, surface)
      Tt->give_bonodes (eid,i,bnodes);
      
      //  loop over boundary nodes
      for (j=0;j<nnbo;j++){
	//  number of DOFs of node
	ndofn = Tt->give_ndofn (bnodes[j]);
	//  nodal values of computed quantities
        reallocv(RSTCKVEC(ndofn, nval));
	//  nodal values on boundary
	nodalval(bnodes[j], nval);
	
	if(nvid[i] < 0)
	  fe[k]=0.0;//for growing_np_problem
	else
	  fe[k]=Tb->lc[lcid].climcond2[l].external_nodval(cid);
	
	k++;
      }
    }//  end of the statement if (bc[i]==gen_climcond){
  }
}


/**
   function assembles nodal values of transmission coefficients
   there are values for all boundary objects on elements
   
   @param lcid - load case id
   @param trc - array of actual transmission coefficients
   
   JK, 24.11.2008, corrected by TKr 04/05/2011
*/
void loadelt::give_trc (double time,bnodvalt *nodval,long lcid,long cid,vector &trc)
{
  long i, j, k, l, ndofn;
  vector nval, av;
  ivector bnodes;
  
  //  all components are set to zero
  fillv (0.0,trc);
  
  k=0;
  //  loop over boundary objects - end nodes, edges or surfaces
  for (i=0;i<nbo;i++){
    
    if (bc[i]==no_bc){
      k+=nnbo;
    }
    
    if (bc[i]>=presc_trmiss){
      //  prescribed transfer on boundary (Cauchy's boundary condition)
      
      //  id of transmisson coefficient
      l=trcid[i];
      
      if (l<0)
	print_err("negative id of transmission coefficients",__FILE__,__LINE__,__func__);
      
      reallocv (RSTCKVEC(nnbo,av));
      
      //  nodal values for required boundary object
      //Tb->lc[lcid].nodval[l].give_val (Tp->time,av);
      nodval[l].give_val (time,av);
      
      for (j=0;j<nnbo;j++){
	trc[k]=av[j];
	k++;
      }
    }
    
    if (bc[i]==gen_climcond){
      
      //  id of nodal values
      l=nvid[i];
      
      if (l<0)
	print_err("negative id of nodal values",__FILE__,__LINE__,__func__);
      
      if(l > Tb->lc[lcid].ncc2-1){
	print_err("id of climatic conditions is greater than the number of clim. conditions", __FILE__, __LINE__, __func__);
	abort ();
      }
      
      //  array of node numbers on boundary object
      reallocv(RSTCKIVEC(nnbo, bnodes));
      //  determination of nodes on required boundary object (end node, edge, surface)
      Tt->give_bonodes (eid,i,bnodes);
      
      //  loop over boundary nodes
      for (j=0;j<nnbo;j++){
	//  number of DOFs of node
	ndofn = Tt->give_ndofn (bnodes[j]);
	//  nodal values of computed quantities
        reallocv(RSTCKVEC(ndofn, nval));
	//  nodal values on boundary
	nodalval(bnodes[j], nval);
	
	if(nvid[i] < 0)
	  trc[k]=0.0;//for growing_np_problem
	else
	  trc[k]=Tb->lc[lcid].climcond2[l].transmission_coeff (lcid,cid,nval[1]);
	
	k++;
      }
    }//  end of the statement if (bc[i]==gen_climcond){
  }
}



/**
   function assembles nodal values of transmission/radiation coefficients
   there are values for all boundary objects on elements
   
   @param lcid - load case id
   @param trr - array of actual transmission/radiation coefficients
   
   JK, 24.11.2008, corrected by TKr 04/05/2011
*/
void loadelt::give_trr (double time,bnodvalt *nodval,vector &trr)
{
  long i,j,k,l;
  vector av;
  
  //  all components are set to zero
  fillv (0.0,trr);
  
  k=0;
  //  loop over boundary objects - end nodes, edges or surfaces
  for (i=0;i<nbo;i++){
    
    if (bc[i]==no_bc){
      k+=nnbo;
    }
    
    if (bc[i]>=presc_trmiss){
      //  prescribed transfer on boundary (Cauchy's boundary condition)
      
      //  id of transmission coefficient emulating the radiation
      l=trrid[i];
      
      if (l<0)
	print_err("negative id of transmission coefficients",__FILE__,__LINE__,__func__);
      
      reallocv (RSTCKVEC(nnbo,av));
      
      //  nodal values for required boundary object
      //Tb->lc[lcid].nodval[l].give_val (Tp->time,av);
      nodval[l].give_val (time,av);
      
      for (j=0;j<nnbo;j++){
	trr[k]=av[j];
	k++;
      }
    }
  }
}

