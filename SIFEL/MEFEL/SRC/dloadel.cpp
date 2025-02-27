#include "dloadel.h"
#include "vector.h"
#include "intools.h"
#include "gfunct.h"
#include "iotools.h"
#include "global.h"
#include "mechtop.h"
#include "globmat.h"
#include "matrix.h"
#include "elemhead.h"
#include "elemswitch.h"



/**
  The constructor inializes attributes to zero values.

  Created by JK,
*/
dloadel::dloadel()
{
  //  element id
  eid = 0L;
  //  the number of DOFs on element
  ndofe = 0L;
  //  the type of element
  tel = elloadtype(0);
  //  the number of edges
  ned = 0;
  //  the number of nodes on edge
  nned = 0;
  //  the number of surfaces
  nsurf = 0;
  //  the number of nodes on surface
  nnsurf = 0;
  //  the number of approximated functions
  napfun = 0L;
  //  the number of components of array nodvale
  nnve = 0;
  //  the number of components of array nodvals
  nnvs = 0;
  //  the number of components of array nodvalv
  nnvv = 0L;
  //  load case id
  nlc = 0L;
  //  subload case id
  nslc = 0L;
  
  //  components of edge load
  nodvale = NULL;
  //  components of surface load
  nodvals = NULL;
  //  components of volume load
  nodvalv = NULL;
  //  edge load indicators
  le = NULL;
  //  surface load indicators
  ls = NULL;
}



/**
  The destructor deallocates used memory.

  Created by JK,
*/
dloadel::~dloadel()
{
  delete [] nodvale;
  delete [] nodvals;
  delete [] nodvalv;
  delete [] le;
  delete [] ls;
}



/**
  The function reads element load from the opened text file
  given by the parameter in.

  @param in - pointer to the opened text file

  @retval 0 - on success

  Created by JK,
  Modified by Tomas Koudelka, 06.2009
*/
long dloadel::read (XFILE *in)
{
  //  element id
  xfscanf(in,"%ld",&eid);
  if (eid<1)
    print_err("nonpositive number of loaded element", __FILE__, __LINE__, __func__);
  if (eid>Mt->ne)
    print_err("number of loaded element (%ld) > total number of elements (%ld)", __FILE__, __LINE__, __func__, eid, Mt->ne);
  eid--;
  
  //  type of element load
  //  1 - volume
  //  2 - edge
  //  3 - surface
  //  4 - edge_surface
  //  5 - edge_volume
  //  6 - surface_volume
  //  7 - edge_surface_volume
  xfscanf (in, "%m", &elloadtype_kwdset, (int*)&tel);

  switch (tel){
  case edge:
    readedgeload (in);
    break;  
  case surface:
    readsurfaceload (in);
    break;  
  case volume:
    readvolumeload (in);
    break;
  case edge_surface:
    readedgeload (in);
    readsurfaceload (in);
    break;
  case edge_volume:
    readedgeload (in);
    readvolumeload (in);
    break;
  case surface_volume:
    readsurfaceload (in);
    readvolumeload (in);
    break;
  case edge_surface_volume:
    readedgeload (in);
    readsurfaceload (in);
    readvolumeload (in);
    break;
  default:
    print_err("unknown type of element load is required", __FILE__, __LINE__, __func__);
  }
  return(0);
}



/**
  The function reads element volume load from the opened text file
  given by the parameter in.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by JK,
*/
void dloadel::readvolumeload (XFILE *in)
{
  long i;
  //  the number of DOFs on element
  ndofe=Mt->give_ndofe (eid);
  nodvalv = new gfunct [ndofe];
  nnvv = ndofe;

  for (i=0;i<ndofe;i++){
    nodvalv[i].read(in);
  }
}



/**
  The function reads element edge load from the opened text file
  given by the parameter in.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by JK,
*/
void dloadel::readedgeload (XFILE *in)
{
  long i,j,k,l;
  
  //  the number of edges
  ned=Mt->give_ned (eid);
  //  the number of nodes on edge
  nned=Mt->give_nned (eid);
  //  the number of approximated functions
  napfun=Mt->give_napfun (eid);
  nnve = ned*nned*napfun;
  // array of load components
  nodvale = new gfunct [nnve];
  //  edge load indicator
  le = new long [ned];

  l=0;
  for (i=0;i<ned;i++){
    xfscanf (in,"%ld",le+i);
    if (le[i]>0){
      for (j=0;j<nned;j++){
	for (k=0;k<napfun;k++){
	  nodvale[l].read(in);
	  l++;
	}
      }
    }
    else{  l+=nned*napfun; }
  }
}



/**
  The function reads element surface load from the opened text file
  given by the parameter in.
   
  @param in - pointer to the input file

  @return The function does not return anything.

  Created by JK,
*/
void dloadel::readsurfaceload (XFILE *in)
{
  long i,j,k,l;

  nsurf=Mt->give_nsurf (eid);
  nnsurf=Mt->give_nnsurf (eid);
  napfun=Mt->give_napfun (eid);
  nnvs = nsurf*nnsurf*napfun;
  // array of load components
  nodvals = new gfunct [nnvs];
  //  surface load indicator
  ls = new long [nsurf];
  
  l=0;
  for (i=0;i<nsurf;i++){
    xfscanf (in,"%ld",ls+i);
    //  ls=0 - surface is not loaded
    //  ls=1 - surface is loaded by load defined in global system
    //  ls=2 - surface is loaded by load defined in local system
    if (ls[i]==1 || ls[i]==2){
      for (j=0;j<nnsurf;j++){
	for (k=0;k<napfun;k++){
	  nodvals[l].read(in);
	  l++;
	}
      }
    }
    else{  l+=nnsurf*napfun; }
  }
}



/**
  The function computes nodal values of load caused by time dependent 
  volume load. Volume load defined at nodes is integrated over element and nodal
  values are obtained and stored in the array nodval.

  @param t - actual time for which the load components should be evaluated
   
  @return The function does not return anything.

  Created by Tomas Koudelka, 25.4.2016
*/
void dloadel::volumeload (double t)
{
  //  number of DOFs on element
  long i, ndofe = Mt->give_ndofe (eid);
  matrix lm(ASTCKMAT(ndofe, ndofe));
  vector nv(ASTCKVEC(nnvv));

  // compute nodal values of load components for the given actual time
  for (i=0; i<nnvv; i++)
    nv(i) = nodvalv[i].getval(t);

  // allocate vector of nodal forces which is data member of the dloadel class
  reallocv(ndofe,nf);
  
  //  load matrix
  loadmat (eid,lm);

  //  load matrix multiplied by nodal values
  mxv (lm,nv,nf);
}



/**
  Function computes nodal forces caused by edge load.
  Nodal forces are stored in vectors nf.

  @param t - actual time for which the load components should be evaluated

  @return The function does not return anything.

  Created by JK,
*/
void dloadel::edgeload (double t)
{
  long i;
  vector nv(ASTCKVEC(nnve));
  ndofe=Mt->give_ndofe (eid);
  reallocv(ndofe,nf);

  // compute nodal values of load components for the given actual time
  for (i=0; i<nnve; i++)
    nv(i) = nodvale[i].getval(t);

  switch (Mt->give_elem_type (eid)){
  case bar2d:{
    break;
  }
  case beam2d:{
    Beam2d->nodeforces (eid,le,nv.a,nf);
    break;
  }
  case beam3d:{
    Beam3d->nodeforces (eid,le,nv.a,nf);
    break;
  }
  case planeelementlt:{
    Pelt->nodeforces (eid,le,nv.a,nf);
    break;
  }
  case planeelementrotlt:{
    Perlt->nodeforces (eid,le,nv.a,nf);
    break;
  }
  case planeelementqt:{
    Peqt->nodeforces (eid,le,nv.a,nf);
    break;
  }
  case planeelementsubqt:{
    Pesqt->nodeforces (eid,le,nv.a,nf);
    break;
  }
  case planeelementlq:{
    Pelq->nodeforces (eid,le,nv.a,nf);
    break;
  }
  case planeelementrotlq:{
    Perlq->nodeforces (eid,le,nv.a,nf);
    break;
  }
  case planeelementqq:{
    Peqq->res_nodeforces (eid,le,nv.a,nf);
    break;
  }

  case axisymmlt:{
    Asymlt->edgeload (eid,le,nv.a,nf);
    break;
  }
  case axisymmlq:{
    Asymlq->edgeload (eid,le,nv.a,nf);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  Function computes nodal forces caused by surface load
  nodal forces are stored in vectors nf

  @param t - actual time for which the load components should be evaluated

  @return The function does not return anything.

  Created by Tomas Koudelka, 25.4.2016
*/
void dloadel::surfaceload (double t)
{
  long i;
  vector nv(ASTCKVEC(nnvs));

  ndofe=Mt->give_ndofe (eid);
  allocv (ndofe,nf); // this MUST be allocated at heap

  // compute nodal values of load components for the given actual time
  for (i=0; i<nnvs; i++)
    nv(i) = nodvals[i].getval(t);
  
  switch (Mt->give_elem_type (eid)){
  case dktel:{
    Dkt->areaforces(eid, nv.a, nf);
    break;
  }
  case lineartet:{
    Ltet->node_forces_surf (0,eid,ls,nv.a,nf);
    break;
  }
  case quadrtet:{
    Qtet->node_forces_surf (0,eid,ls,nv.a,nf);
    break;
  }
  case linearhex:{
    Lhex->node_forces_surf (0,eid,ls,nv.a,nf);
    break;
  }
  case quadrhex:{
    Qhex->node_forces_surf (0,eid,ls,nv.a,nf);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  The function reads element load data from the opened text file
  given by the parameter in. It is used in the mechprep preprocessor.

  @param in - pointer to the opened text file
  @param lc  - total number of load cases

  @retval 0 - on success
  @retval 1 - wrong load case id
  @retval 2 - unknown type of load is required
  @retval 3 - unknown load indicator is required

  Created by Tomas Koudelka
*/
long dloadel::read_prep(XFILE *in, long lc)
{
  long i, j, ncomp;

  xfscanf (in, "%k%ld", "lc_id", &nlc);
  if ((nlc < 1) || (nlc > lc))
  {
    print_err("Line %ld, col %ld, file %s: load case id %ld out of range <1,%ld>",
              __FILE__, __LINE__, __func__, in->line, in->col, in->fname, nlc, lc);
    return 1;
  }
  //  type of element load
  xfscanf (in, "%k%m", "load_type", &elloadtype_kwdset, &tel);
  switch (tel)
  {
    case edge:
      xfscanf (in, "%k%ld", "nedge", &ned);
      xfscanf (in, "%k%ld", "ncomp", &ncomp);
      nnve = ncomp*ned;
      le = new long[ned];
      memset(le, 0, sizeof(*le)*ned);
      nodvale = new gfunct [nnve];
      for(i=0; i<ned; i++)
      {      
        xfscanf (in,"%k%ld","coord_sys", le+i);
        if ((le[i] > -1) && (le[i] < 2))
        {
          if (le[i])
          {
            xfscanf(in, "%k", "load_comp");    
            for(j=0; j<ncomp; j++)
              nodvale[i*ncomp+j].read(in);
          }
        }
        else
        {
          print_err("Line %ld, col %ld, file %s: unknown type %ld of coordinate system is required on %ld. edge of element.", 
                    __FILE__, __LINE__, __func__, in->line, in->col, in->fname, ls[i], i+1);
          return 3;
        }
      }
      break;
    case surface:
      xfscanf (in, "%k%ld", "nsurf", &nsurf);
      xfscanf (in, "%k%ld", "ncomp", &ncomp);
      nnvs = ncomp*nsurf;
      ls = new long[nsurf];
      memset(ls, 0, sizeof(*ls)*nsurf);
      nodvals = new gfunct [nnvs];
      for(i=0; i<nsurf; i++)
      {
        xfscanf (in,"%k%ld","coord_sys", ls+i);
        if ((ls[i] > -1) && (ls[i] < 3))
        {
          if (ls[i])
          {
            xfscanf(in, "%k", "load_comp");
            for(j=0; j<ncomp; j++)
              nodvals[i*ncomp+j].read(in);
          }
        }
        else
        {
          print_err("Line %ld, col %ld, file %s: unknown type %ld of coordinate system is required on %ld. surface of element.", 
                    __FILE__, __LINE__, __func__, in->line, in->col, in->fname, ls[i], i+1);
          return 3;
        }
      }
      break;
    case volume:
      xfscanf (in, "%k%ld", "ncomp", &ncomp);
      nnvv = ncomp;
      nodvalv = new gfunct [nnvv];
      xfscanf(in, "%k", "load_comp");
      for(i=0; i<nnvv; i++)
        nodvalv[i].read(in);
      break;
    default:
      print_err("Line %ld, col %ld, file %s: unknown type %d of element load is required", 
                __FILE__, __LINE__, __func__, in->line, in->col, in->fname, int(tel));
      return 2;
  }
  return(0);
}



/**
  The function prints element load data to the opened text file
  given by the parameter out. It is used in the mechprep preprocessor.

  @param out   - pointer to the opened text file
  @param ident - number of spaces printed on the beginning of each new line except of last

  @retval 0 - on success
  @retval 1 - unknown type of load is required

  Created by Tomas Koudelka, 06.2009
*/
long dloadel::print(FILE *out, int ident=0)
{
  long i, j, k, nnv1;

  //  type of element load
  fprintf (out, " %d", tel);
  if (nodvale)
    {
      nnv1 = nnve/ned;
      k = 0;
      fprintf(out, " ");
      for(i=0; i<ned; i++)
	{
	  fprintf(out, "%ld", le[i]);
	  if (le[i])
	    {
	      for (j=0; j<nnv1; j++)
		{
		  nodvale[k].print(out);
		  k++;
		}
	    }
	  else
	    k += nnv1;
	  if (i == ned-1)
	    fprintf(out, "\n");
	  else
	    fprintf(out, "\n%*c", ident+3,' ');
	}
    }
  if (nodvals)
    {
      nnv1 = nnvs/nsurf;
      k = 0;
      fprintf(out, " ");
      for(i=0; i<nsurf; i++)
	{
	  if (nodvale && (i == 0))
	    fprintf(out, "%*c", ident+2, ' ');
	  
	  fprintf(out, "%ld", ls[i]);
	  if (ls[i])
	    {
	      for (j=0; j<nnv1; j++)
		{
		  nodvals[k].print(out);
		  k++;
		}
	    }
	  else
	    k += nnv1;
	  if (i == nsurf-1)
	    fprintf(out, "\n");
	  else
	    fprintf(out, "\n%*c", ident+3, ' ');
	}
    }
  if (nodvalv)
    {
      if (nodvale || nodvals)
	fprintf(out, "%*c", ident+2, ' ');
      for (i=0; i<nnvv; i++)
	nodvalv[i].print(out);
      fprintf(out, "\n");
    }
  return(0);
}



/**
  The function computes nodal values of element load for the given time.
  The resulting values are stored in the dloadel vector nf.

  @param t - actual time for which the nodal values of element load should be calculated

  @return The function does not return anything, results are stored in the dloadel class member nf.

  Created by Tomas Koudelka, 27.4.2016
*/
void dloadel::compute_load(double t)
{
  switch (tel)
  {
    case noelload:
      break;
    case edge:
      edgeload(t);
      break;
    case surface:
      surfaceload(t);
      break;
    case volume:
      volumeload(t);
      break;
    case edge_surface:
      edgeload(t);
      surfaceload(t);
      break;
    case edge_volume:
      edgeload(t);
      volumeload(t);
      break;
    case surface_volume:
      surfaceload(t);
      volumeload(t);
      break;
    case edge_surface_volume:
      edgeload(t);
      surfaceload(t);
      volumeload(t);
      break;
    case beam_load:
      print_err("beam_load type has not been implemented yet", __FILE__, __LINE__, __func__);
      abort();
  }
}

/**
  The function merges load defined by the parameter lel to the given object.

  @param lel - merged element load

  @retval 0 - on success
  @retval 1 - incompatible number of load components
  @retval 2 - incompatible type of load components 

  Created by Tomas Koudelka 25.4.2016
*/
long dloadel::merge(dloadel &lel)
{
  long i;
  if (lel.nodvale)
  {
    if (nodvale == NULL)
    {
      ned  = lel.ned;
      nnve = lel.nnve;
      le = new long[ned];
      memset(le, 0, sizeof(*le)*ned);
      nodvale = new gfunct [nnve];
      for (i=0; i<ned; i++)
        le[i] = lel.le[i];
      for (i=0; i<nnve; i++)
        nodvale[i].copy(lel.nodvale[i]);
    }
    else
    {
      // check for incompatible number of load components
      if ((ned != lel.ned) || (nnve != lel.nnve))
        return 1;
      // check for incompatible type of load components
      for (i=0; i<ned; i++)
      {
        if (lel.le[i])
        {
          if ((le[i]) && (le[i] != lel.le[i]))
            return 2;
          else
            le[i] = lel.le[i];
        }
      }
      for (i=0; i<nnve; i++)
        nodvale[i].merge(lel.nodvale[i]);
    }
  }

  if (lel.nodvals)  
  {
    if (nodvals == NULL)
    {
      nsurf  = lel.nsurf;
      nnvs = lel.nnvs;
      ls = new long[nsurf];
      memset(ls, 0, sizeof(*ls)*nsurf);
      nodvals = new gfunct [nnvs];
      for (i=0; i<nsurf; i++)
        ls[i] = lel.ls[i];
      for (i=0; i<nnvs; i++)
        nodvals[i].copy(lel.nodvals[i]);
    }
    else
    {
      // check for incompatible number of load components
      if ((nsurf != lel.nsurf) || (nnvs != lel.nnvs))
        return 1;
      // check for incompatible type of load components
      for (i=0; i<nsurf; i++)
      {
        if (lel.ls[i])
        {
          if ((ls[i]) && (ls[i] != lel.ls[i]))
            return 2;
          else
            ls[i] = lel.ls[i];
        }
      }
      for (i=0; i<nnvs; i++)
        nodvals[i].merge(lel.nodvals[i]);
    }
  }

  if (lel.nodvalv)  
  {
    if (nodvalv == NULL)
    {
      nnvv = lel.nnvv;
      nodvalv = new gfunct [nnvv];
      for (i=0; i<nnvv; i++)
        nodvalv[i].copy(lel.nodvalv[i]);
    }
    else
    {
      // check for incompatible number of load components
      if (nnvv != lel.nnvv)
        return 1;
      for (i=0; i<nnvv; i++)
        nodvalv[i].merge(lel.nodvalv[i]);
    }
  }
  set_load_type();
  return 0;
}



/**
  The function detects the actual load type depending on 
  the allocated nodal load values.

  Created by Tomas Koudelka, 25.4.2016
*/
void dloadel::set_load_type()
{
  tel = elloadtype(0);
  if (nodvale)
  {    
    if (nodvals)
    {
      if (nodvalv)
        tel = edge_surface_volume;
      else 
        tel = edge_surface;
    }
    else
    {
      if (nodvalv)
        tel = edge_volume;
      else
        tel = edge;
    }
  }
  else
  {
    if (nodvals)
    {
      if (nodvalv)
        tel = surface_volume;
      else
        tel = surface;
    }
    else
    {
      if (nodvalv)
        tel = volume;
    }
  }
}
