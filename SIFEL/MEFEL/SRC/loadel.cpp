#include "loadel.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "globmat.h"
#include "matrix.h"
#include "intools.h"
#include "elemhead.h"
#include "elemswitch.h"



/**
  The constructor inializes attributes to zero values.

  Created by JK,
*/
loadel::loadel()
{
  //  element id
  eid = 0L;
  //  the number of DOFs on element
  ndofe = 0L;
  //  the type of element
  tel = noelload;
  // indicator of beam load meaning
  elm = noelloadmean;
  //  the number of edges
  ned = 0L;
  //  the number of nodes on edge
  nned = 0L;
  //  the number of surfaces
  nsurf = 0L;
  //  the number of nodes on surface
  nnsurf = 0L;
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
  // number of point load on beam element
  npnt = 0L;
  
  //  components of edge load
  nodvale = NULL;
  //  components of surface load
  nodvals = NULL;
  //  components of volume load
  nodvalv = NULL;
  //  offset of beam continuous load from the first node
  la = NULL;
  //  length of applied beam load
  lf = NULL;
  //  edge load indicators
  le = NULL;
  //  surface load indicators
  ls = NULL;
}



/**
  The destructor deallocates used memory.

  Created by JK,
*/
loadel::~loadel()
{
  delete [] nodvale;
  delete [] nodvals;
  delete [] nodvalv;
  delete [] le;
  delete [] ls;
  delete [] la;
  delete [] lf;
}



/**
  The function reads element load from the opened text file
  given by the parameter in.

  @param in - pointer to the opened text file

  @retval 0 - on success

  Created by JK,
  Modified by Tomas Koudelka, 06.2009
*/
long loadel::read (XFILE *in)
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
  //  8 - beam_load
  xfscanf (in, "%m", &elloadtype_kwdset, (int*)&tel);

  switch (tel){
  case edge:
    readedgeload (in);
    edgeload ();
    break;  
  case surface:
    readsurfaceload (in);
    surfaceload ();
    break;  
  case volume:
    readvolumeload (in);
    volumeload ();
    break;
  case edge_surface:
    readedgeload (in);
    edgeload ();
    readsurfaceload (in);
    surfaceload ();
    break;
  case edge_volume:
    readedgeload (in);
    edgeload ();
    readvolumeload (in);
    volumeload ();
    break;
  case surface_volume:
    readsurfaceload (in);
    surfaceload ();
    readvolumeload (in);
    volumeload ();
    break;
  case edge_surface_volume:
    readedgeload (in);
    edgeload ();
    readsurfaceload (in);
    surfaceload ();
    readvolumeload (in);
    volumeload ();
    break;
  case beam_load:
    readbeamload(in);
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
void loadel::readvolumeload (XFILE *in)
{
  long i;
  //  the number of DOFs on element
  ndofe=Mt->give_ndofe (eid);
  nodvalv = new double [ndofe];
  nnvv = ndofe;

  for (i=0;i<ndofe;i++){
    xfscanf (in,"%lf",nodvalv+i);
  }
}



/**
  The function reads element edge load from the opened text file
  given by the parameter in.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by JK,
*/
void loadel::readedgeload (XFILE *in)
{
  long i,j,k,l;
  
  //  the number of edges
  ned=Mt->give_ned (eid);
  //  the number of nodes on edge
  nned=Mt->give_nned (eid);
  //  the number of approximated functions
  napfun=Mt->give_napfun (eid);
  nnve = ned*nned*napfun;
  nodvale = new double [nnve];
  for (i=0;i<nnve;i++){
    nodvale[i]=0.0;
  }
  //  edge load indicator
  le = new long [ned];

  l=0;
  for (i=0;i<ned;i++){
    xfscanf (in,"%ld",le+i);
    if (le[i]>0){
      for (j=0;j<nned;j++){
	for (k=0;k<napfun;k++){
	  xfscanf (in,"%lf",nodvale+l);
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
void loadel::readsurfaceload (XFILE *in)
{
  long i,j,k,l;

  nsurf=Mt->give_nsurf (eid);
  nnsurf=Mt->give_nnsurf (eid);
  napfun=Mt->give_napfun (eid);
  nnvs = nsurf*nnsurf*napfun;
  nodvals = new double [nnvs];
  for (i=0;i<nnvs;i++){
    nodvals[i]=0.0;
  }
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
	  xfscanf (in,"%lf",nodvals+l);
	  l++;
	}
      }
    }
    else{  l+=nnsurf*napfun; }
  }
}



void loadel::readbeamload(XFILE *in)
{
  long i, j, k;
  double l, aux;

  //  the number of edges
  ned=1;
  //  the number of nodes on edge
  nned=Mt->give_nned(eid);
  //  the number of approximated functions, i.e. the number of DOFs at element nodes
  napfun=Mt->give_napfun(eid);
  // element length
  l = Mt->give_length(eid);

  // read load meaning indicator
  xfscanf(in, "%m", &elloadmeaning_kwdset, (int*)&elm);

  // set load meaning to edge load indicator
  le=new long[ned];
  le[0] = elm;

  switch(elm)
  {
    case load_contmech_glob:      
    case load_contmech_loc:
      // definition of linear continuous mechanical load on a part of beam element 
      //        fa ||||||||| fb
      //           vvvvvvvvv
      //   a +------------------+ b
      //     | la  |   lf  |
      //     +-----+-------+
      //
      //   a = the first node of the beam topology
      //   b = the last node of the beam topology
      //   la = distance of the load beginning from node a
      //   lf = continuous load length
      //   fa = value of continuous load on side from the node a
      //   fb = value of continuous load on side from the node b     
      nnve = 2*napfun;
      nodvale = new double[nnve];
      la = new double[napfun];
      lf = new double[napfun];
      memset(nodvale, 0, sizeof(*nodvale)*nnve);
      memset(la, 0, sizeof(*la)*napfun);
      memset(lf, 0, sizeof(*lf)*napfun);
      k=0;
      for (i=0; i<napfun; i++){
        // read distance of continous load beginning from the first beam node and length of continuos load
        xfscanf(in, "%le%le", la+i, lf+i);
        if ((la[i] < 0.0) || (lf[i] < 0.0)){
          print_err("invalid offset of the load beginning or load length on element %ld (la=%le, lf = %le)",
                    __FILE__, __LINE__, __func__, eid+1, la[i], lf[i]);
          abort();
        }
        aux = la[i] + lf[i];
        if (aux > l){
          print_err("length %le of element %ld is less than the total length of applied load with offset (la+lf = %le)",
                    __FILE__, __LINE__, __func__, l, eid+1, aux);
          abort();
        }
        // read initial and final value of continuous load on the given interval
        xfscanf(in, "%le%le", nodvale+k, nodvale+k+1);
        k += 2;
      }
      break;
    case load_pointmech_glob:      
    case load_pointmech_loc:
      // definition of mechanical point load on a beam element 
      //        Fa_i |
      //             v
      //   a +------------------+ b
      //     | la_i  |
      //     +-------+
      //
      //   a = the first node of the beam topology
      //   b = the last node of the beam topology
      //   la_i = distance of the i-th point load from the node a
      //   Fa_i = value of point load
      xfscanf(in, "%ld", &npnt);
      nnve = npnt*napfun;
      nodvale = new double[nnve];
      la = new double[nnve];
      memset(nodvale, 0, sizeof(*nodvale)*nnve);
      memset(la, 0, sizeof(*la)*nnve);
      k=0;
      for (i=0; i<npnt; i++){
        for(j=0; j<napfun; j++){
          xfscanf(in, "%le%le", la+k, nodvale+k);
          if (la[k] < 0.0 || la[k] > l)
          {
            print_err("position la=%le of beam point load on element %ld is out of range <0.0; %le>",
                      __FILE__, __LINE__, __func__, la[k], eid+1, l);
            abort();
          }
          k++;
        }
      }
      break;
    case load_thermgrad_y:
      // definition of load by thermal gradient in y direction
      //   dT/dy = \frac{\Delta T}/{h_y} =
      //         = \frac{\Delta T_d - \Delta T_h}/{h_y} =
      //
      //     ~~~~~~~~~~~~~~~~~~~~~~~      \Delta T_h  -----+
      //  a +-----------------------+ b                h_y |
      //     ~~~~~~~~~~~~~~~~~~~~~~~      \Delta T_d  -----+
      //    |
      //    |
      //    v y
      nnve = 2;
      nodvale = new double[nnve];      
      memset(nodvale, 0, sizeof(*nodvale)*nnve);
      // read temperature change dT over the profile height (i=0) and 
      // profile height h_y (i=1)
      for(i=0; i<nnve; i++){
        xfscanf(in, "%le", nodvale+i);
      }
      break;
    case load_thermgrad_z:
      // definition of load by thermal gradient in z direction
      //   dT/dz = \frac{\Delta T}/{h_z} =
      //         = \frac{\Delta T_d - \Delta T_h}/{h_z} =
      //
      //     ~~~~~~~~~~~~~~~~~~~~~~~      \Delta T_h  -----+
      //  a +-----------------------+ b                h_z |
      //     ~~~~~~~~~~~~~~~~~~~~~~~      \Delta T_d  -----+
      //    |
      //    |
      //    v z
      nnve = 2;
      nodvale = new double[nnve];      
      memset(nodvale, 0, sizeof(*nodvale)*nnve);
      // read temperature change dT over the profile height (i=0) and 
      // profile height h_z (i=1)
      for(i=0; i<nnve; i++){
        xfscanf(in, "%le", nodvale+i);
      }
      break;
    default:
      print_err("unknown beam load meaning %d is required on element %ld", __FILE__, __LINE__, __func__,  int(elm), eid+1);
      abort();
  }
}



/**
  The function computes nodal values of load caused by volume load.
  Volume load defined at nodes is integrated over element and nodal
   values are obtained and stored in the array nodval
   
  @return The function does not return anything.

  Created by JK,
*/
void loadel::volumeload ()
{
  //  number of DOFs on element
  long ndofe = Mt->give_ndofe (eid);
  matrix lm(ndofe, ndofe);

  // allocate vector of nodal forces which is data member of the loadel class
  reallocv (ndofe,nf);
  
  //  load matrix
  loadmat (eid,lm);

  //  load matrix multiplied by nodal values
  mxv (lm.a,nodvalv,nf.a,ndofe,ndofe);
}



/**
  Function computes nodal forces caused by edge load.
  Nodal forces are stored in vectors nf.

  @return The function does not return anything.

  Created by JK,
*/
void loadel::edgeload ()
{
  ndofe=Mt->give_ndofe (eid);
  allocv (ndofe,nf); // this MUST be allocated on heap because the allocation of member array 
                     // loadel::nf must be persistent even after return from this function

  switch (Mt->give_elem_type (eid)){
  case bar2d:{
    break;
  }
  case beam2d:{
    Beam2d->nodeforces (eid,le,nodvale,nf);
    break;
  }
  case beam3d:{
    Beam3d->nodeforces (eid,le,nodvale,nf);
    break;
  }
  case planeelementlt:{
    Pelt->nodeforces (eid,le,nodvale,nf);
    break;
  }
  case planeelementrotlt:{
    Perlt->nodeforces (eid,le,nodvale,nf);
    break;
  }
  case planeelementqt:{
    Peqt->nodeforces (eid,le,nodvale,nf);
    break;
  }
  case planeelementsubqt:{
    Pesqt->nodeforces (eid,le,nodvale,nf);
    break;
  }
  case planeelementlq:{
    Pelq->nodeforces (eid,le,nodvale,nf);
    break;
  }
  case planeelementrotlq:{
    Perlq->nodeforces (eid,le,nodvale,nf);
    break;
  }
  case planeelementqq:{
    Peqq->res_nodeforces (eid,le,nodvale,nf);
    break;
  }
    
  case dkqel:{
    Dkqelem->node_forces_edge (eid,le,nodvale,nf);
    break;
  }


  case axisymmlt:{
    Asymlt->edgeload (eid,le,nodvale,nf);
    break;
  }
  case axisymmlq:{
    Asymlq->edgeload (eid,le,nodvale,nf);
    break;
  }
  case axisymmqq:{
    Asymqq->edgeload (eid,le,nodvale,nf);
    break;
  }
  case axisymmcq:{
    Asymcq->edgeload (eid,le,nodvale,nf);
    break;
  }
  default:{
    print_err("unknown element type is required (eid=%ld)", __FILE__, __LINE__, __func__, eid+1);
  }
  }
  
}



/**
  Function computes nodal forces caused by surface load
  nodal forces are stored in vectors nf

  @return The function does not return anything.

  Created by JK,
*/
void loadel::surfaceload ()
{
  ndofe=Mt->give_ndofe (eid);
  allocv (ndofe,nf); // this MUST be allocated on heap
  
  switch (Mt->give_elem_type (eid)){
  case planeelementrotlq:{
    Perlq->node_forces_surf (eid,nodvals,nf);
    break;
  }
  case dktel:{
    Dkt->areaforces(eid, nodvals, nf);
    break;
  }

  case quadkirch:{
    Qkirch->surfload (0,eid,nodvals,nf);
    break;
  }
    
  case dkqel:{
    Dkqelem->node_forces_surf (eid,nodvals,nf);
    break;
  }

  case shellqelem:{
    Shq->node_forces_surf (eid,nodvals,nf);
    break;
  }
    
  case lineartet:{
    Ltet->node_forces_surf (0,eid,ls,nodvals,nf);
    break;
  }
  case quadrtet:{
    Qtet->node_forces_surf (0,eid,ls,nodvals,nf);
    break;
  }
  case linearhex:{
    Lhex->node_forces_surf (0,eid,ls,nodvals,nf);
    break;
  }
  case quadrhex:{
    Qhex->node_forces_surf (0,eid,ls,nodvals,nf);
    break;
  }
  default:{
    print_err("unknown element type is required (eid=%ld)", __FILE__, __LINE__, __func__, eid+1);
  }
  }
}


void loadel::beamload()
{
  // allocation of nodal force vector due to element load
  ndofe=Mt->give_ndofe (eid);
  allocv(ndofe, nf); // this MUST be allocated on heap
  
  switch (Mt->give_elem_type (eid)){
    case beam2d:
      Beam2d->beamnodeforces (eid, elm, nnve, la, lf, nodvale, nf);
      break;
      /*    case beam3d:
      Beam3d->beamnodeforces (eid,le,nodvale,nf);
      break;*/
    default:
      print_err("unknown element type is required (eid=%ld)", __FILE__, __LINE__, __func__, eid+1);
      abort();
  }
}



/**
  The function reads element load data from the opened text file
  given by the parameter in. It is used in the mechprep preprocessor.

  @param in - pointer to the opened text file
  @param lc  - total number of load cases
  @param slc - pointer to array with number of subloadcases for individual load cases, defalut value is NULL

  @retval 0 - on success
  @retval 1 - wrong load case id
  @retval 2 - unknown type of load is required
  @retval 3 - unknown load indicator is required

  Created by Tomas Koudelka
*/
long loadel::read_prep(XFILE *in, long lc, long *slc=NULL)
{
  long i, j, ncomp;

  xfscanf (in, "%k%ld", "lc_id", &nlc);
  if ((nlc < 1) || (nlc > lc))
  {
    print_err("Line %ld, col %ld, file %s: load case id %ld out of range <1,%ld>",
              __FILE__, __LINE__, __func__, in->line, in->col, in->fname, nlc, lc);
    return 1;
  }
  if (slc)
  {
    if (Mp->tprob == forced_dynamics){
      if (xfscanf(in, "%+k%ld", "slc_id", &nslc) == 2){
        if ((nslc < 1) || (nslc > slc[nlc-1])){
          print_err("Line %ld, col %ld, file %s: subload case id %ld out of range <1,%ld>",
                    __FILE__, __LINE__, __func__, in->line, in->col, in->fname, nslc, slc[nlc-1]);
          return 1;
        }
      }
      else
        nslc = 0; // no slc_id keyword => static loadcase is being specified -> hidden subload case id is 0
    }
    else{
      xfscanf (in, "%k%ld", "slc_id", &nslc);
      if ((nslc < 1) || (nslc > slc[nlc-1])){
        print_err("Line %ld, col %ld, file %s: subload case id %ld out of range <1,%ld>",
                  __FILE__, __LINE__, __func__, in->line, in->col, in->fname, nslc, slc[nlc-1]);
        return 1;
      }
    }
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
      nodvale = new double [nnve];
      memset(nodvale, 0, sizeof(*nodvale)*nnve);
      for(i=0; i<ned; i++)
      {      
        xfscanf (in,"%k%ld","coord_sys", le+i);
        if ((le[i] > -1) && (le[i] < 2))
        {
          if (le[i])
          {
            xfscanf(in, "%k", "load_comp");    
            for(j=0; j<ncomp; j++)
              xfscanf(in, "%le", nodvale+i*ncomp+j);
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
      nodvals = new double [nnvs];
      memset(nodvals, 0, sizeof(*nodvals)*nnvs);
      for(i=0; i<nsurf; i++)
      {
        xfscanf (in,"%k%ld","coord_sys", ls+i);
        if ((ls[i] > -1) && (ls[i] < 3))
        {
          if (ls[i])
          {
            xfscanf(in, "%k", "load_comp");
            for(j=0; j<ncomp; j++)
              xfscanf(in, "%le", nodvals+i*ncomp+j);
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
      nodvalv = new double [nnvv];
      memset(nodvalv, 0, sizeof(*nodvalv)*nnvv);
      xfscanf(in, "%k", "load_comp");
      for(i=0; i<nnvv; i++)
        xfscanf(in, "%le", nodvalv+i);
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
long loadel::print(FILE *out, int ident=0)
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
		  fprintf(out, " % le", nodvale[k]);
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
		  fprintf(out, " % le", nodvals[k]);
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
	fprintf(out, " % le", nodvalv[i]);
      fprintf(out, "\n");
    }
  return(0);
}



/**
  The function merges load defined by the parameter le to the given object.

  @param lel - merged element load

  @retval 0 - on success
  @retval 1 - incompatible number of load components
  @retval 2 - incompatible type of load components 

  Created by Tomas Koudelka 06.2009
*/
long loadel::merge(loadel &lel)
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
      nodvale = new double [nnve];
      memset(nodvale, 0, sizeof(*nodvale)*nnve);
      for (i=0; i<ned; i++)
        le[i] = lel.le[i];
      for (i=0; i<nnve; i++)
        nodvale[i] = lel.nodvale[i];
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
        nodvale[i] += lel.nodvale[i];
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
      nodvals = new double [nnvs];
      memset(nodvals, 0, sizeof(*nodvals)*nnvs);
      for (i=0; i<nsurf; i++)
        ls[i] = lel.ls[i];
      for (i=0; i<nnvs; i++)
        nodvals[i] = lel.nodvals[i];
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
        nodvals[i] += lel.nodvals[i];
    }
  }

  if (lel.nodvalv)  
  {
    if (nodvalv == NULL)
    {
      nnvv = lel.nnvv;
      nodvalv = new double [nnvv];
      memset(nodvalv, 0, sizeof(*nodvalv)*nnvv);
      for (i=0; i<nnvv; i++)
        nodvalv[i] = lel.nodvalv[i];
    }
    else
    {
      // check for incompatible number of load components
      if (nnvv != lel.nnvv)
        return 1;
      for (i=0; i<nnvv; i++)
        nodvalv[i] += lel.nodvalv[i];
    }
  }
  set_load_type();
  return 0;
}



/**
  The function detects actual type of load depending on 
  the allocated nodal load values.

  Created by Tomas Koudelka, 07.2009
*/
void loadel::set_load_type()
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
