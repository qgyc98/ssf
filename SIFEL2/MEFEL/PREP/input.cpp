#include "input.h"
#include "iotools.h"
#include "kwdset.h"
#include "vector.h"
#include "siftop.h"
#include "gfunct.h"
#include "alias.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "inicd.h"
#include "globprep.h"
#include "loadel.h"
#include "prepalias.h"
#include "hangnode.h"
#include "hngen.h"
#include "bocon.h"
#include "dbcrs.h"
#include "dbmat.h"
#include "entityload.h"
#include "entitytdload.h"
#include "pointset.h"
#include "tempload.h"
#include "aggregator.h"
#include "ipmap.h"
#include "mathem.h"
#include "difcalc.h"
#include <float.h>


/**
 Function   reads problem data from file in a writes it to file out in SIFEL format.

 @param in   - pointer to opened input file
 @param d - pointer to structure with description of preprocessor setup 

 @retval 0 - on succes
 @retval 1 - unable to open temporary file
 @retval 2 - if fails reading number of loading cases
 @retval 3 - if fails reading topology
 @retval 4 - if fails reading materials
 @retval 5 - if fails reading cross-sections
 @retval 6 - if fails reading property file
 @retval 7 - if fails writing output file
*/
long input(XFILE *in, descrip *d)
{
  long err;
  XFILE *itop, *ihn;

  // reading of load cases
  err = input_lc(in);
  if (err)
    return(2);
  // reading of topology
  itop = xfopen(d->topf, "rt");
  if (itop == NULL)
    return(3);
  itop->warning = 1;
  itop->kwdmode = ignore_kwd;
  //in->kwdmode = sequent_mode;
  itop->ignorecase = 1;
  
  
  if (Mp->ssle->prec.pt==boss){
    //  iterative method is preconditioned by the BOSS algorithm
    if (Mp->ssle->prec.agg->impl==2){
      //  the BOSS method is based on the METIS decomposition
      d->paral=2;
    }
  }
  
  if (Mp->ssle->tlinsol==sfeti){
    //  the problem is solved by the FETI method on a single processor
    d->paral=2;
  }


  err = input_siftop(itop, d);
  xfclose(itop);
  if (err)
    return(3);

  // reading of hanging nodes
  Numhn = 0;
  Nod_hang = NULL;
  if (d->hangnf[0])
  {
    ihn = xfopen(d->hangnf, "rt");
    if (ihn == NULL)
      return(3);
    ihn->warning = 1;
    ihn->kwdmode = ignore_kwd;
    //ihn->kwdmode = sequent_mode;
    ihn->ignorecase = 1;
    err = input_hang_nodes(ihn);
    xfclose(ihn);
    if (err)
      return(3);
  }

  // reading material database
  if (d->matsec == no)
    err = input_materials(d->matf, d);
  else
    err = input_materials(in, d);
  if (err)
    return(4);

  // reading cross-section elements
  if (d->crssec == no)
    err = input_crs(d->crf, d);
  else
    err = input_crs(in, d);
  if (err)
    return(5);

  // reading nodal properties
  err = input_nodprop(in);
  switch (err)
  {
    case 0 :
      break;
    case 1 :
      print_err("reading of number of dofs at nodes failed", __FILE__, __LINE__, __func__);
      return(6);
    case 2 :
      print_err("reading of nodal boundary conditions failed", __FILE__, __LINE__, __func__);
      return(6);
    case 3 :
      print_err("reading of coupled dofs at nodes failed", __FILE__, __LINE__, __func__);
      return(6);
    case 4 :
      print_err("reading of springs at nodes failed", __FILE__, __LINE__, __func__);
      return(6);
    case 5 :
      print_err("reading of dof time functions at nodes failed", __FILE__, __LINE__, __func__);
      return(6);
    case 6 :
      print_err("reading of cross-sections at nodes failed", __FILE__, __LINE__, __func__);
      return(6);
    case 7 :
      print_err("reading of local coordinate systems at nodes failed", __FILE__, __LINE__, __func__);
      return(6);
    case 8 :
      print_err("reading of nodal load failed", __FILE__, __LINE__, __func__);
      return(6);
    case 9 :
      print_err("reading of nodal time dependent load failed", __FILE__, __LINE__, __func__);
      return(6);
    case 10 :
      print_err("reading of nodal initial conditions failed", __FILE__, __LINE__, __func__);
      return(6);
    case 11 :
      print_err("reading of nodal initial displacements due to rotation failed", __FILE__, __LINE__, __func__);
      return(6);
    case 12 :
      print_err("reading of nodal temperature load failed", __FILE__, __LINE__, __func__);
      return(6);
    case 13 :
      print_err("reading of periodic boundary condition at nodes failed", __FILE__, __LINE__, __func__);
      return(6);
    case 14 :
      print_err("reading of generated hanging nodes failed", __FILE__, __LINE__, __func__);
      return(6);
    default :
      print_err("unknown error code in nodal sections", __FILE__, __LINE__, __func__);
      return(6);
  }
  // reading nodal properties
  err = input_elemprop(in);
  switch (err)
  {
    case 0 :
      break;
    case 1 :
      print_err("reading of element types failed", __FILE__, __LINE__, __func__);
      return(7);
    case 2 :
      print_err("reading of element material types failed", __FILE__, __LINE__, __func__);
      return(7);
    case 3 :
      print_err("reading of element cross-section failed", __FILE__, __LINE__, __func__);
      return(7);
    case 4 :
      print_err("reading of local coordinate systems at elements failed", __FILE__, __LINE__, __func__);
      return(7);
    case 5 :
      print_err("reading of element load record failed", __FILE__, __LINE__, __func__);
      return(7);
    case 6 :
      print_err("reading of edge load of elements failed", __FILE__, __LINE__, __func__);
      return(7);
    case 7 :
      print_err("reading of surface load of elements failed", __FILE__, __LINE__, __func__);
      return(7);
    case 8 :
      print_err("reading of volume load of elements failed", __FILE__, __LINE__, __func__);
      return(7);
    case 9 :
      print_err("reading of eigenstrain/eigenstresses for elements failed", __FILE__, __LINE__, __func__);
      return(7);
    case 10 :
      print_err("reading of time functions for elements failed", __FILE__, __LINE__, __func__);
      return(7);
    case 11 :
      print_err("reading of time dependent element load record failed", __FILE__, __LINE__, __func__);
      return(7);
    case 12 :
      print_err("reading of time dependent edge load of elements failed", __FILE__, __LINE__, __func__);
      return(7);
    case 13 :
      print_err("reading of time dependent surface load of elements failed", __FILE__, __LINE__, __func__);
      return(7);
    case 14 :
      print_err("reading of time dependent volume load of elements failed", __FILE__, __LINE__, __func__);
      return(7);
    default :
      print_err("unknown error code in element sections", __FILE__, __LINE__, __func__);
      return(7);
  }
  return(0);
}



/**
  The function reads section with files used for the generation of MEFEL input 
  file e.g., topology, material or cros ssection files. It also reads setup of 
  preprocessor.

  @param in - pointer to the opened XFILE structure
  @param d  - structure with input data format description

  @retval 0 - on success
  
  Created by Tomas Koudelka, 7.7.2014

*/
long input_files(XFILE *in, descrip &d)
{
  // reading of line with topology file name
  xfscanf(in, " %1024a", d.topf);
  // reading of line with material database file name
  if (d.matsec == no)
    xfscanf(in, " %1024a", d.matf);
  // reading of line with cross-section database file name
  if (d.crssec == no)
    xfscanf(in, " %1024a", d.crf);
  // reading of line with topology file format indicator,
  xfscanf(in, "%k%m", "mesh_format", &meshform_kwdset, &d.meshfmt);
  // reading of line with edge number indicator
  xfscanf(in, "%k%ld", "edge_numbering", &d.redgn);
  // optional reading of hanging nodes
  d.hangnf[0] = 0;
  xfscanf(in, "%+k %1024a", "hanging_nodes_file", &d.hangnf);
  // optional reading of materials via strings or mechmat
  d.matstr = yes;
  xfscanf(in, "%+k %m", "read_mat_strings", &answertype_kwdset, &d.matstr);
  if (d.matstr == no)
    xfscanf(in, "%k %m", "read_mat_kwd", &answertype_kwdset, &d.matkwd);
  else 
    d.matkwd = no;

  // optional reading of cross section parameters via strings or mechcrsec procedures
  d.crsstr = yes;
  xfscanf(in, "%+k %m", "read_crs_strings", &answertype_kwdset, &d.crsstr);
  if (d.crsstr == no)
    xfscanf(in, "%k %m", "read_crs_kwd", &answertype_kwdset, &d.crskwd);
  else 
    d.crskwd = no;

  return 0;
}



/**
  Input of data about load cases and layers

  @param in - pointer to the opened XFILE structure
  
  Returns:
  @retval 1 - in case of wrong number of load cases
  @retval 2 - in case of wrong number of layers
  @retval 3 - in case of wrong load|subload case index
  @retval 4 - in case of wrong temperature load type 
  @retval 5 - in case of wrong number of prescribed displacements (growing mechanical problem only)
  @retval 6 - in case of multiple assignment of different values of prescribed displacements 
              to the same subload case (growing mechanical problem only)

  created 04.2008 by Tomas Koudelka koudelka@cml.fsv.cvut.cz
*/
long input_lc(XFILE *in)
{
  long i, j, n, check, tnlc;
  long lc_id, slc_id, id;
  char errmsg[1001];
  double tmp;
  

  fprintf(stdout, "\n\nReading of loadcase section . . .");
  memset(errmsg, 0, sizeof(*errmsg)*1001);
  xf_setsec(in, bsec_str[begsec_loadcase]);
  xfscanf(in, "%k%ld", "num_loadcases",  &Nlc);
  tnlc = Nlc / 2;
  Nlay = 0;
  Nslc = NULL;
  Nslc_cum = NULL;
  Tnslc = 0;
  Tlt = NULL;
  Npd = NULL;
  Spd = NULL;
  Nmstrc = 0;
  Mstrc = NULL;  
  switch (Mp->tprob)
  {
    case mat_nonlinear_statics :
    {
      if (tnlc * 2 != Nlc)
      {
        sprintf(errmsg, "Wrong number of load cases\n"  
                        " Number of load cases is %ld\n"  
                        " Number of load cases has to be even\n"  
                        " because you request to solve mat_nonlinear statics.", Nlc);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return (1);
      }
      // reading of temperature load type for particular load cases
      Tlt = new long [Nlc];
      memset(Tlt, 0, sizeof(*Tlt)*Nlc);
      for(i=0; i<Nlc; i++)
      {
        xfscanf(in, "%k%ld", "lc_id", &lc_id);
        if (lc_id < 1)
        {
          sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"  
                          " load case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, lc_id, tnlc);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
        lc_id--;
        xfscanf(in, "%k%ld", "temp_load_type", Tlt+lc_id);
        if ((Tlt[lc_id] < 0) || (Tlt[lc_id] > 3))
        {
          sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"  
                          " temperature load type %ld is out of range <0,3>", in->line, in->col, in->fname, Tlt[lc_id]);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 4;
        }
        // read macro stress components for homogenization problems
        if (Mp->homog == 3)
        {
          if (Mstrc == NULL)
          {
            Mstrc = new double*[Nlc];
            memset(Mstrc, 0, sizeof(*Mstrc)*Nlc);
          }
          xfscanf(in, "%k%ld", "num_macro_stress_comp", &Nmstrc);
          if (Nmstrc > 0)
          {
            Mstrc[i] = new double[Nmstrc];
            memset(Mstrc[i], 0, sizeof(*Mstrc[i])*Nmstrc);
            xfscanf(in, "%k", "macro_stress_comp");
            for (j=0; j<Nmstrc; j++)
              xfscanf (in, "%lf", Mstrc[i]+j);
          }
        }
        // read macro strain components for homogenization problems
        if (Mp->homog == 4)
        {
          if (Mstrc == NULL)
          {
            Mstrc = new double*[Nlc];
            memset(Mstrc, 0, sizeof(*Mstrc)*Nlc);
          }
          xfscanf(in, "%k%ld", "num_macro_strain_comp", &Nmstrc);
          if (Nmstrc > 0)
          {
            Mstrc[i] = new double[Nmstrc];
            memset(Mstrc[i], 0, sizeof(*Mstrc[i])*Nmstrc);
            xfscanf(in, "%k", "macro_strain_comp");
            for (j=0; j<Nmstrc; j++)
              xfscanf (in, "%lf", Mstrc[i]+j);
          }
        }
        // read macro stress/strain components for homogenization problems
        if (Mp->homog == 9)
        {
          if (Mstrc == NULL)
          {
            Mstrc = new double*[Nlc];
            memset(Mstrc, 0, sizeof(*Mstrc)*Nlc);
            Mstrct = new strastre*[Nlc];
            memset(Mstrct, 0, sizeof(*Mstrct)*Nlc);
          }
          xfscanf(in, "%k%ld", "num_macro_comp", &Nmstrc);
          if (Nmstrc > 0)
          {
            Mstrc[i] = new double[Nmstrc];
            memset(Mstrc[i], 0, sizeof(*Mstrc[i])*Nmstrc);
            Mstrct[i] = new strastre[Nmstrc];
            memset(Mstrct[i], 0, sizeof(*Mstrct[i])*Nmstrc);
            xfscanf(in, "%k", "macro_comp");
            for (j=0; j<Nmstrc; j++)
              xfscanf (in, "%m%lf", &strastre_kwdset, Mstrct[i]+j, Mstrc[i]+j);
          }
        }
      }
      break;
    }
    case eigen_dynamics :
    {
      if (Nlc > 0)
      {
        sprintf(errmsg, "Wrong number of load cases\n"
                        " Number of load cases is %ld\n"  
                        " Number of load cases has to be zero\n"  
                        " you request to solve eigen dynamics.", Nlc);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return (1);
      }
      break;
    }
    case forced_dynamics :
    {
      xfscanf(in, "%k%m", "dload_type", &dynload_kwdset, &Tdload);
      if (Tdload == timeindload)
      {
        Nslc = new long[Nlc];
        memset(Nslc, 0, sizeof(*Nslc)*Nlc);
        Nslc_cum = new long[Nlc];
        memset(Nslc_cum, 0, sizeof(*Nslc_cum)*Nlc);
        // reading number of loading cases
        for (i=0;i<Nlc;i++)
        {
          // resolve static subloadcase which is hidden subloadcase in the internal implementation of mechprep
          // it is referenced just by lc_id WITHOUT given slc_id in the mechprep commands
          // no temperature load type is considered in this static load case
          // temperature load can be defined in the time dependent subloadcase with the help of constant time function
          Tnslc++;
          // resolve time dependent subloadcases
          // it is referenced just by lc_id WITH given slc_id in the mechprep commands
          xfscanf(in, "%k%ld", "lc_id", &lc_id);
          if ((lc_id < 1) && (lc_id > Nlc))
          {
            sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                            " load case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, lc_id, Nlc);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 3;
          }
          lc_id--;
          xfscanf(in, "%k%ld", "num_sublc", Nslc+lc_id);
          if (Nslc[lc_id] < 1)
          {
            sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                            " number of subload cases is < 0", in->line, in->col, in->fname);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 3;
          }
          Tnslc+=Nslc[lc_id];
          if (i<Nlc-1)
            Nslc_cum[2*i+1] = Tnslc;
        }
        //  allocation of array of time functions
        Tf = new gfunct [Tnslc];
        Tlt = new long [Tnslc];
        memset(Tlt, 0, sizeof(*Tlt)*Tnslc);
        for (i=0; i<Tnslc-Nlc; i++)
        {
          xfscanf(in, "%k%ld", "tfunc_lc_id", &lc_id);
          if ((lc_id < 1) || (lc_id > Nlc))
          {
            sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                            " load case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, lc_id, Nlc);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 3;
          }
          lc_id--;
          xfscanf(in, "%k%ld", "tfunc_slc_id", &slc_id);
          if ((slc_id < 1) || (slc_id > Nslc[lc_id]))
          {
            sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                            " subload case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, slc_id, Nslc[lc_id]);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 3;
          }
          slc_id--;
          Tf[Nslc_cum[lc_id]+slc_id+1].read(in);
        }
        // reading of temperature load type for particular subload cases
        Tlt = new long [Tnslc];
        memset(Tlt, 0, sizeof(*Tlt)*Tnslc);
        for(i=0; i<Tnslc-Nlc; i++)
        {
          xfscanf(in, "%k%ld", "tempr_type_lc_id", &lc_id);
          if (lc_id < 1)
          {
            sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                            " load case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, lc_id, Nlc);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 3;
          }
          lc_id--;
          xfscanf(in, "%k%ld", "tempr_type_slc_id", &slc_id);
          if ((slc_id < 1) || (slc_id > Nslc[lc_id]))
          {
            sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                            " subload case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, slc_id, Nslc[lc_id]);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 3;
          }
          // do not increment slc_id at this point because there is the hidden first subloadcase in the load case considereded for the static load
          // and it is considered without temperature load. Type of temperature load can be specified in other nonstatic subloadcases 
          xfscanf(in, "%k%le", "temp_load_type", Tlt+Nslc_cum[lc_id]+slc_id);
          if ((Tlt[Nslc_cum[lc_id]+slc_id] < 0) || (Tlt[Nslc_cum[lc_id]+slc_id] > 3))
          {
            sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                            " temperature load type %ld is out of range <0,3>", in->line, in->col, in->fname, Tlt[Nslc_cum[lc_id]+slc_id]);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 4;
          }
          slc_id--;
        }
      }
      else
      {
        if (tnlc * 2 != Nlc)
        {
          sprintf(errmsg, "Wrong number of load cases\n"
                          " Number of load cases is %ld\n"
                          " Number of load cases has to be odd\n"
                          " because you request to solve forced dynamics.", Nlc);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return (1);
        }
        // reading of temperature load type for particular load cases
        Tlt = new long [Nlc];
        memset(Tlt, 0, sizeof(*Tlt)*Nlc);
        for(i=0; i<Nlc; i++)
        {
          xfscanf(in, "%k%ld", "lc_id", &lc_id);
          if (lc_id < 1)
          {
            sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                            " load case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, lc_id, tnlc);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 3;
          }
          lc_id--;
          xfscanf(in, "%k%le", "temp_load_type", Tlt+lc_id);
          if ((Tlt[lc_id] < 0) || (Tlt[lc_id] > 3))
          {
            sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                            " temperature load type %ld is out of range <0,3>", in->line, in->col, in->fname, Tlt[lc_id]);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 4;
          }
        }
      }
      break;
    }

    case growing_mech_structure:
    case mech_timedependent_prob:
    {
      Nslc = new long[Nlc];
      memset(Nslc, 0, sizeof(*Nslc)*Nlc);
      Nslc_cum = new long[Nlc];
      memset(Nslc_cum, 0, sizeof(*Nslc_cum)*Nlc);
      Tnslc = 0;
      for (i=0;i<Nlc;i++)
      { // reading number of loading cases
        xfscanf(in, "%k%ld", "lc_id", &lc_id);
        if (lc_id < 1)
        {
          sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"  
                          " load case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, lc_id, tnlc);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
        lc_id--;
        xfscanf(in, "%k%ld", "num_sublc", Nslc+lc_id);
        if (Nslc[lc_id] < 0)
        {
          sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                          " number of subload cases is < 0", in->line, in->col, in->fname);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 2;
        }
        Tnslc+=Nslc[lc_id];
        if (i<Nlc-1)
          Nslc_cum[i+1] = Tnslc;
      }
      //  allocation of array of time functions
      Tf = new gfunct [Tnslc];
      for (i=0; i<Tnslc; i++)
      {
        xfscanf(in, "%k%ld", "tfunc_lc_id", &lc_id);
        if ((lc_id < 1) || (lc_id > Nlc))
        {
          sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                          " load case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, lc_id, Nlc);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
        lc_id--;
        xfscanf(in, "%k%ld", "tfunc_slc_id", &slc_id);
        if ((slc_id < 1) || (slc_id > Nslc[lc_id]))
        {
          sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                          " subload case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, slc_id, Nslc[lc_id]);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
        slc_id--;
        Tf[Nslc_cum[lc_id]+slc_id].read(in);
      }
      // reading of temperature load type for particular subload cases and 
      Tlt = new long [Tnslc];
      memset(Tlt, 0, sizeof(*Tlt)*Tnslc);
      if (Mp->tprob == growing_mech_structure)
      {
        Npd = new long [Tnslc];
        memset(Npd, 0, sizeof(*Npd)*Tnslc);
        Spd = new double*[Tnslc];
        memset(Spd, 0, sizeof(*Spd)*Tnslc);
      }
      for(i=0; i<Tnslc; i++)
      {
        xfscanf(in, "%k%ld", "tempr_type_lc_id", &lc_id);
        if (lc_id < 1)
        {
          sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                          " load case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, lc_id, tnlc);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
        lc_id--;
        xfscanf(in, "%k%ld", "tempr_type_slc_id", &slc_id);
        if ((slc_id < 1) || (slc_id > Nslc[lc_id]))
        {
          sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                          " subload case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, slc_id, Nslc[lc_id]);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
        slc_id--;
        xfscanf(in, "%k%ld", "temp_load_type", Tlt+Nslc_cum[lc_id]+slc_id);
        if ((Tlt[Nslc_cum[lc_id]+slc_id] < 0) || (Tlt[Nslc_cum[lc_id]+slc_id] > 3))
        {
          sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                          " temperature load type %ld is out of range <0,3>", in->line, in->col, in->fname, Tlt[Nslc_cum[lc_id]+slc_id]);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 4;
        }
        // in case of growing mechanical problem, the prescribed displacements are read      
        if (Mp->tprob == growing_mech_structure)
        {
          xfscanf(in, "%k%ld", "num_pres_displ_lc_id", &lc_id);
          if (lc_id < 1)
          {
            sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                            " load case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, lc_id, tnlc);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 3;
          }
          lc_id--;
          xfscanf(in, "%k%ld", "num_pres_displ_slc_id", &slc_id);
          if ((slc_id < 1) || (slc_id > Nslc[lc_id]))
          {
            sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                            " subload case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, slc_id, Nslc[lc_id]);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 3;
          }
          slc_id--;
          id = Nslc_cum[lc_id]+slc_id;
          xfscanf(in, "%k%ld", "num_presc_displ", &n);
          if (n < 0)
          {
            sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                            " number of prescribed displacements < 0", in->line, in->col, in->fname);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 5;
          }
          if ((Npd[id] != 0) && (Npd[id] != n))
          {
            sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                            " different numbers of prescribed displacements are required", in->line, in->col, in->fname);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 5;
          }
          Npd[id] = n;
          check = 1;
          if (Spd[id] == NULL) 
          {
            check = 0;
            Spd[id] = new double[Npd[id]];
            memset(Spd[id], 0, sizeof(*Spd[id])*Npd[id]);
          }
          if (Npd[id] > 0)
            xfscanf(in, "%k", "presc_displ_val");
          for (j=0; j<Npd[id]; j++)
          {
            xfscanf(in, "%le", &tmp);
            if (check) 
            {
              if (Spd[id][j] != tmp)
              {
                sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                                " different values of prescribed displacements have been assigned", in->line, in->col, in->fname);
                print_err(errmsg, __FILE__, __LINE__, __func__);
                return 6;
              }
            }
            Spd[id][j] = tmp;
          }
        }
      }
      break;
    }
    case layered_linear_statics :
    {
      xfscanf(in, "%k%ld", "num_layers", &Nlay);
      if (Nlay < 1)
      {
        sprintf(errmsg, "Wrong number of layers (=%ld)", Nlay);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return(2);
      }
      // reading of temperature load type for particular load cases
      Tlt = new long [Nlc];
      memset(Tlt, 0, sizeof(*Tlt)*Nlc);
      for(i=0; i<Nlc; i++)
      {
        xfscanf(in, "%k%ld", "lc_id", &lc_id);
        if (lc_id < 1)
        {
          sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                          " load case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, lc_id, tnlc);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
        lc_id--;
        xfscanf(in, "%k%le", "temp_load_type", Tlt+lc_id);
        if ((Tlt[lc_id] < 0) || (Tlt[lc_id] > 3))
        {
          sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"
                          " temperature load type %ld is out of range <0,3>", in->line, in->col, in->fname, Tlt[lc_id]);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 4;
        }
      }
      break;
    }
    case linear_statics :
    {
      // reading of temperature load type for particular load cases
      Tlt = new long [Nlc];
      memset(Tlt, 0, sizeof(*Tlt)*Nlc);
      for(i=0; i<Nlc; i++)
      {
        xfscanf(in, "%k%ld", "lc_id", &lc_id);
        if (lc_id < 1)
        {
          sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"  
                          " load case id %ld is out of range <1,%ld>", in->line, in->col, in->fname, lc_id, tnlc);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
        lc_id--;
        xfscanf(in, "%k%ld", "temp_load_type", Tlt+lc_id);
        if ((Tlt[lc_id] < 0) || (Tlt[lc_id] > 3))
        {
          sprintf(errmsg, "Error : line=%ld, col=%ld, file %s\n"  
                          " temperature load type %ld is out of range <0,3>", in->line, in->col, in->fname, Tlt[lc_id]);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 4;
        }
        // read macro stress components for homogenization problems
        if (Mp->homog == 3)
        {
          if (Mstrc == NULL)
          {
            Mstrc = new double*[Nlc];
            memset(Mstrc, 0, sizeof(*Mstrc)*Nlc);
          }
          xfscanf(in, "%k%ld", "num_macro_stress_comp", &Nmstrc);
          if (Nmstrc > 0)
          {
            Mstrc[i] = new double[Nmstrc];
            memset(Mstrc[i], 0, sizeof(*Mstrc[i])*Nmstrc);
            xfscanf(in, "%k", "macro_stress_comp");
            for (j=0; j<Nmstrc; j++)
              xfscanf (in, "%lf", Mstrc[i]+j);
          }
        }
        // read macro strain components for homogenization problems
        if (Mp->homog == 4)
        {
          if (Mstrc == NULL)
          {
            Mstrc = new double*[Nlc];
            memset(Mstrc, 0, sizeof(*Mstrc)*Nlc);
          }
          xfscanf(in, "%k%ld", "num_macro_strain_comp", &Nmstrc);
          if (Nmstrc > 0)
          {
            Mstrc[i] = new double[Nmstrc];
            memset(Mstrc[i], 0, sizeof(*Mstrc[i])*Nmstrc);
            xfscanf(in, "%k", "macro_strain_comp");
            for (j=0; j<Nmstrc; j++)
              xfscanf (in, "%lf", Mstrc[i]+j);
          }
        }
        // read macro stress/strain components for homogenization problems
        if (Mp->homog == 9)
        {
          if (Mstrc == NULL)
          {
            Mstrc = new double*[Nlc];
            memset(Mstrc, 0, sizeof(*Mstrc)*Nlc);
            Mstrct = new strastre*[Nlc];
            memset(Mstrct, 0, sizeof(*Mstrct)*Nlc);
          }
          xfscanf(in, "%k%ld", "num_macro_comp", &Nmstrc);
          if (Nmstrc > 0)
          {
            Mstrc[i] = new double[Nmstrc];
            memset(Mstrc[i], 0, sizeof(*Mstrc[i])*Nmstrc);
            Mstrct[i] = new strastre[Nmstrc];
            memset(Mstrct[i], 0, sizeof(*Mstrct[i])*Nmstrc);
            xfscanf(in, "%k", "macro_comp");
            for (j=0; j<Nmstrc; j++)
              xfscanf (in, "%m%lf", &strastre_kwdset, Mstrct[i]+j, Mstrc[i]+j);
          }
        }
      }
      break;
    }
    default :
      break;
  }
  // reading of time functions controlling dofs and switching of elements on/off
  if (Mp->tprob == growing_mech_structure)
  {
    fprintf(stdout, "\n\nReading of section with time functions . . .");
    xf_setsec(in, bsec_str[begsec_loadcase]);
    in->kwdmode = sect_mode_full;
    xfscanf(in, "%k", "time_functions");
    in->kwdmode = sequent_mode;
    Gtm->read_gf (in);
  }
  return 0;
}



/**
  Function reads topology file in, format of the file is described by structure d.

  @param in - pointer to the opened XFILE structure
  @param d  - structure with input data format description


  Function returns:
  @retval 0 : on success
  @retval 4 : in case of unknown mesh format

  In case of reading alternative file format :
  @retval 1 : on error in reading of node
  @retval 2 : on error in reading of element
  @retval 3 : on error in reading of global node numbers
*/
long input_siftop(XFILE *in, descrip *d)
{
  long ret;

  fprintf(stdout, "\n\nReading of mesh topology . . .");
  switch (d->meshfmt)
  {
    case t3d:
      Top->import_t3d(in, d->paral);
      break;
    case sifel:
      ret = Top->read(in, d->paral, d->redgn);
      return(ret);
    default:
      print_err("unknown mesh format is required", __FILE__, __LINE__, __func__);
      return(4);
  }
  return(0);
}



/**
  Function reads hanging nodes data from the separated topology file in.

  @param in - pointer to the opened XFILE structure

  Function returns:
  @retval 0 : on success
  @retval 1 : on error in reading of number of master nodes
  @retval 2 : on error in reading of hanging node number
  @retval 3 : on error in reading of hanging node data
*/
long input_hang_nodes(XFILE *in)
{
  long ret;
  long i, id;

  Numhn = -1;
  fprintf(stdout, "\n\nReading of hanging nodes . . .");
  xfscanf(in, "%ld", &Numhn);
  if ((Numhn < 0) || Numhn > Top->nn)
  {
    print_err("Wrong number of hanging nodes", __FILE__, __LINE__, __func__);
    return 1;
  }
  
  Nod_hang = new hangnode*[Top->nn];
  memset (Nod_hang, 0, sizeof(*Nod_hang)*Top->nn);
  for (i=0; i<Numhn; i++)
  {
    xfscanf(in, "%ld", &id);
    if ((id < 1) || Numhn > Top->nn)
    {
      print_err("Wrong hanging node number", __FILE__, __LINE__, __func__);
      return 2;
    }
    id--;
    Nod_hang[id] = new hangnode;
    ret = Nod_hang[id]->read(in);
    if (ret)
      return 3;
  }
  return(0);
}



/**
  Function reads material database from file with name given by the fname

  @param fname  - string with database file name
  @param d - pointer to structure with description of preprocessor setup 

  Returns :
  @retval 0 - on succes
  @retval 1 - if fails opening file
  @retval 2 - if fails reading material parameters
*/
long input_materials(char *fname, descrip *d)
{
  long ret;
  long mpb = Mespr;
  XFILE *in = NULL;

  fprintf(stdout, "\n\nReading of material database . . .");
  in = xfopen(fname, "rt");
  if (in == NULL)
    return(1);

  in->warning = 1;
  in->kwdmode = sequent_mode;
  in->ignorecase = 1;
  Dbmat = new dbmat;  // allocating new dbmat class
  if (d->matstr == yes)
    // material parameters are read as strings
    ret = Dbmat->read(in);
  else
  {
    // material parameters are read with help of mechmat procedures
    Mespr = 0;
    ret = Dbmat->readmm(in, Mm, d);
    Mespr = mpb;
  }
  xfclose(in);
  if (ret)
    return(2);

  return(0);
}



/**
  Function reads material database from corresponding XFILE section.

  @param in - pointer to the opened XFILE file
  @param d - pointer to structure with description of preprocessor setup 

  Returns :
  @retval 0 - on succes
  @retval 1 - if fails opening file
  @retval 2 - if fails reading material parameters

  Created by TKo, 06.2014
*/
long input_materials(XFILE *in, descrip *d)
{
  long ret = 0;
  long mpb = Mespr;

  fprintf(stdout, "\n\nReading of material database . . .");
  if (xf_setsec(in, bsec_str[begsec_mater]))
    return(1);

  in->warning = 1;
  in->ignorecase = 1;
  in->kwdmode = sect_mode_seq;
  Dbmat = new dbmat;  // allocating new dbmat class
  if (d->matstr == yes)
    // material parameters are read as strings
    ret = Dbmat->read(in);
  else
  {
    // material parameters are read with help of mechmat procedures
    Mespr = 0;
    ret = Dbmat->readmm(in, Mm, d);
    Mespr = mpb;
  }

  if (ret)
    return(2);

  return(0);
}



/**
  Function reads cross-section database from file with name given by the fname.

  @param fname - string with database file name
  @param d - pointer to structure with description of preprocessor setup 

  Returns :
   @retval 0 - on succes
   @retval 1 - if fails opening file
   @retval 2 - if fails reading cross-section parameters

  Created by TKo
*/
long input_crs(char *fname, descrip *d)
{
  long ret;
  long mpb = Mespr;
  XFILE *in = NULL;

  fprintf(stdout, "\n\nReading of cross-section database . . .");
  in = xfopen(fname, "rt");
  if (in == NULL)
    return(1);

  in->warning = 1;
  in->ignorecase = 1;
  in->kwdmode = sequent_mode;
  Dbcrs = new dbcrs; // allocating new dbcrs class
  if (d->crsstr == yes)
    // cross section parameters are read as strings
    ret = Dbcrs->read(in);
  else
  {
    // cross section parameters are read with help of mechcrsec procedures
    Mespr = 0;
    ret = Dbcrs->readmc(in, Mc, d);
    Mespr = mpb;
  }
  xfclose(in);
  if (ret)
    return(2);

  return(0);
}



/**
  Function reads cross-section database from file with name given by the fname.

  @param in - pointer to the opened XFILE file
  @param d - pointer to structure with description of preprocessor setup 

  Returns :
   @retval 0 - on succes
   @retval 1 - if fails opening file
   @retval 2 - if fails reading cross-section parameters

  Created by TKo, 06.2014
*/
long input_crs(XFILE *in, descrip *d)
{
  long ret;
  long mpb = Mespr;

  fprintf(stdout, "\n\nReading of cross-section database . . .");
  if (xf_setsec(in, bsec_str[begsec_crsec]))
    return(1);

  in->warning = 1;
  in->kwdmode = sect_mode_seq;
  in->ignorecase = 1;
  Dbcrs = new dbcrs; // allocating new dbcrs class
  if (d->crsstr == yes)
    // cross section parameters are read as strings
    ret = Dbcrs->read(in);
  else
  {
    // cross section parameters are read with help of mechcrsec procedures
    Mespr = 0;
    ret = Dbcrs->readmc(in, Mc, d);
    Mespr = mpb;
  }
  if (ret)
    return(2);

  return(0);
}



/**
  Function reads nodal properties from file in
 
  @param in  - poinetr to opened input file with property description

  Returns :
  @retval 0  - on succes
  @retval 1  - fails input of ndofs
  @retval 2  - fails input of boundary conditions
  @retval 3  - fails input of common code numbers
  @retval 4  - fails input of springs at nodes
  @retval 5  - fails input of time functions of dofs
  @retval 6  - fails input of cross-sections
  @retval 7  - fails input of local coordinate systems
  @retval 8  - fails input of load
  @retval 9  - fails input of time dependent load
  @retval 10 - fails input of initial conditions
  @retval 11 - fails input of temperatures

  Created 04.2008 by Tomas Koudelka koudelka@cml.fsv.cvut.cz
*/
long input_nodprop(XFILE *in)
{
  long err;
  const enumstr nodsects[] = {{"begsec_nodvertpr",3}, {"begsec_nodedgpr",4}, 
                        {"begsec_nodsurfpr", 5}, {"begsec_nodvolpr",6}};
  long nsect = sizeof(nodsects)/sizeof(*nodsects);
  
  Nod_ccn = NULL;

  fprintf(stdout, "\n\nReading of nodal properties . . .");
  in->kwdmode = sequent_mode;
  // reading of nodal dofs
  err = input_nod_ndof(in, nodsects, nsect);
  if (err)
    return 1;
  // reading of boundary conditions at nodes 
  // (i.e. supports or prescribed displacements)
  err = input_nod_bocon(in, nodsects, nsect);
  if (err)
    return 2;
  // reading of coupled dofs at nodes
  err = input_nod_coupl_dofs(in, nodsects, nsect);
  if (err)
    return 3;
  // reading of springs at nodes
  err = input_nod_springs(in, nodsects, nsect);
  if (err)
    return 4;
  // reading of time functions of dofs
  err = input_nod_dof_tfunc(in, nodsects, nsect);
  if (err)
    return 5;
  // reading of cross-sections at nodes
  err = input_nod_crsec(in, nodsects, nsect);   
  if (err)
    return 6;
  // reading of local coordinate system
  err = input_nod_lcs(in, nodsects, nsect);
  if (err)
    return 7;
  // reading of nodal load
  err = input_nod_load(in, nodsects, nsect);   
  if (err)
    return 8;
  // reading of nodal time dependent load
  err = input_nod_tdload(in, nodsects, nsect);
  if (err)
    return 9;
  // reading of initial conditions
  err = input_nod_initcond(in, nodsects, nsect);
  if (err)
    return 10;
  // reading of initial displacements due to rotation about the given axis
  err = input_nod_rotinidispl(in, nodsects, nsect);
  if (err)
    return 11;
  // reading of nodal temperatures  
  err = input_nod_temper(in, nodsects, nsect);
  if (err)
    return 12;
  // reading of periodic nodal boundary condition
  err = input_nod_periodic_bc(in, nodsects, nsect-1);
  if (err)
    return 13;
  err = input_nod_genhn(in, nodsects, nsect);
  if (err)
    return 14;
  return 0;
}



/**
  Function reads element properties from file in
 
  @param in  - poinetr to opened input file with property description

  Returns :
  @retval 0  - on succes
  @retval 1  - fails input of element types
  @retval 2  - fails input of material
  @retval 3  - fails input of cross-sections
  @retval 4  - fails input of local coordinate systems
  @retval 5  - fails input of load record
  @retval 6  - fails input of edge load
  @retval 7  - fails input of surface load
  @retval 8  - fails input of volume load
  @retval 9  - fails input of eigenstrains
  @retval 10 - fails input of time functions for elements

  Created 04.2008 by Tomas Koudelka koudelka@cml.fsv.cvut.cz
*/
long input_elemprop(XFILE *in)
{
  long err;

  // all sections with element properties
  const enumstr elemsects[] = {{bsec_str[begsec_eledgpr].alias, bsec_str[begsec_eledgpr].id},
                                {bsec_str[begsec_elsurfpr].alias, bsec_str[begsec_elsurfpr].id},
                                {bsec_str[begsec_elvolpr].alias, bsec_str[begsec_elvolpr].id}};
  long nsect = sizeof(elemsects)/sizeof(*elemsects);

  // sections with edge properties for elements
  const enumstr edgesect[] = {{bsec_str[begsec_eledgpr].alias, bsec_str[begsec_eledgpr].id}};
  long nedgesect = sizeof(edgesect)/sizeof(*edgesect);

  // sections with surface properties for elements
  const enumstr surfsect[] = {{bsec_str[begsec_elsurfpr].alias, bsec_str[begsec_elsurfpr].id}};
  long nsurfsect = sizeof(surfsect)/sizeof(*surfsect);

  // sections with volume properties for elements
  //  const enumstr volsect[]  = {{bsec_str[begsec_elvolpr].alias, bsec_str[begsec_elvolpr].id}};
  //  long nvolsect = sizeof(volsect)/sizeof(*volsect);

  fprintf(stdout, "\n\nReading of element properties . . .");
  in->kwdmode = sequent_mode;
  // reading of element type
  err = input_elem_type(in, elemsects, nsect);
  if (err)
    return 1;
  Mt->alloc_prep(Top->nn, Top->ne, El_type);
  // reading of element material
  err = input_elem_mat(in, elemsects, nsect);
  if (err)
    return 2;
  // reading of element cross-section
  err = input_elem_crsec(in, elemsects, nsect);   
  if (err)
    return 3;
  // reading of local coordinate system for integartion points
  err = input_elem_lcs(in, elemsects, nsect);
  if (err)
    return 4;
  // reading of load records assigned to elements
  err = input_elem_load(in, elemsects, nsect);   
  if (err)
    return 5;
  // reading of edge loads assigned to elements
  err = input_elem_loadedge(in, edgesect, nedgesect);
  if (err)
    return 6;
  // reading of surface loads assigned to elements
  err = input_elem_loadsurf(in, surfsect, nsurfsect);
  if (err)
    return 7;
  // reading of volume loads assigned to elements
  //  err = input_elem_loadvol(in, volsect, nvolsect);
  err = input_elem_loadvol(in, elemsects, nsect);
  if (err)
    return 8;
  // reading of eigenstrains assigned to elements
  err = input_elem_eigstr(in, elemsects, nsect);
  if (err)
    return 9;
  // reading of element switching time function
  err = input_elem_eltimefunc(in, elemsects, nsect);
  if (err)
    return 10;
  // reading of time dependent element load
  err = input_elem_tdload(in, elemsects, nsect);
  if (err)
    return 11;
  // reading of edge loads assigned to elements
  err = input_elem_tdloadedge(in, edgesect, nedgesect);
  if (err)
    return 12;
  // reading of surface loads assigned to elements
  err = input_elem_tdloadsurf(in, surfsect, nsurfsect);
  if (err)
    return 13;
  // reading of volume loads assigned to elements
  err = input_elem_tdloadvol(in, elemsects, nsect);
  if (err)
    return 14;
  return 0;
}



/**
  The function assign number of dofs at each node. It scans 
  sections nodvertpr, nodedgpr, nodsurfpr, nodvolpr for keyword "ndofn".
  The record which assigns dofs to node looks like follows:
  "ndofn" ndof "propid" prop
  where ndof and prop are positive integer numbers.
  
  The asssigning starts at section for volumes and finishes at
  section for vertices. If a node belongs to several entities and dofs are
  prescribed for these entities, the assigned dofs are rewritten 
  and in case of multiple assigning, the message is written into the log file.
  Finally, the function checks whether all nodes have assigned ndof
 
  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - wrong number of dofs has been assigned to nodes
  @retval 2 - different numbers of dofs have been assigned to node 
  @retval 3 - no number of dofs has been assigned to nodes
  @retval 4 - nodes with required property and entity type cannot be found
*/
long input_nod_ndof(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->nn), col(Top->nn);
  long aline, acol;
  long prop, ndof;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k;
  long prop_used;

  fprintf(stdout, "\n reading of number of dofs at nodes");
  Nod_ndof = new long[Top->nn];
  memset(Nod_ndof, 0, sizeof(*Nod_ndof)*Top->nn);

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "ndofn"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "ndofn", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "ndofn", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%ld%k%ld", "ndofn", &ndof, "propid", &prop);
      in->kwdmode = bkwdmode;
      if (ndof < 1)
      {
        sprintf(errmsg, "Number of assigned dof at node should be > 0\n"
                        "see input file %s, line %ld, column %ld", in->fname, aline, acol);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 1;
      }
      prop_used = 0;
      // assign given ndof to nodes with given property id
      for(k=0; k<Top->nn; k++)
      {    
        if (Top->nodes[k].searchprop(prop, ent, k, Top->edges, Top->surfaces))
        {
          if (Nod_ndof[k] != 0)
            fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already assigned ndofn at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
          if ((Nod_ndof[k] != 0) && (Nod_ndof[k] != ndof))
          {
            sprintf(errmsg, "Different number of dofs is assigned at node %ld\n"
                    "see input file %s, line %ld, column %ld", k+1, in->fname, aline, acol);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 2;
          }
          Nod_ndof[k] = ndof;
          prop_used = 1;
          // backup of line and column for eventual log message
          line[k]     = aline;
          col[k]      = acol;
        }        
      }
      if (Check_unused_prop && (prop_used == 0))
      {
        // no nodes were found with the given property and checking of unused properties is off
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Number of nodal dofs (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 4;
      }
    } // end of nkwd loop
  }   // end of section loop

  // check whether all nodes have assigned a ndof
  for(i=0; i<Top->nn; i++)
  {
    if (Nod_ndof[i] == 0)
    {
      sprintf(errmsg, "Node %ld has not assigned number of dof", i+1);
      print_err(errmsg, __FILE__, __LINE__, __func__);
      return 3;
    }
  }

  return 0;
}



/**
  Function reads and assigns boundary conditions at each node.
  It scans sections nodvertpr, nodedgpr, nodsurfpr, nodvolpr for keyword "bocon".
  The record which assigns dofs to node looks like follows:
  "bocon" "propid" prop "num_bc" nbc {"dir" dir "cond" cond [func]}[nbc] [lc_id nlc [slcid slc] [expr]]}xnbc

  The asssigning starts at section for volumes and finishes at
  section for vertices. If a node belongs to several entities and boundary conditions are
  prescribed for these entities, the assigned boundary conditions are merged and a message 
  is written into the log file. In case of multiple assigning of different boundary conditions 
  in same direction, error of merging is signalized.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - assignment of different number of dofs at one node
  @retval 2 - expression could not be parsed
  @retval 3 - two boundary conditons cannot be merged
  @retval 4 - growing mechanical problem type is required (dofs have to be controlled by time functions)
  @retval 5 - nodes with required property and entity type cannot be found
*/
long input_nod_bocon(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->nn), col(Top->nn);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  ivector setnodes(Top->nn);
  bocon tbc, *aux;
  long i, j, k;
  long prop, ndof;
  long errcode;

  Nod_bocon = new bocon*[Top->nn];
  memset(Nod_bocon, 0, sizeof(*Nod_bocon)*Top->nn);

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of nodal boundary conditions");
  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "bocon"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "bocon", 1);
    if (nkwd && (Mp->tprob == growing_mech_structure))
    {
      print_err("Problem type 'growing_mech_structure' is required.\n"
                " Nodal dofs have to be controlled by time functions.\n"
                " Use 'nod_tfunc' preprocessor keyword for dof control.\n"
                " Invalid keyword 'bocon' found at line %ld, column %ld, file %s",
                __FILE__, __LINE__, __func__, in->line, in->col, in->fname);
      return 4;
    }
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "bocon", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "bocon", "propid", &prop);
      in->kwdmode = bkwdmode;
      // investigation of ndof for nodes with property prop of entity ent
      ndof = Top->get_ndofn(prop, ent, Nod_ndof, setnodes);
      if (ndof == 0)
      {
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s\n"
                        " have different number of dofs. Boundary condition\n" 
                        " (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 1;
      }
      if (ndof == -1)
      {
        if (Check_unused_prop == 0)
        // no nodes were found with the given property and checking of unused properties is off
          continue;
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Boundary condition (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 5;
      }
      // rest of record of boundary condition is read
      errcode = tbc.read(in, ndof);
      if (errcode)
      {
        if (errcode < 2)
        {
          sprintf(errmsg, "Boundary condition  or direction is invalid\n"
                          " in file %s, line %ld, col %ld", in->fname, aline, acol);
        }
        else
        {
          if (errcode < 5)
          {
            sprintf(errmsg, "Wrong (sub)load case id\n"
                            " in file %s, line %ld, col %ld", in->fname, aline, acol);
          }
          else
            sprintf(errmsg, "Expression for boundary condition could not be parsed\n"
                            " in file %s, line %ld, col %ld", in->fname, aline, acol);
        }
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 2;
      }
      // assigning of boundary condition
      for (k=0; k<Top->nn; k++)
      {        
        if (setnodes[k] < 0) // node has not required property
          continue;
        if (Nod_bocon[k])
        {
          fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already assigned boundary condition at line %ld, col %ld, file: %s\n", 
                  aline, acol, in->fname, k+1, line[k], col[k], in->fname);
          aux = Nod_bocon[k];
          Nod_bocon[k] = tbc.merge(Nod_bocon[k]);
          delete aux;
          if (Nod_bocon[k] == NULL)
          {
            sprintf(errmsg, "Conflict in merging of boundary condition at line %ld, col %ld\n"
                            "with previous boundary condition at line %ld, col %ld in file %s", aline, acol, line[k], col[k], in->fname);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 3;
          }
        }
        else
          Nod_bocon[k] = tbc.copy();

        // backup of line and column for eventual log message
        line[k]      = aline;
        col[k]       = acol;
      }  // end of loop for nodes
    }    // end of nkwd loop 
  }      // end of section loop
  return 0;
}



/**
  Function reads and assigns coupled dofs at each node.
  It scans sections nodvertpr, nodedgpr, nodsurfpr, nodvolpr for keyword "dof_coupl".
  The record which assigns coupled dofs to nodes looks like follows:
  "dof_coupl" "propid" prop "num_coup_dir" ndir {"dir" dir_index}[ndir]

  The asssigning starts at section for volumes and finishes at
  section for vertices. If a node belongs to several entities and boundary conditions are
  prescribed for these entities, the assigned boundary conditions are merged and a message 
  is written into the log file.In case of multiple assigning of different boundary conditions 
  in same direction, error of merging is signalized.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - assignment of different number of dofs at one node
  @retval 2 - invalid number of coupled dofs (out of <1;ndof>)  
  @retval 3 - invalid coupled dof number (out of <1;ndof>)
  @retval 4 - growing mechanical problem type is required (dofs have to be controlled by time functions)
  @retval 5 - nodes with required property and entity type cannot be found
*/
long input_nod_coupl_dofs(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->nn), col(Top->nn);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  ivector setnodes(Top->nn);
  long i, j, k, l;
  long prop, ndof, dir, ndir;
  long ccn, accn;

  Nod_ccn = NULL;
  ccn = accn = 1;

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;
  fprintf(stdout, "\n reading of coupled dofs at nodes");

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "dof_coupl"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "dof_coupl", 1);
    if (nkwd && (Mp->tprob == growing_mech_structure))
    {
      print_err("Problem type 'growing_mech_structure' is required.\n"
                " Nodal dofs have to be controlled by time functions.\n"
                " Use 'nod_tfunc' preprocessor keyword for dof control.\n"
                " Invalid keyword 'dof_coupl' found at line %ld, column %ld, file %s",
                __FILE__, __LINE__, __func__, in->line, in->col, in->fname);
      return 4;
    }
    if (Nod_ccn == NULL)
    {
      Nod_ccn = new long*[Top->nn];
      memset(Nod_ccn, 0, sizeof(*Nod_ccn)*Top->nn);
    }
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "dof_coupl", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "dof_coupl", "propid", &prop);
      in->kwdmode = bkwdmode;
      // investigation of ndof for nodes with property prop of entity ent
      ndof = Top->get_ndofn(prop, ent, Nod_ndof, setnodes);
      if (ndof == 0)
      {
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s\n"
                        " have different number of dofs. Coupled dofs \n" 
                        " (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 1;
      }
      if (ndof == -1)
      {
        if (Check_unused_prop == 0)
        // no nodes were found with the given property and checking of unused properties is off
          continue;
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Coupled dofs (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 5;
      }
      xfscanf(in, "%k%ld", "ndir", &ndir);
      if ((ndir < 1) || (ndir > ndof))
      {
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s\n"
                        " have number of coupled dof dircetions out of range <1,%ld>" 
                        " (line %ld, column %ld, file %s)", 
                prop, entitypstr[ent-1].alias, ndof, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 2;
      }
      for (k = 0; k < ndir; k++)
      {
        xfscanf(in, "%k%ld", "dir", &dir);
        if ((dir < 1) || (dir > ndof))
        {
          sprintf(errmsg, "Nodes with property %ld belonging to entity %s\n"
                          " have coupled dof number out of range <1,%ld>" 
                          " (line %ld, column %ld, file %s)", 
                  prop, entitypstr[ent-1].alias, ndof, aline, acol, in->fname);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
        ccn++;
        accn = ccn;
        // assigning of coupled dofs
        for (l=0; l<Top->nn; l++)
        {        
          // investigation for previously assigned common code numbers
          if (setnodes[l] < 0) // node has not required property
            continue;
          if (Nod_ccn[l] == NULL)
            continue;
          if (Nod_ccn[l][dir-1] != 0)
          {
            accn = Nod_ccn[l][dir-1];
            ccn--;
            fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already had coupled dofs at line %ld, col %ld, file: %s\n", aline, acol, in->fname, l+1, line[l], col[l], in->fname);
            break;
          }
        }           
        for (l=0; l<Top->nn; l++)
        {        
          if (setnodes[l] < 0) // node has not required property
            continue;
          if (Nod_ccn[l] == NULL)
          {
            Nod_ccn[l] = new long[ndof];
            memset(Nod_ccn[l], 0, sizeof(*Nod_ccn[l])*ndof);
            // backup of line and column for eventual log message
            line[l] = aline;
            col[l]  = acol;
          }
          Nod_ccn[l][dir-1] = accn;
        } // end of loop for nodes
      }   // end of loop for ndir
    }  // end of nkwd loop 
  }    // end of section loop 
  return 0;
}



/**
  Function reads time functions for switching dofs on/off at each node.
  It scans sections nodvertpr, nodedgpr, nodsurfpr, nodvolpr for keyword "nod_tfunc".
  The record which assigns dofs to node looks like follows:
  "nod_tfunc" "propid" prop {"tfunc_id" tfunc_id}[ndofn]

  The asssigning starts at section for volumes and finishes at
  section for vertices. In case of multiple assigning of different time functions,
  error of assignment is signalized.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - cannot determine unique number of dofs for the given property and entity type
  @retval 2 - number of prescribed time functions is out of range <1,ndofn>
  @retval 3 - dof id is out of range <1,ndofn>
  @retval 4 - nodes with required property and entity type cannot be found
  @retval 5 - invalid time function id

  Created by TKo, 09.2009
*/
long input_nod_dof_tfunc(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->nn), col(Top->nn);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  ivector setnodes(Top->nn);
  long tf_id;
  long i, j, k, l;
  long prop, ndof, ndir, dir;

  if (Nod_ccn == NULL)
  {
    Nod_ccn = new long*[Top->nn];
    memset(Nod_ccn, 0, sizeof(*Nod_ccn)*Top->nn);
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of dof time functions at nodes");

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "nod_tfunc"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "nod_tfunc", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "nod_tfunc", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "nod_tfunc", "propid", &prop);
      in->kwdmode = bkwdmode;
      // investigation of ndof for nodes with property prop of entity ent
      ndof = Top->get_ndofn(prop, ent, Nod_ndof, setnodes);
      if (ndof == 0)
      {
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s\n"
                        " have different number of dofs. Nodal time functions \n" 
                        " (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 1;
      }
      if (ndof == -1)
      {
        if (Check_unused_prop == 0)
        // no nodes were found with the given property and checking of unused properties is off
          continue;
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Nodal time functions (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 4;
      }
      xfscanf(in, "%k%ld", "ndir", &ndir);
      if ((ndir < 1) || (ndir > ndof))
      {
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s\n"
                        " have number of dof dircetions out of range <1,%ld>" 
                        " (line %ld, column %ld, file %s)", 
                prop, entitypstr[ent-1].alias, ndof, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 2;
      }
      for (k = 0; k < ndir; k++)
      {
        xfscanf(in, "%k%ld", "dir", &dir);
        if ((dir < 1) || (dir > ndof))
        {
          sprintf(errmsg, "Nodes with property %ld belonging to entity %s\n"
                          " have time function dof number out of range <1,%ld>" 
                          " (line %ld, column %ld, file %s)", 
                  prop, entitypstr[ent-1].alias, ndof, aline, acol, in->fname);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
        xfscanf(in, "%k%ld", "tfunc_id", &tf_id);
        if ((tf_id < 1) && (tf_id > Gtm->ngf))
        {
          sprintf(errmsg, "Invalid time function index");
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 5;
        }
        // assigning of a time function index to the given dof
        for (l=0; l<Top->nn; l++)
        {        
          if (setnodes[l] < 0) // node has not required property
            continue;
          if (Nod_ccn[l] == NULL)
          { // no time function has been assigned to the given node
            Nod_ccn[l] = new long[ndof];
            memset(Nod_ccn[l], 0, sizeof(*Nod_ccn[l])*ndof);
            Nod_ccn[l][dir-1] = tf_id;
            // backup of line and column for eventual log message
            line[l] = aline;
            col[l]  = acol;
          }
          else
          { // The node had already assigned time function which will be rewritten with actual one.
            // This concept is necessary due to effective assignment 'default' time function for the model body
            // and different time functions for the selected surfaces|edges|vertices.
            if ((Nod_ccn[l][k] != 0) && (Nod_ccn[l][k] != tf_id))
            {
              fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already had assigned different "
                           "time function at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              l+1, line[l], col[l], in->fname);
            }
            Nod_ccn[l][dir-1] = tf_id;
          }
        }  // end of loop for nodes
      }    // end of loop for ndir
    } // end of nkwd loop
  }   // end of section loop
  return 0;
}



/**
  Function reads cross-sections at each node.
  It scans sections nodvertpr, nodedgpr, nodsurfpr, nodvolpr for keyword "nod_crsec".
  The record which assigns a cross-section to node looks like follows:
  "nod_crsec" "propid" prop "type" type_keyword "type_id" type_index

  The asssigning starts at section for volumes and finishes at
  section for vertices. In case of multiple assigning of different cross-sections,
  error of assignment is signalized.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - cannot find required cross-section in database
  @retval 2 - different cross-section type has been assigned at one node
  @retval 3 - nodes with required property and entity type cannot be found
*/
long input_nod_crsec(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->nn), col(Top->nn);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  ivector setnodes(Top->nn);
  crsectype cst;
  long      csti;
  long i, j, k, ic;
  long prop;
  long prop_used;

  Nod_cst = new crsectype[Top->nn];
  memset(Nod_cst, 0, sizeof(*Nod_cst)*Top->nn);
  Nod_csti = new long[Top->nn];
  memset(Nod_csti, 0, sizeof(*Nod_csti)*Top->nn);
  Nod_cstdbi = new long[Top->nn];
  memset(Nod_cstdbi, 0, sizeof(*Nod_cstdbi)*Top->nn);

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of cross-sections at nodes");

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "nod_crsec"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "nod_crsec", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "nod_crsec", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      // read ndof and property id
      cst = csplanestr;
      csti = 2;
      in->kwdmode = sect_mode_seq;
      xfscanf(in, "%k%k%ld", "nod_crsec", "propid", &prop);
      xfscanf(in, "%k%m%k%ld", "type", &crsectype_kwdset, &cst, "type_id", &csti);
      in->kwdmode = bkwdmode;
      // searching of prescribed cross-section type and index in the cross-section database
      ic = Dbcrs->search_crs(cst, csti-1);
      if (ic < 0)
      {
        sprintf(errmsg, "Line %ld, col %ld, file %s: cannot find required cross-section in the database", in->line, in->col, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 1;
      }
      else
        Dbcrs->mark_used(ic, csti-1);
          
      Top->get_propent_nodes(prop, ent, setnodes);
      prop_used = 0;
      // assigning of a time function index to the given dof
      for (k=0; k<Top->nn; k++)
      {        
        if (setnodes[k] < 0) // node has not required property
          continue;
        if ((Nod_cst[k] != nocrosssection) && ((Nod_cst[k] != cst) || (Nod_csti[k] != csti)))        
        {  // node has assigned different cross-section types
          sprintf(errmsg, "Line %ld, col %ld, file %s: Node %ld has already had assigned different"
                          "cross-section type at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                          k+1, line[k], col[k], in->fname);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 2;
        }
        Nod_cst[k]  = cst;
        Nod_csti[k] = csti;
        Nod_cstdbi[k] = ic;
        prop_used = 1;
        // backup of line and column for eventual log message
        line[k] = aline;
        col[k]  = acol;
      }
      if (Check_unused_prop && (prop_used == 0))
      {
        // no nodes were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Nodal cross-section type (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 3;
      }
    }  // end of nkwd loop
  }    // end of section loop
  return 0;
}



/**
  Function reads spring supports at each node.
  It scans sections nodvertpr, nodedgpr, nodsurfpr, nodvolpr for keyword "nod_crsec".
  The record which assigns a springs to node looks like follows:
  "nod_spring" "propid" prop "dir" dir_index "num_mat" nmat nmat*{"type" material_type "type_id" material_id}

  The asssigning starts at section for volumes and finishes at
  section for vertices. In case of multiple assigning of different springs,
  error of assignment is signalized.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - cannot find required material in database
  @retval 2 - different material type has been assigned to spring at one node
  @retval 3 - nodes with required property and entity type cannot be found
  @retval 4 - direction of spring is out of range <1,6> or <1, ndofn>
  @retval 5 - nodes with the given property have different number of dofs 
*/
long input_nod_springs(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->nn), col(Top->nn);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  ivector setnodes(Top->nn);
  long i, j, k, l;
  mattype *tmp_type;
  long *tmp_id, *ic;
  long prop, dir, nmat, ndof, err;
  long prop_used;

  Nod_nsprmat = new long*[Top->nn];
  memset(Nod_nsprmat, 0, sizeof(*Nod_nsprmat)*Top->nn);
  Nod_sprmattype = new mattype**[Top->nn];
  memset(Nod_sprmattype, 0, sizeof(*Nod_sprmattype)*Top->nn);
  Nod_sprmatid = new long**[Top->nn];
  memset(Nod_sprmatid, 0, sizeof(*Nod_sprmatid)*Top->nn);
  Nod_sprmatdbi = new long**[Top->nn];
  memset(Nod_sprmatdbi, 0, sizeof(*Nod_sprmatdbi)*Top->nn);    

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of springs at nodes");

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "nod_crsec"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "nod_spring", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "nod_spring", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read property id, dir and num_mat
      xfscanf(in, "%k%k%ld%k%ld%k%ld", "nod_spring", "propid", &prop, "dir", &dir, "num_mat", &nmat);
      dir--;

      if (nmat<=0)
      {
        print_err("number of material types is <= 0", __FILE__, __LINE__, __func__);
        return 1;
      }
      if ((dir < 0) || (dir > 5))
      {
        print_err("direction of spring is out of range <1,6>", __FILE__, __LINE__, __func__);
        return 4;
      }

      tmp_type = new mattype[nmat];
      memset(tmp_type, 0, sizeof(*tmp_type)*nmat);
      tmp_id   = new long[nmat];
      memset(tmp_id, 0, sizeof(*tmp_id)*nmat);
      ic       = new long[nmat];
      memset(ic, 0, sizeof(*ic)*nmat);

      for (k=0; k<nmat; k++)
      {
        xfscanf(in, "%k%m%k%ld", "type", &mattype_kwdset, tmp_type+k, "type_id", tmp_id+k);
        // searching of prescribed material type and index in the material database
        ic[k] = Dbmat->search_mat(tmp_type[k], tmp_id[k]-1);
        if (ic[k] < 0)
        {
          sprintf(errmsg, "Line %ld, col %ld, file %s: cannot find %ld. required material in the database", in->line, in->col, in->fname, k+1);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          delete [] tmp_type;
          delete [] tmp_id;
          delete [] ic;
          return 1;
        }
        else
          Dbmat->mark_used(ic[k], tmp_id[k]-1);
      }
      in->kwdmode = bkwdmode;
          
      Top->get_propent_nodes(prop, ent, setnodes);
      prop_used = 0;
      // assigning of a time function index to the given dof
      for (k=0; k<Top->nn; k++)
      {        
        if (setnodes[k] < 0) // node has not required property
          continue;

        ndof = Top->get_ndofn(prop, ent, Nod_ndof, setnodes);
        if (ndof == 0)
        {
          sprintf(errmsg, "Nodes with property %ld belonging to entity %s\n"
                          " have different number of dofs. Nodal spring \n" 
                          " (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                          prop, entitypstr[ent-1].alias, aline, acol, in->fname);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          delete [] tmp_type;
          delete [] tmp_id;
          delete [] ic;
          return 5;
        }

        if (Nod_nsprmat[k] == NULL)
        // springs have not been assigned yet
        {
          Nod_nsprmat[k] = new long[ndof];
          memset(Nod_nsprmat[k], 0, sizeof(*Nod_nsprmat[k])*ndof);
          Nod_sprmattype[k] = new mattype*[ndof];
          memset(Nod_sprmattype[k], 0, sizeof(*Nod_sprmattype[k])*ndof);
          Nod_sprmatid[k] = new long*[ndof];
          memset(Nod_sprmatid[k], 0, sizeof(*Nod_sprmatid[k])*ndof);
          Nod_sprmatdbi[k] = new long*[ndof];
          memset(Nod_sprmatdbi[k], 0, sizeof(*Nod_sprmatdbi[k])*ndof);    
        }

        
        prop_used = 1;
        if (Nod_nsprmat[k][dir] != 0)
        {// check for multiple assignment of different material types to one element
          err = 0;
          if (nmat != Nod_nsprmat[k][dir])
            err = 1;
          for (l=0; (l<nmat) && (err==0); l++)
          {
            if ((tmp_type[l] != Nod_sprmattype[k][dir][l]) || (tmp_id[l] != Nod_sprmatid[k][dir][l]))
                err = 1;
          }        
          if (err)           
          {
            sprintf(errmsg, "Different material types have been assigned to spring at node %ld\n"
                    "see input file %s, line %ld, column %ld", k+1, in->fname, aline, acol);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            delete [] tmp_type;
            delete [] tmp_id;
            delete [] ic;
            return 2;
          }
          else
            fprintf(Log, "Line %ld, col %ld, file %s: Spring at node %ld has already assigned material type at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
        }
        else
        {
          Nod_sprmattype[k][dir] = new mattype[nmat];
          Nod_sprmatid[k][dir]   = new long[nmat];
          Nod_sprmatdbi[k][dir]  = new long[nmat];
          Nod_nsprmat[k][dir]    = nmat;
          for(l=0; l<nmat; l++)
          {
            Nod_sprmattype[k][dir][l] = tmp_type[l];
            Nod_sprmatid[k][dir][l]   = tmp_id[l];
            Nod_sprmatdbi[k][dir][l]  = ic[l];
          }
          // backup of line and column for eventual log message
          line[k]     = aline;
          col[k]      = acol;
        }
      } // end of loop for all nodes

      delete [] tmp_type;
      delete [] tmp_id;
      delete [] ic;

      if (Check_unused_prop && (prop_used == 0))
      {
        // no nodes were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Nodal spring (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 3;
      }
    }  // end of nkwd loop
  }    // end of section loop
  return 0;
}



/**
  Function reads local coordinate systems defined at nodes.
  It scans sections nodvertpr, nodedgpr, nodsurfpr, nodvolpr for keyword "nod_lcs".
  The record which assigns local coordinate system to node looks like follows:
  "nod_lcs" "propid" prop "dim" vect_dim {"basevec" base_vector_components}[vect_dim]

  The asssigning starts at section for volumes and finishes at
  section for vertices. In case of multiple assigning of different local coordinate systems,
  error of assignment is signalized.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - invalid dimension of base vectors
  @retval 2 - nodes with required property and entity type cannot be found
*/
long input_nod_lcs(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->nn), col(Top->nn);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  ivector setnodes(Top->nn);
  long dim;
  long i, j, k, l, m;
  long prop, check;
  vector *basev = NULL;
  long prop_used;

  Nod_lcs = new vector*[Top->nn];
  memset(Nod_lcs, 0, sizeof(*Nod_lcs)*Top->nn);

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of nodal local coordinate systems");

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "nod_lcs"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "nod_lcs", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "nod_lcs", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "nod_lcs", "propid", &prop);
      xfscanf(in, "%k%ld", "dim", &dim);
      if ((dim < 2) || (dim > 3))
      {
        sprintf(errmsg, "Line %ld, col %ld, file %s: Wrong dimension of base vectors (%ld)\n"
                "Dimension has to be in range  <2,3>\n", 
                aline, acol, in->fname, dim);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        delete [] basev;
        return 1;
      }
      // reading of base vectors
      basev = new vector[dim];
      for (k=0; k<dim; k++)
      {
        allocv(dim, basev[k]);
        xfscanf(in, "%k", "basevec");
        readv(in, basev[k]);
      }      
      in->kwdmode = bkwdmode;

      // assigning of lcs to the given nodes
      Top->get_propent_nodes(prop, ent, setnodes);
      prop_used = 0;
      for (k=0; k<Top->nn; k++)
      {        
        if (setnodes[k] < 0) // node has not required property
          continue;
        if (Nod_lcs[k] != NULL)
        { 
          check = 1; 
          for (l=0; l<dim; l++)
          {
            for(m=0; m<dim; m++)
            {
              if (Nod_lcs[k][l][m] != basev[l][m])
              {
                check = 0;
                l = dim;
                break;
              }
            }
          }
          if (check == 0)
          {
            // node has assigned different local coordinet systems
            sprintf(errmsg, "Line %ld, col %ld, file %s: Node %ld has already had assigned different"
                            "local coordinate system at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                            k+1, line[k], col[k], in->fname);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            delete [] basev;
            return 1;
          }
        }
        Nod_lcs[k]  = new vector[dim];
        for (l=0; l<dim; l++)
        {
          allocv(dim, Nod_lcs[k][l]);
          copyv(basev[l], Nod_lcs[k][l]);
        } 
        prop_used = 1;
        // backup of line and column for eventual log message
        line[k] = aline;
        col[k]  = acol;
      }
      delete [] basev;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no nodes were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Nodal lcs (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 2;
      }
    } // end of nkwd loop
  }   // end of section loop
  return 0;
}



/**
  Function reads nodal load. It scans sections nodvertpr, nodedgpr, nodsurfpr, 
  nodvolpr for keyword "nod_load". 
  The record which assigns nodal load looks like follows:
  "nod_load" "propid" prop "lc_id" lcid ["slc_id" slcid] "load_comp" {load_components}[ndofn]

  The asssigning starts at section for volumes and finishes at
  section for vertices. If a node belongs to several entities and loads are
  prescribed for these entities, the assigned loads are merged and a message 
  is written into the log file.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  Returns :
  @retval 0 - on succes 
  @retval 1 - load case index out of range 
  @retval 2 - assignment of different number of dofs at slected nodes
  @retval 3 - cannot find nodes with required property and entity type
*/
long input_nod_load(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->nn);
  ivector col(Top->nn);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  ivector setnodes(Top->nn);
  double *force;
  long i, j, k, l;
  long prop, ndof, lcid;
  long tnlc, slcid=0;

  tnlc = Nlc;
  if (Tnslc)
    tnlc = Tnslc;

  Nod_load = new double**[Top->nn];
  memset(Nod_load, 0, sizeof(*Nod_load)*Top->nn);
  for(i=0; i<Top->nn; i++)
  {
    Nod_load[i] = new double*[tnlc];
    memset(Nod_load[i], 0, sizeof(*Nod_load[i])*tnlc);
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of nodal load");

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "nod_load"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "nod_load", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "nod_load", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "nod_load", "propid", &prop);
      xfscanf(in, "%k%ld", "lc_id", &lcid);
      if ((lcid < 1) || (lcid > Nlc))
      {
        sprintf(errmsg, "Line %ld, col %ld, file %s: loadcase index out of range <1,%ld>", in->line, in->col, in->fname, Nlc);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 1;
      }      
      lcid--;
      if (Tnslc)
      {
        if (Mp->tprob == forced_dynamics){
          if (xfscanf(in, "%+k%ld", "slc_id", &slcid) == 2){
            if ((slcid < 1) || (slcid > Nslc[lcid]))
            {
              sprintf(errmsg, "Line %ld, col %ld, file %s: subloadcase index out of range <1,%ld>", in->line, in->col, in->fname, Nslc[lcid]);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              return 1;
            }        
            // do NOT decrement slcid, because the first hidden subloadcase with index 0 is intended for the static loadcase
          }
          else
            slcid = 0; // no slc_id keyword => static loadcase is being specified -> hidden subload case id is 0
        }
        else{
          xfscanf(in, "%k%ld", "slc_id", &slcid);
          if ((slcid < 1) || (slcid > Nslc[lcid])){
            sprintf(errmsg, "Line %ld, col %ld, file %s: subloadcase index out of range <1,%ld>", in->line, in->col, in->fname, Nslc[lcid]);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 1;
          }        
          slcid--;
        }
        lcid = Nslc_cum[lcid]+slcid;
      }
      xfscanf(in, "%k", "load_comp");
      // investigation of ndof for nodes with property prop of entity ent
      ndof = Top->get_ndofn(prop, ent, Nod_ndof, setnodes);
      if (ndof == 0)
      {
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s\n"
                        " have different number of dofs. Nodal load \n" 
                        " (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 2;
      }
      if (ndof == -1)
      {
        if (Check_unused_prop == 0)
        // no nodes were found with the given property and checking of unused properties is off
          continue;
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Nodal load (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 3;
      }
      force = new double[ndof];
      for (k = 0; k < ndof; k++)
        xfscanf(in, "%le", force+k);

      in->kwdmode = bkwdmode;

      for (k=0; k<Top->nn; k++)
      {        
        if (setnodes[k] < 0) // node has not required property
          continue;
        if (Nod_load[k][lcid] == NULL)
        { // node has not assigned load yet
          Nod_load[k][lcid]  = new double[ndof];
          memset(Nod_load[k][lcid], 0, sizeof(*Nod_load[k][lcid])*ndof);
          for(l=0; l<ndof; l++)
            Nod_load[k][lcid][l] += force[l];
          // backup of line and column for eventual log message
          line[k] = aline;
          col[k]  = acol;
        }
        else
        { // node has already assigned load
          fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already had assigned load at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
          for(l=0; l<ndof; l++)
            Nod_load[k][lcid][l] += force[l];
        }
      }
      delete [] force;
    } // end of nkwd loop
  }   // end of section loop
  return 0;
}



/**
  Function reads time dependent nodal load. It scans sections nodvertpr, nodedgpr, nodsurfpr, 
  nodvolpr for keyword "nod_tdload". 
  The record which assigns nodal load looks like follows:
  "nod_tdload" "propid" prop "load_comp" {time_functions_of_load_components}[ndofn]

  The asssigning starts at section for volumes and finishes at
  section for vertices. If a node belongs to several entities and loads are
  prescribed for these entities, the assigned loads are merged and a message 
  is written into the log file.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - subloadcases cannot be required with this type of load
  @retval 2 - index of load case is out of range
  @retval 3 - assignment of different number of dofs at one node
  @retval 4 - cannot find nodes with required property and entity type
*/
long input_nod_tdload(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->nn), col(Top->nn);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  ivector setnodes(Top->nn);
  gfunct *force;
  long i, j, k, l;
  long prop, ndof, lcid;

  Nod_tdload = new gfunct**[Top->nn];
  memset(Nod_tdload, 0, sizeof(*Nod_tdload)*Top->nn);
  for(i=0; i<Top->nn; i++)
  {
    Nod_tdload[i] = new gfunct*[Nlc];
    memset(Nod_tdload[i], 0, sizeof(*Nod_tdload[i])*Nlc);
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of time dependent nodal load");

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "nod_tdload"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "nod_tdload", 1);
    if (nkwd && Tnslc)
    {
      xfscanf(in,"%k", "nod_tdload");
      sprintf(errmsg, "Line %ld, col %ld, file %s: invalid type of load is required in case of subloadcases", in->line, in->col, in->fname);
      print_err(errmsg, __FILE__, __LINE__, __func__);
      return 1;
    }
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "nod_tdload", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "nod_tdload", "propid", &prop);
      xfscanf(in, "%k%ld", "lc_id", &lcid);
      if ((lcid < 1) || (lcid > Nlc))
      {
        sprintf(errmsg, "Line %ld, col %ld, file %s: loadcase index out of range <1,%ld>", in->line, in->col, in->fname, Nlc);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 2;
      }
      lcid--;
      xfscanf(in, "%k", "load_comp");
      // investigation of ndof for nodes with property prop of entity ent
      ndof = Top->get_ndofn(prop, ent, Nod_ndof, setnodes);
      if (ndof == 0)
      {
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s\n"
                        " have different number of dofs. Nodal time dependent load \n" 
                        " (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 3;
      }
      if (ndof == -1)
      {
        if (Check_unused_prop == 0)
        // no nodes were found with the given property and checking of unused properties is off
          continue;
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Nodal time dependent load (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 4;
      }
      force = new gfunct[ndof];
      for (k = 0; k < ndof; k++)
        force[k].read(in);

      in->kwdmode = bkwdmode;

      for (k=0; k<Top->nn; k++)
      {        
        if (setnodes[k] < 0) // node has not required property
          continue;
        if (Nod_tdload[k][lcid] == NULL)
        { // node has not assigned load yet
          Nod_tdload[k][lcid]  = new gfunct[ndof];
          for(l=0; l<ndof; l++)
            Nod_tdload[k][lcid][l].copy(force[l]);
          // backup of line and column for eventual log message
          line[k] = aline;
          col[k]  = acol;
        }
        else
        { // node has already assigned load
          fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already had assigned load at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
          for(l=0; l<ndof; l++)
            Nod_tdload[k][lcid][l].merge(force[l]);
        }
      }
      delete [] force;
    } // end of nkwd loop
  }   // end of section loop
  return 0;
}



/**
  Function reads initial conditions at nodes. It scans sections nodvertpr, nodedgpr, nodsurfpr, 
  nodvolpr for keyword "nod_inicond". 
  The record which assigns nodal initial conditions looks like follows:
  "nod_inicond" "propid" prop "cond" {initial_condition}

  The asssigning starts at section for volumes and finishes at
  section for vertices. If a node belongs to several entities and initial conditions are
  prescribed for these entities, the assigned conditions are merged and a message 
  is written into the log file.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - invalid load case index
  @retval 2 - assignment of different types of initial conditons at one node
  @retval 3 - nodes with required property and entity type cannot be found
*/
long input_nod_initcond(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char emsg[1001];
  ivector line(Top->nn), col(Top->nn);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  ivector setnodes(Top->nn);
  inicd *tic;
  long check;
  long i, j, k, lcid, prop;
  long nini;
  long prop_used;

  Nod_inicd = new inicd**[Top->nn];
  memset(Nod_inicd, 0, sizeof(*Nod_inicd)*Top->nn);
  for(i=0; i<Top->nn; i++)
  {
    Nod_inicd[i] = new inicd*[Nlc];
    memset(Nod_inicd[i], 0, sizeof(*Nod_inicd[i])*Nlc);
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of nodal initial conditions");

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "nod_inicond"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "nod_inicond", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "nod_inicond", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld%k%ld", "nod_inicond", "propid", &prop, "lc_id", &lcid);
      if ((lcid < 1) || (lcid > Nlc))
      {
        sprintf(emsg, "Line %ld, col %ld, file %s: loadcase index out of range <1,%ld>", in->line, in->col, in->fname, Nlc);
        print_err(emsg, __FILE__, __LINE__, __func__);
        return 1;
      }
      lcid--;
      xfscanf(in, "%k", "cond");
      tic = new inicd;
      tic->read(in);

      in->kwdmode = bkwdmode;

      // assigning of initial conditions to the given nodes
      Top->get_propent_nodes(prop, ent, setnodes);
      prop_used = 0;
      for (k=0; k<Top->nn; k++)
      {        
        if (setnodes[k] < 0) // node has not required property
          continue;
        if (Nod_inicd[k][lcid] == NULL)
        { // node has not assigned load yet
          Nod_inicd[k][lcid] = new inicd;
          Nod_inicd[k][lcid]->copy(*tic); 
          // backup of line and column for eventual log message
          line[k] = aline;
          col[k]  = acol;
        }
        else
        { // node has already assigned initial condition
          check =  Nod_inicd[k][lcid]->merge(*tic);
          if (check)
          { // error occured while merging
            sprintf(emsg, "Line %ld, col %ld, file %s: Node %ld has assigned two different types of "
                            "initial condition at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
            print_err(emsg, __FILE__, __LINE__, __func__);
            return 2;
          }
          else
            fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already had assigned initial condition at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
        }
        prop_used = 1;
      }
      delete tic;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no nodes were found with the given property and checking of unused properties is on
        sprintf(emsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                      " Nodal intial conditions (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(emsg, __FILE__, __LINE__, __func__);
        return 3;
      }
    }  // end of nkwd loop
  }    // end of section loop

  // check for completeness of intial condition 
  // in particular load cases, number of initial conditions can be either 0 or Top->nn
  if (Mp->tprob != growing_mech_structure)
  {
    for (i=0; i<Nlc; i++)
    {
      nini = 0;
      for (j=0; j<Top->nn; j++)
      {
        if (Nod_inicd[j][i])
          nini++;
      }
      if ((nini > 0) && (nini < Top->nn))
      {
        for (j=0; j<Top->nn; j++)
        {
          if (Nod_inicd[j][i] == NULL)
            break;
        }
        sprintf(emsg, "number of initial conditions (%ld) for loadcase %ld is less then number of nodes (%ld).\n"
                "The first node without initial condition is %ld", nini, i+1, Top->nn, j+1);
        print_err(emsg, __FILE__, __LINE__, __func__);
        return 3;
      }
    }
  }
  return 0;
}



/**
  Function reads prescribed initial displacements by rotation of selected nodes. It scans sections nodvertpr, nodedgpr, nodsurfpr, 
  nodvolpr for keyword "nod_inipd_rot". 
  The record which assigns nodal initial prescribed displacements looks like follows:
  "nod_inipd_rot" "propid" prop "axis_rec" {rotation_record}.
  This function is applied only for the growing mechanical problems.

  The asssigning starts at section for volumes and finishes at
  section for vertices. If a node belongs to several entities and initial prescribed displacements are
  prescribed for these entities, the both assigned conditions are applied and a message 
  is written into the log file.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - invalid load case index
  @retval 2 - assignment of different types of initial conditons at one node
  @retval 3 - nodes with required property and entity type cannot be found
*/
long input_nod_rotinidispl(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  ivector line(Top->nn), col(Top->nn);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  ivector setnodes(Top->nn);
  long i, j, k, prop;
  long nl, nr;
  long aaxis;
  long prop_used;

  Nod_rot_ipd_num  = 0L;
  Nod_rot_ipd_axis = NULL;

  // command "nod_inipd_rot" is intended for the gradual contruction problems only
  if (Mp->tprob != growing_mech_structure)
    return 0;

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of rotation axis for precribed initial displacements");

  nkwd = 0L;
  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "nod_inicond"
    in->kwdmode = sect_mode_full;
    nkwd += getkwd_sect(in, NULL, aptr, "nod_inipd_rot", 1);
  }
  Nod_rot_ipd_num = nkwd;
  Nod_rot_ipd_axis = new axisrotrec[nkwd];

  aaxis = 0L;
  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "nod_inicond"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "nod_inipd_rot", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "nod_inipd_rot", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "nod_inipd_rot", "propid", &prop);
      xfscanf(in, "%k", "axis_rec");
      Nod_rot_ipd_axis[aaxis].read_prep(in);
      if ((Nod_rot_ipd_axis[aaxis].stime < Mp->timecon.starttime()) || 
          (Nod_rot_ipd_axis[aaxis].stime > Mp->timecon.endtime()) ||
          (Nod_rot_ipd_axis[aaxis].etime > Mp->timecon.endtime()) ||
          (Nod_rot_ipd_axis[aaxis].etime < Mp->timecon.starttime()))
      {
        print_err("Time period <%le,%le> of %ld. axis is out of range <%le, %le>", 
                  __FILE__, __LINE__, __func__, Nod_rot_ipd_axis[aaxis].stime, Nod_rot_ipd_axis[aaxis].etime, 
                  aaxis+1, Mp->timecon.starttime(), Mp->timecon.endtime());
      }
      in->kwdmode = bkwdmode;

      // assigning of initial prescribed displacements to the given nodes
      Top->get_propent_nodes(prop, ent, setnodes);
      prop_used = 0;
      for (k=0; k<Top->nn; k++)
      {        
        if (setnodes[k] < 0) // node has not required property
          setnodes[k] = 0L;
        else
        {
          prop_used = 1L;
          setnodes[k] = 1L;
        }
      }
      if (prop_used)
      {
        nl = Nod_rot_ipd_axis[aaxis].seln.give_num_lst_items(Top->nn, setnodes.a);   // number of  selected objects 
        nr = Nod_rot_ipd_axis[aaxis].seln.give_num_range_items(Top->nn, setnodes.a); // number of ranges in the selection of objects
        if (2*nr < nl)  // range needs two numbers for the specification
          Nod_rot_ipd_axis[aaxis].seln.conv2range(nr, Top->nn, setnodes);
        else
          Nod_rot_ipd_axis[aaxis].seln.conv2lst(nl, Top->nn, setnodes);
      }
      else
      {
        // no nodes were found with the given property and checking of unused properties is on
        print_err("Nodes with property %ld belonging to entity %s have not been found.\n"
                  " Nodal prescribed intial displacements (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                  __FILE__, __LINE__, __func__, prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        return 3;
      }
      aaxis++;
    }
  }
  return 0;
}



/**
  Function reads nodal temperatures. It scans sections nodvertpr, nodedgpr, nodsurfpr, 
  nodvolpr for keyword "nod_temper". 
  The record which assigns nodal temperatures looks like follows:
  "nod_temper" "propid" prop "lc_id" load_case_index "temperature" temperature_value 

  The asssigning starts at section for volumes and finishes at
  section for vertices. If a node belongs to several entities and temperatures are
  prescribed for these entities, the assigned temperatures are OVERWRITTEN and a message 
  is written into the log file.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - wrong loadcase/subloadcase id
  @retval 2 - wrong temperature load type is required (it should be 1 or 2)
  @retval 3 - nodes with required property and entity type cannot be found
*/
long input_nod_temper(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  char emsg[1001];
  long nkwd;
  ivector line(Top->nn), col(Top->nn);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  ivector setnodes(Top->nn);
  tempload *tl;
  long prop;
  long check, tnlc;
  long i, j, k, id;
  long prop_used;

  tnlc = Nlc;
  if (Tnslc)
    tnlc = Tnslc;

  Ntempl = 0;
  Nod_temper = new tempload**[Top->nn];
  memset(Nod_temper, 0, sizeof(*Nod_temper)*Top->nn);
  for(i=0; i<Top->nn; i++)
  {
    Nod_temper[i] = new tempload*[tnlc];
    memset(Nod_temper[i], 0, sizeof(*Nod_temper[i])*tnlc);
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of temperature load at nodes");

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "nod_temper"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "nod_temper", 1);
    in->kwdmode = bkwdmode;
    if (nkwd) // set indicator of temperature load in the given problem
      Ntempl = 1;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "nod_temper", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "nod_temper", "propid", &prop);
      tl = new tempload;
      check = tl->read(in, Nlc, Nslc);
      if (check)
      {
        sprintf(emsg, "Line %ld, col %ld, file %s: wrong loadcase/subloadcase index", in->line, in->col, in->fname);
        print_err(emsg, __FILE__, __LINE__, __func__);
        return 1;
      }

      id = tl->nlc-1;
      if (Tnslc){
        if (Mp->tprob == forced_dynamics)
          id = Nslc_cum[tl->nlc-1]+tl->nslc;          
        else
          id = Nslc_cum[tl->nlc-1]+tl->nslc-1;
      }
    
      if ((Tlt[id] < 1) || (Tlt[id] > 2))
      {
        sprintf(emsg, "Line %ld, col %ld, file %s: wrong temperature load type (%ld) is required\n for load case %ld", in->line, in->col, in->fname, Tlt[id], tl->nlc);
        if (Nslc[tl->nlc-1])
          sprintf(emsg+strlen(emsg), ", subloadcase %ld (it should be 1 or 2)", tl->nslc);
        print_err(emsg, __FILE__, __LINE__, __func__);
        return 2;
      }

      in->kwdmode = bkwdmode;

      // assigning of temperatures to the given nodes
      Top->get_propent_nodes(prop, ent, setnodes);
      prop_used = 0;
      for (k=0; k<Top->nn; k++)
      {        
        if (setnodes[k] < 0) // node has not required property
          continue;
        if (Nod_temper[k][id] == NULL)
        { // node has not assigned temperature yet
          Nod_temper[k][id] = new tempload;
          Nod_temper[k][id]->copy(*tl); 
          prop_used = 1;
          // backup of line and column for eventual log message
          line[k] = aline;
          col[k]  = acol;
        }
        else
        {  // node has already assigned temperature
          Nod_temper[k][id]->copy(*tl); 
          fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already had assigned temperature at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
        }
      }
      delete tl;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no nodes were found with the given property and checking of unused properties is on
        sprintf(emsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Nodal temperatures (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(emsg, __FILE__, __LINE__, __func__);
        return 3;
      }
    }  // end of nkwd loop
  }    // end of section loop
  return 0;
}



/**
  Function reads periodic boundary conditions at nodes. It scans sections nodvertpr, 
  nodedgpr, nodsurfpr,  nodvolpr for keyword "nod_periodbc". 
  The record which assigns nodal periodic boundary conditions looks like follows:
  "nod_periodbc" "propid_master prop "propid_slave" propsl "max_nit" ni "comp_err" cerr

  The asssigning starts at section for volumes and finishes at
  section for vertices. If a node belongs to several entities and periodic boundary conditions
  are  prescribed for these entities, the assigned boundary conditions are compared if they are 
  identical. If the master element/node id is different, the periodic condition is marked
  for removal, message is written into the log file and regular boundary conditions (i.e. support) 
  are being expected in those nodes.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - wrong loadcase/subloadcase id
  @retval 2 - wrong temperature load type is required (it should be 1 or 2)
  @retval 3 - nodes with required property and entity type cannot be found
*/
long input_nod_periodic_bc(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  ivector line(Top->nn), col(Top->nn);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  ivector setmnodes(Top->nn), setsnodes(Top->nn), nodf(Top->nn), umnodes(Top->nn), melem(Top->ne);
  long prop, propsl;
  long i, j, k, l;
  vector bv1(ASTCKVEC(3)), bv2(ASTCKVEC(3)), pv1(ASTCKVEC(3)), pv2(ASTCKVEC(3)), nv(ASTCKVEC(3));
  long nne, eldim;
  vector cg(ASTCKVEC(3));
  vector x, y, z, p(ASTCKVEC(3));
  double norm2, d2, r2, xi, eta, zeta, cerr, aerr, max_err;
  long nsn, nmn, snn, mnn, nme, nfsn, max_errn, mei;
  long ni, nrerr, ncall;
  long cmdid = 0;
  //long nid;
  //double dmin;

  Nperbc = 0;
  Nod_periodbc = new msmap[Top->nn];
  for (i=0; i<Top->nn; i++)
    Nod_periodbc[i].slmas = 0;

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of periodic boundary conditions at nodes");
  fflush(stdout);

  ncall = 0;
  for (i=0, ent=gvertex; i<nsect; i++, ent++) {
  //  for (i=nsect-1, ent=gsurface; i>=0; i--, ent--) {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "nod_periodicbc"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "nod_periodbc", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++) {
      nullv(melem);
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "nod_periodbc", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k %k%ld %k%ld %k%ld %k%le",
              "nod_periodbc", "propid_master", &prop, "propid_slave", &propsl,
              "max_nit", &ni, "comp_err", &cerr);
      cmdid++;
      
      
      in->kwdmode = bkwdmode;

      //
      // assigning of periodic boundary conditions to the given nodes
      //
      // search all master nodes of the given property
      nmn = Top->get_propent_nodes_compact(prop, ent, setmnodes);
      nsn = Top->get_propent_nodes_compact(propsl, ent, setsnodes);
      nme = Top->get_propent_elems_compact(prop, ent, melem);
      if (Check_unused_prop && ((nmn == 0) || (nsn == 0))) {
        // no nodes were found with the given property and checking of unused properties is on
        print_err("Nodes with property %ld belonging to entity %s have not been found.\n"
                  " Nodal periodic boundary conditions (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                  __FILE__, __LINE__, __func__, prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        return 1;
      }
      if ((ent == gcurve) && (nmn < 2)) {
        print_err("the number of master nodes (%ld) with curve property %ld is less than 2.\n"
                  " Nodal periodic boundary conditions (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                  __FILE__, __LINE__, __func__, prop, aline, acol, in->fname);
        return 1;
      }
      if ((ent == gsurface) && (nmn < 3)) {
        print_err("the number of master nodes (%ld) with surface property %ld is less than 3.\n"
                  " Nodal periodic boundary conditions (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                  __FILE__, __LINE__, __func__, prop, aline, acol, in->fname);
        return 1;
      }
      
      if (ent == gcurve) {
        bv1(0) = Top->nodes[setmnodes(1)].x - Top->nodes[setmnodes(0)].x;
        bv1(1) = Top->nodes[setmnodes(1)].y - Top->nodes[setmnodes(0)].y;
        bv1(2) = Top->nodes[setmnodes(1)].z - Top->nodes[setmnodes(0)].z;        
        normalize(bv1);
        for(k=0; k<nsn; k++){
          snn = setsnodes(k);
          pv1(0) = Top->nodes[snn].x - Top->nodes[setmnodes(0)].x;
          pv1(1) = Top->nodes[snn].y - Top->nodes[setmnodes(0)].y;
          pv1(2) = Top->nodes[snn].z - Top->nodes[setmnodes(0)].z;
          scprd(pv1, bv1, norm2);
          cmulv(norm2, bv1, pv2);
          subv(pv1, pv2, pv2);
          if (Nod_periodbc[snn].slmas == 0)// || // no previous periodic bc has been assigned to the node with snn id
            //              (Nod_periodbc[snn].slmas && (Nod_periodbc[snn].ent != ent)))  // periodic bc have been already assigned to the node with snn id but on the different type of entity
          {
            Nod_periodbc[snn].x = Top->nodes[snn].x - pv2(0);
            Nod_periodbc[snn].y = Top->nodes[snn].y - pv2(1);
            Nod_periodbc[snn].z = Top->nodes[snn].z - pv2(2);
            Nod_periodbc[snn].slmas = 1;
            Nod_periodbc[snn].ent = ent;
            line[snn] = aline;
            col[snn] = acol;
          }
          else {
            fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already had a periodic boundary condition assigned at line %ld, col %ld, file: %s."
                         " All periodic boundary conditions in the above node will be ignored.\n",
                    aline, acol, in->fname, snn+1, line[snn], col[snn], in->fname);
            Nod_periodbc[snn].slmas = 2;
          }
        }
      }
      if (ent == gsurface) {
        k = 1;
        //
        // compute normal vector to the plane formed by the selected master nodes
        //
        // compute the first base vector in the plane
        do{
          bv1(0) = Top->nodes[setmnodes(k)].x - Top->nodes[setmnodes(0)].x;
          bv1(1) = Top->nodes[setmnodes(k)].y - Top->nodes[setmnodes(0)].y;
          bv1(2) = Top->nodes[setmnodes(k)].z - Top->nodes[setmnodes(0)].z;
          norm2 = normv(bv1);
          if (norm2 > cerr) {
            cmulv(1.0/norm2, bv1);
            break;
          }
          k++;
        } while(k < nsn);
        if (k == nsn) {
          print_err("master nodes of surface with propery %ld are concentrated in one point,\n cannot compute normal vector of the surface\n"
                    " Nodal periodic boundary conditions (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                    __FILE__, __LINE__, __func__, prop, aline, acol, in->fname);
          return 1;
        }
        k = 2;
        // compute the second base vector in the plane
        do {
          bv2(0) = Top->nodes[setmnodes(k)].x - Top->nodes[setmnodes(0)].x;
          bv2(1) = Top->nodes[setmnodes(k)].y - Top->nodes[setmnodes(0)].y;
          bv2(2) = Top->nodes[setmnodes(k)].z - Top->nodes[setmnodes(0)].z;
          norm2 = normv(bv2);
          if (norm2 <= cerr) {
            k++;
            continue;
          }
          cmulv(1.0/norm2, bv2);
          scprd(bv1, bv2, norm2);
          norm2 = sqrt(norm2);
          norm2 -= 1.0;
          if (fabs(norm2) > cerr)
            break;
          else
            k++;
        } while(k < nsn);
        if (k == nsn){
          print_err("master nodes of surface with propery %ld form stright line, cannot compute normal vector of the surface\n"
                    " Nodal periodic boundary conditions (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                    __FILE__, __LINE__, __func__, prop, aline, acol, in->fname);
          return 1;
        }
        // compute surface normal vector with the help of cross product
        crprd(bv1, bv2, nv);
        norm2 = normv(nv);
        cmulv(1.0/norm2, nv);
        
        for(k=0; k<nsn; k++) {
          snn = setsnodes(k);
          if (Nod_periodbc[snn].slmas == 0)// || // no previous periodic bc has been assigned to the node with snn id
            //              (Nod_periodbc[snn].slmas && (Nod_periodbc[snn].ent != ent)))  // periodic bc have been already assigned to the node with snn id but on the different type of entity
          {   
            pv1(0) = Top->nodes[setmnodes(0)].x - Top->nodes[snn].x;
            pv1(1) = Top->nodes[setmnodes(0)].y - Top->nodes[snn].y;
            pv1(2) = Top->nodes[setmnodes(0)].z - Top->nodes[snn].z;
            scprd(nv, pv1, norm2);
            cmulv(norm2, nv, pv2);
            Nod_periodbc[snn].x = Top->nodes[setsnodes(k)].x + pv2(0);
            Nod_periodbc[snn].y = Top->nodes[setsnodes(k)].y + pv2(1);
            Nod_periodbc[snn].z = Top->nodes[setsnodes(k)].z + pv2(2);
            Nod_periodbc[snn].slmas = 1;
            Nod_periodbc[snn].ent = ent;            
            line[snn] = aline;
            col[snn] = acol;
          }
          else {
            fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already had a periodic boundary condition assigned for entity type %s at line %ld, col %ld, file: %s."
                         " All other assignments of periodic BC will be ignored.\n",
                    //                    " All assignments of periodic BC will be ignored for this entity type.\n",
                    aline, acol, in->fname, snn+1, gentity_kwdset.get_str(ent), line[snn], col[snn], in->fname);
              Nod_periodbc[snn].slmas = 2;
          }
        }
      }
      for(k=0; k<setmnodes.n; k++){
        mnn = setmnodes(k);
        if (Nod_periodbc[mnn].slmas == 0){
          Nod_periodbc[mnn].slmas = -1;
          Nod_periodbc[mnn].ent = ent;            
          line(mnn) = aline;
          col(mnn) = acol;
        }
        if (Nod_periodbc[mnn].slmas > 0){
          fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already had a periodic boundary condition assigned for entity type %s at line %ld, col %ld, file: %s."
                       " All other assignments of periodic BC will be ignored.\n",
                  //                       " All assignments of periodic BC will be ignored for this entity type.\n",
                  aline, acol, in->fname, mnn+1, gentity_kwdset.get_str(ent), line[mnn], col[mnn], in->fname);
          Nod_periodbc[mnn].slmas = 2;
        }        
      }
      nfsn = 0;
      for(k=0; k<nme; k++) {
        mei = melem(k);
        Top->centroid(mei, cg);
        r2 = Top->max_sqrdist_nod_pt(mei, cg);
        r2 = sqr(sqrt(r2)+cerr); // extend the square of radius by the tolerance
        for(l=0; l<nsn; l++){
          snn = setsnodes[l];          
          if (nodf(l) || (Nod_periodbc[snn].slmas == 2)) // periodic BC has been already assigned in previous step or even the node has already assigned multiple periodic BC
            continue;
          d2  = sqr(Nod_periodbc[snn].x-cg(0)) + sqr(Nod_periodbc[snn].y-cg(1)) + sqr(Nod_periodbc[snn].z-cg(2));
          if (d2 <= r2){
            /*
            nid = Top->give_closest_node_coord(mei, Nod_periodbc[snn].x, Nod_periodbc[snn].y, Nod_periodbc[snn].z, dmin);
            if (Nod_periodbc[snn].cerr > dmin){
              if ((Nod_periodbc[snn].nid < 0) && (Nod_periodbc[snn].eid < 0)){
                Nperbc++;
                nfsn++;
              }
              Nod_periodbc[snn].nid = nid;
              Nod_periodbc[snn].cerr = dmin;
              xi = 0.0;
              eta = 0.0;
              zeta = 0.0;
              if (dmin < cerr)
                nodf(l) = 1;
              if (Nod_periodbc[snn].eid >= 0)
                Nod_periodbc[snn].eid = -1;
            }*/
            eldim = Top->give_dimension(mei);
            nne = Top->elements[mei].nne;
            reallocv(RSTCKVEC(nne, x));
            reallocv(RSTCKVEC(nne, y));
            reallocv(RSTCKVEC(nne, z));
            switch (eldim){
              case 1:
                Top->give_node_coord3d(x, y, z, mei);
                xi = eta = zeta = 0.0;
                nrerr = point_natcoord_1d_3d(Nod_periodbc[snn].x, Nod_periodbc[snn].y, Nod_periodbc[snn].z, x, y, z, ni, cerr, xi, aerr);
                corr_nat_coord_bounds(Top->elements[mei].type, xi, eta, zeta);
                bf_1d_3d(p, x, y, z, xi);
                aerr = length(p, Nod_periodbc[snn].x, Nod_periodbc[snn].y, Nod_periodbc[snn].z);
                break;
              case 2:
                Top->give_node_coord2d(x, y, mei);
                xi = eta  = zeta = 0.0;
                nrerr = point_natcoord_2d(Nod_periodbc[snn].x, Nod_periodbc[snn].y, x, y, ni, cerr, xi, eta, aerr);
                corr_nat_coord_bounds(Top->elements[mei].type, xi, eta, zeta);
                bf_2d(p, x, y, xi, eta);
                aerr = length(p, Nod_periodbc[snn].x, Nod_periodbc[snn].y, Nod_periodbc[snn].z);
                break;
              case 3:
                Top->give_node_coord3d(x, y, z, mei);
                xi = eta  = zeta = 0.0;
                nrerr = point_natcoord_3d(Nod_periodbc[snn].x, Nod_periodbc[snn].y, Nod_periodbc[snn].z, x, y, z, ni, cerr, xi, eta, zeta, aerr);                
                corr_nat_coord_bounds(Top->elements[mei].type, xi, eta, zeta);
                bf_3d(p, x, y, z, xi, eta, zeta);
                aerr = length(p, Nod_periodbc[snn].x, Nod_periodbc[snn].y, Nod_periodbc[snn].z);
                ncall++;
                break;
              default:
                print_err("unknown dimension (dim=%ld) of the element is required", __FILE__, __LINE__, __func__, eldim);
                return 2;
            }
            if (nrerr > ni){
              print_err("solution of natural coordinates of slave node %ld did not converge in %ld steps on element %ld ", __FILE__, __LINE__, __func__, snn+1, nrerr, mei+1);
              return 2;
            }
            if (Nod_periodbc[snn].cerr > aerr){
              if ((Nod_periodbc[snn].nid < 0) && (Nod_periodbc[snn].eid < 0)){
                Nperbc++;
                nfsn++;
              }
              Nod_periodbc[snn].xi = xi;
              Nod_periodbc[snn].eta = eta;
              Nod_periodbc[snn].zeta = zeta;
              Nod_periodbc[snn].cerr = aerr;
              Nod_periodbc[snn].eid = mei;
              Nod_periodbc[snn].nid = -1;
              if (aerr < cerr)
                nodf(l) = 1;
            }
          } // node is member of neighbourhood of the mei-th element center 
        } // end of loop over slave nodes
      }   // end of loop over all elements
      if (nsn != nfsn){
        for (l=0; l<nsn; l++){
          snn = setsnodes[l];
          if ((Nod_periodbc[snn].slmas == 1) && (Nod_periodbc[snn].eid < 0) && (Nod_periodbc[snn].nid < 0)){
            print_err("Line %ld, col %ld, file %s: cannot find master node (prop=%ld) for the slave node %ld (prop=%ld) on entity %s.\n"
                      " Check planarity, parallelism and shape correspondence of selected entities of master and slave nodes\n"
                      " in the periodic boundary condition assignment.",
                      __FILE__, __LINE__, __func__,  aline, acol, in->fname, prop, snn+1, propsl, gentity_kwdset.get_str(ent));
            abort();
          }
        }
      }
      periodbc_log(cmdid, in, aline, acol, setsnodes, setmnodes, cerr);
      nullv(nodf);
    }  // end of nkwd loop
  }    // end of section loop

  max_err = 0.0;
  aerr = 0.0;
  nrerr = 0;
  for (k=0; k<Top->nn; k++) {
    if ((Nod_periodbc[k].slmas > 0) && (Nod_periodbc[k].cerr > cerr)) {
      nrerr++;
      aerr += Nod_periodbc[k].cerr;
      if (Nod_periodbc[k].cerr > max_err){
        max_errn = k+1;  
        max_err = Nod_periodbc[k].cerr;
      }
      if (Nod_periodbc[k].eid >= 0)
        fprintf(Log, "Line %ld, col %ld, file %s: master node of node %ld cannot be found in the given tolerance (err=%le > cerr=%le on eid=%ld).\n",
                line[k], col[k], in->fname, k+1, Nod_periodbc[k].cerr, cerr, Nod_periodbc[k].eid+1);
      else {
        if (Nod_periodbc[k].nid >= 0)
          fprintf(Log, "Line %ld, col %ld, file %s: master node of node %ld cannot be found in the given tolerance (err=%le > cerr=%le on nid=%ld).\n",
                  line[k], col[k], in->fname, k+1, Nod_periodbc[k].cerr, cerr, Nod_periodbc[k].nid+1);
        else{
          if (Nod_periodbc[k].cerr != DBL_MAX)
            fprintf(Log, "Line %ld, col %ld, file %s: master node of node %ld cannot be found in the given tolerance (err=%le > cerr=%le master element/node was not determined).\n",
                    line[k], col[k], in->fname, k+1, Nod_periodbc[k].cerr, cerr);
        }
      }
    }
  }
  if (nrerr){
    fprintf(stderr, "\n\nWarning:\n The total number of nodes with periodic boundary conditions: %ld\n", Nperbc);
    fprintf(stderr, " The number of nodes whose coordinates were computed beyond the given tolerance: %ld\n", nrerr);
    fprintf(stderr, "  - required tolerance: %le\n", cerr);
    fprintf(stderr, "  - maximum error: %le at node %ld\n", max_err, max_errn);
    fprintf(stderr, "  - average error: %le\n", aerr/nrerr);
    fprintf(stderr, "Total number of calls of nat. coord. computation: %ld\n", ncall);
  }  
  return 0;
}



/**
  The function searches and generates hanging nodes from a set with a given property.
  Searched elements are given by an element set with given properties. 

  @param in[in]  - pointer to opened input file with property description
  @nodsects[in]  - array with descriptors of sections, which will be searched
  @nsect[in]     - number of searched sections i.e number of elements in array nodsects

  The record which assigns nodal periodic boundary conditions looks like follows:
  "nod_genhn" "propid" prop "num_prop_master_elem" npme "propid_set" {prop_me_i} x npme
  "max_nit" ni "comp_err" cerr

  Returns :
  @retval 0 - on succes  
  @retval 1 - an error
*/
long input_nod_genhn(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  ivector line(Top->nn), col(Top->nn);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  ivector setnodes(Top->nn);
  long prop, npme;
  long i, j, k, l;
  ivector melpr;
  ivector melset(Top->ne);
  long nhn=0, nme=0;

  long mei, eldim, nne, hnf, thnf=0, snn, ni, nrerr;
  vector cg(ASTCKVEC(3)), p(ASTCKVEC(3)), x, y, z;
  double xi, eta, zeta, aerr, cerr, r2, d2;
  ivector hnodf(Top->nn);
  hngen *nod_hngen = new hngen[Top->nn];



  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of setup of hanging node generation");

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;

    // detect number of occurences of keyword "nod_temper"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "nod_genhn", 1);
    in->kwdmode = bkwdmode;
    
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "nod_genhn", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read property id of hanging node set
      xfscanf(in, "%k%k%ld", "nod_genhn", "propid", &prop);
      // read set of master elements
      xfscanf(in, "%k%ld%k", "num_prop_master_elem", &npme, "propid_set");
      reallocv(npme, melpr);
      readv(in, melpr);
      xfscanf(in, "%k%ld%k%le", "max_nit", &ni, "comp_err", &cerr);
      in->kwdmode = bkwdmode;

      // search for set hanging nodes
      nullv(setnodes);
      nhn = Top->get_propent_nodes_compact(prop, ent, setnodes);
      if (Check_unused_prop && (nhn <= 0)){
        print_err("Nodes with propeprty %ld belonging to entity %s have not been found.\n"
                  "Command 'nod_genhn' (line=%ld, col=%ld, file=%s) cannot be performed correctly",
                  __FILE__, __LINE__, __func__, prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        return 1;
      }

      // unused properties were not checked and no hanging nodes were found
      // for the given prop and ent => proceed to the next keyword occurence
      if (nhn == 0)  continue;

      fillv(-1, melset);
      nme = 0;
      // search for set of master elements
      for(k=0; k<npme; k++){
        nme += Top->get_propent_elems(melpr(k), ent, melset, 1);
      }      
      if (nme <= 0){
        print_err("no master element were found for the given property id set.\n"
                  "Command 'nod_genhn' (line=%ld, col=%ld, file=%s) cannot be performed correctly",
                  __FILE__, __LINE__, __func__, prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        return 1;
      }

      hnf = 0; // number of found hanging nodes, it must be equal to nhn after the hanging node generation
      for(k=0; k<Top->ne; k++) {
        mei = melset(k);
        if (mei < 0)   continue;
        Top->centroid(mei, cg);
        r2 = Top->max_sqrdist_nod_pt(mei, cg);
        r2 = sqr(sqrt(r2)+cerr); // extend the square of radius by the tolerance
        for(l=0; l<nhn; l++){
          snn = setnodes[l];          
          if (hnodf(snn)) // hanging node has been already assigned in previous step
            continue;
          d2  = sqr(Top->nodes[snn].x-cg(0)) + sqr(Top->nodes[snn].y-cg(1)) + sqr(Top->nodes[snn].z-cg(2));
          if (d2 <= r2){
            eldim = Top->give_dimension(mei);
            nne = Top->elements[mei].nne;
            reallocv(RSTCKVEC(nne, x));
            reallocv(RSTCKVEC(nne, y));
            reallocv(RSTCKVEC(nne, z));
            switch (eldim){
              case 1:
                Top->give_node_coord3d(x, y, z, mei);
                xi = eta = zeta = 0.0;
                nrerr = point_natcoord_1d_3d(Top->nodes[snn].x, Top->nodes[snn].y, Top->nodes[snn].z, x, y, z, ni, cerr, xi, aerr);
                corr_nat_coord_bounds(Top->elements[mei].type, xi, eta, zeta);
                bf_1d_3d(p, x, y, z, xi);
                aerr = length(p, Top->nodes[snn].x, Top->nodes[snn].y, Top->nodes[snn].z);
                break;
              case 2:
                Top->give_node_coord2d(x, y, mei);
                xi = eta  = zeta = 0.0;
                nrerr = point_natcoord_2d(Top->nodes[snn].x, Top->nodes[snn].y, x, y, ni, cerr, xi, eta, aerr);
                corr_nat_coord_bounds(Top->elements[mei].type, xi, eta, zeta);
                bf_2d(p, x, y, xi, eta);
                aerr = length(p, Top->nodes[snn].x, Top->nodes[snn].y, Top->nodes[snn].z);
                break;
              case 3:
                Top->give_node_coord3d(x, y, z, mei);
                xi = eta  = zeta = 0.0;
                nrerr = point_natcoord_3d(Top->nodes[snn].x, Top->nodes[snn].y, Top->nodes[snn].z, x, y, z, ni, cerr, xi, eta, zeta, aerr);                
                corr_nat_coord_bounds(Top->elements[mei].type, xi, eta, zeta);
                bf_3d(p, x, y, z, xi, eta, zeta);
                aerr = length(p, Top->nodes[snn].x, Top->nodes[snn].y, Top->nodes[snn].z);
                break;
              default:
                print_err("unknown dimension (dim=%ld) of the element is required", __FILE__, __LINE__, __func__, eldim);
                return 2;
            }
            if (nrerr > ni){
              print_err("solution of natural coordinates of slave node %ld did not converge in %ld steps on element %ld ", __FILE__, __LINE__, __func__, snn+1, nrerr, mei+1);
              return 2;
            }
            if (nod_hngen[snn].cerr > aerr){
              nod_hngen[snn].x = Top->nodes[snn].x;
              nod_hngen[snn].y = Top->nodes[snn].y;
              nod_hngen[snn].z = Top->nodes[snn].z;
              nod_hngen[snn].xi   = xi;
              nod_hngen[snn].eta  = eta;
              nod_hngen[snn].zeta = zeta;
              nod_hngen[snn].cerr = aerr;
              nod_hngen[snn].eid  = mei;
              if (aerr < cerr){
                hnodf(snn) = 1;
                hnf++;
                thnf++;
              }
            }
          } // node is member of neighbourhood of the mei-th element center 
        } // end of loop over slave nodes
      }   // end of loop over all elements
      if (nhn != hnf){
        print_err("some of hanging nodes with property %ld on enity %s could not be found on the given element set.\n"
                  " See log file for more detailed desription.\n"
                  "Command 'nod_genhn' (line=%ld, col=%ld, file=%s) cannot be performed correctly.",
                  __FILE__, __LINE__, __func__, prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        fprintf(Log, "Command 'nod_genhn' (line=%ld, col=%ld, file=%s) cannot be performed correctly.\n", aline, acol, in->fname);
        for (k=0; k<setnodes.n; k++){
          if ((nod_hngen[k].eid >= 0) && (nod_hngen[k].cerr > cerr)){
            fprintf(Log, "|-> Hanging node %ld (x=%le, y=%le, z=%le) was found with tolerance"
                         "%le > %le on element %ld.\n",
                    setnodes(k)+1,
                    nod_hngen[k].x, nod_hngen[k].y, nod_hngen[k].z,
                    nod_hngen[k].cerr, cerr, nod_hngen[k].eid+1);
          }
        }
        return 1;
      }
    }  // end of nkwd loop
  }    // end of section loop

  // process, reduce and set array of master nodes
  get_masternodes(Top->nn, nod_hngen);
  
  // generate array of hanging nodes
  long merge = 1;  
  if ((thnf) && (Numhn == 0)){ // allocate array of hanging node records
    Numhn = thnf;
    Nod_hang = new hangnode*[Top->nn];
    memset (Nod_hang, 0, sizeof(*Nod_hang)*Top->nn);
    merge = 0;
  }
  for (i=0; i<Top->nn; i++){
    if (hnodf(i)){
      if (Nod_hang[i] == NULL){
        Nod_hang[i] = new hangnode;
        Nod_hang[i]->nmn = nod_hngen[i].mnodes.n;
        Nod_hang[i]->mnodes = new long[nod_hngen[i].mnodes.n];
        ivector mnod;
        makerefv(mnod, Nod_hang[i]->mnodes, Nod_hang[i]->nmn);
        copyv(nod_hngen[i].mnodes, mnod);
        Nod_hang[i]->natcoord = new double[3];
        Nod_hang[i]->natcoord[0] = nod_hngen[i].xi;
        Nod_hang[i]->natcoord[1] = nod_hngen[i].eta;
        Nod_hang[i]->natcoord[2] = nod_hngen[i].zeta;
        Nod_hang[i]->maset       = nod_hngen[i].et;
        if (merge)  Numhn++;
      }
      else{
        long cmp = Nod_hang[i]->compare(nod_hngen[i]);
        if (cmp == 0){
          print_err("hanging node %ld has been alread defined in the hanging node input file and\n"
                    "  the generated one has different definition.\n",
                    __FILE__, __LINE__, __func__, i+1);
          return 2;
        }
      }
    }
  }
  delete [] nod_hngen;
  
  return 0;
}




/**
  The function assigns element type to the elements. It scans 
  sections eledgpr, elsurfpr, elvolpr for keyword "el_type".
  The record which assigns type to the element looks like follows:
  "el_type" "propid" prop element_type
  where prop is positive integer number.
  
  The asssigning starts at section for volumes and finishes at
  section for edges. If an element  belongs to several entities and types are
  prescribed for these entities, error of assignment is signalized.
  Finally, the function checks whether all elements have assigned type.
 
  @param in  - pointer to opened input file with property description
  @elemsects - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - different types have been assigned to one element
  @retval 2 - different strain/stress states have been assigned to one element 
  @retval 3 - no type has been assigned to elements
  @retval 4 - property numbers of the required entity type have not been read on elements
  @retval 5 - property number of the required entity type have not been found on elements
*/
long input_elem_type(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->ne), col(Top->ne);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k;
  elemtype tmp_type;
  strastrestate tmp_ssst;
  long prop;
  long sres;
  long prop_used;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of element types");
  El_type = new elemtype[Top->ne];
  memset(El_type, 0, sizeof(*El_type)*Top->ne);
  El_ssst = new strastrestate[Top->ne];
  memset(El_ssst, 0, sizeof(*El_ssst)*Top->ne);

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "el_type"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "el_type", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "el_type", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld%m", "el_type", "propid", &prop, &elemtype_kwdset, &tmp_type);
      if ((tmp_type >= planeelementlt) && (tmp_type <= planequadinterface))
        xfscanf(in, "%k%m", "strastrestate", &strastrestate_kwdset, &tmp_ssst);
      
      in->kwdmode = bkwdmode;

      prop_used = 0;
      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, NULL, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { // check for multiple assignment of different types to one element
          if ((El_type[k] != 0) && (El_type[k] != tmp_type))
          {
            sprintf(errmsg, "Different types have been assigned at element %ld\n"
                    "see input file %s, line %ld, column %ld", k+1, in->fname, aline, acol);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 1;
          }
          if (El_type[k] == tmp_type)
            fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned type at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
          El_type[k] = tmp_type;
          prop_used = 1;
          if ((tmp_type >= planeelementlt) && (tmp_type <= planequadinterface))
          {// check for multiple assignment of different ssst to one element
            if ((El_ssst[k] != 0) && (El_ssst[k] != tmp_ssst))
            {
              sprintf(errmsg, "Different strain/stress states have been assigned at element %ld\n"
                              "see input file %s, line %ld, column %ld", k+1, in->fname, aline, acol);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              return 2;
            }
            if (El_ssst[k] == tmp_ssst)
              fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned ssst at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
            El_ssst[k] = tmp_ssst;
            // backup of line and column for eventual log message
            line[k]     = aline;
            col[k]      = acol;
          }
        }        
        if (sres < 0)
        {
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 4;
        }
      }
      if (Check_unused_prop && (prop_used == 0))
      {
        // no elements were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Elements with property %ld belonging to entity %s have not been found.\n"
                        " Element type (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 5;
      }
    }
  }

  // check whether all elements have assigned a type
  for(i=0; i<Top->ne; i++)
  {
    if (El_type[i] == 0)
    {
      sprintf(errmsg, "Element %ld has not assigned a type", i+1);
      print_err(errmsg, __FILE__, __LINE__, __func__);
      return 3;
    }
  }

  return 0;
}



/**
  The function assigns material type to the elements. It scans 
  sections eledgpr, elsurfpr, elvolpr for keyword "el_type".
  The record which assigns material type to the element looks like follows:
  "el_mat" "propid" prop "num_mat" nmat nmat*{"type" material_type "type_id" material_id}
  where prop is positive integer number.
  
  The asssigning starts at section for volumes and finishes at
  section for edges. If an element  belongs to several entities and different 
  material types are prescribed for these entities, error of assignment is signalized.
  Finally, the function checks whether all elements have assigned a material type.
 
  @param in  - pointer to opened input file with property description
  @elemsects - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - different material types have been assigned to one element
  @retval 2 - no type has been assigned to elements
  @retval 3 - property numbers of the required entity type have not been read on elements
  @retval 4 - property number of the required entity type have not been found on elements
*/
long input_elem_mat(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->ne), col(Top->ne);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k, l, err;
  mattype *tmp_type;
  long *tmp_id, *ic;
  long prop, nmat;
  long sres;
  long prop_used;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of element material types");
  El_nmat = new long[Top->ne];
  memset(El_nmat, 0, sizeof(*El_nmat)*Top->ne);
  El_mattype = new mattype*[Top->ne];
  memset(El_mattype, 0, sizeof(*El_mattype)*Top->ne);
  El_matid = new long*[Top->ne];
  memset(El_matid, 0, sizeof(*El_matid)*Top->ne);
  El_matdbi = new long*[Top->ne];
  memset(El_matdbi, 0, sizeof(*El_matdbi)*Top->ne);

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "el_mat"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "el_mat", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    { 
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "el_mat", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld%k%ld", "el_mat", "propid", &prop, "num_mat", &nmat);
      if (nmat<=0)
      {
        print_err("number of material types is <= 0", __FILE__, __LINE__, __func__);
        return 1;
      }
      else
      {
        tmp_type = new mattype[nmat];
        memset(tmp_type, 0, sizeof(*tmp_type)*nmat);
        tmp_id   = new long[nmat];
        memset(tmp_id, 0, sizeof(*tmp_id)*nmat);
        ic       = new long[nmat];
        memset(ic, 0, sizeof(*ic)*nmat);
      }
      for (k=0; k<nmat; k++)
      {
        xfscanf(in, "%k%m%k%ld", "type", &mattype_kwdset, tmp_type+k, "type_id", tmp_id+k);
        // searching of prescribed material type and index in the material database
        ic[k] = Dbmat->search_mat(tmp_type[k], tmp_id[k]-1);
        if (ic[k] < 0)
        {
          sprintf(errmsg, "Line %ld, col %ld, file %s: cannot find %ld. required material in the database", in->line, in->col, in->fname, k+1);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          delete [] tmp_type;
          delete [] tmp_id;
          delete [] ic;
          return 1;
        }
        else
          Dbmat->mark_used(ic[k], tmp_id[k]-1);
      }

      in->kwdmode = bkwdmode;
      
      prop_used = 0;
      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, NULL, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { 
          prop_used = 1;
          if (El_nmat[k] != 0)
          {// check for multiple assignment of different material types to one element
            err = 0;
            if (nmat != El_nmat[k])
              err = 1;
            for (l=0; (l<nmat) && (err==0); l++)
            {
              if ((tmp_type[l] != El_mattype[k][l]) || (tmp_id[l] != El_matid[k][l]))
                err = 1;
            }        
            if (err)           
            {
              sprintf(errmsg, "Different material types have been assigned at element %ld\n"
                      "see input file %s, line %ld, column %ld", k+1, in->fname, aline, acol);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete [] tmp_type;
              delete [] tmp_id;
              delete [] ic;
              return 2;
            }
            else
              fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned material type at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
          }
          else
          {
            El_mattype[k] = new mattype[nmat];
            El_matid[k]   = new long[nmat];
            El_matdbi[k] = new long[nmat];
            El_nmat[k]    = nmat;
            for(l=0; l<nmat; l++)
            {
              El_mattype[k][l] = tmp_type[l];
              El_matid[k][l]   = tmp_id[l];
              El_matdbi[k][l]  = ic[l];
            }
            // backup of line and column for eventual log message
            line[k]     = aline;
            col[k]      = acol;
          }
        }
        if (sres < 0)
        {
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          delete [] tmp_type;
          delete [] tmp_id;
          delete [] ic;
          return 3;
        }
      }
      delete [] tmp_type;
      delete [] tmp_id;
      delete [] ic;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no elements were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Elements with property %ld belonging to entity %s have not been found.\n"
                        " Element material type (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 4;
      }
    }
  }

  // check whether all elements have assigned a type
  for(i=0; i<Top->ne; i++)
  {
    if (El_mattype[i] == NULL)
    {
      sprintf(errmsg, "Element %ld has not assigned a material type", i+1);
      print_err(errmsg, __FILE__, __LINE__, __func__);
      return 2;
    }
  }

  return 0;
}



/**
  The function assigns cross-section type to the elements. It scans 
  sections eledgpr, elsurfpr, elvolpr for keyword "el_type".
  The record which assigns cross-section type to the element looks like follows:
  "el_crsec" "propid" prop crsec_type crsec_id
  where prop is positive integer number.
  
  The asssigning starts at section for volumes and finishes at
  section for edges. If an element  belongs to several entities and different 
  cross-section types are prescribed for these entities, error of assignment is signalized.
 
  @param in  - pointer to opened input file with property description
  @elemsects - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  Returns :
  @retval 0 - on succes 
  @retval 1 - cannot find required cross-section in database 
  @retval 2 - different cross-section types have been assigned to one element
  @retval 3 - property numbers of the required entity type have not been read on elements
  @retval 4 - property number of the required entity type have not been found on elements
*/
long input_elem_crsec(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->ne), col(Top->ne);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k, ic;
  crsectype tmp_type;
  long prop, tmp_id;
  long sres;
  long prop_used;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of element cross-section types");
  El_cst = new crsectype[Top->ne];
  memset(El_cst, 0, sizeof(*El_cst)*Top->ne);
  El_csti = new long[Top->ne];
  memset(El_csti, 0, sizeof(*El_csti)*Top->ne);
  El_cstdbi = new long[Top->ne];
  memset(El_cstdbi, 0, sizeof(*El_cstdbi)*Top->ne);

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "el_crsec"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "el_crsec", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "el_crsec", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld%k%m%k%ld", "el_crsec", "propid", &prop, "type", &crsectype_kwdset, &tmp_type, "type_id", &tmp_id);      
      // searching of prescribed cross-section type and index in the cross-section database
      ic = Dbcrs->search_crs(tmp_type, tmp_id-1);
      if (ic < 0)
      {
        sprintf(errmsg, "Line %ld, col %ld, file %s: cannot find required cross-section in the database", in->line, in->col, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 1;
      }
      else
        Dbcrs->mark_used(ic, tmp_id-1);

      in->kwdmode = bkwdmode;

      prop_used = 0;
      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, NULL, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { // check for multiple assignment of different cross-section types to one element
          if ((El_cst[k] != 0) && ((El_cst[k] != tmp_type) || (El_csti[k] != tmp_id)))
          {
            sprintf(errmsg, "Different cross-section types have been assigned at element %ld\n"
                    "see input file %s, line %ld, column %ld", k+1, in->fname, aline, acol);
            print_err(errmsg, __FILE__, __LINE__, __func__);
            return 2;
          }
          prop_used = 1;
          if (El_cst[k] == tmp_type)
            fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned cross-section type at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
          else
          {
            El_cst[k] = tmp_type;
            El_csti[k]   = tmp_id;
            El_cstdbi[k]  = ic;
            // backup of line and column for eventual log message
            line[k]     = aline;
            col[k]      = acol;
          }
        }
        if (sres < 0)
        {
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }

      }
      if (Check_unused_prop && (prop_used == 0))
      {
        // no elements were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Elements with property %ld belonging to entity %s have not been found.\n"
                        " Element cross-section type (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 4;
      }
    }
  }

  return 0;
}



/**
  Function reads local coordinate systems defined at elements.
  It scans sections elemedgpr, elemsurfpr, elemvolpr for keyword "el_lcs".
  The record which assigns local coordinate system to elements looks like follows:
  "el_lcs" "propid" prop "dim" vect_dim {"basevec" base_vector_components}[vect_dim]

  The asssigning starts at section for volumes and finishes at
  section for edges. In case of multiple assigning of different local coordinate systems,
  error of assignment is signalized.

  @param in  - pointer to opened input file with property description
  @elemsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - invalid dimension of base vectors
  @retval 2 - property numbers of the required entity type have not been read on elements
  @retval 3 - property number of the required entity type have not been found on elements
*/
long input_elem_lcs(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->ne), col(Top->ne);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long dim;
  long i, j, k, l, m;
  long check, prop;
  vector *basev;
  long sres;
  long prop_used;
  long *entid;
  long nentid;

  El_lcs = new vector*[Top->ne];
  memset(El_lcs, 0, sizeof(*El_lcs)*Top->ne);

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of local coordinate systems at elements");

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "el_lcs"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "el_lcs", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "el_lcs", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "el_lcs", "propid", &prop);
      xfscanf(in, "%k%ld", "dim", &dim);
      if ((dim < 2) || (dim > 3))
      {
        sprintf(errmsg, "Line %ld, col %ld, file %s: Wrong dimension of base vectors (%ld)\n"
                        "Dimension has to be in range  <2,3>\n", 
                aline, acol, in->fname, dim);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 1;
      }
      // reading of base vectors
      basev = new vector[dim];
      for (k=0; k<dim; k++)
      {
        allocv(dim, basev[k]);
        xfscanf(in, "%k", "basevec");
        readv(in, basev[k]);
      }     

      in->kwdmode = bkwdmode;

      prop_used = 0; 
      // assigning of lcs to the element
      for (k=0; k<Top->ne; k++)
      {        
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, NULL, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { 
          prop_used = 1;
          // check for multiple assignment of different local coordinate systems to one element
          if (El_lcs[k] != NULL)
          { 
            check = 1; 
            for (l=0; l<dim; l++)
            {
              for(m=0; m<dim; m++)
              {
                if (El_lcs[k][l][m] != basev[l][m])
                {
                  check = 0;
                  l = dim;
                  break;
                }
              }
            }
            if (check == 0)
            {
              // element has assigned different local coordinet systems
              sprintf(errmsg, "Line %ld, col %ld, file %s: Element %ld has already had assigned different"
                              "local coordinate system at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, line[k], col[k], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete [] basev;
              return 1;
            }
          }
          else
          {
            El_lcs[k]  = new vector[dim];
            for (l=0; l<dim; l++)
            {
              allocv(dim, El_lcs[k][l]);
              copyv(basev[l], El_lcs[k][l]);
            } 
            // backup of line and column for eventual log message
            line[k] = aline;
            col[k]  = acol;
          }
        }
        if (sres < 0)
        {
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 2;
        }
      }
      if (Check_unused_prop && (prop_used == 0))
      {
        // no elements were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Elements with property %ld belonging to entity %s have not been found.\n"
                        " Element lcs (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 3;
      }
      delete [] basev;
    }
  }
  return 0;
}



/**
  The function assigns different types of load defined at particular elements to the elements 
  with the required property.  It scans  sections eledgpr, elsurfpr, elvolpr for keyword "el_load".
  The record which assigns volume load to the element looks like follows:
  "el_load" "propid" prop "lc_id" load_case_index "load_type" type_of_load 
  "nedge"|"nsurf"|"" number_of_elem_edges|surfaces| "ncomp" total_number_of_load_componets {load_components}[ncomp] 
  where prop is positive integer number. For more details see loadel.cpp|h
  
  The asssigning starts at section for volumes and finishes at
  section for edges. If an element  belongs to several entities and different 
  load types are prescribed for these entities, merging of load is performed.
 
  @param in  - pointer to opened input file with property description
  @elemsects - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - invalid load/subload case index
  @retval 2 - unknown type of load or load cannot be read
  @retval 3 - merging of load failured
  @retval 4 - property numbers of the required entity type have not been read on elements
  @retval 5 - property number of the required entity type have not been found on elements
*/
long input_elem_load(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k;
  loadel *tmp;
  long lt_id;
  long check;
  long prop, lcid;
  long sres, tnlc;
  long prop_used;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of element load");
  tnlc = Nlc;
  if (Tnslc)
    tnlc = Tnslc;
  El_load = new loadel**[Top->ne];
  memset(El_load, 0, sizeof(*El_load)*Top->ne);
  El_loadln = new long **[Top->ne];
  memset(El_loadln, 0, sizeof(*El_loadln)*Top->ne);
  El_loadcol = new long **[Top->ne];
  memset(El_loadcol, 0, sizeof(*El_loadcol)*Top->ne);
  for(i=0; i<Top->ne; i++)
  {
    El_load[i] = new loadel*[tnlc];
    memset(El_load[i], 0, sizeof(*El_load[i])*tnlc);
    El_loadln[i] = new long*[tnlc];
    memset(El_loadln[i], 0, sizeof(*El_loadln[i])*tnlc);
    El_loadcol[i] = new long*[tnlc];
    memset(El_loadcol[i], 0, sizeof(*El_loadcol[i])*tnlc);
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "el_load"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "el_load", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "el_load", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "el_load", "propid", &prop);
      tmp = new loadel;
      check = tmp->read_prep(in, Nlc, Nslc);
      if (check)
      {
        delete tmp;
        return 1; 
      }
      lt_id = give_lt_id(tmp->tel);
      if (lt_id < 0)
      {
        print_err("Line %ld, col %ld, file %s: unknown type of element load is required",
                  __FILE__, __LINE__, __func__, aline, acol, in->fname);
        delete tmp;
        return 2;
      }
      lcid = tmp->nlc-1;
      if (Tnslc){        
        if (Mp->tprob == forced_dynamics)
          lcid = Nslc_cum[lcid]+tmp->nslc;
        else
          lcid = Nslc_cum[lcid]+tmp->nslc-1;
      }

      in->kwdmode = bkwdmode;

      prop_used = 0;            
      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, NULL, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { // check for multiple assignment of loads to one element
          if (El_load[k][lcid])
          {
            if (El_loadln[k][lcid][lt_id]) 
              fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned load at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
            else
            { // element has assigned different type of load
              El_loadln[k][lcid][lt_id] = aline;
              El_loadcol[k][lcid][lt_id]  = acol;
            }
          }
          else
          {
            El_load[k][lcid] = new loadel;
            // backup of line and column for eventual log message
            El_loadln[k][lcid] = new long[3];
            memset(El_loadln[k][lcid], 0, sizeof(*El_loadln[k][lcid])*3);
            El_loadcol[k][lcid] = new long[3];
            memset(El_loadcol[k][lcid], 0, sizeof(*El_loadcol[k][lcid])*3);
            El_loadln[k][lcid][lt_id] = aline;
            El_loadcol[k][lcid][lt_id]  = acol;
          }
          check = El_load[k][lcid]->merge(*tmp);
          prop_used = 1;
          switch(check)
          {
            case 0:
              break;
            case 1:
              // element has already assigned load with different number of components
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "with incompatible number of components at line %ld, col %ld, file: %s\n", 
                              aline, acol, in->fname, k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              return 3;
            case 2:
              // element has already assigned different type of load
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned different "
                              "type of load at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              return 3;
            default:
              // element has already assigned load which cannot be merged
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "which cannot be merged at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              return 3;
          }
        }        
        if (sres < 0)
        {
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 4;
        }
      }
      delete tmp;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no elements were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Elements with property %ld belonging to entity %s have not been found.\n"
                        " Element general load (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 5;
      }
    }
  }

  return 0;
}



/**
  The function assigns edge load to the elements with with required edge property.
  It scans  sections eledgpr for keyword "edge_load".
  The record which assigns edge load to the element looks like follows:
  "edge_load" "propid" prop "lc_id" load_case_index "ncomp" total_number_of_load_componets {load_components}[ncomp] 
  where prop is positive integer number. For more details see entityload.cpp|h
  
  The asssigning scans only section for edges. If an element  belongs to several entities and different 
  load types are prescribed for these entities, merging of load is performed.
 
  @param in  - pointer to opened input file with property description
  @elemsects - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - read error
  @retval 2 - merging of load failured
  @retval 3 - property numbers of the required entity type have not been read on elements
  @retval 4 - property number of the required entity type have not been found on elements
*/
long input_elem_loadedge(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k;
  entityload *tmp;
  loadel *tmp2;
  long lt_id=0;
  long check;
  long prop, lcid;
  long sres, tnlc;
  long prop_used;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of edge load applied on elements");
  tnlc = Nlc;
  if (Tnslc)
    tnlc = Tnslc;
  if (El_load == NULL)
  {
    El_load = new loadel**[Top->ne];
    memset(El_load, 0, sizeof(*El_load)*Top->ne);
    El_loadln = new long **[Top->ne];
    memset(El_loadln, 0, sizeof(*El_loadln)*Top->ne);
    El_loadcol = new long **[Top->ne];
    memset(El_loadcol, 0, sizeof(*El_loadcol)*Top->ne);
    for(i=0; i<Top->ne; i++)
    {
      El_load[i] = new loadel*[tnlc];
      memset(El_load[i], 0, sizeof(*El_load[i])*tnlc);
      El_loadln[i] = new long*[tnlc];
      memset(El_loadln[i], 0, sizeof(*El_loadln[i])*tnlc);
      El_loadcol[i] = new long*[tnlc];
      memset(El_loadcol[i], 0, sizeof(*El_loadcol[i])*tnlc);
    }
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gcurve; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "el_load"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "edge_load", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "edge_load", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "edge_load", "propid", &prop);
      tmp = new entityload;
      // read load record
      check = tmp->read(in, Nlc, Nslc);
      if (check)
      {
        delete tmp;
        return 1;
      }
      lcid = tmp->nlc-1;
      if (Tnslc){        
        if (Mp->tprob == forced_dynamics)
          lcid = Nslc_cum[lcid]+tmp->nslc;
        else
          lcid = Nslc_cum[lcid]+tmp->nslc-1;
      }

      in->kwdmode = bkwdmode;

      prop_used = 0;            
      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, NULL, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { // check for multiple assignment of loads to one element
          if (El_load[k][lcid])
          {
            if (El_loadln[k][lcid][0])
              fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned load at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
            else
            { // element has assigned different type of load
              El_loadln[k][lcid][lt_id] = aline;
              El_loadcol[k][lcid][lt_id]  = acol;
            }
          }
          else
          {
            El_load[k][lcid] = new loadel;
            // backup of line and column for eventual log message
            El_loadln[k][lcid] = new long[3];
            memset(El_loadln[k][lcid], 0, sizeof(*El_loadln[k][lcid])*3);
            El_loadcol[k][lcid] = new long[3];
            memset(El_loadcol[k][lcid], 0, sizeof(*El_loadcol[k][lcid])*3);
            El_loadln[k][lcid][lt_id] = aline;
            El_loadcol[k][lcid][lt_id]  = acol;
          }
          tmp2 = tmp->edge2loadel(Top->elements[k], prop, Top->nodes, entid, nentid);
          if (tmp2 == NULL)
          {
            print_err("unknown type of load function is required", __FILE__, __LINE__, __func__);
            delete tmp;
            return 1;
          }
          // load merging
          check = El_load[k][lcid]->merge(*tmp2);
          prop_used = 1;
          switch(check)
          {
            case 0:
              break;
            case 1:
              // element has already assigned load with different number of components
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "with incompatible number of components at line %ld, col %ld, file: %s\n", 
                              aline, acol, in->fname, k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
            case 2:
              // element has already assigned different type of load
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned different "
                              "type of load at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
            default:
              // element has already assigned load which cannot be merged
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "which cannot be merged at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
          }
          delete tmp2;
        }        
        if (sres < 0)
        {
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
      }
      delete tmp;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no elements were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Elements with property %ld belonging to entity %s have not been found.\n"
                        " Element edge load (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 4;
      }
    }
  }

  return 0;
}



/**
  The function assigns surface load to the elements with with required surface property.
  It scans  sections elsurfpr for keyword "surf_load".
  The record which assigns surface load to the element looks like follows:
  "surf_load" "propid" prop "lc_id" load_case_index "ncomp" total_number_of_load_componets {load_components}[ncomp] 
  where prop is positive integer number. For more details see entityload.cpp|h
  
  The asssigning scans only section for surfaces. If an element  belongs to several entities and different 
  load types are prescribed for these entities, merging of load is performed.
 
  @param in  - pointer to opened input file with property description
  @elemsects - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - read error
  @retval 2 - merging of load failured
  @retval 3 - property numbers of the required entity type have not been read on elements
  @retval 4 - property number of the required entity type have not been found on elements
*/
long input_elem_loadsurf(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k;
  entityload *tmp;
  loadel *tmp2;
  long lt_id=1;
  long check;
  long prop, lcid;
  long sres, tnlc;
  long prop_used;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of surface load applied on elements");
  tnlc = Nlc;
  if (Tnslc)
    tnlc = Tnslc;
  if (El_load == NULL)
  {
    El_load = new loadel**[Top->ne];
    memset(El_load, 0, sizeof(*El_load)*Top->ne);
    El_loadln = new long **[Top->ne];
    memset(El_loadln, 0, sizeof(*El_loadln)*Top->ne);
    El_loadcol = new long **[Top->ne];
    memset(El_loadcol, 0, sizeof(*El_loadcol)*Top->ne);
    for(i=0; i<Top->ne; i++)
    {
      El_load[i] = new loadel*[tnlc];
      memset(El_load[i], 0, sizeof(*El_load[i])*tnlc);
      El_loadln[i] = new long*[tnlc];
      memset(El_loadln[i], 0, sizeof(*El_loadln[i])*tnlc);
      El_loadcol[i] = new long*[tnlc];
      memset(El_loadcol[i], 0, sizeof(*El_loadcol[i])*tnlc);
    }
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gsurface; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "el_load"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "surf_load", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "surf_load", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "surf_load", "propid", &prop);
      tmp = new entityload;
      // read load record
      check = tmp->read(in, Nlc, Nslc);
      if (check)
      {
        delete tmp;
        return 1;
      }
      lcid = tmp->nlc-1;
      if (Tnslc){        
        if (Mp->tprob == forced_dynamics)
          lcid = Nslc_cum[lcid]+tmp->nslc;
        else
          lcid = Nslc_cum[lcid]+tmp->nslc-1;
      }

      in->kwdmode = bkwdmode;
            
      prop_used = 0;
      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, NULL, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { // check for multiple assignment of load to one element
          if (El_load[k][lcid])
          {
            if (El_loadln[k][lcid][1])
              fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned load at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
            else
            { // element has assigned different type of load
              El_loadln[k][lcid][lt_id] = aline;
              El_loadcol[k][lcid][lt_id]  = acol;
            }
          }
          else
          {
            El_load[k][lcid] = new loadel;
            // backup of line and column for eventual log message
            El_loadln[k][lcid] = new long[3];
            memset(El_loadln[k][lcid], 0, sizeof(*El_loadln[k][lcid])*3);
            El_loadcol[k][lcid] = new long[3];
            memset(El_loadcol[k][lcid], 0, sizeof(*El_loadcol[k][lcid])*3);
            El_loadln[k][lcid][lt_id] = aline;
            El_loadcol[k][lcid][lt_id]  = acol;
          }
          tmp2 = tmp->surface2loadel(Top->elements[k], prop, Top->nodes, entid, nentid);
          if (tmp2 == NULL)
          {
            print_err("unknown type of load function is required", __FILE__, __LINE__, __func__);
            delete tmp;
            return 1;
          }
          // load merging
          check = El_load[k][lcid]->merge(*tmp2);
          prop_used = 1;
          switch(check)
          {
            case 0:
              break;
            case 1:
              // element has already assigned load with different number of components
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "with incompatible number of components at line %ld, col %ld, file: %s\n", 
                              aline, acol, in->fname, k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
            case 2:
              // element has already assigned different type of load
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned different "
                              "type of load at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
            default:
              // element has already assigned load which cannot be merged
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "which cannot be merged at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
          }
          delete tmp2;
        }        
        if (sres < 0)
        {
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
      }
      delete tmp;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no elements were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Elements with property %ld belonging to entity %s have not been found.\n"
                        " Element surface load (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 4;
      }
    }
  }

  return 0;
}



/**
  The function assigns volume load to the elements with with required volume property.
  It scans  sections elsurfpr for keyword "volume_load".
  The record which assigns surface load to the element looks like follows:
  "volume_load" "propid" prop "lc_id" load_case_index "ncomp" total_number_of_load_componets {load_components}[ncomp] 
  where prop is positive integer number. For more details see entityload.cpp|h
  
  The asssigning scans only section for surfaces. If an element  belongs to several entities and different 
  load types are prescribed for these entities, merging of load is performed.
 
  @param in  - pointer to opened input file with property description
  @param elemsects - array with descriptors of sections, which will be searched
  @nsect nsect     - number of searched sections i.e number of elements in array elemsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - read error
  @retval 2 - merging of load failured
  @retval 3 - property numbers of the required entity type have not been read on elements
  @retval 4 - property number of the required entity type have not been found on elements
*/
long input_elem_loadvol(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k;
  entityload *tmp;
  loadel *tmp2;
  long check;
  long lt_id = 2;
  long prop, lcid;
  long sres, tnlc;
  long prop_used;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of volume load applied on elements");
  tnlc = Nlc;
  if (Tnslc)
    tnlc = Tnslc;
  if (El_load == NULL)
  {
    El_load = new loadel**[Top->ne];
    memset(El_load, 0, sizeof(*El_load)*Top->ne);
    El_loadln = new long **[Top->ne];
    memset(El_loadln, 0, sizeof(*El_loadln)*Top->ne);
    El_loadcol = new long **[Top->ne];
    memset(El_loadcol, 0, sizeof(*El_loadcol)*Top->ne);
    for(i=0; i<Top->ne; i++)
    {
      El_load[i] = new loadel*[tnlc];
      memset(El_load[i], 0, sizeof(*El_load[i])*tnlc);
      El_loadln[i] = new long*[tnlc];
      memset(El_loadln[i], 0, sizeof(*El_loadln[i])*tnlc);
      El_loadcol[i] = new long*[tnlc];
      memset(El_loadcol[i], 0, sizeof(*El_loadcol[i])*tnlc);
    }
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "volume_load"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "volume_load", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "volume_load", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "volume_load", "propid", &prop);
      tmp = new entityload;
      // read load record
      check = tmp->read(in, Nlc, Nslc);
      if (check)
      {
        delete tmp;
        return 1;
      }
      lcid = tmp->nlc-1;
      if (Tnslc){        
        if (Mp->tprob == forced_dynamics)
          lcid = Nslc_cum[lcid]+tmp->nslc;
        else
          lcid = Nslc_cum[lcid]+tmp->nslc-1;
      }

      in->kwdmode = bkwdmode;
            
      prop_used = 0;
      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, NULL, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { // check for multiple assignment of loads to one element
          if (El_load[k][lcid])
          {
            if (El_loadln[k][lcid][2])
              fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned load at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
            else
            { // element has assigned different type of load
              El_loadln[k][lcid][lt_id] = aline;
              El_loadcol[k][lcid][lt_id]  = acol;
            }
          }
          else
          {
            El_load[k][lcid] = new loadel;
            // backup of line and column for eventual log message
            El_loadln[k][lcid] = new long[3];
            memset(El_loadln[k][lcid], 0, sizeof(*El_loadln[k][lcid])*3);
            El_loadcol[k][lcid] = new long[3];
            memset(El_loadcol[k][lcid], 0, sizeof(*El_loadcol[k][lcid])*3);
            El_loadln[k][lcid][lt_id] = aline;
            El_loadcol[k][lcid][lt_id]  = acol;
          }
          tmp2 = tmp->vol2loadel(Top->elements[k], prop, Top->nodes);
          if (tmp2 == NULL)
          {
            print_err("unknown type of load function is required", __FILE__, __LINE__, __func__);
            delete tmp;
            return 1;
          }
          // load merging
          check = El_load[k][lcid]->merge(*tmp2);
          prop_used = 1;
          switch(check)
          {
            case 0:
              break;
            case 1:
              // element has already assigned load with different number of components
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "with incompatible number of components at line %ld, col %ld, file: %s\n", 
                              aline, acol, in->fname, k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
            case 2:
              // element has already assigned different type of load
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned different "
                              "type of load at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
            default:
              // element has already assigned load which cannot be merged
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "which cannot be merged at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_loadln[k][lcid][lt_id], El_loadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
          }
          delete tmp2;
        }
        if (sres < 0)
        {
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
      }
      delete tmp;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no elements were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Elements with property %ld belonging to entity %s have not been found.\n"
                        " Element volume load (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 4;
      }
    }
  }

  return 0;
}



/**
  The function assigns different types of time dependent load defined at particular elements to the elements 
  with the required property.  It scans  sections eledgpr, elsurfpr, elvolpr for keyword "el_tdload".
  The record which assigns volume load to the element looks like follows:
  "el_tdload" "propid" prop "lc_id" load_case_index "load_type" type_of_tdload 
  "nedge"|"nsurf"|"" number_of_elem_edges|surfaces| "ncomp" total_number_of_tdload_componets {load_components}[ncomp] 
  where prop is positive integer number. For more details see dloadel.cpp|h
  
  The asssigning starts at section for volumes and finishes at
  section for edges. If an element  belongs to several entities and different 
  load types are prescribed for these entities, merging of load is performed.
 
  @param in  - pointer to opened input file with property description
  @elemsects - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - invalid load/subload case index
  @retval 2 - unknown type of load or load cannot be read
  @retval 3 - merging of load failured
  @retval 4 - property numbers of the required entity type have not been read on elements
  @retval 5 - property number of the required entity type have not been found on elements
*/
long input_elem_tdload(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k;
  dloadel *tmp;
  long lt_id;
  long check;
  long prop, lcid;
  long sres, tnlc;
  long prop_used;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of time dependent element load");
  tnlc = Nlc;
  El_tdload = new dloadel**[Top->ne];
  memset(El_tdload, 0, sizeof(*El_tdload)*Top->ne);
  El_tdloadln = new long **[Top->ne];
  memset(El_tdloadln, 0, sizeof(*El_tdloadln)*Top->ne);
  El_tdloadcol = new long **[Top->ne];
  memset(El_tdloadcol, 0, sizeof(*El_tdloadcol)*Top->ne);
  for(i=0; i<Top->ne; i++)
  {
    El_tdload[i] = new dloadel*[tnlc];
    memset(El_tdload[i], 0, sizeof(*El_tdload[i])*tnlc);
    El_tdloadln[i] = new long*[tnlc];
    memset(El_tdloadln[i], 0, sizeof(*El_tdloadln[i])*tnlc);
    El_tdloadcol[i] = new long*[tnlc];
    memset(El_tdloadcol[i], 0, sizeof(*El_tdloadcol[i])*tnlc);
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "el_tdload"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "el_tdload", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "el_tdload", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "el_tdload", "propid", &prop);
      tmp = new dloadel;
      check = tmp->read_prep(in, Nlc);
      if (check)
      {
        delete tmp;
        return 1; 
      }
      lt_id = give_lt_id(tmp->tel);
      if (lt_id < 0)
      {
        print_err("Line %ld, col %ld, file %s: unknown type of element load is required",
                  __FILE__, __LINE__, __func__, aline, acol, in->fname);
        delete tmp;
        return 2;
      }
      lcid = tmp->nlc-1;

      in->kwdmode = bkwdmode;

      prop_used = 0;            
      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, NULL, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { // check for multiple assignment of loads to one element
          if (El_tdload[k][lcid])
          {
            if (El_tdloadln[k][lcid][lt_id]) 
              fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned load at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
            else
            { // element has assigned different type of load
              El_tdloadln[k][lcid][lt_id] = aline;
              El_tdloadcol[k][lcid][lt_id]  = acol;
            }
          }
          else
          {
            El_tdload[k][lcid] = new dloadel;
            // backup of line and column for eventual log message
            El_tdloadln[k][lcid] = new long[3];
            memset(El_tdloadln[k][lcid], 0, sizeof(*El_tdloadln[k][lcid])*3);
            El_tdloadcol[k][lcid] = new long[3];
            memset(El_tdloadcol[k][lcid], 0, sizeof(*El_tdloadcol[k][lcid])*3);
            El_tdloadln[k][lcid][lt_id] = aline;
            El_tdloadcol[k][lcid][lt_id]  = acol;
          }
          check = El_tdload[k][lcid]->merge(*tmp);
          prop_used = 1;
          switch(check)
          {
            case 0:
              break;
            case 1:
              // element has already assigned load with different number of components
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "with incompatible number of components at line %ld, col %ld, file: %s\n", 
                              aline, acol, in->fname, k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              return 3;
            case 2:
              // element has already assigned different type of load
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned different "
                              "type of load at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              return 3;
            default:
              // element has already assigned load which cannot be merged
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "which cannot be merged at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              return 3;
          }
        }        
        if (sres < 0)
        {
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 4;
        }
      }
      delete tmp;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no elements were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Elements with property %ld belonging to entity %s have not been found.\n"
                        " Element general load (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 5;
      }
    }
  }

  return 0;
}



/**
  The function assigns time dependent edge load to the elements with with required edge property.
  It scans  sections eledgpr for keyword "edge_tdload".
  The record which assigns edge load to the element looks like follows:
  "edge_tdload" "propid" prop "lc_id" load_case_index "ncomp" total_number_of_tdload_components {load_components}[ncomp] 
  where prop is positive integer number. For more details see entitytdload.cpp|h
  
  The asssigning scans only section for edges. If an element  belongs to several entities and different 
  load types are prescribed for these entities, merging of load is performed.
 
  @param in  - pointer to opened input file with property description
  @elemsects - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - read error
  @retval 2 - merging of load failured
  @retval 3 - property numbers of the required entity type have not been read on elements
  @retval 4 - property number of the required entity type have not been found on elements
*/
long input_elem_tdloadedge(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k;
  entitytdload *tmp;
  dloadel *tmp2;
  long lt_id=0;
  long check;
  long prop, lcid;
  long sres, tnlc;
  long prop_used;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of edge load applied on elements");
  tnlc = Nlc;
  if (El_tdload == NULL)
  {
    El_tdload = new dloadel**[Top->ne];
    memset(El_tdload, 0, sizeof(*El_tdload)*Top->ne);
    El_tdloadln = new long **[Top->ne];
    memset(El_tdloadln, 0, sizeof(*El_tdloadln)*Top->ne);
    El_tdloadcol = new long **[Top->ne];
    memset(El_tdloadcol, 0, sizeof(*El_tdloadcol)*Top->ne);
    for(i=0; i<Top->ne; i++)
    {
      El_tdload[i] = new dloadel*[tnlc];
      memset(El_tdload[i], 0, sizeof(*El_tdload[i])*tnlc);
      El_tdloadln[i] = new long*[tnlc];
      memset(El_tdloadln[i], 0, sizeof(*El_tdloadln[i])*tnlc);
      El_tdloadcol[i] = new long*[tnlc];
      memset(El_tdloadcol[i], 0, sizeof(*El_tdloadcol[i])*tnlc);
    }
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gcurve; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "el_tdload"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "edge_tdload", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "edge_tdload", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "edge_tdload", "propid", &prop);
      tmp = new entitytdload;
      // read load record
      check = tmp->read(in, Nlc);
      if (check)
      {
        delete tmp;
        return 1;
      }
      lcid = tmp->nlc-1;

      in->kwdmode = bkwdmode;

      prop_used = 0;            
      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, NULL, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { // check for multiple assignment of loads to one element
          if (El_tdload[k][lcid])
          {
            if (El_tdloadln[k][lcid][0])
              fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned load at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
            else
            { // element has assigned different type of load
              El_tdloadln[k][lcid][lt_id] = aline;
              El_tdloadcol[k][lcid][lt_id]  = acol;
            }
          }
          else
          {
            El_tdload[k][lcid] = new dloadel;
            // backup of line and column for eventual log message
            El_tdloadln[k][lcid] = new long[3];
            memset(El_tdloadln[k][lcid], 0, sizeof(*El_tdloadln[k][lcid])*3);
            El_tdloadcol[k][lcid] = new long[3];
            memset(El_tdloadcol[k][lcid], 0, sizeof(*El_tdloadcol[k][lcid])*3);
            El_tdloadln[k][lcid][lt_id] = aline;
            El_tdloadcol[k][lcid][lt_id]  = acol;
          }
          tmp2 = tmp->tdedge2dloadel(Top->elements[k], prop, Top->nodes, entid, nentid);
          if (tmp2 == NULL)
          {
            print_err("unknown type of load function is required", __FILE__, __LINE__, __func__);
            delete tmp;
            return 1;
          }
          // load merging
          check = El_tdload[k][lcid]->merge(*tmp2);
          prop_used = 1;
          switch(check)
          {
            case 0:
              break;
            case 1:
              // element has already assigned load with different number of components
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "with incompatible number of components at line %ld, col %ld, file: %s\n", 
                              aline, acol, in->fname, k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
            case 2:
              // element has already assigned different type of load
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned different "
                              "type of load at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
            default:
              // element has already assigned load which cannot be merged
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "which cannot be merged at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
          }
          delete tmp2;
        }        
        if (sres < 0)
        {
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
      }
      delete tmp;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no elements were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Elements with property %ld belonging to entity %s have not been found.\n"
                        " Element edge load (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 4;
      }
    }
  }

  return 0;
}



/**
  The function assigns surface load to the elements with with required surface property.
  It scans  sections elsurfpr for keyword "surface_tdload".
  The record which assigns surface load to the element looks like follows:
  "surface_tdload" "propid" prop "lc_id" load_case_index "ncomp" total_number_of_tdload_componets {load_components}[ncomp] 
  where prop is positive integer number. For more details see entitytdload.cpp|h
  
  The asssigning scans only section for surfaces. If an element  belongs to several entities and different 
  load types are prescribed for these entities, merging of load is performed.
 
  @param in  - pointer to opened input file with property description
  @elemsects - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - read error
  @retval 2 - merging of load failured
  @retval 3 - property numbers of the required entity type have not been read on elements
  @retval 4 - property number of the required entity type have not been found on elements
*/
long input_elem_tdloadsurf(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k;
  entitytdload *tmp;
  dloadel *tmp2;
  long lt_id=1;
  long check;
  long prop, lcid;
  long sres, tnlc;
  long prop_used;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of surface load applied on elements");
  tnlc = Nlc;
  if (El_tdload == NULL)
  {
    El_tdload = new dloadel**[Top->ne];
    memset(El_tdload, 0, sizeof(*El_tdload)*Top->ne);
    El_tdloadln = new long **[Top->ne];
    memset(El_tdloadln, 0, sizeof(*El_tdloadln)*Top->ne);
    El_tdloadcol = new long **[Top->ne];
    memset(El_tdloadcol, 0, sizeof(*El_tdloadcol)*Top->ne);
    for(i=0; i<Top->ne; i++)
    {
      El_tdload[i] = new dloadel*[tnlc];
      memset(El_tdload[i], 0, sizeof(*El_tdload[i])*tnlc);
      El_tdloadln[i] = new long*[tnlc];
      memset(El_tdloadln[i], 0, sizeof(*El_tdloadln[i])*tnlc);
      El_tdloadcol[i] = new long*[tnlc];
      memset(El_tdloadcol[i], 0, sizeof(*El_tdloadcol[i])*tnlc);
    }
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gsurface; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "el_tdload"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "surf_tdload", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "surf_tdload", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "surf_tdload", "propid", &prop);
      tmp = new entitytdload;
      // read load record
      check = tmp->read(in, Nlc);
      if (check)
      {
        delete tmp;
        return 1;
      }
      lcid = tmp->nlc-1;

      in->kwdmode = bkwdmode;
            
      prop_used = 0;
      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, NULL, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { // check for multiple assignment of load to one element
          if (El_tdload[k][lcid])
          {
            if (El_tdloadln[k][lcid][1])
              fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned load at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
            else
            { // element has assigned different type of load
              El_tdloadln[k][lcid][lt_id] = aline;
              El_tdloadcol[k][lcid][lt_id]  = acol;
            }
          }
          else
          {
            El_tdload[k][lcid] = new dloadel;
            // backup of line and column for eventual log message
            El_tdloadln[k][lcid] = new long[3];
            memset(El_tdloadln[k][lcid], 0, sizeof(*El_tdloadln[k][lcid])*3);
            El_tdloadcol[k][lcid] = new long[3];
            memset(El_tdloadcol[k][lcid], 0, sizeof(*El_tdloadcol[k][lcid])*3);
            El_tdloadln[k][lcid][lt_id] = aline;
            El_tdloadcol[k][lcid][lt_id]  = acol;
          }
          tmp2 = tmp->tdsurface2dloadel(Top->elements[k], prop, Top->nodes, entid, nentid);
          if (tmp2 == NULL)
          {
            print_err("unknown type of load function is required", __FILE__, __LINE__, __func__);
            delete tmp;
            return 1;
          }
          // load merging
          check = El_tdload[k][lcid]->merge(*tmp2);
          prop_used = 1;
          switch(check)
          {
            case 0:
              break;
            case 1:
              // element has already assigned load with different number of components
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "with incompatible number of components at line %ld, col %ld, file: %s\n", 
                              aline, acol, in->fname, k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
            case 2:
              // element has already assigned different type of load
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned different "
                              "type of load at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
            default:
              // element has already assigned load which cannot be merged
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "which cannot be merged at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
          }
          delete tmp2;
        }        
        if (sres < 0)
        {
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
      }
      delete tmp;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no elements were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Elements with property %ld belonging to entity %s have not been found.\n"
                        " Element surface load (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 4;
      }
    }
  }

  return 0;
}



/**
  The function assigns volume load to the elements with with required volume property.
  It scans  sections elsurfpr for keyword "volume_tdload".
  The record which assigns surface load to the element looks like follows:
  "volume_tdload" "propid" prop "lc_id" load_case_index "ncomp" total_number_of_tdload_componets {load_components}[ncomp] 
  where prop is positive integer number. For more details see entityload.cpp|h
  
  The asssigning scans only section for surfaces. If an element  belongs to several entities and different 
  load types are prescribed for these entities, merging of load is performed.
 
  @param in  - pointer to opened input file with property description
  @param elemsects - array with descriptors of sections, which will be searched
  @nsect nsect     - number of searched sections i.e number of elements in array elemsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - read error
  @retval 2 - merging of load failured
  @retval 3 - property numbers of the required entity type have not been read on elements
  @retval 4 - property number of the required entity type have not been found on elements
*/
long input_elem_tdloadvol(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k;
  entitytdload *tmp;
  dloadel *tmp2;
  long check;
  long lt_id = 2;
  long prop, lcid;
  long sres, tnlc;
  long prop_used;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of volume load applied on elements");
  tnlc = Nlc;
  if (El_tdload == NULL)
  {
    El_tdload = new dloadel**[Top->ne];
    memset(El_tdload, 0, sizeof(*El_tdload)*Top->ne);
    El_tdloadln = new long **[Top->ne];
    memset(El_tdloadln, 0, sizeof(*El_tdloadln)*Top->ne);
    El_tdloadcol = new long **[Top->ne];
    memset(El_tdloadcol, 0, sizeof(*El_tdloadcol)*Top->ne);
    for(i=0; i<Top->ne; i++)
    {
      El_tdload[i] = new dloadel*[tnlc];
      memset(El_tdload[i], 0, sizeof(*El_tdload[i])*tnlc);
      El_tdloadln[i] = new long*[tnlc];
      memset(El_tdloadln[i], 0, sizeof(*El_tdloadln[i])*tnlc);
      El_tdloadcol[i] = new long*[tnlc];
      memset(El_tdloadcol[i], 0, sizeof(*El_tdloadcol[i])*tnlc);
    }
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "volume_tdload"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "volume_tdload", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "volume_tdload", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "volume_tdload", "propid", &prop);
      tmp = new entitytdload;
      // read load record
      check = tmp->read(in, Nlc);
      if (check)
      {
        delete tmp;
        return 1;
      }
      lcid = tmp->nlc-1;

      in->kwdmode = bkwdmode;
            
      prop_used = 0;
      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, NULL, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { // check for multiple assignment of loads to one element
          if (El_tdload[k][lcid])
          {
            if (El_tdloadln[k][lcid][2])
              fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned load at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
            else
            { // element has assigned different type of load
              El_tdloadln[k][lcid][lt_id] = aline;
              El_tdloadcol[k][lcid][lt_id]  = acol;
            }
          }
          else
          {
            El_tdload[k][lcid] = new dloadel;
            // backup of line and column for eventual log message
            El_tdloadln[k][lcid] = new long[3];
            memset(El_tdloadln[k][lcid], 0, sizeof(*El_tdloadln[k][lcid])*3);
            El_tdloadcol[k][lcid] = new long[3];
            memset(El_tdloadcol[k][lcid], 0, sizeof(*El_tdloadcol[k][lcid])*3);
            El_tdloadln[k][lcid][lt_id] = aline;
            El_tdloadcol[k][lcid][lt_id]  = acol;
          }
          tmp2 = tmp->tdvol2dloadel(Top->elements[k], prop, Top->nodes);
          if (tmp2 == NULL)
          {
            print_err("unknown type of load function is required", __FILE__, __LINE__, __func__);
            delete tmp;
            return 1;
          }
          // load merging
          check = El_tdload[k][lcid]->merge(*tmp2);
          prop_used = 1;
          switch(check)
          {
            case 0:
              break;
            case 1:
              // element has already assigned load with different number of components
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "with incompatible number of components at line %ld, col %ld, file: %s\n", 
                              aline, acol, in->fname, k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
            case 2:
              // element has already assigned different type of load
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned different "
                              "type of load at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
            default:
              // element has already assigned load which cannot be merged
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned load "
                              "which cannot be merged at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_tdloadln[k][lcid][lt_id], El_tdloadcol[k][lcid][lt_id], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tmp;
              delete tmp2;
              return 1;
          }
          delete tmp2;
        }
        if (sres < 0)
        {
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 3;
        }
      }
      delete tmp;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no elements were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Elements with property %ld belonging to entity %s have not been found.\n"
                        " Element volume load (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 4;
      }
    }
  }

  return 0;
}



/**
  The function assigns eigenstrains/eigenstresses to the elements with with required property.
  It scans sections elemedgpr, elemsurfpr, elemvolpr for keyword "el_eigstr".
  The record which assigns eigenstrains to elements looks like follows:
  "el_eigstr" "propid" prop "str_type" strt "ncomp" nc "eigstr_comp" {time_functions_of_eigenstr_components}xnc

  The asssigning starts at section for volumes and finishes at
  section for edges. In case of multiple assigning of different eigen strains,
  error of assignment is signalized.
   
  @param in  - pointer to opened input file with property description
  @param elemsects - array with descriptors of sections, which will be searched
  @nsect nsect     - number of searched sections i.e number of elements in array elemsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - read error of gfunct
  @retval 2 - multiple assignment of nonidentical eigenstrains
  @retval 3 - incompatible number of assigned eigenstrain componennts with element
  @retval 4 - property numbers of the required entity type have not been read on elements
  @retval 5 - property number of the required entity type have not been found on elements
  @retval 6 - assignment eigenstresses with eigenstrains at one problem is not supported
*/
long input_elem_eigstr(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  ivector line(Top->ne), col(Top->ne);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k, l;
  gfunct *tmp;
  gfunct *tmp2;
  long ncomp, rncomp;
  long check;
  long *eigstrgf_id;
  long prop;
  long sres;
  long prop_used;
  long *entid;
  long nentid;
  strastre strt, strtp;

  fprintf(stdout, "\n reading of eigenstrains applied on elements");

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;
  El_eigstr = NULL;
  El_eigstrt = NULL;
  strtp = strastre(-1);

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "el_eigstr"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "el_eigstr", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      if (El_eigstr == NULL)
      {
        El_eigstr = new long*[Top->ne];
        memset(El_eigstr, 0, sizeof(*El_eigstr)*Top->ne);
        El_eigstrt = new strastre[Top->ne];
        memset(El_eigstrt, 0, sizeof(*El_eigstrt)*Top->ne);
      }
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "el_eigstr", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      // read property id and number of eigenstrain/eigenstress components
      xfscanf(in, "%k%k%ld%k%m%k%ld%k", "el_eigstr", "propid", &prop, "str_type", &strastre_kwdset, &strt, "ncomp", &ncomp, "eigstr_comp");
      if (strtp < 0)
        strtp = strt;
      if (strt != strtp)
      {
        print_err("Line %ld, col %ld, file %s: Assignment of eigen%s is required\n"
                  " but eigen%s has been already defined at line %ld, col %ld, file: %s\n"
                  "The compound usage of eigenstresses and eigenstrains is not supported", 
                  __FILE__, __LINE__, __func__, in->line, in->col, in->fname, strastre_kwdset.get_str(strt), 
                  strastre_kwdset.get_str(strtp), aline, acol, in->fname);
        return 6;
      }
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;

      eigstrgf_id = new long[ncomp];
      memset(eigstrgf_id, 0, sizeof(*eigstrgf_id)*ncomp);
      for (k=0; k<ncomp; k++)
      {
        tmp = new gfunct;
        check = tmp->read(in);
        if (check != 0)
        {
          delete tmp;
          return 1;
        }
        eigstrgf_id[k] = El_eigstrgf_lst.count();
        El_eigstrgf_lst.append(tmp);
      }
            
      in->kwdmode = bkwdmode;

      prop_used = 0;
      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, NULL, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { // check for multiple assignment of eigenstrains/eigenstresses to one element
          if (El_eigstr[k])
          {
            fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned the same eigenstrains"
                         "at line %ld, col %ld, file: %s\n", 
                         aline, acol, in->fname, k+1, line[k], col[k], in->fname);

            // check for nonidentical eigenstrains/eigenstresses
            if (El_eigstrt[k] != strt)
            {
              print_err("Multiple assignment of nonidentical eigen%s component %ld for element %ld\n" 
                        "at line %ld, col %ld of file '%s' (original component assigned at line %ld, col %ld)", 
                        __FILE__, __LINE__, __func__, strastrestate_kwdset.get_str(El_eigstrt[k]), l+1, k+1, aline, acol, in->fname, line[k], col[k]);
              return 2;
            }
            for(l=0; l<ncomp; l++)
            {
              tmp  = (gfunct *)El_eigstrgf_lst.at(l);              // original gfunct for the eigenstrain component
              tmp2 = (gfunct *)El_eigstrgf_lst.at(eigstrgf_id[l]); // new gfunct for the eigenstrain component

              if (tmp->compare(*tmp2)) 
              {// gfunct objects are not identical
                print_err("Multiple assignment of nonidentical eigen%s component %ld for element %ld\n" 
                          "at line %ld, col %ld of file '%s' (original component assigned at line %ld, col %ld)", 
                          __FILE__, __LINE__, __func__, strastrestate_kwdset.get_str(strt), l+1, k+1, aline, acol, in->fname, line[k], col[k]);
                return 2;
              }
            }
          }
          else
          { 
            rncomp = Mt->give_tncomp(El_type[k]);
            if (El_ssst[k])    rncomp = 4;
            // check equality of number of strain components on the given element and 
            // number of assigned eigenstrain components
            if (rncomp != ncomp)
            {
              print_err("Element %ld has different number of %s components (%ld)\n" 
                        " than the specified at line %ld, col %ld of file '%s'",
                        __FILE__, __LINE__, __func__, k+1, strastrestate_kwdset.get_str(strt), rncomp, aline, acol, in->fname);
              return 3;
            }
            // element has no assigned eigenstrains
            El_eigstr[k] = new long[ncomp];
            memcpy(El_eigstr[k], eigstrgf_id, sizeof(*El_eigstr[k])*ncomp);
            El_eigstrt[k] = strt;
            line[k]  = aline;
            col[k] = acol;
            prop_used = 1;
          }
        }
        if (sres < 0)
        {
          print_err("Element %ld has not read property numbers of %s,\n"
                    " required property cannot be assigned.\n"
                    " see input file %s, line %ld, column %ld", 
                    __FILE__, __LINE__, __func__, k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          return 4;
        }
      }
      delete [] eigstrgf_id;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no elements were found with the given property and checking of unused properties is on
        print_err("Elements with property %ld belonging to entity %s have not been found.\n"
                  " Element eigen%s (line %ld, column %ld, file %s) cannot be assigned correctly.",
                  __FILE__, __LINE__, __func__, prop, entitypstr[ent-1].alias, strastrestate_kwdset.get_str(strt), aline, acol, in->fname);
        return 5;
      }
    }
  }
  return 0;
}



/**
  Function reads time functions for switching on/off of given elements.
  It scans  sections eledgpr, elsurfpr, elvolpr for keyword "el_tfunct".
  The record which assigns volume load to the element looks like follows:
  "el_tfunct" "propid" prop "el_tfunc" time_function_id
  where prop is positive integer number.

  The asssigning starts at section for volumes and finishes at
  section for vertices. In case of multiple assigning of different time functions,
  error of assignment is signalized.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  Returns :
  @retval 0 - on succes  
  @retval 1 - assignment of different time function at one element
  @retval 2 - property numbers of the required entity type have not been read on elements
  @retval 3 - element have not assigned a time function
  @retval 4 - property number of the required entity type have not been found on elements
*/
long input_elem_eltimefunc(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->ne), col(Top->ne);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long tf_id, prop;
  long i, j, k;
  long sres;
  long prop_used;
  long *entid;
  long nentid;

  El_tfunc = new long[Top->ne];
  memset(El_tfunc, 0, sizeof(*El_tfunc)*Top->ne);

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of dof time functions at elements");

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "nod_tfunc"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "el_tfunc", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "el_tfunc", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "el_tfunc", "propid", &prop);
      xfscanf(in, "%k%ld", "tfunc_id", &tf_id);
      if ((tf_id <= 0) || (tf_id > Gtm->ngf))
      {
        sprintf(errmsg, "Line %ld, col %ld, file %s:\n Time function index has to be in range <1,%ld>\n"
                        "time function at line %ld, col %ld, file: %s\n", aline, acol, in->fname, Gtm->ngf, 
                        aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 3;
      }
      in->kwdmode = bkwdmode;
     
      prop_used = 0;
      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, NULL, k, Top->edges, Top->surfaces, entid, nentid);    
        if (sres > 0)
        { // check for multiple assignment of different time functions to one element
          if (El_tfunc[k] && (El_tfunc[k]!=tf_id))
          {
            print_err("Line %ld, col %ld, file %s:\n Element %ld has already had assigned different "
                      "time function at line %ld, col %ld, file: %s\n", __FILE__, __LINE__, __func__,
                      aline, acol, in->fname, k+1, line[k], col[k], in->fname);
            return 1;
          }
          else
          {
            El_tfunc[k] = tf_id;
            prop_used = 1;
            // backup of line and column for eventual log message
            line[k] = aline;
            col[k]  = acol;
          }
        }
        if (sres < 0)
        {
          print_err("Element %ld has not read property numbers of %s,\n"
                    " required property cannot be assigned.\n"
                    " see input file %s, line %ld, column %ld",
                    __FILE__, __LINE__, __func__,
                    k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          return 2;
        }
      }
      if (Check_unused_prop && (prop_used == 0))
      {
        // no elements were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Elements with property %ld belonging to entity %s have not been found.\n"
                        " Element time function (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 4;
      }
    }
  }
  // in case of growing mechanical problem, all elements have to have 
  // assigned a time function
  if (Mp->tprob == growing_mech_structure)
  {
    for (i=0; i< Top->ne; i++)
    {
      if (El_tfunc[i] < 1)
      {
        sprintf(errmsg, "Element %ld has not assigned a time function.\n"
                        "All elements have to have assigned some time function\n"
                        "in case of growing mechanical problem", i+1);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 3;
      }
    }
  }
  return 0;
}



/**
  The overloaded postfix decrementation operator for enum gentity.
  
  Created by TKo, 09.2010
*/
gentity operator --(gentity &a, int)
{
  gentity old;
  switch (a)
  {
    case gvertex:
      old = a;
      a = gregion;
      return old;
    case gcurve:
      old = a;
      a = gvertex;
      return old;
    case gsurface:
      old = a;
      a = gcurve;
      return old;
    case gregion:
      old = a;
      a = gsurface;
      return old;
    default:
      print_err("decreasing of invalid gentity value is required", __FILE__, __LINE__, __func__);
  }
  return gentity(0);
}



/**
  The overloaded postfix incrementation operator for enum gentity.
  
  Created by TKo, 06.2021
*/
gentity operator ++(gentity &a, int)
{
  gentity old;
  switch (a)
  {
    case gvertex:
      old = a;
      a = gcurve;
      return old;
    case gcurve:
      old = a;
      a = gsurface;
      return old;
    case gsurface:
      old = a;
      a = gregion;
      return old;
    case gregion:
      old = a;
      a = gvertex;
      return old;
    default:
      print_err("increasing of invalid gentity value is required", __FILE__, __LINE__, __func__);
  }
  return gentity(0);
}



/**
  Function converts load type alias to long index with different 
  sequence of load types. Tel has sequence volume=1, edge=2, surface=3 and
  resulting index has following values:
  edge    = 0,  surface = 1,  volume  = 2. 
  The function returns -1 for rest values of tel.

  @param tel - alias of converted load type.
  
  Returns:
  @retval -1 - for unknown or unsupported alias
  @retval  0 - for edge load type
  @retval  1 - for surface load type
  @retval  2 - for volume load type
*/
long give_lt_id(elloadtype tel)
{
  long ret;
  switch(tel)
  {
    case edge:
      ret = 0;
      break;
    case surface:
      ret = 1;
      break;
    case volume:
      ret = 2;
      break;
    default:
      ret = -1;
  }
  return ret;
}


void periodbc_log(long cmdid, XFILE *in, long aline, long acol, ivector &setsnodes, ivector &setmnodes, double rerr)
{
  long i,j;
  char fname[50];
  sprintf(fname, "periodbc.%ld.log", cmdid);
  FILE *out = fopen(fname, "wt");
  fprintf(out, "%ld. nod_periodbc command at file %s, line %ld, col %ld:\n\n\n", cmdid, in->fname, aline, acol);

  // print set of master nodes given by property id in nod_periodbc command
  long nmn = setmnodes.n;
  fprintf (out, "Set of master nodes by prop (n=%ld):\n", nmn);
  for(i=0; i<nmn/10; i++){
    for(j=0; j<10; j++)
      fprintf(out, "%6ld ", setmnodes(i*10+j)+1);
    fprintf(out, "\n");
  }
  for(j=0; j<nmn%10; j++)
    fprintf(out, "%6ld ", setmnodes(i*10+j)+1);
  fprintf(out, "\n\n");

  // print set of slave nodes given by property id in nod_periodbc command
  long nsn = setmnodes.n;
  fprintf (out, "Set of slave nodes by prop (n=%ld):\n", nsn);
  for(i=0; i<nsn/10; i++){
    for(j=0; j<10; j++)
      fprintf(out, "%6ld ", setsnodes(i*10+j)+1);
    fprintf(out, "\n");
  }
  for(j=0; j<nsn%10; j++)
    fprintf(out, "%6ld ", setsnodes(i*10+j)+1);
  fprintf(out, "\n\n");

  long nmultn=0, nfn=0;
  nmn = nsn = 0;
  for (i=0; i<Top->nn; i++){
    if (Nod_periodbc[i].slmas == 2){
      nmultn++;
      if (Nod_periodbc[i].cerr > rerr)
        nfn++;
      continue;
    }
    if (Nod_periodbc[i].slmas == 1){
      nsn++;
    }
    if (Nod_periodbc[i].slmas == -1){
      nmn++;
    }
  }
  fprintf(out, "The number of master nodes detected in Nod_periodbc array: %ld\n", nmn);
  fprintf(out, "The number of slave nodes detected in Nod_periodbc array: %ld\n", nsn);
  fprintf(out, "The number of multiple slave/master nodes detected in Nod_periodbc array: %ld\n", nmultn);
  fprintf(out, "The number of failed multiple slave/master nodes with detected in Nod_periodbc array: %ld\n", nfn);
  
  
  // print set of master nodes detected in nod_periodbc command
  fprintf (out, "Set of master nodes detected in Nod_periodbc array (n=%ld):\n", nmn);
  nmn = 0;
  for(i=0; i<Top->nn; i++){
    if (Nod_periodbc[i].slmas == -1){
      fprintf(out, "%6ld", i+1);
      nmn++;  
      if (nmn%10 == 0)
        fprintf(out, "\n");
    }
  }
  fprintf(out, "\n\n");

  // print set of slave nodes detected in nod_periodbc command
  fprintf (out, "Set of slave nodes detected in Nod_periodbc array (n=%ld):\n", nsn);
  nsn = 0;
  for(i=0; i<Top->nn; i++){
    if (Nod_periodbc[i].slmas == 1){
      long eid = Nod_periodbc[i].eid;
      fprintf(out, "%6ld, elem %ld, ", i+1, eid+1);
      long nne = Top->elements[eid].nne;
      for(j=0; j<nne; j++){
        fprintf(out, "%6ld ", Top->elements[eid].nodes[j]+1);
      }
      fprintf(out, ".\n");
    }
  }
  fprintf(out, "\n\n");

  // print set of multiple master/slave nodes 
  fprintf (out, "Set of multiple slave/master nodes detected in Nod_periodbc array (n=%ld):\n", nmultn);
  for(i=0; i<Top->nn; i++){
    if ((Nod_periodbc[i].slmas == 2) && (Nod_periodbc[i].cerr < rerr)){
      long eid = Nod_periodbc[i].eid;
      fprintf(out, "%6ld, elem %ld, ", i+1, eid+1);
      long nne = Top->elements[eid].nne;
      for(j=0; j<nne; j++){
        fprintf(out, "%6ld ", Top->elements[eid].nodes[j]+1);
      }
      fprintf(out, ".\n");
    }
  }
  fprintf(out, "\n\n");


  // print set of failed multiple master/slave nodes 
  fprintf (out, "Set of multiple slave/master nodes detected in Nod_periodbc array (n=%ld):\n", nmultn);
  for(i=0; i<Top->nn; i++){
    if ((Nod_periodbc[i].slmas == 2) && (Nod_periodbc[i].cerr > rerr)){
      long eid = Nod_periodbc[i].eid;
      fprintf(out, "%6ld, eid %ld, err %le\n, ", i+1, eid+1, Nod_periodbc[i].cerr);
    }
  }
  fprintf(out, "\n\n");
  fclose(out);
}



/**
  The function allocates and fills array of master nodes on the given hanging node record array.
  It also attempts to reduce master node list of the given hanging nodes.
  If some of the natural coordinate equals the limit value, i.e. +/- 1.0 then
  the master nodes with the oposite natural coodinate do not contribute to the
  to the hanging node quantity approximation.

  @param nn[in] - number of hanging nodes in the array hn
  @param hn[in,out] - array of hanging node records

  @return The function does not return anything but may modify records in the array hn
          in the case of hanging nodes with limit natural coordinate values.

  Created by Tomas Koudelka, 05.2023
*/
void get_masternodes(long nn, hngen *hn)
{
  long i, j, k, eid, nmn, tnmn;
  ivector mnodes, aux, enod;
  double d, tol;
  vector enc(ASTCKVEC(3)), lnc(ASTCKVEC(3));
  
  for (i=0; i<nn; i++){
    eid = hn[i].eid;
    if (eid < 0)  continue; // it is not hanging node
    nmn = Top->elements[eid].nne;
    reallocv(RSTCKIVEC(nmn, mnodes));
    makerefv(enod, Top->elements[eid].nodes, nmn); 
    copyv(enod, mnodes);
    reallocv(RSTCKIVEC(nmn, aux));
    //    tol = hn[i].cerr;
    tol = 1.0e-4;
    // search nodes with the oposite coordinate and remove them from the master node list
    d = fabs(fabs(hn[i].xi)-1.0);    
    if (d < tol){
      fillv(-1, aux);
      // search nodes with oposite xi coordinate
      Top->give_enod_xicoord(eid, -sgn(hn[i].xi), aux);
      // mark them in the list for removal
      for(k=0; k<aux.n; k++){
        for(j=0; j<mnodes.n; j++){
          if (aux[k] == mnodes[j]){
            mnodes[j] = -1;
            break;
          }
        }
      }
    }    
    d = fabs(fabs(hn[i].eta)-1.0);
    if (d < tol){
      // search nodes with oposite eta coordinate
      Top->give_enod_etacoord(eid, -sgn(hn[i].eta), aux);
      // mark them in the list for removal
      for(k=0; k<aux.n; k++){
        for(j=0; j<mnodes.n; j++){
          if (aux[k] == mnodes[j]){
            mnodes[j] = -1;
            break;
          }
        }
      }
    }
    d = fabs(fabs(hn[i].zeta)-1.0);
    if (d < tol){
      fillv(-1, aux);
      // search nodes with oposite zeta coordinate
      Top->give_enod_zetacoord(eid, -sgn(hn[i].zeta), aux);
      // mark them in the list for removal
      for(k=0; k<aux.n; k++){
        for(j=0; j<mnodes.n; j++){
          if (aux[k] == mnodes[j]){
            mnodes[j] = -1;
            break;
          }
        }
      }
    }
    tnmn = 0;
    for(j=0; j<nmn; j++){
      if (mnodes[j] >= 0) tnmn++;
    }
    reallocv(tnmn, hn[i].mnodes);
    k=0;
    for(j=0; j<nmn; j++){
      if (mnodes[j] >= 0){
        hn[i].mnodes[k] = mnodes[j];
        k++;
      }
    }
    if (tnmn != nmn){
      // master node reduction takes a place
      gentity ent = gentity(0);
      long sid, edid;
      hn[i].et = Top->give_elem_ent(eid, tnmn, ent);
      switch (ent){
        case gvertex:
          break;
        case gcurve:
          edid = Top->elements[eid].compare_edg(hn[i].mnodes.a, hn[i].mnodes.n);
          enc(0) = hn[i].xi;  enc(1) = hn[i].eta;  enc(2) = hn[i].zeta;
          Top->transform_natcoord_edg(hn[i].eid, edid, enc, lnc);
          hn[i].xi = lnc(0);  hn[i].eta = lnc(1);  hn[i].zeta = lnc(2);          
          break;
        case gsurface:
          sid = Top->elements[eid].compare_surf(hn[i].mnodes.a, hn[i].mnodes.n);
          enc(0) = hn[i].xi;  enc(1) = hn[i].eta;  enc(2) = hn[i].zeta;
          Top->transform_natcoord_surf(hn[i].eid, sid, enc, lnc);
          hn[i].xi = lnc(0);  hn[i].eta = lnc(1);  hn[i].zeta = lnc(2);          
          break;
        case gregion:
          print_err("invalid entity type (%d) is required in master node reduction (tnmn=%ld, nmn=%ld).\n"
                    " Error detected at hanging node %ld connetcetd with element %ld.\n",
                    __FILE__, __LINE__, __func__, ent, tnmn, nmn, i+1, hn[i].eid+1);
          break;
        default:
          print_err("unknown type of general entity (%d) is required on hanging node %ld connected with element %ld.\n",
                    __FILE__, __LINE__, __func__, ent, i+1, hn[i].eid+1);
          abort();
      }
    }
    else
      // original master nodes will be used
      hn[i].et = Top->elements[hn[i].eid].type;
  }
}
