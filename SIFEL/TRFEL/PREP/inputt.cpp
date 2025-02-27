#include <string.h>
#include "iotools.h"
#include "kwdset.h"
#include "vector.h"
#include "siftop.h"
#include "gfunct.h"
#include "aliast.h"
#include "globalt.h"

#include "globprept.h"
#include "prepalias.h"
#include "inputt.h"
#include "hangnodet.h"
#include "bocont.h"
#include "dbcrst.h"
#include "dbmatt.h"
#include "entitybocon.h"
#include "intdomt.h"
#include "advectvel.h"
#include "aggregator.h"


/**
 Function   reads problem data from file in a writes it to file out in SIFEL format.

 @param in   - pointer to opened input file
 @param outf - pointer to string with name of output file

 @retval 0 - on succes
 @retval 1 - unable to open temporary file
 @retval 2 - if fails reading number of loading cases
 @retval 3 - if fails reading topology
 @retval 4 - if fails reading materials
 @retval 5 - if fails reading cross-sections
 @retval 6 - if fails reading property file
 @retval 7 - if fails writing output file

 Created by TKo, 09.2010
*/
long inputt(XFILE *in, descript *d)
{
  long err;
  XFILE *itop, *ihn;

  // reading of load cases
  err = input_lct(in);
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
  
  
  if (Tp->ssle->prec.pt==boss){
    //  iterative method is preconditioned by the BOSS algorithm
    if (Tp->ssle->prec.agg->impl==2){
      //  the BOSS method is based on the METIS decomposition
      d->paral=2;
    }
  }
  
  if (Tp->ssle->tlinsol==sfeti){
    //  the problem is solved by the FETI method on a single processor
    d->paral=2;
  }


  err = input_siftopt(itop, d);
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
    err = input_hang_nodest(ihn);
    xfclose(ihn);
    if (err)
      return(3);
  }

  // reading material database
  if (d->matsec == no)
    err = input_materialst(d->matf, d);
  else
    err = input_materialst(in, d);
  if (err)
    return(4);

  // reading cross-section elements
  if (d->crssec == no)
    err = input_crst(d->crf, d);
  else
    err = input_crst(in, d);
  if (err)
    return(5);

  // reading of time functions for growing transport problems
  input_lct(in);

  // reading nodal properties
  err = input_nodpropt(in);
  switch (err)
  {
    case 0 :
      break;
    case 1 :
      print_err("reading of nodal boundary conditions failed", __FILE__, __LINE__, __func__);
      return(6);
    case 2 :
      print_err("reading of coupled dofs at nodes failed", __FILE__, __LINE__, __func__);
      return(6);
    case 3 :
      print_err("reading of dof time functions at nodes failed", __FILE__, __LINE__, __func__);
      return(6);
    case 4 :
      print_err("reading of cross-sections at nodes failed", __FILE__, __LINE__, __func__);
      return(6);
    case 5 :
      print_err("reading of nodal initial conditions failed", __FILE__, __LINE__, __func__);
      return(6);
    case 6 :
      print_err("reading of nodal sources failed", __FILE__, __LINE__, __func__);
      return(6);
    case 7 :
      print_err("reading of nodal advection velocities failed", __FILE__, __LINE__, __func__);
      return(6);
    default :
      print_err("unknown error code in nodal sections", __FILE__, __LINE__, __func__);
      return(6);
  }
  // reading nodal properties
  err = input_elempropt(in);
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
      print_err("reading of element source record failed", __FILE__, __LINE__, __func__);
      return(7);
    case 5 :
      print_err("reading of element vertex boundary condition record failed", __FILE__, __LINE__, __func__);
      return(7);
    case 6 :
      print_err("reading of element edge boundary condition record failed", __FILE__, __LINE__, __func__);
      return(7);
    case 7 :
      print_err("reading of element surface boundary condition record failed", __FILE__, __LINE__, __func__);
      return(7);
    case 8 :
      print_err("reading of time functions for elements failed", __FILE__, __LINE__, __func__);
      return(7);
    default :
      print_err("unknown error code in element sections", __FILE__, __LINE__, __func__);
      return(7);
  }
  return(0);
}



/**
  The function reads section with files used for the generation of TRFEL input 
  file e.g., topology, material or cross section files. It also reads setup of 
  preprocessor.

  @param in - pointer to the opened XFILE structure
  @param d  - structure with input data format description
  
  @retval 0 - on success

  Created by Tomas Koudelka, 7.7.2014
*/
long input_filest(XFILE *in, descript &d)
{
  // reading of line with topology file name
  xfscanf(in, " %a", d.topf);
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

  // initial condition file will not be used by default
  // optional flag for input file with initial conditions at nodes
  d.inicdf = no;
  xfscanf(in, "%+k %m", "inicd_file", &answertype_kwdset, &d.inicdf);
  if (d.inicdf == yes)
    // initial condition file will be used
    xfscanf(in, " %a", d.icf);

  // optional reading of hanging nodes
  d.hangnf[0] = 0;
  xfscanf(in, "%+k %1024a", "hanging_nodes_file", &d.hangnf);

  // optional reading of materials via strings or transmat
  d.matstr = yes;
  xfscanf(in, "%+k %m", "read_mat_strings", &answertype_kwdset, &d.matstr);
  if (d.matstr == no)
    xfscanf(in, "%k %m", "read_mat_kwd", &answertype_kwdset, &d.matkwd);
  else 
    d.matkwd = no;

  // optional reading of cross section parameters via strings or transcrsec procedures
  d.crsstr = yes;
  xfscanf(in, "%+k %m", "read_crs_strings", &answertype_kwdset, &d.crsstr);
  if (d.crsstr == no)
    xfscanf(in, "%k %m", "read_crs_kwd", &answertype_kwdset, &d.crskwd);
  else 
    d.crskwd = no;

  return 0;
}



/**
  Input of data about load cases (time functions for elements in 
  growing transport problems).

  @param in - pointer to the opened XFILE structure
  
  @retval 0 - on success
  
  Created 10.2010 by Tomas Koudelka koudelka@cml.fsv.cvut.cz
*/
long input_lct(XFILE *in)
{
  // reading list of time functions controlling dofs and switching of elements on/off
  if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
  {
    fprintf(stdout, "\n\nReading of section with time functions . . .");
    xf_setsec(in, bsec_str[begsec_loadcase]);
    in->kwdmode = sect_mode_full;
    xfscanf(in, "%k", "time_functions");
    in->kwdmode = sequent_mode;
    Gtt->read_gf (in);
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

  Created by TKo, 09.2010
*/
long input_siftopt(XFILE *in, descript *d)
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
long input_hang_nodest(XFILE *in)
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
  
  Nod_hang = new hangnodet*[Top->nn];
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
    Nod_hang[id] = new hangnodet;
    ret = Nod_hang[id]->read(in);
    if (ret)
      return 3;
  }
  return(0);
}



/**
  Function reads material database from file with name given by the fname

  @param fname - string with database file name
  @param d - pointer to structure with description of preprocessor setup 

  Returns :
  @retval 0 - on success
  @retval 1 - if fails opening file
  @retval 2 - if fails reading material parameters

  Created by TKo, 09.2010
*/
long input_materialst(char *fname, descript *d)
{
  long ret;
  long mpb = Mesprt;
  XFILE *in = NULL;

  fprintf(stdout, "\n\nReading of material database form %s . . .", fname);
  in = xfopen(fname, "rt");
  if (in == NULL)
    return(1);

  in->warning = 1;
  in->kwdmode = sequent_mode;
  in->ignorecase = 1;
  Dbmatt = new dbmatt;  // allocating new dbmatt class
  if (d->matstr == yes)
    // material parameters are read as strings
    ret = Dbmatt->read(in);
  else
  {
    // material parameters are read with help of transmat procedures
    Mesprt = 0;
    ret = Dbmatt->readtm(in, Tm, d);
    Mesprt = mpb;
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
  @retval 0 - on success
  @retval 1 - if fails detecting of material XFILE section
  @retval 2 - if fails reading material parameters

  Created by TKo, 06.2014
*/
long input_materialst(XFILE *in, descript *d)
{
  long ret;
  long mpb = Mesprt;

  fprintf(stdout, "\n\nReading of material database . . .");
  if (xf_setsec(in, bsec_str[begsec_mater]))
    return(1);

  in->warning = 1;
  in->ignorecase = 1;
  in->kwdmode = sect_mode_seq;
  Dbmatt = new dbmatt;  // allocating new dbmatt class
  if (d->matstr == yes)
    // material parameters are read as strings
    ret = Dbmatt->read(in);
  else
  {
    // material parameters are read with help of transmat procedures
    Mesprt = 0;
    ret = Dbmatt->readtm(in, Tm, d);
    Mesprt = mpb;
  }
  if (ret)
    return(2);

  return(0);
}



/**
  Function reads cross-section database from file with name given by the fname

  @param fname - string with database file name
  @param d - pointer to structure with description of preprocessor setup 

  @retval 0 - on success
  @retval 1 - if fails opening file
  @retval 2 - if fails reading cross-section parameters

  Created by TKo, 09.2010
*/
long input_crst(char *fname, descript *d)
{
  long ret;
  long mpb = Mesprt;
  XFILE *in = NULL;

  fprintf(stdout, "\n\nReading of cross-section database from %s . . .", fname);
  in = xfopen(fname, "rt");
  if (in == NULL)
    return(1);

  in->warning = 1;
  in->ignorecase = 1;
  in->kwdmode = sequent_mode;
  Dbcrst = new dbcrst; // allocating new dbcrst class
  if (d->crsstr == yes)
    // cross section parameters are read as strings
    ret = Dbcrst->read(in);
  else
  {
    // cross section parameters are read with help of transcrsec procedures
    Mesprt = 0;
    ret = Dbcrst->readtc(in, Tc, d);
    Mesprt = mpb;
  }
  xfclose(in);
  if (ret)
    return(2);

  return(0);
}



/**
  Function reads cross-section database from corresponding XFILE  section 

  @param in - pointer to the opened XFILE
  @param d - pointer to structure with description of preprocessor setup 

  Returns:
   @retval 0 - on success
   @retval 1 - if fails opening file
   @retval 2 - if fails reading cross-section parameters

  Created by TKo, 06.2014
*/
long input_crst(XFILE *in, descript *d)
{
  long ret;
  long mpb = Mesprt;

  fprintf(stdout, "\n\nReading of cross-section database . . .");
  if (xf_setsec(in, bsec_str[begsec_crsec]))
    return(1);

  in->warning = 1;
  in->kwdmode = sect_mode_seq;
  in->ignorecase = 1;
  Dbcrst = new dbcrst; // allocating new dbcrst class
  if (d->crsstr == yes)
    // cross section parameters are read as strings
    ret = Dbcrst->read(in);
  else
  {
    // cross section parameters are read with help of transcrsec procedures
    Mesprt = 0;
    ret = Dbcrst->readtc(in, Tc, d);
    Mesprt = mpb;
  }
  if (ret)
    return(2);

  return(0);
}



/**
  Function reads nodal properties from file in
 
  @param in  - poinetr to opened input file with property description

  @retval 0  - on success
  @retval 1  - fails input of boundary conditions
  @retval 2  - fails input of common code numbers
  @retval 3  - fails input of time functions of dofs
  @retval 4  - fails input of cross-sections
  @retval 5 -  fails input of initial conditions
  @retval 6 -  fails input of nodal sources
  @retval 7 -  fails input of advection velocities

  Created 09.2010 by Tomas Koudelka koudelka@cml.fsv.cvut.cz
*/
long input_nodpropt(XFILE *in)
{
  long err;
  const enumstr nodsects[] = {{"begsec_nodvertpr",3}, {"begsec_nodedgpr",4}, 
                        {"begsec_nodsurfpr", 5}, {"begsec_nodvolpr",6}};
  long nsect = sizeof(nodsects)/sizeof(*nodsects);
  
  Nod_ccn = NULL;

  fprintf(stdout, "\n\nReading of nodal properties . . .");
  in->kwdmode = sequent_mode;
  // reading of boundary conditions at nodes 
  // (i.e. prescribed values of unknowns)
  err = input_nod_bocont(in, nodsects, nsect);
  if (err)
    return 1;
  // reading of coupled dofs at nodes
  err = input_nod_coupl_dofst(in, nodsects, nsect);
  if (err)
    return 2;
  // reading of time functions of dofs
  err = input_nod_dof_tfunct(in, nodsects, nsect);
  if (err)
    return 3;
  // reading of cross-sections at nodes
  err = input_nod_crsect(in, nodsects, nsect);   
  if (err)
    return 4;
  // reading of initial conditions
  err = input_nod_initcondt(in, nodsects, nsect);
  if (err)
    return 5;
  // reading of sources at nodes
  err = input_nod_sourcet(in, nodsects, nsect);
  if (err)
    return 6;
  // reading of advection velocities at nodes
  if (Tp->advect)
  {
    err = input_nod_advect_vel(in, nodsects, nsect);
    if (err)
      return 7;
  }
  return 0;
}



/**
  Function reads element properties from file in
 
  @param in  - poinetr to opened input file with property description

  @retval 0  - on success
  @retval 1  - fails input of element types
  @retval 2  - fails input of material
  @retval 3  - fails input of cross-sections
  @retval 4  - fails input of element sources
  @retval 5  - fails input of vertex boundary conditions
  @retval 6  - fails input of edge boundary conditions
  @retval 7  - fails input of surface boundary conditions
  @retval 8  - fails input of time functions for elements

  Created 09.2010 by Tomas Koudelka koudelka@cml.fsv.cvut.cz
*/
long input_elempropt(XFILE *in)
{
  long i, err;

  // all sections with element properties
  const enumstr elemsects[] = {{bsec_str[begsec_elvertpr].alias, bsec_str[begsec_elvertpr].id},
                               {bsec_str[begsec_eledgpr].alias, bsec_str[begsec_eledgpr].id},
                               {bsec_str[begsec_elsurfpr].alias, bsec_str[begsec_elsurfpr].id},
                               {bsec_str[begsec_elvolpr].alias, bsec_str[begsec_elvolpr].id}};
  long nsect = sizeof(elemsects)/sizeof(*elemsects);

  // sections with vertex properties for elements
  const enumstr vertsect[] = {{bsec_str[begsec_elvertpr].alias, bsec_str[begsec_elvertpr].id}};
  long nvertsect = sizeof(vertsect)/sizeof(*vertsect);

  // sections with edge properties for elements
  const enumstr edgesect[] = {{bsec_str[begsec_eledgpr].alias, bsec_str[begsec_eledgpr].id}};
  long nedgesect = sizeof(edgesect)/sizeof(*edgesect);

  // sections with surface properties for elements
  const enumstr surfsect[] = {{bsec_str[begsec_elsurfpr].alias, bsec_str[begsec_elsurfpr].id}};
  long nsurfsect = sizeof(surfsect)/sizeof(*surfsect);

  El_loadt = NULL;
  fprintf(stdout, "\n\nReading of element properties . . .");
  in->kwdmode = sequent_mode;
  // reading of element type
  err = input_elem_typet(in, elemsects, nsect);
  if (err)
    return 1;
  // reading of element material
  err = input_elem_matt(in, elemsects, nsect);
  if (err)
    return 2;
  // reading of element cross-section
  err = input_elem_crsect(in, elemsects, nsect);   
  if (err)
    return 3;
  // reading of sources on elements
  err = input_elem_sourcet(in, elemsects, nsect);
  if (err)
    return 4;
  // reading of vertex boundary conditions assigned to elements
  err = input_elem_vertbc(in, vertsect, nvertsect);
  if (err)
    return 5;
  // reading of edge boundary conditions assigned to elements
  err = input_elem_edgebc(in, edgesect, nedgesect);
  if (err)
    return 6;
  // reading of surface boundary conditions assigned to elements
  err = input_elem_surfbc(in, surfsect, nsurfsect);
  for (i=0; i<Tp->ntm; i++) // remove auxiliary arrays used just in reading element BC
  {
    delete [] El_loadtln[i];
    delete [] El_loadtcol[i];
  }
  delete [] El_loadtln;
  delete [] El_loadtcol;

  if (err)
    return 7;
  // reading of element switching time function
  err = input_elem_eltimefunct(in, elemsects, nsect);
  if (err)
    return 8;
  // reading of element integration domain descriptors
  err = input_elem_fluxint(in, surfsect, nsurfsect);
  if (err)
    return 9;

  return 0;
}



/**
  Function reads and assigns boundary conditions at each node.
  It scans sections nodvertpr, nodedgpr, nodsurfpr, nodvolpr for keyword "bocon".
  The record which assigns dofs to node looks like follows:
  "bocon" "propid" prop "lc_id" lcid "ini_cd" ini_val {"cond" cond_val | gfunct}

  The asssigning starts at section for volumes and finishes at
  section for vertices. If a node belongs to several entities and boundary conditions are
  prescribed for these entities, the assigned boundary conditions are merged and a message 
  is written into the log file. In case of multiple assigning of different boundary conditions 
  in same direction, error of merging is signalized.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  @retval 0 - on success  
  @retval 1 - assignment of different number of dofs at one node
  @retval 2 - invalid load case id or expression could not be parsed
  @retval 3 - two boundary conditons cannot be merged
  @retval 4 - growing transport problem type is required (dofs have to be controlled by time functions)
  @retval 5 - nodes with required property and entity type cannot be found
*/
long input_nod_bocont(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  long *line;
  long *col;
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long *setnodes;
  bocont *tbc;
  bocont *aux;
  long i, j, k, lcid, bcid;
  long prop, ndof, nnp;
  long errcode;

  // number of dofs for particualr nodes is the same and it equals to number of transported media
  ndof = Tp->ntm;
  Nod_bocon = new long*[Top->nn];
  memset(Nod_bocon, 0, sizeof(*Nod_bocon)*Top->nn);
  setnodes = new long[Top->nn];
  memset(setnodes, 0, sizeof(*setnodes)*Top->nn);
  line = new long[Top->nn]; 
  memset(line, 0, sizeof(*line)*Top->nn);
  col  = new long[Top->nn]; 
  memset(col, 0, sizeof(*col)*Top->nn);

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
    if (nkwd && ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin)))
    {
      sprintf(errmsg, "Problem type 'growing_np_problem|(_nonlin)' is required.\n"
                      " Nodal dofs have to be controlled by time functions.\n"
                      " Use 'nod_tfunc' preprocessor keyword for dof control.\n"
                      " Invalid keyword 'bocon' found at line %ld, column %ld, file %s",
                      in->line, in->col, in->fname);
      print_err(errmsg, __FILE__, __LINE__, __func__);
      delete [] setnodes;
      delete [] line;
      delete [] col;
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
      
      // investigation of ndof for nodes with property prop of entity ent
      nnp = Top->get_propent_nodes(prop, ent, setnodes);
      if (nnp == 0)
      {
        if (Check_unused_prop == 0)
        // no nodes were found with the given property and checking of unused properties is off
          continue;
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Boundary condition (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        delete [] setnodes;
        delete [] line;
        delete [] col;
        return 5;
      }
      // rest of record of boundary condition is read
      tbc = new bocont;
      errcode = tbc->read(in, ndof);
      if (errcode)
      {
        if (errcode < 3)
        {
          sprintf(errmsg, "Boundary condition or load case id is invalid\n"
                          " in file %s, line %ld, col %ld", in->fname, aline, acol);
        }
        else
          sprintf(errmsg, "Expression for boundary condition could not be parsed\n"
                          " in file %s, line %ld, col %ld", in->fname, aline, acol);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        delete [] setnodes;
        delete tbc;
        delete [] line;
        delete [] col;
        return 2;
      }

      in->kwdmode = bkwdmode;

      lcid = tbc->lcid;
      bcid = 0;
      // assigning of boundary condition
      for (k=0; k<Top->nn; k++)
      {        
        if (setnodes[k] < 0) // node has not required property
          continue;
        if (Nod_bocon[k])  // node has assigned BC from the previous steps
        {
          if (Nod_bocon[k][lcid])  // node has already assigned BC for the lcid
          {
            aux = (bocont *)(Nod_bclst.at(Nod_bocon[k][lcid]-1));
            if (tbc->compare(*aux))
            {
              sprintf(errmsg, "Conflict in merging of boundary condition at line %ld, col %ld\n"
                              "with previous boundary condition at line %ld, col %ld in file %s", aline, acol, line[k], col[k], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete [] setnodes;
              delete tbc;
              delete [] line;
              delete [] col;
              return 3;
            }
            else
              fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already had assigned boundary condition at line %ld, col %ld, file: %s\n",
                      aline, acol, in->fname, k+1, line[k], col[k], in->fname);
          }
          else // no BC in the lcid
          {
            if (bcid == 0)
            {
              Nod_bclst.append(tbc);
              bcid = Nod_bclst.count();
            }
            Nod_bocon[k][lcid] = bcid;
            // backup of line and column for eventual log message
            line[k]      = aline;
            col[k]       = acol;
          }
        }
        else // no BC assigned at node
        {
          Nod_bocon[k] = new long[ndof];
          memset(Nod_bocon[k], 0, sizeof(Nod_bocon[k])*ndof);
          if (bcid == 0)
          {
            Nod_bclst.append(tbc);
            bcid = Nod_bclst.count();
          }
          Nod_bocon[k][lcid] = bcid;
          // backup of line and column for eventual log message
          line[k]      = aline;
          col[k]       = acol;
        }
      }  // end of loop for nodes
    }    // end of nkwd loop 
  }      // end of section loop
  delete [] setnodes;
  delete [] line;
  delete [] col;
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

  @retval 0 - on success  
  @retval 1 - assignment of different number of dofs at one node
  @retval 2 - invalid number of coupled dofs (out of <1;ndof>)  
  @retval 3 - invalid coupled dof number (out of <1;ndof>)
  @retval 4 - growing transport problem type is required (dofs have to be controlled by time functions)
  @retval 5 - nodes with required property and entity type cannot be found

  Created by TKo, 09.2010
*/
long input_nod_coupl_dofst(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  long *line;
  long *col;
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long *setnodes;
  long i, j, k, l;
  long prop, ndof, lcid, nnp, nbc;
  long ccn, accn;

  Nod_ccn = NULL;
  setnodes = new long[Top->nn];
  memset(setnodes, 0, sizeof(*setnodes)*Top->nn);
  line = new long[Top->nn]; 
  memset(line, 0, sizeof(*line)*Top->nn);
  col  = new long[Top->nn]; 
  memset(col, 0, sizeof(*col)*Top->nn);
  ccn = accn = 1;
  // number of dofs for particualr nodes is the same and it equals to number of transported media
  ndof = Tp->ntm;

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
    if (nkwd && ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin)))
    {
      sprintf(errmsg, "Problem type 'growing_np_problem|(_nonlin)' is required.\n"
                      " Nodal dofs have to be controlled by time functions.\n"
                      " Use 'nod_tfunc' preprocessor keyword for dof control.\n"
                      " Invalid keyword 'dof_coupl' found at line %ld, column %ld, file %s",
                      in->line, in->col, in->fname);
      print_err(errmsg, __FILE__, __LINE__, __func__);
      delete [] setnodes;
      delete [] line;
      delete [] col;
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
      // investigation of ndof for nodes with property prop of entity ent
      nnp = Top->get_propent_nodes(prop, ent, setnodes);
      if (nnp == 0)
      {
        if (Check_unused_prop == 0)
        // no nodes were found with the given property and checking of unused properties is off
          continue;
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Coupled dofs (line %ld, column %ld, file %s) cannot be assigned correctly.",
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        delete [] setnodes;
        delete [] line;
        delete [] col;
        return 5;
      }
      xfscanf(in, "%k%ld", "nbc", &nbc);
      if ((nbc < 1) || (nbc > ndof))
      {
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s\n"
                        " have number of coupled dofs out of range <1,%ld>"
                        " (line %ld, column %ld, file %s)", 
                prop, entitypstr[ent-1].alias, ndof, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        delete [] setnodes;
        delete [] line;
        delete [] col;
        return 2;
      }
      for (k = 0; k < nbc; k++)
      {
        xfscanf(in, "%k%ld", "lc_id", &lcid);
        if ((lcid < 1) || (lcid > ndof))
        {
          sprintf(errmsg, "Nodes with property %ld belonging to entity %s\n"
                          " have coupled dof number out of range <1,%ld>" 
                          " (line %ld, column %ld, file %s)", 
                  prop, entitypstr[ent-1].alias, ndof, aline, acol, in->fname);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          delete [] setnodes;
          delete [] line;
          delete [] col;
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
          if (Nod_ccn[l][lcid-1] != 0)
          {
            accn = Nod_ccn[l][lcid-1];
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
          Nod_ccn[l][lcid-1] = accn;
        } // end of loop for nodes
      }   // end of loop for nbc
      in->kwdmode = bkwdmode;
    }  // end of nkwd loop 
  }    // end of section loop 
  delete [] setnodes;
  delete [] line;
  delete [] col;
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

  @retval 0 - on success  
  @retval 1 - cannot determine unique number of dofs for the given property and entity type
  @retval 2 - number of prescribed time functions is out of range <1,ntm>
  @retval 3 - dof id is out of range <1,ntm>
  @retval 4 - nodes with required property and entity type cannot be found
  @retval 5 - invalid time function id


  Created by TKo, 09.2010
*/
long input_nod_dof_tfunct(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  long *line;
  long *col;
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long *setnodes;
  long tf_id;
  long i, j, k, l;
  long prop, ndof, nbc, nnp, lcid;

  if (Nod_ccn == NULL)
  {
    Nod_ccn = new long*[Top->nn];
    memset(Nod_ccn, 0, sizeof(*Nod_ccn)*Top->nn);
  }
  setnodes = new long[Top->nn];
  memset(setnodes, 0, sizeof(*setnodes)*Top->nn);
  line = new long[Top->nn]; 
  memset(line, 0, sizeof(*line)*Top->nn);
  col  = new long[Top->nn]; 
  memset(col, 0, sizeof(*col)*Top->nn);
  // number of dofs for particualr nodes is the same and it equals to number of transported media
  ndof = Tp->ntm;

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
      // investigation of ndof for nodes with property prop of entity ent
      nnp = Top->get_propent_nodes(prop, ent, setnodes);
      if (nnp == 0)
      {
        if (Check_unused_prop == 0)
        // no nodes were found with the given property and checking of unused properties is off
          continue;
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Nodal time functions (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        delete [] setnodes;
        delete [] line;
        delete [] col;
        return 4;
      }
      xfscanf(in, "%k%ld", "nbc", &nbc);
      if ((nbc < 1) || (nbc > ndof))
      {
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s\n"
                        " have number of dof out of range <1,%ld>"
                        " (line %ld, column %ld, file %s)", 
                prop, entitypstr[ent-1].alias, ndof, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        delete [] setnodes;
        delete [] line;
        delete [] col;
        return 2;
      }
      for (k = 0; k < nbc; k++)
      {
        xfscanf(in, "%k%ld", "lc_id", &lcid);
        if ((lcid < 1) || (lcid > ndof))
        {
          sprintf(errmsg, "Nodes with property %ld belonging to entity %s\n"
                          " have time function dof number out of range <1,%ld>" 
                          " (line %ld, column %ld, file %s)", 
                  prop, entitypstr[ent-1].alias, ndof, aline, acol, in->fname);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          delete [] setnodes;
          delete [] line;
          delete [] col;
          return 3;
        }
        xfscanf(in, "%k%ld", "tfunc_id", &tf_id);
        if (tf_id < 1)
        {
          sprintf(errmsg, "Invalid time function index");
          print_err(errmsg, __FILE__, __LINE__, __func__);
          delete [] line;
          delete [] col;
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
            Nod_ccn[l][lcid-1] = tf_id;
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
            Nod_ccn[l][lcid-1] = tf_id;
          }
        }  // end of loop for nodes
      }    // end of loop for nbc
      in->kwdmode = bkwdmode;
    } // end of nkwd loop
  }   // end of section loop
  delete [] setnodes;
  delete [] line;
  delete [] col;
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

  @retval 0 - on success  
  @retval 1 - cannot find required cross-section in database
  @retval 2 - different cross-section type has been assigned at one node
  @retval 3 - nodes with required property and entity type cannot be found

  Created by TKo, 09.2010
*/
long input_nod_crsect(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  long *line;
  long *col;
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long *setnodes;
  crsectypet cst;
  long       csti;
  long i, j, k, ic;
  long prop;
  long nnp;

  Nod_cst = new crsectypet[Top->nn];
  memset(Nod_cst, 0, sizeof(*Nod_cst)*Top->nn);
  Nod_csti = new long[Top->nn];
  memset(Nod_csti, 0, sizeof(*Nod_csti)*Top->nn);
  Nod_cstdbi = new long[Top->nn];
  memset(Nod_cstdbi, 0, sizeof(*Nod_cstdbi)*Top->nn);
  setnodes = new long[Top->nn];
  memset(setnodes, 0, sizeof(*setnodes)*Top->nn);
  line = new long[Top->nn]; 
  memset(line, 0, sizeof(*line)*Top->nn);
  col  = new long[Top->nn]; 
  memset(col, 0, sizeof(*col)*Top->nn);

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
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "nod_crsec", "propid", &prop);
      xfscanf(in, "%k%m%k%ld", "type", &crsectypet_kwdset, &cst, "type_id", &csti);
      // searching of prescribed cross-section type and index in the cross-section database
      ic = Dbcrst->search_crs(cst, csti-1);
      if (ic < 0)
      {
        sprintf(errmsg, "Line %ld, col %ld, file %s: cannot find required cross-section in the database", in->line, in->col, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        delete [] setnodes;
        delete [] line;
        delete [] col;
        return 1;
      }
      else
        Dbcrst->mark_used(ic, csti-1);
          
      in->kwdmode = bkwdmode;

      nnp = Top->get_propent_nodes(prop, ent, setnodes);
      if (nnp == 0)
      {
        if (Check_unused_prop == 0)
        // no nodes were found with the given property and checking of unused properties is off
          continue;
        // no nodes were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Nodal cross-section type (line %ld, column %ld, file %s) cannot be assigned correctly.",
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        delete [] line;
        delete [] col;
        return 3;
      }
      
      // assigning of a time function index to the given dof
      for (k=0; k<Top->nn; k++)
      {        
        if (setnodes[k] < 0) // node has not required property
          continue;
        if ((Nod_cst[k] != nocrosssectiont) && ((Nod_cst[k] != cst) || (Nod_csti[k] != csti)))        
        {  // node has had assigned different cross-section types
          sprintf(errmsg, "Line %ld, col %ld, file %s: Node %ld has already had assigned different"
                          "cross-section type at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                          k+1, line[k], col[k], in->fname);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          delete [] setnodes;
          delete [] line;
          delete [] col;
          return 2;
        }
        Nod_cst[k]  = cst;
        Nod_csti[k] = csti;
        Nod_cstdbi[k] = ic;
        // backup of line and column for eventual log message
        line[k] = aline;
        col[k]  = acol;
      }
    }  // end of nkwd loop
  }    // end of section loop
  delete [] setnodes;
  delete [] line;
  delete [] col;
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

  @retval 0 - on success  
  @retval 1 - cannot read initial conditions from file
  @retval 2 - assignment of different types of initial conditons at one node
  @retval 3 - nodes with required property and entity type cannot be found

  Created by TKo, 09.2010
*/
long input_nod_initcondt(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char emsg[1001];
  long *line;
  long *col;
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long *setnodes;
  long check;
  long i, j, k, l, prop;
  long nini;
  long ndof, nnp;
  double *con;
  prepcondtype cndt;


  Nod_inicd = new double*[Top->nn];
  memset(Nod_inicd, 0, sizeof(*Nod_inicd)*Top->nn);
  // number of dofs for particualr nodes is the same and it equals to number of transported media
  ndof = Tp->ntm;
  if (Inicdf)
  {
    if (read_inicd_file(Inicdf, ndof))
      return 1;
  }

  setnodes = new long[Top->nn];
  memset(setnodes, 0, sizeof(*setnodes)*Top->nn);
  line = new long[Top->nn]; 
  memset(line, 0, sizeof(*line)*Top->nn);
  col  = new long[Top->nn]; 
  memset(col, 0, sizeof(*col)*Top->nn);
  con = new double[ndof];

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of nodal initial conditions");

  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "nod_tdload"
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
      xfscanf(in, "%k%k%ld", "nod_inicond", "propid", &prop);
      xfscanf(in, "%m", &prepcondtype_kwdset, &cndt);
      switch (cndt)
      {
        case cond:
          for (l=0; l<ndof; l++)
            xfscanf(in, "%le", con+l);
          break;
        case file:
          break;
        default:
          print_err("Unknown type of initial condition setup.", __FILE__, __LINE__, __func__);
      }

      in->kwdmode = bkwdmode;

      // assigning of initial conditions to the given nodes
      nnp = Top->get_propent_nodes(prop, ent, setnodes);
      if (nnp == 0)
      {
        if (Check_unused_prop == 0)
          // no nodes were found with the given property and checking of unused properties is off
          continue;
        else
        {
          // no nodes were found with the given property and checking of unused properties is on
          sprintf(emsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Nodal intial conditions (line %ld, column %ld, file %s) cannot be assigned correctly.",
                  prop, entitypstr[ent-1].alias, aline, acol, in->fname);
          print_err(emsg, __FILE__, __LINE__, __func__);
          delete [] con;
          delete [] line;
          delete [] col;
          return 3;
        }
      }        
      for (k=0; k<Top->nn; k++)
      {        
        if (setnodes[k] < 0) // node has not required property
          continue;

        if (Nod_hang)
        {
          if (Nod_hang[k]) // hanging node has defined initial conditions in the master nodes
          {
            fprintf(Log, "Line %ld, col %ld, file %s: Node %ld is a hangining node, initial conditions are defined by its master nodes, prescribed initial values will be ignored at this node\n", aline, acol, in->fname, k+1);
            continue;
          }
        }

        if (Nod_inicd[k] == NULL)
        { // node has not assigned load yet
          Nod_inicd[k] = new double[ndof];
          for(l=0; l<ndof; l++)
            Nod_inicd[k][l] = con[l];
          // backup of line and column for eventual log message
          line[k] = aline;
          col[k]  = acol;
        }
        else
        { // node has already had assigned initial condition
          check = 0;
          if (cndt == cond)
          {
            for(l=0; l<ndof; l++)
            {
              if (con[l] != Nod_inicd[k][l])
                check ++;
            }
          }
          if ((cndt == cond) && (line[k] != 0) && check) // two different conditions were specified in both nodal section and initial condition file
          { // error occured while merging
            sprintf(emsg, "Line %ld, col %ld, file %s: Node %ld has already had assigned different "
                          "initial condition at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
            print_err(emsg, __FILE__, __LINE__, __func__);
            delete [] setnodes;
            delete [] line;
            delete [] col;
            return 2;
          }
          else
          {
            if (line[k])// the same initial conditions were specified at two different lines
              fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already had assigned initial condition at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
            else 
            {           
              if (check)// initial condition was specified in both nodal section and initial condition file
              {
                fprintf(Log, "Line %ld, col %ld, file %s: Initial condition for node %ld specified in file '%s'\n"
                             "                            has been replaced by new one from the above mentioned line\n", 
                             aline, acol, in->fname, k+1, Inicdf);
              }
              else// the same initial conditions were specified at the given line and initial condition file
                fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already had assigned the same initial condition in file: %s\n", aline, acol, in->fname, k+1, Inicdf);
              line[k] = aline;
              col[k] = acol;
            }
          }
        }
      }
    }  // end of nkwd loop
  }    // end of section loop

  // check for completeness of intial condition 
  // number of initial conditions has to be either Top->nn or 0
  nini = 0;
  for (j=0; j<Top->nn; j++)
  {
    if (Nod_inicd[j])
      nini++;
  }
  if ((nini > 0) && (nini < Top->nn-Numhn))
  {
    sprintf(emsg, "number of initial conditions is less then number of nodes");
    print_err(emsg, __FILE__, __LINE__, __func__);
    delete [] con;
    delete [] setnodes;
    delete [] line;
    delete [] col;
    return 3;
  }
  delete [] con;
  delete [] setnodes;
  delete [] line;
  delete [] col;
  return 0;
}



/**
  The function reads initial conditions at nodes stored in the file given by the cndfname.
  The file must be according to the format of the text output file from the TRFEL/PARTRFEL.
  Values of all unknown quanitities at nodes must be included.
  The function stores the initial values directly in the Nod_inicd variable
  which is either allocated and filled with the read values or it is set to the NULL.

  @param cndfname - file name with initial condition values at nodes (created by TRFEL/PARTRF)
  @param ndof - number of DOFs, i.e. number of transported quantities in the problem.
  
  @retval 0 - on success
  @retval 1 - initial condition file could not be opened
  @retval 2 - node number out of range <1, nn>

  Created by Tomas Koudelka, 03.2012.
*/
long read_inicd_file(char *cndfname, long ndof)
{
  long i, j, nid;
  XFILE *inif;
  long nn = Top->nn;
  char label[10];

  inif = xfopen(cndfname, "rt");
  if (inif == NULL)
  {
    print_err("Cannot open initial condition file '%s'", 
              __FILE__, __LINE__, __func__, cndfname);
    return (1);
  }
  inif->warning = 1;
  //inif->kwdmode = ignore_kwd;
  inif->kwdmode = sequent_mode;
  inif->ignorecase = 1;
  Nod_inicd = new double*[nn];
  memset(Nod_inicd, 0, sizeof(*Nod_inicd)*nn);
  for (i=0; i<nn; i++)
  {
    Nod_inicd[i] = new double[ndof];
    memset(Nod_inicd[i], 0, sizeof(*Nod_inicd[i])*ndof);
  }
  // skipping header of the text output file
  xfscanf(inif,"%* a%* a%* a%* a%* a");
  // skipping problem label
  xfscanf(inif,"%* a");
  // skipping step label
  xfscanf(inif,"%* a%* a%* a");
  // skipping "nodal unknowns" label
  xfscanf(inif,"%* a");
  for (i=0; i<nn; i++)
  {
    xfscanf(inif, " %k%ld", "Node", &nid);
    if ((nid < 1) && (nid > nn))
    {
      print_err("Node number is out of range <1,%ld>.\n"
                "Initial condition file '%s', line=%ld, col=%ld", 
                __FILE__, __LINE__, __func__, nn, cndfname, inif->line, inif->col); 
      return(1);
    }
    for (j=0; j<ndof; j++)
    {
      snprintf(label, 10, "r_%ld=", j+1);
      xfscanf(inif, " %k%le", label, Nod_inicd[i]+j);
    }
  }
  return(0);
}
 


/**
  Function reads and assigns nodal sources of given quantity at each node.
  It scans sections nodvertpr, nodedgpr, nodsurfpr, nodvolpr for keyword "nod_src".
  The record which assigns sources to node looks like follows:
  "nod_src" "propid" prop "lc_id" lcid {"src_val" v}|{"src_type" srt {gfunct|hydrh|cemhyd}}
                                      |___________| |____________________________________|
                                            |                          |
                                        stat_prob                all_other_prob

  The record contains call of sourcet::read function over the reading of lcid,  For more details
  see sourcet.cpp.

  The asssigning starts at section for volumes and finishes at
  section for vertices. If a node belongs to several entities and a source of quantity is
  prescribed for these entities, the sources of quantity are compared and if they differ, the error 
  of assignment is signalized.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  @retval 0 - on success  
  @retval 1 - invalid load case id
  @retval 2 - two boundary conditons cannot be merged
  @retval 3 - nodes with required property and entity type cannot be found
  
  Created by TKo, 09.2010
*/
long input_nod_sourcet(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->nn), col(Top->nn);
  ivector setnodes(Top->nn);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  sourcet *tsrc, *aux;
  long lcid;
  long i, j, k;
  long prop, ndof, nnp;
  long srcid;

  Nod_sourcet = new long*[Top->nn];
  memset(Nod_sourcet, 0, sizeof(*Nod_sourcet)*Top->nn);
  // number of dofs for particualr nodes is the same and it equals to number of transported media
  ndof = Tp->ntm;

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of nodal sources");
  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "nod_src"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "nod_src", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "nod_src", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "nod_src", "propid", &prop);
      // investigation for nodes with property prop of entity ent
      nnp = Top->get_propent_nodes(prop, ent, setnodes);
      if (nnp == 0)
      {
        if (Check_unused_prop == 0)
        // no nodes were found with the given property and checking of unused properties is off
          continue;
        sprintf(errmsg, "Nodes with property %ld belonging to entity %s have not been found.\n"
                        " Nodal source (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 3;
      }
      // read load case id
      xfscanf(in, "%k%ld", "lc_id", &lcid);
      if ((lcid < 1) || (lcid > ndof))
      {
        sprintf(errmsg, "Load case index is out of range <1, %ld>\n"
                        " in file %s, line %ld, col %ld", ndof, in->fname, aline, acol);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 1;
      }
      lcid--;

      tsrc = new sourcet;
      // rest of record of nodal source is read
      tsrc->read(in);
      srcid = 0;

      in->kwdmode = bkwdmode;

      // assigning of nodal sources
      for (k=0; k<Top->nn; k++)
      {        
        if (setnodes[k] < 0) // node has not required property
          continue;
        if (Nod_sourcet[k])
        {
          if (Nod_sourcet[k][lcid])
          {
            // node has already had assigned source of quantity 
            fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already had assigned nodal source at line %ld, col %ld, file: %s\n",
                    aline, acol, in->fname, k+1, line[k], col[k], in->fname);
            aux = (sourcet *)(Nod_srclst.at(Nod_sourcet[k][lcid]-1));
            if (tsrc->compare(*aux))
            {
              sprintf(errmsg, "Conflict in assignment of nodal source at line %ld, col %ld\n"
                              "with previous nodal source at line %ld, col %ld in file %s", aline, acol, line[k], col[k], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              delete tsrc;
              return 2;
            }
          }
        }
        // node has had assigned no nodal sources of any quantity
        if (Nod_sourcet[k] == NULL)
        {
          Nod_sourcet[k] = new long[ndof];
          memset(Nod_sourcet[k], 0, sizeof(*Nod_sourcet[k])*ndof);
        }
        // store source in the list of source, it will be used for memory delete
        // (Nod_sourcet contains pointers to one object of tsrc)
        if (srcid == 0)
        {
          Nod_srclst.append(tsrc);
          srcid = Nod_srclst.count();
        }
        Nod_sourcet[k][lcid] = srcid;
        // backup of line and column for eventual log message
        line[k]      = aline;
        col[k]       = acol;
      }  // end of loop for nodes
    }    // end of nkwd loop 
  }      // end of section loop
  return 0;
}






/**
  Function reads and assigns nodal velocities of medium for advection contribution of given quantity at each node.
  It scans sections nodvertpr, nodedgpr, nodsurfpr, nodvolpr for keyword "nod_advect_vel".
  The record which assigns velocities to node looks like follows:
  "nod_advect_vel" "propid" prop "lc_id" lcid "ncomp" nc "velocity_comp" {gfunct}xnc

  The asssigning starts at section for volumes and finishes at  section for vertices. If a node belongs to 
  several entities and a velocities are prescribed for these entities, the velocity values are compared 
  and if they differ, the error of assignment is signalized.

  @param in  - pointer to opened input file with property description
  @nodsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  @retval 0 - on success  
  @retval 1 - invalid load case id
  @retval 2 - invalid number of velocity components
  @retval 3 - wrong definition of advection velocity components
  @retval 4 - two velocities cannot be merged
  @retval 5 - nodes with required property and entity type cannot be found
  @retval 6 - advection velocity vector has not been prescribed at all nodes
  
  Created by TKo, 05.2016
*/
long input_nod_advect_vel(XFILE *in, const enumstr nodsects[], long nsect)
{
  char *aptr;
  long nkwd;
  ivector line(Top->nn), col(Top->nn);
  ivector setnodes(Top->nn);
  long aline, acol, pline;
  gentity ent;
  kwd_handling bkwdmode;
  advectvel *v = NULL;
  double auxv;
  long nc = 0;
  long err;
  long i, j, k, l;
  long prop, ndof, nnp;

  Nod_advect = new double*[Top->nn];
  memset(Nod_advect, 0, sizeof(*Nod_advect)*Top->nn);
  Nod_advect_nc = new long[Top->nn];
  memset(Nod_advect_nc, 0, sizeof(*Nod_advect_nc)*Top->nn);
  // number of dofs for particualr nodes is the same and it equals to number of transported media
  ndof = Tp->ntm;
  aline = 0;

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;
  fprintf(stdout, "\n reading of nodal advection velocities");
  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, nodsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "nod_advect_vel"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "nod_advect_vel", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "nod_advect_vel", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      pline = aline;
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "nod_advect_vel", "propid", &prop);
      // investigation for nodes with property prop of entity ent
      nnp = Top->get_propent_nodes(prop, ent, setnodes);
      if (nnp == 0)
      {
        if (Check_unused_prop == 0)
        // no nodes were found with the given property and checking of unused properties is off
          continue;
        print_err("Nodes with property %ld belonging to entity %s have not been found.\n"
                  " Nodal source (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                  __FILE__, __LINE__, __func__, prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        return 5;
      }

      // read advection velocity record
      if (v)   delete v;
      else nc = -1;
      v = new advectvel;
      err = v->read(in, ndof);
      if (nc < 0)
        nc = v->ncomp;
      else
      {
        if (nc != v->ncomp)
        {
          print_err("Different number of components of advection velocity vector is being defined.\n"
                    "ncomp=%ld on line %ld but ncomp=%ld has already been defined at line %ld of file %s",
                    __FILE__, __LINE__, __func__, v->ncomp, aline, nc, pline, in->fname);
          return 2;
        }
      }
      if (err)
        return err;

      in->kwdmode = bkwdmode;

      // assigning of nodal sources
      for (k=0; k<Top->nn; k++)
      {        
        if (setnodes[k] < 0) // node has not required property
          continue;
        if (Nod_advect[k])
        {         
          if (Nod_advect_nc[k] != v->ncomp)
          {
            print_err("Conflict in assignment of nodal advection velocities (ncomp=%ld) at line %ld, col %ld\n"
                      "with previous nodal advection velocities (ncomp=%ld) at line %ld, col %ld in file %s",
                      __FILE__, __LINE__, __func__, v->ncomp, aline, acol, Nod_advect_nc[k], line[k], col[k], in->fname);
            return 4;
          }
          if (Nod_advect[k])
          {
            // node has already had assigned source of quantity 
            fprintf(Log, "Line %ld, col %ld, file %s: Node %ld has already had assigned nodal advection velocities at line %ld, col %ld, file: %s\n",
                    aline, acol, in->fname, k+1, line[k], col[k], in->fname);
            // check particular components that they has got assigned the same values
            for(l=0; l<v->ncomp; l++) 
            {            
              // compute l-th velocity component value
              err = v->getval(Top->nodes[k], l, auxv);
              if (err) 
              {  // wrong component id passed to getval or wrong function type defined in v
                return 3;
              }
              if (auxv == Nod_advect[k][l]) // components are the same, check next component
                continue;
              else
              { // component has got assigned different value originally => error
                print_err("Conflict in assignment of nodal advection velocity at line %ld, col %ld\n"
                          "with previous nodal advection velocities at line %ld, col %ld in file %s at node %ld.\n"
                          "New value of %ld. velocity component %le != %le assigned originally.",
                          __FILE__, __LINE__, __func__, aline, acol, line[k], col[k], in->fname, k+1, l+1, auxv, Nod_advect[k][l]);
                return 4;
              }
            }
          }
        }
        else
        {
          Nod_advect_nc[k] = v->ncomp;
          Nod_advect[k] = new double[v->ncomp];
          for(l=0; l<v->ncomp; l++) 
          {            
            // compute l-th velocity component value
            err = v->getval(Top->nodes[k], l, Nod_advect[k][l]);
            if (err) 
            {  // wrong component id passed to getval or wrong function type defined in v
              return 3;
            }
          }
        }
      }  // end of loop for nodes
    }    // end of nkwd loop 
  }      // end of section loop

  // check whether all nodes have assigned an advection velocity vector
  for(i=0; i<Top->nn; i++)
  {    
    if (Nod_advect[i] == NULL)
    {
      print_err("Node %ld has not assigned advection velocity vector", __FILE__, __LINE__, __func__, i+1);
      return 6;
    }
  }
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

  @retval 0 - on success  
  @retval 1 - different types have been assigned to one element
  @retval 2 - no type has been assigned to elements
  @retval 3 - property numbers of the required entity type have not been read on elements
  @retval 4 - property number of the required entity type have not been found on elements

  Created by TKo, 09.2010
*/
long input_elem_typet(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->ne), col(Top->ne);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k;
  elemtypet tmp_type;
  long prop, prop_used;
  long sres;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of element types");
  El_type = new elemtypet[Top->ne];
  memset(El_type, 0, sizeof(*El_type)*Top->ne);

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
      xfscanf(in, "%k%k%ld%m", "el_type", "propid", &prop, &elemtypet_kwdset, &tmp_type);
      
      in->kwdmode = bkwdmode;

      prop_used = 0;
      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, Top->nodes, k, Top->edges, Top->surfaces, entid, nentid);
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
            fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already had assigned type at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
          El_type[k] = tmp_type;
          prop_used = 1;
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
                        " Element type (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 4;
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
      return 2;
    }
  }
  return 0;
}



/**
  The function assigns material type to the elements. It scans 
  sections eledgpr, elsurfpr, elvolpr for keyword "el_type".
  The record which assigns material type to the element looks like follows:
  "el_mat" "propid" prop {"type" material_type "type_id" material_id}xnmat
  where prop is positive integer number and nmat=ndof*ndof.
  
  The asssigning starts at section for volumes and finishes at
  section for edges. If an element  belongs to several entities and different 
  material types are prescribed for these entities, error of assignment is signalized.
  Finally, the function checks whether all elements have assigned a material type.
 
  @param in  - pointer to opened input file with property description
  @elemsects - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  @retval 0 - on success  
  @retval 1 - different material types have been assigned to one element
  @retval 2 - no type has been assigned to elements
  @retval 3 - property numbers of the required entity type have not been read on elements
  @retval 4 - property number of the required entity type have not been found on elements
 
  Created by TKo, 09.2010
*/
long input_elem_matt(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->ne), col(Top->ne);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k, l, err;
  mattypet *tmp_type;
  long *tmp_id, *ic;
  long prop, ndof, nmat;
  long sres;
  long prop_used;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of element material types");
  El_nmat = new long[Top->ne];
  memset(El_nmat, 0, sizeof(*El_nmat)*Top->ne);
  El_mattype = new mattypet*[Top->ne];
  memset(El_mattype, 0, sizeof(*El_mattype)*Top->ne);
  El_matid = new long*[Top->ne];
  memset(El_matid, 0, sizeof(*El_matid)*Top->ne);
  El_matdbi = new long*[Top->ne];
  memset(El_matdbi, 0, sizeof(*El_matdbi)*Top->ne);
  // number of dofs for particualr nodes is the same and it equals to number of transported media
  ndof = Tp->ntm;

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
      xfscanf(in, "%k%k%ld", "el_mat", "propid", &prop);
      nmat = ndof*ndof;
      tmp_type = new mattypet[nmat];
      memset(tmp_type, 0, sizeof(*tmp_type)*nmat);
      tmp_id   = new long[nmat];
      memset(tmp_id, 0, sizeof(*tmp_id)*nmat);
      ic       = new long[nmat];
      memset(ic, 0, sizeof(*ic)*nmat);

      for (k=0; k<nmat; k++)
      {
        xfscanf(in, "%k%m%k%ld", "type", &mattypet_kwdset, tmp_type+k, "type_id", tmp_id+k);
        // searching of prescribed material type and index in the material database
        ic[k] = Dbmatt->search_mat(tmp_type[k], tmp_id[k]-1);
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
          Dbmatt->mark_used(ic[k], tmp_id[k]-1);
      }
      
      in->kwdmode = bkwdmode;

      prop_used = 0;

      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, Top->nodes, k, Top->edges, Top->surfaces, entid, nentid);
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
              fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already had assigned material type at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
          }
          else
          {
            El_mattype[k] = new mattypet[nmat];
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
  "el_crsec" "propid" prop "crsec_type" crst "crsec_id" crsid
  where prop is positive integer number, crst is the type of cross section according to aliast and 
  propid is positive integer number.
  
  The asssigning starts at section for volumes and finishes at
  section for edges. If an element  belongs to several entities and different 
  cross-section types are prescribed for these entities, error of assignment is signalized.
 
  @param in  - pointer to opened input file with property description
  @elemsects - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  @retval 0 - on success 
  @retval 1 - cannot find required cross-section in database 
  @retval 2 - different cross-section types have been assigned to one element
  @retval 3 - property numbers of the required entity type have not been read on elements
  @retval 4 - property number of the required entity type have not been found on elements

  Created by TKo, 09.2010
*/
long input_elem_crsect(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->ne), col(Top->ne);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k, ic;
  crsectypet tmp_type;
  long prop, tmp_id;
  long sres;
  long prop_used;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of element cross-section types");
  El_cst = new crsectypet[Top->ne];
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
      xfscanf(in, "%k%k%ld%k%m%k%ld", "el_crsec", "propid", &prop, "type", &crsectypet_kwdset, &tmp_type, "type_id", &tmp_id);
      // searching of prescribed cross-section type and index in the cross-section database
      ic = Dbcrst->search_crs(tmp_type, tmp_id-1);
      if (ic < 0)
      {
        sprintf(errmsg, "Line %ld, col %ld, file %s: cannot find required cross-section in the database", in->line, in->col, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 1;
      }
      else
        Dbcrst->mark_used(ic, tmp_id-1);

      in->kwdmode = bkwdmode;

      prop_used = 0;
      // assign given type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, Top->nodes, k, Top->edges, Top->surfaces, entid, nentid);
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
            fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already had assigned cross-section type at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[k], col[k], in->fname);
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
  Function reads and assigns sources of given quantity on elements.
  It scans sections eledgpr, elsurfpr, elvolpr for keyword "el_src".
  The record which assigns sources to element looks like follows:
  "el_src" "propid" prop "lc_id" lcid {"src_val" v}|{"src_type" srt {gfunct|hydrh|cemhyd}}
                                      |___________| |____________________________________|
                                            |                          |
                                        stat_prob                all_other_prob

  The record contains call of sourcet::read function over the reading of lcid,  For more details
  see sourcet.cpp.

  The asssigning starts at section for volumes and finishes at
  section for vertices. If an element  belongs to several entities and a source of quantity is
  prescribed for these entities, the sources of quantity are compared and if they differ, the error 
  of assignment is signalized.

  @param in  - pointer to opened input file with property description
  @elemsects  - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array nodsects

  @retval 0 - on success  
  @retval 1 - invalid load case id
  @retval 2 - two boundary conditons cannot be merged
  @retval 3 - growing transport problem type is required (dofs have to be controlled by time functions)
  @retval 4 - nodes with required property and entity type cannot be found

  Created by TKo, 09.2010
*/
long input_elem_sourcet(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  ivector line(Top->ne), col(Top->ne);
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  sourcet *tsrc, *aux;
  long lcid;
  long i, j, k;
  long prop, ndof;
  long prop_used;
  long sres, srcid;
  long *entid;
  long nentid;

  El_sourcet = new long*[Top->ne];
  memset(El_sourcet, 0, sizeof(*El_sourcet)*Top->ne);
  // number of dofs for particular nodes is the same and it equals to number of transported media
  ndof = Tp->ntm;

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  fprintf(stdout, "\n reading of sources on elements");
  for (i=nsect-1, ent=gregion; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "el_src"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "el_src", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "el_src", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld", "el_src", "propid", &prop);
      // read load case id
      xfscanf(in, "%k%ld", "lc_id", &lcid);
      if ((lcid < 1) || (lcid > ndof))
      {
        sprintf(errmsg, "Load case index is out of range <1, %ld>\n"
                        " in file %s, line %ld, col %ld", ndof, in->fname, aline, acol);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 1;
      }
      lcid--;

      tsrc = new sourcet;
      // rest of record of source on element is read
      tsrc->read(in);

      in->kwdmode = bkwdmode;

      prop_used = 0;
      srcid = 0;

      // assigning of sources of quantity to elements
      for (k=0; k<Top->ne; k++)
      {        
        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, Top->nodes, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres == 0) // element has not required property
          continue;
        if (El_sourcet[k])
        {
          if (El_sourcet[k][lcid])
          {
            // element has already had assigned source of quantity
            fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already had assigned source at line %ld, col %ld, file: %s\n",
                    aline, acol, in->fname, k+1, line[k], col[k], in->fname);
            aux = (sourcet *)(El_srclst.at(El_sourcet[k][lcid]-1));
            if (tsrc->compare(*aux))
            {
              sprintf(errmsg, "Conflict in assignment of source on element at line %ld, col %ld\n"
                              "with previous source on element at line %ld, col %ld in file %s", aline, acol, line[k], col[k], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__); 
              return 2;
            }
          }
        }
        if (El_sourcet[k] == NULL)
        // element has had assigned no sources of any quantity
        {
          El_sourcet[k] = new long[ndof];
          memset(El_sourcet[k], 0, sizeof(*El_sourcet[k])*ndof);
        }
        if (srcid == 0)
        {
          // store source in the list of source
          El_srclst.append(tsrc);
          srcid = El_srclst.count();
        }
        El_sourcet[k][lcid] = srcid;
        // backup of line and column for eventual log message
        line[k] = aline;
        col[k]  = acol;
      }  // end of loop for element


      if (Check_unused_prop && (prop_used == 0))
      {
        // no elements were found with the given property and checking of unused properties is off
          continue;
        sprintf(errmsg, "Elements with property %ld belonging to entity %s have not been found.\n"
                        " Element source (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        delete tsrc;
        return 4;
      }
    }    // end of nkwd loop 
  }      // end of section loop
  return 0;
}



/**
  The function assigns boundary condition (BC) to the element nodes with required vertex property.
  It scans only section elvertpr for keyword "vertex_bc".
  The record which assigns edge BC to the element looks like follows:
  "vertex_bc" "propid" prop "lc_id" load_case_index "bc_type" boundary_condition_type {entity_bocon} x nbce 
  where prop is positive integer number, entity_bocon is record for entitybocon class and nbce is the
  number of prescribed quantities (nv, trc, trr) in the given bc_type. 
  For more details see entityload.cpp|h and loadelt.cpp|h
  
  It scans only section elvertpr. If an element belongs to several entities and different 
  BC types are prescribed for these entities, merging of BC is performed.
 
  @param in  - pointer to opened input file with property description
  @elemsects - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  @retval 0 - on success  
  @retval 1 - load case id is out of range
  @retval 2 - error in BC reading
  @retval 3 - failure of merging of BC or ceration of bnodvalt objects
  @retval 4 - unknown type of entity or problem type is required
  @retval 5 - property numbers of the required entity type have not been read on elements
  @retval 6 - property number of the required entity type have not been found on elements

  Created by TKo, 09.2010
*/
long input_elem_vertbc(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  const char *coeff_desc[3] = {"nodal values", "transmission coeff.", "radiation coeff."};
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k, l, m;
  entitybocon **tmp;
  climatcond *tmpcc;
  climatcond2 *tmpcc2;
  char *tmpccf;
  long appendcc;
  answertype *trcc;
  list *cclst = NULL;

  loadelt *tmp2;
  long check;
  long prop, lcid;
  long sres;
  long ndof;
  long prop_used, err;
  bocontypet bc;
  long **bnid;
  long aux, nbo, nbce;
  long minid_nv, minid_trc, minid_trr;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of vertex boundary conditions applied on elements");

  // number of dofs for particular nodes is the same and it equals to number of transported media
  ndof = Tp->ntm;
  if (El_loadt == NULL)
  {
    El_loadt = new loadelt**[ndof];
    memset(El_loadt, 0, sizeof(*El_loadt)*ndof);
    El_loadtln = new long *[ndof];
    memset(El_loadtln, 0, sizeof(*El_loadtln)*ndof);
    El_loadtcol = new long *[ndof];
    memset(El_loadtcol, 0, sizeof(*El_loadtcol)*ndof);
    for(i=0; i<ndof; i++)
    {
      El_loadt[i] = new loadelt*[Top->ne];
      memset(El_loadt[i], 0, sizeof(*El_loadt[i])*Top->ne);
      El_loadtln[i] = new long[Top->ne];
      memset(El_loadtln[i], 0, sizeof(*El_loadtln[i])*Top->ne);
      El_loadtcol[i] = new long[Top->ne];
      memset(El_loadtcol[i], 0, sizeof(*El_loadtcol[i])*Top->ne);
    }
    El_nv_lst   = new list[ndof];
    El_trc_lst  = new list[ndof];
    El_trr_lst  = new list[ndof];
    El_cc_lst   = new list[ndof];
    El_ccf_lst  = new list[ndof];
    El_trcc_lst = new list[ndof];
    El_gcc_lst   = new list[ndof];
    El_gccf_lst  = new list[ndof];
    El_gtrcc_lst = new list[ndof];
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gvertex; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "vertex_bc"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "vertex_bc", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "vertex_bc", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld%k%ld%k%m", "vertex_bc", "propid", &prop, "lc_id", &lcid, "bc_type", &bocontype_kwdset, &bc);
      if ((lcid < 1) && (lcid > ndof))
      {
        print_err("load case id is out range <1,%ld>", __FILE__, __LINE__, __func__, ndof);
        return 1;
      }
      lcid--;

      // determine number of necessary coefficients depending on the BC type
      if (bc < presc_trmiss)
        nbce = 1;
      else
        nbce = 3;

      appendcc = 0;
      trcc = NULL;
      tmpcc = NULL;
      tmpcc2 = NULL;
      tmpccf = NULL;
      tmp = NULL;
      if ((bc == det_climcond) || (bc == gen_climcond))
      {
        trcc = new answertype;
        xfscanf(in, "%k%m", "file_climcond", &answertype_kwdset, trcc);
        if (*trcc == no) 
        {
          if (bc == det_climcond)
          {
            tmpcc = new climatcond;
            tmpcc->read(in, 0);
            cclst = El_cc_lst+lcid;
          }
          else
          {
            tmpcc2 = new climatcond2;
            tmpcc2->read(in, 0);
            cclst = El_gcc_lst+lcid;
          }
        }
        else
        {
          tmpccf = new char[in->give_maxlnsize()+1];
          xfscanf(in, "%a", &tmpccf);
          if (bc == det_climcond)
            cclst = El_ccf_lst+lcid;
          else
            cclst = El_gccf_lst+lcid;
        }
        appendcc = 1;
      }
      else
      {
        tmp = new entitybocon*[nbce];
        memset(tmp, 0, sizeof(*tmp)*nbce);

        for(k=0; k<nbce; k++)
        {
          tmp[k] = new entitybocon;
          // read bocon record
          check = tmp[k]->read(in);
          switch(check)
          {
            case 0:
              break;
            case 1:
              print_err("cannot detect coordinates nor time in the expression\n for %s of the given BC", 
                        __FILE__, __LINE__, __func__, coeff_desc[k]);
              for(l=0; l<=k; l++)
                delete tmp[l];
              delete [] tmp;
              return 2;
            default:
              print_err("unknown error of BC reading", __FILE__, __LINE__, __func__);
              return 2;
          }
        }
      }

      in->kwdmode = bkwdmode;

      bnid = new long*[nbce];
      memset(bnid, 0, sizeof(*bnid)*nbce);
      nbo = 0;
      prop_used = 0;
      // storage of number of items in bnodvalt lists from the previous steps
      minid_nv  = El_nv_lst[lcid].count();
      minid_trc = El_trc_lst[lcid].count();
      minid_trr = El_trr_lst[lcid].count();

      // assign given BC type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        // BC can be applied to 1D elements only
        if (Top->elements[k].type >= trianglelinear)
          continue;

        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, Top->nodes, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { 
          // In this code, it is supposed that the given k-th element has the same number of nodes 
          // at each boundary object
          aux = get_nbo(Top->elements[k], gvertex);
          if (nbo != aux)
          {
            // prepare array of indices of bnodvalts for all type of prescribed coefficients
            nbo = aux;
            for (l=0; l<nbce; l++)
            {
              if (bnid[l])
                delete [] bnid[l];
              bnid[l] = new long[nbo];
            }
          }

          // check for multiple assignment of BC to one element
          if (El_loadt[lcid][k])
          {
            if (El_loadtln[lcid][k])
              fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned load at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, El_loadtln[lcid][k], El_loadtcol[lcid][k], in->fname);
            else
            { // element has had assigned different type of BC
              El_loadtln[lcid][k]  = aline;
              El_loadtcol[lcid][k] = acol;
            }
          }
          else
          {
            El_loadt[lcid][k] = new loadelt(k, nbo);
            // backup of line and column for eventual log message
            El_loadtln[lcid][k]  = aline;
            El_loadtcol[lcid][k] = acol;
          }

          // loop over number of required coefficients
          for (l=0; l< nbce; l++)
          {           
            // assembling of nodal values          
            if (l==0)
            {
              if ((bc == det_climcond) || (bc == gen_climcond))
              {
                if (appendcc)
                { 
                  if (bc == det_climcond)
                    El_trcc_lst[lcid].append(trcc);
                  else
                    El_gtrcc_lst[lcid].append(trcc);

                  if (*trcc == no)
                  {
                    if (bc == det_climcond)
                    {
                      El_cc_lst[lcid].append(tmpcc);
                      El_ccf_lst[lcid].append(NULL);
                    }
                    else
                    {
                      El_gcc_lst[lcid].append(tmpcc2);
                      El_gccf_lst[lcid].append(NULL);
                    }
                  }
                  else
                  {
                    if (bc == det_climcond)
                    {
                      El_ccf_lst[lcid].append(tmpccf);
                      El_cc_lst[lcid].append(NULL);
                    }
                    else
                    {
                      El_gccf_lst[lcid].append(tmpccf);
                      El_gcc_lst[lcid].append(NULL);
                    }
                  }
                  appendcc = 0;
                }
                assemble_bclimcond(Top->elements[k], Top->nodes, prop, *cclst, gvertex, bnid[0], entid, nentid);
              }
              else
              {
                err = assemble_bnodvalt(Top->elements[k], Top->nodes, prop, *tmp[0], El_nv_lst[lcid], gvertex, minid_nv, bnid[0], entid, nentid);
                if (err)
                {
                  print_err("cannot assemble nodal values for element %ld", __FILE__, __LINE__, __func__, k+1);
                  for (m=0; m<nbce; m++)
                  {
                    delete tmp[m];
                    delete [] bnid[m];
                  }
                  delete [] bnid;
                  delete [] tmp;
                  return 3;
                }
              }
            }
            // assembling of transmission coefficient
            if (l==1)
            {
              err = assemble_bnodvalt(Top->elements[k], Top->nodes, prop, *tmp[1], El_trc_lst[lcid], gvertex, minid_trc, bnid[1], entid, nentid);
              if (err)
              {
                print_err("cannot assemble transmission coefficients for element %ld", __FILE__, __LINE__, __func__, k+1);
                for (m=0; m<nbce; m++)
                {
                  delete tmp[m];
                  delete [] bnid[m];
                }
                delete [] bnid;
                delete [] tmp;
                return 3;
              }
            }
            // assembling of radiation coefficient
            if (l==2)
            {
              err = assemble_bnodvalt(Top->elements[k], Top->nodes, prop, *tmp[2], El_trr_lst[lcid], gvertex, minid_trr, bnid[2], entid, nentid);
              if (err)
              {
                print_err("cannot assemble radiation coefficients for element %ld", __FILE__, __LINE__, __func__, k+1);
                for (m=0; m<nbce; m++)
                {
                  delete tmp[m];
                  delete [] bnid[m];
                }
                delete [] bnid;
                delete [] tmp;
                return 3;
              }
            }
          }

          // create new object of loadelt from the actual BC
          tmp2 = bc2loadelt(k, Top->elements[k], Top->nodes, prop, bc, bnid, gvertex, entid, nentid);
          if (tmp2 == NULL)
          {
            print_err("unknown type of entity or problem is required", __FILE__, __LINE__, __func__);
            for (l=0; l<nbce; l++)
            {
              if (tmp)
                delete tmp[l];
              delete [] bnid[l];
            }
            delete [] bnid;
            if (tmp)
              delete [] tmp;
            return 3;
          }
          // BC merging
          check = El_loadt[lcid][k]->merge(*tmp2);
          prop_used = 1;
          switch(check)
          {
            case 0:
              break;
            case 1:
              // element has already assigned BC with different values
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned boundary condition"
                              "with different values at line %ld, col %ld, file: %s\n", 
                              aline, acol, in->fname, k+1, El_loadtln[lcid][k], El_loadtcol[lcid][k], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              for (l=0; l<nbce; l++)
              {
                if (tmp)
                  delete tmp[l];
                delete [] bnid[l];
              }
              delete [] bnid;
              if (tmp)
                delete [] tmp;
              delete tmp2;
              return 4;
            case 2:
              // element has already assigned different type of BC
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned different "
                              "type of boundary condition at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_loadtln[lcid][k], El_loadtcol[lcid][k], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              for (l=0; l<nbce; l++)
              {
                if (tmp)
                  delete tmp[l];
                delete [] bnid[l];
              }
              delete [] bnid;
              if (tmp)
                delete [] tmp;
              delete tmp2;
              return 4;
            default:
              // element has already assigned BC which cannot be merged
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned boundary condition "
                              "which cannot be merged at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_loadtln[lcid][k], El_loadtcol[lcid][k], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              for (l=0; l<nbce; l++)
              {
                if (tmp)
                  delete tmp[l];
                delete [] bnid[l];
              }
              delete [] bnid;
              if (tmp)
                delete [] tmp;
              delete tmp2;
              return 4;
          }
          delete tmp2;
        } // end of if (sres > 0)
        if (sres < 0)
        {
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 5;
        }
      } // end of loop k
      for (l=0; l<nbce; l++)
      {
        if (tmp)
          delete tmp[l];
        delete [] bnid[l];
      }
      delete [] bnid;
      if (tmp)
        delete [] tmp;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no nodes were found with the given property and checking of unused properties is on
        sprintf(errmsg, "There are no 1D elements with property %ld belonging to entity %s in the topology.\n"
                        " Element edge boundary condition (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 6;
      }
    }
  }

  return 0;
}



/**
  The function assigns edge boundary condition (BC) to the elements with with required edge property.
  It scans section eledgpr for keyword "edge_bc".
  The record which assigns edge BC to the element looks like follows:
  "edge_bc" "propid" prop "lc_id" load_case_index "bc_type" boundary_condition_type {entity_bocon} x nbce 
  where prop is positive integer number, entity_bocon is record for entitybocon class and nbce is the
  number of prescribed quantities (nv, trc, trr) in the given bc_type. 
  For more details see entityload.cpp|h and loadelt.cpp|h
  
  It scans only section for edges. If an element  belongs to several entities and different 
  BC types are prescribed for these entities, merging of BC is performed.
 
  @param in  - pointer to opened input file with property description
  @elemsects - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  @retval 0 - on success  
  @retval 1 - load case id is out of range
  @retval 2 - error in BC reading
  @retval 3 - failure of merging of BC or ceration of bnodvalt objects
  @retval 4 - unknown type of entity or problem type is required
  @retval 5 - property numbers of the required entity type have not been read on elements
  @retval 6 - property number of the required entity type have not been found on elements

  Created by TKo, 09.2010
*/
long input_elem_edgebc(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  const char *coeff_desc[3] = {"nodal values", "transmission coeff.", "radiation coeff."};
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k, l, m;
  entitybocon **tmp;
  climatcond *tmpcc;
  climatcond2 *tmpcc2;
  char *tmpccf;
  long appendcc;
  answertype *trcc;
  list *cclst = NULL;

  loadelt *tmp2;
  long check;
  long prop, lcid;
  long sres;
  long ndof;
  long prop_used, err;
  bocontypet bc;
  long **bnid;
  long aux, nbo, nbce;
  long minid_nv, minid_trc, minid_trr;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of edge boundary conditions applied on elements");

  // number of dofs for particular nodes is the same and it equals to number of transported media
  ndof = Tp->ntm;


  if (El_loadt == NULL)
  {
    El_loadt = new loadelt**[ndof];
    memset(El_loadt, 0, sizeof(*El_loadt)*ndof);
    El_loadtln = new long *[ndof];
    memset(El_loadtln, 0, sizeof(*El_loadtln)*ndof);
    El_loadtcol = new long *[ndof];
    memset(El_loadtcol, 0, sizeof(*El_loadtcol)*ndof);
    for(i=0; i<ndof; i++)
    {
      El_loadt[i] = new loadelt*[Top->ne];
      memset(El_loadt[i], 0, sizeof(*El_loadt[i])*Top->ne);
      El_loadtln[i] = new long[Top->ne];
      memset(El_loadtln[i], 0, sizeof(*El_loadtln[i])*Top->ne);
      El_loadtcol[i] = new long[Top->ne];
      memset(El_loadtcol[i], 0, sizeof(*El_loadtcol[i])*Top->ne);
    }
    El_nv_lst   = new list[ndof];
    El_trc_lst  = new list[ndof];
    El_trr_lst  = new list[ndof];
    El_cc_lst   = new list[ndof];
    El_ccf_lst  = new list[ndof];
    El_trcc_lst = new list[ndof];
    El_gcc_lst   = new list[ndof];
    El_gccf_lst  = new list[ndof];
    El_gtrcc_lst = new list[ndof];
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gcurve; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "edge_bc"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "edge_bc", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "edge_bc", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld%k%ld%k%m", "edge_bc", "propid", &prop, "lc_id", &lcid, "bc_type", &bocontype_kwdset, &bc);
      if ((lcid < 1) && (lcid > ndof))
      {
        print_err("load case id is out range <1,%ld>", __FILE__, __LINE__, __func__, ndof);
        return 1;
      }
      lcid--;

      // determine number of necessary coefficients depending on the BC type
      if (bc < presc_trmiss)
        nbce = 1;
      else
        nbce = 3;

      appendcc = 0;
      trcc = NULL;
      tmpcc = NULL;
      tmpcc2 = NULL;
      tmpccf = NULL;
      tmp = NULL;
      if ((bc == det_climcond) || (bc == gen_climcond))
      {
        trcc = new answertype;
        xfscanf(in, "%k%m", "file_climcond", &answertype_kwdset, trcc);
        if (*trcc == no) 
        {
          if (bc == det_climcond)
          {
            tmpcc = new climatcond;
            tmpcc->read(in, 0);
            cclst = El_cc_lst+lcid;
          }
          else
          {
            tmpcc2 = new climatcond2;
            tmpcc2->read(in, 0);
            cclst = El_gcc_lst+lcid;
          }
        }
        else
        {
          tmpccf = new char[in->give_maxlnsize()+1];
          xfscanf(in, "%a", &tmpccf);
          if (bc == det_climcond)
            cclst = El_ccf_lst+lcid;
          else
            cclst = El_gccf_lst+lcid;
        }
        appendcc = 1;
      }
      else
      {
        tmp = new entitybocon*[nbce];
        memset(tmp, 0, sizeof(*tmp)*nbce);
      
        for(k=0; k<nbce; k++)
        {
          tmp[k] = new entitybocon;
          // read bocon record
          check = tmp[k]->read(in);
          switch(check)
          {
            case 0:
              break;
            case 1:
              print_err("cannot detect coordinates nor time in the expression\n for %s of the given BC", 
                        __FILE__, __LINE__, __func__, coeff_desc[k]);
              for(l=0; l<=k; l++)
                delete tmp[l];
              delete [] tmp;
              return 2;
            default:
              print_err("unknown error of BC reading", __FILE__, __LINE__, __func__);
              return 2;
          }
        }
      }

      in->kwdmode = bkwdmode;
    
      bnid = new long*[nbce];
      memset(bnid, 0, sizeof(*bnid)*nbce);
      nbo = 0;
      prop_used = 0;
      // storage of number of items in bnodvalt lists from the previous steps
      minid_nv  = El_nv_lst[lcid].count();
      minid_trc = El_trc_lst[lcid].count();
      minid_trr = El_trr_lst[lcid].count();

      // assign given BC type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        // BC can be applied to 2D elements only
        if (Top->elements[k].type >= tetrahedronlinear)
          continue;

        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, Top->nodes, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { 
          // In this code, it is supposed that the given k-th element has the same number of nodes 
          // at each boundary object
          aux = get_nbo(Top->elements[k], gcurve);
          if (nbo != aux)
          {
            // prepare array of indices of bnodvalts for all type of prescribed coefficients
            nbo = aux;
            for (l=0; l<nbce; l++)
            {
              if (bnid[l])
                delete [] bnid[l];
              bnid[l] = new long[nbo];
            }
          }

          // check for multiple assignment of BC to one element
          if (El_loadt[lcid][k])
          {
            if (El_loadtln[lcid][k])
              fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned load at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, El_loadtln[lcid][k], El_loadtcol[lcid][k], in->fname);
            else
            { // element has had assigned different type of BC
              El_loadtln[lcid][k]  = aline;
              El_loadtcol[lcid][k] = acol;
            }
          }
          else
          {
            El_loadt[lcid][k] = new loadelt(k, nbo);
            // backup of line and column for eventual log message
            El_loadtln[lcid][k]  = aline;
            El_loadtcol[lcid][k] = acol;
          }

          // loop over number of required coefficients
          for (l=0; l< nbce; l++)
          {           
            // assembling of nodal values          
            if (l==0)
            {
              if ((bc == det_climcond) || (bc == gen_climcond))
              {
                if (appendcc)
                { 
                  if (bc == det_climcond)
                    El_trcc_lst[lcid].append(trcc);
                  else
                    El_gtrcc_lst[lcid].append(trcc);
                
                  if (*trcc == no)
                  {
                    if (bc == det_climcond)
                    {
                      El_cc_lst[lcid].append(tmpcc);
                      El_ccf_lst[lcid].append(NULL);
                    }
                    else
                    {
                      El_gcc_lst[lcid].append(tmpcc2);
                      El_gccf_lst[lcid].append(NULL);
                    }
                  }
                  else
                  {
                    if (bc == det_climcond)
                    {
                      El_ccf_lst[lcid].append(tmpccf);
                      El_cc_lst[lcid].append(NULL);
                    }
                    else
                    {
                      El_gccf_lst[lcid].append(tmpccf);
                      El_gcc_lst[lcid].append(NULL);
                    }
                  }
                  appendcc = 0;
                }
                assemble_bclimcond(Top->elements[k], Top->nodes, prop, *cclst, gcurve, bnid[0], entid, nentid);
              }
              else
              {
                err = assemble_bnodvalt(Top->elements[k], Top->nodes, prop, *tmp[0], El_nv_lst[lcid], gcurve, minid_nv, bnid[0], entid, nentid);
                if (err)
                {
                  print_err("cannot assemble nodal values for element %ld", __FILE__, __LINE__, __func__, k+1);
                  for (m=0; m<nbce; m++)
                  {
                    delete tmp[m];
                    delete [] bnid[m];
                  }
                  delete [] bnid;
                  delete [] tmp;
                  return 3;
                }
              }
            }
            // assembling of transmission coefficient
            if (l==1)
            {
              err = assemble_bnodvalt(Top->elements[k], Top->nodes, prop, *tmp[1], El_trc_lst[lcid], gcurve, minid_trc, bnid[1], entid, nentid);
              if (err)
              {
                print_err("cannot assemble transmission coefficients for element %ld", __FILE__, __LINE__, __func__, k+1);
                for (m=0; m<nbce; m++)
                {
                  delete tmp[m];
                  delete [] bnid[m];
                }
                delete [] bnid;
                delete [] tmp;
                return 3;
              }
            }
            // assembling of radiation coefficient
            if (l==2)
            {
              err = assemble_bnodvalt(Top->elements[k], Top->nodes, prop, *tmp[2], El_trr_lst[lcid], gcurve, minid_trr, bnid[2], entid, nentid);
              if (err)
              {
                print_err("cannot assemble radiation coefficients for element %ld", __FILE__, __LINE__, __func__, k+1);
                for (m=0; m<nbce; m++)
                {
                  delete tmp[m];
                  delete [] bnid[m];
                }
                delete [] bnid;
                delete [] tmp;
                return 3;
              }
            }
          }

          // create new object of loadelt from the actual BC
          tmp2 = bc2loadelt(k, Top->elements[k], Top->nodes, prop, bc, bnid, gcurve, entid, nentid);
          if (tmp2 == NULL)
          {
            print_err("unknown type of entity or problem is required", __FILE__, __LINE__, __func__);
            for (l=0; l<nbce; l++)
            {
              if (tmp)
                delete tmp[l];
              delete [] bnid[l];
            }
            delete [] bnid;
            if (tmp)
              delete [] tmp;
            return 3;
          }
          // BC merging
          check = El_loadt[lcid][k]->merge(*tmp2);
          prop_used = 1;
          switch(check)
          {
            case 0:
              break;
            case 1:
              // element has already assigned BC with different values
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned boundary condition"
                              "with different values at line %ld, col %ld, file: %s\n", 
                              aline, acol, in->fname, k+1, El_loadtln[lcid][k], El_loadtcol[lcid][k], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              for (l=0; l<nbce; l++)
              {
                if (tmp)
                  delete tmp[l];
                delete [] bnid[l];
              }
              delete [] bnid;
              if (tmp)
                delete [] tmp;
              delete tmp2;
              return 4;
            case 2:
              // element has already assigned different type of BC
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned different "
                              "type of boundary condition at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_loadtln[lcid][k], El_loadtcol[lcid][k], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              for (l=0; l<nbce; l++)
              {
                if (tmp)
                  delete tmp[l];
                delete [] bnid[l];
              }
              delete [] bnid;
              if (tmp)
                delete [] tmp;
              delete tmp2;
              return 4;
            default:
              // element has already assigned BC which cannot be merged
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned boundary condition "
                              "which cannot be merged at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_loadtln[lcid][k], El_loadtcol[lcid][k], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              for (l=0; l<nbce; l++)
              {
                if (tmp)
                  delete tmp[l];
                delete [] bnid[l];
              }
              delete [] bnid;
              if (tmp)
                delete [] tmp;
              delete tmp2;
              return 4;
          }
          delete tmp2;
        } // end of if (sres > 0)
        if (sres < 0)
        {
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 5;
        }
      } // end of loop k
      for (l=0; l<nbce; l++)
      {
        if (tmp)
          delete tmp[l];
        delete [] bnid[l];
      }
      delete [] bnid;
      if (tmp)
        delete [] tmp;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no nodes were found with the given property and checking of unused properties is on
        sprintf(errmsg, "There are no 2D elements with property %ld belonging to entity %s in the topology.\n"
                        " Element edge boundary condition (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 6;
      }
    }
  }

  return 0;
}



/**
  The function assigns surface boundary condition (BC) to the elements with with required surface property.
  It scans  sections elsurfpr for keyword "surf_bc".
  The record assigns surface BC to the element looks like follows:
  "surf_bc" "propid" prop "lc_id" load_case_index "bc_type" boundary_condition_type {entity_bocon} x nbce 
  where prop is positive integer number, entity_bocon is record for entitybocon class and nbce is the
  number of prescribed quantities (nv, trc, trr) in the given bc_type. 
  For more details see entityload.cpp|h and loadelt.cpp|h
  
  It scans only section for surfaces. If an element  belongs to several entities and different 
  BC types are prescribed for these entities, merging of BC is performed.
 
  @param in  - pointer to opened input file with property description
  @elemsects - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  @retval 0 - on success  
  @retval 1 - load case id is out of range
  @retval 2 - error in BC reading
  @retval 3 - failure of merging of BC or bnodvalt object
  @retval 4 - unknown type of entity or problem type is required
  @retval 5 - property numbers of the required entity type have not been read on elements
  @retval 6 - property number of the required entity type have not been found on elements

  Created by TKo, 09.2010
*/
long input_elem_surfbc(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  char errmsg[1001];
  const char *coeff_desc[3] = {"nodal values", "transmission coeff.", "radiation coeff."};
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k, l, m;
  entitybocon **tmp;
  climatcond  *tmpcc;
  climatcond2 *tmpcc2;
  char *tmpccf;
  long appendcc;
  answertype *trcc;
  list *cclst = NULL;

  loadelt *tmp2;
  long check;
  long prop, lcid;
  long sres;
  long ndof;
  long prop_used;
  bocontypet bc;
  long **bnid;
  long aux, nbo, nbce, err;
  long minid_nv, minid_trc, minid_trr;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of surface boundary conditions applied on elements");

  // number of dofs for particular nodes is the same and it equals to number of transported media
  ndof = Tp->ntm;
  if (El_loadt == NULL)
  {
    El_loadt = new loadelt**[ndof];
    memset(El_loadt, 0, sizeof(*El_loadt)*ndof);
    El_loadtln = new long *[ndof];
    memset(El_loadtln, 0, sizeof(*El_loadtln)*ndof);
    El_loadtcol = new long *[ndof];
    memset(El_loadtcol, 0, sizeof(*El_loadtcol)*ndof);
    for(i=0; i<ndof; i++)
    {
      El_loadt[i] = new loadelt*[Top->ne];
      memset(El_loadt[i], 0, sizeof(*El_loadt[i])*Top->ne);
      El_loadtln[i] = new long[Top->ne];
      memset(El_loadtln[i], 0, sizeof(*El_loadtln[i])*Top->ne);
      El_loadtcol[i] = new long[Top->ne];
      memset(El_loadtcol[i], 0, sizeof(*El_loadtcol[i])*Top->ne);
    }
    El_nv_lst   = new list[ndof];
    El_trc_lst  = new list[ndof];
    El_trr_lst  = new list[ndof];
    El_cc_lst   = new list[ndof];
    El_ccf_lst  = new list[ndof];
    El_trcc_lst = new list[ndof];
    El_gcc_lst   = new list[ndof];
    El_gccf_lst  = new list[ndof];
    El_gtrcc_lst = new list[ndof];
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gsurface; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "surf_bc"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "surf_bc", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "surf_bc", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld%k%ld%k%m", "surf_bc", "propid", &prop, "lc_id", &lcid, "bc_type", &bocontype_kwdset, &bc);
      if ((lcid < 1) && (lcid > ndof))
      {
        print_err("load case id is out range <1,%ld>", __FILE__, __LINE__, __func__, ndof);
        return 1;
      }
      lcid--;

      // determine number of necessary coefficients depending on the BC type
      if (bc < presc_trmiss)
        nbce = 1;
      else
        nbce = 3;

      appendcc = 0;
      trcc = NULL;
      tmpcc = NULL;
      tmpcc2 = NULL;
      tmpccf = NULL;
      tmp = NULL;

      if ((bc == det_climcond) || (bc == gen_climcond))
      {
        tmp = NULL;
        trcc = new answertype;
        xfscanf(in, "%k%m", "file_climcond", &answertype_kwdset, trcc);
        if (*trcc == no) 
        {
          if (bc == det_climcond)
          {
            tmpcc = new climatcond;
            tmpcc->read(in, 0);
            cclst = &El_cc_lst[lcid];
          }
          else
          {
            tmpcc2 = new climatcond2;
            tmpcc2->read(in, 0);
            cclst = &El_gcc_lst[lcid];
          }
        }
        else
        {
          tmpccf = new char[in->give_maxlnsize()];
          xfscanf(in, "%a", &tmpccf);
          if (bc == det_climcond)
            cclst = El_ccf_lst+lcid;
          else
            cclst = El_gccf_lst+lcid;
        }
        appendcc = 1;
      }
      else
      {
        tmp = new entitybocon*[nbce];
        memset(tmp, 0, sizeof(*tmp)*nbce);

        for(k=0; k<nbce; k++)
        {
          tmp[k] = new entitybocon;
          // read bocon record
          check = tmp[k]->read(in);
          switch(check)
          {
            case 0:
              break;
            case 1:
              print_err("cannot detect coordinates nor time in the expression\n for %s of the given BC", 
                        __FILE__, __LINE__, __func__, coeff_desc[k]);
              for(l=0; l<=k; l++)
                delete tmp[l];
              delete [] tmp;
              return 2;
            default:
              print_err("unknown error of BC reading", __FILE__, __LINE__, __func__);
              return 2;
          }
        }
      }

      in->kwdmode = bkwdmode;

      bnid = new long*[nbce];
      memset(bnid, 0, sizeof(*bnid)*nbce);
      nbo = 0;
      prop_used = 0;
      // storage of number of items in bnodvalt lists from the previous steps
      minid_nv  = El_nv_lst[lcid].count();
      minid_trc = El_trc_lst[lcid].count();
      minid_trr = El_trr_lst[lcid].count();

      // assign given BC type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        // BC can be applied to 3D elements only
        if (Top->elements[k].type < tetrahedronlinear)
          continue;

        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, Top->nodes, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { 
          // In this code, it is supposed that the given k-th element has the same number of nodes 
          // at each boundary object

          // prepare array of indeces of bnodvalts for all type of prescribed coefficients
          aux = get_nbo(Top->elements[k], gsurface);
          if (nbo != aux)
          {
            // prepare array of indeces of bnodvalts for all type of prescribed coefficients
            nbo = aux;
            for (l=0; l<nbce; l++)
            {
              if (bnid[l])
                delete [] bnid[l];
              bnid[l] = new long[nbo];
            }
          }
	  // check for multiple assignment of BC to one element
          if (El_loadt[lcid][k])
          {
            if (El_loadtln[lcid][k])
              fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has already assigned load at line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, El_loadtln[lcid][k], El_loadtcol[lcid][k], in->fname);
            else
            { // element has had assigned different type of BC
              El_loadtln[lcid][k]  = aline;
              El_loadtcol[lcid][k] = acol;
            }
          }
          else
          {
            El_loadt[lcid][k] = new loadelt(k, nbo);
            // backup of line and column for eventual log message
            El_loadtln[lcid][k]  = aline;
            El_loadtcol[lcid][k] = acol;
          }

          // loop over all required coefficients
          for (l=0; l<nbce; l++)
          {
            // assembling of nodal values          
            if (l==0)
            {
              if ((bc == det_climcond) || (bc == gen_climcond))
              {
                if (appendcc)
                {
                  if (bc == det_climcond)
                    El_trcc_lst[lcid].append(trcc);
                  else
                    El_gtrcc_lst[lcid].append(trcc);

                  if (*trcc == no)
                  {
                     
                    if (bc == det_climcond)
                    {
                      El_cc_lst[lcid].append(tmpcc);
                      El_ccf_lst[lcid].append(NULL);
                    }
                    else
                    {
                      El_gcc_lst[lcid].append(tmpcc2);
                      El_gccf_lst[lcid].append(NULL);
                    }
                  }
                  else
                  {
                    if (bc == det_climcond)
                    {
                      El_ccf_lst[lcid].append(tmpccf);
                      El_cc_lst[lcid].append(NULL);
                    }
                    else
                    {
                      El_gccf_lst[lcid].append(tmpccf);
                      El_gcc_lst[lcid].append(NULL);
                    }
                  }
                  appendcc = 0;
                }
                assemble_bclimcond(Top->elements[k], Top->nodes, prop, *cclst, gsurface, bnid[0], entid, nentid);
              }
              else
              {
                err = assemble_bnodvalt(Top->elements[k], Top->nodes, prop, *tmp[0], El_nv_lst[lcid], gsurface, minid_nv, bnid[0], entid, nentid);
                if (err)
                {
                  print_err("cannot assemble nodal values for element %ld", __FILE__, __LINE__, __func__, k+1);
                  for (m=0; m<nbce; m++)
                  {
                    delete tmp[m];
                    delete [] bnid[m];
                  }
                  delete [] bnid;
                  delete [] tmp;
                  return 3;
                }
              }
            }
            // assembling of transmission coefficient
            if (l==1)
            {
              err = assemble_bnodvalt(Top->elements[k], Top->nodes, prop, *tmp[1], El_trc_lst[lcid], gsurface, minid_trc, bnid[1], entid, nentid);
              if (err)
              {
                print_err("cannot assemble transmission coefficients for element %ld", __FILE__, __LINE__, __func__, k+1);
                for (m=0; m<nbce; m++)
                {
                  delete tmp[m];
                  delete [] bnid[m];
                }
                delete [] bnid;
                delete [] tmp;
                return 3;
              }
            }
            // assembling of radiation coefficient
            if (l==2)
            {
              err = assemble_bnodvalt(Top->elements[k], Top->nodes, prop, *tmp[2], El_trr_lst[lcid], gsurface, minid_trr, bnid[2], entid, nentid);
              if (err)
              {
                print_err("cannot assemble radiation coefficients for element %ld", __FILE__, __LINE__, __func__, k+1);
                for (m=0; m<nbce; m++)
                {
                  delete tmp[m];
                  delete [] bnid[m];
                }
                delete [] bnid;
                delete [] tmp;
                return 3;
              }
            }
          }

          tmp2 = bc2loadelt(k, Top->elements[k], Top->nodes, prop, bc, bnid, gsurface, entid, nentid);
          if (tmp2 == NULL)
          {
            print_err("unknown type of entity or problem is required", __FILE__, __LINE__, __func__);
            for (l=0; l<nbce; l++)
            {
              if (tmp)
                delete tmp[l];
              delete [] bnid[l];
            }
            delete [] bnid;
            if (tmp)
              delete [] tmp;
            return 3;
          }
          // BC merging
          check = El_loadt[lcid][k]->merge(*tmp2);
          prop_used = 1;
          switch(check)
          {
            case 0:
              break;
            case 1:
              // element has already assigned BC with different values
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned boundary condition"
                              "with different values at line %ld, col %ld, file: %s\n", 
                              aline, acol, in->fname, k+1, El_loadtln[lcid][k], El_loadtcol[lcid][k], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              for (l=0; l<nbce; l++)
              {
                if (tmp)
                  delete tmp[l];
                delete [] bnid[l];
              }
              delete [] bnid;
              if (tmp)
                delete [] tmp;
              delete tmp2;
              return 4;
            case 2:
              // element has already assigned different type of BC
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned different "
                              "type of boundary condition at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_loadtln[lcid][k], El_loadtcol[lcid][k], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              for (l=0; l<nbce; l++)
              {
                if (tmp)
                  delete tmp[l];
                delete [] bnid[l];
              }
              delete [] bnid;
              if (tmp)
                delete [] tmp;
              delete tmp2;
              return 4;
            default:
              // element has already assigned BC which cannot be merged
              sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned boundary condition "
                              "which cannot be merged at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                              k+1, El_loadtln[lcid][k], El_loadtcol[lcid][k], in->fname);
              print_err(errmsg, __FILE__, __LINE__, __func__);
              for (l=0; l<nbce; l++)
              {
                if (tmp)
                  delete tmp[l];
                delete [] bnid[l];
              }
              delete [] bnid;
              if (tmp)
                delete [] tmp;
              delete tmp2;
              return 4;
          }
          delete tmp2;
        }        
        if (sres < 0)
        {
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 5;
        }
      }
      for (l=0; l<nbce; l++)
      {
        if (tmp)
          delete tmp[l];
        delete [] bnid[l];
      }
      delete [] bnid;
      if (tmp)
        delete [] tmp;
      if (Check_unused_prop && (prop_used == 0))
      {
        // no nodes were found with the given property and checking of unused properties is on
        sprintf(errmsg, "There are no 3D elements with property %ld belonging to entity %s in the topology.\n"
                        " Element surface boundary condition (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 6;
      }
    }
  }

  return 0;
}



/**
  The function marks element surfaces where integration of selected quantity fluxes will be integrated over the given surface area
  It scans  sections elsurfpr for keyword "flux_int" whose record looks like follows:
  "flux_int" "propid" prop "lc_id" load_case_index "intdom_id" integration_domain_index
  where prop, loda_caes_index and integration_domain_index are positive integer number. The quantity, whose flux will be integrated, 
  is given by load_case_index, surface integral bounds are given by surface with the given property prop. Only section for surfaces 
  are scanned by this function. 
 
  @param in  - pointer to opened input file with property description
  @elemsects - array with descriptors of sections, which will be searched
  @nsect     - number of searched sections i.e number of elements in array elemsects

  @retval 0 - on success  
  @retval 1 - load case id is out of range
  @retval 2 - error in the merging of integration domain indeces
  @retval 3 - property numbers of the required entity type have not been read on elements
  @retval 4 - property number of the required entity type have not been found on elements

  Created by TKo, 02.2018
*/
long input_elem_fluxint(XFILE *in, const enumstr elemsects[], long nsect)
{
  char *aptr;
  long nkwd;
  ivector *line = NULL;
  ivector *col  = NULL;
  long aline, acol;
  gentity ent;
  kwd_handling bkwdmode;
  long i, j, k;
  intdomt *idom;

  long prop, lcid, idomid;
  long sres, ret;
  long ndof;
  long prop_used;
  long *entid;
  long nentid;

  fprintf(stdout, "\n reading of domain description of flux integrals on elements");

  // number of dofs for particular nodes is the same and it equals to number of transported media
  ndof = Tp->ntm;
  Nidomid = 0;
  El_fluxint = new list**[ndof];
  memset(El_fluxint, 0, sizeof(*El_fluxint)*ndof);
  line = new ivector[ndof];
  col  = new ivector[ndof];
  for(i=0; i<ndof; i++)
  {
    El_fluxint[i] = new list*[Top->ne];
    memset(El_fluxint[i], 0, sizeof(*El_fluxint[i])*Top->ne);
    reallocv(Top->ne, line[i]);
    reallocv(Top->ne, col[i]);
  }

  bkwdmode       = in->kwdmode;
  in->ignorecase = 1;

  for (i=nsect-1, ent=gsurface; i>=0; i--, ent--)
  {
    // set searched section
    if (xf_setsec(in, elemsects[i]) == 1)
      continue;
    // detect number of occurences of keyword "flux_int"
    in->kwdmode = sect_mode_full;
    nkwd = getkwd_sect(in, NULL, aptr, "flux_int", 1);
    in->kwdmode = bkwdmode;
    // process each keyword occurence
    for(j=0; j < nkwd; j++)
    {
      // go to the next keyword occurence
      in->kwdmode = sect_mode_fwd;
      getkwd_sect(in, NULL, aptr, "flux_int", 0);
      in->kwdmode = bkwdmode;
      // backup actual line and column
      aline = in->line;
      acol  = in->col;
      in->kwdmode = sect_mode_seq;
      // read ndof and property id
      xfscanf(in, "%k%k%ld%k%ld%k%ld", "flux_int", "propid", &prop, "lc_id", &lcid, "intdom_id", &idomid);
      if ((lcid < 1) && (lcid > ndof))
      {
        print_err("load case id is out range <1,%ld>", __FILE__, __LINE__, __func__, ndof);
        delete [] line;
        delete [] col;
        return 1;
      }
      lcid--;
      if (idomid > Nidomid)
        Nidomid = idomid;

      in->kwdmode = bkwdmode;
      prop_used = 0;

      // assign given BC type to elements with given property id
      for(k=0; k<Top->ne; k++)
      {    
        // BC can be applied to 3D elements only
        if (Top->elements[k].type < tetrahedronlinear)
          continue;

        entid = NULL;
        nentid = 0;
        sres = Top->elements[k].searchprop(prop, ent, Top->nodes, k, Top->edges, Top->surfaces, entid, nentid);
        if (sres > 0)
        { 
	  // check for multiple assignment of flux_int to one element
          if (El_fluxint[lcid][k])
          {
            fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has been already involved into flux integration line %ld, col %ld, file: %s\n", aline, acol, in->fname, k+1, line[lcid][k], col[lcid][k], in->fname);
            idom = new intdomt;
            assemble_int_domain(Top->elements[k], Top->nodes, prop, gsurface, entid, nentid, idomid, *idom);
            ret = merge_intdom(El_fluxint[lcid][k], idom);
            if (ret == 0)
            {
              // all surface integration indeces in idom could be merged with the previous integration descriptors (items of list El_fluxint[lcid][k])
              fprintf(Log, "Line %ld, col %ld, file %s: Element %ld has been already involved into flux integration on domain id=%ld\n", aline, acol, in->fname, k+1, idomid);
              delete idom;
            }
            if (ret == 2)
            {
              print_err("incompatible numbers of boundary objects of merged integration domain descriptors on element %ld.\n" 
                        "see line %ld, col %ld, file: %s\n", __FILE__, __LINE__, __func__, k+1, aline, acol, in->fname);
              delete [] line;
              delete [] col;
              return 2;
            }
            prop_used = 1;
          }
          else
          {
            El_fluxint[lcid][k] = new list;
            idom = new intdomt;
            assemble_int_domain(Top->elements[k], Top->nodes, prop, gsurface, entid, nentid, idomid, *idom);
            El_fluxint[lcid][k]->append(idom);
            // backup of line and column for eventual log message
            line[lcid][k]  = aline;
            col[lcid][k] = acol;
            prop_used = 1;
          }
        }
        if (sres < 0)
        {
          print_err("Property numbers of %s are not read on element %ld,\n"
                    " required property cannot be assigned.\n"
                    " see input file %s, line %ld, column %ld", __FILE__, __LINE__, __func__, entitypstr[ent-1].alias, k+1, in->fname, aline, acol);
          delete [] line;
          delete [] col;
          return 3;
        }
      }
      if (Check_unused_prop && (prop_used == 0))
      {
        // no nodes were found with the given property and checking of unused properties is on
        print_err("Elements with property %ld belonging to entity %s have not been found.\n"
                  " Element flux integration on surface (line %ld, column %ld, file %s) cannot be assigned correctly.",
                   __FILE__, __LINE__, __func__, 
                   prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        delete [] line;
        delete [] col;
        return 4;
      }
    }
  }

  delete [] line;
  delete [] col;

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
  @retval 0 - on success  
  @retval 1 - assignment of different time function at one element
  @retval 2 - property numbers of the required entity type have not been read on elements
  @retval 3 - element have not assigned a time function
  @retval 4 - property number of the required entity type have not been found on elements
*/
long input_elem_eltimefunct(XFILE *in, const enumstr elemsects[], long nsect)
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

  fprintf(stdout, "\n reading of dof time functions at nodes");

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
      if ((tf_id <= 0) || (tf_id > Gtt->ngf))
      {
        sprintf(errmsg, "Line %ld, col %ld, file %s:\n Time function index has to be in range <1,%ld>"
                        "time function at line %ld, col %ld, file: %s\n", aline, acol, in->fname, Gtt->ngf, 
                        line[k], col[k], in->fname);
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
        sres = Top->elements[k].searchprop(prop, ent, Top->nodes, k, Top->edges, Top->surfaces, entid, nentid);    
        if (sres > 0)
        { // check for multiple assignment of different time functions to one element
          if (El_tfunc[k] && (El_tfunc[k]!=tf_id))
          {
            sprintf(errmsg, "Line %ld, col %ld, file %s:\n Element %ld has already had assigned different"
                    "time function at line %ld, col %ld, file: %s\n", aline, acol, in->fname,
                            k+1, line[k], col[k], in->fname);
            print_err(errmsg, __FILE__, __LINE__, __func__);
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
          sprintf(errmsg, "Element %ld has not read property numbers of %s,\n"
                          " required property cannot be assigned.\n"
                          " see input file %s, line %ld, column %ld", k+1, entitypstr[ent-1].alias, in->fname, aline, acol);
          print_err(errmsg, __FILE__, __LINE__, __func__);
          return 2;
        }
      }
      if (Check_unused_prop && (prop_used == 0))
      {
        // no nodes were found with the given property and checking of unused properties is on
        sprintf(errmsg, "Elements with property %ld belonging to entity %s have not been found.\n"
                        " Element time function (line %ld, column %ld, file %s) cannot be assigned correctly.", 
                prop, entitypstr[ent-1].alias, aline, acol, in->fname);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 4;
      }
    }
  }
  // in case of growing transport problems, all elements have to have 
  // assigned a time function
  if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
  {
    for (i=0; i< Top->ne; i++)
    {
      if (El_tfunc[i] < 1)
      {
        sprintf(errmsg, "Element %ld has not assigned a time function.\n"
                        "All elements have to have assigned some time function\n"
                        "in case of growing transport problem", i+1);
        print_err(errmsg, __FILE__, __LINE__, __func__);
        return 3;
      }
    }
  }
  return 0;
}



/**
  Overloaded postfix decrementation operator -- for enum gentity
  
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
  The function fills array of indeces bnid by the indeces in the list of bnodvalt object for the given element and boundary
  condition (BC).

  @param el    - element to which the BC will be assigned
  @param nodes - array of all nodes in the given topology
  @param prop  - entity property id of the element el where the entity BC will be applied
  @param ebc   - object with BC for the given entity ent
  @param bnvl  - list of assigned bnodvalt objects
  @param ent   - entity type for which the the BC will be assigned
  @param minid - the minimum index from which the list bnvl is searched for the identical bnodvalt objects
                 of the given BC. It is used for reduction of ammount of bnodvalt objects.
  @param bnid - array with indeces of bnodval objects for the given coefficients (numbering of indeces must be from zero)
                bnid[i] returns index for i-th boundary object of the given coefficient prescribed by ebc
  @param entid - array of entity indeces for the given property number detected on the additional edges or surfaces
  @param nentid - the length of array entid

  @return The function stores resulting indeces of bnodvalt objects in the parameter bnid.
  @retval 0 - on success
  @retval 1 - cannot create bnodvalt object - parser error

  Created by Tomas Koudelka, 10.2010
*/
long assemble_bnodvalt(selement &el, snode *nodes, long prop, entitybocon &ebc, list &bnvl, gentity ent, 
                       long minid, long *bnid, long *entid, long nentid)
{
  long i, j, l, nbo, nnbo, *propent, **entnodes;
  bnodvalt *tmp;
  long err, rm_ent = 0;
     
  switch (ent)
  {
    case gvertex:
      nbo = el.nne;
      nnbo = 1;
      rm_ent = 1;
      propent = new long[nbo];
      el.searchnodprop(prop, ent, nodes, NULL, NULL, propent);
      entnodes = new long*[nbo];
      for (i=0; i<el.nne; i++)
      {
        entnodes[i] = new long[1];
        entnodes[i][0] = i;
      }
      break;
    case gcurve:
      nbo = el.ned;
      if (nbo == 1)  // 1D element -> take end nodes
      {
        nnbo = 2;
        rm_ent = 1;
        propent = new long[nbo];
        // all nodes of element must lay on the given edge
        propent[0] = prop;
        entnodes    = new long*[nbo];
        entnodes[0] = new long[nnbo];
        // the first and the last nodes are end nodes
        entnodes[0][0] = el.nodes[0];
        entnodes[0][1] = el.nodes[el.nne];
      }
      else  // standard 2D element
      {
        nnbo = el.nned[0];
        propent = el.propedg;
        // entnodes = el.edgenod;
        entnodes = new long*[nbo];
        for (i=0; i<nbo; i++)
          entnodes[i] = new long[nnbo];
        for (i=0; i<nbo; i++)
        {
          for (j=0; j<nnbo; j++)
            entnodes[i][j] = el.edgenod[i][j];
        }
      }
      break;
    case gsurface:
      nbo = el.nsurf;
      nnbo = el.nnsurf[0];
      propent = el.propsurf;
      // entnodes = el.surfnod;
      entnodes = new long*[nbo];
      for (i=0; i<nbo; i++)
        entnodes[i] = new long[nnbo];
      for (i=0; i<nbo; i++)
      {
        for (j=0; j<nnbo; j++)
          entnodes[i][j] = el.surfnod[i][j];
      }
      break;
    default:
      print_err("unknown boundary object is required", __FILE__, __LINE__, __func__);
      abort();
  }
  if ((propent == NULL) && (entid == NULL))
  {
    print_err("boundary properties of elements are not read", __FILE__, __LINE__, __func__);
    abort();
  }

  for (i=0; i<nbo; i++)
    bnid[i] = -1;

  if (propent)
  {
    for (i=0; i<nbo; i++)
    {
      if (propent[i] != prop) // boundary object does not have required property
        continue;
    
      if (ebc.dc == 0) // BC does not depend on coordinates of nodes
      {
        // check for object of bnodvalt with identical nnbo
        for (j=minid; j<bnvl.count(); j++)
        {
          tmp = (bnodvalt *)(bnvl.at(j));
          if (tmp->nsc == nnbo)
          {
            // object with identical nnob has been already created
            bnid[i] = j;
            break;
          }
        }
        if (bnid[i] >= 0)
          continue; // go to next boundary object (i.e. next i)
      }
      if (ebc.dc || bnid[i] < 0)
      // no object of bnodval with identical nnbo has been found
      // or
      // ebc depends on coordinates
      {
        tmp = new bnodvalt;
        tmp->nsc = nnbo;
        tmp->nodval = new gfunct[nnbo];
        for (j=0; j<nnbo; j++)
        {
          err = ebc.getval(nodes[el.nodes[entnodes[i][j]]], tmp->nodval[j]);
          if (err)
            return 1;
        }
        bnvl.append(tmp);
        bnid[i] = bnvl.count()-1;
      }
    }
  }
  
  if (entid)
  {
    // there are entities comming from additional edges or surfaces
    l = 0;
    for (i=0; i<nbo; i++)
    {
      if (i != entid[l]) 
        continue; // boundary object does not have required property
    
      if (ebc.dc == 0) // BC does not depend on coordinates of nodes
      {
        // check for object of bnodvalt with identical nnbo
        for (j=minid; j<bnvl.count(); j++)
        {
          tmp = (bnodvalt *)(bnvl.at(j));
          if (tmp->nsc == nnbo)
          {
            // object with identical nnob has been already created
            bnid[i] = j;
            break;
          }
        }
      }
      if (ebc.dc || bnid[i] < 0)
      // no object of bnodval with identical nnbo has been found
      // or
      // ebc depends on coordinates
      {
        tmp = new bnodvalt;
        tmp->nsc = nnbo;
        tmp->nodval = new gfunct[nnbo];
        for (j=0; j<nnbo; j++)
        {
          err = ebc.getval(nodes[el.nodes[entnodes[i][j]]], tmp->nodval[j]);
          if (err)
            return 1;
        }
        bnvl.append(tmp);
        bnid[i] = bnvl.count()-1;
      }
      if (l < nentid-1)
        l++;
      else
        break;
    }
  }

  if (rm_ent)
    delete [] propent;

  for(i=0; i<nbo; i++)
    delete [] entnodes[i];
  delete [] entnodes;

  return 0;
}



/**
  The function fills array of indeces bnid by the index of the last climatcond object in the list bclimc 
  for the given element and boundary condition (BC).

  @param el      - element to which the BC will be assigned
  @param nodes   - array of all nodes in the given topology
  @param prop    - entity property id of the element el where the entity BC will be applied
  @param bclimc  - list of assigned climatic condition
  @param ent     - entity type for which the the BC will be assigned
  @param bnid    - array with indeces of bnodval objects for the given coefficients (numbering of indeces must be from zero)
                   bnid[i] returns index for i-th boundary object of the given coefficient prescribed by bclimc
  @param entid   - array of entity indeces for the given property number detected on the additional edges or surfaces
  @param nentid  - the length of array entid

  @return The function stores resulting indeces of climatcond objects in the parameter bnid.

  Created by Tomas Koudelka, 09.2011
*/
void assemble_bclimcond(selement &el, snode *nodes, long prop, list &bclimc, gentity ent, 
                        long *bnid, long *entid, long nentid)
{
  long i, l, nbo, *propent;
  long rm_ent = 0;
     
  switch (ent)
  {
    case gvertex:
      nbo = el.nne;
      rm_ent = 1;
      propent = new long[nbo];
      el.searchnodprop(prop, ent, nodes, NULL, NULL, propent);
      break;
    case gcurve:
      nbo = el.ned;
      if (nbo == 1)  // 1D element -> take end nodes
      {
        nbo = 2;
        rm_ent = 1;
        propent = new long[nbo];
        // all nodes of element must lay on the given edge
        propent[0] = prop;
        propent[1] = prop;
      }
      else  // standard 2D element
      {
        propent = el.propedg;
      }
      break;
    case gsurface:
      nbo = el.nsurf;
      propent = el.propsurf;
      break;
    default:
      print_err("unknown boundary object is required", __FILE__, __LINE__, __func__);
      abort();
  }
  if ((propent == NULL) && (entid == NULL))
  {
    print_err("boundary properties of elements are not read", __FILE__, __LINE__, __func__);
    abort();
  }

  for (i=0; i<nbo; i++)
    bnid[i] = -1;

  if (propent)
  {
    for (i=0; i<nbo; i++)
    {
      if (propent[i] != prop) // boundary object does not have required property
        continue;
      else
        bnid[i] = bclimc.count()-1;
    }
  }
  
  if (entid)
  {
    // there are entities comming from additional edges or surfaces
    l = 0;
    for (i=0; i<nbo; i++)
    {
      if (i != entid[l]) 
        continue; // boundary object does not have required property
      else
        bnid[i] = bclimc.count()-1;
    
      if (l < nentid-1)
        l++;
      else
        break;
    }
  }

  if (rm_ent)
    delete [] propent;

  return;
}



/**
  The function creates new object of the loadelt type 
  which contains edge or surface boundary condition (BC) generated for given element and 
  its entities with given property.

  @param eid - element id (numbering must be from zero)
  @param el - element whose BC will be generated
  @param nodes   - array of all nodes in the given topology
  @param prop - entity property id of the element el where the entity BC will be applied
  @param bc - type of assigned BC
  @param bnid - array with indices of bnodval objects for all prescribed coefficients (numbering of indices must be from zero)
                bnid[i][j] returns index for i-th coefficient of j-th boundary object
  @param ent - entity type 
  @param entid - array of entity indeces for the given property number detected on the additional edges or surfaces
  @param nentid - the length of array entid

  @return The function returns pointer to the new allocated object of the loadelt type.

  Created by TKo, 10.2010
*/
loadelt *bc2loadelt(long eid, selement &el, snode *nodes, long prop, bocontypet bc, long **bnid, gentity ent, long *entid, long nentid)
{
  long i, l, nbo, nnbo;
  loadelt *ret = NULL;
  long    *propent, rm_propent;

  if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
    return NULL;

  // array of entity properties from element is used by default
  // -> no delete
  rm_propent = 0;

  switch (ent)
  {
    case gvertex:
      nbo = el.nne;
      nnbo = 1;
      rm_propent = 1;
      propent = new long[nbo];
      el.searchnodprop(prop, ent, nodes, NULL, NULL, propent);
      break;
    case gcurve:
      nbo = el.ned;
      if (nbo == 1)  // 1D element -> take end nodes
      {
        nbo = 2;
        nnbo = 1;
        propent = new long[nbo];
        rm_propent = 1;
        for(i=0; i<nbo; i++)
        {
          // all nodes of element must lay on the given edge
          propent[i] = prop;
        }
      }
      else  // standard 2D element
      {
        nnbo = el.nned[0];
        propent = el.propedg;
      }
      break;
    case gsurface:
      nbo = el.nsurf;
      nnbo = el.nnsurf[0];
      propent = el.propsurf;
      break;
    default:
      return NULL;
  }
  
  if ((propent == NULL) && (entid == NULL)) // no properties were read
    return NULL;

  ret = new loadelt;
  ret->eid = eid; // loadelt uses element numbering from zero
  ret->nbo = nbo;
  ret->nnbo = nnbo;

  //  indicators of boundary conditions
  ret->bc = new bocontypet[nbo];
  memset(ret->bc, 0, sizeof(*ret->bc)*nbo);
  //  id of nodal values or climatic conditions
  ret->nvid = new long[nbo];
  //  id of nodal values describing transmission coefficients
  ret->trcid = new long[nbo];
  //  id of nodal values describing transmission/radiation coefficients
  ret->trrid = new long[nbo];
  
  for(i=0; i<nbo; i++)
    ret->nvid[i] = ret->trcid[i] = ret->trrid[i] = -1;

  if (propent)
  {
    for(i=0; i<nbo; i++)
    {
      if (propent[i] == prop)
      {
        ret->bc[i] = bc;
        ret->nvid[i] = bnid[0][i]; 
        if (bc >= presc_trmiss)
        {
          ret->trcid[i] = bnid[1][i];
          ret->trrid[i] = bnid[2][i];
        }
      }
    }
  }

  if (entid)
  {
    // there are entities comming from additional edges or surfaces
    l = 0;
    for(i=0; i<nbo; i++)
    {
      if (i != entid[l]) 
        continue;  // boundary object does not have required property

      ret->bc[i] = bc;
      ret->nvid[i] = bnid[0][i]; 
      if (bc >= presc_trmiss)
      {
        ret->trcid[i] = bnid[1][i];
        ret->trrid[i] = bnid[2][i];
      }
      if (l < nentid-1)
        l++;
      else
        break;
    }
  }

  if (rm_propent)
    delete [] propent;

  return ret;
}



/**
  The function determines the number of boundary objects of the given 
  entity type on the given element.
  
  @param el - structure of element
  @param ent - entity type

  @return The function returns the number of boundary objects.
  
  Created by TKo, 10.2010
*/
long get_nbo(selement &el, gentity ent)
{
  long nbo=-1;

  switch (ent)
  {
    case gvertex:
      nbo = el.nne;
      break;
    case gcurve:
      nbo = el.ned;
      if (nbo == 1)  // 1D element -> take end nodes
        nbo = 2;
      break;
    case gsurface:
      nbo = el.nsurf;
      break;
    default:
      break;
  }
  return nbo;
}



/**
  The function assembles integration domain description and stores it in the argument idom.
  It is used for the markup of element edges or surfaces, i.e. integration domains, for the computation of 
  flux quantity domain integral.

  @param el    - element to which the BC will be assigned [in]
  @param nodes - array of all nodes in the given topology [in]
  @param prop  - entity property id of the element el where the domain integral will be applied [in]
  @param ent   - entity type for which the the integrated domain will be applied [in]
  @param entid - array of entity indeces for the given property number detected on the additional edges or surfaces [in]
  @param nentid - the length of array entid [in]
  @param idomid - assigned integration domain index [in]
  @param idom   - description of integration domain [out]

  @return The function stores description of integrated domain in idom

  Created by Tomas Koudelka, 02.2018
*/
void assemble_int_domain(selement &el, snode */*nodes*/, long prop, gentity ent, long *entid, long nentid, long idomid, intdomt &idom)
{
  long i, l, nbo, *propent;
  long rm_ent = 0;

  switch (ent)
  {
    case gcurve:
      nbo = el.ned;
      if (nbo == 1)  // 1D element -> take end nodes
      {
        rm_ent = 1;
        propent = new long[nbo];
        // all nodes of element must lay on the given edge
        propent[0] = prop;
      }
      else  // standard 2D element
        propent = el.propedg;
      break;
    case gsurface:
      nbo = el.nsurf;
      propent = el.propsurf;
      break;
    default:
      print_err("unknown boundary object is required", __FILE__, __LINE__, __func__);
      abort();
  }
  if ((propent == NULL) && (entid == NULL))
  {
    print_err("boundary properties of elements are not read", __FILE__, __LINE__, __func__);
    abort();
  }
  idom.idid = new long[nbo];
  idom.n = nbo;
  memset(idom.idid, 0, sizeof(*idom.idid)*nbo);

  if (propent)
  {
    for (i=0; i<nbo; i++)
    {
      if (propent[i] == prop) // boundary object does have the required property
        idom.idid[i] = idomid;
    }
  }

  if (entid)
  {
    // there are entities comming from additional edges or surfaces
    l = 0;
    for (i=0; i<nbo; i++)
    {
      if (i == entid[l]) 
      {
        idom.idid[i] = idomid;
        if (l < nentid-1)
          l++;
        else
          break;
      }
    }
  }

  if (rm_ent)
    delete [] propent;
}



/**
  The function merges new integration domain descriptor with previous ones stored in the given list 
  on the given element. If the merging is not completed fully, the idom is appended inthe descriptor list idoml.
  In such the case, the appended idom has zero domain indeces that was found identical with the indeces of previously stored
  descriptors.

  @param idoml - pointer to the list of previously assigned integration domain descriptors
  @param idom - pointer to the new assigned integration domain descriptor

  @retval 0 - the merging of the new descriptor was merged fully with previous ones, idom was NOT appended in the list 
  @retval 1 - the merging of the new descriptor with the ones in the list could not be done for all boundary objects, 
              idom with some modifications was appended in the list
  @retval 2 - error of merging due to incompatible numbers of boudnary objects in the descriptors

  Created by Tomas Koudelka, 02.2018
*/
long merge_intdom(list *idoml, intdomt *idom)
{
  long i, j, nbo = idom->n;
  long nidoml = idoml->count();
  intdomt *aid;
  long ret = 0;
  long aret;
  
  for (i=0; i<nidoml; i++)
  {
    aret = 0;
    aid = (intdomt *)idoml->at(i);
    if (aid->n == nbo)
    {
      for (j=0; j<nbo; j++)
      {
        if (idom->idid[j] > 0)
        {
          if (aid->idid[j] == idom->idid[j]) // the same integration domain indeces have been detected
            idom->idid[j] = 0;              //  => zero integration domain index in new the integral descriptor       
          else
          {
            if (aid->idid[j] == 0) // there is no active integration domain index in the previous descriptor on the given boundary object
            {
              // copy new integration domain id to the previous descriptor
              aid->idid[j] = idom->idid[j]; 
              idom->idid[j] = 0;
            }
            else  // different integration domain indeces in new and previous descriptors
              aret = 1;
          }
        }
      }
      ret = aret; // set the return value to the actual one
    }
    else // error - incompatible numbers of boundary objects in the merged descriptors
      return 2;
  }

  if (ret)
  {
    // some integration domain indeces could not be merged with the previous ones =>
    // it is necessary to add new descriptor to the list
    idoml->append(idom); 
  }

  return ret;
}
