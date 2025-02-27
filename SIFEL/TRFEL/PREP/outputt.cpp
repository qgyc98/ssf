#include "globprept.h"
#include "outputt.h"
#include "elementt.h"
#include "intdomt.h"
#include <string.h>



/**
  This function writes read data from the preprocessor files to the file given by
  the out parameter in the MEFEL input format.

  @param out - pointer to the opened text file where the data will be written
  @param d   - pointer to the description structure with preprocessor data

  @retval 0 - on success
  @retval 7 - in case of an error

  Created by Tomas Koudelka, 10.2010
*/
long outputt(FILE *out, descript *d)
{
  // Output data in mechtop format
  long err;

  fprintf(stdout, "\nProblem model output :\n");
  Dbmatt->renumber_id();
  Dbcrst->renumber_id();
  fprintf(stdout, " writing of nodes . . .");
  err = wr_nodest(out); // writing of nodes, crsecs in the nodes and local systems
  fflush(out);
  if (err)
    return(7);
  fprintf(stdout, " O.K.\n");

  fprintf(stdout, " writing of boundary conditions . . .");
  err = wr_bocont(out); // writing of boundary conditions at nodes
  fflush(out);
  if (err)
    return(7);
  fprintf(stdout, " O.K.\n");

  fprintf(stdout, " writing of elements . . .");
  err = wr_elementst(out, d); // writing of elements and code numbers on elements
  fflush(out);
  if (err)
    return(7);
  fprintf(stdout, " O.K.\n");

  if (d->paral)
  {
    fprintf(stdout, " writing of global node numbers . . .");
    err = wr_globnodnumt(out);
    fflush(out);
    if (err)
      return(7);
    fprintf(stdout, " O.K.\n");
  }

  fprintf(stdout, " writing of materials . . .");
  err = wr_materialst(out, d); // writing of material parameters
  fflush(out);
  if (err)
    return(7);
  fprintf(stdout, " O.K.\n");

  fprintf(stdout, " writing of cross-sections . . .");
  err = wr_crsecst(out, d);    // writing of cross-section parameters
  fflush(out);
  if (err)
    return(7);
  fprintf(stdout, " O.K.\n");

  fprintf(stdout, " writing of load cases . . .");
  err = wr_loadcaset(out); // writing of load cases
  fflush(out);
  if (err)
    return(7);
  fprintf(stdout, " O.K.\n");

  fprintf(stdout, " writing of flux integration domains . . .");
  err =  wr_fluxint(out);
  fflush(out);
  if (err)
    return(7);
  fprintf(stdout, " O.K.\n");

  fprintf(stdout, " writing of initial conditions . . .");
  err = wr_initcondt(out); // writing of initial conditions
  fflush(out);
  if (err)
    return(7);
  fprintf(stdout, " O.K.\n");

  fprintf(stdout, " writing of advection velocities . . .");
  err = wr_advect(out); // writing of advection velocities
  fflush(out);
  if (err)
    return(7);
  fprintf(stdout, " O.K.\n");

  long i, j, k;
  bocont   *bc;
  sourcet  *src;
  bnodvalt *bnv;
  intdomt  *idom;
  climatcond *cc;
  climatcond2 *cc2;
  char *ccf;
  answertype *trcc;
  long ndof = Tp->ntm;

  // deallocating used memory for nodes
  for (i=0; i<Top->nn; i++)
  {
    if (Nod_hang)
      delete Nod_hang[i];
    delete [] Nod_bocon[i];
    delete [] Nod_ccn[i];
    delete [] Nod_inicd[i];
    if (Nod_advect)
      delete [] Nod_advect[i];
    if (Nod_sourcet)
      delete [] Nod_sourcet[i];
  }
  delete [] Nod_hang;  
  delete [] Nod_bocon; 
  delete [] Nod_ccn;   
  delete [] Nod_cst;   
  delete [] Nod_csti;  
  delete [] Nod_cstdbi;  
  delete [] Nod_inicd; 
  if (Nod_sourcet)
    delete [] Nod_sourcet;
  if (Nod_advect)
  {
    delete [] Nod_advect;
    delete [] Nod_advect_nc;
  }

  for(i=0; i < Nod_bclst.count(); i++)
  {
    bc = (bocont *)(Nod_bclst.at(i));
    delete bc;
  }
  for(i=0; i < Nod_srclst.count(); i++)
  {
    src = (sourcet *)(Nod_srclst.at(i));
    delete src;
  }
  
  // deallocating used memory for elements
  for (i=0; i<Top->ne; i++)
  {
    delete [] El_mattype[i];
    delete [] El_matid[i];
    delete [] El_matdbi[i];
    if (El_sourcet)
      delete [] El_sourcet[i];   
  }
  delete [] El_type;   
  delete [] El_nmat;   
  delete [] El_mattype;
  delete [] El_matid;  
  delete [] El_matdbi;  
  delete [] El_cst;    
  delete [] El_csti;   
  delete [] El_cstdbi;   
  delete [] El_sourcet;

  for (i=0; i<ndof; i++)
  {
    for(j=0; j<Top->ne; j++)
    {
      delete El_loadt[i][j];
      if (El_fluxint[i][j])
      {
        for (k=0; k<El_fluxint[i][j]->count(); k++)
        {
          idom = (intdomt *)El_fluxint[i][j]->at(k);
          delete idom;
        }
      }
      delete El_fluxint[i][j];
    }

    delete [] El_loadt[i];
    delete [] El_fluxint[i];
  }

  if (El_loadt)
  {
    delete [] El_loadt;   

    for (i=0; i<ndof; i++)
    {
      for(j=0; j < El_nv_lst[i].count(); j++)
      {
        bnv = (bnodvalt *)(El_nv_lst[i].at(j));
        delete bnv;
      }

      for(j=0; j < El_trc_lst[i].count(); j++)
      {
        bnv = (bnodvalt *)(El_trc_lst[i].at(j));
        delete bnv;
      }

      for(j=0; j < El_trr_lst[i].count(); j++)
      {
        bnv = (bnodvalt *)(El_trr_lst[i].at(j));
        delete bnv;
      }

      for(j=0; j < El_trcc_lst[i].count(); j++)
      {
        trcc = (answertype *)(El_trcc_lst[i].at(j));
        delete trcc;
      }
      for(j=0; j < El_cc_lst[i].count(); j++)
      {
        cc = (climatcond *)(El_cc_lst[i].at(j));
        delete cc;
      }
      for(j=0; j < El_ccf_lst[i].count(); j++)
      {
        ccf = (char *)(El_ccf_lst[i].at(j));
        delete ccf;
      }
      for(j=0; j < El_gtrcc_lst[i].count(); j++)
      {
        trcc = (answertype *)(El_gtrcc_lst[i].at(j));
        delete trcc;
      }
      for(j=0; j < El_gcc_lst[i].count(); j++)
      {
        cc2 = (climatcond2 *)(El_gcc_lst[i].at(j));
        delete cc2;
      }
      for(j=0; j < El_gccf_lst[i].count(); j++)
      {
        ccf = (char *)(El_gccf_lst[i].at(j));
        delete ccf;
      }
    }
    delete [] El_nv_lst;
    delete [] El_trc_lst;
    delete [] El_trr_lst;
    delete [] El_trcc_lst;
    delete [] El_cc_lst;
    delete [] El_ccf_lst;
    delete [] El_gtrcc_lst;
    delete [] El_gcc_lst;
    delete [] El_gccf_lst;
  }

  delete [] El_fluxint;   

  delete [] El_tfunc;  
  

  for(i=0; i < El_srclst.count(); i++)
  {
    src = (sourcet *)(El_srclst.at(i));
    delete src;
  }
  
  // deallocating rest of objects
  if (Dbmatt)   delete Dbmatt;
  if (Dbcrst)   delete Dbcrst;
  if (Tft)      delete [] Tft;

  return (0);
}



/**
  Function writes section with description of nodes to the text file given by out.

  @param out - pointer to the opened text file, where the data will be written
  @param nprop - number of entries in the input property file

  @retval 0 - if succes
  @retval 1 - if node hasn't assigned ndof
  @retval 2 - if node hasn't assigned cross-section

  Created by Tomas Koudelka, 10.2010
*/
long wr_nodest(FILE *out)
{
  long i, ndofn;

  ndofn = Tp->ntm;
  fprintf(out, "# sections of nodes\n");

  fprintf(out, "%ld\n", Top->nn);
  for (i = 0; i < Top->nn; i++)
  {
    // writing node number and coordinates
    fprintf(out, "%6ld % .10e % .10e % .10e", i+1, Top->nodes[i].x, Top->nodes[i].y, Top->nodes[i].z);
    if (Numhn > 0)
    {
      if (Nod_hang[i])
        Nod_hang[i]->print(out);
      else
        fprintf(out, " %ld", ndofn);
    }
    else
      fprintf(out, " %ld", ndofn);

    fprintf(out, " %d", Nod_cst[i]);

    if (Nod_cst[i])
      fprintf(out, " %ld\n", Dbcrst->crs[Nod_cstdbi[i]].ridx[Nod_csti[i]-1]);
    else
      fprintf(out, "\n");
  }
  fprintf(out, "\n");
  return(0);
}



/**
  Function writes section with description of boundary conditions of nodes
  to the text file given by out.

  @param out - pointer to the opened text file, where the data will be written

  @retval 0 - on succes
  @retval 1 - node has assigned boundary condition and coupled dof at the same direction

  Created by Tomas Koudelka, 10.2010
*/
long wr_bocont(FILE *out)
{
  long i, k, nbc, ndof;
  char emsg[200];

  nbc = 0;
  ndof = Tp->ntm;
  
  fprintf(out, "# boundary conditions at nodes\n");
  for (i = 0; i < Top->nn; i++)
  {
    if (Nod_bocon[i])
    // test if any boundary condition is prescribed
      nbc++;
    if (Nod_ccn[i])
    // test if any common code numbers are prescribed
      nbc++;

  }
  fprintf(out, "%ld\n", nbc);
  for (i = 0; i < Top->nn; i++)
  {
    if ((Nod_bocon[i]) || (Nod_ccn[i]))
    // test if any boundary condition is prescribed or any common code numbers are prescribed
    {
      fprintf(out, "%ld", i+1);
      for (k = 0; k < ndof; k++)
      { // loop for writing bc
        // if load case has not bc, write 1
        if (Nod_bocon[i] && Nod_ccn[i])
        {
          if ((Nod_bocon[i][k]) && (Nod_ccn[i][k]))
          {
            sprintf(emsg, "node %ld has prescribed boundary condition and common code number\n"                
                          " for the same loadcase\n.", i+1);
            print_err(emsg, __FILE__, __LINE__, __func__);
            return(1);
          }
          if ((Nod_bocon[i][k] == 0) && (Nod_ccn[i][k] == 0))
          // no boundary condition nor common code number
          {
            fprintf(out, " 1");
            continue;
          }
        }
        if (Nod_bocon[i])
        {
          if (Nod_bocon[i][k] == 0)
          // no boundary condition
          {
            fprintf(out, " 1");
            continue;
          }
          else
          // prescribed value
          {
            fprintf(out, " %ld", -Nod_bocon[i][k]);
            continue;
          }
        }
        if (Nod_ccn[i])
        // common code number
        {
          if (Nod_ccn[i][k])
          {
            fprintf(out, " %ld", Nod_ccn[i][k]);
            continue;
          }
          else
          {
            if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
              fprintf(out, " 0");
            else
              fprintf(out, " 1");
            continue;
          }
        }
      }
      fprintf(out, "\n");
    }
  }
  fprintf(out, "\n");
  return(0);
}



/**
  Function writes section with description of elements to the text file given by out.

  @param out - pointer to the opened text file, where the data will be written
  @param d   - pointer to the description structure with preprocessor data

  @retval 0 : on succes
  @retval 1 : if element has not assigned a cross-section
  @retval 2 : if element has not assigned a time function

  Created by Tomas Koudelka, 10.2010
*/
long wr_elementst(FILE *out, descript */*d*/)
{
  long i, j, k;
  char emsg[200];
  elemtypet et;

  
  fprintf(out, "# section of elements\n");
  fprintf(out, "%ld\n", Top->ne);
  for (i = 0; i < Top->ne; i++)
  {
    if (((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
        && (El_tfunc[i] == 0))
    {
      sprintf(emsg, "element number %ld has not assigned a time function", i+1);
      print_err(emsg, __FILE__, __LINE__, __func__);
      return(2);
    }

    // writing type of element, element number
    et = El_type[i];
    fprintf(out, "%6ld %3d", i+1, et);

    // writing element nodes
    for (j = 0; j < Top->elements[i].nne; j++)
    {
/*      if ((et == planeelementqt) && (d->t3d) && (j == 3))
      // nodes in T3d file has in this case folowing order 1 2 3 5 6 4 in view of
      // SIFEL order 1 2 3 4 5 6
        continue;*/
      fprintf(out, " %6ld", Top->elements[i].nodes[j]+1);
    }
/*    if ((et == planeelementqt) && (d->t3d))
    // nodes in T3d file has in this case folowing order 1 2 3 5 6 4 in view of
    // SIFEL order 1 2 3 4 5 6
    // so it is need to write fourth node
      fprintf(out, " %6ld", Top->elements[i].nodes[3]+1);*/

    // no code numbers on element (code numbers not supported at this time)
    fprintf(out, " 0");
  
    // time function for growing structure problem
    if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
      fprintf(out, " %ld", El_tfunc[i]);
      
    fprintf(out, " %d", El_cst[i]);
    if (El_cst[i] > 0)
      fprintf(out, " %ld", Dbcrst->crs[El_cstdbi[i]].ridx[El_csti[i]-1]); // if element has cs, write cs index

    for (k = 0; k < El_nmat[i]; k++)
      fprintf(out, " %d %ld", El_mattype[i][k], Dbmatt->mat[El_matdbi[i][k]].ridx[El_matid[i][k]-1]); // writing element material type and index
    fprintf(out, "\n");
  }
  fprintf(out, "\n");
  return(0);
}



/**
  Function writes section with description of elements to the text file given by out.

  @param out - pointer to the opened text file, where the data will be written

  @retval 0 : on succes

  Created by Tomas Koudelka, 10.2010
*/
long wr_globnodnumt(FILE *out)
{
  fprintf(out, "# definition of global node numbers\n");
  for (long i = 0; i < Top->nn; i++)
    fprintf(out, "%ld\n", Top->gnn[i]);
  return(0);
}



/**
 Function writes materials and thier parameters to the file given by the out.

  @param out - pointer to the opened text file, where the data will be written
  @param d   - pointer to the description structure with preprocessor data

  @retval 0 : on succes

  Created by Tomas Koudelka, 10.2010
*/
long wr_materialst(FILE *out, descript *d)
{
  long i, j, wrt;

  fprintf(out, "# definition of materials\n");
  fprintf(out, "%ld\n", Dbmatt->nmatu);  // writing number of material types
  for (i = 0; i < Dbmatt->numt; i++)
  {
    wrt = 1;
    for (j = 0; j < Dbmatt->mat[i].ninst; j++)
    {
      if (Dbmatt->mat[i].instu[j])  // material type index is used (<> 0)
      {
        if (wrt)  // for the first time write material type number and number of used indeces
        {
          fprintf(out, "%d %ld\n", Dbmatt->mat[i].type, Dbmatt->mat[i].ninstu);
          wrt = 0;  // next time don't write
        }
        if (d->matstr == yes)
          // write material type index and property string
          fprintf(out, "%ld %s\n", Dbmatt->mat[i].ridx[j], Dbmatt->mat[i].inst[j]);
        else
        {
          // write material type index and parameters from mechmat
          fprintf(out, "%ld ", Dbmatt->mat[i].ridx[j]);
          Tm->printmatchar(out, Dbmatt->mat[i].type, j);
          fprintf(out, "\n");
        }
      }
    }
  }
  fprintf(out, "\n");
  return(0);
}



/**
 Function writes cross-sections and thier parameters to the file given by the out.

  @param out - pointer to the opened text file, where the data will be written
  @param d   - pointer to the description structure with preprocessor data

  @retval 0 : on succes

  Created by Tomas Koudelka, 10.2010
*/
long wr_crsecst(FILE *out, descript *d)
{
  long i, j, wrt;

  fprintf(out, "# definition of cross-sections\n");
  fprintf(out, "%ld\n", Dbcrst->ncrsu);    // writing number of cs types
  for (i = 0; i < Dbcrst->numt; i++)
  {
    wrt = 1;
    for (j = 0; j < Dbcrst->crs[i].ninst; j++)
    {
      if (Dbcrst->crs[i].instu[j])  // cs type index is used (<> 0)
      {
        if (wrt) // for the first time write cs type number and number of used indeces
        {
          fprintf(out, "%d %ld\n", Dbcrst->crs[i].type, Dbcrst->crs[i].ninstu);
          wrt = 0;  // next time don't write
        }
        if (d->crsstr == yes)
          // write cs type index and property string
          fprintf(out, "%ld %s\n", Dbcrst->crs[i].ridx[j], Dbcrst->crs[i].inst[j]);
        else
        {
          // write cs type index and parameters from mechcrsec
          fprintf(out, "%ld ", Dbcrst->crs[i].ridx[j]);
          Tc->printcrschar(out, Dbcrst->crs[i].type, j);
          fprintf(out, "\n");
        }
      }
    }
  }
  fprintf(out, "\n");
  return(0);
}



/**
 Function writes section with description of load to the text file given by out

  @param out - pointer to the opened text file, where the data will be written

  @retval 0 : on success
  @retval 1 : error in the printing of climatcond objects

  Created by TKo, 10.2010
*/
long wr_loadcaset(FILE *out)
{
  long i, err;

  fprintf(out, "# section with load ceses\n");
  fprintf(out, "%ld\n\n", Tp->ntm);

  for (i = 0; i < Tp->ntm; i++)
  {
    fprintf(out, "#\n# loadcase number %ld\n#\n", i+1);
    switch (Tp->tprob)
    {
      case stationary_problem:
      case nonlinear_stationary_problem:
        wr_prescquant(out); // writing of prescribed unknowns (Dirichlets BC) [+ initial conditions for nonstationary problems] 
        wr_sources(out, i); // writing of nodal and element sources
        wr_loadelt(out, i); // writing of element BC
        fflush(out);
        break;
      case nonstationary_problem:
      case growing_np_problem:
      case nonlinear_nonstationary_problem:
      case growing_np_problem_nonlin:
        wr_prescquant(out); // writing of prescribed unknowns (Dirichlets BC) [+ initial conditions for nonstationary problems] 
        wr_sources(out, i); // writing of nodal and element sources
        wr_loadelt(out, i); // writing of element BC
        fflush(out);
        break;
      case discont_nonstat_problem:
      case discont_nonlin_nonstat_problem:
      case hermes:
        wr_prescquant(out); // writing of prescribed unknowns (Dirichlets BC) [+ initial conditions for nonstationary problems] 
        wr_sources(out, i); // writing of nodal and element sources
        wr_loadelt(out, i); // writing of element BC
        fflush(out);
        break;
      default:
        print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
    }
    err = wr_climatcond(out, i);
    if (err)
      return 1;
    err = wr_climatcond2(out, i);
    if (err)
      return 1;
    //  the number of scaling time functions
    fprintf (out,"# scaling time functions\n");
    fprintf (out,"0 \n");
  }
  
  return(0);
}



/**
  Function writes section with description of nodal prescribed
  quantities to the text file given by out.

  @param out   - pointer to the opened text file, where the data will be written

  @retval 0 : on succes
  @retval 1 : wrong number of prescribed displacements (error in code)

  Created by TKo, 10.2010
*/
long wr_prescquant(FILE *out)
{
  long i;
  long npq = Nod_bclst.count();
  bocont *bc;
  
  fprintf(out, "# prescribed quantities\n");
  fprintf(out, "%ld\n", npq);
  if (npq == 0)
    return(0);

  for (i = 0; i < npq; i++)
  {
    bc = (bocont*)(Nod_bclst.at(i));
    bc->print(out);
  }
  fprintf(out, "\n");
  return(0);
}



/**
 Function writes section with description of nodal load to the text file given by out

  @param out   - pointer to the opened text file, where the data will be written
  @param nlc   - number of loading case which will be written

  @retval 0 : on succes

  Created by Tomas Koudelka, 10.2010
*/
long wr_sources(FILE *out, long nlc)
{
  long i, j, nsrc;
  sourcet *aux;
  long nns = Nod_srclst.count(); 
  long nes = El_srclst.count();
  long *nsl, *esl;
  long nuns, nues;

  nsl = new long[nns];
  memset(nsl, 0, sizeof(*nsl)*nns);
  esl = new long[nes];
  memset(esl, 0, sizeof(*esl)*nes);

  // mark sources at nodes used in the given load case nlc
  nuns = 0;
  if (Nod_sourcet)
  {
    for (i = 0; i < Top->nn; i++)
    {
      if (Nod_sourcet[i])
      {
	if (Nod_sourcet[i][nlc])
	{
	  nuns++;
	  nsl[Nod_sourcet[i][nlc]-1] = 1;
	}
      }
    }
  }

  // mark sources on elements used in the given load case nlc
  nues = 0;
  if (El_sourcet)
  {
    for (i = 0; i < Top->ne; i++)
    {
      if (El_sourcet[i])
      {
	if (El_sourcet[i][nlc])
	{
	  nues++;
	  esl[El_sourcet[i][nlc]-1] = 1;
	}
      }
    }
  }

  // total number of sources used in the nlc load case
  nsrc = 0L;
  for (i=0; i<nns; i++)
  {
    if (nsl[i])   nsrc++;
  }
  for (i=0; i<nes; i++)
  {
    if (esl[i])    nsrc++;
  }
  
  fprintf(out, "# list of sources\n");
  fprintf(out, "%ld\n\n", nsrc);  // writing overall number of nodal and element sources used in the nlc
  if (nsrc)
  {
    fprintf(out, "# list of sources for nodes\n");
    j = 0L;
    for (i = 0; i < Nod_srclst.count(); i++)
    {
      aux = (sourcet*)(Nod_srclst.at(i));
      // print nodal sources used in the given load case nlc
      if (nsl[i])
      {
        j++;
        fprintf(out, "\n\n%ld # nodal source object number %ld\n", j, j);
	aux->print(out);
      }
    }
    fprintf(out, "# list of sources for elements\n");
    for (i = 0; i < El_srclst.count(); i++)
    {
      aux = (sourcet*)(El_srclst.at(i));
      // print element sources used in the given load case nlc
      if (esl[i])
      {
        j++;
        fprintf(out, "\n\n%ld # element source object number %ld\n", j, j);
	aux->print(out);
      }
    }

    fprintf(out, "\n# sources at nodes\n");
    fprintf(out, "%ld\n", nuns);
    if (nuns)
    {
      for (i = 0; i < Top->nn; i++)
      {
	if (Nod_sourcet[i])
	  fprintf(out, "%ld %ld\n", i+1, Nod_sourcet[i][nlc]);
      }
    }

    fprintf(out, "# sources at elements\n");
    fprintf(out, "%ld\n", nues);
    if (nues)
    {
      for (i = 0; i < Top->ne; i++)
      {
	if (El_sourcet[i])
	  fprintf(out, "%ld %ld\n", i+1, El_sourcet[i][nlc]);
      }
    }
  }

  fprintf(out, "# number of nodes with point sources\n0\n"); // not yet implemented

  fprintf(out, "\n");

  delete [] nsl;
  delete [] esl;

  return(0);
}



/**
  Function writes section with description of element load to the text file given by out

  @param out   - pointer to the opened text file, where the data will be written
  @param nlc   - number of loading case which will be written

  @retval 0 : on succes

  Created by Tomas Koudelka, 10.2010
*/
long wr_loadelt(FILE *out, long nlc)
{
  long nnv, ntrc, ntrr;
  bnodvalt *tmp;
  fprintf(out, "# load of elements\n");
  long i, nle;

  nnv  = El_nv_lst[nlc].count();
  ntrc = El_trc_lst[nlc].count();
  ntrr = El_trr_lst[nlc].count();

  if (El_loadt[nlc])
  {
    for (i = 0, nle = 0; i < Top->ne; i++)
    {  // counting overall number of loaded elements
      if(El_loadt[nlc][i])
      {
        nle++;
        // stored indices of bnodval objects has to be renumbered
        // due to separate storage of bnodval objects for
        // nodal values, transmission coefficients and radiation coefficient
        El_loadt[nlc][i]->renumber_id(nnv, ntrc);
      }
    }
  }

  // writing overall number of element BC entries
  fprintf(out, "%ld\n", nle);
  if (El_loadt[nlc])
  {
    for (i = 0; i < Top->ne; i++)
    {
      if (El_loadt[nlc][i])
      {
        // writing of load
        El_loadt[nlc][i]->print(out, nlc);
      }
    }
  }

  fprintf(out, "\n");

  fprintf(out, "# list of nodal values\n");
  fprintf(out, "%ld\n", nnv+ntrc+ntrr);
  fprintf(out, "# list of general nodal values (%ld obj.)\n", El_nv_lst[nlc].count());
  for (i=0; i<El_nv_lst[nlc].count(); i++)
  {
    tmp = (bnodvalt *)(El_nv_lst[nlc].at(i));
    tmp->print(out);
  }

  fprintf(out, "# list of nodal values for transmission coefficients (%ld obj.)\n", El_trc_lst[nlc].count());
  for (i=0; i<El_trc_lst[nlc].count(); i++)
  {
    tmp = (bnodvalt *)(El_trc_lst[nlc].at(i));
    tmp->print(out);
  }
  
  fprintf(out, "# list of nodal values for radiation coefficients (%ld obj.)\n", El_trr_lst[nlc].count());
  for (i=0; i<El_trr_lst[nlc].count(); i++)
  {
    tmp = (bnodvalt *)(El_trr_lst[nlc].at(i));
    tmp->print(out);
  }

  return(0);
}



/**
  Function writes section with description of flux integration domains on elements to the text file given by out

  @param out   - pointer to the opened text file, where the data will be written

  @retval 0 : on succes

  Created by Tomas Koudelka, 03.2018
*/
long wr_fluxint(FILE *out)
{
  long i, j, k, lcid, ndof, nidd;
  intdomt *idom;

  fprintf(out, "\n# section with the list of flux integration domain descriptors\n");

  ndof = Tp->ntm;

  for(lcid=0; lcid<ndof; lcid++)
  {
    if (El_fluxint[lcid])
    {
      for (j=0, nidd=0; j<Top->ne; j++)
      {
        // add the number of descriptors for domain indeces of flux quantity integrals
        if (El_fluxint[lcid][j])
            nidd += El_fluxint[lcid][j]->count();
      }
    }
    else{
      nidd=0;
      continue;
    }

    fprintf(out, "%ld   # number of elements with the integration domain descriptors in load case %ld\n", nidd, lcid+1);
    if (nidd)
    {
      for (i = 0; i < Top->ne; i++)
      {
        if (El_fluxint[lcid][i])
        {
          for(j=0; j<El_fluxint[lcid][i]->count(); j++) // one element may contain several different integration descriptors
          {
            fprintf (out,"\n %ld ",i+1);
            idom = (intdomt *)El_fluxint[lcid][i]->at(j);
            for (k=0; k<idom->n; k++)
            {
              if (idom->idid[k])
                fprintf(out, " 2 %ld ", idom->idid[k]); // BC type = 2, integration domain index
              else
                fprintf(out, " 0 ");  // BC type = 0 => no BC
            }
          }
        }
      }
    }
    fprintf(out, "\n");
  }

  fprintf(out, "\n");
  return 0;
}



/**
  The function writes data of climatic boundary conditions to the text file given by out

  @param out   - pointer to the opened text file, where the data will be written
  @param nlc   - number of loading case which will be written

  @retval 0 : on succes
  @retval 1 : if the mixed types of climatic conditions appeared
  @retval 2 : if no climatic condition has been found on some position in the lists

  Created by Tomas Koudelka, 09.2011
*/
long wr_climatcond(FILE *out, long nlc)
{
  long i, ncc;
  climatcond *cc;
  char *ccf;
  answertype trcc, aux;
  long trcc_id;

  //  the number of climatic conditions
  ncc = El_trcc_lst[nlc].count();
  fprintf (out,"# climatic conditions (%ld obj.)\n", ncc);
  fprintf (out, "%ld", ncc); // all climatcond objects have trcc

  if (ncc == 0)
  {
    fprintf(out, "\n");
    return 0;
  }

  // flag for reading of climatcond from one file or from the input file and form the file 
  // with measurements
  trcc = *((answertype *)El_trcc_lst[nlc].at(0)); 
  if (trcc == no)
  {
    trcc_id = 1;
    fprintf (out, " %ld\n", trcc_id); 
  }
  else
  {
    trcc_id = 2;
    fprintf (out, " %ld\n", trcc_id);
  }
  for(i=0; i<El_trcc_lst[nlc].count(); i++)
  {
    aux = *((answertype *)El_trcc_lst[nlc].at(i));
    if (trcc != aux) // actually, TRFEL does not support mixed 
                     // types of climatic conditions
    {
      print_err("type %ld of climatcond object %ld does not correspond to the original one %ld\n"
                "(TRFEL cannot still handle mixed types)", __FILE__, __LINE__, __func__, trcc_id, i+1, aux);
      return 1;
    }
    fprintf(out, "\n%ld", i+1);  // climatcond object id
    cc = (climatcond *)El_cc_lst[nlc].at(i);
    ccf = (char*)El_ccf_lst[nlc].at(i);
    if ((cc == NULL) && (ccf == NULL)) // this is logical error in the source code of inputt.cpp 
                                       // it should not happen theoretically
    {
      print_err("no climatcond object nor file name for climatcond was specified\n"
                "for the climatcond id %ld (error is in the inputt.cpp source code)", __FILE__, __LINE__, __func__, i+1);
      return 2;
    }
    if (cc)
    {
      fprintf(out, "\n");
      cc->print(out);
    }
    else
      fprintf(out, " %s\n", ccf);
  }
  return 0;
}



/**
  The function writes data of general climatic boundary conditions to the text file given by out

  @param out   - pointer to the opened text file, where the data will be written
  @param nlc   - number of loading case which will be written

  @retval 0 : on succes
  @retval 1 : if the mixed types of climatic conditions appeared
  @retval 2 : if no climatic condition has been found on some position in the lists

  Created by Tomas Koudelka, 07.2014
*/
long wr_climatcond2(FILE *out, long nlc)
{
  long i, ncc;
  climatcond2 *cc;
  char *ccf;
  answertype trcc, aux;
  long trcc_id;

  //  the number of climatic conditions
  ncc = El_gtrcc_lst[nlc].count();
  fprintf (out,"# general climatic conditions (%ld obj.)\n", ncc);
  fprintf (out, "%ld", ncc); // all climatcond objects have trcc

  if (ncc == 0)
  {
    fprintf(out, "\n");
    return 0;
  }

  // flag for reading of climatcond2 from one file or from the input file and form the file 
  // with measurements
  trcc = *((answertype *)El_gtrcc_lst[nlc].at(0)); 
  if (trcc == no)
  {
    trcc_id = 1;
    fprintf (out, " %ld\n", trcc_id); 
  }
  else
  {
    trcc_id = 2;
    fprintf (out, " %ld\n", trcc_id);
  }
  for(i=0; i<El_gtrcc_lst[nlc].count(); i++)
  {
    aux = *((answertype *)El_gtrcc_lst[nlc].at(i));
    if (trcc != aux) // actually, TRFEL does not support mixed 
                     // types of climatic conditions
    {
      print_err("type %ld of climatcond2 object %ld does not correspond to the original one %ld\n"
                "(TRFEL cannot still handle mixed types)", __FILE__, __LINE__, __func__, trcc_id, i+1, aux);
      return 1;
    }
    fprintf(out, "\n%ld", i+1);  // climatcond2 object id
    cc = (climatcond2 *)El_gcc_lst[nlc].at(i);
    ccf = (char*)El_gccf_lst[nlc].at(i);
    if ((cc == NULL) && (ccf == NULL)) // this is logical error in the source code of inputt.cpp 
                                       // it should not happen theoretically
    {
      print_err("no climatcond2 object nor file name for climatcond was specified\n"
                "for the climatcond2 id %ld (error is in the inputt.cpp source code)", __FILE__, __LINE__, __func__, i+1);
      return 2;
    }
    if (cc)
      cc->print(out);
    else
      fprintf(out, "%s\n", ccf);
  }
  return 0;
}



/**
  Function writes section with description of initial conditions for nonstationary problems
  to the text file given by out.

  @param out   - pointer to the opened text file, where the data will be written
  
  @retval 0 : on succes

  Created by Tomas Koudelka, 10.2010
*/
long wr_initcondt(FILE *out)
{
  long i, j;

  fprintf(out, "# initial conditions\n");
  if (Nod_inicd[0])
  {
    for (i = 0; i < Top->nn; i++)
    {
      // writing node number
      fprintf(out, "%6ld", i+1);
      if (Numhn > 0) 
      {
        if (Nod_hang[i]) // there is a hanging node => do not print values 
        {                // because they are defined on master nodes
          fprintf(out, "\n");
          continue;
        }
      }
      // print initial conditions for regular nodes only
      for (j=0; j<Tp->ntm; j++) // loop over number of degrees of freedom
        fprintf(out, " %le", Nod_inicd[i][j]);
      fprintf(out, "\n");
    }
    fprintf(out, "\n");
  }
  return(0);
}



/**
  Function writes section with advection velocities at nodes to the 
  text file given by out.

  @param out   - pointer to the opened text file, where the data will be written
  
  @retval 0 : on succes

  Created by Tomas Koudelka, 5.2016
*/
long wr_advect(FILE *out)
{
  long i, j;

  if (Nod_advect)
  {
    fprintf(out, "# vectors of advection velocities at nodes\n");
    for (i = 0; i < Top->nn; i++)
    {
      // writing node number
      fprintf(out, "%6ld", i+1);
      for (j=0; j<Nod_advect_nc[i]; j++) // loop over number of velocity components
	fprintf(out, " %le", Nod_advect[i][j]);
      fprintf(out, "\n");
    }
    fprintf(out, "\n");
  }
  return(0);
}

