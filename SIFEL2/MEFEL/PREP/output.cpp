#include "output.h"
#include "globprep.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "element.h"
#include "ipmap.h"



/**
  This function writes read data from the preprocessor files to the file given by
  the out parameter in the MEFEL input format.

  @param out - pointer to the opened text file where the data will be written
  @param d   - pointer to the description structure with preprocessor data

  Created  by Tomas Koudelka 2009
*/
long output(FILE *out, descrip *d)
{
  // Output data in mechtop format
  long err;

  fprintf(stdout, "\nProblem model output :");
  Dbmat->renumber_id();
  Dbcrs->renumber_id();
  if (Nperbc){
    fprintf(stdout, "\n writing periodic boundary conditions to log file . . .");
    err = wr_periodbc(Log);
    if (err)
      return(7);
    fprintf(stdout, " O.K.");
    
    fprintf(stdout, "\n transformation of periodic boundary conditions to hanging nodes . . .");
    transf_periodbc_hangnod();
    fprintf(stdout, " O.K.");
  }
  
  fprintf(stdout, "\n writing nodes . . .");
  if (Nlay == 0)
    err = wr_nodes(out); // writing nodes, crsecs in the nodes and local systems
  else
    err = wr_lay_nodes(out); // writing nodes, crsecs in the nodes and local systems for layered problems
  fflush(out);
  if (err)
    return(7);
  fprintf(stdout, " O.K.");

  fprintf(stdout, "\n writing boundary conditions . . .");
  err = wr_bocon(out); // writing boundary conditions at nodes
  fflush(out);
  if (err)
    return(7);
  fprintf(stdout, " O.K.");

  fprintf(stdout, "\n writing elements . . .");
  if (Nlay == 0)
    err = wr_elements(out, d); // writing elements and code numbers on elements
  else
    err = wr_lay_elements(out, d); // writing elements and code numbers on elements for layered problems
  fflush(out);
  if (err)
    return(7);
  fprintf(stdout, " O.K.");

  if (d->paral)
  {
    fprintf(stdout, "\n writing global node numbers . . .");
    err = wr_globnodnum(out,d);
    if (err)
      return(7);
    fprintf(stdout, " O.K.");
  }

//  err = wr_intpoints(out, nprope); // writing integration points of elements
                                   // with indeces of material parameters
  fflush(out);
  if (err)
    return(7);

  fprintf(stdout, "\n writing materials . . .");
  err = wr_materials(out, d); // writing material parameters
  fflush(out);
  if (err)
    return(7);
  fprintf(stdout, " O.K.");

/*
  if (Mp->stresscomp)
  {
    fprintf(stdout, "\n writing auxiliary points for stress computation . . .");
    err = wr_auxpoint(out, nprope); // writing auxiliary points
    fflush(out);
    if (err)
      return(7);
    fprintf(stdout, " O.K.");
    fprintf(stdout, "\n writing local coordinate systems . . .");
    err = wr_ellcsys(out, nprope); // writing local coordinate systems in the auxiliary points
    fprintf(out, "\n");
    fflush(out);
    if (err)
      return(7);
    fprintf(stdout, " O.K.");
  }
  if (Mp->straincomp)
  {
    fprintf(stdout, "\n writing auxiliary points for strain computation . . .");
    err = wr_auxpoint(out, nprope); // writing auxiliary points
    fflush(out);
    if (err)
      return(7);
    fprintf(stdout, " O.K.");
    fprintf(stdout, "\n writing local coordinate systems . . .");
    err = wr_ellcsys(out, nprope); // writing local coordinate systems in the auxiliary points
    fprintf(out, "\n");
    fflush(out);
    if (err)
      return(7);
    fprintf(stdout, " O.K.");
  }*/

  fprintf(stdout, "\n writing cross-sections . . .");
  err = wr_crsecs(out, d);    // writing cross section parameters
  fflush(out);
  if (err)
    return(7);
  fprintf(stdout, " O.K.");

  fprintf(stdout, "\n writing loads . . .");
  err = wr_load(out); // writing loads
  fflush(out);
  if (err)
    return(7);
  fprintf(stdout, " O.K.");

  fprintf(stdout, "\n writing eigenstrains . . .");
  err = wr_eigenstrains(out); // writing eigenstrains
  fflush(out);
  if (err)
    return(8);
  fprintf(stdout, " O.K.");

  long i, j;
  long tnlc = Nlc;
  if (Tnslc)
    tnlc = Tnslc;

  if (Mstrc)
  {
    for (i=0; i<tnlc; i++)
    {
      delete [] Mstrc[i];
    }
    delete [] Mstrc; 
  }
  // deallocating used memory for nodes
  for (i=0; i<Top->nn; i++)
  {
    if (Nod_hang)
      delete Nod_hang[i];
    delete Nod_bocon[i];
    delete [] Nod_ccn[i];
    delete [] Nod_lcs[i];
    if (Nod_nsprmat[i])
    {
      delete [] Nod_nsprmat[i];
      for (j=0; j< Nod_ndof[i]; j++)
      {
        if (Nod_sprmattype[i][j]) 
          delete [] Nod_sprmattype[i][j];
        if (Nod_sprmatid[i][j])
          delete [] Nod_sprmatid[i][j];
        if (Nod_sprmatdbi[i][j])
          delete [] Nod_sprmatdbi[i][j];
      }
      delete [] Nod_sprmattype[i];
      delete [] Nod_sprmatid[i];
      delete [] Nod_sprmatdbi[i];
    }

    for (j=0; j< tnlc; j++)
    {
      delete [] Nod_load[i][j]; 
      delete Nod_temper[i][j];
    }
    for (j=0; j<Nlc; j++)
    {
      delete [] Nod_tdload[i][j]; 
      delete Nod_inicd[i][j];
    }
    delete [] Nod_load[i];
    delete [] Nod_tdload[i];
    delete [] Nod_inicd[i];
    delete [] Nod_temper[i];
  }
  delete [] Nod_hang;  
  delete [] Nod_ndof;  
  delete [] Nod_bocon; 
  delete [] Nod_ccn;   
  delete [] Nod_cst;   
  delete [] Nod_csti;  
  delete [] Nod_cstdbi;  
  delete [] Nod_nsprmat;
  delete [] Nod_sprmattype;
  delete [] Nod_sprmatid;
  delete [] Nod_sprmatdbi;
  delete [] Nod_lcs;   
  delete [] Nod_load;  
  delete [] Nod_tdload;
  delete [] Nod_inicd; 
  delete [] Nod_rot_ipd_axis;
  delete [] Nod_temper;
  delete [] Nod_periodbc;
  
  // deallocating used memory for elements
  for (i=0; i<Top->ne; i++)
  {
    delete [] El_mattype[i];
    delete [] El_matid[i];
    delete [] El_matdbi[i];
    delete [] El_lcs[i];
    for (j=0; j<tnlc; j++)
    {
      delete El_load[i][j];
      delete [] El_loadln[i][j];
      delete [] El_loadcol[i][j];
    }
    for (j=0; j<Nlc; j++)
    {     
      delete El_tdload[i][j];
      delete [] El_tdloadln[i][j];
      delete [] El_tdloadcol[i][j];
    }
    delete [] El_load[i];
    delete [] El_loadln[i];
    delete [] El_loadcol[i];
    delete [] El_tdload[i];
    delete [] El_tdloadln[i];
    delete [] El_tdloadcol[i];
    if (El_eigstr)
      delete [] El_eigstr[i];
  }
  delete [] El_type;   
  delete [] El_ssst;   
  delete [] El_nmat;   
  delete [] El_mattype;
  delete [] El_matid;  
  delete [] El_matdbi;  
  delete [] El_cst;    
  delete [] El_csti;   
  delete [] El_cstdbi;   
  delete [] El_lcs;    
  delete [] El_load;   
  delete [] El_loadln;
  delete [] El_loadcol;
  delete [] El_tdload;   
  delete [] El_tdloadln;
  delete [] El_tdloadcol;
  delete [] El_eigstr;
  delete [] El_eigstrt;
  delete [] El_tfunc;  
  
  for (i=0; i<El_eigstrgf_lst.count(); i++)
    delete (gfunct *)El_eigstrgf_lst.at(i);

  // deallocating rest of objects
  if (Top)      delete Top;
  if (Dbmat)    delete Dbmat;
  if (Dbcrs)    delete Dbcrs;
  if (Nslc)     delete [] Nslc;
  if (Nslc_cum) delete [] Nslc_cum;
  if (Tlt)      delete [] Tlt;
  if (Npd)
  {
    for (i=0; i<Tnslc; i++)
      delete [] Spd[i];
    delete [] Spd;
    delete [] Npd;
  }
  if (Tf)       delete [] Tf;

  return (0);
}



/**
  Function writes section with description of nodes to the text file given by out.

  @param out - pointer to the opened text file, where the data will be written
  @param nprop - number of entries in the input property file

  @retval 0 - if succes
  @retval 1 - if node hasn't assigned ndof
  @retval 2 - if node hasn't assigned cross-section

  Created  by Tomas Koudelka 2009
*/
long wr_nodes(FILE *out)
{
  long i, j, k;

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
        fprintf(out, " %ld", Nod_ndof[i]);
    }
    else
      fprintf(out, " %ld", Nod_ndof[i]);
    fprintf(out, " %d", Nod_cst[i]);
    if (Nod_cst[i])
      fprintf(out, " %ld", Dbcrs->crs[Nod_cstdbi[i]].ridx[Nod_csti[i]-1]);
    if (Nod_lcs[i])
    {
      fprintf(out, " %ld  ", Nod_lcs[i][0].n);
      for(j=0; j<Nod_lcs[i][0].n; j++)
      { 
        for(k=0; k<Nod_lcs[i][j].n; k++)
          fprintf(out, " %le", Nod_lcs[i][j][k]);
        fprintf(out, "   ");
      }
      fprintf(out, "\n");
    }
    else
      fprintf(out, " 0\n");
  }
  fprintf(out, "\n");
  return(0);
}



/**
  Function writes section with description of nodes for layered problems
  to the text file given by out.

  @param out - pointer to the opened text file, where the data will be written

  Returns:
  @retval 0 - on succes

  Created  by Tomas Koudelka 2009
*/
long wr_lay_nodes(FILE *out)
{
  long i, j, k, l;


  fprintf(out, "# sections of nodes\n");
  fprintf(out, "%ld\n", Top->nn*Nlay);
  fprintf(out, "%ld\n", Top->nn);
  for (i = 0; i < Top->nn; i++)
  {
    fprintf (out, "%6ld\n%6ld\n", i+1, Nlay);
    for (l = 0; l < Nlay; l++)
    {
      // writing node number and coordinates
      fprintf(out, "%6ld % e % e % e", Nlay*i+l+1, Top->nodes[i].x, Top->nodes[i].y, Top->nodes[i].z);
      fprintf(out, " %ld %d", Nod_ndof[i], Nod_cst[i]);
      if (Nod_cst[i])
        fprintf(out, " %ld", Dbcrs->crs[Nod_cstdbi[i]].ridx[Nod_csti[i]-1]);
      if (Nod_lcs[i])
      {
        fprintf(out, " %ld  ", Nod_lcs[i][0].n);
        for(j=0; j<Nod_lcs[i][0].n; j++)
        { 
          for(k=0; k<Nod_lcs[i][j].n; k++)
           fprintf(out, " %le", Nod_lcs[i][j][k]);
          fprintf(out, "   ");
        }
        fprintf(out, "\n");
      }
      else
        fprintf(out, " 0\n");
    }
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

  Created  by Tomas Koudelka 2009
*/
long wr_bocon(FILE *out)
{
  long i, k, n, nl, nbc, ndofn;
  char emsg[200];

  nbc = 0;
  Numpd = 0;
  nl = Nlay;
  if (nl == 0)
    nl = 1;
  fprintf(out, "# boundary conditions at nodes\n");
  for (i = 0; i < Top->nn; i++)
  {
    if ((Nod_bocon[i]) || (Nod_ccn[i]))
    // test if any boundary condition common code numbers are prescribed
      nbc++;
  }
  fprintf(out, "%ld\n", nbc*nl);
  for (i = 0; i < Top->nn; i++)
  {
    if ((Nod_bocon[i]) || (Nod_ccn[i]))
    // test if any boundary condition is prescribed or any common code numbers are prescribed
    {
      for (n = 0; n < nl; n++)
      {
        fprintf(out, "%ld", i*nl+n+1);
        ndofn = Nod_ndof[i];
        for (k = 0; k < ndofn; k++)
        { // loop for writing bc
          // if direction has bc, write 0 else write 1
          if (Nod_bocon[i] && Nod_ccn[i])
          {
            if ((Nod_bocon[i]->dir[k] > nobc) && (Nod_ccn[i][k]))
            {
              sprintf(emsg, "node %ld has prescribed boundary condition and common code number\n"                
                            " in the same direction\n.", i+1);
              print_err(emsg, __FILE__, __LINE__, __func__);
              return(1);
            }
            if ((Nod_bocon[i]->dir[k] == nobc) && (Nod_ccn[i][k] == 0))
            // no boundary condition nor common code number
            {
              fprintf(out, " 1");
              continue;
            }
          }
          if (Nod_bocon[i])
          {
            if (Nod_bocon[i]->dir[k] == nobc)
            // no boundary condition
            {
              if (Nod_ccn[i] == NULL)
                fprintf(out, " 1");
            }
            else
            {
              if ((Nod_bocon[i]->nspd[k]) || (Nod_bocon[i]->ndpd[k])) 
              // static or time dependent prescribed displacement
              {
                Numpd++;
                fprintf(out, " %ld", -Numpd);
                continue;
              }

              if (Nod_bocon[i]->con[k] == 0.0)
              // rigid support
              {
                fprintf(out, " 0");
                continue;
              }
            }
          }
          if (Nod_ccn[i])
          // common code number or indeces of time functions
          {
            if (Nod_ccn[i][k])
            {
              fprintf(out, " %ld", Nod_ccn[i][k]);
              continue;
            }
            else
            {
              if (Mp->tprob == growing_mech_structure)
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

  Created  by Tomas Koudelka 2009
*/
long wr_elements(FILE *out, descrip */*d*/)
{
  long i, j, k, l, nspring, eid;
  char emsg[200];
  elemtype et;

  
  nspring = get_nspring();
  fprintf(out, "# section of elements\n");
  fprintf(out, "%ld\n", Top->ne+nspring);
  for (i = 0; i < Top->ne; i++)
  {
    switch (El_type[i])  //  whether element requires cross-section
    {
      case bar2d:
      case beam2d:
      case bar3d:
      case beam3d:
      case beamg3d:
      case barq2d:
      case barq3d:
      case subsoilbeam:
      case beam2dsp:
      {
        if (El_cst[i] == 0)
        {
          sprintf(emsg, "element %ld has not assigned a cross-section", i+1);
          print_err(emsg, __FILE__, __LINE__, __func__);
          return 1;
        }
        break;
      }
      default :
        break;
    }
    if ((Mp->tprob == growing_mech_structure) && (El_tfunc[i] == 0))
    {
      sprintf(emsg, "element number %ld has not assigned a time function", i+1);
      print_err(emsg, __FILE__, __LINE__, __func__);
      return(2);
    }

    // writing type of element, element number
    et = El_type[i];
    fprintf(out, "%6ld ", i+1);
    if ((et >= planeelementlt) && (et <= planeelementsubqt))
      fprintf(out, "%3d %3d", et, El_ssst[i]);
    else
      fprintf(out, "%3d", et);

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
    if (Mp->tprob == growing_mech_structure)
      fprintf(out, " %ld", El_tfunc[i]);
      
    fprintf(out, " %d", El_cst[i]);
    if (El_cst[i] > 0)
      fprintf(out, " %ld", Dbcrs->crs[El_cstdbi[i]].ridx[El_csti[i]-1]); // if element has cs, write cs index

    fprintf(out, " %ld", El_nmat[i]);  // number of material types
    for (l = 0; l < El_nmat[i]; l++)
      fprintf(out, " %d %ld", El_mattype[i][l], Dbmat->mat[El_matdbi[i][l]].ridx[El_matid[i][l]-1]); // writing element material type and index
    fprintf(out, "\n");
  }
  eid = Top->ne+1;
  // writing of spring elements
  fprintf(out, "# spring elements\n");
  for (i=0; i<Top->nn; i++)
  {
    for (j=0; j<Nod_ndof[i]; j++)
    {
      if (Nod_nsprmat[i])
      {
        if (Nod_nsprmat[i][j])
        {
          fprintf(out, "%6ld %3ld %ld 0 0 %ld", eid, spring_1+j, i+1, Nod_nsprmat[i][j]);
          for (k=0; k<Nod_nsprmat[i][j]; k++)
            fprintf(out, " %d %ld", Nod_sprmattype[i][j][k], Dbmat->mat[Nod_sprmatdbi[i][j][k]].ridx[Nod_sprmatid[i][j][k]-1]); // writing spring material type and index
          fprintf(out, "\n");
          eid++;
        }
      }
    }
  }
  fprintf(out, "\n");
  return(0);
}



/**
  Function writes section with description of elements for layered problems to the text
  file given by out.

  @param out - pointer to the opened text file, where the data will be written
  @param d   - pointer to the description structure with preprocessor data

  @retval 0 : on succes
  @retval 1 : if element hasn't assigned type

  Created  by Tomas Koudelka 2009
*/
long wr_lay_elements(FILE *out, descrip */*d*/)
{
  long i, j, l, ln;
  char emsg[200];
  elemtype et;


  fprintf(out, "# section of elements\n");
  fprintf(out, "%ld\n", Top->ne*Nlay);
  for (i = 0; i < Top->ne; i++)
  {
    switch (El_type[i])  //  whether element requires cross-section
    {
      case bar2d:
      case beam2d:
      case bar3d:
      case beam3d:
      case beamg3d:
      case barq2d:
      case barq3d:
      case subsoilbeam:
      case beam2dsp:
      {
        if (El_cst[i] == 0)
        {
          sprintf(emsg, "element %ld has not assigned a cross-section", i+1);
          print_err(emsg, __FILE__, __LINE__, __func__);
          return 1;
        }
        break;
      }
      default :
        break;
    }
    if ((Mp->tprob == growing_mech_structure) && (El_tfunc[i] == 0))
    {
      sprintf(emsg, "element number %ld has not assigned a time function", i+1);
      print_err(emsg, __FILE__, __LINE__, __func__);
      return(2);
    }
    for (ln = 0; ln < Nlay; ln++)
    {
      // writing type of element, element number
      et = El_type[i];

      // writing type of element, element number
      fprintf(out, "%6ld ", Nlay*i+ln+1);
      if ((et >= planeelementlt) && (et <= planeelementsubqt))
        fprintf(out, "%3d %3d", et, El_ssst[i]);
      else
        fprintf(out, "%3d", et);
      // writing element nodes
      for (j = 0; j < Top->elements[i].nne; j++)
      {
/*        if ((et == planeelementqt) && (d->t3d) && (j == 3))
        // nodes in T3d file has in this case folowing order 1 2 3 5 6 4 in view of
        // SIFEL order 1 2 3 4 5 6
          continue;*/
        fprintf(out, " %6ld", (Top->elements[i].nodes[j])*Nlay+ln+1);
      }
/*      if ((et == planeelementqt) && (d->t3d))
      // nodes in T3d file has in this case folowing order 1 2 3 5 6 4 in view of
      // SIFEL order 1 2 3 4 5 6
      // so it is need to write fourth node
        fprintf(out, " %6ld", (Top->elements[i].node[3])*Nlay+ln+1);*/

      // no code numbers on element (code numbers not supported at this time)
      fprintf(out, " 0");

      fprintf(out, " %d", El_cst[i]);
      if (El_cst[i] > 0)
        fprintf(out, " %ld", Dbcrs->crs[El_cstdbi[i]].ridx[El_csti[i]-1]); // if element has cs, write cs index

      fprintf(out, " %ld", El_nmat[i]);  // number of material types
      for (l = 0; l < El_nmat[i]; l++)
        fprintf(out, " %d %ld", El_mattype[i][l], Dbmat->mat[El_matdbi[i][l]].ridx[El_matid[i][l]-1]); // writing element material type and index

      fprintf(out, "\n");
    }
  }
  fprintf(out, "\n");
  return(0);
}



/**
  Function writes section with description of elements to the text file given by out.

  @param out - pointer to the opened text file, where the data will be written

  @retval 0 : on succes

  Created  by Tomas Koudelka 2009
*/
long wr_globnodnum(FILE *out, descrip *d)
{
  long i;
  
  fprintf(out, "# definition of global node numbers\n");
  if (d->paral==2){
    //  sequential version of methods used usually in parallel
    fprintf (out,"%ld %ld\n",Top->meshtype,Top->nsd);
    for (i=0;i<Top->nsd;i++){
      fprintf (out,"%ld\n",Top->nnsd[i]);
    }
  }
  fprintf (out,"\n");
  for (i = 0; i < Top->nn; i++)
    fprintf(out, "%ld\n", Top->gnn[i]);
  return(0);
}



/**
  Function writes materials and thier parameters to the file given by the out.

  @param out - pointer to the opened text file, where the data will be written
  @param d   - pointer to the description structure with preprocessor data

  @retval 0 : on succes

  Created by Tomas Koudelka
*/
long wr_materials(FILE *out, descrip *d)
{
  long i, j, wrt;

  fprintf(out, "# definition of materials\n");
  fprintf(out, "%ld\n", Dbmat->nmatu);  // writing number of material types
  for (i = 0; i < Dbmat->numt; i++)
  {
    wrt = 1;
    for (j = 0; j < Dbmat->mat[i].ninst; j++)
    {
      if (Dbmat->mat[i].instu[j])  // material type index is used (<> 0)
      {
        if (wrt)  // for the first time write material type number and number of used indeces
        {
          fprintf(out, "\n%d %ld\n", Dbmat->mat[i].type, Dbmat->mat[i].ninstu);
          wrt = 0;  // next time don't write
        }
        if (d->matstr == yes)
          // write material type index and property string
          fprintf(out, "%ld %s\n", Dbmat->mat[i].ridx[j], Dbmat->mat[i].inst[j]);
        else
        {
          // write material type index and parameters from mechmat
          fprintf(out, "%ld ", Dbmat->mat[i].ridx[j]);
          Mm->printmatchar(out, Dbmat->mat[i].type, j);
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

  Created by Tomas Koudelka
*/
long wr_crsecs(FILE *out, descrip *d)
{
  long i, j, wrt;

  fprintf(out, "# definition of cross-sections\n");
  fprintf(out, "%ld\n", Dbcrs->ncrsu);    // writing number of cs types
  for (i = 0; i < Dbcrs->numt; i++)
  {
    wrt = 1;
    for (j = 0; j < Dbcrs->crs[i].ninst; j++)
    {
      if (Dbcrs->crs[i].instu[j])  // cs type index is used (<> 0)
      {
        if (wrt) // for the first time write cs type number and number of used indeces
        {
          fprintf(out, "%d %ld\n", Dbcrs->crs[i].type, Dbcrs->crs[i].ninstu);
          wrt = 0;  // next time don't write
        }
        if (d->crsstr == yes)
          // write cs type index and property string
          fprintf(out, "%ld %s\n", Dbcrs->crs[i].ridx[j], Dbcrs->crs[i].inst[j]);
        else
        {
          // write cs type index and parameters from mechcrsec
          fprintf(out, "%ld ", Dbcrs->crs[i].ridx[j]);
          Mc->printcrschar(out, Dbcrs->crs[i].type, j);
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
  @retval 1 : unknown type of problem is required

  Created by Tomas Koudelka
*/
long wr_load(FILE *out)
{
  long tnlc, i, j;

  switch (Mp->tprob)
  {
    case linear_statics:
    case layered_linear_statics:
    case lin_floating_subdomain:
    {
      fprintf(out, "%ld #number of load cases\n\n\n", Nlc);  // writing of number of load cases
      tnlc = Nlc;
      break;
    }
    case eigen_dynamics:
      tnlc = Nlc;
      break;
    case forced_dynamics:{
      tnlc = Nlc;
      if (Tnslc == 0)   tnlc /= 2;
      fprintf(out, "%ld #number of load cases\n\n\n", tnlc);  // writing of number of load cases
      break;
    }
    case mat_nonlinear_statics:
      fprintf(out, "%ld #number of load cases\n\n\n", Nlc/2);  // writing of number of load cases
      tnlc = Nlc;
      break;
    case growing_mech_structure:
    case mech_timedependent_prob:
      fprintf(out, "%ld #number of load cases\n\n\n", Nlc);  // writing of number of load cases
      tnlc = Nlc;
      break;
    default:
      print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
      return 1;
  }
  

  for (i = 0; i < tnlc; i++)
  {
    fprintf(out, "#\n# loadcase number %ld\n#\n", i+1);
    switch (Mp->tprob)
    {
      case linear_statics:
      case lin_floating_subdomain:
        wr_loadn(out, i);          // writing of nodal load
        wr_loadel(out, i);         // writing of element load
        wr_prescdisp(out, i);      // writing of prescribed displacements
        wr_tempload(out, i);
        if ((Mp->homog == 3) || (Mp->homog == 4) || (Mp->homog == 9))
          wr_macrostrastre(out, i);
        fflush(out);
        break;
      case eigen_dynamics:
        break;
      case mat_nonlinear_statics:
        wr_loadn(out, i);          // writing of nodal load
        wr_loadel(out, i);         // writing of element load
        wr_prescdisp(out, i);      // writing of prescribed displacements
        wr_tempload(out, i);
        if ((Mp->homog == 3) || (Mp->homog == 4) || (Mp->homog == 9))
          wr_macrostrastre(out, i);
        if (i%2 == 1)
          wr_initcond(out, i);   // writing of initial conditions
        fflush(out);
        break;
      case forced_dynamics:{
        if (Tnslc == 0)
	{
          fprintf(out, "#\n# static load case (num. %ld in preprocessor)\n#\n\n", 2*i+1);
          wr_loadn(out, 2*i);          // writing of nodal load
          wr_loadel(out, 2*i);         // writing of element load
          wr_prescdisp(out, 2*i);      // writing of prescribed displacements
          wr_tempload(out, 2*i);
          fprintf(out, "\n");
          fprintf(out, "#\n# time dependent load case (num. %ld in preprocessor)\n#\n\n", 2*i+2);
          fprintf(out, "\n %d # load case type = time dependent load\n", int(Tdload));
          wr_dloadn(out, 2*i+1);     // writing of dynamic nodal load
          wr_dloadel(out, 2*i+1);    // writing of dynamic element load
          wr_dprescdisp(out, 2*i+1); // writing of dynamic prescribed displacements
          wr_initcond(out, 2*i+1);
	}
        else
	{
          fprintf(out, "#\n# static load case\n#\n\n");
          wr_loadn(out, Nslc_cum[i]);          // writing of nodal load
          wr_loadel(out, Nslc_cum[i]);         // writing of element load
          wr_prescdisp(out, Nslc_cum[i]);      // writing of prescribed displacements
          wr_tempload(out, Nslc_cum[i]);
          fprintf(out, "\n");
          //  type of dynamic load; 1-usual load case, 10-seismic load
          fprintf(out, "#\n# time dependent load case\n#\n\n");
          fprintf(out, "\n 1 # load case type = time independent load with time scaling factor\n");
          fprintf (out,"%ld # number of subloadcases\n", Nslc[i]);
          for (j = 0; j < Nslc[i]; j++)
	  {
            fprintf(out, "# subloadcase number %ld\n#\n", j+1);
            // writing nodal load
            wr_loadn(out, Nslc_cum[i]+j+1);
            // writing element load
            wr_loadel(out, Nslc_cum[i]+j+1);
            // writing prescribed displacements
            wr_prescdisp(out, Nslc_cum[i]+j+1);
            // writing temperature load type
            wr_tempload(out, Nslc_cum[i]+j+1);      
            fprintf(out, "# time function of subloadcase coefficient\n");
            Tf[Nslc_cum[i]+j+1].print(out);
	  }
          wr_initcond(out, i);
	}
        fflush(out);
	break;
      }
      case growing_mech_structure:
      case mech_timedependent_prob:
      {
        if (Tnslc > 0)
        {
          //  type of dynamic load; 1-usual load case, 10-seismic load
          fprintf (out,"\n1 # type of load is usual load cases\n");
          fprintf (out,"%ld # number of subloadcases\n", Nslc[i]);
          for (j = 0; j < Nslc[i]; j++)
          {
            fprintf(out, "\n# subloadcase number %ld\n#\n", j+1);
            // writing nodal load
            wr_loadn(out, Nslc_cum[i]+j);
            // writing element load
            wr_loadel(out, Nslc_cum[i]+j);
            // writing prescribed displacements
            wr_prescdisp(out, Nslc_cum[i]+j);
            // writing temperature load type
            wr_tempload(out, Nslc_cum[i]+j);
            fprintf(out, "# time function of subloadcase coefficient\n");      
            Tf[Nslc_cum[i]+j].print(out);
          }
        }
        else
        {
          fprintf (out,"\n 20 # type of load is general time dependent load\n");
          // writing dynamic nodal load
          wr_dloadn(out, i);
          // writing dynamic element load
          wr_dloadel(out, i);
          // writing dynamic prescribed displacements
          wr_dprescdisp(out, i); 
        }
        wr_initcond(out, i);
        if (Mp->tprob == growing_mech_structure)
          wr_rotinidispl(out, i);
        fflush(out);
        break;
      }
      default:
        print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
        return 1;
    }
  }
  
  return(0);
}



/**
  Function writes section with description of nodal load to the text file given by out

  @param out   - pointer to the opened text file, where the data will be written
  @param nlc   - number of loading case which will be written

  @retval 0 : on succes

  Created by Tomas Koudelka
*/
long wr_loadn(FILE *out, long nlc)
{
  long i, k, n, nl, nln;
  nl = Nlay;
  if (nl == 0)
    nl = 1;
  fprintf(out, "# nodal load\n");
  for (i = 0, nln = 0; i < Top->nn; i++)
  // countig total number of nodal loads
  {
    if (Nod_load[i])
    {
      if (Nod_load[i][nlc])
        nln++;
    }
  }
  fprintf(out, "%ld\n", nln*nl);  // writing total number of nodal loading entries
  for (i = 0; i < Top->nn; i++)
  {
    if (Nod_load[i])
    {
      if (Nod_load[i][nlc])
      {
        for (n = 0; n < nl; n++)
        {
          // writing node number
          fprintf(out, "%6ld", i*nl+n+1);
          for (k = 0; k < Nod_ndof[i]; k++)
          // loop for writing values of load in single directions
            fprintf(out, " % e", Nod_load[i][nlc][k]);
          fprintf(out, "\n");
        }
      }
    }
  }
  fprintf(out, "\n");
  return(0);
}



/**
  Function writes section with description of dynamic nodal load to the text file given by out

  @param out   - pointer to the opened text file, where the data will be written
  @param nlc   - number of loading case which will be written

  @retval 0 : on succes

  Created by Tomas Koudelka
*/
long wr_dloadn(FILE *out, long nlc)
{
  long i, j, nln;
  fprintf(out, "# nodal load\n");
  for (i = 0, nln = 0; i < Top->nn; i++)
  // counting of total number of nodal loads
  {
    if (Nod_tdload[i]) 
    {
      if (Nod_tdload[i][nlc])
      {
        if ((Mp->tprob == forced_dynamics) || (Mp->tprob == mech_timedependent_prob) ||
            (Mp->tprob == growing_mech_structure))
          nln++;
      }
    }
  }
  fprintf(out, "%ld\n", nln);  // writing total number of nodal loading entries
  for (i = 0; i < Top->nn; i++)
  {
    if (Nod_tdload[i]) 
    {
      if (Nod_tdload[i][nlc])
      {
        if ((Mp->tprob == forced_dynamics) || (Mp->tprob == mech_timedependent_prob) ||
            (Mp->tprob == growing_mech_structure))
        {
          // writing node number
          fprintf(out, "\n%ld # node number\n", i+1);
          for (j=0; j<Nod_ndof[i]; j++)
          {
            Nod_tdload[i][nlc][j].print(out); 
          }
        }
      }
    }
  }
  fprintf(out, "\n\n");
  return(0);
}



/**
  Function writes section with description of element load to the text file given by out

  @param out   - pointer to the opened text file, where the data will be written
  @param nlc   - number of loading case which will be written

  @retval 0 : on succes

  Created by Tomas Koudelka
*/
long wr_loadel(FILE *out, long nlc)
{
  fprintf(out, "# load of elements\n");
  long i, n, nl, nle;
  nl = Nlay;
  if (nl == 0)
    nl = 1;
  for (i = 0, nle = 0; i < Top->ne; i++)
  {  // counting total number of loaded elements
    if (El_load[i])
    {
      if(El_load[i][nlc])
      {
        if ((El_load[i][nlc]->nodvale) ||
            (El_load[i][nlc]->nodvals) ||
            (El_load[i][nlc]->nodvalv))
          nle++; 
      }
    }
  }

  // writing total number of element loads entries
  fprintf(out, "%ld\n", nle*nl);
  for (i = 0; i < Top->ne; i++)
  {
    for (n = 0; n < nl; n++)
    {
      if (El_load[i])
      {
        if (El_load[i][nlc])
        {
          // writing of element number
          fprintf(out, "%6ld", nl*i+n+1);
          // writing of load
          El_load[i][nlc]->print(out, 6);
        }
      }
    }
  }
  fprintf(out, "\n");
  return(0);
}



/*
  Function writes section with description of time dependent element load to the text file given by out

  @param out   - pointer to the opened text file, where the data will be written
  @param nlc   - number of loading case which will be written

  @retval 0 : on succes

  Created by Tomas Koudelka
*/
long wr_dloadel(FILE *out, long nlc)
{
  fprintf(out, "# time dependent load of elements\n");
  long i, nle;
  for (i = 0, nle = 0; i < Top->ne; i++)
  {  // counting total number of loaded elements
    if (El_tdload[i])
    {
      if(El_tdload[i][nlc])
      {
        if ((El_tdload[i][nlc]->nodvale) ||
            (El_tdload[i][nlc]->nodvals) ||
            (El_tdload[i][nlc]->nodvalv))
          nle++; 
      }
    }
  }

  // writing total number of element loads entries
  fprintf(out, "%ld\n", nle);
  for (i = 0; i < Top->ne; i++)
  {
    if (El_tdload[i])
    {
      if (El_tdload[i][nlc])
      {
        // writing of element number
        fprintf(out, "%6ld", i+1);
        // writing of load
        El_tdload[i][nlc]->print(out, 6);
      }
    }
  }
  fprintf(out, "\n");
  return(0);
}



/**
  Function writes section with description of nodal prescribed
  displacement load to the text file given by out.

  @param out   - pointer to the opened text file, where the data will be written
  @param nlc   - number of loading case which will be written

  @retval 0 : on succes
  @retval 1 : wrong number of prescribed displacements (error in code)

  Created by Tomas Koudelka
*/
long wr_prescdisp(FILE *out, long nlc)
{
  long i, j, ndofn, n;
  double tmp;
  const char *coordn[3];
  coordn[0] = "x";
  coordn[1] = "y";
  coordn[2] = "z";
  vector p(ASTCKVEC(3));
  

  fprintf(out, "# prescribed displacements\n");
  if (Mp->tprob == growing_mech_structure)
  {
    fprintf(out, "%ld\n", Npd[nlc]);
    if (Npd[nlc] == 0)
      return(0);
    for (i=0; i<Npd[nlc]; i++)
      fprintf(out, "% e\n", Spd[nlc][i]);
    return(0);
  }

  // total number of prescribed displacements
  fprintf(out, "%ld\n", Numpd);
  if (Numpd == 0)
    return(0);
  n = 0;
  for (i = 0; i < Top->nn; i++)
  {
    if (Nod_bocon[i])
    // test if any boundary condition is prescribed
    {
      ndofn = Nod_ndof[i];
      for (j = 0; j < ndofn; j++)
      // loop for writing bc
      // if direction has bc, write p.d. number and value
      {
        if ((Nod_bocon[i]->dir[j] > nobc) && (Nod_bocon[i]->nspd[j] || (Nod_bocon[i]->ndpd[j])))
        {
          switch (Nod_bocon[i]->dir[j]){
            case sbc:
            case tdbc:
              fprintf(out, "% e\n", Nod_bocon[i]->con[ndofn*nlc+j]);
              n++;
              break;
            case gbc:
              Top->nodes[i].getcoord(p);
              tmp = Nod_bocon[i]->gf[ndofn*nlc+j]->getval(p, coordn);
              fprintf(out, "% e\n", tmp);
              n++;
              break;
            default:
              print_err("unknown type of boundary condition is required", __FILE__, __LINE__, __func__);
              abort();
          }
        }
      }
    }
  }
  if (Numpd != n)
  {
    print_err("wrong number of prescribed displacements", __FILE__, __LINE__, __func__);
    return 1;
  }
  fprintf(out, "\n");
  return(0);
}



/**
  Function writes section with description of dynamic nodal prescribed
  displacement load to the text file given by out.

  @param out   - pointer to the opened text file, where the data will be written
  @param nlc   - number of loading case which will be written

  Returns:
  @retval 0 : on succes

  Created by Tomas Koudelka
*/
long wr_dprescdisp(FILE *out, long nlc)
{
  long i, j, n, ndofn;

  // total number of prescribed displacements
  fprintf(out, "# time dependent prescribed displacements\n");
  fprintf(out, "%ld\n", Numpd);
  if (Numpd == 0)
    return(0);
  n = 0;
  for (i = 0; i < Top->nn; i++)
  {
    if (Nod_bocon[i])
    // test if any boundary condition is prescribed
    {
      ndofn = Nod_ndof[i];
      for (j = 0; j < ndofn; j++)
      // loop for writing bc
      // if direction has bc, write p.d. number and value
      {
        if ((Nod_bocon[i]->dir[j] > nobc) && Nod_bocon[i]->ndpd[j])
        {
          Nod_bocon[i]->gf[ndofn*nlc+j]->print(out);
          n++;
          continue;
        }
        if ((Nod_bocon[i]->dir[j] > nobc) && Nod_bocon[i]->nspd[j])
        {
          fprintf(out, "0 %e\n", Nod_bocon[i]->con[ndofn*nlc+j]);
          n++;
          continue;
        }
      }
    }
  }
  if (Numpd != n)
  {
    print_err("wrong number of prescribed displacements", __FILE__, __LINE__, __func__);
    return 1;
  }
  fprintf(out, "\n");
  return(0);
}



/**
  Function writes temperature load(temperature changes) at nodes.
  
  @param out   - pointer to the opened text file, where the data will be written
  @param nlc   - number of loading case which will be written

  Returns:
  @retval 0 - on success
  @retval 1 - incompatible tempreature load type and assigned temperature load

  Created by Tomas Koudelka
*/
long wr_tempload(FILE *out, long nlc)
{
  long i, prnzero;

  fprintf(out, "# temperature load\n");
  if (Tlt)
  {
    fprintf(out, "%ld\n", Tlt[nlc]);
    if ((Tlt[nlc] == 1) || (Tlt[nlc] == 2))
    {
      for (i = 0; i < Top->nn; i++)
      {
        prnzero = 1;
        if (Nod_temper[i][nlc])
        {
          fprintf(out, "% g\n", Nod_temper[i][nlc]->val);
          prnzero = 0;
        }
        if (prnzero)
          fprintf(out, "0.0\n");
      }
      fprintf(out, "\n");
    }
  }
  else
  // no temperature load was specified
    fprintf(out, "0\n");

  return 0L;
}



/**
  Function writes section with prescribed macro strains/stresses

  @param out - pointer to the opened text file, where the data will be written
  @param nlc - number of loading case which will be written

  @retval 0 : on success

  Created by Tomas Koudelka, 18.2.2014
*/
long wr_macrostrastre(FILE *out, long nlc)
{
  long i;

  if (Mp->homog == 3)
    fprintf(out, "\n# macro stress components\n");
  if (Mp->homog == 4)
    fprintf(out, "\n# macro strain components\n");
  if (Mp->homog == 9)
    fprintf(out, "\n# macro stress/strain components\n");
  if (Mstrc[nlc])
  {
    for (i=0; i<Nmstrc; i++){
      if (Mp->homog == 9)
        fprintf(out, "%d %le\n", Mstrct[nlc][i], Mstrc[nlc][i]);
      else
        fprintf(out, "%le ", Mstrc[nlc][i]);
    }
  }
  else
  {
    if (Mp->homog == 9){
      print_err("no macro stress/strain components were defined in the load case %ld",
                __FILE__, __LINE__, __func__, nlc+1);
      abort();
    }
    for (i=0; i<Nmstrc; i++)
      fprintf(out, "%le ", 0.0);
  }
  fprintf(out, "\n");

  return (0);
}



/**
  Function writes section with description of initial conditions
  to the text file given by out.

  @param out   - pointer to the opened text file, where the data will be written
  @param nlc   - number of loading case which will be written

  @retval 0 : on succes

  Created by Tomas Koudelka
*/
long wr_initcond(FILE *out, long nlc)
{
  long i, j, nini;

  fprintf(out, "# initial conditions\n");
  if (Nod_inicd[0][nlc])
  {
    nini = 0;
    for (j=0; j<Top->nn; j++)
    {
      if (Nod_inicd[j][nlc])
        nini++;
    }
    fprintf(out, "%ld\n", nini);
    for (i = 0; i < Top->nn; i++)
    {
      if (Nod_inicd[i][nlc])
      {
        // writing node number
        fprintf(out, "%6ld ", i+1);
        Nod_inicd[i][nlc]->print(out);
      }
    }
    fprintf(out, "\n");
  }
  else
    fprintf(out, "0\n\n");
  return(0);
}



/**
  Function writes section with 


  @param out   - pointer to the opened text file, where the data will be written
  @param nlc   - number of loading case which will be written

  @retval 0 : on succes

  Created by Tomas Koudelka
*/
long wr_rotinidispl(FILE *out, long /*nlc*/)
{
  long i, n=Nod_rot_ipd_num;

  if (Mp->rot_inidispl == yes){
    fprintf(out, "\n# prescribed initial displacements by rotation\n\n");
    fprintf(out, "%ld #number of axes defined\n", n);
    for (i=0; i<n; i++){
      fprintf(out, "\n# Axis %ld:\n", i+1);
      Nod_rot_ipd_axis[i].print(out);
    }
  }
  return(0);
}



/**
  Function writes section with eigenstrains/eigenstresses on element to the text file given by out.

  @param out   - pointer to the opened text file, where the data will be written

  @retval 0 : on succes

  Created  by Tomas Koudelka 11.2012
*/
long wr_eigenstrains(FILE *out)
{
  long i, j, add_nullgf, ncomp;
  gfunct *tmp;
  const char *strtn;

  // eigenstrains/eigenstresses
  if (El_eigstr == NULL) // no eigenstrains/eigenstresses defined
  {
    fprintf(out, "\n#eigenstrains/eigenstresses\n"); // eigenstrains are defined for each element
    fprintf(out, "0\n");
    return 0;
  }
  else
  {
    strtn = strastre_kwdset.get_str(El_eigstrt[0]);
    if (El_eigstrt[0] == strain)
      fprintf(out, "\n#eigenstrains\n2\n"); // eigenstrains are defined for each element
    else
      fprintf(out, "\n#eigenstresses\n5\n"); // eigenstresses are defined for each element

    add_nullgf = 0;
    for (i=0; i<Top->ne; i++)
    { // output of gfunct id of eigenstrain/eigenstress components on particular elements
      ncomp = Mt->give_tncomp(El_type[i]);
      if (El_ssst[i])
        ncomp = 4;
      for (j=0; j<ncomp; j++)
      {
        if (El_eigstr[i])
          fprintf(out, "%ld ", El_eigstr[i][j]+1);
        else
        {
          fprintf(out, "%ld ", El_eigstrgf_lst.count()+1);
          add_nullgf = 1;
        }
      }
      fprintf(out,"\n");
    }
    // output time functions for eigenstrain/eigenstress components

    if (add_nullgf) // constant zero returning time function has to be added
      fprintf(out, "\n#eigen%s functions\n\n%ld #number of eigen%s functions\n", strtn, El_eigstrgf_lst.count()+1, strtn);
    else
      fprintf(out, "\n#eigen%s functions\n\n%ld #number of eigen%s functions\n", strtn, El_eigstrgf_lst.count(), strtn);
    for (i=0; i<El_eigstrgf_lst.count(); i++)
    {
      fprintf(out, "# %ld. time function for eigen%s\n", i+1, strtn);
      tmp = (gfunct *)El_eigstrgf_lst.at(i);
      tmp->print(out);
    }
    if (add_nullgf)
    {
      fprintf(out, "# %ld. time function for eigen%s - null time function\n", i+1, strtn);
      fprintf(out, "0 0.0\n");
    }
  }
  return 0;
}



/**
  The function computes number of assigned spring supports at all nodes.

  @returns The function returns the total number of spring supports.

  Created by Tomas Koudelka, 2010
*/
long get_nspring()
{
  long i, j, ret;

  ret = 0;
  for(i=0; i<Top->nn; i++)
  {
    for(j=0; j<Nod_ndof[i]; j++)
    {
      if (Nod_nsprmat[i])
      {
        if (Nod_nsprmat[i][j])
          ret++;
      }
    }
  }
  return ret;
}



/**
  The function writes periodic boundary condition to the given text file.

  @param[in] out - pointer to the opened text file for output

  @retval 0 - on success
  @retval 1 - in the case of output error
*/
long wr_periodbc(FILE *out)
{
  long i, failed = 0;
  ivector mpbc(Top->nn);
  
  fprintf(out, "\n\n# number of nodes with periodic boundary condition:\n%ld\n", Nperbc);
  fprintf(out, "# slave_node_id, master_node_id\n");
  fprintf(out, "# or\n");
  fprintf(out, "# slave_node_id, master_elem_id, xi, eta, zeta\n");
  for (i=0; i<Top->nn; i++){
    if (Nod_periodbc[i].slmas == 2){
      mpbc(i) = 1;
      failed++;
      continue;
    }
    if (Nod_periodbc[i].slmas > 0){
      if (Nod_periodbc[i].nid >= 0){
        fprintf(out, "%ld node %ld\n", i+1, Nod_periodbc[i].nid+1);
        continue;
      }
      if (Nod_periodbc[i].eid >= 0){
        fprintf(out, "%ld elem %ld %le %le %le\n",
                i+1, Nod_periodbc[i].eid+1,
                Nod_periodbc[i].xi, Nod_periodbc[i].eta, Nod_periodbc[i].zeta);
        continue;
      }
      fprintf(out, "%ld no element was found\n", i+1);
      mpbc(i) = -1;
      failed++;
    }
  }
  fprintf(out, "# number of nodes with periodic bc failure: %ld\n", failed);  
  fprintf(out, "# list of nodes with multiple periodic BC:\n");  
  for(i=0; i<Top->nn; i++){
    if (mpbc(i) > 0)
      fprintf(out, "%6ld entity %ld\n", i+1, long(Nod_periodbc[i].ent));
  }
  fprintf(out, "# list of nodes where no master nodes/element were found:\n");  
  for(i=0; i<Top->nn; i++){
    if (mpbc(i) < 0)
      fprintf(out, "%6ld\n", i+1);
  }
  
  return 0;
}


/**
  The function transforms periodic boundary conditions to the hanging node fromat.

  @param[in] out - pointer to the opened text file for the log/warning output

  @retval 0 - on success
  @retval 1 - in the case of output error
*/
void transf_periodbc_hangnod()
{
  long i, j, eid, nne, ndofn;
  gtypel et;
  ivector enod;
  ivector mnod(ASTCKIVEC(1));
  vector ncoord(ASTCKVEC(3));
  bocon tbc, *aux;

  if (Nod_hang == NULL){
    Nod_hang = new hangnode*[Top->nn];
    memset (Nod_hang, 0, sizeof(*Nod_hang)*Top->nn);
    Numhn = 0;
  }
  
  for (i=0; i<Top->nn; i++){
    if (Nod_hang[i]){
      print_err("Multiple hanging node definition at node %ld\n", __FILE__, __LINE__, __func__, i+1);
      abort();
    }
    if (Nod_periodbc[i].slmas == 1){
      eid = Nod_periodbc[i].eid;
      if (eid >= 0){ // master nodes were found on element
        nne = Top->elements[eid].nne;
        makerefv(enod, Top->elements[eid].nodes, nne);
        et = Top->elements[eid].type;
        ncoord(0) = Nod_periodbc[i].xi;
        ncoord(1) = Nod_periodbc[i].eta;
        ncoord(2) = Nod_periodbc[i].zeta;
        Nod_hang[i] = new hangnode(enod, ncoord, et);
        // increase element node numbers due to printing out of the hanging nodes
        for (j=0; j<nne; j++)
          Nod_hang[i]->mnodes[j]++;
        Numhn++;
        continue;
      }
      if (Nod_periodbc[i].nid >= 0){  // master node was found directly
        print_err("direct master node mapping has not yet been implemented for periodic boundary conditions\n",
                  __FILE__, __LINE__, __func__);
        abort();
      }
      else{
        print_err("no element nor master node was found in periodic boundary condition at node %ld", __FILE__, __LINE__, __func__, i+1);
        abort();
      }
      continue;
    }
    if (Nod_periodbc[i].slmas == 2){
      if (Nod_periodbc[i].ent == gcurve){
        ndofn = Nod_ndof[i];
        tbc.init(ndofn);
        // mulitple periodic BC at nodes on edges -> corner nodes -> convert to rigid support
        for (j=0; j<ndofn; j++){
          tbc.dir[j] = sbc;
          tbc.con[j] = 0.0;
        }
        if (Nod_bocon[i])
        {
          aux = Nod_bocon[i];
          Nod_bocon[i] = tbc.merge(Nod_bocon[i]);
          delete aux;
          if (Nod_bocon[i] == NULL)
          {
            print_err("Conflict in merging of periodic boundary condition with the Dirichlet boundary conditions at node %ld\n", __FILE__, __LINE__, __func__, i+1);
            abort();
          }
        }
        else
          Nod_bocon[i] = tbc.copy();
      }
      else{
        print_err("multiple master/slave node assignments were found in periodic boundary condition at node %ld", __FILE__, __LINE__, __func__, i+1);
        return;
      }
    }
  }
}
