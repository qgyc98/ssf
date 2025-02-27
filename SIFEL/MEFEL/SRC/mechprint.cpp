#include "mechprint.h"
#include "global.h"
#include "globmat.h"
#include "alias.h"
#include "gtopology.h"
#include "meshtransfer.h"
#include "loadcase.h"
#include "element.h"
#include "node.h"
#include "intpoints.h"
#include "plelemqq.h"
#include "mathem.h"
#include "siftop_element_types.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechbclc.h"
#include "lhsrhs.h"
#include "outdriverm.h"
#include "siftop.h"
#include "elemhead.h"
#include "intp.h"
#include "vector.h"
#include "matrix.h"
#include "vecttens.h"
#include "elemswitch.h"
#include "gfmatrix.h"


#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <set>



/**
  Function opens data files for each type of outputfile

  @param istep - step id. The parameter enables to open file with name enhanced
                 by the step id - istep >= 0 or it leaves required filename untouched for instance that
                 istep < 0. In case of stochastic calculations and in case istep >= 0, the istep precedes the 
                 stochastic step id.
  @param mode - string with control sequence for the file opening. It enables to open 
                new file (mode = "wt") or to append existing ones (mode = "at").
  @param idn1 - id of the first node for GiD mesh (default is idn1 = 1 -> nodes are numbered from 1)
  @param ide1 - id of the first element for GiD mesh (default is ide1 = 1 -> elements are numbered from 1)

  @return The function does not return anything.

  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz

*/
void print_init(long istep, const char *mode, long idn1, long ide1)
{
  char fname[FNAMELEN+20];
  char *path;
  char *name;
  char *suffix;
  long i;
  long sl;


  if ((Outdm->outf == NULL) && (Outdm->textout==1))
  {
    filename_decomposition (Outdm->outfn,path,name,suffix);
    if (St == NULL)
    {
      if (istep < 0)
        sprintf(fname, "%s%s%s", path, name, suffix);
      else
        sprintf(fname, "%s%s.%ld%s", path, name, istep+1, suffix);
    }
    else
    {
      if (istep < 0)
        sprintf(fname, "%s%s.%ld%s", path, name, Mp->ns+1, suffix);
      else
        sprintf(fname, "%s%s.%ld.%ld%s", path, name, istep+1, Mp->ns+1, suffix);
    }
    Outdm->outf = fopen(fname, mode);
    if (Outdm->outf == NULL)
    {
      print_err("unable to open output text file '%s'", __FILE__, __LINE__, __func__, fname);
      abort();
    }
    fseek(Outdm->outf, 0, SEEK_END); // MS Visual C++ requires that
    if (ftell(Outdm->outf) == 0)
      Outdm->print_header(Outdm->outf);

    delete [] path;
    delete [] name;
    delete [] suffix;
  }

  if ((Outdm->gf != grfmt_no) && (Outdm->outgr == NULL))
  {
    filename_decomposition (Outdm->outgrfn,path,name,suffix);
    if (St == NULL)
    {
      if (istep < 0)
        sprintf(fname, "%s%s%s", path, name, suffix);
      else
        sprintf(fname, "%s%s.%ld%s", path, name, istep+1, suffix);
    }
    else
    {
      if (istep < 0)
        sprintf(fname, "%s%s.%ld%s", path, name, Mp->ns+1, suffix);
      else
        sprintf(fname, "%s%s.%ld.%ld%s", path, name, istep+1, Mp->ns+1, suffix);
    }
    if ((Outdm->gf == grfmt_gid) || (Outdm->gf == grfmt_gid_sep) || (Outdm->gf == grfmt_gid_vtk) || (Outdm->gf == grfmt_gidsep_vtk))
    {
      sl = long(strlen(fname));
      sprintf(fname+sl, ".msh");
    }
    if(Outdm->gf != grfmt_vtk)
      Outdm->outgr = fopen(fname, mode);
    if (Outdm->outgr == NULL && Outdm->gf != grfmt_vtk)
    {
      print_err("unable to open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
      abort();
    }
    switch(Outdm->gf)
    {
      case grfmt_no:
        break;   
      case grfmt_open_dx:
        break;
      case grfmt_femcad:
        export_femcad(Outdm->outgr);
        break;
      case grfmt_gid:
      case grfmt_gid_vtk:
        fseek(Outdm->outgr, 0, SEEK_END); // MS Visual C++ requires that
        if (ftell(Outdm->outgr) == 0)
          export_gid_mesh(Outdm->outgr, idn1, ide1);
        fclose(Outdm->outgr);
        if (Outdm->ncut > 0)
	      {
          sprintf(fname+sl, "2d.msh");
          Outdm->outgr = fopen(fname, mode);
          fseek(Outdm->outgr, 0, SEEK_END); // MS Visual C++ requires that
          if (ftell(Outdm->outgr) == 0)
            export_gid_2dmesh(Outdm->outgr, Outdm->ncut, idn1, ide1);
          fclose(Outdm->outgr);
	      }
        sprintf(fname+sl, ".res");
        Outdm->outgr = fopen(fname, mode);
        if (Outdm->outgr == NULL)
        {
          print_err("unable to open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        if (strchr(mode, 'w') != NULL)
        {
          fseek(Outdm->outgr, 0, SEEK_END); // MS Visual C++ requires that
          if (ftell(Outdm->outgr) == 0)
	  {
            fprintf(Outdm->outgr, "GiD Post Results File 1.0\n");
            export_gid_gauss_pt(Outdm->outgr, ide1);
	  }
        }
        break;
      case grfmt_gid_sep:
      case grfmt_gidsep_vtk:
        fseek(Outdm->outgr, 0, SEEK_END); // MS Visual C++ requires that
        if (ftell(Outdm->outgr) == 0)
          export_gid_mesh(Outdm->outgr, idn1, ide1);
        fclose(Outdm->outgr);
        if (Outdm->ncut > 0)
        {
          sprintf(fname+sl, "2d.msh");
          Outdm->outgr = fopen(fname, mode);
          fseek(Outdm->outgr, 0, SEEK_END); // MS Visual C++ requires that
          if (ftell(Outdm->outgr) == 0)
            export_gid_2dmesh(Outdm->outgr, Outdm->ncut, idn1, ide1);
          fclose(Outdm->outgr);
        }
        strncpy(Outdm->outgrfngs, fname, sl);
        Outdm->outgrfngs[sl] = 0;
        Outdm->outgr = NULL;
        Outdm->create_files_gidsp(mode);
        break;
      case grfmt_vtk:
        break;
      default:
        print_err("unknown type of graphics format is required", __FILE__, __LINE__, __func__);
    }
    delete [] path;
    delete [] name;
    delete [] suffix;
  }
  if ((Outdm->ndiag > 0) && (Outdm->outdiagf[0] == NULL))
  {
    filename_decomposition (Outdm->outdiagfn,path,name,suffix);
    for (i=0; i<Outdm->ndiag; i++)
    {
      if (St == NULL)
      {
        if (Outdm->ndiag > 1)
          sprintf(fname, "%s%s.%ld%s", path, name, i+1, suffix);
        else
          sprintf(fname, "%s%s%s", path, name, suffix);
      }
      else
      {
        if (Outdm->ndiag > 1)
          sprintf(fname, "%s%s.%ld.%ld%s", path, name, i+1, Mp->ns+1, suffix);
        else
          sprintf(fname, "%s%s.%ld%s", path, name, Mp->ns+1, suffix);
      }
      if (Outdm->outdiagf[i] == NULL)
      {
        Outdm->outdiagf[i] = fopen(fname, mode);
        if (Outdm->outdiagf[i] == NULL)
        {
          print_err("unable to open diagram file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
      }
    }
    delete [] path;
    delete [] name;
    delete [] suffix;
  }
  Outdm->idn1 = idn1;
  Outdm->ide1 = ide1;


}



/**
  The function performs all prints requirements for given step.

  @param lcid - load case id
  @param istep - step id
  @param lambda - load coefficient or time (it depends on problem type and solver)
  @param fi - array of additional values at nodes which should be printed.
              It depends on problem and solver type. For example for nonlinear
              statics it contains the internal forces.
  
  @return The function does not return anything.

  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz

*/
void print_step(long lcid, long istep, double lambda, double *fi)
{
  switch(Mp->tprob){
    case linear_statics:
    case load_balancing:
    case layered_linear_statics:
      Outdm->print_out(Outdm->outf, lcid, istep, lambda);
      Outdm->print_graphics(Outdm->outgr, lcid, lambda, istep, fi);
      break;
    case lin_floating_subdomain:
    case mat_nonlinear_statics:
    case geom_nonlinear_statics:
    case earth_pressure:
    case forced_dynamics:
    case mech_timedependent_prob:
    case nonlin_floating_subdomain:
    case growing_mech_structure:
      Outdm->print_newstep(Outdm->outf, lcid, istep, lambda);
      Outdm->print_out(Outdm->outf, lcid, istep, lambda);
      Outdm->print_diags(lcid, lambda, istep, fi);
      Outdm->print_graphics(Outdm->outgr, lcid, lambda, istep, fi);
      break;
    case eigen_dynamics:
      Outdm->print_newstep(Outdm->outf, lcid, istep, lambda);
      Outdm->print_out(Outdm->outf, lcid, istep, lambda);
      Outdm->print_graphics(Outdm->outgr, lcid, lambda, istep, fi);
      break;
    default:
      print_err("unsupported problem type is required", __FILE__, __LINE__, __func__);
      break;
  }
}



/**
  The function performs all prints requirements for given step.

  @param lcid - load case id
  @param istep - step id
  @param lambda - load coefficient or time (it depends on problem type and solver)
  @param fi - array of additional values at nodes which should be printed.
              It depends on problem and solver type. For example for nonlinear
              statics it contains the internal forces.
  @param fr - array of residual %vector components at nodes which should be printed.
  
  @return The function does not return anything.

  Created 01.2024 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void print_step(long lcid, long istep, double lambda, double *fi, double *fr)
{
  switch(Mp->tprob){
    case mat_nonlinear_statics:
    case geom_nonlinear_statics:
    case earth_pressure:
    case mech_timedependent_prob:
    case nonlin_floating_subdomain:
    case growing_mech_structure:
      Outdm->print_newstep(Outdm->outf, lcid, istep, lambda);
      Outdm->print_out(Outdm->outf, lcid, istep, lambda);
      Outdm->print_diags(lcid, lambda, istep, fi, fr);
      Outdm->print_graphics(Outdm->outgr, lcid, lambda, istep, fi, fr);
      break;
    default:
      print_err("unsupported problem type is required", __FILE__, __LINE__, __func__);
      break;
  }
}



/**
  The function performs all prints requirements  without respect to
  to selected step/time (forced printing).

  @param lcid - load case id
  @param istep - step id
  @param lambda - load coefficient or time (it depends on problem type and solver)
  @param fi - array of additional values at nodes which should be printed.
              It depends on problem and solver type. For example for nonlinear
              statics it contains the internal forces.
  
  @return The function does not return anything.

  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz

*/
void print_step_forced(long lcid, long istep, double lambda, double *fi)
{
  switch(Mp->tprob){
    case linear_statics:
    case load_balancing:
    case layered_linear_statics:
      Outdm->print_out_forced(Outdm->outf, lcid, istep, lambda);
      Outdm->print_graphics_forced(Outdm->outgr, lcid, lambda, istep, fi);
      break;
    case lin_floating_subdomain:
    case mat_nonlinear_statics:
    case geom_nonlinear_statics:
    case earth_pressure:
    case forced_dynamics:
    case mech_timedependent_prob:
    case nonlin_floating_subdomain:
    case growing_mech_structure:
      Outdm->print_newstep(Outdm->outf, lcid, istep, lambda);
      Outdm->print_out_forced(Outdm->outf, lcid, istep, lambda);
      Outdm->print_diags_forced(lcid, lambda, istep, fi);
      Outdm->print_graphics_forced(Outdm->outgr, lcid, lambda, istep, fi);
      break;
    case eigen_dynamics:
      Outdm->print_newstep(Outdm->outf, lcid, istep, lambda);
      Outdm->print_out_forced(Outdm->outf, lcid, istep, lambda);
      Outdm->print_graphics_forced(Outdm->outgr, lcid, lambda, istep, fi);
      break;
    default:
      print_err("unsupported problem type is required", __FILE__, __LINE__, __func__);
      break;
  }
}



/**
  The function performs all prints requirements  without respect to
  to selected step/time (forced printing).

  @param lcid - load case id
  @param istep - step id
  @param lambda - load coefficient or time (it depends on problem type and solver)
  @param fi - array of additional values at nodes which should be printed.
              It depends on problem and solver type. For example for nonlinear
              statics it contains the internal forces.
  @param fr - array of residual %vector components at nodes which should be printed.
  
  @return The function does not return anything.

  Created 01.2024 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz

*/
void print_step_forced(long lcid, long istep, double lambda, double *fi, double *fr)
{
  switch(Mp->tprob){
    case mat_nonlinear_statics:
    case geom_nonlinear_statics:
    case earth_pressure:
    case forced_dynamics:
    case mech_timedependent_prob:
    case nonlin_floating_subdomain:
    case growing_mech_structure:
      Outdm->print_newstep(Outdm->outf, lcid, istep, lambda);
      Outdm->print_out_forced(Outdm->outf, lcid, istep, lambda);
      Outdm->print_diags_forced(lcid, lambda, istep, fi, fr);
      Outdm->print_graphics_forced(Outdm->outgr, lcid, lambda, istep, fi);
      break;
    default:
      print_err("unsupported problem type is required", __FILE__, __LINE__, __func__);
      break;
  }
}



/**
  The function performs buffer flush of all opened files managed by the outdriver.
  It is usefull for imadiate output of required values and should be called
  at solver after the function print_step is called. 
  
  @return The function does not return anything.

  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz

*/
void print_flush()
{
  long i;
  if (Outdm->outf)
    fflush(Outdm->outf);
  if (Outdm->outgr)
    fflush(Outdm->outgr);
  for (i=0; i<Outdm->ndiag; i++)
  {
    if (Outdm->outdiagf[i])
      fflush(Outdm->outdiagf[i]);
  }
}



/**
  The function performs closing of all opened files managed by the outdriver.
  
  @return The function does not return anything.

  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void print_close()
{
  long i;
  if (Outdm->outf)
    fclose(Outdm->outf);
  Outdm->outf = NULL;
  if (Outdm->outgr)
    fclose(Outdm->outgr);
  Outdm->outgr = NULL;
  for (i=0; i<Outdm->ndiag; i++)
  {
    if (Outdm->outdiagf[i])
      fclose(Outdm->outdiagf[i]);
    Outdm->outdiagf[i] = NULL;
  }
}



/**
  The function exports sets of used elements to the file given by parameter out in GiD format. 

  @param out  - pointer to the opened text file where the output will be produced
  @param idn1 - id of the first node for GiD mesh (default should be idn1 = 1 -> nodes are numbered from 1)
  @param ide1 - id of the first element for GiD mesh (default should be ide1 = 1 -> elements are numbered from 1)

  @return The function does not return anything.

  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void export_gid_mesh(FILE *out, long idn1, long ide1)
{
  long i, print_header, print_coord = 1;

  print_header = 1;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;
    
    switch(Mt->elements[i].te)
    {
      case bar2d:
      case bar3d:
      case beam2d:
      case beam3d:
      case subsoilbeam:
        if (print_header)
	{
          fprintf(out, "MESH \"%ld-beams2\" dimension 3  Elemtype Linear Nnode 2\n", ide1);
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodes(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_element(out, i, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  print_header = 1;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;
    
    switch(Mt->elements[i].te)
    {
      case planequadinterface:
        if (print_header)
	{
          fprintf(out, "MESH \"%ld-plcontact35\" dimension 3  Elemtype Linear Nnode 2\n", ide1);
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodes(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_contact_element(out, i, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  print_header = 1;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;
    
    switch(Mt->elements[i].te)
    {
      case axisymmlqintface:
        if (print_header)
	{
          fprintf(out, "MESH \"%ld-axisymmlqifc66\" dimension 3  Elemtype Linear Nnode 2\n", ide1);
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodes(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_contact_element(out, i, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  print_header = 1;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    switch(Mt->elements[i].te)
    {
      case barq2d:
      case barq3d:
        if (print_header)
	{
          fprintf(out, "MESH \"%ld-beams3\" dimension 3  Elemtype Linear Nnode 3\n", ide1);
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodes(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_element(out, i, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  print_header = 1;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    switch(Mt->elements[i].te)
    {
      case cctel:
      case dktel:
      case dstel:
      case subsoilplatetr:
      case axisymmlt:
      case planeelementlt:
      case planeelementrotlt:
      case shelltrelem:
      case shelltrmelem:
        if (print_header)
	{
          fprintf(out, "MESH \"%ld-trias3\" dimension 3  Elemtype Triangle Nnode 3\n", ide1);
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodes(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_element(out, i, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  print_header = 1;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    switch(Mt->elements[i].te)
    {
      case axisymmqt:
      case planeelementqt:
        if (print_header)
	{
          fprintf(out, "MESH \"%ld-trias6\" dimension 3  Elemtype Triangle Nnode 6\n", ide1);
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodes(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_element(out, i, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  print_header = 1;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    switch(Mt->elements[i].te)
    {
      case axisymmlq:
      case planeelementlq:
      case planeelementrotlq:
      case q4plateel:
      case subsoilplateq:
        if (print_header)
	{
          fprintf(out, "MESH \"%ld-Quads4\" dimension 3  Elemtype Quadrilateral Nnode 4\n", ide1);
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodes(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_element(out, i, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  print_header = 1;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    switch(Mt->elements[i].te)
    {
      case planeelementqq:
      case axisymmqq:
      case axisymmcq:
        if (print_header)
	{
          fprintf(out, "MESH \"%ld-Quads8\" dimension 3  Elemtype Quadrilateral Nnode 8\n", ide1);
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodes(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_element(out, i, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  print_header = 1;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    switch(Mt->elements[i].te)
    {
      case lineartet:
        if (print_header)
	{
          fprintf(out, "MESH \"%ld-Tetras4\" dimension 3  Elemtype Tetrahedra Nnode 4\n", ide1);
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodes(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_element(out, i, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  print_header = 1;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    switch(Mt->elements[i].te)
    {
      case quadrtet:
        if (print_header)
	{
          fprintf(out, "MESH \"%ld-Tetras10\" dimension 3  Elemtype Tetrahedra Nnode 10\n", ide1);
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodes(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_element(out, i, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  print_header = 1;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    switch(Mt->elements[i].te)
    {
      case linearhex:
        if (print_header)
	{
          fprintf(out, "MESH \"%ld-Brick8\" dimension 3  Elemtype Hexahedra Nnode 8\n", ide1);
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodes(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_element(out, i, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  print_header = 1;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    switch(Mt->elements[i].te)
    {
      case quadrhex:
        if (print_header)
	{
          fprintf(out, "MESH \"%ld-Brick20\" dimension 3  Elemtype Hexahedra Nnode 20\n", ide1);
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodes(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_element(out, i, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  print_header = 1;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;
    
    switch(Mt->elements[i].te)
    {
      case hexintface:
        if (print_header)
	{
          fprintf(out, "MESH \"%ld-hexintface\" dimension 3  Elemtype Quadrilateral Nnode 4\n", ide1);
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodes(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_contact_element(out, i, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  print_header = 1;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;
    
    switch(Mt->elements[i].te)
    {
      case linearwed:
        if (print_header)
	{
          fprintf(out, "MESH \"%ld-Wedge6\" dimension 3  Elemtype Prism Nnode 6\n", ide1);
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodes(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_element(out, i, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
}



/**
  The function exports sets of gausspoints of used elements to the file given by parameter 
  out in GiD format. 

  @param out  - pointer to the opened text file where the output will be produced
  @param ide1 - id of the first element for GiD mesh (default is ide1 = 1 -> elements are numbered from 1)

  @return The function does not return anything.

  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void export_gid_gauss_pt(FILE *out, long ide1)
{
  long i, j, k, ii, jj, brk;
  vector gp1, gp2, gp3, gp, w, wt;


  brk = 0;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    switch(Mt->elements[i].te)
    {
      case bar2d:
      case bar3d:
        fprintf(out, "GaussPoints \"Lin_1D\" Elemtype Linear \"%ld-beams2\"\n", ide1);
        fprintf(out, "Number Of Gauss Points: 1\n");
        fprintf(out, "Nodes not included\n");
        fprintf(out, "Natural coordinates: internal\n");
        fprintf(out, "end GaussPoints\n\n");
        brk = 1;
        break;    
      default:
        break;
    }
    if (brk)
      break;
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == planequadinterface)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Linear \"%ld-plcontact35\"\n", planequadinterface, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Pqifc->tnip);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: internal\n");
      fprintf(out, "end GaussPoints\n\n");
      break;    
    }
  }

  brk = 0;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == axisymmlqintface)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Linear \"%ld-axisymmlqifc66\"\n", axisymmlqintface, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Asymlqifc->tnip);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: internal\n");
      fprintf(out, "end GaussPoints\n\n");
      break;    
    }
  }

  brk = 0;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    switch(Mt->elements[i].te)
    {
      case beam2d:
      case beam3d:
      case subsoilbeam:
        fprintf(out, "GaussPoints \"Beam_1D\" Elemtype Linear \"%ld-beams2\"\n", ide1);
        fprintf(out, "Number Of Gauss Points: 2\n");
        fprintf(out, "Nodes included\n");
        fprintf(out, "Natural coordinates: internal\n");
        fprintf(out, "end GaussPoints\n\n");
        brk = 1;
        break;    
      default:
        break;
    }
    if (brk)
      break;
  }

  brk = 0;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    switch(Mt->elements[i].te)
    {
      case barq2d:
      case barq3d:
        fprintf(out, "GaussPoints \"Quad_1D\" Elemtype Linear \"%ld-beams3\"\n", ide1);
        fprintf(out, "Number Of Gauss Points: 3\n");
        fprintf(out, "Nodes not included\n");
        fprintf(out, "Natural coordinates: internal\n");
        fprintf(out, "end GaussPoints\n\n");
        brk = 1;
        break;    
      default:
        break;
    }
    if (brk)
      break;
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == axisymmlt)
    {
      if (Asymlt->nb != 1)
      {
        print_err("%ld block integration scheme is not supported,\n only one block can be used on axisymmlt elements actually", 
                  __FILE__, __LINE__, __func__, Asymlt->nb);
        abort();
      }
      if (Asymlt->intordsm[0][0] != 3)
      {
        print_err("integration order %ld is not supported,\n only intord=3 can be used on axisymmlt elements actually", 
                  __FILE__, __LINE__, __func__, Asymlt->intordsm[0][0]);
        abort();
      }
      fprintf(out, "GaussPoints \"%d\" Elemtype Triangle \"%ld-trias3\"\n", axisymmlt, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Mt->give_totnip(i));
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      reallocv (RSTCKVEC(Asymlt->intordsm[0][0],w));
      reallocv (RSTCKVEC(Asymlt->intordsm[0][0],gp1));
      reallocv (RSTCKVEC(Asymlt->intordsm[0][0],gp2));
      gauss_points_tr (gp1.a,gp2.a,w.a,Asymlt->intordsm[0][0]);
      for (i=0;i<Asymlt->intordsm[0][0];i++)
        fprintf(out, "%le %le\n", gp2[i], gp1[(i+1)%3]);
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if ((Mt->elements[i].te == cctel) || (Mt->elements[i].te == shelltrmelem))
    {
      if (Cct->intordsm[0][0] != 1)
      {
        print_err("integration order %ld is not supported,\n only intord=1 can be used on cct elements actually", 
                  __FILE__, __LINE__, __func__, Cct->intordsm[0][0]);
        abort();
      }
      if (Mt->elements[i].te == cctel)
        fprintf(out, "GaussPoints \"%d\" Elemtype Triangle \"%ld-trias3\"\n", cctel, ide1);
      else
        fprintf(out, "GaussPoints \"%d\" Elemtype Triangle \"%ld-trias3\"\n", shelltrmelem, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Cct->nip[0][0]);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: internal\n");
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if ((Mt->elements[i].te == dktel) || (Mt->elements[i].te == shelltrelem))
    {
      if (Dkt->nb != 1)
      {
        print_err("%ld block integration scheme is not supported,\n only one block can be used on dkt elements actually", 
                  __FILE__, __LINE__, __func__, Dkt->nb);
        abort();
      }
      if (Dkt->intordsm[0][0] != 3)
      {
        print_err("integration order %ld is not supported,\n only intord=3 can be used on dkt or shells actually", 
                  __FILE__, __LINE__, __func__, Dkt->intordsm[0][0]);
        abort();
      }
      if (Mt->elements[i].te == dktel)
        fprintf(out, "GaussPoints \"%d\" Elemtype Triangle \"%ld-trias3\"\n", dktel, ide1);
      else
        fprintf(out, "GaussPoints \"%d\" Elemtype Triangle \"%ld-trias3\"\n", shelltrelem, ide1);
      //      fprintf(out, "Number Of Gauss Points: %ld\n", Mt->give_totnip(i));
      fprintf(out, "Number Of Gauss Points: %ld\n", Dkt->intordsm[0][0]);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      for (ii=0;ii<Dkt->nb;ii++)
      {
        for (jj=0;jj<Dkt->nb;jj++)
        {
          if (Dkt->intordsm[ii][jj]==0)  
            continue;
          reallocv (RSTCKVEC(Dkt->intordsm[ii][jj],w));
          reallocv (RSTCKVEC(Dkt->intordsm[ii][jj],gp1));
          reallocv (RSTCKVEC(Dkt->intordsm[ii][jj],gp2));
          gauss_points_tr (gp1.a,gp2.a,w.a,Dkt->intordsm[ii][jj]);
          for (i=0;i<Dkt->intordsm[ii][jj];i++)
            fprintf(out, "%le %le\n", gp2[i], gp1[(i+1)%3]);
	}
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == subsoilplatetr)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Triangle \"%ld-trias3\"\n", subsoilplatetr, ide1);
      // skip the first integration point in the middle of element
      fprintf(out, "Number Of Gauss Points: %ld\n", Mt->give_totnip(i)-1);

      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      for (ii=0;ii<Spltr->nb;ii++)
      {
        for (jj=0;jj<Spltr->nb;jj++)
        {
          if (Spltr->intordsm[ii][jj]==0)  
            continue;
          reallocv (RSTCKVEC(Spltr->intordsm[ii][jj],w));
          reallocv (RSTCKVEC(Spltr->intordsm[ii][jj],gp1));
          reallocv (RSTCKVEC(Spltr->intordsm[ii][jj],gp2));
          gauss_points_tr (gp1.a,gp2.a,w.a,Spltr->intordsm[ii][jj]);
          // skip the first integration point in the middle of element
          for (i=1;i<Spltr->intordsm[ii][jj];i++)
            fprintf(out, "%le %le\n", gp1[i], gp2[i]);
	}
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == planeelementlt)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Triangle\n", planeelementlt);
      fprintf(out, "Number Of Gauss Points: %ld\n", Mt->give_totnip(i));
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: internal\n");
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == planeelementrotlt)
    {
      if (Perlt->intordsm[0][0] != 3)
      {
        print_err("integration order %ld is not supported,\n only intord=3 can be used on planeelementrotlt or shells actually", 
                  __FILE__, __LINE__, __func__, Perlt->intordsm[0][0]);
        abort();
      }
      fprintf(out, "GaussPoints \"%d\" Elemtype Triangle \"%ld-trias3\"\n", planeelementrotlt, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Perlt->nip[0][0]);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      reallocv (RSTCKVEC(Perlt->intordsm[0][0],w));
      reallocv (RSTCKVEC(Perlt->intordsm[0][0],gp1));
      reallocv (RSTCKVEC(Perlt->intordsm[0][0],gp2));
      gauss_points_tr (gp1.a,gp2.a,w.a,Perlt->intordsm[0][0]);
      for (i=0;i<Perlt->intordsm[0][0];i++)
        fprintf(out, "%le %le\n", gp2[i], gp1[(i+1)%3]);
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == planeelementqt)
    {
      if (Peqt->intordsm[0][0] != 3)
      {
        print_err("integration order %ld is not supported,\n only intord=3 can be used on planeelementqt actually", 
                  __FILE__, __LINE__, __func__, Peqt->intordsm[0][0]);
        abort();
      }
      fprintf(out, "GaussPoints \"%d\" Elemtype Triangle \"%ld-trias6\"\n", planeelementqt, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Peqt->nip[0][0]);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      reallocv (RSTCKVEC(Peqt->intordsm[0][0],w));
      reallocv (RSTCKVEC(Peqt->intordsm[0][0],gp1));
      reallocv (RSTCKVEC(Peqt->intordsm[0][0],gp2));
      gauss_points_tr (gp1.a,gp2.a,w.a,Peqt->intordsm[0][0]);
      for (i=0;i<Peqt->intordsm[0][0];i++)
        fprintf(out, "%le %le\n", gp2[i], gp1[(i+1)%3]);
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }


  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == axisymmqt)
    {
      if (Asymqt->intordsm[0][0] != 3)
      {
        print_err("integration order %ld is not supported,\n only intord=3 can be used on axisymmqt actually", 
                  __FILE__, __LINE__, __func__, Asymqt->intordsm[0][0]);
        abort();
      }
      fprintf(out, "GaussPoints \"%d\" Elemtype Triangle \"%ld-trias6\"\n", axisymmqt, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Asymqt->nip[0][0]);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      reallocv (RSTCKVEC(Asymqt->intordsm[0][0],w));
      reallocv (RSTCKVEC(Asymqt->intordsm[0][0],gp1));
      reallocv (RSTCKVEC(Asymqt->intordsm[0][0],gp2));
      gauss_points_tr (gp1.a,gp2.a,w.a,Asymqt->intordsm[0][0]);
      for (i=0;i<Asymqt->intordsm[0][0];i++)
        fprintf(out, "%le %le\n", gp2[i], gp1[(i+1)%3]);
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == axisymmlq)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Quadrilateral \"%ld-Quads4\"\n", axisymmlq, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", 4L);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      reallocv (RSTCKVEC(Asymlq->intordsm[0][0],w));
      reallocv (RSTCKVEC(Asymlq->intordsm[0][0],gp));
      gauss_points (gp.a,w.a,Asymlq->intordsm[0][0]);
      fprintf(out, "%le %le\n", -gp[0], -gp[0]);
      fprintf(out, "%le %le\n", -gp[0],  gp[0]);
      fprintf(out, "%le %le\n",  gp[0], -gp[0]);
      fprintf(out, "%le %le\n",  gp[0],  gp[0]);
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == planeelementlq)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Quadrilateral \"%ld-Quads4\"\n", planeelementlq, ide1);
      fprintf(out, "Number Of Gauss Points: 4\n");
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      for (ii=0;ii<1;ii++)
      {
        for (jj=0;jj<1;jj++)
        {
          if (Pelq->intordsm[ii][jj]==0)  
            continue;
          reallocv (RSTCKVEC(Pelq->intordsm[ii][jj],w));
          reallocv (RSTCKVEC(Pelq->intordsm[ii][jj],gp));
          gauss_points (gp.a,w.a,Pelq->intordsm[ii][jj]);
          for (i=0;i<Pelq->intordsm[ii][jj];i++)
            for (j=0; j<Pelq->intordsm[ii][jj]; j++)
              fprintf(out, "%le %le\n", -gp[i], -gp[j]);
        }
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }
  
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == planeelementrotlq)
    {
      if (Perlq->intordsm[0][0] != 3)
      {
        print_err("integration order %ld is not supported,\n only intord=3 can be used on planeelementrotlq actually", 
                  __FILE__, __LINE__, __func__, Perlq->intordsm[0][0]);
        abort();
      }
      fprintf(out, "GaussPoints \"%d\" Elemtype Quadrilateral \"%ld-Quads4\"\n", planeelementrotlq, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Perlq->nip[0][0]);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      reallocv (RSTCKVEC(Perlq->intordsm[0][0],w));
      reallocv (RSTCKVEC(Perlq->intordsm[0][0],gp));
      gauss_points (gp.a,w.a,Perlq->intordsm[0][0]);
      for (i=0;i<Perlq->intordsm[0][0];i++)
        for (j=0; j<Perlq->intordsm[0][0]; j++)
          fprintf(out, "%le %le\n", -gp[i], -gp[j]);
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == planeelementqq)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Quadrilateral \"%ld-Quads8\"\n", planeelementqq, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Mt->give_totnip(i));
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      for (ii=0;ii<Peqq->nb;ii++)
      {
        for (jj=0;jj<Peqq->nb;jj++)
        {
          if (Peqq->intordsm[ii][jj]==0)  
            continue;
          reallocv (RSTCKVEC(Peqq->intordsm[ii][jj],w));
          reallocv (RSTCKVEC(Peqq->intordsm[ii][jj],gp));
          gauss_points (gp.a,w.a,Peqq->intordsm[ii][jj]);
          for (i=0;i<Peqq->intordsm[ii][jj];i++)
            for (j=0; j<Peqq->intordsm[ii][jj]; j++)
              fprintf(out, "%le %le\n", -gp[i], -gp[j]);
	}
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == axisymmqq)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Quadrilateral \"%ld-Quads8\"\n", axisymmqq, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Mt->give_totnip(i));
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      for (ii=0;ii<Asymqq->nb;ii++)
      {
        for (jj=0;jj<Asymqq->nb;jj++)
        {
          if (Asymqq->intordsm[ii][jj]==0)  
            continue;
          reallocv (RSTCKVEC(Asymqq->intordsm[ii][jj],w));
          reallocv (RSTCKVEC(Asymqq->intordsm[ii][jj],gp));
          gauss_points (gp.a,w.a,Asymqq->intordsm[ii][jj]);
          for (i=0;i<Asymqq->intordsm[ii][jj];i++)
            for (j=0; j<Asymqq->intordsm[ii][jj]; j++)
              fprintf(out, "%le %le\n", -gp[i], -gp[j]);
	}
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;
    if (Mt->elements[i].te == axisymmcq)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Quadrilateral \"%ld-Quads8\"\n", axisymmcq, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Mt->give_totnip(i));
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      for (ii=0;ii<Asymcq->nb;ii++)
      {
        for (jj=0;jj<Asymcq->nb;jj++)
        {
          if (Asymcq->intordsm[ii][jj]==0)  
            continue;
          reallocv (RSTCKVEC(Asymcq->intordsm[ii][jj],w));
          reallocv (RSTCKVEC(Asymcq->intordsm[ii][jj],gp));
          gauss_points (gp.a,w.a,Asymcq->intordsm[ii][jj]);
          for (i=0;i<Asymcq->intordsm[ii][jj];i++)
            for (j=0; j<Asymcq->intordsm[ii][jj]; j++)
              fprintf(out, "%le %le\n", -gp[i], -gp[j]);
	}
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == lineartet)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype TetraHedra \"%ld-Tetras4\"\n", lineartet, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Mt->give_totnip(i));
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      for (ii=0;ii<Ltet->nb;ii++)
      {
        for (jj=0;jj<Ltet->nb;jj++)
        {
          if (Ltet->intordsm[ii][jj]==0)  
            continue;
          reallocv (RSTCKVEC(Ltet->intordsm[ii][jj],w));
          reallocv (RSTCKVEC(Ltet->intordsm[ii][jj],gp1));
          reallocv (RSTCKVEC(Ltet->intordsm[ii][jj],gp2));
          reallocv (RSTCKVEC(Ltet->intordsm[ii][jj],gp3));
          gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,Ltet->intordsm[ii][jj]);
          for (i=0;i<Ltet->intordsm[ii][jj];i++)
            fprintf(out, "%le %le %le\n", gp1[i], gp2[i], gp3[i]);
	}
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == quadrtet)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype TetraHedra \"%ld-Tetras10\"\n", quadrtet, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Mt->give_totnip(i));
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      for (ii=0;ii<Qtet->nb;ii++)
      {
        for (jj=0;jj<Qtet->nb;jj++)
        {
          if (Qtet->intordsm[ii][jj]==0)  
            continue;
          reallocv (RSTCKVEC(Qtet->intordsm[ii][jj],w));
          reallocv (RSTCKVEC(Qtet->intordsm[ii][jj],gp1));
          reallocv (RSTCKVEC(Qtet->intordsm[ii][jj],gp2));
          reallocv (RSTCKVEC(Qtet->intordsm[ii][jj],gp3));
          gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,Qtet->intordsm[ii][jj]);
          for (i=0;i<Qtet->intordsm[ii][jj];i++)
            fprintf(out, "%le %le %le\n", gp1[i], gp2[i], gp3[i]);
	}
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == linearhex)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Hexahedra \"%ld-Brick8\"\n", linearhex, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Mt->give_totnip(i));
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      for (ii=0;ii<Lhex->nb;ii++)
      {
        for (jj=0;jj<Lhex->nb;jj++)
        {
          if (Lhex->intordsm[ii][jj]==0)  
            continue;
          reallocv (RSTCKVEC(Lhex->intordsm[ii][jj],w));
          reallocv (RSTCKVEC(Lhex->intordsm[ii][jj],gp));
          gauss_points(gp.a,w.a,Lhex->intordsm[ii][jj]);
          for (i=0;i<Lhex->intordsm[ii][jj];i++)
          {
            for (j=0;j<Lhex->intordsm[ii][jj];j++)
            {
              for (k=0;k<Lhex->intordsm[ii][jj];k++)
                fprintf(out, "%le %le %le\n", -gp[i], -gp[j], -gp[k]);
	    }
	  }
	}
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == quadrhex)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Hexahedra \"%ld-Brick20\"\n", quadrhex, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Mt->give_totnip(i));
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      for (ii=0;ii<Qhex->nb;ii++)
      {
        for (jj=0;jj<Qhex->nb;jj++)
        {
          if (Qhex->intordsm[ii][jj]==0)  
            continue;
          reallocv (RSTCKVEC(Qhex->intordsm[ii][jj],w));
          reallocv (RSTCKVEC(Qhex->intordsm[ii][jj],gp));
          gauss_points(gp.a,w.a,Qhex->intordsm[ii][jj]);
          for (i=0;i<Qhex->intordsm[ii][jj];i++)
          {
            for (j=0;j<Qhex->intordsm[ii][jj];j++)
            {
              for (k=0;k<Qhex->intordsm[ii][jj];k++)
                fprintf(out, "%le %le %le\n", -gp[i], -gp[j], -gp[k]);
	    }
	  }
	}
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == hexintface)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Quadrilateral \"%ld-hexintface\"\n", hexintface, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Hexifc->tnip);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: internal\n");
      fprintf(out, "end GaussPoints\n\n");
      break;    
    }
  }

  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    if (Mt->elements[i].te == linearwed)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Quadrilateral \"%ld-Wedge6\"\n", linearwed, ide1);
      fprintf(out, "Number Of Gauss Points: %ld\n", Lwed->tnip);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      fprintf(out, "end GaussPoints\n\n");
      for (ii=0;ii<Lwed->nb;ii++){
        for (jj=0;jj<Lwed->nb;jj++){
          if (Lwed->intordsmt[ii][jj]==0)  continue;

          reallocv (RSTCKVEC(Lwed->intordsmz[ii][jj],w));
          reallocv (RSTCKVEC(Lwed->intordsmz[ii][jj],gp));
          reallocv (RSTCKVEC(Lwed->intordsmt[ii][jj],wt));
          reallocv (RSTCKVEC(Lwed->intordsmt[ii][jj],gp1));
          reallocv (RSTCKVEC(Lwed->intordsmt[ii][jj],gp2));
          
          gauss_points (gp.a,w.a,Lwed->intordsmz[ii][jj]);
          gauss_points_tr (gp1.a,gp2.a,wt.a,Lwed->intordsmt[ii][jj]);
          
          for (i=0;i<Lwed->intordsmt[ii][jj];i++){
            for (k=0;k<Lwed->intordsmz[ii][jj];k++){
              fprintf(out, "%le %le %le\n", gp2[i], gp1[(i+1)%3], 0.5*(1.0-gp[i]));
            }
          }
        }
      }
      break;    
    }
  }
}



/**
  The function exports from 3D brick element mesh the 2D plane element cuts to the file given by parameter out in GiD format. 

  @param out  - pointer to the opened text file where the output will be produced
  @param icut - index of the cut - will be used in the element property
  @param idn1 - id of the first node for GiD mesh (default should be idn1 = 1 -> nodes are numbered from 1)
  @param ide1 - id of the first element for GiD mesh (default should be ide1 = 1 -> elements are numbered from 1)

  @return The function does not return anything.
  
  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void export_gid_2dmesh(FILE *out, long icut, long idn1, long ide1)
{
  long i, print_header, print_coord = 1;
  long range = Mt->ne/icut;

  print_header = 1;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    switch(Mt->elements[i].te)
    {
      case linearhex:
        if (print_header)
	{
          fprintf(out, "MESH Quads4 dimension 3  Elemtype Quadrilateral Nnode 4\n");
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodes(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_2delement(out, i, 0, 4, i/range, 0, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  for (i=Mt->ne-range; i < Mt->ne; i++)
  {
    write_gid_2delement(out, i, 4, 8, icut, range, idn1, ide1);
  }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
}



/**
  The function exports nodes to the file given by parameter out in GiD format. 

  @param out  - pointer to the opened text file where the output will be produced
  @param idn1 - id of the first node for GiD mesh (default should be idn1 = 1 -> nodes are numbered from 1)

  @return The function does not return anything.

  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void write_gid_nodes(FILE *out, long idn1)
{
  long i;
  fprintf(out, "Coordinates\n");
  for (i=0; i<Mt->nn; i++)
    fprintf(out, "%ld %e %e %e\n", i+idn1, Gtm->gnodes[i].x, Gtm->gnodes[i].y, Gtm->gnodes[i].z);
  fprintf(out, "end Coordinates\n");
}



/**
  The function exports element to the file given by parameter out in GiD format. 

  @param out  - pointer to the opened text file where the output will be produced
  @param i    - element id
  @param idn1 - id of the first node for GiD mesh (default should be idn1 = 1 -> nodes are numbered from 1)
  @param ide1 - id of the first element for GiD mesh (default should be ide1 = 1 -> elements are numbered from 1)

  @return The function does not return anything.

  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void write_gid_element(FILE *out, long i, long idn1, long ide1)
{
  long j, matid;
  fprintf(out, "%ld ", i+ide1);
  for (j=0; j<Gtm->gelements[i].give_nne(); j++)
    fprintf(out, "%ld ", Gtm->gelements[i].nodes[j]+idn1);
  matid = mattype_kwdset.get_ordid(Mt->elements[i].tm[0])+1;
  fprintf(out, "%ld%ld\n", matid, Mt->elements[i].idm[0]+1);
//  fprintf(out, "%ld\n", Mt->elements[i].idm[0]+1);
/*  for (j=0; j < Mt->elements[i].nm; j++)
    fprintf(out, "%d%ld", Mt->elements[i].tm[j],Mt->elements[i].idm[j]+1);
  fprintf(out, "\n");*/
}



/**
  The function exports contact element to the file given by parameter out in GiD format. 
  Nodes from half of element are exported only.

  @param out  - pointer to the opened text file where the output will be produced
  @param i    - element id
  @param idn1 - id of the first node for GiD mesh (default should be idn1 = 1 -> nodes are numbered from 1)
  @param ide1 - id of the first element for GiD mesh (default should be ide1 = 1 -> elements are numbered from 1)

  @return The function does not return anything.

  Created 2013 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void write_gid_contact_element(FILE *out, long i, long idn1, long ide1)
{
  long j;
  fprintf(out, "%ld ", i+ide1);
  for (j=0; j<Gtm->gelements[i].give_nne()/2; j++)
    fprintf(out, "%ld ", Gtm->gelements[i].nodes[j]+idn1);
  fprintf(out, "%ld\n", Mt->elements[i].idm[0]+1);
}



/**
  The function exports 3D brick element as 2D plane element to the file given by parameter out in GiD format. 
  Purpose of the function is to perform cuts of the domain which consists of the regular mesh
  of brick elements.

  @param out  - pointer to the opened text file where the output will be produced
  @param i    - brick element id
  @param id1  - brick element index of the first node for plane element
  @param nne  - brick element index of the last node for plane element
  @param icut - index of the cut - will be used in the element property
  @param di   - the start index of the plane elements
  @param idn1 - id of the first node for GiD mesh (default should be idn1 = 1 -> nodes are numbered from 1)
  @param ide1 - id of the first element for GiD mesh (default should be ide1 = 1 -> elements are numbered from 1)

  @return The function does not return anything.

  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void write_gid_2delement(FILE *out, long i, long id1, long nne, long icut, long di, long idn1, long ide1)
{
  long j;
  fprintf(out, "%ld ", i+di+ide1);
  for (j=id1; j<nne; j++)
    fprintf(out, "%ld ", Gtm->gelements[i].nodes[j]+idn1);
  fprintf(out, "%ld%ld\n", icut+1, Mt->elements[i].idm[0]+1);
/*  for (j=0; j < Mt->elements[i].nm; j++)
    fprintf(out, "%ld%ld", Mt->elements[i].tm[j],Mt->elements[i].idm[j]+1);
  fprintf(out, "\n");*/
}



/**
  The function prints set of elements of required types to the GiD mesh file.

  @param[in,out] out - pointer to the opened GiD mesh file
  @param[in] reqet - set of required element types that will be printed
  @param[in] header - string with element set header quoting section of elements of one common type
  @param[in] idn1 - shift index of the nodal numbers
  @param[in] ide1 - shift index of the element numbers

  @return The function does not return anything, but it may change the content of the file referred by the argumnet out.

  Created by Tomas Koudelka, 10.2023
*/
void write_gid_elem_set(FILE *out, const std::set<elemtype> &reqet, const char *header, long idn1, long ide1)
{
  long i;
  
  bool pheader = true;
  for (i=0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;
    if (reqet.find(Mt->elements[i].te) != reqet.end()){
      if (pheader){
        fprintf(out, "%s\n", header);
        fprintf(out, "Elements\n");
        pheader = false;
      }
      write_gid_element(out, i, idn1, ide1);
    }
  }
  if (pheader == false)
    fprintf(out, "End Elements\n");
}



/**
  The function exports nodes and sets of used elements to the file given by parameter out in the GiD mesh format. 

  @param out  - pointer to the opened text file where the output will be produced
  @param idn1 - id of the first node for GiD mesh (default should be idn1 = 1 -> nodes are numbered from 1)
  @param ide1 - id of the first element for GiD mesh (default should be ide1 = 1 -> elements are numbered from 1)

  @return The function does not return anything, but it changes the content of file referred by the argumnet out.

  Created by Tomas Koudelka, 10.2023, tomas.koudelka@fsv.cvut.cz
*/
void write_gid_mesh(FILE *out, long idn1, long ide1)
{
  const long hs = 1025;
  char header[hs];
  write_gid_nodes(out, idn1);
 
  snprintf(header, hs, "MESH \"%ld-beams2\" dimension 3  Elemtype Linear Nnode 2\n", ide1);
  write_gid_elem_set(out, {bar2d, bar3d, beam2d, beam3d, subsoilbeam}, header, idn1, ide1);

  snprintf(header, hs, "MESH \"%ld-plcontact2l\" dimension 3  Elemtype Linear Nnode 2\n", ide1);
  write_gid_elem_set(out, {planequadinterface}, header, idn1, ide1);

  snprintf(header, hs, "MESH \"%ld-axisymmlqifc66\" dimension 3  Elemtype Linear Nnode 2\n", ide1);
  write_gid_elem_set(out, {axisymmlqintface}, header, idn1, ide1);

  snprintf(header, hs, "MESH \"%ld-beams3\" dimension 3  Elemtype Linear Nnode 3\n", ide1);
  write_gid_elem_set(out, {barq2d, barq3d}, header, idn1, ide1);

  snprintf(header, hs, "MESH \"%ld-trias3\" dimension 3  Elemtype Triangle Nnode 3\n", ide1);
  write_gid_elem_set(out, {cctel, dktel, dstel, subsoilplatetr, axisymmlt, planeelementlt,
                           planeelementrotlt, shelltrelem, shelltrmelem}, header, idn1, ide1);

  snprintf(header, hs, "MESH \"%ld-trias6\" dimension 3  Elemtype Triangle Nnode 6\n", ide1);
  write_gid_elem_set(out, {axisymmqt, planeelementqt}, header, idn1, ide1);

  snprintf(header, hs, "MESH \"%ld-Quads4\" dimension 3  Elemtype Quadrilateral Nnode 4\n", ide1);
  write_gid_elem_set(out, {axisymmlq, planeelementlq, planeelementrotlq, q4plateel, subsoilplateq}, header, idn1, ide1);

  snprintf(header, hs, "MESH \"%ld-Quads8\" dimension 3  Elemtype Quadrilateral Nnode 8\n", ide1);
  write_gid_elem_set(out, {planeelementqq, axisymmqq, axisymmcq}, header, idn1, ide1);

  snprintf(header, hs, "MESH \"%ld-Tetras4\" dimension 3  Elemtype Tetrahedra Nnode 4\n", ide1);
  write_gid_elem_set(out, {lineartet}, header, idn1, ide1);

  snprintf(header, hs, "MESH \"%ld-Tetras10\" dimension 3  Elemtype Tetrahedra Nnode 10\n", ide1);
  write_gid_elem_set(out, {quadrtet}, header, idn1, ide1);

  snprintf(header, hs, "MESH \"%ld-Brick8\" dimension 3  Elemtype Hexahedra Nnode 8\n", ide1);
  write_gid_elem_set(out, {linearhex}, header, idn1, ide1);

  snprintf(header, hs, "MESH \"%ld-Brick20\" dimension 3  Elemtype Hexahedra Nnode 20\n", ide1);
  write_gid_elem_set(out, {quadrhex}, header, idn1, ide1);

  snprintf(header, hs, "MESH \"%ld-hexintface\" dimension 3  Elemtype Quadrilateral Nnode 4\n", ide1);
  write_gid_elem_set(out, {hexintface}, header, idn1, ide1);

  snprintf(header, hs, "MESH \"%ld-Wedge6\" dimension 3  Elemtype Prism Nnode 6\n", ide1);
  write_gid_elem_set(out, {linearwed}, header, idn1, ide1);
}



/**
  The function writes displacement vector for all nodes to the file given by parameter out in GiD format. 
  The results are printed for all nodes with no dependency on the outdriver selection.

  @param out - pointer to the opened text file where the output will be produced
  @param lcid - load case id
  @param desclcid - string with description of loadcase

  @return The function does not return anything.

  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void write_gid_displ(FILE *out, long lcid, const char *desclcid)
{
  long i, j;
  vector l,g;
  matrix tm(3,3);
  long ndofn, transf;
  long addmacrol = 0;
  long maxncomp = Mt->max_ncompstr;
  strastrestate ssst = guess_ssst(maxncomp);
  vector mstr(ASTCKVEC(Mt->max_ncompstr));
  matrix mstrt(ASTCKMAT(3,3));
  vector ncoord(ASTCKVEC(3));
  vector ml(ASTCKVEC(3));

  fprintf(out, "\nResult \"Displacements\" \"%ld\" %s Vector OnNodes\n", lcid, desclcid);
  
  fprintf(out, "ComponentNames \"r_1\" \"r_2\" ");
  if (Mt->give_maxndofn() > 2)
    fprintf(out, "\"r_3\"");
  fprintf(out, "\nValues\n");

  // prepare macrolevel component of displacements for homogenization problems
  if ((Mp->homog == 3) || (Mp->homog == 4) || (Mp->homog == 9))
  {
    addmacrol = 1;
    for (j=0; j<maxncomp; j++)
      mstr(j) = macrostrains(lcid, j);
    vector_tensor(mstr, mstrt, ssst, strain);
  }

  for (i=0; i < Mt->nn; i++)
  {    if (Gtm->lnso[i]==0)
      continue;
    if (Outdm->nog.selndisp.presence_id(i) == 0)
      continue;

    ndofn = Mt->give_ndofn (i);
    reallocv (RSTCKVEC(ndofn, l));
    noddispl (lcid, l.a, i);

    if (addmacrol)
    {
      // in homogenization problems (macrostrain(4) and macrostress(3) approaches), 
      // calculated displacements represent fluctuation component of the total displacements
      Mt->give_nodal_coord(i, ncoord);
      mxv(mstrt, ncoord, ml);  // macro level component of displacements ml_i = mstrains_{ij}*x_j
      for (j=0; j<min2(ndofn, 3); j++)
        l(j) += ml(j);
    }

    transf = Mt->nodes[i].transf;
    if (transf>0){

      //fprintf (Outm,"\n uzel %ld   %ld",i,Mt->nodes[i].transf);
      
      reallocv(RSTCKVEC(ndofn,g));
      reallocm(RSTCKMAT(ndofn,ndofn,tm));
      nullm(tm);
      
      tm[0][0]=Mt->nodes[i].e1[0];
      tm[1][0]=Mt->nodes[i].e1[1];
      if (transf == 3)
        tm[2][0]=Mt->nodes[i].e1[2];
      
      tm[0][1]=Mt->nodes[i].e2[0];
      tm[1][1]=Mt->nodes[i].e2[1];
      if (transf == 3)
        tm[2][1]=Mt->nodes[i].e2[2];
      
      if (transf == 3)
      {
        tm[0][2]=Mt->nodes[i].e3[0];
        tm[1][2]=Mt->nodes[i].e3[1];
        tm[2][2]=Mt->nodes[i].e3[2];
      }
      else
      {
        if (ndofn == 3) // 2D beam
          tm[2][2]=1.0;
      }
      
      mxv(tm,l,g);
      copyv (g,l);
    }
    
    fprintf(out, "%ld", i+Outdm->idn1);
    for (j = 0; j < ndofn; j++)
    {
      if (Outdm->nog.selndisp.presence_id(Outdm->nog.seldisp,i,j))
        // for the deflection on plates uncomment following line and comment
        // and select just one appropriate component of displacement vector
        //  fprintf(out, " 0.0 0.0 % e", l[j]);
        fprintf(out, " % e", l[j]);
      else
        fprintf(out, " % e", 0.0);
    }
    fprintf(out, "\n");
  }
  fprintf(out, "End Values\n");
}



/**
  The function writes vector of forces for all nodes to the file given by parameter out in GiD format. 
  The results are printed for all nodes with no dependency on the outdriver selection.

  @param out - pointer to the opened text file where the output will be produced
  @param lcid - load case id
  @param desclcid - string with description of loadcase
  @param veclabel - string with label of printed force %vector
  @param ifor - vector of nodal forces
  @param print_react - flag for printing reactions at DOFs with the prescribed displacements (true)
                       or zero value (false) 

  @return The function does not return anything.

  created 06.2007 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
  modified 04.2016 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void write_gid_nforces(FILE *out, long lcid, const char *desclcid, const char *veclabel, double *ifor, bool print_react)
{
  long i, j, ndof, nc;
  vector f,g;
  matrix tm(3,3);
  long ndofn, transf, print_mom = 0;

  fprintf(out, "\nResult \"%s\" \"%ld\" %s Vector OnNodes\n", veclabel, lcid, desclcid);
  fprintf(out, "ComponentNames \"f_1\" \"f_2\" \"f_3\"\nValues\n");
  for (i=0; i < Mt->nn; i++)
  {
    if (Gtm->lnso[i]==0)
      continue;
    if (Outdm->nog.selnforce.presence_id(i) == 0)
      continue;

    ndof = Mt->give_ndofn(i);
    reallocv (ndof, f);
    nodforce (ifor, i, f, print_react);

    transf = Mt->nodes[i].transf;
    if (transf>0){

      //fprintf (Outm,"\n uzel %ld   %ld",i,Mt->nodes[i].transf);

      ndofn = Mt->give_ndofn (i);
      reallocv(ndofn, g);
      reallocm(ndofn,ndofn,tm);
      fillm(0.0, tm);
      
      tm[0][0]=Mt->nodes[i].e1[0];
      tm[1][0]=Mt->nodes[i].e1[1];
      if (transf == 3)
        tm[2][0]=Mt->nodes[i].e1[2];
      
      tm[0][1]=Mt->nodes[i].e2[0];
      tm[1][1]=Mt->nodes[i].e2[1];
      if (transf == 3)
        tm[2][1]=Mt->nodes[i].e2[2];
      
      if (transf == 3)
      {
        tm[0][2]=Mt->nodes[i].e3[0];
        tm[1][2]=Mt->nodes[i].e3[1];
        tm[2][2]=Mt->nodes[i].e3[2];
      }
      else
      {
        if (ndofn == 3) // 2D beam
          tm[2][2]=1.0;
      }
      
      mxv(tm,f,g);
      copyv (g,f);
    }
    
    fprintf(out, "%ld", i+Outdm->idn1);
    if (ndof > 3)
    {
      nc = 3;
      print_mom = 1;
    }
    else 
      nc = ndof;
    for (j = 0; j < nc; j++)
    {
      if (Outdm->nog.selnforce.presence_id(Outdm->nog.seldisp,i,j))
        fprintf(out, " % e", f[j]);
      else
        fprintf(out, " 0.0");
    }

    for (j = nc; j < 3; j++)
      fprintf(out, " 0.0");
    fprintf(out, "\n");
  }
  fprintf(out, "End Values\n");

  if (print_mom == 0)
    return;

  fprintf(out, "\nResult \"Moments\" \"%ld\" %s Vector OnNodes\n", lcid, desclcid);
  fprintf(out, "ComponentNames \"m_1\" \"m_2\" \"m_3\"\nValues\n");
  for (i=0; i < Mt->nn; i++)
  {
    if (Gtm->lnso[i]==0)
      continue;
    if (Outdm->nog.selnforce.presence_id(i) == 0)
      continue;

    ndof = Mt->give_ndofn(i);
    if (ndof < 4)
      continue;

    reallocv (ndof, f);
    nodforce (ifor, i, f);

    transf = Mt->nodes[i].transf;
    if (transf>0){

      //fprintf (Outm,"\n uzel %ld   %ld",i,Mt->nodes[i].transf);

      ndofn = Mt->give_ndofn (i);
      reallocv(ndofn, g);
      reallocm(ndofn,ndofn,tm);
      fillm(0.0, tm);
      
      tm[0][0]=Mt->nodes[i].e1[0];
      tm[1][0]=Mt->nodes[i].e1[1];
      if (transf == 3)
        tm[2][0]=Mt->nodes[i].e1[2];
      
      tm[0][1]=Mt->nodes[i].e2[0];
      tm[1][1]=Mt->nodes[i].e2[1];
      if (transf == 3)
        tm[2][1]=Mt->nodes[i].e2[2];
      
      if (transf == 3)
      {
        tm[0][2]=Mt->nodes[i].e3[0];
        tm[1][2]=Mt->nodes[i].e3[1];
        tm[2][2]=Mt->nodes[i].e3[2];
      }
      else
      {
        if (ndofn == 3) // 2D beam
          tm[2][2]=1.0;
      }
      
      mxv(tm,f,g);
      copyv (g,f);
    }
    
    fprintf(out, "%ld", i+Outdm->idn1);

    for (j = 3; j < ndof; j++)
    {
      if (Outdm->nog.selnforce.presence_id(Outdm->nog.seldisp,i,j))
        fprintf(out, " % e", f[j]);
      else
        fprintf(out, " 0.0");
    }

    for (j = ndof; j < 6; j++)
      fprintf(out, " 0.0");
    fprintf(out, "\n");
  }
  fprintf(out, "End Values\n");
}



/**
  The function writes a scalar value given by parameter scal on all nodes to the file given by 
  parameter out in GiD format.

  @param out - pointer to the opened text file where the output will be produced
  @param scal - specifies type of required scalar quantity (strain/stres/other)
  @param lcid - load case id
  @param dir - specifies which component of the quantity array will be printed
  @param desclcid - string with description of loadcase

  @return The function does not return anything.

  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void write_gid_nodscalar(FILE *out, strastre scal, long lcid, long dir, const char *desclcid)
{
  const char *sig = "";
  long i, ncompstr, ncompother;
  double aux;
  long ir;
  matrix t;
  vector str, p;
  
  if(Outdm->nlcs > 0)
  {
    reallocv(RSTCKVEC(3, p));
    reallocm(RSTCKMAT(3, 3, t));
  }

  if (scal == stress)
  {
    if (Mm->max_ncompstrn < 6)
    {
      switch (dir)
      {
        case 0:
          sig = "sig_n_1";
          break;
        case 1:
          sig = "sig_n_2";
          break;
        case 2:
          sig = "sig_n_3";
          break;
        case 3:
          sig = "sig_n_4";
          break;
        case 4:
          sig = "sig_n_5";
          break;
        default:
          print_err("unknown direction of stresses is required", __FILE__, __LINE__, __func__);
      }
    }
    if (Mm->max_ncompstrn == 6)
    {
      switch (dir)
      {
        case 0:
          sig = "sig_n_x";
          break;
        case 1:
          sig = "sig_n_y";
          break;
        case 2:
          sig = "sig_n_z";
          break;
        case 3:
          sig = "tau_n_yz";
          break;
        case 4:
          sig = "tau_n_xz";
          break;
        case 5:
          sig = "tau_n_xy";
          break;
        default:
          print_err("unknown direction of stresses is required", __FILE__, __LINE__, __func__);
      }
    }
  }
  if (scal == strain)
  {
    if (Mm->max_ncompstrn < 6)
    {
      switch (dir)
      {
        case 0:
          sig = "eps_n_1";
          break;
        case 1:
          sig = "eps_n_2";
          break;
        case 2:
          sig = "eps_n_3";
          break;
        case 3:
          sig = "eps_n_4";
          break;
        case 4:
          sig = "eps_n_5";
          break;
        default:
          print_err("unknown direction of strains is required", __FILE__, __LINE__, __func__);
      }
    }
    if (Mm->max_ncompstrn == 6)
    {
      switch (dir)
      {
        case 0:
          sig = "eps_n_x";
          break;
        case 1:
          sig = "eps_n_y";
          break;
        case 2:
          sig = "eps_n_z";
          break;
        case 3:
          sig = "eps_n_yz";
          break;
        case 4:
          sig = "eps_n_xz";
          break;
        case 5:
          sig = "eps_n_xy";
          break;
        default:
          print_err("unknown direction of strains is required", __FILE__, __LINE__, __func__);
      }
    }
  }
  if (scal == pstress)
  {
    switch (dir)
    {
      case 0:
        sig = "psig_n_1";
        break;
      case 1:
        sig = "psig_n_2";
        break;
      case 2:
        sig = "psig_n_3";
        break;
      case 3:
        sig = "tau_n_max";
        break;
      default:
        print_err("unknown direction of principal stresses is required", __FILE__, __LINE__, __func__);
    }
  }
  if (scal == pstrain)
  {
    switch (dir)
    {
      case 0:
        sig = "peps_n_1";
        break;
      case 1:
        sig = "peps_n_2";
        break;
      case 3:
        sig = "peps_n_3";
        break;
      default:
        print_err("unknown direction of principal strains is required", __FILE__, __LINE__, __func__);
    }
  }
  if (scal == other)
    fprintf(out, "\nResult \"other_n_%ld\" \"%ld\" %s Scalar OnNodes\n", dir+1, lcid, desclcid);
  else     
    fprintf(out, "\nResult \"%s\" \"%ld\" %s Scalar OnNodes\n", sig, lcid, desclcid);
  fprintf(out, "Values\n");
  switch (scal)
  {
    case stress:
      for (i = 0; i < Mt->nn; i++)
      {
        if (Gtm->lnso[i]==0)
          continue;
        // checking required scalar type and selection of the node and selection of the component
        if (Outdm->nog.selnstre.presence_id(Outdm->nog.selstre, i, dir, ir) == 0)
          continue;
        ncompstr = Mt->nodes[i].ncompstr;
        if (Outdm->nog.transtre[ir] > 0)
        {
          Mt->give_nodal_coord(i, p);
          Outdm->lcs[Outdm->nog.transtre[ir]-1].give_transfmat(t, p, Mp->time);
          str.n = ncompstr;
          str.a = Mt->nodes[i].stress+ncompstr*lcid;
          aux = 0.0;
          if (ncompstr == 4)
            gl_comp_engvectortransf(str, aux, dir, t, planestrain, stress);
          if (ncompstr == 6)
            gl_comp_engvectortransf(str, aux, dir, t, spacestress, stress);
          str.a = NULL;
          fprintf(out, "%ld % e\n", i+Outdm->idn1, aux);
        }
        else
          fprintf(out, "%ld % e\n", i+Outdm->idn1, Mt->nodes[i].stress[ncompstr*lcid+dir]);
      }
      break;
    case strain:
      for (i = 0; i < Mt->nn; i++)
      {
        if (Gtm->lnso[i]==0)
          continue;
        // checking required scalar type and selection of the node and selection of the component
        if (Outdm->nog.selnstra.presence_id(Outdm->nog.selstra, i, dir, ir) == 0)
          continue;
        ncompstr = Mt->nodes[i].ncompstr;
        if (Outdm->nog.transtra[ir] > 0)
        {
          Mt->give_nodal_coord(i, p);
          Outdm->lcs[Outdm->nog.transtra[ir]-1].give_transfmat(t, p, Mp->time);
          str.n = ncompstr;
          str.a = Mt->nodes[i].strain+ncompstr*lcid;
          aux = 0.0;
          if (ncompstr == 4)
            gl_comp_engvectortransf(str, aux, dir, t, planestrain, strain);
          if (ncompstr == 6)
            gl_comp_engvectortransf(str, aux, dir, t, spacestress, strain);
          str.a = NULL;
          fprintf(out, "%ld % e\n", i+Outdm->idn1, aux);
        }
        fprintf(out, "%ld % e\n", i+Outdm->idn1, Mt->nodes[i].strain[ncompstr*lcid+dir]);
      }
      break;
    case pstress:
      for (i = 0; i < Mt->nn; i++)
      {
        if (Gtm->lnso[i]==0)
          continue;
        if (Mt->nodes[i].pstre == NULL)
          continue;
        // checking required scalar type and selection of the node and selection of the component
        if (Outdm->nog.selnstre.presence_id(Outdm->nog.selstre, i, dir) == 0)
          continue;
        if (dir == 3)
          fprintf(out, "%ld % e\n", i+Outdm->idn1, (Mt->nodes[i].pstre[3*lcid]+Mt->nodes[i].pstre[3*lcid+2])/2);
        else
          fprintf(out, "%ld % e\n", i+Outdm->idn1, Mt->nodes[i].pstre[3*lcid+dir]);
      }
      break;
    case pstrain:
      for (i = 0; i < Mt->nn; i++)
      {
        if (Gtm->lnso[i]==0)
          continue;
        if (Mt->nodes[i].pstra == NULL)
          continue;
        // checking required scalar type and selection of the node and selection of the component
        if (Outdm->nog.selnstra.presence_id(Outdm->nog.selstra, i, dir) == 0)
          continue;
        fprintf(out, "%ld % e\n", i+Outdm->idn1, Mt->nodes[i].pstra[3*lcid+dir]);
      }
      break;
    case other:
      for (i = 0; i < Mt->nn; i++)
      {
        if (Gtm->lnso[i]==0)
          continue;
        // checking required scalar type and selection of the node and selection of the component
        if (Outdm->nog.selnoth.presence_id(Outdm->nog.seloth, i, dir) == 0)
          continue;
        ncompother = Mt->nodes[i].ncompother;
        fprintf(out, "%ld % e\n", i+Outdm->idn1, Mt->nodes[i].other[ncompother*lcid+dir]);
      }
      break;
    default:
      print_err("unknown value type is required", __FILE__, __LINE__, __func__);
  }
  fprintf(out, "End Values\n");
}



/**
  The function writes a vector value given by parameter scal for nodes slected by selection with index sid
  to the file given by parameter out in GiD format.

  @param out - pointer to the opened text file where the output will be produced
  @param q - specifies type of required quantity (strain/stres/other)
  @param lcid - load case id
  @param sid - index of nodal selection which will be printed
  @param desclcid - string with description of loadcase

  @return The function does not return anything.

  created 3.10.2007 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void write_gid_nodvector(FILE *out, strastre q, long lcid, long sid, const char *desclcid)
{
  char sig[70];
  long i, j, ncompstr, ncompother;
  long asid, q_id1, q_n;

  switch (q)
  {
    case strain:
      q_id1 = Outdm->nog.selstra[sid].id1[0];
      q_n = Outdm->nog.selstra[sid].ncomp[0];
      sprintf(sig, "eps_n_v%ld-%ld_s%ld", q_id1+1, q_n, sid+1);
      break;
    case stress:
      q_id1 = Outdm->nog.selstre[sid].id1[0];
      q_n = Outdm->nog.selstre[sid].ncomp[0];
      sprintf(sig, "sig_n_v%ld-%ld_s%ld", q_id1+1, q_n, sid+1);
      break;
    case other:
      q_id1 = Outdm->nog.seloth[sid].id1[0];
      q_n = Outdm->nog.seloth[sid].ncomp[0];
      sprintf(sig, "other_n_v%ld-%ld_s%ld", q_id1+1, q_n, sid+1);
      break;
    default:
      print_err("unsupported type of quantity (%d) is required", __FILE__, __LINE__, __func__, q);
      abort();
  }
  fprintf(out, "\nResult \"%s\" \"%ld\" %s Vector OnNodes\n", sig, lcid, desclcid);
  fprintf(out, "ComponentNames \"%s_1\" \"%s_2\" \"%s_3\"\nValues\n", sig, sig, sig);
  switch (q)
  {
    case strain:
      for (i = 0; i < Mt->nn; i++)
      {
        if (Gtm->lnso[i]==0)
          continue;
        // checking required scalar type and selection of the node and selection of the component
        if (Outdm->nog.selnstra.presence_id(i, sid, asid) == 0)
          continue;
        ncompstr = Mt->nodes[i].ncompstr;
        fprintf(out, "%ld", i+Outdm->idn1);
        if (q_id1+q_n > ncompstr)
        {
          print_err("index of selected vector component exceeds total number of strain components\n"
                    " (node=%ld, ncompstr=%ld, id1=%ld, n=%ld)\n", __FILE__, __LINE__, __func__, i+1, ncompstr, q_id1+1, q_n);
          abort();
        }
        for (j = q_id1; j < q_id1+q_n; j++)
          fprintf(out, " % e",  Mt->nodes[i].strain[ncompstr*lcid+j]);
        for (j = 0; j < 3-q_n; j++)
          fprintf(out, " 0.0");
        fprintf(out, "\n");
      }
      break;
    case stress:
      for (i = 0; i < Mt->nn; i++)
      {
        if (Gtm->lnso[i]==0)
          continue;
        // checking required scalar type and selection of the node and selection of the component
        if (Outdm->nog.selnstre.presence_id(i, sid, asid) == 0)
          continue;
        ncompstr = Mt->nodes[i].ncompstr;
        fprintf(out, "%ld", i+Outdm->idn1);
        if (q_id1+q_n > ncompstr)
        {
          print_err("index of selected vector component exceeds total number of stress component\n"
                    " (node=%ld, ncompstr=%ld, id1=%ld, n=%ld)\n", __FILE__, __LINE__, __func__, i+1, ncompstr, q_id1+1, q_n);
          abort();
        }
        for (j = q_id1; j < q_id1+q_n; j++)
          fprintf(out, " % e",  Mt->nodes[i].stress[ncompstr*lcid+j]);
        for (j = 0; j < 3-q_n; j++)
          fprintf(out, " 0.0");
        fprintf(out, "\n");
      }
      break;
    case other:
      for (i = 0; i < Mt->nn; i++)
      {
        if (Gtm->lnso[i]==0)
          continue;
        // checking required scalar type and selection of the node and selection of the component
        if (Outdm->nog.selnoth.presence_id(i, sid, asid) == 0)
          continue;
        ncompother = Mt->nodes[i].ncompother;
        fprintf(out, "%ld", i+Outdm->idn1);
        if (q_id1+q_n > ncompother)
        {
          print_err("index of selected vector component exceeds total number of other values component\n"
                    " (node=%ld, ncompother=%ld, id1=%ld, n=%ld)\n", __FILE__, __LINE__, __func__, i+1, ncompother, q_id1+1, q_n);
          abort();
        }
        for (j = q_id1; j < q_id1+q_n; j++)
          fprintf(out, " % e",  Mt->nodes[i].other[ncompother*lcid+j]);
        for (j = 0; j < 3-q_n; j++)
          fprintf(out, " 0.0");
        fprintf(out, "\n");
      }
      break;
    default:
      print_err("unknown quantity type is required", __FILE__, __LINE__, __func__);
      abort();
  }
  fprintf(out, "End Values\n");
}



/**
  The function writes a tensor(matrix) quantity given by parameter q on elements of given type (parameter te)
  to the file given by parameter out in GiD format. The results are printed at each  element integration 
  points.

  @param out - pointer to the opened text file where the output will be produced
  @param q - specifies type of required scalar quantity (strain/stres/other)
  @param lcid - load case id
  @param sid - index of selection of elements which will be printed
  @param desclcid - string with description of loadcase
  @param te - required element type

  @return The function does not return anything.
 
  created 4.10.2007 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void write_gid_nodtensor(FILE *out, strastre q, long lcid, long sid, const char *desclcid)
{
  char sig[70];
  long i, ncompstr, ncompother;
  long asid, q_id1, q_n;
  matrix t(ASTCKMAT(3,3)); // tensor represenation of the given quantity
  matrix tmat;   // transformation matrix to the local coordinate system
  vector str, p;
  vector v;
  strastrestate ssst;
  
  if(Outdm->nlcs > 0)
  {
    reallocv(RSTCKVEC(3, p));
    reallocm(RSTCKMAT(3, 3, tmat));
  }

  switch (q)
  {
    case strain:
      q_id1 = Outdm->nog.selstra[sid].id1[0];
      q_n = 0;
      sprintf(sig, "eps_n_m_s%ld", sid+1);
      break;
    case stress:
      q_id1 = Outdm->nog.selstre[sid].id1[0];
      q_n = 0;
      sprintf(sig, "sig_n_m_s%ld", sid+1);
      break;
    case other:
      q_id1 = Outdm->nog.seloth[sid].id1[0];
      q_n = Outdm->nog.seloth[sid].ncomp[0];
      sprintf(sig, "other_n_m_s%ld", sid+1);
      break;
    default:
      print_err("unsupported type of quantity (%d) is required", __FILE__, __LINE__, __func__, q);
      abort();
  }
  // detection of maximal number of stress/strain components in the selected elements
  if (q != other)
  {
    for (i = 0; i < Mt->nn; i++)
    {
      if (Gtm->lnso[i]==0)
        continue;
      switch (q)
      {
        case strain:
          // checking required scalar type and selection of the node and selection of the component
          if (Outdm->nog.selnstra.presence_id(i, sid, asid) == 0)
            continue;
          if (q_n < Mt->nodes[i].ncompstr)
            q_n = Mt->nodes[i].ncompstr;
          break;
        case stress:
          // checking required scalar type and selection of the node and selection of the component
          if (Outdm->nog.selnstre.presence_id(i, sid, asid) == 0)
            continue;
          if (q_n < Mt->nodes[i].ncompstr)
            q_n = Mt->nodes[i].ncompstr;
          break;
        default:
          break;
      } 
    }
    if (q_n == 0)
      return;
  }
  if ((q_n != 6) && (q_n != 4) && (q_n != 1))
  {
    print_err("wrong number of required tensor components (%ld)\n"
              " only 1 or 4 or 6 components are allowed.", __FILE__, __LINE__, __func__, q_n);
    abort();
  }
  if (q_n == 6)
    fprintf(out, "\nResult \"%s\" \"%ld\" %s  Matrix OnNodes \n", sig, lcid, desclcid);
  else 
    fprintf(out, "\nResult \"%s\" \"%ld\" %s  PlainDeformationMatrix OnNodes \n", sig, lcid, desclcid);
  fprintf(out, "Values\n");
  
  switch (q)
  {
    case strain:
      for (i = 0; i < Mt->nn; i++)
      {
        if (Gtm->lnso[i]==0)
          continue;
        // checking required scalar type and selection of the node and selection of the component
        if (Outdm->nog.selnstra.presence_id(i, sid, asid) == 0)
          continue;
        ncompstr = Mt->nodes[i].ncompstr;
        fprintf(out, "%ld", i+Outdm->idn1);
        reallocv(RSTCKVEC(ncompstr, v));
        copyv(Mt->nodes[i].strain+ncompstr*lcid+q_id1, v);
        ssst = guess_ssst(ncompstr);
        vector_tensor (v, t, ssst, strain);
        if (Outdm->nog.transtra[sid] > 0)
        {
          Mt->give_nodal_coord(i, p);
          Outdm->lcs[Outdm->nog.transtra[sid]-1].give_transfmat(tmat, p, Mp->time);
          glmatrixtransf(t, tmat);
        }
        if (q_n == 6)
          fprintf(out, " % e % e % e % e % e % e\n", t[0][0], t[1][1], t[2][2], t[0][1], t[1][2], t[0][2]);
        else
          fprintf(out, " % e % e % e % e\n", t[0][0], t[1][1], t[0][1], t[2][2]);
      }
      break;
    case stress:
      for (i = 0; i < Mt->nn; i++)
      {
        if (Gtm->lnso[i]==0)
          continue;
        // checking required scalar type and selection of the node and selection of the component
        if (Outdm->nog.selnstre.presence_id(i, sid, asid) == 0)
          continue;
        ncompstr = Mt->nodes[i].ncompstr;
        fprintf(out, "%ld", i+Outdm->idn1);
        reallocv(RSTCKVEC(ncompstr, v));
        copyv(Mt->nodes[i].stress+ncompstr*lcid+q_id1, v);
        ssst = guess_ssst(ncompstr);
        vector_tensor (v, t, ssst, strain);
        if (Outdm->nog.transtre[sid] > 0)
        {
          Mt->give_nodal_coord(i, p);
          Outdm->lcs[Outdm->nog.transtre[sid]-1].give_transfmat(tmat, p, Mp->time);
          glmatrixtransf(t, tmat);
        }
        if (q_n == 6)
          fprintf(out, " % e % e % e % e % e % e\n", t[0][0], t[1][1], t[2][2], t[0][1], t[1][2], t[0][2]);
        else
          fprintf(out, " % e % e % e % e\n", t[0][0], t[1][1], t[0][1], t[2][2]);
      }
      break;
    case other:
      for (i = 0; i < Mt->nn; i++)
      {
        if (Gtm->lnso[i]==0)
          continue;
        // checking required scalar type and selection of the node and selection of the component
        if (Outdm->nog.selnoth.presence_id(i, sid, asid) == 0)
          continue;
        ncompother = Mt->nodes[i].ncompother;
        fprintf(out, "%ld", i+Outdm->idn1);
        reallocv(RSTCKVEC(q_n, v));
	copyv(Mt->nodes[i].other+ncompother*lcid+q_id1, v);
        ssst = guess_ssst(q_n);
        vector_tensor (v, t, ssst, strain);
        if (q_n == 6)
          fprintf(out, " % e % e % e % e % e % e\n", t[0][0], t[1][1], t[2][2], t[0][1], t[1][2], t[0][2]);
        else
          fprintf(out, " % e % e % e % e\n", t[0][0], t[1][1], t[0][1], t[2][2]);
      }
      break;
    default:
      print_err("unknown quantity type is required", __FILE__, __LINE__, __func__);
      abort();
  }
  fprintf(out, "End Values\n");
}



/**
  The function writes a scalar quantity given by parameters scal and dir on all elements to the file given by parameter out in GiD format. 
  The results are printed at each integration point on given element.

  @param out - pointer to the opened text file where the output will be produced
  @param scal - specifies type of required scalar quantity (strain/stres/other)
  @param lcid - load case id
  @param dir - specifies which component of the quantity array will be printed
  @param desclcid - string with description of loadcase

  @return The function does not return anything.

  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void write_gid_elemscalar(FILE *out, strastre scal, long lcid, long dir, const char *desclcid)
{
  if (Bar2d)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, bar2d);

  if (Bar3d)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, bar3d);

  if (Barq2d)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, barq2d);

  if (Barq3d)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, barq3d);

  if (Beam2d)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, beam2d);

  if (Beam3d)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, beam3d);

  if (Sbeam)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, subsoilbeam);

  if (Spring)
  {
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, spring_1);
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, spring_2);
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, spring_3);
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, spring_4);
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, spring_5);
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, spring_6);
  }

  if (Pqifc)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, planequadinterface);

  if (Pelt)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, planeelementlt);

  if (Peqt)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, planeelementqt);

  if (Asymqt)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, axisymmqt);

  if (Perlt)
  {
    if (Shtr)
    {}
    else
      write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, planeelementrotlt);
  }

  if (Pelq)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, planeelementlq);

  if (Peqq)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, planeelementqq);

  if (Perlq)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, planeelementrotlq);

  if (Pesqt)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, planeelementsubqt);

  if (Cct)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, cctel);

  if (Dkt)
  {
    if (Shtr)
    {}
    else
      write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, dktel);
  }

  if (Dst)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, dstel);

  if (Q4pl)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, q4plateel);

  if (Spltr)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, subsoilplatetr);

  if (Splq)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, subsoilplateq);

  if (Asymlq)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, axisymmlq);

  if (Asymlt)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, axisymmlt);

  if (Asymqq)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, axisymmqq);

  if (Asymcq)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, axisymmcq);

  if (Asymlqifc)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, axisymmlqintface);

  if (Shtr)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, shelltrelem);

  if (Shtrm)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, shelltrmelem);

  if (Shq)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, shellqelem);

  if (Ltet)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, lineartet);

  if (Qtet)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, quadrtet);

  if (Lhex)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, linearhex);

  if (Qhex)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, quadrhex);

  if (Lwed)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, linearwed);

  if (Qwed)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, quadrwed);

  if (Hexifc)
    write_gid_elem_type_scalar(out, scal, lcid, dir, desclcid, hexintface);
}



/**
  The function writes a scalar value given by parameter scal on elements of given type (parameter te)
  to the file given by parameter out in GiD format. The results are printed at each  element integration 
  points.

  @param out - pointer to the opened text file where the output will be produced
  @param scal - specifies type of required scalar quantity (strain/stres/other)
  @param lcid - load case id
  @param dir - specifies which component of the quantity array will be printed
  @param desclcid - string with description of loadcase
  @param te - required element type

  @return The function does not return anything.

  Created 2004 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
  modified 04.2016 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void write_gid_elem_type_scalar(FILE *out, strastre scal, long lcid, long dir, const char *desclcid, elemtype te)
{
  const char *sig = "";
  char gpname[1000];
  long i, j, ipp,ncompstr,ncompother,tnip, ir;
  long print_header = 1;
  long num_copy_ip;
  matrix t;
  vector str, p;
  double aux;
  
  if(Outdm->nlcs > 0)
  {
    reallocv(RSTCKVEC(3, p));
    reallocm(RSTCKMAT(3, 3, t));
  }



  for (i = 0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    tnip = Mt->give_totnip(i);
    if (te == planeelementlq)
      tnip = 4;
    if (te == axisymmlq)
      tnip = 4;
    if (te == shelltrelem)
      tnip = 3;
    if (te == shelltrmelem)
      tnip = 1;
    if (te == cctel)
      tnip = Cct->nip[0][0];
    if (te == planeelementqt)
      tnip = Peqt->nip[0][0];
    if (te == axisymmqt)
      tnip = Asymqt->nip[0][0];
    if (te == planeelementrotlt)
      tnip = Perlt->nip[0][0];
    if (te == planeelementrotlq)
      tnip = Perlq->nip[0][0];
    if (te == beam2d)
      tnip = 2;
    num_copy_ip = 0;
    switch (te)
    {
      case beam3d:
      case subsoilbeam:
	num_copy_ip = 1;
        break;
      default:
        break;
    }

    // checking required element type
    if (Mt->elements[i].te != te)
      continue; 

    // checking required scalar type and selection of the element and selection of the component
    if ((scal == strain) && (Outdm->eog.selestra.presence_id(Outdm->eog.selstra, i, dir, ir) == 0))
      continue;
    if ((scal == stress) && (Outdm->eog.selestre.presence_id(Outdm->eog.selstre, i, dir, ir) == 0))
      continue;
    if ((scal == other) && (Outdm->eog.seleoth.presence_id(Outdm->eog.seloth, i, dir, ir) == 0))
      continue;

    // checking required scalar type and selection of the element and selection of the component
    // in case of sel_all option used for quantity component
    
    ipp = Mt->elements[i].ipp[0][0];
    if (te == shelltrelem)
    {
      ipp = Mt->elements[i].ipp[3][3];
      Mm->ip[ipp].ncompstr=6;
    }
    if (te == shelltrmelem)
    {
      ipp = Mt->elements[i].ipp[2][2];
    }
    if ((scal == strain) && (Outdm->eog.selstra[ir].st == sel_all) && (dir >= Mm->ip[ipp].ncompstr))
      continue;
    if ((scal == stress) && (Outdm->eog.selstre[ir].st == sel_all) && (dir >= Mm->ip[ipp].ncompstr))
      continue;
    if ((scal == other) && (Outdm->eog.seloth[ir].st == sel_all) && (dir >= Mm->ip[ipp].ncompeqother))
      continue;

    if (print_header)
    {
      if (scal == stress)
      {
        switch (Mm->ip[ipp].ssst)
        {
          case bar:
            switch (dir)
            {
              case 0:
                sig = "sig_e_1";
                break;
              case 1:
                sig = "sig_e_2";
                break;
              case 2:
                sig = "sig_e_3";
                break;
              case 3:
                sig = "sig_e_4";
                break;
              case 4:
                sig = "sig_e_5";
                break;
              case 5:
                sig = "sig_e_6";
                break;
              default:
                print_err("unknown direction of stresses is required", __FILE__, __LINE__, __func__);
            }
            break;
          case plbeam:
            switch (dir)
            {
              case 0:
                sig = "N";
                break;
              case 1:
                sig = "V";
                break;
              case 2:
                sig = "M";
                break;
              default:
                print_err("unknown direction of stresses is required", __FILE__, __LINE__, __func__);
            }
            break;
          case planestress:
          case planestrain: 
            switch (dir)
            {
              case 0:
                sig = "sig_e_x";
                break;
              case 1:
                sig = "sig_e_y";
                break;
              case 2:
                sig = "tau_e_xy";
                break;
              case 3:
                sig = "sig_e_z";
                break;
              default:
                print_err("unknown direction of stresses is required", __FILE__, __LINE__, __func__);
            }
            break;
          case planecontact:
            switch(dir){
              case 0:
                sig = "tau_1";
                break;
              case 1:
                if ((te == hexintface) || (te == axisymmlqintface))
                  sig = "tau_2";
                else
                  sig = "sig_n";
                break;
              case 2:
                sig = "sig_n";
                break;
              default:
                print_err("unknown direction of stresses is required", __FILE__, __LINE__, __func__);
            }
            break;
	case platek:{
	  switch (dir){
	  case 0:
	    sig = "m_x";
	    break;
	  case 1:
	    sig = "m_y";
	    break;
	  case 2:
	    sig = "m_xy";
	    break;
	  default:
	    print_err("unknown direction of stress is required", __FILE__, __LINE__, __func__);
	  }
	  break;
	}
	  
	case plates:{
	  switch (dir){
	  case 0:
	    sig = "m_x";
	    break;
	  case 1:
	    sig = "m_y";
	    break;
	  case 2:
	    sig = "m_xy";
	    break;
	  case 3:
	    sig = "q_x";
	    break;
	  case 4:
	    sig = "q_y";
	    break;
	  default:
	    print_err("unknown direction of stress is required", __FILE__, __LINE__, __func__);
	  }
	  break;
	}
	  
	case axisymm:
            switch (dir)
            {
              case 0:
                sig = "sig_e_r";
                break;
              case 1:
                sig = "sig_e_y";
                break;
              case 2:
                sig = "sig_e_phi";
                break;
              case 3:
                sig = "sig_e_ry";
                break;
              default:
                print_err("unknown direction of stresses is required", __FILE__, __LINE__, __func__);
	    }
            break;

	case shell:{
	  switch (dir){
	  case 0:
	    sig = "n_x";
	    break;
	  case 1:
	    sig = "n_y";
	    break;
	  case 2:
	    sig = "n_xy";
	    break;
	  case 3:
	    sig = "m_x";
	    break;
	  case 4:
	    sig = "m_y";
	    break;
	  case 5:
	    sig = "m_xy";
	    break;
	  default:
	    print_err("unknown direction of strains is required", __FILE__, __LINE__, __func__);
	  }
	  break;
	}

          case spacestress:
            switch (dir)
            {
              case 0:
                sig = "sig_e_x";
                break;
              case 1:
                sig = "sig_e_y";
                break;
              case 2:
                sig = "sig_e_z";
                break;
              case 3:
                sig = "tau_e_yz";
                break;
              case 4:
                sig = "tau_e_xz";
                break;
              case 5:
                sig = "tau_e_xy";
                break;
              default:
                print_err("unknown direction of stresses is required", __FILE__, __LINE__, __func__);
            }
            break;
          default:
            print_err("unknown stress/strain state is required", __FILE__, __LINE__, __func__);
        }
      }
      if (scal == strain)
      {
        switch (Mm->ip[ipp].ssst)
        {
          case bar:
            switch (dir)
            {
              case 0:
                sig = "eps_e_1";
                break;
              case 1:
                sig = "eps_e_2";
                break;
              case 2:
                sig = "eps_e_3";
                break;
              case 3:
                sig = "eps_e_4";
                break;
              case 4:
                sig = "eps_e_5";
                break;
              case 5:
                sig = "eps_e_6";
                break;
              default:
                print_err("unknown direction of strains is required", __FILE__, __LINE__, __func__);
            }
	    break;
          case plbeam:
            switch (dir)
            {
              case 0:
                sig = "u";
                break;
              case 1:
                sig = "w";
                break;
              case 2:
                sig = "phi";
                break;
              default:
                print_err("unknown direction of strains is required", __FILE__, __LINE__, __func__);
            }
            break;
          case planestress:
          case planestrain: 
            switch (dir)
            {
              case 0:
                sig = "eps_e_x";
                break;
              case 1:
                sig = "eps_e_y";
                break;
              case 2:
                sig = "eps_e_xy";
                break;
              case 3:
                sig = "eps_e_z";
                break;
              default:
                print_err("unknown direction of strains is required", __FILE__, __LINE__, __func__);
	    }
            break;
	    
          case planecontact:
            switch(dir){
              case 0:
                sig = "slip_1";
                break;
              case 1:
                if ((te == hexintface) || (te == axisymmlqintface))
                  sig = "slip_2";
                else
                  sig = "norm_d";
                break;
              case 2:
                sig = "norm_d";
                break;
              default:
                print_err("unknown direction of strains is required", __FILE__, __LINE__, __func__);
            }
            break;
	case plates:{
	  switch (dir){
	  case 0:{
	    sig = "kappa_x";
	    break;
	  }
	  case 1:{
	    sig = "kappa_y";
	    break;
	  }
	  case 2:{
	    sig = "kappa_xy";
	    break;
	  }
	  case 3:{
	    sig = "gamma_x";
	    break;
	  }
	  case 4:{
	    sig = "gamma_y";
	    break;
	  }
	    
	  default:{
	    print_err("unknown direction of strain is required", __FILE__, __LINE__, __func__);
	  }
	  }
	  break;
	}
	  
	case platek:{
	  switch (dir){
	  case 0:{
	    sig = "kappa_x";
	    break;
	  }
	  case 1:{
	    sig = "kappa_y";
	    break;
	  }
	  case 2:{
	    sig = "kappa_xy";
	    break;
	  }
	  default:{
	    print_err("unknown direction of strain is required", __FILE__, __LINE__, __func__);
	  }
	  }
	  break;
	}
	  
          case axisymm:
            switch (dir)
            {
              case 0:
                sig = "eps_e_r";
                break;
              case 1:
                sig = "eps_e_y";
                break;
              case 2:
                sig = "eps_e_phi";
                break;
              case 3:
                sig = "eps_e_ry";
                break;
              default:
                print_err("unknown direction of strains is required", __FILE__, __LINE__, __func__);
	    }
            break;

	case shell:{
	  switch (dir){
	  case 0:
	    sig = "eps_x";
	    break;
	  case 1:
	    sig = "eps_y";
	    break;
	  case 2:
	    sig = "gamma_xy";
	    break;
	  case 3:
	    sig = "kappa_x";
	    break;
	  case 4:
	    sig = "kappa_y";
	    break;
	  case 5:
	    sig = "kappa_xy";
	    break;
	  default:
	    print_err("unknown direction of strains is required", __FILE__, __LINE__, __func__);
	  }
	  break;
	}
	  
          case spacestress:
            switch (dir)
            {
              case 0:
                sig = "eps_e_x";
                break;
              case 1:
                sig = "eps_e_y";
                break;
              case 2:
                sig = "eps_e_z";
                break;
              case 3:
                sig = "eps_e_yz";
                break;
              case 4:
                sig = "eps_e_xz";
                break;
              case 5:
                sig = "eps_e_xy";
                break;
              default:
                print_err("unknown direction of strains is required", __FILE__, __LINE__, __func__);
            }
            break;
          default:
            print_err("unknown stress/strain state is required", __FILE__, __LINE__, __func__);
        }
      }
      switch (te)
      {
        case bar2d:
        case bar3d:
          sprintf(gpname, "Lin_1D");
          break;
        case beam2d:
        case beam3d:
        case subsoilbeam:
          sprintf(gpname, "Beam_1D");
          break;
        case barq2d:
        case barq3d:
          sprintf(gpname, "Quad_1D");
          break;
        default:
          sprintf(gpname, "%d", te);
      }
      if (scal == other)
        fprintf(out, "\nResult \"other_e_%ld\" \"%ld\" %s Scalar OnGaussPoints \"%s\"\n", dir+1, lcid, desclcid, gpname);
      else     
        fprintf(out, "\nResult \"%s\" \"%ld\" %s Scalar OnGaussPoints \"%s\"\n", sig, lcid, desclcid, gpname);
      fprintf(out, "Values\n");
      print_header = 0;
    }
    switch (scal)
    {
      case strain:
        fprintf(out, "%7ld", i+Outdm->ide1);
        for (j = 0; j < tnip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0]+j;
          if (te == axisymmlq && Mt->give_totnip(i)==16)
          {
            if(j==1)   ipp=Mt->elements[i].ipp[0][0]+3;
            if(j==2)   ipp=Mt->elements[i].ipp[0][0]+12;
            if(j==3)   ipp=Mt->elements[i].ipp[0][0]+15;
          }
          ncompstr = Mm->ip[ipp].ncompstr;
          if (te == shelltrelem)
          {
	    ipp = Mt->elements[i].ipp[3][3]+j;
            ncompstr = 6;          
            reallocv(ncompstr, str);
            rcopyv(Mm->ip[ipp].strain, 0, ncompstr, str.a, 0, ncompstr, ncompstr);	    
          }
          else
            makerefv(str, Mm->ip[ipp].strain+ncompstr*lcid, ncompstr);

          if (j == 0)
          {
            // on soilplatetr skip the first integration point in the middle of element
            if (te == subsoilplatetr)
              continue;

            if (Outdm->eog.transtra[ir] > 0)
            {
              ipcoord(Mm->elip[ipp], ipp, 0, 0, p);
              Outdm->lcs[Outdm->eog.transtra[ir]-1].give_transfmat(t, p, Mp->time);
              gl_comp_engvectortransf(str, aux, dir, t, Mm->ip[ipp].ssst, strain);
              fprintf(out, " % e\n", aux);
            }
            else
              fprintf(out, " % e\n", str[dir]);
          }
          else
          {
            if (Outdm->eog.transtra[ir] > 0)
            {
              ipcoord(Mm->elip[ipp], ipp, 0, 0, p);
              Outdm->lcs[Outdm->eog.transtra[ir]-1].give_transfmat(t, p, Mp->time);
              gl_comp_engvectortransf(str, aux, dir, t, Mm->ip[ipp].ssst, strain);
              fprintf(out, "%7c % e\n", ' ', aux);
            }
            else
              fprintf(out, "%7c % e\n", ' ', str[dir]);
          }
        }
        ipp = Mt->elements[i].ipp[0][0];
        ncompstr = Mm->ip[ipp].ncompstr;
        if (Outdm->eog.transtra[ir] > 0)
        {
          ipcoord(Mm->elip[ipp], ipp, 0, 0, p);
          Outdm->lcs[Outdm->eog.transtra[ir]-1].give_transfmat(t, p, Mp->time);
          gl_comp_engvectortransf(str, aux, dir, t, Mm->ip[ipp].ssst, strain);
        }
        for (j = 0; j < num_copy_ip; j++)
        {
          if (Outdm->eog.transtra[ir] > 0)
            fprintf(out, "%7c % e\n", ' ', aux);
          else
            fprintf(out, "%7c % e\n", ' ', Mm->ip[ipp].strain[ncompstr*lcid+dir]);
        }
        fprintf(out, "\n");
        break;
      case stress:
        fprintf(out, "%7ld", i+Outdm->ide1);
        for (j = 0; j < tnip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0]+j;
          ncompstr = Mm->ip[ipp].ncompstr;
          if (te == axisymmlq && Mt->give_totnip(i)==16)
          {
            if(j==1)   ipp=Mt->elements[i].ipp[0][0]+3;
            if(j==2)   ipp=Mt->elements[i].ipp[0][0]+12;
            if(j==3)   ipp=Mt->elements[i].ipp[0][0]+15;
          }
          if (te == shelltrelem)
          {
	    ipp = Mt->elements[i].ipp[3][3]+j;
            ncompstr = 6;          
            reallocv(ncompstr, str);
            rcopyv(Mm->ip[ipp].stress, 0, ncompstr, str.a, 0, ncompstr, ncompstr);            
          }
          else
            makerefv(str, Mm->ip[ipp].stress+ncompstr*lcid, ncompstr);
          if (j == 0)
          {
            // on soilplatetr skip the first integration point in the middle of element
            if (te == subsoilplatetr)
              continue;

            if (Outdm->eog.transtre[ir] > 0)
            {
              ipcoord(Mm->elip[ipp], ipp, 0, 0, p);
              Outdm->lcs[Outdm->eog.transtre[ir]-1].give_transfmat(t, p, Mp->time);
              gl_comp_engvectortransf(str, aux, dir, t, Mm->ip[ipp].ssst, stress);
              fprintf(out, " % e\n", aux);
            }
            else
            {
              //              fprintf(out, " % e\n", Mm->ip[ipp].stress[ncompstr*lcid+dir]);
              fprintf(out, " % e\n", str[dir]);
            }
          }
          else
          {
            if (Outdm->eog.transtre[ir] > 0)
            {
              ipcoord(Mm->elip[ipp], ipp, 0, 0, p);
              Outdm->lcs[Outdm->eog.transtre[ir]-1].give_transfmat(t, p, Mp->time);
              gl_comp_engvectortransf(str, aux, dir, t, Mm->ip[ipp].ssst, stress);
              fprintf(out, "%7c % e\n", ' ', aux);
            }
            else
              //              fprintf(out, "%7c % e\n", ' ', Mm->ip[ipp].stress[ncompstr*lcid+dir]);
              fprintf(out, "%7c % e\n", ' ', str[dir]);
          }
        }
        ipp = Mt->elements[i].ipp[0][0];
        if (Outdm->eog.transtre[ir] > 0)
        {
          ipcoord(Mm->elip[ipp], ipp, 0, 0, p);
          Outdm->lcs[Outdm->eog.transtre[ir]-1].give_transfmat(t, p, Mp->time);
          gl_comp_engvectortransf(str, aux, dir, t, Mm->ip[ipp].ssst, stress);
        }
        for (j = 0; j < num_copy_ip; j++)
        {
          if (Outdm->eog.transtre[ir] > 0)
            fprintf(out, "%7c % e\n", ' ', aux);
          else
            //            fprintf(out, "%7c % e\n", ' ', Mm->ip[ipp].stress[ncompstr*lcid+dir]);
            fprintf(out, "%7c % e\n", ' ', str[dir]);
        }
        fprintf(out, "\n");
        break;
      case other:
        fprintf(out, "%7ld", i+Outdm->ide1);
        for (j = 0; j < tnip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0]+j;
	  
	  if (te == shelltrelem)
            ipp = Mt->elements[i].ipp[3][3]+j;
	  
          if (te == axisymmlq && Mt->give_totnip(i)==16)
          {
            if(j==1)   ipp=Mt->elements[i].ipp[0][0]+3;
            if(j==2)   ipp=Mt->elements[i].ipp[0][0]+12;
            if(j==3)   ipp=Mt->elements[i].ipp[0][0]+15;
          }
          ncompother = Mm->ip[ipp].ncompeqother;
          if (j == 0)
          {
            // on soilplatetr skip the first integration point in the middle of element
            if (te == subsoilplatetr)
              continue;

            fprintf(out, " % e\n", Mm->ip[ipp].other[ncompother*lcid+dir]);
          }
          else
            fprintf(out, "%7c % e\n", ' ', Mm->ip[ipp].other[ncompother*lcid+dir]);
        }
        for (j = 0; j < num_copy_ip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0];
          ncompother = Mm->ip[ipp].ncompeqother;
          fprintf(out, "%7c % e\n", ' ', Mm->ip[ipp].other[ncompother*lcid+dir]);
        }
/*        for (j = 0; j < tnip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0]+j;
          ncompother = Mm->ip[ipp].ncompother;
          if (j == 0)
            fprintf(out, " % e\n", Mm->ip[ipp].other[ncompother*lcid+dir]);
          else
            fprintf(out, "%7c % e\n", ' ', Mm->ip[ipp].other[ncompother*lcid+dir]);
        }
        for (j = 0; j < num_copy_ip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0];
          ncompother = Mm->ip[ipp].ncompother;
          fprintf(out, "%7c % e\n", ' ', Mm->ip[ipp].other[ncompother*lcid+dir]);
        }*/
        fprintf(out, "\n");
        break;
      default:
        print_err("unknown value type is required",__FILE__, __LINE__, __func__);
    }
  }
  if (print_header == 0)
    fprintf(out, "End Values\n");
}



/**
  The function writes a vector quantity given by parameters q and sid on all selected elements to the file 
  given by parameter out in GiD format. 
  The results are printed at each integration point on given element.

  @param out  - pointer to the opened text file where the output will be produced
  @param q    - specifies type of required vector quantity (strain/stres/other)
  @param lcid - load case id
  @param sid  - index of element selection whose selected quantity array will be printed
  @param desclcid - string with description of loadcase

  @return The function does not return anything.
 
  created 4.10.2007 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void write_gid_elemvector(FILE *out, strastre q, long lcid, long sid, const char *desclcid)
{
  if (Bar2d)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, bar2d);

  if (Bar3d)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, bar3d);

  if (Barq2d)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, barq2d);

  if (Barq3d)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, barq3d);

  if (Beam2d)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, beam2d);

  if (Beam3d)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, beam3d);

  if (Sbeam)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, subsoilbeam);

  if (Spring)
  {
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, spring_1);
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, spring_2);
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, spring_3);
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, spring_4);
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, spring_5);
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, spring_6);
  }

  if (Pqifc)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, planequadinterface);

  if (Pelt)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, planeelementlt);

  if (Peqt)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, planeelementqt);

  if (Asymqt)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, axisymmqt);

  if (Perlt)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, planeelementrotlt);

  if (Pelq)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, planeelementlq);

  if (Peqq)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, planeelementqq);

  if (Perlq)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, planeelementrotlq);

  if (Pesqt)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, planeelementsubqt);

  if (Cct)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, cctel);

  if (Dkt)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, dktel);

  if (Dst)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, dstel);

  if (Q4pl)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, q4plateel);

  if (Spltr)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, subsoilplatetr);

  if (Splq)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, subsoilplateq);

  if (Asymlq)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, axisymmlq);

  if (Asymlt)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, axisymmlt);

  if (Asymqq)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, axisymmqq);

  if (Asymcq)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, axisymmcq);

  if (Asymlqifc)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, axisymmlqintface);

  if (Shtr)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, shelltrelem);

  if (Shq)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, shellqelem);

  if (Ltet)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, lineartet);

  if (Qtet)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, quadrtet);

  if (Lhex)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, linearhex);

  if (Qhex)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, quadrhex);

  if (Lwed)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, linearwed);

  if (Qwed)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, quadrwed);

  if (Hexifc)
    write_gid_elem_type_vector(out, q, lcid, sid, desclcid, hexintface);

}



/**
  The function writes a vector quantity given by parameter q on elements of given type (parameter te)
  to the file given by parameter out in GiD format. The results are printed at each  element integration 
  points.

  @param out - pointer to the opened text file where the output will be produced
  @param q - specifies type of required vector quantity (strain/stres/other)
  @param lcid - load case id
  @param sid - index of selection of elements which will be printed
  @param desclcid - string with description of loadcase
  @param te - required element type

  @return The function does not return anything.

  created 4.10.2007 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void write_gid_elem_type_vector(FILE *out, strastre q, long lcid, long sid, const char *desclcid, elemtype te)
{
  char sig[70];
  char gpname[1000];
  long i, j, k, ipp,ncompstr,ncompother,tnip;
  long print_header = 1;
  long asid, q_id1, q_n;
  long num_copy_ip;


  switch (q)
  {
    case strain:
      q_id1 = Outdm->eog.selstra[sid].id1[0];
      q_n = Outdm->eog.selstra[sid].ncomp[0];
      sprintf(sig, "eps_e_v%ld-%ld_s%ld", q_id1+1, q_n, sid+1);
      break;
    case stress:
      q_id1 = Outdm->eog.selstre[sid].id1[0];
      q_n = Outdm->eog.selstre[sid].ncomp[0];
      sprintf(sig, "sig_e_v%ld-%ld_s%ld", q_id1+1, q_n, sid+1);
      break;
    case other:
      q_id1 = Outdm->eog.seloth[sid].id1[0];
      q_n = Outdm->eog.seloth[sid].ncomp[0];
      sprintf(sig, "other_e_v%ld-%ld_s%ld", q_id1+1, q_n, sid+1);
      break;
    default:
      print_err("unsupported type of quantity (%d) is required in function\n", __FILE__, __LINE__, __func__, q);
      abort();
  }

  for (i = 0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    // checking required element type
    if (Mt->elements[i].te != te)
      continue; 

    // checking required scalar type and selection of the element and selection id
    if (q == strain)
    {
      if (Outdm->eog.selestra.presence_id(i, sid, asid) == 0)
        continue;
    }
    if (q == stress)
    { 
      if (Outdm->eog.selestre.presence_id(i, sid, asid) == 0)
        continue;
    }
    if (q == other)
    { 
      if (Outdm->eog.seleoth.presence_id(i, sid, asid) == 0)
        continue;
    }

    tnip = Mt->give_totnip(i);
    if ((te == planeelementlq) || (te == axisymmlq))
      tnip = 4;
    if (te == shelltrelem)
      tnip = 3;
    if (te == shelltrmelem)
      tnip = 1;
    if (te == cctel)
      tnip = Cct->nip[0][0];
    if (te == planeelementqt)
      tnip = Peqt->nip[0][0];
    if (te == axisymmqt)
      tnip = Asymqt->nip[0][0];
    if (te == planeelementrotlq)
      tnip = Perlq->nip[0][0];
    num_copy_ip = 0;
    switch (te)
    {
      case beam3d:
      case subsoilbeam:
	num_copy_ip = 1;
        break;
      default:
        break;
    }

    if (print_header)
    {
      switch (te)
      {
        case bar2d:
        case bar3d:
        case beam2d:
        case beam3d:
        case subsoilbeam:
          sprintf(gpname, "Lin_1D");
          break;
        case barq2d:
        case barq3d:
          sprintf(gpname, "Quad_1D");
          break;
        default:
          sprintf(gpname, "%d", te);
      }
      fprintf(out, "\nResult \"%s\" \"%ld\" %s Vector OnGaussPoints \"%s\"\n", sig, lcid, desclcid, gpname);
      fprintf(out, "Values\n");
      print_header = 0;
    }
    switch (q)
    {
      case stress:
        for (j = 0; j < tnip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0]+j;
          if ((te == axisymmlq) && (tnip == 16)) 
          {
            if(j==1)   ipp=Mt->elements[i].ipp[0][0]+3;
            if(j==2)   ipp=Mt->elements[i].ipp[0][0]+12;
            if(j==3)   ipp=Mt->elements[i].ipp[0][0]+15;
          }
          if (j == 0)
            fprintf(out, "%7ld", i+Outdm->ide1);
          else
            fprintf(out, "%7c", ' ');
          for (k = q_id1; k < q_id1+q_n; k++)
            //            fprintf(out, " % e", Mm->ip[ipp].stress[ncompstr*lcid+k]);
            fprintf(out, " % e", Mm->ip[ipp].stress[k]);
          for (k = 0; k < 3-q_n; k++)
             fprintf(out, " 0.0");
          fprintf(out, "\n");
        }
        for (j = 0; j < num_copy_ip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0];
          fprintf(out, "%7c", ' ');
          for (k = q_id1; k < q_id1+q_n; k++)
            //            fprintf(out, " % e", Mm->ip[ipp].stress[ncompstr*lcid+k]);
            fprintf(out, " % e", Mm->ip[ipp].stress[k]);
          for (k = 0; k < 3-q_n; k++)
             fprintf(out, " 0.0");
          fprintf(out, "\n");
        }
        fprintf(out, "\n");
        break;
      case strain:
        for (j = 0; j < tnip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0]+j;
          if ((te == axisymmlq) && (tnip == 16))
          {
            if(j==1)   ipp=Mt->elements[i].ipp[0][0]+3;
            if(j==2)   ipp=Mt->elements[i].ipp[0][0]+12;
            if(j==3)   ipp=Mt->elements[i].ipp[0][0]+15;
          }
          ncompstr = Mm->ip[ipp].ncompstr;
          if (j == 0)
            fprintf(out, "%7ld", i+Outdm->ide1);
          else
            fprintf(out, "%7c", ' ');
          for (k = q_id1; k < q_id1+q_n; k++)
            fprintf(out, " % e", Mm->ip[ipp].strain[ncompstr*lcid+k]);
          for (k = 0; k < 3-q_n; k++)
             fprintf(out, " 0.0");
          fprintf(out, "\n");
        }
        for (j = 0; j < num_copy_ip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0];
          ncompstr = Mm->ip[ipp].ncompstr;
          fprintf(out, "%7c", ' ');
          for (k = q_id1; k < q_id1+q_n; k++)
            fprintf(out, " % e", Mm->ip[ipp].strain[ncompstr*lcid+k]);
          for (k = 0; k < 3-q_n; k++)
             fprintf(out, " 0.0");
          fprintf(out, "\n");
        }
        fprintf(out, "\n");
        break;
      case other:
        for (j = 0; j < tnip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0]+j;
          if ((te == axisymmlq) && (tnip == 16))
          {
            if(j==1)   ipp=Mt->elements[i].ipp[0][0]+3;
            if(j==2)   ipp=Mt->elements[i].ipp[0][0]+12;
            if(j==3)   ipp=Mt->elements[i].ipp[0][0]+15;
          }
          ncompother = Mm->ip[ipp].ncompeqother;
          if (j == 0)
            fprintf(out, "%7ld", i+Outdm->ide1);
          else
            fprintf(out, "%7c", ' ');
          for (k = q_id1; k < q_id1+q_n; k++)
            fprintf(out, " % e", Mm->ip[ipp].eqother[ncompother*lcid+k]);
          for (k = 0; k < 3-q_n; k++)
             fprintf(out, " 0.0");
          fprintf(out, "\n");
        }
        for (j = 0; j < num_copy_ip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0];
          ncompother = Mm->ip[ipp].ncompeqother;
          fprintf(out, "%7c", ' ');
          for (k = q_id1; k < q_id1+q_n; k++)
            fprintf(out, " % e", Mm->ip[ipp].eqother[ncompother*lcid+k]);
          for (k = 0; k < 3-q_n; k++)
             fprintf(out, " 0.0");
          fprintf(out, "\n");
        }
        fprintf(out, "\n");
        break;
      default:
        print_err("unknown value type is rquired", __FILE__, __LINE__, __func__);
    }
  }
  if (print_header == 0)
    fprintf(out, "End Values\n");
}



/**
  The function writes a tensor quantity given by parameters q and sid on all selected elements to the file 
  given by parameter out in GiD format. 
  The results are printed at each integration point on given element.

  @param out  - pointer to the opened text file where the output will be produced
  @param q    - specifies type of required tensor quantity (strain/stres/other)
  @param lcid - load case id
  @param sid  - index of element selection whose selected quantity array will be printed
  @param desclcid - string with description of loadcase

  @return The function does not return anything.
 
  created 4.10.2007 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void write_gid_elemtensor(FILE *out, strastre q, long lcid, long sid, const char *desclcid)
{
  if (Bar2d)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, bar2d);

  if (Bar3d)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, bar3d);

  if (Barq2d)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, barq2d);

  if (Barq3d)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, barq3d);

  if (Beam2d)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, beam2d);

  if (Beam3d)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, beam3d);

  if (Sbeam)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, subsoilbeam);

  if (Spring)
  {
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, spring_1);
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, spring_2);
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, spring_3);
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, spring_4);
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, spring_5);
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, spring_6);
  }

  if (Pqifc)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, planequadinterface);

  if (Pelt)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, planeelementlt);

  if (Peqt)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, planeelementqt);

  if (Asymqt)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, axisymmqt);

  if (Perlt)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, planeelementrotlt);

  if (Pelq)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, planeelementlq);

  if (Peqq)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, planeelementqq);

  if (Perlq)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, planeelementrotlq);

  if (Pesqt)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, planeelementsubqt);

  if (Cct)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, cctel);

  if (Dkt)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, dktel);

  if (Dst)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, dstel);

  if (Q4pl)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, q4plateel);

  if (Spltr)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, subsoilplatetr);

  if (Splq)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, subsoilplateq);

  if (Asymlq)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, axisymmlq);

  if (Asymlt)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, axisymmlt);

  if (Asymqq)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, axisymmqq);

  if (Asymcq)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, axisymmcq);

  if (Asymlqifc)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, axisymmlqintface);

  if (Shtr)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, shelltrelem);

  if (Shq)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, shellqelem);

  if (Ltet)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, lineartet);

  if (Qtet)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, quadrtet);

  if (Lhex)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, linearhex);

  if (Qhex)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, quadrhex);

  if (Lwed)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, linearwed);

  if (Qwed)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, quadrwed);

  if (Hexifc)
    write_gid_elem_type_tensor(out, q, lcid, sid, desclcid, hexintface);
}



/**
  The function writes a tensor(matrix) quantity given by parameter q on elements of given type (parameter te)
  to the file given by parameter out in GiD format. The results are printed at each  element integration 
  points.

  @param out - pointer to the opened text file where the output will be produced
  @param q - specifies type of required scalar quantity (strain/stres/other)
  @param lcid - load case id
  @param sid - index of selection of elements which will be printed
  @param desclcid - string with description of loadcase
  @param te - required element type

  @return The function does not return anything.

  created 4.10.2007 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
void write_gid_elem_type_tensor(FILE *out, strastre q, long lcid, long sid, const char *desclcid, elemtype te)
{
  char sig[70];
  char gpname[1000];
  long i, j, ipp, ncompstr, ncompother, tnip;
  long print_header = 1;
  long asid, q_id1, q_n;
  long num_copy_ip;
  vector v;
  matrix t(ASTCKMAT(3,3)); // tensor represenation of the given quantity
  matrix tmat;   // transformation matrix to the local coordinate system
  vector str, p;
  
  if(Outdm->nlcs > 0)
  {
    reallocv(RSTCKVEC(3, p));
    reallocm(RSTCKMAT(3, 3, tmat));
  }



  switch (q)
  {
    case strain:
      q_id1 = 0;
      q_n = 0;
      sprintf(sig, "eps_e_m_s%ld", sid+1);
      break;
    case stress:
      q_id1 = 0;
      q_n = 0;
      sprintf(sig, "sig_e_m_s%ld", sid+1);
      break;
    case other:
      q_id1 = Outdm->eog.seloth[sid].id1[0];
      q_n = Outdm->eog.seloth[sid].ncomp[0];
      sprintf(sig, "other_e_m%ld-%ld_s%ld", q_id1+1, q_n, sid+1);
      break;
    default:
      print_err("unsupported type of quantity (%d) is required", __FILE__, __LINE__, __func__, q);
      abort();
  }

  // detection of maximal number of stress/strain components in the selected elements
  if (q != other)
  {
    for (i = 0; i < Mt->ne; i++)
    {
      if ((Gtm->leso[i]==0) || (Mt->elements[i].te != te))
        continue;
      ipp = Mt->elements[i].ipp[0][0];
      switch (q)
      {
        case strain:
          if (Outdm->eog.selestra.presence_id(i, sid, asid) == 0)
            continue;
          if (q_n < Mm->ip[ipp].ncompstr)
            q_n = Mm->ip[ipp].ncompstr;
          break;
        case stress:
          if (Outdm->eog.selestre.presence_id(i, sid, asid) == 0)
            continue;
          if (q_n < Mm->ip[ipp].ncompstr)
            q_n = Mm->ip[ipp].ncompstr;
          break;
        default:
          break;
      } 
    }
    if (q_n == 0)
      return;
  }
  if ((q_n != 6) && (q_n != 4) && (q_n != 2) && (q_n != 1))
  {
    print_err("invalid number of required tensor components (%ld),\n"
              " only 1, 2, 4 or 6 components are allowed.", __FILE__, __LINE__, __func__, q_n);
    abort();
  }
 
  for (i = 0; i < Mt->ne; i++)
  {
    // checking whether the element is switched on
    if (Gtm->leso[i]==0)
      continue;

    // checking required element type
    if (Mt->elements[i].te != te)
      continue; 

    // checking required scalar type and selection of the element and selection id
    if ((q == strain) && (Outdm->eog.selestra.presence_id(i, sid, asid) == 0))
      continue;
    if ((q == stress) && (Outdm->eog.selestre.presence_id(i, sid, asid) == 0))
      continue;
    if ((q == other) && (Outdm->eog.seleoth.presence_id(i, sid, asid) == 0))
      continue;

    tnip = Mt->give_totnip(i);
    if ((te == planeelementlq) || (te == axisymmlq))
      tnip = 4;
    if (te == shelltrelem)
      tnip = 3;
    if (te == shelltrmelem)
      tnip = 1;
    if (te == cctel)
      tnip = Cct->nip[0][0];
    if (te == planeelementqt)
      tnip = Peqt->nip[0][0];
    if (te == axisymmqt)
      tnip = Asymqt->nip[0][0];
    if (te == planeelementrotlq)
      tnip = Perlq->nip[0][0];
    num_copy_ip = 0;
    switch (te)
    {
      case beam3d:
      case subsoilbeam:
	num_copy_ip = 1;
        break;
      default:
        break;
    }

    if (print_header)
    {
      switch (te)
      {
        case bar2d:
        case bar3d:
        case beam2d:
        case beam3d:
        case subsoilbeam:
          sprintf(gpname, "Lin_1D");
          break;
        case barq2d:
        case barq3d:
          sprintf(gpname, "Quad_1D");
          break;
        default:
          sprintf(gpname, "%d", te);
      }
      if (q_n == 6)
        fprintf(out, "\nResult \"%s\" \"%ld\" %s  Matrix OnGaussPoints \"%s\"\n", sig, lcid, desclcid, gpname);
      else 
        fprintf(out, "\nResult \"%s\" \"%ld\" %s  PlainDeformationMatrix OnGaussPoints \"%s\"\n", sig, lcid, desclcid, gpname);
      fprintf(out, "Values\n");
      print_header = 0;
    }
    switch (q)
    {
      case strain:
        for (j = 0; j < tnip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0]+j;
          if ((te == axisymmlq) && (tnip == 16))
          {
            if(j==1)   ipp=Mt->elements[i].ipp[0][0]+3;
            if(j==2)   ipp=Mt->elements[i].ipp[0][0]+12;
            if(j==3)   ipp=Mt->elements[i].ipp[0][0]+15;
          }
          ncompstr = Mm->ip[ipp].ncompstr;
          if (j == 0)
            fprintf(out, "%7ld", i+Outdm->ide1);
          else
            fprintf(out, "%7c", ' ');
          reallocv(RSTCKVEC(ncompstr, v));
          copyv(Mm->ip[ipp].strain+ncompstr*lcid+q_id1, v);
          vector_tensor (v, t, Mm->ip[ipp].ssst, strain);
          if (Outdm->eog.transtra[sid] > 0)
          {
            ipcoord(Mm->elip[ipp], ipp, 0, 0, p);
            Outdm->lcs[Outdm->eog.transtra[sid]-1].give_transfmat(tmat, p, Mp->time);
            glmatrixtransf(t, tmat);
          }
          if (q_n == 6)
            fprintf(out, " % e % e % e % e % e % e\n", t[0][0], t[1][1], t[2][2], t[0][1], t[1][2], t[0][2]);
          else
            fprintf(out, " % e % e % e % e\n", t[0][0], t[1][1], t[0][1], t[2][2]);
        }
        for (j = 0; j < num_copy_ip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0];
          ncompstr = Mm->ip[ipp].ncompstr;
          fprintf(out, "%7c", ' ');
          reallocv(RSTCKVEC(ncompstr, v));
          copyv(Mm->ip[ipp].strain+ncompstr*lcid+q_id1, v);
          vector_tensor (v, t, Mm->ip[ipp].ssst, strain);
          if (Outdm->eog.transtra[sid] > 0)
          {
            ipcoord(Mm->elip[ipp], ipp, 0, 0, p);
            Outdm->lcs[Outdm->eog.transtra[sid]-1].give_transfmat(tmat, p, Mp->time);
            glmatrixtransf(t, tmat);
          }
          if (q_n == 6)
            fprintf(out, " % e % e % e % e % e % e\n", t[0][0], t[1][1], t[2][2], t[0][1], t[1][2], t[0][2]);
          else
            fprintf(out, " % e % e % e % e\n", t[0][0], t[1][1], t[0][1], t[2][2]);
        }
        fprintf(out, "\n");
        break;
      case stress:
        for (j = 0; j < tnip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0]+j;
          if ((te == axisymmlq) && (tnip == 16))
          {
            if(j==1)   ipp=Mt->elements[i].ipp[0][0]+3;
            if(j==2)   ipp=Mt->elements[i].ipp[0][0]+12;
            if(j==3)   ipp=Mt->elements[i].ipp[0][0]+15;
          }
          ncompstr = Mm->ip[ipp].ncompstr;
          if (j == 0)
            fprintf(out, "%7ld", i+Outdm->ide1);
          else
            fprintf(out, "%7c", ' ');
          reallocv(RSTCKVEC(ncompstr, v));
          copyv(Mm->ip[ipp].stress+ncompstr*lcid+q_id1, v);
          vector_tensor (v, t, Mm->ip[ipp].ssst, stress);
          if (Outdm->eog.transtre[sid] > 0)
          {
            ipcoord(Mm->elip[ipp], ipp, 0, 0, p);
            Outdm->lcs[Outdm->eog.transtre[sid]-1].give_transfmat(tmat, p, Mp->time);
            glmatrixtransf(t, tmat);
          }
          if (q_n == 6)
            fprintf(out, " % e % e % e % e % e % e\n", t[0][0], t[1][1], t[2][2], t[0][1], t[1][2], t[0][2]);
          else
            fprintf(out, " % e % e % e % e\n", t[0][0], t[1][1], t[0][1], t[2][2]);
        }
        for (j = 0; j < num_copy_ip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0];
          ncompstr = Mm->ip[ipp].ncompstr;
          fprintf(out, "%7c", ' ');
          reallocv(RSTCKVEC(ncompstr, v));
          copyv(Mm->ip[ipp].stress+ncompstr*lcid+q_id1, v);
          vector_tensor (v, t, Mm->ip[ipp].ssst, stress);
          if (Outdm->eog.transtre[sid] > 0)
          {
            ipcoord(Mm->elip[ipp], ipp, 0, 0, p);
            Outdm->lcs[Outdm->eog.transtre[sid]-1].give_transfmat(tmat, p, Mp->time);
            glmatrixtransf(t, tmat);
          }
          if (q_n == 6)
            fprintf(out, " % e % e % e % e % e % e\n", t[0][0], t[1][1], t[2][2], t[0][1], t[1][2], t[0][2]);
          else
            fprintf(out, " % e % e % e % e\n", t[0][0], t[1][1], t[0][1], t[2][2]);
        }
        fprintf(out, "\n");
        break;
      case other:
        for (j = 0; j < tnip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0]+j;
          if ((te == axisymmlq) && (tnip == 16))
          {
            if(j==1)   ipp=Mt->elements[i].ipp[0][0]+3;
            if(j==2)   ipp=Mt->elements[i].ipp[0][0]+12;
            if(j==3)   ipp=Mt->elements[i].ipp[0][0]+15;
          }
          ncompother = Mm->ip[ipp].ncompeqother;
          if (j == 0)
            fprintf(out, "%7ld", i+Outdm->ide1);
          else
            fprintf(out, "%7c", ' ');
          reallocv(RSTCKVEC(q_n, v));
          copyv(Mm->ip[ipp].eqother+ncompother*lcid+q_id1, v);
          vector_tensor (v, t, guess_ssst(q_n), stress);
          if (q_n == 6)
            fprintf(out, " % e % e % e % e % e % e\n", t[0][0], t[1][1], t[2][2], t[0][1], t[1][2], t[0][2]);
          else
            fprintf(out, " % e % e % e % e\n", t[0][0], t[1][1], t[0][1], t[2][2]);
        }
        for (j = 0; j < num_copy_ip; j++)
        {
          ipp = Mt->elements[i].ipp[0][0];
          ncompother = Mm->ip[ipp].ncompeqother;
          fprintf(out, "%7c", ' ');
          reallocv(RSTCKVEC(q_n, v));
          copyv(Mm->ip[ipp].eqother+ncompother*lcid+q_id1, v);
          vector_tensor (v, t, guess_ssst(q_n), stress);
          if (q_n == 6)
            fprintf(out, " % e % e % e % e % e % e\n", t[0][0], t[1][1], t[2][2], t[0][1], t[1][2], t[0][2]);
          else
            fprintf(out, " % e % e % e % e\n", t[0][0], t[1][1], t[0][1], t[2][2]);
        }
        fprintf(out, "\n");
        break;
      default:
        print_err("unknown value type is required", __FILE__, __LINE__, __func__);
    }
  }
  if (print_header == 0)
    fprintf(out, "End Values\n");
}



/**
  The function exports mesh topology to the opened text file in FemCAD format.

  @param out - pointer to the opened file

  @return The function does not return anything.
  
  Created by Tomas Koudelka 2003
*/
void export_femcad(FILE *out)
{
  write_nodes(out);
  write_elements(out);
}



/**
  The function exports nodal coordinates to the opened text file in FemCAD format.

  @param out - pointer to the opened file

  @return The function does not return anything.
  
  Created by Tomas Koudelka 2003
*/
void write_nodes(FILE *out)
{
  long i;
  fprintf(out, "\nNODES\n");
  for (i=0; i < Gtm->nn; i++)
    fprintf(out, "%ld %e %e %e\n", i, Gtm->gnodes[i].x, Gtm->gnodes[i].y, Gtm->gnodes[i].z);
}



/**
  The function exports elelemnts to the opened text file in FemCAD format.

  @param out - pointer to the opened file

  @return The function does not return anything.
  
  Created by Tomas Koudelka 2003
*/
void write_elements(FILE *out)
{
  long i,j;

  fprintf(out, "\nELEMENTS\n");
  for (i = 0; i < Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;

    switch (Mt->elements[i].te)
    {
      case bar2d:
      case bar3d:
      case beam2d:
      case beam3d:
      case subsoilbeam:
/*        fprintf(out, "%ld %d", i, ISOLinear1D);
        for (j = 0; j < Gtm->gelements[i].nne; j++)
          fprintf(out, " %ld", Gtm->gelements[i].nodes[j]);*/
        break;    
      case barq2d:
      case barq3d:
/*        fprintf(out, "%ld %d", i, ISOQuadratic1D);
        for (j = 0; j < Gtm->gelements[i].nne; j++)
          fprintf(out, " %ld", Gtm->gelements[i].nodes[j]);*/
        break;    
      case spring_1:
      case spring_2:
      case spring_3:
      case spring_4:
      case spring_5:
      case spring_6:
        break;
      case cctel:
      case dktel:
      case dstel:
      case subsoilplatetr:
      case planeelementlt:
      case planeelementrotlt:
        fprintf(out, "%ld %d", i, TriangleLinear);
        for (j = 0; j < Gtm->gelements[i].nne; j++)
          fprintf(out, " %ld", Gtm->gelements[i].nodes[j]);
        break;
      case planeelementqt:
        fprintf(out, "%ld %d", i, TriangleQuadratic);
        for (j = 0; j < Gtm->gelements[i].nne; j++)
          fprintf(out, " %ld", Gtm->gelements[i].nodes[j]);
        break;
      case axisymmqt:
        fprintf(out, "%ld %d", i, TriangleQuadratic);
        for (j = 0; j < Gtm->gelements[i].nne; j++)
          fprintf(out, " %ld", Gtm->gelements[i].nodes[j]);
        break;
      case planeelementlq:
      case planeelementrotlq:
      case q4plateel:
      case subsoilplateq:
        fprintf(out, "%ld %d", i, ISOLinear2D);
        for (j = 0; j < Gtm->gelements[i].nne; j++)
          fprintf(out, " %ld", Gtm->gelements[i].nodes[j]);
        break;
      case planeelementqq:
        fprintf(out, "%ld %d", i, ISOQuadratic2D);
        for (j = 0; j < Gtm->gelements[i].nne; j++)
          fprintf(out, " %ld", Gtm->gelements[i].nodes[j]);
        break;


      case linearhex:
        fprintf(out, "%ld %d", i, ISOLinear2D);
        for (j = 0; j < 4; j++)
          fprintf(out, " %ld", Gtm->gelements[i].nodes[4+j]);
        break;
      default:
        print_err("unsupported element type is required",__FILE__, __LINE__, __func__);
    }
    fprintf(out, "\n");
  }
}



/**
  The function exports nodal coordinates to the opened text file in FemCAD format.
  The source of data is initialized siftop object.

  @param out - pointer to the opened file
  @param st  - pointer to the siftop object with required topology

  @return The function does not return anything.
  
  Created by Tomas Koudelka 2003
*/
void write_nodes_prep(FILE *out, siftop *st)
{
  long i;
  fprintf(out, "\nNODES\n");
  for (i=0; i < st->nn; i++)
    fprintf(out, "%ld %e %e %e\n", i, st->nodes[i].x, st->nodes[i].y, st->nodes[i].z);
}



/**
  The function exports elelemnts to the opened text file in FemCAD format.
  The source of data is initialized siftop object.

  @param out - pointer to the opened file
  @param st  - pointer to the siftop object with required topology

  @return The function does not return anything.
  
  Created by Tomas Koudelka 2003
*/
void write_elements_prep(FILE *out, siftop *st)
{
  long i,j;

  fprintf(out, "\nELEMENTS\n");
  for (i = 0; i < st->ne; i++)
  {
    fprintf(out, "%ld %d", i, TriangleLinear);
    for (j = 0; j < st->elements[i].nne; j++)
      fprintf(out, " %ld", st->elements[i].nodes[j]-1);
    fprintf(out, "\n");
  }
}



/**
  The function writes nodal displacements from the given load case 
  to the opened text file in FemCAD format.
  
  @param out - pointer to the opened file
  @param lcid  - load case id
  @param desclcid - string with load case description

  @return The function does not return anything.
  
  Created by Tomas Koudelka 2003
*/
void write_displ(FILE *out, long lcid, const char *desclcid)
{
  long i, j;
  vector r;
  fprintf (out,"\nVALUES");
  fprintf (out,"\nVECTOR3D NODES DEFORMATION %s\n", desclcid);
  for (i=0; i < Mt->nn; i++)
  {
    reallocv(RSTCKVEC(Gtm->gnodes[i].ndofn, r));
    noddispl (lcid, r.a, i);
    fprintf(out, "%ld", i);
    for (j = 0; j < Gtm->gnodes[i].ndofn; j++)
      fprintf(out, " % e", r[j]);
    fprintf(out, "  0.0 \n");
  }
}



/**
  The function writes nodal forces from the given load case 
  to the opened text file in FemCAD format.
  
  @param out - pointer to the opened file
  @param lcid  - load case id
  @param desclcid - string with load case description
  @param veclabel - string with label of the force %vector
  @param ifor - %vector of nodal force

  @return The function does not return anything.
  
  Created by Tomas Koudelka 2003
*/
void write_nforces(FILE *out, long /*lcid*/, const char *desclcid, const char *veclabel, double *ifor, bool print_react)
{
  long i, j, ii, ndof;
  vector f, g;
  matrix tm(ASTCKMAT(3,3));
  long ndofn, transf;

  fprintf (out,"\nVALUES");
  fprintf (out,"\nVECTOR3D NODES %s %s\n", veclabel, desclcid);
  for (i=0; i < Mt->nn; i++)
  {
    ndof = Mt->give_ndofn(i);
    reallocv(ndof, f);
    for (j=0;j<ndof;j++)
    {
      ii=Mt->give_dof(i,j);
      if ((ii<1) && print_react)  f[j]=Mt->nodes[i].r[j];
      if (ii>0)  f[j]=ifor[ii-1];
    }

    transf = Mt->nodes[i].transf;
    if (transf>0){

      //fprintf (Outm,"\n uzel %ld   %ld",i,Mt->nodes[i].transf);

      ndofn = Gtm->gnodes[i].ndofn;
      reallocv(ndofn, g);
      reallocm(ndofn,ndofn,tm);
      fillm(0.0, tm);
      
      tm[0][0]=Mt->nodes[i].e1[0];
      tm[1][0]=Mt->nodes[i].e1[1];
      if (transf == 3)
        tm[2][0]=Mt->nodes[i].e1[2];
      
      tm[0][1]=Mt->nodes[i].e2[0];
      tm[1][1]=Mt->nodes[i].e2[1];
      if (transf == 3)
        tm[2][1]=Mt->nodes[i].e2[2];
      
      if (transf == 3)
      {
        tm[0][2]=Mt->nodes[i].e3[0];
        tm[1][2]=Mt->nodes[i].e3[1];
        tm[2][2]=Mt->nodes[i].e3[2];
      }
      else
      {
        if (ndofn == 3) // 2D beam
          tm[2][2]=1.0;
      }
      
      mxv(tm,f,g);
      copyv (g,f);
    }
    fprintf(out, "%ld", i);
    if (ndof > 3)   ndof = 3;
    for (j = 0; j < ndof; j++)
      fprintf(out, " % e", f[j]);
    for (j=ndof; j < 3; j++)
      fprintf(out, " 0.0");
    fprintf(out, "\n");
  }
}



/**
  The function writes value of required nodal quantity from the given load case 
  to the opened text file in FemCAD format.
  
  @param out - pointer to the opened file
  @param scal - type of required quantity
  @param lcid  - load case id
  @param dir - component of the required quantity
  @param desclcid - string with load case description

  @return The function does not return anything.
  
  Created by Tomas Koudelka 2003
*/
void write_nodscalar(FILE *out, strastre scal, long lcid, long dir, const char *desclcid)
{
  const char *sig = "";
  long i, ncompstr, ncompother;

  fprintf(out, "\nVALUES\nSCALAR NODES ");
  if (scal == stress)
  {
    if (Mt->nodes[0].ncompstr < 6)
    {
      switch (dir)
      {
        case 0:
          sig = "sig_x";
          break;
        case 1:
          sig = "sig_y";
          break;
        case 2:
          sig = "tau_xy";
          break;
        case 3:
          sig = "sig_z";
          break;
        default:
          print_err("unknown direction of stresses is required", __FILE__, __LINE__, __func__);
      }
    }
    if (Mt->nodes[0].ncompstr == 6)
    {
      switch (dir)
      {
        case 0:
          sig = "sig_x";
          break;
        case 1:
          sig = "sig_y";
          break;
        case 2:
          sig = "sig_z";
          break;
        case 3:
          sig = "tau_yz";
          break;
        case 4:
          sig = "tau_xz";
          break;
        case 5:
          sig = "tau_xy";
          break;
        default:
          print_err("unknown direction of stresses is required", __FILE__, __LINE__, __func__);
      }
    }
  }
  if (scal == strain)
  {
    if (Mt->nodes[0].ncompstr < 6)
    {
      switch (dir)
      {
        case 0:
          sig = "eps_x";
          break;
        case 1:
          sig = "eps_y";
          break;
        case 2:
          sig = "eps_xy";
          break;
        case 3:
          sig = "eps_z";
          break;
        default:
          print_err("unknown direction of strains is required", __FILE__, __LINE__, __func__);
      }
    }
    if (Mt->nodes[0].ncompstr == 6)
    {
      switch (dir)
      {
        case 0:
          sig = "eps_x";
          break;
        case 1:
          sig = "eps_y";
          break;
        case 2:
          sig = "eps_z";
          break;
        case 3:
          sig = "eps_yz";
          break;
        case 4:
          sig = "eps_xz";
          break;
        case 5:
          sig = "eps_xy";
          break;
        default:
          print_err("unknown direction of strains is required", __FILE__, __LINE__, __func__);
      }
    }
  }
  switch (scal)
  {
    case stress:
      fprintf(out, "%s %s\n", sig, desclcid);
      for (i = 0; i < Mt->nn; i++)
      {
        ncompstr = Mt->nodes[i].ncompstr;
        fprintf(out, "%ld % e\n", i, Mt->nodes[i].stress[ncompstr*lcid+dir]);
      }
      break;
    case strain:
      fprintf(out, "%s %s\n", sig, desclcid);
      for (i = 0; i < Mt->nn; i++)
      {
        ncompstr = Mt->nodes[i].ncompstr;
        fprintf(out, "%ld % e\n", i, Mt->nodes[i].strain[ncompstr*lcid+dir]);
      }
      break;
    case other:
      fprintf(out, "other_n%ld %s\n", dir, desclcid);
      for (i = 0; i < Mt->nn; i++)
      {
        ncompother = Mt->nodes[i].ncompother;
        fprintf(out, "%ld % e\n", i, Mt->nodes[i].other[ncompother*lcid+dir]);
      }
      break;
    default:
      print_err("unknown value type is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function writes value of required quantity on elements from the given load case 
  to the opened text file in FemCAD format.
  
  @param out - pointer to the opened file
  @param scal - type of required quantity
  @param lcid  - load case id
  @param dir - component of the required quantity
  @param desclcid - string with load case description

  @return The function does not return anything.
  
  Created by Tomas Koudelka 2003
*/
void write_elemscalar(FILE *out, strastre scal, long lcid, long dir, const char * desclcid)
{
  const char *sig = "";
  long i, ipp,ncompstr,ncompother;

  fprintf(out, "\nVALUES\nSCALAR ELEMENTS ");
  if (scal == stress)
  {
    if (Mm->ip[0].ncompstr < 6)
    {
      switch (dir)
      {
        case 0:
          sig = "sig_x";
          break;
        case 1:
          sig = "sig_y";
          break;
        case 2:
          sig = "tau_xy";
          break;
        case 3:
          sig = "sig_z";
          break;
        default:
          print_err("unknown direction of stresses is required", __FILE__, __LINE__, __func__);
      }
    }
    if (Mm->ip[0].ncompstr == 6)
    {
      switch (dir)
      {
        case 0:
          sig = "sig_x";
          break;
        case 1:
          sig = "sig_y";
          break;
        case 2:
          sig = "sig_z";
          break;
        case 3:
          sig = "tau_yz";
          break;
        case 4:
          sig = "tau_xz";
          break;
        case 5:
          sig = "tau_xy";
          break;
        default:
          print_err("unknown direction of stresses is required", __FILE__, __LINE__, __func__);
      }
    }
  }
  if (scal == strain)
  {
    if (Mm->ip[0].ncompstr < 6)
    {
      switch (dir)
      {
        case 0:
          sig = "eps_x";
          break;
        case 1:
          sig = "eps_y";
          break;
        case 2:
          sig = "eps_xy";
          break;
        case 3:
          sig = "eps_z";
          break;
        default:
          print_err("unknown direction of strains is required", __FILE__, __LINE__, __func__);
      }
    }
    if (Mm->ip[0].ncompstr == 6)
    {
      switch (dir)
      {
        case 0:
          sig = "eps_x";
          break;
        case 1:
          sig = "eps_y";
          break;
        case 2:
          sig = "eps_z";
          break;
        case 3:
          sig = "eps_yz";
          break;
        case 4:
          sig = "eps_xz";
          break;
        case 5:
          sig = "eps_xy";
          break;
        default:
          print_err("unknown direction of strains is required", __FILE__, __LINE__, __func__);
      }
    }
  }
  switch (scal)
  {
    case stress:
      fprintf(out, "%s %s\n", sig, desclcid);
      for (i = 0; i < Mt->ne; i++)
      {
	if (Gtm->leso[i]==0)
	  continue;

        ipp = Mt->elements[i].ipp[0][0];
        //        fprintf(out, "%ld % e\n", i, Mm->ip[ipp].stress[ncompstr*lcid+dir]);
        fprintf(out, "%ld % e\n", i, Mm->ip[ipp].stress[dir]);
      }
      break;
    case strain:
      fprintf(out, "%s %s\n", sig, desclcid);
      for (i = 0; i < Mt->nn; i++)
      {
        ipp = Mt->elements[i].ipp[0][0];
        ncompstr = Mm->ip[ipp].ncompstr;
        fprintf(out, "%ld % e\n", i, Mm->ip[ipp].strain[ncompstr*lcid+dir]);
      }
      break;
    case other:
      fprintf(out, "other_e%ld %s\n", dir, desclcid);
      for (i = 0; i < Mt->nn; i++)
      {
        ipp = Mt->elements[i].ipp[0][0];
        ncompother = Mm->ip[ipp].ncompother;
        fprintf(out, "%ld % e\n", i, Mm->ip[ipp].eqother[ncompother*lcid+dir]);
      }
      break;
    default:
      print_err("unknown value type is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function writes required values at nodes from the given load case 
  to the opened text file in FemCAD format.

  @param out - pointer to the opened file
  @param val - %vector of nodal values
  @param descr - string with quantity description
  @param desclcid - load case id

  @return The function does not return anything.
  
  Created by Tomas Koudelka 2003
*/
void write_nodscalar(FILE *out,double *val, const char *descr, long desclcid)
{
  long i;
  fprintf(out, "\nVALUES\nSCALAR NODES ");
  fprintf(out, "%s %ld\n", descr, desclcid);
  for (i=0; i < Mt->nn; i++)
    fprintf(out, "%ld %e\n", i, val[i]);  
}



/**
  The function writes required values of principal strains/stresses at nodes from the given load case 
  to the opened text file in FemCAD format.

  @param out - pointer to the opened file
  @param dir - indicator of required values 
               (-1, -2, -3 => eps1, eps2, eps3; 0, 1, 2, 3 => psig1, psig2, psig3, tau_max)
  @param desclcid - string with load case description

  @return The function does not return anything.
  
  Created by Tomas Koudelka 2003
*/
void write_nodscalar(FILE *out, long dir, const char *desclcid)
{
  long i;
  fprintf(out, "\nVALUES\nSCALAR NODES ");
  switch (dir)
  {
    case -1:
      fprintf(out, "%s %s\n", "peps_1", desclcid);
      for (i = 0; i < Mt->nn; i++)
      {
        if (Mt->nodes[i].pstra)
          fprintf(out, "%ld % e\n", i, Mt->nodes[i].pstra[-dir-1]);
        else
          fprintf(out, "%ld 0.0\n",i);
      }
      break;
    case -2:
      fprintf(out, "%s %s\n", "peps_2", desclcid);
      for (i = 0; i < Mt->nn; i++)
      {
        if (Mt->nodes[i].pstra)
          fprintf(out, "%ld % e\n", i, Mt->nodes[i].pstra[-dir-1]);
        else
          fprintf(out, "%ld 0.0\n",i);
      }
      break;
    case -3:
      fprintf(out, "%s %s\n", "peps_3", desclcid);
      for (i = 0; i < Mt->nn; i++)
      {
        if (Mt->nodes[i].pstra)
          fprintf(out, "%ld % e\n", i, Mt->nodes[i].pstra[-dir-1]);
        else
          fprintf(out, "%ld 0.0\n",i);
      }
      break;
    case 0:
      fprintf(out, "%s %s\n", "psig_1", desclcid);
      for (i = 0; i < Mt->nn; i++)
      {
        if (Mt->nodes[i].pstre)
          fprintf(out, "%ld % e\n", i, Mt->nodes[i].pstre[dir]);
        else
          fprintf(out, "%ld 0.0\n",i);
      }
      break;
    case 1:
      fprintf(out, "%s %s\n", "psig_2", desclcid);
      for (i = 0; i < Mt->nn; i++)
      {
        if (Mt->nodes[i].pstre)
          fprintf(out, "%ld % e\n", i, Mt->nodes[i].pstre[dir]);
        else
          fprintf(out, "%ld 0.0\n",i);
      }
      break;
    case 2:
      fprintf(out, "%s %s\n", "psig_3", desclcid);
      for (i = 0; i < Mt->nn; i++)
      {
        if (Mt->nodes[i].pstre)
          fprintf(out, "%ld % e\n", i, Mt->nodes[i].pstre[dir]);
        else
          fprintf(out, "%ld 0.0\n",i);
      }
      break;
    case 3 :
      fprintf(out, "%s %s\n", "tau_max", desclcid);
      for (i = 0; i < Mt->nn; i++)
      {
        if (Mt->nodes[i].pstre)
          fprintf(out, "%ld % e\n", i, (Mt->nodes[i].pstre[2]-Mt->nodes[i].pstre[0])/2);
        else
          fprintf(out, "%ld 0.0\n",i);
      }
      break;
    default:
      break;
  }
}



/**
  The function writes required values on elements from the given load case 
  to the opened text file in FemCAD format.

  @param out - pointer to the opened file
  @param val - %vector of values on elements
  @param descr - string with quantity description
  @param desclcid - string with load case description

  @return The function does not return anything.
  
  Created by Tomas Koudelka 2003
*/
void write_elemscalar(FILE *out,double *val, const char *descr, const char *desclcid)
{
  long i;
  fprintf(out, "\nVALUES\nSCALAR ELEMENTS ");
  fprintf(out, "%s %s\n", descr, desclcid);
  for (i=0; i < Mt->ne; i++){
    if (Gtm->leso[i]==0)
      continue;

    fprintf(out, "%ld %e\n", i, val[i]);  
  }
}



/**
  The function writes required component of nodal displacements from the given load case 
  to the opened text file in FemCAD format.

  @param out - pointer to the opened file
  @param lcid  - load case id
  @param dir - component of the required quantity
  @param desclcid - string with load case description

  @return The function does not return anything.
  
  Created by Tomas Koudelka 2003
*/
void write_deflection (FILE *out,long lcid,long dir,const char *desclcid)
{
  long i;
  vector r;
  
  fprintf(out, "\nVALUES\nSCALAR NODES ");
  fprintf(out, "w %s\n", desclcid);
  for (i=0; i < Mt->nn; i++){
    reallocv(RSTCKVEC(Gtm->gnodes[i].ndofn, r));
    noddispl (lcid, r.a, i);
    fprintf(out, "%ld", i);
    fprintf(out, " % e\n", r[dir]);
  }
}


/**
   function prints displacements
   
   @param out - output stream
   @param lcid - load case id
*/
void print_displacements (FILE *out,long lcid)
{
  long i,j,k,ndofn;
  double *r,*d=NULL;
  
  r=Lsrs->lhs+lcid*Ndofm;
  switch (Mp->tprob){
  case linear_statics:{
    d=Mb->lc[lcid].pd;
    break;
  }
  case mat_nonlinear_statics:{
    d=Mb->lc[lcid].pd;
    break;
  }
  case geom_nonlinear_statics:{
    d=Mb->lc[lcid].pd;
    break;
  }
  case earth_pressure:{
    d=Mb->lc[lcid].pd;
    break;
  }
  // case eigen_dynamics:{
  //   break;
  // }
  // case forced_dynamics:{
  //   d=NULL;
  //   break;
  // }
  // case growing_mech_structure:
  // case mech_timedependent_prob:{
  //   d=NULL;
  //   break;
  // }
  case layered_linear_statics:{
    d=Mb->lc[lcid].pd;
    break;
  }
  case lin_floating_subdomain:{
    d=Mb->lc[lcid].pd;
    break;
  }
  case nonlin_floating_subdomain:{
    d=Mb->lc[lcid].pd;
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__, __func__);
    abort();
  }
  }
  
  
  fprintf (out,"\n\n\n NODE DISPLACEMENTS \n");
  for (i=0;i<Mt->nn;i++){
    fprintf (out,"\n\n node %ld",i+1);
    ndofn=Mt->give_ndofn (i);
    for (j=0;j<ndofn;j++){
      k=Gtm->give_dof (i,j);
      if (k>0)   fprintf (out,"\n gen. displ.    %15.12e",r[k-1]);
      if (k==0)  fprintf (out,"\n gen. displ.    %15.12e",0.0);
      if (k<0)   fprintf (out,"\n gen. displ.    %15.12e",d[0-k-1]);
    }
  }

  /*
  r=Lsrs->lhs+lcid*Ndofm;
  d=B->lc[lcid].pd;
  
  k=T->nodes[11].cn[2];
  if (k>0)   fprintf (out,"\n %f",r[k-1]);
  if (k==0)  fprintf (out,"\n %f",0.0);
  if (k<0)   fprintf (out,"\n %f",d[0-k-1]);
  */
}

void print_multipliers (FILE *out)
{
  long i,j,k,l,nl,nmult;
  double *d;
  ivector cn;
  
  d=Lsrs->lhs;
  
  if (Mp->tprob!=layered_linear_statics){
    print_err("there are no multipliers in your problem.", __FILE__, __LINE__, __func__);
    abort ();
  }
  
  if (Mp->tprob==layered_linear_statics){
    fprintf (out,"\n\n");
    for (i=0;i<Mt->nln;i++){
      nl=Gtm->lgnodes[i].nl;
      nmult=Gtm->give_nmult (i);
      reallocv(RSTCKIVEC(nmult, cn));
      fprintf (out,"\n\n layered node %ld",i);
      for (j=1;j<nl;j++){
        Gtm->give_mult_code_numbers (i,j,cn.a);
        for (k=0;k<nmult;k++){
          l=cn[k];
          if (l==0)  fprintf (out,"\n mult.   %15.12e",0.0);
          if (l>1)   fprintf (out,"\n mult.   %15.12e",d[l-1]);
        }
      }
    }
  }
}



void print_strains_old (FILE */*out*/,long /*lcid*/)
{
  /*
  long i,j,k,ii,te,ne,nip,ncomp;
  ne=Mt->ne;
  
  fprintf (out,"\n\n\n ELEMENT STRAINS \n");
  for (i=0;i<ne;i++){
    te = Mt->give_elem_type (i);
    ii=Mt->elements[i].ipp[0];
    nip=Mt->give_nip(i,0);
    ncomp=Mt->give_ncomp(i);

    switch (te){
    case bar2d:{
      fprintf (out,"\n\n 2D bar element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].strain[k]);
        }
      }
      break;
    }
    case beam2d:{
      fprintf (out,"\n\n 2D beam element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].strain[k]);
        }
      }
      break;
    }

    case planeelementlt:{
      fprintf (out,"\n\n plane lt element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].strain[k]);
        }
      }
      break;
    }
    case planeelementqt:{
      fprintf (out,"\n\n plane qt element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].strain[k]);
        }
      }
      break;
    }
    case planeelementrotlt:{
      fprintf (out,"\n\n plane lt element with rotational degrees of freedom");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].strain[k]);
        }
      }
      break;
    }
    case planeelementlq:{
      fprintf (out,"\n\n plane lq element %ld",i+1);
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].strain[k]);
        }
      }
      break;
    }
    case planeelementqq:{
      fprintf (out,"\n\n plane qq element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].strain[k]);
        }
      }
      break;
    }
    case planeelementrotlq:{
      fprintf (out,"\n\n plane lq element with rotational degree of freedom");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].strain[k]);
        }
      }
      break;
    }

    case cctel:{
      fprintf (out,"\n\n cct element");
      for (j=ii;j<ii+1;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<3;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].strain[k]);
        }
      }
      for (j=ii+1;j<ii+4;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<2;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].strain[k+3]);
        }
      }
      break;
    }
    case axisymmlq:{
      fprintf (out,"\n\n linear quadrilateral element for axisymmetric problems");
      for (j=ii;j<ii+4;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].strain[k]);
        }
      }
      fprintf (out,"\n %4ld",ii+1);
      for (k=0;k<ncomp;k++){
        fprintf (out,"   %20.10f",Mm->ip[ii+4].strain[k]);
      }
      break;
    }
    case lineartet:{
      fprintf (out,"\n\n linear tetrahedral element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].strain[k]);
        }
      }
      break;
    }
    case linearhex:{
      fprintf (out,"\n\n linear hexagonal element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].strain[k]);
        }
      }
      break;
    }
    case quadrhex:{
      fprintf (out,"\n\n quadratic hexagonal element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].strain[k]);
        }
      }
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown element type is required in function");
      fprintf (stderr,"\n print_strains (%s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
  */
}




/**
   function prints strains computed in nodes
   
   @param out - output file
   @param eid - element id
   
   19.5.2002
*/
void print_strains_nodes (FILE *out,long lcid)
{
  long i,j,nn,ncomp;
  nn = Mt->nn;
  
  fprintf (out,"\n\n\n\n STRAINS IN NODES IN LOAD CASE   %ld\n",lcid+1);
  for (i=0;i<nn;i++){
    ncomp=Mt->nodes[i].ncompstr;
    fprintf (out,"\n node %4ld",i+1);
    for (j=0;j<ncomp;j++){
      fprintf (out," %e",Mt->nodes[i].strain[lcid*ncomp+j]);
    }
  }
}

/**
   function prints strains computed in user defined points
   
   @param out - output file
   @param eid - element id
   
   23.2.2002
*/
void print_strains_udp (FILE *out,long eid)
{
  long i,j,nape,ncomp;
  
  nape=Mm->stra.nape[eid];
  ncomp=Mt->give_tncomp(eid);
  fprintf (out,"\n\n strains on element number %ld",eid);
  for (i=0;i<nape;i++){
    fprintf (out,"\n point number %ld\n",i);
    for (j=0;j<ncomp;j++){
      fprintf (out," %20.10f",Mm->stra.ev[eid][i][j]);
    }
  }
}

/**
   function prints strains computed in integration points
   
   @param out - output file
   @param lcid - load case id
   
   19.5.2002
*/
/*
void print_strains_intp (FILE *out,long lcid)
{
  long i,j,k,ii,ipp,ne,nb,ncomp;

  fprintf (out,"\n\n\n\n STRAINS IN INTEGRATION POINTS IN LOAD CASE   %ld\n",lcid+1);
  ne=Mt->ne;
  for (i=0;i<ne;i++){
    fprintf (out,"\n element %ld",i+1);
    nb = Mt->give_nb (i);
    for (ii=0;ii<nb;ii++){
      fprintf (out,"\n block %ld",ii+1);
      ipp=Mt->elements[i].ipp[ii][ii];
      for (j=0;j<Mt->elements[i].nip[ii][ii];j++){
        fprintf (out,"\n int. point %ld",j+1);
        ncomp=Mm->ip[ipp].ncompstr;
        for (k=0;k<ncomp;k++){
          fprintf (out," %e",Mm->ip[ipp].strain[lcid*ncomp+k]);
        }
        ipp++;
      }
    }
  }
}
*/

/**
   function prints strains on elements
   
   @param out - output stream
   
   23.2.2002
*/
/*
void print_strains (FILE *out,long lcid)
{
  print_strains_nodes (out,lcid);
  print_strains_intp (out,lcid);
}
*/

/**
   function prints stresses computed in nodes
   
   @param out - output file
   @param eid - element id
   
   19.5.2002
*/
void print_stresses_nodes (FILE *out,long lcid)
{
  long i,j,nn,ncomp;
  nn = Mt->nn;
  
  fprintf (out,"\n\n\n\n STRESSES IN NODES IN LOAD CASE   %ld\n",lcid+1);
  for (i=0;i<nn;i++){
    ncomp=Mt->nodes[i].ncompstr;
    fprintf (out,"\n node %4ld",i+1);
    for (j=0;j<ncomp;j++){
      fprintf (out," %e",Mt->nodes[i].stress[lcid*ncomp+j]);
    }
  }
}

/**
   function prints stresses computed in user defined points
   
   @param out - output file
   @param eid - element id
   
   23.2.2002
*/
void print_stresses_udp (FILE *out,long eid)
{
  long i,j,nape,ncomp;
  
  nape=Mm->stra.nape[eid];
  ncomp=Mt->give_tncomp(eid);
  fprintf (out,"\n\n stresses on element number %ld",eid);
  for (i=0;i<nape;i++){
    fprintf (out,"\n point number %ld\n",i);
    for (j=0;j<ncomp;j++){
      fprintf (out," %20.10f",Mm->stre.ev[eid][i][j]);
    }
  }
}

/**
   function prints stresses computed in integration points
   
   @param out - output file
   @param eid - element id
   
   23.2.2002
*/
/*
void print_stresses_intp (FILE *out,long lcid)
{
  long i,j,k,ii,ipp,ne,nb,ncomp;
  
  fprintf (out,"\n\n\n\n STRESSES IN INTEGRATION POINTS IN LOAD CASE   %ld\n",lcid+1);
  ne=Mt->ne;
  for (i=0;i<ne;i++){
    fprintf (out,"\n element %ld",i+1);
    nb = Mt->give_nb (i);
    for (ii=0;ii<nb;ii++){
      fprintf (out,"\n block %ld",ii+1);
      ipp=Mt->elements[i].ipp[ii][ii];
      for (j=0;j<Mt->elements[i].nip[ii][ii];j++){
        fprintf (out,"\n int. point %ld",j+1);
        ncomp=Mm->ip[ipp].ncompstr;
        for (k=0;k<ncomp;k++){
          fprintf (out," %e",Mm->ip[ipp].stress[lcid*ncomp+k]);
        }
        ipp++;
      }
    }
  }
}
*/

/**
   function prints stresses on elements
   
   @param out - output stream
   
   19.5.2002
*/
/*
void print_stresses (FILE *out,long lcid)
{
  print_stresses_nodes (out,lcid);
  print_stresses_intp (out,lcid);
}
*/






/**
   function prints stresses computed in integration points
   
   @param out - output stream
   @param lcid - load case id
*/
void print_stresses_old (FILE */*out*/,long /*lcid*/)
{
  /*
  long i,j,k,ii,te,ne,nip,ncomp;
  ne=Mt->ne;
  
  fprintf (out,"\n\n\n ELEMENT STRESSES \n");
  for (i=0;i<ne;i++){
    te = Mt->give_elem_type (i);
    ii=Mt->elements[i].ipp[0];
    nip=Mt->give_nip (i,0);
    ncomp=Mt->give_ncomp (i);

    switch (te){
    case bar2d:{
      fprintf (out,"\n\n 2D bar element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].stress[k]);
        }
      }
      break;
    }
    case beam2d:{
      fprintf (out,"\n\n 2D beam element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].stress[k]);
        }
      }
      break;
    }

    case planeelementlt:{
      fprintf (out,"\n\n plane lt element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].stress[k]);
        }
      }
      break;
    }
    case planeelementqt:{
      fprintf (out,"\n\n plane qt element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].stress[k]);
        }
      }
      break;
    }
    case planeelementrotlt:{
      fprintf (out,"\n\n plane lt element with rotational degrees of freedom");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].stress[k]);
        }
      }
      break;
    }
    case planeelementlq:{
      fprintf (out,"\n\n plane lq element %ld",i+1);
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].stress[k]);
        }
      }
      break;
    }
    case planeelementqq:{
      fprintf (out,"\n\n plane qq element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].stress[k]);
        }
      }
      break;
    }
    case planeelementrotlq:{
      fprintf (out,"\n\n plane lq element with rotational degree of freedom");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].stress[k]);
        }
      }
      break;
    }

    case cctel:{
      fprintf (out,"\n\n cct element");
      for (j=ii;j<ii+1;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<3;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].stress[k]);
        }
      }
      for (j=ii+1;j<ii+4;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<2;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].stress[k+3]);
        }
      }
      break;
    }
    case axisymmlq:{
      fprintf (out,"\n\n linear quadrilateral element for axisymmetric problems");
      for (j=ii;j<ii+4;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].stress[k]);
        }
      }
      fprintf (out,"\n %4ld",ii+1);
      for (k=0;k<ncomp;k++){
        fprintf (out,"   %20.10f",Mm->ip[ii+4].stress[k]);
      }
      break;
    }
    case lineartet:{
      fprintf (out,"\n\n linear tetrahedral element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].stress[k]);
        }
      }
      break;
    }
    case linearhex:{
      fprintf (out,"\n\n linear hexagonal element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].stress[k]);
        }
      }
      break;
    }
    case quadrhex:{
      fprintf (out,"\n\n quadratic hexagonal element");
      for (j=ii;j<ii+nip;j++){
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].stress[k]);
        }
      }
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown element type is required in function");
      fprintf (stderr,"\n print_stresses (%s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
  */
}

void print_other (FILE *out,long /*lcid*/)
  //  nedodelano, spatny pocet slozek ncomp!
{
  long i,j;
  
  fprintf (out,"\n\n\n OTHER \n\n");
  for (i=0;i<Mm->tnip;i++){
    fprintf (out,"\n");
    for (j=0;j<Mm->ip[i].ncompother;j++){
      fprintf (out,"  %f",Mm->ip[i].eqother[j]);
    }
  }
  /*
  long i,j,k,ii,te,ne,nip,ncomp;
  ne=Mt->ne;
  
  fprintf (out,"\n\n\n ELEMENT OTHER ARRAY \n");
  for (i=0;i<ne;i++){
    te = Mt->give_elem_type (i);
    ii=Mt->elements[i].ipp[0];
    nip=Mt->give_nip (i,0);
    switch (te){
    case bar2d:{
      fprintf (out,"\n\n bar 2D element");
      for (j=ii;j<ii+nip;j++){
        ncomp=Mm->ip[j].ncompother;
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].eqother[k]);
        }
      }
      break;
    }
    case beam2d:{
      fprintf (out,"\n\n beam 2D element");
      for (j=ii;j<ii+nip;j++){
        ncomp=Mm->ip[j].ncompother;
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].eqother[k]);
        }
      }
      break;
    }

    case planeelementlq:{
      fprintf (out,"\n\n plane lq element %ld",i+1);
      for (j=ii;j<ii+nip;j++){
        ncomp=Mm->ip[j].ncompother;
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].eqother[k]);
        }
      }
      break;
    }
    case linearhex:{
      fprintf (out,"\n\n linear hexagonal element");
      for (j=ii;j<ii+nip;j++){
        ncomp=Mm->ip[j].ncompother;
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].eqother[k]);
        }
      }
      break;
    }
    case quadrhex:{
      fprintf (out,"\n\n quadratic hexagonal element");
      for (j=ii;j<ii+nip;j++){
        ncomp=Mm->ip[j].ncompother;
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].eqother[k]);
        }
      }
      break;
    }
    case lineartet:{
      fprintf (out,"\n\n linear tetrahedral element");
      for (j=ii;j<ii+nip;j++){
        ncomp=Mm->ip[j].ncompother;
        fprintf (out,"\n %4ld",j+1);
        for (k=0;k<ncomp;k++){
          fprintf (out,"   %20.10f",Mm->ip[j].eqother[k]);
        }
      }
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown element type is required in function");
      fprintf (stderr,"\n print_other (%s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
  */
}


void print_intforces (FILE *out,double *fi)
{
  long i,j,k,ndofn;
  
  fprintf (out,"\n\n\n NODE INTERNAL FORCES \n");
  for (i=0;i<Mt->nn;i++){
    fprintf (out,"\n\n node %ld  ",i+1);
    ndofn=Mt->give_ndofn (i);
    for (j=0;j<ndofn;j++){
      k=Gtm->give_dof (i,j);
      if (k>0)   fprintf (out," %15.12e",fi[k-1]);
      else fprintf (out,"%15s", "reakce");
    }
    fprintf (out,"\n");
  }

}


void print_forces (FILE *out,double *f)
{
  long i, j, mnid, ndofn;
  vector nf;
  
  for (i=0;i<Mt->nn;i++){
    fprintf (out," %ld  ",i+1);
    if (Gtm->lnso[i] == 0)
    {
      fprintf (out,"\n");
      continue;
    }
    ndofn=Mt->give_ndofn (i);
    if (ndofn<0){
      mnid   = Gtm->gnodes[i].mnodes[0];
      ndofn = Gtm->give_ndofn(mnid);
    }
    reallocv(RSTCKVEC(ndofn, nf));
    nodforce (f, nf.a, i);    
    for (j=0;j<ndofn;j++)
      fprintf (out," % 15.12e",nf[j]);
    fprintf (out,"\n");
  }
}


void print_reactions (FILE *out,long /*lcid*/)
{
  long i,j,nn,ndofn;
  nn=Mt->nn;

  fprintf (out,"\n\n\n REACTIONS \n");
  for (i=0;i<nn;i++){
    if (Mt->nodes[i].react==1){
      ndofn=Mt->give_ndofn (i);
      fprintf (out,"\n%4ld",i+1);
      for (j=0;j<ndofn;j++){
        fprintf (out,"   %20.10f",Mt->nodes[i].r[j]);
      }
    }
  }
}

void print_eigenvalues (double *w)
{
  long i,n;
  double pi=3.14159265358979;

  n=Mp->eigsol.neigv;
  
  fprintf (Out,"\n\n\n EIGENVALUES \n");
  //fprintf (Out,"\n");
  for (i=0;i<n;i++){
    //fprintf (Out,"\n omega %4ld      %30.20f   %30.20f",i+1,w[i],sqrt(w[i]));
    fprintf (Out,"\n %4ld      %30.15f   %30.15f   %30.15f",i+1,w[i],sqrt(w[i]),sqrt(w[i])/2.0/pi);
  }

}

void print_eigenvectors ()
{
  long i,j,k,n,m;
  double *v;

  n=Mp->eigsol.neigv;
  m=Ndofm;

  v=Lsrs->give_lhs (0);
  

  fprintf (Out,"\n\n\n EIGENVECTORS \n");
  fprintf (Out,"\n\n");
  k=0;
  for (i=0;i<n;i++){
    fprintf (Out,"\n\n Eigenvector number %ld",i+1);
    for (j=0;j<m;j++){
      fprintf (Out,"\n %e",v[k]);
      k++;
    }
  }


  
  /*
  fprintf (Out,"\n\n\n EIGENVECTORS \n");
  fprintf (Out,"\n\n");
  k=0;
  for (i=0;i<n;i++){
    fprintf (Out,"\n\n Eigenvector number %ld  x components",i+1);
    for (j=0;j<Mt->nn;j++){
      k=Mt->give_dof (j,0);
      if (k<1)  
	fprintf (Out,"\n %e 0.0",Gtm->gnodes[j].z);
      else
	fprintf (Out,"\n %e %e",Gtm->gnodes[j].z,v[i*m+k-1]);
    }

    fprintf (Out,"\n\n Eigenvector number %ld  z components",i+1);
    for (j=0;j<Mt->nn;j++){
      k=Mt->give_dof (j,1);
      if (k<1)  
	fprintf (Out,"\n %e 0.0",Gtm->gnodes[j].z);
      else
	fprintf (Out,"\n %e %e",Gtm->gnodes[j].z,v[i*m+k-1]);
    }
  }
  

  long ndofe;
  vector ifor;
  fprintf (Out,"\n\n\n EIGENMOMENTS \n");
  fprintf (Out,"\n\n");
  for (j=0;j<n;j++){
    fprintf (Out,"\n\n Eigennormalforces number %ld",j+1);
    for (i=0;i<Mt->ne;i++){
      ndofe=Mt->give_ndofe (i);
      reallocv (RSTCKVEC(ndofe,ifor));
      
      Beam2d->res_internal_forces (j,i,ifor);
      
      fprintf (Out,"\n %e %e",Gtm->gnodes[i].z,ifor[0]);
      
    }
    ndofe=Mt->give_ndofe (i-1);
    reallocv (RSTCKVEC(ndofe,ifor));
    
    Beam2d->res_internal_forces (j,i-1,ifor);
    
    fprintf (Out,"\n %e %e",Gtm->gnodes[i].z,-1.0*ifor[3]);
    
    
    fprintf (Out,"\n\n Eigenmoments number %ld",j+1);
    for (i=0;i<Mt->ne;i++){
      ndofe=Mt->give_ndofe (i);
      reallocv (RSTCKVEC(ndofe,ifor));
      
      Beam2d->res_internal_forces (j,i,ifor);
      
      fprintf (Out,"\n %e %e",Gtm->gnodes[i].z,ifor[2]);
      
    }
    ndofe=Mt->give_ndofe (i-1);
    reallocv (RSTCKVEC(ndofe,ifor));
    
    Beam2d->res_internal_forces (j,i-1,ifor);
    
    fprintf (Out,"\n %e %e",Gtm->gnodes[i].z,ifor[5]);
    
  }
  */
  
  /*
  double koef,*av,*u;
  av = new double [Ndofm];
  u = new double [Ndofm];
  nullv (av,Ndofm);
  
  for (i=0;i<Mt->nn;i++){
    j=Mt->give_dof (i,0);
    if (j>0)
      av[j-1]=1.0;
  }
  
  fprintf (Out,"\n\n\n koeficienty");
  for (i=0;i<Mp->eigsol.neigv;i++){
    v=Lsrs->give_lhs (i);
    
    Mmat->gmxv (v,u);
    koef = ss (av,u,Ndofm);
    
    fprintf (Out,"\n koeficient  %ld   %e",i,koef);
  }
  
  delete [] av;
  delete [] u;
  
  */
}


void print_eigenvect_martin (FILE *out)
{
  long i,j,n,ndofn;
  double *v;
  
  /*
  n=Mp->neigv;

  for (i=0;i<n;i++){
    v=Lsrs->give_lhs (i);
    fprintf (Out,"\n\n Eigenvector number %ld",i+1);
    for (j=0;j<Mt->nn;j++){
      fprintf (Out,"\n%ld",j);
      ndofn=Mt->give_ndofn (j);
      if (ndofn==2){
        if (Gtm->give_dof (j,0)>0)  fprintf (Out," %e",v[Gtm->give_dof (j,0)-1]);
        else  fprintf (Out," 0.0");
        if (Gtm->give_dof (j,1)>0)  fprintf (Out," %e",v[Gtm->give_dof (j,1)-1]);
        else  fprintf (Out," 0.0");
        fprintf (Out," 0.0");
      }
      if (ndofn==3){
        if (Gtm->give_dof (j,0)>0)  fprintf (Out," %e",v[Gtm->give_dof (j,0)-1]);
        else  fprintf (Out," 0.0");
        if (Gtm->give_dof (j,1)>0)  fprintf (Out," %e",v[Gtm->give_dof (j,1)-1]);
        else  fprintf (Out," 0.0");
        if (Gtm->give_dof (j,2)>0)  fprintf (Out," %e",v[Gtm->give_dof (j,2)-1]);
        else  fprintf (Out," 0.0");
      }
    }
  }
  */
  
  n=Mp->eigsol.neigv;

  for (i=0;i<n;i++){
    v=Lsrs->give_lhs (i);
    fprintf (out,"\n\n Eigenvector number %ld",i+1);
    for (j=0;j<Mt->nn;j++){
      fprintf (out,"\n%ld",j);
      ndofn=Mt->give_ndofn (j);
      if (ndofn==2){
        if (Gtm->give_dof (j,0)>0)  fprintf (out," %e",v[Gtm->give_dof (j,0)-1]);
        else  fprintf (out," 0.0");
        if (Gtm->give_dof (j,1)>0)  fprintf (out," %e",v[Gtm->give_dof (j,1)-1]);
        else  fprintf (out," 0.0");
        fprintf (out," 0.0");
      }
      if (ndofn==3){
        if (Gtm->give_dof (j,0)>0)  fprintf (out," %e",v[Gtm->give_dof (j,0)-1]);
        else  fprintf (out," 0.0");
        if (Gtm->give_dof (j,1)>0)  fprintf (out," %e",v[Gtm->give_dof (j,1)-1]);
        else  fprintf (out," 0.0");
        if (Gtm->give_dof (j,2)>0)  fprintf (out," %e",v[Gtm->give_dof (j,2)-1]);
        else  fprintf (out," 0.0");
      }
    }
  }
  
}

/*
void print_edges(FILE *out)
{
  long i, j , ned;
  fprintf(out, "GLOBAL NUMBERS OF EDGES ON ELEMENTS\n");
  for (i = 0; i < Gtm->ne; i++)
  {
    ned = Mt->give_ned(i);
    fprintf(out, "Element id %5ld has edges:", i+1);
    for (j = 0; j < ned; j++){
      //fprintf(out, " %5ld", Gtm->gelements[i].edgn[j]+1);
      fprintf(out, "\n chyba v mechprint, line 4972\n");
    }
    fprintf(out, "\n");
  }
  fprintf(out, "\n\n***************************");

  fprintf(out, "TOTAL NUMBER OF EDGES %ld\n", Mt->tned);
  for (i = 0; i < Mt->tned; i++)
    fprintf(out, "Edge id %5ld connects elements %5ld and %5ld\n", i+1, Mt->eadjelem[i][0]+1, Mt->eadjelem[i][1]+1);

  fprintf(out, "\n\n***************************");
  fprintf(out, "ADJACENT ELEMENTS\n");
  for (i = 0; i < Mt->ne; i++)
  {
    ned = Mt->give_ned(i);
    fprintf(out, "Element id %5ld connects these elements:", i+1);
    for (j = 0; j < ned; j++)
      fprintf(out, " %5ld", Mt->adjelem[i][j]+1);
    fprintf(out, "\n");
  }


}
*/

//  ////////////////////       /* termitovo */       ////////////////////////////////////

/**
   function compiles data(from array 'valel') for printing file.dx(.ex)
   
   @param gt - gtopology
   @param mp - probdesc
   @param file - name of "dx"("ex") file; if file=NULL file is generated automatically from 'caption'
   @param caption - name=caption of 'value'
   @param valel - array of 'values' on elements, one element = one value, => dimension is gt->ne
   @param flag - type of file - "d" = file.dx(for OpenDx) , "e" = file.ex(for Elixir)
   
   created  20.11.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void print_valel (gtopology *gt, probdesc *mp, const char *file, char *caption, double *valel,char flag)
{
  long i,j;
  double *valnod;
  char *fff;
  char *tcaption;

  // due to literals passed in caption
  tcaption = new char[strlen(caption)];
  i = 0;
  while (caption[i] != '\0')
  {
    tcaption[i] = caption[i];
    i++;
  }
  tcaption[i] = '\0';
  
  fff = NULL;
  valnod = new double [gt->nn];
  
  if (file == NULL){
    fff = new char [255];
    if (flag == 'e')  sprintf (fff,"%s%s.ex",mp->path,tcaption);
    if (flag == 'd')  sprintf (fff,"%s%s.dx",mp->path,tcaption);
    
    file = fff;
  }
  
  if (gt->nadjelnod == NULL)
    gt->adjacelem (Out);
  
  for (i=0;i<gt->nn;i++){
    valnod[i] = 0.0;
    for (j=0;j<gt->nadjelnod[i];j++){
      valnod[i] += valel[gt->adjelnod[i][j]];
    }
    valnod[i] /= (double)gt->nadjelnod[i];
  }
  
  if (flag == 'e')
    print_ex (gt,file,valnod,valel);
  
  if (flag == 'd'){
    long dimindex = 1;
    print_dx (gt,file,&valnod,&valel,'n',&dimindex,&tcaption,1);
  }
  
  if (fff != NULL){ file = NULL; delete [] fff; }
  delete [] valnod;
  delete [] tcaption;
}
//  
//  /**
//     function compiles data(from array 'valnod') for printing file.dx(.ex)
//     
//     @param gt - gtopology
//     @param mp - probdesc
//     @param mt - mechtop
//     @param file - name of "dx"("ex") file; if file=NULL file is generated automatically from 'cap'
//     @param caption - name=caption of 'value'
//     @param valnod - array of 'values' on nodes, one node = one value, => dimension is gt->nn
//     @param flag - type of file - "d" = file.dx(for OpenDx) , "e" = file.ex(for Elixir)
//     
//     created  5.2.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  void print_valnod (gtopology *gt, probdesc *mp, mechtop *mt, const char *file, char *caption, double *valnod,char flag)
//  {
//    long i,j;
//    char *fff;
//    double *valel = NULL;
//    vector nodval;
//    char *tcaption;
//  
//    // due to literals passed in caption
//    tcaption = new char[strlen(caption)];
//    i = 0;
//    while (caption[i] != '\0')
//    {
//      tcaption[i] = caption[i];
//      i++;
//    }
//    tcaption[i] = '\0';
//    
//    
//    if (file == NULL){
//      fff = new char [255];
//      if (flag == 'e')
//        sprintf (fff,"%s%s.ex",mp->path,tcaption);
//      if (flag == 'd')
//        sprintf (fff,"%s%s.dx",mp->path,tcaption);
//      
//      file = fff;
//    }
//    
//    for (i=0;i<gt->ne;i++)
//      if (mt->elements[i].te == planeelementqq){
//        if (valel==NULL){
//          valel = new double [gt->ne];
//          reallocv (RSTCKVEC(8,nodval));
//        }
//        
//        for (j=0;j<8;j++)
//          nodval[j] = valnod[gt->gelements[i].nodes[j]];
//        
//        valel[i] = Peqq->approx (0.0,0.0,nodval);
//      }
//    
//    if (flag == 'e')
//      print_ex (gt,file,valnod,valel);
//    
//    if (flag == 'd'){
//      long dimindex = 1;
//      print_dx (gt,file,&valnod,&valel,'n',&dimindex,&tcaption,1);
//    }
//    
//    if (valel!=NULL){
//      delete [] valel;
//    }
//  }
//  
/**
   function compiles data for printing file.dx
   data = deformation(from Lsrs), all components of strain and stress(from mt->nodes)
   
   @param gt - gtopology
   @param mp - probdesc
   @param mt - mechtop
   @param mm - mechmat
   @param lcid - id of load case
   @param file - name of "dx" file
   
   created  20.11.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
 */
void print_default_dx (gtopology *gt,probdesc *mp,mechtop *mt,mechmat */*mm*/,long lcid,const char *file)
{
  long i,j,nn,ne,vali,capi,ndofn=0,ncomp,dd;
  double **valnod,**valel,*r;
  char **caption;
  long d[6],*dimindex;
  
  dd = 1;  //  Decomposition of Deformation 0=no 1=yes
  
  if (!mp->detnodstrain){
    //mm->compute_nodestrains (lcid);
    //mm->compute_nodestresses (lcid);
  }
  
  switch (mt->elements[0].te){    //    epsxx    epsyy    epszz    epsyz    epsxz    epsxy
  case planeelementlt:
  case planeelementlq:
  case planeelementqt: { ndofn = 2; d[0]= 1; d[1]= 1; d[2]= 0; d[3]= 0; d[4]= 0; d[5]= 1; break;}
  case planeelementqq: { ndofn = 2; d[0]= 1; d[1]= 1; d[2]= 0; d[3]= 0; d[4]= 0; d[5]= 1; break;}
  case lineartet:
  case quadrtet:
  case linearhex:      { ndofn = 3; d[0]= 1; d[1]= 1; d[2]= 1; d[3]= 1; d[4]= 1; d[5]= 1; break;}
  default:{
    print_err("unknown element type is required", __FILE__,__LINE__, __func__);
    abort();
  }
  }
  
  vali=0;
  nn = gt->nn;
  ne = gt->ne;
  ncomp = mt->give_tncomp (0);
  
  r = new double [ndofn];
  valel = new double*  [ne];
  valnod = new double* [(1+dd)*ndofn+2*ncomp];
  caption = new char*  [ 1+dd *ndofn+2*ncomp];
  dimindex = new long  [ 1+dd *ndofn+2*ncomp];
  
  
  // deformation
  for (i=0;i<ndofn;i++){
    valnod[vali+i] = new double [nn];
    if (dd) valnod[vali+i+ndofn] = new double [nn];
  }
  
  for (i=0;i<nn;i++){
    noddispl (lcid,r,i);
    for (j=0;j<ndofn;j++){
      valnod[vali+j][i] = r[j];
      if (dd) valnod[vali+j+ndofn][i] = r[j];
    }
  }
  vali+=(1+dd)*ndofn;
  
  // strain components
  for (i=0;i<ncomp;i++){
    valnod[vali] = new double [nn];
    for (j=0;j<nn;j++)
      valnod[vali][j] = mt->nodes[j].strain[lcid*ncomp+i];
    
    vali++;
  }
  
  // stress components
  for (i=0;i<ncomp;i++){
    valnod[vali] = new double [nn];
    for (j=0;j<nn;j++)
      valnod[vali][j] = mt->nodes[j].stress[lcid*ncomp+i];
    
    vali++;
  }
  
  // valel
  for (i=0;i<ne;i++)
    if (planeelementqq == mt->elements[i].te){
      valel[i] = new double[vali];
      Peqq->midpoints  (i,dd,valel[i]);
    }
    else
      valel[i] = NULL;
  
  
  capi=0;          dimindex[capi]=ndofn;  caption[capi] = new char[4];   sprintf (caption[capi],"def");  capi++;
  
  if (dd && ndofn>0) { dimindex[capi]=1;  caption[capi] = new char[5];   sprintf (caption[capi],"defx");  capi++; }
  if (dd && ndofn>1) { dimindex[capi]=1;  caption[capi] = new char[5];   sprintf (caption[capi],"defy");  capi++; }
  if (dd && ndofn>2) { dimindex[capi]=1;  caption[capi] = new char[5];   sprintf (caption[capi],"defz");  capi++; }
  
  if (d[0])  { dimindex[capi]=1;  caption[capi] = new char[6];   sprintf (caption[capi],"epsxx");  capi++; }
  if (d[1])  { dimindex[capi]=1;  caption[capi] = new char[6];   sprintf (caption[capi],"epsyy");  capi++; }
  if (d[2])  { dimindex[capi]=1;  caption[capi] = new char[6];   sprintf (caption[capi],"epszz");  capi++; }
  if (d[3])  { dimindex[capi]=1;  caption[capi] = new char[6];   sprintf (caption[capi],"epsyz");  capi++; }
  if (d[4])  { dimindex[capi]=1;  caption[capi] = new char[6];   sprintf (caption[capi],"epsxz");  capi++; }
  if (d[5])  { dimindex[capi]=1;  caption[capi] = new char[6];   sprintf (caption[capi],"epsxy");  capi++; }
  
  if (d[0])  { dimindex[capi]=1;  caption[capi] = new char[6];   sprintf (caption[capi],"sigxx");  capi++; }
  if (d[1])  { dimindex[capi]=1;  caption[capi] = new char[6];   sprintf (caption[capi],"sigyy");  capi++; }
  if (d[2])  { dimindex[capi]=1;  caption[capi] = new char[6];   sprintf (caption[capi],"sigzz");  capi++; }
  if (d[3])  { dimindex[capi]=1;  caption[capi] = new char[6];   sprintf (caption[capi],"sigyz");  capi++; }
  if (d[4])  { dimindex[capi]=1;  caption[capi] = new char[6];   sprintf (caption[capi],"sigxz");  capi++; }
  if (d[5])  { dimindex[capi]=1;  caption[capi] = new char[6];   sprintf (caption[capi],"sigxy");  capi++; }
  
  print_dx (gt,file,valnod,valel,'t',dimindex,caption,capi);
  
  
  delete [] r;
  delete [] dimindex;
  for (i=0;i<vali;i++)
    delete [] valnod[i];
  delete [] valnod;
  for (i=0;i<capi;i++)
    delete [] caption[i];
  delete [] caption;
  for (i=0;i<ne;i++){
    if (valel[i]!=NULL)
      delete [] valnod[i];
  }
  delete [] valel;
}

//  /**
//     function compiles data for printing file.dx
//     data = deformation(from Lsrs), all components from arrays Mt->strain, Mt->stress and Mt->eqother
//     
//     @param gt - gtopology
//     @param mp - probdesc
//     @param mt - mechtop
//     @param lcid - id of load case
//     @param file - name of "dx" file
//     
//     created  20.4.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
//   */
//  void print_default_2_dx (gtopology *gt,probdesc *mp,mechtop *mt,long lcid,const char *file)
//  {
//    long i,j,nn,ne,capi,ndofn,dim,nipval,dd;
//    long *dimindex,nipcomp[3];
//    double *nodvals,**valnod,*r;
//    char **caption;
//    
//    dd = 1;  //  Decomposition of Deformation 0=no 1=yes
//    
//    // interpolation values from int. points to nodes
//    nipcomp[0] = nipcomp[1] = nipcomp[2] = -1;
//    dim = mt->give_dimension (0);
//    
//    nipval = give_nodvals_ip (dim,nipcomp,nodvals);
//    
//    // 
//    nn = gt->nn;
//    ne = gt->ne;
//    ndofn = mt->give_ndofn (0);
//    
//    r = new double [ndofn];
//    dimindex = new long  [ 1+dd *ndofn+nipval];
//    caption = new char*  [ 1+dd *ndofn+nipval];
//    valnod = new double* [(1+dd)*ndofn+nipval];
//    
//    
//    // deformation
//    for (i=0;i<ndofn;i++){
//      valnod[i] = new double [nn];
//      if (dd) valnod[i+ndofn] = new double [nn];
//    }
//    
//    for (i=0;i<nn;i++){
//      noddispl (lcid,r,i);
//      for (j=0;j<ndofn;j++){
//        valnod[j][i] = r[j];
//        if (dd) valnod[j+ndofn][i] = r[j];
//      }
//    }
//    
//    // on int. point components
//    for (i=0;i<nipval;i++)
//      valnod[(1+dd)*ndofn+i] = nodvals + i*nn;
//    
//    
//    capi=0;          dimindex[capi]=ndofn;  caption[capi] = new char[4];   sprintf (caption[capi],"def");  capi++;
//    
//    if (dd && ndofn>0) { dimindex[capi]=1;  caption[capi] = new char[5];   sprintf (caption[capi],"defx");  capi++; }
//    if (dd && ndofn>1) { dimindex[capi]=1;  caption[capi] = new char[5];   sprintf (caption[capi],"defy");  capi++; }
//    if (dd && ndofn>2) { dimindex[capi]=1;  caption[capi] = new char[5];   sprintf (caption[capi],"defz");  capi++; }
//    
//    for (i=0;i<nipcomp[0];i++) { dimindex[capi]=1;  caption[capi] = new char[11];   sprintf (caption[capi],"strain[%ld]",i);  capi++;  }
//    for (i=0;i<nipcomp[1];i++) { dimindex[capi]=1;  caption[capi] = new char[11];   sprintf (caption[capi],"stress[%ld]",i);  capi++;  }
//    for (i=0;i<nipcomp[2];i++) { dimindex[capi]=1;  caption[capi] = new char[12];   sprintf (caption[capi],"eqother[%ld]",i); capi++;  }
//    
//    
//    print_dx (gt,file,valnod,NULL,'n',dimindex,caption,capi);
//    
//    
//    delete [] r;
//    delete [] dimindex;
//    delete [] nodvals;
//    for (i=0;i<capi;i++) delete [] caption[i];  delete [] caption;
//    for (i=0;i<(1+dd)*ndofn;i++) delete [] valnod[i];
//    for (i=0;i<nipval;i++) valnod[i+(1+dd)*ndofn] = NULL;
//    delete [] valnod;
//  }
//  ////////////////////       /* termitovo */       ////////////////////////////////////

void aux_mech_nonlin_print (FILE */*aux*/,double */*r*/,double /*l*/)
{
/*  long i,j,k,m,n;
  
  n=Mp->npun;
  
  for (i=0;i<n;i++){
    //  node number
    j=Mp->requn[2*i+0];
    //  unknown number
    k=Mp->requn[2*i+1];
    //  code number
    m=Mt->give_dof (j,k)-1;
    fprintf (aux,"%e  ",r[m]);
  }
  fprintf (aux,"%e\n",l);
  fflush (aux);*/
}

void aux_mech_time_print (FILE */*aux*/,double */*r*/,double /*l*/)
{
/*  long i,j,k,m,n;
  
  n=Mp->npun;
  
  fprintf (aux,"%e",l);
  for (i=0;i<n;i++){
    //  node number
    j=Mp->requn[2*i+0];
    //  unknown number
    k=Mp->requn[2*i+1];
    //  code number
    m=Mt->give_dof (j,k) - 1;
    fprintf (aux,"  %e",r[m]);
  }
  fprintf (aux,"\n");
  fflush (aux);*/
}



/**
  The function prints header, points and used elements to the VTK file

  @param[in,out] out - pointer to the opened VTK file
  @param[in] istep - integration step
  @param[in] time - actual time or load coefficient

  @return The function does not return anything, but changes content of file given by the argument out.

  Created by Vit Smilauer, modified by Tomas Koudelka, 10.2023
*/
void print_vtk_header(FILE *out, long istep, double time)
{
  long i, j, nn;
  // array of nodal correspondence on linear hexahedron element
  long trlh[8] = {6, 7, 4, 5,  //xy, z=0, xi=x,nu=y, theta=z
                  2, 3, 0, 1}; //xy, z<>0
  // array of nodal correspondence on quadratic hexahedron element
  long trqh[20] = { 6,  7,  4,  5,//xy, z=0, xi=x,nu=y, theta=z
                    2,  3,  0,  1,
                   18, 19, 16, 17,
                   10, 11,  8,  9,
                   14, 15, 12, 13};

  fprintf(out, "# vtk DataFile Version 1.0\n");
  fprintf(out, "Time step %ld time %lf\n", istep, time);
  fprintf(out, "ASCII\n\n");
  fprintf(out, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(out, "POINTS %ld float\n", Mt->nn);
  for (i=0; i<Mt->nn; i++){
    fprintf(out, "%e %e %e\n", Gtm->gnodes[i].x, Gtm->gnodes[i].y, Gtm->gnodes[i].z);
  }

  long int ncells = 0;
  for (i=0; i < Mt->ne; i++){
    ncells += Mt->give_nne (i);
  }
  fprintf(out, "\nCELLS %ld %ld\n", Mt->ne, ncells+Mt->ne);
  for (i=0; i < Mt->ne; i++){
    fprintf(out, "%ld ", Mt->give_nne(i));//number of nodes on element
    for(j=0;j<Mt->give_nne(i);j++){
      switch(Mt->elements[i].te){ //due to node renumbering
        case linearhex:
          nn=trlh[j];
          break;
        case quadrhex:
          nn=trqh[j];
          break;
        default:
          nn=j;
          break;
      }
      fprintf(out, "%ld ", Gtm->gelements[i].nodes[nn]);
    }
    fprintf(out, "\n");
  }
  
  //output cell types
  fprintf(out, "\nCELL_TYPES %ld\n", Mt->ne);
  for (i=0; i < Mt->ne; i++){
    switch(Mt->elements[i].te){
      case bar2d:
      case beam2d:
      case bar3d:
      case beam3d:
        fprintf(out, "3\n");//VTK_LINE
        break;
      case barq2d:
      case barq3d:
        fprintf(out, "21\n");//VTK_QUADRATIC_EDGE
        break;
      case planeelementlt:
      case axisymmlt:
        fprintf(out, "5\n");//VTK_TRIANGLE
        break;
      case planeelementlq:
      case axisymmlq:
        fprintf(out, "9\n");//VTK_QUAD
        break;
      case planeelementqq:
      case axisymmqq:
        fprintf(out, "23\n");//VTK_QUADRATIC_QUAD
        break;
      case lineartet:
        fprintf(out, "10\n");//VTK_TETRA
        break;
      case linearhex:
        fprintf(out, "12\n");//VTK_HEXAHEDRON - Renumbered nodes
        break;
      case quadrhex:
        fprintf(out, "25\n");//VTK_QUADRATIC_HEXAHEDRON - Renumbered nodes
        break;
      default:
        print_err("unknown element type in VTK export", __FILE__, __LINE__, __func__);
        abort();
        break;
    }
  }
}



/**
  The function prints data to the VTK file.
  VTK can not handle integration points directly, results in nodes only

  @param[in,out] out - pointer to the opened VTK file
  @param[in] lcid  - load case id

  @return The function does not return anything, but changes content of file given by the argument out.

  Created by Vit Smilauer, modified by Tomas Koudelka, 10.2023
*/
void write_vtk_unkn(FILE *out, long lcid)
{
  vector r;
  long i, j, k, ncomp, max_ncomp;
  fprintf(out, "\nPOINT_DATA %ld\n", Mt->nn);
  //displacements at nodes
  if (Outdm->nog.selndisp.st != sel_no){
    fprintf(out, "VECTORS Nodal_displacements float\n");
    reallocv(RSTCKVEC(3, r));
    for (i=0; i<Mt->nn; i++){
      r[0]=r[1]=r[2]=0;
      noddispl (lcid, r.a, i);
      fprintf(out, "%e %e %e\n", r[0], r[1], r[2]);
    }
    fprintf(out, "\n");
  }
    // strains at nodes
  if (Outdm->nog.selnstra.st != sel_no){
    fprintf(out, "TENSORS Strain float\n");
    for (i=0; i<Mt->nn; i++){
      ncomp = Mt->nodes[i].ncompstr;
      vector eps(ASTCKVEC(4));
      matrix epst(ASTCKMAT(3,3));
      nullv(eps);
      nullm(epst);
      for (j=0; j<ncomp; j++)
        eps[j] = Mt->nodes[i].strain[j];
      switch (ncomp){
        case 1:{
          vector_tensor (eps, epst, bar, strain);
          break;
        }
        case 4:{
          vector_tensor (eps, epst, planestrain, strain);//or planestress
          break;
        }
        case 6:{
          vector_tensor (eps, epst, spacestress, strain);
        }
      }
      for (k=0; k<3; k++){
        fprintf(out, "%e %e %e\n", epst[k][0], epst[k][1], epst[k][2]);
      }
      fprintf(out, "\n");
    }
  }
  // stresses at nodes
  if (Outdm->nog.selnstre.st != sel_no){
    fprintf(out, "TENSORS Stress float\n");
    for (i=0; i<Mt->nn; i++){
      ncomp = Mt->nodes[i].ncompstr;
      vector sig(ASTCKVEC(4));
      matrix sigt(ASTCKMAT(3,3));
      nullv(sig);
      nullm(sigt);
      for (j=0; j<ncomp; j++)
        sig[j] = Mt->nodes[i].stress[j];
      switch (ncomp){
        case 1:{
          vector_tensor (sig, sigt, bar, stress);
          break;
        }
        case 4:{
          vector_tensor (sig, sigt, planestrain, stress);//or planestress
          break;
        }
        case 6:{
          vector_tensor (sig, sigt, spacestress, stress);
        }
      }
      for (k=0; k<3; k++){
        fprintf(out, "%e %e %e\n", sigt[k][0], sigt[k][1], sigt[k][2]);
      }
      fprintf(out, "\n");
    }
  }
  // other values at nodes
  if (Outdm->nog.selnoth.st != sel_no){
    max_ncomp = 0;
    for (k=0; k<Mm->tnip; k++){//find maximum eqother components
      max_ncomp > Mm->ip[k].ncompother ? max_ncomp = Mm->ip[k].ncompother : 0;
    }
    
    for(i=0; i<max_ncomp; i++){// other[i]
      for (j=0; j<Outdm->nog.selnoth.n; j++){
        if (Outdm->nog.seloth[j].presence_id(i)){
          fprintf(out, "\nSCALARS Eq_other_s%ld_%ld double\nLOOKUP_TABLE default\n", j+1, i+1);
          for (k=0; k < Mt->nn; k++)
            fprintf(out, "%e\n", Mt->nodes[k].giveother(i));
        }
      }
    }
/*    fprintf(out, "TENSORS Other float\n");
    for (i=0; i<Mt->nn; i++){
      ncomp = Mt->nodes[i].ncompother;
      vector oth(ASTCKVEC(4));
      matrix otht(ASTCKMAT(3,3));
      nullv(oth);
      nullm(otht);
      for (j=0; j<ncomp; j++)
        oth[j] = Mt->nodes[i].other[j];
      switch (ncomp){
        case 1:{
          vector_tensor (oth, otht, bar, stress);
          break;
        }
        case 4:{
          vector_tensor (oth, otht, planestrain, stress);//or planestress
          break;
        }
        case 6:{
          vector_tensor (oth, otht, spacestress, stress);
        }
      }
      for (k=0; k<3; k++){
        fprintf(out, "%e %e %e\n", otht[k][0], otht[k][1], otht[k][2]);
      }
      fprintf(out, "\n");
    }*/
  }
  // reactions - not implemented
}
