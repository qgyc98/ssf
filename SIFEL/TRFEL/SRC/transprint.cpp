#include "transprint.h"
#include "globalt.h"
#include "globmatt.h"
#include "aliast.h"
#include "mathem.h"
#include "onemedium.h"
#include "twomedia.h"
#include "threemedia.h"


/**
  Function opens data files for each type of outputfile

  @param istep - step id. The parameter enables to open file with name enhanced
                 by the step id - istep >= 0 or it leaves required filename untouched for instance that
                 istep < 0. In case of stochastic calculations and in case istep >= 0 the istep precedes the 
                 stochastic step id.
  @param mode - string with control sequence for the file opening. It enables to open 
                new file (mode = "wt") or to append existing ones (mode = "at").
  @param idn1 - id of the first node for GiD mesh (default is idn1 = 1 -> nodes are numbered from 1)
  @param ide1 - id of the first element for GiD mesh (default is ide1 = 1 -> elements are numbered from 1)
  Return :
  The function does not return anything
*/
void print_initt(long istep, const char *mode, long idn1, long ide1)
{
  char fname[FNAMELEN+30];
  char *path;
  char *name;
  char *suffix;
  long i;
  long sl;
  
  
  Outdt->idn1 = idn1;
  Outdt->ide1 = ide1;
  if ((Outdt->outf == NULL) && (Outdt->textout==1))
  {
    filename_decomposition(Outdt->outfn,path,name,suffix);
    if (Stt == NULL)
    {
      if (istep < 0)
        sprintf(fname, "%s%s%s", path, name, suffix);
      else
        sprintf(fname, "%s%s.%ld%s", path, name, istep+1, suffix);
    }
    else
    {
      if (istep < 0)
        sprintf(fname, "%s%s.%ld%s", path, name, Tp->ns+1, suffix);
      else
        sprintf(fname, "%s%s.%ld.%ld%s", path, name, istep+1, Tp->ns+1, suffix);
    }
    Outdt->outf = fopen(fname, mode);
    if (Outdt->outf == NULL)
    {
      fprintf(stderr, "\n\nUnable to open output text file '%s' in function print_init()\n", fname);
      fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
      abort();
    }
    fseek(Outdt->outf, 0, SEEK_END); // MS Visual C++ requires that
    if (ftell(Outdt->outf) == 0)
      Outdt->print_header(Outdt->outf);
      
    delete [] path; delete [] name; delete [] suffix;
  }
  
  if ((Outdt->gf != grftt_no) && (Outdt->outgr == NULL))
  {
    filename_decomposition(Outdt->outgrfn,path,name,suffix);
    if (Stt == NULL)
    {
      if (istep < 0)
        sprintf(fname, "%s%s%s", path, name, suffix);
      else
        sprintf(fname, "%s%s.%ld%s", path, name, istep+1, suffix);
    }
    else
    {
      if (istep < 0)
        sprintf(fname, "%s%s.%ld%s", path, name, Tp->ns+1, suffix);
      else
        sprintf(fname, "%s%s.%ld.%ld%s", path, name, istep+1, Tp->ns+1, suffix);
    }
    if ((Outdt->gf == grftt_gid) || (Outdt->gf == grftt_gid_sep) || (Outdt->gf == grftt_gid_vtk) || (Outdt->gf == grftt_gidsep_vtk))
    {
      sl = long(strlen(fname));
      sprintf(fname+sl, ".msh");
    }
    if(Outdt->gf != grftt_vtk)
      Outdt->outgr = fopen(fname, mode);
    
    if (Outdt->outgr == NULL && Outdt->gf != grftt_vtk)
    {
      fprintf(stderr, "\n\nUnable to open graphics file '%s' in function print_init()\n", fname);
      fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
      abort();
    }
    switch(Outdt->gf)
    {
      case grftt_no:
        break;   
      case grftt_open_dx:
        break;
      case grftt_femcad:
        //export_femcad(Outdt->outgr);
        break;
      case grftt_gid:
      case grftt_gid_vtk:
        fseek(Outdt->outgr, 0, SEEK_END); // MS Visual C++ requires that
        if ((ftell(Outdt->outgr) == 0) || Adat)
          export_gid_mesht(Outdt->outgr, idn1, ide1);
        fclose(Outdt->outgr);
        if (Outdt->ncut > 0)
        {
          sprintf(fname+sl, "2d.msh");
          Outdt->outgr = fopen(fname, mode);
          fseek(Outdt->outgr, 0, SEEK_END); // MS Visual C++ requires that
          if ((ftell(Outdt->outgr) == 0) || Adat)
            export_gid_2dmesht(Outdt->outgr, Outdt->ncut, idn1, ide1);
          fclose(Outdt->outgr);
        }
        sprintf(fname+sl, ".res");
        Outdt->outgr = fopen(fname, mode);
        if (Outdt->outgr == NULL)
	{
          fprintf(stderr, "\n\nUnable to open graphics file '%s' in function print_init()\n", fname);
          fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
          abort();
        }
        fseek(Outdt->outgr, 0, SEEK_END); // MS Visual C++ requires that
        if (ftell(Outdt->outgr) == 0)
        {
          fprintf(Outdt->outgr, "GiD Post Results File 1.0\n");
          export_gid_gauss_ptt(Outdt->outgr);
        }
        if (Adat)
          fprintf(Outdt->outgr, "\nOnGroup adapt_step_%ld\n\n", Adat->istep+1);
        break;
      case grftt_gid_sep:
      case grftt_gidsep_vtk:
        fseek(Outdt->outgr, 0, SEEK_END); // MS Visual C++ requires that
        if ((ftell(Outdt->outgr) == 0) || Adat)
          export_gid_mesht(Outdt->outgr, idn1, ide1);
        fclose(Outdt->outgr);
        if (Outdt->ncut > 0)
	{
          sprintf(fname+sl, "2d.msh");
          Outdt->outgr = fopen(fname, mode);
          fseek(Outdt->outgr, 0, SEEK_END); // MS Visual C++ requires that
          if ((ftell(Outdt->outgr) == 0) || Adat)
            export_gid_2dmesht(Outdt->outgr, Outdt->ncut, idn1, ide1);
          fclose(Outdt->outgr);
        }
        strncpy(Outdt->outgrfngs, fname, sl);
        Outdt->outgrfngs[sl]=0;
        Outdt->outgr = NULL;
        Outdt->create_files_gidsp(mode);
        break;
      case grftt_vtk:
        break;
      default:
        fprintf(stderr, "\n\nUnknown type of graphics format is required in function print_init\n");
        fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
    }
    delete [] path; delete [] name; delete [] suffix;
  }
  if ((Outdt->ndiag > 0) && (Outdt->outdiagf[0] == NULL))
  {
    filename_decomposition(Outdt->outdiagfn,path,name,suffix);
    for (i=0; i<Outdt->ndiag; i++)
    {
      if (Stt == NULL)
      {
        if (Outdt->ndiag > 1)
          sprintf(fname, "%s%s.%ld%s", path, name, i+1, suffix);
        else
          sprintf(fname, "%s%s%s", path, name, suffix);
      }
      else
      {
        if (Outdt->ndiag > 1)
          sprintf(fname, "%s%s.%ld.%ld%s", path, name, i+1, Tp->ns+1, suffix);
        else
          sprintf(fname, "%s%s.%ld%s", path, name, Tp->ns+1, suffix);
      }
      if (Outdt->outdiagf[i] == NULL)
      {
        Outdt->outdiagf[i] = fopen(fname, mode);               
        if (Outdt->outdiagf[i] == NULL)
        {
          fprintf(stderr, "\n\nUnable to open diagram file '%s' in function print_init()\n", fname);
          fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
          abort();
        }
        Outdt->odiag[i].print_header(Outdt->outdiagf[i]);
      }
    }
    delete [] path; delete [] name; delete [] suffix;
  }
  //    Outdt->idn1 = idn1;
  //    Outdt->ide1 = ide1;
}



/**
  The function performs all prints requirements for given step.

  @param lcid - load case id
  @param istep - step id
  @param lambda - load coefficient or time (it depends on problem type and solver)
  @param fi - array of additional values at nodes which should be printed.
              It depends on problem and solver type. For example for nonlinear
              statics it contains the internal forces.
  
  Return :
  The function does not return anything.

  created 2004 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz

*/
void print_stept(long lcid, long istep, double lambda, double *fi)
{

  if(Tp->homogt == 0)//in case of no homogenization
    compute_req_valt (lcid);
  switch(Tp->tprob)
    {
    case stationary_problem:
    case nonlinear_stationary_problem:
      Outdt->print_out(Outdt->outf, lcid, istep, lambda);
      Outdt->print_graphics(Outdt->outgr, lcid, lambda, istep, fi);
      break;
    case nonstationary_problem:
    case nonlinear_nonstationary_problem:
    case discont_nonstat_problem:
    case discont_nonlin_nonstat_problem:
    case growing_np_problem:
    case growing_np_problem_nonlin:
      Outdt->print_newstep(Outdt->outf, lcid, istep, lambda);
      Outdt->print_out(Outdt->outf, lcid, istep, lambda);
      Outdt->print_diags(lcid, lambda, istep, fi);
      Outdt->print_graphics(Outdt->outgr, lcid, lambda, istep, fi);
      break;
    default:
      fprintf(stderr, "\n\nUnsupported problem type is required in function print_step\n");
      fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
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

  Created 2014 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz

*/
void print_stept_forced(long lcid, long istep, double lambda, double *fi)
{
  
  compute_req_valt (lcid);
  switch(Tp->tprob)
    {
    case stationary_problem:
    case nonlinear_stationary_problem:
      Outdt->print_out_forced(Outdt->outf, lcid, istep, lambda);
      Outdt->print_graphics_forced(Outdt->outgr, lcid, lambda, istep, fi);
      break;
    case nonstationary_problem:
    case nonlinear_nonstationary_problem:
    case discont_nonstat_problem:
    case discont_nonlin_nonstat_problem:
    case growing_np_problem:
    case growing_np_problem_nonlin:
      Outdt->print_newstep(Outdt->outf, lcid, istep, lambda);
      Outdt->print_out_forced(Outdt->outf, lcid, istep, lambda);
      Outdt->print_diags_forced(lcid, lambda, istep, fi);
      Outdt->print_graphics_forced(Outdt->outgr, lcid, lambda, istep, fi);
      break;
    default:
      fprintf(stderr, "\n\nUnsupported problem type is required in function print_step\n");
      fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
      break;
  }
}



void print_flusht()
{
  long i;
  if (Outdt->outf)
    fflush(Outdt->outf);
  if (Outdt->outgr)
    fflush(Outdt->outgr);
  for (i=0; i<Outdt->ndiag; i++)
  {
    if (Outdt->outdiagf[i])
      fflush(Outdt->outdiagf[i]);
  }
}



void print_closet()
{
  long i;
  if (Outdt->outf)
    fclose(Outdt->outf);
  Outdt->outf = NULL;
  if (Outdt->outgr)
  {
    if (Adat)
      fprintf(Outdt->outgr, "\nEnd OnGroup\n\n");
    fclose(Outdt->outgr);
  }
  Outdt->outgr = NULL;
  if ((Outdt->gf == grftt_gid_sep) && Adat)
    Outdt->close_files_gidsp();

  for (i=0; i<Outdt->ndiag; i++)
  {
    if (Outdt->outdiagf[i])
      fclose(Outdt->outdiagf[i]);
    Outdt->outdiagf[i] = NULL;
  }
}



/**
  The function exports sets of used elements to the file given by parameter out in GiD format. 

  Parameters :
  @param out  - pointer to the opened text file where the output will be produced
  @param idn1 - id of the first node for GiD mesh (default should be idn1 = 1 -> nodes are numbered from 1)
  @param ide1 - id of the first element for GiD mesh (default should be ide1 = 1 -> elements are numbered from 1)

  Return :
  The function does not return anything
*/void export_gid_mesht(FILE *out, long idn1, long ide1)
{
  long i, print_header, print_coord = 1;
  
  if (Adat)
    fprintf(out, "Group adapt_step_%ld\n", Adat->istep+1);

  print_header = 1;
  for (i=0; i < Tt->ne; i++)
    {
      switch(Tt->elements[i].te)
	{
	case barlintax:
	case barlint:
	case barlint3d:
	  if (print_header)
	    {
	      fprintf(out, "MESH beams2 dimension 3  Elemtype Linear Nnode 2\n");
	      //print_header = 0;
	    }
	  if (print_coord)
	    {
	      write_gid_nodest(out, idn1);
	      print_coord = 0;
	    }
	  if (print_header == 1)
	    fprintf(out, "Elements\n");
	  write_gid_elementt(out, i, idn1, ide1);
	  print_header = 0;
	  break;
	default:
	  break;
	}
    }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  
  
  print_header = 1;
  for (i=0; i < Tt->ne; i++)
    {
      switch(Tt->elements[i].te)
	{
	case barquadt:
	case barquadtax:
	  if (print_header)
	    {
	      fprintf(out, "MESH beams3 dimension 3  Elemtype Linear Nnode 3\n");
	      //print_header = 0;
	    }
	  if (print_coord)
	    {
	      write_gid_nodest(out, idn1);
	      print_coord = 0;
	    }
	  if (print_header == 1)
	    fprintf(out, "Elements\n");
	  write_gid_elementt(out, i, idn1, ide1);
	  print_header = 0;
	  break;    
	default:
	  break;
	}
    }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  
  
  print_header = 1;
  for (i=0; i < Tt->ne; i++)
    {
      switch(Tt->elements[i].te)
	{
	case trlint:
	case trlaxisym:
	  if (print_header)
	    {
	      fprintf(out, "MESH trias3 dimension 3  Elemtype Triangle Nnode 3\n");
	      //print_header = 0;
	    }
	  if (print_coord)
	    {
	      write_gid_nodest(out, idn1);
	      print_coord = 0;
	    }
	  if (print_header == 1)
	    fprintf(out, "Elements\n");
	  write_gid_elementt(out, i, idn1, ide1);
	  print_header = 0;
	  break;    
	default:
	  break;
	}
    }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  

  print_header = 1;
  for (i=0; i < Tt->ne; i++)
    {
      switch(Tt->elements[i].te)
	{
	case trquadt://prepared for quadratic elements
	case trqaxisym:
	  if (print_header)
	    {
	      fprintf(out, "MESH trias6 dimension 3  Elemtype Triangle Nnode 6\n");
	      //print_header = 0;
	    }
	  if (print_coord)
	    {
	      write_gid_nodest(out, idn1);
	      print_coord = 0;
	    }
	  if (print_header == 1)
	    fprintf(out, "Elements\n");
	  write_gid_elementt(out, i, idn1, ide1);
	  print_header = 0;
	  break;    
	default:
	  break;
	}
    }
  if (print_header == 0)
    fprintf(out, "end Elements\n");
  
  
  print_header = 1;
  for (i=0; i < Tt->ne; i++)
    {
      switch(Tt->elements[i].te)
	{
	case quadlint:
	case quadlaxisym:
	  if (print_header)
	    {
	      fprintf(out, "MESH Quads4 dimension 3  Elemtype Quadrilateral Nnode 4\n");
	      //print_header = 0;
	    }
	  if (print_coord)
	    {
	      write_gid_nodest(out, idn1);
	      print_coord = 0;
	    }
	  if (print_header == 1)
	    fprintf(out, "Elements\n");
	  write_gid_elementt(out, i, idn1, ide1);
	  print_header = 0;
	  break;    
	default:
	  break;
	}
    }
  if (print_header == 0)
    fprintf(out, "end Elements\n");

  
  print_header = 1;
  for (i=0; i < Tt->ne; i++)
    {
      switch(Tt->elements[i].te)
	{
	case quadquadtax:
	case quadquadt:
	  if (print_header)
	    {
	      fprintf(out, "MESH Quads8 dimension 3  Elemtype Quadrilateral Nnode 8\n");
	      //print_header = 0;
	    }
	  if (print_coord)
	    {
	      write_gid_nodest(out, idn1);
	      print_coord = 0;
	    }
	  if (print_header == 1)
	    fprintf(out, "Elements\n");
	  write_gid_elementt(out, i, idn1, ide1);
	  print_header = 0;
	  break;    
	default:
	  break;
	}
    }
  if (print_header == 0)
    fprintf(out, "end Elements\n");

  
  print_header = 1;
  for (i=0; i < Tt->ne; i++)
    {
      switch(Tt->elements[i].te)
	{
	case lineartett:
	  if (print_header)
	    {
	      fprintf(out, "MESH Tetras4 dimension 3  Elemtype Tetrahedra Nnode 4\n");
	      //print_header = 0;
	    }
	  if (print_coord)
	    {
	      write_gid_nodest(out, idn1);
	      print_coord = 0;
	    }        
	  if (print_header == 1)
	    fprintf(out, "Elements\n");
	  write_gid_elementt(out, i, idn1, ide1);
	  print_header = 0;
	  break;    
	default:
	  break;
	}
    }
  if (print_header == 0)
    fprintf(out, "end Elements\n");

  
  print_header = 1;
  for (i=0; i < Tt->ne; i++)
    {
      switch(Tt->elements[i].te)
	{
	case linearhext:
	  if (print_header)
	    {
	      fprintf(out, "MESH Brick8 dimension 3  Elemtype Hexahedra Nnode 8\n");
	      //print_header = 0;
	    }
	  if (print_coord)
	    {
	      write_gid_nodest(out, idn1);
	      print_coord = 0;
	    }
	  if (print_header == 1)
	    fprintf(out, "Elements\n");
	  write_gid_elementt(out, i, idn1, ide1);
	  print_header = 0;
	  break;    
	default:
	  break;
	}
    }
  if (print_header == 0)
    fprintf(out, "end Elements\n");

  
  print_header = 1;
  for (i=0; i < Tt->ne; i++)
    {
      switch(Tt->elements[i].te)
	{
	case quadratichext:
	  if (print_header)
	    {
	      fprintf(out, "MESH Brick20 dimension 3  Elemtype Hexahedra Nnode 20\n");
	      //print_header = 0;
	    }
	  if (print_coord)
	    {
	      write_gid_nodest(out, idn1);
	      print_coord = 0;
	    }
	  if (print_header == 1)	  
	    fprintf(out, "Elements\n");
	  write_gid_elementt(out, i, idn1, ide1);
	  print_header = 0;
	  break;    
	default:
	  break;
	}
    }
  if (print_header == 0)
    fprintf(out, "end Elements\n");

  print_header = 1;
  for (i=0; i < Tt->ne; i++)
  {
    if (Gtt->leso[i]==0)
      continue;
    
    switch(Tt->elements[i].te)
    {
      case linearwedget:
        if (print_header)
	{
          fprintf(out, "MESH \"Wedge6\" dimension 3  Elemtype Prism Nnode 6\n");
          print_header = 0;
          if (print_coord)
	  {
            write_gid_nodest(out, idn1);
            print_coord = 0;
	  }
          fprintf(out, "Elements\n");
	}
        write_gid_elementt(out, i, idn1, ide1);
        break;    
      default:
        break;
    }
  }
  if (print_header == 0)
  fprintf(out, "end Elements\n");
  if (Adat)
    fprintf(out, "End Group\n\n\n");
}


/**
  The function exports from 3D brick element mesh the 2D plane element cuts to the file given by parameter out in GiD format. 

  Parameters :
  @param out  - pointer to the opened text file where the output will be produced
  @param icut - index of the cut - will be used in the element property
  @param idn1 - id of the first node for GiD mesh (default should be idn1 = 1 -> nodes are numbered from 1)
  @param ide1 - id of the first element for GiD mesh (default should be ide1 = 1 -> elements are numbered from 1)

  Return :
  The function does not return anything
*/
void export_gid_2dmesht(FILE *out, long icut, long idn1, long ide1)
{
  long i, print_header, print_coord = 1;
  long range = Tt->ne/icut;
  
  if (Adat)
    fprintf(out, "Group adapt_step_%ld\n", Adat->istep+1);

  print_header = 1;
  for (i=0; i < Tt->ne; i++)
    {
      switch(Tt->elements[i].te)
	{
	case linearhext:
	  if (print_header)
	    {
	      fprintf(out, "MESH Quads4 dimension 3  Elemtype Quadrilateral Nnode 4\n");
	      //print_header = 0;
	    }
	  if (print_coord)
	    {
	      write_gid_nodest(out, idn1);
	      print_coord = 0;
	    }
	  if (print_header == 1)	    
	    fprintf(out, "Elements\n");
	  write_gid_2delementt(out, i, 0, 4, i/range, 0, idn1, ide1);
	  print_header = 0;
	  break;    
	default:
	  break;
	}
    }
  for (i=Tt->ne-range; i < Tt->ne; i++)
    {
      write_gid_2delementt(out, i, 4, 8, icut, range, idn1, ide1);
    }
  if (print_header == 0)
    fprintf(out, "end Elements\n");

  if (Adat)
    fprintf(out, "End Group\n\n\n");
}


/**
  The function exports sets of gausspoints of used elements to the file given by parameter 
  out in GiD format. 

  Parameters :
  @param out  - pointer to the opened text file where the output will be produced

  Return :
  The function does not return anything
*/
void export_gid_gauss_ptt(FILE *out)
{
  long i, j, k, ii, jj, brk,ngp;
  vector gp1, gp2, gp3, gp, w, wt;

  ii = 0;//only for the first block
  jj = 0;//only for the first block

  // linear bar element //
  brk = 0;
  for (i=0; i < Tt->ne; i++)
  {
    if (Gtt->leso[i]==0)
      continue;

    switch(Tt->elements[i].te)
    {
      case barlint:
      case barlint3d:
      case barlintax :
        fprintf(out, "GaussPoints \"Lin_1D\" Elemtype Linear\n");
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

  // quadratic bar element //
  brk = 0;
  for (i=0; i < Tt->ne; i++)
  {
    if (Gtt->leso[i]==0)
      continue;

    switch(Tt->elements[i].te)
    {
      case barquadt :
      case barquadtax :
        fprintf(out, "GaussPoints \"Quad_1D\" Elemtype Linear\n");
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

  // linear triangular element //
  for (i=0; i < Tt->ne; i++)
  {
    if (Gtt->leso[i]==0)
      continue;

    if (Tt->elements[i].te == trlint)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Triangle\n", trlint);
      ngp=1;
      fprintf(out, "Number Of Gauss Points: %ld\n",ngp);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");

      if (Ltt->intordkm[ii][jj]!=0){//only for the first block
	reallocv(RSTCKVEC(Ltt->intordkm[ii][jj],w));
	reallocv(RSTCKVEC(Ltt->intordkm[ii][jj],gp1));
	reallocv(RSTCKVEC(Ltt->intordkm[ii][jj],gp2));
	gauss_points_tr (gp1.a,gp2.a,w.a,Ltt->intordkm[ii][jj]);
	for (i=0;i<Ltt->intordkm[ii][jj];i++)
	  fprintf(out, "%le %le\n", gp1[i], gp2[i]);
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }


  // linear triangular element - axisymmetric//
  for (i=0; i < Tt->ne; i++)
  {
    if (Gtt->leso[i]==0)
      continue;

    if (Tt->elements[i].te == trlaxisym)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Triangle\n", trlaxisym);
      ngp=1;
      fprintf(out, "Number Of Gauss Points: %ld\n",ngp);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");

      if (Ltat->intordkm[ii][jj]!=0){//only for the first block 
	reallocv(RSTCKVEC(Ltat->intordkm[ii][jj],w));
	reallocv(RSTCKVEC(Ltat->intordkm[ii][jj],gp1));
	reallocv(RSTCKVEC(Ltat->intordkm[ii][jj],gp2));
	gauss_points_tr (gp1.a,gp2.a,w.a,Ltat->intordkm[ii][jj]);
	for (i=0;i<Ltat->intordkm[ii][jj];i++)
	  fprintf(out, "%le %le\n", gp1[i], gp2[i]);
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }


  // linear quadrilateral element //
  for (i=0; i < Tt->ne; i++)
  {
    if (Gtt->leso[i]==0)
      continue;

    if (Tt->elements[i].te == quadlint)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Quadrilateral\n", quadlint);
      ngp=4;
      fprintf(out, "Number Of Gauss Points: %ld\n",ngp);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      
      if (Lqt->intordkm[ii][jj]!=0) {//only for the first block  
	reallocv(RSTCKVEC(Lqt->intordkm[ii][jj],w));
	reallocv(RSTCKVEC(Lqt->intordkm[ii][jj],gp));
	gauss_points (gp.a,w.a,Lqt->intordkm[ii][jj]);
	for (i=0;i<Lqt->intordkm[ii][jj];i++)
	  for (j=0; j<Lqt->intordkm[ii][jj]; j++)
	    fprintf(out, "%le %le\n", -gp[i], -gp[j]);
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }


  // quadratic quadrilateral element //
  for (i=0; i < Tt->ne; i++)
  {
    if (Gtt->leso[i]==0)
      continue;

    if (Tt->elements[i].te == quadquadt)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Quadrilateral\n", quadquadt);
      ngp=9;
      fprintf(out, "Number Of Gauss Points: %ld\n",ngp);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      
      if (Qqt->intordkm[ii][jj]!=0){//only for the first block   
	reallocv(RSTCKVEC(Qqt->intordkm[ii][jj],w));
	reallocv(RSTCKVEC(Qqt->intordkm[ii][jj],gp));
	gauss_points (gp.a,w.a,Qqt->intordkm[ii][jj]);
	for (i=0;i<Qqt->intordkm[ii][jj];i++)
	  for (j=0; j<Qqt->intordkm[ii][jj]; j++)
	    fprintf(out, "%le %le\n", -gp[i], -gp[j]);
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }


  // linear quadrilateral element - axisymmetric//
  for (i=0; i < Tt->ne; i++)
  {
    if (Gtt->leso[i]==0)
      continue;

    if (Tt->elements[i].te == quadlaxisym)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Quadrilateral\n", quadlaxisym);
      ngp=4;
      fprintf(out, "Number Of Gauss Points: %ld\n",ngp);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");

      if (Lqat->intordkm[ii][jj]!=0){//only for the first block   
	reallocv(RSTCKVEC(Lqat->intordkm[ii][jj],w));
	reallocv(RSTCKVEC(Lqat->intordkm[ii][jj],gp));
	gauss_points (gp.a,w.a,Lqat->intordkm[ii][jj]);
	for (i=0;i<Lqat->intordkm[ii][jj];i++)
	  for (j=0; j<Lqat->intordkm[ii][jj]; j++)
	    fprintf(out, "%le %le\n", -gp[i], -gp[j]);
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }


  // quadratic quadrilateral element - axisymmetric//
  for (i=0; i < Tt->ne; i++)
  {
    if (Gtt->leso[i]==0)
      continue;

    if (Tt->elements[i].te == quadquadtax)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Quadrilateral\n", quadquadtax);
      ngp=9;
      fprintf(out, "Number Of Gauss Points: %ld\n",ngp);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      
      if (Qqat->intordkm[ii][jj]!=0) {//only for the first block  
	reallocv(RSTCKVEC(Qqat->intordkm[ii][jj],w));
	reallocv(RSTCKVEC(Qqat->intordkm[ii][jj],gp));
	gauss_points (gp.a,w.a,Qqat->intordkm[ii][jj]);
	for (i=0;i<Qqat->intordkm[ii][jj];i++)
	  for (j=0; j<Qqat->intordkm[ii][jj]; j++)
	    fprintf(out, "%le %le\n", -gp[i], -gp[j]);
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  // linear tetrahedral element //
  for (i=0; i < Tt->ne; i++)
  {
    if (Gtt->leso[i]==0)
      continue;

    if (Tt->elements[i].te == lineartett)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype TetraHedra\n", lineartett);
      ngp=1;
      fprintf(out, "Number Of Gauss Points: %ld\n",ngp);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");

      if (Ltett->intordkm[ii][jj]!=0)  {//only for the first block 
	reallocv (RSTCKVEC(Ltett->intordkm[ii][jj],w));
        reallocv (RSTCKVEC(Ltett->intordkm[ii][jj],gp1));
        reallocv (RSTCKVEC(Ltett->intordkm[ii][jj],gp2));
        reallocv (RSTCKVEC(Ltett->intordkm[ii][jj],gp3));
	gauss_points_tet (gp1.a,gp2.a,gp3.a,w.a,Ltett->intordkm[ii][jj]);
	for (i=0;i<Ltett->intordkm[ii][jj];i++)
	  fprintf(out, "%le %le %le\n", gp1[i], gp2[i], gp3[i]);
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  // linear hexahedral element //
  for (i=0; i < Tt->ne; i++)
  {
    if (Gtt->leso[i]==0)
      continue;

    if (Tt->elements[i].te == linearhext)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Hexahedra\n", linearhext);
      ngp=8;
      fprintf(out, "Number Of Gauss Points: %ld\n",ngp);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");

      if (Lht->intordkm[ii][jj]!=0)  {//only for the first block 
	reallocv (RSTCKVEC(Lht->intordkm[ii][jj],w));
	reallocv (RSTCKVEC(Lht->intordkm[ii][jj],gp));
	gauss_points(gp.a,w.a,Lht->intordkm[ii][jj]);
	for (i=0;i<Lht->intordkm[ii][jj];i++){
          for (j=0;j<Lht->intordkm[ii][jj];j++){
            for (k=0;k<Lht->intordkm[ii][jj];k++)
              fprintf(out, "%le %le %le\n", -gp[i], -gp[j], -gp[k]);
          }
        }
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  // quadratic hexahedral element //
  for (i=0; i < Tt->ne; i++)
  {
    if (Gtt->leso[i]==0)
      continue;

    if (Tt->elements[i].te == quadratichext)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Hexahedra\n", quadratichext);
      ngp=27;
      fprintf(out, "Number Of Gauss Points: %ld\n",ngp);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");
      
      if (Qht->intordkm[0][0]!=0){  //only for the first block
	reallocv (RSTCKVEC(Qht->intordkm[0][0],w));
	reallocv (RSTCKVEC(Qht->intordkm[0][0],gp));
	gauss_points(gp.a,w.a,Qht->intordkm[0][0]);
	for (i=0;i<Qht->intordkm[0][0];i++){
          for (j=0;j<Qht->intordkm[0][0];j++){
            for (k=0;k<Qht->intordkm[0][0];k++)
              fprintf(out, "%le %le %le\n", -gp[i], -gp[j], -gp[k]);
          }
        }
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }
  // linear wedge element //
  for (i=0; i < Tt->ne; i++)
  {
    if (Gtt->leso[i]==0)
      continue;

    if (Tt->elements[i].te == linearwedget)
    {
      fprintf(out, "GaussPoints \"%d\" Elemtype Prism\n", linearwedget);
      ngp=6;
      fprintf(out, "Number Of Gauss Points: %ld\n",ngp);
      fprintf(out, "Nodes not included\n");
      fprintf(out, "Natural coordinates: given\n");

      if (Lwt->intordkmt[ii][jj]!=0)  {//only for the first block 
        reallocv (RSTCKVEC(Lwt->intordkmz[ii][jj],w));
        reallocv (RSTCKVEC(Lwt->intordkmz[ii][jj],gp));
        reallocv (RSTCKVEC(Lwt->intordkmt[ii][jj],wt));
        reallocv (RSTCKVEC(Lwt->intordkmt[ii][jj],gp1));
        reallocv (RSTCKVEC(Lwt->intordkmt[ii][jj],gp2));
          
        gauss_points (gp.a,w.a,Lwt->intordkmz[ii][jj]);
        gauss_points_tr (gp1.a,gp2.a,wt.a,Lwt->intordkmt[ii][jj]);
        
        for (i=0;i<Lwt->intordkmt[ii][jj];i++){
          for (k=0;k<Lwt->intordkmz[ii][jj];k++){
            fprintf(out, "%le %le %le\n", gp2[i], gp1[(i+1)%3], 0.5*(1.0-gp[k]));
          }
        }
      }
      fprintf(out, "end GaussPoints\n\n");
      break;
    }
  }

  
}




/**
  The function exports nodes to the file given by parameter out in GiD format. 

  Parameters :
  @param out  - pointer to the opened text file where the output will be produced
  @param idn1 - id of the first node for GiD mesh (default should be idn1 = 1 -> nodes are numbered from 1)

  Return :
  The function does not return anything
*/
void write_gid_nodest(FILE *out, long idn1)
{
  long i;
  fprintf(out, "Coordinates\n");
  for (i=0; i<Tt->nn; i++)
    fprintf(out, "%ld %e %e %e\n", i+idn1, Gtt->gnodes[i].x, Gtt->gnodes[i].y, Gtt->gnodes[i].z);
  fprintf(out, "end Coordinates\n");
}

/**
  The function exports element to the file given by parameter out in GiD format. 

  Parameters :
  @param out  - pointer to the opened text file where the output will be produced
  @param i    - element id
  @param idn1 - id of the first node for GiD mesh (default should be idn1 = 1 -> nodes are numbered from 1)
  @param ide1 - id of the first element for GiD mesh (default should be ide1 = 1 -> elements are numbered from 1)

  Return :
  The function does not return anything
*/
void write_gid_elementt(FILE *out, long i, long idn1, long ide1)
{
  long j;
  fprintf(out, "%ld ", i+ide1);
  for (j=0; j<Tt->give_nne(i); j++)
    fprintf(out, "%ld ", Gtt->gelements[i].nodes[j]+idn1);
  fprintf(out, "%ld\n", Tt->elements[i].idm[0]+1);
}

/**
  The function exports 3D brick element as 2D plane element to the file given by parameter out in GiD format. 
  Purpose of the function is to perform cuts of the domain which consists of the regular mesh
  of brick elements.

  Parameters :
  @param out  - pointer to the opened text file where the output will be produced
  @param i    - brick element id
  @param id1  - brick element index of the first node for plane element
  @param nne  - brick element index of the last node for plane element
  @param icut - index of the cut - will be used in the element property
  @param di   - the start index of the plane elements
  @param idn1 - id of the first node for GiD mesh (default should be idn1 = 1 -> nodes are numbered from 1)
  @param ide1 - id of the first element for GiD mesh (default should be ide1 = 1 -> elements are numbered from 1)

  Return :
  The function does not return anything
*/
void write_gid_2delementt(FILE *out, long i, long id1, long nne, long icut, long di, long idn1, long ide1)
{
  long j;
  fprintf(out, "%ld ", i+ide1+di);
  for (j=id1; j<nne; j++)
    fprintf(out, "%ld ", Gtt->gelements[i].nodes[j]+idn1);
  fprintf(out, "%ld%ld\n", icut+1, Tt->elements[i].idm[0]+1);
}

void write_gid_unkn(FILE *out, long lcid, const char *desclcid)
{
  long i, j;
  double r;

  for (i=0; i < Tp->ntm; i++)
  {
    fprintf(out, "\nResult \"%ld-%s\" \"%ld\" %s Scalar OnNodes\n", i+1,
            namevartstr[int(Tp->dofname[i])-1].alias, lcid, desclcid);
    fprintf(out, "Values\n");

    for (j=0; j < Tt->nn; j++)
    {
      r = nodalval(j, i);
      fprintf(out, "%ld % e\n", j+Outdt->idn1, r);
    }
    fprintf(out, "End Values\n");
  }
}


void write_gid_nodvectort(FILE *out, strastret scal, long lcid, long unkn, const char *desclcid)
{
  long i, j;
  const char *sig = "";
  const char *unknown = "";

  if (scal == grad)
    sig = "gradients";
  if (scal == flux)
    sig = "fluxes";

  switch (unkn)
    {
    case 0:
      unknown = "1st";
      break;
    case 1:
      unknown = "2nd";
      break;
    case 2:
      unknown = "3rd";
      break;
    default :
      fprintf(stderr, "\n\n Error - unknown unknown of gradient/flux");
      fprintf(stderr, "\n in function write_gid_nodvectort, file %s, line %d\n", __FILE__, __LINE__);
    }

  switch (Tt->nodes[0].ncompgrad){//  transported matter
  case 1:{    
    fprintf(out, "\nResult \"%s - %s\" \"%ld\" %s Vector OnNodes\n",sig,unknown,lcid, desclcid);
    fprintf(out, "ComponentNames \"x\"\nValues\n");
    break;
  }
  case 2:{
    fprintf(out, "\nResult \"%s - %s\" \"%ld\" %s Vector OnNodes\n",sig,unknown,lcid, desclcid);
    fprintf(out, "ComponentNames \"x\" \"y\"\nValues\n");
    break;
  }
  case 3:{
    fprintf(out, "\nResult \"%s - %s\" \"%ld\" %s Vector OnNodes\n",sig,unknown,lcid, desclcid);
    fprintf(out, "ComponentNames \"x\" \"y\" \"z\"\nValues\n");
    break;
  }    
  default:{
    fprintf (stderr,"\n Error - unknown direction is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  

  if (scal == grad){
    for (i=0; i < Tt->nn; i++)
      {
	fprintf(out, "%ld", i+Outdt->idn1);
	
	for (j = 0; j < Tt->nodes[i].ncompgrad; j++)
	  fprintf(out, " % e",Tt->nodes[i].gradient[unkn][j]);

	fprintf(out, "\n");
      }
    fprintf(out, "End Values\n");
  }

  if (scal == flux){
    for (i=0; i < Tt->nn; i++)
      {
	fprintf(out, "%ld", i+Outdt->idn1);
	
	for (j = 0; j < Tt->nodes[i].ncompgrad; j++)
	  fprintf(out, " % e",Tt->nodes[i].flux[unkn][j]);

	fprintf(out, "\n");
      }
    fprintf(out, "End Values\n");
  }
}



/**
  The function writes vector of fluxes for all nodes to the file given by parameter out in GiD format. 
  The results are printed for all nodes with no dependency on the outdriver selection.

  Parameters :
  @param out - pointer to the opened text file where the output will be produced
  @param lcid - load case id
  @param desclcid - string with description of loadcase
  @param ifor - vector of nodal forces=fluxes

  Return :
  The function does not return anything

  created 06.2007 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  modified 10/01/2010 by Tomas Krejci, krejci@cml.fsv.cvut.cz
 */
void write_gid_nforcest(FILE *out, long lcid, const char *desclcid, double *ifor)
{
  long i, j, ii;
  double f;

  for (i=0; i < Tp->ntm; i++)
  {
    fprintf(out, "\nResult \"Fluxes %ld-%s\" \"%ld\" %s Scalar OnNodes\n", i+1,
            namevartstr[int(Tp->dofname[i])-1].alias, lcid, desclcid);
    fprintf(out, "Values\n");

    for (j=0; j < Tt->nn; j++)
    {
      ii=Tt->give_dof(j,i);
      if (ii<0)   f=0.0;
      if (ii==0)  f=0.0;
      if (ii>0)   f=ifor[ii-1];
      fprintf(out, "%ld % e\n", j+Outdt->idn1, f);
    }
    fprintf(out, "End Values\n");
  }
}



void write_gid_nodscalart(FILE *out, strastret scal, long lcid, long dir, const char *desclcid)
{
  const char *sig = "";
  long i, ncompother, ncompeqother;
  
  switch (scal)
    {
    case othert :
      {
	fprintf(out, "\nResult \"other_n_%ld - ", dir+1);
	Tm->give_othervalue_name(out,0,dir);
	fprintf(out, "\" \"%ld\" %s Scalar OnNodes\n", lcid, desclcid);
      }
      break;
    case eqothert :
      {
	fprintf(out, "\nResult \"eqother_n_%ld - ", dir+1);
	Tm->give_eqothervalue_name(out,0,dir);
	fprintf(out, "\" \"%ld\" %s Scalar OnNodes\n", lcid, desclcid);
      }
      break;
    default :
      fprintf(out, "\nResult \"%s\" \"%ld\" %s Scalar OnNodes\n", sig, lcid, desclcid);
    }
  fprintf(out, "Values\n");
  

  switch (scal)
    {
    case othert :
      for (i = 0; i < Tt->nn; i++)
	{
	  ncompother = Tt->nodes[i].ncompother;
	  fprintf(out, "%ld % e\n", i+Outdt->idn1, Tt->nodes[i].other[ncompother*lcid+dir]);
	}
      break;
    case eqothert :
      for (i = 0; i < Tt->nn; i++)
	{
	  ncompeqother = Tt->nodes[i].ncompeqother;
	  fprintf(out, "%ld % e\n", i+Outdt->idn1, Tt->nodes[i].eqother[ncompeqother*lcid+dir]);
	}
      break;
    default :
      fprintf(stderr, "\n\n Error - unknown value type in function write_gid_nodscalart()\n");
      fprintf(stderr, "  in file %s, line %d\n", __FILE__, __LINE__);
    }
  fprintf(out, "End Values\n");
}



/**
  The function writes a scalar quantity given by parameters scal and dir on all elements to the file given by parameter out in GiD format. 
  The results are printed at each integration point on given element.

  Parameters :
  @param out - pointer to the opened text file where the output will be produced
  @param scal - specifies type of required scalar quantity (grad/flux/other/eqother)
  @param lcid - load case id = in this case = number of unknown
  @param dir - specifies which component of the quantity array will be printed
  @param desclcid - string with description of loadcase

  Return :
  The function does not return anything
*/
void write_gid_elemscalart(FILE *out, strastret scal, long lcid, long dir, const char *desclcid)
{
  if (Lbt)
    write_gid_elem_type_scalart(out, scal, lcid, dir, desclcid, barlint);
  if (Lbt3d)
    write_gid_elem_type_scalart(out, scal, lcid, dir, desclcid, barlint3d);
  if (Lbat)
    write_gid_elem_type_scalart(out, scal, lcid, dir, desclcid, barlintax);
  if (Qbt)
    write_gid_elem_type_scalart(out, scal, lcid, dir, desclcid, barquadt);
  if (Qbat)
    write_gid_elem_type_scalart(out, scal, lcid, dir, desclcid, barquadtax);
  if (Ltt)
    write_gid_elem_type_scalart(out, scal, lcid, dir, desclcid, trlint);
  if (Ltat)
    write_gid_elem_type_scalart(out, scal, lcid, dir, desclcid, trlaxisym);
  if (Lqt)
    write_gid_elem_type_scalart(out, scal, lcid, dir, desclcid, quadlint);
  if (Lqat)
    write_gid_elem_type_scalart(out, scal, lcid, dir, desclcid, quadlaxisym);
  if (Qqt)
    write_gid_elem_type_scalart(out, scal, lcid, dir, desclcid, quadquadt);
  if (Qqat)
    write_gid_elem_type_scalart(out, scal, lcid, dir, desclcid, quadquadtax);
  if (Ltett)
    write_gid_elem_type_scalart(out, scal, lcid, dir, desclcid, lineartett);
  if (Lht)
    write_gid_elem_type_scalart(out, scal, lcid, dir, desclcid, linearhext);
  if (Qht)
    write_gid_elem_type_scalart(out, scal, lcid, dir, desclcid, quadratichext);
}



/**
  The function writes a scalar value given by parameter scal on elements of given type (parameter te)
  to the file given by parameter out in GiD format. The results are printed at each  element integration 
  points.

  Parameters :
  @param out - pointer to the opened text file where the output will be produced
  @param scal - specifies type of required scalar quantity (grad/flux/other/eqother)
  @param lcid - load case id = in this case = number of unknown
  @param dir - specifies which component of the quantity array will be printed
  @param desclcid - string with description of loadcase
  @param te - required element type

  Return :
  The function does not return anything
*/
void write_gid_elem_type_scalart(FILE *out, strastret scal, long lcid, long dir, const char *desclcid, elemtypet te)
{
  const char *sig = "";
  char gpname[1000];
  long i, j, tnip, ipp,ncompother,ncompeqother;
  long print_header = 1;


  for (i = 0; i < Tt->ne; i++)
  {
    if (Gtt->leso[i]==0)
      continue;


    // checking required element type
    if (Tt->elements[i].te != te)
      continue; 

    // checking required scalar type and selection of the element and selection of the component
    if ((scal == grad) && (Outdt->eog.selegrad.presence_id(Outdt->eog.selgrad, i, dir) == 0))
      continue;
    if ((scal == flux) && (Outdt->eog.seleflux.presence_id(Outdt->eog.selflux, i, dir) == 0))
      continue;
    if ((scal == othert) && (Outdt->eog.seleoth.presence_id(Outdt->eog.seloth, i, dir) == 0))
      continue;
    if ((scal == eqothert) && (Outdt->eog.seleeqoth.presence_id(Outdt->eog.seleqoth, i, dir) == 0))
      continue;

    if (print_header)
    {
      ipp = Tt->elements[i].ipp[0][0];
      if (scal == flux)
      {
        switch (te)
	  {
	  case barlintax:
	  case barlint:
	  case barlint3d:
	  case barquadt:
	  case barquadtax:
            switch (dir)
            {
              case 0:
		if(lcid==0)
		  sig = "flux_e_1_x";
		if(lcid==1)
		  sig = "flux_e_2_x";		
		if(lcid==2)
		  sig = "flux_e_3_x";		
		if(lcid==3)
		  sig = "flux_e_4_x";		
                break;
              default :
                fprintf(stderr, "\n\n Error - unknown direction of flux");
                fprintf(stderr, "\n in function  write_gid_elemscalart, file %s, line %d\n", __FILE__, __LINE__);
            }
            break;
	  case trlint:
	  case trlaxisym:
	  case trquadt://prepared for quadratic elements
	  case trqaxisym:
	  case quadlint:
	  case quadlaxisym:
	  case quadquadtax:
	  case quadquadt:
            switch (dir)
            {
              case 0:
		if(lcid==0)
		  sig = "flux_e_1_x";
		if(lcid==1)
		  sig = "flux_e_2_x";		
		if(lcid==2)
		  sig = "flux_e_3_x";		
		if(lcid==3)
		  sig = "flux_e_4_x";		
                break;
              case 1:
		if(lcid==0)
		  sig = "flux_e_1_y";
		if(lcid==1)
		  sig = "flux_e_2_y";		
		if(lcid==2)
		  sig = "flux_e_3_y";		
		if(lcid==3)
		  sig = "flux_e_4_y";		
                break;
              default :
                fprintf(stderr, "\n\n Error - unknown direction of flux");
                fprintf(stderr, "\n in function  write_gid_elemscalart, file %s, line %d\n", __FILE__, __LINE__);
            }
            break;
	  case lineartett:
	  case linearhext:
	  case quadratichext:
            switch (dir)
	      {
              case 0:
		if(lcid==0)
		  sig = "flux_e_1_x";
		if(lcid==1)
		  sig = "flux_e_2_x";		
		if(lcid==2)
		  sig = "flux_e_3_x";		
		if(lcid==3)
		  sig = "flux_e_4_x";		
                break;
              case 1:
		if(lcid==0)
		  sig = "flux_e_1_y";
		if(lcid==1)
		  sig = "flux_e_2_y";		
		if(lcid==2)
		  sig = "flux_e_3_y";		
		if(lcid==3)
		  sig = "flux_e_4_y";		
                break;
              case 2:
		if(lcid==0)
		  sig = "flux_e_1_z";
		if(lcid==1)
		  sig = "flux_e_2_z";		
		if(lcid==2)
		  sig = "flux_e_3_z";		
		if(lcid==3)
		  sig = "flux_e_4_z";		
                break;
              default :
                fprintf(stderr, "\n\n Error - unknown direction of flux");
                fprintf(stderr, "\n in function  write_gid_elemscalart, file %s, line %d\n", __FILE__, __LINE__);
            }
            break;
          default:
            fprintf(stderr, "\n\n Error - unknown flux/grad state is required");
            fprintf(stderr, "\n in function  write_gid_elemscalart, file %s, line %d\n", __FILE__, __LINE__);
        }
      }

      if (scal == grad)
      {
        switch (te)
	  {
	  case barlintax:
	  case barlint:
	  case barlint3d:
	  case barquadt:
	  case barquadtax:
            switch (dir)
            {
              case 0:
		if(lcid==0)
		  sig = "grad_e_1_x";
		if(lcid==1)
		  sig = "grad_e_2_x";		
		if(lcid==2)
		  sig = "grad_e_3_x";		
		if(lcid==3)
		  sig = "grad_e_4_x";		
                break;
              default :
                fprintf(stderr, "\n\n Error - unknown direction of grad");
                fprintf(stderr, "\n in function  write_gid_elemscalart, file %s, line %d\n", __FILE__, __LINE__);
            }
            break;
	  case trlint:
	  case trlaxisym:
	  case trquadt://prepared for quadratic elements
	  case trqaxisym:
	  case quadlint:
	  case quadlaxisym:
	  case quadquadtax:
	  case quadquadt:
            switch (dir)
            {
              case 0:
		if(lcid==0)
		  sig = "grad_e_1_x";
		if(lcid==1)
		  sig = "grad_e_2_x";		
		if(lcid==2)
		  sig = "grad_e_3_x";		
		if(lcid==3)
		  sig = "grad_e_4_x";		
                break;
              case 1:
		if(lcid==0)
		  sig = "grad_e_1_y";
		if(lcid==1)
		  sig = "grad_e_2_y";		
		if(lcid==2)
		  sig = "grad_e_3_y";		
		if(lcid==3)
		  sig = "grad_e_4_y";		
                break;
              default :
                fprintf(stderr, "\n\n Error - unknown direction of grad");
                fprintf(stderr, "\n in function  write_gid_elemscalart, file %s, line %d\n", __FILE__, __LINE__);
            }
            break;
	  case lineartett:
	  case linearhext:
	  case quadratichext:
            switch (dir)
	      {
              case 0:
		if(lcid==0)
		  sig = "grad_e_1_x";
		if(lcid==1)
		  sig = "grad_e_2_x";		
		if(lcid==2)
		  sig = "grad_e_3_x";		
		if(lcid==3)
		  sig = "grad_e_4_x";		
                break;
              case 1:
		if(lcid==0)
		  sig = "grad_e_1_y";
		if(lcid==1)
		  sig = "grad_e_2_y";		
		if(lcid==2)
		  sig = "grad_e_3_y";		
		if(lcid==3)
		  sig = "grad_e_4_y";		
                break;
              case 2:
		if(lcid==0)
		  sig = "grad_e_1_z";
		if(lcid==1)
		  sig = "grad_e_2_z";		
		if(lcid==2)
		  sig = "grad_e_3_z";		
		if(lcid==3)
		  sig = "grad_e_4_z";		
                break;
              default :
                fprintf(stderr, "\n\n Error - unknown direction of grad");
                fprintf(stderr, "\n in function  write_gid_elemscalart, file %s, line %d\n", __FILE__, __LINE__);
            }
            break;
          default:
            fprintf(stderr, "\n\n Error - unknown flux/grad state is required");
            fprintf(stderr, "\n in function  write_gid_elemscalart, file %s, line %d\n", __FILE__, __LINE__);
        }
      }

      // element name //
      switch (te)
	{
	case barlintax:
	case barlint:
	case barlint3d:
          sprintf(gpname, "Lin_1D");
          break;
	case barquadt:
	case barquadtax:
          sprintf(gpname, "Quad_1D");
          break;
        default :
          sprintf(gpname, "%d", te);
	}
      

      // number of integration points on elements only for the first block
      switch (te)
	{
	case barlintax:
	case barlint:
	case barlint3d:
          tnip = 2;
          break;
	case barquadt:
	case barquadtax:
          tnip = 3;
          break;
	case trlint:
	case trlaxisym:
          tnip = 1;
          break;
	case quadlint:
	case quadlaxisym:
	  tnip = 4;
          break;
	case quadquadtax:
	case quadquadt:
	  tnip = 9;
          break;
	case lineartett:
	  tnip = 1;
          break;
	case linearhext:
	  tnip = 8;
          break;
	case quadratichext:
	  tnip = 27;
          break;
        default :
          fprintf(stderr, "\n Unknown element type in function write_gid_elem_type_scalart\n");
	  fprintf(stderr, "  in file %s, line %d\n", __FILE__, __LINE__);
	}
      

      switch (scal)
	{
	case othert:
	  fprintf(out, "\nResult \"other_e_%ld\" \"%ld\" %s Scalar OnGaussPoints \"%s\"\nValues\n", dir+1, lcid, desclcid, gpname);
	  break;
	case eqothert:
	  fprintf(out, "\nResult \"eqother_e_%ld\" \"%ld\" %s Scalar OnGaussPoints \"%s\"\nValues\n", dir+1, lcid, desclcid, gpname);
	  break;
	default:
	  fprintf(out, "\nResult \"%s\" \"%ld\" %s Scalar OnGaussPoints \"%s\"\nValues\n", sig, lcid, desclcid, gpname);
	}

      print_header = 0;
    }

    switch (scal)
    {

      case flux :

        fprintf(out, "%7ld", i+Outdt->ide1);
        for (j = 0; j < tnip; j++)
	  {
	    ipp = Tt->elements[i].ipp[0][0]+j;
	    if (j == 0)
	      fprintf(out, " % e\n", Tm->ip[ipp].fluxes[lcid][dir]);
	    else
	      fprintf(out, "%7c % e\n", ' ', Tm->ip[ipp].fluxes[lcid][dir]);
	  }
        fprintf(out, "\n");
        break;
	
    case grad :
      
      fprintf(out, "%7ld", i+Outdt->ide1);
      for (j = 0; j < tnip; j++)
        {
          ipp = Tt->elements[i].ipp[0][0]+j;
          if (j == 0)
	    fprintf(out, " % e\n", Tm->ip[ipp].grad[lcid][dir]);
	    //fprintf (Outt,"\n ipp %ld   ncompstr %ld  lcid %ld  dir %ld  [] %ld",ipp,ncompgrad,lcid,dir,ncompgrad*lcid+dir);
          else
	    fprintf(out, "%7c % e\n", ' ', Tm->ip[ipp].grad[lcid][dir]);
        }
      fprintf(out, "\n");
      break;
      
      case othert :
        fprintf(out, "%7ld", i+Outdt->ide1);
        for (j = 0; j < tnip; j++)
        {
          ipp = Tt->elements[i].ipp[0][0]+j;
          ncompother = Tm->ip[ipp].ncompother;
          if (j == 0)
            fprintf(out, " % e\n", Tm->ip[ipp].other[ncompother*lcid+dir]);
          else
            fprintf(out, "%7c % e\n", ' ', Tm->ip[ipp].other[ncompother*lcid+dir]);
        }
        fprintf(out, "\n");
        break;

      case eqothert :
        fprintf(out, "%7ld", i+Outdt->ide1);
        for (j = 0; j < tnip; j++)
        {
          ipp = Tt->elements[i].ipp[0][0]+j;
          ncompeqother = Tm->ip[ipp].ncompeqother;
          if (j == 0)
            fprintf(out, " % e\n", Tm->ip[ipp].eqother[ncompeqother*lcid+dir]);
          else
            fprintf(out, "%7c % e\n", ' ', Tm->ip[ipp].eqother[ncompeqother*lcid+dir]);
        }
        fprintf(out, "\n");
        break;

      default :
        fprintf(stderr, "\n\n Error - unknown value type in function write_gid_elemscalart()\n");
        fprintf(stderr, "  in file %s, line %d\n", __FILE__, __LINE__);
    }
  }
  if (print_header == 0)
    fprintf(out, "End Values\n");
}



/**
  The function writes a vector quantity given by parameters scal and dir on all elements to the file given by parameter out in GiD format. 
  The results are printed at each integration point on given element.

  Parameters :
  @param out - pointer to the opened text file where the output will be produced
  @param q - specifies type of required vector quantity (grad/flux/other/eqother)
  @param lcid - load case id = in this case = number of unknown
  @param dir - specifies which component of the quantity array will be printed
  @param desclcid - string with description of loadcase

  Return :
  The function does not return anything
*/
void write_gid_elemvectort(FILE *out, strastret q, long lcid, long dir, const char *desclcid)
{
  if (Lbt)
    write_gid_elem_type_vectort(out, q, lcid, dir, desclcid, barlint);
  if (Lbt3d)
    write_gid_elem_type_vectort(out, q, lcid, dir, desclcid, barlint3d);
  if (Lbat)
    write_gid_elem_type_vectort(out, q, lcid, dir, desclcid, barlintax);
  if (Qbt)
    write_gid_elem_type_vectort(out, q, lcid, dir, desclcid, barquadt);
  if (Qbat)
    write_gid_elem_type_vectort(out, q, lcid, dir, desclcid, barquadtax);
  if (Ltt)
    write_gid_elem_type_vectort(out, q, lcid, dir, desclcid, trlint);
  if (Ltat)
    write_gid_elem_type_vectort(out, q, lcid, dir, desclcid, trlaxisym);
  if (Lqt)
    write_gid_elem_type_vectort(out, q, lcid, dir, desclcid, quadlint);
  if (Lqat)
    write_gid_elem_type_vectort(out, q, lcid, dir, desclcid, quadlaxisym);
  if (Qqt)
    write_gid_elem_type_vectort(out, q, lcid, dir, desclcid, quadquadt);
  if (Qqat)
    write_gid_elem_type_vectort(out, q, lcid, dir, desclcid, quadquadtax);
  if (Ltett)
    write_gid_elem_type_vectort(out, q, lcid, dir, desclcid, lineartett);
  if (Lht)
    write_gid_elem_type_vectort(out, q, lcid, dir, desclcid, linearhext);
  if (Qht)
    write_gid_elem_type_vectort(out, q, lcid, dir, desclcid, quadratichext);
}



/**
  The function writes a scalar value given by parameter scal on elements of given type (parameter te)
  to the file given by parameter out in GiD format. The results are printed at each  element integration 
  points.

  Parameters :
  @param out - pointer to the opened text file where the output will be produced
  @param q - specifies type of required vector quantity (grad/flux/other/eqother)
  @param lcid - load case id = in this case = number of unknown
  @param sid - index of selection of elements which will be printed
  @param desclcid - string with description of loadcase
  @param te - required element type

  Return :
  The function does not return anything
*/
void write_gid_elem_type_vectort(FILE *out, strastret q, long lcid, long sid, const char *desclcid, elemtypet te)
{
  char sig[70];
  long q_id1, q_n;
  char gpname[1000];
  long i, j, k, tnip, ipp;
  long print_header = 1;


  switch (q)
  {
    case grad:
      q_id1 = Outdt->eog.selgrad[sid].id1[0];
      q_n = Tp->gdim;
      sprintf(sig, "grad_e_v%ld_s%ld", q_id1+1, sid+1);
      break;
    case flux:
      q_id1 = Outdt->eog.selflux[sid].id1[0];
      q_n = Tp->gdim;
      sprintf(sig, "flux_e_v%ld_s%ld", q_id1+1, sid+1);
      break;
    case othert:
      q_id1 = Outdt->eog.seloth[sid].id1[0];
      q_n = Outdt->eog.seloth[sid].ncomp[0];
      sprintf(sig, "other_e_v%ld-%ld_s%ld", q_id1+1, q_n, sid+1);
      break;
    case eqothert:
      q_id1 = Outdt->eog.seloth[sid].id1[0];
      q_n = Outdt->eog.seloth[sid].ncomp[0];
      sprintf(sig, "other_e_v%ld-%ld_s%ld", q_id1+1, q_n, sid+1);
      break;
    default:
      print_err("unsupported type of quantity (%d) is required in function\n", __FILE__, __LINE__, __func__, q);
      abort();
  }

  for (i = 0; i < Tt->ne; i++)
  {
    if (Gtt->leso[i]==0)
      continue;


    // checking required element type
    if (Tt->elements[i].te != te)
      continue; 

    // checking required quantity type, selection of the element and selection of the vector
    if ((q == grad) && (Outdt->eog.selegrad.presence_id(Outdt->eog.selgrad, i, sid) == 0))
      continue;
    if ((q == flux) && (Outdt->eog.seleflux.presence_id(Outdt->eog.selflux, i, sid) == 0))
      continue;
    if ((q == othert) && (Outdt->eog.seleoth.presence_id(Outdt->eog.seloth, i, sid) == 0))
      continue;
    if ((q == eqothert) && (Outdt->eog.seleeqoth.presence_id(Outdt->eog.seleqoth, i, sid) == 0))
      continue;

    if (print_header)
    {
      ipp = Tt->elements[i].ipp[0][0];

      // element name //
      switch (te)
      {
	case barlintax:
	case barlint:
	case barlint3d:
          sprintf(gpname, "Lin_1D");
          break;
	case barquadt:
	case barquadtax:
          sprintf(gpname, "Quad_1D");
          break;
        default :
          sprintf(gpname, "%d", te);
      }
      
      tnip = Tt->give_tnip(i)/Tp->ntm;      
      fprintf(out, "\nResult \"%s\" \"%ld\" %s Vector OnGaussPoints \"%s\"\n", sig, lcid, desclcid, gpname);
      fprintf(out, "Values\n");
      print_header = 0;
    }

    switch (q)
    {

      case flux :
        for (j = 0; j < tnip; j++)
	{
          ipp = Tt->elements[i].ipp[0][0]+j;
          if (j == 0)
            fprintf(out, "%7ld", i+Outdt->ide1);
          else
            fprintf(out, "%7c", ' ');
          for (k = 0; k < q_n; k++)
            fprintf(out, " % e", Tm->ip[ipp].fluxes[q_id1][k]);
          for (k = 0; k < 3-q_n; k++)
             fprintf(out, " 0.0");
          fprintf(out, "\n");
        }
        fprintf(out, "\n");
        break;
	
      case grad :
        for (j = 0; j < tnip; j++)
	{
          ipp = Tt->elements[i].ipp[0][0]+j;
          if (j == 0)
            fprintf(out, "%7ld", i+Outdt->ide1);
          else
            fprintf(out, "%7c", ' ');
          for (k = 0; k < q_n; k++)
            fprintf(out, " % e", Tm->ip[ipp].grad[q_id1][k]);
          for (k = 0; k < 3-q_n; k++)
             fprintf(out, " 0.0");
          fprintf(out, "\n");
        }
        fprintf(out, "\n");
        break;

      case othert :
        for (j = 0; j < tnip; j++)
	{
          ipp = Tt->elements[i].ipp[0][0]+j;
          if (j == 0)
            fprintf(out, "%7ld", i+Outdt->ide1);
          else
            fprintf(out, "%7c", ' ');
          for (k = q_id1; k < q_id1+q_n; k++)
            fprintf(out, " % e", Tm->ip[ipp].other[k]);
          for (k = 0; k < 3-q_n; k++)
             fprintf(out, " 0.0");
          fprintf(out, "\n");
        }
        fprintf(out, "\n");
        break;

      case eqothert :
        for (j = 0; j < tnip; j++)
	{
          ipp = Tt->elements[i].ipp[0][0]+j;
          if (j == 0)
            fprintf(out, "%7ld", i+Outdt->ide1);
          else
            fprintf(out, "%7c", ' ');
          for (k = q_id1; k < q_id1+q_n; k++)
            fprintf(out, " % e", Tm->ip[ipp].eqother[k]);
          for (k = 0; k < 3-q_n; k++)
             fprintf(out, " 0.0");
          fprintf(out, "\n");
        }
        fprintf(out, "\n");
        break;

      default :
        fprintf(stderr, "\n\n Error - unknown value type in function write_gid_elemscalart()\n");
        fprintf(stderr, "  in file %s, line %d\n", __FILE__, __LINE__);
    }
  }
  if (print_header == 0)
    fprintf(out, "End Values\n");
}



void write_nforcest(FILE *out, long /*lcid*/, const char *desclcid, double *ifor)
{
  long i, j, ii, ndof;
  vector f, g;

  fprintf (out,"\nVALUES");
  fprintf (out,"\nVECTOR3D NODES FLUXES %s\n", desclcid);
  for (i=0; i < Tt->nn; i++)
  {
    ndof = Tt->give_ndofn(i);
    reallocv(RSTCKVEC(ndof, f));
    for (j=0;j<ndof;j++)
    {
      ii=Tt->give_dof(i,j);
      if (ii<0)   f[j]=0.0;
      if (ii==0)  f[j]=0.0;
      if (ii>0)   f[j]=ifor[ii-1];

      //if (ii<1)  f[j]=Tt->nodes[i].nodval[j];
      //if (ii>0)  f[j]=ifor[ii-1];
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
  The function prints header, points and used elements to the VTK file.

  Parameters :
  @param fname - pointer to the opened VTK file
  @param lcid  - load case id
  @param istep - integration step
  @param lambda - actual time in either s,h,days

  @return The function does not return anything.

  Created by V. Smilauer, 2006
  Modified by Tomas Koudelka, 02.2018  
*/
void print_default_vtk (FILE *fname, long /*lcid*/, long istep, double lambda)
{
  long int i,j,nn;
  int tr[20];
  
  fprintf(fname, "# vtk DataFile Version 1.0\n");
  fprintf(fname, "Time step %ld time %lf\n", istep, lambda);
  fprintf(fname, "ASCII\n\n");
  fprintf(fname, "DATASET UNSTRUCTURED_GRID\n");
  
  fprintf(fname, "POINTS %ld double\n", Tt->nn);
  for (i=0; i<Tt->nn; i++) {
    fprintf(fname, "%e %e %e\n", Gtt->gnodes[i].x, Gtt->gnodes[i].y, Gtt->gnodes[i].z);
  }
  
  long int ncells = 0;
  for (i=0; i < Tt->ne; i++){
    ncells += Tt->give_nne (i);
  }
  fprintf(fname, "\nCELLS %ld %ld\n", Tt->ne, ncells+Tt->ne);
  for (i=0; i < Tt->ne; i++){
    fprintf(fname, "%ld ", Tt->give_nne(i));//number of nodes on element
    for(j=0;j<Tt->give_nne(i);j++){
      switch(Tt->elements[i].te){//for nodes renumbering
        case lineartett:
          tr[0] = 0;  tr[1] = 1;  tr[2] = 3;  tr[3] = 2;
          nn=tr[j];
          break;
        case linearhext:
          tr[0] = 6; tr[1] = 7; tr[2] = 4; tr[3] = 5;//xy, z=0, xi=x,nu=y, theta=z
          tr[4] = 2; tr[5] = 3; tr[6] = 0; tr[7] = 1;//xy, z<>0
          nn=tr[j];
          break;
        case quadratichext:
          tr[0] = 6; tr[1] = 7; tr[2] = 4; tr[3] = 5;//xy, z=0, xi=x,nu=y, theta=z
          tr[4] = 2; tr[5] = 3; tr[6] = 0; tr[7] = 1;
          tr[8] = 18; tr[9] = 19; tr[10] = 16; tr[11] = 17;
          tr[12] = 10; tr[13] = 11; tr[14] = 8; tr[15] = 9;
          tr[16] = 14; tr[17] = 15; tr[18] = 12; tr[19] = 13;
          nn=tr[j];
          break;
        default:
          nn=j;
          break;
      }
      fprintf(fname, "%ld ", Gtt->gelements[i].nodes[nn]);
    }
    fprintf(fname, "\n");
  }
  
  //output cell types
  fprintf(fname, "\nCELL_TYPES %ld\n", Tt->ne);
  for (i=0; i < Tt->ne; i++){
    switch(Tt->elements[i].te){
      case barlint:
      case barlint3d:
      case barlintax :
        fprintf(fname, "3\n");//VTK_LINE
        break;
      case barquadt :
      case barquadtax :
        fprintf(fname, "21\n");//VTK_QUADRATIC_EDGE
        break;
      case trlint:
      case trlaxisym:
        fprintf(fname, "5\n");//VTK_TRIANGLE
        break;
      case quadlint:
      case  quadlaxisym:
        fprintf(fname, "9\n");//VTK_QUAD
        break;
      case quadquadt:
      case quadquadtax:
        fprintf(fname, "23\n");//VTK_QUADRATIC_QUAD
        break;
      case lineartett:
        fprintf(fname, "10\n");//VTK_TETRA
        break;
      case linearhext:
        fprintf(fname, "12\n");//VTK_HEXAHEDRON - Renumbered nodes
        break;
      case quadratichext:
        fprintf(fname, "25\n");//VTK_QUADRATIC_HEXAHEDRON - Renumbered nodes
        break;
      default:
        fprintf(stderr, "\n\n Unknown element type in VTK export, file %s, line %d\n", __FILE__, __LINE__);
        abort();
        break;
    }
  }
}



/**
  The function prints datasets in VTK format to the given file.
  VTK works smoothly only on nodes and cells, NOT integration points.

  Parameters :
  @param fname - pointer to the opened VTK file [in]
  @param lcid  - load case id [in]
  @param flag_material - flag for printing material types at particular elements
  @param fi - array of nodal flux resultants (load %vector)

  @return The function does not return anything.

  Created by V. Smilauer, 2006
  Modified by Tomas Koudelka, 02.2018  
*/
void write_vtk_unkn(FILE *fname,long /*lcid*/, int flag_material, double *fi)
{
  vector r;
  long i, j, k, ipp;
  char VarName[255];
  //  double VTK_no_value = -0.11111;//cell or point has no value
  double VTK_no_value = -0.0;//cell or point has no value
  long max_ncomp;
  
  fprintf(fname, "\nPOINT_DATA %ld\n", Tt->nn);
  
  for (j = 0; j < Tp->ntm; j++)
  {
    fprintf(fname, "SCALARS %s double\nLOOKUP_TABLE default\n", namevart_kwdset.get_str(Tp->dofname[j]));
    for (i=0; i < Tt->nn; i++)
    {
      reallocv(RSTCKVEC(Tt->give_ndofn(i), r));
      nodalval(i, r);
      
      //fprintf(fname, "%e", r[j]-273.15);
      if (Outdt->nog.selnunkn.presence_id(Outdt->nog.selunkn, i, j))
        fprintf(fname, "%e\n", r[j]);
      else
        fprintf(fname, "%e\n", VTK_no_value);
    }
  }

  //Point=nodal gradients
  if (Outdt->nog.selngrad.st != sel_no){
    for(i=0; i<Tp->ntm; i++)
    {
      fprintf(fname, "\nVECTORS Grad_q%ld double\n", i+1);
      for (j=0; j<Tt->nn; j++){
        if (Outdt->nog.selngrad.presence_id(j)){
          for (k=0; k<Tt->nodes[j].ncompgrad; k++){
            if (Outdt->nog.selngrad.presence_id(Outdt->nog.selgrad, j, k))
              fprintf(fname, "%e ", Tt->nodes[j].givegrad(i, k));
            else
              fprintf(fname, "%e ", VTK_no_value);
          }
          if (Tt->nodes[j].ncompgrad < 3)
          {
            for (k=0; k<3-Tt->nodes[j].ncompgrad; k++)
              fprintf(fname, "%e ", VTK_no_value);
          }
          fprintf(fname, "\n");
        }
        else
          fprintf(fname, "%e %e %e\n", VTK_no_value, VTK_no_value, VTK_no_value);
      }
    }
  }
  
  //Point=nodal fluxes
  if (Outdt->nog.selnflux.st != sel_no){
    for(i=0; i<Tp->ntm; i++)
    {
      fprintf(fname, "\nVECTORS Flux_q%ld double\n", i+1);
      for (j=0; j<Tt->nn; j++){
        if (Outdt->nog.selnflux.presence_id(j)){
          for (k=0; k<Tt->nodes[j].ncompgrad; k++){
            if (Outdt->nog.selnflux.presence_id(Outdt->nog.selflux, j, k))
              fprintf(fname, "%e ", Tt->nodes[j].giveflux(i, k));
            else
              fprintf(fname, "%e ", VTK_no_value);
          }
          if (Tt->nodes[j].ncompgrad < 3)
          {
            for (k=0; k<3-Tt->nodes[j].ncompgrad; k++)
              fprintf(fname, "%e ", VTK_no_value);
          }
          fprintf(fname, "\n");
        }
        else
          fprintf(fname, "%e %e %e\n", VTK_no_value, VTK_no_value, VTK_no_value);
      }
    }
  }
  
  //print other components on all points=nodes
  if (Outdt->nog.selnoth.st != sel_no){
    for (i=0; i<Tm->givencompother(); i++){//other[i]
      fprintf(fname, "\nSCALARS Other_%ld double\nLOOKUP_TABLE default\n", i+1);
      for (k=0; k<Tt->nn; k++){
        if (Outdt->nog.selnoth.presence_id(k))
          fprintf(fname, "%e\n", Tt->nodes[k].giveother(i));
        else
          fprintf(fname, "%e\n", VTK_no_value);
      }
    }
  }
  
  //print eqother components on all points=nodes
  if (Outdt->nog.selneqoth.st != sel_no){
    max_ncomp = 0;
    for (k=0; k<Tm->tnip; k++){//find maximum eqother components
      max_ncomp > Tm->ip[k].ncompeqother ? max_ncomp = Tm->ip[k].ncompeqother : 0;
    }
    
    for(i=0; i<max_ncomp; i++){// eqother[i]
      for (j=0; j<Outdt->nog.selneqoth.n; j++){
        if (Outdt->nog.seleqoth[j].presence_id(i)){
          fprintf(fname, "\nSCALARS Eq_other_s%ld_%ld double\nLOOKUP_TABLE default\n", j+1, i+1);
          for (k=0; k < Tt->nn; k++)
            fprintf(fname, "%e\n", Tt->nodes[k].giveeqother(i));
        }
      }
    }
  }

  if (Outdt->nog.selnforcet.st != sel_no)
  {
    for (j = 0; j < Tp->ntm; j++)
    {
      fprintf(fname, "SCALARS %s_flux_resultant double\nLOOKUP_TABLE default\n", Tp->dof_mednam[j]);
      for (i=0; i < Tt->nn; i++)
      {
        reallocv(RSTCKVEC(Gtt->gnodes[i].ndofn, r));
        gen_gvnodval(fi, i, r);
        //fprintf(fname, "%e", r[j]-273.15);
        if (Outdt->nog.selnforcet.presence_id(Outdt->nog.selforcet, i, j))
          fprintf(fname, "%e\n", r[j]);
        else
          fprintf(fname, "%e\n", VTK_no_value);
      }
    }
  }

  //
  // CELL DATA printing
  //

  if(flag_material){
    //print materials (usually in the first VTK file), need to extend from tm[0] and idm[0]
    fprintf(fname, "\nCELL_DATA %ld\n", Tt->ne);
    fprintf(fname, "SCALARS Mat_type int\nLOOKUP_TABLE default\n");
    for (i=0; i < Tt->ne; i++){
      fprintf(fname, "%d\n", Tt->elements[i].tm[0]);
    }

    fprintf(fname, "\nSCALARS Mat_id int\nLOOKUP_TABLE default\n");
    for (i=0; i < Tt->ne; i++){
      fprintf(fname, "%ld\n", Tt->elements[i].idm[0]+1);
    }
  }
  else{
    fprintf(fname, "\nCELL_DATA %ld\n", Tt->ne);
  }

  //Cell=element gradients
  if (Outdt->eog.selegrad.st != sel_no){
    fprintf(fname, "\nVECTORS Gradient double\n");
    for (i=0; i<Tt->ne; i++){
      if (Outdt->eog.selegrad.presence_id(i)){
        ipp = Tt->elements[i].ipp[0][0];
        fprintf(fname, "%e %e %e\n", Outdt->eog.selegrad.presence_id(Outdt->eog.selgrad, i, 0)!=0?Tm->ip[ipp].grad[0][0]:0, Outdt->eog.selegrad.presence_id(Outdt->eog.selgrad, i, 1)!=0?Tm->ip[ipp].grad[1][0]:0, Outdt->eog.selegrad.presence_id(Outdt->eog.selgrad, i, 2)!=0?Tm->ip[ipp].grad[2][0]:0);
      }
      else
        fprintf(fname, "%e %e %e\n", VTK_no_value, VTK_no_value, VTK_no_value);
    }
  }
  
  //Cell=element fluxes
  if (Outdt->eog.seleflux.st != sel_no){
    fprintf(fname, "\nVECTORS Flux double\n");
    for (i=0; i<Tt->ne; i++){
      if (Outdt->eog.seleflux.presence_id(i)){
        ipp = Tt->elements[i].ipp[0][0];
        fprintf(fname, "%e %e %e\n", Outdt->eog.seleflux.presence_id(Outdt->eog.selflux, i, 0)!=0?Tm->ip[ipp].fluxes[0][0]:0, Outdt->eog.seleflux.presence_id(Outdt->eog.selflux, i, 1)!=0?Tm->ip[ipp].fluxes[1][0]:0, Outdt->eog.seleflux.presence_id(Outdt->eog.selflux, i, 2)!=0?Tm->ip[ipp].fluxes[2][0]:0);
      }
      else
        fprintf(fname, "%e %e %e\n", VTK_no_value, VTK_no_value, VTK_no_value);
    }
  }
  
  //print other components on all cells=elements
  if (Outdt->eog.seleoth.st != sel_no){
    for (i=0; i<Tm->givencompother(); i++){//other[i]
      fprintf(fname, "\nSCALARS Other_%ld double\nLOOKUP_TABLE default\n", i);
        for (k=0; k<Tt->ne; k++){
          if (Outdt->eog.seleoth.presence_id(k)){
            ipp = Tt->elements[k].ipp[0][0];
            fprintf(fname, "%e\n", Tm->ip[ipp].other[i]);
          }
          else{
            fprintf(fname, "%e\n", VTK_no_value);
          }
        }
    }
  }
  
  
  //print eqother components on all cells=elements
  if (Outdt->eog.seleeqoth.st != sel_no){
    for (k=0; k<Tm->tnip; k++){//find maximum eqother components
      max_ncomp > Tm->ip[k].ncompeqother ? max_ncomp = Tm->ip[k].ncompeqother : 0;
    }
    
    for(i=0; i<max_ncomp; i++){// eqother[i]
      for (j=0; j<Outdt->eog.seleeqoth.n; j++){
        if (Outdt->eog.seleqoth[j].presence_id(i)){
          
          if(Tt->elements[0].tm[0] == cementhydrmat){
            if (i==0) strcpy(VarName, "Rel_heat");
            if (i==1) strcpy(VarName, "DoH");
            if (i==2) strcpy(VarName, "E_paste");
            if (i==3) strcpy(VarName, "nu_paste");
            if (i==4) strcpy(VarName, "E_concrete");
            if (i==5) strcpy(VarName, "nu_concrete");
          }
          else snprintf(VarName, 254, "Eq_other_s%ld", j+1);
          
          fprintf(fname, "\nSCALARS %s_%ld double\nLOOKUP_TABLE default\n", VarName, i+1);
          for (k=0; k < Tt->ne; k++){
            ipp = Tt->elements[k].ipp[0][0];
            fprintf(fname, "%e\n", Tm->ip[ipp].eqother[i]);
          }
        }
      }
    }
  }
}


