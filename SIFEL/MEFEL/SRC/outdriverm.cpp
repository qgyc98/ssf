#include "outdriverm.h"
#include "outquantm.h"
#include "gfmatrix.h"
#include "iotools.h"
#include "galias.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechbclc.h"
#include "lhsrhs.h"
#include "mechprint.h"
#include "node.h"
#include "element.h"
#include "globmat.h"
#include "intpoints.h"
#include "elemswitch.h"
#include "vecttens.h"
#include "list.h"
#include "outresfilem.h"
#include "gnodvalvm.h"
#include "outquantm.h"
#include <math.h>



static int prdisp = 15;
static int prstra = 15;
static int prstre = 15;
static int proth  = 15;
static int prreac = 15;



/**
  Constructor initializes data members to zero values.

  Created by Tomas Koudelka
*/
outdriverm::outdriverm()
{
  memset(outfn, 0, sizeof(*outfn)*FNAMELEN);
  memset(outdiagfn, 0, sizeof(*outdiagfn)*FNAMELEN);
  memset(outgrfn, 0, sizeof(*outgrfn)*FNAMELEN);
  outf = outgr = NULL;
  outdiagf = NULL;
  textout = off;
  gf = grfmt_no;
  ncut = 0;
  idn1 = 1;
  ide1 = 1;
  ndiag = 0;
  odiag = NULL;
  vtk_num = 0;
  nlcs = 0L;
  lcs = NULL;
  nofq = 0;
  ofq = NULL;
}



/**
  Destructor deallocates used memory.

  Created by Tomas Koudelka
*/
outdriverm::~outdriverm()
{
  delete [] outdiagf;
  delete [] odiag;
  delete [] lcs;
  delete [] ofq;
}



/**
  Function reads description of required output from the opened text file in.

  @param in - pointer to opened text file

  @retval 0 - on success
  @retval 1 - error reading output text filename
  @retval 2 - error reading nodal output description
  @retval 3 - error reading element output description
  @retval 4 - error reading user defined point output description
  @retval 5 - error reading output graphics format
  @retval 6 - error reading nodal graphics output description
  @retval 7 - error reading element graphics output description
  @retval 8 - error reading output diagram description
  @retval 9 - error reading local coordinate system section

  Created by Tomas Koudelka
*/
long outdriverm::read(XFILE *in)
{
  long i;
  const char *str;

  xfscanf(in, "%k%m", "textout", &flagsw_kwdset, &textout);
  if (textout)
  {
    if (xfscanf(in," %a", outfn) != 1)
    {
      print_err("cannot read filename for output", __FILE__, __LINE__, __func__);
      return 1;
    }
    fprintf(stdout, "\n Output file name (mechanics): %s", outfn);

    if (no.read(in))
      return 2;
    if (no.dstep.st == sel_no)
      str = " not";
    else
      str = "";
    fprintf(stdout, "\n Output of mechanical nodal values will%s be performed", str);

    if (eo.read(in))
      return 3;
    if (eo.dstep.st == sel_no)
      str = " not";
    else
      str = "";
    fprintf(stdout, "\n Output of mechanical element values will%s be performed", str);
    if (po.read(in))
      return 4;
    if (po.dstep.st == sel_no)
      str = " not";
    else
      str = "";
    fprintf(stdout, "\n Output of mechanical UDP values will%s be performed", str);
  }
  if (xfscanf(in, "%k%m", "outgr_format", &graphfmt_kwdset, (int *)&gf) != 2)
  {
    print_err("cannot read type of grahics output", __FILE__, __LINE__, __func__);
    return 5;
  }
  if (gf == grfmt_open_dx)
  {
    Mp->straincomp = 1;
    Mp->strainpos  = 2;
    Mp->strainaver = 1;
    Mp->stresscomp = 1;
    Mp->stresspos  = 2;
    Mp->stressaver = 1;
    Mp->othercomp = 1;
    Mp->otherpos  = 2;
    Mp->otheraver = 1;
    Mp->reactcomp = 1;
  }
  if (gf == grfmt_sep_files_quant){
    xfscanf(in, "%k%ld", "num_gr_files", &nofq);
    if (nofq < 0){
      print_err("wrong number of output result files (%ld < 0)", __FILE__, __LINE__, __func__, nofq);
      abort();
    }
    if (nofq){
      ofq = new outresfilem[nofq];
      //for(i=0; i<nofq; i++)
      //ofq[i].read(in);
    }
  }
  else if (gf != grfmt_no){
    if (xfscanf(in, " %a", outgrfn) != 1){
      print_err("cannot read filename for graphics output", __FILE__, __LINE__, __func__);
      return 5;
    }
  }
  switch (gf)
  {
    case grfmt_no:
      fprintf(stdout, "\n Graphic output will not be performed");
      break;
    case grfmt_open_dx:
      fprintf(stdout, "\n Graphic output will be in OpenDX format");
      fprintf(stdout, "\n Filename of graphic output: %s", outgrfn);
      break;
    case grfmt_femcad:
      fprintf(stdout, "\n Graphic output will be in FemCAD format");
      fprintf(stdout, "\n Filename of graphic output: %s", outgrfn);
      break;
    case grfmt_gid:
      fprintf(stdout, "\n Graphic output will be in GiD format");
      fprintf(stdout, "\n Filename of graphic output: %s", outgrfn);
      break;
    case grfmt_gid_sep:
      fprintf(stdout, "\n Graphic output will be in GiD format separated by quantities");
      fprintf(stdout, "\n Filename of graphic output: %s", outgrfn);
      break;
    case grfmt_vtk:
      fprintf(stdout, "\n Graphic output will be in VTK format");
      fprintf(stdout, "\n Filename of graphic output: %s", outgrfn);
      break;
    case grfmt_gid_vtk:
      fprintf(stdout, "\n Graphic output will be in GiD and VTK format");
      fprintf(stdout, "\n Filename of graphic output : %s", outgrfn);
      break;
    case grfmt_gidsep_vtk:
      fprintf(stdout, "\n Graphic output will be in GiD format separated to quantities and VTK");
      fprintf(stdout, "\n Filename of graphic output : %s", outgrfn);
      break;
    case grfmt_sep_files_quant:
      fprintf(stdout, "\n Output will be in sparated files with selected quantities");
      fprintf(stdout, "\n Number of files: %ld", nofq);
      for(i=0; i<nofq; i++){
        ofq[i].read(in);
        switch (ofq[i].rffmt){
          case resfmt_open_dx:
            fprintf(stdout, "\n  %ld. filename of graphic output in %s format: %s", i+1, "OpenDX", ofq[i].outfbn);
            break;
          case resfmt_gid:
            fprintf(stdout, "\n  %ld. filename of graphic output in %s format: %s", i+1, "GiD", ofq[i].outfbn);
            break;
          case resfmt_vtk:
            fprintf(stdout, "\n  %ld. filename of graphic output in %s format: %s", i+1, "VTK", ofq[i].outfbn);
            break;
          case resfmt_plain_out:
            fprintf(stdout, "\n  %ld. filename of graphic output in %s format: %s", i+1, "PlainOut", ofq[i].outfbn);
            break;
          case resfmt_diag_dat:
            fprintf(stdout, "\n  %ld. filename of graphic output in %s format: %s", i+1, "DiagramDat", ofq[i].outfbn);
            break;
          default:
            print_err("unknown/unsupported file format (%d) is required in %ld. output file %s",
                      __FILE__, __LINE__, __func__, ofq[i].rffmt, i+1, ofq[i].outfbn);
            abort();
        }
      }
      break;
    default:
      print_err("unknown/unsupported file format (%d) is required",
                __FILE__, __LINE__, __func__, gf);
      abort();
  }

  if ((gf != grfmt_no) && (gf != grfmt_sep_files_quant))
  {
    if (nog.read(in))
      return 6;
    if (nog.dstep.st == sel_no)
      str = " not";
    else
      str = "";
    fprintf(stdout, "\n Graphical output of mechanical nodal values will%s be performed", str);
    
    if (eog.read(in))
      return 7;
    if (eog.dstep.st == sel_no)
      str = " not";
    else
      str = "";
    fprintf(stdout, "\n Graphical output of mechanical element values will%s be performed", str);
  }
  switch (Mp->tprob)
  {
    case lin_floating_subdomain:
    case mat_nonlinear_statics:
    case eigen_dynamics:
    case forced_dynamics:
    case mech_timedependent_prob:
    case nonlin_floating_subdomain:
    case growing_mech_structure:
    case earth_pressure:
      if (xfscanf(in, "%k%ld", "numdiag", &ndiag) != 2)
      {
        print_err("cannot read number of diagrams", __FILE__, __LINE__, __func__);
        return 8;
      }
      if (ndiag == 0)
      {
        fprintf (stdout, "\n Number of diagrams : %ld", ndiag);
        return 0;
      }
      if (xfscanf(in, " %a", outdiagfn) != 1)
      {
        print_err("cannot read filename for diagrams", __FILE__, __LINE__, __func__);
        return 8;
      }
      fprintf (stdout, "\n Number of diagrams: %ld", ndiag);
      fprintf (stdout, "\n Diagram filename: %s",outdiagfn);
      odiag = new outdiagm[ndiag];
      outdiagf = new FILE* [ndiag];
      memset(outdiagf, 0, sizeof(*outdiagf)*ndiag);
      for (i=0; i<ndiag; i++)
      {
        if (odiag[i].read(in))
          return 8;
      }      
      break;
    default:
      break;
  }
  if (testlcs() == 1)
  {
    if (xfscanf(in, "%k%ld", "numlcs", &nlcs) != 2)
    {
      print_err("cannot read number of local coordinate systems", __FILE__, __LINE__, __func__);
      return 9;
    }
    if (nlcs <= 0)
    {
      print_err("wrong number of local coordinate systems", __FILE__, __LINE__, __func__);
      return 9;
    }
    lcs = new lcoordsys[nlcs];
    for (i=0; i<nlcs; i++)
    {
      lcs[i].read(in);
    }
  }
  return 0;
}


/**
  Function prints data with output description to the text file given by out.
  
  @param out - pointer to opened text file for output    

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void outdriverm::print(FILE *out)
{
  long i;
  fprintf(out, "\n#\n# outdriver section\n#\n");
  fprintf(out, "%d\n", (int)textout);
  if (textout == on)
  {
    fprintf(out, "%s\n", outfn);
    no.print(out);
    eo.print(out);
    po.print(out);
  }
  fprintf(out, "%d\n", (int)gf);
  if (gf != grfmt_no)
  {
    if (gf == grfmt_sep_files_quant){
      fprintf(out, "%ld", nofq);
      if (nofq){
        for(i=0; i<nofq; i++)
          ofq[i].print(out);
      }
    }
    fprintf(out, "%s\n", outgrfn);
    nog.print(out);
    eog.print(out);
  }
  switch (Mp->tprob)
  {
    case lin_floating_subdomain:
    case mat_nonlinear_statics:
    case eigen_dynamics:
    case forced_dynamics:
    case mech_timedependent_prob:
    case growing_mech_structure:
    case nonlin_floating_subdomain:
    case earth_pressure:
      fprintf(out, "%ld\n", ndiag);
      if (ndiag == 0)
        break;
      fprintf(out, "%s\n", outdiagfn);
      for (i=0; i<ndiag; i++)
        odiag[i].print(out);
      break;
    default:
      break;
  } 
  if (testlcs() == 1)
  {
    fprintf(out, "\n#\n# Definition of local coordinate systems (lcs)\n#\n\n");
    fprintf(out, "%ld # number of local coordinate systems defined", nlcs);
    for (i=0; i<nlcs; i++)
    {
      fprintf(out, "\n# Transformation matrix of %ld. lcs\n", i+1);
      lcs[i].print(out);
    }
  }
}



/**
  Function prints header to the output text file.
  
  @param out - pointer to the opened text output file

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void outdriverm::print_header(FILE *out)
{
  fprintf(out, "%15s ****  *  ****  ****  *\n", " ");
  fprintf(out, "%15s *     *  *     *     *\n", " ");
  fprintf(out, "%15s  *    *  ***   ***   *\n", " ");
  fprintf(out, "%15s   *   *  *     *     *\n", " ");
  fprintf(out, "%15s****   *  *     ****  ****  MEFEL OUTPUT\n", " ");

  fprintf(out, "\n%s\n", Mp->name);
  fprintf(out, "\n\n\n\n\n");
}



/**
  Function prints step number to the output text file.
  
  @param out   - pointer to the opened text output file
  @param lcid  - load case id
  @param istep - integer step id
  @param time  - time or load step

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void outdriverm::print_newstep(FILE *out, long lcid, long istep, double time)
{
  long i;
  double lambda;

  if (textout)
  {
    if(Mp->tpr == seconds)
      lambda = time;
    if(Mp->tpr == minutes)
      lambda = time/60.0; 
    if(Mp->tpr == hours)
      lambda = time/3600.0;
    if(Mp->tpr == days)
      lambda = time/86400.0;
    
    if ((no.sellc.presence_id(lcid+1) && no.dstep.presence_id(istep, time, Mp->timecon)) ||
        (eo.sellc.presence_id(lcid+1) && eo.dstep.presence_id(istep, time, Mp->timecon)) ||
        (po.sellc.presence_id(lcid+1) && po.dstep.presence_id(istep, time, Mp->timecon)))
    { 
      for (i=0; i<53; i++)
        fprintf(out, "*");
      fprintf(out, "\n%10sStep number=% ld, time/load step=% g\n", " ", istep, lambda);
      for (i=0; i<53; i++)
        fprintf(out, "*");
      fprintf(out, "\n\n\n\n");
    }
  }
}


/**
  Function prints required output values to the output text file.
  
  @param out - pointer to the opened text file 
  @param lcid - load case id
  @param istep - step id
  
  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void outdriverm::print_out(FILE *out, long lcid, long istep, double time)
{
  if (textout)
  {
    fprintf(out, "* LOAD CASE NUMBER %ld\n\n\n", lcid+1);
  
    if (no.sellc.presence_id(lcid+1) && (no.dstep.presence_id(istep, time, Mp->timecon)))
      no.print_out(out, lcid);
    if (eo.sellc.presence_id(lcid+1) && (eo.dstep.presence_id(istep, time, Mp->timecon)))
      eo.print_out(out, lcid);
    if (po.sellc.presence_id(lcid+1) && (po.dstep.presence_id(istep, time, Mp->timecon))) 
      po.print_out(out, lcid);
  }
}



/**
  Function prints required output values to the output text file without respect to
  to selected step/time (forced printing).
  
  @param out - pointer to the opened text file 
  @param lcid - load case id
  @param istep - step id
  
  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void outdriverm::print_out_forced(FILE *out, long lcid, long /*istep*/, double /*time*/)
{
  if (textout)
  {
    fprintf(out, "* LOAD CASE NUMBER %ld\n\n\n", lcid+1);
  
    if (no.sellc.presence_id(lcid+1))
      no.print_out(out, lcid);
    if (eo.sellc.presence_id(lcid+1))
      eo.print_out(out, lcid);
    if (po.sellc.presence_id(lcid+1)) 
      po.print_out(out, lcid);
  }
}



/**
  Function prints diagrams.
  
  @param lcid - required load case id
  @param lambda - actual load parameter/actual time
  @param istep - actual step id
  @param fi - array with values of load vector
  
  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void outdriverm::print_diags(long lcid, double lambda, long istep, double *fi)
{
  long i;
  for(i=0; i<ndiag; i++)
    odiag[i].printval(outdiagf[i], lcid, lambda, istep, fi);
}



/**
  Function prints diagrams without respect to
  to selected step/time (forced printing).
  
  @param lcid - required load case id
  @param lambda - actual load parameter/actual time
  @param istep - actual step id
  @param fi - array with values of load vector
  
  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void outdriverm::print_diags_forced(long lcid, double lambda, long istep, double *fi)
{
  long i;
  for(i=0; i<ndiag; i++)
    odiag[i].printval_forced(outdiagf[i], lcid, lambda, istep, fi);
}



/**
  Function prints diagrams.
  
  @param lcid - required load case id
  @param lambda - actual load parameter/actual time
  @param istep - actual step id
  @param fi - array with values of load vector
  @param fr - array of residual %vector components at nodes which should be printed.
  
  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void outdriverm::print_diags(long lcid, double lambda, long istep, double *fi, double *fr)
{
  long i;
  for(i=0; i<ndiag; i++)
    odiag[i].printval(outdiagf[i], lcid, lambda, istep, fi, fr);
}



/**
  Function prints diagrams without respect to
  to selected step/time (forced printing).
  
  @param lcid - required load case id
  @param lambda - actual load parameter/actual time
  @param istep - actual step id
  @param fi - array with values of load vector
  @param fr - array of residual %vector components at nodes which should be printed.
  
  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void outdriverm::print_diags_forced(long lcid, double lambda, long istep, double *fi, double *fr)
{
  long i;
  for(i=0; i<ndiag; i++)
    odiag[i].printval_forced(outdiagf[i], lcid, lambda, istep, fi, fr);
}



/**
  Function prints required value to the graphics files .
  
  @param out - pointer to the opened text file 
  @param lcid - required load case id
  @param lambda - actual load parameter/actual time
  @param istep - actual step id
  @param fi - array with values of load vector
  
  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void outdriverm::print_graphics(FILE *out, long lcid, double lambda, long istep, double *fi)
{
  long i;
  char dlcid[50];

  switch (gf)
  {
    case grfmt_no:
      break;;
    case grfmt_femcad:
      break;
    case grfmt_open_dx:
    {
      for (i=0; i<Mb->nlc; i++)
      {
        if (gf == grfmt_open_dx)
          print_default_dx (Gtm, Mp, Mt, Mm, i, outgrfn);
      }
      break;
    }
    case grfmt_gid:
    {
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1) && (nog.dstep.presence_id(istep, lambda, Mp->timecon)))
        nog.print_graphics(out, lcid, dlcid, gf, fi);
      if (eog.sellc.presence_id(lcid+1) && (eog.dstep.presence_id(istep, lambda, Mp->timecon)))
        eog.print_graphics(out, lcid, dlcid, gf, ide1);
      break;
    }
    case grfmt_gid_sep:
    {      
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1) && (nog.dstep.presence_id(istep, lambda, Mp->timecon)))
        nog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, fi);
      if (eog.sellc.presence_id(lcid+1) && (eog.dstep.presence_id(istep, lambda, Mp->timecon)))
        eog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, ide1);
      break;
    }
    case grfmt_vtk://not implemented for growing structure
    {
      if (nog.sellc.presence_id(lcid+1) && (nog.dstep.presence_id(istep, lambda, Mp->timecon)))//need to include eog separately
      {
        char fname[FNAMELEN+70];
        FILE *vtk_file;
        char *path;
        char *name;
        char *suffix;
        
        filename_decomposition (outgrfn,path,name,suffix);
        sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
        delete [] path; delete [] name; delete [] suffix;
        
        if((vtk_file = fopen(fname, "w"))==NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        switch (Mp->tprob)
        {
          case mat_nonlinear_statics:
          case mech_timedependent_prob:
            if(Mp->tpr == seconds)
              sprintf(dlcid, "%.15e", lambda);
            if(Mp->tpr == minutes)
              sprintf(dlcid, "%.15e", lambda/60.0);
            if(Mp->tpr == hours)
              sprintf(dlcid, "%.15e", lambda/3600.0);
            if(Mp->tpr == days)
              sprintf(dlcid, "%.15e", lambda/86400.0);
            print_vtk_header(vtk_file, istep, (double) atof(dlcid));//print header,points,nodes,elements
            write_vtk_unkn(vtk_file, lcid);//print datapoints in sections
            fclose(vtk_file);
            vtk_num++;
            break;
          default:
          {
            //not implemented yet
          }
        }
      }
      break;
    }
    case grfmt_gid_vtk:
    {
      //
      // output in GiD format
      //
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1) && (nog.dstep.presence_id(istep, lambda, Mp->timecon)))
        nog.print_graphics(out, lcid, dlcid, gf, fi);
      if (eog.sellc.presence_id(lcid+1) && (eog.dstep.presence_id(istep, lambda, Mp->timecon)))
        eog.print_graphics(out, lcid, dlcid, gf, ide1);

      //
      // output in VTK format
      //
      if (nog.sellc.presence_id(lcid+1) && (nog.dstep.presence_id(istep, lambda, Mp->timecon)))//need to include eog separately
      {
        char fname[FNAMELEN+70];
        FILE *vtk_file;
        char *path;
        char *name;
        char *suffix;
        
        filename_decomposition (outgrfn,path,name,suffix);
        sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
        delete [] path; delete [] name; delete [] suffix;
        
        if((vtk_file = fopen(fname, "w"))==NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        switch (Mp->tprob)
        {
          case mat_nonlinear_statics:
          case mech_timedependent_prob:
            print_vtk_header(vtk_file, istep, (double) atof(dlcid)); //print header,points,nodes,elements
            write_vtk_unkn(vtk_file, lcid); //print datapoints in sections
            fclose(vtk_file);
            vtk_num++;
            break;
          default:
          {
            //not implemented yet
          }
        }
      }
      break;
    }
    case grfmt_gidsep_vtk:
    {      
      //
      // output in GiD format
      //
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1) && (nog.dstep.presence_id(istep, lambda, Mp->timecon)))
        nog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, fi);
      if (eog.sellc.presence_id(lcid+1) && (eog.dstep.presence_id(istep, lambda, Mp->timecon)))
        eog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, ide1);

      //
      // output in VTK format
      //
      if (nog.sellc.presence_id(lcid+1) && (nog.dstep.presence_id(istep, lambda, Mp->timecon)))//need to include eog separately
      {
        char fname[FNAMELEN+70];
        FILE *vtk_file;
        char *path;
        char *name;
        char *suffix;
        
        filename_decomposition (outgrfn,path,name,suffix);
        sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
        delete [] path; delete [] name; delete [] suffix;
        
        if((vtk_file = fopen(fname, "w"))==NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        switch (Mp->tprob)
        {
          case mat_nonlinear_statics:
          case mech_timedependent_prob:
            print_vtk_header(vtk_file, istep, (double) atof(dlcid));//print header,points,nodes,elements
            write_vtk_unkn(vtk_file, lcid);//print datapoints in sections
            fclose(vtk_file);
            vtk_num++;
            break;
          default:
          {
            //not implemented yet
          }
        }
      }
      break;
    }
    case grfmt_sep_files_quant:
      print_sep_files_quant(lcid, lambda, istep, Lsrs->give_lhs(lcid), NULL, fi, NULL);
      break;
    default:
      print_err("unknown type of graphics format is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function prints required value to the graphics files without respect to
  to selected step/time (forced printing).
  
  @param out - pointer to the opened text file 
  @param lcid - required load case id
  @param lambda - actual load parameter/actual time
  @param istep - actual step id
  @param fi - array with values of load vector
  
  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void outdriverm::print_graphics_forced(FILE *out, long lcid, double lambda, long istep, double *fi)
{
  long i;
  char dlcid[50];

  switch (gf)
  {
    case grfmt_no:
      break;;
    case grfmt_femcad:
      break;
    case grfmt_open_dx:
    {
      for (i=0; i<Mb->nlc; i++)
      {
        if (gf == grfmt_open_dx)
          print_default_dx (Gtm, Mp, Mt, Mm, i, outgrfn);
      }
      break;
    }
    case grfmt_gid:
    {
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1))
        nog.print_graphics(out, lcid, dlcid, gf, fi);
      if (eog.sellc.presence_id(lcid+1))
        eog.print_graphics(out, lcid, dlcid, gf, ide1);
      break;
    }
    case grfmt_gid_sep:
    {      
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1))
        nog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, fi);
      if (eog.sellc.presence_id(lcid+1))
        eog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, ide1);
      break;
    }
    case grfmt_vtk:
    {
      if (nog.sellc.presence_id(lcid+1))//need to include eog separately
      {
        char fname[FNAMELEN+70];
        FILE *vtk_file;
        char *path;
        char *name;
        char *suffix;
        
        filename_decomposition (outgrfn,path,name,suffix);
        sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
        delete [] path; delete [] name; delete [] suffix;
        
        if((vtk_file = fopen(fname, "w"))==NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        switch (Mp->tprob)
        {
          case mat_nonlinear_statics:
          case mech_timedependent_prob:
            print_vtk_header(vtk_file, istep, (double) atof(dlcid));//print header,points,nodes,elements
            write_vtk_unkn(vtk_file, lcid);//print datapoints in sections
            fclose(vtk_file);
            vtk_num++;
            break;
          default:
          {
            //not implemented yet
          }
        }
      }
      break;
    }
    case grfmt_gid_vtk:
    {
      //
      // forced output in GiD format
      //
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1))
        nog.print_graphics(out, lcid, dlcid, gf, fi);
      if (eog.sellc.presence_id(lcid+1))
        eog.print_graphics(out, lcid, dlcid, gf, ide1);
      //
      // forced output in VTK format
      //
      if (nog.sellc.presence_id(lcid+1))//need to include eog separately
      {
        char fname[FNAMELEN+70];
        FILE *vtk_file;
        char *path;
        char *name;
        char *suffix;
        
        filename_decomposition (outgrfn,path,name,suffix);
        sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
        delete [] path; delete [] name; delete [] suffix;
        
        if((vtk_file = fopen(fname, "w"))==NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        switch (Mp->tprob)
        {
          case mat_nonlinear_statics:
          case mech_timedependent_prob:
            print_vtk_header(vtk_file, istep, (double) atof(dlcid));//print header,points,nodes,elements
            write_vtk_unkn(vtk_file, lcid);//print datapoints in sections
            fclose(vtk_file);
            vtk_num++;
            break;
          default:
          {
            //not implemented yet
          }
        }
      }
      break;
    }
    case grfmt_gidsep_vtk:
    {      
      //
      // forced output in GiD format
      //
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1))
        nog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, fi);
      if (eog.sellc.presence_id(lcid+1))
        eog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, ide1);

      //
      // forced output in VTK format
      //
      if (nog.sellc.presence_id(lcid+1))//need to include eog separately
      {
        char fname[FNAMELEN+70];
        FILE *vtk_file;
        char *path;
        char *name;
        char *suffix;
        
        filename_decomposition (outgrfn,path,name,suffix);
        sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
        delete [] path; delete [] name; delete [] suffix;
        
        if((vtk_file = fopen(fname, "w"))==NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        switch (Mp->tprob)
        {
          case mat_nonlinear_statics:
          case mech_timedependent_prob:
            print_vtk_header(vtk_file, istep, (double) atof(dlcid));//print header,points,nodes,elements
            write_vtk_unkn(vtk_file, lcid);//print datapoints in sections
            fclose(vtk_file);
            vtk_num++;
            break;
          default:
          {
            //not implemented yet
          }
        }
      }
      break;
    }
    default:
      print_err("unknown type of graphics format is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function prints required value to the graphics files .
  
  @param out - pointer to the opened text file 
  @param lcid - required load case id
  @param lambda - actual load parameter/actual time
  @param istep - actual step id
  @param fi - array with values of load vector
  @param fr - array of residual %vector components at nodes which should be printed.
  
  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void outdriverm::print_graphics(FILE *out, long lcid, double lambda, long istep, double *fi, double *fr)
{
  long i;
  char dlcid[50];

  switch (gf)
  {
    case grfmt_no:
      break;;
    case grfmt_femcad:
      break;
    case grfmt_open_dx:
    {
      for (i=0; i<Mb->nlc; i++)
      {
        if (gf == grfmt_open_dx)
          print_default_dx (Gtm, Mp, Mt, Mm, i, outgrfn);
      }
      break;
    }
    case grfmt_gid:
    {
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1) && (nog.dstep.presence_id(istep, lambda, Mp->timecon)))
        nog.print_graphics(out, lcid, dlcid, gf, fi, fr);
      if (eog.sellc.presence_id(lcid+1) && (eog.dstep.presence_id(istep, lambda, Mp->timecon)))
        eog.print_graphics(out, lcid, dlcid, gf, ide1);
      break;
    }
    case grfmt_gid_sep:
    {      
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1) && (nog.dstep.presence_id(istep, lambda, Mp->timecon)))
        nog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, fi, fr);
      if (eog.sellc.presence_id(lcid+1) && (eog.dstep.presence_id(istep, lambda, Mp->timecon)))
        eog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, ide1);
      break;
    }
    case grfmt_vtk://not implemented for growing structure
    {
      if (nog.sellc.presence_id(lcid+1) && (nog.dstep.presence_id(istep, lambda, Mp->timecon)))//need to include eog separately
      {
        char fname[FNAMELEN+70];
        FILE *vtk_file;
        char *path;
        char *name;
        char *suffix;
        
        filename_decomposition (outgrfn,path,name,suffix);
        sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
        delete [] path; delete [] name; delete [] suffix;
        
        if((vtk_file = fopen(fname, "w"))==NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        switch (Mp->tprob)
        {
          case mat_nonlinear_statics:
          case mech_timedependent_prob:
            if(Mp->tpr == seconds)
              sprintf(dlcid, "%.15e", lambda);
            if(Mp->tpr == minutes)
              sprintf(dlcid, "%.15e", lambda/60.0);
            if(Mp->tpr == hours)
              sprintf(dlcid, "%.15e", lambda/3600.0);
            if(Mp->tpr == days)
              sprintf(dlcid, "%.15e", lambda/86400.0);
            print_vtk_header(vtk_file, istep, (double) atof(dlcid));//print header,points,nodes,elements
            write_vtk_unkn(vtk_file, lcid);//print datapoints in sections
            fclose(vtk_file);
            vtk_num++;
            break;
          default:
          {
            //not implemented yet
          }
        }
      }
      break;
    }
    case grfmt_gid_vtk:
    {
      //
      // output in GiD format
      //
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1) && (nog.dstep.presence_id(istep, lambda, Mp->timecon)))
        nog.print_graphics(out, lcid, dlcid, gf, fi, fr);
      if (eog.sellc.presence_id(lcid+1) && (eog.dstep.presence_id(istep, lambda, Mp->timecon)))
        eog.print_graphics(out, lcid, dlcid, gf, ide1);

      //
      // output in VTK format
      //
      if (nog.sellc.presence_id(lcid+1) && (nog.dstep.presence_id(istep, lambda, Mp->timecon)))//need to include eog separately
      {
        char fname[FNAMELEN+70];
        FILE *vtk_file;
        char *path;
        char *name;
        char *suffix;
        
        filename_decomposition (outgrfn,path,name,suffix);
        sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
        delete [] path; delete [] name; delete [] suffix;
        
        if((vtk_file = fopen(fname, "w"))==NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        switch (Mp->tprob)
        {
          case mat_nonlinear_statics:
          case mech_timedependent_prob:
            print_vtk_header(vtk_file, istep, (double) atof(dlcid)); //print header,points,nodes,elements
            write_vtk_unkn(vtk_file, lcid); //print datapoints in sections
            fclose(vtk_file);
            vtk_num++;
            break;
          default:
          {
            //not implemented yet
          }
        }
      }
      break;
    }
    case grfmt_gidsep_vtk:
    {      
      //
      // output in GiD format
      //
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1) && (nog.dstep.presence_id(istep, lambda, Mp->timecon)))
        nog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, fi, fr);
      if (eog.sellc.presence_id(lcid+1) && (eog.dstep.presence_id(istep, lambda, Mp->timecon)))
        eog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, ide1);

      //
      // output in VTK format
      //
      if (nog.sellc.presence_id(lcid+1) && (nog.dstep.presence_id(istep, lambda, Mp->timecon)))//need to include eog separately
      {
        char fname[FNAMELEN+70];
        FILE *vtk_file;
        char *path;
        char *name;
        char *suffix;
        
        filename_decomposition (outgrfn,path,name,suffix);
        sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
        delete [] path; delete [] name; delete [] suffix;
        
        if((vtk_file = fopen(fname, "w"))==NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        switch (Mp->tprob)
        {
          case mat_nonlinear_statics:
          case mech_timedependent_prob:
            print_vtk_header(vtk_file, istep, (double) atof(dlcid));//print header,points,nodes,elements
            write_vtk_unkn(vtk_file, lcid);//print datapoints in sections
            fclose(vtk_file);
            vtk_num++;
            break;
          default:
          {
            //not implemented yet
          }
        }
      }
      break;
    }
    case grfmt_sep_files_quant:
      print_sep_files_quant(lcid, lambda, istep, Lsrs->give_lhs(lcid), NULL, fi, fr);
      break;
    default:
      print_err("unknown type of graphics format is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function prints required value to the graphics files without respect to
  to selected step/time (forced printing).
  
  @param out - pointer to the opened text file 
  @param lcid - required load case id
  @param lambda - actual load parameter/actual time
  @param istep - actual step id
  @param fi - array with values of load vector
  @param fr - array of residual %vector components at nodes which should be printed.
  
  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void outdriverm::print_graphics_forced(FILE *out, long lcid, double lambda, long istep, double *fi, double *fr)
{
  long i;
  char dlcid[50];

  switch (gf)
  {
    case grfmt_no:
      break;;
    case grfmt_femcad:
      break;
    case grfmt_open_dx:
    {
      for (i=0; i<Mb->nlc; i++)
      {
        if (gf == grfmt_open_dx)
          print_default_dx (Gtm, Mp, Mt, Mm, i, outgrfn);
      }
      break;
    }
    case grfmt_gid:
    {
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1))
        nog.print_graphics(out, lcid, dlcid, gf, fi, fr);
      if (eog.sellc.presence_id(lcid+1))
        eog.print_graphics(out, lcid, dlcid, gf, ide1);
      break;
    }
    case grfmt_gid_sep:
    {      
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1))
        nog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, fi, fr);
      if (eog.sellc.presence_id(lcid+1))
        eog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, ide1);
      break;
    }
    case grfmt_vtk:
    {
      if (nog.sellc.presence_id(lcid+1))//need to include eog separately
      {
        char fname[FNAMELEN+70];
        FILE *vtk_file;
        char *path;
        char *name;
        char *suffix;
        
        filename_decomposition (outgrfn,path,name,suffix);
        sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
        delete [] path; delete [] name; delete [] suffix;
        
        if((vtk_file = fopen(fname, "w"))==NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        switch (Mp->tprob)
        {
          case mat_nonlinear_statics:
          case mech_timedependent_prob:
            print_vtk_header(vtk_file, istep, (double) atof(dlcid));//print header,points,nodes,elements
            write_vtk_unkn(vtk_file, lcid);//print datapoints in sections
            fclose(vtk_file);
            vtk_num++;
            break;
          default:
          {
            //not implemented yet
          }
        }
      }
      break;
    }
    case grfmt_gid_vtk:
    {
      //
      // forced output in GiD format
      //
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1))
        nog.print_graphics(out, lcid, dlcid, gf, fi, fr);
      if (eog.sellc.presence_id(lcid+1))
        eog.print_graphics(out, lcid, dlcid, gf, ide1);
      //
      // forced output in VTK format
      //
      if (nog.sellc.presence_id(lcid+1))//need to include eog separately
      {
        char fname[FNAMELEN+70];
        FILE *vtk_file;
        char *path;
        char *name;
        char *suffix;
        
        filename_decomposition (outgrfn,path,name,suffix);
        sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
        delete [] path; delete [] name; delete [] suffix;
        
        if((vtk_file = fopen(fname, "w"))==NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        switch (Mp->tprob)
        {
          case mat_nonlinear_statics:
          case mech_timedependent_prob:
            print_vtk_header(vtk_file, istep, (double) atof(dlcid));//print header,points,nodes,elements
            write_vtk_unkn(vtk_file, lcid);//print datapoints in sections
            fclose(vtk_file);
            vtk_num++;
            break;
          default:
          {
            //not implemented yet
          }
        }
      }
      break;
    }
    case grfmt_gidsep_vtk:
    {      
      //
      // forced output in GiD format
      //
      switch (Mp->tprob)
      {
        case mech_timedependent_prob:
        case growing_mech_structure:
        case eigen_dynamics:
        case forced_dynamics:
          if(Mp->tpr == seconds)
            sprintf(dlcid, "%.15e", lambda);
          if(Mp->tpr == minutes)
            sprintf(dlcid, "%.15e", lambda/60.0); 
          if(Mp->tpr == hours)
            sprintf(dlcid, "%.15e", lambda/3600.0);
          if(Mp->tpr == days)
            sprintf(dlcid, "%.15e", lambda/86400.0);
          break;
        default:
          sprintf(dlcid, "%ld", istep);
      }
      if (nog.sellc.presence_id(lcid+1))
        nog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, fi, fr);
      if (eog.sellc.presence_id(lcid+1))
        eog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, ide1);

      //
      // forced output in VTK format
      //
      if (nog.sellc.presence_id(lcid+1))//need to include eog separately
      {
        char fname[FNAMELEN+70];
        FILE *vtk_file;
        char *path;
        char *name;
        char *suffix;
        
        filename_decomposition (outgrfn,path,name,suffix);
        sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
        delete [] path; delete [] name; delete [] suffix;
        
        if((vtk_file = fopen(fname, "w"))==NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        switch (Mp->tprob)
        {
          case mat_nonlinear_statics:
          case mech_timedependent_prob:
            print_vtk_header(vtk_file, istep, (double) atof(dlcid));//print header,points,nodes,elements
            write_vtk_unkn(vtk_file, lcid);//print datapoints in sections
            fclose(vtk_file);
            vtk_num++;
            break;
          default:
          {
            //not implemented yet
          }
        }
      }
      break;
    }
    default:
      print_err("unknown type of graphics format is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function creates output graphics files and their headers for GiD separated format (grfmt_gidsp).
  For this format (grfmt_gidsp) each required quantity is printed to the sparated file named 
  in following way: filename.{elem|nodal}_{quanitity}{indexOfQuantity}.res
  quatnitity can be: eps, sig, pesp, psig, other.
  indexOfQuantity can be: 
    -# index number for scalars
    -# _v{firstIndex}-{numCompVec}_s{selId} for vectors, where 
      - firstIndex is the index of the first vector component in the given array (eps, sig, other)
      - numCompVec is the number of vector component
      - selId is the index of selection for nodes/elements which the selected vector belongs to
    -# _m[{firstIndex}-{numCompMtx}]_s{selId} for tensors/matrices, where
      - firstIndex is the index of the first tensor/matrix component in the other array 
      - numCompMtx is the number of tensor/matrix component in the other array
      - selId is the index of selection for nodes/elements which the selected tensor/matrix belongs to

  @param mode - opening mode for opened text files

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void outdriverm::create_files_gidsp(const char *mode)
{
  //  char dlcid[1];
  //  dlcid[0] = '\0';
  char *dlcid = NULL;

  if ((gf != grfmt_gid_sep) && (gf != grfmt_gidsep_vtk))
  {
    print_err("wrong graphics format is required", __FILE__, __LINE__, __func__);
    abort();
  }
  if (nog.dstep.st != sel_no)
    nog.print_graphics(outgrfngs, mode, 0, dlcid, gf, NULL);
  if (eog.dstep.st != sel_no)
    eog.print_graphics(outgrfngs, mode, 0, dlcid, gf, ide1);
}

#ifdef _WIN32
char* realpath(const char*, const char*) { return NULL; };
#endif

/**
  Function creates output graphics files and their headers for particular quantitiest (grfmt_quant).
  For this format (grfmt_quant) each required quantity is printed to the specific file but several quantities
  may use the same file name.

  @param mode[in] - opening mode for opened text files
  @param lcid[in] - load case id (used only for VTK format)
  @param istep[in] - time step id (used only for VTK format)
  @param dlcid[in] - string with descriptor of load case (used only for VTK format), if NULL, no VTK file will be created - function
                     is called from print_init but separate VTK files must be generated for particular time steps

  @return The function does not return anything.

  Created by Tomas Koudelka, 08/2021
*/
void outdriverm::create_files_quant(const char *mode, long lcid, long istep, const char *dlcid)
{
  char *rp1 = new char[FILENAME_MAX];
  char *rp2 = new char[FILENAME_MAX];
  char *p1  = new char[FILENAME_MAX];
  char *p2  = new char[FILENAME_MAX];
  ivector idl(nofq);
  list fnlist;
  long i, j, fndiff;
  char *aux = NULL;

  if (gf != grfmt_sep_files_quant){
    print_err("wrong graphics format is required", __FILE__, __LINE__, __func__);
    abort();
  }
  reallocv(nofq, cf);
  // search for identical file names and formats of the graphic files
  for (i=0; i<nofq; i++){
    get_path(ofq[i].outfbn, p1);
    if (realpath(p1, rp1) == NULL){
      print_err("required path'%s' does not exist\n", __FILE__, __LINE__, __func__, p1);
      abort();
    }
    // append format specifier to the path 
    sprintf(p1, "%s%d", rp1, int(ofq[i].rffmt));
    fndiff = 1;
    for(j=0; j<fnlist.count(); j++){
      get_path((char*)fnlist.at(j), p2);      
      if (realpath(p2, rp2) == NULL){
        print_err("required path'%s' does not exist\n", __FILE__, __LINE__, __func__, p2);
        abort();
      }
      // append format specifier to the path 
      sprintf(p2, "%s%d", rp2, int(ofq[i].rffmt));
      if (strcmp(rp1, rp2) == 0){
        fndiff = 0;       // do not create new file
        cf[i] = idl[j];  // store index of file which has identical name and format
        break;
      }        
    }    
    if (fndiff){  // new file will be created
      idl[fnlist.count()] = i; // store index of new file
      cf[i] = -1;  // set negative indicator for the file creation
      aux = new char[FILENAME_MAX];
      strcpy(aux, p1);
      fnlist.append(aux); // store path with format identfier in the list of created files
    }
  }
  for (i=0; i<nofq; i++){
    if ((cf[i] > 0) && (ofq[i].outf == NULL)){
      // open file
      // create file name
      long stochid = -1;
      if(St) stochid = Mp->ns+1;        
      ofq[i].get_modified_fname(stochid, istep);
      ofq[i].outf = fopen(ofq[i].outfmn, mode);
      if (ofq[i].outf == NULL){
        print_err("cannot open output file '%s'", __FILE__, __LINE__, __func__, ofq[i].outfmn);
        abort();
      }
      fseek(ofq[i].outf, 0, SEEK_END); // MS Visual C++ requires that
    }
    // print header if the file is new
    if (ftell(ofq[i].outf) == 0){
      // new file -> only print header
      ofq[i].print_header(istep, Mp->time);
    }
    else
      ofq[i].outf = ofq[-(cf[i]+1)].outf;
  }
  // release memory of the list of file names
  for (i=0; i<fnlist.count(); i++){
    delete [] (char*)fnlist.at(i);
  }
  delete [] rp1;
  delete [] rp2;
  delete [] p1;
  delete [] p2;
  delete [] aux;
}



/**
  The function tests whether the output setting requires transformation of 
  quantities to the user defined local coordinate system.

  @retval 0 - no local coordinate system is required 
  @retval 1 - some local coordinate system is required 

  Created by Tomas Koudelka 27.11.2014
*/
long outdriverm::testlcs()
{
  long i, j;

  if (no.transtra)
  {
    for (i=0; i<no.selnstra.n; i++)
    {
      if (no.transtra[i] > 0)
        return 1;
    }
  }
  if (no.transtre)
  {
    for (i=0; i<no.selnstre.n; i++)
    {
      if (no.transtre[i] > 0)
        return 1;
    }
  }
  if (eo.transtra)
  {
    for (i=0; i<eo.selestra.n; i++)
    {
      if (eo.transtra[i] > 0)
        return 1;
    }
  }
  if (eo.transtre)
  {
    for (i=0; i<eo.selestre.n; i++)
    {
      if (eo.transtre[i] > 0)
        return 1;
    }
  }
  if (po.transtra)
  {
    for (i=0; i<po.npnt; i++)
    {
      if (po.transtra[i] > 0)
        return 1;
    }
  }
  if (po.transtre)
  {
    for (i=0; i<po.npnt; i++)
    {
      if (po.transtre[i] > 0)
        return 1;
    }
  }
  if (nog.transtra)
  {
    for (i=0; i<nog.selnstra.n; i++)
    {
      if (nog.transtra[i] > 0)
        return 1;
    }
  }
  if (nog.transtre)
  {
    for (i=0; i<nog.selnstre.n; i++)
    {
      if (nog.transtre[i] > 0)
        return 1;
    }
  }
  if (eog.transtra)
  {
    for (i=0; i<eog.selestra.n; i++)
    {
      if (eog.transtra[i] > 0)
        return 1;
    }
  }
  if (eog.transtre)
  {
    for (i=0; i<eog.selestre.n; i++)
    {
      if (eog.transtre[i] > 0)
        return 1;
    }
  }
  for(i=0; i<nofq; i++){
    for(j=0; j<ofq[i].nqnt; j++){
      if (ofq[i].qnt[j].lcs){
        return 1;
      }
    }
  }
  return 0;
}



/** 
  Converts selections given by property id to list or range type
  in all output collections according nodal and element properties 
  defined in the topology top.

  @param top - mesh topology with defined properties of nodes and elements

  @return The function does not return anything but it changes internal 
          representation selections of all nodes and elements that were givne by property id.

  Created by Tomas Koudelka, 9.1.2015
*/
void outdriverm::conv_sel_prop(siftop *top)
{
  no.conv_sel_prop(top);
  eo.conv_sel_prop(top);
  nog.conv_sel_prop(top);
  eog.conv_sel_prop(top);
}



/**
  Constructor initializes data to zero values

  Created by Tomas Koudelka
*/
nodeoutm::nodeoutm()
{
  react = 0;
  seldisp = selstra = selstre = seloth = NULL;
  transtra = transtre = NULL; 
}



/**
  Destructor deallocates used memory

  Created by Tomas Koudelka
*/
nodeoutm::~nodeoutm()
{
  delete [] seldisp;
  delete [] selstra;
  delete [] selstre; 
  delete [] seloth;
  delete [] transtra;
  delete [] transtre;
}



/**
  Function reads data with description for output of nodal values for the text output file.

  @param in - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - error reading step or load case selection
  @retval 2 - error reading displacement selection
  @retval 3 - error reading strain selection
  @retval 4 - error reading stress selection
  @retval 5 - error reading other values selection
  @retval 6 - error reading reaction output flag

  Created by Tomas Koudelka
*/
long nodeoutm::read(XFILE *in)
{
  long i;
  // step and loadcases
  xfscanf(in, "%k", "sel_nodstep");
  dstep.read(in);
  if (dstep.st == sel_no)
    return 0;
  xfscanf(in, "%k", "sel_nodlc");
  sellc.read(in);

  // displacements
  xfscanf(in, "%k", "displ_nodes");
  selndisp.read(in);
  switch (selndisp.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
    case sel_prop:
      xfscanf(in, "%k", "displ_comp");
      seldisp = new sel[selndisp.n];
      for (i=0; i<selndisp.n; i++)
        seldisp[i].read(in);
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 2;
  }

  // strains
  xfscanf(in, "%k", "strain_nodes");
  selnstra.read(in);
  switch (selnstra.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
    case sel_prop:
      xfscanf(in, "%k", "nodstrain_comp");
      Mp->straincomp = 1;
      Mp->strainpos  = 2;
      Mp->strainaver = 1;
      transtra = new long[selnstra.n];
      memset(transtra, 0, sizeof(*transtra)*selnstra.n);
      selstra = new sel[selnstra.n];
      for(i=0; i<selnstra.n; i++)
        selstra[i].read(in);
      xfscanf(in, "%k", "nodstra_transfid");
      for(i=0; i<selnstra.n; i++)
      {
        if (xfscanf(in, "%ld", transtra+i) != 1)
        {
          print_err("cannot read strain selection", __FILE__, __LINE__, __func__);
          return 3;
        }
      }
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 3;
  }


  // stresses
  xfscanf(in, "%k", "stress_nodes");
  selnstre.read(in);
  switch (selnstre.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
    case sel_prop:
      xfscanf(in, "%k", "nodstress_comp");
      Mp->straincomp = 1;
      Mp->strainpos  = 2;
      Mp->strainaver = 1;
      Mp->stresscomp = 1;
      Mp->stresspos  = 2;
      Mp->stressaver = 1;
      transtre = new long[selnstre.n];
      memset(transtre, 0, sizeof(*transtre)*selnstre.n);
      selstre = new sel[selnstre.n];
      for(i=0; i<selnstre.n; i++)
        selstre[i].read(in);
      xfscanf(in, "%k", "nodstre_transfid");
      for(i=0; i<selnstre.n; i++)
      {
        if (xfscanf(in, "%ld", transtre+i) != 1)
        {
          print_err("cannot read stress selection", __FILE__, __LINE__, __func__);
          return 3;
        }
      }
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 4;
  }

  // other values
  xfscanf(in, "%k", "other_nodes");
  selnoth.read(in);
  switch (selnoth.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
    case sel_prop:
      xfscanf(in, "%k", "nodother_comp");
      if(Mp->othercomp == 0)
        Mp->othercomp = 1;
      if(Mp->otheraver == 0)
      {
        Mp->otherpos  = 2;
        Mp->otheraver = 1;
      }
      seloth = new sel[selnoth.n];
      for (i=0; i < selnoth.n; i++)
        seloth[i].read(in);
      break; 
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 5;
  }

  // reactions
  if (xfscanf(in, "%k%ld", "reactions", &react) != 2)
  {
    print_err("cannot read selection of reactions", __FILE__, __LINE__, __func__);
    return 6;
  }
  if (react)
    Mp->reactcomp = 1;
  return 0;
}



/**
  Function prints data with description for output of nodal values to the text file.

  @param out - pointer to the opened text file

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutm::print(FILE *out)
{
  long i;

  // step and loadcases
  dstep.print(out);
  if (dstep.st == sel_no)
    return;
  sellc.print(out);
  fprintf(out, "\n");

  // displacements
  selndisp.print(out);
  switch(selndisp.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selndisp.n; i++)
        seldisp[i].print(out);
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }


  // strains
  selnstra.print(out);
  switch (selnstra.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selnstra.n; i++)
        selstra[i].print(out);
      for(i=0; i<selnstra.n; i++)
        fprintf(out, "%ld ", transtra[i]);
      fprintf(out, "\n");
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }

  // stresses
  selnstre.print(out);
  switch (selnstre.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selnstre.n; i++)
        selstre[i].print(out);
      for(i=0; i<selnstre.n; i++)
        fprintf(out, "%ld ", transtre[i]);
      fprintf(out, "\n");
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }

  // other values
  selnoth.print(out);
  switch (selnoth.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selnoth.n; i++)
        seloth[i].print(out);
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }

  // reactions
  fprintf(out, "%ld\n", react);
}



/**
  Function prints required output values for selected nodes and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

  @return The function does not return anything.

  Created by Tomas Koudelka, modified by JK, 17. 4. 2014
*/
void nodeoutm::print_out(FILE *out, long lcid)
{
  if (selndisp.st != sel_no)
  {
    fprintf(out, "** Nodal displacements:\n\n");
    print_disp(out, lcid);
  }
  if ((Mp->homog == 3) || (Mp->homog == 9))
  {
    long j,ncomp = Mm->max_ncompstre;
    double *lhs = Lsrs->give_lhs(0);
    strastre *mstrastre = Mb->give_mstrastre(lcid);
    fprintf(out, "** Calculated macro-strains for the homogenization problem:\n\n");
    j=Ndofm-ncomp;
    for (long i=0; i<ncomp; i++){
      if (mstrastre[i] == stress){
        fprintf(out, "EPS_%ld=% .*e\n", i+1, prstra, lhs[j]);
        j++;
      }
      else
        fprintf(out, "EPS_%ld= --\n", i+1);
    }
    fprintf(out, "\n");
  }

  if (selnstra.st != sel_no)
  {
    fprintf(out, "** Nodal averaged strains:\n\n");
    print_stra(out, lcid);
  }
  if (selnstre.st != sel_no)
  {
    fprintf(out, "** Nodal averaged stresses:\n\n");
    print_stre(out, lcid);
  }
  if (selnoth.st != sel_no)
  {
    fprintf(out, "** Nodal averaged other values:\n\n");
    print_other(out);
  }
  if (react)
  {
    fprintf(out, "** Reactions:\n\n");
    print_react(out);
  }  
}



/**
  Function prints required displacements for selected nodes and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutm::print_disp(FILE *out, long lcid)
{
  long i, j, ndofn;
  double *r;
  
  //  Lahovice
  // vector l,g;
  // matrix tm(ASTCKMAT(3,3));
  // long transf;
  //  konec, Lahovice

  for (i=0; i<Mt->nn; i++)
  {
    if (Gtm->lnso[i]==0)
      continue;
    if (selndisp.presence_id(i))
    {
      fprintf(out, " Node %7ld", i+1);
      //fprintf(out, "  %7ld", i+1);
      //  Lahovice
      // fprintf(out, "%7ld", i+1);
      //  konec, Lahovice
      ndofn = Mt->give_ndofn(i);
      r = new double [ndofn];
      noddispl (lcid, r, i);
      
      /*
      // uprava pro Lahovice, patek 25.6.2010
      fprintf(out, "  % 12.8le % 12.8le % 12.8le",Gtm->gnodes[i].x,Gtm->gnodes[i].y,Gtm->gnodes[i].z);
      
      reallocv (RSTCKVEC(Gtm->gnodes[i].ndofn,l));
      noddispl(lcid, l.a, i);
      
      transf = Mt->nodes[i].transf;
      if (transf>0){
	
	//fprintf (Outm,"\n uzel %ld   %ld",i,Mt->nodes[i].transf);
	
	ndofn = Gtm->gnodes[i].ndofn;
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
      //  konec upravy pro Lahovice
      */


      for (j=0; j<ndofn; j++)
      {
        if (selndisp.presence_id(seldisp,i,j)){
          //  Lahovice
          // fprintf(out, "   %14.11le", r[j]);
          //  konec, Lahovice

          fprintf(out, "   r_%ld=% .*e", j+1, prdisp, r[j]);
          //fprintf(out, "   %15.12le", r[j]);
	}
      }
      delete [] r;

      fprintf(out, "\n");
    } 
  }
  fprintf(out, "\n");
}



/**
  Function prints required strains for selected nodes and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutm::print_stra(FILE *out, long lcid)
{
  long i, j, ncomp, id, ir;
  vector aux;
  matrix t;
  vector str, p;

  if(Outdm->nlcs > 0)
  {
    reallocv(RSTCKVEC(3, p));
    reallocm(RSTCKMAT(3, 3, t));
  }

  for (i=0; i<Mt->nn; i++)
  {
    if (Gtm->lnso[i]==0)
      continue;
    if (selnstra.presence_id(i, ir))
    {
      fprintf(out, " Node %7ld", i+1);
      ncomp = Mt->nodes[i].ncompstr;
      id = lcid*ncomp;
      if (transtra[ir] > 0)
      {
        reallocv(RSTCKVEC(ncomp, aux));
        str.n = ncomp;
        str.a = Mt->nodes[i].strain+id;
        Mt->give_nodal_coord(i, p);
        Outdm->lcs[transtra[ir]-1].give_transfmat(t, p, Mp->time);
        if (ncomp == 4)
          gl_engvectortransf(str, aux, t, planestrain, strain);
        if (ncomp == 6)
          gl_engvectortransf(str, aux, t, spacestress, strain);
        str.a = NULL;
      }
      for (j=0; j<ncomp; j++)
      {
        if (selnstra.presence_id(selstra,i,j))
        {
          if (transtra[ir] > 0)
            fprintf(out, "   leps_%ld=% .*e", j+1, prstra, aux[j]);
          else
            fprintf(out, "   eps_%ld=% .*e", j+1, prstra, Mt->nodes[i].strain[id+j]);
        }
      }
      if (transtra[ir] < 0)
      {
        vector eps(ASTCKVEC(4)), peps(ASTCKVEC(3)); 
        matrix epst(ASTCKMAT(3,3)), pvect(ASTCKMAT(3,3));
        nullv(eps);
        for (j=0; j<ncomp; j++)
          eps[j] = Mt->nodes[i].strain[j];
        if (ncomp == 4)
          vector_tensor (eps, epst, planestrain, strain);
        if (ncomp == 6)
          vector_tensor (eps, epst, spacestress, strain);
        princ_val (epst, peps, pvect, 20, 1.0e-4, Mp->zero, 3,1);
        if (Mt->nodes[i].pstra == NULL)
          Mt->nodes[i].pstra = new double [9];
        for (j=0; j<3; j++)
        {
          Mt->nodes[i].pstra[j] = peps[j];
          fprintf(out, "   peps_%ld=% .*e", j+1, prstra, peps[j]);
        }
        for (j=0; j<3; j++)
        {
          Mt->nodes[i].pstra[3+j] = pvect[0][j];
          fprintf(out, "   a_1%ld=% .*e", j+1, prstra, pvect[0][j]);
        }
        for (j=0; j<3; j++)
        {
          Mt->nodes[i].pstra[6+j] = pvect[2][j];
          fprintf(out, "   a_3%ld=% .*e", j+1, prstra, pvect[2][j]);
        }  
      }      
      fprintf(out, "\n");
    } 
  }
  fprintf(out, "\n");
}



/**
  Function prints required stresses for selected nodes and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutm::print_stre(FILE *out, long lcid)
{
  long i, ir, j, ncomp, id;
  vector aux;
  matrix t;
  vector str, p;

  if(Outdm->nlcs > 0)
  {
    reallocv(RSTCKVEC(3, p));
    reallocm(RSTCKMAT(3, 3, t));
  }

  for (i=0; i<Mt->nn; i++)
  {
    if (Gtm->lnso[i]==0)
      continue;
    if (selnstre.presence_id(i, ir))
    {
      fprintf(out, " Node %7ld", i+1);
      ncomp = Mt->nodes[i].ncompstr;
      id = lcid*ncomp;
      if (transtre[ir] > 0)
      {
        reallocv(RSTCKVEC(ncomp, aux));
        str.n = ncomp;
        str.a = Mt->nodes[i].strain+id;
        Mt->give_nodal_coord(i, p);
        Outdm->lcs[transtre[ir]-1].give_transfmat(t, p, Mp->time);
        if (ncomp == 4)
          gl_engvectortransf(str, aux, t, planestrain, stress);
        if (ncomp == 6)
          gl_engvectortransf(str, aux, t, spacestress, stress);
        str.a = NULL;
      }
      for (j=0; j<ncomp; j++)
      {
        if (selnstre.presence_id(selstre, i, j))
        {
          if (transtre[ir] > 0)
            fprintf(out, "   lsig_%ld=% .*e", j+1, prstre, aux[j]);
          else
            fprintf(out, "   sig_%ld=% .*e", j+1, prstre, Mt->nodes[i].stress[id+j]);
        }
      }
      if (transtre[ir] < 0)
      {
        vector sig(ASTCKVEC(ncomp)), psig(ASTCKVEC(3)); 
        matrix sigt(ASTCKMAT(3,3)), pvect(ASTCKMAT(3,3));
        nullv(sig);
        for (j=0; j<ncomp; j++)
          sig[j] = Mt->nodes[i].stress[j];
        if (ncomp == 4)
          vector_tensor (sig, sigt, planestrain, stress);
        if (ncomp == 6)
          vector_tensor (sig, sigt, spacestress, stress);
        
        princ_val (sigt, psig, pvect, 20, 1.0e-4, Mp->zero, 3,1);
        if (Mt->nodes[i].pstre == NULL)
          Mt->nodes[i].pstre = new double [9];
        for (j=0; j<3; j++)
        {
          Mt->nodes[i].pstre[j] = psig[j];
          fprintf(out, "   psig_%ld=% .*e", j+1, prstre, psig[j]);
        }
        fprintf(out, "   tau_max=% .*e", prstre, (psig[2]-psig[0])/2);
        for (j=0; j<3; j++)
        {
          Mt->nodes[i].pstre[3+j] = pvect[0][j];
          fprintf(out, "   a_1%ld=% .*e", j+1, prstre, pvect[0][j]);
        }
        for (j=0; j<3; j++)
        {
          Mt->nodes[i].pstre[6+j] = pvect[2][j];
          fprintf(out, "   a_3%ld=% .*e", j+1, prstre, pvect[2][j]);
        }  
      }
      fprintf(out, "\n");
   }
 }
  fprintf(out, "\n");
}



/**
  Function prints required other values for selected nodes and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutm::print_other(FILE *out)
{
  long i, j, ncomp;

  for (i=0; i<Mt->nn; i++)
  {
    if (Gtm->lnso[i]==0)
      continue;
    if (selnoth.presence_id(i))
    {
      fprintf(out, " Node %7ld", i+1);
      ncomp = Mt->nodes[i].ncompother;
      for (j=0; j<ncomp; j++)
      {
        if (selnoth.presence_id(seloth, i, j))
          fprintf(out, "   other_%ld=% .*e", j+1, proth, Mt->nodes[i].other[j]);
      }
      fprintf(out, "\n");
    } 
  }
  fprintf(out, "\n");
}



/**
  Function prints all reactions for for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutm::print_react(FILE *out)
{
  long i, j, ndofn;

  for (i=0; i<Mt->nn; i++)
  {
    if (Gtm->lnso[i]==0)
      continue;
    if (Mt->nodes[i].react)
    {
      //fprintf(out, " Node %7ld", i+1);
      fprintf(out, "  %7ld", i+1);
      ndofn = Mt->give_ndofn(i);
      for (j=0; j<ndofn; j++)
        fprintf(out, "   R_%ld=% .*e", j+1, prreac, Mt->nodes[i].r[j]);
      fprintf(out, "\n");
    } 
  }
  fprintf(out, "\n");
}



/**
  Converts selections given by property id to list or range type
  in all output collections according nodal properties defined in the topology top.

  @param top - mesh topology with defined properties of nodes and elements

  @return The function does not return anything but it changes internal 
          representation selections of all nodes that were givne by property id.

  Created by Tomas Koudelka, 9.1.2015
*/
void nodeoutm::conv_sel_prop(siftop *top)
{
  sel *newsel;
  long *newtrans;

  if (selndisp.st == sel_prop)
  {
    if (selndisp.conv_selprop(top, gnod, seldisp, newsel, NULL, newtrans) == 0){
      delete [] seldisp;
      seldisp = newsel;
    }
  }

  if (selnstra.st == sel_prop)
  {
    if (selnstra.conv_selprop(top, gnod, selstra, newsel, transtra, newtrans) == 0){
      delete [] selstra;
      selstra = newsel;
      delete [] transtra;
      transtra = newtrans;
    }
  }

  if (selnstre.st == sel_prop)
  {
    if (selnstre.conv_selprop(top, gnod, selstre, newsel, transtre, newtrans) == 0){
      delete [] selstre;
      selstre = newsel;
      delete [] transtre;
      transtre = newtrans;
    }
  }

  if (selnoth.st == sel_prop)
  {
    if (selnoth.conv_selprop(top, gnod, seloth, newsel, NULL, newtrans) == 0){
      delete [] seloth;
      seloth = newsel;
    }
  }
}








/**
  Constructor initializes data to zero values

  Created by Tomas Koudelka
*/
nodeoutgm::nodeoutgm()
{
  seldisp = selstra = selstre = seloth = selforce = NULL;
  transtra = transtre = NULL; 
}



/**
  Destructor deallocates used memory

  Created by Tomas Koudelka
*/
nodeoutgm::~nodeoutgm()
{
  delete [] seldisp;
  delete [] selstra;
  delete [] selstre; 
  delete [] seloth;
  delete [] selforce;
  delete [] transtra;
  delete [] transtre;
}



/**
  Function reads data with description for output of nodal values from the graphic file.

  @param in - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - error reading displacement selection 
  @retval 2 - error reading strain selection
  @retval 3 - error reading stress selection
  @retval 4 - error reading other values selection
  @retval 5 - error reading force selection

  Created by Tomas Koudelka
*/
long nodeoutgm::read(XFILE *in)
{
  long i;
  // step and loadcases
  xfscanf(in, "%k", "sel_nodstep");
  dstep.read(in);
  if (dstep.st == sel_no)
    return 0;
  xfscanf(in, "%k", "sel_nodlc");
  sellc.read(in);

  // displacements
  xfscanf(in, "%k", "displ_nodes");
  selndisp.read(in);
  switch (selndisp.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
    case sel_prop:
      xfscanf(in, "%k", "displ_comp");
      seldisp = new sel[selndisp.n];
      for (i=0; i<selndisp.n; i++)
        seldisp[i].read(in);
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 1;
  }

  // strains
  xfscanf(in, "%k", "strain_nodes");
  selnstra.read(in);
  switch (selnstra.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
    case sel_prop:
      xfscanf(in, "%k", "nodstrain_comp");
      Mp->straincomp = 1;
      Mp->strainpos  = 2;
      Mp->strainaver = 1;
      transtra = new long[selnstra.n];
      memset(transtra, 0, sizeof(*transtra)*selnstra.n);
      selstra = new sel[selnstra.n];
      for(i=0; i<selnstra.n; i++)
        selstra[i].read(in);
      xfscanf(in, "%k", "nodstra_transfid");
      for(i=0; i<selnstra.n; i++)
      {
        if (xfscanf(in, "%ld", transtra+i) != 1)
        {
          print_err("cannot read strain selection", __FILE__, __LINE__, __func__);
          return 2;
        }
      }
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 2;
  }


  // stresses
  xfscanf(in, "%k", "stress_nodes");
  selnstre.read(in);
  switch (selnstre.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
    case sel_prop:
      xfscanf(in, "%k", "nodstress_comp");
      Mp->straincomp = 1;
      Mp->strainpos  = 2;
      Mp->strainaver = 1;
      Mp->stresscomp = 1;
      Mp->stresspos  = 2;
      Mp->stressaver = 1;
      transtre = new long[selnstre.n];
      memset(transtre, 0, sizeof(*transtre)*selnstre.n);
      selstre = new sel[selnstre.n];
      for(i=0; i<selnstre.n; i++)
        selstre[i].read(in);
      xfscanf(in, "%k", "nodstre_transfid");
      for(i=0; i<selnstre.n; i++)
      {
        if (xfscanf(in, "%ld", transtre+i) != 1)
        {
          print_err("cannot read stress selection", __FILE__, __LINE__, __func__);
          return 3;
        }
      }
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 3;
  }

  // other values
  xfscanf(in, "%k", "other_nodes");
  selnoth.read(in);
  switch (selnoth.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
    case sel_prop:
      xfscanf(in, "%k", "nodother_comp");
      if(Mp->othercomp == 0)
        Mp->othercomp = 1;
      if(Mp->otheraver == 0)
      {
        Mp->otherpos  = 2;
        Mp->otheraver = 1;
      }
      seloth = new sel[selnoth.n];
      for (i=0; i < selnoth.n; i++)
        seloth[i].read(in);
      break; 
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 4;
  }

  // 
  xfscanf(in, "%k", "force_nodes");
  selnforce.read(in);
  switch (selnforce.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
    case sel_prop:
      xfscanf(in, "%k", "force_comp");
      selforce = new sel[selnforce.n];
      Mp->reactcomp = 1;
      for (i=0; i<selndisp.n; i++)
        selforce[i].read(in);
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 5;
  }

  return 0;
}



/**
  Function prints data with description for output of nodal values to the text file.

  @param out - pointer to the opened text file

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print(FILE *out)
{
  long i;

  // step and loadcases
  dstep.print(out);
  if (dstep.st == sel_no)
    return;
  sellc.print(out);
  fprintf(out, "\n");

  // displacements
  selndisp.print(out);
  switch(selndisp.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selndisp.n; i++)
        seldisp[i].print(out);
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }


  // strains
  selnstra.print(out);
  switch (selnstra.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selnstra.n; i++)
        selstra[i].print(out);
      for(i=0; i<selnstra.n; i++)
        fprintf(out, "%ld ", transtra[i]);
      fprintf(out, "\n");
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }

  // stresses
  selnstre.print(out);
  switch (selnstre.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selnstre.n; i++)
        selstre[i].print(out);
      for(i=0; i<selnstre.n; i++)
        fprintf(out, "%ld ", transtre[i]);
      fprintf(out, "\n");
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }

  // other values
  selnoth.print(out);
  switch (selnoth.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selnoth.n; i++)
        seloth[i].print(out);
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }
  selnforce.print(out);
  switch (selnforce.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for (i=0; i<selndisp.n; i++)
        selforce[i].print(out);
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }
}



/**
  Function prints required output values for selected nodes and for given load case and step
  to the output grahics file.
  
  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format
  @param ifor - vector of nodal forces

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_graphics(FILE *out, long lcid, const char *desclcid, graphfmt gf, double *ifor)
{
  if (gf == grfmt_open_dx)
    return;

  if (selndisp.st != sel_no)
  {
    if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
      write_gid_displ(out, lcid, desclcid);
    else
      write_displ(out, lcid, desclcid);
  }
  print_gr_stra_scal(out, lcid, desclcid, gf);
  print_gr_stre_scal(out, lcid, desclcid, gf);
  print_gr_oth_scal(out, lcid, desclcid, gf);

  print_gr_stra_vec(out, lcid, desclcid, gf);
  print_gr_stre_vec(out, lcid, desclcid, gf);
  print_gr_oth_vec(out, lcid, desclcid, gf);

  print_gr_stra_mtx(out, lcid, desclcid, gf);
  print_gr_stre_mtx(out, lcid, desclcid, gf);
  print_gr_oth_mtx(out, lcid, desclcid, gf);

  if (selnforce.st != sel_no)
  {
    if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))    
      write_gid_nforces(out, lcid, desclcid, "Forces", ifor, true);
    else
      write_nforces(out, lcid, desclcid, "Forces", ifor, true);
  }
}



/**
  Function prints required output values for selected nodes and for given load case and step
  to the output grahics file.
  
  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format
  @param ifor - %vector of nodal forces
  @param fr - residual %vector at nodes

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_graphics(FILE *out, long lcid, const char *desclcid, graphfmt gf, double *ifor, double *fr)
{
  if (gf == grfmt_open_dx)
    return;

  if (selndisp.st != sel_no)
  {
    if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
      write_gid_displ(out, lcid, desclcid);
    else
      write_displ(out, lcid, desclcid);
  }
  print_gr_stra_scal(out, lcid, desclcid, gf);
  print_gr_stre_scal(out, lcid, desclcid, gf);
  print_gr_oth_scal(out, lcid, desclcid, gf);

  print_gr_stra_vec(out, lcid, desclcid, gf);
  print_gr_stre_vec(out, lcid, desclcid, gf);
  print_gr_oth_vec(out, lcid, desclcid, gf);

  print_gr_stra_mtx(out, lcid, desclcid, gf);
  print_gr_stre_mtx(out, lcid, desclcid, gf);
  print_gr_oth_mtx(out, lcid, desclcid, gf);

  if (selnforce.st != sel_no)
  {
    if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk)){
      write_gid_nforces(out, lcid, desclcid, "Forces", ifor, true);
      write_gid_nforces(out, lcid, desclcid, "Residual", fr, false);
    }
    else{
      write_nforces(out, lcid, desclcid, "Forces", ifor, true);
      write_nforces(out, lcid, desclcid, "Residual", fr, false);
    }
  }
}



/**
  Function prints required output values for selected nodes and for given load case and step
  to the several output grahics files named by printed quantity component.
  
  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format
  @param ifor - vector of nodal forces

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_graphics(const char *outfn, const char *mode, long lcid, const char *desclcid, graphfmt gf, double *ifor)
{
  char fname[FNAMELEN+70];
  FILE *out;

  if ((gf != grfmt_gid_sep) &&  (gf != grfmt_gidsep_vtk))
  {
    print_err("invalid graphics format is required", __FILE__, __LINE__, __func__);
    return;
  }

  if (selndisp.st != sel_no)
  {
    sprintf(fname, "%s.displ.res", outfn);
    out = fopen(fname, mode);
    if (out == NULL)
    {
      print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
      abort();
    }
    fseek(out, 0, SEEK_END); // MS Visual C++ requires that
    if (ftell(out) == 0)
      fprintf(out, "GiD Post Results File 1.0\n");
    else
    {
      if (desclcid) // regular step output
        write_gid_displ(out, lcid, desclcid);
      else{
          // do nothing because restorage from backup is just performed
      }
    }
    fclose(out);
  }

  print_gr_stra_scal(outfn, mode, lcid, desclcid);
  print_gr_stre_scal(outfn, mode, lcid, desclcid);
  print_gr_oth_scal(outfn, mode, lcid, desclcid);

  print_gr_stra_vec(outfn, mode, lcid, desclcid);
  print_gr_stre_vec(outfn, mode, lcid, desclcid);
  print_gr_oth_vec(outfn, mode, lcid, desclcid);

  print_gr_stra_mtx(outfn, mode, lcid, desclcid);
  print_gr_stre_mtx(outfn, mode, lcid, desclcid);
  print_gr_oth_mtx(outfn, mode, lcid, desclcid);

  if (selnforce.st != sel_no)
  {
    sprintf(fname, "%s.force.res", outfn);
    out = fopen(fname, mode);
    if (out == NULL)
    {
      print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
      abort();
    }
    fseek(out, 0, SEEK_END); // MS Visual C++ requires that
    if (ftell(out) == 0)
      fprintf(out, "GiD Post Results File 1.0\n");
    else{
      if (desclcid) // regular step output
        write_gid_nforces(out, lcid, desclcid, "Forces", ifor, true);
      else{
        // do nothing because restorage from backup is just performed
      }
    }
    fclose(out);
  }
}



/**
  Function prints required output values for selected nodes and for given load case and step
  to the several output grahics files named by printed quantity component.
  
  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format
  @param ifor - vector of nodal forces
  @param fr - residual %vector at nodes

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_graphics(const char *outfn, const char *mode, long lcid, const char *desclcid, graphfmt gf,
                               double *ifor, double *fr)
{
  char fname[FNAMELEN+70];
  FILE *out;

  if ((gf != grfmt_gid_sep) &&  (gf != grfmt_gidsep_vtk))
  {
    print_err("invalid graphics format is required", __FILE__, __LINE__, __func__);
    return;
  }

  if (selndisp.st != sel_no)
  {
    sprintf(fname, "%s.displ.res", outfn);
    out = fopen(fname, mode);
    if (out == NULL)
    {
      print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
      abort();
    }
    fseek(out, 0, SEEK_END); // MS Visual C++ requires that
    if (ftell(out) == 0)
      fprintf(out, "GiD Post Results File 1.0\n");
    else
    {
      if (desclcid) // regular step output
        write_gid_displ(out, lcid, desclcid);
      else{
          // do nothing because restorage from backup is just performed
      }
    }
    fclose(out);
  }

  print_gr_stra_scal(outfn, mode, lcid, desclcid);
  print_gr_stre_scal(outfn, mode, lcid, desclcid);
  print_gr_oth_scal(outfn, mode, lcid, desclcid);

  print_gr_stra_vec(outfn, mode, lcid, desclcid);
  print_gr_stre_vec(outfn, mode, lcid, desclcid);
  print_gr_oth_vec(outfn, mode, lcid, desclcid);

  print_gr_stra_mtx(outfn, mode, lcid, desclcid);
  print_gr_stre_mtx(outfn, mode, lcid, desclcid);
  print_gr_oth_mtx(outfn, mode, lcid, desclcid);

  if (selnforce.st != sel_no)
  {
    sprintf(fname, "%s.force.res", outfn);
    out = fopen(fname, mode);
    if (out == NULL)
    {
      print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
      abort();
    }
    fseek(out, 0, SEEK_END); // MS Visual C++ requires that
    if (ftell(out) == 0)
      fprintf(out, "GiD Post Results File 1.0\n");
    else{
      if (desclcid){
        // regular step output
        write_gid_nforces(out, lcid, desclcid, "Forces", ifor, true);
        write_gid_nforces(out, lcid, desclcid, "Residual", fr, false);
      }
      else{
        // do nothing because restorage from backup is just performed
      }
    }
    fclose(out);
  }
}



/** 
  Function prints values of selected strains as scalars for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_stra_scal(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long i, j, k;

  if (selnstra.st != sel_no)
  {
    for(i=0; i<Mm->max_ncompstrn; i++)
    {
      for (j=0; j<selnstra.n; j++)
      {
        if ((selstra[j].st == sel_mtx) || (selstra[j].st == sel_range_mtx) ||
            (selstra[j].st == sel_range_vec) || (selstra[j].st == sel_vec))
          continue;
        if (selstra[j].presence_id(i))
        {
          if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))    
            write_gid_nodscalar(out, strain, lcid, i, desclcid);
          else
            write_nodscalar(out, strain, lcid, i, desclcid);
          break;
        }
      }
    }
    for (j=0; j<selnstra.n; j++)
    {
      if ((transtra[j] < 0) && (selstra[j].st != sel_mtx) && 
          (selstra[j].st != sel_range_mtx) && (selstra[j].st != sel_range_vec) && 
          (selstra[j].st != sel_vec))
      {
        for(k=0; k<3; k++)
        {
          if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
            write_gid_nodscalar(out, pstrain, lcid, k, desclcid);
          else
            write_nodscalar(out, -k-1, desclcid);
        }
        break;
      }
    }
  }
  return;
}



/**
  Function prints values of selected stresses as scalars for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_stre_scal(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long i, j, k;

  if (selnstre.st != sel_no)
  {
    for(i=0; i<Mm->max_ncompstrn; i++)
    {
      for (j=0; j<selnstre.n; j++)
      {
        if ((selstre[j].st == sel_mtx) || (selstre[j].st == sel_range_mtx) ||
            (selstre[j].st == sel_range_vec) || (selstre[j].st == sel_vec))
          continue;
        if (selstre[j].presence_id(i))
        {
          if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
            write_gid_nodscalar(out, stress, lcid, i, desclcid);
          else
            write_nodscalar(out, stress, lcid, i, desclcid);
          break;
        }
      }
    }
    for (j=0; j<selnstre.n; j++)
    {
      if ((transtre[j] < 0) && (selstre[j].st != sel_mtx) && (selstre[j].st != sel_range_mtx) &&
          (selstre[j].st != sel_range_vec) && (selstre[j].st != sel_vec))
      {
        for(k=0; k<4; k++)
        {
          if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
            write_gid_nodscalar(out, pstress, lcid, k, desclcid);
          else 
            write_nodscalar(out, k, desclcid);
        }
        break;
      }
    }
  }
  return;
}



/** 
  Function prints values of selected other items as scalars for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_oth_scal(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long i, j;

  if (selnoth.st != sel_no)
  {
    for(i=0; i<Mm->max_ncompothern; i++)
    {
      for (j=0; j<selnoth.n; j++)
      {
        if ((seloth[j].st == sel_mtx) || (seloth[j].st == sel_range_mtx) ||
            (seloth[j].st == sel_range_vec) || (seloth[j].st == sel_vec))
          continue;
        if (seloth[j].presence_id(i))
        {
          if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
            write_gid_nodscalar(out, other, lcid, i, desclcid);
          else
            write_nodscalar(out, other, lcid, i, desclcid);
          break;
        }
      }
    }
  }
  return;
}



/** 
  Function prints values of selected strains as vectors for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_stra_vec(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long j, id1, n;
  if (selnstra.st != sel_no)
  {
    for (j=0; j<selnstra.n; j++)
    {
      if (selstra[j].st != sel_range_vec)
        continue;
      id1 = selstra[j].id1[0];
      n  = selstra[j].ncomp[0];
      if (id1+n-1 < Mm->max_ncompstrn)
      {
        if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
          write_gid_nodvector(out, strain, lcid, j, desclcid);
        else
        {
          print_err("required export format of a strain vector is not supported", __FILE__, __LINE__ , __func__);
          abort();
        }
      }
      else
      {
        print_err("required index of strain vector exceeds limits\n"
                  " selected range is id1=%ld, n=%ld, max ncompst=%ld",
                  __FILE__, __LINE__, __func__, id1, n, Mm->max_ncompstrn);
        abort();
      }
    }
  }
  return;
}



/** 
  Function prints values of selected stresses as vectors for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_stre_vec(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long j, id1, n;
  if (selnstre.st != sel_no)
  {
    for (j=0; j<selnstre.n; j++)
    {
      if (selstre[j].st != sel_range_vec)
        continue;
      id1 = selstre[j].id1[0];
      n  = selstre[j].ncomp[0];
      if (id1+n-1 < Mm->max_ncompstrn)
      {
        if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))   
          write_gid_nodvector(out, stress, lcid, j, desclcid);
        else
        {
          print_err("required export format of a stress vector is not supported\n", __FILE__, __LINE__, __func__);
          abort();
        }
      }
      else
      {
        print_err("required index of stress vector exceeds limits\n"
                  " selected range is id1=%ld, n=%ld, max ncompst=%ld\n",
                  __FILE__, __LINE__, __func__,id1, n, Mm->max_ncompstrn);
        abort();
      }
    }
  }
  return;
}



/** 
  Function prints values of selected other items as vectors for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_oth_vec(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long j, id1, n;
  if (selnoth.st != sel_no)
  {
    for (j=0; j<selnoth.n; j++)
    {
      if (seloth[j].st != sel_range_vec)
        continue;
      id1 = seloth[j].id1[0];
      n  = seloth[j].ncomp[0];
      if (id1+n-1 < Mm->max_ncompothern)
      {
        if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
          write_gid_nodvector(out, other, lcid, j, desclcid);
        else
        {
          print_err("required export format of a other values vector is not supported\n", __FILE__, __LINE__, __func__);
          abort();
        }
        break;
      }
      else
      {
        print_err("required index of other values vector exceeds limits\n"
                  " selected range is id1=%ld, n=%ld, max nncompo=%ld\n",
                  __FILE__, __LINE__, __func__, id1, n, Mm->max_ncompothern);
        abort();
      }
    }
  }
  return;
}



/** 
  Function prints values of selected strains as a tensor for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_stra_mtx(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long j;
  if (selnstra.st != sel_no)
  {
    for (j=0; j<selnstra.n; j++)
    {
      if ((selstra[j].st != sel_mtx) && (selstra[j].st != sel_range_mtx))
        continue;
      if (selstra[j].st == sel_range_mtx)
      {
        print_err("sel_range_mtx selection type of a strain tensor is not supported, use sel_mtx(%d)", __FILE__, __LINE__, __func__, sel_mtx);
        abort();
      }
      if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
        write_gid_nodtensor(out, strain, lcid, j, desclcid);
      else
      {
        print_err("required export format of a strain tensor is not supported", __FILE__, __LINE__, __func__);
        abort();
      }
    }
  }
  return;
}



/** 
  Function prints values of selected stresses as a tensor for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_stre_mtx(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long j;
  if (selnstre.st != sel_no)
  {
    for (j=0; j<selnstre.n; j++)
    {
      if ((selstre[j].st != sel_mtx) && (selstre[j].st != sel_range_mtx))
        continue;
      if (selstre[j].st == sel_range_mtx)
      {
        print_err("sel_range_mtx selection type of a stress tensor is not supported, use sel_mtx(%d)",
                  __FILE__, __LINE__, __func__, sel_mtx);
        abort();
      }
      if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
        write_gid_nodtensor(out, stress, lcid, j, desclcid);
      else
      {
        print_err("required export format of a stress tensor is not supported", __FILE__, __LINE__, __func__);
        abort();
      }
    }
  }
  return;
}



/**
  Function prints values of selected other values as a tensor for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_oth_mtx(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long j;
  if (selnoth.st != sel_no)
  {
    for (j=0; j<selnoth.n; j++)
    {
      if ((seloth[j].st != sel_mtx) && (seloth[j].st != sel_range_mtx))
        continue;
      if (seloth[j].st == sel_mtx)
      {
        print_err("sel_mtx selection type of a other values tensor is not supported, use sel_range_mtx(%d)", 
                  __FILE__, __LINE__, __func__, sel_range_mtx);
        abort();
      }
      if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
        write_gid_nodtensor(out, other, lcid, j, desclcid);
      else
      {
        print_err("required export format of an other values is not supported", __FILE__, __LINE__, __func__);
        abort();
      }
    }
  }
  return;
}



//-------------------------------------------------------------------------------------------------------------------

// Functions for output to separated graphic files for nodes

//-------------------------------------------------------------------------------------------------------------------



/** 
  Function prints values of selected strains as scalars for selected load case to the separated 
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_stra_scal(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long i, j, k;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selnstra.st != sel_no)
  {
    for(i=0; i<Mm->max_ncompstrn; i++)
    {
      for (j=0; j<selnstra.n; j++)
      {
        if ((selstra[j].st == sel_mtx) || (selstra[j].st == sel_range_mtx) ||
            (selstra[j].st == sel_range_vec) || (selstra[j].st == sel_vec))
          continue;
        if (selstra[j].presence_id(i))
        {
          sprintf(fname, "%s.nodal_eps%ld.res", outfn, i+1);
          out = fopen(fname, mode);
          if (out == NULL)
          {
            print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
            abort();
          }
          fseek(out, 0, SEEK_END); // MS Visual C++ requires that
          if (ftell(out) == 0)  // new file -> only print header
            fprintf(out, "GiD Post Results File 1.0\n");
          else
          {
            if (desclcid) // regular step output
              write_gid_nodscalar(out, strain, lcid, i, desclcid);
            else{
              // do nothing because restorage from backup is just performed
            }
          }
          fclose(out);
          break;
        }
      }
    }
    for (j=0; j<selnstra.n; j++)
    {
      if ((transtra[j] < 0) && (selstra[j].st != sel_mtx) && 
          (selstra[j].st != sel_range_mtx) && (selstra[j].st != sel_range_vec) &&
          (selstra[j].st != sel_vec))
      {
        for(k=0; k<3; k++)
        {
          sprintf(fname, "%s.nodal_peps%ld.res", outfn, k+1);
          out = fopen(fname, mode);
          if (out == NULL)
          {
            print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
            abort();
          }
          fseek(out, 0, SEEK_END); // MS Visual C++ requires that
          if (ftell(out) == 0)  // new -> only print header
            fprintf(out, "GiD Post Results File 1.0\n");
          else
          {
            if (desclcid) // regular step output
              write_gid_nodscalar(out, pstrain, lcid, k, desclcid);
            else{
              // do nothing because restorage from backup is just performed
            }
          }
          fclose(out);
        }
        break;
      }
    }
  }
  return;
}



/** 
  Function prints values of selected stresses as scalars for selected load case to the
  separated graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_stre_scal(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long i, j, k;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selnstre.st != sel_no)
  {
    for(i=0; i<Mm->max_ncompstrn; i++)
    {
      for (j=0; j<selnstre.n; j++)
      {
        if ((selstre[j].st == sel_mtx) || (selstre[j].st == sel_range_mtx) ||
            (selstre[j].st == sel_range_vec) || (selstre[j].st == sel_vec))
          continue;
        if (selstre[j].presence_id(i))
        {
          sprintf(fname, "%s.nodal_sig%ld.res", outfn, i+1);
          out = fopen(fname, mode);
          if (out == NULL)
          {
            print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
            abort();
          }
          fseek(out, 0, SEEK_END); // MS Visual C++ requires that
          if (ftell(out) == 0)  // new file -> only print header
            fprintf(out, "GiD Post Results File 1.0\n");
          else
          {
            if (desclcid) // regular step output
              write_gid_nodscalar(out, stress, lcid, i, desclcid);
            else{
              // do nothing because restorage from backup is just performed
            }
          }
          fclose(out);
          break;
        }
      }
    }
    for (j=0; j<selnstre.n; j++)
    {
      if ((transtre[j] < 0) && (selstre[j].st != sel_mtx) && (selstre[j].st != sel_range_mtx) &&
          (selstre[j].st != sel_range_vec) && (selstre[j].st != sel_vec))
      {
        for(k=0; k<4; k++)
        {
          sprintf(fname, "%s.nodal_psig%ld.res", outfn, k+1);
          out = fopen(fname, mode);
          if (out == NULL)
          {
            print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
            abort();
          }
          fseek(out, 0, SEEK_END); // MS Visual C++ requires that
          if (ftell(out) == 0)  // new file -> only print header
            fprintf(out, "GiD Post Results File 1.0\n");
          else
          {
            if (desclcid) // regular step output
              write_gid_nodscalar(out, pstress, lcid, k, desclcid);
            else{
              // do nothing because restorage from backup is just performed
            }
          }
          fclose(out);
        }
        break;
      }
    }
  }
  return;
}



/**
  Function prints values of selected other items as scalars for selected load case to separated 
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_oth_scal(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long i, j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selnoth.st != sel_no)
  {
    for(i=0; i<Mm->max_ncompothern; i++)
    {
      for (j=0; j<selnoth.n; j++)
      {
        if ((seloth[j].st == sel_mtx) || (seloth[j].st == sel_range_mtx) ||
            (seloth[j].st == sel_range_vec) || (seloth[j].st == sel_vec))
          continue;
        if (seloth[j].presence_id(i))
        {
          sprintf(fname, "%s.nodal_other%ld.res", outfn, i+1);
          out = fopen(fname, mode);
          if (out == NULL)
          {
            print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
            abort();
          }
          fseek(out, 0, SEEK_END); // MS Visual C++ requires that
          if (ftell(out) == 0)
            fprintf(out, "GiD Post Results File 1.0\n");
          else
          {
            if (desclcid) // regular step output
              write_gid_nodscalar(out, other, lcid, i, desclcid);
            else{
              // do nothing because restorage from backup is just performed
            }
          }
          fclose(out);
          break;
        }
      }
    }
  }
  return;
}



/** 
  Function prints values of selected strains as vectors for selected load case to separated graphics files
  named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_stra_vec(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long j, id1, n;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selnstra.st != sel_no)
  {
    for (j=0; j<selnstra.n; j++)
    {
      if (selstra[j].st != sel_range_vec)
        continue;
      id1 = selstra[j].id1[0];
      n  = selstra[j].ncomp[0];
      if (id1+n-1 < Mm->max_ncompstrn)
      {
        sprintf(fname, "%s.nodal_eps_v%ld-%ld_s%ld.res", outfn, id1+1, n, j+1);
        out = fopen(fname, mode);
        if (out == NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        fseek(out, 0, SEEK_END); // MS Visual C++ requires that
        if (ftell(out) == 0)
          fprintf(out, "GiD Post Results File 1.0\n");
        else
        {
          if (desclcid) // regular step output
            write_gid_nodvector(out, strain, lcid, j, desclcid);
          else{
            // do nothing because restorage from backup is just performed
          }
        }
        fclose(out);
      }
      else
      {
        print_err("required index of strain vector exceeds limits\n"
                  " selected range is id1=%ld, n=%ld, max ncompst=%ld\n",
                  __FILE__, __LINE__, __func__, id1, n, Mm->max_ncompstrn);
        abort();
      }
    }
  }
  return;
}



/**
  Function prints values of selected stresses as vectors for selected load case to separated 
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_stre_vec(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long j, id1, n;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selnstre.st != sel_no)
  {
    for (j=0; j<selnstre.n; j++)
    {
      if (selstre[j].st != sel_range_vec)
        continue;
      id1 = selstre[j].id1[0];
      n  = selstre[j].ncomp[0];
      if (id1+n-1 < Mm->max_ncompstrn)
      {
        sprintf(fname, "%s.nodal_sig_v%ld-%ld_s%ld.res", outfn, id1+1, n, j+1);
        out = fopen(fname, mode);
        if (out == NULL)
        {
          print_err("cannot open graphics file '%s'",__FILE__, __LINE__, __func__, fname);
          abort();
        }
        fseek(out, 0, SEEK_END); // MS Visual C++ requires that
        if (ftell(out) == 0)
          fprintf(out, "GiD Post Results File 1.0\n");
        else
        {
          if (desclcid) // regular step output
            write_gid_nodvector(out, stress, lcid, j, desclcid);
          else{
            // do nothing because restorage from backup is just performed
          }
        }
        fclose(out);
      }
      else
      {
        print_err("required index of stress vector exceeds limits\n"
                  " selected range is id1=%ld, n=%ld, max ncompst=%ld\n",
                  __FILE__, __LINE__, __func__, id1, n, Mm->max_ncompstrn);
        abort();
      }
    }
  }
  return;
}



/**
  Function prints values of selected other items as vectors for selected load case to 
  separated graphics file named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_oth_vec(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long j, id1, n;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selnoth.st != sel_no)
  {
    for (j=0; j<selnoth.n; j++)
    {
      if (seloth[j].st != sel_range_vec)
        continue;
      id1 = seloth[j].id1[0];
      n  = seloth[j].ncomp[0];
      if (id1+n-1 < Mm->max_ncompothern)
      {
        sprintf(fname, "%s.nodal_other_v%ld-%ld_s%ld.res", outfn, id1+1, n, j+1);
        out = fopen(fname, mode);
        if (out == NULL)
        {
          print_err("cannot open graphics file '%s'",__FILE__, __LINE__, __func__, fname);
          abort();
        }
        fseek(out, 0, SEEK_END); // MS Visual C++ requires that
        if (ftell(out) == 0)
          fprintf(out, "GiD Post Results File 1.0\n");
        else
        {
          if (desclcid) // regular step output
            write_gid_nodvector(out, other, lcid, j, desclcid);
          else{
            // do nothing because restorage from backup is just performed
          }
        }
        fclose(out);
      }
      else
      {
        print_err("required index of other values vector exceeds limits\n"
                  " selected range is id1=%ld, n=%ld, max nncompo=%ld\n",
                  __FILE__, __LINE__, __func__, id1, n, Mm->max_ncompothern);
        abort();
      }
    }
  }
  return;
}



/** 
  Function prints values of selected strains as a tensor for selected load case to graphics file.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_stra_mtx(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selnstra.st != sel_no)
  {
    for (j=0; j<selnstra.n; j++)
    {
      if ((selstra[j].st != sel_mtx) && (selstra[j].st != sel_range_mtx))
        continue;
      if (selstra[j].st == sel_range_mtx)
      {
        print_err("sel_range_mtx selection type of a strain tensor is not supported, use sel_mtx(%d)\n",
                  __FILE__, __LINE__, __func__, sel_mtx);
        abort();
      }
      sprintf(fname, "%s.nodal_eps_m_s%ld.res", outfn, j+1);
      out = fopen(fname, mode);
      if (out == NULL)
      {
        print_err("cannot open graphics file '%s'",__FILE__, __LINE__, __func__, fname);
        abort();
      }
      fseek(out, 0, SEEK_END); // MS Visual C++ requires that
      if (ftell(out) == 0)
        fprintf(out, "GiD Post Results File 1.0\n");
      else
      {
        if (desclcid) // regular step output
          write_gid_nodtensor(out, strain, lcid, j, desclcid);
        else{
          // do nothing because restorage from backup is just performed
        }
      }
      fclose(out);
    }
  }
  return;
}



/** 
  Function prints values of selected stresses as a tensor for selected load case to graphics file.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_stre_mtx(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selnstre.st != sel_no)
  {
    for (j=0; j<selnstre.n; j++)
    {
      if ((selstre[j].st != sel_mtx) && (selstre[j].st != sel_range_mtx))
        continue;
      if (selstre[j].st == sel_range_mtx)
      {
        print_err("sel_range_mtx selection type of a stress tensor is not supported, use sel_mtx(%d)\n",
                  __FILE__, __LINE__, __func__, sel_mtx);
        abort();
      }
      sprintf(fname, "%s.nodal_sig_m_s%ld.res", outfn, j+1);
      out = fopen(fname, mode);
      if (out == NULL)
      {
        print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
        abort();
      }
      fseek(out, 0, SEEK_END); // MS Visual C++ requires that
      if (ftell(out) == 0)
        fprintf(out, "GiD Post Results File 1.0\n");
      else
      {
        if (desclcid) // regular step output
          write_gid_nodtensor(out, stress, lcid, j, desclcid);
        else{
          // do nothing because restorage from backup is just performed
        }
      }
      fclose(out);
    }
  }
  return;
}



/**
  Function prints values of selected other values as a tensor for selected load case to graphics file.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void nodeoutgm::print_gr_oth_mtx(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selnoth.st != sel_no)
  {
    for (j=0; j<selnoth.n; j++)
    {
      if ((seloth[j].st != sel_mtx) && (seloth[j].st != sel_range_mtx))
        continue;
      if (seloth[j].st == sel_mtx)
      {
        print_err("sel_mtx selection type of an other values tensor is not supported, use sel_range_mtx(%d)\n", 
                  __FILE__, __LINE__, __func__, sel_range_mtx);
        abort();
      }
      sprintf(fname, "%s.nodal_other_m%ld-%ld_s%ld.res", outfn, seloth[j].id1[0]+1, seloth[j].ncomp[0], j+1);
      out = fopen(fname, mode);
      if (out == NULL)
      {
        print_err("cannot open graphics file '%s'",__FILE__, __LINE__, __func__, fname);
        abort();
      }
      fseek(out, 0, SEEK_END); // MS Visual C++ requires that
      if (ftell(out) == 0)
        fprintf(out, "GiD Post Results File 1.0\n");
      else
      {
        if (desclcid) // regular step output
          write_gid_nodtensor(out, other, lcid,j, desclcid);
        else{
          // do nothing because restorage from backup is just performed
        }
      }
      fclose(out);
    }
  }
  return;
}



/**
  Converts selections given by property id to list or range type
  in all output collections according nodal properties defined in the topology top.

  @param top - mesh topology with defined properties of nodes and elements

  @return The function does not return anything but it changes internal 
          representation selections of all nodes that were givne by property id.

  Created by Tomas Koudelka, 9.1.2015
*/
void nodeoutgm::conv_sel_prop(siftop *top)
{
  sel *newsel;
  long *newtrans;

  if (selndisp.st == sel_prop)
  {
    if (selndisp.conv_selprop(top, gnod, seldisp, newsel, NULL, newtrans) == 0){
      delete [] seldisp;
      seldisp = newsel;
    }
  }

  if (selnstra.st == sel_prop)
  {
    if (selnstra.conv_selprop(top, gnod, selstra, newsel, transtra, newtrans) == 0){
      delete [] selstra;
      selstra = newsel;
      delete [] transtra;
      transtra = newtrans;
    }
  }

  if (selnstre.st == sel_prop)
  {
    if (selnstre.conv_selprop(top, gnod, selstre, newsel, transtre, newtrans) == 0){
      delete [] selstre;
      selstre = newsel;
      delete [] transtre;
      transtre = newtrans;
    }
  }

  if (selnoth.st == sel_prop)
  {
    if (selnoth.conv_selprop(top, gnod, seloth, newsel, NULL, newtrans) == 0){
      delete [] seloth;
      seloth = newsel;
    }
  }

  if (selnforce.st == sel_prop)
  {
    if (selnforce.conv_selprop(top, gnod, selforce, newsel, NULL, newtrans) == 0){
      delete [] selforce;
      selforce = newsel;
    }
  }
}



/**
  Constructor initializes data to zero values

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
elemoutm::elemoutm()
{
  selstra = selstre = seloth = NULL;
  transtra = transtre = NULL; 
}



/**
  Destructor deallocates used memory

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
elemoutm::~elemoutm()
{
  delete [] selstra;
  delete [] selstre;
  delete [] seloth;
  delete [] transtra;
  delete [] transtre;
}



/**
  Function reads data with description for output of element values from the text file.

  @param in - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - error reading step 
  @retval 2 - error reading strain selection
  @retval 3 - error reading stress selection
  @retval 4 - error reading other values selection

  Created by Tomas Koudelka
*/
long elemoutm::read(XFILE *in)
{
  long i;
  // step and loadcases
  xfscanf(in, "%k", "sel_elemstep");
  dstep.read(in);
  if (dstep.st == sel_no)
    return 0;
  xfscanf(in, "%k", "sel_elemlc");
  sellc.read(in);
  if (sellc.st == sel_no)
    return 0;


  // strains
  xfscanf(in, "%k", "strain_elems");
  selestra.read(in);
  switch (selestra.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
    case sel_prop:
      xfscanf(in, "%k", "elemstrain_comp");
      Mp->straincomp = 1;
      if (Mp->strainpos == 0)
	Mp->strainpos = 1;
      transtra = new long[selestra.n];
      memset(transtra, 0, sizeof(*transtra)*selestra.n);
      selstra = new sel[selestra.n];
      for(i=0; i<selestra.n; i++)
        selstra[i].read(in);
      xfscanf(in, "%k", "elemstra_transfid");
      for(i=0; i<selestra.n; i++)
      {
        if (xfscanf(in, "%ld", transtra+i) != 1)
        {
          print_err("cannot read strain selection", __FILE__, __LINE__, __func__);
          return 2;
        }
      }
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 2;
  }


  // stresses
  xfscanf(in, "%k", "stress_elems");
  selestre.read(in);
  switch (selestre.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
    case sel_prop:
      xfscanf(in, "%k", "elemstress_comp");
      Mp->straincomp = 1;
      if (Mp->strainpos == 0)
	Mp->strainpos = 1;
      Mp->stresscomp = 1;
      if (Mp->stresspos == 0)
	Mp->stresspos = 1;
      transtre = new long[selestre.n];
      memset(transtre, 0, sizeof(*transtre)*selestre.n);
      selstre = new sel[selestre.n];
      for(i=0; i<selestre.n; i++)
        selstre[i].read(in);
      xfscanf(in, "%k", "elemstre_transfid");
      for(i=0; i<selestre.n; i++)
      {
        if (xfscanf(in, "%ld", transtre+i) != 1)
        {
          print_err("cannot read stress selection", __FILE__, __LINE__, __func__);
          return 3;
        }
      }
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 3;
  }

  // other values
  xfscanf(in, "%k", "other_elems");
  seleoth.read(in);
  switch (seleoth.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
    case sel_prop:
      xfscanf(in, "%k", "elemother_comp");
      Mp->othercomp = 1;
      if (Mp->otherpos == 0)
	Mp->otherpos = 1; 
      seloth = new sel[seleoth.n]; 
      for (i=0; i<seleoth.n; i++)
        seloth[i].read(in);
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 4;
  }
  return 0;
}



/**
  Function prints data with description for output of element values to the text file.

  @param out - pointer to the opened text file

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutm::print(FILE *out)
{
  long i;

  // step and loadcases
  dstep.print(out);
  if (dstep.st == sel_no)
    return;
  sellc.print(out);
  if (sellc.st == sel_no)
    return;

  // strains
  selestra.print(out);
  switch (selestra.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selestra.n; i++)
        selstra[i].print(out);
      for(i=0; i<selestra.n; i++)
        fprintf(out, "%ld ", transtra[i]);
      fprintf(out, "\n");
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }

  // stresses
  selestre.print(out);
  switch (selestre.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selestre.n; i++)
        selstre[i].print(out);
      for(i=0; i<selestre.n; i++)
        fprintf(out, "%ld ", transtre[i]);
      fprintf(out, "\n");
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }

  // other values
  seleoth.print(out);
  switch (seleoth.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<seleoth.n; i++)
        seloth[i].print(out);
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }
}



/**
  Function prints required output values for selected elements and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutm::print_out(FILE *out, long lcid)
{
  if (selestra.st != sel_no)
  {
    fprintf(out, "** Element strains:\n\n");
    print_stra(out, lcid);
  }
  if (selestre.st != sel_no)
  {
    fprintf(out, "** Element stresses:\n\n");
    print_stre(out, lcid);
  }
  if (seleoth.st != sel_no)
  {
    fprintf(out, "** Element other values:\n\n");
    print_other(out);
  }
}



/**
  Function prints required strains for selected elements and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutm::print_stra(FILE *out, long lcid)
{
  long i, j, k, ncomp, id, ipp, tnipe, ir;
  vector aux;
  matrix t;
  vector str, p;

  if(Outdm->nlcs > 0)
  {
    reallocv(RSTCKVEC(3, p));
    reallocm(RSTCKMAT(3, 3, t));
  }

  for (i=0; i<Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;
    if (selestra.presence_id(i, ir))
    {
      ipp = Mt->elements[i].ipp[0][0];
      tnipe = Mt->give_totnip(i);
      fprintf(out, "   Element %7ld, integration points %ld - %ld\n  ", i+1, ipp+1, ipp+tnipe);
      for (j=0; j<tnipe; j++)
      {
	ncomp = Mm->ip[ipp+j].ncompstr;
	id = lcid*ncomp;
        if (transtra[ir] > 0)
        {
          reallocv(RSTCKVEC(ncomp, aux));
          str.n = ncomp;
          str.a = Mm->ip[ipp+j].strain+id;
          ipcoord(Mm->elip[ipp+j], ipp+j, 0, 0, p);
          Outdm->lcs[transtra[ir]-1].give_transfmat(t, p, Mp->time);
          gl_engvectortransf(str, aux, t, planestrain, strain);
          str.a = NULL;
        }
        for (k=0; k<ncomp; k++)
        {
          if (selestra.presence_id(selstra, i, k))
          {
            if (transtra[ir] > 0)
              fprintf(out, "   leps_%ld=% .*e", k+1, prstra, aux[k]);
            else
              fprintf(out, "   eps_%ld=% .*e", k+1, prstra, Mm->ip[ipp+j].strain[id+k]);
          }
        }
        if (j < tnipe-1)
          fprintf(out, "\n  ");
        else
          fprintf(out, "\n");
      }
    } 
  }
  fprintf(out, "\n");
}



/**
  Function prints required stresses for selected elements and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

  Created by Tomas Koudelka
*/
void elemoutm::print_stre(FILE *out, long lcid)
{
  long i, j, k, ncomp, ipp, tnipe, id, ir;
  vector aux;
  matrix t;
  vector str, p;
  
  if(Outdm->nlcs > 0)
  {
    reallocv(RSTCKVEC(3, p));
    reallocm(RSTCKMAT(3, 3, t));
  }


  for (i=0; i<Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;
    if (selestre.presence_id(i, ir))
    {
      ipp = Mt->elements[i].ipp[0][0];
      tnipe = Mt->give_totnip(i);
      fprintf(out, "   Element %7ld, integration points %ld - %ld\n  ", i+1, ipp+1, ipp+tnipe);
      for (j=0; j<tnipe; j++)
      {
	ncomp = Mm->ip[ipp+j].ncompstr;
	id = lcid*ncomp;
        if (transtre[ir] > 0)
        {
          reallocv(RSTCKVEC(ncomp, aux));
          str.n = ncomp;
          str.a = Mm->ip[ipp+j].stress+id;
          ipcoord(Mm->elip[ipp+j], ipp+j, 0, 0, p);
          Outdm->lcs[transtre[ir]-1].give_transfmat(t, p, Mp->time);
          gl_engvectortransf(str, aux, t, planestrain, strain);
          str.a = NULL;
        }
        for (k=0; k<ncomp; k++)
        {
          if (selestre.presence_id(selstre, i, k))
          {
            if (transtre[ir] > 0)
              fprintf(out, "   lsig_%ld=% .*e", k+1, prstre, aux[k]);
            else
              //              fprintf(out, "   sig_%ld=% .*e", k+1, prstre, Mm->ip[ipp+j].stress[id+k]);
              fprintf(out, "   sig_%ld=% .*e", k+1, prstre, Mm->ip[ipp+j].stress[k]);
          }
        }
        if (j < tnipe-1)
          fprintf(out, "\n  ");
        else
          fprintf(out, "\n");
      }
    } 
  }
  fprintf(out, "\n");
}



/**
  Function prints required other values for selected elements and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutm::print_other(FILE *out)
{
  long i, j, k, ncomp, ipp, tnipe;

  for (i=0; i<Mt->ne; i++)
  {
    if (Gtm->leso[i]==0)
      continue;
    if (seleoth.presence_id(i))
    {
      ipp = Mt->elements[i].ipp[0][0];
      tnipe = Mt->give_totnip(i);
      fprintf(out, "   Element %7ld, integration points %ld - %ld\n  ", i+1, ipp+1, ipp+tnipe);
      ncomp = Mm->ip[ipp].ncompeqother;
      for (j=0; j<tnipe; j++)
      {
        for (k=0; k<ncomp; k++)
        {
          if (seleoth.presence_id(seloth, i, k))
            fprintf(out, "   other_%ld=% .*e", k+1, proth, Mm->ip[ipp+j].eqother[k]);
        }
        if (j < tnipe-1)
          fprintf(out, "\n  ");
        else
          fprintf(out, "\n");
      }
    } 
  }
  fprintf(out, "\n");
}



/**
  Converts selections given by property id to list or range type
  in all output collections according element properties defined in the topology top.

  @param top - mesh topology with defined properties of elements

  @return The function does not return anything but it changes internal 
          representation selections of all elements that were givne by property id.

  Created by Tomas Koudelka, 9.1.2015
*/
void elemoutm::conv_sel_prop(siftop *top)
{
  sel *newsel;
  long *newtrans;

  if (selestra.st == sel_prop)
  {
    if (selestra.conv_selprop(top, gelem, selstra, newsel, transtra, newtrans) == 0){
      delete [] selstra;
      selstra = newsel;
      delete [] transtra;
      transtra = newtrans;
    }
  }

  if (selestre.st == sel_prop)
  {
    if (selestre.conv_selprop(top, gelem, selstre, newsel, transtre, newtrans) == 0){
      delete [] selstre;
      selstre = newsel;
      delete [] transtre;
      transtre = newtrans;
    }
  }

  if (seleoth.st == sel_prop)
  {
    if (seleoth.conv_selprop(top, gelem, seloth, newsel, NULL, newtrans) == 0){
      delete [] seloth;
      seloth = newsel;
    }
  }
}





/**
  Constructor initializes data to zero values

  Created by Tomas Koudelka
*/
elemoutgm::elemoutgm()
{
  selstra = selstre = seloth = NULL;
  transtra = transtre = NULL; 
  ide1 = 1;
}



/**
  Destructor deallocates used memory

  Created by Tomas Koudelka
*/
elemoutgm::~elemoutgm()
{
  delete [] selstra;
  delete [] selstre;
  delete [] seloth;
  delete [] transtra;
  delete [] transtre;
}



/**
  Function reads data with description for output of element values from the text file.

  @param in - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - error reading step 
  @retval 2 - error reading strain selection
  @retval 3 - error reading stress selection
  @retval 4 - error reading other values selection

  Created by Tomas Koudelka
*/
long elemoutgm::read(XFILE *in)
{
  long i;
  // step and loadcases
  xfscanf(in, "%k", "sel_elemstep");
  dstep.read(in);
  if (dstep.st == sel_no)
    return 0;
  xfscanf(in, "%k", "sel_elemlc");
  sellc.read(in);
  if (sellc.st == sel_no)
    return 0;


  // strains
  xfscanf(in, "%k", "strain_elems");
  selestra.read(in);
  switch (selestra.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
    case sel_prop:
      xfscanf(in, "%k", "elemstrain_comp");
      Mp->straincomp = 1;
      if (Mp->strainpos == 0)
        Mp->strainpos = 1;
      transtra = new long[selestra.n];
      memset(transtra, 0, sizeof(*transtra)*selestra.n);
      selstra = new sel[selestra.n];
      for(i=0; i<selestra.n; i++)
        selstra[i].read(in);
      xfscanf(in, "%k", "elemstra_transfid");
      for(i=0; i<selestra.n; i++)
      {
        if (xfscanf(in, "%ld", transtra+i) != 1)
        {
          print_err("cannot read strain selection", __FILE__, __LINE__, __func__);
          return 2;
        }
      }
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 2;
  }


  // stresses
  xfscanf(in, "%k", "stress_elems");
  selestre.read(in);
  switch (selestre.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
    case sel_prop:
      xfscanf(in, "%k", "elemstress_comp");
      Mp->straincomp = 1;
      if (Mp->strainpos == 0)
        Mp->strainpos = 1;
      Mp->stresscomp = 1;
      if (Mp->stresspos == 0)
	Mp->stresspos = 1;
      transtre = new long[selestre.n];
      memset(transtre, 0, sizeof(*transtre)*selestre.n);
      selstre = new sel[selestre.n];
      for(i=0; i<selestre.n; i++)
        selstre[i].read(in);
      xfscanf(in, "%k", "elemstre_transfid");
      for(i=0; i<selestre.n; i++)
      {
        if (xfscanf(in, "%ld", transtre+i) != 1)
        {
          print_err("cannot read stress selection", __FILE__, __LINE__, __func__);
          return 3;
        }
      }
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 3;
  }

  // other values
  xfscanf(in, "%k", "other_elems");
  seleoth.read(in);
  switch (seleoth.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
    case sel_prop:
      xfscanf(in, "%k", "elemother_comp");
      Mp->othercomp = 1;
      if (Mp->otherpos == 0)
	Mp->otherpos = 1;
      seloth = new sel[seleoth.n]; 
      for (i=0; i<seleoth.n; i++)
        seloth[i].read(in);
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 4;
  }
  return 0;
}



/**
  Function prints data with description for output of element values to the text file.

  @param out - pointer to the opened text file

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print(FILE *out)
{
  long i;

  // step and loadcases
  dstep.print(out);
  if (dstep.st == sel_no)
    return;
  sellc.print(out);
  if (sellc.st == sel_no)
    return;

  // strains
  selestra.print(out);
  switch (selestra.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selestra.n; i++)
        selstra[i].print(out);
      for(i=0; i<selestra.n; i++)
        fprintf(out, "%ld ", transtra[i]);
      fprintf(out, "\n");
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }

  // stresses
  selestre.print(out);
  switch (selestre.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selestre.n; i++)
        selstre[i].print(out);
      for(i=0; i<selestre.n; i++)
        fprintf(out, "%ld ", transtre[i]);
      fprintf(out, "\n");
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }

  // other values
  seleoth.print(out);
  switch (seleoth.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<seleoth.n; i++)
        seloth[i].print(out);
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }
}



/**  
  Function prints required output values for selected elements and for given load case and step
  to the output grahics file.
  
  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format
  @param idelem1 - number of the first element for GiD output (normally should be 1), 
                   set by print_graphics function

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_graphics(FILE *out, long lcid, const char *desclcid, graphfmt gf, long idelem1)
{
  ide1 = idelem1;
  if (gf == grfmt_open_dx)
    return;
    
  print_gr_stra_scal(out, lcid, desclcid, gf);
  print_gr_stre_scal(out, lcid, desclcid, gf);
  print_gr_oth_scal(out, lcid, desclcid, gf);

  print_gr_stra_vec(out, lcid, desclcid, gf);
  print_gr_stre_vec(out, lcid, desclcid, gf);
  print_gr_oth_vec(out, lcid, desclcid, gf);

  print_gr_stra_mtx(out, lcid, desclcid, gf);
  print_gr_stre_mtx(out, lcid, desclcid, gf);
  print_gr_oth_mtx(out, lcid, desclcid, gf);

}



/** 
  Function prints values of selected strains as scalars for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_stra_scal(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long i, j;

  if (selestra.st != sel_no)
  {
    for(i=0; i<Mm->max_ncompstre; i++)
    {
      for (j=0; j<selestra.n; j++)
      {
        if ((selstra[j].st == sel_mtx) || (selstra[j].st == sel_range_mtx) ||
            (selstra[j].st == sel_range_vec) || (selstra[j].st == sel_vec))
          continue;
        if (selstra[j].presence_id(i))
        {
          if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
            write_gid_elemscalar(out, strain, lcid, i, desclcid);
          if (gf == grfmt_femcad)
            write_elemscalar(out, strain, lcid, i, desclcid);
          break;
        }
      }
    }
  }
}



/** 
  Function prints values of selected stresses as scalars for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_stre_scal(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long i, j;

  if (selestre.st != sel_no)
  {
    for(i=0; i<Mm->max_ncompstre; i++)
    {
      for (j=0; j<selestre.n; j++)
      {
        if ((selstre[j].st == sel_mtx) || (selstre[j].st == sel_range_mtx) ||
            (selstre[j].st == sel_range_vec) || (selstre[j].st == sel_vec))
          continue;
        if (selstre[j].presence_id(i))
        {
          if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
            write_gid_elemscalar(out, stress, lcid, i, desclcid);
          if (gf == grfmt_femcad)
            write_elemscalar(out, stress, lcid, i, desclcid);
          break;
        }
      }
    }
  }
}



/** 
  Function prints values of selected other items as scalars for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_oth_scal(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long i, j;

  if (seleoth.st != sel_no)
  {
    for(i=0; i<Mm->max_ncompothere; i++)
    {
      for (j=0; j<seleoth.n; j++)
      {
        if ((seloth[j].st == sel_mtx) || (seloth[j].st == sel_range_mtx) ||
            (seloth[j].st == sel_range_vec) || (seloth[j].st == sel_vec))
          continue;
        if (seloth[j].presence_id(i))
        {
          if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
            write_gid_elemscalar(out, other, lcid, i, desclcid);
          if (gf == grfmt_femcad)
            write_elemscalar(out, other, lcid, i, desclcid);
          break;
        }
      }
    }
  }
}



/** 
  Function prints values of selected strains as vectors for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_stra_vec(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long j, id1, n;

  if (selestra.st != sel_no)
  {
    for (j=0; j<selestra.n; j++)
    {
      if (selstra[j].st != sel_range_vec)
        continue;
      id1 = selstra[j].id1[0];
      n  = selstra[j].ncomp[0];
      if (id1+n-1 < Mm->max_ncompstre)
      {
        if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
          write_gid_elemvector(out, strain, lcid, j, desclcid);
        else
        {
          print_err("required export format of a strain vector is not supported", __FILE__, __LINE__, __func__);
          abort();
        }
      }
      else
      {
        print_err("required index of strain vector exceeds limits\n"
                  " selected range is id1=%ld, n=%ld, max ncompst=%ld",
                  __FILE__, __LINE__, __func__, id1, n, Mm->max_ncompstrn);
        abort();
      }
    }
  }
}



/** 
  Function prints values of selected stresses as vectors for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_stre_vec(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long j, id1, n;

  if (selestre.st != sel_no)
  {
    for (j=0; j<selestre.n; j++)
    {
      if (selstre[j].st != sel_range_vec)
        continue;
      id1 = selstre[j].id1[0];
      n  = selstre[j].ncomp[0];
      if (id1+n-1 < Mm->max_ncompstre)
      {
        if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
          write_gid_elemvector(out, stress, lcid, j, desclcid);
        else
        {
          print_err("required export format of a stress vector is not supported", __FILE__, __LINE__, __func__);
          abort();
        }
      }
      else
      {
        print_err("required index of stress vector exceeds limits\n"
                  " selected range is id1=%ld, n=%ld, max ncompst=%ld",
                  __FILE__, __LINE__, __func__, id1, n, Mm->max_ncompstrn);
        abort();
      }
    }
  }
}



/** 
  Function prints values of selected other items as vectors for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_oth_vec(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long j;
  long id1, n;

  if (seleoth.st != sel_no)
  {
    for (j=0; j<seleoth.n; j++)
    {
      if (seloth[j].st != sel_range_vec)
        continue;
      id1 = seloth[j].id1[0];
      n  = seloth[j].ncomp[0];
      if (id1+n-1 < Mm->max_ncompothere)
      {
        if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
          write_gid_elemvector(out, other, lcid, j, desclcid);
        else
        {
          print_err("unsupported export format of an other values vector is required",__FILE__, __LINE__, __func__);
          abort();
        }
      }
      else
      {
        print_err("required index of other values vector exceeds limits\n"
                  " selected range is id1=%ld, n=%ld, max encompo=%ld",
                  __FILE__, __LINE__, __func__, id1, n, Mm->max_ncompothere);
        abort();
      }
    }
  }
}



/**
  Function prints values of selected strains as tensors for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_stra_mtx(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long j;

  if (selestra.st != sel_no)
  {
    for (j=0; j<selestra.n; j++)
    {
      if ((selstra[j].st != sel_mtx) && (selstra[j].st != sel_range_mtx))
        continue;
      if (selstra[j].st == sel_range_mtx)
      {
        print_err("sel_range_mtx selection type of a strain tensor is not supported, use sel_mtx(%d)",
                  __FILE__, __LINE__, __func__, sel_mtx);
        abort();
      }
      if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
        write_gid_elemtensor(out, strain, lcid, j, desclcid);
      else
      {
        print_err("required export format of a strain tensor is not supported", __FILE__, __LINE__, __func__);
        abort();
      }
    }
  }
}



/** 
  Function prints values of selected stresses as tensors for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_stre_mtx(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long j;

  if (selestre.st != sel_no)
  {
    for (j=0; j<selestre.n; j++)
    {
      if ((selstre[j].st != sel_mtx) && (selstre[j].st != sel_range_mtx))
        continue;
      if (selstre[j].st == sel_range_mtx)
      {
        print_err("sel_range_mtx selection type of a stress tensor is not supported, use sel_mtx(%d)", 
                  __FILE__, __LINE__, __func__, sel_mtx);
        abort();
      }
      if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
        write_gid_elemtensor(out, stress, lcid, j, desclcid);
      else
      {
        print_err("required export format of a stress tensor is not supported", __FILE__, __LINE__, __func__);
        abort();
      }
    }
  }
}



/** 
  Function prints values of selected other items as tensors for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_oth_mtx(FILE *out, long lcid, const char *desclcid, graphfmt gf)
{
  long j, id1, n;

  if (seleoth.st != sel_no)
  {
    for (j=0; j<seleoth.n; j++)
    {
      if ((seloth[j].st != sel_mtx) && (seloth[j].st != sel_range_mtx))
        continue;
      if (seloth[j].st == sel_mtx)
      {
        print_err("sel_mtx selection type of an other values tensor is not supported, use sel_range_mtx(%d)",
                  __FILE__, __LINE__, __func__, sel_range_mtx);
        abort();
      }
      id1 = seloth[j].id1[0];
      n  = seloth[j].ncomp[0];
      if (id1+n-1 < Mm->max_ncompothere)
      {
        if ((gf == grfmt_gid) || (gf == grfmt_gid_vtk))
          write_gid_elemtensor(out, other, lcid, j, desclcid);
        else
        {
          print_err("required export format of an other values tensor is not supported", __FILE__, __LINE__, __func__);
          abort();
        }
      }
      else
      {
        print_err("required index of other values vector exceeds limits\n"
                  " selected range is id1=%ld, n=%ld, max encompo=%ld\n",
                  __FILE__, __LINE__, __func__, id1, n, Mm->max_ncompothere);
        abort();
      }
    }
  }
}



/**
  Function prints required output values for selected elements and for given load case and step
  to the several output grahics files named by printed quantity component. If the opening mode is 
  "wt", than the output graphics files are only opened and header is printed, no data are written.
  
  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format
  @param idelem1 - number of the first element for GiD output (normally should be 1), 
                   set by print_graphics function

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_graphics(const char *outfn, const char *mode, long lcid, const char *desclcid, graphfmt gf, long idelem1)
{
  ide1 = idelem1;
  if ((gf != grfmt_gid_sep) && (gf != grfmt_gidsep_vtk))
  {
    print_err("invalid graphics format is required", __FILE__, __LINE__, __func__);
    return;
  }
  print_gr_stra_scal(outfn, mode, lcid, desclcid);
  print_gr_stre_scal(outfn, mode, lcid, desclcid);
  print_gr_oth_scal(outfn, mode, lcid, desclcid);

  print_gr_stra_vec(outfn, mode, lcid, desclcid);
  print_gr_stre_vec(outfn, mode, lcid, desclcid);
  print_gr_oth_vec(outfn, mode, lcid, desclcid);

  print_gr_stra_mtx(outfn, mode, lcid, desclcid);
  print_gr_stre_mtx(outfn, mode, lcid, desclcid);
  print_gr_oth_mtx(outfn, mode, lcid, desclcid);
  return;
}



/** 
  Function prints values of selected strains as scalars for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_stra_scal(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long i, j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selestra.st != sel_no)
  {
    for(i=0; i<Mm->max_ncompstre; i++)
    {
      for (j=0; j<selestra.n; j++)
      {
        if ((selstra[j].st == sel_mtx) || (selstra[j].st == sel_range_mtx) ||
            (selstra[j].st == sel_range_vec) || (selstra[j].st == sel_vec))
          continue;
        if (selstra[j].presence_id(i))
        {
          sprintf(fname, "%s.elem_eps%ld.res", outfn, i+1);
          out = fopen(fname, mode);
          if (out == NULL)
          {
            print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
            abort();
          }
          fseek(out, 0, SEEK_END); // MS Visual C++ requires that
          if (ftell(out) == 0)
          {
            fprintf(out, "GiD Post Results File 1.0\n");
            export_gid_gauss_pt(out, ide1);
          }
          else
          {
            if (desclcid) // regular step output
              write_gid_elemscalar(out, strain, lcid, i, desclcid);
            else{
              // do nothing because restorage from backup is just performed
            }
          }
          fclose(out);
          break;
        }
      }
    }
  }
}



/**  
  Function prints values of selected stresses as scalars for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_stre_scal(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long i, j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selestre.st != sel_no)
  {
    for(i=0; i<Mm->max_ncompstre; i++)
    {
      for (j=0; j<selestre.n; j++)
      {
        if ((selstre[j].st == sel_mtx) || (selstre[j].st == sel_range_mtx) ||
            (selstre[j].st == sel_range_vec) || (selstre[j].st == sel_vec))
          continue;
        if (selstre[j].presence_id(i))
        {
          sprintf(fname, "%s.elem_sig%ld.res", outfn, i+1);
          out = fopen(fname, mode);
          if (out == NULL)
          {
            print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
            abort();
          }
          fseek(out, 0, SEEK_END); // MS Visual C++ requires that
          if (ftell(out) == 0)
          {
            fprintf(out, "GiD Post Results File 1.0\n");
            export_gid_gauss_pt(out, ide1);
          }
          else
          {
            if (desclcid) // regular step output
              write_gid_elemscalar(out, stress, lcid, i, desclcid);
            else{
              // do nothing because restorage from backup is just performed
            }
          }
          fclose(out);
          break;
        }
      }
    }
  }
}



/**  
  Function prints values of selected other items as scalars for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_oth_scal(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long i, j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (seleoth.st != sel_no)
  {
    for(i=0; i<Mm->max_ncompothere; i++)
    {
      for (j=0; j<seleoth.n; j++)
      {
        if ((seloth[j].st == sel_mtx) || (seloth[j].st == sel_range_mtx) ||
            (seloth[j].st == sel_range_vec) || (seloth[j].st == sel_vec))
          continue;
        if (seloth[j].presence_id(i))
        {
          sprintf(fname, "%s.elem_other%ld.res", outfn, i+1);
          out = fopen(fname, mode);
          if (out == NULL)
          {
            print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
            abort();
          }
          fseek(out, 0, SEEK_END); // MS Visual C++ requires that
          if (ftell(out) == 0)
          {
            fprintf(out, "GiD Post Results File 1.0\n");
            export_gid_gauss_pt(out, ide1);
          }
          else
          {
            if (desclcid) // regular step output
              write_gid_elemscalar(out, other, lcid, i, desclcid);
            else{
              // do nothing because restorage from backup is just performed
            }
          }
          fclose(out);
          break;
        }
      }
    }
  }
}



/**
  Function prints values of selected strains as vectors for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_stra_vec(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long j, id1, n;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selestra.st != sel_no)
  {
    for (j=0; j<selestra.n; j++)
    {
      if (selstra[j].st != sel_range_vec)
        continue;
      id1 = selstra[j].id1[0];
      n  = selstra[j].ncomp[0];
      if (id1+n-1 < Mm->max_ncompstre)
      {
        sprintf(fname, "%s.elem_eps_v%ld-%ld_s%ld.res", outfn, id1+1, n, j+1);
        out = fopen(fname, mode);
        if (out == NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        fseek(out, 0, SEEK_END); // MS Visual C++ requires that
        if (ftell(out) == 0)
        {
          fprintf(out, "GiD Post Results File 1.0\n");
          export_gid_gauss_pt(out, ide1);
        }
        else
        {
          if (desclcid) // regular step output
            write_gid_elemvector(out, strain, lcid, j, desclcid);
          else{
            // do nothing because restorage from backup is just performed
          }
        }
        fclose(out);
      }
      else
      {
        print_err("required index of strain vector exceeds limits\n"
                  " selected range is id1=%ld, n=%ld, max ncompst=%ld\n",
                  __FILE__, __LINE__, __func__, id1, n, Mm->max_ncompstre);
        abort();
      }
    }
  }
}



/**
  Function prints values of selected stresses as vectors for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_stre_vec(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long j, id1, n;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selestre.st != sel_no)
  {
    for (j=0; j<selestre.n; j++)
    {
      if (selstre[j].st != sel_range_vec)
        continue;
      id1 = selstre[j].id1[0];
      n  = selstre[j].ncomp[0];
      if (id1+n-1 < Mm->max_ncompstre)
      {
        sprintf(fname, "%s.elem_sig_v%ld-%ld_s%ld.res", outfn, id1+1, n, j+1);
        out = fopen(fname, mode);
        if (out == NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        fseek(out, 0, SEEK_END); // MS Visual C++ requires that
        if (ftell(out) == 0)
        {
          fprintf(out, "GiD Post Results File 1.0\n");
          export_gid_gauss_pt(out, ide1);
        }
        else
        {
          if (desclcid) // regular step output
            write_gid_elemvector(out, stress, lcid, j, desclcid);
          else{
            // do nothing because restorage from backup is just performed
          }
        }
        fclose(out);
      }
      else
      {
        print_err("required index of stress vector exceeds limits\n"
                  " selected range is id1=%ld, n=%ld, max ncompst=%ld\n",
                  __FILE__, __LINE__, __func__, id1, n, Mm->max_ncompstre);
        abort();
      }
    }
  }
}



/**
  Function prints values of selected other items as vectors for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_oth_vec(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long id1, n, j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (seleoth.st != sel_no)
  {
    for (j=0; j<seleoth.n; j++)
    {
      if (seloth[j].st != sel_range_vec)
        continue;
      id1 = seloth[j].id1[0];
      n  = seloth[j].ncomp[0];
      if (id1+n-1 < Mm->max_ncompothere)
      {
        sprintf(fname, "%s.elem_other_v%ld-%ld_s%ld.res", outfn, id1+1, n, j+1);
        out = fopen(fname, mode);
        if (out == NULL)
        {
          print_err("cannot open graphics file '%s'",__FILE__, __LINE__, __func__, fname);
          abort();
        }
        fseek(out, 0, SEEK_END); // MS Visual C++ requires that
        if (ftell(out) == 0)
        {
          fprintf(out, "GiD Post Results File 1.0\n");
          export_gid_gauss_pt(out, ide1);
        }
        else
        {
          if (desclcid) // regular step output
            write_gid_elemvector(out, other, lcid, j, desclcid);
          else{
            // do nothing because restorage from backup is just performed
          }
        }
        fclose(out);
      }
      else
      {
        print_err("required index of other values vector exceeds limits\n"
                  " selected range is id1=%ld, n=%ld, max encompo=%ld\n",
                  __FILE__, __LINE__, __func__, id1, n, Mm->max_ncompothere);
        abort();
      }
    }
  }
}



/**
  Function prints values of selected strains as tensors for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_stra_mtx(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selestra.st != sel_no)
  {
    for (j=0; j<selestra.n; j++)
    {
      if ((selstra[j].st != sel_mtx) && (selstra[j].st != sel_range_mtx))
        continue;
      if (selstra[j].st == sel_range_mtx)
      {
        print_err("sel_range_mtx selection type of a strain tensor is required, use sel_mtx(%d)\n",
                  __FILE__, __LINE__, __func__, sel_mtx);
        abort();
      }
      sprintf(fname, "%s.elem_eps_m_s%ld.res", outfn, j+1);
      out = fopen(fname, mode);
      if (out == NULL)
      {
        print_err("cannot open graphics file '%s'",__FILE__, __LINE__, __func__, fname);
        abort();
      }
      fseek(out, 0, SEEK_END); // MS Visual C++ requires that
      if (ftell(out) == 0)
      {
        fprintf(out, "GiD Post Results File 1.0\n");
        export_gid_gauss_pt(out, ide1);
      }
      else
      {
        if (desclcid) // regular step output
          write_gid_elemtensor(out, strain, lcid, j, desclcid);
        else{
          // do nothing because restorage from backup is just performed
        }
      }
      fclose(out);
    }
  }
}



/**
  Function prints values of selected stresses as tensors for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_stre_mtx(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selestre.st != sel_no)
  {
    for (j=0; j<selestre.n; j++)
    {
      if ((selstre[j].st != sel_mtx) && (selstre[j].st != sel_range_mtx))
        continue;
      if (selstre[j].st == sel_range_mtx)
      {
        print_err("sel_range_mtx selection type of a stress tensor is not supported, use sel_mtx(%d)\n",
                  __FILE__, __LINE__, __func__, sel_mtx);
        abort();
      }
      sprintf(fname, "%s.elem_sig_m_s%ld.res", outfn, j+1);
      out = fopen(fname, mode);
      if (out == NULL)
      {
        print_err("cannot open graphics file '%s'",__FILE__, __LINE__, __func__, fname);
        abort();
      }
      fseek(out, 0, SEEK_END); // MS Visual C++ requires that
      if (ftell(out) == 0)
      {
        fprintf(out, "GiD Post Results File 1.0\n");
        export_gid_gauss_pt(out, ide1);
      }
      else
      {
        if (desclcid) // regular step output
          write_gid_elemtensor(out, stress, lcid, j, desclcid);
        else{
          // do nothing because restorage from backup is just performed
        }
      }
      fclose(out);
    }
  }
}



/**  
  Function prints values of selected other items as tensors for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgm::print_gr_oth_mtx(const char *outfn, const char *mode, long lcid, const char *desclcid)
{
  long j, id1, n;
  char fname[FNAMELEN+70];
  FILE *out;

  if (seleoth.st != sel_no)
  {
    for (j=0; j<seleoth.n; j++)
    {
      if ((seloth[j].st != sel_mtx) && (seloth[j].st != sel_range_mtx))
        continue;
      if (seloth[j].st == sel_mtx)
      {
        print_err("sel_mtx selection type of an other values tensor is not supported, use sel_range_mtx(%d)",
                  __FILE__, __LINE__, __func__, sel_range_mtx);
        abort();
      }
      id1 = seloth[j].id1[0];
      n  = seloth[j].ncomp[0];
      if (id1+n-1 < Mm->max_ncompothere)
      {
        sprintf(fname, "%s.elem_other_m%ld-%ld_s%ld.res", outfn, id1+1, n, j+1);
        out = fopen(fname, mode);
        if (out == NULL)
        {
          print_err("cannot open graphics file '%s'", __FILE__, __LINE__, __func__, fname);
          abort();
        }
        fseek(out, 0, SEEK_END); // MS Visual C++ requires that
        if (ftell(out) == 0)
        {
          fprintf(out, "GiD Post Results File 1.0\n");
          export_gid_gauss_pt(out, ide1);
        }
        else
        {
          if (desclcid) // regular step output
            write_gid_elemtensor(out, other, lcid, j, desclcid);
          else{
            // do nothing because restorage from backup is just performed
          }
        }
        fclose(out);
      }
      else
      {
        print_err("required index of other values tensor exceeds limits\n"
                  " selected range is id1=%ld, n=%ld, max encompo=%ld\n",
                  __FILE__, __LINE__, __func__, id1, n, Mm->max_ncompothere);
        abort();
      }
    }
  }
}



/**
  Converts selections given by property id to list or range type
  in all output collections according element properties defined in the topology top.

  @param top - mesh topology with defined properties of elements

  @return The function does not return anything but it changes internal 
          representation selections of all elements that were givne by property id.

  Created by Tomas Koudelka, 9.1.2015
*/
void elemoutgm::conv_sel_prop(siftop *top)
{
  sel *newsel;
  long *newtrans;

  if (selestra.st == sel_prop)
  {
    if (selestra.conv_selprop(top, gelem, selstra, newsel, transtra, newtrans) == 0){
      delete [] selstra;
      selstra = newsel;
      delete [] transtra;
      transtra = newtrans;
    }
  }

  if (selestre.st == sel_prop)
  {
    if (selestre.conv_selprop(top, gelem, selstre, newsel, transtre, newtrans) == 0){
      delete [] selstre;
      selstre = newsel;
      delete [] transtre;
      transtre = newtrans;
    }
  }

  if (seleoth.st == sel_prop)
  {
    if (seleoth.conv_selprop(top, gelem, seloth, newsel, NULL, newtrans) == 0){
      delete [] seloth;
      seloth = newsel;
    }
  }
}



/**
  Constructor initializes data to zero values

  Created by Tomas Koudelka
*/
pointoutm::pointoutm()
{
  npnt = 0;
  ksi = eta = zeta = NULL;
  selpnt = selstra = selstre = seloth = NULL;
  transtra = transtre = NULL; 
}



/**
  Destructor deallocates used memory

  Created by Tomas Koudelka
*/
pointoutm::~pointoutm()
{
  delete [] ksi;
  delete [] eta;
  delete [] zeta;
  delete [] selpnt;
  delete [] selstra;
  delete [] selstre;
  delete [] seloth;
  delete [] transtra;
  delete [] transtre;
}



/**
  Function reads data with description for output of user defined point values from the text file.

  @param in - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - error reading step or number of points
  @retval 2 - error reading natural coordinates
  @retval 3 - error reading sets of points on elements
  @retval 4 - error reading transformation id

  Created by Tomas Koudelka
*/
long pointoutm::read(XFILE *in)
{
  long i;
  // step and number of points
  xfscanf(in, "%k", "sel_pointstep");
  dstep.read(in);
  if (dstep.st == sel_no)
    return 0;
  if (xfscanf(in, "%k%ld", "numpoints", &npnt) != 2)
  {
    print_err("cannot read number of points", __FILE__, __LINE__, __func__);
    return 1;
  }
  if (npnt == 0)
    return 0;

  ksi  = new double [npnt];
  eta  = new double [npnt];
  zeta = new double [npnt];
  memset(ksi, 0, sizeof(*ksi)*npnt);
  memset(eta, 0, sizeof(*eta)*npnt);
  memset(zeta, 0, sizeof(*zeta)*npnt);
  for (i=0; i<npnt; i++)
  {
    if (xfscanf(in, "%k %le %k %le %k %le", "ksi", ksi+i, "eta", eta+i, "zeta", zeta+i) != 6)
    {
      print_err("cannot read coordinates of UDP", __FILE__, __LINE__, __func__);
      return 2;
    }
  }

  selelem.read(in);
  switch(selelem.st)
  {
    case sel_no:
      return 0;
    case sel_all:
    case sel_range:
    case sel_list:
      selpnt = new sel[selelem.n];
      for (i=0; i<selelem.n; i++)
        selpnt[i].read(in);
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return 3;
  }

  sellc.read(in);  
  if (sellc.st == sel_no)
    return 0;

  selstra = new sel[npnt];
  selstre = new sel[npnt];
  seloth  = new sel[npnt];
  transtra = new long[npnt];
  memset(transtra, 0, sizeof(transtra)*npnt);
  transtre = new long[npnt];
  memset(transtre, 0, sizeof(transtre)*npnt);

  for (i=0; i<npnt; i++)
  {
    // strains
    selstra[i].read(in);
    if (xfscanf(in, "%k%ld", "transfid", transtra+i) != 2)
    {
      print_err("cannot read strain transformation id", __FILE__, __LINE__, __func__);
      return 4;
    }
    // stresses
    selstre[i].read(in);
    if (xfscanf(in, "%k%ld", "transfid", transtre+i) != 2)
    {
      print_err("cannot read stress transformation id", __FILE__, __LINE__, __func__);
      return 4;
    }
    // other values
    seloth[i].read(in);

    // setup of flags
    if (selstra[i].st != sel_no)
      Mp->straincomp = 1;
    if (selstre[i].st != sel_no)
      Mp->stresscomp = 1;
    //if (seloth[i].st != sel_no)
    //Mp->othercomp = 1;
  }
  return 0;
}



/**
  Function prints data with description for output of element values to the text file.

  @param out - pointer to the opened tetx file

  @return The function does not return anything.
 
  Created by Tomas Koudelka
*/
void pointoutm::print(FILE *out)
{
  long i;
  dstep.print(out);
  if ((dstep.st == sel_no) || (npnt == 0))
  {
    fprintf(out, "\n");
    return;
  }
  fprintf(out, " %ld\n", npnt);
  if (npnt == 0)
    return;
  
  for (i=0; i<npnt; i++)
    fprintf(out, "%e %e %e\n", ksi[i], eta[i], zeta[i]);


  selelem.print(out);
  switch(selelem.st)
  {
    case sel_no:
      return;
    case sel_all:
    case sel_range:
    case sel_list:
      for (i=0; i<selelem.n; i++)
        selpnt[i].print(out);
      break;
    default:
      print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
      return;
  }

  sellc.print(out);  
  if (sellc.st == sel_no)
    return;

  for (i=0; i<npnt; i++)
  {
    // strains
    selstra[i].print(out);
    fprintf(out, "%ld\n", transtra[i]);
    // stresses
    selstre[i].print(out);
    fprintf(out, "%ld\n", transtre[i]);
    // other values
    seloth[i].print(out);
  }
}



/**
  Function prints required output values for user defined point and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id  

  @return The function does not return anything.
 
  Created by Tomas Koudelka
*/
void pointoutm::print_out(FILE *out, long lcid)
{
  long i;

  for (i=0; i<npnt; i++)
  {
    if (selstra[i].st != sel_no)
      fprintf(out, "** UDP strains:\n\n");
  }
  print_stra(out, lcid);

  for (i=0; i<npnt; i++)
  {
    if (selstre[i].st != sel_no)
      fprintf(out, "** UDP stresses:\n\n");
  }
  print_stre(out, lcid);

  for (i=0; i<npnt; i++)
  {
    if (seloth[i].st != sel_no)
      fprintf(out, "** UDP other values:\n\n");
    print_other(out, i);
  }
}



/**
  Function prints required strains for user defined point and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

  @return The function does not return anything.
 
  Created by Tomas Koudelka
*/
void pointoutm::print_stra(FILE *out, long lcid)
{
  long i,j,k,l,n, ncomp, id;

  for (i=0; i<Mt->ne; i++)
  {
    if (selelem.presence_id(i))
    {
      fprintf(out, "   Element %7ld:\n",i+1);
      ncomp = Mt->give_tncomp(i);
      id = lcid*ncomp;
      switch(selelem.st)
      {
        case sel_no:
          n=0;
          break;
        case sel_all:
        case sel_range:
        case sel_list:
          n=selelem.n;
          break;
        default:
          print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
          return;
      }
      for (j=0; j<n; j++)
      {
        for (k=0; k<npnt; k++)
        {
          if (selpnt[j].presence_id(k))
          {
            fprintf(out, "    point %ld:", k+1);
            for (l=0; l<ncomp; l++)
            {
              if(selstra[k].presence_id(l))
                fprintf(out, "   eps_%ld=% *e", l+1, prstra, Mm->stra.ev[i][k][id+l]);
            }
            fprintf(out, "\n"); 
          }
        }
      }
      fprintf(out, "\n");
    }
  } 
}



/**
  Function prints required stresses for user defined point and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id
  @param pid  - user defined point id

  @return The function does not return anything.
 
  Created by Tomas Koudelka
*/
void pointoutm::print_stre(FILE *out, long lcid)
{
  long i,j,k,l,n, ncomp, id;

  for (i=0; i<Mt->ne; i++)
  {
    if (selelem.presence_id(i))
    {
      fprintf(out, "   Element %7ld:\n", i+1);
      ncomp = Mt->give_tncomp(i);
      id = lcid*ncomp;
      switch(selelem.st)
      {
        case sel_no:
          n=0;
          break;
        case sel_all:
        case sel_range:
        case sel_list:
          n=selelem.n;
          break;
        default:
          print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
          return;
      }
      for (j=0; j<n; j++)
      {
        for (k=0; k<npnt; k++)
        {
          if (selpnt[j].presence_id(k))
          {
            fprintf(out, "    point %ld:", k+1);
            for (l=0; l<ncomp; l++)
            {
              if(selstra[k].presence_id(l))
                fprintf(out, "   sig_%ld=% *e", l+1, prstre, Mm->stre.ev[i][k][id+l]);
            }
            fprintf(out, "\n"); 
          }
        }
      }
      fprintf(out, "\n");
    }
  } 
}



/**
  Function prints required other values for user defined point and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param pid  - user defined point id

  @return The function does not return anything.
 
  Created by Tomas Koudelka
*/
void pointoutm::print_other(FILE */*out*/, long /*pid*/)
{
  long i,j,k,l,n, ncomp;

  for (i=0; i<Mt->ne; i++)
  {
    if (selelem.presence_id(i))
    {
//      fprintf(out, "   Element %7ld:\n"i+1);
      ncomp = Mm->ip[Mt->elements[i].ipp[0][0]].ncompother;
      switch(selelem.st)
      {
        case sel_no:
          n=0;
          break;
        case sel_all:
        case sel_range:
        case sel_list:
          n=selelem.n;
          break;
        default:
          print_err("unknown type of selection is required", __FILE__, __LINE__, __func__);
          return;
      }
      for (j=0; j<n; j++)
      {
        for (k=0; k<npnt; k++)
        {
          if (selpnt[j].presence_id(k))
          {
//            fprintf(out, "    point %ld:", k+1);
            for (l=0; l<ncomp; l++)
            {
//              if(seloth[k].presence_id(l))
//                fprintf(out, "   other_%ld=% *e", l+1, proth, Mm->other.ev[i][k][l]);
            }
//            fprintf(out, "\n"); 
          }
        }
      }
//      fprintf(out, "\n");
    }
  } 
}



/**
  The function prints required quantities to the separated files.

  @param[in] lcid - required load case id
  @param[in] time - actual time or load parameter
  @param[in] istep - actual step id
  @param[in] r  - array of the nodal values
  @param[in] fl - array with values of the load vector
  @param[in] fi - array with values of the internal force vector (it should be almost identical with fl)
  @param[in] fr - array with values of the residual vector (fl-fi)
  
  @return The function does not return anything, but modifies the content of the output files.

  Created by Tomas Koudelka
*/
void outdriverm::print_sep_files_quant(long lcid, double time, long istep, double *r, double *fl, double *fi, double *fr)
{
  long i;
  gnodvalvm gnv(r, fl, fi, fr);
  
  for (i=0; i<nofq; i++){ // loop over all files
    ofq[i].print_quants(lcid, time, istep, lcs, gnv, idn1, ide1);
  }
}




/**
  The function prints header for the plain text output file format.

  @param[in,out] out - pointer to the opened text file.

  @return The function does not return anything, but changes the content of the file out.

  Created by Tomas Koudelka, 10.2023
*/
void outdriverm::print_header_plain_out(FILE *out)
{
  fprintf(out, "%15s ****  *  ****  ****  *\n", " ");
  fprintf(out, "%15s *     *  *     *     *\n", " ");
  fprintf(out, "%15s  *    *  ***   ***   *\n", " ");
  fprintf(out, "%15s   *   *  *     *     *\n", " ");
  fprintf(out, "%15s****   *  *     ****  ****  MEFEL OUTPUT\n", " ");

  fprintf(out, "\n%s\n", Mp->name);
  fprintf(out, "\n\n\n\n\n");
}
