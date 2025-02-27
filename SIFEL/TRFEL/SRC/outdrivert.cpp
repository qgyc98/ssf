#include <string.h>
#include <stdlib.h>
#include "outdrivert.h"
#include "intools.h"
#include "galias.h"
#include "globalt.h"
#include "probdesct.h"
#include "transbclc.h"
#include "transtop.h"
#include "transprint.h"
#include "nodet.h"
#include "elementt.h"
#include "globmatt.h"
#include "intpointst.h"
#include <math.h>

static int prunkn  = 15;
static int prgrad  = 15;
static int prflux  = 15;
static int proth   = 15;
static int preqoth = 15;
//static int preforcet = 15;

/**
  Constructor initializes data members to zero values.
*/
outdrivert::outdrivert()
{
  memset(outfn, 0, sizeof(*outfn)*FNAMELEN);
  memset(outdiagfn, 0, sizeof(*outdiagfn)*FNAMELEN);
  memset(outgrfn, 0, sizeof(*outgrfn)*FNAMELEN);
  outf = outgr = NULL;
  outdiagf = NULL;
  textout = off;
  nlcs = 0;
  nclcs = NULL;
  bvlcs = NULL;
  gf    = grftt_no;
  ncut = 0;
  idn1 = 1;
  ide1 = 1;
  ndiag = 0;
  odiag = NULL;
  vtk_num = 0;
}



/**
  Destructor deallocates used memory.
*/
outdrivert::~outdrivert()
{
  long i;

  delete [] outdiagf;
  delete [] nclcs;
  for (i=0; i<nlcs; i++)
    delete [] bvlcs[i];
  delete [] odiag;
}



/**
  Function reads description of required output from the opened text file in.

  @param in - pointer to opened text file

  @retval 0 - on success
  @retval 1 - error reading output text filename
  @retval 2 - error reading nodal output description
  @retval 3 - error reading element output description
  @retval 4 - error reading user defined point output description
  @retval 5 - error reading output graphics description
  @retval 6 - error reading output diagram description
*/
long outdrivert::read(XFILE *in)
{
  long i;
  const char *str;

  // text output
  xfscanf(in, "%k%m", "textout", &flagsw_kwdset, &textout);
  if (textout)
    {
      if (xfscanf(in," %a", outfn) != 1)
	{
	  fprintf(stderr, "\n\nError reading filename for output\n");
	  fprintf(stderr, " function outdrivert::read (file %s, line %d)\n", __FILE__, __LINE__);
	  return 1;
	}
      fprintf(stdout, "\n Output file name (transport): %s", outfn);

      if (no.read(in))
	return 2; 
      if (no.dstep.st == sel_no)
	str = " not";
      else
	str = "";
      fprintf(stdout, "\n Output of transport nodal values will%s be performed", str);

      if (eo.read(in))
	return 3;
      if (eo.dstep.st == sel_no)
	str = " not";
      else
	str = "";
      fprintf(stdout, "\n Output of transport element values will%s be performed", str);

      xfscanf(in, "%k", "user_defined_points");
      if (po.read(in))
	return 4;
      if (po.dstep.st == sel_no)
	str = " not";
      else
	str = "";
      fprintf(stdout, "\n Output of transport UDP values will%s be performed", str);
    }

  
  // graphic output
  if (xfscanf(in, "%k%m", "outgr_format", &graphftt_kwdset, (int *)&gf) != 2)
    {
      fprintf(stderr, "\n\nError reading type of grahics output\n");
      fprintf(stderr, " function outdrivert::read (file %s, line %d)\n", __FILE__, __LINE__);
      return 5;
    }
  
  if (gf == grftt_open_dx)
    {
      Tp->gradcomp = 1;
      Tp->gradpos  = 2;
      Tp->gradaver = 1;
      Tp->fluxcomp = 1;
      Tp->fluxpos  = 2;
      Tp->fluxaver = 1;
      Tp->othercomp = 1;
      Tp->otherpos  = 2;
      Tp->otheraver = 1;
      Tp->eqothercomp = 1;
      Tp->eqotherpos  = 2;
      Tp->eqotheraver = 1;
    }
  
  if (gf != grftt_no)
    {
      if (xfscanf(in, " %a", outgrfn) != 1)
	{
	  fprintf(stderr, "\n\nError reading filename for graphics output\n");
	  fprintf(stderr, " function outdrivert::read (file %s, line %d)\n", __FILE__, __LINE__);
	  return 5;
	}
    }

  if (Mesprt)
  {
    switch (gf)
    {
    case grftt_no:
      fprintf(stdout, "\n Graphic output will not be performed");
      break;
    case grftt_open_dx:
      fprintf(stdout, "\n Graphic output will be in OpenDX format");
      fprintf(stdout, "\n Filename of graphic output : %s", outgrfn);
      break;
    case grftt_femcad:
      fprintf(stdout, "\n Graphic output will be in FemCAD format");
      fprintf(stdout, "\n Filename of graphic output : %s", outgrfn);
      break;
    case grftt_gid:
      fprintf(stdout, "\n Graphic output will be in GiD format");
      fprintf(stdout, "\n Filename of graphic output : %s", outgrfn);
      break;
    case grftt_gid_sep:
      fprintf(stdout, "\n Graphic output will be in GiD format separated to quantities");
      fprintf(stdout, "\n Filename of graphic output : %s", outgrfn);
      break;
    case grftt_vtk:
      fprintf(stdout, "\n Graphic output will be in VTK format");
      fprintf(stdout, "\n Filename of graphic output : %s", outgrfn);
      break;
    case grftt_gid_vtk:
      fprintf(stdout, "\n Graphic output will be in GiD and VTK format");
      fprintf(stdout, "\n Filename of graphic output : %s", outgrfn);
      break;
    case grftt_gidsep_vtk:
      fprintf(stdout, "\n Graphic output will be in GiD format separated to quantities and VTK");
      fprintf(stdout, "\n Filename of graphic output : %s", outgrfn);
      break;
    }
  }

  if (gf != grftt_no)
    {
      if (nog.read(in))
	return 6; 
      if (nog.dstep.st == sel_no)
	str = " not";
      else
	str = "";
      if (Mesprt)
        fprintf(stdout, "\n Graphical output of transport nodal values will%s be performed", str);
      
      if (eog.read(in))
	return 7;
      if (eog.dstep.st == sel_no)
	str = " not";
      else
	str = "";
      if (Mesprt)
        fprintf(stdout, "\n Graphical output of transport element values will%s be performed", str);
    }

  //diagrams output
  switch (Tp->tprob)
    {
    case nonstationary_problem:
    case nonlinear_nonstationary_problem:
    case discont_nonstat_problem:
    case discont_nonlin_nonstat_problem:
    case growing_np_problem:
    case growing_np_problem_nonlin:
      if (xfscanf(in, "%k%ld", "numdiag", &ndiag) != 2)
	{
	  fprintf(stderr, "\n\nError reading number of diagrams\n");
	  fprintf(stderr, " function outdrivert::read (file %s, line %d)\n", __FILE__, __LINE__);
	  return 8;
	}
      if (ndiag == 0)
      {
        if (Mesprt)
	  fprintf (stdout, "\n Number of diagrams : %ld", ndiag);
        return 0;
      }
      if (xfscanf(in, " %a", outdiagfn) != 1)
	{
	  fprintf(stderr, "\n\nError reading filename for diagrams\n");
	  fprintf(stderr, " function outdrivert::read (file %s, line %d)\n", __FILE__, __LINE__);
	  return 9;
	}
      if (Mesprt)
      {
        fprintf (stdout, "\n Number of diagrams : %ld", ndiag);
        fprintf (stdout, "\n Diagram filename : %s",outdiagfn);
      }
      odiag = new outdiagt[ndiag];
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
  return 0;
}


/**
  Function prints data with output description to the text file given by out.
  
  @param out - pointer to opened text file for output    
*/
void outdrivert::print(FILE *out)
{
  long i;

  fprintf(out, "\n\n#\n# outdriver section\n#\n");
  fprintf(out, "%d\n", (int)textout);
  
  if (textout == on){
    fprintf(out, "\n%s\n", outfn);
    no.print(out);
    eo.print(out);
    po.print(out);
  }
  
  fprintf(out, "%d\n", (int)gf);
  if (gf != grftt_no){
    fprintf(out, "%s\n", outgrfn);
    nog.print(out);
    eog.print(out);
  }

  switch (Tp->tprob)
  {
  case nonstationary_problem:
  case nonlinear_nonstationary_problem:
  case discont_nonstat_problem:
  case discont_nonlin_nonstat_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:
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
}


/**
  Function prints header to the output text file.
  
  @param out - pointer to the opened text output file
*/
void outdrivert::print_header(FILE *out)
{
  fprintf(out, "%15s ****  *  ****  ****  *\n", " ");
  fprintf(out, "%15s *     *  *     *     *\n", " ");
  fprintf(out, "%15s  *    *  ***   ***   *\n", " ");
  fprintf(out, "%15s   *   *  *     *     *\n", " ");
  fprintf(out, "%15s****   *  *     ****  ****  TRFEL OUTPUT\n", " ");

  fprintf(out, "\n%s\n", Tp->name);
  fprintf(out, "\n\n\n\n\n");
}



/**
  Function prints step number to the output text file.
  
  @param out  - pointer to the opened text output file
  @param lcid - load case id
  @param step - integer step id
  @param time - time or load step
*/
void outdrivert::print_newstep(FILE *out, long /*lcid*/, long istep, double time)
{
  long i;
  double lambda;

  if (textout == on)
  {

    switch (Tp->tprt)
      {
      case secondst:
	lambda = time;
	break;
      case minutest:
	lambda = time/60.0;
	break;
      case hourst:
	lambda = time/3600.0;
	break;
      case dayst:
	lambda = time/86400.0;
	break;
      default:
	{
	  print_err("unknown type of time printing",__FILE__,__LINE__,__func__);
	}	
      }
    
    if (no.dstep.presence_id(istep, time, Tp->timecont) ||
        eo.dstep.presence_id(istep, time, Tp->timecont) ||
        po.dstep.presence_id(istep, time, Tp->timecont))
    {
      for (i=0; i<53; i++)
        fprintf(out, "*");
      fprintf(out, "\n%10sStep number=% ld, time/load step=% .15g\n", " ", istep, lambda);
      for (i=0; i<53; i++)
        fprintf(out, "*");
      fprintf(out, "\n\n\n\n");
    }
  }
}


/**
  Function prints required output values to the output text file.
  
  @param out  - pointer to the opened text file 
  @param lcid  - load case id
  @param istep - step id
  @param time  - actual time step
  
*/
void outdrivert::print_out(FILE *out, long lcid, long istep, double time)
{
  if (textout == on)
  {
    if (no.dstep.presence_id(istep, time, Tp->timecont))
      no.print_out(out, lcid);
    if (eo.dstep.presence_id(istep, time, Tp->timecont))
      eo.print_out(out, lcid);
    if (po.dstep.presence_id(istep, time, Tp->timecont)) 
      po.print_out(out, lcid);
  }
}



/**
  Function prints required output values to the output text file.
  
  @param out  - pointer to the opened text file 
  @param lcid  - load case id
  @param istep - step id
  @param time  - actual time step

  Created by Tomas Koudelka, 2014
  
*/
void outdrivert::print_out_forced(FILE *out, long lcid, long /*istep*/, double /*time*/)
{
  if (textout == on)
  {
    no.print_out(out, lcid);
    eo.print_out(out, lcid);
    po.print_out(out, lcid);
  }
}



/**
  Function prints diagrams.
  
  @param lcid - required load case id
  @param lambda - actual load parameter/actual time
  @param istep - actual step id
  @param fi - array with values of load vector
  
*/
void outdrivert::print_diags(long lcid, double lambda, long istep, double *fi)
{
  long i;
  for(i=0; i<ndiag; i++)
    odiag[i].printval(outdiagf[i], lcid, lambda, istep, fi);
}



/**
  Function prints diagrams.
  
  @param lcid - required load case id
  @param lambda - actual load parameter/actual time
  @param istep - actual step id
  @param fi - array with values of load vector
  
  Created by Tomas Koudelka, 2014
*/
void outdrivert::print_diags_forced(long lcid, double lambda, long istep, double *fi)
{
  long i;
  for(i=0; i<ndiag; i++)
    odiag[i].printval_forced(outdiagf[i], lcid, lambda, istep, fi);
}



/**
  Function prints required value to the graphics files.
  
  @param out - pointer to the opened text file 
  @param lcid - required load case id
  @param lambda - actual load parameter/actual time
  @param istep - actual step id
  @param fi - array with values of load vector
  
*/
void outdrivert::print_graphics(FILE *out, long lcid, double lambda, long istep, double *fi)
{
  long i;
  char dlcid[50];


  switch (gf)
  {
    case grftt_no:
      break;;
    case grftt_femcad:
      break;
    case grftt_open_dx:
    {
      for (i=0; i<1; i++)//for (i=0; i<Tb->nlc; i++)
      {
        //if (gf == grftt_open_dx)
        //print_default_dx (Gtt,i,outgrfn);
      }
      break;
    }
    case grftt_gid:
    {
      if (Tp->tprob == stationary_problem)
        sprintf(dlcid, "%ld", istep);
      else
      {
        switch (Tp->tprt){
          case secondst:
            sprintf(dlcid, "%.15e", lambda);
            break;
          case minutest:
            sprintf(dlcid, "%.15e", lambda/60.0);
            break;
          case hourst:
            sprintf(dlcid, "%.15e", lambda/3600.0);
            break;
          case dayst:
            sprintf(dlcid, "%.15e", lambda/86400.0);
            break;
          default:
            print_err("unknown type of time printing",__FILE__,__LINE__,__func__);
        }
      }
      if ((nog.dstep.presence_id(istep, lambda, Tp->timecont)))
        nog.print_graphics(out, lcid, dlcid, gf, fi);
	
      if ((eog.dstep.presence_id(istep, lambda, Tp->timecont)))
        eog.print_graphics(out, lcid, dlcid, gf, ide1);
      break;
    }
    case grftt_gid_sep:{
      if (Tp->tprob == stationary_problem)
        sprintf(dlcid, "%ld", istep);
      else{
        switch (Tp->tprt){
          case secondst:
            sprintf(dlcid, "%.15e", lambda);
            break;
          case minutest:
            sprintf(dlcid, "%.15e", lambda/60.0);
            break;
          case hourst:
            sprintf(dlcid, "%.15e", lambda/3600.0);
            break;
          case dayst:
            sprintf(dlcid, "%.15e", lambda/86400.0);
            break;
          default:
            print_err("unknown type of time printing",__FILE__,__LINE__,__func__);
        }
      }
      if ((nog.dstep.presence_id(istep, lambda, Tp->timecont)))
        nog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, fi, hf_off);

      if ((eog.dstep.presence_id(istep, lambda, Tp->timecont)))
        eog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, hf_off, ide1);
      break;
    }
    case grftt_vtk:{
      //works only for nonstationary problem
      if (nog.dstep.presence_id(istep, lambda, Tp->timecont))
      {
        char fname[FNAMELEN+70];
        FILE *vtk_file;
        char *path;
        char *name;
        char *suffix;

        filename_decomposition(outgrfn,path,name,suffix);
        sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
        delete [] path; delete [] name; delete [] suffix;

        if((vtk_file = fopen(fname, "w"))==NULL){
          fprintf(stderr, "\n\nUnable to open graphics file '%s' in function print_graphics()\n", fname);
          fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
          abort();
        }
        print_default_vtk(vtk_file,lcid,istep,lambda);//print header,points,nodes,elements
        write_vtk_unkn(vtk_file,lcid,vtk_num==0?1:0,fi);//print datapoints in sections
        fclose(vtk_file);
        vtk_num++;
      }
      break;
    }
    case grftt_gid_vtk:{
      if (Tp->tprob == stationary_problem)
        sprintf(dlcid, "%ld", istep);
      else{
        switch (Tp->tprt){
          case secondst:
            sprintf(dlcid, "%.15e", lambda);
            break;
          case minutest:
            sprintf(dlcid, "%.15e", lambda/60.0);
            break;
          case hourst:
            sprintf(dlcid, "%.15e", lambda/3600.0);
            break;
          case dayst:
            sprintf(dlcid, "%.15e", lambda/86400.0);
            break;
          default:
            print_err("unknown type of time printing",__FILE__,__LINE__,__func__);
        }
      }
      if ((nog.dstep.presence_id(istep, lambda, Tp->timecont)))
      {
        // output of nodal values in GiD format
        nog.print_graphics(out, lcid, dlcid, gf, fi);


        // output of nodal values in VTK format
        char fname[FNAMELEN+70];
        FILE *vtk_file;
        char *path;
        char *name;
        char *suffix;

        filename_decomposition(outgrfn,path,name,suffix);
        sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
        delete [] path; delete [] name; delete [] suffix;

        if((vtk_file = fopen(fname, "w"))==NULL){
          fprintf(stderr, "\n\nUnable to open graphics file '%s' in function print_graphics()\n", fname);
          fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
          abort();
        }
        print_default_vtk(vtk_file,lcid,istep,lambda);//print header,points,nodes,elements
        write_vtk_unkn(vtk_file,lcid,vtk_num==0?1:0,fi);//print datapoints in sections
        fclose(vtk_file);
        vtk_num++;
      }	
      if ((eog.dstep.presence_id(istep, lambda, Tp->timecont)))
        eog.print_graphics(out, lcid, dlcid, gf, ide1);

      break;
    }
    case grftt_gidsep_vtk:{
      if (Tp->tprob == stationary_problem)
        sprintf(dlcid, "%ld", istep);
      else{
        switch (Tp->tprt){
          case secondst:
            sprintf(dlcid, "%.15e", lambda);
            break;
          case minutest:
            sprintf(dlcid, "%.15e", lambda/60.0);
            break;
          case hourst:
            sprintf(dlcid, "%.15e", lambda/3600.0);
            break;
          case dayst:
            sprintf(dlcid, "%.15e", lambda/86400.0);
            break;
          default:
            print_err("unknown type of time printing",__FILE__,__LINE__,__func__);
        }
      }
      if ((nog.dstep.presence_id(istep, lambda, Tp->timecont)))
      {
        nog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, fi, hf_off);
        // output of nodal values in VTK format
        char fname[FNAMELEN+70];
        FILE *vtk_file;
        char *path;
        char *name;
        char *suffix;

        filename_decomposition(outgrfn,path,name,suffix);
        sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
        delete [] path; delete [] name; delete [] suffix;

        if((vtk_file = fopen(fname, "w"))==NULL){
          fprintf(stderr, "\n\nUnable to open graphics file '%s' in function print_graphics()\n", fname);
          fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
          abort();
        }
        print_default_vtk(vtk_file,lcid,istep,lambda);//print header,points,nodes,elements
        write_vtk_unkn(vtk_file,lcid,vtk_num==0?1:0,fi);//print datapoints in sections
        fclose(vtk_file);
        vtk_num++;
      }
      if ((eog.dstep.presence_id(istep, lambda, Tp->timecont)))
        eog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, hf_off, ide1);
      break;
    }
    default:
      print_err("unknown type of graphics format is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function prints required value to the graphics files.
  
  @param out - pointer to the opened text file 
  @param lcid - required load case id
  @param lambda - actual load parameter/actual time
  @param istep - actual step id
  @param fi - array with values of load vector
  
  Created by Tomas Koudelka, 2014
*/
void outdrivert::print_graphics_forced(FILE *out, long lcid, double lambda, long istep, double *fi)
{
  long i;
  char dlcid[50];


  switch (gf)
    {
    case grftt_no:
      break;;
    case grftt_femcad:
      break;
    case grftt_open_dx:
      {
	for (i=0; i<1; i++)//for (i=0; i<Tb->nlc; i++)
	  {
	    //if (gf == grftt_open_dx)
	    //print_default_dx (Gtt,i,outgrfn);
	  }
	break;
      }
    case grftt_gid:
    {
      if (Tp->tprob == stationary_problem)
        sprintf(dlcid, "%ld", istep);
      else
      {
        switch (Tp->tprt)
        {
          case secondst:
            sprintf(dlcid, "%.15e", lambda);
            break;
          case minutest:
            sprintf(dlcid, "%.15e", lambda/60.0);
            break;
          case hourst:
            sprintf(dlcid, "%.15e", lambda/3600.0);
            break;
          case dayst:
            sprintf(dlcid, "%.15e", lambda/86400.0);
            break;
          default:
            print_err("unknown type of time printing",__FILE__,__LINE__,__func__);
        }
      }
      nog.print_graphics(out, lcid, dlcid, gf, fi);      
      eog.print_graphics(out, lcid, dlcid, gf, ide1);
      break;
    }
    case grftt_gid_sep:
    {
      if (Tp->tprob == stationary_problem)
        sprintf(dlcid, "%ld", istep);
      else
      {
        switch (Tp->tprt)
        {
          case secondst:
            sprintf(dlcid, "%.15e", lambda);
            break;
          case minutest:
            sprintf(dlcid, "%.15e", lambda/60.0);
            break;
          case hourst:
            sprintf(dlcid, "%.15e", lambda/3600.0);
            break;
          case dayst:
            sprintf(dlcid, "%.15e", lambda/86400.0);
            break;
          default:
            print_err("unknown type of time printing",__FILE__,__LINE__,__func__);
        }
      }
      nog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, fi, hf_off);        
      eog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, hf_off, ide1);      
      break;
    }
    case grftt_vtk:
    {//works only for nonstationary problem
      char fname[FNAMELEN+70];
      FILE *vtk_file;
      char *path;
      char *name;
      char *suffix;
      
      filename_decomposition(outgrfn,path,name,suffix);
      sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
      delete [] path; delete [] name; delete [] suffix;

      if((vtk_file = fopen(fname, "w"))==NULL) {
        fprintf(stderr, "\n\nUnable to open graphics file '%s' in function print_graphics()\n", fname);
        fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
        abort();
      }
      print_default_vtk(vtk_file,lcid,istep,lambda);//print header,points,nodes,elements
      write_vtk_unkn(vtk_file,lcid,vtk_num==0?1:0,fi);//print datapoints in sections
      fclose(vtk_file);
      vtk_num++;

      break;
    }
    case grftt_gid_vtk:
    {
      //
      // forced output of nodal values in the GiD format
      //
      if (Tp->tprob == stationary_problem)
        sprintf(dlcid, "%ld", istep);
      else
      {
        switch (Tp->tprt)
        {
          case secondst:
            sprintf(dlcid, "%.15e", lambda);
            break;
          case minutest:
            sprintf(dlcid, "%.15e", lambda/60.0);
            break;
          case hourst:
            sprintf(dlcid, "%.15e", lambda/3600.0);
            break;
          case dayst:
            sprintf(dlcid, "%.15e", lambda/86400.0);
            break;
          default:
            print_err("unknown type of time printing",__FILE__,__LINE__,__func__);
        }
      }
      nog.print_graphics(out, lcid, dlcid, gf, fi);      
      eog.print_graphics(out, lcid, dlcid, gf, ide1);

      //
      // forced output of nodal values in the VTK format
      //
      //works only for nonstationary problem
      char fname[FNAMELEN+70];
      FILE *vtk_file;
      char *path;
      char *name;
      char *suffix;
      
      filename_decomposition(outgrfn,path,name,suffix);
      sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
      delete [] path; delete [] name; delete [] suffix;

      if((vtk_file = fopen(fname, "w"))==NULL) {
        fprintf(stderr, "\n\nUnable to open graphics file '%s' in function print_graphics()\n", fname);
        fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
        abort();
      }
      print_default_vtk(vtk_file,lcid,istep,lambda);//print header,points,nodes,elements
      write_vtk_unkn(vtk_file,lcid,vtk_num==0?1:0,fi);//print datapoints in sections
      fclose(vtk_file);
      vtk_num++;
      break;
    }
    case grftt_gidsep_vtk:
    {
      //
      // forced output of nodal values in the GiD format
      //
      if (Tp->tprob == stationary_problem)
        sprintf(dlcid, "%ld", istep);
      else
      {
        switch (Tp->tprt)
        {
          case secondst:
            sprintf(dlcid, "%.15e", lambda);
            break;
          case minutest:
            sprintf(dlcid, "%.15e", lambda/60.0);
            break;
          case hourst:
            sprintf(dlcid, "%.15e", lambda/3600.0);
            break;
          case dayst:
            sprintf(dlcid, "%.15e", lambda/86400.0);
            break;
          default:
            print_err("unknown type of time printing",__FILE__,__LINE__,__func__);
        }
      }
      nog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, fi, hf_off);        
      eog.print_graphics(outgrfngs, "at", lcid, dlcid, gf, hf_off, ide1);      

      //
      // forced output of nodal values in the VTK format
      //
      //works only for nonstationary problem
      char fname[FNAMELEN+70];
      FILE *vtk_file;
      char *path;
      char *name;
      char *suffix;
      filename_decomposition(outgrfn,path,name,suffix);
      sprintf(fname,"%s%s%04ld%s", path, name, istep, suffix);//increase index in filename
      delete [] path; delete [] name; delete [] suffix;

      if((vtk_file = fopen(fname, "w"))==NULL) {
        fprintf(stderr, "\n\nUnable to open graphics file '%s' in function print_graphics()\n", fname);
        fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
        abort();
      }
      print_default_vtk(vtk_file,lcid,istep,lambda);//print header,points,nodes,elements
      write_vtk_unkn(vtk_file,lcid,vtk_num==0?1:0,fi);//print datapoints in sections
      fclose(vtk_file);
      vtk_num++;
      break;
    }
    default:
      print_err("unknown type of graphics format is required", __FILE__, __LINE__, __func__);
  }
}



/**
  Function creates output graphics files and their headers for GiD separated format (grftt_gidsp).
  For this format (grftt_gidsp) each required quantity is printed to the sparated file named 
  in following way : filename.{elem|nodal}_{quanitity}{indexOfQuantity}.flavia.res
  quatnitity can be : uknowns, grad, flux, pgrad, pflux, other, eqother.
  @param mode - opening mode for opened text files

*/
void outdrivert::create_files_gidsp(const char *mode)
{
  const char *dlcid = "";

  if ((gf != grftt_gid_sep) && (gf != grftt_gidsep_vtk))
  {
    fprintf(stderr, "\nWrong graphics format is required in function outdrivert::create_files_gidsp\n");
    fprintf(stderr, "(file %s, line %d)\n", __FILE__, __LINE__);
    abort();
  }
  if (nog.dstep.st != sel_no)
    nog.print_graphics(outgrfngs, mode, 0, dlcid, gf, NULL, header);
  if (eog.dstep.st != sel_no)
    eog.print_graphics(outgrfngs, mode, 0, dlcid, gf, header, ide1);
}




/**
  Function prints footers for GiD separated format (grftt_gidsp) in case of adaptivity.
  For this format (grftt_gidsp) each required quantity is printed to the sparated file named 
  in following way : filename.{elem|nodal}_{quanitity}{indexOfQuantity}.flavia.res
  quatnitity can be : uknowns, grad, flux, pgrad, pflux, other, eqother.

*/
void outdrivert::close_files_gidsp()
{
  const char *dlcid = "";

  if ((gf != grftt_gid_sep) && (gf != grftt_gidsep_vtk))
  {
    fprintf(stderr, "\nWrong graphics format is required in function outdrivert::create_files_gidsp\n");
    fprintf(stderr, "(file %s, line %d)\n", __FILE__, __LINE__);
    abort();
  }
  if (nog.dstep.st != sel_no)
    nog.print_graphics(outgrfngs, "at", 0, dlcid, gf, NULL, footer);
  if (eog.dstep.st != sel_no)
    eog.print_graphics(outgrfngs, "at", 0, dlcid, gf, footer, ide1);
}



/** 
  Converts selections given by property id to list or range type
  in all output collections according nodal and element properties 
  defined in the topology top.

  @param top - mesh topology with defined properties of nodes and elements

  @return The function does not return anything but it changes internal 
          representation selections of all nodes and elements that were givne by property id.

  Created by Tomas Koudelka, 08.2016
*/
void outdrivert::conv_sel_prop(siftop *top)
{
  no.conv_sel_prop(top);
  eo.conv_sel_prop(top);
  nog.conv_sel_prop(top);
  eog.conv_sel_prop(top);
}



/**
  Constructor initializes data to zero values
*/
nodeoutt::nodeoutt()
{
  selunkn = selgrad = selflux = seloth = seleqoth = NULL;
}



/**
  Destructor deallocates used memory
*/
nodeoutt::~nodeoutt()
{
  delete [] selunkn;
  delete [] selgrad;
  delete [] selflux; 
  delete [] seloth;
  delete [] seleqoth;
}



/**
  Function reads data with description for output of nodal values from the text file.

  @param in - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - error reading step or load case selection
  @retval 2 - error reading unknown selection
  @retval 3 - error reading gradient selection
  @retval 4 - error reading flux selection
  @retval 5 - error reading other values selection
  @retval 6 - error reading eqother values selection
 */
long nodeoutt::read(XFILE *in)
{
  long i;
  // step
  xfscanf(in, "%k", "sel_nodstep");
  dstep.read(in);
  if (dstep.st == sel_no)
    return 0;

  // unknowns
  xfscanf(in, "%k", "unknowns_nodes");
  selnunkn.read(in);
  switch (selnunkn.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "unknowns_comp");
      selunkn = new sel[selnunkn.n];
      for (i=0; i<selnunkn.n; i++)
        selunkn[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 2;
  }

  // gradients
 xfscanf(in, "%k", "grad_nodes");
 selngrad.read(in);
  switch (selngrad.st)
    {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "nodgrad_unkn");
      if(Tp->gradcomp == 0)
	Tp->gradcomp = 1;
      if(Tp->gradaver == 0)
      {
	Tp->gradaver = 1;
        Tp->gradpos  = 2;
      }
      selgrad = new sel[selngrad.n];
      for(i=0; i<selngrad.n; i++)
        selgrad[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 3;
    }
  
  
  // fluxes
  xfscanf(in, "%k", "flux_nodes");
 selnflux.read(in);
  switch (selnflux.st)
    {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "nodflux_unkn");
      if(Tp->gradcomp == 0)
	Tp->gradcomp = 1;
      if(Tp->fluxcomp == 0)
	Tp->fluxcomp = 1;
      if(Tp->gradaver == 0)
      {
        Tp->gradpos  = 2;
	Tp->gradaver = 1;
      }
      if(Tp->fluxaver == 0)
      {
        Tp->fluxpos  = 2;
	Tp->fluxaver = 1;
      }
      selflux = new sel[selnflux.n];
      for(i=0; i<selnflux.n; i++)
        selflux[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
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
      xfscanf(in, "%k", "nodother_comp");
      if(Tp->othercomp == 0)
	Tp->othercomp = 1;
      if(Tp->otheraver == 0)
      {
        Tp->otherpos  = 2;
	Tp->otheraver = 1;
      }
      seloth = new sel[selnoth.n];
      for (i=0; i < selnoth.n; i++)
        seloth[i].read(in);
      break; 
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 5;
  }

  // eq_other values
  xfscanf(in, "%k", "eqother_nodes");
  selneqoth.read(in);
  switch (selneqoth.st)
    {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "nodeqother_comp");
      if(Tp->eqothercomp == 0)
	Tp->eqothercomp = 1;
      if(Tp->eqotheraver == 0)
      {
        Tp->eqotherpos  = 2;
	Tp->eqotheraver = 1;
      }
      seleqoth = new sel[selneqoth.n];
      for (i=0; i < selneqoth.n; i++)
        seleqoth[i].read(in);
      break; 
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 6;
  }

  // 
  xfscanf(in, "%k", "flux_nodes");
  selnforcet.read(in);
  switch (selnforcet.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "flux_comp");
      selforcet = new sel[selnforcet.n];
      for (i=0; i<selnunkn.n; i++)
        selforcet[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 2;
  }

  return 0;
}



/**
  Function prints data with description for output of nodal values to the text file.

  @param out - pointer to the opened text file
*/
void nodeoutt::print(FILE *out)
{
  long i;

  // step and loadcases
  dstep.print(out);
  if (dstep.st == sel_no)
    return;
  fprintf(out, "\n");

  // unknowns
  selnunkn.print(out);
  switch(selnunkn.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selnunkn.n; i++)
        selunkn[i].print(out);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
  }


  // gradients
  selngrad.print(out);
  switch (selngrad.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selngrad.n; i++)
        selgrad[i].print(out);
      fprintf(out, "\n");
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
  }

  // fluxes
  selnflux.print(out);
  switch (selnflux.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selnflux.n; i++)
        selflux[i].print(out);
      fprintf(out, "\n");
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
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
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
  }

  // eq_other values
  selneqoth.print(out);
  switch (selneqoth.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selneqoth.n; i++)
        seleqoth[i].print(out);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
  }

  //forces = fluxes
  selnforcet.print(out);
  switch (selnforcet.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for (i=0; i<selnunkn.n; i++)
        selforcet[i].print(out);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
  }

}


/**
  Function prints required output values for selected nodes and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

 */
void nodeoutt::print_out(FILE *out, long lcid)
{
  if (selnunkn.st != sel_no)
  {
    fprintf(out, "** Nodal unknowns :\n\n");
    print_unkn(out, lcid);
  }
  if (selngrad.st != sel_no)
  {
    fprintf(out, "** Nodal averaged gradients of unknowns :\n\n");
    print_grad(out, lcid);
  }
  if (selnflux.st != sel_no)
  {
    fprintf(out, "** Nodal averaged fluxes of unknowns :\n\n");
    print_flux(out, lcid);
  }
  if (selnoth.st != sel_no)
  {
    fprintf(out, "** Nodal averaged other values :\n\n");
    print_other(out);
  }
  if (selneqoth.st != sel_no)
  {
    fprintf(out, "** Nodal averaged eq_other values :\n\n");
    print_eqother(out);
  }

  /* if (selnforcet.st != sel_no) //not used
     {
     fprintf(out, "** Nodal optional variables (fluxes, etc.) :\n\n");
     print_forcet(out);
     }
  */
}




/**
  Function prints required unknowns for selected nodes and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

 */
void nodeoutt::print_unkn(FILE *out, long /*lcid*/)
{
  long i, j, ndofn;
  vector r;

  for (i=0; i<Tt->nn; i++)
  {
    if (selnunkn.presence_id(i))
    {
      fprintf(out, " Node %7ld", i+1);
      ndofn = Tt->give_ndofn(i);
      reallocv(RSTCKVEC(ndofn, r));
      nodalval(i, r);
      for (j=0; j<ndofn; j++)
      {
        if (selnunkn.presence_id(selunkn,i,j)){
          // fprintf(out, "%20.15le",r[j]);
          fprintf(out, "   r_%ld=% .*e", j+1, prunkn, r[j]);
	}
      }
      fprintf(out, "\n");
    } 
  }
  fprintf(out, "\n");
}



/**
  Function prints required gradients for selected nodes and for given load case
  to the output text file - not completed
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

 */
void nodeoutt::print_grad(FILE *out, long /*lcid*/)
{
  long i, j, k, ncomp, ir;//co to je ir??!!
  
  ncomp = Tt->nodes[0].ncompgrad;
  
  for (i=0; i<Tt->nn; i++)
    {
      if (selngrad.presence_id(i, ir))
	{
	  fprintf(out, " Node %7ld\n", i+1);
	  ncomp = Tt->nodes[i].ncompgrad;
	  for (j=0; j<Tp->ntm; j++){
	    if (selngrad.presence_id(selgrad,i,j))
	      for (k=0; k<ncomp; k++)
		{
		  fprintf(out, "   grad_%ld_%ld=% .*e", j+1, k+1, prgrad, Tt->nodes[i].gradient[j][k]);
		}
	    fprintf(out, "\n"); 
	  }
	} 
    }
  fprintf(out, "\n");
}



/**
  Function prints required flux for selected nodes and for given load case
  to the output text file - not completed
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

 */
void nodeoutt::print_flux(FILE *out, long /*lcid*/)
{
  long i, ir, j, k, ncomp;//co to je ir??!!
  
  ncomp = Tt->nodes[0].ncompgrad;
  for (i=0; i<Tt->nn; i++)
    {
      if (selnflux.presence_id(i, ir))
	{
	  fprintf(out, " Node %7ld\n", i+1);
	  ncomp = Tt->nodes[i].ncompgrad;
	  for (j=0; j<Tp->ntm; j++){
	    if (selnflux.presence_id(selflux, i, j))
	      for (k=0; k<ncomp; k++)
		{
		  fprintf(out, "   flux_%ld_%ld=% .*e", j+1, k+1, prflux, Tt->nodes[i].flux[j][k]);
		}
	    fprintf(out, "\n");
	  }
	}
    }
  fprintf(out, "\n");
}



/**
  Function prints required other values for selected nodes and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

 */
void nodeoutt::print_other(FILE *out)
{
  long i, j, ncomp;
  
  ncomp = Tt->nodes[0].ncompother;
  for (i=0; i<Tt->nn; i++)
    {
      if (selnoth.presence_id(i))
	{
	  fprintf(out, " Node %7ld", i+1);
	  ncomp = Tt->nodes[i].ncompother;
	  for (j=0; j<ncomp; j++)
	    {
	      if (selnoth.presence_id(seloth, i, j)){
		fprintf(out, "\n");
		fprintf(out, "   other_%ld ", j+1);
		Tm->give_othervalue_name(out,0,j);
		fprintf(out, " =% .*e", proth, Tt->nodes[i].other[j]);
	      }
	    }
	  fprintf(out, "\n");
	} 
    }
  fprintf(out, "\n");
}




/**
  Function prints required eq_other values for selected nodes and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

 */
void nodeoutt::print_eqother(FILE *out)
{
  long i, j, ncomp;
  
  ncomp = Tt->nodes[0].ncompeqother;
  for (i=0; i<Tt->nn; i++)
    {
      if (selneqoth.presence_id(i))
	{
	  fprintf(out, " Node %7ld", i+1);
	  ncomp = Tt->nodes[i].ncompeqother;
	  for (j=0; j<ncomp; j++)
	    {
	      if (selneqoth.presence_id(seleqoth, i, j)){
		fprintf(out, "\n");
		fprintf(out, "   eqother_%ld ", j+1);
		Tm->give_eqothervalue_name(out,0,j);
		fprintf(out, " =% .*e", preqoth, Tt->nodes[i].eqother[j]);
	      }
	    }
	  fprintf(out, "\n");
	} 
    }
  fprintf(out, "\n");
}



/**
  Function prints all fluxes for for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

 */
void nodeoutt::print_forcet(FILE *out)
{
  /*
  long i, ii, j, ndof;

  for (i=0; i<Tt->nn; i++)
  {
     fprintf(out, " Node %7ld", i+1);
     ndof = Tt->give_ndofn(i);
     for (j=0; j<ndof; j++){
       ii=Tt->give_dof(i,j);
       if (ii<0)   f[j]=0.0;
       if (ii==0)  f[j]=0.0;
       if (ii>0)   f[j]=ifor[ii-1];
       
       fprintf(out, "   flux_%ld=% .*e", j+1, prforcet,f[j]);//not completed and not used
       
       
     }
     fprintf(out, "\n");
  }
  */
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
void nodeoutt::conv_sel_prop(siftop *top)
{
  sel *newsel;
  long *newtrans;

  if (selnunkn.st == sel_prop)
  {
    selnunkn.conv_selprop(top, gnod, selunkn, newsel, NULL, newtrans);
    delete [] selunkn;
    selunkn = newsel;
  }

  if (selngrad.st == sel_prop)
  {
    selngrad.conv_selprop(top, gnod, selgrad, newsel, NULL, newtrans);
    delete [] selgrad;
    selgrad = newsel;
  }

  if (selnflux.st == sel_prop)
  {
    selnflux.conv_selprop(top, gnod, selflux, newsel, NULL, newtrans);
    delete [] selflux;
    selflux = newsel;
  }

  if (selnoth.st == sel_prop)
  {
    selnoth.conv_selprop(top, gnod, seloth, newsel, NULL, newtrans);
    delete [] seloth;
    seloth = newsel;
  }
}








/**
  Constructor initializes data to zero values
*/
nodeoutgt::nodeoutgt()
{
  selunkn = selgrad = selflux = seloth = seleqoth = NULL;
}



/**
  Destructor deallocates used memory
*/
nodeoutgt::~nodeoutgt()
{
  delete [] selunkn;
  delete [] selgrad;
  delete [] selflux; 
  delete [] seloth;
  delete [] seleqoth;
}



/**
  Function reads data with description for output of nodal values from the text file.

  @param in - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - error reading step or load case selection
  @retval 2 - error reading unknown selection
  @retval 3 - error reading gradient selection
  @retval 4 - error reading flux selection
  @retval 5 - error reading other values selection
  @retval 6 - error reading eqother values selection
 */
long nodeoutgt::read(XFILE *in)
{
  long i;
  // step
  xfscanf(in, "%k", "sel_nodstep");
  dstep.read(in);
  if (dstep.st == sel_no)
    return 0;

  // unknowns
  xfscanf(in, "%k", "unknowns_nodes");
  selnunkn.read(in);
  switch (selnunkn.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "unknowns_comp");
      selunkn = new sel[selnunkn.n];
      for (i=0; i<selnunkn.n; i++)
        selunkn[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 2;
  }

  // gradients
  xfscanf(in, "%k", "grad_nodes");
  selngrad.read(in);
  switch (selngrad.st)
    {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "nodgrad_unkn");
      if(Tp->gradcomp == 0)
	Tp->gradcomp = 1;
      if(Tp->gradaver == 0)
      {
        Tp->gradpos  = 2;
	Tp->gradaver = 1;
      }
      selgrad = new sel[selngrad.n];
      for(i=0; i<selngrad.n; i++)
        selgrad[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 3;
    }
  
  
  // fluxes
  xfscanf(in, "%k", "flux_nodes");
 selnflux.read(in);
  switch (selnflux.st)
    {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "nodflux_unkn");
      if(Tp->gradcomp == 0)
	Tp->gradcomp = 1;
      if(Tp->fluxcomp == 0)
	Tp->fluxcomp = 1;
      if(Tp->gradaver == 0)
      {
        Tp->gradpos  = 2;
	Tp->gradaver = 1;
      }
      if(Tp->fluxaver == 0)
      {
        Tp->fluxpos  = 2;
	Tp->fluxaver = 1;
      }
      selflux = new sel[selnflux.n];
      for(i=0; i<selnflux.n; i++)
        selflux[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
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
      xfscanf(in, "%k", "nodother_comp");
      if(Tp->othercomp == 0)
	Tp->othercomp = 1;
      if(Tp->otheraver == 0)
      {
        Tp->otherpos  = 2;
	Tp->otheraver = 1;
      }
      seloth = new sel[selnoth.n];
      for (i=0; i < selnoth.n; i++)
        seloth[i].read(in);
      break; 
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 5;
  }

  // eq_other values
  xfscanf(in, "%k", "eqother_nodes");
  selneqoth.read(in);
  switch (selneqoth.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "nodeqother_comp");
      if(Tp->eqothercomp == 0)
	Tp->eqothercomp = 1;
      if(Tp->eqotheraver == 0)
      {
        Tp->eqotherpos  = 2;
	Tp->eqotheraver = 1;
      }
      seleqoth = new sel[selneqoth.n];
      for (i=0; i < selneqoth.n; i++)
        seleqoth[i].read(in);
      break; 
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 6;
  }

  // 
  xfscanf(in, "%k", "fluxres_nodes");
  selnforcet.read(in);
  switch (selnforcet.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "flux_comp");
      selforcet = new sel[selnforcet.n];
      for (i=0; i<selnforcet.n; i++)
        selforcet[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 2;
  }



  return 0;
}



/**
  Function prints data with description for output of nodal values to the text file.

  @param out - pointer to the opened text file
*/
void nodeoutgt::print(FILE *out)
{
  long i;

  // step and loadcases
  dstep.print(out);
  if (dstep.st == sel_no)
    return;
  fprintf(out, "\n");

  // unknowns
  selnunkn.print(out);
  switch(selnunkn.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selnunkn.n; i++)
        selunkn[i].print(out);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
  }


  // gradients
  selngrad.print(out);
  switch (selngrad.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selngrad.n; i++)
        selgrad[i].print(out);
      fprintf(out, "\n");
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
  }

  // fluxes
  selnflux.print(out);
  switch (selnflux.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selnflux.n; i++)
        selflux[i].print(out);
      fprintf(out, "\n");
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
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
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
  }

  // eq_other values
  selneqoth.print(out);
  switch (selneqoth.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selneqoth.n; i++)
        seleqoth[i].print(out);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
  }

  // 
  selnforcet.print(out);
  switch (selnforcet.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for (i=0; i<selnunkn.n; i++)
        selforcet[i].print(out);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function nodeoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
  }

}

/**
  Function prints required output values for selected nodes and for given load case and step
  to the output grahics file.
  
  @param out  - pointer to the opened grahics file
  @param lcid - load case id = 0
  @param desclcid - load case description
  @param gf - graphics format
  @param ifor - vector of nodal forces (fluxes)

 */
void nodeoutgt::print_graphics(FILE *out, long lcid, const char *desclcid, graphftt gf, double *ifor)
{
  long i, j;

  if (gf == grftt_open_dx)
    return;
  
  if (selnunkn.st != sel_no)
  {
    if ((gf == grftt_gid) || (gf == grftt_gid_vtk))
      write_gid_unkn(out, lcid, desclcid);
  }
  
  if (selnflux.st != sel_no)
  {
    for(i=0; i<Tp->ntm; i++)
    {
      for (j=0; j<selnflux.n; j++)
      {
        if (selflux[j].presence_id(i))
        {
          if ((gf == grftt_gid) || (gf == grftt_gid_vtk)) //flux is vector
            write_gid_nodvectort(out, flux, lcid, i, desclcid);
          break;
        }
      }
    }
  }
  
  
  if (selngrad.st != sel_no)
  {
    for(i=0; i<Tp->ntm; i++)
    {
      for (j=0; j<selngrad.n; j++)
      {
        if (selgrad[j].presence_id(i))
        {
          if ((gf == grftt_gid) || (gf == grftt_gid_vtk)) //gradient is vector
            write_gid_nodvectort(out, grad, lcid, i, desclcid);
          break;
        }
      }
    }
  }
  
  
  if (selnoth.st != sel_no)
  {
    for(i=0; i<Tt->nodes[0].ncompother; i++)
    {
      for (j=0; j<selnoth.n; j++)
      {
        if (seloth[j].presence_id(i))
        {
          if ((gf == grftt_gid) || (gf == grftt_gid_vtk)) //eqother is scalar
            write_gid_nodscalart(out, othert, lcid, i, desclcid);
          break;
        }
      }
    }
  }
  
  
  if (selneqoth.st != sel_no)
  {
    for(i=0; i<Tt->nodes[0].ncompeqother; i++)
    {
      for (j=0; j<selneqoth.n; j++)
      {
        if (seleqoth[j].presence_id(i))
        {
          if ((gf == grftt_gid)  || (gf == grftt_gid_vtk)) //eqother is scalar
            write_gid_nodscalart(out, eqothert, lcid, i, desclcid);
          break;
        }
      }
    }
  }
  
  if (selnforcet.st != sel_no)
  {
    if ((gf == grftt_gid)  || (gf == grftt_gid_vtk))   
      write_gid_nforcest(out, lcid, desclcid, ifor);
    else
      write_nforcest(out, lcid, desclcid, ifor);
  }
}





/**
  Function prints required output values for selected nodes and for given load case and step
  to the several output grahics GiD files named by printed quantity component.
  
  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format
  @param ifor - vector of nodal forces (fluxes)
  @param hf - flag for writing of GiD result file header/footer,
              hf_off(default value) = no header/footer is written (standard time step)
              header                = header and/or group beginning is written (inital time step)
              footer                = footer of group is written (final time step)
*/
void nodeoutgt::print_graphics(char *outfn, const char *mode, long lcid, const char *desclcid, graphftt gf, double *ifor, hflagsw hf)
{
  long i, j;
  char fname[FNAMELEN+70];
  FILE *out;
  
  if ((gf != grftt_gid_sep) && (gf != grftt_gidsep_vtk))
  {
    print_err("Invalid graphics format %d is required", __FILE__, __LINE__, __func__, int(gf));
    return;
  }
  
  if (selnunkn.st != sel_no)
  {
    sprintf(fname, "%s.unkn.flavia.res", outfn);
    out = fopen(fname, mode);
    if (out == NULL)
    {
      fprintf(stderr, "\n\nUnable to open graphics file '%s' in function nodeoutt::print_graphics\n", fname);
      fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
      abort();
    }
    switch (hf)
    {
    case header:{
      fseek(out, 0, SEEK_END); // MS Visual C++ requires that
      if (ftell(out) == 0)
        {
          fprintf(out, "GiD Post Results File 1.0\n");
          export_gid_gauss_ptt(out);
        }
      if (Adat)
	fprintf(out, "\nOnGroup adapt_step_%ld\n\n", Adat->istep+1);
      break;
    }
    case hf_off:
      write_gid_unkn(out, lcid, desclcid);
      break;
    case footer:
      if (Adat)
	fprintf(out, "End OnGroup\n\n");
      break;
    default:
      {
        print_err("unknown type of header/footer flag is required", __FILE__, __LINE__, __func__);
        abort();
      }
    }
    fclose(out);
  }
  if (selngrad.st != sel_no)
  {
    for(i=0; i<Tp->ntm; i++)
    {
      for (j=0; j<selngrad.n; j++)
      {
        if (selgrad[j].presence_id(i))
        {
          sprintf(fname, "%s.nodal_grad%ld.flavia.res", outfn, i+1);
          out = fopen(fname, mode);
          if (out == NULL)
          {
            fprintf(stderr, "\n\nUnable to open graphics file '%s' in function nodeoutt::print_graphics\n", fname);
            fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
            abort();
          }
          if (hf == footer)
          {
            if (Adat)
              fprintf(out, "End OnGroup\n\n");
            fclose(out);
            break;
          }
          if (hf == header)
          {
            fseek(out, 0, SEEK_END); // MS Visual C++ requires that
            if (ftell(out) == 0){
              fprintf(out, "GiD Post Results File 1.0\n");
              export_gid_gauss_ptt(out);
            }
            if (Adat)
              fprintf(out, "\nOnGroup adapt_step_%ld\n\n", Adat->istep+1);
          }
          else
            write_gid_nodvectort(out, grad, lcid, i, desclcid);//gradient is vector
          fclose(out);
          break;
        }
      }
    }
  }
  if (selnflux.st != sel_no)
  {
    for(i=0; i<Tp->ntm; i++)
    {
      for (j=0; j<selnflux.n; j++)
      {
        if (selflux[j].presence_id(i))
        {
          sprintf(fname, "%s.nodal_flux%ld.flavia.res", outfn, i+1);
          out = fopen(fname, mode);
          if (out == NULL)
          {
            fprintf(stderr, "\n\nUnable to open graphics file '%s' in function nodeoutt::print_graphics\n", fname);
            fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
            abort();
          }
          if (hf == footer)
          {
            if (Adat)
              fprintf(out, "End OnGroup\n\n");
            fclose(out);
            break;
          }
          if (hf == header)
          {
            fseek(out, 0, SEEK_END); // MS Visual C++ requires that
            if (ftell(out) == 0){
              fprintf(out, "GiD Post Results File 1.0\n");
              export_gid_gauss_ptt(out);
            }
            if (Adat)
              fprintf(out, "\nOnGroup adapt_step_%ld\n\n", Adat->istep+1);
          }
          else
            write_gid_nodvectort(out, flux, lcid, i, desclcid);//flux is vector
          fclose(out);
          break;
        }
      }
    }
  }
  if (selnoth.st != sel_no)
  {
    for(i=0; i<Tt->nodes[0].ncompother; i++)
    {
      for (j=0; j<selnoth.n; j++)
      {
        if (seloth[j].presence_id(i))
        {
          sprintf(fname, "%s.nodal_other%ld.flavia.res", outfn, i+1);
          out = fopen(fname, mode);
          if (out == NULL)
          {
            fprintf(stderr, "\n\nUnable to open graphics file '%s' in function nodeoutt::print_graphics\n", fname);
            fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
            abort();
          }
          if (hf == footer)
          {
            if (Adat)
              fprintf(out, "End OnGroup\n\n");
            fclose(out);
            break;
          }
          if (hf == header)
          {
            fseek(out, 0, SEEK_END); // MS Visual C++ requires that
            if (ftell(out) == 0){
              fprintf(out, "GiD Post Results File 1.0\n");
              export_gid_gauss_ptt(out);
            }
            if (Adat)
              fprintf(out, "\nOnGroup adapt_step_%ld\n\n", Adat->istep+1);
          }
          else
            write_gid_nodscalart(out, othert, lcid, i, desclcid);//other is scalar
          fclose(out);
          break;
        }
      }
    }
  }
  if (selneqoth.st != sel_no)
  {
    for(i=0; i<Tt->nodes[0].ncompeqother; i++)
    {
      for (j=0; j<selneqoth.n; j++)
      {
        if (seleqoth[j].presence_id(i))
        {
          sprintf(fname, "%s.nodal_eqother%ld.flavia.res", outfn, i+1);
          out = fopen(fname, mode);
          if (out == NULL)
          {
            fprintf(stderr, "\n\nUnable to open graphics file '%s' in function nodeoutt::print_graphics\n", fname);
            fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
            abort();
          }
          if (hf == footer)
          {
            if (Adat)
              fprintf(out, "End OnGroup\n\n");
            fclose(out);
            break;
          }
          if (hf == header)
          {
            fseek(out, 0, SEEK_END); // MS Visual C++ requires that
            if (ftell(out) == 0){
              fprintf(out, "GiD Post Results File 1.0\n");
              export_gid_gauss_ptt(out);
            }
            if (Adat)
              fprintf(out, "\nOnGroup adapt_step_%ld\n\n", Adat->istep+1);
          }
          else
            write_gid_nodscalart(out, eqothert, lcid, i, desclcid);//eqother is scalar
          fclose(out);
          break;
        }
      }
    }
  }
  if (selnforcet.st != sel_no)
  {
    sprintf(fname, "%s.flux.res", outfn);
    out = fopen(fname, mode);
    if (out == NULL)
    {
      fprintf(stderr, "\n\nUnable to open graphics file '%s' in function nodeoutgt::print_graphics\n", fname);
      fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
      abort();
    }
    switch (hf)
    {
      case header:
        fseek(out, 0, SEEK_END); // MS Visual C++ requires that
        if (ftell(out) == 0)
        {
          fprintf(out, "GiD Post Results File 1.0\n");
          export_gid_gauss_ptt(out);
        }
        if (Adat)
          fprintf(out, "\nOnGroup adapt_step_%ld\n\n", Adat->istep+1);
        break;
      case hf_off:
        write_gid_nforcest(out, lcid, desclcid, ifor);
        break;
      case footer:
        if (Adat)
          fprintf(out, "End OnGroup\n\n");
        break;
      default:
      {
        print_err("unknown type of header/footer flag is required", __FILE__, __LINE__, __func__);
        abort();
      }
    }
    fclose(out);
  }
}



/**
  Converts selections given by property id to list or range type
  in all output collections according nodal properties defined in the topology top.

  @param top - mesh topology with defined properties of nodes and elements

  @return The function does not return anything but it changes internal 
          representation selections of all nodes that were givne by property id.

  Created by Tomas Koudelka, 08.2016
*/
void nodeoutgt::conv_sel_prop(siftop *top)
{
  sel *newsel;
  long *newtrans;

  if (selnunkn.st == sel_prop)
  {
    selnunkn.conv_selprop(top, gnod, selunkn, newsel, NULL, newtrans);
    delete [] selunkn;
    selunkn = newsel;
  }

  if (selngrad.st == sel_prop)
  {
    selngrad.conv_selprop(top, gnod, selgrad, newsel, NULL, newtrans);
    delete [] selgrad;
    selgrad = newsel;
  }

  if (selnflux.st == sel_prop)
  {
    selnflux.conv_selprop(top, gnod, selflux, newsel, NULL, newtrans);
    delete [] selflux;
    selflux = newsel;
  }

  if (selnoth.st == sel_prop)
  {
    selnoth.conv_selprop(top, gnod, seloth, newsel, NULL, newtrans);
    delete [] seloth;
    seloth = newsel;
  }

  if (selnforcet.st == sel_prop)
  {
    selnforcet.conv_selprop(top, gnod, selforcet, newsel, NULL, newtrans);
    delete [] selforcet;
    selforcet = newsel;
  }
}





/**
  Constructor initializes data to zero values
*/
elemoutt::elemoutt()
{
  selgrad = selflux = seloth = seleqoth = NULL;
}



/**
  Destructor deallocates used memory
*/
elemoutt::~elemoutt()
{
  delete [] selgrad;
  delete [] selflux;
  delete [] seloth;
  delete [] seleqoth;
}



/**
  Function reads data with description for output of element values from the text file
   - not completed

  @param in - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - error reading step 
  @retval 2 - error reading gradient selection
  @retval 3 - error reading flux selection
  @retval 4 - error reading other values selection
  @retval 5 - error reading eqother values selection
 */
long elemoutt::read(XFILE *in)
{
  long i;
  // step
  xfscanf(in, "%k", "sel_elemstep");
  dstep.read(in);
  if (dstep.st == sel_no)
    return 0;

  // gradients
  xfscanf(in, "%k", "grad_elems");
  selegrad.read(in);
  switch (selegrad.st)
    {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "elemgrad_unkn");
      Tp->gradcomp = 1;
      Tp->gradpos = 1;
      selgrad = new sel[selegrad.n];
      for(i=0; i<selegrad.n; i++)
        selgrad[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 2;
    }
  
  
  // fluxes
  xfscanf(in, "%k", "flux_elems");
  seleflux.read(in);
  switch (seleflux.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "elemflux_comp");
      Tp->gradcomp = 1;
      Tp->gradpos = 1;
      Tp->fluxcomp = 1;
      Tp->fluxpos = 1;
      selflux = new sel[seleflux.n];
      for(i=0; i<seleflux.n; i++)
        selflux[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
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
      xfscanf(in, "%k", "elemother_comp");
      Tp->othercomp = 1;
      Tp->otherpos = 1;
      seloth = new sel[seleoth.n]; 
      for (i=0; i<seleoth.n; i++)
        seloth[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 4;
  }

  // eq_other values
  xfscanf(in, "%k", "eqother_elems");
  seleeqoth.read(in);
  switch (seleeqoth.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "elemeqother_comp");
      Tp->eqothercomp = 1;
      Tp->eqotherpos = 1;
      seleqoth = new sel[seleeqoth.n]; 
      for (i=0; i<seleeqoth.n; i++)
        seleqoth[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 5;
  }
  return 0;
}



/**
  Function prints data with description for output of element values to the text file
 - not completed

  @param out - pointer to the opened text file
*/
void elemoutt::print(FILE *out)
{
  long i;

  // step and loadcases
  dstep.print(out);
  if (dstep.st == sel_no)
    return;

  // gradients
  selegrad.print(out);
  switch (selegrad.st)
    {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selegrad.n; i++)
	selgrad[i].print(out);
      fprintf(out, "\n");
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
    }
  
  // fluxes
  seleflux.print(out);
  switch (seleflux.st)
    {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<seleflux.n; i++)
        selflux[i].print(out);
      fprintf(out, "\n");
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
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
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
  }
  
  // eq_other values
  seleeqoth.print(out);
  switch (seleeqoth.st)
    {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<seleeqoth.n; i++)
        seleqoth[i].print(out);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
  }
}


/**
  Function prints required output values for selected elements and for given load case
  to the output text file  - not completed
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

 */
void elemoutt::print_out(FILE *out, long lcid)
{
  if (selegrad.st != sel_no)
  {
    fprintf(out, "** Element gradients :\n\n");
    print_grad(out, lcid);
  }
  if (seleflux.st != sel_no)
  {
    fprintf(out, "** Element fluxes :\n\n");
    print_flux(out, lcid);
  }
  if (seleoth.st != sel_no)
  {
    fprintf(out, "** Element other values :\n\n");
    print_other(out);
  }
  if (seleeqoth.st != sel_no)
  {
    fprintf(out, "** Element eq_other values :\n\n");
    print_eqother(out);
  }
}


/**
  Function prints required grads for selected elements and for given load case
  to the output text file - not completed
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

 */
void elemoutt::print_grad(FILE *out, long /*lcid*/)
{
  long i, ii, j, k, ncomp, ipp, tnipe;
  
  for (i=0; i<Tt->ne; i++)
    {
      if (selegrad.presence_id(i))
	{
	  ipp = Tt->elements[i].ipp[0][0];
	  tnipe = Tt->give_tnip(i);
	  fprintf(out, "   Element %7ld, integration points %ld - %ld\n  ", i+1, ipp+1, ipp+tnipe);
	  ncomp = Tt->give_ncomp (i);
	  for (j=0; j<tnipe; j++)
	    {
	      for (ii=0; ii<Tp->ntm; ii++)
		for (k=0; k<ncomp; k++)
		  {
		    if (selegrad.presence_id(selgrad, i, ii))
		      fprintf(out, "   grad_%ld_%ld=% .*e", ii+1, k+1, prgrad, Tm->ip[ipp+j].grad[ii][k]);
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
   Function prints required fluxes for selected elements and for given load case
   to the output text file - not completed
   
   @param out  - pointer to the opened text file
   @param lcid - load case id
   
*/
void elemoutt::print_flux(FILE *out, long /*lcid*/)
{
  long i, ii, j, k, ncomp, ipp, tnipe;
  
  for (i=0; i<Tt->ne; i++)
    {
      if (seleflux.presence_id(i))
	{
	  ipp = Tt->elements[i].ipp[0][0];
	  tnipe = Tt->give_tnip(i);
	  fprintf(out, "   Element %7ld, integration points %ld - %ld\n  ", i+1, ipp+1, ipp+tnipe);
	  ncomp = Tt->give_ncomp (i);
	  for (j=0; j<tnipe; j++)
	    {
	      for (ii=0; ii<Tp->ntm; ii++)
		for (k=0; k<ncomp; k++)
		  {
		    if (seleflux.presence_id(selflux, i, ii))
		      fprintf(out, "   flux_%ld_%ld=% .*e", ii+1, k+1, prflux, Tm->ip[ipp+j].fluxes[ii][k]);
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
   to the output text file - not completed
   
   @param out  - pointer to the opened text file
   
*/
void elemoutt::print_other(FILE *out)
{
  long i, j, k, ncomp, ipp, tnipe;
  
  for (i=0; i<Tt->ne; i++)
    {
      if (seleoth.presence_id(i))
	{
	  ipp = Tt->elements[i].ipp[0][0];
	  tnipe = Tt->give_tnip(i);
	  fprintf(out, "   Element %7ld, integration points %ld - %ld\n  ", i+1, ipp+1, ipp+tnipe);
	  ncomp = Tm->givencompother();
	  for (j=0; j<tnipe; j++)
	    {
	      for (k=0; k<ncomp; k++)
		{
		  if (seleoth.presence_id(seloth, i, k)){
		    fprintf(out, "\n");
		    fprintf(out, "   other_%ld ", k+1);
		    Tm->give_othervalue_name(out,0,k);
		    fprintf(out, " =% .*e", proth, Tm->ip[ipp+j].other[k]);
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
   to the output text file - not completed
   
   @param out  - pointer to the opened text file
   
*/
void elemoutt::print_eqother(FILE *out)
{
  long i, j, k, ncomp, ipp, tnipe;
  
  for (i=0; i<Tt->ne; i++)
    {
      if (seleeqoth.presence_id(i))
	{
	  ipp = Tt->elements[i].ipp[0][0];
	  tnipe = Tt->give_tnip(i);
	  fprintf(out, "   Element %7ld, integration points %ld - %ld\n  ", i+1, ipp+1, ipp+tnipe);
	  ncomp = Tm->ip[ipp].ncompeqother;
	  for (j=0; j<tnipe; j++)
	    {
	      for (k=0; k<ncomp; k++)
		{
		  if (seleeqoth.presence_id(seleqoth, i, k)){
		    fprintf(out, "\n");
		    fprintf(out, "   eqother_%ld ", k+1);
		    Tm->give_eqothervalue_name(out,0,k);
		    fprintf(out, " =% .*e", preqoth, Tm->ip[ipp+j].eqother[k]);
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
  Converts selections given by property id to list or range type
  in all output collections according element properties defined in the topology top.

  @param top - mesh topology with defined properties of elements

  @return The function does not return anything but it changes internal 
          representation selections of all elements that were givne by property id.

  Created by Tomas Koudelka, 08.2016
*/
void elemoutt::conv_sel_prop(siftop *top)
{
  sel *newsel;
  long *newtrans;

  if (selegrad.st == sel_prop)
  {
    selegrad.conv_selprop(top, gelem, selgrad, newsel, NULL, newtrans);
    delete [] selgrad;
    selgrad = newsel;
  }

  if (seleflux.st == sel_prop)
  {
    seleflux.conv_selprop(top, gelem, selflux, newsel, NULL, newtrans);
    delete [] selflux;
    selflux = newsel;
  }

  if (seleoth.st == sel_prop)
  {
    seleoth.conv_selprop(top, gelem, seloth, newsel, NULL, newtrans);
    delete [] seloth;
    seloth = newsel;
  }

  if (seleeqoth.st == sel_prop)
  {
    seleeqoth.conv_selprop(top, gelem, seleqoth, newsel, NULL, newtrans);
    delete [] seleqoth;
    seleqoth = newsel;
  }
}






/**
  Constructor initializes data to zero values
*/
elemoutgt::elemoutgt()
{
  selgrad = selflux = seloth = seleqoth = NULL;
  ide1 = 1;
}



/**
  Destructor deallocates used memory
*/
elemoutgt::~elemoutgt()
{
  delete [] selgrad;
  delete [] selflux;
  delete [] seloth;
  delete [] seleqoth;
}



/**
  Function reads data with description for output of element values from the text file
   - not completed

  @param in - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - error reading step 
  @retval 2 - error reading gradient selection
  @retval 3 - error reading flux selection
  @retval 4 - error reading other values selection
  @retval 5 - error reading eqother values selection
 */
long elemoutgt::read(XFILE *in)
{
  long i;
  // step
  xfscanf(in, "%k", "sel_elemstep");
  dstep.read(in);
  if (dstep.st == sel_no)
    return 0;

  // gradients
  xfscanf(in, "%k", "grad_elems");
  selegrad.read(in);
  switch (selegrad.st)
    {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "elemgrad_unkn");
      Tp->gradcomp = 1;
      Tp->gradpos = 1;
      selgrad = new sel[selegrad.n];
      for(i=0; i<selegrad.n; i++)
        selgrad[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 2;
    }
  
  
  // fluxes
  xfscanf(in, "%k", "flux_elems");
  seleflux.read(in);
  switch (seleflux.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "elemflux_unkn");
      Tp->gradcomp = 1;
      Tp->fluxcomp = 1;
      Tp->gradpos = 1;
      Tp->fluxpos = 1;
      selflux = new sel[seleflux.n];
      for(i=0; i<seleflux.n; i++)
        selflux[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
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
      xfscanf(in, "%k", "elemother_comp");
      Tp->othercomp = 1;
      Tp->otherpos = 1;
      seloth = new sel[seleoth.n]; 
      for (i=0; i<seleoth.n; i++)
        seloth[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 4;
  }

  // eq_other values
  xfscanf(in, "%k", "eqother_elems");
  seleeqoth.read(in);
  switch (seleeqoth.st)
  {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      xfscanf(in, "%k", "elemeqother_comp");
      Tp->eqothercomp = 1;
      Tp->eqotherpos = 1;
      seleqoth = new sel[seleeqoth.n]; 
      for (i=0; i<seleeqoth.n; i++)
        seleqoth[i].read(in);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 5;
  }
  return 0;
}



/**
  Function prints data with description for output of element values to the text file
 - not completed

  @param out - pointer to the opened text file
*/
void elemoutgt::print(FILE *out)
{
  long i;

  // step and loadcases
  dstep.print(out);
  if (dstep.st == sel_no)
    return;

  // gradients
  selegrad.print(out);
  switch (selegrad.st)
    {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<selegrad.n; i++)
	selgrad[i].print(out);
      fprintf(out, "\n");
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
    }
  
  // fluxes
  seleflux.print(out);
  switch (seleflux.st)
    {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<seleflux.n; i++)
        selflux[i].print(out);
      fprintf(out, "\n");
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
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
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
  }
  
  // eq_other values
  seleeqoth.print(out);
  switch (seleeqoth.st)
    {
    case sel_no:
      break;
    case sel_all:
    case sel_range:
    case sel_list:
      for(i=0; i<seleeqoth.n; i++)
        seleqoth[i].print(out);
      break;
    default:
      fprintf(stderr, "\n\nUnknown type of selection is required in function elemoutt::print\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
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
void elemoutgt::print_graphics(FILE *out, long lcid, const char *desclcid, graphftt gf, long idelem1)
{
  ide1 = idelem1;
  if (gf == grftt_open_dx)
    return;
    
  print_gr_grad_scal(out, desclcid, gf);
  print_gr_flux_scal(out, desclcid, gf);
  print_gr_oth_scal(out, lcid, desclcid, gf);
  print_gr_eqoth_scal(out, lcid, desclcid, gf);

  print_gr_grad_vec(out, lcid, desclcid, gf);
  print_gr_flux_vec(out, lcid, desclcid, gf);
  print_gr_oth_vec(out, lcid, desclcid, gf);
  print_gr_eqoth_vec(out, lcid, desclcid, gf);
}



/** 
  Function prints values of selected gradients as scalars for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgt::print_gr_grad_scal(FILE *out, const char *desclcid, graphftt gf)
{
  long i,ii,j;
 
  if (selegrad.st != sel_no)
  {
    for(ii=0; ii<Tp->ntm; ii++)
    {
      for(i=0; i<Tp->gdim; i++)
      {
	for (j=0; j<selegrad.n; j++)
	{
          if ((selgrad[j].st == sel_mtx) || (selgrad[j].st == sel_range_mtx) ||
              (selgrad[j].st == sel_range_vec) || (selgrad[j].st == sel_vec))
            continue;
          if (selgrad[j].presence_id(ii))
          {
            if ((gf == grftt_gid) || (gf == grftt_gid_vtk))
              write_gid_elemscalart(out, grad, ii, i, desclcid);
            break;
          }
        }
      }
    }
  }
}



/** 
  Function prints values of selected fluxes as scalars for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgt::print_gr_flux_scal(FILE *out, const char *desclcid, graphftt gf)
{
  long i,ii,j;
 
  if (seleflux.st != sel_no)
  {
    for(ii=0; ii<Tp->ntm; ii++)
    {
      for(i=0; i<Tp->gdim; i++)
      {
        for (j=0; j<seleflux.n; j++)
        {
          if ((selflux[j].st == sel_mtx) || (selflux[j].st == sel_range_mtx) ||
              (selflux[j].st == sel_range_vec) || (selflux[j].st == sel_vec))
            continue;
          if (selflux[j].presence_id(ii))
          {
            if ((gf == grftt_gid) || (gf == grftt_gid_vtk))
              write_gid_elemscalart(out, flux, ii, i, desclcid);
            break;
          }
        }
      }
    }
  }
}



/** 
  Function prints selected other values as scalars for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgt::print_gr_oth_scal(FILE *out, long lcid, const char *desclcid, graphftt gf)
{
  long i,j;
 
  if (seleoth.st != sel_no)
  {
    for(i=0; i<Tm->ip[0].ncompother; i++)
    {
      for (j=0; j<seleoth.n; j++)
      {
        if ((seloth[j].st == sel_mtx) || (seloth[j].st == sel_range_mtx) ||
            (seloth[j].st == sel_range_vec) || (seloth[j].st == sel_vec))
          continue;
        if (seloth[j].presence_id(i))
        {
          if ((gf == grftt_gid) || (gf == grftt_gid_vtk))
            write_gid_elemscalart(out, othert, lcid, i, desclcid);
          break;
        }
      }
    }
  }
}



/** 
  Function prints selected eqother values as scalars for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgt::print_gr_eqoth_scal(FILE *out, long lcid, const char *desclcid, graphftt gf)
{
  long i,j;
 
  if (seleeqoth.st != sel_no)
  {
    for(i=0; i<Tm->ip[0].ncompeqother; i++)
    {
      for (j=0; j<seleeqoth.n; j++)
      {
        if ((seleqoth[j].st == sel_mtx) || (seleqoth[j].st == sel_range_mtx) ||
            (seleqoth[j].st == sel_range_vec) || (seleqoth[j].st == sel_vec))
          continue;
        if (seleqoth[j].presence_id(i))
        {
          if ((gf == grftt_gid) || (gf == grftt_gid_vtk))
            write_gid_elemscalart(out, eqothert, lcid, i, desclcid);
          break;
        }
      }
    }
  }
}



/** 
  Function prints values of selected gradients as vectors for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param lcid - load case id
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgt::print_gr_grad_vec(FILE *out, long lcid, const char *desclcid, graphftt gf)
{
  long j;

  if (selegrad.st != sel_no)
  {
    for (j=0; j<selegrad.n; j++)
    {
      if (selgrad[j].st != sel_vec)
        continue;
      if ((gf == grftt_gid) || (gf == grftt_gid_vtk))
        write_gid_elemvectort(out, grad, lcid, j, desclcid);
      else
      {
        print_err("required export format of a strain vector is not supported", __FILE__, __LINE__, __func__);
        abort();
      }
    }
  }
}



/** 
  Function prints values of selected fluxes as vectors for selected load case to graphics file.

  @param out  - pointer to the opened grahics file
  @param desclcid - load case description
  @param gf - graphics format

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgt::print_gr_flux_vec(FILE *out, long lcid, const char *desclcid, graphftt gf)
{
  long j;

  if (seleflux.st != sel_no)
  {
    for (j=0; j<seleflux.n; j++)
    {
      if (selflux[j].st != sel_vec)
        continue;
      if ((gf == grftt_gid) || (gf == grftt_gid_vtk))
        write_gid_elemvectort(out, flux, lcid, j, desclcid);
      else
      {
        print_err("required export format of a strain vector is not supported", __FILE__, __LINE__, __func__);
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
void elemoutgt::print_gr_oth_vec(FILE *out, long lcid, const char *desclcid, graphftt gf)
{
  long j;

  if (seleoth.st != sel_no)
  {
    for (j=0; j<seleoth.n; j++)
    {
      if (seloth[j].st != sel_range_vec)
        continue;
      if ((gf == grftt_gid) || (gf == grftt_gid_vtk))
        write_gid_elemvectort(out, othert, lcid, j, desclcid);
      else
      {
        print_err("unsupported export format of an other values vector is required",__FILE__, __LINE__, __func__);
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
void elemoutgt::print_gr_eqoth_vec(FILE *out, long lcid, const char *desclcid, graphftt gf)
{
  long j;

  if (seleeqoth.st != sel_no)
  {
    for (j=0; j<seleeqoth.n; j++)
    {
      if (seleqoth[j].st != sel_range_vec)
        continue;
      if ((gf == grftt_gid) || (gf == grftt_gid_vtk))
        write_gid_elemvectort(out, eqothert, lcid, j, desclcid);
      else
      {
        print_err("unsupported export format of an eqother values vector is required",__FILE__, __LINE__, __func__);
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
  @param hf - flag for writing of GiD result file header/footer,
              hf_off(default value) = no header/footer is written (standard time step)
              header                = header and/or group beginning is written (inital time step)
              footer                = footer of group is written (final time step)
              
*/
void elemoutgt::print_graphics(const char *outfn, const char *mode, long lcid, const char *desclcid, 
                               graphftt gf, hflagsw hf, long idelem1)
{
  ide1 = idelem1;
  if ((gf != grftt_gid_sep) && (gf != grftt_gidsep_vtk))
  {
    print_err("invalid graphics format is required", __FILE__, __LINE__, __func__);
    return;
  }
  print_gr_grad_scal(outfn, mode, desclcid, hf);
  print_gr_flux_scal(outfn, mode, desclcid, hf);
  print_gr_oth_scal(outfn, mode, lcid, desclcid, hf);
  print_gr_eqoth_scal(outfn, mode, lcid, desclcid, hf);

  print_gr_grad_vec(outfn, mode, lcid, desclcid, hf);
  print_gr_flux_vec(outfn, mode, lcid, desclcid, hf);
  print_gr_oth_vec(outfn, mode, lcid, desclcid, hf);
  print_gr_eqoth_vec(outfn, mode, lcid, desclcid, hf);

  return;
}



/** 
  Function prints values of selected gradients as scalars for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param desclcid - load case description
  @param hf - flag for writing of GiD result file header/footer,
              hf_off(default value) = no header/footer is written (standard time step)
              header                = header and/or group beginning is written (inital time step)
              footer                = footer of group is written (final time step)

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgt::print_gr_grad_scal(const char *outfn, const char *mode, const char *desclcid, hflagsw hf)
{
  long i,ii,j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (selegrad.st != sel_no)
  {
    for(ii=0; ii<Tp->ntm; ii++)
    {
      for(i=0; i<Tm->ip[0].ncompgrad; i++)
      {
        for (j=0; j<selegrad.n; j++)
        {
          if (selgrad[j].presence_id(ii))
          {
            sprintf(fname, "%s.elem_grad%ld_%ld.flavia.res", outfn, ii+1, i+1);
            out = fopen(fname, mode);
            if (out == NULL)
            {
              fprintf(stderr, "\n\nUnable to open graphics file '%s' in function elemoutt::print_graphics\n", fname);
              fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
              abort();
            }
            if (hf == footer)
            {
              if (Adat)
                fprintf(out, "End OnGroup\n\n");
              fclose(out);
              break;
            }
            if (hf == header)
            {
              fseek(out, 0, SEEK_END); // MS Visual C++ requires that
              if (ftell(out) == 0)
              {
                fprintf(out, "GiD Post Results File 1.0\n");
                export_gid_gauss_ptt(out);
              }
              if (Adat)
                fprintf(out, "\nOnGroup adapt_step_%ld\n\n", Adat->istep+1);
            }
            else
              write_gid_elemscalart(out, grad, ii, i, desclcid);
            fclose(out);
            break;
          }
        }
      }
    }
  }
}



/** 
  Function prints values of selected fluxes as scalars for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param desclcid - load case description
  @param hf - flag for writing of GiD result file header/footer,
              hf_off(default value) = no header/footer is written (standard time step)
              header                = header and/or group beginning is written (inital time step)
              footer                = footer of group is written (final time step)

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgt::print_gr_flux_scal(const char *outfn, const char *mode, const char *desclcid, hflagsw hf)
{
  long i,ii,j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (seleflux.st != sel_no)
  {
    for(ii=0; ii<Tp->ntm; ii++)
    {
      for(i=0; i<Tm->ip[0].ncompgrad; i++)
      {
        for (j=0; j<seleflux.n; j++)
        {
          if (selflux[j].presence_id(ii))
          {
            sprintf(fname, "%s.elem_flux%ld_%ld.flavia.res", outfn, ii+1, i+1);
            out = fopen(fname, mode);
            if (out == NULL)
            {
              fprintf(stderr, "\n\nUnable to open graphics file '%s' in function elemoutt::print_graphics\n", fname);
              fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
              abort();
            }
            if (hf == footer)
            {
              if (Adat)
                fprintf(out, "End OnGroup\n\n");
              fclose(out);
              break;
            }
            if (hf == header)
            {
              fseek(out, 0, SEEK_END); // MS Visual C++ requires that
              if (ftell(out) == 0)
              {
                fprintf(out, "GiD Post Results File 1.0\n");
                export_gid_gauss_ptt(out);
              }
              if (Adat)
                fprintf(out, "\nOnGroup adapt_step_%ld\n\n", Adat->istep+1);
            }
            else
              write_gid_elemscalart(out, flux, ii, i, desclcid);
            fclose(out);
            break;
          }
        }
      }
    }
  }
}



/** 
  Function prints selected other values as scalars for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description
  @param hf - flag for writing of GiD result file header/footer,
              hf_off(default value) = no header/footer is written (standard time step)
              header                = header and/or group beginning is written (inital time step)
              footer                = footer of group is written (final time step)

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgt::print_gr_oth_scal(const char *outfn, const char *mode, long lcid, const char *desclcid, hflagsw hf)
{
  long i, j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (seleoth.st != sel_no)
  {
    for(i=0; i<Tm->ip[0].ncompother; i++)
    {
      for (j=0; j<seleoth.n; j++)
      {
        if (seloth[j].presence_id(i))
        {
          sprintf(fname, "%s.elem_other%ld.flavia.res", outfn, i+1);
          out = fopen(fname, mode);
          if (out == NULL)
          {
            fprintf(stderr, "\n\nUnable to open graphics file '%s' in function elemoutt::print_graphics\n", fname);
            fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
            abort();
          }
          if (hf == footer)
          {
            if (Adat)
              fprintf(out, "End OnGroup\n\n");
            fclose(out);
            break;
          }
          if (hf == header)
          {
            fseek(out, 0, SEEK_END); // MS Visual C++ requires that
            if (ftell(out) == 0)
            {
              fprintf(out, "GiD Post Results File 1.0\n");
              export_gid_gauss_ptt(out);
            }
            if (Adat)
              fprintf(out, "\nOnGroup adapt_step_%ld\n\n", Adat->istep+1);
          }
          else
            write_gid_elemscalart(out, othert, lcid, i, desclcid);
          fclose(out);
          break;
        }
      }
    }
  }
}



/** 
  Function prints selected eqother values as scalars for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description
  @param hf - flag for writing of GiD result file header/footer,
              hf_off(default value) = no header/footer is written (standard time step)
              header                = header and/or group beginning is written (inital time step)
              footer                = footer of group is written (final time step)

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgt::print_gr_eqoth_scal(const char *outfn, const char *mode, long lcid, const char *desclcid, hflagsw hf)
{
  long i, j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (seleeqoth.st != sel_no)
  {
    for(i=0; i<Tm->ip[0].ncompeqother; i++)
    {
      for (j=0; j<seleeqoth.n; j++)
      {
        if (seleqoth[j].presence_id(i))
        {
          sprintf(fname, "%s.elem_eqother%ld.flavia.res", outfn, i+1);
          out = fopen(fname, mode);
          if (out == NULL)
          {
            fprintf(stderr, "\n\nUnable to open graphics file '%s' in function elemoutt::print_graphics\n", fname);
            fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
            abort();
          }
          if (hf == footer)
          {
            if (Adat)
              fprintf(out, "End OnGroup\n\n");
            fclose(out);
            break;
          }
          if (hf == header)
          {
            fseek(out, 0, SEEK_END); // MS Visual C++ requires that
            if (ftell(out) == 0)
            {
              fprintf(out, "GiD Post Results File 1.0\n");
              export_gid_gauss_ptt(out);
            }
            if (Adat)
              fprintf(out, "\nOnGroup adapt_step_%ld\n\n", Adat->istep+1);
          }
          else
            write_gid_elemscalart(out, eqothert, lcid, i, desclcid);
          fclose(out);
          break;
        }
      }
    }
  }
}



/** 
  Function prints values of selected gradients as vectors for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description
  @param hf - flag for writing of GiD result file header/footer,
              hf_off(default value) = no header/footer is written (standard time step)
              header                = header and/or group beginning is written (inital time step)
              footer                = footer of group is written (final time step)

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgt::print_gr_grad_vec(const char *outfn, const char *mode, long lcid, const char *desclcid, hflagsw hf)
{
  long j;
  char fname[FNAMELEN+70];
  FILE *out;


  if (selegrad.st != sel_no)
  {
    for (j=0; j<selegrad.n; j++)
    {
      if (selgrad[j].st != sel_vec)
        continue;
      sprintf(fname, "%s.elem_grad%ld_v_s%ld.flavia.res", outfn, selgrad[j].id1[0]+1, j+1);
      out = fopen(fname, mode);
      if (out == NULL)
      {
        fprintf(stderr, "\n\nUnable to open graphics file '%s' in function elemoutt::print_graphics\n", fname);
        fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
        abort();
      }
      if (hf == footer)
      {
        if (Adat)
          fprintf(out, "End OnGroup\n\n");
        fclose(out);
        break;
      }
      if (hf == header)
      {
        fseek(out, 0, SEEK_END); // MS Visual C++ requires that
        if (ftell(out) == 0)
        {
          fprintf(out, "GiD Post Results File 1.0\n");
          export_gid_gauss_ptt(out);
        }
        if (Adat)
          fprintf(out, "\nOnGroup adapt_step_%ld\n\n", Adat->istep+1);
      }
      else
        write_gid_elemvectort(out, grad, lcid, j, desclcid);
      fclose(out);
      break;
    }
  }
}



/** 
  Function prints values of selected fluxes as vectors for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description
  @param hf - flag for writing of GiD result file header/footer,
              hf_off(default value) = no header/footer is written (standard time step)
              header                = header and/or group beginning is written (inital time step)
              footer                = footer of group is written (final time step)

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgt::print_gr_flux_vec(const char *outfn, const char *mode, long lcid, const char *desclcid, hflagsw hf)
{
  long j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (seleflux.st != sel_no)
  {
    for (j=0; j<seleflux.n; j++)
    {
      if (selflux[j].st != sel_vec)
        continue;
      sprintf(fname, "%s.elem_flux%ld_v_s%ld.flavia.res", outfn, selflux[j].id1[0]+1, j+1);
      out = fopen(fname, mode);
      if (out == NULL)
      {
        fprintf(stderr, "\n\nUnable to open graphics file '%s' in function elemoutt::print_graphics\n", fname);
        fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
        abort();
      }
      if (hf == footer)
      {
        if (Adat)
          fprintf(out, "End OnGroup\n\n");
        fclose(out);
        break;
      }
      if (hf == header)
      {
        fseek(out, 0, SEEK_END); // MS Visual C++ requires that
        if (ftell(out) == 0)
        {
          fprintf(out, "GiD Post Results File 1.0\n");
          export_gid_gauss_ptt(out);
        }
        if (Adat)
          fprintf(out, "\nOnGroup adapt_step_%ld\n\n", Adat->istep+1);
      }
      else
        write_gid_elemvectort(out, flux, lcid, j, desclcid);
      fclose(out);
      break;
    }
  }
}



/** 
  Function prints selected other values as vectors for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description
  @param hf - flag for writing of GiD result file header/footer,
              hf_off(default value) = no header/footer is written (standard time step)
              header                = header and/or group beginning is written (inital time step)
              footer                = footer of group is written (final time step)

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgt::print_gr_oth_vec(const char *outfn, const char *mode, long lcid, const char *desclcid, hflagsw hf)
{
  long j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (seleoth.st != sel_no)
  {
    for (j=0; j<seleoth.n; j++)
    {
      if (seloth[j].st != sel_range_vec)
        continue;
      sprintf(fname, "%s.elem_other_v%ld-%ld_s%ld.flavia.res", outfn, seloth[j].id1[0]+1, seloth[j].ncomp[0], j+1);
      out = fopen(fname, mode);
      if (out == NULL)
      {
        fprintf(stderr, "\n\nUnable to open graphics file '%s' in function elemoutt::print_graphics\n", fname);
        fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
        abort();
      }
      if (hf == footer)
      {
        if (Adat)
          fprintf(out, "End OnGroup\n\n");
        fclose(out);
        break;
      }
      if (hf == header)
      {
        fseek(out, 0, SEEK_END); // MS Visual C++ requires that
        if (ftell(out) == 0)
        {
          fprintf(out, "GiD Post Results File 1.0\n");
          export_gid_gauss_ptt(out);
        }
        if (Adat)
          fprintf(out, "\nOnGroup adapt_step_%ld\n\n", Adat->istep+1);
      }
      else
        write_gid_elemvectort(out, flux, lcid, j, desclcid);
      fclose(out);
      break;
    }
  }
}



/** 
  Function prints selected eqother values as scalars for selected load case to separated
  graphics files named by printed quantity component.

  @param outfn  - string with file name part
  @param mode - opening mode for graphics files
  @param lcid - load case id
  @param desclcid - load case description
  @param hf - flag for writing of GiD result file header/footer,
              hf_off(default value) = no header/footer is written (standard time step)
              header                = header and/or group beginning is written (inital time step)
              footer                = footer of group is written (final time step)

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void elemoutgt::print_gr_eqoth_vec(const char *outfn, const char *mode, long lcid, const char *desclcid, hflagsw hf)
{
  long j;
  char fname[FNAMELEN+70];
  FILE *out;

  if (seleeqoth.st != sel_no)
  {
    for (j=0; j<seleeqoth.n; j++)
    {
      if (seleqoth[j].st != sel_range_vec)
        continue;
      sprintf(fname, "%s.elem_eqother_v%ld-%ld_s%ld.flavia.res", outfn, seleqoth[j].id1[0]+1, seleqoth[j].ncomp[0], j+1);
      out = fopen(fname, mode);
      if (out == NULL)
      {
        fprintf(stderr, "\n\nUnable to open graphics file '%s' in function elemoutt::print_graphics\n", fname);
        fprintf(stderr, " (file %s, line %d)\n", __FILE__, __LINE__);
        abort();
      }
      if (hf == footer)
      {
        if (Adat)
          fprintf(out, "End OnGroup\n\n");
        fclose(out);
        break;
      }
      if (hf == header)
      {
        fseek(out, 0, SEEK_END); // MS Visual C++ requires that
        if (ftell(out) == 0)
        {
          fprintf(out, "GiD Post Results File 1.0\n");
          export_gid_gauss_ptt(out);
        }
        if (Adat)
          fprintf(out, "\nOnGroup adapt_step_%ld\n\n", Adat->istep+1);
      }
      else
        write_gid_elemvectort(out, flux, lcid, j, desclcid);
      fclose(out);
      break;
    }
  }
}



/**
  Converts selections given by property id to list or range type
  in all output collections according element properties defined in the topology top.

  @param top - mesh topology with defined properties of elements

  @return The function does not return anything but it changes internal 
          representation selections of all elements that were givne by property id.

  Created by Tomas Koudelka, 08.2016
*/
void elemoutgt::conv_sel_prop(siftop *top)
{
  sel *newsel;
  long *newtrans;

  if (selegrad.st == sel_prop)
  {
    selegrad.conv_selprop(top, gelem, selgrad, newsel, NULL, newtrans);
    delete [] selgrad;
    selgrad = newsel;
  }

  if (seleflux.st == sel_prop)
  {
    seleflux.conv_selprop(top, gelem, selflux, newsel, NULL, newtrans);
    delete [] selflux;
    selflux = newsel;
  }

  if (seleoth.st == sel_prop)
  {
    seleoth.conv_selprop(top, gelem, seloth, newsel, NULL, newtrans);
    delete [] seloth;
    seloth = newsel;
  }

  if (seleeqoth.st == sel_prop)
  {
    seleeqoth.conv_selprop(top, gelem, seleqoth, newsel, NULL, newtrans);
    delete [] seleqoth;
    seleqoth = newsel;
  }
}



/************* This below is not finished **************/

/**
  Constructor initializes data to zero values
*/
pointoutt::pointoutt()
{
  npnt = 0;
  ksi = eta = zeta = NULL;
  selgrad = selflux = seloth = seleqoth = NULL;
}



/**
  Destructor deallocates used memory
*/
pointoutt::~pointoutt()
{
  delete [] ksi;
  delete [] eta;
  delete [] zeta;
  delete [] selgrad;
  delete [] selflux;
  delete [] seloth;
  delete [] seleqoth;
}



/**
  Function reads data with description for output of user defined point values from the text file.

  @param in - pointer to the opened text file

  @retval 0 - on success
  @retval 1 - error reading step or number of points
  @retval 2 - error reading natural coordinates
  @retval 3 - error reading sets of points on elements
  @retval 4 - error reading transformation id
 */
long pointoutt::read(XFILE *in)
{
  long i;
  // step and number of points
  dstep.read(in);
  if (dstep.st == sel_no)
    return 0;
  if (xfscanf(in, "%k%ld", "numpoints", &npnt) != 2)
  {
    fprintf(stderr, "\n\nError reading number of points\n");
    fprintf(stderr, " in function pointoutt::read, (file %s, line %d)\n", __FILE__, __LINE__);
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
      fprintf(stderr, "\n\nError reading coordinates of UDP\n");
      fprintf(stderr, " in function pointoutt::read, (file %s, line %d)\n", __FILE__, __LINE__);
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
      fprintf(stderr, "\n\nUnknown type of selection is required in function pointoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return 3;
  }

  sellc.read(in);  
  if (sellc.st == sel_no)
    return 0;

  selgrad = new sel[npnt];
  selflux = new sel[npnt];
  seloth  = new sel[npnt];
  seleqoth  = new sel[npnt];

  for (i=0; i<npnt; i++)
  {
    // gradients
    selgrad[i].read(in);

    // fluxes
    selflux[i].read(in);

    // other values
    seloth[i].read(in);

    // other values
    seleqoth[i].read(in);
    
    // setup of flags
    if (selgrad[i].st != sel_no)
      Tp->gradcomp = 1;
    if (selflux[i].st != sel_no)
      Tp->fluxcomp = 1;
    if (seloth[i].st != sel_no)
      Tp->othercomp = 1;
    if (seleqoth[i].st != sel_no)
      Tp->eqothercomp = 1;
  }
  return 0;
}



/**
  Function prints data with description for output of element values to the text file.

  @param out - pointer to the opened tetx file
*/
void pointoutt::print(FILE *out)
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
      fprintf(stderr, "\n\nUnknown type of selection is required in function pointoutt::read\n");
      fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
      return;
  }

  sellc.print(out);  
  if (sellc.st == sel_no)
    return;

  for (i=0; i<npnt; i++)
  {
    // gradients
    selgrad[i].print(out);
    // fluxes
    selflux[i].print(out);
    // other values
    seloth[i].print(out);
    // eq_other values
    seleqoth[i].print(out);
  }
}



/**
  Function prints required output values for user defined point and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id  

*/
void pointoutt::print_out(FILE *out, long lcid)
{
  long i;
  
  for (i=0; i<npnt; i++)
    {
      if (selgrad[i].st != sel_no)
	fprintf(out, "** UDP gradients :\n\n");
    }
  print_grad(out, lcid);
  
  for (i=0; i<npnt; i++)
    {
      if (selflux[i].st != sel_no)
	fprintf(out, "** UDP fluxes :\n\n");
    }
  print_flux(out, lcid);
  
  for (i=0; i<npnt; i++)
    {
      if (seloth[i].st != sel_no)
	fprintf(out, "** UDP other values :\n\n");
      print_other(out, i);
    }
  for (i=0; i<npnt; i++)
    {
      if (seleqoth[i].st != sel_no)
	fprintf(out, "** UDP eqother values :\n\n");
      print_eqother(out, i);
    }
}



/**
  Function prints required gradients for user defined point and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id

 */
void pointoutt::print_grad(FILE *out, long lcid)
{
  long i,j,k,l,n, ncomp, id;

  for (i=0; i<Tt->ne; i++)
  {
    if (selelem.presence_id(i))
    {
      fprintf(out, "   Element %7ld:\n",i+1);
      ncomp = Tt->give_ncomp(i);
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
          fprintf(stderr, "\n\nUnknown type of selection is required in function pointoutt::print_grad\n");
          fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
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
              if(selgrad[k].presence_id(l))
                fprintf(out, "   eps_%ld=% *e", l+1, prgrad, Tm->grad.ev[i][k][id+l]);//not completed
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
  Function prints required fluxes for user defined point and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param lcid - load case id
  @param pid  - user defined point id

 */
void pointoutt::print_flux(FILE *out, long lcid)
{
  long i,j,k,l,n, ncomp, id;

  for (i=0; i<Tt->ne; i++)
  {
    if (selelem.presence_id(i))
    {
      fprintf(out, "   Element %7ld:\n", i+1);
      ncomp = Tt->give_ncomp(i);
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
          fprintf(stderr, "\n\nUnknown type of selection is required in function pointoutt::print_flux\n");
          fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
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
              if(selgrad[k].presence_id(l))
                fprintf(out, "   sig_%ld=% *e", l+1, prgrad, Tm->flux.ev[i][k][id+l]);//not completed
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

 */
void pointoutt::print_other(FILE */*out*/, long /*pid*/)
{
  long i,j,k,l,n, ncomp;

  for (i=0; i<Tt->ne; i++)
  {
    if (selelem.presence_id(i))
    {
//      fprintf(out, "   Element %7ld:\n"i+1);
      ncomp = Tm->ip[Tt->elements[i].ipp[0][0]].ncompother;
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
          fprintf(stderr, "\n\nUnknown type of selection is required in function pointoutt::print_other\n");
          fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
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
//                fprintf(out, "   other_%ld=% *e", l+1, proth, Mm->other.ev[i][k][l]);//not completed
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
  Function prints required eq_other values for user defined point and for given load case
  to the output text file.
  
  @param out  - pointer to the opened text file
  @param pid  - user defined point id

 */
void pointoutt::print_eqother(FILE */*out*/, long /*pid*/)
{
  long i,j,k,l,n, ncomp;

  for (i=0; i<Tt->ne; i++)
  {
    if (selelem.presence_id(i))
    {
//      fprintf(out, "   Element %7ld:\n"i+1);
      ncomp = Tm->ip[Tt->elements[i].ipp[0][0]].ncompeqother;
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
          fprintf(stderr, "\n\nUnknown type of selection is required in function pointoutt::print_eqother\n");
          fprintf(stderr, " (file %s, line %d)", __FILE__, __LINE__);
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
//              if(seleqoth[k].presence_id(l))
//                fprintf(out, "   eqother_%ld=% *e", l+1, preqoth, Mm->eqother.ev[i][k][l]);//not completed
            }
//            fprintf(out, "\n"); 
          }
        }
      }
//      fprintf(out, "\n");
    }
  } 
}
