#include "outdiagm.h"
#include "iotools.h"
#include "gtopology.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "meshtransfer.h"
#include "element.h"
#include "node.h"
#include "intpoints.h"
#include "globmat.h"



/**
  The constructor initializes data to zero values

  Created by Tomas Koudelka,
*/
outdiagm::outdiagm()
{
  pid = NULL; eid = NULL; ipeid = NULL; nif = NULL;
  pu = NULL;
  ipu = NULL;
  x = y = z = NULL;
}


/**
  The destructor deallocates used memory

  Created by Tomas Koudelka,
*/
outdiagm::~outdiagm()
{
  delete [] nif;
  delete [] pid;
  delete [] eid;
  delete [] ipeid;
  delete [] pu;
  delete [] ipu;
  delete [] x;
  delete [] y;
  delete [] z;
}



/** 
  The function reads data from the file given by the pointer to the 
  opened text file.
   
  @param  in - pointer to the opened text file where the data will be read from

  @retval 0 - on success
  @retval 1 - in case read error   

  Created by Tomas Koudelka,
*/
long outdiagm::read(XFILE *in)
{
  long i;
  if (xfscanf (in, "%k%ld", "numunknowns", &npun) != 2)
  {
    print_err("cannot read number of printed unknowns", __FILE__,__LINE__, __func__);
    return 1;
  }
  xfscanf(in, "%k", "sel_diagstep");
  dstep.read(in);
  if (dstep.st == sel_no)
    return 0;
  nif   = new nodip[npun];
  pid   = new long [npun];
  eid   = new long [npun];
  ipeid = new long [npun];
  pu    = new prunk[npun];
  ipu   = new long [npun];
  x     = new double [npun];
  y     = new double [npun];
  z     = new double [npun];
  memset(nif, 0, sizeof(*nif)*npun);
  memset(pu, 0, sizeof(*pu)*npun);
  memset(ipu, 0, sizeof(*ipu)*npun);
  memset(x, 0, sizeof(*x)*npun);
  memset(y, 0, sizeof(*y)*npun);
  memset(z, 0, sizeof(*z)*npun);

  for (i=0; i<npun; i++)
  {
    pid[i] = eid[i] = ipeid[i] = -1;
    if (xfscanf (in, "%k%m", "point", &nodip_kwdset, (int *)(nif+i)) != 2)
    {
      print_err("cannot read type of point", __FILE__,__LINE__, __func__);
      return 1;
    }
    switch (nif[i])
    {
      case atnode:
        if (xfscanf (in, "%k%ld", "node", pid+i) != 2)
        {
          print_err("cannot read node id", __FILE__,__LINE__, __func__);
          return 1;
        }
        pid[i]--;
        break;
      case atip:
        if (xfscanf (in, "%k %ld %k %ld", "elem", eid+i, "ip", ipeid+i) != 4)
        {
          print_err("cannot read element or ip number", __FILE__,__LINE__, __func__);
          return 1;
        }
        eid[i]--; ipeid[i]--;
        break;
      case atxyz:
        if (xfscanf (in, "%k %le %k %le %k %le", "x", x+i, "y", y+i, "z", z+i) != 6)
        {
          print_err("cannot read point coordinates", __FILE__,__LINE__, __func__);
          return 1;
        }
        if (Gtm->gnodes) // due to preprocessor where pid[i] need not to be set
          pid[i] = Gtm->give_nearest_node(x[i],y[i],z[i]);
        break;
      default:
        print_err("unknown type of point is required", __FILE__, __LINE__, __func__); 
        return 1;
    }
  
    if (xfscanf (in, "%k%m", "quant_type", &prunk_kwdset, (int *)(pu+i)) != 2)
    {
      print_err("cannot read type of printed unknown", __FILE__,__LINE__, __func__);
      return 1;
    }
    if ((pu[i] == pr_stepid) || (pu[i] == pr_appload) || (pu[i] == pr_time))
      continue;
    if (xfscanf (in, "%k%ld", "compid", ipu+i) != 2)
      {
	print_err("cannot read type of printed unknown", __FILE__,__LINE__, __func__);
	return 1;
      }
    ipu[i]--;
    if (nif[i] == atnode)
      {
	switch (pu[i])
	  {
	  case pr_strains:
	    if(Mp->straincomp == 0)
	      Mp->straincomp = 1;
	    if(Mp->strainaver == 0)
	      Mp->strainaver = 1;
	    if(Mp->strainpos == 0)
	      Mp->strainpos = 2;
	    break;
	  case pr_stresses:
	    if(Mp->straincomp == 0)
	      Mp->straincomp = 1;
	    if(Mp->strainaver == 0)
	      Mp->strainaver = 1;
	    if(Mp->strainpos == 0)
	      Mp->strainpos = 2;
	    if(Mp->stresscomp == 0)
	      Mp->stresscomp = 1;
	    if(Mp->stressaver == 0)
	      Mp->stressaver = 1;
	    if(Mp->stresspos == 0)
	      Mp->stresspos = 2;	  
	    break;
	  case pr_react:
	    Mp->reactcomp = 1;
	    break;
	  case pr_other:
	    if(Mp->othercomp == 0)
	      Mp->othercomp = 1;
	    if(Mp->otheraver == 0)
	      Mp->otheraver = 1;
	    if(Mp->otherpos == 0)
	      Mp->otherpos = 2;
	    break;
	  default:
	    break;
	  }
      }
    
    if ((nif[i] == atip) || (nif[i] == atxyz))
      {
	switch (pu[i])
	  {
	  case pr_strains:
	    if(Mp->straincomp == 0)
	      Mp->straincomp = 1;
	    if(Mp->strainaver == 0)
	      Mp->strainaver = 1;
	    if(Mp->strainpos == 0)
	      Mp->strainpos = 1;
	    break;
	  case pr_stresses:
	    if(Mp->straincomp == 0)
	      Mp->straincomp = 1;
	    if(Mp->strainaver == 0)
	      Mp->strainaver = 1;
	    if(Mp->strainpos == 0)
	      Mp->strainpos = 1;
	    if(Mp->stresscomp == 0)
	      Mp->stresscomp = 1;
	    if(Mp->stressaver == 0)
	      Mp->stressaver = 1;
	    if(Mp->stresspos == 0)
	      Mp->stresspos = 1;	  
	    break;
	  case pr_react:
	    Mp->reactcomp = 1;
	    break;
	  case pr_other:
	    if(Mp->othercomp == 0)
	      Mp->othercomp = 1;
	    if(Mp->otheraver == 0)
	      Mp->otheraver = 1;
	    if(Mp->otherpos == 0)
	      Mp->otherpos = 1;
	    break;
	  default:
	    break;
	  }
      }
  }
  return 0;
}



/** 
  The function prints data to the file given by the pointer out.
   
  @param  in - pointer to the opened text file where the data will be read from

  @retval 0 - on success

  Created by Tomas Koudelka,
*/
long outdiagm::print(FILE *out)
{
  long i;
  fprintf(out, "%ld ", npun);
  dstep.print(out);
  if (dstep.st == sel_no)
    return 0;

  for (i=0; i<npun; i++)
  {
    fprintf(out, "%d ", int(nif[i]));
    switch (nif[i])
    {
      case atnode:
        fprintf (out, "%ld ", pid[i]+1);
        break;
      case atip:
        fprintf (out, "%ld %ld ", eid[i]+1, ipeid[i]+1);
        break;
      case atxyz:
        fprintf (out, "%e %e %e ", x[i], y[i], z[i]);
        break;
      default:
        print_err("unknown type of point is required", __FILE__, __LINE__, __func__); 
        return 1;
    }
    fprintf (out, "%d", int(pu[i]));

    if ((pu[i] == pr_stepid) || (pu[i] == pr_appload) || (pu[i] == pr_time))
    {
      fprintf (out, "\n");
      continue;
    }
    fprintf (out, " %ld\n", ipu[i]+1);
  }

  return 0;
}



/**
  Function prints header for diagram file

  @param out - pointer to opened output text file

  @retval 0 - on success
  @retval 1 - in case of unknown required value type

  Created by Tomas Koudelka,
*/
long outdiagm::print_header(FILE *out)
{
  long i;
  for (i=0; i<npun; i++)
  {
    switch (pu[i])
    {
      case pr_displ:
        fprintf(out, "r_%ld ", ipu[i]);
        break;
      case pr_strains:
        fprintf(out, "eps_%ld ", ipu[i]);
        break;
      case pr_stresses:
        fprintf(out, "sig_%ld ", ipu[i]);
        break;
      case pr_forces:
        fprintf(out, "f_%ld ", ipu[i]);
        break;
      case pr_react:
        fprintf(out, "R_%ld ", ipu[i]);
        break;
      case pr_stepid:
        fprintf(out, "step ");
        break;
      case pr_appload:
        fprintf(out, "lambda ");
        break;
      case pr_time:
        fprintf(out, "time ");
        break;
      case pr_other:
        fprintf(out, "other_%ld ", ipu[i]);
        break;
      default:
        print_err("unknown type of value for diagram is required", __FILE__, __LINE__, __func__, pid[i]+1);
        return 1;
    }
  }
  fprintf(out, "\n");
  return 0;
}



/**
  The function prints one row with required values to the output file out.
  It supposes that values in the integration points and nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param lambda - coefficient of applied load
  @param istep - index of current step
  @param fi - pointer to array load vector


  @retval 0 - on success
  @retval 1 - in case of wrong ip identification
  @retval 2 - in case of unknown required value type

  Created by Tomas Koudelka,
*/
long outdiagm::printval(FILE *out, long lcid, double lambda, long istep, double *fi)
{
  long i;

  // check if we should print at given step
  if (dstep.presence_id(istep, lambda, Mp->timecon) == 0)
    return 0;

  for (i=0; i<npun; i++)
  {
    if ((pid[i] < 0) && (nif[i] == atip))
    {
      if ((eid[i] < Mt->ne) && (ipeid[i] >= 0) && (ipeid[i] < Mt->give_tnip(eid[i])))
        pid[i] = Mt->elements[eid[i]].ipp[0][0]+ipeid[i];
      else
      {
        print_err("invalid element number or element ip number is required\n"
                  "(element %ld, ip %ld)", __FILE__, __LINE__, __func__, eid[i]+1, ipeid[i]+1); 
        return 1;
      }
    }
    switch (pu[i])
    {
      case pr_displ:
        print_displacements(out, lcid, i);
        break;
      case pr_strains:
        print_strains(out, lcid, i);
        break;
      case pr_stresses:
        print_stresses(out, lcid, i);
        break;
      case pr_macrostrain:
        print_macrostrain(out, lcid, i);
        break;
      case pr_macrostress:
        print_macrostress(out, lcid, i);
        break;
      case pr_forces:
        print_forces(out, i, fi);
        break;
      case pr_react:
        print_reactions(out, lcid, i);
        break;
      case pr_stepid :
        fprintf(out, "%ld ", istep);
        break;
      case pr_time:
	if(Mp->tpr == seconds)
	  fprintf(out, "% .15e ", lambda);
	if(Mp->tpr == minutes)
	  fprintf(out, "% .15e ", lambda/60.0);
	if(Mp->tpr == hours)
	  fprintf(out, "% .15e ", lambda/3600.0);
	if(Mp->tpr == days)
	  fprintf(out, "% .15e ", lambda/86400.0);
        break;
      case pr_appload:
      case pr_eigval:
        fprintf(out, "% .15e ", lambda);
        break;
      case pr_other:
        print_others(out, lcid, i);
        break;
      default:
        print_err("unknown type of value for diagram is required at node or ip number %ld",
                  __FILE__, __LINE__, __func__, pid[i]+1);
        return 2;
    }
  }
  fprintf(out, "\n");

  return 0;
}



/**
  The function prints one row with required values to the output file out.
  It supposes that values in the integration points and nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param lambda - coefficient of applied load
  @param istep - index of current step
  @param fi - pointer to array load vector
  @param fr - array of residual %vector components at nodes which should be printed.


  @retval 0 - on success
  @retval 1 - in case of wrong ip identification
  @retval 2 - in case of unknown required value type

  Created by Tomas Koudelka,
*/
long outdiagm::printval(FILE *out, long lcid, double lambda, long istep, double *fi, double *fr)
{
  long i;

  // check if we should print at given step
  if (dstep.presence_id(istep, lambda, Mp->timecon) == 0)
    return 0;

  for (i=0; i<npun; i++)
  {
    if ((pid[i] < 0) && (nif[i] == atip))
    {
      if ((eid[i] < Mt->ne) && (ipeid[i] >= 0) && (ipeid[i] < Mt->give_tnip(eid[i])))
        pid[i] = Mt->elements[eid[i]].ipp[0][0]+ipeid[i];
      else
      {
        print_err("invalid element number or element ip number is required\n"
                  "(element %ld, ip %ld)", __FILE__, __LINE__, __func__, eid[i]+1, ipeid[i]+1); 
        return 1;
      }
    }
    switch (pu[i])
    {
      case pr_displ:
        print_displacements(out, lcid, i);
        break;
      case pr_strains:
        print_strains(out, lcid, i);
        break;
      case pr_stresses:
        print_stresses(out, lcid, i);
        break;
      case pr_macrostrain:
        print_macrostrain(out, lcid, i);
        break;
      case pr_macrostress:
        print_macrostress(out, lcid, i);
        break;
      case pr_forces:
        print_forces(out, i, fi);
        break;
      case pr_residual:
        print_forces(out, i, fr);
        break;
      case pr_react:
        print_reactions(out, lcid, i);
        break;
      case pr_stepid :
        fprintf(out, "%ld ", istep);
        break;
      case pr_time:
	if(Mp->tpr == seconds)
	  fprintf(out, "% .15e ", lambda);
	if(Mp->tpr == minutes)
	  fprintf(out, "% .15e ", lambda/60.0);
	if(Mp->tpr == hours)
	  fprintf(out, "% .15e ", lambda/3600.0);
	if(Mp->tpr == days)
	  fprintf(out, "% .15e ", lambda/86400.0);
        break;
      case pr_appload:
      case pr_eigval:
        fprintf(out, "% .15e ", lambda);
        break;
      case pr_other:
        print_others(out, lcid, i);
        break;
      default:
        print_err("unknown type of value for diagram is required at node or ip number %ld",
                  __FILE__, __LINE__, __func__, pid[i]+1);
        return 2;
    }
  }
  fprintf(out, "\n");

  return 0;
}



/**
  The function prints one row with required values to the output file out without respect to
  to selected step/time (forced printing). It supposes that values in the integration points 
  and nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param lambda - coefficient of applied load
  @param istep - index of current step
  @param fi - pointer to array load vector


  @retval 0 - on success
  @retval 1 - in case of wrong ip identification
  @retval 2 - in case of unknown required value type

  Created by Tomas Koudelka,
*/
long outdiagm::printval_forced(FILE *out, long lcid, double lambda, long istep, double *fi)
{
  long i;


  for (i=0; i<npun; i++)
  {
    if ((pid[i] < 0) && (nif[i] == atip))
    {
      if ((eid[i] < Mt->ne) && (ipeid[i] >= 0) && (ipeid[i] < Mt->give_tnip(eid[i])))
        pid[i] = Mt->elements[eid[i]].ipp[0][0]+ipeid[i];
      else
      {
        print_err("invalid element number or element ip number is required\n"
                  "(element %ld, ip %ld)", __FILE__, __LINE__, __func__, eid[i]+1, ipeid[i]+1); 
        return 1;
      }
    }
    switch (pu[i])
    {
      case pr_displ:
        print_displacements(out, lcid, i);
        break;
      case pr_strains:
        print_strains(out, lcid, i);
        break;
      case pr_stresses:
        print_stresses(out, lcid, i);
        break;
      case pr_forces:
        print_forces(out, i, fi);
        break;
      case pr_react:
        print_reactions(out, lcid, i);
        break;
      case pr_stepid:
        fprintf(out, "%ld ", istep);
        break;
      case pr_time:
	if(Mp->tpr == seconds)
	  fprintf(out, "% .15e ", lambda);
	if(Mp->tpr == minutes)
	  fprintf(out, "% .15e ", lambda/60.0);
	if(Mp->tpr == hours)
	  fprintf(out, "% .15e ", lambda/3600.0);
	if(Mp->tpr == days)
	  fprintf(out, "% .15e ", lambda/86400.0);
        break;
      case pr_appload:
      case pr_eigval:
        fprintf(out, "% .15e ", lambda);
        break;
      case pr_other:
        print_others(out, lcid, i);
        break;
      default:
        print_err("unknown type of value for diagram is required at node or ip number %ld",
                  __FILE__, __LINE__, __func__, pid[i]+1);
        return 2;
    }
  }
  fprintf(out, "\n");

  return 0;
}



/**
  The function prints one row with required values to the output file out without respect to
  to selected step/time (forced printing). It supposes that values in the integration points 
  and nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param lambda - coefficient of applied load
  @param istep - index of current step
  @param fi - pointer to array load vector
  @param fr - array of residual %vector components at nodes which should be printed.


  @retval 0 - on success
  @retval 1 - in case of wrong ip identification
  @retval 2 - in case of unknown required value type

  Created by Tomas Koudelka,
*/
long outdiagm::printval_forced(FILE *out, long lcid, double lambda, long istep, double *fi, double *fr)
{
  long i;


  for (i=0; i<npun; i++)
  {
    if ((pid[i] < 0) && (nif[i] == atip))
    {
      if ((eid[i] < Mt->ne) && (ipeid[i] >= 0) && (ipeid[i] < Mt->give_tnip(eid[i])))
        pid[i] = Mt->elements[eid[i]].ipp[0][0]+ipeid[i];
      else
      {
        print_err("invalid element number or element ip number is required\n"
                  "(element %ld, ip %ld)", __FILE__, __LINE__, __func__, eid[i]+1, ipeid[i]+1); 
        return 1;
      }
    }
    switch (pu[i])
    {
      case pr_displ:
        print_displacements(out, lcid, i);
        break;
      case pr_strains:
        print_strains(out, lcid, i);
        break;
      case pr_stresses:
        print_stresses(out, lcid, i);
        break;
      case pr_forces:
        print_forces(out, i, fi);
        break;
      case pr_residual:
        print_forces(out, i, fr);
        break;
      case pr_react:
        print_reactions(out, lcid, i);
        break;
      case pr_stepid:
        fprintf(out, "%ld ", istep);
        break;
      case pr_time:
	if(Mp->tpr == seconds)
	  fprintf(out, "% .15e ", lambda);
	if(Mp->tpr == minutes)
	  fprintf(out, "% .15e ", lambda/60.0);
	if(Mp->tpr == hours)
	  fprintf(out, "% .15e ", lambda/3600.0);
	if(Mp->tpr == days)
	  fprintf(out, "% .15e ", lambda/86400.0);
        break;
      case pr_appload:
      case pr_eigval:
        fprintf(out, "% .15e ", lambda);
        break;
      case pr_other:
        print_others(out, lcid, i);
        break;
      default:
        print_err("unknown type of value for diagram is required at node or ip number %ld",
                  __FILE__, __LINE__, __func__, pid[i]+1);
        return 2;
    }
  }
  fprintf(out, "\n");

  return 0;
}



/**
  The function prints displacement at required point to the output file out.
  It supposes that values in the nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param idp - point index

  @retval 0 - on success
  @retval 1 - in case that value at integration point is required

  Created by Tomas Koudelka,
*/
long outdiagm::print_displacements(FILE *out, long lcid, long idp)
{
  double *r;
  switch (nif[idp])
  {
    case atxyz:
    case atnode:
      r = new double [Gtm->gnodes[pid[idp]].ndofn];
      noddispl (lcid, r, pid[idp]);
      fprintf(out, "% .15e ", r[ipu[idp]]);
      delete [] r;
      break;
    case atip:
      print_err("unsupported combination of printed values is required -\n"
                " displacements at ip number %ld", __FILE__, __LINE__, __func__, pid[idp]+1);
      return 1;
    default:
      print_err("unknown type of point is required (point number %ld)", __FILE__, __LINE__, __func__, pid[idp]+1); 
      return 1;
  }
  return  0;
}



/**
  The function prints strain at required point to the output file out.
  It supposes that values in the integration points and nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param idp - point index

  @retval 0 - on success
  @retval 1 - in case that unknown type of point is required
  @retval 2 - in case that unknown component of strain is required

  Created by Tomas Koudelka,
*/
long outdiagm::print_strains(FILE *out, long /*lcid*/, long idp)
{
  switch (nif[idp])
  {
    case atxyz:
    case atnode:
      if (ipu[idp] < Mt->nodes[pid[idp]].ncompstr)
        fprintf(out, "% .15e ", Mt->nodes[pid[idp]].strain[ipu[idp]]);
      else
      {
        print_err("invalid component number %ld is required at node number %ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1); 
        return 2;
      }
      break;
    case atip:
      if (ipu[idp] < Mm->ip[pid[idp]].ncompstr)
        fprintf(out, "% .15e ", Mm->ip[pid[idp]].strain[ipu[idp]]);
      else
      {
        print_err("invalid component number %ld is required at integration point number %ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1);
        return 2;
      }
      break;
    default:
      print_err("unknown type of point is required (point number %ld)", __FILE__, __LINE__, __func__, pid[idp]+1); 
      return 1;
  }
  return 0;
}



/**
  The function prints stress at required point to the output file out.
  It supposes that values in the integration points and nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param idp - point index

  @retval 0 - on success
  @retval 1 - in case that unknown type of point is required
  @retval 2 - in case that unknown component of strain is required

  Created by Tomas Koudelka,
*/
long outdiagm::print_stresses(FILE *out, long /*lcid*/, long idp)
{
  switch (nif[idp])
  {
    case atxyz:
    case atnode:
      if (ipu[idp] < Mt->nodes[pid[idp]].ncompstr)
        fprintf(out, "% .15e ", Mt->nodes[pid[idp]].stress[ipu[idp]]);
      else
      {
        print_err("invalid component number %ld is required at node number %ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1); 
        return 2;
      }
      break;
    case atip:
      if (ipu[idp] < Mm->ip[pid[idp]].ncompstr)
        fprintf(out, "% .15e ", Mm->ip[pid[idp]].stress[ipu[idp]]);
      else
      {
        print_err("invalid component number %ld is required at integration point number %ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1);
        return 2;
      }
      break;
    default:
      print_err("unknown type of point is required (point number %ld)", __FILE__, __LINE__, __func__, pid[idp]+1); 
      return 1;
  }
  return 0;
}



/**
  The function prints macro-strain components at required point to the output file out.
  It supposes that values in the integration points and nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param idp - point index

  @retval 0 - on success
  @retval 1 - in case that unknown type of point is required
  @retval 2 - in case that unknown component of strain is required
  @retval 3 - in case that macro-stress/strain components are not defined in the problem

  Created by Tomas Koudelka, 9.1.2015
*/
long outdiagm::print_macrostrain(FILE *out, long lcid, long idp)
{
  if ((Mp->homog == 3) || (Mp->homog == 4) || (Mp->homog == 9))
  {
    long ncomp = Mm->max_ncompstre;
    switch (nif[idp])
    {
      case atxyz:
      case atnode:
        if (ipu[idp] < ncomp)
          fprintf(out, "% .15e ", macrostrains(lcid, ipu[idp]));
        else
        {
          print_err("invalid component number %ld is required at node number %ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1); 
          return 2;
        }
        break;
      case atip:
        if (ipu[idp] < ncomp)
          fprintf(out, "% .15e ", macrostrains(lcid, ipu[idp]));
        else
        {
          print_err("invalid component number %ld is required at integration point number %ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1);
          return 2;
        }
        break;
      default:
        print_err("unknown type of point is required (point number %ld)", __FILE__, __LINE__, __func__, pid[idp]+1); 
        return 1;
    }
  }
  else
  {
    print_err("macro-stress/strain is required but they are not defined\n in the problem (point number %ld)",
               __FILE__, __LINE__, __func__, pid[idp]+1); 
    return 3;
  }

  return 0;
}



/**
  The function prints macro-strain components at required point to the output file out.
  It supposes that values in the integration points and nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param idp - point index

  @retval 0 - on success
  @retval 1 - in case that unknown type of point is required
  @retval 2 - in case that unknown component of strain is required
  @retval 3 - in case that macro-stress/strain components are not defined in the problem

  Created by Tomas Koudelka, 9.1.2015
*/
long outdiagm::print_macrostress(FILE *out, long lcid, long idp)
{
  if ((Mp->homog == 3) || (Mp->homog == 4) || (Mp->homog == 9))
  {
    long ncomp = Mm->max_ncompstre;    

    switch (nif[idp])
    {
      case atxyz:
      case atnode:
        if (ipu[idp] < ncomp)
          fprintf(out, "% .15e ", macrostresses(lcid, ipu[idp]));
        else
        {
          print_err("invalid component number %ld is required at node number %ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1); 
          return 2;
        }
        break;
      case atip:
        if (ipu[idp] < ncomp)
          fprintf(out, "% .15e ", macrostresses(lcid, ipu[idp]));
        else
        {
          print_err("invalid component number %ld is required at integration point number %ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1);
          return 2;
        }
        break;
      default:
        print_err("unknown type of point is required (point number %ld)", __FILE__, __LINE__, __func__, pid[idp]+1); 
        return 1;
    }
  }
  else
  {
    print_err("macro-stress/strain is required but they are not defined\n in the problem (point number %ld)",
               __FILE__, __LINE__, __func__, pid[idp]+1); 
    return 3;
  }

  return 0;
}



/**
  The function prints forces at required point to the output file out.
  Value from load vector fi is used and printed.

  @param out - pointer to opened output text file
  @param idp - point index

  @retval 0 - on success
  @retval 1 - in case that point type is atip
  @retval 2 - in case that unknown component of force is required

  Created by Tomas Koudelka,
*/
long outdiagm::print_forces(FILE *out, long idp, double *fi)
{
  long ii;
  switch (nif[idp])
  {
    case atxyz:
    case atnode:
      ii=Mt->give_dof(pid[idp], ipu[idp]);
      if (ii > 0)
        fprintf(out, "% .15e ", fi[ii-1]);
      else      
      {
        print_err("invalid component number is required at node number %ld", __FILE__, __LINE__, __func__, pid[idp]+1); 
        return 2;
      }
      break;
    case atip:
      print_err("unsupported combination of printed values is required -\n"
                " forces in ip number %ld",
                __FILE__, __LINE__, __func__, pid[idp]+1);
      return 1;
    default:
      print_err("unknown type of point is required (point number %ld)", __FILE__, __LINE__, __func__, pid[idp]+1); 
      return 1;
  }
  return 0;
}



/**
  The function prints reaction at required node to the output file out.
  It supposes that values in the nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param idp - required component of reaction

  @retval 0 - on success
  @retval 1 - in case that unknown type of point is required

  Created by Tomas Koudelka,
*/
long outdiagm::print_reactions(FILE *out, long /*lcid*/, long idp)
{
  switch (nif[idp])
  {
    case atxyz:
    case atnode:
      if (Mt->nodes[pid[idp]].react)
        fprintf(out, "% .15e ", Mt->nodes[pid[idp]].r[ipu[idp]]);
      break;
    case atip:
      print_err("unsupported combination of printed values is required -\n" 
                " forces in ip number %ld", __FILE__, __LINE__, __func__, pid[idp]+1);
      return 1;
    default:
      print_err("unknown type of point is required (point number %ld)", __FILE__, __LINE__, __func__, pid[idp]+1); 
      return 1;
  }
  return 0;
}



/**
  The function prints other value at required point to the output file out.
  It supposes that values in the integration points and nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param idp - required component of other array

  @retval 0 - on success
  @retval 1 - in case that unknown type of point is required
  @retval 2 - in case that unknown component of strain is required

  Created by Tomas Koudelka,
*/
long outdiagm::print_others(FILE *out, long /*lcid*/, long idp)
{
  switch (nif[idp])
  {
    case atxyz:
    case atnode:
      if (ipu[idp] < Mt->nodes[pid[idp]].ncompother)
        fprintf(out, "% .15e ", Mt->nodes[pid[idp]].other[ipu[idp]]);
      else
      {
        print_err("invalid component number %ld is required at node number %ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1); 
        return 2;
      }
      break;
    case atip:
      if (ipu[idp] < Mm->ip[pid[idp]].ncompeqother)
        fprintf(out, "% .15e ", Mm->ip[pid[idp]].other[ipu[idp]]);
      else
      {
        print_err("invalid component number %ld is required at integration point number %ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1); 
        return 2;
      }
      break;
    default:
      print_err("unknown type of point is required (point number %ld)", __FILE__, __LINE__, __func__, pid[idp]+1); 
      return 1;
  }
  return 0;
}
