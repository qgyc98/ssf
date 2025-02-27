#include <string.h>
#include "outdiagt.h"
#include "gtopology.h"
#include "globalt.h"
#include "probdesct.h"
//#include "meshtransfert.h"
#include "transtop.h"
#include "elementt.h"
#include "nodet.h"
#include "intpointst.h"
#include "globmatt.h"
#include "selection.h"



/**
  The constructor initializes data to zero values

  Created by TKr, TKo
*/
outdiagt::outdiagt()
{
  pid = NULL; eid = NULL; ipeid = NULL; nif = NULL;
  pu = NULL;
  ipu = NULL;
  x = y = z = NULL;
}


/**
  The destructor deallocates used memory

  Created by TKr, TKo
*/
outdiagt::~outdiagt()
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

  Created by TKr, TKo
*/
long outdiagt::read(XFILE *in)
{
  long i;
  if (xfscanf (in, "%k%ld", "numunknowns", &npun) != 2)
  {
    print_err("cannot read number of printed unknowns", __FILE__,__LINE__,__func__);
    return 1;
  }
  xfscanf (in, "%k", "sel_diagstep");
  dstep.read(in);
  if (dstep.st == sel_no)
    return 0;
  nif   = new nodip[npun];
  pid   = new long [npun];
  eid   = new long [npun];
  ipeid = new long [npun];
  pu    = new prunkt[npun];
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
        print_err("cannot read type of point", __FILE__,__LINE__,__func__);
        return 1;
      }
      switch (nif[i])
	{
	case atnode:
	  if (xfscanf (in, "%k%ld", "node", pid+i) != 2)
	  {
	    print_err("cannot read node id", __FILE__,__LINE__,__func__);
            return 1;
          }
	  pid[i]--;
	  break;
	case atip:
	  if (xfscanf (in, "%k %ld %k %ld", "elem", eid+i, "ip", ipeid+i) != 4)
	  {
            print_err("cannot read element or ip number", __FILE__,__LINE__,__func__);
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
          if (Gtt->gnodes) // due to preprocessor where pid[i] need not to be set
            pid[i] = Gtt->give_nearest_node(x[i],y[i],z[i]);
          break;
	default:
	  print_err("unknown type of point is required", __FILE__, __LINE__, __func__); 
	  return 1;
	}
      
      if (xfscanf (in, "%k%m", "unknowntype", &prunkt_kwdset, (int *)(pu+i)) != 2)
    	{
	  print_err("cannot read type of printed unknown", __FILE__,__LINE__, __func__);
	  return 1;
	}
      if ((pu[i] == pr_stepidt) || (pu[i] == pr_apploadt) || (pu[i] == pr_timet))
	continue;
      if (xfscanf (in, "%k%ld", "compid", ipu+i) != 2)
	{
	  print_err("cannot read component of printed unknown", __FILE__,__LINE__, __func__);
	  return 1;
	}
      ipu[i]--;
      if (nif[i] == atnode)
	{
	  switch (pu[i])
	    {
	    case pr_gradients :
	      if(Tp->gradcomp == 0)
		Tp->gradcomp = 1;
	      if(Tp->gradaver == 0)
		Tp->gradaver = 1;
	      Tp->gradpos = 2;
	      break;
	    case pr_fluxes :
	      if(Tp->gradcomp == 0)
		Tp->gradcomp = 1;
	      if(Tp->fluxcomp == 0)
		Tp->fluxcomp = 1;
	      if(Tp->gradaver == 0)
		Tp->gradaver = 1;
	      Tp->gradpos = 2;
	      if(Tp->fluxaver == 0)
		Tp->fluxaver = 1;
	      Tp->fluxpos = 2;
	      break;
	    case pr_othert :
	      if(Tp->othercomp == 0)
		Tp->othercomp = 1;
	      if(Tp->otheraver == 0)
		Tp->otheraver = 1;
	      if(Tp->otherpos == 0)
		Tp->otherpos = 3;
	      break;
	    case pr_eqothert :
	      if(Tp->eqothercomp == 0)
		Tp->eqothercomp = 1;
	      if(Tp->eqotheraver == 0)
		Tp->eqotheraver = 1;
	      Tp->eqotherpos = 2;
	      break;
	    default:
	      break;
	    }
	}
      if ((nif[i] == atip) || (nif[i] == atxyz))
	{
	  switch (pu[i])
	    {
	    case pr_gradients :
	      if(Tp->gradcomp == 0)
		Tp->gradcomp = 1;
	      if(Tp->gradaver == 0)
		Tp->gradaver = 1;
	      if(Tp->gradpos == 0)
		Tp->gradpos = 1;
	      break;
	    case pr_fluxes :
	      if(Tp->gradcomp == 0)
		Tp->gradcomp = 1;
	      if(Tp->fluxcomp == 0)
		Tp->fluxcomp = 1;
	      if(Tp->gradaver == 0)
		Tp->gradaver = 1;
	      if(Tp->gradpos == 0)
		Tp->gradpos = 1;
	      if(Tp->fluxaver == 0)
		Tp->fluxaver = 1;
	      if(Tp->fluxpos == 0)
		Tp->fluxpos = 1;
	      break;
	    case pr_othert :
	      if(Tp->othercomp == 0)
		Tp->othercomp = 1;
	      if(Tp->otheraver == 0)
		Tp->otheraver = 1;
	      if(Tp->otherpos == 0)
		Tp->otherpos = 1;
	      break;
	    case pr_eqothert :
	      if(Tp->eqothercomp == 0)
		Tp->eqothercomp = 1;
	      if(Tp->eqotheraver == 0)
		Tp->eqotheraver = 1;
	      if(Tp->eqotherpos == 0)
		Tp->eqotherpos = 1;
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

  Created by TKr, TKo
*/
long outdiagt::print(FILE *out)
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

    if ((pu[i] == pr_stepidt) || (pu[i] == pr_apploadt) || (pu[i] == pr_timet))
    {
      fprintf (out, "\n");
      continue;
    }
    fprintf (out, " %ld\n", ipu[i]+1);
  }

  return 0;
}



/**
  Function prints header for diagram file.

  @param out - pointer to opened output text file

  @retval 0 - on success
  @retval 1 - in case of unknown required value type

  Created by TKr, TKo, this should be controled again 
*/
long outdiagt::print_header(FILE *out)
{
  long i;
  char emsg[200];

  fprintf(out, "# ");
  for (i=0; i<npun; i++)
    {
      switch (pu[i])
	{
	case pr_unknowns :
	  fprintf(out, "unknown_%ld ", ipu[i]);
	  break;
	case pr_gradients :
	  fprintf(out, "gradient_%ld ", ipu[i]);
	  break;
	case pr_fluxes :
	  fprintf(out, "flux_%ld ", ipu[i]);
	  break;
	case pr_surffluxes :
	  fprintf(out, "surfflux_%ld ", ipu[i]);
	  break;
	case pr_stepidt :
	  fprintf(out, "step ");
	  break;
	case pr_apploadt :
	  fprintf(out, "lambda ");
	  break;
	case pr_timet :
	  fprintf(out, "time ");
	  break;
	case pr_othert :
	  fprintf(out, "other_%ld ", ipu[i]);
	  break;
	case pr_eqothert :
	  fprintf(out, "eqother_%ld ", ipu[i]);
	  break;
	default:
	  print_err("unknown type of value for diagram is required in node or ip number %ld", __FILE__, __LINE__, __func__, pid[i]+1);
	  return 1;
	}
    }
  fprintf(out, "\n");
  fprintf(out, "#");
  for (i=0; i<npun; i++)
  {
    if ((pid[i] < 0) && (nif[i] == atip))
    {
      if ((eid[i] < Tt->ne) && (ipeid[i] >= 0) && (ipeid[i] < Tt->give_tnip(eid[i])))
        pid[i] = Tt->elements[eid[i]].ipp[0][0]+ipeid[i];
      else
      {
        sprintf(emsg, "invalid element number (%ld) or element ip number (%ld) is required", eid[i]+1, ipeid[i]+1);
        print_err(emsg, __FILE__, __LINE__, __func__);
        return 1;
      }
    }
    switch(nif[i])
    {
      case atnode:
      case atxyz:
        fprintf(out, " node=%ld", pid[i]+1);
        break;
      case atip:
        fprintf(out, " elem=%ld", eid[i]+1);
        break;
      default:
        print_err("unknown type of value for diagram is required in node or ip number %ld", __FILE__, __LINE__, __func__, pid[i]+1);
        return 1;
    }
  }
  fprintf(out, "\n");        

  fprintf(out, "#");
  for (i=0; i<npun; i++)
  {
    if ((pid[i] < 0) && (nif[i] == atip))
    {
      if ((eid[i] < Tt->ne) && (ipeid[i] >= 0) && (ipeid[i] < Tt->give_tnip(eid[i])))
        pid[i] = Tt->elements[eid[i]].ipp[0][0]+ipeid[i];
      else
      {
        sprintf(emsg, "invalid element number (%ld) or element ip number (%ld) is required", eid[i]+1, ipeid[i]+1);
        print_err(emsg, __FILE__, __LINE__, __func__);
        return 1;
      }
    }
    switch(nif[i])
    {
      case atnode:
      case atxyz:
        fprintf(out, " %le", Gtt->gnodes[pid[i]].x);
        break;
      case atip:
	//        fprintf(out, " %le", Gtt->gnodes[pid[i]].x);
        break;
      // case atxyz:
      //   fprintf(out, "  ---");
      //   break;
      default:
        print_err("unknown point type is required", __FILE__, __LINE__, __func__);
        return 1;
    }
  }
  fprintf(out, "\n");

  fprintf(out, "# ");
  for (i=0; i<npun; i++)
  {
    if ((pid[i] < 0) && (nif[i] == atip))
    {
      if ((eid[i] < Tt->ne) && (ipeid[i] >= 0) && (ipeid[i] < Tt->give_tnip(eid[i])))
        pid[i] = Tt->elements[eid[i]].ipp[0][0]+ipeid[i];
      else
      {
        sprintf(emsg, "invalid element number (%ld) or element ip number (%ld) is required", eid[i]+1, ipeid[i]+1);
        print_err(emsg, __FILE__, __LINE__, __func__);
        return 1;
      }
    }
    switch(nif[i])
    {
      case atnode:
      case atxyz:
        fprintf(out, " %le", Gtt->gnodes[pid[i]].y);
        break;
      case atip:
	//        fprintf(out, " %le", Gtt->gnodes[pid[i]].y);
        break;
      // case atxyz:
      //   fprintf(out, "  ---");
      //   break;
      default:
        print_err("unknown point type is required", __FILE__, __LINE__, __func__);
        return 1;
    }
  }
  fprintf(out, "\n");

  fprintf(out, "# ");
  for (i=0; i<npun; i++)
  {
    if ((pid[i] < 0) && (nif[i] == atip))
    {
      if ((eid[i] < Tt->ne) && (ipeid[i] >= 0) && (ipeid[i] < Tt->give_tnip(eid[i])))
        pid[i] = Tt->elements[eid[i]].ipp[0][0]+ipeid[i];
      else
      {
        sprintf(emsg, "invalid element number (%ld) or element ip number (%ld) is required", eid[i]+1, ipeid[i]+1);
        print_err(emsg, __FILE__, __LINE__, __func__);
        return 1;
      }
    }
    switch(nif[i])
    {
      case atnode:
      case atxyz:
        fprintf(out, " %le", Gtt->gnodes[pid[i]].z);
        break;
      case atip:
	//        fprintf(out, " %le", Gtt->gnodes[pid[i]].z);
        break;
      // case atxyz:
      //   fprintf(out, "  ---");
      //   break;
      default:
        print_err("unknown point type is required", __FILE__, __LINE__, __func__);
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

  Created by TKr, TKo
*/
long outdiagt::printval(FILE *out, long lcid, double lambda, long istep, double */*fi*/)
{
  long i;

  // check if we should print at given step
  if (dstep.presence_id(istep, lambda, Tp->timecont) == 0)
    return 0;

  for (i=0; i<npun; i++)
  {
    if ((pid[i] < 0) && (nif[i] == atip))
    {
      if ((eid[i] < Tt->ne) && (ipeid[i] >= 0) && (ipeid[i] < Tt->give_tnip(eid[i])))
        pid[i] = Tt->elements[eid[i]].ipp[0][0]+ipeid[i];
      else
      {
        print_err("invalid element number or element ip number is required on element %ld, ip %ld", __FILE__, __LINE__, __func__, eid[i]+1, ipeid[i]+1); 
        return 1;
      }
    }
    switch (pu[i])
      {
      case pr_unknowns :
        print_unknowns(out,lcid,i);
        break;
      case pr_gradients:
	print_gradients(out,lcid,i);
       	break;
      case pr_fluxes:
	print_fluxes(out,lcid,i);
	break;
      case pr_surffluxes:
	print_surffluxes(out,lcid,i);
	break;
      case pr_othert:
	print_others(out,lcid,i);
	break;
      case pr_eqothert:
	print_eqothers(out,lcid,i);
	break;
      case pr_stepidt :
        fprintf(out, "%ld ", istep);
        break;
      case pr_apploadt :
      case pr_timet:
	switch (Tp->tprt)
	  {
	  case secondst:
	    fprintf(out, " %.15e ", lambda);
	    break;
	  case minutest:
	    fprintf(out, " %.15e ", lambda/60.0);
	    break;
	  case hourst:
	    fprintf(out, " %.15e ", lambda/3600.0);
	    break;
	  case dayst:
	    fprintf(out, " %.15e ", lambda/86400.0);
	    break;
	  default:
	    {
	      print_err("unknown type of time printing",__FILE__,__LINE__,__func__);
	    }	
	  }
	
        break;
      default:
        print_err("unknown type of value for diagram is required in node or ip number %ld", __FILE__, __LINE__, __func__, pid[i]+1);
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


  @retval 0 - on success
  @retval 1 - in case of wrong ip identification
  @retval 2 - in case of unknown required value type

  Created by TKo
*/
long outdiagt::printval_forced(FILE *out, long lcid, double lambda, long istep, double */*fi*/)
{
  long i;

  for (i=0; i<npun; i++)
  {
    if ((pid[i] < 0) && (nif[i] == atip))
    {
      if ((eid[i] < Tt->ne) && (ipeid[i] >= 0) && (ipeid[i] < Tt->give_tnip(eid[i])))
        pid[i] = Tt->elements[eid[i]].ipp[0][0]+ipeid[i];
      else
      {
        print_err("invalid element number or element ip number is required on element %ld, ip %ld", __FILE__, __LINE__, __func__, eid[i]+1, ipeid[i]+1); 
        return 1;
      }
    }
    switch (pu[i])
      {
      case pr_unknowns :
        print_unknowns(out,lcid,i);
        break;
      case pr_gradients:
	print_gradients(out,lcid,i);
       	break;
      case pr_fluxes:
	print_fluxes(out,lcid,i);
	break;
      case pr_surffluxes:
	print_surffluxes(out,lcid,i);
	break;
      case pr_othert:
	print_others(out,lcid,i);
	break;
      case pr_eqothert:
	print_eqothers(out,lcid,i);
	break;
      case pr_stepidt :
        fprintf(out, "%ld ", istep);
        break;
      case pr_apploadt :
      case pr_timet:
	switch (Tp->tprt)
	  {
	  case secondst:
	    fprintf(out, " %.15e ", lambda);
	    break;
	  case minutest:
	    fprintf(out, " %.15e ", lambda/60.0);
	    break;
	  case hourst:
	    fprintf(out, " %.15e ", lambda/3600.0);
	    break;
	  case dayst:
	    fprintf(out, " %.15e ", lambda/86400.0);
	    break;
	  default:
	    {
	      print_err("unknown type of time printing",__FILE__,__LINE__,__func__);
	    }	
	  }
	
        break;
      default:
        print_err("unknown type of value for diagram is required in node or ip number %ld", __FILE__, __LINE__, __func__, pid[i]+1);
        return 2;
    }
  }
  fprintf(out, "\n");

  return 0;
}



/**

  The function prints unknowns at required point to the output file out.
  It supposes that values in the nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param idp - index of printed unknown

  @retval 0 - on success
  @retval 1 - in case that value at integration point is required
  
  Created by TKr, TKo
*/
long outdiagt::print_unknowns(FILE *out, long /*lcid*/, long idp)
{
  double r;
  switch (nif[idp])
    {
    case atnode:
    case atxyz:
      r = nodalval(pid[idp], ipu[idp]);
      fprintf(out, "% e ", r);
      break;
    case atip:      
      r = Tm->ip[pid[idp]].av[ipu[idp]];
      fprintf(out, "% e ", r);
      break;
    default:
      print_err("unknown type of point is required in point number %ld", __FILE__, __LINE__, __func__, pid[idp]+1); 
      return 1;
  }
  return  0;
}

/**

  This function prints gradient value at required point to the output file out.
  It supposes that values in the integration points and nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param idp - required component of gradient array

  @retval 0 - on success
  @retval 1 - in case that unknown type of point is required
  @retval 2 - in case that unknown component of strain is required
*/
long outdiagt::print_gradients(FILE *out, long /*lcid*/, long idp)
{
  long i;

  switch (nif[idp])
    {
    case atnode:
    case atxyz:
      if (ipu[idp] < Tp->ntm){
	for(i=0;i<Tt->nodes[pid[idp]].ncompgrad;i++)
	  fprintf(out, "% e ", Tt->nodes[pid[idp]].gradient[ipu[idp]][i]);
      }
      else
      {
        print_err("invalid component number %ld is required at node %ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1); 
        return 2;
      }
      break;
    case atip:
      if (ipu[idp] < Tp->ntm){
	for(i=0;i<Tm->ip[pid[idp]].ncompgrad;i++)
	  fprintf(out, "% e ", Tm->ip[pid[idp]].grad[ipu[idp]][i]);
      }
      else
      {
        print_err("invalid component number %ld is required at ip=%ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1); 
        return 2;
      }
      break;
    default:
      print_err("unknown type of point is required in point number %ld", __FILE__, __LINE__, __func__, pid[idp]+1); 
      return 1;
  }
  return 0;
}


/**
  This function prints flux value at required point to the output file out.
  It supposes that values in the integration points and nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param idp - required component of flux array

  @retval 0 - on success
  @retval 1 - in case that unknown type of point is required
  @retval 2 - in case that unknown component of strain is required
*/
long outdiagt::print_fluxes(FILE *out, long /*lcid*/, long idp)
{
  long i;

  switch (nif[idp])
    {
    case atnode:
    case atxyz:
      if (ipu[idp] < Tp->ntm){
	for(i=0;i<Tt->nodes[pid[idp]].ncompgrad;i++)
	  fprintf(out, "% e ", Tt->nodes[pid[idp]].flux[ipu[idp]][i]);
      }
      else
      {
        print_err("invalid component number %ld is required at node %ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1); 
        return 2;
      }
      break;
    case atip:
      if (ipu[idp] < Tp->ntm){
	for(i=0;i<Tm->ip[pid[idp]].ncompgrad;i++)
	  fprintf(out, "% e ", Tm->ip[pid[idp]].fluxes[ipu[idp]][i]);
      }
      else
      {
        print_err("invalid component number %ld is required at ip=%ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1); 
        return 2;
      }
      break;
    default:
      print_err("unknown type of point is required in point number %ld", __FILE__, __LINE__, __func__, pid[idp]+1); 
      return 1;
  }
  return 0;
}




/**
  The function prints computed surface fluxes to the output file out.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param idp - point index

  @retval 0 - on success
  @retval 1 - in case that unknown type of point is required
  @retval 2 - in case that unknown component of surface flux is required

  Created by Tomas Krejci according to Tomas Koudelka, 26/2/2021
*/
long outdiagt::print_surffluxes(FILE *out, long lcid, long idp)
{
  switch (nif[idp])
    {
    case atxyz:
    case atnode:
      if (ipu[idp] <  Tb->nbf)
	fprintf(out, "% .15e ", surface_fluxes(lcid, ipu[idp]));
      else
        {
          print_err("invalid component number is required at node number %ld", __FILE__, __LINE__, __func__, pid[idp]+1); 
          return 2;
        }
      break;
    case atip:
      if (ipu[idp] < Tb->nbf)
	fprintf(out, "% .15e ", surface_fluxes(lcid, ipu[idp]));
      else
        {
          print_err("invalid component number is required at integration point number %ld", __FILE__, __LINE__, __func__, pid[idp]+1);
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
  The function prints other value at required point to the output file out.
  It supposes that values in the integration points and nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param idp - required component of other array

  @retval 0 - on success
  @retval 1 - in case that unknown type of point is required
  @retval 2 - in case that unknown component of strain is required

  Created by TKr, TKo
*/
long outdiagt::print_others(FILE *out, long lcid, long idp)
{
  long ncompother = 0;
  switch (nif[idp])
    {
    case atnode:
    case atxyz:
      if (ipu[idp] < Tt->nodes[pid[idp]].ncompother)
	fprintf(out, "% e ", Tt->nodes[pid[idp]].other[ipu[idp]]);
      else
      {
        print_err("invalid component number %ld is required at node %ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1); 
        return 2;
      }
      break;
    case atip:
      ncompother = Tm->ip[pid[idp]].ncompother;
      if (ipu[idp] < ncompother)
	fprintf(out, "% e ", Tm->ip[pid[idp]].other[ncompother*lcid+ipu[idp]]);
      else
      {
        print_err("invalid component number %ld is required at ip=%ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1); 
        return 2;
      }
      break;
    default:
      print_err("unknown type of point is required in point number %ld", __FILE__, __LINE__, __func__, pid[idp]+1); 
      return 1;
  }
  return 0;
}


/**
  This function prints eqother value at required point to the output file out.
  It supposes that values in the integration points and nodes are actualized and valid.

  @param out - pointer to opened output text file
  @param lcid - load case id
  @param idp - required component of eqother array

  @retval 0 - on success
  @retval 1 - in case that unknown type of point is required
  @retval 2 - in case that unknown component of strain is required
*/
long outdiagt::print_eqothers(FILE *out, long lcid, long idp)
{
  long ncompeqother = 0;
  switch (nif[idp])
    {
    case atnode:
    case atxyz:
      if (ipu[idp] < Tt->nodes[pid[idp]].ncompeqother)
	fprintf(out, "% e ", Tt->nodes[pid[idp]].eqother[ipu[idp]]);
      else
      {
        print_err("invalid component number %ld is required at node %ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1); 
        return 2;
      }
      break;
    case atip:
      ncompeqother = Tm->ip[pid[idp]].ncompeqother;
      if (ipu[idp] < ncompeqother)
	fprintf(out, "% e ", Tm->ip[pid[idp]].eqother[ncompeqother*lcid+ipu[idp]]);
      else
      {
        print_err("invalid component number %ld is required at ip=%ld", __FILE__, __LINE__, __func__, ipu[idp]+1, pid[idp]+1); 
        return 2;
      }
      break;
    default:
      print_err("unknown type of point is required in point number %ld", __FILE__, __LINE__, __func__, pid[idp]+1); 
      return 1;
  }
  return 0;
}
