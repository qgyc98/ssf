#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "slipsurf.h"
#include "global.h"
#include "gtopology.h"
#include "mechmat.h"
#include "mechtop.h"
#include "element.h"




/**
  The function detects evolved plastic zones in the integration points. For each detected plastic zone, 
  new id is generated starting from 1. If the consistency parameter (gamma) is zero in the given integration point 
  then the plastic zone id = 0 will be assigned in such the case. Resulting generated id are stored in the 
  array plast_ip.

  @param plast_ip  - allocated array where the id of the detected plastic araeas will be stored for each 
                     integration point, i.e. plast_ip[ipp] = plastic zone id at the integration point ipp.
  @param gamma_min - the minimum value of the consistency parameter which will be taken into account

  @return The function returns the number of plastic zones detected.

  Created by Tomas Koudelka, 03.2013
*/
long detect_plastic_zones(long *plast_ip, double gamma_min)
{
  long i, j, k;
  long max_gamma_ip;  // integration point number where the maximum of consistency parameter was found
  long max_gamma_eid; // element number where the maximum of consistency parameter was found
  double max_gamma = 0.0; // maximum value of consistency parameter
  long *adjelem;      // list of adjacent elements - adjelem[i] = i-th adjacent element to the given element
  long nadjelem;      // number of adjacent elements to the given element
  long  new_nadjelem; // new actual number of adjacent elements
  long *gpl_adjelem;  // array of flags for all elements: 
                      //  gpl_adjelem[i] = -1 => i-th element has been already searched for plastic state
                      //  gpl_adjelem[i] =  0 => i-th element has been never searched for plastic state
                      //  gpl_adjelem[i] =  1 => i-th element will be searched for plastic state
                                                          
  long *plast_adjelem;// array of plastic flags for adjacent elements - 
                      // plast_adjelem[i] = 1 => i-th element has at least one int. point in plastic state
  long nplast_ip;     // number of integration points in plastic state
  long  zone_id = 0;  // number of 


  for (i=0; i<Mm->tnip; i++)
    plast_ip[i] = -1;


  
  Gtm->adjacelem(NULL);
  gpl_adjelem = new long[Mt->ne];
  do
  {
    // find new integration points with plastic strains
    max_gamma = detect_max_gamma(plast_ip, max_gamma_ip, max_gamma_eid, gamma_min);
  
    if (max_gamma_ip == -1)// no plastic zones were found
      break;

    zone_id++;
    nadjelem = Gtm->nadjelel[max_gamma_eid]; // number of adjacent elements at the ip with actual maximum cons. param
    adjelem = new long [nadjelem];           // array of adjacent element numbers
    memcpy(adjelem, Gtm->adjelel[max_gamma_eid], sizeof(*Gtm->adjelel[max_gamma_eid])*nadjelem);
    plast_adjelem = new long [nadjelem];
    memset(plast_adjelem, 0, sizeof(*plast_adjelem)*nadjelem);
    memset(gpl_adjelem, 0, sizeof(*gpl_adjelem)*Mt->ne);
    do 
    {
      nplast_ip = detect_plast_ip(nadjelem, adjelem, plast_ip, zone_id, plast_adjelem);
      if (nplast_ip == 0)
      {
        nadjelem = 0;
        delete [] adjelem;
        delete [] plast_adjelem;
      }
      else
      {
        for(i=0; i<nadjelem; i++)
          gpl_adjelem[adjelem[i]]=-1;

        new_nadjelem = 0;
        for(i=0; i<nadjelem; i++)
        {
          if (plast_adjelem[i])
          {
            for (j=0; j<Gtm->nadjelel[adjelem[i]]; j++)
            {
              k = Gtm->adjelel[adjelem[i]][j]; // id of j-th adjacent element to adjelem[i]
              if (gpl_adjelem[k] == 0)
              {
                new_nadjelem++;
                gpl_adjelem[k] = 1;
              }
            }
          }
        }
        if (new_nadjelem > nadjelem)
        { 
          delete [] plast_adjelem;
          plast_adjelem = new long[new_nadjelem];
          delete [] adjelem;
          adjelem = new long[new_nadjelem];
        }
        j=0;
        for (i=0; i<Mt->ne; i++)
        {
          if (gpl_adjelem[i] == 1)
          {
            adjelem[j] = i;
            j++;
          }
        }
        nadjelem = new_nadjelem;
        memset(plast_adjelem, 0, sizeof(*plast_adjelem)*nadjelem);        
      }

    }while (nadjelem);

  }while (1);

  delete [] gpl_adjelem;
  return zone_id;
}



/**
  The function searches elements whose numbers are in the array adjelem for integration points
  with evolved plasticity. The integration points which have nonzero consistency parameter (gamma) 
  are marked in the array plast_ip by the number stored in the zone_id. If the integration point of the
  i-th adjacent element was newly marked than the plast_adjelem[i] is set to 1.  

  @param nadjelem - number of elements which will be  searched for int. points in plastic state
  @param adjelem  - array with adjacent element numbers (dimension is nadjelem)
  @param plast_ip - array with plastic zone numbers of int. points,
                    plast_ip[ipp] > 0  - plastic zone number for the given ipp
                    plast_ip[ipp] = 0  - plasticity was detected but plastic zone id has not been assigned yet
                    plast_ip[ipp] = -1 - no plasticity was detected in the given ipp
  @param zone_id  - actual plastic zone id. It will be assigned to the searched int. points that are in the plastic state
  @param plast_adjelem - array of element numbers whose int. points were found being in the plastic state

  @return The function returns number of int. points with nonzero consistency parameter 

  Created by Tomas Koudelka 03.2013
*/
long detect_plast_ip(long nadjelem, long *adjelem, long *plast_ip, long zone_id, long *plast_adjelem)
{
  long i, j, nplast_ip, tnip, ipp;

  nplast_ip = 0;
  for(i=0; i<nadjelem; i++)
  {
    if (Gtm->leso[adjelem[i]]==0)   //  only active elements are searched
      continue;

    tnip = Mt->give_tnip(adjelem[i]);
    ipp  = Mt->elements[adjelem[i]].ipp[0][0];
    for(j=0; j<tnip; j++)
    {
      if (plast_ip[ipp] == 0) // given ip is in plastic state but no zone id was assigned
      { 
        plast_ip[ipp] = zone_id;
        plast_adjelem[i] = 1;
        nplast_ip++;
      }
      ipp++;
    }
  }
  return nplast_ip;
}



/**
  The function searches integration points for the maximum value of the consistency parameter.
  Only those points are searched that the plastic zone identifier have not yet been assigned to
  (i.e. plast_ip[i] is zero or negative). The function returns the maximum value of consistency
  parameter gamma and there are also passed resulting max_gamma_ip and max_gamma_eid by the 
  function arguments.

  @param plast_ip - array with plastic zone numbers of int. points,
                    plast_ip[ipp] > 0  - plastic zone number for the given ipp
                    plast_ip[ipp] = 0  - plasticity was detected but plastic zone id has not been assigned yet
                    plast_ip[ipp] = -1 - no plasticity was detected in the given ipp
  @param max_gamma_ip  - integration point id where the maximum value of the consistency parameter was found (output)
  @param max_gamma_eid - element id where the max_gamma_ip was found (output) 
  @param gamma_min - the minimum value of the consistency parameter which will be taken into account

  @return The function returns the maximum value of the consistency parameter (gamma) which was found in the
          remaining integration points. 

  Created by Tomas Koudelka 03.2013
*/
double detect_max_gamma(long *plast_ip, long &max_gamma_ip, long &max_gamma_eid, double gamma_min)
{
  long ipp = 0;
  long nplast_ip = 0;
  long i, j;
  double gamma;
  double max_gamma = gamma_min; // maximum value of consistency parameter
  max_gamma_ip  = -1;
  max_gamma_eid = -1;

  for(i=0; i<Mt->ne; i++)
  {
    if (Gtm->leso[i]==1)
    {
      //  only active elements are searched
      ipp =  Mt->elements[i].ipp[0][0];
      for (j=0; j<Mt->give_tnip(i); j++)
      {        
        if (plast_ip[ipp] < 1)
        {
          gamma = Mm->give_consparam(ipp);
          if (gamma > max_gamma)
          {
            max_gamma     = gamma;
            max_gamma_ip  = ipp;
            max_gamma_eid = i;
          }
          if (gamma > gamma_min )
          {
            plast_ip[ipp] = 0;
            nplast_ip++;
          }
        }
        ipp++;
      }
    }
  }
  
  return max_gamma;
}


/*
  The function interpolates the plastic zone by the slip surface with help of 
  weight least square method. The exponential approximation of the slip surface 
  is assumed in the form y = a*exp(b*x).

  @param id[in] - plastic zone id
  @param plast_ip[in] - array of plastic zone id for particular int. points (plast[ipp] = id of plastic zone)
  @param a[out] - parameter a of the approximation function y = a*exp(b*x)
  @param b[out] - parameter b of the approximation function y = a*exp(b*x)
  @param minx[out] - minumum x coordinate of the points in the given plastic zone
  @param maxx[out] - maximum x coordinate of the points in the given plastic zone

  @retval 0 - approximation was calculated sucessfully
  @retval 1 - plastic zone contains less than 2 points => no slip surface can be determined

  Created by Tomas Koudelka 03.2013
*/
long approx_slip_surf(long id, long *plast_ip, double &a, double &b, double &minx, double &maxx)
{
  long i,j, ipp, npip, tnip;
  double sw2_x_lny = 0.0;  // \sum (w_i)^2 (x_i)^2 \ln(y_i)
  double sw2_x     = 0.0;  // \sum (w_i)^2 (x_i) 
  double sw2_lny   = 0.0;  // \sum (w_i)^2 \ln(y_i)
  double sw2_x2    = 0.0;  // \sum (w_i)^2 (x_i)^2 
  double sw2       = 0.0;  // \sum (w_i)^2
  double w;
  matrix coord;
  vector auxc(3);
  double phi, r;

  // !!!!!!
  // detection of the minimum negative y coordinate of the domain must be provided
  // in order to avoid problems with negativ or zero logarithm argument
  
  npip=0;
  minx = DBL_MAX;
  maxx = DBL_MIN;
  for(i=0; i<Mt->ne; i++)
  {
    if (Gtm->leso[i]==1)
    {
      //  only active elements are searched
      ipp =  Mt->elements[i].ipp[0][0];
      tnip = Mt->give_tnip(i);
      reallocm(tnip, 3, coord);
      Mt->give_ipcoord_elem(i, coord);
      for (j=0; j<tnip; j++)
      {        
        if (plast_ip[ipp] == id)
        {
          npip++;
          w = Mm->give_consparam(ipp);
          /* 
          if (coord(j, 0) < minx)
            minx = coord(j, 0);
          if (coord(j, 0) > maxx)
            maxx = coord(j, 0);
          sw2_x_lny += w*w*coord[j][0]*log(coord[j][1]);
          sw2_x     += w*w*coord[j][0];
          sw2_lny   += w*w*log(coord[j][1]);
          sw2_x2    += w*w*coord[j][0]*coord[j][0];
          sw2       += w*w;*/
          auxc(0) = coord[j][0] - 10.0;
          auxc(1) = coord[j][1] - 10.0;
          r = normv(auxc);
          phi = acos(auxc(0)/r);
          if (phi < minx)
            minx = phi;
          if (phi > maxx)
            maxx = phi;
          fprintf (Out, "%ld %le %le %le\n", Mm->elip[ipp]+1, phi, r, w*w);
          sw2_x_lny += w*w*phi*log(r);
          sw2_x     += w*w*phi;
          sw2_lny   += w*w*log(r);
          sw2_x2    += w*w*phi*phi;
          sw2       += w*w;
        }
        ipp++;
      }
    }
  }
  if (npip < 2)
  {
    print_err("Cannot determine the slip surface for the given plastic zone %ld\n."
              "At least two different points in the plastic state must be given.", __FILE__, __LINE__, __func__, id);
    return 1;
  }

  a = (sw2_x_lny*sw2_x - sw2_lny*sw2_x2)/(sw2_x*sw2_x-sw2*sw2_x2);
  a = exp(a);
  b = (sw2_lny*sw2_x - sw2_x_lny*sw2)/(sw2_x*sw2_x-sw2*sw2_x2);


  return 0;
}

