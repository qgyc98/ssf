#include "sequent.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "matrix.h"
#include "vector.h"
#include "gtopology.h"
#include "mathem.h"
#include "element.h"
#include "elemhead.h"
#include "intpoints.h"
#include "elemswitch.h"
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <vector>

#ifndef M_PI
 #define M_PI 3.14159265358979323846
#endif



/**
  Function creates list of adjacent integration points and their distances. 
  Results are stored in object Mt (class mechtop):
  - Mt->nadjip[i] - number of adjacent integration points to the i-th integration point
  - Mt->adjip[i] - array of numbers of adjacent integration points to the i-th integration point
  - Mt->dist[i] - array of distances of adjacent integration points
                  to the i-th integration point
  @return The function does not return anything.

  Created by JK,
  Modified by TKo
*/ 
void adjacip (void)
{
  time_t bt,et;
  bt = time (NULL);

  long i,j,k,l,ii,jj,nb,ipp,nip,totnip,nae,nnae,nte,nlmid;
  long *adjelel,*adjelelback,*treatedel,*change;
  double rlim;
  //  long   *tmp_adjip;
  //  double *tmp_dist;
  vector coord(3),auxcoord(3);
  double *ipc; // array of ip coordinates, ipc+3*i = pointer to the array of 3 coordinates of the i-th int. point
  // maximum reached sizes of often reallocated arrays
  long adjelel_max     = 0;
  long adjelelback_max = 0;
  long treatedel_max   = 0;
  long change_max      = 0;
  long max_nadjip      = 0;
  long nipp = 0;

  
  adjelel     = NULL;
  adjelelback = NULL;
  treatedel   = NULL;
  change      = NULL;
  //  total number of integration points
  totnip=Mt->tnip=Mm->tnip;
  Mt->nadjip = new long [totnip];
  Mt->adjip = new long* [totnip];
  Mt->dist = new double* [totnip];
  ipc = new double[3*totnip];
  memset(Mt->nadjip, 0, sizeof(*Mt->nadjip)*totnip);
  memset(Mt->adjip, 0, sizeof(*Mt->adjip)*totnip);
  memset(Mt->dist, 0, sizeof(*Mt->dist)*totnip);
  memset(ipc, 0, sizeof(*ipc)*3*totnip);

  // coordinates of integration points
  for(i=0; i<Gtm->ne; i++){
    nb = Mt->give_nb (i);
    for (ii=0;ii<nb;ii++){
      for (jj=0;jj<nb;jj++){
	ipp=Mt->elements[i].ipp[ii][jj];
	nip = Mt->give_nip (i,ii,jj);	
	for (j=0;j<nip;j++){
          makerefv(coord, ipc+3*ipp, 3);
          ipcoord(i, ipp, ii, jj, coord);
          ipp++;
        }
      }
    }
  }
  
  fprintf(stdout, "\n Adjacent integration points\n");
  //  max_nadjip = give_max_adjacip();
  //  fprintf(stdout, " guess of maximum number of integration point is %ld\n", max_nadjip);

  /*
  struct cmpfunc{
    static int compare (const void * a, const void * b) {
      return ( *(int*)a - *(int*)b );
    }
  };
  FILE *ftrel = fopen("treatedel.log", "wt");
  fprintf(ftrel, "\n\nList of treated elements of adjac. int. points:\n");
  */
  
  for (i=0;i<Gtm->ne;i++){
    nb = Mt->give_nb (i);
    if(i%100 == 0)
      fprintf(stdout, "\n Adjacent integration points on element %ld, max_nadjip=%ld", i+1,max_nadjip);
    for (ii=0;ii<nb;ii++){
      for (jj=0;jj<nb;jj++){
	
	ipp=Mt->elements[i].ipp[ii][jj];
        if (!(Mm->ip[ipp].hmt & 2)){
	  //  integration point number ipp does not contain nonlocal material model
	  Mt->nadjip[ipp]=0;
	  continue;
	}
	nip = Mt->give_nip (i,ii,jj);
	
	for (j=0;j<nip;j++){
	  
	  if (!(Mm->ip[ipp].hmt & 2)){
	    //  integration point number j does not contain nonlocal material model
	    continue;
	  }
	  
	  nipp++;
	  
          nlmid = Mm->givenonlocid(ipp);
	  //  radius of nonlocal volume
	  rlim = Mm->nonlocradius (ipp,nlmid);
	  
	  //  coordinates of integration points
          makerefv(coord, ipc+3*ipp, 3);
	  //ipcoord (i,ipp,ii,jj,coord);
	  
	  //  number of adjacent integration points determination
	  Mt->nadjip[ipp]=0;
	  //  number of adjacent elements to required element
	  nae=Gtm->nadjelel[i];
	  nte=nae;
	  //  list of adjacent elements to required element
          reallocate(adjelel, nae, adjelel_max);
          reallocate(treatedel, nae, treatedel_max);
	  for (k=0;k<nae;k++){
	    adjelel[k]=Gtm->adjelel[i][k];
	    treatedel[k]=adjelel[k];
	  }
	  
	  //  computes the number of integration points in a circle (ball) with the radies rlim
	  //  these integration points will be used in nonlocal models
	  //  the number of integration points are stored in the array Mt->nadjip[ipp]
	  do{
	    in_dist (ipp,coord,ipc,adjelel,nae,rlim);
	    
            reallocate(adjelelback, nae, adjelelback_max);
	    
	    newnadjelel (ipp,adjelelback,adjelel,nae,nnae);
	    if (nnae>0){
	      
              reallocate(adjelel, nnae, adjelel_max);
	      
	      newadjelel (adjelelback,adjelel,nae,nnae,treatedel,nte);
	      
              reallocate(change, nte, change_max);
	      
	      for (k=0;k<nte;k++){
		change[k]=treatedel[k];
	      }
	      
              reallocate(treatedel, nte+nnae, treatedel_max);
	      
	      for (k=0;k<nte;k++){
		treatedel[k]=change[k];
	      }
	      l=0;
	      for (k=nte;k<nte+nnae;k++){
		treatedel[k]=adjelel[l];  l++;
	      }
	      nte+=nnae;
	      nae=nnae;
	    }
	    else{}
	  }
	  while (nnae!=0);

          /*
          qsort(treatedel, nte, sizeof(*treatedel), cmpfunc::compare);
          fprintf(ftrel, "ipp=%ld, nte=%ld:", ipp, nte);
          for(long kk=0; kk<nte; kk++){
            fprintf(ftrel, " %ld", treatedel[kk]);
          }
          fprintf(ftrel, "\n");
          */
          
	  //  computes the number of integration points in a circle (ball) with the radies rlim
	  //  these integration points will be used in nonlocal models
	  //  the number of integration points are stored in the array Mt->nadjip[ipp]
	  
	  //  adjacent integration points determination
          if (max_nadjip < Mt->nadjip[ipp])
            max_nadjip = Mt->nadjip[ipp];
	  Mt->adjip[ipp] = new long [Mt->nadjip[ipp]];
	  Mt->dist[ipp] = new double [Mt->nadjip[ipp]];
	  //	  Mt->adjip[ipp] = new long [max_nadjip];
	  //	  Mt->dist[ipp]  = new double [max_nadjip];
	  Mt->nadjip[ipp]=0;
	  nae=Gtm->nadjelel[i];
	  nte=nae;
          reallocate(adjelel, nae, adjelel_max);
          reallocate(treatedel, nae, treatedel_max);
	  for (k=0;k<nae;k++){
	    adjelel[k]=Gtm->adjelel[i][k];
	    treatedel[k]=adjelel[k];
	  }
	  
	  do{
	    dist (ipp,max_nadjip,coord,ipc,adjelel,nae,rlim);
	    
            reallocate(adjelelback, nae, adjelelback_max);
	    
	    newnadjelel (ipp,adjelelback,adjelel,nae,nnae);
	    if (nnae>0){
	      
              reallocate(adjelel, nnae, adjelel_max);
	      
	      newadjelel (adjelelback,adjelel,nae,nnae,treatedel,nte);
	      
              reallocate(change, nte, change_max);
	      
	      for (k=0;k<nte;k++){
		change[k]=treatedel[k];
	      }
	      
              reallocate(treatedel, nte+nnae, treatedel_max);
	      
	      for (k=0;k<nte;k++){
		treatedel[k]=change[k];
	      }
	      l=0;
	      for (k=nte;k<nte+nnae;k++){
		treatedel[k]=adjelel[l];  l++;
	      }
	      nte+=nnae;
	      nae=nnae;
	    }
	  }
	  while (nnae!=0);
	  ipp++;
	}
      }
    }
  }

  delete [] adjelelback;
  delete [] adjelel;
  delete [] treatedel;
  delete [] change;
  delete [] ipc;
  //fclose(ftrel);  
  
  // compression of arrays with adjacent integration points and distances
  /*  for (ipp=0; ipp<totnip; ipp++)
      {	  
      tmp_dist  = Mt->dist[ipp];
      Mt->dist[ipp]  = new double [Mt->nadjip[ipp]];
      memcpy(Mt->dist[ipp], tmp_dist, sizeof(*tmp_dist)*Mt->nadjip[ipp]);
      delete [] tmp_dist;
      tmp_adjip = Mt->adjip[ipp];
      Mt->adjip[ipp] = new long [Mt->nadjip[ipp]];
      memcpy(Mt->adjip[ipp], tmp_adjip, sizeof(*tmp_adjip)*Mt->nadjip[ipp]);
      delete [] tmp_adjip;
      }
  */

  /*
  for (ipp=0; ipp<Mm->tnip; ipp++){
    long nadjip = Mt->nadjip[ipp];
    qsort(Mt->adjip[ipp], nadjip, sizeof(*Mt->adjip[ipp]), cmpfunc::compare);
  }*/

  et = time (NULL);
  fprintf (Out,"\n\n time of adajcip             %ld", long(et-bt));
  fprintf (stdout,"\n\n time of adajcip             %ld", long(et-bt));
  fflush(Out);
}



/**
  The function saves adjacent integration points on the elements.
  The backup is controlled by the hdb controller.

  @return The function does not return anything.

  Created by Tomas Koudelka, 2008 
*/
void save_adjacip()
{
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];

  switch(Mp->hdbcont.hdbtype)
  {
    case hdbs_nonloc:
    case hdbrs_nonloc:
    case hdbrs_single:
    case hdbrs_multiple:
    case hdbs_single:
    case hdbs_multiple:
      sprintf(name, "%s.adjacip.bac",Mp->hdbcont.hdbnames);
      if (Mp->hdbcont.hdbfmts == text)
        aux = fopen(name,"wt");
      else
        aux = fopen(name,"wb");
      if (aux==NULL)
      {
        sprintf(emsg, "cannot open backup file %s", name);
        print_err(emsg, __FILE__, __LINE__, __func__);
        abort();
      }
      if (Mp->hdbcont.hdbfmts == text)
        save_adjacip_txt(aux);
      else 
        save_adjacip_bin(aux);
      break;
    default:
      break;
  }    
}



/**
  The function restores adjacent integration points on the elements.
  The process is controlled by the hdb controller.

  @retval 0 - restorage has not been performed
  @retval 1 - successfull restorage

  Created by Tomas Koudelka, 2008 
*/
long restore_adjacip()
{
  FILE *aux;
  char name[FNAMELEN+20];
  char emsg[FNAMELEN+100];

  switch(Mp->hdbcont.hdbtype)
  {
    case hdbr_nonloc:
    case hdbrs_nonloc:
    case hdbrs_single:
    case hdbrs_multiple:
    case hdbr_single:
    case hdbr_multiple:
      sprintf(name, "%s.adjacip.bac",Mp->hdbcont.hdbnamer);
      if (Mp->hdbcont.hdbfmtr == text)
        aux = fopen(name,"rt");
      else
        aux = fopen(name,"rb");
      if (aux==NULL)
      {
        sprintf(emsg, "cannot open backup file %s", name);
        print_err(emsg, __FILE__, __LINE__, __func__);
        abort();
      }
      if (Mp->hdbcont.hdbfmtr == text)
        restore_adjacip_txt(aux);
      else 
        restore_adjacip_bin(aux);
      break;
    default:
      return 0;
  }    
  return 1;
}



/**
  The function saves adjacent integration points on the elements
  in the text format. 
  
  @param aux - pointer to the opened file
  
  @return The function does not return anything.

  Created by Tomas Koudelka, 2008 
*/
void save_adjacip_txt(FILE *aux)
{
  long i,j;
  if (Mespr)
    fprintf(stdout, "\n Saving adjacent integration points to backup file");
  for (i=0;i<Mm->tnip;i++)
  {
    fprintf (aux,"%ld %ld",i,Mt->nadjip[i]);
    for (j=0;j<Mt->nadjip[i];j++)
      fprintf (aux," %ld %le",Mt->adjip[i][j],Mt->dist[i][j]);
    fprintf(aux, "\n");
  }
}



/**
  The function saves adjacent integration points on the elements
  in the binary format. 
  
  @param aux - pointer to the opened file

  @return The function does not return anything.

  Created by Tomas Koudelka, 2008 
*/
void save_adjacip_bin(FILE *aux)
{
  long i,j;
  if (Mespr)
    fprintf(stdout, "\n Saving adjacent integration points to backup file");
  for (i=0;i<Mm->tnip;i++)
  {
    fwrite(&i, sizeof(i), 1, aux);
    fwrite(Mt->nadjip+i, sizeof(*Mt->nadjip), 1, aux);
    for (j=0;j<Mt->nadjip[i];j++)
    {
      fwrite(Mt->adjip[i]+j, sizeof(*Mt->adjip[i]), 1, aux);
      fwrite(Mt->dist[i]+j, sizeof(*Mt->dist[i]), 1, aux);
    }
  }
}



/**
  The function restores adjacent integration points on the elements
  from the text file.

  @param aux - pointer to the opened file

  @return The function does not return anything.

  Created by Tomas Koudelka, 2008 
*/
void restore_adjacip_txt(FILE *aux)
{
  long i,j,ipp;
  Mt->nadjip = new long[Mm->tnip];
  memset(Mt->nadjip, 0, sizeof(*Mt->nadjip)*Mm->tnip);
  Mt->adjip = new long*[Mm->tnip];
  memset(Mt->adjip, 0, sizeof(*Mt->adjip)*Mm->tnip);
  Mt->dist = new double*[Mm->tnip];
  memset(Mt->dist, 0, sizeof(*Mt->dist)*Mm->tnip);
  if (Mespr)
    fprintf(stdout, "\n Restoring adjacent integration points from backup file");
  for (i=0;i<Mm->tnip;i++)
  {
    fscanf (aux, "%ld%ld", &ipp, Mt->nadjip+i);
    if ((ipp < 0) || (ipp >= Mm->tnip))
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }    
    Mt->adjip[i] = new long[Mt->nadjip[i]];
    memset(Mt->adjip[i], 0, sizeof(*Mt->adjip[i])*Mt->nadjip[i]);
    Mt->dist[i]   = new double[Mt->nadjip[i]];
    memset(Mt->dist[i], 0, sizeof(*Mt->dist[i])*Mt->nadjip[i]);
    for (j=0; j<Mt->nadjip[i]; j++)
      fscanf (aux, "%ld%le", Mt->adjip[i]+j, Mt->dist[i]+j);
  }
}



/**
  The function restores adjacent integration points on the elements
  from the binary file.

  @param aux - pointer to the opened file

  @return The function does not return anything.

  Created by Tomas Koudelka, 2008 
*/
void restore_adjacip_bin(FILE *aux)
{
  long i,j,ipp;
  Mt->nadjip = new long[Mm->tnip];
  memset(Mt->nadjip, 0, sizeof(*Mt->nadjip)*Mm->tnip);
  Mt->adjip = new long*[Mm->tnip];
  memset(Mt->adjip, 0, sizeof(*Mt->adjip)*Mm->tnip);
  Mt->dist = new double*[Mm->tnip];
  memset(Mt->dist, 0, sizeof(*Mt->dist)*Mm->tnip);
  if (Mespr)
    fprintf(stdout, "\n Restoring adjacent integration points from backup file");
  for (i=0;i<Mm->tnip;i++)
  {
    fread(&ipp, sizeof(i), 1, aux);
    fread(Mt->nadjip+i, sizeof(*Mt->nadjip), 1, aux);
    if ((ipp < 0) || (ipp >= Mm->tnip))
    {
      print_err("invalid integration point number", __FILE__, __LINE__, __func__);
      abort();
    }
    Mt->adjip[i] = new long[Mt->nadjip[i]];
    memset(Mt->adjip[i], 0, sizeof(*Mt->adjip[i])*Mt->nadjip[i]);
    Mt->dist[i]   = new double[Mt->nadjip[i]];
    memset(Mt->dist[i], 0, sizeof(*Mt->dist[i])*Mt->nadjip[i]);
    for (j=0;j<Mt->nadjip[i];j++)
    {
      fread(Mt->adjip[i]+j, sizeof(*Mt->adjip[i]), 1, aux);
      fread(Mt->dist[i]+j, sizeof(*Mt->dist[i]), 1, aux);
    }
  }
}



/**
  Function estimates maximum number of adjacent integration points.

  @retval Function returns estimated maximum number of adjacent integration points. 

  Created by Tomas Koudelka 8.2008
*/
long give_max_adjacip()
{
  long i,tnip,max_tnip,ipp,ret,nlmid;
  long max_1d, max_2d, max_3d;
  double vol, area, length, min_vol, min_area, min_length;
  double r, max_r;

  ret = 0;
  max_tnip = 0;
  max_r = 0.0;
  min_vol = min_area = min_length = DBL_MAX;
  for (i=0; i<Mt->ne; i++)
  {
    ipp=Mt->elements[i].ipp[0][0];
    
    if (!(Mm->ip[ipp].hmt & 2))
    //  integration point number ipp does not contain nonlocal material model
      continue;
    
    nlmid = Mm->givenonlocid(ipp);
    //  maximum radius of nonlocal volume
    r = Mm->nonlocradius (ipp,nlmid);
    if (r > max_r)
      max_r = r;
    // maximum number of integration points
    tnip = Mt->give_tnip(i);
    if (tnip > max_tnip)
      max_tnip = tnip;
    // minimum 
    switch (Mt->give_dimension(i))
    {
      case 1:
        length = fabs(Mt->give_length(i));
        if (length < min_length)
          min_length = length;
        break;
      case 2:
        area = fabs(Mt->give_area(i));
        if (area < min_area)
          min_area = area;
        break;
      case 3:
        vol = fabs(Mt->give_volume(i));
        if (vol < min_vol)
          min_vol = vol;
        break;
      default:
        print_err("unknown dimension of element", __FILE__, __LINE__, __func__);
    }
  }
  if (min_length == 0.0)
    print_err("minimum length of element is zero ", __FILE__, __LINE__, __func__);
  if (min_area == 0.0)
    print_err("minimum area of element is zero ", __FILE__, __LINE__, __func__);
  if (min_vol == 0.0)
    print_err("minimum volume of element is zero ", __FILE__, __LINE__, __func__);
    
  // maximum number of 1d elements
  max_1d = 15*(long(2.0*max_r/min_length));
  // maximum number of 2d elements
  max_2d = 15*(long(M_PI*max_r*max_r/min_area));
  // maximum number of 3d elements
  max_3d = 15*(long(4.0/3.0*M_PI*max_r*max_r*max_r/min_vol));
  if (max_1d > max_2d)
    ret = max_1d;
  else
    ret = max_2d;
  if (max_3d > ret)
    ret = max_3d;

  return ret;
}



/**
   Function performs "smart" allocation and reallocation of memory for
   array of long integer values.

   @param array - pointer to array which should be allocated
   @param req_size - required size of array which should be allocated
   @param capacity - maximum size of allocated array in history

   @return The function does not return anything.

   Created by Tomas Koudelka 8.2008
*/
void reallocate(long *&array, long req_size, long &capacity)
{
  if (req_size == 0)
    return;
  if (req_size > capacity)
  {
    //    fprintf(Out,"\n Reallocation req_size=%ld, capacity=%ld", req_size, capacity);
    capacity = 2*req_size;
    if (array)
      delete [] array;
    array = new long [capacity];
  }
}



/**
  Function calculates number of non-eliminated adjacent elements.
  Function is used in function adjacip.
   
  @param ipp - integration point pointer
  @param adjelelback - auxiliary array of numbers of adjacent elements
  @param adjelel - array of numbers of adjacent elements
  @param nae - number of adjacent elements
  @param nnae - new number of adjacent elements

  @return The function does not return anything.
 
  Created by JK,
*/
void newnadjelel (long /*ipp*/,long *adjelelback,long *adjelel,long &nae,long &nnae)
{
  long i,nee;
  
  //  number of eliminated elements
  nee=0;
  for (i=0;i<nae;i++){
    if (adjelel[i]<0)  nee++;
  }
  
  if (nae!=nee){
    nnae=0;
    for (i=0;i<nae;i++){
      adjelelback[i]=adjelel[i];
      if (adjelel[i]>-1){
	nnae+=Gtm->nadjelel[adjelel[i]];
      }
    }
  }
  else{
    nnae=0;
  }
}



/**
  Function creates list of adjacent elements to required element.
  Function is used in function adjacip.
   
  @param adjelelback - auxiliary array of adjacent elements
  @param adjelel - array of adjacent elements
  @param nae - number of adjacent elements
  @param nnae - new number of adjacent elements
  @param treatedel - array of treated elements
  @param nte - number of treated elements
   
  @return The function does not return anything.

  Created by JK,
*/
void newadjelel (long *adjelelback,long *adjelel, long &nae, long &nnae, long *treatedel, long &nte)
{
  long i,j,k,l,prev,min;
  
  j=0;
  for (i=0;i<nae;i++){
    if (adjelelback[i]>-1){
      for (k=0;k<Gtm->nadjelel[adjelelback[i]];k++){
	adjelel[j]=Gtm->adjelel[adjelelback[i]][k];
	j++;
      }
    }
  }
  
  //  sorting
  prev=-1;
  for (i=0;i<nnae;i++){
    min=LONG_MAX;
    for (j=i;j<nnae;j++){
      if (min>adjelel[j]){
	min=adjelel[j];  k=j;
      }
    }
    if (min==prev){
      nnae--;  i--;
      adjelel[k]=adjelel[nnae];
    }
    else{
      l=adjelel[i];
      adjelel[i]=min;
      adjelel[k]=l;
      prev=min;
    }
  }
  
  for (i=0;i<nnae;i++){
    k=adjelel[i];
    for (j=0;j<nte;j++){
      if (treatedel[j]==k){
	nnae--;
	adjelel[i]=adjelel[nnae];
	i--;
	break;
      }
    }
  }

}



/**
  Function computes number of integration points which are closer
  to required integration point (ipp) than rlim.
   
  Function is used in adjacip.
   
  @param[in] ipp - integration point pointer
  @param[in] coord - coordinates of one of integration points
  @param[in] ipc - array of coordinates of all integration points, 
                   ipc+3*i = pointer to array of 3 coordinates of i-th integration point
  @param[in,out] adjelel - array of adjacent elements
  @param[in] nae - number of adjacent elements
  @param[in] rlim - limit distance
   
  @return The function does not return anything.
  
  Created by JK,
  Modified by TKo, 04.2018
*/
void in_dist (long ipp, vector &coord, double *ipc, long *adjelel, long nae, double rlim)
{
  long i,j,ri,ci,lj,uj,ii,eid,nip,nb;
  double r;
  vector auxcoord(3);
  
  for (i=0;i<nae;i++){
    eid = adjelel[i];
    nb = Mt->give_nb(eid);
    ii=0;
    for (ri=0; ri<nb; ri++){
      for (ci=0; ci<nb; ci++){
        nip = Mt->give_nip (eid, ri, ci);
        lj = Mt->elements[eid].ipp[ri][ci];
        uj = lj+nip;
        for (j=lj;j<uj;j++){
          if (Mm->ip[j].hmt & 2){
            //  integration point number j contains nonlocal material model
            makerefv(auxcoord, ipc+j*3, 3); 
            //ipcoord (eid, j, ri, ci, auxcoord);
            r = length (coord, auxcoord);
            if (r<rlim){
              Mt->nadjip[ipp]++;
              ii++;
            }
          }
        }
      }
    }
    if (ii==0)  
      adjelel[i]=-1-adjelel[i];
  }
}



/**
  Function computes distances of integration points which are closer
  to required integration point (ipp) than rlim.
   
  Function is used in adjacip.
   
  @param[in] ipp - integration point pointer
  @param[in] max_nadjip - maximum number of adjacent integration points
  @param[in] coord - coordinates of one of integration points
  @param[in] ipc - array of coordinates of all integration points, 
                   ipc+3*i = pointer to array of 3 coordinates of i-th integration point  
  @param[in,out] adjelel - array of adjacent elements
  @param[in] nae - number of adjacent elements
  @param[in] rlim - limit distance
   
  @return The function does not return anything.

  Created by JK,
  Modified by TKo, 04.2018
*/
void dist (long ipp,long /*max_nadjip*/,vector &coord, double *ipc, long *adjelel, long nae, double rlim)
{
  long i,j,ri,ci,lj,uj,ii,eid,nip,nb;
  double r;
  vector auxcoord(3);
  
  for (i=0;i<nae;i++){
    eid = adjelel[i];
    nb = Mt->give_nb(eid);
    ii=0;    
    for (ri=0; ri<nb; ri++){
      for (ci=0; ci<nb; ci++){
        nip = Mt->give_nip (eid,ri,ci);
        lj = Mt->elements[eid].ipp[ri][ci];
        uj = lj+nip;
    
        for (j=lj;j<uj;j++){
          if (Mm->ip[j].hmt & 2){
            //  integration point number j contains nonlocal material model
            makerefv(auxcoord, ipc+j*3, 3); 
            //ipcoord (eid,j,ri,ci,auxcoord);
            r = length (coord,auxcoord);
            if (r<rlim){
              Mt->adjip[ipp][Mt->nadjip[ipp]]=j;
              Mt->dist[ipp][Mt->nadjip[ipp]]=r;
              Mt->nadjip[ipp]++;
              ii++;
            }
          }
        }
      }
    }
    if (ii==0)  adjelel[i]=-1-adjelel[i];
  }
}



/**
  The function creates list of adjacent integration points and their distances. 
  Results are stored in arrays of object Mt (class mechtop):
  - Mt->nadjip[i] - number of adjacent integration points to the i-th integration point
  - Mt->adjip[i] - array of numbers of adjacent integration points to the i-th integration point
  - Mt->dist[i] - array of distances of adjacent integration points
                  to the i-th integration point
  @return The function does not return anything, but allocates and changes arrays of Mt object.

  Created by TKo, 12.2023
*/ 
void adjacip2 (void)
{
  time_t bt,et;
  bt = time (NULL);

  long i, j, ii, jj, nb, ipp, nip, totnip, nae, nnae, nlmid, nadjip;

  ivector adjelel; // list of adjacent elements used in the actual lookup pass for adjacent int. points

  ivector adjelelback; // backup of the adjacent element list from the passed lookup

  // list of adjacent int. points, std::vector allows for automatic growing while pushing back new elements
  std::vector<long> ipl;

  // list of distances of adjacent int. points, std::vector allows for automatic growing while pushing back new elements
  std::vector<double> ipd;

  // list of already treated elements in one adjacent int. point lookup,
  // std::set allows for automatic duplicate element removal
  std::set<long> treatedel;  

  double rlim; // limit radius of the adjacent int. point  search
  vector coord; // coordinates of the actual int. point whose neighbours are looked for

  vector ipc(3*Mm->tnip); // array of int. point coordinates, ipc+3*i = pointer to the array of 3 coordinates of the i-th int. point

  // maximum reached sizes of often reallocated arrays
  long max_nadjip = 0;

  // total number of integration points involved in the neighbour search, i.e. number of int. points with a nonlocal material model
  long nipp = 0;


  //  total number of integration points
  totnip = Mt->tnip = Mm->tnip;
  Mt->nadjip = new long [totnip];
  Mt->adjip = new long* [totnip];
  Mt->dist = new double* [totnip];
  memset(Mt->nadjip, 0, sizeof(*Mt->nadjip)*totnip);
  memset(Mt->adjip, 0, sizeof(*Mt->adjip)*totnip);
  memset(Mt->dist, 0, sizeof(*Mt->dist)*totnip);

  // calculate coordinates of integration points
  for(i=0; i<Gtm->ne; i++){
    nb = Mt->give_nb (i);
    for (ii=0;ii<nb;ii++){
      for (jj=0;jj<nb;jj++){
	ipp=Mt->elements[i].ipp[ii][jj];
	nip = Mt->give_nip (i,ii,jj);	
	for (j=0;j<nip;j++){
          makerefv(coord, &ipc[3*ipp], 3);
          ipcoord(i, ipp, ii, jj, coord);
          ipp++;
        }
      }
    }
  }
  
  fprintf(stdout, "\n Adjacent integration points\n");

  // struct cmpfunc{
  //   static int compare (const void * a, const void * b) {
  //     return ( *(int*)a - *(int*)b );
  //   }
  // };
  // FILE *ftrel = fopen("treatedel.log", "wt");
  // fprintf(ftrel, "\n\nList of treated elements of adjac. int. points:\n");

  for (i=0;i<Gtm->ne;i++){
    if(i%100 == 0)
      fprintf(stdout, "\n Adjacent integration points on element %ld, max_nadjip=%ld", i+1,max_nadjip);
    ipp=Mt->elements[i].ipp[0][0];
    if (!(Mm->ip[ipp].hmt & 2)){
      //  integration point number ipp does not contain nonlocal material model
      Mt->nadjip[ipp]=0;
      continue;
    }
    nip = Mt->give_tnip(i);
	
    for (j=0;j<nip;j++){
      nipp++;

      // index of the nonlocal material model
      nlmid = Mm->givenonlocid(ipp);
      //  radius of nonlocal volume
      rlim = Mm->nonlocradius (ipp,nlmid);
	  
      // reference to coordinates of the actual integration point
      makerefv(coord, &ipc[3*ipp], 3);
	  
      //  number of adjacent integration points determination
      Mt->nadjip[ipp]=0;
      
      //  number of adjacent elements to required element
      nae=Gtm->nadjelel[i];

      //  create initial list of adjacent elements to required element including the given element
      reallocv(nae, adjelel);
      copyv(Gtm->adjelel[i], adjelel.a, nae);
      // initialize list of treated elements
      treatedel.clear(); // erase all elements
      treatedel.insert(adjelel.a, adjelel.a + adjelel.n);  // copy actual list of adjacent elements
      // empty list of adjacent int. points
      ipl.clear();
      // empty list of distances of adjacent int. point
      ipd.clear();
	  
      //  computes the number of integration points in a circle (ball) with the radius rlim
      //  these integration points will be used in nonlocal models
      //  the number of adjacent integration points are stored in the array Mt->nadjip[ipp]
      do{
        in_dist2(ipp, coord, ipc, rlim, adjelel, ipl, ipd);
        ipl.capacity();
        ipd.capacity();
	    
        reallocv(nae, adjelelback);
	    
        newnadjelel2(adjelel, adjelelback, nnae);
        if (nnae>0){
          newadjelel2(adjelelback, treatedel, adjelel, nnae);
          nae=nnae;
        }
      }
      while (nnae!=0);
	  
      // storage of found adjacent integration points and their distances in mechtop arrays
      nadjip = Mt->nadjip[ipp];
      if (max_nadjip < nadjip)
        max_nadjip = nadjip;
      if (Mt->adjip[ipp] == NULL)
        Mt->adjip[ipp] = new long[nadjip];
      else{
        print_err("multiple allocation of adjip[%ld], eid=%ld", __FILE__, __LINE__, __func__, ipp, i+1);
        abort();
      }
        
      Mt->dist[ipp]  = new double[nadjip];
      copyv(ipl.data(), Mt->adjip[ipp], nadjip); // ipl.data() returns the pointer to the array of std::vector elements
      copyv(ipd.data(), Mt->dist[ipp],  nadjip);
      
      // fprintf(ftrel, "ipp=%ld, nte=%ld:", ipp, treatedel.size());
      // for (std::set<long>::iterator it = treatedel.begin(); it != treatedel.end(); it++){
      //   fprintf(ftrel, " %ld", *it);
      // }
      // fprintf(ftrel, "\n");

      ipp++;
    }
  }

  /*
  for (ipp=0; ipp<Mm->tnip; ipp++){
    long nadjip = Mt->nadjip[ipp];
    qsort(Mt->adjip[ipp], nadjip, sizeof(*Mt->adjip[ipp]), cmpfunc::compare);
  }*/
  
  et = time (NULL);
  fprintf (Out,"\n\n time of adajcip             %ld", long(et-bt));
  fprintf (stdout,"\n\n time of adajcip             %ld", long(et-bt));
  fflush(Out);
  // fclose(ftrel);  
}



/**
  The function calculates a number of non-eliminated adjacent elements.
  It is used in function adjacip.
   
  @param[in] adjelel - array of identfiers of adjacent elements
  @param[out] adjelelback - auxiliary array of identfiers of adjacent elements, it will contain backup of adjelel array
  @param[out] nnae - new number of adjacent elements

  @return The function does not return anything.
 
  Created by TKo 12.2023, according to the original JK algorithm
*/
void newnadjelel2(const ivector &adjelel, ivector &adjelelback, long &nnae)
{
  long i, nee, nae = adjelel.n;
  
  // number of eliminated elements, i.e. no adjacent int. points were found on these elements
  nee=0;
  for(i=0; i<nae; i++){
    if (adjelel[i] < 0)  nee++;
  }
  nnae = 0;
  if (nae!=nee){
    for (i=0; i<adjelel.n; i++){
      adjelelback[i] = adjelel[i];
      if (adjelel[i] > -1){
	nnae += Gtm->nadjelel[adjelel[i]];
      }
    }
  }
}



/**
  Function creates list of adjacent elements to required element.
  Function is used in function adjacip.
   
  @param[in]  adjelelback - auxiliary array of adjacent elements from the previous lookup pass
  @param[in,out]  treatedel - array of already treated elements, it may be actualized by new set of adjacent elements
  @param[out] adjelel - array of new adjacent elements for the neighbour lookup
  @param[out] nnae - new number of adjacent elements
   
  @return The function does not return anything.

  Created by TKo 10.2023, according to original JK algorithm
*/
void newadjelel2(const ivector &adjelelback, std::set<long> &treatedel, ivector &adjelel, long &nnae)
{
  long i, j, eid, nadjel, aeid;
  std::set<long> new_adjelel;

  // create list of new adjacent elements to already passed (adjacent) elements
  for (i=0; i<adjelelback.n; i++){
    eid = adjelelback[i];
    if (eid > -1){
      nadjel = Gtm->nadjelel[eid];
      // insert all adjacent elements of eid-th element, duplicates are removed automatically in std::set
      for (j=0; j<nadjel; j++){
        aeid = Gtm->adjelel[eid][j];
        if (treatedel.count(aeid) == 0){ // std::set.count returns a number of its elements matching aeid, i.e. 0 or 1
          new_adjelel.insert(aeid); // insert just non-treated elements
          treatedel.insert(aeid);   // add new element for the treatment
        }
      }
    }
  }

  nnae = new_adjelel.size();
  reallocv(nnae, adjelel);
  // copy new adjacent elements to array adjelel
  decltype(new_adjelel)::iterator it = new_adjelel.begin(); // get iterator (pointer) to the first element of new_adjelel
  for(i=0; i<nnae; it++, i++)
    adjelel[i] = *it;
}



/**
  Function computes number of integration points which are closer
  to required integration point (ipp) than rlim.
   
  Function is used in adjacip2.
   
  @param[in] ipp - integration point pointer
  @param[in] coord - coordinates of one of integration points
  @param[in] ipc - array of coordinates of all integration points, 
                   ipc+3*i = pointer to array of 3 coordinates of i-th integration point
  @param[in] rlim - limit distance
  @param[in,out] adjelel - array of adjacent elements
  @param[in,out] - list of adjacent integration points.
  @param[in,out] - list of distances of adjacent integration points.
   
  @return The function does not return anything.
  
  Created by JK,
  Modified by TKo, 04.2018
*/
void in_dist2 (long ipp, const vector &coord, const vector &ipc, double rlim, ivector &adjelel, std::vector<long> &ipl, std::vector<double> &ipd)
{
  long i, j, lj, uj, eid, nip, count;
  double r;
  vector auxcoord;
  
  for (i=0;i<adjelel.n;i++){
    eid = adjelel[i];
    count = 0; 
    nip = Mt->give_tnip(eid);
    lj  = Mt->elements[eid].ipp[0][0];
    uj  = lj+nip;
    for (j=lj; j<uj; j++){
      if (Mm->ip[j].hmt & 2){// integration point number j has a nonlocal material model
        // make reference to coordinates of the j-th int. point
        makerefv(auxcoord, &ipc[j*3], 3); 
        r = length(coord, auxcoord); // calculate distance of the j-th int. point
        if (r<rlim){         // j-th int. point is closer than the limit distance
          Mt->nadjip[ipp]++; // increase the number of adjacent int. points
          ipl.push_back(j);  // store the adjacent int. point number
          ipd.push_back(r);  // store the adjacent int. point distance
          count++;           // increase number of found adjacent int. point
        }
      }
    }
    if (count==0)  // if no adjacent int. points were found on the actual element, mark it by negative value
      adjelel[i]=-1-adjelel[i];
  }
}



void print_adjacip(FILE *out)
{
  fprintf(out, "\nADJACENT IP ACCORDING TO NONLOCAL MODELS\n\n");
  long eid = -1;
  for(long i=0; i<Mm->tnip; i++){
    long nadjip = Mt->nadjip[i];
    if (Mm->elip[i] != eid){
      eid = Mm->elip[i];
      fprintf(out, "Element %ld, ipp[0][0]=%ld: ----------------------------------\n", eid+1, Mt->elements[eid].ipp[0][0]);
    }
    fprintf(out, "IPP %7ld, nadjip %4ld:", i, nadjip);
    for(long j=0; j<nadjip; j++){
      fprintf(out, " %ld", Mt->adjip[i][j]);
    }
    fprintf(out, "\n");
  }
}
