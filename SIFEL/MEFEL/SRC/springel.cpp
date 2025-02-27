#include "springel.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "globmat.h"
#include "element.h"
#include "matrix.h"
#include "vector.h"
#include "gtopology.h"
#include <math.h>



/**
  The constructor inializes attributes to zero values and allocates
  internal data arrays.
  
  Created by Tomas Koudelka, 12.8.2001
*/
springel::springel (void)
{
  long i;

  nne=1;  ndofe=-1;  tncomp=1;  napfun=0;
  ssst = bar;

  nb=1;

  ncomp = new long [nb];
  ncomp[0]=1;

  cncomp = new long [nb];
  cncomp[0]=0;

  nip = new long* [nb];
  intordsm = new long* [nb];
  for (i=0;i<nb;i++){
    nip[i] = new long [nb];
    intordsm[i] = new long [nb];
  }

  nip[0][0]=1;
  tnip=1;
  
  intordsm[0][0]=1;
}



/**
  The destructor is defined only for the formal purposes.
  
  Created by Tomas Koudelka, 12.8.2001
*/
springel::~springel (void)
{
  long i;

  for (i=0;i<nb;i++){
    delete [] nip[i];
    delete [] intordsm[i];
  }
  delete [] nip;
  delete [] intordsm;

  delete [] cncomp;
  delete [] ncomp;
}





/**
  The function returns correct number of dofs on the spring for given element eid and also setups 
  ndofe attribute to this value.
  
  @param eid - element id

  @return The function returns number of DOFs.

  Created by Tomas Koudelka, 12.8.2001
*/
long springel::give_ndofe (long eid)
{  
  if ((Gtm->give_ndofe(eid) == 0) && (Gtm->give_nne(eid) > 0))
  {  
    ivector enodes(1);
    Mt->give_elemnodes (eid, enodes);
    Gtm->gelements[eid].ndofe = Mt->give_ndofn(enodes[0]);
  }

  return Gtm->give_ndofe(eid);
}



/**
  The function computes stiffness %matrix of given block. The type
  of element determines direction of the spring support and thus the
  stiffness contribution is stored to the appropriate position of the
  matrix sm.

  @param eid - element id.
  @param ri  - block row id
  @param ci  - block column id
  @param sm  - stiffness %matrix where the results are stored (output)

  @return The function returns assmebled block of stiffness %matrix in the parameter sm.

  Created by Tomas Koudelka, 12.8.2001
*/
void springel::stiffness_matrix (long eid,long ri,long ci,matrix &sm)
{
  matrix d(1, 1);
  Mm->matstiff (d,Mt->elements[eid].ipp[ri+0][ci+0]);
  switch (Mt->elements[eid].te)
  {
    case spring_1:
      sm[0][0] = d[0][0];
      break;
    case spring_2:
      sm[1][1] = d[0][0];
      break;
    case spring_3:
      sm[2][2] = d[0][0];
      break;
    case spring_4:
      sm[3][3] = d[0][0];
      break;
    case spring_5:
      sm[4][4] = d[0][0];
      break;
    case spring_6:
      sm[5][5] = d[0][0];
      break;
    default:{break;}
  }
}



/**
  The function computes resulting stiffness %matrix

  @param eid - element id.
  @param sm  - stiffness %matrix where the results are stored (output)

  @return The function returns assmebled resulting stiffness %matrix in the parameter sm.

  Created by Tomas Koudelka, 12.8.2001
*/
void springel::res_stiffness_matrix (long eid,matrix &sm)
{
  stiffness_matrix (eid,0,0,sm);
}



/**
  The function computes mass matrix of given element.

  @param eid - element id.
  @param mm  - mass %matrix where the results are stored (output)

  @return The function returns assmebled %matrix in the parameter mm.

  Created by Tomas Koudelka, 12.8.2001
*/
void springel::mass_matrix (long /*eid*/,matrix &mm)
{
  fillm (0.0,mm);
}



/**
  The function computes strains on the required spring element for 
  the given load case.

  @param eid - element id
  @param lcid - load case id 
 
  @return The function does not return anything.

  Created by Tomas Koudelka, 12.8.2001
*/
void springel::strains (long eid,long lcid)
{
  long n = Mt->give_ndofe(eid);
  vector  r(n);
  vector  sig(n);
  vector  eps(1);
  long ii;

  eldispl (0,eid,r.a);
  ii = Mt->elements[eid].ipp[0][0];
  switch (Mt->elements[eid].te)
  {
    case spring_1:
      eps[0] = r[0];
      break;
    case spring_2:
      eps[0] = r[1];
      break;
    case spring_3:
      eps[0] = r[2];
      break;
    case spring_4:
      eps[0] = r[3];
      break;
    case spring_5:
      eps[0] = r[4];
      break;
    case spring_6:
      eps[0] = r[5];
      break;
    default:{break;}
  }
  Mm->storestrain (lcid,ii,eps);
}



/**
  The function computes stresses.

  @param eid - element id
  @param lcid - load case id 
 
  @return The function does not return anything.

  Created by Tomas Koudelka, 12.8.2001
*/
void springel::stresses (long eid,long lcid)
{
  long n = Mt->give_ndofe(eid);
  ivector cn(n);
  vector  r(n);
  vector  sig(1);
  vector  eps(1);
  matrix  d(1,1);
  long ii;

  ii = Mt->elements[eid].ipp[0][0];
  Mm->givestrain(lcid,ii,eps);
  Mm->matstiff(d, ii);
  mxv(d, eps, sig);
  Mm->storestress(lcid, ii, sig);
}



/**
  The function computes internal forces of given block.

  @param lcid - load case id
  @param eid - element id
  @param ri - block row id
  @param ci - block column id
  @param ifor - vector of internal forces

  @return The function returns %vector of internal forces in the parameter ifor

  Created by Tomas Koudelka, 12.8.2001
*/
void springel::internal_forces (long lcid,long eid,long ri,long ci,vector &ifor)
{
  long ii;
  matrix d(1,1);
  vector sig(1);

  ii = Mt->elements[eid].ipp[ri+0][ci+0];
  if (Mp->strcomp==1)
    Mm->computenlstresses (ii,Mm->ip[ii]);
  Mm->givestress (lcid,ii, sig);
  fillv(0.0, ifor);
  switch (Mt->elements[eid].te)
  {
    case spring_1:
      ifor[0] = sig[0];
      break;
    case spring_2:
      ifor[1] = sig[0];
      break;
    case spring_3:
      ifor[2] = sig[0];
      break;
    case spring_4:
      ifor[3] = sig[0];
      break;
    case spring_5:
      ifor[4] = sig[0];
      break;
    case spring_6:
      ifor[5] = sig[0];
      break;
    default:{break;}
  }
/*  fillv(0.0, ifor);
  switch (Mt->elements[eid].te)
  {
    case spring_1:
      ifor[0] = d[0][0]*r[0];
      break;
    case spring_2:
      ifor[1] = d[0][0]*r[1];
      break;
    case spring_3:
      ifor[2] = d[0][0]*r[2];
      break;
    case spring_4:
      ifor[3] = d[0][0]*r[3];
      break;
    case spring_5:
      ifor[4] = d[0][0]*r[4];
      break;
    case spring_6:
      ifor[5] = d[0][0]*r[5];
      break;
    default:{break;}
  }*/
}



/**
  The function computes resulting internal forces.

  @param lcid - load case id
  @param eid - element id
  @param ifor - %vector of internal forces (output)

  @return The function returns %vector of internal forces in the parameter ifor

  Created by Tomas Koudelka, 12.8.2001
*/
void springel::res_internal_forces (long lcid,long eid,vector &ifor)
{
  internal_forces (lcid,eid,0,0,ifor);
}



/**
  The function computes initial values of the given quantities at each integration point of the
  element from the nodal values given by the parameter nodval. Initial condition types must be 
  the same for all nodes of the element.

  @param eid - element id
  @param ri  - block row index
  @param ci  - block column index
  @param nodval - nodal values of particular initial conditions.
                  nodval[i][j] represents value of j-th initial condition at i-th node of the given element.
  @param ictn - array of types of initial condition for each node of element.
                The type of initial condition determines which values are being specified in the node. 
                (ictn[i] & inistrain) returns nonzero if nodal values of initial strains are specified
                (ictn[i] & inistress) returns nonzero if nodal values of initial stresses are specified
                (ictn[i] & iniother)  returns nonzero if nodal values of initial values of eqother array are specified
                (ictn[i] & inicond)   returns nonzero if nodal values of other initial conditions are specified

  @retur The function does not return anything.

  Created by Tomas Koudelka 0.5.2023
*/
void springel::inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn)
{
  long i, j, k, ipp;
  long ii, jj, nv = nodval.n;
  double ipval;
  vector anv(ASTCKVEC(nne));
  long nstra, nstre, ncompstr, ncompeqother;
  long idstra, idstre, idoth, idic;
  inictype ict;
  int aux;

  nstra = idstra = nstre = idstre = idoth = idic = 0;

  ict = ictn[0];
  for (i=0; i<nne; i++)
  {
    aux = int(ictn[i])-int(ict);
    if (aux < 0)  aux = -aux;    
    aux &= ~(inidisp);
    aux &= ~(inidisp_x);
    aux &= ~(inidisp_y);
    aux &= ~(inidisp_z);
    if ((ictn[i] != ict) && aux)
    {
      print_err("Incompatible types of initial conditions on element %ld\n"
                " at %ld. and %ld. nodes", __FILE__, __LINE__, __func__, eid+1, 1, i+1);
      abort();
    }
  }
  for (j = 0; j < nv; j++) // for all initial values
  {
    for(i = 0; i < nne; i++) // for all nodes on element
      anv[i] = nodval[i][j];
    for (ii = 0; ii < nb; ii++)
    {
      for (jj = 0; jj < nb; jj++)
      {
        ipp=Mt->elements[eid].ipp[ri+ii][ci+jj];
        if (intordsm[ii][jj] == 0)
          continue;
        for (k = 0; k < intordsm[ii][jj]; k++)
        {
          //  value in integration point
          ipval = anv(0);
          ncompstr =  Mm->ip[ipp].ncompstr;
          ncompeqother = Mm->ip[ipp].ncompeqother;
          if ((ictn[0] & inistrain) && (j < ncompstr))
          {
            Mm->ip[ipp].strain[idstra] += ipval;
            ipp++;
            continue;
          }
          if ((ictn[0] & inistress) && (j < nstra + ncompstr))
          {
            Mm->ip[ipp].stress[idstre] += ipval;
            ipp++;
            continue;
          }
          if ((ictn[0] & iniother) && (j < nstra+nstre+ncompeqother))
          {
            Mm->ip[ipp].eqother[idoth] += ipval;
            ipp++;
            continue;
          }
          if ((ictn[0] & inicond) && (j < nv))
          {
            if (Mm->ic[ipp] == NULL)
            {
              Mm->ic[ipp] = new double[nv-j];
              memset(Mm->ic[ipp], 0, sizeof(*Mm->ic[ipp])*(nv-j)); 
	    }
            Mm->ic[ipp][idic] += ipval;
            ipp++;
            continue;
          }
          ipp++;
        }
      }
    }
    ipp=Mt->elements[eid].ipp[ri][ci];
    ncompstr =  Mm->ip[ipp].ncompstr;
    ncompeqother = Mm->ip[ipp].ncompeqother;
    if ((ictn[0] & inistrain) && (j < ncompstr))
    {
      nstra++;
      idstra++;
      continue;
    }
    if ((ictn[0] & inistress) && (j < nstra + ncompstr))
    {      
      nstre++;
      idstre++;
      continue;
    }  
    if ((ictn[0] & iniother)  && (j < nstra + nstre + ncompeqother))
    {
      idoth++;
      continue;
    }
    if ((ictn[0] & inicond) && (j < nv))
    {
      idic++;
      continue;
    }
  }
}



/**
  Function interpolates the nodal values to the integration points on the element.
   
  @param eid - element id
  @param nodval - array of nodal values
  @param ipval - array of values at integration points
   
  @return The function returns %vector of interpolated values in the parameter ipval

  Created by Tomas Koudelka, 12.8.2001
*/
void springel::intpointval (long /*eid*/,vector &nodval,vector &ipval)
{
  copyv(nodval, ipval);
}


void springel::ipcoord(long eid, long /*ipp*/, vector &ipcoord)
{
  long nne = Mt->give_nne(eid);
  ivector enod(ASTCKIVEC(nne));
  Mt->give_nodal_coord(enod[0], ipcoord);
}
