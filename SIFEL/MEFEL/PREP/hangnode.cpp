#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hangnode.h"
#include "hngen.h"
#include "vector.h"



/**
  This constructor inializes attributes to zero values.

  Created by TKo, 10.2012
*/
hangnode::hangnode()
{
  nmn = 0;
  mnodes = NULL;
  natcoord = NULL;
  maset = noel;
}



/**
  This constructor intializes attributes from the array of master node numbers 
  and natural coordinates given by vectors mn a nd natc.
*/
hangnode::hangnode(ivector &mn, vector &natc, gtypel met)
{
  nmn = mn.n;
  mnodes = new long[nmn];
  copyv(mn, mnodes);
  natcoord = new double[3];
  copyv(natc, natcoord);
  maset = met;
}



/**
  This destructor releases allocated dynamic memory

  Created by TKo, 10.2012
*/
hangnode::~hangnode()
{
  delete [] mnodes;
  delete [] natcoord;
}



/**
  The function reads additional data about hanging node
  from the opened text file in.

  @param in - pointer to the opned text file

  @retval 0 - on success
  @retval 1 - wrong number of master nodes

  Created by TKo, 10.2012
*/
long hangnode::read(XFILE *in)
{
  long i;
  
  nmn = 0;
  xfscanf(in, "%ld", &nmn);
  nmn = -nmn;  // number of master nodes is in negative format
  if (nmn < 1)
    return 1;
  mnodes = new long[nmn];
  memset(mnodes, 0, sizeof(*mnodes)*nmn);
  natcoord = new double[3];
  for (i=0; i<nmn; i++){
    xfscanf(in, "%ld", mnodes+i);
    mnodes[i]--;
  }
  xfscanf(in, "%le %le %le", natcoord, natcoord+1, natcoord+2);
  xfscanf(in, "%m", &gtypel_kwdset, &maset);
  return 0;
}



/**
  The function prints additional data about hanging node
  to the opened text file out.

  @param out - pointer to the opened text file

  @retval 0 - on success

  Created by TKo, 10.2012
*/
long hangnode::print(FILE *out)
{
  long i;
  
  fprintf(out, " %ld", -nmn);
  for (i=0; i<nmn; i++)
    fprintf(out, " %ld", mnodes[i]+1);
  fprintf(out, " %le %le %le", natcoord[0], natcoord[1], natcoord[2]);
  fprintf(out, " %d", maset);
  return 0;
}



/**
  The function compares record of generated hanging nodes with the 
  given hanging node record for equality.

  @param hng[in] - record of generated hanging node to be compared

  @retval 0 - if the hanging node records are different
  @retval 1 - if the hanging node records are equal

  Created by TKo, 30.4.2023
*/
long hangnode::compare(const hngen &hng)
{
  if (nmn != hng.mnodes.n)
    return 1;
  for (long i=0; i<nmn; i++){
    if (mnodes[i] != hng.mnodes[i])
      return 1;
  }
  if (natcoord[0] != hng.xi)
    return 1;
  if (natcoord[1] != hng.eta)
    return 1;
  if (natcoord[1] != hng.zeta)
    return 1;
  if (maset != hng.et)
    return 1;

  return 0;
}

