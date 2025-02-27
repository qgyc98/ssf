#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hangnodet.h"



/**
  This constructor inializes attributes to zero values.

  Created by TKo, 07.2018
*/
hangnodet::hangnodet()
{
  nmn = 0;
  mnodes = NULL;
  natcoord = NULL;
  masentity = noel;
}



/**
  This destructor releases allocated dynamic memory

  Created by TKo, 07.2018
*/
hangnodet::~hangnodet()
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

  Created by TKo, 07.2018
*/
long hangnodet::read(XFILE *in)
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
  for (i=0; i<nmn; i++)
    xfscanf(in, "%ld", mnodes+i);
  xfscanf(in, "%le %le %le", natcoord, natcoord+1, natcoord+2);
  xfscanf(in, "%m", &gtypel_kwdset, &masentity);
  return 0;
}



/**
  The function prints additional data about hanging node
  to the opened text file out.

  @param out - pointer to the opened text file

  @retval 0 - on success

  Created by TKo, 07.2018
*/
long hangnodet::print(FILE *out)
{
  long i;
  
  fprintf(out, " %ld", -nmn);
  for (i=0; i<nmn; i++)
    fprintf(out, " %ld", mnodes[i]);
  fprintf(out, " %le %le %le", natcoord[0], natcoord[1], natcoord[2]);
  fprintf(out, " %d", masentity);
  return 0;
}
