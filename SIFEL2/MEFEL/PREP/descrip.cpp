#include "descrip.h"
#include <stdio.h>


descrip::descrip()
{
  topf[0] = matf[0] = crf[0] = hangnf[0] = '\0';

  meshfmt = sifel;
  paral = 0;
  redgn = 0;
  matsec = crssec = no;
  matstr = crsstr = no;
  matkwd = crskwd = no;
}



descrip::~descrip()
{
}



long descrip::print(FILE *out)
{
  fprintf(out, "%s\n", topf);
  if (matsec == no)
    fprintf(out, "%s\n", matf);
  if (crssec == no)
    fprintf(out, "%s\n", crf);

  fprintf(out, "mesh_format %d  # mesh format indicator\n", int(meshfmt));
  fprintf(out, "edge_numbering %ld  # edge/surface property on elements indicator\n", redgn);

  if (hangnf[0])
    fprintf(out, "hanging_nodes_file %s", hangnf);

  if ((matstr == no) || (matkwd == yes))
  {
    fprintf(out, "read_mat_strings %d\n", matstr);
    fprintf(out, "read_mat_kwd %d\n", matkwd);
  }

  if ((crsstr == no) || (crskwd == yes))
  {
    fprintf(out, "read_crs_strings %d\n", crsstr);
    fprintf(out, "read_crs_kwd %d\n", crskwd);
  }

  return 0;
}
