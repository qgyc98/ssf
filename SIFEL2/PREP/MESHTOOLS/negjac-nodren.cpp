#include "iotools.h"
#include "xfile.h"
#include "siftop.h"


int main(int argc, char*argv[])
{
  XFILE *in;
  FILE  *out;
  siftop top;
  long i, aux;
  

  fprintf (stdout,"\n\n *** NODE RENUMBERING FOR ELEMENTS WITH NEGATIVE JACOBIAN ***\n");
  fprintf (stdout," --------------------------------------------------------------------\n");
  if (argc < 4){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : negjac-nodren mesh_input_file_name   mesh_output_file_name  read_elem_boundary_prop(0=no,1=yes)\n\n");
    return(1);
  }
  in = xfopen(argv[1],"r");
  if (in == NULL){
    fprintf (stderr,"\n Cannot open input mesh file '%s'.\n", argv[1]);
    return(2);
  }
  
  out = fopen(argv[2],"w");
  if (out == NULL){
    fprintf (stderr,"\n Cannot open output mesh file '%s'.\n", argv[2]);
    return(2);
  }
  long rebp;
  sscanf(argv[3], "%ld", &rebp);
  if ((rebp < 0) || (rebp > 1)){
    fprintf(stderr, "\n Wrong indicator for the reading of element boundary properties, rebp=%ld but it must be 0 or 1\n", rebp);
    return (3);
  }

  top.read(in, 0, rebp);
  for(i=0; i<top.ne; i++){
    if (top.elements[i].type != isolinear3d)
      continue;
    // swap nodes 1 and 4
    aux = top.elements[i].nodes[0];
    top.elements[i].nodes[0] = top.elements[i].nodes[3];
    top.elements[i].nodes[3] = aux;
    // swap nodes 2 and 3
    aux = top.elements[i].nodes[1];
    top.elements[i].nodes[1] = top.elements[i].nodes[2];
    top.elements[i].nodes[2] = aux;
    // swap nodes 5 and 8
    aux = top.elements[i].nodes[4];
    top.elements[i].nodes[4] = top.elements[i].nodes[7];
    top.elements[i].nodes[7] = aux;
    // swap nodes 6 and 7
    aux = top.elements[i].nodes[5];
    top.elements[i].nodes[5] = top.elements[i].nodes[6];
    top.elements[i].nodes[6] = aux;
    
    // swap surface properties on surfaces 2 and 4
    aux = top.elements[i].propsurf[1];
    top.elements[i].propsurf[1] = top.elements[i].propsurf[3];
    top.elements[i].propsurf[3] = aux;
    
    // swap edge properties on edges 1 and 3
    aux = top.elements[i].propedg[0];
    top.elements[i].propedg[0] = top.elements[i].propedg[2];
    top.elements[i].propedg[2] = aux;
    // swap edge properties on edges 6 and 7
    aux = top.elements[i].propedg[5];
    top.elements[i].propedg[5] = top.elements[i].propedg[6];
    top.elements[i].propedg[6] = aux;
    // swap edge properties on edges 5 and 8
    aux = top.elements[i].propedg[4];
    top.elements[i].propedg[4] = top.elements[i].propedg[7];
    top.elements[i].propedg[7] = aux;
    // swap edge properties on edges 9 and 11
    aux = top.elements[i].propedg[8];
    top.elements[i].propedg[8] = top.elements[i].propedg[10];
    top.elements[i].propedg[10] = aux;    
  }

  top.print(out);
  
  xfclose(in);
  fclose(out);
}
