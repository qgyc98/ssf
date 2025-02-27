#include "iotools.h"
#include "xfile.h"
#include "siftop.h"


int main(int argc, char*argv[])
{
  XFILE *in;
  FILE  *out;
  long i, j;
  long nn;
  long ne, eid;
  snode    *nodes;
  selement *elements;
  long nedges, nfaces;
  sedges   *edges = NULL;
  sfaces   *surfaces = NULL;
  long maxorignid;
  long     *new2orignid, *orig2newnid;

  fprintf (stdout,"\n\n *** NODE RENUMBERING AFTER DELETED ELEMENTS ***\n");
  fprintf (stdout," --------------------------------------------------\n");
  if (argc < 4){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : mshed-nodren mesh_input_file_name   mesh_output_file_name  read_elem_boundary_prop(0=no,1=yes,2=yes for 2D elements)\n\n");
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

  //  number of nodes
  xfscanf(in, "%k%ld", "num_nodes", &nn);
  nodes = new snode[nn];
  new2orignid = new long[nn];
  memset(new2orignid, 0, sizeof(*new2orignid)*nn);
  maxorignid = -1L;
  for (i = 0; i < nn; i++)
  {
    xfscanf(in, "%k%ld", "node_id", new2orignid+i);
    new2orignid[i]--;
    if (new2orignid[i] > maxorignid)
      maxorignid = new2orignid[i]; 
    xfscanf(in, "%k%le %k%le %k%le %k%ld", "x", &nodes[i].x, "y", &nodes[i].y, "z", &nodes[i].z, "numprop", &nodes[i].nprop);
    if (nodes[i].nprop > 0)
    {
      nodes[i].alloc(nodes[i].nprop);
      for (j=0; j < nodes[i].nprop; j++)
        xfscanf(in, "%k%m%ld", "prop", &entityp_kwdset, nodes[i].entid+j, nodes[i].prop+j);
    }
  }
  orig2newnid = new long[maxorignid+1];
  memset(orig2newnid, 0, sizeof(*orig2newnid)*(maxorignid+1));
  for (i=0; i<nn; i++){
    orig2newnid[new2orignid[i]] = i;
  }

  
  //  number of elements
  xfscanf(in, "%k%ld", "num_elements", &ne);
  elements = new selement[ne];
  for (i = 0; i < ne; i++)
  {
    xfscanf(in, "%k%ld", "elem_id", &eid);
    xfscanf(in, "%k%m", "eltype", &gtypel_kwdset, &elements[i].type);
    if (elements[i].alloc(rebp))
    {
      print_err("unknown type on element %ld is required (file %s, line=%ld, col=%ld)",
                __FILE__, __LINE__, __func__, eid+1, in->fname, in->line, in->col);
      return (2);
    }
    xfscanf(in, "%k", "enodes");
    for (j = 0; j < elements[i].nne; j++)
    {
      xfscanf(in, "%ld", &elements[i].nodes[j]);
      elements[i].nodes[j]--;
      elements[i].nodes[j] = orig2newnid[elements[i].nodes[j]];
      if ((elements[i].nodes[j] < 0) || (elements[i].nodes[j] >= nn))
      {
        print_err("node number %ld on element %ld is out of range <1, %ld> (file %s, line=%ld, col=%ld)",
                  __FILE__, __LINE__, __func__, elements[i].nodes[j]+1, i+1, nn, in->fname, in->line, in->col);
        return (2);
      }
    }
    xfscanf(in, "%k%ld", "eprop", &elements[i].prop);
    if (rebp > 0)
    {
      if (elements[i].ned)
        xfscanf(in, "%k", "propedg");
      for(j = 0; j < elements[i].ned; j++)
        xfscanf(in, "%ld", elements[i].propedg+j);
      if (rebp == 1) // due to compatibility with Mesheditor output which cannot print surface property for 2D elements
      {              // in the case of modified file by Mesheditor, readege flag can be set to 2 and zero surface properties will be used on all 2D elements
        if (elements[i].nsurf)
          xfscanf(in, "%k", "propsurf");
        for(j = 0; j < elements[i].nsurf; j++)
          xfscanf(in, "%ld", elements[i].propsurf+j);
      }
    }
    else{ // skip rest of line
      skipline(in);
    }
  }  

  
  fprintf(out, "%ld\n", nn);
  for(i=0; i<nn; i++){
    fprintf (out,"%ld %15.10le %15.10le %15.10le %ld",i+1, nodes[i].x, nodes[i].y, nodes[i].z, nodes[i].nprop);
    for(j=0; j<nodes[i].nprop; j++)
    {
      fprintf(out, " %d %ld", nodes[i].entid[j], nodes[i].prop[j]);
    }
    fprintf(out, "\n");
  }

  fprintf(out, "\n%ld\n", ne);
  for(i=0; i<ne; i++){
    fprintf (out,"%ld %d",i+1, elements[i].type);
    for (j = 0; j < elements[i].nne; j++)
      fprintf (out," %ld",elements[i].nodes[j]+1);
    fprintf (out," %ld",elements[i].prop);
    if (elements[i].propedg)
    {
      fprintf(out, "  ");
      for (j=0; j < elements[i].ned; j++)
        fprintf (out," %ld",elements[i].propedg[j]);
    }
    if (elements[i].propsurf)
    {
      fprintf(out, "  ");
      for (j=0; j < elements[i].nsurf; j++)
        fprintf (out," %ld",elements[i].propsurf[j]);
    }
    fprintf(out, "\n");
  }

  delete [] nodes;
  delete [] elements;

  kwd_handling bmode = in->kwdmode;
  in->kwdmode = sequent_mode;
  // read additional data about surfaces created in the MeshEditor
  if (xfscanf(in, " %+k", "faces"))
  {
    in->kwdmode = bmode;
    xfscanf(in, "%ld", &nfaces);  
    surfaces = new sfaces(nfaces);
    surfaces->read(in);
    fprintf(out, "faces %ld\n", nfaces);
    for (i=0; i<nfaces; i++){
      for (j=0; j<surfaces->faces[i].nnod; j++)
        surfaces->faces[i].nodes[j] = orig2newnid[surfaces->faces[i].nodes[j]];
      surfaces->faces[i].print(out, NULL);
    }    
  }
  delete surfaces;
  
  in->kwdmode = sequent_mode;
  // read additional data about edges created in the MeshEditor
  if (xfscanf(in, " %+k", "edges"))
  {
    in->kwdmode = bmode;
    xfscanf(in, "%ld", &nedges);
    edges = new sedges(nedges);
    edges->read(in); 
    fprintf(out, "edges %ld\n", nedges);
    for (i=0; i<nedges; i++){
      edges->edges[i].n1 = orig2newnid[edges->edges[i].n1];
      edges->edges[i].n2 = orig2newnid[edges->edges[i].n2];
      edges->edges[i].print(out, NULL);
    }
  }
  delete edges;

  delete [] orig2newnid;
  delete [] new2orignid;

  in->kwdmode = bmode;
  xfclose(in);
  fclose(out);
}
