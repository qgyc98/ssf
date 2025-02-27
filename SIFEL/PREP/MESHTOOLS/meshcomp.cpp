#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

#include "iotools.h"
#include "siftop.h"



/**
  Program compares two meshes and searches identical nodes and elements.
  Mesh files can be either in sifel format or t3d format.
  The format of mesh files and presence of global node numbering are specified interactively.
  Input of required maximum of distance tolerance is also specified interactively.

  Result of the program is a message at terminal window. One of the following results can 
  be obtained:
  Meshes are identical (but they can have different node and element numbering)
  One mesh is subdomain of the other (but they can have different node and element numbering)
  Meshes are different (all nodes of one mesh could not be paired with nodes of the other mesh or
                        connectivity of all elements of one mesh are not identical with
                        conncetivity of elements of the other mesh)

  Created by TKo, 5.2010
*/
int main (int argc,char *argv[])
{
  XFILE *in1;
  XFILE *in2;
  siftop top1, top2;
  meshform meshfmt1, meshfmt2;
  long paral1, paral2;
  double err;

  fprintf (stdout,"\n\n *** MESH COMPARATOR ***\n");
  fprintf (stdout," --------------------------\n");
  if (argc < 3){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : meshcomp mesh1_file_name mesh2_file_name\n\n");
    return(1);
  }
  in1 = xfopen(argv[1],"r");
  if (in1 == NULL){
    fprintf (stderr,"\n Cannot open mesh file 1.");
    return(2);
  }
  in2 = xfopen(argv[2],"r");
  if (in2 == NULL){
    fprintf (stderr,"\n Cannot open mesh file 2.");
    return(2);
  }

  printf("Input format of mesh 1 (0=sifel/1=t3d): ");
  scanf("%d", (int *)(&meshfmt1));
  printf("\nInput format of mesh 2 (0=sifel/1=t3d): ");
  scanf("%d", (int *)(&meshfmt2));
  printf("\nInput global node numbering indicator for mesh 1 (0=no/1=yes): ");
  scanf("%ld", &paral1);
  printf("\nInput global node numbering indicator for mesh 2 (0=no/1=yes): ");
  scanf("%ld", &paral2);
  printf("\nInput required maximum distance of two identical nodes (err): ");
  scanf("%lf", &err);
  printf("\n");

  printf("\n Reading of mesh 1 (%s) ...", argv[1]);
  fflush(stdout);
  switch (meshfmt1)
  {
    case t3d:
      top1.import_t3d(in1, paral1);
      break;
    case sifel:
      top1.read(in1, paral1, 0);
      break;
    default:
      print_err("unknown mesh format is required", __FILE__, __LINE__, __func__);
  }
  printf(" O.K.\n");
  xfclose(in1);
  printf("\n Reading of mesh 2 (%s) ...", argv[2]);
  fflush(stdout);
  switch (meshfmt2)
  {
    case t3d:
      top2.import_t3d(in2, paral2);
      break;
    case sifel:
      top2.read(in2, paral2, 0);
      break;
    default:
      print_err("unknown mesh format is required", __FILE__, __LINE__, __func__);
  }
  printf(" O.K.\n");
  xfclose(in2);
  fflush(stdout);
  
  printf("\n Comparing nodes ... ");


  double dx, dy, dz, dr;
  double x1, y1, z1;
  long *n1, *n2, sn1, sn2;
  long i, j, k, l;

  n1 = new long[top1.nn];
  memset(n1, 0, sizeof(*n1)*top1.nn);
  n2 = new long[top2.nn];
  memset(n2, 0, sizeof(*n2)*top2.nn);

  for(i=0; i<top1.nn; i++)
  {
    x1 = top1.nodes[i].x; 
    y1 = top1.nodes[i].y; 
    z1 = top1.nodes[i].z; 
    for(j=0; j<top2.nn; j++)
    {
      if (n2[j])
        continue;
      dx = x1 - top2.nodes[j].x;
      dy = y1 - top2.nodes[j].y;
      dz = z1 - top2.nodes[j].z;
      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      if (dr < err)
      {
        n1[i] = j+1;
        n2[j] = i+1;
        break;
      }
    }
    if ((i%10 == 0) || (i==top1.nn-1))
      printf("\r Comparing nodes ... %.1lf%%", double(i)*100/top1.nn);
  }
  printf("\n\n");
  for(i=0, sn1=0; i<top1.nn; i++)
  {
    if (n1[i])  sn1++;
  }
  for(i=0, sn2=0; i<top2.nn; i++)
  {
    if (n2[i])  sn2++;
  }

  printf("Mesh 1 (%s):\n total number of nodes: %ld\n number of paired nodes: %ld\n\n", argv[1], top1.nn, sn1);
  printf("Mesh 2 (%s):\n total number of nodes: %ld\n number of paired nodes: %ld\n\n", argv[2], top2.nn, sn2);
  printf("Result of comparison:\n");
  if ((sn1 == top1.nn) && (sn2 == top2.nn))
    printf(" Nodes are identical.\n");
  if ((sn1 == top1.nn) && (sn2 != top2.nn))
    printf(" Nodes of mesh1 are contained in mesh2.\n");
  if ((sn2 == top2.nn) && (sn1 != top1.nn))
    printf(" Nodes of mesh2 are contained in mesh1.\n");
  if ((sn2 != top2.nn) && (sn1 != top1.nn))
  {
    printf(" Meshes have different nodes.\n");
    printf("\nConclusion of comparison:\n");
    printf(" MESHES ARE DIFFERENT.\n\n");
    delete [] n1;
    delete [] n2;
    return 0;
  }

  printf("\n Comparing elements ... ");


  long nne1, fn, tfn;
  long *e1, *e2, se1, se2;

  e1 = new long[top1.ne];
  memset(e1, 0, sizeof(*e1)*top1.ne);
  e2 = new long[top2.ne];
  memset(e2, 0, sizeof(*e2)*top2.ne);
  for(i=0; i<top1.ne; i++)
  {
    nne1 = top1.elements[i].nne;

    // checking of paired nodes on element
    for (k=0; k<nne1; k++)
    {
      if (n1[top1.elements[i].nodes[k]] == 0) // node has not been paired
      {
        k = -1;
        break;
      }
    }
    if ((i%10 == 0) || (i==top1.ne-1))
      printf("\r Comparing elements ... %.1lf%%", double(i)*100/top1.ne);
    if (k < 0)  // one of element nodes has not been paired
      continue;

    for (j=0; j<top2.ne; j++)
    {
      if (e2[j])
        continue;
      if (nne1 == top2.elements[j].nne)
      {
        tfn = 0;
        for (k=0; k<nne1; k++)
        {
          fn = 0;
          for (l=0; l<nne1; l++)
          {
            if (n1[top1.elements[i].nodes[k]]-1 == top2.elements[j].nodes[l])
            {
              fn = 1;
              tfn++;
              break;  
            }
          }
          if (fn == 0)
            break;
        }
        if (tfn == nne1)
        {  
          e1[i] = j+1;
          e2[j] = i+1;
          break;
        }
      }
    }
  }
  printf("\n\n");

  for(i=0, se1=0; i<top1.ne; i++)
  {
    if (e1[i])  se1++;
  }
  for(i=0, se2=0; i<top2.ne; i++)
  {
    if (e2[i])  se2++;
  }

  printf("Mesh 1 (%s):\n total number of elements: %ld\n number of paired elements: %ld\n\n", argv[1], top1.ne, se1);
  printf("Mesh 2 (%s):\n total number of elements: %ld\n number of paired elements: %ld\n\n", argv[2], top2.ne, se2);
  printf("Conclusion of comparison:\n");
  if ((se1 == top1.ne) && (se2 == top2.ne))
    printf(" MESHES ARE IDENTICAL.\n\n");
  if ((se1 == top1.ne) && (se2 != top2.ne))
    printf(" MESH1 IS SUBDOMAIN OF MESH2.\n\n");
  if ((se2 == top2.ne) && (se1 != top1.ne))
    printf(" MESH2 IS SUBDOMAIN OF MESH1.\n\n");
  if ((se2 != top2.ne) && (se1 != top1.ne))
  {
    printf(" Meshes have different elements.\n\n");
    printf(" MESHES ARE DIFFERENT.\n\n");
  }
  delete [] n1;
  delete [] n2;
  delete [] e1;
  delete [] e2;
  return 0;
}
