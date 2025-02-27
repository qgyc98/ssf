#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#include "iotools.h"
#include "siftop.h"

siftop    *Top;     ///< Sifel topology

/**
  This structure holds data about input files of preprocessor
*/
struct descrip
{
  char topf[1025];   ///< topology file name
  meshform meshfmt;  ///< format of topology file
  long paral;        ///< indicator whether sequential version of preprocessor is used (= 0) or paralell (= 1)
  long redgn;        ///< indicator whether edge numbers of element should be read in topology file (=1)
};


/**
  Function reads topology file in, format of the file is described by structure d.

  @param in - pointer to the opened XFILE structure


  Function returns:
  @retval 0 : on success
  @retval 4 : in case of unknown mesh format

  In case of reading alternative file format :
  @retval 1 : on error in reading of node
  @retval 2 : on error in reading of element
  @retval 3 : on error in reading of global node numbers
*/
long input_siftop(XFILE *in,descrip *d)
{
  long ret;

  fprintf(stdout, "\n\nReading of mesh topology . . .");
  switch (d->meshfmt)
  {
    case t3d:
      //  it works only for one mesh, mesh decomposed into submeshes cannot be used
      Top->import_t3d(in, 0);
      break;
    case sifel:
      //  it works only for one mesh, mesh decomposed into submeshes cannot be used
      ret = Top->read(in, 0, d->redgn);
      return(ret);
    default:
      print_err("unknown mesh format is required", __FILE__, __LINE__, __func__);
      return(4);
  }
  return(0);
}


int main (int argc,char *argv[])
{
  long i;
  double xk,xp,r;
  double alpha,alphac,ys,xl,yl,x,y,xn;
  descrip d;
  FILE *out;
  //FILE *gid;
  XFILE *in,*itop;
  Top = new siftop;

  if (argc < 4){
    fprintf (stderr,"\n Wrong number of command line parameters.");
    fprintf (stderr,"\n Use : ohybak input_file_name output_file_name gid_file_name\n\n");
    return(1);
  }

  //  input file
  in = xfopen(argv[1],"r");
  if (in == NULL){
    fprintf (stderr,"\n Input file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }
  
  //  output file
  out = fopen(argv[2],"w");
  if (out == NULL){
    fprintf (stderr,"\n Output file has not been specified.");
    fprintf (stderr,"\n Try it again!\n\n");
    return(2);
  }
  
  // reading of line with topology file name
  xfscanf(in, " %1024s", d.topf);
  // reading of line with topology file format indicator,
  xfscanf(in, "%k%m", "mesh_format", &meshform_kwdset, &d.meshfmt);
  // reading of line with edge number indicator
  xfscanf(in, "%k%ld", "edge_numbering", &d.redgn);

  xfscanf (in,"%lf %lf %lf",&xp,&xk,&r);

  alphac=(xk-xp)/r;
  ys=r*(cos(alphac-1.0));
      
  //  input file
  itop = xfopen(d.topf,"r");
  if (in == NULL){
    fprintf (stderr,"\n Topology file has not been specified.");
    return(2);
  }
  input_siftop(itop,&d);
  
  double ox,oy;
  for (i=0;i<Top->nn;i++){
    x=Top->nodes[i].x;
    y=Top->nodes[i].y;
    ox=x;
    oy=y;

    //  prvni rovna cast
    if (ox<=xp)
      continue;
    
    //  smerovy oblouk, pro laiky zatacka
    if (xp<ox && ox<=xk){
      alpha=(x-xp)/r;
      ys=r*(cos(alpha)-1.0);
      x=xp+(r+y)*sin(alpha);
      y=ys+y*cos(alpha);

      Top->nodes[i].x=x;
      Top->nodes[i].y=y;
    }
    
    //  treti rovna cast, pro odborniky smerova prima
    if (xk<ox){
      xl=xp+(r+y)*sin(alphac);
      yl=ys+y*cos(alphac);
      xn=xl+(x-xk)*cos(alphac);
      y=yl-(x-xk)*sin(alphac);
      Top->nodes[i].x=xn;
      Top->nodes[i].y=y;
    }
  }

  Top->print (out);
  // gid = fopen(argv[3],"w");
//   Top->export_gid_mesh(gid);
//   fclose (gid);
//   xfclose (in);
//   xfclose (itop);
//   fclose (out);
//   fprintf (stdout,"\n");
}
