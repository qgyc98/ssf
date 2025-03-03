#include "import_T3d.h"
#include "siftop_element_types.h"
#include "iotools.h"
#include <stdlib.h>
#include <stdio.h>



siftop *S=NULL; ///< Pointer to siftop structure, which will contain topology from t3d file
XFILE  *f;      ///< Pointer to the input file



long type,degree,renum,mode ;
long nodes,edges,trias,quads,tetras,pyrams,wedges,bricks,elements ;
long modeI[6] ;
long nodeI[20] ;

#define Vertex 1
#define Curve 2
#define Surface 3
#define Region 4
#define Patch 5
#define Shell 6

#define Parameters 0
#define Tangent 1
#define Normal 2
#define BoundaryEntities 3
#define AssociatedElements 4
#define Neighbourhood 5



/**
  Function sets topology pointer of global variable S to value oS.

   @param oS structure of siftop

   @retval 0 always

   created  29.5.2000, Ondrej Hrstka
   modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long import_T3d_set_topology_pointer ( siftop &oS )
{
    S=&oS ;
    return 0 ;
}



/**
  This function decodes mode variable and stores results in modeI array

   created  29.5.2000, Ondrej Hrstka
   modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
void decode_mode ( void )
{
  long i ;
  for ( i=0 ; i<6 ; i++ )
    modeI[i]=( mode>>i ) & 1 ;
}



/**
  This function imports node coordinates and property number from t3d file f.
  Results are stored in S variable.

  @retval 0 always

   created  29.5.2000, Ondrej Hrstka
   modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long import_T3d_nodes ( long paral )
{
  long i,id,property,entity,eid ;
  double x,y,z ;

  if (paral)
    S->gnn = new long[nodes];
  for ( i=0 ; i<nodes ; i++ )
  {
    xfscanf(f, "%ld", &id);
    if (paral)
      xfscanf(f, "%ld", &S->gnn[id-1]);

    xfscanf(f, "%le", &x);
    xfscanf(f, "%le", &y);
    xfscanf(f, "%le", &z);
    xfscanf(f, "%ld", &entity);
    xfscanf(f, "%ld", &eid);
    xfscanf(f, "%ld", &property);

    S->nodes[id-1].alloc(1);
    S->nodes[id-1].x = x;
    S->nodes[id-1].y = y;
    S->nodes[id-1].z = z;
    S->nodes[id-1].entid[0] = evertex;
    S->nodes[id-1].prop[0]  = property;


    switch ( entity )
    {
      case Vertex:
      case Region:
        break ;
      case Curve:
        if ( modeI[Parameters] || modeI[Tangent] )
          skipline(f); // skip rest of line
        break ;
      case Surface:
      case Patch:
      case Shell:
        if ( modeI[Parameters] || modeI[Normal] )
          skipline(f); // skip rest of line
        break ;
    }
  }
  return(0);
}



/**
  This function imports 1D elements from t3d file f.
  Results are stored in S variable.

  @retval 0 always

   created  29.5.2000, Ondrej Hrstka
   modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long import_T3d_edges ( void )
{
  long id,esize,i,j,dummy,property;
  gtypel etype;
  
  switch ( degree )
  {
    case 1:
      etype=isolinear1d;
      esize=2;
      break ;
    case 2:
      etype=isoquadratic1d;
      esize=3;
      break ;
    default:
      print_err("unknown degree of element is required", __FILE__, __LINE__, __func__);
  }

  for ( i=0 ; i<edges ; i++ )
  {
    xfscanf(f, "%ld", &id);
    for ( j=0 ; j<esize ; j++ )
      xfscanf(f, "%ld", nodeI+j);

    xfscanf(f, "%ld", &dummy);
    xfscanf(f, "%ld", &dummy);
    xfscanf(f, "%ld", &property);

    S->elements[id-1].type = etype;
    S->elements[id-1].alloc(1);
    for(j=0; j<esize; j++)
      S->elements[id-1].nodes[j] = nodeI[j]; 
    S->elements[id-1].prop = property;
  }
  return(0);
}



/**
  This function imports 2D triangle elements from t3d file f.
  Results are stored in S variable.

  @retval 0 always

   created  29.5.2000, Ondrej Hrstka
   modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long import_T3d_trias ( void )
{
  long id,esize,i,j,dummy,property ;
  gtypel etype;
  
  switch (degree)
  {
    case 1:
      etype=trianglelinear;
      esize=3;
      break ;
    case 2:
      etype=trianglequadratic;
      esize=6;
      break ;
    default:
      print_err("unknown degree of element is required", __FILE__, __LINE__, __func__);
  }

  for ( i=0 ; i<trias ; i++ )
  {
    xfscanf(f, "%ld", &id);
    for ( j=0 ; j<esize ; j++ )
      xfscanf(f, "%ld", nodeI+j);

    xfscanf(f, "%ld", &dummy);
    xfscanf(f, "%ld", &dummy);
    xfscanf(f, "%ld", &property);

    S->elements[id-1].type = etype;
    S->elements[id-1].alloc(1);
    for(j=0; j<esize; j++)
      S->elements[id-1].nodes[j] = nodeI[j]; 
    S->elements[id-1].prop = S->elements[id-1].propsurf[0] = property;

    if ( modeI[Neighbourhood] )
    {
        
      skipline(f); // skip rest of line
    }
    if (modeI[BoundaryEntities])
    {
      for (j = 0; j < S->elements[id-1].ned; j++)
      {
        xfscanf(f, "%ld", &dummy);
      }
      for (j = 0; j < S->elements[id-1].ned; j++)
      {
        xfscanf(f, "%ld", &S->elements[id-1].propedg+j);
      }
    }
    else
      skipline(f); // skip rest of line
  }
  return(0);
}



/**
  This function imports 2D quad elements from t3d file f.
  Results are stored in S variable.

  @retval 0 always

   created  29.5.2000, Ondrej Hrstka
   modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long import_T3d_quads ( void )
{
  long id,esize,i,j,dummy,property ;
  gtypel etype;
  
  switch (degree)
  {
    case 1:
      etype=isolinear2d;
      esize=4;
      break ;
    case 2:
      etype=isoquadratic2d;
      esize=8;
      break ;
    default:
      print_err("unknown degree of element is required", __FILE__, __LINE__, __func__);
  }

  for ( i=0 ; i<quads ; i++ )
  {
    xfscanf(f, "%ld", &id);
    for ( j=0 ; j<esize ; j++ )
      xfscanf(f, "%ld", nodeI+j);

    xfscanf(f, "%ld", &dummy);
    xfscanf(f, "%ld", &dummy);
    xfscanf(f, "%ld", &property);

    S->elements[id-1].type = etype;
    S->elements[id-1].alloc(1);
    for(j=0; j<esize; j++)
      S->elements[id-1].nodes[j] = nodeI[j]; 
    S->elements[id-1].prop = S->elements[id-1].propsurf[0] = property;

    if ( modeI[Neighbourhood] )
    {    
      skipline(f); // skip rest of line
    }
    if (modeI[BoundaryEntities])
    {
      for (j = 0; j < S->elements[id-1].ned; j++)
      {
        xfscanf(f, "%ld", &dummy);
      }
      for (j = 0; j < S->elements[id-1].ned; j++)
      {
        xfscanf(f, "%ld", &S->elements[id-1].propedg+j);
      }
    }
    else
      skipline(f); // skip rest of line
  }
  return(0);
}



/**
  This function imports 3D tetrahedron elements from t3d file f.
  Results are stored in S variable.

  @retval 0 always

   created  29.5.2000, Ondrej Hrstka
   modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified 10.11.2005, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long import_T3d_tetras ( void )
{
  long id,esize,i,j,dummy,property ;
  gtypel etype;
  
  switch (degree)
  {
    case 1:
      etype=tetrahedronlinear;
      esize=4;
      break ;
    case 2:
      etype=tetrahedronquadratic;
      esize=10;
      break ;
    default:
      print_err("unknown degree of element is required", __FILE__, __LINE__, __func__);
  }

  for ( i=0 ; i<tetras ; i++ )
  {
    xfscanf(f, "%ld", &id);
    for ( j=0 ; j<esize ; j++ )
      xfscanf(f, "%ld", nodeI+j);

    xfscanf(f, "%ld", &dummy);
    xfscanf(f, "%ld", &dummy);
    xfscanf(f, "%ld", &property);

    S->elements[id-1].type = etype;
    S->elements[id-1].alloc(1);
    for(j=0; j<esize; j++)
      S->elements[id-1].nodes[j] = nodeI[j]; 
    S->elements[id-1].prop = property;

    if ( modeI[Neighbourhood] || modeI[AssociatedElements] )
    {
      skipline(f); // skip rest of line
    }
    if (modeI[BoundaryEntities])
    {
      for (j = 0; j < S->elements[id-1].nsurf; j++)
      {
        xfscanf(f, "%ld", &dummy);
        xfscanf(f, "%ld", &dummy);
      }
      for (j = 0; j < S->elements[id-1].nsurf; j++)
      {
        xfscanf(f, "%ld", &S->elements[id-1].propsurf+j);
      }
      skipline(f); // skip rest of line
    }
    else
      skipline(f); // skip rest of line
  }
  return(0);
}



/**
  This function imports 3D pyramidal elements from t3d file f.
  Results are stored in S variable.

  @retval 0 always

   created  29.5.2000, Ondrej Hrstka
   modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long import_T3d_pyrams ( void )
{
  long id,esize,i,j,dummy,property ;
  gtypel etype;
  
  switch (degree)
  {
    case 1:
      etype=pyramidelinear;
      esize=5;
      break ;
    case 2:
      etype=pyramidequadratic;
      esize=13;
      break ;
    default:
      print_err("unknown degree of element is required", __FILE__, __LINE__, __func__);
  }

  for ( i=0 ; i<pyrams ; i++ )
  {
    xfscanf(f, "%ld", &id);
    for ( j=0 ; j<esize ; j++ )
      xfscanf(f, "%ld", nodeI+j);

    xfscanf(f, "%ld", &dummy);
    xfscanf(f, "%ld", &dummy);
    xfscanf(f, "%ld", &property);

    S->elements[id-1].type = etype;
    S->elements[id-1].alloc(1);
    for(j=0; j<esize; j++)
      S->elements[id-1].nodes[j] = nodeI[j]; 
    S->elements[id-1].prop = property;

    if ( modeI[BoundaryEntities] || modeI[Neighbourhood] || modeI[AssociatedElements] )
    {
      skipline(f); // skip rest of line
    }
/*        if ( modeI[Neighbourhood] || modeI[AssociatedElements] )
        {
          skipline(f); // skip rest of line
        }
        if (modeI[BoundaryEntities])
        {
          for (j = 0; j < S->elements[id-1].nsurf; j++)
          {
            xfscanf(f, "%ld", &S->elements[id-1].propsurf+j);
          }
        }*/
  }
  return(0);
}



/**
  This function imports 3D wedge elements from t3d file f.
  Results are stored in S variable.

  @retval 0 always

   created  29.5.2000, Ondrej Hrstka
   modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long import_T3d_wedges ( void )
{
  long id,esize,i,j,dummy,property ;
  gtypel etype;
  
  switch (degree)
  {
    case 1:
      etype=wedgelinear;
      esize=6;
      break ;
    case 2:
      etype=wedgequadratic;
      esize=15;
      break ;
    default:
      print_err("unknown degree of element is required", __FILE__, __LINE__, __func__);
  }

  for ( i=0 ; i<wedges ; i++ )
  {
    xfscanf(f, "%ld", &id);
    for ( j=0 ; j<esize ; j++ )
      xfscanf(f, "%ld", nodeI+j);

    xfscanf(f, "%ld", &dummy);
    xfscanf(f, "%ld", &dummy);
    xfscanf(f, "%ld", &property);

    S->elements[id-1].type = etype;
    S->elements[id-1].alloc(1);
    for(j=0; j<esize; j++)
      S->elements[id-1].nodes[j] = nodeI[j]; 
    S->elements[id-1].prop = property;

    if ( modeI[BoundaryEntities] || modeI[Neighbourhood] || modeI[AssociatedElements] )
    {
      skipline(f); // skip rest of line
    }
/*        if ( modeI[Neighbourhood] || modeI[AssociatedElements] )
        {
          skipline(f); // skip rest of line
        }
        if (modeI[BoundaryEntities])
        {
          for (j = 0; j < S->elements[id-1].nsurf; j++)
          {
            xfscanf(f, "%ld", &S->elements[id-1].propsurf+j);
          }
        }*/
  }
  return(0);
}



/**
  This function imports 3D brick elements from t3d file f.
  Results are stored in S variable.

  @retval 0 always

   created  29.5.2000, Ondrej Hrstka
   modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz

*/
long import_T3d_bricks ( void )
{
  long id,esize,i,j,dummy,property ;
  gtypel etype;
  
  switch (degree)
  {
    case 1:
      etype=isolinear3d;
      esize=8;
      break ;
    case 2:
      etype=isoquadratic3d;
      esize=20;
      break ;
    default:
      print_err("unknown degree of element is required", __FILE__, __LINE__, __func__);
  }

  for ( i=0 ; i<bricks ; i++ )
  {
    xfscanf(f, "%ld", &id);
    for ( j=0 ; j<esize ; j++ )
      xfscanf(f, "%ld", nodeI+j);

    xfscanf(f, "%ld", &dummy);
    xfscanf(f, "%ld", &dummy);
    xfscanf(f, "%ld", &property);

    S->elements[id-1].type = etype;
    S->elements[id-1].alloc(1);
    for(j=0; j<esize; j++)
      S->elements[id-1].nodes[j] = nodeI[j]; 
    S->elements[id-1].prop = property;


    if ( modeI[Neighbourhood] || modeI[AssociatedElements] )
    {
      skipline(f); // skip rest of line
    }

    if (modeI[BoundaryEntities])
    {
      for (j = 0; j < S->elements[id-1].nsurf; j++)
      {
        xfscanf(f, "%ld", &dummy);
        xfscanf(f, "%ld", &dummy);
      }
      for (j = 0; j < S->elements[id-1].nsurf; j++)
      {
        xfscanf(f, "%ld", &S->elements[id-1].propsurf+j);
      }
    }
    else
      skipline(f); // skip rest of line
  }
  return(0);
}



/**
  This function imports 3D pyramidal elements from t3d file f.
  Results are stored in S variable.

  @param ofname - string with valid name and path of T3d output file

  @retval -1 if topology pointer S is not set
  @retval -2 if file specified by ofname cannot be opened
  @retval  1

   created  29.5.2000, Ondrej Hrstka
   modified 23.1.2002, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long import_T3d ( char *ofname, long paral )
{
  long temp;
  char emsg[1100];
  
  if ( !S )
  {
    print_err("topology pointer is not set",__FILE__, __LINE__, __func__ ) ;
    return -1 ;
  }
  f=xfopen( ofname,"rt" ) ;  
  if ( !f )
  {
    sprintf(emsg, "cannot open T3d output file %s for reading",ofname);
    print_err(emsg, __FILE__, __LINE__, __func__);
    return -2 ;
  }
  f->warning = 1;
  f->kwdmode = ignore_kwd;
  f->ignorecase = 1;

  xfscanf(f, "%ld", &type);   // reads type of generator
  xfscanf(f, "%ld", &degree); // reads degrees of element aproximation
  xfscanf(f, "%ld", &renum);  // reads indicator of renumbering
  xfscanf(f, "%ld", &mode);   // reads mode for elements
    
  if (paral==1){
    xfscanf(f, "%ld", &temp);
  }
    
  decode_mode() ;

  xfscanf(f, "%ld", &nodes);
//    nodes++ ;

  switch ( type )
  {
    case 3:
      xfscanf(f, "%ld", &edges);
      xfscanf(f, "%ld", &trias);
      xfscanf(f, "%ld", &tetras);
      quads=0 ;
      pyrams=0 ;
      wedges=0 ;
      bricks=0 ;
      break ;
    case 4:
      xfscanf(f, "%ld", &edges);
      xfscanf(f, "%ld", &quads);
      xfscanf(f, "%ld", &bricks);
      trias=0 ;
      tetras=0 ;
      pyrams=0 ;
      wedges=0 ;
      break ;
    case 7:
      xfscanf(f, "%ld", &edges);
      xfscanf(f, "%ld", &trias);
      xfscanf(f, "%ld", &quads);
      xfscanf(f, "%ld", &tetras);
      xfscanf(f, "%ld", &pyrams);
      xfscanf(f, "%ld", &wedges);
      xfscanf(f, "%ld", &bricks);
      break ;
  }

  S->nn = nodes;
  import_T3d_nodes(paral) ;

//    elements=edges+trias+quads+tetras+pyrams+wedges+bricks+1 ;
  elements=edges+trias+quads+tetras+pyrams+wedges+bricks ;
  S->ne = elements;

  import_T3d_edges();
  import_T3d_trias();
  import_T3d_quads();
  import_T3d_tetras();
  import_T3d_pyrams();
  import_T3d_wedges();
  import_T3d_bricks();

  return 1 ;
}










