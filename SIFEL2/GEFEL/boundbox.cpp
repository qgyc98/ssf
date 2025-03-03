#include "boundbox.h"
#include "iotools.h"
#include "vector.h"
#include "siftop.h"
#include <stdlib.h>
#include <limits>



/**
  Constructor initializes bounding box limit coordinates.

  Created by TKo, 12.2023
*/
boundbox::boundbox()
{
  xmin = xmax = ymin = ymax = zmin = zmax = 0.0;
}



/**
  The constructor initializes bounding box limit coordinates from the givne topology.

  Created by TKo, 03.2024
*/
boundbox::boundbox(const siftop &top)
{
  init(top);
}



/**
  The function initializes bounding box limit coordinates from the givne topology.

  @param[in] top - topology which is used for the bounding box initialization.

  @return The function does not return antyhthing but changes coordinate limits of 
          the given bounding box instance.

  Created by TKo, 03.2024
*/
void boundbox::init(const siftop &top)
{
  long i;
  xmin = ymin = zmin =  std::numeric_limits<double>::max();
  xmax = ymax = zmax = -std::numeric_limits<double>::max();
  for (i=0; i<top.nn; i++){
    if (top.nodes[i].x < xmin)
      xmin = top.nodes[i].x;
    if (top.nodes[i].y < ymin)
      ymin = top.nodes[i].y;
    if (top.nodes[i].z < zmin)
      zmin = top.nodes[i].z;
    if (top.nodes[i].x > xmax)
      xmax = top.nodes[i].x;
    if (top.nodes[i].y > ymax)
      ymax = top.nodes[i].y;
    if (top.nodes[i].z > zmax)
      zmax = top.nodes[i].z;
  }
}


/**
  Initializes bounding box limit coordinates.

  @param[in] xcmin - minimum x-coordinate
  @param[in] xcmax - maximum x-coordinate
  @param[in] ycmin - minimum y-coordinate
  @param[in] ycmax - maximum y-coordinate
  @param[in] zcmin - minimum z-coordinate
  @param[in] zcmax - maximum z-coordinate

  Created by TKo, 12.2023
*/
void boundbox::init(double xcmin, double xcmax, double ycmin, double ycmax, double zcmin, double zcmax)
{
  if ((xcmax < xcmin) || (ycmax < ycmin) || (zcmax < zcmin)){
    print_err("invalid limits of the bounding box"
              "(xcmax < xcmin)=%d, (ycmax < ycmin)=%d, (zcmax < zcmin)=%d", __FILE__, __LINE__, __func__,
              (xcmax < xcmin), (ycmax < ycmin), (zcmax < zcmin));
    abort();
  }
    
  xmin = xcmin;  xmax = xcmax;
  ymin = ycmin;  ymax = ycmax;
  zmin = zcmin;  zmax = zcmax;
}



/**
  The function returns whether the given point is inside of the bounding box.

  @param[in] p - point coordinates collected in three-component %vector.

  @retval true  : if the point is inside the given bounding box.
  @retval false : if the point is outside the given bounding box.

  Created by TKo, 12.2023.
*/
bool boundbox::is_inside(const vector &p)
{
  if ((p[0] > xmax) || (p[0] < xmin) ||
      (p[1] > ymax) || (p[1] < ymin) ||
      (p[2] > zmax) || (p[2] < zmin))
    return false;
                        
  return true;
}



/**
  The function returns whether the given point is inside of the bounding box.

  @param[in] x - x-coordinate of the investigated point.
  @param[in] y - y-coordinate of the investigated point.
  @param[in] z - z-coordinate of the investigated point.

  @retval true  : if the point is inside the given bounding box.
  @retval false : if the point is outside the given bounding box.

  Created by TKo, 12.2023.
*/
bool boundbox::is_inside(double x, double y, double z)
{
  if ((x > xmax) || (x < xmin) ||
      (y > ymax) || (y < ymin) ||
      (z > zmax) || (z < zmin))
    return false;
                        
  return true;
}
