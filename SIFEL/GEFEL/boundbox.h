#ifndef BOUNDBOX_H
#define BOUNDBOX_H

struct vector;
class  siftop;

/**
  Class for handling of a bounding box.
  
  Created by TKo, 12.2023
*/
class boundbox
{
 public:
  boundbox();
  boundbox(const siftop &top);
  
  /// initializes limit bound. box coordinates to the given values
  void init(double xcmin, double xcmax, double ycmin, double ycmax, double zcmin, double zcmax);
  /// initializes limit bound. box coordinates to the given values
  void init(const siftop &top);
  /// test if the given point is inside of the bounding box
  bool is_inside(const vector &p);
  /// test if the given point is inside of the bounding box
  bool is_inside(double x, double y, double z);

  /// minimum x coordinate of bbox
  double xmin;
  /// maximum x coordinate of bbox
  double xmax;
  /// minimum y coordinate of bbox
  double ymin;
  /// maximum y coordinate of bbox
  double ymax;
  /// minimum z coordinate of bbox
  double zmin;
  /// maximum z coordinate of bbox
  double zmax;    
};
#endif
