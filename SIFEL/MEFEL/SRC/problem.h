#ifndef PROBLEM_H
#define PROBLEM_H

class probdesc;
class gtopology;
class mechtop;
class mechmat;
class mechcrsec;
class mechbclc;
class lhsrhs;
class gmatrix;


class problem
{
 public:
  problem (void);
  ~problem (void);
  
  void globinic (void);
  void dealoc (void);
  void deinic (void);
  
  
  probdesc *mp;
  gtopology *gt;
  
  mechtop *mt;
  mechmat *mm;
  mechcrsec *mc;
  mechbclc *mb;
  
  lhsrhs *lsrs;
  gmatrix *smat;
  
};

#endif
