#ifndef CEBFIPCONTACTMAT_H
#define CEBFIPCONTACTMAT_H

#include "xfile.h"
struct matrix;

class cebfipcontactmat
{
  public:
    cebfipcontactmat (void);
    ~cebfipcontactmat (void);
    void read (XFILE *in);

    void matstiff (matrix &d,long ipp);

    void nlstresses (long ipp, long im, long ido);

    void updateval (long ipp, long im, long ido);

    long sgn (double a);


    /// concrete compressive strength
    double fcc;
    
    /// concrete tensile strength
    double fct;
    
    /// confining conditions /// 1 - unconfined, 2 - confined
    long conf;
    
    /// bond conditions      /// 1 - good, 2 - poor
    long bond;
    
    /// normal stiffnes
    double normal;
    
    /// length differention on graph
    double s;
    double ss;
    double sss;

    /// shape of the first branch
    double alfa;

    /// maximum shear stress
    double taumax;

    /// residual strength
    double tauf;
};

#endif
