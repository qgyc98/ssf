#ifndef RICHARDS_H
#define RICHARDS_H

#include "genfile.h"
#include "aliast.h"

class richards
{
  public:
    richards();    //constructor
    ~richards();   //destructor
    
    void read(XFILE *in);
    void print(FILE *out);
    
    ///  dimension of problem solved
    long dim;
    
    
    ///  function computes hydraulic conductivity
    double kk_value (double h);
    ///  function computes the derivative of K with respect to h (dK/dh van Genuchten).
    double dkkdh_value (double h);
    ///  function computes second derivative of K with respect to h (d^2K/dh^2 van Genuchten).
    double ddkkddh_value (double h);
    ///  function computes capacity (van Genuchten).
    double c_value (double h);
    ///  function computes derivative of capacity with respect to hydraulic head (van Genuchten).
    double dcdh_value (double h);
    ///  function computes water content (retention curve) (van Genuchten).
    double theta_val (double h);
    ///  function computes Darcian velocity (Dracy-Buckingham law)
    vector darcian_flux(long ipp);
    /// This function computes first order Taylor serie error of the solution over time. 
    double taylor_error(long ipp) ;
    
    void matcond (matrix &d,long ri,long ci,long ipp);
    void matcap (double &c,long ri,long ci,long ipp);
    
    void matcond1d (matrix &d,long ri,long ci,long ipp);
    void matcond2d (matrix &d,long ri,long ci,long ipp);
    void matcond3d (matrix &d,long ri,long ci,long ipp);

    void matcond2 (matrix &d,long ri,long ci,long ipp);
    void matcond1d_2 (matrix &d,long ri,long ci,long ipp);
    void matcond2d_2 (matrix &d,long ri,long ci,long ipp);
    void matcond3d_2 (matrix &d,long ri,long ci,long ipp);
    
    void give_dof_names(namevart *dofname, long ntm);
   

  private:
    ///  parameters of van Genuchten model
    double alpha;
    double n;
    double m;
    

    ///  residual water content
    double thetar;
    ///  specific storage
    double storage;
    
  public:
    ///  saturated hydraulic conductivities
    double kksxx,kksxy,kksxz,kksyy,kksyz,kkszz;
    ///  saturated water content
    double thetas;
};

#endif
