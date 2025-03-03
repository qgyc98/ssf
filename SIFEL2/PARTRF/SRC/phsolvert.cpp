#include "phsolvert.h"
#include "pglobalt.h"
#include "seqfilest.h"
#include <string.h>
#include <math.h>
#include "mpi.h"

void par_homogenization ()
{
  long i,j,n;
  double s,zero,alpha,*d,*p,*lhs,*lhsi,*tdlhs,*rhs;

  double dt,end_time;
  

  n=Ndoft;
  
  lhs = Lsrst->give_lhs (0);
  tdlhs = Lsrst->give_tdlhs (0);
  lhsi = Lsrst->lhsi;
  rhs = Lsrst->give_rhs (0);

  d = new double [n];
  p = new double [n];
  //tdlhs = new double [n];

  //  initial values
  nullv (lhs,n);
  nullv (tdlhs,n);
  nullv (d,n);
  nullv (p,n);
  


  
  alpha=Tp->alpha;

  zero=Tp->zero;
  



  // **************************
  //  main iteration loop  ***
  // **************************

  //  loop counter
  i=0;
  
  //  starting time
  Tp->time = Tp->timecont.starttime ();
  //  time increment
  dt = Tp->timecont.initialtimeincr ();
  //  end time
  end_time = Tp->timecont.endtime ();
  

  //  nodes - integration points interpolation
  approximation ();
  actual_previous_change ();

  Tm->initmaterialmodels();
  

  print_initt(-1, "wt",Ptp->fni,Ptp->fei);
  print_stept(0,i,Tp->time,NULL);
  print_flusht();


  do{
    

    if (Myrank==0) fprintf (stdout,"\n iteration number %ld",i);
    
    //fprintf (Outt,"\n\n\n\n\n iteration number %ld",i);
    
    //par_aux_nonstat_print (gr,lhsi,d,time);
    
    
    //fprintf (Out,"\n\n\n\n\n\n kontrola integracnich bodu");
    //for (long ijk=0;ijk<Tm->nip;ijk++){
    //fprintf (Out,"\n %ld  %lf %lf",ijk,Tm->ip[ijk].av[0],Tm->ip[ijk].av[1]);
    //}


    //  zde bude rozeslani hodnot z mastera na slaves

    //  na slaves se resi jednotlive ulohy
    if (Myrank!=0){
      //  slaves solve problems on microscale
      // homog_trans ();
    }
    //  sber aktualnich dat ze slaves na mastera
    
    conductivity_matrix (0);

    capacity_matrix (0);
    
    if (i==0){
      Psolt->computation_stat (Gtt,Kmat);
    }
    
    //Kmat->printdiag (Outt);
    //fprintf (Outt,"\n\n MATICE K \n");
    //Kmat->printmat (Outt);
    //fprintf (Outt,"\n\n MATICE C \n");
    //Cmat->printmat (Outt);

    for (j=0;j<n;j++){
      p[j]=lhs[j]+(1.0-alpha)*dt*tdlhs[j];
    }
    
    //  auxiliary vector
    Kmat->gmxv (p,d);
    
    

    //fprintf (Outt,"\n\n\n\n\n kontrola K.(d+dt*(1-a)*v)");
    //for (long ijk=0;ijk<n;ijk++){
    //fprintf (Outt,"\n %ld %le",ijk,d[ijk]);
    //}
    

    
    Kmat->scalgm (alpha*dt);
    Kmat->addgm (1.0,*Cmat);
    
    
    
    //fprintf (Out,"\n\n\n\n\n kontrola d+dt*(1-a)*v");
    //for (long ijk=0;ijk<n;ijk++){
    //fprintf (Out,"\n %ld %lf",ijk,p[ijk]);
    //}

    
    
    
    //fprintf (Out,"\n\n\n\n\n kontrola K.(d+dt*(1-a)*v)");
    //for (long ijk=0;ijk<n;ijk++){
    //fprintf (Out,"\n %ld %lf",ijk,r[ijk]);
    //}

    
    //  right hand side
    trfel_right_hand_side (0,rhs,n);
    
    //fprintf (Outt,"\n\n\n\n\n kontrola prave strany");
    //for (long ijk=0;ijk<n;ijk++){
    //fprintf (Outt,"\n %ld %le",ijk,rhs[ijk]);
    //}
    
    
    
    for (j=0;j<n;j++){
      rhs[j] = rhs[j] - d[j];
      d[j]=tdlhs[j];
    }
    
    
    
    //fprintf (Outt,"\n\n\n\n\n kontrola prave strany");
    //for (long ijk=0;ijk<n;ijk++){
    //fprintf (Outt,"\n %ld %le",ijk,rhs[ijk]);
    //}

 
    //fprintf (Outt,"\n\n MATICE SOUSTAVY \n");
    //Kmat->printmat (Outt);




    //par_linear_solver ();
    Psolt->par_linear_solver (Gtt,Kmat,tdlhs,rhs,Outt,Mesprt);
    



    //fprintf (Outt,"\n\n MATICE SOUSTAVY PO RESENI \n");
    //Kmat->printmat (Outt);
    
    //fprintf (Outt,"\n\n\n\n\n kontrola reseni soustavy");
    //for (long ijk=0;ijk<n;ijk++){
    //fprintf (Outt,"\n tdlhs %ld      %le",ijk,tdlhs[ijk]);
    //}
    
    
    
    for (j=0;j<n;j++){
      s=(1.0-alpha)*d[j]+alpha*tdlhs[j];
      lhs[j]+=dt*s;
    }
    
    
    //fprintf (Outt,"\n\n\n\n\n kontrola obnovenych promennych");
    //for (long ijk=0;ijk<n;ijk++){
    //fprintf (Outt,"\n lhs %ld   %le",ijk,lhs[ijk]);
    //}

    //  nodes - integration points interpolation
    approximation ();
    



    Tp->time = Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    i++;
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();

    printf ("\n Tp->time   %e",Tp->time);
    printf ("\n i          %ld",i);
    
  }while(Tp->time<end_time);



  delete [] p;  delete [] d;
  //fclose (gr); 

  print_closet();
  
  Psolt->computation_stat_print (Outt);
  
  
}

