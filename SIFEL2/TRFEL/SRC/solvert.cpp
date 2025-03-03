#include "solvert.h"
#include "spsolvert.h"
#include "nspsolvert.h"
#include "npsolvert.h"
#include "nnpsolvert.h"
#include "cpnpsolvert.h"
#include "cpnnpsolvert.h"
#include "dnpsolvert.h"
#include "dnnpsolvert.h"
#include "homogtrans.h"
#include "globalt.h"
#include "globmatt.h"
#include "libtrace.h"
#include "trfelinit.h"


void solve_trfel_problem ()
{
  /// termitovo zkontrolovat-smazat
  if (Tp->stochasticcalc==0 && Tp->adaptivityflag==0)
    solve_trfel_deterministic_problem ();
  if (Tp->stochasticcalc    && Tp->adaptivityflag==0)
    solve_trfel_stochastic_problem ();
  if (Tp->stochasticcalc==0 && Tp->adaptivityflag)
    solve_trfel_adaptivity_problem ();
  if (Tp->stochasticcalc    && Tp->adaptivityflag)
    print_err("stochastic and adaptivity mix not supported",__FILE__,__LINE__,__func__);
}

void solve_trfel_deterministic_problem ()
{
  switch (Tp->tprob){
  case stationary_problem:{
    if(Tp->homogt == 1)
      transport_homogenization();
      //transport_homogenization_old_old();
    else{
      solve_stationary_problem ();
      //solve_radiation_stationary_problem ();
    }
    break;
  }
  case nonlinear_stationary_problem:{
    if(Tp->homogt == 1)
      transport_homogenization();
      //transport_homogenization_old();
    else{
      solve_nonlinear_stationary_problem_pokus ();
      //solve_nonlinear_stationary_problem ();
      //solve_nonlinear_stationary_problem_old ();
    }
    break;
  }
  case nonstationary_problem:{
    solve_nonstationary_problem ();
    break;
  }
  case nonlinear_nonstationary_problem:{
    //solve_nonlinear_nonstationary_problem ();
    solve_nonstationary_problem ();
    break;
  }
  case growing_np_problem:{
    //solve_nonstationary_growing_problem ();
    //solve_nonstationary_growing_vform ();
    solve_nonstationary_problem ();
    break;
  }
  case growing_np_problem_nonlin:{
    //solve_nonstationary_growing_problem_nonlin ();
    solve_nonstationary_problem ();
    break;
  }
  case discont_nonstat_problem:{
    solve_discont_nonstationary_problem ();
    break;
  }
  case discont_nonlin_nonstat_problem:{
    solve_discont_nonlin_nonstationary_problem ();
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}

void solve_trfel_stochastic_problem ()
{
  print_err("stochastic transport problems are not implemented at this time",__FILE__,__LINE__,__func__);
}




/// **********************************************************************************************************
/// *** TERMITOVO
/// **********************************************************************************************************

void print_VTK_nodedata (char *filename, gtopology *gt, const double *data)
{
  long i;
  
  FILE *out = fopen (filename,"w");   if (out==NULL)  { print_err("test file has not been opened", __FILE__, __LINE__, __func__);  exit (1); }
  
  /// HEAD
  fprintf (out, "# vtk DataFile Version 3.0\n");
  fprintf (out, "testing VTK file written by Termitovo fci in SIFEL\n");
  fprintf (out, "ASCII\n");
  fprintf (out, "DATASET POLYDATA\n");
  fprintf (out, "POINTS %ld float\n", gt->nn);
  
  for (i=0; i<gt->nn; i++)
    fprintf (out, "%13.6f %13.6f %13.6f\n", gt->gnodes[i].x, gt->gnodes[i].y, gt->gnodes[i].z);
  
  fprintf (out, "POLYGONS %ld %ld\n", gt->ne, 4*gt->ne);
  
  for (i=0; i<gt->ne; i++)
    fprintf (out, "3 %ld %ld %ld\n", gt->gelements[i].nodes[0], gt->gelements[i].nodes[1], gt->gelements[i].nodes[2]);
  
  fprintf (out, "POINT_DATA %ld\n", gt->nn);
  fprintf (out, "SCALARS displacement float 2\n");
  fprintf (out, "LOOKUP_TABLE default\n");
  
  for (i=0; i<gt->nn; i++)
    fprintf (out, "%13.6f %13.6f\n", data[i], data[i+gt->nn]);
  
  fclose (out);
}

void newmeshgen (long i)
{
  int ret;
  char procname[1023];
  char genbin[255];    sprintf (genbin, "/home/dr/Bin/T3d");
  char prepbin[255];   sprintf (prepbin, "./transprep");
  double d = 0.6;         // velikost prvku, zatim natvrdo, bude se zadavat nebo bude nejlepe primo v this.t3d.in
  
  const char *filename = Adat->give_filename();
  const char *ni       = Adat->give_ni();
  int w                = Adat->give_niwidth();
  
  /// MESH GENERATION
  if (Mesprt) fprintf (stdout," ADAPTIVITY: T3D generation\n");
  sprintf (procname, "%s  -i %s.T3d.in  -o %s.T3d.out  -m %s.T3d.bgm%s  -d %g  -p 264 -r 1 -k %d", genbin, filename, filename, filename, ni, d, 1); // Tt->give_degree(0));
  fprintf (stdout,"\n %s\n",procname);
  if (Mesprt<2)  strcat(procname, "  >> adaptivity.log 2>&1");
  if (Mesprt>1)  fprintf (stdout, "%s\n", procname);
  ret = system (procname);  if (ret != 0) { print_err("Mesh generation failed", __FILE__, __LINE__, __func__);  exit (1); }
  
  /// SIFEL.IN GENERATION
  if (Mesprt) fprintf (stdout," ADAPTIVITY: sifel.in generation\n");
  sprintf (procname, "%s %s.sifel.pr %s.sifel.in.%0*ld", prepbin, filename, filename, w, i);
  if (Mesprt<2)  strcat(procname, "  >> adaptivity.log 2>&1");
  if (Mesprt>1)  fprintf (stdout, "%s\n", procname);
  ret = system (procname);  if (ret != 0) { print_err("TRFEL PREP failed", __FILE__, __LINE__, __func__);  exit (1); }
  
  
  /// FILE BACKUP
  sprintf (procname,"mv  %s.T3d.out  %s.T3d.out.%0*ld", filename, filename, w, i);     system (procname);
  //sprintf (procname,"mv  %s.T3d.bgm  %s.T3d.bgm.%0*ld", filename, filename, w, i-1);   system (procname);
}

void solve_trfel_adaptivity_problem ()
{
  /// FIRST ADAPTIVE STEP
  fprintf (stdout,"\n\n ADAPTIVITY: start with mesh 0\n");
  solve_trfel_deterministic_problem ();
  
  switch (Tp->tprob) {
  /// STATIONARY PROBLEM
  case stationary_problem:
    if (Adat->answer)  fprintf (stdout,"\n ADAPTIVITY: new mesh required\n");
    else               fprintf (stdout,"\n ADAPTIVITY: new mesh NOT required\n");
    break;
  /// NONSTATIONARY_PROBLEM
  case nonstationary_problem: {
    gtopology *Gtt_old;
    adaptivityt *Adat_old;
    
    int argc = 2;
    const char *argv[2];
    argv[0] = "trfel";
//    argv[0] = get_program_name();
    char infile[255];
    argv[1] = infile;
    
    int i = 1;
    while (Adat->answer) {
      fprintf (stdout,"\n ADAPTIVITY: switch form mesh %d to mesh %d\n", i-1, i);
      
      /// OLD MESH procedures: mesh generation, ...
      newmeshgen (i);
      sprintf (infile, "%s.sifel.in.%0*d", Adat->give_filename(), Adat->give_niwidth(), i);
      if (Gtt->nadjelnod == NULL)
	Gtt->adjacelem (Outt);
      
      
      /// BACKUP
      Adat->statedata_backup();
      //
      Gtt_old    = Gtt;     Gtt   = NULL;
      Adat_old   = Adat;    Adat  = NULL;
      
      
      /// DELETE OLD PROBLEM
      fclose(Outt);
      delete_globt ();
      
      
      /// NEW MESH procedures
      trfel_init (argc, argv);
      if (Gtt->nadjelnod == NULL)
	Gtt->adjacelem (Outt);
      
      /// set up
      Tp->timecont.be_copy_of(*Adat_old->tctrl);   //  nemam to dat do npsolver ??
      Tp->time = Tp->timecont.time;               //  -- "" --
      Adat->initialize(i);
      
      
      /// DATA TRANSFER
      Adat->statedata_transfer(Adat_old, Gtt_old);
      
      /// check print
      char tmpfile[255+12];
      sprintf (tmpfile, "%s.r_old.vtk",   infile);   print_VTK_nodedata (tmpfile, Gtt_old, Adat_old->give_r());
      sprintf (tmpfile, "%s.r_new.vtk",   infile);   print_VTK_nodedata (tmpfile, Gtt,     Adat    ->give_r());
      sprintf (tmpfile, "%s.rdr_old.vtk", infile);   print_VTK_nodedata (tmpfile, Gtt_old, Adat_old->give_rdr());
      sprintf (tmpfile, "%s.rdr_new.vtk", infile);   print_VTK_nodedata (tmpfile, Gtt,     Adat    ->give_rdr());
      
      /// delete
      delete Gtt_old;   Gtt_old = NULL;
      delete Adat_old;  Adat_old = NULL;
      
      
      // ??????????????????????????????????????????????????????????????????
      // toto jsem mel v mefelu, je to v trfelu treba??
      /// mefel_right_hand_side (lcid,Lsrs->rhs);
      
      
      /// SOLVE
      solve_trfel_deterministic_problem ();
    
      
      i++;
      
      /// docasne aby to delalo jen jeden krok
      ///Adat->answer = 0;
    }    
    
    break;
  }
  default:
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);  exit (1);
  }
}
