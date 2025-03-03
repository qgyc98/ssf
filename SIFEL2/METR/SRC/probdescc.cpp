#include "probdescc.h"
#include "globalc.h"
#include "probdesc.h"
#include "gmatrix.h"
#include "intools.h"



probdescc::probdescc (void)
{
  kwdsw = nokwd;

  niilnr=0;  errnr=0.0;
  alpha=0.0;  time=0.0;
  zero = 1.0e-20;

  
  bb = notdefbb;

  path = filename = suffix = NULL;
  //  smaller number of integration points
  savemode=1;
  
  // clean matrices
  cleanmatrix = (coupcleanmatrices) 0;

  //  type of nonlinear solver
  fcsolv = (coupsolver) 0;
  
  //  type of residuum
  restype = (residuumtype) 0;
  
  //  solver of system of linear algebraic equations
  ssle = new slesolv();

  cgravity = grc_no;
  memset(gr, 0, sizeof(*gr)*3);
  // type of the passing data between modules (TRFEL, MEFEL) in coupled problems
  // pass_by_closest_ip=1 -> data are copied to nodes from the closest int. point, requires the same number of elements in both modules - minimum additional memory requirements, fastest way, may be inaccurate for certain combination of elements
  // pass_by_nodes_comp=2 -> data are calculated at nodes with the help of state varibales from the closest int. point, requires the first n elements in both modules to be the same in shape and same ordering of the first m nodes,
  //                         where n is the minimum number of elements in modules, m is the minimum number of nodes on elements - minimum additional memory requirements, slower than type 1, better accuracy than type 1
  // pass_by_aux_ip = 3 -> data are calculated in auxiliary int. points, the number of nodes and elements in meshes of particular modules can be independent - may leed to large memory requirements, slower than type 1, exact.
  // pass_by_copy_ip = 4 -> data are copied directly between int. points, the SAME int. point must be used on both meshes for the first Tt->ne elements - fast, exact, minimu memory requirements.
  //dpt = pass_by_copy_ip;
  dpt = pass_by_aux_ip;
  // the maximum number of itartions in the computation of auxiliary int. point coordinates
  aip_ni = 100;
  // the maximum distance of two points at which they may be assumed to be identical in the computation of auxiliary int. point coordinates
  aip_err = 1.0e-6;

  incr_step=0;
  inner_step=0;


}

probdescc::~probdescc (void)
{
  delete [] path; delete [] filename; delete [] suffix;
  delete ssle;
}

/**
   function reads basic data about solved problem
   
   @param in - input stream
   
   23.11.2002
*/
void probdescc::read (XFILE *in)
{
  char str_tmp[1001];//auxiliary string
  const char *pkwd;
  
  // problem description
  xfscanf (in," %a",name);
  
  //  problem type
  xfscanf (in,"%k%m","problemtype",&problemtypec_kwdset,(int*)&tprob);

  if (tprob !=fully_coupled_material){
    // line with MEFEL input filename and optional -kwd switch
    xfscanf (in," %a",str_tmp);
    // processing of optional keyword switch -kwd={0,1,2,3,4}
    pkwd = strstrcis(str_tmp, "-kwd=");
    if (pkwd)
    {
      // pickup MEFEL input filename from the beginning of line until the first whitespace occurence
      memset(minfile, 0, 1001*sizeof(*minfile));
      strncpy (minfile, str_tmp, (strchr(str_tmp,' ')-str_tmp)/sizeof(char));
      Mp->kwdsw = pkwd_sw(-1);
      sscanf(pkwd+5, "%d", (int *)&Mp->kwdsw);
      if ((Mp->kwdsw < 0) || (Mp->kwdsw > 4))
      {
	fflush(stdout);
	print_err("unknown -kwd value %d is obtained\n", __FILE__, __LINE__, __func__, Mp->kwdsw);
	abort();
      }
    }
    else
    {
      strcpy (minfile, str_tmp);
      Mp->kwdsw = nokwd;
    }
    
    // line with TRFEL input filename and optional -kwd switch
    xfscanf (in," %a",str_tmp);
    // processing of optional keyword switch -kwd={0,1,2,3,4}
    pkwd = strstrcis(str_tmp, "-kwd=");
    if (pkwd)
    {
      memset(tinfile, 0, 1001*sizeof(*tinfile));
      // pickup TRFEL input filename from the beginning of line until the first whitespace occurence
      strncpy (tinfile, str_tmp, (strchr(str_tmp,' ')-str_tmp)/sizeof(char));
      Tp->kwdsw = pkwd_sw(-1);
      sscanf(pkwd+5, "%d", (int *)&Tp->kwdsw);
      if ((Tp->kwdsw < 0) || (Tp->kwdsw > 4))
      {
	fflush(stdout);
	print_err("unknown -kwd value %d is obtained\n", __FILE__, __LINE__, __func__, Tp->kwdsw);
	abort();
      }
    }
    else
    {
      strcpy (tinfile, str_tmp);
      Tp->kwdsw = nokwd;
    }
  }
  
  //  detail messages
  Mesprc=0;
  xfscanf (in,"%k%ld","mesprc",&Mesprc);
  if (Mesprc==1)  fprintf (stdout,"\n detail information will be printed");
  else fprintf (stdout,"\n only important messages will be printed");
 
  //  Babuska-Brezzi
  xfscanf (in,"%k%m","BB",&babuskabrezzi_kwdset,(int*)&bb);
  //  type of the passing data between modules
  xfscanf (in,"%k%m","Datapasstype",&datapasstype_kwdset,(int*)&dpt);
  if (Mesprc)
    fprintf(stdout, "\n type of the data passing between TRFEL and MEFEL: %s", datapasstype_kwdset.get_str(dpt));
  //  transported matter
  xfscanf (in,"%k%m","transmatter",&transmatterc_kwdset,(int*)&tmatt);
  //  names of transported media
  xfscanf (in,"%k%m","mednames",&mednamesc_kwdset,(int*)&mednam);
  //  passing data type
  // xfscanf(in, "%k%m", "data_passing_type", &datapasstype_kwdset, (int*)&dpt);
  // if (dpt == pass_by_aux_ip)
  //   xfscanf(in, "%k%ld%k%le", "aux_ip_ni", &aip_ni, "aux_ip_err" , &aip_err);
  
  switch (tprob){
  case fully_coupled_material:
  case fully_coupled_mech_trans:{
    
    //  gravity acceleration is taken into account
    xfscanf (in,"%k%m","gravityacceleration",&gravityaccelerationc_kwdset,(int*)&cgravity);
    if (Mesprc==1) fprintf (stdout,"\n Gravity acceleration is taken into account in METR");

    if(cgravity == grc_yes){
      xfscanf (in,"%lf %lf %lf",gr+0,gr+1,gr+2);
    }

    //  time controller
    //  start time, end time, number of important times
    //  list of important times, time increment function
    timecon.read (in);
    
    //  alpha coefficient
    xfscanf (in,"%lf",&alpha);
    
    //  type of nonlinear solver
    //  tnlinsol = 1 - arc-length
    //  tnlinsol = 2 - Newton-Raphson method
    xfscanf (in,"%d",(int*)&tnlinsol);
    //  maximum number of iteration in one increment
    //  maximum norm of residual vector
    xfscanf (in,"%ld %lf",&niilnr,&errnr);
    
    //  type of nonlinear solver
    //  fcsolv = 1 - full Newton
    //  fcsolv = 2 - modified Newton
    //  restype = 1 - fluxes
    //  restype = 2 - comparison of lhs and rhs
    xfscanf (in,"%k%m%k%m", "coupledsolver", &coupsolver_kwdset, (int*)&fcsolv, "resid_type", &residuumtype_kwdset, (int*)&restype);
    
    //  type of storage of zero-order and first order matrices
    //  zero-order matrix = stiffness and conductivity matrices
    //  first-order matrix = capacity matrix
    xfscanf (in,"%k%m%k%m","zeroordermat",&storagetype_kwdset,(int*)&tstord0,"firstordermat",&storagetype_kwdset,(int*)&tstord1);
    
    //  data about solver of system of linear equations
    ssle->read (Gtu,in,Mesprc);
    
    break;
  }
  case growing_par_coupl_mech_trans:
  case par_coupl_mech_trans:{
    
    //  time controller
    //  start time, end time, number of important times
    //  list of important times, time increment function
    timecon.read (in);
    
    //  linearc=0
    //  fullnewtonc=1
    //  modnewtonc=2
    xfscanf (in,"%k%m","nonlinsolvertype",&nonlinsolvertypec_kwdset,(int*)&tnlinsol);

    xfscanf (in,"%ld %lf",&niilnr,&errnr);
    
    //  type of coupled solver
    xfscanf (in,"%k%m","cleanmatrices",&coupcleanmatrices_kwdset,(int*)&cleanmatrix);
    
    //  this problem requires storage of actual nodal values on nodes
    Tp->nvs=1;
    //  this problem requires storage of actual nodal values on nodes
    Tp->pnvs=1;
    //  this problem requires storage of initial nodal values on nodes
    Tp->invs=1;
    //  this problem requires storage of time derivatives of nodal values on nodes
    //tdnvs=1;

    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
    abort ();
  }
  }
}


/**  function prints basic data about solved problem
     JK, 25.9.2001
     
     @param in - output file
*/
void probdescc::print (FILE *out)
{
  //  problem name
  fprintf (out, "%s\n", name);

  //  problem type
  fprintf (out,"%d\n", int(tprob));
  
  if (Cp->tprob!=fully_coupled_material){
    //  MEFEL input file name
    fprintf (out,"%s\n",minfile);

    //  TRFEL input file name
    fprintf (out,"%s\n",tinfile);
  }
  
  //  detail messages
  fprintf (out,"%ld\n", Mesprc);

  // BB condition
  fprintf (out,"%d\n", bb);

  //  type of the passing data between modules
  fprintf (out,"%d\n",dpt);
  
  //  transported matter
  fprintf (out,"%d\n", tmatt);

  //  names of transported media
  fprintf (out, "%d\n\n", mednam);
     
  //  type of the passing values between modules
  // fprintf (out, " %d\n", dpt);
  // if (dpt == pass_by_aux_ip)
  //  fprintf(out, "%ld %le\n", aip_ni, aip_err);


  switch (tprob){
  case fully_coupled_mech_trans:{
    
    //  gravity acceleration
    fprintf (out,"%d # gravity\n",cgravity);
    if(cgravity == grc_yes){
      fprintf (out,"%le %le %le\n",gr[0],gr[1],gr[2]);
    }

    timecon.print (out);

    fprintf (out," %d",tnlinsol);

    fprintf (out, " %lf \n", alpha);

    fprintf (out," %ld %lf\n",niilnr,errnr);
    
    fprintf (out," %d %d\n",fcsolv,restype);

    fprintf (out," %d %d\n",tstord0,tstord1);
    
    ssle->print (out);
    
    fprintf (out,"\n");
    
    break;
  }

  case growing_par_coupl_mech_trans:
  case par_coupl_mech_trans:{
    
    timecon.print (out);
    
    
    fprintf (out," %d",tnlinsol);

    fprintf (out," %ld %lf\n",niilnr,errnr);
    
    fprintf (out," %d\n",fcsolv);

    break;
  }

  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }  
}
