#include "outresfilem.h"
#include "iotools.h"
#include "outquantm.h"
#include "mechprint.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "tensorcname.h"
#include "vectorcnamem.h"

#include "xfile.h"



outresfilem::outresfilem() : outfbn{}, outfmn{}
{
  outfbnl = 0;
  outfbn[0] = '\0';
  rffmt = resfmt_no;
  nqnt = 0;
  qnt = NULL;
  outf = NULL;
  outfmn[0] = '\0';
}



outresfilem::~outresfilem()
{
  delete [] qnt;
}



/**
  The function reads the given output file description.

  @param[in] in - pointer to the opened text file

  @retval 0 - on success,
  @retval 1 - error in the reading of base output file name,
  @retval 2 - error in the reading of output file format,
  @retval 3 - error in the reading of time step selection,
  @retval 4 - error in the reading of load case selection,
  @retval 5 - error in the reading of number of required quantities,
  @retval 6 - error in the reading of some quantity description.

  Created by Tomas Koudelka, 10.2023  
*/
long outresfilem::read(XFILE *in)
{
  long i, ret;
  
  // read name of the given result output file
  ret = xfscanf(in, "%s", outfbn);
  if (ret != 1){
    print_err("cannot read base file name from file", __FILE__, __LINE__, __func__);
    return 1;
  }
  outfbnl = strlen(outfbn);

  // read format of the given result file
  ret = xfscanf(in, "%k%m", "res_file_fmt", &resfilefmt_kwdset, &rffmt);
  if (ret != 2){
    print_err("cannot read result file format", __FILE__, __LINE__, __func__);
    return 2;
  }

  // read selection of the required time steps in which the output to the given file will be performed
  ret = selstep.read(in);
  if (ret){
    print_err("cannot read selection of time step", __FILE__, __LINE__, __func__);
    return 3;
  }

  // read selection of the required load cases in which the output to the given file will be performed
  ret = sellc.read(in);
  if (ret){
    print_err("cannot read load case selection", __FILE__, __LINE__, __func__);
    return 4;
  }

  // read number of required quantities 
  ret = xfscanf(in, "%k%ld", "nquant", &nqnt);
  if (ret != 2){
    print_err("cannot read number of required quantities", __FILE__, __LINE__, __func__);
    return 5;
  }

  // read particular quantities
  qnt = new outquantm[nqnt];
  for(i=0; i<nqnt; i++){
    ret = qnt[i].read(in);
    if (ret != 2){
      print_err("cannot read %ld. quantity description", __FILE__, __LINE__, __func__, i+1);
      return 6;
    }
  }
  
  return 0;
}



/**
  The function prints the given output file description.

  @param[in,out] out - pointer to the opened text file

  @return The function does not return anything, but changes the content of the output file given by the argument out.

  Created by Tomas Koudelkam, 10.2023
*/
void outresfilem::print(FILE *out)
{
  long i;
  
  // print name of the given result output file
  fprintf(out, "%s\n", outfbn);

  // print format of the given result file
  fprintf(out, "%d\n", int(rffmt));

  // print selection of the required time steps in which the output to the given file will be performed
  selstep.print(out);

  // read selection of the required load cases in which the output to the given file will be performed
  sellc.print(out);

  // read number of required quantities 
  fprintf(out, "%ld\n", nqnt);

  // read particular quantities
  for(i=0; i<nqnt; i++)
    qnt[i].print(out);
}



/**
  The function prints initial header to the output file for 
  particular result formats.

  @param[in] istep - time/load step identifier
  @param[in] time - actual time or load coefficient

  @return The function does not return anything, but changes content of the output file outf.

  Created by Tomas Koudelka, 10.2023
*/
void outresfilem::print_header(long istep, double time)
{
  long i, j;
  char *label = NULL;
  strastrestate ssst = strastrestate(0);
  
  switch(rffmt){
    case resfmt_gid:      
      fprintf(outf, "GiD Post Results File 1.0\n");
      break;
    case resfmt_vtk:
      print_vtk_header(outf, istep, time); //print header, points, nodes, elements
      break;
    case resfmt_plain_out:
      fprintf(outf, "%15s ****  *  ****  ****  *\n", " ");
      fprintf(outf, "%15s *     *  *     *     *\n", " ");
      fprintf(outf, "%15s  *    *  ***   ***   *\n", " ");
      fprintf(outf, "%15s   *   *  *     *     *\n", " ");
      fprintf(outf, "%15s****   *  *     ****  ****  MEFEL OUTPUT\n", " ");
      // print specific header string
      fprintf(outf, "\n%s\n", Mp->name);
      fprintf(outf, "\n\n\n\n\n");
      break;
    case resfmt_diag_dat:
      fprintf(outf, "#");
      for (i=0; i<nqnt; i++){
        if ((qnt[i].defqlabel == true) && (label == NULL))
          label = new char[outquantm::maxl_qlabel];
        ssst = qnt[i].give_pnt_ssst(qnt[i].selid.id1[0]);
        for(j=0; j<qnt[i].ncmp; j++){
          if (qnt[i].seliq.presence_id(j)){
            if (qnt[i].defqlabel == false){
              fprintf(outf, " %s", qnt[i].qlabel);
              if (qnt[i].qr == tensq){
                fprintf(outf, "_%s", tensor_cname::tens_indstr(ssst, j));
              }
              if ((qnt[i].qr == vectq) && (qnt[i].mqn < tempr_strain_q)) // global nodal vectors                
                fprintf(outf, "_%s", vector_cnamem::vect_indstr(ssst, j));
              if ((qnt[i].qr == vectq) && (qnt[i].mqn >= tempr_strain_q)) // other, nonmechq
                fprintf(outf, "_%ld", j+1);
            }
            else{
              ssst = qnt[i].give_pnt_ssst(qnt[i].selid.id1[0]);
              qnt[i].give_quant_comp_deflabel(qnt[i].mqn, qnt[i].reftensq, ssst, j, outquantm::maxl_qlabel, label);
              fprintf(outf, " %s", label);
            }
          }
        }
      }
      fprintf(outf, "\n");
      break;
    default:
      print_err("unknown type of file format (%d) is required in the output file '%s'.",
                __FILE__, __LINE__, __func__, int(rffmt), outfmn);
      abort();
  }
  delete [] label;
}



/**
  The function creates modifed file name of the given output file
  with the help of base file name outfbn decorated optionally by the additional
  filename modifiers according to the following rules:
  - if the argument step id given by the argument stepid is non-negative then 
    this id is appended to the base file name
  - if the parameter set identfier in the stochastic calculations given by the argument stochid is 
    non-negative then this id is  appended to the base file name
  - if the base file name has no extension then the default one is appended

  @param[in] stochid - index of the parameter set in the case of stochastic calculation, otherwise stochid=-1
  @param[in] stepid - time step index, the base file name is not decorated by stepid if stepid=-1

  Created by Tomas Koudelka, 10.2023
*/
void outresfilem::get_modified_fname(long stochid, long stepid)
{
  char *path, *name, *ext;
  long pathl, namel;
  long pc;

  // decompose base file name into path, name and extension
  filename_decomp_ptr(outfbn, path, name, ext);
  if (path)
    pathl = name - path;
  else
    pathl = 0;

  if (ext){
    namel = ext - name;
    //extl = strlen(ext);
  }
  else{
    namel = outfbnl;
    //extl = 0;
  }

  // take path and file name
  pc = sprintf(outfmn, "%*s", int(pathl+namel), outfbn);
  // add optional identifier of stochastic parameter set
  if (stochid >= 0)
    pc += snprintf(outfmn+pc, fnl-pc, ".%ld", stochid);
  // add optional identifier of stochastic parameter set
  if (stepid >= 0)
    pc += snprintf(outfmn+pc, fnl-pc, ".%ld", stepid);
  
  if (ext == NULL){
    // no file extension is specified -> add default one
    const char *dext = give_fmt_dext(rffmt);
    pc += snprintf(outfmn+pc, fnl-pc, ".%s", dext);
  }
}



/**
  The function prints all required quantities to the given output file.

  @param[in] lcid - required load case id,
  @param[in] time - actual time or actual load coefficient,
  @param[in] istep - actual step id,
  @param[in] lcsarr - array of the local coordinte systems,
  @param[in] gnv - actual nodal values of quantities related to nodes only
  @param[in] idn1 - shift index of node numbers
  @param[in] ide1 - shift index of element numbers

  @return The function does not return anything but changes content of the file given by outf data member.

  Created by Tomas Koudelka, 10.2023
*/
void outresfilem::print_quants(long lcid, double time, long istep, lcoordsys *lcsarr,
                               gnodvalvm &gnv, long idn1, long ide1)
{
  long i;
  
  for(i=0; i<nqnt; i++){ // loop over all quantities in the file
    qnt[i].outf = outf;
    qnt[i].outfname = outfmn;
    qnt[i].idn1 = idn1;      
    qnt[i].ide1 = ide1;
    qnt[i].print_quant(lcid, time, istep, lcsarr, gnv, rffmt);
  }
}



/**
  The function returns default file extension for the required result file format.

  @param[in] fmt - required result file format

  @return The function returns pointer to the constant string with result file extension.

  Created by Tomas Koudelka, 10.2023
*/
const char*outresfilem::give_fmt_dext(resfilefmt fmt)
{
  switch(fmt){
    case resfmt_no:
      return "";
    case resfmt_open_dx:
      return "";
    case resfmt_femcad:
      return "";
    case resfmt_gid:
      return ".res";
    case resfmt_vtk:
      return ".vtk";
    case resfmt_plain_out:
      return ".out";
    case resfmt_diag_dat:
      return ".dat";
    default:
      print_err("unknown type %d of the result file format", __FILE__, __LINE__, __func__);
      abort();
  }
  return "";
}
