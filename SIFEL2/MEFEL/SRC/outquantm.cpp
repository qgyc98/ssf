#include "outquantm.h"
#include "global.h"
#include "mechmat.h"
#include "mechtop.h"
#include "probdesc.h"
#include "elemswitch.h"
#include "node.h"
#include "element.h"
#include "vector.h"
#include "matrix.h"
#include "vecttens.h"
#include "tensorcname.h"
#include "vectorcnamem.h"
#include "alias.h"
#include "galias.h"
#include "lcoordsys.h"

#include <stdio.h>
#include <stdlib.h>

double **outquantm::qvv = NULL;
double  *outquantm::qvs = NULL;
long outquantm::maxnpnt = 0;
long outquantm::maxncmp = 0;
const char *outquantm::namepar[] = {"x", "y", "z", "t", "q"};

// setup of the Jacobi method for the determination of principal tensor values
long outquantm::ni_jac=30;
double outquantm::err_jac=1.0e-10;
double outquantm::zero=1.0e-15;;


outquantm::outquantm()
{
  mqn = nomech_q;
  reftensq = nomech_q;
  qlabel[0] = '\0';
  tpnt = atnode;
  tnpnt = 0;
  nudp = 0;
  ksi = eta = zeta = NULL;
  qr = qrf = scalq;
  qstf = no;
  lcs = 0;
  ncmp = ncmpr = 0;
  outfname = NULL;
  outf = NULL;
  ide1 = 0;
}



outquantm::~outquantm()
{
  delete [] ksi;
  delete [] eta;
  delete [] zeta;
}



/**
  The function reads description of printed quantity to the file for the graphic postprocessor.
  @param in[in] - pointer to the opened input file

  @retval 0 - on success
  @retval 1 - in the case of reading error

  Created by Tomas Koudelka, 07.2021
*/
long outquantm::read(XFILE *in)
{
  long i, n;
  
  // read name id of required quantity
  xfscanf(in, "%k%m","qname", &mechquant_kwdset, &mqn);
  switch(mqn){
    case first_inv:
    case second_inv:
    case third_inv:
    case tensor_norm:
    case tensdeviator:
      xfscanf(in, "%k%m","reftens_qname", &mechquant_kwdset, &reftensq);
      if ((reftensq == other_q) || (reftensq == nonmech_q)){
        xfscanf(in, "%k%m","reftens_strastre", &mechquant_kwdset, &other_tens_strastre);
      }
      break;
    default:
      break;
  }

  // read string with label of required quantity, use "@" for automatic label
  xfscanf(in, "%k%255s", "qlabel", qlabel);
  if (qlabel[0] == '@')    defqlabel = true;

  // read type of points where quantity values willl be printed out
  xfscanf(in, "%k%m", "pnt_type", &nodip_kwdset, &tpnt);

  // read selection of node/elem numbers where the quantity will be printed out 
  xfscanf(in, "%k", "pnt_sel");
  selid.read(in);

  switch (tpnt){
    case atnode: // no additional data are necessary in the case of nodes
      tnpnt = n = Mt->nn;
      break;
    case atip: // read selection of integration points on selected elements (local ordering numbers)
      tnpnt = n = Mm->tnip;
      seleip.read(in);
      break;
    case atxyz: // read set of user defined points (UDP) for selected elements
      xfscanf(in, "%k%ld", "num_udp", &nudp); // number of UDPs 
      ksi  = new double[nudp];
      eta  = new double[nudp];
      zeta = new double[nudp];
      for(i=0; i<nudp; i++){ // reading of natrual coordinates of UDPs
        xfscanf(in, "%le %le", ksi+i, eta+i);
        if (Mp->dim == 3)
          xfscanf(in, "%le", zeta+i);
        else
          zeta[i] = 0.0;
      }
      tnpnt = n = Mt->ne*nudp;
      break;
    default:
      print_err("unknown type of point is required (%d)", __FILE__, __LINE__, __func__, int(tpnt));
      abort();
  }
  // set the first dimension of auxiliary arrays qvs and qvv
  if (n > maxnpnt)  maxnpnt = n;

  // get quantity represenation and its number of components according to mechmat and number of quantity components
  qr = Mm->give_quant_rep(mqn, reftensq, ncmp, tsn, tti);
  ncmpr = ncmp;
  switch(qr){
    case scalq:
      seliq.st = sel_all;
      break;
    case vectq:
      xfscanf(in, "%k", "vect_comp_sel");
      // read component selection for vector types of quantities
      seliq.read(in);
      ncmpr = seliq.give_nselcomp(ncmp); // reduce the maximum according to the number of selected components      
      // read quatity representation in which will be printed out (scalar/vector/tensor)
      xfscanf(in, "%k%m", "quant_out_rep", &quantrep_kwdset, &qrf);
      break;
    case tensq:
      xfscanf(in, "%k", "vect_comp_sel");
      // read component selection for tensor types of quantities
      seliq.read(in);
      ncmpr = seliq.give_nselcomp(ncmp); // reduce the maximum according to the number of selected components      
      // read quatity representation in which will be printed out (scalar/vector/tensor)
      xfscanf(in, "%k%m", "quant_out_rep", &quantrep_kwdset, &qrf);
      break;
    default:
      print_err("unknown type of quantity (%d) is required", __FILE__, __LINE__, __func__, int(qr));
      abort();      
  }
  // set the second dimension of auxiliary array qvv
  if (ncmpr > maxncmp)  maxncmp = ncmpr;

  // function for scaling/tranformation of the quantity (constant function means the scaling factor, pars
  xfscanf(in, "%k%m", "scaling", answertype_kwdset, &qstf);
  if (qstf == yes)
    qst.read(in); // read scaling function, const type of function means scaling factor, parser is used for other transformations

  // read id of transformation for vector or tensor quantity representation
  switch(qr){
    case scalq:
      break;
    case vectq:
    case tensq:   // read id of local coordinate system (LCS)
      xfscanf(in, "%k%ld", "lcs", &lcs);
      break;
    default:
      print_err("unknown representation of quantity (%d) is required", __FILE__, __LINE__, __func__, int(qr));
      abort();      
  }
  return 0;
}



/**
  The function prints description of printed quantity to the text file.
  @param out[in] - pointer to the opened output file

  @return The function does not return anything.

  Created by Tomas Koudelka, 07.2021
*/
void outquantm::print(FILE *out)
{
  long i;
  
  // print name id of required quantity
  fprintf(out, "%d ", int(mqn));

  // print string with label of required quantity
  fprintf(out, "%s ", qlabel);

  // print type of points where quantity values willl be printed out
  fprintf(out, "%d ", int(tpnt));

  // print selection of node/elem numbers where the quantity will be printed out 
  selid.print(out);

  switch (tpnt){
    case atnode: // no additional data are necessary in the case of nodes
      break;
    case atip: // read selection of integration points on selected elements (local ordering numbers)
      seleip.print(out);
      break;
    case atxyz: // read set of user defined points (UDP) for selected elements
      fprintf(out, "%ld\n", nudp); // number of UDPs 
      for(i=0; i<nudp; i++){ // reading of natrual coordinates of UDPs
        fprintf(out, "%le %le", ksi[i], eta[i]);
        if (Mp->dim == 3)
          fprintf(out, " %le", zeta[i]);
        fprintf(out, "\n");
      }      
      break;
    default:
      print_err("unknown type of point is required (%d)", __FILE__, __LINE__, __func__, int(tpnt));
      abort();
  }

  // print quatity format in which will be printed out (scalar/vector/tensor)
  fprintf(out, "%d ", int(qrf));
  switch(qr){
    case scalq:
      break;
    case vectq:
    case tensq: // print component selection for vector and tensor types of quantities
      seliq.print(out);
      break;
    default:
      print_err("unknown type of quantity (%d) is required", __FILE__, __LINE__, __func__, int(qr));
      abort();      
  }

  // function for scaling/tranformation of the quantity (constant function means the scaling factor, pars
  fprintf(out, "%d ", int(qstf));
  if (qstf == yes)
    qst.print(out); // print scaling function, const type of function means scaling factor, parser is used for the complex transformations

  // print id of transformation for vector or tensor quantity represenation
  switch(qr){
    case scalq:
      break;
    case vectq:
    case tensq:   // read id of local coordinate system (LCS)
      fprintf(out, "%ld", lcs);
      break;
    default:
      print_err("unknown type of quantity (%d) is required", __FILE__, __LINE__, __func__, int(qr));
      abort();      
  }
  fprintf(out, "\n");
}



/**
  The function converts selection of points with the required quantity output to the regular selection
  in the cases where the selction is given by the property number. It is used in preprocessor especially.
   
  @param[in] top - pointer to the general topology of the problem solved
  
  @return The function does not return anything but it changes the selid data members.

  Created by TKO, 09.2023
*/
void outquantm::conv_sel_prop(siftop *top)
{
  if (selid.st == sel_prop){
    if (tpnt == atnode)
      selid.conv_selprop(top, gnod);
    else
      selid.conv_selprop(top, gelem);
  }
}



/**
  The function prints the reuired quantity to the output file given by argument out.
  
  @param[in] lcid - load case id.
  @param[in] time - actual time or load coefficient
  @param[in] stepid - time/load step identfier
  @param[in] lcsarray - array of the local coordinate system definitions
  @param[in] nv - actual nodal values of quantities related to nodes only
  @param[in] fmt - format of the output result file.

  @return The function does not return anything but it makes output to the file given by argument out.

  Created by TKo, 09.2023
*/
void outquantm::print_quant(long lcid, double time, long stepid, lcoordsys *lcsarray, gnodvalvm &nv, resfilefmt fmt)
{ 
  collect_quant_val(lcid, lcsarray, nv);
  
  switch (fmt){
    case resfmt_open_dx:
    case resfmt_femcad:
    case resfmt_gid:
    case resfmt_vtk:
    case resfmt_plain_out:
    case resfmt_diag_dat:
      //print_diag_dat(outf, selid, seleip, tnpnt, selip, ncmpr, 
      break;
    default:
      print_err("unknown format of the output file (%d) is required.\n", __FILE__ , __LINE__, __func__, int(fmt));
      abort();
  }
}



/**
  The function collects quantity values in arrays qvv or qvs according to the natural quantity 
  representation. The transformations or scaling of the selected quantity is being performed 
  according to the appropriate options.
  
  @param[in] lcid - load case id
  @param[in] lcsarray - array of the local coordinate system defintions
  
  @return The function does not return anything but changes values either in array qvs or qvv.

  Created by Tomas Koudelka, 09.2023
*/
void outquantm::collect_quant_val(long lcid, lcoordsys *lcsarray, gnodvalvm &nv)
{
  vector auxv;
  matrix auxm;

  switch(qr){
    case scalq:
      clean_qvs();
      break;
    case vectq:
      clean_qvv();
      break;
    case tensq:
      clean_qvv();
      break;
    default:
      print_err("unknown cahracter type %d of the required quantity", __FILE__, __LINE__, __func__, int(qr));
      abort();
  }
  // collect values from nodes or elements
  switch (tpnt){
    case atnode:
      collect_nodal_quant_val(lcid, nv);
      break;
    case atip:
      collect_elem_quant_val(lcid);
      break;
    default:
      print_err("unknown point type %d is required", __FILE__, __LINE__, __func__, int(tpnt));
      abort();
  }

  // make transformation to the required coordinate system
  make_coord_transf(lcsarray);

  
  // make quantity scaling or other transformation of the quantity values
  make_scale_transf();
}



/**
  The function collects nodal quantity values in arrays qvv or qvs according to the natural 
  quantity representation.
  
  @param[in] lcid - load case id
  
  @return The function does not return anything but changes values either in array qvs or qvv.

  Created by Tomas Koudelka, 09.2023
*/
void outquantm::collect_nodal_quant_val(long lcid, gnodvalvm &nv)
{
  long i;
  vector auxv;
  matrix auxm;

  if ((selid.st == sel_list) || (selid.st == sel_single)){ // list of nodal indices or single node are given
    for(i=0; i<selid.n; i++){
      switch(qr){
        case scalq:
          Mt->give_nodal_quant(selid.id1[i], mqn, reftensq, lcid, qvs[i]);
          break;
        case vectq:
          makerefv(auxv, qvv[selid.id1[i]], maxncmp);
          Mt->give_nodal_quant(selid.id1[i], mqn, lcid, seliq, nv, auxv);
          break;
        case tensq:
          makerefm(auxm, qvv[selid.id1[i]], 3, 3);
          Mt->give_nodal_quant(selid.id1[i], mqn, reftensq, lcid, auxm);
          break;
        default:
          print_err("unknown type %d of the required quantity", __FILE__, __LINE__, __func__, int(qr));
          abort();
      }
    }
  }
  else{ // all nodes or range of nodal indices are given
    for(i=0; i<Mt->nn; i++){
      if(selid.presence_id(i)){
        switch(qr){
          case scalq:
            Mt->give_nodal_quant(i, mqn, reftensq, lcid, qvs[i]);
            break;
          case vectq:
            makerefv(auxv, qvv[i], maxncmp);
            Mt->give_nodal_quant(i, mqn, lcid, seliq, nv, auxv);
            break;
          case tensq:
            makerefm(auxm, qvv[i], 3, 3);
            Mt->give_nodal_quant(i, mqn, reftensq, lcid, auxm);
            break;
          default:
            print_err("unknown type %d of the required quantity", __FILE__, __LINE__, __func__, int(qr));
          abort();
        }
      }
    }
  }
}



/**
  The function collects element quantity values in arrays qvv or qvs according to the natural 
  quantity representation.
  
  @param[in] compid - component id
  
  @return The function does not return anything but changes values either in array qvs or qvv.

  Created by Tomas Koudelka, 09.2023
*/
void outquantm::collect_elem_quant_val(long lcid)
{
  long i, j, ipp, tnip;
  vector auxv;
  matrix auxm;

  if ((selid.st == sel_list) || (selid.st == sel_single)){ // list of element indices is given
    for(i=0; i<selid.n; i++){
      ipp=Mt->elements[selid.id1[i]].ipp[0][0];
      tnip = Mt->give_tnip(selid.id1[i]);
      switch(qr){
        case scalq:
          for(j=0; j<tnip; j++, ipp++){
            if (seleip.presence_id(j))
              Mm->give_quant(ipp, mqn, reftensq, lcid, seliq, other_tens_strastre, qvs[ipp]);
          }
          break;
        case vectq:
          for(j=0; j<tnip; j++, ipp++){
            if (seleip.presence_id(j)){
              makerefv(auxv, qvv[ipp], maxncmp);
              Mm->give_quant(i, mqn, lcid, seliq, auxv);
            }
          }
          break;
        case tensq:
          for(j=0; j<tnip; j++, ipp++){
            if (seleip.presence_id(j)){
              makerefm(auxm, qvv[ipp], 3, 3);
              Mm->give_quant(ipp, mqn, reftensq, lcid, seliq, other_tens_strastre, auxm);
            }
          }
          break;
        default:
          print_err("unknown type %d of the required quantity", __FILE__, __LINE__, __func__, int(qr));
          abort();
      }
    }
  }
  else{ // all elements or index ranges are selected
    for(i=0; i<Mt->ne; i++){
      if(selid.presence_id(i)){
        ipp=Mt->elements[i].ipp[0][0];
        tnip = Mt->give_tnip(i);
        switch(qr){
          case scalq:
            for(j=0; j<tnip; j++, ipp++){
              if (seleip.presence_id(j))
                Mm->give_quant(ipp, mqn, reftensq, lcid, seliq, other_tens_strastre, qvs[ipp]);
            }
            break;
          case vectq:
            for(j=0; j<tnip; j++, ipp++){
              if (seleip.presence_id(j)){
                makerefv(auxv, qvv[ipp], maxncmp);
                Mm->give_quant(i, mqn, lcid, seliq, auxv);
              }
            }
            break;
          case tensq:
            for(j=0; j<tnip; j++, ipp++){
              if (seleip.presence_id(j)){
                makerefm(auxm, qvv[ipp], 3, 3);
                Mm->give_quant(ipp, mqn, reftensq, lcid, seliq, other_tens_strastre, auxm);
              }
            }
            break;
          default:
            print_err("unknown type %d of the required quantity", __FILE__, __LINE__, __func__, int(qr));
            abort();
        }
      }
    }
  }
}



/**
  The function transforms %vector or tensor quantities stored in auxiliary array qvv
  to the coordinate system defined by lcs.

  @param[in] pid - point id, where the transformation %matrix will be assmebled
  @param[in] lcsarray - array of the local coordinate system definitions

  @return The function does not return anything, but changes values in the array qvv.
*/
void outquantm::make_coord_transf(lcoordsys *lcsarray)
{
  long i, j, tmdim;
  matrix auxm, tmat;
  vector auxv;
  vector pc(ASTCKVEC(3));

  
  // make transformation to the local coordinate system
  if (lcs > 0){
    tmdim = lcsarray[lcs-1].dim;
    reallocm(RSTCKMAT(tmdim, tmdim, tmat));
    if (qr == vectq){
      for(i=0; i<tnpnt; i++){
        give_pnt_coord(i, pc);
        lcsarray[lcs-1].give_transfmat(tmat, pc, Mp->time);
        makerefv(auxv, qvv[i], ncmpr);
        glvectortransfblock(auxv, tmat);
      }
    }
    if (qr == tensq){
      for(i=0; i<tnpnt; i++){
        give_pnt_coord(i, pc);
        lcsarray[lcs-1].give_transfmat(tmat, pc, Mp->time);
        makerefm(auxm, qvv[i], 3, 3);
        glmatrixtransfblock(auxm, tmat);        
      }
    }
  }

  // make transformation of the nodal quantitiy vector values to the global coordinate system
  // if local coordinate system is defined at nodes
  if ((lcs == 0) && (qr == vectq) && (tpnt == atnode)){    
    for (i=0; i<tnpnt; i++){
      tmdim = Mt->nodes[i].transf;
      if (tmdim > 0){
        give_pnt_coord(i, pc);
        lcsarray[lcs-1].give_transfmat(tmat, pc, Mp->time);
        if ((tmdim == 2) && (Gtm->gnodes[i].ndofn == 3)){
          // 2D beam problem
          tmdim = 3;
          tmat(2,2) = 1.0;
        }
        reallocm(RSTCKMAT(tmdim, tmdim, tmat));
        switch(mqn){
          case displ_q:        // displacement vector
          case react_q:        // reactions vector at nodes     
          case load_vect_q:    // load vector at nodes
          case int_force_q:    // internal force vector
          case resid_vect_q:   // residual vector
            makerefv(auxv, qvv[i], ncmpr);
            lgvectortransfblock(auxv, tmat);
            break;
          default:
            break;
        }
      }
    }
  }
  
  // make transformation to the principal coordinate system
  if ((lcs < 0) && (qr == tensq)){
    long ret = 1;
    matrix tv, pvect(ASTCKMAT(3,3));
    vector auxv, pval(ASTCKVEC(3));
    strastrestate ssst;
    ivector id(ASTCKIVEC(3));

    for(i=0; i<tnpnt; i++){
      ssst = give_pnt_ssst(i);
      switch (tsn){
        case voigtred:
          makerefv(auxv, qvv[i], ncmpr);
          reallocm(RSTCKMAT(3, 3, tv));
          vector_tensor(auxv, tv, ssst, tti);
          break;
        default:
          break;
      }
      if (ncmpr <= 6){
        makerefv(auxv, qvv[i], ncmpr);
        reallocm(RSTCKMAT(3, 3, tv));
        vector_tensor(auxv, tv, ssst, tti);
        ret = 0;
      }
      if (ncmpr == 9){
        makerefm(tv, qvv[i], 3, 3);
        ret = 0;
      }
      if (ret){
        print_err("cannot assemble the second order tensor from %ld component array",
                  __FILE__, __LINE__, __func__, ncmpr);
        abort();
      }
      princ_val(tv, pval, pvect, ni_jac, err_jac, zero, Mp->dim, 1);
      give_all_normal_indices(ssst, id);
      memset(qvv[i], 0, sizeof(*qvv[i]*ncmpr));
      for (j=0; j<3; j++){
        if (id(j) >= 0)
          qvv[i][id(j)] = pval(j);
      }
    }
  }
}



/**
  The function performs additional scaling/transformation of the qunatities in the
  auxiliary arrays qvs or qvv according to the given quantity representation. It can be used 
  for the unit conversion as well as for the other more complex transformations. The 
  scaling/transformations is made with the help of general function qst. 

  If the qst is a constant function then the simple qunatity scaling is performed.
  For another (complex) tansformations a parser function type can be used. The parser 
  expression can refer to point coordinates denoted as x, y, z, time denoted as t 
  and quantity value denoted as q. For vector/tensor qunatities, the same transformation
  formula is applied on all components.

  @return The function does not return anything, but it changes values in arrays qvs and qvv.

  Created by Tomas Koudelka, 09.2023
*/
void outquantm::make_scale_transf()
{
  long i, j, nc;
  
  if (qstf == yes){
    vector pc(ASTCKVEC(5));
    double *val;
    for(i=0; i<tnpnt; i++){
      give_pnt_coord(i, pc);
      pc(3) = Mp->time;
      if (qr == scalq){
        nc = 1;
        val = qvs+i;
      }
      else{
        nc = ncmpr;
        val = qvv[i];
      }
      for(j=0; j<nc; j++){
        if (qst.tfunc == constant){
          // the constant is considered to be a qunatity scaling factor
          val[j] *= qst.getval(Mp->time);
        }
        else{
          // other (complex) transformation is made with the help of parser,
          // e.g. temperature conversion from K to deg C can be made with the help
          // of the following parser expression 'q-273.15'
          pc(4) = val[j];
          val[j] = qst.getval(pc, namepar);
        }
      }
    }
  }
}



/**
  The function returns position vector of the given pid-th point.

  @param[in] pid - point id, i.e. nodal or element ip identifier
  @param[out] pc - resulting position %vector, it will contain point coordinates x, y, and z on the output

  @return The function returns position %vector in the argument pc.

  Created by Tomas Koudelka, 09.2023
*/
void outquantm::give_pnt_coord(long pid, vector &pc)
{
  nullv(pc);
  
  if (tpnt == atnode)  Mt->give_nodal_coord(pid, pc);
  if (tpnt == atip)    ipcoord(Mm->elip[pid], pid, 0, 0, pc);  
}



/**
  The function returns stress/strain state at the given pid-th point.

  @param[in] pid - point id, i.e. nodal or element ip identifier

  @return The function returns stress/strain state of the given point.

  Created by Tomas Koudelka, 09.2023
*/
strastrestate outquantm::give_pnt_ssst(long pid)
{
  strastrestate ssst=strastrestate(0);
  
  if (tpnt == atnode)  ssst = guess_ssst(Mt->nodes[pid].ncompstr);
  if (tpnt == atip)    ssst = Mm->ip[pid].ssst;

  return ssst;
}



/**
  The function sets all values in the auxiliary array qs to zero.

  @return The function does not return anything but it zeroes array qs.

  Created by Tomas Koudelka, 09.2023
*/
void outquantm::clean_qvs()
{
  memset(qvs, 0, sizeof(*qvs)*tnpnt);
}



/**
  The function sets all component values in the auxiliary array qvv to zero.

  @return The function does not return anything but it zeroes array qvv.

  Created by Tomas Koudelka, 09.2023
*/
void outquantm::clean_qvv()
{
  long i;
  
  for (i=0; i<tnpnt; i++){
    memset(qvv[i], 0, sizeof(*qvv[i])*maxncmp);
  }
}


/**
 The function creates string with default label of the required quantity.
  
  @param[in] mq - required mechanical quantity identifier
  @param[in] rmq - referenced mechanical qunatity for tensor calculations, otherwise it is ignored
  @param[in] maxll - maximum length of label
  @param[out] label - output string with default label of the required quantity, 
                      it must be allocated to the maxll length at least  
  
  @return The function returns the default quantity label in the argument label.

  Created by Tomas Koudelka, 10.2023
*/
void outquantm::give_quant_deflabel(mechquant mq, mechquant rmq, long maxll, char *label)
{
  long l;
  
  label[0] = '\0';
  switch (mq){
    case displ_q:        // displacements vector (r)      
      snprintf(label, maxll, "r");
      break;
    case strain_q:       // strain tensor
      snprintf(label, maxll, "eps");
      break;
    case stress_q:       // stress tensor (sig)
      snprintf(label, maxll, "sig");
      break;
    case other_q:          // other array values as a vector (other)
      snprintf(label, maxll, "other");
      break;
    case react_q:        // reactions vector at nodes  (R)
      snprintf(label, maxll, "R");
      break;
    case load_vect_q:    // load vector at nodes (f_l)
      snprintf(label, maxll, "f_l");
      break;
    case int_force_q:    // internal force vector (f_i)
      snprintf(label, maxll, "f_i");
      break;
    case resid_vect_q:   // residual vector (f_r)
      snprintf(label, maxll, "f_r");
      break;
    case tempr_strain_q: // temperature strain tensor (epst)
      snprintf(label, maxll, "epst");
      break;
    case eig_strain_q:   // eigenstrain tensor (eps0)
      snprintf(label, maxll, "eps0");
      break;
    case eig_stress_q:   // eigenstress tensor (sig0)
      snprintf(label, maxll, "sig0");
      break;
    case time_q:         // time (constant scalar quantity at whole domain)
      snprintf(label, maxll, "t");
      break;
    case step_id_q:      // step id (constant scalar quantity at whole domain)
      snprintf(label, maxll, "step_id");
      break;
    case load_fact_q:    // load factor in nonlinear statics problem type (constant scalar quantity at whole domain)
      snprintf(label, maxll, "lambda");
      break;
    case eigval_q:       // eigen values in eigenvalue problem type (constant scalar quantity at whole domain)
      snprintf(label, maxll, "lambda");
      break;
    case macrostrain_q:  // macrostrain (constant tensor quantity at whole domain) (Sig)
      snprintf(label, maxll, "Eps");
      break;
    case macrostress_q:  // macrostress (constant tensor quantity at whole domain)  (Eps)
      snprintf(label, maxll, "Sig");
      break;
    case nonmech_q:      // nonmechanical quantities at ip
      snprintf(label, maxll, "nonmechq");
      break;
    case first_inv:       // A11 + A22 + A33 (I_1)
      snprintf(label, maxll, "I_1(");
      give_quant_deflabel(rmq, nomech_q, maxll, label+4);
      l = strlen(label);
      snprintf(label+l, maxll-l, ")");
      break;
    case second_inv:      // A11*A22 + A11*A33 + A22*A33 - A12^2 - A13^2 - A23^2 (I_2)
      snprintf(label, maxll, "I_2(");
      give_quant_deflabel(rmq, nomech_q, maxll, label+4);
      l = strlen(label);
      snprintf(label+l, maxll-l, ")");
      break;
    case third_inv:       // det|A| (I_3)
      snprintf(label, maxll, "I_3(");
      give_quant_deflabel(rmq, nomech_q, maxll, label+4);
      l = strlen(label);
      snprintf(label+l, maxll-l, ")");
      break;
    case tensor_norm:     // ||A|| = sqrt(a_ij*a_ij) (||.||)
      snprintf(label, maxll, "||");
      give_quant_deflabel(rmq, nomech_q, maxll, label+2);
      l = strlen(label);
      snprintf(label+l, maxll-l, "||");
      break;
    case tensdeviator:    // D_ij = A_ij - delta_ij*A_kk/3 (d)
      snprintf(label, maxll, "d(");
      give_quant_deflabel(rmq, nomech_q, maxll, label+2);
      l = strlen(label);
      snprintf(label+l, maxll-l, ")");
      break;
    case strain_vol:      // eps_v = eps_x + eps_y + eps_z (eps_v)
      snprintf(label, maxll, "eps_v");
      break;
    case mean_stress:     // sig_m = (sig_x + sig_y + sig_z)/3 (sig_m)
      snprintf(label, maxll, "sig_m");
      break;
    case j2inv:           // negative value of the second invariant of stress deviator, i.e. J2 = 1/2 s_ij s_ij (J2)
      snprintf(label, maxll, "J_2s");
      break;
    case von_mises_stress:// sig_eff = sqrt(3*J2)  (sig_vM)
      snprintf(label, maxll, "sig_vM");
      break;
    case strain_deviator: // e_ij = eps_ij - delta_ij*eps_v/3  (d(eps))
      snprintf(label, maxll, "d(eps)");
      break;
    case stress_deviator: // s_ij = sig_ij - delta_ij*sig_m    (d(sig))
      snprintf(label, maxll, "d(sig)");
      break;
    case strain_pl:       // plastic strain tensor              (epsp)
      snprintf(label, maxll, "epsp");
      break;
    case cons_param:      // consistency parameter gamma
      snprintf(label, maxll, "gamma"); 
      break;
    case damage_scal:     // scalar damage (omega)
      snprintf(label, maxll, "omega");
      break;
    case damaget_scal:   // scalar damage (omega_t) in tension
      snprintf(label, maxll, "omega_t");
      break;
    case damagec_scal:   // scalar damage (omega_c) in compression
      snprintf(label, maxll, "omega_c");
      break;
    case damage_tens:    // damage tensor (Omega)
      snprintf(label, maxll, "Omega");
      break;
    case damaget_tens:   // damage tensor (Omega_t) in tension
      snprintf(label, maxll, "Omega_t");
      break;
    case damagec_tens:   // damage tensor (Omega_c) in compression
      snprintf(label, maxll, "Omega_c");
    default:
      snprintf(label, maxll, "%s", mechquant_kwdset.get_str(mq));
  }
  return;
}



/**
  The function returns pointer to the string with default label of the required quantity component.

  @param[in] mq - required mechanical quantity identifier
  @param[in] rmq - referenced mechanical qunatity for tensor calculations, otherwise it is ignored
  @param[in] ssst - stress/strain state of the tensor quantity, not used for other representations (scalar/vector) 
                    of teh required mq.
  @param[in] compid - required component id
  @param[in] maxll - maximum length of label
  @param[out] label - output string with default label of the required quantity, 
                      it must be allocated to the maxll length at least  
  
  
  @return The function returns the default quantity component label in the argument label.

  Created by Tomas Koudelka, 10.2023
*/
void outquantm::give_quant_comp_deflabel(mechquant mq, mechquant rmq, strastrestate ssst, long compid, long maxll, char *label)
{
  long l;
  label[0] = '\0';
  give_quant_deflabel(mq, rmq, maxll, label);
  switch(mq){
    case displ_q:        // displacements vector (r)      
      snprintf(label, maxll, "r_%s", vector_cnamem::displ_cmpstr(ssst, compid));      
      break;
    case strain_q:       // strain tensor
      snprintf(label, maxll, "%s", tensor_cname::strain_cmpstr(ssst, compid));
      break;
    case stress_q:       // stress tensor (sig)
      snprintf(label, maxll, "%s", tensor_cname::stress_cmpstr(ssst, compid));
      break;
    case other_q:          // other array values as a vector (other)
      snprintf(label, maxll, "other_%ld", compid+1);
      break;
    case react_q:        // reactions vector at nodes  (R)
      snprintf(label, maxll, "%s", vector_cnamem::react_cmpstr(ssst, compid));
      break;
    case load_vect_q:    // load vector at nodes (f_l)
      snprintf(label, maxll, "f_l_%s", vector_cnamem::vect_indstr(ssst, compid));
      break;
    case int_force_q:    // internal force vector (f_i)
      snprintf(label, maxll, "f_i_%s", vector_cnamem::vect_indstr(ssst, compid));
      break;
    case resid_vect_q:   // residual vector (f_r)
      snprintf(label, maxll, "f_r_%s", vector_cnamem::vect_indstr(ssst, compid));
      break;
    case tempr_strain_q: // temperature strain tensor (epst)
      snprintf(label, maxll, "epst_%s", tensor_cname::tens_indstr(ssst, compid));
      break;
    case eig_strain_q:   // eigenstrain tensor (eps0)
      snprintf(label, maxll, "eps0_%s", tensor_cname::tens_indstr(ssst, compid));
      break;
    case eig_stress_q:   // eigenstress tensor (sig0)
      snprintf(label, maxll, "sig0_%s", tensor_cname::tens_indstr(ssst, compid));
      break;
    case macrostrain_q:  // macrostrain (constant tensor quantity at whole domain) (Sig)
      snprintf(label, maxll, "Eps_%s", tensor_cname::tens_indstr(ssst, compid));
      break;
    case macrostress_q:  // macrostress (constant tensor quantity at whole domain)  (Eps)
      snprintf(label, maxll, "Sig_%s", tensor_cname::tens_indstr(ssst, compid));
      break;
    case nonmech_q:      // nonmechanical quantities at ip
      snprintf(label, maxll, "nonmechq_%ld", compid+1);
      break;
    case tensdeviator:    // D_ij = A_ij - delta_ij*A_kk/3 (d)
      snprintf(label, maxll, "d(");
      give_quant_deflabel(rmq, nomech_q, maxll-2, label+2);
      l = strlen(label);
      snprintf(label+l, maxll-l, ")_%s", tensor_cname::tens_indstr(ssst, compid));
      break;
    case strain_deviator: // e_ij = eps_ij - delta_ij*eps_v/3  (d(eps))
      snprintf(label, maxll, "d(eps)_%s", tensor_cname::tens_indstr(ssst, compid));
      break;
    case stress_deviator: // s_ij = sig_ij - delta_ij*sig_m    (d(sig))
      snprintf(label, maxll, "d(sig)_%s", tensor_cname::tens_indstr(ssst, compid));
      break;
    case strain_pl:       // plastic strain tensor              (epsp)
      snprintf(label, maxll, "epsp_%s", tensor_cname::tens_indstr(ssst, compid));
      break;
    case cons_param:      // consistency parameter gamma
      snprintf(label, maxll, "gamma"); 
      break;
    case damage_scal:     // scalar damage (omega)
      snprintf(label, maxll, "omega");
      break;
    case damaget_scal:   // scalar damage (omega_t) in tension
      snprintf(label, maxll, "omega_t");
      break;
    case damagec_scal:   // scalar damage (omega_c) in compression
      snprintf(label, maxll, "omega_c");
      break;
    case damage_tens:    // damage tensor (Omega)
      snprintf(label, maxll, "Omega_%s", tensor_cname::tens_indstr(ssst, compid));
      break;
    case damaget_tens:   // damage tensor (Omega_t) in tension
      snprintf(label, maxll, "Omega_t_%s", tensor_cname::tens_indstr(ssst, compid));
      break;
    case damagec_tens:   // damage tensor (Omega_c) in compression
      snprintf(label, maxll, "Omega_c_%s", tensor_cname::tens_indstr(ssst, compid));
    default:
      snprintf(label, maxll, "%s", mechquant_kwdset.get_str(mq));
  }
  return;
}
