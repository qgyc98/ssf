#ifndef AGGREGATOR_H
#define AGGREGATOR_H

class aggregator
{
 public:
  
  //   JK -> PM

  //  defines number of aggregates
  void define_na (long i);
  //  defines matrix size
  //  i - number of unknowns (n)
  //  j - number of nonzero entries in matrix (nz)
  void define_matsize (long i,long j);
  //  assembles matrix data
  //  i - array of row indices
  //  j - array of column indices
  //  v - array of matrix entries
  void matrix_assemble (long *i,long *j,double *v);
  //  definice stupne finalniho prolongatoru
  //  i - pozadovany stupen zhlazeneho prolongatoru
  void define_degree (long i);

  
  //  PM -> JK

  //  assembles list of numbers of nodes in aggregates
  //  pocet uzlu v agregatech
  //void assemble_lnnagr (long *i);
  void assemble_lnnagr ();
  //  assembles node numbers in aggregates
  //  seznamy uzlu v jednotlivych agregatech
  //  lnagr[i][j] = k - j-ty uzel v i-tem agregatu ma cislo k
  //void assemble_lnagr (long **lnagr);
  void assemble_lnagr ();
  
  //  sestavuje pocet nenulovych prvku ve sloupcich tentativniho prolongatoru
  //  i[j]=k - v j-tem sloupci je k nenulovych prvku
  void assemble_lnuagr_tp (long *i);

  //  sestavuje tentative prolongator
  //  agrid - cislo pozadovaneho agregatu
  //  cid - cislo pozadovaneho sloupce v agregatu
  //  i - pole indexu
  //  a - pole hodnot (v pripade rigid body motions v mechanice)
  void assemble_tentative_prol (long agrid,long cid,long *i,double *a);
  
  //  vraci skutecny stupen prolongatoru
  long give_deg ();

  //  sestavuje pocet nenulovych prvku ve sloupcich smoothed prolongatoru
  //  i[j]=k - v j-tem sloupci je k nenulovych prvku
  void assemble_lnuagr_sp (long *i);

  //  uprava, sobota 6.1.2007
  //  agrid - cislo pozadovaneho agregatu
  //  i obsahuje pole dotcenych uzly po zhlazeni na agregatu agrid
  void assemble_lnuagr_sp (long agrid,long *i);
  //  konec upravy


  //  sestavuje smoothed prolongator
  //  agrid - cislo pozadovaneho agregatu
  //  cid - cislo pozadovaneho sloupce v agregatu
  //  i - pole indexu
  //  a - pole hodnot (v pripade rigid body motions v mechanice)
  void assemble_smoothed_prol (long agrid,long cid,long *i,double *a);
  
  
  //  number of aggregates
  long na;

  //  list of numbers of nodes in aggregates
  long *lnnagr;
  //  list of node numbers in aggregates
  long **lnagr;
  
  //  number of nonzero entries in matrix
  long nz;
  //  number of unknowns
  long n;
  //  row indices
  long *ri;
  //  column indices
  long *ci;
  //  matrix entries
  double *val;

};

#endif
