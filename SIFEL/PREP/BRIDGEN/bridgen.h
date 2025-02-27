#ifndef BRIDGEN_H
#define BRIDGEN_H

struct section        // struktura usek, nesouci informace o krajnich rezech, vlastni delce a delelni
{
  long   a;           // cislo prvniho rezu v useku
  long   b;           // cislo druheho rezu v useku
  double length;      // delka useku
  long   division;    // pocet rezu v useku
  long   quadlin;     // urcuje zda ma usek linearni nebo kvadraticky prubeh; dva pro kvadraticky, jedna pro linearni
  long   prop;        // cislo objemove property generovanych 3D prvku
  section(){a=b=division=0;length=0.0;quadlin=0;prop=-1;};
};

void gennodes(long neh, long nev, snode *nodes, long *non, snode *uzly, long &j, long **edges, long **corner,long ***edgenod, long **elemdiv, long k);
void searchsimedg(selement *&elements, long i, long **edges, long **corner);
void genelements (long neh, long nev, long prope, long *propedg, long *propsurf, selement *elements, long **edgenod, long **edges, long **corner, long l, long &k);
void copy_similar_edgenod (long ***edgenod, long **edges, long a, long k, long b);
void read (XFILE *in, long *&edgdiv, section *&arsec, siftop *&cuts, long &numcuts, siftop *&topcuts, long &numsec);
void div_elements (long *edgdiv, long **elemdiv, long *propedg, long **edges, long k);
void gensection (section &sec, siftop *&topcuts, siftop &top, long &numel, long &numnod, double &zetko);
void gensection2 (long pcut_id, section &sec, siftop *&topcuts, siftop &top, long &numel, long &numnod, double &zetko);
void alloc (long **&edges, long ***&elemdiv, long numcuts, siftop *&cuts, long **&corner);
void dealloc (long **edges, long ***elemdiv, long numcuts, siftop *cuts, long **corner);
void copy_section_nodes2top(siftop &sect, long ini_nodesc_id, long num_nodes, long &anid, siftop &top);
void copy_shift_set_section_nodes2top(siftop &sect, long ini_nodesc_id, long num_nodes, double dx, double dy, double dz,
                                      long &anid, siftop &top);
#endif
