#ifndef ABC_H
#define ABC_H

void gennodes(long neh, long nev, snode *nodes, long *non, snode *uzly, long &j, long **edges,long ***edgenod,long k, selement &el);
void searchsimedg(selement *elements, long i, long **edges);
void genelements (long neh, long nev, selement &el, selement *elements, long **edgenod, long **edges, long l, long &k);
void copy_similar_edgenod (long ***edgenod, long **edges, long a, long k, long b);
void read_edgdiv (XFILE *in, long *&hrany, long &a);
void div_elements (long *edgdiv, long **elemdiv, long *propedg, long **edges, long k);


#endif
