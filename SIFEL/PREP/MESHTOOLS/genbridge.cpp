#include <stdio.h>
#include <string.h>
#include <stdlib.h>


int main(int argc, char* argv[])
{
  FILE *top,*out;
  long i,j;
  // pocet zadanych rezu
  long ncut;
  //cislo vertexu
  long ver;
  // x souradnice
  double **x;
  // y souradnice
  double **y;
  // z souradnice
  double **z;
  // typ pricneho rezu
  long *typecut;
  // deleni pricneho rezu na prvky
  long *cscount;
  // deleni podelne
  long *lcount;
  
  cscount = new long[9];
  top = fopen(argv[1],"r");
  // cteni deleni krivek
  // pocet pro tloustku desky
  fscanf(top,"%ld",&cscount[0]);
  // konzola 1
  fscanf(top,"%ld",&cscount[1]);
  // konzola 2
  fscanf(top,"%ld",&cscount[2]);
  // konzola 3
  fscanf(top,"%ld",&cscount[3]);
  // konzola 4
  fscanf(top,"%ld",&cscount[4]);
  // prostredek komurky
  fscanf(top,"%ld",&cscount[5]);
  // tloustka komurky
  fscanf(top,"%ld",&cscount[6]);
  // bok komurky
  fscanf(top,"%ld",&cscount[7]);
  // ostatni v komurce
  fscanf(top,"%ld",&cscount[8]);
  // cteni poctu rezu
  fscanf(top,"%ld",&ncut);
  typecut = new long[ncut];
  lcount = new long[ncut];
  x = new double*[ncut];
  y = new double*[ncut];
  z = new double*[ncut];
  // cteni jednotlivych rezu
  for(i = 0; i < ncut; i++){
    x[i] = new double[32];
    y[i] = new double[32];
    z[i] = new double[32];
    if(i == 0){
      fscanf(top,"%ld",&typecut[i]);
    }
    else{
      fscanf(top,"%ld %ld",&typecut[i],&lcount[i]);
    }
    for(j = 0; j < 32; j++){
      fscanf(top,"%lf",&x[i][j]);
      fscanf(top,"%lf",&y[i][j]);
      fscanf(top,"%lf",&z[i][j]);
    }
  }
  fclose(top);
  out = fopen(argv[2],"w");
  for(i = 0; i < ncut; i++){
    if(typecut[i] == 1){
      if(i == 0){
	// prvni rez
	fprintf(out,"# %ld rez - typ %ld\n",i+1,typecut[i]);
	for(j = 0; j < 32; j++){
	  ver=j+1;
	  fprintf(out,"vertex %ld xyz %lf %lf %lf\n",ver,x[0][j],y[0][j],z[0][j]);
	}
	fprintf(out,"curve  1 vertex  1   2 count %ld\n",cscount[0]);
	fprintf(out,"curve  2 vertex  2   3 count %ld\n",cscount[1]);
	fprintf(out,"curve  3 vertex  3   4 count %ld\n",cscount[2]);
	fprintf(out,"curve  4 vertex  4   5 count %ld\n",cscount[7]);
	fprintf(out,"curve  5 vertex  5   6 count %ld\n",cscount[7]); 
	fprintf(out,"curve  6 vertex  6   7 count %ld\n",cscount[8]);
	fprintf(out,"curve  7 vertex  7   8 count %ld\n",cscount[5]);
	fprintf(out,"curve  8 vertex  8   9 count %ld\n",cscount[8]);
	fprintf(out,"curve  9 vertex  9  10 count %ld\n",cscount[7]);
	fprintf(out,"curve 10 vertex 10 11 count %ld\n",cscount[7]);
	fprintf(out,"curve 11 vertex 11 12 count %ld\n",cscount[3]);
	fprintf(out,"curve 12 vertex 12 13 count %ld\n",cscount[4]);
	fprintf(out,"curve 13 vertex 13 14 count  %ld\n",cscount[0]);
	fprintf(out,"curve 14 vertex 14 15 count %ld\n",cscount[4]);
	fprintf(out,"curve 15 vertex 15 16 count %ld\n",cscount[3]);
	fprintf(out,"curve 16 vertex 16 17 count %ld\n",cscount[6]);
	fprintf(out,"curve 17 vertex 17 18 count %ld\n",cscount[8]);
	fprintf(out,"curve 18 vertex 18 19 count %ld\n",cscount[5]);
	fprintf(out,"curve 19 vertex 19 20 count %ld\n",cscount[8]);
	fprintf(out,"curve 20 vertex 20 21 count %ld\n",cscount[6]);
	fprintf(out,"curve 21 vertex 21 22 count %ld\n",cscount[2]);
	fprintf(out,"curve 22 vertex 22  1 count %ld\n",cscount[1]);
	fprintf(out,"curve 23 vertex 23 24 count %ld\n",cscount[8]);
	fprintf(out,"curve 24 vertex 24 25 count %ld\n",cscount[7]);
	fprintf(out,"curve 25 vertex 25 26 count %ld\n",cscount[7]);
	fprintf(out,"curve 26 vertex 26 27 count %ld\n",cscount[8]);
	fprintf(out,"curve 27 vertex 27 28 count %ld\n",cscount[5]);
	fprintf(out,"curve 28 vertex 28 29 count %ld\n",cscount[8]);
	fprintf(out,"curve 29 vertex 29 30 count %ld\n",cscount[7]);
	fprintf(out,"curve 30 vertex 30 31 count %ld\n",cscount[7]);
	fprintf(out,"curve 31 vertex 31 32 count %ld\n",cscount[8]);
	fprintf(out,"curve 32 vertex 32 23 count %ld\n",cscount[5]);
	fprintf(out,"\n");
	fprintf(out,"curve 33 vertex  3 22 count %ld \n",cscount[0]);
	fprintf(out,"curve 34 vertex  4 21 count %ld \n",cscount[0]);
	fprintf(out,"curve 35 vertex 24 20 count %ld \n",cscount[0]);
	fprintf(out,"curve 36 vertex 23 19 count %ld \n",cscount[0]);
	fprintf(out,"curve 37 vertex  6 26 count %ld\n",cscount[6]);
	fprintf(out,"curve 38 vertex  7 27 count %ld\n",cscount[6]);
	fprintf(out,"curve 39 vertex  8 28 count %ld\n",cscount[6]);
	fprintf(out,"curve 40 vertex  9 29 count %ld\n",cscount[6]);
	fprintf(out,"curve 41 vertex 32 18 count %ld \n",cscount[0]);
	fprintf(out,"curve 42 vertex 31 17 count %ld \n",cscount[0]);
	fprintf(out,"curve 43 vertex 11 16 count %ld \n",cscount[0]);
	fprintf(out,"curve 44 vertex 12 15 count %ld \n",cscount[0]);
	fprintf(out,"\n");
	fprintf(out,"curve 45 vertex  4 24 count %ld\n",cscount[6]);
	fprintf(out,"curve 46 vertex 31 11 count %ld\n",cscount[6]);
	fprintf(out,"curve 47 vertex  5 25 count %ld\n",cscount[6]);
	fprintf(out,"curve 48 vertex 30 10 count %ld\n",cscount[6]);
	fprintf(out,"surface 1 curve 1 2 33 22\n");
	fprintf(out,"surface 2 curve 33 3 34 21\n");
	fprintf(out,"surface 3 curve 20 34 45 35\n");
	fprintf(out,"surface 4 curve 45  4 47 24\n");
	fprintf(out,"surface 5 curve 47  5 37 25\n");
	fprintf(out,"surface 6 curve 35 23 36 19\n");
	fprintf(out,"surface 7 curve 37  6 38 26\n");
	fprintf(out,"surface 8 curve 38  7 39 27\n");
	fprintf(out,"surface 9 curve 36 32 41 18\n");
	fprintf(out,"surface 10 curve 39 8 40 28\n");
	fprintf(out,"surface 11 curve 40  9 48 29\n");
	fprintf(out,"surface 12 curve 48 10 46 30\n");
	fprintf(out,"surface 13 curve 41 31 42 17\n");
	fprintf(out,"surface 14 curve 46 43 16 42\n");
	fprintf(out,"surface 15 curve 43 11 44 15\n");
	fprintf(out,"surface 16 curve 44 12 13 14\n");
	fprintf(out,"\n");

      }
      if(i != 0){
	fprintf(out,"# %ld rez - typ %ld\n",i+1,typecut[i]);
	for(j = 0; j < 32; j++){
	  ver=j+1+i*100;
	  fprintf(out,"vertex %ld xyz %lf %lf %lf\n",ver,x[i][j],y[i][j],z[i][j]);
	}
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",1+i*100,1+i*100,2+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  2+i*100,  2+i*100,  3+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  3+i*100,  3+i*100,  4+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  4+i*100,  4+i*100,  5+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  5+i*100,  5+i*100,  6+i*100,cscount[7]); 
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  6+i*100,  6+i*100,  7+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  7+i*100,  7+i*100,  8+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  8+i*100,  8+i*100,  9+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  9+i*100,  9+i*100, 10+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 10+i*100, 10+i*100,11+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 11+i*100, 11+i*100,12+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 12+i*100, 12+i*100,13+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 13+i*100, 13+i*100,14+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 14+i*100, 14+i*100,15+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 15+i*100, 15+i*100,16+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 16+i*100, 16+i*100,17+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 17+i*100, 17+i*100,18+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 18+i*100, 18+i*100,19+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 19+i*100, 19+i*100,20+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 20+i*100, 20+i*100,21+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 21+i*100, 21+i*100,22+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 22+i*100, 22+i*100, 1+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 23+i*100, 23+i*100,24+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 24+i*100, 24+i*100,25+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 25+i*100, 25+i*100,26+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 26+i*100, 26+i*100,27+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 27+i*100, 27+i*100,28+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 28+i*100, 28+i*100,29+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 29+i*100, 29+i*100,30+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 30+i*100, 30+i*100,31+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 31+i*100, 31+i*100,32+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 32+i*100, 32+i*100,23+i*100,cscount[5]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 33+i*100,  3+i*100, 22+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 34+i*100,  4+i*100, 21+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 35+i*100, 24+i*100, 20+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 36+i*100, 23+i*100, 19+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 37+i*100,  6+i*100, 26+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 38+i*100,  7+i*100, 27+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 39+i*100,  8+i*100, 28+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 40+i*100,  9+i*100, 29+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 41+i*100, 32+i*100, 18+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 42+i*100, 31+i*100, 17+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 43+i*100, 11+i*100, 16+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 44+i*100, 12+i*100, 15+i*100,cscount[0]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 45+i*100,  4+i*100, 24+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 46+i*100, 31+i*100, 11+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 47+i*100,  5+i*100, 25+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 48+i*100, 30+i*100, 10+i*100,cscount[6]);
	fprintf(out,"\n");
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  1+i*100, 1+i*100, 2+i*100,33+i*100,22+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  2+i*100,33+i*100, 3+i*100,34+i*100,21+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  3+i*100,20+i*100,34+i*100,45+i*100,35+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  4+i*100,45+i*100, 4+i*100,47+i*100,24+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  5+i*100,47+i*100, 5+i*100,37+i*100,25+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  6+i*100,35+i*100,23+i*100,36+i*100,19+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  7+i*100,37+i*100, 6+i*100,38+i*100,26+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  8+i*100,38+i*100, 7+i*100,39+i*100,27+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  9+i*100,36+i*100,32+i*100,41+i*100,18+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 10+i*100,39+i*100, 8+i*100,40+i*100,28+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 11+i*100,40+i*100, 9+i*100,48+i*100,29+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 12+i*100,48+i*100,10+i*100,46+i*100,30+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 13+i*100,41+i*100,31+i*100,42+i*100,17+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 14+i*100,46+i*100,43+i*100,16+i*100,42+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 15+i*100,43+i*100,11+i*100,44+i*100,15+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 16+i*100,44+i*100,12+i*100,13+i*100,14+i*100);
	fprintf(out,"\n");
	fprintf(out,"# spojeni rezu %ld (typ %ld) a %ld (typ %ld)\n",i,typecut[i-1],i+1,typecut[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10001+i*100,1+(i-1)*100,1+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10002+i*100,2+(i-1)*100,2+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10003+i*100,3+(i-1)*100,3+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10004+i*100,4+(i-1)*100,4+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10005+i*100,5+(i-1)*100,5+i*100,lcount[i]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10006+i*100,6+(i-1)*100,6+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][5]+(x[i][5]-x[i-1][5])/2,y[i-1][5]+(y[i][5]-y[i-1][5])/2,z[i-1][5]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10007+i*100,7+(i-1)*100,7+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][6]+(x[i][6]-x[i-1][6])/2,y[i-1][6]+(y[i][6]-y[i-1][6])/2,z[i-1][6]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10008+i*100,8+(i-1)*100,8+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][7]+(x[i][7]-x[i-1][7])/2,y[i-1][7]+(y[i][7]-y[i-1][7])/2,z[i-1][7]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10009+i*100,9+(i-1)*100,9+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][8]+(x[i][8]-x[i-1][8])/2,y[i-1][8]+(y[i][8]-y[i-1][8])/2,z[i-1][8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10010+i*100,10+(i-1)*100,10+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10011+i*100,11+(i-1)*100,11+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10012+i*100,12+(i-1)*100,12+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10013+i*100,13+(i-1)*100,13+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10014+i*100,14+(i-1)*100,14+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10015+i*100,15+(i-1)*100,15+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10016+i*100,16+(i-1)*100,16+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10017+i*100,17+(i-1)*100,17+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10018+i*100,18+(i-1)*100,18+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10019+i*100,19+(i-1)*100,19+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10020+i*100,20+(i-1)*100,20+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10021+i*100,21+(i-1)*100,21+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10022+i*100,22+(i-1)*100,22+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10023+i*100,23+(i-1)*100,23+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10024+i*100,24+(i-1)*100,24+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10025+i*100,25+(i-1)*100,25+i*100,lcount[i]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10026+i*100,26+(i-1)*100,26+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][25]+(x[i][25]-x[i-1][25])/2,y[i-1][25]+(y[i][25]-y[i-1][25])/2,z[i-1][25]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10027+i*100,27+(i-1)*100,27+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][26]+(x[i][26]-x[i-1][26])/2,y[i-1][26]+(y[i][26]-y[i-1][26])/2,z[i-1][26]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10028+i*100,28+(i-1)*100,28+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][27]+(x[i][27]-x[i-1][27])/2,y[i-1][27]+(y[i][27]-y[i-1][27])/2,z[i-1][27]);      
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10029+i*100,29+(i-1)*100,29+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][28]+(x[i][28]-x[i-1][28])/2,y[i-1][28]+(y[i][28]-y[i-1][28])/2,z[i-1][28]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10030+i*100,30+(i-1)*100,30+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10031+i*100,31+(i-1)*100,31+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10032+i*100,32+(i-1)*100,32+i*100,lcount[i]);
	for(j = 0; j < 32; j++){
	  if(j != 21 && j != 31){
	    fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,j+10002+i*100,j+1+i*100,j+10001+i*100);
	  }
	  else{
	    if(j == 21){
	      fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,10001+i*100,j+1+i*100,j+10001+i*100);
	    }
	    if(j == 31){
	      fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,10023+i*100,j+1+i*100,j+10001+i*100);
	    }
	  }
	}
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10033+i*100,(i-1)*100+33,10003+i*100,i*100+33,10022+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10034+i*100,(i-1)*100+34,10004+i*100,i*100+34,10021+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10035+i*100,(i-1)*100+35,10024+i*100,i*100+35,10020+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10036+i*100,(i-1)*100+36,10023+i*100,i*100+36,10019+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10037+i*100,(i-1)*100+37,10006+i*100,i*100+37,10026+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10038+i*100,(i-1)*100+38,10007+i*100,i*100+38,10027+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10039+i*100,(i-1)*100+39,10008+i*100,i*100+39,10028+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10040+i*100,(i-1)*100+40,10009+i*100,i*100+40,10029+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10041+i*100,(i-1)*100+41,10032+i*100,i*100+41,10018+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10042+i*100,(i-1)*100+42,10031+i*100,i*100+42,10017+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10043+i*100,(i-1)*100+43,10011+i*100,i*100+43,10016+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10044+i*100,(i-1)*100+44,10012+i*100,i*100+44,10015+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10045+i*100,(i-1)*100+45,10004+i*100,i*100+45,10024+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10046+i*100,(i-1)*100+46,10031+i*100,i*100+46,10011+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10047+i*100,(i-1)*100+47,10005+i*100,i*100+47,10025+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10048+i*100,(i-1)*100+48,10030+i*100,i*100+48,10010+i*100);
	fprintf(out,"\n# definice regionu\n");
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",1+i*100,10001+i*100,10002+i*100,10033+i*100,10022+i*100,(i-1)*100+1,i*100+1);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",2+i*100,10033+i*100,10003+i*100,10034+i*100,10021+i*100,(i-1)*100+2,i*100+2);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",3+i*100,10034+i*100,10045+i*100,10035+i*100,10020+i*100, 3+(i-1)*100,3+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",4+i*100,10004+i*100,10047+i*100,10024+i*100,10045+i*100, 4+(i-1)*100,4+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld  hexa map yes equidistant size def\n",5+i*100,10005+i*100,10037+i*100,10025+i*100,10047+i*100,5+(i-1)*100,5+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld  hexa map yes equidistant size def\n",6+i*100,10035+i*100,10023+i*100,10036+i*100,10019+i*100,6+(i-1)*100,6+i*100);
	fprintf(out,"region %ld boundary surface %ld %ld -%ld -%ld %ld -%ld  hexa map yes equidistant size def\n",7+i*100,10038+i*100,10026+i*100,10037+i*100,10006+i*100,7+(i-1)*100,7+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld %ld %ld -%ld hexa map yes equidistant size def\n",8+i*100,10038+i*100,10007+i*100,10039+i*100,10027+i*100,8+(i-1)*100,8+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",9+i*100,10036+i*100,10032+i*100,10041+i*100,10018+i*100,9+(i-1)*100,9+i*100 );
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld %ld %ld -%ld hexa map yes equidistant size def\n",10+i*100,10039+i*100,10008+i*100,10040+i*100,10028+i*100,10+(i-1)*100,10+i*100);
	fprintf(out,"region %ld boundary surface %ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",11+i*100,10029+i*100,10040+i*100,10009+i*100,10048+i*100,11+(i-1)*100,11+i*100);
	fprintf(out,"region %ld boundary surface %ld %ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",12+i*100,10030+i*100,10048+i*100,10010+i*100,10046+i*100,12+(i-1)*100,12+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",13+i*100,10041+i*100,10031+i*100,10042+i*100,10017+i*100,13+(i-1)*100,13+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",14+i*100,10042+i*100,10046+i*100,10043+i*100,10016+i*100,14+(i-1)*100,14+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",15+i*100,10043+i*100,10011+i*100,10044+i*100,10015+i*100,15+(i-1)*100,15+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",16+i*100,10044+i*100,10012+i*100,10013+i*100,10014+i*100,16+(i-1)*100,16+i*100);
      }
    }
    // zacatek ztuzidla
    if(typecut[i] == 2){
      if(i == 0){
	fprintf(out,"# %ld rez - typ %ld\n",i+1,typecut[i]);
	for(j = 0; j < 32; j++){
	  ver=j+1+i*100;
	  fprintf(out,"vertex %ld xyz %lf %lf %lf\n",ver,x[i][j],y[i][j],z[i][j]);
	}
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",1+i*100,1+i*100,2+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  2+i*100,  2+i*100,  3+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  3+i*100,  3+i*100,  4+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  4+i*100,  4+i*100,  5+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  5+i*100,  5+i*100,  6+i*100,cscount[7]); 
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  6+i*100,  6+i*100,  7+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  7+i*100,  7+i*100,  8+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  8+i*100,  8+i*100,  9+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  9+i*100,  9+i*100, 10+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 10+i*100, 10+i*100,11+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 11+i*100, 11+i*100,12+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 12+i*100, 12+i*100,13+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 13+i*100, 13+i*100,14+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 14+i*100, 14+i*100,15+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 15+i*100, 15+i*100,16+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 16+i*100, 16+i*100,17+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 17+i*100, 17+i*100,18+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 18+i*100, 18+i*100,19+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 19+i*100, 19+i*100,20+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 20+i*100, 20+i*100,21+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 21+i*100, 21+i*100,22+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 22+i*100, 22+i*100, 1+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 23+i*100, 23+i*100,24+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 24+i*100, 24+i*100,25+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 25+i*100, 25+i*100,26+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 26+i*100, 26+i*100,27+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 27+i*100, 27+i*100,28+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 28+i*100, 28+i*100,29+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 29+i*100, 29+i*100,30+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 30+i*100, 30+i*100,31+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 31+i*100, 31+i*100,32+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 32+i*100, 32+i*100,23+i*100,cscount[5]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 33+i*100,  3+i*100, 22+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 34+i*100,  4+i*100, 21+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 35+i*100, 24+i*100, 20+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 36+i*100, 23+i*100, 19+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 37+i*100,  6+i*100, 26+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 38+i*100,  7+i*100, 27+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 39+i*100,  8+i*100, 28+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 40+i*100,  9+i*100, 29+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 41+i*100, 32+i*100, 18+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 42+i*100, 31+i*100, 17+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 43+i*100, 11+i*100, 16+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 44+i*100, 12+i*100, 15+i*100,cscount[0]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 45+i*100,  4+i*100, 24+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 46+i*100, 31+i*100, 11+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 47+i*100,  5+i*100, 25+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 48+i*100, 30+i*100, 10+i*100,cscount[6]);
	fprintf(out,"\n");
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  1+i*100, 1+i*100, 2+i*100,33+i*100,22+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  2+i*100,33+i*100, 3+i*100,34+i*100,21+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  3+i*100,20+i*100,34+i*100,45+i*100,35+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  4+i*100,45+i*100, 4+i*100,47+i*100,24+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  5+i*100,47+i*100, 5+i*100,37+i*100,25+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  6+i*100,35+i*100,23+i*100,36+i*100,19+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  7+i*100,37+i*100, 6+i*100,38+i*100,26+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  8+i*100,38+i*100, 7+i*100,39+i*100,27+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  9+i*100,36+i*100,32+i*100,41+i*100,18+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 10+i*100,39+i*100, 8+i*100,40+i*100,28+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 11+i*100,40+i*100, 9+i*100,48+i*100,29+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 12+i*100,48+i*100,10+i*100,46+i*100,30+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 13+i*100,41+i*100,31+i*100,42+i*100,17+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 14+i*100,46+i*100,43+i*100,16+i*100,42+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 15+i*100,43+i*100,11+i*100,44+i*100,15+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 16+i*100,44+i*100,12+i*100,13+i*100,14+i*100);
	fprintf(out,"\n");
	fprintf(out,"vertex %ld xyz %lf %lf %lf\n",33+i*100,x[i][24],y[i][22],z[i][24]);
	fprintf(out,"vertex %ld xyz %lf %lf %lf\n",34+i*100,x[i][29],y[i][31],z[i][29]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",49+i*100,23+i*100,33+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",50+i*100,33+i*100,25+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",51+i*100,33+i*100,27+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",52+i*100,32+i*100,34+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",53+i*100,34+i*100,30+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",54+i*100,34+i*100,28+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",55+i*100,33+i*100,34+i*100,cscount[5]);
	fprintf(out,"\n");
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",17+i*100,23+i*100,24+i*100,50+i*100,49+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",18+i*100,25+i*100,26+i*100,51+i*100,50+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",19+i*100,32+i*100,49+i*100,55+i*100,52+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",20+i*100,55+i*100,51+i*100,27+i*100,54+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",21+i*100,52+i*100,53+i*100,30+i*100,31+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",22+i*100,53+i*100,54+i*100,28+i*100,29+i*100);
	fprintf(out,"\n");
      }
      else{
	if(typecut[i-1]  == 1){
	for(j = 0; j < 32; j++){
	  ver=j+1+i*100;
	  fprintf(out,"vertex %ld xyz %lf %lf %lf\n",ver,x[i][j],y[i][j],z[i][j]);
	}
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",1+i*100,1+i*100,2+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  2+i*100,  2+i*100,  3+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  3+i*100,  3+i*100,  4+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  4+i*100,  4+i*100,  5+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  5+i*100,  5+i*100,  6+i*100,cscount[7]); 
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  6+i*100,  6+i*100,  7+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  7+i*100,  7+i*100,  8+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  8+i*100,  8+i*100,  9+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  9+i*100,  9+i*100, 10+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 10+i*100, 10+i*100,11+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 11+i*100, 11+i*100,12+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 12+i*100, 12+i*100,13+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 13+i*100, 13+i*100,14+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 14+i*100, 14+i*100,15+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 15+i*100, 15+i*100,16+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 16+i*100, 16+i*100,17+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 17+i*100, 17+i*100,18+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 18+i*100, 18+i*100,19+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 19+i*100, 19+i*100,20+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 20+i*100, 20+i*100,21+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 21+i*100, 21+i*100,22+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 22+i*100, 22+i*100, 1+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 23+i*100, 23+i*100,24+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 24+i*100, 24+i*100,25+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 25+i*100, 25+i*100,26+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 26+i*100, 26+i*100,27+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 27+i*100, 27+i*100,28+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 28+i*100, 28+i*100,29+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 29+i*100, 29+i*100,30+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 30+i*100, 30+i*100,31+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 31+i*100, 31+i*100,32+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 32+i*100, 32+i*100,23+i*100,cscount[5]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 33+i*100,  3+i*100, 22+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 34+i*100,  4+i*100, 21+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 35+i*100, 24+i*100, 20+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 36+i*100, 23+i*100, 19+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 37+i*100,  6+i*100, 26+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 38+i*100,  7+i*100, 27+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 39+i*100,  8+i*100, 28+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 40+i*100,  9+i*100, 29+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 41+i*100, 32+i*100, 18+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 42+i*100, 31+i*100, 17+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 43+i*100, 11+i*100, 16+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 44+i*100, 12+i*100, 15+i*100,cscount[0]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 45+i*100,  4+i*100, 24+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 46+i*100, 31+i*100, 11+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 47+i*100,  5+i*100, 25+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 48+i*100, 30+i*100, 10+i*100,cscount[6]);
	fprintf(out,"\n");
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  1+i*100, 1+i*100, 2+i*100,33+i*100,22+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  2+i*100,33+i*100, 3+i*100,34+i*100,21+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  3+i*100,20+i*100,34+i*100,45+i*100,35+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  4+i*100,45+i*100, 4+i*100,47+i*100,24+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  5+i*100,47+i*100, 5+i*100,37+i*100,25+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  6+i*100,35+i*100,23+i*100,36+i*100,19+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  7+i*100,37+i*100, 6+i*100,38+i*100,26+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  8+i*100,38+i*100, 7+i*100,39+i*100,27+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  9+i*100,36+i*100,32+i*100,41+i*100,18+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 10+i*100,39+i*100, 8+i*100,40+i*100,28+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 11+i*100,40+i*100, 9+i*100,48+i*100,29+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 12+i*100,48+i*100,10+i*100,46+i*100,30+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 13+i*100,41+i*100,31+i*100,42+i*100,17+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 14+i*100,46+i*100,43+i*100,16+i*100,42+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 15+i*100,43+i*100,11+i*100,44+i*100,15+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 16+i*100,44+i*100,12+i*100,13+i*100,14+i*100);
	fprintf(out,"\n");
	fprintf(out,"vertex %ld xyz %lf %lf %lf\n",33+i*100,x[i][24],y[i][22],z[i][24]);
	fprintf(out,"vertex %ld xyz %lf %lf %lf\n",34+i*100,x[i][29],y[i][31],z[i][29]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",49+i*100,23+i*100,33+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",50+i*100,33+i*100,25+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",51+i*100,33+i*100,27+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",52+i*100,32+i*100,34+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",53+i*100,34+i*100,30+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",54+i*100,34+i*100,28+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",55+i*100,33+i*100,34+i*100,cscount[5]);
	fprintf(out,"\n");
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",17+i*100,23+i*100,24+i*100,50+i*100,49+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",18+i*100,25+i*100,26+i*100,51+i*100,50+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",19+i*100,32+i*100,49+i*100,55+i*100,52+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",20+i*100,55+i*100,51+i*100,27+i*100,54+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",21+i*100,52+i*100,53+i*100,30+i*100,31+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",22+i*100,53+i*100,54+i*100,28+i*100,29+i*100);
	fprintf(out,"\n");
	fprintf(out,"# spojeni rezu 0 a 1 nebo 3 a 1 - linearni spojeni\n");
	for(j = 0; j < 32; j++){
	  fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10001+j+i*100,j+1+(i-1)*100,j+1+i*100,lcount[i]);
	}
	fprintf(out,"\n");
	for(j = 0; j < 32; j++){
	  if(j != 21 && j != 31){
	    fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,j+10002+i*100,j+1+i*100,j+10001+i*100);
	  }
	  else{
	    if(j == 21){
	      fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,10001+i*100,j+1+i*100,j+10001+i*100);
	    }
	    if(j == 31){
	      fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,10023+i*100,j+1+i*100,j+10001+i*100);
	    }
	  }
	}
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10033+i*100,(i-1)*100+33,10003+i*100,i*100+33,10022+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10034+i*100,(i-1)*100+34,10004+i*100,i*100+34,10021+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10035+i*100,(i-1)*100+35,10024+i*100,i*100+35,10020+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10036+i*100,(i-1)*100+36,10023+i*100,i*100+36,10019+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10037+i*100,(i-1)*100+37,10006+i*100,i*100+37,10026+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10038+i*100,(i-1)*100+38,10007+i*100,i*100+38,10027+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10039+i*100,(i-1)*100+39,10008+i*100,i*100+39,10028+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10040+i*100,(i-1)*100+40,10009+i*100,i*100+40,10029+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10041+i*100,(i-1)*100+41,10032+i*100,i*100+41,10018+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10042+i*100,(i-1)*100+42,10031+i*100,i*100+42,10017+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10043+i*100,(i-1)*100+43,10011+i*100,i*100+43,10016+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10044+i*100,(i-1)*100+44,10012+i*100,i*100+44,10015+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10045+i*100,(i-1)*100+45,10004+i*100,i*100+45,10024+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10046+i*100,(i-1)*100+46,10031+i*100,i*100+46,10011+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10047+i*100,(i-1)*100+47,10005+i*100,i*100+47,10025+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10048+i*100,(i-1)*100+48,10030+i*100,i*100+48,10010+i*100);
	fprintf(out,"\n");
	fprintf(out,"# definice regionu\n");
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",1+i*100,10001+i*100,10002+i*100,10033+i*100,10022+i*100,(i-1)*100+1,i*100+1);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",2+i*100,10033+i*100,10003+i*100,10034+i*100,10021+i*100,(i-1)*100+2,i*100+2);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",3+i*100,10034+i*100,10045+i*100,10035+i*100,10020+i*100, 3+(i-1)*100,3+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",4+i*100,10004+i*100,10047+i*100,10024+i*100,10045+i*100, 4+(i-1)*100,4+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld  hexa map yes equidistant size def\n",5+i*100,10005+i*100,10037+i*100,10025+i*100,10047+i*100,5+(i-1)*100,5+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld  hexa map yes equidistant size def\n",6+i*100,10035+i*100,10023+i*100,10036+i*100,10019+i*100,6+(i-1)*100,6+i*100);
	fprintf(out,"region %ld boundary surface %ld %ld -%ld -%ld %ld -%ld  hexa map yes equidistant size def\n",7+i*100,10038+i*100,10026+i*100,10037+i*100,10006+i*100,7+(i-1)*100,7+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld %ld %ld -%ld hexa map yes equidistant size def\n",8+i*100,10038+i*100,10007+i*100,10039+i*100,10027+i*100,8+(i-1)*100,8+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",9+i*100,10036+i*100,10032+i*100,10041+i*100,10018+i*100,9+(i-1)*100,9+i*100 );
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld %ld %ld -%ld hexa map yes equidistant size def\n",10+i*100,10039+i*100,10008+i*100,10040+i*100,10028+i*100,10+(i-1)*100,10+i*100);
	fprintf(out,"region %ld boundary surface %ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",11+i*100,10029+i*100,10040+i*100,10009+i*100,10048+i*100,11+(i-1)*100,11+i*100);
	fprintf(out,"region %ld boundary surface %ld %ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",12+i*100,10030+i*100,10048+i*100,10010+i*100,10046+i*100,12+(i-1)*100,12+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",13+i*100,10041+i*100,10031+i*100,10042+i*100,10017+i*100,13+(i-1)*100,13+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",14+i*100,10042+i*100,10046+i*100,10043+i*100,10016+i*100,14+(i-1)*100,14+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",15+i*100,10043+i*100,10011+i*100,10044+i*100,10015+i*100,15+(i-1)*100,15+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",16+i*100,10044+i*100,10012+i*100,10013+i*100,10014+i*100,16+(i-1)*100,16+i*100);
	}
	if(typecut[i-1] == 3 || typecut[i-1] == 5){
	  fprintf(out,"# %ld rez - typ %ld\n",i+1,typecut[i]);
	for(j = 0; j < 32; j++){
	  ver=j+1+i*100;
	  fprintf(out,"vertex %ld xyz %lf %lf %lf\n",ver,x[i][j],y[i][j],z[i][j]);
	}
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",1+i*100,1+i*100,2+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  2+i*100,  2+i*100,  3+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  3+i*100,  3+i*100,  4+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  4+i*100,  4+i*100,  5+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  5+i*100,  5+i*100,  6+i*100,cscount[7]); 
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  6+i*100,  6+i*100,  7+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  7+i*100,  7+i*100,  8+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  8+i*100,  8+i*100,  9+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  9+i*100,  9+i*100, 10+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 10+i*100, 10+i*100,11+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 11+i*100, 11+i*100,12+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 12+i*100, 12+i*100,13+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 13+i*100, 13+i*100,14+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 14+i*100, 14+i*100,15+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 15+i*100, 15+i*100,16+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 16+i*100, 16+i*100,17+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 17+i*100, 17+i*100,18+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 18+i*100, 18+i*100,19+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 19+i*100, 19+i*100,20+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 20+i*100, 20+i*100,21+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 21+i*100, 21+i*100,22+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 22+i*100, 22+i*100, 1+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 23+i*100, 23+i*100,24+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 24+i*100, 24+i*100,25+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 25+i*100, 25+i*100,26+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 26+i*100, 26+i*100,27+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 27+i*100, 27+i*100,28+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 28+i*100, 28+i*100,29+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 29+i*100, 29+i*100,30+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 30+i*100, 30+i*100,31+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 31+i*100, 31+i*100,32+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 32+i*100, 32+i*100,23+i*100,cscount[5]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 33+i*100,  3+i*100, 22+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 34+i*100,  4+i*100, 21+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 35+i*100, 24+i*100, 20+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 36+i*100, 23+i*100, 19+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 37+i*100,  6+i*100, 26+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 38+i*100,  7+i*100, 27+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 39+i*100,  8+i*100, 28+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 40+i*100,  9+i*100, 29+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 41+i*100, 32+i*100, 18+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 42+i*100, 31+i*100, 17+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 43+i*100, 11+i*100, 16+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 44+i*100, 12+i*100, 15+i*100,cscount[0]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 45+i*100,  4+i*100, 24+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 46+i*100, 31+i*100, 11+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 47+i*100,  5+i*100, 25+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 48+i*100, 30+i*100, 10+i*100,cscount[6]);
	fprintf(out,"\n");
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  1+i*100, 1+i*100, 2+i*100,33+i*100,22+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  2+i*100,33+i*100, 3+i*100,34+i*100,21+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  3+i*100,20+i*100,34+i*100,45+i*100,35+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  4+i*100,45+i*100, 4+i*100,47+i*100,24+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  5+i*100,47+i*100, 5+i*100,37+i*100,25+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  6+i*100,35+i*100,23+i*100,36+i*100,19+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  7+i*100,37+i*100, 6+i*100,38+i*100,26+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  8+i*100,38+i*100, 7+i*100,39+i*100,27+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  9+i*100,36+i*100,32+i*100,41+i*100,18+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 10+i*100,39+i*100, 8+i*100,40+i*100,28+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 11+i*100,40+i*100, 9+i*100,48+i*100,29+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 12+i*100,48+i*100,10+i*100,46+i*100,30+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 13+i*100,41+i*100,31+i*100,42+i*100,17+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 14+i*100,46+i*100,43+i*100,16+i*100,42+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 15+i*100,43+i*100,11+i*100,44+i*100,15+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 16+i*100,44+i*100,12+i*100,13+i*100,14+i*100);
	fprintf(out,"\n");
	fprintf(out,"vertex %ld xyz %lf %lf %lf\n",33+i*100,x[i][24],y[i][22],z[i][24]);
	fprintf(out,"vertex %ld xyz %lf %lf %lf\n",34+i*100,x[i][29],y[i][31],z[i][29]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",49+i*100,23+i*100,33+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",50+i*100,33+i*100,25+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",51+i*100,33+i*100,27+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",52+i*100,32+i*100,34+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",53+i*100,34+i*100,30+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",54+i*100,34+i*100,28+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",55+i*100,33+i*100,34+i*100,cscount[5]);
	fprintf(out,"\n");
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",17+i*100,23+i*100,24+i*100,50+i*100,49+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",18+i*100,25+i*100,26+i*100,51+i*100,50+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",19+i*100,32+i*100,49+i*100,55+i*100,52+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",20+i*100,55+i*100,51+i*100,27+i*100,54+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",21+i*100,52+i*100,53+i*100,30+i*100,31+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",22+i*100,53+i*100,54+i*100,28+i*100,29+i*100);
	fprintf(out,"\n");
	fprintf(out,"# spojeni rezu %ld (typ %ld) a %ld (typ %ld)\n",i,typecut[i-1],i+1,typecut[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10001+i*100,1+(i-1)*100,1+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10002+i*100,2+(i-1)*100,2+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10003+i*100,3+(i-1)*100,3+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10004+i*100,4+(i-1)*100,4+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10005+i*100,5+(i-1)*100,5+i*100,lcount[i]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10006+i*100,6+(i-1)*100,6+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][5]+(x[i][5]-x[i-1][5])/2,y[i-1][5]+(y[i][5]-y[i-1][5])/2,z[i-1][5]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10007+i*100,7+(i-1)*100,7+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][6]+(x[i][6]-x[i-1][6])/2,y[i-1][6]+(y[i][6]-y[i-1][6])/2,z[i-1][6]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10008+i*100,8+(i-1)*100,8+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][7]+(x[i][7]-x[i-1][7])/2,y[i-1][7]+(y[i][7]-y[i-1][7])/2,z[i-1][7]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10009+i*100,9+(i-1)*100,9+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][8]+(x[i][8]-x[i-1][8])/2,y[i-1][8]+(y[i][8]-y[i-1][8])/2,z[i-1][8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10010+i*100,10+(i-1)*100,10+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10011+i*100,11+(i-1)*100,11+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10012+i*100,12+(i-1)*100,12+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10013+i*100,13+(i-1)*100,13+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10014+i*100,14+(i-1)*100,14+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10015+i*100,15+(i-1)*100,15+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10016+i*100,16+(i-1)*100,16+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10017+i*100,17+(i-1)*100,17+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10018+i*100,18+(i-1)*100,18+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10019+i*100,19+(i-1)*100,19+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10020+i*100,20+(i-1)*100,20+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10021+i*100,21+(i-1)*100,21+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10022+i*100,22+(i-1)*100,22+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10023+i*100,23+(i-1)*100,23+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10024+i*100,24+(i-1)*100,24+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10025+i*100,25+(i-1)*100,25+i*100,lcount[i]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10026+i*100,26+(i-1)*100,26+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][25]+(x[i][25]-x[i-1][25])/2,y[i-1][25]+(y[i][25]-y[i-1][25])/2,z[i-1][25]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10027+i*100,27+(i-1)*100,27+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][26]+(x[i][26]-x[i-1][26])/2,y[i-1][26]+(y[i][26]-y[i-1][26])/2,z[i-1][26]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10028+i*100,28+(i-1)*100,28+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][27]+(x[i][27]-x[i-1][27])/2,y[i-1][27]+(y[i][27]-y[i-1][27])/2,z[i-1][27]);      
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10029+i*100,29+(i-1)*100,29+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][28]+(x[i][28]-x[i-1][28])/2,y[i-1][28]+(y[i][28]-y[i-1][28])/2,z[i-1][28]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10030+i*100,30+(i-1)*100,30+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10031+i*100,31+(i-1)*100,31+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10032+i*100,32+(i-1)*100,32+i*100,lcount[i]);
	for(j = 0; j < 32; j++){
	  if(j != 21 && j != 31){
	    fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,j+10002+i*100,j+1+i*100,j+10001+i*100);
	  }
	  else{
	    if(j == 21){
	      fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,10001+i*100,j+1+i*100,j+10001+i*100);
	    }
	    if(j == 31){
	      fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,10023+i*100,j+1+i*100,j+10001+i*100);
	    }
	  }
	}
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10033+i*100,(i-1)*100+33,10003+i*100,i*100+33,10022+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10034+i*100,(i-1)*100+34,10004+i*100,i*100+34,10021+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10035+i*100,(i-1)*100+35,10024+i*100,i*100+35,10020+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10036+i*100,(i-1)*100+36,10023+i*100,i*100+36,10019+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10037+i*100,(i-1)*100+37,10006+i*100,i*100+37,10026+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10038+i*100,(i-1)*100+38,10007+i*100,i*100+38,10027+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10039+i*100,(i-1)*100+39,10008+i*100,i*100+39,10028+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10040+i*100,(i-1)*100+40,10009+i*100,i*100+40,10029+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10041+i*100,(i-1)*100+41,10032+i*100,i*100+41,10018+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10042+i*100,(i-1)*100+42,10031+i*100,i*100+42,10017+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10043+i*100,(i-1)*100+43,10011+i*100,i*100+43,10016+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10044+i*100,(i-1)*100+44,10012+i*100,i*100+44,10015+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10045+i*100,(i-1)*100+45,10004+i*100,i*100+45,10024+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10046+i*100,(i-1)*100+46,10031+i*100,i*100+46,10011+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10047+i*100,(i-1)*100+47,10005+i*100,i*100+47,10025+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10048+i*100,(i-1)*100+48,10030+i*100,i*100+48,10010+i*100);
	fprintf(out,"\n# definice regionu\n");
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",1+i*100,10001+i*100,10002+i*100,10033+i*100,10022+i*100,(i-1)*100+1,i*100+1);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",2+i*100,10033+i*100,10003+i*100,10034+i*100,10021+i*100,(i-1)*100+2,i*100+2);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",3+i*100,10034+i*100,10045+i*100,10035+i*100,10020+i*100, 3+(i-1)*100,3+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",4+i*100,10004+i*100,10047+i*100,10024+i*100,10045+i*100, 4+(i-1)*100,4+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld  hexa map yes equidistant size def\n",5+i*100,10005+i*100,10037+i*100,10025+i*100,10047+i*100,5+(i-1)*100,5+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld  hexa map yes equidistant size def\n",6+i*100,10035+i*100,10023+i*100,10036+i*100,10019+i*100,6+(i-1)*100,6+i*100);
	fprintf(out,"region %ld boundary surface %ld %ld -%ld -%ld %ld -%ld  hexa map yes equidistant size def\n",7+i*100,10038+i*100,10026+i*100,10037+i*100,10006+i*100,7+(i-1)*100,7+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld %ld %ld -%ld hexa map yes equidistant size def\n",8+i*100,10038+i*100,10007+i*100,10039+i*100,10027+i*100,8+(i-1)*100,8+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",9+i*100,10036+i*100,10032+i*100,10041+i*100,10018+i*100,9+(i-1)*100,9+i*100 );
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld %ld %ld -%ld hexa map yes equidistant size def\n",10+i*100,10039+i*100,10008+i*100,10040+i*100,10028+i*100,10+(i-1)*100,10+i*100);
	fprintf(out,"region %ld boundary surface %ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",11+i*100,10029+i*100,10040+i*100,10009+i*100,10048+i*100,11+(i-1)*100,11+i*100);
	fprintf(out,"region %ld boundary surface %ld %ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",12+i*100,10030+i*100,10048+i*100,10010+i*100,10046+i*100,12+(i-1)*100,12+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",13+i*100,10041+i*100,10031+i*100,10042+i*100,10017+i*100,13+(i-1)*100,13+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",14+i*100,10042+i*100,10046+i*100,10043+i*100,10016+i*100,14+(i-1)*100,14+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",15+i*100,10043+i*100,10011+i*100,10044+i*100,10015+i*100,15+(i-1)*100,15+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",16+i*100,10044+i*100,10012+i*100,10013+i*100,10014+i*100,16+(i-1)*100,16+i*100);
	}
      }
    }
    // konec ztuzidla
    if(typecut[i] == 3){
      fprintf(out,"# %ld rez - typ %ld\n",i+1,typecut[i]);
	for(j = 0; j < 32; j++){
	  ver=j+1+i*100;
	  fprintf(out,"vertex %ld xyz %lf %lf %lf\n",ver,x[i][j],y[i][j],z[i][j]);
	}
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",1+i*100,1+i*100,2+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  2+i*100,  2+i*100,  3+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  3+i*100,  3+i*100,  4+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  4+i*100,  4+i*100,  5+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  5+i*100,  5+i*100,  6+i*100,cscount[7]); 
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  6+i*100,  6+i*100,  7+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  7+i*100,  7+i*100,  8+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  8+i*100,  8+i*100,  9+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  9+i*100,  9+i*100, 10+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 10+i*100, 10+i*100,11+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 11+i*100, 11+i*100,12+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 12+i*100, 12+i*100,13+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 13+i*100, 13+i*100,14+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 14+i*100, 14+i*100,15+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 15+i*100, 15+i*100,16+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 16+i*100, 16+i*100,17+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 17+i*100, 17+i*100,18+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 18+i*100, 18+i*100,19+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 19+i*100, 19+i*100,20+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 20+i*100, 20+i*100,21+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 21+i*100, 21+i*100,22+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 22+i*100, 22+i*100, 1+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 23+i*100, 23+i*100,24+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 24+i*100, 24+i*100,25+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 25+i*100, 25+i*100,26+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 26+i*100, 26+i*100,27+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 27+i*100, 27+i*100,28+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 28+i*100, 28+i*100,29+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 29+i*100, 29+i*100,30+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 30+i*100, 30+i*100,31+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 31+i*100, 31+i*100,32+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 32+i*100, 32+i*100,23+i*100,cscount[5]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 33+i*100,  3+i*100, 22+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 34+i*100,  4+i*100, 21+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 35+i*100, 24+i*100, 20+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 36+i*100, 23+i*100, 19+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 37+i*100,  6+i*100, 26+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 38+i*100,  7+i*100, 27+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 39+i*100,  8+i*100, 28+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 40+i*100,  9+i*100, 29+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 41+i*100, 32+i*100, 18+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 42+i*100, 31+i*100, 17+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 43+i*100, 11+i*100, 16+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 44+i*100, 12+i*100, 15+i*100,cscount[0]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 45+i*100,  4+i*100, 24+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 46+i*100, 31+i*100, 11+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 47+i*100,  5+i*100, 25+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 48+i*100, 30+i*100, 10+i*100,cscount[6]);
	fprintf(out,"\n");
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  1+i*100, 1+i*100, 2+i*100,33+i*100,22+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  2+i*100,33+i*100, 3+i*100,34+i*100,21+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  3+i*100,20+i*100,34+i*100,45+i*100,35+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  4+i*100,45+i*100, 4+i*100,47+i*100,24+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  5+i*100,47+i*100, 5+i*100,37+i*100,25+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  6+i*100,35+i*100,23+i*100,36+i*100,19+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  7+i*100,37+i*100, 6+i*100,38+i*100,26+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  8+i*100,38+i*100, 7+i*100,39+i*100,27+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  9+i*100,36+i*100,32+i*100,41+i*100,18+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 10+i*100,39+i*100, 8+i*100,40+i*100,28+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 11+i*100,40+i*100, 9+i*100,48+i*100,29+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 12+i*100,48+i*100,10+i*100,46+i*100,30+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 13+i*100,41+i*100,31+i*100,42+i*100,17+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 14+i*100,46+i*100,43+i*100,16+i*100,42+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 15+i*100,43+i*100,11+i*100,44+i*100,15+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 16+i*100,44+i*100,12+i*100,13+i*100,14+i*100);
	fprintf(out,"\n");
	fprintf(out,"vertex %ld xyz %lf %lf %lf\n",33+i*100,x[i][24],y[i][22],z[i][24]);
	fprintf(out,"vertex %ld xyz %lf %lf %lf\n",34+i*100,x[i][29],y[i][31],z[i][29]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",49+i*100,23+i*100,33+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",50+i*100,33+i*100,25+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",51+i*100,33+i*100,27+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",52+i*100,32+i*100,34+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",53+i*100,34+i*100,30+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",54+i*100,34+i*100,28+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",55+i*100,33+i*100,34+i*100,cscount[5]);
	fprintf(out,"\n");
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",17+i*100,23+i*100,24+i*100,50+i*100,49+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",18+i*100,25+i*100,26+i*100,51+i*100,50+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",19+i*100,32+i*100,49+i*100,55+i*100,52+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",20+i*100,55+i*100,51+i*100,27+i*100,54+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",21+i*100,52+i*100,53+i*100,30+i*100,31+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",22+i*100,53+i*100,54+i*100,28+i*100,29+i*100);
	fprintf(out,"\n");
	for(j = 0; j < 32; j++){
	  fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10001+j+i*100,j+1+(i-1)*100,j+1+i*100,lcount[i]);
	}
	fprintf(out,"\n");
	for(j = 0; j < 32; j++){
	  if(j != 21 && j != 31){
	    fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,j+10002+i*100,j+1+i*100,j+10001+i*100);
	  }
	  else{
	    if(j == 21){
	      fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,10001+i*100,j+1+i*100,j+10001+i*100);
	    }
	    if(j == 31){
	      fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,10023+i*100,j+1+i*100,j+10001+i*100);
	    }
	  }
	}
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10033+i*100,33+(i-1)*100,33+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10034+i*100,34+(i-1)*100,34+i*100,lcount[i]);
	
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10033+i*100,(i-1)*100+33,10003+i*100,i*100+33,10022+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10034+i*100,(i-1)*100+34,10004+i*100,i*100+34,10021+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10035+i*100,(i-1)*100+35,10024+i*100,i*100+35,10020+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10036+i*100,(i-1)*100+36,10023+i*100,i*100+36,10019+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10037+i*100,(i-1)*100+37,10006+i*100,i*100+37,10026+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10038+i*100,(i-1)*100+38,10007+i*100,i*100+38,10027+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10039+i*100,(i-1)*100+39,10008+i*100,i*100+39,10028+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10040+i*100,(i-1)*100+40,10009+i*100,i*100+40,10029+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10041+i*100,(i-1)*100+41,10032+i*100,i*100+41,10018+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10042+i*100,(i-1)*100+42,10031+i*100,i*100+42,10017+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10043+i*100,(i-1)*100+43,10011+i*100,i*100+43,10016+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10044+i*100,(i-1)*100+44,10012+i*100,i*100+44,10015+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10045+i*100,(i-1)*100+45,10004+i*100,i*100+45,10024+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10046+i*100,(i-1)*100+46,10031+i*100,i*100+46,10011+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10047+i*100,(i-1)*100+47,10005+i*100,i*100+47,10025+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10048+i*100,(i-1)*100+48,10030+i*100,i*100+48,10010+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10049+i*100,(i-1)*100+49,10033+i*100,i*100+49,10023+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10050+i*100,(i-1)*100+50,10033+i*100,i*100+50,10025+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10051+i*100,(i-1)*100+51,10033+i*100,i*100+51,10027+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10052+i*100,(i-1)*100+52,10034+i*100,i*100+52,10032+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10053+i*100,(i-1)*100+53,10034+i*100,i*100+53,10030+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10054+i*100,(i-1)*100+54,10034+i*100,i*100+54,10028+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10055+i*100,(i-1)*100+55,10033+i*100,i*100+55,10034+i*100);
	fprintf(out,"\n# definice regionu\n");
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",1+i*100,10001+i*100,10002+i*100,10033+i*100,10022+i*100,(i-1)*100+1,i*100+1);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",2+i*100,10033+i*100,10003+i*100,10034+i*100,10021+i*100,(i-1)*100+2,i*100+2);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",3+i*100,10034+i*100,10045+i*100,10035+i*100,10020+i*100, 3+(i-1)*100,3+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",4+i*100,10004+i*100,10047+i*100,10024+i*100,10045+i*100, 4+(i-1)*100,4+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld  hexa map yes equidistant size def\n",5+i*100,10005+i*100,10037+i*100,10025+i*100,10047+i*100,5+(i-1)*100,5+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld  hexa map yes equidistant size def\n",6+i*100,10035+i*100,10023+i*100,10036+i*100,10019+i*100,6+(i-1)*100,6+i*100);
	fprintf(out,"region %ld boundary surface %ld %ld -%ld -%ld %ld -%ld  hexa map yes equidistant size def\n",7+i*100,10038+i*100,10026+i*100,10037+i*100,10006+i*100,7+(i-1)*100,7+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld %ld %ld -%ld hexa map yes equidistant size def\n",8+i*100,10038+i*100,10007+i*100,10039+i*100,10027+i*100,8+(i-1)*100,8+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",9+i*100,10036+i*100,10032+i*100,10041+i*100,10018+i*100,9+(i-1)*100,9+i*100 );
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld %ld %ld -%ld hexa map yes equidistant size def\n",10+i*100,10039+i*100,10008+i*100,10040+i*100,10028+i*100,10+(i-1)*100,10+i*100);
	fprintf(out,"region %ld boundary surface %ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",11+i*100,10029+i*100,10040+i*100,10009+i*100,10048+i*100,11+(i-1)*100,11+i*100);
	fprintf(out,"region %ld boundary surface %ld %ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",12+i*100,10030+i*100,10048+i*100,10010+i*100,10046+i*100,12+(i-1)*100,12+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",13+i*100,10041+i*100,10031+i*100,10042+i*100,10017+i*100,13+(i-1)*100,13+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",14+i*100,10042+i*100,10046+i*100,10043+i*100,10016+i*100,14+(i-1)*100,14+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",15+i*100,10043+i*100,10011+i*100,10044+i*100,10015+i*100,15+(i-1)*100,15+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",16+i*100,10044+i*100,10012+i*100,10013+i*100,10014+i*100,16+(i-1)*100,16+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",17+i*100,10024+i*100,10050+i*100,10049+i*100,10023+i*100,17+(i-1)*100,17+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld -%ld %ld %ld -%ld hexa map yes equidistant size def\n",18+i*100,10025+i*100, 10026+i*100, 10051+i*100, 10050+i*100,18+(i-1)*100, 18+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",19+i*100,10049+i*100, 10055+i*100, 10052+i*100, 10032+i*100, 19+(i-1)*100, 19+i*100);
	fprintf(out,"region %ld boundary surface %ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",20+i*100,10051+i*100, 10027+i*100, 10054+i*100, 10055+i*100, 20+(i-1)*100, 20+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",21+i*100,10052+i*100, 10053+i*100, 10030+i*100, 10031+i*100, 21+(i-1)*100, 21+i*100);
	fprintf(out,"region %ld boundary surface %ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",22+i*100,10054+i*100, 10028+i*100, 10029+i*100, 10053+i*100, 22+(i-1)*100, 22+i*100);
    }
    // konec zarodku - zacatek noveho pole
    if(typecut[i] == 4){
      if(i == 0){
	// prvni rez
	fprintf(out,"# %ld rez - typ %ld\n",i+1,typecut[i]);
	for(j = 0; j < 32; j++){
	  ver=j+1;
	  fprintf(out,"vertex %ld xyz %lf %lf %lf\n",ver,x[0][j],y[0][j],z[0][j]);
	}
	fprintf(out,"\n");
	fprintf(out,"curve  1 vertex  1   2 count %ld\n",cscount[0]);
	fprintf(out,"curve  2 vertex  2   3 count %ld\n",cscount[1]);
	fprintf(out,"curve  3 vertex  3   4 count %ld\n",cscount[2]);
	fprintf(out,"curve  4 vertex  4   5 count %ld\n",cscount[7]);
	fprintf(out,"curve  5 vertex  5   6 count %ld\n",cscount[7]); 
	fprintf(out,"curve  6 vertex  6   7 count %ld\n",cscount[8]);
	fprintf(out,"curve  7 vertex  7   8 count %ld\n",cscount[5]);
	fprintf(out,"curve  8 vertex  8   9 count %ld\n",cscount[8]);
	fprintf(out,"curve  9 vertex  9  10 count %ld\n",cscount[7]);
	fprintf(out,"curve 10 vertex 10 11 count %ld\n",cscount[7]);
	fprintf(out,"curve 11 vertex 11 12 count %ld\n",cscount[3]);
	fprintf(out,"curve 12 vertex 12 13 count %ld\n",cscount[4]);
	fprintf(out,"curve 13 vertex 13 14 count  %ld\n",cscount[0]);
	fprintf(out,"curve 14 vertex 14 15 count %ld\n",cscount[4]);
	fprintf(out,"curve 15 vertex 15 16 count %ld\n",cscount[3]);
	fprintf(out,"curve 16 vertex 16 17 count %ld\n",cscount[6]);
	fprintf(out,"curve 17 vertex 17 18 count %ld\n",cscount[8]);
	fprintf(out,"curve 18 vertex 18 19 count %ld\n",cscount[5]);
	fprintf(out,"curve 19 vertex 19 20 count %ld\n",cscount[8]);
	fprintf(out,"curve 20 vertex 20 21 count %ld\n",cscount[6]);
	fprintf(out,"curve 21 vertex 21 22 count %ld\n",cscount[2]);
	fprintf(out,"curve 22 vertex 22  1 count %ld\n",cscount[1]);
	fprintf(out,"curve 23 vertex 23 24 count %ld\n",cscount[8]);
	fprintf(out,"curve 24 vertex 24 25 count %ld\n",cscount[7]);
	fprintf(out,"curve 25 vertex 25 26 count %ld\n",cscount[7]);
	fprintf(out,"curve 26 vertex 26 27 count %ld\n",cscount[8]);
	fprintf(out,"curve 27 vertex 27 28 count %ld\n",cscount[5]);
	fprintf(out,"curve 28 vertex 28 29 count %ld\n",cscount[8]);
	fprintf(out,"curve 29 vertex 29 30 count %ld\n",cscount[7]);
	fprintf(out,"curve 30 vertex 30 31 count %ld\n",cscount[7]);
	fprintf(out,"curve 31 vertex 31 32 count %ld\n",cscount[8]);
	fprintf(out,"curve 32 vertex 32 23 count %ld\n",cscount[5]);
	fprintf(out,"\n");
	fprintf(out,"curve 33 vertex  3 22 count %ld \n",cscount[0]);
	fprintf(out,"curve 34 vertex  4 21 count %ld \n",cscount[0]);
	fprintf(out,"curve 35 vertex 24 20 count %ld \n",cscount[0]);
	fprintf(out,"curve 36 vertex 23 19 count %ld \n",cscount[0]);
	fprintf(out,"curve 37 vertex  6 26 count %ld\n",cscount[6]);
	fprintf(out,"curve 38 vertex  7 27 count %ld\n",cscount[6]);
	fprintf(out,"curve 39 vertex  8 28 count %ld\n",cscount[6]);
	fprintf(out,"curve 40 vertex  9 29 count %ld\n",cscount[6]);
	fprintf(out,"curve 41 vertex 32 18 count %ld \n",cscount[0]);
	fprintf(out,"curve 42 vertex 31 17 count %ld \n",cscount[0]);
	fprintf(out,"curve 43 vertex 11 16 count %ld \n",cscount[0]);
	fprintf(out,"curve 44 vertex 12 15 count %ld \n",cscount[0]);
	fprintf(out,"\n");
	fprintf(out,"curve 45 vertex  4 24 count %ld\n",cscount[6]);
	fprintf(out,"curve 46 vertex 31 11 count %ld\n",cscount[6]);
	fprintf(out,"curve 47 vertex  5 25 count %ld\n",cscount[6]);
	fprintf(out,"curve 48 vertex 30 10 count %ld\n",cscount[6]);
	fprintf(out,"\n");
	fprintf(out,"surface 1 curve 1 2 33 22\n");
	fprintf(out,"surface 2 curve 33 3 34 21\n");
	fprintf(out,"surface 3 curve 20 34 45 35\n");
	fprintf(out,"surface 4 curve 45  4 47 24\n");
	fprintf(out,"surface 5 curve 47  5 37 25\n");
	fprintf(out,"surface 6 curve 35 23 36 19\n");
	fprintf(out,"surface 7 curve 37  6 38 26\n");
	fprintf(out,"surface 8 curve 38  7 39 27\n");
	fprintf(out,"surface 9 curve 36 32 41 18\n");
	fprintf(out,"surface 10 curve 39 8 40 28\n");
	fprintf(out,"surface 11 curve 40  9 48 29\n");
	fprintf(out,"surface 12 curve 48 10 46 30\n");
	fprintf(out,"surface 13 curve 41 31 42 17\n");
	fprintf(out,"surface 14 curve 46 43 16 42\n");
	fprintf(out,"surface 15 curve 43 11 44 15\n");
	fprintf(out,"surface 16 curve 44 12 13 14\n");
	fprintf(out,"\n");
	}
      if(i != 0){
	fprintf(out,"# %ld rez - typ %ld\n",i+1,typecut[i]);
	for(j = 0; j < 32; j++){
	  ver=j+1+i*100;
	  fprintf(out,"vertex %ld xyz %lf %lf %lf\n",ver,x[i][j],y[i][j],z[i][j]);
	}
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",1+i*100,1+i*100,2+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  2+i*100,  2+i*100,  3+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  3+i*100,  3+i*100,  4+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  4+i*100,  4+i*100,  5+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  5+i*100,  5+i*100,  6+i*100,cscount[7]); 
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  6+i*100,  6+i*100,  7+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  7+i*100,  7+i*100,  8+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  8+i*100,  8+i*100,  9+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  9+i*100,  9+i*100, 10+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 10+i*100, 10+i*100,11+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 11+i*100, 11+i*100,12+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 12+i*100, 12+i*100,13+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 13+i*100, 13+i*100,14+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 14+i*100, 14+i*100,15+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 15+i*100, 15+i*100,16+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 16+i*100, 16+i*100,17+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 17+i*100, 17+i*100,18+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 18+i*100, 18+i*100,19+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 19+i*100, 19+i*100,20+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 20+i*100, 20+i*100,21+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 21+i*100, 21+i*100,22+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 22+i*100, 22+i*100, 1+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 23+i*100, 23+i*100,24+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 24+i*100, 24+i*100,25+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 25+i*100, 25+i*100,26+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 26+i*100, 26+i*100,27+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 27+i*100, 27+i*100,28+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 28+i*100, 28+i*100,29+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 29+i*100, 29+i*100,30+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 30+i*100, 30+i*100,31+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 31+i*100, 31+i*100,32+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 32+i*100, 32+i*100,23+i*100,cscount[5]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 33+i*100,  3+i*100, 22+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 34+i*100,  4+i*100, 21+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 35+i*100, 24+i*100, 20+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 36+i*100, 23+i*100, 19+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 37+i*100,  6+i*100, 26+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 38+i*100,  7+i*100, 27+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 39+i*100,  8+i*100, 28+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 40+i*100,  9+i*100, 29+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 41+i*100, 32+i*100, 18+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 42+i*100, 31+i*100, 17+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 43+i*100, 11+i*100, 16+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 44+i*100, 12+i*100, 15+i*100,cscount[0]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 45+i*100,  4+i*100, 24+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 46+i*100, 31+i*100, 11+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 47+i*100,  5+i*100, 25+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 48+i*100, 30+i*100, 10+i*100,cscount[6]);
	fprintf(out,"\n");
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  1+i*100, 1+i*100, 2+i*100,33+i*100,22+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  2+i*100,33+i*100, 3+i*100,34+i*100,21+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  3+i*100,20+i*100,34+i*100,45+i*100,35+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  4+i*100,45+i*100, 4+i*100,47+i*100,24+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  5+i*100,47+i*100, 5+i*100,37+i*100,25+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  6+i*100,35+i*100,23+i*100,36+i*100,19+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  7+i*100,37+i*100, 6+i*100,38+i*100,26+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  8+i*100,38+i*100, 7+i*100,39+i*100,27+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  9+i*100,36+i*100,32+i*100,41+i*100,18+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 10+i*100,39+i*100, 8+i*100,40+i*100,28+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 11+i*100,40+i*100, 9+i*100,48+i*100,29+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 12+i*100,48+i*100,10+i*100,46+i*100,30+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 13+i*100,41+i*100,31+i*100,42+i*100,17+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 14+i*100,46+i*100,43+i*100,16+i*100,42+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 15+i*100,43+i*100,11+i*100,44+i*100,15+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 16+i*100,44+i*100,12+i*100,13+i*100,14+i*100);
	fprintf(out,"\n");
	fprintf(out,"# spojeni rezu 0 a 1 nebo 3 a 1 - linearni spojeni\n");
	for(j = 0; j < 32; j++){
	  fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10001+j+i*100,j+1+(i-1)*100,j+1+i*100,lcount[i]);
	}
	fprintf(out,"\n");
 	for(j = 0; j < 32; j++){
	  if(j != 21 && j != 31){
	    fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,j+10002+i*100,j+1+i*100,j+10001+i*100);
	  }
	  else{
	    if(j == 21){
	      fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,10001+i*100,j+1+i*100,j+10001+i*100);
	    }
	    if(j == 31){
	      fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,10023+i*100,j+1+i*100,j+10001+i*100);
	    }
	  }
	}
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10033+i*100,(i-1)*100+33,10003+i*100,i*100+33,10022+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10034+i*100,(i-1)*100+34,10004+i*100,i*100+34,10021+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10035+i*100,(i-1)*100+35,10024+i*100,i*100+35,10020+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10036+i*100,(i-1)*100+36,10023+i*100,i*100+36,10019+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10037+i*100,(i-1)*100+37,10006+i*100,i*100+37,10026+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10038+i*100,(i-1)*100+38,10007+i*100,i*100+38,10027+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10039+i*100,(i-1)*100+39,10008+i*100,i*100+39,10028+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10040+i*100,(i-1)*100+40,10009+i*100,i*100+40,10029+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10041+i*100,(i-1)*100+41,10032+i*100,i*100+41,10018+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10042+i*100,(i-1)*100+42,10031+i*100,i*100+42,10017+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10043+i*100,(i-1)*100+43,10011+i*100,i*100+43,10016+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10044+i*100,(i-1)*100+44,10012+i*100,i*100+44,10015+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10045+i*100,(i-1)*100+45,10004+i*100,i*100+45,10024+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10046+i*100,(i-1)*100+46,10031+i*100,i*100+46,10011+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10047+i*100,(i-1)*100+47,10005+i*100,i*100+47,10025+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10048+i*100,(i-1)*100+48,10030+i*100,i*100+48,10010+i*100);
	fprintf(out,"\n# definice regionu\n");
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",1+i*100,10001+i*100,10002+i*100,10033+i*100,10022+i*100,(i-1)*100+1,i*100+1);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",2+i*100,10033+i*100,10003+i*100,10034+i*100,10021+i*100,(i-1)*100+2,i*100+2);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",3+i*100,10034+i*100,10045+i*100,10035+i*100,10020+i*100, 3+(i-1)*100,3+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",4+i*100,10004+i*100,10047+i*100,10024+i*100,10045+i*100, 4+(i-1)*100,4+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld  hexa map yes equidistant size def\n",5+i*100,10005+i*100,10037+i*100,10025+i*100,10047+i*100,5+(i-1)*100,5+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld  hexa map yes equidistant size def\n",6+i*100,10035+i*100,10023+i*100,10036+i*100,10019+i*100,6+(i-1)*100,6+i*100);
	fprintf(out,"region %ld boundary surface %ld %ld -%ld -%ld %ld -%ld  hexa map yes equidistant size def\n",7+i*100,10038+i*100,10026+i*100,10037+i*100,10006+i*100,7+(i-1)*100,7+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld %ld %ld -%ld hexa map yes equidistant size def\n",8+i*100,10038+i*100,10007+i*100,10039+i*100,10027+i*100,8+(i-1)*100,8+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",9+i*100,10036+i*100,10032+i*100,10041+i*100,10018+i*100,9+(i-1)*100,9+i*100 );
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld %ld %ld -%ld hexa map yes equidistant size def\n",10+i*100,10039+i*100,10008+i*100,10040+i*100,10028+i*100,10+(i-1)*100,10+i*100);
	fprintf(out,"region %ld boundary surface %ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",11+i*100,10029+i*100,10040+i*100,10009+i*100,10048+i*100,11+(i-1)*100,11+i*100);
	fprintf(out,"region %ld boundary surface %ld %ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",12+i*100,10030+i*100,10048+i*100,10010+i*100,10046+i*100,12+(i-1)*100,12+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",13+i*100,10041+i*100,10031+i*100,10042+i*100,10017+i*100,13+(i-1)*100,13+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",14+i*100,10042+i*100,10046+i*100,10043+i*100,10016+i*100,14+(i-1)*100,14+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",15+i*100,10043+i*100,10011+i*100,10044+i*100,10015+i*100,15+(i-1)*100,15+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",16+i*100,10044+i*100,10012+i*100,10013+i*100,10014+i*100,16+(i-1)*100,16+i*100);
     }
    }
    // rez uprostred pole
    if(typecut[i] == 5){
      fprintf(out,"# \n");
      fprintf(out,"# %ld rez - typ %ld\n",i+1,typecut[i]);
      for(j = 0; j < 32; j++){
	ver=j+1+i*100;
	fprintf(out,"vertex %ld xyz %lf %lf %lf\n",ver,x[i][j],y[i][j],z[i][j]);
      }
      fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",1+i*100,1+i*100,2+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  2+i*100,  2+i*100,  3+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  3+i*100,  3+i*100,  4+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  4+i*100,  4+i*100,  5+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  5+i*100,  5+i*100,  6+i*100,cscount[7]); 
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  6+i*100,  6+i*100,  7+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  7+i*100,  7+i*100,  8+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  8+i*100,  8+i*100,  9+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  9+i*100,  9+i*100, 10+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 10+i*100, 10+i*100,11+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 11+i*100, 11+i*100,12+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 12+i*100, 12+i*100,13+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 13+i*100, 13+i*100,14+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 14+i*100, 14+i*100,15+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 15+i*100, 15+i*100,16+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 16+i*100, 16+i*100,17+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 17+i*100, 17+i*100,18+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 18+i*100, 18+i*100,19+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 19+i*100, 19+i*100,20+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 20+i*100, 20+i*100,21+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 21+i*100, 21+i*100,22+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 22+i*100, 22+i*100, 1+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 23+i*100, 23+i*100,24+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 24+i*100, 24+i*100,25+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 25+i*100, 25+i*100,26+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 26+i*100, 26+i*100,27+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 27+i*100, 27+i*100,28+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 28+i*100, 28+i*100,29+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 29+i*100, 29+i*100,30+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 30+i*100, 30+i*100,31+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 31+i*100, 31+i*100,32+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 32+i*100, 32+i*100,23+i*100,cscount[5]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 33+i*100,  3+i*100, 22+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 34+i*100,  4+i*100, 21+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 35+i*100, 24+i*100, 20+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 36+i*100, 23+i*100, 19+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 37+i*100,  6+i*100, 26+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 38+i*100,  7+i*100, 27+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 39+i*100,  8+i*100, 28+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 40+i*100,  9+i*100, 29+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 41+i*100, 32+i*100, 18+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 42+i*100, 31+i*100, 17+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 43+i*100, 11+i*100, 16+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 44+i*100, 12+i*100, 15+i*100,cscount[0]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 45+i*100,  4+i*100, 24+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 46+i*100, 31+i*100, 11+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 47+i*100,  5+i*100, 25+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 48+i*100, 30+i*100, 10+i*100,cscount[6]);
	fprintf(out,"\n");
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  1+i*100, 1+i*100, 2+i*100,33+i*100,22+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  2+i*100,33+i*100, 3+i*100,34+i*100,21+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  3+i*100,20+i*100,34+i*100,45+i*100,35+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  4+i*100,45+i*100, 4+i*100,47+i*100,24+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  5+i*100,47+i*100, 5+i*100,37+i*100,25+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  6+i*100,35+i*100,23+i*100,36+i*100,19+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  7+i*100,37+i*100, 6+i*100,38+i*100,26+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  8+i*100,38+i*100, 7+i*100,39+i*100,27+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  9+i*100,36+i*100,32+i*100,41+i*100,18+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 10+i*100,39+i*100, 8+i*100,40+i*100,28+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 11+i*100,40+i*100, 9+i*100,48+i*100,29+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 12+i*100,48+i*100,10+i*100,46+i*100,30+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 13+i*100,41+i*100,31+i*100,42+i*100,17+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 14+i*100,46+i*100,43+i*100,16+i*100,42+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 15+i*100,43+i*100,11+i*100,44+i*100,15+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 16+i*100,44+i*100,12+i*100,13+i*100,14+i*100);
	fprintf(out,"\n");
	fprintf(out,"\n");
	fprintf(out,"# spojeni rezu %ld (typ %ld) a %ld (typ %ld)\n",i,typecut[i-1],i+1,typecut[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10001+i*100,1+(i-1)*100,1+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10002+i*100,2+(i-1)*100,2+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10003+i*100,3+(i-1)*100,3+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10004+i*100,4+(i-1)*100,4+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10005+i*100,5+(i-1)*100,5+i*100,lcount[i]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10006+i*100,6+(i-1)*100,6+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][5]+(x[i][5]-x[i-1][5])/2,y[i-1][5]+(y[i][5]-y[i-1][5])/2,z[i][5]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10007+i*100,7+(i-1)*100,7+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][6]+(x[i][6]-x[i-1][6])/2,y[i-1][6]+(y[i][6]-y[i-1][6])/2,z[i][6]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10008+i*100,8+(i-1)*100,8+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][7]+(x[i][7]-x[i-1][7])/2,y[i-1][7]+(y[i][7]-y[i-1][7])/2,z[i][7]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10009+i*100,9+(i-1)*100,9+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][8]+(x[i][8]-x[i-1][8])/2,y[i-1][8]+(y[i][8]-y[i-1][8])/2,z[i][8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10010+i*100,10+(i-1)*100,10+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10011+i*100,11+(i-1)*100,11+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10012+i*100,12+(i-1)*100,12+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10013+i*100,13+(i-1)*100,13+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10014+i*100,14+(i-1)*100,14+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10015+i*100,15+(i-1)*100,15+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10016+i*100,16+(i-1)*100,16+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10017+i*100,17+(i-1)*100,17+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10018+i*100,18+(i-1)*100,18+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10019+i*100,19+(i-1)*100,19+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10020+i*100,20+(i-1)*100,20+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10021+i*100,21+(i-1)*100,21+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10022+i*100,22+(i-1)*100,22+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10023+i*100,23+(i-1)*100,23+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10024+i*100,24+(i-1)*100,24+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10025+i*100,25+(i-1)*100,25+i*100,lcount[i]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10026+i*100,26+(i-1)*100,26+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][25]+(x[i][25]-x[i-1][25])/2,y[i-1][25]+(y[i][25]-y[i-1][25])/2,z[i][25]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10027+i*100,27+(i-1)*100,27+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][26]+(x[i][26]-x[i-1][26])/2,y[i-1][26]+(y[i][26]-y[i-1][26])/2,z[i][26]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10028+i*100,28+(i-1)*100,28+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][27]+(x[i][27]-x[i-1][27])/2,y[i-1][27]+(y[i][27]-y[i-1][27])/2,z[i][27]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10029+i*100,29+(i-1)*100,29+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][28]+(x[i][28]-x[i-1][28])/2,y[i-1][28]+(y[i][28]-y[i-1][28])/2,z[i][28]);      
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10030+i*100,30+(i-1)*100,30+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10031+i*100,31+(i-1)*100,31+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10032+i*100,32+(i-1)*100,32+i*100,lcount[i]);
	for(j = 0; j < 32; j++){
	  if(j != 21 && j != 31){
	    fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,j+10002+i*100,j+1+i*100,j+10001+i*100);
	  }
	  else{
	    if(j == 21){
	      fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,10001+i*100,j+1+i*100,j+10001+i*100);
	    }
	    if(j == 31){
	      fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,10023+i*100,j+1+i*100,j+10001+i*100);
	    }
	  }
	}
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10033+i*100,(i-1)*100+33,10003+i*100,i*100+33,10022+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10034+i*100,(i-1)*100+34,10004+i*100,i*100+34,10021+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10035+i*100,(i-1)*100+35,10024+i*100,i*100+35,10020+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10036+i*100,(i-1)*100+36,10023+i*100,i*100+36,10019+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10037+i*100,(i-1)*100+37,10006+i*100,i*100+37,10026+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10038+i*100,(i-1)*100+38,10007+i*100,i*100+38,10027+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10039+i*100,(i-1)*100+39,10008+i*100,i*100+39,10028+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10040+i*100,(i-1)*100+40,10009+i*100,i*100+40,10029+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10041+i*100,(i-1)*100+41,10032+i*100,i*100+41,10018+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10042+i*100,(i-1)*100+42,10031+i*100,i*100+42,10017+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10043+i*100,(i-1)*100+43,10011+i*100,i*100+43,10016+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10044+i*100,(i-1)*100+44,10012+i*100,i*100+44,10015+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10045+i*100,(i-1)*100+45,10004+i*100,i*100+45,10024+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10046+i*100,(i-1)*100+46,10031+i*100,i*100+46,10011+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10047+i*100,(i-1)*100+47,10005+i*100,i*100+47,10025+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10048+i*100,(i-1)*100+48,10030+i*100,i*100+48,10010+i*100);
	fprintf(out,"\n# definice regionu\n");
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",1+i*100,10001+i*100,10002+i*100,10033+i*100,10022+i*100,(i-1)*100+1,i*100+1);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",2+i*100,10033+i*100,10003+i*100,10034+i*100,10021+i*100,(i-1)*100+2,i*100+2);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",3+i*100,10034+i*100,10045+i*100,10035+i*100,10020+i*100, 3+(i-1)*100,3+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",4+i*100,10004+i*100,10047+i*100,10024+i*100,10045+i*100, 4+(i-1)*100,4+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld  hexa map yes equidistant size def\n",5+i*100,10005+i*100,10037+i*100,10025+i*100,10047+i*100,5+(i-1)*100,5+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld  hexa map yes equidistant size def\n",6+i*100,10035+i*100,10023+i*100,10036+i*100,10019+i*100,6+(i-1)*100,6+i*100);
	fprintf(out,"region %ld boundary surface %ld %ld -%ld -%ld %ld -%ld  hexa map yes equidistant size def\n",7+i*100,10038+i*100,10026+i*100,10037+i*100,10006+i*100,7+(i-1)*100,7+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld %ld %ld -%ld hexa map yes equidistant size def\n",8+i*100,10038+i*100,10007+i*100,10039+i*100,10027+i*100,8+(i-1)*100,8+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",9+i*100,10036+i*100,10032+i*100,10041+i*100,10018+i*100,9+(i-1)*100,9+i*100 );
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld %ld %ld -%ld hexa map yes equidistant size def\n",10+i*100,10039+i*100,10008+i*100,10040+i*100,10028+i*100,10+(i-1)*100,10+i*100);
	fprintf(out,"region %ld boundary surface %ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",11+i*100,10029+i*100,10040+i*100,10009+i*100,10048+i*100,11+(i-1)*100,11+i*100);
	fprintf(out,"region %ld boundary surface %ld %ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",12+i*100,10030+i*100,10048+i*100,10010+i*100,10046+i*100,12+(i-1)*100,12+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",13+i*100,10041+i*100,10031+i*100,10042+i*100,10017+i*100,13+(i-1)*100,13+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",14+i*100,10042+i*100,10046+i*100,10043+i*100,10016+i*100,14+(i-1)*100,14+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",15+i*100,10043+i*100,10011+i*100,10044+i*100,10015+i*100,15+(i-1)*100,15+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",16+i*100,10044+i*100,10012+i*100,10013+i*100,10014+i*100,16+(i-1)*100,16+i*100);
    }
    if(typecut[i] == 6){
      fprintf(out,"# %ld rez - typ %ld\n",i+1,typecut[i]);
	for(j = 0; j < 32; j++){
	  ver=j+1+i*100;
	  fprintf(out,"vertex %ld xyz %lf %lf %lf\n",ver,x[i][j],y[i][j],z[i][j]);
	}
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",1+i*100,1+i*100,2+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  2+i*100,  2+i*100,  3+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  3+i*100,  3+i*100,  4+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  4+i*100,  4+i*100,  5+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  5+i*100,  5+i*100,  6+i*100,cscount[7]); 
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  6+i*100,  6+i*100,  7+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  7+i*100,  7+i*100,  8+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  8+i*100,  8+i*100,  9+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",  9+i*100,  9+i*100, 10+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 10+i*100, 10+i*100,11+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 11+i*100, 11+i*100,12+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 12+i*100, 12+i*100,13+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 13+i*100, 13+i*100,14+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 14+i*100, 14+i*100,15+i*100,cscount[4]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 15+i*100, 15+i*100,16+i*100,cscount[3]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 16+i*100, 16+i*100,17+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 17+i*100, 17+i*100,18+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 18+i*100, 18+i*100,19+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 19+i*100, 19+i*100,20+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 20+i*100, 20+i*100,21+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 21+i*100, 21+i*100,22+i*100,cscount[2]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 22+i*100, 22+i*100, 1+i*100,cscount[1]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 23+i*100, 23+i*100,24+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 24+i*100, 24+i*100,25+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 25+i*100, 25+i*100,26+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 26+i*100, 26+i*100,27+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 27+i*100, 27+i*100,28+i*100,cscount[5]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 28+i*100, 28+i*100,29+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 29+i*100, 29+i*100,30+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 30+i*100, 30+i*100,31+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 31+i*100, 31+i*100,32+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 32+i*100, 32+i*100,23+i*100,cscount[5]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 33+i*100,  3+i*100, 22+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 34+i*100,  4+i*100, 21+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 35+i*100, 24+i*100, 20+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 36+i*100, 23+i*100, 19+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 37+i*100,  6+i*100, 26+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 38+i*100,  7+i*100, 27+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 39+i*100,  8+i*100, 28+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 40+i*100,  9+i*100, 29+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 41+i*100, 32+i*100, 18+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 42+i*100, 31+i*100, 17+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 43+i*100, 11+i*100, 16+i*100,cscount[0]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 44+i*100, 12+i*100, 15+i*100,cscount[0]);
	fprintf(out,"\n");
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 45+i*100,  4+i*100, 24+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 46+i*100, 31+i*100, 11+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 47+i*100,  5+i*100, 25+i*100,cscount[6]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n", 48+i*100, 30+i*100, 10+i*100,cscount[6]);
	fprintf(out,"\n");
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  1+i*100, 1+i*100, 2+i*100,33+i*100,22+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  2+i*100,33+i*100, 3+i*100,34+i*100,21+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  3+i*100,20+i*100,34+i*100,45+i*100,35+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  4+i*100,45+i*100, 4+i*100,47+i*100,24+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  5+i*100,47+i*100, 5+i*100,37+i*100,25+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  6+i*100,35+i*100,23+i*100,36+i*100,19+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  7+i*100,37+i*100, 6+i*100,38+i*100,26+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  8+i*100,38+i*100, 7+i*100,39+i*100,27+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",  9+i*100,36+i*100,32+i*100,41+i*100,18+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 10+i*100,39+i*100, 8+i*100,40+i*100,28+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 11+i*100,40+i*100, 9+i*100,48+i*100,29+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 12+i*100,48+i*100,10+i*100,46+i*100,30+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 13+i*100,41+i*100,31+i*100,42+i*100,17+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 14+i*100,46+i*100,43+i*100,16+i*100,42+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 15+i*100,43+i*100,11+i*100,44+i*100,15+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n", 16+i*100,44+i*100,12+i*100,13+i*100,14+i*100);
	fprintf(out,"\n");
	fprintf(out,"vertex %ld xyz %lf %lf %lf\n",33+i*100,x[i][24],y[i][22],z[i][24]);
	fprintf(out,"vertex %ld xyz %lf %lf %lf\n",34+i*100,x[i][29],y[i][31],z[i][29]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",49+i*100,23+i*100,33+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",50+i*100,33+i*100,25+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",51+i*100,33+i*100,27+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",52+i*100,32+i*100,34+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",53+i*100,34+i*100,30+i*100,cscount[8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",54+i*100,34+i*100,28+i*100,cscount[7]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",55+i*100,33+i*100,34+i*100,cscount[5]);
	fprintf(out,"\n");
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",17+i*100,23+i*100,24+i*100,50+i*100,49+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",18+i*100,25+i*100,26+i*100,51+i*100,50+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",19+i*100,32+i*100,49+i*100,55+i*100,52+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",20+i*100,55+i*100,51+i*100,27+i*100,54+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",21+i*100,52+i*100,53+i*100,30+i*100,31+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",22+i*100,53+i*100,54+i*100,28+i*100,29+i*100);
	fprintf(out,"\n");
	fprintf(out,"\n");
	fprintf(out,"# spojeni rezu %ld (typ %ld) a %ld (typ %ld)\n",i,typecut[i-1],i+1,typecut[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10001+i*100,1+(i-1)*100,1+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10002+i*100,2+(i-1)*100,2+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10003+i*100,3+(i-1)*100,3+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10004+i*100,4+(i-1)*100,4+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10005+i*100,5+(i-1)*100,5+i*100,lcount[i]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10006+i*100,6+(i-1)*100,6+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][5]+(x[i][5]-x[i-1][5])/2,y[i-1][5]+(y[i][5]-y[i-1][5])/2,z[i][5]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10007+i*100,7+(i-1)*100,7+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][6]+(x[i][6]-x[i-1][6])/2,y[i-1][6]+(y[i][6]-y[i-1][6])/2,z[i][6]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10008+i*100,8+(i-1)*100,8+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][7]+(x[i][7]-x[i-1][7])/2,y[i-1][7]+(y[i][7]-y[i-1][7])/2,z[i][7]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10009+i*100,9+(i-1)*100,9+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][8]+(x[i][8]-x[i-1][8])/2,y[i-1][8]+(y[i][8]-y[i-1][8])/2,z[i][8]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10010+i*100,10+(i-1)*100,10+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10011+i*100,11+(i-1)*100,11+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10012+i*100,12+(i-1)*100,12+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10013+i*100,13+(i-1)*100,13+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10014+i*100,14+(i-1)*100,14+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10015+i*100,15+(i-1)*100,15+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10016+i*100,16+(i-1)*100,16+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10017+i*100,17+(i-1)*100,17+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10018+i*100,18+(i-1)*100,18+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10019+i*100,19+(i-1)*100,19+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10020+i*100,20+(i-1)*100,20+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10021+i*100,21+(i-1)*100,21+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10022+i*100,22+(i-1)*100,22+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10023+i*100,23+(i-1)*100,23+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10024+i*100,24+(i-1)*100,24+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10025+i*100,25+(i-1)*100,25+i*100,lcount[i]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10026+i*100,26+(i-1)*100,26+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][25]+(x[i][25]-x[i-1][25])/2,y[i-1][25]+(y[i][25]-y[i-1][25])/2,z[i][25]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10027+i*100,27+(i-1)*100,27+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][26]+(x[i][26]-x[i-1][26])/2,y[i-1][26]+(y[i][26]-y[i-1][26])/2,z[i][26]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10028+i*100,28+(i-1)*100,28+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][27]+(x[i][27]-x[i-1][27])/2,y[i-1][27]+(y[i][27]-y[i-1][27])/2,z[i][27]);
	fprintf(out,"curve %ld order 3 vertex %ld %ld count %ld\n",10029+i*100,29+(i-1)*100,29+i*100,lcount[i]);
	fprintf(out,"polygon 1 xyz %lf %lf %lf weight 1\n",x[i-1][28]+(x[i][28]-x[i-1][28])/2,y[i-1][28]+(y[i][28]-y[i-1][28])/2,z[i][28]);      
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10030+i*100,30+(i-1)*100,30+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10031+i*100,31+(i-1)*100,31+i*100,lcount[i]);
	fprintf(out,"curve %ld vertex %ld %ld count %ld\n",10032+i*100,32+(i-1)*100,32+i*100,lcount[i]);
	for(j = 0; j < 32; j++){
	  if(j != 21 && j != 31){
	    fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,j+10002+i*100,j+1+i*100,j+10001+i*100);
	  }
	  else{
	    if(j == 21){
	      fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,10001+i*100,j+1+i*100,j+10001+i*100);
	    }
	    if(j == 31){
	      fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",j+10001+i*100,j+1+(i-1)*100,10023+i*100,j+1+i*100,j+10001+i*100);
	    }
	  }
	}
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10033+i*100,(i-1)*100+33,10003+i*100,i*100+33,10022+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10034+i*100,(i-1)*100+34,10004+i*100,i*100+34,10021+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10035+i*100,(i-1)*100+35,10024+i*100,i*100+35,10020+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10036+i*100,(i-1)*100+36,10023+i*100,i*100+36,10019+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10037+i*100,(i-1)*100+37,10006+i*100,i*100+37,10026+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10038+i*100,(i-1)*100+38,10007+i*100,i*100+38,10027+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10039+i*100,(i-1)*100+39,10008+i*100,i*100+39,10028+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10040+i*100,(i-1)*100+40,10009+i*100,i*100+40,10029+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10041+i*100,(i-1)*100+41,10032+i*100,i*100+41,10018+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10042+i*100,(i-1)*100+42,10031+i*100,i*100+42,10017+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10043+i*100,(i-1)*100+43,10011+i*100,i*100+43,10016+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10044+i*100,(i-1)*100+44,10012+i*100,i*100+44,10015+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10045+i*100,(i-1)*100+45,10004+i*100,i*100+45,10024+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10046+i*100,(i-1)*100+46,10031+i*100,i*100+46,10011+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10047+i*100,(i-1)*100+47,10005+i*100,i*100+47,10025+i*100);
	fprintf(out,"surface %ld curve %ld %ld %ld %ld\n",10048+i*100,(i-1)*100+48,10030+i*100,i*100+48,10010+i*100);
	fprintf(out,"\n# definice regionu\n");
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",1+i*100,10001+i*100,10002+i*100,10033+i*100,10022+i*100,(i-1)*100+1,i*100+1);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",2+i*100,10033+i*100,10003+i*100,10034+i*100,10021+i*100,(i-1)*100+2,i*100+2);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",3+i*100,10034+i*100,10045+i*100,10035+i*100,10020+i*100, 3+(i-1)*100,3+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",4+i*100,10004+i*100,10047+i*100,10024+i*100,10045+i*100, 4+(i-1)*100,4+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld  hexa map yes equidistant size def\n",5+i*100,10005+i*100,10037+i*100,10025+i*100,10047+i*100,5+(i-1)*100,5+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld  hexa map yes equidistant size def\n",6+i*100,10035+i*100,10023+i*100,10036+i*100,10019+i*100,6+(i-1)*100,6+i*100);
	fprintf(out,"region %ld boundary surface %ld %ld -%ld -%ld %ld -%ld  hexa map yes equidistant size def\n",7+i*100,10038+i*100,10026+i*100,10037+i*100,10006+i*100,7+(i-1)*100,7+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld %ld %ld -%ld hexa map yes equidistant size def\n",8+i*100,10038+i*100,10007+i*100,10039+i*100,10027+i*100,8+(i-1)*100,8+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",9+i*100,10036+i*100,10032+i*100,10041+i*100,10018+i*100,9+(i-1)*100,9+i*100 );
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld %ld %ld -%ld hexa map yes equidistant size def\n",10+i*100,10039+i*100,10008+i*100,10040+i*100,10028+i*100,10+(i-1)*100,10+i*100);
	fprintf(out,"region %ld boundary surface %ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",11+i*100,10029+i*100,10040+i*100,10009+i*100,10048+i*100,11+(i-1)*100,11+i*100);
	fprintf(out,"region %ld boundary surface %ld %ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",12+i*100,10030+i*100,10048+i*100,10010+i*100,10046+i*100,12+(i-1)*100,12+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",13+i*100,10041+i*100,10031+i*100,10042+i*100,10017+i*100,13+(i-1)*100,13+i*100);
	fprintf(out,"region %ld boundary surface -%ld %ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",14+i*100,10042+i*100,10046+i*100,10043+i*100,10016+i*100,14+(i-1)*100,14+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld %ld -%ld %ld -%ld hexa map yes equidistant size def\n",15+i*100,10043+i*100,10011+i*100,10044+i*100,10015+i*100,15+(i-1)*100,15+i*100);
	fprintf(out,"region %ld boundary surface -%ld -%ld -%ld -%ld %ld -%ld hexa map yes equidistant size def\n",16+i*100,10044+i*100,10012+i*100,10013+i*100,10014+i*100,16+(i-1)*100,16+i*100);
    }
  }
}
