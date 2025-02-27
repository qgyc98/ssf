#include <stdlib.h>

long max(long i,long j)
{
  if (i>j) return i;
  else     return j;
}


long min(long i,long j)
{
  if (j>i) return i;
  else     return j;
}

long bandwidth (long *adr,long n)
/*
  funkce zjistuje sirku pasu ve skyline
  
  adr - pole adres diagonalnich prvku
  n - pocet slozek pole adr
  
  5.8.1998
*/
{
  long i,j,k;
  
  j=0;
  for (i=0;i<n;i++){
    k=adr[i+1]-adr[i];
    if (k>j)  j=k;
  }
  return j;
}

void eliminuj_4i_rev(double *a,long n,long s)
{
  double sum1,sum2,sum3,sum4,x;
  long i,j,k;
  
  for(j=0;j<n;j++)
  {
  	sum1=sum2=0.0;  	
  	for(k=0;k<(j-1);k+=2)
  	{
  	  sum1+=(a[j*s+k]*a[j*s+k]);
  	  sum2+=(a[j*s+k+1]*a[j*s+k+1]);
  	}    	
  	if (j&1) sum1+=(a[j*s+j-1]*a[j*s+j-1]);
  	//if ((a[j*s+j]-sum1-sum2)<=0) { printf ("Chyba \n"); }
  	a[j*s+j]=sqrt(a[j*s+j]-sum1-sum2);
  if(j&1)
  {	
  	for(i=j+1;i<(n-3);i+=4)
  	{
  	  sum1=sum2=sum3=sum4=0.0;  	  
  	  for(k=0;k<j;k++)
  	  {
  	    x=a[j*s+k];		
  	    sum1+=a[i*s+k]*x;
  	    sum2+=a[(i+1)*s+k]*x;	
  	    sum3+=a[(i+2)*s+k]*x;	
  	    sum4+=a[(i+3)*s+k]*x;	
  	  }	
  	  x=a[j*s+j];
  	  a[i*s+j]=(a[i*s+j]-sum1)/x;
  	  a[(i+1)*s+j]=(a[(i+1)*s+j]-sum2)/x;  	  
  	  a[(i+2)*s+j]=(a[(i+2)*s+j]-sum3)/x;  	  
  	  a[(i+3)*s+j]=(a[(i+3)*s+j]-sum4)/x;  	  
        }
        for(;i<n;i++)
        {
  	  sum1=0.0;  	  
  	  for(k=0;k<j;k++)
  	  {
  	    sum1+=a[i*s+k]*a[j*s+k];	
  	  }	
  	  a[i*s+j]=(a[i*s+j]-sum1)/a[j*s+j];        	
        }	
}
else
{
  	for(i=n-1;i>(j+3);i-=4)
  	{
  	  sum1=sum2=sum3=sum4=0.0;  	  
  	  for(k=0;k<j;k++)
  	  {
  	    x=a[j*s+k];		
  	    sum1+=a[i*s+k]*x;
  	    sum2+=a[(i-1)*s+k]*x;	
  	    sum3+=a[(i-2)*s+k]*x;	
  	    sum4+=a[(i-3)*s+k]*x;	
  	  }	
  	  x=a[j*s+j];
  	  a[i*s+j]=(a[i*s+j]-sum1)/x;
  	  a[(i-1)*s+j]=(a[(i-1)*s+j]-sum2)/x;  	  
  	  a[(i-2)*s+j]=(a[(i-2)*s+j]-sum3)/x;  	  
  	  a[(i-3)*s+j]=(a[(i-3)*s+j]-sum4)/x;  	  
        }
        for(;i>j;i--)
        {
  	  sum1=0.0;  	  
  	  for(k=0;k<j;k++)
  	  {
  	    sum1+=a[i*s+k]*a[j*s+k];	
  	  }	
  	  a[i*s+j]=(a[i*s+j]-sum1)/a[j*s+j];        	
        }	
}
    }		
}	
  

void faze1_opt_sp(double *a,double *b,long n,long s,long n2)
{
   long i,j,k;	
   double *sumy,x,sum1,sum2,sum3,sum4;	
   sumy = new double [n];       
   for(i=0;i<n;i++)
   {
   	for(j=0;j<(n-3);j+=4)
   	{
   		sum1=sum2=sum3=sum4=0.0;
   		for(k=0;k<j;k++)
   		{
   			x=sumy[k];
   			sum1+=x*a[j*s+k];//L
   			sum2+=x*a[(j+1)*s+k];//L
   			sum3+=x*a[(j+2)*s+k];//L   			
   			sum4+=x*a[(j+3)*s+k];//L
   		}
   		x=b[i*s+j]-sum1;
   		x=sumy[j]=x/a[j*s+j];//A,L	
   		//sum2+=(x*a[(j+1)*n+j])/a[j*n+j];
   		sum2+=x*a[(j+1)*s+j];
   		sum3+=x*a[(j+2)*s+j];
   		sum4+=x*a[(j+3)*s+j];
   		x=b[i*s+j+1]-sum2;
   		x=sumy[j+1]=x/a[(j+1)*s+(j+1)];//A,L	
   		//sum3+=(x*a[(j+2)*n+j+1])/a[(j+1)*n+(j+1)];
   		sum3+=x*a[(j+2)*s+j+1];
   		sum4+=x*a[(j+3)*s+j+1];
   		x=b[i*s+j+2]-sum3;
   		sumy[j+2]=x/a[(j+2)*s+(j+2)];//A,L

   		//sum4+=(x*a[(j+3)*n+j+2])/a[(j+2)*n+(j+2)];	
   		sum4+=sumy[j+2]*a[(j+3)*s+j+2];
   		x=b[i*s+j+3]-sum4;
   		sumy[j+3]=x/a[(j+3)*s+(j+3)];//A,L	   		
   	}	
   	for(;j<n;j++)
   	{
   		sum1=sum2=0.0;
   		for(k=0;k<(j-1);k+=2)
   		{
   			sum1+=sumy[k]*a[j*s+k];//L
   			sum2+=sumy[k+1]*a[j*s+k+1];//L
   		}
   		if (j&1) sum1+=sumy[j-1]*a[j*s+j-1];//L
   		sumy[j]=(b[i*s+j]-sum1-sum2)/a[j*s+j];//A,L	
   	}	   	
   	for(k=0;k<n;k++) b[i*s+k]=sumy[k];   		
   }	
   delete [] sumy;
}

   
void faze2_rev_sp(double *a,double *b,long n,long s,long n2)
{
   double sum,x,sum1,sum2,sum3,sum4;	
   long i,j,k;	
   for(i=0;i<n2;i++)
   {
   	if (i&1)
   	{
   	for(j=0;j<=(i-3);j+=4)
   	{
   		sum1=sum2=sum3=sum4=0.0;
   		for(k=0;k<n;k++)
   		{
   			//sum+=a[(i+n2)*n+k]*a[(k+n2)*n+j];//A
   			x=b[i*s+k];
   			sum1+=x*b[j*s+k];//A
   			sum2+=x*b[(j+1)*s+k];//A
   			sum3+=x*b[(j+2)*s+k];//A
   			sum4+=x*b[(j+3)*s+k];//A
   		}
   		a[i*s+j]-=sum1;//B
   		a[i*s+j+1]-=sum2;//B
   		a[i*s+j+2]-=sum3;//B
   		a[i*s+j+3]-=sum4;//B
   	}
   	for(;j<=i;j++)
   	{
   		sum=0.0;
   		for(k=0;k<n;k++)
   		{
   			//sum+=a[(i+n2)*n+k]*a[(k+n2)*n+j];//A
   			sum+=b[i*s+k]*b[j*s+k];//A
   		}
   		a[i*s+j]-=sum;//B
   	}
	}
	else
	{
   	for(j=i;j>=3;j-=4)
   	{
   		sum1=sum2=sum3=sum4=0.0;
   		for(k=0;k<n;k++)
   		{
   			//sum+=a[(i+n2)*n+k]*a[(k+n2)*n+j];//A
   			x=b[i*s+k];
   			sum1+=x*b[j*s+k];//A
   			sum2+=x*b[(j-1)*s+k];//A
   			sum3+=x*b[(j-2)*s+k];//A
   			sum4+=x*b[(j-3)*s+k];//A
   		}
   		a[i*s+j]-=sum1;//B
   		a[i*s+j-1]-=sum2;//B
   		a[i*s+j-2]-=sum3;//B
   		a[i*s+j-3]-=sum4;//B
   	}
   	for(;j>=0;j--)
   	{
   		sum=0.0;
   		for(k=0;k<n;k++)
   		{
   			//sum+=a[(i+n2)*n+k]*a[(k+n2)*n+j];//A
   			sum+=b[i*s+k]*b[j*s+k];//A
   		}
   		a[i*s+j]-=sum;//B
   	}	
	}	   	
   }
}


void faze1(double *a,double *b,long n,long s,long n2)
{
   long i,j,k;	
   double sum1;
   for(i=0;i<n;i++)
   {
   	for(j=0;j<n;j++)
   	{
   		sum1=0.0;
   		for(k=0;k<j;k++) sum1+=b[i*s+k]*a[j*s+k];//L   		   		
   		b[i*s+j]=(b[i*s+j]-sum1)/a[j*s+j];//A,L	
   	}	
   }	
}


void faze2(double *a,double *b,long n,long s,long n2)
{
   double sum1;
   long i,j,k;	
   for(i=0;i<n2;i++)
   {
   	for(j=0;j<=i;j++)
   	{
   		sum1=0.0;
   		for(k=0;k<n;k++)  sum1+=b[i*s+k]*b[j*s+k];//A
   		a[i*s+j]-=sum1;//B
   	}
   }
}


void faze1_block(double *a,double *b,long n,long n2,long block)
{
   long i,j,j2,step1,k;	
   double *sumy,x,sum1,sum2,sum3,sum4;	
   sumy=new double [n];       
   step1=(n+block-1)/block;
   for(j2=0;j2<n;j2+=step1)
   {
   
   for(i=0;i<n2;i++)
   {
	if (j2)
	{
	   for(k=0;k<j2;k++) sumy[k]=b[i*n+k];   		
	}	
   	for(j=j2;j<min(n,j2+step1)-3;j+=4)
   	{
   		sum1=sum2=sum3=sum4=0.0;
   		for(k=0;k<j;k++)
   		{
   			x=sumy[k];
   			sum1+=x*a[j*n+k];//L
   			sum2+=x*a[(j+1)*n+k];//L
   			sum3+=x*a[(j+2)*n+k];//L   			
   			sum4+=x*a[(j+3)*n+k];//L
   		}
   		x=b[i*n+j]-sum1;
   		x=sumy[j]=x/a[j*n+j];//A,L	
   		//sum2+=(x*a[(j+1)*n+j])/a[j*n+j];
   		sum2+=x*a[(j+1)*n+j];
   		sum3+=x*a[(j+2)*n+j];
   		sum4+=x*a[(j+3)*n+j];
   		x=b[i*n+j+1]-sum2;
   		x=sumy[j+1]=x/a[(j+1)*n+(j+1)];//A,L	
   		//sum3+=(x*a[(j+2)*n+j+1])/a[(j+1)*n+(j+1)];
   		sum3+=x*a[(j+2)*n+j+1];
   		sum4+=x*a[(j+3)*n+j+1];
   		x=b[i*n+j+2]-sum3;
   		sumy[j+2]=x/a[(j+2)*n+(j+2)];//A,L
   		//sum4+=(x*a[(j+3)*n+j+2])/a[(j+2)*n+(j+2)];	
   		sum4+=sumy[j+2]*a[(j+3)*n+j+2];
   		x=b[i*n+j+3]-sum4;
   		sumy[j+3]=x/a[(j+3)*n+(j+3)];//A,L	   			   		
   	}	
   	for(;j<min(n,j2+step1);j++)
   	{
   		sum1=sum2=0.0;
   		for(k=0;k<(j-1);k+=2)
   		{
   			sum1+=sumy[k]*a[j*n+k];//L
   			sum2+=sumy[k+1]*a[j*n+k+1];//L
   		}
   		if (j&1) sum1+=sumy[j-1]*a[j*n+j-1];//L
   		sumy[j]=(b[i*n+j]-sum1-sum2)/a[j*n+j];//A,L	
   	}	   	   	
   	for(k=j2;k<min(n,j2+step1);k++) b[i*n+k]=sumy[k];   		
   }	
   }
   delete [] sumy;
}


void faze2_block(double *a,double *b,long n,long n2,long block)
{
   double sum,x,sum1,sum2,sum3,sum4;	
   long i,j,j2,k,step2;	

   step2=(n+block-1)/block;
   for(j2=0;j2<n;j2+=step2)
   {
   for(i=0;i<n2;i++)
   {
   	for(j=j2;j<=min((i-3),j2+step2-1);j+=4)
   	{
   		sum1=sum2=sum3=sum4=0.0;
   		for(k=0;k<n;k++)
   		{
   			//sum+=a[(i+n2)*n+k]*a[(k+n2)*n+j];//A
   			x=b[i*n+k];
   			sum1+=x*b[j*n+k];//A
   			sum2+=x*b[(j+1)*n+k];//A
   			sum3+=x*b[(j+2)*n+k];//A
   			sum4+=x*b[(j+3)*n+k];//A
   		}
   		a[i*n+j]-=sum1;//B
   		a[i*n+j+1]-=sum2;//B
   		a[i*n+j+2]-=sum3;//B
   		a[i*n+j+3]-=sum4;//B
   	}
   	for(;j<=min(i,j2+step2-1);j++)
   	{
   		sum=0.0;
   		for(k=0;k<n;k++)
   		{
   			//sum+=a[(i+n2)*n+k]*a[(k+n2)*n+j];//A
   			sum+=b[i*n+k]*b[j*n+k];//A
   		}
   		a[i*n+j]-=sum;//B
   	}   	
   }
   }			
}


void napln_a(long i1,long i2,long band,double *pole,double *a,long *adr)
{
  long i,j,a1,a2,kolik;
  double *x;	
  for(i=i1;i<i2;i++)	
  {
    a1=adr[i];
    a2=adr[i+1];
    kolik=min(a2-a1,i+1-i1);
    x=pole+a1;
    for(j=0;j<kolik;j++) {a[(i-i1)*band+(i-i1)-j]=*x;   x++; } 
    for(j=kolik;j<=(i-i1);j++) a[(i-i1)*band+(i-i1)-j]=0.0;
  }	

}	


void uloz_a(long i1,long i2,long band,double *pole,double *a,long *adr)
{
  long i,j,a1,a2,kolik;
  double *x;	
  for(i=i1;i<i2;i++)	
  {
    a1=adr[i];
    a2=adr[i+1];
    kolik=min(a2-a1,i+1-i1);
    x=pole+adr[i];
    for(j=0;j<kolik;j++) { *x=a[(i-i1)*band+(i-i1)-j];    x++; }
  }		

}	


void napln_b(long i1,long i2,long band,double *pole,double *b,long *adr)
{
  long i,j,a1,a2,kolik;
  double *x;	
  for(i=i1;i<i2;i++)	
  {
    a1=adr[i];
    a2=adr[i+1];
    kolik=i+1-i1;    
    x=pole+adr[i]+kolik;
    for(j=0;j<((a2-a1)-kolik);j++) { b[(i-i1)*band+(band-1)-j]=*x;    x++; }
    for(j=((a2-a1)-kolik);j<band;j++) b[(i-i1)*band+(band-1)-j]=0.0;
  }		

}


void uloz_b(long i1,long i2,long band,double *pole,double *b,long *adr)
{
  long i,j,a1,a2,kolik;
  double *x;	
  for(i=i1;i<i2;i++)	
  {
    a1=adr[i];
    a2=adr[i+1];
    //kolik=min(a2-a1,i+1-i1)+1;    
    kolik=i+1-i1;    
    x=pole+adr[i]+kolik;
    for(j=0;j<((a2-a1)-kolik);j++) { *x=b[(i-i1)*band+(band-1)-j];    x++; }
  }	

}


void blokove(double *pole, long *adr, long n, long band)
{
  long i,i1,i2;
  double *a,*b;	
  
  printf("Band = %li \n",band);
  a=new double [band*band];
  b=new double [band*band];
  
  i=0;  
  //if ((n-i)<(2*band)) krok=n-i;
  i1=i2=band;
  napln_a(i,i1,band,pole,a,adr);
  /*
  printf("*************************A**********************\n");  
  for(k=0;k<band;k++)
  {
    for(j=0;j<band;j++)  printf("  %g  ",a[j+k*band]);
    printf("\n");
  } */   

  /*
  printf("Pred elim \n");
  fflush(stdout);*/  
  eliminuj_4i_rev(a,band,band);
  /*
  printf("*************************sqrt(A)***********************\n");  
  for(k=0;k<band;k++)
  {
    for(j=0;j<band;j++)  printf("  %g  ",a[j+k*band]);
    printf("\n");
  } */       

  while(i1<n)
  {
    i2=min(i1+band,n); 	      
    //printf("Cyklus   i= %li i1= %li i2= %li\n",i,i1,i2);
    
    //napln_b(i1,i2,band,pole,b,adr);
    napln_b(i1,i2,i1-i,pole,b,adr);        
    /*
    printf("**********************************\n");
    printf("Po napln B \n");
    for(k=0;k<(i2-i1);k++)
    {
      for(j=0;j<(i2-i1);j++)  printf("  %g  ",b[j+k*(i2-i1)]);
      printf("\n");
    } */       
                
    //faze1(a,b,i1-i,i1-i,i2-i1);
    faze1_block(a,b,i1-i,i2-i1,4);
    /*
    printf("**********************************\n");
    printf("Po faze1 \n");
    for(k=0;k<(i2-i1);k++)
    {
      for(j=0;j<(i2-i1);j++)  printf("  %g  ",b[j+k*(i2-i1)]);
      printf("\n");
    } */       
        
    //uloz_a(i,i1,band,pole,a,adr);
    uloz_a(i,i1,i1-i,pole,a,adr);
    
    //napln_a(i1,i2,band,pole,a,adr);
    napln_a(i1,i2,i2-i1,pole,a,adr);
    /*
    printf("**********************************\n");
    printf("Po napln A \n");
    for(k=0;k<(i2-i1);k++)
    {
      for(j=0;j<(i2-i1);j++)  printf("  %g  ",a[j+k*(i2-i1)]);
      printf("\n");
    } */       
    
    //faze2(a,b,i1-i,i1-i,i2-i1);
    faze2(a,b,i1-i,i2-i1,4);
    /*
    printf("**********************************\n");
    printf("Po faze2 \n");
    for(k=0;k<(i2-i1);k++)
    {
      for(j=0;j<(i2-i1);j++)  printf("  %g  ",a[j+k*(i2-i1)]);
      printf("\n");
    } */           
    //uloz_b(i1,i2,band,pole,b,adr);
    uloz_b(i1,i2,i1-i,pole,b,adr);
    
    eliminuj_4i_rev(a,i2-i1,i2-i1);
    /*
    printf("**********************************\n");
    printf("Po elim A \n");
    for(k=0;k<(i2-i1);k++)
    {
      for(j=0;j<(i2-i1);j++)  printf("  %g  ",a[j+k*(i2-i1)]);
      printf("\n");
    } */       
    
    i=i1;
    i1=i2;
  }
  //uloz_a(i,n,band,pole,a,adr);
  uloz_a(i,n,i1-i,pole,a,adr);
  //tisk_sky(adr,pole,n);
  delete [] a;
  delete [] b;
}


void ldl_sky2 (double *a,double *x,double *y,long *adr,long n,long tc)
/*
  funkce resi soustavu linearnich algebraickych rovnic
  reseni se provadi rozkladem LDL
  matice soustavy je ulozena ve skylinu
  
  a - matice soustavy
  x - vektor reseni
  y - vektor prave strany
  adr - pole adres diagonalnich prvku
  n - pocet neznamych
  tc - typ vypoctu  tc=1 - provede se eliminace i zpetny chod
                    tc=2 - provede se pouze eliminace
		    tc=3 - provede se pouze zpetny chod
							
  10.7.1996
  funkce je totozna s procedurou ve fortranu, ktera je stejne
  rychla jako colsol od Batheho
  funkce byla testovana s vysledky od Batheho
	                              ldl ve fortranu
				      stare ldl v c
*/
{
  long i,k,ac,ac1,ack,ack1;
  double s,g;
  
  if (tc==1 || tc==2)
  {
    i=bandwidth(adr,n);
    blokove(a,adr,n,i);
  }
  if (tc==1 || tc==3){
    /********************/
    /*  vypocet reseni  */
    /********************/
    /*  vypocet  Lz=y => z (y se prepisuji na z) */
    for (k=1;k<n;k++){
      /*  smycka pres nezname  */
      ack=adr[k]+1;  ack1=adr[k+1];
      ac1=k-1;  s=0.0;
      for (i=ack;i<ack1;i++){
	s+=a[i]*y[ac1];
	ac1--;
      }
      y[k]-=s;
    }
    /*  deleni zbyle soustavy diagonalnimi prvky (DLx=z => Lx=1/Dz) */
    for (k=0;k<n;k++){
      ac=adr[k];  y[k]/=a[ac];
    }
    /*  vypocet Lx=1/Dz => x  */
    for (k=n-1;k>-1;k--){
      /*  smycka pres nezname  */
      ack=adr[k]+1;  ack1=adr[k+1];
      x[k]=y[k];  g=x[k];
      ac=k-1;
      for (i=ack;i<ack1;i++){
	y[ac]-=a[i]*g;
	ac--;
      }
    }
  }
}

