#include <stdlib.h>

#if !defined ( __SADE__ )
#include "sade.h"
#endif

//#include "../VILOG/vilog.c"

#if defined ( __PARALLEL__ )
#include "parallel.c"
#endif

sade::sade ( void )
{
//------------------variables of SADE:

  //pool_rate = 25;
  // radioactivity = 0.05;
  //local_radioactivity = 0.05;
  //gradient = 0;
  //cross_rate = 0.1;
  //mutation_rate=0.5;
  //mutagen_rate=200;

    ActualSize=0;
    generation=0;
    fitness_call=0;
    PoolSize=0;
    SelectedSize=0;
    Force=NULL;
    mutagen=NULL;
    bsf=NULL;
    bsf_value=MinDouble;
    bsf_birth=0;
    btg=0;
    btg_value=MinDouble;
}

sade::~sade ( void )
{
}

// -------------------------------------- SADE TECHNOLOGY --------------------------------------------------

void sade::new_point ( double* p )
{
    int i;

    for ( i=0; i<F->Dim; i++ )
    {
        p[i] = RND.give_double ( F->Domain[i][0], F->Domain[i][1] );
    }
}

void sade::configuration ( void )
{
  fprintf ( stderr, "\n SADE algorithm is used: \n\n" ) ;
  fprintf ( stderr, " pool_rate: %ld \n", pool_rate ) ;
  fprintf ( stderr, " radioactivity: %f \n", radioactivity ) ;
  fprintf ( stderr, " local_radioactivity: %f \n", local_radioactivity ) ;
  fprintf ( stderr, " gradient: %ld \n", gradient ) ;
  fprintf ( stderr, " cross_rate: %f \n", cross_rate ) ;
  fprintf ( stderr, " mutation_rate: %f \n", mutation_rate ) ;
  fprintf ( stderr, " mutagen_rate: %f \n", mutagen_rate ) ;

  if ( !F ) {
      fprintf (stderr,"\n\n You must allocate obj. function first! Exit.");
      fprintf (stderr,"\n sade::configuration (file %s, line %d).\n",
	       __FILE__,__LINE__);
      exit ( 1 ) ;
  }    
    
    PoolSize = 2*pool_rate*F->Dim;
    SelectedSize = pool_rate*F->Dim;

    long i;
    mutagen = new double [F->Dim];
    for ( i=0; i<F->Dim; i++ )
    {
        mutagen[i] = (F->Domain[i][1]-F->Domain[i][0])/mutagen_rate;
        if ( F->Domain[i][2]>mutagen[i] ) mutagen[i] = F->Domain[i][2];
    }
    RND.init();
}

#if !defined ( __PARALLEL__ )
void sade::EVALUATE_GENERATION ( long start )
{
    long i, k;


    for ( i=start*SelectedSize; i<ActualSize; i++ )
    {
        if ( F->Return_to_domain )
        {
            for ( k=0; k<F->Dim; k++ )
            {
                if ( CH[i][k]<F->Domain[k][0] ) CH[i][k]=F->Domain[k][0];
                if ( CH[i][k]>F->Domain[k][1] ) CH[i][k]=F->Domain[k][1];
            }
        }
        Force[i] = F->value ( CH[i] );
        fitness_call++;
        if ( Force[i]>btg_value )
        {
            btg_value = Force[i];
            btg=i;
            btg_or = origine[i];
        }
    }
    if ( btg_value>bsf_value )
    {
        memcpy ( bsf, CH[btg], F->Dim*sizeof ( double ) );
        bsf_value = btg_value;
        bsf_birth = fitness_call;
        //F->evaluate( bsf, fitness_call, bsf_value );
    }
}	
#endif

void sade::SELECT ( void )
{
    long i1=0, i2=0, dead, last;

    while ( ActualSize > SelectedSize )
    {
        i1 = RND.give_long ( 0,ActualSize-1 );
        i2 = RND.give_long ( 1,ActualSize-1 );
        if ( i1==i2 ) i2--;
        if ( Force[i1] >= Force[i2] ) dead = i2;
        else dead = i1;
        last = ActualSize-1;
        memcpy( CH[dead], CH[last], F->Dim*sizeof(double) );
	if ( btg==last ) btg=dead;
        Force[dead] = Force[last];
        origine[dead] = origine[last];
        ActualSize--;
    }
}	

void sade::MUTATE ( void )
{
    double p=0.0, *x;
    long i, j, index=0;

    for ( i=0; i<SelectedSize; i++ )
    {
        if ( ActualSize == PoolSize ) break;
        p = RND.give_double ( 0.0, 1.0 );
        if ( p<=radioactivity )
        {
            index = RND.give_long ( 0,SelectedSize-1 );
            x = new double[F->Dim] ;
	    new_point ( x ) ;
            for ( j=0;  j<F->Dim; j++ )
            {
                CH[ActualSize][j] = CH[index][j]+mutation_rate*( x[j]-CH[index][j] );
            }
            delete [] x;
            origine[ActualSize] = 2;
            ActualSize++;
        }
    }
}

void sade::LOCAL_MUTATE ( void )
{
    double p=0.0,dCH=0.0 ;
    long i,j,index=0 ;
    for ( i=0; i<SelectedSize; i++ )
    {
        if ( ActualSize == PoolSize ) break;
        p = RND.give_double( 0.0, 1.0 );
        if ( p <= local_radioactivity )
        {
            index = RND.give_long( 0, SelectedSize-1 );

            for ( j=0; j<F->Dim; j++ )
            {
                dCH = RND.give_double( -mutagen[j],mutagen[j] );
                CH[ActualSize][j] = CH[index][j]+dCH;
            }
            origine[ActualSize] = 3;
            ActualSize++;
        }
    }
}

void sade::CROSS ( void )
{
    long i1=0, i2=0, i3=0, j=0 ;
    long direction = 1;
    while ( ActualSize < PoolSize )
    {
        i1 = RND.give_long( 0,SelectedSize-1 );
        i2 = RND.give_long( 1,SelectedSize-1 );
        if ( i1==i2 ) i2--;
        i3 = RND.give_long( 0,SelectedSize-1 );
        if ( gradient && (Force[i1]>Force[i2]) ) direction=-1;

        for ( j=0 ; j<F->Dim ; j++ )
            CH[ActualSize][j] = CH[i3][j]+cross_rate*(double)(direction)*( CH[i2][j]-CH[i1][j] ) ;
        origine[ActualSize] = 1;
        ActualSize++;
    }
}

void sade::FIRST_GENERATION ( void)
{
    long i;

    Force = new double [PoolSize];
    allocm ( PoolSize, F->Dim, CH ) ;
    origine = new long [PoolSize];
    for ( i=0; i<PoolSize; i++ )
    {
        new_point ( CH[i] ) ;
        origine[i] = 0;
    }
    bsf = new double [F->Dim];
    ActualSize = PoolSize;
    EVALUATE_GENERATION ( 0 );
    SELECT ();
    //vilog_news ();
}	

void sade::clear_pool ( void )
{
    if ( Force ) delete [] Force;
    if ( bsf ) delete [] bsf;
    if ( mutagen ) delete [] mutagen;
    if ( origine ) delete [] origine;
}

// -------------------------------------- STOPPING CRITERIA ------------------------------------------------

long sade::to_continue ()
{
    if ( fitness_call > fitness_calls_limit ) return 0;
    if ( F->optimum ) if ( *F->optimum-bsf_value<=F->precision ) return 0;
    return 1;
}

// -------------------------------------- MAIN -------------------------------------------------------------

void sade::run ( double *assessment )
{
    long stop=0 ;

    ceraf=1;
    configuration ();
    //init_vilog ();

    FIRST_GENERATION ();

    while ( !stop )
    {
        MUTATE ();
        LOCAL_MUTATE ();
        CROSS ();
        EVALUATE_GENERATION ( 1 );
        SELECT ();
        //if ( !vilog_news () ) stop=1;
        if ( !to_continue () ) stop=1;
        generation++;
    }
    //vilog_results ();
    //vilog_make_out ();

    F->user_evaluate ( bsf );
    if ( assessment ) *assessment = bsf_value/(double)fitness_call;

    //close_vilog ();
    clear_pool ();
}
