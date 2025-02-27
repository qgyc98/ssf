#if !defined ( __SADE__ )
#define        __SADE__

#if !defined ( __general_h_ )
#include "general.h"
#endif

#if !defined ( __OBJECTIVE_FUNCTION__ )
#include "obj_funct.h"
#endif

#include <math.h>
#include <string.h>

#include "randy.h"
#include "matrix.h"

#if defined ( __HAVE_KRESLITKO__ )

#include "kreslitko.h"

#define CloseSignal 101
#define NextSignal 102
#define StartStopSignal 103
#define PostScriptSignal 104

#endif

#if defined ( __PARALLEL__ )
#include "my_mpi.h"
#endif

// Will need deklaration of randy somewhere
extern randy RND ;

class sade
/** This file declares optimization algorithm SADE, which is genetic algorithm
 working on real domains.
 created: 7.11.2003, Anicka Kucerova, anicka@cml.fsv.cvut.cz
 */
{
    public:
        sade ( void );
        ~sade ( void );
#if defined ( __PARALLEL__ )
        void run_master ( my_mpi *oM, double *assessment=NULL ) ;
	void run_slave ( my_mpi *oM );
        void CONFIGURATION_SLAVE ( void ) ;
	void EVALUATE_SLAVE ( void ) ;
	void EVALUATE_MASTER ( int start_id ) ;
        int to_continue_slave ( void ) ;
#endif
        void run ( double *assessment=NULL );
        long vilog;  ///<enable visualisation and logging
        long fitness_calls_limit;  ///<maximal numbre of calculation of objective function

    private:
// -------------------------------------- SADE TECHNOLOGY --------------------------------------------------
        void new_point ( double* p ); 
        void configuration ( void );  
        void EVALUATE_GENERATION ( long start ); 
        void SELECT ( void ); 
        void MUTATE ( void ); 
        void LOCAL_MUTATE ( void ); 
        void CROSS ( void ); 
        void FIRST_GENERATION ( void); 
        void clear_pool ( void ); 
// -------------------------------------- STOPPING CRITERIA ------------------------------------------------
        long to_continue (); 

//--------------------------------------- VISUALISATION AND LOGGING ----------------------------------------
        void init_vilog ( void ); 
        long vilog_news ( void ); 
        void vilog_results ();
        void vilog_make_out ();
        void close_vilog ( void );
#if defined ( __HAVE_KRESLITKO__ )
        kreslitko *Ko ;
        long running ;
        void open_graphics ( void );
        void draw_chromos ( double **och , long on );
        void draw_top ( void );
        long until_next_step ( void );
        void close_graphics ( void );
        double xmin, xmax, dx, ymin, ymax, dy;
#endif

//------------------variables of SADE:
    public:
        long pool_rate;
        double radioactivity, local_radioactivity;
        double mutation_rate, mutagen_rate;
        double cross_rate;
        long gradient;
        obj_funct *F;
        randy RND;

    private:
        long ActualSize, generation, fitness_call, PoolSize, SelectedSize;
        double *Force, *mutagen;
        long *origine, btg_or;
        matrix CH;

        double *bsf, bsf_value, btg_value;
        int btg;
        long bsf_birth;

	long ceraf ;

#if defined ( __PARALLEL__ )
	int mpi_mode ;
	my_mpi *M ;
#endif

};

#endif // __SADE__

