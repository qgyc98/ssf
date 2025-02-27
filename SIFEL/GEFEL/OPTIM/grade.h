#if !defined ( __GRADE__ )
#define        __GRADE__

#if !defined ( __general_h_ )
#include "general.h"
#endif

#if !defined ( __OBJECTIVE_FUNCTION__ )
#include "obj_funct.h"
#endif

#include <math.h>
#include <string.h>

#include <randy.h>
#include "matrix.h"

#if defined ( __HAVE_KRESLITKO__ )

#include <kreslitko.h>

#define CloseSignal 101
#define NextSignal 102
#define StartStopSignal 103
#define PostScriptSignal 104

#endif

#if defined ( __PARALLEL__ )
#include "my_mpi.h"
#endif

// Will need declaration of randy somewhere
extern randy RND ;

class grade
{
    public:
        grade ( void );
        ~grade ( void );
        void run ( double *assessment=NULL );
        long vilog;
        long vilog_gati;
        long fitness_calls_limit;
        double cross_limit;

#if defined ( __PARALLEL__ )
        void run_master ( my_mpi *oM, double *assessment=NULL ) ;
	void run_slave ( my_mpi *oM );
        void CONFIGURATION_SLAVE ( void ) ;
	void EVALUATE_SLAVE ( void ) ;
	void EVALUATE_MASTER ( int start_id ) ;
        int to_continue_slave ( void ) ;
#endif


    private:
	void new_point ( double* p );
        void configuration ( void );
        void EVALUATE_GENERATION ( long start );
        void SELECT ( void );
        void MUTATE ( void );
        void LOCAL_MUTATE ( void );
        void CROSS ( void );
        void FIRST_GENERATION ( void);
        void clear_pool ( void );
        long to_continue ();

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

};

#endif // __GRADE__
