#if !defined ( __CML_FILE__ )
#define        __CML_FILE__

#if !defined ( __GENERIC_CONSTANTS__ )
#define        __GENERIC_CONSTANTS__
 
#define MinDouble -1.0e+99
#define MaxDouble +1.0e+99
#define NearZero +1.0e-16
#define NoDouble -7.87813047009E+177
#define NoInt -9279187

 
#endif

#include <stdio.h>

typedef FILE* PFILE ;
struct keyword ;
struct section ;

#define LABEL label
#define SECTION section_name

class cmlfile
{
    public:
        cmlfile ( void ) ;
        cmlfile ( const char *ofname ) ;
       ~cmlfile ( void ) ;
        // modified by Matej
    public:
        void set_names ( long onames ) ;
        void set_name_string ( long oname , char *oname_string ) ;
        void get_value_from_name ( long olabel , long &ovalue ) ;
        void get_value_from_name ( long olabel , long oorder , long &ovalue ) ;
	void get_value_with_index_from_name ( long olabel , long oindex , long &ovalue ) ;
	void get_value_with_index_from_name ( long olabel , long oindex , long oorder , long &ovalue ) ;
	void get_section_value ( long olabel, long oline, long oorder, long &ovalue ) ;
	void get_section_value ( long olabel, long oline, long oorder, double &ovalue ) ;
	void get_section_value ( long olabel, long oline, long oorder, char *ovalue ) ;
	long name ( const char *obuffer ) ;
    public:
        keyword *get_name ;
	long names ;
	// end of modification
	// start modification by smilauer 20.11.2006
    public:
	void get_next_line_in_section ( long olabel, long &ovalue );
	void get_next_line_in_section ( long olabel, double &ovalue );
	void get_next_line_in_section ( long olabel, char *ovalue );
	// end of modification
    public:
        void open ( const char *ofname ) ;
        void open_for_output ( const char *ofname ) ;
        void close ( void ) ;
        void set_labels ( long olabels ) ;
        void set_sections ( long osections ) ;
        void set_label_string ( long olabel , const char *olabel_string ) ;
        void set_section_string ( long olabel , const char *osection_string ) ;
        void set_double_format ( const char *odouble_format ) ;
        void compile ( void ) ;
        void check_requirements ( void ) ;
        void require ( long olabel ) ;
        void minimal ( long olabel , long ominimal ) ;
        void optionalize ( long olabel ) ;
        void refuse ( long olabel ) ;
        void check_label ( long olabel ) ;
        void set_minimal_rows ( long osection , long ominimal ) ;
        void check_section ( long osection ) ;
        void load_labels ( const char *orcname ) ;
        void inherit_labels ( cmlfile &of ) ;
        void get_value_count ( long olabel , long &ocount ) ;
        void get_value ( long olabel , long &ovalue ) ;
        void get_value ( long olabel , double &ovalue ) ;
        void get_value ( long olabel , char *ovalue ) ;
        void get_value ( long olabel , long oorder , long &ovalue ) ;
        void get_value ( long olabel , long oorder , double &ovalue ) ;
        void get_value ( long olabel , long oorder , char *ovalue ) ;
        void get_value_with_key ( long olabel , long okey , long &ovalue ) ;
        void get_value_with_key ( long olabel , long okey , double &ovalue ) ;
        void get_value_with_key ( long olabel , long okey , char *ovalue ) ;
        void get_value_with_key ( long olabel , long okey , long oorder , long &ovalue ) ;
        void get_value_with_key ( long olabel , long okey , long oorder , double &ovalue ) ;
        void get_value_with_key ( long olabel , long okey , long oorder , char *ovalue ) ;
        void get_value_with_key ( long olabel , const char *okey , long &ovalue ) ;
        void get_value_with_key ( long olabel , const char *okey , double &ovalue ) ;
        void get_value_with_key ( long olabel , const char *okey , char *ovalue ) ;
        void get_value_with_key ( long olabel , const char *okey , long oorder , long &ovalue ) ;
        void get_value_with_key ( long olabel , const char *okey , long oorder , double &ovalue ) ;
        void get_value_with_key ( long olabel , const char *okey , long oorder , char *ovalue ) ;
        void get_value_with_index ( long olabel , long oindex , long &ovalue ) ;
        void get_value_with_index ( long olabel , long oindex , double &ovalue ) ;
        void get_value_with_index ( long olabel , long oindex , char *ovalue ) ;
        void get_value_with_index ( long olabel , long oindex , long oorder , long &ovalue ) ;
        void get_value_with_index ( long olabel , long oindex , long oorder , double &ovalue ) ;
        void get_value_with_index ( long olabel , long oindex , long oorder , char *ovalue ) ;
        void get_percentage ( long olabel , double &ovalue ) ;
        void get_percentage ( long olabel , long oorder , double &ovalue ) ;
        void get_percentage_with_key ( long olabel , long okey , double &ovalue ) ;
        void get_percentage_with_key ( long olabel , long okey , long oorder , double &ovalue ) ;
        void get_percentage_with_key ( long olabel , const char *okey , double &ovalue ) ;
        void get_percentage_with_key ( long olabel , const char *okey , long oorder , double &ovalue ) ;
        void get_percentage_with_index ( long olabel , long okey , double &ovalue ) ;
        void get_percentage_with_index ( long olabel , long okey , long oorder , double &ovalue ) ;
        void get_line ( long olabel , char *oline ) ;
        void get_line_with_key ( long olabel , long okey , char *oline ) ;
        void get_line_with_key ( long olabel , const char *okey , char *oline ) ;
        void get_line_with_index ( long olabel , long oindex , char *oline ) ;
        void find_section ( long olabel , long &osection_rows ) ;
        void get_section_line ( long olabel , char *oline ) ;
        void end_reading_section ( void ) ;
        void out ( long olabel , long ovalue ) ;
        void out ( long olabel , double ovalue ) ;
        void out ( long olabel , const char *ovalue ) ;
        void out_label ( long olabel ) ;
        void out_long ( long ovalue ) ;
        void out_double ( double ovalue ) ;
        void out_string ( const char *ovalue ) ;
        void out_free_string ( const char *ovalue ) ;
        void out_eoln ( void ) ;
        void out_section ( long olabel ) ;
        void out_data ( const char *oline ) ;
        void out_section_end ( void ) ;
        void clear_errors ( void ) ;
        long error_opening_file ( void ) ;
        long error_reading_file ( void ) ;
        long error_writing_file ( void ) ;
        long error_in_requirements ( void ) ;
        long any_error ( void ) ;
        long label ( const char *obuffer ) ;
        long section_name ( const char *obuffer ) ;
    private:
        void open_next_input ( const char *ofname ) ;
        void stat_section ( void ) ;
        void jump_over_section ( void ) ;
        PFILE *input ;
        FILE *output ;
        keyword *stop ;
        section *block ;
        char *buffer,*double_format ;
        long files,labels,sections,actual ;
        long normal_read_disable,requirement_insufficient,compiled ;
        long open_error,read_error,write_error ;
} ;

typedef cmlfile* Pcmlfile ;

#endif

