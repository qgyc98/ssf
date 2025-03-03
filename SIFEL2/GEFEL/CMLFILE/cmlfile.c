#include "cmlfile.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BufferLength 1024
#define KeywordLength 64
#define MaxFiles 16
#define DefaultDoubleFormat "%10.10E"

#define KeywordSection "SECTION"
#define KeywordEndSection "END"
#define KeywordFile "FILE"

#define OK 1

#define None -1
#define Comment -2
#define File -3
#define Section -4
#define EndSection -5

#define Required 1
#define Optional 0
#define Refused -1

long read_to_eoln ( FILE *ofile , char *obuffer )
{
    if ( !fgets( obuffer,BufferLength,ofile ))
    {
        return None ;
    }
    long i=0 ;
    while ( obuffer[i] )
    {
        i++ ;
    }
    obuffer[i-1]=0 ;
    return OK ;
}

long str_equal ( char *opattern , char *ostring )
{
    long i=0 ;
    while ( opattern[i] )
    {
        if ( opattern[i]!=ostring[i] )
        {
            return 0 ;
        }
        i++ ;
    }
    return 1 ;
}

long read_string ( FILE *ofile , char *obuffer )
{
    long i,j ;
    j=fscanf( ofile,"%s",obuffer ) ;

    if ( j!=1 )
    {
        return j ;
    }
    if ( obuffer[0]!='"' )
    {
        return 1 ;
    }
    i=0 ;
    while ( obuffer[i] )
    {
        i++ ;
    }
    if ( obuffer[i-1]=='"' )
    {
        obuffer[i-1]=0 ;
        memcpy( obuffer,obuffer+1,i ) ;
        return 1 ;
    }
    char c ;
    do
    {
        c=fgetc( ofile ) ;
        if ( c==EOF )
        {
            return -1 ;
        }
        obuffer[i]=c ;
        i++ ;
    }
    while ( c!='"' ) ;
    obuffer[i-1]=0 ;
    memcpy( obuffer,obuffer+1,i ) ;
    return 1 ;
}

struct keyword
{
    long *file,*position ;
    long i,n,req ;
    char *word ;
    keyword ( void ) ;
   ~keyword ( void ) ;
    void copy ( keyword &ok ) ;
} ;

keyword::keyword ( void )
{
    file=NULL ;
    position=NULL ;
    word=new char[KeywordLength] ;
    i=0 ;
    n=0 ;
    req=Optional ;
}

keyword::~keyword ( void )
{
    if ( file )
    {
        delete [] file ;
    }
    if ( position )
    {
        delete [] position ;
    }
    if ( word )
    {
        delete [] word ;
    }
}

void keyword::copy ( keyword &ok )
{
    req=ok.req ;
    memcpy( word,ok.word,KeywordLength ) ;
}

struct section
{
    char *word ;
    long *pos,*rows,*file,*end ;
    long files,rows_all,actual,i,req ;
    section ( void ) ;
   ~section ( void ) ;
    void copy ( section &os ) ;
} ;

section::section ( void )
{
    pos=new long[MaxFiles] ;
    rows=new long[MaxFiles] ;
    file=new long[MaxFiles] ;
    end=new long[MaxFiles] ;
    word=new char[KeywordLength] ;
    files=0 ;
    rows_all=0 ;
    actual=0 ;
    i=0 ;
    req=0 ;
}

section::~section ( void )
{
    if ( pos )
    {
        delete [] pos ;
    }
    if ( rows )
    {
        delete [] rows ;
    }
    if ( file )
    {
        delete [] file ;
    }
    if ( end )
    {
        delete [] end ;
    }
    if ( word )
    {
        delete [] word ;
    }
}

void section::copy ( section &os )
{
    req=os.req ;
    memcpy( word,os.word,KeywordLength ) ;
}

cmlfile::cmlfile ( void )
{
    input=new PFILE[MaxFiles] ;
    long i ;
    for ( i=0 ; i<MaxFiles ; i++ )
    {
        input[i]=NULL ;
    }
    files=0 ;
    labels=0 ;
	// modified by matej
	names=0 ;
	get_name=NULL ;
	// end of modification
    sections=0 ;
    stop=NULL ;
    block=NULL ;
    output=NULL ;

    buffer=new char[BufferLength] ;
    double_format=new char[16] ;
    strcpy( double_format,DefaultDoubleFormat ) ;

    clear_errors() ;
    normal_read_disable=0 ;
    requirement_insufficient=0 ;
    compiled=0 ;
}

cmlfile::cmlfile ( char *ofname )
{
    input=new PFILE[MaxFiles] ;
    long i ;
    for ( i=0 ; i<MaxFiles ; i++ )
    {
        input[i]=NULL ;
    }
	// modified by matej
	names=0 ;
	get_name=NULL ;
	// end of modification
	labels=0 ;
    sections=0 ;
    stop=NULL ;
    block=NULL ;
    output=NULL ;

    buffer=new char[BufferLength] ;
    double_format=new char[16] ;
    strcpy( double_format,DefaultDoubleFormat ) ;

    clear_errors() ;
    normal_read_disable=0 ;
    requirement_insufficient=0 ;
    compiled=0 ;

    input[0]=fopen( ofname,"rt" ) ;
    if ( !input[0] )
    {
        printf( "Cannot open file %s.\n",ofname ) ;
        open_error=1 ;
        return ;
    }
    files=1 ;
}

cmlfile::~cmlfile ( void )
{
    long i ;
    for ( i=0 ; i<files ; i++ )
    {
        fclose( input[i] ) ;
    }
    delete [] input ;
    if ( output )
    {
        fclose( output ) ;
    }
    if ( buffer )
    {
        delete [] buffer ;
    }
    if ( double_format )
    {
        delete [] double_format ;
    }
    if ( stop )
    {
        delete [] stop ;
    }
    if ( block )
    {
        delete [] block ;
    }
//  modified by Matej
	if ( get_name )
	  {
		delete [] get_name ;
	  }
// end of modification
}

//  modified by Matej
void cmlfile::set_names ( long onames )
{
  names=onames ;
  get_name=new keyword[onames] ;
}

void cmlfile::set_name_string ( long oname , char *oname_string )
{
  strcpy( get_name[oname].word,oname_string ) ;
}

void cmlfile::get_value_from_name ( long olabel , long &ovalue )
{
  char s[KeywordLength] ;

  get_value( olabel,s ) ;
  ovalue=name( s ) ;
}

void cmlfile::get_value_from_name ( long olabel , long oorder , long &ovalue )
{
  char s[KeywordLength] ;

  get_value( olabel,oorder,s ) ;
  ovalue=name( s ) ;
}

void cmlfile::get_value_with_index_from_name ( long olabel , long oindex , long &ovalue )
{
  char s[KeywordLength] ;

  get_value_with_index( olabel, oindex, s ) ;
  ovalue=name( s ) ;
}

void cmlfile::get_value_with_index_from_name ( long olabel , long oindex , long oorder , long &ovalue )
{
  char s[KeywordLength] ;

  get_value_with_index( olabel, oindex, oorder, s ) ;
  ovalue=name( s ) ;
}

void cmlfile::get_section_value ( long olabel, long oline, long oorder, long &ovalue )
{
  if ( !compiled )
    {
      compile() ;
    }

  actual=block[olabel].file[0] ;
  fseek( input[actual],block[olabel].pos[0],SEEK_SET ) ;
  read_to_eoln( input[actual],buffer ) ;

  long i, j ;
  for ( i=0 ; ( i<oline ) && ( i<block[olabel].rows_all ) ; i++ )
    {
      read_to_eoln( input[actual],buffer ) ;
    }

  for ( i=0 ; i<oorder ; i++ )
    {
      read_string( input[actual],buffer ) ;
    }
  j=fscanf( input[actual],"%ld",&ovalue ) ;
  
  if ( j!=1 )
    {
      printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
      read_error=1 ;
    }
}

void cmlfile::get_section_value ( long olabel, long oline, long oorder, double &ovalue )
{
  if ( !compiled )
    {
      compile() ;
    }

  actual=block[olabel].file[0] ;
  fseek( input[actual],block[olabel].pos[0],SEEK_SET ) ;
  read_to_eoln( input[actual],buffer ) ;

  long i, j ;
  for ( i=0 ; ( i<oline ) && ( i<block[olabel].rows_all ) ; i++ )
    {
      read_to_eoln( input[actual],buffer ) ;
    }

  for ( i=0 ; i<oorder ; i++ )
    {
      read_string( input[actual],buffer ) ;
    }
  j=fscanf( input[actual],"%lf",&ovalue ) ;
  
  if ( j!=1 )
    {
      printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
      read_error=1 ;
    }
}

void cmlfile::get_section_value ( long olabel, long oline, long oorder, char *ovalue )
{
  if ( !compiled )
    {
      compile() ;
    }

  actual=block[olabel].file[0] ;
  fseek( input[actual],block[olabel].pos[0],SEEK_SET ) ;
  read_to_eoln( input[actual],buffer ) ;

  long i, j ;
  for ( i=0 ; ( i<oline ) && ( i<block[olabel].rows_all ) ; i++ )
    {
      read_to_eoln( input[actual],buffer ) ;
    }

  for ( i=0 ; i<oorder ; i++ )
    {
      read_string( input[actual],buffer ) ;
    }
  j=read_string( input[actual],ovalue ) ;
  
  if ( j!=1 )
    {
      printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
      read_error=1 ;
    }
}

long cmlfile::name ( char *obuffer )
{
  long i ;
  for ( i=0 ; i<names ; i++ )
    {
	  if ( str_equal( obuffer,get_name[i].word ))
        {
		  return i ;
        }
    }
  return None ;
}

// end of modification
void cmlfile::open ( char *ofname )
{
    long i ;
    for ( i=0 ; i<files ; i++ )
    {
        fclose( input[i] ) ;
    }
    input[0]=fopen( ofname,"rt" ) ;
    if ( !input[0] )
    {
        printf( "Cannot open file %s.\n",ofname ) ;
        open_error=1 ;
        return ;
    }
    files=1 ;
    compiled=0 ;
}

void cmlfile::open_for_output ( char *ofname )
{
    output=fopen( ofname,"wt" ) ;
    if ( !output )
    {
        printf( "For some reason cannot open file %s.\n",ofname ) ;
        open_error=1 ;
    }
}

void cmlfile::close ( void )
{
    long i ;
    for ( i=0 ; i<files ; i++ )
    {
        fclose( input[i] ) ;
        input[i]=NULL ;
    }
    files=0 ;
    if ( output )
    {
        fclose( output ) ;
        output=NULL ;
    }
}

void cmlfile::set_labels ( long olabels )
{
    labels=olabels ;
    stop=new keyword[olabels] ;
    compiled=0 ;
}
void cmlfile::set_sections ( long osections )
{
    sections=osections ;
    block=new section[osections] ;
    compiled=0 ;
}

void cmlfile::set_label_string ( long olabel , char *olabel_string )
{
    strcpy( stop[olabel].word,olabel_string ) ;
    compiled=0 ;
}

void cmlfile::set_section_string ( long olabel , char *osection_string )
{
    strcpy( block[olabel].word,osection_string ) ;
    compiled=0 ;
}

void cmlfile::set_double_format ( char *odouble_format )
{
    strcpy( double_format,odouble_format ) ;
}

void cmlfile::compile ( void )
{
    requirement_insufficient=0 ;
    if ( !input[0] )
    {
        printf( "Nothing to compile.\n" ) ;
        return ;
    }
    long i,j ;

    actual=0 ;

    do
    {
        while ( !feof( input[actual] ))
        {
            fscanf( input[actual],"%s",buffer ) ;
            j=label( buffer ) ;
            switch ( j )
            {
                case None:
                    break ;
                case Comment:
                    break ;
                case File:
                    fscanf( input[actual],"%s",buffer ) ;
                    open_next_input( buffer ) ;
                    break ;
                case Section:
                    stat_section() ;
                    break ;
                case EndSection:
                    break ;
                default:
                    stop[j].n++ ;
                    break ;
            }
        }
        actual++ ;
    }
    while ( input[actual] ) ;

    for ( i=0 ; i<labels ; i++ )
    {
        stop[i].file=new long[stop[i].n] ;
        stop[i].position=new long[stop[i].n] ;
    }

    for ( actual=0 ; actual<files ; actual++ )
    {
        rewind( input[actual] ) ;
        while ( !feof( input[actual] ))
        {
            fscanf( input[actual],"%s",buffer ) ;
            j=label( buffer ) ;
            switch ( j )
            {
                case None:
                    break ;
                case Comment:
                    break ;
                case File:
                    read_to_eoln( input[actual],buffer ) ;
                    break ;
                case Section:
                    jump_over_section() ;
                    break ;
                case EndSection:
                    break ;
                default:
                    stop[j].file[stop[j].i]=actual ;
                    stop[j].position[stop[j].i]=ftell( input[actual] ) ;
                    stop[j].i++ ;
                    break ;
            }
        }
    }
    compiled=1 ;
}

void cmlfile::check_requirements ( void )
{
    if ( !compiled )
    {
        compile() ;
    }
    long i ;
    for ( i=0 ; i<labels ; i++ )
    {
        if ( stop[i].n<stop[i].req )
        {
            printf( "Insufficient input for label %s. This means fatal error.\n",stop[i].word ) ;
            requirement_insufficient=1 ;
        }
    }
    for ( i=0 ; i<sections ; i++ )
    {
        if ( block[i].rows_all<block[i].req )
        {
            printf( "Insufficient input for section %s. This means fatal error.\n",block[i].word ) ;
            requirement_insufficient=1 ;
        }
    }
}

void cmlfile::require ( long olabel )
{
    stop[olabel].req=Required ;
}

void cmlfile::minimal ( long olabel , long ominimal )
{
    stop[olabel].req=ominimal ;
}

void cmlfile::optionalize ( long olabel )
{
    stop[olabel].req=Optional ;
}

void cmlfile::refuse ( long olabel )
{
    stop[olabel].req=Refused ;
}

void cmlfile::check_label ( long olabel )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( stop[olabel].n<stop[olabel].req )
    {
        printf( "Insufficient input for label %s. This means fatal error.\n",stop[olabel].word ) ;
        requirement_insufficient=1 ;
    }
}

void cmlfile::set_minimal_rows ( long osection , long ominimal )
{
    block[osection].req=ominimal ;
}

void cmlfile::check_section ( long osection )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( block[osection].rows_all<block[osection].req )
    {
        printf( "Insufficient input for section %s. This means fatal error.\n",block[osection].word ) ;
        requirement_insufficient=1 ;
    }
}

void cmlfile::load_labels ( char *orcname )
{
    FILE *rc ;
    rc=fopen( orcname,"rt" ) ;
    if ( !rc )
    {
        printf( "Cannot open rc file %s.\n",orcname ) ;
        open_error=1 ;
        return ;
    }
    long j,l=0,s=0 ;
    while ( !feof( rc ))
    {
        fscanf( rc,"%s",buffer ) ;
        if ( str_equal( KeywordSection,buffer ))
        {
            s++ ;
            read_to_eoln( rc,buffer ) ;
        }
        else
        {
            if (( buffer[0]>='A' ) && ( buffer[0]<='z' ))
            {
                l++ ;
                read_to_eoln( rc,buffer ) ;
            }
        }
    }
    rewind( rc ) ;
    set_labels( l ) ;
    set_sections( s ) ;

    l=0 ; s=0 ;

    while ( !feof( rc ))
    {
        j=fscanf( rc,"%s",buffer ) ;
        if ( j!=1 )
        {
            break ;
        }
        if ( str_equal( KeywordSection,buffer ))
        {
            fscanf( rc,"%s",buffer ) ;
            fscanf( rc,"%ld",&j ) ;
            set_section_string( s,buffer ) ;
            set_minimal_rows( s,j ) ;
            s++ ;
        }
        else
        {
            if (( buffer[0]>='A' ) && ( buffer[0]<='z' ))
            {
                fscanf( rc,"%ld",&j ) ;
                set_label_string( l,buffer ) ;
                minimal( l,j ) ;
                l++ ;
            }
        }
    }
    compiled=0 ;
}

void cmlfile::inherit_labels ( cmlfile &of )
{
    if ( stop )
    {
        delete [] stop ;
    }
    if ( block )
    {
        delete [] block ;
    }
    labels=of.labels ;
    sections=of.sections ;

    stop=new keyword[labels] ;
    block=new section[sections] ;

    long i ;
    for ( i=0 ; i<labels ; i++ )
    {
        stop[i].copy( of.stop[i] ) ;
    }
    for ( i=0 ; i<sections ; i++ )
    {
        block[i].copy( of.block[i] ) ;
    }
    compiled=0 ;
}

void cmlfile::get_value_count ( long olabel , long &ocount )
{
    if ( !compiled )
    {
        compile() ;
    }
    ocount=stop[olabel].n ;
}

void cmlfile::get_value ( long olabel , long &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoInt ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoInt ;
        read_error=1 ;
        return ;
    }
    actual=stop[olabel].file[0] ;
    fseek( input[actual],stop[olabel].position[0],SEEK_SET ) ;
    long j ;
    j=fscanf( input[actual],"%ld",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value ( long olabel , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }
    actual=stop[olabel].file[0] ;
    fseek( input[actual],stop[olabel].position[0],SEEK_SET ) ;
    long j ;
    j=fscanf( input[actual],"%lf",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value ( long olabel , char *ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue[0]=0 ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue[0]=0 ;
        read_error=1 ;
        return ;
    }
    actual=stop[olabel].file[0] ;
    fseek( input[actual],stop[olabel].position[0],SEEK_SET ) ;
    long j ;
    j=read_string( input[actual],ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value ( long olabel , long oorder , long &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoInt ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoInt ;
        read_error=1 ;
        return ;
    }
    actual=stop[olabel].file[0] ;
    fseek( input[actual],stop[olabel].position[0],SEEK_SET ) ;
    long i,j ;
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],buffer ) ;
    }
    j=fscanf( input[actual],"%ld",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value ( long olabel , long oorder , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }
    actual=stop[olabel].file[0] ;
    fseek( input[actual],stop[olabel].position[0],SEEK_SET ) ;
    long i,j ;
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],buffer ) ;
    }
    j=fscanf( input[actual],"%lf",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value ( long olabel , long oorder , char *ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue[0]=0 ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue[0]=0 ;
        read_error=1 ;
        return ;
    }
    actual=stop[olabel].file[0] ;
    fseek( input[actual],stop[olabel].position[0],SEEK_SET ) ;
    long i,j ;
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],ovalue ) ;
    }
    j=read_string( input[actual],ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_key ( long olabel , long okey , long &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoInt ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoInt ;
        read_error=1 ;
        return ;
    }
    long i=0,pos,key ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        fscanf( input[actual],"%ld",&key ) ;
        if ( key!=okey )
        {
            i++ ;
        }
    }
    while (( key!=okey ) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %ld does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue=NoInt ;
        return ;
    }
    long j ;
    j=fscanf( input[actual],"%ld",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_key ( long olabel , long okey , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }
    long i=0,pos,key ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        fscanf( input[actual],"%ld",&key ) ;
        if ( key!=okey )
        {
            i++ ;
        }
    }
    while (( key!=okey ) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %ld does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue=NoDouble ;
        return ;
    }
    long j ;
    j=fscanf( input[actual],"%lf",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_key ( long olabel , long okey , char *ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue[0]=0 ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue[0]=0 ;
        read_error=1 ;
        return ;
    }
    long i=0,pos,key ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        fscanf( input[actual],"%ld",&key ) ;
        if ( key!=okey )
        {
            i++ ;
        }
    }
    while (( key!=okey ) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %ld does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue[0]=0 ;
        return ;
    }
    long j ;
    j=read_string( input[actual],ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_key ( long olabel , long okey , long oorder , long &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoInt ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoInt ;
        read_error=1 ;
        return ;
    }
    long i=0,pos,key ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        fscanf( input[actual],"%ld",&key ) ;
        if ( key!=okey )
        {
            i++ ;
        }
    }
    while (( key!=okey ) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %ld does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue=NoInt ;
        return ;
    }
    long j ;
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],buffer ) ;
    }
    j=fscanf( input[actual],"%ld",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_key ( long olabel , long okey , long oorder , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }
    long i=0,pos,key ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        fscanf( input[actual],"%ld",&key ) ;
        if ( key!=okey )
        {
            i++ ;
        }
    }
    while (( key!=okey ) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %ld does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue=NoDouble ;
        return ;
    }
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],buffer ) ;
    }
    long j ;
    j=fscanf( input[actual],"%lf",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_key ( long olabel , long okey , long oorder , char *ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue[0]=0 ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue[0]=0 ;
        read_error=1 ;
        return ;
    }
    long i=0,pos,key ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        fscanf( input[actual],"%ld",&key ) ;
        if ( key!=okey )
        {
            i++ ;
        }
    }
    while (( key!=okey ) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %ld does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue[0]=0 ;
        return ;
    }
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],ovalue ) ;
    }
    long j ;
    j=read_string( input[actual],ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_key ( long olabel , char *okey , long &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoInt ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoInt ;
        read_error=1 ;
        return ;
    }
    long i=0,pos ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        read_string( input[actual],buffer ) ;
        if ( !str_equal( okey,buffer ))
        {
            i++ ;
        }
    }
    while (( !str_equal( okey,buffer )) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %s does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue=NoInt ;
        return ;
    }
    long j ;
    j=fscanf( input[actual],"%ld",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_key ( long olabel , char *okey , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }
    long i=0,pos ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        read_string( input[actual],buffer ) ;
        if ( !str_equal( okey,buffer ))
        {
            i++ ;
        }
    }
    while (( !str_equal( okey,buffer )) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %s does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue=NoDouble ;
        return ;
    }
    long j ;
    j=fscanf( input[actual],"%lf",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_key ( long olabel , char *okey , char *ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue[0]=0 ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue[0]=0 ;
        read_error=1 ;
        return ;
    }
    long i=0,pos ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        read_string( input[actual],buffer ) ;
        if ( !str_equal( okey,buffer ))
        {
            i++ ;
        }
    }
    while (( !str_equal( okey,buffer )) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %s does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue[0]=0 ;
        return ;
    }
    long j ;
    j=read_string( input[actual],ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_key ( long olabel , char *okey , long oorder , long &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoInt ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoInt ;
        read_error=1 ;
        return ;
    }
    long i=0,pos ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        read_string( input[actual],buffer ) ;
        if ( !str_equal( okey,buffer ))
        {
            i++ ;
        }
    }
    while (( !str_equal( okey,buffer )) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %s does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue=NoInt ;
        return ;
    }
    long j ;
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],buffer ) ;
    }
    j=fscanf( input[actual],"%ld",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_key ( long olabel , char *okey , long oorder , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }
    long i=0,pos ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        read_string( input[actual],buffer ) ;
        if ( !str_equal( okey,buffer ))
        {
            i++ ;
        }
    }
    while (( !str_equal( okey,buffer )) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %s does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue=NoDouble ;
        return ;
    }
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],buffer ) ;
    }
    long j ;
    j=fscanf( input[actual],"%lf",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_key ( long olabel , char *okey , long oorder , char *ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue[0]=0 ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue[0]=0 ;
        read_error=1 ;
        return ;
    }
    long i=0,pos ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        read_string( input[actual],buffer ) ;
        if ( !str_equal( okey,buffer ))
        {
            i++ ;
        }
    }
    while (( !str_equal( okey,buffer )) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %s does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue[0]=0 ;
        return ;
    }
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],ovalue ) ;
    }
    long j ;
    j=read_string( input[actual],ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_index ( long olabel , long oindex , long &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoInt ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoInt ;
        read_error=1 ;
        return ;
    }
    if ( oindex>=stop[olabel].n )
    {
        printf( "Index over possible value.\n" ) ;
        ovalue=NoInt ;
        read_error=1 ;
        return ;
    }

    actual=stop[olabel].file[oindex] ;
    fseek( input[actual],stop[olabel].position[oindex],SEEK_SET ) ;

    long j ;
    j=fscanf( input[actual],"%ld",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_index ( long olabel , long oindex , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }
    if ( oindex>=stop[olabel].n )
    {
        printf( "Index over possible value.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }

    actual=stop[olabel].file[oindex] ;
    fseek( input[actual],stop[olabel].position[oindex],SEEK_SET ) ;

    long j ;
    j=fscanf( input[actual],"%lf",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_index ( long olabel , long oindex , char *ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue[0]=0 ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue[0]=0 ;
        read_error=1 ;
        return ;
    }
    if ( oindex>=stop[olabel].n )
    {
        printf( "Index over possible value.\n" ) ;
        ovalue[0]=0 ;
        read_error=1 ;
        return ;
    }

    actual=stop[olabel].file[oindex] ;
    fseek( input[actual],stop[olabel].position[oindex],SEEK_SET ) ;

    long j ;
    j=read_string( input[actual],ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_index ( long olabel , long oindex , long oorder , long &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoInt ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoInt ;
        read_error=1 ;
        return ;
    }
    if ( oindex>=stop[olabel].n )
    {
        printf( "Index over possible value.\n" ) ;
        ovalue=NoInt ;
        read_error=1 ;
        return ;
    }

    actual=stop[olabel].file[oindex] ;
    fseek( input[actual],stop[olabel].position[oindex],SEEK_SET ) ;

    long i,j ;
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],buffer ) ;
    }
    j=fscanf( input[actual],"%ld",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_index ( long olabel , long oindex , long oorder , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }
    if ( oindex>=stop[olabel].n )
    {
        printf( "Index over possible value.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }

    actual=stop[olabel].file[oindex] ;
    fseek( input[actual],stop[olabel].position[oindex],SEEK_SET ) ;

    long i,j ;
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],buffer ) ;
    }
    j=fscanf( input[actual],"%lf",&ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_value_with_index ( long olabel , long oindex , long oorder , char *ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue[0]=0 ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue[0]=0 ;
        read_error=1 ;
        return ;
    }
    if ( oindex>=stop[olabel].n )
    {
        printf( "Index over possible value.\n" ) ;
        ovalue[0]=0 ;
        read_error=1 ;
        return ;
    }

    actual=stop[olabel].file[oindex] ;
    fseek( input[actual],stop[olabel].position[oindex],SEEK_SET ) ;

    long i,j ;
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],ovalue ) ;
    }
    j=read_string( input[actual],ovalue ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_percentage ( long olabel , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        read_error=1 ;
        return ;
    }
    actual=stop[olabel].file[0] ;
    fseek( input[actual],stop[olabel].position[0],SEEK_SET ) ;
    long j ;
    j=fscanf( input[actual],"%s",buffer ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
    sscanf( buffer,"%lf",&ovalue ) ;
    j=0 ;
    while ( buffer[j] )
    {
        if ( buffer[j]=='%' )
        {
            ovalue/=100.0 ;
            return ;
        }
        j++ ;
    }
}

void cmlfile::get_percentage ( long olabel , long oorder , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }
    actual=stop[olabel].file[0] ;
    fseek( input[actual],stop[olabel].position[0],SEEK_SET ) ;
    long i,j ;
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],buffer ) ;
    }
    j=fscanf( input[actual],"%s",buffer ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
    sscanf( buffer,"%lf",&ovalue ) ;
    j=0 ;
    while ( buffer[j] )
    {
        if ( buffer[j]=='%' )
        {
            ovalue/=100.0 ;
            return ;
        }
        j++ ;
    }
}

void cmlfile::get_percentage_with_key ( long olabel , long okey , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }
    long i=0,pos,key ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        fscanf( input[actual],"%ld",&key ) ;
        if ( key!=okey )
        {
            i++ ;
        }
    }
    while (( key!=okey ) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %ld does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue=NoDouble ;
        return ;
    }
    long j ;
    j=fscanf( input[actual],"%s",buffer ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
    sscanf( buffer,"%lf",&ovalue ) ;
    j=0 ;
    while ( buffer[j] )
    {
        if ( buffer[j]=='%' )
        {
            ovalue/=100.0 ;
            return ;
        }
        j++ ;
    }
}

void cmlfile::get_percentage_with_key ( long olabel , long okey , long oorder , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }
    long i=0,pos,key ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        fscanf( input[actual],"%ld",&key ) ;
        if ( key!=okey )
        {
            i++ ;
        }
    }
    while (( key!=okey ) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %ld does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue=NoDouble ;
        return ;
    }
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],buffer ) ;
    }
    long j ;
    j=fscanf( input[actual],"%s",buffer ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
    sscanf( buffer,"%lf",&ovalue ) ;
    j=0 ;
    while ( buffer[j] )
    {
        if ( buffer[j]=='%' )
        {
            ovalue/=100.0 ;
            return ;
        }
        j++ ;
    }
}

void cmlfile::get_percentage_with_key ( long olabel , char *okey , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }
    long i=0,pos ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        fscanf( input[actual],"%s",buffer ) ;
        if ( !str_equal( okey,buffer ))
        {
            i++ ;
        }
    }
    while (( !str_equal( okey,buffer )) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %s does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue=NoDouble ;
    }
    long j ;
    j=fscanf( input[actual],"%s",buffer ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
    sscanf( buffer,"%lf",&ovalue ) ;
    j=0 ;
    while ( buffer[j] )
    {
        if ( buffer[j]=='%' )
        {
            ovalue/=100.0 ;
            return ;
        }
        j++ ;
    }
}

void cmlfile::get_percentage_with_key ( long olabel , char *okey , long oorder , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }
    long i=0,pos ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        fscanf( input[actual],"%s",buffer ) ;
        if ( !str_equal( okey,buffer ))
        {
            i++ ;
        }
    }
    while (( !str_equal( okey,buffer )) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %s does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        ovalue=NoDouble ;
    }
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],buffer ) ;
    }
    long j ;
    j=fscanf( input[actual],"%s",buffer ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
    sscanf( buffer,"%lf",&ovalue ) ;
    j=0 ;
    while ( buffer[j] )
    {
        if ( buffer[j]=='%' )
        {
            ovalue/=100.0 ;
            return ;
        }
        j++ ;
    }
}

void cmlfile::get_percentage_with_index ( long olabel , long oindex , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }
    if ( oindex>=stop[olabel].n )
    {
        printf( "Index over possible value.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }

    actual=stop[olabel].file[oindex] ;
    fseek( input[actual],stop[olabel].position[oindex],SEEK_SET ) ;

    long j ;
    j=fscanf( input[actual],"%s",buffer ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
    sscanf( buffer,"%lf",&ovalue ) ;
    j=0 ;
    while ( buffer[j] )
    {
        if ( buffer[j]=='%' )
        {
            ovalue/=100.0 ;
            return ;
        }
        j++ ;
    }
}

void cmlfile::get_percentage_with_index ( long olabel , long oindex , long oorder , double &ovalue )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        ovalue=NoDouble ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }
    if ( oindex>=stop[olabel].n )
    {
        printf( "Index over possible value.\n" ) ;
        ovalue=NoDouble ;
        read_error=1 ;
        return ;
    }

    actual=stop[olabel].file[oindex] ;
    fseek( input[actual],stop[olabel].position[oindex],SEEK_SET ) ;

    long i,j ;
    for ( i=0 ; i<oorder ; i++ )
    {
        read_string( input[actual],buffer ) ;
    }
    j=fscanf( input[actual],"%s",buffer ) ;

    if ( j!=1 )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
    sscanf( buffer,"%lf",&ovalue ) ;
    j=0 ;
    while ( buffer[j] )
    {
        if ( buffer[j]=='%' )
        {
            ovalue/=100.0 ;
            return ;
        }
        j++ ;
    }
}

void cmlfile::get_line ( long olabel , char *oline )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        oline[0]=0 ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        oline[0]=0 ;
        read_error=1 ;
        return ;
    }
    actual=stop[olabel].file[0] ;
    fseek( input[actual],stop[olabel].position[0],SEEK_SET ) ;
    long j ;
    j=read_to_eoln( input[actual],oline ) ;
    if ( j==None )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_line_with_key ( long olabel , long okey , char *oline )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        oline[0]=0 ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        oline[0]=0 ;
        read_error=1 ;
        return ;
    }
    long i=0,pos,key ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        fscanf( input[actual],"%ld",&key ) ;
        if ( key!=okey )
        {
            i++ ;
        }
    }
    while (( key!=okey ) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %ld does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        oline[0]=0 ;
        return ;
    }
    long j ;
    j=read_to_eoln( input[actual],oline ) ;

    if ( j==None )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_line_with_key ( long olabel , char *okey , char *oline )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        oline[0]=0 ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        oline[0]=0 ;
        read_error=1 ;
        return ;
    }
    long i=0,pos ;
    do
    {
        actual=stop[olabel].file[i] ;
        pos=stop[olabel].position[i] ;
        fseek( input[actual],pos,SEEK_SET ) ;
        fscanf( input[actual],"%s",buffer ) ;
        if ( !str_equal( okey,buffer ))
        {
            i++ ;
        }
    }
    while (( !str_equal( okey,buffer )) && ( i<stop[olabel].n )) ;
    if ( i==stop[olabel].n )
    {
        printf( "Key value %s does not occur in data for label %s.\n",okey,stop[olabel].word ) ;
        oline[0]=0 ;
        return ;
    }
    long j ;
    j=read_to_eoln( input[actual],oline ) ;

    if ( j==None )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::get_line_with_index ( long olabel , long oindex , char *oline )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( !stop[olabel].n )
    {
        oline[0]=0 ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        oline[0]=0 ;
        read_error=1 ;
        return ;
    }
    if ( oindex>=stop[olabel].n )
    {
        printf( "Index over possible value.\n" ) ;
        oline[0]=0 ;
        read_error=1 ;
        return ;
    }

    actual=stop[olabel].file[oindex] ;
    fseek( input[actual],stop[olabel].position[oindex],SEEK_SET ) ;

    long j ;
    j=read_to_eoln( input[actual],oline ) ;

    if ( j==None )
    {
        printf( "Seems that some error occured reading data for label %s.\n",stop[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::find_section ( long olabel , long &osection_rows )
{
    if ( !compiled )
    {
        compile() ;
    }
    if ( block[olabel].files==0 )
    {
        osection_rows=0 ;
        return ;
    }
    if ( normal_read_disable )
    {
        printf( "Cannot read line input now because now reading SECTION data.\n" ) ;
        read_error=1 ;
    }
    block[olabel].actual=0 ;
    block[olabel].i=0 ;
    actual=block[olabel].file[0] ;
    osection_rows=block[olabel].rows_all ;
    fseek( input[actual],block[olabel].pos[0],SEEK_SET ) ;
    read_to_eoln( input[actual],buffer ) ;
    normal_read_disable=1 ;
}

void cmlfile::get_section_line ( long olabel , char *oline )
{
    if ( !compiled )
    {
        compile() ;
    }
    long i=block[olabel].actual ;
    if ( block[olabel].i>=block[olabel].rows[i] )
    {
        i++ ;
        if ( i>=files )
        {
            normal_read_disable=0 ;
            return ;
        }
        block[olabel].actual=i ;
        block[olabel].i=0 ;
        actual=block[olabel].file[i] ;
        fseek( input[actual],block[olabel].pos[i],SEEK_SET ) ;
        read_to_eoln( input[actual],buffer ) ;
    }
    long j ;
    j=read_to_eoln( input[actual],oline ) ;
    block[olabel].i++ ;

    if ( j==None )
    {
        printf( "Seems that some error occured reading data for label %s.\n",block[olabel].word ) ;
        read_error=1 ;
    }
}

void cmlfile::end_reading_section ( void )
{
    normal_read_disable=0 ;
}

void cmlfile::out ( long olabel , long ovalue )
{
    if ( !output )
    {
        printf( "Cannot write to output file.\n" ) ;
        write_error=1 ;
    }
    fprintf( output,"%s %ld\n",stop[olabel].word,ovalue ) ;
}

void cmlfile::out ( long olabel , double ovalue )
{
    if ( !output )
    {
        printf( "Cannot write to output file.\n" ) ;
        write_error=1 ;
    }
    fprintf( output,"%s ",stop[olabel].word ) ;
    fprintf( output,double_format,ovalue ) ;
    fprintf( output,"\n" ) ;
}


void cmlfile::out ( long olabel , char *ovalue )
{
    if ( !output )
    {
        printf( "Cannot write to output file.\n" ) ;
        write_error=1 ;
    }
    fprintf( output,"%s ",stop[olabel].word ) ;
    long i=0,k=0 ;
    while ( ovalue[i] )
    {
        if ( ovalue[i]==' ' )
        {
            k=1 ;
            break ;
        }
        i++ ;
    }
    if ( k )
    {
        fputc( 34,output ) ;
    }
    fprintf( output,"%s",ovalue ) ;
    if ( k )
    {
        fputc( 34,output ) ;
    }
    fprintf( output,"%c",' ' ) ;
    fprintf( output,"\n" ) ;
}

void cmlfile::out_label ( long olabel )
{
    if ( !output )
    {
        printf( "Cannot write to output file.\n" ) ;
        write_error=1 ;
    }
    fprintf( output,"%s ",stop[olabel].word ) ;
}

void cmlfile::out_long ( long ovalue )
{
    if ( !output )
    {
        printf( "Cannot write to output file.\n" ) ;
        write_error=1 ;
    }
    fprintf( output,"%ld ",ovalue ) ;
}

void cmlfile::out_double ( double ovalue )
{
    if ( !output )
    {
        printf( "Cannot write to output file.\n" ) ;
        write_error=1 ;
    }
    fprintf( output,double_format,ovalue ) ;
    fprintf( output," " ) ;
}

void cmlfile::out_string ( char *ovalue )
{
    if ( !output )
    {
        printf( "Cannot write to output file.\n" ) ;
        write_error=1 ;
    }
    long i=0,k=0 ;
    while ( ovalue[i] )
    {
        if ( ovalue[i]==' ' )
        {
            k=1 ;
            break ;
        }
        i++ ;
    }
    if ( k )
    {
        fputc( 34,output ) ;
    }
    fprintf( output,"%s",ovalue ) ;
    if ( k )
    {
        fputc( 34,output ) ;
    }
    fprintf( output,"%c",' ' ) ;
}

void cmlfile::out_free_string ( char *ovalue )
{
    if ( !output )
    {
        printf( "Cannot write to output file.\n" ) ;
        write_error=1 ;
    }
    fprintf( output,"%s",ovalue ) ;
}

void cmlfile::out_eoln ( void )
{
    if ( !output )
    {
        printf( "Cannot write to output file.\n" ) ;
        write_error=1 ;
    }
    fprintf( output,"\n" ) ;
}

void cmlfile::out_section ( long olabel )
{
    if ( !output )
    {
        printf( "Cannot write to output file.\n" ) ;
        write_error=1 ;
    }
    fprintf( output,"%s %s\n",KeywordSection,block[olabel].word ) ;
}

void cmlfile::out_data ( char *oline )
{
    if ( !output )
    {
        printf( "Cannot write to output file.\n" ) ;
        write_error=1 ;
    }
    fprintf( output,"%s\n",oline ) ;
}

void cmlfile::out_section_end ( void )
{
    if ( !output )
    {
        printf( "Cannot write to output file.\n" ) ;
        write_error=1 ;
    }
    fprintf( output,"%s\n",KeywordEndSection ) ;
}

void cmlfile::clear_errors ( void )
{
    open_error=0 ;
    read_error=0 ;
    write_error=0 ;
}

long cmlfile::error_opening_file ( void )
{
    return open_error ;
}

long cmlfile::error_reading_file ( void )
{
    return read_error ;
}

long cmlfile::error_writing_file ( void )
{
    return write_error ;
}

long cmlfile::error_in_requirements ( void )
{
    return requirement_insufficient ;
}

long cmlfile::any_error ( void )
{
    if ( open_error | read_error | write_error | requirement_insufficient )
    {
        return 1 ;
    }
    return 0 ;
}

void cmlfile::open_next_input ( char *ofname )
{
    input[files]=fopen( ofname,"rt" ) ;
    if ( !input[files] )
    {
        printf( "Cannot open file %s.\n",ofname ) ;
        open_error=1 ;
        return ;
    }
    files++ ;
}

void cmlfile::stat_section ( void )
{
    fscanf( input[actual],"%s",buffer ) ;
    long i,j,b ;
    j=section_name( buffer ) ;
    if ( j==None )
    {
        printf( "Invalid section %s.\n",buffer ) ;
        return ;
    }
    i=block[j].files ;
    block[j].pos[i]=ftell( input[actual] ) ;
    block[j].end[i]=ftell( input[actual] ) ;
    block[j].file[i]=actual ;
    block[j].rows[i]=0 ;

    b=read_to_eoln( input[actual],buffer ) ;
    b=read_to_eoln( input[actual],buffer ) ;
    if ( b==None )
    {
        printf( "Seems that some error occured reading data for section %s.\n",block[j].word ) ;
        read_error=1 ;
        return ;
    }

    while ( !str_equal( KeywordEndSection,buffer ))
    {
        block[j].end[i]=ftell( input[actual] ) ;
        block[j].rows[i]++ ;
        read_to_eoln( input[actual],buffer ) ;
        if ( b==None )
        {
            printf( "Seems that some error occured reading data for section %s.\n",block[j].word ) ;
            read_error=1 ;
            return ;
        }

    }
    block[j].rows_all+=block[j].rows[i] ;
    block[j].files++ ;
}

void cmlfile::jump_over_section ( void )
{
    fscanf( input[actual],"%s",buffer ) ;
    long j ;
    j=section_name( buffer ) ;
    if ( j==None )
    {
        printf( "Invalid section %s.\n",buffer ) ;
        return ;
    }
    long i=0 ;
    while ( block[j].file[i]!=actual )
    {
        i++ ;
    }
    fseek( input[actual],block[j].end[i],SEEK_SET ) ;
}

long cmlfile::label ( char *obuffer )
{
    if ( obuffer[0]=='#' )
    {
        return Comment ;
    }
    if ( str_equal( KeywordSection,obuffer ))
    {
        return Section ;
    }
    if ( str_equal( KeywordEndSection,obuffer ))
    {
        return EndSection ;
    }
    if ( str_equal( KeywordFile,obuffer ))
    {
        return File ;
    }
    long i ;
    for ( i=0 ; i<labels ; i++ )
    {
        if ( str_equal( obuffer,stop[i].word ))
        {
            return i ;
        }
    }
    return None ;
}

long cmlfile::section_name ( char *obuffer )
{
    long i ;
    for ( i=0 ; i<sections ; i++ )
    {
        if ( str_equal( block[i].word,obuffer ))
        {
            return i ;
        }
    }
    return None ;
}
