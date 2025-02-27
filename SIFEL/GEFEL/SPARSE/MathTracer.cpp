// MathTracer.cpp

#include "MathTracer.h"

DSS_NAMESPASE_BEGIN
////////////////// ////////////////////////////////////////////////////////////////////////////////////
// You can modify this class, or inherit a new one . If you want to redirect the output somewhere else.
MathTracer::MathTracer()
{
  act_block = 0;
  min_pivot = 0;
  stabil_pivot = 1.0;
  break_flag = 0;
}

void MathTracer::Write(double a)
{
  sprintf(m_string,"%e",a);
  Write(m_string);
}

void MathTracer::Write(int a)
{
  sprintf(m_string,"%d",a);
  Write(m_string);
}

void MathTracer::Writeln()
{
  printf("\n");
}

void MathTracer::Writeln(const char* str)
{
  printf("%s", str);
  printf("\n");
}

void MathTracer::Write(const char* str)
{
  printf("%s", str);
}

/*void MathTracer::DrawProgress(double e)
  {
  //e;
  }*/

// true - continue factorization
// false - break factorization
bool MathTracer::CallUnstableDialog()
{
  return true;
}

/*void MathTracer::PrintUnstablePivot(long pivot)
  {
  //pivot;
  }
*/

char* MathTracer::NowString()
{
  memcpy(m_string,"\0",1);
  struct tm *today;
  time_t ltime;
  time(&ltime);
  today = localtime( &ltime );
  if (today)
    //return ctime(&ltime);
    //strftime( m_string, 128,"%C", today );
    sprintf( m_string,"%s", ctime( &ltime ) );
  return m_string;
}

clock_t MathTracer::ClockStart(void)
{
  return clock();
}
char* MathTracer::MeasureClock(clock_t& clock_start)
{
  clock_t end = clock();
  double duration = (double)(end - clock_start) / CLOCKS_PER_SEC;
  sprintf(m_string, "%0.3f s", duration );
  return m_string;
}

void MathTracer::CS()
{
  //temporary_measure_start = CTime::GetCurrentTime();
  time(&m_temporary_measure_start);
  m_clock_start = clock();
}

char* MathTracer::MC_()
{
  clock_t end = clock();
  double duration = (double)(end - m_clock_start) / CLOCKS_PER_SEC;
  sprintf(m_string, "%0.3f s", duration );
  return m_string;
}



MathTracerFile::~MathTracerFile() 
{
  if (out)
  {
    if ((out == stdout) && (out == stderr))
      out = NULL;
    else
      fclose(out);
  }
}

void MathTracerFile::SetOutFile(FILE *outfile) 
{
  if (out)
  {
    if ((out == stdout) && (out == stderr))
      out = NULL;
    else
      fclose(out);
  }
  out = outfile;
}

void MathTracerFile::Writeln()
{
  fprintf(out, "\n");
}

void MathTracerFile::Writeln(const char* str)
{
  fprintf(out, "%s\n", str);
}

void MathTracerFile::Write(const char* str)
{
  fprintf(out, "%s", str);
}

void MathTracerFile::Flush()
{
  fflush(out);
}


DSS_NAMESPASE_END
