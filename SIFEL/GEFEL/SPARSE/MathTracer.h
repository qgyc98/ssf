// MathTracer.h

#ifndef _MATH_TRACER_H__
#define _MATH_TRACER_H__

#include <stdio.h>
#include "DSSAfx.h"

DSS_NAMESPASE_BEGIN
////////////////// ////////////////////////////////////////////////////////////////////////////////////
// You can modify this class, or inherit a new one . If you want to redirect the output somewhere else.
class MathTracer 
{
private:
	char m_string[128];

public:
	double min_pivot;
	double stabil_pivot;
	int break_flag;
	long act_block;
	long act_row;

	MathTracer();

	virtual void Write(double a);
	virtual void Write(int a);
	virtual void Writeln();
	virtual void Writeln(const char* str);
	virtual void Write(const char* str);

//	virtual void DrawProgress(double e);

	// true - continue factorization
	// false - break factorization
	virtual bool CallUnstableDialog();

//	virtual void PrintUnstablePivot(long pivot);

	char* NowString();
	void CS(void);
	char* MC_();

	clock_t ClockStart(void);
	char* MeasureClock(clock_t& clock);

protected:
	time_t m_temporary_measure_start;
	clock_t m_clock_start;
};


/**
  Derived class from MathTracer which redirects the solver output to the text file
  By default, the output is sent to the standard error file.
 */
class MathTracerFile : public MathTracer
{
 private:
  FILE *out; // pointer to the opened text file where the output should be sent 

	// silence compiler warning about the listed method hiding
	using MathTracer::Write;
	using MathTracer::Writeln;

 public:
  MathTracerFile() : MathTracer() {out = stderr;};
  MathTracerFile(FILE *outfile) : MathTracer() {out = outfile;};
  ~MathTracerFile();

  virtual void SetOutFile(FILE *outfile);
  virtual void Writeln();
  virtual void Writeln(const char* str);
  virtual void Write(const char* str);
  virtual void Flush();
};


class MathTracerNull : public MathTracer
{
 private:
	// silence compiler warning about the listed method hiding
	using MathTracer::Write;
	using MathTracer::Writeln;

public:
  MathTracerNull() : MathTracer() {};

  virtual void Writeln() {};
  virtual void Writeln(char*) {};
  virtual void Write(char*) {};
};
DSS_NAMESPASE_END

#endif //_MATH_TRACER_H__
