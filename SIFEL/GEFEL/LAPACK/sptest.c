/*
 *  TEST MODULE for the sparse matrix routines
 *
 *  Author:                     Advisor:
 *     Kenneth S. Kundert           Alberto Sangiovanni-Vincentelli
 *     UC Berkeley
 *
 *  This file contains the test routine for the sparse matrix routines.
 *  They are able to read matrices from files and solve them.
 *
 *  >>> Functions contained in this file:
 *  main
 *  ReadMatrixFromFile
 *  PrintMatrixErrorMessage
 *  InterpretCommandLine
 *  GetProgramName
 *  Usage
 
 *  Revision and copyright information.
 *
 *  Copyright (c) 1985,86,87,88
 *  by Kenneth S. Kundert and the University of California.
 *
 *  Permission to use, copy, modify, and distribute this software and
 *  its documentation for any purpose and without fee is hereby granted,
 *  provided that the copyright notices appear in all copies and
 *  supporting documentation and that the authors and the University of
 *  California are properly credited.  The authors and the University of
 *  California make no representations as to the suitability of this
 *  software for any purpose.  It is provided `as is', without express
 *  or implied warranty.
 *
 *  IMPORTS
 *
 *  >>> Import descriptions:
 *  stdio.h math.h ctype.h
 *      Standard C libraries.
 *  spConfig.h
 *      Macros that customize the sparse matrix package. It is not normally
 *      necessary, nor is normally particularly desirable to include this
 *      file into the calling routines.  Nor should spINSIDE_SPARSE be defined.
 *      It is done in this test file so that the complex test routines may be
 *      removed when they don't exist in Sparse.
 *  spmatrix.h
 *      Macros and declarations to be imported by the user.
 */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#define spINSIDE_SPARSE
#include "spConfig.h"
#undef spINSIDE_SPARSE
#include "spmatrix.h"

static struct Bin {
    ComplexVector Array;
    struct Bin     *Next;
}  *FirstArray, *CurrentArray = NULL;
/*
 *  COMPLEX ARRAY ALLOCATION
 *
 *  These routines are used to check out and in arrays of complex numbers.
 */
static void CheckInAllComplexArrays()
{
    if (CurrentArray != (void*)0)
        CurrentArray = FirstArray;
    return;
}

static int ReadMatrixFromFile();
static void PrintMatrixErrorMessage(int  Error);
static void EnlargeVectors(int NewSize, int Clear );
/*
 *   GLOBAL DECLARATIONS
 *
 *   These variables, types, and macros are used throughout this file.
 *
 *   >>> Macros
 *   PRINT_LIMIT
 *       The maximum number of terms to be printed form the solution vector.
 *       May be overridden by the -n option.
 *
 *   >>> Variable types:
 *   ElementRecord
 *       A structure for holding data for each matrix element.
 *
 *   >>> Global variables:
 *  ProgramName  (char *)
 *      The name of the test program with any path stripped off.
 *  FileName  (char *)
 *      The name of the current file name.
 *  Complex  (BOOLEAN)
 *      The flag that indicates whether the matrix is complex or real.
 *  Element  (ElementRecord[])
 *      Array used to hold information about all the elements.
 *  Matrix  (char *)
 *      The pointer to the matrix.
 *  Size  (int)
 *      The size of the matrix.
 *  RHS  (double)
 *      The right-hand side vector (b in Ax = b).
 *  Solution  (double)
 *      The solution vector (x in Ax = b).
 *  RHS_Verif  (double)
 *      The calculated RHS vector.
 */

/* Begin global declarations. */
#define  PRINT_LIMIT  9

static char *ProgramName, *FileName;
static int  Size, MaxSize = 0, PrintLimit = PRINT_LIMIT, Iterations = 1, ColumnAsRHS = 1;
static BOOLEAN  Complex, SolutionOnly = 0, RealAsComplex = 0, Transposed = 0;
static BOOLEAN  CreatePlotFiles = 0, UseColumnAsRHS = 0, PrintLimitSet = 0;
static BOOLEAN ExpansionWarningGiven, DiagPivoting = DIAGONAL_PIVOTING;
static char Message[BUFSIZ];
MatrixPtr Matrix = NULL;
static double AbsThreshold = 0.0, RelThreshold = 0.0;
static double *RHS, *Solution, *RHS_Verif, *iRHS, *iSolution, *iRHS_Verif;
static ComplexVector cRHS, cSolution, cRHS_Verif;
int Init();
static void Usage(char *sBadOpt)
{
    if (sBadOpt != NULL)
		fprintf(stderr, "%s: unknown or deformed option `%s'.\n", ProgramName, sBadOpt);
    fprintf(stderr, "\nUsage: %s [options] [file names]\n\nOptions:\n", ProgramName);
    fprintf(stderr, "    -s      Print solution rather than run statistics.\n");
    fprintf(stderr, "    -r x    Use x as relative threshold.\n");
    fprintf(stderr, "    -a x    Use x as absolute threshold.\n");
    fprintf(stderr, "    -n n    Print first n terms of solution vector.\n");
    fprintf(stderr, "    -i n    Repeat build/factor/solve n times.\n");
    fprintf(stderr, "    -b n    Use n'th column of matrix as b in Ax=b.\n");
#if DIAGONAL_PIVOTING
    fprintf(stderr, "    -c      Use complete (as opposed to diagonal) pivoting.\n");
#endif
#if DOCUMENTATION
    fprintf(stderr, "    -p      Create plot files `filename.bef' and `filename.aft'.\n");
#endif
#if spCOMPLEX
    fprintf(stderr, "    -x      Treat real matrix as complex with imaginary part zero.\n");
#endif
#if TRANSPOSE
    fprintf(stderr, "    -t      Solve transposed system.\n");
#endif
    fprintf(stderr, "    -u      Print usage message.\n");
    exit(1);
}

static int InterpretCommandLine(int  ArgCount, char* Args[])
{
	int I, FileCount = 0;
	char *GetProgramName();

/* Begin `InterpretCommandLine'. */

/* Determine the name of the program. */
    char* ProgramName = GetProgramName(Args[0]);

/* Step through the argument list, interpreting and deleting the options. */
    for (I = 1; I < ArgCount; I++)
    {   if (!strcmp(Args[I], "-a"))
        {
			if (ArgCount == I || (sscanf(Args[I+1],"%lf", &AbsThreshold) != 1))
            {
				AbsThreshold = 0.0;
                Usage(Args[I]);
            }
            else I++;
        }
        else if (!strcmp(Args[I], "-r"))
        {
			if (ArgCount == I || (sscanf(Args[I+1], "%lf", &RelThreshold) != 1))
            {
				RelThreshold = 0.0;
                Usage(Args[I]);
            }
            else I++;
        }
        else if (!strcmp(Args[I], "-x"))
        {
#if spCOMPLEX
            RealAsComplex = 1;
#else
            fprintf(stderr, "%s: Sparse is not configured to solve complex matrices.\n",
                    ProgramName);
            fprintf(stderr,"    Enable spCOMPLEX in `spConfig.h'.\n");
#endif
        }
        else if (!strcmp(Args[I], "-s"))
            SolutionOnly = 1;
        else if (!strcmp(Args[I], "-c"))
            DiagPivoting = 0;
        else if (!strcmp(Args[I], "-t"))
        {
#if TRANSPOSE
            Transposed = 1;
#else
            fprintf(stderr, "%s: Sparse is not configured to solve transposed system.\n",
                    ProgramName);
            fprintf(stderr,"    Enable TRANSPOSE in `spConfig.h'.\n");
#endif
        }
        else if (!strcmp(Args[I], "-n"))
        {
			if (ArgCount == I || (sscanf(Args[I+1],"%ld", &PrintLimit) != 1))
            {
				PrintLimit = PRINT_LIMIT;
                Usage(Args[I]);
            }
            else
            {
				PrintLimitSet = 1; I++;
            }
        }
        else if (!strcmp(Args[I], "-i"))
        {
			if (ArgCount == I || (sscanf(Args[I+1],"%ld", &Iterations) != 1))
            {
				Iterations = 1; Usage(Args[I]);
            }
            else I++;
        }
        else if (!strcmp(Args[I], "-b"))
        {
			if (ArgCount == I || (sscanf(Args[I+1],"%ld", &ColumnAsRHS) != 1))
            {
				ColumnAsRHS = 1;
                UseColumnAsRHS = 1;
            }
            else
            {
				UseColumnAsRHS = 1;
                ColumnAsRHS = MAX( ColumnAsRHS, 1 );
                I++;
            }
        }
        else if (!strcmp(Args[I], "-p"))
        {
#if DOCUMENTATION
            CreatePlotFiles = 1;
#else
            fprintf(stderr, "%s: Sparse is not configured to generate plot files.\n",
                    ProgramName);
            fprintf(stderr,"    Enable DOCUMENTATION in `spConfig.h'.\n");
#endif
        }
        else if (!strcmp(Args[I], "-u"))
            Usage( (char *)NULL );
        else if (Args[I][0] == '-')
            Usage(Args[I]);
        else Args[++FileCount] = Args[I];
    }
    return FileCount;
}
/*
 *  MAIN TEST ROUTINE for sparse matrix routines
 *
 *  This routine reads a matrix from a file and solves it several
 *  times.  The solution is printed along with some statistics.
 *  The format expected for the input matrix is the same as what is output by
 *  spFileMatrix and spFileVector.
 *
 *  >>> Local variables:
 *  Determinant  (double)
 *      Real portion of the determinant of the matrix.
 *  iDeterminant  (double)
 *      Imaginary portion of the determinant of the matrix.
 *  Error  (int)
 *      Holds error status.
 *  Last  (int)
 *      The index of the last term to be printed of the solution.
 *  Iterations  (int)
 *      The number of times that the matrix will be factored and solved.
 *  MaxRHS  (double)
 *      The largest term in the given RHS vector.
 *  Residual  (double)
 *      The sum of the magnitude of the differences in each corresponding
 *      term of the given and calculated RHS vector.
 *  Exponent  (int)
 *      Exponent for the determinant.
 *  Growth  (double)
 *      The growth that has occurred in the matrix during the factorization.
 *
 *  >>> Timing variables:
 *  BuildTime  (double)
 *      The time required to build up the matrix including the time to clear
 *      the matrix.
 *  ConditionTime  (double)
 *      The time required to compute the condition number.
 *  DeterminantTime  (double)
 *      The time required to compute the determinant.
 *  FactorTime  (double)
 *      The time required to factor the matrix without ordering.
 *  InitialFactorTime  (double)
 *      The time required to factor the matrix with ordering.
 *  SolveTime  (double)
 *      The time required to perform the forward and backward elimination.
 *  StartTime  (double)
 *      The time that a timing interval was started.
 *
 *  >>> Global variables:
 *  Complex  <used>
 *  Matrix  <set>
 *  Size  <used>
 *  RHS  <used>
 *  iRHS  <used>
 */
int main(int ArgCount, char* Args[])
{
	long I;
	int  Error, Last, Elements, Fillins, Exponent;
	double  MaxRHS, Residual, Determinant, iDeterminant,
		ConditionNumber, PseudoCondition,
		LargestBefore, LargestAfter, Roundoff, InfNorm,
		BuildTime, FactorTime, SolveTime, PartitionTime,
		InitialFactorTime, ConditionTime, DeterminantTime;
	BOOLEAN StandardInput;
	char j, PlotFile[BUFSIZ], ErrMsg[BUFSIZ];
	clock_t StartTime, BeginTime;
//	extern char *sbrk();
    BeginTime = clock();
    ArgCount = InterpretCommandLine( ArgCount, Args );
/* Assure that the Sparse is compatible with this test program.*/
/*#if !EXPANDABLE || !INITIALIZE || !ARRAY_OFFSET
        fprintf(stderr, "%s: Sparse is configured inappropriately for test program.\n", ProgramName);
        fprintf(stderr, "    Enable EXPANDABLE, INITIALIZE, and ARRAY_OFFSET in `spConfig.h'.\n");
        exit(1);
#   endif*/
    do
    {
/* Initialization. */
        BuildTime = FactorTime = SolveTime = 0.0;
        ExpansionWarningGiven = 0;
/* Create matrix. */
        Matrix = spCreate(0, spCOMPLEX, &Error );
        if (Matrix == NULL)
        {
			fprintf(stderr, "%s: insufficient memory available.\n", ProgramName);
            exit(1);
        }
        if( Error >= spFATAL ) goto End;

/* Read matrix. */
        if (ArgCount == 0)
            FileName = NULL;
        else
            FileName = *(++Args);
        Error = ReadMatrixFromFile();
        if (Error) goto End;
        StandardInput = (FileName == NULL);

/* Clear solution vector if row and column numbers are not densely packed. */
        if (spGetSize(Matrix, 1) != spGetSize(Matrix, 0))
        {   if (Complex || RealAsComplex)
            {   for (I = Size; I > 0; I--)
                    cSolution[I].Real = cSolution[I].Imag = 0.0;
            }
            else
            {   for (I = Size; I > 0; I--)
                    Solution[I] = 0.0;
            }
        }

/* Perform initial build, factor, and solve. */
        (void)spInitialize(Matrix, Init);

#if MODIFIED_NODAL
        spMNA_Preorder( Matrix );
#endif
#if DOCUMENTATION
        if (CreatePlotFiles)
        {   if (StandardInput)
                (void)sprintf(PlotFile, "bef");
            else
                (void)sprintf(PlotFile, "%s.bef", FileName);
            if (! spFileMatrix( Matrix, PlotFile, FileName, 0, 0, 0 ))
            {   (void)sprintf(ErrMsg,"%s: plotfile `%s'",ProgramName,PlotFile);
                perror( ErrMsg );
            }
        }
#if 0
        spPrint( Matrix, 0 /*reodered*/, 0 /*data*/, 0 /*header*/ );
#endif
#endif /* DOCUMENTATION */
#if STABILITY
        if (! SolutionOnly) LargestBefore = spLargestElement(Matrix);
#endif
#if CONDITION
        if (! SolutionOnly) InfNorm = spNorm(Matrix);
#endif

        StartTime = clock();
        Error = spOrderAndFactor(Matrix, RHS, RelThreshold, AbsThreshold, DiagPivoting);
        InitialFactorTime = clock() - StartTime;
        PrintMatrixErrorMessage( Error );
        if( Error >= spFATAL )
            goto End;

#if DOCUMENTATION
        if (CreatePlotFiles)
        {
			if (StandardInput)
                (void)sprintf(PlotFile, "aft");
            else
                (void)sprintf(PlotFile, "%s.aft", FileName);
            if (! spFileMatrix( Matrix, PlotFile, FileName, 1, 0, 0 ))
            {
				(void)sprintf(ErrMsg,"%s: plotfile `%s'",ProgramName,PlotFile);
                perror( ErrMsg );
            }
        }
#if 0
        spFileStats( Matrix, FileName, "stats" );
#endif
#endif /* DOCUMENTATION */

/*
 * IMAG_VECTORS is a macro that replaces itself with `, iRHS, iSolution'
 * if the options spCOMPLEX and spSEPARATED_COMPLEX_VECTORS are set,
 * otherwise it disappears without a trace.
 */
#if TRANSPOSE
        if (Transposed)
            spSolveTransposed( Matrix, RHS, Solution IMAG_VECTORS );
        else
#endif
            spSolve(Matrix, RHS, Solution IMAG_VECTORS);
        if (SolutionOnly)
            Iterations = 0;
        else
        {
#if STABILITY
            LargestAfter = spLargestElement(Matrix);
            Roundoff = spRoundoff(Matrix, LargestAfter);
#endif
#if CONDITION
            StartTime = clock();
            ConditionNumber = spCondition(Matrix, InfNorm, &Error);
            ConditionTime = clock() - StartTime;
            PrintMatrixErrorMessage(Error);
#endif
#if PSEUDOCONDITION
            PseudoCondition = spPseudoCondition(Matrix);
#endif

            StartTime = clock();
            spPartition(Matrix, spDEFAULT_PARTITION);
            PartitionTime = clock() - StartTime;
        }

/* Solve system of equations Iterations times. */
        for(I = 1; I <= Iterations; I++)
        {
			StartTime = clock();
            (void)spInitialize(Matrix, Init);
            BuildTime += clock() - StartTime;

            StartTime = clock();
            Error = spFactor(Matrix);
            FactorTime += clock() - StartTime;
            if( Error != spOKAY ) PrintMatrixErrorMessage( Error );
            if( Error >= spFATAL ) goto End;
            StartTime = clock();
/*
 * IMAG_VECTORS is a macro that replaces itself with `, iRHS, iSolution'
 * if the options spCOMPLEX and spSEPARATED_COMPLEX_VECTORS are set,
 * otherwise it disappears without a trace.
 */
#if TRANSPOSE
            if (Transposed)
                spSolveTransposed(Matrix, RHS, Solution IMAG_VECTORS );
            else
#endif
                spSolve( Matrix, RHS, Solution IMAG_VECTORS );
            SolveTime += clock() - StartTime;
        }
/* Print Solution. */
        if (SolutionOnly)
        {
			if (PrintLimitSet)
                Last = MIN( PrintLimit, Size );
            else
                Last = Size;
            j = ' ';
        }
        else
        {
			Last = (PrintLimit > Size) ? Size : PrintLimit;
            if (Last > 0) printf("Solution:\n");
            j = 'j';
        }

        if (Complex || RealAsComplex)
        {
#if spSEPARATED_COMPLEX_VECTORS
            for (I = 1; I <= Last; I++)
				printf("%-16.9lg   %-.9lg%c\n", (double)Solution[I], (double)iSolution[I], j);
#else
            for (I = 1; I <= Last; I++)
				printf("%-16.9lg   %-.9lg%c\n", (double)cSolution[I].Real, (double)cSolution[I].Imag, j);
#endif
        }
        else
			for (I = 1; I <= Last; I++)
                printf("%-.9lg\n", (double)Solution[I]);
        if (Last < Size && Last != 0)
            printf("Solution list truncated.\n");
        printf("\n");

#if DETERMINANT
/* Calculate determinant. */
        if (!SolutionOnly)
        {
			StartTime = clock();
#if spCOMPLEX
            spDeterminant( Matrix, &Exponent, &Determinant, &iDeterminant );
#else
            spDeterminant( Matrix, &Exponent, &Determinant );
#endif
            DeterminantTime = clock() - StartTime;
            if (Complex || RealAsComplex)
            {
				Determinant = hypot(Determinant, iDeterminant);
                while (Determinant >= 10.0)
                {
					Determinant *= 0.1; Exponent++;
                }
            }
        }
#else
        Determinant = 0.0; Exponent = 0;
#endif

#if MULTIPLICATION
        if (!SolutionOnly)
        {
/* Calculate difference between actual RHS vector and RHS vector calculated from solution. */
/* Find the largest element in the given RHS vector. */
            MaxRHS = 0.0;
            if (Complex || RealAsComplex)
            {
#if spSEPARATED_COMPLEX_VECTORS
                for (I = 1; I <= Size; I++)
                {
					if (ABS(RHS[I]) > MaxRHS)
                        MaxRHS = ABS(RHS[I]);
                    if (ABS(iRHS[I]) > MaxRHS)
                        MaxRHS = ABS(iRHS[I]);
                }
#else
                for (I = 1; I <= Size; I++)
                {
					if (ABS(cRHS[I].Real) > MaxRHS)
                        MaxRHS = ABS(cRHS[I].Real);
                    if (ABS(cRHS[I].Imag) > MaxRHS)
                        MaxRHS = ABS(cRHS[I].Imag);
                }
#endif
            }
            else
				for (I = 1; I <= Size; I++)
					if (ABS(RHS[I]) > MaxRHS)
						MaxRHS = ABS(RHS[I]);

/* Rebuild matrix. */
            (void)spInitialize(Matrix, Init);
#if spCOMPLEX && spSEPARATED_COMPLEX_VECTORS
            if (Transposed)
				spMultTransposed(Matrix, RHS_Verif, Solution, iRHS_Verif, iSolution);
            else
				spMultiply(Matrix, RHS_Verif, Solution, iRHS_Verif, iSolution);
#else
            if (Transposed)
                spMultTransposed(Matrix, RHS_Verif, Solution);
            else
                spMultiply(Matrix, RHS_Verif, Solution);
#endif

/* Calculate residual. */
            Residual = 0.0;
            if (Complex || RealAsComplex)
            {
#if spSEPARATED_COMPLEX_VECTORS
                for (I = 1; I <= Size; I++)
					Residual += ABS(RHS[I] - RHS_Verif[I]) + ABS(iRHS[I] - iRHS_Verif[I]);
#else
				for (I = 1; I <= Size; I++)
					Residual += ABS(cRHS[I].Real - cRHS_Verif[I].Real) + ABS(cRHS[I].Imag - cRHS_Verif[I].Imag);
#endif
            }
			else
				for (I = 1; I <= Size; I++)
					Residual += ABS(RHS[I] - RHS_Verif[I]);
        }
#endif

/* Print statistics. */
        if (! SolutionOnly)
        {
			Elements = spElementCount(Matrix);
            Fillins = spFillinCount(Matrix);
            
            printf("Initial factor time = %.2lf.\n", InitialFactorTime);
            printf("Partition time = %.2lf.\n", PartitionTime);
            if (Iterations > 0)
            {
				printf("Build time = %.3lf.\n", (BuildTime/Iterations));
                printf("Factor time = %.3lf.\n",(FactorTime/Iterations));
                printf("Solve time = %.3lf.\n", (SolveTime/Iterations));
            }
#if STABILITY
            printf("Condition time = %.2lf.\n", ConditionTime);
#endif
#if DETERMINANT
            printf("Determinant time = %.2lf.\n", DeterminantTime);
#endif
            printf("\nTotal number of elements = %d.\n", Elements);
            printf("Average number of elements per row initially = %.2lf.\n",
                        (double)(Elements - Fillins) / (double)spGetSize(Matrix, 0));
            printf("Total number of fill-ins = %d.\n", Fillins);
#if DETERMINANT || MULTIPLICATION || PSEUDOCONDITION || CONDITION || STABILITY
            putchar('\n');
#endif
#if STABILITY
            if (LargestBefore != 0.0)
                printf("Growth = %.2lg.\n", LargestAfter / LargestBefore);
            printf("Max error in matrix = %.2lg.\n", Roundoff);
#endif
#if STABILITY
            if(ABS(ConditionNumber) > 10 * DBL_MIN);
                printf("Condition number = %.2lg.\n", 1.0 / ConditionNumber);
#endif
#if CONDITION && STABILITY
            printf("Estimated upper bound of error in solution = %.2lg.\n",
                    Roundoff / ConditionNumber);
#endif
#if PSEUDOCONDITION
            printf("PseudoCondition = %.2lg.\n", PseudoCondition);
#endif
#if DETERMINANT
            printf("Determinant = %.3lg", (double)Determinant );
            if (Determinant != 0.0 && Exponent != 0)
                printf("e%d", Exponent);
            putchar('.'); putchar('\n');
#endif
#if MULTIPLICATION
            if (MaxRHS != 0.0)
                printf("Normalized residual = %.2lg.\n", (Residual / MaxRHS));
#endif
        }

End:;
        (void)fflush(stdout);
        CheckInAllComplexArrays();
        spDestroy(Matrix);
        Matrix = NULL;
    } while( --ArgCount > 0 );

    if (!SolutionOnly)
    {
		printf("\nAggregate resource usage:\n");
        printf("    Time required = %.2lf seconds.\n", (double)(clock() - BeginTime) / CLOCKS_PER_SEC );
        //printf("    Virtual memory used = %d kBytes.\n\n", ((int)sbrk(0))/1000);
    }
    return 0;
}
/*
 *   READ MATRIX FROM FILE
 *
 *   This function reads the input file for the matrix and the RHS vector.
 *   If no RHS vector exists, one is created.  If there is an error in the
 *   file, the appropriate error messages are delivered to standard output.
 *
 *   >>> Returned:
 *   The error status is returned.  If no error occurred, a zero is returned.
 *   Otherwise, a one is returned.
 *
 *   >>> Local variables:
 *   pMatrixFile  (FILE *)
 *       The pointer to the file that holds the matrix.
 *   InputString  (char [])
 *       String variable for holding input from the matrix file.
 *   Message  (char [])
 *       String variable that contains a one line descriptor of the matrix.
 *       Descriptor is taken from matrix file.
 *
 *   >>> Global variables:
 *   Complex  <set>
 *   Size  <set>
 *   Element  <set>
 *   RHS  <set>
 *   iRHS  <set>
 */

static int ReadMatrixFromFile()
{
	long I, Reads;
	FILE *pMatrixFile;
	char  sInput[BUFSIZ], sType[BUFSIZ], *p;
	int Error, Row, Col, Count = 0, LineNumber, RHS_Col, IntSize;
	double Real, Imag = 0.0;
	double* pElement;
	ComplexNumber *pValue, *pInitInfo, *CheckOutComplexArray();
	static char *EndOfFile = "%s: unexpected end of file `%s' at line %d.\n",
		*Syntax = "%s: syntax error in file `%s' at line %d.\n";

/* Open matrix file in read mode. */
    if (!SolutionOnly) putchar('\n');

    if (FileName == NULL)
    {
		FileName = "standard input";
        pMatrixFile = stdin;
    }
    else
    {
		pMatrixFile = fopen(FileName, "r");
        if (pMatrixFile == NULL)
        {
			fprintf(stderr, "%s: file %s was not found.\n", ProgramName, FileName);
            return 1;
        }
    }
    Complex = 0;
    LineNumber = 1;

/* Read and print label. */
    if (NULL == fgets(Message, BUFSIZ, pMatrixFile))
    {
		fprintf(stderr, EndOfFile, ProgramName, FileName, LineNumber); return 1;
    }

/* For compatibility with the old file syntax. */
    if (!strncmp( Message, "Starting", 8 ))
    {   /* Test for complex matrix. */
        if (strncmp( Message, "Starting complex", 15 ) == 0)
            Complex = 1;
        LineNumber++;
        if (NULL == fgets( Message, BUFSIZ, pMatrixFile ))
        {
			fprintf(stderr, EndOfFile, ProgramName, FileName, LineNumber);
            return 1;
        }
    }
    if (!SolutionOnly) printf("%-s\n", Message);

/* Read size of matrix and determine type of matrix. */
    LineNumber++;
    if (NULL == fgets( sInput, BUFSIZ, pMatrixFile ))
    {
		fprintf(stderr, EndOfFile, ProgramName, FileName, LineNumber); return 1;
    }
    if ((Reads = sscanf( sInput,"%d %s", &Size, sType )) < 1)
    {
		fprintf(stderr, Syntax, ProgramName, FileName, LineNumber); return 1;
    }
    if (Reads == 2)
    {
		for (p = sType; *p != '\0'; p++)
            if (isupper(*p)) *p += 'a'-'A';
        if (strncmp( sType, "complex", 7 ) == 0)
            Complex = 1;
        else if (strncmp( sType, "real", 7 ) == 0)
            Complex = 0;
        else
        {
			fprintf(stderr, Syntax, ProgramName, FileName, LineNumber); return 1;
        }
    }
    EnlargeVectors( Size, 1 );
    RHS_Col = MIN( Size, ColumnAsRHS );

#if ! spCOMPLEX
    if (Complex)
    {
		fprintf(stderr, "%s: Sparse is not configured to solve complex matrices.\n",
                    ProgramName);
        fprintf(stderr,"    Enable spCOMPLEX in `spConfig.h'.\n");
        return 1;
    }
#endif
#if ! REAL
	if (!(Complex || RealAsComplex))
	{
		fprintf(stderr, "%s: Sparse is not configured to solve real matrices.\n",
				ProgramName);
		fprintf(stderr,"    Enable REAL in `spConfig.h'.\n");
		return 1;
	}
#endif

/* Read matrix elements. */
    do
    {   if (Count == 0)
            pValue = CheckOutComplexArray(Count = 1000);
        LineNumber++;
        if (NULL == fgets( sInput, BUFSIZ, pMatrixFile ))
        {
			fprintf(stderr, "%s: unexpected end of file `%s' at line %d.\n",
                    ProgramName, FileName, LineNumber);
            return 1;
        }
        if (Complex)
        {
			Reads = sscanf( sInput,"%d%d%lf%lf", &Row, &Col, &Real, &Imag );
            if (Reads != 4)
            {   fprintf(stderr, "%s: syntax error in file `%s' at line %d.\n",
                        ProgramName, FileName, LineNumber);
                return 1;
            }
        }
        else
        {
			Reads = sscanf( sInput,"%d%d%lf", &Row, &Col, &Real );
            if (Reads != 3)
            {
				fprintf(stderr, "%s: syntax error in file `%s' at line %d.\n",
                        ProgramName, FileName, LineNumber);
                return 1;
            }
        }
        if(Row < 0 || Col < 0)
        {
			fprintf(stderr, "%s: index not positive in file `%s' at line %d.\n",
                        ProgramName, FileName, LineNumber);
            return 1;
        }
        if(Row > Size || Col > Size)
        {
			if (!ExpansionWarningGiven)
            {
				fprintf( stderr,
         "%s: computed and given matrix size differ in file `%s' at line %d.\n",
                        ProgramName, FileName, LineNumber);
                ExpansionWarningGiven = 1;
            }
            Size = MAX(Row, Col);
            EnlargeVectors( Size, 0 );
        }
        pElement = spGetElement(Matrix, Row, Col);
        if (pElement == NULL)
        {
			fprintf(stderr, "%s: insufficient memory available.\n", ProgramName);
            exit(1);
        }
        pInitInfo = (ComplexNumber*)spGetInitInfo(pElement);
        if (pInitInfo == NULL)
        {
			pValue[--Count].Real = Real;
            pValue[Count].Imag = Imag;
            spInstallInitInfo(pElement, (char *)(pValue + Count));
        }
        else
        {
			pInitInfo->Real += Real;
            pInitInfo->Imag += Imag;
        }

/* Save into RHS vector if in desired column. */
        if (Col == RHS_Col)
        {
			if (Complex || RealAsComplex)
            {
#if spSEPARATED_COMPLEX_VECTORS
                RHS[Row] = Real;
                iRHS[Row] = Imag;
#else
                cRHS[Row].Real = Real;
                cRHS[Row].Imag = Imag;
#endif
            }
            else RHS[Row] = Real;
        }
    } while (Row != 0 && Col != 0);

    Size = spGetSize( Matrix, 1 );
    if (Error = spError( Matrix ) != spOKAY)
    {
		PrintMatrixErrorMessage( Error );
        if(Error >= spFATAL) return 1;
    }

/* Read RHS vector. */
    if (! UseColumnAsRHS && (NULL != fgets( sInput, BUFSIZ, pMatrixFile )))
    {
/* RHS vector exists, read it. */
        LineNumber++;
        for (I = 1; I <= Size; I++)
        {
			if (I != 1 || (strncmp( sInput, "Beginning", 8 ) == 0))
            {
				LineNumber++;
                if (NULL == fgets( sInput, BUFSIZ, pMatrixFile ))
                {
					fprintf(stderr, "%s: unexpected end of file `%s' at line %d.\n",
                            ProgramName, FileName, LineNumber);
                    return 1;
                }
            }
			if (Complex)
            {
#if spSEPARATED_COMPLEX_VECTORS
                Reads = sscanf( sInput,"%lf%lf", &RHS[I], &iRHS[I] );
#else
                Reads = sscanf( sInput, "%lf%lf", &cRHS[I].Real, &cRHS[I].Imag );
#endif
                if (Reads != 2)
                {   fprintf(stderr,
                            "%s: syntax error in file `%s' at line %d.\n",
                            ProgramName, FileName, LineNumber);
                    return 1;
                }
            }
            else /* Not complex. */
            {
				Reads = sscanf( sInput, "%lf", &RHS[I] );
                if (Reads != 1)
                {
					fprintf(stderr, "%s: syntax error in file `%s' at line %d.\n",
                            ProgramName, FileName, LineNumber);
                    return 1;
                }
            }
        }
        if (RealAsComplex && ! Complex)
        {
#if spSEPARATED_COMPLEX_VECTORS
            for (I = 1; I <= Size; I++) iRHS[I] = 0.0;
#else
            for (I = Size; I > 0; I--)
            {
				cRHS[I].Real = RHS[I];
                cRHS[I].Imag = 0.0;
            }
#endif
        }
    }

/* Print out the size and the type of the matrix. */
    if (! SolutionOnly)
    {
		IntSize = spGetSize( Matrix, 0 );
        printf("Matrix is %d x %d ", IntSize, IntSize);
        if (IntSize != Size)
            printf("(external size is %d x %d) ", Size, Size);
        if (Complex || RealAsComplex)
            printf("and complex.\n",Size,Size);
        else
            printf("and real.\n",Size,Size);
    }
	if (Complex || RealAsComplex)
		spSetComplex(Matrix);
	else
		spSetReal(Matrix);
    (void)fclose(pMatrixFile);
    return 0;
}
/*
 *   INITIALIZE MATRIX ELEMENT
 *
 *   This function copys the InitInfo to the Real and Imag matrix element
 *   values.
 *
 *   >>> Returned:
 *   A zero is returns, signifying that no error can be made.
 */

/*ARGSUSED*/

static int Init(double *pElement, char *pInitInfo)
{
    *pElement = *(double*)pInitInfo;
    if (Complex || RealAsComplex)
        *(pElement+1) = *(double *)pInitInfo+1;
    return 0;
}

static ComplexVector CheckOutComplexArray(int Count)
{
	struct Bin Bin, *pBin;
	ComplexVector Temp;

    if (CurrentArray == NULL || CurrentArray->Next == NULL)
    {
		pBin = ALLOC( struct Bin, 1);
        Temp = ALLOC( ComplexNumber, Count );
        if (pBin == NULL || Temp == NULL)
        {
			fprintf( stderr, "%s: insufficient memory available.\n", ProgramName);
            exit(1);
        }
        Bin.Array = Temp; Bin.Next = NULL;
        *pBin = Bin;
        if (CurrentArray == NULL)
            FirstArray = CurrentArray = pBin;
        else
			CurrentArray->Next = CurrentArray = pBin;
    }
    else
    {
		pBin = CurrentArray;
        CurrentArray = pBin->Next;
    }
    return pBin->Array;
}
/*
 *  PRINT ERROR MESSAGE
 *
 *  Prints error message for matrix routines if one is required.
 *
 *  >>> Argument:
 *  Error  <input>  (int)
 *      The error status of the matrix routines.
 */

static void PrintMatrixErrorMessage(int  Error)
{
	int Row, Col;

    if (Error == spOKAY)
        return;
    if (Error >= spFATAL)
		fprintf(stderr, "%s: fatal error detected in file `%s'.\n", ProgramName, FileName );
    else
		fprintf(stderr, "%s: warning in `%s'.\n", ProgramName, FileName );
    if (Error == spNO_MEMORY)
        fprintf(stderr, "    Insufficient memory available.\n");
    else if (Error == spSINGULAR)
    {
		spWhereSingular( Matrix, &Row, &Col );
        printf("    Singular matrix (detected at row %d, column %d).\n", Row, Col );
    }
    else if (Error == spZERO_DIAG)
    {
		spWhereSingular( Matrix, &Row, &Col );
        printf("    Zero diagonal detected at row %d, column %d.\n", Row, Col );
    }
    else if (Error == spPANIC)
        fprintf(stderr, "    Matrix routines being used improperly.\n");
    else if (Error == spSMALL_PIVOT)
        fprintf(stderr, "    A pivot was chosen that is smaller than the absolute threshold.\n");
    return;
}
/*
 *  PROGRAM NAME
 *
 *  Removes path from argv[0] and returns program name.
 *  Assumes UNIX style path names.
 */

static char* GetProgramName(char *Arg0)
{
	char *pTail = Arg0 + strlen(Arg0)-1;
    while (pTail != Arg0 && *(pTail-1) != '/')
        --pTail;
    return pTail;
}

/*
 *  ENLARGE VECTORS
 *
 *  Allocate or enlarge vectors.
 */
static void EnlargeVectors(int  NewSize, int Clear)
{
	int I, PrevSize = MaxSize;
	double* OldRHS = RHS, *iOldRHS = iRHS;
	ComplexVector cOldRHS = cRHS;
/****#ifdef spCOMPLEX
#   define SCALE 2
#else
#   define SCALE 1
#endif *****/
#define SCALE 1

    if (NewSize > PrevSize)
    {
		if (MaxSize != 0)
        {
			free( (char *)Solution );
            free( (char *)RHS_Verif );
        }
        RHS = ALLOC( double, SCALE*(NewSize+1) );
        Solution = ALLOC( double, SCALE*(NewSize+1) );
        RHS_Verif = ALLOC( double, SCALE*(NewSize+1) );
        if (! RHS || ! Solution || ! RHS_Verif)
        {
			fprintf(stderr, "%s: insufficient memory available.\n", ProgramName);
            exit(1);
        }
        cRHS = (ComplexVector)RHS;
        cSolution = (ComplexVector)Solution;
        cRHS_Verif = (ComplexVector)RHS_Verif;
        iRHS = RHS + NewSize + 1;
        iSolution = Solution + NewSize + 1;
        iRHS_Verif = RHS_Verif + NewSize + 1;

/* Copy data from old RHS to new RHS. */
        if (!Clear)
/* Copy old RHS vector to new. */
            if (Complex || RealAsComplex)
            {
				for (I = PrevSize; I > 0; I--)
                {
#if spSEPARATED_COMPLEX_VECTORS || LINT
                    RHS[I] = OldRHS[I];
                    iRHS[I] = iOldRHS[I];
#endif
#if spSEPARATED_COMPLEX_VECTORS || LINT
                    cRHS[I] = cOldRHS[I];
#endif
                }
            }
            else
				for (I = PrevSize; I > 0; I--)
					RHS[I] = OldRHS[I];
        if (MaxSize != 0) free( (char *)OldRHS );
        MaxSize = NewSize;
    }

/* Either completely clear or clear remaining portion of RHS vector. */
    if ((NewSize > PrevSize) || Clear)
    {   if (Clear)
        {
			NewSize = MAX( NewSize, PrevSize );
            PrevSize = 0;
        }
        if (Complex || RealAsComplex)
			for (I = NewSize; I > PrevSize; I--)
            {
#if spSEPARATED_COMPLEX_VECTORS || LINT
                RHS[I] = 0.0;
                iRHS[I] = 0.0;
#endif
#if ! spSEPARATED_COMPLEX_VECTORS || LINT
                cRHS[I].Real = 0.0;
                cRHS[I].Imag = 0.0;
#endif
            }
        else
			for (I = NewSize; I > PrevSize; I--)
                RHS[I] = 0.0;
    }
    return;
}

