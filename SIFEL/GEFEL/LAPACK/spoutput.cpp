/*
 *  MATRIX OUTPUT MODULE
 *
 *  Author:                     Advisor:
 *      Kenneth S. Kundert          Alberto Sangiovanni-Vincentelli
 *      UC Berkeley
 *
 *  This file contains the output-to-file and output-to-screen routines for
 *  the matrix package.
 *
 *  >>> User accessible functions contained in this file:
 *  spPrint
 *  spFileMatrix
 *  spFileVector
 *  spFileStats
 *
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
 */
#define spINSIDE_SPARSE
#include "spconfig.h"
#include "spmatrix.h"
/*
 *  PRINT MATRIX
 *
 *  Formats and send the matrix to standard output.  Some elementary
 *  statistics are also output.  The matrix is output in a format that is
 *  readable by people.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to matrix.
 *  PrintReordered  <input>  (int)
 *      Indicates whether the matrix should be printed out in its original
 *      form, as input by the user, or whether it should be printed in its
 *      reordered form, as used by the matrix routines.  A zero indicates that
 *      the matrix should be printed as inputed, a one indicates that it
 *      should be printed reordered.
 *  Data  <input>  (int)
 *      Boolean flag that when false indicates that output should be
 *      compressed such that only the existence of an element should be
 *      indicated rather than giving the actual value.  Thus 11 times as
 *      many can be printed on a row.  A zero signifies that the matrix
 *      should be printed compressed. A one indicates that the matrix
 *      should be printed in all its glory.
 *  Header  <input>  (int)
 *      Flag indicating that extra information should be given, such as row
 *      and column numbers.
 *
 *  >>> Local variables:
 *  Col  (int)
 *      Column being printed.
 *  ElementCount  (int)
 *      Variable used to count the number of nonzero elements in the matrix.
 *  LargestElement  (double)
 *      The magnitude of the largest element in the matrix.
 *  LargestDiag  (double)
 *      The magnitude of the largest diagonal in the matrix.
 *  Magnitude  (double)
 *      The absolute value of the matrix element being printed.
 *  PrintOrdToIntColMap  (int [])
 *      A translation array that maps the order that columns will be
 *      printed in (if not PrintReordered) to the internal column numbers.
 *  PrintOrdToIntRowMap  (int [])
 *      A translation array that maps the order that rows will be
 *      printed in (if not PrintReordered) to the internal row numbers.
 *  pElement  (ElementPtr)
 *      Pointer to the element in the matrix that is to be printed.
 *  pImagElements  (ElementPtr [ ])
 *      Array of pointers to elements in the matrix.  These pointers point
 *      to the elements whose real values have just been printed.  They are
 *      used to quickly access those same elements so their imaginary values
 *      can be printed.
 *  Row  (int)
 *      Row being printed.
 *  Size  (int)
 *      The size of the matrix.
 *  SmallestDiag  (double)
 *      The magnitude of the smallest diagonal in the matrix.
 *  SmallestElement  (double)
 *      The magnitude of the smallest element in the matrix excluding zero
 *      elements.
 *  StartCol  (int)
 *      The column number of the first column to be printed in the group of
 *      columns currently being printed.
 *  StopCol  (int)
 *      The column number of the last column to be printed in the group of
 *      columns currently being printed.
 *  Top  (int)
 *      The largest expected external row or column number.
 */
void spPrint(MatrixPtr Matrix, int PrintReordered, int Data, int Header)
{
	int  J = 0, I, Row, Col, Size, Top, StartCol = 1, StopCol, Columns, ElementCount = 0;
	int  *PrintOrdToIntRowMap, *PrintOrdToIntColMap;
	double  Magnitude, SmallestDiag, SmallestElement,
		LargestElement = 0.0, LargestDiag = 0.0;
	ElementPtr  pElement, pImagElements[PRINTER_WIDTH/10+1];

    ASSERT( IS_SPARSE( Matrix ) );
    Size = Matrix->Size;

/* Create a packed external to internal row and column translation array. */
# if TRANSLATE
    Top = Matrix->AllocatedExtSize;
#else
    Top = Matrix->AllocatedSize;
#endif
    CALLOC( PrintOrdToIntRowMap, int, Top + 1 );
    CALLOC( PrintOrdToIntColMap, int, Top + 1 );
    if ( PrintOrdToIntRowMap == NULL || PrintOrdToIntColMap == NULL)
    {
		Matrix->Error = spNO_MEMORY; return;
    }
    for (I = 1; I <= Size; I++)
    {
		PrintOrdToIntRowMap[ Matrix->IntToExtRowMap[I] ] = I;
        PrintOrdToIntColMap[ Matrix->IntToExtColMap[I] ] = I;
    }

/* Pack the arrays. */
    for (J = 1, I = 1; I <= Top; I++)
		if (PrintOrdToIntRowMap[I] != 0)
            PrintOrdToIntRowMap[ J++ ] = PrintOrdToIntRowMap[ I ];
    for (J = 1, I = 1; I <= Top; I++)
		if (PrintOrdToIntColMap[I] != 0)
            PrintOrdToIntColMap[ J++ ] = PrintOrdToIntColMap[ I ];
/* Print header. */
    if (Header)
    {
		printf("MATRIX SUMMARY\n\n");
        printf("Size of matrix = %1u x %1u.\n", Size, Size);
        if ( Matrix->Reordered && PrintReordered )
            printf("Matrix has been reordered.\n");
        putchar('\n');
		if ( Matrix->Factored )
            printf("Matrix after factorization:\n");
        else
            printf("Matrix before factorization:\n");
        SmallestElement = DBL_MAX;
        SmallestDiag = SmallestElement;
    }
/* Determine how many columns to use. */
    Columns = PRINTER_WIDTH;
    if (Header) Columns -= 5;
    if (Data) Columns = (Columns+1) / 10;

/*
 * Print matrix by printing groups of complete columns until all the columns
 * are printed.
 */
    J = 0;
    while ( J <= Size )

/* Calculate index of last column to printed in this group. */
    {
		StopCol = StartCol + Columns - 1;
        if (StopCol > Size)
            StopCol = Size;

/* Label the columns. */
        if (Header)
        {
			if (Data)
            {
				printf("    ");
                for (I = StartCol; I <= StopCol; I++)
                {
					if (PrintReordered)
                        Col = I;
                    else
                        Col = PrintOrdToIntColMap[I];
                    printf(" %9d", Matrix->IntToExtColMap[ Col ]);
                }
                printf("\n\n");
            }
            else
            {
				if (PrintReordered)
                    printf("Columns %1d to %1d.\n",StartCol,StopCol);
                else
                {
					printf("Columns %1d to %1d.\n",
                        Matrix->IntToExtColMap[ PrintOrdToIntColMap[StartCol] ],
                        Matrix->IntToExtColMap[ PrintOrdToIntColMap[StopCol] ]);
                }
            }
        }

/* Print every row ...  */
        for (I = 1; I <= Size; I++)
        {
			if (PrintReordered)
                Row = I;
            else
                Row = PrintOrdToIntRowMap[I];

            if (Header)
            {
				if (PrintReordered && !Data)
                    printf("%4d", I);
                else
                    printf("%4d", Matrix->IntToExtRowMap[ Row ]);
                if (!Data) putchar(' ');
            }

/* ... in each column of the group. */
            for (J = StartCol; J <= StopCol; J++)
            {
				if (PrintReordered)
                    Col = J;
                else
                    Col = PrintOrdToIntColMap[J];

                pElement = Matrix->FirstInCol[Col];
                while(pElement != NULL && pElement->Row != Row)
                    pElement = pElement->NextInCol;

                if (Data)
                    pImagElements[J - StartCol] = pElement;

                if (pElement != NULL)

/* Case where element exists */
                {
					if (Data)
					  printf(" %9.3g", (double)pElement->Real);
                    else
                        putchar('x');

/* Update status variables */
                    if ((Magnitude = ELEMENT_MAG(pElement)) > LargestElement)
                        LargestElement = Magnitude;
                    if ((Magnitude < SmallestElement) && (Magnitude != 0.0))
                        SmallestElement = Magnitude;
                    ElementCount++;
                }

/* Case where element is structurally zero */
                else
                {
					if (Data)
                        printf("       ...");
                    else
                        putchar('.');
                }
            }
            putchar('\n');
        }

/* Calculate index of first column in next group. */
        StartCol = StopCol;
        StartCol++;
        putchar('\n');
    }
    if (Header)
    {
		printf("\nLargest element in matrix = %-1.4g.\n", LargestElement);
        printf("Smallest element in matrix = %-1.4g.\n", SmallestElement);

/* Search for largest and smallest diagonal values */
        for (I = 1; I <= Size; I++)
        {
			if (Matrix->Diag[I] != NULL)
            {
				Magnitude = ELEMENT_MAG( Matrix->Diag[I] );
                if ( Magnitude > LargestDiag ) LargestDiag = Magnitude;
                if ( Magnitude < SmallestDiag ) SmallestDiag = Magnitude;
            }
        }

    /* Print the largest and smallest diagonal values */
        if ( Matrix->Factored )
        {
			printf("\nLargest diagonal element = %-1.4g.\n", LargestDiag);
            printf("Smallest diagonal element = %-1.4g.\n", SmallestDiag);
        }
        else
        {
			printf("\nLargest pivot element = %-1.4g.\n", LargestDiag);
            printf("Smallest pivot element = %-1.4g.\n", SmallestDiag);
        }

    /* Calculate and print sparsity and number of fill-ins created. */
        printf("\nDensity = %2.2f%%.\n", ((double)(ElementCount * 100)) /
                                                       ((double)(Size * Size)));
        if (!Matrix->NeedsOrdering)
            printf("Number of fill-ins = %1d.\n", Matrix->Fillins);
    }
    putchar('\n');
    (void)fflush(stdout);

    FREE(PrintOrdToIntColMap);
    FREE(PrintOrdToIntRowMap);
    return;
}
/*
 *  OUTPUT MATRIX TO FILE
 *
 *  Writes matrix to file in format suitable to be read back in by the
 *  matrix test program.
 *
 *  >>> Returns:
 *  One is returned if routine was successful, otherwise zero is returned.
 *  The calling function can query errno (the system global error variable)
 *  as to the reason why this routine failed.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to matrix.
 *  File  <input>  (char *)
 *      Name of file into which matrix is to be written.
 *  Label  <input>  (char *)
 *      String that is transferred to file and is used as a label.
 *  Reordered  <input> (BOOLEAN)
 *      Specifies whether matrix should be output in reordered form,
 *      or in original order.
 *  Data  <input> (BOOLEAN)
 *      Indicates that the element values should be output along with
 *      the indices for each element.  This parameter must be true if
 *      matrix is to be read by the sparse test program.
 *  Header  <input> (BOOLEAN)
 *      Indicates that header is desired.  This parameter must be true if
 *      matrix is to be read by the sparse test program.
 *
 *  >>> Local variables:
 *  Col  (int)
 *      The original column number of the element being output.
 *  pElement  (ElementPtr)
 *      Pointer to an element in the matrix.
 *  pMatrixFile  (FILE *)
 *      File pointer to the matrix file.
 *  Row  (int)
 *      The original row number of the element being output.
 *  Size  (int)
 *      The size of the matrix.
 */
int spFileMatrix(MatrixPtr Matrix, char *File, char *Label, int Reordered, int Data, int Header)
{
	int I, Size, Row, Col, Err;
	ElementPtr pElement;
	FILE *pMatrixFile;
    if ((pMatrixFile = fopen(File, "w")) == NULL)
        return 0;
    Size = Matrix->Size;
    if (Header)
    {
		if (Matrix->Factored && Data)
			Err = fprintf(pMatrixFile, "Warning : The following matrix is factored in to LU form.\n");
        if (Err < 0) return 0;
        if (fprintf(pMatrixFile, "%s\n", Label) < 0) return 0;
        Err = fprintf( pMatrixFile, "%d\t%s\n", Size, (Matrix->Complex ? "complex" : "real"));
        if (Err < 0) return 0;
    }
    if (!Data)
    {   for (I = 1; I <= Size; I++)
        {   pElement = Matrix->FirstInCol[I];
            while (pElement != NULL)
            {   if (Reordered)
                {   Row = pElement->Row;
                    Col = I;
                }
                else
                {   Row = Matrix->IntToExtRowMap[pElement->Row];
                    Col = Matrix->IntToExtColMap[I];
                }
                pElement = pElement->NextInCol;
                if (fprintf(pMatrixFile, "%d\t%d\n", Row, Col) < 0) return 0;
            }
        }
/* Output terminator, a line of zeros. */
        if (Header)
            if (fprintf(pMatrixFile, "0\t0\n") < 0) return 0;
    }
    if (Data && !Matrix->Complex)
    {   for (I = 1; I <= Size; I++)
        {   pElement = Matrix->FirstInCol[I];
            while (pElement != NULL)
            {   Row = Matrix->IntToExtRowMap[pElement->Row];
                Col = Matrix->IntToExtColMap[I];
                Err = fprintf
                (   pMatrixFile,"%d\t%d\t%-.15g\n",
                    Row, Col, (double)pElement->Real
                );
                if (Err < 0) return 0;
                pElement = pElement->NextInCol;
            }
        }
/* Output terminator, a line of zeros. */
        if (Header)
            if (fprintf(pMatrixFile,"0\t0\t0.0\n") < 0) return 0;
    }
    if (fclose(pMatrixFile) < 0) return 0;
    return 1;
}
/*
 *  OUTPUT SOURCE VECTOR TO FILE
 *
 *  Writes vector to file in format suitable to be read back in by the
 *  matrix test program.  This routine should be executed after the function
 *  spFileMatrix.
 *
 *  >>> Returns:
 *  One is returned if routine was successful, otherwise zero is returned.
 *  The calling function can query errno (the system global error variable)
 *  as to the reason why this routine failed.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to matrix.
 *  File  <input>  (char *)
 *      Name of file into which matrix is to be written.
 *  RHS  <input>  (double [])
 *      Right-hand side vector. This is only the real portion if
 *      spSEPARATED_COMPLEX_VECTORS is true.
 *  iRHS  <input>  (double [])
 *      Right-hand side vector, imaginary portion.  Not necessary if matrix
 *      is real or if spSEPARATED_COMPLEX_VECTORS is set false.
 *
 *  >>> Local variables:
 *  pMatrixFile  (FILE *)
 *      File pointer to the matrix file.
 *  Size  (int)
 *      The size of the matrix.
 *
 *  >>> Obscure Macros
 *  IMAG_RHS
 *      Replaces itself with `, iRHS' if the options spCOMPLEX and
 *      spSEPARATED_COMPLEX_VECTORS are set, otherwise it disappears
 *      without a trace.
 */

int spFileVector(MatrixPtr  Matrix,char * File, double* RHS)
{
	int  I, Size/*, Err*/;
	FILE  *pMatrixFile;

    ASSERT( IS_SPARSE( Matrix ) && RHS != NULL)

/* Open File in append mode. */
    if ((pMatrixFile = fopen(File,"a")) == NULL)
        return 0;

/* Correct array pointers for ARRAY_OFFSET. */
#if !ARRAY_OFFSET
#if spCOMPLEX
    if (Matrix->Complex)
    {
#if spSEPARATED_COMPLEX_VECTORS
        ASSERT(iRHS != NULL)
        --RHS; --iRHS;
#else
        RHS -= 2;
#endif
    }
    else
#endif /* spCOMPLEX */
        --RHS;
#endif /* !ARRAY_OFFSET */


/* Output vector. */
    Size = Matrix->Size;
#if REAL && spCOMPLEX
    else
#endif
#if REAL
		for (I = 1; I <= Size; I++)
        if (fprintf(pMatrixFile, "%-.15g\n", (double)RHS[I]) < 0)
                return 0;
#endif /* REAL */

/* Close file. */
    if (fclose(pMatrixFile) < 0) return 0;
    return 1;
}
/*
 *  OUTPUT STATISTICS TO FILE
 *
 *  Writes useful information concerning the matrix to a file.  Should be
 *  executed after the matrix is factored.
 * 
 *  >>> Returns:
 *  One is returned if routine was successful, otherwise zero is returned.
 *  The calling function can query errno (the system global error variable)
 *  as to the reason why this routine failed.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to matrix.
 *  File  <input>  (char *)
 *      Name of file into which matrix is to be written.
 *  Label  <input>  (char *)
 *      String that is transferred to file and is used as a label.
 *
 *  >>> Local variables:
 *  Data  (double)
 *      The value of the matrix element being output.
 *  LargestElement  (double)
 *      The largest element in the matrix.
 *  NumberOfElements  (int)
 *      Number of nonzero elements in the matrix.
 *  pElement  (ElementPtr)
 *      Pointer to an element in the matrix.
 *  pStatsFile  (FILE *)
 *      File pointer to the statistics file.
 *  Size  (int)
 *      The size of the matrix.
 *  SmallestElement  (double)
 *      The smallest element in the matrix excluding zero elements.
 */

int spFileStats(MatrixPtr Matrix, char *File, char *Label)
{
	int  Size, I, NumberOfElements;
	ElementPtr  pElement;
	double Data, LargestElement, SmallestElement;
	FILE  *pStatsFile;
    ASSERT( IS_SPARSE( Matrix ) );

    if ((pStatsFile = fopen(File, "a")) == NULL)
        return 0;

/* Output statistics. */
    Size = Matrix->Size;
    if (!Matrix->Factored)
        fprintf(pStatsFile, "Matrix has not been factored.\n");
    fprintf(pStatsFile, "|||  Starting new matrix  |||\n");
    fprintf(pStatsFile, "%s\n", Label);
    if (Matrix->Complex)
        fprintf(pStatsFile, "Matrix is complex.\n");
    else
        fprintf(pStatsFile, "Matrix is real.\n");
    fprintf(pStatsFile,"     Size = %d\n",Size);

/* Search matrix. */
    NumberOfElements = 0;
    LargestElement = 0.0;
    SmallestElement = DBL_MAX;

    for (I = 1; I <= Size; I++)
    {
		pElement = Matrix->FirstInCol[I];
        while (pElement != NULL)
        {
			NumberOfElements++;
            Data = ELEMENT_MAG(pElement);
            if (Data > LargestElement)
                LargestElement = Data;
            if (Data < SmallestElement && Data != 0.0)
                SmallestElement = Data;
            pElement = pElement->NextInCol;
        }
    }
    SmallestElement = MIN( SmallestElement, LargestElement );
/* Output remaining statistics. */
    fprintf(pStatsFile, "     Initial number of elements = %d\n",
            NumberOfElements - Matrix->Fillins);
    fprintf(pStatsFile,
            "     Initial average number of elements per row = %f\n",
            (double)(NumberOfElements - Matrix->Fillins) / (double)Size);
    fprintf(pStatsFile, "     Fill-ins = %d\n",Matrix->Fillins);
    fprintf(pStatsFile, "     Average number of fill-ins per row = %f%%\n",
            (double)Matrix->Fillins / (double)Size);
    fprintf(pStatsFile, "     Total number of elements = %d\n",
            NumberOfElements);
    fprintf(pStatsFile, "     Average number of elements per row = %f\n",
            (double)NumberOfElements / (double)Size);
    fprintf(pStatsFile,"     Density = %f%%\n",
            (double)(100.0*NumberOfElements)/(double)(Size*Size));
    fprintf(pStatsFile,"     Relative Threshold = %e\n", Matrix->RelThreshold);
    fprintf(pStatsFile,"     Absolute Threshold = %e\n", Matrix->AbsThreshold);
    fprintf(pStatsFile,"     Largest Element = %e\n", LargestElement);
    fprintf(pStatsFile,"     Smallest Element = %e\n\n\n", SmallestElement);
/* Close file. */
    (void)fclose(pStatsFile);
    return 1;
}
