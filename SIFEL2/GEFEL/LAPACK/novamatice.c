#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define spINSIDE_SPARSE
#include "spConfig.h"
#include "spmatrix.h"
char ProgramName[] = "Program";

MatrixPtr Matrix = NULL;
int cteniMatice(int *ind, int *jnd, double* nenul, double *prava, char *nazev, int *pocet)
{
	int n, i, j, *iind = ind;
	double d;
	char radek[80];
	FILE *soubor = fopen(nazev, "r");
	fgets(radek, 80, soubor);
	fgets(radek, 80, soubor);
	sscanf(radek, "%d%", &n);
	while(1)
	{
		fscanf(soubor, "%d%d%lf", &i, &j, &d);
		if (i == 0 && j == 0)
			break;
		else
		{
			*iind++ = i; *jnd++ = j; *nenul++ = d;
		}
	}
	for (i = 0; i < n; i++)
		fscanf(soubor, "%lf", prava++);
	fclose(soubor);
	*pocet = (int)(iind - ind);
	return n;
}


int novaMatice(int n, int *ind, int *jnd, double* nenul, int rozm)
{
	int Error, *lastInd = ind + n;
	ElementPtr pElement;
	do
    {
		if(*ind > rozm || *jnd > rozm)
		{
			fprintf( stderr, "Špatné indexy %d %d \n", *ind, *jnd);
			exit(1);
        }
        pElement = spGetElement(Matrix, *ind, *jnd);
        if (pElement == NULL)
        {
			fprintf(stderr, "Insufficient memory available.\n");
            exit(1);
        }
		pElement->Real = *nenul++;
//		pElement->pInitInfo = nenul++;
		jnd++;
    } while (++ind <lastInd);

/*    rozm = spGetSize(Matrix, 1); */
    if ((Error = spError(Matrix)) != spOKAY)
    {
//		PrintMatrixErrorMessage( Error );
        if(Error >= spFATAL) return 1;
    }
    return 0;
}
static int Init(ElementPtr pElement, double *pInitInfo)
{
	pElement->Real = *pInitInfo;
    return 0;
}
void main()
{
	int pocet, Error, ret, m, n;
	int *indi = (int*)malloc(60 * sizeof(int));
	int *indj = (int*)malloc(60 * sizeof(int));
	double *hodn = (double*)malloc(60 * sizeof(double));
	double *prava = (double*)malloc(25 * sizeof(double));
	double *Solution = (double*)malloc(25 * sizeof(double));
	int rozm = cteniMatice(indi, indj, hodn, prava, "mat2.txt", &pocet);
	Matrix = spCreate(rozm, &Error);
	ret = novaMatice(pocet, indi, indj, hodn, rozm);
    if ((m = spGetSize(Matrix, 1)) != (n = spGetSize(Matrix, 0)))
		return;
   //spInitialize(Matrix, Init);
   Error = spFactor(Matrix);
   spSolve( Matrix, prava, Solution);
}