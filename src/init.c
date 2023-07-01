#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(klhsc)(double *Kmat, int *nobs, double *y, int *nlam, double *ulam, double *tol, int *maxit, double *gam, int *anlam, int *npass, int *jerr, double *alpmat);
extern void F77_NAME(llhsc)(double *Xmat, int *nobs, int *np, double *y, int *nlam, double *ulam, double *tol, int *maxit, double *gam, int *anlam, int *npass, int *jerr, double *btmat);


static const R_FortranMethodDef FortranEntries[] = {
    {"klhsc",  (DL_FUNC) &F77_NAME(klhsc),  12},
    {"llhsc",  (DL_FUNC) &F77_NAME(llhsc),  13},
    {NULL, NULL, 0}
};

void R_init_lkysvm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
