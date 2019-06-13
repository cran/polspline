#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void heftpq(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void logcensor(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void nlogcensor(void *, void *, void *, void *, void *, void *);
extern void nlogcensorx(void *);
extern void polymarsF(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rpqlsd(void *, void *, void *, void *, void *, void *, void *);
extern void share(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sharex(void *, void *);
extern void sheft(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sheftx(void *);
extern void sphare(void *, void *, void *, void *, void *, void *, void *, void *);
extern void spoly(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void spolyx(void *);
extern void ssumm(void *, void *, void *, void *, void *, void *, void *, void *);
extern void tspsps(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void tspspsx(void *);

static const R_CMethodDef CEntries[] = {
    {"heftpq",      (DL_FUNC) &heftpq,      10},
    {"logcensor",   (DL_FUNC) &logcensor,   13},
    {"nlogcensor",  (DL_FUNC) &nlogcensor,   6},
    {"nlogcensorx", (DL_FUNC) &nlogcensorx,  1},
    {"polymarsF",   (DL_FUNC) &polymarsF,   37},
    {"rpqlsd",      (DL_FUNC) &rpqlsd,       7},
    {"share",       (DL_FUNC) &share,       17},
    {"sharex",      (DL_FUNC) &sharex,       2},
    {"sheft",       (DL_FUNC) &sheft,       16},
    {"sheftx",      (DL_FUNC) &sheftx,       1},
    {"sphare",      (DL_FUNC) &sphare,       8},
    {"spoly",       (DL_FUNC) &spoly,       16},
    {"spolyx",      (DL_FUNC) &spolyx,       1},
    {"ssumm",       (DL_FUNC) &ssumm,        8},
    {"tspsps",      (DL_FUNC) &tspsps,       9},
    {"tspspsx",     (DL_FUNC) &tspspsx,      1},
    {NULL, NULL, 0}
};

void R_init_polspline(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
