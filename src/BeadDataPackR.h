#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


SEXP roundLocsFileValues(SEXP inputVector);
SEXP composeIntensityFlags(SEXP neg, SEXP large);
SEXP int2Bits(SEXP flags);
SEXP decodeInd(SEXP indices);
SEXP adjustValues(SEXP mat);
SEXP returnTrueIndex(SEXP predX, SEXP predY, SEXP nrow);

