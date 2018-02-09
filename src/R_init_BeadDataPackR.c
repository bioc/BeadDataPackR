#include "BeadDataPackR.h"

static const R_CallMethodDef callMethods[] = 
{
    {"roundLocsFileValues", (DL_FUNC)&roundLocsFileValues, 1},
    {"composeIntensityFlags", (DL_FUNC)&composeIntensityFlags, 2},
    {"int2Bits", (DL_FUNC)&int2Bits, 1},
    {"decodeInd", (DL_FUNC)&decodeInd, 1},
    {"adjustValues", (DL_FUNC)&adjustValues, 1},
    {"returnTrueIndex", (DL_FUNC)&returnTrueIndex, 3},
    {"bitsToInts", (DL_FUNC)&bitsToInts, 1},
    {NULL, NULL, 0}
    
};

void R_init_BeadDataPackR(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, TRUE);
}
