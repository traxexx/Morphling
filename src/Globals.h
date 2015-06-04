#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>

extern bool DEBUG_MODE; // if set 1 print detailed message in functions
extern bool SINGLE_SIDE; // if include single-end results in printing vcf
extern bool PSEUDO_CHR; // if 1, include pseudo chromosomes
extern bool PRINT_NON_VARIANT; // if 1, print GT=0 windows. For debug use
extern bool APPLY_BOTH_END;
extern bool APPLY_DEPTH_FILTER; // filter variants by depth
extern bool REFINE_BREAK_POINT; // refine break point?
extern bool REF_ALLELE; // print ref base out?
extern bool PASS_ONLY; // only print pass variant?

extern bool PARALLEL; // if TURE, print re-genotype command.

extern int WIN;  // win length
extern int STEP; // step length
extern int MIN_VARIANT_QUALITY; // min variant quality
extern int LEVEL; // minimum level
extern int NON_OFFSET; // minimum # non-evidence reads = LEVEL + NON_OFFSET
extern int MIN_EVENT_QUALITY;

extern std::string REF_CHR;
extern std::string MPATH; // path of Morphling ( with / at last )

#endif