#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>
#include <fstream>

extern bool CLEAN;

//extern int MIN_EVIDENCE; // min #evidence reads in Sites.cpp when doing Export

extern int NMEI;

extern int WIN;

extern int MIN_CLIP;

extern int MATCH;

extern int MIN_QUALITY;

extern int NFIELD;

extern std::string REF_CHR;

extern std::string PATH;

extern std::ofstream LOG;

extern int N_NEG_REF;

#endif