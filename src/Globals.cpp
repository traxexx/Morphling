#include "Globals.h"

bool DEBUG_MODE = 0; // default not in debug mode

bool SINGLE_SIDE = 0; // default not printing single-anchor resutls

bool PSEUDO_CHR = 0; // default not discovering pseudo chr

bool PRINT_NON_VARIANT = 0; // do not print

bool APPLY_BOTH_END = 1;

bool APPLY_DEPTH_FILTER = 1; // default use depth filter

bool REFINE_BREAK_POINT = 1; // refine break point?

bool PASS_ONLY = 0;

bool REF_ALLELE = 1; // print ref base out?

int MIN_EVENT_QUALITY = 10;

int WIN = 0;

int STEP = 0;

int MIN_VARIANT_QUALITY = 10;

std::string REF_CHR;

std::string MPATH;

int LEVEL;

int NON_OFFSET;


