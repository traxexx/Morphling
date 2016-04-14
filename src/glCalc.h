#ifndef GLCALC_H
#define GLCALC_H

#include <vector>
#include "Likelihood.h"

using std::vector;

// given log10(gl) and prior, return dosage
int getGtFromGl( vector<float> & log10_prior, vector<float> & gl );

// EM to calculate non-HWE genotype frequency
bool setAFfromGL( vector<float> & gf, vector<lhRec> & vl );

// convert vector of log10(x) to probablity, normalized to sum = 1
void convertToProb( vector<float>  & lp );

#endif
